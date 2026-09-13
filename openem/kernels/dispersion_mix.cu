// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Mingdian Liu
// Mixed-averaging cells: Tidy3D's PolarizedAveraging is a **weighted sum** of the arithmetic and
// harmonic averages on a slanted interface,
//
//     eps_eff = (1-beta) * sum_i f_i eps_i + beta / (sum_i f_i/eps_i)
//
// where beta = n_i^2, the squared projection of the interface normal onto that component.
//
// Measured by fitting cell by cell on a sphere at five frequencies: median residual 7.9e-06, with
// beta correlating to n_i^2 at 0.9995. Neither pure arithmetic nor pure harmonic fits any slanted
// interface cell at all.
//
// That **sum** cannot be expressed on either path alone: the arithmetic eps has fixed poles in eps
// (E form), the harmonic eps has fixed poles only in 1/eps (D form), and their sum has fixed poles
// in neither. The way out is to make the harmonic part's D a separate unknown d_h, which leaves
// **one scalar linear equation** per cell per step:
//
//   E        = k2·d_h + SQ,            k2 = ζ_∞ + Σ Re(B_Q)
//   d_arith  = (ε_∞ + Σ2Re(B_P))·E + SP
//   d_total  = (1−β)·d_arith + β·d_h
//     ⇒ d_h = (d_total − k1·SQ − (1−β)·SP) / (k1·k2 + β),   k1 = (1−β)(ε_∞+Σ2Re(B_P))
//
// SP and SQ are the parts of the two pole recursions that **do not involve step n+1**; the pre pass
// computes them from E^n and d_h^n.
//
// **beta=0 degenerates exactly into the E form and beta=1 exactly into the D form.** Both of those
// paths were verified independently, so this kernel inherits their checks (the two in
// tests/test_dispersion.py).
//
// d always stores D/eps0, which keeps both eps and zeta dimensionless.

__device__ __forceinline__ float pick_m(
    const float* __restrict__ Ex, const float* __restrict__ Ey,
    const float* __restrict__ Ez, const int c, const int id)
{
    return (c == 0) ? Ex[id] : ((c == 1) ? Ey[id] : Ez[id]);
}

// First pass: advance both sets of poles to the "half step" using E^n and d_h^n, producing SP and SQ.
extern "C" __global__ void mix_pre(
    float* __restrict__ p_re, float* __restrict__ p_im,     // poles of eps (arithmetic branch)
    float* __restrict__ q_re, float* __restrict__ q_im,     // poles of 1/eps (harmonic branch)
    float* __restrict__ sp, float* __restrict__ sq,         // (n_entry,) outputs
    const float* __restrict__ pa_re, const float* __restrict__ pa_im,  // A_P − 1
    const float* __restrict__ pb_re, const float* __restrict__ pb_im,  // B_P
    const float* __restrict__ qa_re, const float* __restrict__ qa_im,  // A_Q − 1
    const float* __restrict__ qb_re, const float* __restrict__ qb_im,  // B_Q
    const int* __restrict__ p_ofs, const int* __restrict__ q_ofs,
    const int* __restrict__ comp, const int* __restrict__ cell,
    const float* __restrict__ d_h,
    const float* __restrict__ Ex, const float* __restrict__ Ey,
    const float* __restrict__ Ez,
    const int n_entry)
{
    const int e = blockIdx.x * blockDim.x + threadIdx.x;
    if (e >= n_entry) return;

    const float en = pick_m(Ex, Ey, Ez, comp[e], cell[e]);
    // Arithmetic branch: P_half = P + (A_P-1)P + B_P*E^n, SP = sum 2*Re(P_half)
    float s = 0.0f;
    for (int t = p_ofs[e]; t < p_ofs[e + 1]; ++t) {
        const float pr = p_re[t], pi = p_im[t];
        const float dr = pa_re[t] * pr - pa_im[t] * pi + pb_re[t] * en;
        const float di = pa_re[t] * pi + pa_im[t] * pr + pb_im[t] * en;
        p_re[t] = pr + dr;
        p_im[t] = pi + di;
        s += 2.0f * (pr + dr);
    }
    sp[e] = s;
    // Harmonic branch: Q_half = Q + (A_Q-1)Q + B_Q*d_h^n, SQ = sum Re(Q_half) over all roots
    const float dh = d_h[e];
    s = 0.0f;
    for (int t = q_ofs[e]; t < q_ofs[e + 1]; ++t) {
        const float qr = q_re[t], qi = q_im[t];
        const float dr = qa_re[t] * qr - qa_im[t] * qi + qb_re[t] * dh;
        const float di = qa_re[t] * qi + qa_im[t] * qr + qb_im[t] * dh;
        q_re[t] = qr + dr;
        q_im[t] = qi + di;
        s += qr + dr;
    }
    sq[e] = s;
}

// Second pass: the E array now holds the stretched curl plus the source term. Solve for d_h^{n+1}
// and E^{n+1}, then finish both sets of poles.
extern "C" __global__ void mix_post(
    float* __restrict__ Ex, float* __restrict__ Ey, float* __restrict__ Ez,
    float* __restrict__ d_tot, float* __restrict__ d_h,
    float* __restrict__ p_re, float* __restrict__ p_im,
    float* __restrict__ q_re, float* __restrict__ q_im,
    const float* __restrict__ pb_re, const float* __restrict__ pb_im,
    const float* __restrict__ qb_re, const float* __restrict__ qb_im,
    const int* __restrict__ p_ofs, const int* __restrict__ q_ofs,
    const int* __restrict__ comp, const int* __restrict__ cell,
    const float* __restrict__ sp, const float* __restrict__ sq,
    const float* __restrict__ k1, const float* __restrict__ k2,
    const float* __restrict__ one_mb,      // 1 − β
    const float* __restrict__ inv_den,     // 1/(k1·k2 + β)
    const float dt_over_eps0,
    const int n_entry)
{
    const int e = blockIdx.x * blockDim.x + threadIdx.x;
    if (e >= n_entry) return;

    const int c = comp[e], id = cell[e];
    const float curl = pick_m(Ex, Ey, Ez, c, id);
    const float dt_new = d_tot[e] + dt_over_eps0 * curl;
    d_tot[e] = dt_new;

    const float s_p = sp[e], s_q = sq[e];
    const float dh = (dt_new - k1[e] * s_q - one_mb[e] * s_p) * inv_den[e];
    d_h[e] = dh;
    const float en = k2[e] * dh + s_q;

    for (int t = p_ofs[e]; t < p_ofs[e + 1]; ++t) {
        p_re[t] += pb_re[t] * en;
        p_im[t] += pb_im[t] * en;
    }
    for (int t = q_ofs[e]; t < q_ofs[e + 1]; ++t) {
        q_re[t] += qb_re[t] * dh;
        q_im[t] += qb_im[t] * dh;
    }
    if (c == 0)      Ex[id] = en;
    else if (c == 1) Ey[id] = en;
    else             Ez[id] = en;
}
