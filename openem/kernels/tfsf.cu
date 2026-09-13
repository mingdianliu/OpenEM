// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Mingdian Liu
// Incident-field correction on the six faces of a TFSF box: a **generic face kernel**. The
// injection axis and the in-plane polarization are folded into the descriptor on the Python side
// (solver._tfsf_box):
//
//   correction = F[base + p*sp + q*sq] += (cb[id] *) coef * row[c0 + cp*p + cq*q]
//
// The table column positions c0/cp/cq are **double**. At oblique incidence, the projection of a
// shell point onto the 1D auxiliary line is still linear in the in-plane scan coordinates (p, q),
// but the coefficients are no longer integers, so the lookup becomes a linear interpolation
// between two adjacent columns. At normal incidence c0/cp/cq are all integral, so pos is exactly
// an integer and w == 0.0, and `row[i0] + 0.0*(row[i1]-row[i0])` is **exactly equal** to row[i0]
// in IEEE754. Hence a single code path, with normal incidence bitwise unchanged; see
// tests/test_tfsf.py.
//
// F is one of the six field component arrays. (p, q) are the two scan axes on the face,
// generalizing the older i/j/k spelling. row is the current row of the 1D auxiliary-grid incident
// table, either the E or the H table, precomputed by tfsf1d.py. c0/cp/cq map a point on the face
// to a table column; cp and cq are non-zero only when a scan axis happens to be the injection
// axis. Every coefficient, sign, 1/dl and polarization component beta included, is derived term by
// term on the Python side from yee.cu's difference expressions (see the comments in
// solver._tfsf_box), so neither direction nor polarization creates a code branch.
//
// At normal incidence with arbitrary in-plane polarization the incident field is
// E = (beta_u*u_hat + beta_v*v_hat)*E1D and H = (beta_u*v_hat - beta_v*u_hat)*H1D, with (u, v, a)
// in cyclic right-handed order: 16 kinds of face correction in all, of which the half with beta=0
// never launch. These kernels launch only when the scene carries a TFSF source, so a scene without
// one is bitwise unchanged on a single code path.

// Table lookup: a real-valued column position with linear interpolation. Out of range is clamped
// to the endpoint with w=0. Normal incidence with c0 = row_n-1 lands exactly here and returns the
// end column itself, the same as an integer index would.
__device__ __forceinline__ double tfsf_tab(
    const double* __restrict__ row, const int row_n,
    const double c0f, const double cpf, const double cqf,
    const int p, const int q)
{
    double pos = c0f + cpf * (double)p + cqf * (double)q;
    int i0 = (int)floor(pos);
    double w = pos - (double)i0;
    if (i0 < 0) { i0 = 0; w = 0.0; }
    else if (i0 >= row_n - 1) { i0 = row_n - 1; w = 0.0; }
    const int i1 = (i0 + 1 < row_n) ? (i0 + 1) : i0;
    return row[i0] + w * (row[i1] - row[i0]);
}

// ---- E-step correction, multiplied by that component's per-cell cb = dt/(eps0*eps) ----
extern "C" __global__ void tfsf_corr_e(
    float* __restrict__ F,
    const float* __restrict__ cbF,
    const float coef,
    const double* __restrict__ tab,     // the whole incident table; current row = tab + (*step)*row_n
    const int* __restrict__ step, const int row_n,
    const double c0f, const double cpf, const double cqf,
    const long long base, const long long sp, const long long sq,
    const int np, const int nq)
{
    const int q = blockIdx.x * blockDim.x + threadIdx.x;
    const int p = blockIdx.y * blockDim.y + threadIdx.y;
    if (p >= np || q >= nq) return;
    const double* row = tab + (size_t)(*step) * row_n;
    const long long id = base + (long long)p * sp + (long long)q * sq;
    F[id] += cbF[id] * coef * (float)tfsf_tab(row, row_n, c0f, cpf, cqf, p, q);
}

// ---- H-step correction; the coefficient already includes dt/mu0 ----
extern "C" __global__ void tfsf_corr_h(
    float* __restrict__ F,
    const float coef,
    const double* __restrict__ tab,
    const int* __restrict__ step, const int row_n,
    const double c0f, const double cpf, const double cqf,
    const long long base, const long long sp, const long long sq,
    const int np, const int nq)
{
    const int q = blockIdx.x * blockDim.x + threadIdx.x;
    const int p = blockIdx.y * blockDim.y + threadIdx.y;
    if (p >= np || q >= nq) return;
    const double* row = tab + (size_t)(*step) * row_n;
    const long long id = base + (long long)p * sp + (long long)q * sq;
    F[id] += coef * (float)tfsf_tab(row, row_n, c0f, cpf, cqf, p, q);
}

