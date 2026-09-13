// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Mingdian Liu
// One-way plane wave injection, plus the runtime DFT for frequency-domain monitors.
//
// ---- Injection ----------------------------------------------------------
// Single-face TF/SF injection travelling only in +z. Modelled on the FDTDX inheritance chain
// (UniformPlaneSource -> LinearlyPolarizedPlaneSource -> TFSFPlaneSource, where the source is
// added as a correction after the curl update; see update_E/update_H in
// objects/sources/source.py) and on the one-way injection formulas in Meep, src/sources.cpp.
//
// The source plane sits at Ex's z position, edges[ks]. For +z propagation, k >= ks is the
// total-field region and k < ks the scattered-field region; -z propagation mirrors the whole
// thing. The two correction terms are
//   Hy[kh] += coef_h * Ex_inc(t^n)             kh = ks-1 for +z, ks for -z
//   Ex[ks] += cex * coef_e * Hy_inc(t^{n+1/2})
// **Direction creates no code branch**: kh and the signs of both coefficients are computed on the
// Python side and passed in.
// The incident waveform table is sampled on the Python side with source_time.amp_time at the
// **exact retarded times**, folding the half-cell spatial offset into the time argument, with no
// interpolation.
//
// ── DFT ───────────────────────────────────────────────────────────────────
// Modelled on the runtime accumulation in Meep, src/dft.cpp. E lives on whole steps and H on half
// steps, and each enters the phase factor at **its own true time**, so their relative phase comes
// out right on its own with no extra time-colocation correction.
// Spatially H is half a cell off from E along z, so before accumulating it is colocated onto the E
// plane as (H[k]+H[k-1])/2.
// The accumulators are float64: accumulating 68,400 steps in fp32 loses precision.
//
// The phase factor cos/sin(2*pi*f*t) **depends only on (f, t)** and is identical at every cell of
// the monitor. Having each thread recompute nf double-precision sincos pairs is pure waste (270
// million of them per step on the largest scene), so it is hoisted into a small (nf, 4) table:
// dft_phase_table computes it once per step, with the apodization window multiplied in.
// **Bitwise equivalent**: the same sincos, the same multiplications, the same operand order, so
// computing it once gives exactly what computing it a hundred thousand times would.

#define IDX(i, j, k) (((i) * ny + (j)) * nz + (k))

// ---- Injection: H half step ----
extern "C" __global__ void inject_h(
    float* __restrict__ Hy,
    const float coef,          // +-dt/mu0 * inv_dz_p[kh]; the sign follows the direction
    const double* __restrict__ ex_tab,   // (n_total,) the whole Ex_inc table
    const int* __restrict__ step,        // device-side step counter (set up before the graph)
    const int kh,              // H correction plane: ks-1 for +z, ks for -z
    const int nx, const int ny, const int nz)
{
    const int j = blockIdx.x * blockDim.x + threadIdx.x;
    const int i = blockIdx.y * blockDim.y + threadIdx.y;
    if (i >= nx || j >= ny) return;
    // The (float) conversion rounds exactly as the host-side np.float32(table[n]) did, so this is
    // bitwise identical
    const float ex_inc = (float)ex_tab[*step];
    Hy[IDX(i, j, kh)] += coef * ex_inc;
}

// ---- Injection: E whole step ----
extern "C" __global__ void inject_e(
    float* __restrict__ Ex,
    const float* __restrict__ cex,   // dt/(eps0*eps_r), per cell
    const float coef_e,              // +-inv_dz_d[ks]; the sign follows the direction
    const double* __restrict__ hy_tab,   // (n_total,) the whole Hy_inc table
    const int* __restrict__ step,
    const int ks,
    const int nx, const int ny, const int nz)
{
    const int j = blockIdx.x * blockDim.x + threadIdx.x;
    const int i = blockIdx.y * blockDim.y + threadIdx.y;
    if (i >= nx || j >= ny) return;
    const int id = IDX(i, j, ks);
    const float hy_inc = (float)hy_tab[*step];
    Ex[id] += cex[id] * coef_e * hy_inc;
}

// ---- Mode source: H half-step correction ----
// The same thing as the plane wave's inject_h, differing in two ways: the incident field is a
// **two-dimensional complex profile** rather than a constant, and there are two tangential
// components rather than one. The complex multiplication is expanded here as Re(prof*amp), with
// prof varying per point and amp per step, so only two scalars go to the device each step and the
// profile is uploaded once.
// H_t1 and H_t2 are the two tangential components of the normal axis, in cyclic order
// t1 = (axis+1)%3. The profile is (n1, n2) with row width w2; with x as the normal the second
// dimension is z, already padded to nzp = w2.
extern "C" __global__ void inject_mode_h(
    float* __restrict__ Ht1, float* __restrict__ Ht2,
    const float* __restrict__ e1_re, const float* __restrict__ e1_im,
    const float* __restrict__ e2_re, const float* __restrict__ e2_im,
    const double* __restrict__ amp_tab,       // (n,2) per-step (re, im)
    const int* __restrict__ step,
    const float coef_1, const float coef_2,   // carry the sign and dt/mu0/dx
    const int kh, const int axis, const int n1, const int n2, const int w2,
    const int nx, const int ny, const int nz)
{
    const int o2 = blockIdx.x * blockDim.x + threadIdx.x;   // t2, the contiguous dimension
    const int o1 = blockIdx.y * blockDim.y + threadIdx.y;   // t1
    if (o1 >= n1 || o2 >= n2) return;
    const float amp_re = (float)amp_tab[2 * (*step)];
    const float amp_im = (float)amp_tab[2 * (*step) + 1];
    const int o  = o1 * w2 + o2;
    const int id = axis == 0 ? IDX(kh, o1, o2)
                 : axis == 1 ? IDX(o2, kh, o1) : IDX(o1, o2, kh);
    Ht1[id] += coef_1 * (e2_re[o] * amp_re - e2_im[o] * amp_im);
    Ht2[id] += coef_2 * (e1_re[o] * amp_re - e1_im[o] * amp_im);
}

// ---- Mode source: E whole-step correction ----
extern "C" __global__ void inject_mode_e(
    float* __restrict__ Et1, float* __restrict__ Et2,
    const float* __restrict__ ce1, const float* __restrict__ ce2,
    const float* __restrict__ h1_re, const float* __restrict__ h1_im,
    const float* __restrict__ h2_re, const float* __restrict__ h2_im,
    const double* __restrict__ amp_tab,       // (n,2) per-step (re, im)
    const int* __restrict__ step,
    const float coef_1, const float coef_2,   // carry the sign and 1/dx_dual
    const int ks, const int axis, const int n1, const int n2, const int w2,
    const int nx, const int ny, const int nz)
{
    const int o2 = blockIdx.x * blockDim.x + threadIdx.x;
    const int o1 = blockIdx.y * blockDim.y + threadIdx.y;
    if (o1 >= n1 || o2 >= n2) return;
    const float amp_re = (float)amp_tab[2 * (*step)];
    const float amp_im = (float)amp_tab[2 * (*step) + 1];
    const int o  = o1 * w2 + o2;
    const int id = axis == 0 ? IDX(ks, o1, o2)
                 : axis == 1 ? IDX(o2, ks, o1) : IDX(o1, o2, ks);
    Et1[id] += ce1[id] * coef_1 * (h2_re[o] * amp_re - h2_im[o] * amp_im);
    Et2[id] += ce2[id] * coef_2 * (h1_re[o] * amp_re - h1_im[o] * amp_im);
}

// ---- Frequency-domain accumulation ----
// Phasors for 4 components (Ex, Ey, colocated Hx, colocated Hy) times nf frequencies on one z
// plane. Output layout is (comp, f, i, j), with comp in the order Ex, Ey, Hx, Hy.
// The phase is computed inside the kernel with sincos, to avoid shipping four small arrays to the
// device every step: at 68,400 steps and 2 monitors that would be half a million tiny transfers.
// Point current source. Each thread handles one stencil point.
//
//   eps0*eps dE/dt = curl H - J   =>   E -= (dt/(eps0 eps)) * J
//
// J = amp/dV * w, where amp is the time waveform of the current moment, w the multilinear weight
// and dV that point's Yee dual volume; coef[] holds w/dV. cex/cey/cez are dt/(eps0*eps), already
// including 1/(1+k) when lossy, matching how the reference implementation treats dipole injection.
//
// The amplitude is sampled at (n+0.5)*dt, following the reference implementation: electric and
// dipole samples are evaluated at (n+0.5)*dt.
//
// atomicAdd is used because the stencils of two dipoles may land on the same cell.
extern "C" __global__ void inject_dipole(
    float* __restrict__ Ex, float* __restrict__ Ey, float* __restrict__ Ez,
    const float* __restrict__ cex, const float* __restrict__ cey,
    const float* __restrict__ cez,
    const int* __restrict__ flat,      // (n,) flat indices
    const int* __restrict__ comp,      // (n,) 0/1/2
    const float* __restrict__ coef,    // (n,) w/dV
    const int* __restrict__ src_of,    // (n,) which source each point belongs to
    const double* __restrict__ amp,    // (n_src, nsteps1) amplitudes at half steps
    const int* __restrict__ step_p, const int nsteps1, const int n)
{
    const int t = blockIdx.x * blockDim.x + threadIdx.x;
    const int step = *step_p;
    // The step guard is a **backstop**: solver.run already refuses, on the host, any step count
    // beyond the amplitude table. It stays because an out-of-range device read raises nothing and
    // silently injects garbage. What that looked like in practice was the whole field turning NaN a
    // few hundred steps later, which reads as "dispersion instability" and cost a full round of
    // investigation.
    if (t >= n || step >= nsteps1) return;

    float* Ep[3] = {Ex, Ey, Ez};
    const float* cp[3] = {cex, cey, cez};
    const int c = comp[t];
    const int id = flat[t];
    const float a = (float)amp[(size_t)src_of[t] * nsteps1 + step];
    atomicAdd(&Ep[c][id], -cp[c][id] * coef[t] * a);
}

// Point magnetic current source: the H dual of inject_dipole.
//
//   mu0 dH/dt = -curl E - M   =>   H -= (dt/mu0) * M
//
// mu0 is uniform and there is no magnetic loss, so the coefficient is a **scalar**, ch = dt/mu0,
// the same value update_h uses; no per-cell arrays like cex/cey/cez are needed. M = amp/dV * w,
// and coef[] is still w/dV.
//
// The amplitude is sampled at n*dt: the H step advances H^{n-1/2} to H^{n+1/2}, whose midpoint is
// the whole step, matching the convention by which inject_dipole samples at (n+0.5)*dt, the
// midpoint of the E step. The reference implementation states the same: magnetic samples use
// n*dt. The extra phase offset it mentions afterwards belongs only to the E/H pairing of a planar
// Huygens source and does not apply to a point source.
extern "C" __global__ void inject_dipole_h(
    float* __restrict__ Hx, float* __restrict__ Hy, float* __restrict__ Hz,
    const float ch,                    // dt / mu0
    const int* __restrict__ flat,      // (n,) flat indices
    const int* __restrict__ comp,      // (n,) 0/1/2
    const float* __restrict__ coef,    // (n,) w/dV
    const int* __restrict__ src_of,    // (n,) which source each point belongs to
    const double* __restrict__ amp,    // (n_src, nsteps1) amplitudes at whole steps
    const int* __restrict__ step_p, const int nsteps1, const int n)
{
    const int t = blockIdx.x * blockDim.x + threadIdx.x;
    const int step = *step_p;
    if (t >= n || step >= nsteps1) return;   // backstop guard, for the same reason as inject_dipole

    float* Hp[3] = {Hx, Hy, Hz};
    const int c = comp[t];
    const float a = (float)amp[(size_t)src_of[t] * nsteps1 + step];
    atomicAdd(&Hp[c][flat[t]], -ch * coef[t] * a);
}

// Phasors for all 6 components within one index box, **stored at their native Yee positions, not
// colocated**. Colocation is a linear operation that commutes with the DFT, so it is done on the
// host.
// Output layout is (6, nf, ni, nj, nk), with components in the order Ex, Ey, Ez, Hx, Hy, Hz.
//
// apod_e and apod_h are this step's apodization window values, E at t_e and H at t_h (see
// apodization.py). Without apodization they are always 1.0, and multiplying by 1.0 is exact in
// IEEE, so no second code path is needed.
extern "C" __global__ void accumulate_dft_box(
    double* __restrict__ re, double* __restrict__ im,
    const float* __restrict__ Ex, const float* __restrict__ Ey,
    const float* __restrict__ Ez, const float* __restrict__ Hx,
    const float* __restrict__ Hy, const float* __restrict__ Hz,
    const double* __restrict__ tab,     // (nf, 4) slice of the phase table; see dft_phase_table
    const double* __restrict__ win_e,   // (n_total,) apodization window; always 1.0 without it
    const double* __restrict__ win_h,
    const int* __restrict__ stp,
    const int i0, const int j0, const int k0,
    const int ni, const int nj, const int nk, const int nf,
    const int nx, const int ny, const int nz,
    const int f0, const int nf_all, const int tmod)   // ring-table addressing
{
    // Thread layout: flatten the bounding box (j, k) with **k innermost** onto blockIdx.x, put
    // **frequency on blockIdx.y**, and i on blockIdx.z.
    //
    // One thread per cell with nf serialized inside does not fill the GPU: one mode plane had only
    // 4,836 cells against nf=1001. Running threadIdx.x along k instead wastes lanes whenever k is
    // narrow (a box of (3,52,31) uses 31 of 64 lanes, 48%). Flattening covers both: lane
    // utilization is about 90% and the accumulator offset off = a*nj*nk + tid is fully contiguous
    // in tid.
    //
    // Each (component, frequency, cell) bucket is still accumulated in time order by a single
    // thread, so **the summation order is unchanged and the result is bitwise equivalent**. grid.y
    // and grid.z cap at 65535, which bounds nf and ni; the host asserts both.
    const int tid = blockIdx.x * blockDim.x + threadIdx.x;
    const int f = blockIdx.y;
    const int a = blockIdx.z;
    if (a >= ni || tid >= nj * nk || f >= nf) return;
    const int b = tid / nk;
    const int cc = tid - b * nk;

    const int id = IDX(i0 + a, j0 + b, k0 + cc);
    const double vals[6] = {
        (double)Ex[id], (double)Ey[id], (double)Ez[id],
        (double)Hx[id], (double)Hy[id], (double)Hz[id],
    };

    const int cells = ni * nj * nk;
    const int off = (a * nj + b) * nk + cc;

    // The multiplication grouping matches the original code, vals * (cos * apod), so this is
    // bitwise equivalent
    const double apod_e = win_e[*stp], apod_h = win_h[*stp];
    const double* ph = tab + (((size_t)((*stp) % tmod) * nf_all) + f0 + f) * 4;
    const double cE = ph[0] * apod_e, sE = ph[1] * apod_e;
    const double cH = ph[2] * apod_h, sH = ph[3] * apod_h;
    for (int c = 0; c < 6; ++c) {
        const double cr = (c < 3) ? cE : cH;   // E takes the whole-step phase, H the half-step one
        const double si = (c < 3) ? sE : sH;
        const size_t o = (size_t)(c * nf + f) * cells + off;
        re[o] += vals[c] * cr;
        im[o] += vals[c] * si;
    }
}

// Time-domain sampling: write the requested components within the index box into slot, weighted by
// weights. Stored at native Yee positions; colocation and interpolation happen on the host. Output
// layout is (nc, nslots, ni, nj, nk).
//
// weights and add exist for **the time average of H**: E lives natively on whole steps m*dt and H
// on half steps, while the H that Tidy3D reports at t = m*dt is the average of the two adjacent
// half steps (established by measurement). So the solver fires two shots:
//   step m   : weights = [E:1, H:0.5], add=0 (overwrite)
//   step m+1 : weights = [E:0, H:0.5], add=1 (accumulate)
extern "C" __global__ void sample_time_box(
    float* __restrict__ out,
    const float* __restrict__ Ex, const float* __restrict__ Ey,
    const float* __restrict__ Ez, const float* __restrict__ Hx,
    const float* __restrict__ Hy, const float* __restrict__ Hz,
    const int* __restrict__ comps, const float* __restrict__ weights, const int nc,
    const int* __restrict__ step_p,     // device step counter; the sampling time is m = *step+1
    const int beg, const int iv, const int end_, const int nslots,
    const int add,                      // 0 = first shot (E plus half of H), 1 = second (the other half)
    const int i0, const int j0, const int k0,
    const int ni, const int nj, const int nk,
    const int nx, const int ny, const int nz)
{
    // The slot logic moved here from the host: inside the graph the kernel launches unconditionally
    // every step and guards itself. The semantics match the old host version word for word:
    // m_step is the current step plus 1, and the second shot looks at m_step-1.
    const int m_ref = (*step_p + 1) - add;
    const int soff = m_ref - beg;
    if (soff < 0 || m_ref >= end_ || soff % iv) return;
    const int slot = soff / iv;
    if (slot >= nslots) return;
    const int cc = blockIdx.x * blockDim.x + threadIdx.x;
    const int b = blockIdx.y * blockDim.y + threadIdx.y;
    const int a = blockIdx.z;
    if (a >= ni || b >= nj || cc >= nk) return;

    const int id = IDX(i0 + a, j0 + b, k0 + cc);
    const float* F[6] = {Ex, Ey, Ez, Hx, Hy, Hz};
    const int cells = ni * nj * nk;
    const int off = (a * nj + b) * nk + cc;
    for (int c = 0; c < nc; ++c) {
        const float w = weights[c];
        if (w == 0.0f) continue;
        const size_t o = ((size_t)c * nslots + slot) * cells + off;
        const float v = w * F[comps[c]][id];
        out[o] = add ? out[o] + v : v;
    }
}

// Phasors for the 4 tangential components times nf frequencies on one plane, with **any normal axis**.
//
// The two transverse axes are taken in **ascending order**: axis=0 gives (y,z), 1 gives (x,z),
// 2 gives (x,y). The 4 output components are always (E_t1, E_t2, H_t1, H_t2), matching
// term_x = Re(c0*conj(c3)) and term_y = Re(c1*conj(c2)) in flux.plane_flux.
//
// **The Poynting vector carries an extra minus sign when axis=1**: S_y = E_z H_x* - E_x H_z*,
// while ascending (t1,t2) = (x,z) produces E_x H_z* - E_z H_x*. That sign is handled on the host
// in normal_dir, not here.
//
// The phase factors come from dft_phase_table's (nf, 4) table; see the note at the top of this file.
extern "C" __global__ void dft_phase_table(
    double* __restrict__ tab,           // (nf, 4) = [ce, se, ch, sh], **without apodization**
    const double* __restrict__ freqs,   // (nf,) Hz; every DFT monitor's frequencies concatenated
    const double* __restrict__ t_e_tab, // (n_total,) E times, (n+1)*dt
    const double* __restrict__ t_h_tab, // (n_total,) H times, (n+0.5)*dt
    const int* __restrict__ step,
    const int nf, const int tmod,      // ring table; the row is *step % tmod
    float* __restrict__ tab32)         // fp32 copy, used by the batched flush kernels
{
    const int f = blockIdx.x * blockDim.x + threadIdx.x;
    if (f >= nf) return;
    const double t_e = t_e_tab[*step], t_h = t_h_tab[*step];
    const double two_pi = 6.283185307179586476925286766559;
    double ce, se, ch, sh;
    sincos(two_pi * freqs[f] * t_e, &se, &ce);
    sincos(two_pi * freqs[f] * t_h, &sh, &ch);
    // Apodization is **not** multiplied in here. Doing so would tie the table to a monitor and drop
    // back to one launch per monitor per step, which measured 1.76 times slower on scenes with many
    // monitors and few frequencies. It is multiplied inside accumulate instead, with the grouping
    // unchanged.
    double* out = tab + ((size_t)((*step) % tmod) * nf + f) * 4;
    out[0] = ce;
    out[1] = se;
    out[2] = ch;
    out[3] = sh;
    float* out32 = tab32 + ((size_t)((*step) % tmod) * nf + f) * 4;
    out32[0] = (float)ce;
    out32[1] = (float)se;
    out32[2] = (float)ch;
    out32[3] = (float)sh;
}

// A plane monitor has far fewer cells than a volume one (66,000 on one plane of the largest scene),
// so one thread per cell does not fill the GPU. Frequency therefore goes on **blockIdx.z** and each
// thread handles a single (cell, frequency) pair.
// Each (component, frequency, cell) bucket is still accumulated in time order by a single thread,
// so **the summation order is unchanged and the result is bitwise equivalent**. The thread x
// dimension still runs along t2, the contiguous one, so accumulator writes stay coalesced.
// grid.z caps at 65535, which bounds nf; the host asserts it.
extern "C" __global__ void accumulate_dft(
    double* __restrict__ re, double* __restrict__ im,   // (4, nf, n1, n2)
    const float* __restrict__ Ex, const float* __restrict__ Ey,
    const float* __restrict__ Ez,
    const float* __restrict__ Hx, const float* __restrict__ Hy,
    const float* __restrict__ Hz,
    const double* __restrict__ tab,     // (nf, 4) slice of the phase table
    const double* __restrict__ win_e,   // (n_total,) apodization window
    const double* __restrict__ win_h,
    const int* __restrict__ stp,
    const int km,                       // monitor plane: the E index along the normal axis
    const int axis,                     // normal axis, 0/1/2
    const int nf,
    const int nz_log,                   // logical nz, for the transverse extent; nz is the pitch
    const int nx, const int ny, const int nz,
    const int f0, const int nf_all, const int tmod)   // ring-table addressing
{
    const int t1 = (axis == 0) ? 1 : 0;
    const int t2 = (axis == 2) ? 1 : 2;
    const int n1 = (t1 == 0) ? nx : ny;
    const int n2 = (t2 == 1) ? ny : nz_log;

    const int b = blockIdx.x * blockDim.x + threadIdx.x;   // along t2, the contiguous dimension
    const int a = blockIdx.y * blockDim.y + threadIdx.y;   // along t1
    const int f = blockIdx.z;                              // frequency
    if (a >= n1 || b >= n2 || f >= nf) return;

    int c[3];
    c[axis] = km;  c[t1] = a;  c[t2] = b;
    const int id = IDX(c[0], c[1], c[2]);
    c[axis] = km - 1;
    const int idm = IDX(c[0], c[1], c[2]);

    const float* Ep[3] = {Ex, Ey, Ez};
    const float* Hp[3] = {Hx, Hy, Hz};
    // H is half a cell off from E along the normal: colocate it onto the E plane
    const double vals[4] = {
        (double)Ep[t1][id],
        (double)Ep[t2][id],
        0.5 * ((double)Hp[t1][id] + (double)Hp[t1][idm]),
        0.5 * ((double)Hp[t2][id] + (double)Hp[t2][idm]),
    };

    // The multiplication grouping matches the original code, vals * (cos * apod), so this is
    // bitwise equivalent
    const double apod_e = win_e[*stp], apod_h = win_h[*stp];
    const double* ph = tab + (((size_t)((*stp) % tmod) * nf_all) + f0 + f) * 4;
    const double cE = ph[0] * apod_e, sE = ph[1] * apod_e;
    const double cH = ph[2] * apod_h, sH = ph[3] * apod_h;
    const size_t plane = (size_t)n1 * n2;
    const size_t off = (size_t)a * n2 + b;

    for (int cc = 0; cc < 4; ++cc) {
        const double cr = (cc < 2) ? cE : cH;   // E takes the whole-step phase, H the half-step one
        const double si = (cc < 2) ? sE : sH;
        const size_t o = (size_t)(cc * nf + f) * plane + off;
        re[o] += vals[cc] * cr;
        im[o] += vals[cc] * si;
    }
}


// Advance the device-side step counter: after graph capture, zero host work per step.
extern "C" __global__ void step_advance(int* __restrict__ step)
{
    if (blockIdx.x == 0 && threadIdx.x == 0) ++*step;
}


// ---- Batched DFT accumulation ----
// The snapshot kernel copies the **raw f32 field values** each monitor needs into a ring of T rows
// every step. The flush kernel does real work only at a batch boundary, when (*stp+1) % T == 0: it
// replays the multiply-adds of the T steps in the batch in step order into an fp32 partial sum,
// then adds that into the main accumulator in fp64 once at the end of the batch. fp64 throughput
// on an H800 is cut to about 1 TFLOPS and the fp64 multiply-add used to be the arithmetic
// bottleneck of the DFT section. The error bound is about T*eps32, measured at roughly 5e-7
// relative to the spectral peak, so this path is accepted under tolerance rather than bitwise.
// The fp64 accumulator is still read and written once per T steps.
// The tail batch, when a pinned step count is not a multiple of T, is flushed by the host outside
// the loop with force_cnt>0. By then *stp has been advanced one step too far, which the kernel
// subtracts back.

extern "C" __global__ void dft_snap_box(
    float* __restrict__ snap,           // (T, 6, cells)
    const float* __restrict__ Ex, const float* __restrict__ Ey,
    const float* __restrict__ Ez, const float* __restrict__ Hx,
    const float* __restrict__ Hy, const float* __restrict__ Hz,
    const int* __restrict__ stp, const int T,
    const int i0, const int j0, const int k0,
    const int ni, const int nj, const int nk,
    const int nx, const int ny, const int nz)
{
    const int tid = blockIdx.x * blockDim.x + threadIdx.x;
    const int a = blockIdx.z;
    if (a >= ni || tid >= nj * nk) return;
    const int b = tid / nk;
    const int cc = tid - b * nk;
    const int id = IDX(i0 + a, j0 + b, k0 + cc);
    const size_t cells = (size_t)ni * nj * nk;
    const size_t off = (size_t)(a * nj + b) * nk + cc;
    float* row = snap + (size_t)((*stp) % T) * 6 * cells;
    row[0 * cells + off] = Ex[id];
    row[1 * cells + off] = Ey[id];
    row[2 * cells + off] = Ez[id];
    row[3 * cells + off] = Hx[id];
    row[4 * cells + off] = Hy[id];
    row[5 * cells + off] = Hz[id];
}

extern "C" __global__ void dft_flush_box(
    double* __restrict__ re, double* __restrict__ im,
    const float* __restrict__ snap,
    const float* __restrict__ tab,
    const float* __restrict__ win_e, const float* __restrict__ win_h,
    const int* __restrict__ stp,
    const int T, const int force_cnt,
    const int f0, const int nf_all, const int tmod,
    const int cells, const int nf)
{
    // Threads are (cell, group of 8 frequencies). Occupancy is full; the snapshot is streamed once
    // per frequency group as a coalesced read (nf/8 times rather than nf); and the fp64 accumulator
    // is read and written exactly once per batch, which is the lower bound.
    // fp32 partial sums within the batch, added into the main accumulator in fp64 once at the end;
    // accepted under tolerance.
    const int tile = (blockIdx.y + blockIdx.z * gridDim.y);
    const int off = tile * blockDim.x + threadIdx.x;
    const int fc = blockIdx.x * 8;
    if (off >= cells || fc >= nf) return;
    const int s_now = (force_cnt > 0) ? (*stp - 1) : (*stp);
    int cnt = force_cnt;
    if (cnt == 0) {
        if ((s_now + 1) % T != 0) return;
        cnt = T;
    }
    const int s_lo = s_now - cnt + 1;
    const int fhi = (fc + 8 < nf) ? fc + 8 : nf;
    for (int f = fc; f < fhi; ++f) {
        float aR[6], aI[6];
        for (int c = 0; c < 6; ++c) { aR[c] = 0.f; aI[c] = 0.f; }
        for (int q = 0; q < cnt; ++q) {
            const int s = s_lo + q;
            const float apod_e = win_e[s], apod_h = win_h[s];
            const float* ph = tab
                + (((size_t)(s % tmod) * nf_all) + f0 + f) * 4;
            const float cE = ph[0] * apod_e, sE = ph[1] * apod_e;
            const float cH = ph[2] * apod_h, sH = ph[3] * apod_h;
            const float* row = snap + (size_t)((s - s_lo + (s_lo % T)) % T)
                                       * 6 * cells;
            for (int c = 0; c < 6; ++c) {
                const float v = row[c * cells + off];
                aR[c] = fmaf(v, (c < 3) ? cE : cH, aR[c]);
                aI[c] = fmaf(v, (c < 3) ? sE : sH, aI[c]);
            }
        }
        for (int c = 0; c < 6; ++c) {
            const size_t o = (size_t)(c * nf + f) * cells + off;
            re[o] += (double)aR[c]; im[o] += (double)aI[c];
        }
    }
}

extern "C" __global__ void dft_snap_plane(
    float* __restrict__ snap,           // (T, 6, n1, n2)
    const float* __restrict__ Ex, const float* __restrict__ Ey,
    const float* __restrict__ Ez, const float* __restrict__ Hx,
    const float* __restrict__ Hy, const float* __restrict__ Hz,
    const int* __restrict__ stp, const int T,
    const int km, const int axis, const int nz_log,
    const int nx, const int ny, const int nz)
{
    const int t1 = (axis == 0) ? 1 : 0;
    const int t2 = (axis == 2) ? 1 : 2;
    const int n1 = (t1 == 0) ? nx : ny;
    const int n2 = (t2 == 1) ? ny : nz_log;
    const int b = blockIdx.x * blockDim.x + threadIdx.x;
    const int a = blockIdx.y * blockDim.y + threadIdx.y;
    if (a >= n1 || b >= n2) return;
    int c[3];
    c[axis] = km;  c[t1] = a;  c[t2] = b;
    const int id = IDX(c[0], c[1], c[2]);
    c[axis] = km - 1;
    const int idm = IDX(c[0], c[1], c[2]);
    const float* Ep[3] = {Ex, Ey, Ez};
    const float* Hp[3] = {Hx, Hy, Hz};
    const size_t pl = (size_t)n1 * n2;
    const size_t off = (size_t)a * n2 + b;
    float* row = snap + (size_t)((*stp) % T) * 6 * pl;
    row[0 * pl + off] = Ep[t1][id];
    row[1 * pl + off] = Ep[t2][id];
    row[2 * pl + off] = Hp[t1][id];
    row[3 * pl + off] = Hp[t1][idm];
    row[4 * pl + off] = Hp[t2][id];
    row[5 * pl + off] = Hp[t2][idm];
}

extern "C" __global__ void dft_flush_plane(
    double* __restrict__ re, double* __restrict__ im,   // (4, nf, n1, n2)
    const float* __restrict__ snap,
    const float* __restrict__ tab,
    const float* __restrict__ win_e, const float* __restrict__ win_h,
    const int* __restrict__ stp,
    const int T, const int force_cnt,
    const int f0, const int nf_all, const int tmod,
    const int n1, const int n2, const int nf)
{
    // As in dft_flush_box, with the plane points flattened to one dimension; the snapshot and
    // accumulator offsets use only off.
    // fp32 partial sums within the batch, including the 0.5 colocation average, added into the main
    // accumulator in fp64 once at the end.
    const size_t pl = (size_t)n1 * n2;
    const int tile = (blockIdx.y + blockIdx.z * gridDim.y);
    const int off = tile * blockDim.x + threadIdx.x;
    const int fc = blockIdx.x * 8;
    if (off >= (int)pl || fc >= nf) return;
    const int s_now = (force_cnt > 0) ? (*stp - 1) : (*stp);
    int cnt = force_cnt;
    if (cnt == 0) {
        if ((s_now + 1) % T != 0) return;
        cnt = T;
    }
    const int s_lo = s_now - cnt + 1;
    const int fhi = (fc + 8 < nf) ? fc + 8 : nf;
    for (int f = fc; f < fhi; ++f) {
        float aR[4], aI[4];
        for (int c = 0; c < 4; ++c) { aR[c] = 0.f; aI[c] = 0.f; }
        for (int q = 0; q < cnt; ++q) {
            const int s = s_lo + q;
            const float apod_e = win_e[s], apod_h = win_h[s];
            const float* ph = tab
                + (((size_t)(s % tmod) * nf_all) + f0 + f) * 4;
            const float cE = ph[0] * apod_e, sE = ph[1] * apod_e;
            const float cH = ph[2] * apod_h, sH = ph[3] * apod_h;
            const float* row = snap + (size_t)(s % T) * 6 * pl;
            const float vals[4] = {
                row[0 * pl + off],
                row[1 * pl + off],
                0.5f * (row[2 * pl + off] + row[3 * pl + off]),
                0.5f * (row[4 * pl + off] + row[5 * pl + off]),
            };
            for (int c = 0; c < 4; ++c) {
                aR[c] = fmaf(vals[c], (c < 2) ? cE : cH, aR[c]);
                aI[c] = fmaf(vals[c], (c < 2) ? sE : sH, aI[c]);
            }
        }
        for (int c = 0; c < 4; ++c) {
            const size_t o = (size_t)(c * nf + f) * pl + off;
            re[o] += (double)aR[c]; im[o] += (double)aI[c];
        }
    }
}

extern "C" __global__ void dft_flush_box_pf(
    double* __restrict__ re, double* __restrict__ im,
    const float* __restrict__ snap,
    const float* __restrict__ tab,
    const float* __restrict__ win_e, const float* __restrict__ win_h,
    const int* __restrict__ stp,
    const int T, const int force_cnt,
    const int f0, const int nf_all, const int tmod,
    const int cells, const int nf)
{
    // fp32 partial sums within the batch, added into the main accumulator in fp64 once at the end.
    const int off = blockIdx.x * blockDim.x + threadIdx.x;
    const int f = blockIdx.y;
    if (off >= cells || f >= nf) return;
    const int s_now = (force_cnt > 0) ? (*stp - 1) : (*stp);
    int cnt = force_cnt;
    if (cnt == 0) {
        if ((s_now + 1) % T != 0) return;
        cnt = T;
    }
    float aR[6], aI[6];
    for (int c = 0; c < 6; ++c) { aR[c] = 0.f; aI[c] = 0.f; }
    for (int q = 0; q < cnt; ++q) {
        const int s = s_now - cnt + 1 + q;
        const float apod_e = win_e[s], apod_h = win_h[s];
        const float* ph = tab + (((size_t)(s % tmod) * nf_all) + f0 + f) * 4;
        const float cE = ph[0] * apod_e, sE = ph[1] * apod_e;
        const float cH = ph[2] * apod_h, sH = ph[3] * apod_h;
        const float* row = snap + (size_t)(s % T) * 6 * cells;
        for (int c = 0; c < 6; ++c) {
            const float v = row[c * cells + off];
            aR[c] = fmaf(v, (c < 3) ? cE : cH, aR[c]);
            aI[c] = fmaf(v, (c < 3) ? sE : sH, aI[c]);
        }
    }
    for (int c = 0; c < 6; ++c) {
        const size_t o = (size_t)(c * nf + f) * cells + off;
        re[o] += (double)aR[c]; im[o] += (double)aI[c];
    }
}

extern "C" __global__ void dft_flush_plane_pf(
    double* __restrict__ re, double* __restrict__ im,   // (4, nf, n1, n2)
    const float* __restrict__ snap,
    const float* __restrict__ tab,
    const float* __restrict__ win_e, const float* __restrict__ win_h,
    const int* __restrict__ stp,
    const int T, const int force_cnt,
    const int f0, const int nf_all, const int tmod,
    const int n1, const int n2, const int nf)
{
    // fp32 partial sums within the batch, including the 0.5 colocation average, added into the main
    // accumulator in fp64 once at the end.
    const int b = blockIdx.x * blockDim.x + threadIdx.x;
    const int a = blockIdx.y * blockDim.y + threadIdx.y;
    const int f = blockIdx.z;
    if (a >= n1 || b >= n2 || f >= nf) return;
    const int s_now = (force_cnt > 0) ? (*stp - 1) : (*stp);
    int cnt = force_cnt;
    if (cnt == 0) {
        if ((s_now + 1) % T != 0) return;
        cnt = T;
    }
    const size_t pl = (size_t)n1 * n2;
    const size_t off = (size_t)a * n2 + b;
    float aR[4], aI[4];
    for (int c = 0; c < 4; ++c) { aR[c] = 0.f; aI[c] = 0.f; }
    for (int q = 0; q < cnt; ++q) {
        const int s = s_now - cnt + 1 + q;
        const float apod_e = win_e[s], apod_h = win_h[s];
        const float* ph = tab + (((size_t)(s % tmod) * nf_all) + f0 + f) * 4;
        const float cE = ph[0] * apod_e, sE = ph[1] * apod_e;
        const float cH = ph[2] * apod_h, sH = ph[3] * apod_h;
        const float* row = snap + (size_t)(s % T) * 6 * pl;
        const float vals[4] = {
            row[0 * pl + off],
            row[1 * pl + off],
            0.5f * (row[2 * pl + off] + row[3 * pl + off]),
            0.5f * (row[4 * pl + off] + row[5 * pl + off]),
        };
        for (int c = 0; c < 4; ++c) {
            aR[c] = fmaf(vals[c], (c < 2) ? cE : cH, aR[c]);
            aI[c] = fmaf(vals[c], (c < 2) ? sE : sH, aI[c]);
        }
    }
    for (int c = 0; c < 4; ++c) {
        const size_t o = (size_t)(c * nf + f) * pl + off;
        re[o] += (double)aR[c]; im[o] += (double)aI[c];
    }
}


// Shared-memory tiled version of dft_flush_box. The block is 16 cells by 16 frequency threads. The
// threads first cooperatively load a (cnt, 6, 16) snapshot strip into shared memory (30 KB at
// T=80), then each consumes its own (cell, frequency + 16k).
// The multiply-adds for each (cell, frequency) are in exactly the same order as in dft_flush_box,
// so this is bitwise identical, while snapshot DRAM traffic drops from ceil(nf/8) passes to one;
// measured at 2.17x on the benchmark. The signature matches dft_flush_box; the host only changes
// the launch geometry and the dynamic shared-memory size.
#define P49_TILE 16
extern "C" __global__ void dft_flush_box_sm2(
    double* __restrict__ re, double* __restrict__ im,
    const float* __restrict__ snap,
    const float* __restrict__ tab,
    const float* __restrict__ win_e, const float* __restrict__ win_h,
    const int* __restrict__ stp,
    const int T, const int force_cnt,
    const int f0, const int nf_all, const int tmod,
    const int cells, const int nf)
{
    const int s_now = (force_cnt > 0) ? (*stp - 1) : (*stp);
    int cnt = force_cnt;
    if (cnt == 0) { if ((s_now + 1) % T != 0) return; cnt = T; }
    const int s_lo = s_now - cnt + 1;
    const int off0 = blockIdx.x * P49_TILE;
    if (off0 >= cells) return;
    const int tile_n = min(P49_TILE, cells - off0);

    extern __shared__ float sm[];              // (cnt, 6, P49_TILE)
    const int nthr = blockDim.x * blockDim.y;
    const int tid = threadIdx.y * blockDim.x + threadIdx.x;
    const int total = cnt * 6 * P49_TILE;
    for (int i = tid; i < total; i += nthr) {
        const int q = i / (6 * P49_TILE);
        const int c = (i / P49_TILE) % 6;
        const int o = i % P49_TILE;
        const int s = s_lo + q;
        const float* row = snap + (size_t)((s - s_lo + (s_lo % T)) % T)
                                   * 6 * cells;
        sm[i] = (o < tile_n) ? row[(size_t)c * cells + off0 + o] : 0.f;
    }
    __syncthreads();

    const int o = threadIdx.x;
    if (o >= tile_n) return;
    const int off = off0 + o;
    for (int f = threadIdx.y; f < nf; f += blockDim.y) {
        float aR[6], aI[6];
        for (int c = 0; c < 6; ++c) { aR[c] = 0.f; aI[c] = 0.f; }
        for (int q = 0; q < cnt; ++q) {
            const int s = s_lo + q;
            const float apod_e = win_e[s], apod_h = win_h[s];
            const float* ph = tab + (((size_t)(s % tmod) * nf_all) + f0 + f) * 4;
            const float cE = ph[0] * apod_e, sE = ph[1] * apod_e;
            const float cH = ph[2] * apod_h, sH = ph[3] * apod_h;
            const float* v6 = sm + (q * 6) * P49_TILE;
            for (int c = 0; c < 6; ++c) {
                const float v = v6[c * P49_TILE + o];
                aR[c] = fmaf(v, (c < 3) ? cE : cH, aR[c]);
                aI[c] = fmaf(v, (c < 3) ? sE : sH, aI[c]);
            }
        }
        for (int c = 0; c < 6; ++c) {
            const size_t oo = (size_t)(c * nf + f) * cells + off;
            re[oo] += (double)aR[c]; im[oo] += (double)aI[c];
        }
    }
}
