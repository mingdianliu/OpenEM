// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Mingdian Liu
// A near-to-far surface monitor accumulates only the **tangential** components.
//
// Basis (checked against Tidy3D source, field_projection/projector.py from line 730): the
// equivalent surface currents use only E1, E2, H1, H2, that is two tangential E and two tangential
// H. The normal components are never read, so a third of a six-component store is pure waste. On
// the largest scattering case this took re/im from 57 GB down to 38 GB and the total requirement
// from 85 GB to 66 GB, which fits in 80 GB.
//
// Layout: re/im are (ncomp, nf, cells) and snap is (T, ncomp, cells). cmap[slot], in [0,6), gives
// the physical component for that slot (0..2 = E, 3..5 = H). With ncomp=6 and cmap=[0..5] this is
// word for word equivalent to the original kernel.

extern "C" __global__ void dft_snap_box_sub(
    float* __restrict__ snap,           // (T, ncomp, cells)
    const float* __restrict__ Ex, const float* __restrict__ Ey,
    const float* __restrict__ Ez, const float* __restrict__ Hx,
    const float* __restrict__ Hy, const float* __restrict__ Hz,
    const int* __restrict__ stp, const int T,
    const int* __restrict__ cmap, const int ncomp,
    const int i0, const int j0, const int k0,
    const int ni, const int nj, const int nk,
    const int nx, const int ny, const int nz)
{
    const int tid = blockIdx.x * blockDim.x + threadIdx.x;
    const int a = blockIdx.z;
    if (a >= ni || tid >= nj * nk) return;
    const int b = tid / nk;
    const int cc = tid - b * nk;
    const int id = ((i0 + a) * ny + (j0 + b)) * nz + (k0 + cc);
    const size_t cells = (size_t)ni * nj * nk;
    const size_t off = (size_t)(a * nj + b) * nk + cc;
    const float* F[6] = {Ex, Ey, Ez, Hx, Hy, Hz};
    float* row = snap + (size_t)((*stp) % T) * ncomp * cells;
    for (int s = 0; s < ncomp; ++s)
        row[(size_t)s * cells + off] = F[cmap[s]][id];
}

extern "C" __global__ void dft_flush_box_pf_sub(
    double* __restrict__ re, double* __restrict__ im,
    const float* __restrict__ snap,
    const float* __restrict__ tab,
    const float* __restrict__ win_e, const float* __restrict__ win_h,
    const int* __restrict__ stp,
    const int T, const int force_cnt,
    const int f0, const int nf_all, const int tmod,
    const int cells, const int nf,
    const int* __restrict__ cmap, const int ncomp)
{
    // fp32 partial sums within the batch, added into the fp64 accumulator once at the batch end.
    const int off = blockIdx.x * blockDim.x + threadIdx.x;
    const int f = blockIdx.y;
    if (off >= cells || f >= nf) return;
    const int s_now = (force_cnt > 0) ? (*stp - 1) : (*stp);
    int cnt = force_cnt;
    if (cnt == 0) {
        if ((s_now + 1) % T != 0) return;
        cnt = T;
    }
    // ncomp <= 6, so accumulate in a register array as the original kernel does, avoiding a global
    // read and write every step
    float aR[6], aI[6];
    for (int c = 0; c < ncomp; ++c) { aR[c] = 0.f; aI[c] = 0.f; }
    for (int q = 0; q < cnt; ++q) {
        const int s = s_now - cnt + 1 + q;
        const float apod_e = win_e[s], apod_h = win_h[s];
        const float* ph = tab + (((size_t)(s % tmod) * nf_all) + f0 + f) * 4;
        const float cE = ph[0] * apod_e, sE = ph[1] * apod_e;
        const float cH = ph[2] * apod_h, sH = ph[3] * apod_h;
        const float* row = snap + (size_t)(s % T) * ncomp * cells;
        for (int c = 0; c < ncomp; ++c) {
            const float v = row[(size_t)c * cells + off];
            aR[c] = fmaf(v, (cmap[c] < 3) ? cE : cH, aR[c]);
            aI[c] = fmaf(v, (cmap[c] < 3) ? sE : sH, aI[c]);
        }
    }
    for (int c = 0; c < ncomp; ++c) {
        const size_t o = (size_t)(c * nf + f) * cells + off;
        re[o] += (double)aR[c]; im[o] += (double)aI[c];
    }
}
