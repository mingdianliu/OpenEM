// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Mingdian Liu
// Batched time-monitor sampling: one launch handles every monitor.
//
// The problem: one scene had 13 time monitors times 2 phases, i.e. 26 launches per step. A probe
// measured 530 us of which only 86 us was real work; the other 444 us was launch overhead. Five
// point monitors of 8 to 12 cells each accounted for 14 launches and 220 us of that.
//
// This kernel packs the parameters into device arrays, one row per monitor. blockIdx.y selects the
// monitor and blockIdx.x * blockDim.x flattens the cells. What each thread computes is **word for
// word** what sample_time_box computes: the same slot guard, the same w * F[comp][id], the same
// add semantics.

extern "C" __global__ void sample_time_batch(
    float* const* __restrict__ outs,        // [nmon] output buffer per monitor
    const float* __restrict__ Ex, const float* __restrict__ Ey,
    const float* __restrict__ Ez, const float* __restrict__ Hx,
    const float* __restrict__ Hy, const float* __restrict__ Hz,
    const int* const* __restrict__ comps_p, // [nmon] component table per monitor
    const float* const* __restrict__ w_p,   // [nmon] weights per monitor, for this phase
    const int* __restrict__ nc_a,           // [nmon]
    const int* __restrict__ step_p,
    const int* __restrict__ beg_a, const int* __restrict__ iv_a,
    const int* __restrict__ end_a, const int* __restrict__ nslots_a,
    const int add,
    const int* __restrict__ i0_a, const int* __restrict__ j0_a,
    const int* __restrict__ k0_a,
    const int* __restrict__ ni_a, const int* __restrict__ nj_a,
    const int* __restrict__ nk_a,
    const int* __restrict__ cells_a,        // [nmon] ni*nj*nk, precomputed
    const int nmon, const int nx, const int ny, const int nz)
{
    const int mi = blockIdx.y;
    if (mi >= nmon) return;

    // Slot guard, word for word the same as sample_time_box
    const int m_ref = (*step_p + 1) - add;
    const int soff = m_ref - beg_a[mi];
    if (soff < 0 || m_ref >= end_a[mi] || soff % iv_a[mi]) return;
    const int slot = soff / iv_a[mi];
    if (slot >= nslots_a[mi]) return;

    const int cells = cells_a[mi];
    const int off = blockIdx.x * blockDim.x + threadIdx.x;
    if (off >= cells) return;

    // Unflatten off into (a, b, cc), the same (a*nj + b)*nk + cc convention as the single-launch version
    const int nj = nj_a[mi], nk = nk_a[mi];
    const int a = off / (nj * nk);
    const int rem = off - a * nj * nk;
    const int b = rem / nk;
    const int cc = rem - b * nk;

    const int id = ((i0_a[mi] + a) * ny + (j0_a[mi] + b)) * nz + (k0_a[mi] + cc);
    const float* F[6] = {Ex, Ey, Ez, Hx, Hy, Hz};
    const int nc = nc_a[mi];
    const int* comps = comps_p[mi];
    const float* weights = w_p[mi];
    float* out = outs[mi];
    const int nslots = nslots_a[mi];
    for (int c = 0; c < nc; ++c) {
        const float w = weights[c];
        if (w == 0.0f) continue;
        const size_t o = ((size_t)c * nslots + slot) * cells + off;
        const float v = w * F[comps[c]][id];
        out[o] = add ? out[o] + v : v;
    }
}
