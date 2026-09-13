// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Mingdian Liu
// Batched absorber decay: one launch covers the union of every slab.
//
// The problem: one launch per slab (10 per step on one scene, 8 on another), and together those
// slabs cover nearly the whole domain (16.2M of 16.8M cells on one of them). Measured at 521 us to
// move 389 MB, i.e. 0.75 TB/s, only 27% of roofline, from the poor shape of a thin box plus the
// tail effect of many launches.
//
// This kernel maps threads over the union box and loops inside over the slabs covering that cell,
// multiplying the decay factors **in the original order**. The multiplication order is exactly the
// same as with per-slab launches, so the result is bitwise identical. A slab's decay array is tiny
// (nlay floats) and stays in L1/L2.

extern "C" __global__ void absorb_batch(
    float* __restrict__ Fx, float* __restrict__ Fy, float* __restrict__ Fz,
    const float* const* __restrict__ dec_int,   // [ns] whole-cell profile per slab
    const float* const* __restrict__ dec_half,  // [ns] half-cell profile per slab
    const int* __restrict__ axis_a,             // [ns]
    const int* __restrict__ i0_a, const int* __restrict__ j0_a,
    const int* __restrict__ k0_a,
    const int* __restrict__ bi_a, const int* __restrict__ bj_a,
    const int* __restrict__ bk_a,
    const int is_e, const int ns,
    const int ux0, const int uy0, const int uz0,   // union box origin
    const int ux1, const int uy1, const int uz1,   // union box end, exclusive
    const int nx, const int ny, const int nz)
{
    const int k = blockIdx.x * blockDim.x + threadIdx.x + uz0;
    const int j = blockIdx.y * blockDim.y + threadIdx.y + uy0;
    const int i = blockIdx.z + ux0;
    if (i >= ux1 || j >= uy1 || k >= uz1) return;

    float fx = Fx[((i * ny) + j) * nz + k];
    float fy = Fy[((i * ny) + j) * nz + k];
    float fz = Fz[((i * ny) + j) * nz + k];
    bool hit = false;

    for (int s = 0; s < ns; ++s) {
        const int i0 = i0_a[s], j0 = j0_a[s], k0 = k0_a[s];
        const int bi = bi_a[s], bj = bj_a[s], bk = bk_a[s];
        if (i < i0 || i >= i0 + bi || j < j0 || j >= j0 + bj
            || k < k0 || k >= k0 + bk) continue;
        hit = true;
        const int ax = axis_a[s];
        const int off = (ax == 0) ? (i - i0) : (ax == 1) ? (j - j0) : (k - k0);
        const float di = dec_int[s][off];
        const float dh = dec_half[s][off];
        // Same convention as absorb_slab: the normal component takes the half-cell profile and the
        // tangential ones the whole-cell profile; is_e picks the staggering for the E or H family
        if (is_e) {
            if (ax == 0) { fx *= dh; fy *= di; fz *= di; }
            else if (ax == 1) { fx *= di; fy *= dh; fz *= di; }
            else { fx *= di; fy *= di; fz *= dh; }
        } else {
            if (ax == 0) { fx *= di; fy *= dh; fz *= dh; }
            else if (ax == 1) { fx *= dh; fy *= di; fz *= dh; }
            else { fx *= dh; fy *= dh; fz *= di; }
        }
    }
    if (!hit) return;
    Fx[((i * ny) + j) * nz + k] = fx;
    Fy[((i * ny) + j) * nz + k] = fy;
    Fz[((i * ny) + j) * nz + k] = fz;
}


// Version with 1D pre-multiplied profiles. Slabs on the same axis never intersect (the host checks
// this), so each cell has at most one factor per axis; the six profile values are the original
// single factors and the multiplication order follows the slab order, so this is bitwise identical.
// A cell whose six g values are all 1 returns before it ever reads F.
extern "C" __global__ void absorb_batch_1d(
    float* __restrict__ Fx, float* __restrict__ Fy, float* __restrict__ Fz,
    const float* __restrict__ gx_i, const float* __restrict__ gx_h,
    const float* __restrict__ gy_i, const float* __restrict__ gy_h,
    const float* __restrict__ gz_i, const float* __restrict__ gz_h,
    const int is_e,
    const int ux0, const int uy0, const int uz0,
    const int ux1, const int uy1, const int uz1,
    const int nx, const int ny, const int nz)
{
    const int k = blockIdx.x * blockDim.x + threadIdx.x + uz0;
    const int j = blockIdx.y * blockDim.y + threadIdx.y + uy0;
    const int i = blockIdx.z + ux0;
    if (i >= ux1 || j >= uy1 || k >= uz1) return;
    const float xi = gx_i[i], xh = gx_h[i];
    const float yi = gy_i[j], yh = gy_h[j];
    const float zi = gz_i[k], zh = gz_h[k];
    if (xi == 1.0f && xh == 1.0f && yi == 1.0f && yh == 1.0f
        && zi == 1.0f && zh == 1.0f) return;
    const long id = ((long)i * ny + j) * nz + k;
    float fx = Fx[id], fy = Fy[id], fz = Fz[id];
    if (is_e) {
        fx *= xh; fy *= xi; fz *= xi;
        fx *= yi; fy *= yh; fz *= yi;
        fx *= zi; fy *= zi; fz *= zh;
    } else {
        fx *= xi; fy *= xh; fz *= xh;
        fx *= yh; fy *= yi; fz *= yh;
        fx *= zh; fy *= zh; fz *= zi;
    }
    Fx[id] = fx; Fy[id] = fy; Fz[id] = fz;
}
