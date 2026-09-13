// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Mingdian Liu
// update_h fused with the H-family absorber decay.
// Applicability, decided at runtime on the solver side: there is an absorbing layer **and** no
// source injection on the H side, so nothing sits between update_h and absorb(H) and the fusion
// cannot miss the source's share of the decay.
#define IDX(i, j, k) (((i) * ny + (j)) * nz + (k))

__device__ __forceinline__ float stretch(
    float* __restrict__ psi, const int id, const float d,
    const float a, const float b, const float invk, const int on)
{
    // on: 0 = no PML on this axis, 1 = normal, 2 = take the identity-region fast path
    if (!on) return fmaf(invk, d, 0.0f);
    // In the identity region (a==0 and b==1) psi is always 0, so that read and write are wasted.
    // Only enable it when the hit rate is high (on==2): at a low hit rate this branch blocks the
    // early issue of the psi load and measured 6% slower. on is a scalar parameter, uniform across
    // threads, so the branch itself is nearly free.
    if (on == 2 && a == 0.0f && b == 1.0f) return fmaf(invk, d, 0.0f);
    const float p = b * psi[id] + a * d;
    psi[id] = p;
    return invk * d + p;
}

extern "C" __global__ void update_h_absorb(
    float* __restrict__ Hx, float* __restrict__ Hy, float* __restrict__ Hz,
    const float* __restrict__ Ex, const float* __restrict__ Ey,
    const float* __restrict__ Ez,
    float* __restrict__ psi_hx_y, float* __restrict__ psi_hx_z,
    float* __restrict__ psi_hy_z, float* __restrict__ psi_hy_x,
    float* __restrict__ psi_hz_x, float* __restrict__ psi_hz_y,
    const int on_x, const int on_y, const int on_z,     // whether this derivative axis has a PML
    const float* __restrict__ idx_p, const float* __restrict__ idy_p,
    const float* __restrict__ idz_p,                    // 1 / primal spacing
    const float* __restrict__ ax, const float* __restrict__ bx,
    const float* __restrict__ kx,
    const float* __restrict__ ay, const float* __restrict__ by,
    const float* __restrict__ ky,
    const float* __restrict__ az, const float* __restrict__ bz,
    const float* __restrict__ kz,                       // CPML coefficients at H positions
    const int* __restrict__ nxt_x, const int* __restrict__ nxt_y,
    const int* __restrict__ nxt_z,
    const float* __restrict__ mnx_x, const float* __restrict__ mnx_y,
    const float* __restrict__ mnx_z,
    const float ch,                                     // dt / mu0
    // ---- absorber decay for the H family, fused in place ----
    const float* const* __restrict__ dec_int,
    const float* const* __restrict__ dec_half,
    const int* __restrict__ axis_a,
    const int* __restrict__ i0_a, const int* __restrict__ j0_a,
    const int* __restrict__ k0_a,
    const int* __restrict__ bi_a, const int* __restrict__ bj_a,
    const int* __restrict__ bk_a,
    const int ns,
    const int nx, const int ny, const int nz)
{
    const int k = blockIdx.x * blockDim.x + threadIdx.x;
    const int j = blockIdx.y * blockDim.y + threadIdx.y;
    const int i = blockIdx.z;
    if (i >= nx || j >= ny || k >= nz) return;

    const int id = IDX(i, j, k);
    const int ip = nxt_x[i], jp = nxt_y[j], kp = nxt_z[k];
    const float mi = mnx_x[i], mj = mnx_y[j], mk = mnx_z[k];
    const float ex = Ex[id], ey = Ey[id], ez = Ez[id];

    // mi/mj/mk each apply to the E components tangential to the high wall of that axis:
    //   x wall -> Ey, Ez     y wall -> Ex, Ez     z wall -> Ex, Ey
    // which are exactly the components the forward differences below read.

    // Hx: mu0 dHx/dt = -(dEz/dy - dEy/dz)
        const float dzy_x = (Ez[IDX(i, jp, k)] * mj - ez) * idy_p[j];
        const float dyz_x = (Ey[IDX(i, j, kp)] * mk - ey) * idz_p[k];
        const float sy_x = stretch(psi_hx_y, id, dzy_x, ay[j], by[j], ky[j], on_y);
        const float sz_x = stretch(psi_hx_z, id, dyz_x, az[k], bz[k], kz[k], on_z);
        const float vx = Hx[id] - ch * (sy_x - sz_x);
    // Hy: mu0 dHy/dt = -(dEx/dz - dEz/dx)
        const float dxz_y = (Ex[IDX(i, j, kp)] * mk - ex) * idz_p[k];
        const float dzx_y = (Ez[IDX(ip, j, k)] * mi - ez) * idx_p[i];
        const float sz_y = stretch(psi_hy_z, id, dxz_y, az[k], bz[k], kz[k], on_z);
        const float sx_y = stretch(psi_hy_x, id, dzx_y, ax[i], bx[i], kx[i], on_x);
        const float vy = Hy[id] - ch * (sz_y - sx_y);
    // Hz: mu0 dHz/dt = -(dEy/dx - dEx/dy)
        const float dyx_z = (Ey[IDX(ip, j, k)] * mi - ey) * idx_p[i];
        const float dxy_z = (Ex[IDX(i, jp, k)] * mj - ex) * idy_p[j];
        const float sx_z = stretch(psi_hz_x, id, dyx_z, ax[i], bx[i], kx[i], on_x);
        const float sy_z = stretch(psi_hz_y, id, dxy_z, ay[j], by[j], ky[j], on_y);
        const float vz = Hz[id] - ch * (sx_z - sy_z);

    // H-family absorber decay, multiplied onto the field **one factor at a time**, in the same
    // association order as absorb_batch's ((v*d1)*d2)*d3. Multiplying the factors together first
    // and then the field differs by 1 ulp, measured at 3.4e-08.
    float hx_ = vx, hy_ = vy, hz_ = vz;
    for (int s = 0; s < ns; ++s) {
        const int i0 = i0_a[s], j0 = j0_a[s], k0 = k0_a[s];
        if (i < i0 || i >= i0 + bi_a[s] || j < j0 || j >= j0 + bj_a[s]
            || k < k0 || k >= k0 + bk_a[s]) continue;
        const int ax = axis_a[s];
        const int off = (ax == 0) ? (i - i0) : (ax == 1) ? (j - j0) : (k - k0);
        const float di = dec_int[s][off], dh = dec_half[s][off];
        // is_e=0: half = (comp != axis)
        hx_ *= (ax != 0) ? dh : di;
        hy_ *= (ax != 1) ? dh : di;
        hz_ *= (ax != 2) ? dh : di;
    }
    Hx[id] = hx_;  Hy[id] = hy_;  Hz[id] = hz_;
}
