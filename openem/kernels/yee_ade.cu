// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Mingdian Liu
// Fully fused ADE: the pole update is inlined into update_e and the hist array disappears entirely.
//
// The structure is update_e_lut with the hist read removed and replaced by an inlined post(n) plus
// pre(n+1) on that cell's own poles, using the same instruction sequence as dispersion_step. The
// justification is a deferral identity: the E that dispersion_step reads at the end of step n is
// the same memory that update_e reads at step n+1, with nothing writing it in between. That holds
// when the cell is outside any absorbing layer and there is no second tensor pass. So deferring it
// into update_e reads a bitwise identical en, and w held in a register equals the value that would
// have been written to hist and read back. The first step is free: E^0 = 0 gives w = 0, which is
// the initial hist.
//
// Entry layout: at solver initialization the entries are reordered into canonical (comp, cell)
// order. That is a pure permutation which leaves the pole order within each entry untouched, so
// the P values stay bitwise identical. estart/ecount then index a slice per cell directly, and the
// comp, cell and ofs streams disappear.
//
// Applicability, decided on the solver side: the LUT hits, there is no modulation and no tensor.
// Entries inside an absorbing layer stay on the original dispersion_step_dec, which is a thin
// shell; entries outside it come here.

#define IDX(i, j, k) (((i) * ny + (j)) * nz + (k))

__device__ __forceinline__ float stretch(
    float* __restrict__ psi, const int id, const float d,
    const float a, const float b, const float invk, const int on)
{
    if (!on) return fmaf(invk, d, 0.0f);
    const float p = b * psi[id] + a * d;
    psi[id] = p;
    return invk * d + p;
}

// Inlined pole update for this cell and component: post followed by pre, returning w = sum 2*Re(dP).
__device__ __forceinline__ float ade_w(
    float* __restrict__ p_re, float* __restrict__ p_im,
    const float* __restrict__ lut_am1_re, const float* __restrict__ lut_am1_im,
    const float* __restrict__ lut_b_re, const float* __restrict__ lut_b_im,
    const int* __restrict__ cidx,
    const unsigned int s, const int n, const float en)
{
    float w = 0.0f;
    for (int t = 0; t < n; ++t) {
        const int ci = cidx[s + t];
        const float br = lut_b_re[ci], bi = lut_b_im[ci];
        const float pr = p_re[s + t] + br * en;
        const float pi = p_im[s + t] + bi * en;
        const float ar = lut_am1_re[ci], ai = lut_am1_im[ci];
        const float dr = ar * pr - ai * pi + br * en;
        const float di = ar * pi + ai * pr + bi * en;
        p_re[s + t] = pr + dr;
        p_im[s + t] = pi + di;
        w += 2.0f * dr;
    }
    return w;
}

extern "C" __global__ void update_e_ade(
    float* __restrict__ Ex, float* __restrict__ Ey, float* __restrict__ Ez,
    const float* __restrict__ Hx, const float* __restrict__ Hy,
    const float* __restrict__ Hz,
    float* __restrict__ psi_ex_y, float* __restrict__ psi_ex_z,
    float* __restrict__ psi_ey_z, float* __restrict__ psi_ey_x,
    float* __restrict__ psi_ez_x, float* __restrict__ psi_ez_y,
    const int on_x, const int on_y, const int on_z,
    const unsigned short* __restrict__ qx,
    const unsigned short* __restrict__ qy,
    const unsigned short* __restrict__ qz,
    const float* __restrict__ lca_x, const float* __restrict__ lcb_x,
    const float* __restrict__ lca_y, const float* __restrict__ lcb_y,
    const float* __restrict__ lca_z, const float* __restrict__ lcb_z,
    const float* __restrict__ idx_d, const float* __restrict__ idy_d,
    const float* __restrict__ idz_d,
    const float* __restrict__ ax, const float* __restrict__ bx,
    const float* __restrict__ kx,
    const float* __restrict__ ay, const float* __restrict__ by,
    const float* __restrict__ ky,
    const float* __restrict__ az, const float* __restrict__ bz,
    const float* __restrict__ kz,
    const int* __restrict__ prv_x, const int* __restrict__ prv_y,
    const int* __restrict__ prv_z,
    const float* __restrict__ mpv_x, const float* __restrict__ mpv_y,
    const float* __restrict__ mpv_z,
    const float* __restrict__ pec_x, const float* __restrict__ pec_y,
    const float* __restrict__ pec_z,
    // ---- fused ADE section, replacing three hist reads plus the bounding box ----
    float* __restrict__ p_re, float* __restrict__ p_im,
    const float* __restrict__ lut_am1_re, const float* __restrict__ lut_am1_im,
    const float* __restrict__ lut_b_re, const float* __restrict__ lut_b_im,
    const int* __restrict__ cidx,
    const unsigned int* __restrict__ estart_x,   // (ncell,) start of this cell's pole slice
    const unsigned int* __restrict__ estart_y,
    const unsigned int* __restrict__ estart_z,
    const unsigned char* __restrict__ ecount_x,  // (ncell,) number of poles; 0 means no entry
    const unsigned char* __restrict__ ecount_y,
    const unsigned char* __restrict__ ecount_z,
    const int eb_x0, const int eb_x1, const int eb_y0, const int eb_y1,
    const int eb_z0, const int eb_z1,            // bounding box of dispersive entries outside the layer
    const float cw,
    const int nx, const int ny, const int nz)
{
    const int k = blockIdx.x * blockDim.x + threadIdx.x;
    const int j = blockIdx.y * blockDim.y + threadIdx.y;
    const int i = blockIdx.z;
    if (i >= nx || j >= ny || k >= nz) return;

    const int id = IDX(i, j, k);
    const int im = prv_x[i], jm = prv_y[j], km = prv_z[k];
    const float mi = mpv_x[i], mj = mpv_y[j], mk = mpv_z[k];
    const float hx = Hx[id], hy = Hy[id], hz = Hz[id];

    const bool ein = (i >= eb_x0) & (i < eb_x1) & (j >= eb_y0) & (j < eb_y1)
                   & (k >= eb_z0) & (k < eb_z1);

    // Ex
    {
        const float en = Ex[id];
        float w = 0.0f;
        if (ein) {
            const int n = ecount_x[id];
            if (n) w = ade_w(p_re, p_im, lut_am1_re, lut_am1_im,
                             lut_b_re, lut_b_im, cidx, estart_x[id], n, en);
        }
        const float dzy = (hz - Hz[IDX(i, jm, k)] * mj) * idy_d[j];
        const float dyz = (hy - Hy[IDX(i, j, km)] * mk) * idz_d[k];
        const float sy = stretch(psi_ex_y, id, dzy, ay[j], by[j], ky[j], on_y);
        const float sz = stretch(psi_ex_z, id, dyz, az[k], bz[k], kz[k], on_z);
        const int ix_ = (int)qx[id];
        Ex[id] = (lca_x[ix_] * en
                  + lcb_x[ix_] * (sy - sz - cw * w)) * pec_y[j] * pec_z[k];
    }
    // Ey
    {
        const float en = Ey[id];
        float w = 0.0f;
        if (ein) {
            const int n = ecount_y[id];
            if (n) w = ade_w(p_re, p_im, lut_am1_re, lut_am1_im,
                             lut_b_re, lut_b_im, cidx, estart_y[id], n, en);
        }
        const float dxz = (hx - Hx[IDX(i, j, km)] * mk) * idz_d[k];
        const float dzx = (hz - Hz[IDX(im, j, k)] * mi) * idx_d[i];
        const float sz = stretch(psi_ey_z, id, dxz, az[k], bz[k], kz[k], on_z);
        const float sx = stretch(psi_ey_x, id, dzx, ax[i], bx[i], kx[i], on_x);
        const int iy_ = (int)qy[id];
        Ey[id] = (lca_y[iy_] * en
                  + lcb_y[iy_] * (sz - sx - cw * w)) * pec_z[k] * pec_x[i];
    }
    // Ez
    {
        const float en = Ez[id];
        float w = 0.0f;
        if (ein) {
            const int n = ecount_z[id];
            if (n) w = ade_w(p_re, p_im, lut_am1_re, lut_am1_im,
                             lut_b_re, lut_b_im, cidx, estart_z[id], n, en);
        }
        const float dyx = (hy - Hy[IDX(im, j, k)] * mi) * idx_d[i];
        const float dxy = (hx - Hx[IDX(i, jm, k)] * mj) * idy_d[j];
        const float sx = stretch(psi_ez_x, id, dyx, ax[i], bx[i], kx[i], on_x);
        const float sy = stretch(psi_ez_y, id, dxy, ay[j], by[j], ky[j], on_y);
        const int iz_ = (int)qz[id];
        Ez[id] = (lca_z[iz_] * en
                  + lcb_z[iz_] * (sx - sy - cw * w)) * pec_x[i] * pec_y[j];
    }
}

// The natively fused lean kernel: it runs only over the domain interior, leaving the skin (PML, the
// absorber shell, the domain faces) to the full kernel.
// Every identity value in the interior is skipped: psi is neither read nor written (it is always 0,
// and fmaf(invk,d,0) is bitwise the same as invk*d+0), PEC is not multiplied (times 1.0f is a
// bitwise identity), mpv is not multiplied (times 1.0f), and prv is simply i-1 since the interior
// never wraps. Curl and ADE are done in one pass. The launch covers only the interior box [o, e).
