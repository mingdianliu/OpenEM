// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Mingdian Liu
// Splitting the domain into an interior and a boundary shell.
// This file does not change a single byte of yee.cu: the old path is bit-frozen, and only lean
// mode uses the two kernels here.
//
// update_e_ade_lean is the interior lean kernel. Every piece of identity machinery is dropped:
//   psi is always 0 so it is neither read nor written, PEC and mpv are times 1.0f so they are not
//   multiplied, and the neighbour is simply i-1. Dispersive entries inline post and pre (the
//   deferral identity, on the same instruction sequence as dispersion_step). A non-dispersive
//   scene passes an empty eb box. The interior box must start 32-aligned in z; breaking that
//   alignment measured 28% slower.
// update_e_lut_box is a mechanical transformation of update_e_lut that takes a box origin, used to
//   run the six non-overlapping boxes of the boundary shell.

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

// The natively fused lean kernel: it runs only over the domain interior, leaving the skin (PML, the
// absorber shell, the domain faces) to the full kernel.
// Every identity value in the interior is skipped: psi is neither read nor written (it is always 0,
// and fmaf(invk,d,0) is bitwise the same as invk*d+0), PEC is not multiplied (times 1.0f is a
// bitwise identity), mpv is not multiplied (times 1.0f), and prv is simply i-1 since the interior
// never wraps. Curl and ADE are done in one pass. The launch covers only the interior box [o, e).
extern "C" __global__ void update_e_ade_lean(
    float* __restrict__ Ex, float* __restrict__ Ey, float* __restrict__ Ez,
    const float* __restrict__ Hx, const float* __restrict__ Hy,
    const float* __restrict__ Hz,
    const unsigned short* __restrict__ qx,
    const unsigned short* __restrict__ qy,
    const unsigned short* __restrict__ qz,
    const float* __restrict__ lca_x, const float* __restrict__ lcb_x,
    const float* __restrict__ lca_y, const float* __restrict__ lcb_y,
    const float* __restrict__ lca_z, const float* __restrict__ lcb_z,
    const float* __restrict__ idx_d, const float* __restrict__ idy_d,
    const float* __restrict__ idz_d,
    const float* __restrict__ kx, const float* __restrict__ ky,
    const float* __restrict__ kz,          // invk; always 1 in the interior, and fmaf keeps it bitwise
    float* __restrict__ p_re, float* __restrict__ p_im,
    const float* __restrict__ lut_am1_re, const float* __restrict__ lut_am1_im,
    const float* __restrict__ lut_b_re, const float* __restrict__ lut_b_im,
    const int* __restrict__ cidx,
    const unsigned int* __restrict__ estart_x,
    const unsigned int* __restrict__ estart_y,
    const unsigned int* __restrict__ estart_z,
    const unsigned char* __restrict__ ecount_x,
    const unsigned char* __restrict__ ecount_y,
    const unsigned char* __restrict__ ecount_z,
    const int eb_x0, const int eb_x1, const int eb_y0, const int eb_y1,
    const int eb_z0, const int eb_z1,
    const float cw,
    const int ox, const int oy, const int oz,     // launch origin (z already 32-aligned)
    const int oz_real,                            // true interior z start, used as a guard
    const int ex_, const int ey_, const int ez_,  // interior box end, exclusive
    const int nx, const int ny, const int nz)
{
    const int k = blockIdx.x * blockDim.x + threadIdx.x + oz;
    const int j = blockIdx.y * blockDim.y + threadIdx.y + oy;
    const int i = blockIdx.z + ox;
    if (i >= ex_ || j >= ey_ || k >= ez_ || k < oz_real) return;

    const int id = IDX(i, j, k);
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
        const float dzy = (hz - Hz[IDX(i, j - 1, k)]) * idy_d[j];
        const float dyz = (hy - Hy[IDX(i, j, k - 1)]) * idz_d[k];
        const float sy = fmaf(ky[j], dzy, 0.0f);
        const float sz = fmaf(kz[k], dyz, 0.0f);
        const int ix_ = (int)qx[id];
        Ex[id] = lca_x[ix_] * en + lcb_x[ix_] * (sy - sz - cw * w);
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
        const float dxz = (hx - Hx[IDX(i, j, k - 1)]) * idz_d[k];
        const float dzx = (hz - Hz[IDX(i - 1, j, k)]) * idx_d[i];
        const float sz = fmaf(kz[k], dxz, 0.0f);
        const float sx = fmaf(kx[i], dzx, 0.0f);
        const int iy_ = (int)qy[id];
        Ey[id] = lca_y[iy_] * en + lcb_y[iy_] * (sz - sx - cw * w);
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
        const float dyx = (hy - Hy[IDX(i - 1, j, k)]) * idx_d[i];
        const float dxy = (hx - Hx[IDX(i, j - 1, k)]) * idy_d[j];
        const float sx = fmaf(kx[i], dyx, 0.0f);
        const float sy = fmaf(ky[j], dxy, 0.0f);
        const int iz_ = (int)qz[id];
        Ez[id] = lca_z[iz_] * en + lcb_z[iz_] * (sx - sy - cw * w);
    }
}



// Interior lean kernel on the H side: psi is neither read nor written (fmaf identity), mnx is not
// multiplied (times 1.0f), and the neighbour is simply i+1. The launch origin is already 32-aligned
// in z, with a k < oz_real guard, as on the E side.
extern "C" __global__ void update_h_lean(
    float* __restrict__ Hx, float* __restrict__ Hy, float* __restrict__ Hz,
    const float* __restrict__ Ex, const float* __restrict__ Ey,
    const float* __restrict__ Ez,
    const float* __restrict__ idx_p, const float* __restrict__ idy_p,
    const float* __restrict__ idz_p,
    const float* __restrict__ kx, const float* __restrict__ ky,
    const float* __restrict__ kz,
    const float ch,
    const int ox, const int oy, const int oz,
    const int oz_real,
    const int ex_, const int ey_, const int ez_,
    const int nx, const int ny, const int nz)
{
    const int k = blockIdx.x * blockDim.x + threadIdx.x + oz;
    const int j = blockIdx.y * blockDim.y + threadIdx.y + oy;
    const int i = blockIdx.z + ox;
    if (i >= ex_ || j >= ey_ || k >= ez_ || k < oz_real) return;

    const int id = IDX(i, j, k);
    const float ex = Ex[id], ey = Ey[id], ez = Ez[id];
    // Hx
    {
        const float dzy = (Ez[IDX(i, j + 1, k)] - ez) * idy_p[j];
        const float dyz = (Ey[IDX(i, j, k + 1)] - ey) * idz_p[k];
        const float sy = fmaf(ky[j], dzy, 0.0f);
        const float sz = fmaf(kz[k], dyz, 0.0f);
        Hx[id] -= ch * (sy - sz);
    }
    // Hy
    {
        const float dxz = (Ex[IDX(i, j, k + 1)] - ex) * idz_p[k];
        const float dzx = (Ez[IDX(i + 1, j, k)] - ez) * idx_p[i];
        const float sz = fmaf(kz[k], dxz, 0.0f);
        const float sx = fmaf(kx[i], dzx, 0.0f);
        Hy[id] -= ch * (sz - sx);
    }
    // Hz
    {
        const float dyx = (Ey[IDX(i + 1, j, k)] - ey) * idx_p[i];
        const float dxy = (Ex[IDX(i, j + 1, k)] - ex) * idy_p[j];
        const float sx = fmaf(kx[i], dyx, 0.0f);
        const float sy = fmaf(ky[j], dxy, 0.0f);
        Hz[id] -= ch * (sx - sy);
    }
}


// Single boundary-shell kernel: launch over the whole domain and return early in the interior,
// replacing six fragmented box launches.
extern "C" __global__ void update_e_lut_shell(
    float* __restrict__ Ex, float* __restrict__ Ey, float* __restrict__ Ez,
    const float* __restrict__ Hx, const float* __restrict__ Hy,
    const float* __restrict__ Hz,
    float* __restrict__ psi_ex_y, float* __restrict__ psi_ex_z,
    float* __restrict__ psi_ey_z, float* __restrict__ psi_ey_x,
    float* __restrict__ psi_ez_x, float* __restrict__ psi_ez_y,
    const int on_x, const int on_y, const int on_z,     // whether this derivative axis has a PML
    const unsigned short* __restrict__ qx,              // uint16 index into the
    const unsigned short* __restrict__ qy,              //   (ca, cb) pair table
    const unsigned short* __restrict__ qz,
    const float* __restrict__ lca_x, const float* __restrict__ lcb_x,
    const float* __restrict__ lca_y, const float* __restrict__ lcb_y,
    const float* __restrict__ lca_z, const float* __restrict__ lcb_z,
    const float* __restrict__ idx_d, const float* __restrict__ idy_d,
    const float* __restrict__ idz_d,                    // 1 / dual spacing
    const float* __restrict__ ax, const float* __restrict__ bx,
    const float* __restrict__ kx,
    const float* __restrict__ ay, const float* __restrict__ by,
    const float* __restrict__ ky,
    const float* __restrict__ az, const float* __restrict__ bz,
    const float* __restrict__ kz,                       // CPML coefficients at E positions
    const int* __restrict__ prv_x, const int* __restrict__ prv_y,
    const int* __restrict__ prv_z,
    const float* __restrict__ mpv_x, const float* __restrict__ mpv_y,
    const float* __restrict__ mpv_z,
    const float* __restrict__ pec_x, const float* __restrict__ pec_y,
    const float* __restrict__ pec_z,
    const float* __restrict__ hist_x, const float* __restrict__ hist_y,
    const float* __restrict__ hist_z,                   // Σ 2Re(P_half - P^n)
    const float cw,                                     // eps0 / dt
    const int hb_x0, const int hb_x1, const int hb_y0, const int hb_y1,
    const int hb_z0, const int hb_z1,                   // bounding box of the dispersion hist
    const int ox, const int oy, const int oz,
    const int ex_, const int ey_, const int ez_,
    const int zseg0, const int zseg1,   // z range of this segment, [zseg0, zseg1)
    const int nx, const int ny, const int nz)
{
    const int k = blockIdx.x * blockDim.x + threadIdx.x + zseg0;
    const int j = blockIdx.y * blockDim.y + threadIdx.y;
    const int i = blockIdx.z;
    if (i >= nx || j >= ny || k >= zseg1 || k >= nz) return;
    // Interior cells belong to the lean kernel, so skip them here: the shell is the whole domain
    // minus the interior
    if (i >= ox && i < ex_ && j >= oy && j < ey_
        && k >= oz && k < ez_) return;

    const int id = IDX(i, j, k);
    const int im = prv_x[i], jm = prv_y[j], km = prv_z[k];
    const float mi = mpv_x[i], mj = mpv_y[j], mk = mpv_z[k];
    const float hx = Hx[id], hy = Hy[id], hz = Hz[id];

    // Outside the dispersion bounding box hist is always 0: dispersion writes only entry cells and
    // this kernel reads only its own cell id. So outside the box a constant 0.0f replaces three
    // loads and is bitwise identical to the +0.0f that would have been read. A non-dispersive scene
    // passes an empty box (all zeros) and skips the loads everywhere, saving the traffic of an
    // all-zero array.
    const bool hin = (i >= hb_x0) & (i < hb_x1) & (j >= hb_y0) & (j < hb_y1)
                   & (k >= hb_z0) & (k < hb_z1);
    const float ph_x = hin ? hist_x[id] : 0.0f;
    const float ph_y = hin ? hist_y[id] : 0.0f;
    const float ph_z = hin ? hist_z[id] : 0.0f;

    // pec applies to the components tangential to the low wall of that axis, i.e. the two other
    // than its own:
    //   Ex -> pec_y * pec_z     Ey -> pec_z * pec_x     Ez -> pec_x * pec_y
    // The normal component is not constrained by PEC, so Ex is not multiplied by pec_x.

    // Ex: eps dEx/dt = dHz/dy - dHy/dz
    {
        const float dzy = (hz - Hz[IDX(i, jm, k)] * mj) * idy_d[j];
        const float dyz = (hy - Hy[IDX(i, j, km)] * mk) * idz_d[k];
        const float sy = stretch(psi_ex_y, id, dzy, ay[j], by[j], ky[j], on_y);
        const float sz = stretch(psi_ex_z, id, dyz, az[k], bz[k], kz[k], on_z);
        const int ix_ = (int)qx[id];
        Ex[id] = (lca_x[ix_] * Ex[id]
                  + lcb_x[ix_] * (sy - sz - cw * ph_x)) * pec_y[j] * pec_z[k];
    }
    // Ey: eps dEy/dt = dHx/dz - dHz/dx
    {
        const float dxz = (hx - Hx[IDX(i, j, km)] * mk) * idz_d[k];
        const float dzx = (hz - Hz[IDX(im, j, k)] * mi) * idx_d[i];
        const float sz = stretch(psi_ey_z, id, dxz, az[k], bz[k], kz[k], on_z);
        const float sx = stretch(psi_ey_x, id, dzx, ax[i], bx[i], kx[i], on_x);
        const int iy_ = (int)qy[id];
        Ey[id] = (lca_y[iy_] * Ey[id]
                  + lcb_y[iy_] * (sz - sx - cw * ph_y)) * pec_z[k] * pec_x[i];
    }
    // Ez: eps dEz/dt = dHy/dx - dHx/dy
    {
        const float dyx = (hy - Hy[IDX(im, j, k)] * mi) * idx_d[i];
        const float dxy = (hx - Hx[IDX(i, jm, k)] * mj) * idy_d[j];
        const float sx = stretch(psi_ez_x, id, dyx, ax[i], bx[i], kx[i], on_x);
        const float sy = stretch(psi_ez_y, id, dxy, ay[j], by[j], ky[j], on_y);
        const int iz_ = (int)qz[id];
        Ez[id] = (lca_z[iz_] * Ez[id]
                  + lcb_z[iz_] * (sx - sy - cw * ph_z)) * pec_x[i] * pec_y[j];
    }
}

extern "C" __global__ void update_h_shell(
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
    const int ox, const int oy, const int oz,
    const int ex_, const int ey_, const int ez_,
    const int zseg0, const int zseg1,   // z range of this segment, [zseg0, zseg1)
    const int nx, const int ny, const int nz)
{
    const int k = blockIdx.x * blockDim.x + threadIdx.x + zseg0;
    const int j = blockIdx.y * blockDim.y + threadIdx.y;
    const int i = blockIdx.z;
    if (i >= nx || j >= ny || k >= zseg1 || k >= nz) return;
    // Interior cells belong to the lean kernel, so skip them here: the shell is the whole domain
    // minus the interior
    if (i >= ox && i < ex_ && j >= oy && j < ey_
        && k >= oz && k < ez_) return;

    const int id = IDX(i, j, k);
    const int ip = nxt_x[i], jp = nxt_y[j], kp = nxt_z[k];
    const float mi = mnx_x[i], mj = mnx_y[j], mk = mnx_z[k];
    const float ex = Ex[id], ey = Ey[id], ez = Ez[id];

    // mi/mj/mk each apply to the E components tangential to the high wall of that axis:
    //   x wall -> Ey, Ez     y wall -> Ex, Ez     z wall -> Ex, Ey
    // which are exactly the components the forward differences below read.

    // Hx: mu0 dHx/dt = -(dEz/dy - dEy/dz)
    {
        const float dzy = (Ez[IDX(i, jp, k)] * mj - ez) * idy_p[j];
        const float dyz = (Ey[IDX(i, j, kp)] * mk - ey) * idz_p[k];
        const float sy = stretch(psi_hx_y, id, dzy, ay[j], by[j], ky[j], on_y);
        const float sz = stretch(psi_hx_z, id, dyz, az[k], bz[k], kz[k], on_z);
        Hx[id] -= ch * (sy - sz);
    }
    // Hy: mu0 dHy/dt = -(dEx/dz - dEz/dx)
    {
        const float dxz = (Ex[IDX(i, j, kp)] * mk - ex) * idz_p[k];
        const float dzx = (Ez[IDX(ip, j, k)] * mi - ez) * idx_p[i];
        const float sz = stretch(psi_hy_z, id, dxz, az[k], bz[k], kz[k], on_z);
        const float sx = stretch(psi_hy_x, id, dzx, ax[i], bx[i], kx[i], on_x);
        Hy[id] -= ch * (sz - sx);
    }
    // Hz: mu0 dHz/dt = -(dEy/dx - dEx/dy)
    {
        const float dyx = (Ey[IDX(ip, j, k)] * mi - ey) * idx_p[i];
        const float dxy = (Ex[IDX(i, jp, k)] * mj - ex) * idy_p[j];
        const float sx = stretch(psi_hz_x, id, dyx, ax[i], bx[i], kx[i], on_x);
        const float sy = stretch(psi_hz_y, id, dxy, ay[j], by[j], ky[j], on_y);
        Hz[id] -= ch * (sx - sy);
    }
}
