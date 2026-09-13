// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Mingdian Liu
// General non-uniform Yee update with CPML on all three axes. Each axis is independently periodic
// or absorbing, and the medium may be lossy.
//
// **There is only one code path.** No topology variants, no lossy/lossless variants. The trick is
// pushed all the way: fill the coefficients with identity values.
//
//   * Periodic axis: fill that axis's CPML coefficients with (a=0, b=1, invk=1), so psi stays 0,
//     the derivative passes through unchanged and the PML term vanishes on its own.
//   * Inside vs outside the PML: one coefficient array of length n, filled with identity values
//     outside the PML too.
//   * Lossless medium: ca filled with 1.0, so the semi-implicit update degenerates to explicit.
//
// The kernel therefore contains no powf, no "is this cell inside the PML", no "is this axis
// periodic or absorbing" and no "is this medium lossy" branch. Boundary topology is expressed by
// tables precomputed on the Python side (Grid.Axis.index_tables):
//
//   nxt[i] / prv[i]     neighbour indices (periodic wraparound; clamped inside the domain at an
//                       absorbing end)
//   mnx[i] / mpv[i]     whether that neighbour really exists; an out-of-domain neighbour at an
//                       absorbing end is multiplied by 0
//   pec[i]              zeroes tangential E at i=0 of an absorbing axis (PEC closes off the PML)
//
// mpv does more than select a value: it makes the backward difference at i=0 of an absorbing end
// read a finite number, so psi cannot accumulate an Inf on cells that PEC is going to zero anyway,
// which would then give 0*Inf = NaN.
//
// Yee layout (measured from sim.grid.yee, aligned component by component with Tidy3D):
//   Ex(ctr_x, edge_y, edge_z)   Ey(edge_x, ctr_y, edge_z)   Ez(edge_x, edge_y, ctr_z)
//   Hx(edge_x, ctr_y, ctr_z)    Hy(ctr_x, edge_y, ctr_z)    Hz(ctr_x, ctr_y, edge_z)
//
// Therefore:
//   H update = forward difference / primal spacing, with CPML coefficients at H positions
//              (create_sfactor_f)
//   E update = backward difference / dual spacing, with CPML coefficients at E positions
//              (create_sfactor_b)
//
// Each field component has two derivative directions, so there are 12 psi arrays in all. The ones
// on periodic axes are identically zero yet still take part in the arithmetic and occupy memory;
// that is recorded as known debt, correctness first.
//
// Memory layout is C order (nx, ny, nz) with **z contiguous**. The thread x dimension maps to k so
// accesses coalesce.

#define IDX(i, j, k) (((i) * ny + (j)) * nz + (k))

// CPML psi recursion plus coordinate stretching. Outside a PML, where (a=0, b=1, invk=1), it
// returns d unchanged.
//
// ``on == 0`` means **this derivative axis has no PML anywhere** (a is identically 0, b identically
// 1): psi starts at 0 and stays there, and the recursion degenerates to invk*d. Skipping the psi
// read and write then saves 8 B per cell per call. Across 12 psi arrays that is 96 B per cell per
// step, and a scene whose six faces are all absorbers with no PML face at all was spending every
// one of those bytes for nothing.
//
// ``fmaf(invk, d, 0.0f)`` is bitwise identical to the original ``invk*d + psi`` with psi=0, down to
// the sign of a zero: the compiler emits exactly that FMA for the original expression.
// ``on`` has one value for the whole launch, so the branch is warp-uniform and never diverges.
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

extern "C" __global__ void update_h(
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

// E update. Semi-implicit in the conductivity, and carrying the dispersion history term:
//
//   E^{n+1} = ca * E^n + cb * (stretched curl H - cw * hist)
//   ca = (eps_inf - s) / (eps_inf + s + G)
//   cb = (dt/eps0)     / (eps_inf + s + G)
//   cw = eps0 / dt
//
// s = sigma*dt/(2*eps0) is the conductivity term and G = sum_p 2*Re(B_p) is the implicit coupling
// of the dispersion into E^{n+1} (see dispersion.cu). **When the medium is lossless and
// non-dispersive, s = G = 0 and hist is identically zero**, and the expression above degenerates
// word for word into the explicit update with ca=1 and cb=dt/(eps0*eps). Still one code path.
//
// hist is written by dispersion_pre. In a non-dispersive scene all three pointers point at **the
// same** all-zero array, so the multiply-add here runs as usual and the result is unchanged. That
// costs one extra nx*ny*nz float32 to avoid a branch, and is recorded as known debt.
extern "C" __global__ void update_e(
    float* __restrict__ Ex, float* __restrict__ Ey, float* __restrict__ Ez,
    const float* __restrict__ Hx, const float* __restrict__ Hy,
    const float* __restrict__ Hz,
    float* __restrict__ psi_ex_y, float* __restrict__ psi_ex_z,
    float* __restrict__ psi_ey_z, float* __restrict__ psi_ey_x,
    float* __restrict__ psi_ez_x, float* __restrict__ psi_ez_y,
    const int on_x, const int on_y, const int on_z,     // whether this derivative axis has a PML
    const float* __restrict__ cax, const float* __restrict__ cay,
    const float* __restrict__ caz,                      // (1-k)/(1+k)
    const float* __restrict__ cbx, const float* __restrict__ cby,
    const float* __restrict__ cbz,                      // (dt/(eps0 eps))/(1+k)
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
        Ex[id] = (cax[id] * Ex[id]
                  + cbx[id] * (sy - sz - cw * ph_x)) * pec_y[j] * pec_z[k];
    }
    // Ey: eps dEy/dt = dHx/dz - dHz/dx
    {
        const float dxz = (hx - Hx[IDX(i, j, km)] * mk) * idz_d[k];
        const float dzx = (hz - Hz[IDX(im, j, k)] * mi) * idx_d[i];
        const float sz = stretch(psi_ey_z, id, dxz, az[k], bz[k], kz[k], on_z);
        const float sx = stretch(psi_ey_x, id, dzx, ax[i], bx[i], kx[i], on_x);
        Ey[id] = (cay[id] * Ey[id]
                  + cby[id] * (sz - sx - cw * ph_y)) * pec_z[k] * pec_x[i];
    }
    // Ez: eps dEz/dt = dHy/dx - dHx/dy
    {
        const float dyx = (hy - Hy[IDX(im, j, k)] * mi) * idx_d[i];
        const float dxy = (hx - Hx[IDX(i, jm, k)] * mj) * idy_d[j];
        const float sx = stretch(psi_ez_x, id, dyx, ax[i], bx[i], kx[i], on_x);
        const float sy = stretch(psi_ez_y, id, dxy, ay[j], by[j], ky[j], on_y);
        Ez[id] = (caz[id] * Ez[id]
                  + cbz[id] * (sx - sy - cw * ph_z)) * pec_x[i] * pec_y[j];
    }
}


// Lookup-table version of update_e for (ca, cb): the same expression as update_e, with the
// coefficients coming from a uint16 index into a deduplicated table instead. The gather bit
// pattern is identical.
extern "C" __global__ void update_e_lut(
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
