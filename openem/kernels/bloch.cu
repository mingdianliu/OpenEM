// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Mingdian Liu
// **Complex-field** kernel variants for Bloch scenes, where k is non-zero.
//
// Relationship to kernels/yee.cu: the same difference, CPML and mask structure, with each field
// replaced by a pair of real arrays, re and im. The CPML coefficients (a, b, invk), the spacings,
// ca/cb and pec are all **real**, so applying them to a complex field just means computing re and
// im separately; by linearity the psi recursion holds independently for each part.
//
// The complex phase of k appears **only** in neighbour reads that wrap around periodically: the
// mnx/mpv masks become complex (grid.Axis.index_tables(bloch_k)), with exp(+-i*2*pi*k) on the
// wrapping elements, 1+0i elsewhere, and 0+0i outside an absorbing end. As on the real path, the
// mask expresses both topology and phase, so the kernel stays branch-free.
//
// These kernels **launch only when Scene.any_bloch holds** (or a test forces it). A scene without
// Bloch never compiles or launches them and is bitwise unchanged, on a single code path, the same
// way the absorber and the modulation kernels do it.
//
// The complex subset excludes dispersion, tensors, modulation and loss (the host fails closed
// otherwise), so update_e_c has no hist term and ca/cb are the explicit coefficients of a pure eps.
//
// Memory layout and launch match yee.cu: C order (nx, ny, nz), with the thread x dimension mapped
// to k.

#define IDX(i, j, k) (((i) * ny + (j)) * nz + (k))

// psi recursion plus coordinate stretching, the same as stretch in yee.cu; called once for re and
// once for im, since the coefficients are real.
__device__ __forceinline__ float stretch_c(
    float* __restrict__ psi, const int id, const float d,
    const float a, const float b, const float invk)
{
    const float p = b * psi[id] + a * d;
    psi[id] = p;
    return invk * d + p;
}

extern "C" __global__ void update_h_c(
    float* __restrict__ Hxr, float* __restrict__ Hxi,
    float* __restrict__ Hyr, float* __restrict__ Hyi,
    float* __restrict__ Hzr, float* __restrict__ Hzi,
    const float* __restrict__ Exr, const float* __restrict__ Exi,
    const float* __restrict__ Eyr, const float* __restrict__ Eyi,
    const float* __restrict__ Ezr, const float* __restrict__ Ezi,
    float* __restrict__ psi_hx_y_r, float* __restrict__ psi_hx_y_i,
    float* __restrict__ psi_hx_z_r, float* __restrict__ psi_hx_z_i,
    float* __restrict__ psi_hy_z_r, float* __restrict__ psi_hy_z_i,
    float* __restrict__ psi_hy_x_r, float* __restrict__ psi_hy_x_i,
    float* __restrict__ psi_hz_x_r, float* __restrict__ psi_hz_x_i,
    float* __restrict__ psi_hz_y_r, float* __restrict__ psi_hz_y_i,
    const float* __restrict__ idx_p, const float* __restrict__ idy_p,
    const float* __restrict__ idz_p,
    const float* __restrict__ ax, const float* __restrict__ bx,
    const float* __restrict__ kx,
    const float* __restrict__ ay, const float* __restrict__ by,
    const float* __restrict__ ky,
    const float* __restrict__ az, const float* __restrict__ bz,
    const float* __restrict__ kz,
    const int* __restrict__ nxt_x, const int* __restrict__ nxt_y,
    const int* __restrict__ nxt_z,
    const float* __restrict__ mnxr_x, const float* __restrict__ mnxi_x,
    const float* __restrict__ mnxr_y, const float* __restrict__ mnxi_y,
    const float* __restrict__ mnxr_z, const float* __restrict__ mnxi_z,
    const float ch,
    const int nx, const int ny, const int nz)
{
    const int k = blockIdx.x * blockDim.x + threadIdx.x;
    const int j = blockIdx.y * blockDim.y + threadIdx.y;
    const int i = blockIdx.z;
    if (i >= nx || j >= ny || k >= nz) return;

    const int id = IDX(i, j, k);
    const int ip = nxt_x[i], jp = nxt_y[j], kp = nxt_z[k];
    const float mir = mnxr_x[i], mii = mnxi_x[i];
    const float mjr = mnxr_y[j], mji = mnxi_y[j];
    const float mkr = mnxr_z[k], mki = mnxi_z[k];

    // Hx: mu0 dHx/dt = -(dEz/dy - dEy/dz)
    {
        const int ja = IDX(i, jp, k), ka = IDX(i, j, kp);
        const float zr = Ezr[ja], zi = Ezi[ja], yr = Eyr[ka], yi = Eyi[ka];
        const float dzy_r = (zr * mjr - zi * mji - Ezr[id]) * idy_p[j];
        const float dzy_i = (zr * mji + zi * mjr - Ezi[id]) * idy_p[j];
        const float dyz_r = (yr * mkr - yi * mki - Eyr[id]) * idz_p[k];
        const float dyz_i = (yr * mki + yi * mkr - Eyi[id]) * idz_p[k];
        const float syr = stretch_c(psi_hx_y_r, id, dzy_r, ay[j], by[j], ky[j]);
        const float syi = stretch_c(psi_hx_y_i, id, dzy_i, ay[j], by[j], ky[j]);
        const float szr = stretch_c(psi_hx_z_r, id, dyz_r, az[k], bz[k], kz[k]);
        const float szi = stretch_c(psi_hx_z_i, id, dyz_i, az[k], bz[k], kz[k]);
        Hxr[id] -= ch * (syr - szr);
        Hxi[id] -= ch * (syi - szi);
    }
    // Hy: mu0 dHy/dt = -(dEx/dz - dEz/dx)
    {
        const int ka = IDX(i, j, kp), ia = IDX(ip, j, k);
        const float xr = Exr[ka], xi = Exi[ka], zr = Ezr[ia], zi = Ezi[ia];
        const float dxz_r = (xr * mkr - xi * mki - Exr[id]) * idz_p[k];
        const float dxz_i = (xr * mki + xi * mkr - Exi[id]) * idz_p[k];
        const float dzx_r = (zr * mir - zi * mii - Ezr[id]) * idx_p[i];
        const float dzx_i = (zr * mii + zi * mir - Ezi[id]) * idx_p[i];
        const float szr = stretch_c(psi_hy_z_r, id, dxz_r, az[k], bz[k], kz[k]);
        const float szi = stretch_c(psi_hy_z_i, id, dxz_i, az[k], bz[k], kz[k]);
        const float sxr = stretch_c(psi_hy_x_r, id, dzx_r, ax[i], bx[i], kx[i]);
        const float sxi = stretch_c(psi_hy_x_i, id, dzx_i, ax[i], bx[i], kx[i]);
        Hyr[id] -= ch * (szr - sxr);
        Hyi[id] -= ch * (szi - sxi);
    }
    // Hz: mu0 dHz/dt = -(dEy/dx - dEx/dy)
    {
        const int ia = IDX(ip, j, k), ja = IDX(i, jp, k);
        const float yr = Eyr[ia], yi = Eyi[ia], xr = Exr[ja], xi = Exi[ja];
        const float dyx_r = (yr * mir - yi * mii - Eyr[id]) * idx_p[i];
        const float dyx_i = (yr * mii + yi * mir - Eyi[id]) * idx_p[i];
        const float dxy_r = (xr * mjr - xi * mji - Exr[id]) * idy_p[j];
        const float dxy_i = (xr * mji + xi * mjr - Exi[id]) * idy_p[j];
        const float sxr = stretch_c(psi_hz_x_r, id, dyx_r, ax[i], bx[i], kx[i]);
        const float sxi = stretch_c(psi_hz_x_i, id, dyx_i, ax[i], bx[i], kx[i]);
        const float syr = stretch_c(psi_hz_y_r, id, dxy_r, ay[j], by[j], ky[j]);
        const float syi = stretch_c(psi_hz_y_i, id, dxy_i, ay[j], by[j], ky[j]);
        Hzr[id] -= ch * (sxr - syr);
        Hzi[id] -= ch * (sxi - syi);
    }
}

// Complex variant of the E update. The subset is lossless and non-dispersive, so ca/cb are the
// explicit coefficients of a pure eps, real and per cell, and there is no hist term. pec is a real
// mask and behaves as in yee.cu.
extern "C" __global__ void update_e_c(
    float* __restrict__ Exr, float* __restrict__ Exi,
    float* __restrict__ Eyr, float* __restrict__ Eyi,
    float* __restrict__ Ezr, float* __restrict__ Ezi,
    const float* __restrict__ Hxr, const float* __restrict__ Hxi,
    const float* __restrict__ Hyr, const float* __restrict__ Hyi,
    const float* __restrict__ Hzr, const float* __restrict__ Hzi,
    float* __restrict__ psi_ex_y_r, float* __restrict__ psi_ex_y_i,
    float* __restrict__ psi_ex_z_r, float* __restrict__ psi_ex_z_i,
    float* __restrict__ psi_ey_z_r, float* __restrict__ psi_ey_z_i,
    float* __restrict__ psi_ey_x_r, float* __restrict__ psi_ey_x_i,
    float* __restrict__ psi_ez_x_r, float* __restrict__ psi_ez_x_i,
    float* __restrict__ psi_ez_y_r, float* __restrict__ psi_ez_y_i,
    const float* __restrict__ cax, const float* __restrict__ cay,
    const float* __restrict__ caz,
    const float* __restrict__ cbx, const float* __restrict__ cby,
    const float* __restrict__ cbz,
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
    const float* __restrict__ mpvr_x, const float* __restrict__ mpvi_x,
    const float* __restrict__ mpvr_y, const float* __restrict__ mpvi_y,
    const float* __restrict__ mpvr_z, const float* __restrict__ mpvi_z,
    const float* __restrict__ pec_x, const float* __restrict__ pec_y,
    const float* __restrict__ pec_z,
    const int nx, const int ny, const int nz)
{
    const int k = blockIdx.x * blockDim.x + threadIdx.x;
    const int j = blockIdx.y * blockDim.y + threadIdx.y;
    const int i = blockIdx.z;
    if (i >= nx || j >= ny || k >= nz) return;

    const int id = IDX(i, j, k);
    const int im = prv_x[i], jm = prv_y[j], km = prv_z[k];
    const float mir = mpvr_x[i], mii = mpvi_x[i];
    const float mjr = mpvr_y[j], mji = mpvi_y[j];
    const float mkr = mpvr_z[k], mki = mpvi_z[k];

    // Ex: eps dEx/dt = dHz/dy - dHy/dz
    {
        const int ja = IDX(i, jm, k), ka = IDX(i, j, km);
        const float zr = Hzr[ja], zi = Hzi[ja], yr = Hyr[ka], yi = Hyi[ka];
        const float dzy_r = (Hzr[id] - (zr * mjr - zi * mji)) * idy_d[j];
        const float dzy_i = (Hzi[id] - (zr * mji + zi * mjr)) * idy_d[j];
        const float dyz_r = (Hyr[id] - (yr * mkr - yi * mki)) * idz_d[k];
        const float dyz_i = (Hyi[id] - (yr * mki + yi * mkr)) * idz_d[k];
        const float syr = stretch_c(psi_ex_y_r, id, dzy_r, ay[j], by[j], ky[j]);
        const float syi = stretch_c(psi_ex_y_i, id, dzy_i, ay[j], by[j], ky[j]);
        const float szr = stretch_c(psi_ex_z_r, id, dyz_r, az[k], bz[k], kz[k]);
        const float szi = stretch_c(psi_ex_z_i, id, dyz_i, az[k], bz[k], kz[k]);
        const float p = pec_y[j] * pec_z[k];
        Exr[id] = (cax[id] * Exr[id] + cbx[id] * (syr - szr)) * p;
        Exi[id] = (cax[id] * Exi[id] + cbx[id] * (syi - szi)) * p;
    }
    // Ey: eps dEy/dt = dHx/dz - dHz/dx
    {
        const int ka = IDX(i, j, km), ia = IDX(im, j, k);
        const float xr = Hxr[ka], xi = Hxi[ka], zr = Hzr[ia], zi = Hzi[ia];
        const float dxz_r = (Hxr[id] - (xr * mkr - xi * mki)) * idz_d[k];
        const float dxz_i = (Hxi[id] - (xr * mki + xi * mkr)) * idz_d[k];
        const float dzx_r = (Hzr[id] - (zr * mir - zi * mii)) * idx_d[i];
        const float dzx_i = (Hzi[id] - (zr * mii + zi * mir)) * idx_d[i];
        const float szr = stretch_c(psi_ey_z_r, id, dxz_r, az[k], bz[k], kz[k]);
        const float szi = stretch_c(psi_ey_z_i, id, dxz_i, az[k], bz[k], kz[k]);
        const float sxr = stretch_c(psi_ey_x_r, id, dzx_r, ax[i], bx[i], kx[i]);
        const float sxi = stretch_c(psi_ey_x_i, id, dzx_i, ax[i], bx[i], kx[i]);
        const float p = pec_z[k] * pec_x[i];
        Eyr[id] = (cay[id] * Eyr[id] + cby[id] * (szr - sxr)) * p;
        Eyi[id] = (cay[id] * Eyi[id] + cby[id] * (szi - sxi)) * p;
    }
    // Ez: eps dEz/dt = dHy/dx - dHx/dy
    {
        const int ia = IDX(im, j, k), ja = IDX(i, jm, k);
        const float yr = Hyr[ia], yi = Hyi[ia], xr = Hxr[ja], xi = Hxi[ja];
        const float dyx_r = (Hyr[id] - (yr * mir - yi * mii)) * idx_d[i];
        const float dyx_i = (Hyi[id] - (yr * mii + yi * mir)) * idx_d[i];
        const float dxy_r = (Hxr[id] - (xr * mjr - xi * mji)) * idy_d[j];
        const float dxy_i = (Hxi[id] - (xr * mji + xi * mjr)) * idy_d[j];
        const float sxr = stretch_c(psi_ez_x_r, id, dyx_r, ax[i], bx[i], kx[i]);
        const float sxi = stretch_c(psi_ez_x_i, id, dyx_i, ax[i], bx[i], kx[i]);
        const float syr = stretch_c(psi_ez_y_r, id, dxy_r, ay[j], by[j], ky[j]);
        const float syi = stretch_c(psi_ez_y_i, id, dxy_i, ay[j], by[j], ky[j]);
        const float p = pec_x[i] * pec_y[j];
        Ezr[id] = (caz[id] * Ezr[id] + cbz[id] * (sxr - syr)) * p;
        Ezi[id] = (caz[id] * Ezi[id] + cbz[id] * (sxi - syi)) * p;
    }
}

// Complex variant of the point current source: it injects the full complex amplitude
// exp(i*phi)*g(t) (waveform.amp_half_complex). Coefficients and timing match inject_dipole in
// source_dft.cu; only re and im each get their own shot.
extern "C" __global__ void inject_dipole_c(
    float* __restrict__ Exr, float* __restrict__ Exi,
    float* __restrict__ Eyr, float* __restrict__ Eyi,
    float* __restrict__ Ezr, float* __restrict__ Ezi,
    const float* __restrict__ cex, const float* __restrict__ cey,
    const float* __restrict__ cez,
    const int* __restrict__ flat, const int* __restrict__ comp,
    const float* __restrict__ coef, const int* __restrict__ src_of,
    const double* __restrict__ amp_re,   // (n_src, nsteps1) complex amplitude at half steps
    const double* __restrict__ amp_im,
    const int step, const int nsteps1, const int n)
{
    const int t = blockIdx.x * blockDim.x + threadIdx.x;
    if (t >= n || step >= nsteps1) return;   // backstop guard, for the same reason as inject_dipole

    float* Er[3] = {Exr, Eyr, Ezr};
    float* Ei[3] = {Exi, Eyi, Ezi};
    const float* cp[3] = {cex, cey, cez};
    const int c = comp[t];
    const int id = flat[t];
    const size_t o = (size_t)src_of[t] * nsteps1 + step;
    const float s = -cp[c][id] * coef[t];
    atomicAdd(&Er[c][id], s * (float)amp_re[o]);
    atomicAdd(&Ei[c][id], s * (float)amp_im[o]);
}

// Complex variant of the point magnetic current source, using the whole-step complex amplitude
// waveform.amp_int_complex; otherwise as inject_dipole_h.
extern "C" __global__ void inject_dipole_h_c(
    float* __restrict__ Hxr, float* __restrict__ Hxi,
    float* __restrict__ Hyr, float* __restrict__ Hyi,
    float* __restrict__ Hzr, float* __restrict__ Hzi,
    const float ch,
    const int* __restrict__ flat, const int* __restrict__ comp,
    const float* __restrict__ coef, const int* __restrict__ src_of,
    const double* __restrict__ amp_re,   // (n_src, nsteps1) complex amplitude at whole steps
    const double* __restrict__ amp_im,
    const int step, const int nsteps1, const int n)
{
    const int t = blockIdx.x * blockDim.x + threadIdx.x;
    if (t >= n || step >= nsteps1) return;

    float* Hr[3] = {Hxr, Hyr, Hzr};
    float* Hi[3] = {Hxi, Hyi, Hzi};
    const int c = comp[t];
    const int id = flat[t];
    const size_t o = (size_t)src_of[t] * nsteps1 + step;
    const float s = -ch * coef[t];
    atomicAdd(&Hr[c][id], s * (float)amp_re[o]);
    atomicAdd(&Hi[c][id], s * (float)amp_im[o]);
}

// Complex variant of time-domain sampling: layout, the two-shot H time average, and the weights
// and add semantics all match sample_time_box.
extern "C" __global__ void sample_time_box_c(
    float* __restrict__ out_re, float* __restrict__ out_im,
    const float* __restrict__ Exr, const float* __restrict__ Exi,
    const float* __restrict__ Eyr, const float* __restrict__ Eyi,
    const float* __restrict__ Ezr, const float* __restrict__ Ezi,
    const float* __restrict__ Hxr, const float* __restrict__ Hxi,
    const float* __restrict__ Hyr, const float* __restrict__ Hyi,
    const float* __restrict__ Hzr, const float* __restrict__ Hzi,
    const int* __restrict__ comps, const float* __restrict__ weights, const int nc,
    const int slot, const int nslots, const int add,
    const int i0, const int j0, const int k0,
    const int ni, const int nj, const int nk,
    const int nx, const int ny, const int nz)
{
    const int cc = blockIdx.x * blockDim.x + threadIdx.x;
    const int b = blockIdx.y * blockDim.y + threadIdx.y;
    const int a = blockIdx.z;
    if (a >= ni || b >= nj || cc >= nk) return;

    const int id = IDX(i0 + a, j0 + b, k0 + cc);
    const float* Fr[6] = {Exr, Eyr, Ezr, Hxr, Hyr, Hzr};
    const float* Fi[6] = {Exi, Eyi, Ezi, Hxi, Hyi, Hzi};
    const int cells = ni * nj * nk;
    const int off = (a * nj + b) * nk + cc;
    for (int c = 0; c < nc; ++c) {
        const float w = weights[c];
        if (w == 0.0f) continue;
        const size_t o = ((size_t)c * nslots + slot) * cells + off;
        const float vr = w * Fr[comps[c]][id];
        const float vi = w * Fi[comps[c]][id];
        out_re[o] = add ? out_re[o] + vr : vr;
        out_im[o] = add ? out_im[o] + vi : vi;
    }
}


// ---- Complex single-face TF/SF injection (oblique plane wave, oblique GaussianBeam, or a mode
// source in a Bloch scene) ----
// Same layout as the real inject_mode_h/e, except that the incident profile is **genuinely
// complex** and so is the target field: the real path injects Re(profile*amp), while here the full
// complex product profile*amp is added into both re and im.
// profile = e_re + i*e_im, at the respective Yee tangential positions, with the transverse phase
// e^{i k_t . r_t} already folded in. amp = amp_re + i*amp_im, an analytic signal stored as an
// interleaved (n,2) per-step table. The tangential pair follows the cyclic order t1/t2, and the
// three axis branches select the injection face, as in the real version in source_dft.cu. The
// coefficient coef carries the sign and dt/mu0/dl.
extern "C" __global__ void inject_mode_h_c(
    float* __restrict__ Ht1r, float* __restrict__ Ht1i,
    float* __restrict__ Ht2r, float* __restrict__ Ht2i,
    const float* __restrict__ e1_re, const float* __restrict__ e1_im,
    const float* __restrict__ e2_re, const float* __restrict__ e2_im,
    const double* __restrict__ amp_tab, const int step,
    const float coef_1, const float coef_2,
    const int kh, const int axis, const int n1, const int n2, const int w2,
    const int nx, const int ny, const int nz)
{
    const int o2 = blockIdx.x * blockDim.x + threadIdx.x;   // t2
    const int o1 = blockIdx.y * blockDim.y + threadIdx.y;   // t1
    if (o1 >= n1 || o2 >= n2) return;
    const float ar = (float)amp_tab[2 * step];
    const float ai = (float)amp_tab[2 * step + 1];
    const int o  = o1 * w2 + o2;
    const int id = axis == 0 ? IDX(kh, o1, o2)
                 : axis == 1 ? IDX(o2, kh, o1) : IDX(o1, o2, kh);
    // Ht1 ← coef_1 · (e2 · amp)
    Ht1r[id] += coef_1 * (e2_re[o] * ar - e2_im[o] * ai);
    Ht1i[id] += coef_1 * (e2_re[o] * ai + e2_im[o] * ar);
    // Ht2 ← coef_2 · (e1 · amp)
    Ht2r[id] += coef_2 * (e1_re[o] * ar - e1_im[o] * ai);
    Ht2i[id] += coef_2 * (e1_re[o] * ai + e1_im[o] * ar);
}

extern "C" __global__ void inject_mode_e_c(
    float* __restrict__ Et1r, float* __restrict__ Et1i,
    float* __restrict__ Et2r, float* __restrict__ Et2i,
    const float* __restrict__ ce1, const float* __restrict__ ce2,
    const float* __restrict__ h1_re, const float* __restrict__ h1_im,
    const float* __restrict__ h2_re, const float* __restrict__ h2_im,
    const double* __restrict__ amp_tab, const int step,
    const float coef_1, const float coef_2,
    const int ks, const int axis, const int n1, const int n2, const int w2,
    const int nx, const int ny, const int nz)
{
    const int o2 = blockIdx.x * blockDim.x + threadIdx.x;
    const int o1 = blockIdx.y * blockDim.y + threadIdx.y;
    if (o1 >= n1 || o2 >= n2) return;
    const float ar = (float)amp_tab[2 * step];
    const float ai = (float)amp_tab[2 * step + 1];
    const int o  = o1 * w2 + o2;
    const int id = axis == 0 ? IDX(ks, o1, o2)
                 : axis == 1 ? IDX(o2, ks, o1) : IDX(o1, o2, ks);
    const float c1 = ce1[id] * coef_1, c2 = ce2[id] * coef_2;
    Et1r[id] += c1 * (h2_re[o] * ar - h2_im[o] * ai);
    Et1i[id] += c1 * (h2_re[o] * ai + h2_im[o] * ar);
    Et2r[id] += c2 * (h1_re[o] * ar - h1_im[o] * ai);
    Et2i[id] += c2 * (h1_re[o] * ai + h1_im[o] * ar);
}


// ---- Frequency-domain DFT accumulation for complex fields (FieldMonitor and FluxMonitor on the
// Bloch path) ----
// The real version accumulates sum val*e^{i*omega*t}*apod = re + i*im. For a complex field
// u = u_re + i*u_im:
//   phasor = Σ (u_re+i·u_im)(cos+i·sin)·apod
//     re += (u_re·cos − u_im·sin)·apod
//     im += (u_re·sin + u_im·cos)·apod
// Layout matches accumulate_dft_box: (j,k) flattened onto grid.x, frequency onto grid.y, i onto
// grid.z. step is a scalar, since the Bloch path does no graph capture. The phase table tab is
// (nf,4), i.e. [cE, sE, cH, sH] per frequency.
extern "C" __global__ void accumulate_dft_box_c(
    double* __restrict__ re, double* __restrict__ im,
    const float* __restrict__ Exr, const float* __restrict__ Exi,
    const float* __restrict__ Eyr, const float* __restrict__ Eyi,
    const float* __restrict__ Ezr, const float* __restrict__ Ezi,
    const float* __restrict__ Hxr, const float* __restrict__ Hxi,
    const float* __restrict__ Hyr, const float* __restrict__ Hyi,
    const float* __restrict__ Hzr, const float* __restrict__ Hzi,
    const double* __restrict__ tab,     // (nf,4) this step's phase times apod
    const int step,
    const int i0, const int j0, const int k0,
    const int ni, const int nj, const int nk, const int nf,
    const int nx, const int ny, const int nz)
{
    const int tid = blockIdx.x * blockDim.x + threadIdx.x;
    const int f = blockIdx.y;
    const int a = blockIdx.z;
    if (a >= ni || tid >= nj * nk || f >= nf) return;
    const int b = tid / nk;
    const int cc = tid - b * nk;
    const int id = IDX(i0 + a, j0 + b, k0 + cc);
    const double vr[6] = {(double)Exr[id], (double)Eyr[id], (double)Ezr[id],
                          (double)Hxr[id], (double)Hyr[id], (double)Hzr[id]};
    const double vi[6] = {(double)Exi[id], (double)Eyi[id], (double)Ezi[id],
                          (double)Hxi[id], (double)Hyi[id], (double)Hzi[id]};
    const int cells = ni * nj * nk;
    const int off = (a * nj + b) * nk + cc;
    const double* ph = tab + (size_t)f * 4;
    const double cE = ph[0], sE = ph[1], cH = ph[2], sH = ph[3];
    for (int c = 0; c < 6; ++c) {
        const double cr = (c < 3) ? cE : cH;
        const double si = (c < 3) ? sE : sH;
        const size_t o = (size_t)(c * nf + f) * cells + off;
        re[o] += vr[c] * cr - vi[c] * si;
        im[o] += vr[c] * si + vi[c] * cr;
    }
}


// ---- Colocated flux DFT for complex fields (FluxMonitor and DiffractionMonitor on the Bloch
// path) ----
// As in the real accumulate_dft: E_t1 and E_t2 are taken on the plane, while H_t1 and H_t2 are
// averaged over half a cell along the normal to colocate them onto the E plane.
// For a complex field: phasor = sum u*e^{i*omega*t}*apod, so
//   re += u_re*cos - u_im*sin,  im += u_re*sin + u_im*cos.
// Output is a (4,nf,n1,n2) complex phasor, one array each for re and im, in the order
// E_t1, E_t2, H_t1, H_t2. step is a scalar; tab is (nf,4) holding [cE, sE, cH, sH] for this step.
extern "C" __global__ void accumulate_dft_c(
    double* __restrict__ re, double* __restrict__ im,
    const float* __restrict__ Exr, const float* __restrict__ Exi,
    const float* __restrict__ Eyr, const float* __restrict__ Eyi,
    const float* __restrict__ Ezr, const float* __restrict__ Ezi,
    const float* __restrict__ Hxr, const float* __restrict__ Hxi,
    const float* __restrict__ Hyr, const float* __restrict__ Hyi,
    const float* __restrict__ Hzr, const float* __restrict__ Hzi,
    const double* __restrict__ tab, const int step,
    const int km, const int axis, const int nf,
    const int nz_log, const int nx, const int ny, const int nz)
{
    const int t1 = (axis == 0) ? 1 : 0;
    const int t2 = (axis == 2) ? 1 : 2;
    const int n1 = (t1 == 0) ? nx : ny;
    const int n2 = (t2 == 1) ? ny : nz_log;
    const int b = blockIdx.x * blockDim.x + threadIdx.x;
    const int a = blockIdx.y * blockDim.y + threadIdx.y;
    const int f = blockIdx.z;
    if (a >= n1 || b >= n2 || f >= nf) return;
    int c[3]; c[axis] = km; c[t1] = a; c[t2] = b;
    const int id = IDX(c[0], c[1], c[2]);
    c[axis] = km - 1;
    const int idm = IDX(c[0], c[1], c[2]);
    const float* Er[3] = {Exr, Eyr, Ezr};
    const float* Ei[3] = {Exi, Eyi, Ezi};
    const float* Hr[3] = {Hxr, Hyr, Hzr};
    const float* Hi[3] = {Hxi, Hyi, Hzi};
    const double vr[4] = {(double)Er[t1][id], (double)Er[t2][id],
        0.5*((double)Hr[t1][id]+(double)Hr[t1][idm]),
        0.5*((double)Hr[t2][id]+(double)Hr[t2][idm])};
    const double vi[4] = {(double)Ei[t1][id], (double)Ei[t2][id],
        0.5*((double)Hi[t1][id]+(double)Hi[t1][idm]),
        0.5*((double)Hi[t2][id]+(double)Hi[t2][idm])};
    const double* ph = tab + (size_t)f * 4;
    const double cE = ph[0], sE = ph[1], cH = ph[2], sH = ph[3];
    const size_t plane = (size_t)n1 * n2;
    const size_t off = (size_t)a * n2 + b;
    for (int cc = 0; cc < 4; ++cc) {
        const double cr = (cc < 2) ? cE : cH;
        const double si = (cc < 2) ? sE : sH;
        const size_t o = (size_t)(cc * nf + f) * plane + off;
        re[o] += vr[cc] * cr - vi[cc] * si;
        im[o] += vr[cc] * si + vi[cc] * cr;
    }
}
