// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Mingdian Liu
// Adiabatic absorber (Tidy3D's `Absorber`): a **matched** lossy medium graded layer by layer.
//
// Under the matching condition sigma_m/mu = sigma_e/eps = Gamma, Maxwell's equations are
// equivalent to lossless evolution multiplied by exp(-Gamma*t); substituting E,H = e^{-Gamma t}
// (E',H') recovers the lossless equations. So the implementation is: once per step, multiply the E
// family and the H family inside the absorbing layer by **the same** exp(-Gamma*dt) profile.
// Matching holds by construction, no permeability array is needed, and the main update, the ADE
// and the mix kernels are untouched. Those were exactly the coupling points where the old psi
// recursion and an ADE mismatch used to diverge.
//
// Staggering: along the absorbing axis, the **normal E component sits at the half cell** (Ex is
// offset by half a cell along x) and the tangential ones at the whole cell. H is the opposite (Hx
// is at the whole cell along x, Hy and Hz at the half cell). Both profiles are sampled under the
// same convention as cpml.steps.
//
// decay_* is stored in **global index order** within the slab (l = g - g0), with g0=0 at the lo end
// and g0 = n_ax - nlay at the hi end, matching the output order of cpml.profile.
#define IDX(i, j, k) (((i) * ny + (j)) * nz + (k))

extern "C" __global__ void absorb_slab(
    float* __restrict__ Fx, float* __restrict__ Fy, float* __restrict__ Fz,
    const float* __restrict__ decay_int,    // (nlay,) whole-cell positions
    const float* __restrict__ decay_half,   // (nlay,) half-cell positions
    const int is_e,                         // 1 = E family, 0 = H family
    const int axis,                         // absorbing axis, 0/1/2
    const int i0, const int j0, const int k0,   // slab bounding-box origin
    const int bi, const int bj, const int bk,   // slab bounding-box size
    const int nx, const int ny, const int nz)
{
    // Thread layout: flatten the bounding box (j, k) into one dimension with **k innermost**, and
    // give i to blockIdx.y.
    //
    // The old spelling mapped threadIdx.x to t2 = (axis==2) ? 1 : 2. With z as the absorbing axis
    // that makes t2 = j, so adjacent threads are nz floats apart and one warp touches 32 cache
    // lines, a 32x traffic amplification. On one scene the z faces held 63% of the slab cells and
    // it measured 12.68 ms per step against a roofline of 0.59 ms, a factor of 21. An independent
    // reference implementation had recorded the same class of bug on a different axis, fixed the
    // same way by remapping threads per axis.
    //
    // Why not simply run threadIdx.x along k: the slab is only nlay wide in k (40 in the larger
    // scenes, 12 in a small test), far less than a warp, leaving lane utilization at 19 to 62%,
    // which measured as a net loss on small scenes. Flattened, lane utilization is about 100% and a
    // warp's addresses are contiguous within a k row (only crossing rows changes cache line), which
    // at nlay=40 averages about 1.25 cache lines per warp.
    //
    // The layer number is the offset along the absorbing axis within the bounding box; i0/j0/k0 on
    // that axis is exactly g0.
    const int tid = blockIdx.x * blockDim.x + threadIdx.x;
    const int ii = blockIdx.y;
    if (ii >= bi || tid >= bj * bk) return;
    const int jj = tid / bk;
    const int kk = tid - jj * bk;

    const int l = (axis == 0) ? ii : ((axis == 1) ? jj : kk);
    const int id = IDX(i0 + ii, j0 + jj, k0 + kk);

    const float di = decay_int[l];
    const float dh = decay_half[l];
    float* F[3] = {Fx, Fy, Fz};
    #pragma unroll
    for (int comp = 0; comp < 3; ++comp) {
        // E: the normal component is at the half cell; H: the tangential ones are
        const int half = is_e ? (comp == axis) : (comp != axis);
        F[comp][id] *= half ? dh : di;
    }
}

// Decay of the ADE polarization history P inside the absorbing layer is **not** done here. It is
// folded into dispersion_post in kernels/dispersion.cu, which saves one read-modify-write pass
// over p; see the note there.
