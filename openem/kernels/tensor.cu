// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Mingdian Liu
// E update for full-tensor eps/sigma cells (FullyAnisotropicMedium), on the "D form plus inverse
// tensor" path.
//
// The continuous equations, semi-implicit and on the same time discretization as the scalar path,
// except that eps and sigma are both 3x3:
//
//   ε₀ ε (E^{n+1}−E^n)/Δt + σ (E^{n+1}+E^n)/2 = curl H − J
//   ⇒ A E^{n+1} = B E^n + (Δt/ε₀)(curl H − J),   A = ε + s,  B = ε − s,  s = σΔt/(2ε₀)
//   ⇒ E^{n+1} = M1 E^n + M2 (curl H − J),        M1 = A⁻¹B,  M2 = (Δt/ε₀)·A⁻¹
//
// sigma need not be symmetric: its antisymmetric part is exactly the gyrotropic (magneto-optic)
// term, A remains invertible and the formulas are unchanged. For diagonal eps and sigma, M1 and M2
// degenerate word for word into the scalar path's ca and cb; see tests/test_tensor.py.
//
// On a Yee grid the three E components sit at different positions (see the comment at the top of
// yee.cu), so an off-position component is taken as a **4-neighbour average**: to move component j
// to the position of component i, take {0,+1} along axis i (edge to center) and {-1,0} along axis j
// (center to edge), each multiplied by the boundary mask (periodic wraparound, or 0 at an
// absorbing end), with a constant weight of 1/4. That is exact on a uniform grid; the weighting on
// a non-uniform grid is a documented concession.
//
// The hook into the main loop reuses the "curl passthrough" track from dispersion_mix: a tensor
// cell has ca=0 and cb=1, so update_e writes the stretched curl plus source term into the E array
// unchanged. This file stores M1*<E^n> **before** update_e using E^n (tensor_dot, first pass),
// reads back c = curl - J **after** every E-side injection to form M2*<c> (second pass), and
// finally writes back with tensor_scatter. The neighbours that <c_j> reads must themselves be
// passthrough cells, which the entry table guarantees by carrying a ring of diagonal halo
// (model.TensorEps.validate locks the closure). A scene with no tensor cells never launches these
// three kernels, so being bitwise unchanged holds by construction, as with the absorber.
//
// The gather (second pass) reads the E array and writes only the private enew; only scatter writes
// back into E. The two passes are separate because the neighbour average reads cells that other
// entries are about to write, and writing in place would race.

#define IDX3(i, j, k) (((i) * ny + (j)) * nz + (k))

__device__ __forceinline__ float pick3(
    const float* __restrict__ Ex, const float* __restrict__ Ey,
    const float* __restrict__ Ez, const int c, const int id)
{
    return (c == 0) ? Ex[id] : ((c == 1) ? Ey[id] : Ez[id]);
}

// out[e] = wbase*base[e] + sum_b m[e][b] * <F_b> at the Yee point of component comp[e]
// First pass:  out=rhs,  base=rhs, wbase=0, m=M1, F=E^n
// Second pass: out=enew, base=rhs, wbase=1, m=M2, F=(stretched curl - J)
extern "C" __global__ void tensor_dot(
    float* __restrict__ out,
    const float* __restrict__ base, const float wbase,
    const float* __restrict__ Fx, const float* __restrict__ Fy,
    const float* __restrict__ Fz,
    const float* __restrict__ m,             // (n_entry, 3) row-major
    const int* __restrict__ comp, const int* __restrict__ cell,
    const int* __restrict__ nxt_x, const int* __restrict__ nxt_y,
    const int* __restrict__ nxt_z,
    const int* __restrict__ prv_x, const int* __restrict__ prv_y,
    const int* __restrict__ prv_z,
    const float* __restrict__ mnx_x, const float* __restrict__ mnx_y,
    const float* __restrict__ mnx_z,
    const float* __restrict__ mpv_x, const float* __restrict__ mpv_y,
    const float* __restrict__ mpv_z,
    const int ny, const int nz, const int n_entry)
{
    const int e = blockIdx.x * blockDim.x + threadIdx.x;
    if (e >= n_entry) return;

    const int a = comp[e];
    const int id = cell[e];
    const int k = id % nz;
    const int j = (id / nz) % ny;
    const int i = id / (nz * ny);

    const int* nxt[3] = {nxt_x, nxt_y, nxt_z};
    const int* prv[3] = {prv_x, prv_y, prv_z};
    const float* mnx[3] = {mnx_x, mnx_y, mnx_z};
    const float* mpv[3] = {mpv_x, mpv_y, mpv_z};
    const float* F[3] = {Fx, Fy, Fz};
    const int c0[3] = {i, j, k};

    float acc = wbase * base[e] + m[e * 3 + a] * pick3(Fx, Fy, Fz, a, id);
    for (int b = 0; b < 3; ++b) {
        if (b == a) continue;
        const float mb = m[e * 3 + b];
        // axis a: edge to center, {0, +1}; axis b: center to edge, {-1, 0}
        int ca[3] = {c0[0], c0[1], c0[2]};
        ca[a] = nxt[a][c0[a]];
        int cb[3] = {c0[0], c0[1], c0[2]};
        cb[b] = prv[b][c0[b]];
        int cab[3] = {ca[0], ca[1], ca[2]};
        cab[b] = cb[b];
        const float wa = mnx[a][c0[a]];
        const float wb = mpv[b][c0[b]];
        const float avg = 0.25f * (
            F[b][IDX3(c0[0], c0[1], c0[2])]
            + wa * F[b][IDX3(ca[0], ca[1], ca[2])]
            + wb * F[b][IDX3(cb[0], cb[1], cb[2])]
            + wa * wb * F[b][IDX3(cab[0], cab[1], cab[2])]);
        acc += mb * avg;
    }
    out[e] = acc;
}

// Write back E^{n+1}. The pec factors follow the same convention as update_e: component a is
// multiplied by the pec of the **other two axes**. Tangential E is identically 0 on the wall at
// i=0 of an absorbing axis, and a passthrough cell is no exception; otherwise the tensor path
// would lift the PML's floor.
extern "C" __global__ void tensor_scatter(
    float* __restrict__ Ex, float* __restrict__ Ey, float* __restrict__ Ez,
    const float* __restrict__ enew,
    const int* __restrict__ comp, const int* __restrict__ cell,
    const float* __restrict__ pec_x, const float* __restrict__ pec_y,
    const float* __restrict__ pec_z,
    const int ny, const int nz, const int n_entry)
{
    const int e = blockIdx.x * blockDim.x + threadIdx.x;
    if (e >= n_entry) return;

    const int a = comp[e];
    const int id = cell[e];
    const int k = id % nz;
    const int j = (id / nz) % ny;
    const int i = id / (nz * ny);

    const float p[3] = {pec_x[i], pec_y[j], pec_z[k]};
    const float v = enew[e] * p[(a + 1) % 3] * p[(a + 2) % 3];
    if (a == 0)      Ex[id] = v;
    else if (a == 1) Ey[id] = v;
    else             Ez[id] = v;
}
