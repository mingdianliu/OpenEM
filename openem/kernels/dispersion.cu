// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Mingdian Liu
// Auxiliary differential equations (ADE) for dispersive materials: one first-order equation per
// complex pole, discretized with the trapezoidal rule.
//
// Tidy3D reduces every dispersive material to a PoleResidue (``_pole_residue_dict`` is an abstract
// method on ``DispersiveMedium`` that every subclass must implement), so only one model has to be
// understood here:
//
//   ε(ω) = ε_∞ − Σ_i [ c_i/(jω + a_i) + c_i*/(jω + a_i*) ]
//
// Mapping to the "pole q, residue r" form gives q = a_i and r = c_i, with no sign flip; measured
// point by point, the two sides agree on eps to 8.6e-15 relative. The time-domain equation for
// each conjugate pair is **one complex first-order equation**,
//
//   dP/dt = q P + r E,        physical polarization = 2 Re(P)
//
// rather than the real second-order equation from the textbooks. The reason is accuracy: in the
// real second-order form the dE/dt coupling term (non-zero for gold, silver and silicon alike)
// uses a forward difference centered at n+1/2 instead of n, which drops the whole scheme to
// **first order**; gold at dt=1e-17 then has 3.3% error in eps. The complex-pole form with the
// trapezoidal rule is **second order**, giving 0.084% under the same conditions. That was measured
// as a convergence-order table across four materials.
//
// Trapezoidal rule, (P^{n+1}-P^n)/dt = q(P^{n+1}+P^n)/2 + r(E^{n+1}+E^n)/2:
//
//   P^{n+1} = A P^n + B (E^{n+1} + E^n)
//   A = (1 + qΔt/2)/(1 − qΔt/2),   B = (rΔt/2)/(1 − qΔt/2)
//
// The coefficient arrays store **A - 1**, not A; see the note inside dispersion_pre below.
//
// **Storage is packed by dispersive cell-component**, not a dense array. Dispersive materials
// usually occupy a minority of cells (45.8% in the densest case in the validation set, far less in
// the plasmonic ones), so a dense allocation would spend memory on vacuum. The layout is CSR-like:
//
//   comp[e] / cell[e]     which component (0/1/2) and which flat cell index entry e refers to
//   ofs[e] .. ofs[e+1]    the range of entry e's poles in the coefficient arrays
//
// A given (comp, cell) **may appear in only one entry**; a mixed cell merges the poles of both
// materials into a single entry. So writing hist needs no atomics, and the build side guarantees it.
//
// The E update splits into two passes with update_e in between:
//
//   dispersion_pre    advance P to P_half = A P^n + B E^n using E^n, and compute the history term
//                     W = sum 2*Re(P_half - P^n) that update_e needs
//   update_e          E^{n+1} = ca E^n + cb (curl H - (eps0/dt) W)
//   dispersion_post   add B E^{n+1} to obtain the complete P^{n+1}
//
// This needs only **one** P array, with no P^{n-1}, half the storage of the real second-order form.

// Fetch component c at cell id. All three pointers are passed in so the kernel does no pointer
// arithmetic.
__device__ __forceinline__ float pick(
    const float* __restrict__ Ex, const float* __restrict__ Ey,
    const float* __restrict__ Ez, const int c, const int id)
{
    return (c == 0) ? Ex[id] : ((c == 1) ? Ey[id] : Ez[id]);
}

// First pass: P <- A P + B E^n, also writing out the history term update_e needs.
// The coefficients are deduplicated **losslessly** into a table of (am1, b) pairs (lut_*), with one
// int32 index stored per pole. Measured across the whole set: am1 has at most 7 distinct values and
// b at most 526,812, so the table is at most 8.4 MB and stays in L2, while the coefficient traffic
// per pole drops from 16 B to 4 B. The table holds the original values themselves (np.unique does
// not change a bit), so a lookup returns the same float and the result is **bitwise equivalent**.
// This is a different thing entirely from a "quantized LUT", which would change the physics and is
// not allowed.
extern "C" __global__ void dispersion_pre(
    float* __restrict__ p_re, float* __restrict__ p_im,
    float* __restrict__ hist_x, float* __restrict__ hist_y,
    float* __restrict__ hist_z,
    const float* __restrict__ lut_am1_re, const float* __restrict__ lut_am1_im,
    const float* __restrict__ lut_b_re, const float* __restrict__ lut_b_im,
    const int* __restrict__ cidx,
    const int* __restrict__ comp, const int* __restrict__ cell,
    const int* __restrict__ ofs,
    const float* __restrict__ Ex, const float* __restrict__ Ey,
    const float* __restrict__ Ez,
    const int n_entry)
{
    const int e = blockIdx.x * blockDim.x + threadIdx.x;
    if (e >= n_entry) return;

    const int c = comp[e], id = cell[e];
    const float en = pick(Ex, Ey, Ez, c, id);
    const int lo = ofs[e], hi = ofs[e + 1];

    float w = 0.0f;
    for (int t = lo; t < hi; ++t) {
        const float pr = p_re[t], pi = p_im[t];
        // Compute the **increment** dP = (A-1)P + B E directly; do not write it as (A P + B E) - P.
        // For a high-Q pole |A| approaches 1 and P is larger than E by 1/(1-|A|), about 275 times
        // for the steepest silicon pole, so subtracting two large numbers loses roughly 2.5
        // significant digits in float32.
        //
        // For the record: this is genuinely better conditioned, but **measurement showed it is not
        // the source of the 4e-05 residual** in the equivalence check (4.048e-05 before, 4.042e-05
        // after). That residual turned out to be the transient of a finite DFT window. The form is
        // kept only because it costs nothing.
        const int ci = cidx[t];
        const float ar = lut_am1_re[ci], ai = lut_am1_im[ci];
        const float br = lut_b_re[ci],  bi = lut_b_im[ci];
        const float dr = ar * pr - ai * pi + br * en;
        const float di = ar * pi + ai * pr + bi * en;
        p_re[t] = pr + dr;
        p_im[t] = pi + di;
        w += 2.0f * dr;                 // the physical polarization is 2 Re(P)
    }
    // Each entry owns its (comp, cell) exclusively, so this is a plain write rather than an
    // accumulation and needs no atomics.
    if (c == 0)      hist_x[id] = w;
    else if (c == 1) hist_y[id] = w;
    else             hist_z[id] = w;
}

// Second pass: E^{n+1} is now known, so add B E^{n+1}.
//
// **The decay of P inside the absorbing layer also happens here.** Each step P is multiplied by
// the same decay profile as its host E. This keeps the decaying frame consistent: if E is in the
// decaying frame and P is not, the ADE inside the layer evolves a mismatched E/P system.
// Multiplying by the same profile is equivalent to damping the oscillators in the layer,
// dP/dt = (q-Gamma)P + rE, which is passive and entirely allowed by the physics of an absorber.
// Note this is **not** the fix for the slow divergence seen in a ring resonator; that one was bulk
// gain from a negative-weight pole. The factors are precomputed on the host: P hangs off an E
// component, so it takes the E-family profile value at that component's staggered position, with
// multiple slabs multiplying together in a corner.
//
// This used to be a separate kernel. It walked a sparse entry list through a second indirection
// (entry[i] -> ofs[e]) and then read and wrote the whole of p again, which on the largest scene
// meant 77.68M entries and 1.94 ms per step. Folded in here it costs only one extra dense factor
// array of length n_entry, whose value outside the layer is 1.0f. **Multiplying by 1.0f is a
// bitwise identity for any value**, so a scene with no absorber, or with all its dispersion
// outside the layer, produces exactly the same result. With has_dec=0 even the load is skipped,
// and the branch is warp-uniform.
extern "C" __global__ void dispersion_post(
    float* __restrict__ p_re, float* __restrict__ p_im,
    const float* __restrict__ lut_b_re, const float* __restrict__ lut_b_im,
    const int* __restrict__ cidx,
    const int* __restrict__ comp, const int* __restrict__ cell,
    const int* __restrict__ ofs,
    const float* __restrict__ Ex, const float* __restrict__ Ey,
    const float* __restrict__ Ez,
    const float* __restrict__ dec,   // (n_entry,) per-step decay factor; 1 outside the layer
    const int has_dec,               // 0 = no dispersive entries inside the layer, dec is not read
    const int n_entry)
{
    const int e = blockIdx.x * blockDim.x + threadIdx.x;
    if (e >= n_entry) return;

    const int c = comp[e], id = cell[e];
    const float en = pick(Ex, Ey, Ez, c, id);
    const float dc = has_dec ? dec[e] : 1.0f;
    for (int t = ofs[e]; t < ofs[e + 1]; ++t) {
        const int ci = cidx[t];
        p_re[t] = (p_re[t] + lut_b_re[ci] * en) * dc;
        p_im[t] = (p_im[t] + lut_b_im[ci] * en) * dc;
    }
}


// Fusion of post(n) and pre(n+1). Nothing modifies E on these cells between the two passes (the
// applicability test lives in _disp_fused on the solver side), so against the same en we first add
// B*E^{n+1} and then advance to P_half. The instruction sequence is exactly that of running them
// separately, hence **bitwise identical**, but P is read and written only once, saving 16 B per
// pole per step.
// Note this kernel carries **no** absorber decay: its applicability test already excludes
// dispersion inside the layer.
extern "C" __global__ void dispersion_step(
    float* __restrict__ p_re, float* __restrict__ p_im,
    float* __restrict__ hist_x, float* __restrict__ hist_y,
    float* __restrict__ hist_z,
    const float* __restrict__ lut_am1_re, const float* __restrict__ lut_am1_im,
    const float* __restrict__ lut_b_re, const float* __restrict__ lut_b_im,
    const int* __restrict__ cidx,
    const int* __restrict__ comp, const int* __restrict__ cell,
    const int* __restrict__ ofs,
    const float* __restrict__ Ex, const float* __restrict__ Ey,
    const float* __restrict__ Ez,
    const int n_entry)
{
    const int e = blockIdx.x * blockDim.x + threadIdx.x;
    if (e >= n_entry) return;

    const int c = comp[e], id = cell[e];
    const float en = pick(Ex, Ey, Ez, c, id);
    const int lo = ofs[e], hi = ofs[e + 1];

    float w = 0.0f;
    for (int t = lo; t < hi; ++t) {
        const int ci = cidx[t];
        const float br = lut_b_re[ci], bi = lut_b_im[ci];
        // post: P += B E^{n+1}. The dc that the original post multiplied by is always 1.0f here, a
        // bitwise identity, so it is omitted.
        const float pr = p_re[t] + br * en;
        const float pi = p_im[t] + bi * en;
        // pre(n+1): dP = (A-1)P + B E, the same expression as dispersion_pre
        const float ar = lut_am1_re[ci], ai = lut_am1_im[ci];
        const float dr = ar * pr - ai * pi + br * en;
        const float di = ar * pi + ai * pr + bi * en;
        p_re[t] = pr + dr;
        p_im[t] = pi + di;
        w += 2.0f * dr;
    }
    if (c == 0)      hist_x[id] = w;
    else if (c == 1) hist_y[id] = w;
    else             hist_z[id] = w;
}


// Fusion of post(n) and pre(n+1) for scenes with dispersion inside an absorbing layer. Beyond the
// plain applicability test, the only thing touching E between the two passes is the E-family
// absorber decay: at most three factors for that cell, in slab application order (f1/f2/f3, with
// 1.0f where that slab does not apply), multiplied in sequence. ((en*f1)*f2)*f3 is **bitwise
// identical** to absorb_slab multiplying once per slab, since times 1.0f is an identity. The decay
// of P still uses the dense dec, a single multiplication, the same expression as dispersion_post.
// Tensor scenes remain excluded, because the second tensor pass modifies E after post.
extern "C" __global__ void dispersion_step_dec(
    float* __restrict__ p_re, float* __restrict__ p_im,
    float* __restrict__ hist_x, float* __restrict__ hist_y,
    float* __restrict__ hist_z,
    const float* __restrict__ lut_am1_re, const float* __restrict__ lut_am1_im,
    const float* __restrict__ lut_b_re, const float* __restrict__ lut_b_im,
    const int* __restrict__ cidx,
    const int* __restrict__ comp, const int* __restrict__ cell,
    const int* __restrict__ ofs,
    const float* __restrict__ Ex, const float* __restrict__ Ey,
    const float* __restrict__ Ez,
    const float* __restrict__ dec,   // (n_entry,) per-step decay factor for P; 1 outside the layer
    const float* __restrict__ f1, const float* __restrict__ f2,
    const float* __restrict__ f3,    // E-family decay factors, in slab application order
    // Recomputing these three from the slab profiles was tried, saving 1.2 GB per step, and
    // measured 5% *slower* on the largest scene: its absorbing layers are 60/40/40 cells thick
    // with 77% of the cells inside them, so the six-iteration loop ran nearly full for almost
    // every entry and the saved traffic did not pay for the loop and register pressure.
    const int n_entry)
{
    const int e = blockIdx.x * blockDim.x + threadIdx.x;
    if (e >= n_entry) return;

    const int c = comp[e], id = cell[e];
    const float en = pick(Ex, Ey, Ez, c, id);
    const float en2 = ((en * f1[e]) * f2[e]) * f3[e];
    const float dc = dec[e];
    const int lo = ofs[e], hi = ofs[e + 1];

    float w = 0.0f;
    for (int t = lo; t < hi; ++t) {
        const int ci = cidx[t];
        const float br = lut_b_re[ci], bi = lut_b_im[ci];
        // post: P = (P + B E^{n+1}) * dec, the same expression as dispersion_post
        const float pr = (p_re[t] + br * en) * dc;
        const float pi = (p_im[t] + bi * en) * dc;
        // pre(n+1): by now the absorber has decayed E to en2
        const float ar = lut_am1_re[ci], ai = lut_am1_im[ci];
        const float dr = ar * pr - ai * pi + br * en2;
        const float di = ar * pi + ai * pr + bi * en2;
        p_re[t] = pr + dr;
        p_im[t] = pi + di;
        w += 2.0f * dr;
    }
    if (c == 0)      hist_x[id] = w;
    else if (c == 1) hist_y[id] = w;
    else             hist_z[id] = w;
}
// Variant that additionally writes en2 back into E, replacing the absorber's E pass when coverage
// is complete.
extern "C" __global__ void dispersion_step_dec_we(
    float* __restrict__ p_re, float* __restrict__ p_im,
    float* __restrict__ hist_x, float* __restrict__ hist_y,
    float* __restrict__ hist_z,
    const float* __restrict__ lut_am1_re, const float* __restrict__ lut_am1_im,
    const float* __restrict__ lut_b_re, const float* __restrict__ lut_b_im,
    const int* __restrict__ cidx,
    const int* __restrict__ comp, const int* __restrict__ cell,
    const int* __restrict__ ofs,
    float* __restrict__ Ex, float* __restrict__ Ey,
    float* __restrict__ Ez,          // write en2 back, replacing the absorber's E pass
    const float* __restrict__ dec,   // (n_entry,) per-step decay factor for P; 1 outside the layer
    const float* __restrict__ f1, const float* __restrict__ f2,
    const float* __restrict__ f3,    // E-family decay factors, in slab application order
    // Recomputing these three from the slab profiles was tried, saving 1.2 GB per step, and
    // measured 5% *slower* on the largest scene: its absorbing layers are 60/40/40 cells thick
    // with 77% of the cells inside them, so the six-iteration loop ran nearly full for almost
    // every entry and the saved traffic did not pay for the loop and register pressure.
    const int n_entry)
{
    const int e = blockIdx.x * blockDim.x + threadIdx.x;
    if (e >= n_entry) return;

    const int c = comp[e], id = cell[e];
    const float en = pick(Ex, Ey, Ez, c, id);
    const float en2 = ((en * f1[e]) * f2[e]) * f3[e];
    const float dc = dec[e];
    const int lo = ofs[e], hi = ofs[e + 1];

    float w = 0.0f;
    for (int t = lo; t < hi; ++t) {
        const int ci = cidx[t];
        const float br = lut_b_re[ci], bi = lut_b_im[ci];
        // post: P = (P + B E^{n+1}) * dec, the same expression as dispersion_post
        const float pr = (p_re[t] + br * en) * dc;
        const float pi = (p_im[t] + bi * en) * dc;
        // pre(n+1): by now the absorber has decayed E to en2
        const float ar = lut_am1_re[ci], ai = lut_am1_im[ci];
        const float dr = ar * pr - ai * pi + br * en2;
        const float di = ar * pi + ai * pr + bi * en2;
        p_re[t] = pr + dr;
        p_im[t] = pi + di;
        w += 2.0f * dr;
    }
    // Write the decayed E back, the same expression as the absorb_batch E pass. Each (comp, cell)
    // has exactly one entry (checked on the host), so the writes never conflict.
    if (c == 0)      { hist_x[id] = w; Ex[id] = en2; }
    else if (c == 1) { hist_y[id] = w; Ey[id] = en2; }
    else             { hist_z[id] = w; Ez[id] = en2; }
}


// ---- fast path for the two-pole bucket ----
// On the host, _pole_bucket_setup has already sorted the entries into buckets by pole count, in a
// stable ascending order (1-pole | 2-pole | 3-or-more). In the two-pole bucket [e0, e0+n) the slots
// run contiguously as slot0 + 2j with slot0 even, because the orphan slots ahead of the bucket are
// aligned; an orphan slot belongs to no entry range and is never read or written.
// cellc = comp<<30 | cell, packed. The per-pole arithmetic is in exactly the same order as the
// general kernel, so it is bitwise equivalent; float2 merges the stride-2 P and cidx accesses into
// coalesced loads and stores, measured at 1.35x on the benchmark.
// This kernel receives cellc (and dec/f) pointers already sliced at e0, with j starting from 0.

extern "C" __global__ void dispersion_step_p2(
    float* __restrict__ p_re, float* __restrict__ p_im,
    float* __restrict__ hist_x, float* __restrict__ hist_y,
    float* __restrict__ hist_z,
    const float* __restrict__ lut_am1_re, const float* __restrict__ lut_am1_im,
    const float* __restrict__ lut_b_re, const float* __restrict__ lut_b_im,
    const int* __restrict__ cidx,
    const int* __restrict__ cellc,
    const float* __restrict__ Ex, const float* __restrict__ Ey,
    const float* __restrict__ Ez,
    const int slot0, const int n)
{
    const int j = blockIdx.x * blockDim.x + threadIdx.x;
    if (j >= n) return;
    const int cc = cellc[j];
    const int c = ((unsigned)cc) >> 30, id = cc & 0x3FFFFFFF;
    const float en = pick(Ex, Ey, Ez, c, id);
    const int t0 = slot0 + 2 * j;
    float2 Pr = *reinterpret_cast<float2*>(p_re + t0);
    float2 Pi = *reinterpret_cast<float2*>(p_im + t0);
    const int2 ci = *reinterpret_cast<const int2*>(cidx + t0);
    float w = 0.0f;
    {
        const float br = lut_b_re[ci.x], bi = lut_b_im[ci.x];
        const float pr = Pr.x + br * en;
        const float pi = Pi.x + bi * en;
        const float ar = lut_am1_re[ci.x], ai = lut_am1_im[ci.x];
        const float dr = ar * pr - ai * pi + br * en;
        const float di = ar * pi + ai * pr + bi * en;
        Pr.x = pr + dr; Pi.x = pi + di; w += 2.0f * dr;
    }
    {
        const float br = lut_b_re[ci.y], bi = lut_b_im[ci.y];
        const float pr = Pr.y + br * en;
        const float pi = Pi.y + bi * en;
        const float ar = lut_am1_re[ci.y], ai = lut_am1_im[ci.y];
        const float dr = ar * pr - ai * pi + br * en;
        const float di = ar * pi + ai * pr + bi * en;
        Pr.y = pr + dr; Pi.y = pi + di; w += 2.0f * dr;
    }
    *reinterpret_cast<float2*>(p_re + t0) = Pr;
    *reinterpret_cast<float2*>(p_im + t0) = Pi;
    if (c == 0)      hist_x[id] = w;
    else if (c == 1) hist_y[id] = w;
    else             hist_z[id] = w;
}

extern "C" __global__ void dispersion_step_dec_p2(
    float* __restrict__ p_re, float* __restrict__ p_im,
    float* __restrict__ hist_x, float* __restrict__ hist_y,
    float* __restrict__ hist_z,
    const float* __restrict__ lut_am1_re, const float* __restrict__ lut_am1_im,
    const float* __restrict__ lut_b_re, const float* __restrict__ lut_b_im,
    const int* __restrict__ cidx,
    const int* __restrict__ cellc,
    const float* __restrict__ Ex, const float* __restrict__ Ey,
    const float* __restrict__ Ez,
    const float* __restrict__ dec,
    const float* __restrict__ f1, const float* __restrict__ f2,
    const float* __restrict__ f3,
    const int slot0, const int n)
{
    const int j = blockIdx.x * blockDim.x + threadIdx.x;
    if (j >= n) return;
    const int cc = cellc[j];
    const int c = ((unsigned)cc) >> 30, id = cc & 0x3FFFFFFF;
    const float en = pick(Ex, Ey, Ez, c, id);
    const float en2 = ((en * f1[j]) * f2[j]) * f3[j];
    const float dc = dec[j];
    const int t0 = slot0 + 2 * j;
    float2 Pr = *reinterpret_cast<float2*>(p_re + t0);
    float2 Pi = *reinterpret_cast<float2*>(p_im + t0);
    const int2 ci = *reinterpret_cast<const int2*>(cidx + t0);
    float w = 0.0f;
    {
        const float br = lut_b_re[ci.x], bi = lut_b_im[ci.x];
        const float pr = (Pr.x + br * en) * dc;
        const float pi = (Pi.x + bi * en) * dc;
        const float ar = lut_am1_re[ci.x], ai = lut_am1_im[ci.x];
        const float dr = ar * pr - ai * pi + br * en2;
        const float di = ar * pi + ai * pr + bi * en2;
        Pr.x = pr + dr; Pi.x = pi + di; w += 2.0f * dr;
    }
    {
        const float br = lut_b_re[ci.y], bi = lut_b_im[ci.y];
        const float pr = (Pr.y + br * en) * dc;
        const float pi = (Pi.y + bi * en) * dc;
        const float ar = lut_am1_re[ci.y], ai = lut_am1_im[ci.y];
        const float dr = ar * pr - ai * pi + br * en2;
        const float di = ar * pi + ai * pr + bi * en2;
        Pr.y = pr + dr; Pi.y = pi + di; w += 2.0f * dr;
    }
    *reinterpret_cast<float2*>(p_re + t0) = Pr;
    *reinterpret_cast<float2*>(p_im + t0) = Pi;
    if (c == 0)      hist_x[id] = w;
    else if (c == 1) hist_y[id] = w;
    else             hist_z[id] = w;
}

// Writeback variant: as above, but also writes en2 back into E, replacing the absorber's E pass
// when coverage is complete.
extern "C" __global__ void dispersion_step_dec_we_p2(
    float* __restrict__ p_re, float* __restrict__ p_im,
    float* __restrict__ hist_x, float* __restrict__ hist_y,
    float* __restrict__ hist_z,
    const float* __restrict__ lut_am1_re, const float* __restrict__ lut_am1_im,
    const float* __restrict__ lut_b_re, const float* __restrict__ lut_b_im,
    const int* __restrict__ cidx,
    const int* __restrict__ cellc,
    float* __restrict__ Ex, float* __restrict__ Ey, float* __restrict__ Ez,
    const float* __restrict__ dec,
    const float* __restrict__ f1, const float* __restrict__ f2,
    const float* __restrict__ f3,
    const int slot0, const int n)
{
    const int j = blockIdx.x * blockDim.x + threadIdx.x;
    if (j >= n) return;
    const int cc = cellc[j];
    const int c = ((unsigned)cc) >> 30, id = cc & 0x3FFFFFFF;
    const float en = pick(Ex, Ey, Ez, c, id);
    const float en2 = ((en * f1[j]) * f2[j]) * f3[j];
    const float dc = dec[j];
    const int t0 = slot0 + 2 * j;
    float2 Pr = *reinterpret_cast<float2*>(p_re + t0);
    float2 Pi = *reinterpret_cast<float2*>(p_im + t0);
    const int2 ci = *reinterpret_cast<const int2*>(cidx + t0);
    float w = 0.0f;
    {
        const float br = lut_b_re[ci.x], bi = lut_b_im[ci.x];
        const float pr = (Pr.x + br * en) * dc;
        const float pi = (Pi.x + bi * en) * dc;
        const float ar = lut_am1_re[ci.x], ai = lut_am1_im[ci.x];
        const float dr = ar * pr - ai * pi + br * en2;
        const float di = ar * pi + ai * pr + bi * en2;
        Pr.x = pr + dr; Pi.x = pi + di; w += 2.0f * dr;
    }
    {
        const float br = lut_b_re[ci.y], bi = lut_b_im[ci.y];
        const float pr = (Pr.y + br * en) * dc;
        const float pi = (Pi.y + bi * en) * dc;
        const float ar = lut_am1_re[ci.y], ai = lut_am1_im[ci.y];
        const float dr = ar * pr - ai * pi + br * en2;
        const float di = ar * pi + ai * pr + bi * en2;
        Pr.y = pr + dr; Pi.y = pi + di; w += 2.0f * dr;
    }
    *reinterpret_cast<float2*>(p_re + t0) = Pr;
    *reinterpret_cast<float2*>(p_im + t0) = Pi;
    if (c == 0)      { hist_x[id] = w; Ex[id] = en2; }
    else if (c == 1) { hist_y[id] = w; Ey[id] = en2; }
    else             { hist_z[id] = w; Ez[id] = en2; }
}
