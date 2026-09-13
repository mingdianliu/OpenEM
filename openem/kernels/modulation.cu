// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Mingdian Liu
// Time-modulated media (modulation_spec): each step, replace ca/cb of the modulated cells with
// that step's values.
//
// Physics and discretization. A modulated cell is lossless and non-dispersive (the build side
// fails closed otherwise), with the constitutive relation D = eps0 * eps(t) * E and
// eps(t) = eps_s + delta_eps(t), where
//
//   delta_eps(t) = amp * cos(omega*t - phi)
//
// The convention is taken from the tidy3d client, components/time_modulation.py:
//   delta_eps = Re[A_t e^{i phi_t - i omega t} * A_r e^{i phi_r}],
//   amp = A_t * A_r, phi = phi_t + phi_r.
//
// Discretized in **D form**, the standard charge-conserving approach: advance D first, then divide
// by that step's eps:
//
//   eps0 * eps^{n+1} * E^{n+1} = eps0 * eps^n * E^n + dt * curl H^{n+1/2}
//   => E^{n+1} = [eps^n/eps^{n+1}] E^n + [(dt/eps0)/eps^{n+1}] * curl H
//
// That is exactly the shape of update_e's ca/cb (ca = eps^n/eps^{n+1},
// cb = (dt/eps0)/eps^{n+1}), so **update_e needs no change**: this kernel runs just before
// update_e each step and rewrites only the modulated cells' ca/cb. A scene with no modulation
// never launches it at all and stays bitwise unchanged, on a single code path, the same way the
// absorber does it.
//
// The trigonometry is not done in the kernel. Each step the host computes two pairs of scalars in
// float64, (cos omega*t, sin omega*t) at t_n and at t_{n+1}, and passes them in. A cell's own
// phase is expanded with the angle-sum formula:
//   cos(omega*t - phi) = cos(omega*t) * cos(phi) + sin(omega*t) * sin(phi)
// with cos(phi) and sin(phi) precomputed into a table at build time. This avoids the precision
// loss of a float32 large-angle remainder, and spares the common globally-in-phase case from
// paying for a per-cell cosine.
//
// Storage uses the same packed layout as dispersion.cu: comp[e] and cell[e], with each
// (comp, cell) appearing exactly once, locked by validate on the build side, so the writes need no
// atomics.

extern "C" __global__ void modulate_coeffs(
    float* __restrict__ cax, float* __restrict__ cay, float* __restrict__ caz,
    float* __restrict__ cbx, float* __restrict__ cby, float* __restrict__ cbz,
    const float* __restrict__ amp,        // A_t * A_r
    const float* __restrict__ cph,        // cos(phi)
    const float* __restrict__ sph,        // sin(phi)
    const float* __restrict__ eps_stat,   // static eps_s, subpixel averaging included
    const int* __restrict__ comp, const int* __restrict__ cell,
    const float* __restrict__ cs_tab,     // (n_total, 4) = [c0, s0, c1, s1], computed by the host
                                          // in float64 then narrowed to f32
    const float* __restrict__ sd,         // s = sigma*dt/(2*eps0); all zero when lossless
    const float* __restrict__ gg,          // implicit dispersion coupling G; all zero when non-dispersive
    const int* __restrict__ step,
    const float dte,                      // dt / ε₀
    const int n_entry)
{
    const int e = blockIdx.x * blockDim.x + threadIdx.x;
    if (e >= n_entry) return;

    const float c0 = cs_tab[4 * (*step) + 0], s0 = cs_tab[4 * (*step) + 1];
    const float c1 = cs_tab[4 * (*step) + 2], s1 = cs_tab[4 * (*step) + 3];
    const int c = comp[e], id = cell[e];
    const float a = amp[e], cp = cph[e], sp = sph[e];
    const float e0 = eps_stat[e] + a * (c0 * cp + s0 * sp);   // eps(t_n)
    const float e1 = eps_stat[e] + a * (c1 * cp + s1 * sp);   // eps(t_{n+1})
    // Combined update form: when dispersion or loss coexists with modulation, the numerator takes
    // eps_inf at step n and the denominator at step n+1. With s = G = 0 this collapses back to
    // exactly e0/e1, since subtracting and adding zero are exact in IEEE754, so a lossless
    // non-dispersive scene is bitwise identical to before this was added.
    const float s = sd[e], g = gg[e];
    const float den = e1 + s + g;
    const float ca = (e0 - s) / den;
    const float cb = dte / den;
    if (c == 0)      { cax[id] = ca; cbx[id] = cb; }
    else if (c == 1) { cay[id] = ca; cby[id] = cb; }
    else             { caz[id] = ca; cbz[id] = cb; }
}
