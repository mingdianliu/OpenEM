# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Mingdian Liu
"""Boundary types, CPML coefficients, Absorber conductivity, and failing closed on boundary
topology.
"""

from __future__ import annotations

import numpy as np
import tidy3d as td

from openem import cpml, knobs
from openem.grid import ABSORBING, Axis, Grid, unravel_cells
from openem.scene._util import E_KEYS
from openem.units import UM

def boundary_types(sim: td.Simulation) -> list[tuple[str, str]]:
    """The boundary type names of each axis, as (low end, high end).

    The ``bloch_vec`` of a ``BlochBoundary`` is extracted and checked separately by
    :func:`bloch_vector`. For non-zero k the solver takes the complex-field path rather than
    failing closed.
    """
    return [
        (type(getattr(sim.boundary_spec, ax).minus).__name__,
         type(getattr(sim.boundary_spec, ax).plus).__name__)
        for ax in "xyz"
    ]


def bloch_vector(sim: td.Simulation) -> tuple[float, float, float]:
    """The Bloch wave vector of each axis, in Tidy3D's ``bloch_vec`` convention with units of 2*pi/L.

    A non-Bloch axis gets 0. Both ends must be a ``BlochBoundary`` with the same ``bloch_vec``:
    "phase on one end but not the other" has no physical meaning, and silently taking one of them
    would produce a solution at the wrong k point.
    """
    out = []
    for ax in "xyz":
        edge = getattr(sim.boundary_spec, ax)
        ks = [float(side.bloch_vec)
              for side in (edge.minus, edge.plus)
              if type(side).__name__ == "BlochBoundary"]
        if not ks:
            out.append(0.0)
            continue
        if len(ks) != 2 or ks[0] != ks[1]:
            raise NotImplementedError(
                f"the two BlochBoundary ends of axis {ax} disagree (bloch_vec={ks}); the phase of "
                "the periodic wraparound must be unique")
        out.append(ks[0])
    return tuple(out)


def _boundary(sim: td.Simulation, ax: int, side: str):
    """The tidy3d boundary object at the ``side`` end (``"lo"`` or ``"hi"``) of axis ``ax``."""
    edge = getattr(sim.boundary_spec, "xyz"[ax])
    return edge.minus if side == "lo" else edge.plus


def _pml_params(sim: td.Simulation, ax: int, side: str) -> cpml.PMLParams | None:
    """The PML parameters of that end; None when it is not a PML."""
    b = _boundary(sim, ax, side)
    if type(b).__name__ not in ("PML", "StablePML", "Absorber"):
        return None
    p = b.parameters
    # The `AbsorberParams` an `Absorber` uses has only the three sigma_* fields: it is a simple
    # absorbing layer made of a pure conductivity gradient, with no coordinate stretching and no
    # CFS. So kappa identically 1 and alpha identically 0 is its definition, not a guessed default.
    # The normalization of sigma follows the same formula as the PML's. **An Absorber does not go
    # through the CPML psi recursion**: its sigma is folded into the E update coefficients as an
    # ordinary conductivity; see :func:`fold_absorber_sigma`.
    if not hasattr(p, "kappa_min"):
        return cpml.PMLParams(
            num_layers=b.num_layers,
            sigma_min=p.sigma_min,
            sigma_max=p.sigma_max,
            sigma_order=p.sigma_order,
            kappa_min=1.0,
            kappa_max=1.0,
            kappa_order=1,
            alpha_min=0.0,
            alpha_max=0.0,
            alpha_order=1,
        )
    return cpml.PMLParams(
        num_layers=b.num_layers,
        sigma_min=p.sigma_min,
        sigma_max=p.sigma_max,
        sigma_order=p.sigma_order,
        kappa_min=p.kappa_min,
        kappa_max=p.kappa_max,
        kappa_order=p.kappa_order,
        alpha_min=p.alpha_min,
        alpha_max=p.alpha_max,
        alpha_order=p.alpha_order,
    )


def _avg_speed(eps: np.ndarray, ax: int, side: str, n_layers: int) -> float:
    """Mean wave speed inside the PML relative to c0, i.e. 1/n; Tidy3D's ``s_value`` uses it to
    scale sigma.
    """
    sl = [slice(None)] * 3
    sl[ax] = slice(0, n_layers) if side == "lo" else slice(-n_layers, None)
    return float(1.0 / np.sqrt(np.mean(eps[tuple(sl)])))


#: The boundary topologies the kernel actually implements. Each axis is either periodic or
#: absorbing, and the two ends must agree.
PERIODIC = ("Periodic", "BlochBoundary")


def refuse_unsupported_topology(grid: Grid) -> None:
    """Fail closed on any boundary topology the kernel cannot support.

    Each axis must be entirely periodic, or any combination of absorbing, PEC and PMC walls.
    ``kernels/yee.cu`` handles periodic and absorbing alike (a periodic axis simply fills its CPML
    coefficients with identity values), and PEC or PMC walls go through the index tables and masks
    (index_tables in grid.py). Periodic at one end and absorbing at the other is not implemented.

    The trap this comes from: an early kernel applied the PML only to z and hardcoded x and y as
    periodic, while `Scene` dutifully recorded PML on all six faces and computed the x and y
    coefficients too. Those were then silently discarded, the wave wrapped around through x and y,
    the field never decayed, and nothing raised anything.
    """
    # The test is by **category**, not by type name: PML at one end and Absorber at the other are
    # differently named but both absorbing, with identical index tables and masks, which is legal.
    # (Tidy3D itself rejects a PML plus Periodic mixture on the Boundary, so that never reaches
    # here.)
    # A PEC wall shares the masks of an absorbing end: at the low end pec[i=0]=0 pins tangential E
    # to zero (a PML is closed off by PEC anyway), and at the high end tangential E falls on edge n,
    # which is not stored, with mnx[n-1]=0. index_tables, the ghost convention
    # (GHOST_BY_BOUNDARY["PECBoundary"]=1) and _pml_params (which returns None and builds no psi)
    # are all already in place, so any axis combining absorbing and PEC ends is admitted directly.
    # A PMC wall works, but the high end carries a half-cell geometric error (measured; see
    # tests/test_pmc.py):
    #   the low wall is exact, since mpv[0]=-1 makes H(-1/2) = -H(+1/2) and H=0 lands exactly on
    #     node 0;
    #   the high wall is pulled in by dl/2, because tangential E on the wall is a **free** degree of
    #     freedom yet is not stored (drop-last), so all that can be done is to force the difference
    #     in the H update to zero, which puts dE_t/dn = 0 at the last H node z_{n-1/2} rather than
    #     at the wall, z_n.
    # The test: the |Ex|/max of the last point of a standing wave hits cos(pi/N) at both N=40 and
    # N=80, i.e. the wall pulled in by half a cell, and rules out cos(2*pi/N), the wall at its
    # nominal position. The consequence is an O(dl) geometric error: a PEC-plus-PMC cavity is
    # effectively dl/2 short. Fixing it needs one more stored column of tangential E at node n on
    # that axis, or a special case in the kernel.
    wall_ok = (*ABSORBING, "PECBoundary", "PMCBoundary")
    bad = []
    for ax, a in zip("xyz", grid.axes):
        ends = (a.boundary_lo, a.boundary_hi)
        if all(e in PERIODIC for e in ends) or all(e in wall_ok for e in ends):
            continue
        bad.append(f"axis {ax} is {a.boundary_lo} / {a.boundary_hi}")
    if bad:
        raise NotImplementedError(
            "unsupported boundary topology: " + "; ".join(bad)
            + f". Each axis must be periodic ({'/'.join(PERIODIC)}), or a combination of "
            f"absorbing and PEC/PMC walls ({'/'.join(ABSORBING)}/PECBoundary/PMCBoundary)"
        )


#: Former name.
refuse_unsupported = refuse_unsupported_topology



def build_grid(sim: td.Simulation) -> Grid:
    """``sim.grid.boundaries`` -> :class:`~openem.grid.Grid`, converting units to metres on the way.

    Raises:
        AssertionError: the cell count we derive disagrees with ``sim.grid.num_cells``, which means
            Tidy3D's grid convention has been misunderstood and we must not continue.
    """
    b = sim.grid.boundaries
    bnd = boundary_types(sim)
    flat = [float(s) == 0.0 for s in sim.size]      # in Tidy3D only size=0 marks a flat 2D axis
    grid = Grid(
        x=Axis(np.asarray(b.x, dtype=np.float64) * UM, *bnd[0], flat=flat[0]),
        y=Axis(np.asarray(b.y, dtype=np.float64) * UM, *bnd[1], flat=flat[1]),
        z=Axis(np.asarray(b.z, dtype=np.float64) * UM, *bnd[2], flat=flat[2]),
    )
    refuse_unsupported_topology(grid)
    if tuple(grid.shape) != tuple(sim.grid.num_cells):
        raise AssertionError(
            f"cell count mismatch: computed {grid.shape} locally, tidy3d reports {tuple(sim.grid.num_cells)}"
        )
    return grid


def _is_absorber(sim: td.Simulation, ax: int, side: str) -> bool:
    return type(_boundary(sim, ax, side)).__name__ == "Absorber"


def _faces(sim: td.Simulation, grid: Grid, want: str):
    """Yield the absorbing faces one at a time as ``(ax, side, params, dl)``, ordered by axis and
    then lo before hi.

    ``want="psi"`` yields the psi-PML faces (PML and StablePML; an Absorber does not count);
    ``want="absorber"`` yields the Absorber faces, raising ValueError when the layer count exceeds
    the total number of cells on that axis.
    ``dl`` is the representative step at that end; Tidy3D's ``create_sfactor_*`` likewise passes a
    single value, ``dls[0]`` at the min end and ``dls[-1]`` at the max end.
    """
    for ax in range(3):
        axis = grid.axes[ax]
        for side in ("lo", "hi"):
            params = _pml_params(sim, ax, side)
            if params is None or _is_absorber(sim, ax, side) != (want == "absorber"):
                continue
            if want == "absorber" and int(params.num_layers) > axis.n:
                raise ValueError(
                    f"axis {'xyz'[ax]}, {side} end: the Absorber has "
                    f"{int(params.num_layers)} layers, more than the {axis.n} cells on that axis")
            dl = float(axis.dl[0] if side == "lo" else axis.dl[-1])
            yield ax, side, params, dl


def absorber_bands(sim: td.Simulation, grid: Grid) -> list[tuple[int, int, int]]:
    """``(axis, starting global cell, layer count)`` for each Absorber face."""
    out = []
    for ax, side, params, _dl in _faces(sim, grid, "absorber"):
        nlay = int(params.num_layers)
        out.append((ax, 0 if side == "lo" else grid.axes[ax].n - nlay, nlay))
    return out


def fold_absorber_sigma(sim: td.Simulation, grid: Grid, eps: dict,
                        sigma: dict) -> dict:
    """Fold an Absorber's graded conductivity into the sigma arrays of the three E components, in
    S/m, returning the same dict.

    Tidy3D's ``Absorber`` is a graded lossy layer with **conductivity only and no magnetic loss**,
    described as a multilayer system with gradually increasing conductivity, closed off by PEC.
    The normalization of sigma uses the same formula as the PML::

        σ(x) = sigma_max · avg_speed / (η₀ · dl) · (x/d)^sigma_order

    ``avg_speed = 1/sqrt(mean eps)`` over the layer band, and ``dl`` is the representative step at
    that end.

    Evidence: the reference reflectance of an absorbing boundary, as two curves against layer count
    and against resolution, agrees point by point with a 1D transfer-matrix model that is E-only
    under this same normalization. The older matched implementation, decaying E and H by the same
    Gamma, came out about 10 times stronger than the reference.

    Folding it into sigma rather than a separate decay kernel: sigma goes through ``ca/cb``
    (coeffs._e_coeffs) so each cell uses **its own eps**. A waveguide core at eps=12 and cladding at
    2.1 differ in decay rate by a factor of 6, and a slab-averaged eps would be wrong. Dispersive
    cells are handled by the implicit coupling on the same path and are unconditionally stable.

    Staggering: along the absorbing axis the normal component (``comp == axis``) sits at the half
    cell and samples depth at the H positions, while the tangential ones sit at the whole cell and
    sample at the E positions, as in ``cpml.steps``.
    """
    for ax, side, params, dl in _faces(sim, grid, "absorber"):
        nlay = int(params.num_layers)
        g0 = 0 if side == "lo" else grid.axes[ax].n - nlay
        speed = _avg_speed(eps["ez"], ax, side, nlay)
        prof = {at: cpml.profile(params, dl, side, at, speed)[0]  # type: ignore[arg-type]
                for at in ("E", "H")}
        shape = [1, 1, 1]
        shape[ax] = nlay
        sl = [slice(None)] * 3
        sl[ax] = slice(g0, g0 + nlay)
        for c, key in enumerate(E_KEYS):
            if sigma[key] is None:
                sigma[key] = np.zeros(grid.shape, dtype=np.float64)
            sigma[key][tuple(sl)] += prof["H" if c == ax else "E"].reshape(shape)
    return sigma


#: The CFS alpha coefficient injected into the psi recursion when the PML band contains cells with
#: ADE poles. Dimensionless, normalized like sigma: ``alpha(layer) = coeff * avg_speed/(eta0*dl) *
#: step``, on the same depth profile as sigma at order 1.
#:
#: **Why it is injected**: a CPML psi recursion with alpha=0 diverges at the **scheme level** over a
#: long run against ADE poles in the same cells. A term-by-term whole-domain CPU replay of one case
#: (identical trajectories in fp32 and fp64) measured exponential growth of about 6.7e-4 per step,
#: with the mode sitting on the silicon rows inside the x PML band at an oscillation frequency of
#: omega*dt about 0.0084. That is one seventh of the working band at 0.059, exactly the
#: pathological region where 1/s(omega -> 0) -> 0 at alpha=0, and removing the pole entries inside
#: the PML band made it entirely stable. Sweeping the pole damping gamma from 0 to 1e-2 barely
#: changed the growth rate, because the low-frequency region uses the **quasi-static** response of
#: the pole rather than its resonance; kappa identically 1 and a raw-alpha StablePML did not help
#: either. A CFS alpha clamps the low-frequency stretch to a finite real number (s -> kappa +
#: sigma/alpha), the same mechanism as tidy3d's own StablePML "for dispersive materials". But our
#: profile copies tidy3d's s_value faithfully, where alpha is not normalized and stays in raw S/m,
#: so 0.9 S/m amounts to nothing; here it is converted to an effective value under the sigma
#: normalization.
#:
#: **Magnitude**: 0.05 gives alpha_max about 3.7e3 S/m on that grid, screening omega below
#: alpha/eps0 = 4.2e14 rad/s (omega*dt about 0.02, a factor of 2.4 above the diverging mode at
#: 0.0084). At the lower edge of the working band alpha/(omega*eps0) is at most 0.35, and since it
#: grows as step to the first power it only acts deep in the layer, perturbing the in-band
#: absorption by a few percent, with the residual PML reflection still far below any observable
#: tolerance. Both 0.05 and 0.9 were entirely stable over an 18k-step replay, more than an order of
#: magnitude of margin.
#:
#: **Gating**: injected only into a face whose layer band genuinely contains pole cells and whose
#: user-supplied alpha_max is 0. A non-dispersive scene, a scene whose poles stay out of the PML,
#: and any face where the user gave an explicit CFS are all bitwise unchanged.
#: ``OPENEM_PML_DISP_ALPHA`` overrides it; setting 0 restores the old behaviour.
PML_DISP_ALPHA = 0.05


def _pml_disp_alpha() -> float:
    return float(knobs.env("PML_DISP_ALPHA"))     # default registered in knobs.KNOBS, equal to PML_DISP_ALPHA (locked by tests/test_knobs_scene.py)


def _band_has_poles(disp_cells: np.ndarray | None, grid: Grid, ax: int,
                    side: str, n_layers: int) -> bool:
    """Whether the layer band of this PML face contains any cell carrying poles."""
    if disp_cells is None or disp_cells.size == 0 or n_layers <= 0:
        return False
    cells = np.asarray(disp_cells, dtype=np.int64)
    coord = unravel_cells(cells, grid.shape)[ax]
    n_ax = grid.shape[ax]
    if side == "lo":
        return bool(np.any(coord < n_layers))
    return bool(np.any(coord >= n_ax - n_layers))


def _cfs_coeffs(params: cpml.PMLParams, dl: float, side: str, dt: float,
                speed: float, acoef: float) -> cpml.PMLCoeffs:
    """psi recursion coefficients after injecting a CFS alpha under the same normalization as sigma;
    sigma and kappa are untouched.

    The only difference from :func:`cpml.build` is the alpha profile,
    ``acoef * speed/(eta0*dl) * step`` at order 1. cpml.profile stays bitwise faithful to tidy3d's
    ``s_value``, which is why this variant lives here rather than in the cpml module.
    """
    out = {}
    for at in ("E", "H"):
        sigma, kappa, _alpha = cpml.profile(params, dl, side, at, speed)  # type: ignore[arg-type]
        alpha = acoef * speed / (cpml.ETA_0 * dl) * cpml.steps(
            params.num_layers, side, at)  # type: ignore[arg-type]
        a, b = cpml.recursion_coeffs(sigma, kappa, alpha, dt)
        out[f"a_{at}"] = a
        out[f"b_{at}"] = b
        out[f"inv_kappa_{at}"] = 1.0 / kappa
    return cpml.PMLCoeffs(**out)  # type: ignore[arg-type]


def build_pml(
    sim: td.Simulation, grid: Grid, eps_normal: np.ndarray,
    disp_cells: np.ndarray | None = None,
) -> dict[tuple[int, str], cpml.PMLCoeffs]:
    """Recursion coefficients for each psi-PML face.

    ``disp_cells`` is the set of flat cell indices carrying ADE poles anywhere in the scene
    (Dispersion together with DispersionMix). A psi-PML face whose layer band contains pole cells
    and whose user alpha_max is 0 gets a CFS alpha injected (:data:`PML_DISP_ALPHA`, the fix for
    the scheme-level instability); every other face is bitwise unchanged.

    ``dl`` is the representative step at that end (see :func:`_faces`).

    **Absorber faces are not handled here**: an Absorber is a graded lossy medium rather than a
    coordinate stretch, and its sigma is folded into the E update coefficients by
    :func:`fold_absorber_sigma`.
    """
    pml: dict[tuple[int, str], cpml.PMLCoeffs] = {}
    for ax, side, params, dl in _faces(sim, grid, "psi"):
        speed = _avg_speed(eps_normal, ax, side, params.num_layers)
        acoef = _pml_disp_alpha()
        if (acoef != 0.0 and params.alpha_max == 0.0 and params.alpha_min == 0.0
                and _band_has_poles(disp_cells, grid, ax, side,
                                    params.num_layers)):
            # Pole cells inside the band, so inject a CFS alpha; this is the scheme-level
            # instability between an alpha=0 psi recursion and the ADE
            pml[(ax, side)] = _cfs_coeffs(
                params, dl, side, float(sim.dt), speed, acoef)
            continue
        pml[(ax, side)] = cpml.build(
            params, dl=dl, side=side, dt=sim.dt, avg_speed=speed)  # type: ignore[arg-type]
    return pml


#: Depth-graded artificial damping gamma_max applied to ADE poles inside a psi-PML band.
#: Dimensionless, relative to |q|: zero at the interface and growing outward as step cubed, the
#: same order as the sigma profile.
#:
#: **Why it was added**: the CFS alpha injection cured the scheme-level family of alpha=0 CPML
#: against ADE (about 6.7e-4 per step) on the whole-domain CPU replay, but GPU probes showed a
#: faster residual mode of the same family still present on the complete solve path: one
#: configuration went NaN at about 36k steps and two others at about 13.5k, implying growth rates
#: of 3e-3 to 1e-2 per step. That residual mode is not reproducible on CPU, even though the replay
#: was aligned term by term with yee.cu and yee_ade.cu, including the PEC closure, the P43 fusion,
#: the LUT and the source injection timing, and stayed entirely stable. But it is sensitive to
#: alpha (one configuration was cured by alpha alone), so it belongs to the same psi-by-pole
#: coupling inside the band. At alpha=0.05 the psi-side damping mid-layer is about 1e-2 per step,
#: the same order as the residual growth rate, which leaves no margin. This mechanism adds a second
#: line from the **pole side**: gamma_d(layer) = gamma_max * step^3 lowers the pole spectral radius
#: deep in the layer by gamma_d*|q|*dt of about 0.034 per step (gamma_max=0.3, the dominant silicon
#: pole at omega*dt about 0.9, step=0.5), an order of magnitude above the implied growth rate.
#: Controlled experiments had already shown this mechanism alone cures the CPU-reproducible family;
#: it is now productized and stacked with alpha.
#:
#: **Cost**: only the pole coefficients of cells inside a psi-PML band change, and the interior of
#: a PML is not a physical region. At the interface layer gamma_d = gamma_max/1728, about 1.7e-4,
#: perturbing in-band reflection far less than the discretization error of the sigma profile
#: itself. Absorber faces do not count as a band; they have their own decaying-frame consistency
#: mechanism (see kernels/absorber.cu).
#:
#: **Gating**: no dispersion, no psi-PML face, poles outside the band, or the environment variable
#: set to 0, all leave everything bitwise unchanged. ``OPENEM_PML_POLE_DAMP`` overrides it.
PML_POLE_DAMP = 0.3


def _pml_pole_damp() -> float:
    return float(knobs.env("PML_POLE_DAMP"))      # default registered in knobs.KNOBS, equal to PML_POLE_DAMP


def _pml_depth_profiles(sim: td.Simulation, grid: Grid) -> list[np.ndarray] | None:
    """Normalized depth profile of the psi-PML on each axis, zero outside the band and ignoring
    Absorber faces. None when there is none at all.
    """
    prof = [np.zeros(n_ax) for n_ax in grid.shape]
    any_face = False
    for ax, side, params, _dl in _faces(sim, grid, "psi"):
        n_ax, v = grid.shape[ax], prof[ax]
        m = min(int(params.num_layers), n_ax)
        if m <= 0:
            continue
        any_face = True
        i = np.arange(m, dtype=np.float64)
        if side == "lo":
            v[:m] = np.maximum(v[:m], (m - i) / m)
        else:
            v[n_ax - m:] = np.maximum(v[n_ax - m:], (i + 1) / m)
    return prof if any_face else None


def _damped_pole_arrays(am1: np.ndarray, b: np.ndarray, gd: np.ndarray,
                        dt: float) -> tuple[np.ndarray, np.ndarray]:
    """Shift the real part of a whole pole array by a per-pole gamma_d; entries with gamma_d=0 are
    bitwise untouched.

    ``(A-1, B)`` is modified directly in the discrete domain, by the same formula as the controlled
    experiment:
      q = (2/dt) * am1/(am1+2), the inverse of the trapezoidal map; q' = q - gamma_d*|q|;
      am1' = (1 + q'*dt/2)/(1 - q'*dt/2) - 1; b' = b * (1 - q*dt/2)/(1 - q'*dt/2).
    """
    m = gd > 0.0
    if not np.any(m):
        return am1, b
    am1 = am1.copy()
    b = b.copy()
    a_ = am1[m]
    with np.errstate(divide="ignore", invalid="ignore"):
        q = (2.0 / dt) * a_ / (a_ + 2.0)
        q2 = q - gd[m] * np.abs(q)
        den0 = 1.0 - q * dt / 2.0
        den2 = 1.0 - q2 * dt / 2.0
        a_new = (1.0 + q2 * dt / 2.0) / den2 - 1.0
        b_new = b[m] * den0 / den2
    # The inverse map is ill-conditioned when A is near -1 (omega*dt approaching pi), so that pole
    # is left untouched, failing closed to the old behaviour rather than writing an untrustworthy
    # coefficient
    ok = np.isfinite(a_new) & np.isfinite(b_new)
    am1[m] = np.where(ok, a_new, a_)
    b[m] = np.where(ok, b_new, b[m])
    return am1, b


def _per_pole_gamma(gd_e: np.ndarray, ofs: np.ndarray) -> np.ndarray:
    """Expand a per-entry gamma_d onto the per-pole slots; the CSR offsets ``ofs`` give the pole
    count of each entry.
    """
    return np.repeat(gd_e, np.diff(ofs).astype(np.int64))


def damp_pml_poles(sim: td.Simulation, grid: Grid, disp, disp_mix, dt: float):
    """Inject depth-graded damping (:data:`PML_POLE_DAMP`) into the poles inside a psi-PML band,
    returning a new table.

    Entries outside the band, gamma_max=0, no psi-PML face, and no dispersion all leave everything
    bitwise unchanged. ``Dispersion.g`` (= 2*sum Re(B)) is recomputed for the entries that changed.
    """
    import dataclasses

    gmax = _pml_pole_damp()
    if gmax == 0.0 or disp is None:
        return disp, disp_mix
    prof = _pml_depth_profiles(sim, grid)
    if prof is None:
        return disp, disp_mix
    nx, ny, nz = grid.shape
    # gamma_d per cell: take the max depth across axes, then cube it
    depth = np.maximum.reduce([
        prof[0][:, None, None] * np.ones((1, ny, nz)),
        prof[1][None, :, None] * np.ones((nx, 1, nz)),
        prof[2][None, None, :] * np.ones((nx, ny, 1))])
    gd_cell = (gmax * depth ** 3).ravel()

    if disp.n_entry:
        gd_e = gd_cell[disp.cell.astype(np.int64)]
        if np.any(gd_e > 0.0):
            counts = np.diff(disp.pole_ofs).astype(np.int64)
            am1, b = _damped_pole_arrays(disp.am1, disp.b, _per_pole_gamma(gd_e, disp.pole_ofs), dt)
            g = disp.g.copy()
            me = gd_e > 0.0
            # g = 2 sum Re(B): recompute only the entries that changed
            seg = np.add.reduceat(2.0 * b.real, disp.pole_ofs[:-1])
            seg = np.where(counts > 0, seg, 0.0)
            g[me] = seg[me]
            disp = dataclasses.replace(disp, am1=am1, b=b, g=g)

    if disp_mix is not None and disp_mix.n_entry:
        gd_e = gd_cell[disp_mix.cell.astype(np.int64)]
        if np.any(gd_e > 0.0):
            pa, pb = _damped_pole_arrays(
                disp_mix.pa, disp_mix.pb, _per_pole_gamma(gd_e, disp_mix.p_ofs), dt)
            qa, qb = _damped_pole_arrays(
                disp_mix.qa, disp_mix.qb, _per_pole_gamma(gd_e, disp_mix.q_ofs), dt)
            disp_mix = dataclasses.replace(disp_mix, pa=pa, pb=pb, qa=qa, qb=qb)

    return disp, disp_mix

