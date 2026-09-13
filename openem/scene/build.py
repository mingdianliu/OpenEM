# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Mingdian Liu
"""Assemble media, boundaries, sources and monitors into one :class:`Scene`.

This module only does the assembly; the conventions and the fail-closed checks of each kind of
physics live in their own modules.
"""

from __future__ import annotations

from typing import NamedTuple

import numpy as np
import tidy3d as td

from openem.grid import unravel_cells
from openem.model import Scene
from openem.scene import boundaries, dispersion, media, modulation, monitors, sources
from openem.scene.poles import any_dispersive


def _inside_tfsf_box(ii, jj, kk, t) -> np.ndarray:
    """Per-cell test of whether a cell sits at least 2 cells inside the TFSF box (clear of the
    correction band on the box faces); a bool array."""
    (ilo, jlo, klo), (ihi, jhi, khi) = t.box_lo, t.box_hi
    return ((ii >= ilo + 2) & (ii <= ihi - 3) & (jj >= jlo + 2)
            & (jj <= jhi - 3) & (kk >= klo + 2) & (kk <= khi - 3))


def _check_pec_placement(pec, shape, srcs) -> None:
    """Fail closed on every PEC placement whose convention is not settled, or where an error would
    be hidden silently.

    - **An injection plane crossing PEC cells** (plane-wave or mode source): the injection on those
      cells is quietly swallowed by cb=0, and the reference implementation's behavior for such a
      placement has never been verified.
    - **TFSF**: the PEC has to sit at least 2 cells inside the box. A PEC outside the box (the
      scattered-field region) is physically wrong to begin with, since what it constrains should be
      the total field and not the scattered field; a PEC touching the box faces or the box columns
      makes the ε uniformity check read values the mask had already scrubbed to 1.0 and pass
      silently.
    - **PermittivityMonitor**: the ε of masked cells has been refilled with 1.0, so reporting it
      would mislead.
    """
    cells = np.unique(np.concatenate([pec[k] for k in ("ex", "ey", "ez")]))
    if cells.size == 0:
        return
    ii, jj, kk = unravel_cells(cells, shape)
    for src in srcs.planes:
        if np.any(kk == src.plane_index):
            raise NotImplementedError(
                f"the injection plane k={src.plane_index} of a plane-wave source crosses PEC "
                "cells; the convention is not settled")
    for m in srcs.mode_srcs:
        n_hit = int(np.count_nonzero(ii == m.plane_index))
        if n_hit:
            # Let through since 2026-09-04 (FieldProjections: the mode source plane crosses a PEC
            # screen). PEC cells have ca=cb=0 so the injection is swallowed, but the mode profile
            # itself is ≈0 inside PEC (the mode solver treats PEC as a boundary too), so what gets
            # swallowed was zero anyway. The remaining risk is a half-cell error in the first cell
            # against the PEC face, the same order as the accuracy concession staircased PEC
            # already makes.
            print(f"[openem] the injection plane i={m.plane_index} of a mode source crosses "
                  f"{n_hit} PEC cells; injection on PEC cells is ignored as E≡0")
    for t in srcs.tfsf_srcs:
        inside = _inside_tfsf_box(ii, jj, kk, t)
        if not bool(inside.all()):
            (ilo, jlo, klo), (ihi, jhi, khi) = t.box_lo, t.box_hi
            raise NotImplementedError(
                f"{int((~inside).sum())} PEC cells are not at least 2 cells inside the TFSF box "
                f"(box [{ilo},{ihi}]x[{jlo},{jhi}]x[{klo},{khi}])")
    # PermittivityMonitor coexisting with PEC: the ε array of masked cells is refilled with 1.0
    # (the kernel zeroes them from the pec mask and never looks at ε), and the readout side
    # (nb.backend._perm_data) uses the pec_ex/ey/ez masks of the Scene to restore those cells to
    # tidy3d's convention, constants.pec_val = -1e8 (measured on the client: sampling ε on a PEC
    # cell always gives -1e8, with no intermediate value).


def _check_tensor_placement(tensor, shape, srcs, absorbers) -> None:
    """Fail closed on every combination of tensor cells with a source or an absorber whose
    convention has never been settled.

    - **An injection plane crossing tensor cells** (plane-wave or mode source): the cb the
      injection is multiplied by is 1 on pass-through cells, so the whole amplitude convention is
      wrong.
    - **TFSF**: the tensor cells (halo included) all have to sit at least 2 cells inside the box,
      since the face correction is multiplied by cb as well and the 1D incident background assumes
      the box shell is a scalar medium.
    - **Absorber**: the convention for stacking the in-layer attenuation on top of a pass-through
      value has never been verified.
    """
    ii, jj, kk = unravel_cells(tensor.cell, shape)
    coords = (ii, jj, kk)
    for src in srcs.planes:
        if np.any(coords[src.axis] == src.plane_index):
            raise NotImplementedError(
                f"the injection plane {'xyz'[src.axis]}={src.plane_index} of a plane-wave source "
                "crosses tensor cells; the convention is not settled")
    for m in srcs.mode_srcs:
        if np.any(ii == m.plane_index):
            raise NotImplementedError(
                f"the injection plane i={m.plane_index} of a mode source crosses tensor cells; "
                "the convention is not settled")
    for t in srcs.tfsf_srcs:
        inside = _inside_tfsf_box(ii, jj, kk, t)
        if not bool(inside.all()):
            raise NotImplementedError(
                f"{int((~inside).sum())} tensor cells (halo included) are not at least 2 cells "
                "inside the TFSF box")
    for ax, g0, nlay in absorbers:
        if np.any((coords[ax] >= g0) & (coords[ax] < g0 + nlay)):
            raise NotImplementedError(
                f"tensor cells fall inside the Absorber layer of the {'xyz'[ax]} axis; the "
                "convention is not settled")


class Materials(NamedTuple):
    """The seven things the material side builds. ``pec``, ``disp``, ``disp_mix``, ``mod`` and
    ``tensor`` are None when the scene has nothing of that kind."""

    eps: dict
    sigma: dict
    pec: dict | None
    disp: object
    disp_mix: object
    mod: object
    tensor: object


def _materials(sim: td.Simulation, grid) -> Materials:
    """The :class:`Materials` of the material side.

    A dispersive scene takes a different extraction path: eps stores ε_∞ and the poles come out
    separately (scene/dispersion.py). The non-dispersive scene keeps its original path untouched;
    it was verified against 9 reference cases, so swapping the implementation would gain nothing
    and only risk a regression. The kernel side has a single path, locked down by the
    bitwise-identical test in test_dispersion. The three unsupported combinations fail closed here.
    """
    disp = disp_mix = pec = mod = tensor = None
    # Each of these booleans is computed once: every call has to walk sim.volumetric_structures (a
    # plain tidy3d property that re-converts Medium2D each time)
    has_tensor = media.has_tensor(sim)
    is_disp = any_dispersive(sim)
    with_pec = media.has_pec(sim)
    is_mod = modulation.any_modulated(sim)
    if has_tensor:
        # The three pass-through tracks (dispersion mix, tensor) cannot occupy the same cell; the
        # ordering of the PEC mask against the tensor, and how modulation rewrites ca/cb on
        # pass-through cells, have no settled convention either. Unverified means fail closed.
        if is_disp:
            raise NotImplementedError(
                "FullyAnisotropicMedium together with a dispersive medium is not supported yet")
        if with_pec:
            raise NotImplementedError(
                "FullyAnisotropicMedium together with PEC is not supported yet")
        if is_mod:
            raise NotImplementedError(
                "FullyAnisotropicMedium together with a time-varying medium is not supported yet")
    if is_disp:
        if with_pec:
            raise NotImplementedError(
                "PEC together with a dispersive medium is not supported yet: the dispersion "
                "extraction path has no PEC mask")
        eps, sigma, disp, disp_mix = dispersion.build(sim, grid)
        # Inject depth-graded damping into the poles inside the ψ-PML band (a second margin on top
        # of the CFS α, boundaries.PML_POLE_DAMP; outside the band, and at γ=0, bitwise unchanged)
        disp, disp_mix = boundaries.damp_pml_poles(
            sim, grid, disp, disp_mix, float(sim.dt))
        if is_mod:
            # Dispersive + time-varying: the modulation rides on ε_∞. What dispersion.build
            # returns is exactly ε_∞ at the three Yee component positions, and ca/cb are
            # recomputed every step by kmod; the recurrence coefficients of the ADE depend only on
            # the pole parameters and dt, are unaffected by the modulation, and stay as they are.
            mod = modulation.build(sim, grid, eps, sigma)
    else:
        eps, sigma, pec = media.epsilon_and_sigma(
            media.tensor_free_sim(sim) if has_tensor else sim)
        if is_mod:
            if pec is not None:
                raise NotImplementedError(
                    "a time-varying medium together with PEC is not supported yet: ca/cb would be "
                    "rewritten by two paths")
            mod = modulation.build(sim, grid, eps, sigma)
    if has_tensor:
        tensor = media.tensor_entries(sim, grid, eps, sigma)
    return Materials(eps, sigma, pec, disp, disp_mix, mod, tensor)


def _dispersive_cell_sets(disp, disp_mix):
    """The sets of dispersive cells: ``(disp_cells, disp_any, disp_cell_arr)``.

    - ``disp_cells``: flat cell sets per component, ``{comp: set}``. The σ at the source plane
      decides whether the incident-wave impedance is complex, and a plane wave looks it up on its
      own polarization component (an anisotropic medium differs component by component);
    - ``disp_any``: the flat cell set where any component carries a pole (used by the TFSF box
      shell check);
    - ``disp_cell_arr``: the cell index array of Dispersion ∪ DispersionMix (faces with poled cells
      inside the PML band need a CFS α injected, boundaries.PML_DISP_ALPHA).
    """
    disp_cell_arr = None
    if disp is not None:
        disp_cell_arr = (disp.cell if disp_mix is None
                         else np.concatenate([disp.cell, disp_mix.cell]))
    disp_cells = ({c: set(np.unique(disp.cell[disp.comp == c]).tolist()) for c in range(3)}
                  if disp is not None else None)
    disp_any = set(np.unique(disp.cell).tolist()) if disp is not None else set()
    if disp_mix is not None:
        disp_any |= set(np.unique(disp_mix.cell).tolist())
    return disp_cells, disp_any, disp_cell_arr


def _check_placements(sim, grid, mats: Materials, srcs) -> None:
    """Fail closed on every placement of the PEC mask, modulated cells or tensor cells against
    sources and absorbers whose convention has never been settled."""
    if mats.pec is not None:
        _check_pec_placement(mats.pec, grid.shape, srcs)
    if mats.mod is not None:
        modulation.check_placement(mats.mod, grid.shape, srcs)
    if mats.tensor is not None and mats.tensor.n_entry:
        _check_tensor_placement(mats.tensor, grid.shape, srcs,
                                boundaries.absorber_bands(sim, grid))


def from_simulation(sim: td.Simulation) -> Scene:
    """Build a Scene from a ``td.Simulation``.

    Symmetry boundaries **run the full domain, with no symmetry reduction**: the client's
    ``sim.grid`` is itself the full-domain grid. A localized source has to have its mirror images
    added explicitly; a plane wave does not.

    Raises:
        NotImplementedError: physics that is not supported yet. **Always fail closed**: degrading
            silently would leave later deviations impossible to attribute to anything.
    """
    media.refuse_unsupported(sim)
    # This has to be aligned explicitly, the default will not do: both directions fail silently.
    # **Set it and do not restore it** (it must not be wrapped in local_subpixel): within one
    # cloud_emulation scope, the readout code that runs after the scene is built (_perm_data in
    # nb.backend, modes.amplitudes, mode_solver_data) calls sim.epsilon and the local ModeSolver
    # directly and relies on this knob still being equal to bool(sim.subpixel). Restoring it to
    # None makes tidy3d follow "use extras if it is there" and go through the C++ engine, which
    # changes the readout convention of a subpixel=False simulation (2026-09-10 review). The
    # cloud_emulation scope restores it for everyone when it ends.
    media.set_local_subpixel(bool(sim.subpixel))
    return _build(sim)


def _build(sim: td.Simulation) -> Scene:
    grid = boundaries.build_grid(sim)
    mats = _materials(sim, grid)
    eps, sigma = mats.eps, mats.sigma
    disp_cells, disp_any, disp_cell_arr = _dispersive_cell_sets(mats.disp,
                                                                mats.disp_mix)
    pml = boundaries.build_pml(sim, grid, eps["ez"], disp_cell_arr)
    srcs = sources.build(sim, grid,
                         sources.MediaTables(eps, sigma, disp_cells, disp_any))
    # Absorber = a lossy layer with graded conductivity: σ is folded into the E update coefficients.
    #
    # **It must come after the sources.** A source computes the incident-wave impedance from the σ
    # at the source plane, and the sample point is hard-coded at the (0, 0) corner of that plane,
    # which is exactly where the absorber layers of the other two axes sit. Folding first makes the
    # plane wave treat the absorber conductivity as material loss, and the whole field of Metalens
    # came out 81 times too large (measured 2026-09-07). The absorber is a boundary, not a medium,
    # and the source is supposed not to see it.
    sigma = boundaries.fold_absorber_sigma(sim, grid, eps, sigma)
    mons = monitors.build(sim, grid)
    _check_placements(sim, grid, mats, srcs)

    return Scene(
        grid=grid,
        dt=float(sim.dt),
        num_time_steps=int(sim.num_time_steps),
        shutoff=float(sim.shutoff),
        symmetry=tuple(int(v) for v in sim.symmetry),
        bloch_k=boundaries.bloch_vector(sim),
        eps_ex=eps["ex"],
        eps_ey=eps["ey"],
        eps_ez=eps["ez"],
        sigma_ex=sigma["ex"],
        sigma_ey=sigma["ey"],
        sigma_ez=sigma["ez"],
        pec_ex=mats.pec["ex"] if mats.pec is not None else None,
        pec_ey=mats.pec["ey"] if mats.pec is not None else None,
        pec_ez=mats.pec["ez"] if mats.pec is not None else None,
        dispersion=mats.disp,
        dispersion_mix=mats.disp_mix,
        modulation=mats.mod,
        tensor=mats.tensor,
        pml=pml,
        sources=srcs.planes,
        tfsf_sources=srcs.tfsf_srcs,
        dipoles=srcs.dipoles,
        mode_sources=srcs.mode_srcs,
        flux_monitors=mons.flux,
        field_monitors=mons.field,
        field_time_monitors=mons.field_time,
        flux_time_monitors=mons.flux_time,
        permittivity_monitors=mons.permittivity,
        projection_monitors=mons.projection,
        box_flux_monitors=mons.box_flux,
        mode_monitors=mons.mode,
    )


def from_file(path: str) -> Scene:
    """Build a Scene straight from ``simulation.json``.

    Warning:
        A job containing ``CustomMedium`` or ``TriangleMesh`` **must not** use JSON, because JSON
        drops that data silently; those need ``simulation.hdf5.gz``.
    """
    return from_simulation(td.Simulation.from_file(path))
