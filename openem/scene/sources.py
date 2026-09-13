# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Mingdian Liu
"""One-sided plane wave injection, point dipoles, and mirror sources for symmetry planes."""

from __future__ import annotations

from typing import NamedTuple

import numpy as np
import tidy3d as td

from openem import cpml, tfsf1d, tfsf_oblique, waveform
from openem.grid import Grid, cyclic_axes, nearest_edge_index, normal_axis, transverse_axes
from openem.model import PlaneWaveSource, PointDipole, TFSFSource
from openem.scene import modes as modes_mod
from openem.scene._util import E_KEYS, plane_index, sign_of, upstream_h_index, yee_coords
from openem.units import UM

def incident_tables(
    source_time, dt: float, num_steps: int, dl_half: float, eps_r: float,
    direction: int = +1, sigma: float = 0.0, freq0: float | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """The two incident waveform tables for one-face TF/SF injection, timed from the source
    plane ``z_s``.

    - ``Ex_inc`` sits exactly on ``z_s`` and is sampled at whole steps ``n*dt``
    - ``Hy_inc`` sits on the neighboring Hy cell center (``z_s - dl/2`` for +z, ``z_s + dl/2``
      for -z). In both cases that point is *upstream*, i.e. ``dl/(2v)`` *earlier* than the
      source plane, so the time argument is ``(n+0.5)*dt + dl/(2v)`` either way: the direction
      only decides which cell ``dl`` is taken from, it does not change this sign.
    - The *sign* of ``Hy_inc`` does flip with direction: E × H has to point along the
      propagation direction, so if Ex is unchanged, Hy changes sign.

    The half-cell offset is folded into the time argument and ``amp_time`` is called directly,
    with *no interpolation*: that delay is not a whole or half number of steps.

    **Lossy media**: in the wave impedance ``η = η₀/n``, ``n = sqrt(ε_r + i σ/(ω ε₀))`` is
    complex; using only its real part leaves the injected H mismatched to E and the source is
    no longer one-sided. Directionality was measured to break down as ``Im ε`` grows
    (0 -> 5.5e-08, 0.043 -> 4.8e-05, 0.258 -> 2.5e-03, 1.0 -> 9.9e-02).

    How it is done: ``amp_time`` returns a *complex analytic signal*, so multiply by the complex
    ``n/η₀`` and then take the real part. The phase velocity uses ``Re(n)``; the upstream point
    has not been absorbed yet, so the amplitude gets a factor ``exp(+κ ω dl/(2c))``. When
    ``Im(n) = 0`` the whole expression degenerates *term for term* into the original real
    version (the complex constant becomes real, ``Re(z·a) = a·Re(z)``), so there is still only
    one code path.

    The complex coefficients are evaluated at ``freq0``: exact at a single frequency, a
    first-order approximation across the band. Bandwidths in the validation set are all around
    ±10%, and that ``exp`` factor is of order 1e-05 over half a cell, so the only thing that
    really acts here is the impedance.

    Args:
        dl_half: Distance from the source plane to the neighboring Hy cell center, in meters.
        eps_r: Real part of the relative permittivity at the source plane.
        direction: +1 / -1, the propagation direction of the wave.
        sigma: Conductivity at the source plane [S/m]. 0 means lossless.
        freq0: Frequency at which the complex coefficients are evaluated; not needed when
            ``sigma=0``.
    """
    eps_c = complex(eps_r)
    if sigma:
        if freq0 is None:
            raise ValueError("a plane wave in a lossy medium needs freq0 to fix the complex "
                             "wave impedance")
        eps_c += 1j * sigma / (2 * np.pi * freq0 * cpml.EPSILON_0)
    n_c = np.sqrt(eps_c)
    v = cpml.C_0 / n_c.real                      # phase velocity
    n = np.arange(num_steps + 1, dtype=np.float64)
    ex_inc = np.real(np.asarray(source_time.amp_time(n * dt)))
    t_h = (n + 0.5) * dt + dl_half / v
    # n_c/η₀ is 1/η; the upstream point is not yet absorbed, so scale by exp(+κ ω dl_half/c)
    gain = 1.0
    if sigma:
        gain = float(np.exp(n_c.imag * 2 * np.pi * freq0 * dl_half / cpml.C_0))
    coef = n_c / cpml.ETA_0 * gain
    hy_inc = direction * np.real(np.asarray(source_time.amp_time(t_h)) * coef)
    return np.ascontiguousarray(ex_inc), np.ascontiguousarray(hy_inc)


#: Old name: model.lossless (model.py) and a few tests still call it this way.
_incident_tables = incident_tables


def _linear_weights(coord: np.ndarray, x: float) -> list[tuple[int, float]]:
    """Write ``x`` as a linear combination of the two neighboring points on ``coord``, with the
    weights summing to 1.

    An exact hit returns a single entry, so a coordinate that lands exactly on a Yee point does
    not become two entries because of floating point. The flat axis of a 2D simulation (a single
    sample point) also returns a single entry.
    """
    if coord.size == 1:
        # **Flat axis of a 2D simulation**: the axis has a single Yee sample point and the field
        # is constant along it, so any target coordinate maps onto that one point. This used to
        # raise "outside the range": measured on the z axis of Autograd3InverseDesign, where the
        # only point is -0.03125 µm while the monitor/source asks for 0.0.
        return [(0, 1.0)]
    i = int(np.argmin(np.abs(coord - x)))
    if abs(float(coord[i]) - x) <= 1e-12 * max(abs(x), 1.0):
        return [(i, 1.0)]
    j = i + 1 if coord[i] < x else i - 1
    if j < 0 or j >= coord.size:
        raise ValueError(
            f"coordinate {x!r} lies outside this component's Yee coordinate range "
            f"[{coord[0]!r}, {coord[-1]!r}], cannot interpolate"
        )
    w = (x - float(coord[i])) / (float(coord[j]) - float(coord[i]))
    return [(i, 1.0 - w), (j, w)]


def _stencil_product(per_axis, coef_fn) -> tuple[list, list]:
    """Tensor product of the three axes' ``[(index, weight)]`` lists -> ``(idx, coef)``; the
    coefficient comes from ``coef_fn(i, j, k, wi, wj, wk)`` (the multiplication order lives in
    the caller's lambda and is bitwise identical to the original triple loop)."""
    idx, coef = [], []
    for i, wi in per_axis[0]:
        for j, wj in per_axis[1]:
            for k, wk in per_axis[2]:
                idx.append((i, j, k))
                coef.append(coef_fn(i, j, k, wi, wj, wk))
    return idx, coef


#: Relative tolerance for deciding that a source lies on a symmetry plane. The actual values in
#: the validation set are either exactly 0 or off by more than 0.005 µm.
ON_PLANE_RTOL = 1e-12


def reflection_sign(comp: int, dim: int, sym: int, magnetic: bool = False) -> int:
    """The sign a mirror source has to be multiplied by: ``symmetry[dim] × the eigenvalue of
    this component under reflection about that axis``.

    The eigenvalues come from the Tidy3D source, ``components/data/em_fields.py:20-25``: an E
    component along the reflection axis gets -1, otherwise +1, and **H is the other way round**
    (+1 along the axis, -1 for the rest, because H is a pseudovector). They are applied in
    ``monitor_data.py:450-452`` as ``sym_val * sym_eigenvalue``. Electric current density is a
    polar vector like E and magnetic current density is a pseudovector like H, so the two
    tables are reused as they stand.

    Physical sanity check: next to a PEC plane (``symmetry=-1``), the image of a *tangential*
    electric dipole takes -1 (which is what makes the tangential E vanish on the plane) and the
    image of a *normal* one takes +1. On a PMC plane (``symmetry=+1``), the image of a *normal*
    magnetic dipole (Bandstructure's Hz all sit on the z=0 symmetry plane) takes +1: PMC allows
    a normal H, and the source is doubled by the on-plane rule (see :func:`mirror_images`).
    """
    along = comp == dim
    eig = (1 if along else -1) if magnetic else (-1 if along else 1)
    return sym * eig


def mirror_images(
    center: tuple[float, float, float],
    comp: int,
    sym: tuple[int, int, int],
    sym_center: tuple[float, float, float],
    extent: tuple[float, float, float],
    magnetic: bool = False,
) -> list[tuple[list[float], float]]:
    """``[(position, sign), ...]``, including the original position (sign +1).

    We run the full domain while Tidy3D solves in the reduced one: a *local* source declared in
    the positive octant is equivalent, in the full domain, to itself plus its mirror images. So
    the images have to be filled in explicitly.

    **An axis on which the source sits right on the symmetry plane still produces an image.**
    That image lands back on the original position, the coefficients add, and the source is
    doubled, which is the literal meaning of symmetrizing the source distribution. Confirmed by
    measurement: all 7 of ResonanceFinder's dipoles sit on the z plane, and skipping the z image
    left our field at exactly **0.5 times** Tidy3D's (the ratio over the 7 monitors lands
    between 0.4988 and 0.5011); adding it back makes it right.

    A self-image with sign -1 would cancel the source out. That is physically correct (a source
    that lies on the plane yet is required to be odd under it has to be zero), but it is always
    a contradictory input, so fail closed instead of silently returning a zero field.

    Plane waves do not come through here: they already fill the whole domain and are consistent
    under reflection in the plane, so there is no image to add.
    """
    images: list[tuple[list[float], float]] = [([float(v) for v in center], 1.0)]
    for dim in range(3):
        if sym[dim] == 0:
            continue
        s = float(reflection_sign(comp, dim, sym[dim], magnetic))
        on_plane = abs(center[dim] - sym_center[dim]) <= ON_PLANE_RTOL * max(extent[dim], 1.0)
        if on_plane and s < 0:
            raise NotImplementedError(
                f"source lies on the {'xyz'[dim]} symmetry plane, but the mirror sign for this "
                "component is -1, which would require the source to be identically zero: "
                "contradictory input"
            )
        flipped = []
        for pos, w in images:
            q = list(pos)
            q[dim] = 2.0 * sym_center[dim] - q[dim]
            flipped.append((q, w * s))
        images = images + flipped
    return images


def dipole(sim: td.Simulation, src, grid: Grid) -> PointDipole:
    """``td.PointDipole`` -> :class:`PointDipole`, including symmetry-plane mirror images.

    The component coordinates are taken from **Tidy3D's own** ``sim.grid.yee`` rather than
    rebuilt from the convention. All images go into **one** stencil (the sign is folded into
    ``coef``), so one declared source still corresponds to a single waveform table.
    """
    pol = str(src.polarization)
    magnetic = pol[0] == "H"
    if not src.interpolate:
        raise NotImplementedError(
            "interpolate=False (snapping to the nearest Yee point) is not supported yet; "
            "the validation set uses the default True throughout"
        )
    comp = "xyz".index(pol[1])
    axes_coords = yee_coords(sim, pol[1], magnetic)
    sym = tuple(int(v) for v in sim.symmetry)
    sym_center = tuple(float(v) for v in sim.center)
    extent = tuple(abs(float(v)) for v in sim.size)
    images = mirror_images(
        tuple(float(v) for v in src.center), comp, sym, sym_center, extent,
        magnetic)

    # Tidy3D's amplitude=1 means a **current moment of 1 A·µm**, which in SI is 1e-6 A·m; the
    # current density on the grid is therefore amplitude * UM / dV.
    # By duality (J<->M, ε₀<->μ₀) a magnetic dipole means a **magnetic current moment of
    # 1 V·µm**. Like the electric case, this one is derived rather than documented: the check is
    # the vacuum radiated power against ε₀ω²|K|²/(12πc) (tests/test_dipole_h.py), confirmed end
    # to end by the FieldTimeMonitor of Bandstructure sim_0 against the reference
    # implementation's output.
    vol = (grid.yee_h_dual_volumes() if magnetic else grid.yee_dual_volumes())[comp]
    idx, coef = [], []
    for pos, sign in images:
        per_axis = [_linear_weights(axes_coords[ax], pos[ax] * UM) for ax in range(3)]
        i_, c_ = _stencil_product(
            per_axis, lambda i, j, k, wi, wj, wk: sign * wi * wj * wk * UM / float(vol[i, j, k]))
        idx += i_
        coef += c_
    return PointDipole(
        component=comp,
        indices=np.asarray(idx, dtype=np.int32),
        coef=np.asarray(coef, dtype=np.float64),
        waveform=waveform.sample(src.source_time, sim.dt, sim.num_time_steps),
        magnetic=magnetic,
    )



def _axis_weights_uniform(grid: Grid, src, ax: int, comp: int,
                          coords: np.ndarray) -> list[tuple[int, float]]:
    """``[(cell index, weight)]`` for a ``UniformCurrentSource`` along one axis (one branch per
    size case).

    - Zero thickness: reverse-interpolate onto the two neighboring planes; the δ is discretized
      as ``UM·w/μ`` using this component's Yee dual measure.
    - Filling the axis (inf): 1.0 per cell.
    - Finite non-zero size: the current density is uniform over [c-s/2, c+s/2]. Under the
      density convention the weight of each Yee point is the overlap between its own cell and
      that interval, divided by the cell length (a cell entirely inside the interval gives 1, an
      edge cell gives its covered fraction, outside gives 0). This is continuous with the
      "reverse interpolation over cell length" convention of the size=0 branch and the "1.0 per
      cell" of the inf branch: s->0 degenerates into s/dl (a point source converted to a
      density), s->inf into 1.0 per cell.
    """
    size_ax = float(src.size[ax])
    if size_ax == 0.0:
        mu = grid.axes[ax].dl if ax == comp else grid.axes[ax].dl_dual
        pairs = _linear_weights(coords, float(src.center[ax]) * UM)
        return [(i, UM * w / float(mu[i])) for i, w in pairs]
    if np.isinf(size_ax):
        return [(i, 1.0) for i in range(grid.axes[ax].n)]
    lo = (float(src.center[ax]) - 0.5 * size_ax) * UM
    hi = (float(src.center[ax]) + 0.5 * size_ax) * UM
    edges = np.asarray(grid.axes[ax].edges, dtype=np.float64)
    cell_lo, cell_hi = edges[:-1], edges[1:]     # edges.size == n+1 always holds (grid.py)
    pairs = []
    for i in range(grid.axes[ax].n):
        ov = min(hi, float(cell_hi[i])) - max(lo, float(cell_lo[i]))
        if ov > 0:
            pairs.append((i, ov / (float(cell_hi[i]) - float(cell_lo[i]))))
    if not pairs:
        raise ValueError(
            f"the {'xyz'[ax]} interval [{lo:.3e},{hi:.3e}] m of UniformCurrentSource "
            "covers no cell")
    return pairs


def uniform_current(sim: td.Simulation, src, grid: Grid) -> PointDipole:
    """``td.UniformCurrentSource`` -> current spread over Yee points (reusing
    :class:`PointDipole`).

    Only the forms that actually occur in the validation set are supported: E polarization,
    ``interpolate=True``, amplitude convention ``density``, and each axis either zero thickness
    (reverse-interpolated onto the two neighboring planes) or filling the axis (inf).

    **Amplitude convention** (the same derivation chain as the dipole above): under the
    ``density`` convention ``amplitude=1`` is a current density of 1 A/µm², and a zero-thickness
    axis carries a δ. Discretizing that δ onto the Yee sample line of the axis (weight w over
    the dual measure μ at that point), the SI current density is::

        J_cell = amp · 1e12 · ∏_{zero-thickness axes d} (UM · w_d / μ_d)

    Checked term by term: with all three axes at zero this is a point dipole and the formula
    gives ``amp·1e-6·w³/dV``, exactly :func:`dipole`'s ``coef = w·UM/dV``; with a single zero
    axis (a current sheet) it gives ``amp·1e6·w/μ``, i.e. the surface current density
    ``K = amp A/µm`` laid into the dual cell. The check is the analytic solution for a current
    sheet in vacuum, ``|E| = η|K|/2`` (tests/test_tensor.py, no free parameters).

    μ is this component's Yee dual measure: axis == component axis -> primal ``dl``, otherwise
    dual ``dl_dual`` (the same convention as ``grid.yee_dual_volumes``).
    """
    pol = str(src.polarization)
    if pol[0] != "E":
        raise NotImplementedError(
            f"magnetic-current polarization ({pol}) of UniformCurrentSource is not supported "
            "yet: it would have to be injected into the H update")
    if not src.interpolate:
        raise NotImplementedError("UniformCurrentSource interpolate=False is not supported yet")
    if str(getattr(src, "current_amplitude_definition", "density")) != "density":
        raise NotImplementedError(
            "UniformCurrentSource only supports current_amplitude_definition='density'")
    if any(int(v) != 0 for v in sim.symmetry):
        raise NotImplementedError(
            "UniformCurrentSource with symmetry boundaries is not supported yet: the mirror "
            "rules have never been checked")

    comp = "xyz".index(pol[1])
    axes_coords = yee_coords(sim, pol[1])
    per_axis = [_axis_weights_uniform(grid, src, ax, comp, axes_coords[ax])
                for ax in range(3)]

    idx, coef = _stencil_product(per_axis, lambda i, j, k, wi, wj, wk: 1.0e12 * wi * wj * wk)
    return PointDipole(
        component=comp,
        indices=np.asarray(idx, dtype=np.int32),
        coef=np.asarray(coef, dtype=np.float64),
        waveform=waveform.sample(src.source_time, sim.dt, sim.num_time_steps),
    )


def _split_complex_waveform(wf: waveform.Waveform) -> tuple[waveform.Waveform, waveform.Waveform]:
    """``(Re(g) table, Im(g) table)``: complex-amplitude injection split into two real tables
    (see :func:`custom_current`)."""
    wf_im = waveform.Waveform(
        amp_int=np.imag(wf.amp_int_complex).astype(np.float64),
        amp_half=np.imag(wf.amp_half_complex).astype(np.float64),
        amp_int_complex=wf.amp_int_complex, dt=wf.dt,
        amp_half_complex=wf.amp_half_complex)
    return wf, wf_im


def _targets_for_axis(grid: Grid, src, ax: int, comp: int, coords_yee: np.ndarray,
                      cds: np.ndarray) -> tuple[list[tuple[int, float]], np.ndarray]:
    """Target Yee coordinates plus weights along one axis:
    ``([(cell index, weight)], target physical coordinates in m)``.

    - Zero-thickness injection axis (src.size==0, more than 1 grid cell): the target is the
      source center plane and the δ is discretized as UM·w/μ against the dual measure; when the
      dataset has more than one coordinate along that axis it is first interpolated onto the
      plane (a probe with three z points was measured to be off by a factor of 403 = (1/dz)²,
      which is exactly this missing δ).
    - 2D flat axis: the axis has a single cell, that cell is the target, weight 1.
    - Finite transverse axis: the targets are every Yee point inside the source box, weight 1;
      outside the dataset's range the value is filled with 0.
    """
    if float(src.size[ax]) == 0.0 and grid.axes[ax].n > 1:
        mu = grid.axes[ax].dl if ax == comp else grid.axes[ax].dl_dual
        plane = float(src.center[ax]) * UM
        pairs = _linear_weights(coords_yee, plane)
        return [(i, UM * w / float(mu[i])) for i, w in pairs], np.full(len(pairs), plane)
    if grid.axes[ax].n == 1:
        return [(0, 1.0)], np.array([float(cds[0]) if cds.size == 1 else float(coords_yee[0])])
    half = 0.5 * float(src.size[ax]) * UM + 1e-9 * UM
    c0 = float(src.center[ax]) * UM
    sel = np.flatnonzero(np.abs(coords_yee - c0) <= half)
    return [(int(i), 1.0) for i in sel], coords_yee[sel]


def _interp_dataset(vals: np.ndarray, ds_coords: list, tgt_pos: list) -> np.ndarray:
    """Interpolate the dataset onto the target points: an axis with a single coordinate is not
    interpolated (that one value is taken), the others are linear with 0 outside the range."""
    from scipy.interpolate import RegularGridInterpolator as _RGI

    interp_axes = [ax for ax in range(3) if ds_coords[ax].size > 1]
    squeeze = tuple(ax for ax in range(3) if ds_coords[ax].size == 1)
    v_src = vals
    for ax in sorted(squeeze, reverse=True):
        v_src = np.take(v_src, 0, axis=ax)
    P1, P2, P3 = np.meshgrid(*tgt_pos, indexing="ij")
    # Point-source dataset (a single coordinate on all three axes; the adjoint
    # CustomCurrentSource of Autograd15Antenna's point-field objective): no axis is left to
    # interpolate over and np.stack([]) blows up with "need at least one array to stack", so
    # broadcast that one value directly (2026-09-05)
    if interp_axes:
        pts = np.stack([np.ravel(P) for P, ax in zip((P1, P2, P3), range(3)) if ax in interp_axes], axis=-1)
        grids = [ds_coords[ax] for ax in interp_axes]
        fr = _RGI(grids, v_src.real, method="linear", bounds_error=False, fill_value=0.0)
        fi = _RGI(grids, v_src.imag, method="linear", bounds_error=False, fill_value=0.0)
        vg = (fr(pts) + 1j * fi(pts)).reshape(P1.shape)
    else:
        vg = np.broadcast_to(np.asarray(v_src, dtype=np.complex128).reshape(()), P1.shape)
    return np.asarray(vg, dtype=np.complex128)


def _component_dipoles(comp: int, magnetic: bool, tgt_idx: list, vals: np.ndarray,
                       wf, wf_im) -> list[PointDipole]:
    """Target points of one component -> at most two :class:`PointDipole` objects (the Re(g)
    table with coefficient +Re(a), the Im(g) table with -Im(a)). The coefficients
    ``g = 1.0e12*wi*wj*wk; cre = g*v.real; cim = -g*v.imag`` are kept as they are."""
    idx, cre, cim = [], [], []
    for a0, (i, wi) in enumerate(tgt_idx[0]):
        for a1, (j, wj) in enumerate(tgt_idx[1]):
            for a2, (kk, wk) in enumerate(tgt_idx[2]):
                v = complex(vals[a0, a1, a2])
                if v == 0:
                    continue
                idx.append((i, j, kk))
                g = 1.0e12 * wi * wj * wk
                cre.append(g * v.real)
                cim.append(-g * v.imag)
    out: list[PointDipole] = []
    if not idx:
        return out
    ai = np.asarray(idx, dtype=np.int32)
    for coefs, w in ((cre, wf), (cim, wf_im)):
        c = np.asarray(coefs, dtype=np.float64)
        if not np.any(c):
            continue
        out.append(PointDipole(component=comp, indices=ai, coef=c,
                               waveform=w, magnetic=magnetic))
    return out


def _dataset_values(key: str, arr) -> np.ndarray | None:
    """Complex amplitude ``(x, y, z)`` of one component of the dataset; all zeros returns None
    (no dipole is generated).

    Raises:
        NotImplementedError: There is a frequency axis with more than one frequency; a single
            source can only carry a current pattern at one frequency.
    """
    vals = np.asarray(arr.values, dtype=np.complex128)
    if vals.ndim == 4:                     # (x, y, z, f)
        if vals.shape[3] != 1:
            raise NotImplementedError(
                f"component {key} of CustomCurrentSource carries {vals.shape[3]} frequencies: "
                "a single source can only have a current pattern at one frequency")
        vals = vals[..., 0]
    return vals if np.any(vals) else None


def custom_current(sim: td.Simulation, src, grid: Grid) -> list[PointDipole]:
    """``td.CustomCurrentSource`` -> a set of :class:`PointDipole`.

    This is exactly the type tidy3d uses for its **adjoint source**: a point-by-point complex
    amplitude current over the monitor plane. In the dataset ``Ex/Ey/Ez`` is the electric
    current J and ``Hx/Hy/Hz`` the magnetic current M, with coordinates **relative to
    ``center``** (simulation coordinate = dataset coordinate + center).

    **Complex amplitudes need no kernel change**: the real signal being injected is
    ``Re(a·g(t)) = Re(a)·Re(g) - Im(a)·Im(g)``, where g is the analytic signal from
    ``amp_time``. So each component is split into **two** PointDipoles, one carrying the
    ``Re(g)`` table with coefficient ``+Re(a)``, one carrying the ``Im(g)`` table with
    ``-Im(a)``. Both are real, so the existing batch injection takes them as they are.

    The **amplitude convention** follows :func:`uniform_current`'s derivation chain: dataset
    values are current *densities* (A/µm² for electric current, V/µm² for magnetic current) and
    are multiplied by ``1e12`` to get SI. When the dataset has only one coordinate along some
    axis (the normal of a planar source) that axis carries a δ, discretized as ``UM·w/μ``
    against this component's Yee dual measure, term for term the same as the zero-thickness axis
    of ``uniform_current``. Non-degenerate axes snap to the nearest Yee sample point (weight 1).
    **The flat axis of a 2D simulation does not count as zero thickness**: it has a single
    coordinate only because the whole simulation is one cell thick, and treating it as a δ would
    multiply in a spurious ``UM/dl``.

    Raises:
        NotImplementedError: ``interpolate=False``, ``confine_to_bounds=True``, no
            ``current_dataset``, a dataset with several frequencies, or an all-zero dataset.
    """
    if not getattr(src, "interpolate", True):
        raise NotImplementedError("CustomCurrentSource interpolate=False is not supported yet")
    if getattr(src, "confine_to_bounds", False):
        raise NotImplementedError(
            "CustomCurrentSource confine_to_bounds=True is not supported yet")
    ds = getattr(src, "current_dataset", None)
    if ds is None:
        raise NotImplementedError("CustomCurrentSource has no current_dataset")

    wf, wf_im = _split_complex_waveform(waveform.sample(src.source_time, sim.dt, sim.num_time_steps))

    out: list[PointDipole] = []
    for key in ("Ex", "Ey", "Ez", "Hx", "Hy", "Hz"):
        arr = getattr(ds, key, None)
        if arr is None:
            continue
        vals = _dataset_values(key, arr)
        if vals is None:
            continue
        magnetic = key[0] == "H"
        comp = "xyz".index(key[1])
        axes_yee = yee_coords(sim, key[1], magnetic)
        ds_coords = [np.asarray(arr.coords[a].values, dtype=np.float64) * UM
                     + float(src.center[i]) * UM for i, a in enumerate("xyz")]

        # per_axis[ax][n] = the [(cell index, weight), ...] that data point n of this axis
        # spreads over. Dataset -> Yee grid is done by **interpolation**, not by "snap each data
        # point to the nearest cell". Snapping skips cells when the dataset is coarser than the
        # grid (CustomFieldSource tutorial: field amplitude down to 1/11, flux to 1/127) and
        # piles several points onto one cell when it is finer. Tidy3D's semantics are exactly a
        # linear interpolation of the dataset onto its own grid (interpolate=True).
        # Each axis yields "target Yee coordinates + weights":
        #   Zero-thickness injection axis (src.size==0, more than 1 grid cell): the target is
        #     the source center plane and the δ is discretized as UM·w/μ against the dual
        #     measure; when the dataset has more than one coordinate along that axis it is first
        #     interpolated onto the plane (a probe with three z points was measured to be off by
        #     a factor of 403 = (1/dz)², which is exactly this missing δ).
        #   2D flat axis: the axis has a single cell, that cell is the target, weight 1.
        #   Finite transverse axis: the targets are every Yee point inside the source box,
        #     weight 1; outside the dataset's range the value is filled with 0.
        tgt_idx: list[list[tuple[int, float]]] = []      # per axis: targets [(cell idx, weight)]
        tgt_pos: list[np.ndarray] = []                    # per axis: target coords (m), for interp
        for ax in range(3):
            ti_, tp_ = _targets_for_axis(grid, src, ax, comp, axes_yee[ax], ds_coords[ax])
            tgt_idx.append(ti_)
            tgt_pos.append(tp_)
        if any(len(t) == 0 for t in tgt_idx):
            continue
        # tgt_idx[ax][n] = (cell index, weight), vals[a0,a1,a2] = value at that target point
        vals = _interp_dataset(vals, ds_coords, tgt_pos)
        if tuple(len(t) for t in tgt_idx) != vals.shape:
            raise NotImplementedError(
                f"after interpolation, {key} of CustomCurrentSource has shape {vals.shape}, "
                f"which does not match the target point counts "
                f"{tuple(len(t) for t in tgt_idx)} of the three axes")
        out.extend(_component_dipoles(comp, magnetic, tgt_idx, vals, wf, wf_im))
    if not out:
        raise NotImplementedError("the dataset of CustomCurrentSource is all zeros")
    return out


def custom_field(sim: td.Simulation, src, grid: Grid) -> list[PointDipole]:
    """``td.CustomFieldSource`` -> a set of :class:`PointDipole` (equivalent currents).

    ``CustomFieldSource`` supplies the **tangential complex E/H fields on a plane**; the
    injection axis is the one whose size is 0 (``src.injection_axis``). The surface equivalence
    principle turns them into equivalent surface electric / magnetic currents:

        J_s = n̂ × H,   M_s = -n̂ × E

    n̂ is the normal of the injection plane (``+injection_axis``). Giving both E and H makes a
    Huygens source, which is directional by itself; if only one of them is given, OpenEM goes
    along with that too (it is then not one-sided, which matches tidy3d's semantics).

    The resulting J (electric current) and M (magnetic current) are packed into a
    ``current_dataset`` (Ex/Ey/Ez = J, Hx/Hy/Hz = M) and handed straight to
    :func:`custom_current`: the amplitude convention, the split of a complex amplitude into two
    tables, the δ discretization of the zero-thickness axis and magnetic current injection are
    all reused from there instead of being written a second time.

    Raises:
        NotImplementedError: No ``field_dataset``, a leftover tangential component on the
            injection axis, or an all-zero dataset after the conversion (reported by
            :func:`custom_current`).
    """
    ds = getattr(src, "field_dataset", None)
    if ds is None:
        raise NotImplementedError("CustomFieldSource has no field_dataset")
    n = int(src.injection_axis)              # normal axis of the injection plane
    # n̂ × v componentwise: cross[i] = Σ ε(n,j,i)·v[j]. For the unit normal n̂ = ê_n,
    # (n̂ × v)_i = ε_{n j i} v_j, i.e. components (n+1)%3 and (n+2)%3 of v are moved around in
    # cyclic order.
    t1, t2 = cyclic_axes(n)                  # the two tangential axes (cyclic order)
    Ec = [getattr(ds, f"E{a}", None) for a in "xyz"]
    Hc = [getattr(ds, f"H{a}", None) for a in "xyz"]

    # n̂ × F: the result lands on the tangential axes. (n̂×F)_{t1} = +F_{t2}, (n̂×F)_{t2} = -F_{t1}
    #   (cyclic order n->t1->t2->n, ê_n × ê_t1 = ê_t2, ê_n × ê_t2 = -ê_t1, hence
    #    ê_n × (F_t1 ê_t1 + F_t2 ê_t2) = F_t1 ê_t2 - F_t2 ê_t1)
    def _cross_n(F):
        out = [None, None, None]
        out[t1] = _neg(F[t2])                 # -F_{t2}
        out[t2] = F[t1]                        # +F_{t1}
        return out

    J = _cross_n(Hc)                           # J_s = n̂ × H
    M = [_neg(x) for x in _cross_n(Ec)]        # M_s = -n̂ × E

    comps = {}
    for ax, key in enumerate("xyz"):
        if J[ax] is not None:
            comps[f"E{key}"] = J[ax]           # electric current goes into E* of the current src
        if M[ax] is not None:
            comps[f"H{key}"] = M[ax]           # magnetic current goes into H*
    if not comps:
        raise NotImplementedError("the field_dataset of CustomFieldSource has no tangential "
                                  "component")
    cur_ds = td.FieldDataset(**comps)
    equiv = td.CustomCurrentSource(
        center=src.center, size=src.size, source_time=src.source_time,
        current_dataset=cur_ds,
        interpolate=getattr(src, "interpolate", True))
    return custom_current(sim, equiv, grid)


class MediaTables(NamedTuple):
    """The medium-side inputs needed to build a source: four things that get passed all the way
    down, and spelling them out one by one in the signatures would drown the real parameters.

    ``eps`` / ``sigma`` are ``{"ex": ..., "ey": ..., "ez": ...}`` (the values of ``sigma`` may
    all be None); ``disp_cells`` is the per-component set of flat indices of dispersive cells,
    ``{comp: set}`` (a plane wave source checks against its own polarization component), and
    ``disp_any`` is the union over all components (used when checking the TFSF box shell).
    """

    eps: dict[str, np.ndarray]
    sigma: dict[str, np.ndarray | None]
    disp_cells: dict[int, set] | None = None
    disp_any: set[int] | None = None


def _neg(arr):
    """Negate a DataArray; ``None`` is returned unchanged."""
    if arr is None:
        return None
    return -arr


def plane_wave(sim: td.Simulation, src, grid: Grid, med: MediaTables) -> PlaneWaveSource:
    """``td.PlaneWave`` -> :class:`PlaneWaveSource` (one-face TF/SF injection).

    The incident wave impedance takes the ε/σ of the **polarization component**: in an
    anisotropic medium it is the ε along the polarization direction that the wave sees. In an
    isotropic scene all three components agree at the source plane, so nothing changes.

    Polarization: with ``pol_angle=0``, E lies along the first transverse axis; with ``±π/2``,
    along the second. ``-π/2`` and ``+π/2`` differ only by an overall sign of the waveform and
    the injection tables do not tell them apart: flux, and everything normalized between the
    structure run and the reference run, is independent of an overall sign.
    """
    if abs(src.angle_theta) > 0:
        raise NotImplementedError(
            f"only normal incidence is supported, angle_theta={src.angle_theta}")
    if abs(float(getattr(src, "angle_phi", 0.0))) > 0:
        raise NotImplementedError(f"only angle_phi=0 is supported, got {src.angle_phi}")

    ax = normal_axis(src.size)
    trans = transverse_axes(ax)
    pol = float(src.pol_angle)
    if abs(pol) < 1e-12:
        pol_axis = trans[0]
    elif abs(abs(pol) - np.pi / 2) < 1e-12:
        pol_axis = trans[1]
    else:
        raise NotImplementedError(f"only pol_angle ∈ {{0, ±π/2}} is supported, got {src.pol_angle}")
    eps_pol = med.eps[E_KEYS[pol_axis]]
    sigma_pol = med.sigma[E_KEYS[pol_axis]]
    direction = sign_of(src.direction)
    ks = plane_index(grid, ax, src.center[ax])
    # Neighboring Hy cell center: ks-1 for +z, ks for -z
    dl_h = float(grid.axes[ax].dl[upstream_h_index(ks, direction)])
    if ax != 2:
        # build routes every plane wave with normal_axis != 2 to plane_wave_beam; the 1D
        # auxiliary grid kernel only supports propagation along z (the src.axis != 2 guard in
        # sources_setup._plane_source).
        raise NotImplementedError(
            f"the 1D auxiliary grid plane wave only supports a z normal, got axis {ax}")
    at = (0, 0, ks)
    if med.disp_cells and _flat(at, grid.shape) in med.disp_cells.get(pol_axis, ()):
        raise NotImplementedError(
            f"the source plane {at} of the plane wave lies in a **dispersive** medium: the "
            "incident wave impedance would need ε(ω), and the only complex ε supported so far "
            "is ε_r + iσ/(ωε₀)."
        )
    # Keep the raw parameters around: :func:`openem.model.lossless` recomputes from them
    inc_kw = dict(
        source_time=src.source_time, dt=sim.dt, num_steps=sim.num_time_steps,
        dl_half=dl_h / 2.0,
        eps_r=float(eps_pol[at]),
        direction=direction,
        sigma=0.0 if sigma_pol is None else float(sigma_pol[at]),
        freq0=float(src.source_time.freq0),
    )
    ex_inc, hy_inc = incident_tables(**inc_kw)
    return PlaneWaveSource(
        axis=ax,
        plane_index=ks,
        direction=direction,
        pol_axis=pol_axis,
        waveform=waveform.sample(src.source_time, sim.dt, sim.num_time_steps),
        angle_theta=float(src.angle_theta),
        num_freqs=int(src.num_freqs or 1),
        ex_inc=ex_inc,
        hy_inc=hy_inc,
        center=float(src.center[ax]) * UM,
        incident=inc_kw,
    )


def _flat(idx: tuple[int, int, int], shape: tuple[int, int, int]) -> int:
    """3D index -> flat index, the same convention as ``cell`` in :class:`Dispersion`."""
    return (idx[0] * shape[1] + idx[1]) * shape[2] + idx[2]


#: Extra 3D cells the TFSF box column carries beyond the box (the buffer between the injection
#: plane and the box face).
TFSF_GAP = 4

#: Relative tolerance for calling two ε values numerically equal. Subpixel averaging gives
#: column-by-column identical values for a transversely uniform layered structure, so the only
#: difference should come from the floating-point representation.
TFSF_EPS_RTOL = 1e-9

#: Relative tolerance for calling the grid uniform inside the box in oblique-incidence TFSF
#: (``np.ptp(dl)/mean(dl)``).
TFSF_DL_RTOL = 1e-12


def _column_segment(col: np.ndarray, c0: int, c1: int) -> np.ndarray:
    """The segment of the box column on the 1D edges, ``col[c0 : c1+2]``; if the tail runs out
    (``c1+1`` is already the last cell) one cell is padded with the end value, because eps1 and
    sigma1 are edge quantities and there is one more edge than cell."""
    if c1 + 1 < col.size:
        return col[c0:c1 + 2]
    return np.concatenate([col[c0:], [col[-1]]])


def _tfsf_ring(ulo: int, uhi: int, vlo: int, vhi: int) -> tuple[np.ndarray, np.ndarray]:
    """The **perimeter** cells ``(uu, vv)`` of the box footprint grown by one cell (in the
    transverse plane), used to check transverse uniformity.

    Only the perimeter, not the whole footprint: the scatterer inside the box (a sphere) sits in
    the middle of the footprint and taking the whole patch would sweep it in too. The perimeter
    is closest to the box faces, which is where a structure cutting into a box face shows up
    first.
    """
    uu, vv = [], []
    for i in range(ulo - 1, uhi + 1):
        for j in (vlo - 1, vlo, vhi - 1, vhi):
            uu.append(i)
            vv.append(j)
    for j in range(vlo + 1, vhi - 1):
        for i in (ulo - 1, ulo, uhi - 1, uhi):
            uu.append(i)
            vv.append(j)
    return np.asarray(uu), np.asarray(vv)


class _TfsfBox(NamedTuple):
    """Integer geometry of the TFSF box, shared by the normal- and oblique-incidence paths.

    ``a`` is the injection axis; ``lo`` / ``hi`` are the primal edge indices of the six faces
    (in xyz order, the two lists returned by :func:`tfsf_box_edges`). The derived names below
    are the same indices rewritten in (u, v, a) axis order: the transverse axes are
    ``u = (a+1)%3`` and ``v = (a+2)%3``, taken in cyclic right-handed order.
    ``c0`` / ``c1`` are the range of 3D cells covered by the 1D box column (:data:`TFSF_GAP`
    extra cells beyond each box face); only the normal-incidence path reads them, because on an
    open axis ``lo``/``hi`` are the two ends of the whole axis and there is no box column at all.
    """

    a: int
    lo: list[int]
    hi: list[int]

    @property
    def u(self) -> int:
        return cyclic_axes(self.a)[0]

    @property
    def v(self) -> int:
        return cyclic_axes(self.a)[1]

    @property
    def alo(self) -> int:
        return self.lo[self.a]

    @property
    def ahi(self) -> int:
        return self.hi[self.a]

    @property
    def ulo(self) -> int:
        return self.lo[self.u]

    @property
    def uhi(self) -> int:
        return self.hi[self.u]

    @property
    def vlo(self) -> int:
        return self.lo[self.v]

    @property
    def vhi(self) -> int:
        return self.hi[self.v]

    @property
    def c0(self) -> int:
        return self.alo - TFSF_GAP

    @property
    def c1(self) -> int:
        return self.ahi + TFSF_GAP - 1

    def uva(self, arr):
        """Reorder an array to (u, v, a) axis order. ``moveaxis`` returns a view, so the values
        are unchanged."""
        return np.moveaxis(arr, (self.u, self.v, self.a), (0, 1, 2))

    def ring(self) -> tuple[np.ndarray, np.ndarray]:
        """The perimeter cells of the box footprint grown by one cell (see :func:`_tfsf_ring`)."""
        return _tfsf_ring(self.ulo, self.uhi, self.vlo, self.vhi)


class _Aux1D(NamedTuple):
    """The 1D auxiliary line together with its two injection tables: the group of fields inside
    :class:`TFSFSource` that describes the incident field."""

    dl1: np.ndarray
    eps1: np.ndarray
    sigma1: np.ndarray | None
    inj: int
    col0: int
    ex_inc: np.ndarray
    hy_inc: np.ndarray


def tfsf_box_edges(src, grid):
    """The primal edge indices of the six TFSF box faces, plus the list of open axes.

    tidy3d allows a TFSF size of ``td.inf`` (tfsf2 in TFSF.ipynb is ``size=(inf, inf, 4)``, with
    periodic/Bloch boundaries): the box then fills the whole computational domain along that
    axis, those two faces do not exist, and there is no TF/SF correction there. Such an axis is
    recorded as open and its indices are the two ends of the whole axis.

    Raises:
        NotImplementedError: The injection axis itself is infinite; there is then no box at all,
            so TF/SF makes no sense.
        ValueError: Only one end of an axis is infinite (tidy3d never generates this, so if it
            shows up the data is broken).
    """
    lo, hi = src.bounds
    a = int(src.injection_axis)
    lo_e, hi_e, open_axes = [], [], []
    for q in range(3):
        f_lo, f_hi = np.isfinite(float(lo[q])), np.isfinite(float(hi[q]))
        if f_lo and f_hi:
            lo_e.append(nearest_edge_index(grid.axes[q].edges, float(lo[q]) * UM))
            hi_e.append(nearest_edge_index(grid.axes[q].edges, float(hi[q]) * UM))
            continue
        if f_lo != f_hi:
            raise ValueError(
                f"the TFSF box is infinite at only one end of the {'xyz'[q]} axis "
                f"({lo[q]}, {hi[q]})")
        if q == a:
            raise NotImplementedError(
                "the TFSF source has infinite size along the **injection axis**: there is then "
                "no total-field/scattered-field interface and TF/SF makes no sense. Either give "
                "this notebook a finite box or run the whole notebook on Tidy3D's solver.")
        open_axes.append(q)
        lo_e.append(0)
        hi_e.append(grid.axes[q].edges.size - 1)
    return lo_e, hi_e, tuple(open_axes)


#: Old name (tests/test_tfsf.py refers to it this way).
_tfsf_box_edges = tfsf_box_edges


def _oblique_uniform_axes(ax, bx: _TfsfBox, open_axes) -> None:
    """Each of the three axes must have a uniform grid inside the box, otherwise the projected
    position is no longer a linear function of the coordinate on the face: fail closed. On an
    open axis the box fills the whole axis, so the whole axis is checked."""
    for q in range(3):
        seg = ax[q].dl if q in open_axes else ax[q].dl[max(bx.lo[q] - 2, 0):bx.hi[q] + 2]
        if float(np.ptp(seg)) > TFSF_DL_RTOL * float(seg.mean()):
            raise NotImplementedError(
                f"oblique-incidence TFSF requires a uniform grid on the {'xyz'[q]} axis inside "
                f"the box, got dl ∈ [{seg.min():.4e}, {seg.max():.4e}]: the projected position "
                "is then no longer a linear function of the coordinate on the face")


def _oblique_shell(bx: _TfsfBox, shape, open_axes) -> np.ndarray:
    """Mask of the **box shell** whose background has to be checked: the box grown by two cells,
    with the interior carved out.

    The TF/SF correction only acts on the shell, and a scatterer inside the box has nothing to
    do with the incident field (FullyAnisotropic's sphere sits right in the middle of the box).
    This is the same criterion as the perimeter ring of the normal-incidence path, except that
    oblique incidence has to check all six faces. An open axis has no faces, so the whole axis
    counts as inside the box: otherwise the two cells at each end of the computational domain
    would be taken for shell and the background there would be required to be uniform for
    nothing (under periodic boundaries that is exactly where a structure may sit flush).
    """
    shell = np.zeros(shape, dtype=bool)
    shell[tuple(slice(max(bx.lo[q] - 2, 0), bx.hi[q] + 3) for q in range(3))] = True
    inner = tuple(slice(0, shape[q]) if q in open_axes
                  else slice(bx.lo[q] + 2, bx.hi[q] - 1) for q in range(3))
    if all(s.stop > s.start for s in inner):
        shell[inner] = False
    return shell


def _oblique_background(med: MediaTables, shell: np.ndarray) -> float:
    """The background ``eps_bg`` on the box shell. Non-uniform, lossy or dispersive backgrounds
    all fail to define a unique incident field, so fail closed."""
    ref = None
    for c, e in med.eps.items():
        blk = e[shell]
        if ref is None:
            ref = float(blk.flat[0])
        dev = float(np.abs(blk - ref).max())
        if dev > TFSF_EPS_RTOL * max(1.0, abs(ref)):
            raise NotImplementedError(
                f"eps_{c} is not uniform on the box shell of the oblique-incidence TFSF "
                f"(largest deviation {dev:.3e}): the 1D line runs along k̂, and a layered or "
                "anisotropic background does not define a unique incident field")
    for s in med.sigma.values():
        if s is not None and float(np.abs(s[shell]).max()) > 0.0:
            raise NotImplementedError("oblique-incidence TFSF does not support a lossy box "
                                      "shell yet")
    if med.disp_any:
        raise NotImplementedError("oblique-incidence TFSF does not support a dispersive "
                                  "background yet")
    return float(ref)


def _tfsf_oblique(sim, src, grid: Grid, med: MediaTables) -> TFSFSource:
    """Oblique-incidence TFSF (``angle_theta ≠ 0``). The 1D auxiliary line runs along k̂; see
    :mod:`openem.tfsf_oblique`.

    The restrictions all fail closed rather than guessing: inside the box **each of the three
    axes must have a uniform grid** (otherwise the projected position is not linear in the
    coordinate on the face), and on the **box shell** the background must be isotropic, uniform,
    lossless and non-dispersive (otherwise one 1D line cannot define a unique incident field).
    *Inside* the box any structure is allowed, since the TF/SF correction only acts on the
    shell.
    """
    a = int(src.injection_axis)
    k_hat = np.asarray(src._dir_vector, dtype=np.float64)
    e_hat = np.asarray(src._pol_vector, dtype=np.float64)
    ax = [grid.axes[q] for q in range(3)]
    lo_e, hi_e, open_axes = tfsf_box_edges(src, grid)
    bx = _TfsfBox(a, lo_e, hi_e)
    shp = grid.shape
    # On an open axis the box already touches both ends of the domain; with no face there is
    # no margin to require
    if not all(2 <= lo_e[q] < hi_e[q] <= shp[q] - 3
               for q in range(3) if q not in open_axes):
        raise ValueError(
            f"the oblique-incidence TFSF box {lo_e}..{hi_e} is too close to the grid boundary "
            f"({shp}): each of the six face corrections needs two cells of margin")
    _oblique_uniform_axes(ax, bx, open_axes)
    eps_bg = _oblique_background(med, _oblique_shell(bx, shp, open_axes))
    dt = float(sim.dt)
    n_steps = int(sim.num_time_steps)
    freq0 = float(src.source_time.freq0)
    dl_axes = [float(ax[q].dl[lo_e[q]]) for q in range(3)]
    dl1v = tfsf1d.matched_dl(k_hat, dl_axes, dt, freq0, eps_r=eps_bg)
    pmin, pmax = tfsf_oblique.projection_range(grid, k_hat, lo_e, hi_e)
    proj0 = pmin - 2.0 * dl1v
    n_col_e = int(np.ceil((pmax - proj0) / dl1v)) + 3
    n_run = tfsf1d.runway_cells(n_steps, dt, dl1v)
    dl1 = np.full(2 * n_run + n_col_e, dl1v, dtype=np.float64)
    eps1 = np.full(dl1.size, eps_bg, dtype=np.float64)
    col0 = n_run                       # the projection proj0 lands on 1D edge col0
    inj = col0 - TFSF_GAP // 2         # injection plane, a few cells upstream of the box
    ex_inc, hy_inc = incident_tables(
        source_time=src.source_time, dt=dt, num_steps=n_steps,
        dl_half=dl1v / 2.0, eps_r=eps_bg,
        direction=+1,                  # the line runs along +k̂, direction is folded into k̂
        sigma=0.0, freq0=freq0)
    # Tidy3D's amplitude convention, in the words of the td.TFSF docs: it "injects 1 W of power
    # per µm² of source area **along the injection_axis**", and under angled incidence the power
    # injected along the injection axis stays the same. The power density normal to the
    # injection axis is I·cosθ, so to make it equal the normal-incidence value the plane wave
    # intensity has to be larger by 1/cosθ and the amplitude by 1/√cosθ. At normal incidence
    # cosθ=1 and this is the identity.
    #
    # **The factor must not be multiplied into ex_inc**: that table is at the same time the
    # normalization divisor (normalize.source_table -> plane_flux_divisor), so multiplying it in
    # would scale the divisor too and cancel exactly (measured: flux_theta/flux_diag was still
    # 1.0035, unmoved). Tidy3D normalizes by the **source waveform spectrum**, which carries no
    # such geometric factor, so the table keeps the reference amplitude, the factor is recorded
    # on its own, and it is applied only where sources_setup._tfsf_box builds the 1D tables.
    amp_scale = 1.0 / np.sqrt(abs(float(k_hat[a])))
    return TFSFSource(
        **_tfsf_common(src, bx,
                       _Aux1D(dl1, eps1, None, inj, col0, ex_inc, hy_inc),
                       dt, n_steps),
        direction=+1,
        pol_u=float(e_hat[(a + 1) % 3]), pol_v=float(e_hat[(a + 2) % 3]),
        k_hat=tuple(float(x) for x in k_hat),
        e_hat=tuple(float(x) for x in e_hat),
        dl1_step=float(dl1v), proj0=float(proj0), n_col_e=int(n_col_e),
        amp_scale=float(amp_scale), open_axes=open_axes)


def _tfsf_common(src, bx: _TfsfBox, aux: _Aux1D, dt: float, n_steps: int) -> dict:
    """The :class:`TFSFSource` fields shared by the normal- and oblique-incidence paths."""
    return {
        "box_lo": tuple(bx.lo), "box_hi": tuple(bx.hi),
        "dl1": aux.dl1, "eps1": aux.eps1, "sigma1": aux.sigma1,
        "inj": int(aux.inj), "col0": int(aux.col0),
        "ex_inc": aux.ex_inc, "hy_inc": aux.hy_inc,
        "waveform": waveform.sample(src.source_time, dt, n_steps),
        "num_freqs": int(src.num_freqs or 1), "axis": bx.a,
    }


def _tfsf_check_background(med: MediaTables, bx: _TfsfBox, grid) -> None:
    """Background checks for normal-incidence TFSF: ε and σ transversely uniform over the
    perimeter ring of the box and over the footprints of the two a faces, and neither the box
    shell nor the box column lying in a dispersive medium. Otherwise the incident field cannot
    be defined, so fail closed."""
    a, u, v = bx.a, bx.u, bx.v
    ulo, uhi, vlo, vhi = bx.ulo, bx.uhi, bx.vlo, bx.vhi
    c0, c1 = bx.c0, bx.c1
    alo, ahi = bx.alo, bx.ahi
    uu, vv = bx.ring()
    em_all = {c: bx.uva(e) for c, e in med.eps.items()}   # moveaxis is a view, values unchanged
    for c, em in em_all.items():
        ring = em[uu, vv, :]                           # (n_ring, na)
        p = em[ulo, vlo, :].astype(np.float64)
        dev = np.abs(ring[:, c0:c1 + 1] - p[c0:c1 + 1]).max()
        if dev > TFSF_EPS_RTOL * max(1.0, np.abs(p).max()):
            raise NotImplementedError(
                f"eps_{c} is transversely non-uniform on the perimeter ring of the TFSF box "
                f"(largest deviation {dev:.3e}): a structure cuts into a box face or crosses "
                "the box column, and the incident field cannot be defined by a 1D background")
    for kk in (alo - 1, alo, ahi, ahi + 1):            # the whole footprint of the a faces
        inner = kk in (alo, ahi)
        for c, em in em_all.items():
            face = em[ulo - 1:uhi + 1, vlo - 1:vhi + 1, kk]
            dev = float(np.abs(face - em[ulo, vlo, kk]).max())
            if dev > TFSF_EPS_RTOL:
                if inner:
                    # The layer just inside the box (total-field region): placing a scatterer
                    # flush against the box face is a legal layout (PhaseChangeAntennas: the
                    # antenna's bottom face = the box bottom = the substrate surface). The
                    # incident field is defined by the 1D background, so non-uniformity here
                    # goes into the scattered field and does not break the TF/SF split.
                    print(f"[openem] inside the TFSF injection-axis face ({'xyz'[a]}={kk}): "
                          f"eps_{c} is non-uniform (deviation {dev:.2e}), treating it as a "
                          f"flush scatterer")
                    continue
                raise NotImplementedError(
                    f"outside the TFSF injection-axis face ({'xyz'[a]}={kk}), eps_{c} is "
                    "non-uniform: a structure cuts into the injection plane (scattered-field "
                    "region)")
    if med.sigma["ex"] is not None:
        # A lossy background works now (the 1D line uses the two-coefficient ca/cb form), but
        # the σ profile along the injection axis has to be **transversely uniform**: there is
        # only one 1D line, and without transverse uniformity there is no single incident field.
        # The criterion is the same as for ε: compare the perimeter ring against a reference
        # column.
        for c, s in med.sigma.items():
            sm = bx.uva(s)
            # **Only the perimeter ring of the shell** (the uu/vv from _tfsf_ring), not the
            # whole box volume: the TF/SF correction acts only on the shell, and however large
            # σ is inside the box has nothing to do with the incident field.
            ring = sm[uu, vv, :][:, c0:c1 + 1]
            ref = sm[uu[0], vv[0], c0:c1 + 1]
            dev = float(np.abs(ring - ref).max())
            scale = max(float(np.abs(ref).max()), 1.0)
            if dev / scale > TFSF_EPS_RTOL:
                raise NotImplementedError(
                    f"sigma_{c} is transversely non-uniform on the TFSF box shell (relative "
                    f"deviation {dev / scale:.3e}): one 1D line cannot define a unique "
                    "incident field")
    if med.disp_any:
        shell = set()
        idx3 = [0, 0, 0]
        for i, j in zip(uu, vv):
            for kk in range(c0, c1 + 1):
                idx3[u], idx3[v], idx3[a] = int(i), int(j), int(kk)
                shell.add(_flat(tuple(idx3), grid.shape))
        if shell & med.disp_any:
            raise NotImplementedError("the TFSF box shell or box column lies in a dispersive "
                                      "medium: not supported yet")


def _tfsf_1d_line(bx: _TfsfBox, eps_col, sig_col, dl_col, runway):
    """The 1D auxiliary grid ``(dl1, eps1, sigma1)``: the box column copies the 3D dl, with a
    long enough runway attached at each end.

    eps1[q] is the ε of edge q; the box column segment copies 3D edges c0..c1+1 and the runway
    extends the end values. The σ profile uses the same slicing and end-value extension as eps1;
    in a lossless scene it is None and the 1D line takes the original path.
    ``runway`` is the number of runway cells at the two ends, ``(n_lo, n_hi)``.
    """
    n_lo, n_hi = runway
    dl1 = np.concatenate([np.full(n_lo, dl_col[0]), dl_col, np.full(n_hi, dl_col[-1])])
    eps_seg = _column_segment(eps_col, bx.c0, bx.c1)
    eps1 = np.concatenate([np.full(n_lo, eps_seg[0]), eps_seg,
                           np.full(n_hi - 1, eps_seg[-1])])
    sigma1 = None
    if sig_col is not None:
        sig_seg = _column_segment(sig_col, bx.c0, bx.c1)
        if np.abs(sig_seg).max() > 0.0:
            sigma1 = np.concatenate([np.full(n_lo, sig_seg[0]), sig_seg,
                                     np.full(n_hi - 1, sig_seg[-1])])
    return dl1, eps1, sigma1


def _tfsf_injection(bx: _TfsfBox, direction, col0, n_lo, eps_col):
    """The injection plane ``inj`` (a 1D index) and the ε segment ``seg`` required to be uniform.

    The injection plane sits TFSF_GAP//2 cells outside the box face, inside the uniform stretch
    on the injection side. Uniformity only has to cover **the injection point and one cell on
    either side**, which is where the analytic injection tables are seeded; the remaining cells
    between the injection point and the box face are propagated by the 1D line itself along the
    true ε profile, so an interface (substrate/air) is allowed to fall there. The box bottom of
    PhaseChangeAntennas is exactly the substrate surface, and Tidy3D likewise relies on the 1D
    auxiliary solution to handle a layered background. Requiring the whole stretch from the end
    of the box column to the box face to be uniform, as it used to, was too strict.
    """
    if direction > 0:
        inj = col0 - TFSF_GAP // 2
        seg = eps_col[bx.c0:bx.alo - TFSF_GAP // 2 + 2]          # c0 .. injection cell+1
    else:
        inj = n_lo + (bx.ahi - bx.c0) + TFSF_GAP // 2
        seg = eps_col[bx.ahi + TFSF_GAP // 2 - 1:bx.c1 + 2]      # injection cell-1 .. c1+1
    return inj, seg


def _tfsf_normal_box(src, grid: Grid) -> _TfsfBox:
    """The box geometry for normal-incidence TFSF, together with the two margin checks.

    Raises:
        NotImplementedError: The box has infinite size on some axis (an open axis). Normal
            incidence goes through the 16 hand-written terms in
            ``sources_setup._tfsf_face_corrections`` (integer constants kept so the results stay
            bitwise identical), which have not been split up per face yet; the oblique path
            (angle_theta != 0) already supports it.
        ValueError: The box touches the grid boundary (each of the six face corrections needs
            one cell of margin), or the 1D box column runs past the injection axis.
    """
    a = int(src.injection_axis)
    lo_e, hi_e, open_axes = tfsf_box_edges(src, grid)
    if open_axes:
        raise NotImplementedError(
            f"the normal-incidence TFSF box has infinite size on the "
            f"{'/'.join('xyz'[q] for q in open_axes)} axis, which is not supported yet: the 16 "
            "hand-written face corrections have not been split up per face. The oblique path "
            "(angle_theta != 0) already supports it.")
    shp = grid.shape
    if not all(1 <= lo_e[q] < hi_e[q] <= shp[q] - 1 for q in range(3)):
        raise ValueError(
            f"the TFSF box {lo_e}..{hi_e} touches or exceeds the grid boundary {shp}: each of "
            "the six face corrections needs one cell of margin")
    bx = _TfsfBox(a, lo_e, hi_e)
    if bx.c0 < 0 or bx.c1 > shp[a] - 1:                # 3D cells covered by the 1D box column
        raise ValueError(
            f"the TFSF box column [{bx.c0}, {bx.c1}] runs past the {'xyz'[a]} axis "
            f"(n={shp[a]})")
    return bx


def tfsf(sim: td.Simulation, src, grid: Grid, med: MediaTables) -> TFSFSource:
    """``td.TFSF`` -> :class:`TFSFSource`.

    Any injection axis, any in-plane polarization (``pol_angle``/``angle_phi`` are folded into
    ``(β_u, β_v)`` by tidy3d's ``_pol_vector``). The incident field goes through the 1D
    auxiliary grid (:mod:`openem.tfsf1d`), which supports a substrate layered along the
    injection axis and cutting through the box; the background has to be **transversely
    uniform** on the perimeter ring of the box, otherwise the incident field cannot be defined
    and this fails closed.

    Below, the injection axis is written a and the transverse axes are ``u=(a+1)%3`` and
    ``v=(a+2)%3`` (cyclic right-handed order); for array indexing, ``np.moveaxis`` moves
    (u, v, a) to (0, 1, 2) so one single set of formulas applies.
    """
    a = int(src.injection_axis)
    if abs(float(src.angle_theta)) > 0:
        return _tfsf_oblique(sim, src, grid, med)
    pol = np.asarray(src._pol_vector, dtype=np.float64)
    if abs(pol[a]) > 1e-12:
        raise AssertionError(
            f"at normal incidence the polarization vector has an injection-axis component: "
            f"{pol}")
    u, v = cyclic_axes(a)
    beta_u, beta_v = float(pol[u]), float(pol[v])
    direction = sign_of(src.direction)

    ax = [grid.axes[q] for q in range(3)]
    bx = _tfsf_normal_box(src, grid)
    na = grid.shape[a]

    # ---- Background checks: the perimeter ring plus the footprints of the two a faces have to
    # be transversely uniform ----
    _tfsf_check_background(med, bx, grid)

    # ---- 1D auxiliary grid: the box column copies the 3D dl, with a long enough runway at
    # each end ----
    # The 1D line solves the E_u polarization, so ε is taken from the profile of component u.
    # When β_v != 0 the two polarizations share one line, which requires the ε_u and ε_v
    # profiles to agree; if they do not (an anisotropic background), the incident field would
    # need two lines, which is not implemented, so fail closed.
    key_u, key_v = E_KEYS[u], E_KEYS[v]
    eps_col = bx.uva(med.eps[key_u])[bx.ulo, bx.vlo, :].astype(np.float64)
    if abs(beta_v) > 0:
        eps_col_v = bx.uva(med.eps[key_v])[bx.ulo, bx.vlo, :].astype(np.float64)
        dev = np.abs(eps_col_v[bx.c0:bx.c1 + 2 if bx.c1 + 1 < na else None]
                     - eps_col[bx.c0:bx.c1 + 2 if bx.c1 + 1 < na else None]).max()
        if dev > TFSF_EPS_RTOL:
            raise NotImplementedError(
                f"on the box column of a skew-polarized TFSF, the ε_{key_u} and ε_{key_v} "
                f"profiles disagree (largest deviation {dev:.3e}): each polarization would "
                "need its own 1D incident line, which is not implemented")
    dl_col = ax[a].dl[bx.c0:bx.c1 + 1].astype(np.float64)
    n_steps = int(sim.num_time_steps)
    dt = float(sim.dt)
    freq0 = float(src.source_time.freq0)
    n_lo = tfsf1d.runway_cells(n_steps, dt, float(dl_col[0]))
    n_hi = tfsf1d.runway_cells(n_steps, dt, float(dl_col[-1]))
    sig_col = (bx.uva(med.sigma[key_u])[bx.ulo, bx.vlo, :].astype(np.float64)
               if med.sigma.get(key_u) is not None else None)
    dl1, eps1, sigma1 = _tfsf_1d_line(bx, eps_col, sig_col, dl_col, (n_lo, n_hi))
    col0 = n_lo + (bx.alo - bx.c0)                      # 1D index of 3D edge alo
    inj, seg = _tfsf_injection(bx, direction, col0, n_lo, eps_col)
    if np.abs(seg - seg[0]).max() > TFSF_EPS_RTOL:
        raise NotImplementedError(
            "the ε near the TFSF injection point (from injection cell ±1 to the end of the box "
            f"column) is non-uniform along the injection axis: the analytic injection tables "
            f"require local uniformity at the injection point, got {seg.min():.6g}.."
            f"{seg.max():.6g}")
    kh1 = inj - 1 if direction > 0 else inj
    ex_inc, hy_inc = incident_tables(
        source_time=src.source_time, dt=dt, num_steps=n_steps,
        dl_half=float(dl1[kh1]) / 2.0, eps_r=float(eps1[inj]),
        direction=direction, sigma=0.0, freq0=freq0)

    return TFSFSource(
        **_tfsf_common(src, bx, _Aux1D(dl1, eps1, sigma1, inj, col0, ex_inc, hy_inc),
                       dt, n_steps),
        direction=direction, pol_u=beta_u, pol_v=beta_v,
    )


def _beamish(src) -> bool:
    """Whether this source is injected through the beam carrier (a ModeSource-style 2D complex
    field): ModeSource, GaussianBeam, and any PlaneWave at oblique incidence or with a non-z
    normal."""
    if type(src).__name__ == "ModeSource" or isinstance(src, td.GaussianBeam):
        return True
    if isinstance(src, td.PlaneWave):
        return (abs(float(src.angle_theta)) > 0
                or abs(float(getattr(src, "angle_phi", 0.0))) > 0
                or normal_axis(src.size) != 2)
    return False


class BuiltSources(NamedTuple):
    """The four kinds of source from :func:`build`, already grouped by Scene's fields. The
    placement check (``_check_placements`` in scene/build.py) looks at all four, so the whole
    group is passed together."""

    planes: list[PlaneWaveSource]
    dipoles: list[PointDipole]
    mode_srcs: list
    tfsf_srcs: list[TFSFSource]


def build(sim: td.Simulation, grid: Grid, med: MediaTables) -> BuiltSources:
    """All sources, dispatched by type.

    Raises:
        ValueError: The scene contains no source at all.
        NotImplementedError: A source type that is not supported yet.
    """
    planes: list[PlaneWaveSource] = []
    dipoles: list[PointDipole] = []
    mode_srcs: list = []
    tfsf_srcs: list[TFSFSource] = []

    # BUGFIX(2026-09-09): when one scene mixes both injection conventions, the 1D auxiliary
    # grid plane wave and the beam carrier, the readout side has only one normalize.field_norm
    # divisor and the two conventions cannot both be right (the beam carrier is normalized to
    # unit power, divisor D; the plane wave to unit E amplitude, divisor D·√(nA/2η₀)).
    # On top of that, normalize.source_table takes sc.sources[0] (which holds only the
    # normal-incidence ones) while nb.backend._adjoint_coefficient multiplies back
    # sim.sources[normalize_index]; in a multi-source scene those are simply not the same
    # source, and they differ by an arbitrary complex constant.
    # So as soon as a beam carrier can appear in the scene, normal-incidence plane waves are
    # sent through the beam carrier as well: the convention becomes uniformly the mode branch,
    # and mode_sources[0] corresponds one to one with sim.sources[0].
    # (If the current source is at oblique incidence or has a non-z normal it is a beam carrier
    # itself and _mixed is necessarily True, so the PlaneWave branch below only needs to test
    # _mixed.)
    _mixed = any(_beamish(s) for s in sim.sources)
    for src in sim.sources:
        if isinstance(src, td.TFSF):
            tfsf_srcs.append(tfsf(sim, src, grid, med))
        elif isinstance(src, td.PointDipole):
            dipoles.append(dipole(sim, src, grid))
        elif isinstance(src, td.UniformCurrentSource):
            dipoles.append(uniform_current(sim, src, grid))
        elif type(src).__name__ == "CustomCurrentSource":
            # this is tidy3d's adjoint source (built by web/api/autograd)
            dipoles.extend(custom_current(sim, src, grid))
        elif type(src).__name__ == "CustomFieldSource":
            # tangential E/H on a plane -> equivalent currents (J = n̂ × H, M = -n̂ × E),
            # reusing custom_current
            dipoles.extend(custom_field(sim, src, grid))
        elif isinstance(src, td.PlaneWave):
            # The kernel of the TF/SF 1D auxiliary grid path only supports propagation along z
            # (the src.axis != 2 guard in sources_setup._plane_source). Plane waves with an x/y
            # normal, which is exactly what the adjoint source of a diffraction objective looks
            # like (Autograd4MultiObjective: an x-normal diffraction monitor -> a backward plane
            # wave propagating along x), go through the same ModeSource-style carrier, which
            # picks its axis from normal_axis(src.size) throughout and holds for all three axes.
            if _mixed:
                # Oblique incidence: the transverse phase goes through the ModeSource-style 2D
                # complex field carrier (PlaneWaveBeamProfile), which needs Bloch boundaries and
                # the complex solver. Normal incidence on a z normal still uses the 1D auxiliary
                # grid, which is cheaper.
                mode_srcs.append(modes_mod.plane_wave_beam(sim, src, grid))
            else:
                planes.append(plane_wave(sim, src, grid, med))
        elif type(src).__name__ == "ModeSource":
            mode_srcs.append(modes_mod.mode_source(sim, src, grid))
        elif isinstance(src, td.GaussianBeam):
            mode_srcs.append(modes_mod.gaussian_beam(sim, src, grid))
        else:
            raise NotImplementedError(
                f"only PlaneWave / PointDipole / UniformCurrentSource / "
                f"CustomCurrentSource / CustomFieldSource / "
                f"ModeSource / GaussianBeam / TFSF are supported, got {type(src).__name__}"
            )
    if not planes and not dipoles and not mode_srcs and not tfsf_srcs:
        raise ValueError("the scene contains no sources")
    return BuiltSources(planes, dipoles, mode_srcs, tfsf_srcs)
