# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Mingdian Liu
"""Data model for the solver. **Does not import tidy3d.**

This is a separate file in order to hold that boundary: ``scene/`` is the builder and must
import tidy3d, while ``Scene`` itself is nothing but numpy arrays, and the solve path
(solver / flux / serialize) does not need tidy3d. Putting both in one file would make
``import openem.serialize`` crash inside a solve job container, where tidy3d is not
necessarily installable. ``tests/test_no_tidy3d_leak.py`` locks this down.
"""

from __future__ import annotations

from dataclasses import dataclass, field, replace

import numpy as np

from openem import cpml, waveform
from openem.grid import Grid, ravel_cells, unravel_cells


@dataclass(frozen=True)
class PlaneWaveSource:
    """One-way plane wave source.

    Attributes:
        axis: Propagation axis (0/1/2).
        plane_index: Cell index of the source plane.
        direction: ``+1`` or ``-1``.
        pol_axis: Polarization axis of E (0/1/2).
        waveform: Pre-sampled time waveform.
        angle_theta: Polar angle of incidence, radians. Only 0 (normal incidence) is supported.
        num_freqs: Number of frequency points Tidy3D uses to approximate **the variation of the
            source's spatial profile with frequency** (default 1). It is only needed at oblique
            incidence, where the transverse wavenumber is ∝ ω; **at normal incidence it has no
            effect** (any value gives the same result), so we record it but never use it.
        ex_inc: ``(N+1,)`` ``Ex_inc(n*dt)`` on the source plane, for the TF/SF H correction.
        hy_inc: ``(N+1,)`` ``Hy_inc((n+0.5)*dt + half-cell delay) / eta``, for the E correction.
        incident: The raw parameters needed to recompute ``ex_inc/hy_inc``. **Only for
            :func:`lossless`**: both tables are pinned at scene build time to the ``(ε, σ)`` at
            the source plane, so changing the medium forces a recompute. Otherwise the reference
            run injects an H matched to a lossy medium into a lossless one.
    """

    axis: int
    plane_index: int
    direction: int
    pol_axis: int
    waveform: waveform.Waveform
    angle_theta: float
    num_freqs: int
    ex_inc: np.ndarray
    hy_inc: np.ndarray
    incident: dict | None = None
    #: Source center coordinate along the propagation axis as given by the user (SI m). Tidy3D
    #: references the plane wave phase to it; once the source plane snaps to the grid line
    #: edges[plane_index], the readout side has to add e^{ik*dir*(center-z_e)}
    #: (normalize._src_snap_phase). Old scene.npz files lack this field -> nan -> no correction.
    center: float = float("nan")


@dataclass(frozen=True)
class TFSFSource:
    """TFSF box: the incident plane wave exists only inside the box; all six faces get a
    correction.

    The incident field comes from a **1D auxiliary grid** (:mod:`openem.tfsf1d`): a 1D Yee line
    with the same dt, the same dl and the same ε profile as the 3D injection axis. The incident
    field satisfies the 3D discrete equations bitwise, so leakage outside the box is down to
    float32 rounding alone, and when a substrate crosses the box (MultipoleExpansion) the
    transmission and reflection at the interface also come out right on their own. The injection
    axis is arbitrary, as is the in-plane polarization (``pol_u/pol_v``). At **oblique
    incidence** (``angle_theta ≠ 0``) the 1D line runs along k̂ instead and the table column
    positions become real-valued; see :mod:`openem.tfsf_oblique`.

    Attributes:
        axis: Injection axis 0/1/2. The in-plane polarization basis follows the cyclic
            right-handed order ``u=(axis+1)%3``, ``v=(axis+2)%3``.
        pol_u / pol_v: Components of the E polarization vector on ``(û, v̂)`` (tidy3d's
            ``_pol_vector``, which includes ``pol_angle`` and ``angle_phi``).
        direction: ``+1`` / ``-1``, propagation direction of the incident wave.
        box_lo / box_hi: Primal edge indices of the box's six faces, ``(ilo, jlo, klo)`` /
            ``(ihi, jhi, khi)``. The interior (total-field region) is cells ``ilo..ihi-1``, and
            so on.
        dl1: ``(n1,)`` primal spacing of the 1D grid [m]. Over the box columns it copies the dl
            of the 3D injection axis; both ends are extended by a long enough runway
            (tfsf1d.runway_cells).
        eps1: ``(n1,)`` relative permittivity at the 1D edge positions (the ε of component u on
            the box shell columns, the runway extended with the end values; when β_v≠0 the
            builder requires the ε_u and ε_v profiles to agree).
        inj: 1D injection plane (E edge index).
        col0: 1D index corresponding to the 3D edge ``box_lo[axis]``.
        ex_inc / hy_inc: ``(N+1,)`` 1D injection tables (sources._incident_tables, with the wave
            impedance set by the ε at the injection plane).
        waveform: Time waveform (used by shutoff's source_end_step and by source spectrum
            normalization).
        num_freqs: Recorded but unused: it has no effect at normal incidence.
    """

    direction: int
    box_lo: tuple[int, int, int]
    box_hi: tuple[int, int, int]
    dl1: np.ndarray
    eps1: np.ndarray
    inj: int
    col0: int
    ex_inc: np.ndarray
    hy_inc: np.ndarray
    waveform: waveform.Waveform
    num_freqs: int = 1
    axis: int = 2
    pol_u: float = 1.0
    pol_v: float = 0.0
    #: ``(n1,)`` conductivity at the 1D edge positions [S/m]; None for a lossless background.
    #: When present, the 1D line uses the two-coefficient lossy update (tfsf1d.tables).
    sigma1: np.ndarray | None = None
    #: The next five are nonzero only at **oblique incidence** (``angle_theta ≠ 0``):
    #: ``k_hat``/``e_hat`` are the unit propagation direction and the E polarization direction
    #: in xyz (tidy3d's ``_dir_vector`` / ``_pol_vector``), ``dl1_step`` is the 1D cell spacing
    #: [m] obtained by matching the discrete dispersion, ``proj0`` is the projected value [m]
    #: mapped onto 1D edge ``col0``, and ``n_col_e`` is the number of columns in the E table.
    #: At normal incidence ``k_hat is None`` and we take the hand-written 16 terms in
    #: sources_setup (bitwise identical results).
    k_hat: tuple[float, float, float] | None = None
    e_hat: tuple[float, float, float] | None = None
    dl1_step: float = 0.0
    proj0: float = 0.0
    n_col_e: int = 0
    #: "Open axes": the source is infinite along such an axis, the box spans the whole domain,
    #: and those two faces do not exist (tfsf2 in TFSF.ipynb has size=(inf, inf, 4), so both x
    #: and y are open). An empty tuple means all six faces are present.
    open_axes: tuple[int, ...] = ()
    #: Factor between the actual injected amplitude and ``ex_inc/hy_inc``. At oblique incidence
    #: it is ``1/√cosθ`` (Tidy3D's solver normalizes power along the **injection axis normal**;
    #: see scene/sources._tfsf_oblique); at normal incidence it is always 1.0 (multiplying by
    #: 1.0 is exact in IEEE754, so results stay bitwise identical). The two tables keep the
    #: reference amplitude, because they are also the normalization divisor
    #: (normalize.source_table).
    amp_scale: float = 1.0


@dataclass(frozen=True)
class ModeSource:
    """Guided mode source, TF/SF injection on a single plane.

    Same TF/SF machinery as :class:`PlaneWaveSource`, with only two differences: the incident
    field goes from a **constant** to a **2D complex mode profile**, and the tangential
    components go from one pair (Ex, Hy) to two pairs.

    The profile comes from the client-side ``ModeSolver`` (its ``n_eff`` is bitwise identical to
    the reference implementation's, and modes are normalized to unit power). It already covers
    the whole transverse plane, zero-padded outside the mode box. Each of the four components
    sits at **its own Yee transverse position**; nothing is interpolated.

    The tangential pair follows the **cyclic order** ``(t1, t2) = ((axis+1)%3, (axis+2)%3)``:
    the Yee curl keeps its form under a cyclic permutation, so a single sign table for the four
    correction terms (``sources_setup._MODE_SIGNS``) serves all three normals. Field names such
    as ``ey_inc``/``ez_inc`` keep the x-normal naming but really mean **E_t1/E_t2** (for an x
    normal, t1 happens to be y and t2 z).

    Attributes:
        axis: Normal axis 0/1/2.
        plane_index: Index of the source plane along the normal axis.
        direction: ``+1`` / ``-1``.
        ey_inc, ez_inc: ``(n_t1, n_t2)`` complex profile (E_t1/E_t2), on the source plane,
            sampled at whole steps.
        hy_inc, hz_inc: ``(n_t1, n_t2)`` complex profile (H_t1/H_t2), on the neighboring H
            plane (half a cell away, with the phase already folded in), sampled at half steps.
        amp_e: ``(N+1,)`` ``amp_time(n·dt)``, complex.
        amp_h: ``(N+1,)`` ``amp_time((n+0.5)·dt)``, complex.
        n_eff: Effective index of the mode; needed both for the half-cell phase and for the
            acceptance criterion.
    """

    plane_index: int
    direction: int
    ey_inc: np.ndarray
    ez_inc: np.ndarray
    hy_inc: np.ndarray
    hz_inc: np.ndarray
    amp_e: np.ndarray
    amp_h: np.ndarray
    n_eff: complex
    axis: int = 0
    #: Broadband correction terms (present only when ``num_freqs > 1`` and a single profile is
    #: not enough; all None otherwise). The profiles are ``(K, n_t1, n_t2)`` arrays of
    #: **Pₖ − P₀**, the coefficients ``(K, N+1)`` arrays of aₖ, so the injection is
    #: P₀·S(t) + Σₖ (Pₖ−P₀)·aₖ(t). See the derivation in scene/modes.py.
    bb_ey: np.ndarray | None = None
    bb_ez: np.ndarray | None = None
    bb_hy: np.ndarray | None = None
    bb_hz: np.ndarray | None = None
    bb_amp_e: np.ndarray | None = None
    bb_amp_h: np.ndarray | None = None


@dataclass(frozen=True)
class PointDipole:
    """Point current source, spread over the surrounding Yee points.

    Tidy3D's ``PointDipole`` is a ``ReverseInterpolatedSource``: with ``interpolate=True``, the
    "equivalent source data is applied on the surrounding Yee grid points ... using linear
    interpolation" (``source/current.py:73-80``). In the validation set, **not one** dipole
    lands exactly on the Yee point of its own component, so the interpolation has to be done;
    but at most two axes need it (z always hits exactly), so the number of weights is 2 or 4,
    and OpenEM's "four-edge" stencil is just enough for these scenes. Written here as general
    trilinear interpolation, with no such limit baked in.

    Attributes:
        component: 0/1/2 -> x/y/z component (E or H depending on ``magnetic``).
        indices: ``(nw, 3)`` cell indices of the stencil points.
        coef: ``(nw,)`` ``weight / dV``: the multilinear weight divided by the Yee dual volume
            at that point. In Tidy3D ``amplitude=1`` is a **current moment** of 1 A·µm, so the
            current density on the grid is ``amp/dV``. (This one is derived, not copied.)
        waveform: The source's own waveform. **Every source must be sampled separately**: in the
            multi-source scenes of the validation set (ResonanceFinder / Bandstructure) the 7
            sources share freq0/fwidth but each has a different ``phase``.
        magnetic: True means a magnetic dipole (polarization=Hx/Hy/Hz): it is injected into the
            H update, uses the Yee coordinates and dual volume of H, ``amplitude=1`` is by
            duality a **magnetic current moment** of 1 V·µm, and the amplitude is sampled at
            whole steps ``n·dt`` (the midpoint of an H step).
    """

    component: int
    indices: np.ndarray
    coef: np.ndarray
    waveform: waveform.Waveform
    magnetic: bool = False


@dataclass(frozen=True)
class FieldMonitor:
    """Frequency-domain field monitor. **Stores phasors at their Yee positions; no colocation
    inside the kernel.**

    Colocation is a linear operation and commutes with the DFT, so it is done on the host: the
    convention is still being cross-checked, and in Python it is cheap to change how it is done.

    Attributes:
        name: Monitor name in Tidy3D.
        freqs: Frequencies [Hz].
        origin: Starting cell index of the storage box, ``(i0, j0, k0)``.
        interval_space: Tidy3D's spatial downsampling stride. **The solver does not
            downsample**: the box is accumulated at full resolution as usual (a superset, so
            nothing is lost), and the downsampling is applied at write-back / comparison time
            using Tidy3D's own choice of points. What Tidy3D keeps is a subset of those same Yee
            points, so taking the subset afterwards is exact.
        box: Size of the storage box in cells, ``(ni, nj, nk)``. One spare cell is already
            reserved at the low end of each axis, because colocating onto a grid point needs
            both ``F[i-1]`` and ``F[i]``.
        apodization: Apodization window ``(start, end, width)`` (seconds; start/end may be
            None), ``None`` means no apodization. For the window function see
            ``apodization.window``.
    """

    name: str
    freqs: np.ndarray
    origin: tuple[int, int, int]
    box: tuple[int, int, int]
    apodization: tuple[float | None, float | None, float] | None = None
    interval_space: tuple[int, int, int] = (1, 1, 1)


#: Order of the field components in the arrays, shared by kernel and host.
COMPONENTS = ("Ex", "Ey", "Ez", "Hx", "Hy", "Hz")

#: Name suffix of the internal FieldMonitor that a ModeMonitor expands into (the name is built
#: by scene/modes.plane_monitor; td_readout.finalize uses it to skip those monitors, and
#: projection / mode decomposition use it to recognize a plane holding fields at their Yee
#: positions). It lives here because td_readout cannot import the scene package (no tidy3d in
#: the job container); scene/modes.PLANE_SUFFIX and td_readout._PLANE_SUFFIX are aliases of it.
MODE_PLANE_SUFFIX = "__mode_plane"


@dataclass(frozen=True)
class FieldTimeMonitor:
    """Time-domain field monitor. Like :class:`FieldMonitor`, it stores raw values **at their
    Yee positions**.

    The set of sampled steps is ``m = step_begin, +interval, ...`` (``m < step_end``), at times
    ``m·dt``. The interval comes from Tidy3D's own ``TimeMonitor.time_inds(sim.tmesh)``
    (``monitor.py:331``); we do not reimplement that half-open interval logic.

    Attributes:
        name: Monitor name in Tidy3D.
        comps: Indices into :data:`COMPONENTS` of the components to record. Tidy3D's ``fields``
            is a subset: ResonanceFinder only asks for ``Ey``.
        origin / box: Starting cell index and size in cells of the storage box.
        step_begin / step_end / interval: Sampling window.
        num_slots: Nominal number of slots. Early termination leaves it partly unfilled; for the
            actual count see :attr:`solver.Result.time_slots`.
    """

    name: str
    comps: tuple[int, ...]
    origin: tuple[int, int, int]
    box: tuple[int, int, int]
    step_begin: int
    step_end: int
    interval: int
    num_slots: int

    def slot_of(self, m: int) -> int | None:
        """Which slot time step ``m`` writes into; None if it is not a sampled step.
        (The solver's own slot logic lives in the kernel and in readout._readout; this helper is
        only used by the tests.)"""
        if m < self.step_begin or m >= self.step_end:
            return None
        off = m - self.step_begin
        if off % self.interval:
            return None
        s = off // self.interval
        return s if s < self.num_slots else None


@dataclass(frozen=True)
class FluxTimeMonitor:
    """Time-domain flux monitor: **one scalar** per sampled step, the instantaneous Poynting
    flux through the plane.

    The spatial convention is **colocated integration** (``use_colocated_integration=True``,
    hard-coded for every FluxTimeMonitor): the four tangential components are linearly
    interpolated onto "tangential primal boundary points x the monitor's exact normal
    coordinate", and the area element is the dual cell spanned by the midpoints of neighboring
    samples, with trapezoidal weights truncated exactly at the monitor bounds. This is **not**
    the same convention as the staggered Yee area element (the ``flux.area_weights`` path used
    by :class:`FluxMonitor`). For the full geometry see ``flux.flux_time_geometry``. The time
    convention matches :class:`FieldTimeMonitor` (the ``time_inds`` window plus the whole-step
    time average of H, ``H̄ = ½(H^{m−1/2}+H^{m+1/2})``). The instantaneous flux has no ½ like the
    frequency domain does; that factor comes from the CW time average.

    Its relation to the frequency-domain flux is **Parseval**, not equality frequency by
    frequency: ``Σ_m s(m) = 4·dt·∫₀^∞ P(f) df`` (the DFT of s is a convolution of the E and H
    spectra and has no pointwise relation to ½Re(E×H*) at the monitor frequencies).

    Attributes:
        name: Monitor name in Tidy3D.
        axis / plane_index / normal_dir: As in :class:`FluxMonitor`.
        frac: Linear interpolation weight ∈ [0, 1) of the plane between ``edges[plane_index]``
            and ``edges[plane_index+1]``. What gets interpolated is the **flux scalar** (each of
            the two neighboring primal planes is integrated first, then
            ``(1−frac)·S_km + frac·S_{km+1}``), not the field; that was settled by a
            self-consistency check against the reference outputs. All 6 planes of
            NanobeamCavity lie off the Yee boundaries.
        t_bounds: The monitor's ``((lo, hi), (lo, hi))`` on the two tangential axes (ascending),
            in **meters**. The area element is truncated exactly at these bounds, and samples
            within half a dual cell outside them still get a partial weight, which is why
            coordinates are stored rather than cell indices.
        step_begin / step_end / interval / num_slots: As in :class:`FieldTimeMonitor`. For how
            many slots were actually filled, see the length of ``solver.Result.flux_time``.
    """

    name: str
    axis: int
    plane_index: int
    normal_dir: int
    step_begin: int
    step_end: int
    interval: int
    num_slots: int
    t_bounds: tuple[tuple[float, float], tuple[float, float]] = ((0.0, 0.0), (0.0, 0.0))
    frac: float = 0.0


@dataclass(frozen=True)
class PermittivityMonitor:
    """Frequency-domain monitor for the diagonal entries of the ε tensor.

    Tidy3D's docstring states explicitly that the values live at **Yee grid locations**
    (``monitor.py:1270-1274``), so it follows the same convention as our internal
    ``eps_ex/ey/ez``: no need to run the solver, just slice it out of :meth:`Scene.eps_box`.

    Attributes:
        name: Monitor name in Tidy3D.
        freqs: Frequencies [Hz]. For a non-dispersive material ε does not depend on frequency,
            so this value is only used for reporting.
        origin / box: Starting cell index and size in cells of the storage box.
    """

    name: str
    freqs: np.ndarray
    origin: tuple[int, int, int]
    box: tuple[int, int, int]


@dataclass(frozen=True)
class FluxMonitor:
    """Planar flux monitor.

    Attributes:
        name: Monitor name in Tidy3D; results must be written back under the same name.
        axis: Normal axis of the plane.
        plane_index: Cell index of the plane.
        normal_dir: ``+1`` / ``-1``, Tidy3D's ``normal_dir``.
        freqs: Frequencies, Hz.
        source_spectrum: Source spectrum, same length as ``freqs``, used to normalize the flux.
        transverse: Cell index ranges on the two transverse axes, ``((j0, j1), (k0, k1))``
            (closed intervals); ``None`` means the monitor covers the whole plane. **Required
            whenever the monitor does not span the plane transversely**: without it the whole
            plane is integrated. The two BraggGratings monitors both sit at
            ``plane_index=1680`` and differ only in their y center, so with this field missing
            they return **exactly the same** numbers.
        apodization: Apodization window ``(start, end, width)``, same convention as
            :class:`FieldMonitor`.
    """

    name: str
    axis: int
    plane_index: int
    normal_dir: int
    freqs: np.ndarray
    source_spectrum: np.ndarray
    transverse: tuple[tuple[int, int], tuple[int, int]] | None = None
    apodization: tuple[float | None, float | None, float] | None = None
    #: Tangential physical bounds (meters) ``((lo1, hi1), (lo2, hi2))``, used by the colocated
    #: integration convention (2026-09-06); None means an old result or old scene, which can
    #: only integrate at the Yee positions.
    t_bounds: tuple[tuple[float, float], tuple[float, float]] | None = None


def _check_csr(cls_name: str, shape: tuple[int, int, int], comp: np.ndarray,
               cell: np.ndarray, ofs_list: list[tuple[str, np.ndarray, int]]) -> np.ndarray:
    """Structural check shared by the dispersion tables; returns the joint key of
    ``(comp, cell)``.

    Better to stop with an error here than to let the kernel run with out-of-range indices or
    duplicate entries: a duplicate lets a later entry overwrite the cell that writes hist / E,
    and the result is a run that completes but is physically wrong.
    """
    n = int(comp.size)
    ncell = int(np.prod(shape))
    if comp.shape != (n,) or cell.shape != (n,):
        raise ValueError(f"{cls_name}: comp/cell lengths differ")
    if n == 0:
        return np.zeros(0, dtype=np.int64)
    if comp.min() < 0 or comp.max() > 2:
        raise ValueError(f"{cls_name}: comp out of range {comp.min()}..{comp.max()}")
    if cell.min() < 0 or cell.max() >= ncell:
        raise ValueError(f"{cls_name}: cell out of range {cell.max()} >= {ncell}")
    for name, ofs, npole in ofs_list:
        if ofs.shape != (n + 1,):
            raise ValueError(f"{cls_name}.{name} has shape {ofs.shape}, expected ({n + 1},)")
        if ofs[0] != 0 or ofs[-1] != npole:
            raise ValueError(
                f"{cls_name}.{name} endpoints {ofs[0]}..{ofs[-1]}, expected 0..{npole}")
        if np.any(np.diff(ofs) < 0):
            raise ValueError(f"{cls_name}.{name} is not monotonic")
    key = comp.astype(np.int64) * ncell + cell
    # O(n) scatter to find duplicates. np.unique used to hash-deduplicate hundreds of millions
    # of entries: measured at 57 s on GroupDelay, the bulk of the initialization time.
    # key ∈ [0, 3*ncell) ⇒ one bool marker array of size 3*ncell is enough (about 56 MB).
    seen = np.zeros(3 * ncell, dtype=bool)
    seen[key] = True
    n_uniq = int(seen.sum())
    del seen
    if n_uniq != n:
        raise ValueError(f"{cls_name}: {n - n_uniq} duplicate (comp, cell) entries")
    return key


@dataclass(frozen=True)
class Dispersion:
    """Pole data for dispersive materials, stored compressed by "dispersive cell component".

    Tidy3D reduces every dispersive material to a PoleResidue, so there is only one model here.
    In the time domain each conjugate pair is one complex first-order equation
    ``dP/dt = qP + rE``, and the physical polarization is ``2Re(P)``. For the discretization and
    the derivation see ``kernels/dispersion.cu``.

    The layout is CSR-like: entry ``e`` occupies ``pole_ofs[e] : pole_ofs[e+1]`` of the
    coefficient arrays. **A given (comp, cell) appears in one entry only**: a mixed cell merges
    the poles of both materials into the same entry, so the kernel needs no atomics when writing
    hist. :meth:`validate` locks this down.

    Attributes:
        comp: ``(n_entry,)`` component, 0/1/2 -> Ex/Ey/Ez.
        cell: ``(n_entry,)`` flat cell index ``(i*ny + j)*nz + k``.
        pole_ofs: ``(n_entry + 1,)`` CSR offsets; the last entry equals the total pole count.
        am1: ``(n_pole,)`` complex **``A − 1``**, where ``A = (1+qΔt/2)/(1−qΔt/2)``. Storing the
            difference rather than ``A`` avoids float32 cancellation; see
            ``kernels/dispersion.cu``.
        b: ``(n_pole,)`` complex ``B = (rΔt/2)/(1 − qΔt/2)``, **already multiplied by the
            geometric weight of the cell** (the arithmetic average in a mixed cell scales the
            residues by the fill fraction and leaves the pole positions unchanged).
        g: ``(n_entry,)`` ``Σ_p 2Re(B_p)``, the implicit coupling term of the E update. It is
            stored separately because it has to be scattered back into a dense array by
            ``(comp, cell)`` before being merged with ``ca``/``cb``.
    """

    comp: np.ndarray
    cell: np.ndarray
    pole_ofs: np.ndarray
    am1: np.ndarray
    b: np.ndarray
    g: np.ndarray

    @property
    def n_entry(self) -> int:
        return int(self.comp.size)

    @property
    def n_pole(self) -> int:
        return int(self.am1.size)

    def validate(self, shape: tuple[int, int, int]) -> None:
        """Shapes, index ranges, uniqueness of (comp, cell), monotonicity of the CSR offsets.

        Better to stop with an error here than to let the kernel run with out-of-range indices
        or duplicate entries: a duplicate lets the entry written last overwrite hist, and the
        result is a run that completes but is physically wrong.
        """
        n = self.n_entry
        for name, arr in (("g", self.g, ), ("am1", self.am1,), ("b", self.b,)):
            want = n if name == "g" else self.n_pole
            if arr.shape != (want,):
                raise ValueError(f"Dispersion.{name} has shape {arr.shape}, expected ({want},)")
        _check_csr("Dispersion", shape, self.comp, self.cell,
                   [("pole_ofs", self.pole_ofs, self.n_pole)])

    def g_dense(self, shape: tuple[int, int, int]) -> tuple[np.ndarray, ...]:
        """Scatter ``g`` into three dense ``shape`` arrays, ready to merge into ``ca``/``cb``."""
        out = [np.zeros(int(np.prod(shape)), dtype=np.float64) for _ in range(3)]
        for c in range(3):
            m = self.comp == c
            out[c][self.cell[m]] = self.g[m]
        return tuple(o.reshape(shape) for o in out)


@dataclass(frozen=True)
class DispersionMix:
    """Cells on a slanted interface: a weighted sum of the arithmetic and harmonic averages
    (``β = n_i²``).

    ``ε_eff = (1−β)·Σfᵢεᵢ + β/(Σfᵢ/εᵢ)``. This **sum** cannot be expressed through a fixed set
    of poles of either ε or 1/ε, so the D of the harmonic part becomes a separate unknown and
    one scalar linear equation is solved per cell per step (derivation and kernel in
    ``kernels/dispersion_mix.cu``).

    **β=0 degenerates exactly into :class:`Dispersion`** (poles of ε), **β=1 into the pure
    harmonic branch** (poles of 1/ε). Each of the two limits has its own independent acceptance
    criterion, so this class needs no new reference of its own.

    Attributes:
        comp / cell: As in :class:`Dispersion`.
        p_ofs / pa / pb: CSR offsets and ``(A−1, B)`` of the arithmetic branch, from the poles
            of **ε**.
        q_ofs / qa / qb: CSR offsets and ``(A−1, B)`` of the harmonic branch, from the poles of
            **1/ε**.
        beta: ``(n_entry,)`` mixing coefficient.
        eps_inf: ``(n_entry,)`` ``Σfᵢε_∞ᵢ`` (high-frequency limit of the arithmetic branch).
        zeta_inf: ``(n_entry,)`` ``Σfᵢ/ε_∞ᵢ`` (high-frequency limit of the harmonic branch).
    """

    comp: np.ndarray
    cell: np.ndarray
    p_ofs: np.ndarray
    pa: np.ndarray
    pb: np.ndarray
    q_ofs: np.ndarray
    qa: np.ndarray
    qb: np.ndarray
    beta: np.ndarray
    eps_inf: np.ndarray
    zeta_inf: np.ndarray

    @property
    def n_entry(self) -> int:
        return int(self.comp.size)

    def coefficients(self) -> tuple[np.ndarray, ...]:
        """``(k1, k2, 1−β, 1/(k1k2+β))``: the four coefficients the kernel uses directly.

        ``G_P = Σ2Re(B_P)`` and ``G_Q = ΣRe(B_Q)`` are the implicit coupling of the two pole
        sets to time n+1.
        """
        gp = (np.add.reduceat(2.0 * self.pb.real, self.p_ofs[:-1])
              if self.pb.size else np.zeros(self.n_entry))
        gq = (np.add.reduceat(self.qb.real, self.q_ofs[:-1])
              if self.qb.size else np.zeros(self.n_entry))
        # For entries with no poles reduceat spills into the next segment; zero them explicitly
        gp = np.where(np.diff(self.p_ofs) > 0, gp, 0.0)
        gq = np.where(np.diff(self.q_ofs) > 0, gq, 0.0)
        one_mb = 1.0 - self.beta
        k1 = one_mb * (self.eps_inf + gp)
        k2 = self.zeta_inf + gq
        den = k1 * k2 + self.beta
        if np.any(np.abs(den) < 1e-30):
            raise ValueError("mixed-cell denominator k1·k2+β hit zero: coefficients unusable")
        return k1, k2, one_mb, 1.0 / den

    def validate(self, shape: tuple[int, int, int], *others) -> None:
        n = self.n_entry
        for name, arr in (("beta", self.beta), ("eps_inf", self.eps_inf),
                          ("zeta_inf", self.zeta_inf)):
            if arr.shape != (n,):
                raise ValueError(f"DispersionMix.{name} has shape {arr.shape}, expected ({n},)")
        key = _check_csr("DispersionMix", shape, self.comp, self.cell,
                         [("p_ofs", self.p_ofs, int(self.pa.size)),
                          ("q_ofs", self.q_ofs, int(self.qa.size))])
        if n == 0:
            return
        if np.any((self.beta < -1e-6) | (self.beta > 1 + 1e-6)):
            raise ValueError(
                f"β outside [0,1]: {self.beta.min():.4f}..{self.beta.max():.4f}"
                " (the square of the normal projection; out of range means an unreliable fit)")
        ncell = int(np.prod(shape))
        for o in others:
            if o is None or not o.n_entry:
                continue
            ok = o.comp.astype(np.int64) * ncell + o.cell
            if np.intersect1d(key, ok).size:
                raise ValueError("the same (comp, cell) appears in more than one dispersion path")
        self.coefficients()


@dataclass(frozen=True)
class TimeModulation:
    """Time-varying medium: ``δε(t) = amp · cos(2π·freq·t − phase)`` on the modulated cells.

    The mathematical form is copied from the tidy3d client's ``components/time_modulation.py``
    (not guessed): ``δε(r,t) = Re[amp_time(t)·amp_space(r)]``, where
    ``amp_time = A_t·exp(i·φ_t − i·2π·f·t)`` and ``amp_space = A_r·exp(i·φ_r)``, which expands
    to ``A_t·A_r·cos(2πft − φ_t − φ_r)``. Stored here are the combined
    ``amp = A_t·A_r`` and ``phase = φ_t + φ_r``.

    The update equation uses the **D form** (the standard charge-conserving discretization, see
    kernels/modulation.cu)::

        ε^{n+1} E^{n+1} = ε^n E^n + (dt/ε₀)·curl H
        ⇒ ca = ε(t_n)/ε(t_{n+1}),  cb = (dt/ε₀)/ε(t_{n+1})

    which folds straight into the ``ca/cb`` update_e already has: every step the kernel only
    rewrites the coefficients of the modulated cells, and in a scene without modulation it is
    never launched at all (a single path, the same approach as the absorber).

    Attributes:
        freq: Modulation frequency [Hz] (the validation set only contains cases with one global
            frequency; multiple frequencies fail closed).
        comp: ``(n_entry,)`` component, 0/1/2 -> Ex/Ey/Ez.
        cell: ``(n_entry,)`` flat cell index ``(i*ny + j)*nz + k``.
        amp: ``(n_entry,)`` modulation amplitude ``A_t·A_r`` (an absolute change of ε,
            dimensionless).
        phase: ``(n_entry,)`` phase ``φ_t + φ_r`` [rad].
    """

    freq: float
    comp: np.ndarray
    cell: np.ndarray
    amp: np.ndarray
    phase: np.ndarray

    @property
    def n_entry(self) -> int:
        return int(self.comp.size)

    def validate(self, shape: tuple[int, int, int]) -> None:
        """Shapes, index ranges, uniqueness of (comp, cell), finiteness of the amplitudes.

        A duplicate (comp, cell) lets the coefficient written later overwrite the earlier one,
        giving a run that completes but is physically wrong; better to stop here.
        """
        n = self.n_entry
        for name, arr in (("amp", self.amp), ("phase", self.phase)):
            if arr.shape != (n,):
                raise ValueError(f"TimeModulation.{name} has shape {arr.shape}, expected ({n},)")
            if not np.all(np.isfinite(arr)):
                raise ValueError(f"TimeModulation.{name} contains non-finite values")
        if not (self.freq > 0):
            raise ValueError(f"TimeModulation.freq = {self.freq}, must be positive")
        _check_csr("TimeModulation", shape, self.comp, self.cell, [])


@dataclass(frozen=True)
class TensorEps:
    """Cells of a fully anisotropic medium (``FullyAnisotropicMedium``).

    Each entry is "the 3×3 material tensor of component ``comp`` at its own Yee position". From
    it the solver computes row ``comp`` of ``M1 = A⁻¹B`` and of ``M2 = (Δt/ε₀)A⁻¹``
    (``A = ε + σΔt/(2ε₀)``, ``B = ε − σΔt/(2ε₀)``), averaging the components that do not sit at
    this position over 4 neighbors; for the discretization see ``kernels/tensor.cu``.

    The entry table has to be **closed**: the 4-neighbor average of a true tensor entry (one
    with off-diagonal coupling) reads the "stretched curl" of the neighboring cells, so those
    cells must be in the table too, passing their value through; otherwise what is read is the
    E^{n+1} of an ordinary cell. The builder (scene/media.py) guarantees this with a ring of
    **diagonal halo** entries, and :meth:`validate` enforces it.

    Attributes:
        comp: ``(n,)`` 0/1/2 -> Ex/Ey/Ez.
        cell: ``(n,)`` flat cell index, same convention as :class:`Dispersion`.
        eps: ``(n, 3, 3)`` relative permittivity tensor (symmetric). Halo entries are diagonal
            (the diagonal element takes the ε of that cell's scalar path, the remaining diagonal
            entries are set to 1).
        sigma: ``(n, 3, 3)`` conductivity [S/m]. **Symmetry is not required**: the antisymmetric
            part is exactly the gyrotropic (magneto-optic) term.
    """

    comp: np.ndarray
    cell: np.ndarray
    eps: np.ndarray
    sigma: np.ndarray

    @property
    def n_entry(self) -> int:
        return int(self.comp.size)

    def rows(self, dt: float) -> tuple[np.ndarray, np.ndarray]:
        """``(m1, m2)``, each ``(n, 3)`` float64: the two coefficient rows the kernel uses."""
        s = self.sigma * (dt / (2.0 * cpml.EPSILON_0))
        a_inv = np.linalg.inv(self.eps + s)
        m1_full = a_inv @ (self.eps - s)
        m2_full = (dt / cpml.EPSILON_0) * a_inv
        rows = np.arange(self.n_entry)
        return m1_full[rows, self.comp, :], m2_full[rows, self.comp, :]

    def coupled_mask(self) -> np.ndarray:
        """``(n,)`` bool: the "true tensor" entries that have off-diagonal coupling (halo
        entries are False)."""
        off = ~np.eye(3, dtype=bool)
        return (np.abs(self.eps[:, off]).max(axis=1) > 0.0) \
            | (np.abs(self.sigma[:, off]).max(axis=1) > 0.0)

    def validate(self, shape: tuple[int, int, int],
                 tables: list[dict] | None = None) -> None:
        """Shapes, uniqueness, and (when the neighbor tables are given) closure of the
        4-neighbor stencil.

        With closure broken the kernel takes a neighboring cell's E^{n+1} for the "stretched
        curl": the run completes but is physically wrong, so better to stop here.
        """
        n = self.n_entry
        for name, arr in (("eps", self.eps), ("sigma", self.sigma)):
            if arr.shape != (n, 3, 3):
                raise ValueError(f"TensorEps.{name} has shape {arr.shape}, expected ({n}, 3, 3)")
        key = _check_csr("TensorEps", shape, self.comp, self.cell, [])
        if n == 0 or tables is None:
            return
        ncell = int(np.prod(shape))
        have = set(key.tolist())
        coords = unravel_cells(self.cell, shape)
        cpl = self.coupled_mask()
        missing = 0
        for a in range(3):
            for b in range(3):
                if a == b:
                    continue
                m = cpl & (self.comp == a) & (
                    (np.abs(self.eps[:, a, b]) > 0) | (np.abs(self.sigma[:, a, b]) > 0))
                if not m.any():
                    continue
                c = [coords[0][m], coords[1][m], coords[2][m]]
                # Axis a takes {0, nxt}, axis b takes {prv, 0}: same stencil as kernels/tensor.cu
                for da in (False, True):
                    for db in (False, True):
                        q = [c[0].copy(), c[1].copy(), c[2].copy()]
                        if da:
                            q[a] = tables[a]["nxt"][c[a]]
                        if db:
                            q[b] = tables[b]["prv"][c[b]]
                        flat = ravel_cells(q[0].astype(np.int64), q[1], q[2], shape)
                        want = np.int64(b) * ncell + flat
                        missing += sum(1 for w in want.tolist() if w not in have)
        if missing:
            raise ValueError(
                f"TensorEps closure is broken: {missing} of the 4-neighbor stencil points are "
                "not in the entry table (halo not wide enough); the kernel would take the E of "
                "an ordinary cell for the stretched curl")


@dataclass(frozen=True)
class AbsorberSlab:
    """One ``Absorber`` face: a **matched** lossy medium graded layer by layer (an adiabatic
    absorber).

    Tidy3D's Absorber is a "multilayer system with graded conductivity terminated by a PEC" (its
    own docstring), not a coordinate stretch; OpenEM's reference implementation carries a
    matched magnetic loss (``matkh.bin``). Our implementation: every step, E and H inside the
    layers are each multiplied once by the same ``exp(−Γ·dt)``.

    Attributes:
        axis: Absorbing axis.
        g0: Global cell index where this slab starts along the absorbing axis (0 at the lo end,
            n−nlay at the hi end).
        decay_int: ``(nlay,)`` per-step decay factor at the whole-cell positions.
        decay_half: ``(nlay,)`` the same at the half-cell positions.
    """

    axis: int
    g0: int
    decay_int: np.ndarray
    decay_half: np.ndarray

    def validate(self) -> None:
        for nm in ("decay_int", "decay_half"):
            d = getattr(self, nm)
            if d.ndim != 1 or not np.all((d > 0.0) & (d <= 1.0)):
                raise ValueError(
                    f"AbsorberSlab.{nm} must be (nlay,) with every factor in (0, 1]; got "
                    f"shape={d.shape}, min={d.min():.3e}, max={d.max():.3e}")


@dataclass(frozen=True)
class ModeMonitorSpec:
    """One ``td.ModeMonitor``: the solver only computes the fields on the plane, the
    decomposition happens in post-processing.

    Attributes:
        name: Name of the original monitor, used to line results up with the reference output.
        plane_name: Name of the ``FieldMonitor`` it expands into, which the phasors come from.
    """

    name: str
    plane_name: str


@dataclass(frozen=True)
class BoxFluxMonitor:
    """Total flux through a closed box: split into six face :class:`FluxMonitor` objects and
    summed along the outward normals.

    Tidy3D allows a ``FluxMonitor`` to be three-dimensional, meaning the net power through a
    closed surface (outward is positive). The solver only computes planar flux, so the box is
    split up, computed face by face and summed; each face carries its sign in ``normal_dir``
    (the outward normal of a ``-`` face points in the negative direction).

    Attributes:
        name: Monitor name in Tidy3D.
        face_names: Names of the six face monitors.
    """

    name: str
    face_names: tuple[str, ...]


@dataclass(frozen=True)
class ProjectionMonitor:
    """Far-field projection monitor: the solver does not compute it directly, it only records
    which faces it rests on.

    The projection box is expanded into six face :class:`FieldMonitor` objects
    (``scene/projection.py``) and the solver computes the phasors on those faces as usual; the
    near-to-far transform is handed to the Tidy3D client's ``FieldProjector`` in
    post-processing. **Anything the client can give us, we do not write ourselves**, the same
    reasoning as for subpixel averaging.

    Attributes:
        name: Monitor name in Tidy3D; the same name must be used when writing results back.
        surface_names: Names of the six face monitors, in the order matching ``normal_dirs``.
        normal_dirs: The **outward normal** of each face, ``"+"`` or ``"-"``.
    """

    name: str
    surface_names: tuple[str, ...]
    normal_dirs: tuple[str, ...]


@dataclass(frozen=True)
class Scene:
    """Everything the solver takes as input. Pure numpy, no tidy3d objects."""

    grid: Grid
    dt: float
    num_time_steps: int
    shutoff: float

    #: Relative permittivity at the Yee component positions, each ``(nx, ny, nz)``
    eps_ex: np.ndarray
    eps_ey: np.ndarray
    eps_ez: np.ndarray

    #: Coefficients of each PML face; the key is ``(axis, "lo"|"hi")``
    pml: dict[tuple[int, str], cpml.PMLCoeffs]

    #: Conductivity [S/m], at the same points and with the same shape as ``eps_*``. ``None``
    #: means lossless. Either all three are None or none of them are; :attr:`any_loss` decides.
    sigma_ex: np.ndarray | None = None
    sigma_ey: np.ndarray | None = None
    sigma_ez: np.ndarray | None = None

    #: Flat indices ``(i*ny + j)*nz + k`` of the PEC cells, per component; E is identically 0 in
    #: these cells (ca=cb=0 in the solver, and every source injection is multiplied by cb, so it
    #: is zeroed right along with it). Staircasing convention; the mask comes from exactly one
    #: place, the ε verdict in ``scene/media.py`` (pec_val). ``None`` means the scene has no
    #: PEC; all three are None together or all three are arrays together (possibly empty).
    pec_ex: np.ndarray | None = None
    pec_ey: np.ndarray | None = None
    pec_ez: np.ndarray | None = None

    #: Dispersion pole data (E form). ``None`` means nothing in the domain is dispersive.
    dispersion: Dispersion | None = None

    #: Pole data for cells on a slanted interface (mixed form). ``None`` means there are none.
    dispersion_mix: DispersionMix | None = None

    #: Time-varying medium (modulation_spec). ``None`` means no modulation, and then the solver
    #: never launches the modulation kernel even once, leaving results bitwise identical.
    modulation: TimeModulation | None = None

    #: Cells of a fully anisotropic medium (FullyAnisotropicMedium). ``None`` means there are
    #: none.
    tensor: TensorEps | None = None

    #: Tidy3D's ``symmetry``, per axis -1 (PEC plane) / 0 (none) / +1 (PMC plane).
    #: The solver itself never reads it; it only serves to reproduce the cell-count convention
    #: of Tidy3D's solver. For folding onto a half domain with it see fold.py (used by the speed
    #: benchmark script).
    symmetry: tuple[int, int, int] = (0, 0, 0)

    #: Bloch wave vector per axis (Tidy3D's ``bloch_vec`` convention, in units of 2π/L; a
    #: periodic boundary is 0). If any component is nonzero the solver takes the complex-field
    #: path (kernels/bloch.cu); for the subset of physics that supports it see
    #: solver._refuse_unsupported_bloch.
    bloch_k: tuple[float, float, float] = (0.0, 0.0, 0.0)

    sources: list[PlaneWaveSource] = field(default_factory=list)
    tfsf_sources: list[TFSFSource] = field(default_factory=list)
    dipoles: list[PointDipole] = field(default_factory=list)
    mode_sources: list[ModeSource] = field(default_factory=list)
    flux_monitors: list[FluxMonitor] = field(default_factory=list)
    field_monitors: list[FieldMonitor] = field(default_factory=list)
    field_time_monitors: list[FieldTimeMonitor] = field(default_factory=list)
    flux_time_monitors: list[FluxTimeMonitor] = field(default_factory=list)
    permittivity_monitors: list[PermittivityMonitor] = field(default_factory=list)
    projection_monitors: list[ProjectionMonitor] = field(default_factory=list)
    box_flux_monitors: list[BoxFluxMonitor] = field(default_factory=list)
    mode_monitors: list[ModeMonitorSpec] = field(default_factory=list)
    absorbers: list[AbsorberSlab] = field(default_factory=list)

    @property
    def shape(self) -> tuple[int, int, int]:
        return self.grid.shape

    def eps_box(self, mon: PermittivityMonitor) -> tuple[np.ndarray, ...]:
        """``(eps_xx, eps_yy, eps_zz)`` sliced to the monitor box, at their Yee positions.

        For a lossy medium only the real part of the complex ε is given here; the imaginary part
        would have to be rebuilt from σ and the frequency, and every case in the validation set
        that carries a PermittivityMonitor is lossless.
        """
        i0, j0, k0 = mon.origin
        ni, nj, nk = mon.box
        sl = (slice(i0, i0 + ni), slice(j0, j0 + nj), slice(k0, k0 + nk))
        return tuple(e[sl] for e in (self.eps_ex, self.eps_ey, self.eps_ez))

    @property
    def any_dispersion_mix(self) -> bool:
        """Whether any cell takes the mixed form."""
        return self.dispersion_mix is not None and self.dispersion_mix.n_entry > 0

    @property
    def any_modulation(self) -> bool:
        """Whether any cell is modulated."""
        return self.modulation is not None and self.modulation.n_entry > 0

    @property
    def any_dispersion(self) -> bool:
        """Whether there is any dispersive material. Decides whether P and hist are
        allocated."""
        return self.dispersion is not None and self.dispersion.n_entry > 0

    @property
    def any_tensor(self) -> bool:
        """Whether there are fully anisotropic cells. Decides whether the three tensor kernels
        are launched."""
        return self.tensor is not None and self.tensor.n_entry > 0

    @property
    def any_bloch(self) -> bool:
        """Whether the Bloch wave vector is nonzero. Decides between the real and complex field
        path."""
        return any(k != 0.0 for k in self.bloch_k)

    @property
    def any_loss(self) -> bool:
        """Whether there is any conductivity. Decides which kernel the solver uses."""
        return self.sigma_ex is not None

    @property
    def any_pec(self) -> bool:
        """Whether there are any PEC cells."""
        return self.pec_ex is not None

    @property
    def num_cells_cloud(self) -> int:
        """The convention behind "computational grid points" in Tidy3D's log: on a symmetric
        axis its solver computes only half the domain and adds one ghost layer at each end,
        i.e. ``∏(N/2 + 2)`` (matched exactly on 5 cases).

        Only used to report the same number Tidy3D's solver does. **We run the full domain**;
        for what is actually allocated see ``grid.shape_with_ghost``.
        """
        total = 1
        for axis, sym in zip(self.grid.axes, self.symmetry):
            if sym == 0:
                total *= axis.n_with_ghost
            else:
                # Symmetric axis: half the domain plus one ghost layer at each end
                total *= axis.n // 2 + 2
        return total


def homogeneous(sc: Scene) -> Scene:
    """Drop every structure and fill the whole domain with the medium at the source plane.

    Composed with :func:`lossless`, this is the reference run that measures the incident power.
    The PML coefficients are left unchanged.
    """
    if sc.any_tensor:
        raise NotImplementedError(
            "homogeneous() does not support tensor-medium scenes yet: the medium at the source "
            "plane is a 3x3 tensor, and the reference run convention for filling the domain "
            "with it has never been defined")
    src = sc.sources[0]
    at = (0, 0, src.plane_index)
    const = np.full(sc.eps_ex.shape, float(sc.eps_ex[at]), dtype=np.float64)
    out = replace(sc, eps_ex=const, eps_ey=const.copy(), eps_ez=const.copy(),
                  pec_ex=None, pec_ey=None, pec_ez=None,
                  dispersion=_homogeneous_dispersion(sc, at))
    if sc.sigma_ex is not None:
        sig = np.full(sc.eps_ex.shape, float(sc.sigma_ex[at]), dtype=np.float64)
        out = replace(out, sigma_ex=sig, sigma_ey=sig.copy(), sigma_ez=sig.copy())
    return out


def lossless(sc: Scene) -> Scene:
    """Zero the conductivity **and recompute the TF/SF incident tables**.

    Measuring the incident power requires a lossless reference: a uniform lossless medium
    neither reflects nor absorbs, so the flux at the monitor equals the injected power.
    BeerLambert's loss lives in the **background** medium (it has no structures at all), so
    :func:`homogeneous` alone cannot produce a reference run; that way the reference and the
    structured run would be identical.

    **The incident tables must be recomputed along with it.** ``hy_inc`` contains the complex
    wave impedance and the half-cell gain of the medium at the source plane. Clearing σ without
    recomputing leaves the reference run injecting an H matched to a lossy medium into a
    lossless one, with the emitted amplitude off by a factor of ``gain × |n_c|/n``. On
    BeerLambert this factor was measured at **1.00353**, and all of it lands in
    ``T = flux_lossy/flux_lossless``: T comes out 0.64% low, disguised as half a cell of extra
    optical path.
    """
    if sc.any_tensor and np.any(np.abs(sc.tensor.sigma) > 0):
        raise NotImplementedError(
            "lossless() does not support scenes with a σ tensor yet: clearing only the scalar σ "
            "would silently leave the tensor loss behind")
    out = replace(sc, sigma_ex=None, sigma_ey=None, sigma_ez=None)
    if sc.sigma_ex is None or not sc.sources:
        return out                       # already lossless, or no plane wave source
    from openem.scene import sources as _sources

    new = []
    for src in sc.sources:
        if src.incident is None:
            raise ValueError(
                "a plane wave source in a lossy scene carries no incident parameters, so "
                "lossless() cannot recompute the incident tables. Clearing only σ makes the "
                "emitted amplitude of the reference run wrong by a complex wave impedance "
                "factor, and it is **invisible**: the flux is still a perfectly normal number."
            )
        kw = dict(src.incident)
        kw["sigma"] = 0.0
        ex, hy = _sources._incident_tables(**kw)
        new.append(replace(src, ex_inc=ex, hy_inc=hy))
    return replace(out, sources=new)


def _homogeneous_dispersion(sc: Scene, at: tuple[int, int, int]) -> Dispersion | None:
    """Copy the poles of the cell at the source plane over the whole domain, for
    :func:`homogeneous`.

    What the reference run needs is an infinite uniform medium identical to the medium at the
    source. Simply dropping the dispersion here would change the material of the reference run
    whenever the medium at the source is itself dispersive: the incident power would come out
    wrong, and **the error would be invisible** (the flux is still a perfectly normal number).
    So it is carried over as is, and if it cannot be, this raises.
    """
    d = sc.dispersion
    if d is None or d.n_entry == 0:
        return None
    nx, ny, nz = sc.shape
    flat = ravel_cells(at[0], at[1], at[2], sc.shape)
    ncell = nx * ny * nz
    src = {int(d.comp[e]): int(e) for e in np.flatnonzero(d.cell == flat)}
    if not src:
        return None                      # the medium at the source is not dispersive: common
    if len(src) != 3:
        raise NotImplementedError(
            f"only {len(src)}/3 components of the cell at the source plane are dispersive, so "
            "no uniform reference run can be built")
    counts = {c: int(d.pole_ofs[e + 1] - d.pole_ofs[e]) for c, e in src.items()}
    comp = np.repeat(np.arange(3, dtype=np.int32), ncell)
    cell = np.tile(np.arange(ncell, dtype=np.int32), 3)
    per = np.array([counts[c] for c in range(3)], dtype=np.int64).repeat(ncell)
    ofs = np.concatenate([[0], np.cumsum(per)]).astype(np.int32)
    a_parts, b_parts = [], []
    for c in range(3):
        e = src[c]
        sl = slice(int(d.pole_ofs[e]), int(d.pole_ofs[e + 1]))
        a_parts.append(np.tile(d.am1[sl], ncell))
        b_parts.append(np.tile(d.b[sl], ncell))
    g = np.concatenate([np.full(ncell, d.g[src[c]]) for c in range(3)])
    return Dispersion(comp=comp, cell=cell, pole_ofs=ofs,
                      am1=np.concatenate(a_parts), b=np.concatenate(b_parts), g=g)
