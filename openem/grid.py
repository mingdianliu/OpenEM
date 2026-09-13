# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Mingdian Liu
"""Non-uniform Yee grid.

We do not generate the grid: the cell boundary coordinates are taken straight from
``sim.grid.boundaries`` (see docs/architecture.md). This module only turns them into the spacings
the update equations need, and counts the ghost layers.

Yee convention (the same as Tidy3D, FDTDX ``core/physics/curl.py:29-32``):

- updating **H** divides by the **primal** spacing ``dl[i] = edges[i+1] - edges[i]``
- updating **E** divides by the **dual** spacing ``(dl[i] + dl[i-1]) / 2``

Ghost rule: add one layer at every end where a neighbour **cannot be reached**.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

#: How many ghost layers each boundary type needs at that end (see docs/architecture.md).
#: A periodic end can reach its neighbour by wrapping round, so it needs none.
GHOST_BY_BOUNDARY: dict[str, int] = {
    "PML": 1,
    "StablePML": 1,
    "Absorber": 1,
    "PECBoundary": 1,
    "PMCBoundary": 1,
    "Periodic": 0,
    "BlochBoundary": 0,
    # Symmetry plane (the low end after fold.py has folded): the quantities on the plane are
    # inside the domain and the neighbour outside it goes through the mirror mask, so one ghost
    # layer is needed, as for PEC.
    "SymmetryPEC": 1,
    "SymmetryPMC": 1,
}

#: Boundary types that absorb or truncate (the ghost values at those ends follow the
#: absorbing-end rule, with tangential E set to 0)
ABSORBING = ("PML", "StablePML", "Absorber")

#: Whether each of the six field components sits on the primal edge (True) or at the cell centre
#: (False) along each axis.
#: Taken from the measured ``sim.grid.yee``: Ex(ctr_x, edge_y, edge_z), Hx(edge_x, ctr_y, ctr_z)
#: and so on (the layout table at the top of kernels/yee.cu). This is the single source and
#: everywhere else derives from it: colocate.AT_CENTER is its negation; tfsf_oblique looks it up
#: by component number through :func:`yee_on_edge`; fold's mirror index mapping uses it directly.
YEE_ON_EDGE: dict[str, tuple[bool, bool, bool]] = {
    "Ex": (False, True, True),
    "Ey": (True, False, True),
    "Ez": (True, True, False),
    "Hx": (True, False, False),
    "Hy": (False, True, False),
    "Hz": (False, False, True),
}

#: Mirror eigenvalue coefficients of the six field components along axis ax: σ = sign * s (the
#: derivation is at the top of fold.py).
#: Tangential E = +s, normal E = -s; tangential H = -s, normal H = +s.
#: The key is the component name, and the value gives, per axis, the sign before multiplying by
#: s: MIRROR_SIGN[comp][ax].
MIRROR_SIGN: dict[str, tuple[int, int, int]] = {
    "Ex": (-1, +1, +1),
    "Ey": (+1, -1, +1),
    "Ez": (+1, +1, -1),
    "Hx": (+1, -1, -1),
    "Hy": (-1, +1, -1),
    "Hz": (-1, -1, +1),
}


def yee_on_edge(comp: int, magnetic: bool) -> tuple[bool, bool, bool]:
    """Look up :data:`YEE_ON_EDGE` by component number (0/1/2 = x/y/z)."""
    return YEE_ON_EDGE[("H" if magnetic else "E") + "xyz"[comp]]


def transverse_axes(axis: int) -> tuple[int, int]:
    """The two axes other than the normal ``axis``, in **ascending** order (the (t1, t2)
    convention of monitors and flux)."""
    return tuple(i for i in range(3) if i != axis)  # type: ignore[return-value]


def cyclic_axes(axis: int) -> tuple[int, int]:
    """The two axes other than the normal ``axis``, in **cyclic right-handed order**
    ``((axis+1)%3, (axis+2)%3)`` (the (u, v) convention of mode sources and TFSF; the Yee curl
    keeps its form under a cyclic permutation)."""
    return (axis + 1) % 3, (axis + 2) % 3


def unravel_cells(cell, shape):
    """Flat cell index ``(i*ny + j)*nz + k`` -> ``(ii, jj, kk)``. ``shape`` uses the last two
    dimensions ``(ny, nz)`` (the pitched version passes ``(..., ny, nzp)``). Integer arithmetic,
    written identically for numpy and cupy, and scalars work too."""
    ny, nz = int(shape[-2]), int(shape[-1])
    return cell // (ny * nz), (cell // nz) % ny, cell % nz


def ravel_cells(ii, jj, kk, shape):
    """The inverse of :func:`unravel_cells`: ``(i*ny + j)*nz + k``."""
    ny, nz = int(shape[-2]), int(shape[-1])
    return (ii * ny + jj) * nz + kk


def outer3(a, b, c):
    """Outer product of three 1D tables, ``a[i]·b[j]·c[k]``, with the multiplication order fixed
    at (a·b)·c: the three Yee dual volumes and the solver's GPU version all go through this one
    expression, written identically for numpy and cupy."""
    return a[:, None, None] * b[None, :, None] * c[None, None, :]


def ghost_count(boundary_type: str, has_symmetry: bool = False) -> int:
    """How many ghost layers this end needs.

    Args:
        boundary_type: Tidy3D's boundary type name, e.g. ``"PML"`` / ``"Periodic"``.
        has_symmetry: a periodic axis that also carries a symmetry plane needs one layer as well.

    Raises:
        KeyError: unknown boundary type. Better to raise than to treat it silently as 0: one
            layer too few would put the cell count noticeably out of step with the reference
            implementation, and nothing would warn you.
    """
    if boundary_type not in GHOST_BY_BOUNDARY:
        raise KeyError(
            f"unknown boundary type {boundary_type!r}; the known ones are "
            f"{sorted(GHOST_BY_BOUNDARY)}. The number of ghost layers has to be defined "
            "explicitly, never guessed."
        )
    n = GHOST_BY_BOUNDARY[boundary_type]
    if has_symmetry and n == 0:
        # one layer when the same end is both periodic and a symmetry plane
        return 1
    return n


@dataclass(frozen=True)
class Axis:
    """The grid of one axis.

    Attributes:
        edges: ``(n+1,)`` cell boundary coordinates in **metres**, including the PML region.
            (Tidy3D's unit is µm; the conversion happens in scene/.)
        boundary_lo: boundary type at the low end.
        boundary_hi: boundary type at the high end.
    """

    edges: np.ndarray
    boundary_lo: str
    boundary_hi: str
    #: Whether this is a 2D flat axis in Tidy3D's sense (``sim.size==0``). ``None`` = unknown,
    #: falling back to the old rule "one cell means flat" (axes rebuilt by fold/serialize land
    #: here). An axis with size!=0 that happens to have exactly one cell (y=0.01 µm in
    #: GratingEfficiency) is **genuinely 3D**, and its area has to use the actual dl.
    flat: bool | None = None

    def __post_init__(self) -> None:
        if self.edges.ndim != 1 or self.edges.size < 2:
            raise ValueError(
                f"edges must be a 1D array of length >= 2, got shape={self.edges.shape}")
        if not np.all(np.diff(self.edges) > 0):
            raise ValueError("edges must be strictly increasing")

    @property
    def n(self) -> int:
        """Number of physical cells, ghost layers excluded."""
        return self.edges.size - 1

    @property
    def is_flat(self) -> bool:
        """2D flat axis: the explicit flag wins, and without it fall back to "just one cell"."""
        return bool(self.flat) if self.flat is not None else self.n == 1

    @property
    def dl(self) -> np.ndarray:
        """Primal spacing ``(n,)``, used when updating H."""
        return np.diff(self.edges)

    @property
    def ghost_lo(self) -> int:
        return ghost_count(self.boundary_lo)

    @property
    def ghost_hi(self) -> int:
        return ghost_count(self.boundary_hi)

    @property
    def n_with_ghost(self) -> int:
        """Array length including ghost layers, which is the size actually allocated."""
        return self.n + self.ghost_lo + self.ghost_hi

    @property
    def is_periodic(self) -> bool:
        return self.boundary_lo in ("Periodic", "BlochBoundary")

    @property
    def dl_dual(self) -> np.ndarray:
        """Dual spacing ``(n,)``, used when updating E.

        ``dl_dual[i] = (dl[i] + dl[i-1]) / 2``. At i=0 this needs a ``dl[-1]`` from outside the
        domain: a periodic axis wraps round and takes ``dl[n-1]``, while an absorbing end takes
        ``dl[0]``, as if the grid extended one cell outwards.
        """
        dl = self.dl
        prev = np.empty_like(dl)
        prev[1:] = dl[:-1]
        # Symmetry plane at the low end: the mirrored cell spacing is dl[0] (fold.py has already
        # verified that the grid is symmetric about the plane), the same value as the absorbing
        # end's "extend one cell outwards", so no branch is needed here.
        prev[0] = dl[-1] if self.is_periodic else dl[0]
        return 0.5 * (dl + prev)

    @property
    def centers(self) -> np.ndarray:
        """Cell centre coordinates ``(n,)``."""
        return 0.5 * (self.edges[:-1] + self.edges[1:])

    def index_tables(self, bloch_k: float = 0.0) -> dict[str, np.ndarray]:
        """Neighbour indices and masks: the kernels express the boundary topology through these
        tables instead of through branches.

        - ``nxt`` / ``prv``: neighbour indices. A periodic axis wraps round; an absorbing axis is
          clamped inside the domain (to be used together with the masks).
        - ``mnx`` / ``mpv``: whether that neighbour really exists. An out-of-domain neighbour on
          an absorbing axis is multiplied by 0, so both "tangential E outside the wall = 0" and
          "the backward difference cannot read H" become multiplications rather than branches.
        - ``pec``: 0 at ``i=0`` on an absorbing axis, because the PML is backed by PEC and the
          tangential E there is identically zero.

        On a periodic axis all three masks are 1, so one and the same kernel code is correct for
        both topologies.

        Args:
            bloch_k: Bloch wave vector (Tidy3D's ``bloch_vec`` convention, in units of 2π/L).
                When it is non-zero only a periodic axis is allowed and ``mnx``/``mpv`` turn
                **complex128**: the Bloch condition is ``u(x+L) = u(x)·exp(+i·2π·k)``, so the
                forward difference wrapping round at ``i=n-1`` picks up ``exp(+i·2π·k)`` and the
                backward difference wrapping round at ``i=0`` picks up ``exp(-i·2π·k)``. The
                wrap-around phase appears **only** in those two mask elements, so the complex
                kernel (kernels/bloch.cu) is still branch-free.
                With ``bloch_k=0`` the real tables returned are bitwise identical to before.
        """
        n = self.n
        i = np.arange(n)
        if self.is_periodic:
            mnx = np.ones(n)
            mpv = np.ones(n)
            if bloch_k != 0.0:
                mnx = mnx.astype(np.complex128)
                mpv = mpv.astype(np.complex128)
                mnx[n - 1] = np.exp(+2j * np.pi * bloch_k)
                mpv[0] = np.exp(-2j * np.pi * bloch_k)
            return {
                "nxt": ((i + 1) % n).astype(np.int32),
                "prv": ((i - 1) % n).astype(np.int32),
                "mnx": mnx,
                "mpv": mpv,
                "pec": np.ones(n),
            }
        if bloch_k != 0.0:
            raise ValueError(
                f"bloch_k={bloch_k:g} != 0 only makes sense on a periodic axis, but this axis "
                f"is {self.boundary_lo}/{self.boundary_hi}")
        pec = np.ones(n)
        mpv = (i > 0).astype(np.float64)
        if self.boundary_lo == "PMCBoundary":
            mpv[0] = -1.0        # PMC low wall: tangential H odd (H(-½)=-H(+½) -> H=0 on wall)
        elif self.boundary_lo in ("SymmetryPEC", "SymmetryPMC"):
            # Symmetry plane at the low end (fold.py). The backward difference of a tangential
            # quantity on the plane (i=0) has to read the mirror outside it: E on the edge plane
            # reads H at the centre, and the mirrored cell is ctr[0] itself, so the mask carries
            # the mirror eigenvalue. Tangential H is -s under either symmetry (s=+1 PMC -> -1;
            # s=-1 PEC -> +1, but then tangential E is zeroed by pec and the mask value does not
            # matter). The derivation is at the top of fold.py.
            s = -1.0 if self.boundary_lo == "SymmetryPEC" else +1.0
            mpv[0] = -s
            if s < 0:
                pec[0] = 0.0        # PEC symmetry plane: tangential E on it is always zero
        else:
            pec[0] = 0.0            # absorbing end: the PML is backed by PEC
        nxt = np.minimum(i + 1, n - 1).astype(np.int32)
        mnx = (i + 1 < n).astype(np.float64)
        if self.boundary_hi == "PMCBoundary":
            nxt[n - 1] = n - 1
            mnx[n - 1] = 1.0     # PMC high wall: tangential E even (∂E_t/∂n=0)
        elif self.boundary_hi in ("SymmetryPEC", "SymmetryPMC"):
            # High end of the double mirror wall (fold.py). The last row is the wrap-around cell
            # (identical to full-domain row 0), and its forward neighbour is full-domain row 1 =
            # the mirror about the high wall = local row n-2; the mask carries the mirror
            # eigenvalue +s of the component being read (tangential E). With s=-1 the tangential
            # E on the wall is always zero and pec is zeroed, just as at the low end.
            s_hi = -1.0 if self.boundary_hi == "SymmetryPEC" else +1.0
            nxt[n - 1] = n - 2
            mnx[n - 1] = s_hi
            if s_hi < 0:
                pec[n - 1] = 0.0
        return {
            "nxt": nxt,
            "prv": np.maximum(i - 1, 0).astype(np.int32),
            "mnx": mnx,
            "mpv": mpv,
            "pec": pec,
        }


@dataclass(frozen=True)
class Grid:
    """The grid of all three axes."""

    x: Axis
    y: Axis
    z: Axis

    @property
    def axes(self) -> tuple[Axis, Axis, Axis]:
        return (self.x, self.y, self.z)

    @property
    def shape(self) -> tuple[int, int, int]:
        """Physical cell counts ``(nx, ny, nz)``, matching Tidy3D's ``grid.num_cells``."""
        return (self.x.n, self.y.n, self.z.n)

    @property
    def shape_with_ghost(self) -> tuple[int, int, int]:
        """Array shape including ghost layers, which is the size actually allocated."""
        return tuple(a.n_with_ghost for a in self.axes)  # type: ignore[return-value]

    @property
    def num_cells(self) -> int:
        """Total number of physical cells, matching ``prod(sim.grid.num_cells)``."""
        nx, ny, nz = self.shape
        return nx * ny * nz

    @property
    def num_cells_cloud(self) -> int:
        """Cell count including ghost layers, i.e. the "computational grid points" figure in the
        reference implementation's log.

        On a symmetry axis the reference implementation counts only the half domain; that
        convention is in :attr:`model.Scene.num_cells_cloud`.
        """
        nx, ny, nz = self.shape_with_ghost
        return nx * ny * nz

    def cell_volumes(self) -> np.ndarray:
        """Primal cell volumes ``(nx, ny, nz)``. (Currently only tests/test_grid.py uses it.)"""
        return outer3(self.x.dl, self.y.dl, self.z.dl)

    def yee_dual_volumes(self) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """The dual volume ``(nx, ny, nz)`` of each of ``Ex``/``Ey``/``Ez``.

        The three E components sit at different places on the Yee grid, so their dual cell
        volumes differ too:

        - ``Ex`` at (ctr_x, edge_y, edge_z) -> ``dx_primal · dy_dual · dz_dual``
        - ``Ey`` at (edge_x, ctr_y, edge_z) -> ``dx_dual · dy_primal · dz_dual``
        - ``Ez`` at (edge_x, edge_y, ctr_z) -> ``dx_dual · dy_dual · dz_primal``

        The volume weighting of shutoff has to use this and not one shared primal volume; the
        basis is ``grid_weighting=component-yee-dual`` in the reference solver.
        On a uniform grid all three are equal, so BiosensorGrating cannot tell them apart.
        """
        xp, yp, zp = (a.dl for a in self.axes)
        xd, yd, zd = (a.dl_dual for a in self.axes)
        return outer3(xp, yd, zd), outer3(xd, yp, zd), outer3(xd, yd, zp)

    def yee_h_dual_volumes(self) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """The dual volume of each of ``Hx``/``Hy``/``Hz``, dual to :meth:`yee_dual_volumes`.

        The Yee positions of the H components swap primal and dual exactly (the layout table at
        the top of kernels/yee.cu):

        - ``Hx`` at (edge_x, ctr_y, ctr_z) -> ``dx_dual · dy_primal · dz_primal``
        - ``Hy`` at (ctr_x, edge_y, ctr_z) -> ``dx_primal · dy_dual · dz_primal``
        - ``Hz`` at (ctr_x, ctr_y, edge_z) -> ``dx_primal · dy_primal · dz_dual``

        This is what a magnetic dipole divides by when spreading a magnetic moment into a
        magnetic current density (scene/sources.py).
        """
        xp, yp, zp = (a.dl for a in self.axes)
        xd, yd, zd = (a.dl_dual for a in self.axes)
        return outer3(xd, yp, zp), outer3(xp, yd, zp), outer3(xp, yp, zd)


def normal_axis(size: tuple[float, float, float]) -> int:
    """The one axis whose ``size`` is 0 is the plane normal.

    Both sources and monitors use it to recognise the zero-thickness dimension.
    """
    zeros = [i for i, s in enumerate(size) if s == 0]
    if len(zeros) != 1:
        raise ValueError(f"expected exactly one zero dimension (a plane), got size={size}")
    return zeros[0]


def nearest_edge_index(edges: np.ndarray, coord: float) -> int:
    """Snap a plane coordinate to the nearest primal boundary index.

    What sits on a source plane or a FluxMonitor plane is tangential E, exactly at ``edges[k]``,
    so the snap has to go to a **boundary**. The "which cell does it fall in" criterion of
    :func:`nearest_index` will not do: in floating point ``edges[k]`` can be a hair larger than
    the target, and that would be off by a whole cell.

    The "half cell" has to be the **local spacing on either side of the snap point**, not the
    smallest spacing on the whole axis: an Auto grid can go down to ~1 nm inside a metal layer
    while the source plane lands in a coarse ~9.5 nm region, so taking the global minimum as the
    threshold would reject a legitimate snap (0.14 local cells off) as "not on a boundary" (both
    the plane-wave source at z=-0.7 µm and the T flux plane at z=+0.7 µm of MIMResonator ran into
    this). Tidy3D's own semantics is to snap a plane to the nearest grid position, and being
    within half a **local** cell is exactly what "nearest" means.

    Raises:
        ValueError: the offset exceeds half a local cell, which means the plane is nowhere near a
            Yee boundary.
    """
    i = int(np.argmin(np.abs(edges - coord)))
    resid = abs(float(edges[i]) - coord)
    dl = np.diff(edges)
    local = max(float(dl[i]) if i < dl.size else 0.0,
                float(dl[i - 1]) if i > 0 else 0.0)
    half = 0.5 * local
    if resid > half:
        raise ValueError(
            f"coordinate {coord!r} is off the nearest cell boundary {edges[i]!r} by {resid!r}, "
            f"more than the local half cell {half!r}"
        )
    return i


def nearest_index(edges: np.ndarray, coord: float) -> int:
    """Which cell a physical coordinate falls in: returns the ``i`` satisfying
    ``edges[i] <= coord < edges[i+1]``.

    A coordinate outside the domain raises instead of being clamped: silent clamping would make a
    monitor record on the wrong plane.
    (Currently only tests/test_grid.py uses it.)
    """
    if coord < edges[0] or coord > edges[-1]:
        raise ValueError(
            f"coordinate {coord!r} falls outside the grid [{edges[0]}, {edges[-1]}]")
    i = int(np.searchsorted(edges, coord, side="right") - 1)
    return min(max(i, 0), edges.size - 2)
