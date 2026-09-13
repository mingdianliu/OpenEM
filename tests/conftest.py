# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Mingdian Liu
"""Some extraction-side tests read an external archive at
``$OPENEM_TABLEA/<example>/Tidy3D/simulation.json``. This turns a missing archive into a skip
instead of a FileNotFoundError.

Without the archive those tests simply cannot run, and an explicit skip carrying the path in its
reason reflects the real state far better than forty red crosses. Patched at module level rather
than through a fixture, so module-scoped fixtures go through it too.
"""
import os
import pytest

_TABLEA = os.environ.get("OPENEM_TABLEA", "/nonexistent")  # external reference archive; without
#                                                           OPENEM_TABLEA these cases skip

try:
    import tidy3d as td

    _orig_from_file = td.Simulation.from_file.__func__

    def _patched_from_file(cls, fname, *a, **k):
        try:
            return _orig_from_file(cls, fname, *a, **k)
        except FileNotFoundError:
            if str(fname).startswith(_TABLEA):
                pytest.skip(f"external archive missing: {fname} (point OPENEM_TABLEA at the "
                            "archive directory to run this)")
            raise

    td.Simulation.from_file = classmethod(_patched_from_file)

    _orig_sd_from_file = td.SimulationData.from_file.__func__

    def _patched_sd_from_file(cls, fname, *a, **k):
        try:
            return _orig_sd_from_file(cls, fname, *a, **k)
        except (FileNotFoundError, OSError) as e:
            if str(fname).startswith(_TABLEA) and ("No such file" in str(e) or "Unable to" in str(e)):
                pytest.skip(f"archive missing: {fname}")
            raise

    td.SimulationData.from_file = classmethod(_patched_sd_from_file)
except Exception:  # do nothing when tidy3d is absent
    pass
