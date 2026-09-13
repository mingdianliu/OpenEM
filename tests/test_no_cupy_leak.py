# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Mingdian Liu
"""Lock the cupy boundary: the numpy-only modules must not import cupy, even transitively.

``__init__.py`` states it: the solver chain (solver, device_tables, sources_setup, monitors_setup,
dispersion_setup, pitch, coeffs, fusion, readout, device, plus the optional GPU path in tfsf1d;
setup_tables is only a re-export shim over those) may import cupy. Every other module sees numpy
only, because they have to run on a machine without a GPU, for writeback, scoring and export.

Same approach as tests/test_no_tidy3d_leak.py: import in a subprocess and inspect ``sys.modules``.
"""

import subprocess
import sys

#: The numpy-only modules; not one of them may touch cupy
NUMPY_ONLY = ["openem.model", "openem.grid", "openem.cpml", "openem.waveform",
              "openem.serialize", "openem.flux", "openem.fold", "openem.normalize",
              "openem.colocate", "openem.td_readout", "openem.shutoff",
              "openem.apodization", "openem.tfsf_oblique", "openem.tfsf1d",
              "openem.results", "openem.units", "openem.knobs"]


def test_numpy_only_modules_do_not_import_cupy():
    code = (
        "import " + ", ".join(NUMPY_ONLY) + "\n"
        "import sys\n"
        "leaked = sorted(m for m in sys.modules if m.split('.')[0] == 'cupy')\n"
        "assert not leaked, 'these modules dragged cupy in: ' + repr(leaked)\n"
        "print('OK')\n"
    )
    r = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True)
    assert r.returncode == 0, f"stdout={r.stdout}\nstderr={r.stderr}"
