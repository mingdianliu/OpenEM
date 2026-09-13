# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Mingdian Liu
"""Table-building helpers. They now live in four modules, split by what they are about; this file
only re-exports them (a shim):

- ``sources_setup.py``    tables for plane waves, TFSF boxes, mode sources, dipoles, modulated media
- ``monitors_setup.py``   buffers for the four monitor kinds, DFT phase tables, emission geometry
                          (including the three table builders that used to sit in readout)
- ``dispersion_setup.py`` reordering of dispersion entries (dedup, bucketing, deferral, partitioning)
                          and the absorber tables
- ``pitch.py``            array utilities: z pitch, padding, upload, pinned copies

``_pml_axis_arrays`` is in coeffs.py and ``_refuse_unsupported_bloch`` is in solver.py.
Existing spellings such as ``from openem.setup_tables import _mode_source`` keep working, but new
code should import from the four modules above directly.
"""

from openem.coeffs import _pml_axis_arrays  # noqa: F401
from openem.pitch import (  # noqa: F401
    _d2h_pinned,
    _f32,
    _f64,
    _pad32_dev_to,
    _pad32_host_to,
    _pad_z_to,
    _padz1_to,
    _pitch,
    _plain_neighbors,
    _repitch_cells,
)
from openem.sources_setup import (  # noqa: F401
    _MODE_SIGNS,
    _dipole_table,
    _dipole_table_c,
    _mode_source,
    _mode_source_entry,
    _modulation_setup,
    _plane_source,
    _plane_source_one,
    _tfsf_box,
)
from openem.monitors_setup import (  # noqa: F401
    _DFT_TMOD,
    _N2F_TAN,
    _dft_table,
    _field_mon_entry,
    _flush_launch,
    _flux_mon_entry,
    _flux_time_buffers,
    _monitor_buffers,
    _n2f_cmap,
    _pick_T,
    _tmon_entry,
    _use_sm,
)
from openem.dispersion_setup import (  # noqa: F401
    _absorber_disp_decay,
    _absorber_slabs,
    _ade_defer_setup,
    _coef_lut,
    _entry_tables,
    _gather_poles,
    _identity_span,
    _lean_disp_partition,
    _lean_geometry,
    _lean_setup,
    _lean_side_out,
    _pole_bucket_setup,
    _shell_boxes,
    _shell_segments,
)
