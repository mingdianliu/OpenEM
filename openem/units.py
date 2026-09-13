# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Mingdian Liu
"""Length units.

Tidy3D works in **um**; internally we work in **metres** everywhere. The conversion happens only
in the ``scene/`` layer, so units never come up again anywhere else.
"""

#: um -> m
UM = 1e-6
