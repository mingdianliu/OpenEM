# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Mingdian Liu
"""Sits at the repo root so pytest puts this directory on sys.path, which is what makes
``from openem import ...`` work.

(pytest defaults to importmode=prepend, which inserts the directory holding conftest.py at the
front of sys.path.)
"""

import importlib.util

import pytest

#: Without tidy3d-extras installed, every test case that needs subpixel averaging skips instead of
#: failing. tidy3d-extras is Flexcompute's add-on and public environments, CI included, usually lack
#: it. Its absence only affects the "subpixel-averaged epsilon" path; the solver itself and the rest
#: of the tests do not depend on it.
HAS_TIDY3D_EXTRAS = importlib.util.find_spec("tidy3d_extras") is not None
_EXTRAS_HINT = "tidy3d-extras"


@pytest.hookimpl(hookwrapper=True)
def pytest_runtest_makereport(item, call):
    report = (yield).get_result()
    if HAS_TIDY3D_EXTRAS or not report.failed or call.excinfo is None:
        return
    # Covers both the directly raised case and the one where a test meant to assert some other
    # exception fails because the missing-extras error was raised first.
    exc = call.excinfo.value
    if _EXTRAS_HINT in str(exc) or _EXTRAS_HINT in str(getattr(exc, "__cause__", "")):
        report.outcome = "skipped"
        report.longrepr = (str(item.fspath), item.location[1] or 0,
                           "skipped: needs tidy3d-extras (subpixel averaging)")
