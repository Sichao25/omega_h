"""
Pytest configuration for PyOmega_h tests.

The lib fixture is defined here (conftest.py) so that exactly the same
Omega_h::Library instance exists for the entire test session. Each
test file must NOT define its own lib fixture — doing so creates duplicate
Library instances, and whichever one is destroyed first will call
MPI_Finalize, causing the other Library's Comm destructors to fail with:
"Attempting to use an MPI routine after finalizing MPICH".
"""

import atexit
import os
import sys

import pytest
import PyOmega_h as omega_h

_global_lib = None


@pytest.fixture(scope="session")
def lib():
    """Session-scoped library instance shared across all tests.

    The Library must outlive all Mesh and Comm objects. The pybind11
    bindings use shared_ptr with py::keep_alive to ensure this ordering
    through normal Python GC reference counting.
    """
    global _global_lib
    if _global_lib is None:
        _global_lib = omega_h.OmegaHLibrary()
    return _global_lib


@pytest.fixture(scope="session")
def world(lib):
    """Session-scoped communicator."""
    return lib.world()
