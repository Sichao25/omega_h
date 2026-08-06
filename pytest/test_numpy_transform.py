"""
Direct unit tests for numpy_to_omega_h_read and omega_h_write_to_numpy.

Verifies that data is properly deep-copied through host memory
(not aliased from numpy buffer or device pointer).

Run with: pytest test_numpy_transform.py -v
"""

import numpy as np
import PyOmega_h as omega_h


def test_numpy_to_read_copies_data():
    """numpy_to_omega_h_read must copy data through host, not alias."""
    size = 10
    original = np.arange(size, dtype=np.float64)
    read_view = omega_h.numpy_to_read_float64(original)

    # Access Read data via HostRead (buffer protocol)
    host_read = omega_h.HostRead_float64(read_view)
    result = np.array(host_read, copy=False)
    np.testing.assert_array_equal(result, original)

    # Mutate original — result must be unchanged (deep copy)
    original[:] = -999.0
    result_after = np.array(host_read, copy=False)
    expected = np.arange(size, dtype=np.float64)
    np.testing.assert_array_equal(result_after, expected)


def test_write_to_numpy_copies_data():
    """omega_h_write_to_numpy must copy device data to host numpy array."""
    size = 10
    expected = np.arange(size, dtype=np.float64)

    host_write = omega_h.HostWrite_float64(size)
    # Fill via buffer protocol
    view = np.array(host_write, copy=False)
    view[:] = expected

    # HostWrite.write() returns Write; convert to numpy
    write_view = host_write.write()
    result = omega_h.write_to_numpy_float64(write_view)
    np.testing.assert_array_equal(result, expected)

    # Mutate host — result must be unchanged (deep copy)
    view[:] = -999.0
    np.testing.assert_array_equal(result, expected)

