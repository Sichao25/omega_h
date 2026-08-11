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


def test_numpy_to_read_sliced_step():
    """numpy_to_omega_h_read must handle non-contiguous sliced arrays."""
    base = np.arange(20, dtype=np.float64)
    sliced = base[::3]  # [0, 3, 6, 9, 12, 15, 18] — non-contiguous
    assert not sliced.flags["C_CONTIGUOUS"]

    read_view = omega_h.numpy_to_read_float64(sliced)
    host_read = omega_h.HostRead_float64(read_view)
    result = np.array(host_read, copy=False)
    np.testing.assert_array_equal(result, sliced)

    # Verify deep copy — mutating base doesn't affect result
    base[:] = -999.0
    result_after = np.array(host_read, copy=False)
    expected = np.array([0, 3, 6, 9, 12, 15, 18], dtype=np.float64)
    np.testing.assert_array_equal(result_after, expected)


def test_numpy_to_read_sliced_reversed():
    """numpy_to_omega_h_read must handle reversed (negative stride) arrays."""
    base = np.array([1, 2, 3, 4, 5], dtype=np.int32)
    sliced = base[::-1]  # [5, 4, 3, 2, 1] — negative stride
    assert not sliced.flags["C_CONTIGUOUS"]

    read_view = omega_h.numpy_to_read_int32(sliced)
    host_read = omega_h.HostRead_int32(read_view)
    result = np.array(host_read, copy=False)
    np.testing.assert_array_equal(result, sliced)

    # Verify deep copy
    base[:] = 0
    result_after = np.array(host_read, copy=False)
    expected = np.array([5, 4, 3, 2, 1], dtype=np.int32)
    np.testing.assert_array_equal(result_after, expected)


def test_numpy_to_read_sliced_offset():
    """numpy_to_omega_h_read must handle sliced arrays with offset + step."""
    base = np.array([10, 11, 12, 13, 14, 15, 16], dtype=np.int64)
    sliced = base[2:6:2]  # [12, 14] — offset by 2, step 2
    assert not sliced.flags["C_CONTIGUOUS"]

    read_view = omega_h.numpy_to_read_int64(sliced)
    host_read = omega_h.HostRead_int64(read_view)
    result = np.array(host_read, copy=False)
    np.testing.assert_array_equal(result, sliced)

    # Verify deep copy
    base[2] = 999
    result_after = np.array(host_read, copy=False)
    expected = np.array([12, 14], dtype=np.int64)
    np.testing.assert_array_equal(result_after, expected)

