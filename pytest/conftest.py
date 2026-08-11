"""
Pytest configuration for PyOmega_h tests.

The PyOmega_h bindings now use shared_ptr<Library> with py::keep_alive
policies on Mesh and Comm objects. This ensures the Library outlives all
dependent objects, so normal Python GC cleanup works correctly without
needing os._exit().
"""
