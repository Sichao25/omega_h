cmake_minimum_required(VERSION 3.10)
project(TestOmegaHConfig CXX)

set(CMAKE_CXX_EXTENSIONS OFF)

# Enable testing
enable_testing()

# Find Omega_h package
find_package(Omega_h REQUIRED)

# Print found Omega_h configuration
message(STATUS "Omega_h_VERSION: ${Omega_h_VERSION}")
message(STATUS "Omega_h_USE_MPI: ${Omega_h_USE_MPI}")
message(STATUS "Omega_h_USE_Kokkos: ${Omega_h_USE_Kokkos}")
message(STATUS "Omega_h_USE_ZLIB: ${Omega_h_USE_ZLIB}")
message(STATUS "Omega_h_USE_libMeshb: ${Omega_h_USE_libMeshb}")
message(STATUS "Omega_h_USE_Gmsh: ${Omega_h_USE_Gmsh}")
message(STATUS "Omega_h_USE_Gmodel: ${Omega_h_USE_Gmodel}")
message(STATUS "Omega_h_USE_SEACASExodus: ${Omega_h_USE_SEACASExodus}")
message(STATUS "Omega_h_USE_ADIOS2: ${Omega_h_USE_ADIOS2}")
message(STATUS "Omega_h_USE_pybind11: ${Omega_h_USE_pybind11}")
message(STATUS "Omega_h_USE_DOLFIN: ${Omega_h_USE_DOLFIN}")
message(STATUS "Omega_h_USE_EGADS: ${Omega_h_USE_EGADS}")

# Create test executable
add_executable(test_export test_export.cpp)

# Link with Omega_h
target_link_libraries(test_export Omega_h::omega_h)

# Set C++ standard to match Omega_h (C++17)
set_property(TARGET test_export PROPERTY CXX_STANDARD 17)
set_property(TARGET test_export PROPERTY CXX_STANDARD_REQUIRED ON)

# Add test
add_test(NAME test_export COMMAND test_export)
set_tests_properties(test_export PROPERTIES 
  PASS_REGULAR_EXPRESSION "Test completed successfully!")
