PETSC_DIR=/lore/seols/develop/petsc-3.21.4
PETSC_ARCH=real-gcc12.3.0-mpich4.1.1
PREFIX=/lore/seols/romulus-install
#module use /opt/scorec/spack/rhel9/v0201_4/lmod/linux-rhel9-x86_64/Core/
#module load gcc/12.3.0-iil3lno cmake/3.26.3-2duxfcd mpich/4.1.1-xpoyz4t zlib
#module load zstd/1.5.5-jmx6qld zlib/1.2.13-mvjz5oi libtirpc
#module load zlib/1.2.13-mjocrm2
cmake .. \
  -DCMAKE_C_COMPILER=mpicc \
  -DCMAKE_CXX_COMPILER=mpicxx \
  -DOmega_h_USE_CUDA=OFF \
  -DOmega_h_USE_Kokkos=OFF \
  -DOmega_h_USE_MPI=ON \
  -DBUILD_TESTING=ON \
  -DOmega_h_USE_ADIOS2=ON \
  -DADIOS2_INCLUDE_DIR=$PREFIX/include \
  -DADIOS2_LIB_DIR==$PREFIX/lib64 \
  -DCMAKE_INSTALL_PREFIX=$PREFIX \
  -DCMAKE_BUILD_TYPE=Debug
