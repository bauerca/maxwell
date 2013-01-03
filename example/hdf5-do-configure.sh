#!/bin/bash

rm CMakeCache.txt
rm -rf CMakeFiles/

cmake \
  -D CMAKE_INSTALL_PREFIX:PATH=/scr_verus/bauerca/prog \
  -D CMAKE_BUILD_TYPE:STRING=Release \
  -D HDF5_ENABLE_Z_LIB_SUPPORT:BOOL=ON \
  $@ \
  ../
  
# Add to above if using MPI
#  -D HDF5_ENABLE_PARALLEL:BOOL=OFF \
#  -D CMAKE_CXX_COMPILER:FILEPATH=/usr/local/mpi/bin/mpicxx \
#  -D CMAKE_C_COMPILER:FILEPATH=/usr/local/mpi/bin/mpicc \
