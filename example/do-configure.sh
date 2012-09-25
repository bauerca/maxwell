#!/bin/bash

rm CMakeCache.txt
rm -rf CMakeFiles

cmake \
  -D CMAKE_BUILD_TYPE:STRING=RELEASE \
  -D CMAKE_INSTALL_PREFIX:FILEPATH=/scr_verus/bauerca/prog \
  -D MAXWELL_USE_MPI:BOOL=ON \
  -D MAXWELL_MPI_BASE_DIR:FILEPATH=/usr/local/mpi \
  $@ \
  ../
