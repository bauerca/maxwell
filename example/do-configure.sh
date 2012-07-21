#!/bin/bash

rm CMakeCache.txt
rm -rf CMakeFiles

cmake \
  -D MAXWELL_USE_MPI:BOOL=ON \
  -D MAXWELL_MPI_BASE_DIR:FILEPATH=/usr/local/mpi \
  $@ \
  ../
