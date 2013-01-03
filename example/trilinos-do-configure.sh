rm ./CMakeCache.txt
rm -rf ./CMakeFiles/

cmake \
  -D CMAKE_INSTALL_PREFIX:PATH=/usr/local/trilinos-11.0.3 \
  -D CMAKE_BUILD_TYPE:STRING=RELEASE \
  -D TPL_ENABLE_MPI:BOOL=OFF \
  -D MPI_BASE_DIR:PATH=/usr/local/mpi \
  -D Trilinos_ENABLE_ALL_PACKAGES:BOOL=OFF \
  -D Trilinos_ENABLE_TESTS:BOOL=OFF \
  -D Trilinos_ENABLE_Anasazi:BOOL=ON \
  -D Trilinos_ENABLE_ML:BOOL=ON \
  -D Trilinos_ENABLE_AztecOO:BOOL=ON \
  -D Trilinos_ENABLE_Ifpack:BOOL=ON \
  -D Trilinos_ENABLE_Amesos:BOOL=ON \
  -D Trilinos_ENABLE_Komplex:BOOL=ON \
  -D Trilinos_ENABLE_EpetraExt:BOOL=ON \
  -D Trilinos_ENABLE_Epetra:BOOL=ON \
  -D Trilinos_ENABLE_Zoltan:BOOL=ON \
  -D Trilinos_ENABLE_Teuchos:BOOL=ON \
  -D Teuchos_ENABLE_COMPLEX:BOOL=ON \
  -D Teuchos_ENABLE_EXTENDED:BOOL=ON \
  -D Trilinos_ENABLE_Triutils:BOOL=ON \
  -D Trilinos_ENABLE_Fortran:BOOL=OFF \
  ../
