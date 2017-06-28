maxwell
=======

Distributed-memory parallel eigensolver for Maxwell's equations.

See the following publications for an idea of what this software can do.

- [A fast multigrid-based electromagnetic eigensolver for curved metal boundaries on the Yee mesh](http://carlbauer.org/static/pdf/dey-mittra-multigrid-solver.pdf)
- [A second-order 3D electromagnetics algorithm for curved interfaces between anisotropic dielectrics on a Yee mesh](http://carlbauer.org/static/pdf/second-order-yee-dielectrics.pdf)

EXAMPLE CONFIGURATION
=====================

::

  '/opt/contrib-clangcxx11/cmake-3.8.1-ser/bin/cmake' \
    -DCMAKE_INSTALL_PREFIX:PATH=/opt/volatile-clangcxx11/maxwell-r19-ser \
    -DCMAKE_BUILD_TYPE:STRING=RelWithDebInfo \
    -DCMAKE_COLOR_MAKEFILE:BOOL=FALSE \
    -DCMAKE_VERBOSE_MAKEFILE:BOOL=TRUE \
    -DCMAKE_C_COMPILER:FILEPATH='/usr/bin/clang' \
    -DCMAKE_CXX_COMPILER:FILEPATH='/usr/bin/clang++' \
    -DCMAKE_Fortran_COMPILER:FILEPATH='/opt/homebrew/bin/gfortran-4.9' \
    -DCMAKE_C_FLAGS:STRING='-fPIC -pipe' \
    -DCMAKE_CXX_FLAGS:STRING='-std=c++11 -stdlib=libc++ -fPIC -pipe' \
    -DCMAKE_Fortran_FLAGS:STRING='-fPIC -pipe' \
    -DTRILINOS_BASE_DIR:FILEPATH='/opt/contrib-clangcxx11/trilinos-11.14.3-sercomm' \
    -DUSE_PREBUILT_HDF5:BOOL=TRUE \
    -DHDF5_BASE_DIR:FILEPATH='/opt/contrib-clangcxx11/hdf5-1.8.18-ser' \
    /Users/userdir/projects/vorpalall/maxwell

LICENSE
=======

Copyright (c) 2012 Regents of the University of Colorado
This software was developed by the Plasmas and Beams Group at the Center for
Integrated Plasma Studies under U. S. Department of Energy grant
DE-FG02-04ER41317.

Permission is hereby granted, free of charge, to any person obtaining a copy of
this software and associated documentation files (the "Software"), to deal in
the Software without restriction, including without limitation the rights to
use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies
of the Software, and to permit persons to whom the Software is furnished to do
so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
