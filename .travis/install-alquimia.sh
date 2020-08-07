#!/bin/sh

cd $ALQUIMIA_DIR
mkdir build
cd build
cmake .. -DCMAKE_INSTALL_PREFIX=$PWD/../alquimia-install \
         -DCMAKE_C_COMPILER=`which mpicc` \
         -DCMAKE_CXX_COMPILER=`which mpicxx` \
         -DCMAKE_Fortran_COMPILER=`which mpif90` \
         -DCMAKE_BUILD_TYPE=$BUILD_TYPE \
         -DXSDK_WITH_PFLOTRAN=ON \
         -DTPL_PFLOTRAN_LIBRARIES=$PFLOTRAN_DIR/src/pflotran/libpflotranchem.a \
         -DTPL_PFLOTRAN_INCLUDE_DIRS=$PFLOTRAN_DIR/src/pflotran \
         -DXSDK_WITH_CRUNCHFLOW=ON \
         -DTPL_CRUNCHFLOW_LIBRARIES=$CRUNCHTOPE_DIR/source/libcrunchchem.a \
         -DTPL_CRUNCHFLOW_INCLUDE_DIRS=$CRUNCHTOPE_DIR/source
make -j4
make test
make install
make test_install
