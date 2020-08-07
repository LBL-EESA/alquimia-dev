#!/bin/sh

git clone https://gitlab.com/petsc/petsc.git $PETSC_DIR
cd $PETSC_DIR
git checkout v3.11.3
./configure --with-mpi=1 --with-debug=1 --with-shared-libraries=1 --download-fblaslapack=1 --with-debug=$DEBUG
make
