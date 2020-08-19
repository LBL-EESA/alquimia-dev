#!/bin/sh

git clone https://gitlab.com/petsc/petsc.git $PETSC_DIR
cd $PETSC_DIR
git checkout balay/fix-make-f90-target/maint
./configure --with-mpi=1 --with-debug=1 --with-shared-libraries=1 --download-fblaslapack=1 --with-debug=$DEBUG
make
