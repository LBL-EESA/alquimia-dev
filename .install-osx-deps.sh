# Install required software
brew update
brew install openmpi wget

# Make sure the weird gfortran library links are in place.
ln -s /usr/local/lib/gcc/5/libgfortran.dylib /usr/local/lib/libgfortran.dylib
ln -s /usr/local/lib/gcc/5/libgfortran.a /usr/local/lib/libgfortran.a

# Go get PETSc 3.6.x and build it.
wget http://ftp.mcs.anl.gov/pub/petsc/release-snapshots/petsc-lite-3.6.3.tar.gz
tar xzf petsc-lite-3.6.3.tar.gz
pushd $PETSC_DIR
./configure --with-mpi=1 --with-debug=$DEBUG --with-shared-libraries=1
make
ln -s $PETSC_DIR/lib/petsc/conf $PETSC_DIR/conf
ln -s $PETSC_DIR/include/petsc/finclude $PETSC_DIR/include/finclude
popd

# Go get pflotran and build it.
wget https://bitbucket.org/pflotran/pflotran-dev/get/3fe478242357.zip
unzip -q 3fe478242357.zip
mv pflotran-pflotran-dev-3fe478242357 $PFLOTRAN_DIR
pushd $PFLOTRAN_DIR/src/pflotran
make pflotran_rxn
popd

