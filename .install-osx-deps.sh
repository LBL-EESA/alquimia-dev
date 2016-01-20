# Install required software
brew update
brew install openmpi wget

# Make sure the weird gfortran library links are in place.
ln -s /usr/local/lib/gcc/5/libgfortran.dylib /usr/local/lib/libgfortran.dylib
ln -s /usr/local/lib/gcc/5/libgfortran.a /usr/local/lib/libgfortran.a

# Go get PETSc 3.5.x and build it.
wget http://ftp.mcs.anl.gov/pub/petsc/release-snapshots/petsc-lite-3.5.4.tar.gz
tar xzvf petsc-lite-3.5.4.tar.gz
pushd $PETSC_DIR
./configure --with-mpi=1 --with-debug=$DEBUG --with-shared-libraries=0
make
popd

# Go get pflotran and build it.
wget https://bitbucket.org/pflotran/pflotran-dev/get/611092f80ddb.zip
unzip 611092f80ddb.zip
mv pflotran-pflotran-dev-611092f80ddb $PFLOTRAN_DIR
pushd $PFLOTRAN_DIR/src/pflotran
make pflotran_rxn
popd

