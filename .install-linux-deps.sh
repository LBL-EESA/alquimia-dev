# Install necessary software.
sudo apt-get update -qq
sudo apt-get install -y cmake gcc libopenmpi-dev openmpi-bin liblapack-dev gfortran mercurial wget

# Go get PETSc 3.5.x and build it.
wget http://ftp.mcs.anl.gov/pub/petsc/release-snapshots/petsc-lite-3.5.4.tar.gz
tar xzvf petsc-lite-3.5.4.tar.gz
pushd $PETSC_DIR
./configure --with-mpi=1 --with-debug=$DEBUG --with-shared-libraries=1
make
popd

# Go get pflotran and build it.
wget https://bitbucket.org/pflotran/pflotran-dev/get/611092f80ddb.zip
unzip 611092f80ddb.zip
mv pflotran-pflotran-dev-611092f80ddb $PFLOTRAN_DIR
pushd $PFLOTRAN_DIR/src/pflotran
make pflotran_rxn
popd
