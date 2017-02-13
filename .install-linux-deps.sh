# Install necessary software.
sudo apt-get update -qq
sudo apt-get install -y cmake gcc libopenmpi-dev openmpi-bin liblapack-dev gfortran mercurial wget

# Go get PETSc and build it.
wget https://bitbucket.org/petsc/petsc.zip
unzip -q petsc.zip
pushd $PETSC_DIR
./configure --with-mpi=1 --with-debug=$DEBUG --with-shared-libraries=1 --download-pflotran
make
popd

# Go get pflotran and build it.
# -already installed under PETSc
