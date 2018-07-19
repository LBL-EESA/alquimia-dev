# Install necessary software.
sudo apt-get update -qq
sudo apt-get install -y cmake gcc libopenmpi-dev openmpi-bin liblapack-dev gfortran git
export TMPDIR=/tmp

# Go get PETSc and build it.
git clone https://bitbucket.org/petsc/petsc petsc
pushd $PETSC_DIR
./configure --with-mpi=1 --with-debug=$DEBUG --with-shared-libraries=1 --download-pflotran --download-pflotran-commit=origin/master
make
popd

# Go get pflotran and build it.
# -already installed under PETSc
