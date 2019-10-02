# Install necessary software.
sudo apt-get update -qq
sudo apt-get install -y cmake gcc libopenmpi-dev openmpi-bin liblapack-dev gfortran git
export TMPDIR=/tmp

# Go get PETSc and build it.
#git clone https://bitbucket.org/petsc/petsc petsc
git clone https://gitlab.com/petsc/petsc.git petsc
pushd $PETSC_DIR
git checkout v3.11.3
./configure --with-mpi=1 --download-hdf5 --with-debug=$DEBUG --with-shared-libraries=1 --download-pflotran --download-pflotran-commit=origin/master
make
##popd

# Get and build CrunchTope, inside external packages for convenience
cd $PETSC_ARCH/externalpackages
git clone https://bitbucket.org/crunchflow/crunchtope-dev.git --branch petsc-upgrade
cd crunchtope-dev/source
git apply MakefileForAlquimia.patch
make libcrunchchem.a
popd

# Go get pflotran and build it.
# -already installed under PETSc
