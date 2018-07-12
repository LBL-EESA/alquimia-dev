# Install required software
brew update
brew install gcc
brew install open-mpi
export TMPDIR=/tmp

# Make sure the weird gfortran library links are in place.
ln -s /usr/local/lib/gcc/5/libgfortran.dylib /usr/local/lib/libgfortran.dylib
ln -s /usr/local/lib/gcc/5/libgfortran.a /usr/local/lib/libgfortran.a

# Go get PETSc 3.6.x and build it.
git clone https://bitbucket.org/petsc/petsc petsc
pushd $PETSC_DIR
./configure --with-mpi=1 --with-debugging=$DEBUG --with-shared-libraries=1 --download-pflotran --with-fc=`which mpif90`
make
make test
popd

# Go get pflotran and build it.
# -already installed under PETSc
