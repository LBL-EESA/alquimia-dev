cd $ALQUIMIA_DIR
if [ -d "${ALQUIMIA_DIR}/build" ]; then
  rm -rf "${ALQUIMIA_DIR}/build"
fi
ORIGIN_DIR=${pwd}
mkdir ${ALQUIMIA_DIR}/build ; cd ${ALQUIMIA_DIR}/build
cmake .. \
  -DCMAKE_INSTALL_PREFIX=$ALQUIMIA_DIR \
  -DCMAKE_C_COMPILER=`which mpicc` \
  -DCMAKE_CXX_COMPILER=`which mpicxx` \
  -DCMAKE_Fortran_COMPILER=`which mpif90` \
  -DCMAKE_BUILD_TYPE=Debug \
  -DXSDK_WITH_PFLOTRAN=OFF \
  -DXSDK_WITH_CRUNCHFLOW=ON \
  -DTPL_CRUNCHFLOW_LIBRARIES=$CRUNCHFLOW_DIR/source/libcrunchchem.a \
  -DTPL_CRUNCHFLOW_INCLUDE_DIRS=$CRUNCHFLOW_DIR/source \
  -DCMAKE_BUILD_TYPE=DEBUG \
  -DCMAKE_C_FLAGS="-W -Wall -Wextra" \
  -DCMAKE_CXX_FLAGS="-W -Wall -Wextra" \
  -DCMAKE_Fortran_FLAGS="-W -Wall -Wextra"
make -j 6 VERBOSE=1
make test
make install
cd $ORIGIN_DIR

