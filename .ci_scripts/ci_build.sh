set -ex

cmake --version
cmake . -B build \
      -DGHOST_ROOT=${GHOST_INSTALL_DIR} \
      -DGHOST_DEBUG=ON \
      -DnetCDFFortran_DEBUG=ON \
      -DnetCDF_Fortran_INCLUDE_DIR=/usr/lib64/gfortran/modules/ \
      -DBUILD_TESTING=ON

cmake --build build

cmake --build build --target check
