set -ex

cmake --version
cmake . -B build \
      -DGHOST_ROOT=${GHOST_INSTALL_DIR} \
      -DGHOST_DEBUG=ON \
      -DnetCDFFortran_DEBUG=ON \
      -DnetCDF_Fortran_INCLUDE_DIR=/usr/lib64/gfortran/modules/

cmake --build build
