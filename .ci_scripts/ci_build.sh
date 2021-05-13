set -ex

cmake --version
cmake . -B build \
      -DGHOST_ROOT=${GHOST_INSTALL_DIR} \
      -DGHOST_DEBUG=ON

cmake --build build
