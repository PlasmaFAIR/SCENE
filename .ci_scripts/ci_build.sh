set -ex

cmake . -B build \
      -DGHOST_DIR=${GHOST_INSTALL_DIR}
cmake --build build
