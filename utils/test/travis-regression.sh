#!/usr/bin/env bash
#
# Expects following env-vars:
#
# FC, CC, WITH_MPI ("false"/"true")
# SOURCE_DIR, BUILD_DIR, INSTALL_DIR
#
set -ex

if [ "$TRAVIS_PULL_REQUEST" != "false" ]; then
  BUILD_TYPE="Debug"
else
  BUILD_TYPE="Release"
fi

cmake_options=(
  "-DCMAKE_INSTALL_PREFIX=${INSTALL_DIR}"
  "-DWITH_MPI=${WITH_MPI}"
  "-DWITH_API=true"
  "-DFYPP_FLAGS='-DTRAVIS'"
  "-DHYBRID_CONFIG_METHODS='Submodule'"
  "-DCMAKE_BUILD_TYPE=${BUILD_TYPE}"
)

cmake -B ${BUILD_DIR} "${cmake_options[@]}" ${SOURCE_DIR}
cmake --build ${BUILD_DIR} -- -j
pushd ${BUILD_DIR}
ctest -j --output-on-failure
popd
cmake --install ${BUILD_DIR}
