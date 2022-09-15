#!/usr/bin/env bash
#
# Build external components first and uses them as external dependencies during build
#
# Expects following env-vars:
#
# WITH_MPI ("false"/"true"), SOURCE_DIR, BUILD_DIR, INSTALL_DIR
#
set -ex

if [ "${WITH_MPI}" == "true" ]; then

  cmake -DCMAKE_INSTALL_PREFIX=${INSTALL_DIR} \
        -DBUILD_EXPORTED_TARGETS_ONLY=True -DCMAKE_BUILD_TYPE=Debug \
        -B ${BUILD_DIR}/mpifx ${SOURCE_DIR}/external/mpifx/origin
  cmake --build ${BUILD_DIR}/mpifx -- -j
  cmake --install ${BUILD_DIR}/mpifx

  cmake -DWITH_MPI=True \
        -DHYBRID_CONFIG_METHODS='Find' \
        -DCMAKE_INSTALL_PREFIX=${INSTALL_DIR} \
        -DCMAKE_BUILD_TYPE=Debug \
        -B ${BUILD_DIR}/fortnet ${SOURCE_DIR}
  cmake --build ${BUILD_DIR}/fortnet -- VERBOSE=1 fnet

else

  cmake -DWITH_MPI=False \
        -DHYBRID_CONFIG_METHODS='Find' \
        -DBUILD_SHARED_LIBS=False \
        -DCMAKE_BUILD_TYPE=Debug \
        -DCMAKE_INSTALL_PREFIX=${INSTALL_DIR} \
        -B ${BUILD_DIR}/fortnet ${SOURCE_DIR}
  cmake --build ${BUILD_DIR}/fortnet -- VERBOSE=1 fnet

fi
