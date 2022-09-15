#!/usr/bin/env bash
#
# Must be called, after Fortnet had been installed
#
# Expects following env-vars:
#
# FC, CC, ("false"/"true")
# SOURCE_DIR (Fortnet source dir), BUILD_DIR,
# INSTALL_DIR (Fortnet install dir with already installed Fortnet)
#
set -ex

# Integration test for CMake builds
CMAKE_PREFIX_PATH="${INSTALL_DIR}:${CMAKE_PREFIX_PATH}" \
    ${SOURCE_DIR}/test/prog/fortnet/integration/cmake/runtest.sh ${BUILD_DIR}/cmake

# Integration test for PKG-CONFIG builds
PKG_CONFIG_PATH="${INSTALL_DIR}/lib/pkgconfig:$PKG_CONFIG_PATH" \
    ${SOURCE_DIR}/test/prog/fortnet/integration/pkgconfig/runtest.sh ${BUILD_DIR}/pkgconfig
