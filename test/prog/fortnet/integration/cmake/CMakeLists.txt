cmake_minimum_required(VERSION 3.16)

project(TestFortnetBuild LANGUAGES Fortran)

find_package(Fortnet REQUIRED)

add_executable(test_build test_build.f90)
target_link_libraries(test_build Fortnet::Fortnet)


