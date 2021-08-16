set(builddir ${CMAKE_CURRENT_BINARY_DIR})
set(srcdir ${CMAKE_CURRENT_SOURCE_DIR})
file(TO_NATIVE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/../../../" projectdir)

# execute_process is verbatim, therefore ";" must be replaced with " "
string(REGEX REPLACE ";" " " releasename "${RELEASE}")


execute_process(
  COMMAND ${projectdir}/utils/test/testlist_to_fypp
  INPUT_FILE ${srcdir}/tests
  OUTPUT_FILE ${builddir}/_fortnet_tests.fypp)

set(fypp_flags ${FYPP_CONFIG_FLAGS})
list(APPEND fypp_flags -I${projectdir}/prog/fortnet/include -DRELEASE="'${releasename}'")

execute_process(
  COMMAND ${FYPP} ${fypp_flags} -DMPI_PROCS=${TEST_MPI_PROCS}
  INPUT_FILE ${builddir}/_fortnet_tests.fypp
  OUTPUT_FILE ${builddir}/_fortnet_tests)

file(STRINGS ${builddir}/_fortnet_tests fortnet_tests_raw)
foreach(line IN LISTS fortnet_tests_raw)
  string(STRIP "${line}" testname)
  if(NOT "${testname}" STREQUAL "")
    list(APPEND tests ${testname})
  endif()
endforeach()

# copy datasets to build folder
execute_process(
  COMMAND ${CMAKE_COMMAND} -E copy_directory
  ${srcdir}/datasets ${builddir}/datasets)

foreach(test IN LISTS tests)
  add_test(
    NAME fortnet_${test}
    COMMAND ${PYTHON_INTERPRETER} ${srcdir}/bin/testwithworkdir ${test} ${builddir}
        -r ${srcdir} -P ${TEST_RUNNER} -p $<TARGET_FILE:fnet>)
endforeach()
