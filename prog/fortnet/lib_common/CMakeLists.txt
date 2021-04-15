set(curdir "lib_common")

set(sources-fpp
  ${curdir}/loss.F90
  ${curdir}/nestedtypes.F90
  ${curdir}/parallel.F90
  ${curdir}/random.F90
  ${curdir}/workarounds.F90)

set(ALL-SOURCES-FPP ${ALL-SOURCES-FPP} ${sources-fpp} PARENT_SCOPE)
