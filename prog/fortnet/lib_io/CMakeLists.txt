set(curdir "lib_io")

set(sources-fpp
  ${curdir}/fnetdata.F90
  ${curdir}/fnetout.F90
  ${curdir}/iterout.F90
  ${curdir}/netstat.F90)

set(ALL-SOURCES-FPP ${ALL-SOURCES-FPP} ${sources-fpp} PARENT_SCOPE)
