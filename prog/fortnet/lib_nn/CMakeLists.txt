set(curdir "lib_nn")

set(sources-fpp
  ${curdir}/bpnn.F90
  ${curdir}/layer.F90
  ${curdir}/network.F90
  ${curdir}/transfer.F90)

set(ALL-SOURCES-FPP ${ALL-SOURCES-FPP} ${sources-fpp} PARENT_SCOPE)
