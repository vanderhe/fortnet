set(curdir "lib_extlibs")

set(sources-fpp
  ${curdir}/mpifx.F90)

if(WITH_SOCKETS)
  list(APPEND sources-fpp ${curdir}/fsockets.F90)
endif(WITH_SOCKETS)

set(ALL-SOURCES-FPP ${ALL-SOURCES-FPP} ${sources-fpp} PARENT_SCOPE)
