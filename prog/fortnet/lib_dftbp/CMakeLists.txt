set(curdir "lib_dftbp")

set(sources-fpp
  ${curdir}/accuracy.F90
  ${curdir}/assert.F90
  ${curdir}/blas.F90
  ${curdir}/blasroutines.F90
  ${curdir}/charmanip.F90
  ${curdir}/conjgrad.F90
  ${curdir}/constants.F90
  ${curdir}/dynneighlist.F90
  ${curdir}/fire.F90
  ${curdir}/globalenv.F90
  ${curdir}/hsdparser.F90
  ${curdir}/hsdutils.F90
  ${curdir}/hsdutils2.F90
  ${curdir}/latpointiter.F90
  ${curdir}/lbfgs.F90
  ${curdir}/linemin.F90
  ${curdir}/linkedlist.F90
  ${curdir}/linkedlisti0.F90
  ${curdir}/linkedlisti1.F90
  ${curdir}/linkedlistlc0.F90
  ${curdir}/linkedlistmc0.F90
  ${curdir}/linkedlistr0.F90
  ${curdir}/linkedlistr1.F90
  ${curdir}/linkedlistr2.F90
  ${curdir}/linkedlists0.F90
  ${curdir}/message.F90
  ${curdir}/mpifx.F90
  ${curdir}/ranlux.F90
  ${curdir}/simplealgebra.F90
  ${curdir}/sorting.F90
  ${curdir}/steepdesc.F90
  ${curdir}/tokenreader.F90
  ${curdir}/typegeometry.F90
  ${curdir}/typegeometryhsd.F90
  ${curdir}/unitconversion.F90
  ${curdir}/xmlf90.F90
  ${curdir}/xmlutils.F90)

if(WITH_MPI)
  list(APPEND sources-fpp
    ${curdir}/mpienv.F90)
endif()

set(ALL-SOURCES-FPP ${ALL-SOURCES-FPP} ${sources-fpp} PARENT_SCOPE)
