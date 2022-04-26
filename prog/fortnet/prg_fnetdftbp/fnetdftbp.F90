!--------------------------------------------------------------------------------------------------!
!  FORTNET: A Behler-Parrinello-Neural-Network Implementation                                      !
!  Copyright (C) 2020 - 2022  T. W. van der Heide                                                  !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Demonstrates possible usage of Fortnet within DFTB+.
program fnetdftbp

  use dftbp_accuracy, only : dp
  use dftbp_message, only : error
  use dftbp_globalenv, only : initGlobalEnv, destructGlobalEnv, stdOut
  use dftbp_typegeometry, only : TGeometry, normalize
  use dftbp_simplealgebra, only : invert33, determinant33

  use fnet_initprogram, only : TProgramVariables, printFortnetHeader, printDateAndTime,&
      & readFromNetstat, copyrightYear, version
  use fnet_nestedtypes, only : TEnv, TEnv_init, TPredicts, TPredicts_init, TIntArray1D,&
      & TRealArray1D, TRealArray2D, TJacobians, TJacobian
  use fnet_forces, only : TForces, forceAnalysis
  use fnet_dftbpfeatures, only : TFeatures_init_dftbp, TFeatures_collect_dftbp
  use fnet_acsf, only : TAcsf
  use fnet_bpnn, only : TBpnn
  use fnet_netstat, only : readExtFeaturesConfig, inquireExtFeatures, inquireAcsf
  use fnet_intmanip, only : getNumberOfUniqueInt

#:if WITH_MPI
  use fnet_mpifx
#:endif

  implicit none


  !> instance containing program variables
  type(TProgramVariables) :: prog

  !> Behler-Parrinello-Neural-Network instance
  type(TBpnn) :: bpnn

  !> representation of ACSF mappings and derivatives
  type(TAcsf) :: acsf, acsfPrime

  !> determinant to calculate inverse lattice vectors
  real(dp) :: det

  !> true, if current mpi process is the lead
  logical :: tLead

  !> runs over atoms
  integer :: iAtom

  !! #######################Input from DFTB+########################################################

  !> program version
  character(len=*), parameter :: netstatpath = 'fortnet.hdf5'

  !> geometry to calculate (should later be provided by DFTB+)
  type(TGeometry) :: geo

  !> index mapping local atom --> atomic number
  integer, allocatable :: atomicNumbers(:)

  !> Mulliken populations for each atom
  !! shape: mullikenPopulations(nAtoms)%array(nMullikenPopulations)
  type(TRealArray1D), allocatable :: mullikenPopulations(:)

  !! #######################Output for DFTB+########################################################

  !> obtained global prediction (i.e. total energy contribution)
  real(dp) :: globalPrediction

  !> obtained atomic forces, shape: [nAtom]
  real(dp), allocatable :: atomicPredictions(:)

  !> obtained atomic forces, shape: [3, nAtom]
  real(dp), allocatable :: atomicForces(:,:)

  !> derivative of global prediction (in our case the total energy) w.r.t.
  !! external input features (in our case Mulliken populations)
  !! shape: mullikenPopulations(nUniqueElements)%array(nMullikenPopulations)
  type(TRealArray1D), allocatable :: dEnergydMullikenpop(:)

  ! write standard output header
  call printFortnetHeader(version, copyrightYear)
  call printDateAndTime()

  ! initialise global environment
  call initGlobalEnv()

  ! initialise program variables
  call TEnv_init(prog%env)

#:if WITH_MPI
  tLead = prog%env%globalMpiComm%lead
#:else
  tLead = .true.
#:endif

  ! set up example Si2 geometry (geo is of standard DFTB+ type)
  geo%tPeriodic = .true.
  ! no support for helical boundary conditions
  geo%tHelical = .false.
  geo%tFracCoord = .true.
  geo%nSpecies = 1
  allocate(geo%species(2))
  geo%species(:) = 1
  allocate(geo%coords(3, 2))
  geo%coords(:, 1) = [0.00_dp, 0.00_dp, 0.00_dp]
  geo%coords(:, 2) = [0.25_dp, 0.25_dp, 0.25_dp]
  geo%nAtom = size(geo%coords, dim=2)
  if (geo%tPeriodic) then
    allocate(geo%latVecs(3, 3))
    geo%latvecs(:, 1) = [0.0_dp, 1.0_dp, 1.0_dp]
    geo%latvecs(:, 1) = geo%latvecs(:, 1) * 5.3_dp
    geo%latvecs(:, 2) = [1.0_dp, 0.0_dp, 1.0_dp]
    geo%latvecs(:, 2) = geo%latvecs(:, 2) * 5.3_dp
    geo%latvecs(:, 3) = [1.0_dp, 1.0_dp, 0.0_dp]
    geo%latvecs(:, 3) = geo%latvecs(:, 3) * 5.3_dp
    geo%origin = [0.0_dp, 0.0_dp, 0.0_dp]
    ! useless shift:
    geo%coords(:,:) = geo%coords - spread(geo%origin, 2, geo%nAtom)
    if (geo%tFracCoord) then
      geo%coords(:,:) = matmul(geo%latVecs, geo%coords)
      geo%origin(:) = matmul(geo%latVecs, geo%origin)
    end if
    allocate(geo%recVecs2p(3, 3))
    det = determinant33(geo%latVecs)
    if (abs(det) < 1e-12_dp) then
      call error('Linear dependent lattice vectors (RIP).')
    end if
    call invert33(geo%recVecs2p, geo%latVecs, det)
    geo%recVecs2p(:,:) = reshape(geo%recVecs2p, [3, 3], order=[2, 1])
  end if
  call normalize(geo)

  ! we also need the atomic number of the atoms in geo
  atomicNumbers = [14, 14]

  ! handle possible external features (Mulliken populations)
  ! case 1: ACSF only (hand over Mulliken populations anyway, Fortnet will not use them)
  ! case 2: ACSF + Mulliken populations as additional atomic input features
  allocate(mullikenPopulations(geo%nAtom))
  ! insert some dummy Mulliken populations (s, p, d) for testing
  do iAtom = 1, geo%nAtom
    mullikenPopulations(iAtom)%array = [2.0_dp, 2.0_dp, 0.0_dp]
  end do

  ! initializes the BPNN and ACSF instance (only call once)
  call handleInitialisation(prog, netstatpath, bpnn, acsf)

  ! calculate ACSF mappings for the current geometry (call everytime the geometry changed)
  call calculateMappings(prog%env, geo, atomicNumbers, mullikenPopulations, acsf, acsfPrime)

  ! evaluate network and Jacobian/forces
  call predict(prog, bpnn, acsf, acsfPrime, geo, mullikenPopulations, globalPrediction,&
      & atomicPredictions, atomicForces, dEnergydMullikenpop)

  write(stdOut, '(1F20.10)') globalPrediction
  write(stdOut, '(2F20.10)') atomicPredictions
  write(stdOut, '(3F20.10)') atomicForces
  ! write(stdOut, '(3F20.10)') dEnergydMullikenpop(1)%array

  ! #############################Geometry-Change####################################################

  geo%coords(:, 1) = [0.00_dp, 0.00_dp, 0.00_dp]
  geo%coords(:, 2) = [0.25_dp, 0.25_dp, 0.25_dp]
  if (geo%tPeriodic) then
    geo%latvecs(:, 1) = [0.0_dp, 1.0_dp, 1.0_dp]
    geo%latvecs(:, 1) = geo%latvecs(:, 1) * 7.3_dp
    geo%latvecs(:, 2) = [1.0_dp, 0.0_dp, 1.0_dp]
    geo%latvecs(:, 2) = geo%latvecs(:, 2) * 7.3_dp
    geo%latvecs(:, 3) = [1.0_dp, 1.0_dp, 0.0_dp]
    geo%latvecs(:, 3) = geo%latvecs(:, 3) * 7.3_dp
    geo%origin = [0.0_dp, 0.0_dp, 0.0_dp]
    ! useless shift:
    geo%coords(:,:) = geo%coords - spread(geo%origin, 2, geo%nAtom)
    if (geo%tFracCoord) then
      geo%coords(:,:) = matmul(geo%latVecs, geo%coords)
      geo%origin(:) = matmul(geo%latVecs, geo%origin)
    end if
    det = determinant33(geo%latVecs)
    if (abs(det) < 1e-12_dp) then
      call error('Linear dependent lattice vectors (RIP).')
    end if
    call invert33(geo%recVecs2p, geo%latVecs, det)
    geo%recVecs2p(:,:) = reshape(geo%recVecs2p, [3, 3], order=[2, 1])
  end if
  call normalize(geo)

  ! if geometry changed, re-evaluate...
  call calculateMappings(prog%env, geo, atomicNumbers, mullikenPopulations, acsf, acsfPrime)
  call predict(prog, bpnn, acsf, acsfPrime, geo, mullikenPopulations, globalPrediction,&
      & atomicPredictions, atomicForces, dEnergydMullikenpop)

  write(stdOut, '(1F20.10)') globalPrediction
  write(stdOut, '(2F20.10)') atomicPredictions
  write(stdOut, '(3F20.10)') atomicForces
  ! write(stdOut, '(3F20.10)') dEnergydMullikenpop(1)%array

  call destructGlobalEnv()


contains

  !> Initialises several derived types needed for running the core.
  subroutine handleInitialisation(prog, netstatpath, bpnn, acsf)

    !> representation of program variables
    type(TProgramVariables), intent(inout) :: prog

    !> filename or path to netstat file
    character(len=*), intent(in) :: netstatpath

    !> Behler-Parrinello-Neural-Network instance
    type(TBpnn), intent(out) :: bpnn

    !> representation of ACSF mappings
    type(TAcsf), intent(out) :: acsf

    !! true, if current process is the lead
    logical :: tLead

  #:if WITH_MPI
    tLead = prog%env%globalMpiComm%lead
  #:else
    tLead = .true.
  #:endif

    write(stdOut, '(/A,/)') 'Initialisation:'
    write(stdOut, '(A)') 'running in prediction mode'

    ! path to Fortnet's configuration file
    prog%inp%data%netstatpath = netstatpath

    if (tLead) then
      prog%inp%features%nFeatures = 0
      call inquireExtFeatures(prog%inp%data%netstatpath, prog%inp%features%tExtFeatures)
      write(stdOut, '(A,L1,/)') 'Network was trained with additional external features: ',&
          & prog%inp%features%tExtFeatures
      if (prog%inp%features%tExtFeatures) then
        call readExtFeaturesConfig(prog%inp%data%netstatpath, prog%inp%features%ext)
        prog%inp%features%nFeatures = prog%inp%features%nFeatures&
            & + prog%inp%features%ext%nExtFeatures
      end if
      call inquireAcsf(prog%inp%data%netstatpath, prog%inp%features%tMappingFeatures)
      if (prog%inp%features%tMappingFeatures) then
        call acsf%fromFile(prog%inp%data%netstatpath,&
            & tReduce=prog%inp%features%mapping%tReduce,&
            & tStandardize=prog%inp%features%mapping%tStandardize,&
            & nRadial=prog%inp%features%mapping%nRadial,&
            & nAngular=prog%inp%features%mapping%nAngular)
        if (prog%inp%features%mapping%nRadial + prog%inp%features%mapping%nAngular < 1) then
          call error('Parsed netstat file does not contain ACSF inputs. Aborting.')
        end if
        prog%inp%features%nFeatures = prog%inp%features%nFeatures&
            & + size(acsf%gFunctions%func)
      end if
      call bpnn%fromFile(prog%inp%data%netstatpath)
    end if
  #:if WITH_MPI
    call mpifx_bcast(prog%env%globalMpiComm, prog%inp%features%nFeatures)
    call mpifx_bcast(prog%env%globalMpiComm, prog%inp%features%tMappingFeatures)
    call mpifx_bcast(prog%env%globalMpiComm, prog%inp%features%tExtFeatures)
    if (prog%inp%features%tMappingFeatures) then
      call acsf%syncConfig(prog%env%globalMpiComm)
    end if
    if (prog%inp%features%tExtFeatures) then
      call prog%inp%features%ext%syncConfig(prog%env%globalMpiComm)
    end if
    call bpnn%sync(prog%env%globalMpiComm)
  #:endif

    write(stdout, '(A,/)') repeat('-', 80)

  end subroutine handleInitialisation


  !> Carries out the calculation of structural mappings, if desired.
  subroutine calculateMappings(env, geo, atomicNumbers, mullikenPopulations, acsf, acsfPrime)

    !> if compiled with mpi enabled, contains mpi communicator
    type(TEnv), intent(in) :: env

    !> geometry to calculate ACSF for
    type(TGeometry), intent(in) :: geo

    !> index mapping local atom --> atomic number
    integer, intent(in) :: atomicNumbers(:)

    !> Mulliken populations for each atom
    !! shape: mullikenPopulations(nAtoms)%array(nMullikenPopulations)
    type(TRealArray1D), intent(in) :: mullikenPopulations(:)

    !> representation of acsf mapping information
    type(TAcsf), intent(inout) :: acsf, acsfPrime

    !! contains geometries to calculate ACSF for (in our case with length of one)
    type(TGeometry), allocatable :: geos(:)

    !! index mapping local atom --> atomic number
    type(TIntArray1D), allocatable :: localAtToAtNum(:)

    !! additional external, atomic features for list of geometries
    type(TRealArray2D), allocatable :: extFeatures(:)

    !! runs over atoms
    integer :: iAtom

    ! create Fortnet compatible fields
    allocate(geos(1))
    geos(1) = geo
    allocate(localAtToAtNum(1))
    localAtToAtNum(1)%array = atomicNumbers

    if (prog%inp%features%tExtFeatures) then
      allocate(extFeatures(1))
      ! warning: this usage of mullikenPopulations only works for single-species systems
      allocate(extFeatures(1)%array(size(mullikenPopulations(1)%array), geo%nAtom))
      do iAtom = 1, geo%nAtom
        extFeatures(1)%array(:, iAtom) = mullikenPopulations(iAtom)%array
      end do
    end if

    write(stdOut, '(A)', advance='no') 'Calculating ACSF...'

    ! copy ACSF configurations in advance of any calculation
    acsfPrime = acsf

    call acsf%calculate(geos, env, localAtToAtNum, extFeaturesInp=extFeatures)
    call acsfPrime%calculatePrime(geos, env, localAtToAtNum, extFeaturesInp=extFeatures)

    write(stdout, '(A)') 'done'

  end subroutine calculateMappings


  !> Runs Fortnet's core routines, depending on the obtained configuration.
  subroutine predict(prog, bpnn, acsf, acsfPrime, geo, mullikenPopulations,&
      & globalPrediction, atomicPredictions, atomicForces, dEnergydMullikenpop)

    !> representation of program variables
    type(TProgramVariables), intent(inout) :: prog

    !> Behler-Parrinello-Neural-Network instance
    type(TBpnn), intent(inout) :: bpnn

    !> representation of acsf mapping information
    type(TAcsf), intent(inout) :: acsf, acsfPrime

    !> geometry to calculate
    type(TGeometry), intent(in) :: geo

    !> Mulliken populations for each atom
    !! shape: mullikenPopulations(nAtoms)%array(nMullikenPopulations)
    type(TRealArray1D), intent(in) :: mullikenPopulations(:)

    !> obtained global prediction (i.e. total energy contribution)
    real(dp), intent(out) :: globalPrediction

    !> obtained atomic forces, shape: [nAtom]
    real(dp), intent(out), allocatable :: atomicPredictions(:)

    !> obtained atomic forces, shape: [3, nAtom]
    real(dp), intent(out), allocatable :: atomicForces(:,:)

    !> derivative of global prediction (in our case the total energy) w.r.t.
    !! external input features (in our case Mulliken populations)
    !! shape: dEnergydMullikenpop(nUniqueElements)%array(nMullikenPopulations)
    type(TRealArray1D), intent(out), allocatable :: dEnergydMullikenpop(:)

    !! obtained atomic forces
    type(TForces) :: forces

    !! contains Jacobians of a single systems
    type(TJacobians) :: jacobian

    !! index mapping local atom --> global species index
    type(TIntArray1D), allocatable :: localAtToGlobalSp(:)

    !! additional external, atomic features of the dataset
    real(dp), allocatable :: extFeatures(:,:)

    !! true, if current process is the lead
    logical :: tLead

    !! number of datapoints
    integer :: nDatapoints

    !! number of ACSF mappings
    integer :: nAcsf

    !! number of additional external atomic input features
    integer :: nExtFeatures

    !! total number of atomic input features
    integer :: nFeatures

    !! number of unique elements in geometry
    integer :: nElements

    ! !! maximum ACSF cutoff
    ! real(dp) :: maxAcsfCutoff

  #:if WITH_MPI
    tLead = prog%env%globalMpiComm%lead
  #:else
    tLead = .true.
  #:endif

    ! we only want to calculate a single structure
    nDatapoints = 1

    if (bpnn%nAtomicTargets /= 0 .or. bpnn%nGlobalTargets /= 1) then
      call error('DFTB+Fortnet interfacing only working for exactly one global target.')
    end if

    ! create Fortnet compatible fields
    if (prog%inp%features%tExtFeatures) then
      ! warning: this usage of mullikenPopulations only works for single-species systems
      allocate(extFeatures(size(mullikenPopulations(1)%array), geo%nAtom))
      do iAtom = 1, geo%nAtom
        extFeatures(:, iAtom) = mullikenPopulations(iAtom)%array
      end do
      nExtFeatures = size(extFeatures, dim=1)
    else
      nExtFeatures = 0
    end if

    ! count total number of features
    nAcsf = size(acsf%gFunctions%func)
    nFeatures = nAcsf + nExtFeatures

    ! since there is only one geometry, no difference between local/global species indexing
    allocate(localAtToGlobalSp(1))
    localAtToGlobalSp(1)%array = geo%species

    call TFeatures_init_dftbp(prog%features, acsf, geo, nExtFeatures)

    if (tLead) then
      call TFeatures_collect_dftbp(prog%features, acsf, extFeatures=extFeatures)
      if (bpnn%dims(1) /= nFeatures) then
        call error('Mismatch in number of features and BPNN input nodes.')
      end if
    end if

  #:if WITH_MPI
    call prog%features%sync(prog%env%globalMpiComm)
  #:endif

    if (allocated(acsf%vals%vals)) deallocate(acsf%vals%vals)
    if (allocated(acsfPrime%vals%vals)) deallocate(acsfPrime%vals%vals)

    write(stdOut, '(A)', advance='no') 'Feeding networks...'
    atomicPredictions = reshape(bpnn%iPredict(prog%features%trainFeatures(1)%array, geo%species),&
        & [geo%nAtom])
    write(stdOut, '(A)') 'done'
    globalPrediction = sum(atomicPredictions)

    write(stdOut, '(A)', advance='no') 'Calculating Jacobian...'
    allocate(jacobian%sys(nDatapoints))
    jacobian%sys(1) = bpnn%iJacobian(prog%features%trainFeatures(1)%array, geo%species)
    write(stdOut, '(A)') 'done'
    write(stdOut, '(A)', advance='no') 'Calculating forces...'
    forces = forceAnalysis(bpnn, prog%features%trainFeatures, jacobian, acsfPrime, prog%env,&
        & localAtToGlobalSp, nExtFeaturesToIgnore=nExtFeatures)
    write(stdOut, '(A,/)') 'done'
    atomicForces = forces%geos(1)%array

    write(stdout, '(A,/)') repeat('-', 80)

    ! maxAcsfCutoff = acsf%getMaxCutoff()
    ! write(stdOut, '(1F20.10,A)') maxAcsfCutoff, ' Bohr'

    if (prog%inp%features%tExtFeatures) then
      call getNumberOfUniqueInt(geo%species, nElements)
      allocate(dEnergydMullikenpop(nElements))
      ! warning: this is only valid for single species systems
      dEnergydMullikenpop(1)%array = jacobian%sys(1)%atom(1)%array(1, nAcsf+1:nFeatures)
    end if

  end subroutine predict

end program fnetdftbp
