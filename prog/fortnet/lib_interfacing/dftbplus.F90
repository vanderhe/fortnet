!--------------------------------------------------------------------------------------------------!
!  FORTNET: A Behler-Parrinello-Neural-Network Implementation                                      !
!  Copyright (C) 2020 - 2022  T. W. van der Heide                                                  !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Contains tweaks of some functionality, tailored towards DFTB+.
module fnet_dftbpfeatures

  use dftbp_accuracy, only: dp
  use dftbp_message, only : error
  use dftbp_typegeometry, only : TGeometry

  use fnet_acsf, only : TAcsf
  use fnet_features, only : tFeatures

#:if WITH_MPI
  use fnet_mpifx
#:endif

  implicit none

  private

  public :: TFeatures_init_dftbp, TFeatures_collect_dftbp


contains

  !> Initialises a feature instance.
  subroutine TFeatures_init_dftbp(features, acsf, geo, nExtFeatures)

    !> collected features
    type(TFeatures), intent(inout) :: features

    !> representation of acsf mapping informations
    type(TAcsf), intent(in) :: acsf

    !> geometry instance
    type(TGeometry), intent(in) :: geo

    !> number of additional external atomic input features
    integer, intent(in) :: nExtFeatures

    !! total number of input features
    integer :: nFeatures

    ! potentially dangerous, assumes that acsf are present/allocated
    nFeatures = size(acsf%gFunctions%func) + nExtFeatures

    if (nFeatures == 0) then
      call error('No features present to collect. Aborting.')
    end if

    ! we only want to calculate a single structure
    if (allocated(features%trainFeatures)) deallocate(features%trainFeatures)
    allocate(features%trainFeatures(1))
    allocate(features%trainFeatures(1)%array(nFeatures, geo%nAtom))
    features%trainFeatures(1)%array(:,:) = 0.0_dp

  end subroutine TFeatures_init_dftbp


  !> Collects features from ACSF calculations and external sources.
  subroutine TFeatures_collect_dftbp(features, acsf, extFeatures)

    !> collected features of data and mapping block
    type(TFeatures), intent(inout) :: features

    !> representation of acsf mapping informations
    type(TAcsf), intent(in) :: acsf

    !> additional external, atomic features of the dataset
    real(dp), intent(in), optional :: extFeatures(:,:)

    !! number of additional external atomic input features
    integer :: nExtFeatures

    !! number of ACSF mappings
    integer :: nAcsf

    !! total number of input features
    integer :: nFeatures

    if (present(extFeatures)) then
      nExtFeatures = size(extFeatures, dim=1)
    else
      nExtFeatures = 0
    end if
    nAcsf = size(acsf%gFunctions%func)
    nFeatures = nAcsf + nExtFeatures

    if ((nAcsf > 0) .and. (nExtFeatures > 0)) then
      features%trainFeatures(1)%array(1:nAcsf, :) = acsf%vals%vals(1)%array
      features%trainFeatures(1)%array(nAcsf+1:nFeatures, :) = extFeatures
    elseif ((nAcsf > 0) .and. (.not. (nExtFeatures > 0))) then
      features%trainFeatures(1)%array = acsf%vals%vals(1)%array
    elseif ((.not. (nAcsf > 0)) .and. (nExtFeatures > 0)) then
      features%trainFeatures(1)%array = extFeatures
    else
      call error('No features present to collect. Aborting.')
    end if

  end subroutine TFeatures_collect_dftbp

end module fnet_dftbpfeatures
