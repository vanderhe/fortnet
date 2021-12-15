#:include 'common.fypp'

submodule (fnet_acsf) acsf_derivs

  implicit none
    
  contains

  module subroutine GeoAcsfGrad(this, geo, gFunctions, localAtToAtNum, extFeatures)
    type(TRealArray4D), intent(inout) :: this
    type(TGeometry)   , intent(in)    :: geo               !> system geometry container
    type(TGFunction), intent(in)      :: gFunctions(:)     !> list of multiple G-functions
    integer, intent(in)               :: localAtToAtNum(:) !> index mapping local atom --> atomic number
    real(dp), intent(in), optional    :: extFeatures(:,:)  !> atom dependent scaling parameters for cutoff function

    !> geometries reduced to atoms of a single species
    type(TGeometry) :: geo1, geo2

    !> atomic prefactors of central atom
    real(dp) :: atomId1, atomId2

    !> atomic prefactors of neighboring atoms
    real(dp), allocatable :: atomIds1(:), atomIds2(:)

    !> neighbor coordinates and squared distances
    real(dp), allocatable :: neighDists1(:), neighDists2(:), neighCoords1(:,:), neighCoords2(:,:)

    !> (?) global indices of neghbouring atoms in the reduced geometry 
    integer, allocatable  :: atomIndices1(:), atomIndices2(:)

    integer :: iAtom, iAcsf, iAtomOut1, i
    
    do iAtom = 1, geo%nAtom
      do iAcsf = 1, size(gFunctions)

        call buildGFunctionNeighborlists(iAtom, geo, gFunctions(iAcsf), localAtToAtNum , &
                      & extFeatures, geo1, geo2, iAtomOut1, atomId1, atomId2, atomIds1, atomIds2, &
                      & neighDists1, neighDists2, neighCoords1, neighCoords2, atomIndices1, atomIndices2)
        
        this%array(:, iAcsf, iAtom, atomIndices1) = gFuncGrad(gFunctions(iAcsf), geo%coords(:,iAtom), AtomId1, &
                                                                        & neighCoords1, neighDists1, atomIds1, &
                                                                        & neighCoords2, neighDists2, atomIds2)
            
        this%array(:, iAcsf, iAtom, iAtom) = -sum(this%array(:, iAcsf, iAtom, atomIndices1), dim=2)

        if (gFunctions(iAcsf)%tAngular) then
          if (gFunctions(iAcsf)%atomicNumbers(1) /= gFunctions(iAcsf)%atomicNumbers(2)) then
             this%array(:, iAcsf, iAtom, atomIndices2) = gFuncGrad(gFunctions(iAcsf), geo%coords(:,iAtom), AtomId1, &
                                                                             & neighCoords2, neighDists2, atomIds2, &
                                                                             & neighCoords1, neighDists1, atomIds1)

             this%array(:, iAcsf, iAtom, iAtom) = this%array(:, iAcsf, iAtom, iAtom) &
                                                 & - sum(this%array(:, iAcsf, iAtom, atomIndices2), dim=2)
          end if
        end if
        
      end do
    end do

  end subroutine GeoAcsfGrad

  function gFuncGrad(gFunction, iCoords, iAtomId, neighCoords1, neighDists1, atomIds1, neighCoords2, neighDists2, atomIds2)
    type(TGFunction), intent(in)   :: gFunction
    real(dp), intent(in)           :: iCoords(3)
    real(dp), intent(in)           :: iAtomId
    real(dp), intent(in)           :: neighCoords1(:,:), neighDists1(:), atomIds1(:)
    real(dp), intent(in), optional :: neighCoords2(:,:), neighDists2(:), atomIds2(:)

    real(dp) :: gFuncGrad(3, size(neighDists1))

    real(dp) :: df_ij, dfc_ik, dfc_jk, fc_ik, fc_jk, dist_jk, a_ijk
    real(dp) :: T_ik(3,3), uvec_jk(3), uvec_ik(3)
    integer  :: jAtom, kAtom, i
    real(dp) :: rCut, xi
    logical  :: is_g5 
    real(dp), allocatable :: uvecs2(:,:), fc2(:), sqDists2(:)
    logical , allocatable  :: mask(:)
    
    gFuncGrad = 0.0_dp
    rCut = gFunction%rCut
    
    if (gFunction%type == 'g1') then
      do jAtom = 1, size(neighDists1)
        gFuncGrad(:,jAtom) = df_cutoff(neighDists1(jAtom), iAtomId, atomIds1(jAtom), rCut, 1) & 
                             & * (neighCoords1(:,jAtom) - iCoords) / neighDists1(jAtom)
      end do
      return
    end if

    if (gFunction%type == 'g2') then
      do jAtom = 1, size(neighDists1)
        gFuncGrad(:,jAtom) = ( df_cutoff(neighDists1(jAtom), iAtomId, atomIds1(jAtom), rCut, 1)  &
                            & - f_cutoff(neighDists1(jAtom), iAtomId, atomIds1(jAtom), rCut)     &
                            &       * 2.0_dp*gFunction%eta*(neighDists1(jAtom) - gFunction%rs) ) & 
                            & * exp(-gFunction%eta*(neighDists1(jAtom) - gFunction%rs)**2) &
                            & * (neighCoords1(:,jAtom) - iCoords) / neighDists1(jAtom)
      end do
      return
    end if

    if (gFunction%type == 'g3') then
      associate (kappa => gFunction%kappa)
        do jAtom = 1, size(neighDists1)
          gFuncGrad(:,jAtom) = &
            &( cos(kappa* neighDists1(jAtom))* df_cutoff(neighDists1(jAtom), iAtomId, atomIds1(jAtom), rCut, 1) &
            & -sin(kappa* neighDists1(jAtom))*  f_cutoff(neighDists1(jAtom), iAtomId, atomIds1(jAtom), rCut)* kappa )  &
            &* (neighCoords1(:,jAtom) - iCoords) / neighDists1(jAtom)
        end do
      end associate
      return
    end if
    
    is_g5 = (gFunction%type == 'g5')
    xi    = gFunction%xi
    
    allocate(uvecs2  , source=neighCoords2)
    allocate(sqDists2, source=neighDists2)
    allocate(fc2, mold=neighDists2); fc2(:) = 0.0_dp
    do jAtom = 1, size(neighDists2)
        uvecs2(:,jAtom) = (uvecs2(:,jAtom) - iCoords) / neighDists2(jAtom)
        sqDists2(jAtom) = neighDists2(jAtom)**2
        if (neighDists2(jAtom) > rCut) cycle
        fc2(jAtom)      = f_cutoff(neighDists2(jAtom), iAtomId, atomIds2(jAtom), rCut)
    end do
    
    ! first handle the case where both neighbourhoods are the same
    ! this happens when there's no species resolution or when Z1 == Z2 
    if (gFunction%atomicNumbers(1) == gFunction%atomicNumbers(2)) then
      allocate(mask(size(neighDists2)))
      
      do kAtom = 1, size(neighDists1)
        T_ik   = T_matr(uvecs2(:,kAtom), neighDists1(kAtom)) * gFunction%lambda * xi
        dfc_ik = -2.0_dp * gFunction%eta * neighDists1(kAtom) * fc2(kAtom) & 
                & + df_cutoff(neighDists1(kAtom), iAtomId, atomIds1(kAtom), rCut, 1)
        
        ! case kAtom == jAtom
        gFuncGrad(:,kAtom) = gFuncGrad(:,kAtom) + & 
                         & (1.0_dp + gFunction%lambda)**xi * fc2(kAtom)* dfc_ik * &
                         & exp(-2.0_dp*gFunction%eta*sqDists2(kAtom)) * uvecs2(:,kAtom)
        
        ! TODO: perhaps "pack" function would be better here                 
        mask(:) = .false.
        mask(kAtom) = .true.
                
        do jAtom = 1, size(neighDists2)
          if (mask(jAtom)) cycle  ! skip jAtom == kAtom
          
          a_ijk = 1.0_dp + gFunction%lambda*dot_product(uvecs2(:,jAtom), uvecs2(:,kAtom))
          
          if (is_g5) then
            gFuncGrad(:,kAtom) = gFuncGrad(:,kAtom) + &
                               & fc2(jAtom) * exp(-gFunction%eta*(sqDists2(kAtom) + sqDists2(jAtom))) &
                               & *a_ijk**(xi - 1.0_dp)* (fc2(kAtom)*matmul(T_ik, uvecs2(:, jAtom)) + a_ijk*dfc_ik*uvecs2(:, kAtom))
            cycle
          end if

          uvec_jk = neighCoords1(:,kAtom) - neighCoords2(:,jAtom)
          dist_jk = norm2(uvec_jk)
          
          if (dist_jk > rCut) cycle
          
          uvec_jk = uvec_jk / dist_jk
          fc_jk   = f_cutoff(dist_jk, AtomIds1(kAtom), atomIds2(jAtom), rCut)
          dfc_jk  = -2.0_dp*gFunction%eta*dist_jk*fc_jk + df_cutoff(dist_jk, AtomIds1(kAtom), atomIds2(jAtom), rCut, 1)

          gFuncGrad(:,kAtom) = gFuncGrad(:,kAtom) + &
                           & fc2(jAtom) * exp(-gFunction%eta*(sqDists2(kAtom) + sqDists2(jAtom) + dist_jk**2)) &
                           & * a_ijk**(xi - 1.0_dp)*( fc2(kAtom)*fc_jk * matmul(T_ik, uvecs2(:, jAtom)) + &
                                                     & a_ijk *(fc_jk*dfc_ik*uvecs2(:, kAtom) + fc2(kAtom)*dfc_jk*uvec_jk) ) 
        end do

      end do
      gFuncGrad = gFuncGrad*2.0_dp**(2.0_dp - xi)
      return
    end if
    
    ! finally, case Z1 /= Z2
    do kAtom = 1, size(neighDists1)
      
      fc_ik   = f_cutoff(neighDists1(kAtom), iAtomId, atomIds1(kAtom), rCut)
      dfc_ik  = -2.0_dp * gFunction%eta * neighDists1(kAtom) * fc_ik &
                 & + df_cutoff(neighDists1(kAtom), iAtomId, atomIds1(kAtom), rCut, 1)
      
      uvec_ik = (neighCoords1(:,kAtom) - iCoords) / neighDists1(kAtom)
      T_ik    = T_matr(uvec_ik, neighDists1(kAtom)) * gFunction%lambda * xi
      
      do jAtom = 1, size(neighDists2)
        a_ijk = 1.0_dp + gFunction%lambda* dot_product(uvecs2(:,jAtom), uvec_ik)
        
        if (is_g5) then
          gFuncGrad(:,kAtom) = gFuncGrad(:,kAtom) + &
                          & fc2(jAtom)* exp(-gFunction%eta*(neighDists1(kAtom)**2 + sqDists2(jAtom))) &
                          & * a_ijk**(xi - 1.0_dp)*( fc_ik*matmul(T_ik, uvecs2(:, jAtom)) + a_ijk*dfc_ik*uvec_ik )
          cycle
        end if
        
        uvec_jk = neighCoords1(:,kAtom) - neighCoords2(:,jAtom)
        dist_jk = norm2(uvec_jk)
        
        if (dist_jk > rCut) cycle
        
        uvec_jk = uvec_jk / dist_jk
        fc_jk   = f_cutoff(dist_jk, AtomIds1(kAtom), atomIds2(jAtom), rCut)
        dfc_jk  = -2.0_dp*gFunction%eta*dist_jk*fc_jk + df_cutoff(dist_jk, AtomIds1(kAtom), atomIds2(jAtom), rCut, 1)

        gFuncGrad(:,kAtom) = gFuncGrad(:,kAtom) + &
                         & fc2(jAtom) * exp(-gFunction%eta*(neighDists1(kAtom)**2 + sqDists2(jAtom) + dist_jk**2)) &
                         & * a_ijk**(xi - 1.0_dp)*( fc_ik * fc_jk * matmul(T_ik, uvecs2(:, jAtom))  + &
                                                  & a_ijk *(fc_jk*dfc_ik*uvec_ik + fc_ik*dfc_jk*uvec_jk) ) 
      end do
    end do
    gFuncGrad = gFuncGrad*2.0_dp**(1.0_dp - xi)

  end function gFuncGrad


  pure function T_matr(uvec, nrm)
    real(dp), intent(in)  :: uvec(3), nrm
    real(dp) :: T_matr(3,3)
    
    integer  :: i, j

    T_matr = 0.0_dp
    do i = 1, 3
      T_matr(i,i) = 1.0_dp - uvec(i)*uvec(i)
      do j = i+1, 3
        T_matr(i,j) = -uvec(i)*uvec(j)
        T_matr(j,i) = T_matr(i,j)
      end do
    end do

    T_matr = T_matr / nrm
  end function
  
! *** N.B.: in f_cutoff and df_cutoff there is no checking whether rr < rcut ***
  pure function f_cutoff(rr, atomId1, atomId2, rcut) 
    !> atom distance (in cutoff range)
    real(dp), intent(in) :: rr

    !> atom ID of center atom
    real(dp), intent(in) :: atomId1

    !> atom ID of neighbor
    real(dp), intent(in) :: atomId2

    !> cutoff radius
    real(dp), intent(in) :: rcut

    !> corresponding activation function values
    real(dp) :: f_cutoff

    f_cutoff = 0.5_dp * atomId1 * atomId2 * (cos(pi * rr / rcut) + 1.0_dp)
  end function f_cutoff

  pure function df_cutoff(rr, atomId1, atomId2, rcut, deriv) 
    !> atom distance (in cutoff range)
    real(dp), intent(in) :: rr

    !> atom ID of center atom
    real(dp), intent(in) :: atomId1

    !> atom ID of neighbor
    real(dp), intent(in) :: atomId2

    !> cutoff radius
    real(dp), intent(in) :: rcut
    
    !> derivative order
    integer, intent(in)  :: deriv
    
    real(dp) :: df_cutoff
    
    ! here we use the formula for the n-th derivative of cos(ax)
    df_cutoff = 0.5_dp * atomId1 * atomId2 * (pi/rcut)**deriv * cos((pi*rr/rcut) + 0.5_dp*(deriv*pi))
  end function df_cutoff
  
end submodule acsf_derivs
