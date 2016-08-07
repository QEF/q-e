  !                                                                            
  ! Copyright (C) 2010-2016 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino 
  ! Copyright (C) 2007-2009 Jesse Noffsinger, Brad Malone, Feliciano Giustino  
  !                                                                            
  ! This file is distributed under the terms of the GNU General Public         
  ! License. See the file `LICENSE' in the root directory of the               
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .             
  !                                                                            
  !                                                                            
  ! 
  !-----------------------------------------------------------------------
  SUBROUTINE ephbloch2wane ( nbnd, nbndsub, nks, nkstot, xk, &
       cu, cuq, epmatk, nrr, irvec, wslen, epmatw)
  !-----------------------------------------------------------------------
  !!
  !!  From the electron-phonon matrix elements in Bloch representation (coarse 
  !!  mesh), find the corresponding matrix elements in Wannier representation
  !!
  !!
  !-----------------------------------------------------------------------
  !
  USE kinds,     ONLY : DP
  USE pwcom,     ONLY : at, bg, celldm
  USE constants_epw, ONLY : bohr2ang, twopi, ci, czero, cone
  USE io_epw,    ONLY : iuwane
  USE io_global, ONLY : ionode_id
  USE mp_global, ONLY : inter_pool_comm
  USE mp       , ONLY : mp_sum 
  USE mp_world,  ONLY : mpime
  implicit none
  !
  !  input variables
  !
  INTEGER, INTENT (in) :: nbnd
  !! Number of bands
  INTEGER, INTENT (in) :: nbndsub
  !! Number of bands in the optimal subspace
  INTEGER, INTENT (in) :: nks
  !! Number of kpoints
  INTEGER, INTENT (in) :: nrr
  !! Number of kpoint blocks, in the pool
  INTEGER, INTENT (in) :: nkstot
  !! Number of kpoint blocks, total
  INTEGER, INTENT (in) :: irvec (3, nrr)
  !! Number of WS points and coordinates
  !
  REAL(kind=DP), INTENT (in) :: xk (3, nks)
  !! kpoint coordinates (cartesian in units of 2piba)
  REAL(kind=DP), INTENT (in) :: wslen (nrr)
  !! WS vectors length (alat units)
  COMPLEX(kind=DP), INTENT (in) :: cu (nbnd, nbndsub, nks)
  !! Rotation matrix from wannier code
  COMPLEX(kind=DP), INTENT (in) :: cuq (nbnd, nbndsub, nks)
  !! Rotation matrix from wannier code
  COMPLEX(kind=DP), INTENT (in) :: epmatk ( nbnd, nbnd, nks)
  !! e-p matrix in bloch representation, coarse mesh
  !
  ! output variables 
  !
  COMPLEX(kind=DP), INTENT(out) :: epmatw ( nbndsub, nbndsub, nrr)
  !!  e-p matrix  in wannier basis 
  !
  ! Work variables 
  !
  complex(kind=DP) :: epmats (nbndsub, nbndsub, nks), eptmp(nbndsub, nbnd)
  !  e-p matrix  in smooth Bloch basis, coarse mesh
  !  e-p matrix, temporary
  !
  integer :: ik, ir
  real(kind=DP) :: rdotk, tmp
  complex(kind=DP) :: cfac
  !
  !
  !----------------------------------------------------------
  !  STEP 1: rotation to optimally smooth Bloch states
  !----------------------------------------------------------
  !
  !  g~ = U_k+q^\dagger g U_k
  !
  !  g   is epmatk (ibnd, jbnd, ik)
  !  g~  is epmats (ibnd, jbnd, ik)
  !
  CALL start_clock ( 'ep: step 1' )
  !
  !
  DO ik = 1, nks
     !
     ! the two zgemm calls perform the following ops:
     ! epmats  = [ cu(ikq)^\dagger * epmatk ] * cu(ikk)
     ! [here we have a size-reduction from nbnd*nbnd to nbndsub*nbndsub] 
     !
     CALL zgemm ('c', 'n', nbndsub, nbnd, nbnd, cone, cuq(:,:,ik),  &
          nbnd, epmatk(:,:,ik), nbnd, czero, eptmp, nbndsub)
     CALL zgemm ('n', 'n', nbndsub, nbndsub, nbnd, cone, eptmp,     &
          nbndsub, cu(:,:,ik), nbnd, czero, epmats(:,:,ik), nbndsub)
     !
  ENDDO
  !
  CALL stop_clock ( 'ep: step 1' )
  !
  !----------------------------------------------------------------------
  !  STEP 2: Fourier transform to obtain matrix elements in wannier basis
  !----------------------------------------------------------------------
  !
  !  g (R) = (1/nkc) sum_k e^{-ikR} g~(k)
  !
  !  epmatw (nbndsub,nbndsub,ir) is g(R)
  !
  CALL start_clock ( 'ep: step 2' )
  !
  epmatw (:, :, :) = czero
  !
  ! bring xk in crystal coordinates
  !
  CALL cryst_to_cart (nks, xk, at, -1)
  !
  DO ir = 1, nrr
     !
     DO ik = 1, nks
       !
       !
       rdotk = twopi * dot_product( xk ( :, ik), dble(irvec( :, ir) ))
       cfac = exp( -ci*rdotk ) / dble(nkstot)
       epmatw ( :, :, ir) = epmatw ( :, :, ir) + cfac * epmats ( :, :, ik)
       !
     ENDDO
     !
  ENDDO
  !
  CALL mp_sum(epmatw,inter_pool_comm)  
  !
  ! bring xk back into cart coord
  !
  CALL cryst_to_cart (nks, xk, bg, 1)
  !
  !
  !  Sheck spatial decay of matrix elements in Wannier basis
  !  the unit in r-space is angstrom, and I am plotting 
  !  the matrix for the first mode only
  !
  IF (mpime.eq.ionode_id) THEN
    OPEN (unit=iuwane,file='decay.epwane')
    WRITE(iuwane, '(a)') '# Spatial decay of e-p matrix elements in Wannier basis'
    DO ir = 1, nrr
      ! 
      tmp =  maxval ( abs(epmatw(:,:,ir)) ) 
      WRITE(iuwane, *) wslen(ir) * celldm (1) * bohr2ang, tmp
      !
    ENDDO
    CLOSE(iuwane)
  ENDIF
  !
  CALL stop_clock ( 'ep: step 2' )
  !
  END SUBROUTINE ephbloch2wane

  
