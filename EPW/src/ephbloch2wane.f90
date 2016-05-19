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
  SUBROUTINE ephbloch2wane ( nbnd, nbndsub, nks, nkstot, lgamma, xk, &
       cu, cuq, lwin, lwinq, epmatk, nrr, irvec, wslen, epmatw)
  !-----------------------------------------------------------------------
  !
  !  From the electron-phonon matrix elements in Bloch representation (coarse 
  !  mesh), find the corresponding matrix elements in Wannier representation
  !
  !
  !-----------------------------------------------------------------------
  !
  USE kinds,     ONLY : DP
  USE pwcom,     ONLY : at, bg, celldm
  USE constants_epw, ONLY : bohr2ang, twopi, ci, czero, cone
#ifdef __PARA
  USE io_global, ONLY : ionode_id
  USE mp_global, ONLY : inter_pool_comm,my_pool_id
  USE mp       , ONLY : mp_sum 
  USE mp_world,  ONLY : mpime
#endif
  implicit none
  !
  !  input variables
  !
  integer :: nbnd, nbndsub, nks, nrr, irvec (3, nrr), nkstot
  ! number of bands
  ! number of bands in the optimal subspace
  ! number of kpoints
  ! number of kpoint blocks, in the pool
  ! number of kpoint blocks, total
  ! number of WS points and coordinates
  logical :: lgamma
  ! true if we are working with q=0
  real(kind=DP) :: xk (3, nks), wslen (nrr)
  ! kpoint coordinates (cartesian in units of 2piba)
  ! WS vectors length (alat units)
  complex(kind=DP) :: cu (nbnd, nbndsub, nks), cuq (nbnd, nbndsub, nks), epmatk ( nbnd, nbnd, nks)
  ! rotation matrix from wannier code
  ! e-p matrix in bloch representation, coarse mesh
  logical :: lwin( nbnd, nks ), lwinq( nbnd, nks )
  ! identify bands within outer energy window (for disentanglement)
  !
  ! output variables 
  !
  complex(kind=DP) :: epmatw ( nbndsub, nbndsub, nrr)
  !  e-p matrix  in wannier basis 
  !
  ! work variables 
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
#ifdef __PARA
  CALL mp_sum(epmatw,inter_pool_comm)  
#endif
  !
  ! bring xk back into cart coord
  !
  CALL cryst_to_cart (nks, xk, bg, 1)
  !
  !
  !  check spatial decay of matrix elements in Wannier basis
  !  the unit in r-space is angstrom, and I am plotting 
  !  the matrix for the first mode only
  !
#ifdef __PARA
    IF (mpime.eq.ionode_id) THEN
#endif
      OPEN (unit=301,file='decay.epwane')
      WRITE(301, '(/3x,a/)') '#Spatial decay of e-p matrix elements in Wannier basis'
      DO ir = 1, nrr
        ! 
        tmp =  maxval ( abs(epmatw(:,:,ir)) ) 
        WRITE(301, *) wslen(ir) * celldm (1) * bohr2ang, tmp
        !
      ENDDO
      CLOSE(301)
#ifdef __PARA
    ENDIF
#endif
  !
  CALL stop_clock ( 'ep: step 2' )
  !
  END SUBROUTINE ephbloch2wane

  
