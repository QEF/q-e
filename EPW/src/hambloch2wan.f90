  !                                                                            
  ! Copyright (C) 2010-2016 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino
  ! Copyright (C) 2007-2009 Jesse Noffsinger, Brad Malone, Feliciano Giustino  
  !                                                                            
  ! This file is distributed under the terms of the GNU General Public         
  ! License. See the file `LICENSE' in the root directory of the               
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .             
  !                                                                            
  !
  !--------------------------------------------------------------------------
  SUBROUTINE hambloch2wan ( nbnd, nbndsub, nks, nkstot, et, xk, cu, &
     lwin, nrr, irvec, wslen, chw )
  !--------------------------------------------------------------------------
  !!
  !!  From the Hamiltonian in Bloch representationi (coarse mesh), 
  !!  find the corresponding Hamiltonian in Wannier representation 
  !!
  !--------------------------------------------------------------------------
  !
  USE kinds,     ONLY : DP
  USE pwcom,     ONLY : at, bg, celldm
  USE constants_epw, ONLY : bohr2ang, twopi, ci, czero
  USE io_global, ONLY : ionode_id
  USE mp_global, ONLY : inter_pool_comm
  USE mp,        ONLY : mp_barrier,mp_sum
  USE mp_world,  ONLY : mpime
  !
  implicit none
  !
  LOGICAL, INTENT (in) :: lwin( nbnd, nks )
  !! identify bands within outer energy window (for disentanglement)
  !
  INTEGER, INTENT (in) :: nbnd
  !! number of bands 
  INTEGER, INTENT (in) :: nbndsub
  !! number of bands in the optimal subspace 
  INTEGER, INTENT (in) :: nks
  !! number of kpoints
  INTEGER, INTENT (in) :: nkstot
  !! number of kpoint blocks, in the pool
  INTEGER, INTENT (in) ::  nrr
  !! number of kpoint blocks, total 
  INTEGER, INTENT (in) :: irvec (3, nrr)
  !! number of WS points and coordinates
  !
  REAL(kind=DP), INTENT (in) :: et (nbnd, nks)
  !! hamiltonian eigenvalues, coarse mesh
  REAL(kind=DP), INTENT (in) :: xk (3, nks)
  !! kpoint coordinates (cartesian in units of 2piba)
  REAL(kind=DP), INTENT (in) :: wslen (nrr)
  !! WS vectors length (alat units)
  !
  COMPLEX(kind=DP), INTENT (in) :: cu (nbnd, nbndsub, nks)
  !! rotation matrix from wannier code
  COMPLEX(kind=DP), INTENT(OUT) :: chw ( nbndsub, nbndsub, nrr)
  !! Hamiltonian in smooth Bloch basis, coarse mesh 
  !
  ! work variables 
  !
  INTEGER :: ik
  !! Counter of k-point index
  INTEGER :: ibnd
  !! Counter on band index
  INTEGER :: jbnd
  !! Counter on band index
  INTEGER :: ir
  !! Counter on real-space index
  INTEGER :: mbnd
  !! Counter on band index 
  INTEGER :: j
  !! Counter on band index
  REAL(kind=DP) :: rdotk
  !! $$\mathbf{r}\cdot\mathbf{k}
  REAL(kind=DP) :: tmp
  !! Maximum value of the real space Hamiltonian
  REAL(kind=dp)    :: et_opt(nbnd,nks)
  !! hamiltonian eigenvalues within the outer window in the first ndimwin(ik) entries
  !
  COMPLEX(kind=DP) :: chs(nbndsub, nbndsub, nks) 
  !! Hamiltonian
  COMPLEX(kind=DP) :: cfac
  !! $$e^{-i\mathbf{r}\cdot\mathbf{k}}$$
  COMPLEX(kind=DP) :: ctmp
  !! Temporary variable to store the Hamiltonian
  !
  !----------------------------------------------------------
  !    STEP 1: rotation to optimally smooth Bloch states
  !----------------------------------------------------------
  !
  CALL start_clock ( 'Ham: step 1' )
  !
  et_opt=0.d0
  !
  ! slim down et to contain states within the outer window
  !
  do ik=1,nks
     ibnd=0
     do j=1,nbnd
        if(lwin(j,ik)) then
           ibnd=ibnd+1
           et_opt(ibnd,ik)=et(j,ik)
        end if
     end do
  end do
  !
  !  H~ (k) = U(k)^\dagger * H(k) * U(k)
  !  H~ (k) is chs( nbndsub, nbndsub, ik )
  !
  DO ik = 1, nks
   !
   !
   DO jbnd = 1, nbndsub
    DO ibnd = 1, jbnd
       !
       ctmp = czero
       !
       DO mbnd = 1, nbnd
         ctmp = ctmp + conjg(cu (mbnd,ibnd,ik)) * et_opt (mbnd,ik) * cu (mbnd,jbnd,ik)
       ENDDO
       !
       chs (ibnd , jbnd , ik) = ctmp 
       chs (jbnd , ibnd , ik) = conjg(ctmp)
       !
    ENDDO
   ENDDO
  ENDDO
  !
  CALL stop_clock ( 'Ham: step 1' )
  !
  !----------------------------------------------------------
  !  STEP 2: Fourier transform to go into Wannier basis
  !----------------------------------------------------------
  !
  !  H (R) = (1/nk) sum_k e^{-ikR} H~ (k)
  !  chw (nbndsub, nbndsub, ir) is H (R)
  !
  CALL start_clock ( 'Ham: step 2' )
  !
  ! bring xk in crystal coordinates
  !
  CALL cryst_to_cart (nks, xk, at, -1)
  !
  chw ( :, :, :) = czero 
  ! 
  DO ir = 1, nrr
    !
    DO ik = 1, nks
       !
       rdotk = twopi * dot_product( xk ( :, ik), dble(irvec( :, ir) ))
       cfac = exp( -ci*rdotk ) / dble(nkstot)
       chw ( :, :, ir ) = chw ( :, :, ir ) + cfac * chs ( :, :, ik )
       !
    ENDDO
    !
  ENDDO
  CALL mp_sum(chw,inter_pool_comm) 
  !
  ! bring xk back into cart coord
  !
  CALL cryst_to_cart (nks, xk, bg, 1)
  !
    !
    !  check spatial decay of Hamiltonian in Wannier basis
    !  the unit in r-space is angstrom
    !
    IF (mpime.eq.ionode_id) THEN
       open(unit=300,file='decay.H')
       WRITE(300, '(/3x,a/)') '#Spatial decay of Hamiltonian in Wannier basis'
      DO ir = 1, nrr
        !
        tmp =  maxval ( abs( chw (:,:,ir)) )
        WRITE(300, *) wslen(ir) * celldm (1) * bohr2ang, tmp
        !
      ENDDO
      close(300)
    ENDIF
    CALL mp_barrier(inter_pool_comm)
  !
  CALL stop_clock ( 'Ham: step 2' )
  !
  END SUBROUTINE hambloch2wan

