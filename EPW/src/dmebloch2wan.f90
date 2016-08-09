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
  SUBROUTINE dmebloch2wan ( nbnd, nbndsub, nks, nkbl, dme, xk, cu, &
     nrr, irvec, wslen )
  !--------------------------------------------------------------------------
  !!
  !!  From the Dipole in Bloch representationi (coarse mesh), 
  !!  find the corresponding Dipole in Wannier representation 
  !!
  !
  USE kinds,         ONLY : DP
  USE pwcom,         ONLY : at, bg, celldm
  USE elph2,         ONLY : cdmew
  USE constants_epw, ONLY : bohr2ang, twopi, ci, czero
  USE io_global,     ONLY : ionode_id
  USE mp_global,     ONLY : inter_pool_comm
  USE mp_world,      ONLY : mpime
  USE mp,            ONLY : mp_barrier,mp_sum
  implicit none
  !
  !  input variables
  !
  INTEGER, INTENT (in) :: nbnd
  !! number of bands
  INTEGER, INTENT (in) :: nbndsub
  !! number of bands in the optimal subspace
  INTEGER, INTENT (in) :: nks
  !! number of kpoints
  INTEGER, INTENT (in) :: nkbl
  !! number of kpoint blocks, in the pool
  INTEGER, INTENT (in) :: nrr
  !! number of WS points 
  INTEGER, INTENT (in) :: irvec (3, nrr) 
  !! Coordinate of Wannier space points
  ! 
  REAL(kind=DP), INTENT (in) :: xk (3, nks)
  !! kpoint coordinates (cartesian in units of 2piba) 
  REAL(kind=DP), INTENT (in) :: wslen (nrr)
  !! WS vectors length (alat units)
  ! 
  COMPLEX(kind=DP), INTENT (in) :: dme (3,nbnd, nbnd,nks)
  !! Dipole matrix elements on coarse mesh
  COMPLEX(kind=DP), INTENT (in) :: cu (nbnd, nbndsub, nks)
  !! rotation matrix from wannier code
  !
  ! Local variables
  INTEGER :: ipol
  !! Counter on polarization
  INTEGER :: ik
  !! Counter on k-point
  INTEGER :: ir
  !! Counter on WS points
  REAL(kind=DP) :: rdotk
  !! $$ mathbf{r}\cdot\mathbf{k} $$
  REAL(kind=DP) :: tmp
  !! Temporary variables 
  COMPLEX(kind=DP) :: cps(3, nbndsub, nbndsub, nks)
  !! Hamiltonian in smooth Bloch basis, coarse mesh 
  COMPLEX(kind=DP) :: cfac
  !! $$ e^{-i\mathbf{r}\cdot\mathbf{k}} $$
  COMPLEX(kind=DP) :: dme_utmp(nbnd,nbndsub)
  !!
  !
  !----------------------------------------------------------
  !    STEP 1: rotation to optimally smooth Bloch states
  !----------------------------------------------------------
  !
  CALL start_clock ( 'Dipole: step 1' )
  !
  !  p~ (k) = U(k)^\dagger * p(k) * U(k)
  !  p~ (k) is cps( ipol, nbndsub, nbndsub, ik )
  !
  DO ik = 1, nks
     DO ipol = 1, 3
        ! copied from ephbloch2wane.  produce equivalent results
        !
        dme_utmp(:,:) = matmul( dme(ipol,:,:,ik), cu(:,:,ik) )
        cps (ipol, :,:, ik) = matmul ( conjg(transpose( cu(:,:,ik))), dme_utmp (:,:) )
        !
    ENDDO
  ENDDO
  !
  !
  !----------------------------------------------------------
  !  STEP 2: Fourier transform to go into Wannier basis
  !----------------------------------------------------------
  !
  !  p (R) = (1/nk) sum_k e^{-ikR} p~ (k)
  !  cdmew (ipol, nbndsub, nbndsub, ir) is p (R)
  !
  CALL start_clock ( 'Dipole: step 2' )
  !
  ! bring xk in crystal coordinates
  !
  CALL cryst_to_cart (nks, xk, at, -1)
  !
  cdmew ( :, :, :, :) = czero 
  ! 
  DO ir = 1, nrr
    !
    DO ik = 1, nks
       !
       rdotk = twopi * dot_product( xk ( :, ik), dble(irvec( :, ir) ))
       cfac = exp( -ci*rdotk ) / dble(nkbl)
       DO ipol = 1, 3
          cdmew ( ipol, :, :, ir ) = cdmew ( ipol, :, :, ir ) + cfac * cps ( ipol, :, :, ik )
       ENDDO
       !
    ENDDO
    !
  ENDDO
  CALL mp_sum(cdmew,inter_pool_comm) 
  !
  ! bring xk back into cart coord
  !
  CALL cryst_to_cart (nks, xk, bg, 1)
  !
  !
  !  check spatial decay of Dipole in Wannier basis
  !  the unit in r-space is angstrom
  !
  IF (mpime.eq.ionode_id) then
     OPEN(unit=300,file='decay.P')
     WRITE(300, '(/3x,a/)') '#Spatial decay of Dipole in Wannier basis'
     DO ir = 1, nrr
        !
        tmp =  maxval ( abs( cdmew (:, :,:,ir)) )
        WRITE(300, *) wslen(ir) * celldm (1) * bohr2ang, tmp
        !
     ENDDO
     close(300)
  ENDIF
  CALL mp_barrier(inter_pool_comm)
  !
  CALL stop_clock ( 'Dipole: step 2' )
  !
  END SUBROUTINE dmebloch2wan

