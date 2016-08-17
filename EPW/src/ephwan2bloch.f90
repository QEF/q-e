  !                                                                            
  ! Copyright (C) 2010-2016 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino 
  ! Copyright (C) 2007-2009 Jesse Noffsinger, Brad Malone, Feliciano Giustino  
  !                                                                            
  ! This file is distributed under the terms of the GNU General Public         
  ! License. See the file `LICENSE' in the root directory of the               
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .             
  !                                                                            
  ! 
  !---------------------------------------------------------------------------
  SUBROUTINE ephwan2bloch ( nbnd, nrr, irvec, ndegen, epmatw, &
         xk, cufkk, cufkq, epmatf, nmodes)
  !---------------------------------------------------------------------------
  !!
  !! Interpolation from Wannier to the fine Bloch grid of the electron-phonon 
  !! matrix elements
  !!
  USE kinds,         ONLY : DP
  USE constants_epw, ONLY : twopi, ci, czero, cone
  implicit none
  !
  INTEGER, INTENT (in) :: nbnd
  !! number of bands (possibly in the optimal subspace)
  INTEGER, INTENT (in) :: nrr
  !! Number of Wigner-Size points
  INTEGER, INTENT (in) :: irvec ( 3, nrr)
  !! Coordinates of WS points
  INTEGER, INTENT (in) :: ndegen (nrr)
  !! Degeneracy of WS points
  INTEGER, INTENT (in) :: nmodes
  !! number of phonon modes
  !
  REAL(kind=DP), INTENT (in) :: xk(3)
  !! kpoint for the interpolation (WARNING: this must be in crystal coord!)
  !
  COMPLEX(kind=DP), INTENT (in) :: epmatw ( nbnd, nbnd, nrr, nmodes)
  !! e-p matrix in Wannier representation
  COMPLEX(kind=DP), INTENT (in) :: cufkk (nbnd, nbnd)
  !! rotation matrix U(k)
  COMPLEX(kind=DP), INTENT (in) :: cufkq (nbnd, nbnd)
  !! rotation matrix U(k+q)
  COMPLEX(kind=DP), INTENT (out) :: epmatf (nbnd, nbnd, nmodes)
  !! e-p matrix in Bloch representation, fine grid
  !
  ! work variables 
  integer :: ir, imode
  real(kind=DP) :: rdotk
  complex(kind=DP) :: cfac, eptmp( nbnd, nbnd)
  !
  !----------------------------------------------------------
  !  STEP 3: inverse Fourier transform of g to fine k mesh
  !----------------------------------------------------------
  !
  !  g~ (k') = sum_R 1/ndegen(R) e^{-ik'R} g (R)
  !
  !  g~(k') is epmatf (nbnd, nbnd, ik )
  !  every pool works with its own subset of k points on the fine grid
  !
  epmatf = czero
  !
  DO ir = 1, nrr
     !
     ! note xk is assumed to be already in cryst coord
     !
     rdotk = twopi * dot_product ( xk, dble(irvec(:, ir)) )
     cfac = exp( ci*rdotk ) / dble( ndegen(ir) )
     !
     DO imode = 1, nmodes
       epmatf (:, :, imode) = epmatf (:, :, imode) + cfac * epmatw ( :, :, ir, imode)
     ENDDO
     !
  ENDDO
  !
  !----------------------------------------------------------
  !  STEP 4: un-rotate to Bloch space, fine grid
  !----------------------------------------------------------
  !
  !  g (k') = U_q^\dagger (k') g~ (k') U_k (k')
  !
  ! the two zgemm calls perform the following ops:
  !  epmatf  = [ cufkq * epmatf ] * cufkk^\dagger
  !
  DO imode = 1, nmodes
    !
    CALL zgemm ('n', 'n', nbnd, nbnd, nbnd, cone, cufkq, &
         nbnd, epmatf (:,:,imode), nbnd, czero, eptmp, nbnd)
    CALL zgemm ('n', 'c', nbnd, nbnd, nbnd, cone, eptmp, &
         nbnd, cufkk, nbnd, czero, epmatf(:,:,imode), nbnd)
    !
  ENDDO
  !
  !
  END SUBROUTINE ephwan2bloch
