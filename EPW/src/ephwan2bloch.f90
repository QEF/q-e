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
  !
#include "f_defs.h"
  USE kinds,         ONLY : DP
  USE constants_epw, ONLY : twopi, ci, czero, cone
  implicit none
  !
  !  input variables
  !
  integer :: nbnd, nrr, irvec ( 3, nrr), ndegen (nrr), nmodes
  ! number of bands (possibly in tyhe optimal subspace)
  ! number of WS points
  ! coordinates of WS points
  ! degeneracy of WS points
  ! number of phonon modes
  complex(kind=DP) :: epmatw ( nbnd, nbnd, nrr, nmodes), cufkk (nbnd, nbnd), &
    cufkq (nbnd, nbnd)
  ! e-p matrix in Wannier representation
  ! rotation matrix U(k)
  ! rotation matrix U(k+q)
  real(kind=DP) :: xk(3)
  ! kpoint for the interpolation (WARNING: this must be in crystal coord!)
  !
  !  output variables
  !
  complex(kind=DP) :: epmatf (nbnd, nbnd, nmodes)
  ! e-p matrix in Bloch representation, fine grid
  !
  ! work variables 
  !
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
