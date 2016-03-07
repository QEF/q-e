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
  subroutine ephwan2blochp ( nmodes, xxq, irvec, ndegen, nrr_q, cuf, epmatf, nbnd, nrr_k )
  !---------------------------------------------------------------------------
  !
  ! even though this is for phonons, I use the same notations
  ! adopted for the electronic case (nmodes->nmodes etc)
  !
  !
#include "f_defs.h"
  USE kinds,         only : DP
  USE epwcom,        only : parallel_k, parallel_q
  USE io_epw,        only : iunepmatwp
  use elph2,         only : epmatwp
  USE constants_epw, ONLY : twopi, ci, czero
#ifdef __PARA 
  USE mp_global,     ONLY : inter_pool_comm, intra_pool_comm, mp_sum
#endif
  implicit none
  !
  !  input variables
  !
  integer :: nmodes, nrr_q, irvec ( 3, nrr_q), ndegen (nrr_q), nbnd, nrr_k
  ! number of bands (possibly in tyhe optimal subspace)
  ! number of WS points
  ! coordinates of WS points
  ! degeneracy of WS points
  ! n of bands
  ! n of electronic WS points
  complex(kind=DP) :: epmatw ( nbnd, nbnd, nrr_k, nmodes), cuf (nmodes, nmodes)
  ! e-p matrix in Wanner representation
  ! rotation matrix U(k)
  real(kind=DP) :: xxq(3)
  ! kpoint for the interpolation (WARNING: this must be in crystal coord!)
  !
  !  output variables
  !
  complex(kind=DP) :: epmatf (nbnd, nbnd, nrr_k, nmodes)
  ! e-p matrix in Bloch representation, fine grid
  !
  ! work variables 
  !
  integer :: ibnd, jbnd, ir, ire, ir_start, ir_stop, imode
  real(kind=DP) :: rdotk
  complex(kind=DP) :: eptmp( nbnd, nbnd, nrr_k, nmodes)
  complex(kind=DP) :: cfac(nrr_q)
  !
  CALL start_clock('ephW2Bp')
  !----------------------------------------------------------
  !  STEP 3: inverse Fourier transform of g to fine k mesh
  !----------------------------------------------------------
  !
  !  g~ (k') = sum_R 1/ndegen(R) e^{-ik'R} g (R)
  !
  !  g~(k') is epmatf (nmodes, nmodes, ik )
  !  every pool works with its own subset of k points on the fine grid
  !
  IF (parallel_k) THEN
     CALL para_bounds(ir_start, ir_stop, nrr_q)
  ELSEIF (parallel_q) THEN
     ir_start = 1
     ir_stop  = nrr_q
  ELSE 
     CALL errore ('ephwan2blochp', 'Problem with parallel_k/q scheme', nrr_q)
  ENDIF
  !
  eptmp = czero
  cfac(:) = czero
  !
  DO ir = ir_start, ir_stop
     !   
     ! note xxq is assumed to be already in cryst coord
     !
     rdotk = twopi * dot_product ( xxq, dble(irvec(:, ir)) )
     cfac(ir) = exp( ci*rdotk ) / dble( ndegen(ir) )
  ENDDO
  ! 
  DO ir = ir_start, ir_stop
    eptmp(:,:,:,:) = eptmp(:,:,:,:) +&
      cfac(ir)*epmatwp( :, :, :, :, ir)
  ENDDO
  !
#ifdef __PARA
  IF (parallel_k) CALL mp_sum(eptmp, inter_pool_comm)
#endif
  !
  !----------------------------------------------------------
  !  STEP 4: un-rotate to Bloch space, fine grid
  !----------------------------------------------------------
  !
  ! epmatf(j) = sum_i eptmp(i) * uf(i,j)
  !
  Call zgemm( 'n', 'n', nbnd * nbnd * nrr_k, nmodes, nmodes, ( 1.d0, 0.d0 ),eptmp  , nbnd * nbnd * nrr_k, &
                                                                                    cuf, nmodes         , &
                                                             ( 0.d0, 0.d0 ),epmatf, nbnd * nbnd * nrr_k )
 ! DO ibnd = 1, nbnd
 !  DO jbnd = 1, nbnd
 !    !
 !    CALL zgemm ('n','n',nrr_k, nmodes, nmodes,(1.d0,0.d0),eptmp(ibnd,jbnd,:,:),nrr_k,&
 !        cuf,nmodes,(0.d0,0.d0), epmatf(ibnd,jbnd,:,:), nrr_k )
 !    !
 !  ENDDO
 ! ENDDO
  !
  CALL stop_clock('ephW2Bp')
  !
  end subroutine ephwan2blochp

