!
! Copyright (C) 2001-2007 Quantum-ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

!-----------------------------------------------------------------------
subroutine lr_ch_psi_all (n, h, ah, e, ik, m)
  !-----------------------------------------------------------------------
  !
  ! This routine applies the operator ( H - \epsilon S + alpha_pv P_v)
  ! to a vector h. The result is given in Ah.
  !
  !OBM!!! 19jul2009 modified from phonon version according to Brent walker's guideline
  !
#include "f_defs.h"

  USE kinds, only : DP
  USE wvfct, ONLY : npwx, nbnd
  USE uspp, ONLY: nkb, vkb
  USE noncollin_module, ONLY : noncolin, npol

  USE control_ph, ONLY : alpha_pv, nbnd_occ
  !USE eqv,        ONLY : evq
  !USE qpoint,     ONLY : ikqs !directly use ik instead

  USE mp_global, ONLY: intra_pool_comm
  USE mp,        ONLY: mp_sum

  use control_flags,                only : gamma_only
  use wavefunctions_module, only : evc !evq is replaced by evc
  !use lr_variables,         only : lr_alpha_pv, nbnd_occ, 
  use lr_variables,         only : lr_verbosity
  use io_global,            only : stdout
  implicit none

  integer :: n, m, ik
  ! input: the dimension of h
  ! input: the number of bands
  ! input: the k point

  real(DP) :: e (m)
  ! input: the eigenvalue

  complex(DP) :: h (npwx*npol, m), ah (npwx*npol, m)
  ! input: the vector
  ! output: the operator applied to the vector
  !
  !   local variables
  !
  integer :: ibnd, ikq, ig
  ! counter on bands
  ! the point k+q
  ! counter on G vetors

  complex(DP), allocatable :: hpsi (:,:), spsi (:,:)
  ! scalar products
  ! the product of the Hamiltonian and h
  ! the product of the S matrix and h
  !obm debug
  !real(DP) :: obm_debug
  !complex(kind=dp), external :: ZDOTC
  
  call start_clock ('ch_psi')
  
  If (lr_verbosity > 5) WRITE(stdout,'("<lr_ch_psi_all>")')
  
  allocate (hpsi( npwx*npol , m))    
  allocate (spsi( npwx*npol , m))    
  hpsi (:,:) = (0.d0, 0.d0)
  spsi (:,:) = (0.d0, 0.d0)
  !
  !   compute the product of the hamiltonian with the h vector
  !!OBM debug
  !     obm_debug=0
  !     do ibnd=1,m
  !        !
  !        obm_debug=obm_debug+ZDOTC(npwx*npol,h(:,ibnd),1,h(:,ibnd),1)
  !        !
  !     enddo
  !     print *, "lr_ch_psi_all h", obm_debug
  !!obm_debug
 !
  call lr_h_psiq (npwx, n, m, h, hpsi, spsi)
  !!OBM debug
  !     obm_debug=0
  !     do ibnd=1,m
  !        !
  !        obm_debug=obm_debug+ZDOTC(npwx*npol,hpsi(:,ibnd),1,hpsi(:,ibnd),1)
  !        !
  !     enddo
  !     print *, "lr_ch_psi_all hpsi", obm_debug
  !!obm_debug

  call start_clock ('last')
  !
  !   then we compute the operator H-epsilon S
  !
  ah=(0.d0,0.d0)
  do ibnd = 1, m
     do ig = 1, n
        ah (ig, ibnd) = hpsi (ig, ibnd) - e (ibnd) * spsi (ig, ibnd)
     enddo
  enddo
  IF (noncolin) THEN
     do ibnd = 1, m
        do ig = 1, n
           ah (ig+npwx,ibnd)=hpsi(ig+npwx,ibnd)-e(ibnd)*spsi(ig+npwx,ibnd)
        end do
     end do
  END IF
  !
  !   Here we compute the projector in the valence band
  !
  if(gamma_only) then                     
       call lr_ch_psi_all_gamma()        
    else                                      
       call lr_ch_psi_all_k()                 
  endif                                    
  !!OBM debug
  !     obm_debug=0
  !     do ibnd=1,m
  !        !
  !        obm_debug=obm_debug+ZDOTC(npwx*npol,ah(:,ibnd),1,ah(:,ibnd),1)
  !        !
  !     enddo
  !     print *, "lr_ch_psi_all ah", obm_debug
  !!obm_debug
 !              
  deallocate (spsi)
  deallocate (hpsi)
  call stop_clock ('last')
  call stop_clock ('ch_psi')
  return
contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!K-point part
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine lr_ch_psi_all_k()
    
     USE becmod, ONLY : becp, becp_nc, calbec

     IMPLICIT NONE
 
     complex(kind=dp), allocatable :: ps(:,:)

     allocate (ps  ( nbnd , m))    
     !ikq = ikqs(ik)
     ps (:,:) = (0.d0, 0.d0)
    
     IF (noncolin) THEN
        call ZGEMM ('C', 'N', nbnd_occ (ik) , m, npwx*npol, (1.d0, 0.d0) , evc, &
          npwx*npol, spsi, npwx*npol, (0.d0, 0.d0) , ps, nbnd)
     ELSE
        call ZGEMM ('C', 'N', nbnd_occ (ik) , m, n, (1.d0, 0.d0) , evc, &
          npwx, spsi, npwx, (0.d0, 0.d0) , ps, nbnd)
     ENDIF
     ps (:,:) = ps(:,:) * alpha_pv
#ifdef __PARA
     call mp_sum ( ps, intra_pool_comm )
#endif

     hpsi (:,:) = (0.d0, 0.d0)
     IF (noncolin) THEN
        call ZGEMM ('N', 'N', npwx*npol, m, nbnd_occ (ik) , (1.d0, 0.d0) , evc, &
             npwx*npol, ps, nbnd, (1.d0, 0.d0) , hpsi, npwx*npol)
     ELSE
        call ZGEMM ('N', 'N', n, m, nbnd_occ (ik) , (1.d0, 0.d0) , evc, &
             npwx, ps, nbnd, (1.d0, 0.d0) , hpsi, npwx)
     END IF
     spsi(:,:) = hpsi(:,:)
     !
     !    And apply S again
     !
     IF (noncolin) THEN
        call calbec (n, vkb, hpsi, becp_nc, m)
     ELSE
        call calbec (n, vkb, hpsi, becp, m)
     END IF
     call s_psi (npwx, n, m, hpsi, spsi)
     do ibnd = 1, m
        do ig = 1, n
           ah (ig, ibnd) = ah (ig, ibnd) + spsi (ig, ibnd)
        enddo
     enddo
     IF (noncolin) THEN
        do ibnd = 1, m
           do ig = 1, n
              ah (ig+npwx, ibnd) = ah (ig+npwx, ibnd) + spsi (ig+npwx, ibnd)
           enddo
        enddo
     END IF
     deallocate (ps)
  end subroutine lr_ch_psi_all_k
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!gamma part
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine lr_ch_psi_all_gamma()
    
     USE becmod, ONLY : rbecp,  calbec

     IMPLICIT NONE
 
     real(kind=dp), allocatable :: ps(:,:)

     allocate (ps  ( nbnd , m))    
     !ikq = ikqs(ik)
     ps (:,:) = 0.d0
    
     IF (noncolin) THEN
       call errore('lr_ch_psi_all', 'non collin in gamma point not implemented',1)
     ELSE
       CALL DGEMM( 'C', 'N', nbnd, m, n, 2.D0,evc, 2*npwx*npol, spsi, 2*npwx*npol, 0.D0, ps, nbnd )
     ENDIF
     ps (:,:) = ps(:,:) * alpha_pv
#ifdef __PARA
     call mp_sum ( ps, intra_pool_comm )
#endif

     hpsi (:,:) = (0.d0, 0.d0)
     IF (noncolin) THEN
        call ZGEMM ('N', 'N', npwx*npol, m, nbnd_occ (ik) , (1.d0, 0.d0) , evc, &
             npwx*npol, ps, nbnd, (1.d0, 0.d0) , hpsi, npwx*npol)
     ELSE
        call DGEMM ('N', 'N', 2*n, m, nbnd_occ (ik) , 1.d0 , evc, &
             2*npwx, ps, nbnd, 1.d0 , hpsi, 2*npwx)
     END IF
     spsi(:,:) = hpsi(:,:)
     !
     !    And apply S again
     !
     IF (noncolin) THEN
        !call calbec (n, vkb, hpsi, becp_nc, m)
     ELSE
        call calbec (n, vkb, hpsi, rbecp, m)
     END IF
     call s_psi (npwx, n, m, hpsi, spsi)
     do ibnd = 1, m
        do ig = 1, n
           ah (ig, ibnd) = ah (ig, ibnd) + spsi (ig, ibnd)
        enddo
     enddo
     IF (noncolin) THEN
        do ibnd = 1, m
           do ig = 1, n
              ah (ig+npwx, ibnd) = ah (ig+npwx, ibnd) + spsi (ig+npwx, ibnd)
           enddo
        enddo
     END IF
     deallocate (ps)
  end subroutine lr_ch_psi_all_gamma

end subroutine lr_ch_psi_all
