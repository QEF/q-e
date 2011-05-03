!
! Copyright (C) 2001-2007 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

!-----------------------------------------------------------------------
SUBROUTINE lr_ch_psi_all (n, h, ah, e, ik, m)
  !-----------------------------------------------------------------------
  !
  ! This routine applies the operator ( H - \epsilon S + alpha_pv P_v)
  ! to a vector h. The result is given in Ah.
  !
  ! Modified by Osman Baris Malcioglu (2009)
  !
#include "f_defs.h"

  USE kinds, ONLY : DP
  USE wvfct, ONLY : npwx, nbnd
  USE uspp, ONLY: nkb, vkb
  USE noncollin_module, ONLY : noncolin, npol

  USE control_ph, ONLY : alpha_pv, nbnd_occ
  !USE eqv,        ONLY : evq
  !USE qpoint,     ONLY : ikqs !directly use ik instead

  USE mp_global, ONLY: intra_pool_comm
  USE mp,        ONLY: mp_sum

  USE control_flags,                ONLY : gamma_only
  USE wavefunctions_module, ONLY : evc !evq is replaced by evc
  !use lr_variables,         only : lr_alpha_pv, nbnd_occ,
  USE lr_variables,         ONLY : lr_verbosity
  USE io_global,            ONLY : stdout
  IMPLICIT NONE

  INTEGER :: n, m, ik
  ! input: the dimension of h
  ! input: the number of bands
  ! input: the k point

  real(DP) :: e (m)
  ! input: the eigenvalue

  COMPLEX(DP) :: h (npwx*npol, m), ah (npwx*npol, m)
  ! input: the vector
  ! output: the operator applied to the vector
  !
  !   local variables
  !
  INTEGER :: ibnd, ikq, ig
  ! counter on bands
  ! the point k+q
  ! counter on G vetors

  COMPLEX(DP), ALLOCATABLE :: hpsi (:,:), spsi (:,:)
  ! scalar products
  ! the product of the Hamiltonian and h
  ! the product of the S matrix and h
  !OBM debug
  !real(DP) :: obm_debug
  !complex(kind=dp), external :: ZDOTC

  CALL start_clock ('ch_psi')

  IF (lr_verbosity > 5) WRITE(stdout,'("<lr_ch_psi_all>")')

  ALLOCATE (hpsi( npwx*npol , m))
  ALLOCATE (spsi( npwx*npol , m))
  hpsi (:,:) = (0.d0, 0.d0)
  spsi (:,:) = (0.d0, 0.d0)
  !
  !   compute the product of the hamiltonian with the h vector
  !!OBM debug
  !    obm_debug=0
  !    do ibnd=1,m
  !       !
  !       obm_debug=obm_debug+ZDOTC(npwx*npol,h(:,ibnd),1,h(:,ibnd),1)
  !       !
  !    enddo
  !    print *, "lr_ch_psi_all h", obm_debug
  !!obm_debug
 !
  CALL lr_h_psiq (npwx, n, m, h, hpsi, spsi)
  !!OBM debug
  !    obm_debug=0
  !    do ibnd=1,m
  !       !
  !       obm_debug=obm_debug+ZDOTC(npwx*npol,hpsi(:,ibnd),1,hpsi(:,ibnd),1)
  !       !
  !    enddo
  !    print *, "lr_ch_psi_all hpsi", obm_debug
  !!obm_debug

  CALL start_clock ('last')
  !
  !   then we compute the operator H-epsilon S
  !
  ah=(0.d0,0.d0)
  DO ibnd = 1, m
     DO ig = 1, n
        ah (ig, ibnd) = hpsi (ig, ibnd) - e (ibnd) * spsi (ig, ibnd)
     ENDDO
  ENDDO
  IF (noncolin) THEN
     DO ibnd = 1, m
        DO ig = 1, n
           ah (ig+npwx,ibnd)=hpsi(ig+npwx,ibnd)-e(ibnd)*spsi(ig+npwx,ibnd)
        ENDDO
     ENDDO
  ENDIF
  !
  !   Here we compute the projector in the valence band
  !
  IF(gamma_only) THEN
       CALL lr_ch_psi_all_gamma()
    ELSE
       CALL lr_ch_psi_all_k()
  ENDIF
  !!OBM debug
  !    obm_debug=0
  !    do ibnd=1,m
  !       !
  !       obm_debug=obm_debug+ZDOTC(npwx*npol,ah(:,ibnd),1,ah(:,ibnd),1)
  !       !
  !    enddo
  !    print *, "lr_ch_psi_all ah", obm_debug
  !!obm_debug
 !
  DEALLOCATE (spsi)
  DEALLOCATE (hpsi)
  CALL stop_clock ('last')
  CALL stop_clock ('ch_psi')
  RETURN
CONTAINS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!K-point part
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE lr_ch_psi_all_k()

     USE becmod, ONLY : becp, calbec

     IMPLICIT NONE

     COMPLEX(kind=dp), ALLOCATABLE :: ps(:,:)

     ALLOCATE (ps  ( nbnd , m))
     !ikq = ikqs(ik)
     ps (:,:) = (0.d0, 0.d0)

     IF (noncolin) THEN
        CALL ZGEMM ('C', 'N', nbnd_occ (ik) , m, npwx*npol, (1.d0, 0.d0) , evc, &
          npwx*npol, spsi, npwx*npol, (0.d0, 0.d0) , ps, nbnd)
     ELSE
        CALL ZGEMM ('C', 'N', nbnd_occ (ik) , m, n, (1.d0, 0.d0) , evc, &
          npwx, spsi, npwx, (0.d0, 0.d0) , ps, nbnd)
     ENDIF
     ps (:,:) = ps(:,:) * alpha_pv
#ifdef __PARA
     CALL mp_sum ( ps, intra_pool_comm )
#endif

     hpsi (:,:) = (0.d0, 0.d0)
     IF (noncolin) THEN
        CALL ZGEMM ('N', 'N', npwx*npol, m, nbnd_occ (ik) , (1.d0, 0.d0) , evc, &
             npwx*npol, ps, nbnd, (1.d0, 0.d0) , hpsi, npwx*npol)
     ELSE
        CALL ZGEMM ('N', 'N', n, m, nbnd_occ (ik) , (1.d0, 0.d0) , evc, &
             npwx, ps, nbnd, (1.d0, 0.d0) , hpsi, npwx)
     ENDIF
     spsi(:,:) = hpsi(:,:)
     !
     !    And apply S again
     !
     CALL calbec (n, vkb, hpsi, becp, m)
     CALL s_psi (npwx, n, m, hpsi, spsi)
     DO ibnd = 1, m
        DO ig = 1, n
           ah (ig, ibnd) = ah (ig, ibnd) + spsi (ig, ibnd)
        ENDDO
     ENDDO
     IF (noncolin) THEN
        DO ibnd = 1, m
           DO ig = 1, n
              ah (ig+npwx, ibnd) = ah (ig+npwx, ibnd) + spsi (ig+npwx, ibnd)
           ENDDO
        ENDDO
     ENDIF
     DEALLOCATE (ps)
  END SUBROUTINE lr_ch_psi_all_k
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!gamma part
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE lr_ch_psi_all_gamma()

     USE becmod, ONLY : becp,  calbec
     USE realus, ONLY : real_space, fft_orbital_gamma, &
                                    bfft_orbital_gamma, calbec_rs_gamma, &
                                    s_psir_gamma,real_space_debug

     IMPLICIT NONE

     real(kind=dp), ALLOCATABLE :: ps(:,:)

     ALLOCATE (ps  ( nbnd , m))
     !ikq = ikqs(ik)
     ps (:,:) = 0.d0

     IF (noncolin) THEN
       CALL errore('lr_ch_psi_all', 'non collin in gamma point not implemented',1)
     ELSE
       CALL DGEMM( 'C', 'N', nbnd, m, n, 2.D0,evc, 2*npwx*npol, spsi, 2*npwx*npol, 0.D0, ps, nbnd )
     ENDIF
     ps (:,:) = ps(:,:) * alpha_pv
#ifdef __PARA
     CALL mp_sum ( ps, intra_pool_comm )
#endif

     hpsi (:,:) = (0.d0, 0.d0)
     IF (noncolin) THEN
        CALL ZGEMM ('N', 'N', npwx*npol, m, nbnd_occ (ik) , (1.d0, 0.d0) , evc, &
             npwx*npol, ps, nbnd, (1.d0, 0.d0) , hpsi, npwx*npol)
     ELSE
        CALL DGEMM ('N', 'N', 2*n, m, nbnd_occ (ik) , 1.d0 , evc, &
             2*npwx, ps, nbnd, 1.d0 , hpsi, 2*npwx)
     ENDIF
     spsi(:,:) = hpsi(:,:)
     !
     !    And apply S again
     !
   IF (real_space_debug >6 ) THEN
     DO ibnd=1,m,2
      CALL fft_orbital_gamma(hpsi,ibnd,m)
      CALL calbec_rs_gamma(ibnd,m,becp%r)
      CALL s_psir_gamma(ibnd,m)
      CALL bfft_orbital_gamma(spsi,ibnd,m)
     ENDDO
    ELSE
     CALL calbec (n, vkb, hpsi, becp, m)
     CALL s_psi (npwx, n, m, hpsi, spsi)
    ENDIF
     DO ibnd = 1, m
        DO ig = 1, n
           ah (ig, ibnd) = ah (ig, ibnd) + spsi (ig, ibnd)
        ENDDO
     ENDDO
     IF (noncolin) THEN
        DO ibnd = 1, m
           DO ig = 1, n
              ah (ig+npwx, ibnd) = ah (ig+npwx, ibnd) + spsi (ig+npwx, ibnd)
           ENDDO
        ENDDO
     ENDIF
     DEALLOCATE (ps)
  END SUBROUTINE lr_ch_psi_all_gamma

END SUBROUTINE lr_ch_psi_all
