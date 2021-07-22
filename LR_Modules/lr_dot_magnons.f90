!
! Copyright (C) 2021 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
FUNCTION lr_dot_magnons(x,y)
  !---------------------------------------------------------------------
  !
  ! Extension of lr_dot.f90 to magnons
  !
  USE kinds,                ONLY : dp
  USE io_global,            ONLY : stdout
  USE klist,                ONLY : nks, xk, wk, ngk
  USE lsda_mod,             ONLY : nspin
  USE wvfct,                ONLY : npwx,nbnd,wg
  USE gvecw,                ONLY : gcutw
  USE control_flags,        ONLY : gamma_only
  USE gvect,                ONLY : gstart, ngm, g
  USE mp,                   ONLY : mp_sum
  USE mp_global,            ONLY : inter_pool_comm, intra_bgrp_comm
  USE noncollin_module,     ONLY : noncolin, npol
  USE control_lr,           ONLY : nbnd_occ, nbnd_occx
  USE qpoint,               ONLY : nksq
  !
  IMPLICIT NONE
  !
  COMPLEX(kind=dp) :: x(npwx*npol,nbnd_occx,nksq,2), &
                      y(npwx*npol,nbnd_occx,nksq,2)
  COMPLEX(kind=dp) :: lr_dot_magnons
  REAL(kind=dp) :: degspin
  INTEGER :: ibnd, ik
  !
  CALL start_clock ('lr_dot_magnons')
  !
  lr_dot_magnons = (0.0d0,0.0d0)
  !
  IF (nspin==2) THEN
      degspin = 1.0d0
  ELSE
      degspin = 2.0d0
  ENDIF
  IF (noncolin) degspin = 1.0d0
  !
  CALL lr_dot_k_magnons()
  !
!  lr_dot_magnons = lr_dot_magnons/degspin
  !
  CALL stop_clock ('lr_dot_magnons')
  !
  RETURN
  !
CONTAINS
  !
  SUBROUTINE lr_dot_k_magnons
    !
    ! MAGNONS
    ! Noncollinear case is implemented
    !
    USE qpoint,      ONLY : ikks, ikqs
    !
    IMPLICIT NONE
    INTEGER :: ios
    INTEGER :: ik,  &
               ikk, & ! index of the point k
               ikq, & ! index of the point k+q
               npwq, &! number of the plane-waves at point k+q
               imk
    ! 
    DO ik = 1, nksq
       !
       ikk  = ikks(ik)
       ikq  = ikqs(ik)
       npwq = ngk(ikq)
       !
       IF ( mod(ik,2) == 0) THEN
          imk  = ikk - 3      ! position of -k
       ELSE
          imk  = ikk + 3      ! position of -k
       ENDIF
       !
       ! Resonant part (upper batch)
       !
       DO ibnd = 1, nbnd_occ(ikk)
          !
          lr_dot_magnons = lr_dot_magnons + wk(ikk) * &
                  dot_product(x(:,ibnd,ik,1),y(:,ibnd,ik,1))
          !
       ENDDO
       !
       ! Anti - Resonant part (lower batch)
       !
       DO ibnd = 1, nbnd_occ(imk)
          !
          lr_dot_magnons = lr_dot_magnons + wk(imk) * &
                  dot_product(x(:,ibnd,ik,2),y(:,ibnd,ik,2))
          !
       ENDDO
       !
    ENDDO
    !
#if defined(__MPI)
    CALL mp_sum(lr_dot_magnons, intra_bgrp_comm)
    CALL mp_sum(lr_dot_magnons, inter_pool_comm)
#endif
    !
    RETURN
    !
  END SUBROUTINE lr_dot_k_magnons
  !
END FUNCTION lr_dot_magnons
