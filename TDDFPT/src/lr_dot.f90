!
! Copyright (C) 2001-2015 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
FUNCTION lr_dot(x,y)
  !---------------------------------------------------------------------
  !
  ! This subroutine calculates a dot product of the conjugate 
  ! of a complex vector x and a complex vector y 
  ! (sums over the bands and k-points).
  !
  ! Brent Walker, ICTP, 2004
  ! Modified by Osman Baris Malcioglu, SISSA, 2009
  ! Modified by Iurii Timrov, SISSA, 2013 (extension to EELS)
  !
  USE kinds,                ONLY : dp
  USE io_global,            ONLY : stdout
  USE klist,                ONLY : nks, xk, wk
  USE realus,               ONLY : npw_k
  USE lsda_mod,             ONLY : nspin
  USE wvfct,                ONLY : npwx,nbnd,wg,npw,igk,ecutwfc,g2kin
  USE control_flags,        ONLY : gamma_only
  USE gvect,                ONLY : gstart, ngm, g
  USE cell_base,            ONLY : tpiba2
  USE mp,                   ONLY : mp_sum
  USE mp_global,            ONLY : inter_pool_comm, intra_bgrp_comm
  USE lr_variables,         ONLY : lr_verbosity, lr_periodic, eels
  USE noncollin_module,     ONLY : noncolin, npol
  USE control_ph,           ONLY : nbnd_occ
  USE qpoint,               ONLY : npwq, igkq, ikks, ikqs, nksq
  !
  IMPLICIT NONE
  !
  COMPLEX(kind=dp) :: x(npwx*npol,nbnd,nksq), y(npwx*npol,nbnd,nksq)
  COMPLEX(kind=dp) :: lr_dot
  COMPLEX(kind=dp) :: temp_k
  REAL(kind=dp) :: temp_gamma, degspin
  INTEGER :: ibnd, ik
  REAL(kind=dp), EXTERNAL    :: DDOT
  COMPLEX(kind=dp), EXTERNAL :: ZDOTC
  !
  IF (lr_verbosity > 5) THEN
     WRITE(stdout,'("<lr_dot>")')
  ENDIF
  !
  CALL start_clock ('lr_dot')
  !
  lr_dot = (0.0d0,0.0d0)
  temp_gamma = 0.0d0
  temp_k = (0.0d0,0.0d0)
  !
  IF (nspin==2) THEN
      degspin = 1.0d0
  ELSE
      degspin = 2.0d0
  ENDIF
  IF (noncolin) degspin = 1.0d0
  !
  IF (eels) THEN
     !
     CALL lr_dot_k_eels()
     !
  ELSE
     !
     IF (gamma_only) THEN
        !
        CALL lr_dot_gamma()
        lr_dot = cmplx(temp_gamma,0.0d0,dp)
        !
     ELSE
        !
        CALL lr_dot_k()
        lr_dot = temp_k
        !
     ENDIF
     !
  ENDIF
  !
  lr_dot = lr_dot/degspin
  !
  IF (lr_verbosity > 5) WRITE(stdout,'("<end of lr_dot>")')
  !
  CALL stop_clock ('lr_dot')
  !
  RETURN
  !
CONTAINS
  !
  SUBROUTINE lr_dot_gamma
    !
    ! Optical case: gamma_only
    ! Noncollinear case is not implemented
    !
    DO ibnd=1,nbnd
       !
       temp_gamma = temp_gamma + 2.D0*wg(ibnd,1)*DDOT(2*npw_k(1),x(:,ibnd,1),1,y(:,ibnd,1),1)
       !
       ! G=0 has been accounted twice, so we subtract one contribution.
       !
       IF (gstart==2) temp_gamma = temp_gamma - wg(ibnd,1)*dble(x(1,ibnd,1))*dble(y(1,ibnd,1))
       !
    ENDDO
    !
#ifdef __MPI
    CALL mp_sum(temp_gamma, intra_bgrp_comm)
#endif
    !
    RETURN
    !
  END SUBROUTINE lr_dot_gamma
  !
  SUBROUTINE lr_dot_k
    !
    ! Optical case: generalized k point case 
    ! Noncollinear case is not implemented
    !
    DO ik=1,nks   
       DO ibnd=1,nbnd
          !
          temp_k = temp_k + wg(ibnd,ik) * ZDOTC(npw_k(ik),x(1,ibnd,ik),1,y(1,ibnd,ik),1)
          !
       ENDDO
    ENDDO
    !
#ifdef __MPI
    CALL mp_sum(temp_k, inter_pool_comm)
    CALL mp_sum(temp_k, intra_bgrp_comm)
#endif
    !
    RETURN
    !
  END SUBROUTINE lr_dot_k
  !
  SUBROUTINE lr_dot_k_eels
    !
    ! EELS
    ! Noncollinear case is implemented
    !
    IMPLICIT NONE
    INTEGER :: ios, ikk, ikq
    !
    !IF (nksq > 1) rewind (unit = iunigk)
    !
    DO ik = 1, nksq
       !
       IF (lr_periodic) THEN
          ikk = ik
          ikq = ik
       ELSE
          ikk = ikks(ik)
          ikq = ikqs(ik)
       ENDIF
       !
       ! Determination of the number of plane waves npwq at point k+q (ikq).
       !
       CALL gk_sort( xk(1,ikq), ngm, g, ( ecutwfc / tpiba2 ), npwq, igkq, g2kin)
       !
!      IF (nksq > 1) THEN
!         read (iunigk, err = 100, iostat = ios) npw, igk
!100      call errore ('lr_dot', 'reading igk', abs (ios) )
!         read (iunigk, err = 200, iostat = ios) npwq, igkq
!200      call errore ('lr_dot', 'reading igkq', abs (ios) )
!      ENDIF
       !
       DO ibnd = 1, nbnd_occ(ikk)
          !
          IF (noncolin) THEN
             lr_dot = lr_dot + wk(ikk) * ZDOTC(npwx*npol,x(1,ibnd,ik),1,y(1,ibnd,ik),1)
          ELSE
             lr_dot = lr_dot + wk(ikk) * ZDOTC(npwq,x(1,ibnd,ik),1,y(1,ibnd,ik),1)
          ENDIF
          !
       ENDDO
       !
    ENDDO
    !
#ifdef __MPI
    CALL mp_sum(lr_dot, inter_pool_comm)
    CALL mp_sum(lr_dot, intra_bgrp_comm)
#endif
    !
    RETURN
    !
  END SUBROUTINE lr_dot_k_eels
  !
END FUNCTION lr_dot
!-----------------------------------------------------------------------
