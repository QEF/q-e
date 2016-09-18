!
! Copyright (C) 2001-2016 Quantum ESPRESSO group
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
  USE klist,                ONLY : nks, xk, wk, ngk
  USE lsda_mod,             ONLY : nspin
  USE wvfct,                ONLY : npwx,nbnd,wg
  USE gvecw,                ONLY : gcutw
  USE control_flags,        ONLY : gamma_only
  USE gvect,                ONLY : gstart, ngm, g
  USE mp,                   ONLY : mp_sum
  USE mp_global,            ONLY : inter_pool_comm, intra_bgrp_comm
  USE lr_variables,         ONLY : lr_verbosity, eels
  USE noncollin_module,     ONLY : noncolin, npol
  USE control_lr,           ONLY : nbnd_occ
  USE qpoint,               ONLY : nksq
  !
  IMPLICIT NONE
  !
  COMPLEX(kind=dp) :: x(npwx*npol,nbnd,nksq), &
                      y(npwx*npol,nbnd,nksq)
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
       temp_gamma = temp_gamma + 2.D0*wg(ibnd,1)*DDOT(2*ngk(1),x(:,ibnd,1),1,y(:,ibnd,1),1)
       !
       ! G=0 has been accounted twice, so we subtract one contribution.
       !
       IF (gstart==2) temp_gamma = temp_gamma - wg(ibnd,1)*dble(x(1,ibnd,1))*dble(y(1,ibnd,1))
       !
    ENDDO
    !
#if defined(__MPI)
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
          temp_k = temp_k + wg(ibnd,ik) * ZDOTC(ngk(ik),x(1,ibnd,ik),1,y(1,ibnd,ik),1)
          !
       ENDDO
    ENDDO
    !
#if defined(__MPI)
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
    USE qpoint,      ONLY : ikks, ikqs
    !
    IMPLICIT NONE
    INTEGER :: ios
    INTEGER :: ik,  &
               ikk, & ! index of the point k
               ikq, & ! index of the point k+q
               npwq   ! number of the plane-waves at point k+q
    ! 
    DO ik = 1, nksq
       !
       ikk  = ikks(ik)
       ikq  = ikqs(ik)
       npwq = ngk(ikq)
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
#if defined(__MPI)
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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Debugging subroutines
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE check_vector_gamma (x)
   !----------------------------------------------------------------------------
   ! Checks the inner product for a given vector, and its imaginary and real component
   ! input: evc
   ! output : screen output
   !
   USE kinds,                ONLY : dp
   USE mp_global,            ONLY : inter_pool_comm, intra_bgrp_comm
   USE mp,                   ONLY : mp_sum
   USE klist ,               ONLY : ngk
   USE gvect,                ONLY : gstart
   USE io_global,            ONLY : stdout
   !
   IMPLICIT NONE
   !input/output
   COMPLEX(kind=dp),INTENT(in)  :: x(:)
   !
   ! local variables
   !
   REAL(kind=dp) :: temp_gamma
   REAL(kind=dp), EXTERNAL :: DDOT
   !
   temp_gamma = 2.D0*DDOT(2*ngk(1),x(:),1,x(:),1)
   !    
   IF (gstart==2) temp_gamma = temp_gamma - dble(x(1))*dble(x(1))
   ! 
#if defined(__MPI)
   CALL mp_sum(temp_gamma, intra_bgrp_comm)
#endif
   !    
   WRITE(stdout,'("<x> = ",E15.8)') temp_gamma
   !
   RETURN
   !  
END SUBROUTINE check_vector_gamma

SUBROUTINE check_vector_f (x)
   !-----------------------------------------------------------------------
   !
   ! Checks the inner product for a given vector, and its imaginary and real component
   ! input: evc
   ! output: screen output
   ! 
   USE kinds,                ONLY : dp
   USE mp_global,            ONLY : inter_pool_comm, intra_bgrp_comm
   USE mp,                   ONLY : mp_sum
   USE klist ,               ONLY : ngk
   USE gvect,                ONLY : gstart
   USE io_global,            ONLY : stdout
   !
   IMPLICIT NONE
   !input/output
   COMPLEX(kind=dp),INTENT(in)  :: x(:)
   !
   ! local variables
   !
   COMPLEX(kind=dp) :: temp_f
   COMPLEX(kind=dp), EXTERNAL :: ZDOTC
   !
   temp_f = ZDOTC(ngk(1),x(:),1,x(:),1)
   !
#if defined(__MPI)
   CALL mp_sum(temp_f, intra_bgrp_comm)
#endif
   !
   WRITE(stdout,'("<x> = ",2E15.8,1X)') temp_f
   !
   RETURN
   ! 
END SUBROUTINE check_vector_f

SUBROUTINE check_all_bands_gamma (x,sx,nbnd1,nbnd2)
  !----------------------------------------------------------------------
  !
  ! Checks all bands of given KS states for orthoganilty
  ! input: evc and sevc
  ! output : screen output
  !
  USE kinds,                ONLY : dp 
  USE mp_global,            ONLY : inter_pool_comm, intra_bgrp_comm
  USE mp,                   ONLY : mp_sum
  USE klist ,               ONLY : ngk
  USE io_global,            ONLY : stdout
  USE gvect,                ONLY : gstart
  !
  IMPLICIT NONE
  !input/output
  INTEGER, INTENT(in) :: nbnd1,nbnd2 !Total number of bands for x and sx
  COMPLEX(kind=dp),INTENT(in) :: x(:,:), sx(:,:)
  !
  ! local variables
  !
  INTEGER :: ibnd, jbnd
  REAL(kind=dp) :: temp_gamma
  REAL(kind=dp), EXTERNAL :: DDOT
  !
  DO ibnd=1,nbnd1
     DO jbnd=ibnd,nbnd2
        !
        temp_gamma = 2.D0*DDOT(2*ngk(1),x(:,ibnd),1,sx(:,jbnd),1)
        !
        IF (gstart==2) temp_gamma = temp_gamma - dble(x(1,ibnd))*dble(sx(1,jbnd))
        !
#if defined(__MPI)
        CALL mp_sum(temp_gamma, intra_bgrp_comm)
#endif
        !
        WRITE(stdout,'("<x,",I02,"|S|x,",I02,"> =",E15.8)') ibnd,jbnd,temp_gamma
     ENDDO
  ENDDO
  !
  RETURN
  !
END SUBROUTINE check_all_bands_gamma

SUBROUTINE check_density_gamma (rx,nbnd)
  !---------------------------------------------------------------------------------
  !
  ! Checks the contirbution of a given function transformed into real space
  ! input: revc
  ! output : screen output
  !
  USE kinds,                ONLY : dp
  USE mp_global,            ONLY : inter_pool_comm, intra_bgrp_comm
  USE mp,                   ONLY : mp_sum
  USE wvfct,                ONLY : wg
  USE fft_base,             ONLY : dfftp
  USE io_global,            ONLY : stdout
  USE cell_base,            ONLY : omega
  !
  IMPLICIT NONE
  !input/output
  INTEGER, INTENT(in)          :: nbnd !Total number of bands for x and sx
  COMPLEX(kind=dp),INTENT(in)  :: rx(:,:)
  ! 
  ! local variables
  !
  INTEGER :: ibnd
  REAL(kind=dp) :: temp_gamma,w1,w2
  !
  DO ibnd=1,nbnd,2
     w1 = wg(ibnd,1)/omega
     !
     IF (ibnd<nbnd) THEN
        w2 = wg(ibnd+1,1)/omega
     ELSE
        w2 = w1
     ENDIF
     temp_gamma = sum(w1*dble(rx(1:dfftp%nnr,ibnd))*dble(rx(1:dfftp%nnr,ibnd))&
                + w2*aimag(rx(1:dfftp%nnr,ibnd))*aimag(rx(1:dfftp%nnr,ibnd)))
#if defined(__MPI)
     CALL mp_sum(temp_gamma, intra_bgrp_comm)
#endif
     WRITE(stdout,'("Contribution of bands ",I02," and ",I02," to total density",E15.8)') ibnd,ibnd+1,temp_gamma
     !    
  ENDDO
  !
  RETURN
  !
END SUBROUTINE check_density_gamma
