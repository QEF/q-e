!-----------------------------------------------------------------------
!
FUNCTION lr_dot(x,y)
  !---------------------------------------------------------------------
  ! Brent Walker, ICTP, 2004
  !---------------------------------------------------------------------
  ! ... wrapper for PWSCF linear response inner product routines
  ! ... sums over the bands
  ! ... call for each k-point with arguments:
  ! ... call lr_dot(npw_k(ik),evc1(1,1,ik,1),1,evc1(1,1,ik,2),1)
  !---------------------------------------------------------------------
  ! Modified by Osman Baris Malcioglu (2009)
  !
  USE io_global,            ONLY : stdout
  USE kinds,                ONLY : dp
  USE klist,                ONLY : nks
  !use lr_variables,         only : npw_k
  USE realus,               ONLY : npw_k
  USE lsda_mod,             ONLY : nspin
  USE wvfct,                ONLY : npwx,nbnd,wg
  USE control_flags,        ONLY : gamma_only
  USE gvect,                ONLY : gstart
  USE mp,                   ONLY : mp_sum
  USE mp_global,            ONLY : inter_pool_comm, intra_bgrp_comm
  USE lr_variables,   ONLY : lr_verbosity
   USE io_global,      ONLY : stdout
 !
  IMPLICIT NONE
  !
  COMPLEX(kind=dp) :: x(npwx,nbnd,nks),y(npwx,nbnd,nks)
  COMPLEX(kind=dp) :: lr_dot
  COMPLEX(kind=dp) :: temp_k
  real(kind=dp) :: temp_gamma
  real(kind=dp) :: degspin
  INTEGER :: ibnd
  INTEGER :: ik
  real(kind=dp), EXTERNAL    :: DDOT
  COMPLEX(kind=dp), EXTERNAL :: zdotc
  !
  IF (lr_verbosity > 5) THEN
     WRITE(stdout,'("<lr_dot>")')
  ENDIF
  CALL start_clock ('lr_dot')
  !
  lr_dot=(0.0d0,0.0d0)
  !
  temp_gamma=0.0d0
  temp_k=(0.0d0,0.0d0)
  !
  degspin=2.0d0
  IF(nspin==2) degspin=1.0d0
  !
  IF(gamma_only) THEN
     CALL lr_dot_gamma()
     lr_dot=cmplx(temp_gamma,0.0d0,dp)
  ELSE
     CALL lr_dot_k()
     lr_dot=temp_k
  ENDIF
  !
  lr_dot=lr_dot/degspin
  !
  IF (lr_verbosity > 5) WRITE(stdout,'("<end of lr_dot>")')
  CALL stop_clock ('lr_dot')
  !
  RETURN
  !
CONTAINS
  !
  SUBROUTINE lr_dot_gamma
    !
    DO ibnd=1,nbnd
       !
       temp_gamma = temp_gamma + 2.D0*wg(ibnd,1)*DDOT(2*npw_k(1),x(:,ibnd,1),1,y(:,ibnd,1),1)
       IF (gstart==2) temp_gamma = temp_gamma - wg(ibnd,1)*dble(x(1,ibnd,1))*dble(y(1,ibnd,1))
       !
    ENDDO
    !
#ifdef __MPI
    !call reduce(1,temp_gamma)
    CALL mp_sum(temp_gamma, intra_bgrp_comm)
#endif
    !
    RETURN
  END SUBROUTINE lr_dot_gamma
  !
  SUBROUTINE lr_dot_k
    !
    DO ik=1,nks
       DO ibnd=1,nbnd
          !
          temp_k=temp_k+wg(ibnd,ik)*zdotc(npw_k(ik),x(1,ibnd,ik),1,y(1,ibnd,ik),1)
          !
       ENDDO
    ENDDO
#ifdef __MPI
    !call poolreduce(2,temp_k)
    CALL mp_sum(temp_k,inter_pool_comm)
    !call reduce(2,temp_k)
    CALL mp_sum(temp_k, intra_bgrp_comm)
#endif
    !
    RETURN
  END SUBROUTINE lr_dot_k
  !
END FUNCTION lr_dot
!-----------------------------------------------------------------------
