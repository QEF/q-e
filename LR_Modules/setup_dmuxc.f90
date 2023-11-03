!
! Copyright (C) 2001-2016 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
SUBROUTINE setup_dmuxc
  !-----------------------------------------------------------------------
  !! This subroutine computes dmuxc (derivative of the XC potential - LDA case).
  !
  USE kinds,            ONLY : DP
  USE eqv,              ONLY : dmuxc
  USE lsda_mod,         ONLY : lsda
  USE fft_base,         ONLY : dfftp
  USE scf,              ONLY : rho, rho_core
  USE noncollin_module, ONLY : noncolin, domag
  USE xc_lib,           ONLY : dmxc
  !
  IMPLICIT NONE
  !
  REAL(DP), ALLOCATABLE :: rho_aux(:,:)
  ! auxiliary array for density
  INTEGER  :: ir, is, js, ns, dfftp_nnr
  !
  CALL start_clock ('setup_dmuxc')
  !
  dfftp_nnr = dfftp%nnr !to avoid unnecessary copies in acc loop
  !
  !$acc data copyin( rho ) copyout( dmuxc )
  !$acc data copyin( rho%of_r )
  !
  ns = 1
  IF ( lsda ) ns = 2
  IF ( (.NOT. lsda) .AND. noncolin .AND. domag ) ns = 4
  !
  ALLOCATE( rho_aux(dfftp%nnr,ns) )
  !$acc data create( rho_aux )
  !
  IF ( lsda ) THEN
     !
     !$acc parallel loop present(rho) copyin(rho_core)
     DO ir = 1, dfftp_nnr
       rho_aux(ir,1) = ( rho%of_r(ir,1) + rho%of_r(ir,2) + rho_core(ir) )*0.5_DP
       rho_aux(ir,2) = ( rho%of_r(ir,1) - rho%of_r(ir,2) + rho_core(ir) )*0.5_DP
     ENDDO
     !
     CALL dmxc( dfftp_nnr, 2, rho_aux, dmuxc, gpu_args_=.TRUE. )
     !
  ELSE
     !
     IF ( noncolin .AND. domag ) THEN
        !
        !$acc parallel loop present(rho) copyin(rho_core)
        DO ir = 1, dfftp_nnr
          rho_aux(ir,1)   = rho%of_r(ir,1) + rho_core(ir)
          rho_aux(ir,2:4) = rho%of_r(ir,2:4)
        ENDDO
        !
        CALL dmxc( dfftp_nnr, 4, rho_aux, dmuxc, gpu_args_=.TRUE. )
        !
     ELSE
        !
        !$acc parallel loop present(rho) copyin(rho_core)
        DO ir = 1, dfftp_nnr
          rho_aux(ir,1) = rho%of_r(ir,1) + rho_core(ir)
        ENDDO
        !
        CALL dmxc( dfftp_nnr, 1, rho_aux, dmuxc, gpu_args_=.TRUE. )
        !
     ENDIF
     !
  ENDIF
  !
  !$acc end data
  DEALLOCATE( rho_aux )
  !
  !$acc end data
  !$acc end data
  !
  CALL stop_clock ('setup_dmuxc')
  !
  RETURN
  !
END SUBROUTINE setup_dmuxc
