!
! Copyright (C) 2001-2023 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
SUBROUTINE force_corr( forcescc )
  !-----------------------------------------------------------------------
  !! This routine calculates the force term vanishing at full
  !! self-consistency. It follows the suggestion of Chan-Bohnen-Ho
  !! (PRB 47, 4771 (1993)). The true charge density is approximated
  !! by means of a free atom superposition.
  !  (alessio f.)
  !
  !! It uses superposition of atomic charges interpolated in tab_rhoat
  !
  USE kinds,                ONLY : DP
  USE constants,            ONLY : tpi
  USE uspp_param,           ONLY : upf
  USE ions_base,            ONLY : nat, ntyp => nsp, ityp, tau
  USE cell_base,            ONLY : tpiba, omega
  USE fft_base,             ONLY : dfftp
  USE fft_rho,              ONLY : rho_r2g
  USE gvect,                ONLY : ngm, gstart, g, ngl, gl, igtongl
  USE lsda_mod,             ONLY : nspin
  USE scf,                  ONLY : vnew
  USE control_flags,        ONLY : gamma_only
  USE mp_bands,             ONLY : intra_bgrp_comm
  USE mp,                   ONLY : mp_sum, mp_max
  USE rhoat_mod,            ONLY : init_tab_rhoat, interp_rhoat
  !  !
  IMPLICIT NONE
  !
  REAL(DP) :: forcescc(3,nat)
  !
  ! ... local variables
  !
  ! work space
  REAL(DP), ALLOCATABLE :: rhoatg(:), vauxr(:)
  COMPLEX(DP), ALLOCATABLE :: vauxg(:,:)
  ! temp factors
  REAL(DP) :: qmax, arg, fact, forcesccx, forcesccy, forcesccz
  ! counters
  INTEGER :: ir, isup, isdw, ig, nt, na, ierr
  !
  !! Initialize the interpolation table if not done previously
  qmax = tpiba * SQRT ( MAXVAL ( gl ) )
  CALL mp_max (qmax, intra_bgrp_comm)
  CALL init_tab_rhoat (qmax, omega, intra_bgrp_comm, ierr)
  !  
  ! ... vnew is V_out - V_in
  !
  ALLOCATE( vauxr(dfftp%nnr), vauxg(dfftp%nnr,1) )
  !
  IF (nspin==1 .OR. nspin==4) THEN
     vauxr(:) = vnew%of_r(:,1)
  ELSE
     isup = 1
     isdw = 2
     vauxr(:) = (vnew%of_r(:,isup) + vnew%of_r(:,isdw)) * 0.5d0
  ENDIF
  !
  ALLOCATE( rhoatg(ngl) )
  !
  !$acc data copyin(tau,gl,ityp) create(rhoatg,vauxg)
  !
  CALL rho_r2g( dfftp, vauxr, vauxg )
  !
  IF (gamma_only) THEN
     fact = 2.d0*omega
  ELSE
     fact = 1.d0*omega
  ENDIF
  !
  DO nt = 1, ntyp
     !
     CALL interp_rhoat( nt, ngl, gl, tpiba**2, rhoatg )
     !
#if defined(_OPENACC)
     !$acc parallel loop gang copy(forcescc)
#else
     !$omp parallel do private(forcesccx,forcesccy,forcesccz,ig,arg)
#endif
     DO na = 1, nat
        IF ( nt == ityp(na) ) THEN
           forcescc(1:3,na) = 0.0_DP
           forcesccx = 0.0d0
           forcesccy = 0.0d0
           forcesccz = 0.0d0
           !$acc loop vector reduction(+:forcesccx,forcesccy,forcesccz)
           DO ig = gstart, ngm
              arg = ( g(1,ig)*tau(1,na) + g(2,ig)*tau(2,na) + g(3,ig)*tau(3,na) ) * tpi
              forcesccx = forcesccx + DBLE(fact * &
                          rhoatg(igtongl(ig)) * CMPLX(SIN(arg),COS(arg),kind=DP) *  &
                          g(1,ig) * tpiba * CONJG(vauxg(ig,1)))
              forcesccy = forcesccy + DBLE(fact * &
                          rhoatg(igtongl(ig)) * CMPLX(SIN(arg),COS(arg),kind=DP) *  &
                          g(2,ig) * tpiba * CONJG(vauxg(ig,1)))
              forcesccz = forcesccz + DBLE(fact * &
                          rhoatg(igtongl(ig)) * CMPLX(SIN(arg),COS(arg),kind=DP) * &
                          g(3,ig) * tpiba * CONJG(vauxg(ig,1)))
           ENDDO
           forcescc(1,na) = forcesccx
           forcescc(2,na) = forcesccy
           forcescc(3,na) = forcesccz
        ENDIF
     ENDDO
  ENDDO
  !
  CALL mp_sum( forcescc, intra_bgrp_comm )
  !
  !$acc end data 
  !
  DEALLOCATE( rhoatg )
  DEALLOCATE( vauxr, vauxg )
  !
  RETURN
  !
END SUBROUTINE force_corr
