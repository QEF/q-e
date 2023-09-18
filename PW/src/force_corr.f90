!
! Copyright (C) 2001-2007 Quantum ESPRESSO group
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
  !! It uses superposition of atomic charges contained in the array rho_at
  !! and read from pseudopotential files.
  !
  USE kinds,                ONLY : DP
  USE constants,            ONLY : tpi
  USE uspp_param,           ONLY : upf
  USE uspp_data,            ONLY : tab_rho, dq
  USE ions_base,            ONLY : nat, ntyp => nsp, ityp, tau
  USE cell_base,            ONLY : tpiba, omega
  USE fft_base,             ONLY : dfftp
  USE fft_rho,              ONLY : rho_r2g
  USE gvect,                ONLY : ngm, gstart, g, ngl, gl, igtongl
  USE lsda_mod,             ONLY : nspin
  USE scf,                  ONLY : vnew
  USE control_flags,        ONLY : gamma_only
  USE mp_bands,             ONLY : intra_bgrp_comm
  USE mp,                   ONLY : mp_sum
  !
  IMPLICIT NONE
  !
  REAL(DP) :: forcescc(3,nat)
  !
  ! ... local variables
  !
  ! work space
  REAL(DP), ALLOCATABLE :: rhocgnt(:), vauxr(:)
  COMPLEX(DP), ALLOCATABLE :: vauxg(:,:)
  ! temp factors
  REAL(DP) :: gx, arg, fact, forcesccx, forcesccy, forcesccz, &
       ux, px, wx, vx
  ! counters
  INTEGER :: ir, isup, isdw, ig, nt, na, i0,i1,i2,i3
  !
  ! ... vnew is V_out - V_in, psic is the temp space
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
  ALLOCATE( rhocgnt(ngl) )
  !
  !$acc data copyin(tau,gl,ityp) create(rhocgnt,vauxg) present(tab_rho)
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
     !$acc parallel loop 
     DO ig = 1, ngl
        !
        gx = SQRT(gl(ig)) * tpiba
        px = gx / dq - int (gx/dq)
        ux = 1.d0 - px
        vx = 2.d0 - px
        wx = 3.d0 - px
        i0 = INT(gx/dq) + 1
        i1 = i0 + 1
        i2 = i0 + 2
        i3 = i0 + 3
        rhocgnt (ig) = &
                     tab_rho(i0, nt) * ux * vx * wx / 6.d0 + &
                     tab_rho(i1, nt) * px * vx * wx / 2.d0 - &
                     tab_rho(i2, nt) * px * ux * wx / 2.d0 + &
                     tab_rho(i3, nt) * px * ux * vx / 6.d0
     ENDDO
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
                          rhocgnt(igtongl(ig)) * CMPLX(SIN(arg),COS(arg),kind=DP) *  &
                          g(1,ig) * tpiba * CONJG(vauxg(ig,1)))
              forcesccy = forcesccy + DBLE(fact * &
                          rhocgnt(igtongl(ig)) * CMPLX(SIN(arg),COS(arg),kind=DP) *  &
                          g(2,ig) * tpiba * CONJG(vauxg(ig,1)))
              forcesccz = forcesccz + DBLE(fact * &
                          rhocgnt(igtongl(ig)) * CMPLX(SIN(arg),COS(arg),kind=DP) * &
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
  DEALLOCATE( rhocgnt )
  DEALLOCATE( vauxr, vauxg )
  !
  RETURN
  !
END SUBROUTINE force_corr
