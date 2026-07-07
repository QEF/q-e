!
! Copyright (C) 2001-2025 Quantum ESPRESSO Foundation
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
SUBROUTINE vhpsi_U( ldap, np, mps, psip, hpsi )
  !-----------------------------------------------------------------------
  !! This routine computes the Hubbard potential applied to the electronic
  !! structure of the current k-point. The result is added to hpsi.
  !! GPU version with OpenACC
  !
  USE kinds,         ONLY : DP
  USE ldaU,          ONLY : Hubbard_lmax, Hubbard_l, is_hubbard,   &
                            nwfcU, wfcU, offsetU, lda_plus_u_kind, &
                            is_hubbard_back, offsetU_back, backall, &
                            offsetU_back1, ldim_back, ldmx_b, &
                            Hubbard_l2, Hubbard_l3
  USE lsda_mod,      ONLY : current_spin
  USE scf,           ONLY : v
  USE ions_base,     ONLY : nat, ntyp => nsp, ityp
  USE control_flags, ONLY : gamma_only, offload_type
  USE becmod,        ONLY : calbec
  USE noncollin_module, ONLY: noncolin, npol
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: ldap
  !! leading dimension of arrays psip, hpsi
  INTEGER, INTENT(IN) :: np
  !! true dimension of psip, hpsi
  INTEGER, INTENT(IN) :: mps
  !! number of states psip
  COMPLEX(DP), INTENT(IN) :: psip(npol*ldap,mps)
  !! the wavefunction
  COMPLEX(DP), INTENT(INOUT) :: hpsi(npol*ldap,mps)
  !! Hamiltonian dot psi
  !
  IF ( .NOT. ANY(is_hubbard(:)) .AND. .NOT.ANY(is_hubbard_back(:)) ) RETURN
  !
  !$acc data present(wfcU)
  !
  IF (gamma_only) THEN
     CALL vhpsi_gamma_acc ()
  ELSE IF ( noncolin ) THEN
     CALL vhpsi_nc_acc( ldap, np, mps, psip, hpsi )
  ELSE
     CALL vhpsi_k_acc ()
  ENDIF
  !
  !$acc end data
  !
  RETURN
  !
CONTAINS
  !
SUBROUTINE vhpsi_gamma_acc()
  !
  ! Gamma-only version
  !
  IMPLICIT NONE
  !
  ! ... local variables
  !
  REAL(DP), ALLOCATABLE :: proj_r(:,:)
  REAL(DP), ALLOCATABLE :: rtemp(:,:), vns_r(:,:,:),  vnsb_r(:,:,:)
  !
  INTEGER :: na, nt, ldim, ldim0, ldimax, ldimaxt
  !
  ALLOCATE( proj_r(nwfcU, mps) )
  !$acc enter data create(proj_r)
  CALL calbec(offload_type, np, wfcU, psip, proj_r)
  !
  ldimax = 2*Hubbard_lmax+1
  ldimaxt = MAX(ldimax, ldmx_b)
  ALLOCATE( rtemp(ldimaxt,mps) )
  !
  !$acc enter data create(rtemp)
  IF (ANY(is_hubbard(:))) THEN
     ALLOCATE( vns_r(ldimax,ldimax,nat) )
     !$acc enter data create(vns_r)
     vns_r = v%ns(:,:,current_spin,:)
     !$acc update device(vns_r)
  ENDIF
  IF (ANY(is_hubbard_back(:))) THEN
     ALLOCATE( vnsb_r(ldmx_b,ldmx_b,nat) )
     !$acc enter data create(vnsb_r)
     vnsb_r = v%nsb(:,:,current_spin,:)
     !$acc update device(vnsb_r)
  ENDIF
  !
  DO nt = 1, ntyp
     !
     ! Compute the action of the Hubbard potential on the KS wave functions:
     ! V_Hub |psip > = \sum v%ns |wfcU> <wfcU|psip>
     ! where v%ns = U ( delta/2 - rho%ns ) is computed in v_of_rho
     !
     IF ( is_hubbard(nt) ) THEN
        !  
        ldim = 2*Hubbard_l(nt) + 1
        !
        DO na = 1, nat
           IF ( nt == ityp(na) ) THEN
              !
              !$acc host_data use_device(proj_r,vns_r,rtemp)
              CALL MYDGEMM( 'N','N', ldim,mps,ldim, 1.0_dp, &
                   vns_r(1,1,na), ldimax, &
                   proj_r(offsetU(na)+1,1), nwfcU, 0.0_dp, rtemp, ldimaxt )
              !$acc end host_data
              !
              !$acc host_data use_device(wfcU,rtemp,hpsi)
              CALL MYDGEMM( 'N','N', 2*np, mps, ldim, 1.0_dp, &
                   wfcU(1,offsetU(na)+1), 2*ldap, rtemp, ldimaxt, &
                   1.0_dp, hpsi, 2*ldap )
              !$acc end host_data
              !
           ENDIF
        ENDDO
        !
     ENDIF
     !
     ! If the background is used then compute extra 
     ! contribution to the Hubbard potential
     !
     IF ( is_hubbard_back(nt) ) THEN
        !
        ldim = ldim_back(nt)
        !
        DO na = 1, nat
           IF ( nt == ityp(na) ) THEN
              !
              ldim = 2*Hubbard_l2(nt)+1
              !
              !$acc host_data use_device(proj_r,vnsb_r,rtemp)
              CALL MYDGEMM( 'N','N', ldim,mps,ldim, 1.0_dp, &
                   vnsb_r(1,1,na),ldmx_b, &
                   proj_r(offsetU_back(na)+1,1), &
                   nwfcU, 0.0_dp, rtemp, ldimaxt )
              !$acc end host_data
              !
              !$acc host_data use_device(wfcU, rtemp,hpsi)
              CALL MYDGEMM( 'N','N', 2*np, mps, ldim, 1.0_dp, &
                   wfcU(1,offsetU_back(na)+1), 2*ldap, rtemp, &
                   ldimaxt, 1.0_dp, hpsi, 2*ldap )
              !$acc end host_data
              !
              IF (backall(nt)) THEN
                 !
                 ldim0 = 2*Hubbard_l2(nt)+1
                 ldim  = 2*Hubbard_l3(nt)+1
                 !
                 !$acc host_data use_device(proj_r,vnsb_r,rtemp)
                 CALL MYDGEMM( 'N', 'N', ldim,mps,ldim, 1.0_dp,     &
                      vnsb_r(ldim0+1,ldim0+1,na),                       &
                      ldim_back(nt), proj_r(offsetU_back1(na)+1,1), &
                      nwfcU, 0.0_dp, rtemp, ldimaxt )
                 !$acc end host_data
                 !
                 !$acc host_data use_device(wfcU, rtemp,hpsi)
                 CALL MYDGEMM( 'N', 'N', 2*np, mps, ldim, 1.0_dp, &
                      wfcU(1,offsetU_back1(na)+1), 2*ldap, rtemp, &
                      ldimaxt, 1.0_dp, hpsi, 2*ldap )
                 !$acc end host_data
                 !
              ENDIF
           ENDIF
        ENDDO
        !
     ENDIF
     !
  ENDDO
  !
  IF (ANY(is_hubbard(:))) THEN
     !$acc exit data delete(vns_r)
     DEALLOCATE( vns_r )
  END IF
  IF (ANY(is_hubbard_back(:))) THEN
     !$acc exit data delete(vnsb_r)
     DEALLOCATE( vnsb_r )
  END IF
  !$acc exit data delete(rtemp)
  DEALLOCATE( rtemp )
  !$acc exit data delete(proj_r)
  DEALLOCATE ( proj_r )
  !
END SUBROUTINE vhpsi_gamma_acc
!
SUBROUTINE vhpsi_k_acc()
  !
  ! k-point version
  !
  IMPLICIT NONE
  !
  ! ... local variables
  !
  COMPLEX(DP), ALLOCATABLE :: proj_k(:,:)
  COMPLEX(DP), ALLOCATABLE :: ctemp(:,:), vns_c(:,:,:), vnsb_c(:,:,:)
  INTEGER :: na, nt, ldim, ldim0, ldimax, ldimaxt
  !
  ALLOCATE( proj_k(nwfcU, mps) )
  !$acc enter data create(proj_k)
  CALL calbec(offload_type, np, wfcU, psip, proj_k)
  !
  ldimax = 2*Hubbard_lmax+1
  ldimaxt = MAX(ldimax, ldmx_b)
  ALLOCATE( ctemp(ldimaxt,mps) )
  !$acc enter data create(ctemp)
  IF (ANY(is_hubbard(:))) THEN
     ALLOCATE( vns_c(ldimax,ldimax,nat) )
     !$acc enter data create(vns_c)
     vns_c = CMPLX(v%ns(:,:,current_spin,:),KIND=DP)
     !$acc update device(vns_c)
  ENDIF
  IF (ANY(is_hubbard_back(:))) THEN
     ALLOCATE( vnsb_c(ldmx_b,ldmx_b,nat) )
     !$acc enter data create(vnsb_c)
     vnsb_c = CMPLX(v%nsb(:,:,current_spin,:),KIND=DP)
     !$acc update device(vnsb_c)
  ENDIF
  !
  DO nt = 1, ntyp
     !
     ! Compute the action of the Hubbard potential on the KS wave functions:
     ! V_Hub |psip > = \sum v%ns |wfcU> <wfcU|psip>
     ! where v%ns = U ( delta/2 - rho%ns ) is computed in v_of_rho
     !
     IF ( is_hubbard(nt) ) THEN
        !  
        ldim = 2*Hubbard_l(nt) + 1
        !
        DO na = 1, nat
           IF ( nt == ityp(na) ) THEN
              !
              !$acc host_data use_device(proj_k,vns_c,ctemp)
              CALL MYZGEMM( 'N', 'N', ldim, mps, ldim, (1.0_dp,0.0_dp), &
                   vns_c(:,:,na), ldimax, proj_k(offsetU(na)+1,1), nwfcU, &
                   (0.0_dp,0.0_dp), ctemp, ldimaxt )
              !$acc end host_data
              !
              !$acc host_data use_device(wfcU, ctemp,hpsi)
              CALL MYZGEMM( 'N', 'N', np, mps, ldim, (1.0_dp,0.0_dp), &
                   wfcU(1,offsetU(na)+1), ldap, ctemp, ldimaxt, &
                   (1.0_dp,0.0_dp), hpsi, ldap)
              !$acc end host_data
              !
           ENDIF
        ENDDO
        !
     ENDIF
     !
     ! If the background is used then compute extra 
     ! contribution to the Hubbard potential
     !
     IF ( is_hubbard_back(nt) ) THEN
        !
        ldim = ldim_back(nt)
        !
        DO na = 1, nat
           IF ( nt == ityp(na) ) THEN
              !
              !
              ldim = 2*Hubbard_l2(nt)+1
              !
              !$acc host_data use_device(proj_k,vnsb_c,ctemp)
              CALL MYZGEMM( 'N', 'N', ldim,mps,ldim, (1.0_dp,0.0_dp),     &
                   vnsb_c(:,:,na), ldmx_b, proj_k(offsetU_back(na)+1,1), &
                   nwfcU, (0.0_dp,0.0_dp), ctemp, ldimaxt )
              !$acc end host_data
              !
              !$acc host_data use_device(wfcU,ctemp,hpsi)
              CALL MYZGEMM( 'N', 'N', np, mps, ldim, (1.0_dp,0.0_dp), &
                   wfcU(1,offsetU_back(na)+1), ldap, ctemp,           &
                   ldimaxt, (1.0_dp,0.0_dp), hpsi, ldap )
              !$acc end host_data
              !
              IF (backall(nt)) THEN
                 !
                 ldim0 = 2*Hubbard_l2(nt)+1
                 ldim  = 2*Hubbard_l3(nt)+1
                 !
                 !$acc host_data use_device(proj_k,vnsb_c,ctemp)
                 CALL MYZGEMM( 'N', 'N', ldim,mps,ldim,(1.0_dp,0.0_dp), &
                      vnsb_c(ldim0+1,ldim0+1,na),ldmx_b,                   &
                      proj_k(offsetU_back1(na)+1,1), nwfcU,             &
                      (0.0_dp,0.0_dp), ctemp, ldimaxt )
                 !$acc end host_data
                 ! 
                 !$acc host_data use_device(wfcU, ctemp,hpsi)
                 CALL MYZGEMM( 'N', 'N', np, mps, ldim, (1.0_dp,0.0_dp), &
                      wfcU(1,offsetU_back1(na)+1), ldap, ctemp,          &
                      ldimaxt, (1.0_dp,0.0_dp), hpsi, ldap )
                 !$acc end host_data
                 !
              ENDIF
              !
           ENDIF
        ENDDO
        !
     ENDIF
     !
  ENDDO
  !
  IF (ANY(is_hubbard(:))) THEN
     !$acc exit data delete(vns_c)
     DEALLOCATE( vns_c )
  END IF
  IF (ANY(is_hubbard_back(:))) THEN
     !$acc exit data delete(vnsb_c)
     DEALLOCATE( vnsb_c )
  END IF
  !$acc exit data delete(ctemp)
  DEALLOCATE( ctemp )
  !$acc exit data delete(proj_k)
  DEALLOCATE( proj_k )
  !
END SUBROUTINE vhpsi_k_acc
!
!-------------------------------------------------------------------------
END SUBROUTINE vhpsi_U
!-------------------------------------------------------------------------

!-------------------------------------------------------------------------
SUBROUTINE vhpsi_nc_acc( lda, np, mps, psi, hpsi )
  !-----------------------------------------------------------------------
  !! Noncollinear version of \(\texttt{vhpsi} routine (A. Smogunov).
  !! Ported to GPU (DFT+U only)
  !
  USE kinds,            ONLY: dp
  USE control_flags,    ONLY: offload_type
  USE ldaU,             ONLY: Hubbard_lmax, Hubbard_l, is_Hubbard, nwfcU, &
                              wfcU, offsetU, lda_plus_u_kind
  USE scf,              ONLY: v
  USE ions_base,        ONLY: nat, ntyp => nsp, ityp
  USE noncollin_module, ONLY: npol
  USE mp_bands,         ONLY: intra_bgrp_comm
  USE mp,               ONLY: mp_sum
  USE lsda_mod,         ONLY: nspin
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: lda
  !! leading dimension of arrays psi, hpsi
  INTEGER, INTENT(IN) :: np
  !! true dimension of psi, hpsi
  INTEGER, INTENT(IN) :: mps
  !! number of states psi
  COMPLEX(dp), INTENT(IN) :: psi(lda*npol,mps)
  !! the wavefunction
  COMPLEX(dp), INTENT(INOUT) :: hpsi(lda*npol,mps)
  !! Hamiltonian dot psi
  !
  ! ... local variables
  !
  INTEGER :: ibnd, na, nwfc, is1, is2, nt, m1, m2, ldim
  COMPLEX(dp), ALLOCATABLE :: proj(:,:)
  COMPLEX(dp), ALLOCATABLE :: ctemp(:,:), vns(:,:)
  !
  !
  IF ( lda_plus_u_kind == 2 ) CALL errore('vhpsi_nc','incorrectly called',1)
  !
  !$acc data present(wfcU,psi,hpsi)
  ALLOCATE( proj(nwfcU, mps) )
  !$acc enter data create(proj)
  !
  ! calculate proj=<psi_at | psi_k> 
  !$acc host_data use_device(wfcU,psi,proj)
  CALL MYZGEMM ('C', 'N', nwfcU, mps, lda*npol, (1.0_dp, 0.0_dp), wfcU, &
                    lda*npol, psi, lda*npol, (0.0_dp, 0.0_dp),  proj, nwfcU)
  CALL mp_sum ( proj, intra_bgrp_comm )
  !$acc end host_data
  !
  DO nt = 1, ntyp
     !
     ! Compute the action of the Hubbard potential on the KS wave functions:
     ! V_Hub |psi > = \sum v%ns |wfcU> <wfcU|psi>
     ! where v%ns = U ( delta/2 - rho%ns ) is computed in v_of_rho
     !
     IF ( is_hubbard(nt) ) THEN
        !  
        ldim = 2*Hubbard_l(nt) + 1
        !
        ALLOCATE ( ctemp(ldim*npol,mps) )
        ALLOCATE ( vns (ldim*npol,ldim*npol) )
        !$acc data create(vns, ctemp)
        !
        DO na = 1, nat
           IF ( nt == ityp(na) ) THEN
              !
              do is1 = 1, npol
                 do is2 = 1, npol
                    DO m2 = 1, ldim
                       DO m1 = 1, ldim
                          vns (m1+ldim*(is1-1), m2+ldim*(is2-1))  &
                               = v%ns_nc(m1,m2,npol*(is1-1)+is2,na)
                       ENDDO
                    ENDDO
                 enddo
              enddo
              !$acc update device(vns)
              !
              !$acc host_data use_device(vns,proj,ctemp)
              CALL MYZGEMM ('n','n', ldim*npol, mps, ldim*npol,(1.0_dp,0.0_dp),&
                   vns, ldim*npol, proj(offsetU(na)+1,1),&
                   nwfcU,(0.0_dp,0.0_dp),ctemp,ldim*npol)
              !$acc end host_data
              !
              !$acc host_data use_device(wfcU,ctemp,hpsi)
              CALL MYZGEMM ('n','n', lda*npol, mps, ldim*npol, (1.0_dp,0.0_dp),&
                   wfcU(1,offsetU(na)+1), lda*npol, ctemp, ldim*npol, &
                   (1.0_dp,0.0_dp), hpsi, lda*npol)
              !$acc end host_data
              !
           ENDIF
        ENDDO
        !
        !$acc end data
        DEALLOCATE (vns)
        DEALLOCATE (ctemp)
        !
     ENDIF
     !
  ENDDO
  !
  !$acc exit data delete(proj)
  deallocate (proj)
  !$acc end data
  !
  RETURN
  !
END SUBROUTINE vhpsi_nc_acc

