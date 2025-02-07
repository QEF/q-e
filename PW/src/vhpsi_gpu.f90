!
! Copyright (C) 2001-2025 Quantum ESPRESSO Foundation
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
SUBROUTINE vhpsi_gpu( ldap, np, mps, psip, hpsi )
  !-----------------------------------------------------------------------
  !! This routine computes the Hubbard potential applied to the electronic
  !! structure of the current k-point. The result is added to hpsi.
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
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: ldap
  !! leading dimension of arrays psip, hpsi
  INTEGER, INTENT(IN) :: np
  !! true dimension of psip, hpsi
  INTEGER, INTENT(IN) :: mps
  !! number of states psip
  COMPLEX(DP), INTENT(IN) :: psip(ldap,mps)
  !! the wavefunction
  COMPLEX(DP), INTENT(INOUT) :: hpsi(ldap,mps)
  !! Hamiltonian dot psi
  !
  IF ( lda_plus_u_kind  == 2 ) &
       CALL errore('vhpsi', 'DFT+U+V case not implemented for GPU', 1 )
  IF ( lda_plus_u_kind /= 0 .AND. lda_plus_u_kind /= 1 ) RETURN
  IF ( .NOT. ANY(is_hubbard(:)) .AND. .NOT.ANY(is_hubbard_back(:)) ) RETURN
  !
  CALL start_clock( 'vhpsi' )
  !
  ! Offset of atomic wavefunctions initialized in setup and stored in offsetU
  !
  !$acc data present(wfcU)
  !
  ! proj = <wfcU|psip>
  IF (gamma_only) THEN
     CALL vhpsi_gamma_acc ()
  ELSE
     CALL vhpsi_k_acc ()
  ENDIF
  !
  !$acc end data
  !
  CALL stop_clock( 'vhpsi' )
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
END SUBROUTINE vhpsi_gpu
!-------------------------------------------------------------------------
