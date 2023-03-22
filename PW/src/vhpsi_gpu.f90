!
! Copyright (C) 2001-2022 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#if !defined(__CUDA)
#define cublasDgemm dgemm
#define cublasZgemm zgemm
#endif
!-----------------------------------------------------------------------
SUBROUTINE vhpsi_gpu( ldap, np, mps, psip_d, hpsi_d )
  !-----------------------------------------------------------------------
  !! This routine computes the Hubbard potential applied to the electronic
  !! structure of the current k-point. The result is added to hpsi.
  !
  USE kinds,         ONLY : DP
  USE becmod,        ONLY : bec_type, calbec, allocate_bec_type, &
                            deallocate_bec_type
  USE ldaU,          ONLY : Hubbard_lmax, Hubbard_l, is_hubbard,   &
                            nwfcU, wfcU, offsetU, lda_plus_u_kind, &
                            is_hubbard_back, Hubbard_l2, offsetU_back, &
                            backall, offsetU_back1
  USE lsda_mod,      ONLY : current_spin
  USE scf,           ONLY : v
  USE ions_base,     ONLY : nat, ntyp => nsp, ityp
  USE control_flags, ONLY : gamma_only
  USE mp,            ONLY : mp_sum
  !
  USE becmod_gpum,      ONLY : bec_type_d
  USE becmod_subs_gpum, ONLY : allocate_bec_type_gpu, deallocate_bec_type_gpu, &
                               calbec_gpu
  USE device_memcpy_m,    ONLY: dev_memcpy, dev_memset
  !
#if defined(__CUDA)
  USE cudafor
  USE cublas
#endif  
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: ldap
  !! leading dimension of arrays psip, hpsi
  INTEGER, INTENT(IN) :: np
  !! true dimension of psip, hpsi
  INTEGER, INTENT(IN) :: mps
  !! number of states psip
  COMPLEX(DP), INTENT(IN) :: psip_d(ldap,mps)
  !! the wavefunction
  COMPLEX(DP), INTENT(INOUT) :: hpsi_d(ldap,mps)
  !! Hamiltonian dot psi
  !
  ! ... local variables
  !
  INTEGER :: dimwf1
  TYPE(bec_type_d), TARGET :: proj_d
  COMPLEX(DP), ALLOCATABLE :: wfcU_d(:,:)
  !
#if defined(__CUDA)
  attributes(DEVICE) :: wfcU_d, psip_d, hpsi_d
#endif
  !
  CALL start_clock_gpu( 'vhpsi' )
  !
  ! Offset of atomic wavefunctions initialized in setup and stored in offsetU
  !
  ! Allocate the array proj
  CALL allocate_bec_type_gpu( nwfcU, mps, proj_d ) 
  !
  dimwf1 = SIZE(wfcU(:,1))
  ALLOCATE( wfcU_d(dimwf1,nwfcU) )
  CALL dev_memcpy(wfcU_d, wfcU)
  !
  ! proj = <wfcU|psip>
  CALL calbec_gpu( np, wfcU_d, psip_d, proj_d )
  !
  IF ( lda_plus_u_kind.EQ.0 .OR. lda_plus_u_kind.EQ.1 ) THEN
     CALL vhpsi_U_gpu()  ! DFT+U
  ELSEIF ( lda_plus_u_kind.EQ.2 ) THEN
     CALL errore('vhpsi', 'DFT+U+V case not implemented for GPU', 1 )
  ENDIF
  !
  CALL deallocate_bec_type_gpu( proj_d )
  !
  DEALLOCATE( wfcU_d )
  !
  CALL stop_clock_gpu( 'vhpsi' )
  !
  RETURN
  !
CONTAINS
  !
SUBROUTINE vhpsi_U_gpu()
  !
  ! This routine applies the Hubbard potential with U_I
  ! to the KS wave functions. 
  !
  USE ldaU,      ONLY : ldim_back, ldmx_b, Hubbard_l3
  !
  IMPLICIT NONE
  !
  INTEGER :: na, nt, ldim, ldim0, ldimax, ldimaxt
  !
  REAL(DP),    ALLOCATABLE :: rtemp_d(:,:), vns_d(:,:,:),  vnsb_d(:,:,:)
  COMPLEX(DP), ALLOCATABLE :: ctemp_d(:,:), vaux_d(:,:,:), vauxb_d(:,:,:)
  !
#if defined(__CUDA)
  attributes(DEVICE) :: vns_d, vnsb_d, rtemp_d, ctemp_d, vaux_d, vauxb_d
#endif   
  !
  ldimax = 2*Hubbard_lmax+1
  ldimaxt = MAX(ldimax, ldmx_b)
  !
  IF ( ANY(is_hubbard(:)) .OR. ANY(is_hubbard_back(:)) ) THEN
    IF (gamma_only) THEN
      ALLOCATE( rtemp_d(ldimaxt,mps) )
      IF (ANY(is_hubbard(:))) THEN
        ALLOCATE( vns_d(ldimax,ldimax,nat) )
        vns_d = v%ns(:,:,current_spin,:)
      ENDIF
      IF (ANY(is_hubbard_back(:))) THEN
        ALLOCATE( vnsb_d(ldmx_b,ldmx_b,nat) )
        vnsb_d = v%nsb(:,:,current_spin,:)
      ENDIF
    ELSE
      ALLOCATE( ctemp_d(ldimaxt,mps) )
      IF (ANY(is_hubbard(:))) THEN
        ALLOCATE( vaux_d(ldimax,ldimax,nat) )
        vaux_d = CMPLX(v%ns(:,:,current_spin,:))
      ENDIF
      IF (ANY(is_hubbard_back(:))) THEN
         ALLOCATE( vauxb_d(ldmx_b,ldmx_b,nat) )
         vauxb_d = CMPLX(v%nsb(:,:,current_spin,:))
      ENDIF
    ENDIF
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
              IF (gamma_only) THEN
                 !
                 CALL cublasDgemm( 'N','N', ldim,mps,ldim, 1.0_dp, &
                      vns_d(1,1,na), ldimax, &
                      proj_d%r_d(offsetU(na)+1,1), nwfcU, 0.0_dp, rtemp_d, ldimaxt )
                 !
                 CALL cublasDgemm( 'N','N', 2*np, mps, ldim, 1.0_dp, &
                      wfcU_d(1,offsetU(na)+1), 2*ldap, rtemp_d, ldimaxt, &
                      1.0_dp, hpsi_d, 2*ldap )
                 !
              ELSE
                 !
                 CALL cublasZgemm( 'N', 'N', ldim, mps, ldim, (1.0_dp,0.0_dp), &
                      vaux_d(:,:,na), ldimax, proj_d%k_d(offsetU(na)+1,1), nwfcU, &
                      (0.0_dp,0.0_dp), ctemp_d, ldimaxt )
                 !
                 CALL cublasZgemm( 'N', 'N', np, mps, ldim, (1.0_dp,0.0_dp), &
                      wfcU_d(1,offsetU(na)+1), ldap, ctemp_d, ldimaxt, &
                      (1.0_dp,0.0_dp), hpsi_d, ldap)
                 !
              ENDIF
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
              IF (gamma_only) THEN
                 !
                 ldim = 2*Hubbard_l2(nt)+1
                 !
                 CALL cublasDgemm( 'N','N', ldim,mps,ldim, 1.0_dp, &
                      vnsb_d(1,1,na),ldmx_b, &
                      proj_d%r_d(offsetU_back(na)+1,1), &
                      nwfcU, 0.0_dp, rtemp_d, ldimaxt )
                 !
                 CALL cublasDgemm( 'N','N', 2*np, mps, ldim, 1.0_dp, &
                      wfcU_d(1,offsetU_back(na)+1), 2*ldap, rtemp_d, &
                      ldimaxt, 1.0_dp, hpsi_d, 2*ldap )
                 !
                 IF (backall(nt)) THEN
                    !
                    ldim0 = 2*Hubbard_l2(nt)+1
                    ldim  = 2*Hubbard_l3(nt)+1
                    !
                    CALL cublasDgemm( 'N', 'N', ldim,mps,ldim, 1.0_dp,     &
                         vnsb_d(ldim0+1,ldim0+1,na),                       &
                         ldim_back(nt), proj_d%r_d(offsetU_back1(na)+1,1), &
                         nwfcU, 0.0_dp, rtemp_d, ldimaxt )
                    !
                    CALL cublasDgemm( 'N', 'N', 2*np, mps, ldim, 1.0_dp, &
                         wfcU_d(1,offsetU_back1(na)+1), 2*ldap, rtemp_d, &
                         ldimaxt, 1.0_dp, hpsi_d, 2*ldap )
                    !
                 ENDIF
                 !
              ELSE
                 !
                 ldim = 2*Hubbard_l2(nt)+1
                 !
                 CALL cublasZgemm( 'N', 'N', ldim,mps,ldim, (1.0_dp,0.0_dp),     &
                      vauxb_d(:,:,na), ldmx_b, proj_d%k_d(offsetU_back(na)+1,1), &
                      nwfcU, (0.0_dp,0.0_dp), ctemp_d, ldimaxt )
                 !
                 CALL cublasZgemm( 'N', 'N', np, mps, ldim, (1.0_dp,0.0_dp), &
                      wfcU_d(1,offsetU_back(na)+1), ldap, ctemp_d,           &
                      ldimaxt, (1.0_dp,0.0_dp), hpsi_d, ldap )
                 !
                 IF (backall(nt)) THEN
                    !
                    ldim0 = 2*Hubbard_l2(nt)+1
                    ldim  = 2*Hubbard_l3(nt)+1
                    !
                    CALL cublasZgemm( 'N', 'N', ldim,mps,ldim,(1.0_dp,0.0_dp), &
                         vauxb_d(ldim0+1,ldim0+1,na),ldmx_b,                   &
                         proj_d%k_d(offsetU_back1(na)+1,1), nwfcU,             &
                         (0.0_dp,0.0_dp), ctemp_d, ldimaxt )
                    ! 
                    CALL cublasZgemm( 'N', 'N', np, mps, ldim, (1.0_dp,0.0_dp), &
                         wfcU_d(1,offsetU_back1(na)+1), ldap, ctemp_d,          &
                         ldimaxt, (1.0_dp,0.0_dp), hpsi_d, ldap )
                    !
                 ENDIF
                 !
              ENDIF
           ENDIF
        ENDDO
        !
     ENDIF
     !
  ENDDO
  !
  IF ( ANY(is_hubbard(:)) .OR. ANY(is_hubbard_back(:)) ) THEN
    IF (gamma_only) THEN
      DEALLOCATE( rtemp_d )
      IF (ANY(is_hubbard(:))) DEALLOCATE( vns_d )
      IF (ANY(is_hubbard_back(:))) DEALLOCATE( vnsb_d )
    ELSE
      DEALLOCATE( ctemp_d )
      IF (ANY(is_hubbard(:))) DEALLOCATE( vaux_d )
      IF (ANY(is_hubbard_back(:))) DEALLOCATE( vauxb_d )
    ENDIF
  ENDIF
  !
  !
  RETURN
  !
END SUBROUTINE vhpsi_U_gpu
!
END SUBROUTINE vhpsi_gpu
!-------------------------------------------------------------------------
