!
! Copyright (C) 2001-2025 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
SUBROUTINE v_hubbard ( noncolin, rho, v, eth )
  !---------------------------------------------------------------------
  !
  USE kinds, ONLY : dp
  USE ldaU,  ONLY : lda_plus_u, lda_plus_u_kind, ldmx_b, v_nsg, &
                    Hubbard_l, Hubbard_lmax, apply_U, orbital_resolved 
  USE scf,   ONLY : scf_type
#if defined (__OSCDFT)
  USE plugin_flags,     ONLY : use_oscdft
  USE oscdft_base,      ONLY : oscdft_ctx
  USE oscdft_functions, ONLY : oscdft_v_constraint
#endif  
  !
  IMPLICIT NONE
  !
  TYPE(scf_type), INTENT(INOUT) :: rho
  !! the valence charge
  TYPE(scf_type), INTENT(INOUT) :: v
  !! the scf (Hxc) potential
  REAL(dp), INTENT(OUT) :: eth
  !! the hubbard energy
  LOGICAL, INTENT(IN) :: noncolin
  !! noncollinear calculation
  !
  REAL(dp) :: eth1
  !! the hubbard energy coming from the background states
  !
  IF ( lda_plus_u ) THEN
     !
     IF (lda_plus_u_kind == 0) THEN
        !
        ! DFT+U (simplified)
        !
        IF (noncolin) THEN
           IF (orbital_resolved) THEN
              CALL v_hubbard_resolved_nc (rho%ns_nc, v%ns_nc, eth)
           ELSE
              CALL v_hubbard_nc (rho%ns_nc, v%ns_nc, eth)
           ENDIF
        ELSE
           IF (orbital_resolved) THEN
              CALL v_hubbard_resolved(rho%ns, v%ns, eth)
           ELSE
              CALL v_hubbard_plain (rho%ns, v%ns, eth)
           ENDIF
        ENDIF
        !
        ! Background
        IF (ldmx_b.GT.0) THEN
           CALL v_hubbard_b (rho%nsb, v%nsb, eth1)
           eth = eth + eth1
        ENDIF
        !
     ELSEIF (lda_plus_u_kind == 1) THEN
        !
        ! DFT+U (full)
        !
        IF (noncolin) THEN
           CALL v_hubbard_full_nc (rho%ns_nc, v%ns_nc, eth)
        ELSE
           CALL v_hubbard_full (rho%ns, v%ns, eth)
        ENDIF
        !
     ELSEIF (lda_plus_u_kind == 2) THEN
        !
        ! DFT+U+V (simplified)
        !
        IF (noncolin) THEN
           CALL v_hubbard_extended_nc (rho%nsg, v_nsg, eth)
        ELSE
           CALL v_hubbard_extended (rho%nsg, v_nsg, eth)
        ENDIF
     ELSE
        !
        CALL errore('v_of_rho', 'Not allowed value of lda_plus_u_kind',1)
        !
     ENDIF
     !
  ENDIF
  !
#if defined (__OSCDFT)
  IF (use_oscdft .AND. (oscdft_ctx%inp%oscdft_type==2)) THEN
     IF (lda_plus_u_kind == 0) THEN
        CALL oscdft_v_constraint (oscdft_ctx, Hubbard_lmax, Hubbard_l, rho%ns, v%ns, eth)
     ELSEIF (lda_plus_u_kind == 2) THEN
        CALL oscdft_v_constraint_extended (rho%nsg, v_nsg, eth)
     ENDIF
  ENDIF
#endif
  !
END SUBROUTINE v_hubbard
!-----------------------------------------------------------------------
SUBROUTINE v_hubbard_plain ( ns, v_hub, eth )
  !---------------------------------------------------------------------
  !
  !! Computes Hubbard potential and Hubbard energy.
  !! DFT+U: Simplified rotationally-invariant formulation by
  !! Dudarev et al., Phys. Rev. B 57, 1505 (1998).
  !! DFT+U+J0: B. Himmetoglu et al., Phys. Rev. B 84, 115108 (2011).
  !
  USE kinds,                ONLY : DP
  USE ions_base,            ONLY : nat, ityp
  USE ldaU,                 ONLY : Hubbard_lmax, Hubbard_l, Hubbard_U, &
                                   Hubbard_alpha, Hubbard_J0, Hubbard_beta, &
                                   dfpt_hub
  USE lsda_mod,             ONLY : nspin
  USE control_flags,        ONLY : iverbosity
  USE io_global,            ONLY : stdout
  !
  IMPLICIT NONE
  !
  REAL(DP), INTENT(IN)  :: ns(2*Hubbard_lmax+1,2*Hubbard_lmax+1,nspin,nat)
  !! Occupation matrix
  REAL(DP), INTENT(OUT) :: v_hub(2*Hubbard_lmax+1,2*Hubbard_lmax+1,nspin,nat)
  !! Hubbard potential
  REAL(DP), INTENT(OUT) :: eth
  !! Hubbard energy
  COMPLEX(DP) :: check
  !
  !  ... local variables
  !
  REAL(DP) :: effU, sgn(2) 
  INTEGER  :: is, isop, na, nt, m1, m2
  !
  eth    = 0.d0
  sgn(1) =  1.d0  
  sgn(2) = -1.d0
  !
  v_hub(:,:,:,:) = 0.d0
  check = (0.0,0.0)
  !
  DO na = 1, nat
     !
     nt = ityp (na)
     !
     IF (Hubbard_U(nt) /= 0.d0 .OR. Hubbard_alpha(nt) /= 0.d0) THEN
        !
        IF (Hubbard_J0(nt) /= 0.d0) THEN
           effU = Hubbard_U(nt) - Hubbard_J0(nt)
        ELSE
           effU = Hubbard_U(nt)
        ENDIF 
        ! 
        DO is = 1, nspin
           DO m1 = 1, 2 * Hubbard_l(nt) + 1
              eth = eth + ( Hubbard_alpha(nt) + 0.5D0*effU )*ns(m1,m1,is,na)
              v_hub(m1,m1,is,na) = v_hub(m1,m1,is,na) + &
                                   Hubbard_alpha(nt)  + 0.5D0*effU
               check = check + Hubbard_alpha(nt)  + 0.5D0*effU
              DO m2 = 1, 2 * Hubbard_l(nt) + 1
                 eth = eth - 0.5D0 * effU * ns(m2,m1,is,na)* ns(m1,m2,is,na)
                 v_hub(m1,m2,is,na) = v_hub(m1,m2,is,na) - &
                                        effU * ns(m2,m1,is,na)
                 check = check -effU * ns(m2,m1,is,na)
              ENDDO
           ENDDO
        ENDDO
        !
     ENDIF
     !
     IF (Hubbard_J0(nt) /= 0.d0 .OR. Hubbard_beta(nt) /= 0.d0) THEN
        !
        DO is = 1, nspin
           isop = 1
           IF ( nspin == 2 .AND. is == 1) isop = 2
           DO m1 = 1, 2 * Hubbard_l(nt) + 1
              eth = eth + sgn(is)*Hubbard_beta(nt) * ns(m1,m1,is,na)
              v_hub(m1,m1,is,na) = v_hub(m1,m1,is,na) + sgn(is)*Hubbard_beta(nt)
              DO m2 = 1, 2*Hubbard_l(nt)+1
                 eth = eth + 0.5D0*Hubbard_J0(nt)*ns(m2,m1,is,na)*ns(m1,m2,isop,na)
                 v_hub(m1,m2,is,na) = v_hub(m1,m2,is,na) + Hubbard_J0(nt) * &
                                                            ns(m2,m1,isop,na)
              ENDDO
           ENDDO
        ENDDO
        !
     END IF
     !   
  ENDDO
  !
  IF (nspin==1) eth = 2.d0 * eth
  !
  ! Hubbard energy
  !
  IF ( iverbosity > 0 .AND. .NOT.dfpt_hub ) THEN
     WRITE(stdout,'(/5x,"HUBBARD ENERGY = ",f9.4,1x," (Ry)")') eth
     !write(stdout,*) "check coll U", check
  ENDIF
  !
  RETURN
  !
END SUBROUTINE v_hubbard_plain
!-----------------------------------------------------------------------

!----------------------------------------------------------------------
SUBROUTINE v_hubbard_nc( ns, v_hub, eth )
  !---------------------------------------------------------------------
  !
  !! Computes Hubbard potential and Hubbard energy in the \textbf{noncollinear formulation}.
  !! DFT+U: Simplified rotationally-invariant formulation by
  !! Dudarev et al., Phys. Rev. B 57, 1505 (1998).
  !! DFT+U+J0: B. Himmetoglu et al., Phys. Rev. B 84, 115108 (2011).
  !
  USE kinds,                ONLY : DP
  USE ions_base,            ONLY : nat, ityp
  USE ldaU,                 ONLY : Hubbard_lmax, Hubbard_l, Hubbard_U, &
                                   Hubbard_alpha, Hubbard_J0, Hubbard_beta, &
                                   dfpt_hub
  USE lsda_mod,             ONLY : nspin
  USE control_flags,        ONLY : iverbosity
  USE io_global,            ONLY : stdout
  !
  IMPLICIT NONE
  !
  COMPLEX(DP) :: ns(2*Hubbard_lmax+1,2*Hubbard_lmax+1,nspin,nat)
  !! Occupation matrix
  COMPLEX(DP) :: v_hub(2*Hubbard_lmax+1,2*Hubbard_lmax+1,nspin,nat)
  !! Hubbard potential
  REAL(DP), INTENT(OUT) :: eth
  !! Hubbard energy
  !
  !  ... local variables
  !
  INTEGER  :: is, is1, na, nt, m1, m2
  !
  eth    = 0.d0
  v_hub(:,:,:,:) = 0.d0
  !
  DO na = 1, nat
     !
     nt = ityp (na)
     !
     IF (Hubbard_U(nt) /= 0.d0) THEN
        DO is = 1, nspin
           !
           IF (is == 2) THEN
            is1 = 3
           ELSEIF (is == 3) THEN
            is1 = 2
           ELSE
            is1 = is
           ENDIF
           !
           ! Non spin-flip contribution
           ! (diagonal [spin indexes] occupancy matrices)
           IF (is1 == is) THEN
              !     
              ! diagonal part [spin indexes]     
              DO m1 = 1, 2*Hubbard_l(nt) + 1
                 ! Hubbard energy
                 eth = eth + ( Hubbard_alpha(nt) + 0.5D0*Hubbard_U(nt) )&
                             * ns(m1,m1,is,na)
                 ! Hubbard potential
                 v_hub(m1,m1,is,na) = v_hub(m1,m1,is,na) + & 
                                      Hubbard_alpha(nt) + 0.5D0*Hubbard_U(nt)
                 ! 
                 ! NON-diagonal part [spin indexes]
                 DO m2 = 1, 2 * Hubbard_l(nt) + 1
                    ! Hubbard energy
                    eth = eth - 0.5D0 * Hubbard_U(nt) * ns(m1,m2,is,na)*ns(m2,m1,is,na)
                    ! Hubbard potential
                    v_hub(m1,m2,is,na) = v_hub(m1,m2,is,na) - &
                                         Hubbard_U(nt) * ns(m2,m1,is,na)
                 ENDDO
              ENDDO  
           !
           ! Spin-flip contribution
           ! (NON-diagonal [spin indexes] occupancy matrices)   
           ELSE
              DO m1 = 1, 2*Hubbard_l(nt) + 1
                 DO m2 = 1, 2 * Hubbard_l(nt) + 1 
                    ! Hubbard energy
                    eth = eth - 0.5D0 * Hubbard_U(nt) &
                          * ns(m1,m2,is,na)*ns(m2,m1,is1,na)
                    ! Hubbard potential
                    v_hub(m1,m2,is,na) = v_hub(m1,m2,is,na) - &
                                         Hubbard_U(nt) * ns(m2,m1,is1,na)
                 ENDDO
              ENDDO         
           ENDIF
        ENDDO  
     ENDIF
     !
  ENDDO
  !
  ! Hubbard energy
  !
  IF ( iverbosity > 0 ) THEN
     WRITE(stdout,'(/5x,"HUBBARD ENERGY = ",f9.4,1x," (Ry)")') eth
  ENDIF
  !
  RETURN
  !
END SUBROUTINE v_hubbard_nc
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
SUBROUTINE v_hubbard_b (ns, v_hub, eth)
  !-------------------------------------------------------------------------
  !
  !! Computes Hubbard potential and Hubbard energy for background states.
  !! DFT+U: Simplified rotationally-invariant formulation by
  !! Dudarev et al., Phys. Rev. B 57, 1505 (1998).
  !
  USE kinds,                ONLY : DP
  USE ions_base,            ONLY : nat, ityp
  USE ldaU,                 ONLY : Hubbard_J0, Hubbard_beta, Hubbard_U2,  &
                                   ldim_back, ldmx_b, Hubbard_alpha_back, &
                                   is_hubbard_back, dfpt_hub
  USE lsda_mod,             ONLY : nspin
  USE control_flags,        ONLY : iverbosity
  USE io_global,            ONLY : stdout

  IMPLICIT NONE
  !
  REAL(DP), INTENT(IN)  :: ns(ldmx_b,ldmx_b,nspin,nat)
  REAL(DP), INTENT(OUT) :: v_hub(ldmx_b,ldmx_b,nspin,nat)
  REAL(DP), INTENT(OUT) :: eth
  REAL(DP) :: effU
  INTEGER :: is, is1, na, nt, m1, m2, m3, m4
  !
  eth = 0.d0
  !
  v_hub(:,:,:,:) = 0.d0
  !
  DO na = 1, nat
     !
     nt = ityp (na)
     !
     IF (Hubbard_J0(nt).NE.0.d0) &
          CALL errore('v_hubbard_b', 'J0 is not supported in DFT+U with multiple channels per atomic type',1)
     !
     IF (Hubbard_beta(nt).NE.0.d0) &
     CALL errore('v_hubbard_b', 'Hubbard_beta is not supported in DFT+U with multiple channels per atomic type',1) 
     !
     IF (is_hubbard_back(nt)) THEN
        !
        effU = Hubbard_U2(nt)
        !
        DO is = 1, nspin
           !
           DO m1 = 1, ldim_back(nt) 
              !
              eth = eth + ( Hubbard_alpha_back(nt) + 0.5D0 * effU ) * &
                              ns(m1,m1,is,na)
              v_hub(m1,m1,is,na) = v_hub(m1,m1,is,na) + &
                            ( Hubbard_alpha_back(nt) + 0.5D0 * effU )
              !
              DO m2 = 1, ldim_back(nt)
                 !
                 eth = eth - 0.5D0 * effU * &
                            ns(m2,m1,is,na)* ns(m1,m2,is,na)
                 v_hub(m1,m2,is,na) = v_hub(m1,m2,is,na) - &
                                       effU * ns(m2,m1,is,na)
                 !
              ENDDO
              !
           ENDDO
           !
        ENDDO
        !
     ENDIF
     !
  ENDDO
  !
  IF (nspin.EQ.1) eth = 2.d0 * eth
  !
  ! Hubbard energy
  !
  IF ( iverbosity > 0 .AND. .NOT.dfpt_hub ) THEN
     WRITE(stdout,'(/5x,"HUBBARD BACKGROUND ENERGY = ",f9.4,1x," (Ry)")') eth
  ENDIF
  !
  RETURN
  !
END SUBROUTINE v_hubbard_b
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
SUBROUTINE v_hubbard_full( ns, v_hub, eth )
  !---------------------------------------------------------------------
  !
  !! Computes Hubbard potential and Hubbard energy.
  !! DFT+U(+J) : Formulation by Liechtenstein et al., Phys. Rev. B 52, R5467 (1995).
  !
  USE kinds,                ONLY : DP
  USE ions_base,            ONLY : nat, ityp
  USE ldaU,                 ONLY : Hubbard_lmax, Hubbard_l, Hubbard_U, &
                                   Hubbard_J, Hubbard_alpha
  USE lsda_mod,             ONLY : nspin
  USE control_flags,        ONLY : iverbosity
  USE io_global,            ONLY : stdout
  !
  IMPLICIT NONE
  !
  REAL(DP), INTENT(IN)  :: ns(2*Hubbard_lmax+1,2*Hubbard_lmax+1,nspin,nat)
  !! Occupation matrix
  REAL(DP), INTENT(OUT) :: v_hub(2*Hubbard_lmax+1,2*Hubbard_lmax+1,nspin,nat)
  !! Hubbard potential
  REAL(DP), INTENT(OUT) :: eth
  !! Hubbard energy
  !
  !  ... local variables
  !
  REAL(DP) :: n_tot, n_spin, eth_dc, eth_u, mag2
  INTEGER  :: is, isop, is1, na, nt, m1, m2, m3, m4
  REAL(DP), ALLOCATABLE :: u_matrix(:,:,:,:)
  !
  ALLOCATE( u_matrix(2*Hubbard_lmax+1, 2*Hubbard_lmax+1, 2*Hubbard_lmax+1, 2*Hubbard_lmax+1) )
  !
  eth    = 0.d0
  eth_dc = 0.d0
  eth_u  = 0.d0
  !
  v_hub(:,:,:,:) = 0.d0
  !
  DO na = 1, nat
     !
     nt = ityp (na)
     !
     IF (Hubbard_U(nt)/=0.d0) THEN
        !
        ! Initialize U(m1,m2,m3,m4) matrix 
        !
        CALL hubbard_matrix( Hubbard_lmax, Hubbard_l(nt), Hubbard_U(nt), &
                               Hubbard_J(1,nt), u_matrix )
        !
        ! Total N and M^2 for DC (double counting) term
        !
        n_tot = 0.d0
        !
        DO is = 1, nspin
           DO m1 = 1, 2 * Hubbard_l(nt) + 1
              n_tot = n_tot + ns(m1,m1,is,na)
           ENDDO
        ENDDO
        !
        IF (nspin==1) n_tot = 2.d0 * n_tot
        !
        mag2  = 0.d0
        ! 
        IF (nspin==2) THEN
           DO m1 = 1, 2 * Hubbard_l(nt) + 1
              mag2 = mag2 + ns(m1,m1,1,na) - ns(m1,m1,2,na)
           ENDDO
        ENDIF
        mag2  = mag2**2
        !
        ! Hubbard energy: DC term
        !
        eth_dc = eth_dc + 0.5d0*( Hubbard_U(nt)*n_tot*(n_tot-1.d0) - &
                                  Hubbard_J(1,nt)*n_tot*(0.5d0*n_tot-1.d0) - &
                                  0.5d0*Hubbard_J(1,nt)*mag2 )
        !
        DO is = 1, nspin
           !
           ! n_spin = up/down N
           !
           n_spin = 0.d0
           DO m1 = 1, 2 * Hubbard_l(nt) + 1
              n_spin = n_spin + ns(m1,m1,is,na)
           ENDDO
           !
           DO m1 = 1, 2 * Hubbard_l(nt) + 1
              !
              ! Hubbard potential: DC contribution  
              !
              v_hub(m1,m1,is,na) = v_hub(m1,m1,is,na) + Hubbard_J(1,nt)*n_spin + &
                         0.5d0*(Hubbard_U(nt)-Hubbard_J(1,nt)) - Hubbard_U(nt)*n_tot
              !
              ! +U contributions 
              !
              DO m2 = 1, 2 * Hubbard_l(nt) + 1
                 DO m3 = 1, 2 * Hubbard_l(nt) + 1
                    DO m4 = 1, 2 * Hubbard_l(nt) + 1
                       !
                       DO is1 = 1, nspin
                          v_hub(m1,m2,is,na) = v_hub(m1,m2,is,na) + (MOD(nspin,2)+1) * &
                                               u_matrix(m1,m3,m2,m4) * ns(m3,m4,is1,na)
                       ENDDO
                       !
                       v_hub(m1,m2,is,na) = v_hub(m1,m2,is,na) - &
                                            u_matrix(m1,m3,m4,m2) * ns(m3,m4,is,na)
                       !
                       eth_u = eth_u + 0.5d0*( ( u_matrix(m1,m2,m3,m4)-u_matrix(m1,m2,m4,m3) ) * &
                                       ns(m1,m3,is,na)*ns(m2,m4,is,na)+u_matrix(m1,m2,m3,m4)   * &
                                       ns(m1,m3,is,na)*ns(m2,m4,nspin+1-is,na) )
                    ENDDO ! m4
                 ENDDO ! m3
              ENDDO ! m2
              !
           ENDDO ! m1
           !
        ENDDO ! is
        !
     ENDIF
     !
  ENDDO ! na
  ! 
  IF (nspin==1) eth_u = 2.d0 * eth_u
  eth = eth_u - eth_dc
  !
  ! Hubbard energy
  !
  IF ( iverbosity > 0 ) THEN
     WRITE(stdout,'(/5x,"HUBBARD ENERGIES (dc, U, total) ",3f9.4,1x," (Ry)")') eth_dc, eth_u, eth
  ENDIF
  !
  DEALLOCATE (u_matrix)
  !
  RETURN
  !
END SUBROUTINE v_hubbard_full
!---------------------------------------------------------------

!---------------------------------------------------------------
SUBROUTINE v_hubbard_full_nc( ns, v_hub, eth )
  !-------------------------------------------------------------
  !
  !! Computes Hubbard potential and Hubbard energy (noncollinear case).
  !! DFT+U(+J) : Formulation by Liechtenstein et al., Phys. Rev. B 52, R5467 (1995).
  !
  USE kinds,                ONLY : DP
  USE ions_base,            ONLY : nat, ityp
  USE noncollin_module,     ONLY : noncolin
  USE ldaU,                 ONLY : Hubbard_lmax, Hubbard_l, &
                                   Hubbard_U, Hubbard_J, Hubbard_alpha
  USE lsda_mod,             ONLY : nspin
  USE control_flags,        ONLY : iverbosity
  USE io_global,            ONLY : stdout
  !
  IMPLICIT NONE
  !
  COMPLEX(DP) :: ns(2*Hubbard_lmax+1,2*Hubbard_lmax+1,nspin,nat)
  !! Occupation matrix
  COMPLEX(DP) :: v_hub(2*Hubbard_lmax+1,2*Hubbard_lmax+1,nspin,nat)
  !! Hubbard potential
  REAL(DP) :: eth
  !! Hubbard matrix
  !
  !  ... local variables
  !
  REAL(DP) :: eth_dc, eth_noflip, eth_flip, mx, my, mz, mag2
  INTEGER :: is, is1, js, i, j, na, nt, m1, m2, m3, m4
  COMPLEX(DP) :: n_tot, n_aux
  REAL(DP), ALLOCATABLE :: u_matrix(:,:,:,:)
  !
  ALLOCATE( u_matrix(2*Hubbard_lmax+1, 2*Hubbard_lmax+1, 2*Hubbard_lmax+1, 2*Hubbard_lmax+1) )
  !
  eth        = 0.d0
  eth_dc     = 0.d0  
  eth_noflip = 0.d0  
  eth_flip   = 0.d0  
  !
  v_hub(:,:,:,:) = 0.d0
  !
  DO na = 1, nat  
     !
     nt = ityp (na)  
     !
     IF (Hubbard_U(nt) /= 0.d0) THEN  
        !
        ! Initialize U(m1,m2,m3,m4) matrix 
        !
        CALL hubbard_matrix( Hubbard_lmax, Hubbard_l(nt), Hubbard_U(nt), &
                             Hubbard_J(1,nt), u_matrix )
        !
        ! Total N and M^2 for DC (double counting) term
        !
        n_tot = 0.d0
        mx    = 0.d0
        my    = 0.d0
        mz    = 0.d0
        DO m1 = 1, 2*Hubbard_l(nt)+1
          n_tot = n_tot + ns(m1,m1,1,na) + ns(m1,m1,4,na)
          mx = mx + DBLE( ns(m1, m1, 2, na) + ns(m1, m1, 3, na) )
          my = my + 2.d0 * AIMAG( ns(m1, m1, 2, na) )
          mz = mz + DBLE( ns(m1, m1, 1, na) - ns(m1, m1, 4, na) )
        ENDDO  
        mag2 = mx**2 + my**2 + mz**2  
        !
        ! Hubbard energy: DC term
        !
        mx = REAL(n_tot)
        eth_dc = eth_dc + 0.5d0*( Hubbard_U(nt)*mx*(mx-1.d0) - &
                                  Hubbard_J(1,nt)*mx*(0.5d0*mx-1.d0) - &
                                  0.5d0*Hubbard_J(1,nt)*mag2 )   
        !
        DO is = 1, nspin  
           !
           IF (is == 2) THEN
            is1 = 3
           ELSEIF (is == 3) THEN
            is1 = 2
           ELSE
            is1 = is
           ENDIF
           !
           ! Hubbard energy:
           !
           IF (is1 == is) THEN
             !
             ! Non spin-flip contribution
             !
             DO m1 = 1, 2*Hubbard_l(nt)+1
              DO m2 = 1, 2*Hubbard_l(nt)+1
                DO m3 = 1, 2*Hubbard_l(nt)+1
                 DO m4 = 1, 2*Hubbard_l(nt)+1
                   eth_noflip = eth_noflip + 0.5d0*(                            &
                              ( u_matrix(m1,m2,m3,m4)-u_matrix(m1,m2,m4,m3) )*  & 
                              ns(m1,m3,is,na)*ns(m2,m4,is,na) +                 &
                      u_matrix(m1,m2,m3,m4)*ns(m1,m3,is,na)*ns(m2,m4,nspin+1-is,na) )
                 ENDDO
                ENDDO
              ENDDO
             ENDDO
             !
           ELSE
             ! 
             ! Spin-flip contribution
             !
             DO m1 = 1, 2*Hubbard_l(nt)+1
              DO m2 = 1, 2*Hubbard_l(nt)+1
               DO m3 = 1, 2*Hubbard_l(nt)+1
                DO m4 = 1, 2*Hubbard_l(nt)+1
                   eth_flip = eth_flip - 0.5d0*u_matrix(m1,m2,m4,m3)* &
                                     ns(m1,m3,is,na)*ns(m2,m4,is1,na) 
                ENDDO
               ENDDO
              ENDDO
             ENDDO
             !
           ENDIF
           !
           ! Hubbard potential: non spin-flip contribution 
           !
           IF (is1 == is) THEN
             !
             DO m1 = 1, 2*Hubbard_l(nt)+1
              DO m2 = 1, 2*Hubbard_l(nt)+1
               DO m3 = 1, 2*Hubbard_l(nt)+1
                DO m4 = 1, 2*Hubbard_l(nt)+1
                  v_hub(m1,m2,is,na) = v_hub(m1,m2,is,na) + &
                    u_matrix(m1,m3,m2,m4)*( ns(m3,m4,1,na)+ns(m3,m4,4,na) ) 
                ENDDO 
               ENDDO
              ENDDO
             ENDDO
             !
           ENDIF
           !
           ! n_aux = /sum_{i} n_{i,i}^{sigma2, sigma1} for DC term
           !
           n_aux = 0.d0
           DO m1 = 1, 2*Hubbard_l(nt)+1
              n_aux = n_aux + ns(m1,m1,is1,na)  
           ENDDO
           !
           DO m1 = 1, 2*Hubbard_l(nt)+1
             ! 
             ! Hubbard potential: DC contribution  
             !
             v_hub(m1,m1,is,na) = v_hub(m1,m1,is,na) + Hubbard_J(1,nt)*n_aux
             !
             IF (is1 == is) THEN
                v_hub(m1,m1,is,na) = v_hub(m1,m1,is,na) + &
                     0.5d0*(Hubbard_U(nt)-Hubbard_J(1,nt)) - Hubbard_U(nt)*n_tot  
             ENDIF
             !
             ! Hubbard potential: spin-flip contribution
             !
             DO m2 = 1, 2*Hubbard_l(nt)+1  
              DO m3 = 1, 2*Hubbard_l(nt)+1
               DO m4 = 1, 2*Hubbard_l(nt)+1
                  v_hub(m1,m2,is,na) = v_hub(m1,m2,is,na) - &
                             u_matrix(m1,m3,m4,m2) * ns(m3,m4,is1,na) 
               ENDDO 
              ENDDO
             ENDDO
             !
           ENDDO
           !
        ENDDO ! is
        !
     ENDIF
     !
  ENDDO ! na
  !
  eth = eth_noflip + eth_flip - eth_dc
  !
  ! Hubbard energies
  !
  IF ( iverbosity > 0 ) THEN
    WRITE(stdout,*) '--- in v_hubbard ---'
    WRITE(stdout,'("Hub. E (dc, noflip, flip, total) ",4f9.4)') &
                                 eth_dc, eth_noflip, eth_flip, eth 
    WRITE(stdout,*) '-------'
  ENDIF
  !
  DEALLOCATE (u_matrix)
  !
  RETURN
  !
END SUBROUTINE v_hubbard_full_nc
!----------------------------------------------------------------------------

!------------------------------------------------------------------------------------
SUBROUTINE v_hubbard_extended (nsg, v_hub, eth)
  !-----------------------------------------------------------------------------------
  !
  !! Computes extended Hubbard potential and Hubbard energy.
  !! DFT+U+V: Simplified rotationally-invariant formulation by
  !! V.L. Campo Jr and M. Cococcioni, J. Phys.: Condens. Matter 22, 055602 (2010).
  !
  USE kinds,         ONLY : DP
  USE ions_base,     ONLY : nat, ityp
  USE ldaU,          ONLY : Hubbard_l, Hubbard_alpha, Hubbard_J0, Hubbard_beta,&
                            ldim_u, ldmx_tot, max_num_neighbors, at_sc,        &
                            neighood, Hubbard_V, Hubbard_alpha_back,           &
                            is_hubbard, is_hubbard_back, dfpt_hub
  USE lsda_mod,      ONLY : nspin
  USE control_flags, ONLY : iverbosity
  USE io_global,     ONLY : stdout
  !
  IMPLICIT NONE
  !
  COMPLEX(DP), INTENT(IN)  :: nsg  (ldmx_tot, ldmx_tot, max_num_neighbors, nat, nspin)
  COMPLEX(DP), INTENT(OUT) :: v_hub(ldmx_tot, ldmx_tot, max_num_neighbors, nat, nspin)
  REAL(DP),    INTENT(OUT) :: eth
  ! 
  ! Local variables 
  !
  INTEGER :: is, isop, na, na1, na2, nt, nt1, nt2, m1, m2, viz, equiv_na2, i_type
  COMPLEX(DP) :: check, check_en
  INTEGER, EXTERNAL :: type_interaction, find_viz
  !
  eth  = 0.d0
  v_hub(:,:,:,:,:) = (0.d0, 0.d0)
  check = (0.d0, 0.d0)
  check_en = (0.d0, 0.d0)
  !
  DO na1 = 1, nat
     !
     nt1 = ityp(na1)
     !
     IF ( is_hubbard(nt1) .OR. is_hubbard_back(nt1) ) THEN
        !
        DO is = 1, nspin
           !
           DO viz = 1, neighood(na1)%num_neigh
              !
              na2 = neighood(na1)%neigh(viz)
              equiv_na2 = at_sc(na2)%at
              nt2 = ityp(equiv_na2)
              !
              IF ((is_hubbard(nt2).OR.is_hubbard_back(nt2)) .AND. &
                  (Hubbard_V(na1,na2,1).NE.0.d0 .OR. &
                   Hubbard_V(na1,na2,2).NE.0.d0 .OR. &
                   Hubbard_V(na1,na2,3).NE.0.d0 .OR. &
                   Hubbard_V(na1,na2,4).NE.0.d0) ) THEN
                  !
                  ! For both standard and background states of a center atom
                  DO m1 = 1, ldim_u(nt1)
                     ! For both standard and background states of the neighbor atom
                     DO m2 = 1, ldim_u(nt2)
                        !
                        i_type = type_interaction(na1,m1,equiv_na2,m2)
                        !
                        v_hub(m2,m1,viz,na1,is) = &
                           - CONJG(nsg(m2,m1,viz,na1,is)) * Hubbard_V(na1,na2,i_type)
                        check = check  - CONJG(nsg(m2,m1,viz,na1,is)) * Hubbard_V(na1,na2,i_type)

                        !
                        eth = eth - nsg(m2,m1,viz,na1,is) * CONJG(nsg(m2,m1,viz,na1,is)) &
                                    * Hubbard_V(na1,na2,i_type) * 0.5d0
                        check_en = check_en -nsg(m2,m1,viz,na1,is) * CONJG(nsg(m2,m1,viz,na1,is)) &
                                    * Hubbard_V(na1,na2,i_type) * 0.5d0
                        !
                     ENDDO
                  ENDDO
                  !
                  IF ( na1.EQ.na2 ) THEN
                     !
                     na = find_viz(na1,na1)
                     !
                     ! This is the diagonal term (like in the DFT+U only case)
                     ! 
                     DO m1 = 1, ldim_u(nt1)
                        !
                        i_type = type_interaction(na1,m1,equiv_na2,m1)
                        !
                        v_hub(m1,m1,na,na1,is) = v_hub(m1,m1,na,na1,is) &
                                       + Hubbard_V(na1,na1,i_type) * 0.5d0
                        check = check  + Hubbard_V(na1,na1,i_type) * 0.5d0
                        ! 
                        eth = eth + nsg(m1,m1,na,na1,is) &
                                       * Hubbard_V(na1,na1,i_type) * 0.5d0
                        check_en = check_en +nsg(m1,m1,na,na1,is) &
                                       * Hubbard_V(na1,na1,i_type) * 0.5d0
                        !
                     ENDDO
                     !
                     ! Hubbard_J0 (only on-site)
                     !
                     IF ( nspin.EQ.2 .AND. &
                          (Hubbard_J0(nt1).NE.0.d0 .OR. Hubbard_beta(nt1).NE.0.d0) ) THEN
                          !
                          IF (is.EQ.1) THEN
                             isop = 2
                          ELSE
                             isop = 1
                          ENDIF
                          !
                          DO m1 = 1, 2*Hubbard_l(nt1)+1
                             !
                             v_hub(m1,m1,na,na1,is) = v_hub(m1,m1,na,na1,is) - &
                                                      Hubbard_J0(nt1) * 0.5d0
                             eth = eth - 0.5d0 * Hubbard_J0(nt1) * nsg(m1,m1,na,na1,is)
                             !
                             IF (is.EQ.1) THEN
                                v_hub(m1,m1,na,na1,is) = v_hub(m1,m1,na,na1,is) + &
                                                         Hubbard_beta(nt1)
                                eth = eth + Hubbard_beta(nt1) * nsg(m1,m1,na,na1,is)
                             ELSE
                                v_hub(m1,m1,na,na1,is) = v_hub(m1,m1,na,na1,is) - &
                                                         Hubbard_beta(nt1)
                                eth = eth - Hubbard_beta(nt1) * nsg(m1,m1,na,na1,is)
                             ENDIF
                             !
                             DO m2 = 1, 2*Hubbard_l(nt1)+1
                                v_hub(m2,m1,na,na1,is) = v_hub(m2,m1,na,na1,is) + &
                                   CONJG(nsg(m2,m1,na,na1,is) + nsg(m2,m1,na,na1,isop)) * &
                                   Hubbard_J0(nt1)
                                eth = eth + nsg(m2,m1,na,na1,is) *    &
                                   CONJG(nsg(m2,m1,na,na1,is) + nsg(m2,m1,na,na1,isop)) * &
                                   Hubbard_J0(nt1) * 0.5d0
                             ENDDO
                             !
                          ENDDO ! m1
                          !
                     ENDIF
                     !
                  ENDIF
                  !
              ENDIF
              !
           ENDDO ! viz
           !
        ENDDO ! is
        !  
     ENDIF
     !
     ! Hubbard_alpha or Hubbard_alpha_back
     !
     IF ( ldim_u(nt1).GT.0 .AND. &
          (Hubbard_alpha(nt1).NE.0.0d0 .OR. Hubbard_alpha_back(nt1).NE.0.0d0) ) THEN
          !
          na = find_viz(na1,na1)
          !
          DO is = 1, nspin
             !
             DO m1 = 1, 2*Hubbard_l(nt1)+1
                v_hub(m1,m1,na,na1,is) = v_hub(m1,m1,na,na1,is) + Hubbard_alpha(nt1)
                eth = eth + nsg(m1,m1,na,na1,is)* Hubbard_alpha(nt1)
             ENDDO
             !
             IF ( ldim_u(nt1).GT.2*Hubbard_l(nt1)+1 .AND. &
                  Hubbard_alpha_back(nt1).NE.0.0d0 ) THEN
                ! Background states
                DO m1 = 2*Hubbard_l(nt1)+2, ldim_u(nt1)
                   v_hub(m1,m1,na,na1,is) = v_hub(m1,m1,na,na1,is) + Hubbard_alpha_back(nt1)
                   eth = eth + nsg(m1,m1,na,na1,is)* Hubbard_alpha_back(nt1)
                ENDDO
             ENDIF
             !
          ENDDO
          !
     ENDIF
     !
  ENDDO ! na1
  !
  IF (nspin.EQ.1) eth = eth * 2.d0
  !
  ! Hubbard energy
  !
  IF ( iverbosity > 0 .AND. .NOT.dfpt_hub ) THEN
     WRITE(stdout,'(/5x,"HUBBARD ENERGY = ",f9.4,1x," (Ry)")') eth
     !write(stdout,*) "check col UV",  check
     !write(stdout,*) "check_en col UV",  check_en
  ENDIF
  !
  RETURN
  !
END SUBROUTINE v_hubbard_extended
!---------------------------------------------------------------------

SUBROUTINE v_hubbard_extended_nc (nsg, v_hub, eth)
   !-----------------------------------------------------------------------------------
   !
   !! Computes extended Hubbard potential and Hubbard energy.
   !! DFT+U+V: Simplified rotationally-invariant formulation by
   !! V.L. Campo Jr and M. Cococcioni, J. Phys.: Condens. Matter 22, 055602 (2010).
   !
   USE kinds,             ONLY : DP
   USE ions_base,         ONLY : nat, ityp
   USE ldaU,              ONLY : Hubbard_l, Hubbard_alpha, Hubbard_J0, Hubbard_beta,   &
                                 ldim_u, ldmx_tot, max_num_neighbors, at_sc, neighood, &
                                 Hubbard_V, Hubbard_alpha_back, is_hubbard, &
                                 is_hubbard_back, dfpt_hub
   USE lsda_mod,          ONLY : nspin
   USE control_flags,     ONLY : iverbosity
   USE io_global,         ONLY : stdout
   USE noncollin_module,  ONLY : npol
   !
   IMPLICIT NONE
   !
   COMPLEX(DP), INTENT(IN)  :: nsg  (ldmx_tot, ldmx_tot, max_num_neighbors, nat, nspin)
   COMPLEX(DP), INTENT(OUT) :: v_hub(ldmx_tot, ldmx_tot, max_num_neighbors, nat, nspin)
   REAL(DP),    INTENT(OUT) :: eth
   COMPLEX(DP) :: check, check_en
   ! 
   ! Local variables 
   !
   INTEGER :: is, is1, isop, na, na1, na2, nt, nt1, nt2, m1, m2, viz, equiv_na2
   INTEGER, EXTERNAL :: find_viz
   !
   eth  = 0.d0
   v_hub(:,:,:,:,:) = (0.d0, 0.d0)
   check = (0.0,0.0)
   check_en = (0.0,0.0)
   !
   !write(stdout,*) nsg
   DO na1 = 1, nat
      !
      nt1 = ityp(na1)
      !
      IF ( is_hubbard(nt1) ) THEN
         !
         DO is = 1, nspin
            !
            IF (is == 2) THEN
               is1 = 3
            ELSEIF (is == 3) THEN
               is1 = 2
            ELSE
               is1 = is
            ENDIF
            DO viz = 1, neighood(na1)%num_neigh
               !
               na2 = neighood(na1)%neigh(viz)
               equiv_na2 = at_sc(na2)%at
               nt2 = ityp(equiv_na2)
               !
               IF (is_hubbard(nt2) .AND. &
                   (Hubbard_V(na1,na2,1).NE.0.d0) ) THEN
                   !
                   ! Here no need to use is1: complex conjugation is enough
                   ! For both standard and background states of a center atom
                   DO m1 = 1, ldim_u(nt1)
                      ! For both standard and background states of the neighbor atom
                      DO m2 = 1, ldim_u(nt2)
                         !
                         v_hub(m2,m1,viz,na1,is) = - CONJG(nsg(m2,m1,viz,na1,is)) * Hubbard_V(na1,na2,1)
                        check = check - CONJG(nsg(m2,m1,viz,na1,is)) * Hubbard_V(na1,na2,1)
                         !
                         eth = eth - nsg(m2,m1,viz,na1,is) * CONJG(nsg(m2,m1,viz,na1,is)) &
                                     * Hubbard_V(na1,na2,1) * 0.5d0
                        check_en = check_en - nsg(m2,m1,viz,na1,is) * CONJG(nsg(m2,m1,viz,na1,is)) &
                                     * Hubbard_V(na1,na2,1) * 0.5d0
                         !
                      ENDDO
                   ENDDO
                   !
                   IF ( na1.EQ.na2 .AND. is1.EQ.is) THEN
                      !
                      na = find_viz(na1,na1)
                      !
                      ! This is the diagonal term (like in the DFT+U only case)
                      ! 
                      DO m1 = 1, ldim_u(nt1)
                         !
                         v_hub(m1,m1,na,na1,is) = v_hub(m1,m1,na,na1,is) &
                                        + Hubbard_V(na1,na1,1) * 0.5d0
                        check = check + Hubbard_V(na1,na1,1) * 0.5d0
                         ! 
                         eth = eth + nsg(m1,m1,na,na1,is) &
                                        * Hubbard_V(na1,na1,1) * 0.5d0
                        check_en = check_en + nsg(m1,m1,na,na1,is) &
                                        * Hubbard_V(na1,na1,1) * 0.5d0
                         !
                      ENDDO
                      !
                   ENDIF
                   !
               ENDIF
               !
            ENDDO ! viz
            !
         ENDDO ! is
         !  
      ENDIF
      !
      ! Hubbard_alpha
      !!
      IF ( ldim_u(nt1).GT.0 .AND. (Hubbard_alpha(nt1).NE.0.0d0 ) ) THEN
         !
         na = find_viz(na1,na1)
         !
         DO is = 1, npol
            !
            DO m1 = 1, 2*Hubbard_l(nt1)+1
               v_hub(m1,m1,na,na1,is**2) = v_hub(m1,m1,na,na1,is**2) + Hubbard_alpha(nt1)
               eth = eth + nsg(m1,m1,na,na1,is**2)* Hubbard_alpha(nt1)
            ENDDO
            !
         ENDDO
         !
      ENDIF
      !
   ENDDO ! na1
   !
   !
   IF (nspin.EQ.1) eth = eth * 2.d0
   !
   ! Hubbard energy
   !
   IF ( iverbosity > 0 .AND. .NOT.dfpt_hub ) THEN
      WRITE(stdout,'(/5x,"HUBBARD ENERGY = ",f9.4,1x," (Ry)")') eth
   ENDIF
   !
   RETURN
   !
 END SUBROUTINE v_hubbard_extended_nc
!
!-----------------------------------------------------------------------
SUBROUTINE v_hubbard_resolved( ns, v_hub, eth )
!---------------------------------------------------------------------
!
!! Computes Hubbard potential and Hubbard energy
!! for a manifold of selected spin- or magnetic quantum orbitals.
!! The Hubbard potential is first calculated in the diagonal representation
!! based on the eigenvalues of the occupation matrix and then 
!! re-rotated using the eigenvectors in order to retain compatiblity with
!! other parts of the code.
!! See Macke et al., arXiv:2312.13580 (2023).
!! Uses the simplified rotationally-invariant formulation by
!! Dudarev et al., Phys. Rev. B 57, 1505 (1998).
!
USE kinds,                ONLY : DP
USE ions_base,            ONLY : nat, ityp
USE ldaU,                 ONLY : Hubbard_lmax, Hubbard_l, Hubbard_Um, &
                                 Hubbard_alpha_m, lambda_ns, dfpt_hub, &
                                 eigenvecs_ref, order_um, apply_U, hub_pot_fix
USE lsda_mod,             ONLY : nspin
USE constants,            ONLY : eps16, RYTOEV
USE control_flags,        ONLY : iverbosity
USE io_global,            ONLY : stdout
!
IMPLICIT NONE
!
REAL(DP), INTENT(IN)  :: ns(2*Hubbard_lmax+1,2*Hubbard_lmax+1,nspin,nat)
!! occupation matrix
REAL(DP), INTENT(OUT) :: v_hub(2*Hubbard_lmax+1,2*Hubbard_lmax+1,nspin,nat)
!! Hubbard potential
REAL(DP), INTENT(OUT) :: eth
!! Hubbard energy
!
!  ... local variables
!
! Hubbard potential in the diagonal representation
REAL(DP)                 :: v_hub_diag(2*Hubbard_lmax+1,nspin,nat), temp
REAL(DP)                 :: effU, effalpha
! eigenvectors of the ns occupation matrix in the current iteration
COMPLEX(DP)              :: eigenvecs_current(2*Hubbard_lmax+1,2*Hubbard_lmax+1,nspin)
!
INTEGER                  :: is, na, nt, m1, m2, m3, ldim, m_order
! the ordering vector for the orbital-tracking routine
INTEGER                  :: order(2*Hubbard_lmax+1)
LOGICAL                  :: is_first
IF (.NOT. ALLOCATED(order_um)) THEN 
  IF (nspin == 2 ) THEN 
     ALLOCATE(order_um(2*Hubbard_lmax+1,nspin, nat))  
  ELSE 
     ALLOCATE(order_um(2*Hubbard_lmax+1,1,nat)) 
  END IF 
  order_um = 0
END IF 

!
!
eth    = 0.d0
lambda_ns(:,:,:) = 0.d0
! orbital occupations (=eigenvalues of rho%ns)
!
v_hub(:,:,:,:) = 0.d0
v_hub_diag(:,:,:) = 0.d0 
!
IF ( .NOT. apply_U ) RETURN
! Do not apply corrections before the eigenstates have stabilized
! The eigenstates are considered stable a) when starting from a
! converged charge density or b) when the treshold of iterative 
! diagonalization gets small enough (typically within 3-6 iterations).
! This is controlled by electrons.f90.
! 
DO na = 1, nat
   !
   nt = ityp (na)
   !
   IF ( ANY(ABS(Hubbard_Um(:,:,nt)) .GT. eps16) .OR. &
        ANY(ABS(Hubbard_alpha_m(:,:,nt)) .GT. eps16) ) THEN       
      !
      ldim = 2 * Hubbard_l(nt) + 1
      eigenvecs_current(:,:,:) = CMPLX(0.d0,0.d0, kind=dp)
      is_first = ALL(eigenvecs_ref(:,:,:,na) == eigenvecs_current)  

      !
      effU = 0.0
      effalpha = 0.0
      !
      CALL diag_ns( ldim, ns(1:ldim,1:ldim,:,na), lambda_ns(1:ldim,:,na), &
                    eigenvecs_current(1:ldim,1:ldim,:) )
      !
      DO is = 1, nspin
         !
         ! sort eigenvectors with respect to the (reference) order established in eigvecs_first
         IF (is_first)  THEN  
             order(1:ldim) = order_um(1:ldim,is,na)
             IF (ALL(order(1:ldim)==0)) order(1:ldim) = [(m1,m1=1,ldim)] 
             DO m1 =1, ldim 
               eigenvecs_ref(1:ldim, order(m1), is, na) = eigenvecs_current(1:ldim, m1, is) 
             END DO 
         END IF 
         order(:) = 0
         CALL order_eigenvecs( order(1:ldim), eigenvecs_current(1:ldim,1:ldim,is), &
                                 eigenvecs_ref(1:ldim,1:ldim,is,na), ldim )
         !
         
         order_um(1:ldim,is,na) = order(1:ldim) 
         DO m1 = 1, ldim
            !
            ! calculate Hubbard potential and -energy
            ! using the eigenvalues and their ordering
            m_order = order(m1)
            effU = Hubbard_Um(m_order,is,nt)
            effalpha = Hubbard_alpha_m(m_order,is,nt)
            !
            ! linear Hubbard U terms:
            eth = eth + ( effalpha + 0.5D0*effU ) * lambda_ns(m1,is,na)
            v_hub_diag(m1,is,na) = v_hub_diag(m1,is,na) + &
                                       effalpha + 0.5D0*effU
            ! quadratic Hubbard U terms:
            eth = eth - 0.5D0 * effU * lambda_ns(m1,is,na) * lambda_ns(m1,is,na)
            v_hub_diag(m1,is,na) = v_hub_diag(m1,is,na) - &
                                       effU * lambda_ns(m1,is,na)
            !
            ! The following can be eliminated once code is approved
#if defined(__DEBUG)
            WRITE( stdout,'(/5x,"v_of_rho:")')
            WRITE( stdout,'(/5x,"AT:",i2," ,m_order:",i2,", is: ",i2,", LAMBDA:",f7.5," U:",f7.5)') &
                  na,m_order,is,lambda_ns(m1,is,na),effU*RYTOEV
            IF ( effalpha /= 0.d0) &
               WRITE( stdout,'(/5x,"m_order:",i2,", is: ",i2,", LAMBDA:",f8.5," alpha:",f8.5)') &
                  m_order,is,lambda_ns(m1,is,na),effalpha*RYTOEV
            WRITE( stdout,'(/5x,"Hubbard potential of this state (v_hub_diag) = ",f8.5)') &
                  v_hub_diag(m1,is,na)
            WRITE( stdout,'(/5x,"Cumulative running total of Hubbard energies (eth) = ",f8.5)') eth
#endif
            !
         ENDDO
         !
         ! backrotation of v_hub_diag to the non-diagonal v_hub
         DO m1 = 1, ldim
            DO m2 = 1, ldim
               temp = CMPLX(0.d0,0.d0, kind=dp)
               DO m3 = 1, ldim
                  temp = temp + CONJG(eigenvecs_current(m1,m3,is))* &
                     v_hub_diag(m3,is,na)*eigenvecs_current(m2,m3,is)
               ENDDO
               v_hub(m1,m2,is,na) = DBLE(temp)
            ENDDO
         ENDDO
      ENDDO ! is
      !
      IF ( ALL(eigenvecs_ref(:,:,:,na) .EQ. 0.d0) ) THEN
         !
         ! if this routine is executed for the first time,
         ! save the current eigenvecs as reference eigenvecs
         eigenvecs_ref(1:ldim,1:ldim,:,na) = eigenvecs_current(1:ldim,1:ldim,:)
      ENDIF
      !
   ENDIF
   !
ENDDO ! nt
!
IF (nspin==1) eth = 2.d0 * eth
!
! print Hubbard energy
!
IF ( iverbosity > 0 .AND. .NOT.dfpt_hub ) THEN
   WRITE(stdout,'(/5x,"HUBBARD ENERGY = ",f9.5,1x," (Ry)")') eth
ENDIF
!
RETURN
!
END SUBROUTINE v_hubbard_resolved
!----------------------------------------------------------------------------
SUBROUTINE v_hubbard_resolved_nc( ns, v_hub, eth )
!----------------------------------------------------------------------------
!
!! Computes Hubbard potential and Hubbard energy
!! for a manifold of selected orbitals (noncollinear formulation).
!! The Hubbard potential is first calculated in the diagonal representation
!! based on the eigenvalues of the occupation matrix and then 
!! re-rotated using the eigenvectors in order to retain compatiblity with
!! other parts of the code.
!! Uses the simplified rotationally-invariant formulation by
!! Dudarev et al., Phys. Rev. B 57, 1505 (1998).
!
USE kinds,                ONLY : DP
USE ions_base,            ONLY : nat, ityp
USE ldaU,                 ONLY : Hubbard_lmax, Hubbard_l, Hubbard_Um_nc, &
                                 Hubbard_alpha_m_nc, lambda_ns, order_um,&
                                 eigenvecs_ref, apply_U, hub_pot_fix, dfpt_hub
USE lsda_mod,             ONLY : nspin
USE constants,            ONLY : eps16, RYTOEV
USE control_flags,        ONLY : iverbosity
USE io_global,            ONLY : stdout
!
IMPLICIT NONE
!
COMPLEX(DP), INTENT(IN) :: ns(2*Hubbard_lmax+1,2*Hubbard_lmax+1,nspin,nat)
!! occupation matrix
COMPLEX(DP), INTENT(OUT):: v_hub(2*Hubbard_lmax+1,2*Hubbard_lmax+1,nspin,nat)
!! Hubbard potential
REAL(DP), INTENT(OUT) :: eth
!! Hubbard energy
!
!  ... local variables
!
REAL(DP)                 :: effU, effalpha
COMPLEX(DP)              :: v_hub_diag(4*Hubbard_lmax+2,nat), temp
COMPLEX(DP)              :: v_hub_temp(4*Hubbard_lmax+2,4*Hubbard_lmax+2)
COMPLEX(DP)              :: eigenvecs_current(4*Hubbard_lmax+2,4*Hubbard_lmax+2)
!
INTEGER                  :: is, na, nt, m1, m2, m3, ldim, m_order
!
INTEGER                  :: order(4*Hubbard_lmax+2) 
LOGICAL                  :: is_first
!! ordering vector
!
IF (.NOT. ALLOCATED(order_um)) THEN 
  ALLOCATE(order_um(4*Hubbard_lmax+2,1,nat)) 
  order_um = 0
END IF 

eth    = 0.d0
lambda_ns(:,:,:) = 0.d0
! orbital occupations (=eigenvalues of rho%ns)
!
v_hub(:,:,:,:) = 0.d0
v_hub_diag(:,:) = 0.d0 
! Hubbard potential in the diagonal representation
!
IF ( .NOT. apply_U ) RETURN
! Do not apply corrections before the eigenstates have stabilized
! The eigenstates are considered stable a) when starting from a
! converged charge density or b) when the treshold of iterative 
! diagonalization gets small enough (typically within 3-6 iterations)
! this is checked in electrons.f90
! 
DO na = 1, nat
   !
   nt = ityp (na)
   !
   IF ( ANY(ABS(Hubbard_Um_nc(:,nt)) .GT. eps16) .OR. &
        ANY(ABS(Hubbard_alpha_m_nc(:,nt)) .GT. eps16) ) THEN       
      !
      ldim = 2 * Hubbard_l(nt) + 1
      eigenvecs_current(:,:) = CMPLX(0.d0,0.d0, kind=dp)
      is_first = ALL(eigenvecs_ref(:,1,:,na) == eigenvecs_current)  
      !
      effU = 0.0
      effalpha = 0.0
      !
      CALL diag_ns_nc( ldim, ns(1:ldim,1:ldim,:,na), lambda_ns(1:2*ldim,1,na), &
                     eigenvecs_current(1:2*ldim,1:2*ldim) )
                              !
      ! sort eigenvectors with respect to the (reference) order established in eigvecs_first
      IF (is_first)  THEN  
        order(1:2*ldim) = order_um(1:2*ldim,1,na)
        IF (ALL(order(1:2*ldim)==0)) order(1:2*ldim) = [(m1,m1=1,2*ldim)] 
        DO m1 =1, 2*ldim 
          eigenvecs_ref(1:2*ldim, order(m1), 1, na) = eigenvecs_current(1:2*ldim, m1) 
        END DO 
      END IF 
      order(:) = 0
      !
      CALL order_eigenvecs( order(1:2*ldim), eigenvecs_current(1:2*ldim,1:2*ldim), &
                                 eigenvecs_ref(1:2*ldim,1:2*ldim,1,na), 2*ldim )
      !
      ! No need to iterate over is (all done in diag_ns_nc)
      order_um(1:2*ldim,1,na) = order(1:2*ldim) 
      DO m1 = 1, 2*ldim
         !
         ! calculate Hubbard potential and -energy
         ! using the eigenvalues and their ordering
         m_order = order(m1)
         effU = Hubbard_Um_nc(m_order,nt)
         effalpha = Hubbard_alpha_m_nc(m_order,nt)
         !
         ! linear Hubbard U terms:
         eth = eth + ( effalpha + 0.5D0*effU ) * lambda_ns(m1,1,na)
         v_hub_diag(m1,na) = v_hub_diag(m1,na) + &
                                    effalpha + 0.5D0*effU
         ! quadratic Hubbard U terms:
         eth = eth - 0.5D0 * effU * lambda_ns(m1,1,na) * lambda_ns(m1,1,na)
         v_hub_diag(m1,na) = v_hub_diag(m1,na) - &
                                    effU * lambda_ns(m1,1,na)
         !
         ! The following block can be eliminated once code is approved
#if defined(__DEBUG)
         WRITE( stdout,'(/5x,"v_of_rho (NC):")')
         WRITE( stdout,'(/5x,"AT:",i2," ,m_order:",i2,", LAMBDA:",f7.5," U:",f7.5)') &
               na,m_order,lambda_ns(m1,1,na),effU*RYTOEV
         IF ( effalpha /= 0.d0) &
            WRITE( stdout,'(/5x,"m_order:",i2,", LAMBDA:",f8.5," alpha:",f8.5)') &
               m_order,lambda_ns(m1,1,na),effalpha*RYTOEV
         WRITE( stdout,'(/5x,"Hubbard potential of this state (v_hub_diag) = ",f8.5)') &
               v_hub_diag(m1,na)
         WRITE( stdout,'(/5x,"Cumulative running total of Hubbard energies (eth) = ",f8.5)') eth
#endif
         !
      ENDDO
         !
      ! backrotation of v_hub_diag to the non-diagonal v_hub
      ! first, go back to the 2ldim*2ldim matrix
      v_hub_temp(:,:) = 0.d0
      !
      DO m1 = 1, 2*ldim
         DO m2 = 1, 2*ldim
            temp = CMPLX(0.d0,0.d0, kind=dp)
            DO m3 = 1, 2*ldim
               temp = temp + CONJG(eigenvecs_current(m1,m3))* &
                  v_hub_diag(m3,na)*eigenvecs_current(m2,m3)
            ENDDO
            v_hub_temp(m1,m2) = temp
         ENDDO
      ENDDO
      ! now, sort the different quadrants of v_hub_temp into the actual
      ! v_hub. upper left quadrant is spin 1, upper right spin 2 etc.
      !
      DO m1 = 1, ldim
         DO m2 = 1,ldim
            v_hub(m1,m2,1,na) = v_hub_temp(m1,m2)
            v_hub(m1,m2,2,na) = v_hub_temp(m1,ldim+m2)
            v_hub(m1,m2,3,na) = v_hub_temp(ldim+m1,m2)
            v_hub(m1,m2,4,na) = v_hub_temp(ldim+m1,ldim+m2)
         ENDDO
      ENDDO
      !
      IF ( ALL(eigenvecs_ref(:,:,1,na) .EQ. 0.d0) ) THEN
         ! if this routine is executed for the first time,
         ! save the current eigenvecs as reference eigenvecs
         eigenvecs_ref(1:2*ldim,1:2*ldim,1,na) = &
                     & eigenvecs_current(1:2*ldim,1:2*ldim)
      ENDIF
      !
   ENDIF
   !
ENDDO ! nt
!
! print Hubbard energy
!
IF ( iverbosity > 0 .AND. .NOT.dfpt_hub ) THEN
   WRITE(stdout,'(/5x,"HUBBARD ENERGY = ",f9.5,1x," (Ry)")') eth
ENDIF
!
RETURN
!
END SUBROUTINE v_hubbard_resolved_nc
!----------------------------------------------------------------------------
