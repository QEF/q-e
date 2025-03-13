!!-------------------------------------------------------------------------
!
MODULE sic_mod
   !
   ! ... Module for self-interaction-corrected calculations
   !     Written by Stefano Falletta
   !
   USE starting_scf,     ONLY : starting_pot
   USE cell_base,        ONLY : omega
   USE ener,             ONLY : esic
   USE fft_base,         ONLY : dfftp, dffts
   USE gvect,            ONLY : ngm
   USE io_global,        ONLY : stdout
   USE kinds,            ONLY : DP
   USE klist,            ONLY : nkstot, degauss, nks
   USE lsda_mod,         ONLY : nspin, isk, current_spin
   USE noncollin_module, ONLY : noncolin
   USE scf,              ONLY : scf_type, create_scf_type, &
                                scf_type_copy, rhoz_or_updw, rho
   USE xc_lib,           ONLY : xclib_dft_is
   USE uspp,             ONLY : okvan
   !
   USE mp_bands,         ONLY : inter_bgrp_comm, intra_bgrp_comm
   USE mp,               ONLY : mp_sum
   USE electrons_base,   ONLY : set_nelup_neldw
   USE klist,            ONLY : nelec, tot_magnetization, nelup, neldw
   USE control_flags,    ONLY : use_gpu, lbfgs
   ! 
   IMPLICIT NONE
   SAVE
   !
   CHARACTER(len=20) :: pol_type      ! 'e' for electron, 'h' for hole
   REAL(DP)          :: sic_gamma     ! prefactor for localizing potential
   LOGICAL           :: sic_energy    ! calculate sic energy
   LOGICAL           :: sic_first     ! calculate sic energy
   INTEGER           :: isp           ! spin with/without the polaron
   TYPE(scf_type), ALLOCATABLE :: rho_n ! neutral density
   INTEGER           :: fp, fn        ! occupations charged and neutral states
   INTEGER           :: qp, qn        ! supercell charges polaron and neutral states
   !
   CONTAINS
   !
   !--------------------
   SUBROUTINE init_sic()
   !--------------------
      !
      !   ... initialization SIC module 
      !
      IMPLICIT NONE
      !
      IF (pol_type /= 'e' .and. pol_type /= 'h') CALL errore('sic_init', 'error in pol_type', 1)
      IF (starting_pot /= 'atomic') CALL errore('sic_init', 'only atomic starting_pot supported', 1)
      IF (degauss /= 0.D0)          CALL errore('sic_init', 'gaussian smearing not allowed', 1)
      IF (nspin .ne. 2)             CALL errore('sic_init', 'spin polarized calculation required',1)
      IF (nkstot < 2)               CALL errore('sic_init', 'error in the value of nkstot',1)
      IF (dffts%has_task_groups)    CALL errore('sic_init', 'task groups not implemented',1)
      IF (noncolin)                 CALL errore('sic_init', 'non-collinear spin calculations not implemented',1)
      IF (okvan)                    CALL errore('sic_init', 'norm-conserving pseudopotentials required',1)
      IF (xclib_dft_is('meta'))     CALL errore('sic_init', 'meta-GGA not implemented',1)
      IF (xclib_dft_is('hybrid'))   CALL errore('sic_init', 'hybrid not implemented',1)
      IF (lbfgs .AND. .NOT. sic_energy) CALL errore('sic_init', 'use damped ion dynamics when sic_energy = .false.',1)
      IF (pol_type == 'e') THEN
         isp = 1
         fp = 1
         fn = 0
      END IF
      IF (pol_type == 'h') THEN
         isp = 2
         fp = 0
         fn = 1
      END IF
      sic_first = .true.
      esic = 0.d0
      !
   END SUBROUTINE init_SIC
   !
   !--------------------------------------------------
   SUBROUTINE add_vsic(rho, rho_core, rhog_core, v)
   !--------------------------------------------------
      !
      ! ... apply self-interaction correction to the KS potential 
      !
      IMPLICIT NONE
      !
      TYPE(scf_type), INTENT(INOUT) :: rho                   ! electron density
      REAL(DP),       INTENT(IN)    :: rho_core(dfftp%nnr)   ! core charge in real space
      COMPLEX(DP),    INTENT(IN)    :: rhog_core(ngm)        ! core charge in reciprocal space
      TYPE(scf_type), INTENT(INOUT) :: v                     ! Hxc potential of rho
      TYPE(scf_type), ALLOCATABLE   :: rho_aux               ! auxiliary density
      REAL(DP),       ALLOCATABLE   :: vxc(:,:)              ! xc potential of rho
      REAL(DP),       ALLOCATABLE   :: vxc_aux(:,:)          ! xc potential of rho_aux
      REAL(DP),       ALLOCATABLE   :: vh_aux(:,:)           ! hartree potential
      REAL(DP) :: vtxc_aux, etxc_aux, eh_aux, charge_aux, sgn
      REAL(DP), PARAMETER :: ry2ev = 13.605698066
      INTEGER :: is
      !
      IF(tot_magnetization /= 0) THEN
         ALLOCATE(rho_aux)
         ALLOCATE(vxc(dfftp%nnr,nspin))
         ALLOCATE(vxc_aux(dfftp%nnr,nspin))
         ALLOCATE(vh_aux(dfftp%nnr,nspin))
         !
         IF(pol_type == 'e') sgn = +1
         IF(pol_type == 'h') sgn = -1
         ! 
         ! ... rho_aux : density without polaron
         !
         CALL create_scf_type(rho_aux)
         CALL scf_type_copy(rho, rho_aux)
         CALL rhoz_or_updw(rho_aux, 'r_and_g', '->updw')
         rho_aux%of_r(:,isp) = rho_aux%of_r(:,isp) - sgn*rho%pol_r(:,1)  
         rho_aux%of_g(:,isp) = rho_aux%of_g(:,isp) - sgn*rho%pol_g(:,1)
         CALL rhoz_or_updw(rho_aux, 'r_and_g', '->rhoz')
         !
         ! ... add local potential to KS Hamiltonian
         !
         vxc(:,:)     = 0.d0
         vxc_aux(:,:) = 0.d0
         CALL v_xc(rho,      rho_core, rhog_core, etxc_aux, vtxc_aux, vxc)
         CALL v_xc(rho_aux,  rho_core, rhog_core, etxc_aux, vtxc_aux, vxc_aux)
         v%of_r(:,:) = v%of_r(:,:) + sic_gamma*(vxc(:,:) - vxc_aux(:,:))
         !
         ! ... SIC energy
         !
         IF(sic_energy) THEN
            esic = 0
            CALL rhoz_or_updw(rho, 'r_and_g', '->updw')
            DO is = 1, nspin
               esic = esic + SUM( sic_gamma*(vxc(:,is) - vxc_aux(:,is))*(rho%of_r(:,is) - rho_n%of_r(:,is)) )
            END DO
            CALL rhoz_or_updw(rho, 'r_and_g', '->rhoz')
            CALL mp_sum(esic, intra_bgrp_comm)
            esic = 0.5*esic*omega/(dfftp%nr1*dfftp%nr2*dfftp%nr3)
         END IF
         !
         ! ... monitor localization through Hartree energy
         ! 
         vh_aux(:,:) = 0.d0
         CALL v_h(rho%pol_g(:,1), eh_aux, charge_aux, vh_aux)
         WRITE(stdout,'(5X,"EH[np] =",f13.8," eV")') eh_aux*ry2ev
         !
         ! ... deallocate
         !
         DEALLOCATE(rho_aux)
         DEALLOCATE(vxc)
         DEALLOCATE(vxc_aux)
         DEALLOCATE(vh_aux)
      END IF
      ! 
   END SUBROUTINE add_vsic
   !
   !--------------------
   SUBROUTINE occ_f2fn()
   !--------------------
      !
      !   ... switch occupation of polaron level from charged to neutral 
      !
      nelec = nelec - fp + fn
      tot_magnetization = 0
      CALL set_nelup_neldw ( tot_magnetization, nelec, nelup, neldw )
      !
   END SUBROUTINE occ_f2fn 
   !
   !--------------------
   SUBROUTINE occ_fn2f()
   !--------------------
      !
      !   ... switch occupation of polaron level from neutral to charged
      !
      nelec = nelec - fn + fp
      tot_magnetization = 1
      CALL set_nelup_neldw ( tot_magnetization, nelec, nelup, neldw )
      !
   END SUBROUTINE occ_fn2f
   !
   !--------------------------
   SUBROUTINE deallocate_sic()
   !--------------------------
      !
      !   ... deallocate SIC module 
      !
      IMPLICIT NONE
      !
      IF(sic_energy .AND. ALLOCATED(rho_n)) DEALLOCATE(rho_n)
      !
   END SUBROUTINE deallocate_sic
   !
   !--------------------------
   SUBROUTINE save_rhon(rho)
   !--------------------------
      !
      !   ... save density rho_n neutral calculation
      !
      IMPLICIT NONE
      !
      TYPE(scf_type), INTENT(INOUT) :: rho ! electron density
      !
      IF(.NOT. ALLOCATED(rho_n)) ALLOCATE(rho_n)
      CALL scf_type_copy(rho,rho_n)
      CALL rhoz_or_updw(rho_n, 'r_and_g', '->updw')
      !
   END SUBROUTINE save_rhon
   !
END MODULE sic_mod
