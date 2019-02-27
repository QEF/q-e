!
! Copyright (C) 2001-2011 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
SUBROUTINE hinit1()
  !----------------------------------------------------------------------------
  !
  ! ... Atomic configuration dependent hamiltonian initialization
  ! ... Important note: does not recompute structure factors,
  ! ... they must be computed before this routine is called
  !
  USE ions_base,     ONLY : nat, nsp, ityp, tau
  USE cell_base,     ONLY : at, bg, omega, tpiba2
  USE fft_base,      ONLY : dfftp
  USE gvect,         ONLY : ngm, g
  USE gvecs,         ONLY : doublegrid
  USE ldaU,          ONLY : lda_plus_u
  USE lsda_mod,      ONLY : nspin
  USE scf,           ONLY : vrs, vltot, v, kedtau
  USE control_flags, ONLY : tqr, treinit_gvecs 
  USE realus,        ONLY : generate_qpointlist,betapointlist,init_realspace_vars,real_space
  USE wannier_new,   ONLY : use_wannier
  USE martyna_tuckerman, ONLY : tag_wg_corr_as_obsolete
  USE scf,           ONLY : rho
  USE paw_variables, ONLY : okpaw, ddd_paw
  USE paw_onecenter, ONLY : paw_potential
  USE paw_symmetry,  ONLY : paw_symmetrize_ddd
  USE dfunct,        ONLY : newd
  USE fft_base,   ONLY : dfftp
  USE fft_base,   ONLY : dffts
  USE fft_types,  ONLY : fft_type_deallocate
  USE funct,         ONLY : dft_is_hybrid
  USE exx_base,      ONLY : exx_grid_init, exx_mp_init, exx_div_check, coulomb_fac, coulomb_done 
  USE exx,           ONLY : exx_fft_initialized, dfftt, exx_fft_create, deallocate_exx 
  USE exx_band,     ONLY : igk_exx 

  !
  IMPLICIT NONE
  !
  !
  ! ... calculate the total local potential
  !
  CALL setlocal()
  !
  ! these routines can be used to patch quantities that are dependent
  ! on the ions and cell parameters
  !
  CALL plugin_init_ions()
  CALL plugin_init_cell()
  !
  ! ... plugin contribution to local potential
  !
  CALL plugin_scf_potential(rho,.FALSE.,-1.d0,vltot)
  !
  ! ... define the total local potential (external+scf)
  !
  !
  ! ... if tereinit_gvecs is true re-set FFT grids
  !
        IF (  treinit_gvecs  ) THEN  

           CALL clean_pw( .FALSE. ) 
           CALL close_files(.TRUE.)
           dfftp%nr1=0; dfftp%nr2=0; dfftp%nr3=0
           dffts%nr1=0; dffts%nr2=0; dffts%nr3=0
           !
           CALL init_run()
           IF ( dft_is_hybrid() ) THEN
              IF (ALLOCATED(coulomb_fac) .AND. dft_is_hybrid() ) & 
                        DEALLOCATE (coulomb_fac, coulomb_done) 

                CALL deallocate_exx
                if (allocated(igk_exx))  &
                        deallocate(igk_exx) 

                dfftt%nr1=0; dfftt%nr2=0; dfftt%nr3=0 
                CALL fft_type_deallocate(dfftt) 
                CALL exx_grid_init(REINIT = .TRUE.)
                CALL exx_mp_init()
                exx_fft_initialized = .FALSE. 
                CALL exx_div_check()
           ENDIF
        ELSE 
                IF (ALLOCATED(coulomb_fac) .AND. dft_is_hybrid() ) & 
                        DEALLOCATE (coulomb_fac, coulomb_done) 

        END IF  
  CALL set_vrs( vrs, vltot, v%of_r, kedtau, v%kin_r, dfftp%nnr, nspin, doublegrid )
  !
  IF ( tqr ) CALL generate_qpointlist()

  IF (real_space ) then
   call betapointlist()
   call init_realspace_vars()
  endif
  !
  ! ... update the D matrix and the PAW coefficients
  !
  IF (okpaw) THEN
     CALL compute_becsum(1)
     CALL PAW_potential(rho%bec, ddd_paw)
     CALL PAW_symmetrize_ddd(ddd_paw)
  ENDIF
  ! 
  CALL newd()
  !
  ! ... and recalculate the products of the S with the atomic wfcs used 
  ! ... in LDA+U calculations
  !
  IF ( lda_plus_u ) CALL orthoUwfc () 
  IF ( use_wannier ) CALL orthoatwfc( .true. )
  !
  call tag_wg_corr_as_obsolete
  !
  RETURN
  !
END SUBROUTINE hinit1

