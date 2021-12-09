!
! Copyright (C) 2001-2013 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
SUBROUTINE do_phonon(auxdyn)
  !-----------------------------------------------------------------------
  !! This is the main driver of the phonon code. It assumes that the 
  !! preparatory stuff has been already done.  
  !! When the code calls this routine it has already read input
  !! decided which irreducible representations have to be calculated
  !! and it has set the variables that decide which work this routine
  !! will do. The parallel stuff has been already setup by the calling
  !! codes. This routine makes the two loops over
  !! the q-points and the irreps and does only the calculations
  !! that have been decided by the driver routine.  
  !! At a generic q-point, if necessary, it recalculates the band structure 
  !! calling pwscf again. Then it can calculate the response to an atomic
  !! displacement, the dynamical matrix at that q-point, and the 
  !! electron-phonon interaction at that q. At q=0 it can calculate 
  !! the linear response to an electric field perturbation and hence the
  !! dielectric constant, the Born effective charges and the polarizability
  !! at imaginary frequencies.  
  !! At q=0, from the second order response to an electric field,
  !! it can calculate also the electro-optic and the raman tensors.
  !

  USE disp,            ONLY : nqs
  USE control_ph,      ONLY : epsil, trans, qplot, only_init, &
                              only_wfc, rec_code, where_rec, reduce_io
  USE el_phon,         ONLY : elph, elph_mat, elph_simple, elph_epa
  !
  ! YAMBO >
  USE YAMBO,           ONLY : elph_yambo
  ! YAMBO <
  !
  USE elph_tetra_mod, ONLY : elph_tetra, elph_tetra_lambda, elph_tetra_gamma
  USE elph_scdft_mod, ONLY : elph_scdft
  USE io_global,      ONLY : stdout
  ! FIXME: see below setup_pw
  USE noncollin_module, ONLY : noncolin, domag
  USE ahc,            ONLY : elph_ahc, elph_do_ahc
  USE io_files,       ONLY : iunwfc
  USE buffers,        ONLY : close_buffer

  IMPLICIT NONE
  !
  CHARACTER (LEN=256), INTENT(IN) :: auxdyn
  INTEGER :: iq, qind
  LOGICAL :: do_band, do_iq, setup_pw
  !
  qind = 0
  !
  DO iq = 1, nqs
     !
     CALL prepare_q(auxdyn, do_band, do_iq, setup_pw, iq)
     !
     !  If this q is not done in this run, cycle
     !
     IF (.NOT.do_iq) CYCLE
     qind = qind + 1
     !
     !  If necessary the bands are recalculated
     !
     ! Note (A. Urru): This has still to be cleaned (setup_pw 
     ! should be correctly set by prepare_q: here we force it 
     ! to be .true. in order for the code to work properly in 
     ! the case SO-MAG).
     !
     setup_pw=setup_pw .OR. (noncolin .AND. domag)
     IF (setup_pw) THEN
        IF (reduce_io .AND. (qind == 1)) THEN
           CALL close_buffer( iunwfc, 'DELETE' )
        ENDIF
        CALL run_nscf(do_band, iq)
     ENDIF
     !
     !  If only_wfc=.TRUE. the code computes only the wavefunctions 
     !
     IF (only_wfc) THEN
        where_rec='only_wfc'
        rec_code=-1000
        GOTO 100
     ENDIF
     !
     !  Initialize the quantities which do not depend on
     !  the linear response of the system
     !
     CALL initialize_ph()
     !
     !  electric field perturbation
     !
     IF (epsil) CALL phescf()
     !
     !  IF only_init is .true. the code computes only the 
     !  initialization parts.
     !
     IF (only_init) THEN
        where_rec='only_init'
        rec_code=-1000
        GOTO 100
     ENDIF
     !
     !  phonon perturbation
     !
     IF ( trans ) THEN
        !
        CALL phqscf()
        CALL dynmatrix_new(iq)
        !
     END IF
     !
     CALL rotate_dvscf_star(iq)
     !
     !  electron-phonon interaction
     !
     IF ( elph ) THEN
        !
        IF ( .NOT. trans ) THEN
           !
           CALL dvanqq()
           IF ( elph_mat ) THEN
              CALL ep_matrix_element_wannier()
           ELSE
              CALL elphon()
           END IF
           !
        END IF
        !
        IF ( elph_mat ) THEN
           CALL elphsum_wannier(iq)
        ELSEIF( elph_simple ) THEN
           CALL elphsum_simple()
        ELSEIF( elph_epa ) THEN
           CALL elphfil_epa(iq)
        ELSEIF( elph_yambo ) THEN
           CALL elph_yambo_eval_and_IO()
        ELSEIF(elph_tetra == 1) THEN
           CALL elph_tetra_lambda()
        ELSEIF(elph_tetra == 2) THEN
           CALL elph_tetra_gamma()
        ELSEIF(elph_tetra == 3) THEN
           CALL elph_scdft()
        ELSEIF( elph_ahc ) THEN
           CALL elph_do_ahc()
        ELSE 
           CALL elphsum()
        END IF
        !
     END IF
     !
     ! ... cleanup of the variables for the next q point
     !
100  CALL clean_pw_ph(iq)
     !
  END DO

END SUBROUTINE do_phonon
