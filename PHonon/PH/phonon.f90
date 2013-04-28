!
! Copyright (C) 2001-2009 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
PROGRAM phonon
  !-----------------------------------------------------------------------
  !
  ! ... This is the main driver of the phonon code.
  ! ... It reads all the quantities calculated by pwscf, it
  ! ... checks if some recover file is present and determines
  ! ... which calculation needs to be done. Finally, it makes
  ! ... a loop over the q points. At a generic q, if necessary it
  ! ... recalculates the band structure calling pwscf again.
  ! ... Then it can calculate the response to an atomic displacement,
  ! ... the dynamical matrix at that q, and the electron-phonon
  ! ... interaction at that q. At q=0 it can calculate the linear response
  ! ... to an electric field perturbation and hence the dielectric
  ! ... constant, the Born effective charges and the polarizability
  ! ... at imaginary frequencies.
  ! ... At q=0, from the second order response to an electric field,
  ! ... it can calculate also the electro-optic and the raman tensors.
  ! ... Presently implemented:
  ! ... dynamical matrix (q/=0)   NC [4], US [4], PAW [4]
  ! ... dynamical matrix (q=0)    NC [5], US [5], PAW [4]
  ! ... dielectric constant       NC [5], US [5], PAW [3]
  ! ... born effective charges    NC [5], US [5], PAW [3]
  ! ... polarizability (iu)       NC [2], US [2]
  ! ... electron-phonon           NC [3], US [3]
  ! ... electro-optic             NC [1]
  ! ... raman tensor              NC [1]
  !
  ! NC = norm conserving pseudopotentials
  ! US = ultrasoft pseudopotentials
  ! PAW = projector augmented-wave
  ! [1] LDA, 
  ! [2] [1] + GGA, 
  ! [3] [2] + LSDA/sGGA, 
  ! [4] [3] + Spin-orbit/nonmagnetic,
  ! [5] [4] + Spin-orbit/magnetic (experimental when available)
  !
  ! Not implemented in ph.x:
  ! [6] [5] + constraints on the magnetization
  ! [7] [6] + Hubbard U
  ! [8] [7] + Hybrid functionals
  ! [9] ? + External Electric field
  ! [10] ? + nonperiodic boundary conditions.

  USE disp,            ONLY : nqs
  USE control_ph,      ONLY : epsil, trans, bands_computed, qplot, only_init, &
                              only_wfc
  USE el_phon,         ONLY : elph, elph_mat, elph_simple
  USE output,          ONLY : fildrho
  USE check_stop,      ONLY : check_stop_init
  USE ph_restart,      ONLY : ph_writefile
  USE mp_global,       ONLY: mp_startup, nimage
  USE environment,     ONLY: environment_start

  !
  IMPLICIT NONE
  !
  INTEGER :: iq, ierr
  LOGICAL :: do_band, do_iq, setup_pw
  CHARACTER (LEN=9)   :: code = 'PHONON'
  CHARACTER (LEN=256) :: auxdyn
  !
  ! Initialize MPI, clocks, print initial messages
  !
  CALL mp_startup ( start_images=.true. )
  CALL environment_start ( code )
  !
  ! ... and begin with the initialization part
  !
  CALL phq_readin()
  !
  CALL check_stop_init()
  !
  ! ... Checking the status of the calculation and if necessary initialize
  ! ... the q mesh
  !
  CALL check_initial_status(auxdyn)
  !
  DO iq = 1, nqs
     !
     CALL prepare_q(auxdyn, do_band, do_iq, setup_pw, iq)
     !
     !  If this q is not done in this run, cycle
     !
     IF (.NOT.do_iq) CYCLE
     !
     !  If necessary the bands are recalculated
     !
     IF (setup_pw) CALL run_nscf(do_band, iq)
     !
     !  If only_wfc=.TRUE. the code computes only the wavefunctions 
     !
     IF (only_wfc) GOTO 100
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
     !  IF only_init is .true. the code computes only the initialization parts.
     !
     IF (only_init) GOTO 100
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
     call rotate_dvscf_star(iq)
     !
     !  electron-phonon interaction
     !
     IF ( elph ) THEN
        !
        IF ( .NOT. trans ) THEN
           !
           CALL dvanqq()
           IF ( elph_mat ) then
              call ep_matrix_element_wannier()
           ELSE
              CALL elphon()
           END IF
           !
        END IF
        !
        IF ( elph_mat ) then
           call elphsum_wannier(iq)
        ELSEIF( elph_simple ) then
           CALL elphsum_simple()
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
  !
  !  reset the status of the recover files
  !
  CALL ph_writefile('status_ph',1,0,ierr)
  !
  IF (qplot) CALL write_qplot_data(auxdyn)
  !
  IF (bands_computed) CALL print_clock_pw()
  !
  CALL stop_smoothly_ph( .TRUE. )
  !
  STOP
  !
END PROGRAM phonon
