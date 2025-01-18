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
  !! This is the main driver of the phonon code.
  !! It reads all the quantities calculated by \(\texttt{pwscf}\), it
  !! checks if some recover file is present and determines which 
  !! calculation needs to be done. Finally, it calls \(\texttt{do_phonon}\)
  !! that does the loop over the q points.  
  !! Presently implemented:
  !
  !! * dynamical matrix (\(q\neq 0\))   NC [4], US [4], PAW [4]
  !! * dynamical matrix (\(q=0\))       NC [5], US [5], PAW [4]
  !! * dielectric constant              NC [5], US [5], PAW [3]
  !! * Born effective charges           NC [5], US [5], PAW [3]
  !! * polarizability (iu)              NC [2], US [2]
  !! * electron-phonon                  NC [3], US [3]
  !! * electro-optic                    NC [1]
  !! * Raman tensor                     NC [1]
  !
  !! NC = norm conserving pseudopotentials  
  !! US = ultrasoft pseudopotentials  
  !! PAW = projector augmented-wave  
  !! [1] LDA,  
  !! [2] [1] + GGA,  
  !! [3] [2] + LSDA/sGGA,  
  !! [4] [3] + Spin-orbit/nonmagnetic, non-local vdW functionals, DFT-D2  
  !! [5] [4] + Spin-orbit/magnetic (experimental when available)
  !
  !! Not implemented in \(\texttt{ph.x}\):  
  !! [6] [5] + constraints on the magnetization  
  !! [7] Tkatchenko-Scheffler, DFT-D3  
  !! [8] Hybrid and meta-GGA functionals  
  !! [9] External Electric field  
  !! [10] nonperiodic boundary conditions.
  !
  USE control_flags,   ONLY : use_para_diag
  USE control_ph,      ONLY : bands_computed, qplot
  USE check_stop,      ONLY : check_stop_init
  USE ph_restart,      ONLY : ph_writefile
  USE environment,     ONLY : environment_start
  USE mp_global,       ONLY : mp_startup
  ! YAMBO >
  USE YAMBO,           ONLY : elph_yambo,dvscf_yambo
  ! YAMBO <
  !
  IMPLICIT NONE
  !
  INTEGER :: iq, ierr
  LOGICAL :: do_band, do_iq
  CHARACTER (LEN=9)   :: code = 'PHONON'
  CHARACTER (LEN=256) :: auxdyn
  !
  ! Initialize MPI, clocks, print initial messages
  !
  CALL mp_startup ( start_images=.true., images_only=.true. )
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
  ! ... Do the loop over the q points and irreps.
  !
  CALL do_phonon(auxdyn)
  !
  !  reset the status of the recover files
  !
  CALL ph_writefile('status_ph',1,0,ierr)
  !
  ! YAMBO >
  IF (.not.elph_yambo.and..not.dvscf_yambo) then
    ! YAMBO <
    !
    IF (qplot) CALL write_qplot_data(auxdyn)
    !
    IF (bands_computed) CALL print_clock_pw()
    !
    ! YAMBO >
  ENDIF
  ! YAMBO <
  !
  IF ( use_para_diag ) CALL laxlib_end()
  CALL stop_smoothly_ph( .TRUE. )
  !
  STOP
  !
END PROGRAM phonon
