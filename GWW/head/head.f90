!
! Copyright (C) 2001-2013 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!

!
!-----------------------------------------------------------------------
PROGRAM head
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
  ! ... dynamical matrix (q/=0)   NC [4], US [4], PAW [3]
  ! ... dynamical matrix (q=0)    NC [5], US [5], PAW [3]
  ! ... dielectric constant       NC [5], US [5], PAW [3] 
  ! ... born effective charges    NC [5], US [5], PAW [3]
  ! ... polarizability (iu)       NC [2], US [2]
  ! ... elctron-phonon            NC [3], US [3]
  ! ... electro-optic             NC [1]
  ! ... raman tensor              NC [1]
  !
  ! NC = norm conserving pseudopotentials
  ! US = ultrasoft pseudopotentials
  ! PAW = projector augmented-wave
  ! [1] LDA, [2] [1]+GGA, [3] [2]+LSDA/sGGA, [4] [3]+Spin-orbit/nonmagnetic,
  ! [5] [4]+Spin-orbit/magnetic
  !
  USE io_global,       ONLY : stdout
  USE disp,            ONLY : nqs
  USE control_ph,      ONLY : epsil, trans,  bands_computed, ldisp
  USE output,          ONLY : fildrho
  USE check_stop,      ONLY : check_stop_init
  USE ph_restart,      ONLY : ph_writefile, destroy_status_run
  USE save_ph,         ONLY : clean_input_variables
  USE mp_global,       ONLY: mp_startup, nimage
  !USE path_io_routines, ONLY : io_path_start
  USE environment,     ONLY: environment_start
  USE wannier_gw,     ONLY : l_head
  USE control_ph,      ONLY : epsil, trans, qplot, only_init, &
                              only_wfc
  USE el_phon,         ONLY : elph, elph_mat, elph_simple
  !
  IMPLICIT NONE
  !
  INTEGER :: iq,ierr
  LOGICAL :: do_band, do_iq, setup_pw
  CHARACTER (LEN=9)   :: code = 'PHONON'
  CHARACTER (LEN=256) :: auxdyn
  !
  ! Initialize MPI, clocks, print initial messages
  !
#if defined(__MPI)
  CALL mp_startup ( )
  !IF (nimage>1) CALL io_path_start()
#endif
  CALL environment_start ( code )
  !
  WRITE( stdout, '(/5x,"Ultrasoft (Vanderbilt) Pseudopotentials")' )
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

  ldisp=.false.
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
     !IF (setup_pw) CALL run_pwscf(do_band)
     IF (setup_pw) CALL run_nscf(do_band, iq)
     !
     !  Initialize the quantities which do not depend on
     !  the linear response of the system
 
     CALL initialize_ph()
     !
     !  electric field perturbation
     !

     IF (epsil) CALL phescf()
     
     if(l_head) then

        call solve_head
     endif
     !
     !  phonon perturbation
     !
     !IF ( trans ) THEN
     !   !
     !   CALL phqscf()
     !   CALL dynmatrix()
     !   !
     !   IF ( fildrho /= ' ' ) CALL punch_plot_ph()
     !   !
     !END IF
     !
     !  electron-phonon interaction
     !
     !IF ( elph ) THEN
     !   !
     !   IF ( .NOT. trans ) THEN
     !      ! 
     !      CALL dvanqq()
     !      CALL elphon()
     !      !
     !   END IF
     !   !
     !   CALL elphsum()
     !   !
     !END IF
     !
     ! ... cleanup of the variables for the next q point
     !
     CALL clean_pw_ph(iq)
     !
  END DO

  CALL ph_writefile('init',0,0,ierr)
  CALL collect_grid_files()
  CALL destroy_status_run()
  !
  IF (bands_computed) CALL print_clock_pw()
  !
  CALL stop_ph( .TRUE. )
  !
  STOP
  !
END PROGRAM head
