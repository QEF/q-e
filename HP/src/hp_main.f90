!
! Copyright (C) 2001-2018 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
PROGRAM hp_main
  !-----------------------------------------------------------------------
  !
  ! This is the main driver of the HP code.
  !
  USE kinds,             ONLY : DP
  USE io_global,         ONLY : stdout, ionode
  USE check_stop,        ONLY : check_stop_init
  USE mp_global,         ONLY : mp_startup, mp_global_end
  USE environment,       ONLY : environment_start, environment_end
  USE ions_base,         ONLY : nat, ityp, atm, tau, amass
  USE io_files,          ONLY : tmp_dir
  USE ldaU_hp,           ONLY : perturbed_atom, start_q, last_q, nqs, code, &
                                compute_hp, sum_pertq, perturb_only_atom,   &
                                determine_num_pert_only, tmp_dir_save
  !
  IMPLICIT NONE
  !
  INTEGER :: iq, na, ipol
  LOGICAL :: do_iq, setup_pw
  !
  ! Initialize MPI, clocks, print initial messages
  !
#if defined (__MPI)
  CALL mp_startup()
#endif
  !
  CALL environment_start(code)
  !
  ! Print the preamble
  !
  CALL hp_print_preamble()
  !
  ! Read the input parameters and the data produced by PWscf
  !
  CALL hp_readin()
  !
  ! Initialization
  !
  CALL check_stop_init()
  !
  ! Allocate and initialize various arrays needed
  ! for the HP calculation
  !
  CALL hp_init()
  !
  ! Print a summary
  !
  CALL hp_summary()
  !
  IF (compute_hp) GO TO 100
  !
  IF (determine_num_pert_only) GO TO 103
  !
  ! Loop over all atoms in the primitive cell
  !
  DO na = 1, nat
     !
     ! Check if the current atom must be perturbed 
     !
     CALL hp_check_pert(na)
     ! 
     IF (.NOT.perturbed_atom(na) ) CYCLE
     !
     ! The following reinitialization is needed in order to
     ! perform a linear-response calculation for the next atom.
     !
     IF ( na > 1 ) THEN
        CALL clean_pw(.true.)
        CALL close_files(.true.)
        tmp_dir=tmp_dir_save
        CALL read_file()
     ENDIF
     !
     ! Print information about the perturb atom
     !
     WRITE( stdout, * )
     WRITE( stdout, 3336) 
     WRITE( stdout, '(/26x,"PERTURBED ATOM #",2x,i3,/)') na
     WRITE( stdout, '(5x,"site n.  atom      mass           positions (alat units)")')
     WRITE( stdout, '(7x,i2,3x,a6,3x,f8.4,"   tau(",i2,") = (",3f9.5,"  )")')  &
                      & na, atm(ityp(na)), amass(ityp(na)), na, (tau(ipol,na),ipol=1,3)
     WRITE( stdout, * )
     WRITE( stdout, 3336)
     !
     ! Check if the type of the perturbed atom is unique
     ! or there are other atoms of the same type
     !
     CALL hp_check_type(na)
     ! 
     ! Check the status of the calculation and 
     ! initialize the q mesh (q_points) and R mesh (R_points)
     !
     CALL hp_generate_grids()
     !
     IF (sum_pertq) GO TO 102 
     !
     ! Loop over the q points
     !
     DO iq = 1, nqs
        !
        ! Setup a calculation for a specific q point
        !
        CALL hp_prepare_q(iq, do_iq, setup_pw)
        !
        !  If this q is not done in this run, cycle
        !
        IF (.NOT.do_iq) CYCLE
        !
        ! If necessary the bands are recalculated
        !
        IF (setup_pw) CALL hp_run_nscf(.true.) 
        !
        ! Initialize the quantities which do not depend on
        ! the linear response of the system
        !
        CALL hp_load_q()
        !
        ! SCF solution of the linearized Kohn-Sham equation
        !
        CALL hp_solve_linear_system (na, iq)
        !
        ! Write dns0 and dnsscf for the current q to file
        !
        CALL hp_write_dnsq(iq) 
        !
        ! Cleanup of the variables for the next q point
        !
        CALL hp_clean_q(.TRUE.)
        !
     ENDDO
     !
     IF (start_q > 1 .OR. last_q < nqs) THEN
        !
        WRITE( stdout, '(/6x,"Not all q points were considered. Stopping smoothly...",/)')   
        CALL hp_dealloc_1()
        GO TO 103
        !
     ENDIF
     !
102 CONTINUE
     !
     ! Read dns0 and dnsscf for all q from file
     !
     IF (sum_pertq) CALL hp_read_dnsq()
     !
     ! Sum over q of the response occupation matrices 
     !
     CALL hp_dnstot_sum_q()
     !
     ! Calculate the response functions chi0 and chi 
     !
     CALL hp_calc_chi()
     !
     ! Write one column of chi0 and chi to file
     !
     CALL hp_write_chi()
     ! 
     ! Deallocate some arrays
     !
     CALL hp_dealloc_1()
     !
     ! If perturb_only_atom(na)=.true., then this is not a full calculation
     ! but a calculation for only one Hubbard atom na. Hence, stop smoothly.
     !
     IF (perturb_only_atom(na)) GO TO 103
     !
     ! last_q must be recomputed for the next perturbation,
     ! therefore we need to reset it back to -1.
     !
     last_q = -1
     !
  ENDDO
  !
100 CONTINUE
  !
  IF (ionode) THEN
     !
     ! Collect various pieces (columns) of chi0 and chi
     ! (this is needed when various perturbations were
     ! considered not in one single run)
     !
     IF (compute_hp) CALL hp_read_chi()
     !
     ! Write full chi0 and chi to file
     !
     CALL hp_write_chi_full()
     !
  ENDIF
  !
101 CONTINUE
  !
  ! Calculation of Hubbard U (serial post-processing) 
  !
  IF (ionode) CALL hp_postproc()
  !
103 CONTINUE
  !
  ! Deallocate some arrays
  !
  CALL hp_dealloc_2()
  !
  ! Print clocks
  !
  IF (.NOT.compute_hp .AND. .NOT.sum_pertq .AND. .NOT.determine_num_pert_only) THEN
     WRITE( stdout, * )
     WRITE( stdout, * )  '    PRINTING TIMING FROM PWSCF ROUTINES: '
     CALL print_clock_pw()
     CALL hp_print_clock()
  ENDIF
  !
  CALL environment_end(code)
  !
#if defined (__MPI)
  CALL mp_global_end()
#endif
  !
3336 FORMAT('     ',69('='))
  !
  STOP
  !
CONTAINS
  !
SUBROUTINE hp_print_preamble()
  !
  IMPLICIT NONE
  !
  WRITE( stdout, '(/5x,"=--------------------------------------------------------------------------=")')
  WRITE( stdout, '(/7x,"Calculation of Hubbard parameters from DFPT; please cite this program as")')
  WRITE( stdout, '(/7x,"I. Timrov, N. Marzari, and M. Cococcioni, Phys. Rev. B 98, 085127 (2018)")')
  WRITE( stdout, '(/5x,"=--------------------------------------------------------------------------=")')
  !
  RETURN
  !
END SUBROUTINE hp_print_preamble
  !
END PROGRAM hp_main
