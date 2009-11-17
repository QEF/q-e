!-----------------------------------------------------------------------
PROGRAM lr_main
  !---------------------------------------------------------------------
  ! Brent Walker, ICTP, 2004
  ! Dario Rocca, SISSA, 2006
  ! O. Baris Malcioglu, SISSA, 2008
  !---------------------------------------------------------------------
  ! ... overall driver routine for applying lanczos algorithm
  ! ... to the matrix of equations coming from tddft
  ! ... spectrum version
  !---------------------------------------------------------------------
  !
#include "f_defs.h"
  !
  USE lr_lanczos,            ONLY : one_lanczos_step
  USE io_global,             ONLY : stdout
  USE kinds,                 ONLY : dp
  USE lr_variables,          ONLY : restart, restart_step,&
                                    itermax, lr_verbosity,&
                                    evc1, norm0, charge_response,&
                                    n_ipol, d0psi, rho_1_tot, &
                                    LR_iteration, LR_polarization, &
                                    plot_type
  USE io_files,              ONLY : nd_nmbr
  USE global_version,        ONLY : version_number
  USE charg_resp,            ONLY : lr_calc_w_T, read_wT_beta_gamma_z, &
                                    lr_dump_rho_tot_compat1, lr_dump_rho_tot_cube,&
                                    lr_dump_rho_tot_xyzd,lr_dump_rho_tot_xcrys,&
                                    lr_dump_rho_tot_pxyd
  USE ions_base,             ONLY : tau,nat,atm,ityp
  USE environment,           ONLY: environment_start
  USE mp_global,             ONLY : nimage, mp_startup
  USE control_flags,         ONLY : use_task_groups, ortho_para
  !
  IMPLICIT NONE
  !
  ! Local variables
  !
  !INTEGER           :: iter, ip, op
  INTEGER           :: ip,pol_index
  INTEGER           :: iter_restart
  LOGICAL           :: rflag
  CHARACTER (len=9) :: code = 'TDDFPT'
  !
  pol_index=1
  !CALL init_clocks( .TRUE. )
  !
  !
  !CALL startup (nd_nmbr, code, version_number)
#ifdef __PARA
  CALL mp_startup ( use_task_groups, ortho_para )
#endif
  CALL environment_start ( 'TDDFPT' )
  !
  CALL start_clock('lr_main')
  !
  WRITE( stdout, '(/5x,"Ultrasoft (Vanderbilt) Pseudopotentials")' )
  If (lr_verbosity > 5) THEN
   WRITE(stdout,'("<lr_main>")')
  endif
  !
  !   Reading input file
  !
  CALL lr_readin ( )
  !
  WRITE(stdout,'(/,5X,"Lanczos linear response spectrum calculation")')
  WRITE(stdout,'(5x,"Number of Lanczos iterations = ",i6)') itermax
  !
  !   Allocate and zero lr variables
  !
  !OBM_DEBUG
  If (lr_verbosity > 6) THEN
     WRITE(stdout,'(/,5X,"Step-main1")')
  endif
  !OBM_DEBUG
  !
  if (restart)  call test_restart()
  !
  CALL lr_alloc_init()  
  !
  !  Charge response: initialisation
  !
  !
  !   Read in ground state wavefunctions
  ! 
  !OBM_DEBUG
  If (lr_verbosity > 6) THEN
   WRITE(stdout,'(/,5X,"Step-main2")')
  endif
  !OBM_DEBUG
  CALL lr_read_wf()
  !
  !   Set up initial response orbitals
  !
  !OBM_DEBUG
  If (lr_verbosity > 6) THEN
   WRITE(stdout,'(/,5X,"Step-main3")')
  endif
  !OBM_DEBUG
  !
  !
  !Initialisation of degauss/openshell related stuff
  !
  call lr_init_nfo()
  !
  IF ( restart )      CALL lr_read_d0psi()
  IF ( .NOT.restart ) CALL lr_solve_e()
  !
  !   Set up initial stuff for derivatives
  !
  !OBM_DEBUG
  If (lr_verbosity > 6) THEN
   WRITE(stdout,'(/,5X,"Step-main4")')
  endif
  !OBM_DEBUG
  CALL lr_dv_setup()
  !OBM_DEBUG
  If (lr_verbosity > 6) THEN
   WRITE(stdout,'(/,5X,"Step-main5")')
  endif
  !OBM_DEBUG
  !Coordinates of the read atom, just in case
  IF (lr_verbosity > 0) THEN
   WRITE(stdout,'(/,5X,"Positions of atoms in internal coordinates")')
   DO ip=1,nat ! I am using ip here as a counter over atoms
      WRITE(stdout,'(A3,2X,3F15.8)') atm(ityp(ip)),tau(1:3,ip)
   ENDDO
  ENDIF
  !
  !   Lanczos loop where the real work happens
  !
   DO ip=1, n_ipol
    if (n_ipol/=1) then 
      LR_polarization=ip
      pol_index=LR_polarization
    endif
    !WRITE(stdout,'(5x,"Lanczos iterations for polarization direction ",i2)') LR_polarization
    ! 
    if (charge_response == 2 ) then  
         !
         ! Read precalculated beta gamma z
         !
         rho_1_tot(:)=0 !zero the response charge at the beginning of a loop
         call read_wT_beta_gamma_z() 
         call lr_calc_w_T()
     endif 
     !
     !
     IF (restart) then 
      !
        !
        CALL lr_restart(iter_restart,rflag)
        CALL lr_read_d0psi()
        !
       !
     else
      !  
        !
        CALL lr_read_d0psi()
        evc1(:,:,:,1) = d0psi(:,:,:,pol_index)
        CALL lr_normalise( evc1(:,:,:,1), norm0(pol_index) )
        evc1(:,:,:,2) = evc1(:,:,:,1)
        !
        iter_restart=1
        !
      ! 
     END IF
     !
     CALL sd0psi() 
     !
     IF ( .NOT.restart ) write(stdout,'(/5x,"Starting Lanczos loop",1x,i8)')   LR_polarization
     IF ( restart )      write(stdout,'(/5x,"Restarting Lanczos loop",1x,i8)') LR_polarization
     !
     lancz_loop1 : DO LR_iteration = iter_restart, itermax
        !
        write(stdout,'(/5x,"Lanczos iteration:",1x,i8)') LR_iteration
        !
        call one_lanczos_step()
        !
        IF ( mod(LR_iteration,restart_step)==0 .OR. LR_iteration==itermax ) CALL lr_write_restart()
        !
     END DO lancz_loop1
     ! 
    if (charge_response == 1 .and. lr_verbosity > 0) then
         call lr_calc_w_T()
    endif
    if (charge_response == 2 ) then 
      !call lr_dump_rho_tot_compat1()
      if (plot_type == 1 .or. plot_type == 5) call lr_dump_rho_tot_xyzd(rho_1_tot,"summed-rho")
      if (plot_type == 2 .or. plot_type == 5) call lr_dump_rho_tot_xcrys(rho_1_tot,"summed-rho")
      if (plot_type == 3 .or. plot_type == 5) call lr_dump_rho_tot_cube(rho_1_tot,"summed-rho")
      !call lr_dump_rho_tot_pxyd(rho_1_tot,"summed-rho")
      !call lr_dump_rho_tot_xcrys(rho_1_tot,"summed-rho")
    endif 
     !
  END DO
  !
  !  if(lr_verbosity>2) call lr_diagonalise(itermax-1)
  !
  WRITE(stdout,'(5x,"End of Lanczos iterations")') 
  !
  !
  !   Calculate the spectrum !This is now done by a post processing program
  !
  !If (lr_verbosity > 3) THEN
  !  CALL lr_calculate_spectrum()
  !endif
  !
  !   Deallocate pw variables
  !
  CALL clean_pw( .FALSE. )  
  !
  WRITE(stdout,'(5x,"Finished linear response calculation...")')
  !
  CALL stop_clock('lr_main')
  !
  CALL print_clock_lr()
  !
  CALL stop_lr()
  !
  If (lr_verbosity > 5) THEN
   WRITE(stdout,'("<end of lr_main>")')
  endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Additional small-time subroutines
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
CONTAINS
 SUBROUTINE test_restart()
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!This tests whether the restart flag was used correctly
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 use lr_variables,     only : n_ipol,LR_polarization,restart
 use io_files,         only: prefix, tmp_dir, nd_nmbr
 USE mp,               ONLY : mp_bcast, mp_barrier
 USE io_global,        ONLY : ionode, ionode_id

 IMPLICIT NONE
  character(len=256) :: tempfile, filename
  logical :: exst
  character(len=6), external :: int_to_char
 
 if ( n_ipol == 1 ) then
  filename = trim(prefix)//'.d0psi.'//trim(int_to_char(1))
  tempfile = trim(tmp_dir) // trim(filename) //nd_nmbr
 else 
  filename = trim(prefix)//'.d0psi.'//trim(int_to_char(LR_polarization))
  tempfile = trim(tmp_dir) // trim(filename) //nd_nmbr
 endif
 inquire (file = tempfile, exist = exst)
 !print *, tempfile," exst=",exst
 if (.not. exst) then
    WRITE(stdout,'("Warning: Some files are missing, unable to restart.")')
    restart=.false.
 endif
 !print *,"a"
 if ( n_ipol == 1 ) then
  filename = trim(prefix)//'.restart_lanczos.'//trim(int_to_char(1))
  tempfile = trim(tmp_dir) // trim(filename) //nd_nmbr
 else 
  filename = trim(prefix)//'.restart_lanczos.'//trim(int_to_char(LR_polarization))
  tempfile = trim(tmp_dir) // trim(filename)//nd_nmbr
 endif
 inquire (file = tempfile, exist = exst)
 !print *, tempfile," exst=",exst
 if (.not. exst) then
    WRITE(stdout,'("Warning: Some files are missing, unable to restart.")')
    restart=.false.
 endif


 if ( n_ipol == 1 ) then
  filename = trim(prefix) // ".beta_gamma_z." // trim(int_to_char(1))
  tempfile = trim(tmp_dir) // trim(filename)
 else 
  filename = trim(prefix) // ".beta_gamma_z." // trim(int_to_char(LR_polarization))
  tempfile = trim(tmp_dir) // trim(filename)
 endif

 inquire (file = tempfile, exist = exst)
 !print *, tempfile," exst=",exst
 if (.not. exst) then
    WRITE(stdout,'("Warning: Some files are missing, unable to restart.")')
    restart=.false.
 endif

#ifdef __PARA
 call mp_bcast(restart, ionode_id)
#endif

 END SUBROUTINE test_restart
END PROGRAM lr_main
!-----------------------------------------------------------------------
