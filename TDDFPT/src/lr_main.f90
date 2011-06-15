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
                                    n_ipol, d0psi, rho_1_tot, rho_1_tot_im,&
                                    LR_iteration, LR_polarization, &
                                    plot_type, no_hxc, nbnd_total, project, F,R, &
                                    itermax_int, revc0, lr_io_level
  USE io_files,              ONLY : nd_nmbr
  USE global_version,        ONLY : version_number
  USE charg_resp,            ONLY : lr_calc_w_T, read_wT_beta_gamma_z, &
                                    lr_dump_rho_tot_compat1, lr_dump_rho_tot_cube,&
                                    lr_dump_rho_tot_xyzd,lr_dump_rho_tot_xcrys,&
                                    lr_dump_rho_tot_pxyd,chi,lr_calc_R,w_t_norm0_store,&
                                    resonance_condition
  USE ions_base,             ONLY : tau,nat,atm,ityp
  USE environment,           ONLY: environment_start
  USE mp_global,             ONLY : nimage, mp_startup

  USE control_ph,            ONLY : nbnd_occ
  USE wvfct,                 ONLY : nbnd
  USE wavefunctions_module,  ONLY : psic
  USE control_flags,          ONLY : tddfpt
  !Debugging
  USE lr_variables, ONLY: check_all_bands_gamma, check_density_gamma,check_vector_gamma
  !
  IMPLICIT NONE
  !
  ! Local variables
  !
  !INTEGER           :: iter, ip, op
  INTEGER           :: ip,pol_index,ibnd_occ,ibnd_virt,ibnd
  INTEGER           :: iter_restart,iteration
  LOGICAL           :: rflag
  CHARACTER (len=9) :: code = 'TDDFPT'
  COMPLEX(kind=dp)     :: sum_F,sum_c
  !
  !
  pol_index=1
  !CALL init_clocks( .TRUE. )
  !
  !
  !CALL startup (nd_nmbr, code, version_number)
#ifdef __PARA
  CALL mp_startup ( )
#endif
  CALL environment_start ( 'TDDFPT' )
  WRITE( stdout, '(/5x,"----------------------------------------")' )
  WRITE( stdout, '(/5x,"This is TDDFPT (Time Dependent Density Functional Perturbation Theory)")' )
  WRITE( stdout, '(/5x,"Sub Version: 0.9 ")' )
!  WRITE( stdout, '(/5x,"")' )
!  WRITE( stdout, '(/5x,"Please cite this project as:  ")' )
!  WRITE( stdout, '(/5x,"O.B.Malcioglu, R. Gebauer, D. Rocca, S. Baroni,")' )
!  WRITE( stdout, '(/5x,"""turboTDDFT – a code for the simulation of molecular")' )
!  WRITE( stdout, '(/5x,"spectra using the Liouville-Lanczos approach to")' )
!  WRITE( stdout, '(/5x,"time-dependent density-functional perturbation theory""")' )
!  WRITE( stdout, '(/5x,"CPC, (in press)")' )

  WRITE( stdout, '(/5x,"----------------------------------------")' )
  !
  CALL start_clock('lr_main')
  !
  WRITE( stdout, '(/5x,"Ultrasoft (Vanderbilt) Pseudopotentials")' )
  IF (lr_verbosity > 5) THEN
   WRITE(stdout,'("<lr_main>")')
  ENDIF
  !Let the phonon routines know that they are doing tddfpt.
  tddfpt=.TRUE.
  !
  !   Reading input file and PWSCF xml, some initialisation
  !
  CALL lr_readin ( )
  !
  WRITE(stdout,'(/,5X,"Lanczos linear response spectrum calculation")')
  WRITE(stdout,'(5x,"Number of Lanczos iterations = ",i6)') itermax
  !
  !
  IF (no_hxc)  WRITE(stdout,'(5x,"No Hartree/Exchange/Correlation")')
  !
  !   Allocate and zero lr variables
  !
  !
  !Initialisation of degauss/openshell related stuff
  !
  CALL lr_init_nfo()
  IF (project) THEN
   IF(nbnd > nbnd_occ(1)) THEN
    WRITE(stdout,'(/,5X,"Virtual states in ground state run will be used in projection analysis")')
   ELSE
    WRITE(stdout,'(/,5X,"No virtual states for projection found")')
    project=.false.
   ENDIF
  ENDIF

  IF (nbnd>nbnd_occ(1)) THEN
   WRITE(stdout,'(/,5X,"Warning: There are virtual states in the input file, trying to disregard in response calculation")')
   nbnd_total=nbnd
   nbnd=nbnd_occ(1)
  ELSE
   nbnd_total=nbnd
  ENDIF
  !
  CALL lr_alloc_init()
  !
  !  Charge response: initialisation
  !
  !
  !   Read in ground state wavefunctions
  !
  CALL lr_read_wf()
  !
  !   Set up initial response orbitals
  !
  !
  IF ( test_restart(1) ) THEN
    CALL lr_read_d0psi()
  ELSE
    CALL lr_solve_e()
  ENDIF

  DEALLOCATE( psic )

  IF(project) THEN
     CALL sd0psi() !after this d0psi is Sd0psi the d0psi is read afterwards again...
     CALL lr_calc_R()
     DO ip=1, n_ipol
      WRITE(stdout,'(/,/5x,"Oscillator strengths for polarization direction",1x,i8)') ip
      WRITE(stdout,'(5x,"occ",1x,"con",8x,"Re(R)",14x,"Im(R)")')
      DO ibnd_occ=1,nbnd
       DO ibnd_virt=1,(nbnd_total-nbnd)
       WRITE(stdout,'(5x,i3,1x,i3,3x,E16.8,2X,E16.8)') &
       ibnd_occ,ibnd_virt,dble(R(ibnd_occ,ibnd_virt,ip)),aimag(R(ibnd_occ,ibnd_virt,ip))
       ENDDO
      ENDDO
     ENDDO
  ENDIF
  !
  !   Set up initial stuff for derivatives
  !
  CALL lr_dv_setup()

  !Coordinates of the read atom, just in case
  IF (lr_verbosity > 1) THEN
   WRITE(stdout,'(/,5X,"Positions of atoms in internal coordinates")')
   DO ip=1,nat ! I am using ip here as a counter over atoms
      WRITE(stdout,'(5X,A3,2X,3F15.8)') atm(ityp(ip)),tau(1:3,ip)
   ENDDO
  ENDIF
  !
  !   Lanczos loop where the real work happens
  !
   DO ip=1, n_ipol
    IF (n_ipol/=1) THEN
      LR_polarization=ip
      pol_index=LR_polarization
    ENDIF
    !
    IF (charge_response == 1 ) THEN
         !
         ! Read precalculated beta gamma z
         !
         CALL read_wT_beta_gamma_z()
         CALL lr_calc_w_T()
     ENDIF
     !
     !
     IF (test_restart(2)) THEN
      !
        !
        CALL lr_restart(iter_restart,rflag)
        CALL lr_read_d0psi()
        !
        WRITE(stdout,'(/5x,"Restarting Lanczos loop",1x,i8)') LR_polarization
       !
     ELSE
      !
        !
        CALL lr_read_d0psi()
        evc1(:,:,:,1) = d0psi(:,:,:,pol_index)
        CALL lr_normalise( evc1(:,:,:,1), norm0(pol_index) )
        evc1(:,:,:,2) = evc1(:,:,:,1)
        !
        iter_restart=1
        !
        WRITE(stdout,'(/5x,"Starting Lanczos loop",1x,i8)')   LR_polarization
      !
     ENDIF

     !
     CALL sd0psi() !after this d0psi is Sd0psi !OBM:Check if this is really necessary
     !
     !

     !
     lancz_loop1 : DO iteration = iter_restart, itermax
        !
        LR_iteration=iteration
        WRITE(stdout,'(/5x,"Lanczos iteration:",1x,i8)') LR_iteration
        !
        CALL one_lanczos_step()
        !
        IF ( lr_io_level > 0 .and. (mod(LR_iteration,restart_step)==0 .or. &
                              LR_iteration==itermax .or. LR_iteration==1) )&
           CALL lr_write_restart()
     ENDDO lancz_loop1
     !

    IF (charge_response == 1 ) THEN
      IF (resonance_condition) THEN
       !response charge density, absorbtive
       IF (plot_type == 1 .or. plot_type == 5) &
        CALL lr_dump_rho_tot_xyzd(aimag(rho_1_tot_im(:,1)),"absorbtive")
       IF (plot_type == 2 .or. plot_type == 5) &
        CALL lr_dump_rho_tot_xcrys(aimag(rho_1_tot_im(:,1)),"absorbtive")
       IF (plot_type == 3 .or. plot_type == 5) &
        CALL lr_dump_rho_tot_cube(aimag(rho_1_tot_im(:,1)),"absorbtive")
       !response charge density, dispersive
       IF (plot_type == 1 .or. plot_type == 5) &
        CALL lr_dump_rho_tot_xyzd(dble(rho_1_tot_im(:,1)),"dispersive")
       IF (plot_type == 2 .or. plot_type == 5) &
        CALL lr_dump_rho_tot_xcrys(dble(rho_1_tot_im(:,1)),"dispersive")
       IF (plot_type == 3 .or. plot_type == 5) &
        CALL lr_dump_rho_tot_cube(dble(rho_1_tot_im(:,1)),"dispersive")
      ELSE
      IF (plot_type == 1 .or. plot_type == 5) CALL lr_dump_rho_tot_xyzd(rho_1_tot(:,1),"summed-rho")
      IF (plot_type == 2 .or. plot_type == 5) CALL lr_dump_rho_tot_xcrys(rho_1_tot(:,1),"summed-rho")
      IF (plot_type == 3 .or. plot_type == 5) CALL lr_dump_rho_tot_cube(rho_1_tot(:,1),"summed-rho")
      ENDIF
    ENDIF
     IF (project) THEN
      WRITE(stdout,'(/,/5x,"Projection of virtual states for polarization direction",1x,i8)') LR_polarization
      WRITE(stdout,'(2x,"occ",1x,"vir",8x,"Re(F)",14x,"Im(F)",8x, &
    & " Frac. pres. in Re(chi_",I1,"_",I1,") and Im(chi_",I1,"_",I1,")")') &
    &  ip,ip,ip,ip
       sum_f=cmplx(0.0d0,0.0d0,dp)
      DO ibnd_occ=1,nbnd
       DO ibnd_virt=1,(nbnd_total-nbnd)
       F(ibnd_occ,ibnd_virt,ip)=F(ibnd_occ,ibnd_virt,ip)*cmplx(w_T_norm0_store,0.0d0,dp)
       sum_f=F(ibnd_occ,ibnd_virt,ip)*conjg(R(ibnd_occ,ibnd_virt,ip))
       WRITE(stdout,'(2x,i3,1x,i3,3x,E16.8,2X,E16.8,17X,F8.5,2x,F8.5)') &
       ibnd_occ,ibnd_virt,dble(F(ibnd_occ,ibnd_virt,ip)),aimag(F(ibnd_occ,ibnd_virt,ip)),&
       (dble(sum_f)/dble(chi(ip,ip))), (aimag(sum_f)/aimag(chi(ip,ip)))
       ENDDO
      ENDDO
     ENDIF
     !
  ENDDO
  DEALLOCATE( revc0 )
  !
  !
  WRITE(stdout,'(5x,"End of Lanczos iterations")')
  !
  !Final reports
  !
  IF (project .and. n_ipol == 3) THEN
      WRITE(stdout,'(/,/5x,"Participation of virtual states to absorbtion coefficent")')
      WRITE(stdout,'(5x,"occ",1x,"vir",5x,"Re(Tr(F.R))",6x,"Im(TR(F.R))",5x,"fraction in alpha")')
      DO ibnd_occ=1,nbnd
       DO ibnd_virt=1,(nbnd_total-nbnd)
       sum_f=cmplx(0.0d0,0.0d0,dp)
       sum_c=cmplx(0.0d0,0.0d0,dp)
        DO ip=1,n_ipol
         sum_f=sum_f+F(ibnd_occ,ibnd_virt,ip)*conjg(R(ibnd_occ,ibnd_virt,ip))
         sum_c=sum_c+chi(ip,ip)
        ENDDO
        WRITE(stdout,'(5x,i3,1x,i3,3x,E16.8,2X,E16.8,2X,F8.5)') &
       ibnd_occ,ibnd_virt,dble(sum_F),aimag(sum_F),(aimag(sum_f)/aimag(sum_c))
       ENDDO
      ENDDO
  ENDIF

  !
  !   Deallocate pw variables
  !
  CALL clean_pw( .false. )
  !
  WRITE(stdout,'(5x,"Finished linear response calculation...")')
  !
  CALL stop_clock('lr_main')
  !
  CALL print_clock_lr()
  !
  CALL stop_lr()
  !
  IF (lr_verbosity > 5) THEN
   WRITE(stdout,'("<end of lr_main>")')
  ENDIF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Additional small-time subroutines
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
CONTAINS
 LOGICAL FUNCTION test_restart(test_this)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!This tests whether the restart flag is applicable
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 USE lr_variables,     ONLY : n_ipol,LR_polarization,restart,bgz_suffix
 USE io_files,         ONLY: prefix, tmp_dir, nd_nmbr, wfc_dir
 USE mp,               ONLY : mp_bcast, mp_barrier,mp_sum
 USE io_global,        ONLY : ionode, ionode_id

 IMPLICIT NONE
  INTEGER, INTENT(in) :: test_this
  CHARACTER(len=256) :: tempfile, filename, tmp_dir_saved
  LOGICAL :: exst
  CHARACTER(len=6), EXTERNAL :: int_to_char
  INTEGER :: i, temp_restart

 !
 !test_this= 1 : d0psi files
 !test_this= 2 : lanczos restart files

 temp_restart=0
 !print *, "test_restart with restart=",restart
 IF (.not.restart) THEN
  test_restart = .false.
  RETURN
 ENDIF
 test_restart=.true.
 IF (test_this == 1) THEN
 !
 !Check for parallel i/o files that are in wfc_dir
 tmp_dir_saved = tmp_dir
 IF ( wfc_dir /= 'undefined' ) tmp_dir = wfc_dir
 !
 IF ( n_ipol == 1 ) THEN
  filename = trim(prefix)//'.d0psi.'//trim(int_to_char(1))
  tempfile = trim(tmp_dir) // trim(filename) //nd_nmbr
  INQUIRE (file = tempfile, exist = exst)
  !print *, tempfile," exst=",exst
  IF (.not. exst) THEN
    temp_restart=1
  ENDIF
 ELSE
  DO i=1, n_ipol
   filename = trim(prefix)//'.d0psi.'//trim(int_to_char(i))
   tempfile = trim(tmp_dir) // trim(filename) //nd_nmbr
   INQUIRE (file = tempfile, exist = exst)
   !print *, tempfile," exst=",exst
   IF (.not. exst) THEN
     temp_restart=1
   ENDIF
  ENDDO
 ENDIF

 tmp_dir = tmp_dir_saved

 IF ( wfc_dir /= 'undefined' ) THEN
 ! check if these files can be read from outdir instead of wfcdir
 !
 IF ( n_ipol == 1 ) THEN
  filename = trim(prefix)//'.d0psi.'//trim(int_to_char(1))
  tempfile = trim(tmp_dir) // trim(filename) //nd_nmbr
  INQUIRE (file = tempfile, exist = exst)
  IF (exst) THEN
    temp_restart=0
  ENDIF
 ELSE
  DO i=1, n_ipol
   filename = trim(prefix)//'.d0psi.'//trim(int_to_char(i))
   tempfile = trim(tmp_dir) // trim(filename) //nd_nmbr
   INQUIRE (file = tempfile, exist = exst)
   IF (exst) THEN
     temp_restart=0
   ENDIF
  ENDDO
 ENDIF
 ENDIF
 ENDIF !for test_this = 1
 IF (test_this == 2) THEN

 !Restart files are always written in outdir
 IF ( n_ipol == 1 ) THEN
  filename = trim(prefix)//'.restart_lanczos.'//trim(int_to_char(1))
  tempfile = trim(tmp_dir) // trim(filename) //nd_nmbr
 ELSE
  filename = trim(prefix)//'.restart_lanczos.'//trim(int_to_char(LR_polarization))
  tempfile = trim(tmp_dir) // trim(filename)//nd_nmbr
 ENDIF
 INQUIRE (file = tempfile, exist = exst)
 !print *, tempfile," exst=",exst
 IF (.not. exst) THEN
    temp_restart=1
 ENDIF
 !
 !End of parallel file i/o
 !
 IF ( n_ipol == 1 ) THEN
  filename = trim(prefix) // trim(bgz_suffix) // trim(int_to_char(1))
  tempfile = trim(tmp_dir) // trim(filename)
 ELSE
  filename = trim(prefix) // trim(bgz_suffix) // trim(int_to_char(LR_polarization))
  tempfile = trim(tmp_dir) // trim(filename)
 ENDIF
 INQUIRE (file = tempfile, exist = exst)
 !print *, tempfile," exst=",exst
 IF (.not. exst) THEN
    temp_restart=1
 ENDIF
 ENDIF !for test_this = 2

  !print *,"temp_restart",temp_restart
#ifdef __PARA
    CALL mp_sum(temp_restart)
#endif
  !print *, "current temp_restart", temp_restart
  IF (temp_restart > 0 ) THEN
   !print *,"restart falsified",nd_nmbr
   !WRITE(stdout,'(5X,A,3X,"is missing, unable to restart.")') offender
   WRITE(stdout,'(5X,"There are missing files!")')
   IF (test_this==1) WRITE(stdout,'(5X,"d0psi files can not be found,trying to recompansate")')
   IF (test_this==2) WRITE(stdout,'(5X,"lanczos restart files can not be found, starting run from scratch")')
   test_restart=.false.
  ENDIF

  RETURN
 END FUNCTION test_restart
END PROGRAM lr_main
!-----------------------------------------------------------------------
