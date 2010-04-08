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
                                    plot_type, no_hxc, nbnd_total, project, F,R, &
                                    itermax_int
  USE io_files,              ONLY : nd_nmbr
  USE global_version,        ONLY : version_number
  USE charg_resp,            ONLY : lr_calc_w_T, read_wT_beta_gamma_z, &
                                    lr_dump_rho_tot_compat1, lr_dump_rho_tot_cube,&
                                    lr_dump_rho_tot_xyzd,lr_dump_rho_tot_xcrys,&
                                    lr_dump_rho_tot_pxyd,chi,lr_calc_R
  USE ions_base,             ONLY : tau,nat,atm,ityp
  USE environment,           ONLY: environment_start
  USE mp_global,             ONLY : nimage, mp_startup
  
  USE control_ph,            ONLY : nbnd_occ
  use wvfct,                 only : nbnd
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
  complex(kind=dp)     :: sum_F,sum_c
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
  WRITE( stdout, '(/5x,"Sub Version: BETA 1 ")' )
  WRITE( stdout, '(/5x,"")' )
  WRITE( stdout, '(/5x,"Please cite this project as:  ")' )
  WRITE( stdout, '(/5x,"----------------------------------------")' )
  !
  CALL start_clock('lr_main')
  !
  WRITE( stdout, '(/5x,"Ultrasoft (Vanderbilt) Pseudopotentials")' )
  If (lr_verbosity > 5) THEN
   WRITE(stdout,'("<lr_main>")')
  endif
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
  !OBM_DEBUG
  If (lr_verbosity > 6) THEN
     WRITE(stdout,'(/,5X,"Step-main1")')
  endif
  !OBM_DEBUG
  !
  !
  !Initialisation of degauss/openshell related stuff
  !
  call lr_init_nfo()
  if (project) then
   if(nbnd > nbnd_occ(1)) then 
    WRITE(stdout,'(/,5X,"Virtual states in ground state run will be used in projection analysis")')
   else
    WRITE(stdout,'(/,5X,"No virtual states for projection found")')
    project=.false.
   endif
  endif

  IF (nbnd>nbnd_occ(1)) then
   WRITE(stdout,'(/,5X,"Warning: There are virtual states in the input file, trying to disregard in response calculation")')
   nbnd_total=nbnd
   nbnd=nbnd_occ(1)
  else
   nbnd_total=nbnd
  endif
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
  IF ( test_restart() ) then 
    CALL lr_read_d0psi()
  else
    CALL lr_solve_e()
  endif
  if(project) then
     CALL sd0psi() !after this d0psi is Sd0psi the d0psi is read afterwards again...
     call lr_calc_R()
     DO ip=1, n_ipol
      write(stdout,'(/,/5x,"Oscillator strengths for polarization direction",1x,i8)') ip
      write(stdout,'(5x,"occ",1x,"con",8x,"Re(R)",14x,"Im(R)")')
      do ibnd_occ=1,nbnd
       do ibnd_virt=1,(nbnd_total-nbnd)
       write(stdout,'(5x,i3,1x,i3,3x,E16.8,2X,E16.8)') & 
       ibnd_occ,ibnd_virt,DBLE(R(ibnd_occ,ibnd_virt,ip)),AIMAG(R(ibnd_occ,ibnd_virt,ip)) 
       enddo
      enddo
     enddo
  endif
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
         rho_1_tot(:,:)=0.0D0 !zero the response charge at the beginning of a loop
         call read_wT_beta_gamma_z() 
         call lr_calc_w_T()
     endif 
     !
     !
     IF (test_restart()) then 
      !
        !
        CALL lr_restart(iter_restart,rflag)
        CALL lr_read_d0psi()
        !
        write(stdout,'(/5x,"Restarting Lanczos loop",1x,i8)') LR_polarization
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
        write(stdout,'(/5x,"Starting Lanczos loop",1x,i8)')   LR_polarization
      ! 
     END IF 
     if (lr_verbosity >10) then
      write(stdout,'("d0psi")')
      do ibnd=1,nbnd
             call check_vector_gamma(d0psi(:,ibnd,1,pol_index))
      enddo
     endif

     ! 
     CALL sd0psi() !after this d0psi is Sd0psi !OBM:Check if this is really necessary
     !
     !
     if (lr_verbosity >10) then
      write(stdout,'("initial evc1")')
      do ibnd=1,nbnd
             call check_vector_gamma(evc1(:,ibnd,1,1))
      enddo
      write(stdout,'("initial sd0psi")')
      do ibnd=1,nbnd
             call check_vector_gamma(d0psi(:,ibnd,1,pol_index))
      enddo
     endif

     !
     lancz_loop1 : DO iteration = iter_restart, itermax
        !
        LR_iteration=iteration
        write(stdout,'(/5x,"Lanczos iteration:",1x,i8)') LR_iteration
        !
        call one_lanczos_step()
        !
        IF ( mod(LR_iteration,restart_step)==0 .OR. LR_iteration==itermax .OR. LR_iteration==1 ) CALL lr_write_restart()
     END DO lancz_loop1
     ! 
    if (charge_response == 1 .and. lr_verbosity > 0) then
         call lr_calc_w_T()
    endif
    if (charge_response == 2 ) then 
      !call lr_dump_rho_tot_compat1()
      if (plot_type == 1 .or. plot_type == 5) call lr_dump_rho_tot_xyzd(rho_1_tot(:,1),"summed-rho")
      if (plot_type == 2 .or. plot_type == 5) call lr_dump_rho_tot_xcrys(rho_1_tot(:,1),"summed-rho")
      if (plot_type == 3 .or. plot_type == 5) call lr_dump_rho_tot_cube(rho_1_tot(:,1),"summed-rho")
      !call lr_dump_rho_tot_pxyd(rho_1_tot,"summed-rho")
      !call lr_dump_rho_tot_xcrys(rho_1_tot,"summed-rho")
    endif
     if (project) then
      write(stdout,'(/,/5x,"Projection of virtual states for polarization direction",1x,i8)') LR_polarization
      write(stdout,'(5x,"occ",1x,"con",8x,"Re(F)",14x,"Im(F)",5x,"%chi_",I1,"_",I1)') ip,ip
      do ibnd_occ=1,nbnd
       do ibnd_virt=1,(nbnd_total-nbnd)
       sum_f=F(ibnd_occ,ibnd_virt,ip)*R(ibnd_occ,ibnd_virt,ip)
       sum_f=sum_f/chi(ip,ip)
       write(stdout,'(5x,i3,1x,i3,3x,E16.8,2X,E16.8,2X,"%",F5.2)') & 
       ibnd_occ,ibnd_virt,DBLE(F(ibnd_occ,ibnd_virt,ip)),AIMAG(F(ibnd_occ,ibnd_virt,ip)),100d0*abs(sum_f)       
       enddo
      enddo
     endif
     !
  END DO
  !
  !  if(lr_verbosity>2) call lr_diagonalise(itermax-1)
  !
  WRITE(stdout,'(5x,"End of Lanczos iterations")') 
  !
  !Final reports
  !
  if (project .and. n_ipol == 3) then
      write(stdout,'(/,/5x,"Participation of virtual states to absorbtion coefficent")') 
      write(stdout,'(5x,"occ",1x,"con",5x,"Re(Tr(F.R))",6x,"Im(TR(F.R))",5x,"% in alpha")')
      do ibnd_occ=1,nbnd
       do ibnd_virt=1,(nbnd_total-nbnd)
       sum_f=cmplx(0.0d0,0.0d0,dp)
       sum_c=cmplx(0.0d0,0.0d0,dp)
        do ip=1,n_ipol
         sum_f=sum_f+F(ibnd_occ,ibnd_virt,ip)*R(ibnd_occ,ibnd_virt,ip)
         sum_c=sum_c+chi(ip,ip)
        enddo
        write(stdout,'(5x,i3,1x,i3,3x,E16.8,2X,E16.8,2X,"%",F5.2)') & 
       ibnd_occ,ibnd_virt,DBLE(sum_F),AIMAG(sum_F),100d0*AIMAG(sum_f)/AIMAG(sum_c)
       enddo
      enddo
  endif

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
 LOGICAL FUNCTION test_restart()
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!This tests whether the restart flag is applicable
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 use lr_variables,     only : n_ipol,LR_polarization,restart
 use io_files,         only: prefix, tmp_dir, nd_nmbr
 USE mp,               ONLY : mp_bcast, mp_barrier,mp_sum
 USE io_global,        ONLY : ionode, ionode_id

 IMPLICIT NONE
  character(len=256) :: tempfile, filename
  logical :: exst
  character(len=6), external :: int_to_char
  integer :: i, temp_restart

 !
 temp_restart=0
 !print *, "test_restart with restart=",restart
 if (.not.restart) then 
  test_restart = .false.
  return  
 endif
 test_restart=.true.
 if ( n_ipol == 1 ) then
  filename = trim(prefix)//'.d0psi.'//trim(int_to_char(1))
  tempfile = trim(tmp_dir) // trim(filename) //nd_nmbr 
  inquire (file = tempfile, exist = exst)
  !print *, tempfile," exst=",exst
  if (.not. exst) then
    temp_restart=1
  endif
 else 
  DO i=1, n_ipol
   filename = trim(prefix)//'.d0psi.'//trim(int_to_char(i))
   tempfile = trim(tmp_dir) // trim(filename) //nd_nmbr
   inquire (file = tempfile, exist = exst)
   !print *, tempfile," exst=",exst
   if (.not. exst) then
     temp_restart=1
   endif
  END DO
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
    temp_restart=1
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
    temp_restart=1
 endif
 
  !print *,"temp_restart",temp_restart
#ifdef __PARA
    call mp_sum(temp_restart)
#endif
  !print *, "current temp_restart", temp_restart
  if (temp_restart > 0 ) then 
   !print *,"restart falsified",nd_nmbr
   !WRITE(stdout,'(5X,A,3X,"is missing, unable to restart.")') offender
   WRITE(stdout,'(5X,"Some files are missing, unable to restart current loop")')
   test_restart=.false.
  endif
  
  RETURN 
 END FUNCTION test_restart
END PROGRAM lr_main
!-----------------------------------------------------------------------
