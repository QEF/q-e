subroutine read_input_and_bcast(filerecon, r_paw)
  USE kinds, ONLY : DP
  use xspectra
  USE cell_base,       ONLY :  at
  USE io_files,      ONLY : prefix, tmp_dir
  USE cut_valence_green, ONLY :&
       cut_ierror, &    ! convergence tolerance for one step in the integral
       cut_stepu , &    ! integration initial step, upper side
       cut_stepl , &    ! integration initial step, lower side
       cut_startt, &    ! integration start value of the t variable
       cut_tinf  , &    ! maximum value of the lower integration boundary
       cut_tsup  , &    ! minimum value of the upper integration boudary
       cut_desmooth,&   ! size of the interval near the fermi energy
                        ! in which cross section is smoothed
       cut_nmemu,&      ! size of the memory of the values of the green function, upper side
       cut_nmeml,&      ! size of the memory of the values of the green function, lower side
       cut_occ_states  ! true if you want tou remove occupied states from the spectrum
  USE gamma_variable_mod, ONLY : gamma_value, gamma_energy, &
                                 gamma_lines, gamma_tab, gamma_points, &
                                 gamma_mode, gamma_file
  USE io_global,       ONLY : stdout,ionode,ionode_id   ! Modules/io_global.f90
  USE mp,              ONLY : mp_bcast, mp_sum             !parallelization
  USE mp_world,        ONLY : nproc, world_comm
  USE parameters,      ONLY : ntypx,lmaxx,lqmax
  USE control_flags, ONLY : twfcollect
  USE klist, ONLY : nelup, neldw, nelec

  IMPLICIT NONE

  INTEGER :: nargs, iiarg, ierr, ios, i
  LOGICAL :: found ! input_file found or not ?
  REAL(DP) :: norm, xeps_dot_xk
  REAL(DP) :: r_paw(0:lmaxx)
  CHARACTER (LEN=256) :: input_file
  CHARACTER (LEN=256) :: filerecon(ntypx)


  ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  ! $   Namelists Definition
  ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  namelist / input_xspectra / &
       calculation,&       !
       verbosity, &        ! high/low
       prefix, &           ! prefix of the pwscf output files
       outdir, &           ! directory tmp_dir or where the files are
       xiabs,&
       xkvec,&
       xepsilon,&
       xcoordcrys,&
       ef_r,&              ! obsolete since June 2014
       xe0,&            ! Zero of energy for cross section plot in eV
       xonly_plot,&
       xread_wf,&
       x_save_file,&
       xniter,&
       xerror,&
       xcheck_conv, &
       show_status, &
       nelup,neldw, &
       wf_collect,&
       U_projection_type,&
       time_limit,&
       restart_mode,&
       edge,   &            ! 'K', 'L2' or 'L3'
       lplus,   &            !  if true only the l+1 transition is calculated for L23
       lminus            !  if true only the l-1 transition is calculated for L23


  namelist / plot / &
       xnepoint,&
       xgamma,&
       xemax,&
       xemin,&
       cut_occ_states,&
       terminator,&
       gamma_mode,&
       gamma_file,&
       gamma_energy,&
       gamma_value

  namelist / pseudos /&
       filerecon,&
       filecore,&
       r_paw

  namelist / cut_occ /&
       cut_ierror, cut_stepu, cut_stepl, cut_startt, cut_tinf, cut_tsup,&
       cut_desmooth, cut_nmemu, cut_nmeml


  CALL set_xspectra_namelists_defaults()

  ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  ! $   Check if the input is from file or from stdin and read it
  ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  IF ( ionode ) THEN

     ! This part is similar to subroutine input_from_file (in flib/inpfile.f90)

     nargs = command_argument_count()
     found = .FALSE.
     input_file = ' '
    
     DO iiarg = 1, (nargs-1)
       !
       CALL get_command_argument( iiarg, input_file )
       IF ( TRIM( input_file ) == '-input' .OR. &
            TRIM( input_file ) == '-inp'   .OR. &
            TRIM( input_file ) == '-in'    .OR. &
            TRIM( input_file ) == '-i' ) THEN
          !
          CALL get_command_argument( ( iiarg + 1 ) , input_file )
          found = .TRUE.
          EXIT
       ENDIF
       !
     ENDDO

     IF (found) THEN
       OPEN ( UNIT = 5, FILE = input_file, FORM = 'FORMATTED', &
              STATUS = 'OLD', IOSTAT = ierr )
       IF ( ierr > 0 ) THEN
         CALL errore('iosys', '    input file ' // TRIM( input_file ) // &
                & ' not found' , ierr )
!        WRITE (stdout, '(/,5x,"*** STOP: input file ",A," not found ***",/)' ) &
!          TRIM( input_file )
       ENDIF
     ELSE
       ierr = -1
     END IF

     ! ... Reading namelists
     WRITE(stdout,1000) ! return+line
     WRITE(stdout,'(5x,a,a)')  &
       '                           Reading ','input_file'
     WRITE(stdout,1001) ! line+return

     READ(5, input_xspectra, err = 200, iostat = ios)
200  CALL errore ('input_xspectra', 'reading input_xspectra namelist', abs (ios) )

     READ(5, plot, err = 300, iostat = ios)
300  CALL errore ('plot', 'reading plot namelist', abs (ios) )

     READ(5, pseudos, err = 400, iostat = ios)
400  CALL errore ('pseudos', 'reading pseudos namelist', abs (ios) )

     READ(5, cut_occ, err = 500, iostat = ios)
500  CALL errore ('cut_occ', 'reading cut_occ namelist', abs (ios) )


     tmp_dir = TRIM(outdir)

  ENDIF

  ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  ! $   Variables broadcasting
  ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  CALL mp_bcast( calculation, ionode_id, world_comm )
  CALL mp_bcast( edge, ionode_id, world_comm )
  CALL mp_bcast( two_edges, ionode_id, world_comm )
  CALL mp_bcast( lplus, ionode_id, world_comm )
  CALL mp_bcast( lminus, ionode_id, world_comm )
  CALL mp_bcast( tmp_dir, ionode_id, world_comm )
  CALL mp_bcast( prefix,  ionode_id, world_comm )
  CALL mp_bcast( outdir,  ionode_id, world_comm )
  CALL mp_bcast( xnepoint,  ionode_id, world_comm )
  CALL mp_bcast( xniter,  ionode_id, world_comm )
  CALL mp_bcast( xcheck_conv,  ionode_id, world_comm )
  CALL mp_bcast( xang_mom,  ionode_id, world_comm )
  CALL mp_bcast( xgamma,  ionode_id, world_comm )
  CALL mp_bcast( xerror,  ionode_id, world_comm )
  CALL mp_bcast( xemin,  ionode_id, world_comm )
  CALL mp_bcast( xemax,  ionode_id, world_comm )
  CALL mp_bcast( show_status, ionode_id, world_comm)
  CALL mp_bcast( verbosity, ionode_id, world_comm)

  CALL mp_bcast( xkvec,  ionode_id, world_comm )
  CALL mp_bcast( xepsilon,  ionode_id, world_comm )

  CALL mp_bcast( xonly_plot,  ionode_id, world_comm )
  CALL mp_bcast( filerecon,  ionode_id, world_comm )
  CALL mp_bcast( filecore,  ionode_id, world_comm )
  CALL mp_bcast( xiabs,  ionode_id, world_comm )
  CALL mp_bcast( r_paw,  ionode_id, world_comm )
  CALL mp_bcast( xread_wf,  ionode_id, world_comm )
  CALL mp_bcast( x_save_file,  ionode_id, world_comm )
  CALL mp_bcast( xcoordcrys,  ionode_id, world_comm )
  CALL mp_bcast( ef_r,  ionode_id, world_comm )
  CALL mp_bcast( xe0,  ionode_id, world_comm )
  CALL mp_bcast( cut_occ_states, ionode_id, world_comm )
  CALL mp_bcast( terminator, ionode_id, world_comm )
  CALL mp_bcast( wf_collect, ionode_id, world_comm )
  CALL mp_bcast( twfcollect, ionode_id, world_comm )

  CALL mp_bcast( U_projection_type, ionode_id, world_comm )

  CALL mp_bcast( gamma_mode, ionode_id, world_comm )
  CALL mp_bcast( gamma_energy, ionode_id, world_comm )
  CALL mp_bcast( gamma_value, ionode_id, world_comm )

  CALL mp_bcast( cut_ierror, ionode_id, world_comm )
  CALL mp_bcast( cut_stepu, ionode_id, world_comm )
  CALL mp_bcast( cut_stepl, ionode_id, world_comm )
  CALL mp_bcast( cut_startt, ionode_id, world_comm )
  CALL mp_bcast( cut_tinf, ionode_id, world_comm )
  CALL mp_bcast( cut_tsup, ionode_id, world_comm )
  CALL mp_bcast( cut_desmooth, ionode_id, world_comm )
  CALL mp_bcast( cut_nmemu, ionode_id, world_comm )
  CALL mp_bcast( cut_nmeml, ionode_id, world_comm )

  !... restart
  CALL mp_bcast( time_limit, ionode_id, world_comm )
  CALL mp_bcast( restart_mode, ionode_id, world_comm )


 1000 FORMAT(/,5x,&
  '-------------------------------------------------------------------------')
 1001 FORMAT(5x,&
  '-------------------------------------------------------------------------',&
  /)


  ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  ! $ check on input variables          $
  ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$


  IF(TRIM(ADJUSTL(calculation)).EQ.'xanes_dipole') THEN
     xang_mom=1                    !so it is not necessary to specify xang_mom
     calculation='xanes'
  ELSEIF(TRIM(ADJUSTL(calculation)).EQ.'xanes_quadrupole') THEN
     xang_mom=2                    !so it is not necessary to specify xang_mom
     calculation='xanes'
  ENDIF




  ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  ! $   check on wfcollect
  ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  IF(xread_wf.AND.wf_collect) THEN
     CALL errore ('main','incompatibility xread_wf and wf_collect',1)
  ENDIF

  twfcollect=wf_collect

  IF(trim(adjustl(edge)).NE.'K'.AND. &
       trim(adjustl(edge)).NE.'L1'.AND. &
       trim(adjustl(edge)).NE.'L2'.AND. &
       trim(adjustl(edge)).NE.'L3'.AND. &
       trim(adjustl(edge)).NE.'L23') then
     write(stdout,*) 'Calculation for this edge not implemented!'
     write(stdout,*) 'Program is stopped'
     call stop_xspectra()
  ENDIF

  IF(trim(adjustl(edge)).eq.'L23') then
     write(stdout,*) 'Calculation of either L2 or L3'
     write(stdout,*) 'Please choose the one you want'
     call stop_xspectra()
  ENDIF

  IF(xang_mom.eq.2.and.         &
       (trim(adjustl(edge)).eq.'L2'.or.trim(adjustl(edge)).eq.'L3'&
       .or.trim(adjustl(edge)).eq.'L23') ) then
     write(stdout,*) 'Quadrupolar cross section for L23 edges not implemented.'
     write(stdout,*) 'Program is stopped'
     call stop_xspectra()
  ENDIF

 ! 
 !   lplus and lminus cannot be both positive
 !
  if( lplus .and. lminus ) then
    lplus  = .false.
    lminus = .false.
  end if

  



end subroutine read_input_and_bcast
