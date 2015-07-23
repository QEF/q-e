!
! Copyright (C) 2009-2010 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
PROGRAM X_Spectra
  USE kinds, ONLY : DP
  USE constants,       ONLY : rytoev,pi,fpi
  USE io_global,       ONLY : stdout,ionode,ionode_id   ! Modules/io_global.f90
  USE io_files,        ONLY : prefix, iunwfc, nwordwfc, tmp_dir, diropn
  USE cell_base,       ONLY : bg, at, celldm
  USE parameters,      ONLY : ntypx,lmaxx,lqmax
  USE ions_base,       ONLY : nat, ntyp => nsp, ityp, tau
  USE ktetra,          ONLY : ltetra, ntetra, tetra
  USE start_k,         ONLY : nk1, nk2, nk3, k1, k2, k3
  USE wvfct,           ONLY : npwx,nbnd,npw,igk,et, wg! et(nbnd,nkstot)
  USE radial_grids,    ONLY : ndmx
  USE atom,            ONLY : rgrid
  USE becmod,          ONLY : becp
  USE uspp,            ONLY : vkb, nkb, okvan 
  USE uspp_param,      ONLY : upf
  USE xspectra
  USE ener,            ONLY : ef, ef_up, ef_dw !Fermi energy (ef in Ry)
  USE symm_base,       ONLY : nsym,s
  USE paw_gipaw,       ONLY : read_recon,  &
       paw_vkb,             & ! |p> projectors
       paw_becp,            & ! product of projectors and wf.
       paw_nkb,             & ! total number of beta functions, with st.fact.
       paw_lmaxkb,          &
       paw_recon,           &
       set_paw_upf
  USE klist,           ONLY : &
       nkstot,            & ! total number of k-points
       nks,               & ! number of k-points for local pool
       nelec,nelup,neldw,             & !number of electrons
       xk,                & ! k-points coordinates
       wk ,               & ! k-points weight
       npk,               &
       degauss,lgauss,ngauss,    &
       two_fermi_energies
  USE lsda_mod,        ONLY : nspin,lsda,isk,current_spin
  USE noncollin_module,ONLY : noncolin
  USE mp,              ONLY : mp_bcast, mp_sum             !parallelization
  USE mp_global,       ONLY : mp_startup, mp_global_end
  USE mp_pools,        ONLY : intra_pool_comm, npool
  USE mp_world,        ONLY : nproc, world_comm
  USE control_flags,   ONLY : gamma_only
  USE environment,     ONLY : environment_start

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

  USE control_flags,   ONLY : twfcollect
  !<CG>
  USE gamma_variable_mod, ONLY : gamma_value, gamma_energy, &
                                 gamma_lines, gamma_tab, gamma_points, &
                                 gamma_mode, gamma_file
  USE xspectra_paw_variables, ONLY : xspectra_paw_nhm, init_xspectra_paw_nhm
  USE edge_energy, ONLY: getE
  !</CG>

  IMPLICIT NONE 
  !
  ! ... local variables
  !
  !<DC> clean the declarations
  INTEGER :: nargs,iiarg,ierr,ios,il,ibnd,ibnd_up,ibnd_dw !,xm_r,nc_r,ncomp_max
  INTEGER :: iargc
  !INTEGER :: nt,nb,na,i,j,k,nrest,nline
  INTEGER :: nt,na,i,j,k !nrest,nline
  INTEGER, ALLOCATABLE :: ncalcv(:,:)
  INTEGER, ALLOCATABLE :: paw_iltonhb(:,:,:)
              ! corresp l, projector, type <--> cumulative over all the species
  REAL (DP) :: norm, core_energy
  REAL (DP) :: rc(ntypx,0:lmaxx),r_paw(0:lmaxx)
  REAL (DP) :: core_wfn(ndmx)
  REAL (DP) :: ehomo, elumo, middle_gap ! in eV 
  REAL (DP) :: e_core, e_core_ryd
  REAL (DP), ALLOCATABLE :: a(:,:,:),b(:,:,:),xnorm(:,:)      !lanczos vectors
  !REAL (DP), EXTERNAL   :: efermig,efermit
 
  LOGICAL :: exst !, loc_set
  LOGICAL :: terminator, show_status, wf_collect 
  LOGICAL :: found ! input_file found or not ?

  CHARACTER (LEN=256) :: input_file, filerecon(ntypx), filecore 
  CHARACTER (LEN=256) :: outdir
  CHARACTER (LEN=25)  :: calculation 
  CHARACTER (LEN=4)   :: verbosity
  CHARACTER (LEN=10)  :: dummy_char
  !REAL (DP) :: gamma_energy(2), gamma_value(2)
  REAL (DP) :: xeps_dot_xk ! scalar product between xepsilon and xkvec
  !REAL(dp) :: auxrpol(3,2) non used variable
  !</DC>

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


  ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  ! $   Initialize MPI environment, clocks and a few other things
  ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

#ifdef __MPI
  CALL mp_startup ( )
#endif
  CALL environment_start ( 'XSpectra' )
  CALL banner_xspectra()

  ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  ! $   Set default values for some namelist variables
  ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  ! ... Namelist input_xspectra
  !calculation='xanes'
  calculation='xanes_dipole'
  prefix=' '
  verbosity='low'
  x_save_file='xanes.sav'
  CALL get_env( 'ESPRESSO_TMPDIR', outdir )
  IF ( TRIM( outdir ) == ' ' ) outdir = './'
  xniter=50
  xcheck_conv=5
  xonly_plot=.FALSE.
  xread_wf=.FALSE.
  xerror=1.d-2
  xiabs=1               !identify the adsorbing atom
  !<DC>
  !DO i=2,3
  !   xkvec(i)=0.d0
  !   xepsilon(i)=0.d0
  !ENDDO
  xkvec(1)=0.d0       
  xkvec(2)=1.d0
  xkvec(3)=0.d0
  xepsilon(1)=1.d0
  xepsilon(2:3)=0.d0
  !ef_r=0.d0
  xe0=xe0_default 
  !</DC>
  xcoordcrys=.true.
  show_status=.false.
  wf_collect=.false.
  U_projection_type='atomic'
  restart_mode='from_scratch'
  time_limit=1.d8
  edge='K'
  lplus=.false.
  lminus=.false.

 
  ! ... Namelist plot
  xnepoint=100
  xemin=0.d0
  xemax=10.d0
  xgamma=0.1d0
  cut_occ_states=.FALSE.
  terminator=.false.
  gamma_mode='constant'
  gamma_file='gamma.dat'
  

  ! ... Namelist pseudos 
  filecore='Core.wfc'

  ! ... Namelist cut_occ (for cutting the occupied states, paste_fermi function)
  cut_ierror=1.d-7
  cut_stepu=1.d-2
  cut_stepl=1.d-3
  cut_startt=1.d0
  cut_tinf=1.d-6
  cut_tsup=100.d0
  cut_desmooth=1.d-2
  cut_nmemu=100000
  cut_nmeml=100000

  ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  ! $   Check if the input is from file or from stdin and read it
  ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  IF ( ionode ) THEN

     ! This part is similar to subroutine input_from_file (in flib/inpfile.f90)

     nargs = iargc()
     found = .FALSE.
     input_file = ' '
     
     DO iiarg = 1, (nargs-1)
       !
       CALL getarg( iiarg, input_file )
       IF ( TRIM( input_file ) == '-input' .OR. &
            TRIM( input_file ) == '-inp'   .OR. &
            TRIM( input_file ) == '-in'    .OR. &
            TRIM( input_file ) == '-i' ) THEN
          !
          CALL getarg( ( iiarg + 1 ) , input_file )
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


 !  
 !   lplus and lminus cannot be both positive
 !
  if( lplus .and. lminus ) then
    lplus  = .false.
    lminus = .false.
  end if


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

  ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  ! $   Writes the relevant parameters of the XSpectra calculation
  ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  WRITE(stdout, '(5x,a,a,/)') 'calculation: ', TRIM(ADJUSTL(calculation))

!<NM> Write xepsilon and (if needed) xkvec
 IF(.NOT.xonly_plot) THEN
     IF ( xcoordcrys ) THEN
         WRITE(stdout,'(5x,a,3(f10.6,1x),/)') &
         'xepsilon  [crystallographic coordinates]: ', (xepsilon(i),i=1,3)
     ELSE
        WRITE(stdout,'(5x,a,3(f10.6,1x),/)') &
         'xepsilon  [cartesian coordinates]: ', (xepsilon(i),i=1,3)
     ENDIF

     IF ( TRIM(ADJUSTL(calculation)).EQ.'xanes_quadrupole' )  THEN
        IF ( xcoordcrys ) THEN
            WRITE(stdout,'(5x,a,3(f10.6,1x),/)') &
            'xkvec  [crystallographic coordinates]: ', (xkvec(i),i=1,3)
        ELSE
           WRITE(stdout,'(5x,a,3(f10.6,1x),/)') &
           'xkvec [cartesian coordinates]: ', (xkvec(i),i=1,3)
        ENDIF
     ENDIF
!<NM>
ENDIF


  ! ... Writes xonly_plot, its meaning and plot parameters
 
  IF (xonly_plot.EQV..FALSE.) then
     WRITE(stdout,'(5x,a)') 'xonly_plot: FALSE'
     WRITE(stdout,'(8x,a,/)') &  
          '=> complete calculation: Lanczos + spectrum plot'
     WRITE(stdout,'(5x,a,a20)') 'filecore (core-wavefunction file): ', &
                                 filecore
  ELSE
     WRITE(stdout,'(5x,a)') 'xonly_plot: TRUE'
     WRITE(stdout,'(8x,a)') & 
          '=> only the spectrum plot'
  ENDIF
  WRITE(stdout,*) 

  WRITE(stdout,'(5x,a)') 'main plot parameters:'
  IF (cut_occ_states) THEN
     WRITE(stdout,'(8x,a)') 'cut_occ_states: TRUE'
  ELSE
     WRITE(stdout,'(8x,a)') 'cut_occ_states: FALSE'
  ENDIF
  WRITE(stdout,'(8x,a,a8)')   'gamma_mode:  ', gamma_mode 
  IF (TRIM(ADJUSTL(gamma_mode)).EQ.'constant') THEN
    WRITE(stdout,'(8x,a,f5.2)') '-> using xgamma [eV]: ', xgamma
  ELSEIF (TRIM(ADJUSTL(gamma_mode)).EQ.'file') THEN
    WRITE(stdout,'(8x,a,a50)')  '-> using gamma_file: ', gamma_file
  ELSEIF (TRIM(ADJUSTL(gamma_mode)).EQ.'variable') THEN
    WRITE(stdout,'(8x,a,f5.2,a1,f5.2,a)') &
     '-> first, constant up to point (', &
                                     gamma_energy(1),',',gamma_value(1),') [eV]' 
    WRITE(stdout,'(8x,a,f5.2,a1,f5.2,a)') &
     '-> then, linear up to point (',gamma_energy(2),',',gamma_value(2),') [eV]' 
    WRITE(stdout,'(8x,a)') '-> finally, constant up to xemax'
  ENDIF
  WRITE(stdout,'(8x,a,f6.2)') 'xemin [eV]: ', xemin
  WRITE(stdout,'(8x,a,f6.2)') 'xemax [eV]: ', xemax
  WRITE(stdout,'(8x,a,i4)')   'xnepoint: ', xnepoint
  IF (abs(xe0-xe0_default)<1.d-3) THEN
     WRITE(stdout,'(8x,a,/)') &
          'energy zero automatically set to the Fermi level'
     IF (xonly_plot) THEN
        WRITE(stdout,'(5x,a)') 'Fermi level read in x_save_file'
     ELSE
        WRITE(stdout,'(5x,3a)') &
          'Fermi level determined from SCF save directory (', &
          trim(prefix)//'.save',')'
     ENDIF
     WRITE(stdout,'(5x,a)') &
         'NB: For an insulator (SCF calculated with occupations="fixed")'
     WRITE(stdout,'(5x,a)') &
         '    the Fermi level will be placed at the position of HOMO.'
  ELSE
     WRITE(stdout,'(8x,a,f10.6,3a)') 'xe0 [eV]: ', xe0, &
         ' (energy zero read in ','input file',')'
  ENDIF
  WRITE(stdout,*) 
  WRITE(stdout,'(5x,a)') 'WARNING: variable ef_r is obsolete'
  !</DC>

  ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  ! $   Writing the status of the code (working features and things to do)
  ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  IF(show_status) CALL WRITE_status_of_the_code

  ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  ! $   Initialising calculation and clock
  ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  IF(edge.EQ.'L1') edge='K'

  IF(TRIM(ADJUSTL(calculation)).EQ.'xanes_dipole') THEN
     xang_mom=1                    !so it is not necessary to specify xang_mom
     calculation='xanes'
  ELSEIF(TRIM(ADJUSTL(calculation)).EQ.'xanes_quadrupole') THEN
     xang_mom=2                    !so it is not necessary to specify xang_mom
     calculation='xanes'
  ENDIF

  call select_nl_init(edge, nl_init, two_edges, n_lanczos)     

  CALL start_clock( calculation  )

  ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  ! $   check on wfcollect
  ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  IF(xread_wf.AND.wf_collect) THEN
     CALL errore ('main','incompatibility xread_wf and wf_collect',1)
  ENDIF

  twfcollect=wf_collect

  ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  ! $   Reads, initializes and writes several things in two cases :
  ! $   case 1: xonlyplot=.false.  (complete calc., i.e. Lanczos + plot)
  ! $   case 2: xonlyplot=.true.   (the plot only)
  ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  IF(.NOT.xonly_plot) THEN

     WRITE(stdout,1000) ! return+line
     WRITE(stdout,'(5x,2a)')  &
     '                 Reading SCF save directory: ',trim(prefix)//'.save'
     WRITE(stdout,1001) ! line+return

     CALL read_file()
     ehomo=0.d0
     elumo=0.d0
     CALL get_homo_lumo(ehomo,elumo)
     ehomo = ehomo*rytoev
     elumo = elumo*rytoev

     call reset_k_points_and_reinit_nscf()

 

     IF(xread_wf) THEN
        IF (okvan) THEN
           WRITE(stdout,'(/,5x,a,f10.6,/)') &
                        'Approx. ram memory needed per proc in MB = ',&
                        (16.0*4.0*npwx*npool)/(nproc*1000000.0)
        ELSE
           WRITE(stdout,'(/,5x,a,f10.6,/)') &
                        'Approx. ram memory needed per proc in MB = ',&
                        (16.0*3.0*npwx*npool)/(nproc*1000000.0)
        ENDIF
     ENDIF

 
        ! ... normalize xepsilon 

     IF ( xcoordcrys ) CALL cryst_to_cart(1,xepsilon,at,1)
     norm=DSQRT(xepsilon(1)**2+xepsilon(2)**2+xepsilon(3)**2)
     DO i=1,3
       xepsilon(i)=xepsilon(i)/norm
     ENDDO
 
     ! ... If needed normalize xkvec 
     !     and check orthogonality with xepsilon 

     IF(xang_mom.EQ.2) THEN
     
        IF ( xcoordcrys ) CALL cryst_to_cart(1,xkvec,at,1)
        norm=DSQRT(xkvec(1)**2+xkvec(2)**2+xkvec(3)**2)
        xkvec(:)=xkvec(:)/norm
        xeps_dot_xk=xkvec(1)*xepsilon(1)+&
                    xkvec(2)*xepsilon(2)+& 
                    xkvec(3)*xepsilon(3) 
        IF ((ABS(xeps_dot_xk)) >= 1.0d-6) THEN
           WRITE(stdout,'(5x,a)') &
              'ERROR: xkvec and xepsilon are not orthogonal'
           WRITE(stdout,'(12x,a,f10.6,/)') 'scalar product=', xeps_dot_xk
           WRITE(stdout,'(5x,a)') 'STOP'
           CALL stop_xspectra ()
        ENDIF

     ENDIF


     !... Is the type associated to xiabs existing ?
     
     i=0
     DO na=1,nat
        IF(ityp(na).EQ.xiabs) i=i+1
     ENDDO
     IF(i.NE.1) THEN
        CALL errore( 'main program', 'Wrong xiabs!!!',i)
     ENDIF

     !...  Reads reconstruction files
     WRITE(stdout,1000) ! return+line 
     WRITE(stdout,'(5x,a)')  &
     '          Reading core wavefunction file for the absorbing atom'
     WRITE(stdout,1001) ! line+return
     DO nt = 1, ntyp
        call set_paw_upf(nt, upf(nt))
     ENDDO

     CALL read_core_abs(filecore,core_wfn, nl_init)

     WRITE(stdout,'(5x,a," successfully read")') TRIM(ADJUSTL(filecore))
     
     IF ( .NOT. paw_recon(xiabs)%gipaw_data_in_upf_file ) &
          CALL read_recon ( filerecon(xiabs), xiabs, paw_recon(xiabs) ) !*apsi

     !... Assign paw radii to species (this will become soon obsolete)
     ! <DC> Why? </DC>
     !
     !  read recon should be parallelized
     !
     WRITE(stdout,1000) ! return+line 
     WRITE(stdout,'(5x,a)')  &
     '                         Attributing the PAW radii '
     WRITE(stdout,'(5x,a)')  &
     '                for the absorbing atom [units: Bohr radius]'
     WRITE(stdout,1001) ! line+return

     DO nt=1,ntyp
        IF ((.NOT.paw_recon(nt)%gipaw_data_in_upf_file)) then
           paw_recon(nt)%paw_nbeta=0
        ENDIF
        DO  j=1,paw_recon(nt)%paw_nbeta
           il=paw_recon(nt)%psphi(j)%label%l
           ! The change made by DC is wrong, this has to be reset as
           ! it was before.
           !<DC> changed the following block (writing format)
           !IF(xiabs.EQ.nt.AND.DABS(r_paw(il)).lt.1.d-6) THEN
           !   !*apsi  if(r(kkbeta(nt),nt).GT.1.d-3) THEN
           !   !*apsi     rc(nt,il)=  r(kkbeta(nt),nt)
           !   !*apsi  ELSE
           !   !*apsi     WRITE(stdout,*) 'Warning-no-kkbeta r_paw(',il,')=1.0'
           !   !*apsi     rc(nt,il)=1.0
           !   !*apsi  ENDIF
           !   !<CG>  to be verified
           !   IF (paw_recon(nt)%psphi(j)%label%rc > 1.d-3) THEN
           !      WRITE(stdout,*) 'warning, r_paw(', il,' ) set to ', &
           !           paw_recon(nt)%psphi(j)%label%rc
           !      rc(nt, il)= paw_recon(nt)%psphi(j)%label%rc*3.0/2.0
           !   ELSE
           !      WRITE(stdout,*) 'Warning, no rc'
           !      WRITE(stdout,*) 'warning, r_paw(', il,' ) set to 1.5'
           !      rc(nt, il)= 1.5d0
           !   ENDIF
           !   !</CG>
           !ELSEIF(xiabs.EQ.nt.AND.DABS(r_paw(il)).GT.1.d-6) THEN
           !   rc(nt,il)=r_paw(il)
           !ELSEIF(nt.NE.xiabs) THEN
           !   !*apsi  IF(r(kkbeta(nt),nt).GT.1.d-3) THEN
           !   !*apsi     rc(nt,il)=r(kkbeta(nt),nt)
           !   !*apsi  ELSE
           !   !*apsi     rc(nt,il)=1.0
           !   !*apsi  ENDIF
           !   !<CG> to be verified
           !   IF(paw_recon(nt)%psphi(j)%label%rc.GT.1.d-3) THEN
           !      rc(nt,il)=paw_recon(nt)%psphi(j)%label%rc*3.0/2.0
           !   ELSE
           !      rc(nt,il)=1.5
           !   ENDIF
           !   !</CG> 
           !ENDIF
           !... assigns rc for the absorbing atom
           IF ( xiabs.EQ.nt ) then
              IF ( DABS(r_paw(il)) > 1.d-6 ) THEN
                 rc(nt,il)=r_paw(il)
                 WRITE(stdout,'(8x,a,i2,a,i2,a,f5.2,3a)') &
                 'PAW proj', j,': r_paw(l=',il,')=',rc(nt,il), &
                 '  (from ','input file)',')'
              ELSE
                 !<CG>  to be verified
                 IF (paw_recon(nt)%psphi(j)%label%rc > 1.d-3) THEN
                    rc(nt, il)= paw_recon(nt)%psphi(j)%label%rc*3.0/2.0
                    WRITE(stdout,'(8x,a,i2,a,i2,a,f5.2,a)') &
                    'PAW proj', j,': r_paw(l=',il,')=',rc(nt,il),'  (1.5*r_cut)'
                 ELSE
                    rc(nt, il)= 1.5d0
                    WRITE(stdout,'(8x,a,i2,a,i2,a,f5.2,a)') &
                    'PAW proj',j,': r_paw(l=',il,')=',rc(nt,il),'  (set to 1.5)'
                 ENDIF
                 !</CG>
              ENDIF
           !... assigns rc for the other atoms (not further used)
           ELSE
              IF (paw_recon(nt)%psphi(j)%label%rc > 1.d-3) THEN
                 rc(nt,il)=paw_recon(nt)%psphi(j)%label%rc*3.0d0/2.0d0
              ELSE
                 rc(nt,il)=1.5d0
              ENDIF
           ENDIF
           !</DC>
        ENDDO
     ENDDO
     WRITE(stdout,*)
     WRITE(stdout,'(8x,a)')&
           'NB: The calculation will not necessary use all these r_paw values.'
     WRITE(stdout,'(8x,a)')&
           '    - For a edge in the electric-dipole approximation,'
     WRITE(stdout,'(8x,a)')&
           '      only the r_paw(l=1) values are used.'
     WRITE(stdout,'(8x,a)')&
           '    - For a K edge in the electric-quadrupole approximation,'
     WRITE(stdout,'(8x,a,/)')'      only the r_paw(l=2) values are used.'
     WRITE(stdout,'(8x,a,/)')     '    - For a L2 or L3 edge in the electric-quadrupole approximation,'
     WRITE(stdout,'(8x,a,/)')'      all projectors (s, p and d) are used.'

     !<CG>
     DO nt=1,ntyp
        DO j = 1,paw_recon(nt)%paw_nbeta
           paw_recon(nt)%psphi(j)%label%rc=rc(nt,paw_recon(nt)%psphi(j)%label%l)
           paw_recon(nt)%aephi(j)%label%rc=rc(nt,paw_recon(nt)%aephi(j)%label%l)
        ENDDO
     ENDDO
     !</CG>


     !...  write band energies if xread_wf=true

     IF(xread_wf.AND.TRIM(verbosity).EQ.'high') THEN

        WRITE(stdout,'(5x,a)') &
           'Band energies read from scf save file [units: eV]'
        WRITE(stdout,'(5x,a)') &
           '-------------------------------------------------'
        DO i=1,nkstot
           WRITE(stdout,'(5x,"k=[",3f14.8,"]   spin=",1i2)') &
                xk(1,i),xk(2,i),xk(3,i),isk(i)
           WRITE(stdout, '(8f9.4)') (et(j,i)*rytoev,j=1,nbnd)
        ENDDO

        WRITE(stdout,*)

     ENDIF

     CALL init_gipaw_1
   
     !
     WRITE(stdout,1000) ! return+line
     WRITE(stdout,'(5x,a)') & 
     '                      Getting the Fermi energy '
     WRITE(stdout,1001) ! line+return
     !
     IF (lsda) THEN
       WRITE(stdout,'(5x,a,a)') 'From SCF save directory',&
                         ' (spin polarized work):' 
       !
       IF (abs(ehomo)<1.e+6) THEN ! insulator => HOMO exists
         WRITE(stdout,'(8x,a,f9.4,a)') 'ehomo [eV]: ', ehomo,&
                        ' (highest occupied level:max of up and down)'
         ef=ehomo
         IF (abs(elumo)<1.e+6) THEN  ! insulator and LUMO exists 
            WRITE(stdout,'(8x,a,f9.4,a)') 'elumo [eV]: ', elumo,&
                        ' (lowest occupied level:min of up and down)'
         ELSE
            WRITE(stdout,'(8x,a)') 'No LUMO values in SCF calculation'
         ENDIF
       ELSE IF (abs(ef)>1.e-4) THEN
         ef=ef*rytoev !ef in eV
       ELSE
         WRITE(stdout,'(8x,a,f9.4)') 'ef_up [eV]: ', ef_up*rytoeV
         WRITE(stdout,'(8x,a,f9.4)') 'ef_dw [eV]: ', ef_dw*rytoeV
         ef=max(ef_dw,ef_up)*rytoeV
         WRITE(stdout,'(8x,a,f9.4)') '-> ef set to the max of ef_up and ef_dw '
       ENDIF
       WRITE(stdout,'(8x,a,f9.4)') 'ef    [eV]: ', ef   
      
       WRITE(stdout,'(/,5x,a)') &
          '-> ef (in eV) will be written in x_save_file'
     ELSE
       WRITE(stdout,'(5x,a)') 'From SCF save directory:'
       !
       ef=ef*rytoev
       IF (abs(ehomo)<1.e+6) THEN ! insulator => HOMO exists
         WRITE(stdout,'(8x,a,f9.4,a)') 'ehomo [eV]: ', ehomo,&
                                      ' (highest occupied level)'
         ef=ehomo
         IF (abs(elumo)<1.e+6) THEN  ! insulator and LUMO exists 
           WRITE(stdout,'(8x,a,f9.4,a)') 'elumo [eV]: ', elumo,&
                                      ' (lowest occupied level)'
         ELSE
           WRITE(stdout,'(8x,a)') 'No LUMO value in SCF calculation'
         ENDIF
       ENDIF
       WRITE(stdout,'(8x,a,f9.4)') 'ef    [eV]: ', ef   
       
       WRITE(stdout,'(/,5x,a)') &
          '-> ef (in eV) will be written in x_save_file'
     
     ENDIF

     !
     WRITE(stdout,1000) ! return+line
     WRITE(stdout,'(5x,a)') &
     '                      Energy zero of the spectrum '
     WRITE(stdout,1001) ! line+return
     !
     IF (abs(xe0-xe0_default)<1.d-3) THEN ! no xe0 in input
          WRITE(stdout,'(5x,a)') &
          '-> ef will be used as energy zero of the spectrum'
     ELSE
          WRITE(stdout,'(5x,a,/,7x,3a)') &  
          '-> ef will NOT be used as energy zero of the spectrum',&
          '(because xe0 read in ', 'input file',')'  
     ENDIF
    

     !        CALL mp_bcast( ef, ionode_id )  !Why should I need this ?


     !... Allocation of variables for paw projectors
     !... and initialization of Vanderbilt and Paw projectors
     !
     ! CG : becp allocated in the xanes_dipole and xanes_quadrupole subroutines
     !     call allocate_bec_type ( nkb, nbnd, becp )

     !CALL init_gipaw_1 ! Already called above

     !...  Definition of a specific indexation to avoid M. Profeta's crazy one 
     CALL init_xspectra_paw_nhm

     ALLOCATE (paw_iltonhb(0:paw_lmaxkb,xspectra_paw_nhm, ntyp))

     CALL define_index_arrays(paw_iltonhb)

     !... Allocates PAW projectors
     ALLOCATE (paw_vkb( npwx,  paw_nkb))
     ALLOCATE (paw_becp(paw_nkb, nbnd))

  ELSEIF(xonly_plot) THEN  

     CALL read_header_save_file(x_save_file)
     nks = nkstot
     IF(lsda) THEN
        isk(1:nkstot/2)=1
        isk(nkstot/2+1:nkstot)=2
     !  wk(1:nkstot)=2.d0/nkstot
     ELSEIF(.NOT.lsda) THEN
        isk(1:nkstot)=1
     !  wk(1:nkstot)=2.d0/nkstot
     ENDIF
     wk(1:nkstot)=2.d0/nkstot
     CALL divide_et_impera( xk, wk, isk, lsda, nkstot, nks )

  ENDIF   

  !<CG>
  ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  ! $ Computing gamma tabulated values
  ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  !<MCB> THIS MUST BE CHANGED
  !<DC> why ? </DC>

  IF (TRIM(gamma_mode).EQ.'file') THEN
     CALL read_gamma_file
  ELSEIF (TRIM(gamma_mode).EQ.'variable') THEN
     gamma_lines=2
     ALLOCATE(gamma_points(2,2))
     gamma_points(1,1)=gamma_energy(1)
     gamma_points(2,1)=gamma_energy(2)
     gamma_points(1,2)=gamma_value(1)
     gamma_points(2,2)=gamma_value(2)
  ENDIF

  IF ((TRIM(gamma_mode).EQ.'file').OR.(TRIM(gamma_mode).EQ.'variable')) THEN
     ALLOCATE( gamma_tab(xnepoint))
     CALL initialize_gamma_tab
     DEALLOCATE(gamma_points)
  ENDIF

  !</CG>

  xnitermax=xniter

  !<DC> Block moved below into the previous "if xonlyplot" block
  !IF(xonly_plot) THEN
  !   CALL read_header_save_file(x_save_file)
  !   nks = nkstot
  !   WRITE(6,*) 'nks=',nks
  !   IF(lsda) THEN
  !      isk(1:nkstot/2)=1
  !      isk(nkstot/2+1:nkstot)=2
  !      wk(1:nkstot)=2.d0/nkstot
  !   ELSEIF(.NOT.lsda) THEN
  !      isk(1:nkstot)=1
  !      wk(1:nkstot)=2.d0/nkstot
  !   ENDIF
  !   CALL divide_et_impera( xk, wk, isk, lsda, nkstot, nks )
  !ENDIF

  ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  ! $  Checks PAW relations between pseudo partial waves and projector 
  ! $  (radial parts)
  ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  IF(.NOT.xonly_plot.AND.TRIM(verbosity).EQ.'high') &
      CALL check_paw_projectors(xiabs)

  ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  ! $   Allocates and initializes Lanczos variables
  ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  ALLOCATE(a(xnitermax,n_lanczos,nks))
  ALLOCATE(b(xnitermax,n_lanczos,nks))
  ALLOCATE(xnorm(n_lanczos,nks))
  ALLOCATE(ncalcv(n_lanczos,nks))

  a(:,:,:)=0.d0
  b(:,:,:)=0.d0
  xnorm(:,:)=0.d0
  ncalcv(:,:)=0

  ! for restart
  ALLOCATE(calculated(n_lanczos,nks))
  calculated(:,:)=0

  ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  ! $  And now we go...  XANES CALCULATION
  ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  IF(TRIM(calculation).EQ.'xanes') THEN
     IF(.NOT.xonly_plot) THEN ! complte calculation
        WRITE(stdout, 1000) ! line 
        WRITE(stdout,'(5x,a)')&
        '                     Starting XANES calculation'
        IF(nl_init(2).eq.0) then
           IF (xang_mom==1) WRITE(stdout,'(5x,a)')&
           '                in the electric dipole approximation'
           IF (xang_mom==2) WRITE(stdout,'(5x,a)')&
           '              in the electric quadrupole approximation'
           WRITE(stdout,1001)  ! line 
        ELSEIF(nl_init(2).eq.1) then
           WRITE(stdout,'(5x,a)')&
                '                in the electric dipole approximation'
        ENDIF
           
        IF(TRIM(ADJUSTL(restart_mode)).eq.'restart') THEN
          CALL read_header_save_file(x_save_file)
          CALL read_save_file(a,b,xnorm,ncalcv,x_save_file,core_energy)
        ENDIF

        !...  Writes information about the method

        WRITE(stdout,'(7(5x,a,/))') &
        "Method of calculation based on the Lanczos recursion algorithm",&
        "--------------------------------------------------------------",&
        "   - STEP 1: Construction of a kpoint-dependent Lanczos basis,",&
        "     in which the Hamiltonian is tridiagonal (each 'iter' ",&
        "     corresponds to the calculation of one more Lanczos vector)",&
        "   - STEP 2: Calculation of the cross-section as a continued fraction",&
        "     averaged over the k-points."

        WRITE(stdout,'(5x,"... Begin STEP 1 ...",/)')

        save_file_version = 2 ! adds ef after core_energy

        IF(nl_init(2).EQ.0.AND.xang_mom .eq. 1) THEN

           save_file_kind='xanes_dipole'
           CALL xanes_dipole(a,b,ncalcv,xnorm,core_wfn,&
                             paw_iltonhb,terminator,verbosity)

        ELSEIF(nl_init(2).EQ.0.AND.xang_mom .eq. 2) THEN

           save_file_kind='xanes_quadrupole'
           CALL xanes_quadrupole(a,b,ncalcv,xnorm,core_wfn,&
                                 paw_iltonhb,terminator,verbosity)

        ELSEIF(nl_init(2).EQ.1.AND.xang_mom .eq. 1) then
           call xanes_dipole_general_edge(a,b,ncalcv,nl_init,xnorm,core_wfn,paw_iltonhb,terminator, verbosity)
        ENDIF
        !
        ! write_save_file should be changed for L2,3
        !
        CALL write_save_file(a,b,xnorm,ncalcv,x_save_file)
        WRITE(stdout,'(/,5x,a)') &
              'Results of STEP 1 successfully written in x_save_file'
        WRITE(stdout,'(5x,a18,/,5x,a2,2x,a65)') 'x_save_file name: ',&
                      '->', x_save_file 
        WRITE(stdout,'(5x,a21,i2)') 'x_save_file version: ', save_file_version

        WRITE(stdout,'(/,5x,"... End STEP 1 ...",/)')

     ELSE ! only the spectrum plot
        WRITE(stdout,1000) ! return+line
        WRITE(stdout,'(5x,a)')&
        '                          Reading x_save_file'
        WRITE(stdout,1001) ! line+return

        CALL read_save_file(a,b,xnorm,ncalcv,x_save_file,core_energy)

        IF (save_file_version .eq. 1 .AND. abs(xe0-xe0_default)<1.e-3) THEN
          WRITE(stdout,'(5x,3a)') &
            "STOP: the variable 'xe0' must be assigned in ", &
             'input file', ' (since save_file_version = 1)'
          CALL stop_xspectra()
        ENDIF
     ENDIF

     IF (TRIM(save_file_kind).eq.'unfinished') CALL stop_xspectra ()       
 
     IF (.NOT.xonly_plot) THEN
        !
        WRITE(stdout,'(5x,"... Begin STEP 2 ...",/)')
        WRITE(stdout,'(5x,a)') &
         'The spectrum is calculated using the following parameters:'
        !
        e_core = getE(upf(xiabs)%psd,edge)

     ELSE
        WRITE(stdout,1000) ! return+line
        WRITE(stdout,'(5x,a,a)')  &
        '               Starting the calculation of the spectrum'
        WRITE(stdout,1001) !line+return
        WRITE(stdout,'(5x,a)') 'Using the following parameters:'
        !
        e_core = core_energy
        !
     ENDIF

     e_core_ryd = e_core/rytoev

     IF (abs(xe0-xe0_default)<1.d-3) THEN ! xe0 not in input_file
        write(stdout,'(8x,a,f9.4,a)') 'energy-zero of the spectrum [eV]: ',&
          ef 
        ! ef is either read in x_save_file (version 2) if xonlyplot=T
        ! or determined from the scf calculation
        xe0_ry = ef/rytoev    ! ef est en eV ici
     ELSE ! xe0 is read in the input
        WRITE(stdout,'(8x,a,f9.4,a)') 'xe0 [eV]: ', xe0
        xe0_ry=xe0/rytoev
     ENDIF

     IF (cut_occ_states) THEN
        WRITE(stdout,'(8x,a)') 'the occupied states are cut'
     ELSE
        WRITE(stdout,'(8x,a)') 'the occupied states are NOT cut'
     ENDIF
     WRITE(stdout,'(8x,a,f6.2)') 'xemin [eV]: ', xemin
     WRITE(stdout,'(8x,a,f6.2)') 'xemax [eV]: ', xemax
     WRITE(stdout,'(8x,a,i4)')   'xnepoint: ', xnepoint
     IF (TRIM(ADJUSTL(gamma_mode)).EQ.'constant') THEN
        WRITE(stdout,'(8x,a,f8.3)')'constant broadening parameter [eV]: ', xgamma
     ELSE
        WRITE(stdout,'(8x,a)') 'energy-dependent broadening parameter:'
        IF (TRIM(ADJUSTL(gamma_mode)).EQ.'file') THEN
           WRITE(stdout,'(8x,a,a30)')' -> using gamma_file: ', gamma_file
        ELSEIF (TRIM(ADJUSTL(gamma_mode)).EQ.'variable') THEN
           WRITE(stdout,'(8x,a,f5.2,a1,f5.2,a)')                     &
                 ' -> first, constant up to point (', gamma_energy(1), &
                 ',', gamma_value(1), ') [eV]'
           WRITE(stdout,'(8x,a,f5.2,a1,f5.2,a)')                     &
                 ' -> then, linear up to point (', gamma_energy(2),   &
                 ',', gamma_value(2), ') [eV]'
           WRITE(stdout,'(8x,a)') ' -> finally, constant up to xemax'
        ENDIF
     ENDIF
     WRITE(stdout,'(8x,"Core level energy [eV]:",1x,g11.4)') -e_core
     WRITE(stdout,'(8x,a,/)') &
     ' (from electron binding energy of neutral atoms in X-ray data booklet)'
     
     IF(xang_mom.EQ.1) THEN
        CALL plot_xanes_dipole(a,b,xnorm,ncalcv,terminator,e_core_ryd,1)
        if(two_edges) CALL plot_xanes_dipole(a,b,xnorm,ncalcv,terminator,e_core_ryd,2)
     ELSEIF(xang_mom.EQ.2) THEN
        CALL plot_xanes_quadrupole(a,b,xnorm,ncalcv,terminator,e_core_ryd)
     ELSEIF(nl_init(2).eq.1) then
        
     ENDIF
     !
     WRITE(stdout,'(5x,"Cross-section successfully written in ",a,/)') &
      'xanes.dat'
     IF (.NOT. xonly_plot) WRITE(stdout,'(5x,"... End STEP 2 ...",/)')

  ELSEIF(TRIM(calculation).EQ.'rxes') THEN
     CALL errore( 'Main', 'rxes Not yet implemented',1)
  ELSEIF(TRIM(calculation).EQ.'bethe_salpeter') THEN
     CALL errore( 'Main', 'bethe_salpeter Not yet implemented',1)
  ENDIF

  ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  !  Deallocation
  ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  !  IF(.NOT.xonly_plot) THEN
  !     call deallocate_bec_type ( becp )
  !  ENDIF

  DEALLOCATE(a)
  DEALLOCATE(b)
  DEALLOCATE(xnorm)
  DEALLOCATE(ncalcv)


  !WRITE (stdout,*) 'End program ', TRIM(calculation)

  CALL stop_clock( calculation  )
  CALL print_clock( calculation )

  WRITE (stdout, 1000)
  WRITE (stdout,'(5x,a)') '                           END JOB XSpectra'
  WRITE (stdout, 1001)
  CALL stop_xspectra () 

  ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  ! Formats 
  ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

 1000 FORMAT(/,5x,&
  '-------------------------------------------------------------------------')
 1001 FORMAT(5x,&
  '-------------------------------------------------------------------------',&
  /)

END program X_Spectra

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
SUBROUTINE stop_xspectra
  !----------------------------------------------------------------------------
  !
  ! Synchronize processes before stopping. This is a copy of stop_pp.
  !
  USE control_flags, ONLY: twfcollect
  USE io_files, ONLY: iunwfc
  USE mp_global, ONLY: mp_global_end
  USE parallel_include
  !
#ifdef __MPI

  INTEGER :: info
  LOGICAL :: op

  INQUIRE ( iunwfc, opened = op )

  IF ( op ) THEN
     IF (twfcollect) THEN
        CLOSE (unit = iunwfc, status = 'delete')
     ELSE
        CLOSE (unit = iunwfc, status = 'keep')
     ENDIF
  ENDIF

  CALL mp_global_end()

#endif

  STOP
END SUBROUTINE stop_xspectra


!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!  XANES K-edge calculation in the electric dipole approximation
!------------------------------------------------------------------------------
SUBROUTINE xanes_dipole(a,b,ncalcv,xnorm,core_wfn,paw_iltonhb,&
                        terminator,verbosity)
  !----------------------------------------------------------------------------
  USE constants,       ONLY : fpi
  USE io_global,       ONLY : stdout     ! Modules/io_global.f90
  USE kinds,           ONLY : DP
  USE parameters,      ONLY : ntypx
  USE radial_grids,    ONLY : ndmx
  USE ions_base,       ONLY : nat, ntyp => nsp, ityp
  USE wvfct,           ONLY : npwx, nbndx, nbnd, npw, igk, g2kin, et,&
                              !current_k, ecutwfc
                              ecutwfc
  USE symm_base,       ONLY : d1,d2,d3
  USE noncollin_module,ONLY : noncolin
  USE lsda_mod,        ONLY : nspin,lsda,isk,current_spin
  USE cell_base,       ONLY : tpiba2, bg
  USE wavefunctions_module, ONLY: evc
  USE klist,           ONLY : &
       nkstot,                & ! total number of k-points
       nks,                   & ! number of k-points per pool
       xk,                    & ! k-points coordinates
       wk                       ! k-points weight
  USE gvect,           ONLY: g, ngm, ngl
  USE fft_base,        ONLY: dfftp
  USE paw_gipaw,       ONLY : &
       paw_vkb,               & ! |p> projectors
       paw_becp,              & ! product of projectors and wf.
       paw_nkb,               & ! total number of beta functions, with st.fact.
       paw_lmaxkb,paw_recon
  USE becmod,          ONLY : becp, allocate_bec_type, deallocate_bec_type !CG
  USE scf,             ONLY : vltot, vrs, v, kedtau
  USE gvecs,           ONLY : doublegrid
  USE mp_world,        ONLY : world_comm
  USE mp_pools,        ONLY : intra_pool_comm, root_pool, npool
  USE mp,              ONLY : mp_sum, mp_bcast, mp_barrier !CG
  USE io_global,       ONLY : ionode

  USE xspectra,        ONLY : xiabs, xanes_dip, xang_mom, xniter,&
                              xnitermax, xepsilon,time_limit,calculated,&
                              save_file_kind
  USE atom,            ONLY : rgrid, msh
  !  use atom,        ONLY : &
  !       mesh,     &!mesh(ntypx) number of mesh points
  !       msh ,     &!msh(ntypx)the point at rcut=end of radial integration
  !       r   
  USE radin_mod
  USE basis,           ONLY : natomwfc
  USE uspp,            ONLY : vkb, nkb, okvan !CG
  USE uspp_param,      ONLY : upf
  USE ldaU,            ONLY : lda_plus_u, init_lda_plus_u 
  !<CG>
  USE xspectra_paw_variables, ONLY : xspectra_paw_nhm
  !</CG>

  IMPLICIT NONE
  !
  REAL(dp), INTENT(INOUT) :: a(xnitermax,1,nks)
  REAL(dp), INTENT(INOUT) :: b(xnitermax,1,nks)     
  REAL(dp), INTENT(INOUT) :: xnorm(1,nks)
  REAL(dp), INTENT(IN)    :: core_wfn(ndmx)
  INTEGER, INTENT(INOUT)  :: ncalcv(1,nks)
  INTEGER, INTENT(IN)     :: paw_iltonhb(0:paw_lmaxkb,xspectra_paw_nhm,ntyp)
  LOGICAL, INTENT(IN)     :: terminator
  !
  !... Local variables
  !
  INTEGER  :: is,ik,iabso,nr,ip,jp,l,j,icrd,ip_l,nrc,nt,na
  INTEGER  :: ipx,ipx_0,ipy,ipz,nline,nrest,npw_partial
  INTEGER  :: nunfinished
  LOGICAL  :: recalc
  REAL(dp) :: pref,prefb,v_of_0,xnorm_partial
  REAL(dp) :: norm, normps
  REAL(dp), ALLOCATABLE :: aux(:)
  COMPLEX(KIND=DP), EXTERNAL :: zdotc
  COMPLEX(dp), ALLOCATABLE :: paw_vkb_cplx(:,:)
  COMPLEX(dp), ALLOCATABLE :: psiwfc(:), spsiwfc(:)
  CHARACTER(LEN=4) :: verbosity

  REAL(dp) :: timenow 
  REAL(DP), EXTERNAL ::  get_clock
  EXTERNAL :: zdscal

  timenow=0
  nunfinished=0

  ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  !  Constant definitions
  ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  pref=SQRT(3.d0/2.d0)
  prefb=SQRT(3.d0)

  ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  !  Variable allocation and initialization
  ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  ALLOCATE(aux(rgrid(xiabs)%mesh)) !allocation too big, it needs only up to msh
  ALLOCATE (paw_vkb_cplx( npwx,  paw_nkb))
  ALLOCATE(xanes_dip(paw_recon(xiabs)%paw_nl(xang_mom)))
  ALLOCATE(psiwfc(npwx))
  xanes_dip(:)=0.d0

  ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  ! Radial Dipole Matrix element
  ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  !  radial part:  <\phi_n|r|\psi_i>=int r^3 \phi_n \psi_i dr.
  !  Here only ps and ae wavefunctions from XX.recon
  !
  WRITE(stdout,'(8x,a)') &
       'Radial transition matrix element(s) used in the calculation of the'
  WRITE(stdout,'(8x,a)') &
       'initial vector of the Lanczos basis (|tilde{phi}_abs> normalized)'
  !WRITE(stdout,'(5x,a,i2,a,i2,a,i3,a)') 'There are ',&
  !                                      paw_recon(xiabs)%paw_nl(xang_mom),&
  !                                      ' projector(s)/channel for l=',&
  !                                      xang_mom,' and atom type',&
  !                                      xiabs,'.'

  ! ... Checks that the core wf is correctly normalized
  !
  nr=msh(xiabs)  ! extended up to all the NON ZERO points in the mesh.
  IF(TRIM(verbosity).EQ.'high') THEN
     aux(1:nr)=core_wfn(1:nr)*core_wfn(1:nr)
     WRITE (stdout,'(8x,"Norm of core wfc = ",f10.6)') &
           SQRT(para_radin(aux(1:nr),rgrid(xiabs)%r(1:nr),nr))
  ENDIF

  ! ... Calculates the radial integral
  !
  ip_l=0

  DO ip=1,paw_recon(xiabs)%paw_nbeta
     !  IF(psphi(xiabs,ip)%label%l.EQ.xang_mom) THEN
     IF(paw_recon(xiabs)%aephi(ip)%label%l.EQ.xang_mom) THEN
        nrc=paw_recon(xiabs)%aephi(ip)%label%nrc
        ip_l=ip_l+1
        ! here below, recall that psi is r*psi and you have a Jacobian=r^2
        aux(1:nrc) = rgrid(xiabs)%r(1:nrc) * &
                     paw_recon(xiabs)%aephi(ip)%psi(1:nrc) * &
                     core_wfn(1:nrc)
        ! here we have to integrate only inside the augmentation region.
        xanes_dip(ip_l)=para_radin(aux(1:nrc),rgrid(xiabs)%r(1:nrc),nrc)
     ENDIF
  ENDDO

  !... Writes the radial transition matrix element(s)
  ! 
  DO ip=1,paw_recon(xiabs)%paw_nl(xang_mom) !
     WRITE(stdout,'(8x,a,i1,a,i1,a,f14.9)') &
            '| For PAW proj. (l=',xang_mom,') #',&
            ip, ': radial matrix element =', xanes_dip(ip)
  ENDDO
  WRITE(stdout,*)
  DEALLOCATE(aux)

  ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  !  Determines the index of the first projector of the absorbing atom
  ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  !
  ipx_0=0
  IF(xiabs.NE.1) THEN
     DO nt=1, xiabs-1
        DO na=1, nat
           IF (ityp(na).EQ.nt) THEN
              ipx_0=ipx_0+paw_recon(nt)%paw_nh
           ENDIF
        ENDDO
     ENDDO
  ENDIF

  ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  ! Starts the loop over the k-points
  ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  !*apsi  call set_vrs(vrs,vltot,vr,nrxx,nspin,doublegrid)
  !<CG>
  CALL set_vrs(vrs,vltot,v%of_r,kedtau, v%kin_r,dfftp%nnr,nspin,doublegrid)
  !</CG>

  IF (npool /= 1) WRITE(stdout,'(a,i5,a,i3)') 'NB: the ', nks,&
     ' k-point are not all listed below because npool=', npool

  DO ik=1,nks

     WRITE(stdout,'(8x,a)')&
           '|-------------------------------------------------------------'
     WRITE(stdout,'(8x,a ,i5,a,3(f7.4,a),f7.4,a,i3)') '! k-point # ',ik, &
        ':  (', xk(1,ik),', ',xk(2,ik),', ',xk(3,ik),'), ',wk(ik),', ', isk(ik)
     WRITE(stdout,'(8x,a)')&
           '|-------------------------------------------------------------'

     IF (calculated(1,ik).EQ.1) CYCLE

     timenow=get_clock( 'xanes' )
     CALL mp_bcast(timenow,root_pool, intra_pool_comm)     

     IF( timenow.GT.time_limit) THEN
       nunfinished=1
       EXIT
     ENDIF

     IF (verbosity.eq.'high') &
       WRITE(stdout,'(8x,"|   Total cpu time spent up to now: ",F9.2," s")')&
             timenow 

     ! <DC> current_k=ik </DC>

     IF(lsda) current_spin=isk(ik)

     !... gk_sort: sort k-points and exit kinetic energies 
     !
     CALL gk_sort(xk (1,ik),ngm,g,ecutwfc/tpiba2,npw,igk,g2kin)  !CHECK
     g2kin=g2kin*tpiba2                                          !CHECK
     
     npw_partial = npw
     CALL mp_sum( npw_partial, intra_pool_comm )

     IF(xniter.ge.npw_partial) THEN
        xniter = npw_partial
        WRITE(stdout,'(8x,a)') '|   Hilbert space is saturated'
        WRITE(stdout,'(8x,a,i10)') '|   xniter is set equal to ',npw_partial
        WRITE(stdout,'(8x,a)') &
               '|   Increase kinetic-energy cutoff in your SCF calculation!'
     ENDIF

     !<CG>        
     CALL init_gipaw_2(npw,igk,xk(1,ik),paw_vkb)
     !</CG>
     IF (.NOT.lda_plus_u) CALL init_us_2(npw,igk,xk(1,ik),vkb)
     IF (lda_plus_u) CALL orthoUwfc_k(ik)

     ! Angular Matrix element
     !
     !... Calculates the complex PAW projectors, paw_vkb_cplx, from
     !     LC of paw_vkb, i.e., the real PAW projectors expressed using
     !     real spherical harmonics)
     !
     !*************************************************************************
     ! Here I define human projectors <CG>
     !
     ! paw_vkb(1:npw,ipx) are real spherical harmonics. 
     ! The real spherical harmonics are defined as
     ! 
     !     y_{l,2m}  =[Y_{l,m}+(-1)^m Y_{l,-m}]/SQRT(2)
     !     y_{l,2m+1}=[Y_{l,m}-(-1)^m Y_{l,-m}]/(i*SQRT(2))
     !
     ! (remember Y_{l,m}=(-1)^m Y_{l,-m)^*  )
     !
     ! The complex spherical harmonics can be written has a function of the real
     ! ones as:
     !
     !     Y_{l,m}  =        [y_{l,2m}+iy_{l,2m+1}]/SQRT(2)     
     !     Y_{l,-m} = (-1)^m [y_{l,2m}-iy_{l,2m+1}]/SQRT(2)     
     !
     ! The paw_vkb_cplx are the Y_{l,m} so the usual spherical harmonics.
     !
     ! Rotational invariance has been checked.
     !*************************************************************************

     DO ip=1,paw_recon(xiabs)%paw_nl(xang_mom)

        ipx=ipx_0+paw_iltonhb(xang_mom,ip,xiabs)
        ! m= 0
        paw_vkb_cplx(1:npw,ipx)= paw_vkb(1:npw,ipx)
        ! m=+1 
        paw_vkb_cplx(1:npw,ipx+1)=       &
             (paw_vkb(1:npw,ipx+1)+(0.d0,1.d0)*paw_vkb(1:npw,ipx+2))/SQRT(2.0)
        ! m=-1 
        paw_vkb_cplx(1:npw,ipx+2)=-       &
             (paw_vkb(1:npw,ipx+1)-(0.d0,1.d0)*paw_vkb(1:npw,ipx+2))/SQRT(2.0)

     ENDDO



     !... Calculates the initial state for the Lanczos procedure,
     !    stored in psiwfc.
     !    It includes the radial matrix element, the angular matrix element,
     !    and the associated PAW projector.
     !
     psiwfc(1:npw)=(0.d0,0.d0)
     DO ip=1,paw_recon(xiabs)%paw_nl(xang_mom)   
        ipx=ipx_0+paw_iltonhb(xang_mom,ip,xiabs)

        !      WARNING, storage of spherical harmonics is the following:
        !         given Y_{l,m}=P_{lm}exp(im\phi) one has
        !                                                         counter
        !   l, m=0  -------> Y_{l,0}                              1     z
        !   l, m=+1 -------> Y_{l,1}+Y_{l,-1} or cos(m\phi)       2     
        !   l, m=-1 -------> Y_{l,1}-Y_{l,-1} or sin(m\phi)       3     
        !  .....
        !   l, m=+l -------> Y_{l,1}+Y_{l,-1} or cos(m\phi)      2*l
        !   l, m=-l -------> Y_{l,1}-Y_{l,-1} or sin(m\phi)      2*l+1


        psiwfc(1:npw)=psiwfc(1:npw) + ( pref                                   &
                                        * paw_vkb_cplx(1:npw,ipx+2)            &
                                        * (xepsilon(1)+(0.d0,1.d0)*xepsilon(2))&
                                       +pref                                   &
                                        * paw_vkb_cplx(1:npw,ipx+1)            &
                                        *(-xepsilon(1)+(0.d0,1.d0)*xepsilon(2))&
                                       +prefb                                  &
                                        *paw_vkb_cplx(1:npw,ipx)               &
                                        *xepsilon(3)                           &
                                      )                                        &
                                    * xanes_dip(ip)/SQRT(fpi)

     ENDDO
     psiwfc(1:npw)=psiwfc(1:npw)*SQRT(fpi)/3.0

     !... Normalizes the wavefunction psiwfc(1:npw)
     !

     !<CG>
     CALL allocate_bec_type(nkb,1,becp)

     write(6,*) 'okvan=',okvan
     IF (okvan) THEN
        ALLOCATE(spsiwfc(npwx))
        spsiwfc(:)=(0.d0,0.d0)
        recalc=.true.
        CALL sm1_psi(recalc,npwx, npw, 1, psiwfc, spsiwfc)
        xnorm_partial=zdotc(npw,psiwfc,1,spsiwfc,1)
        DEALLOCATE(spsiwfc)
     ELSE
!        xnorm_partial=0.d0
!        do ip=1,npw
!          xnorm_partial=xnorm_partial+conjg(psiwfc(ip))*psiwfc(ip)
!       enddo
        xnorm_partial=real(zdotc(npw,psiwfc,1,psiwfc,1),dp)

     ENDIF
     !</CG>

     CALL mp_sum( xnorm_partial, intra_pool_comm )

     xnorm(1,ik)=SQRT(xnorm_partial)
     WRITE(stdout,'(8x,a,e15.8)') '|   Norm of the initial Lanczos vector:',&
                                     xnorm(1,ik)
     norm=1.d0/xnorm(1,ik)
     
     CALL zdscal(npw,norm,psiwfc,1)
    
     !... Starts the Lanczos procedure
     !
     IF (okvan) THEN
        CALL lanczos_uspp(a(:,1,ik),b(:,1,ik),psiwfc,ncalcv(1,ik),terminator)
     ELSE
        CALL lanczos(a(:,1,ik),b(:,1,ik),psiwfc,ncalcv(1,ik),terminator)
     ENDIF

     !!      Then I write small report of the lanczos results
     !!
     !IF(TRIM(verbosity).EQ.'high') THEN
     !   WRITE( stdout,*) '-----------------------------------------'
     !   WRITE( stdout,*) 'k-point number =',ik
     !   WRITE( stdout,*) 'k-point coordinate, isk'
     !   WRITE( stdout,'(3f12.6,1i2)') xk(1,ik),xk(2,ik),xk(3,ik),isk(ik)
     !   WRITE( stdout,*) 'Norm of the initial vector =',xnorm(1,ik)
     !   WRITE( stdout,*) 'Number of iterations =',ncalcv(1,ik)
     !   !        nline=ncalcv(icrd,ik)/6
     !   !        nrest=ncalcv(icrd,ik)-nline*6
     !   !        WRITE( stdout,*) 'a vectors:'
     !   !        DO ip=1,nline
     !   !           WRITE( stdout,"(6(f10.6,3x))") (a((ip-1)*6+j,icrd,ik),j=1,6)
     !   !        ENDDO
     !   !        WRITE( stdout,"(6(f10.6,3x))") (a(nline*6+j,icrd,ik),j=1,nrest)
     !   !        WRITE( stdout,*) 'b vectors:'
     !   !        DO ip=1,nline
     !   !           WRITE( stdout,"(6(f10.6,3x))") (b((ip-1)*6+j,icrd,ik),j=1,6)
     !   !        ENDDO
     !   !        WRITE( stdout,"(6(f10.6,3x))") (b(nline*6+j,icrd,ik),j=1,nrest)
     !   WRITE( stdout,*) '-----------------------------------------'
     !ENDIF

     CALL deallocate_bec_type ( becp ) ! CG
     calculated(1,ik)=1  

  ENDDO  ! LOOP on k-points
 
  CALL mp_barrier(world_comm)
  CALL mp_sum(nunfinished, world_comm)
  IF (nunfinished >= 1) THEN 
    save_file_kind='unfinished'
    WRITE(stdout,'(5x,a)') 'Calculation not finished'
  ENDIF

  DEALLOCATE(psiwfc)
  DEALLOCATE(xanes_dip)
  DEALLOCATE (paw_vkb_cplx)

END SUBROUTINE xanes_dipole

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!  XANES calculation in the electric quadrupole approximation
!------------------------------------------------------------------------------
SUBROUTINE xanes_quadrupole(a,b,ncalcv,xnorm,core_wfn,paw_iltonhb,&
                            terminator,verbosity)
  !----------------------------------------------------------------------------
  USE io_global,       ONLY: stdout     ! Modules/io_global.f90
  USE kinds,           ONLY: DP
  USE constants,       ONLY: pi
  USE parameters,      ONLY: ntypx
  USE radial_grids,    ONLY: ndmx
  USE ions_base,       ONLY: nat, ntyp => nsp, ityp
  USE wvfct,           ONLY: npwx,nbndx,nbnd,npw,igk,&
       !g2kin,et, current_k, ecutwfc
                              g2kin,et, ecutwfc
  !       ,igk_l2g
  USE lsda_mod,        ONLY: nspin,lsda,isk,current_spin
  USE cell_base,       ONLY: tpiba2, bg
  USE wavefunctions_module, ONLY: evc
  USE klist,           ONLY: &
       nkstot,               & ! total number of k-points
       nks,                  & ! number of k-points per pool
       xk,                   & ! k-points coordinates
       wk                      ! k-points weight
  USE gvect,           ONLY: g,ngm,ngl
  USE fft_base,        ONLY: dfftp
  USE paw_gipaw,       ONLY: &
       paw_vkb,              & ! |p> projectors
       paw_becp,             & ! product of projectors and wf.
       paw_nkb,              & ! total number of beta functions, with st.fact.
       paw_lmaxkb,           & 
       paw_recon
  USE becmod,          ONLY: becp, allocate_bec_type, deallocate_bec_type ! CG
  USE scf,             ONLY: vltot, v, vrs, kedtau !CG
  USE gvecs,           ONLY: doublegrid
  USE mp_pools,        ONLY: intra_pool_comm, root_pool, npool
  USE mp_world,        ONLY: world_comm
  USE mp,              ONLY: mp_sum,mp_barrier, mp_bcast !CG
  USE xspectra,        ONLY: xiabs, xanes_qua, xang_mom, xniter, xnitermax,&
                             xkvec, xepsilon, save_file_kind,              &
                             calculated, time_limit
  USE atom,            ONLY: rgrid, msh
  !  use atom,        ONLY : &
  !       mesh,     &!mesh(ntypx) number of mesh points
  !       msh ,     &!msh(ntypx)the point at rcut=end of radial integration
  !       r   
  USE radin_mod
  USE uspp,            ONLY: vkb, nkb, okvan !CG
  USE ldaU,            ONLY: lda_plus_u
  USE basis,           ONLY: natomwfc
  !<CG>
  USE xspectra_paw_variables, ONLY: xspectra_paw_nhm
  !</CG>

  IMPLICIT NONE
  !
  REAL(dp), INTENT(INOUT) :: a(xnitermax,1,nks)
  REAL(dp), INTENT(INOUT) :: b(xnitermax,1,nks)     
  REAL(dp), INTENT(INOUT) :: xnorm(1,nks)
  REAL(dp), INTENT(IN)    :: core_wfn(ndmx)
  INTEGER, INTENT(INOUT)  :: ncalcv(1,nks)
  INTEGER, INTENT(IN)     :: paw_iltonhb(0:paw_lmaxkb,xspectra_paw_nhm, ntyp)
  LOGICAL, INTENT(IN)     :: terminator
  !
  !... Local variables
  !
  INTEGER :: is,ik,iabso,nr,ip,jp,l,j,icrd,ip_l,nrc,nt,na
  INTEGER :: ipx,ipx_0,ipy,ipz,nline,nrest,npw_partial
  INTEGER :: nunfinished
  LOGICAL :: recalc
  REAL (dp) :: pref,prefb,v_of_0,xnorm_partial,prefm2,prefm1,prefm0
  REAL (dp) :: norm,normps


  REAL (dp), ALLOCATABLE :: aux(:)
  COMPLEX(KIND=DP), EXTERNAL :: zdotc
  COMPLEX(dp), ALLOCATABLE :: paw_vkb_cplx(:,:)
  COMPLEX(KIND=dp), ALLOCATABLE :: psi(:)
  COMPLEX(dp), ALLOCATABLE :: psiwfc(:), spsiwfc(:)
  CHARACTER(LEN=4) :: verbosity

  REAL(DP), EXTERNAL ::  get_clock
  REAL(dp) :: timenow
  EXTERNAL zdscal

  nunfinished=0
  timenow=0

  ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  !  Constant Definitions
  ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$


  prefm2=SQRT(3.0/40.0)/3.0
  prefm1=prefm2
  prefm0=2.0*SQRT(1.0/40.0)/3.0

  pref=SQRT(2.d0)
  prefb=1.0/pref


  ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  !  Variable allocation and initialization
  ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ 

  ALLOCATE(aux(rgrid(xiabs)%mesh))!overdimensionated, necessary only up to msh
  ALLOCATE(psi(npwx))
  ALLOCATE (paw_vkb_cplx( npwx,  paw_nkb))
  ALLOCATE(xanes_qua(paw_recon(xiabs)%paw_nl(xang_mom)))
  ALLOCATE(psiwfc(npwx))
  xanes_qua(:)=0.d0
  psi(:)=(0.d0)

  ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  ! Radial Quadrupole Matrix element
  ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  !  Radial part :   <\phi_n|r|\psi_i>=int r^3 \phi_n \psi_i dr.
  !  Here only ps and ae wavefunctions from XX.recon
  !

  

  !WRITE( stdout,*) 'Calculation Quadrupole matrix element'
  !WRITE( stdout,*) 'There are ',paw_recon(xiabs)%paw_nl(xang_mom),'
  !projectors/channel'
  !WRITE( stdout,*) 'xang_mom=',xang_mom,' xiabs=',xiabs

  ! ... Checks that the core wf is correctly normalized

  nr=msh(xiabs)  ! extended up to all the NON ZERO points in the mesh.
  IF(TRIM(verbosity).EQ.'high') THEN
     aux(1:nr)=core_wfn(1:nr)*core_wfn(1:nr)
     WRITE (stdout,'(8x,"Norm of core wfc =",f10.6)') &
            SQRT(para_radin(aux(1:nr),rgrid(xiabs)%r(1:nr),nr))
  ENDIF

  ! ... Calculate the radial integral

  ip_l=0

  DO ip=1,paw_recon(xiabs)%paw_nbeta
     IF(paw_recon(xiabs)%psphi(ip)%label%l.EQ.xang_mom) THEN
        nrc=paw_recon(xiabs)%aephi(ip)%label%nrc
        ip_l=ip_l+1
        !    here below, recall that psi is r*psi and you have a Jacobian=r^2
        !        aux(1:nr)=r(1:nr,xiabs)*r(1:nr,xiabs)* &
        !             psphi(xiabs,ip)%psi(1:nr)*core_wfn(1:nr)
        aux(1:nrc)=rgrid(xiabs)%r(1:nrc)*rgrid(xiabs)%r(1:nrc)*&
                   paw_recon(xiabs)%aephi(ip)%psi(1:nrc)*core_wfn(1:nrc)
        !    here we have to integrate only inside the augmentation region.
        xanes_qua(ip_l)=para_radin(aux(1:nrc),rgrid(xiabs)%r(1:nrc),nrc)
     ENDIF
  ENDDO
 
  ! ... Write the radial transition matrix element(s)


  DO ip=1,paw_recon(xiabs)%paw_nl(xang_mom) !
     WRITE(stdout,'(8x,a,i1,a,i1,a,f14.9)') &
            '| For PAW proj. (l=',xang_mom,') #',&
            ip, ': radial matrix element =', xanes_qua(ip)
  ENDDO
  WRITE(stdout,*)
  DEALLOCATE(aux)
  DEALLOCATE(psi)

  ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  !  Determines the index of the first projector of the absorbing atom
  ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$  

  ipx_0=0
  IF(xiabs.NE.1) THEN
     DO nt=1, xiabs-1
        DO na=1, nat
           IF (ityp(na).EQ.nt) THEN
              ipx_0=ipx_0+paw_recon(nt)%paw_nh
           ENDIF
        ENDDO
     ENDDO
  ENDIF


  ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  ! Starts the loop over the k-points
  ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$  


  !*apsi  call set_vrs(vrs,vltot,vr,nrxx,nspin,doublegrid)
  !  CALL set_vrs(vrs,vltot,v%of_r,nrxx,nspin,doublegrid)
  !  CALL newd

  !<CG>
  ! set_vrs (vrs, vltot, vr, kedtau, kedtaur,nrxx, nspin, doublegrid)
  CALL set_vrs(vrs,vltot,v%of_r,kedtau, v%kin_r,dfftp%nnr,nspin,doublegrid)
  !</CG>


  DO ik=1,nks

     WRITE(stdout,'(8x,a)')&
           '|-------------------------------------------------------------'
     WRITE(stdout,'(8x,a ,i5,a,3(f7.4,a))') '! k-point # ',ik, &
        ':  (', xk(1,ik),', ',xk(2,ik),', ',xk(3,ik),') '
     WRITE(stdout,'(8x,a ,f7.4,a,i3)') '! weight:',wk(ik), &
        ' spin state:', isk(ik)
     WRITE(stdout,'(8x,a)')&
           '|-------------------------------------------------------------'
     
     IF (calculated(1,ik).EQ.1) CYCLE

     timenow=get_clock( 'xanes' )
     CALL mp_bcast(timenow,root_pool, intra_pool_comm)

     IF( timenow.GT.time_limit) THEN
       nunfinished=1
       EXIT
     ENDIF

     IF (verbosity.eq.'high') &
     WRITE(stdout,'(8x, "| Total cpu time spent up to now is ",F9.2," s")')&
           timenow

     !<DC> current_k=ik </DC>
     
     IF(lsda) current_spin=isk(ik)

     !... gk_sort : sort k-points and exit kinetic energies 
     CALL gk_sort(xk (1,ik),ngm,g,ecutwfc/tpiba2,npw,igk,g2kin)  !CHECK
     g2kin=g2kin*tpiba2                                          !CHECK

     npw_partial = npw
     CALL mp_sum( npw_partial, intra_pool_comm )

     IF(xniter.ge.npw_partial) THEN
        xniter = npw_partial

        WRITE(stdout,'(8x,a)') '|   Hilbert space is saturated'
        WRITE(stdout,'(8x,a,i10)') '|   xniter is set equal to ',npw_partial
        WRITE(stdout,'(8x,a)') &
               '|   Increase kinetic-energy cutoff in your SCF calculation!'
        !        CALL stop_pp
     ENDIF

     !<CG>
     CALL init_gipaw_2(npw,igk,xk(1,ik),paw_vkb)
     !</CG>
     if(.not.lda_plus_u) CALL init_us_2(npw,igk,xk(1,ik),vkb)
     IF (lda_plus_u) CALL orthoUwfc_k(ik)

     ! Angular Matrix element
     !
     !... Calculates the complex PAW projectors, paw_vkb_cplx, from
     !     LC of paw_vkb, i.e., the real PAW projectors expressed using
     !     real spherical harmonics)
     !
     !*************************************************************************
     ! Here I define human projectors
     !
     ! paw_vkb(1:npw,ipx) are real spherical harmonics. The real spherical 
     ! harmonics are defined as
     ! 
     !     y_{l,2m}  =[(-)^{m} Y_{l,m}+Y_{l,-m}]/SQRT(2)
     !     y_{l,2m+1}=[(-)^{m} Y_{l,m}-Y_{l,-m}]/(i*SQRT(2))
     !
     ! The complex spherical harmonics can be written has a function of the real
     ! ones as :
     !
     !     Y_{l,m}  =       [y_{l,2m}+iy_{l,2m+1}]/SQRT(2)     
     !     Y_{l,-m} = (-)^m [y_{l,2m}-iy_{l,2m+1}]/SQRT(2)     
     !
     !  The paw_vkb_cplx are the Y_{l,m}, so the usual spherical harmonics
     !*************************************************************************



     DO ip=1,paw_recon(xiabs)%paw_nl(xang_mom)   

        ipx=ipx_0+paw_iltonhb(xang_mom,ip,xiabs)
        paw_vkb_cplx(1:npw,ipx)= paw_vkb(1:npw,ipx)   !m=0
        paw_vkb_cplx(1:npw,ipx+1)=       &
             (paw_vkb(1:npw,ipx+1)+(0.d0,1.d0)*paw_vkb(1:npw,ipx+2))/SQRT(2.0)  
             !m=+1
        paw_vkb_cplx(1:npw,ipx+2)=       &
             -(paw_vkb(1:npw,ipx+1)-(0.d0,1.d0)*paw_vkb(1:npw,ipx+2))/SQRT(2.0)
             !m=-1
        paw_vkb_cplx(1:npw,ipx+3)=       &
             (paw_vkb(1:npw,ipx+3)+(0.d0,1.d0)*paw_vkb(1:npw,ipx+4))/SQRT(2.0)
             !m=+2
        paw_vkb_cplx(1:npw,ipx+4)=       &
             (paw_vkb(1:npw,ipx+3)-(0.d0,1.d0)*paw_vkb(1:npw,ipx+4))/SQRT(2.0)
             !m=-2
     ENDDO


     !... Calculates the initial state for the Lanczos procedure,
     !    stored in psiwfc.
     !    It includes the radial matrix element, the angular matrix element,
     !    and the associated PAW projector.
     ! 

     psiwfc(:)=(0.d0,0.d0)

     DO ip=1,paw_recon(xiabs)%paw_nl(xang_mom)   
        ipx=ipx_0+paw_iltonhb(xang_mom,ip,xiabs)

        !      WARNING, storage of spherical harmonics is the following:
        !         given Y_{l,m}=P_{lm}exp(im\phi) one has
        !                                                         counter
        !   l, m=0  -------> Y_{l,0}                              1     z
        !   l, m=+1 -------> Y_{l,1}+Y_{l,-1} or cos(m\phi)       2     
        !   l, m=-1 -------> Y_{l,1}-Y_{l,-1} or sin(m\phi)       3     
        !  .....
        !   l, m=+l -------> Y_{l,1}+Y_{l,-1} or cos(m\phi)      2*l
        !   l, m=-l -------> Y_{l,1}-Y_{l,-1} or sin(m\phi)      2*l+1

        psiwfc(1:npw)=psiwfc(1:npw)+prefm2*&
             (xepsilon(1)-(0.d0,1.d0)*xepsilon(2))*&
             (xkvec(1)-(0.d0,1.d0)*xkvec(2))*&
             xanes_qua(ip)*&
             paw_vkb_cplx(1:npw,ipx+3)

        psiwfc(1:npw)=psiwfc(1:npw)+prefm2*&
             (xepsilon(1)+(0.d0,1.d0)*xepsilon(2))*&
             (xkvec(1)+(0.d0,1.d0)*xkvec(2))*&
             xanes_qua(ip)*&
             paw_vkb_cplx(1:npw,ipx+4)

        psiwfc(1:npw)=psiwfc(1:npw)-prefm1*( &
             (xepsilon(1)-(0.d0,1.d0)*xepsilon(2))* &
             xkvec(3)+ &
             (xkvec(1)-(0.d0,1.d0)*xkvec(2))* &
             xepsilon(3) )* &
             xanes_qua(ip)*&
             paw_vkb_cplx(1:npw,ipx+1)

        psiwfc(1:npw)=psiwfc(1:npw)+prefm1*( &
             (xepsilon(1)+(0.d0,1.d0)*xepsilon(2))* &
             xkvec(3)+ &
             (xkvec(1)+(0.d0,1.d0)*xkvec(2))* &
             xepsilon(3) )* &
             xanes_qua(ip)*&
             paw_vkb_cplx(1:npw,ipx+2)

        psiwfc(1:npw)=psiwfc(1:npw)+prefm0*( &
             pref*xkvec(3)*xepsilon(3)-&
             prefb*(xepsilon(1)*xkvec(1)+xepsilon(2)*xkvec(2)))*&
             xanes_qua(ip)*&
             paw_vkb_cplx(1:npw,ipx)

     ENDDO


      !... Normalizes the wavefunction psiwfc(1:npw)

!     CALL allocate_bec_type(nkb,1, becp) ! CG
     CALL allocate_bec_type(nkb,natomwfc, becp) 
     !<CG>
     IF (okvan) THEN
        ALLOCATE(spsiwfc(npwx))
        spsiwfc(:)=(0.d0,0.d0)
        recalc=.true.
        CALL sm1_psi(recalc,npwx, npw, 1, psiwfc, spsiwfc)
        xnorm_partial=zdotc(npw,psiwfc,1,spsiwfc,1)
        DEALLOCATE(spsiwfc)
     ELSE
        xnorm_partial=zdotc(npw,psiwfc,1,psiwfc,1)
     ENDIF
     !</CG>

     CALL mp_sum( xnorm_partial, intra_pool_comm )

     xnorm(1,ik)=SQRT(xnorm_partial)
     WRITE(stdout,'(8x,a,e15.8)') '|   Norm of the initial Lanczos vector:',&
                                     xnorm(1,ik)
     norm=1.d0/xnorm(1,ik)

     CALL zdscal(npw,norm,psiwfc,1)
    
     !... Starts the Lanczos procedure

     IF (okvan) THEN
        CALL lanczos_uspp(a(:,1,ik),b(:,1,ik),psiwfc,ncalcv(1,ik), terminator)
     ELSE
        CALL lanczos(a(:,1,ik),b(:,1,ik),psiwfc,ncalcv(1,ik), terminator)
     ENDIF
     
     !!      Then I write small report of the lanczos results
     !IF(TRIM(verbosity).EQ.'high') THEN
      ! WRITE( stdout,*) '-----------------------------------------'
      ! WRITE( stdout,*) 'k-point number =',ik
      ! WRITE( stdout,*) 'Norm of the initial vector =',xnorm(1,ik)
      ! WRITE( stdout,*) 'Number of iterations =',ncalcv(1,ik)
      ! WRITE( stdout,*) 'k-point coordinate, isk'
      ! WRITE( stdout,'(3f12.6,1i2)') xk(1,ik),xk(2,ik),xk(3,ik),isk(ik)
        !     nline=ncalcv(1,ik)/6
        !     nrest=ncalcv(1,ik)-nline*6
        !     WRITE( stdout,*) 'a vectors:'
        !     DO ip=1,nline
        !        WRITE( stdout,"(6(f10.6,3x))") (a((ip-1)*6+j,1,ik),j=1,6)
        !     ENDDO
        !     WRITE( stdout,"(6(f10.6,3x))") (a(nline*6+j,1,ik),j=1,nrest)
        !     WRITE( stdout,*) 'b vectors:'
        !     DO ip=1,nline
        !        WRITE( stdout,"(6(f10.6,3x))") (b((ip-1)*6+j,1,ik),j=1,6)
        !     ENDDO
        !     WRITE( stdout,"(6(f10.6,3x))") (b(nline*6+j,1,ik),j=1,nrest)
      ! WRITE( stdout,*) '-----------------------------------------'
    !ENDIF

     CALL deallocate_bec_type (becp) ! CG
     calculated(1,ik)=1

  ENDDO   !LOOP on k-points

  CALL mp_barrier(world_comm)
  CALL mp_sum(nunfinished, world_comm)
  IF (nunfinished >= 1) THEN
    save_file_kind='unfinished'
    write(stdout,'(5x,a)') 'Calculation not finished'
  ENDIF

  DEALLOCATE(psiwfc)
  DEALLOCATE (paw_vkb_cplx)
  DEALLOCATE(xanes_qua)

END SUBROUTINE xanes_quadrupole


!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
SUBROUTINE lanczos (a,b,psi,ncalcv,terminator)
  !----------------------------------------------------------------------------
  ! subroutine written by CG
  USE kinds,     ONLY: DP
  USE constants, ONLY: rytoev
  USE wvfct,     ONLY: npwx,nbndx, nbnd,npw,igk,g2kin
  USE becmod,    ONLY: becp
  USE uspp,      ONLY: vkb, nkb
  !USE cell_base, ONLY: omega
  USE xspectra,  ONLY: xniter, &
                       xnepoint, &
                       xcheck_conv,&
                       xnitermax,&
                       xemin,&
                       xemax,&
                       xgamma,&
                       xerror
  USE mp_global, ONLY: intra_pool_comm  
  USE mp,        ONLY: mp_sum
  USE io_global, ONLY: stdout

  IMPLICIT NONE
  ! 
  REAL(dp), DIMENSION (xnitermax), INTENT(INOUT) :: a, b
  COMPLEX(dp), DIMENSION (npwx),   INTENT(INOUT) :: psi
  INTEGER, INTENT(INOUT) ::  ncalcv
  LOGICAL, INTENT(IN)    :: terminator

  !... Local variables
  LOGICAL  :: converge
  LOGICAL  :: iconv
  INTEGER  :: ibnd, j, i, m
  REAL(dp) :: norm, error, xemin_ry, xemax_ry, xgamma_ry
  COMPLEX (dp) :: ac,bc
  REAL(dp), ALLOCATABLE :: comp(:)
  COMPLEX (dp), ALLOCATABLE :: hpsi(:), u(:)

  REAL (dp) :: ddot
  COMPLEX (DP) :: zdotc
  EXTERNAL :: zdotc,ddot
  EXTERNAL :: h_psi


  ALLOCATE(hpsi(npwx))
  ALLOCATE(u(npwx))
  ALLOCATE(comp(xnepoint))

  hpsi(:)=(0.d0,0.d0)
  u(:)=(0.d0,0.d0)
  a(:)=0.d0
  b(:)=0.d0

  xemax_ry=xemax/rytoev
  xemin_ry=xemin/rytoev
  xgamma_ry=xgamma/rytoev

  iconv=.false.

  ! ------------------------  1st iteration --------------------------
  !
  ! -- computes H*Psi  (|u0>=|Psi>)

  CALL h_psi( npwx, npw,1, psi, hpsi )


  ! -- computes a_(1)=<Psi|HPsi>

  a(1)=dble(zdotc(npw,psi,1,hpsi,1))

  CALL mp_sum(a(1), intra_pool_comm)

  ac=-a(1)*(1.d0,0.d0)
  !
  ! -- computes t3=H*Psi-a_1*Psi

  CALL zaxpy(npw,ac,psi,1,hpsi,1)

  !
  ! -- computes the norm of t3

  b(1) = dble(zdotc(npw,hpsi,1,hpsi,1))
  CALL mp_sum( b(1), intra_pool_comm )
  b(1) = SQRT( b(1) )
  !
  ! -- computes the vector |u1>

  CALL zdscal(npw,1.d0/b(1),hpsi,1)

  !
  ! -- saves |u1>


  u(1:npw)=hpsi(1:npw)

  hpsi(:)=(0.d0,0.d0)


  !
  ! -------------------------- Next iterations -----------------------
  !


  comp(:)=0.d0
  comp(1)=1.d0



  DO i=2,xniter


     !From here below we have:
     !   u=psi_j
     !   hpsi=empty
     !   psi=psi_j-1

     CALL h_psi( npwx, npw, 1, u, hpsi )


     !
     ! I compute hpsi=hpsi-b_{j-1}*psi_{j-1} (j is actually i-1)

     bc=-b(i-1)*(1.d0,0.d0)
     CALL zaxpy(npw,bc,psi,1,hpsi,1)

     !
     ! computes a(i)=<t2|t3>=<t2|H|t2>


     a(i)=REAL(zdotc(npw,hpsi,1,u,1),dp)
     CALL mp_sum( a(i), intra_pool_comm )
     !
     ! computes t3=t3-a(i)*t2
     !

     ac=-a(i)*(1.d0,0.d0)
     CALL zaxpy(npw,ac,u,1,hpsi,1)


     !
     ! computes b(i) the norm of t3
     !

     b(i)=REAL(zdotc(npw,hpsi,1,hpsi,1),dp)
     CALL mp_sum( b(i), intra_pool_comm )
     b(i) = SQRT( b(i) )


     !
     ! saves initial vector in t1 (t1 = t2)

     psi(1:npw)=u(1:npw)
     !
     ! computes vector t3/norm(t3) and saves it in t2

     CALL zdscal(npw,1.d0/b(i),hpsi,1)
     u(1:npw)=hpsi(1:npw)
     hpsi(1:npw)=(0.d0,0.d0)


     !
     !   I should gather all the ncalcv,error at the end and write only
     !   then.
     !
     IF(mod(i,xcheck_conv).EQ.0) THEN
        IF(converge(a,b,i,comp,error,xemin_ry,xemax_ry,&
                    xgamma_ry,xnepoint,xerror,terminator)) THEN
           WRITE(stdout,'(8x,a,i6,a,f12.8)') '!   => CONVERGED at iter ',i,&
                                             ' with error=',error
           ncalcv=i
           iconv=.true.
           EXIT
        ELSE
           WRITE(stdout,'(8x,a,i6,a,f12.8)') '|   Estimated error at iter ',i,&
                                             ': ', error
        ENDIF
     ENDIF

  ENDDO

  IF(.NOT.iconv) THEN
     WRITE(stdout,'(8x,a,i6,a)') '!   XANES not converged after', i-1,&
                                 ' iterations'
     WRITE(stdout,'(8x,a,i6,a,f12.8)') '!   Estimated final error after ',&
          i-1,'iterations: ', &
          converge(a,b,i-1,comp,error,xemin_ry,xemax_ry,&
                   xgamma_ry,xnepoint,xerror,terminator)
     ncalcv=i-1
  ENDIF

  DEALLOCATE(hpsi)
  DEALLOCATE(u)
  DEALLOCATE(comp)

END SUBROUTINE lanczos


!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
SUBROUTINE lanczos_uspp (a,b,psi,ncalcv,terminator)
  !----------------------------------------------------------------------------
  ! subroutine written by CG
  USE kinds,     ONLY: DP
  USE wvfct,     ONLY: npwx,nbndx, nbnd,npw,igk,g2kin
  USE becmod,    ONLY: becp, calbec
  USE constants, ONLY: rytoev
  USE uspp,      ONLY: vkb, nkb
  !USE cell_base, ONLY: omega
  USE xspectra,  ONLY: xniter,&
                       xnepoint,&
                       xcheck_conv,&
                       xnitermax,&
                       xemin,&
                       xemax,&
                       xgamma,&
                       xerror
  USE mp_global, ONLY: intra_pool_comm
  USE mp,        ONLY: mp_sum
  USE io_global, ONLY: stdout

  IMPLICIT NONE
  ! 
  REAL (dp), DIMENSION (xnitermax), INTENT (inout) :: a, b
  COMPLEX (dp), DIMENSION (npwx),   INTENT (inout) :: psi
  INTEGER, INTENT (inout) :: ncalcv
  LOGICAL, INTENT (in)    :: terminator

  !... local variables
  LOGICAL :: converge
  LOGICAL :: iconv
  LOGICAL :: recalc
  INTEGER :: ibnd,j,i,m
  REAL(dp) :: norm,error,xemin_ry,xemax_ry,xgamma_ry
  COMPLEX(dp):: ac,bc
  REAL(dp), ALLOCATABLE :: comp(:)
  COMPLEX(dp), ALLOCATABLE :: u(:), v1(:), v2(:), v3(:)
  COMPLEX(dp) :: vecteuraux1(npwx,1), vecteuraux2(npwx,1)

  REAL (dp) :: ddot
  COMPLEX(dp) :: zdotc
  EXTERNAL :: zdotc,ddot
  EXTERNAL :: h_psi

  ALLOCATE(v1(npwx))
  ALLOCATE(u(npwx))
  ALLOCATE(comp(xnepoint))
  ALLOCATE(v2(npwx))
  ALLOCATE(v3(npwx))

  v1(:) = (0.d0,0.d0)
  v2(:) = (0.d0,0.d0)
  u(:)  = (0.d0,0.d0)
  a(:)  = 0.d0
  b(:)  = 0.d0

  xemax_ry = xemax/rytoev
  xemin_ry = xemin/rytoev
  xgamma_ry= xgamma/rytoev

  iconv=.false.

  ! ------------------------  1st iteration --------------------------
  !
  recalc=.true.
  CALL sm1_psi(recalc,npwx, npw, 1, psi, v1) ! -- computes v1= S^-1 psi
  CALL h_psi( npwx, npw,1, v1, u )           ! -- computes u = H v1

  ! -- computes a_(1)=<v1|u>

  a(1)=dble(zdotc(npw,v1,1,u,1))
  CALL mp_sum(a(1), intra_pool_comm)
  ac=-a(1)*(1.d0,0.d0)
  !
  ! -- computes u= S^-1 (u-a_1*psi)

  CALL zaxpy(npw,ac,psi,1,u,1)               ! -- computes u=u-a_1*psi
  recalc=.false.
  CALL sm1_psi(recalc,npwx, npw, 1, u,v1 )  ! -- computes u=S^-1 u
  !
  ! -- computes the norm

  b(1) = zdotc(npw,u,1,v1,1)                 ! -- computes b_1 =sqrt(<u|v1>)
  CALL mp_sum( b(1), intra_pool_comm )
  b(1) = SQRT( b(1) )
  !
  ! -- computes the vector |u1>
  CALL zdscal(npw,1.d0/b(1),u,1)             ! -- computes u=u/b1
  CALL zdscal(npw,1.d0/b(1),v1,1)            ! -- computes v1=v1/b1
  !
  !
  ! -------------------------- Next iterations -----------------------
  !
  comp(:)=0.d0
  comp(1)=1.d0

  DO i=2,xniter

     CALL h_psi( npwx, npw,1, v1,v2)         ! -- computes v2= H v1

     a(i)=REAL(zdotc(npw,v1,1,v2,1),dp)      ! -- computes a_i=<v1|v2>
     CALL mp_sum( a(i), intra_pool_comm )
     !
     ! I compute hpsi=hpsi-b_{j-1}*psi_{j-1} (j is actually i-1)

     bc=-b(i-1)*(1.d0,0.d0)
     CALL zaxpy(npw,bc,psi,1,v2,1)           ! -- computes v2=v2-b_{i-1}*psi

     ac=-a(i)*(1.d0,0.d0)
     CALL zaxpy(npw,ac,u,1,v2,1)             ! -- computes v2=v2-a_i*u

     v1(:)=(0.d0,0.d0)
     recalc=.false.
     CALL sm1_psi(recalc,npwx, npw, 1,v2 ,v1 ) ! computes v1= S^-1 v2

     b(i)=REAL(zdotc(npw,v2,1,v1,1),dp)      ! -- computes b_i=sqrt(<v2|v1>)
     CALL mp_sum( b(i), intra_pool_comm )
     b(i) = SQRT( b(i) )

     !
     ! saves initial vector in t1 (t1 = t2)

     psi(1:npwx)=u(1:npwx)                   ! Psi -> u
     !
     CALL zdscal(npw,1.d0/b(i),v2,1)         ! -- computes v2=v2/b_i
     CALL zdscal(npw,1.d0/b(i),v1,1)         ! -- computes v1=v1/b_i
     u(1:npwx)=v2(1:npwx)                    ! v2 -> u
     !
     !   I should gather all the ncalcv,error at the end and write only
     !   then.
     !
     IF(mod(i,xcheck_conv).EQ.0) THEN
        IF(converge(a, b, i, comp,error, xemin_ry, xemax_ry,&
                    xgamma_ry, xnepoint, xerror, terminator)) THEN
           WRITE(stdout,'(8x,a,i6,a,f12.8)') '!   => CONVERGED at iter ',i,&
                                             ' with error=',error
           ncalcv=i
           iconv=.true.
           EXIT
        ELSE
           WRITE(stdout,'(8x,a,i6,a,f12.8)') '|   Estimated error at iter ',i,&
                                             ': ', error
        ENDIF
     ENDIF

  ENDDO

  IF(.NOT.iconv) THEN
     WRITE(stdout,'(8x,a,i6,a)') '!   XANES not converged after', i-1,&
                                 ' iterations'
     WRITE(stdout,'(8x,a,i6,a,l1)') '!   Estimated final error after ',&
          i-1,'iterations: ', &
          converge(a,b,i-1,comp,error,xemin_ry,xemax_ry,&
                   xgamma_ry,xnepoint,xerror,terminator)
     ncalcv=i-1
  ENDIF

  DEALLOCATE(v1,v2,u)
  DEALLOCATE(comp)

END SUBROUTINE lanczos_uspp

!</CG>

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
LOGICAL FUNCTION converge(a,b,m,comp,estimated_error,xemin,xemax,xgamma,&
                          xnepoint,xerror,use_term)
  !----------------------------------------------------------------------------
  USE kinds, ONLY : dp
  USE xspectra, ONLY : xnitermax
  IMPLICIT NONE

  !
  ! INPUT:
  ! -----
  !
  INTEGER :: xnepoint
  REAL(dp) :: xemin,xemax,xgamma,xerror

  !
  ! INPUT / OUTPUT:
  ! ----------------
  !
  REAL(dp) :: a(xnitermax),b(xnitermax)
  INTEGER :: m        ! number of calculated vectors
  REAL(dp) :: comp(xnepoint)
  REAL(dp) :: estimated_error  
  !
  ! ---------------------- local variables ------------------------------
  !
  REAL(dp) :: deltae,tmp,area,e,err
  REAL(dp)  :: continued_fraction
  INTEGER  :: n
  LOGICAL :: use_term
  !
  deltae = (xemax-xemin) / xnepoint
  err = 0.d0
  area = 0.d0
  n = 0
  e = xemin
  DO n = 1, xnepoint
     e = e + deltae
     tmp = continued_fraction(a,b,e,xgamma,m,use_term)
     err = err + abs(comp(n)-tmp)
     area = area + abs(tmp)
     comp(n) = tmp
  ENDDO
  err = err / area
  estimated_error = err
  IF(err < xerror) THEN
     converge = .true.
  ELSE
     converge = .false.
  ENDIF
END FUNCTION converge
 
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
FUNCTION continued_fraction(a,b,e,gamma,m, term)
  !----------------------------------------------------------------------------
  USE kinds,    ONLY: dp
  USE xspectra, ONLY: xnitermax, xcheck_conv
  USE io_global,ONLY: stdout
  IMPLICIT NONE

  ! computes the continued fraction.
  !
  ! INPUT:
  ! -----
  !
  REAL(dp) :: continued_fraction
  INTEGER  :: m
  REAL(dp) :: a(xnitermax)
  REAL(dp) :: b(xnitermax)
  REAL(dp) :: gamma
  REAL(dp) :: e
  LOGICAL :: term
  !
  ! ---------------------- local variables ------------------------------
  !
  INTEGER :: i, p,q
  COMPLEX(dp) :: res ,lastterm ! intermediate variable
  REAL(dp) :: aa, bb

  q=xcheck_conv/2
  IF (term) THEN
     aa=0.0
     bb=0.0
     DO p=1, q
        aa=aa+a(m-p)
        bb=bb+b(m-p)
     ENDDO
     aa=aa/q
     bb=bb/q
     res=lastterm(aa-e,bb*bb,gamma)
  ELSE
     res = CMPLX(a(m)-e,gamma,kind=DP)
  ENDIF
  DO i = 1, m -1
     res = CMPLX(a(m-i)-e, -gamma,kind=DP) -b(m-i)*b(m-i)/res
  ENDDO
  continued_fraction = AIMAG(1/res)
END FUNCTION continued_fraction

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
SUBROUTINE read_core_abs(filename,core_wfn,nl_init)
  !--------------------------------------------------------------------------
  USE kinds,   ONLY: DP
  USE atom,    ONLY: rgrid
  !USE atom,        ONLY : mesh     !mesh(ntypx) number of mesh points
  USE xspectra,ONLY: xiabs
  USE io_global,       ONLY : ionode, stdout

  IMPLICIT NONE

  INTEGER :: i, ierr, nbp, iblind
  INTEGER, dimension(2) :: nl_init
  CHARACTER (LEN=80) :: filename
  REAL(KIND=dp):: x
  REAL(KIND=dp):: core_wfn(*)

  !WRITE(6,*) 'xmesh=',rgrid(xiabs)%mesh
  open(unit=33,file=filename,form='formatted',iostat=ierr,status='old',err=123)

123   if( ierr == 29 )  then
          if( ionode ) write(stdout, &
          '("ERROR: core wavefunction file ",A,&
         & " does not exist or is not located in the right folder !")') &
           trim(filename)
          call stop_xspectra
        else if( ierr /= 0 ) then
          if( ionode ) write(stdout, '("ERROR reading the core wavefunction file")')
          call stop_xspectra
        else
          continue
        end if

  rewind(33)

  !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  !<OB>           Brute force identification of the core state:
  !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  nbp = rgrid(xiabs)%mesh
! Determine how many lines to skip before reading the good core function
  select case( nl_init(1) )
    case(1)
      iblind = 1                       !1s
    case(2)
      if( nl_init(2) == 0 ) then       !2s
        iblind = nbp + 2
      else                             !2p
        iblind = 2*nbp + 3
      end if
    case(3)
      if( nl_init(2) == 0 ) then       !3s
        iblind = 3*nbp + 4
      else if( nl_init(2) == 1 ) then  !3p
        iblind = 4*nbp + 5
      else
        iblind = 5*nbp + 6             !3d
      end if
  end select

  do i = 1, iblind
    READ(33,*) 
  end do

  DO i=1,nbp
   READ(33 ,*) x,core_wfn(i)
   write(277,*) x,core_wfn(i)
  ENDDO

  close(33)

END SUBROUTINE read_core_abs

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
SUBROUTINE plot_xanes_dipole(a,b,xnorm,ncalcv,terminator,e1s_ry,ispectra)
  !--------------------------------------------------------------------------
  ! Calculates and plots the electric-dipole absorption K-edge cross section 
  ! as a continued fraction,
  ! from the a_i and b_i coefficients calculated for each k-point 
  !--------------------------------------------------------------------------
  USE kinds,      ONLY: DP
  USE constants,  ONLY: rytoev, fpi
  USE xspectra,   ONLY: xang_mom, xemax, xemin, xiabs, xnepoint,n_lanczos, &
                        xgamma, xonly_plot, xnitermax, xe0_ry,edge, two_edges
  !*apsi  USE uspp_param, ONLY : psd  !psd(ntypx) label for the atoms 
  USE klist,      ONLY: nkstot, & ! total number of k-points
                        nks,    & ! number of k-points per pool
                        xk,     & ! k-points coordinates
                        wk        ! k-points weight
  !USE ener,       ONLY: ef
  USE io_global,  ONLY: stdout, ionode  
  USE mp_global,  ONLY: inter_pool_comm !CG
  USE lsda_mod,   ONLY: nspin,isk
  USE mp,         ONLY: mp_sum
  USE uspp_param, ONLY: upf
  USE gamma_variable_mod, ONLY: gamma_tab, gamma_mode, &
                                gamma_value, gamma_energy, gamma_file
  USE cut_valence_green,  ONLY: cut_occ_states, cut_desmooth, &
                                memu, meml, cut_nmemu, cut_nmeml

  IMPLICIT NONE

  REAL(dp), INTENT (in) :: a(xnitermax,1,nks),&
                           b(xnitermax,1,nks),&
                           xnorm(1,nks)
  REAL(dp), INTENT (in) :: e1s_ry 
  INTEGER,  INTENT (in) :: ncalcv(1,nks), ispectra
  LOGICAL,  INTENT (in) :: terminator
  
  !... Local variables
  INTEGER  :: i,ik,n,icoord           !loops
  INTEGER  :: lmax, lanczos_i, lanczos_f
  INTEGER  :: iestart, i_lanczos
  REAL(dp) :: alpha2
  REAL(dp) :: energy,de,mod_xgamma,xemax_ryd,xemin_ryd,xgamma_ryd
  REAL(dp) :: e0 ! in Ry
  REAL(dp) :: tmp_var
  REAL(dp) :: Intensity_coord(1,xnepoint,nspin)
  REAL(dp) :: continued_fraction 
  REAL(dp) :: paste_fermi, desmooth,t1,t2,f1,f2,df1,df2,poly(4) !CG 
  LOGICAL  :: first

  ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  !  constant and initialization
  ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  alpha2 = fpi/137.04

  desmooth = cut_desmooth/rytoev  ! This is in Rydberg

  e0 = xe0_ry 

  !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  !... Output file for the cross section 
  !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  IF( ionode ) THEN
     if(.not.two_edges.or.(two_edges.and.ispectra.eq.1)) then
        OPEN (unit=277,file='xanes.dat',form='formatted',status='unknown')
        REWIND(277) 
        !... writes input parameters in file 277
        WRITE(277,"('# Final state angular momentum:',1x,i3)") xang_mom

        IF (TRIM(ADJUSTL(gamma_mode)).EQ.'constant') THEN
           WRITE(277,"('# Broadening parameter (in eV):',1x,f8.3)") xgamma
        ELSE 
           WRITE(277,'("# Energy-dependent broadening parameter:")')
           IF (TRIM(ADJUSTL(gamma_mode)).EQ.'file') THEN
              WRITE(277,"('# -> using gamma_file:',1x,a50)") gamma_file
           ELSEIF (TRIM(ADJUSTL(gamma_mode)).EQ.'variable') THEN
              WRITE(277,"('# -> first, constant up to point (',f5.2,a1,f5.2,a)") &
                   gamma_energy(1),',',gamma_value(1),') [eV]'
              WRITE(277,"('# -> then, linear up to point (',f5.2,a1,f5.2,a)") &
                   gamma_energy(2),',',gamma_value(2),') [eV]'
              WRITE(277,"('# -> finally, constant up to xemax')")
           ENDIF
        ENDIF

        WRITE(277,"('# Absorbing atom type (xiabs):',i4)") xiabs

        IF(nspin.GT.1) THEN
           WRITE(277,"('# Energy (eV)   sigma_tot   sigma_up    sigma_down ')")
        ELSE
           WRITE(277,"('# Energy (eV)   sigma')")
        ENDIF
     endif
  ENDIF

  !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  !... Converts in Rydberg most of the relevant quantities
  !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  xemax_ryd = xemax/rytoeV + e0
  xemin_ryd = xemin/rytoeV + e0
  xgamma_ryd= xgamma/rytoeV

  !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  !... Calculates the continued fraction
  !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  Intensity_coord(:,:,:) = 0.d0

  de = (xemax_ryd-xemin_ryd) / REAL(xnepoint-1)

!<OB>
  if( .not. two_edges ) then
     lanczos_i = 1
     lanczos_f = n_lanczos
  else
     if( ispectra == 1 ) then
        lanczos_i = 1
        if( edge == 'L23' ) then
           lanczos_f = 2
        else if( edge == 'M45' ) then
           lanczos_f = 4
        else
           write(stdout,*) 'Output not yet programmed...'
        end if
     else
        lanczos_f = n_lanczos
        if( edge == 'L23' ) then
           lanczos_i = 3
        else if( edge == 'M45' ) then
           lanczos_i = 5
        else
           write(stdout,*) 'Output not yet programmed...'
        end if
     end if
  end if
  !<OB>
  


  IF (TRIM(gamma_mode).EQ.'constant') THEN

     IF(cut_occ_states) THEN
        ALLOCATE(memu(cut_nmemu,2))
        ALLOCATE(meml(cut_nmeml,2))
        iestart=(e0-xemin_ryd)/de
        do i_lanczos = lanczos_i, lanczos_f
          DO ik=1,nks

             first=.true.  ! to erase the memory of paste_fermi
             !<CG>
             t1=e0-desmooth
             f1=paste_fermi(t1,e0,a(1,i_lanczos,ik),b(1,i_lanczos,ik),xgamma_ryd,ncalcv(i_lanczos,ik)-1,terminator, first)
             df1=paste_fermi(t1-de,e0,a(1,i_lanczos,ik),b(1,i_lanczos,ik),xgamma_ryd,ncalcv(i_lanczos,ik)-1,terminator, first)
             df1=(f1-df1)/de
             t2=e0+desmooth
             f2=continued_fraction(a(1,i_lanczos,ik),b(1,i_lanczos,ik),t2,xgamma_ryd,ncalcv(i_lanczos,ik)-1,terminator)&
                  +paste_fermi(t2,e0,a(1,i_lanczos,ik),b(1,i_lanczos,ik),xgamma_ryd,ncalcv(i_lanczos,ik)-1,terminator, first)
             df2=continued_fraction(a(1,i_lanczos,ik),b(1,i_lanczos,ik),t2+de,xgamma_ryd,ncalcv(i_lanczos,ik)-1,terminator)&
                  +paste_fermi(t2+de,e0,a(1,i_lanczos,ik),b(1,i_lanczos,ik),xgamma_ryd,ncalcv(i_lanczos,ik)-1,terminator, first)
             df2=(df2-f2)/de
             CALL determine_polycut(t1,t2,f1,f2,df1,df2,poly) ! calculates interpolation polynome
  
             DO n=1,xnepoint
                energy=xemin_ryd+de*(n-1)
                IF ((energy-e0<desmooth).AND.(energy-e0>-desmooth)) THEN  ! interpolation 
                   tmp_var=poly(1)+poly(2)*energy+poly(3)*energy**2+poly(4)*energy**3
                   tmp_var=tmp_var*xnorm(i_lanczos,ik)*xnorm(i_lanczos,ik)
                   !</CG>
                ELSE
                   tmp_var=0.d0
                   IF (n>iestart) THEN
                      tmp_var=  &
                           continued_fraction(a(1,i_lanczos,ik),b(1,i_lanczos,ik),energy, &
                           xgamma_ryd,ncalcv(i_lanczos,ik)-1,terminator)*  &
                           xnorm(i_lanczos,ik)*xnorm(i_lanczos,ik)
                   ENDIF
                   tmp_var = tmp_var + paste_fermi(energy,e0,a(1,i_lanczos,ik),&
                        b(1,i_lanczos,ik),xgamma_ryd,ncalcv(i_lanczos,ik)-1,terminator, first) &
                        *xnorm(i_lanczos,ik)*xnorm(i_lanczos,ik)
                ENDIF
                Intensity_coord(i_lanczos,n,isk(ik)) = Intensity_coord(i_lanczos,n,isk(ik))+tmp_var*wk(ik)
             ENDDO
          ENDDO
        end do
        DEALLOCATE(memu)
        DEALLOCATE(meml)


     ELSE
        do i_lanczos = lanczos_i, lanczos_f
          DO ik=1,nks

             DO n=1,xnepoint
                energy=xemin_ryd+de*(n-1)
                tmp_var=  &
                     continued_fraction(a(1,i_lanczos,ik),b(1,i_lanczos,ik),&
                     energy,xgamma_ryd,ncalcv(i_lanczos,ik)-1, terminator)*  &
                     xnorm(i_lanczos,ik)*xnorm(i_lanczos,ik)
                Intensity_coord(i_lanczos,n,isk(ik)) = Intensity_coord(i_lanczos,n,isk(ik))+tmp_var*wk(ik)
             ENDDO
          ENDDO
        end do
     ENDIF

  ELSE ! nonconstant gamma

     IF(cut_occ_states) THEN
        ALLOCATE(memu(cut_nmemu,2))
        ALLOCATE(meml(cut_nmeml,2))
        iestart=(e0-xemin_ryd)/de
        do i_lanczos = lanczos_i, lanczos_f 
          DO ik=1,nks

             first=.true.  ! to erase the memory of paste_fermi
  
  
             xgamma_ryd=gamma_tab(iestart)
             t1=e0-desmooth
             f1=paste_fermi(t1,e0,a(1,i_lanczos,ik),b(1,i_lanczos,ik),xgamma_ryd,ncalcv(i_lanczos,ik)-1,terminator, first)
             df1=paste_fermi(t1-de,e0,a(1,i_lanczos,ik),b(1,i_lanczos,ik),xgamma_ryd,ncalcv(i_lanczos,ik)-1,terminator, first)
             df1=(f1-df1)/de
             t2=e0+desmooth
             f2=continued_fraction(a(1,i_lanczos,ik),b(1,i_lanczos,ik),t2,xgamma_ryd,ncalcv(i_lanczos,ik)-1,terminator)&
                  +paste_fermi(t2,e0,a(1,i_lanczos,ik),b(1,i_lanczos,ik),xgamma_ryd,ncalcv(i_lanczos,ik)-1,terminator, first)
             df2=continued_fraction(a(1,i_lanczos,ik),b(1,i_lanczos,ik),t2+de,xgamma_ryd,ncalcv(i_lanczos,ik)-1,terminator)&
                  +paste_fermi(t2+de,e0,a(1,i_lanczos,ik),b(1,i_lanczos,ik),xgamma_ryd,ncalcv(i_lanczos,ik)-1,terminator, first)
             df2=(df2-f2)/de
             CALL determine_polycut(t1,t2,f1,f2,df1,df2,poly)
  
             DO n=1,xnepoint
                energy=xemin_ryd+de*(n-1)
                xgamma_ryd=gamma_tab(n)
                IF ((energy-e0<desmooth).AND.(energy-e0>-desmooth)) THEN  ! interpolation
                   tmp_var=poly(1)+poly(2)*energy+poly(3)*energy**2+poly(4)*energy**3
                   tmp_var=tmp_var*xnorm(i_lanczos,ik)*xnorm(i_lanczos,ik)
                ELSE
                   tmp_var=0.d0
                   IF (n>iestart) THEN
                      tmp_var=  &
                           continued_fraction(a(1,i_lanczos,ik),b(1,i_lanczos,ik),&
                           energy,xgamma_ryd,ncalcv(i_lanczos,ik)-1,terminator)*  &
                           xnorm(i_lanczos,ik)*xnorm(i_lanczos,ik)
                   ENDIF
                   tmp_var = tmp_var + paste_fermi(energy,e0,a(1,i_lanczos,ik),&
                        b(1,i_lanczos,ik),xgamma_ryd,ncalcv(i_lanczos,ik)-1,terminator, first) &
                        *xnorm(i_lanczos,ik)*xnorm(i_lanczos,ik)
                ENDIF
                !            Intensity_tot(n)=Intensity_tot(n)+tmp_var*wk(ik)
                Intensity_coord(i_lanczos,n,isk(ik)) = Intensity_coord(i_lanczos,n,isk(ik))+tmp_var*wk(ik)
             ENDDO
          ENDDO
        end do
        DEALLOCATE(memu)
        DEALLOCATE(meml)


     ELSE
        do i_lanczos = lanczos_i, lanczos_f 
          DO ik=1,nks

             DO n=1,xnepoint
                energy=xemin_ryd+de*(n-1)
                xgamma_ryd=gamma_tab(n)
                tmp_var=  &
                     continued_fraction(a(1,i_lanczos,ik),b(1,i_lanczos,ik),&
                     energy,xgamma_ryd,ncalcv(i_lanczos,ik)-1, terminator)*  &
                     xnorm(i_lanczos,ik)*xnorm(i_lanczos,ik)
                Intensity_coord(i_lanczos,n,isk(ik)) = Intensity_coord(i_lanczos,n,isk(ik))+tmp_var*wk(ik)
             ENDDO
          ENDDO
       end do
     ENDIF

  ENDIF ! gamma_mode


  !... Considers the two cases of constant and non-constant broadening parameter
  !... Case 1: gamma is constant
!  IF (TRIM(gamma_mode).EQ.'constant') THEN
!
!     IF(cut_occ_states) THEN
!
!        ALLOCATE(memu(cut_nmemu,2))
!        ALLOCATE(meml(cut_nmeml,2))
!        iestart = (e0-xemin_ry)/de
!
!        DO ik = 1, nks
!           first = .true.  ! to erase the memory of paste_fermi
!           !<CG>
!           t1 = e0 - desmooth
!           f1 = paste_fermi(t1, e0, a(1,1,ik), b(1,1,ik),&
!                            xgamma_ry, ncalcv(1,ik)-1,   &
!                            terminator, first)
!           df1 = paste_fermi(t1-de, e0, a(1,1,ik), b(1,1,ik),&
!                             xgamma_ry, ncalcv(1,ik)-1,      &
!                             terminator, first)
!           df1 = (f1-df1)/de
!           t2 = e0 + desmooth
!           f2 = continued_fraction(a(1,1,ik), b(1,1,ik),         &
!                                   t2, xgamma_ry, ncalcv(1,ik)-1,&
!                                   terminator)                   &
!                + paste_fermi(t2, e0, a(1,1,ik), b(1,1,ik),      &
!                              xgamma_ry, ncalcv(1,ik)-1,         &
!                              terminator, first)
!           df2 = continued_fraction(a(1,1,ik), b(1,1,ik),             &
!                                    t2+de, xgamma_ry, ncalcv(1,ik)-1, &
!                                    terminator)                       &
!                + paste_fermi(t2+de, e0, a(1,1,ik), b(1,1,ik),        &
!                              xgamma_ry, ncalcv(1,ik)-1,              &
!                              terminator, first)
!           df2 = (df2-f2)/de
!
!           !... Calculates interpolation polynome
!           CALL determine_polycut(t1,t2,f1,f2,df1,df2,poly) 
!
!           DO n = 1, xnepoint
!              energy = xemin_ry + de*(n-1)
!              IF ((energy-e0<desmooth).AND.(energy-e0>-desmooth)) THEN  
!                 ! interpolation 
!                 tmp_var = poly(1) +           &
!                           poly(2)*energy +    &
!                           poly(3)*energy**2 + &
!                           poly(4)*energy**3
!                 tmp_var = tmp_var*xnorm(1,ik)*xnorm(1,ik)
!                 !</CG>
!              ELSE
!                 tmp_var = 0.d0
!                 IF (n > iestart) &
!                   tmp_var = continued_fraction(a(1,1,ik), b(1,1,ik),       &
!                                                energy, xgamma_ry,          &
!                                                ncalcv(1,ik)-1, terminator) &
!                              *xnorm(1,ik)*xnorm(1,ik)
!                 tmp_var = tmp_var +                                     &
!                           paste_fermi(energy, e0, a(1,1,ik), b(1,1,ik), &
!                                        xgamma_ry, ncalcv(1,ik)-1,       &
!                                        terminator, first)               &
!                           *xnorm(1,ik)*xnorm(1,ik)
!              ENDIF
!              Intensity_coord(1,n,isk(ik)) = Intensity_coord(1,n,isk(ik)) &
!                                             + tmp_var*wk(ik)
!           ENDDO
!
!        ENDDO
!
!        DEALLOCATE(memu)
!        DEALLOCATE(meml)
!
!     ELSE ! if occupied states are not cut
!        DO ik = 1, nks
!           DO n = 1, xnepoint
!              energy = xemin_ry + de*(n-1)
!              tmp_var= continued_fraction(a(1,1,ik), b(1,1,ik), &
!                                          energy, xgamma_ry,    &
!                                          ncalcv(1,ik)-1,       &
!                                          terminator)           &
!                       *xnorm(1,ik)*xnorm(1,ik)
!              Intensity_coord(1,n,isk(ik)) = Intensity_coord(1,n,isk(ik)) &
!                                             + tmp_var*wk(ik)
!           ENDDO
!        ENDDO
!     ENDIF
!
  !... Case 2: gamma is not constant (energy-dependent)
!  ELSE 
!
!     IF(cut_occ_states) THEN
!
!        ALLOCATE(memu(cut_nmemu,2))
!        ALLOCATE(meml(cut_nmeml,2))
!        iestart=(e0-xemin_ry)/de
!        DO ik=1,nks
!           first=.true.  ! to erase the memory of paste_fermi
!           xgamma_ry = gamma_tab(iestart)
!           t1 = e0 - desmooth
!           f1 = paste_fermi(t1, e0, a(1,1,ik), b(1,1,ik), &
!                            xgamma_ry, ncalcv(1,ik)-1,    &
!                            terminator, first)
!           df1 = paste_fermi(t1-de, e0, a(1,1,ik), b(1,1,ik), &
!                             xgamma_ry, ncalcv(1,ik)-1,       &
!                             terminator, first)
!           df1 = (f1-df1)/de
!           t2 = e0 + desmooth
!           f2 = continued_fraction(a(1,1,ik), b(1,1,ik),          &
!                                   t2, xgamma_ry, ncalcv(1,ik)-1, &
!                                   terminator)                    &
!                + paste_fermi(t2, e0, a(1,1,ik), b(1,1,ik),       &
!                              xgamma_ry, ncalcv(1,ik)-1,          &
!                              terminator, first)
!           df2 = continued_fraction(a(1,1,ik), b(1,1,ik),             &
!                                    t2+de, xgamma_ry, ncalcv(1,ik)-1, &
!                                    terminator)                       &
!                 + paste_fermi(t2+de, e0, a(1,1,ik), b(1,1,ik),       &
!                               xgamma_ry, ncalcv(1,ik)-1,             &
!                               terminator, first)
!           df2 = (df2 - f2)/de
!           CALL determine_polycut(t1,t2,f1,f2,df1,df2,poly)
!
!           DO n=1,xnepoint
!              energy = xemin_ry + de*(n-1)
!              xgamma_ry = gamma_tab(n)
!              IF ((energy-e0<desmooth).AND.(energy-e0>-desmooth)) THEN  
!                 ! interpolation
!                 tmp_var = poly(1) + &
!                           poly(2)*energy + &
!                           poly(3)*energy**2 + &
!                           poly(4)*energy**3
!                 tmp_var = tmp_var*xnorm(1,ik)*xnorm(1,ik)
!              ELSE
!                 tmp_var=0.d0
!                 IF (n>iestart) tmp_var = &
!                                   continued_fraction(a(1,1,ik), b(1,1,ik), &
!                                                      energy, xgamma_ry,    &
!                                                      ncalcv(1,ik)-1,       &
!                                                      terminator)           &
!                                   *xnorm(1,ik)*xnorm(1,ik)
!                 tmp_var = tmp_var + &
!                           paste_fermi(energy, e0, a(1,1,ik), b(1,1,ik),&
!                                       xgamma_ry, ncalcv(1,ik)-1,       &
!                                       terminator, first)               &
!                           *xnorm(1,ik)*xnorm(1,ik)
!              ENDIF
!              Intensity_coord(1,n,isk(ik)) = Intensity_coord(1,n,isk(ik)) &
!                                             + tmp_var*wk(ik)
!           ENDDO
!
!        ENDDO
!        DEALLOCATE(memu)
!        DEALLOCATE(meml)
!
!     ELSE ! if occupied states are not cut
!        DO ik=1,nks
!           DO n=1,xnepoint
!              energy = xemin_ry + de*(n-1)
!              xgamma_ry = gamma_tab(n)
!              tmp_var =  continued_fraction(a(1,1,ik), b(1,1,ik), &
!                                            energy, xgamma_ry,    &
!                                            ncalcv(1,ik)-1,       &
!                                            terminator)           &
!                         *xnorm(1,ik)*xnorm(1,ik)
!              Intensity_coord(1,n,isk(ik)) = Intensity_coord(1,n,isk(ik))&
!                                             + tmp_var*wk(ik)
!           ENDDO
!        ENDDO
!     ENDIF
!
!  ENDIF ! gamma_mode

  !  CALL poolreduce( nspin*xnepoint, Intensity_coord )

  !<CG>  replaces poolreduce
#ifdef __MPI
  CALL mp_sum ( Intensity_coord, inter_pool_comm )
#endif
  !</CG>

  !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  !... Writes the final cross section in file 277
  !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  IF(ionode) THEN
     IF(nspin == 1) THEN
        DO n=1,xnepoint
           energy = xemin_ryd + de*(n-1)
           Intensity_coord(:,n,:) = Intensity_coord(:,n,:) * &
                                    (energy+e1s_ry) *          &
                                    alpha2 
           WRITE(277,'(2f14.8)') (energy-e0)*rytoeV, Intensity_coord(1,n,1)
        ENDDO
     ELSEIF(nspin == 2) THEN
        DO n=1,xnepoint
           energy = xemin_ryd + de*(n-1)
           Intensity_coord(:,n,:) = Intensity_coord(:,n,:) * &
                                    (energy+e1s_ry)          * &
                                     alpha2 !
           WRITE(277,'(4f14.8)') (energy-e0)*rytoev, &
                Intensity_coord(1,n,1)+Intensity_coord(1,n,2),&
                Intensity_coord(1,n,1),Intensity_coord(1,n,2)
        ENDDO
     ENDIF

     CLOSE(277)
  ENDIF

END SUBROUTINE plot_xanes_dipole

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
SUBROUTINE plot_xanes_quadrupole(a,b,xnorm,ncalcv,terminator,e1s_ry)
  !--------------------------------------------------------------------------
  ! Calculates and plots the electric-quadrupole absorption cross section 
  ! as a continued fraction,
  ! from the a_i and b_i coefficients calculated for each k-point 
  !--------------------------------------------------------------------------
  USE kinds,      ONLY: DP
  USE constants,  ONLY: pi, rytoev
  USE xspectra,   ONLY: xang_mom, xemax, xemin, xiabs, xnepoint, &
                        xgamma, xonly_plot, xnitermax, xe0_ry
  !*apsi  USE uspp_param, ONLY : psd  !psd(ntypx) label for the atoms 
  USE klist,      ONLY: nkstot,& ! total number of k-points
                        nks,   & ! number of k-points per pool
                        xk,    & ! k-points coordinates
                        wk       ! k-points weight
  !USE ener,       ONLY: ef
  USE io_global,  ONLY: stdout,ionode  
  USE mp_global,  ONLY: inter_pool_comm
  USE lsda_mod,   ONLY: nspin,isk
  USE mp,         ONLY: mp_sum
  USE uspp_param, ONLY: upf
  USE gamma_variable_mod, ONLY : gamma_tab, gamma_mode, &
                                 gamma_file, gamma_value, gamma_energy
  USE cut_valence_green , ONLY : cut_occ_states, cut_desmooth, &
                                 cut_nmemu, cut_nmeml, memu, meml

  IMPLICIT NONE

  REAL(dp), INTENT (in) :: a(xnitermax,1,nks),&
                           b(xnitermax,1,nks),&
                           xnorm(1,nks)
  REAL(dp), INTENT (in) :: e1s_ry 
  INTEGER,  INTENT (in) :: ncalcv(1,nks)
  LOGICAL,  INTENT (in) :: terminator

  !... Local variables
  INTEGER  :: i,ik,n,icoord           !loops
  INTEGER  :: lmax
  INTEGER  :: iestart
  REAL(dp) :: alpha2,constantqua
  REAL(dp) :: energy,de,mod_xgamma,xemax_ry,xemin_ry,xgamma_ry
  REAL(dp) :: e0
  REAL(dp) :: tmp_var
  REAL(dp) :: Intensity_tot(xnepoint,nspin)
  REAL(dp) :: continued_fraction
  REAL(dp) :: paste_fermi,desmooth,t1,t2,f1,f2,df1,df2,poly(4)
  LOGICAL  :: first

  ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  !  constant and initialization
  ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  constantqua = pi/(137.04*137.04*137.04)

  alpha2=4.d0*pi/137.04

  desmooth=cut_desmooth/rytoev  ! This is in Rydberg
  
  e0 = xe0_ry 

  !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  !... Output file for the cross section 
  !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  IF( ionode ) THEN

     open (unit=277,file='xanes.dat',form='formatted',status='unknown')

     REWIND(277) ! to be at the initial point of file 277
     !... writes input parameters in file 277
     WRITE(277,"('# final state angular momentum:',1x,i3)") xang_mom

     IF (TRIM(ADJUSTL(gamma_mode)).EQ.'constant') THEN
        WRITE(277,"('# Broadening parameter (in eV):',1x,f8.3)") xgamma
     ELSE
        WRITE(277,'("# Energy-dependent broadening parameter:")')
        IF (TRIM(ADJUSTL(gamma_mode)).EQ.'file') THEN
           WRITE(277,"('# -> using gamma_file:',1x,a50)") gamma_file
        ELSEIF (TRIM(ADJUSTL(gamma_mode)).EQ.'variable') THEN
           WRITE(277,"('# -> first, constant up to point (',f5.2,a1,f5.2,a)") &
                 gamma_energy(1),',',gamma_value(1),') [eV]'
           WRITE(277,"('# -> then, linear up to point (',f5.2,a1,f5.2,a)") &
                 gamma_energy(2),',',gamma_value(2),') [eV]'
           WRITE(277,"('# -> finally, constant up to xemax')")
        ENDIF
     ENDIF

     WRITE(277,"('# Absorbing atom type (xiabs):',i4)") xiabs

     IF(nspin.GT.1) THEN
        WRITE(277,"('# Energy (eV)   sigma_tot   sigma_up    sigma_down ')")
     ELSE
        WRITE(277,"('# Energy (eV)   sigma')")
     ENDIF

  ENDIF

  !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  !... Converts in Rydberg most of the relevant quantities
  !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  xemax_ry = xemax/rytoev + e0
  xemin_ry = xemin/rytoev + e0
  xgamma_ry= xgamma/rytoev

  !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  !... Calculates the continued fraction
  !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  Intensity_tot(:,:)=0.d0

  de = (xemax_ry-xemin_ry) / REAL(xnepoint-1)

  !... Considers the two cases of constant and non-constant broadening parameter
  !... Case 1: gamma is constant
  IF (TRIM(gamma_mode).EQ.'constant') THEN

     IF(cut_occ_states) THEN

        ALLOCATE(memu(cut_nmemu,2))
        ALLOCATE(meml(cut_nmeml,2))
        DO ik=1,nks
           iestart = (e0 - xemin_ry)/de
           first = .true.
           t1 = e0 - desmooth
           f1 = paste_fermi(t1, e0, a(1,1,ik), b(1,1,ik), &
                            xgamma_ry, ncalcv(1,ik)-1,    &
                            terminator, first)
           df1 = paste_fermi(t1-de, e0, a(1,1,ik), b(1,1,ik), &
                             xgamma_ry, ncalcv(1,ik)-1,       &
                             terminator, first)
           df1 = (f1 - df1)/de
           t2 = e0 + desmooth
           f2 = continued_fraction(a(1,1,ik), b(1,1,ik), t2,  &
                                   xgamma_ry, ncalcv(1,ik)-1, &
                                   terminator)                &
                + paste_fermi(t2, e0, a(1,1,ik), b(1,1,ik),   &
                              xgamma_ry, ncalcv(1,ik)-1,      &
                              terminator, first)
           df2 = continued_fraction(a(1,1,ik), b(1,1,ik), t2+de, &
                                    xgamma_ry, ncalcv(1,ik)-1,   &
                                    terminator)                  &
                 + paste_fermi(t2+de, e0, a(1,1,ik), b(1,1,ik),  &
                               xgamma_ry,ncalcv(1,ik)-1,terminator, first)
           df2=(df2-f2)/de
           CALL determine_polycut(t1,t2,f1,f2,df1,df2,poly)

           DO n=1,xnepoint
              energy = xemin_ry + de*(n-1)
              IF ((energy-e0<desmooth).AND.(energy-e0>-desmooth)) THEN  
              !interpolation
                 tmp_var = poly(1) + &
                           poly(2)*energy + &
                           poly(3)*energy**2 + &
                           poly(4)*energy**3
                 tmp_var = tmp_var*xnorm(1,ik)*xnorm(1,ik)
              ELSE
                 tmp_var=0.d0
                 IF (n>iestart) tmp_var = &
                                   continued_fraction(a(1,1,ik), b(1,1,ik),&
                                                      energy, xgamma_ry,   &
                                                      ncalcv(1,ik)-1,      &
                                                      terminator)          &
                                   *xnorm(1,ik)*xnorm(1,ik)
                 tmp_var = tmp_var + &
                           paste_fermi(energy, e0, a(1,1,ik), b(1,1,ik), &
                                       xgamma_ry, ncalcv(1,ik)-1,        &
                                       terminator, first)                &
                           *xnorm(1,ik)*xnorm(1,ik)
              ENDIF
              Intensity_tot(n,isk(ik)) = Intensity_tot(n,isk(ik)) &
                                         + tmp_var*wk(ik)
           ENDDO
        ENDDO
        DEALLOCATE(memu)
        DEALLOCATE(meml)

     ELSE ! if occupied states are not cut
        DO ik=1,nks
           DO n=1,xnepoint
              energy = xemin_ry + de*(n-1)
              tmp_var= continued_fraction(a(1,1,ik), b(1,1,ik), &
                                          energy, xgamma_ry,    &
                                          ncalcv(1,ik)-1,       &
                                          terminator)           &
                       *xnorm(1,ik)*xnorm(1,ik)
              Intensity_tot(n,isk(ik)) = Intensity_tot(n,isk(ik)) &
                                         + tmp_var*wk(ik)
           ENDDO
        ENDDO
     ENDIF
 
  !... Case 2: gamma is not constant (energy-dependent)
  ELSE 

     IF(cut_occ_states) THEN ! if occupied states are cut
        ALLOCATE(memu(cut_nmemu,2))
        ALLOCATE(meml(cut_nmeml,2))
        iestart = (e0-xemin_ry)/de
        DO ik=1,nks
           first = .true. ! to erase memory of paste_fermi
           xgamma_ry = gamma_tab(iestart)
           t1 = e0 - desmooth
           f1 = paste_fermi(t1, e0, a(1,1,ik), b(1,1,ik), &
                            xgamma_ry, ncalcv(1,ik)-1,    &
                            terminator, first)
           df1 = paste_fermi(t1-de, e0, a(1,1,ik), b(1,1,ik), &
                             xgamma_ry, ncalcv(1,ik)-1,       &
                             terminator, first)
           df1 = (f1-df1) / de
           t2 = e0 + desmooth
           f2 = continued_fraction(a(1,1,ik), b(1,1,ik), t2,  &
                                   xgamma_ry, ncalcv(1,ik)-1, &
                                   terminator)                &
                + paste_fermi(t2, e0, a(1,1,ik), b(1,1,ik),   &
                              xgamma_ry, ncalcv(1,ik)-1,      &
                              terminator, first)
           df2 = continued_fraction(a(1,1,ik), b(1,1,ik), t2+de, &
                                    xgamma_ry, ncalcv(1,ik)-1,   &
                                    terminator)                  &
                 + paste_fermi(t2+de, e0, a(1,1,ik), b(1,1,ik),  &
                               xgamma_ry, ncalcv(1,ik)-1,        &
                               terminator, first)
           df2 = (df2-f2) / de
           CALL determine_polycut(t1,t2,f1,f2,df1,df2,poly)

           DO n=1,xnepoint
              energy = xemin_ry + de*(n-1)
              xgamma_ry = gamma_tab(n)
              IF ((energy-e0<desmooth).AND.(energy-e0>-desmooth)) THEN  
                 ! interpolation
                 tmp_var = poly(1) + &
                           poly(2)*energy + &
                           poly(3)*energy**2 + &
                           poly(4)*energy**3
                 tmp_var = tmp_var*xnorm(1,ik)*xnorm(1,ik)
              ELSE
                 tmp_var = 0.d0
                 IF (n>iestart) tmp_var = &
                                  continued_fraction(a(1,1,ik), b(1,1,ik), &
                                                     energy, xgamma_ry,    &
                                                     ncalcv(1,ik)-1,       &
                                                     terminator)           &
                                  *xnorm(1,ik)*xnorm(1,ik)
                 tmp_var = tmp_var + &
                           paste_fermi(energy, e0, a(1,1,ik),b(1,1,ik),&
                                       xgamma_ry, ncalcv(1,ik)-1,      &
                                       terminator, first)              &
                           *xnorm(1,ik)*xnorm(1,ik)
              ENDIF
              Intensity_tot(n,isk(ik)) = Intensity_tot(n,isk(ik)) &
                                         + tmp_var*wk(ik)
           ENDDO
        ENDDO
        DEALLOCATE(memu)
        DEALLOCATE(meml)

     ELSE ! if occupied states are not cut
        DO ik=1,nks
           DO n=1,xnepoint
              energy = xemin_ry + de*(n-1)
              xgamma_ry = gamma_tab(n)
              tmp_var = continued_fraction(a(1,1,ik), b(1,1,ik), &
                                           energy, xgamma_ry,    &
                                           ncalcv(1,ik)-1,       &
                                           terminator)           &
                        *xnorm(1,ik)*xnorm(1,ik)
              Intensity_tot(n,isk(ik)) = Intensity_tot(n,isk(ik)) &
                                         + tmp_var*wk(ik)
           ENDDO
        ENDDO
     ENDIF

  ENDIF ! gammma_mode

  !  CALL poolreduce( nspin*xnepoint, Intensity_tot )

  !<CG>  replaces poolreduce
#ifdef __MPI
  CALL mp_sum ( Intensity_tot, inter_pool_comm )
#endif
  !</CG>

  !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  !... Writes the final cross section in file 277
  !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  IF(ionode) THEN
     IF(nspin == 1) THEN
        DO n=1,xnepoint
           energy = xemin_ry + de*(n-1)
           !Intensity_tot(n,:)=Intensity_tot(n,:) &
           !                   *(energy+e_1s)*(energy+e_1s)*(energy+e_1s) &
           !                   * constantqua      !normalized
           Intensity_tot(n,:) = Intensity_tot(n,:) &
                                * (energy+e1s_ry)**3 &
                                * constantqua     
           WRITE(277,'(2f14.8)') (energy-e0)*rytoev, Intensity_tot(n,:)
        ENDDO
     ELSEIF(nspin == 2) THEN
        DO n=1,xnepoint
           energy = xemin_ry + de*(n-1)
           !Intensity_tot(n,:) = Intensity_tot(n,:) &
           !                     *(energy+e_1s)*(energy+e_1s)*(energy+e_1s) &
           !                     *constantqua     !normalized
           Intensity_tot(n,:) = Intensity_tot(n,:) &
                                * (energy+e1s_ry)**3 &
                                * constantqua     
           WRITE(277,'(4f14.8)') (energy-e0)*rytoev, &
                                  Intensity_tot(n,1)+Intensity_tot(n,2),&
                                  Intensity_tot(n,1), Intensity_tot(n,2)
        ENDDO
     ENDIF

     CLOSE(277)
  ENDIF
  

END SUBROUTINE plot_xanes_quadrupole

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
SUBROUTINE define_index_arrays(paw_iltonhb)
  !----------------------------------------------------------------------------
  USE paw_gipaw,   ONLY: paw_lmaxkb, paw_recon
  USE ions_base,   ONLY: ntyp => nsp
  USE parameters,  ONLY: lmaxx
  USE xspectra_paw_variables, ONLY : xspectra_paw_nhm ! CG

  IMPLICIT NONE
  ! Arguments
  INTEGER :: paw_iltonhb(0:paw_lmaxkb,xspectra_paw_nhm,ntyp) ! CG

  ! Local
  INTEGER :: nt,ih,nb,l
  INTEGER :: ip_per_l(0:lmaxx) 

  DO nt = 1, ntyp
     ih = 1
     ip_per_l(:) = 0
     DO nb = 1, paw_recon(nt)%paw_nbeta
        l = paw_recon(nt)%aephi(nb)%label%l
        ip_per_l(l) = ip_per_l(l) + 1
        paw_iltonhb(l,ip_per_l(l),nt) = ih
        ih = ih + 2*l + 1
     ENDDO
  ENDDO

END SUBROUTINE define_index_arrays

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
FUNCTION lastterm(a,b,g)
  !----------------------------------------------------------------------------
  USE kinds, ONLY: dp
  IMPLICIT NONE
  REAL(dp)    :: a, b, g, y1, y2, z1, z2, r
  COMPLEX(dp) :: lastterm

  y1 = a*a - g*g - 4*b
  y2 = -2*a*g
  r  = 0.5*SQRT(y1*y1 + y2*y2)

  IF (g<0) THEN
     z1 =  a/2 + 0.5*SIGN(SQRT(y1/2 + r),y2)
     z2 = -g/2 + 0.5*SQRT(-y1/2 + r)
  ELSE
     z1 =  a/2 - 0.5*SIGN(SQRT(y1/2 + r),y2)
     z2 = -g/2 - 0.5*SQRT(-y1/2+r)
  ENDIF
 
 lastterm = CMPLX(z1,z2,kind=DP)

END FUNCTION lastterm

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
FUNCTION paste_fermi(e,e0,a,b,gamma,m,term,first)
  !----------------------------------------------------------------------------
  USE kinds,     ONLY: dp
  USE constants, ONLY: tpi
  USE xspectra,  ONLY: xnitermax, xcheck_conv
  USE cut_valence_green, ONLY: cut_ierror, cut_stepu, cut_stepl, &
                               cut_startt, cut_tinf, cut_tsup,   &
                               cut_nmemu, cut_nmeml, memu, meml
  IMPLICIT NONE

  REAL(dp) :: paste_fermi
  ! Arguments
  INTEGER  :: m
  REAL(dp) :: a(xnitermax)
  REAL(dp) :: b(xnitermax)
  REAL(dp) :: gamma
  REAL(dp) :: e, e0
  LOGICAL  :: term, first
  ! Local
  COMPLEX(dp) :: green, y, dy, c1, c2, e1, e2
  REAL(dp)    :: t, dt, t1, ta, tb
  INTEGER, save :: n1
  INTEGER, save :: n2
  INTEGER :: nn1, nn2

  IF (first) THEN
     memu(:,:) = (0.d0,0.d0)
     meml(:,:) = (0.d0,0.d0) 
     n1 = 0
     n2 = 0
     first = .false.
  ENDIF

  dy = cut_ierror + 1.0
  y = 0.d0

  nn1 = 1
  nn2 = 1

  t1 = 0.5773502692

  t = cut_startt

  DO WHILE ((abs(dy)>cut_ierror).OR.(t<cut_tsup))
     dt = cut_stepu*t
     ta = t + dt*(1-t1)/2
     tb = t + dt*(1+t1)/2
     e1 = CMPLX(e0,ta,kind=DP)
     e2 = CMPLX(e0,tb,kind=DP)

     IF (nn1>n1) THEN
        c1 = green(a,b,e1,m,term)
        c2 = green(a,b,e2,m,term)
        IF (nn1<cut_nmemu) THEN
           memu(nn1,1) = c1
           memu(nn1,2) = c2
           n1 = nn1
        ENDIF
     ELSE
        c1 = memu(nn1,1)
        c2 = memu(nn1,2)
     ENDIF

     dy = (dt/2) * &
          ( c1/CMPLX(e0-e,ta-gamma,kind=DP)          &
           + CONJG(c1)/CMPLX(e0-e,-ta-gamma,kind=DP) &
           + c2/CMPLX(e0-e,tb-gamma,kind=DP)         &
           + CONJG(c2)/CMPLX(e0-e,-tb-gamma,kind=DP) )
     y = y + dy
     t = t + dt
     nn1 = nn1 + 1
  ENDDO

  t = cut_startt
  dy = cut_ierror + 1

  DO WHILE((abs(dy)>cut_ierror).OR.(t>cut_tinf))
     dt = cut_stepl * t
     ta = t - dt*(1-t1)/2
     tb = t - dt*(1+t1)/2
     e1 = CMPLX(e0,ta,kind=DP)
     e2 = CMPLX(e0,tb,kind=DP)

     IF (nn2>n2) THEN
        c1 = green(a,b,e1,m,term)
        c2 = green(a,b,e2,m,term)
        IF (nn2<cut_nmeml) THEN
           meml(nn2,1) = c1
           meml(nn2,2) = c2
           n2 = nn2
        ENDIF
     ELSE
        c1 = meml(nn2,1)
        c2 = meml(nn2,2)
     ENDIF

     dy = (dt/2) * &
          ( c1/CMPLX(e0-e,ta-gamma,kind=DP)           &
            + CONJG(c1)/CMPLX(e0-e,-ta-gamma,kind=DP) &
            + c2/CMPLX(e0-e,tb-gamma,kind=DP)         &
            + CONJG(c2)/CMPLX(e0-e,-tb-gamma,kind=DP))
     y = y + dy
     t = t - dt
     nn2 = nn2 + 1
  ENDDO

  paste_fermi = AIMAG(y)/tpi

END FUNCTION paste_fermi

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
FUNCTION green(a,b,e,m,term)
  !----------------------------------------------------------------------------
  USE kinds,    ONLY : dp
  USE xspectra, ONLY : xnitermax, xcheck_conv

  IMPLICIT NONE

  COMPLEX(dp) :: green
  ! Arguments
  INTEGER  :: m
  REAL(dp) :: a(xnitermax)
  REAL(dp) :: b(xnitermax)
  COMPLEX(dp) :: e
  LOGICAL :: term
  ! Local
  INTEGER :: i, p, q
  COMPLEX(dp) :: res, lastterm 
  REAL(dp) :: aa, bb

  q = xcheck_conv/2
  IF (term) THEN
     aa = 0.0
     bb = 0.0
     DO p = 1, q
        aa = aa + a(m-p)
        bb = bb + b(m-p)
     ENDDO
     aa = aa/q
     bb = bb/q

     res = lastterm(aa-REAL(e), bb*bb, AIMAG(e))
  ELSE
     res = CMPLX(a(m)-REAL(e),AIMAG(e),kind=DP)
  ENDIF
  DO i = 1, m-1
     res = a(m-i) - e - b(m-i)*b(m-i)/res
  ENDDO

  green = 1/res

END FUNCTION green

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
SUBROUTINE check_paw_projectors(xiabs)
  !----------------------------------------------------------------------------
  USE kinds,           ONLY: DP
  USE constants,       ONLY: pi
  USE paw_gipaw,       ONLY: paw_lmaxkb, paw_recon
  USE xspectra_paw_variables, ONLY: xspectra_paw_nhm
  USE atom,            ONLY: rgrid, msh
  !  USE atom,  ONLY : &
  !       mesh,     &!mesh(ntypx) number of mesh points              
  !       msh ,     &!msh(ntypx)the point at rcut=end of radial integration 
  !       r, rab
  USE ions_base,       ONLY: ntyp => nsp
  USE io_global,       ONLY: stdout
  USE radin_mod

  IMPLICIT NONE
  INTEGER, INTENT(IN) :: xiabs
  ! Local  variables
  INTEGER  :: nr, nrc, ip, jp, lmax, l, ip_l, jtyp, n1, n2, nrs, ndm, ih, jh
  REAL(dp) :: overlap, rexx, overlap2
  REAL (dp), ALLOCATABLE :: aux(:), f(:,:)
  REAL(dp) , ALLOCATABLE :: s(:,:), e(:), v(:,:)

  ALLOCATE(aux(rgrid(xiabs)%mesh)) !allocation too big, it needs only up to msh
  ALLOCATE(f(rgrid(xiabs)%mesh,2)) !allocation too big, it needs only up to msh

  WRITE(stdout,'(/,5x,a)')&
  '-------------------------------------------------------------------------'
  WRITE(stdout,'(5x,a)')  &
  '                         Verification of the PAW relations '
  WRITE(stdout,'(5x,a,/)')&
  '-------------------------------------------------------------------------'

  WRITE(stdout,'(8x,a)') 'atom type   total number of projectors'
  DO jtyp=1,ntyp
     WRITE (stdout,'(13x,i4,3x,i4)') jtyp, paw_recon(jtyp)%paw_nbeta
  ENDDO
  WRITE(stdout,*)
  
  !... I calculate maximum l
  !
  lmax=0
  DO ip = 1, paw_recon(xiabs)%paw_nbeta
     IF(paw_recon(xiabs)%psphi(ip)%label%l.GT.lmax) &
          lmax = paw_recon(xiabs)%psphi(ip)%label%l
  ENDDO

  WRITE(stdout,'(8x,a)') 'atom type    l   number of projectors per ang. mom.'
  DO jtyp = 1, ntyp
     DO l = 0, lmax
        WRITE(stdout,'(13x,i4,3x,i2,3x,i4)') jtyp, l, paw_recon(jtyp)%paw_nl(l)
     ENDDO
  ENDDO
  WRITE(stdout,*)

  !... We calculate the overlaps between partial waves and projectors
  !    to see if they are equal to the Croneker delta.

  nr = msh(xiabs)  ! extended up to all the NON ZERO points in the mesh.

  WRITE(stdout,'(8x,a)')'Overlaps between partial waves and projectors (radial)'
  WRITE(stdout,'(8x,a)')'------------------------------------------------------'
  WRITE(stdout,*)
  WRITE(stdout,'(8x,a)') &
                    '< \tilde{phi}_{l,n} | \tilde{p}_{l,nn} > = delta_{n,nn}  ?'
  WRITE(stdout,*) 

  DO ip = 1, paw_recon(xiabs)%paw_nbeta
     DO jp = 1, paw_recon(xiabs)%paw_nbeta
        IF(paw_recon(xiabs)%psphi(ip)%label%l .EQ. &
                        paw_recon(xiabs)%psphi(jp)%label%l) THEN
           nrc=Count(rgrid(xiabs)%r(1:nr)<=paw_recon(xiabs)%psphi(ip)%label%rc)
           IF(nrc > nr) THEN
              WRITE(stdout,'(8x,a,i8,a,i8)') 'STOP: nrc=', nrc,' > nr=', nr
              CALL errore ( "nrc > nr", "xanes_dipole", 0 )
           ENDIF
           aux(1:nrc) = paw_recon(xiabs)%psphi(ip)%psi(1:nrc) &
                        * paw_recon(xiabs)%paw_betar(1:nrc,jp)
           aux(nrc+1:nr) = 0.d0
           WRITE(stdout,'(8x,"<tilde{phi}_",2i2,10X,"|tilde{p}_",2i2,">=",1f14.8)')  &
                ip,paw_recon(xiabs)%psphi(ip)%label%l,jp, &
                paw_recon(xiabs)%psphi(jp)%label%l, &
                para_radin(aux(1:nr),rgrid(xiabs)%r(1:nr),nr)
        ENDIF
     ENDDO
  ENDDO

  WRITE(stdout,*)

  WRITE(stdout,'(8x,a)')'Checking normalization of pseudo,ae wfc and projectors'
  WRITE(stdout,'(8x,a)')'(radial part only, integral up to r_c)'
  WRITE(stdout,'(8x,a)')'------------------------------------------------------'
  WRITE(stdout,*)
  WRITE(stdout,'(8x,a)') 'l,   n, |proj|^2, |pswf|^2 , |aewf|^2'
  DO l = 0, lmax
     DO ip = 1, paw_recon(xiabs)%paw_nbeta
        IF(paw_recon(xiabs)%psphi(ip)%label%l.EQ.l) THEN
           nrc =Count(rgrid(xiabs)%r(1:nr)<=paw_recon(xiabs)%psphi(ip)%label%rc)
           aux(1:nrc) = paw_recon(xiabs)%paw_betar(1:nrc,ip) &
                      * paw_recon(xiabs)%paw_betar(1:nrc,ip)
           overlap = para_radin(aux(1:nrc),rgrid(xiabs)%r(1:nrc),nrc)
           aux(1:nrc) = paw_recon(xiabs)%aephi(ip)%psi(1:nrc) &
                      * paw_recon(xiabs)%aephi(ip)%psi(1:nrc)
           overlap2 = para_radin(aux(1:nrc),rgrid(xiabs)%r(1:nrc),nrc)
           aux(1:nrc) = paw_recon(xiabs)%psphi(ip)%psi(1:nrc) &
                      * paw_recon(xiabs)%psphi(ip)%psi(1:nrc)
           rexx = para_radin(aux(1:nrc),rgrid(xiabs)%r(1:nrc),nrc)
           WRITE(stdout,'(8x,2i4,3f14.8)')l,ip,overlap,overlap2,rexx
        ENDIF
     ENDDO
  ENDDO
  WRITE(stdout,*)
  
  GOTO 323

  WRITE(stdout,'(8x,a)') '<phi|chi> = \sum_nl <phi|phi_l> <p_l|chi>_nrc ?'
  WRITE(stdout,'(8x,a)') '-----------------------------------------------'
  WRITE(stdout,'(8x,a)') &
        'WARNING: this test assumes a form of the phi/chi function'
  
  !  DO l=0,lmax
  DO l = 1, 1
     ip_l = 0
     DO ip = 1, paw_recon(xiabs)%paw_nbeta
        IF(ip_l.EQ.0.AND.paw_recon(xiabs)%psphi(ip)%label%l.EQ.l) ip_l = ip
     ENDDO
     
     f(:,:) = 0.d0
     DO ip = 1, paw_recon(xiabs)%paw_nbeta
        IF(paw_recon(xiabs)%psphi(ip)%label%l.EQ.l) THEN
           f(1:nr,1) = f(1:nr,1) + &
                       paw_recon(xiabs)%psphi(ip)%psi(1:nr)/REAL(ip,dp)
           f(1:nr,2) = f(1:nr,2) + &
                       1.123*paw_recon(xiabs)%psphi(ip)%psi(1:nr)/REAL(ip,dp)
        ENDIF
     ENDDO
     rexx = 0.d0
     DO ip = 1, paw_recon(xiabs)%paw_nbeta
        IF(paw_recon(xiabs)%psphi(ip)%label%l.EQ.l) THEN
           nrc=Count(rgrid(xiabs)%r(1:nr)<=paw_recon(xiabs)%psphi(ip)%label%rc)
           IF(nrc > nr) THEN
              WRITE(stdout,'(8x,a,i8,a,i8)') 'STOP: nrc=', nrc,' > nr=', nr
              CALL errore ( "nrc > nr", "xanes_dipole", 0 )
           ENDIF
           aux(1:nrc) = f(1:nrc,1)*paw_recon(xiabs)%paw_betar(1:nrc,ip)
           overlap = para_radin(aux(1:nrc),rgrid(xiabs)%r(1:nrc),nrc)
           aux(1:nrc) = f(1:nrc,2)*paw_recon(xiabs)%psphi(ip)%psi(1:nrc)
           overlap = overlap*para_radin(aux(1:nrc),rgrid(xiabs)%r(1:nrc),nrc)
           WRITE(stdout,'(8x,"overlap(l=",1i2,",n=",1i2,")= ",1f14.8)')&
                          l, ip, overlap
           rexx = rexx + overlap
        ENDIF
     ENDDO
     aux(1:nr) = f(1:nr,1)*f(1:nr,2)
     WRITE(stdout,'(8x,"sum/overlap=",1f14.8)') &
                    rexx/para_radin(aux,rgrid(xiabs)%r(1:nr),nrc)
     WRITE(stdout,'(8x,"sum projectors=",1f14.8," overlap=",1f14.8)') &
                    rexx,para_radin(aux,rgrid(xiabs)%r(1:nr),nrc)
     WRITE(stdout,*)
  ENDDO
  ! ENDDO
  !WRITE(stdout,*)
  !WRITE(stdout,*) '================================================================'
  
323 CONTINUE
  !
  !  Check linear dependence of projectors
  !

  WRITE(stdout,'(8x,a)') 'Checking linear dependence of projectors'
  WRITE(stdout,'(8x,a)') '----------------------------------------'
  WRITE(stdout,*)

  DEALLOCATE(aux)

  ndm = MAXVAL (msh(1:ntyp))

  ALLOCATE(aux(ndm))

  DO l = 0, paw_lmaxkb
     IF (paw_recon(xiabs)%paw_nl(l)>0) THEN
        ALLOCATE (s(paw_recon(xiabs)%paw_nl(l),paw_recon(xiabs)%paw_nl(l)))
        ALLOCATE (e(paw_recon(xiabs)%paw_nl(l)),v(paw_recon(xiabs)%paw_nl(l),&
                    paw_recon(xiabs)%paw_nl(l)))
        DO ih = 1, paw_recon(xiabs)%paw_nl(l)
           n1 = paw_recon(xiabs)%paw_iltonh(l,ih)
           nrc = paw_recon(xiabs)%psphi(n1)%label%nrc
           nrs = paw_recon(xiabs)%psphi(n1)%label%nrs
           DO jh = 1, paw_recon(xiabs)%paw_nl(l)
              n2 = paw_recon(xiabs)%paw_iltonh(l,jh)
              CALL step_f(aux,paw_recon(xiabs)%psphi(n1)%psi(1:msh(xiabs)) * &
                   paw_recon(xiabs)%psphi(n2)%psi(1:msh(xiabs)), &
                   rgrid(xiabs)%r(1:msh(xiabs)),nrs,nrc, 1.d0, msh(xiabs) )
              CALL simpson( msh(xiabs), aux, rgrid(xiabs)%rab(1), s(ih,jh))
           ENDDO
        ENDDO
  
        WRITE(stdout,'(8x,"atom type:",1i4)') xiabs
        WRITE(stdout,'(8x,"number of projectors projector  =",1i3, " angular momentum=",1i4)') &
                    paw_recon(xiabs)%paw_nl(l),l

        DO ih = 1, paw_recon(xiabs)%paw_nl(l)
           WRITE(stdout,'(8x,10f14.8)') (s(ih,jh),jh=1,paw_recon(xiabs)%paw_nl(l)) 
        ENDDO
        WRITE(stdout,'(8x,a)') 'Eigenvalues S matrix:'

        IF(paw_recon(xiabs)%paw_nl(l).EQ.1) THEN
           WRITE(stdout,'(10x,1i4,1f14.8)') 1,s(1,1)
        ELSE 
           CALL rdiagh(paw_recon(xiabs)%paw_nl(l), s, &
                       paw_recon(xiabs)%paw_nl(l) , e, v )
           DO ih=1,paw_recon(xiabs)%paw_nl(l)
              WRITE(stdout,'(10x,1i4,1f14.8)') ih,e(ih)
           ENDDO
        ENDIF
        WRITE(stdout,*)
        DEALLOCATE(s,e,v)
     ENDIF
  ENDDO
  !WRITE(stdout,*) '================================================================'

  DEALLOCATE(aux,f)

END SUBROUTINE check_paw_projectors

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
SUBROUTINE read_save_file(a,b,xnorm,ncalcv,x_save_file,core_energy)
  !----------------------------------------------------------------------------
  ! This routine reads the x_save_file (default name: xanes.sav)
  ! (automatically named prefix_K.save for a K edge, in the next version). 
  ! This routine is used only if xonly_plot=.true.
  !----------------------------------------------------------------------------
  USE kinds,       ONLY: DP
  USE constants,   ONLY: rytoev
  USE klist,       ONLY: nks, nkstot
  USE xspectra,    ONLY: xnitermax, xang_mom, xiabs,   &
                         n_lanczos, save_file_version, &
                         save_file_kind, calculated,   &
                         xe0, xe0_default, xe0_ry
  USE ener,        ONLY: ef
  USE io_global,   ONLY: stdout, ionode
  USE lsda_mod,    ONLY: nspin, lsda

  IMPLICIT NONE
  ! Arguments
  REAL(dp), INTENT (INOUT) ::  a(xnitermax,n_lanczos,nks)
  REAL(dp), INTENT (INOUT) ::  b(xnitermax,n_lanczos,nks)     
  REAL(dp), INTENT (INOUT) ::  xnorm(n_lanczos,nks) 
  INTEGER, INTENT (INOUT)  ::  ncalcv(n_lanczos,nks)
  CHARACTER(LEN=256), INTENT (IN) :: x_save_file
  REAL(dp), INTENT (OUT)   ::  core_energy
  ! Local variables
  ! NB: The '_r' suffix means 'read in x_save_file'
  INTEGER  :: ierr, nkstot_r
  INTEGER  :: xm_r, nc_r, ncomp_max
  INTEGER  :: i, j, k, ncalcv_max
  INTEGER  :: calculated_all(n_lanczos,nkstot)
  INTEGER, ALLOCATABLE :: ncalcv_all(:,:)
  REAL(dp) :: xepsilon_r(3), xkvec_r(3)
  REAL(dp), ALLOCATABLE :: a_all(:,:), b_all(:,:), xnorm_all(:,:), aux(:,:)

  ALLOCATE(a_all(xnitermax,nkstot))
  ALLOCATE(b_all(xnitermax,nkstot))
  ALLOCATE(xnorm_all(n_lanczos,nkstot))
  ALLOCATE(ncalcv_all(n_lanczos,nkstot))
  ALLOCATE(aux(xnitermax,nks))

  a(:,:,:) = 0.d0
  b(:,:,:) = 0.d0
  xnorm(:,:) = 0.d0
  ncalcv(:,:) = 0
  calculated_all(:,:) = 0


  OPEN ( UNIT = 10, FILE = x_save_file, FORM = 'FORMATTED', &
       STATUS = 'UNKNOWN', IOSTAT = ierr )
  CALL errore( 'iosys', 'x_save_file ' // TRIM( x_save_file ) // &
       & ' not found' , ierr )
  WRITE(stdout,'(5x,"x_save_file name: ",a)') TRIM( x_save_file )
  REWIND(10)

  IF (save_file_version == 0) then
     WRITE(stdout,'(5x,a)') 'x_save_file version: old'
  !ELSEIF(save_file_version == 1) then
  !   WRITE(stdout,'(5x,a)') 'x_save_file version: 1'
  ELSE
     WRITE(stdout,'(5x,a,i3)') 'x_save_file version: ', save_file_version
     DO i = 1, 6 
        READ(10,*)      
     ENDDO
  ENDIF
  WRITE(stdout,*)

  READ(10,*) lsda, nspin
  WRITE(stdout,'(5x,a,i2)') 'nspin:',nspin
  READ(10,*) xm_r, nkstot_r, xnitermax
  WRITE(stdout,'(5x,a,i4)') 'number of k-points:',nkstot
  WRITE(stdout,*)
  WRITE(stdout,'(5x,a,i4)') 'final-state angular momentum (xm_r): ',xm_r
  IF (xm_r==1) THEN
     WRITE(stdout,'(5x,a)') ' => electric-dipole approximation'
  ELSEIF (xm_r==2) THEN
     WRITE(stdout,'(5x,a)') ' => electric-quadrupole approximation'
  ELSE
     WRITE(stdout,'(5x,a)') 'Wrong value of xm_r: STOP'
     CALL stop_xspectra
  ENDIF
  IF(xm_r.NE.xang_mom) & 
     CALL errore('read_save_file','xm_r is different from xang_mom=',xang_mom)

  READ(10,*) ncalcv_max
  IF(ncalcv_max.GT.xnitermax) THEN
     WRITE(stdout,'(5x,a,i5)') 'ncalcv_max=',ncalcv_max
     CALL errore('read_save_file','ncalcv_max is grater than xnitermax=', &
                 xnitermax)
  ENDIF

  WRITE(stdout,*)
  IF (save_file_version < 2) THEN
     READ(10,*) core_energy
  ELSE
     READ(10,*) core_energy, ef
     WRITE(stdout,'(5x,a,f9.4)') 'Fermi level [eV]:', ef 
  ENDIF
  WRITE(stdout,'(5x,a,f10.3,/)') 'core energy [eV]:', core_energy

  READ(10,*) (xkvec_r(i),i=1,3)
  !WRITE(stdout,*) '---------------------------------------------------------'
  !WRITE(stdout,*) 'xkvec read from savefile'
  !WRITE(stdout,*) (xkvec_r(i),i=1,3)
  READ(10,*) (xepsilon_r(i),i=1,3)
  !WRITE(stdout,*) 'xepsilon read from file'
  !WRITE(stdout,*) (xepsilon_r(i),i=1,3)
  !WRITE(stdout,*) '---------------------------------------------------------'
  !write(stdout,*) 'n_lanczos=', n_lanczos
  WRITE(stdout,'(5x,a,1x,3(f10.6,1x))') 'xepsilon [Cartesian frame]:',&
       (xepsilon_r(i),i=1,3)
  IF (xm_r==2) WRITE(stdout,'(5x,a,1x,3(f10.6,1x))') &
       'xkvec [Cartesian frame]:', (xkvec_r(i),i=1,3)
  WRITE(stdout,*)

  DO i = 1, n_lanczos
     IF (TRIM(save_file_kind).eq.'unfinished') THEN
       READ(10,*) (calculated_all(i,j),j=1,nkstot)
     ENDIF
     READ(10,*) (xnorm_all(i,j),j=1,nkstot)
     READ(10,*) (ncalcv_all(i,j),j=1,nkstot)
     READ(10,*) ((a_all(j,k),j=1,ncalcv_max),k=1,nkstot)
     READ(10,*) ((b_all(j,k),j=1,ncalcv_max),k=1,nkstot)
#ifdef __MPI
     CALL poolscatter(xnitermax,nkstot,a_all,nks,aux)
     a(1:xnitermax,i,1:nks)=aux(1:xnitermax,1:nks)
     CALL poolscatter(xnitermax,nkstot,b_all,nks,aux)
     b(1:xnitermax,i,1:nks)=aux(1:xnitermax,1:nks)
#else
     a(1:xnitermax,i,1:nkstot)=a_all(1:ncalcv_max,1:nkstot)
     b(1:xnitermax,i,1:nkstot)=b_all(1:ncalcv_max,1:nkstot)
#endif

  ENDDO
  CLOSE(10)

#ifdef __MPI 
  CALL poolscatter(n_lanczos,nkstot,xnorm_all,nks,xnorm)
  CALL ipoolscatter(n_lanczos,nkstot,ncalcv_all,nks,ncalcv)
  CALL ipoolscatter(n_lanczos,nkstot,calculated_all,nks,calculated)
#else
  IF(nks.NE.nkstot) THEN
     CALL errore('read_save_file','nks\=nkstot',1)
  ENDIF
  xnorm(1:n_lanczos,1:nkstot) = xnorm_all(1:n_lanczos,1:nkstot)
  ncalcv(1:n_lanczos,1:nkstot) = ncalcv_all(1:n_lanczos,1:nkstot)
  calculated(1:n_lanczos,1:nkstot) = calculated_all(1:n_lanczos,1:nkstot)
#endif
  DEALLOCATE(a_all)
  DEALLOCATE(b_all)
  DEALLOCATE(xnorm_all)
  DEALLOCATE(ncalcv_all)
  DEALLOCATE(aux)

END SUBROUTINE read_save_file

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
SUBROUTINE read_header_save_file(x_save_file)
  !----------------------------------------------------------------------------
  USE kinds,    ONLY: DP
  USE klist,    ONLY: nkstot
  USE lsda_mod, ONLY: nspin,lsda
  USE xspectra, ONLY: save_file_version, save_file_kind, n_lanczos
  USE io_global,ONLY: stdout

  IMPLICIT NONE
  CHARACTER(LEN=256), INTENT (IN) :: x_save_file
  INTEGER    :: ierr, nkstot_r
  INTEGER    :: xm_r
  CHARACTER  :: c

  OPEN ( UNIT = 10, FILE = x_save_file, FORM = 'FORMATTED', &
       STATUS = 'UNKNOWN', IOSTAT = ierr )
  CALL errore( 'iosys', 'x_save_file ' // TRIM( x_save_file ) // &
       & ' not found' , ierr )
  REWIND(10)

  READ(10, '(a1)') c
  REWIND(10)
  IF (c == '#') then
     READ(10, '(20x,i8)') save_file_version
     READ(10, '(20x,a32)') save_file_kind
     READ(10,*) 
     READ(10,'(27x,i4)') n_lanczos
     READ(10,*) 
     READ(10,*) 
  ELSE
     save_file_version=0
     save_file_kind='xanes_old'
     n_lanczos=1
  ENDIF

  READ(10,*) lsda,nspin
  READ(10,*) xm_r,nkstot_r
  nkstot=nkstot_r
  CLOSE(10)

  RETURN 
END SUBROUTINE read_header_save_file

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
SUBROUTINE write_save_file(a,b,xnorm,ncalcv,x_save_file)
  !----------------------------------------------------------------------------
  USE kinds,      ONLY: DP
  USE constants,  ONLY: rytoev
  USE klist,      ONLY: nks, nkstot
  USE xspectra,   ONLY: xnitermax, xang_mom, xkvec, xepsilon, xiabs, &
                        save_file_version, save_file_kind,           &
                        n_lanczos, calculated, edge
  USE ener,       ONLY: ef
  USE io_global,  ONLY: ionode
  !*apsi  USE uspp_param, ONLY : psd
  USE lsda_mod,   ONLY: nspin,lsda
  USE uspp_param, ONLY: upf
  USE edge_energy, ONLY: getE
  !
  IMPLICIT NONE
  ! Arguments
  REAL(dp), INTENT (IN) :: a(xnitermax,n_lanczos,nks)
  REAL(dp), INTENT (IN) :: b(xnitermax,n_lanczos,nks)     
  REAL(dp), INTENT (IN) :: xnorm(n_lanczos,nks)
  INTEGER,  INTENT (IN) :: ncalcv(n_lanczos,nks)
  CHARACTER(LEN=256), INTENT (IN) :: x_save_file ! could be removed, 
                                                 ! exists in module
  ! Local variables
  INTEGER :: ierr
  INTEGER :: i, j, k 
  INTEGER :: ncalcv_max
  INTEGER :: calculated_all(n_lanczos,nkstot)
  INTEGER, ALLOCATABLE :: ncalcv_all(:,:)
  REAL(dp), ALLOCATABLE :: a_all(:,:), b_all(:,:), xnorm_all(:,:)
  CHARACTER(LEN=8) :: dte

  IF (ionode)  CALL DATE_AND_TIME(date=dte)

  ALLOCATE(a_all(xnitermax,nkstot))
  ALLOCATE(b_all(xnitermax,nkstot))
  ALLOCATE(xnorm_all(n_lanczos,nkstot))
  ALLOCATE(ncalcv_all(n_lanczos,nkstot))

  ncalcv_all(:,:) = 0
  xnorm_all(:,:)  = 0.d0

  ncalcv_all(1:n_lanczos,1:nks) = ncalcv(1:n_lanczos,1:nks)
  xnorm_all(1:n_lanczos,1:nks)  = xnorm(1:n_lanczos,1:nks)
  calculated_all(1:n_lanczos,1:nks) = calculated(1:n_lanczos,1:nks)

#ifdef __MPI
  CALL poolrecover(xnorm_all,n_lanczos,nkstot,nks)
  CALL ipoolrecover(ncalcv_all,n_lanczos,nkstot,nks)
  CALL ipoolrecover(calculated_all,n_lanczos,nkstot,nks)
#endif

  ncalcv_max = 0
  DO j = 1, n_lanczos
     DO i = 1, nkstot       
        IF (ncalcv_all(j,i).GT.ncalcv_max) ncalcv_max = ncalcv_all(j,i)
     ENDDO
  ENDDO

  IF ( ionode ) THEN
     OPEN ( UNIT = 10, FILE = x_save_file, FORM = 'FORMATTED', &
            STATUS = 'UNKNOWN', IOSTAT = ierr )
     REWIND(10)
     WRITE(10, '(a20,i8)') '# save_file_version=',save_file_version
     WRITE(10, '(a20,a32)') '# save_file_kind   =',save_file_kind
     WRITE(10,'(a7,a8)') '# date=', dte
     WRITE(10,'(a27,i4)') '# number of lanczos stored=',n_lanczos
     WRITE(10,'(a1)') '#'
     WRITE(10,'(a1)') '#'
     WRITE(10,*) lsda,nspin
     WRITE(10,*) xang_mom,nkstot,xnitermax
     WRITE(10,*) ncalcv_max
!     WRITE(10,*) mygetK(upf(xiabs)%psd), ef ! ef here is in eV
     WRITE(10,*) getE(upf(xiabs)%psd,edge), ef ! ef here is in eV
     WRITE(10,*) (xkvec(i),i=1,3)
     WRITE(10,*) (xepsilon(i),i=1,3)
  ENDIF
  DO i=1, n_lanczos
     a_all(:,:) = 0.d0
     b_all(:,:) = 0.d0
     a_all(1:xnitermax,1:nks) =  a(1:xnitermax,i,1:nks)
     b_all(1:xnitermax,1:nks) =  b(1:xnitermax,i,1:nks)
#ifdef __MPI
     CALL poolrecover(a_all,xnitermax,nkstot,nks)
     CALL poolrecover(b_all,xnitermax,nkstot,nks)
#endif
     IF ( ionode) THEN
        IF (TRIM(save_file_kind).EQ.'unfinished') THEN
          WRITE(10,*) (calculated_all(i,j),j=1,nkstot)
        ENDIF
        WRITE(10,*) (xnorm_all(i,j),j=1,nkstot)
        WRITE(10,*) (ncalcv_all(i,j),j=1,nkstot)
        WRITE(10,*) ((a_all(j,k),j=1,ncalcv_max),k=1,nkstot)
        WRITE(10,*) ((b_all(j,k),j=1,ncalcv_max),k=1,nkstot)
     ENDIF
  ENDDO
  CLOSE(10)

  DEALLOCATE(a_all)
  DEALLOCATE(b_all)
  DEALLOCATE(xnorm_all)
  DEALLOCATE(ncalcv_all)

END SUBROUTINE write_save_file

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
SUBROUTINE write_status_of_the_code
  !----------------------------------------------------------------------------
  USE io_global,       ONLY : stdout
  IMPLICIT NONE
  WRITE(stdout,'(5x,a)')&
   '-------------------------------------------------------------------------'
  WRITE (stdout,'(5x,a)') &
   '                      STATUS OF THE CODE (22/04/2009) '
  WRITE(stdout,'(5x,a)')&
   '-------------------------------------------------------------------------'
  WRITE (stdout,'(5x,a)') 'Working features (22/04/2009)'
  WRITE (stdout,'(5x,a)') '-----------------------------'
  WRITE (stdout,'(5x,a)') &
   '- XANES works both in the electric-dipole and -quadrupole approximation,'
  WRITE (stdout,'(5x,a)') '- Spin polarized works'
  WRITE (stdout,'(5x,a)') '- DFT+U implemented, validated'
  WRITE (stdout,'(5x,a)') '- Ultrasoft pseudo works'
  WRITE (stdout,'(5x,a)') '- Cut occupied states working, improved'
  WRITE (stdout,'(5x,a)') '- Terminator working'
  WRITE (stdout,'(5x,a)') '- Multiprojectors TM+USPP working (MCB,CG)'
  WRITE (stdout,'(5x,a)') '- New save file format, with version numbering'
  WRITE (stdout,'(5x,a)') &
   '- Time limit implemented, with restart, seems to work'
  WRITE (stdout,'(5x,a)') &
   '- DFT+U tested ONLY for non ortho wfc, but implemented'
  WRITE(stdout,*) 
  WRITE (stdout,'(5x,a)') 'TO DO'
  WRITE (stdout,'(5x,a)') '-----'
  WRITE (stdout,'(5x,a)') '- L2,3 edges [OB]'
  WRITE (stdout,'(5x,a)') '- Generalization to all edges [OB]'
  WRITE (stdout,'(5x,a)') '- XMCD [?]'
  WRITE (stdout,'(5x,a)') '- IXS [DC]'
  WRITE (stdout,'(5x,a)') '- EELS [DC]'
  WRITE (stdout,'(5x,a)') '- REXS [DC]'
  WRITE (stdout,'(5x,a)') '- Bethe-Salpeter [?] '
  WRITE (stdout,'(5x,a)') '- RXES [?]' 

END SUBROUTINE write_status_of_the_code

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
SUBROUTINE read_gamma_file
  !----------------------------------------------------------------------------
  USE gamma_variable_mod

  IMPLICIT NONE
  INTEGER :: nl, ierr, i

  OPEN ( UNIT = 21, FILE = gamma_file, FORM = 'formatted',&
         STATUS = 'unknown', IOSTAT = ierr)
  CALL errore('io ', 'gamma file '//TRIM(gamma_file)//' not found', abs (ierr))
  REWIND(21)

  nl = 0

  DO
     READ (21,'(a1)',iostat=ierr)
     IF (ierr.NE.0) EXIT
     nl = nl + 1
  ENDDO
  CLOSE(21)

  gamma_lines = nl
  ALLOCATE(gamma_points(nl,2))

  OPEN ( UNIT = 21, FILE = gamma_file, FORM = 'formatted',&
         STATUS = 'unknown', IOSTAT = ierr)
  REWIND(21)

  DO i=1,nl
     READ(21,*) gamma_points(i,1), gamma_points(i,2)
  ENDDO

  CLOSE(21)

END SUBROUTINE read_gamma_file

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
SUBROUTINE initialize_gamma_tab
  !----------------------------------------------------------------------------
  USE kinds,     ONLY: dp
  USE xspectra,  ONLY: xemin, xemax, xnepoint
  USE io_global, ONLY: stdout
  USE constants, ONLY: rytoev
  USE gamma_variable_mod

  IMPLICIT NONE
  REAL(dp) :: e, x, y, dx
  INTEGER  :: i, j, n

  dx = (xemax-xemin)/dfloat(xnepoint)

  DO n = 1, xnepoint
     x = xemin + (n-1)*dx
     i = 1
     DO j = 1, gamma_lines
        IF(x > gamma_points(j,1)) i = i + 1
     ENDDO

     IF (i == 1) THEN
        y = gamma_points(1,2)
     ELSEIF (i == (gamma_lines+1)) THEN
        y = gamma_points(gamma_lines,2)
     ELSE
        y = (  gamma_points(i-1,2) * (gamma_points(i,1)-x)     &
             + gamma_points(i,2)   * (x-gamma_points(i-1,1))  )&
            / ( gamma_points(i,1) - gamma_points(i-1,1) )
     ENDIF
     gamma_tab(n) = y/rytoev
  ENDDO

END SUBROUTINE initialize_gamma_tab

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
SUBROUTINE determine_polycut(t1,t2,f1,f2,df1,df2,poly)
  !----------------------------------------------------------------------------
  ! Calculates the interpolation polynome between 2 points - CG
  !----------------------------------------------------------------------------
  USE kinds, ONLY: dp

  IMPLICIT NONE
  REAL(dp) :: t1, t2, f1, f2, df1, df2, poly(4)

  poly(4) = ( (t2-t1) * (df2+df1) - 2*(f2-f1) )/( (t2-t1)**3 )
  poly(3) = (df2-df1)/(2*(t2-t1)) - &
            1.5d0*(t2+t1)*( (t2-t1)*(df2+df1) - 2*(f2-f1) )/( (t2-t1)**3 )
  poly(2) = df1 - 2*t1*poly(3) - 3*t1*t1*poly(4)
  poly(1) = f1 - poly(2)*t1 - poly(3)*t1**2 - poly(4)*t1**3

END SUBROUTINE determine_polycut
