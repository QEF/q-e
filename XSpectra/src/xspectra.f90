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
  USE uspp_param,           ONLY : upf
  USE xspectra
  USE ener,            ONLY : ef, ef_up, ef_dw !Fermi energy
  USE symm_base,       ONLY : nsym,s
  USE paw_gipaw,              ONLY : read_recon,  &
       paw_vkb,             & ! |p> projectors
       paw_becp,            & ! product of projectors and wf.
       paw_nkb,             & ! total number of beta functions, with st.fact.
       paw_lmaxkb,          &
       paw_recon,           &
       set_paw_upf
  USE klist,       ONLY : &
       nkstot,            & ! total number of k-points
       nks,               & ! number of k-points for local pool
       nelec,nelup,neldw,             & !number of electrons
       xk,                & ! k-points coordinates
       wk ,               & ! k-points weight
       npk,               &
       degauss,lgauss,ngauss,    &
       two_fermi_energies
  USE lsda_mod,    ONLY : nspin,lsda,isk,current_spin
  USE noncollin_module,     ONLY : noncolin
  USE mp,         ONLY : mp_bcast, mp_sum             !parallelization
  USE mp_global,  ONLY : mp_startup, mp_global_end
  USE mp_pools,   ONLY : intra_pool_comm, npool
  USE mp_world,   ONLY : nproc, world_comm
  USE control_flags, ONLY : gamma_only
  USE environment,   ONLY : environment_start

  USE cut_valence_green, ONLY :&
       cut_ierror, &    ! convergence tolerance for one step in the integral
       cut_stepu , &    ! integration initial step, upper side
       cut_stepl , &    ! integration initial step, lower side
       cut_startt, &    ! integration start value of the t variable
       cut_tinf  , &    ! maximum value of the lower integration boundary
       cut_tsup  , &    ! minimum value of the upper integration boudary
       cut_desmooth,&   ! size of the interval near the fermi energy in which cross section is smoothed
       cut_nmemu,&      ! size of the memory of the values of the green function, upper side
       cut_nmeml,&      ! size of the memory of the values of the green function, lower side
       cut_occ_states  ! true if you want tou remove occupied states from the spectrum

  USE control_flags, ONLY : twfcollect
  !<CG>
  USE gamma_variable_mod, ONLY : gamma_lines, gamma_tab, gamma_points, gamma_mode, gamma_file
  USE xspectra_paw_variables, ONLY : xspectra_paw_nhm, init_xspectra_paw_nhm
  !</CG>

  IMPLICIT NONE 
  !
  ! ... local variables
  !
  LOGICAL terminator, show_status, wf_collect
  INTEGER :: nargs,iiarg,ierr,ios,il,ibnd,ibnd_up,ibnd_dw,xm_r,nc_r,ncomp_max
  INTEGER :: iargc
  INTEGER nt,nb,na,i,j,k,nrest,nline
  INTEGER, ALLOCATABLE :: ncalcv(:,:)
  INTEGER, ALLOCATABLE ::&
       paw_iltonhb(:,:,:)      ! corresp l, projector, type <--> cumulative over all the species
  REAL (DP) ehomo, elumo,norm, core_energy
  REAL(KIND=DP) :: core_wfn(ndmx)
  REAL(dp), ALLOCATABLE:: a(:,:,:),b(:,:,:),xnorm(:,:)      !lanczos vectors
  REAL (DP), EXTERNAL :: efermig,efermit
  REAL(DP) :: rc(ntypx,0:lmaxx),r_paw(0:lmaxx)
  LOGICAL :: exst, loc_set

  CHARACTER (LEN=256)  :: input_file, filerecon(ntypx),filecore
  CHARACTER(LEN=256) :: outdir 
  CHARACTER(LEN=25) :: calculation
  CHARACTER(LEN=4) :: verbosity
  CHARACTER(LEN=10) :: dummy_char
  REAL(dp) :: gamma_energy(2), gamma_value(2)
  REAL(dp) :: auxrpol(3,2)

  ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  ! $   Namelists Definition
  ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$


  namelist / input_xspectra / &
       calculation,&       ! calulation type : xanes...
       verbosity, &        ! high/low
       prefix, &           ! prefix of the pwscf output files
       outdir, &           ! directory tmp_dir or where the files are
       xiabs,&
       xkvec,&
       xepsilon,&
       xcoordcrys,&
       ef_r,&
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
       restart_mode   

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


  ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  !    initialising MPI environment, clocks, a few other things
  ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

#ifdef __MPI
  CALL mp_startup ( )
#endif
  CALL environment_start ( 'XSPECTRA' )

  ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  ! $   Default values for namelists
  ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$


  !
  ! Set default values for namelist:
  !

  calculation='xanes'
  prefix=' '
  verbosity='low'
  x_save_file='xanes.sav'
  CALL get_env( 'ESPRESSO_TMPDIR', outdir )
  IF ( TRIM( outdir ) == ' ' ) outdir = './'
  xnepoint=100
  xniter=50
  xcheck_conv=5
  xonly_plot=.FALSE.
  xread_wf=.FALSE.
  xemin=0.d0
  xemax=10.d0
  xgamma=0.1d0
  xerror=1.d-2
  xiabs=1               !identify the adsorbing atom
  DO i=2,3
     xkvec(i)=0.d0
     xepsilon(i)=0.d0
  ENDDO
  xkvec(1)=1.d0
  xepsilon(1)=1.d0
  xcoordcrys=.true.
  ef_r=0.d0
  cut_occ_states=.FALSE.
  terminator=.false.
  show_status=.false.
  wf_collect=.false.
  gamma_mode='constant'
  gamma_file='gamma.dat'
  U_projection_type='atomic'

  ! default values for cutting the occupied states (paste_fermi function)
  cut_ierror=1.d-7
  cut_stepu=1.d-2
  cut_stepl=1.d-3
  cut_startt=1.d0
  cut_tinf=1.d-6
  cut_tsup=100.d0
  cut_desmooth=1.d-2
  cut_nmemu=100000
  cut_nmeml=100000

  ! Set default values for other variables

  filecore='Core.wfc'

  restart_mode='from_scratch'
  time_limit=1.d8


  ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  ! $   Check if the input is from file or from stdin
  ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  !
  ! ... Input from file ?
  !
  IF ( ionode ) THEN

     nargs = iargc()
!     nargs = command_argument_count() ! if iargc does not work

     !
     DO iiarg = 1, ( nargs - 1 )
        !
        CALL getarg( iiarg, input_file )
        IF ( TRIM( input_file ) == '-input' .OR. &
             TRIM( input_file ) == '-inp'   .OR. &
             TRIM( input_file ) == '-in' ) THEN
           !
           CALL getarg( ( iiarg + 1 ) , input_file )
           OPEN ( UNIT = 5, FILE = input_file, FORM = 'FORMATTED', &
                STATUS = 'OLD', IOSTAT = ierr )
           CALL errore( 'iosys', 'input file ' // TRIM( input_file ) // &
                & ' not found' , ierr )
           !
        END IF
        !
     END DO



     ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
     ! $   Reading namelists
     ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$


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


  ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  ! $   Variables broadcasting 
  ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  CALL mp_bcast( calculation, ionode_id, world_comm )
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

! restart
  CALL mp_bcast( time_limit, ionode_id, world_comm )
  CALL mp_bcast( restart_mode, ionode_id, world_comm )


  ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  ! $ Writing the status of the code (working features and things to do)
  ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  IF(show_status) CALL WRITE_status_of_the_code

  ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  !    Initialising clocks
  ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  !

  IF(TRIM(ADJUSTL(calculation)).EQ.'xanes_dipole') THEN
     n_lanczos=1
     xang_mom=1                       !so it is not necessary to specify xang_mom
     calculation='xanes'
  ELSEIF(TRIM(ADJUSTL(calculation)).EQ.'xanes_quadrupole') THEN
     n_lanczos=1
     xang_mom=2                       !so it is not necessary to specify xang_mom
     calculation='xanes'
  ENDIF

  CALL start_clock( calculation  )
  !CALL stop_clock( code )


  ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  !  check on wfcollect
  ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  IF(xread_wf.AND.wf_collect) THEN
     CALL errore ('main','incompatibility xread_wf and wf_collect',1)
  ENDIF


  twfcollect=wf_collect


  ! $$$$$$$$$  notstart if on  xonlyplot  $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  IF(.NOT.xonly_plot) THEN

     ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
     !    read pwscf structural and k-points infos, also ditributes across the pools
     ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

     CALL read_file()

     ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
     !    initialize everything as in a nscf calculation
     ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

     call reset_k_points_and_reinit_nscf()

     IF(xread_wf) THEN
        WRITE(stdout,*) ' '
        IF (okvan) THEN
           WRITE(stdout,*) 'Approx. ram memory needed per proc in MB = ',(16.0*4.0*npwx*npool)/(nproc*1000000.0)
        ELSE
           WRITE(stdout,*) 'Approx. ram memory needed per proc in MB = ', (16.0*3.0*npwx*npool)/(nproc*1000000.0)
        ENDIF
        WRITE(stdout,*)
     ENDIF
     WRITE(stdout,*) 'k-points : nkstot=', nkstot


     ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$  
     ! normalize xkvec and xepsilon
     ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$  

     IF(xcoordcrys) CALL cryst_to_cart(1,xepsilon,at,1)
     IF(xang_mom.EQ.2) THEN
        IF(xcoordcrys) CALL cryst_to_cart(1,xkvec,at,1)
        norm=DSQRT(xkvec(1)**2+xkvec(2)**2+xkvec(3)**2)
        DO i=1,3
           xkvec(i)=xkvec(i)/norm
        ENDDO
     ENDIF
     norm=DSQRT(xepsilon(1)**2+xepsilon(2)**2+xepsilon(3)**2)
     DO i=1,3
        xepsilon(i)=xepsilon(i)/norm
     ENDDO

     ! check orthogonality

     WRITE (stdout,*) '--- Polarisation and k vector [cartesian coordinates]----'
     WRITE (stdout,'(a,1x,3(f10.8, 1x))') 'xepsilon(:)=', (xepsilon(i),i=1,3)
     WRITE (stdout,'(a,1x,3(f10.8, 1x))') 'xkvec(:)=', (xkvec(i),i=1,3)
     IF(xang_mom.EQ.2) THEN
        IF ((abs(xkvec(1)*xepsilon(1)+xkvec(2)*xepsilon(2)+xkvec(3)*xepsilon(3))).ge.1.0d-6) THEN
           WRITE(stdout,*) 'WARNING, xkvec and xepsilon are not orthogonal'
           WRITE(stdout,*) 'scalar product=',xkvec(1)*xepsilon(1)+xkvec(2)*xepsilon(2)+xkvec(3)*xepsilon(3)
        ENDIF
     ENDIF



     ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
     !  is the type associated to xiabs existing ?
     ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$  

     
     i=0
     DO na=1,nat
        IF(ityp(na).EQ.xiabs) i=i+1
     ENDDO
     IF(i.NE.1) THEN
        CALL errore( 'main program', 'Wrong xiabs!!!',i)
     ENDIF

     ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
     !  Reads reconstruction files
     ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$  

     DO nt = 1, ntyp
        call set_paw_upf(nt, upf(nt))
     ENDDO


     CALL read_core_abs(filecore,core_wfn)
     

     IF ( .NOT. paw_recon(xiabs)%gipaw_data_in_upf_file ) &
          CALL read_recon ( filerecon(xiabs), xiabs, paw_recon(xiabs) ) !*apsi


     ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
     !  Assign paw radii to species (this will become soon obsolete)
     ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$  
     !
     !  read recon should be parallelized
     !


     DO nt=1,ntyp
        IF ((.NOT.paw_recon(nt)%gipaw_data_in_upf_file)) then
           paw_recon(nt)%paw_nbeta=0
        ENDIF
        DO  j=1,paw_recon(nt)%paw_nbeta
           il=paw_recon(nt)%psphi(j)%label%l
           IF(xiabs.EQ.nt.AND.DABS(r_paw(il)).lt.1.d-6) THEN
              !*apsi              if(r(kkbeta(nt),nt).GT.1.d-3) THEN
              !*apsi                rc(nt,il)=  r(kkbeta(nt),nt)
              !*apsi              ELSE
              !*apsi                 WRITE(stdout,*) 'Warning-no-kkbeta r_paw(',il,')=1.0'
              !*apsi                 rc(nt,il)=1.0
              !*apsi              ENDIF
              !<CG>  to be verified
              IF (paw_recon(nt)%psphi(j)%label%rc > 1.d-3) THEN
                 WRITE(stdout,*) 'warning, r_paw(', il,' ) set to ', &
                      paw_recon(nt)%psphi(j)%label%rc
                 rc(nt, il)= paw_recon(nt)%psphi(j)%label%rc*3.0/2.0
              ELSE
                 WRITE(stdout,*) 'Warning, no rc'
                 WRITE(stdout,*) 'warning, r_paw(', il,' ) set to 1.5'
                 rc(nt, il)= 1.5d0
              ENDIF
              !</CG>

           ELSEIF(xiabs.EQ.nt.AND.DABS(r_paw(il)).GT.1.d-6) THEN
              rc(nt,il)=r_paw(il)
           ELSEIF(nt.NE.xiabs) THEN
              !*apsi              IF(r(kkbeta(nt),nt).GT.1.d-3) THEN
              !*apsi                 rc(nt,il)=r(kkbeta(nt),nt)
              !*apsi              ELSE
              !*apsi                 rc(nt,il)=1.0
              !*apsi              ENDIF
              !<CG> to be verified
              IF(paw_recon(nt)%psphi(j)%label%rc.GT.1.d-3) THEN
                 rc(nt,il)=paw_recon(nt)%psphi(j)%label%rc*3.0/2.0
              ELSE
                 rc(nt,il)=1.5
              ENDIF
              !<CG>
           ENDIF
        ENDDO
     ENDDO

     !<CG>
     DO nt=1,ntyp
        DO j = 1,paw_recon(nt)%paw_nbeta
           paw_recon(nt)%psphi(j)%label%rc = rc(nt,paw_recon(nt)%psphi(j)%label%l)
           paw_recon(nt)%aephi(j)%label%rc = rc(nt,paw_recon(nt)%aephi(j)%label%l)
        ENDDO
     ENDDO
     !</CG>


     ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
     !  write band energies if xread_wf=true
     ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$  


     IF(xread_wf.AND.TRIM(verbosity).EQ.'high') THEN

        WRITE(stdout,*) '------- Band energies read from file -----------'

        DO i=1,nkstot
           WRITE(stdout,'("k=[",3f14.8,"]   spin=",1i2)') &
                xk(1,i),xk(2,i),xk(3,i),isk(i)
           WRITE(stdout, '(8f9.4)') (et(j,i)*rytoev,j=1,nbnd)
        ENDDO


        WRITE(stdout,*) '------------------------------------------------'
     ENDIF


     ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
     !  Calculate the Fermi level if possible
     ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$  

     WRITE(stdout,*)
     WRITE(stdout,*) '=========================================================='

     IF(TRIM(calculation).EQ.'fermi_level'.AND..NOT.xread_wf) THEN
        WRITE(stdout,*) 'Impossible to calculate the Fermi level'
        WRITE(stdout,*) '   without reading wavefunctions (xread_wf=.true.)'
        CALL errore( 'main program', &
             'xread_wf incompatible with calculation type', 1 )
     ENDIF


     IF(.NOT.xread_wf) THEN        !Fermi level must be read from input
        ef=ef_r
        WRITE( stdout,*)
        WRITE( stdout,'(">> Fermi level from input (Ry)=",1f14.8)') ef
     ENDIF

     IF(TRIM(calculation).EQ.'fermi_level') THEN
        IF(ltetra.OR.lgauss) THEN
           WRITE( stdout,*)
           WRITE( stdout,*) 'Metallic case'
           WRITE( stdout,*)
           CALL errore('input','Read fermi level from scf/nscf output',1)
        ELSE
           WRITE(stdout,*)
           WRITE(stdout,*)  'Insulating case:'
           WRITE(stdout,*)
           !   SPINLESS CASE
           IF ( nspin.EQ.1) THEN
              ibnd =  NINT (nelec) / 2
              ef = MAXVAL ( et( ibnd  , 1:nkstot) )
              WRITE( stdout, '(">> Fermi level (Ry) = ",1f14.8)') ef
              !SPIN POLARIZED CALCULATION
           ELSEIF(nspin.EQ.2) THEN
              ibnd    = NINT( nelec )
              IF(nelup.EQ.0.OR.nelup.EQ.0) THEN
                 WRITE(stdout,*) 'WARNING, nelup=0 or neldw=0'
              ENDIF
              WRITE(stdout,*) 'nel=',nelec,' nelup=',nelup, 'neldw=',neldw

              ibnd_up = NINT( nelup )
              ibnd_dw = NINT( neldw )
              IF ( ibnd_up == 0 ) THEN
                 !
                 ef = MAXVAL( et(ibnd_dw,1:nkstot/2) )
                 ef_dw = ef
                 WRITE( stdout,&
                      '(">> Fermi level down (Ry)= ",1f14.8)')&
                      ef_dw
                 !
              ELSE IF ( ibnd_dw == 0 ) THEN
                 !
                 ef = MAXVAL( et(ibnd_up,1:nkstot/2) )
                 ef_up = ef
                 WRITE( stdout,&
                      '(">> Fermi level up (Ry)= ",1f14.8)')&
                      ef_up
                 !
              ELSE
                 !
                 WRITE(stdout,*) 'nkstot=',nkstot
                 ef    = MAX( MAXVAL( et(ibnd_up,1:nkstot/2) ), &
                      MAXVAL( et(ibnd_dw,nkstot/2+1:nkstot) ) )
                 ef_up =  MAXVAL( et(ibnd_up,1:nkstot/2) )
                 ef_dw =  MAXVAL( et(ibnd_dw,nkstot/2+1:nkstot) )
                 WRITE( stdout,&
                      '(">> Fermi level up (Ry)= ",1f14.8,&
                      &" Fermi level down (Ry)= ",1f14.8)') &
                      ef_up,ef_dw
                 !
              END IF
              WRITE( stdout,'(">> Fermi level (Ry)=",1f14.8)') ef
           END IF

        ENDIF
        WRITE( stdout,*)
        WRITE (stdout,*) &
             '============================================================'

        call stop_xspectra () 

     ENDIF
     !        CALL mp_bcast( ef, ionode_id )  !Why should I need this ?


     ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
     !  Allocation of variables for paw projectors
     ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$  

     ! CG : becp allocated in the xanes_dipole and xanes_quadrupole subroutines
     !     call allocate_bec_type ( nkb, nbnd, becp )
     ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
     !  Initialise Vanderbilt and Paw projectors
     ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$  
     !-----------

     !     CALL init_us_1  ! CG

     !<CG>     
     CALL init_gipaw_1
     !</CG>     

     ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
     !  Definition of a particular indexation to avoid Mickael Profeta's crazy indexation
     ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$  
     ! <CG>
     CALL init_xspectra_paw_nhm
     !<\CG>
     ALLOCATE (paw_iltonhb(0:paw_lmaxkb,xspectra_paw_nhm, ntyp))
     CALL define_index_arrays(paw_iltonhb)

     ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
     !  Allocate paw projectors
     ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$  

     ALLOCATE (paw_vkb( npwx,  paw_nkb))
     ALLOCATE (paw_becp(paw_nkb, nbnd))

  ELSEIF(xonly_plot) THEN  !$$$$$$$$$$$$$$$$$$  xonly_plot if structure
     ! Fermi level read from file
     ef=ef_r
     WRITE( stdout,*)
     WRITE( stdout,'(">> Fermi level read from file (Ry)=",1f14.8)') ef
     WRITE( stdout,*)
  ENDIF   !$$$$$$$$$$$$$$$$$$  xonly_plot if structure

  !<CG>
  ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  ! $ Computing gamma tabulated values
  ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  !<MCB> THIS MUST BE CHANGED

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

  ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  !</CG>

  xnitermax=xniter



  IF(xonly_plot) THEN
     CALL read_header_save_file(x_save_file)
     nks = nkstot
     WRITE(6,*) 'nks=',nks
     IF(lsda) THEN
        isk(1:nkstot/2)=1
        isk(nkstot/2+1:nkstot)=2
        wk(1:nkstot)=2.d0/nkstot
     ELSEIF(.NOT.lsda) THEN
        isk(1:nkstot)=1
        wk(1:nkstot)=2.d0/nkstot
     ENDIF
     CALL divide_et_impera( xk, wk, isk, lsda, nkstot, nks )
  ENDIF

  ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  !  Verification of paw relations between pseudo partial waves and projector (radial parts)
  ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$  

  IF(.NOT.xonly_plot.AND.TRIM(verbosity).EQ.'high')  CALL check_paw_projectors(xiabs)


  ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  !  Allocate and initialise lanczosvectors
  ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$  

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

  ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  !  And now we go...  XANES CALCULATION
  ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$  

  !  IF (okvan) THEN
  !    WRITE(stdout,*) 'Ultrasoft not implemented'
  !    CALL stop_pp
  !  ENDIF


  IF(TRIM(calculation).EQ.'xanes') THEN
     IF(.NOT.xonly_plot) THEN
        IF(TRIM(ADJUSTL(restart_mode)).eq.'restart') THEN
          CALL read_header_save_file(x_save_file)
          CALL read_save_file(a,b,xnorm,ncalcv,x_save_file,core_energy)
        ENDIF
        IF(xang_mom.EQ.1) THEN
           save_file_version=1
           save_file_kind='xanes_dipole'

           CALL xanes_dipole(a,b,ncalcv,xnorm,core_wfn,paw_iltonhb,terminator,verbosity)
           ! open save file and write everything
           CALL write_save_file(a,b,xnorm,ncalcv,x_save_file)
        ELSEIF(xang_mom.EQ.2) THEN
           save_file_version=1
           save_file_kind='xanes_quadrupole'
           CALL xanes_quadrupole(a,b,ncalcv,xnorm,core_wfn,paw_iltonhb,terminator,verbosity)
           CALL write_save_file(a,b,xnorm,ncalcv,x_save_file)
        ENDIF
     ELSE
        CALL read_save_file(a,b,xnorm,ncalcv,x_save_file,core_energy)
     ENDIF

     IF (TRIM(save_file_kind).eq.'unfinished') CALL stop_xspectra ()       

     IF(xang_mom.EQ.1) THEN
        CALL plot_xanes_dipole(a,b,xnorm,ncalcv,terminator,core_energy)
     ELSEIF(xang_mom.EQ.2) THEN
        CALL plot_xanes_quadrupole(a,b,xnorm,ncalcv,terminator,core_energy)
     ENDIF

  ELSEIF(TRIM(calculation).EQ.'rxes') THEN
     CALL errore( 'Main', 'rxes Not yet implemented',1)
  ELSEIF(TRIM(calculation).EQ.'bethe_salpeter') THEN
     CALL errore( 'Main', 'bethe_salpeter Not yet implemented',1)
  ENDIF

  ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  !  Deallocation
  ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$  

  !  IF(.NOT.xonly_plot) THEN
  !     call deallocate_bec_type ( becp )
  !  ENDIF

  DEALLOCATE(a)
  DEALLOCATE(b)
  DEALLOCATE(xnorm)
  DEALLOCATE(ncalcv)


  WRITE (stdout,*) 'End program ', TRIM(calculation)

    CALL stop_clock( calculation  )
    CALL print_clock( calculation )

  CALL stop_xspectra () 

END program X_Spectra

!--------------------------------------------------------------------
SUBROUTINE stop_xspectra
  !--------------------------------------------------------------------
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


! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!  Dipolar Calculation
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$  


SUBROUTINE xanes_dipole(a,b,ncalcv,xnorm,core_wfn,paw_iltonhb,terminator,verbosity)
  USE constants,       ONLY : fpi
  USE io_global,       ONLY : stdout     ! Modules/io_global.f90
  USE kinds,           ONLY : DP
  USE parameters,      ONLY : ntypx
  USE radial_grids,    ONLY : ndmx
  USE ions_base,       ONLY : nat, ntyp => nsp, ityp
  USE wvfct,           ONLY : npwx, nbndx, nbnd, npw, igk, g2kin, et,&
                              current_k, ecutwfc
  USE symm_base,       ONLY : d1,d2,d3
  USE noncollin_module,     ONLY : noncolin
  USE lsda_mod,        ONLY : nspin,lsda,isk,current_spin
  USE cell_base,       ONLY: tpiba2, bg
  USE wavefunctions_module, ONLY: evc
  USE klist,           ONLY : &
       nkstot,            & ! total number of k-points
       nks,               & ! number of k-points per pool
       xk,                & ! k-points coordinates
       wk                   ! k-points weight
  USE gvect,            ONLY: g, ngm, ngl
  USE fft_base,         ONLY: dfftp
  USE paw_gipaw,        ONLY : &
       paw_vkb,             & ! |p> projectors
       paw_becp,            & ! product of projectors and wf.
       paw_nkb,             & ! total number of beta functions, with st.fact.
       paw_lmaxkb,paw_recon
  USE becmod,          ONLY : becp, allocate_bec_type, deallocate_bec_type !CG
  USE scf,             ONLY : vltot, vrs, v, kedtau
  USE gvecs,           ONLY : doublegrid
  USE mp_world,        ONLY : world_comm
  USE mp_pools,        ONLY : intra_pool_comm, root_pool
  USE mp,              ONLY : mp_sum, mp_bcast, mp_barrier !CG
  USE io_global,       ONLY : ionode

  USE xspectra,        ONLY : xiabs, xanes_dip, xang_mom, xniter,&
                              xnitermax, xepsilon,time_limit,calculated,&
                            save_file_kind

  USE atom,       ONLY : rgrid, msh
  !  use atom,        ONLY : &
  !       mesh,     &!mesh(ntypx) number of mesh points
  !       msh ,     &!msh(ntypx)the point at rcut=end of radial integration
  !       r   
  USE radin_mod
  USE basis,                ONLY : natomwfc

  USE uspp,   ONLY : vkb, nkb, okvan !CG
  USE uspp_param,       ONLY : upf

  USE ldaU,   ONLY : lda_plus_u, init_lda_plus_u 
  !<CG>
  USE xspectra_paw_variables, ONLY : xspectra_paw_nhm
  !</CG>

  IMPLICIT NONE
  REAL(dp) core_wfn(ndmx)
  REAL(dp) a(xnitermax,1,nks),b(xnitermax,1,nks)     
  REAL (dp)  xnorm(1,nks)
  INTEGER :: is,ik,iabso,nr,ip,jp,l,j,icrd,ip_l,nrc,nt,na
  INTEGER :: ipx,ipx_0,ipy,ipz,nline,nrest,npw_partial
  INTEGER :: ncalcv(1,nks)
  INTEGER :: paw_iltonhb(0:paw_lmaxkb,xspectra_paw_nhm, ntyp)
  REAL (dp) pref,prefb,v_of_0,xnorm_partial
  REAL (dp) norm
  REAL (dp), ALLOCATABLE :: aux(:)
  COMPLEX(KIND=DP), EXTERNAL :: zdotc
  COMPLEX(dp), ALLOCATABLE :: paw_vkb_cplx(:,:)
  LOGICAL :: terminator
  REAL(dp) :: normps
  CHARACTER(LEN=4) :: verbosity

  EXTERNAL zdscal

  LOGICAL :: recalc
  COMPLEX(dp), ALLOCATABLE :: psiwfc(:), spsiwfc(:)
  REAL(DP), EXTERNAL ::  get_clock
  REAL(dp) :: timenow=0 
  INTEGER :: nunfinished=0

  ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  !  Constant Definitions
  ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$  

  pref=SQRT(3.d0/2.d0)
  prefb=SQRT(3.d0)

  ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  !  Variable allocation and initialization
  ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$  

  ALLOCATE(aux(rgrid(xiabs)%mesh)) !allocation too big, it needs only up to msh
  ALLOCATE (paw_vkb_cplx( npwx,  paw_nkb))
  ALLOCATE(xanes_dip(paw_recon(xiabs)%paw_nl(xang_mom)))
  ALLOCATE(psiwfc(npwx))
  xanes_dip(:)=0.d0


  ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  ! Radial Dipole Matrix element
  ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$  

  !  I compute the radial part, 
  !          <\phi_n|r|\psi_i>=int r^3 \phi_n \psi_i dr.
  !  Here only ps and ae wavefunctions from XX.recon
  !

  WRITE( stdout,*) 'Calculation dipole matrix element'
  WRITE( stdout,*) 'There are ',paw_recon(xiabs)%paw_nl(xang_mom),' projectors/channel'
  WRITE( stdout,*) 'for angular moment ',xang_mom,' and atom type ',xiabs

  ! I check that the core wf is correctly normalized

  nr=msh(xiabs)  ! extended up to all the NON ZERO points in the mesh.

 IF(TRIM(verbosity).EQ.'high') THEN
     aux(1:nr)=core_wfn(1:nr)*core_wfn(1:nr)
     WRITE (stdout,'(" norm of core wfc =",1f14.8)') SQRT(para_radin(aux(1:nr),rgrid(xiabs)%r(1:nr),nr))
 ENDIF



  ! and I calculate the radial integral

  ip_l=0

  DO ip=1,paw_recon(xiabs)%paw_nbeta
     !     IF(psphi(xiabs,ip)%label%l.EQ.xang_mom) THEN
     IF(paw_recon(xiabs)%aephi(ip)%label%l.EQ.xang_mom) THEN

        nrc=paw_recon(xiabs)%aephi(ip)%label%nrc
        ip_l=ip_l+1
        !    here below, recall that psi is r*psi and you have a Jacobian=r^2
        aux(1:nrc)=rgrid(xiabs)%r(1:nrc)*paw_recon(xiabs)%aephi(ip)%psi(1:nrc)*core_wfn(1:nrc)
        !    here we have to integrate only inside the augmentation region.
        xanes_dip(ip_l)=para_radin(aux(1:nrc),rgrid(xiabs)%r(1:nrc),nrc)
     ENDIF
  ENDDO


  WRITE(6,*) '----------------------------------------------------------------'
  DO ip=1,paw_recon(xiabs)%paw_nl(xang_mom) !
     WRITE( stdout,'("dipole radial matrix element proj. (",i2,")=",f14.8)') ip,xanes_dip(ip)
  ENDDO
  WRITE(6,*) '----------------------------------------------------------------'
  DEALLOCATE(aux)

  !
  !  We count the projectors for all the atoms
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


  ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  ! Beginning the loop over the k-points
  ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$  

  !*apsi  call set_vrs(vrs,vltot,vr,nrxx,nspin,doublegrid)
  !<CG>
  CALL set_vrs(vrs,vltot,v%of_r,kedtau, v%kin_r,dfftp%nnr,nspin,doublegrid)
  !</CG>


  
  DO ik=1,nks

     IF (calculated(1,ik).EQ.1) CYCLE

     timenow=get_clock( 'xanes' )
     CALL mp_bcast(timenow,root_pool, intra_pool_comm)     

     IF( timenow.GT.time_limit) THEN
       nunfinished=1
       EXIT
     ENDIF
 
     WRITE(stdout,*)
     WRITE(stdout,*) 'Starting k-point : ', ik
     WRITE( stdout,'(" total cpu time spent up to now is ",F9.2," secs")') timenow

     current_k=ik

     IF(lsda) current_spin=isk(ik)

     !gk_sort  sort k-points and exit kinetic energies 
     CALL gk_sort(xk (1,ik),ngm,g,ecutwfc/tpiba2,npw,igk,g2kin)  !CHECK
     g2kin=g2kin*tpiba2                                          !CHECK

     npw_partial = npw
     CALL mp_sum( npw_partial, intra_pool_comm )

     IF(xniter.ge.npw_partial) THEN
        xniter = npw_partial
        WRITE(stdout,*) 'Hilbert space is saturated'
        WRITE(stdout,*) 'xniter is set equal to ',npw_partial
        WRITE(stdout,*) 'Hint: Increase Kinetic Energy cutoff in your SCF simulation'
     ENDIF


     !<CG>        
     CALL init_gipaw_2(npw,igk,xk(1,ik),paw_vkb)
     !</CG>
     if(.not.lda_plus_u) CALL init_us_2(npw,igk,xk(1,ik),vkb)
     IF (lda_plus_u) CALL orthoUwfc_k(ik)


     ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
     ! Angular Matrix element
     ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$  
     !
     ! Here I define human projectors
     !
     ! paw_vkb(1:npw,ipx) are real spherical harmonics. The real spherical harmonics are
     ! defined as
     ! 
     !     y_{l,2m}  =[Y_{l,m}+(-1)^m Y_{l,-m}]/SQRT(2)
     !     y_{l,2m+1}=[Y_{l,m}-(-1)^m Y_{l,-m}]/(i*SQRT(2))
     !
     ! (remember Y_{l,m}=(-1)^m Y_{l,-m)^*  )
     !
     ! The complex spherical harmonics can be written has a function of the real
     ! ones as :
     !
     !     Y_{l,m}  =        [y_{l,2m}+iy_{l,2m+1}]/SQRT(2)     
     !     Y_{l,-m} = (-1)^m [y_{l,2m}-iy_{l,2m+1}]/SQRT(2)     
     !
     !  The paw_vkb_cplx are the Y_{l,m} so the usual spherical harmonics
     !
     ! rotational invariance has been checked


     DO ip=1,paw_recon(xiabs)%paw_nl(xang_mom)

        ipx=ipx_0+paw_iltonhb(xang_mom,ip,xiabs)


        paw_vkb_cplx(1:npw,ipx)= paw_vkb(1:npw,ipx)   !m=0

        paw_vkb_cplx(1:npw,ipx+1)=       &
             (paw_vkb(1:npw,ipx+1)+(0.d0,1.d0)*paw_vkb(1:npw,ipx+2))/SQRT(2.0)    !m=+1

        paw_vkb_cplx(1:npw,ipx+2)=-       &
             (paw_vkb(1:npw,ipx+1)-(0.d0,1.d0)*paw_vkb(1:npw,ipx+2))/SQRT(2.0)    !m=-1



     ENDDO

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


        psiwfc(1:npw)=psiwfc(1:npw)&
             +(                    &
             pref*paw_vkb_cplx(1:npw,ipx+2)*(xepsilon(1)+(0.d0,1.d0)*xepsilon(2))&
             + pref*paw_vkb_cplx(1:npw,ipx+1)*(-xepsilon(1)+(0.d0,1.d0)*xepsilon(2))&
             +prefb*paw_vkb_cplx(1:npw,ipx)*xepsilon(3) &
             )*xanes_dip(ip)/SQRT(fpi)


     ENDDO
     psiwfc(1:npw)=psiwfc(1:npw)*SQRT(fpi)/3.0





     ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
     ! Starting Lanczos
     ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$  


     !   
     !   I normalize the wavefunction psiwfc(1:npw)
     !

     !<CG>
     CALL allocate_bec_type(nkb,1,becp)
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
     WRITE( stdout,*) 'norm initial vector=',xnorm(1,ik)
     norm=1.d0/xnorm(1,ik)

     CALL zdscal(npw,norm,psiwfc,1)
     !
     !      Then I call the lanczos routine
     !
     WRITE(stdout,*) 'Starting lanczos'

     IF (okvan) THEN
        CALL lanczos_uspp(a(:,1,ik),b(:,1,ik),psiwfc,ncalcv(1,ik), terminator)
     ELSE
        CALL lanczos(a(:,1,ik),b(:,1,ik),psiwfc,ncalcv(1,ik), terminator)
     ENDIF


     !
     !      Then I write small report of the lanczos results
     !

     IF(TRIM(verbosity).EQ.'high') THEN
        WRITE( stdout,*) '-----------------------------------------'
        WRITE( stdout,*) 'k-point number =',ik
        WRITE( stdout,*) 'k-point coordinate, isk'
        WRITE( stdout,'(3f12.6,1i2)') xk(1,ik),xk(2,ik),xk(3,ik),isk(ik)
        WRITE( stdout,*) 'Norm of the initial vector =',xnorm(1,ik)
        WRITE( stdout,*) 'Number of iterations =',ncalcv(1,ik)

        !        nline=ncalcv(icrd,ik)/6
        !        nrest=ncalcv(icrd,ik)-nline*6
        !        WRITE( stdout,*) 'a vectors:'
        !        DO ip=1,nline
        !           WRITE( stdout,"(6(f10.6,3x))") (a((ip-1)*6+j,icrd,ik),j=1,6)
        !        ENDDO
        !        WRITE( stdout,"(6(f10.6,3x))") (a(nline*6+j,icrd,ik),j=1,nrest)
        !        WRITE( stdout,*) 'b vectors:'
        !        DO ip=1,nline
        !           WRITE( stdout,"(6(f10.6,3x))") (b((ip-1)*6+j,icrd,ik),j=1,6)
        !        ENDDO
        !        WRITE( stdout,"(6(f10.6,3x))") (b(nline*6+j,icrd,ik),j=1,nrest)
        WRITE( stdout,*) '-----------------------------------------'
     ENDIF

     CALL deallocate_bec_type ( becp ) ! CG
     calculated(1,ik)=1

  ENDDO  !on k points

  ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  ! Array deallocation
  ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$  

  CALL mp_barrier(world_comm)
  CALL mp_sum(nunfinished, world_comm)
  IF (nunfinished >= 1) THEN 
    save_file_kind='unfinished'
    write(stdout,*) 'calculation not finished'
  ENDIF

  DEALLOCATE(psiwfc)
  DEALLOCATE(xanes_dip)
  DEALLOCATE (paw_vkb_cplx)

END SUBROUTINE xanes_dipole

! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!  Quadrupolar Calculation
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$  


SUBROUTINE xanes_quadrupole(a,b,ncalcv,xnorm,core_wfn,paw_iltonhb,terminator,verbosity)
  USE io_global,       ONLY : stdout     ! Modules/io_global.f90
  USE kinds, ONLY : DP
  USE constants, ONLY : pi
  USE parameters,       ONLY : ntypx
  USE radial_grids,     ONLY : ndmx
  USE ions_base,   ONLY : nat, ntyp => nsp, ityp
  USE wvfct,            ONLY : npwx,nbndx,nbnd,npw,igk,&
       g2kin,et, current_k, ecutwfc
  !       ,igk_l2g
  USE lsda_mod,    ONLY : nspin,lsda,isk,current_spin
  USE cell_base, ONLY: tpiba2, bg
  USE wavefunctions_module, ONLY: evc
  USE klist,       ONLY : &
       nkstot,            & ! total number of k-points
       nks,               & ! number of k-points per pool
       xk,                & ! k-points coordinates
       wk                   ! k-points weight
  USE gvect, ONLY: g,ngm,ngl
  USE fft_base,        ONLY : dfftp
  USE paw_gipaw,     ONLY : &
       paw_vkb,             & ! |p> projectors
       paw_becp,            & ! product of projectors and wf.
       paw_nkb,             & ! total number of beta functions, with st.fact.
       paw_lmaxkb, & 
       paw_recon
  USE becmod, ONLY:becp, allocate_bec_type, deallocate_bec_type ! CG
  USE scf, ONLY: vltot,v,vrs, kedtau !CG
  USE gvecs, ONLY : doublegrid
  USE mp_pools,   ONLY : intra_pool_comm, root_pool
  USE mp_world,   ONLY : world_comm
  USE mp,         ONLY : mp_sum,mp_barrier, mp_bcast !CG
  USE xspectra,  ONLY:  xiabs,xanes_qua,xang_mom,xniter,xnitermax,xkvec,xepsilon,&
                        save_file_kind, calculated, time_limit
  USE atom,       ONLY : rgrid, msh
  !  use atom,        ONLY : &
  !       mesh,     &!mesh(ntypx) number of mesh points
  !       msh ,     &!msh(ntypx)the point at rcut=end of radial integration
  !       r   
  USE radin_mod
  USE uspp,   ONLY : vkb, nkb, okvan !CG
  USE ldaU,   ONLY : lda_plus_u
  USE basis,                ONLY : natomwfc
  !<CG>
  USE xspectra_paw_variables, ONLY : xspectra_paw_nhm
  !</CG>

  IMPLICIT NONE
  REAL(dp) core_wfn(ndmx)
  REAL(dp) a(xnitermax,1,nks),b(xnitermax,1,nks)     
  REAL (dp)  xnorm(1,nks)
  INTEGER is,ik,iabso,nr,ip,jp,l,j, ip_l,nrc,nt,na,ipx_0
  INTEGER ipx,ipy,ipz,nline,nrest,npw_partial,icrd
  INTEGER ncalcv(1,nks)
  INTEGER paw_iltonhb(0:paw_lmaxkb,xspectra_paw_nhm, ntyp) !CG

  REAL (dp) pref,prefb,v_of_0,xnorm_partial
  REAL (dp) norm,prefm2,prefm1,prefm0
  REAL (dp), ALLOCATABLE :: aux(:)
  COMPLEX(KIND=dp), ALLOCATABLE :: psi(:)
  COMPLEX(KIND=DP), EXTERNAL :: zdotc
  COMPLEX(dp), ALLOCATABLE :: paw_vkb_cplx(:,:)
  LOGICAL terminator
  REAL(dp) :: normps

  EXTERNAL zdscal

  COMPLEX(dp), ALLOCATABLE :: psiwfc(:), spsiwfc(:)
  LOGICAL :: recalc
  CHARACTER(LEN=4) :: verbosity
  REAL(DP), EXTERNAL ::  get_clock
  REAL(DP) :: timenow=0
  INTEGER :: nunfinished=0


  ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  !  Constant Definitions
  ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$  


  prefm2=SQRT(3.0/40.0)/3.0
  prefm1=prefm2
  prefm0=2.0*SQRT(1.0/40.0)/3.0

  pref=SQRT(2.d0)
  prefb=1.0/pref


  ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  !  Variable allocation and initialization
  ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$  

  ALLOCATE(aux(rgrid(xiabs)%mesh))!overdimensionated, necessary only up to msh
  ALLOCATE(psi(npwx))
  ALLOCATE (paw_vkb_cplx( npwx,  paw_nkb))
  ALLOCATE(xanes_qua(paw_recon(xiabs)%paw_nl(xang_mom)))
  ALLOCATE(psiwfc(npwx))

  xanes_qua(:)=0.d0
  psi(:)=(0.d0)



  ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  ! Radial Quadrupole Matrix element
  ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$  

  !
  !  First of all I build the quadrupole matrix element (radial part).
  !  Namely I compute the radial part, 
  !          <\phi_n|r|\psi_i>=int r^3 \phi_n \psi_i dr.
  !  Here only ps and ae wavefunctions from XX.recon
  !

  WRITE( stdout,*) 'Calculation Quadrupole matrix element'
  WRITE( stdout,*) 'There are ',paw_recon(xiabs)%paw_nl(xang_mom),' projectors/channel'
  WRITE( stdout,*) 'xang_mom=',xang_mom,' xiabs=',xiabs

  ! I check that the fondamental orthogonality condition of paw is satisfied:

  nr=msh(xiabs)  ! extended up to all the NON ZERO points in the mesh.


  ! I check that the core wf is correctly normalized

  IF(TRIM(verbosity).EQ.'high') THEN
     aux(1:nr)=core_wfn(1:nr)*core_wfn(1:nr)
     WRITE (stdout,'("norm of core wfc =",1f14.8)') SQRT(para_radin(aux(1:nr),rgrid(xiabs)%r(1:nr),nr))
  ENDIF

  ! and I calculate the radial integral

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


  DO ip=1,paw_recon(xiabs)%paw_nl(xang_mom) !
     WRITE( stdout,*)  'XANES quadrupole matrix element proj',ip,': ',xanes_qua(ip)
  ENDDO

  DEALLOCATE(aux)
  DEALLOCATE(psi)


  !
  !  We count the projectors for all the atoms
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


  ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  ! Beginning the loop over the k-points
  ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$  


  !*apsi  call set_vrs(vrs,vltot,vr,nrxx,nspin,doublegrid)
  !  CALL set_vrs(vrs,vltot,v%of_r,nrxx,nspin,doublegrid)
  !  CALL newd

  !<CG>
  ! set_vrs (vrs, vltot, vr, kedtau, kedtaur,nrxx, nspin, doublegrid)
  CALL set_vrs(vrs,vltot,v%of_r,kedtau, v%kin_r,dfftp%nnr,nspin,doublegrid)
  !</CG>




  DO ik=1,nks
     IF (calculated(1,ik).EQ.1) CYCLE

     timenow=get_clock( 'xanes' )
     CALL mp_bcast(timenow,root_pool, intra_pool_comm)

     IF( timenow.GT.time_limit) THEN
       nunfinished=1
       EXIT
     ENDIF

     WRITE(stdout,*)
     WRITE(stdout,*) 'Starting k-point : ', ik
     WRITE( stdout,'(" total cpu time spent up to now is ",F9.2," secs")') timenow

     current_k=ik
     IF(lsda) current_spin=isk(ik)

     !gk_sort  sort k-points and exit kinetic energies 
     CALL gk_sort(xk (1,ik),ngm,g,ecutwfc/tpiba2,npw,igk,g2kin)  !CHECK
     g2kin=g2kin*tpiba2                                          !CHECK

     npw_partial = npw
     CALL mp_sum( npw_partial, intra_pool_comm )
     IF(xniter.ge.npw_partial) THEN
        xniter = npw_partial
        WRITE(stdout,*) 'Hilbert space is saturated'
        WRITE(stdout,*) 'xniter is set equal to ',npw_partial
        WRITE(stdout,*) 'Hint: Increase Kinetic Energy cutoff in your SCF simulation'
        !        CALL stop_pp
     ENDIF

     !<CG>
     CALL init_gipaw_2(npw,igk,xk(1,ik),paw_vkb)
     !</CG>
     if(.not.lda_plus_u) CALL init_us_2(npw,igk,xk(1,ik),vkb)
     IF (lda_plus_u) CALL orthoUwfc_k(ik)


     !
     ! Here I define human projectors
     !
     ! paw_vkb(1:npw,ipx) are real spherical harmonics. The real spherical harmonics are
     ! defined as
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
     !



     DO ip=1,paw_recon(xiabs)%paw_nl(xang_mom)   

        ipx=ipx_0+paw_iltonhb(xang_mom,ip,xiabs)
        paw_vkb_cplx(1:npw,ipx)= paw_vkb(1:npw,ipx)   !m=0
        paw_vkb_cplx(1:npw,ipx+1)=       &
             (paw_vkb(1:npw,ipx+1)+(0.d0,1.d0)*paw_vkb(1:npw,ipx+2))/SQRT(2.0)   !m=+1
        paw_vkb_cplx(1:npw,ipx+2)=       &
             -(paw_vkb(1:npw,ipx+1)-(0.d0,1.d0)*paw_vkb(1:npw,ipx+2))/SQRT(2.0)    !m=-1
        paw_vkb_cplx(1:npw,ipx+3)=       &
             (paw_vkb(1:npw,ipx+3)+(0.d0,1.d0)*paw_vkb(1:npw,ipx+4))/SQRT(2.0)   !m=+2
        paw_vkb_cplx(1:npw,ipx+4)=       &
             (paw_vkb(1:npw,ipx+3)-(0.d0,1.d0)*paw_vkb(1:npw,ipx+4))/SQRT(2.0)    !m=-2


     ENDDO

     !      WARNING, storage of spherical harmonics is the following:
     !         given Y_{l,m}=P_{lm}exp(im\phi) one has
     !                                                         counter
     !   l, m=0  -------> Y_{l,0}                              1     z
     !   l, m=+1 -------> Y_{l,1}+Y_{l,-1} or cos(m\phi)       2     
     !   l, m=-1 -------> Y_{l,1}-Y_{l,-1} or sin(m\phi)       3     
     !  .....
     !   l, m=+l -------> Y_{l,1}+Y_{l,-1} or cos(m\phi)      2*l
     !   l, m=-l -------> Y_{l,1}-Y_{l,-1} or sin(m\phi)      2*l+1


     ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
     ! Angular Matrix element
     ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$  

     psiwfc(:)=(0.d0,0.d0)

     DO ip=1,paw_recon(xiabs)%paw_nl(xang_mom)   
        ipx=ipx_0+paw_iltonhb(xang_mom,ip,xiabs)


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


     ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
     ! Starting Lanczos
     ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$  

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

     !   
     !   I normalize the wavefunction pwswfc(1:npw)
     !

     CALL mp_sum( xnorm_partial, intra_pool_comm )

     xnorm(1,ik)=SQRT(xnorm_partial)
     WRITE( stdout,*) 'norm initial vector=',xnorm(1,ik)
     norm=1.d0/xnorm(1,ik)

     CALL zdscal(npw,norm,psiwfc,1)
     !
     !      Then I call the lanczos routine
     !
     IF (okvan) THEN
        CALL lanczos_uspp(a(:,1,ik),b(:,1,ik),psiwfc,ncalcv(1,ik), terminator)
     ELSE
        CALL lanczos(a(:,1,ik),b(:,1,ik),psiwfc,ncalcv(1,ik), terminator)
     ENDIF
     !
     !      Then I write small report of the lanczos results
     !

     IF(TRIM(verbosity).EQ.'high') THEN
        WRITE( stdout,*) '-----------------------------------------'
        WRITE( stdout,*) 'k-point number =',ik
        WRITE( stdout,*) 'Norm of the initial vector =',xnorm(1,ik)
        WRITE( stdout,*) 'Number of iterations =',ncalcv(1,ik)
        WRITE( stdout,*) 'k-point coordinate, isk'
        WRITE( stdout,'(3f12.6,1i2)') xk(1,ik),xk(2,ik),xk(3,ik),isk(ik)
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
        WRITE( stdout,*) '-----------------------------------------'
     ENDIF

     CALL deallocate_bec_type (becp) ! CG

     calculated(1,ik)=1

  ENDDO   !en loop over k points

  CALL mp_barrier(world_comm)
  CALL mp_sum(nunfinished, world_comm)
  IF (nunfinished >= 1) THEN
    save_file_kind='unfinished'
    write(stdout,*) 'calculation not finished'
  ENDIF


  ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  ! Array deallocation
  ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$  

  DEALLOCATE(psiwfc)
  DEALLOCATE (paw_vkb_cplx)
  DEALLOCATE(xanes_qua)



END SUBROUTINE xanes_quadrupole

!<CG>
SUBROUTINE lanczos (a,b,psi,ncalcv, terminator)
  USE kinds, ONLY : DP
  USE constants, ONLY : rytoev
  USE wvfct,  ONLY:  npwx,nbndx, nbnd,npw,igk,g2kin
  !USE wavefunctions_module, ONLY : psic
  USE becmod, ONLY:becp
  USE uspp,   ONLY : vkb, nkb
  USE cell_base, ONLY:omega
  USE xspectra, ONLY : xniter, xnepoint, xcheck_conv,xnitermax,xemin,xemax,xgamma,xerror
  USE mp_global,  ONLY : intra_pool_comm  
  USE mp,         ONLY : mp_sum
  USE io_global,       ONLY : stdout

  IMPLICIT NONE
  ! INPUT
  INTEGER ncalcv
  REAL(dp) :: a(xnitermax),b(xnitermax)
  COMPLEX(KIND=DP) :: psi(npwx)
  LOGICAL terminator

  !INTERNAL
  LOGICAL converge
  LOGICAL iconv
  INTEGER ibnd,j,i,m
  REAL(KIND=dp)  norm,error,xemin_ryd,xemax_ryd,xgamma_ryd
  COMPLEX(KIND=dp):: ac,bc
  REAL(KIND=dp), ALLOCATABLE :: comp(:)
  COMPLEX(KIND=dp), ALLOCATABLE :: hpsi(:),u(:)
  REAL (KIND=dp) :: ddot
  COMPLEX(KIND=DP) :: zdotc
  EXTERNAL zdotc,ddot
  External h_psi

  ALLOCATE(hpsi(npwx))
  ALLOCATE(u(npwx))
  ALLOCATE(comp(xnepoint))

  hpsi(:)=(0.d0,0.d0)
  u(:)=(0.d0,0.d0)
  a(:)=0.d0
  b(:)=0.d0

  xemax_ryd=xemax/rytoev
  xemin_ryd=xemin/rytoev
  xgamma_ryd=xgamma/rytoev

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

  b(1) = zdotc(npw,hpsi,1,hpsi,1)
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
        IF(converge(a,b,i,comp,error,xemin_ryd,xemax_ryd,xgamma_ryd,xnepoint,xerror,terminator)) THEN
           WRITE( stdout,*) 'CONVERGED at iter ',i,' with error=',error
           ncalcv=i
           iconv=.true.
           EXIT
        ELSE
           WRITE( stdout,*) 'Estimated error at iter ',i,' is ', error
        ENDIF
     ENDIF


  ENDDO

  IF(.NOT.iconv) THEN
     WRITE( stdout,*) 'XANES not converged after', i-1, ' iterations'
     WRITE( stdout,*) 'Estimated final error after ',i-1,'iterations is ', &
          converge(a,b,i-1,comp,error,xemin_ryd,xemax_ryd,xgamma_ryd,xnepoint,xerror,terminator)
     ncalcv=i-1
  ENDIF



  DEALLOCATE(hpsi)
  DEALLOCATE(u)
  DEALLOCATE(comp)


END SUBROUTINE lanczos


!<CG>
SUBROUTINE lanczos_uspp (a,b,psi,ncalcv, terminator)
  USE kinds, ONLY : DP
  USE wvfct,  ONLY:  npwx,nbndx, nbnd,npw,igk,g2kin
  USE becmod, ONLY: becp, calbec
  USE constants, ONLY: rytoev
  USE uspp,   ONLY : vkb, nkb
  USE cell_base, ONLY:omega
  USE xspectra, ONLY : xniter, xnepoint,xcheck_conv,xnitermax,xemin,xemax,xgamma,xerror
  USE mp_global,  ONLY : intra_pool_comm
  USE mp,         ONLY : mp_sum
  USE io_global,       ONLY : stdout

  IMPLICIT NONE
  ! INPUT
  INTEGER ncalcv
  REAL(dp) :: a(xnitermax),b(xnitermax)
  COMPLEX(KIND=DP) :: psi(npwx)
  LOGICAL terminator

  !INTERNAL
  LOGICAL converge
  LOGICAL iconv
  INTEGER ibnd,j,i,m
  REAL(KIND=dp)  norm,error,xemin_ryd,xemax_ryd,xgamma_ryd
  COMPLEX(KIND=dp):: ac,bc
  REAL(KIND=dp), ALLOCATABLE :: comp(:)
  COMPLEX(KIND=dp), ALLOCATABLE :: u(:), v1(:), v2(:), v3(:)
  REAL (KIND=dp) :: ddot
  COMPLEX(KIND=DP) :: zdotc
  EXTERNAL zdotc,ddot
  External h_psi
  LOGICAL recalc

  COMPLEX(KIND=DP) :: vecteuraux1(npwx,1), vecteuraux2(npwx,1)

  ALLOCATE(v1(npwx))
  ALLOCATE(u(npwx))
  ALLOCATE(comp(xnepoint))
  ALLOCATE(v2(npwx))
  ALLOCATE(v3(npwx))

  v1(:)=(0.d0,0.d0)
  v2(:)=(0.d0,0.d0)
  u(:)=(0.d0,0.d0)
  a(:)=0.d0
  b(:)=0.d0

  xemax_ryd=xemax/rytoev
  xemin_ryd=xemin/rytoev
  xgamma_ryd=xgamma/rytoev

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
        IF(converge(a,b,i,comp,error,xemin_ryd,xemax_ryd,xgamma_ryd,xnepoint,xerror,terminator)) THEN
           WRITE( stdout,*) 'CONVERGED at iter ',i,' with error=',error
           ncalcv=i
           iconv=.true.
           EXIT
        ELSE
           WRITE( stdout,*) 'Estimated error at iter ',i,' is ', error
        ENDIF
     ENDIF

  ENDDO

  IF(.NOT.iconv) THEN
     WRITE( stdout,*) 'XANES not converged after', i-1, ' iterations'
     WRITE( stdout,*) 'Estimated final error after ',i-1,'iterations is ', &
          converge(a,b,i-1,comp,error,xemin_ryd,xemax_ryd,xgamma_ryd,xnepoint,xerror,terminator)
     ncalcv=i-1
  ENDIF

  DEALLOCATE(v1,v2,u)
  DEALLOCATE(comp)

END SUBROUTINE lanczos_uspp

!</CG>



LOGICAL FUNCTION converge(a,b,m,comp,estimated_error,xemin,xemax,xgamma,xnepoint,xerror,use_term)
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
!
! ------------------------------------------------------------------------
!
FUNCTION continued_fraction(a,b,e,gamma,m, term)
  USE kinds, ONLY : dp
  USE xspectra, ONLY : xnitermax, xcheck_conv
  USE io_global,       ONLY : stdout
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
!

SUBROUTINE read_core_abs(filename,core_wfn)
  USE kinds, ONLY : DP
  USE atom,        ONLY : rgrid
  !USE atom,        ONLY : mesh     !mesh(ntypx) number of mesh points
  USE xspectra,       ONLY : xiabs
  IMPLICIT NONE
  INTEGER :: i
  CHARACTER (LEN=80) :: filename
  REAL(KIND=dp):: x
  REAL(KIND=dp):: core_wfn(*)

  WRITE(6,*) 'xmesh=',rgrid(xiabs)%mesh
  open(unit=33,file=filename,form='formatted',status='old')
  rewind(33)
  READ(33,*)
  DO i=1, rgrid(xiabs)%mesh
     !   read(33,'(2f16.8)') x,core_wfn(i)
     READ(33 ,*) x,core_wfn(i)
  ENDDO
  close(33)
END SUBROUTINE read_core_abs


! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!  Continuum fraction and plotting dipolar part
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$  


SUBROUTINE plot_xanes_dipole(a,b,xnorm,ncalcv, terminator, core_energy)
  USE kinds, ONLY : DP
  USE constants, ONLY : rytoev, pi
  USE xspectra, ONLY : xang_mom,xemax,xemin,xiabs,xnepoint,xgamma,xonly_plot
  !*apsi  USE uspp_param, ONLY : psd  !psd(ntypx) label for the atoms 
  USE klist, ONLY : nelec, &
       nkstot,            & ! total number of k-points
       nks,               & ! number of k-points per pool
       xk,                & ! k-points coordinates
       wk                   ! k-points weight
  USE ener, ONLY : ef
  USE io_global,       ONLY : stdout,ionode  
  USE mp_global,  ONLY : intra_pool_comm, inter_pool_comm !CG
  USE xspectra, ONLY : xnitermax,xkvec,xepsilon
  USE cut_valence_green, ONLY : cut_occ_states, cut_desmooth, memu, meml, cut_nmemu, cut_nmeml
  USE lsda_mod,    ONLY : nspin,isk
  !<CG>
  USE mp,                   ONLY : mp_sum
  USE uspp_param,           ONLY : upf
  USE gamma_variable_mod  , ONLY : gamma_tab, gamma_mode
  !</CG>
  IMPLICIT NONE
  LOGICAL terminator
  INTEGER iestart
  INTEGER ncalcv(1,nks)
  REAL(dp) a(xnitermax,1,nks),b(xnitermax,1,nks)     
  REAL (dp)  xnorm(1,nks)
  REAL(KIND=dp) alpha2,core_energy
  INTEGER :: i,ik,n,icoord           !loops
  INTEGER :: lmax
  REAL (dp) :: mygetK,e_1s,energy,de,mod_xgamma,xemax_ryd,xemin_ryd,xgamma_ryd
  REAL (dp) :: tmp_var
  REAL (dp) :: Intensity_coord(1,xnepoint,nspin)
  !   , Intensity_tot(xnepoint)
  REAL(dp)  :: continued_fraction, paste_fermi, desmooth,t1,t2,f1,f2,df1,df2,poly(4) !CG 
  LOGICAL :: first

  ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  !  constant and initialization
  ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$  

  IF(xonly_plot) THEN
     e_1s=core_energy           !the K-edge in eV.
  ELSE
     !<CG>
     e_1s=mygetK(upf(xiabs)%psd) !mygetK gets the K-edge in eV.
     !</CG>
  ENDIF


  e_1s=e_1s/rytoeV            !This is in Rydberg

  alpha2=4.d0*pi/137.04

  desmooth=cut_desmooth/rytoev  ! This is in rydberg

  !
  ! small output
  !

  WRITE( stdout,"(/'CALCULATION XANES SPECTRA')")
  WRITE( stdout,"('---------------------')")
  WRITE( stdout,"(/' Final state angular momentum:',1x,i3)") xang_mom
  IF (TRIM(gamma_mode).EQ.'constant')  WRITE( stdout,"(' Broadening parameter (in eV):',1x,g20.12)") xgamma     !check
  WRITE( stdout,"(' Nb of energy points:',1x,i5)") xnepoint
  WRITE( stdout,"(' Maximum energy (in eV):',1x,g20.12)") xemax            !check
  WRITE( stdout,"(' Minimum energy (in eV):',1x,g20.12)") xemin            !check
  WRITE( stdout,"(' Binding energy of the 1s level (in eV):',1x,g20.12)") -e_1s*rytoeV

  IF( ionode ) THEN
     open (unit=277,file='xanes.dat',form='formatted',status='unknown')
     rewind(277)
     WRITE(277,"('# final state angular momentum:',1x,i3)") xang_mom
     IF (TRIM(gamma_mode).EQ.'constant')     WRITE(277,"('# broadening parameter (in eV):',1x,g20.12)") xgamma      !

     WRITE(277,"('# absorbing atom type:',i4)") xiabs

     IF(nspin.GT.1) THEN
        WRITE(277,"('# Energy (eV)   sigmatot   sigmaup    sigmadown ')")
     ELSE
        WRITE(277,"('# Energy (eV)   sigma')")
     ENDIF

  ENDIF

  !
  ! I convert in Rydberg most of the relevant quantities
  !

  xemax_ryd=xemax/rytoeV+ef
  xemin_ryd=xemin/rytoeV+ef
  xgamma_ryd=xgamma/rytoeV
  WRITE(stdout,*) 'xemin(ryd)=',xemin_ryd
  WRITE(stdout,*) 'xemax(ryd)=',xemax_ryd
  IF (TRIM(gamma_mode).EQ.'constant') THEN
     WRITE(stdout,*) 'gamma in rydberg=',xgamma_ryd
  ELSE
     WRITE(stdout,*) 'nonconstant gamma'
  ENDIF


  ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  !  Continuum fraction
  ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$  

  Intensity_coord(:,:,:)=0.d0


  de = (xemax_ryd-xemin_ryd) / REAL(xnepoint-1)

  IF (TRIM(gamma_mode).EQ.'constant') THEN

     IF(cut_occ_states) THEN
        ALLOCATE(memu(cut_nmemu,2))
        ALLOCATE(meml(cut_nmeml,2))
        iestart=(ef-xemin_ryd)/de
        DO ik=1,nks
           first=.true.  ! to erase the memory of paste_fermi
           !<CG>
           t1=ef-desmooth
           f1=paste_fermi(t1,ef,a(1,1,ik),b(1,1,ik),xgamma_ryd,ncalcv(1,ik)-1,terminator, first)
           df1=paste_fermi(t1-de,ef,a(1,1,ik),b(1,1,ik),xgamma_ryd,ncalcv(1,ik)-1,terminator, first)
           df1=(f1-df1)/de
           t2=ef+desmooth
           f2=continued_fraction(a(1,1,ik),b(1,1,ik),t2,xgamma_ryd,ncalcv(1,ik)-1,terminator)&
                +paste_fermi(t2,ef,a(1,1,ik),b(1,1,ik),xgamma_ryd,ncalcv(1,ik)-1,terminator, first)
           df2=continued_fraction(a(1,1,ik),b(1,1,ik),t2+de,xgamma_ryd,ncalcv(1,ik)-1,terminator)&
                +paste_fermi(t2+de,ef,a(1,1,ik),b(1,1,ik),xgamma_ryd,ncalcv(1,ik)-1,terminator, first)
           df2=(df2-f2)/de
           CALL determine_polycut(t1,t2,f1,f2,df1,df2,poly) ! calculates interpolation polynome

           DO n=1,xnepoint
              energy=xemin_ryd+de*(n-1)
              IF ((energy-ef<desmooth).AND.(energy-ef>-desmooth)) THEN  ! interpolation 
                 tmp_var=poly(1)+poly(2)*energy+poly(3)*energy**2+poly(4)*energy**3
                 tmp_var=tmp_var*xnorm(1,ik)*xnorm(1,ik)
                 !</CG>
              ELSE
                 tmp_var=0.d0
                 IF (n>iestart) THEN
                    tmp_var=  &
                         continued_fraction(a(1,1,ik),b(1,1,ik),energy,xgamma_ryd,ncalcv(1,ik)-1,terminator)*  &
                         xnorm(1,ik)*xnorm(1,ik)
                 ENDIF
                 tmp_var = tmp_var + paste_fermi(energy,ef,a(1,1,ik),b(1,1,ik),xgamma_ryd,ncalcv(1,ik)-1,terminator, first) &
                      *xnorm(1,ik)*xnorm(1,ik)
              ENDIF
              !            Intensity_tot(n)=Intensity_tot(n)+tmp_var*wk(ik)
              Intensity_coord(1,n,isk(ik)) = Intensity_coord(1,n,isk(ik))+tmp_var*wk(ik)
           ENDDO
        ENDDO
        DEALLOCATE(memu)
        DEALLOCATE(meml)


     ELSE
        DO ik=1,nks
           DO n=1,xnepoint
              energy=xemin_ryd+de*(n-1)
              tmp_var=  &
                   continued_fraction(a(1,1,ik),b(1,1,ik),energy,xgamma_ryd,ncalcv(1,ik)-1, terminator)*  &
                   xnorm(1,ik)*xnorm(1,ik)
              Intensity_coord(1,n,isk(ik)) = Intensity_coord(1,n,isk(ik))+tmp_var*wk(ik)
           ENDDO
        ENDDO
     ENDIF

  ELSE ! nonconstant gamma

     IF(cut_occ_states) THEN
        ALLOCATE(memu(cut_nmemu,2))
        ALLOCATE(meml(cut_nmeml,2))
        iestart=(ef-xemin_ryd)/de
        DO ik=1,nks
           first=.true.  ! to erase the memory of paste_fermi


           xgamma_ryd=gamma_tab(iestart)
           t1=ef-desmooth
           f1=paste_fermi(t1,ef,a(1,1,ik),b(1,1,ik),xgamma_ryd,ncalcv(1,ik)-1,terminator, first)
           df1=paste_fermi(t1-de,ef,a(1,1,ik),b(1,1,ik),xgamma_ryd,ncalcv(1,ik)-1,terminator, first)
           df1=(f1-df1)/de
           t2=ef+desmooth
           f2=continued_fraction(a(1,1,ik),b(1,1,ik),t2,xgamma_ryd,ncalcv(1,ik)-1,terminator)&
                +paste_fermi(t2,ef,a(1,1,ik),b(1,1,ik),xgamma_ryd,ncalcv(1,ik)-1,terminator, first)
           df2=continued_fraction(a(1,1,ik),b(1,1,ik),t2+de,xgamma_ryd,ncalcv(1,ik)-1,terminator)&
                +paste_fermi(t2+de,ef,a(1,1,ik),b(1,1,ik),xgamma_ryd,ncalcv(1,ik)-1,terminator, first)
           df2=(df2-f2)/de
           CALL determine_polycut(t1,t2,f1,f2,df1,df2,poly)

           DO n=1,xnepoint
              energy=xemin_ryd+de*(n-1)
              xgamma_ryd=gamma_tab(n)
              IF ((energy-ef<desmooth).AND.(energy-ef>-desmooth)) THEN  ! interpolation
                 tmp_var=poly(1)+poly(2)*energy+poly(3)*energy**2+poly(4)*energy**3
                 tmp_var=tmp_var*xnorm(1,ik)*xnorm(1,ik)
              ELSE
                 tmp_var=0.d0
                 IF (n>iestart) THEN
                    tmp_var=  &
                         continued_fraction(a(1,1,ik),b(1,1,ik),energy,xgamma_ryd,ncalcv(1,ik)-1,terminator)*  &
                         xnorm(1,ik)*xnorm(1,ik)
                 ENDIF
                 tmp_var = tmp_var + paste_fermi(energy,ef,a(1,1,ik),b(1,1,ik),xgamma_ryd,ncalcv(1,ik)-1,terminator, first) &
                      *xnorm(1,ik)*xnorm(1,ik)
              ENDIF
              !            Intensity_tot(n)=Intensity_tot(n)+tmp_var*wk(ik)
              Intensity_coord(1,n,isk(ik)) = Intensity_coord(1,n,isk(ik))+tmp_var*wk(ik)
           ENDDO
        ENDDO
        DEALLOCATE(memu)
        DEALLOCATE(meml)


     ELSE
        DO ik=1,nks
           DO n=1,xnepoint
              energy=xemin_ryd+de*(n-1)
              xgamma_ryd=gamma_tab(n)
              tmp_var=  &
                   continued_fraction(a(1,1,ik),b(1,1,ik),energy,xgamma_ryd,ncalcv(1,ik)-1, terminator)*  &
                   xnorm(1,ik)*xnorm(1,ik)
              Intensity_coord(1,n,isk(ik)) = Intensity_coord(1,n,isk(ik))+tmp_var*wk(ik)
           ENDDO
        ENDDO
     ENDIF

  ENDIF ! gamma_mode

  !  CALL poolreduce( nspin*xnepoint, Intensity_coord )

  !<CG>  replaces poolreduce
#ifdef __MPI
  CALL mp_sum ( Intensity_coord, inter_pool_comm )
#endif
  !</CG>

  ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  !  Plotting xanes spectrum
  ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$  

  IF(ionode) THEN
     IF(nspin.EQ.1) THEN
        DO n=1,xnepoint
           energy=xemin_ryd+de*(n-1)
           Intensity_coord(:,n,:)=Intensity_coord(:,n,:)*(energy+e_1s)*alpha2 !
           WRITE(277,'(2f14.8)') (energy-ef)*rytoeV, Intensity_coord(1,n,1)
        ENDDO
     ELSEIF(nspin.EQ.2) THEN
        DO n=1,xnepoint
           energy=xemin_ryd+de*(n-1)
           Intensity_coord(:,n,:)=Intensity_coord(:,n,:)*(energy+e_1s)*alpha2 !
           WRITE(277,'(4f14.8)') (energy-ef)*rytoev, &
                Intensity_coord(1,n,1)+Intensity_coord(1,n,2),&
                Intensity_coord(1,n,1),Intensity_coord(1,n,2)
        ENDDO
     ENDIF

     close(277)
  ENDIF

END SUBROUTINE plot_xanes_dipole



! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!  Continuum fraction and plotting quadrupolar part
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$  



SUBROUTINE plot_xanes_quadrupole(a,b,xnorm,ncalcv, terminator,core_energy)
  USE kinds, ONLY : DP
  USE constants, ONLY: pi, rytoev
  USE xspectra, ONLY : xang_mom,xemax,xemin,xiabs,xnepoint,xgamma,xonly_plot
  !*apsi  USE uspp_param, ONLY : psd  !psd(ntypx) label for the atoms 
  USE klist, ONLY : nelec, &
       nkstot,            & ! total number of k-points
       nks,               & ! number of k-points per pool
       xk,                & ! k-points coordinates
       wk                   ! k-points weight
  USE ener, ONLY : ef
  USE io_global,       ONLY : stdout,ionode  
  USE mp_global,  ONLY : intra_pool_comm, inter_pool_comm
  USE xspectra, ONLY : xnitermax,xepsilon
  USE cut_valence_green , ONLY : cut_occ_states, cut_desmooth, cut_nmemu, cut_nmeml, memu, meml
  USE lsda_mod,    ONLY : nspin,isk
  !<CG>
  USE mp,                   ONLY : mp_sum
  USE uspp_param,           ONLY : upf
  USE gamma_variable_mod, ONLY : gamma_tab, gamma_mode
  !</CG>
  IMPLICIT NONE
  LOGICAL terminator
  INTEGER iestart
  INTEGER ncalcv(1,nks)
  REAL(dp) a(xnitermax,1,nks),b(xnitermax,1,nks)     
  REAL (dp)  xnorm(1,nks)
  REAL(KIND=dp) alpha2,constantqua
  INTEGER :: i,ik,n,icoord           !loops
  INTEGER :: lmax
  REAL (dp) :: mygetK,e_1s,energy,de,mod_xgamma,xemax_ryd,xemin_ryd,xgamma_ryd
  REAL (dp) :: tmp_var,core_energy
  REAL (dp) :: Intensity_tot(xnepoint,nspin)
  REAL(dp)  :: continued_fraction, paste_fermi, desmooth,t1,t2,f1,f2,df1,df2,poly(4)
  LOGICAL :: first

  ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  !  constant and initialization
  ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$  

  constantqua=pi/(137.04*137.04*137.04)
  IF(xonly_plot) THEN
     e_1s=core_energy           !the K-edge in eV.
  ELSE
     !<CG>
     e_1s=mygetK(upf(xiabs)%psd) !mygetK gets the K-edge in eV.
     !</CG>
  ENDIF
  e_1s=e_1s/rytoeV            !This is in Rydberg

  desmooth=cut_desmooth/rytoev  ! This is in rydberg

  alpha2=4.d0*pi/137.04

  !
  ! small output
  !

  WRITE( stdout,"(/'CALCULATION XANES SPECTRA')")
  WRITE( stdout,"('---------------------')")
  WRITE( stdout,"(/' Final state angular momentum:',1x,i3)") xang_mom
  IF (TRIM(gamma_mode).EQ.'constant')   WRITE( stdout,"(' Broadening parameter (in eV):',1x,g20.12)") xgamma     !check
  WRITE( stdout,"(' Nb of energy points:',1x,i5)") xnepoint
  WRITE( stdout,"(' Maximum energy (in eV):',1x,g20.12)") xemax            !check
  WRITE( stdout,"(' Minimum energy (in eV):',1x,g20.12)") xemin            !check
  WRITE( stdout,"(' Binding energy of the 1s level (in eV):',1x,g20.12)") -e_1s*rytoev


  IF( ionode ) THEN
     open (unit=277,file='xanes.dat',form='formatted',status='unknown')
     rewind(277)
     WRITE(277,"('# final state angular momentum:',1x,i3)") xang_mom
     IF (TRIM(gamma_mode).EQ.'constant')  WRITE(277,"('# broadening parameter (in eV):',1x,g20.12)") xgamma      

     WRITE(277,"('# absorbing atom type:',i4)") xiabs
     IF(nspin.EQ.1) THEN
        WRITE(277,"('# Energy (eV)   sigmatot')")
     ELSEIF(nspin.EQ.2) THEN
        WRITE(277,"('# Energy (eV)   sigmatot, sigmaup, sigmadown')")
     ENDIF
  ENDIF

  !
  ! I convert in Rydberg most of the relevant quantities
  !

  xemax_ryd=xemax/rytoev+ef
  xemin_ryd=xemin/rytoev+ef
  xgamma_ryd=xgamma/rytoev

  WRITE(stdout,*) 'xemin(ryd)=',xemin_ryd
  WRITE(stdout,*) 'xemax(ryd)=',xemax_ryd
  IF (TRIM(gamma_mode).EQ.'constant') THEN
     WRITE(stdout,*) 'gamma in rydberg=',xgamma_ryd
  ELSE
     WRITE(stdout,*) 'gamma from file'
  ENDIF


  !  WRITE(stdout,*) 'terminator=', terminator


  ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  !  Continuum fraction
  ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$  


  Intensity_tot(:,:)=0.d0

  de = (xemax_ryd-xemin_ryd) / REAL(xnepoint-1)

  IF (TRIM(gamma_mode).EQ.'constant') THEN
     IF(cut_occ_states) THEN
        ALLOCATE(memu(cut_nmemu,2))
        ALLOCATE(meml(cut_nmeml,2))
        DO ik=1,nks
           iestart=(ef-xemin_ryd)/de
           first=.true.
           !<CG>
           t1=ef-desmooth
           f1=paste_fermi(t1,ef,a(1,1,ik),b(1,1,ik),xgamma_ryd,ncalcv(1,ik)-1,terminator, first)
           df1=paste_fermi(t1-de,ef,a(1,1,ik),b(1,1,ik),xgamma_ryd,ncalcv(1,ik)-1,terminator, first)
           df1=(f1-df1)/de
           t2=ef+desmooth
           f2=continued_fraction(a(1,1,ik),b(1,1,ik),t2,xgamma_ryd,ncalcv(1,ik)-1,terminator)&
                +paste_fermi(t2,ef,a(1,1,ik),b(1,1,ik),xgamma_ryd,ncalcv(1,ik)-1,terminator, first)
           df2=continued_fraction(a(1,1,ik),b(1,1,ik),t2+de,xgamma_ryd,ncalcv(1,ik)-1,terminator)&
                +paste_fermi(t2+de,ef,a(1,1,ik),b(1,1,ik),xgamma_ryd,ncalcv(1,ik)-1,terminator, first)
           df2=(df2-f2)/de
           CALL determine_polycut(t1,t2,f1,f2,df1,df2,poly)

           DO n=1,xnepoint
              energy=xemin_ryd+de*(n-1)
              IF ((energy-ef<desmooth).AND.(energy-ef>-desmooth)) THEN  !interpolation
                 tmp_var=poly(1)+poly(2)*energy+poly(3)*energy**2+poly(4)*energy**3
                 tmp_var=tmp_var*xnorm(1,ik)*xnorm(1,ik)
                 !</CG>
              ELSE
                 tmp_var=0.d0
                 IF (n>iestart) THEN
                    tmp_var=  &
                         continued_fraction(a(1,1,ik),b(1,1,ik),energy,xgamma_ryd,ncalcv(1,ik)-1,terminator)*  &
                         xnorm(1,ik)*xnorm(1,ik)
                 ENDIF
                 tmp_var=tmp_var+paste_fermi(energy,ef,a(1,1,ik),b(1,1,ik),xgamma_ryd,ncalcv(1,ik)-1,terminator, first)&
                      *xnorm(1,ik)*xnorm(1,ik)
              ENDIF
              Intensity_tot(n,isk(ik))=Intensity_tot(n,isk(ik))+tmp_var*wk(ik)
           ENDDO
        ENDDO
        DEALLOCATE(memu)
        DEALLOCATE(meml)

     ELSE
        !     WRITE(stdout,*) 'in plor_
        DO ik=1,nks
           DO n=1,xnepoint
              energy=xemin_ryd+de*(n-1)
              tmp_var=  &
                   continued_fraction(a(1,1,ik),b(1,1,ik),energy,xgamma_ryd,ncalcv(1,ik)-1, terminator)*  &
                   xnorm(1,ik)*xnorm(1,ik)
              Intensity_tot(n,isk(ik))=Intensity_tot(n,isk(ik))+tmp_var*wk(ik)
           ENDDO
        ENDDO
     ENDIF

  ELSE ! nonconstant gamma

     IF(cut_occ_states) THEN
        ALLOCATE(memu(cut_nmemu,2))
        ALLOCATE(meml(cut_nmeml,2))
        DO ik=1,nks
           iestart=(ef-xemin_ryd)/de
           first=.true.
           !<CG>
           xgamma_ryd=gamma_tab(iestart)
           t1=ef-desmooth
           f1=paste_fermi(t1,ef,a(1,1,ik),b(1,1,ik),xgamma_ryd,ncalcv(1,ik)-1,terminator, first)
           df1=paste_fermi(t1-de,ef,a(1,1,ik),b(1,1,ik),xgamma_ryd,ncalcv(1,ik)-1,terminator, first)
           df1=(f1-df1)/de
           t2=ef+desmooth
           f2=continued_fraction(a(1,1,ik),b(1,1,ik),t2,xgamma_ryd,ncalcv(1,ik)-1,terminator)&
                +paste_fermi(t2,ef,a(1,1,ik),b(1,1,ik),xgamma_ryd,ncalcv(1,ik)-1,terminator, first)
           df2=continued_fraction(a(1,1,ik),b(1,1,ik),t2+de,xgamma_ryd,ncalcv(1,ik)-1,terminator)&
                +paste_fermi(t2+de,ef,a(1,1,ik),b(1,1,ik),xgamma_ryd,ncalcv(1,ik)-1,terminator, first)
           df2=(df2-f2)/de
           CALL determine_polycut(t1,t2,f1,f2,df1,df2,poly)

           DO n=1,xnepoint
              energy=xemin_ryd+de*(n-1)
              xgamma_ryd=gamma_tab(n)
              IF ((energy-ef<desmooth).AND.(energy-ef>-desmooth)) THEN  !interpolation
                 tmp_var=poly(1)+poly(2)*energy+poly(3)*energy**2+poly(4)*energy**3
                 tmp_var=tmp_var*xnorm(1,ik)*xnorm(1,ik)
                 !</CG>
              ELSE
                 tmp_var=0.d0
                 IF (n>iestart) THEN
                    tmp_var=  &
                         continued_fraction(a(1,1,ik),b(1,1,ik),energy,xgamma_ryd,ncalcv(1,ik)-1,terminator)*  &
                         xnorm(1,ik)*xnorm(1,ik)
                 ENDIF
                 tmp_var=tmp_var+paste_fermi(energy,ef,a(1,1,ik),b(1,1,ik),xgamma_ryd,ncalcv(1,ik)-1,terminator, first)&
                      *xnorm(1,ik)*xnorm(1,ik)
              ENDIF
              Intensity_tot(n,isk(ik))=Intensity_tot(n,isk(ik))+tmp_var*wk(ik)
           ENDDO
        ENDDO
        DEALLOCATE(memu)
        DEALLOCATE(meml)

     ELSE
        !     WRITE(stdout,*) 'in plor_
        DO ik=1,nks
           DO n=1,xnepoint
              energy=xemin_ryd+de*(n-1)
              xgamma_ryd=gamma_tab(n)
              tmp_var=  &
                   continued_fraction(a(1,1,ik),b(1,1,ik),energy,xgamma_ryd,ncalcv(1,ik)-1, terminator)*  &
                   xnorm(1,ik)*xnorm(1,ik)
              Intensity_tot(n,isk(ik))=Intensity_tot(n,isk(ik))+tmp_var*wk(ik)
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


  ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  !  Plotting xanes spectrum
  ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$  

  IF(ionode) THEN
     IF(nspin.EQ.1) THEN
        DO n=1,xnepoint
           energy=xemin_ryd+de*(n-1)
           Intensity_tot(n,:)=Intensity_tot(n,:)*(energy+e_1s)*(energy+e_1s)*(energy+e_1s)*constantqua     !normalized
           WRITE(277,'(2f14.8)') (energy-ef)*rytoev,Intensity_tot(n,:)
        ENDDO
     ELSEIF(nspin.EQ.2) THEN
        DO n=1,xnepoint
           energy=xemin_ryd+de*(n-1)
           Intensity_tot(n,:)=Intensity_tot(n,:)*(energy+e_1s)*(energy+e_1s)*(energy+e_1s)*constantqua     !normalized
           WRITE(277,'(4f14.8)') (energy-ef)*rytoev,Intensity_tot(n,1)+Intensity_tot(n,2),&
                Intensity_tot(n,1),Intensity_tot(n,2)
        ENDDO
     ENDIF
     close(277)
  ENDIF



END SUBROUTINE plot_xanes_quadrupole


SUBROUTINE define_index_arrays(paw_iltonhb)
  USE paw_gipaw,   ONLY : paw_lmaxkb, paw_recon
  USE ions_base,   ONLY : ntyp => nsp
  USE parameters,  ONLY : lmaxx
  USE xspectra_paw_variables, ONLY : xspectra_paw_nhm ! CG
  ! Arguments
  INTEGER paw_iltonhb(0:paw_lmaxkb,xspectra_paw_nhm,ntyp) ! CG

  ! Local
  INTEGER nt,ih,nb,l
  INTEGER ip_per_l(0:lmaxx) 

  DO nt = 1, ntyp
     ih = 1
     ip_per_l(:)=0
     DO nb = 1, paw_recon(nt)%paw_nbeta
        l = paw_recon(nt)%aephi(nb)%label%l
        ip_per_l(l)=ip_per_l(l)+1
        paw_iltonhb(l,ip_per_l(l),nt)= ih
        ih=ih+2*l+1
     ENDDO
  ENDDO


END SUBROUTINE define_index_arrays



FUNCTION lastterm(a,b,g)
  USE kinds, ONLY : dp
  IMPLICIT NONE
  REAL(dp)  :: a,b, g,y1,y2,z1,z2,r
  COMPLEX(dp) :: lastterm

y1 =a*a-g*g-4*b
y2=-2*a*g
r=SQRT(y1*y1+y2*y2)/2

IF (g<0) THEN
z1=a/2+0.5*SIGN(SQRT(y1/2+r),y2)
z2=-g/2+SQRT(-y1/2+r)/2
ELSE
z1=a/2-0.5*SIGN(SQRT(y1/2+r),y2)
z2=-g/2-SQRT(-y1/2+r)/2
ENDIF
 
 lastterm=CMPLX(z1,z2,kind=DP)

END FUNCTION lastterm



FUNCTION paste_fermi(e,ef,a,b,gamma,m,term, first)
  USE kinds, ONLY : dp
  USE xspectra, ONLY : xnitermax, xcheck_conv
  USE cut_valence_green, ONLY :&
       cut_ierror,cut_stepu, cut_stepl,&
       cut_startt,cut_tinf, cut_tsup, cut_nmemu, cut_nmeml, memu, meml
  IMPLICIT NONE

  REAL(dp) :: paste_fermi
  INTEGER  :: m
  REAL(dp) :: a(xnitermax)
  REAL(dp) :: b(xnitermax)
  REAL(dp) :: gamma
  REAL(dp) :: e,ef
  LOGICAL :: term
  COMPLEX(dp) :: green,y,dy,c1,c2, e1, e2
  REAL(dp) ::twopi, t,dt, t1, ta, tb
  INTEGER, save :: n1
  INTEGER, save :: n2
  INTEGER :: nn1,nn2
  LOGICAL :: first

  IF (first) THEN
     memu(:,:)=(0.d0,0.d0)
     meml(:,:)=(0.d0,0.d0) 
     n1=0
     n2=0
     first=.false.
  ENDIF

  dy=cut_ierror+1.0
  y=0.d0
  twopi=6.28318530

  nn1=1
  nn2=1

  t1=0.5773502692

  t=cut_startt

  DO WHILE ((abs(dy)>cut_ierror).OR.(t<cut_tsup))
     dt=cut_stepu*t
     ta=t+dt*(1-t1)/2
     tb=t+dt*(1+t1)/2
     e1=CMPLX(ef,ta,kind=DP)
     e2=CMPLX(ef,tb,kind=DP)

     IF (nn1>n1) THEN
        c1=green(a,b,e1,m,term)
        c2=green(a,b,e2,m,term)
        IF (nn1<cut_nmemu) THEN
           memu(nn1,1)=c1
           memu(nn1,2)=c2
           n1=nn1
        ENDIF
     ELSE
        c1=memu(nn1,1)
        c2=memu(nn1,2)
     ENDIF

     dy=(dt/2)*(c1/CMPLX(ef-e,ta-gamma,kind=DP)+CONJG(c1)/CMPLX(ef-e,-ta-gamma,kind=DP)&
        +c2/CMPLX(ef-e,tb-gamma,kind=DP)+CONJG(c2)/CMPLX(ef-e,-tb-gamma,kind=DP))
     y=y+dy
     t=t+dt
     nn1=nn1+1
  ENDDO

  t=cut_startt
  dy=cut_ierror+1

  DO WHILE((abs(dy)>cut_ierror).OR.(t>cut_tinf))
     dt=cut_stepl*t
     ta=t-dt*(1-t1)/2
     tb=t-dt*(1+t1)/2
     e1=CMPLX(ef,ta,kind=DP)
     e2=CMPLX(ef,tb,kind=DP)

     IF (nn2>n2) THEN
        c1=green(a,b,e1,m,term)
        c2=green(a,b,e2,m,term)
        IF (nn2<cut_nmeml) THEN
           meml(nn2,1)=c1
           meml(nn2,2)=c2
           n2=nn2
        ENDIF
     ELSE
        c1=meml(nn2,1)
        c2=meml(nn2,2)
     ENDIF

     dy=(dt/2)*(c1/CMPLX(ef-e,ta-gamma,kind=DP)+CONJG(c1)/CMPLX(ef-e,-ta-gamma,kind=DP)&
        +c2/CMPLX(ef-e,tb-gamma,kind=DP)+CONJG(c2)/CMPLX(ef-e,-tb-gamma,kind=DP))
     y=y+dy
     t=t-dt
     nn2=nn2+1
  ENDDO

  paste_fermi=AIMAG(y)/twopi

END FUNCTION paste_fermi


FUNCTION green(a,b,e,m, term)
  USE kinds, ONLY : dp
  USE xspectra, ONLY : xnitermax, xcheck_conv
  IMPLICIT NONE

  COMPLEX(dp) :: green
  INTEGER  :: m
  REAL(dp) :: a(xnitermax)
  REAL(dp) :: b(xnitermax)
  COMPLEX(dp) :: e
  LOGICAL :: term
  INTEGER :: i, p,q
  COMPLEX(dp) :: res ,lastterm 
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

     res=lastterm(aa-REAL(e),bb*bb,AIMAG(e))
  ELSE
     res = CMPLX(a(m)-REAL(e),AIMAG(e),kind=DP)
  ENDIF
  DO i = 1, m -1
     res = a(m-i)-e -b(m-i)*b(m-i)/res
  ENDDO

  green = 1/res

END FUNCTION green


SUBROUTINE check_paw_projectors(xiabs)
  USE kinds,           ONLY : DP
  USE constants,       ONLY : pi
  USE paw_gipaw,       ONLY : &
       paw_lmaxkb, &
       paw_recon
  USE xspectra_paw_variables, ONLY : xspectra_paw_nhm
  USE atom,            ONLY : rgrid, msh
  !  USE atom,  ONLY : &
  !       mesh,     &!mesh(ntypx) number of mesh points              
  !       msh ,     &!msh(ntypx)the point at rcut=end of radial integration 
  !       r, rab
  USE ions_base,       ONLY : ntyp => nsp
  USE io_global,       ONLY : stdout
  USE radin_mod
  IMPLICIT NONE
  INTEGER xiabs
  ! internal
  INTEGER :: nr,nrc,ip,jp,lmax,l,ip_l,jtyp,n1,n2,nrs,ndm,ih,jh
  REAL(dp) :: overlap,rexx,overlap2
  REAL (dp), ALLOCATABLE :: aux(:),f(:,:)
  REAL(dp) , ALLOCATABLE :: s(:,:),e(:),v(:,:)

  ALLOCATE(aux(rgrid(xiabs)%mesh)) !allocation too big, it needs only up to msh
  ALLOCATE(f(rgrid(xiabs)%mesh,2)) !allocation too big, it needs only up to msh


  WRITE(stdout,*) '----  PAW projectors from reconstruction files -----'
  WRITE(stdout,*)
  WRITE(stdout,*) 'atom type,  total   number of projectors'
  DO jtyp=1,ntyp
     WRITE (stdout,'(2i4)') jtyp,paw_recon(jtyp)%paw_nbeta
  ENDDO
  WRITE(stdout,*)

  !
  ! I calculate maximum l
  !

  lmax=0
  DO ip=1,paw_recon(xiabs)%paw_nbeta
     IF(paw_recon(xiabs)%psphi(ip)%label%l.GT.lmax) &
          lmax = paw_recon(xiabs)%psphi(ip)%label%l
  ENDDO


  WRITE(stdout,*) 'atom type,  l,   number of projectors per ang. mom.'

  DO jtyp=1,ntyp
     DO l=0,lmax
        WRITE(stdout,'(3i4)') jtyp,l,paw_recon(jtyp)%paw_nl(l)
     ENDDO
  ENDDO



  ! We calculate the overlaps between partial waves and projectors
  ! to see if they are equal to the croneker delta.

  nr=msh(xiabs)  ! extended up to all the NON ZERO points in the mesh.


  WRITE(stdout,*) '----  Overlaps between partial waves and projectors (radial) -----'
  WRITE(stdout,*)
  WRITE(stdout,*) '<tilde{phi} l,n|tilde{p} l,nn>=delta_{n,nn}'
  WRITE(stdout,*) 

  
  DO ip=1,paw_recon(xiabs)%paw_nbeta
     DO jp=1,paw_recon(xiabs)%paw_nbeta
        IF(paw_recon(xiabs)%psphi(ip)%label%l.EQ.paw_recon(xiabs)%psphi(jp)%label%l) THEN
           nrc=Count(rgrid(xiabs)%r(1:nr).le.paw_recon(xiabs)%psphi(ip)%label%rc)
           IF(nrc.GT.nr) THEN
              WRITE(stdout,*) 'nrc=',nrc,' > ',nr,' = nr' 
              CALL errore ( "nrc > nr", "xanes_dipole", 0 )
           ENDIF
           aux(1:nrc)=paw_recon(xiabs)%psphi(ip)%psi(1:nrc)*paw_recon(xiabs)%paw_betar(1:nrc,jp)
           aux(nrc+1:nr)=0.d0
           WRITE(stdout,'("<tilde{phi}_",2i2,10X,"|tilde{p}_",2i2,">=",1f14.8)')  &
                ip,paw_recon(xiabs)%psphi(ip)%label%l,jp, &
                paw_recon(xiabs)%psphi(jp)%label%l, &
                para_radin(aux(1:nr),rgrid(xiabs)%r(1:nr),nr)
        ENDIF
     ENDDO
  ENDDO

  WRITE(stdout,*)

  !
  !
  !

  WRITE(stdout,*) '---- Check normalization pseudo,ae wf and projectors -----------'
  WRITE(stdout,*) '----    (radial part only, integral up to r_c)    -----------'
  WRITE(stdout,*)
  WRITE(stdout,*) 'l,   n, |proj|^2, |pswf|^2 , |aewf|^2'
  WRITE(stdout,*)
  DO l=0,lmax
     DO ip=1,paw_recon(xiabs)%paw_nbeta
        IF(paw_recon(xiabs)%psphi(ip)%label%l.EQ.l) THEN
           nrc=Count(rgrid(xiabs)%r(1:nr).le.paw_recon(xiabs)%psphi(ip)%label%rc)
           aux(1:nrc) = paw_recon(xiabs)%paw_betar(1:nrc,ip) &
                      * paw_recon(xiabs)%paw_betar(1:nrc,ip)
           overlap=para_radin(aux(1:nrc),rgrid(xiabs)%r(1:nrc),nrc)
           aux(1:nrc)=paw_recon(xiabs)%aephi(ip)%psi(1:nrc)*paw_recon(xiabs)%aephi(ip)%psi(1:nrc)
           overlap2=para_radin(aux(1:nrc),rgrid(xiabs)%r(1:nrc),nrc)
           aux(1:nrc)=paw_recon(xiabs)%psphi(ip)%psi(1:nrc)*paw_recon(xiabs)%psphi(ip)%psi(1:nrc)
           rexx=para_radin(aux(1:nrc),rgrid(xiabs)%r(1:nrc),nrc)
           WRITE(stdout,'(2i4,3f14.8)')l,ip,overlap,overlap2,rexx
        ENDIF
     ENDDO
  ENDDO
  WRITE(stdout,*)
  
  goto 323
  WRITE(stdout,*) '---- <phi|chi>= \sum_nl <phi |phi_l><p_l| chi>_nrc  --------'
  WRITE(stdout,*) 'WARNING : this test assumes a form of the phi/chi function'
  
  !  DO l=0,lmax
  DO l=1,1
     ip_l=0
     DO ip=1,paw_recon(xiabs)%paw_nbeta
        IF(ip_l.EQ.0.AND.paw_recon(xiabs)%psphi(ip)%label%l.EQ.l) ip_l=ip
     ENDDO
     
     f(:,:)=0.d0
     DO ip=1,paw_recon(xiabs)%paw_nbeta
        IF(paw_recon(xiabs)%psphi(ip)%label%l.EQ.l) THEN
           f(1:nr,1)=f(1:nr,1)+paw_recon(xiabs)%psphi(ip)%psi(1:nr)/REAL(ip,dp)
           f(1:nr,2)=f(1:nr,2)+1.123*paw_recon(xiabs)%psphi(ip)%psi(1:nr)/REAL(ip,dp)
        ENDIF
     ENDDO
     rexx=0.d0
     DO ip=1,paw_recon(xiabs)%paw_nbeta
        IF(paw_recon(xiabs)%psphi(ip)%label%l.EQ.l) THEN
           nrc=Count(rgrid(xiabs)%r(1:nr).le.paw_recon(xiabs)%psphi(ip)%label%rc)
           IF(nrc.GT.nr) THEN
              WRITE(stdout,*) 'nrc=',nrc,' > ',nr,' = nr' 
              CALL errore ( "nrc > nr", "xanes_dipole", 0 )
           ENDIF
           aux(1:nrc)=f(1:nrc,1)*paw_recon(xiabs)%paw_betar(1:nrc,ip)
           overlap=para_radin(aux(1:nrc),rgrid(xiabs)%r(1:nrc),nrc)
           aux(1:nrc)=f(1:nrc,2)*paw_recon(xiabs)%psphi(ip)%psi(1:nrc)
           overlap=overlap*para_radin(aux(1:nrc),rgrid(xiabs)%r(1:nrc),nrc)
           WRITE(stdout,'("overlap(l=",1i2,",n=",1i2,")= ",1f14.8)')l,ip,overlap
           rexx=rexx+overlap
        ENDIF
     ENDDO
     aux(1:nr)=f(1:nr,1)*f(1:nr,2)
     WRITE(stdout,'("sum/overlap=",1f14.8)') rexx/para_radin(aux,rgrid(xiabs)%r(1:nr),nrc)
     WRITE(stdout,'("sum projectors=",1f14.8," overlap=",1f14.8)') rexx,para_radin(aux,rgrid(xiabs)%r(1:nr),nrc)
     WRITE(stdout,*)'---+++++----'
  ENDDO
  !           aux(1:nrc)=paw_recon(xiabs)%psphi(ip)%psi(1:nrc)*paw_recon(xiabs)%paw_betar(1:nrc,jp)

  !           
  !
  !
  !        ENDIF
  !          ENDDO
  WRITE(stdout,*)
  WRITE(stdout,*) '================================================================'
  
323 continue
  !
  !  Check linear dependence of projectors
  !

  WRITE(stdout,*) '================================================================'
  WRITE(stdout,*) '           Checking linear dipendence of projectors             '
  WRITE(stdout,*)

  DEALLOCATE(aux)

  ndm = MAXVAL (msh(1:ntyp))

  ALLOCATE(aux(ndm))

  DO l=0,paw_lmaxkb
     IF (paw_recon(xiabs)%paw_nl(l)>0) THEN
        ALLOCATE (s(paw_recon(xiabs)%paw_nl(l),paw_recon(xiabs)%paw_nl(l)))
        ALLOCATE (e(paw_recon(xiabs)%paw_nl(l)),v(paw_recon(xiabs)%paw_nl(l),paw_recon(xiabs)%paw_nl(l)))
        DO ih=1,paw_recon(xiabs)%paw_nl(l)
           n1=paw_recon(xiabs)%paw_iltonh(l,ih)
           nrc=paw_recon(xiabs)%psphi(n1)%label%nrc
           nrs=paw_recon(xiabs)%psphi(n1)%label%nrs
           DO jh=1,paw_recon(xiabs)%paw_nl(l)
              n2=paw_recon(xiabs)%paw_iltonh(l,jh)
              CALL step_f(aux,paw_recon(xiabs)%psphi(n1)%psi(1:msh(xiabs)) * &
                   paw_recon(xiabs)%psphi(n2)%psi(1:msh(xiabs)), &
                   rgrid(xiabs)%r(1:msh(xiabs)),nrs,nrc, 1.d0, msh(xiabs) )
              CALL simpson ( msh(xiabs), aux, rgrid(xiabs)%rab(1), s(ih,jh))
           ENDDO
        ENDDO
  
        WRITE(stdout,'("atom type=",1i4)') xiabs
        WRITE(stdout,'("number of projectors projector  =",1i4," angular momentum=",1i4)') &
             paw_recon(xiabs)%paw_nl(l),l
        DO ih=1,paw_recon(xiabs)%paw_nl(l)
           WRITE(stdout,'(10f14.8)') (s(ih,jh),jh=1,paw_recon(xiabs)%paw_nl(l)) 
        ENDDO
        WRITE(stdout,*) 'Eigenvalues S matrix:'
        WRITE(stdout,*)

        IF(paw_recon(xiabs)%paw_nl(l).EQ.1) THEN
           WRITE(stdout,'(1i4,1f14.8)') 1,s(1,1)
        ELSE 
           CALL rdiagh(paw_recon(xiabs)%paw_nl(l),s,paw_recon(xiabs)%paw_nl(l) , e, v )
           DO ih=1,paw_recon(xiabs)%paw_nl(l)
              WRITE(stdout,'(1i4,1f14.8)') ih,e(ih)
           ENDDO
        ENDIF
        WRITE(stdout,*)
        DEALLOCATE(s,e,v)
     ENDIF
  ENDDO
  WRITE(stdout,*) '================================================================'


  DEALLOCATE(aux,f)
END SUBROUTINE check_paw_projectors


SUBROUTINE read_save_file(a,b,xnorm,ncalcv,x_save_file,core_energy)
  USE kinds,       ONLY : DP
  USE klist,       ONLY : nks,nkstot
  USE xspectra,    ONLY : xnitermax,xang_mom,xiabs,&
       n_lanczos, save_file_version, save_file_kind, calculated     
  USE io_global,   ONLY : stdout,ionode
  USE lsda_mod,    ONLY : nspin,lsda
  IMPLICIT NONE
  CHARACTER(LEN=256) :: x_save_file
  INTEGER            :: ierr,nkstot_r
  INTEGER            :: xm_r,nc_r,ncomp_max
  INTEGER ncalcv(n_lanczos,nks)
  REAL(dp) core_energy
  REAL(dp) a(xnitermax,n_lanczos,nks),b(xnitermax,n_lanczos,nks)     
  REAL (dp)  xnorm(n_lanczos,nks)
  INTEGER, ALLOCATABLE :: ncalcv_all(:,:)
  REAL(dp), ALLOCATABLE :: a_all(:,:),b_all(:,:),xnorm_all(:,:),aux(:,:)
  REAL(dp) xkvec_r(3)
  REAL(dp) :: xepsilon_r(3)
  INTEGER i,j,k,ncalcv_max
  INTEGER :: calculated_all(n_lanczos,nkstot)

  ALLOCATE(a_all(xnitermax,nkstot))
  ALLOCATE(b_all(xnitermax,nkstot))
  ALLOCATE(xnorm_all(n_lanczos,nkstot))
  ALLOCATE(ncalcv_all(n_lanczos,nkstot))
  ALLOCATE(aux(xnitermax,nks))

  a(:,:,:)=0.d0
  b(:,:,:)=0.d0
  xnorm(:,:)=0.d0
  ncalcv(:,:)=0
  calculated_all(:,:)=0


  OPEN ( UNIT = 10, FILE = x_save_file, FORM = 'FORMATTED', &
       STATUS = 'UNKNOWN', IOSTAT = ierr )
  CALL errore( 'iosys', 'x_save_file ' // TRIM( x_save_file ) // &
       & ' not found' , ierr )
  REWIND(10)

  IF (save_file_version.eq.0) then
     WRITE(stdout,*) 'save file version : old'
  ELSEIF(save_file_version.eq.1) then
     WRITE(stdout,*) 'save file version : 1'
     DO i=1,6 
        READ(10,*)      
     ENDDO
  ENDIF


  READ(10,*) lsda,nspin
  READ(10,*) xm_r,nkstot_r,xnitermax

  IF(xm_r.NE.xang_mom) THEN
     WRITE(stdout,*) 'xm_r=',xm_r
     CALL errore('read_save_file','xm_r is different from xang_mom=',xang_mom)
  ENDIF

  READ(10,*) ncalcv_max
  IF(ncalcv_max.GT.xnitermax) THEN
     WRITE(stdout,*) 'ncalcv_max=',ncalcv_max
     CALL errore('read_save_file','ncalcv_max is grater than xnitermax=',xnitermax)
  ENDIF

  READ(10,*) core_energy

  READ(10,*) (xkvec_r(i),i=1,3)
  WRITE(stdout,*) '---------------------------------------------------------'
  WRITE(stdout,*) 'xkvec read from savefile'
  WRITE(stdout,*) (xkvec_r(i),i=1,3)
  READ(10,*) (xepsilon_r(i),i=1,3)
  WRITE(stdout,*) 'xepsilon read from file'
  WRITE(stdout,*) (xepsilon_r(i),i=1,3)
  WRITE(stdout,*) '---------------------------------------------------------'
  write(stdout,*) 'n_lanczos=', n_lanczos
  DO i=1, n_lanczos
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
  close(10)

#ifdef __MPI 
  CALL poolscatter(n_lanczos,nkstot,xnorm_all,nks,xnorm)
  CALL ipoolscatter(n_lanczos,nkstot,ncalcv_all,nks,ncalcv)
  CALL ipoolscatter(n_lanczos,nkstot,calculated_all,nks,calculated)
#else
  IF(nks.NE.nkstot) THEN
     CALL errore('read_save_file','nks\=nkstot',1)
  ENDIF
  xnorm(1:n_lanczos,1:nkstot)=xnorm_all(1:n_lanczos,1:nkstot)
  ncalcv(1:n_lanczos,1:nkstot)=ncalcv_all(1:n_lanczos,1:nkstot)
  calculated(1:n_lanczos,1:nkstot)=calculated_all(1:n_lanczos,1:nkstot)
#endif
  DEALLOCATE(a_all)
  DEALLOCATE(b_all)
  DEALLOCATE(xnorm_all)
  DEALLOCATE(ncalcv_all)
  DEALLOCATE(aux)


END SUBROUTINE read_save_file

SUBROUTINE read_header_save_file(x_save_file)
  USE kinds, ONLY : DP
  USE klist,      ONLY : nkstot
  USE lsda_mod,    ONLY : nspin,lsda
  USE xspectra, ONLY : save_file_version, save_file_kind, n_lanczos
  USE io_global, ONLY : stdout
  IMPLICIT NONE
  CHARACTER(LEN=256) :: x_save_file
  INTEGER            :: ierr,nkstot_r
  INTEGER            :: xm_r
  CHARACTER          :: c

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


SUBROUTINE write_save_file(a,b,xnorm,ncalcv,x_save_file)
  USE kinds, ONLY : DP
  USE klist,      ONLY : nks,nkstot
  USE xspectra,      ONLY : xnitermax,xang_mom,xkvec,xepsilon,xiabs,&
       save_file_version, save_file_kind, n_lanczos, calculated
  USE io_global,       ONLY : ionode
  !*apsi  USE uspp_param, ONLY : psd
  USE lsda_mod,    ONLY : nspin,lsda
  USE uspp_param, ONLY : upf
  IMPLICIT NONE
  CHARACTER(LEN=256) :: x_save_file
  INTEGER            :: ierr
  INTEGER ncalcv(n_lanczos,nks)
  REAL(dp) a(xnitermax,n_lanczos,nks),b(xnitermax,n_lanczos,nks)     
  REAL (dp)  xnorm(n_lanczos,nks)
  INTEGER, ALLOCATABLE :: ncalcv_all(:,:)
  REAL(dp), ALLOCATABLE :: a_all(:,:),b_all(:,:),xnorm_all(:,:)
  REAL (dp) :: mygetK
  INTEGER i,j,k,ncalcv_max
  CHARACTER(LEN=8) :: dte
  INTEGER :: calculated_all(n_lanczos,nkstot)

  IF (ionode)  CALL DATE_AND_TIME(date=dte)

  ALLOCATE(a_all(xnitermax,nkstot))
  ALLOCATE(b_all(xnitermax,nkstot))
  ALLOCATE(xnorm_all(n_lanczos,nkstot))
  ALLOCATE(ncalcv_all(n_lanczos,nkstot))

  ncalcv_all(:,:)=0
  xnorm_all(:,:)=0.d0

  ncalcv_all(1:n_lanczos,1:nks)=ncalcv(1:n_lanczos,1:nks)
  xnorm_all(1:n_lanczos,1:nks)=xnorm(1:n_lanczos,1:nks)
  calculated_all(1:n_lanczos,1:nks)=calculated(1:n_lanczos,1:nks)

#ifdef __MPI
  CALL poolrecover(xnorm_all,n_lanczos,nkstot,nks)
  CALL ipoolrecover(ncalcv_all,n_lanczos,nkstot,nks)
  CALL ipoolrecover(calculated_all,n_lanczos,nkstot,nks)
#endif


  ncalcv_max=0
  DO j=1, n_lanczos
     DO i=1,nkstot
        IF(ncalcv_all(j,i).GT.ncalcv_max) ncalcv_max=ncalcv_all(j,i)
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
     WRITE(10,*) mygetK(upf(xiabs)%psd)
     WRITE(10,*) (xkvec(i),i=1,3)
     WRITE(10,*) (xepsilon(i),i=1,3)
  ENDIF
  DO i=1, n_lanczos
     a_all(:,:)=0.d0
     b_all(:,:)=0.d0
     a_all(1:xnitermax,1:nks)=  a(1:xnitermax,i,1:nks)
     b_all(1:xnitermax,1:nks)=  b(1:xnitermax,i,1:nks)
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
  close(10)

  DEALLOCATE(a_all)
  DEALLOCATE(b_all)
  DEALLOCATE(xnorm_all)
  DEALLOCATE(ncalcv_all)


END SUBROUTINE write_save_file



SUBROUTINE write_status_of_the_code
  USE io_global,       ONLY : stdout
  IMPLICIT NONE

  WRITE (stdout,*) '======= Working features (22/04/2009) ==========='
  WRITE (stdout,*) 'xanes works both in dipolar and quadrupolar part,'
  WRITE (stdout,*) 'spin polarized works'
  WRITE (stdout,*) 'DFT+U implemented, validated'
  WRITE (stdout,*) 'ultrasoft pseudo works'
  WRITE (stdout,*) 'cut occupied states working, improved'
  WRITE (stdout,*) 'terminator working'
  WRITE (stdout,*) 'Multiprojectors TM+USPP working (MCB,CG)'
  WRITE (stdout,*) 'new save file format, with version numbering'
  WRITE (stdout,*) 'time limit implemented, with restart, seems to work'
  WRITE (stdout,*) 'DFT+U tested ONLY for non ortho wfc, but implemented'
  WRITE (stdout,*) '======= TO DO                         ==========='
  WRITE (stdout,*) 'Bethe-Salpeter [CG] '
  WRITE (stdout,*) 'RXES [DC] o [CG]' 
  WRITE (stdout,*) '================================================='

END SUBROUTINE write_status_of_the_code




SUBROUTINE read_gamma_file
  USE gamma_variable_mod
  IMPLICIT NONE
  INTEGER :: nl, ierr, i

  open (unit=21, file=gamma_file, form='formatted',status='unknown', iostat=ierr)
  CALL errore ('io ', 'gamma file '//TRIM(gamma_file)//' not found', abs (ierr) )
  rewind(21)

  nl=0

  DO
     READ (21,'(a1)',iostat=ierr)
     IF (ierr.NE.0) EXIT
     nl=nl+1
  ENDDO
  close(21)

  gamma_lines=nl
  ALLOCATE(gamma_points(nl,2))

  open (unit=21, file=gamma_file, form='formatted',status='unknown', iostat=ierr)
  rewind(21)

  DO i=1,nl
     READ(21,*) gamma_points(i,1), gamma_points(i,2)
  ENDDO

  close(21)


END SUBROUTINE read_gamma_file

SUBROUTINE initialize_gamma_tab
  USE xspectra, ONLY : xemin, xemax, xnepoint
  USE kinds, ONLY :dp
  USE io_global, ONLY : stdout
  USE constants, ONLY : rytoev
  USE gamma_variable_mod
  IMPLICIT NONE
  REAL(dp) :: e,x,y,dx
  INTEGER :: i,j,n
  dx=(xemax-xemin)/dfloat(xnepoint)

  DO n=1, xnepoint
     x=xemin+(n-1)*dx
     i=1
     DO j=1, gamma_lines
        IF(x>gamma_points(j,1)) i=i+1
     ENDDO

     IF (i.EQ.1) THEN
        y=gamma_points(1,2)
     ELSEIF (i.EQ.(gamma_lines+1)) THEN
        y=gamma_points(gamma_lines,2)
     ELSE
        y=(gamma_points(i-1,2)*(gamma_points(i,1)-x)+gamma_points(i,2)*(x-gamma_points(i-1,1)))&
             /(gamma_points(i,1)-gamma_points(i-1,1))
     ENDIF
     gamma_tab(n)=y/rytoev
  ENDDO

END SUBROUTINE initialize_gamma_tab

SUBROUTINE determine_polycut(t1,t2,f1,f2,df1,df2,poly)
  ! calculates the interpolation polynome betwenn 2 points
  USE kinds, ONLY : dp
  IMPLICIT NONE
  REAL(dp) :: t1,t2,f1,f2,df1,df2,poly(4)

  poly(4)=((t2-t1)*(df2+df1)-2*(f2-f1))/((t2-t1)**3)
  poly(3)=(df2-df1)/(2*(t2-t1))-1.5d0*(t2+t1)*((t2-t1)*(df2+df1)-2*(f2-f1))/((t2-t1)**3)
  poly(2)=df1-2*t1*poly(3)-3*t1*t1*poly(4)
  poly(1)=f1-poly(2)*t1-poly(3)*t1**2-poly(4)*t1**3

END SUBROUTINE determine_polycut

!</CG>

