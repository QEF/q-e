PROGRAM Xxanes
#include "f_defs.h"
  USE kinds, only : DP
  USE constants,          ONLY : degspin
  USE io_global,       ONLY : stdout,ionode,ionode_id   ! Modules/io_global.f90
  USE io_files,        ONLY : nd_nmbr, prefix, tmp_dir
  USE parser,          ONLY :  read_line
  USE cell_base,       ONLY : bg, at, celldm
  USE global_version,  ONLY : version_number
  use parameters,       only : ntypx,lmaxx,npsx,lqmax
  USE ions_base,          ONLY : nat, ntyp => nsp, ityp, tau
  USE ktetra,             ONLY : nk1, nk2, nk3, k1, k2, k3, &
       ltetra, ntetra, tetra
  USE control_flags,    ONLY : gamma_only
  use wvfct,            ONLY : npwx,nbnd,npw,igk,et! et(nbnd,nkstot)
  USE radial_grids,     ONLY : ndmx
  USE atom,             ONLY : rgrid
  use becmod, ONLY:becp,rbecp
  USE uspp,   ONLY : vkb, nkb, okvan 
  USE xanes
  USE ener, ONLY : ef, ef_up, ef_dw !Fermi energy
  USE symme,   ONLY : nsym,s
  use paw_gipaw,              only : read_recon,  &
       paw_vkb,             & ! |p> projectors
       paw_becp,            & ! product of projectors and wf.
       paw_nkb,             & ! total number of beta functions, with st.fact.
       paw_lmaxkb,          &
       paw_recon

  use klist,       ONLY : &
       nkstot,            & ! total number of k-points
       nks,               & ! number of k-points for local pool
       nelec,nelup,neldw,             & !number of electrons
       xk,                & ! k-points coordinates
       wk ,               & ! k-points weight
       npk,               &
       degauss,lgauss,ngauss,    &
       two_fermi_energies

  use lsda_mod,    ONLY : nspin,lsda,isk,current_spin
  USE noncollin_module,     ONLY : noncolin
  use mp,         only : mp_bcast, mp_sum             !parallelization
  USE mp_global,  ONLY : intra_pool_comm, nproc, npool 

  use cut_valence_green, only :&
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

  USE control_flags, only : twfcollect
  !<CG>
  USE gamma_variable_mod, only : gamma_lines, gamma_tab, gamma_points, gamma_mode, gamma_file
  use xanes_paw_variables, only : xanes_paw_nhm, init_xanes_paw_nhm
  !</CG>

  implicit none 
  !
  ! ... local variables
  !
  real(kind=dp) ryd2eV, pi
  parameter(ryd2ev = 13.6058d0, pi = 3.141592653589793d0)
  logical use_paratec_recon, terminator, show_status, wf_collect
  integer :: nargs,iiarg,ierr,ios,il,ibnd,ibnd_up,ibnd_dw,xm_r,nc_r,ncomp_max
  integer, external :: iargc
  integer nt,nb,na,i,j,k,nrest,nline
  integer, allocatable :: ncalcv(:,:)
  INTEGER, ALLOCATABLE ::&
       paw_iltonhb(:,:,:)      ! corresp l, projector, type <--> cumulative over all the species
  real (DP) ehomo, elumo,norm, core_energy
  real(kind=DP) :: core_wfn(ndmx)
  real(dp), allocatable:: a(:,:,:),b(:,:,:),xnorm(:,:)      !lanczos vectors
  REAL (DP), EXTERNAL :: efermig,efermit
  real(DP) :: rc(ntypx,0:lmaxx),r_paw(0:lmaxx)
  LOGICAL :: vloc_set

  CHARACTER (LEN=256)  :: input_file, filerecon(ntypx),filecore
  character(len=256) :: outdir 
  character(len=25) :: calculation
  character(len=4) :: verbosity
  character(len=10) :: dummy_char
  real(dp) :: gamma_energy(2), gamma_value(2)


  ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  ! $   Namelists Definition
  ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$


  namelist / input_xspectra / &
       calculation,&       ! calulation type : xanes...
       verbosity, &        ! high/low
       prefix, &           ! prefix of the pwscf output files
       outdir, &           ! directory tmp_dir or where the files are
       xang_mom,&
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
       U_projection_type   

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
       r_paw, &
       use_paratec_recon

  namelist / cut_occ /&
       cut_ierror, cut_stepu, cut_stepl, cut_startt, cut_tinf, cut_tsup,&
       cut_desmooth, cut_nmemu, cut_nmeml


  ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  !    Initialising post processing (This initialise MPI, it has to be left here)
  ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  CALL start_postproc(nd_nmbr)

  ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  ! $   Default values for namelists
  ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$


  !
  ! Set default values for namelist: (fforgot what is kedge)
  !

  calculation='xanes'
  prefix=''
  verbosity='low'
  x_save_file='xanes.sav'
  outdir='./'
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
  do i=2,3
     xkvec(i)=0.d0
     xepsilon(i)=0.d0
  enddo
  xkvec(1)=1.d0
  xepsilon(1)=1.d0
  xcoordcrys=.true.
  ef_r=0.d0
  use_paratec_recon=.FALSE.
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


  ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  ! $   Check if the input is from file or from stdin
  ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  !
  ! ... Input from file ?
  !
  if ( ionode ) then

     nargs = iargc()

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


     read(5, input_xspectra, err = 200, iostat = ios) 
200  call errore ('input_xspectra', 'reading input_xspectra namelist', abs (ios) )

     read(5, plot, err = 300, iostat = ios) 
300  call errore ('plot', 'reading plot namelist', abs (ios) )

     read(5, pseudos, err = 400, iostat = ios) 
400  call errore ('pseudos', 'reading pseudos namelist', abs (ios) )

     read(5, cut_occ, err = 500, iostat = ios) 
500  call errore ('cut_occ', 'reading cut_occ namelist', abs (ios) )


     tmp_dir = trim(outdir) 

  end if


  ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  ! $   Variables broadcasting 
  ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  CALL mp_bcast( calculation, ionode_id)
  CALL mp_bcast( tmp_dir, ionode_id ) 
  CALL mp_bcast( prefix,  ionode_id )
  CALL mp_bcast( outdir,  ionode_id ) 
  CALL mp_bcast( xnepoint,  ionode_id ) 
  CALL mp_bcast( xniter,  ionode_id ) 
  CALL mp_bcast( xcheck_conv,  ionode_id ) 
  CALL mp_bcast( xang_mom,  ionode_id ) 
  CALL mp_bcast( xgamma,  ionode_id ) 
  CALL mp_bcast( xerror,  ionode_id ) 
  CALL mp_bcast( xemin,  ionode_id ) 
  CALL mp_bcast( xemax,  ionode_id ) 
  CALL mp_bcast( show_status, ionode_id)

  CALL mp_bcast( xkvec,  ionode_id ) 
  CALL mp_bcast( xepsilon,  ionode_id ) 

  CALL mp_bcast( xonly_plot,  ionode_id ) 
  CALL mp_bcast( filerecon,  ionode_id ) 
  CALL mp_bcast( filecore,  ionode_id ) 
  CALL mp_bcast( xiabs,  ionode_id ) 
  CALL mp_bcast( r_paw,  ionode_id ) 
  CALL mp_bcast( xread_wf,  ionode_id ) 
  CALL mp_bcast( x_save_file,  ionode_id ) 
  CALL mp_bcast( xcoordcrys,  ionode_id ) 
  CALL mp_bcast( ef_r,  ionode_id )   
  CALL mp_bcast( use_paratec_recon, ionode_id )   
  CALL mp_bcast( cut_occ_states, ionode_id )
  CALL mp_bcast( terminator, ionode_id )
  CALL mp_bcast( wf_collect, ionode_id )
  CALL mp_bcast( twfcollect, ionode_id )

  CALL mp_bcast( U_projection_type, ionode_id )

  CALL mp_bcast( gamma_mode, ionode_id )
  CALL mp_bcast( gamma_energy, ionode_id )
  CALL mp_bcast( gamma_value, ionode_id )

  CALL mp_bcast( cut_ierror, ionode_id )
  CALL mp_bcast( cut_stepu, ionode_id )
  CALL mp_bcast( cut_stepl, ionode_id )
  CALL mp_bcast( cut_startt, ionode_id )
  CALL mp_bcast( cut_tinf, ionode_id )
  CALL mp_bcast( cut_tsup, ionode_id )
  CALL mp_bcast( cut_desmooth, ionode_id )
  CALL mp_bcast( cut_nmemu, ionode_id )
  CALL mp_bcast( cut_nmeml, ionode_id )

  ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  ! $ Writing the status of the code (working features and things to do)
  ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  if(show_status) call write_status_of_the_code

  ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  !    Initialising clocks
  ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  !
  !
  ! ... use ".FALSE." to disable all clocks except the total cpu time clock
  ! ... use ".TRUE."  to enable clocks
  !

  if(trim(adjustl(calculation)).eq.'xanes_dipole') then
     xang_mom=1                       !so it is not necessary to specify xang_mom
     calculation='xanes'
  elseif(trim(adjustl(calculation)).eq.'xanes_quadrupole') then
     xang_mom=2                       !so it is not necessary to specify xang_mom
     calculation='xanes'
  endif


  CALL init_clocks( .TRUE. )
  CALL start_clock( calculation  )
  !CALL stop_clock( code )


  ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  !  check on wfcollect
  ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  if(xread_wf.and.wf_collect) then
     call errore ('main','incompatibility xread_wf and wf_collect',1)
  endif


  twfcollect=wf_collect


  ! $$$$$$$$$  notstart if on  xonlyplot  $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  if(.not.xonly_plot) then

     ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
     !    read pwscf structural and k-points infos, also ditributes across the pools
     ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

     call read_file_xanes(xread_wf)

     write(stdout,*) 'k-points : nkstot=', nkstot

     ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
     !    write out crystal structure, kpoints list etc.
     ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$


     write(stdout,*) '-------------- Crystal Structure ------------ '
     write(stdout,*) 'celldm(1:6)'
     write(stdout,'(3f14.8)') (celldm(i),i=1,3)
     write(stdout,'(3f14.8)') (celldm(i),i=4,6)
     write(stdout,*) 'direct lattice vectors'
     do i=1,3
        write(stdout,'(3f14.8)') (at(i,j),j=1,3)
     enddo
     write(stdout,*) 'reciprocal lattice vectors'
     do i=1,3
        write(stdout,'(3f14.8)') (bg(i,j),j=1,3)
     enddo

     write(stdout,*) 'nks=',nks,' nkstot=',nkstot 
     write(stdout,*) ' ----k-point list [units 2*pi/celldm(1)], weight---------'

     do i=1,nkstot
        write(stdout,'(1i6,4f14.8)') i,(xk(j,i) , j=1,3),wk(i)
     enddo
     write(stdout,*) '-------------------------------------------------'     

     ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$  
     ! normalize xkvec and xepsilon
     ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$  

     if(xcoordcrys) call cryst_to_cart(1,xepsilon,at,1)
     if(xang_mom.eq.2) then
        if(xcoordcrys) call cryst_to_cart(1,xkvec,at,1)
        norm=dsqrt(xkvec(1)**2+xkvec(2)**2+xkvec(3)**2)
        do i=1,3
           xkvec(i)=xkvec(i)/norm
        enddo
     endif
     norm=dsqrt(xepsilon(1)**2+xepsilon(2)**2+xepsilon(3)**2)
     do i=1,3
        xepsilon(i)=xepsilon(i)/norm
     enddo

     ! check orthogonality

     write (stdout,*) '--- Polarisation and k vector [cartesian coordinates]----'
     write (stdout,'(a,1x,3(f10.8, 1x))') 'xepsilon(:)=', (xepsilon(i),i=1,3)
     write (stdout,'(a,1x,3(f10.8, 1x))') 'xkvec(:)=', (xkvec(i),i=1,3)
     if(xang_mom.eq.2) then 
        if ((abs(xkvec(1)*xepsilon(1)+xkvec(2)*xepsilon(2)+xkvec(3)*xepsilon(3))).ge.1.0d-6) then
           write(stdout,*) 'WARNING, xkvec and xepsilon are not orthogonal'
           write(stdout,*) 'scalar product=',xkvec(1)*xepsilon(1)+xkvec(2)*xepsilon(2)+xkvec(3)*xepsilon(3)
        endif
     endif


     ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
     !  is the type associated to xiabs existing ?
     ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$  


     i=0
     do na=1,nat
        if(ityp(na).eq.xiabs) i=i+1
     enddo
     if(i.ne.1) then
        call errore( 'main program', 'Wrong xiabs!!!',i)
     endif

     ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
     !  Reads reconstruction files
     ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$  

     if(use_paratec_recon) then
        write(stdout,*) 'Warning : Paratec recon not implemented'
        call stop_pp
        call read_core_abs_paratec(filecore,core_wfn,xiabs)
        !   IF ( .NOT. paw_recon(xiabs)%gipaw_data_in_upf_file ) &
        !        call read_recon_paratec ( filerecon(xiabs), xiabs, paw_recon(xiabs), &
        !             vloc_set ) !*apsi
     else
        call read_core_abs(filecore,core_wfn)
        IF ( .NOT. paw_recon(xiabs)%gipaw_data_in_upf_file ) &
             call read_recon ( filerecon(xiabs), xiabs, paw_recon(xiabs) ) !*apsi
     endif

     ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
     !  Reads potentials and so on from post processing
     ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$  

     if(.not.twfcollect) call openfil_pp

     ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
     !  Assign paw radii to species (this will become soon obsolete)
     ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$  
     !
     !  read recon should be parallelized
     !

     do nt=1,ntyp
        DO  j=1,paw_recon(nt)%paw_nbeta
           il=paw_recon(nt)%psphi(j)%label%l
           if(xiabs.eq.nt.and.dabs(r_paw(il)).lt.1.d-6) then
              !*apsi              if(r(kkbeta(nt),nt).gt.1.d-3) then
              !*apsi                rc(nt,il)=  r(kkbeta(nt),nt)
              !*apsi              else
              !*apsi                 write(stdout,*) 'Warning-no-kkbeta r_paw(',il,')=1.0'
              !*apsi                 rc(nt,il)=1.0
              !*apsi              endif
              !<CG>  to be verified
              if (paw_recon(nt)%psphi(il)%label%rc > 1.d-3) then
                 write(stdout,*) 'warning, r_paw(', il,' ) set to ', paw_recon(nt)%psphi(il)%label%rc
                 rc(nt, il)= paw_recon(nt)%psphi(il)%label%rc
              else
                 write(stdout,*) 'Warning, no rc'
                 write(stdout,*) 'warning, r_paw(', il,' ) set to 1.0'
                 rc(nt, il)= 1.d0
              endif
              !</CG>

           elseif(xiabs.eq.nt.and.dabs(r_paw(il)).gt.1.d-6) then
              rc(nt,il)=r_paw(il)
           elseif(nt.ne.xiabs) then
              !*apsi              if(r(kkbeta(nt),nt).gt.1.d-3) then
              !*apsi                 rc(nt,il)=r(kkbeta(nt),nt)
              !*apsi              else
              !*apsi                 rc(nt,il)=1.0
              !*apsi              endif
              !<CG> to be verified
              if(paw_recon(nt)%psphi(il)%label%rc.gt.1.d-3) then
                 rc(nt,il)=paw_recon(nt)%psphi(il)%label%rc
              else
                 rc(nt,il)=1.0
              endif
              !<CG>
           endif
        enddo
     enddo

     !<CG>
     do nt=1,ntyp
        do il = 1,paw_recon(nt)%paw_nbeta
           paw_recon(nt)%psphi(il)%label%rc = rc(nt,paw_recon(nt)%psphi(il)%label%l)
           paw_recon(nt)%aephi(il)%label%rc = rc(nt,paw_recon(nt)%aephi(il)%label%l)
        enddo
     enddo
     !</CG>


     ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
     !  write band energies if xread_wf=true
     ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$  


     if(xread_wf.and.trim(verbosity).eq.'high') then

        write(stdout,*) '------- Band energies read from file -----------'

        do i=1,nkstot
           write(stdout,'("k=[",3f14.8,"]   spin=",1i2)') &
                xk(1,i),xk(2,i),xk(3,i),isk(i)
           write(stdout, '(8f9.4)') (et(j,i)*13.605,j=1,nbnd)
        enddo


        write(stdout,*) '------------------------------------------------'
     endif


     ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
     !  Calculate the Fermi level if possible
     ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$  

     write(stdout,*)
     write(stdout,*) '=========================================================='

     if(trim(calculation).eq.'fermi_level'.and..not.xread_wf) then
        write(stdout,*) 'Impossible to calculate the Fermi level'
        write(stdout,*) '   without reading wavefunctions (xread_wf=.true.)'
        CALL errore( 'main program', &
             'xread_wf incompatible with calculation type', 1 )
     endif


     if(.not.xread_wf) then        !Fermi level must be read from input
        ef=ef_r
        write( stdout,*)
        write( stdout,'(">> Fermi level from input (Ryd)=",1f14.8)') ef
     endif

     if(trim(calculation).eq.'fermi_level') then
        if(ltetra.or.lgauss) then
           write( stdout,*)
           write( stdout,*) 'Metallic case'
           write( stdout,*)
           call errore('input','Read fermi level from scf/nscf output',1)
        else
           write(stdout,*)
           write(stdout,*)  'Insulating case:'
           write(stdout,*)
           !   SPINLESS CASE
           IF ( nspin.eq.1) THEN
              ibnd =  nint (nelec) / 2
              ef = MAXVAL ( et( ibnd  , 1:nkstot) )
              WRITE( stdout, '(">> Fermi level (Ryd) = ",1f14.8)') ef
              !SPIN POLARIZED CALCULATION
           ELSEIF(nspin.eq.2) then
              ibnd    = NINT( nelec )
              if(nelup.eq.0.or.nelup.eq.0) then
                 write(stdout,*) 'WARNING, nelup=0 or neldw=0'
              endif
              write(stdout,*) 'nel=',nelec,' nelup=',nelup, 'neldw=',neldw

              ibnd_up = NINT( nelup )
              ibnd_dw = NINT( neldw )
              IF ( ibnd_up == 0 ) THEN
                 !
                 ef = MAXVAL( et(ibnd_dw,1:nkstot/2) )
                 ef_dw = ef
                 write( stdout,&
                      '(">> Fermi level down (Ryd)= ",1f14.8)')&
                      ef_dw
                 !
              ELSE IF ( ibnd_dw == 0 ) THEN
                 !
                 ef = MAXVAL( et(ibnd_up,1:nkstot/2) )
                 ef_up = ef
                 write( stdout,&
                      '(">> Fermi level up (Ryd)= ",1f14.8)')&
                      ef_up
                 !
              ELSE
                 !
                 write(stdout,*) 'nkstot=',nkstot
                 ef    = MAX( MAXVAL( et(ibnd_up,1:nkstot/2) ), &
                      MAXVAL( et(ibnd_dw,nkstot/2+1:nkstot) ) )
                 ef_up =  MAXVAL( et(ibnd_up,1:nkstot/2) )
                 ef_dw =  MAXVAL( et(ibnd_dw,nkstot/2+1:nkstot) )
                 write( stdout,&
                      '(">> Fermi level up (Ryd)= ",1f14.8,&
                      &" Fermi level down (Ryd)= ",1f14.8)') &
                      ef_up,ef_dw
                 !
              END IF
              write( stdout,'(">> Fermi level (Ryd)=",1f14.8)') ef
           END IF

        ENDIF
        write( stdout,*)
        write (stdout,*) &
             '============================================================'
     endif
     !        call mp_bcast( ef, ionode_id )  !Why should I need this ?





     ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
     !  Allocation of variables for paw projectors
     ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$  

! CG : becp allocated in the xanes_dipole and xanes_quadrupole subroutines
!     if( gamma_only ) then
!        ALLOCATE( rbecp( nkb, nbnd ) )
!     else
!        ALLOCATE( becp( nkb, nbnd ) )
!     endif

     ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
     !  Initialise Vanderbilt and Paw projectors
     ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$  
     !-----------

     !     call init_us_1  ! CG

     !<CG>     
     call init_gipaw_1
     !</CG>     

     ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
     !  Definition of a particular indexation to avoid Mickael Profeta's crazy indexation
     ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$  
     ! <CG>
     call init_xanes_paw_nhm
     !<\CG>
     allocate (paw_iltonhb(0:paw_lmaxkb,xanes_paw_nhm, ntyp))
     call define_index_arrays(paw_iltonhb)

     ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
     !  Allocate paw projectors
     ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$  

     allocate (paw_vkb( npwx,  paw_nkb))
     allocate (paw_becp(paw_nkb, nbnd))

  elseif(xonly_plot) then  !$$$$$$$$$$$$$$$$$$  xonly_plot if structure
     ! Fermi level read from file
     ef=ef_r
     write( stdout,*)
     write( stdout,'(">> Fermi level read from file (Ryd)=",1f14.8)') ef
     write( stdout,*)
  endif   !$$$$$$$$$$$$$$$$$$  xonly_plot if structure

  !<CG>
  ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  ! $ Computing gamma tabulated values
  ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  if (trim(gamma_mode).eq.'file') then
     if (ionode) call read_gamma_file
     call mp_bcast( gamma_lines, ionode_id )
     if (.not.ionode) allocate (gamma_points(gamma_lines,2))
     call mp_bcast( gamma_points, ionode_id )
  elseif (trim(gamma_mode).eq.'variable') then
     gamma_lines=2
     allocate(gamma_points(2,2))
     gamma_points(1,1)=gamma_energy(1)
     gamma_points(2,1)=gamma_energy(2)
     gamma_points(1,2)=gamma_value(1)
     gamma_points(2,2)=gamma_value(2)
  endif

  if ((trim(gamma_mode).eq.'file').or.(trim(gamma_mode).eq.'variable')) then
     allocate( gamma_tab(xnepoint))
     call initialize_gamma_tab
     deallocate(gamma_points)
  endif

  ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  !</CG>

  xnitermax=xniter

  if(xonly_plot) then
     call read_header_save_file(x_save_file)
     nks = nkstot
     write(6,*) 'nks=',nks
     if(lsda) then
        isk(1:nkstot/2)=1
        isk(nkstot/2+1:nkstot)=2
        wk(1:nkstot)=2.d0/nkstot
     elseif(.not.lsda) then
        isk(1:nkstot)=1
        wk(1:nkstot)=2.d0/nkstot
     endif
     CALL divide_et_impera( xk, wk, isk, lsda, nkstot, nks )
  endif

  ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  !  Verification of paw relations between pseudo partial waves and projector (radial parts)
  ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$  

  if(.not.xonly_plot.and.trim(verbosity).eq.'high')  call check_paw_projectors(xiabs)


  ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  !  Allocate and initialise lanczosvectors
  ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$  

  allocate(a(xnitermax,1,nks))
  allocate(b(xnitermax,1,nks))
  allocate(xnorm(1,nks))
  allocate(ncalcv(1,nks))

  a(:,:,:)=0.d0
  b(:,:,:)=0.d0
  xnorm(:,:)=0.d0
  ncalcv(:,:)=0

  ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  !  And now we go...  XANES CALCULATION
  ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$  

  if (okvan) then
    write(stdout,*) 'Ultrasoft not implemented'
    call stop_pp
  endif


  if(trim(calculation).eq.'xanes') then
     if(.not.xonly_plot) then

        write(stdout,*) ' '
        write(stdout,*) 'Approx. ram memory needed per proc in MB = ', (16*3*npwx*npool)/(nproc*1000000)
        write(stdout,*) 

        if(xang_mom.eq.1) then
           call xanes_dipole(a,b,ncalcv,xnorm,core_wfn,paw_iltonhb,terminator,verbosity)
           ! open save file and write everything
           call write_save_file(a,b,xnorm,ncalcv,x_save_file)
        elseif(xang_mom.eq.2) then
           call xanes_quadrupole(a,b,ncalcv,xnorm,core_wfn,paw_iltonhb,terminator,verbosity)
           call write_save_file(a,b,xnorm,ncalcv,x_save_file)
        endif
     else
        call read_save_file(a,b,xnorm,ncalcv,x_save_file,core_energy)
     endif

     if(xang_mom.eq.1) then
        call plot_xanes_dipole(a,b,xnorm,ncalcv,terminator,core_energy)
     elseif(xang_mom.eq.2) then
        call plot_xanes_quadrupole(a,b,xnorm,ncalcv,terminator,core_energy)
     endif

  elseif(trim(calculation).eq.'rxes') then
     call errore( 'Main', 'rxes Not yet implemented',1)
  elseif(trim(calculation).eq.'bethe_salpeter') then
     call errore( 'Main', 'bethe_salpeter Not yet implemented',1)
  elseif(trim(calculation).eq.'hpsi') then
     call verify_hpsi
  endif

  ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  !  Deallocation
  ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$  

!  if(.not.xonly_plot) then
!     if( gamma_only ) then
!        DEALLOCATE(rbecp)
!     else
!        DEALLOCATE(becp)
!     endif
!  endif

  deallocate(a)
  deallocate(b)
  deallocate(xnorm)
  deallocate(ncalcv)


  write (stdout,*) 'End program ', trim(calculation)


  call stop_pp

end program Xxanes









! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!  Dipolar Calculation
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$  


subroutine xanes_dipole(a,b,ncalcv,xnorm,core_wfn,paw_iltonhb,terminator,verbosity)
  USE io_files,         ONLY : nd_nmbr, prefix, tmp_dir, &
       nwordwfc, iunwfc
  USE io_global,        ONLY : stdout     ! Modules/io_global.f90
  USE kinds,            only : DP
  use parameters,       ONLY : ntypx
  USE radial_grids,     ONLY : ndmx
  USE ions_base,        ONLY : nat, ntyp => nsp, ityp
  use wvfct,            ONLY : npwx, nbndx, nbnd, npw, igk, g2kin, et
  use lsda_mod,         ONLY : nspin,lsda,isk,current_spin
  use cell_base,        only: tpiba2, bg
  use wavefunctions_module, only: evc
  use klist,            ONLY : &
       nkstot,            & ! total number of k-points
       nks,               & ! number of k-points per pool
       xk,                & ! k-points coordinates
       wk                   ! k-points weight
  use gvect,            ONLY: g,ngm,ecutwfc,ngl,nrxx
  !,ig_l2g(ngm),ngm_l,ngm_g
  use paw_gipaw,        ONLY : &
       paw_vkb,             & ! |p> projectors
       paw_becp,            & ! product of projectors and wf.
       paw_nkb,             & ! total number of beta functions, with st.fact.
       paw_lmaxkb,paw_recon
  use becmod,     ONLY : becp,rbecp,allocate_bec, deallocate_bec !CG
  use scf,        ONLY : vltot, vrs, v, kedtau !CG
  use gsmooth,    ONLY : doublegrid
  USE mp_global,  ONLY : intra_pool_comm
  USE mp,         ONLY : mp_sum

  USE xanes,      only : xiabs, xanes_dip, xang_mom, xniter, xnitermax, xepsilon

  USE atom,       ONLY : rgrid, msh
  !  use atom,        ONLY : &
  !       mesh,     &!mesh(ntypx) number of mesh points
  !       msh ,     &!msh(ntypx)the point at rcut=end of radial integration
  !       r   
  use radin_mod

  USE uspp,   ONLY : vkb, nkb, okvan !CG

  USE ldaU,   ONLY : lda_plus_u
  !<CG>
  use xanes_paw_variables, only : xanes_paw_nhm
  !</CG>

  implicit none
  real(dp) core_wfn(ndmx)
  real(dp) a(xnitermax,1,nks),b(xnitermax,1,nks)     
  real (dp)  xnorm(1,nks)
  integer :: is,ik,iabso,nr,ip,jp,l,j,icrd,ip_l,nrc,nt,na
  integer :: ipx,ipx_0,ipy,ipz,nline,nrest,npw_partial
  integer :: ncalcv(1,nks)
  integer :: paw_iltonhb(0:paw_lmaxkb,xanes_paw_nhm, ntyp)
  real (dp) pref,prefb,pi,arg,v_of_0,xnorm_partial
  real (dp) norm
  real (dp), allocatable :: aux(:)
  complex(kind=DP) :: ZDOTC
  complex(dp), allocatable :: paw_vkb_cplx(:,:)
  logical :: terminator
  real(dp) :: normps
  character(len=4) :: verbosity

  external ZDOTC
  external ZDSCAL

  complex(dp), allocatable :: psiwfc(:)

  ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  !  Constant Definitions
  ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$  

  pi=2.d0*dasin(1.d0)
  pref=sqrt(3.d0/2.d0)
  prefb=sqrt(3.d0)
  arg=sqrt(4.d0*pi)


  ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  !  Variable allocation and initialization
  ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$  

  allocate(aux(rgrid(xiabs)%mesh)) !allocation too big, it needs only up to msh
  allocate (paw_vkb_cplx( npwx,  paw_nkb))
  allocate(xanes_dip(paw_recon(xiabs)%paw_nl(xang_mom)))
  allocate(psiwfc(npwx))
  xanes_dip(:)=0.d0


  ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  ! Radial Dipole Matrix element
  ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$  

  !  I compute the radial part, 
  !          <\phi_n|r|\psi_i>=int r^3 \phi_n \psi_i dr.
  !  Here only ps and ae wavefunctions from XX.recon
  !

  write( stdout,*) 'Calculation dipole matrix element'
  write( stdout,*) 'There are ',paw_recon(xiabs)%paw_nl(xang_mom),' projectors/channel'
  write( stdout,*) 'for angular moment ',xang_mom,' and atom type ',xiabs

  ! I check that the core wf is correctly normalized

  nr=msh(xiabs)  ! extended up to all the NON ZERO points in the mesh.

  if(trim(verbosity).eq.'high') then
     aux(1:nr)=core_wfn(1:nr)*core_wfn(1:nr)
     write (stdout,'(" norm of core wfc =",1f14.8)') sqrt(para_radin(aux(1:nr),rgrid(xiabs)%r(1:nr),nr))
  endif


  ! and I calculate the radial integral

  ip_l=0

  do ip=1,paw_recon(xiabs)%paw_nbeta
     !     if(psphi(xiabs,ip)%label%l.eq.xang_mom) then
     if(paw_recon(xiabs)%aephi(ip)%label%l.eq.xang_mom) then

        nrc=paw_recon(xiabs)%aephi(ip)%label%nrc
        ip_l=ip_l+1
        !    here below, recall that psi is r*psi and you have a Jacobian=r^2
        aux(1:nrc)=rgrid(xiabs)%r(1:nrc)*paw_recon(xiabs)%aephi(ip)%psi(1:nrc)*core_wfn(1:nrc)
        !    here we have to integrate only inside the augmentation region.
        xanes_dip(ip_l)=para_radin(aux(1:nrc),rgrid(xiabs)%r(1:nrc),nrc)

     endif
  enddo


  write(6,*) '----------------------------------------------------------------'
  do ip=1,paw_recon(xiabs)%paw_nl(xang_mom) !
     write( stdout,'("dipole radial matrix element proj. (",i2,")=",f14.8)') ip,xanes_dip(ip)
  enddo
  write(6,*) '----------------------------------------------------------------'
  deallocate(aux)


  !
  !  We count the projectors for all the atoms
  !

  ipx_0=0
  if(xiabs.ne.1) then
     do nt=1, xiabs-1
        do na=1, nat
           if (ityp(na).eq.nt) then
              ipx_0=ipx_0+paw_recon(nt)%paw_nh
           endif
        enddo
     enddo
  endif


  ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  ! Beginning the loop over the k-points
  ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$  

  !*apsi  call set_vrs(vrs,vltot,vr,nrxx,nspin,doublegrid)
  !<CG>
  call set_vrs(vrs,vltot,v%of_r,kedtau, v%kin_r,nrxx,nspin,doublegrid)
  !</CG>

  !  call newd   ! CG

    if (lda_plus_u) call init_xanes_ldau

  do ik=1,nks

     if(lsda) current_spin=isk(ik)

     !gk_sort  sort k-points and exit kinetic energies 
     call gk_sort(xk (1,ik),ngm,g,ecutwfc/tpiba2,npw,igk,g2kin)  !CHECK
     g2kin=g2kin*tpiba2                                          !CHECK

     npw_partial = npw
     call mp_sum( npw_partial, intra_pool_comm )

     if(xniter.ge.npw_partial) then
        xniter = npw_partial
        write(stdout,*) 'Hilbert space is saturated'
        write(stdout,*) 'xniter is set equal to ',npw_partial
        write(stdout,*) 'Hint: Increase Kinetic Energy cutoff in your SCF simulation'
     endif


     !<CG>        
     call init_gipaw_2(npw,igk,xk(1,ik),paw_vkb)
     !</CG>
     call init_us_2(npw,igk,xk(1,ik),vkb)

     ! initialise (not orthogonalized) atomic wfc if lda_plus_u=T 
          if (lda_plus_u) call init_xanes_ldau_2(ik)

     ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
     ! Angular Matrix element
     ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$  
     !
     ! Here I define human projectors
     !
     ! paw_vkb(1:npw,ipx) are real spherical harmonics. The real spherical harmonics are
     ! defined as
     ! 
     !     y_{l,2m}  =[Y_{l,m}+(-1)^m Y_{l,-m}]/sqrt(2)
     !     y_{l,2m+1}=[Y_{l,m}-(-1)^m Y_{l,-m}]/(i*sqrt(2))
     !
     ! (remember Y_{l,m}=(-1)^m Y_{l,-m)^*  )
     !
     ! The complex spherical harmonics can be written has a function of the real
     ! ones as :
     !
     !     Y_{l,m}  =        [y_{l,2m}+iy_{l,2m+1}]/sqrt(2)     
     !     Y_{l,-m} = (-1)^m [y_{l,2m}-iy_{l,2m+1}]/sqrt(2)     
     !
     !  The paw_vkb_cplx are the Y_{l,m} so the usual spherical harmonics
     !
     ! rotational invariance has been checked


     do ip=1,paw_recon(xiabs)%paw_nl(xang_mom)

        ipx=ipx_0+paw_iltonhb(xang_mom,ip,xiabs)


        paw_vkb_cplx(1:npw,ipx)= paw_vkb(1:npw,ipx)   !m=0

        paw_vkb_cplx(1:npw,ipx+1)=       &
             (paw_vkb(1:npw,ipx+1)+(0.d0,1.d0)*paw_vkb(1:npw,ipx+2))/sqrt(2.0)    !m=+1

        paw_vkb_cplx(1:npw,ipx+2)=-       &
             (paw_vkb(1:npw,ipx+1)-(0.d0,1.d0)*paw_vkb(1:npw,ipx+2))/sqrt(2.0)    !m=-1



     enddo

     psiwfc(1:npw)=(0.d0,0.d0)
     do ip=1,paw_recon(xiabs)%paw_nl(xang_mom)   
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
             pref*paw_vkb_cplx(1:npw,ipx+2)*cmplx(xepsilon(1), xepsilon(2))&
             + pref*paw_vkb_cplx(1:npw,ipx+1)*cmplx(-xepsilon(1), xepsilon(2))&
             +prefb*paw_vkb_cplx(1:npw,ipx)*xepsilon(3) &
             )*xanes_dip(ip)/arg


     enddo
     psiwfc(1:npw)=psiwfc(1:npw)*sqrt(4*pi)/3.0




     ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
     ! Starting Lanczos
     ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$  


     !   
     !   I normalize the wavefunction psiwfc(1:npw)
     !

     call allocate_bec(nkb,1) ! CG

     xnorm_partial=ZDOTC(npw,psiwfc,1,psiwfc,1)

     call mp_sum( xnorm_partial, intra_pool_comm )

     xnorm(1,ik)=sqrt(xnorm_partial)
     write( stdout,*) 'norm initial vector=',xnorm(1,ik)
     norm=1.d0/xnorm(1,ik)

     call ZDSCAL(npw,norm,psiwfc,1)
     !
     !      Then I call the lanczos routine
     !
     write(stdout,*) 'Starting lanczos'

     call lanczos(a(:,1,ik),b(:,1,ik),psiwfc,ncalcv(1,ik), terminator)

     

     !
     !      Then I write small report of the lanczos results
     !

     if(trim(verbosity).eq.'high') then
        write( stdout,*) '-----------------------------------------'
        write( stdout,*) 'k-point number =',ik
        write( stdout,*) 'k-point coordinate, isk'
        write( stdout,'(3f12.6,1i2)') xk(1,ik),xk(2,ik),xk(3,ik),isk(ik)
        write( stdout,*) 'Norm of the initial vector =',xnorm(1,ik)
        write( stdout,*) 'Number of iterations =',ncalcv(1,ik)

        !        nline=ncalcv(icrd,ik)/6
        !        nrest=ncalcv(icrd,ik)-nline*6
        !        write( stdout,*) 'a vectors:'
        !        do ip=1,nline
        !           write( stdout,"(6(f10.6,3x))") (a((ip-1)*6+j,icrd,ik),j=1,6)
        !        enddo
        !        write( stdout,"(6(f10.6,3x))") (a(nline*6+j,icrd,ik),j=1,nrest)
        !        write( stdout,*) 'b vectors:'
        !        do ip=1,nline
        !           write( stdout,"(6(f10.6,3x))") (b((ip-1)*6+j,icrd,ik),j=1,6)
        !        enddo
        !        write( stdout,"(6(f10.6,3x))") (b(nline*6+j,icrd,ik),j=1,nrest)
        write( stdout,*) '-----------------------------------------'
     endif

  call deallocate_bec ( ) ! CG

  enddo  !on k points

  ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  ! Array deallocation
  ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$  

  deallocate(psiwfc)
  deallocate(xanes_dip)
  deallocate (paw_vkb_cplx)

end subroutine xanes_dipole

! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!  Quadrupolar Calculation
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$  


subroutine xanes_quadrupole(a,b,ncalcv,xnorm,core_wfn,paw_iltonhb,terminator,verbosity)
  USE io_files,        ONLY : nd_nmbr, prefix, tmp_dir, &
       nwordwfc, iunwfc
  USE io_global,       ONLY : stdout     ! Modules/io_global.f90
  USE kinds, only : DP
  use parameters,       only : ntypx
  USE radial_grids,     ONLY : ndmx
  USE ions_base,   ONLY : nat, ntyp => nsp, ityp
  use wvfct,            ONLY : npwx,nbndx,nbnd,npw,igk,&
       g2kin,et
  !       ,igk_l2g
  use lsda_mod,    ONLY : nspin,lsda,isk,current_spin
  use cell_base, only: tpiba2, bg
  use wavefunctions_module, only: evc
  use klist,       ONLY : &
       nkstot,            & ! total number of k-points
       nks,               & ! number of k-points per pool
       xk,                & ! k-points coordinates
       wk                   ! k-points weight
  use gvect, only: g,ngm,ecutwfc,ngl,nrxx
  !,ig_l2g(ngm),ngm_l,ngm_g
  use paw_gipaw,     only : &
       paw_vkb,             & ! |p> projectors
       paw_becp,            & ! product of projectors and wf.
       paw_nkb,             & ! total number of beta functions, with st.fact.
       paw_lmaxkb, & 
       paw_recon
  use becmod, ONLY:becp,rbecp,allocate_bec, deallocate_bec ! CG
  use scf, ONLY: vltot,v,vrs, kedtau !CG
  use gsmooth, ONLY : doublegrid
  USE mp_global,  ONLY : intra_pool_comm
  USE mp,         ONLY : mp_sum
  USE xanes,  only:  xiabs,xanes_qua,xang_mom,xniter,xnitermax,xkvec,xepsilon
  USE atom,       ONLY : rgrid, msh
!  use atom,        ONLY : &
!       mesh,     &!mesh(ntypx) number of mesh points
!       msh ,     &!msh(ntypx)the point at rcut=end of radial integration
!       r   
  use radin_mod
  USE uspp,   ONLY : vkb, nkb, okvan !CG
  USE ldaU,   ONLY : lda_plus_u
!<CG>
  use xanes_paw_variables, only : xanes_paw_nhm
!</CG>

  implicit none
  real(dp) core_wfn(ndmx)
  real(dp) a(xnitermax,1,nks),b(xnitermax,1,nks)     
  real (dp)  xnorm(1,nks)
  integer is,ik,iabso,nr,ip,jp,l,j, ip_l,nrc,nt,na,ipx_0
  integer ipx,ipy,ipz,nline,nrest,npw_partial,icrd
  integer ncalcv(1,nks)
  integer paw_iltonhb(0:paw_lmaxkb,xanes_paw_nhm, ntyp) !CG
  
  real (dp) pref,prefb,pi,arg,v_of_0,xnorm_partial
  real (dp) norm,prefm2,prefm1,prefm0
  real (dp), allocatable :: aux(:)
  complex(kind=dp), allocatable :: psi(:)
  complex(kind=DP) :: ZDOTC
  complex(dp), allocatable :: paw_vkb_cplx(:,:)
  logical terminator
  real(dp) :: normps

  external ZDOTC
  external ZDSCAL

  complex(dp), allocatable :: psiwfc(:)
  character(len=4) :: verbosity



! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!  Constant Definitions
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$  


  pi=2.d0*dasin(1.d0)
  prefm2=sqrt(3.0/40.0)/3.0
  prefm1=prefm2
  prefm0=2.0*sqrt(1.0/40.0)/3.0

  pref=sqrt(2.d0)
  prefb=1.0/pref


! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!  Variable allocation and initialization
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$  

  allocate(aux(rgrid(xiabs)%mesh))!overdimensionated, necessary only up to msh
  allocate(psi(npwx))
  allocate (paw_vkb_cplx( npwx,  paw_nkb))
  allocate(xanes_qua(paw_recon(xiabs)%paw_nl(xang_mom)))
  allocate(psiwfc(npwx))

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

  write( stdout,*) 'Calculation Quadrupole matrix element'
  write( stdout,*) 'There are ',paw_recon(xiabs)%paw_nl(xang_mom),' projectors/channel'
  write( stdout,*) 'xang_mom=',xang_mom,' xiabs=',xiabs
  
  ! I check that the fondamental orthogonality condition of paw is satisfied:

  !     nr=mesh(xiabs)  ! extended up to all the points in the mesh.
  nr=msh(xiabs)  ! extended up to all the NON ZERO points in the mesh.

        
! I check that the core wf is correctly normalized

  if(trim(verbosity).eq.'high') then
     aux(1:nr)=core_wfn(1:nr)*core_wfn(1:nr)
     write (stdout,'("norm of core wfc =",1f14.8)') sqrt(para_radin(aux(1:nr),rgrid(xiabs)%r(1:nr),nr))
  endif

! and I calculate the radial integral

  ip_l=0

  do ip=1,paw_recon(xiabs)%paw_nbeta
     if(paw_recon(xiabs)%psphi(ip)%label%l.eq.xang_mom) then
        nrc=paw_recon(xiabs)%aephi(ip)%label%nrc
        ip_l=ip_l+1
        !    here below, recall that psi is r*psi and you have a Jacobian=r^2
        !        aux(1:nr)=r(1:nr,xiabs)*r(1:nr,xiabs)* &
        !             psphi(xiabs,ip)%psi(1:nr)*core_wfn(1:nr)
        aux(1:nrc)=rgrid(xiabs)%r(1:nrc)*rgrid(xiabs)%r(1:nrc)*&
             paw_recon(xiabs)%aephi(ip)%psi(1:nrc)*core_wfn(1:nrc)
        !    here we have to integrate only inside the augmentation region.
        xanes_qua(ip_l)=para_radin(aux(1:nrc),rgrid(xiabs)%r(1:nrc),nrc)
     endif
  enddo


  do ip=1,paw_recon(xiabs)%paw_nl(xang_mom) !
     write( stdout,*)  'XANES quadrupole matrix element proj',ip,': ',xanes_qua(ip)
  enddo
  
  deallocate(aux)
  deallocate(psi)


  !
  !  We count the projectors for all the atoms
  !
  

  ipx_0=0
  if(xiabs.ne.1) then
     do nt=1, xiabs-1
        do na=1, nat
           if (ityp(na).eq.nt) then
              ipx_0=ipx_0+paw_recon(nt)%paw_nh
           endif
        enddo
     enddo
  endif

 
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! Beginning the loop over the k-points
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$  


!*apsi  call set_vrs(vrs,vltot,vr,nrxx,nspin,doublegrid)
!  call set_vrs(vrs,vltot,v%of_r,nrxx,nspin,doublegrid)
!  call newd

!<CG>
! set_vrs (vrs, vltot, vr, kedtau, kedtaur,nrxx, nspin, doublegrid)
  call set_vrs(vrs,vltot,v%of_r,kedtau, v%kin_r,nrxx,nspin,doublegrid)
!</CG>

  if (lda_plus_u) call init_xanes_ldau

  do ik=1,nks

     if(lsda) current_spin=isk(ik)

     !gk_sort  sort k-points and exit kinetic energies 
     call gk_sort(xk (1,ik),ngm,g,ecutwfc/tpiba2,npw,igk,g2kin)  !CHECK
     g2kin=g2kin*tpiba2                                          !CHECK

     npw_partial = npw
     call mp_sum( npw_partial, intra_pool_comm )
     if(xniter.ge.npw_partial) then
        xniter = npw_partial
        write(stdout,*) 'Hilbert space is saturated'
        write(stdout,*) 'xniter is set equal to ',npw_partial
        write(stdout,*) 'Hint: Increase Kinetic Energy cutoff in your SCF simulation'
        !        call stop_pp
     endif

!<CG>
     call init_gipaw_2(npw,igk,xk(1,ik),paw_vkb)
!</CG>
     call init_us_2(npw,igk,xk(1,ik),vkb)

  ! initialise orthogonalized atomic wfc if lda_plus_u=T
     if (lda_plus_u) call init_xanes_ldau_2(ik)

     !
     ! Here I define human projectors
     !
     ! paw_vkb(1:npw,ipx) are real spherical harmonics. The real spherical harmonics are
     ! defined as
     ! 
     !     y_{l,2m}  =[(-)^{m} Y_{l,m}+Y_{l,-m}]/sqrt(2)
     !     y_{l,2m+1}=[(-)^{m} Y_{l,m}-Y_{l,-m}]/(i*sqrt(2))
     !
     ! The complex spherical harmonics can be written has a function of the real
     ! ones as :
     !
     !     Y_{l,m}  =       [y_{l,2m}+iy_{l,2m+1}]/sqrt(2)     
     !     Y_{l,-m} = (-)^m [y_{l,2m}-iy_{l,2m+1}]/sqrt(2)     
     !
     !  The paw_vkb_cplx are the Y_{l,m}, so the usual spherical harmonics
     !



     do ip=1,paw_recon(xiabs)%paw_nl(xang_mom)   
        
        ipx=ipx_0+paw_iltonhb(xang_mom,ip,xiabs)
        paw_vkb_cplx(1:npw,ipx)= paw_vkb(1:npw,ipx)   !m=0
        paw_vkb_cplx(1:npw,ipx+1)=       &
             (paw_vkb(1:npw,ipx+1)+(0.d0,1.d0)*paw_vkb(1:npw,ipx+2))/sqrt(2.0)   !m=+1
        paw_vkb_cplx(1:npw,ipx+2)=       &
             -(paw_vkb(1:npw,ipx+1)-(0.d0,1.d0)*paw_vkb(1:npw,ipx+2))/sqrt(2.0)    !m=-1
        paw_vkb_cplx(1:npw,ipx+3)=       &
             (paw_vkb(1:npw,ipx+3)+(0.d0,1.d0)*paw_vkb(1:npw,ipx+4))/sqrt(2.0)   !m=+2
        paw_vkb_cplx(1:npw,ipx+4)=       &
             (paw_vkb(1:npw,ipx+3)-(0.d0,1.d0)*paw_vkb(1:npw,ipx+4))/sqrt(2.0)    !m=-2


     enddo

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

     do ip=1,paw_recon(xiabs)%paw_nl(xang_mom)   
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



     enddo

! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! Starting Lanczos
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$  

     call allocate_bec(nkb,1) ! CG

!   
!   I normalize the wavefunction pwswfc(1:npw)
!
     xnorm_partial=ZDOTC(npw,psiwfc,1,psiwfc,1)

     call mp_sum( xnorm_partial, intra_pool_comm )

     xnorm(1,ik)=sqrt(xnorm_partial)
     write( stdout,*) 'norm initial vector=',xnorm(1,ik)
     norm=1.d0/xnorm(1,ik)

     call ZDSCAL(npw,norm,psiwfc,1)
     !
     !      Then I call the lanczos routine
     !
     call lanczos(a(:,1,ik),b(:,1,ik),psiwfc,ncalcv(1,ik), terminator)
     !
     !      Then I write small report of the lanczos results
     !
     
     if(trim(verbosity).eq.'high') then
        write( stdout,*) '-----------------------------------------'
        write( stdout,*) 'k-point number =',ik
        write( stdout,*) 'Norm of the initial vector =',xnorm(1,ik)
        write( stdout,*) 'Number of iterations =',ncalcv(1,ik)
        write( stdout,*) 'k-point coordinate, isk'
        write( stdout,'(3f12.6,1i2)') xk(1,ik),xk(2,ik),xk(3,ik),isk(ik)
!     nline=ncalcv(1,ik)/6
!     nrest=ncalcv(1,ik)-nline*6
!     write( stdout,*) 'a vectors:'
!     do ip=1,nline
!        write( stdout,"(6(f10.6,3x))") (a((ip-1)*6+j,1,ik),j=1,6)
!     enddo
!     write( stdout,"(6(f10.6,3x))") (a(nline*6+j,1,ik),j=1,nrest)
!     write( stdout,*) 'b vectors:'
!     do ip=1,nline
!        write( stdout,"(6(f10.6,3x))") (b((ip-1)*6+j,1,ik),j=1,6)
!     enddo
!     write( stdout,"(6(f10.6,3x))") (b(nline*6+j,1,ik),j=1,nrest)
        write( stdout,*) '-----------------------------------------'
     endif

     call deallocate_bec() ! CG

  enddo   !en loop over k points

  ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  ! Array deallocation
  ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$  

  deallocate(psiwfc)
  deallocate (paw_vkb_cplx)
  deallocate(xanes_qua)



end subroutine xanes_quadrupole



subroutine lanczos (a,b,psi,ncalcv, terminator)
  USE kinds, only : DP
  use wvfct,  only:  npwx,nbndx, nbnd,npw,igk,g2kin
  !use wavefunctions_module, ONLY : psic
  use becmod, ONLY:becp,rbecp
  USE gsmooth,  ONLY : nls, nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, nrxxs
  USE uspp,   ONLY : vkb, nkb
  USE cell_base, only:omega
  USE xanes, ONLY : xniter, xnepoint, xcheck_conv,xnitermax,xemin,xemax,xgamma,xerror
  USE mp_global,  ONLY : intra_pool_comm  
  USE mp,         ONLY : mp_sum
  USE io_global,       ONLY : stdout

  implicit none
  ! INPUT
  integer ncalcv
  real(dp) :: a(xnitermax),b(xnitermax)
  complex(kind=DP) :: psi(npwx)
  logical terminator

  !INTERNAL
  logical converge
  logical iconv
  integer ibnd,j,i,m
  real(kind=dp)  norm,error,xemin_ryd,xemax_ryd,ryd2eV,xgamma_ryd
  parameter(ryd2ev = 13.6058d0)
  complex(kind=dp):: ac,bc
  real(kind=dp), allocatable :: comp(:)
  complex(kind=dp), allocatable :: hpsi(:),u(:)
  real (kind=dp) :: ddot
  complex(kind=DP) :: ZDOTC
  external ZDOTC,ddot
  External h_psi

  allocate(hpsi(npwx))
  allocate(u(npwx))
  allocate(comp(xnepoint))

  hpsi(:)=(0.d0,0.d0)
  u(:)=(0.d0,0.d0)
  a(:)=0.d0
  b(:)=0.d0

  xemax_ryd=xemax/Ryd2eV
  xemin_ryd=xemin/Ryd2eV
  xgamma_ryd=xgamma/Ryd2eV

  iconv=.false.

  ! ------------------------  1st iteration --------------------------
  !
  ! -- computes H*Psi  (|u0>=|Psi>)

  CALL h_psi( npwx, npw,1, psi, hpsi )


  ! -- computes a_(1)=<Psi|HPsi>

  a(1)=dble(ZDOTC(npw,psi,1,hpsi,1))

  call mp_sum(a(1), intra_pool_comm)

  ac=-a(1)*(1.d0,0.d0)
  !
  ! -- computes t3=H*Psi-a_1*Psi

  call zaxpy(npw,ac,psi,1,hpsi,1)

  !
  ! -- computes the norm of t3

  b(1) = ZDOTC(npw,hpsi,1,hpsi,1)
  call mp_sum( b(1), intra_pool_comm )
  b(1) = sqrt( b(1) )
  !
  ! -- computes the vector |u1>

  call zdscal(npw,1.d0/b(1),hpsi,1)

  !
  ! -- saves |u1>


  u(1:npw)=hpsi(1:npw)

  hpsi(:)=(0.d0,0.d0)


  !
  ! -------------------------- Next iterations -----------------------
  !


  comp(:)=0.d0
  comp(1)=1.d0



  do i=2,xniter


     !From here below we have:
     !   u=psi_j
     !   hpsi=empty
     !   psi=psi_j-1

     CALL h_psi( npwx, npw, 1, u, hpsi )


     !
     ! I compute hpsi=hpsi-b_{j-1}*psi_{j-1} (j is actually i-1)

     bc=-b(i-1)*(1.d0,0.d0)
     call zaxpy(npw,bc,psi,1,hpsi,1)

     !
     ! computes a(i)=<t2|t3>=<t2|H|t2>


     a(i)=real(ZDOTC(npw,hpsi,1,u,1),dp)
     call mp_sum( a(i), intra_pool_comm )
     !
     ! computes t3=t3-a(i)*t2
     !

     ac=-a(i)*(1.d0,0.d0)
     call zaxpy(npw,ac,u,1,hpsi,1)


     !
     ! computes b(i) the norm of t3
     !

     b(i)=real(ZDOTC(npw,hpsi,1,hpsi,1),dp)
     call mp_sum( b(i), intra_pool_comm )
     b(i) = sqrt( b(i) )


     !
     ! saves initial vector in t1 (t1 = t2)

     psi(1:npw)=u(1:npw)
     !
     ! computes vector t3/norm(t3) and saves it in t2

     call zdscal(npw,1.d0/b(i),hpsi,1)
     u(1:npw)=hpsi(1:npw)
     hpsi(1:npw)=(0.d0,0.d0)


     !
     !   I should gather all the ncalcv,error at the end and write only
     !   then.
     !


     if(mod(i,xcheck_conv).eq.0) then
        if(converge(a,b,i,comp,error,xemin_ryd,xemax_ryd,xgamma_ryd,xnepoint,xerror,terminator)) then
           write( stdout,*) 'CONVERGED at iter ',i,' with error=',error
           ncalcv=i
           iconv=.true.
           exit
        else
           write( stdout,*) 'Estimated error at iter ',i,' is ', error
        endif
     endif


  enddo

  if(.not.iconv) then
     write( stdout,*) 'XANES not converged after', i-1, ' iterations'
     write( stdout,*) 'Estimated final error after ',i-1,'iterations is ', &
          converge(a,b,i-1,comp,error,xemin_ryd,xemax_ryd,xgamma_ryd,xnepoint,xerror,terminator)
     ncalcv=i-1
  endif



  deallocate(hpsi)
  deallocate(u)
  deallocate(comp)


end subroutine lanczos


logical function converge(a,b,m,comp,estimated_error,xemin,xemax,xgamma,xnepoint,xerror,use_term)
  USE kinds, ONLY : dp
  USE xanes, ONLY : xnitermax
  implicit none

  !
  ! INPUT:
  ! -----
  !
  integer :: xnepoint
  real(dp) :: xemin,xemax,xgamma,xerror

  !
  ! INPUT / OUTPUT:
  ! ----------------
  !
  real(dp) :: a(xnitermax),b(xnitermax)
  integer :: m        ! number of calculated vectors
  real(dp) :: comp(xnepoint)
  real(dp) :: estimated_error  
  !
  ! ---------------------- local variables ------------------------------
  !
  real(dp) :: deltae,tmp,area,e,err
  real(dp)  :: continued_fraction
  integer  :: n
  logical :: use_term
  !
  deltae = (xemax-xemin) / xnepoint
  err = 0.d0
  area = 0.d0
  n = 0
  e = xemin
  do n = 1, xnepoint
    e = e + deltae
    tmp = continued_fraction(a,b,e,xgamma,m,use_term)
    err = err + abs(comp(n)-tmp)
    area = area + abs(tmp)
    comp(n) = tmp
  enddo
  err = err / area
  estimated_error = err
  if(err < xerror) then
     converge = .true.
  else
     converge = .false.
  endif
end function converge
!
! ------------------------------------------------------------------------
!
function continued_fraction(a,b,e,gamma,m, term)
  USE kinds, ONLY : dp
  USE xanes, ONLY : xnitermax, xcheck_conv
   USE io_global,       ONLY : stdout
  implicit none

  ! computes the continued fraction.
  !
  ! INPUT:
  ! -----
  !
  real(dp) :: continued_fraction
  integer  :: m
  real(dp) :: a(xnitermax)
  real(dp) :: b(xnitermax)
  real(dp) :: gamma
  real(dp) :: e
  logical :: term
  !
  ! ---------------------- local variables ------------------------------
  !
  integer :: i, p,q
  complex(dp) :: res ,lastterm ! intermediate variable
  real(dp) :: aa, bb

  q=xcheck_conv/2
  if (term) then
    aa=0.0
    bb=0.0
    do p=1, q
	aa=aa+a(m-p)
	bb=bb+b(m-p)
    enddo
    aa=aa/q
    bb=bb/q
    res=lastterm(aa-e,bb*bb,gamma)
  else
  res = cmplx(a(m)-e,gamma)
  endif
  do i = 1, m -1
    res = cmplx(a(m-i)-e, -gamma) -b(m-i)*b(m-i)/res
  enddo
  continued_fraction = aimag(1/res)
end function continued_fraction
!

subroutine read_core_abs_paratec(filename,core_wfn,xiabs)
  USE io_global,       ONLY : stdout
  USE kinds,           only : DP
  USE ions_base,       ONLY : ityp,nat
  USE radial_grids,    only : ndmx
  USE atom,            ONLY : rgrid, msh
  !  use atom,        ONLY : mesh, msh, r     !mesh(ntypx) number of mesh points
  use splinelib

  implicit none
  INTEGER, INTENT ( IN ) :: xiabs

  logical :: core_found
  integer :: i,icounter, iostatus,na, jtyp, j
  character (LEN=80) :: filename
  real(kind=dp):: x(ndmx),d1
  real(kind=dp):: core_wfn_para(ndmx),core_wfn(*)
  character (LEN=90)::dummy_str
  real(dp), allocatable :: xdata(:), tab(:), tab_d2y(:)

  core_found=.FALSE.

  open(unit=33,file=filename,form='formatted',status='old')

  rewind(33)
  do
     read(33,'(a)',iostat = iostatus) dummy_str
     if ( iostatus /= 0 ) exit

     if ( index ( dummy_str, "core wavefunctions" ) /= 0 ) then
        core_found=.TRUE.
        read(33,'(a)') dummy_str
        if ( dummy_str(1:7) /= "#n=1l=0" ) then
           call errore ( "read_core_abs_paratec", "missing 1s state", 0 )
        end if
        icounter=0
        do 
           read(33,'(a)') dummy_str
           if ( dummy_str(1:1) == "#" .or. dummy_str(1:1) == "&" ) then
              exit
           else
              icounter=icounter+1
              if(icounter.gt.ndmx) &
                   call errore ( "read_core_abs_paratec", "ndmx too small", 0 )
              read(dummy_str,*) x(icounter),core_wfn_para(icounter)
           end if
        end do
     end if
  end do
  if(.not.core_found) &
       call errore ( "read_core_abs_paratec", "no core state(s) found", 0 )

  close(33)

  jtyp=0
  do na=1,nat
     if(ityp(na).eq.xiabs) jtyp=jtyp+1
  enddo

  ! Interpolate to the grid of the pseudo potential
  allocate( xdata(icounter), tab(icounter), tab_d2y(icounter) )

  xdata = x ( 1:icounter )
  tab = core_wfn_para ( 1:icounter )

  ! initialize spline interpolation; for 3.x, x >= 1
  d1 = (tab(2) - tab(1)) / (xdata(2) - xdata(1))
  call spline(xdata, tab, 0.d0, d1, tab_d2y)
  !  call spline(xdata, tab, tab_d2y)

  do j = 1, msh(jtyp)
     core_wfn(j) = splint(xdata, tab, tab_d2y, rgrid(jtyp)%r(j))
  end do

  deallocate( xdata, tab, tab_d2y )
end subroutine read_core_abs_paratec
   
subroutine read_core_abs(filename,core_wfn)
  USE kinds, only : DP
  USE atom,        ONLY : rgrid
  !use atom,        ONLY : mesh     !mesh(ntypx) number of mesh points
  use xanes,       ONLY : xiabs
  implicit none
  integer :: i
  character (LEN=80) :: filename
  real(kind=dp):: x
  real(kind=dp):: core_wfn(*)
  
  open(unit=33,file=filename,form='formatted',status='old')
  rewind(33)
  read(33,*)
  do i=1, rgrid(xiabs)%mesh
     !   read(33,'(2f16.8)') x,core_wfn(i)
     read(33 ,*) x,core_wfn(i)
  enddo
  close(33)
end subroutine read_core_abs


! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!  Continuum fraction and plotting dipolar part
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$  


subroutine plot_xanes_dipole(a,b,xnorm,ncalcv, terminator, core_energy)
  USE kinds, only : DP
  USE xanes, only : xang_mom,xemax,xemin,xiabs,xnepoint,xgamma,xonly_plot
  !*apsi  USE uspp_param, ONLY : psd  !psd(ntypx) label for the atoms 
  USE klist, ONLY : nelec, &
       nkstot,            & ! total number of k-points
       nks,               & ! number of k-points per pool
       xk,                & ! k-points coordinates
       wk                   ! k-points weight
  USE ener, ONLY : ef
  USE io_global,       ONLY : stdout,ionode  
  USE mp_global,  ONLY : intra_pool_comm, inter_pool_comm !CG
  USE xanes, ONLY : xnitermax,xkvec,xepsilon
  use cut_valence_green, only : cut_occ_states, cut_desmooth, memu, meml, cut_nmemu, cut_nmeml
  use lsda_mod,    ONLY : nspin,isk
  !<CG>
  USE mp,                   ONLY : mp_sum
  USE uspp_param,           ONLY : upf
  use gamma_variable_mod  , only : gamma_tab, gamma_mode
  !</CG>
  implicit none
  logical terminator
  integer iestart
  integer ncalcv(1,nks)
  real(dp) a(xnitermax,1,nks),b(xnitermax,1,nks)     
  real (dp)  xnorm(1,nks)
  real(kind=dp) ryd2eV,alpha2,pi,core_energy
  parameter(ryd2ev = 13.6058d0)
  integer :: i,ik,n,icoord           !loops
  integer :: lmax
  real (dp) :: mygetK,e_1s,energy,de,mod_xgamma,xemax_ryd,xemin_ryd,xgamma_ryd
  real (dp) :: tmp_var
  real (dp) :: Intensity_coord(1,xnepoint,nspin)
  !   , Intensity_tot(xnepoint)
  real(dp)  :: continued_fraction, paste_fermi, desmooth,t1,t2,f1,f2,df1,df2,poly(4) !CG 
  logical :: first

  ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  !  constant and initialization
  ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$  

  pi=2.d0*dasin(1.d0)

  if(xonly_plot) then
     e_1s=core_energy           !the K-edge in eV.
  else
     !<CG>
     e_1s=mygetK(upf(xiabs)%psd) !mygetK gets the K-edge in eV.
     !</CG>
  endif

  e_1s=e_1s/Ryd2eV            !This is in Rydberg

  alpha2=4.d0*pi/137.04

  desmooth=cut_desmooth/ryd2ev  ! This is in rydberg

  !
  ! small output
  !

  write( stdout,"(/'CALCULATION XANES SPECTRA')")
  write( stdout,"('---------------------')")
  write( stdout,"(/' Final state angular momentum:',1x,i3)") xang_mom
  if (trim(gamma_mode).eq.'constant')  write( stdout,"(' Broadening parameter (in eV):',1x,g20.12)") xgamma     !check
  write( stdout,"(' Nb of energy points:',1x,i5)") xnepoint
  write( stdout,"(' Maximum energy (in eV):',1x,g20.12)") xemax            !check
  write( stdout,"(' Minimum energy (in eV):',1x,g20.12)") xemin            !check
  write( stdout,"(' Binding energy of the 1s level (in eV):',1x,g20.12)") -e_1s*Ryd2eV

  if( ionode ) then
     open (unit=277,file='xanes.dat',form='formatted',status='unknown')
     rewind(277)
     write(277,"('# final state angular momentum:',1x,i3)") xang_mom
     if (trim(gamma_mode).eq.'constant')     write(277,"('# broadening parameter (in eV):',1x,g20.12)") xgamma      !

     write(277,"('# absorbing atom type:',i4)") xiabs

     if(nspin.gt.1) then
        write(277,"('# Energy (eV)   sigmatot   sigmaup    sigmadown ')")
     else
        write(277,"('# Energy (eV)   sigma')")
     endif

  endif

  !
  ! I convert in Rydberg most of the relevant quantities
  !

  xemax_ryd=xemax/Ryd2eV+ef
  xemin_ryd=xemin/Ryd2eV+ef
  xgamma_ryd=xgamma/Ryd2eV

  write(stdout,*) 'xemin(ryd)=',xemin_ryd
  write(stdout,*) 'xemax(ryd)=',xemax_ryd
  if (trim(gamma_mode).eq.'constant') then
     write(stdout,*) 'gamma in rydberg=',xgamma_ryd
  else
     write(stdout,*) 'nonconstant gamma'
  endif


  ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  !  Continuum fraction
  ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$  

  Intensity_coord(:,:,:)=0.d0


  de = (xemax_ryd-xemin_ryd) / real(xnepoint-1)

  if (trim(gamma_mode).eq.'constant') then

     if(cut_occ_states) then
        allocate(memu(cut_nmemu,2))
        allocate(meml(cut_nmeml,2))
        iestart=(ef-xemin_ryd)/de
        do ik=1,nks
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
           call determine_polycut(t1,t2,f1,f2,df1,df2,poly) ! calculates interpolation polynome

           do n=1,xnepoint
              energy=xemin_ryd+de*(n-1)
              if ((energy-ef<desmooth).and.(energy-ef>-desmooth)) then  ! interpolation 
                 tmp_var=poly(1)+poly(2)*energy+poly(3)*energy**2+poly(4)*energy**3
                 tmp_var=tmp_var*xnorm(1,ik)*xnorm(1,ik)
                 !</CG>
              else
                 tmp_var=0.d0
                 if (n>iestart) then
                    tmp_var=  &
                         continued_fraction(a(1,1,ik),b(1,1,ik),energy,xgamma_ryd,ncalcv(1,ik)-1,terminator)*  &
                         xnorm(1,ik)*xnorm(1,ik)
                 endif
                 tmp_var = tmp_var + paste_fermi(energy,ef,a(1,1,ik),b(1,1,ik),xgamma_ryd,ncalcv(1,ik)-1,terminator, first) &
                      *xnorm(1,ik)*xnorm(1,ik)
              endif
              !            Intensity_tot(n)=Intensity_tot(n)+tmp_var*wk(ik)
              Intensity_coord(1,n,isk(ik)) = Intensity_coord(1,n,isk(ik))+tmp_var*wk(ik)
           enddo
        enddo
        deallocate(memu)
        deallocate(meml)


     else
        do ik=1,nks
           do n=1,xnepoint
              energy=xemin_ryd+de*(n-1)
              tmp_var=  &
                   continued_fraction(a(1,1,ik),b(1,1,ik),energy,xgamma_ryd,ncalcv(1,ik)-1, terminator)*  &
                   xnorm(1,ik)*xnorm(1,ik)
              Intensity_coord(1,n,isk(ik)) = Intensity_coord(1,n,isk(ik))+tmp_var*wk(ik)
           enddo
        enddo
     endif

  else ! nonconstant gamma

     if(cut_occ_states) then
        allocate(memu(cut_nmemu,2))
        allocate(meml(cut_nmeml,2))
        iestart=(ef-xemin_ryd)/de
        do ik=1,nks
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
           call determine_polycut(t1,t2,f1,f2,df1,df2,poly)

           do n=1,xnepoint
              energy=xemin_ryd+de*(n-1)
              xgamma_ryd=gamma_tab(n)
              if ((energy-ef<desmooth).and.(energy-ef>-desmooth)) then  ! interpolation
                 tmp_var=poly(1)+poly(2)*energy+poly(3)*energy**2+poly(4)*energy**3
                 tmp_var=tmp_var*xnorm(1,ik)*xnorm(1,ik)
              else
                 tmp_var=0.d0
                 if (n>iestart) then
                    tmp_var=  &
                         continued_fraction(a(1,1,ik),b(1,1,ik),energy,xgamma_ryd,ncalcv(1,ik)-1,terminator)*  &
                         xnorm(1,ik)*xnorm(1,ik)
                 endif
                 tmp_var = tmp_var + paste_fermi(energy,ef,a(1,1,ik),b(1,1,ik),xgamma_ryd,ncalcv(1,ik)-1,terminator, first) &
                      *xnorm(1,ik)*xnorm(1,ik)
              endif
              !            Intensity_tot(n)=Intensity_tot(n)+tmp_var*wk(ik)
              Intensity_coord(1,n,isk(ik)) = Intensity_coord(1,n,isk(ik))+tmp_var*wk(ik)
           enddo
        enddo
        deallocate(memu)
        deallocate(meml)


     else
        do ik=1,nks
           do n=1,xnepoint
              energy=xemin_ryd+de*(n-1)
              xgamma_ryd=gamma_tab(n)
              tmp_var=  &
                   continued_fraction(a(1,1,ik),b(1,1,ik),energy,xgamma_ryd,ncalcv(1,ik)-1, terminator)*  &
                   xnorm(1,ik)*xnorm(1,ik)
              Intensity_coord(1,n,isk(ik)) = Intensity_coord(1,n,isk(ik))+tmp_var*wk(ik)
           enddo
        enddo
     endif

  endif ! gamma_mode

  !  call poolreduce( nspin*xnepoint, Intensity_coord )

  !<CG>  replaces poolreduce
#ifdef __PARA
  call mp_sum ( Intensity_coord, inter_pool_comm )
#endif
  !</CG>

  ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  !  Plotting xanes spectrum
  ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$  

  if(ionode) then
     if(nspin.eq.1) then
        do n=1,xnepoint
           energy=xemin_ryd+de*(n-1)
           Intensity_coord(:,n,:)=Intensity_coord(:,n,:)*(energy+e_1s)*alpha2 !
           write(277,'(2f14.8)') (energy-ef)*ryd2eV, Intensity_coord(1,n,1)
        enddo
     elseif(nspin.eq.2) then
        do n=1,xnepoint
           energy=xemin_ryd+de*(n-1)
           Intensity_coord(:,n,:)=Intensity_coord(:,n,:)*(energy+e_1s)*alpha2 !
           write(277,'(4f14.8)') (energy-ef)*ryd2eV, &
                Intensity_coord(1,n,1)+Intensity_coord(1,n,2),&
                Intensity_coord(1,n,1),Intensity_coord(1,n,2)
        enddo
     endif

     close(277)
  endif

end subroutine plot_xanes_dipole



! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!  Continuum fraction and plotting quadrupolar part
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$  



subroutine plot_xanes_quadrupole(a,b,xnorm,ncalcv, terminator,core_energy)
  USE kinds, only : DP
  USE xanes, only : xang_mom,xemax,xemin,xiabs,xnepoint,xgamma,xonly_plot
  !*apsi  USE uspp_param, ONLY : psd  !psd(ntypx) label for the atoms 
  USE klist, ONLY : nelec, &
       nkstot,            & ! total number of k-points
       nks,               & ! number of k-points per pool
       xk,                & ! k-points coordinates
       wk                   ! k-points weight
  USE ener, ONLY : ef
  USE io_global,       ONLY : stdout,ionode  
  USE mp_global,  ONLY : intra_pool_comm, inter_pool_comm
  USE xanes, ONLY : xnitermax,xepsilon
  use cut_valence_green , only : cut_occ_states, cut_desmooth, cut_nmemu, cut_nmeml, memu, meml
  use lsda_mod,    ONLY : nspin,isk
  !<CG>
  USE mp,                   ONLY : mp_sum
  USE uspp_param,           ONLY : upf
  use gamma_variable_mod, only : gamma_tab, gamma_mode
  !</CG>
  implicit none
  logical terminator
  integer iestart
  integer ncalcv(1,nks)
  real(dp) a(xnitermax,1,nks),b(xnitermax,1,nks)     
  real (dp)  xnorm(1,nks)
  real(kind=dp) ryd2eV,alpha2,pi,constantqua
  parameter(ryd2ev = 13.6058d0)
  integer :: i,ik,n,icoord           !loops
  integer :: lmax
  real (dp) :: mygetK,e_1s,energy,de,mod_xgamma,xemax_ryd,xemin_ryd,xgamma_ryd
  real (dp) :: tmp_var,core_energy
  real (dp) :: Intensity_tot(xnepoint,nspin)
  real(dp)  :: continued_fraction, paste_fermi, desmooth,t1,t2,f1,f2,df1,df2,poly(4)
  logical :: first

  ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  !  constant and initialization
  ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$  

  pi=2.d0*dasin(1.d0)
  constantqua=pi/(137.04*137.04*137.04)
  if(xonly_plot) then
     e_1s=core_energy           !the K-edge in eV.
  else
     !<CG>
     e_1s=mygetK(upf(xiabs)%psd) !mygetK gets the K-edge in eV.
     !</CG>
  endif
  e_1s=e_1s/Ryd2eV            !This is in Rydberg

  desmooth=cut_desmooth/ryd2ev  ! This is in rydberg

  alpha2=4.d0*pi/137.04

  !
  ! small output
  !

  write( stdout,"(/'CALCULATION XANES SPECTRA')")
  write( stdout,"('---------------------')")
  write( stdout,"(/' Final state angular momentum:',1x,i3)") xang_mom
  if (trim(gamma_mode).eq.'constant')   write( stdout,"(' Broadening parameter (in eV):',1x,g20.12)") xgamma     !check
  write( stdout,"(' Nb of energy points:',1x,i5)") xnepoint
  write( stdout,"(' Maximum energy (in eV):',1x,g20.12)") xemax            !check
  write( stdout,"(' Minimum energy (in eV):',1x,g20.12)") xemin            !check
  write( stdout,"(' Binding energy of the 1s level (in eV):',1x,g20.12)") -e_1s*Ryd2eV


  if( ionode ) then
     open (unit=277,file='xanes.dat',form='formatted',status='unknown')
     rewind(277)
     write(277,"('# final state angular momentum:',1x,i3)") xang_mom
     if (trim(gamma_mode).eq.'constant')  write(277,"('# broadening parameter (in eV):',1x,g20.12)") xgamma      

     write(277,"('# absorbing atom type:',i4)") xiabs
     if(nspin.eq.1) then
        write(277,"('# Energy (eV)   sigmatot')")
     elseif(nspin.eq.2) then
        write(277,"('# Energy (eV)   sigmatot, sigmaup, sigmadown')")
     endif
  endif

  !
  ! I convert in Rydberg most of the relevant quantities
  !

  xemax_ryd=xemax/Ryd2eV+ef
  xemin_ryd=xemin/Ryd2eV+ef
  xgamma_ryd=xgamma/Ryd2eV

  write(stdout,*) 'xemin(ryd)=',xemin_ryd
  write(stdout,*) 'xemax(ryd)=',xemax_ryd
  if (trim(gamma_mode).eq.'constant') then
     write(stdout,*) 'gamma in rydberg=',xgamma_ryd
  else
     write(stdout,*) 'gamma from file'
  endif


  !  write(stdout,*) 'terminator=', terminator


  ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  !  Continuum fraction
  ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$  


  Intensity_tot(:,:)=0.d0

  de = (xemax_ryd-xemin_ryd) / real(xnepoint-1)

  if (trim(gamma_mode).eq.'constant') then
     if(cut_occ_states) then
        allocate(memu(cut_nmemu,2))
        allocate(meml(cut_nmeml,2))
        do ik=1,nks
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
           call determine_polycut(t1,t2,f1,f2,df1,df2,poly)

           do n=1,xnepoint
              energy=xemin_ryd+de*(n-1)
              if ((energy-ef<desmooth).and.(energy-ef>-desmooth)) then  !interpolation
                 tmp_var=poly(1)+poly(2)*energy+poly(3)*energy**2+poly(4)*energy**3
                 tmp_var=tmp_var*xnorm(1,ik)*xnorm(1,ik)
                 !</CG>
              else
                 tmp_var=0.d0
                 if (n>iestart) then
                    tmp_var=  &
                         continued_fraction(a(1,1,ik),b(1,1,ik),energy,xgamma_ryd,ncalcv(1,ik)-1,terminator)*  &
                         xnorm(1,ik)*xnorm(1,ik)
                 endif
                 tmp_var=tmp_var+paste_fermi(energy,ef,a(1,1,ik),b(1,1,ik),xgamma_ryd,ncalcv(1,ik)-1,terminator, first)&
                      *xnorm(1,ik)*xnorm(1,ik)
              endif
              Intensity_tot(n,isk(ik))=Intensity_tot(n,isk(ik))+tmp_var*wk(ik)
           end do
        enddo
        deallocate(memu)
        deallocate(meml)

     else
        !     write(stdout,*) 'in plor_
        do ik=1,nks
           do n=1,xnepoint
              energy=xemin_ryd+de*(n-1)
              tmp_var=  &
                   continued_fraction(a(1,1,ik),b(1,1,ik),energy,xgamma_ryd,ncalcv(1,ik)-1, terminator)*  &
                   xnorm(1,ik)*xnorm(1,ik)
              Intensity_tot(n,isk(ik))=Intensity_tot(n,isk(ik))+tmp_var*wk(ik)
           enddo
        enddo
     endif

  else ! nonconstant gamma

     if(cut_occ_states) then
        allocate(memu(cut_nmemu,2))
        allocate(meml(cut_nmeml,2))
        do ik=1,nks
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
           call determine_polycut(t1,t2,f1,f2,df1,df2,poly)

           do n=1,xnepoint
              energy=xemin_ryd+de*(n-1)
              xgamma_ryd=gamma_tab(n)
              if ((energy-ef<desmooth).and.(energy-ef>-desmooth)) then  !interpolation
                 tmp_var=poly(1)+poly(2)*energy+poly(3)*energy**2+poly(4)*energy**3
                 tmp_var=tmp_var*xnorm(1,ik)*xnorm(1,ik)
                 !</CG>
              else
                 tmp_var=0.d0
                 if (n>iestart) then
                    tmp_var=  &
                         continued_fraction(a(1,1,ik),b(1,1,ik),energy,xgamma_ryd,ncalcv(1,ik)-1,terminator)*  &
                         xnorm(1,ik)*xnorm(1,ik)
                 endif
                 tmp_var=tmp_var+paste_fermi(energy,ef,a(1,1,ik),b(1,1,ik),xgamma_ryd,ncalcv(1,ik)-1,terminator, first)&
                      *xnorm(1,ik)*xnorm(1,ik)
              endif
              Intensity_tot(n,isk(ik))=Intensity_tot(n,isk(ik))+tmp_var*wk(ik)
           end do
        enddo
        deallocate(memu)
        deallocate(meml)

     else
        !     write(stdout,*) 'in plor_
        do ik=1,nks
           do n=1,xnepoint
              energy=xemin_ryd+de*(n-1)
              xgamma_ryd=gamma_tab(n)
              tmp_var=  &
                   continued_fraction(a(1,1,ik),b(1,1,ik),energy,xgamma_ryd,ncalcv(1,ik)-1, terminator)*  &
                   xnorm(1,ik)*xnorm(1,ik)
              Intensity_tot(n,isk(ik))=Intensity_tot(n,isk(ik))+tmp_var*wk(ik)
           enddo
        enddo
     endif

  endif ! gammma_mode

  !  call poolreduce( nspin*xnepoint, Intensity_tot )

  !<CG>  replaces poolreduce
#ifdef __PARA
  call mp_sum ( Intensity_tot, inter_pool_comm )
#endif
  !</CG>


  ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  !  Plotting xanes spectrum
  ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$  

  if(ionode) then
     if(nspin.eq.1) then
        do n=1,xnepoint
           energy=xemin_ryd+de*(n-1)
           Intensity_tot(n,:)=Intensity_tot(n,:)*(energy+e_1s)*(energy+e_1s)*(energy+e_1s)*constantqua     !normalized
           write(277,'(2f14.8)') (energy-ef)*ryd2eV,Intensity_tot(n,:)
        enddo
     elseif(nspin.eq.2) then
        do n=1,xnepoint
           energy=xemin_ryd+de*(n-1)
           Intensity_tot(n,:)=Intensity_tot(n,:)*(energy+e_1s)*(energy+e_1s)*(energy+e_1s)*constantqua     !normalized
           write(277,'(4f14.8)') (energy-ef)*ryd2eV,Intensity_tot(n,1)+Intensity_tot(n,2),&
                Intensity_tot(n,1),Intensity_tot(n,2)
        enddo
     endif
     close(277)
  endif



end subroutine plot_xanes_quadrupole


subroutine define_index_arrays(paw_iltonhb)
  USE paw_gipaw,   ONLY : paw_lmaxkb, paw_recon
  USE ions_base,   ONLY : ntyp => nsp
  USE parameters,  ONLY : lmaxx
  USE xanes_paw_variables, only : xanes_paw_nhm ! CG
  ! Arguments
  integer paw_iltonhb(0:paw_lmaxkb,xanes_paw_nhm,ntyp) ! CG

  ! Local
  integer nt,ih,nb,l
  integer ip_per_l(0:lmaxx) 

  do nt = 1, ntyp
     ih = 1
     ip_per_l(:)=0
     do nb = 1, paw_recon(nt)%paw_nbeta
        l = paw_recon(nt)%aephi(nb)%label%l
        ip_per_l(l)=ip_per_l(l)+1
        paw_iltonhb(l,ip_per_l(l),nt)= ih
        ih=ih+2*l+1
     enddo
  enddo


end subroutine define_index_arrays



function lastterm(a,b,g)
  USE kinds, ONLY : dp
  implicit none
  real(dp)  :: a,b, g,y1,y2,z1,z2,r
  complex(dp) :: lastterm

y1 =a*a-g*g-4*b
y2=-2*a*g
r=sqrt(y1*y1+y2*y2)/2

if (g<0) then
z1=a/2+0.5*sign(sqrt(y1/2+r),y2)
z2=-g/2+sqrt(-y1/2+r)/2
else
z1=a/2-0.5*sign(sqrt(y1/2+r),y2)
z2=-g/2-sqrt(-y1/2+r)/2
endif
 
 lastterm=cmplx(z1,z2)

end function lastterm



function paste_fermi(e,ef,a,b,gamma,m,term, first)
  USE kinds, ONLY : dp
  USE xanes, ONLY : xnitermax, xcheck_conv
  use cut_valence_green, only :&
         cut_ierror,cut_stepu, cut_stepl,&
         cut_startt,cut_tinf, cut_tsup, cut_nmemu, cut_nmeml, memu, meml
  implicit none

  real(dp) :: paste_fermi
  integer  :: m
  real(dp) :: a(xnitermax)
  real(dp) :: b(xnitermax)
  real(dp) :: gamma
  real(dp) :: e,ef
  logical :: term
  complex(dp) :: green,y,dy,c1,c2, e1, e2
  real(dp) ::twopi, t,dt, t1, ta, tb
  integer, save :: n1
  integer, save :: n2
  integer :: nn1,nn2
  logical :: first

  if (first) then
   memu(:,:)=cmplx(0.d0,0.d0)
   meml(:,:)=cmplx(0.d0,0.d0) 
   n1=0
   n2=0
   first=.false.
  endif

  dy=cut_ierror+1.0
  y=0.d0
  twopi=6.28318530

  nn1=1
  nn2=1

  t1=0.5773502692
  
  t=cut_startt

  do while ((abs(dy)>cut_ierror).or.(t<cut_tsup))
   dt=cut_stepu*t
   ta=t+dt*(1-t1)/2
   tb=t+dt*(1+t1)/2
   e1=cmplx(ef,ta)
   e2=cmplx(ef,tb)

   if (nn1>n1) then
    c1=green(a,b,e1,m,term)
    c2=green(a,b,e2,m,term)
    if (nn1<cut_nmemu) then
     memu(nn1,1)=c1
     memu(nn1,2)=c2
     n1=nn1
    endif 
   else
    c1=memu(nn1,1)
    c2=memu(nn1,2)
   endif

   dy=(dt/2)*(c1/cmplx(ef-e,ta-gamma)+conjg(c1)/cmplx(ef-e,-ta-gamma)+c2/cmplx(ef-e,tb-gamma)+conjg(c2)/cmplx(ef-e,-tb-gamma))
   y=y+dy
   t=t+dt
   nn1=nn1+1
  end do

  t=cut_startt
  dy=cut_ierror+1

  do while((abs(dy)>cut_ierror).or.(t>cut_tinf))
   dt=cut_stepl*t
   ta=t-dt*(1-t1)/2
   tb=t-dt*(1+t1)/2
   e1=cmplx(ef,ta)
   e2=cmplx(ef,tb)

   if (nn2>n2) then
    c1=green(a,b,e1,m,term)
    c2=green(a,b,e2,m,term)
    if (nn2<cut_nmeml) then
     meml(nn2,1)=c1
     meml(nn2,2)=c2
     n2=nn2
    endif
   else
    c1=meml(nn2,1)
    c2=meml(nn2,2)
   endif

   dy=(dt/2)*(c1/cmplx(ef-e,ta-gamma)+conjg(c1)/cmplx(ef-e,-ta-gamma)+c2/cmplx(ef-e,tb-gamma)+conjg(c2)/cmplx(ef-e,-tb-gamma))
   y=y+dy
   t=t-dt
   nn2=nn2+1
  end do

  paste_fermi=aimag(y)/twopi

end function paste_fermi


function green(a,b,e,m, term)
  USE kinds, ONLY : dp
  USE xanes, ONLY : xnitermax, xcheck_conv
  implicit none

  complex(dp) :: green
  integer  :: m
  real(dp) :: a(xnitermax)
  real(dp) :: b(xnitermax)
  complex(dp) :: e
  logical :: term
  integer :: i, p,q
  complex(dp) :: res ,lastterm 
  real(dp) :: aa, bb


  q=xcheck_conv/2
  if (term) then
    aa=0.0
    bb=0.0
    do p=1, q
	aa=aa+a(m-p)
	bb=bb+b(m-p)
    enddo
    aa=aa/q
    bb=bb/q

   res=lastterm(aa-real(e),bb*bb,aimag(e))
  else
  res = cmplx(a(m)-real(e),aimag(e))
  endif
  do i = 1, m -1
    res = a(m-i)-e -b(m-i)*b(m-i)/res
  enddo
  
  green = 1/res
  
end function green


subroutine check_paw_projectors(xiabs)
  USE kinds,           only : DP
  USE paw_gipaw,       only : &
        paw_lmaxkb, &
        paw_recon
  use xanes_paw_variables, only : xanes_paw_nhm
  USE atom,            ONLY : rgrid, msh
!  USE atom,  ONLY : &
!       mesh,     &!mesh(ntypx) number of mesh points              
!       msh ,     &!msh(ntypx)the point at rcut=end of radial integration 
!       r, rab
  USE ions_base,       ONLY : ntyp => nsp
  USE io_global,       ONLY : stdout
  use radin_mod
  implicit none
  integer xiabs
  ! internal
  integer :: nr,nrc,ip,jp,lmax,l,ip_l,jtyp,n1,n2,nrs,ndm,ih,jh
  real(dp) :: pi,overlap,rexx,overlap2
  real (dp), allocatable :: aux(:),f(:,:)
  real(dp) , allocatable :: s(:,:),e(:),v(:,:)

  pi=2.d0*dasin(1.d0)

  allocate(aux(rgrid(xiabs)%mesh)) !allocation too big, it needs only up to msh
  allocate(f(rgrid(xiabs)%mesh,2)) !allocation too big, it needs only up to msh

  
  write(stdout,*) '----  PAW projectors from reconstruction files -----'
  write(stdout,*)
  write(stdout,*) 'atom type,  total   number of projectors'
  do jtyp=1,ntyp
     write (stdout,'(2i4)') jtyp,paw_recon(jtyp)%paw_nbeta
  enddo
  write(stdout,*)

  !
  ! I calculate maximum l
  !
  
  lmax=0
  do ip=1,paw_recon(xiabs)%paw_nbeta
     if(paw_recon(xiabs)%psphi(ip)%label%l.gt.lmax) &
          lmax = paw_recon(xiabs)%psphi(ip)%label%l
  enddo


  write(stdout,*) 'atom type,  l,   number of projectors per ang. mom.'

  do jtyp=1,ntyp
     do l=0,lmax
        write(stdout,'(3i4)') jtyp,l,paw_recon(jtyp)%paw_nl(l)
     enddo
  enddo



  ! We calculate the overlaps between partial waves and projectors
  ! to see if they are equal to the croneker delta.

  nr=msh(xiabs)  ! extended up to all the NON ZERO points in the mesh.


  write(stdout,*) '----  Overlaps between partial waves and projectors (radial) -----'
  write(stdout,*)
  write(stdout,*) '<tilde{phi} l,n|tilde{p} l,nn>=delta_{n,nn}'
  write(stdout,*) 

  
  do ip=1,paw_recon(xiabs)%paw_nbeta
     do jp=1,paw_recon(xiabs)%paw_nbeta
        if(paw_recon(xiabs)%psphi(ip)%label%l.eq.paw_recon(xiabs)%psphi(jp)%label%l) then
           nrc=Count(rgrid(xiabs)%r(1:nr).le.paw_recon(xiabs)%psphi(ip)%label%rc)
           if(nrc.gt.nr) then
              write(stdout,*) 'nrc=',nrc,' > ',nr,' = nr' 
              call errore ( "nrc > nr", "xanes_dipole", 0 )
           endif
           aux(1:nrc)=paw_recon(xiabs)%psphi(ip)%psi(1:nrc)*paw_recon(xiabs)%paw_betar(1:nrc,jp)
           aux(nrc+1:nr)=0.d0
           write(stdout,'("<tilde{phi}_",2i2,10X,"|tilde{p}_",2i2,">=",1f14.8)')  &
                ip,paw_recon(xiabs)%psphi(ip)%label%l,jp, &
                paw_recon(xiabs)%psphi(jp)%label%l, &
                para_radin(aux(1:nr),rgrid(xiabs)%r(1:nr),nr)
        endif
     enddo
  enddo

  write(stdout,*)

  !
  !
  !

  write(stdout,*) '---- Check normalization pseudo,ae wf and projectors -----------'
  write(stdout,*) '----    (radial part only, integral up to r_c)    -----------'
  write(stdout,*)
  write(stdout,*) 'l,   n, |proj|^2, |pswf|^2 , |aewf|^2'
  write(stdout,*)
  do l=0,lmax
     do ip=1,paw_recon(xiabs)%paw_nbeta
        if(paw_recon(xiabs)%psphi(ip)%label%l.eq.l) then
           nrc=Count(rgrid(xiabs)%r(1:nr).le.paw_recon(xiabs)%psphi(ip)%label%rc)
           aux(1:nrc) = paw_recon(xiabs)%paw_betar(1:nrc,ip) &
                      * paw_recon(xiabs)%paw_betar(1:nrc,ip)
           overlap=para_radin(aux(1:nrc),rgrid(xiabs)%r(1:nrc),nrc)
           aux(1:nrc)=paw_recon(xiabs)%aephi(ip)%psi(1:nrc)*paw_recon(xiabs)%aephi(ip)%psi(1:nrc)
           overlap2=para_radin(aux(1:nrc),rgrid(xiabs)%r(1:nrc),nrc)
           aux(1:nrc)=paw_recon(xiabs)%psphi(ip)%psi(1:nrc)*paw_recon(xiabs)%psphi(ip)%psi(1:nrc)
           rexx=para_radin(aux(1:nrc),rgrid(xiabs)%r(1:nrc),nrc)
           write(stdout,'(2i4,3f14.8)')l,ip,overlap,overlap2,rexx
        endif
     enddo
  enddo
  write(stdout,*)
  
  goto 323
  write(stdout,*) '---- <phi|chi>= \sum_nl <phi |phi_l><p_l| chi>_nrc  --------'
  write(stdout,*) 'WARNING : this test assumes a form of the phi/chi function'
  
  !  do l=0,lmax
  do l=1,1
     ip_l=0
     do ip=1,paw_recon(xiabs)%paw_nbeta
        if(ip_l.eq.0.and.paw_recon(xiabs)%psphi(ip)%label%l.eq.l) ip_l=ip
     enddo
     
     f(:,:)=0.d0
     do ip=1,paw_recon(xiabs)%paw_nbeta
        if(paw_recon(xiabs)%psphi(ip)%label%l.eq.l) then
           f(1:nr,1)=f(1:nr,1)+paw_recon(xiabs)%psphi(ip)%psi(1:nr)/real(ip,dp)
           f(1:nr,2)=f(1:nr,2)+1.123*paw_recon(xiabs)%psphi(ip)%psi(1:nr)/real(ip,dp)
        endif
     enddo
     rexx=0.d0
     do ip=1,paw_recon(xiabs)%paw_nbeta
        if(paw_recon(xiabs)%psphi(ip)%label%l.eq.l) then
           nrc=Count(rgrid(xiabs)%r(1:nr).le.paw_recon(xiabs)%psphi(ip)%label%rc)
           if(nrc.gt.nr) then
              write(stdout,*) 'nrc=',nrc,' > ',nr,' = nr' 
              call errore ( "nrc > nr", "xanes_dipole", 0 )
           endif
           aux(1:nrc)=f(1:nrc,1)*paw_recon(xiabs)%paw_betar(1:nrc,ip)
           overlap=para_radin(aux(1:nrc),rgrid(xiabs)%r(1:nrc),nrc)
           aux(1:nrc)=f(1:nrc,2)*paw_recon(xiabs)%psphi(ip)%psi(1:nrc)
           overlap=overlap*para_radin(aux(1:nrc),rgrid(xiabs)%r(1:nrc),nrc)
           write(stdout,'("overlap(l=",1i2,",n=",1i2,")= ",1f14.8)')l,ip,overlap
           rexx=rexx+overlap
        endif
     enddo
     aux(1:nr)=f(1:nr,1)*f(1:nr,2)
     write(stdout,'("sum/overlap=",1f14.8)') rexx/para_radin(aux,rgrid(xiabs)%r(1:nr),nrc)
     write(stdout,'("sum projectors=",1f14.8," overlap=",1f14.8)') rexx,para_radin(aux,rgrid(xiabs)%r(1:nr),nrc)
     write(stdout,*)'---+++++----'
  enddo
  !           aux(1:nrc)=paw_recon(xiabs)%psphi(ip)%psi(1:nrc)*paw_recon(xiabs)%paw_betar(1:nrc,jp)

  !           
  !
  !
  !        endif
  !          enddo
  write(stdout,*)
  write(stdout,*) '================================================================'
  
323 continue
  !
  !  Check linear dependence of projectors
  !

  write(stdout,*) '================================================================'
  write(stdout,*) '           Checking linear dipendence of projectors             '
  write(stdout,*)

  deallocate(aux)

  ndm = MAXVAL (msh(1:ntyp))

  allocate(aux(ndm))

  do l=0,paw_lmaxkb
     if (paw_recon(xiabs)%paw_nl(l)>0) then
        allocate (s(paw_recon(xiabs)%paw_nl(l),paw_recon(xiabs)%paw_nl(l)))
        allocate (e(paw_recon(xiabs)%paw_nl(l)),v(paw_recon(xiabs)%paw_nl(l),paw_recon(xiabs)%paw_nl(l)))
        do ih=1,paw_recon(xiabs)%paw_nl(l)
           n1=paw_recon(xiabs)%paw_iltonh(l,ih)
           nrc=paw_recon(xiabs)%psphi(n1)%label%nrc
           nrs=paw_recon(xiabs)%psphi(n1)%label%nrs
           do jh=1,paw_recon(xiabs)%paw_nl(l)
              n2=paw_recon(xiabs)%paw_iltonh(l,jh)
              call step_f(aux,paw_recon(xiabs)%psphi(n1)%psi(1:msh(xiabs)) * &
                   paw_recon(xiabs)%psphi(n2)%psi(1:msh(xiabs)), &
                   rgrid(xiabs)%r(1:msh(xiabs)),nrs,nrc, 1.d0, msh(xiabs) )
              call simpson ( msh(xiabs), aux, rgrid(xiabs)%rab(1), s(ih,jh))
           enddo
        enddo
  
        write(stdout,'("atom type=",1i4)') xiabs
        write(stdout,'("number of projectors projector  =",1i4," angular momentum=",1i4)') &
             paw_recon(xiabs)%paw_nl(l),l
        do ih=1,paw_recon(xiabs)%paw_nl(l)
           write(stdout,'(10f14.8)') (s(ih,jh),jh=1,paw_recon(xiabs)%paw_nl(l)) 
        enddo
        write(stdout,*) 'Eigenvalues S matrix:'
        write(stdout,*)

        if(paw_recon(xiabs)%paw_nl(l).eq.1) then
           write(stdout,'(1i4,1f14.8)') 1,s(1,1)
        else 
           call rdiagh(paw_recon(xiabs)%paw_nl(l),s,paw_recon(xiabs)%paw_nl(l) , e, v )
           do ih=1,paw_recon(xiabs)%paw_nl(l)
              write(stdout,'(1i4,1f14.8)') ih,e(ih)
           enddo
        endif
        write(stdout,*)
        deallocate(s,e,v)
     endif
  enddo
  write(stdout,*) '================================================================'


  deallocate(aux,f)
end subroutine check_paw_projectors


subroutine read_save_file(a,b,xnorm,ncalcv,x_save_file,core_energy)
  USE kinds,       only : DP
  USE klist,       ONLY : nks,nkstot
  USE xanes,       ONLY : xnitermax,xang_mom,xiabs
  USE io_global,   ONLY : stdout,ionode
  USE lsda_mod,    ONLY : nspin,lsda
  implicit none
  character(LEN=256) :: x_save_file
  integer            :: ierr,nkstot_r
  integer            :: xm_r,nc_r,ncomp_max
  integer ncalcv(1,nks)
  real(dp) core_energy
  real(dp) a(xnitermax,1,nks),b(xnitermax,1,nks)     
  real (dp)  xnorm(1,nks)
  integer, allocatable :: ncalcv_all(:,:)
  real(dp), allocatable :: a_all(:,:),b_all(:,:),xnorm_all(:,:),aux(:,:)
  real(dp) xkvec_r(3),xepsilon_r(3)
  integer i,j,k,ncalcv_max

  allocate(a_all(xnitermax,nkstot))
  allocate(b_all(xnitermax,nkstot))
  allocate(xnorm_all(1,nkstot))
  allocate(ncalcv_all(1,nkstot))
  allocate(aux(xnitermax,nks))


  OPEN ( UNIT = 10, FILE = x_save_file, FORM = 'FORMATTED', &
       STATUS = 'UNKNOWN', IOSTAT = ierr )
  CALL errore( 'iosys', 'x_save_file ' // TRIM( x_save_file ) // &
       & ' not found' , ierr )
  read(10,*) lsda,nspin
  read(10,*) xm_r,nkstot_r,xnitermax
  
  if(xm_r.ne.xang_mom) then
     write(stdout,*) 'xm_r=',xm_r
     call errore('read_save_file','xm_r is different from xang_mom=',xang_mom)
  endif

  read(10,*) ncalcv_max
  if(ncalcv_max.gt.xnitermax) then
     write(stdout,*) 'ncalcv_max=',ncalcv_max
     call errore('read_save_file','ncalcv_max is grater than xnitermax=',xnitermax)
  endif

  read(10,*) core_energy

  read(10,*) (xkvec_r(i),i=1,3)
  write(stdout,*) '---------------------------------------------------------'
  write(stdout,*) 'xkvec read from savefile'
  write(stdout,*) (xkvec_r(i),i=1,3)
  read(10,*) (xepsilon_r(i),i=1,3)
  write(stdout,*) 'xepsilon read from file'
  write(stdout,*) (xepsilon_r(i),i=1,3)
  write(stdout,*) '---------------------------------------------------------'

  read(10,*) (xnorm_all(1,j),j=1,nkstot)
  read(10,*) (ncalcv_all(1,j),j=1,nkstot)
  read(10,*) ((a_all(i,k),i=1,ncalcv_max),k=1,nkstot)
  read(10,*) ((b_all(i,k),i=1,ncalcv_max),k=1,nkstot)             
  close(10)


  call poolscatter(xnitermax,nkstot,a_all,nks,aux)
  a(1:xnitermax,1,1:nks)=aux(1:xnitermax,1:nks)
  call poolscatter(xnitermax,nkstot,b_all,nks,aux)
  b(1:xnitermax,1,1:nks)=aux(1:xnitermax,1:nks)
  call poolscatter(1,nkstot,xnorm_all,nks,xnorm)
  call ipoolscatter(1,nkstot,ncalcv_all,nks,ncalcv)

  deallocate(a_all)
  deallocate(b_all)
  deallocate(xnorm_all)
  deallocate(ncalcv_all)
  deallocate(aux)
end subroutine read_save_file

subroutine read_header_save_file(x_save_file)
  USE kinds, only : DP
  USE klist,      ONLY : nkstot
  USE lsda_mod,    ONLY : nspin,lsda
  implicit none
  character(LEN=256) :: x_save_file
  integer            :: ierr,nkstot_r
  integer            :: xm_r

  OPEN ( UNIT = 10, FILE = x_save_file, FORM = 'FORMATTED', &
       STATUS = 'UNKNOWN', IOSTAT = ierr )
  CALL errore( 'iosys', 'x_save_file ' // TRIM( x_save_file ) // &
       & ' not found' , ierr )
  rewind(10)
  read(10,*) lsda,nspin
  read(10,*) xm_r,nkstot_r
  nkstot=nkstot_r
  close(10)

  return 
end subroutine read_header_save_file

subroutine write_save_file(a,b,xnorm,ncalcv,x_save_file)
  USE kinds, only : DP
  USE klist,      ONLY : nks,nkstot
  USE xanes,      ONLY : xnitermax,xang_mom,xkvec,xepsilon,xiabs
  USE io_global,       ONLY : ionode
!*apsi  USE uspp_param, ONLY : psd
  USE lsda_mod,    ONLY : nspin,lsda
  use uspp_param, only : upf
  implicit none
  character(LEN=256) :: x_save_file
  integer            :: ierr
  integer ncalcv(1,nks)
  real(dp) a(xnitermax,1,nks),b(xnitermax,1,nks)     
  real (dp)  xnorm(1,nks)
  integer, allocatable :: ncalcv_all(:,:)
  real(dp), allocatable :: a_all(:,:),b_all(:,:),xnorm_all(:,:)
  real (dp) :: mygetK
  integer i,j,k,ncalcv_max

  allocate(a_all(xnitermax,nkstot))
  allocate(b_all(xnitermax,nkstot))
  allocate(xnorm_all(1,nkstot))
  allocate(ncalcv_all(1,nkstot))
  
  ncalcv_all(1,1:nks)=ncalcv(1,1:nks)
  xnorm_all(1,1:nks)=xnorm(1,1:nks)
  a_all(1:xnitermax,1:nks)=  a(1:xnitermax,1,1:nks)
  b_all(1:xnitermax,1:nks)=  b(1:xnitermax,1,1:nks)

  call poolrecover(a_all,xnitermax,nkstot,nks)
  call poolrecover(b_all,xnitermax,nkstot,nks)
  call poolrecover(xnorm_all,1,nkstot,nks)
  call ipoolrecover(ncalcv_all,1,nkstot,nks)
  
  ncalcv_max=0
  do i=1,nkstot
     if(ncalcv_all(1,i).gt.ncalcv_max) ncalcv_max=ncalcv_all(1,i)
  enddo
  if ( ionode ) then
     OPEN ( UNIT = 10, FILE = x_save_file, FORM = 'FORMATTED', &
          STATUS = 'UNKNOWN', IOSTAT = ierr )
     write(10,*) lsda,nspin
     write(10,*) xang_mom,nkstot,xnitermax
     write(10,*) ncalcv_max
!*apsi     write(10,*) mygetK(psd(xiabs))
     write(10,*) mygetK(upf(xiabs)%psd)
     write(10,*) (xkvec(i),i=1,3)
     write(10,*) (xepsilon(i),i=1,3)
     write(10,*) (xnorm_all(1,j),j=1,nkstot)
     write(10,*) (ncalcv_all(1,j),j=1,nkstot)
     write(10,*) ((a_all(i,k),i=1,ncalcv_max),k=1,nkstot)
     write(10,*) ((b_all(i,k),i=1,ncalcv_max),k=1,nkstot)             
     close(10)
  endif

  deallocate(a_all)
  deallocate(b_all)
  deallocate(xnorm_all)
  deallocate(ncalcv_all)


end subroutine write_save_file



subroutine write_status_of_the_code
USE io_global,       ONLY : stdout
implicit none

  write (stdout,*) '======= Working features (22/02/2008) ==========='
  write (stdout,*) 'xanes works both in dipolar and quadrupolar part,'
  write (stdout,*) 'spin polarized works'
  write (stdout,*) 'DFT+U implemented, validated for norm-conserving '
  write (stdout,*) 'cut occupied states working, improved'
  write (stdout,*) 'terminator working'
  write (stdout,*) 'Multiprojectors TM+USPP working (MCB,CG)'
  write (stdout,*) 'ultrasoft not implemented'
  write (stdout,*) 'DFT+U tested only for non ortho wfc, but implemented'
  write (stdout,*) '======= TO DO                         ==========='
  write (stdout,*) 'Bethe-Salpeter [CG] '
  write (stdout,*) 'RXES [DC] o [CG]' 
  write (stdout,*) '================================================='

end subroutine write_status_of_the_code

!<CG>
subroutine verify_hpsi
  USE io_files,         ONLY : nd_nmbr, prefix, tmp_dir, &
       nwordwfc, iunwfc
  USE io_global,        ONLY : stdout     ! Modules/io_global.f90
  USE kinds,            only : DP
  use parameters,       ONLY : ntypx
  USE radial_grids,     ONLY : ndmx
  USE ions_base,        ONLY : nat, ntyp => nsp, ityp
  use wvfct,            ONLY : npwx, nbndx, nbnd, npw, igk, g2kin, et
  use lsda_mod,         ONLY : nspin,lsda,isk,current_spin
  use cell_base,        only: tpiba2, bg
  use wavefunctions_module, only: evc
  use klist,            ONLY : &
       nkstot,            & ! total number of k-points
       nks,               & ! number of k-points per pool
       xk,                & ! k-points coordinates
       wk                   ! k-points weight
  use gvect,            ONLY: g,ngm,ecutwfc,ngl,nrxx
  use paw_gipaw,        ONLY : &
       paw_vkb,             & ! |p> projectors
       paw_becp,            & ! product of projectors and wf.
       paw_nkb,             & ! total number of beta functions, with st.fact.
       paw_lmaxkb,paw_recon
  use becmod,     ONLY : becp,rbecp, calbec,allocate_bec, deallocate_bec !CG
  use scf,        ONLY : vltot, vrs, v, kedtau !CG
  use gsmooth,    ONLY : doublegrid
  USE mp_global,  ONLY : intra_pool_comm, mpime,my_pool_id, npool
  USE mp,         ONLY : mp_sum
  USE xanes,      only : xiabs, xanes_dip, xang_mom, xniter, xnitermax, xepsilon
  USE atom,       ONLY : rgrid, msh
  use radin_mod
  USE uspp,   ONLY : vkb, nkb, okvan
  USE ldaU,   ONLY : lda_plus_u
  use xanes_paw_variables, only : xanes_paw_nhm
  use wavefunctions_module, only: evc
  USE uspp_param, ONLY : upf

  implicit none
  integer :: is,ik,iabso,nr,ip,jp,l,j,icrd,ip_l,nrc,nt,na,i
  integer :: ipx,ipx_0,ipy,ipz,nline,nrest,npw_partial
  real (dp) v_of_0
  real (dp) norm
  complex(kind=DP) :: ZDOTC
  complex(dp), allocatable :: paw_vkb_cplx(:,:)
  real(dp) :: normps

  external ZDOTC
  external ZDSCAL

  complex(dp), allocatable :: psiwfc(:)

  integer :: ipw
  complex(dp) :: prodscal, hevc(npwx), vecteur(npwx), prodscal_part, normtemp
  integer :: indice, numk, nkppool, ind2, rest, mpimea, mpimeb
  logical :: exst, opnd
  character ( len=32) :: filehpsi, filenumber
  real(dp) :: maxdiff
  complex(dp) :: psi_h_psi,  psi_psi,  psi_s_psi, psi_sm1s_psi
  real(dp) :: difference

  call set_vrs(vrs,vltot,v%of_r,kedtau, v%kin_r,nrxx,nspin,doublegrid)
  if (lda_plus_u) call init_xanes_ldau

  mpimea=mpime
  filenumber=''


  do  j=1, 3
     mpimeb=mod(mpimea,10)
     filenumber=char(mpimeb+48)//filenumber
     mpimea=mpimea/10
  enddo

  filehpsi='hpsi_'//trim(filenumber)//'.out'
  open(unit=1000+mpime,file=filehpsi,form='formatted',status='unknown')
  rewind(1000+mpime)

  maxdiff=0.d0

  do ik=1,nks
     if(lsda) current_spin=isk(ik)
     !gk_sort  sort k-points and exit kinetic energies
     call gk_sort(xk (1,ik),ngm,g,ecutwfc/tpiba2,npw,igk,g2kin)  !CHECK
     g2kin=g2kin*tpiba2                                          !CHECK
     npw_partial = npw
     call mp_sum( npw_partial, intra_pool_comm )
     call init_gipaw_2(npw,igk,xk(1,ik),paw_vkb)
     call init_us_2(npw,igk,xk(1,ik),vkb)
     if (lda_plus_u) call init_xanes_ldau_2(ik)

     write(1000+mpime,*) 'lda_plus_u=', lda_plus_u
     write(1000+mpime,*) 'mypoolid=', my_pool_id,' ik=', ik
     write(1000+mpime,*) 'npool=', npool, ' nkstot=', nkstot
     write(1000+mpime,*) 'nwordwfc=', nwordwfc, ' iunwfc=', iunwfc

     !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     ! open saved eigenvectors
     !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     inquire (unit = iunwfc, opened = opnd)
     if (.not.opnd)   CALL diropn( iunwfc, 'wfc', nwordwfc, exst )
     call davcio( evc, nwordwfc, iunwfc, ik, -1 )

     numk=0
     nkppool=nkstot/npool
     rest=nkstot-nkppool*npool
     do indice=0, my_pool_id-1
        numk=numk+nkppool
        if (indice<rest) numk=numk+1
     enddo
     numk=numk+ik
     write(1000+mpime,*) 'npw=', npw, 'npwx=', npwx
     !write(1000+mpime,*) 'ns(1,1,1,1)=', ns(1,1,1,1)
     !write(1000+mpime,*) 'natomwfc=', natomwfc
     !write(1000+mpime,*) 'swfcatom=', ((swfcatom(ind2, indice), ind2=1,3), indice=1, natomwfc)

     do indice=1, nbnd
        ! calculating <psi | psi>
        psi_psi=ZDOTC(npw,evc(:,indice),1,evc(:,indice),1)
        call mp_sum( psi_psi, intra_pool_comm )
        ! calculating <psi | H | psi >
        hevc(:)=(0.d0,0.d0)
        CALL h_psi( npwx, npw,1, evc(:,indice), hevc )
        psi_h_psi=ZDOTC(npw,evc(:,indice),1,hevc,1)
        call mp_sum( psi_h_psi, intra_pool_comm )
           difference=abs(psi_h_psi-et(indice,numk)*psi_psi)
           write(1000+mpime,*) 'k-point', ik, 'absolute k =', numk
           write(1000+mpime,*) 'bande : ', indice
           write(1000+mpime,*) 'current spin=', current_spin
           write(1000+mpime,*) 'et(bande,ik)=', et(indice, numk)
           write(1000+mpime,*) '|<psi|psi>|=',abs(psi_psi)
           write(1000+mpime,*) '|<psi|H|psi>|=',abs(psi_h_psi)
           write(1000+mpime,*) '|<psi|H|psi>-E*<psi|psi>|=', difference
        if (difference > maxdiff)  maxdiff=difference
        if (difference > 1.d-4) write(1000+mpime,*) 'warning : difference too big'

        write(1000+mpime,*) ' '
     enddo
     write(1000+mpime,*) '--------------------- end ------------------ '
  enddo  !on k points

  close(1000+mpime)

  write(stdout,*) '> h_psi test : maximum difference is ', maxdiff

end subroutine verify_hpsi







subroutine read_gamma_file
  use gamma_variable_mod
  implicit none
  integer :: nl, ierr, i

  open (unit=21, file=gamma_file, form='formatted',status='unknown', iostat=ierr)
  call errore ('io ', 'gamma file '//trim(gamma_file)//' not found', abs (ierr) )
  rewind(21)

  nl=0

  do
     read (21,'(a1)',iostat=ierr)
     if (ierr.ne.0) exit
     nl=nl+1
  enddo
  close(21)

  gamma_lines=nl
  allocate(gamma_points(nl,2))

  open (unit=21, file=gamma_file, form='formatted',status='unknown', iostat=ierr)
  rewind(21)

  do i=1,nl
     read(21,*) gamma_points(i,1), gamma_points(i,2)
  enddo

  close(21)


end subroutine read_gamma_file

subroutine initialize_gamma_tab
  use xanes, only : xemin, xemax, xnepoint
  use kinds, only :dp
  use io_global, only : stdout
  use gamma_variable_mod
  implicit none
  real(dp) :: e,x,y,dx
  integer :: i,j,n
  real (dp), parameter :: ryd2ev = 13.6058d0
  dx=(xemax-xemin)/dfloat(xnepoint)

  do n=1, xnepoint
     x=xemin+(n-1)*dx
     i=1
     do j=1, gamma_lines
        if(x>gamma_points(j,1)) i=i+1
     enddo

     if (i.eq.1) then
        y=gamma_points(1,2)
     elseif (i.eq.(gamma_lines+1)) then
        y=gamma_points(gamma_lines,2)
     else
        y=(gamma_points(i-1,2)*(gamma_points(i,1)-x)+gamma_points(i,2)*(x-gamma_points(i-1,1)))&
             /(gamma_points(i,1)-gamma_points(i-1,1))
     endif
     gamma_tab(n)=y/ryd2ev
  enddo

end subroutine initialize_gamma_tab

subroutine determine_polycut(t1,t2,f1,f2,df1,df2,poly)
  ! calculates the interpolation polynome betwenn 2 points
  use kinds, only : dp
  implicit none
  real(dp) :: t1,t2,f1,f2,df1,df2,poly(4)

  poly(4)=((t2-t1)*(df2+df1)-2*(f2-f1))/((t2-t1)**3)
  poly(3)=(df2-df1)/(2*(t2-t1))-1.5d0*(t2+t1)*((t2-t1)*(df2+df1)-2*(f2-f1))/((t2-t1)**3)
  poly(2)=df1-2*t1*poly(3)-3*t1*t1*poly(4)
  poly(1)=f1-poly(2)*t1-poly(3)*t1**2-poly(4)*t1**3

end subroutine determine_polycut

!</CG>

