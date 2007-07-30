!
! Copyright (C) 2004 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!---------------------------------------------------------------
subroutine ld1_readin
  !---------------------------------------------------------------
  !
  !     This routine reads the input parameters of the calculation
  !
  USE io_global,  ONLY : ionode, ionode_id
  USE mp,         ONLY : mp_bcast
  use ld1inc
  use funct, only : set_dft_from_name
  use atomic_paw, only : paw_io, paw2us
  implicit none

  integer ::  &
       n,i,   &          ! counter on wavefunctions
       nc,    &          ! counter on configuration
       ns,ns1,&          ! counter on pseudo wavefunctions
       c1,    &          ! counter
       ios               ! I/O control

  real(DP) :: &
       edum(nwfsx), zdum ! auxiliary variables
  character(len=6)  :: zdum_

  character(len=80) :: config, configts(ncmax1)
  character(len=2)  :: atom
  character(len=20) :: dft
  character, external :: atom_name*2
  integer, external :: atomic_number
  logical, external :: matches

  namelist /input/ xmin,    &  ! the minimum x of the linear mesh
       dx,      &  ! parameters of the mesh
       rmax,    &  ! the maximum r of the mesh
       zed,     &  ! the atomic charge
       atom,    &  ! atomic symbol - can be specified instead of zed
       beta,    &  ! the mixing coefficient
       tr2,     &  ! the scf threshold
       iswitch, &  ! the type of calculation
       nld, rlderiv, eminld, emaxld, deld,& ! log derivatives
       config,  &  ! a string with electron configuration
       lsd,     &  ! if 1 lsda is computed      
       rel,     &  ! 0 non-relativistic calculation
                   ! 1 scalar-relativistic calculation
                   ! 2 dirac-relativistic calculation
       dft,     &  ! LDA, GGA, exchange only or Hartree ?
       isic,    &  ! if 1 self-interaction correction
       latt,    &  ! if <> 0 Latter correction is applied
       title,   &  ! the title of the run
       prefix,  &  ! the prefix for file names
       vdw         ! if .true. vdW coefficient in TF+vW will be calculated

  namelist /test/                 &
       nconf,         & ! the number of configurations
       configts,      & ! the configurations of the tests
       file_pseudo,   & ! input file containing the pseudopotential
       file_pseudopw, & ! output file containing the pseudopotential
       ecutmin,       & ! for test with spherical Bessel functions:
       ecutmax,       & ! min and max energy cutoff for j_l(qr),
       decut,         & ! step: ecut = ecutmin, ecutmin+decut, ... , ecutmax
       rm               ! radius of the box

  namelist /inputp/ &
       pseudotype,&! the pseudopotential type
       tm,    &    ! use Troullier-Martins instead of RRKJ
       rho0,  &    ! value of the charge at the origin
       zval,  &    ! the pseudo valence
       lloc,  &    ! l component considered as local 
       nlcc,  &    ! if true nlcc is set
       rcore, &    ! the core radius for nlcc
       rcloc, &    ! the local cut-off for pseudo
       lpaw,  &    ! if true create a PAW dataset
       author, &   ! the author of the PP
       file_pseudopw, & ! output file where the pseudopotential is written
       file_screen,   & ! output file for the screening potential
       file_core,     & ! output file for total and core charge
       file_beta,     & ! output file for the beta functions
       file_chi,      & ! outpu  file for the chi functions
       file_qvan,     & ! output file for the qvan functions
       file_wfcaegen, & ! output file where the all-electron wfc used for 
                        !        pseudo generation are written
       file_wfcncgen, & ! output file where the norm-conserving wfc used for 
                        !        pseudo generation are written
       file_wfcusgen, & ! output file where the ultra-soft wfc used for 
                        !        pseudo generation are written
       file_recon       ! output file needed for the paw reconstruction
       
   !
  prefix       = 'ld1'
  file_pseudo  = ' '
  file_pseudopw= ' '
  file_recon   = ' '
  file_screen  = ' '
  file_core    = ' '
  file_chi     = ' '
  file_beta    = ' '
  file_qvan    = ' '
  file_wfcaegen = ' '
  file_wfcncgen = ' '
  file_wfcusgen = ' '
  !
  !   set default values 
  !
  atom  = '  '
  zed   = 0.0_dp
  xmin  = -7.0_dp
  dx    =  0.0125_dp
  rmax  =100.0_dp

  beta  =  0.2_dp
  tr2   = 1.0e-14_dp
  iswitch=1

  rlderiv=4.0_dp
  eminld=-3.0_dp
  emaxld=3.0_dp
  nld=0
  deld=0.03_dp

  rel = 5 
  lsd = 0
  dft = 'LDA'
  isic= 0
  latt= 0
  title = ' '
  config= ' '

  lpaw = .false.
  author='anonymous'

  vdw  = .false.
  
  ! read the namelist input

  if (ionode) read(5,input,err=100,iostat=ios) 
100  call mp_bcast(ios, ionode_id)
  call errore('ld1_readin','reading input namelist ',abs(ios))
  call bcast_input()
  call mp_bcast(atom, ionode_id )
  call mp_bcast(config, ionode_id )
  call mp_bcast(dft, ionode_id )

  call set_dft_from_name(dft)

  if (zed == 0.0_dp .and. atom /= ' ') then
     zed = DBLE(atomic_number(atom))
  else if (zed /= 0.0_dp .and. atom == ' ') then
     if (DBLE(int(zed)) /= zed .or. zed < 1.0_dp .or. zed > 100) then
        write(zdum_,'(f6.2)') zed
        call errore('ld1_readin','wrong nuclear charge zed: '//zdum_,1)
     end if
     atom = atom_name(nint(zed))
  else
     zdum = DBLE(atomic_number(atom))
     if (zdum /= zed) call errore &
          ('ld1_readin','inconsistent Z/atom specification',nint(zdum))
  end if
  if (iswitch < 1 .or. iswitch > 3) &
       call errore('ld1_readin','wrong iswitch',1)
  if (eminld > emaxld) &
       call errore('ld1_readin','eminld or emaxld wrong',1)
  if (deld < 0.0_dp) &
       call errore('ld1_readin','negative deld',1)
  if (nld > nwfsx) &
       call errore('ld1_readin','too many nld',1)
  if (xmin > -2) call errore('ld1_readin','wrong xmin',1)
  if (dx <=0.0_dp) call errore('ld1_readin','wrong dx',1)

  if (isic == 1 .and. latt == 1) call errore('ld1_readin', &
       &    'isic and latter correction not allowed',1)
  if (isic == 1 .and. iswitch .ne. 1 ) call errore('ld1_readin', &
       &    'SIC available with all-electron only', 1)
 

  zmesh=zed
  if (rel == 5 ) then
     if (zed < 19.0_dp) then
        rel=0
     else
        rel=1
     endif
  endif
  if (rel < 0 .or. rel > 2) call errore('ld1_readin','wrong rel',1)
  !
  !     No lsda with pseudopotential generation
  !
  if (iswitch > 2) lsd = 0
  if (lsd == 0) then
     nspin = 1
  else if(lsd == 1) then
     nspin = 2
     if (rel == 2) call errore('ld1_readin', &
       &    'local spin density and spin-orbit not allowed',1)
  else
     call errore('ld1_readin','lsd not correct',1)
  endif

  if (config == ' ') then
     if (ionode) call read_config (rel, lsd, nwf, el, nn, ll, oc, isw, jj)
     call bcast_config()
  else
     call el_config (config, rel, lsd, .true., nwf, el, nn, ll, oc, isw, jj)
  end if
  !
  !  In the spin polarized or relativistic case adjust the occupations
  !
  if (lsd == 1) then
     call occ_spin(nwf,nwfx,el,nn,ll,oc,isw)
  else if (rel == 2) then
     call occ_spinorb(nwf,nwfx,el,nn,ll,jj,oc,isw)
  endif
  !
  ! generate the radial grid - note that if iswitch = 2 the radial grid
  ! is not generated but read from the pseudopotential file
  !
  if (iswitch /= 2) then
     call do_mesh(rmax,zmesh,xmin,dx,0,ndm,mesh,r,r2,rab,sqr)
     rhoc=0.0_dp
  endif
  !
  if (iswitch == 1) then
     !
     !    no more data needed for AE calculations
     !
     return
     !     
  else if (iswitch == 3) then
     !
     !    reading input for PP generation
     !
     zval=0.0_dp
     lloc=-1
     rcloc=-1_dp
     nlcc=.false.
     rcore=0.0_dp
     rho0=0.0_dp
     tm  = .false.
     pseudotype=0
     jjs=0.0_dp

     if (ionode) read(5,inputp,err=500,iostat=ios)
500  call mp_bcast(ios, ionode_id)
     call errore('ld1_readin','reading inputp',abs(ios))
     call bcast_inputp()

     if (lloc < 0 .and. rcloc <=0.0_dp) &
          call errore('ld1_readin','rcloc must be positive',1)
     if (pseudotype < 1.or.pseudotype > 3) &
          call errore('ld1_readin','specify correct pseudotype',1)
     if (rel==2 .and. pseudotype==1 ) &
          call errore('ld1_readin','Generation of a FR PP with'// & 
                  &     ' pseudotype=1 not allowed',1)

     !
     if (ionode) &
        call read_psconfig (rel, lsd, nwfs, els, nns, lls, ocs, &
             isws, jjs, enls, rcut, rcutus )
     call bcast_psconfig()
     !
     if (rel==2) call occ_spinorbps &
          (nwfs,nwfsx,els,nns,lls,jjs,ocs,rcut,rcutus,enls,isws)
     !
     lmax = maxval(lls(1:nwfs))
     !
     zdum = zed
     do n=1,nwf
        if ( oc(n) > 0.0_dp) zdum = zdum - oc(n)
     end do
     do ns=1,nwfs
        if ( ocs(ns) > 0.0_dp) zdum = zdum + ocs(ns)
     end do
     if (zval == 0) then
        zval = zdum
        if ( abs(nint(zdum)-zdum) > 1.d-8 ) then
           write(zdum_,'(f6.2)') zdum
           call errore ('ld1_readin', 'found noninteger valence ' &
                &//zdum_//', if you want this specify zval in inputp',1)
        end if
     else if ( abs(zval-zdum) > 1.d-8 ) then
        call errore ('ld1_readin',&
             ' supplied and calculated valence charge do not match',1)
     end if
     !
     do ns=1,nwfs
        if (pseudotype < 3) rcutus(ns) = rcut(ns)
        do ns1=1,ns-1
           if (lls(ns) == lls(ns1).and.pseudotype == 1) &
                call errore('ld1_readin','two wavefunctions for same l',1)
        enddo
        !
        if (enls(ns) /= 0.0_dp .and. ocs(ns) > 0.0_dp) &
             call errore('ld1_readin','unbound states must be empty',1)
        if (rcut(ns) /= rcutus(ns)) then
           !
           ! this channel is US. Check that there is at least another energy
           !
          c1=0
          do ns1=1,nwfs
             if (els(ns) == els(ns1) .and. jjs(ns) == jjs(ns1)) c1=c1+1 
          enddo
          if (c1 < 2) call errore('ld1_readin', &
                        'US requires at least two energies per channel',1)
        endif
     enddo
     if (nwfs > 1) then
        if (els(nwfs)==els(nwfs-1) .and. jjs(nwfs)==jjs(nwfs-1) .and. &
            lloc > -1) call errore('ld1_readin','only one local channel',1)
     endif
     nlc=0
     nnl=0

  end if
  !
  !    reading input for PP testing
  !
  jjts=0.0_dp
  jjtsc=0.0_dp
  
  nconf=1
  configts=' '
  ecutmin = 0.0_dp
  ecutmax = 0.0_dp
  decut   = 5.0_dp
  rm      =30.0_dp

  if (ionode) read(5,test,err=300,iostat=ios)
300  call mp_bcast(ios, ionode_id)

  if (iswitch==2) call errore('ld1_readin','reading test',abs(ios))
  call bcast_test()
  call mp_bcast(configts, ionode_id)
  !
  !  PP generation: if namelist test is not found, use defaults
  !
  if (iswitch == 3 .and. ios /= 0 ) then
     !
     ! use for testing the same configuration as for PP generation
     ! (unless a different one is explicitely specified in namelist &test)
     !
     ns1 = 0
     do ns=1,nwfs
        !
        if ( ocs(ns)  > 0.0_dp .or. &
            (ocs(ns) == 0.0_dp .and. enls(ns) == 0.0_dp) ) then
           !
           ! copy states used in the PP generation to testing configuration
           ! Only bound states must be copied. Note that this WILL NOT WORK
           ! if bound states are not used in the generation of the PP
           !
           ns1 = ns1 + 1
           eltsc (ns1,1)= els (ns)
           nntsc (ns1,1)= nns (ns)
           lltsc (ns1,1)= lls (ns)
           octsc (ns1,1)= ocs (ns)
           iswtsc(ns1,1)= isws(ns)
           jjtsc (ns1,1)= jjs (ns)
        end if
     end do
     !
     nwftsc(1) = ns1
     !
     return
     !
  endif
  !
  if (nconf > ncmax1.or.nconf < 1) &
       call errore('ld1_readin','nconf is wrong',1)
  if (iswitch == 3 .and. nconf > 1) &
       call errore('ld1_readin','too many test configurations',1)
  !  
  do nc=1,nconf
     if (configts(nc) == ' ') then
        if (ionode) &
        call read_psconfig (rel, lsd, nwftsc(nc), eltsc(1,nc), &
             nntsc(1,nc), lltsc(1,nc), octsc(1,nc), iswtsc(1,nc), &
             jjtsc(1,nc), edum(1), rcuttsc(1,nc), rcutustsc(1,nc) )
        call mp_bcast(eltsc(:,nc),ionode_id)
     else
        call el_config( configts(nc), rel, lsd, .false., nwftsc(nc), &
             eltsc(1,nc), nntsc(1,nc), lltsc(1,nc), octsc(1,nc), &
             iswtsc(1,nc), jjtsc(1,nc))
     endif
  enddo
  call bcast_pstsconfig()
  do nc=1,nconf
     !
     !  adjust the occupations of the test cases if this is a lsd run
     !
     if (lsd == 1) then
        call occ_spin(nwftsc(nc),nwfsx,eltsc(1,nc),nntsc(1,nc),lltsc(1,nc),&
             octsc(1,nc), iswtsc(1,nc)) 
     else if (rel == 2) then
        call occ_spinorb(nwftsc(nc),nwfsx,eltsc(1,nc), &
             nntsc(1,nc),lltsc(1,nc),jjtsc(1,nc),octsc(1,nc),iswtsc(1,nc))
     endif
  end do
  !
  !    PP testing: reading the pseudopotential
  !
  if (iswitch ==2) then
     lpaw=.false.
     !
     if (file_pseudo == ' ') &
       call errore('ld1_readin','file_pseudo is needed',1)
     if (matches('.upf',file_pseudo) .or. matches('.UPF', file_pseudo)) then
        !
        !    UPF format
        !
        call read_pseudoupf
        !
     else if (matches('.PAW',file_pseudo) .or. matches('.PAW',file_pseudo)) then
        !
        !    PAW dataset
        !
        lpaw=.true.
        open(unit=111, file=trim(file_pseudo), status='unknown',  &
             form='formatted', err=50, iostat=ios)
50      call errore('ld1_readin','open error on file '//file_pseudo,abs(ios))
        call paw_io(pawsetup,111,"INP")
        close(111)
        call paw2us ( pawsetup, zval, mesh, r, r2, sqr, dx, nbeta, lls, &
             ikk, betas, qq, qvan, pseudotype )
        !
     else if ( matches('.rrkj3', file_pseudo) .or. &
               matches('.RRKJ3', file_pseudo)) then
        !
        !    Old RRKJ format
        !
        call read_newpseudo (ios)
        !
        if (ios /= 0) then
           !
           !    try old Norm-Conserving format
           !
           pseudotype = 1
           !
        else
           lmax=0
           do ns=1,nwfs
              lmax=max(lmax,lls(ns))
           enddo
        end if
     else
        !
        !    Old Norm-Conserving format
        !
        pseudotype = 1
        !
     endif
     !
     if (pseudotype == 1) then
        !
        call read_pseudo  &
             (file_pseudo,zed,xmin,rmax,dx,mesh,ndm,r,r2,rab,sqr, &
             dft,lmax,lloc,zval,nlcc,rhoc,vnl,vpsloc,rel)
        call set_dft_from_name(dft)
        !
        do ns=1,lmax+1
           ikk(ns)=mesh
        enddo
     endif
     !
  endif
  !
  if (lpaw) then
     if (pseudotype /= 3) call errore('ld1_readin', &
          'please start from a US for generating a PAW dataset' ,pseudotype)
     if (rel /= 0) call errore('ld1_readin', &
          'relativistic PAW not implemented' ,rel)
     if (latt /= 0) call errore('ld1_readin', &
          'Latter correction not implemented in PAW' ,latt)
     call errore('ld1_readin', &
          'PAW dataset generation and test is experimental' ,-1)
  endif

  return

end subroutine ld1_readin

subroutine bcast_input()
  USE io_global,  ONLY : ionode_id
  USE mp,         ONLY : mp_bcast
  use ld1inc

implicit none
#ifdef __PARA
   call mp_bcast( xmin, ionode_id )
   call mp_bcast( dx, ionode_id )
   call mp_bcast( rmax, ionode_id )
   call mp_bcast( zed, ionode_id )
   call mp_bcast( beta, ionode_id )
   call mp_bcast( tr2, ionode_id )
   call mp_bcast( iswitch, ionode_id )
   call mp_bcast( nld, ionode_id )
   call mp_bcast( rlderiv, ionode_id )
   call mp_bcast( eminld, ionode_id )
   call mp_bcast( emaxld, ionode_id )
   call mp_bcast( deld, ionode_id )
   call mp_bcast( lsd, ionode_id )
   call mp_bcast( rel, ionode_id )
   call mp_bcast( isic, ionode_id )
   call mp_bcast( latt, ionode_id )
   call mp_bcast( title, ionode_id )
   call mp_bcast( prefix, ionode_id )
   call mp_bcast( vdw, ionode_id )
#endif
return
end subroutine bcast_input

subroutine bcast_inputp()
  USE io_global,  ONLY : ionode_id
  USE mp,         ONLY : mp_bcast
  use ld1inc

implicit none
#ifdef __PARA
  call mp_bcast( pseudotype, ionode_id )
  call mp_bcast( tm,  ionode_id ) 
  call mp_bcast( rho0,  ionode_id )
  call mp_bcast( zval,  ionode_id )
  call mp_bcast( lloc,  ionode_id )
  call mp_bcast( nlcc,  ionode_id )
  call mp_bcast( rcore, ionode_id )
  call mp_bcast( rcloc, ionode_id )
  call mp_bcast( lpaw,  ionode_id )
  call mp_bcast( file_pseudopw, ionode_id )
  call mp_bcast( file_screen, ionode_id ) 
  call mp_bcast( file_core, ionode_id )
  call mp_bcast( file_beta, ionode_id )
  call mp_bcast( file_chi, ionode_id )
  call mp_bcast( file_qvan, ionode_id )
  call mp_bcast( file_wfcaegen, ionode_id )
  call mp_bcast( file_wfcncgen, ionode_id )
  call mp_bcast( file_wfcusgen, ionode_id )
  call mp_bcast( file_recon,  ionode_id )
#endif
  return
end subroutine bcast_inputp

subroutine bcast_test()
  USE io_global,  ONLY : ionode_id
  USE mp,         ONLY : mp_bcast
  use ld1inc

implicit none
#ifdef __PARA
   call mp_bcast( nconf, ionode_id ) 
   call mp_bcast( file_pseudo, ionode_id )
   call mp_bcast( ecutmin, ionode_id ) 
   call mp_bcast( ecutmax, ionode_id ) 
   call mp_bcast( decut, ionode_id ) 
   call mp_bcast( rm, ionode_id )      
#endif
return
end subroutine bcast_test

subroutine bcast_config()
  USE io_global,  ONLY : ionode_id
  USE mp,         ONLY : mp_bcast
  use ld1inc

implicit none
#ifdef __PARA
  call mp_bcast( nwf, ionode_id )
  call mp_bcast( el, ionode_id )
  call mp_bcast( nn, ionode_id )
  call mp_bcast( ll, ionode_id )
  call mp_bcast( oc, ionode_id )
  call mp_bcast( isw, ionode_id )
  call mp_bcast( jj, ionode_id )
#endif
return
end subroutine bcast_config

subroutine bcast_psconfig()
  USE io_global,  ONLY : ionode_id
  USE mp,         ONLY : mp_bcast
  use ld1inc

implicit none
#ifdef __PARA
  call mp_bcast( nwfs, ionode_id )
  call mp_bcast( els, ionode_id )
  call mp_bcast( nns, ionode_id )
  call mp_bcast( lls, ionode_id )
  call mp_bcast( ocs, ionode_id )
  call mp_bcast( jjs, ionode_id )
  call mp_bcast( isws, ionode_id )
  call mp_bcast( enls, ionode_id )
  call mp_bcast( rcut, ionode_id )
  call mp_bcast( rcutus, ionode_id )
#endif
return
end subroutine bcast_psconfig

subroutine bcast_pstsconfig()
  USE io_global,  ONLY : ionode_id
  USE mp,         ONLY : mp_bcast
  use ld1inc

implicit none
#ifdef __PARA
  call mp_bcast( nwftsc, ionode_id )
  call mp_bcast( nntsc, ionode_id )
  call mp_bcast( lltsc, ionode_id )
  call mp_bcast( octsc, ionode_id )
  call mp_bcast( jjtsc, ionode_id )
  call mp_bcast( iswtsc, ionode_id )
  call mp_bcast( rcuttsc, ionode_id ) 
  call mp_bcast( rcutustsc, ionode_id )
#endif
return
end subroutine bcast_pstsconfig



