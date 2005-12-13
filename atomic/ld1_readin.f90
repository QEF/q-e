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
       edum(nwfsx), zdum        ! auxiliary

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
       file_pseudo      ! input file containing the pseudopotential

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
       file_pseudopw, & ! output file where the pseudopotential is written
       file_screen,   & ! output file for the screening potential
       file_core,     & ! output file for total and core charge
       file_beta,     & ! output file for the beta functions
       file_chi,      & ! outpu  file for the chi functions
       file_qvan,     & ! output file for the qvan functions
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
  lsd   = 0
  dft= 'LDA'
  latt  = 0
  title = ' '
  config=' '

  lpaw = .false.

  vdw  = .false.

  ! read the namelist input

  read(5,input,err=100,iostat=ios) 
100 call errore('ld1_readin','reading input namelist ',abs(ios))

  call set_dft_from_name(dft)

  if (zed == 0.0_dp .and. atom /= ' ') then
     zed = DBLE(atomic_number(atom))
  else if (zed /= 0.0_dp .and. atom == ' ') then
     if(DBLE(int(zed)) /= zed .or. zed < 1.0_dp .or. zed > 100) &
          call errore('ld1_readin','wrong zed',1)
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
     call read_config (rel, lsd, nwf, el, nn, ll, oc, isw, jj)
  else
     call el_config(config,.true.,nwf,el,nn,ll,oc,isw)
     !
     ! check same labels corresponding to different spin or j value
     !
     jj(1:nwf)=0.d0
     do n=1,nwf 
        do i=n+1,nwf  
           if (el(i) == el(n)) then
              if (rel == 2) then
                 if (ll(n) > 0) then
                    jj(n) = ll(n) + (isw(n)-1.5)
                    jj(i) = ll(i) + (isw(i)-1.5)
                    if ( oc(n) > (2.0_dp*jj(n)+1.0_dp) ) &
                         call errore('ld1_readin','occupation wrong',n)
                    if ( oc(i) > (2.0_dp*jj(i)+1.0_dp) ) &
                         call errore('ld1_readin','occupation wrong',i)
                 else
                    call errore('ld1_readin',el(i)//' appears twice',i)
                 end if
              else if ( lsd==0 ) then
                 call errore('ld1_readin',el(i)//' appears twice',i)
              end if
           end if
        end do
     end do
     if (rel == 2) isw(1:nwf)=1 
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

     read(5,inputp,err=500,iostat=ios)
500  call errore('ld1_readin','reading inputp',abs(ios))

     if (lloc < 0 .and. rcloc <=0.0_dp) &
          call errore('ld1_readin','rcloc must be positive',1)
     if (pseudotype < 1.or.pseudotype > 3) &
          call errore('ld1_readin','specify correct pseudotype',1)
     !
     call read_psconfig (rel, lsd, nwfs, els, nns, lls, ocs, &
          isws, jjs, enls, rcut, rcutus )
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
     if ( abs(nint(zdum)-zdum) > 1.d-8 ) call errore &
          ('ld1_readin',' calculated valence charge not integer?',1)
     if (zval == 0) then
        zval = zdum
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

  read(5,test,err=300,iostat=ios)

300 continue
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
  call errore('ld1_readin','reading test',abs(ios))
  !
  if (nconf > ncmax1.or.nconf < 1) &
       call errore('ld1_readin','nconf is wrong',1)
  if (iswitch == 3 .and. nconf > 1) &
       call errore('ld1_readin','too many test configurations',1)
  !  
  do nc=1,nconf
     if (configts(nc) == ' ') then
        call read_psconfig (rel, lsd, nwftsc(nc), eltsc(1,nc), &
             nntsc(1,nc), lltsc(1,nc), octsc(1,nc), iswtsc(1,nc), &
             jjtsc(1,nc), edum(1), rcuttsc(1,nc), rcutustsc(1,nc) )
     else
        call el_config(configts(nc),.false.,nwftsc(nc),eltsc(1,nc),  &
             &     nntsc(1,nc),lltsc(1,nc),octsc(1,nc),iswtsc(1,nc))
     endif
     !
     !  adjust the occupations of the test cases if this is a lsd run
     !
     if (lsd == 1) then
        call occ_spin(nwftsc(nc),nwfsx,eltsc(1,nc),nntsc(1,nc),lltsc(1,nc),&
             octsc(1,nc), iswtsc(1,nc)) 
     else if (rel == 2) then
        call occ_spinorb(nwftsc(nc),nwfsx,eltsc(1,nc), &
             &  nntsc(1,nc),lltsc(1,nc),jjtsc(1,nc),octsc(1,nc),iswtsc(1,nc))
     else
        jjtsc=0.0_dp
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
  end if
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
  end if

  return

end subroutine ld1_readin
!
!---------------------------------------------------------------
subroutine occ_spin(nwf,nwfx,el,nn,ll,oc,isw)
  !---------------------------------------------------------------
  !
  !  This routine splits the occupations of the states between spin-up
  !  and spin down. If the occupations are lower than 2*l+1 it does
  !  nothing, otherwise 2*l+1 states are assumed with spin up and
  !  the difference with spin down. 
  !
  use kinds, only : DP
  implicit none
  integer :: nwf, nwfx, nn(nwfx), ll(nwfx), isw(nwfx)
  real(DP) :: oc(nwfx)
  character(len=2) :: el(nwfx)

  integer :: nwf0, n, n1
  logical :: ok

  nwf0=nwf
  do n=1,nwf0
     if (oc(n) > (2*ll(n)+1)) then
        !
        !    check that the new state is not already available
        !
        do n1=n+1,nwf0
           if (el(n1)==el(n)) call errore('ld1_readin','wrong occupations',1)
        enddo
        !
        !    and add it
        !
        nwf=nwf+1
        if (nwf > nwfx) call errore('ld1_readin','too many wavefunctions',1)
        el(nwf)=el(n)
        nn(nwf)=nn(n)
        ll(nwf)=ll(n)
        oc(nwf)=oc(n)-2*ll(n)-1
        oc(n)=2*ll(n)+1
        if (isw(n) == 1) isw(nwf)=2 
        if (isw(n) == 2) isw(nwf)=1 
     else
        ok=.true.
        do n1=1,nwf0
           if (n1 /= n) ok=ok.and.(el(n1) /= el(n))  
        enddo
        if (ok) then
           nwf=nwf+1
           if (nwf > nwfx) &
                & call errore('occ_spin','too many wavefunctions',1)
           el(nwf)=el(n)
           nn(nwf)=nn(n)
           ll(nwf)=ll(n)
           oc(nwf)=0.0_dp
           if (isw(n) == 1) isw(nwf)=2 
           if (isw(n) == 2) isw(nwf)=1 
        endif
     endif
  enddo
  return
end subroutine occ_spin
!
!---------------------------------------------------------------
subroutine read_config(rel, lsd, nwf, el, nn, ll, oc, isw, jj)
  !---------------------------------------------------------------
  !
  use kinds, only: dp
  use ld1_parameters, only: nwfx
  implicit none
  ! input
  integer :: rel, lsd 
  ! output: atomic states
  character(len=2) :: el(nwfx)
  integer :: nwf, nn(nwfx), ll(nwfx), isw(nwfx)
  real(DP) :: oc(nwfx), jj(nwfx)
  ! local variables
  integer :: ios, n, ncheck
  character (len=2) :: label
  character (len=1), external :: capital
  !
  !
  read(5,*,err=200,iostat=ios) nwf
200 call errore('read_config','reading nwf ',abs(ios))
  if (nwf <= 0) call errore('read_config','nwf is wrong',1)
  if (nwf > nwfx) call errore('read_config','too many wfcs',1)
  !
  !     read the occupation of the states
  !
  do n=1,nwf  
     if (rel < 2) then
        jj(n) = 0.0_dp
        if (lsd == 0) then
           read(5,*,err=20,end=20,iostat=ios) &
                el(n), nn(n), ll(n), oc(n)
           isw(n)=1
20         call errore('read_config','reading orbital (lda)',abs(ios))
        else  
           read(5,*,err=21,end=21,iostat=ios) &
                el(n), nn(n), ll(n), oc(n), isw(n)
21         call errore('read_config','reading orbital (lsd)',abs(ios))
           if(isw(n) > 2 .or. isw(n) < 1) &
                call errore('read_config','spin variable wrong ',n)
        endif
     else
        read(5,*,err=22,end=22,iostat=ios) &
             el(n), nn(n), ll(n), oc(n), jj(n)
        isw(n)=1
        if ((abs(ll(n)+0.5_dp-jj(n)) > 1.e-3_dp) .and. &
            (abs(ll(n)-0.5_dp-jj(n)) > 1.e-3_dp) .and. abs(jj(n)) > 1.e-3_dp) &
            call errore('read_config','jj wrong',n)
        if (oc(n) > (2.0_dp*jj(n)+1.0_dp) .and. abs(jj(n)) > 1e-3_dp) &
             call errore('read_config','occupations wrong',n)
22      call errore('read_config','reading orbital (rel)',abs(ios))
     endif
     !
     ! Check: no two same wavefunctions
     !
     do ncheck=1,n-1
        if ( el(ncheck) == el(n) .and. isw(ncheck) == isw(n) .and. &
             jj(ncheck) == jj(n) ) then
           call errore('read_config', &
                'same wavefunction '//el(n)//' appears twice',n)
        endif
     enddo
     !
     ! More sanity checks
     !
     write(label,'(a2)') el(n)
     read (label,'(i1)') ncheck
     if (ncheck /= nn(n)  .or. &
         capital(label(2:2)) == 'S' .and. ll(n) /= 0 .or. &
         capital(label(2:2)) == 'P' .and. ll(n) /= 1 .or. &
         capital(label(2:2)) == 'D' .and. ll(n) /= 2 .or. &
         capital(label(2:2)) == 'F' .and. ll(n) /= 3 .or. &
         oc(n) > 2.0_dp*(2*ll(n)+1) .or. nn(n) < ll(n)+1  ) &
         call errore('read_config',label//' wrong?',n)
  enddo
  !
  return
end subroutine read_config
!
!---------------------------------------------------------------
subroutine read_psconfig (rel, lsd, nwfs, els, nns, lls, ocs, &
     isws, jjs, enls, rcut, rcutus )
  !---------------------------------------------------------------
  !
  use kinds, only: dp
  use ld1_parameters, only: nwfsx
  implicit none
  ! input
  integer :: rel, lsd 
  ! output: atomic states
  character(len=2) :: els(nwfsx)
  integer :: nwfs, nns(nwfsx), lls(nwfsx), isws(nwfsx)
  real(DP) :: ocs(nwfsx), jjs(nwfsx), enls(nwfsx), &
       rcut(nwfsx), rcutus(nwfsx)
  ! local variables
  integer :: ios, n
  character (len=2) :: label
  character (len=1), external :: capital

  read(5,*,err=600,iostat=ios) nwfs
600 call errore('read_psconfig','reading nwfs',abs(ios))

  if (nwfs <= 0 .or. nwfs > nwfsx) &
       call errore('read_psconfig','nwfs is wrong',1)

  do n=1,nwfs
     if (rel < 2) then
        if (lsd == 1) then
           read(5,*,err=30,end=30,iostat=ios) &
                els(n), nns(n), lls(n), ocs(n), enls(n), &
                rcut(n), rcutus(n), isws(n)
           if (isws(n) > 2 .or. isws(n) < 1) &
                call errore('read_psconfig', 'spin variable wrong ',n)
           if (ocs(n) > (2.0_dp*lls(n)+1.0_dp))                 &
             call errore('read_psconfig','occupations (ls) wrong',n)
        else
           read(5,*,err=30,end=30,iostat=ios) &
                els(n), nns(n), lls(n), ocs(n), enls(n), &
                rcut(n), rcutus(n)
           isws(n)=1
           if (ocs(n) > 2.0_dp*(2.0_dp*lls(n)+1.0_dp))                 &
             call errore('read_psconfig','occupations (l) wrong',n)
        end if
        jjs(n)=0.0_dp
     else
        read(5,*,err=30,end=30,iostat=ios) &
             els(n), nns(n), lls(n), ocs(n), enls(n),     &
             rcut(n), rcutus(n), jjs(n)
        isws(n)=1
        if ((abs(lls(n)+0.5_dp-jjs(n)) > 1.e-3_dp).and.      &
            (abs(lls(n)-0.5_dp-jjs(n)) > 1.e-3_dp).and. abs(jjs(n)) > 1.e-3_dp) &
             call errore('read_psconfig', 'jjs wrong',n)
        if (ocs(n) > (2.0_dp*jjs(n)+1.0_dp).and. abs(jjs(n)) > 1.e-3_dp) &
             call errore('read_psconfig','occupations (j) wrong',n)
     endif
     write(label,'(a2)') els(n)
     if ( capital(label(2:2)) == 'S'.and.lls(n) /= 0.or.   &
          capital(label(2:2)) == 'P'.and.lls(n) /= 1.or.   &
          capital(label(2:2)) == 'D'.and.lls(n) /= 2.or.   &
          capital(label(2:2)) == 'F'.and.lls(n) /= 3.or.   &
          ocs(n) > 2*(2*lls(n)+1).or.                 &
          nns(n) < lls(n)+1 )                         &
          call errore('read_psconfig','ps-label'//' wrong?',n)
     if (rcut(n) > rcutus(n)) &
          call errore('read_psconfig','rcut or rcutus is wrong',1)
  enddo
30 call errore('read_psconfig','reading pseudo wavefunctions',abs(ios))
  !
  return
end subroutine read_psconfig
!------------------------------------------------------------------------
