!
! Copyright (C) 2004 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!---------------------------------------------------------------
subroutine ld1_readin
  !---------------------------------------------------------------
  !
  !     This routine reads the input parameters of the calculation
  !
  use ld1inc
  use funct
  implicit none

  integer ::  &
       n,i,   &          ! counter on wavefunctions
       nc,    &          ! counter on configuration
       ns,ns1,&          ! counter on pseudo wavefunctions
       c1,    &          ! counter
       ios               ! I/O control

  real(kind=dp) :: &
       edum(nwfsx), zdum        ! auxiliary

  character(len=80) :: config, configts(ncmax1)
  character(len=2) :: atom
  logical, dimension(nwfsx) :: unbound
  logical :: found
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
       file_wavefunctions,& ! file names with wavefunctions
       file_logderae ! file with logder  

  namelist /test/                 &
       nconf,         & ! the number of configurations
       configts,      & ! the configurations of the tests
       file_pseudo,   & ! input file containing the pseudopotential
       file_wavefunctionsps,& ! output file for pseudowfc
       file_tests       ! output file for  transferability test

  namelist /inputp/ &
       pseudotype,&! the pseudopotential type
       tm,    &    ! use Troullier-Martins instead of RRKJ
       rho0,  &    ! value of the charge at the origin
       zval,  &    ! the pseudo valence
       lloc,  &    ! l component considered as local 
       nlcc,  &    ! if true nlcc is set
       rcore, &    ! the core radius for nlcc
       rcloc, &    ! the local cut-off for pseudo
       file_screen,   & ! output file for screening potential
       file_core,     & ! output file for total and core charge
       file_beta,     & ! output file for beta functions
       file_chi,      & ! output file for chi functions
       file_qvan,     & ! output file for qvan functions
       file_pseudopw, & ! output file where the pseudopotential is written
       file_recon,    & ! output file needed for paw reconstruction
       file_logderps    ! output file for pseudo logarithmic derivatives
  !
  !   read the namelist input and set default values 
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
  file_wavefunctions= ' '
  file_recon= ' '
  file_logderae= ' '

  read(5,input,err=100,iostat=ios) 
100 call errore('ld1_readin','reading input namelist ',abs(ios))

  call which_dft(dft)

  if (zed == 0.0_dp .and. atom /= ' ') then
     zed = dble(atomic_number(atom))
  else if (zed /= 0.0_dp .and. atom == ' ') then
     if(dble(int(zed)) /= zed .or. zed < 1.0_dp .or. zed > 100) &
          call errore('ld1_readin','wrong zed',1)
     atom = atom_name(nint(zed))
  else
     zdum = dble(atomic_number(atom))
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
     if (zed >= 19.0_dp) then
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
  else
     call errore('ld1_readin','lsd not correct',1)
  endif
  if (rel == 2 .and. lsd == 1) call errore('ld1_readin', &
       &    'local spin density and spin-orbit not allowed',1)

  if (config == ' ') then
     call read_config (rel, lsd, nwf, el, nn, ll, oc, isw, jj)
  else
     call el_config(config,.true.,nwf,el,nn,ll,oc,isw)
     !
     ! check same labels corresponding to different spin or j value
     !
     do n=1,nwf  
        do i=n+1,nwf  
           if (el(i) == el(n)) then
              if ( lsd==0 ) then
                 call errore('ld1_readin',el(i)//' appears twice',i)
              else if (rel == 2) then
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
              end if
           end if
        end do
     end do
  end if

  if (iswitch /= 2) then
     call do_mesh(rmax,zmesh,xmin,dx,0,ndm,mesh,r,r2,sqr)
     rhoc=0.0_dp
  endif
  !
  !  In the spin polarized case adjust the occupations
  !
  if (lsd == 1)  call occ_spin(nwf,nwfx,el,nn,ll,oc,isw)

  if (iswitch == 1) then
     !
     !    no more data needed for AE calculations
     !
     return
     !
  else if (iswitch == 2) then
     !
     !    reading input for PP testing
     !
     jjs=0.0_dp
     jjts=0.0_dp
     jjtsc=0.0_dp
     
     do n=1,nwf
        oc_old(n)=oc(n)
     enddo

     nconf=1
     configts=' '
     file_wavefunctionsps= ' '
     file_pseudo=' '
     file_tests=' '

     read(5,test,err=300,iostat=ios)
300  call errore('ld1_readin','reading test',abs(ios))
     
     if (nconf > ncmax1.or.nconf < 1) &
          call errore('ld1_readin','nconf is wrong',1)
     if (iswitch == 3 .and. nconf > 1) &
          call errore('ld1_readin','too many test configurations',1)

     do nc=1,nconf
        if (configts(nc) == ' ') then
           call read_psconfig (rel, lsd, nwftsc(nc), eltsc(1,nc), &
                nntsc(1,nc), lltsc(1,nc), octsc(1,nc), iswtsc(1,nc), &
                jjtsc(1,nc), edum(1), rcuttsc(1,nc), rcutustsc(1,nc) )
        else
           call el_config(configts(nc),.false.,nwftsc(nc),eltsc(1,nc),  &
                &     nntsc(1,nc),lltsc(1,nc),octsc(1,nc),iswtsc(1,nc))
        endif
     enddo
35   call errore('ld1_readin','reading test wavefunctions',abs(ios))
     !
     !  adjust the occupations of the test cases if this is a lsd run
     !
     if (lsd == 1) then
        do nc=1,nconf
           call occ_spin(nwftsc(nc),nwfsx,eltsc(1,nc),nntsc(1,nc),lltsc(1,nc),&
                octsc(1,nc), iswtsc(1,nc)) 
        enddo
     endif
     !
     !    reading the pseudopotential
     !
     if (file_pseudo == ' ') &
          call errore('ld1_readin','file_pseudo is needed',1)
     if (matches('upf',file_pseudo) .or. matches('UPF', file_pseudo)) then
        !
        !    UPF format
        !
        call read_pseudoupf
        !
     else if ( matches('rrkj3', file_pseudo) .or. &
               matches('RRKJ3', file_pseudo)) then
        !
        !    Old RRKJ format
        !
        call read_newpseudo (ios)
        !
        if (ios /= 0) then
           !
           !    try old Norm-Conserving format
           !
           call read_pseudo  &
                (file_pseudo,zed,xmin,rmax,dx,mesh,ndm,r,r2,sqr, &
                dft,lmax,lloc,zval,nlcc,rhoc,vnl,vnlo,vpsloc,rel)
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
        call read_pseudo  &
             (file_pseudo,zed,xmin,rmax,dx,mesh,ndm,r,r2,sqr, &
             dft,lmax,lloc,zval,nlcc,rhoc,vnl,vnlo,vpsloc,rel)
        pseudotype = 1
     endif
     !
     ! initialize a few variables that may cause trouble otherwise
     !
     file_pseudopw=' '
     file_logderps=' '
     
  else if (iswitch == 3) then
     !
     !    reading input for PP generation
     !
     file_pseudopw=' '
     file_screen=' '
     file_core=' '
     file_chi=' '
     file_beta=' '
     file_qvan=' '
     file_logderps=' '
     zval=0.0_dp
     lloc=-1
     rcloc=1.5_dp
     nlcc=.false.
     rcore=0.0_dp
     rho0=0.0_dp
     tm  = .false.
     pseudotype=0

     read(5,inputp,err=500,iostat=ios)
500  call errore('ld1_readin','reading inputp',abs(ios))

     if (rcloc <=0.0_dp) &
          call errore('ld1_readin','rcloc is negative',1)
     if (pseudotype < 1.or.pseudotype > 3) &
          call errore('ld1_readin','specify correct pseudotype',1)
     !
     call read_psconfig (rel, lsd, nwfs, els, nns, lls, ocs, &
          isws, jjs, enls, rcut, rcutus )
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
             ' suppied and calculated valence charge do not match',1)
     end if
     !
     do ns=1,nwfs
        do ns1=1,ns-1
           if (lls(ns) == lls(ns1).and.pseudotype == 1) &
                call errore('ld1_readin','two wavefunctions for same l',1)
        enddo
        !
        ! flag unbound states, i.e. those with negative occupancy
        ! and those with zero occupancy and nonzero reference energy
        ! WARNING: this should be done in a cleaner way
        !
        unbound(ns) = (ocs(ns) < 0.0_dp) .or. &
             (ocs(ns) == 0.0_dp .and. enls(ns) /= 0.0_dp) 
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
     !
     ! initialize a few variables used in subsequent PP testing
     !
     file_wavefunctionsps= ' '
     file_tests=' '
     nconf=1
     ns1 = 0
     do n=1,nwfs
        if (.not.unbound(n)) then
           !
           ! copy states used in the PP generation to testing configuration
           ! Only bound states must be copied. Note that this WILL NOT WORK
           ! if bound states are not used in the generation of the PP
           !
           ns1 = ns1 + 1
           eltsc (ns1,1)= els (n)
           nntsc (ns1,1)= nns (n)
           lltsc (ns1,1)= lls (n)
           octsc (ns1,1)= ocs (n)
           iswtsc(ns1,1)= isws(n)
           jjtsc (ns1,1)= jjs (n)
        end if
     end do
     !
     ! add additional valence states for PAW reconstruction, if not used
     ! in PP generation. They must appear last in the list of all-electron
     ! states and must be empty. To be done in a cleaner way !
     !
     do n=nwf,1,-1
        found = .false.
        do ns=1,nwfs
           if ( el(n) == els(ns) ) found = .true.
        end do
        if (found) exit
        !
        ns1 = ns1 + 1
        eltsc (ns1,1)= el (n)
        lltsc (ns1,1)= ll (n)
        if (oc(n) > 0.0_dp) &
             call errore('ld1_readin','state'//el(n)//'should not be there',n)
        octsc (ns1,1)= oc (n)
        iswtsc(ns1,1)= isw(n)
        jjtsc (ns1,1)= jj (n)
        !
        nntsc (ns1,1)= lltsc (ns1,1) + 1
        do ns=1,ns1-1
           if ( lltsc(ns,1) == lltsc(ns1,1) ) &
                nntsc(ns1,1) = nntsc(ns1,1) + 1
        end do
        !
     enddo
     !
     nwftsc(1) = ns1
     !
     do n=1,nwf
        oc_old(n)=oc(n)
     enddo
     !
  endif

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
  real(kind=dp) :: oc(nwfx)
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
  real(kind=dp) :: oc(nwfx), jj(nwfx)
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
        if (lsd == 0) then
           read(5,'(a2,2i3,f6.2)',err=20,end=20,iostat=ios) &
                el(n), nn(n), ll(n), oc(n)
           isw(n)=1
20         call errore('read_config','reading orbital (lda)',abs(ios))
        else  
           read(5,'(a2,2i3,f6.2,i3)',err=21,end=21,iostat=ios) &
                el(n), nn(n), ll(n), oc(n), isw(n)
21         call errore('read_config','reading orbital (lsd)',abs(ios))
           if(isw(n) > 2 .or. isw(n) < 1) &
                call errore('read_config','spin variable wrong ',n)
        endif
     else
        read(5,'(a2,2i3,2f6.2)',err=22,end=22,iostat=ios) &
             el(n), nn(n), ll(n), oc(n), jj(n)
        isw(n)=1
        if ((abs(ll(n)+0.5_dp-jj(n)) > 1.e-3_dp) .and. &
            (abs(ll(n)-0.5_dp-jj(n)) > 1.e-3_dp) .and. abs(jj(n)) > 1.e-3_dp)  &
            call errore('read_config','jj wrong',n)
        if (oc(n) > (2.0_dp*jj(n)+1.0_dp) .and. abs(jj(n)) > 1e-3_dp) &
             call errore('read_config','occupations wrong',n)
22      call errore('read_config','reading orbital (rel)',abs(ios))
     endif
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
  real(kind=dp) :: ocs(nwfsx), jjs(nwfsx), enls(nwfsx), &
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
           read(5,'(a2,2i3,4f6.2,i3)',err=30,end=30,iostat=ios) &
                els(n), nns(n), lls(n), ocs(n), enls(n), &
                rcut(n), rcutus(n), isws(n)
           if (isws(n) > 2 .or. isws(n) < 1) &
                call errore('read_psconfig', 'spin variable wrong ',n)
           if (ocs(n) > (2.0_dp*lls(n)+1.0_dp))                 &
             call errore('read_psconfig','occupations (ls) wrong',n)
        else
           read(5,'(a2,2i3,4f6.2)',err=30,end=30,iostat=ios) &
                els(n), nns(n), lls(n), ocs(n), enls(n), &
                rcut(n), rcutus(n)
           isws(n)=1
           if (ocs(n) > 2.0_dp*(2.0_dp*lls(n)+1.0_dp))                 &
             call errore('read_psconfig','occupations (l) wrong',n)
        end if
        jjs(n)=0.0_dp
     else
        read(5,'(a2,2i3,5f6.2)',err=30,end=30,iostat=ios) &
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

