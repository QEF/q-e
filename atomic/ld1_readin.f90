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
       ios,   &          ! I/O control
       ncheck            ! used to check input

  real(kind=dp) :: &
       rdum, zdum        ! auxiliary

  character(len=80) :: config, configts(ncmax1)
  character(len=2) :: atom, label
  character, external :: capital*1, atom_name*2
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
       pseudotype,    & ! the pseudopotential type
       file_wavefunctionsps,& ! the file with pseudowfc
       file_pseudopw, & ! the file where to write pseudopotential
       file_logderps, & ! file with logarithmic derivatives
       file_pseudo,   & ! filename of the pseudopotential
       file_tests,    & ! file with transferability test
       file_recon       ! file needed for the paw reconstruction

  namelist /inputp/ &
       tm,    &    ! use Troullier-Martins instead of RRKJ
       rho0,  &    ! value of the charge at the origin
       zval,  &    ! the pseudo valence
       lloc,  &    ! l component considered as local 
       nlcc,  &    ! if true nlcc is set
       rcore, &    ! the core radius for nlcc
       rcloc, &    ! the local cut-off for pseudo
       file_screen,   & ! file with screening potential
       file_core,     & ! file with total and core charge
       file_beta,     & ! file with the beta functions
       file_chi,      & ! file with the chi functions
       file_qvan,     & ! file with the qvan functions
       file_pseudopw, & ! the file where to write pseudopotential
       file_logderps    ! file with the pseudo log derivatives
  !
  !   read the namelist input and set default values 
  !
  atom  = '  '
  zed   = 0.d0
  xmin  = -7.d0
  dx    =  0.0125d0
  rmax  =100.0d0

  beta  =  0.2d0
  tr2   = 1.0d-14
  iswitch=1

  rlderiv=4.d0
  eminld=-3.d0
  emaxld=3.d0
  nld=0
  deld=0.03

  rel = 5 
  lsd   = 0
  dft= 'LDA'
  latt  = 0
  title = ' '
  config=' '
  file_wavefunctions= ' '
  file_logderae= ' '

  read(5,input,err=100,iostat=ios) 
100 call errore('ld1_readin','reading input namelist ',abs(ios))

  call which_dft(dft)

  if (zed == 0.d0 .and. atom /= ' ') then
     zed = dble(atomic_number(atom))
  else if (zed /= 0.d0 .and. atom == ' ') then
     if(dble(int(zed)) /= zed .or. zed < 1.d0 .or. zed > 100) &
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
  if (deld < 0.d0) &
       call errore('ld1_readin','negative deld',1)
  if (nld > nwfsx) &
       call errore('ld1_readin','too many nld',1)
  if (xmin > -2) call errore('ld1_readin','wrong xmin',1)
  if (dx <=0.d0) call errore('ld1_readin','wrong dx',1)

  if (isic == 1 .and. latt == 1) call errore('ld1_readin', &
       &    'isic and latter correction not allowed',1)

  zmesh=zed
  if (rel == 5 ) then
     if (zed >= 19.d0) then
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
     read(5,*,err=200,iostat=ios) nwf
200  call errore('ld1_readin','reading nwf ',abs(ios))
     if (nwf <= 0) call errore('ld1_readin','nwf is wrong',1)
     if (nwf > nwfx) call errore('ld1_readin','too many wfcs',1)
     !
     !     read the occupation of the states
     !
     do n=1,nwf  
        if (rel < 2) then
           if (lsd == 0) then
              read(5,'(a2,2i3,f6.2)',err=20,end=20,iostat=ios) &
                   el(n), nn(n), ll(n), oc(n)
              isw(n)=1
20            call errore('ld1','reading wavefunctions',abs(ios))
           else  
              read(5,'(a2,2i3,f6.2,i3)',err=21,end=21,iostat=ios) &
                   el(n), nn(n), ll(n), oc(n), isw(n)
21            call errore('ld1','reading wavefunctions 1',abs(ios))
              if(isw(n) > 2 .or. isw(n) < 1) &
                   call errore('ld1_readin','spin variable wrong ',n)
           endif
        else
           read(5,'(a2,2i3,2f6.2)',err=22,end=22,iostat=ios) &
                el(n), nn(n), ll(n), oc(n), jj(n)
           isw(n)=1
           if ((abs(ll(n)+0.5d0-jj(n)) > 1.d-3) .and. &
               (abs(ll(n)-0.5d0-jj(n)) > 1.d-3) .and. abs(jj(n)) > 1.d-3)  &
                call errore('ld1_readin','jj wrong',n)
           if (oc(n) > (2.d0*jj(n)+1.d0) .and. abs(jj(n)) > 1d-3) &
                call errore('ld1_readin','occupations wrong',n)
22         call errore('ld1','reading wavefunctions 2',abs(ios))
        endif
        write(label,'(a2)') el(n)
        read (label,'(i1)') ncheck
        if (ncheck /= nn(n)  .or. &
             capital(label(2:2)) == 'S' .and. ll(n) /= 0 .or. &
             capital(label(2:2)) == 'P' .and. ll(n) /= 1 .or. &
             capital(label(2:2)) == 'D' .and. ll(n) /= 2 .or. &
             capital(label(2:2)) == 'F' .and. ll(n) /= 3 .or. &
             oc(n) > 2.d0*(2*ll(n)+1) .or. nn(n) < ll(n)+1  ) &
             call errore('ld1_readin',label//' wrong?',n)
     enddo
  else
     call el_config(config,.true.,nwf,el,nn,ll,oc,isw)
  endif

  if (iswitch /= 2) then
     call do_mesh(rmax,zmesh,xmin,dx,0,ndm,mesh,r,r2,sqr)
     rhoc=0.d0
  endif
  !
  !  In the spin polarized case adjust the occupations
  !
  if (lsd == 1)  call occ_spin(nwf,nwfx,el,nn,ll,oc,isw)

  if (iswitch == 1) return

  jjs=0.d0
  jjts=0.d0
  jjtsc=0.d0

  do n=1,nwf
     oc_old(n)=oc(n)
  enddo

  nconf=1
  pseudotype=0
  configts=' '
  file_wavefunctionsps= ' '
  file_logderps=' '
  file_pseudo=' '
  file_pseudopw=' '
  file_tests=' '
  file_recon= ' '

  read(5,test,err=300,iostat=ios)
300 call errore('ld1_readin','reading test',abs(ios))

  if (pseudotype < 1.or.pseudotype > 3) &
       call errore('ld1_readin','wrong pseudotype',1)
  if (nconf > ncmax1.or.nconf < 1) &
       call errore('ld1_readin','nconf is wrong',1)
  if (nconf > 1.and.file_pseudopw /= ' ') &
       call errore('ld1_readin','test or write pseudo?',1)
  if (file_pseudopw == file_pseudo .and. file_pseudo /= ' ') &
       call errore('ld1_readin','rewrite file pseudo?',1)
  if (iswitch == 3 .and. nconf > 1) &
       call errore('ld1_readin','too many test configurations',1)

  do nc=1,nconf
     if (configts(nc) == ' ') then
        read(5,*,err=400,iostat=ios) nwftsc(nc)
400     call errore('ld1_readin','reading nwfts',abs(ios))

        if (nwftsc(nc) <=0.or.nwftsc(nc) > nwfsx) &
             call errore('ld1_readin','nwfts is wrong',1)

        do n=1,nwftsc(nc)
           if (rel < 2) then
              if (lsd == 1) then
                 read(5,'(a2,2i3,4f6.2,i3)',err=23,end=23,iostat=ios)&
                      eltsc(n,nc), nntsc(n,nc), lltsc(n,nc), &
                      octsc(n,nc), rdum, rcuttsc(n,nc),  &
                      rcutustsc(n,nc), iswtsc(n,nc)

                 if (iswtsc(n,nc) > 2.or.iswtsc(n,nc) < 1) &
                      call errore('ld1_readin', &
                      'spin ps variable wrong ',n)
23               call errore('ld1_readin','reading wavefunctions ',abs(ios))

              else
                 read(5,'(a2,2i3,4f6.2)',err=24,end=24,iostat=ios) &
                      eltsc(n,nc), nntsc(n,nc), lltsc(n,nc), &
                      octsc(n,nc), rdum, rcuttsc(n,nc), rcutustsc(n,nc)
                 iswtsc(n,nc)=1
24               call errore('ld1_readin','reading wavefunctions ',abs(ios))
              endif
           else
              read(5,'(a2,2i3,5f6.2)',err=25,end=25,iostat=ios) &
                   eltsc(n,nc), nntsc(n,nc), lltsc(n,nc), &
                   octsc(n,nc), rdum, rcuttsc(n,nc), &
                   rcutustsc(n,nc), jjtsc(n,nc)
25            call errore('ld1_readin','reading wavefunctions ',abs(ios))
              iswtsc(n,nc)=1
              if ((abs(lltsc(n,nc)+0.5d0-jjtsc(n,nc)) > 1.d-3) .and.   &
                  (abs(lltsc(n,nc)-0.5d0-jjtsc(n,nc)) > 1.d-3) .and.   &
                   abs(jjtsc(n,nc)) > 1.d-3) call errore('ld1_readin', &
                   'jjtsc wrong',n)
              if (octsc(n,nc) > (2.d0*jjtsc(n,nc)+1.d0)      &
                   .and.abs(jjtsc(n,nc)) > 1d-3) &
                   call errore('ld1_readin','test occupations wrong',n)
           endif
           do i=1,n-1
              if(eltsc(n,nc) == eltsc(i,nc) .and. &
                 iswtsc(n,nc) == iswtsc(i,nc).and.jjtsc(n,nc) == jjtsc(i,nc)) &
                   call errore('ld1_readin','why many equal states?',-n)
           end do

           write(label,'(a2)') eltsc(n,nc)
           if (capital(label(2:2)) == 'S' .and. lltsc(n,nc) /= 0 .or. &
               capital(label(2:2)) == 'P' .and. lltsc(n,nc) /= 1 .or. &
               capital(label(2:2)) == 'D' .and. lltsc(n,nc) /= 2 .or. &
               capital(label(2:2)) == 'F' .and. lltsc(n,nc) /= 3 .or. &
               octsc(n,nc) > 2.d0**(2*lltsc(n,nc)+1) .or.   &
               nntsc(n,nc) < lltsc(n,nc)+1 )  &
               call errore('ld1_readin','test-label'//' wrong?',n)
        enddo
     else
        call el_config(configts(nc),.false.,nwftsc(nc),eltsc(1,nc),  &
             &     nntsc(1,nc),lltsc(1,nc),octsc(1,nc),iswtsc(1,nc))
     endif
  enddo
35 call errore('ld1_readin','reading test wavefunctions',abs(ios))
  !
  !  adjust the occupations of the test cases if this is a lsd run
  !
  if (lsd == 1) then
     do nc=1,nconf
        call occ_spin(nwftsc(nc),nwfsx,eltsc(1,nc),nntsc(1,nc),lltsc(1,nc), &
             &   octsc(1,nc), iswtsc(1,nc)) 
     enddo
  endif
  !
  !    reading the pseudopotential
  !
  if (iswitch == 2) then
     jjs=0.d0
     if (file_pseudo == ' ') &
          call errore('ld1_readin','iswitch=2 and file_pseudo?',1)
     if (matches('upf',file_pseudo) .or. matches('UPF', file_pseudo)) then
        call read_pseudoupf
     else
        if (pseudotype == 1) then
           call read_pseudo &
                (file_pseudo,zed,xmin,rmax,dx,mesh,ndm,r,r2,sqr, &
                dft,lmax,lloc,zval,nlcc,rhoc,vnl,vnlo,vpsloc,rel)
        elseif (pseudotype == 2 .or. pseudotype == 3) then
           call read_newpseudo
           lmax=0
           do ns=1,nwfs
              lmax=max(lmax,lls(ns))
           enddo
        else
           call errore('ld1_readin','pseudotype not programmed ',1)
        endif
     endif
  endif

  if (iswitch == 3) then
     !
     !    pseudopotential input reading
     !
     file_pseudopw=' '
     file_screen=' '
     file_core=' '
     file_chi=' '
     file_beta=' '
     file_qvan=' '
     zval=0.d0
     lloc=-1
     rcloc=1.5d0
     nlcc=.false.
     rcore=0.d0
     rho0=0.d0
     tm  = .false.

     read(5,inputp,err=500,iostat=ios)
500  call errore('ld1_readin','reading inputp',abs(ios))

     if (file_pseudopw == ' ') &
          call errore('ld1_readin','iswitch=3 and file_pseudopw?',1)
     if (rcloc <=0.d0) &
          call errore('ld1_readin','rcloc is negative',1)
     if (nlcc .and. rcore == 0.d0) &
          call errore('ld1_readin','must give a core radius',1)

     read(5,*,err=600,iostat=ios) nwfs
600  call errore('ld1_readin','reading nwfs',abs(ios))

     if (nwfs <= 0 .or. nwfs > nwfsx) &
          call errore('ld1_readin','nwfs is wrong',1)

     do n=1,nwfs
        if (rel < 2) then
           read(5,'(a2,2i3,4f6.2)',err=30,end=30,iostat=ios) &
                els(n), nns(n), lls(n), ocs(n), enls(n), &
                rcut(n), rcutus(n)
           jjs(n)=0.d0
        else
           read(5,'(a2,2i3,5f6.2)',err=30,end=30,iostat=ios) &
                els(n), nns(n), lls(n), ocs(n), enls(n),       &
                rcut(n), rcutus(n), jjs(n)
           if ((abs(lls(n)+0.5d0-jjs(n)) > 1.d-3).and.      &
                (abs(lls(n)-0.5d0-jjs(n)) > 1.d-3)  )        &
                call errore('ld1_readin', 'jjs wrong',n)
           if (octsc(n,nc) > (2.d0*jjtsc(n,nc)+1.d0))       &
                call errore('ld1_readin','test occupations wrong',n)
        endif
        write(label,'(a2)') els(n)
        if (capital(label(2:2)) == 'S'.and.lls(n) /= 0.or.   &
             capital(label(2:2)) == 'P'.and.lls(n) /= 1.or.   &
             capital(label(2:2)) == 'D'.and.lls(n) /= 2.or.   &
             capital(label(2:2)) == 'F'.and.lls(n) /= 3.or.   &
             ocs(n) > 2*(2*lls(n)+1).or.                 &
             nns(n) < lls(n)+1 )                         &
             call errore('ld1_readin','ps-label'//' wrong?',n)
        if (rcut(n) > rcutus(n)) &
             call errore('ld1_readin','rcut or rcutus is wrong',1)
     enddo
30   call errore('ld1_readin','reading pseudo wavefunctions',abs(ios))

     lmax=0
     do ns=1,nwfs
        do ns1=1,ns-1
           if (lls(ns) == lls(ns1).and.pseudotype == 1) &
                call errore('ld1_readin','too wavefunction for l',1)
        enddo
        lmax=max(lmax,lls(ns))
     enddo
     nlc=0
     nnl=0
  endif

  return
end subroutine ld1_readin

subroutine occ_spin(nwf,nwfx,el,nn,ll,oc,isw)
  !
  !  This routine splits the occupations of the states between spin-up
  !  and spin down. If the occupations are lower than 2*l+1 it does
  !  nothing, otherwise 2*l+1 states are assumed with spin up and
  !  the difference with spin down. 
  !
  implicit none
  integer, parameter :: dp=kind(1.d0)
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
           oc(nwf)=0.d0
           if (isw(n) == 1) isw(nwf)=2 
           if (isw(n) == 2) isw(nwf)=1 
        endif
     endif
  enddo
  return
end subroutine occ_spin
