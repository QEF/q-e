!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------

subroutine vcsmd
  !-----------------------------------------------------------------------
  ! Main (interface) routine between PWSCF and the variable-cell shape
  ! molecular dynamics code by R.M. Wentzcovitch, PRB 44, 2358 (1991).
  !
  ! Molecular and/or cell dynamics is performed according to the value of
  ! the switch variable calc:
  !
  !  calc  = 'md'   : standard molecular dynamics
  !  calc  = 'mm'   : structural minimization by damped dynamics
  !  calc  = 'cd'   : Parrinello-Rahman cell dynamics
  !  calc  = 'cm'   : Parrinello-Rahman cell minimization by damped dynami
  !  calc  = 'nd'   : Wentzcovitch's new cell dynamics
  !  calc  = 'nm'   : Wentzcovitch's new cell minimization by damped dynam
  !
  ! Dynamics performed using Beeman algorithm, J. Comp. Phys. 20, 130 (1976))

  !
#include "machine.h"
  use pwcom
  use io, only : prefix
#ifdef __PARA
  use para
#endif
  implicit none
  !
  ! I/O variable first
  !
  ! PWSCF variables
  !  nat  = total number of atoms
  !  ntyp = total number of atomic types
  !  ityp(na)  = atomic type for na-th atom
  !  tau(i,na) =  position of the na-th atom
  !  at (icar,ivec) = direct Bravais lattice vectors
  !  bg (icar,ivec) = reciprocal lattice vectors
  !  amass(nt) = mass (in atomic ryd units) for atom of nt-th type
  !  cmass = cell mass in ryd units.
  !  press = target pressure in ryd/(a.u.)^3
  !
  ! local variable
  !

  real(kind=DP) :: p,            & ! virial pressure
                   vcell,        & ! cell volume
                   avec (3, 3),  & ! at(3,3) * alat
                   aveci (3, 3), & ! avec at t-dt
                   avecd (3, 3), & ! d(avec)/dt
                   avec2d (3, 3),& ! d2(avec)/dt2
                   avec2di (3, 3),&! d2(avec)/dt2 at t-dt
                   avec0 (3, 3), & ! avec at t = 0
                   sig0 (3, 3),  & ! sigma at t=0
                   v0              ! volume at t=0
  real(kind=DP), allocatable ::      &
                   rat (:, :),   & ! atomic positions (lattice coord)
                   rati (:, :),  & ! rat at previous step
                   ratd (:, :),  & ! rat derivatives at current step
                   rat2d (:, :), & ! rat 2nd derivatives at current step
                   rat2di (:,:), & ! rat 2nd derivatives at previous step
                   tauold (:, :, :)! additional  history variables
  real(kind=DP) :: &
           avmod (3), theta (3, 3), & ! used to monitor cell dynamics
           enew, e_start,           & ! DFT energy at current and first step
           eold,                    & ! DFT energy at previous step
           uta, eka, eta, ekla, utl, etl, ut, ekint, edyn,  & ! other energies
           acu, ack, acp, acpv, avu, avk, avp, avpv,        & ! acc.& avrg. ener
           tnew, pv,                & ! istantaneous temperature and p*vcell
           sigmamet (3, 3),         & ! sigma = avec^-1 * vcell = bg/alat*omega
           vx2 (ntypx), vy2 (ntypx), vz2 (ntypx),     & ! work vectors
           vmean (ntypx), rms (ntypx), ekin (ntypx),  & ! work vectors
           tempo, time_au, epsp

  character(len=3) :: iost  ! status (old or new) for I/O files
  character(len=6) :: ipos  ! status ('append' or 'asis') for I/O files

  logical :: exst
  integer :: idone, na, nst, ipol, i, j, k
  ! last completed time step
  ! counter on atoms
  ! counter on completed moves
  !
  ! I/O units
  !
  integer :: iun_e, iun_eal, iun_ave, iun_p, iun_avec, iun_tv

  parameter (iun_e=21, iun_eal=22, iun_ave=23, iun_p=24, iun_avec=25, iun_tv=26)
  !
  ! Allocate work arrays
  !
  allocate (rat(3,nat), rati(3,nat), ratd(3,nat), rat2d(3,nat), rat2di(3,nat), &
            tauold(3,nat,3))
  !
  ! open MD history file (if not present this is a new run!)
  !
  call seqopn (4, trim(prefix)//'.md', 'formatted', exst)
  if (.not.exst) then
     close (unit = 4, status = 'delete')
     if (istep.ne.1) call errore ('vcsmd', ' previous MD history got lost', 1)
     tnew = 0.d0
     acu = 0.d0
     ack = 0.d0
     acp = 0.d0
     acpv = 0.d0
     avu = 0.d0
     avk = 0.d0
     avp = 0.d0
     avpv = 0.d0
     nzero = 0
     tauold(:,:,:) = 0.d0

     eold = etot + 2*epse ! set value for eold at first iteration
  else
     !
     ! read MD run history (idone is the last completed MD step)
     !
     read (4, * ) rati, ratd, rat2d, rat2di, tauold
     read (4, * ) aveci, avecd, avec2d, avec2di
     read (4, * ) avec0, sig0, v0, e_start, eold
     read (4, * ) acu, ack, acp, acpv, avu, avk, avp, avpv, sigmamet
     read (4, * ) idone, nzero, ntimes
     close (unit = 4, status = 'keep')
     istep = idone+1
     call DCOPY (3 * nat, tauold (1, 1, 2), 1, tauold (1, 1, 3), 1)
     call DCOPY (3 * nat, tauold (1, 1, 1), 1, tauold (1, 1, 2), 1)
  endif

  !
  ! check if convergence for structural minimization is acheived
  !
  if (calc.eq.'mm') then
     conv_ions = eold - etot .lt. epse
     do i = 1, 3
        do na = 1, nat - fixatom
           conv_ions = conv_ions.and.abs (force (i, na) ) .lt. epsf
        enddo
     enddo
     if (conv_ions) then
        write (6,'(/5x,"Damped Dynamics: convergence achieved, Efinal=",&
              &     f15.8)') etot
        write (6,'(/72("-")//5x,"Final estimate of positions")')
        if (ltaucry) write (6, '(/5x,"Cartesian coordinates")')
        do na = 1, nat
           write (6,'(a3,3x,3f14.9)') atm(ityp(na)), (tau(i,na), i=1,3)
        enddo
        if (ltaucry) then
           write (6, '(/5x,"In crystal coordinates")')
           call cryst_to_cart (nat, tau, bg, - 1)
           do na = 1, nat
              write (6,'(a3,3x,3f14.9)') atm(ityp(na)), (tau(i,na), i=1,3)
           enddo
           call cryst_to_cart (nat, tau, at, 1)
        endif
        write (6, '(/)')
        return
     end if
  end if
  if (calc.eq.'nm' .or. calc.eq.'cm') then
     epsp = 0.5  ! kbar
     conv_ions = eold - etot .lt. epse
     do i = 1, 3
        do na = 1, nat - fixatom
           conv_ions = conv_ions.and.abs (force (i, na) ) .lt. epsf
        enddo
     enddo
     do i =1,3
        conv_ions = conv_ions .and. abs (sigma(i,i) - press)*uakbar.lt. epsp
        do j =i+1,3
           conv_ions = conv_ions .and. abs ( sigma(i,j) )*uakbar.lt. epsp
        end do
     end do
     if (conv_ions) then
        if (calc.eq.'cm') write (6,'(/5x,"Parrinello-Rahman Damped Dynamics: convergence achieved, Efinal=",&
              &     f15.8)') etot
        if (calc.eq.'nm') write (6,'(/5x,"Wentzcovitch Damped Dynamics: convergence achieved, Efinal=",&
              &     f15.8)') etot
        write (6,'(/72("-")//5x,"Final estimate of lattice vectors (input alat units)")')
        write (6, '(3f14.9)') ( (at (i, k) , i = 1, 3) , k = 1, 3)
        write (6,'(a,f12.4,a)') '  final unit-cell volume =', omega, ' (a.u.)^3'
        write (6,'(a,f12.4,a)') '  input alat = ', alat, ' (a.u.)'

        write (6,'(//5x,"Final estimate of positions")')
        if (ltaucry) write (6, '(/5x,"Cartesian coordinates (input alat units)")')
        do na = 1, nat
           write (6,'(a3,3x,3f14.9)') atm(ityp(na)), (tau(i,na), i=1,3)
        enddo
        if (ltaucry) then
           write (6, '(/5x,"In crystal coordinates")')
           call cryst_to_cart (nat, tau, bg, - 1)
           do na = 1, nat
              write (6,'(a3,3x,3f14.9)') atm(ityp(na)), (tau(i,na), i=1,3)
           enddo
           call cryst_to_cart (nat, tau, at, 1)
        endif
        write (6, '(/)')
        return
     end if
  end if

  call DCOPY (3 * nat, tau, 1, tauold (1, 1, 1), 1)
  time_au = 0.0000242 * e2
  tempo = (istep - 1) * dt * time_au

  if (istep.eq.1 .and. calc.eq.'mm')  &
     write(6,'(/5x,"Damped Dynamics Minimization", /5x, &
             & "convergence thresholds: EPSE = ", e8.2,"  EPSF = ",e8.2)') &
               epse, epsf
  if (istep.eq.1 .and. calc.eq.'cm')  &
     write(6,'(/5x,"Parrinello-Rahman Damped Cell-Dynamics Minimization", /5x, &
             & "convergence thresholds: EPSE = ", e8.2,"  EPSF = ",e8.2,&
             & "  EPSP = ",e8.2 )') epse, epsf, epsp
  if (istep.eq.1 .and. calc.eq.'nm')  &
     write(6,'(/5x,"Wentzcovitch Damped Cell-Dynamics Minimization", /5x, &
             & "convergence thresholds: EPSE = ", e8.2,"  EPSF = ",e8.2,&
             & "  EPSP = ",e8.2 )') epse, epsf, epsp

  write (6, '(/5x,"Entering Dynamics;  it = ",i5,"   time = ", &
       &                          f8.5," pico-seconds"/)') istep, tempo
  !       write (*,*) ' enter vcsmd ', istep,idone,exst
  !
  ! save cell shape of previous step
  !
  call DCOPY (9, at, 1, at_old, 1)
  omega_old = omega
  !
  ! Translate
  !
  ! define rat as the atomic positions in lattice coordinates
  !
  call DCOPY (3 * nat, tau, 1, rat, 1)
  call cryst_to_cart (nat, rat, bg, - 1)
  !
  call DCOPY (9, at, 1, avec, 1)
  call DSCAL (9, alat, avec, 1)
  !
  ! convert forces to lattice coordinates
  !
  call cryst_to_cart (nat, force, bg, - 1)
  call DSCAL (3 * nat, 1.d0 / alat, force, 1)
  !
  ! scale stress to stress*omega
  !
  call DSCAL (9, omega, sigma, 1)

  vcell = omega
  if (istep.eq.1) then
     e_start = etot

     enew = etot - e_start

#ifdef DEBUG_VCSMD
     write (45,*) 'istep=',istep
     write (45,*) 'ntyp=', ntyp
     write (45,*) 'nat=', nat
     write (45,*) 'rat=',rat
     write (45,*) 'ityp=',ityp
     write (45,*) 'avec=',avec
     write (45,*) 'vcell=',vcell
     write (45,*) 'force=',force
     write (45,*) 'sigma=',sigma
     write (45,*) 'calc=',calc
     write (45,*) 'temperature=',temperature
     write (45,*) 'vx2=',vx2
     write (45,*) 'vy2=',vy2
     write (45,*) 'vz2=',vz2
     write (45,*) 'rms=',rms
     write (45,*) 'vmean=',vmean
     write (45,*) 'ekin=',ekin
     write (45,*) 'avmod=',avmod
     write (45,*) 'theta=',theta
     write (45,*) 'amass=',amass
     write (45,*) 'cmass=',cmass
     write (45,*) 'press=',press
     write (45,*) 'p=',p
     write (45,*) 'dt=',dt
     write (45,*) 'aveci=',aveci
     write (45,*) 'avecd=',avecd
     write (45,*) 'avec2d=',avec2d
     write (45,*) 'avec2di=',avec2di
     write (45,*) 'sigmamet=',sigmamet
     write (45,*) 'sig0=', sig0
     write (45,*) 'avec0=',avec0
     write (45,*) 'v0=',v0
     write (45,*) 'rati=',rati
     write (45,*) 'ratd=',ratd
     write (45,*) 'rat2d=',rat2d
     write (45,*) 'rat2di=',rat2di
     write (45,*) 'enew=',enew
     write (45,*) 'uta=',uta
     write (45,*) 'eka=',eka
     write (45,*) 'eta=',eta
     write (45,*) 'ekla=',ekla
     write (45,*) 'utl=',utl
     write (45,*) 'etl=',etl
     write (45,*) 'ut=',ut
     write (45,*) 'ekint=',ekint
     write (45,*) 'edyn=',edyn
#endif

     call init (ntyp, nat, ntyp, nat-fixatom, rat, ityp, avec, vcell, force, &
          sigma, calc, temperature, vx2, vy2, vz2, rms, vmean, ekin, &
          avmod, theta, amass, cmass, press, p, dt, aveci, avecd, avec2d, &
          avec2di, sigmamet, sig0, avec0, v0, rati, ratd, rat2d, rat2di, &
          enew, uta, eka, eta, ekla, utl, etl, ut, ekint, edyn)

#ifdef DEBUG_VCSMD
     write (46,*) 'istep=',istep
     write (46,*) 'ntyp=', ntyp
     write (46,*) 'nat=', nat
     write (46,*) 'rat=',rat
     write (46,*) 'ityp=',ityp
     write (46,*) 'avec=',avec
     write (46,*) 'vcell=',vcell
     write (46,*) 'force=',force
     write (46,*) 'sigma=',sigma
     write (46,*) 'calc=',calc
     write (46,*) 'temperature=',temperature
     write (46,*) 'vx2=',vx2
     write (46,*) 'vy2=',vy2
     write (46,*) 'vz2=',vz2
     write (46,*) 'rms=',rms
     write (46,*) 'vmean=',vmean
     write (46,*) 'ekin=',ekin
     write (46,*) 'avmod=',avmod
     write (46,*) 'theta=',theta
     write (46,*) 'amass=',amass
     write (46,*) 'cmass=',cmass
     write (46,*) 'press=',press
     write (46,*) 'p=',p
     write (46,*) 'dt=',dt
     write (46,*) 'aveci=',aveci
     write (46,*) 'avecd=',avecd
     write (46,*) 'avec2d=',avec2d
     write (46,*) 'avec2di=',avec2di
     write (46,*) 'sigmamet=',sigmamet
     write (46,*) 'sig0=', sig0
     write (46,*) 'avec0=',avec0
     write (46,*) 'v0=',v0
     write (46,*) 'rati=',rati
     write (46,*) 'ratd=',ratd
     write (46,*) 'rat2d=',rat2d
     write (46,*) 'rat2di=',rat2di
     write (46,*) 'enew=',enew
     write (46,*) 'uta=',uta
     write (46,*) 'eka=',eka
     write (46,*) 'eta=',eta
     write (46,*) 'ekla=',ekla
     write (46,*) 'utl=',utl
     write (46,*) 'etl=',etl
     write (46,*) 'ut=',ut
     write (46,*) 'ekint=',ekint
     write (46,*) 'edyn=',edyn
#endif
  else
     !         write (*,'(3f12.6)') ((ratd(ipol,na),ipol=1,3),na=1,nat)

     enew = etot - e_start

#ifdef DEBUG_VCSMD
     write (45,*) 'xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx'
     write (45,*) 'istep=', istep
     write (45,*) 'ntyp=',ntyp
     write (45,*) 'nat=',nat
     write (45,*) 'ityp=',ityp
     write (45,*) 'rat=',rat
     write (45,*) 'avec=',avec
     write (45,*) 'vcell=',vcell
     write (45,*) 'force=',force
     write (45,*) 'sigma=',sigma
     write (45,*) 'calc=',calc
     write (45,*) 'avmod=',avmod
     write (45,*) 'theta=',theta
     write (45,*) 'amass=',amass
     write (45,*) 'cmass=',cmass
     write (45,*) 'press=',press
     write (45,*) 'p=',p
     write (45,*) 'dt=',dt
     write (45,*) 'avecd=',avecd
     write (45,*) 'avec2d=',avec2d
     write (45,*) 'aveci=',aveci
     write (45,*) 'avec2di=',avec2di
     write (45,*) 'sigmamet=',sigmamet
     write (45,*) 'sig0=',sig0
     write (45,*) 'avec0=',avec0
     write (45,*) 'v0=',v0
     write (45,*) 'ratd=',ratd
     write (45,*) 'rat2d=',rat2d
     write (45,*) 'rati=',rati
     write (45,*) 'rat2di=',rat2di
     write (45,*) 'enew=',enew
     write (45,*) 'uta=',uta
     write (45,*) 'eka=',eka
     write (45,*) 'eta=',eta
     write (45,*) 'ekla=',ekla
     write (45,*) 'utl=',utl
     write (45,*) 'etl=',etl
     write (45,*) 'ut=',ut
     write (45,*) 'ekint=',ekint
     write (45,*) 'edyn=',edyn
     write (45,*) 'temperature=',temperature
     write (45,*) 'ttol=',ttol
     write (45,*) 'ntcheck=',ntcheck
     write (45,*) 'ntimes=',ntimes
     write (45,*) 'istep=',istep
     write (45,*) 'tnew=',tnew
     write (45,*) 'nzero=',nzero
     write (45,*) 'nat=',nat
     write (45,*) 'acu=',acu
     write (45,*) 'ack=',ack
     write (45,*) 'acp=',acp
     write (45,*) 'acpv=',acpv
     write (45,*) 'avu=',avu
     write (45,*) 'avk=',avk
     write (45,*) 'avp=',avp
     write (45,*) 'avpv=',avpv
#endif

     call move (ntyp, nat, ntyp, ityp, rat, avec, vcell, force, &
          sigma, calc, avmod, theta, amass, cmass, press, p, dt, avecd, &
          avec2d, aveci, avec2di, sigmamet, sig0, avec0, v0, ratd, rat2d, &
          rati, rat2di, enew, uta, eka, eta, ekla, utl, etl, ut, ekint, &
          edyn, temperature, ttol, ntcheck, ntimes, istep, tnew, nzero, &
          nat-fixatom, acu, ack, acp, acpv, avu, avk, avp, avpv)

#ifdef DEBUG_VCSMD
     write (46,*) 'xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx'
     write (46,*) 'istep=', istep
     write (46,*) 'ntyp=',ntyp
     write (46,*) 'nat=',nat
     write (46,*) 'ityp=',ityp
     write (46,*) 'rat=',rat
     write (46,*) 'avec=',avec
     write (46,*) 'vcell=',vcell
     write (46,*) 'force=',force
     write (46,*) 'sigma=',sigma
     write (46,*) 'calc=',calc
     write (46,*) 'avmod=',avmod
     write (46,*) 'theta=',theta
     write (46,*) 'amass=',amass
     write (46,*) 'cmass=',cmass
     write (46,*) 'press=',press
     write (46,*) 'p=',p
     write (46,*) 'dt=',dt
     write (46,*) 'avecd=',avecd
     write (46,*) 'avec2d=',avec2d
     write (46,*) 'aveci=',aveci
     write (46,*) 'avec2di=',avec2di
     write (46,*) 'sigmamet=',sigmamet
     write (46,*) 'sig0=',sig0
     write (46,*) 'avec0=',avec0
     write (46,*) 'v0=',v0
     write (46,*) 'ratd=',ratd
     write (46,*) 'rat2d=',rat2d
     write (46,*) 'rati=',rati
     write (46,*) 'rat2di=',rat2di
     write (46,*) 'enew=',enew
     write (46,*) 'uta=',uta
     write (46,*) 'eka=',eka
     write (46,*) 'eta=',eta
     write (46,*) 'ekla=',ekla
     write (46,*) 'utl=',utl
     write (46,*) 'etl=',etl
     write (46,*) 'ut=',ut
     write (46,*) 'ekint=',ekint
     write (46,*) 'edyn=',edyn
     write (46,*) 'temperature=',temperature
     write (46,*) 'ttol=',ttol
     write (46,*) 'ntcheck=',ntcheck
     write (46,*) 'ntimes=',ntimes
     write (46,*) 'istep=',istep
     write (46,*) 'tnew=',tnew
     write (46,*) 'nzero=',nzero
     write (46,*) 'nat=',nat
     write (46,*) 'acu=',acu
     write (46,*) 'ack=',ack
     write (46,*) 'acp=',acp
     write (46,*) 'acpv=',acpv
     write (46,*) 'avu=',avu
     write (46,*) 'avk=',avk
     write (46,*) 'avp=',avp
     write (46,*) 'avpv=',avpv
#endif

  endif

  pv = p * omega
  !
  !  write on output files several control quantities
  !
#ifdef __PARA
  !
  ! only the first processor does this I/O
  !
  if (me.eq.1.and.mypool.eq.1) then
#endif
     !
     !  NB: at the first iteration files should not be present,
     !      for subsequent iterations they should.
     !
     if (istep.eq.1) then
        call delete_if_present ( 'e')
        call delete_if_present ( 'eal')
        call delete_if_present ( 'ave')
        call delete_if_present ( 'p')
        call delete_if_present ( 'avec')
        call delete_if_present ( 'tv')
        iost = 'new'
        ipos = 'asis'
     else
        iost = 'old'
        ipos = 'append'
     endif
     open (unit=iun_e,   file='e',   status=iost,form='formatted',position=ipos)
     open (unit=iun_eal, file='eal', status=iost,form='formatted',position=ipos)
     open (unit=iun_ave, file='ave', status=iost,form='formatted',position=ipos)
     open (unit=iun_p,   file='p',   status=iost,form='formatted',position=ipos)
     open (unit=iun_avec,file='avec',status=iost,form='formatted',position=ipos)
     open (unit=iun_tv,  file='tv',  status=iost,form='formatted',position=ipos)
     nst = istep - 1
     write (iun_e,   101) ut, ekint, edyn, pv, nst
     write (iun_eal, 103) uta, eka, eta, utl, ekla, etl, nst
     write (iun_ave, 104) avu, avk, nst
     write (iun_p,   105) press, p, avp, nst
     if (calc (1:1) .ne.'m') write (iun_avec, 103) (avmod(k), k=1,3),  &
                                    theta(1,2), theta(2,3), theta(3,1), nst

     write (iun_tv, 104) vcell, tnew, nst
     close (unit = iun_e,    status = 'keep')
     close (unit = iun_eal,  status = 'keep')
     close (unit = iun_ave,  status = 'keep')
     close (unit = iun_p,    status = 'keep')
     close (unit = iun_avec, status = 'keep')
     close (unit = iun_tv,   status = 'keep')
#ifdef __PARA
  endif
#endif
  !
  ! update configuration in PWSCF variables
  !
  call DCOPY (9, avec, 1, at, 1)
  call DSCAL (9, 1.d0 / alat, at, 1)
  call volume (alat, at (1, 1), at (1, 2), at (1, 3), omega)

  call recips (at(1,1), at(1,2), at(1,3), bg(1,1), bg(1,2) , bg(1,3) )
  call DCOPY (3 * nat, rat, 1, tau, 1)
  if (lmovecell) then
     write (6, * ) ' new lattice vectors (alat unit) :'
     write (6, '(3f14.9)') ( (at (i, k) , i = 1, 3) , k = 1, 3)
     write (6,'(a,f12.4,a)') '  new unit-cell volume =', omega, ' (a.u.)^3'
  endif
  write (6, * ) ' new positions in cryst coord'
  write (6,'(a3,3x,3f14.9)') (atm(ityp(na)), (tau(i,na), i=1,3),na=1,nat)
  write (6, * ) ' new positions in cart coord (alat unit)'
  call cryst_to_cart (nat, tau, at, + 1)
  write (6,'(a3,3x,3f14.9)') (atm(ityp(na)), (tau(i,na), i=1,3),na=1,nat)
  write (6, '(/5x,"Ekin = ",f14.8," Ryd   T = ",f6.1," K ", &
       &       " Etot = ",f14.8)') ekint, tnew, edyn + e_start
  !
  ! save MD history on file
  !
  call seqopn (4, trim(prefix)//'.md', 'formatted', exst)
  write (4, * ) rati, ratd, rat2d, rat2di, tauold
  write (4, * ) aveci, avecd, avec2d, avec2di
  write (4, * ) avec0, sig0, v0, e_start, etot
  write (4, * ) acu, ack, acp, acpv, avu, avk, avp, avpv, sigmamet
  write (4, * ) istep, nzero, ntimes
  close (unit = 4, status = 'keep')
  !
  ! find alpha & beta for charge-density and wfc extrapolation (uptate_pot
  !
  if (calc(1:1).eq.'m') call find_alpha_and_beta(nat,tau,tauold,alpha0,beta0)
  !
  ! Deallocate
  !
  deallocate (rat, rati, ratd, rat2d, rat2di, tauold)
  !

  return
101 format(1x,4d12.5,i6)
103 format(1x,6d12.5,i6)
104 format(1x,2d12.5,i6)

105 format(1x,3d12.5,i6)
end subroutine vcsmd

subroutine delete_if_present(filename)

   character (len=*) :: filename
   logical           :: exst, opnd
   integer           :: iunit

   inquire (file = filename, exist = exst)

   if (exst) then
      do iunit = 99, 1, - 1
         inquire (unit = iunit, opened = opnd)
         if (.not.opnd) goto 10
      enddo
      call errore ('delete_if_present', 'free unit not found?!?', 1)
10    continue
      open  (unit=iunit, file= filename , status='old')
      close (unit=iunit, status = 'delete')
      write (6,*) 'WARNING: ',filename,' file was present; old file deleted '
   end if

end subroutine
