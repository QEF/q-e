!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!---------------------------------------------------------------------
program matdyn
  !-----------------------------------------------------------------------
  !  this program calculates the phonon frequencies for a list of generic 
  !  q vectors starting from the interatomic force constants generated 
  !  from the dynamical matrices as written by DFPT phonon code through 
  !  the companion program q2r
  !
  !  matdyn can generate a supercell of the original cell for mass
  !  approximation calculation. If supercell data are not specified
  !  in input, the unit cell, lattice vectors, atom types and positions
  !  are read from the force constant file
  !
  !  Input cards: namelist &input
  !     flfrc     file produced by q2r containing force constants (needed)
  !     asr       if .true. use Acoustic Sum Rules (default: .false., no ASR)
  !     dos       if .true. calculate phonon Density of States (DOS)
  !               using tetrahedra and a uniform q-point grid (see below)
  !               NB: may not work properly in noncubic materials
  !               if .false. calculate phonon bands from the list of q-points
  !               supplied in input
  !     nk1,nk2,nk3  uniform q-point grid for DOS calculation (includes q=0)
  !     deltaE    energy step, in cm^(-1), for DOS calculation: from min
  !               to max phonon energy (default: 1 cm^(-1))
  !     fldos     output file for dos (default: 'matdyn.dos')
  !     flfrq     output file for frequencies (default: 'matdyn.freq')
  !     flvec     output file for normal modes (default: 'matdyn.modes')
  !     at        supercell lattice vectors - must form a superlattice of the 
  !               original lattice
  !     l1,l2,l3  supercell lattice vectors are original cell vectors
  !               multiplied by l1, l2, l3 respectively
  !     ntyp      number of atom types in the supercell
  !     amass     masses of atoms in the supercell
  !     readtau   read  atomic positions of the supercell from input
  !               (used to specify different masses)
  !     fltau     write atomic positions of the supercell to file "fltau"
  !               (default: fltau=' ', do not write)
  !
  !  if (readtau) atom types and positions in the supercell follow:
  !     (tau(i,na),i=1,3), ityp(na)
  !  Then, if (.not.dos) :
  !     nq         number of q-points
  !     (q(i,n), i=1,3)    nq q-points in 2pi/a units
  !  If q = 0, the direction qhat (q=>0) for the non-analytic part
  !  is extracted from the sequence of q-points as follows:
  !     qhat = q(n) - q(n-1)  or   qhat = q(n) - q(n+1) 
  !  depending on which one is available and nonzero.
  !  For low-symmetry crystals, specify twice q = 0 in the list
  !  if you want to have q = 0 results for two different directions
  !
  implicit none
  !
  ! variables *_blk refer to the original cell, other variables
  ! to the (super)cell (which may coincide with the original cell)
  !
  integer, parameter:: nax=16, nax_blk=16, nrx=8, nqx=500, nrwsx=200, &
                       nax3=3*nax
  real(kind=8), parameter :: eps=1.0e-6,  rydcm1 = 13.6058*8065.5, &
       amconv = 1.66042e-24/9.1095e-28*0.5
  integer :: nr1, nr2, nr3, nsc, nk1, nk2, nk3, ntetra, ibrav
  character(len=30) :: flfrc, flfrq, flvec, fltau, fldos
  logical :: asr, dos, has_zstar
  complex(kind=8) :: dyn(3,3,nax,nax), dyn_blk(3,3,nax_blk,nax_blk)
  complex(kind=8) :: z(3*nax,3*nax)         ! eigenvalues
  real(kind=8) :: frc(nrx,nrx,nrx,3,3,nax_blk,nax_blk) ! force constants
  real(kind=8) :: at(3,3), bg(3,3), omega,   &! cell parameters and volume
                  alat, tau(3,nax),          &! atomic positions 
                  at_blk(3,3), bg_blk(3,3),  &! original cell
                  omega_blk,                 &! original cell volume
                  tau_blk(3,nax_blk),        &! original atomic positions 
                  epsil(3,3),                &! dielectric tensor
                  zeu(3,3,nax_blk),          &! effective charges
                  amass(nax),                 &! atomic masses
                  amass_blk(nax_blk),         &! original atomic masses
                  q(3,nqx),                   &! list of q-points
                  w2(3*nax,nqx),              &! frequencies (square)
                  freq(3*nax,nqx),            &! frequencies
                  atws(3,3),      &! lattice vector for WS initialization
                  rws(0:3,nrwsx)   ! nearest neighbor list, rws(0,*) = norm^2
  real(kind=8), allocatable:: tetra(:,:)
  !
  integer :: nat, nat_blk,                 & 
             ityp_blk(nax_blk), ityp(nax), &
             ntyp, ntyp_blk,               &
             itau_blk(nax),                &
             l1, l2, l3,                   &! supercell dimensions
             nrws                          ! number of nearest neighbor
  !
  logical :: readtau
  !
  real(kind=8) :: qhat(3), qh, deltaE, Emin, Emax, E, DOSofE(1)
  integer :: n, i, j, it, nq, na, nb, ndos, iout
  namelist /input/ flfrc, amass, asr, flfrq, flvec, at, dos, deltaE,  &
       &           fldos, nk1, nk2, nk3, l1, l2, l3, ntyp, readtau, fltau

  CHARACTER (LEN=80)  :: input_file
  INTEGER             :: nargs, iiarg, ierr

  !
  ! set namelist default
  !
  dos = .false.
  deltaE = 1.0
  nk1 = 0 
  nk2 = 0 
  nk3 = 0 
  asr  =.false.
  readtau=.false.
  flfrc=' '
  fldos='matdyn.dos'
  flfrq='matdyn.freq'
  flvec='matdyn.modes'
  fltau=' '
  amass(:) =0.d0
  amass_blk(:) =0.d0
  at(:,:) = 0.d0
  ntyp=0
  l1=1
  l2=1
  l3=1

  !
  ! ... Input from file ?
  !
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

  !
  read (5,input)
  !
  ! convert masses to atomic units
  !
  amass(:) = amass(:) * amconv
  !
  ! read force constants 
  !
  call readfc ( flfrc, nr1, nr2, nr3, nrx, frc, epsil, zeu, nat_blk, &
       nax_blk, ibrav, alat, tau_blk, at_blk, ntyp_blk, ityp_blk,    &
       amass_blk, omega_blk, has_zstar)
  !
  call recips ( at_blk(1,1),at_blk(1,2),at_blk(1,3),  &
                bg_blk(1,1),bg_blk(1,2),bg_blk(1,3) )
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! set up (super)cell
  !
  ! types of atoms
  ! 
  if (ntyp < 0) then
     call errore ('matdyn','wrong ntyp ', abs(ntyp))
  else if (ntyp == 0) then
     ntyp=ntyp_blk
  end if
  !
  ! masses (for mass approximation)
  ! 
  do it=1,ntyp
     if (amass(it) < 0.d0) then
        call errore ('matdyn','wrong mass in the namelist',it)
     else if (amass(it) == 0.d0) then
        if (it.le.ntyp_blk) then
           write (*,'(a,i3,a,a)') ' mass for atomic type ',it,      &
                &                     ' not given; uses mass from file ',flfrc
           amass(it) = amass_blk(it)
        else
           call errore ('matdyn','missing mass in the namelist',it)
        end if
     end if
  end do
  !
  ! lattice vectors
  !
  if (SUM(abs(at(:,:))) == 0.d0) then
     if (l1.le.0 .or. l2.le.0 .or. l3.le.0) call                    &
          &             errore ('matdyn',' wrong l1,l2 or l3',1)
     at(:,1) = at_blk(:,1)*dfloat(l1)
     at(:,2) = at_blk(:,2)*dfloat(l2)
     at(:,3) = at_blk(:,3)*dfloat(l3)
  end if
  !
  call check_at(at,bg_blk,alat,omega)
  !
  ! the supercell contains "nsc" times the original unit cell
  !
  nsc = nint(omega/omega_blk)
  if (abs(omega/omega_blk-nsc) > eps) &
       call errore ('matdyn', 'volume ratio not integer', 1)
  !
  ! read/generate atomic positions of the (super)cell
  !
  nat = nat_blk * nsc
  if (nat.gt.nax) call errore ('matdyn','nat.gt.nax',nat)
  !
  if (readtau) then
     call read_tau (nat,nat_blk,ntyp,bg_blk,tau,tau_blk,ityp,itau_blk)
  else
     call set_tau  (nat,nat_blk,at,at_blk,tau,tau_blk,ityp,ityp_blk,itau_blk)
  endif
  !
  if (fltau.ne.' ') call write_tau(fltau,nat,tau,ityp)
  !
  ! reciprocal lattice vectors
  !
  call recips (at(1,1),at(1,2),at(1,3),bg(1,1),bg(1,2),bg(1,3))
  !
  ! build the WS cell corresponding to the force constant grid
  !
  atws(:,1) = at_blk(:,1)*dfloat(nr1)
  atws(:,2) = at_blk(:,2)*dfloat(nr2)
  atws(:,3) = at_blk(:,3)*dfloat(nr3)
  ! initialize WS r-vectors
  call wsinit(rws,nrwsx,nrws,atws)
  !
  ! end of (super)cell setup
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (dos) then
     if (nk1 < 1 .or. nk2 < 1 .or. nk3 < 1) &
          call errore  ('matdyn','specify correct q-point grid!',1)
     ntetra = 6 * nk1 * nk2 * nk3
     allocate ( tetra(4,ntetra) )
     call gen_qpoints (ibrav, at, bg, nat, tau, ityp, nk1, nk2, nk3, &
          ntetra, nqx, nq, q, tetra)
  else
  !
  ! read q-point list
  !
     read (5,*) nq
     if (nq.gt.nqx) call errore ('matdyn','too many k-points',nq)
     do n = 1,nq
        read (5,*) (q(i,n),i=1,3)
     end do
  end if
  !
  if(asr) call set_asr(nr1,nr2,nr3,nrx,frc,zeu,nat_blk,nax_blk)
  !
  if (flvec.eq.' ') then
     iout=6
  else
     iout=4
     open (unit=iout,file=flvec,status='unknown',form='formatted')
  end if
  do n=1, nq
     dyn(:,:,:,:) = (0.d0, 0.d0)

     call setupmat (q(1,n),dyn,nat,nax,at,bg,tau,itau_blk,nsc,alat, &
             dyn_blk,nat_blk,nax_blk,at_blk,bg_blk,tau_blk,omega_blk,  &
             epsil,zeu,frc,nr1,nr2,nr3,nrx,has_zstar,rws,nrws)

     if (q(1,n)==0.d0 .and. q(2,n)==0.d0 .and. q(3,n)==0.d0) then
        !
        ! q = 0 : we need the direction q => 0 for the non-analytic part
        !
        if ( (n == 1 .and. nq > 1) .or. &
             (n > 1 .and. n < nq .and.  &
              q(1,n-1)==0.d0.and.q(2,n-1)==0.d0.and.q(3,n-1)==0.d0) ) then
              ! if q is the first point in the list, or
              ! if preceding q is also 0 :
           qhat(:) = q(:,n) - q(:,n+1)
        else if ( n > 1 ) then
              ! if q is not the first point in the list
           qhat(:) = q(:,n) - q(:,n-1)
        else
           qhat(:) = 0.d0
        end if
        qh = sqrt(qhat(1)**2+qhat(2)**2+qhat(3)**2)
        if (qh /= 0.d0) qhat(:) = qhat(:) / qh
        if (qh /= 0.d0 .and. .not. has_zstar) call errore  &
             ('matdyn','non-analytic term for q=0 missing !', -1)
        !
        call nonanal (nax,nat,dyn,qhat,itau_blk,nax_blk,epsil,zeu,omega)
        !
     end if
     !
     call dyndiag(nax,nat,amass,ityp,dyn,w2(1,n),z)
     !
     call writemodes(nax,nat,q(1,n),w2(1,n),z,iout)
     !
  end do
  !
  if(iout .ne. 6) close(unit=iout)
  !
  do n=1,nq
     ! freq(i,n) = frequencies in cm^(-1)
     !             negative sign if omega^2 is negative
     do i=1,3*nat
        freq(i,n)= sqrt(abs(w2(i,n)))*rydcm1
        if (w2(i,n).lt.0.0) freq(i,n) = -freq(i,n)
     end do
  end do
  !
  if(flfrq.ne.' ') then
     open (unit=2,file=flfrq ,status='unknown',form='formatted')
     write(2,*) nq, 3*nat
     do n=1, nq
        write(2,'(6f10.4)') (freq(i,n),i=1,3*nat)
     end do
     close(unit=2)
  end if
  !
  if (dos) then
     Emin = 0.0 
     Emax = 0.0
     do n=1,nq
        do i=1, 3*nat
           Emin = min (Emin, freq(i,n))
           Emax = max (Emax, freq(i,n))
        end do
     end do
     !
     ndos = nint ( (Emax - Emin) / DeltaE+0.500001)  
     open (unit=2,file=fldos,status='unknown',form='formatted')
     do n= 1, ndos  
        E = Emin + (n - 1) * DeltaE  
        call dos_t(freq, 1, 3*nax, 3*nat, nq, ntetra, tetra, E, DOSofE)
        write (2, '(2e12.4)') E, DOSofE (1)
     end do
     close(unit=2)
  end if
  !
  stop
end program matdyn
!
!-----------------------------------------------------------------------
subroutine readfc (flfrc,nr1,nr2,nr3,nrx,frc,epsil,zeu,nat,nax,    &
                   ibrav,alat,tau,at,ntyp,ityp,amass,omega,has_zstar)
  !-----------------------------------------------------------------------
  !
  implicit none
  ! I/O variable
  character(len=30) flfrc
  integer ibrav, nr1,nr2,nr3,nrx, nat, nax, ntyp, ityp(nax)
  real(kind=8) frc(nrx,nrx,nrx,3,3,nax,nax), epsil(3,3),zeu(3,3,nax)
  real(kind=8) alat, at(3,3), tau(3,nax)
  logical has_zstar
  ! local variables
  integer i, j, na, nb, m1,m2,m3
  integer ibid, jbid, nabid, nbbid, m1bid,m2bid,m3bid
  real(kind=8) amass(nax), amass_from_file, celldm(6), omega
  integer nt
  character(len=3) atm
  !
  !
  open (unit=1,file=flfrc,status='old',form='formatted')
  !
  !
  !  read cell data
  !
  read(1,*) ntyp,nat,ibrav,(celldm(i),i=1,6)
  if (nat.gt.nax) call errore ('readfc','too many atoms',nat)
  call latgen(ibrav,celldm,at(1,1),at(1,2),at(1,3),omega)
  alat = celldm(1)
  at = at / alat !  bring at in units of alat
  call volume(alat,at(1,1),at(1,2),at(1,3),omega)
  !
  !  read atomic types, positions and masses
  !
  do nt = 1,ntyp
     read(1,*) i,atm,amass_from_file
     if (i.ne.nt) call errore ('readfc','wrong data read',nt)
     if (amass(nt).eq.0.d0) then
        amass(nt) = amass_from_file
     else
        write(*,*) 'for atomic type',nt,' mass from file not used'
     end if
  end do
  do na=1,nat
     read(1,*) i,ityp(na),(tau(j,na),j=1,3)
     if (i.ne.na) call errore ('readfc','wrong data read',na)
  end do
  !
  !  read macroscopic variable
  !
  read (1,*) has_zstar
  if (has_zstar) then
     read(1,*) ((epsil(i,j),j=1,3),i=1,3)
     do na=1,nat
        read(1,*) 
        read(1,*) ((zeu(i,j,na),j=1,3),i=1,3)
     end do
  else
     zeu  (:,:,:) = 0.d0
     epsil(:,:) = 0.d0
  end if
  !
  !  read real space part
  !
  read (1,*) nr1,nr2,nr3
  if (nr1.gt.nrx)  call errore ('readin','nr1 .gt. nrx ',+1)
  if (nr2.gt.nrx)  call errore ('readin','nr2 .gt. nrx ',+1)
  if (nr3.gt.nrx)  call errore ('readin','nr3 .gt. nrx ',+1)
  !
  if(nat.gt.nax) call errore  ('readfc','nax too small', nat)
  !
  !  read real-space interatomic force constants
  !
  frc(:,:,:,:,:,:,:) = 0.d0
  do i=1,3
     do j=1,3
        do na=1,nat
           do nb=1,nat
              read (1,*) ibid, jbid, nabid, nbbid
              if(i .ne.ibid  .or. j .ne.jbid .or.                   &
                 na.ne.nabid .or. nb.ne.nbbid)                      &
                 call errore  ('readfc','error in reading',1)
              read (1,*) (((m1bid, m2bid, m3bid,                    &
                          frc(m1,m2,m3,i,j,na,nb),                  &
                           m1=1,nr1),m2=1,nr2),m3=1,nr3)
           end do
        end do
     end do
  end do
  !
  close(unit=1)
  !
  return
end subroutine readfc
!
!-----------------------------------------------------------------------
subroutine frc_blk(nax,dyn,q,tau,nat,                             &
     &                   nr1,nr2,nr3,nrx,frc,at,bg,rws,nrws)
  !-----------------------------------------------------------------------
  ! calculates the dynamical matrix at q from the (short-range part of the)
  ! force constants 
  !
  implicit none
  integer nr1, nr2, nr3, nrx, nax, nat, n1, n2, n3, &
          ipol, jpol, na, nb, m1, m2, m3, nint, i,j, nrws
  complex(kind=8) dyn(3,3,nax,nax), cmplx
  real(kind=8) frc(nrx,nrx,nrx,3,3,nax,nax), tau(3,nax), q(3), arg, &
               at(3,3), bg(3,3), r(3), weight, r_ws(3),  &
               total_weight, rws(0:3,nrws)
  real(kind=8), parameter:: tpi = 2.0*3.14159265358979d0
  real(kind=8), external :: wsweight
  !
  do na=1, nat
     do nb=1, nat
        total_weight=0.0d0
        do n1=-2*nrx,2*nrx
           do n2=-2*nrx,2*nrx
              do n3=-2*nrx,2*nrx
                 !
                 ! SUM OVER R VECTORS IN THE SUPERCELL - VERY VERY SAFE RANGE!
                 !
                 do i=1, 3
                    r(i) = n1*at(i,1)+n2*at(i,2)+n3*at(i,3)
                    r_ws(i) = r(i) + tau(i,na)-tau(i,nb)
                 end do
                 weight = wsweight(r_ws,rws,nrws)
                 if (weight .gt. 0.0) then
                    !
                    ! FIND THE VECTOR CORRESPONDING TO R IN THE ORIGINAL CELL
                    !
                    m1 = mod(n1+1,nr1)
                    if(m1.le.0) m1=m1+nr1
                    m2 = mod(n2+1,nr2)
                    if(m2.le.0) m2=m2+nr2
                    m3 = mod(n3+1,nr3)
                    if(m3.le.0) m3=m3+nr3
                    !
                    ! FOURIER TRANSFORM
                    !
                    arg = tpi*(q(1)*r(1) + q(2)*r(2) + q(3)*r(3))
                    do ipol=1, 3
                       do jpol=1, 3
                          dyn(ipol,jpol,na,nb) =                 &
                               dyn(ipol,jpol,na,nb) +            &
                               frc(m1,m2,m3,ipol,jpol,na,nb)     &
                               *cmplx(cos(arg),-sin(arg))*weight
                       end do
                    end do
                 end if
                 total_weight=total_weight + weight
              end do
           end do
        end do
        if (abs(total_weight-nr1*nr2*nr3).gt.1.0d-8) then
           write(*,*) total_weight
           call errore ('frc_blk','wrong total_weight',1)
        end if
     end do
  end do
  !
  return
end subroutine frc_blk
!
!-----------------------------------------------------------------------
subroutine setupmat (q,dyn,nat,nax,at,bg,tau,itau_blk,nsc,alat, &
     &         dyn_blk,nat_blk,nax_blk,at_blk,bg_blk,tau_blk,omega_blk, &
     &                 epsil,zeu,frc,nr1,nr2,nr3,nrx,has_zstar,rws,nrws)
  !-----------------------------------------------------------------------
  ! compute the dynamical matrix (the analytic part only)
  !
  implicit none
  real(kind=8), parameter :: tpi=2.d0*3.14159265358979d0
  !
  ! I/O variables
  !
  integer:: nr1, nr2, nr3, nrx, nax, nat, nat_blk, nax_blk, &
       &    nsc, nrws, itau_blk(nat)
  real(kind=8) :: q(3), tau(3,nax), at(3,3), bg(3,3), alat,      &
                  epsil(3,3), zeu(3,3,nax_blk), rws(0:3,nrws),   &
                  frc(nrx,nrx,nrx,3,3,nax_blk,nax_blk)
  real(kind=8) :: tau_blk(3,nax_blk), at_blk(3,3), bg_blk(3,3), omega_blk
  complex(kind=8) dyn_blk(3,3,nax_blk,nax_blk)
  complex(kind=8) ::  dyn(3,3,nax,nax)
  logical has_zstar
  !
  ! local variables
  !
  real(kind=8) :: arg
  complex(kind=8) :: cfac(nat)
  integer :: i,j,k, na,nb, na_blk, nb_blk, iq
  real(kind=8) qp(3), qbid(3,nsc) ! automatic array
  !
  !
  call q_gen(nsc,qbid,at_blk,bg_blk,at,bg)
  !
  do iq=1,nsc
     !
     do k=1,3
        qp(k)= q(k) + qbid(k,iq)
     end do
     !
     dyn_blk(:,:,:,:) = (0.d0,0.d0)
     call frc_blk (nax_blk,dyn_blk,qp,tau_blk,nat_blk,              &
          &              nr1,nr2,nr3,nrx,frc,at_blk,bg_blk,rws,nrws)
     if (has_zstar) &
          call rgd_blk(nax_blk,nat_blk,dyn_blk,qp,tau_blk,   &
                       epsil,zeu,bg_blk,omega_blk,+1.d0)
     !
     do na=1,nat
        na_blk = itau_blk(na)
        do nb=1,nat
           nb_blk = itau_blk(nb)
           !
           arg=tpi* ( qp(1) * ( (tau(1,na)-tau_blk(1,na_blk)) -   &
                                (tau(1,nb)-tau_blk(1,nb_blk)) ) + &
                      qp(2) * ( (tau(2,na)-tau_blk(2,na_blk)) -   &
                                (tau(2,nb)-tau_blk(2,nb_blk)) ) + &
                      qp(3) * ( (tau(3,na)-tau_blk(3,na_blk)) -   &
                                (tau(3,nb)-tau_blk(3,nb_blk)) ) )
           !
           cfac(nb) = cmplx(cos(arg),sin(arg))/nsc
           !
        end do ! nb
        !
        do i=1,3
           do j=1,3
              !
              do nb=1,nat
                 nb_blk = itau_blk(nb)
                 dyn(i,j,na,nb) = dyn(i,j,na,nb) + cfac(nb) * &
                      dyn_blk(i,j,na_blk,nb_blk)
              end do ! nb
              !
           end do ! j
        end do ! i
     end do ! na
     !
  end do ! iq
  !
  return
end subroutine setupmat
!----------------------------------------------------------------------
subroutine set_asr(nr1,nr2,nr3,nrx,frc,zeu,nat,nax)
  !-----------------------------------------------------------------------
  !
  implicit none
  integer nr1, nr2, nr3, nrx, nr, i, j, na, nb, n1,n2,n3, nat, nax
  real(kind=8) frc(nrx,nrx,nrx,3,3,nax,nax), sum, zeu(3,3,nax)
  !
  ! Acoustic Sum Rule on effective charges
  !
  do i=1,3
     do j=1,3
        sum=0.0
        do na=1,nat
           sum = sum + zeu(i,j,na)
        end do
        do na=1,nat
           zeu(i,j,na) = zeu(i,j,na) - sum/nat
        end do
     end do
  end do
  !
  ! Acoustic Sum Rule on force constants in real space
  !
  do i=1,3
     do j=1,3
        do na=1,nat
           sum=0.0
           do nb=1,nat
              do n1=1,nr1
                 do n2=1,nr2
                    do n3=1,nr3
                       sum=sum+frc(n1,n2,n3,i,j,na,nb)
                    end do
                 end do
              end do
           end do
           frc(1,1,1,i,j,na,na) = frc(1,1,1,i,j,na,na) - sum
           !               write(6,*) ' na, i, j, sum = ',na,i,j,sum
        end do
     end do
  end do
  !
  return
end subroutine set_asr
!
!-----------------------------------------------------------------------
subroutine q_gen(nsc,qbid,at_blk,bg_blk,at,bg)
  !-----------------------------------------------------------------------
  ! generate list of q (qbid) that are G-vectors of the supercell
  ! but not of the bulk
  !
  implicit none
  integer :: nsc
  real(kind=8) qbid(3,nsc), at_blk(3,3), bg_blk(3,3), at(3,3), bg(3,3)
  !
  integer, parameter:: nr1=4, nr2=4, nr3=4, &
                       nrm=(2*nr1+1)*(2*nr2+1)*(2*nr3+1)
  real(kind=8), parameter:: eps=1.0e-7
  integer :: i, j, k,i1, i2, i3, idum(nrm), iq
  real(kind=8) :: qnorm(nrm), qbd(3,nrm) ,qwork(3), delta
  logical lbho
  !
  i = 0
  do i1=-nr1,nr1
     do i2=-nr2,nr2
        do i3=-nr3,nr3
           i = i + 1
           do j=1,3
              qwork(j) = i1*bg(j,1) + i2*bg(j,2) + i3*bg(j,3)
           end do ! j
           !
           qnorm(i)  = qwork(1)**2 + qwork(2)**2 + qwork(3)**2
           !
           do j=1,3
              !
              qbd(j,i) = at_blk(1,j)*qwork(1) + &
                         at_blk(2,j)*qwork(2) + &
                         at_blk(3,j)*qwork(3)
           end do ! j
           !
           idum(i) = 1
           !
        end do ! i3
     end do ! i2
  end do ! i1
  !
  do i=1,nrm-1
     if (idum(i).eq.1) then
        do j=i+1,nrm
           if (idum(j).eq.1) then
              lbho=.true.
              do k=1,3
                 delta = qbd(k,i)-qbd(k,j)
                 lbho = lbho.and. (abs(nint(delta)-delta).lt.eps)
              end do ! k
              if (lbho) then
                 if(qnorm(i).gt.qnorm(j)) then
                    qbd(1,i) = qbd(1,j)
                    qbd(2,i) = qbd(2,j)
                    qbd(3,i) = qbd(3,j)
                    qnorm(i) = qnorm(j)
                 end if
                 idum(j) = 0
              end if
           end if
        end do ! j
     end if
  end do ! i
  !
  iq = 0
  do i=1,nrm
     if (idum(i).eq.1) then
        iq=iq+1
        qbid(1,iq)= bg_blk(1,1)*qbd(1,i) +  &
                    bg_blk(1,2)*qbd(2,i) +  &
                    bg_blk(1,3)*qbd(3,i)
        qbid(2,iq)= bg_blk(2,1)*qbd(1,i) +  &
                    bg_blk(2,2)*qbd(2,i) +  &
                    bg_blk(2,3)*qbd(3,i)
        qbid(3,iq)= bg_blk(3,1)*qbd(1,i) +  &
                    bg_blk(3,2)*qbd(2,i) +  &
                    bg_blk(3,3)*qbd(3,i)
     end if
  end do ! i
  !
  if (iq.ne.nsc) call errore('q_gen',' probably nr1,nr2,nr3 too small ', iq)
  return
end subroutine q_gen
!
!-----------------------------------------------------------------------
subroutine check_at(at,bg_blk,alat,omega)
  !-----------------------------------------------------------------------
  implicit none
  !
  real(kind=8) :: at(3,3), bg_blk(3,3), alat, omega
  real(kind=8) :: work(3,3)
  integer :: i,j
  real(kind=8), parameter :: small=1.d-6
  !
  work(:,:) = at(:,:)
  call cryst_to_cart(3,work,bg_blk,-1)
  !
  do j=1,3
     do i =1,3
        if ( abs(work(i,j)-nint(work(i,j))) > small) then
           write (*,'(3f9.4)') work(:,:)
           call errore ('check_at','at not multiple of at_blk',1)
        end if
     end do
  end do
  !
  omega =alat**3 * abs(at(1,1)*(at(2,2)*at(3,3)-at(3,2)*at(2,3))- &
                       at(1,2)*(at(2,1)*at(3,3)-at(2,3)*at(3,1))+ &
                       at(1,3)*(at(2,1)*at(3,2)-at(2,2)*at(3,1)))
  !
  return
end subroutine check_at
!
!-----------------------------------------------------------------------
subroutine set_tau                                                &
     &        (nat,nat_blk,at,at_blk,tau,tau_blk,ityp,ityp_blk,itau_blk)
  !-----------------------------------------------------------------------
  implicit none
  integer nat, nat_blk,ityp(nat),ityp_blk(nat_blk), itau_blk(nat)
  real(kind=8) at(3,3),at_blk(3,3),tau(3,nat),tau_blk(3,nat_blk)
  !
  real(kind=8) bg(3,3), r(3) ! work vectors
  integer i,i1,i2,i3,na,na_blk
  real(kind=8) small
  integer NN1,NN2,NN3
  parameter (NN1=8, NN2=8, NN3=8, small=1.d-8)
  !
  call recips (at(1,1),at(1,2),at(1,3),bg(1,1),bg(1,2),bg(1,3))
  !
  na = 0
  !
  do i1 = -NN1,NN1
     do i2 = -NN2,NN2
        do i3 = -NN3,NN3
           r(1) = i1*at_blk(1,1) + i2*at_blk(1,2) + i3*at_blk(1,3)
           r(2) = i1*at_blk(2,1) + i2*at_blk(2,2) + i3*at_blk(2,3)
           r(3) = i1*at_blk(3,1) + i2*at_blk(3,2) + i3*at_blk(3,3)
           call cryst_to_cart(1,r,bg,-1)
           !
           if ( r(1).gt.-small .and. r(1).lt.1.d0-small .and.          &
                r(2).gt.-small .and. r(2).lt.1.d0-small .and.          &
                r(3).gt.-small .and. r(3).lt.1.d0-small ) then
              call cryst_to_cart(1,r,at,+1)
              !
              do na_blk=1, nat_blk
                 na = na + 1
                 if (na.gt.nat) call errore('set_tau','too many atoms',na)
                 tau(1,na)    = tau_blk(1,na_blk) + r(1)
                 tau(2,na)    = tau_blk(2,na_blk) + r(2)
                 tau(3,na)    = tau_blk(3,na_blk) + r(3)
                 ityp(na)     = ityp_blk(na_blk)
                 itau_blk(na) = na_blk
              end do
              !
           end if
           !
        end do
     end do
  end do
  !
  if (na.ne.nat) call errore('set_tau','too few atoms: increase NNs',na)
  !
  return
end subroutine set_tau
!
!-----------------------------------------------------------------------
subroutine read_tau(nat,nat_blk,ntyp,bg_blk,tau,tau_blk,ityp,itau_blk)
  !---------------------------------------------------------------------
  implicit none
  integer nat, nat_blk, ntyp, ityp(nat),itau_blk(nat)
  real(kind=8) bg_blk(3,3),tau(3,nat),tau_blk(3,nat_blk)
  !
  real(kind=8) r(3) ! work vectors
  integer i,na,na_blk
  !
  real(kind=8) small
  parameter ( small = 1.d-6 )
  !
  do na=1,nat
     read(*,*) (tau(i,na),i=1,3), ityp(na)
     if (ityp(na).le.0 .or. ityp(na) .gt. ntyp) &
          call errore('read_tau',' wrong atomic type', na)
     do na_blk=1,nat_blk
        r(1) = tau(1,na) - tau_blk(1,na_blk)
        r(2) = tau(2,na) - tau_blk(2,na_blk)
        r(3) = tau(3,na) - tau_blk(3,na_blk)
        call cryst_to_cart(1,r,bg_blk,-1)
        if (abs( r(1)-nint(r(1)) ) .lt. small .and.                 &
            abs( r(2)-nint(r(2)) ) .lt. small .and.                 &
            abs( r(3)-nint(r(3)) ) .lt. small ) then
           itau_blk(na) = na_blk
           go to 999
        end if
     end do
     call errore ('read_tau',' wrong atomic position ', na)
999  continue
  end do
  !
  return
end subroutine read_tau
!
!-----------------------------------------------------------------------
subroutine write_tau(fltau,nat,tau,ityp)
  !-----------------------------------------------------------------------
  implicit none
  integer nat, ityp(nat)
  real(kind=8) tau(3,nat)
  character(len=*) fltau
  !
  integer i,na
  !
  open (unit=4,file=fltau, status='new')
  do na=1,nat
     write(4,'(3(f12.6),i3)') (tau(i,na),i=1,3), ityp(na)
  end do
  close (4)
  !
  return 
end subroutine write_tau
!
!-----------------------------------------------------------------------
subroutine gen_qpoints (ibrav, at, bg, nat, tau, ityp, nk1, nk2, nk3, &
     ntetra, nqx, nq, q, tetra)
  !-----------------------------------------------------------------------
  !
  implicit none
  ! input 
  integer :: ibrav, nat, nqx, nk1, nk2, nk3, ntetra, ityp(*)
  real(kind=8) :: at(3,3), bg(3,3), tau(3,nat)
  ! output 
  integer :: nq, tetra(4,ntetra)
  real(kind=8) :: q(3,nqx)
  ! local
  integer :: nrot, nsym, s(3,3,48), ftau(3,48), irt(48,nat)
  logical :: minus_q, invsym
  real(kind=8) :: xqq(3), wk(nqx), mdum(3,nat)
  character(len=45)   ::  sname(48)
  !
  xqq (:) =0.d0
  if (ibrav == 4 .or. ibrav == 5) then  
     !
     !  hexagonal or trigonal bravais lattice
     !
     call hexsym (at, s, sname, nrot)  
  elseif (ibrav >= 1 .and. ibrav <= 14) then  
     !
     !  cubic bravais lattice
     !
     call cubicsym (at, s, sname, nrot)  
  elseif (ibrav == 0) then  
     call errore ('gen_qpoints', 'assuming cubic symmetry',-1)  
     call cubicsym (at, s, sname, nrot)  
  else  
     call errore ('gen_qpoints', 'wrong ibrav', 1)  
  endif
  !
  call kpoint_grid ( nrot, s, bg, nqx, 0,0,0, nk1,nk2,nk3, nq, q, wk)
  !
  call sgama (nrot, nat, s, sname, at, bg, tau, ityp, nsym, 6, &
       6, 6, irt, ftau, nqx, nq, q, wk, invsym, minus_q, xqq, &
       0, 0, .false., mdum)
  
  if (ntetra /= 6 * nk1 * nk2 * nk3) &
       call errore ('gen_qpoints','inconsistent ntetra',1)

  call tetrahedra (nsym, s, minus_q, at, bg, nqx, 0, 0, 0, &
       nk1, nk2, nk3, nq, q, wk, ntetra, tetra)
  !
  return
end subroutine gen_qpoints
!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!--------------------------------------------------------------------
subroutine dos_t (et, nspin, nbndx, nbnd, nks, ntetra, tetra, e, dost)
  !------------------------------------------------------------------
  !
  use kinds, only : DP
  implicit none  
  integer :: nspin, nbndx, nbnd, nks, ntetra, tetra (4, ntetra)  

  real(kind=DP) :: et (nbndx, nks), e, dost (2)  
  integer :: itetra (4), nk, ns, nt, ibnd, i  

  real(kind=DP) :: etetra (4), e1, e2, e3, e4  
  do ns = 1, nspin  
     dost (ns) = 0.0  
     !
     ! nk is used to select k-points with up (ns=1) or down (ns=2) spin
     !
     if (ns.eq.1) then  
        nk = 0  
     else  
        nk = nks / 2  

     endif
     do nt = 1, ntetra  
        do ibnd = 1, nbnd  
           ! these are the energies at the vertexes of the nt-th tetrahedron
           do i = 1, 4  
              etetra (i) = et (ibnd, tetra (i, nt) + nk)  
           enddo
           itetra (1) = 0  
           call hpsort (4, etetra, itetra)  
           e1 = etetra (1)  
           e2 = etetra (2)  
           e3 = etetra (3)  
           e4 = etetra (4)  
           if (e.lt.e4.and.e.ge.e3) then  
              dost (ns) = dost (ns) + 1.d0 / ntetra * (3.0 * (e4 - e) **2 / &
                   (e4 - e1) / (e4 - e2) / (e4 - e3) )
           elseif (e.lt.e3.and.e.ge.e2) then  
              dost (ns) = dost (ns) + 1.d0 / ntetra / (e3 - e1) / (e4 - e1) &
                   * (3.0 * (e2 - e1) + 6.0 * (e-e2) - 3.0 * (e3 - e1 + e4 - e2) &
                   / (e3 - e2) / (e4 - e2) * (e-e2) **2)
           elseif (e.lt.e2.and.e.gt.e1) then  
              dost (ns) = dost (ns) + 1.d0 / ntetra * 3.0 * (e-e1) **2 / &
                   (e2 - e1) / (e3 - e1) / (e4 - e1)
           endif
        enddo


     enddo
     ! add correct spin normalization : 2 for LDA, 1 for LSDA calculations

     dost (ns) = dost (ns) * 2.d0 / nspin  

  enddo
  return  
end subroutine dos_t
