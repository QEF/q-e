!
! Copyright (C) 2001-2004 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "f_defs.h"
!
!---------------------------------------------------------------------
PROGRAM matdyn
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
  USE kinds,      ONLY : DP
  USE mp,         ONLY : mp_start, mp_env, mp_end, mp_barrier
  USE mp_global,  ONLY : nproc, mpime, mp_global_start  
  !
  IMPLICIT NONE
  !
  INTEGER :: gid
  !
  ! variables *_blk refer to the original cell, other variables
  ! to the (super)cell (which may coincide with the original cell)
  !
  INTEGER, PARAMETER:: nax=16, nax_blk=16, nrx=8, nqx=500, nrwsx=200, &
                       nax3=3*nax
  REAL(KIND=DP), PARAMETER :: eps=1.0e-6,  rydcm1 = 13.6058*8065.5, &
       amconv = 1.66042e-24/9.1095e-28*0.5
  INTEGER :: nr1, nr2, nr3, nsc, nk1, nk2, nk3, ntetra, ibrav
  CHARACTER(LEN=256) :: flfrc, flfrq, flvec, fltau, fldos
  LOGICAL :: asr, dos, has_zstar
  COMPLEX(KIND=DP) :: dyn(3,3,nax,nax), dyn_blk(3,3,nax_blk,nax_blk)
  COMPLEX(KIND=DP) :: z(3*nax,3*nax)         ! eigenvalues
  REAL(KIND=DP) :: frc(nrx,nrx,nrx,3,3,nax_blk,nax_blk) ! force constants
  REAL(KIND=DP) :: at(3,3), bg(3,3), omega,   &! cell parameters and volume
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
                  atws(3,3),      &! lattice vector for WS initialization
                  rws(0:3,nrwsx)   ! nearest neighbor list, rws(0,*) = norm^2
  REAL(KIND=DP), ALLOCATABLE:: tetra(:,:), freq(:,:)
  !
  INTEGER :: nat, nat_blk,                 & 
             ityp_blk(nax_blk), ityp(nax), &
             ntyp, ntyp_blk,               &
             itau_blk(nax),                &
             l1, l2, l3,                   &! supercell dimensions
             nrws                          ! number of nearest neighbor
  !
  LOGICAL :: readtau
  !
  REAL(KIND=DP) :: qhat(3), qh, deltaE, Emin, Emax, E, DOSofE(1)
  INTEGER :: n, i, j, it, nq, na, nb, ndos, iout
  NAMELIST /input/ flfrc, amass, asr, flfrq, flvec, at, dos, deltaE,  &
       &           fldos, nk1, nk2, nk3, l1, l2, l3, ntyp, readtau, fltau

  CHARACTER (LEN=256)  :: input_file
  INTEGER             :: nargs, iiarg, ierr, ILEN
  INTEGER, EXTERNAL   :: iargc
  !
  !
  CALL mp_start()
  !
  CALL mp_env( nproc, mpime, gid )
  !
  IF ( mpime == 0 ) THEN
     !
     ! ... all calculations are done by the first cpu
     !
     ! set namelist default
     !
     dos = .FALSE.
     deltaE = 1.0
     nk1 = 0 
     nk2 = 0 
     nk3 = 0 
     asr  =.FALSE.
     readtau=.FALSE.
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
     READ (5,input)
     !
     ! convert masses to atomic units
     !
     amass(:) = amass(:) * amconv
     !
     ! read force constants 
     !
     CALL readfc ( flfrc, nr1, nr2, nr3, nrx, frc, epsil, zeu, nat_blk, &
          nax_blk, ibrav, alat, tau_blk, at_blk, ntyp_blk, ityp_blk,    &
          amass_blk, omega_blk, has_zstar)
     !
     CALL recips ( at_blk(1,1),at_blk(1,2),at_blk(1,3),  &
          bg_blk(1,1),bg_blk(1,2),bg_blk(1,3) )
     !
     ! set up (super)cell
     !
     ! types of atoms
     ! 
     IF (ntyp < 0) THEN
        CALL errore ('matdyn','wrong ntyp ', ABS(ntyp))
     ELSE IF (ntyp == 0) THEN
        ntyp=ntyp_blk
     END IF
     !
     ! masses (for mass approximation)
     ! 
     DO it=1,ntyp
        IF (amass(it) < 0.d0) THEN
           CALL errore ('matdyn','wrong mass in the namelist',it)
        ELSE IF (amass(it) == 0.d0) THEN
           IF (it.LE.ntyp_blk) THEN
              WRITE (*,'(a,i3,a,a)') ' mass for atomic type ',it,      &
                   &                     ' not given; uses mass from file ',flfrc
              amass(it) = amass_blk(it)
           ELSE
              CALL errore ('matdyn','missing mass in the namelist',it)
           END IF
        END IF
     END DO
     !
     ! lattice vectors
     !
     IF (SUM(ABS(at(:,:))) == 0.d0) THEN
        IF (l1.LE.0 .OR. l2.LE.0 .OR. l3.LE.0) CALL                    &
             &             errore ('matdyn',' wrong l1,l2 or l3',1)
        at(:,1) = at_blk(:,1)*DBLE(l1)
        at(:,2) = at_blk(:,2)*DBLE(l2)
        at(:,3) = at_blk(:,3)*DBLE(l3)
     END IF
     !
     CALL check_at(at,bg_blk,alat,omega)
     !
     ! the supercell contains "nsc" times the original unit cell
     !
     nsc = NINT(omega/omega_blk)
     IF (ABS(omega/omega_blk-nsc) > eps) &
          CALL errore ('matdyn', 'volume ratio not integer', 1)
     !
     ! read/generate atomic positions of the (super)cell
     !
     nat = nat_blk * nsc
     IF (nat.GT.nax) CALL errore ('matdyn','nat.gt.nax',nat)
     !
     IF (readtau) THEN
        CALL read_tau (nat,nat_blk,ntyp,bg_blk,tau,tau_blk,ityp,itau_blk)
     ELSE
        CALL set_tau  (nat,nat_blk,at,at_blk,tau,tau_blk,ityp,ityp_blk,itau_blk)
     ENDIF
     !
     IF (fltau.NE.' ') CALL write_tau(fltau,nat,tau,ityp)
     !
     ! reciprocal lattice vectors
     !
     CALL recips (at(1,1),at(1,2),at(1,3),bg(1,1),bg(1,2),bg(1,3))
     !
     ! build the WS cell corresponding to the force constant grid
     !
     atws(:,1) = at_blk(:,1)*DBLE(nr1)
     atws(:,2) = at_blk(:,2)*DBLE(nr2)
     atws(:,3) = at_blk(:,3)*DBLE(nr3)
     ! initialize WS r-vectors
     CALL wsinit(rws,nrwsx,nrws,atws)
     !
     ! end of (super)cell setup
     !
     IF (dos) THEN
        IF (nk1 < 1 .OR. nk2 < 1 .OR. nk3 < 1) &
             CALL errore  ('matdyn','specify correct q-point grid!',1)
        ntetra = 6 * nk1 * nk2 * nk3
        ALLOCATE ( tetra(4,ntetra) )
        CALL gen_qpoints (ibrav, at, bg, nat, tau, ityp, nk1, nk2, nk3, &
             ntetra, nqx, nq, q, tetra)
     ELSE
        !
        ! read q-point list
        !
        READ (5,*) nq
        IF (nq.GT.nqx) CALL errore ('matdyn','too many k-points',nq)
        DO n = 1,nq
           READ (5,*) (q(i,n),i=1,3)
        END DO
     END IF
     !
     IF(asr) CALL set_asr(nr1,nr2,nr3,nrx,frc,zeu,nat_blk,nax_blk)
     !
     IF (flvec.EQ.' ') THEN
        iout=6
     ELSE
        iout=4
        OPEN (unit=iout,file=flvec,status='unknown',form='formatted')
     END IF
     DO n=1, nq
        dyn(:,:,:,:) = (0.d0, 0.d0)

        CALL setupmat (q(1,n),dyn,nat,nax,at,bg,tau,itau_blk,nsc,alat, &
             dyn_blk,nat_blk,nax_blk,at_blk,bg_blk,tau_blk,omega_blk,  &
             epsil,zeu,frc,nr1,nr2,nr3,nrx,has_zstar,rws,nrws)

        IF (q(1,n)==0.d0 .AND. q(2,n)==0.d0 .AND. q(3,n)==0.d0) THEN
           !
           ! q = 0 : we need the direction q => 0 for the non-analytic part
           !
           IF ( (n == 1 .AND. nq > 1) .OR. &
                (n > 1 .AND. n < nq .AND.  &
                q(1,n-1)==0.d0.AND.q(2,n-1)==0.d0.AND.q(3,n-1)==0.d0) ) THEN
              ! if q is the first point in the list, or
              ! if preceding q is also 0 :
              qhat(:) = q(:,n) - q(:,n+1)
           ELSE IF ( n > 1 ) THEN
              ! if q is not the first point in the list
              qhat(:) = q(:,n) - q(:,n-1)
           ELSE
              qhat(:) = 0.d0
           END IF
           qh = SQRT(qhat(1)**2+qhat(2)**2+qhat(3)**2)
           IF (qh /= 0.d0) qhat(:) = qhat(:) / qh
           IF (qh /= 0.d0 .AND. .NOT. has_zstar) CALL errore  &
                ('matdyn','non-analytic term for q=0 missing !', -1)
           !
           CALL nonanal (nax,nat,dyn,qhat,itau_blk,nax_blk,epsil,zeu,omega)
           !
        END IF
        !
        CALL dyndiag(nax,nat,amass,ityp,dyn,w2(1,n),z)
        !
        CALL writemodes(nax,nat,q(1,n),w2(1,n),z,iout)
        !
     END DO
     !
     IF(iout .NE. 6) CLOSE(unit=iout)
     !
     ALLOCATE (freq(3*nat, nq))
     DO n=1,nq
        ! freq(i,n) = frequencies in cm^(-1)
        !             negative sign if omega^2 is negative
        DO i=1,3*nat
           freq(i,n)= SQRT(ABS(w2(i,n)))*rydcm1
           IF (w2(i,n).LT.0.0) freq(i,n) = -freq(i,n)
        END DO
     END DO
     !
     IF(flfrq.NE.' ') THEN
        OPEN (unit=2,file=flfrq ,status='unknown',form='formatted')
        WRITE(2,*) nq, 3*nat
        DO n=1, nq
           WRITE(2,'(6f10.4)') (freq(i,n),i=1,3*nat)
        END DO
        CLOSE(unit=2)
     END IF
     !
     IF (dos) THEN
        Emin = 0.0 
        Emax = 0.0
        DO n=1,nq
           DO i=1, 3*nat
              Emin = MIN (Emin, freq(i,n))
              Emax = MAX (Emax, freq(i,n))
           END DO
        END DO
        !
        ndos = NINT ( (Emax - Emin) / DeltaE+0.500001)  
        OPEN (unit=2,file=fldos,status='unknown',form='formatted')
        DO n= 1, ndos  
           E = Emin + (n - 1) * DeltaE  
           CALL dos_t(freq, 1, 3*nat, nq, ntetra, tetra, E, DOSofE)
           WRITE (2, '(2e12.4)') E, DOSofE (1)
        END DO
        CLOSE(unit=2)
     END IF
     !
  END IF
  ! 
  CALL mp_barrier()
  !
  CALL mp_end()
  !
  STOP
  !
END PROGRAM matdyn
!
!-----------------------------------------------------------------------
SUBROUTINE readfc (flfrc,nr1,nr2,nr3,nrx,frc,epsil,zeu,nat,nax,    &
                   ibrav,alat,tau,at,ntyp,ityp,amass,omega,has_zstar)
  !-----------------------------------------------------------------------
  !
  USE kinds,      ONLY : DP
  !
  IMPLICIT NONE
  ! I/O variable
  CHARACTER(LEN=256) flfrc
  INTEGER ibrav, nr1,nr2,nr3,nrx, nat, nax, ntyp, ityp(nax)
  REAL(KIND=DP) frc(nrx,nrx,nrx,3,3,nax,nax), epsil(3,3),zeu(3,3,nax)
  REAL(KIND=DP) alat, at(3,3), tau(3,nax)
  LOGICAL has_zstar
  ! local variables
  INTEGER i, j, na, nb, m1,m2,m3
  INTEGER ibid, jbid, nabid, nbbid, m1bid,m2bid,m3bid
  REAL(KIND=DP) amass(nax), amass_from_file, celldm(6), omega
  INTEGER nt
  CHARACTER(LEN=3) atm
  !
  !
  OPEN (unit=1,file=flfrc,status='old',form='formatted')
  !
  !
  !  read cell data
  !
  READ(1,*) ntyp,nat,ibrav,(celldm(i),i=1,6)
  IF (nat.GT.nax) CALL errore ('readfc','too many atoms',nat)
  CALL latgen(ibrav,celldm,at(1,1),at(1,2),at(1,3),omega)
  alat = celldm(1)
  at = at / alat !  bring at in units of alat
  CALL volume(alat,at(1,1),at(1,2),at(1,3),omega)
  !
  !  read atomic types, positions and masses
  !
  DO nt = 1,ntyp
     READ(1,*) i,atm,amass_from_file
     IF (i.NE.nt) CALL errore ('readfc','wrong data read',nt)
     IF (amass(nt).EQ.0.d0) THEN
        amass(nt) = amass_from_file
     ELSE
        WRITE(*,*) 'for atomic type',nt,' mass from file not used'
     END IF
  END DO
  DO na=1,nat
     READ(1,*) i,ityp(na),(tau(j,na),j=1,3)
     IF (i.NE.na) CALL errore ('readfc','wrong data read',na)
  END DO
  !
  !  read macroscopic variable
  !
  READ (1,*) has_zstar
  IF (has_zstar) THEN
     READ(1,*) ((epsil(i,j),j=1,3),i=1,3)
     DO na=1,nat
        READ(1,*) 
        READ(1,*) ((zeu(i,j,na),j=1,3),i=1,3)
     END DO
  ELSE
     zeu  (:,:,:) = 0.d0
     epsil(:,:) = 0.d0
  END IF
  !
  !  read real space part
  !
  READ (1,*) nr1,nr2,nr3
  IF (nr1.GT.nrx)  CALL errore ('readin','nr1 .gt. nrx ',+1)
  IF (nr2.GT.nrx)  CALL errore ('readin','nr2 .gt. nrx ',+1)
  IF (nr3.GT.nrx)  CALL errore ('readin','nr3 .gt. nrx ',+1)
  !
  IF(nat.GT.nax) CALL errore  ('readfc','nax too small', nat)
  !
  !  read real-space interatomic force constants
  !
  frc(:,:,:,:,:,:,:) = 0.d0
  DO i=1,3
     DO j=1,3
        DO na=1,nat
           DO nb=1,nat
              READ (1,*) ibid, jbid, nabid, nbbid
              IF(i .NE.ibid  .OR. j .NE.jbid .OR.                   &
                 na.NE.nabid .OR. nb.NE.nbbid)                      &
                 CALL errore  ('readfc','error in reading',1)
              READ (1,*) (((m1bid, m2bid, m3bid,                    &
                          frc(m1,m2,m3,i,j,na,nb),                  &
                           m1=1,nr1),m2=1,nr2),m3=1,nr3)
           END DO
        END DO
     END DO
  END DO
  !
  CLOSE(unit=1)
  !
  RETURN
END SUBROUTINE readfc
!
!-----------------------------------------------------------------------
SUBROUTINE frc_blk(nax,dyn,q,tau,nat,                             &
     &                   nr1,nr2,nr3,nrx,frc,at,bg,rws,nrws)
  !-----------------------------------------------------------------------
  ! calculates the dynamical matrix at q from the (short-range part of the)
  ! force constants 
  !
  USE kinds,      ONLY : DP
  !
  IMPLICIT NONE
  INTEGER nr1, nr2, nr3, nrx, nax, nat, n1, n2, n3, &
          ipol, jpol, na, nb, m1, m2, m3, nint, i,j, nrws
  COMPLEX(KIND=DP) dyn(3,3,nax,nax), cmplx
  REAL(KIND=DP) frc(nrx,nrx,nrx,3,3,nax,nax), tau(3,nax), q(3), arg, &
               at(3,3), bg(3,3), r(3), weight, r_ws(3),  &
               total_weight, rws(0:3,nrws)
  REAL(KIND=DP), PARAMETER:: tpi = 2.0*3.14159265358979d0
  REAL(KIND=DP), EXTERNAL :: wsweight
  !
  DO na=1, nat
     DO nb=1, nat
        total_weight=0.0d0
        DO n1=-2*nrx,2*nrx
           DO n2=-2*nrx,2*nrx
              DO n3=-2*nrx,2*nrx
                 !
                 ! SUM OVER R VECTORS IN THE SUPERCELL - VERY VERY SAFE RANGE!
                 !
                 DO i=1, 3
                    r(i) = n1*at(i,1)+n2*at(i,2)+n3*at(i,3)
                    r_ws(i) = r(i) + tau(i,na)-tau(i,nb)
                 END DO
                 weight = wsweight(r_ws,rws,nrws)
                 IF (weight .GT. 0.0) THEN
                    !
                    ! FIND THE VECTOR CORRESPONDING TO R IN THE ORIGINAL CELL
                    !
                    m1 = MOD(n1+1,nr1)
                    IF(m1.LE.0) m1=m1+nr1
                    m2 = MOD(n2+1,nr2)
                    IF(m2.LE.0) m2=m2+nr2
                    m3 = MOD(n3+1,nr3)
                    IF(m3.LE.0) m3=m3+nr3
                    !
                    ! FOURIER TRANSFORM
                    !
                    arg = tpi*(q(1)*r(1) + q(2)*r(2) + q(3)*r(3))
                    DO ipol=1, 3
                       DO jpol=1, 3
                          dyn(ipol,jpol,na,nb) =                 &
                               dyn(ipol,jpol,na,nb) +            &
                               frc(m1,m2,m3,ipol,jpol,na,nb)     &
                               *CMPLX(COS(arg),-SIN(arg))*weight
                       END DO
                    END DO
                 END IF
                 total_weight=total_weight + weight
              END DO
           END DO
        END DO
        IF (ABS(total_weight-nr1*nr2*nr3).GT.1.0d-8) THEN
           WRITE(*,*) total_weight
           CALL errore ('frc_blk','wrong total_weight',1)
        END IF
     END DO
  END DO
  !
  RETURN
END SUBROUTINE frc_blk
!
!-----------------------------------------------------------------------
SUBROUTINE setupmat (q,dyn,nat,nax,at,bg,tau,itau_blk,nsc,alat, &
     &         dyn_blk,nat_blk,nax_blk,at_blk,bg_blk,tau_blk,omega_blk, &
     &                 epsil,zeu,frc,nr1,nr2,nr3,nrx,has_zstar,rws,nrws)
  !-----------------------------------------------------------------------
  ! compute the dynamical matrix (the analytic part only)
  !
  USE kinds,      ONLY : DP
  !
  IMPLICIT NONE
  REAL(KIND=DP), PARAMETER :: tpi=2.d0*3.14159265358979d0
  !
  ! I/O variables
  !
  INTEGER:: nr1, nr2, nr3, nrx, nax, nat, nat_blk, nax_blk, &
       &    nsc, nrws, itau_blk(nat)
  REAL(KIND=DP) :: q(3), tau(3,nax), at(3,3), bg(3,3), alat,      &
                  epsil(3,3), zeu(3,3,nax_blk), rws(0:3,nrws),   &
                  frc(nrx,nrx,nrx,3,3,nax_blk,nax_blk)
  REAL(KIND=DP) :: tau_blk(3,nax_blk), at_blk(3,3), bg_blk(3,3), omega_blk
  COMPLEX(KIND=DP) dyn_blk(3,3,nax_blk,nax_blk)
  COMPLEX(KIND=DP) ::  dyn(3,3,nax,nax)
  LOGICAL has_zstar
  !
  ! local variables
  !
  REAL(KIND=DP) :: arg
  COMPLEX(KIND=DP) :: cfac(nat)
  INTEGER :: i,j,k, na,nb, na_blk, nb_blk, iq
  REAL(KIND=DP) qp(3), qbid(3,nsc) ! automatic array
  !
  !
  CALL q_gen(nsc,qbid,at_blk,bg_blk,at,bg)
  !
  DO iq=1,nsc
     !
     DO k=1,3
        qp(k)= q(k) + qbid(k,iq)
     END DO
     !
     dyn_blk(:,:,:,:) = (0.d0,0.d0)
     CALL frc_blk (nax_blk,dyn_blk,qp,tau_blk,nat_blk,              &
          &              nr1,nr2,nr3,nrx,frc,at_blk,bg_blk,rws,nrws)
     IF (has_zstar) &
          CALL rgd_blk(nax_blk,nat_blk,dyn_blk,qp,tau_blk,   &
                       epsil,zeu,bg_blk,omega_blk,+1.d0)
     !
     DO na=1,nat
        na_blk = itau_blk(na)
        DO nb=1,nat
           nb_blk = itau_blk(nb)
           !
           arg=tpi* ( qp(1) * ( (tau(1,na)-tau_blk(1,na_blk)) -   &
                                (tau(1,nb)-tau_blk(1,nb_blk)) ) + &
                      qp(2) * ( (tau(2,na)-tau_blk(2,na_blk)) -   &
                                (tau(2,nb)-tau_blk(2,nb_blk)) ) + &
                      qp(3) * ( (tau(3,na)-tau_blk(3,na_blk)) -   &
                                (tau(3,nb)-tau_blk(3,nb_blk)) ) )
           !
           cfac(nb) = CMPLX(COS(arg),SIN(arg))/nsc
           !
        END DO ! nb
        !
        DO i=1,3
           DO j=1,3
              !
              DO nb=1,nat
                 nb_blk = itau_blk(nb)
                 dyn(i,j,na,nb) = dyn(i,j,na,nb) + cfac(nb) * &
                      dyn_blk(i,j,na_blk,nb_blk)
              END DO ! nb
              !
           END DO ! j
        END DO ! i
     END DO ! na
     !
  END DO ! iq
  !
  RETURN
END SUBROUTINE setupmat
!----------------------------------------------------------------------
SUBROUTINE set_asr(nr1,nr2,nr3,nrx,frc,zeu,nat,nax)
  !-----------------------------------------------------------------------
  !
  USE kinds,      ONLY : DP
  !
  IMPLICIT NONE
  INTEGER nr1, nr2, nr3, nrx, nr, i, j, na, nb, n1,n2,n3, nat, nax
  REAL(KIND=DP) frc(nrx,nrx,nrx,3,3,nax,nax), sum, zeu(3,3,nax)
  !
  ! Acoustic Sum Rule on effective charges
  !
  DO i=1,3
     DO j=1,3
        sum=0.0
        DO na=1,nat
           sum = sum + zeu(i,j,na)
        END DO
        DO na=1,nat
           zeu(i,j,na) = zeu(i,j,na) - sum/nat
        END DO
     END DO
  END DO
  !
  ! Acoustic Sum Rule on force constants in real space
  !
  DO i=1,3
     DO j=1,3
        DO na=1,nat
           sum=0.0
           DO nb=1,nat
              DO n1=1,nr1
                 DO n2=1,nr2
                    DO n3=1,nr3
                       sum=sum+frc(n1,n2,n3,i,j,na,nb)
                    END DO
                 END DO
              END DO
           END DO
           frc(1,1,1,i,j,na,na) = frc(1,1,1,i,j,na,na) - sum
           !               write(6,*) ' na, i, j, sum = ',na,i,j,sum
        END DO
     END DO
  END DO
  !
  RETURN
END SUBROUTINE set_asr
!
!-----------------------------------------------------------------------
SUBROUTINE q_gen(nsc,qbid,at_blk,bg_blk,at,bg)
  !-----------------------------------------------------------------------
  ! generate list of q (qbid) that are G-vectors of the supercell
  ! but not of the bulk
  !
  USE kinds,      ONLY : DP
  !
  IMPLICIT NONE
  INTEGER :: nsc
  REAL(KIND=DP) qbid(3,nsc), at_blk(3,3), bg_blk(3,3), at(3,3), bg(3,3)
  !
  INTEGER, PARAMETER:: nr1=4, nr2=4, nr3=4, &
                       nrm=(2*nr1+1)*(2*nr2+1)*(2*nr3+1)
  REAL(KIND=DP), PARAMETER:: eps=1.0e-7
  INTEGER :: i, j, k,i1, i2, i3, idum(nrm), iq
  REAL(KIND=DP) :: qnorm(nrm), qbd(3,nrm) ,qwork(3), delta
  LOGICAL lbho
  !
  i = 0
  DO i1=-nr1,nr1
     DO i2=-nr2,nr2
        DO i3=-nr3,nr3
           i = i + 1
           DO j=1,3
              qwork(j) = i1*bg(j,1) + i2*bg(j,2) + i3*bg(j,3)
           END DO ! j
           !
           qnorm(i)  = qwork(1)**2 + qwork(2)**2 + qwork(3)**2
           !
           DO j=1,3
              !
              qbd(j,i) = at_blk(1,j)*qwork(1) + &
                         at_blk(2,j)*qwork(2) + &
                         at_blk(3,j)*qwork(3)
           END DO ! j
           !
           idum(i) = 1
           !
        END DO ! i3
     END DO ! i2
  END DO ! i1
  !
  DO i=1,nrm-1
     IF (idum(i).EQ.1) THEN
        DO j=i+1,nrm
           IF (idum(j).EQ.1) THEN
              lbho=.TRUE.
              DO k=1,3
                 delta = qbd(k,i)-qbd(k,j)
                 lbho = lbho.AND. (ABS(NINT(delta)-delta).LT.eps)
              END DO ! k
              IF (lbho) THEN
                 IF(qnorm(i).GT.qnorm(j)) THEN
                    qbd(1,i) = qbd(1,j)
                    qbd(2,i) = qbd(2,j)
                    qbd(3,i) = qbd(3,j)
                    qnorm(i) = qnorm(j)
                 END IF
                 idum(j) = 0
              END IF
           END IF
        END DO ! j
     END IF
  END DO ! i
  !
  iq = 0
  DO i=1,nrm
     IF (idum(i).EQ.1) THEN
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
     END IF
  END DO ! i
  !
  IF (iq.NE.nsc) CALL errore('q_gen',' probably nr1,nr2,nr3 too small ', iq)
  RETURN
END SUBROUTINE q_gen
!
!-----------------------------------------------------------------------
SUBROUTINE check_at(at,bg_blk,alat,omega)
  !-----------------------------------------------------------------------
  !
  USE kinds,      ONLY : DP
  !
  IMPLICIT NONE
  !
  REAL(KIND=DP) :: at(3,3), bg_blk(3,3), alat, omega
  REAL(KIND=DP) :: work(3,3)
  INTEGER :: i,j
  REAL(KIND=DP), PARAMETER :: small=1.d-6
  !
  work(:,:) = at(:,:)
  CALL cryst_to_cart(3,work,bg_blk,-1)
  !
  DO j=1,3
     DO i =1,3
        IF ( ABS(work(i,j)-NINT(work(i,j))) > small) THEN
           WRITE (*,'(3f9.4)') work(:,:)
           CALL errore ('check_at','at not multiple of at_blk',1)
        END IF
     END DO
  END DO
  !
  omega =alat**3 * ABS(at(1,1)*(at(2,2)*at(3,3)-at(3,2)*at(2,3))- &
                       at(1,2)*(at(2,1)*at(3,3)-at(2,3)*at(3,1))+ &
                       at(1,3)*(at(2,1)*at(3,2)-at(2,2)*at(3,1)))
  !
  RETURN
END SUBROUTINE check_at
!
!-----------------------------------------------------------------------
SUBROUTINE set_tau                                                &
     &        (nat,nat_blk,at,at_blk,tau,tau_blk,ityp,ityp_blk,itau_blk)
  !-----------------------------------------------------------------------
  !
  USE kinds,      ONLY : DP
  !
  IMPLICIT NONE
  INTEGER nat, nat_blk,ityp(nat),ityp_blk(nat_blk), itau_blk(nat)
  REAL(KIND=DP) at(3,3),at_blk(3,3),tau(3,nat),tau_blk(3,nat_blk)
  !
  REAL(KIND=DP) bg(3,3), r(3) ! work vectors
  INTEGER i,i1,i2,i3,na,na_blk
  REAL(KIND=DP) small
  INTEGER NN1,NN2,NN3
  PARAMETER (NN1=8, NN2=8, NN3=8, small=1.d-8)
  !
  CALL recips (at(1,1),at(1,2),at(1,3),bg(1,1),bg(1,2),bg(1,3))
  !
  na = 0
  !
  DO i1 = -NN1,NN1
     DO i2 = -NN2,NN2
        DO i3 = -NN3,NN3
           r(1) = i1*at_blk(1,1) + i2*at_blk(1,2) + i3*at_blk(1,3)
           r(2) = i1*at_blk(2,1) + i2*at_blk(2,2) + i3*at_blk(2,3)
           r(3) = i1*at_blk(3,1) + i2*at_blk(3,2) + i3*at_blk(3,3)
           CALL cryst_to_cart(1,r,bg,-1)
           !
           IF ( r(1).GT.-small .AND. r(1).LT.1.d0-small .AND.          &
                r(2).GT.-small .AND. r(2).LT.1.d0-small .AND.          &
                r(3).GT.-small .AND. r(3).LT.1.d0-small ) THEN
              CALL cryst_to_cart(1,r,at,+1)
              !
              DO na_blk=1, nat_blk
                 na = na + 1
                 IF (na.GT.nat) CALL errore('set_tau','too many atoms',na)
                 tau(1,na)    = tau_blk(1,na_blk) + r(1)
                 tau(2,na)    = tau_blk(2,na_blk) + r(2)
                 tau(3,na)    = tau_blk(3,na_blk) + r(3)
                 ityp(na)     = ityp_blk(na_blk)
                 itau_blk(na) = na_blk
              END DO
              !
           END IF
           !
        END DO
     END DO
  END DO
  !
  IF (na.NE.nat) CALL errore('set_tau','too few atoms: increase NNs',na)
  !
  RETURN
END SUBROUTINE set_tau
!
!-----------------------------------------------------------------------
SUBROUTINE read_tau(nat,nat_blk,ntyp,bg_blk,tau,tau_blk,ityp,itau_blk)
  !---------------------------------------------------------------------
  !
  USE kinds,      ONLY : DP
  !
  IMPLICIT NONE
  !
  INTEGER nat, nat_blk, ntyp, ityp(nat),itau_blk(nat)
  REAL(KIND=DP) bg_blk(3,3),tau(3,nat),tau_blk(3,nat_blk)
  !
  REAL(KIND=DP) r(3) ! work vectors
  INTEGER i,na,na_blk
  !
  REAL(KIND=DP) small
  PARAMETER ( small = 1.d-6 )
  !
  DO na=1,nat
     READ(*,*) (tau(i,na),i=1,3), ityp(na)
     IF (ityp(na).LE.0 .OR. ityp(na) .GT. ntyp) &
          CALL errore('read_tau',' wrong atomic type', na)
     DO na_blk=1,nat_blk
        r(1) = tau(1,na) - tau_blk(1,na_blk)
        r(2) = tau(2,na) - tau_blk(2,na_blk)
        r(3) = tau(3,na) - tau_blk(3,na_blk)
        CALL cryst_to_cart(1,r,bg_blk,-1)
        IF (ABS( r(1)-NINT(r(1)) ) .LT. small .AND.                 &
            ABS( r(2)-NINT(r(2)) ) .LT. small .AND.                 &
            ABS( r(3)-NINT(r(3)) ) .LT. small ) THEN
           itau_blk(na) = na_blk
           go to 999
        END IF
     END DO
     CALL errore ('read_tau',' wrong atomic position ', na)
999  CONTINUE
  END DO
  !
  RETURN
END SUBROUTINE read_tau
!
!-----------------------------------------------------------------------
SUBROUTINE write_tau(fltau,nat,tau,ityp)
  !-----------------------------------------------------------------------
  !
  USE kinds,      ONLY : DP
  !
  IMPLICIT NONE
  !
  INTEGER nat, ityp(nat)
  REAL(KIND=DP) tau(3,nat)
  CHARACTER(LEN=*) fltau
  !
  INTEGER i,na
  !
  OPEN (unit=4,file=fltau, status='new')
  DO na=1,nat
     WRITE(4,'(3(f12.6),i3)') (tau(i,na),i=1,3), ityp(na)
  END DO
  CLOSE (4)
  !
  RETURN 
END SUBROUTINE write_tau
!
!-----------------------------------------------------------------------
SUBROUTINE gen_qpoints (ibrav, at, bg, nat, tau, ityp, nk1, nk2, nk3, &
     ntetra, nqx, nq, q, tetra)
  !-----------------------------------------------------------------------
  !
  USE kinds,      ONLY : DP
  !
  IMPLICIT NONE
  ! input 
  INTEGER :: ibrav, nat, nqx, nk1, nk2, nk3, ntetra, ityp(*)
  REAL(KIND=DP) :: at(3,3), bg(3,3), tau(3,nat)
  ! output 
  INTEGER :: nq, tetra(4,ntetra)
  REAL(KIND=DP) :: q(3,nqx)
  ! local
  INTEGER :: nrot, nsym, s(3,3,48), ftau(3,48), irt(48,nat)
  LOGICAL :: minus_q, invsym
  REAL(KIND=DP) :: xqq(3), wk(nqx), mdum(3,nat)
  CHARACTER(LEN=45)   ::  sname(48)
  !
  xqq (:) =0.d0
  IF (ibrav == 4 .OR. ibrav == 5) THEN  
     !
     !  hexagonal or trigonal bravais lattice
     !
     CALL hexsym (at, s, sname, nrot)  
  ELSEIF (ibrav >= 1 .AND. ibrav <= 14) THEN  
     !
     !  cubic bravais lattice
     !
     CALL cubicsym (at, s, sname, nrot)  
  ELSEIF (ibrav == 0) THEN  
     CALL errore ('gen_qpoints', 'assuming cubic symmetry',-1)  
     CALL cubicsym (at, s, sname, nrot)  
  ELSE  
     CALL errore ('gen_qpoints', 'wrong ibrav', 1)  
  ENDIF
  !
  CALL kpoint_grid ( nrot, s, bg, nqx, 0,0,0, nk1,nk2,nk3, nq, q, wk)
  !
  CALL sgama (nrot, nat, s, sname, at, bg, tau, ityp, nsym, 6, &
       6, 6, irt, ftau, nqx, nq, q, wk, invsym, minus_q, xqq, &
       0, 0, .FALSE., mdum)
  
  IF (ntetra /= 6 * nk1 * nk2 * nk3) &
       CALL errore ('gen_qpoints','inconsistent ntetra',1)

  CALL tetrahedra (nsym, s, minus_q, at, bg, nqx, 0, 0, 0, &
       nk1, nk2, nk3, nq, q, wk, ntetra, tetra)
  !
  RETURN
END SUBROUTINE gen_qpoints
