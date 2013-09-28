!
! Copyright (C) 2001-2008 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
PROGRAM q2trans 
  !----------------------------------------------------------------------------
  !
  !  q2r.x:
  !     reads force constant matrices C(q) produced by the phonon code
  !     for a grid of q-points, calculates the corresponding set of
  !     interatomic force constants (IFC), C(R)
  !
  !  Input data: Namelist "input"
  !     fildyn     :  input file name (character, must be specified)
  !                   "fildyn"0 contains information on the q-point grid
  !                   "fildyn"1-N contain force constants C_n = C(q_n)
  !                   for n=1,...N, where N is the number of q-points
  !                   in the irreducible brillouin zone
  !                   Normally this should be the same as specified
  !                   on input to the phonon code
  !                   In the non collinear/spin-orbit case the files 
  !                   produced by ph.x are in .xml format. In this case 
  !                   fildyn is the same as in the phonon code + the .xml 
  !                   extension.
  !     flfrc      :  output file containing the IFC in real space
  !                   (character, must be specified)
  !     zasr       :  Indicates type of Acoustic Sum Rules used for the Born
  !                   effective charges (character):
  !                   - 'no': no Acoustic Sum Rules imposed (default)
  !                   - 'simple':  previous implementation of the asr used
  !                     (3 translational asr imposed by correction of
  !                     the diagonal elements of the force-constants matrix)
  !                   - 'crystal': 3 translational asr imposed by optimized
  !                      correction of the IFC (projection).
  !                   - 'one-dim': 3 translational asr + 1 rotational asr
  !                     imposed by optimized correction of the IFC (the
  !                     rotation axis is the direction of periodicity; it
  !                     will work only if this axis considered is one of
  !                     the cartesian axis).
  !                   - 'zero-dim': 3 translational asr + 3 rotational asr
  !                     imposed by optimized correction of the IFC.
  !                   Note that in certain cases, not all the rotational asr
  !                   can be applied (e.g. if there are only 2 atoms in a
  !                   molecule or if all the atoms are aligned, etc.).
  !                   In these cases the supplementary asr are cancelled
  !                   during the orthonormalization procedure (see below).
  !
  !  If a file "fildyn"0 is not found, the code will ignore variable "fildyn"
  !  and will try to read from the following cards the missing information
  !  on the q-point grid and file names:
  !     nr1,nr2,nr3:  dimensions of the FFT grid formed by the q-point grid
  !     nfile      :  number of files containing C(q_n), n=1,nfile
  !  followed by nfile cards:
  !     filin      :  name of file containing C(q_n)
  !  The name and order of files is not important as long as q=0 is the first
  !
USE iotk_module  
USE kinds,      ONLY : DP
  USE mp,         ONLY : mp_bcast
  USE mp_global,  ONLY : mp_startup, mp_global_end
  USE mp_world,   ONLY : nproc, mpime, world_comm
  USE dynamicalq, ONLY : phiq, tau, ityp, zeu
  USE fft_scalar, ONLY : cfft3d
  USE io_global, ONLY : ionode_id, ionode, stdout
  USE io_dyn_mat, ONLY : read_dyn_mat_param, read_dyn_mat_header, &
                         read_dyn_mat, read_dyn_mat_tail, &
                         write_dyn_mat_header, write_ifc
  USE environment, ONLY : environment_start, environment_end
use constants, only: pi, fpi, e2
  !
  IMPLICIT NONE
  !
  INTEGER,       PARAMETER :: ntypx = 10
  REAL(DP), PARAMETER :: eps=1.D-5, eps12=1.d-12
  INTEGER                  :: nr1, nr2, nr3, nr(3)
  !     dimensions of the FFT grid formed by the q-point grid
  !
  CHARACTER(len=20)  :: crystal
  CHARACTER(len=256) :: fildyn, filin, filj, filf, flfrc, file_ifc
  CHARACTER(len=3)   :: atm(ntypx)
  CHARACTER(LEN=6), EXTERNAL :: int_to_char
  !
  LOGICAL :: lq, lrigid, lrigid1, lnogridinfo, xmldyn, ltrans
  CHARACTER (LEN=10) :: zasr, iasr
  INTEGER :: m1, m2, m3, m(3), l1, l2, l3, j1, j2, na1, na2, ipol, nn
  INTEGER :: nat, nq, ntyp, iq, icar, nfile, ifile, nqs, nq_log
  INTEGER :: na, nt, n1, n2, n3, nrx, ndummy, nb
  !
  INTEGER :: gid, ibrav, ierr, nspin_mag, ios, idir
  !
  INTEGER, ALLOCATABLE ::  nc(:,:,:)
  COMPLEX(DP), ALLOCATABLE :: phid(:,:,:,:,:)
  REAL(DP), ALLOCATABLE :: ifc3(:,:,:,:,:,:,:), ifc(:,:,:,:,:), ifc0(:,:,:,:), frc(:,:,:,:,:,:,:),  kfc(:,:,:), k00(:,:), k01(:,:)
  REAL(DP), ALLOCATABLE :: m_loc(:,:)
  !
  REAL(DP) :: celldm(6), at(3,3), bg(3,3)
  REAL(DP) :: q(3,48), omega, xq, amass(ntypx), resi, sum1, sum2
  REAL(DP) :: epsil(3,3), d1(3), dd1, d2(3), dd2

  REAL(DP) :: amconv = 1.66042d-24/9.1095d-28*0.5d0 !12.0107
  !
  logical           :: la2F, onedim
  LOGICAL, EXTERNAL :: has_xml

INTEGER  :: dimwan
INTEGER  :: nkpts
INTEGER  :: nrtot

CHARACTER(256)       :: fileout 

INTEGER  :: i, j, ik, ir, nsize

LOGICAL   :: have_overlap, htype, noNA, readifc, dielec
REAL ::  fermi_energy

INTEGER, ALLOCATABLE :: nk(:), ivr(:,:)
REAL,    ALLOCATABLE :: wr(:)
COMPLEX, ALLOCATABLE :: rham(:,:,:), ovp(:,:,:)
REAL, ALLOCATABLE :: r_rham(:,:,:), r_ovp(:,:,:)

CHARACTER(600)       :: attr, card
INTEGER, PARAMETER         ::   &
    stdin = 5

  !
  NAMELIST / input / fildyn, flfrc, zasr, la2F, onedim, noNA, idir, fileout, readifc, file_ifc, nr1,nr2,nr3, nat, amass
  !
  CALL mp_startup()
  CALL environment_start('Q2R')
  !
  IF (ionode) CALL input_from_file ( )
     !
  fildyn = ' '
  flfrc = ' '
  zasr = 'no'
  ltrans = .true.
  onedim=.false.
  noNA=.true.
  idir=1
!  ifc=.false.
     !
  la2F=.false.
     !
     !
  IF (ionode)  READ ( 5, input, IOSTAT =ios )
 
  CALL mp_bcast(ios, ionode_id, world_comm )
  CALL errore('q2r','error reading input namelist', abs(ios))

     !
     ! Define IFCs for transport calculation (MBN, April 2009)
     !
     ALLOCATE (ifc3(3,3,nat,nat,nr1,nr2,nr3) )
     ALLOCATE (frc(nr1,nr2,nr3,3,3,nat,nat) )

allo_dir: SELECT CASE (idir)
     CASE(1)
     ALLOCATE (ifc(3,3,nat,nat,nr1) )
     ALLOCATE (kfc(3*nat,3*nat,nr1/2+1) )
     ALLOCATE (k00(3*nat*(nr1/2+1),3*nat*(nr1/2+1) ), k01(3*nat*(nr1/2+1),3*nat*(nr1/2+1) ) )
     CASE(2)
     ALLOCATE (ifc(3,3,nat,nat,nr2) )
     ALLOCATE (kfc(3*nat,3*nat,nr2/2+1) )
     ALLOCATE (k00(3*nat*(nr2/2+1),3*nat*(nr2/2+1) ), k01(3*nat*(nr2/2+1),3*nat*(nr2/2+1) ) )
     CASE(3)
     ALLOCATE (ifc(3,3,nat,nat,nr3) )
     ALLOCATE (kfc(3*nat,3*nat,nr3/2+1) )
     ALLOCATE (k00(3*nat*(nr3/2+1),3*nat*(nr3/2+1) ), k01(3*nat*(nr3/2+1),3*nat*(nr3/2+1) ) )
END SELECT allo_dir

     ALLOCATE (ifc0(3,3,nat,nat) )

     ifc(:,:,:,:,:)=0.0
     ifc0(:,:,:,:)=0.0
     frc(:,:,:,:,:,:,:)=0.0
     ifc3(:,:,:,:,:,:,:)=0.0

     allocate (ityp(nat))

  ! Read the IFCs from file

filin=TRIM(file_ifc)//'.fc'
OPEN(3,FILE=TRIM(filin))

READ(3,*) ntyp, nat
do na=1,ntyp
   READ(3,*) 
end do
do na=1,nat
   READ (3,*) ndummy, ityp(na)
end do

READ(3,*) dielec
IF(dielec) THEN
do i=1,3
   READ(3,*)
end do
do na=1,nat
   READ(3,*)
   do i=1,3
      READ(3,*)
   end do
end do
END IF

READ (3,'(3i4)') ndummy,ndummy,ndummy

do i=1,3
   do j=1,3
      do na=1,nat
         do nb=1,nat
            READ (3,*)
            do m1=1,nr1
              do m2=1,nr2
                 do m3=1,nr3
                    READ (3,'(3i4,2x,1pe18.11)')ndummy,ndummy,ndummy,ifc3(i,j,na,nb,m1,m2,m3)
                 end do
               end do
             end do
         end do
      end do
   end do
end do

CLOSE (3)

  ! Construct the IFC matrix with the correct symmetry for transport calculations

direction: SELECT CASE (idir)

CASE(1)
  DO j1=1,3
     DO j2=1,3
        DO na1=1,nat
           DO na2=1,nat
              DO m1=1,nr1
                 DO m2=1,nr2
                       DO m3=1,nr3
                       ! for transport in one-dim systems
                       IF(onedim) THEN
                          IF(m2.eq.1.and.m3.eq.1) ifc(j1,j2,na1,na2,m1)=ifc3(j1,j2,na1,na2,m1,m2,m3)
                       ENDIF
                       ! for transport in 3-dim systems: sum on the plane perpendicular to the transport direction
                        ifc(j1,j2,na1,na2,m1)=ifc(j1,j2,na1,na2,m1)+ifc3(j1,j2,na1,na2,m1,m2,m3)
                    END DO
                 END DO
              END DO
           END DO
        END DO
     END DO
  END DO
  nrx=nr1

CASE(2)
  DO j1=1,3
     DO j2=1,3
        DO na1=1,nat
           DO na2=1,nat
              DO m2=1,nr2
                 DO m1=1,nr1
                    DO m3=1,nr3
                       ! for transport in one-dim systems
                       IF(onedim) THEN
                          IF(m1.eq.1.and.m3.eq.1) ifc(j1,j2,na1,na2,m2)=ifc3(j1,j2,na1,na2,m1,m2,m3)
                       ENDIF
                       ! for transport in 3-dim systems: sum on the plane perpendicular to the transport direction
                        ifc(j1,j2,na1,na2,m2)=ifc(j1,j2,na1,na2,m2)+ifc3(j1,j2,na1,na2,m1,m2,m3)
                    END DO
                 END DO
              END DO
           END DO
        END DO
     END DO
  END DO
  nrx=nr2

CASE(3)
  DO j1=1,3
     DO j2=1,3
        DO na1=1,nat
           DO na2=1,nat
              DO m3=1,nr3
                 DO m2=1,nr2
                    DO m1=1,nr1
                       ! for transport in one-dim systems
                       IF(onedim) THEN
                          IF(m1.eq.1.and.m2.eq.1) ifc(j1,j2,na1,na2,m3)=ifc3(j1,j2,na1,na2,m1,m2,m3)
                       ENDIF
                       ! for transport in 3-dim systems: sum on the plane perpendicular to the transport direction
                        ifc(j1,j2,na1,na2,m3)=ifc(j1,j2,na1,na2,m3)+ifc3(j1,j2,na1,na2,m1,m2,m3)
                    END DO
                 END DO
              END DO
           END DO
        END DO
     END DO
  END DO
  nrx=nr3

  END SELECT direction

  ! Correction for finite IFC in the center of the real space mesh 
  IF(nrx > 1) THEN
  DO j1=1,3
     DO j2=1,3
        DO na1=1,nat
           DO na2=1,nat
              ifc0(j1,j2,na1,na2)= ifc(j1,j2,na1,na2,nrx/2+1)
           END DO
        END DO
     END DO
  END DO
  
  DO j1=1,3
     DO j2=1,3
        DO na1=1,nat
           DO na2=1,nat
              DO m1=1,nrx
                 ifc(j1,j2,na1,na2,m1)= ifc(j1,j2,na1,na2,m1)-ifc0(j1,j2,na1,na2)
              END DO
           END DO
        END DO
     END DO
  END DO
  ENDIF
  
  ! Impose the acoustic sum rule for the shifted IFC: the interatomic force of the atom on itself should be 
  ! equal to minus the sum of all interatomic forces generated by all others atoms (action-reaction law!)
  ! eq. (82) in Gonze and Lee, PRB 55, 10355 (1997)
  
  DO j1=1,3
     DO j2=1,3
        DO na1=1,nat
           sum1=0.0
           DO na2=1,nat
              IF(na1.ne.na2) sum1=sum1+ifc(j1,j2,na1,na2,1)
           END DO
           sum2=0.0
           DO na2=1,nat
              DO m1=2,nrx
                 sum2=sum2+ifc(j1,j2,na1,na2,m1)
              END DO
           END DO
        END DO
     END DO
  END DO

  ! Check the range of the IFC in the slab

!  DO j1=1,3
!     DO j2=1,3
!        DO na1=1,nat
!           DO m1=1,nrx
!              WRITE(*,'(4I3,1x,1F12.6)') na1, j1, j2, m1, ifc(j1,j2,1,na1,m1)
!           ENDDO
!        ENDDO
!     ENDDO
!  ENDDO
  
  ! Write the IFC for heat transport.
  
  DO m1=1,nrx/2+1
     DO na1=1,nat
        DO na2=1,nat
           DO j1=1,3
              DO j2=1,3
                 kfc(3*(na1-1)+j1,3*(na2-1)+j2,m1) = ifc(j1,j2,na1,na2,m1)
              END DO
           END DO
        END DO
     END DO
  END DO
  
  ! define k00
  
  DO i=1,3*nat*(nrx/2+1)
     DO j=1,3*nat*(nrx/2+1)
        k00(i,j)=0.0d0
     END DO
  END DO
  
  DO m1=1,nrx/2+1
     DO m2=m1,nrx/2+1
        DO j1=1,3*nat
           DO j2=1,3*nat
              k00(3*nat*(m1-1)+j1,3*nat*(m2-1)+j2 ) = &
                   amconv/SQRT(amass(ityp((j1-1)/3+1))*amass(ityp((j2-1)/3+1)))*kfc(j1,j2,m2-m1+1)
           END DO
        END DO
     END DO
  END DO
  
  DO i=1,3*nat*(nrx/2+1)
     DO j=1,i-1
        k00(i,j)=k00(j,i)
     END DO
  END DO
  
  ! define k01
  
  DO i=1,3*nat*(nrx/2+1)
     DO j=1,3*nat*(nrx/2+1)
        k01(i,j)=0.0d0
     END DO
  END DO
  
  DO m1=1,nrx/2+1
     DO m2=1,m1-1
        DO j1=1,3*nat
           DO j2=1,3*nat
              k01(3*nat*(m1-1)+j1,3*nat*(m2-1)+j2 ) = &
                   amconv/SQRT(amass(ityp((j1-1)/3+1))*amass(ityp((j2-1)/3+1)))*kfc(j1,j2,m2+(nrx/2+1)-m1+1)
           END DO
        END DO
     END DO
  END DO

! 
! write to file 
! 
nrtot=2
dimwan=3*nat*(nrx/2+1)
nkpts=1
ALLOCATE (nk(3))

nk(:)  = 1
nr(:)  = 0

ALLOCATE( wr(nkpts) )
wr(:)  = 1

ALLOCATE( ivr(3,nrtot) )
ALLOCATE( rham(dimwan,dimwan,nrtot) )
ALLOCATE( ovp(dimwan,dimwan,nrtot) )

! extract submatrices for leads

rham(:,:,1)=cmplx(K00(:,:),0.0)
rham(:,:,2)=(0.0,0.0)

! 00
do i=1,36
   write(77,'(36(1x,2f12.8))')(rham(i,j,1),j=1,36)
end do

! 01 - 13
do i=13,24
   do j=1,12
      rham(i,j,2)=rham(i-12,j+24,1)
   end do
end do

! 01 - 13
do i=25,36
   do j=13,24
      rham(i,j,2)=rham(i-24,j+12,1)
   end do
end do

! 01 - 12
do i=25,36
   do j=1,12
      rham(i,j,2)=rham(i-24,j+12,1)
   end do
end do

do i=1,36
   write(78,'(24(1x,2f12.8))')(rham(i,j,2),j=1,36)
end do

! end extracting matrices

rham(:,:,1)=cmplx(K00(:,:),0.0)
rham(:,:,2)=cmplx(K01(:,:),0.0)

vectors: SELECT CASE (idir)

CASE(1)
    nr(1)=2
    nr(2)=1
    nr(3)=1
    ivr(1,1)=0
    ivr(2,1)=0
    ivr(3,1)=0
    ivr(1,2)=1
    ivr(2,2)=0
    ivr(3,2)=0

CASE(2)
    nr(1)=1
    nr(2)=2
    nr(3)=1
    ivr(1,1)=0
    ivr(2,1)=0
    ivr(3,1)=0
    ivr(1,2)=0
    ivr(2,2)=1
    ivr(3,2)=0

CASE(3)
    nr(1)=1
    nr(2)=1
    nr(3)=2
    ivr(1,1)=0
    ivr(2,1)=0
    ivr(3,1)=0
    ivr(1,2)=0
    ivr(2,2)=0
    ivr(3,2)=1

END SELECT vectors

fermi_energy=0.0

have_overlap  = .false.

CALL iotk_open_write( stdout, FILE=TRIM(fileout))

CALL iotk_write_begin(stdout,"HAMILTONIAN")

CALL iotk_write_attr( attr, "dimwann", dimwan, FIRST=.TRUE. )
CALL iotk_write_attr( attr, "nkpts", nkpts )
CALL iotk_write_attr( attr, "nk", nk )
CALL iotk_write_attr( attr, "nrtot", nrtot )
CALL iotk_write_attr( attr, "nr", nr )
CALL iotk_write_attr( attr, "have_overlap", have_overlap )
CALL iotk_write_attr( attr, "fermi_energy", fermi_energy )

CALL iotk_write_empty( stdout, "DATA", ATTR=attr)

nsize=3*2
CALL iotk_write_attr( attr, "type", "integer", FIRST=.TRUE. )
CALL iotk_write_attr( attr, "size", nsize )
CALL iotk_write_attr( attr, "columns", 3 )
CALL iotk_write_attr( attr, "units", "crystal" )
CALL iotk_write_dat( stdout, "IVR", ivr, COLUMNS=3, ATTR=attr )

CALL iotk_write_attr( attr, "type", "real", FIRST=.TRUE. )
CALL iotk_write_attr( attr, "size", nkpts )
CALL iotk_write_dat( stdout, "WR", wr, ATTR=attr )

CALL iotk_write_begin(stdout,"RHAM")
DO ir = 1, nrtot
    CALL iotk_write_dat(stdout,"VR"//TRIM(iotk_index(ir)), rham(:,:,ir))
ENDDO
CALL iotk_write_end(stdout,"RHAM")
CALL iotk_write_end(stdout,"HAMILTONIAN")

CALL iotk_close_write( stdout )


  !
  CALL environment_end('Q2R')

  CALL mp_global_end()
  !
END PROGRAM q2trans 
!
!----------------------------------------------------------------------------
SUBROUTINE gammaq2r( nqtot, nat, nr1, nr2, nr3, at )
  !----------------------------------------------------------------------------
  !
  USE kinds, ONLY : DP
  USE fft_scalar, ONLY : cfft3d
  USE io_global, ONLY : ionode, ionode_id, stdout
  USE mp,        ONLY : mp_bcast
  USE mp_world,  ONLY : world_comm
  !
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: nqtot, nat, nr1, nr2, nr3
  REAL(DP), INTENT(IN) :: at(3,3)
  !
  INTEGER, ALLOCATABLE :: nc(:,:,:)
  COMPLEX(DP), ALLOCATABLE :: gaminp(:,:,:,:,:), gamout(:,:,:,:,:)
  !
  REAL(DP), PARAMETER :: eps=1.D-5, eps12=1.d-12
  INTEGER  :: nsig = 10, isig, filea2F, nstar, count_q, nq, nq_log, iq, &
       icar, ipol, m1,m2,m3, m(3), nr(3), j1,j2, na1, na2, nn
  LOGICAL :: lq
  REAL(DP) :: deg, ef, dosscf
  REAL(DP) :: q(3,48), xq, resi
  character(len=14) :: name

  !
  ALLOCATE (gaminp(3,3,nat,nat,48), gamout(nr1*nr2*nr3,3,3,nat,nat) )
  ALLOCATE ( nc (nr1,nr2,nr3) )
  write (stdout,*)
  write (stdout,*) '  Preparing gamma for a2F '
  write (stdout,*)
  !
  nr(1) = nr1
  nr(2) = nr2
  nr(3) = nr3
  !
  DO isig=1, nsig
     filea2F = 50 + isig
     write(name,"(A7,I2)") 'a2Fq2r.',filea2F
     IF (ionode) open(filea2F, file=name, STATUS = 'old', FORM = 'formatted')
     nc = 0
     !
     ! to pass to matdyn, for each isig, we read: degauss, Fermi energy and DOS
     !
     DO count_q=1,nqtot
        !
        IF (ionode) THEN
           READ(filea2F,*) deg, ef, dosscf
           READ(filea2F,*) nstar
        ENDIF
        CALL mp_bcast(deg, ionode_id, world_comm )
        CALL mp_bcast(ef, ionode_id, world_comm )
        CALL mp_bcast(dosscf, ionode_id, world_comm )
        CALL mp_bcast(nstar, ionode_id, world_comm )
        !
        CALL read_gamma ( nstar, nat, filea2F, q, gaminp )
        !
        do nq = 1,nstar
           lq = .true.
           do ipol=1,3
              xq = 0.0d0
              do icar=1,3
                 xq = xq + at(icar,ipol) * q(icar,nq) * nr(ipol)
              end do
              lq = lq .AND. (ABS(NINT(xq) - xq) < eps)
              iq = NINT(xq)
              !
              m(ipol)= mod(iq,nr(ipol)) + 1
              if (m(ipol) < 1) m(ipol) = m(ipol) + nr(ipol)
           end do !ipol
           IF (.NOT.lq) CALL errore('init','q not allowed',1)
           !
           if(nc(m(1),m(2),m(3)) == 0) then
              nc(m(1),m(2),m(3)) = 1
              CALL TRASL( gamout, gaminp, nq, nr1, nr2, nr3, nat, m(1), m(2), m(3) )
           else
              call errore('init',' nc already filled: wrong q grid or wrong nr',1)
           end if
        enddo ! stars for given q-point
     ENDDO ! q-points
     !
     nq_log = SUM (nc)
     if (nq_log == nr1*nr2*nr3) then
        write (stdout,*)
        write (stdout,'(" Broadening = ",F10.3)') deg
        write (stdout,'(5x,a,i4)') ' q-space grid ok, #points = ',nq_log
     else
        call errore('init',' missing q-point(s)!',1)
     end if
     do j1=1,3
        do j2=1,3
           do na1=1,nat
              do na2=1,nat
                 call cfft3d ( gamout(:,j1,j2,na1,na2), &
                      nr1,nr2,nr3, nr1,nr2,nr3, 1 )
              end do
           end do
        end do
     end do
     gamout = gamout / DBLE (nr1*nr2*nr3)
     !
     IF (ionode) close(filea2F)
     !
     filea2F = 60 + isig
     write(name,"(A10,I2)") 'a2Fmatdyn.',filea2F
     IF (ionode) THEN
     open(filea2F, file=name, STATUS = 'unknown')
     !
     WRITE(filea2F,*) deg, ef, dosscf
     write(filea2F,'(3i4)') nr1, nr2, nr3

     do j1=1,3
        do j2=1,3
           do na1=1,nat
              do na2=1,nat
                 write(filea2F,'(4i4)') j1,j2,na1,na2
                 nn=0
                 DO m3=1,nr3
                    DO m2=1,nr2
                       DO m1=1,nr1
                          nn=nn+1
                          write(filea2F,'(3i4,2x,1pe18.11)')   &
                               m1,m2,m3, DBLE(gamout(nn,j1,j2,na1,na2))
                       END DO
                    END DO
                 END DO
              end do  ! na2
           end do  ! na1
        end do   !  j2
     end do   ! j1
     close(filea2F)
     ENDIF  ! ionode

     resi = SUM ( ABS ( AIMAG( gamout ) ) )

     IF (resi > eps12) THEN
        WRITE (stdout,"(/5x,' fft-check warning: sum of imaginary terms = ',e12.7)") resi
     ELSE
        WRITE (stdout,"(/5x,' fft-check success (sum of imaginary terms < 10^-12)')")
     END IF

  ENDDO
  !
  DEALLOCATE (gaminp, gamout )
  !
END SUBROUTINE gammaq2r
!
!-----------------------------------------------------------------------
subroutine read_gamma (nqs, nat, ifn, xq, gaminp)
  !-----------------------------------------------------------------------
  !
  USE kinds, ONLY : DP
  USE io_global, ONLY : ionode, ionode_id, stdout
  USE mp,        ONLY : mp_bcast
  USE mp_world,  ONLY : world_comm
  implicit none
  !
  ! I/O variables
  integer, intent(in) :: nqs, nat, ifn
  real(DP), intent(out) :: xq(3,48)
  complex(DP), intent(out) :: gaminp(3,3,nat,nat,48)
  !
  logical :: lrigid
  integer :: i, j, na, nb, nt, iq
  real(DP) :: phir(3),phii(3)
  CHARACTER(LEN=75) :: line
  !
  !
  Do iq=1,nqs
     IF (ionode) THEN
        READ(ifn,*)
        READ(ifn,*)
        READ(ifn,*)
        READ(ifn,'(11X,3F14.9)')  (xq(i,iq),i=1,3)
     !     write(*,*) 'xq    ',iq,(xq(i,iq),i=1,3)
        READ(ifn,*)
     END IF
     CALL mp_bcast(xq(:,iq), ionode_id, world_comm)
     do na=1,nat
        do nb=1,nat
           IF (ionode) read(ifn,*) i,j
           CALL mp_bcast(i, ionode_id, world_comm)
           CALL mp_bcast(j, ionode_id, world_comm)
           if (i.ne.na) call errore('read_gamma','wrong na read',na)
           if (j.ne.nb) call errore('read_gamma','wrong nb read',nb)
           do i=1,3
              IF (ionode) read (ifn,*) (phir(j),phii(j),j=1,3)
              CALL mp_bcast(phir, ionode_id, world_comm)
              CALL mp_bcast(phii, ionode_id, world_comm)
              do j = 1,3
                 gaminp(i,j,na,nb,iq) = CMPLX(phir(j),phii(j),kind=DP)
              end do
              !           write(*,*) 'gaminp  ',(gaminp(i,j,na,nb,iq),j=1,3)
           end do
        end do
     end do
     !
  ENDDO
  RETURN
  !
end subroutine read_gamma
!
!----------------------------------------------------------------------------
SUBROUTINE trasl( phid, phiq, nq, nr1, nr2, nr3, nat, m1, m2, m3 )
  !----------------------------------------------------------------------------
  !
  USE kinds, ONLY : DP
  !
  IMPLICIT NONE
  INTEGER, intent(in) ::  nr1, nr2, nr3, m1, m2, m3, nat, nq
  COMPLEX(DP), intent(in) :: phiq(3,3,nat,nat,48)
  COMPLEX(DP), intent(out) :: phid(nr1,nr2,nr3,3,3,nat,nat)
  !
  INTEGER :: j1,j2,  na1, na2
  !
  DO j1=1,3
     DO j2=1,3
        DO na1=1,nat
           DO na2=1,nat
              phid(m1,m2,m3,j1,j2,na1,na2) = &
                   0.5d0 * (      phiq(j1,j2,na1,na2,nq) +  &
                          CONJG(phiq(j2,j1,na2,na1,nq)))
           END DO
        END DO
     END DO
  END DO
  !
  RETURN
END SUBROUTINE trasl
!----------------------------------------------------------------------
subroutine set_zasr ( zasr, nr1,nr2,nr3, nat, ibrav, tau, zeu)
  !-----------------------------------------------------------------------
  !
  ! Impose ASR - refined version by Nicolas Mounet
  !
  USE kinds, ONLY : DP
  USE io_global, ONLY : stdout
  implicit none
  character(len=10) :: zasr
  integer ibrav,nr1,nr2,nr3,nr,m,p,k,l,q,r
  integer n,i,j,n1,n2,n3,na,nb,nat,axis,i1,j1,na1
  !
  real(DP) sum, zeu(3,3,nat)
  real(DP) tau(3,nat), zeu_new(3,3,nat)
  !
  real(DP) zeu_u(6*3,3,3,nat)
  ! These are the "vectors" associated with the sum rules on effective charges
  !
  integer zeu_less(6*3),nzeu_less,izeu_less
  ! indices of vectors zeu_u that are not independent to the preceding ones,
  ! nzeu_less = number of such vectors, izeu_less = temporary parameter
  !
  real(DP) zeu_w(3,3,nat), zeu_x(3,3,nat),scal,norm2
  ! temporary vectors and parameters

  ! Initialization.
  ! n is the number of sum rules to be considered (if zasr.ne.'simple')
  ! and 'axis' is the rotation axis in the case of a 1D system
  ! (i.e. the rotation axis is (Ox) if axis='1', (Oy) if axis='2'
  ! and (Oz) if axis='3')
  !
  if((zasr.ne.'simple').and.(zasr.ne.'crystal').and.(zasr.ne.'one-dim') &
                       .and.(zasr.ne.'zero-dim')) then
      call errore('q2r','invalid Acoustic Sum Rulei for Z*:' // zasr, 1)
  endif
  if(zasr.eq.'crystal') n=3
  if(zasr.eq.'one-dim') then
     ! the direction of periodicity is the rotation axis
     ! It will work only if the crystal axis considered is one of
     ! the cartesian axis (typically, ibrav=1, 6 or 8, or 4 along the
     ! z-direction)
     if (nr1*nr2*nr3.eq.1) axis=3
     if ((nr1.ne.1).and.(nr2*nr3.eq.1)) axis=1
     if ((nr2.ne.1).and.(nr1*nr3.eq.1)) axis=2
     if ((nr3.ne.1).and.(nr1*nr2.eq.1)) axis=3
     if (((nr1.ne.1).and.(nr2.ne.1)).or.((nr2.ne.1).and. &
          (nr3.ne.1)).or.((nr1.ne.1).and.(nr3.ne.1))) then
        call errore('q2r','too many directions of &
             &   periodicity in 1D system',axis)
     endif
     if ((ibrav.ne.1).and.(ibrav.ne.6).and.(ibrav.ne.8).and. &
          ((ibrav.ne.4).or.(axis.ne.3)) ) then
        write(stdout,*) 'zasr: rotational axis may be wrong'
     endif
     write(stdout,'("zasr rotation axis in 1D system= ",I4)') axis
     n=4
  endif
  if(zasr.eq.'zero-dim') n=6

  ! Acoustic Sum Rule on effective charges
  !
  if(zasr.eq.'simple') then
     do i=1,3
        do j=1,3
           sum=0.0d0
           do na=1,nat
               sum = sum + zeu(i,j,na)
            end do
            do na=1,nat
               zeu(i,j,na) = zeu(i,j,na) - sum/nat
            end do
         end do
      end do
   else
      ! generating the vectors of the orthogonal of the subspace to project
      ! the effective charges matrix on
      !
      zeu_u(:,:,:,:)=0.0d0
      do i=1,3
         do j=1,3
            do na=1,nat
               zeu_new(i,j,na)=zeu(i,j,na)
            enddo
         enddo
      enddo
      !
      p=0
      do i=1,3
         do j=1,3
            ! These are the 3*3 vectors associated with the
            ! translational acoustic sum rules
            p=p+1
            zeu_u(p,i,j,:)=1.0d0
            !
         enddo
      enddo
      !
      if (n.eq.4) then
         do i=1,3
            ! These are the 3 vectors associated with the
            ! single rotational sum rule (1D system)
            p=p+1
            do na=1,nat
               zeu_u(p,i,MOD(axis,3)+1,na)=-tau(MOD(axis+1,3)+1,na)
               zeu_u(p,i,MOD(axis+1,3)+1,na)=tau(MOD(axis,3)+1,na)
            enddo
            !
         enddo
      endif
      !
      if (n.eq.6) then
         do i=1,3
            do j=1,3
               ! These are the 3*3 vectors associated with the
               ! three rotational sum rules (0D system - typ. molecule)
               p=p+1
               do na=1,nat
                  zeu_u(p,i,MOD(j,3)+1,na)=-tau(MOD(j+1,3)+1,na)
                  zeu_u(p,i,MOD(j+1,3)+1,na)=tau(MOD(j,3)+1,na)
               enddo
               !
            enddo
         enddo
      endif
      !
      ! Gram-Schmidt orthonormalization of the set of vectors created.
      !
      nzeu_less=0
      do k=1,p
         zeu_w(:,:,:)=zeu_u(k,:,:,:)
         zeu_x(:,:,:)=zeu_u(k,:,:,:)
         do q=1,k-1
            r=1
            do izeu_less=1,nzeu_less
               if (zeu_less(izeu_less).eq.q) r=0
            enddo
            if (r.ne.0) then
               call sp_zeu(zeu_x,zeu_u(q,:,:,:),nat,scal)
               zeu_w(:,:,:) = zeu_w(:,:,:) - scal* zeu_u(q,:,:,:)
            endif
         enddo
         call sp_zeu(zeu_w,zeu_w,nat,norm2)
         if (norm2.gt.1.0d-16) then
            zeu_u(k,:,:,:) = zeu_w(:,:,:) / DSQRT(norm2)
         else
            nzeu_less=nzeu_less+1
            zeu_less(nzeu_less)=k
         endif
      enddo
      !
      ! Projection of the effective charge "vector" on the orthogonal of the
      ! subspace of the vectors verifying the sum rules
      !
      zeu_w(:,:,:)=0.0d0
      do k=1,p
         r=1
         do izeu_less=1,nzeu_less
            if (zeu_less(izeu_less).eq.k) r=0
         enddo
         if (r.ne.0) then
            zeu_x(:,:,:)=zeu_u(k,:,:,:)
            call sp_zeu(zeu_x,zeu_new,nat,scal)
            zeu_w(:,:,:) = zeu_w(:,:,:) + scal*zeu_u(k,:,:,:)
         endif
      enddo
      !
      ! Final substraction of the former projection to the initial zeu, to get
      ! the new "projected" zeu
      !
      zeu_new(:,:,:)=zeu_new(:,:,:) - zeu_w(:,:,:)
      call sp_zeu(zeu_w,zeu_w,nat,norm2)
      write(stdout,'("Norm of the difference between old and new effective ", &
           &  "charges: " , F25.20)') SQRT(norm2)
      !
      ! Check projection
      !
      !write(6,'("Check projection of zeu")')
      !do k=1,p
      !  zeu_x(:,:,:)=zeu_u(k,:,:,:)
      !  call sp_zeu(zeu_x,zeu_new,nat,scal)
      !  if (DABS(scal).gt.1d-10) write(6,'("k= ",I8," zeu_new|zeu_u(k)= ",F15.10)') k,scal
      !enddo
      !
      do i=1,3
         do j=1,3
            do na=1,nat
               zeu(i,j,na)=zeu_new(i,j,na)
            enddo
         enddo
      enddo
   endif
   !
   !
   return
 end subroutine set_zasr
!
!----------------------------------------------------------------------
SUBROUTINE set_asr (asr, nr1, nr2, nr3, frc, zeu, nat, ibrav, tau)
  !-----------------------------------------------------------------------
  !
  USE kinds,      ONLY : DP
  USE io_global,  ONLY : stdout
  !
  IMPLICIT NONE
  CHARACTER (LEN=10), intent(in) :: asr
  INTEGER, intent(in) :: nr1, nr2, nr3, nat, ibrav
  REAL(DP), intent(in) :: tau(3,nat)
  REAL(DP), intent(inout) :: frc(nr1,nr2,nr3,3,3,nat,nat), zeu(3,3,nat)
  !
  INTEGER :: axis, n, i, j, na, nb, n1,n2,n3, m,p,k,l,q,r, i1,j1,na1
  REAL(DP) :: zeu_new(3,3,nat)
  REAL(DP), ALLOCATABLE :: frc_new(:,:,:,:,:,:,:)
  type vector
     real(DP),pointer :: vec(:,:,:,:,:,:,:)
  end type vector
  !
  type (vector) u(6*3*nat)
  ! These are the "vectors" associated with the sum rules on force-constants
  !
  integer :: u_less(6*3*nat),n_less,i_less
  ! indices of the vectors u that are not independent to the preceding ones,
  ! n_less = number of such vectors, i_less = temporary parameter
  !
  integer, allocatable :: ind_v(:,:,:)
  real(DP), allocatable :: v(:,:)
  ! These are the "vectors" associated with symmetry conditions, coded by
  ! indicating the positions (i.e. the seven indices) of the non-zero elements (there
  ! should be only 2 of them) and the value of that element. We do so in order
  ! to limit the amount of memory used.
  !
  real(DP), allocatable :: w(:,:,:,:,:,:,:), x(:,:,:,:,:,:,:)
  ! temporary vectors and parameters
  real(DP) :: scal,norm2, sum
  !
  real(DP) :: zeu_u(6*3,3,3,nat)
  ! These are the "vectors" associated with the sum rules on effective charges
  !
  integer :: zeu_less(6*3),nzeu_less,izeu_less
  ! indices of the vectors zeu_u that are not independent to the preceding ones,
  ! nzeu_less = number of such vectors, izeu_less = temporary parameter
  !
  real(DP) :: zeu_w(3,3,nat), zeu_x(3,3,nat)
  ! temporary vectors

  ! Initialization. n is the number of sum rules to be considered (if asr.ne.'simple')
  ! and 'axis' is the rotation axis in the case of a 1D system
  ! (i.e. the rotation axis is (Ox) if axis='1', (Oy) if axis='2' and (Oz) if axis='3')
  !
  if((asr.ne.'simple').and.(asr.ne.'crystal').and.(asr.ne.'one-dim') &
                      .and.(asr.ne.'zero-dim')) then
     call errore('set_asr','invalid Acoustic Sum Rule:' // asr, 1)
  endif
  !
  if(asr.eq.'simple') then
     !
     ! Simple Acoustic Sum Rule on effective charges
     !
     do i=1,3
        do j=1,3
           sum=0.0d0
           do na=1,nat
              sum = sum + zeu(i,j,na)
           end do
           do na=1,nat
              zeu(i,j,na) = zeu(i,j,na) - sum/nat
           end do
        end do
     end do
     !
     ! Simple Acoustic Sum Rule on force constants in real space
     !
     do i=1,3
        do j=1,3
           do na=1,nat
              sum=0.0d0
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
      !
   end if

  if(asr.eq.'crystal') n=3
  if(asr.eq.'one-dim') then
     ! the direction of periodicity is the rotation axis
     ! It will work only if the crystal axis considered is one of
     ! the cartesian axis (typically, ibrav=1, 6 or 8, or 4 along the
     ! z-direction)
     if (nr1*nr2*nr3.eq.1) axis=3
     if ((nr1.ne.1).and.(nr2*nr3.eq.1)) axis=1
     if ((nr2.ne.1).and.(nr1*nr3.eq.1)) axis=2
     if ((nr3.ne.1).and.(nr1*nr2.eq.1)) axis=3
     if (((nr1.ne.1).and.(nr2.ne.1)).or.((nr2.ne.1).and. &
          (nr3.ne.1)).or.((nr1.ne.1).and.(nr3.ne.1))) then
        call errore('set_asr','too many directions of &
             & periodicity in 1D system',axis)
     endif
     if ((ibrav.ne.1).and.(ibrav.ne.6).and.(ibrav.ne.8).and. &
          ((ibrav.ne.4).or.(axis.ne.3)) ) then
        write(stdout,*) 'asr: rotational axis may be wrong'
     endif
     write(stdout,'("asr rotation axis in 1D system= ",I4)') axis
     n=4
  endif
  if(asr.eq.'zero-dim') n=6
  !
  ! Acoustic Sum Rule on effective charges
  !
  ! generating the vectors of the orthogonal of the subspace to project
  ! the effective charges matrix on
  !
  zeu_u(:,:,:,:)=0.0d0
  do i=1,3
     do j=1,3
        do na=1,nat
           zeu_new(i,j,na)=zeu(i,j,na)
        enddo
     enddo
  enddo
  !
  p=0
  do i=1,3
     do j=1,3
        ! These are the 3*3 vectors associated with the
        ! translational acoustic sum rules
        p=p+1
        zeu_u(p,i,j,:)=1.0d0
        !
     enddo
  enddo
  !
  if (n.eq.4) then
     do i=1,3
        ! These are the 3 vectors associated with the
        ! single rotational sum rule (1D system)
        p=p+1
        do na=1,nat
           zeu_u(p,i,MOD(axis,3)+1,na)=-tau(MOD(axis+1,3)+1,na)
           zeu_u(p,i,MOD(axis+1,3)+1,na)=tau(MOD(axis,3)+1,na)
        enddo
        !
     enddo
  endif
  !
  if (n.eq.6) then
     do i=1,3
        do j=1,3
           ! These are the 3*3 vectors associated with the
           ! three rotational sum rules (0D system - typ. molecule)
           p=p+1
           do na=1,nat
              zeu_u(p,i,MOD(j,3)+1,na)=-tau(MOD(j+1,3)+1,na)
              zeu_u(p,i,MOD(j+1,3)+1,na)=tau(MOD(j,3)+1,na)
           enddo
           !
        enddo
     enddo
  endif
  !
  ! Gram-Schmidt orthonormalization of the set of vectors created.
  !
  nzeu_less=0
  do k=1,p
     zeu_w(:,:,:)=zeu_u(k,:,:,:)
     zeu_x(:,:,:)=zeu_u(k,:,:,:)
     do q=1,k-1
        r=1
        do izeu_less=1,nzeu_less
           if (zeu_less(izeu_less).eq.q) r=0
        enddo
        if (r.ne.0) then
           call sp_zeu(zeu_x,zeu_u(q,:,:,:),nat,scal)
           zeu_w(:,:,:) = zeu_w(:,:,:) - scal* zeu_u(q,:,:,:)
        endif
     enddo
     call sp_zeu(zeu_w,zeu_w,nat,norm2)
     if (norm2.gt.1.0d-16) then
        zeu_u(k,:,:,:) = zeu_w(:,:,:) / DSQRT(norm2)
     else
        nzeu_less=nzeu_less+1
        zeu_less(nzeu_less)=k
     endif
  enddo
  !
  ! Projection of the effective charge "vector" on the orthogonal of the
  ! subspace of the vectors verifying the sum rules
  !
  zeu_w(:,:,:)=0.0d0
  do k=1,p
     r=1
     do izeu_less=1,nzeu_less
        if (zeu_less(izeu_less).eq.k) r=0
     enddo
     if (r.ne.0) then
        zeu_x(:,:,:)=zeu_u(k,:,:,:)
        call sp_zeu(zeu_x,zeu_new,nat,scal)
        zeu_w(:,:,:) = zeu_w(:,:,:) + scal*zeu_u(k,:,:,:)
     endif
  enddo
  !
  ! Final substraction of the former projection to the initial zeu, to get
  ! the new "projected" zeu
  !
  zeu_new(:,:,:)=zeu_new(:,:,:) - zeu_w(:,:,:)
  call sp_zeu(zeu_w,zeu_w,nat,norm2)
  write(stdout,'("Norm of the difference between old and new effective ", &
       & "charges: ",F25.20)') SQRT(norm2)
  !
  ! Check projection
  !
  !write(6,'("Check projection of zeu")')
  !do k=1,p
  !  zeu_x(:,:,:)=zeu_u(k,:,:,:)
  !  call sp_zeu(zeu_x,zeu_new,nat,scal)
  !  if (DABS(scal).gt.1d-10) write(6,'("k= ",I8," zeu_new|zeu_u(k)= ",F15.10)') k,scal
  !enddo
  !
  do i=1,3
     do j=1,3
        do na=1,nat
           zeu(i,j,na)=zeu_new(i,j,na)
        enddo
     enddo
  enddo
  !
  ! Acoustic Sum Rule on force constants
  !
  !
  ! generating the vectors of the orthogonal of the subspace to project
  ! the force-constants matrix on
  !
  do k=1,18*nat
     allocate(u(k) % vec(nr1,nr2,nr3,3,3,nat,nat))
     u(k) % vec (:,:,:,:,:,:,:)=0.0d0
  enddo
  ALLOCATE (frc_new(nr1,nr2,nr3,3,3,nat,nat))
  do i=1,3
     do j=1,3
        do na=1,nat
           do nb=1,nat
              do n1=1,nr1
                 do n2=1,nr2
                    do n3=1,nr3
                       frc_new(n1,n2,n3,i,j,na,nb)=frc(n1,n2,n3,i,j,na,nb)
                    enddo
                 enddo
              enddo
           enddo
        enddo
     enddo
  enddo
  !
  p=0
  do i=1,3
     do j=1,3
        do na=1,nat
           ! These are the 3*3*nat vectors associated with the
           ! translational acoustic sum rules
           p=p+1
           u(p) % vec (:,:,:,i,j,na,:)=1.0d0
           !
        enddo
     enddo
  enddo
  !
  if (n.eq.4) then
     do i=1,3
        do na=1,nat
           ! These are the 3*nat vectors associated with the
           ! single rotational sum rule (1D system)
           p=p+1
           do nb=1,nat
              u(p) % vec (:,:,:,i,MOD(axis,3)+1,na,nb)=-tau(MOD(axis+1,3)+1,nb)
              u(p) % vec (:,:,:,i,MOD(axis+1,3)+1,na,nb)=tau(MOD(axis,3)+1,nb)
           enddo
           !
        enddo
     enddo
  endif
  !
  if (n.eq.6) then
     do i=1,3
        do j=1,3
           do na=1,nat
              ! These are the 3*3*nat vectors associated with the
              ! three rotational sum rules (0D system - typ. molecule)
              p=p+1
              do nb=1,nat
                 u(p) % vec (:,:,:,i,MOD(j,3)+1,na,nb)=-tau(MOD(j+1,3)+1,nb)
                 u(p) % vec (:,:,:,i,MOD(j+1,3)+1,na,nb)=tau(MOD(j,3)+1,nb)
              enddo
              !
           enddo
        enddo
     enddo
  endif
  !
  allocate (ind_v(9*nat*nat*nr1*nr2*nr3,2,7), v(9*nat*nat*nr1*nr2*nr3,2) )
  m=0
  do i=1,3
     do j=1,3
        do na=1,nat
           do nb=1,nat
              do n1=1,nr1
                 do n2=1,nr2
                    do n3=1,nr3
                       ! These are the vectors associated with the symmetry constraints
                       q=1
                       l=1
                       do while((l.le.m).and.(q.ne.0))
                          if ((ind_v(l,1,1).eq.n1).and.(ind_v(l,1,2).eq.n2).and. &
                               (ind_v(l,1,3).eq.n3).and.(ind_v(l,1,4).eq.i).and. &
                               (ind_v(l,1,5).eq.j).and.(ind_v(l,1,6).eq.na).and. &
                               (ind_v(l,1,7).eq.nb)) q=0
                          if ((ind_v(l,2,1).eq.n1).and.(ind_v(l,2,2).eq.n2).and. &
                               (ind_v(l,2,3).eq.n3).and.(ind_v(l,2,4).eq.i).and. &
                               (ind_v(l,2,5).eq.j).and.(ind_v(l,2,6).eq.na).and. &
                               (ind_v(l,2,7).eq.nb)) q=0
                          l=l+1
                       enddo
                       if ((n1.eq.MOD(nr1+1-n1,nr1)+1).and.(n2.eq.MOD(nr2+1-n2,nr2)+1) &
                            .and.(n3.eq.MOD(nr3+1-n3,nr3)+1).and.(i.eq.j).and.(na.eq.nb)) q=0
                       if (q.ne.0) then
                          m=m+1
                          ind_v(m,1,1)=n1
                          ind_v(m,1,2)=n2
                          ind_v(m,1,3)=n3
                          ind_v(m,1,4)=i
                          ind_v(m,1,5)=j
                          ind_v(m,1,6)=na
                          ind_v(m,1,7)=nb
                          v(m,1)=1.0d0/DSQRT(2.0d0)
                          ind_v(m,2,1)=MOD(nr1+1-n1,nr1)+1
                          ind_v(m,2,2)=MOD(nr2+1-n2,nr2)+1
                          ind_v(m,2,3)=MOD(nr3+1-n3,nr3)+1
                          ind_v(m,2,4)=j
                          ind_v(m,2,5)=i
                          ind_v(m,2,6)=nb
                          ind_v(m,2,7)=na
                          v(m,2)=-1.0d0/DSQRT(2.0d0)
                       endif
                    enddo
                 enddo
              enddo
           enddo
        enddo
     enddo
  enddo
  !
  ! Gram-Schmidt orthonormalization of the set of vectors created.
  ! Note that the vectors corresponding to symmetry constraints are already
  ! orthonormalized by construction.
  !
  n_less=0
  allocate (w(nr1,nr2,nr3,3,3,nat,nat), x(nr1,nr2,nr3,3,3,nat,nat))
  do k=1,p
     w(:,:,:,:,:,:,:)=u(k) % vec (:,:,:,:,:,:,:)
     x(:,:,:,:,:,:,:)=u(k) % vec (:,:,:,:,:,:,:)
     do l=1,m
        !
        call sp2(x,v(l,:),ind_v(l,:,:),nr1,nr2,nr3,nat,scal)
        do r=1,2
           n1=ind_v(l,r,1)
           n2=ind_v(l,r,2)
           n3=ind_v(l,r,3)
           i=ind_v(l,r,4)
           j=ind_v(l,r,5)
           na=ind_v(l,r,6)
           nb=ind_v(l,r,7)
           w(n1,n2,n3,i,j,na,nb)=w(n1,n2,n3,i,j,na,nb)-scal*v(l,r)
        enddo
     enddo
     if (k.le.(9*nat)) then
        na1=MOD(k,nat)
        if (na1.eq.0) na1=nat
        j1=MOD((k-na1)/nat,3)+1
        i1=MOD((((k-na1)/nat)-j1+1)/3,3)+1
     else
        q=k-9*nat
        if (n.eq.4) then
           na1=MOD(q,nat)
           if (na1.eq.0) na1=nat
           i1=MOD((q-na1)/nat,3)+1
        else
           na1=MOD(q,nat)
           if (na1.eq.0) na1=nat
           j1=MOD((q-na1)/nat,3)+1
           i1=MOD((((q-na1)/nat)-j1+1)/3,3)+1
        endif
     endif
     do q=1,k-1
        r=1
        do i_less=1,n_less
           if (u_less(i_less).eq.q) r=0
        enddo
        if (r.ne.0) then
           call sp3(x,u(q) % vec (:,:,:,:,:,:,:), i1,na1,nr1,nr2,nr3,nat,scal)
           w(:,:,:,:,:,:,:) = w(:,:,:,:,:,:,:) - scal* u(q) % vec (:,:,:,:,:,:,:)
        endif
     enddo
     call sp1(w,w,nr1,nr2,nr3,nat,norm2)
     if (norm2.gt.1.0d-16) then
        u(k) % vec (:,:,:,:,:,:,:) = w(:,:,:,:,:,:,:) / DSQRT(norm2)
     else
        n_less=n_less+1
        u_less(n_less)=k
     endif
  enddo
  !
  ! Projection of the force-constants "vector" on the orthogonal of the
  ! subspace of the vectors verifying the sum rules and symmetry contraints
  !
  w(:,:,:,:,:,:,:)=0.0d0
  do l=1,m
     call sp2(frc_new,v(l,:),ind_v(l,:,:),nr1,nr2,nr3,nat,scal)
     do r=1,2
        n1=ind_v(l,r,1)
        n2=ind_v(l,r,2)
        n3=ind_v(l,r,3)
        i=ind_v(l,r,4)
        j=ind_v(l,r,5)
        na=ind_v(l,r,6)
        nb=ind_v(l,r,7)
        w(n1,n2,n3,i,j,na,nb)=w(n1,n2,n3,i,j,na,nb)+scal*v(l,r)
     enddo
  enddo
  do k=1,p
     r=1
     do i_less=1,n_less
        if (u_less(i_less).eq.k) r=0
     enddo
     if (r.ne.0) then
        x(:,:,:,:,:,:,:)=u(k) % vec (:,:,:,:,:,:,:)
        call sp1(x,frc_new,nr1,nr2,nr3,nat,scal)
        w(:,:,:,:,:,:,:) = w(:,:,:,:,:,:,:) + scal*u(k)%vec(:,:,:,:,:,:,:)
     endif
     deallocate(u(k) % vec)
  enddo
  !
  ! Final substraction of the former projection to the initial frc, to get
  ! the new "projected" frc
  !
  frc_new(:,:,:,:,:,:,:)=frc_new(:,:,:,:,:,:,:) - w(:,:,:,:,:,:,:)
  call sp1(w,w,nr1,nr2,nr3,nat,norm2)
  write(stdout,'("Norm of the difference between old and new force-constants:",&
       &     F25.20)') SQRT(norm2)
  !
  ! Check projection
  !
  !write(6,'("Check projection IFC")')
  !do l=1,m
  !  call sp2(frc_new,v(l,:),ind_v(l,:,:),nr1,nr2,nr3,nat,scal)
  !  if (DABS(scal).gt.1d-10) write(6,'("l= ",I8," frc_new|v(l)= ",F15.10)') l,scal
  !enddo
  !do k=1,p
  !  x(:,:,:,:,:,:,:)=u(k) % vec (:,:,:,:,:,:,:)
  !  call sp1(x,frc_new,nr1,nr2,nr3,nat,scal)
  !  if (DABS(scal).gt.1d-10) write(6,'("k= ",I8," frc_new|u(k)= ",F15.10)') k,scal
  !  deallocate(u(k) % vec)
  !enddo
  !
  do i=1,3
     do j=1,3
        do na=1,nat
           do nb=1,nat
              do n1=1,nr1
                 do n2=1,nr2
                    do n3=1,nr3
                       frc(n1,n2,n3,i,j,na,nb)=frc_new(n1,n2,n3,i,j,na,nb)
                    enddo
                 enddo
              enddo
           enddo
        enddo
     enddo
  enddo
  deallocate (x, w)
  deallocate (v, ind_v)
  deallocate (frc_new)
  !
  return
end subroutine set_asr
!

!----------------------------------------------------------------------
subroutine sp_zeu(zeu_u,zeu_v,nat,scal)
  !-----------------------------------------------------------------------
  !
  ! does the scalar product of two effective charges matrices zeu_u and zeu_v
  ! (considered as vectors in the R^(3*3*nat) space, and coded in the usual way)
  !
  USE kinds, ONLY : DP
  implicit none
  integer i,j,na,nat
  real(DP) zeu_u(3,3,nat)
  real(DP) zeu_v(3,3,nat)
  real(DP) scal
  !
  !
  scal=0.0d0
  do i=1,3
    do j=1,3
      do na=1,nat
        scal=scal+zeu_u(i,j,na)*zeu_v(i,j,na)
      enddo
    enddo
  enddo
  !
  return
  !
end subroutine sp_zeu
!----------------------------------------------------------------------
subroutine sp1(u,v,nr1,nr2,nr3,nat,scal)
  !-----------------------------------------------------------------------
  !
  ! does the scalar product of two force-constants matrices u and v (considered as
  ! vectors in the R^(3*3*nat*nat*nr1*nr2*nr3) space, and coded in the usual way)
  !
  USE kinds, ONLY: DP
  implicit none
  integer nr1,nr2,nr3,i,j,na,nb,n1,n2,n3,nat
  real(DP) u(nr1,nr2,nr3,3,3,nat,nat)
  real(DP) v(nr1,nr2,nr3,3,3,nat,nat)
  real(DP) scal
  !
  !
  scal=0.0d0
  do i=1,3
    do j=1,3
      do na=1,nat
        do nb=1,nat
          do n1=1,nr1
            do n2=1,nr2
              do n3=1,nr3
                scal=scal+u(n1,n2,n3,i,j,na,nb)*v(n1,n2,n3,i,j,na,nb)
              enddo
            enddo
          enddo
        enddo
      enddo
    enddo
  enddo
  !
  return
  !
end subroutine sp1
!
!----------------------------------------------------------------------
subroutine sp2(u,v,ind_v,nr1,nr2,nr3,nat,scal)
  !-----------------------------------------------------------------------
  !
  ! does the scalar product of two force-constants matrices u and v (considered as
  ! vectors in the R^(3*3*nat*nat*nr1*nr2*nr3) space). u is coded in the usual way
  ! but v is coded as explained when defining the vectors corresponding to the
  ! symmetry constraints
  !
  USE kinds, ONLY: DP
  implicit none
  integer nr1,nr2,nr3,i,nat
  real(DP) u(nr1,nr2,nr3,3,3,nat,nat)
  integer ind_v(2,7)
  real(DP) v(2)
  real(DP) scal
  !
  !
  scal=0.0d0
  do i=1,2
    scal=scal+u(ind_v(i,1),ind_v(i,2),ind_v(i,3),ind_v(i,4),ind_v(i,5),ind_v(i,6), &
         ind_v(i,7))*v(i)
  enddo
  !
  return
  !
end subroutine sp2
!
!----------------------------------------------------------------------
subroutine sp3(u,v,i,na,nr1,nr2,nr3,nat,scal)
  !-----------------------------------------------------------------------
  !
  ! like sp1, but in the particular case when u is one of the u(k)%vec
  ! defined in set_asr (before orthonormalization). In this case most of the
  ! terms are zero (the ones that are not are characterized by i and na), so
  ! that a lot of computer time can be saved (during Gram-Schmidt).
  !
  USE kinds, ONLY: DP
  implicit none
  integer nr1,nr2,nr3,i,j,na,nb,n1,n2,n3,nat
  real(DP) u(nr1,nr2,nr3,3,3,nat,nat)
  real(DP) v(nr1,nr2,nr3,3,3,nat,nat)
  real(DP) scal
  !
  !
  scal=0.0d0
  do j=1,3
    do nb=1,nat
      do n1=1,nr1
        do n2=1,nr2
          do n3=1,nr3
            scal=scal+u(n1,n2,n3,i,j,na,nb)*v(n1,n2,n3,i,j,na,nb)
          enddo
        enddo
      enddo
    enddo
  enddo
  !
  return
  !
end subroutine sp3
!

