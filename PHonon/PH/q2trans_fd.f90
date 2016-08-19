!
! Copyright (C) 2001-2008 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
PROGRAM q2trans_fd
  !----------------------------------------------------------------------------
  !
  !     reads force constant matrices C(q) produced by the cwFD phonon code
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
  USE io_global, ONLY : ionode_id, ionode, stdout
  USE io_dyn_mat, ONLY : read_dyn_mat_param, read_dyn_mat_header, &
                         read_dyn_mat, read_dyn_mat_tail, &
                         write_dyn_mat_header, write_ifc
  USE environment, ONLY : environment_start, environment_end
USE constants, ONLY: pi, fpi, e2
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
  CHARACTER(len=6), EXTERNAL :: int_to_char
  !
  LOGICAL :: lq, lrigid, lrigid1, lnogridinfo, xmldyn, ltrans
  CHARACTER (len=10) :: zasr, iasr
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
  LOGICAL           :: la2F, onedim
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

     ALLOCATE (ityp(nat))

  ! Read the IFCs from file

filin=trim(file_ifc)//'.fc'
OPEN(3,FILE=trim(filin))

READ(3,*) ntyp, nat
DO na=1,ntyp
   READ(3,*)
ENDDO
DO na=1,nat
   READ (3,*) ndummy, ityp(na)
ENDDO

READ(3,*) dielec
IF(dielec) THEN
DO i=1,3
   READ(3,*)
ENDDO
DO na=1,nat
   READ(3,*)
   DO i=1,3
      READ(3,*)
   ENDDO
ENDDO
ENDIF

READ (3,'(3i4)') ndummy,ndummy,ndummy

DO i=1,3
   DO j=1,3
      DO na=1,nat
         DO nb=1,nat
            READ (3,*)
            DO m1=1,nr1
              DO m2=1,nr2
                 DO m3=1,nr3
                    READ (3,'(3i4,2x,1pe18.11)')ndummy,ndummy,ndummy,ifc3(i,j,na,nb,m1,m2,m3)
                 ENDDO
               ENDDO
             ENDDO
         ENDDO
      ENDDO
   ENDDO
ENDDO

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
                          IF(m2==1.and.m3==1) ifc(j1,j2,na1,na2,m1)=ifc3(j1,j2,na1,na2,m1,m2,m3)
                       ENDIF
                       ! for transport in 3-dim systems: sum on the plane perpendicular to the transport direction
                        ifc(j1,j2,na1,na2,m1)=ifc(j1,j2,na1,na2,m1)+ifc3(j1,j2,na1,na2,m1,m2,m3)
                    ENDDO
                 ENDDO
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ENDDO
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
                          IF(m1==1.and.m3==1) ifc(j1,j2,na1,na2,m2)=ifc3(j1,j2,na1,na2,m1,m2,m3)
                       ENDIF
                       ! for transport in 3-dim systems: sum on the plane perpendicular to the transport direction
                        ifc(j1,j2,na1,na2,m2)=ifc(j1,j2,na1,na2,m2)+ifc3(j1,j2,na1,na2,m1,m2,m3)
                    ENDDO
                 ENDDO
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ENDDO
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
                          IF(m1==1.and.m2==1) ifc(j1,j2,na1,na2,m3)=ifc3(j1,j2,na1,na2,m1,m2,m3)
                       ENDIF
                       ! for transport in 3-dim systems: sum on the plane perpendicular to the transport direction
                        ifc(j1,j2,na1,na2,m3)=ifc(j1,j2,na1,na2,m3)+ifc3(j1,j2,na1,na2,m1,m2,m3)
                    ENDDO
                 ENDDO
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ENDDO
  nrx=nr3

  END SELECT direction

  ! Correction for finite IFC in the center of the real space mesh
  IF(nrx > 1) THEN
  DO j1=1,3
     DO j2=1,3
        DO na1=1,nat
           DO na2=1,nat
              ifc0(j1,j2,na1,na2)= ifc(j1,j2,na1,na2,nrx/2+1)
           ENDDO
        ENDDO
     ENDDO
  ENDDO

  DO j1=1,3
     DO j2=1,3
        DO na1=1,nat
           DO na2=1,nat
              DO m1=1,nrx
                 ifc(j1,j2,na1,na2,m1)= ifc(j1,j2,na1,na2,m1)-ifc0(j1,j2,na1,na2)
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ENDDO
  ENDIF

  ! Impose the acoustic sum rule for the shifted IFC: the interatomic force of the atom on itself should be
  ! equal to minus the sum of all interatomic forces generated by all others atoms (action-reaction law!)
  ! eq. (82) in Gonze and Lee, PRB 55, 10355 (1997)

  DO j1=1,3
     DO j2=1,3
        DO na1=1,nat
           sum1=0.0
           DO na2=1,nat
              IF(na1/=na2) sum1=sum1+ifc(j1,j2,na1,na2,1)
           ENDDO
           sum2=0.0
           DO na2=1,nat
              DO m1=2,nrx
                 sum2=sum2+ifc(j1,j2,na1,na2,m1)
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ENDDO

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
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ENDDO

  ! define k00

  DO i=1,3*nat*(nrx/2+1)
     DO j=1,3*nat*(nrx/2+1)
        k00(i,j)=0.0d0
     ENDDO
  ENDDO

  DO m1=1,nrx/2+1
     DO m2=m1,nrx/2+1
        DO j1=1,3*nat
           DO j2=1,3*nat
              k00(3*nat*(m1-1)+j1,3*nat*(m2-1)+j2 ) = &
                   amconv/sqrt(amass(ityp((j1-1)/3+1))*amass(ityp((j2-1)/3+1)))*kfc(j1,j2,m2-m1+1)
           ENDDO
        ENDDO
     ENDDO
  ENDDO

  DO i=1,3*nat*(nrx/2+1)
     DO j=1,i-1
        k00(i,j)=k00(j,i)
     ENDDO
  ENDDO

  ! define k01

  DO i=1,3*nat*(nrx/2+1)
     DO j=1,3*nat*(nrx/2+1)
        k01(i,j)=0.0d0
     ENDDO
  ENDDO

  DO m1=1,nrx/2+1
     DO m2=1,m1-1
        DO j1=1,3*nat
           DO j2=1,3*nat
              k01(3*nat*(m1-1)+j1,3*nat*(m2-1)+j2 ) = &
                   amconv/sqrt(amass(ityp((j1-1)/3+1))*amass(ityp((j2-1)/3+1)))*kfc(j1,j2,m2+(nrx/2+1)-m1+1)
           ENDDO
        ENDDO
     ENDDO
  ENDDO

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
DO i=1,36
   WRITE(77,'(36(1x,2f12.8))')(rham(i,j,1),j=1,36)
ENDDO

! 01 - 13
DO i=13,24
   DO j=1,12
      rham(i,j,2)=rham(i-12,j+24,1)
   ENDDO
ENDDO

! 01 - 13
DO i=25,36
   DO j=13,24
      rham(i,j,2)=rham(i-24,j+12,1)
   ENDDO
ENDDO

! 01 - 12
DO i=25,36
   DO j=1,12
      rham(i,j,2)=rham(i-24,j+12,1)
   ENDDO
ENDDO

DO i=1,36
   WRITE(78,'(24(1x,2f12.8))')(rham(i,j,2),j=1,36)
ENDDO

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

CALL iotk_open_write( stdout, FILE=trim(fileout))

CALL iotk_write_begin(stdout,"HAMILTONIAN")

CALL iotk_write_attr( attr, "dimwann", dimwan, FIRST=.true. )
CALL iotk_write_attr( attr, "nkpts", nkpts )
CALL iotk_write_attr( attr, "nk", nk )
CALL iotk_write_attr( attr, "nrtot", nrtot )
CALL iotk_write_attr( attr, "nr", nr )
CALL iotk_write_attr( attr, "have_overlap", have_overlap )
CALL iotk_write_attr( attr, "fermi_energy", fermi_energy )

CALL iotk_write_empty( stdout, "DATA", ATTR=attr)

nsize=3*2
CALL iotk_write_attr( attr, "type", "integer", FIRST=.true. )
CALL iotk_write_attr( attr, "size", nsize )
CALL iotk_write_attr( attr, "columns", 3 )
CALL iotk_write_attr( attr, "units", "crystal" )
CALL iotk_write_dat( stdout, "IVR", ivr, COLUMNS=3, ATTR=attr )

CALL iotk_write_attr( attr, "type", "real", FIRST=.true. )
CALL iotk_write_attr( attr, "size", nkpts )
CALL iotk_write_dat( stdout, "WR", wr, ATTR=attr )

CALL iotk_write_begin(stdout,"RHAM")
DO ir = 1, nrtot
    CALL iotk_write_dat(stdout,"VR"//trim(iotk_index(ir)), rham(:,:,ir))
ENDDO
CALL iotk_write_end(stdout,"RHAM")
CALL iotk_write_end(stdout,"HAMILTONIAN")

CALL iotk_close_write( stdout )


  !
  CALL environment_end('Q2R')

  CALL mp_global_end()
  !
END PROGRAM q2trans_fd
