!
! Copyright (C) 2001-2004 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "f_defs.h"
!
!----------------------------------------------------------------------------
PROGRAM q2r
  !----------------------------------------------------------------------------
  !
  USE kinds,      ONLY : DP
  USE mp,         ONLY : mp_start, mp_env, mp_end, mp_barrier
  USE mp_global,  ONLY : nproc, mpime, mp_global_start
  !
  IMPLICIT NONE
  !
  INTEGER,       PARAMETER :: nax = 16, nrx1 = 8, nrx2 = 8, nrx3 = 8
  REAL(kind=DP), PARAMETER :: eps=1.D-5
  INTEGER                  :: nr1, nr2, nr3, nr(3)
  !
  CHARACTER(len=20) :: crystal
  CHARACTER(len=80) :: title
  CHARACTER(len=80) :: filin,filj,filf,fild
  CHARACTER(len=3)  :: atm(nax)
  !
  LOGICAL :: lq,lrigid,zasr, lrigid_save 
  INTEGER :: m1, m2, m3, l1, l2, l3, i, j, j1, j2, na1, na2, ipol
  INTEGER :: nat, nq, ntyp, iq, icar, nfile, nqtot, ifile, nqs
  INTEGER :: na, nt
  !
  INTEGER :: gid
  !
  INTEGER :: m(3), nc(nrx1,nrx2,nrx3),ibrav,ityp(nax)
  !
  REAL(KIND=DP) :: celldm(6), at(3,3), bg(3,3), tau(3,nax)
  REAL(KIND=DP) :: q(3,48),omega, xq, amass(nax), resi,sum
  REAL(KIND=DP) :: epsil(3,3),zeu(3,3,nax)
  !
  COMPLEX(KIND=DP) :: phiq(3,3,nax,nax,48)
  COMPLEX(KIND=DP) :: phid(nrx1,nrx2,nrx3,3,3,nax,nax)
  !
  NAMELIST / input / nr1, nr2, nr3, fild, zasr
  !
  CHARACTER (LEN=80)  :: input_file
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
     nr1 = 0
     nr2 = 0
     nr3 = 0
     !
     ! ... Input from file ?
     !
     nargs = iargc()
     !
     DO iiarg = 1, ( nargs - 1 )
        !
        CALL getarg( iiarg, input_file )
        !
        IF ( TRIM( input_file ) == '-input' .OR. &
             TRIM( input_file ) == '-inp'   .OR. &
             TRIM( input_file ) == '-in' ) THEN
           !
           CALL getarg( ( iiarg + 1 ) , input_file )
           !
           OPEN ( UNIT = 5, FILE = input_file, FORM = 'FORMATTED', &
                STATUS = 'OLD', IOSTAT = ierr )
           !
           CALL errore( 'iosys', 'input file ' // TRIM( input_file ) // &
                & ' not found' , ierr )
           !
        END IF
        !
     END DO
     !
     READ ( 5, input )
     !
     ! check input
     !
     IF (nr1 > nrx1) CALL errore ('q2r',' nr1 too big, increase nrx1',nrx1)
     IF (nr2 > nrx2) CALL errore ('q2r',' nr2 too big, increase nrx2',nrx2)
     IF (nr3 > nrx3) CALL errore ('q2r',' nr3 too big, increase nrx3',nrx3)
     IF (nr1 < 1) CALL errore ('q2r',' nr1 wrong or missing',1)
     IF (nr2 < 1) CALL errore ('q2r',' nr2 wrong or missing',1)
     IF (nr3 < 1) CALL errore ('q2r',' nr3 wrong or missing',1)
     !
     ! copy nrX -> nr(X)
     !
     nr(1) = nr1
     nr(2) = nr2
     nr(3) = nr3
     !
     ! D matrix (analytical part)
     !
     !
     nqtot = 0
     !
     DO l1=1,nr1
        DO l2=1,nr2
           DO l3=1,nr3
              nc(l1,l2,l3)=0
           END DO
        END DO
     END DO
     !
     ! Reciprocal space dyn.mat. read from file
     !
     READ (5,*) nfile
     DO ifile=1,nfile
        READ(5,'(a)') filin
        WRITE (6,*) ' reading dyn.mat. from file ',filin
        OPEN(unit=1,file=filin,status='old',form='formatted')
        CALL read_file(nqs,q,phiq,nax,epsil,zeu,lrigid,  &
             ntyp,nat,ibrav,celldm,atm,amass,ityp,tau)
        IF (ifile.EQ.1) THEN
           lrigid_save=lrigid
           CALL latgen(ibrav,celldm,at(1,1),at(1,2),at(1,3),omega)
           at = at / celldm(1)  !  bring at in units of alat 
           CALL volume(celldm(1),at(1,1),at(1,2),at(1,3),omega)
           CALL recips(at(1,1),at(1,2),at(1,3),bg(1,1),bg(1,2),bg(1,3))
           IF (lrigid.AND.zasr) THEN
              DO i=1,3
                 DO j=1,3
                    sum=0.d0
                    DO na=1,nat
                       sum=sum+zeu(i,j,na)
                    END DO
                    DO na=1,nat
                       zeu(i,j,na)=zeu(i,j,na)-sum/nat
                    END DO
                 END DO
              END DO
           END IF
        END IF
        IF (lrigid.AND..NOT.lrigid_save) CALL errore('main',            &
             &          'in this case Gamma must be the first file ',1)
        !
        WRITE (6,*) ' nqs= ',nqs
        CLOSE(unit=1)
        DO nq = 1,nqs
           WRITE(6,'(a,3f12.8)') ' q= ',(q(i,nq),i=1,3)
           lq = .TRUE.
           DO ipol=1,3
              xq = 0.0
              DO icar=1,3
                 xq = xq + at(icar,ipol) * q(icar,nq) * nr(ipol)
              END DO
              lq = lq .AND. (ABS(NINT(xq) - xq) .LT. eps)
              iq = NINT(xq)
              !
              m(ipol)= MOD(iq,nr(ipol)) + 1
              IF (m(ipol) .LT. 1) m(ipol) = m(ipol) + nr(ipol) 
           END DO
           IF (.NOT.lq) CALL errore('init','q not allowed',1)
           IF(nc(m(1),m(2),m(3)).EQ.0) THEN
              nc(m(1),m(2),m(3))=1
              IF (lrigid) CALL rgd_blk (nax,nat,phiq(1,1,1,1,nq),q(1,nq), &
                   tau,epsil,zeu,bg,omega,-1.d0)
              CALL trasl(phid,phiq,nq,nrx1,nrx2,nrx3,nat,m(1),m(2),m(3),nax)
              nqtot=nqtot+1
           ELSE
              WRITE (*,'(3i4)') (m(i),i=1,3)
              CALL errore('init',' nc already filled: wrong q grid or wrong nr',1)
           END IF
        END DO
     END DO
     !
     ! Check grid dimension
     !
     IF (nqtot .EQ. nr1*nr2*nr3) THEN
        WRITE (6,'(/5x,a,i4)') ' q-space grid ok, #points = ',nqtot
     ELSE
        CALL errore('init',' missing q-point(s)!',1)
     END IF
     !
     ! dyn.mat. FFT
     !
     DO j1=1,3
        DO j2=1,3
           DO na1=1,nat
              DO na2=1,nat
                 CALL tolerant_cft3(phid(1,1,1,j1,j2,na1,na2), &
                      nr1,nr2,nr3,nrx1,nrx2,nrx3,1)
                 CALL DSCAL(2*nrx1*nrx2*nrx3,1.d0/(nr1*nr2*nr3),       &
                      phid(1,1,1,j1,j2,na1,na2),1)
              END DO
           END DO
        END DO
     END DO
     !
     ! Real space force constants written to file (analytical part)
     !
     resi = 0
     OPEN(unit=2,file=fild,status='unknown',form='formatted')
     WRITE(2,'(i3,i5,i3,6f11.7)') ntyp,nat,ibrav,celldm
     DO nt = 1,ntyp
        WRITE(2,*) nt," '",atm(nt),"' ",amass(nt)
     END DO
     DO na=1,nat
        WRITE(2,'(2i5,3f15.7)') na,ityp(na),(tau(j,na),j=1,3)
     END DO
     WRITE (2,*) lrigid
     IF (lrigid) THEN
        WRITE(2,'(3f15.7)') ((epsil(i,j),j=1,3),i=1,3)
        DO na=1,nat
           WRITE(2,'(i5)') na
           WRITE(2,'(3f15.7)') ((zeu(i,j,na),j=1,3),i=1,3)
        END DO
     END IF
     WRITE (2,'(4i4)') nr1, nr2, nr3 
     DO j1=1,3
        DO j2=1,3
           DO na1=1,nat
              DO na2=1,nat
                 DO m1=1,nr1
                    DO m2=1,nr2
                       DO m3=1,nr3
                          resi = resi + &
                                 dabs(dimag(phid(m1,m2,m3,j1,j2,na1,na2)))
                       END DO
                    END DO
                 END DO
                 WRITE (2,'(4i4)') j1,j2,na1,na2
                 WRITE (2,'(3i4,2x,1pe18.11)')   &
                      (((m1,m2,m3,REAL(phid(m1,m2,m3,j1,j2,na1,na2)), &
                      m1=1,nr1),m2=1,nr2),m3=1,nr3)
              END DO
           END DO
        END DO
     END DO
     CLOSE(2)
     WRITE (6,"(/5x,' fft-check: imaginary sum = ',e12.7)") resi
     !
  END IF
  ! 
  CALL mp_barrier()
  !
  CALL mp_end()
  !
END PROGRAM q2r
!
!----------------------------------------------------------------------------
SUBROUTINE read_file( nqs, xq, phi, nax, epsil, zeu, lrigid, &
                      ntyp, nat, ibrav, celldm, atm, amass, ityp, tau )
  !----------------------------------------------------------------------------
  !
  USE kinds, ONLY : DP
  !
  IMPLICIT NONE
  !
  ! I/O variables
  LOGICAL :: lrigid
  INTEGER :: nqs, nax, ntyp, nat, ibrav, ityp(nax)
  REAL(KIND=DP) :: epsil(3,3),zeu(3,3,nax)
  REAL(KIND=DP) :: xq(3,48), celldm(6), amass(nax), tau(3,nax)
  COMPLEX(KIND=DP) :: phi(3,3,nax,nax,48)
  CHARACTER(LEN=3) atm(nax)
  ! local variables
  INTEGER :: ntyp1,nat1,ibrav1,ityp1
  INTEGER :: i, j, na, nb, nt
  REAL(KIND=DP) :: tau1(3), amass1, celldm1(6),q2
  REAL(KIND=DP) :: phir(3),phii(3)
  COMPLEX(KIND=DP) dcmplx
  CHARACTER(LEN=75) :: line
  CHARACTER(LEN=3)  :: atm1
  LOGICAL :: first
  DATA first/.TRUE./
  SAVE first
  !
  READ(1,*) 
  READ(1,*) 
  IF (first) THEN
     !
     ! read cell information from file
     !
     READ(1,*) ntyp,nat,ibrav,(celldm(i),i=1,6)
     IF (nat.GT.nax) CALL errore('read_f','nax too small',nat)
     IF (ntyp.GT.nat) CALL errore('read_f','ntyp.gt.nat!!',ntyp)
     DO nt = 1,ntyp
        READ(1,*) i,atm(nt),amass(nt)
        IF (i.NE.nt) CALL errore('read_f','wrong data read',nt)
     END DO
     DO na=1,nat
        READ(1,*) i,ityp(na),(tau(j,na),j=1,3)
        IF (i.NE.na) CALL errore('read_f','wrong data read',na)
     END DO
     !
     first=.FALSE.
     lrigid=.FALSE.
     !
  ELSE
     !
     ! check cell information with previous one
     !
     READ(1,*) ntyp1,nat1,ibrav1,(celldm1(i),i=1,6)
     IF (ntyp1.NE.ntyp) CALL errore('read_f','wrong ntyp',1)
     IF (nat1.NE.nat) CALL errore('read_f','wrong nat',1)
     IF (ibrav1.NE.ibrav) CALL errore('read_f','wrong ibrav',1)
     DO i=1,6
        IF(celldm1(i).NE.celldm(i)) CALL errore('read_f','wrong celldm',i)
     END DO
     DO nt = 1,ntyp
        READ(1,*) i,atm1,amass1
        IF (i.NE.nt) CALL errore('read_f','wrong data read',nt)
        IF (atm1.NE.atm(nt)) CALL errore('read_f','wrong atm',nt)
        IF (amass1.NE.amass(nt)) CALL errore('read_f','wrong amass',nt)
     END DO
     DO na=1,nat
        READ(1,*) i,ityp1,(tau1(j),j=1,3)
        IF (i.NE.na) CALL errore('read_f','wrong data read',na)
        IF (ityp1.NE.ityp(na)) CALL errore('read_f','wrong ityp',na)
        IF (tau1(1).NE.tau(1,na)) CALL errore('read_f','wrong tau1',na)
        IF (tau1(2).NE.tau(2,na)) CALL errore('read_f','wrong tau2',na)
        IF (tau1(3).NE.tau(3,na)) CALL errore('read_f','wrong tau3',na)
     END DO
  END IF
  !
  !
  nqs = 0
100 CONTINUE
  READ(1,*)
  READ(1,'(a)') line
  IF (line(6:14).NE.'Dynamical') THEN
     IF (nqs.EQ.0) CALL errore('read',' stop with nqs=0 !!',1)
     q2 = xq(1,nqs)**2 + xq(2,nqs)**2 + xq(3,nqs)**2
     IF (q2.NE.0.d0) RETURN
     DO WHILE (line(6:15).NE.'Dielectric') 
        READ(1,'(a)',err=200, END=200) line
     END DO
     lrigid=.TRUE.
     READ(1,*) ((epsil(i,j),j=1,3),i=1,3)
     READ(1,*)
     READ(1,*)
     READ(1,*)
     WRITE (*,*) 'macroscopic fields =',lrigid
     WRITE (*,'(3f10.5)') ((epsil(i,j),j=1,3),i=1,3)
     DO na=1,nat
        READ(1,*)
        READ(1,*) ((zeu(i,j,na),j=1,3),i=1,3)
        WRITE (*,*) ' na= ', na
        WRITE (*,'(3f10.5)') ((zeu(i,j,na),j=1,3),i=1,3)
     END DO
     RETURN
200  WRITE (*,*) ' Dielectric Tensor not found'
     lrigid=.FALSE.     
     RETURN
  END IF
  !
  nqs = nqs + 1
  READ(1,*) 
  READ(1,'(a)') line
  READ(line(11:75),*) (xq(i,nqs),i=1,3)
  READ(1,*) 
  !
  DO na=1,nat
     DO nb=1,nat
        READ(1,*) i,j
        IF (i.NE.na) CALL errore('read_f','wrong na read',na)
        IF (j.NE.nb) CALL errore('read_f','wrong nb read',nb)
        DO i=1,3
           READ (1,*) (phir(j),phii(j),j=1,3)
           DO j = 1,3
              phi(i,j,na,nb,nqs) = dcmplx(phir(j),phii(j))
           END DO
        END DO
     END DO
  END DO
  !
  go to 100
  !
END SUBROUTINE read_file
!
!----------------------------------------------------------------------------
SUBROUTINE trasl( phi, phiq, nq, nrx1, nrx2, nrx3, nat, m1, m2, m3, nax )
  !----------------------------------------------------------------------------
  !
  USE kinds, ONLY : DP
  !
  IMPLICIT NONE
  INTEGER:: j1,j2, m1, m2, m3, nrx1, nrx2, nrx3, na1, na2, nat, nax, nq
  !
  COMPLEX(KIND=DP) :: phi(nrx1,nrx2,nrx3,3,3,nax,nax)
  COMPLEX(KIND=DP) :: phiq(3,3,nax,nax,48)
  !
  DO j1=1,3
     DO j2=1,3
        DO na1=1,nat
           DO na2=1,nat
              phi(m1,m2,m3,j1,j2,na1,na2) = &
                   0.5 * (      phiq(j1,j2,na1,na2,nq) +  &
                          CONJG(phiq(j2,j1,na2,na1,nq)))
           END DO
        END DO
     END DO
  END DO
  !
  RETURN 
END SUBROUTINE trasl

# if defined __AIX || defined __FFTW || defined __SGI
#  define __FFT_MODULE_DRV
# endif

!----------------------------------------------------------------------------
SUBROUTINE tolerant_cft3( f, nr1, nr2, nr3, nrx1, nrx2, nrx3, iflg )
  !----------------------------------------------------------------------------
  !
  !  cft3 called for vectors with arbitrary maximal dimensions
  !
  USE kinds,      ONLY : DP
#if defined __FFT_MODULE_DRV
  USE fft_scalar, ONLY : cfft3d
#endif

  IMPLICIT NONE
  INTEGER :: nr1,nr2,nr3,nrx1,nrx2,nrx3,iflg
  COMPLEX(KIND=DP) :: f(nrx1,nrx2,nrx3)
#if defined __FFT_MODULE_DRV
  COMPLEX(KIND=DP) :: ftmp(nrx1*nrx2*nrx3)
#endif
  INTEGER :: i0,i1,i2,i3
  !
  i0=0
  DO i3=1,nr3
     DO i2=1,nr2
        DO i1=1,nr1
           i0 = i0 + 1
#if defined __FFT_MODULE_DRV
           ftmp(i0) = f(i1,i2,i3)
#else
           f(i0,1,1) = f(i1,i2,i3)
#endif
        ENDDO
     ENDDO
  ENDDO
  !
  IF (nr1.NE.1 .OR. nr2.NE.1 .OR. nr3.NE.1)                         &
#if defined __FFT_MODULE_DRV
       &    CALL cfft3d(ftmp,nr1,nr2,nr3,nr1,nr2,nr3,1)
#else
       &    CALL cft_3(f,nr1,nr2,nr3,nr1,nr2,nr3,1,iflg)
#endif
  !
  i0=nr1*nr2*nr3
  DO i3=nr3,1,-1
     DO i2=nr2,1,-1
        DO i1=nr1,1,-1
#if defined __FFT_MODULE_DRV
           f(i1,i2,i3) = ftmp(i0)
#else
           f(i1,i2,i3) = f(i0,1,1)
#endif
           i0 = i0 - 1
        ENDDO
     ENDDO
  ENDDO

  !
  RETURN
END SUBROUTINE tolerant_cft3
