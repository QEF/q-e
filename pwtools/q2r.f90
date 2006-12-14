!
! Copyright (C) 2001-2005 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "f_defs.h"
!
Module dynamicalq
  !
  ! All variables read from file that need dynamical allocation
  !
  USE kinds, ONLY: DP
  COMPLEX(DP), ALLOCATABLE :: phiq(:,:,:,:,:) 
  REAL(DP), ALLOCATABLE ::  tau(:,:), zeu(:,:,:)
  INTEGER, ALLOCATABLE ::  ityp(:)
  !
end Module dynamicalq
!
!----------------------------------------------------------------------------
PROGRAM q2r
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
  !                   "fildyn"1-N contain force constans C_n = C(q_n)
  !                   for n=1,...N, where N is the number of q-points
  !                   in the irreducible brillouin zone
  !                   Normally this should be the same as specified
  !                   on input to the phonon code
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
  USE kinds,      ONLY : DP
  USE mp,         ONLY : mp_start, mp_env, mp_end, mp_barrier
  USE mp_global,  ONLY : nproc, mpime, mp_global_start
  USE dynamicalq, ONLY : phiq, tau, ityp, zeu
  USE fft_scalar, ONLY : cfft3d
  !
  IMPLICIT NONE
  !
  INTEGER,       PARAMETER :: ntypx = 10
  REAL(DP), PARAMETER :: eps=1.D-5, eps12=1.d-12
  INTEGER                  :: nr1, nr2, nr3, nr(3)
  !     dimensions of the FFT grid formed by the q-point grid
  !
  CHARACTER(len=20)  :: crystal
  CHARACTER(len=256) :: fildyn, filin, filj, filf, flfrc
  CHARACTER(len=3)   :: atm(ntypx)
  CHARACTER(len=9)   :: symm_type
  CHARACTER(LEN=6), EXTERNAL :: int_to_char
  !
  LOGICAL :: lq, lrigid, lrigid_save, lnogridinfo
  CHARACTER (LEN=10) :: zasr
  INTEGER :: m1, m2, m3, m(3), l1, l2, l3, i, j, j1, j2, na1, na2, ipol, nn
  INTEGER :: nat, nq, ntyp, iq, icar, nfile, ifile, nqs, nq_log
  INTEGER :: na, nt
  !
  INTEGER :: gid, ibrav, ierr
  !
  INTEGER, ALLOCATABLE ::  nc(:,:,:)
  COMPLEX(DP), ALLOCATABLE :: phid(:,:,:,:,:)
  !
  REAL(DP) :: celldm(6), at(3,3), bg(3,3)
  REAL(DP) :: q(3,48),omega, xq, amass(ntypx), resi
  REAL(DP) :: epsil(3,3)
  !
  logical           :: la2F
  !
  NAMELIST / input / fildyn, flfrc, zasr, la2F
  !
  CALL mp_start()
  !
  CALL mp_env( nproc, mpime, gid )
  !
  IF ( mpime == 0 ) THEN
     !
     ! ... all calculations are done by the first cpu
     !
     fildyn = ' '
     flfrc = ' '
     !
     la2F=.false.
     !
     CALL input_from_file ( ) 
     !
     READ ( 5, input )
     !
     ! check input
     !
     IF (flfrc == ' ')  CALL errore ('q2r',' bad flfrc',1)
     !
     OPEN (unit=1, file=TRIM(fildyn)//'0', status='old', form='formatted', &
          iostat=ierr)
     lnogridinfo = ( ierr /= 0 )
     IF (lnogridinfo) THEN
        WRITE (6,*) ' file ',TRIM(fildyn)//'0', ' not found'
        WRITE (6,*) ' reading grid info from input'
        READ (5, *) nr1, nr2, nr3
        READ (5, *) nfile
     ELSE
        WRITE (6,*) ' reading grid info from file ',TRIM(fildyn)//'0'
        READ (1, *) nr1, nr2, nr3
        READ (1, *) nfile
        CLOSE (unit=1, status='keep')
     END IF
     !
     IF (nr1 < 1 .OR. nr1 > 1024) CALL errore ('q2r',' nr1 wrong or missing',1)
     IF (nr2 < 1 .OR. nr2 > 1024) CALL errore ('q2r',' nr2 wrong or missing',1)
     IF (nr3 < 1 .OR. nr2 > 1024) CALL errore ('q2r',' nr3 wrong or missing',1)
     IF (nfile < 1 .OR. nfile > 1024) &
        CALL errore ('q2r','too few or too many file',MAX(1,nfile))
     !
     ! copy nrX -> nr(X)
     !
     nr(1) = nr1
     nr(2) = nr2
     nr(3) = nr3
     !
     ! D matrix (analytical part)
     !
     ntyp = ntypx ! avoids spurious out-of-bound errors
     !
     ALLOCATE ( nc(nr1,nr2,nr3) )
     nc = 0
     !
     ! Force constants in reciprocal space read from file
     !
     DO ifile=1,nfile
        IF (lnogridinfo) THEN
           READ(5,'(a)') filin
        ELSE
           filin = TRIM(fildyn) // TRIM( int_to_char( ifile ) )
        END IF
        WRITE (6,*) ' reading force constants from file ',TRIM(filin)
        OPEN (unit=1, file=filin, status='old', form='formatted', iostat=ierr)
        IF (ierr /= 0) CALL errore('q2r','file '//TRIM(filin)//' missing!',1)
        CALL read_file (nqs, q, epsil, lrigid,  &
             ntyp, nat, ibrav, symm_type, celldm, at, atm, amass)
        IF (ifile == 1) THEN
           ! it must be allocated here because nat is read from file
           ALLOCATE (phid(nr1*nr2*nr3,3,3,nat,nat) )
           !
           lrigid_save=lrigid

           CALL latgen(ibrav,celldm,at(1,1),at(1,2),at(1,3),omega)
           at = at / celldm(1)  !  bring at in units of alat 

           CALL volume(celldm(1),at(1,1),at(1,2),at(1,3),omega)
           CALL recips(at(1,1),at(1,2),at(1,3),bg(1,1),bg(1,2),bg(1,3))
           IF (lrigid .AND. (zasr.NE.'no')) THEN
              CALL set_zasr ( zasr, nr1,nr2,nr3, nat, ibrav, tau, zeu)
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
              xq = 0.0d0
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
              IF (lrigid) THEN
                 CALL rgd_blk (nr1,nr2,nr3,nat,phiq(1,1,1,1,nq),q(1,nq), &
                  tau,epsil,zeu,bg,omega,-1.d0)
              END IF
              CALL trasl ( phid, phiq, nq, nr1,nr2,nr3, nat, m(1),m(2),m(3))
           ELSE
              WRITE (*,'(3i4)') (m(i),i=1,3)
              CALL errore('init',' nc already filled: wrong q grid or wrong nr',1)
           END IF
        END DO
     END DO
     !
     ! Check grid dimension
     !
     nq_log = SUM (nc)
     IF (nq_log == nr1*nr2*nr3) THEN
        WRITE (6,'(/5x,a,i4)') ' q-space grid ok, #points = ',nq_log
     ELSE
        CALL errore('init',' missing q-point(s)!',1)
     END IF
     !
     ! dyn.mat. FFT (use serial version)
     !
     DO j1=1,3
        DO j2=1,3
           DO na1=1,nat
              DO na2=1,nat
                 CALL cfft3d ( phid (:,j1,j2,na1,na2), &
                      nr1,nr2,nr3, nr1,nr2,nr3, 1 )
                 phid(:,j1,j2,na1,na2) = &
                      phid(:,j1,j2,na1,na2) / DBLE(nr1*nr2*nr3)
              END DO
           END DO
        END DO
     END DO
     !
     ! Real space force constants written to file (analytical part)
     !
     OPEN(unit=2,file=flfrc,status='unknown',form='formatted')
     WRITE(2,'(i3,i5,i3,6f11.7)') ntyp,nat,ibrav,celldm
     if (ibrav==0) then
        write (2,'(a)') symm_type
        write (2,'(2x,3f15.9)') ((at(i,j),i=1,3),j=1,3)
     end if
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
                 WRITE (2,'(4i4)') j1,j2,na1,na2
                 nn=0
                 DO m3=1,nr3
                    DO m2=1,nr2
                       DO m1=1,nr1
                          nn=nn+1
                          WRITE (2,'(3i4,2x,1pe18.11)')   &
                               m1,m2,m3, DBLE(phid(nn,j1,j2,na1,na2))
                       END DO
                    END DO
                 END DO
              END DO
           END DO
        END DO
     END DO
     CLOSE(2)
     resi = SUM ( ABS (AIMAG ( phid ) ) )
     IF (resi > eps12) THEN
        WRITE (6,"(/5x,' fft-check warning: sum of imaginary terms = ',e12.7)") resi
     ELSE
        WRITE (6,"(/5x,' fft-check success (sum of imaginary terms < 10^-12)')")
     END IF
     !
     DEALLOCATE(phid, phiq, zeu, nc)
     !
     IF(la2F) CALL gammaq2r ( nfile, nat, nr1, nr2, nr3, at )
     !
     DEALLOCATE (tau, ityp)
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
SUBROUTINE gammaq2r( nqtot, nat, nr1, nr2, nr3, at )
  !----------------------------------------------------------------------------
  !
  USE kinds, ONLY : DP
  USE fft_scalar, ONLY : cfft3d
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
  write (6,*)
  write (6,*) '  Preparing gamma for a2F '
  write (6,*)
  !
  nr(1) = nr1
  nr(2) = nr2
  nr(3) = nr3
  !
  DO isig=1, nsig
     filea2F = 50 + isig
     write(name,"(A7,I2)") 'a2Fq2r.',filea2F
     open(filea2F, file=name, STATUS = 'old', FORM = 'formatted')
     nc = 0
     !
     ! to pass to matdyn, for each isig, we read: degauss, Fermi energy and DOS
     !
     DO count_q=1,nqtot
        !
        READ(filea2F,*) deg, ef, dosscf
        READ(filea2F,*) nstar
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
        write (6,*)
        write (6,'(" Broadening = ",F10.3)') deg
        write (6,'(5x,a,i4)') ' q-space grid ok, #points = ',nq_log
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
     close(filea2F)
     !
     filea2F = 60 + isig
     write(name,"(A10,I2)") 'a2Fmatdyn.',filea2F
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

     resi = SUM ( ABS ( AIMAG( gamout ) ) )
     
     IF (resi > eps12) THEN
        WRITE (6,"(/5x,' fft-check warning: sum of imaginary terms = ',e12.7)") resi
     ELSE
        WRITE (6,"(/5x,' fft-check success (sum of imaginary terms < 10^-12)')")
     END IF

  ENDDO
  !
  DEALLOCATE (gaminp, gamout )
  !
END SUBROUTINE gammaq2r
!
!----------------------------------------------------------------------------
SUBROUTINE read_file( nqs, xq, epsil, lrigid, &
                      ntyp, nat, ibrav, symm_type, celldm, at, atm, amass )
  !----------------------------------------------------------------------------
  !
  USE kinds, ONLY : DP
  USE dynamicalq, ONLY: phiq, tau, ityp, zeu
  !
  IMPLICIT NONE
  !
  REAL(DP), PARAMETER :: eps8=1.D-8
  ! I/O variables
  LOGICAL :: lrigid
  INTEGER :: nqs, ntyp, nat, ibrav
  REAL(DP) :: epsil(3,3)
  REAL(DP) :: xq(3,48), celldm(6), at(3,3), amass(ntyp)
  CHARACTER(LEN=3) atm(ntyp)
  CHARACTER(LEN=9) symm_type
  ! local variables
  INTEGER :: ntyp1,nat1,ibrav1,ityp1
  INTEGER :: i, j, na, nb, nt
  REAL(DP) :: tau1(3), amass1, at1(3,3), celldm1(6), q2
  REAL(DP) :: phir(3),phii(3)
  CHARACTER(LEN=75) :: line
  CHARACTER(LEN=3)  :: atm1
  CHARACTER(LEN=9) symm_type1
  LOGICAL, SAVE :: first =.TRUE.
  !
  READ(1,*) 
  READ(1,*) 
  IF (first) THEN
     !
     ! read cell information from file
     !
     READ(1,*) ntyp,nat,ibrav,(celldm(i),i=1,6)
     if (ibrav==0) then
        read (1,'(a)') symm_type
        read (1,*) ((at(i,j),i=1,3),j=1,3)
     end if

     IF (ntyp.GT.nat) CALL errore('read_f','ntyp.gt.nat!!',ntyp)
     DO nt = 1,ntyp
        READ(1,*) i,atm(nt),amass(nt)
        IF (i.NE.nt) CALL errore('read_f','wrong data read',nt)
     END DO
     ALLOCATE ( ityp(nat), tau(3,nat) )
     DO na=1,nat
        READ(1,*) i,ityp(na),(tau(j,na),j=1,3)
        IF (i.NE.na) CALL errore('read_f','wrong data read',na)
     END DO
     !
     ALLOCATE ( phiq (3,3,nat,nat,48), zeu (3,3,nat) )
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
        IF( abs (celldm1(i)-celldm(i)) > eps8 ) &
             CALL errore('read_f','wrong celldm',i)
     END DO
     if (ibrav==0) then
         read (1,*) symm_type1
         if (symm_type1 /= symm_type) &
            CALL errore('read_f','wrong symm_type for ibrav=0',1)
         read (1,*) ((at1(i,j),i=1,3),j=1,3)
         do i=1,3
            do j=1,3
               if( abs (at1(i,j)-at(i,j)) > eps8) &
                 CALL errore('read_f','wrong at(i,j)',i+3*(j-1))
            end do
         end do
     end if
     DO nt = 1,ntyp
        READ(1,*) i,atm1,amass1
        IF (i.NE.nt) CALL errore('read_f','wrong data read',nt)
        IF (atm1.NE.atm(nt)) CALL errore('read_f','wrong atm',nt)
        IF (abs(amass1-amass(nt)) > eps8 ) &
             CALL errore('read_f','wrong amass',nt)
     END DO
     DO na=1,nat
        READ(1,*) i,ityp1,(tau1(j),j=1,3)
        IF (i.NE.na) CALL errore('read_f','wrong data read',na)
        IF (ityp1.NE.ityp(na)) CALL errore('read_f','wrong ityp',na)
        IF ( abs (tau1(1)-tau(1,na)) > eps8 .OR. &
             abs (tau1(2)-tau(2,na)) > eps8 .OR. &
             abs (tau1(3)-tau(3,na)) > eps8 ) &
             CALL errore('read_f','wrong tau',na)
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
              phiq (i,j,na,nb,nqs) = CMPLX (phir(j),phii(j))
           END DO
        END DO
     END DO
  END DO
  !
  go to 100
  !
END SUBROUTINE read_file
!
!-----------------------------------------------------------------------
subroutine read_gamma (nqs, nat, ifn, xq, gaminp)
  !-----------------------------------------------------------------------
  !
  USE kinds, ONLY : DP
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
     READ(ifn,*)
     READ(ifn,*) 
     READ(ifn,*)
     READ(ifn,'(11X,3F14.9)')  (xq(i,iq),i=1,3)
     !     write(*,*) 'xq    ',iq,(xq(i,iq),i=1,3)
     READ(ifn,*)
     do na=1,nat
        do nb=1,nat
           read(ifn,*) i,j
           if (i.ne.na) call errore('read_f','wrong na read',na)
           if (j.ne.nb) call errore('read_f','wrong nb read',nb)
           do i=1,3
              read (ifn,*) (phir(j),phii(j),j=1,3)
              do j = 1,3
                 gaminp(i,j,na,nb,iq) = CMPLX (phir(j),phii(j))
              end do
              !           write(*,*) 'gaminp  ',(gaminp(i,j,na,nb,iq),j=1,3)
           end do
        end do
     end do
     !
  Enddo
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
  implicit none
  character(len=10) :: zasr
  integer ibrav,nr1,nr2,nr3,nr,m,p,k,l,q,r
  integer n,i,j,n1,n2,n3,na,nb,nat,axis,i1,j1,na1
  !
  real(8) sum, zeu(3,3,nat)
  real(8) tau(3,nat), zeu_new(3,3,nat)
  !
  real(8) zeu_u(6*3,3,3,nat)
  ! These are the "vectors" associated with the sum rules on effective charges
  !
  integer zeu_less(6*3),nzeu_less,izeu_less
  ! indices of vectors zeu_u that are not independent to the preceding ones,
  ! nzeu_less = number of such vectors, izeu_less = temporary parameter
  !
  real(8) zeu_w(3,3,nat), zeu_x(3,3,nat),scal,norm2
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
        write(6,*) 'zasr: rotational axis may be wrong'
     endif
     write(6,'("zasr rotation axis in 1D system= ",I4)') axis
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
      write(6,'("Norm of the difference between old and new effective ", &
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
subroutine sp_zeu(zeu_u,zeu_v,nat,scal)
  !-----------------------------------------------------------------------
  !
  ! does the scalar product of two effective charges matrices zeu_u and zeu_v 
  ! (considered as vectors in the R^(3*3*nat) space, and coded in the usual way)
  !
  implicit none
  integer i,j,na,nat
  real(8) zeu_u(3,3,nat)
  real(8) zeu_v(3,3,nat)
  real(8) scal  
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
