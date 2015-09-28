!
! Copyright (C) 2006 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
PROGRAM plotproj
!
!  This small program is used to select the band eigenvalues whose
!  wavefunctions projected on atomic wavefunctions have projections larger
!  than a given threshold. It requires two input files. The first is a
!  file with the band eigenvalues, written in the output of pw.x.
!  The input file with the bands has the following format:
!  nbnd, nks     ! number of bands, number of k points
!  --- blank line
!  kvector coordinates
!  --- blank line
!  bands eigenvalues
!  ...
!  --- blank line
!  kvector coordinates
!  --- blank line
!  bands eigenvalues
!  ...
!
!  The second file is written by the projwfc.x program with the option
!  lsym=.false.
!
!  The input of this program is:
!  filename     ! name of the file with the band eigenvalues
!  filename1    ! name of the file with the projections
!  fileout      ! name of the output file where the bands are written
!  threshold    ! see below
!  ncri         ! number of criterions for selecting the bands
!  for each criterion
!  first_atomic_wfc, last_atomic_wfc   ! the band is selected if the
!                                        sum of the projections on
!                                        the atomic wavefunctions between
!                                        first_atomic_wfc and
!                                        last_atomic_wfc is larger than
!                                        threshold. The sum is done on
!                                        all criterions.
!

  IMPLICIT NONE
  INTEGER, PARAMETER :: DP = selected_real_kind(14,200)
  REAL(DP), ALLOCATABLE :: e(:,:), k(:,:), kx(:)
  INTEGER :: nks = 0, nbnd = 0, ios, n, i, ibnd, na, idum, nat, &
       natomwfc, nwfc, ntyp, ncri, icri
  LOGICAL, ALLOCATABLE :: toplot(:,:)
  CHARACTER(len=256) :: filename, filename1
  REAL(DP) :: psum, threshold
  REAL(DP), ALLOCATABLE :: proj(:,:,:)
  INTEGER, ALLOCATABLE :: first_atomic_wfc(:), last_atomic_wfc(:)

  CALL get_file ( filename )

  OPEN(UNIT=1,FILE=filename,FORM='formatted',status='old',iostat=ios)
  IF (ios/=0) STOP 'Error opening band file '

  READ(1,*, err=20, iostat=ios) nbnd, nks

  IF (nks <= 0 .or. nbnd <= 0 ) THEN
     STOP 'Error reading file header'
  ELSE
     PRINT '("Reading ",i4," bands at ",i4," k-points")', nbnd, nks
  ENDIF

  ALLOCATE (e(nbnd,nks))
  ALLOCATE (k(3,nks))
  ALLOCATE (kx(nks))
  ALLOCATE (toplot(nbnd,nks))

  DO n=1,nks
     READ(1, *, ERR=20, IOSTAT=ios)
     READ(1, '(13x,3f7.4)', ERR=20, IOSTAT=ios) (k(i,n), i=1,3)
     READ(1, *, ERR=20, IOSTAT=ios)
     READ(1, '(2x,8f9.4)', END=20, ERR=20) (e(i,n),i=1,nbnd)
     IF (n==1) THEN
        kx(n) = sqrt (k(1,1)**2 + k(2,1)**2 + k(3,1)**2)
     ELSE
        kx(n) = kx(n-1) + sqrt ( (k(1,n)-k(1,n-1))**2 + &
             (k(2,n)-k(2,n-1))**2 + &
             (k(3,n)-k(3,n-1))**2 )
     ENDIF
  ENDDO

20 IF (ios/=0) STOP "problem reading files"
  CLOSE(UNIT=1)

  CALL get_file ( filename1 )
  OPEN(UNIT=1, FILE=filename1, FORM='formatted', STATUS='old', IOSTAT=ios)
  IF (ios/=0) STOP 'Error opening projection file '
  READ(1, *, ERR=20, IOSTAT=ios)
  READ (1, '(8i8)', ERR=20, IOSTAT=ios) idum, idum, idum, idum, idum, &
       idum, nat, ntyp
  DO i=1,2+nat+ntyp
     READ(1, *, ERR=20, IOSTAT=ios)
  ENDDO
  READ (1, '(3i8)',ERR=20, IOSTAT=ios) natomwfc, nks, nbnd
  READ (1, *, ERR=20, IOSTAT=ios)

  ALLOCATE( proj(natomwfc,nbnd,nks) )
  DO nwfc = 1, natomwfc
     READ(1, *, ERR=20, IOSTAT=ios)
     DO n=1,nks
        DO ibnd=1,nbnd
           READ(1, '(2i8,f20.10)', ERR=20, IOSTAT=ios) idum,idum,proj(nwfc,ibnd,n)
        ENDDO
     ENDDO
  ENDDO
  CLOSE(1)

  WRITE(*,'("output file > ")', advance="NO")
  READ(5,'(a)', END=25, ERR=25)  filename

  IF (filename == ' ' ) THEN
     PRINT '("skipping ...")'
     GOTO 25
  ENDIF

  OPEN (UNIT=2,FILE=filename,FORM='formatted',STATUS='unknown',IOSTAT=ios)
  IF (ios/=0) STOP "Error opening output file "

  READ(5, *, ERR=20, IOSTAT=ios) threshold
  READ(5, *, ERR=20, IOSTAT=ios) ncri
  IF (ncri<1)  STOP '("no orbital given ...")'
  ALLOCATE(first_atomic_wfc(ncri))
  ALLOCATE(last_atomic_wfc(ncri))
  DO icri=1,ncri
     READ(5, *, ERR=20, IOSTAT=ios) first_atomic_wfc(icri),  &
          last_atomic_wfc(icri)
     IF (first_atomic_wfc(icri)>natomwfc.or.last_atomic_wfc(icri)>natomwfc .or. &
          first_atomic_wfc(icri)<1 .or. &
          last_atomic_wfc(icri)<first_atomic_wfc(icri) ) THEN
        PRINT '("Problem with ...",i5)', icri
        GOTO 25
     ENDIF
  ENDDO

  toplot=.false.
  DO i=1,nbnd
     DO n=1,nks
        psum=0.d0
        DO icri=1,ncri
           DO nwfc=first_atomic_wfc(icri),last_atomic_wfc(icri)
              psum=psum+abs(proj(nwfc,i,n))
           ENDDO
        ENDDO
        toplot(i,n)=toplot(i,n).or.(psum > threshold)
     ENDDO
  ENDDO

  DO i=1,nbnd
     DO n=1,nks
        IF (toplot(i,n))  WRITE (2,'(2f10.4)') kx(n), e(i,n)
     ENDDO
  ENDDO

  CLOSE (UNIT = 2)
25 CONTINUE

END PROGRAM plotproj
