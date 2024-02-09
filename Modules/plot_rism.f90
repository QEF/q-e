!
! Copyright (C) 2016 National Institute of Advanced Industrial Science and Technology (AIST)
! [ This code is written by Satomichi Nishihara. ]
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------------------------------------------
SUBROUTINE plot_rism(filplot, title, nsite, nr1x, nr2x, nr3x, nr1, nr2, nr3, &
                   & nat, ntyp, ibrav, celldm, at, gcutm, laue, &
                   & asite, atm, ityp, zv, tau, plot, iflag)
  !--------------------------------------------------------------------------
  !
  ! ... iflag >0 : write 3D-RISM's data
  ! ... iflag< 0 : read 3D-RISM's data
  !
  USE io_global, ONLY : stdout
  USE kinds,     ONLY : DP
  !
  IMPLICIT NONE
  !
  CHARACTER(LEN=*),  INTENT(IN)    :: filplot
  CHARACTER(LEN=75), INTENT(INOUT) :: title
  INTEGER,           INTENT(INOUT) :: nsite
  INTEGER,           INTENT(INOUT) :: nr1x
  INTEGER,           INTENT(INOUT) :: nr2x
  INTEGER,           INTENT(INOUT) :: nr3x
  INTEGER,           INTENT(INOUT) :: nr1
  INTEGER,           INTENT(INOUT) :: nr2
  INTEGER,           INTENT(INOUT) :: nr3
  INTEGER,           INTENT(INOUT) :: nat
  INTEGER,           INTENT(INOUT) :: ntyp
  INTEGER,           INTENT(INOUT) :: ibrav
  REAL(DP),          INTENT(INOUT) :: celldm(6)
  REAL(DP),          INTENT(INOUT) :: at(3, 3)
  REAL(DP),          INTENT(INOUT) :: gcutm
  LOGICAL,           INTENT(INOUT) :: laue
  CHARACTER(LEN=12), INTENT(INOUT) :: asite(1:*)
  CHARACTER(LEN=*),  INTENT(INOUT) :: atm(1:*)
  INTEGER,           INTENT(INOUT) :: ityp(1:*)
  REAL(DP),          INTENT(INOUT) :: zv(1:*)
  REAL(DP),          INTENT(INOUT) :: tau(3, 1:*)
  REAL(DP),          INTENT(INOUT) :: plot(1:*)
  INTEGER,           INTENT(INOUT) :: iflag
  !
  INTEGER           :: i
  INTEGER           :: ipol
  INTEGER           :: it
  INTEGER           :: ia
  INTEGER           :: is
  INTEGER           :: ir
  INTEGER           :: nfft
  INTEGER           :: ndum
  INTEGER           :: ios
  INTEGER           :: iunplot
  INTEGER, EXTERNAL :: find_free_unit
  !
  ! ... check parameters
  IF (TRIM(filplot) == '') THEN
    CALL errore('plot_rism', 'filename missing', 1)
  END IF
  !
  IF (iflag == 0) THEN
    CALL errore('plot_rism', 'iflag==0 not allowed', 1)
  END IF
  !
  ! ... open file
  iunplot = find_free_unit()
  !
  IF (iflag > 0) THEN
     WRITE(stdout, '(5X,"Writing data to file  ",A)') TRIM(ADJUSTL(filplot))
     OPEN(unit=iunplot, file=TRIM(ADJUSTL(filplot)), form='formatted', status='unknown', iostat=ios)
  ELSE
     WRITE(stdout, '(5X,"Reading data from file  ",A)') TRIM(ADJUSTL(filplot))
     OPEN(unit=iunplot, file=TRIM(ADJUSTL(filplot)), form='formatted', status='old', iostat=ios)
  END IF
  !
  IF (ios /= 0) THEN
    CALL errore('plot_rism', 'opening file ' // TRIM(ADJUSTL(filplot)), ABS(ios))
  END IF
  !
  REWIND(iunplot)
  !
  IF (iflag > 0) THEN
    ! ... write file
    WRITE(iunplot, '(A)') title
    WRITE(iunplot, '(9I8)') nsite, nr1x, nr2x, nr3x, nr1, nr2, nr3, nat, ntyp
    WRITE(iunplot, '(I6,2X,6F16.8)') ibrav, celldm
    IF (ibrav == 0) THEN
      DO i = 1, 3
        WRITE(iunplot, '(3E22.12e3)') (at(ipol, i), ipol = 1, 3)
      END DO
    END IF
    WRITE(iunplot, '(F20.10,L3)')          gcutm, laue
    WRITE(iunplot, '(I4,3X,A12)')          (is, asite(is), is = 1, nsite)
    WRITE(iunplot, '(I4,3X,A2,3X,F5.2)')   (it, atm(it), zv(it), it = 1, ntyp)
    WRITE(iunplot, '(I4,3X,3F15.9,3X,I2)') (ia, (tau(ipol, ia), ipol = 1, 3), ityp(ia), ia = 1, nat)
    !
    nfft = nr1x * nr2x * nr3x
    DO is = 1, nsite
       WRITE(iunplot, '(5(1PE18.9e3))')    (plot(ir), ir = (is - 1) * nfft + 1, is * nfft)
    END DO
    !
  ELSE
    ! ... read file
    READ(iunplot, '(A)') title
    READ(iunplot, *) nsite, nr1x, nr2x, nr3x, nr1, nr2, nr3, nat, ntyp
    READ(iunplot, *) ibrav, celldm
    IF (ibrav == 0) THEN
      DO i = 1, 3
        READ(iunplot, *) (at(ipol, i), ipol = 1, 3)
      END DO
    END IF
    READ(iunplot, *)                    gcutm, laue
    READ(iunplot, '(I4,3X,A12)')        (ndum, asite(is), is = 1, nsite)
    READ(iunplot, '(I4,3X,A2,3X,F5.2)') (ndum, atm(it), zv(it), it = 1, ntyp)
    READ(iunplot, *)                    (ndum, (tau(ipol, ia), ipol = 1, 3), ityp(ia), ia = 1, nat)
    !
    nfft = nr1x * nr2x * nr3x
    DO is = 1, nsite
       READ(iunplot, *)                 (plot(ir), ir = (is - 1) * nfft + 1, is * nfft)
    END DO
    !
  END IF
  !
  ! ... close file
  CLOSE(unit=iunplot)
  !
END SUBROUTINE plot_rism
!
!--------------------------------------------------------------------------
SUBROUTINE read_rism_header(filplot, title, nsite, nr1x, nr2x, nr3x, nr1, nr2, nr3, &
                          & nat, ntyp, ibrav, celldm, at, gcutm, laue)
  !--------------------------------------------------------------------------
  !
  ! ... read header of file "filplot"
  !
  USE io_global, ONLY : stdout
  USE kinds,     ONLY : DP
  !
  IMPLICIT NONE
  !
  CHARACTER(LEN=*),  INTENT(IN)  :: filplot
  CHARACTER(LEN=75), INTENT(OUT) :: title
  INTEGER,           INTENT(OUT) :: nsite
  INTEGER,           INTENT(OUT) :: nr1x
  INTEGER,           INTENT(OUT) :: nr2x
  INTEGER,           INTENT(OUT) :: nr3x
  INTEGER,           INTENT(OUT) :: nr1
  INTEGER,           INTENT(OUT) :: nr2
  INTEGER,           INTENT(OUT) :: nr3
  INTEGER,           INTENT(OUT) :: nat
  INTEGER,           INTENT(OUT) :: ntyp
  INTEGER,           INTENT(OUT) :: ibrav
  REAL(DP),          INTENT(OUT) :: celldm(6)
  REAL(DP),          INTENT(OUT) :: at(3, 3)
  REAL(DP),          INTENT(OUT) :: gcutm
  LOGICAL,           INTENT(OUT) :: laue
  !
  INTEGER           :: i
  INTEGER           :: ipol
  INTEGER           :: ios
  INTEGER           :: iunplot
  INTEGER, EXTERNAL :: find_free_unit
  !
  ! ... check parameters
  IF (TRIM(filplot) == '') THEN
    CALL errore('read_rism_h', 'filename missing', 1)
  END IF
  !
  ! ... open file
  iunplot = find_free_unit()
  !
  WRITE(stdout, '(5X,"Reading header from file  ",A)') TRIM(ADJUSTL(filplot))
  OPEN(unit=iunplot, file=TRIM(ADJUSTL(filplot)), form='formatted', status='old', iostat=ios)
  !
  IF (ios /= 0) THEN
    CALL errore ('read_rism_h', 'opening file ' // TRIM(ADJUSTL(filplot)), ABS(ios))
  END IF
  !
  REWIND(iunplot)
  !
  ! ... read file
  READ(iunplot, '(A)') title
  READ(iunplot, *)     nsite, nr1x, nr2x, nr3x, nr1, nr2, nr3, nat, ntyp
  READ(iunplot, *)     ibrav, celldm
  IF (ibrav == 0) THEN
    DO i = 1, 3
      READ(iunplot, *) (at(ipol, i), ipol = 1, 3)
    END DO
  END IF
  READ(iunplot, *) gcutm, laue
  !
  ! ... close file
  CLOSE(unit=iunplot)
  !
END SUBROUTINE read_rism_header
