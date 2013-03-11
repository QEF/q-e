!
! Copyright (C) 2011 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------
MODULE dfile_autoname
!----------------------------------------------------------------------
  USE kinds,     ONLY : DP
  !
  PUBLIC :: dfile_name, dfile_generate_name, dfile_get_qlist
  !
  PRIVATE
  CHARACTER(len=12),PARAMETER :: dfile_directory_basename='.dfile_dir'
  !
CONTAINS
!----------------------------------------------------------------------
FUNCTION dfile_directory_file(basename, prefix)
  !----------------------------------------------------------------------
  IMPLICIT NONE
  CHARACTER(len=*),INTENT(in) :: basename
  CHARACTER(len=*),INTENT(in) :: prefix
  CHARACTER(len=512) :: dfile_directory_file
  dfile_directory_file = TRIM(prefix)//"."// &
                           TRIM(basename)//dfile_directory_basename
  RETURN
  !----------------------------------------------------------------------
END FUNCTION dfile_directory_file
!----------------------------------------------------------------------
!
!----------------------------------------------------------------------
FUNCTION open_dfile_directory(basename, prefix)
  !----------------------------------------------------------------------
  IMPLICIT NONE
  CHARACTER(len=*),INTENT(in) :: basename
  CHARACTER(len=*),INTENT(in) :: prefix   ! directory where to operate
  INTEGER :: open_dfile_directory, ios
  INTEGER, EXTERNAL :: find_free_unit
  CHARACTER(len=256) :: filename
  !
  filename = dfile_directory_file(basename, prefix)
  open_dfile_directory = find_free_unit()
  !
  OPEN(UNIT  = open_dfile_directory, &
       ACCESS= 'sequential',           &
       POSITION='append',              &
       FILE  = TRIM(filename),         &
       FORM  ='formatted', status='unknown', iostat=ios)
  !
  IF(ios/=0) CALL errore('open_dfile_directory','Cannot open: '//TRIM(filename),ABS(ios))
  !
  RETURN
  !----------------------------------------------------------------------
END FUNCTION open_dfile_directory
!----------------------------------------------------------------------
!
!----------------------------------------------------------------------
FUNCTION scan_dfile_directory(iunit, xq, at, found, equiv) 
  !----------------------------------------------------------------------
  IMPLICIT NONE
  CHARACTER(len=256) :: scan_dfile_directory
  !
  REAL(DP),INTENT(in) :: xq(3)
  REAL(DP),INTENT(in) :: at(3,3)
  INTEGER,INTENT(in)  :: iunit
  LOGICAL,INTENT(out) :: found
  LOGICAL,INTENT(in),OPTIONAL :: equiv ! if .false. only look for exactly q
                                       ! if .true. any q+G is ok (default)
  !
  INTEGER  :: ios, iq_
  REAL(DP) :: xp(3), aq(3), ap(3)
  CHARACTER(len=256) :: xp_name
  REAL(DP),PARAMETER :: gam(3) = (/ 0._dp, 0._dp, 0._dp /), accept = 1.e-5_dp
  !
  LOGICAL :: equiv_
  LOGICAL,EXTERNAL :: eqvect
  !
  found=.false.
  scan_dfile_directory = ''
  !
  equiv_ = .true.
  IF(present(equiv)) equiv_ = equiv
  ! xq in crystal coordinates:
  aq = xq
  CALL cryst_to_cart (1,aq,at,-1)
  !
  REWIND(iunit)
  ios=0
  !
  SCAN_FILE : &
  DO WHILE(ios==0)
      READ(iunit,*,iostat=ios) xp, ap, iq_, xp_name
    ! ap = xp
    ! CALL cryst_to_cart (1,ap,at,-1)
    !
    IF (equiv_) THEN
      IF (eqvect(aq,ap,gam,accept) .and. ios==0) THEN
        found=.true.
        scan_dfile_directory = TRIM(ADJUSTL(xp_name))
        EXIT SCAN_FILE
      ENDIF
    ELSE
      IF ( ALL(ABS(ap-aq)<accept) ) THEN
        found=.true.
        scan_dfile_directory = TRIM(ADJUSTL(xp_name))
        EXIT SCAN_FILE
      ENDIF
    ENDIF
  ENDDO SCAN_FILE
  !
  RETURN
  !----------------------------------------------------------------------
END FUNCTION scan_dfile_directory
!----------------------------------------------------------------------
!
!----------------------------------------------------------------------
FUNCTION dfile_name(xq, at, name, prefix, generate, index_q, equiv)
  !----------------------------------------------------------------------
  ! automatically generate a name for fildrho file
  USE io_global,    ONLY : ionode
  IMPLICIT NONE
  ! function:
  CHARACTER(len=256) :: dfile_name
  ! input variables:
  INTEGER, INTENT(in)         :: index_q      ! index of the q-point  
  REAL(DP),INTENT(in)         :: xq(3)    ! the q point in cartesian axes
  REAL(DP),INTENT(in)         :: at(3,3)  ! the lattice vectors, to transform the q to crystal coords
  CHARACTER(len=*),INTENT(in) :: prefix   ! directory where to operate
  CHARACTER(len=*),INTENT(in) :: name     ! input fildrho
  LOGICAL,INTENT(in)          :: generate ! make a new name if not found
  LOGICAL,INTENT(in),OPTIONAL :: equiv    ! accept an equivalent point q+G instead of q
                                          ! (only used when generate=.false.)
  !
  INTEGER :: iunit = -1, ios
  LOGICAL :: found
  CHARACTER(len=256) :: basename
  REAL(DP) :: aq(3)
  !
  ! Only ionode returns a meaningful name, as only ionode should do the i/o
  IF(.not.ionode) THEN
    dfile_name = ' '
    RETURN
  ENDIF
  !
  IF(name(1:5) /= 'auto:') THEN
    dfile_name = name
    RETURN
  ENDIF
  !
  basename = TRIM(name(6:))
  !
  
  iunit = open_dfile_directory(basename, prefix)
  rewind(iunit)
 
  dfile_name = scan_dfile_directory(iunit, xq, at, found, equiv)

  CLOSE(iunit)
  
  !
  ! Return here if point was found
  !IF(found) print*, "xq found as ", TRIM(dfile_name)
  IF(found) RETURN
  !
  IF(.not.generate) THEN
    WRITE(*,'(7x,"Error: ",3f12.6)') xq
    WRITE(*,'(7x,"Error: ",a,2x,a)') TRIM(name), TRIM(prefix)
    CALL errore('dfile_name','Requested q vector not found @ '//TRIM(basename), 1)
  ENDIF
  !
  ! Make up a new name
  dfile_name = TRIM(dfile_generate_name(xq, at, basename))
  !
  ! Append the new name to the list
  iunit = open_dfile_directory(basename, prefix)
  aq = xq
  CALL cryst_to_cart (1,aq,at,-1)
  !
  WRITE(iunit,*,iostat=ios) xq, aq, index_q, TRIM(dfile_name)

  IF(ios/=0) CALL errore('dfile_name','Cannot write dfile_directory',1)
  CLOSE(iunit)
  !
  RETURN
  !----------------------------------------------------------------------
END FUNCTION dfile_name
!----------------------------------------------------------------------
!
!----------------------------------------------------------------------
SUBROUTINE dfile_get_qlist(xqs, nqs, name, prefix, index_q)
  !----------------------------------------------------------------------
  ! automatically generate a name for fildrho file
  USE io_global,    ONLY : ionode
  IMPLICIT NONE
  ! input variables:
  INTEGER,INTENT(in)          :: nqs        ! max number of points
  INTEGER                     :: index_q(nqs)    ! index of irreducible q
  REAL(DP),INTENT(out)        :: xqs(3,nqs) ! the q point in cartesian axes
  CHARACTER(len=*),INTENT(in) :: prefix     ! directory where to operate
  CHARACTER(len=*),INTENT(in) :: name       ! input fildrho
  !
  INTEGER :: iunit = -1, ios, iq
  REAL(DP) :: xq_dummy(3)
  CHARACTER(len=256) :: basename
  !
  ! Only ionode scans for the filename, and does NOT broadcast. The broadcast
  ! must be done outside, because here we do not know if all processors are
  ! calling this subroutine!
  !
  IF (.not. ionode) THEN
   xqs = 0._dp
   RETURN
  ENDIF
  !
  IF(name(1:5) == 'auto:') THEN
    basename = TRIM(name(6:))
  ELSE
    basename = TRIM(name)
  ENDIF
  !
  iunit = open_dfile_directory(basename, prefix)
  REWIND(iunit)
  !
  GET_Q_LOOP : &
  DO iq = 1,nqs
     READ(iunit,*,iostat=ios) xqs(:,iq), xq_dummy(:),index_q(iq)
     IF(ios/=0) THEN
       CALL errore('dfile_get_qlist', 'Error while reading q point', iq)
     ENDIF
  ENDDO &
  GET_Q_LOOP
  !
  CLOSE(iunit)
  !
  RETURN
  !----------------------------------------------------------------------
END SUBROUTINE dfile_get_qlist
!----------------------------------------------------------------------
!
!----------------------------------------------------------------------
FUNCTION dfile_generate_name(xq, at, name)
  !----------------------------------------------------------------------
  ! automatically generate a name for fildrho file
  IMPLICIT NONE
  ! function:
  CHARACTER(len=256) :: dfile_generate_name
  ! input variables:
  REAL(DP),INTENT(in)         :: xq(3)   ! the q point in cartesian axes
  REAL(DP),INTENT(in)         :: at(3,3) ! lattice vectors
  CHARACTER(len=*),INTENT(in) :: name    ! input fildrho
  ! local variables:
  REAL(DP) :: aq(3) ! xq in crystal axes
  !
  ! take xq to crystalline coordinates
  aq = xq
  CALL cryst_to_cart (1,aq,at,-1)
  !
  WRITE(dfile_generate_name, '(a,".",a,"_",a,"_",a)') TRIM(name), &
        TRIM(real2frac(aq(1))), TRIM(real2frac(aq(2))), TRIM(real2frac(aq(3)))
  !
  dfile_generate_name = TRIM(dfile_generate_name)
  !
  RETURN
  !
  !----------------------------------------------------------------------
END FUNCTION dfile_generate_name
!----------------------------------------------------------------------
!
! Convert a real number to a string containing fraction, using 'o' as a division sign, es.
! 0    --> "0"
! -4   --> "-4"
! 0.25 --> "1o4"
! -1.66666666667 -> "-5/3"
! 
!----------------------------------------------------------------------
FUNCTION real2frac(r) RESULT (f)
  !----------------------------------------------------------------------
  IMPLICIT NONE
  REAL(DP),INTENT(in) :: r
  CHARACTER(len=64) :: f
  !
  INTEGER :: d, n
  INTEGER,PARAMETER :: max_denominator = 48000
  REAL(DP),PARAMETER :: accept = 1.d-6
  CHARACTER(len=64) :: nc,dc
  !
  IF(max_denominator*accept*20>1._dp) &
    CALL errore('real2frac', 'incompatible parameters', 2)
  ! Threat zero and integers separately:
  IF (ABS(r)<accept) THEN
    f = '0'
    RETURN
  ELSE IF (ABS(r-NINT(r))<accept) THEN
    WRITE(nc, '(i16)') NINT(r)
    f = TRIM(ADJUSTL(nc))
    RETURN
  ENDIF
  !
  ! If the number is not an integer, recompose the fraction
  DO d = 1, max_denominator+1
    IF( ABS(r*d-NINT(r*d)) < accept ) EXIT
  ENDDO
  !
!  IF (d > max_denominator) CALL errore('real2frac', 'not a fraction', 1)
  ! 
  IF (d > max_denominator) THEN
     WRITE(*, '("WARNING from real2frac:",e25.15," is not a fraction, falling back to hex." )') r
     WRITE(f,'(Z64)') r
     f='0x'//TRIM(ADJUSTL(f))
     RETURN
  ENDIF
  !
  n = NINT(r*d)
  !
  WRITE(nc, '(i16)') n
  WRITE(dc, '(i16)') d
  !
  f = TRIM(ADJUSTL(nc))//'o'//TRIM(ADJUSTL(dc))
  !
  RETURN
  !
  !----------------------------------------------------------------------
END FUNCTION real2frac
!----------------------------------------------------------------------
!
!----------------------------------------------------------------------
END MODULE dfile_autoname
!----------------------------------------------------------------------
