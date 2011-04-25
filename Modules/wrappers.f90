!
! Copyright (C) 2004-2009 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------------------------------------------
!
MODULE wrappers
  !
  ! these routines are used to pass fortran strings to C routines  in a
  ! safe way. Strings are converted to integer arrays here,  passed to 
  ! C wrappers, converted back to strings. Other ways to pass fortran
  ! strings to C turned out to be non portable and not safe
  !
   USE kinds, ONLY : DP
   USE io_global,  ONLY : stdout
   IMPLICIT NONE
   SAVE
CONTAINS
   !
   FUNCTION feval_infix( ierr, str )
     REAL(DP) :: feval_infix
     INTEGER :: ierr
     CHARACTER(LEN=*) :: str
     INTEGER :: i, ilen
     INTEGER, ALLOCATABLE :: istr(:)
     REAL(DP), EXTERNAL :: eval_infix_wrapper
     ALLOCATE( istr( LEN( str ) ) )
     DO i = 1, LEN( str )
        istr(i) = ICHAR( str(i:i) )
        IF( istr(i) < 0 .OR. istr(i) > 127 ) &
           CALL errore( ' feval_infix ', ' invalid character ', ABS( istr(i) ) )
     END DO 
     ilen = LEN( str )
     feval_infix = eval_infix_wrapper( ierr, istr, ilen )
     DEALLOCATE( istr )
     RETURN
   END FUNCTION
   !
   FUNCTION f_mkdir( dirname )
     INTEGER :: f_mkdir
     CHARACTER(LEN=*) :: dirname
     INTEGER :: i, ilen
     INTEGER, ALLOCATABLE :: istr(:)
     INTEGER, EXTERNAL :: c_mkdir_int
     ALLOCATE( istr( LEN_TRIM( dirname ) ) )
     DO i = 1, LEN_TRIM( dirname )
        istr(i) = ICHAR( dirname(i:i) )
        IF( istr(i) < 0 .OR. istr(i) > 127 ) &
           CALL errore( ' f_mkdir ', ' invalid character ', ABS( istr(i) ) )
     END DO 
     ilen = LEN_TRIM( dirname )
     f_mkdir = c_mkdir_int( istr, ilen )
     DEALLOCATE( istr )
     RETURN
   END FUNCTION
   !
   FUNCTION f_chdir( dirname )
     INTEGER :: f_chdir
     CHARACTER(LEN=*) :: dirname
     INTEGER :: i, ilen
     INTEGER, ALLOCATABLE :: istr(:)
     INTEGER, EXTERNAL :: c_chdir_int
     ALLOCATE( istr( LEN_TRIM( dirname ) ) )
     DO i = 1, LEN_TRIM( dirname )
        istr(i) = ICHAR( dirname(i:i) )
        IF( istr(i) < 0 .OR. istr(i) > 127 ) &
           CALL errore( ' f_chdir ', ' invalid character ', ABS( istr(i) ) )
     END DO
     ilen = LEN_TRIM( dirname )
     f_chdir = c_chdir_int( istr, ilen )
     DEALLOCATE( istr )
     RETURN
   END FUNCTION
   !
   FUNCTION f_rename( oldname, newname )
     INTEGER :: f_rename
     CHARACTER(LEN=*) :: oldname
     CHARACTER(LEN=*) :: newname
     INTEGER :: i, lold, lnew
     INTEGER, ALLOCATABLE :: iold(:)
     INTEGER, ALLOCATABLE :: inew(:)
     INTEGER, EXTERNAL :: c_rename_int
     lold = LEN( oldname )
     lnew = LEN( newname )
     ALLOCATE( iold( lold ) )
     ALLOCATE( inew( lnew ) )
     DO i = 1, lold
        iold(i) = ICHAR( oldname(i:i) )
        IF( iold(i) < 0 .OR. iold(i) > 127 ) &
           CALL errore( ' f_rename ', ' invalid character ', ABS( iold(i) ) )
     END DO
     DO i = 1, lnew
        inew(i) = ICHAR( newname(i:i) )
        IF( inew(i) < 0 .OR. inew(i) > 127 ) &
           CALL errore( ' f_rename ', ' invalid character ', ABS( inew(i) ) )
     END DO
     f_rename = c_rename_int( iold, lold, inew, lnew )
     DEALLOCATE( inew )
     DEALLOCATE( iold )
     RETURN
   END FUNCTION
   !
   FUNCTION f_link( oldname, newname )
     INTEGER :: f_link
     CHARACTER(LEN=*) :: oldname
     CHARACTER(LEN=*) :: newname
     INTEGER :: i, lold, lnew
     INTEGER, ALLOCATABLE :: iold(:)
     INTEGER, ALLOCATABLE :: inew(:)
     INTEGER, EXTERNAL :: c_link_int
     lold = LEN( oldname )
     lnew = LEN( newname )
     ALLOCATE( iold( lold ) )
     ALLOCATE( inew( lnew ) )
     DO i = 1, lold
        iold(i) = ICHAR( oldname(i:i) )
        IF( iold(i) < 0 .OR. iold(i) > 127 ) &
           CALL errore( ' f_link ', ' invalid character ', ABS( iold(i) ) )
     END DO
     DO i = 1, lnew
        inew(i) = ICHAR( newname(i:i) )
        IF( inew(i) < 0 .OR. inew(i) > 127 ) &
           CALL errore( ' f_link ', ' invalid character ', ABS( inew(i) ) )
     END DO
     f_link = c_link_int( iold, lold, inew, lnew )
     DEALLOCATE( inew )
     DEALLOCATE( iold )
     RETURN
   END FUNCTION
   !

   SUBROUTINE md5_from_file (filename, md5)
     
     CHARACTER(LEN=*), INTENT (IN) :: filename
     CHARACTER(len=32), INTENT (OUT) :: md5
     CHARACTER(LEN=len_trim(filename)) :: ftrim

     INTEGER , EXTERNAL :: file_md5
     INTEGER, ALLOCATABLE :: istr(:)
     INTEGER :: getter(32), retval, i, ilen

     ftrim = TRIM(filename)     
     ilen = LEN(ftrim )
     ALLOCATE( istr( ilen ) )
     DO i = 1, ilen
        istr(i) = ichar( ftrim(i:i) )
     ENDDO
     
     retval = file_md5( istr, ilen, getter )

     DO i = 1,32
        md5(i:i) = char(  getter(i) )
     ENDDO

     DEALLOCATE( istr )

  END SUBROUTINE 
  !

END MODULE
