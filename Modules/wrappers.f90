!
! Copyright (C) 2004-2013 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------------------------------------------
!
! This module contains fortran wrappers to POSIX system calls. 
! The wrappers are used to convert the Fortran CHARACTER array to
! null-terminated C *char. The conversion and the interface is done
! with the F95 intrinsic iso_c_binding module.
!
! Additionally, it provides interfaces to the C functions in clib/: 
! eval_infix, md5_from_file, f_mkdir_safe
!
! NOTE: the mkdir function is NOT called directly as it returns error if
!       directory already exists. We use instead a C wrapper c_mkdir_safe
!
MODULE wrappers
  USE kinds, ONLY : DP
  USE io_global, ONLY : stdout
  USE ISO_C_BINDING
  IMPLICIT NONE
  !
  ! C std library functions fortran wrappers:
  PUBLIC  f_remove, f_rename, f_chdir, f_mkdir, f_rmdir, f_getcwd
  ! more stuff:
  PUBLIC  f_copy, feval_infix, md5_from_file, f_mkdir_safe
  !
  ! HELP:
  ! integer f_remove(pathname)
  ! integer f_rename(oldfile, newfile)
  ! integer f_chdir(newdir)
  ! integer f_chmod(mode) i.e. mode=777 (disable)
  ! integer f_mkdir(dirname, mode) mode is optional
  ! integer f_rmdir(dirname)
  ! subroutine f_getcwd(dirname) 
  ! All "*name" are fortran characters of implicit length,
  ! "mode" are integers, all functions return 0 if successful, -1 otherwise
  !
  ! real(dp) :: result = feval_infix(integer:: ierr, character(len=*) :: expression)
  ! subroutine md5_from_file(character(len=*) :: filename, character(len=32) ::md5)
  PRIVATE
  !
  SAVE
  !
  ! Interfaces to the C functions, these are kept private as Fortran
  ! characters have (?) to be converted explicitly to C character arrays.
  ! Use the f_* wrappers instead
  INTERFACE
    FUNCTION remove(pathname) BIND(C,name="remove") RESULT(r)
      USE iso_c_binding
      CHARACTER(kind=c_char),INTENT(in)  :: pathname(*)
      INTEGER(c_int)        :: r
    END FUNCTION
    FUNCTION rename(input,output) BIND(C,name="rename") RESULT(r)
      USE iso_c_binding
      CHARACTER(kind=c_char),INTENT(in) :: input(*)
      CHARACTER(kind=c_char),INTENT(in) :: output(*)
      INTEGER(c_int)        :: r
    END FUNCTION
    FUNCTION chmod(filename,mode) BIND(C,name="chmod") RESULT(r)
      USE iso_c_binding
      CHARACTER(kind=c_char),INTENT(in)  :: filename(*)
      INTEGER(c_int),VALUE  ,INTENT(in)  :: mode
      INTEGER(c_int)        :: r
    END FUNCTION
    FUNCTION chdir(filename) BIND(C,name="chdir") RESULT(r)
      USE iso_c_binding
      CHARACTER(kind=c_char),INTENT(in)  :: filename(*)
      INTEGER(c_int)        :: r
    END FUNCTION
    FUNCTION mkdir(dirname,mode) BIND(C,name="mkdir") RESULT(r)
      USE iso_c_binding
      CHARACTER(kind=c_char),INTENT(in)  :: dirname(*)
      INTEGER(c_int),VALUE  ,INTENT(in)  :: mode
      INTEGER(c_int)        :: r
    END FUNCTION
    FUNCTION rmdir(dirname) BIND(C,name="rmdir") RESULT(r)
      USE iso_c_binding
      CHARACTER(kind=c_char),INTENT(in)  :: dirname(*)
      INTEGER(c_int)        :: r
    END FUNCTION
    FUNCTION getcwd(buffer,size) BIND(C,name="getcwd") RESULT(r)
      USE iso_c_binding
      CHARACTER(kind=c_char) ,INTENT(out) :: buffer(*)
      INTEGER(c_size_t),VALUE,INTENT(in)  :: size
      TYPE(c_ptr)  :: r
    END FUNCTION
  END INTERFACE
  !
  ! ====================================================================
CONTAINS
  ! ====================================================================
  ! fortran wrappers functions that call the C functions after converting
  ! fortran characters to C character arrays
  FUNCTION f_remove(filename) RESULT(r)
    CHARACTER(*),INTENT(in)  :: filename
    INTEGER(c_int) :: r
    r= remove(TRIM(filename)//C_NULL_CHAR)
  END FUNCTION

  FUNCTION f_rename(input,output) RESULT(k)
    CHARACTER(*),INTENT(in)  :: input,output
    INTEGER :: k
    k= rename(TRIM(input)//C_NULL_CHAR,TRIM(output)//C_NULL_CHAR)
  END FUNCTION

  FUNCTION f_chdir(dirname) RESULT(r)
    CHARACTER(*),INTENT(in)  :: dirname
    INTEGER(c_int) :: r
    r= chdir(TRIM(dirname)//C_NULL_CHAR)
  END FUNCTION
  !
  ! f_mkdir, causes an ERROR if dirname already exists: use f_mkdir_safe instead
  FUNCTION f_mkdir(dirname, mode) RESULT(r)
    CHARACTER(*),INTENT(in)  :: dirname
    INTEGER,INTENT(in) :: mode
    INTEGER(c_int) :: r
    INTEGER(c_int) :: c_mode
    c_mode = INT(mode, kind=c_int)
    r= mkdir(TRIM(dirname)//C_NULL_CHAR, c_mode)
  END FUNCTION
  ! Note: permissions are usually in octal, e.g.:
  !       mode = o'640' => rw-r-----
  FUNCTION f_chmod(filename, mode) RESULT(r)
    CHARACTER(*),INTENT(in)  :: filename
    INTEGER,INTENT(in) :: mode
    INTEGER(c_int) :: r
    INTEGER(c_int) :: c_mode
    c_mode = INT(mode, kind=c_int)
    r= chmod(TRIM(filename)//C_NULL_CHAR, c_mode)
  END FUNCTION

  FUNCTION f_rmdir(dirname) RESULT(r)
    CHARACTER(*),INTENT(in)  :: dirname
    INTEGER(c_int) :: r
    r= rmdir(TRIM(dirname)//C_NULL_CHAR)
  END FUNCTION
  
  SUBROUTINE f_getcwd(output)
    CHARACTER(kind=c_char,len=*),INTENT(out) :: output
    TYPE(c_ptr) :: buffer
    INTEGER(C_SIZE_T) :: length,i  ! was kind=C_LONG, which fails on WIN32
    length=LEN(output)
    buffer=getcwd(output,length)
    DO i=1,length
      IF(output(i:i) == C_NULL_CHAR) EXIT
    ENDDO
    output(i:)=' '
  END SUBROUTINE
  ! ==================================================================== 
  ! copy a file, uses clibs/copy.c which currently does a binary copy
  ! using an 8kb buffer
  ! 
  ! returns:
  !  0 : no error
  ! -1 : cannot open source
  ! -2 : cannot open dest
  ! -3 : error while writing
  ! -4 : disk full while writing
  FUNCTION f_copy(source, dest) RESULT(r)
    INTERFACE
    FUNCTION c_copy(source,dest) BIND(C,name="copy") RESULT(r)
      USE iso_c_binding
      CHARACTER(kind=c_char),INTENT(in)  :: source(*), dest(*)
      INTEGER(c_int)         :: r
    END FUNCTION c_copy
    END INTERFACE
    CHARACTER(*),INTENT(in)  :: source, dest
    INTEGER(c_int) :: r
    r= c_copy(TRIM(source)//C_NULL_CHAR, TRIM(dest)//C_NULL_CHAR)
  END FUNCTION
  !
  ! safe mkdir from clib/c_mkdir.c that creates a directory, if necessary, 
  ! and checks permissions. It can be called in parallel.
  ! Returns: 0 = all ok
  !          1 = error
  !         -1 = the directory already existed and is properly writable
  FUNCTION f_mkdir_safe(dirname) RESULT(r)
    INTERFACE
    FUNCTION mkdir_safe(dirname) BIND(C,name="c_mkdir_safe") RESULT(r)
      USE iso_c_binding
      CHARACTER(kind=c_char),INTENT(in)  :: dirname(*)
      INTEGER(c_int)         :: r
    END FUNCTION mkdir_safe
    END INTERFACE
    CHARACTER(*),INTENT(in)  :: dirname
    INTEGER(c_int) :: r
    r= mkdir_safe(TRIM(dirname)//C_NULL_CHAR)
  END FUNCTION
  !
  ! Two more wrappers for eval_infix (simple algebric expression parser)
  ! and for get_md5 which computes the md5 sum of a file.
  !
  FUNCTION feval_infix(fierr, fstr)
    USE ISO_C_BINDING
    IMPLICIT NONE
    REAL(DP) :: feval_infix
    INTEGER :: fierr
    CHARACTER(len=*) :: fstr
    INTEGER :: filen
    !
    INTERFACE
    FUNCTION ceval_infix(cierr, cstr, cilen) BIND(C, name="eval_infix")
    !REAL(kind=c_double) FUNCTION ceval_infix(cierr, cstr, cilen) BIND(C, name="eval_infix")
    !  double eval_infix( int *ierr, const char *strExpression, int len )
      USE ISO_C_BINDING
      REAL(kind=c_double) :: ceval_infix
      INTEGER(kind=c_int)    :: cierr
      CHARACTER(kind=c_char) :: cstr(*)
      INTEGER(kind=c_int),VALUE :: cilen
    END FUNCTION ceval_infix
    END INTERFACE
    !
    INTEGER(kind=c_int) :: cierr
    INTEGER(kind=c_int) :: cilen
    CHARACTER(len=len_trim(fstr)+1,kind=c_char) :: cstr
    !
    INTEGER :: i
    !
    filen = len_trim(fstr)
    cilen = INT(filen, kind=c_int)
    DO i = 1,filen
      cstr(i:i) = fstr(i:i)
    ENDDO
    cstr(filen+1:filen+1)=C_NULL_CHAR
    !
    feval_infix = REAL( ceval_infix(cierr, cstr, cilen), kind=DP)
    fierr = INT(cierr)
    RETURN
  END FUNCTION feval_infix
  !
  !
  SUBROUTINE md5_from_file (ffile, fmd5)
    IMPLICIT NONE
    CHARACTER(LEN=*), INTENT (IN) :: ffile
    CHARACTER(len=32), INTENT (OUT) :: fmd5
    !
    INTERFACE
    SUBROUTINE cget_md5(cfile, cmd5, cierr) BIND(C, name="get_md5")
    ! void get_md5(const char *file, char *md5, int err)
      USE ISO_C_BINDING
      CHARACTER(kind=c_char) :: cfile(*)
      CHARACTER(kind=c_char) :: cmd5(*)
      INTEGER(kind=c_int)    :: cierr
    END SUBROUTINE cget_md5
    END INTERFACE
    !
    INTEGER,PARAMETER :: md5_length = 32
    INTEGER :: i
    !
    CHARACTER(len=len_trim(ffile)+1,kind=c_char) :: cfile!(*)
    CHARACTER(len=(md5_length+1),kind=c_char)    :: cmd5!(*)
    INTEGER(kind=c_int)    :: cierr
    !
    cfile = TRIM(ffile)//C_NULL_CHAR
    !
    CALL cget_md5(cfile, cmd5, cierr)
    !
    DO i = 1,md5_length
       fmd5(i:i) = cmd5(i:i)
    ENDDO
    !
  END SUBROUTINE 
END MODULE
! ==================================================================== 
