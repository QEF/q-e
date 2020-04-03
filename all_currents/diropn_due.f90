subroutine diropn_due (prefix_due, unit, extension, recl, exst, tmp_dir_)
  !-----------------------------------------------------------------------
  !
  !     this routine opens a file named "prefix"."extension" in tmp_dir 
  !     for direct I/O access
  !     If appropriate, the node number is added to the file name
  !
#if defined(__SX6)
#  define DIRECT_IO_FACTOR 1
#else
#  define DIRECT_IO_FACTOR 8 
#endif
  !
  ! the  record length in direct-access I/O is given by the number of
  ! real*8 words times DIRECT_IO_FACTOR (may depend on the compiler)
  !
  use kinds
  use io_files, only : nd_nmbr, tmp_dir
  implicit none
  !
  !    first the input variables
  !
  CHARACTER(len=*), intent(in) :: prefix_due         ! prepended to file names
  character(len=*), intent(in) :: extension
  ! input: name of the file to open
  character(len=256) ::access_test,form_test
  character(len=*), optional :: tmp_dir_
  ! optional variable, if present it is used as tmp_dir
  integer :: unit, recl,recl_test
  ! input: unit of the file to open
  ! input: length of the records
  logical :: exst
  ! output: if true the file exists
  !
  !    local variables
  !
  character(len=256) :: tempfile, filename
  ! complete file name
  integer :: ios
  integer*8 :: unf_recl
  ! used to check I/O operations
  ! length of the record
  logical :: opnd
  ! if true the file is already opened
  !
  if (unit < 0) call errore ('diropn_due', 'wrong unit', 1)
  !
  !    we first check that the file is not already openend
  !
  ios = 0
  
  inquire (unit = unit, opened = opnd)
  if (opnd) call errore ('diropn_due', "can't open a connected unit", abs(unit))
  !
  !      then we check the filename extension
  !
  if (extension == ' ') call errore ('diropn_due','filename extension not given',2)
  filename = trim(prefix_due) // '.' // trim(extension)
  !print *, 'ciao ciao' , trim(extension) // nd_nmbr, trim(filename)
  if (present(tmp_dir_)) then
     tempfile = trim(tmp_dir_) // trim(filename) //nd_nmbr
  else
     tempfile = trim(tmp_dir) // trim(filename) //nd_nmbr
  endif
  !tempfile = trim(tmp_dir) // trim(filename) //trim(nd_nmbr)

  inquire (file = tempfile, exist = exst)
  !
  !      the unit for record length is unfortunately machine-dependent
  !
  unf_recl = DIRECT_IO_FACTOR * int(recl, kind=kind(unf_recl))
  if (unf_recl <= 0) call errore ('diropn_due', 'wrong record length', 3)

  open (unit, file = trim(adjustl(tempfile)), iostat = ios, form = 'unformatted', &
       status = 'unknown', access = 'direct', recl = unf_recl)
  if (ios /= 0) call errore ('diropn_due', 'error opening '//trim(adjustl(tempfile)), unit)
  return
end subroutine diropn_due
