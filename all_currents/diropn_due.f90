subroutine diropn_due (prefix_due, unit, extension, recl, exst, tmp_dir_)
use io_files, only : prefix, diropn
implicit none
  character(len=*), intent(in) :: extension,prefix_due
  ! input: name of the file to open
  character(len=*), optional :: tmp_dir_
  ! optional variable, if present it is used as tmp_dir
  integer :: unit, recl
  ! input: unit of the file to open
  ! input: length of the records
  logical :: exst

character(len=256) :: old_prefix

old_prefix=prefix
prefix=prefix_due
call diropn(unit,extension,recl,exst,tmp_dir_)
prefix=old_prefix

end subroutine


