program main
use iotk_module
use iotk_files_interf
implicit none

character(500) :: source,dest
character(iotk_namlenx) :: root
character(iotk_attlenx) :: attr
character :: key

logical        :: binary_dest,binary_source
integer :: maxsize
#ifdef __MPI
integer :: ierr
#endif

#ifdef __MPI
call MPI_init(ierr)
#endif

write(*,*) "Name of input file:"
read(*,*)  source
!write(*,*) "Textual or Binary:"
!read(*,*)  key
binary_source = .false.
!if(key=="B" .or. key=="b") binary_source = .true.
call iotk_magic(trim(source),binary=binary_source)
write(*,*) "Name of output file:"
read(*,*)  dest
write(*,*) "Textual or Binary:"
read(*,*)  key
binary_dest = .false.
maxsize = -1
if(key=="B" .or. key=="b") then
  binary_dest = .true.
else
  write(*,*) "Maximum data size (-1 for unlimited size)"
  read(*,*)  maxsize
end if

call iotk_open_read(60, trim(source),binary=binary_source,root=root,attr=attr)
call iotk_open_write(61,trim(dest),  binary=binary_dest,  root=root,attr=attr)
call iotk_copy_tag(60,61,maxsize=maxsize)
call iotk_close_write(61)
call iotk_close_read(60)

#ifdef __MPI
call MPI_Finalize(ierr)
#endif
end program main
