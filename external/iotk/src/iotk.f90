
program iotk
  use iotk_module
  use iotk_base
  use iotk_error_interf
  implicit none
  integer :: nargs,ierrl
  character(iotk_linlenx) :: args(iotk_maxargs)

  call iotk_readcmdline(args,nargs,eos=.true.,ierr=ierrl)
  if(ierrl/=0) goto 1
  call iotk_tool(args(1:nargs))
1 continue
  if(ierrl/=0) call iotk_error_handler(ierrl)
end program iotk
