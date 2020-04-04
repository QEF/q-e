program simple_bse

  USE  start_end
  USE input_simple_exc
  
  implicit none

  TYPE(input_options) :: sin


!setup MPI environment
  call startup

  call read_input_simple_exc( sin )  
  select case(sin%task)
     case(0)!solve eigen-problem
         call simple_eigen(sin)
     case(1)!find spectrum lanczos
        call lanczos(sin)
  end select

  call stop_run
  stop

end program simple_bse
