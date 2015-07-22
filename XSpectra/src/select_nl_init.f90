subroutine select_nl_init(edge, nl_init, two_edges, n_lanczos)
  USE kinds, ONLY : DP
  USE io_global, ONLY : stdout
  IMPLICIT NONE
  logical :: two_edges
  integer :: n_lanczos
  integer, dimension(2) :: nl_init
  character(LEN=16) :: edge
  CHARACTER(LEN=10) :: dummy_char, n_0, l_0
 
  if( TRIM(ADJUSTL(edge)) .eq. 'K' ) edge = 'K1'
  dummy_char = TRIM(ADJUSTL(edge))
  n_0 = dummy_char(1:1)
  if( len(TRIM(ADJUSTL(edge))) == 2 ) then
    l_0 = dummy_char(2:2)
    two_edges = .false.
  else
    l_0 = dummy_char(2:3)
    two_edges = .true.
  end if
 
  select case (n_0)
    case('K')
      nl_init(1) = 1
    case('L')
      nl_init(1) = 2
    case('M')
      nl_init(1) = 3
    case('N')
      nl_init(1) = 4
    case default
      write(stdout,*) 'Needs to be extended'
  end select

   select case (l_0)
    case('1')
      nl_init(2) = 0
      n_lanczos = 1
    case('2')
      nl_init(2) = 1
      n_lanczos = 2    ! the manifold
    case('3')
      nl_init(2) = 1
      n_lanczos = 4
    case('4')
      nl_init(2) = 2
      n_lanczos = 4
    case('5')
      nl_init(2) = 2
      n_lanczos = 6
    case('6')
      nl_init(2) = 3
      n_lanczos = 6
    case('7')
      nl_init(2) = 3
      n_lanczos = 6
    case('23')
      nl_init(2) = 1
      n_lanczos = 6
    case('45')
      nl_init(2) = 2
      n_lanczos = 10
    case('67')
      nl_init(2) = 3
      n_lanczos = 14
    case default
      write(stdout,*) 'Needs to be extended'
  end select


end subroutine select_nl_init
