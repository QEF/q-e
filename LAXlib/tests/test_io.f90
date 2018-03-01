MODULE test_io
    !
    IMPLICIT NONE
    !
    INTERFACE read_problem
       MODULE PROCEDURE read_cmplx_problem, read_real_problem
    END INTERFACE
    !
    CONTAINS
    !
    SUBROUTINE read_cmplx_problem(fname, ldh, n, m, h, s, e, v)
        USE la_param, ONLY : DP
        IMPLICIT NONE
        character(len=*), intent(in) :: fname
        integer, intent(out) :: ldh, n, m
        complex(dp), allocatable, intent(inout) :: h(:,:)
        complex(dp), allocatable, intent(inout) :: s(:,:)
        real(dp), allocatable, intent(inout)    :: e(:)
        complex(dp), allocatable, intent(inout) :: v(:,:)
        real(dp) :: aux1, aux2
        !
        integer :: i, j, t
        logical :: exist
        !
        character(len=20):: fname_
    
        !
        fname_ = trim(fname)
        !
        print *, "reading ", fname_
        inquire(file=fname_, exist=exist)
        
        if (.not. exist ) then
            print *, "read_problem, could not read file"
            stop 1
        end if
        
        open (unit=15, file=fname_, status='old', form='unformatted', action='read')
        
        t = 0
        read (15)  n
        read (15)  m
        read (15)  ldh
        print *, n, m, ldh
        t = t + 3
        !
        ALLOCATE(h(ldh, n), s(ldh, n), e(n), v(ldh, m))
        !
        READ(15) h
        READ(15) s
        READ(15) e
        READ(15) v
        !
        close(15)
    END SUBROUTINE read_cmplx_problem
    !
    SUBROUTINE read_real_problem(fname, ldh, n, m, h, s, e, v)
        USE la_param, ONLY : DP
        IMPLICIT NONE
        character(len=*), intent(in) :: fname
        integer, intent(out) :: ldh, n, m
        real(dp), allocatable, intent(inout) :: h(:,:)
        real(dp), allocatable, intent(inout) :: s(:,:)
        real(dp), allocatable, intent(inout)    :: e(:)
        real(dp), allocatable, intent(inout) :: v(:,:)
        real(dp) :: aux1, aux2
        !
        integer :: i, j, t
        logical :: exist
        !
        character(len=20):: fname_
    
        !
        fname_ = trim(fname)
        !
        print *, "reading ", fname_
        inquire(file=fname_, exist=exist)
        
        if (.not. exist ) then
            print *, "read_problem, could not read file"
            stop 1
        end if
        
        open (unit=15, file=fname_, status='old', form='unformatted', action='read')
        
        t = 0
        read (15)  n
        read (15)  m
        read (15)  ldh
        print *, n, m, ldh
        t = t + 3
        !
        ALLOCATE(h(ldh, n), s(ldh, n), e(n), v(ldh, m))
        !
        READ(15) h
        READ(15) s
        READ(15) e
        READ(15) v
        !
        close(15)
    END SUBROUTINE read_real_problem
END MODULE test_io
