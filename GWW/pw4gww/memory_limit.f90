! FOR GWW
!
! This module checks the amount of free memory
! based on the value of the memory per core
! available (provided by the user)
! Such value is defined in mem_per_core
!
! Author: G. Stenuit
!
!----------------------------------------------------------------------------
MODULE memory_limit
  !----------------------------------------------------------------------------
  !
! #ifdef __GWW
  !
  USE io_global,            ONLY : stdout, ionode
  USE kinds,    ONLY : i4b, DP
  USE mp, ONLY : mp_sum, mp_barrier
  USE wannier_gw,  ONLY : mem_per_core
  implicit none
  SAVE
  !
  REAL(kind=DP) :: mem_used, mem_free ! used and free memory in Byte !!

CONTAINS

  subroutine test_allocation_mp(n_row, n_col, choice, mem_u)
        ! This subroutine determines the size of a matrice
        ! that can be allocated on a NODE. In part., choice: 1=row; 2=col
        ! fixes which dimension (1=row,2=col) will be modified
        ! in order to test if the allocation is allowed
        implicit none
        integer, intent(inout) :: n_row, n_col ! the dimension of the matrix to be tested for alloc
        integer, intent(in) :: choice ! choice indicates which dim will be reduced if the allocation fails
        real(kind=DP), intent(in) :: mem_u !mem_u provides the amount (in kB) of already used memory
        integer :: err_alloc ! flag for the allocation (err_alloc == 0 => OK; != 0 => not OK)
        real(kind=DP) :: mem_esti ! estimation of the memory needed to allocation the matrx
                            ! this test is done since the err_alloc flags behaves strangely on nodes
        integer, parameter :: step_div=2
        complex(kind=DP), dimension(:,:), allocatable :: MAT_TEST

        write(stdout,*) 'Routine test_allocation_mp'
        call flush_unit(stdout)

        if ( (choice /= 1) .and. (choice /= 2) ) then
                write(stdout,*) 'Please choose 1 (for row) or 2 (for col) for the third argument !'
                call flush_unit(stdout)
                stop
        endif


        mem_esti=n_row*n_col
        mem_esti=mem_esti*16 ! 16 corresponds to the worst case: allocation of complex
        mem_free=mem_per_core-mem_u
        write(stdout,*) 'mem_limit=', mem_per_core/(1024*1024), 'MB (provides by the user)'
        write(stdout,*) 'mem_esti=', mem_esti/(1024*1024), ' MB'
        write(stdout,*) 'mem free=', mem_free/(1024*1024), ' MB'
        call flush_unit(stdout)
        if ( mem_esti >= mem_free ) then
                write(stdout,*) '!! Warning n_row or n_col too big ! !!'
                call flush_unit(stdout)
                if ( choice == 1 ) then
                        n_row=mem_free/(n_col*16)
                        n_row=n_row
                else
                        n_col=mem_free/(n_row*16)
                        n_col=n_col
                endif
                write(stdout,*) '=> program will use : n_row=', n_row, ' and n_col=', n_col
                call flush_unit(stdout)
        endif


        do
                write(stdout,*) 'Allocation : TESTS'
                call flush_unit(stdout)
                allocate(MAT_TEST(n_row, n_col), stat=err_alloc)

                if ( err_alloc == 0 ) then
                        write(stdout,*) 'Allocation of MAT(', n_row, ',', n_col, ') : OK'
                        call flush_unit(stdout)
                        deallocate(MAT_TEST)
                        exit
                else
                        write(stdout,*) 'Allocation of MAT(', n_row, ',', n_col, ') : FAILED !'
                        if (choice == 1) then
                                n_row=n_row/step_div
                        else
                                n_col=n_col/step_div
                        endif
                        write(stdout,*) 'Try now with n_row=', n_row, 'and n_col=', n_col
                        call flush_unit(stdout)
                endif
                if (n_row == 0 .or. n_col == 0) then
                        exit
                endif
        enddo

  end subroutine test_allocation_mp

!
! #endif __GWW
!
END MODULE memory_limit

