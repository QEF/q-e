!
! Copyright (C) 2001-2013 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!


MODULE para_gww
!this modules contains arrays indicating if the 
!processor should perform the task

SAVE

  LOGICAL, ALLOCATABLE :: is_my_time(:) !for 2n+1 times and frequencies
  LOGICAL  :: is_my_last!for extra 0 time calculation
  LOGICAL, ALLOCATABLE :: is_my_pola(:) !for 0 to n calculations
  LOGICAL, ALLOCATABLE :: is_my_state(:)!for KS states considered 1 to n_max
  LOGICAL, ALLOCATABLE :: is_my_state_range(:)!for KS states considered i_min to i_max

CONTAINS

  subroutine free_memory_para_gww
    
    implicit none

    if(allocated(is_my_time)) deallocate(is_my_time)
    if(allocated(is_my_pola)) deallocate(is_my_pola)
    if(allocated(is_my_state)) deallocate(is_my_state)
    if(allocated(is_my_state_range)) deallocate(is_my_state_range)


    return

  end subroutine free_memory_para_gww
  
  SUBROUTINE setup_para_gww(ntimes,nstates, i_min, i_max)

!this subroutine initialize the para variables for the gww
!calculation, parallelization is achieved on imaginary times
!and frequencies

    USE mp_world,  ONLY : mpime, nproc
    USE io_global, ONLY : stdout

    implicit none

    INTEGER, INTENT(in) :: ntimes!number of time samples
    INTEGER, INTENT(in) :: nstates!max number of states
    INTEGER, INTENT(in) :: i_min!lowest state for which the self-energy is calculated
    INTEGER, INTENT(in) :: i_max!upper  state for which the self-energy is calculated

    INTEGER :: ndelta, it, ip, iqq

!allocates arrays
    allocate(is_my_time(-ntimes:ntimes))
    allocate(is_my_pola(0:ntimes))
    allocate(is_my_state(nstates))
    allocate(is_my_state_range(i_min:i_max))
    
    is_my_time(:)=.false.
    is_my_pola(:)=.false.
    is_my_state(:)=.false.
    is_my_state_range(:)=.false.
    

    ndelta=(2*ntimes+1)/nproc
    if(ndelta*nproc < (2*ntimes+1)) ndelta=ndelta+1

    iqq=-ntimes
    do ip=0,nproc-1
       do it=1,ndelta
          if(iqq <= ntimes.and.(mpime==ip)) then
             is_my_time(iqq)=.true.
          endif
          if(it==1.and.(mpime==ip))  write(stdout,*) 'min', iqq, ip,it,ndelta
          if(it==ndelta.and.(mpime==ip))  write(stdout,*) 'max', iqq, ip,it,ndelta
          iqq=iqq+1
       enddo
    enddo


    if((mpime+1)==nproc) then
       is_my_last=.true.
    else
       is_my_last=.false.
    endif

    ndelta=(ntimes+1)/nproc 
    if(ndelta*nproc < (ntimes+1)) ndelta=ndelta+1

    iqq=0
    do ip=0,nproc-1
       do it=1,ndelta
          if(iqq <= ntimes.and.(mpime==ip)) then
             is_my_pola(iqq)=.true.
          endif
          if(it==1.and.(mpime==ip))  write(stdout,*) 'min pola', iqq
          if(it==ndelta.and.(mpime==ip))  write(stdout,*) 'max pola', iqq
          iqq=iqq+1
       enddo
    enddo


    ndelta=(nstates)/nproc
    if(ndelta*nproc < nstates) ndelta=ndelta+1

    iqq=1
    do ip=0,nproc-1
       do it=1,ndelta
          if(iqq <= nstates.and.(mpime==ip)) then
             is_my_state(iqq)=.true.
          endif
          if(it==1.and.(mpime==ip))  write(stdout,*) 'min state', iqq
          if(it==ndelta.and.(mpime==ip))  write(stdout,*) 'max state', iqq
          iqq=iqq+1
       enddo
    enddo

    ndelta=(i_max-i_min+1)/nproc
    if(ndelta*nproc <  (i_max-i_min+1)) ndelta=ndelta+1

    iqq=1
    do ip=0,nproc-1
       do it=1,ndelta
          if(iqq <= (i_max-i_min+1).and.(mpime==ip)) then
             is_my_state_range(iqq+i_min-1)=.true.
          endif
          if(it==1.and.(mpime==ip))  write(stdout,*) 'min state range', iqq +i_min-1
          if(it==ndelta.and.(mpime==ip))  write(stdout,*) 'max state range', iqq+i_min-1
          iqq=iqq+1
       enddo
    enddo

    


    return

  END SUBROUTINE setup_para_gww

END MODULE para_gww



