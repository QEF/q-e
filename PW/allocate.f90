!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
module allocate
  use error_handler
  use parameters, only : DP
  integer, public             :: maxsize, totsize

  integer, private, parameter :: real_len     = 8,&
       integer_len  = 4,&
       complex_len  = 16

  interface mallocate
     module procedure allocate_real_1d
     module procedure allocate_real_2d
     module procedure allocate_real_3d
     module procedure allocate_real_4d
     module procedure allocate_real_5d
     module procedure allocate_real_6d
     module procedure allocate_real_7d

     module procedure allocate_complex_1d
     module procedure allocate_complex_2d
     module procedure allocate_complex_3d
     module procedure allocate_complex_4d
     module procedure allocate_complex_5d
     module procedure allocate_complex_6d
     module procedure allocate_complex_7d

     module procedure allocate_integer_1d
     module procedure allocate_integer_2d
     module procedure allocate_integer_3d
     module procedure allocate_integer_4d
     module procedure allocate_integer_5d
     module procedure allocate_integer_6d
     module procedure allocate_integer_7d
  end interface

  interface mfree
     module procedure free_real_1d
     module procedure free_real_2d
     module procedure free_real_3d
     module procedure free_real_4d
     module procedure free_real_5d
     module procedure free_real_6d
     module procedure free_real_7d

     module procedure free_complex_1d
     module procedure free_complex_2d
     module procedure free_complex_3d
     module procedure free_complex_4d
     module procedure free_complex_5d
     module procedure free_complex_6d
     module procedure free_complex_7d

     module procedure free_integer_1d
     module procedure free_integer_2d
     module procedure free_integer_3d
     module procedure free_integer_4d
     module procedure free_integer_5d
     module procedure free_integer_6d
     module procedure free_integer_7d
  end interface

contains

  subroutine allocate_real_1d(var,dim1) 
    real(kind = DP), pointer, dimension(:) :: var
    integer                     :: dim1 
    integer                     :: ierr
    call add_name("allocate_real_1d")
    if (associated(var)) call error_mem("Variable already allocated")
    allocate(var(dim1),STAT=ierr)
    if (ierr .ne. 0) call error_mem("Memory allocation failed") 
    totsize = totsize + real_len*dim1
    maxsize = max (maxsize, totsize)
    call chop_name
    return
  end subroutine allocate_real_1d

  subroutine free_real_1d(vin) 
    real(kind = DP), pointer, dimension(:) :: vin
    integer                     :: ierr
    call add_name("free_real_1d")
    totsize = totsize - real_len*product(shape(vin))
    deallocate(vin,STAT=ierr)
    if (ierr .ne. 0) call error_mem("Memory deallocation failed") 
    call chop_name
    return
  end subroutine free_real_1d

  subroutine allocate_real_2d(var,dim1,dim2) 
    real(kind = DP), pointer, dimension(:,:) :: var
    integer                       :: dim1,dim2
    integer                     :: ierr
    call add_name("allocate_real_2d")
    if (associated(var)) call error_mem("Variable already allocated")
    allocate(var(dim1,dim2),STAT=ierr)
    if (ierr .ne. 0) call error_mem("Memory allocation failed") 
    totsize = totsize + real_len*dim1*dim2
    maxsize = max (maxsize, totsize)
    call chop_name
    return
  end subroutine allocate_real_2d

  subroutine free_real_2d(vin) 
    real(kind = DP), pointer, dimension(:,:) :: vin
    integer                     :: ierr
    call add_name("free_real_2d")
    totsize = totsize - real_len*product(shape(vin))
    deallocate(vin,STAT=ierr)
    if (ierr .ne. 0) call error_mem("Memory deallocation failed") 
    call chop_name
    return
  end subroutine free_real_2d
  
  subroutine allocate_real_3d(var,dim1,dim2,dim3) 
    real(kind = DP), pointer, dimension(:,:,:) :: var
    integer                         :: dim1,dim2,dim3
    integer                     :: ierr
    call add_name("allocate_real_3d")
    if (associated(var)) call error_mem("Variable already allocated")
    allocate(var(dim1,dim2,dim3),STAT=ierr)
    if (ierr .ne. 0) call error_mem("Memory allocation failed") 
    totsize = totsize + real_len*dim1*dim2*dim3 
    maxsize = max (maxsize, totsize)
    call chop_name
    return
  end subroutine allocate_real_3d

  subroutine free_real_3d(vin) 
    real(kind = DP), pointer, dimension(:,:,:) :: vin
    integer                     :: ierr
    call add_name("free_real_3d")
    totsize = totsize - real_len*product(shape(vin))
    deallocate(vin,STAT=ierr)
    if (ierr .ne. 0) call error_mem("Memory deallocation failed") 
    call chop_name
    return
  end subroutine free_real_3d

  subroutine allocate_real_4d(var,dim1,dim2,dim3,dim4) 
    real(kind = DP), pointer, dimension(:,:,:,:) :: var
    integer                           :: dim1,dim2,dim3,dim4
    integer                     :: ierr
    call add_name("allocate_real_4d")
    if (associated(var)) call error_mem("Variable already allocated")
    allocate(var(dim1,dim2,dim3,dim4),STAT=ierr)
    if (ierr .ne. 0) call error_mem("Memory allocation failed") 
    totsize = totsize + real_len*dim1*dim2*dim3*dim4
    maxsize = max (maxsize, totsize)
    call chop_name
    return
  end subroutine allocate_real_4d

  subroutine free_real_4d(vin) 
    real(kind = DP), pointer, dimension(:,:,:,:) :: vin
    integer                     :: ierr
    call add_name("free_real_4d")
    totsize = totsize - real_len*product(shape(vin))
    deallocate(vin,STAT=ierr)
    if (ierr .ne. 0) call error_mem("Memory deallocation failed") 
    call chop_name
    return
  end subroutine free_real_4d

  subroutine allocate_real_5d(var,dim1,dim2,dim3,dim4,dim5) 
    real(kind = DP), pointer, dimension(:,:,:,:,:) :: var
    integer                           :: dim1,dim2,dim3,dim4,dim5
    integer                     :: ierr
    call add_name("allocate_real_5d")
    if (associated(var)) call error_mem("Variable already allocated")
    allocate(var(dim1,dim2,dim3,dim4,dim5),STAT=ierr)
    if (ierr .ne. 0) call error_mem("Memory allocation failed") 
    totsize = totsize + real_len*dim1*dim2*dim3*dim4*dim5
    maxsize = max (maxsize, totsize)
    call chop_name
    return
  end subroutine allocate_real_5d

  subroutine free_real_5d(vin) 
    real(kind = DP), pointer, dimension(:,:,:,:,:) :: vin
    integer                     :: ierr
    call add_name("free_real_5d")
    totsize = totsize - real_len*product(shape(vin))
    deallocate(vin,STAT=ierr)
    if (ierr .ne. 0) call error_mem("Memory deallocation failed") 
    call chop_name
    return
  end subroutine free_real_5d


  subroutine allocate_real_6d(var,dim1,dim2,dim3,dim4,dim5,dim6) 
    real(kind = DP), pointer, dimension(:,:,:,:,:,:) :: var
    integer                           :: dim1,dim2,dim3,dim4,dim5,dim6
    integer                     :: ierr
    call add_name("allocate_real_6d")
    if (associated(var)) call error_mem("Variable already allocated")
    allocate(var(dim1,dim2,dim3,dim4,dim5,dim6),STAT=ierr)
    if (ierr .ne. 0) call error_mem("Memory allocation failed") 
    totsize = totsize + real_len*dim1*dim2*dim3*dim4*dim5*dim6
    maxsize = max (maxsize, totsize)
    call chop_name
    return
  end subroutine allocate_real_6d

  subroutine free_real_6d(vin) 
    real(kind = DP), pointer, dimension(:,:,:,:,:,:) :: vin
    integer                     :: ierr
    call add_name("free_real_6d")
    totsize = totsize - real_len*product(shape(vin))
    deallocate(vin,STAT=ierr)
    if (ierr .ne. 0) call error_mem("Memory deallocation failed") 
    call chop_name
    return
  end subroutine free_real_6d


  subroutine allocate_real_7d(var,dim1,dim2,dim3,dim4,dim5,dim6,dim7) 
    real(kind = DP), pointer, dimension(:,:,:,:,:,:,:) :: var
    integer                           :: dim1,dim2,dim3,dim4,dim5,dim6,dim7
    integer                     :: ierr
    call add_name("allocate_real_7d")
    if (associated(var)) call error_mem("Variable already allocated")
    allocate(var(dim1,dim2,dim3,dim4,dim5,dim6,dim7),STAT=ierr)
    if (ierr .ne. 0) call error_mem("Memory allocation failed") 
    totsize = totsize + real_len*dim1*dim2*dim3*dim4*dim5*dim6*dim7
    maxsize = max (maxsize, totsize)
    call chop_name
    return
  end subroutine allocate_real_7d

  subroutine free_real_7d(vin) 
    real(kind = DP), pointer, dimension(:,:,:,:,:,:,:) :: vin
    integer                     :: ierr
    call add_name("free_real_7d")
    totsize = totsize - real_len*product(shape(vin))
    deallocate(vin,STAT=ierr)
    if (ierr .ne. 0) call error_mem("Memory allocation failed") 
    call chop_name
    return
  end subroutine free_real_7d
!  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!  %%%%%%%%%%%%%%%%%%%%%%%%   END REAL   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  subroutine allocate_complex_1d(var,dim1) 
    complex(kind = DP), pointer, dimension(:) :: var
    integer                     :: dim1 
    integer                     :: ierr
    call add_name("allocate_complex_1d")
    if (associated(var)) call error_mem("Variable already allocated")
    allocate(var(dim1),STAT=ierr)
    if (ierr .ne. 0) call error_mem("Memory allocation failed") 
    totsize = totsize + complex_len*dim1
    maxsize = max (maxsize, totsize)
    call chop_name
    return
  end subroutine allocate_complex_1d

  subroutine free_complex_1d(vin) 
    complex(kind=DP), pointer, dimension(:) :: vin
    integer                     :: ierr
    call add_name("free_complex_1d")
    totsize = totsize - complex_len*product(shape(vin))
    deallocate(vin,STAT=ierr)
    if (ierr .ne. 0) call error_mem("Memory deallocation failed") 
    call chop_name
    return
  end subroutine free_complex_1d

  subroutine allocate_complex_2d(var,dim1,dim2) 
    complex(kind=DP), pointer, dimension(:,:) :: var
    integer                       :: dim1,dim2
    integer                     :: ierr
    call add_name("allocate_complex_2d")
    if (associated(var)) call error_mem("Variable already allocated")
    allocate(var(dim1,dim2),STAT=ierr)
    if (ierr .ne. 0) call error_mem("Memory allocation failed") 
    totsize = totsize + complex_len*dim1*dim2
    maxsize = max (maxsize, totsize)
    call chop_name
    return
  end subroutine allocate_complex_2d

  subroutine free_complex_2d(vin) 
    complex(kind=DP), pointer, dimension(:,:) :: vin
    integer                     :: ierr
    call add_name("free_complex_2d")
    totsize = totsize - complex_len*product(shape(vin))
    deallocate(vin,STAT=ierr)
    if (ierr .ne. 0) call error_mem("Memory deallocation failed") 
    call chop_name
    return
  end subroutine free_complex_2d

  subroutine allocate_complex_3d(var,dim1,dim2,dim3) 
    complex(kind=DP), pointer, dimension(:,:,:) :: var
    integer                         :: dim1,dim2,dim3
    integer                     :: ierr
    call add_name("allocate_complex_3d")
    if (associated(var)) call error_mem("Variable already allocated")
    allocate(var(dim1,dim2,dim3),STAT=ierr)
    if (ierr .ne. 0) call error_mem("Memory allocation failed") 
    totsize = totsize + complex_len*dim1*dim2*dim3 
    maxsize = max (maxsize, totsize)
    call chop_name
    return
  end subroutine allocate_complex_3d

  subroutine free_complex_3d(vin) 
    complex(kind=DP), pointer, dimension(:,:,:) :: vin
    integer                     :: ierr
    call add_name("free_complex_3d")
    totsize = totsize - complex_len*product(shape(vin))
    deallocate(vin,STAT=ierr)
    if (ierr .ne. 0) call error_mem("Memory deallocation failed") 
    call chop_name
    return
  end subroutine free_complex_3d

  subroutine allocate_complex_4d(var,dim1,dim2,dim3,dim4) 
    complex(kind=DP), pointer, dimension(:,:,:,:) :: var
    integer                           :: dim1,dim2,dim3,dim4
    integer                     :: ierr
    call add_name("allocate_complex_4d")
    if (associated(var)) call error_mem("Variable already allocated")
    allocate(var(dim1,dim2,dim3,dim4),STAT=ierr)
    if (ierr .ne. 0) call error_mem("Memory allocation failed") 
    totsize = totsize + complex_len*dim1*dim2*dim3*dim4
    maxsize = max (maxsize, totsize)
    call chop_name
    return
  end subroutine allocate_complex_4d

  subroutine free_complex_4d(vin) 
    complex(kind=DP), pointer, dimension(:,:,:,:) :: vin
    integer                     :: ierr
    call add_name("free_complex_4d")
    totsize = totsize - complex_len*product(shape(vin))
    deallocate(vin,STAT=ierr)
    if (ierr .ne. 0) call error_mem("Memory deallocation failed") 
    call chop_name
    return
  end subroutine free_complex_4d

  subroutine allocate_complex_5d(var,dim1,dim2,dim3,dim4,dim5) 
    complex(kind=DP), pointer, dimension(:,:,:,:,:) :: var
    integer                           :: dim1,dim2,dim3,dim4,dim5
    integer                     :: ierr
    call add_name("allocate_complex_5d")
    if (associated(var)) call error_mem("Variable already allocated")
    allocate(var(dim1,dim2,dim3,dim4,dim5),STAT=ierr)
    if (ierr .ne. 0) call error_mem("Memory allocation failed") 
    totsize = totsize + complex_len*dim1*dim2*dim3*dim4*dim5
    maxsize = max (maxsize, totsize)
    call chop_name
    return
  end subroutine allocate_complex_5d

  subroutine free_complex_5d(vin) 
    complex(kind=DP), pointer, dimension(:,:,:,:,:) :: vin
    integer                     :: ierr
    call add_name("free_complex_5d")
    totsize = totsize - complex_len*product(shape(vin))
    deallocate(vin,STAT=ierr)
    if (ierr .ne. 0) call error_mem("Memory deallocation failed") 
    call chop_name
    return
  end subroutine free_complex_5d

  subroutine allocate_complex_6d(var,dim1,dim2,dim3,dim4,dim5,dim6) 
    complex(kind = DP), pointer, dimension(:,:,:,:,:,:) :: var
    integer                           :: dim1,dim2,dim3,dim4,dim5,dim6
    integer                     :: ierr
    call add_name("allocate_complex_6d")
    if (associated(var)) call error_mem("Variable already allocated")
    allocate(var(dim1,dim2,dim3,dim4,dim5,dim6),STAT=ierr)
    if (ierr .ne. 0) call error_mem("Memory allocation failed") 
    totsize = totsize + complex_len*dim1*dim2*dim3*dim4*dim5*dim6
    maxsize = max (maxsize, totsize)
    call chop_name
    return
  end subroutine allocate_complex_6d

  subroutine free_complex_6d(vin) 
    complex(kind = DP), pointer, dimension(:,:,:,:,:,:) :: vin
    integer                     :: ierr
    call add_name("free_complex_6d")
    totsize = totsize - complex_len*product(shape(vin))
    deallocate(vin,STAT=ierr)
    if (ierr .ne. 0) call error_mem("Memory deallocation failed") 
    call chop_name
    return
  end subroutine free_complex_6d


  subroutine allocate_complex_7d(var,dim1,dim2,dim3,dim4,dim5,dim6,dim7) 
    complex(kind = DP), pointer, dimension(:,:,:,:,:,:,:) :: var
    integer                           :: dim1,dim2,dim3,dim4,dim5,dim6,dim7
    integer                     :: ierr
    call add_name("allocate_complex_7d")
    if (associated(var)) call error_mem("Variable already allocated")
    allocate(var(dim1,dim2,dim3,dim4,dim5,dim6,dim7),STAT=ierr)
    if (ierr .ne. 0) call error_mem("Memory allocation failed") 
    totsize = totsize + complex_len*dim1*dim2*dim3*dim4*dim5*dim6*dim7
    maxsize = max (maxsize, totsize)
    call chop_name
    return
  end subroutine allocate_complex_7d

  subroutine free_complex_7d(vin) 
    complex(kind = DP), pointer, dimension(:,:,:,:,:,:,:) :: vin
    integer                     :: ierr
    call add_name("free_complex_7d")
    totsize = totsize - complex_len*product(shape(vin))
    deallocate(vin,STAT=ierr)
    if (ierr .ne. 0) call error_mem("Memory deallocation failed") 
    call chop_name
    return
  end subroutine free_complex_7d
!  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!  %%%%%%%%%%%%%%%%%%%%%%%%   END COMPLEX  %%%%%%%%%%%%%%%%%%%%%%%%%%%%
!  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  subroutine allocate_integer_1d(var,dim1) 
    integer, pointer, dimension(:) :: var
    integer                     :: dim1 
    integer                     :: ierr
    call add_name("allocate_integer_1d")
    if (associated(var)) call error_mem("Variable already allocated")
    allocate(var(dim1),STAT=ierr)
    if (ierr .ne. 0) call error_mem("Memory allocation failed") 
    totsize = totsize + integer_len*dim1
    maxsize = max (maxsize, totsize)
    call chop_name
    return
  end subroutine allocate_integer_1d

  subroutine free_integer_1d(vin) 
    integer, pointer, dimension(:) :: vin
    integer                     :: ierr
    call add_name("free_integer_1d")
    totsize = totsize - integer_len*product(shape(vin))
    deallocate(vin,STAT=ierr)
    if (ierr .ne. 0) call error_mem("Memory deallocation failed") 
    call chop_name
    return
  end subroutine free_integer_1d

  subroutine allocate_integer_2d(var,dim1,dim2) 
    integer, pointer, dimension(:,:) :: var
    integer                       :: dim1,dim2
    integer                     :: ierr
    call add_name("allocate_integer_2d")
    if (associated(var)) call error_mem("Variable already allocated")
    allocate(var(dim1,dim2),STAT=ierr)
    if (ierr .ne. 0) call error_mem("Memory allocation failed") 
    totsize = totsize + integer_len*dim1*dim2
    maxsize = max (maxsize, totsize)
    call chop_name
    return
  end subroutine allocate_integer_2d

  subroutine free_integer_2d(vin) 
    integer, pointer, dimension(:,:) :: vin
    integer                     :: ierr
    call add_name("free_integer_2d")
    totsize = totsize - integer_len*product(shape(vin))
    deallocate(vin,STAT=ierr)
    if (ierr .ne. 0) call error_mem("Memory deallocation failed") 
    call chop_name
    return
  end subroutine free_integer_2d

  subroutine allocate_integer_3d(var,dim1,dim2,dim3) 
    integer, pointer, dimension(:,:,:) :: var
    integer                         :: dim1,dim2,dim3
    integer                     :: ierr
    call add_name("allocate_integer_3d")
    if (associated(var)) call error_mem("Variable already allocated")
    allocate(var(dim1,dim2,dim3),STAT=ierr)
    if (ierr .ne. 0) call error_mem("Memory allocation failed") 
    totsize = totsize + integer_len*dim1*dim2*dim3 
    maxsize = max (maxsize, totsize)
    call chop_name
    return
  end subroutine allocate_integer_3d

  subroutine free_integer_3d(vin) 
    integer, pointer, dimension(:,:,:) :: vin
    integer                     :: ierr
    call add_name("free_integer_3d")
    totsize = totsize - integer_len*product(shape(vin))
    deallocate(vin,STAT=ierr)
    if (ierr .ne. 0) call error_mem("Memory deallocation failed") 
    call chop_name
    return
  end subroutine free_integer_3d

  subroutine allocate_integer_4d(var,dim1,dim2,dim3,dim4) 
    integer, pointer, dimension(:,:,:,:) :: var
    integer                           :: dim1,dim2,dim3,dim4
    integer                     :: ierr
    call add_name("allocate_integer_4d")
    if (associated(var)) call error_mem("Variable already allocated")
    allocate(var(dim1,dim2,dim3,dim4),STAT=ierr)
    if (ierr .ne. 0) call error_mem("Memory allocation failed") 
    totsize = totsize + integer_len*dim1*dim2*dim3*dim4
    maxsize = max (maxsize, totsize)
    call chop_name
    return
  end subroutine allocate_integer_4d

  subroutine free_integer_4d(vin) 
    integer, pointer, dimension(:,:,:,:) :: vin
    integer                     :: ierr
    call add_name("free_integer_4d")
    totsize = totsize - integer_len*product(shape(vin))
    deallocate(vin,STAT=ierr)
    if (ierr .ne. 0) call error_mem("Memory deallocation failed") 
    call chop_name
    return
  end subroutine free_integer_4d

  subroutine allocate_integer_5d(var,dim1,dim2,dim3,dim4,dim5) 
    integer, pointer, dimension(:,:,:,:,:) :: var
    integer                           :: dim1,dim2,dim3,dim4,dim5
    integer                     :: ierr
    call add_name("allocate_integer_5d")
    if (associated(var)) call error_mem("Variable already allocated")
    allocate(var(dim1,dim2,dim3,dim4,dim5),STAT=ierr)
    if (ierr .ne. 0) call error_mem("Memory allocation failed") 
    totsize = totsize + integer_len*dim1*dim2*dim3*dim4*dim5
    maxsize = max (maxsize, totsize)
    call chop_name
    return
  end subroutine allocate_integer_5d

  subroutine free_integer_5d(vin) 
    integer, pointer, dimension(:,:,:,:,:) :: vin
    integer                     :: ierr
    call add_name("free_integer_5d")
    totsize = totsize - integer_len*product(shape(vin))
    deallocate(vin,STAT=ierr)
    if (ierr .ne. 0) call error_mem("Memory deallocation failed") 
    call chop_name
    return
  end subroutine free_integer_5d

  subroutine allocate_integer_6d(var,dim1,dim2,dim3,dim4,dim5,dim6) 
    integer, pointer, dimension(:,:,:,:,:,:) :: var
    integer                           :: dim1,dim2,dim3,dim4,dim5,dim6
    integer                     :: ierr
    call add_name("allocate_integer_6d")
    if (associated(var)) call error_mem("Variable already allocated")
    allocate(var(dim1,dim2,dim3,dim4,dim5,dim6),STAT=ierr)
    if (ierr .ne. 0) call error_mem("Memory allocation failed") 
    totsize = totsize + integer_len*dim1*dim2*dim3*dim4*dim5*dim6
    maxsize = max (maxsize, totsize)
    call chop_name
    return
  end subroutine allocate_integer_6d

  subroutine free_integer_6d(vin) 
    integer, pointer, dimension(:,:,:,:,:,:) :: vin
    integer                     :: ierr
    call add_name("free_integer_6d")
    totsize = totsize - integer_len*product(shape(vin))
    deallocate(vin,STAT=ierr)
    if (ierr .ne. 0) call error_mem("Memory deallocation failed") 
    call chop_name
    return
  end subroutine free_integer_6d


  subroutine allocate_integer_7d(var,dim1,dim2,dim3,dim4,dim5,dim6,dim7) 
    integer, pointer, dimension(:,:,:,:,:,:,:) :: var
    integer                           :: dim1,dim2,dim3,dim4,dim5,dim6,dim7
    integer                     :: ierr
    call add_name("allocate_integer_7d")
    if (associated(var)) call error_mem("Variable already allocated")
    allocate(var(dim1,dim2,dim3,dim4,dim5,dim6,dim7),STAT=ierr)
    if (ierr .ne. 0) call error_mem("Memory allocation failed") 
    totsize = totsize + integer_len*dim1*dim2*dim3*dim4*dim5*dim6*dim7
    maxsize = max (maxsize, totsize)
    call chop_name
    return
  end subroutine allocate_integer_7d

  subroutine free_integer_7d(vin) 
    integer, pointer, dimension(:,:,:,:,:,:,:) :: vin
    integer                     :: ierr
    call add_name("free_integer_7d")
    totsize = totsize - integer_len*product(shape(vin))
    deallocate(vin,STAT=ierr)
    if (ierr .ne. 0) call error_mem("Memory deallocation failed") 
    call chop_name
    return
  end subroutine free_integer_7d
!  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!  %%%%%%%%%%%%%%%%%%%%%%%%   END INTEGER  % %%%%%%%%%%%%%%%%%%%%%%%%%%
!  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end module allocate
