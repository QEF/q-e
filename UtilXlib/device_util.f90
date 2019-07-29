#if defined(__CUDA)
#define __HAVE_DEVICE
#endif
!
module device_util_m
  implicit none
!
interface dev_memcpy
    !
    subroutine dp_dev_memcpy_r1d(array_out, array_in, &
                                            range1, lbound1 )
      implicit none
      !   
      integer, parameter :: PRCSN = selected_real_kind(14,200)
      real(PRCSN), intent(inout) :: array_out(:)
      real(PRCSN), intent(in)    :: array_in(:)
      integer, optional, intent(in) ::  range1(2)
      integer, optional, intent(in) ::  lbound1
#if defined(__CUDA)
      attributes(device) :: array_out, array_in
#endif
       !
    end subroutine dp_dev_memcpy_r1d
    !
    subroutine dp_dev_memcpy_r2d(array_out, array_in, &
                                            range1, lbound1, &
                                            range2, lbound2 )
      implicit none
      !   
      integer, parameter :: PRCSN = selected_real_kind(14,200)
      real(PRCSN), intent(inout) :: array_out(:,:)
      real(PRCSN), intent(in)    :: array_in(:,:)
      integer, optional, intent(in) ::  range1(2), range2(2)
      integer, optional, intent(in) ::  lbound1, lbound2
#if defined(__CUDA)
      attributes(device) :: array_out, array_in
#endif
       !
    end subroutine dp_dev_memcpy_r2d
    !
    subroutine dp_dev_memcpy_r3d(array_out, array_in, &
                                            range1, lbound1, &
                                            range2, lbound2, &
                                            range3, lbound3 )
      implicit none
      !   
      integer, parameter :: PRCSN = selected_real_kind(14,200)
      real(PRCSN), intent(inout) :: array_out(:,:,:)
      real(PRCSN), intent(in)    :: array_in(:,:,:)
      integer, optional, intent(in) ::  range1(2), range2(2), range3(2)
      integer, optional, intent(in) ::  lbound1, lbound2, lbound3
#if defined(__CUDA)
      attributes(device) :: array_out, array_in
#endif
       !
    end subroutine dp_dev_memcpy_r3d
    !
    subroutine dp_dev_memcpy_r4d(array_out, array_in, &
                                            range1, lbound1, &
                                            range2, lbound2, &
                                            range3, lbound3, &
                                            range4, lbound4 )
      implicit none
      !   
      integer, parameter :: PRCSN = selected_real_kind(14,200)
      real(PRCSN), intent(inout) :: array_out(:,:,:,:)
      real(PRCSN), intent(in)    :: array_in(:,:,:,:)
      integer, optional, intent(in) ::  range1(2), range2(2), range3(2), range4(2)
      integer, optional, intent(in) ::  lbound1, lbound2, lbound3, lbound4
#if defined(__CUDA)
      attributes(device) :: array_out, array_in
#endif
       !
    end subroutine dp_dev_memcpy_r4d
    !
    subroutine dp_dev_memcpy_c1d(array_out, array_in, &
                                            range1, lbound1 )
      implicit none
      !   
      integer, parameter :: PRCSN = selected_real_kind(14,200)
      complex(PRCSN), intent(inout) :: array_out(:)
      complex(PRCSN), intent(in)    :: array_in(:)
      integer, optional, intent(in) ::  range1(2)
      integer, optional, intent(in) ::  lbound1
#if defined(__CUDA)
      attributes(device) :: array_out, array_in
#endif
       !
    end subroutine dp_dev_memcpy_c1d
    !
    subroutine dp_dev_memcpy_c2d(array_out, array_in, &
                                            range1, lbound1, &
                                            range2, lbound2 )
      implicit none
      !   
      integer, parameter :: PRCSN = selected_real_kind(14,200)
      complex(PRCSN), intent(inout) :: array_out(:,:)
      complex(PRCSN), intent(in)    :: array_in(:,:)
      integer, optional, intent(in) ::  range1(2), range2(2)
      integer, optional, intent(in) ::  lbound1, lbound2
#if defined(__CUDA)
      attributes(device) :: array_out, array_in
#endif
       !
    end subroutine dp_dev_memcpy_c2d
    !
    subroutine dp_dev_memcpy_c3d(array_out, array_in, &
                                            range1, lbound1, &
                                            range2, lbound2, &
                                            range3, lbound3 )
      implicit none
      !   
      integer, parameter :: PRCSN = selected_real_kind(14,200)
      complex(PRCSN), intent(inout) :: array_out(:,:,:)
      complex(PRCSN), intent(in)    :: array_in(:,:,:)
      integer, optional, intent(in) ::  range1(2), range2(2), range3(2)
      integer, optional, intent(in) ::  lbound1, lbound2, lbound3
#if defined(__CUDA)
      attributes(device) :: array_out, array_in
#endif
       !
    end subroutine dp_dev_memcpy_c3d
    !
    subroutine dp_dev_memcpy_c4d(array_out, array_in, &
                                            range1, lbound1, &
                                            range2, lbound2, &
                                            range3, lbound3, &
                                            range4, lbound4 )
      implicit none
      !   
      integer, parameter :: PRCSN = selected_real_kind(14,200)
      complex(PRCSN), intent(inout) :: array_out(:,:,:,:)
      complex(PRCSN), intent(in)    :: array_in(:,:,:,:)
      integer, optional, intent(in) ::  range1(2), range2(2), range3(2), range4(2)
      integer, optional, intent(in) ::  lbound1, lbound2, lbound3, lbound4
#if defined(__CUDA)
      attributes(device) :: array_out, array_in
#endif
       !
    end subroutine dp_dev_memcpy_c4d
    !
    subroutine sp_dev_memcpy_r1d(array_out, array_in, &
                                            range1, lbound1 )
      implicit none
      !   
      integer, parameter :: PRCSN = selected_real_kind(6, 37)
      real(PRCSN), intent(inout) :: array_out(:)
      real(PRCSN), intent(in)    :: array_in(:)
      integer, optional, intent(in) ::  range1(2)
      integer, optional, intent(in) ::  lbound1
#if defined(__CUDA)
      attributes(device) :: array_out, array_in
#endif
       !
    end subroutine sp_dev_memcpy_r1d
    !
    subroutine sp_dev_memcpy_r2d(array_out, array_in, &
                                            range1, lbound1, &
                                            range2, lbound2 )
      implicit none
      !   
      integer, parameter :: PRCSN = selected_real_kind(6, 37)
      real(PRCSN), intent(inout) :: array_out(:,:)
      real(PRCSN), intent(in)    :: array_in(:,:)
      integer, optional, intent(in) ::  range1(2), range2(2)
      integer, optional, intent(in) ::  lbound1, lbound2
#if defined(__CUDA)
      attributes(device) :: array_out, array_in
#endif
       !
    end subroutine sp_dev_memcpy_r2d
    !
    subroutine sp_dev_memcpy_r3d(array_out, array_in, &
                                            range1, lbound1, &
                                            range2, lbound2, &
                                            range3, lbound3 )
      implicit none
      !   
      integer, parameter :: PRCSN = selected_real_kind(6, 37)
      real(PRCSN), intent(inout) :: array_out(:,:,:)
      real(PRCSN), intent(in)    :: array_in(:,:,:)
      integer, optional, intent(in) ::  range1(2), range2(2), range3(2)
      integer, optional, intent(in) ::  lbound1, lbound2, lbound3
#if defined(__CUDA)
      attributes(device) :: array_out, array_in
#endif
       !
    end subroutine sp_dev_memcpy_r3d
    !
    subroutine sp_dev_memcpy_r4d(array_out, array_in, &
                                            range1, lbound1, &
                                            range2, lbound2, &
                                            range3, lbound3, &
                                            range4, lbound4 )
      implicit none
      !   
      integer, parameter :: PRCSN = selected_real_kind(6, 37)
      real(PRCSN), intent(inout) :: array_out(:,:,:,:)
      real(PRCSN), intent(in)    :: array_in(:,:,:,:)
      integer, optional, intent(in) ::  range1(2), range2(2), range3(2), range4(2)
      integer, optional, intent(in) ::  lbound1, lbound2, lbound3, lbound4
#if defined(__CUDA)
      attributes(device) :: array_out, array_in
#endif
       !
    end subroutine sp_dev_memcpy_r4d
    !
    subroutine sp_dev_memcpy_c1d(array_out, array_in, &
                                            range1, lbound1 )
      implicit none
      !   
      integer, parameter :: PRCSN = selected_real_kind(6, 37)
      complex(PRCSN), intent(inout) :: array_out(:)
      complex(PRCSN), intent(in)    :: array_in(:)
      integer, optional, intent(in) ::  range1(2)
      integer, optional, intent(in) ::  lbound1
#if defined(__CUDA)
      attributes(device) :: array_out, array_in
#endif
       !
    end subroutine sp_dev_memcpy_c1d
    !
    subroutine sp_dev_memcpy_c2d(array_out, array_in, &
                                            range1, lbound1, &
                                            range2, lbound2 )
      implicit none
      !   
      integer, parameter :: PRCSN = selected_real_kind(6, 37)
      complex(PRCSN), intent(inout) :: array_out(:,:)
      complex(PRCSN), intent(in)    :: array_in(:,:)
      integer, optional, intent(in) ::  range1(2), range2(2)
      integer, optional, intent(in) ::  lbound1, lbound2
#if defined(__CUDA)
      attributes(device) :: array_out, array_in
#endif
       !
    end subroutine sp_dev_memcpy_c2d
    !
    subroutine sp_dev_memcpy_c3d(array_out, array_in, &
                                            range1, lbound1, &
                                            range2, lbound2, &
                                            range3, lbound3 )
      implicit none
      !   
      integer, parameter :: PRCSN = selected_real_kind(6, 37)
      complex(PRCSN), intent(inout) :: array_out(:,:,:)
      complex(PRCSN), intent(in)    :: array_in(:,:,:)
      integer, optional, intent(in) ::  range1(2), range2(2), range3(2)
      integer, optional, intent(in) ::  lbound1, lbound2, lbound3
#if defined(__CUDA)
      attributes(device) :: array_out, array_in
#endif
       !
    end subroutine sp_dev_memcpy_c3d
    !
    subroutine sp_dev_memcpy_c4d(array_out, array_in, &
                                            range1, lbound1, &
                                            range2, lbound2, &
                                            range3, lbound3, &
                                            range4, lbound4 )
      implicit none
      !   
      integer, parameter :: PRCSN = selected_real_kind(6, 37)
      complex(PRCSN), intent(inout) :: array_out(:,:,:,:)
      complex(PRCSN), intent(in)    :: array_in(:,:,:,:)
      integer, optional, intent(in) ::  range1(2), range2(2), range3(2), range4(2)
      integer, optional, intent(in) ::  lbound1, lbound2, lbound3, lbound4
#if defined(__CUDA)
      attributes(device) :: array_out, array_in
#endif
       !
    end subroutine sp_dev_memcpy_c4d
    !
end interface dev_memcpy
!
!
interface dev_memset
    !
    subroutine dp_dev_memset_r1d(array_out, val, &
                                            
                                            range1, lbound1 )
       implicit none
       !
       integer, parameter :: PRCSN = selected_real_kind(14,200)
       real(PRCSN), intent(inout) :: array_out(:)
       real(PRCSN), intent(in)    :: val
       integer, optional, intent(in) ::  range1(2)
       integer, optional, intent(in) ::  lbound1
#if defined(__CUDA)
       attributes(device) :: array_out
#endif
       !
    end subroutine dp_dev_memset_r1d
    !
    subroutine dp_dev_memset_r2d(array_out, val, &
                                            
                                            range1, lbound1, &
                                            range2, lbound2 )
       implicit none
       !
       integer, parameter :: PRCSN = selected_real_kind(14,200)
       real(PRCSN), intent(inout) :: array_out(:,:)
       real(PRCSN), intent(in)    :: val
       integer, optional, intent(in) ::  range1(2), range2(2)
       integer, optional, intent(in) ::  lbound1, lbound2
#if defined(__CUDA)
       attributes(device) :: array_out
#endif
       !
    end subroutine dp_dev_memset_r2d
    !
    subroutine dp_dev_memset_r3d(array_out, val, &
                                            
                                            range1, lbound1, &
                                            range2, lbound2, &
                                            range3, lbound3 )
       implicit none
       !
       integer, parameter :: PRCSN = selected_real_kind(14,200)
       real(PRCSN), intent(inout) :: array_out(:,:,:)
       real(PRCSN), intent(in)    :: val
       integer, optional, intent(in) ::  range1(2), range2(2), range3(2)
       integer, optional, intent(in) ::  lbound1, lbound2, lbound3
#if defined(__CUDA)
       attributes(device) :: array_out
#endif
       !
    end subroutine dp_dev_memset_r3d
    !
    subroutine dp_dev_memset_r4d(array_out, val, &
                                            
                                            range1, lbound1, &
                                            range2, lbound2, &
                                            range3, lbound3, &
                                            range4, lbound4 )
       implicit none
       !
       integer, parameter :: PRCSN = selected_real_kind(14,200)
       real(PRCSN), intent(inout) :: array_out(:,:,:,:)
       real(PRCSN), intent(in)    :: val
       integer, optional, intent(in) ::  range1(2), range2(2), range3(2), range4(2)
       integer, optional, intent(in) ::  lbound1, lbound2, lbound3, lbound4
#if defined(__CUDA)
       attributes(device) :: array_out
#endif
       !
    end subroutine dp_dev_memset_r4d
    !
    subroutine dp_dev_memset_c1d(array_out, val, &
                                            
                                            range1, lbound1 )
       implicit none
       !
       integer, parameter :: PRCSN = selected_real_kind(14,200)
       complex(PRCSN), intent(inout) :: array_out(:)
       complex(PRCSN), intent(in)    :: val
       integer, optional, intent(in) ::  range1(2)
       integer, optional, intent(in) ::  lbound1
#if defined(__CUDA)
       attributes(device) :: array_out
#endif
       !
    end subroutine dp_dev_memset_c1d
    !
    subroutine dp_dev_memset_c2d(array_out, val, &
                                            
                                            range1, lbound1, &
                                            range2, lbound2 )
       implicit none
       !
       integer, parameter :: PRCSN = selected_real_kind(14,200)
       complex(PRCSN), intent(inout) :: array_out(:,:)
       complex(PRCSN), intent(in)    :: val
       integer, optional, intent(in) ::  range1(2), range2(2)
       integer, optional, intent(in) ::  lbound1, lbound2
#if defined(__CUDA)
       attributes(device) :: array_out
#endif
       !
    end subroutine dp_dev_memset_c2d
    !
    subroutine dp_dev_memset_c3d(array_out, val, &
                                            
                                            range1, lbound1, &
                                            range2, lbound2, &
                                            range3, lbound3 )
       implicit none
       !
       integer, parameter :: PRCSN = selected_real_kind(14,200)
       complex(PRCSN), intent(inout) :: array_out(:,:,:)
       complex(PRCSN), intent(in)    :: val
       integer, optional, intent(in) ::  range1(2), range2(2), range3(2)
       integer, optional, intent(in) ::  lbound1, lbound2, lbound3
#if defined(__CUDA)
       attributes(device) :: array_out
#endif
       !
    end subroutine dp_dev_memset_c3d
    !
    subroutine dp_dev_memset_c4d(array_out, val, &
                                            
                                            range1, lbound1, &
                                            range2, lbound2, &
                                            range3, lbound3, &
                                            range4, lbound4 )
       implicit none
       !
       integer, parameter :: PRCSN = selected_real_kind(14,200)
       complex(PRCSN), intent(inout) :: array_out(:,:,:,:)
       complex(PRCSN), intent(in)    :: val
       integer, optional, intent(in) ::  range1(2), range2(2), range3(2), range4(2)
       integer, optional, intent(in) ::  lbound1, lbound2, lbound3, lbound4
#if defined(__CUDA)
       attributes(device) :: array_out
#endif
       !
    end subroutine dp_dev_memset_c4d
    !
    subroutine sp_dev_memset_r1d(array_out, val, &
                                            
                                            range1, lbound1 )
       implicit none
       !
       integer, parameter :: PRCSN = selected_real_kind(6, 37)
       real(PRCSN), intent(inout) :: array_out(:)
       real(PRCSN), intent(in)    :: val
       integer, optional, intent(in) ::  range1(2)
       integer, optional, intent(in) ::  lbound1
#if defined(__CUDA)
       attributes(device) :: array_out
#endif
       !
    end subroutine sp_dev_memset_r1d
    !
    subroutine sp_dev_memset_r2d(array_out, val, &
                                            
                                            range1, lbound1, &
                                            range2, lbound2 )
       implicit none
       !
       integer, parameter :: PRCSN = selected_real_kind(6, 37)
       real(PRCSN), intent(inout) :: array_out(:,:)
       real(PRCSN), intent(in)    :: val
       integer, optional, intent(in) ::  range1(2), range2(2)
       integer, optional, intent(in) ::  lbound1, lbound2
#if defined(__CUDA)
       attributes(device) :: array_out
#endif
       !
    end subroutine sp_dev_memset_r2d
    !
    subroutine sp_dev_memset_r3d(array_out, val, &
                                            
                                            range1, lbound1, &
                                            range2, lbound2, &
                                            range3, lbound3 )
       implicit none
       !
       integer, parameter :: PRCSN = selected_real_kind(6, 37)
       real(PRCSN), intent(inout) :: array_out(:,:,:)
       real(PRCSN), intent(in)    :: val
       integer, optional, intent(in) ::  range1(2), range2(2), range3(2)
       integer, optional, intent(in) ::  lbound1, lbound2, lbound3
#if defined(__CUDA)
       attributes(device) :: array_out
#endif
       !
    end subroutine sp_dev_memset_r3d
    !
    subroutine sp_dev_memset_r4d(array_out, val, &
                                            
                                            range1, lbound1, &
                                            range2, lbound2, &
                                            range3, lbound3, &
                                            range4, lbound4 )
       implicit none
       !
       integer, parameter :: PRCSN = selected_real_kind(6, 37)
       real(PRCSN), intent(inout) :: array_out(:,:,:,:)
       real(PRCSN), intent(in)    :: val
       integer, optional, intent(in) ::  range1(2), range2(2), range3(2), range4(2)
       integer, optional, intent(in) ::  lbound1, lbound2, lbound3, lbound4
#if defined(__CUDA)
       attributes(device) :: array_out
#endif
       !
    end subroutine sp_dev_memset_r4d
    !
    subroutine sp_dev_memset_c1d(array_out, val, &
                                            
                                            range1, lbound1 )
       implicit none
       !
       integer, parameter :: PRCSN = selected_real_kind(6, 37)
       complex(PRCSN), intent(inout) :: array_out(:)
       complex(PRCSN), intent(in)    :: val
       integer, optional, intent(in) ::  range1(2)
       integer, optional, intent(in) ::  lbound1
#if defined(__CUDA)
       attributes(device) :: array_out
#endif
       !
    end subroutine sp_dev_memset_c1d
    !
    subroutine sp_dev_memset_c2d(array_out, val, &
                                            
                                            range1, lbound1, &
                                            range2, lbound2 )
       implicit none
       !
       integer, parameter :: PRCSN = selected_real_kind(6, 37)
       complex(PRCSN), intent(inout) :: array_out(:,:)
       complex(PRCSN), intent(in)    :: val
       integer, optional, intent(in) ::  range1(2), range2(2)
       integer, optional, intent(in) ::  lbound1, lbound2
#if defined(__CUDA)
       attributes(device) :: array_out
#endif
       !
    end subroutine sp_dev_memset_c2d
    !
    subroutine sp_dev_memset_c3d(array_out, val, &
                                            
                                            range1, lbound1, &
                                            range2, lbound2, &
                                            range3, lbound3 )
       implicit none
       !
       integer, parameter :: PRCSN = selected_real_kind(6, 37)
       complex(PRCSN), intent(inout) :: array_out(:,:,:)
       complex(PRCSN), intent(in)    :: val
       integer, optional, intent(in) ::  range1(2), range2(2), range3(2)
       integer, optional, intent(in) ::  lbound1, lbound2, lbound3
#if defined(__CUDA)
       attributes(device) :: array_out
#endif
       !
    end subroutine sp_dev_memset_c3d
    !
    subroutine sp_dev_memset_c4d(array_out, val, &
                                            
                                            range1, lbound1, &
                                            range2, lbound2, &
                                            range3, lbound3, &
                                            range4, lbound4 )
       implicit none
       !
       integer, parameter :: PRCSN = selected_real_kind(6, 37)
       complex(PRCSN), intent(inout) :: array_out(:,:,:,:)
       complex(PRCSN), intent(in)    :: val
       integer, optional, intent(in) ::  range1(2), range2(2), range3(2), range4(2)
       integer, optional, intent(in) ::  lbound1, lbound2, lbound3, lbound4
#if defined(__CUDA)
       attributes(device) :: array_out
#endif
       !
    end subroutine sp_dev_memset_c4d
    !
    !
end interface dev_memset
!
!
interface dev_memsync
    !
    subroutine dp_h2d_memsync_r1d(array_out, array_in, &
                                             range1, lbound1  )
#if defined(__CUDA)
       use cudafor
#endif
       implicit none
       !   
       integer, parameter :: PRCSN = selected_real_kind(14,200)
       real(PRCSN), intent(inout) :: array_out(:) 
       real(PRCSN), intent(in)    :: array_in(:) 
       integer, optional, intent(in) ::  range1(2)
       integer, optional, intent(in) ::  lbound1
#if defined(__CUDA)
       attributes(device) :: array_out
#endif
       !
    end subroutine dp_h2d_memsync_r1d
    !
    subroutine dp_h2d_memsync_r2d(array_out, array_in, &
                                             range1, lbound1 , &
                                             range2, lbound2  )
#if defined(__CUDA)
       use cudafor
#endif
       implicit none
       !   
       integer, parameter :: PRCSN = selected_real_kind(14,200)
       real(PRCSN), intent(inout) :: array_out(:,:) 
       real(PRCSN), intent(in)    :: array_in(:,:) 
       integer, optional, intent(in) ::  range1(2), range2(2)
       integer, optional, intent(in) ::  lbound1, lbound2
#if defined(__CUDA)
       attributes(device) :: array_out
#endif
       !
    end subroutine dp_h2d_memsync_r2d
    !
    subroutine dp_h2d_memsync_r3d(array_out, array_in, &
                                             range1, lbound1 , &
                                             range2, lbound2 , &
                                             range3, lbound3  )
#if defined(__CUDA)
       use cudafor
#endif
       implicit none
       !   
       integer, parameter :: PRCSN = selected_real_kind(14,200)
       real(PRCSN), intent(inout) :: array_out(:,:,:) 
       real(PRCSN), intent(in)    :: array_in(:,:,:) 
       integer, optional, intent(in) ::  range1(2), range2(2), range3(2)
       integer, optional, intent(in) ::  lbound1, lbound2, lbound3
#if defined(__CUDA)
       attributes(device) :: array_out
#endif
       !
    end subroutine dp_h2d_memsync_r3d
    !
    subroutine dp_h2d_memsync_r4d(array_out, array_in, &
                                             range1, lbound1 , &
                                             range2, lbound2 , &
                                             range3, lbound3 , &
                                             range4, lbound4  )
#if defined(__CUDA)
       use cudafor
#endif
       implicit none
       !   
       integer, parameter :: PRCSN = selected_real_kind(14,200)
       real(PRCSN), intent(inout) :: array_out(:,:,:,:) 
       real(PRCSN), intent(in)    :: array_in(:,:,:,:) 
       integer, optional, intent(in) ::  range1(2), range2(2), range3(2), range4(2)
       integer, optional, intent(in) ::  lbound1, lbound2, lbound3, lbound4
#if defined(__CUDA)
       attributes(device) :: array_out
#endif
       !
    end subroutine dp_h2d_memsync_r4d
    !
    subroutine dp_h2d_memsync_c1d(array_out, array_in, &
                                             range1, lbound1  )
#if defined(__CUDA)
       use cudafor
#endif
       implicit none
       !   
       integer, parameter :: PRCSN = selected_real_kind(14,200)
       complex(PRCSN), intent(inout) :: array_out(:) 
       complex(PRCSN), intent(in)    :: array_in(:) 
       integer, optional, intent(in) ::  range1(2)
       integer, optional, intent(in) ::  lbound1
#if defined(__CUDA)
       attributes(device) :: array_out
#endif
       !
    end subroutine dp_h2d_memsync_c1d
    !
    subroutine dp_h2d_memsync_c2d(array_out, array_in, &
                                             range1, lbound1 , &
                                             range2, lbound2  )
#if defined(__CUDA)
       use cudafor
#endif
       implicit none
       !   
       integer, parameter :: PRCSN = selected_real_kind(14,200)
       complex(PRCSN), intent(inout) :: array_out(:,:) 
       complex(PRCSN), intent(in)    :: array_in(:,:) 
       integer, optional, intent(in) ::  range1(2), range2(2)
       integer, optional, intent(in) ::  lbound1, lbound2
#if defined(__CUDA)
       attributes(device) :: array_out
#endif
       !
    end subroutine dp_h2d_memsync_c2d
    !
    subroutine dp_h2d_memsync_c3d(array_out, array_in, &
                                             range1, lbound1 , &
                                             range2, lbound2 , &
                                             range3, lbound3  )
#if defined(__CUDA)
       use cudafor
#endif
       implicit none
       !   
       integer, parameter :: PRCSN = selected_real_kind(14,200)
       complex(PRCSN), intent(inout) :: array_out(:,:,:) 
       complex(PRCSN), intent(in)    :: array_in(:,:,:) 
       integer, optional, intent(in) ::  range1(2), range2(2), range3(2)
       integer, optional, intent(in) ::  lbound1, lbound2, lbound3
#if defined(__CUDA)
       attributes(device) :: array_out
#endif
       !
    end subroutine dp_h2d_memsync_c3d
    !
    subroutine dp_h2d_memsync_c4d(array_out, array_in, &
                                             range1, lbound1 , &
                                             range2, lbound2 , &
                                             range3, lbound3 , &
                                             range4, lbound4  )
#if defined(__CUDA)
       use cudafor
#endif
       implicit none
       !   
       integer, parameter :: PRCSN = selected_real_kind(14,200)
       complex(PRCSN), intent(inout) :: array_out(:,:,:,:) 
       complex(PRCSN), intent(in)    :: array_in(:,:,:,:) 
       integer, optional, intent(in) ::  range1(2), range2(2), range3(2), range4(2)
       integer, optional, intent(in) ::  lbound1, lbound2, lbound3, lbound4
#if defined(__CUDA)
       attributes(device) :: array_out
#endif
       !
    end subroutine dp_h2d_memsync_c4d
    !
    subroutine sp_h2d_memsync_r1d(array_out, array_in, &
                                             range1, lbound1  )
#if defined(__CUDA)
       use cudafor
#endif
       implicit none
       !   
       integer, parameter :: PRCSN = selected_real_kind(6, 37)
       real(PRCSN), intent(inout) :: array_out(:) 
       real(PRCSN), intent(in)    :: array_in(:) 
       integer, optional, intent(in) ::  range1(2)
       integer, optional, intent(in) ::  lbound1
#if defined(__CUDA)
       attributes(device) :: array_out
#endif
       !
    end subroutine sp_h2d_memsync_r1d
    !
    subroutine sp_h2d_memsync_r2d(array_out, array_in, &
                                             range1, lbound1 , &
                                             range2, lbound2  )
#if defined(__CUDA)
       use cudafor
#endif
       implicit none
       !   
       integer, parameter :: PRCSN = selected_real_kind(6, 37)
       real(PRCSN), intent(inout) :: array_out(:,:) 
       real(PRCSN), intent(in)    :: array_in(:,:) 
       integer, optional, intent(in) ::  range1(2), range2(2)
       integer, optional, intent(in) ::  lbound1, lbound2
#if defined(__CUDA)
       attributes(device) :: array_out
#endif
       !
    end subroutine sp_h2d_memsync_r2d
    !
    subroutine sp_h2d_memsync_r3d(array_out, array_in, &
                                             range1, lbound1 , &
                                             range2, lbound2 , &
                                             range3, lbound3  )
#if defined(__CUDA)
       use cudafor
#endif
       implicit none
       !   
       integer, parameter :: PRCSN = selected_real_kind(6, 37)
       real(PRCSN), intent(inout) :: array_out(:,:,:) 
       real(PRCSN), intent(in)    :: array_in(:,:,:) 
       integer, optional, intent(in) ::  range1(2), range2(2), range3(2)
       integer, optional, intent(in) ::  lbound1, lbound2, lbound3
#if defined(__CUDA)
       attributes(device) :: array_out
#endif
       !
    end subroutine sp_h2d_memsync_r3d
    !
    subroutine sp_h2d_memsync_r4d(array_out, array_in, &
                                             range1, lbound1 , &
                                             range2, lbound2 , &
                                             range3, lbound3 , &
                                             range4, lbound4  )
#if defined(__CUDA)
       use cudafor
#endif
       implicit none
       !   
       integer, parameter :: PRCSN = selected_real_kind(6, 37)
       real(PRCSN), intent(inout) :: array_out(:,:,:,:) 
       real(PRCSN), intent(in)    :: array_in(:,:,:,:) 
       integer, optional, intent(in) ::  range1(2), range2(2), range3(2), range4(2)
       integer, optional, intent(in) ::  lbound1, lbound2, lbound3, lbound4
#if defined(__CUDA)
       attributes(device) :: array_out
#endif
       !
    end subroutine sp_h2d_memsync_r4d
    !
    subroutine sp_h2d_memsync_c1d(array_out, array_in, &
                                             range1, lbound1  )
#if defined(__CUDA)
       use cudafor
#endif
       implicit none
       !   
       integer, parameter :: PRCSN = selected_real_kind(6, 37)
       complex(PRCSN), intent(inout) :: array_out(:) 
       complex(PRCSN), intent(in)    :: array_in(:) 
       integer, optional, intent(in) ::  range1(2)
       integer, optional, intent(in) ::  lbound1
#if defined(__CUDA)
       attributes(device) :: array_out
#endif
       !
    end subroutine sp_h2d_memsync_c1d
    !
    subroutine sp_h2d_memsync_c2d(array_out, array_in, &
                                             range1, lbound1 , &
                                             range2, lbound2  )
#if defined(__CUDA)
       use cudafor
#endif
       implicit none
       !   
       integer, parameter :: PRCSN = selected_real_kind(6, 37)
       complex(PRCSN), intent(inout) :: array_out(:,:) 
       complex(PRCSN), intent(in)    :: array_in(:,:) 
       integer, optional, intent(in) ::  range1(2), range2(2)
       integer, optional, intent(in) ::  lbound1, lbound2
#if defined(__CUDA)
       attributes(device) :: array_out
#endif
       !
    end subroutine sp_h2d_memsync_c2d
    !
    subroutine sp_h2d_memsync_c3d(array_out, array_in, &
                                             range1, lbound1 , &
                                             range2, lbound2 , &
                                             range3, lbound3  )
#if defined(__CUDA)
       use cudafor
#endif
       implicit none
       !   
       integer, parameter :: PRCSN = selected_real_kind(6, 37)
       complex(PRCSN), intent(inout) :: array_out(:,:,:) 
       complex(PRCSN), intent(in)    :: array_in(:,:,:) 
       integer, optional, intent(in) ::  range1(2), range2(2), range3(2)
       integer, optional, intent(in) ::  lbound1, lbound2, lbound3
#if defined(__CUDA)
       attributes(device) :: array_out
#endif
       !
    end subroutine sp_h2d_memsync_c3d
    !
    subroutine sp_h2d_memsync_c4d(array_out, array_in, &
                                             range1, lbound1 , &
                                             range2, lbound2 , &
                                             range3, lbound3 , &
                                             range4, lbound4  )
#if defined(__CUDA)
       use cudafor
#endif
       implicit none
       !   
       integer, parameter :: PRCSN = selected_real_kind(6, 37)
       complex(PRCSN), intent(inout) :: array_out(:,:,:,:) 
       complex(PRCSN), intent(in)    :: array_in(:,:,:,:) 
       integer, optional, intent(in) ::  range1(2), range2(2), range3(2), range4(2)
       integer, optional, intent(in) ::  lbound1, lbound2, lbound3, lbound4
#if defined(__CUDA)
       attributes(device) :: array_out
#endif
       !
    end subroutine sp_h2d_memsync_c4d
    !
    !
#if defined(__HAVE_DEVICE)
    subroutine dp_d2h_memsync_r1d(array_out, array_in, &
                                             range1, lbound1  )
#if defined(__CUDA)
       use cudafor
#endif 
       implicit none
       !
       integer, parameter :: PRCSN = selected_real_kind(14,200)
       real(PRCSN), intent(inout) :: array_out(:)
       real(PRCSN), intent(in)    :: array_in(:)
       integer, optional, intent(in) ::  range1(2)
       integer, optional, intent(in) ::  lbound1
#if defined(__CUDA)
       attributes(device) :: array_in
#endif
       !
    end subroutine dp_d2h_memsync_r1d
    !
    subroutine dp_d2h_memsync_r2d(array_out, array_in, &
                                             range1, lbound1 , &
                                             range2, lbound2  )
#if defined(__CUDA)
       use cudafor
#endif 
       implicit none
       !
       integer, parameter :: PRCSN = selected_real_kind(14,200)
       real(PRCSN), intent(inout) :: array_out(:,:)
       real(PRCSN), intent(in)    :: array_in(:,:)
       integer, optional, intent(in) ::  range1(2), range2(2)
       integer, optional, intent(in) ::  lbound1, lbound2
#if defined(__CUDA)
       attributes(device) :: array_in
#endif
       !
    end subroutine dp_d2h_memsync_r2d
    !
    subroutine dp_d2h_memsync_r3d(array_out, array_in, &
                                             range1, lbound1 , &
                                             range2, lbound2 , &
                                             range3, lbound3  )
#if defined(__CUDA)
       use cudafor
#endif 
       implicit none
       !
       integer, parameter :: PRCSN = selected_real_kind(14,200)
       real(PRCSN), intent(inout) :: array_out(:,:,:)
       real(PRCSN), intent(in)    :: array_in(:,:,:)
       integer, optional, intent(in) ::  range1(2), range2(2), range3(2)
       integer, optional, intent(in) ::  lbound1, lbound2, lbound3
#if defined(__CUDA)
       attributes(device) :: array_in
#endif
       !
    end subroutine dp_d2h_memsync_r3d
    !
    subroutine dp_d2h_memsync_r4d(array_out, array_in, &
                                             range1, lbound1 , &
                                             range2, lbound2 , &
                                             range3, lbound3 , &
                                             range4, lbound4  )
#if defined(__CUDA)
       use cudafor
#endif 
       implicit none
       !
       integer, parameter :: PRCSN = selected_real_kind(14,200)
       real(PRCSN), intent(inout) :: array_out(:,:,:,:)
       real(PRCSN), intent(in)    :: array_in(:,:,:,:)
       integer, optional, intent(in) ::  range1(2), range2(2), range3(2), range4(2)
       integer, optional, intent(in) ::  lbound1, lbound2, lbound3, lbound4
#if defined(__CUDA)
       attributes(device) :: array_in
#endif
       !
    end subroutine dp_d2h_memsync_r4d
    !
    subroutine dp_d2h_memsync_c1d(array_out, array_in, &
                                             range1, lbound1  )
#if defined(__CUDA)
       use cudafor
#endif 
       implicit none
       !
       integer, parameter :: PRCSN = selected_real_kind(14,200)
       complex(PRCSN), intent(inout) :: array_out(:)
       complex(PRCSN), intent(in)    :: array_in(:)
       integer, optional, intent(in) ::  range1(2)
       integer, optional, intent(in) ::  lbound1
#if defined(__CUDA)
       attributes(device) :: array_in
#endif
       !
    end subroutine dp_d2h_memsync_c1d
    !
    subroutine dp_d2h_memsync_c2d(array_out, array_in, &
                                             range1, lbound1 , &
                                             range2, lbound2  )
#if defined(__CUDA)
       use cudafor
#endif 
       implicit none
       !
       integer, parameter :: PRCSN = selected_real_kind(14,200)
       complex(PRCSN), intent(inout) :: array_out(:,:)
       complex(PRCSN), intent(in)    :: array_in(:,:)
       integer, optional, intent(in) ::  range1(2), range2(2)
       integer, optional, intent(in) ::  lbound1, lbound2
#if defined(__CUDA)
       attributes(device) :: array_in
#endif
       !
    end subroutine dp_d2h_memsync_c2d
    !
    subroutine dp_d2h_memsync_c3d(array_out, array_in, &
                                             range1, lbound1 , &
                                             range2, lbound2 , &
                                             range3, lbound3  )
#if defined(__CUDA)
       use cudafor
#endif 
       implicit none
       !
       integer, parameter :: PRCSN = selected_real_kind(14,200)
       complex(PRCSN), intent(inout) :: array_out(:,:,:)
       complex(PRCSN), intent(in)    :: array_in(:,:,:)
       integer, optional, intent(in) ::  range1(2), range2(2), range3(2)
       integer, optional, intent(in) ::  lbound1, lbound2, lbound3
#if defined(__CUDA)
       attributes(device) :: array_in
#endif
       !
    end subroutine dp_d2h_memsync_c3d
    !
    subroutine dp_d2h_memsync_c4d(array_out, array_in, &
                                             range1, lbound1 , &
                                             range2, lbound2 , &
                                             range3, lbound3 , &
                                             range4, lbound4  )
#if defined(__CUDA)
       use cudafor
#endif 
       implicit none
       !
       integer, parameter :: PRCSN = selected_real_kind(14,200)
       complex(PRCSN), intent(inout) :: array_out(:,:,:,:)
       complex(PRCSN), intent(in)    :: array_in(:,:,:,:)
       integer, optional, intent(in) ::  range1(2), range2(2), range3(2), range4(2)
       integer, optional, intent(in) ::  lbound1, lbound2, lbound3, lbound4
#if defined(__CUDA)
       attributes(device) :: array_in
#endif
       !
    end subroutine dp_d2h_memsync_c4d
    !
    subroutine sp_d2h_memsync_r1d(array_out, array_in, &
                                             range1, lbound1  )
#if defined(__CUDA)
       use cudafor
#endif 
       implicit none
       !
       integer, parameter :: PRCSN = selected_real_kind(6, 37)
       real(PRCSN), intent(inout) :: array_out(:)
       real(PRCSN), intent(in)    :: array_in(:)
       integer, optional, intent(in) ::  range1(2)
       integer, optional, intent(in) ::  lbound1
#if defined(__CUDA)
       attributes(device) :: array_in
#endif
       !
    end subroutine sp_d2h_memsync_r1d
    !
    subroutine sp_d2h_memsync_r2d(array_out, array_in, &
                                             range1, lbound1 , &
                                             range2, lbound2  )
#if defined(__CUDA)
       use cudafor
#endif 
       implicit none
       !
       integer, parameter :: PRCSN = selected_real_kind(6, 37)
       real(PRCSN), intent(inout) :: array_out(:,:)
       real(PRCSN), intent(in)    :: array_in(:,:)
       integer, optional, intent(in) ::  range1(2), range2(2)
       integer, optional, intent(in) ::  lbound1, lbound2
#if defined(__CUDA)
       attributes(device) :: array_in
#endif
       !
    end subroutine sp_d2h_memsync_r2d
    !
    subroutine sp_d2h_memsync_r3d(array_out, array_in, &
                                             range1, lbound1 , &
                                             range2, lbound2 , &
                                             range3, lbound3  )
#if defined(__CUDA)
       use cudafor
#endif 
       implicit none
       !
       integer, parameter :: PRCSN = selected_real_kind(6, 37)
       real(PRCSN), intent(inout) :: array_out(:,:,:)
       real(PRCSN), intent(in)    :: array_in(:,:,:)
       integer, optional, intent(in) ::  range1(2), range2(2), range3(2)
       integer, optional, intent(in) ::  lbound1, lbound2, lbound3
#if defined(__CUDA)
       attributes(device) :: array_in
#endif
       !
    end subroutine sp_d2h_memsync_r3d
    !
    subroutine sp_d2h_memsync_r4d(array_out, array_in, &
                                             range1, lbound1 , &
                                             range2, lbound2 , &
                                             range3, lbound3 , &
                                             range4, lbound4  )
#if defined(__CUDA)
       use cudafor
#endif 
       implicit none
       !
       integer, parameter :: PRCSN = selected_real_kind(6, 37)
       real(PRCSN), intent(inout) :: array_out(:,:,:,:)
       real(PRCSN), intent(in)    :: array_in(:,:,:,:)
       integer, optional, intent(in) ::  range1(2), range2(2), range3(2), range4(2)
       integer, optional, intent(in) ::  lbound1, lbound2, lbound3, lbound4
#if defined(__CUDA)
       attributes(device) :: array_in
#endif
       !
    end subroutine sp_d2h_memsync_r4d
    !
    subroutine sp_d2h_memsync_c1d(array_out, array_in, &
                                             range1, lbound1  )
#if defined(__CUDA)
       use cudafor
#endif 
       implicit none
       !
       integer, parameter :: PRCSN = selected_real_kind(6, 37)
       complex(PRCSN), intent(inout) :: array_out(:)
       complex(PRCSN), intent(in)    :: array_in(:)
       integer, optional, intent(in) ::  range1(2)
       integer, optional, intent(in) ::  lbound1
#if defined(__CUDA)
       attributes(device) :: array_in
#endif
       !
    end subroutine sp_d2h_memsync_c1d
    !
    subroutine sp_d2h_memsync_c2d(array_out, array_in, &
                                             range1, lbound1 , &
                                             range2, lbound2  )
#if defined(__CUDA)
       use cudafor
#endif 
       implicit none
       !
       integer, parameter :: PRCSN = selected_real_kind(6, 37)
       complex(PRCSN), intent(inout) :: array_out(:,:)
       complex(PRCSN), intent(in)    :: array_in(:,:)
       integer, optional, intent(in) ::  range1(2), range2(2)
       integer, optional, intent(in) ::  lbound1, lbound2
#if defined(__CUDA)
       attributes(device) :: array_in
#endif
       !
    end subroutine sp_d2h_memsync_c2d
    !
    subroutine sp_d2h_memsync_c3d(array_out, array_in, &
                                             range1, lbound1 , &
                                             range2, lbound2 , &
                                             range3, lbound3  )
#if defined(__CUDA)
       use cudafor
#endif 
       implicit none
       !
       integer, parameter :: PRCSN = selected_real_kind(6, 37)
       complex(PRCSN), intent(inout) :: array_out(:,:,:)
       complex(PRCSN), intent(in)    :: array_in(:,:,:)
       integer, optional, intent(in) ::  range1(2), range2(2), range3(2)
       integer, optional, intent(in) ::  lbound1, lbound2, lbound3
#if defined(__CUDA)
       attributes(device) :: array_in
#endif
       !
    end subroutine sp_d2h_memsync_c3d
    !
    subroutine sp_d2h_memsync_c4d(array_out, array_in, &
                                             range1, lbound1 , &
                                             range2, lbound2 , &
                                             range3, lbound3 , &
                                             range4, lbound4  )
#if defined(__CUDA)
       use cudafor
#endif 
       implicit none
       !
       integer, parameter :: PRCSN = selected_real_kind(6, 37)
       complex(PRCSN), intent(inout) :: array_out(:,:,:,:)
       complex(PRCSN), intent(in)    :: array_in(:,:,:,:)
       integer, optional, intent(in) ::  range1(2), range2(2), range3(2), range4(2)
       integer, optional, intent(in) ::  lbound1, lbound2, lbound3, lbound4
#if defined(__CUDA)
       attributes(device) :: array_in
#endif
       !
    end subroutine sp_d2h_memsync_c4d
    !
#endif
    !
end interface dev_memsync
!
!
interface dev_conjg
    !

    subroutine dp_dev_conjg_c1d(array_inout, &
                                   
                                   range1, lbound1 )
       implicit none
       !
       integer, parameter :: PRCSN = selected_real_kind(14,200)
       COMPLEX(PRCSN), intent(inout) :: array_inout(:)
       integer, optional, intent(in) ::  range1(2)
       integer, optional, intent(in) ::  lbound1
#if defined(__CUDA)
       attributes(device) :: array_inout
#endif
       !
    end subroutine dp_dev_conjg_c1d
    !
    subroutine dp_dev_conjg_c2d(array_inout, &
                                   
                                   range1, lbound1, &
                                   range2, lbound2 )
       implicit none
       !
       integer, parameter :: PRCSN = selected_real_kind(14,200)
       COMPLEX(PRCSN), intent(inout) :: array_inout(:,:)
       integer, optional, intent(in) ::  range1(2), range2(2)
       integer, optional, intent(in) ::  lbound1, lbound2
#if defined(__CUDA)
       attributes(device) :: array_inout
#endif
       !
    end subroutine dp_dev_conjg_c2d
    !
    subroutine dp_dev_conjg_c3d(array_inout, &
                                   
                                   range1, lbound1, &
                                   range2, lbound2, &
                                   range3, lbound3 )
       implicit none
       !
       integer, parameter :: PRCSN = selected_real_kind(14,200)
       COMPLEX(PRCSN), intent(inout) :: array_inout(:,:,:)
       integer, optional, intent(in) ::  range1(2), range2(2), range3(2)
       integer, optional, intent(in) ::  lbound1, lbound2, lbound3
#if defined(__CUDA)
       attributes(device) :: array_inout
#endif
       !
    end subroutine dp_dev_conjg_c3d
    !
    subroutine dp_dev_conjg_c4d(array_inout, &
                                   
                                   range1, lbound1, &
                                   range2, lbound2, &
                                   range3, lbound3, &
                                   range4, lbound4 )
       implicit none
       !
       integer, parameter :: PRCSN = selected_real_kind(14,200)
       COMPLEX(PRCSN), intent(inout) :: array_inout(:,:,:,:)
       integer, optional, intent(in) ::  range1(2), range2(2), range3(2), range4(2)
       integer, optional, intent(in) ::  lbound1, lbound2, lbound3, lbound4
#if defined(__CUDA)
       attributes(device) :: array_inout
#endif
       !
    end subroutine dp_dev_conjg_c4d
    !
    subroutine sp_dev_conjg_c1d(array_inout, &
                                   
                                   range1, lbound1 )
       implicit none
       !
       integer, parameter :: PRCSN = selected_real_kind(6, 37)
       COMPLEX(PRCSN), intent(inout) :: array_inout(:)
       integer, optional, intent(in) ::  range1(2)
       integer, optional, intent(in) ::  lbound1
#if defined(__CUDA)
       attributes(device) :: array_inout
#endif
       !
    end subroutine sp_dev_conjg_c1d
    !
    subroutine sp_dev_conjg_c2d(array_inout, &
                                   
                                   range1, lbound1, &
                                   range2, lbound2 )
       implicit none
       !
       integer, parameter :: PRCSN = selected_real_kind(6, 37)
       COMPLEX(PRCSN), intent(inout) :: array_inout(:,:)
       integer, optional, intent(in) ::  range1(2), range2(2)
       integer, optional, intent(in) ::  lbound1, lbound2
#if defined(__CUDA)
       attributes(device) :: array_inout
#endif
       !
    end subroutine sp_dev_conjg_c2d
    !
    subroutine sp_dev_conjg_c3d(array_inout, &
                                   
                                   range1, lbound1, &
                                   range2, lbound2, &
                                   range3, lbound3 )
       implicit none
       !
       integer, parameter :: PRCSN = selected_real_kind(6, 37)
       COMPLEX(PRCSN), intent(inout) :: array_inout(:,:,:)
       integer, optional, intent(in) ::  range1(2), range2(2), range3(2)
       integer, optional, intent(in) ::  lbound1, lbound2, lbound3
#if defined(__CUDA)
       attributes(device) :: array_inout
#endif
       !
    end subroutine sp_dev_conjg_c3d
    !
    subroutine sp_dev_conjg_c4d(array_inout, &
                                   
                                   range1, lbound1, &
                                   range2, lbound2, &
                                   range3, lbound3, &
                                   range4, lbound4 )
       implicit none
       !
       integer, parameter :: PRCSN = selected_real_kind(6, 37)
       COMPLEX(PRCSN), intent(inout) :: array_inout(:,:,:,:)
       integer, optional, intent(in) ::  range1(2), range2(2), range3(2), range4(2)
       integer, optional, intent(in) ::  lbound1, lbound2, lbound3, lbound4
#if defined(__CUDA)
       attributes(device) :: array_inout
#endif
       !
    end subroutine sp_dev_conjg_c4d
    !
    !
end interface dev_conjg

end module device_util_m

!
!=======================================
!
subroutine dp_dev_memcpy_r1d(array_out, array_in, &
                                            range1, lbound1 )
    implicit none
    !
    integer, parameter :: PRCSN = selected_real_kind(14,200)
    real(PRCSN), intent(inout) :: array_out(:)
    real(PRCSN), intent(in)    :: array_in(:)
    integer, optional, intent(in) ::  range1(2)
    integer, optional, intent(in) ::  lbound1
#if defined(__CUDA)
    attributes(device) :: array_out, array_in
#endif
    !
    integer :: i1, d1s, d1e
    integer :: lbound1_, range1_(2)
    !
    lbound1_=1
    if (present(lbound1)) lbound1_=lbound1 
    range1_=(/1,size(array_out, 1)/)
    if (present(range1)) range1_=range1 
    !
    ! the lower bound of the assumed shape array passed to the subroutine is 1
    ! lbound and range instead refer to the indexing in the parent caller.
    d1s = range1_(1) -lbound1_ +1
    d1e = range1_(2) -lbound1_ +1
    !
    !
    !$cuf kernel do(1)
    do i1 = d1s, d1e
        array_out(i1 ) = array_in(i1 )
    enddo
    !
end subroutine dp_dev_memcpy_r1d
!
subroutine dp_dev_memcpy_r2d(array_out, array_in, &
                                            range1, lbound1, &
                                            range2, lbound2 )
    implicit none
    !
    integer, parameter :: PRCSN = selected_real_kind(14,200)
    real(PRCSN), intent(inout) :: array_out(:,:)
    real(PRCSN), intent(in)    :: array_in(:,:)
    integer, optional, intent(in) ::  range1(2), range2(2)
    integer, optional, intent(in) ::  lbound1, lbound2
#if defined(__CUDA)
    attributes(device) :: array_out, array_in
#endif
    !
    integer :: i1, d1s, d1e
    integer :: lbound1_, range1_(2)
    integer :: i2, d2s, d2e
    integer :: lbound2_, range2_(2)
    !
    lbound1_=1
    if (present(lbound1)) lbound1_=lbound1 
    range1_=(/1,size(array_out, 1)/)
    if (present(range1)) range1_=range1 
    !
    ! the lower bound of the assumed shape array passed to the subroutine is 1
    ! lbound and range instead refer to the indexing in the parent caller.
    d1s = range1_(1) -lbound1_ +1
    d1e = range1_(2) -lbound1_ +1
    !
    lbound2_=1
    if (present(lbound2)) lbound2_=lbound2 
    range2_=(/1,size(array_out, 2)/)
    if (present(range2)) range2_=range2 
    !
    ! the lower bound of the assumed shape array passed to the subroutine is 1
    ! lbound and range instead refer to the indexing in the parent caller.
    d2s = range2_(1) -lbound2_ +1
    d2e = range2_(2) -lbound2_ +1
    !
    !
    !$cuf kernel do(2)
    do i2 = d2s, d2e
    do i1 = d1s, d1e
        array_out(i1,i2 ) = array_in(i1,i2 )
    enddo
    enddo
    !
end subroutine dp_dev_memcpy_r2d
!
subroutine dp_dev_memcpy_r3d(array_out, array_in, &
                                            range1, lbound1, &
                                            range2, lbound2, &
                                            range3, lbound3 )
    implicit none
    !
    integer, parameter :: PRCSN = selected_real_kind(14,200)
    real(PRCSN), intent(inout) :: array_out(:,:,:)
    real(PRCSN), intent(in)    :: array_in(:,:,:)
    integer, optional, intent(in) ::  range1(2), range2(2), range3(2)
    integer, optional, intent(in) ::  lbound1, lbound2, lbound3
#if defined(__CUDA)
    attributes(device) :: array_out, array_in
#endif
    !
    integer :: i1, d1s, d1e
    integer :: lbound1_, range1_(2)
    integer :: i2, d2s, d2e
    integer :: lbound2_, range2_(2)
    integer :: i3, d3s, d3e
    integer :: lbound3_, range3_(2)
    !
    lbound1_=1
    if (present(lbound1)) lbound1_=lbound1 
    range1_=(/1,size(array_out, 1)/)
    if (present(range1)) range1_=range1 
    !
    ! the lower bound of the assumed shape array passed to the subroutine is 1
    ! lbound and range instead refer to the indexing in the parent caller.
    d1s = range1_(1) -lbound1_ +1
    d1e = range1_(2) -lbound1_ +1
    !
    lbound2_=1
    if (present(lbound2)) lbound2_=lbound2 
    range2_=(/1,size(array_out, 2)/)
    if (present(range2)) range2_=range2 
    !
    ! the lower bound of the assumed shape array passed to the subroutine is 1
    ! lbound and range instead refer to the indexing in the parent caller.
    d2s = range2_(1) -lbound2_ +1
    d2e = range2_(2) -lbound2_ +1
    !
    lbound3_=1
    if (present(lbound3)) lbound3_=lbound3 
    range3_=(/1,size(array_out, 3)/)
    if (present(range3)) range3_=range3 
    !
    ! the lower bound of the assumed shape array passed to the subroutine is 1
    ! lbound and range instead refer to the indexing in the parent caller.
    d3s = range3_(1) -lbound3_ +1
    d3e = range3_(2) -lbound3_ +1
    !
    !
    !$cuf kernel do(3)
    do i3 = d3s, d3e
    do i2 = d2s, d2e
    do i1 = d1s, d1e
        array_out(i1,i2,i3 ) = array_in(i1,i2,i3 )
    enddo
    enddo
    enddo
    !
end subroutine dp_dev_memcpy_r3d
!
subroutine dp_dev_memcpy_r4d(array_out, array_in, &
                                            range1, lbound1, &
                                            range2, lbound2, &
                                            range3, lbound3, &
                                            range4, lbound4 )
    implicit none
    !
    integer, parameter :: PRCSN = selected_real_kind(14,200)
    real(PRCSN), intent(inout) :: array_out(:,:,:,:)
    real(PRCSN), intent(in)    :: array_in(:,:,:,:)
    integer, optional, intent(in) ::  range1(2), range2(2), range3(2), range4(2)
    integer, optional, intent(in) ::  lbound1, lbound2, lbound3, lbound4
#if defined(__CUDA)
    attributes(device) :: array_out, array_in
#endif
    !
    integer :: i1, d1s, d1e
    integer :: lbound1_, range1_(2)
    integer :: i2, d2s, d2e
    integer :: lbound2_, range2_(2)
    integer :: i3, d3s, d3e
    integer :: lbound3_, range3_(2)
    integer :: i4, d4s, d4e
    integer :: lbound4_, range4_(2)
    !
    lbound1_=1
    if (present(lbound1)) lbound1_=lbound1 
    range1_=(/1,size(array_out, 1)/)
    if (present(range1)) range1_=range1 
    !
    ! the lower bound of the assumed shape array passed to the subroutine is 1
    ! lbound and range instead refer to the indexing in the parent caller.
    d1s = range1_(1) -lbound1_ +1
    d1e = range1_(2) -lbound1_ +1
    !
    lbound2_=1
    if (present(lbound2)) lbound2_=lbound2 
    range2_=(/1,size(array_out, 2)/)
    if (present(range2)) range2_=range2 
    !
    ! the lower bound of the assumed shape array passed to the subroutine is 1
    ! lbound and range instead refer to the indexing in the parent caller.
    d2s = range2_(1) -lbound2_ +1
    d2e = range2_(2) -lbound2_ +1
    !
    lbound3_=1
    if (present(lbound3)) lbound3_=lbound3 
    range3_=(/1,size(array_out, 3)/)
    if (present(range3)) range3_=range3 
    !
    ! the lower bound of the assumed shape array passed to the subroutine is 1
    ! lbound and range instead refer to the indexing in the parent caller.
    d3s = range3_(1) -lbound3_ +1
    d3e = range3_(2) -lbound3_ +1
    !
    lbound4_=1
    if (present(lbound4)) lbound4_=lbound4 
    range4_=(/1,size(array_out, 4)/)
    if (present(range4)) range4_=range4 
    !
    ! the lower bound of the assumed shape array passed to the subroutine is 1
    ! lbound and range instead refer to the indexing in the parent caller.
    d4s = range4_(1) -lbound4_ +1
    d4e = range4_(2) -lbound4_ +1
    !
    !
    !$cuf kernel do(4)
    do i4 = d4s, d4e
    do i3 = d3s, d3e
    do i2 = d2s, d2e
    do i1 = d1s, d1e
        array_out(i1,i2,i3,i4 ) = array_in(i1,i2,i3,i4 )
    enddo
    enddo
    enddo
    enddo
    !
end subroutine dp_dev_memcpy_r4d
!
subroutine dp_dev_memcpy_c1d(array_out, array_in, &
                                            range1, lbound1 )
    implicit none
    !
    integer, parameter :: PRCSN = selected_real_kind(14,200)
    complex(PRCSN), intent(inout) :: array_out(:)
    complex(PRCSN), intent(in)    :: array_in(:)
    integer, optional, intent(in) ::  range1(2)
    integer, optional, intent(in) ::  lbound1
#if defined(__CUDA)
    attributes(device) :: array_out, array_in
#endif
    !
    integer :: i1, d1s, d1e
    integer :: lbound1_, range1_(2)
    !
    lbound1_=1
    if (present(lbound1)) lbound1_=lbound1 
    range1_=(/1,size(array_out, 1)/)
    if (present(range1)) range1_=range1 
    !
    ! the lower bound of the assumed shape array passed to the subroutine is 1
    ! lbound and range instead refer to the indexing in the parent caller.
    d1s = range1_(1) -lbound1_ +1
    d1e = range1_(2) -lbound1_ +1
    !
    !
    !$cuf kernel do(1)
    do i1 = d1s, d1e
        array_out(i1 ) = array_in(i1 )
    enddo
    !
end subroutine dp_dev_memcpy_c1d
!
subroutine dp_dev_memcpy_c2d(array_out, array_in, &
                                            range1, lbound1, &
                                            range2, lbound2 )
    implicit none
    !
    integer, parameter :: PRCSN = selected_real_kind(14,200)
    complex(PRCSN), intent(inout) :: array_out(:,:)
    complex(PRCSN), intent(in)    :: array_in(:,:)
    integer, optional, intent(in) ::  range1(2), range2(2)
    integer, optional, intent(in) ::  lbound1, lbound2
#if defined(__CUDA)
    attributes(device) :: array_out, array_in
#endif
    !
    integer :: i1, d1s, d1e
    integer :: lbound1_, range1_(2)
    integer :: i2, d2s, d2e
    integer :: lbound2_, range2_(2)
    !
    lbound1_=1
    if (present(lbound1)) lbound1_=lbound1 
    range1_=(/1,size(array_out, 1)/)
    if (present(range1)) range1_=range1 
    !
    ! the lower bound of the assumed shape array passed to the subroutine is 1
    ! lbound and range instead refer to the indexing in the parent caller.
    d1s = range1_(1) -lbound1_ +1
    d1e = range1_(2) -lbound1_ +1
    !
    lbound2_=1
    if (present(lbound2)) lbound2_=lbound2 
    range2_=(/1,size(array_out, 2)/)
    if (present(range2)) range2_=range2 
    !
    ! the lower bound of the assumed shape array passed to the subroutine is 1
    ! lbound and range instead refer to the indexing in the parent caller.
    d2s = range2_(1) -lbound2_ +1
    d2e = range2_(2) -lbound2_ +1
    !
    !
    !$cuf kernel do(2)
    do i2 = d2s, d2e
    do i1 = d1s, d1e
        array_out(i1,i2 ) = array_in(i1,i2 )
    enddo
    enddo
    !
end subroutine dp_dev_memcpy_c2d
!
subroutine dp_dev_memcpy_c3d(array_out, array_in, &
                                            range1, lbound1, &
                                            range2, lbound2, &
                                            range3, lbound3 )
    implicit none
    !
    integer, parameter :: PRCSN = selected_real_kind(14,200)
    complex(PRCSN), intent(inout) :: array_out(:,:,:)
    complex(PRCSN), intent(in)    :: array_in(:,:,:)
    integer, optional, intent(in) ::  range1(2), range2(2), range3(2)
    integer, optional, intent(in) ::  lbound1, lbound2, lbound3
#if defined(__CUDA)
    attributes(device) :: array_out, array_in
#endif
    !
    integer :: i1, d1s, d1e
    integer :: lbound1_, range1_(2)
    integer :: i2, d2s, d2e
    integer :: lbound2_, range2_(2)
    integer :: i3, d3s, d3e
    integer :: lbound3_, range3_(2)
    !
    lbound1_=1
    if (present(lbound1)) lbound1_=lbound1 
    range1_=(/1,size(array_out, 1)/)
    if (present(range1)) range1_=range1 
    !
    ! the lower bound of the assumed shape array passed to the subroutine is 1
    ! lbound and range instead refer to the indexing in the parent caller.
    d1s = range1_(1) -lbound1_ +1
    d1e = range1_(2) -lbound1_ +1
    !
    lbound2_=1
    if (present(lbound2)) lbound2_=lbound2 
    range2_=(/1,size(array_out, 2)/)
    if (present(range2)) range2_=range2 
    !
    ! the lower bound of the assumed shape array passed to the subroutine is 1
    ! lbound and range instead refer to the indexing in the parent caller.
    d2s = range2_(1) -lbound2_ +1
    d2e = range2_(2) -lbound2_ +1
    !
    lbound3_=1
    if (present(lbound3)) lbound3_=lbound3 
    range3_=(/1,size(array_out, 3)/)
    if (present(range3)) range3_=range3 
    !
    ! the lower bound of the assumed shape array passed to the subroutine is 1
    ! lbound and range instead refer to the indexing in the parent caller.
    d3s = range3_(1) -lbound3_ +1
    d3e = range3_(2) -lbound3_ +1
    !
    !
    !$cuf kernel do(3)
    do i3 = d3s, d3e
    do i2 = d2s, d2e
    do i1 = d1s, d1e
        array_out(i1,i2,i3 ) = array_in(i1,i2,i3 )
    enddo
    enddo
    enddo
    !
end subroutine dp_dev_memcpy_c3d
!
subroutine dp_dev_memcpy_c4d(array_out, array_in, &
                                            range1, lbound1, &
                                            range2, lbound2, &
                                            range3, lbound3, &
                                            range4, lbound4 )
    implicit none
    !
    integer, parameter :: PRCSN = selected_real_kind(14,200)
    complex(PRCSN), intent(inout) :: array_out(:,:,:,:)
    complex(PRCSN), intent(in)    :: array_in(:,:,:,:)
    integer, optional, intent(in) ::  range1(2), range2(2), range3(2), range4(2)
    integer, optional, intent(in) ::  lbound1, lbound2, lbound3, lbound4
#if defined(__CUDA)
    attributes(device) :: array_out, array_in
#endif
    !
    integer :: i1, d1s, d1e
    integer :: lbound1_, range1_(2)
    integer :: i2, d2s, d2e
    integer :: lbound2_, range2_(2)
    integer :: i3, d3s, d3e
    integer :: lbound3_, range3_(2)
    integer :: i4, d4s, d4e
    integer :: lbound4_, range4_(2)
    !
    lbound1_=1
    if (present(lbound1)) lbound1_=lbound1 
    range1_=(/1,size(array_out, 1)/)
    if (present(range1)) range1_=range1 
    !
    ! the lower bound of the assumed shape array passed to the subroutine is 1
    ! lbound and range instead refer to the indexing in the parent caller.
    d1s = range1_(1) -lbound1_ +1
    d1e = range1_(2) -lbound1_ +1
    !
    lbound2_=1
    if (present(lbound2)) lbound2_=lbound2 
    range2_=(/1,size(array_out, 2)/)
    if (present(range2)) range2_=range2 
    !
    ! the lower bound of the assumed shape array passed to the subroutine is 1
    ! lbound and range instead refer to the indexing in the parent caller.
    d2s = range2_(1) -lbound2_ +1
    d2e = range2_(2) -lbound2_ +1
    !
    lbound3_=1
    if (present(lbound3)) lbound3_=lbound3 
    range3_=(/1,size(array_out, 3)/)
    if (present(range3)) range3_=range3 
    !
    ! the lower bound of the assumed shape array passed to the subroutine is 1
    ! lbound and range instead refer to the indexing in the parent caller.
    d3s = range3_(1) -lbound3_ +1
    d3e = range3_(2) -lbound3_ +1
    !
    lbound4_=1
    if (present(lbound4)) lbound4_=lbound4 
    range4_=(/1,size(array_out, 4)/)
    if (present(range4)) range4_=range4 
    !
    ! the lower bound of the assumed shape array passed to the subroutine is 1
    ! lbound and range instead refer to the indexing in the parent caller.
    d4s = range4_(1) -lbound4_ +1
    d4e = range4_(2) -lbound4_ +1
    !
    !
    !$cuf kernel do(4)
    do i4 = d4s, d4e
    do i3 = d3s, d3e
    do i2 = d2s, d2e
    do i1 = d1s, d1e
        array_out(i1,i2,i3,i4 ) = array_in(i1,i2,i3,i4 )
    enddo
    enddo
    enddo
    enddo
    !
end subroutine dp_dev_memcpy_c4d
!
subroutine sp_dev_memcpy_r1d(array_out, array_in, &
                                            range1, lbound1 )
    implicit none
    !
    integer, parameter :: PRCSN = selected_real_kind(6, 37)
    real(PRCSN), intent(inout) :: array_out(:)
    real(PRCSN), intent(in)    :: array_in(:)
    integer, optional, intent(in) ::  range1(2)
    integer, optional, intent(in) ::  lbound1
#if defined(__CUDA)
    attributes(device) :: array_out, array_in
#endif
    !
    integer :: i1, d1s, d1e
    integer :: lbound1_, range1_(2)
    !
    lbound1_=1
    if (present(lbound1)) lbound1_=lbound1 
    range1_=(/1,size(array_out, 1)/)
    if (present(range1)) range1_=range1 
    !
    ! the lower bound of the assumed shape array passed to the subroutine is 1
    ! lbound and range instead refer to the indexing in the parent caller.
    d1s = range1_(1) -lbound1_ +1
    d1e = range1_(2) -lbound1_ +1
    !
    !
    !$cuf kernel do(1)
    do i1 = d1s, d1e
        array_out(i1 ) = array_in(i1 )
    enddo
    !
end subroutine sp_dev_memcpy_r1d
!
subroutine sp_dev_memcpy_r2d(array_out, array_in, &
                                            range1, lbound1, &
                                            range2, lbound2 )
    implicit none
    !
    integer, parameter :: PRCSN = selected_real_kind(6, 37)
    real(PRCSN), intent(inout) :: array_out(:,:)
    real(PRCSN), intent(in)    :: array_in(:,:)
    integer, optional, intent(in) ::  range1(2), range2(2)
    integer, optional, intent(in) ::  lbound1, lbound2
#if defined(__CUDA)
    attributes(device) :: array_out, array_in
#endif
    !
    integer :: i1, d1s, d1e
    integer :: lbound1_, range1_(2)
    integer :: i2, d2s, d2e
    integer :: lbound2_, range2_(2)
    !
    lbound1_=1
    if (present(lbound1)) lbound1_=lbound1 
    range1_=(/1,size(array_out, 1)/)
    if (present(range1)) range1_=range1 
    !
    ! the lower bound of the assumed shape array passed to the subroutine is 1
    ! lbound and range instead refer to the indexing in the parent caller.
    d1s = range1_(1) -lbound1_ +1
    d1e = range1_(2) -lbound1_ +1
    !
    lbound2_=1
    if (present(lbound2)) lbound2_=lbound2 
    range2_=(/1,size(array_out, 2)/)
    if (present(range2)) range2_=range2 
    !
    ! the lower bound of the assumed shape array passed to the subroutine is 1
    ! lbound and range instead refer to the indexing in the parent caller.
    d2s = range2_(1) -lbound2_ +1
    d2e = range2_(2) -lbound2_ +1
    !
    !
    !$cuf kernel do(2)
    do i2 = d2s, d2e
    do i1 = d1s, d1e
        array_out(i1,i2 ) = array_in(i1,i2 )
    enddo
    enddo
    !
end subroutine sp_dev_memcpy_r2d
!
subroutine sp_dev_memcpy_r3d(array_out, array_in, &
                                            range1, lbound1, &
                                            range2, lbound2, &
                                            range3, lbound3 )
    implicit none
    !
    integer, parameter :: PRCSN = selected_real_kind(6, 37)
    real(PRCSN), intent(inout) :: array_out(:,:,:)
    real(PRCSN), intent(in)    :: array_in(:,:,:)
    integer, optional, intent(in) ::  range1(2), range2(2), range3(2)
    integer, optional, intent(in) ::  lbound1, lbound2, lbound3
#if defined(__CUDA)
    attributes(device) :: array_out, array_in
#endif
    !
    integer :: i1, d1s, d1e
    integer :: lbound1_, range1_(2)
    integer :: i2, d2s, d2e
    integer :: lbound2_, range2_(2)
    integer :: i3, d3s, d3e
    integer :: lbound3_, range3_(2)
    !
    lbound1_=1
    if (present(lbound1)) lbound1_=lbound1 
    range1_=(/1,size(array_out, 1)/)
    if (present(range1)) range1_=range1 
    !
    ! the lower bound of the assumed shape array passed to the subroutine is 1
    ! lbound and range instead refer to the indexing in the parent caller.
    d1s = range1_(1) -lbound1_ +1
    d1e = range1_(2) -lbound1_ +1
    !
    lbound2_=1
    if (present(lbound2)) lbound2_=lbound2 
    range2_=(/1,size(array_out, 2)/)
    if (present(range2)) range2_=range2 
    !
    ! the lower bound of the assumed shape array passed to the subroutine is 1
    ! lbound and range instead refer to the indexing in the parent caller.
    d2s = range2_(1) -lbound2_ +1
    d2e = range2_(2) -lbound2_ +1
    !
    lbound3_=1
    if (present(lbound3)) lbound3_=lbound3 
    range3_=(/1,size(array_out, 3)/)
    if (present(range3)) range3_=range3 
    !
    ! the lower bound of the assumed shape array passed to the subroutine is 1
    ! lbound and range instead refer to the indexing in the parent caller.
    d3s = range3_(1) -lbound3_ +1
    d3e = range3_(2) -lbound3_ +1
    !
    !
    !$cuf kernel do(3)
    do i3 = d3s, d3e
    do i2 = d2s, d2e
    do i1 = d1s, d1e
        array_out(i1,i2,i3 ) = array_in(i1,i2,i3 )
    enddo
    enddo
    enddo
    !
end subroutine sp_dev_memcpy_r3d
!
subroutine sp_dev_memcpy_r4d(array_out, array_in, &
                                            range1, lbound1, &
                                            range2, lbound2, &
                                            range3, lbound3, &
                                            range4, lbound4 )
    implicit none
    !
    integer, parameter :: PRCSN = selected_real_kind(6, 37)
    real(PRCSN), intent(inout) :: array_out(:,:,:,:)
    real(PRCSN), intent(in)    :: array_in(:,:,:,:)
    integer, optional, intent(in) ::  range1(2), range2(2), range3(2), range4(2)
    integer, optional, intent(in) ::  lbound1, lbound2, lbound3, lbound4
#if defined(__CUDA)
    attributes(device) :: array_out, array_in
#endif
    !
    integer :: i1, d1s, d1e
    integer :: lbound1_, range1_(2)
    integer :: i2, d2s, d2e
    integer :: lbound2_, range2_(2)
    integer :: i3, d3s, d3e
    integer :: lbound3_, range3_(2)
    integer :: i4, d4s, d4e
    integer :: lbound4_, range4_(2)
    !
    lbound1_=1
    if (present(lbound1)) lbound1_=lbound1 
    range1_=(/1,size(array_out, 1)/)
    if (present(range1)) range1_=range1 
    !
    ! the lower bound of the assumed shape array passed to the subroutine is 1
    ! lbound and range instead refer to the indexing in the parent caller.
    d1s = range1_(1) -lbound1_ +1
    d1e = range1_(2) -lbound1_ +1
    !
    lbound2_=1
    if (present(lbound2)) lbound2_=lbound2 
    range2_=(/1,size(array_out, 2)/)
    if (present(range2)) range2_=range2 
    !
    ! the lower bound of the assumed shape array passed to the subroutine is 1
    ! lbound and range instead refer to the indexing in the parent caller.
    d2s = range2_(1) -lbound2_ +1
    d2e = range2_(2) -lbound2_ +1
    !
    lbound3_=1
    if (present(lbound3)) lbound3_=lbound3 
    range3_=(/1,size(array_out, 3)/)
    if (present(range3)) range3_=range3 
    !
    ! the lower bound of the assumed shape array passed to the subroutine is 1
    ! lbound and range instead refer to the indexing in the parent caller.
    d3s = range3_(1) -lbound3_ +1
    d3e = range3_(2) -lbound3_ +1
    !
    lbound4_=1
    if (present(lbound4)) lbound4_=lbound4 
    range4_=(/1,size(array_out, 4)/)
    if (present(range4)) range4_=range4 
    !
    ! the lower bound of the assumed shape array passed to the subroutine is 1
    ! lbound and range instead refer to the indexing in the parent caller.
    d4s = range4_(1) -lbound4_ +1
    d4e = range4_(2) -lbound4_ +1
    !
    !
    !$cuf kernel do(4)
    do i4 = d4s, d4e
    do i3 = d3s, d3e
    do i2 = d2s, d2e
    do i1 = d1s, d1e
        array_out(i1,i2,i3,i4 ) = array_in(i1,i2,i3,i4 )
    enddo
    enddo
    enddo
    enddo
    !
end subroutine sp_dev_memcpy_r4d
!
subroutine sp_dev_memcpy_c1d(array_out, array_in, &
                                            range1, lbound1 )
    implicit none
    !
    integer, parameter :: PRCSN = selected_real_kind(6, 37)
    complex(PRCSN), intent(inout) :: array_out(:)
    complex(PRCSN), intent(in)    :: array_in(:)
    integer, optional, intent(in) ::  range1(2)
    integer, optional, intent(in) ::  lbound1
#if defined(__CUDA)
    attributes(device) :: array_out, array_in
#endif
    !
    integer :: i1, d1s, d1e
    integer :: lbound1_, range1_(2)
    !
    lbound1_=1
    if (present(lbound1)) lbound1_=lbound1 
    range1_=(/1,size(array_out, 1)/)
    if (present(range1)) range1_=range1 
    !
    ! the lower bound of the assumed shape array passed to the subroutine is 1
    ! lbound and range instead refer to the indexing in the parent caller.
    d1s = range1_(1) -lbound1_ +1
    d1e = range1_(2) -lbound1_ +1
    !
    !
    !$cuf kernel do(1)
    do i1 = d1s, d1e
        array_out(i1 ) = array_in(i1 )
    enddo
    !
end subroutine sp_dev_memcpy_c1d
!
subroutine sp_dev_memcpy_c2d(array_out, array_in, &
                                            range1, lbound1, &
                                            range2, lbound2 )
    implicit none
    !
    integer, parameter :: PRCSN = selected_real_kind(6, 37)
    complex(PRCSN), intent(inout) :: array_out(:,:)
    complex(PRCSN), intent(in)    :: array_in(:,:)
    integer, optional, intent(in) ::  range1(2), range2(2)
    integer, optional, intent(in) ::  lbound1, lbound2
#if defined(__CUDA)
    attributes(device) :: array_out, array_in
#endif
    !
    integer :: i1, d1s, d1e
    integer :: lbound1_, range1_(2)
    integer :: i2, d2s, d2e
    integer :: lbound2_, range2_(2)
    !
    lbound1_=1
    if (present(lbound1)) lbound1_=lbound1 
    range1_=(/1,size(array_out, 1)/)
    if (present(range1)) range1_=range1 
    !
    ! the lower bound of the assumed shape array passed to the subroutine is 1
    ! lbound and range instead refer to the indexing in the parent caller.
    d1s = range1_(1) -lbound1_ +1
    d1e = range1_(2) -lbound1_ +1
    !
    lbound2_=1
    if (present(lbound2)) lbound2_=lbound2 
    range2_=(/1,size(array_out, 2)/)
    if (present(range2)) range2_=range2 
    !
    ! the lower bound of the assumed shape array passed to the subroutine is 1
    ! lbound and range instead refer to the indexing in the parent caller.
    d2s = range2_(1) -lbound2_ +1
    d2e = range2_(2) -lbound2_ +1
    !
    !
    !$cuf kernel do(2)
    do i2 = d2s, d2e
    do i1 = d1s, d1e
        array_out(i1,i2 ) = array_in(i1,i2 )
    enddo
    enddo
    !
end subroutine sp_dev_memcpy_c2d
!
subroutine sp_dev_memcpy_c3d(array_out, array_in, &
                                            range1, lbound1, &
                                            range2, lbound2, &
                                            range3, lbound3 )
    implicit none
    !
    integer, parameter :: PRCSN = selected_real_kind(6, 37)
    complex(PRCSN), intent(inout) :: array_out(:,:,:)
    complex(PRCSN), intent(in)    :: array_in(:,:,:)
    integer, optional, intent(in) ::  range1(2), range2(2), range3(2)
    integer, optional, intent(in) ::  lbound1, lbound2, lbound3
#if defined(__CUDA)
    attributes(device) :: array_out, array_in
#endif
    !
    integer :: i1, d1s, d1e
    integer :: lbound1_, range1_(2)
    integer :: i2, d2s, d2e
    integer :: lbound2_, range2_(2)
    integer :: i3, d3s, d3e
    integer :: lbound3_, range3_(2)
    !
    lbound1_=1
    if (present(lbound1)) lbound1_=lbound1 
    range1_=(/1,size(array_out, 1)/)
    if (present(range1)) range1_=range1 
    !
    ! the lower bound of the assumed shape array passed to the subroutine is 1
    ! lbound and range instead refer to the indexing in the parent caller.
    d1s = range1_(1) -lbound1_ +1
    d1e = range1_(2) -lbound1_ +1
    !
    lbound2_=1
    if (present(lbound2)) lbound2_=lbound2 
    range2_=(/1,size(array_out, 2)/)
    if (present(range2)) range2_=range2 
    !
    ! the lower bound of the assumed shape array passed to the subroutine is 1
    ! lbound and range instead refer to the indexing in the parent caller.
    d2s = range2_(1) -lbound2_ +1
    d2e = range2_(2) -lbound2_ +1
    !
    lbound3_=1
    if (present(lbound3)) lbound3_=lbound3 
    range3_=(/1,size(array_out, 3)/)
    if (present(range3)) range3_=range3 
    !
    ! the lower bound of the assumed shape array passed to the subroutine is 1
    ! lbound and range instead refer to the indexing in the parent caller.
    d3s = range3_(1) -lbound3_ +1
    d3e = range3_(2) -lbound3_ +1
    !
    !
    !$cuf kernel do(3)
    do i3 = d3s, d3e
    do i2 = d2s, d2e
    do i1 = d1s, d1e
        array_out(i1,i2,i3 ) = array_in(i1,i2,i3 )
    enddo
    enddo
    enddo
    !
end subroutine sp_dev_memcpy_c3d
!
subroutine sp_dev_memcpy_c4d(array_out, array_in, &
                                            range1, lbound1, &
                                            range2, lbound2, &
                                            range3, lbound3, &
                                            range4, lbound4 )
    implicit none
    !
    integer, parameter :: PRCSN = selected_real_kind(6, 37)
    complex(PRCSN), intent(inout) :: array_out(:,:,:,:)
    complex(PRCSN), intent(in)    :: array_in(:,:,:,:)
    integer, optional, intent(in) ::  range1(2), range2(2), range3(2), range4(2)
    integer, optional, intent(in) ::  lbound1, lbound2, lbound3, lbound4
#if defined(__CUDA)
    attributes(device) :: array_out, array_in
#endif
    !
    integer :: i1, d1s, d1e
    integer :: lbound1_, range1_(2)
    integer :: i2, d2s, d2e
    integer :: lbound2_, range2_(2)
    integer :: i3, d3s, d3e
    integer :: lbound3_, range3_(2)
    integer :: i4, d4s, d4e
    integer :: lbound4_, range4_(2)
    !
    lbound1_=1
    if (present(lbound1)) lbound1_=lbound1 
    range1_=(/1,size(array_out, 1)/)
    if (present(range1)) range1_=range1 
    !
    ! the lower bound of the assumed shape array passed to the subroutine is 1
    ! lbound and range instead refer to the indexing in the parent caller.
    d1s = range1_(1) -lbound1_ +1
    d1e = range1_(2) -lbound1_ +1
    !
    lbound2_=1
    if (present(lbound2)) lbound2_=lbound2 
    range2_=(/1,size(array_out, 2)/)
    if (present(range2)) range2_=range2 
    !
    ! the lower bound of the assumed shape array passed to the subroutine is 1
    ! lbound and range instead refer to the indexing in the parent caller.
    d2s = range2_(1) -lbound2_ +1
    d2e = range2_(2) -lbound2_ +1
    !
    lbound3_=1
    if (present(lbound3)) lbound3_=lbound3 
    range3_=(/1,size(array_out, 3)/)
    if (present(range3)) range3_=range3 
    !
    ! the lower bound of the assumed shape array passed to the subroutine is 1
    ! lbound and range instead refer to the indexing in the parent caller.
    d3s = range3_(1) -lbound3_ +1
    d3e = range3_(2) -lbound3_ +1
    !
    lbound4_=1
    if (present(lbound4)) lbound4_=lbound4 
    range4_=(/1,size(array_out, 4)/)
    if (present(range4)) range4_=range4 
    !
    ! the lower bound of the assumed shape array passed to the subroutine is 1
    ! lbound and range instead refer to the indexing in the parent caller.
    d4s = range4_(1) -lbound4_ +1
    d4e = range4_(2) -lbound4_ +1
    !
    !
    !$cuf kernel do(4)
    do i4 = d4s, d4e
    do i3 = d3s, d3e
    do i2 = d2s, d2e
    do i1 = d1s, d1e
        array_out(i1,i2,i3,i4 ) = array_in(i1,i2,i3,i4 )
    enddo
    enddo
    enddo
    enddo
    !
end subroutine sp_dev_memcpy_c4d
!
!
!======================
!
subroutine dp_dev_memset_r1d(array_out, val, & 
                                             
                                            range1, lbound1 )
    implicit none
    !
    integer, parameter :: PRCSN = selected_real_kind(14,200)
    real(PRCSN), intent(inout) :: array_out(:)
    real(PRCSN), intent(in)    :: val
    integer, optional, intent(in) ::  range1(2)
    integer, optional, intent(in) ::  lbound1
#if defined(__CUDA)
    attributes(device) :: array_out
#endif
    !
    integer :: i1, d1s, d1e
    integer :: lbound1_, range1_(2)
    !
    lbound1_=1
    if (present(lbound1)) lbound1_=lbound1 
    range1_=(/1,size(array_out, 1)/)
    if (present(range1)) range1_=range1 
    !
    d1s = range1_(1) -lbound1_ +1
    d1e = range1_(2) -lbound1_ +1
    !
    !
    !$cuf kernel do(1)
    do i1 = d1s, d1e
        array_out(i1 ) = val
    enddo
    !
end subroutine dp_dev_memset_r1d
!
subroutine dp_dev_memset_r2d(array_out, val, & 
                                             
                                            range1, lbound1, & 
                                            range2, lbound2 )
    implicit none
    !
    integer, parameter :: PRCSN = selected_real_kind(14,200)
    real(PRCSN), intent(inout) :: array_out(:,:)
    real(PRCSN), intent(in)    :: val
    integer, optional, intent(in) ::  range1(2), range2(2)
    integer, optional, intent(in) ::  lbound1, lbound2
#if defined(__CUDA)
    attributes(device) :: array_out
#endif
    !
    integer :: i1, d1s, d1e
    integer :: lbound1_, range1_(2)
    integer :: i2, d2s, d2e
    integer :: lbound2_, range2_(2)
    !
    lbound1_=1
    if (present(lbound1)) lbound1_=lbound1 
    range1_=(/1,size(array_out, 1)/)
    if (present(range1)) range1_=range1 
    !
    d1s = range1_(1) -lbound1_ +1
    d1e = range1_(2) -lbound1_ +1
    !
    lbound2_=1
    if (present(lbound2)) lbound2_=lbound2 
    range2_=(/1,size(array_out, 2)/)
    if (present(range2)) range2_=range2 
    !
    d2s = range2_(1) -lbound2_ +1
    d2e = range2_(2) -lbound2_ +1
    !
    !
    !$cuf kernel do(2)
    do i2 = d2s, d2e
    do i1 = d1s, d1e
        array_out(i1,i2 ) = val
    enddo
    enddo
    !
end subroutine dp_dev_memset_r2d
!
subroutine dp_dev_memset_r3d(array_out, val, & 
                                             
                                            range1, lbound1, & 
                                            range2, lbound2, & 
                                            range3, lbound3 )
    implicit none
    !
    integer, parameter :: PRCSN = selected_real_kind(14,200)
    real(PRCSN), intent(inout) :: array_out(:,:,:)
    real(PRCSN), intent(in)    :: val
    integer, optional, intent(in) ::  range1(2), range2(2), range3(2)
    integer, optional, intent(in) ::  lbound1, lbound2, lbound3
#if defined(__CUDA)
    attributes(device) :: array_out
#endif
    !
    integer :: i1, d1s, d1e
    integer :: lbound1_, range1_(2)
    integer :: i2, d2s, d2e
    integer :: lbound2_, range2_(2)
    integer :: i3, d3s, d3e
    integer :: lbound3_, range3_(2)
    !
    lbound1_=1
    if (present(lbound1)) lbound1_=lbound1 
    range1_=(/1,size(array_out, 1)/)
    if (present(range1)) range1_=range1 
    !
    d1s = range1_(1) -lbound1_ +1
    d1e = range1_(2) -lbound1_ +1
    !
    lbound2_=1
    if (present(lbound2)) lbound2_=lbound2 
    range2_=(/1,size(array_out, 2)/)
    if (present(range2)) range2_=range2 
    !
    d2s = range2_(1) -lbound2_ +1
    d2e = range2_(2) -lbound2_ +1
    !
    lbound3_=1
    if (present(lbound3)) lbound3_=lbound3 
    range3_=(/1,size(array_out, 3)/)
    if (present(range3)) range3_=range3 
    !
    d3s = range3_(1) -lbound3_ +1
    d3e = range3_(2) -lbound3_ +1
    !
    !
    !$cuf kernel do(3)
    do i3 = d3s, d3e
    do i2 = d2s, d2e
    do i1 = d1s, d1e
        array_out(i1,i2,i3 ) = val
    enddo
    enddo
    enddo
    !
end subroutine dp_dev_memset_r3d
!
subroutine dp_dev_memset_r4d(array_out, val, & 
                                             
                                            range1, lbound1, & 
                                            range2, lbound2, & 
                                            range3, lbound3, & 
                                            range4, lbound4 )
    implicit none
    !
    integer, parameter :: PRCSN = selected_real_kind(14,200)
    real(PRCSN), intent(inout) :: array_out(:,:,:,:)
    real(PRCSN), intent(in)    :: val
    integer, optional, intent(in) ::  range1(2), range2(2), range3(2), range4(2)
    integer, optional, intent(in) ::  lbound1, lbound2, lbound3, lbound4
#if defined(__CUDA)
    attributes(device) :: array_out
#endif
    !
    integer :: i1, d1s, d1e
    integer :: lbound1_, range1_(2)
    integer :: i2, d2s, d2e
    integer :: lbound2_, range2_(2)
    integer :: i3, d3s, d3e
    integer :: lbound3_, range3_(2)
    integer :: i4, d4s, d4e
    integer :: lbound4_, range4_(2)
    !
    lbound1_=1
    if (present(lbound1)) lbound1_=lbound1 
    range1_=(/1,size(array_out, 1)/)
    if (present(range1)) range1_=range1 
    !
    d1s = range1_(1) -lbound1_ +1
    d1e = range1_(2) -lbound1_ +1
    !
    lbound2_=1
    if (present(lbound2)) lbound2_=lbound2 
    range2_=(/1,size(array_out, 2)/)
    if (present(range2)) range2_=range2 
    !
    d2s = range2_(1) -lbound2_ +1
    d2e = range2_(2) -lbound2_ +1
    !
    lbound3_=1
    if (present(lbound3)) lbound3_=lbound3 
    range3_=(/1,size(array_out, 3)/)
    if (present(range3)) range3_=range3 
    !
    d3s = range3_(1) -lbound3_ +1
    d3e = range3_(2) -lbound3_ +1
    !
    lbound4_=1
    if (present(lbound4)) lbound4_=lbound4 
    range4_=(/1,size(array_out, 4)/)
    if (present(range4)) range4_=range4 
    !
    d4s = range4_(1) -lbound4_ +1
    d4e = range4_(2) -lbound4_ +1
    !
    !
    !$cuf kernel do(4)
    do i4 = d4s, d4e
    do i3 = d3s, d3e
    do i2 = d2s, d2e
    do i1 = d1s, d1e
        array_out(i1,i2,i3,i4 ) = val
    enddo
    enddo
    enddo
    enddo
    !
end subroutine dp_dev_memset_r4d
!
subroutine dp_dev_memset_c1d(array_out, val, & 
                                             
                                            range1, lbound1 )
    implicit none
    !
    integer, parameter :: PRCSN = selected_real_kind(14,200)
    complex(PRCSN), intent(inout) :: array_out(:)
    complex(PRCSN), intent(in)    :: val
    integer, optional, intent(in) ::  range1(2)
    integer, optional, intent(in) ::  lbound1
#if defined(__CUDA)
    attributes(device) :: array_out
#endif
    !
    integer :: i1, d1s, d1e
    integer :: lbound1_, range1_(2)
    !
    lbound1_=1
    if (present(lbound1)) lbound1_=lbound1 
    range1_=(/1,size(array_out, 1)/)
    if (present(range1)) range1_=range1 
    !
    d1s = range1_(1) -lbound1_ +1
    d1e = range1_(2) -lbound1_ +1
    !
    !
    !$cuf kernel do(1)
    do i1 = d1s, d1e
        array_out(i1 ) = val
    enddo
    !
end subroutine dp_dev_memset_c1d
!
subroutine dp_dev_memset_c2d(array_out, val, & 
                                             
                                            range1, lbound1, & 
                                            range2, lbound2 )
    implicit none
    !
    integer, parameter :: PRCSN = selected_real_kind(14,200)
    complex(PRCSN), intent(inout) :: array_out(:,:)
    complex(PRCSN), intent(in)    :: val
    integer, optional, intent(in) ::  range1(2), range2(2)
    integer, optional, intent(in) ::  lbound1, lbound2
#if defined(__CUDA)
    attributes(device) :: array_out
#endif
    !
    integer :: i1, d1s, d1e
    integer :: lbound1_, range1_(2)
    integer :: i2, d2s, d2e
    integer :: lbound2_, range2_(2)
    !
    lbound1_=1
    if (present(lbound1)) lbound1_=lbound1 
    range1_=(/1,size(array_out, 1)/)
    if (present(range1)) range1_=range1 
    !
    d1s = range1_(1) -lbound1_ +1
    d1e = range1_(2) -lbound1_ +1
    !
    lbound2_=1
    if (present(lbound2)) lbound2_=lbound2 
    range2_=(/1,size(array_out, 2)/)
    if (present(range2)) range2_=range2 
    !
    d2s = range2_(1) -lbound2_ +1
    d2e = range2_(2) -lbound2_ +1
    !
    !
    !$cuf kernel do(2)
    do i2 = d2s, d2e
    do i1 = d1s, d1e
        array_out(i1,i2 ) = val
    enddo
    enddo
    !
end subroutine dp_dev_memset_c2d
!
subroutine dp_dev_memset_c3d(array_out, val, & 
                                             
                                            range1, lbound1, & 
                                            range2, lbound2, & 
                                            range3, lbound3 )
    implicit none
    !
    integer, parameter :: PRCSN = selected_real_kind(14,200)
    complex(PRCSN), intent(inout) :: array_out(:,:,:)
    complex(PRCSN), intent(in)    :: val
    integer, optional, intent(in) ::  range1(2), range2(2), range3(2)
    integer, optional, intent(in) ::  lbound1, lbound2, lbound3
#if defined(__CUDA)
    attributes(device) :: array_out
#endif
    !
    integer :: i1, d1s, d1e
    integer :: lbound1_, range1_(2)
    integer :: i2, d2s, d2e
    integer :: lbound2_, range2_(2)
    integer :: i3, d3s, d3e
    integer :: lbound3_, range3_(2)
    !
    lbound1_=1
    if (present(lbound1)) lbound1_=lbound1 
    range1_=(/1,size(array_out, 1)/)
    if (present(range1)) range1_=range1 
    !
    d1s = range1_(1) -lbound1_ +1
    d1e = range1_(2) -lbound1_ +1
    !
    lbound2_=1
    if (present(lbound2)) lbound2_=lbound2 
    range2_=(/1,size(array_out, 2)/)
    if (present(range2)) range2_=range2 
    !
    d2s = range2_(1) -lbound2_ +1
    d2e = range2_(2) -lbound2_ +1
    !
    lbound3_=1
    if (present(lbound3)) lbound3_=lbound3 
    range3_=(/1,size(array_out, 3)/)
    if (present(range3)) range3_=range3 
    !
    d3s = range3_(1) -lbound3_ +1
    d3e = range3_(2) -lbound3_ +1
    !
    !
    !$cuf kernel do(3)
    do i3 = d3s, d3e
    do i2 = d2s, d2e
    do i1 = d1s, d1e
        array_out(i1,i2,i3 ) = val
    enddo
    enddo
    enddo
    !
end subroutine dp_dev_memset_c3d
!
subroutine dp_dev_memset_c4d(array_out, val, & 
                                             
                                            range1, lbound1, & 
                                            range2, lbound2, & 
                                            range3, lbound3, & 
                                            range4, lbound4 )
    implicit none
    !
    integer, parameter :: PRCSN = selected_real_kind(14,200)
    complex(PRCSN), intent(inout) :: array_out(:,:,:,:)
    complex(PRCSN), intent(in)    :: val
    integer, optional, intent(in) ::  range1(2), range2(2), range3(2), range4(2)
    integer, optional, intent(in) ::  lbound1, lbound2, lbound3, lbound4
#if defined(__CUDA)
    attributes(device) :: array_out
#endif
    !
    integer :: i1, d1s, d1e
    integer :: lbound1_, range1_(2)
    integer :: i2, d2s, d2e
    integer :: lbound2_, range2_(2)
    integer :: i3, d3s, d3e
    integer :: lbound3_, range3_(2)
    integer :: i4, d4s, d4e
    integer :: lbound4_, range4_(2)
    !
    lbound1_=1
    if (present(lbound1)) lbound1_=lbound1 
    range1_=(/1,size(array_out, 1)/)
    if (present(range1)) range1_=range1 
    !
    d1s = range1_(1) -lbound1_ +1
    d1e = range1_(2) -lbound1_ +1
    !
    lbound2_=1
    if (present(lbound2)) lbound2_=lbound2 
    range2_=(/1,size(array_out, 2)/)
    if (present(range2)) range2_=range2 
    !
    d2s = range2_(1) -lbound2_ +1
    d2e = range2_(2) -lbound2_ +1
    !
    lbound3_=1
    if (present(lbound3)) lbound3_=lbound3 
    range3_=(/1,size(array_out, 3)/)
    if (present(range3)) range3_=range3 
    !
    d3s = range3_(1) -lbound3_ +1
    d3e = range3_(2) -lbound3_ +1
    !
    lbound4_=1
    if (present(lbound4)) lbound4_=lbound4 
    range4_=(/1,size(array_out, 4)/)
    if (present(range4)) range4_=range4 
    !
    d4s = range4_(1) -lbound4_ +1
    d4e = range4_(2) -lbound4_ +1
    !
    !
    !$cuf kernel do(4)
    do i4 = d4s, d4e
    do i3 = d3s, d3e
    do i2 = d2s, d2e
    do i1 = d1s, d1e
        array_out(i1,i2,i3,i4 ) = val
    enddo
    enddo
    enddo
    enddo
    !
end subroutine dp_dev_memset_c4d
!
subroutine sp_dev_memset_r1d(array_out, val, & 
                                             
                                            range1, lbound1 )
    implicit none
    !
    integer, parameter :: PRCSN = selected_real_kind(6, 37)
    real(PRCSN), intent(inout) :: array_out(:)
    real(PRCSN), intent(in)    :: val
    integer, optional, intent(in) ::  range1(2)
    integer, optional, intent(in) ::  lbound1
#if defined(__CUDA)
    attributes(device) :: array_out
#endif
    !
    integer :: i1, d1s, d1e
    integer :: lbound1_, range1_(2)
    !
    lbound1_=1
    if (present(lbound1)) lbound1_=lbound1 
    range1_=(/1,size(array_out, 1)/)
    if (present(range1)) range1_=range1 
    !
    d1s = range1_(1) -lbound1_ +1
    d1e = range1_(2) -lbound1_ +1
    !
    !
    !$cuf kernel do(1)
    do i1 = d1s, d1e
        array_out(i1 ) = val
    enddo
    !
end subroutine sp_dev_memset_r1d
!
subroutine sp_dev_memset_r2d(array_out, val, & 
                                             
                                            range1, lbound1, & 
                                            range2, lbound2 )
    implicit none
    !
    integer, parameter :: PRCSN = selected_real_kind(6, 37)
    real(PRCSN), intent(inout) :: array_out(:,:)
    real(PRCSN), intent(in)    :: val
    integer, optional, intent(in) ::  range1(2), range2(2)
    integer, optional, intent(in) ::  lbound1, lbound2
#if defined(__CUDA)
    attributes(device) :: array_out
#endif
    !
    integer :: i1, d1s, d1e
    integer :: lbound1_, range1_(2)
    integer :: i2, d2s, d2e
    integer :: lbound2_, range2_(2)
    !
    lbound1_=1
    if (present(lbound1)) lbound1_=lbound1 
    range1_=(/1,size(array_out, 1)/)
    if (present(range1)) range1_=range1 
    !
    d1s = range1_(1) -lbound1_ +1
    d1e = range1_(2) -lbound1_ +1
    !
    lbound2_=1
    if (present(lbound2)) lbound2_=lbound2 
    range2_=(/1,size(array_out, 2)/)
    if (present(range2)) range2_=range2 
    !
    d2s = range2_(1) -lbound2_ +1
    d2e = range2_(2) -lbound2_ +1
    !
    !
    !$cuf kernel do(2)
    do i2 = d2s, d2e
    do i1 = d1s, d1e
        array_out(i1,i2 ) = val
    enddo
    enddo
    !
end subroutine sp_dev_memset_r2d
!
subroutine sp_dev_memset_r3d(array_out, val, & 
                                             
                                            range1, lbound1, & 
                                            range2, lbound2, & 
                                            range3, lbound3 )
    implicit none
    !
    integer, parameter :: PRCSN = selected_real_kind(6, 37)
    real(PRCSN), intent(inout) :: array_out(:,:,:)
    real(PRCSN), intent(in)    :: val
    integer, optional, intent(in) ::  range1(2), range2(2), range3(2)
    integer, optional, intent(in) ::  lbound1, lbound2, lbound3
#if defined(__CUDA)
    attributes(device) :: array_out
#endif
    !
    integer :: i1, d1s, d1e
    integer :: lbound1_, range1_(2)
    integer :: i2, d2s, d2e
    integer :: lbound2_, range2_(2)
    integer :: i3, d3s, d3e
    integer :: lbound3_, range3_(2)
    !
    lbound1_=1
    if (present(lbound1)) lbound1_=lbound1 
    range1_=(/1,size(array_out, 1)/)
    if (present(range1)) range1_=range1 
    !
    d1s = range1_(1) -lbound1_ +1
    d1e = range1_(2) -lbound1_ +1
    !
    lbound2_=1
    if (present(lbound2)) lbound2_=lbound2 
    range2_=(/1,size(array_out, 2)/)
    if (present(range2)) range2_=range2 
    !
    d2s = range2_(1) -lbound2_ +1
    d2e = range2_(2) -lbound2_ +1
    !
    lbound3_=1
    if (present(lbound3)) lbound3_=lbound3 
    range3_=(/1,size(array_out, 3)/)
    if (present(range3)) range3_=range3 
    !
    d3s = range3_(1) -lbound3_ +1
    d3e = range3_(2) -lbound3_ +1
    !
    !
    !$cuf kernel do(3)
    do i3 = d3s, d3e
    do i2 = d2s, d2e
    do i1 = d1s, d1e
        array_out(i1,i2,i3 ) = val
    enddo
    enddo
    enddo
    !
end subroutine sp_dev_memset_r3d
!
subroutine sp_dev_memset_r4d(array_out, val, & 
                                             
                                            range1, lbound1, & 
                                            range2, lbound2, & 
                                            range3, lbound3, & 
                                            range4, lbound4 )
    implicit none
    !
    integer, parameter :: PRCSN = selected_real_kind(6, 37)
    real(PRCSN), intent(inout) :: array_out(:,:,:,:)
    real(PRCSN), intent(in)    :: val
    integer, optional, intent(in) ::  range1(2), range2(2), range3(2), range4(2)
    integer, optional, intent(in) ::  lbound1, lbound2, lbound3, lbound4
#if defined(__CUDA)
    attributes(device) :: array_out
#endif
    !
    integer :: i1, d1s, d1e
    integer :: lbound1_, range1_(2)
    integer :: i2, d2s, d2e
    integer :: lbound2_, range2_(2)
    integer :: i3, d3s, d3e
    integer :: lbound3_, range3_(2)
    integer :: i4, d4s, d4e
    integer :: lbound4_, range4_(2)
    !
    lbound1_=1
    if (present(lbound1)) lbound1_=lbound1 
    range1_=(/1,size(array_out, 1)/)
    if (present(range1)) range1_=range1 
    !
    d1s = range1_(1) -lbound1_ +1
    d1e = range1_(2) -lbound1_ +1
    !
    lbound2_=1
    if (present(lbound2)) lbound2_=lbound2 
    range2_=(/1,size(array_out, 2)/)
    if (present(range2)) range2_=range2 
    !
    d2s = range2_(1) -lbound2_ +1
    d2e = range2_(2) -lbound2_ +1
    !
    lbound3_=1
    if (present(lbound3)) lbound3_=lbound3 
    range3_=(/1,size(array_out, 3)/)
    if (present(range3)) range3_=range3 
    !
    d3s = range3_(1) -lbound3_ +1
    d3e = range3_(2) -lbound3_ +1
    !
    lbound4_=1
    if (present(lbound4)) lbound4_=lbound4 
    range4_=(/1,size(array_out, 4)/)
    if (present(range4)) range4_=range4 
    !
    d4s = range4_(1) -lbound4_ +1
    d4e = range4_(2) -lbound4_ +1
    !
    !
    !$cuf kernel do(4)
    do i4 = d4s, d4e
    do i3 = d3s, d3e
    do i2 = d2s, d2e
    do i1 = d1s, d1e
        array_out(i1,i2,i3,i4 ) = val
    enddo
    enddo
    enddo
    enddo
    !
end subroutine sp_dev_memset_r4d
!
subroutine sp_dev_memset_c1d(array_out, val, & 
                                             
                                            range1, lbound1 )
    implicit none
    !
    integer, parameter :: PRCSN = selected_real_kind(6, 37)
    complex(PRCSN), intent(inout) :: array_out(:)
    complex(PRCSN), intent(in)    :: val
    integer, optional, intent(in) ::  range1(2)
    integer, optional, intent(in) ::  lbound1
#if defined(__CUDA)
    attributes(device) :: array_out
#endif
    !
    integer :: i1, d1s, d1e
    integer :: lbound1_, range1_(2)
    !
    lbound1_=1
    if (present(lbound1)) lbound1_=lbound1 
    range1_=(/1,size(array_out, 1)/)
    if (present(range1)) range1_=range1 
    !
    d1s = range1_(1) -lbound1_ +1
    d1e = range1_(2) -lbound1_ +1
    !
    !
    !$cuf kernel do(1)
    do i1 = d1s, d1e
        array_out(i1 ) = val
    enddo
    !
end subroutine sp_dev_memset_c1d
!
subroutine sp_dev_memset_c2d(array_out, val, & 
                                             
                                            range1, lbound1, & 
                                            range2, lbound2 )
    implicit none
    !
    integer, parameter :: PRCSN = selected_real_kind(6, 37)
    complex(PRCSN), intent(inout) :: array_out(:,:)
    complex(PRCSN), intent(in)    :: val
    integer, optional, intent(in) ::  range1(2), range2(2)
    integer, optional, intent(in) ::  lbound1, lbound2
#if defined(__CUDA)
    attributes(device) :: array_out
#endif
    !
    integer :: i1, d1s, d1e
    integer :: lbound1_, range1_(2)
    integer :: i2, d2s, d2e
    integer :: lbound2_, range2_(2)
    !
    lbound1_=1
    if (present(lbound1)) lbound1_=lbound1 
    range1_=(/1,size(array_out, 1)/)
    if (present(range1)) range1_=range1 
    !
    d1s = range1_(1) -lbound1_ +1
    d1e = range1_(2) -lbound1_ +1
    !
    lbound2_=1
    if (present(lbound2)) lbound2_=lbound2 
    range2_=(/1,size(array_out, 2)/)
    if (present(range2)) range2_=range2 
    !
    d2s = range2_(1) -lbound2_ +1
    d2e = range2_(2) -lbound2_ +1
    !
    !
    !$cuf kernel do(2)
    do i2 = d2s, d2e
    do i1 = d1s, d1e
        array_out(i1,i2 ) = val
    enddo
    enddo
    !
end subroutine sp_dev_memset_c2d
!
subroutine sp_dev_memset_c3d(array_out, val, & 
                                             
                                            range1, lbound1, & 
                                            range2, lbound2, & 
                                            range3, lbound3 )
    implicit none
    !
    integer, parameter :: PRCSN = selected_real_kind(6, 37)
    complex(PRCSN), intent(inout) :: array_out(:,:,:)
    complex(PRCSN), intent(in)    :: val
    integer, optional, intent(in) ::  range1(2), range2(2), range3(2)
    integer, optional, intent(in) ::  lbound1, lbound2, lbound3
#if defined(__CUDA)
    attributes(device) :: array_out
#endif
    !
    integer :: i1, d1s, d1e
    integer :: lbound1_, range1_(2)
    integer :: i2, d2s, d2e
    integer :: lbound2_, range2_(2)
    integer :: i3, d3s, d3e
    integer :: lbound3_, range3_(2)
    !
    lbound1_=1
    if (present(lbound1)) lbound1_=lbound1 
    range1_=(/1,size(array_out, 1)/)
    if (present(range1)) range1_=range1 
    !
    d1s = range1_(1) -lbound1_ +1
    d1e = range1_(2) -lbound1_ +1
    !
    lbound2_=1
    if (present(lbound2)) lbound2_=lbound2 
    range2_=(/1,size(array_out, 2)/)
    if (present(range2)) range2_=range2 
    !
    d2s = range2_(1) -lbound2_ +1
    d2e = range2_(2) -lbound2_ +1
    !
    lbound3_=1
    if (present(lbound3)) lbound3_=lbound3 
    range3_=(/1,size(array_out, 3)/)
    if (present(range3)) range3_=range3 
    !
    d3s = range3_(1) -lbound3_ +1
    d3e = range3_(2) -lbound3_ +1
    !
    !
    !$cuf kernel do(3)
    do i3 = d3s, d3e
    do i2 = d2s, d2e
    do i1 = d1s, d1e
        array_out(i1,i2,i3 ) = val
    enddo
    enddo
    enddo
    !
end subroutine sp_dev_memset_c3d
!
subroutine sp_dev_memset_c4d(array_out, val, & 
                                             
                                            range1, lbound1, & 
                                            range2, lbound2, & 
                                            range3, lbound3, & 
                                            range4, lbound4 )
    implicit none
    !
    integer, parameter :: PRCSN = selected_real_kind(6, 37)
    complex(PRCSN), intent(inout) :: array_out(:,:,:,:)
    complex(PRCSN), intent(in)    :: val
    integer, optional, intent(in) ::  range1(2), range2(2), range3(2), range4(2)
    integer, optional, intent(in) ::  lbound1, lbound2, lbound3, lbound4
#if defined(__CUDA)
    attributes(device) :: array_out
#endif
    !
    integer :: i1, d1s, d1e
    integer :: lbound1_, range1_(2)
    integer :: i2, d2s, d2e
    integer :: lbound2_, range2_(2)
    integer :: i3, d3s, d3e
    integer :: lbound3_, range3_(2)
    integer :: i4, d4s, d4e
    integer :: lbound4_, range4_(2)
    !
    lbound1_=1
    if (present(lbound1)) lbound1_=lbound1 
    range1_=(/1,size(array_out, 1)/)
    if (present(range1)) range1_=range1 
    !
    d1s = range1_(1) -lbound1_ +1
    d1e = range1_(2) -lbound1_ +1
    !
    lbound2_=1
    if (present(lbound2)) lbound2_=lbound2 
    range2_=(/1,size(array_out, 2)/)
    if (present(range2)) range2_=range2 
    !
    d2s = range2_(1) -lbound2_ +1
    d2e = range2_(2) -lbound2_ +1
    !
    lbound3_=1
    if (present(lbound3)) lbound3_=lbound3 
    range3_=(/1,size(array_out, 3)/)
    if (present(range3)) range3_=range3 
    !
    d3s = range3_(1) -lbound3_ +1
    d3e = range3_(2) -lbound3_ +1
    !
    lbound4_=1
    if (present(lbound4)) lbound4_=lbound4 
    range4_=(/1,size(array_out, 4)/)
    if (present(range4)) range4_=range4 
    !
    d4s = range4_(1) -lbound4_ +1
    d4e = range4_(2) -lbound4_ +1
    !
    !
    !$cuf kernel do(4)
    do i4 = d4s, d4e
    do i3 = d3s, d3e
    do i2 = d2s, d2e
    do i1 = d1s, d1e
        array_out(i1,i2,i3,i4 ) = val
    enddo
    enddo
    enddo
    enddo
    !
end subroutine sp_dev_memset_c4d
!
!
!======================
!
subroutine dp_h2d_memsync_r1d(array_out, array_in, &
                                             range1, lbound1  )
#if defined(__CUDA)
    use cudafor
#endif
    implicit none
    !
    integer, parameter :: PRCSN = selected_real_kind(14,200)
    real(PRCSN), intent(inout) :: array_out(:)
    real(PRCSN), intent(in)    :: array_in(:)
    integer, optional, intent(in) ::  range1(2)
    integer, optional, intent(in) ::  lbound1
#if defined(__CUDA)
    attributes(device) :: array_out
    integer :: ierr
#endif
    !
    integer :: d1_start, d1_end, d1_size, d1_ld
    integer :: lbound1_, range1_(2)
    !
    lbound1_=1
    if (present(lbound1)) lbound1_=lbound1 
    range1_=(/1,size(array_out, 1)/)
    if (present(range1)) range1_=range1 
    !
    d1_start = range1_(1) -lbound1_ +1
    d1_end   = range1_(2) -lbound1_ +1
    d1_size  = range1_(2) -range1_(1) + 1
    d1_ld    = size(array_out, 1)
    !
#if defined(__CUDA)
    ierr = cudaMemcpy( array_out(d1_start), array_in(d1_start), d1_size, cudaMemcpyHostToDevice )
#else
    array_out(d1_start:d1_end) = &
              array_in(d1_start:d1_end)
#endif
    !
end subroutine dp_h2d_memsync_r1d
!
subroutine dp_h2d_memsync_r2d(array_out, array_in, &
                                             range1, lbound1 , &
                                             range2, lbound2  )
#if defined(__CUDA)
    use cudafor
#endif
    implicit none
    !
    integer, parameter :: PRCSN = selected_real_kind(14,200)
    real(PRCSN), intent(inout) :: array_out(:,:)
    real(PRCSN), intent(in)    :: array_in(:,:)
    integer, optional, intent(in) ::  range1(2), range2(2)
    integer, optional, intent(in) ::  lbound1, lbound2
#if defined(__CUDA)
    attributes(device) :: array_out
    integer :: ierr
#endif
    !
    integer :: d1_start, d1_end, d1_size, d1_ld
    integer :: lbound1_, range1_(2)
    integer :: d2_start, d2_end, d2_size, d2_ld
    integer :: lbound2_, range2_(2)
    !
    lbound1_=1
    if (present(lbound1)) lbound1_=lbound1 
    range1_=(/1,size(array_out, 1)/)
    if (present(range1)) range1_=range1 
    !
    d1_start = range1_(1) -lbound1_ +1
    d1_end   = range1_(2) -lbound1_ +1
    d1_size  = range1_(2) -range1_(1) + 1
    d1_ld    = size(array_out, 1)
    lbound2_=1
    if (present(lbound2)) lbound2_=lbound2 
    range2_=(/1,size(array_out, 2)/)
    if (present(range2)) range2_=range2 
    !
    d2_start = range2_(1) -lbound2_ +1
    d2_end   = range2_(2) -lbound2_ +1
    d2_size  = range2_(2) -range2_(1) + 1
    d2_ld    = size(array_out, 2)
    !
#if defined(__CUDA)
    ierr = cudaMemcpy2D( array_out(d1_start, d2_start) , d1_ld, array_in(d1_start, d2_start), d2_ld, d1_size, d2_size )
#else
    array_out(d1_start:d1_end,d2_start:d2_end) = &
              array_in(d1_start:d1_end,d2_start:d2_end)
#endif
    !
end subroutine dp_h2d_memsync_r2d
!
subroutine dp_h2d_memsync_r3d(array_out, array_in, &
                                             range1, lbound1 , &
                                             range2, lbound2 , &
                                             range3, lbound3  )
#if defined(__CUDA)
    use cudafor
#endif
    implicit none
    !
    integer, parameter :: PRCSN = selected_real_kind(14,200)
    real(PRCSN), intent(inout) :: array_out(:,:,:)
    real(PRCSN), intent(in)    :: array_in(:,:,:)
    integer, optional, intent(in) ::  range1(2), range2(2), range3(2)
    integer, optional, intent(in) ::  lbound1, lbound2, lbound3
#if defined(__CUDA)
    attributes(device) :: array_out
    integer :: ierr
#endif
    !
    integer :: d1_start, d1_end, d1_size, d1_ld
    integer :: lbound1_, range1_(2)
    integer :: d2_start, d2_end, d2_size, d2_ld
    integer :: lbound2_, range2_(2)
    integer :: d3_start, d3_end, d3_size, d3_ld
    integer :: lbound3_, range3_(2)
    !
    lbound1_=1
    if (present(lbound1)) lbound1_=lbound1 
    range1_=(/1,size(array_out, 1)/)
    if (present(range1)) range1_=range1 
    !
    d1_start = range1_(1) -lbound1_ +1
    d1_end   = range1_(2) -lbound1_ +1
    d1_size  = range1_(2) -range1_(1) + 1
    d1_ld    = size(array_out, 1)
    lbound2_=1
    if (present(lbound2)) lbound2_=lbound2 
    range2_=(/1,size(array_out, 2)/)
    if (present(range2)) range2_=range2 
    !
    d2_start = range2_(1) -lbound2_ +1
    d2_end   = range2_(2) -lbound2_ +1
    d2_size  = range2_(2) -range2_(1) + 1
    d2_ld    = size(array_out, 2)
    lbound3_=1
    if (present(lbound3)) lbound3_=lbound3 
    range3_=(/1,size(array_out, 3)/)
    if (present(range3)) range3_=range3 
    !
    d3_start = range3_(1) -lbound3_ +1
    d3_end   = range3_(2) -lbound3_ +1
    d3_size  = range3_(2) -range3_(1) + 1
    d3_ld    = size(array_out, 3)
    !
#if defined(__CUDA)
    !call dxlib_errore('cu_memsync_','3D arrays not implemented yet',1)
    array_out(d1_start:d1_end,d2_start:d2_end,d3_start:d3_end) = &
              array_in(d1_start:d1_end,d2_start:d2_end,d3_start:d3_end)
#else
    array_out(d1_start:d1_end,d2_start:d2_end,d3_start:d3_end) = &
              array_in(d1_start:d1_end,d2_start:d2_end,d3_start:d3_end)
#endif
    !
end subroutine dp_h2d_memsync_r3d
!
subroutine dp_h2d_memsync_r4d(array_out, array_in, &
                                             range1, lbound1 , &
                                             range2, lbound2 , &
                                             range3, lbound3 , &
                                             range4, lbound4  )
#if defined(__CUDA)
    use cudafor
#endif
    implicit none
    !
    integer, parameter :: PRCSN = selected_real_kind(14,200)
    real(PRCSN), intent(inout) :: array_out(:,:,:,:)
    real(PRCSN), intent(in)    :: array_in(:,:,:,:)
    integer, optional, intent(in) ::  range1(2), range2(2), range3(2), range4(2)
    integer, optional, intent(in) ::  lbound1, lbound2, lbound3, lbound4
#if defined(__CUDA)
    attributes(device) :: array_out
    integer :: ierr
#endif
    !
    integer :: d1_start, d1_end, d1_size, d1_ld
    integer :: lbound1_, range1_(2)
    integer :: d2_start, d2_end, d2_size, d2_ld
    integer :: lbound2_, range2_(2)
    integer :: d3_start, d3_end, d3_size, d3_ld
    integer :: lbound3_, range3_(2)
    integer :: d4_start, d4_end, d4_size, d4_ld
    integer :: lbound4_, range4_(2)
    !
    lbound1_=1
    if (present(lbound1)) lbound1_=lbound1 
    range1_=(/1,size(array_out, 1)/)
    if (present(range1)) range1_=range1 
    !
    d1_start = range1_(1) -lbound1_ +1
    d1_end   = range1_(2) -lbound1_ +1
    d1_size  = range1_(2) -range1_(1) + 1
    d1_ld    = size(array_out, 1)
    lbound2_=1
    if (present(lbound2)) lbound2_=lbound2 
    range2_=(/1,size(array_out, 2)/)
    if (present(range2)) range2_=range2 
    !
    d2_start = range2_(1) -lbound2_ +1
    d2_end   = range2_(2) -lbound2_ +1
    d2_size  = range2_(2) -range2_(1) + 1
    d2_ld    = size(array_out, 2)
    lbound3_=1
    if (present(lbound3)) lbound3_=lbound3 
    range3_=(/1,size(array_out, 3)/)
    if (present(range3)) range3_=range3 
    !
    d3_start = range3_(1) -lbound3_ +1
    d3_end   = range3_(2) -lbound3_ +1
    d3_size  = range3_(2) -range3_(1) + 1
    d3_ld    = size(array_out, 3)
    lbound4_=1
    if (present(lbound4)) lbound4_=lbound4 
    range4_=(/1,size(array_out, 4)/)
    if (present(range4)) range4_=range4 
    !
    d4_start = range4_(1) -lbound4_ +1
    d4_end   = range4_(2) -lbound4_ +1
    d4_size  = range4_(2) -range4_(1) + 1
    d4_ld    = size(array_out, 4)
    !
#if defined(__CUDA)
    !call dxlib_errore('cu_memsync_','4D arrays not implemented yet',1)
    array_out(d1_start:d1_end,d2_start:d2_end,d3_start:d3_end,d4_start:d4_end) = &
              array_in(d1_start:d1_end,d2_start:d2_end,d3_start:d3_end,d4_start:d4_end)
#else
    array_out(d1_start:d1_end,d2_start:d2_end,d3_start:d3_end,d4_start:d4_end) = &
              array_in(d1_start:d1_end,d2_start:d2_end,d3_start:d3_end,d4_start:d4_end)
#endif
    !
end subroutine dp_h2d_memsync_r4d
!
subroutine dp_h2d_memsync_c1d(array_out, array_in, &
                                             range1, lbound1  )
#if defined(__CUDA)
    use cudafor
#endif
    implicit none
    !
    integer, parameter :: PRCSN = selected_real_kind(14,200)
    complex(PRCSN), intent(inout) :: array_out(:)
    complex(PRCSN), intent(in)    :: array_in(:)
    integer, optional, intent(in) ::  range1(2)
    integer, optional, intent(in) ::  lbound1
#if defined(__CUDA)
    attributes(device) :: array_out
    integer :: ierr
#endif
    !
    integer :: d1_start, d1_end, d1_size, d1_ld
    integer :: lbound1_, range1_(2)
    !
    lbound1_=1
    if (present(lbound1)) lbound1_=lbound1 
    range1_=(/1,size(array_out, 1)/)
    if (present(range1)) range1_=range1 
    !
    d1_start = range1_(1) -lbound1_ +1
    d1_end   = range1_(2) -lbound1_ +1
    d1_size  = range1_(2) -range1_(1) + 1
    d1_ld    = size(array_out, 1)
    !
#if defined(__CUDA)
    ierr = cudaMemcpy( array_out(d1_start), array_in(d1_start), d1_size, cudaMemcpyHostToDevice )
#else
    array_out(d1_start:d1_end) = &
              array_in(d1_start:d1_end)
#endif
    !
end subroutine dp_h2d_memsync_c1d
!
subroutine dp_h2d_memsync_c2d(array_out, array_in, &
                                             range1, lbound1 , &
                                             range2, lbound2  )
#if defined(__CUDA)
    use cudafor
#endif
    implicit none
    !
    integer, parameter :: PRCSN = selected_real_kind(14,200)
    complex(PRCSN), intent(inout) :: array_out(:,:)
    complex(PRCSN), intent(in)    :: array_in(:,:)
    integer, optional, intent(in) ::  range1(2), range2(2)
    integer, optional, intent(in) ::  lbound1, lbound2
#if defined(__CUDA)
    attributes(device) :: array_out
    integer :: ierr
#endif
    !
    integer :: d1_start, d1_end, d1_size, d1_ld
    integer :: lbound1_, range1_(2)
    integer :: d2_start, d2_end, d2_size, d2_ld
    integer :: lbound2_, range2_(2)
    !
    lbound1_=1
    if (present(lbound1)) lbound1_=lbound1 
    range1_=(/1,size(array_out, 1)/)
    if (present(range1)) range1_=range1 
    !
    d1_start = range1_(1) -lbound1_ +1
    d1_end   = range1_(2) -lbound1_ +1
    d1_size  = range1_(2) -range1_(1) + 1
    d1_ld    = size(array_out, 1)
    lbound2_=1
    if (present(lbound2)) lbound2_=lbound2 
    range2_=(/1,size(array_out, 2)/)
    if (present(range2)) range2_=range2 
    !
    d2_start = range2_(1) -lbound2_ +1
    d2_end   = range2_(2) -lbound2_ +1
    d2_size  = range2_(2) -range2_(1) + 1
    d2_ld    = size(array_out, 2)
    !
#if defined(__CUDA)
    ierr = cudaMemcpy2D( array_out(d1_start, d2_start) , d1_ld, array_in(d1_start, d2_start), d2_ld, d1_size, d2_size )
#else
    array_out(d1_start:d1_end,d2_start:d2_end) = &
              array_in(d1_start:d1_end,d2_start:d2_end)
#endif
    !
end subroutine dp_h2d_memsync_c2d
!
subroutine dp_h2d_memsync_c3d(array_out, array_in, &
                                             range1, lbound1 , &
                                             range2, lbound2 , &
                                             range3, lbound3  )
#if defined(__CUDA)
    use cudafor
#endif
    implicit none
    !
    integer, parameter :: PRCSN = selected_real_kind(14,200)
    complex(PRCSN), intent(inout) :: array_out(:,:,:)
    complex(PRCSN), intent(in)    :: array_in(:,:,:)
    integer, optional, intent(in) ::  range1(2), range2(2), range3(2)
    integer, optional, intent(in) ::  lbound1, lbound2, lbound3
#if defined(__CUDA)
    attributes(device) :: array_out
    integer :: ierr
#endif
    !
    integer :: d1_start, d1_end, d1_size, d1_ld
    integer :: lbound1_, range1_(2)
    integer :: d2_start, d2_end, d2_size, d2_ld
    integer :: lbound2_, range2_(2)
    integer :: d3_start, d3_end, d3_size, d3_ld
    integer :: lbound3_, range3_(2)
    !
    lbound1_=1
    if (present(lbound1)) lbound1_=lbound1 
    range1_=(/1,size(array_out, 1)/)
    if (present(range1)) range1_=range1 
    !
    d1_start = range1_(1) -lbound1_ +1
    d1_end   = range1_(2) -lbound1_ +1
    d1_size  = range1_(2) -range1_(1) + 1
    d1_ld    = size(array_out, 1)
    lbound2_=1
    if (present(lbound2)) lbound2_=lbound2 
    range2_=(/1,size(array_out, 2)/)
    if (present(range2)) range2_=range2 
    !
    d2_start = range2_(1) -lbound2_ +1
    d2_end   = range2_(2) -lbound2_ +1
    d2_size  = range2_(2) -range2_(1) + 1
    d2_ld    = size(array_out, 2)
    lbound3_=1
    if (present(lbound3)) lbound3_=lbound3 
    range3_=(/1,size(array_out, 3)/)
    if (present(range3)) range3_=range3 
    !
    d3_start = range3_(1) -lbound3_ +1
    d3_end   = range3_(2) -lbound3_ +1
    d3_size  = range3_(2) -range3_(1) + 1
    d3_ld    = size(array_out, 3)
    !
#if defined(__CUDA)
    !call dxlib_errore('cu_memsync_','3D arrays not implemented yet',1)
    array_out(d1_start:d1_end,d2_start:d2_end,d3_start:d3_end) = &
              array_in(d1_start:d1_end,d2_start:d2_end,d3_start:d3_end)
#else
    array_out(d1_start:d1_end,d2_start:d2_end,d3_start:d3_end) = &
              array_in(d1_start:d1_end,d2_start:d2_end,d3_start:d3_end)
#endif
    !
end subroutine dp_h2d_memsync_c3d
!
subroutine dp_h2d_memsync_c4d(array_out, array_in, &
                                             range1, lbound1 , &
                                             range2, lbound2 , &
                                             range3, lbound3 , &
                                             range4, lbound4  )
#if defined(__CUDA)
    use cudafor
#endif
    implicit none
    !
    integer, parameter :: PRCSN = selected_real_kind(14,200)
    complex(PRCSN), intent(inout) :: array_out(:,:,:,:)
    complex(PRCSN), intent(in)    :: array_in(:,:,:,:)
    integer, optional, intent(in) ::  range1(2), range2(2), range3(2), range4(2)
    integer, optional, intent(in) ::  lbound1, lbound2, lbound3, lbound4
#if defined(__CUDA)
    attributes(device) :: array_out
    integer :: ierr
#endif
    !
    integer :: d1_start, d1_end, d1_size, d1_ld
    integer :: lbound1_, range1_(2)
    integer :: d2_start, d2_end, d2_size, d2_ld
    integer :: lbound2_, range2_(2)
    integer :: d3_start, d3_end, d3_size, d3_ld
    integer :: lbound3_, range3_(2)
    integer :: d4_start, d4_end, d4_size, d4_ld
    integer :: lbound4_, range4_(2)
    !
    lbound1_=1
    if (present(lbound1)) lbound1_=lbound1 
    range1_=(/1,size(array_out, 1)/)
    if (present(range1)) range1_=range1 
    !
    d1_start = range1_(1) -lbound1_ +1
    d1_end   = range1_(2) -lbound1_ +1
    d1_size  = range1_(2) -range1_(1) + 1
    d1_ld    = size(array_out, 1)
    lbound2_=1
    if (present(lbound2)) lbound2_=lbound2 
    range2_=(/1,size(array_out, 2)/)
    if (present(range2)) range2_=range2 
    !
    d2_start = range2_(1) -lbound2_ +1
    d2_end   = range2_(2) -lbound2_ +1
    d2_size  = range2_(2) -range2_(1) + 1
    d2_ld    = size(array_out, 2)
    lbound3_=1
    if (present(lbound3)) lbound3_=lbound3 
    range3_=(/1,size(array_out, 3)/)
    if (present(range3)) range3_=range3 
    !
    d3_start = range3_(1) -lbound3_ +1
    d3_end   = range3_(2) -lbound3_ +1
    d3_size  = range3_(2) -range3_(1) + 1
    d3_ld    = size(array_out, 3)
    lbound4_=1
    if (present(lbound4)) lbound4_=lbound4 
    range4_=(/1,size(array_out, 4)/)
    if (present(range4)) range4_=range4 
    !
    d4_start = range4_(1) -lbound4_ +1
    d4_end   = range4_(2) -lbound4_ +1
    d4_size  = range4_(2) -range4_(1) + 1
    d4_ld    = size(array_out, 4)
    !
#if defined(__CUDA)
    !call dxlib_errore('cu_memsync_','4D arrays not implemented yet',1)
    array_out(d1_start:d1_end,d2_start:d2_end,d3_start:d3_end,d4_start:d4_end) = &
              array_in(d1_start:d1_end,d2_start:d2_end,d3_start:d3_end,d4_start:d4_end)
#else
    array_out(d1_start:d1_end,d2_start:d2_end,d3_start:d3_end,d4_start:d4_end) = &
              array_in(d1_start:d1_end,d2_start:d2_end,d3_start:d3_end,d4_start:d4_end)
#endif
    !
end subroutine dp_h2d_memsync_c4d
!
subroutine sp_h2d_memsync_r1d(array_out, array_in, &
                                             range1, lbound1  )
#if defined(__CUDA)
    use cudafor
#endif
    implicit none
    !
    integer, parameter :: PRCSN = selected_real_kind(6, 37)
    real(PRCSN), intent(inout) :: array_out(:)
    real(PRCSN), intent(in)    :: array_in(:)
    integer, optional, intent(in) ::  range1(2)
    integer, optional, intent(in) ::  lbound1
#if defined(__CUDA)
    attributes(device) :: array_out
    integer :: ierr
#endif
    !
    integer :: d1_start, d1_end, d1_size, d1_ld
    integer :: lbound1_, range1_(2)
    !
    lbound1_=1
    if (present(lbound1)) lbound1_=lbound1 
    range1_=(/1,size(array_out, 1)/)
    if (present(range1)) range1_=range1 
    !
    d1_start = range1_(1) -lbound1_ +1
    d1_end   = range1_(2) -lbound1_ +1
    d1_size  = range1_(2) -range1_(1) + 1
    d1_ld    = size(array_out, 1)
    !
#if defined(__CUDA)
    ierr = cudaMemcpy( array_out(d1_start), array_in(d1_start), d1_size, cudaMemcpyHostToDevice )
#else
    array_out(d1_start:d1_end) = &
              array_in(d1_start:d1_end)
#endif
    !
end subroutine sp_h2d_memsync_r1d
!
subroutine sp_h2d_memsync_r2d(array_out, array_in, &
                                             range1, lbound1 , &
                                             range2, lbound2  )
#if defined(__CUDA)
    use cudafor
#endif
    implicit none
    !
    integer, parameter :: PRCSN = selected_real_kind(6, 37)
    real(PRCSN), intent(inout) :: array_out(:,:)
    real(PRCSN), intent(in)    :: array_in(:,:)
    integer, optional, intent(in) ::  range1(2), range2(2)
    integer, optional, intent(in) ::  lbound1, lbound2
#if defined(__CUDA)
    attributes(device) :: array_out
    integer :: ierr
#endif
    !
    integer :: d1_start, d1_end, d1_size, d1_ld
    integer :: lbound1_, range1_(2)
    integer :: d2_start, d2_end, d2_size, d2_ld
    integer :: lbound2_, range2_(2)
    !
    lbound1_=1
    if (present(lbound1)) lbound1_=lbound1 
    range1_=(/1,size(array_out, 1)/)
    if (present(range1)) range1_=range1 
    !
    d1_start = range1_(1) -lbound1_ +1
    d1_end   = range1_(2) -lbound1_ +1
    d1_size  = range1_(2) -range1_(1) + 1
    d1_ld    = size(array_out, 1)
    lbound2_=1
    if (present(lbound2)) lbound2_=lbound2 
    range2_=(/1,size(array_out, 2)/)
    if (present(range2)) range2_=range2 
    !
    d2_start = range2_(1) -lbound2_ +1
    d2_end   = range2_(2) -lbound2_ +1
    d2_size  = range2_(2) -range2_(1) + 1
    d2_ld    = size(array_out, 2)
    !
#if defined(__CUDA)
    ierr = cudaMemcpy2D( array_out(d1_start, d2_start) , d1_ld, array_in(d1_start, d2_start), d2_ld, d1_size, d2_size )
#else
    array_out(d1_start:d1_end,d2_start:d2_end) = &
              array_in(d1_start:d1_end,d2_start:d2_end)
#endif
    !
end subroutine sp_h2d_memsync_r2d
!
subroutine sp_h2d_memsync_r3d(array_out, array_in, &
                                             range1, lbound1 , &
                                             range2, lbound2 , &
                                             range3, lbound3  )
#if defined(__CUDA)
    use cudafor
#endif
    implicit none
    !
    integer, parameter :: PRCSN = selected_real_kind(6, 37)
    real(PRCSN), intent(inout) :: array_out(:,:,:)
    real(PRCSN), intent(in)    :: array_in(:,:,:)
    integer, optional, intent(in) ::  range1(2), range2(2), range3(2)
    integer, optional, intent(in) ::  lbound1, lbound2, lbound3
#if defined(__CUDA)
    attributes(device) :: array_out
    integer :: ierr
#endif
    !
    integer :: d1_start, d1_end, d1_size, d1_ld
    integer :: lbound1_, range1_(2)
    integer :: d2_start, d2_end, d2_size, d2_ld
    integer :: lbound2_, range2_(2)
    integer :: d3_start, d3_end, d3_size, d3_ld
    integer :: lbound3_, range3_(2)
    !
    lbound1_=1
    if (present(lbound1)) lbound1_=lbound1 
    range1_=(/1,size(array_out, 1)/)
    if (present(range1)) range1_=range1 
    !
    d1_start = range1_(1) -lbound1_ +1
    d1_end   = range1_(2) -lbound1_ +1
    d1_size  = range1_(2) -range1_(1) + 1
    d1_ld    = size(array_out, 1)
    lbound2_=1
    if (present(lbound2)) lbound2_=lbound2 
    range2_=(/1,size(array_out, 2)/)
    if (present(range2)) range2_=range2 
    !
    d2_start = range2_(1) -lbound2_ +1
    d2_end   = range2_(2) -lbound2_ +1
    d2_size  = range2_(2) -range2_(1) + 1
    d2_ld    = size(array_out, 2)
    lbound3_=1
    if (present(lbound3)) lbound3_=lbound3 
    range3_=(/1,size(array_out, 3)/)
    if (present(range3)) range3_=range3 
    !
    d3_start = range3_(1) -lbound3_ +1
    d3_end   = range3_(2) -lbound3_ +1
    d3_size  = range3_(2) -range3_(1) + 1
    d3_ld    = size(array_out, 3)
    !
#if defined(__CUDA)
    !call dxlib_errore('cu_memsync_','3D arrays not implemented yet',1)
    array_out(d1_start:d1_end,d2_start:d2_end,d3_start:d3_end) = &
              array_in(d1_start:d1_end,d2_start:d2_end,d3_start:d3_end)
#else
    array_out(d1_start:d1_end,d2_start:d2_end,d3_start:d3_end) = &
              array_in(d1_start:d1_end,d2_start:d2_end,d3_start:d3_end)
#endif
    !
end subroutine sp_h2d_memsync_r3d
!
subroutine sp_h2d_memsync_r4d(array_out, array_in, &
                                             range1, lbound1 , &
                                             range2, lbound2 , &
                                             range3, lbound3 , &
                                             range4, lbound4  )
#if defined(__CUDA)
    use cudafor
#endif
    implicit none
    !
    integer, parameter :: PRCSN = selected_real_kind(6, 37)
    real(PRCSN), intent(inout) :: array_out(:,:,:,:)
    real(PRCSN), intent(in)    :: array_in(:,:,:,:)
    integer, optional, intent(in) ::  range1(2), range2(2), range3(2), range4(2)
    integer, optional, intent(in) ::  lbound1, lbound2, lbound3, lbound4
#if defined(__CUDA)
    attributes(device) :: array_out
    integer :: ierr
#endif
    !
    integer :: d1_start, d1_end, d1_size, d1_ld
    integer :: lbound1_, range1_(2)
    integer :: d2_start, d2_end, d2_size, d2_ld
    integer :: lbound2_, range2_(2)
    integer :: d3_start, d3_end, d3_size, d3_ld
    integer :: lbound3_, range3_(2)
    integer :: d4_start, d4_end, d4_size, d4_ld
    integer :: lbound4_, range4_(2)
    !
    lbound1_=1
    if (present(lbound1)) lbound1_=lbound1 
    range1_=(/1,size(array_out, 1)/)
    if (present(range1)) range1_=range1 
    !
    d1_start = range1_(1) -lbound1_ +1
    d1_end   = range1_(2) -lbound1_ +1
    d1_size  = range1_(2) -range1_(1) + 1
    d1_ld    = size(array_out, 1)
    lbound2_=1
    if (present(lbound2)) lbound2_=lbound2 
    range2_=(/1,size(array_out, 2)/)
    if (present(range2)) range2_=range2 
    !
    d2_start = range2_(1) -lbound2_ +1
    d2_end   = range2_(2) -lbound2_ +1
    d2_size  = range2_(2) -range2_(1) + 1
    d2_ld    = size(array_out, 2)
    lbound3_=1
    if (present(lbound3)) lbound3_=lbound3 
    range3_=(/1,size(array_out, 3)/)
    if (present(range3)) range3_=range3 
    !
    d3_start = range3_(1) -lbound3_ +1
    d3_end   = range3_(2) -lbound3_ +1
    d3_size  = range3_(2) -range3_(1) + 1
    d3_ld    = size(array_out, 3)
    lbound4_=1
    if (present(lbound4)) lbound4_=lbound4 
    range4_=(/1,size(array_out, 4)/)
    if (present(range4)) range4_=range4 
    !
    d4_start = range4_(1) -lbound4_ +1
    d4_end   = range4_(2) -lbound4_ +1
    d4_size  = range4_(2) -range4_(1) + 1
    d4_ld    = size(array_out, 4)
    !
#if defined(__CUDA)
    !call dxlib_errore('cu_memsync_','4D arrays not implemented yet',1)
    array_out(d1_start:d1_end,d2_start:d2_end,d3_start:d3_end,d4_start:d4_end) = &
              array_in(d1_start:d1_end,d2_start:d2_end,d3_start:d3_end,d4_start:d4_end)
#else
    array_out(d1_start:d1_end,d2_start:d2_end,d3_start:d3_end,d4_start:d4_end) = &
              array_in(d1_start:d1_end,d2_start:d2_end,d3_start:d3_end,d4_start:d4_end)
#endif
    !
end subroutine sp_h2d_memsync_r4d
!
subroutine sp_h2d_memsync_c1d(array_out, array_in, &
                                             range1, lbound1  )
#if defined(__CUDA)
    use cudafor
#endif
    implicit none
    !
    integer, parameter :: PRCSN = selected_real_kind(6, 37)
    complex(PRCSN), intent(inout) :: array_out(:)
    complex(PRCSN), intent(in)    :: array_in(:)
    integer, optional, intent(in) ::  range1(2)
    integer, optional, intent(in) ::  lbound1
#if defined(__CUDA)
    attributes(device) :: array_out
    integer :: ierr
#endif
    !
    integer :: d1_start, d1_end, d1_size, d1_ld
    integer :: lbound1_, range1_(2)
    !
    lbound1_=1
    if (present(lbound1)) lbound1_=lbound1 
    range1_=(/1,size(array_out, 1)/)
    if (present(range1)) range1_=range1 
    !
    d1_start = range1_(1) -lbound1_ +1
    d1_end   = range1_(2) -lbound1_ +1
    d1_size  = range1_(2) -range1_(1) + 1
    d1_ld    = size(array_out, 1)
    !
#if defined(__CUDA)
    ierr = cudaMemcpy( array_out(d1_start), array_in(d1_start), d1_size, cudaMemcpyHostToDevice )
#else
    array_out(d1_start:d1_end) = &
              array_in(d1_start:d1_end)
#endif
    !
end subroutine sp_h2d_memsync_c1d
!
subroutine sp_h2d_memsync_c2d(array_out, array_in, &
                                             range1, lbound1 , &
                                             range2, lbound2  )
#if defined(__CUDA)
    use cudafor
#endif
    implicit none
    !
    integer, parameter :: PRCSN = selected_real_kind(6, 37)
    complex(PRCSN), intent(inout) :: array_out(:,:)
    complex(PRCSN), intent(in)    :: array_in(:,:)
    integer, optional, intent(in) ::  range1(2), range2(2)
    integer, optional, intent(in) ::  lbound1, lbound2
#if defined(__CUDA)
    attributes(device) :: array_out
    integer :: ierr
#endif
    !
    integer :: d1_start, d1_end, d1_size, d1_ld
    integer :: lbound1_, range1_(2)
    integer :: d2_start, d2_end, d2_size, d2_ld
    integer :: lbound2_, range2_(2)
    !
    lbound1_=1
    if (present(lbound1)) lbound1_=lbound1 
    range1_=(/1,size(array_out, 1)/)
    if (present(range1)) range1_=range1 
    !
    d1_start = range1_(1) -lbound1_ +1
    d1_end   = range1_(2) -lbound1_ +1
    d1_size  = range1_(2) -range1_(1) + 1
    d1_ld    = size(array_out, 1)
    lbound2_=1
    if (present(lbound2)) lbound2_=lbound2 
    range2_=(/1,size(array_out, 2)/)
    if (present(range2)) range2_=range2 
    !
    d2_start = range2_(1) -lbound2_ +1
    d2_end   = range2_(2) -lbound2_ +1
    d2_size  = range2_(2) -range2_(1) + 1
    d2_ld    = size(array_out, 2)
    !
#if defined(__CUDA)
    ierr = cudaMemcpy2D( array_out(d1_start, d2_start) , d1_ld, array_in(d1_start, d2_start), d2_ld, d1_size, d2_size )
#else
    array_out(d1_start:d1_end,d2_start:d2_end) = &
              array_in(d1_start:d1_end,d2_start:d2_end)
#endif
    !
end subroutine sp_h2d_memsync_c2d
!
subroutine sp_h2d_memsync_c3d(array_out, array_in, &
                                             range1, lbound1 , &
                                             range2, lbound2 , &
                                             range3, lbound3  )
#if defined(__CUDA)
    use cudafor
#endif
    implicit none
    !
    integer, parameter :: PRCSN = selected_real_kind(6, 37)
    complex(PRCSN), intent(inout) :: array_out(:,:,:)
    complex(PRCSN), intent(in)    :: array_in(:,:,:)
    integer, optional, intent(in) ::  range1(2), range2(2), range3(2)
    integer, optional, intent(in) ::  lbound1, lbound2, lbound3
#if defined(__CUDA)
    attributes(device) :: array_out
    integer :: ierr
#endif
    !
    integer :: d1_start, d1_end, d1_size, d1_ld
    integer :: lbound1_, range1_(2)
    integer :: d2_start, d2_end, d2_size, d2_ld
    integer :: lbound2_, range2_(2)
    integer :: d3_start, d3_end, d3_size, d3_ld
    integer :: lbound3_, range3_(2)
    !
    lbound1_=1
    if (present(lbound1)) lbound1_=lbound1 
    range1_=(/1,size(array_out, 1)/)
    if (present(range1)) range1_=range1 
    !
    d1_start = range1_(1) -lbound1_ +1
    d1_end   = range1_(2) -lbound1_ +1
    d1_size  = range1_(2) -range1_(1) + 1
    d1_ld    = size(array_out, 1)
    lbound2_=1
    if (present(lbound2)) lbound2_=lbound2 
    range2_=(/1,size(array_out, 2)/)
    if (present(range2)) range2_=range2 
    !
    d2_start = range2_(1) -lbound2_ +1
    d2_end   = range2_(2) -lbound2_ +1
    d2_size  = range2_(2) -range2_(1) + 1
    d2_ld    = size(array_out, 2)
    lbound3_=1
    if (present(lbound3)) lbound3_=lbound3 
    range3_=(/1,size(array_out, 3)/)
    if (present(range3)) range3_=range3 
    !
    d3_start = range3_(1) -lbound3_ +1
    d3_end   = range3_(2) -lbound3_ +1
    d3_size  = range3_(2) -range3_(1) + 1
    d3_ld    = size(array_out, 3)
    !
#if defined(__CUDA)
    !call dxlib_errore('cu_memsync_','3D arrays not implemented yet',1)
    array_out(d1_start:d1_end,d2_start:d2_end,d3_start:d3_end) = &
              array_in(d1_start:d1_end,d2_start:d2_end,d3_start:d3_end)
#else
    array_out(d1_start:d1_end,d2_start:d2_end,d3_start:d3_end) = &
              array_in(d1_start:d1_end,d2_start:d2_end,d3_start:d3_end)
#endif
    !
end subroutine sp_h2d_memsync_c3d
!
subroutine sp_h2d_memsync_c4d(array_out, array_in, &
                                             range1, lbound1 , &
                                             range2, lbound2 , &
                                             range3, lbound3 , &
                                             range4, lbound4  )
#if defined(__CUDA)
    use cudafor
#endif
    implicit none
    !
    integer, parameter :: PRCSN = selected_real_kind(6, 37)
    complex(PRCSN), intent(inout) :: array_out(:,:,:,:)
    complex(PRCSN), intent(in)    :: array_in(:,:,:,:)
    integer, optional, intent(in) ::  range1(2), range2(2), range3(2), range4(2)
    integer, optional, intent(in) ::  lbound1, lbound2, lbound3, lbound4
#if defined(__CUDA)
    attributes(device) :: array_out
    integer :: ierr
#endif
    !
    integer :: d1_start, d1_end, d1_size, d1_ld
    integer :: lbound1_, range1_(2)
    integer :: d2_start, d2_end, d2_size, d2_ld
    integer :: lbound2_, range2_(2)
    integer :: d3_start, d3_end, d3_size, d3_ld
    integer :: lbound3_, range3_(2)
    integer :: d4_start, d4_end, d4_size, d4_ld
    integer :: lbound4_, range4_(2)
    !
    lbound1_=1
    if (present(lbound1)) lbound1_=lbound1 
    range1_=(/1,size(array_out, 1)/)
    if (present(range1)) range1_=range1 
    !
    d1_start = range1_(1) -lbound1_ +1
    d1_end   = range1_(2) -lbound1_ +1
    d1_size  = range1_(2) -range1_(1) + 1
    d1_ld    = size(array_out, 1)
    lbound2_=1
    if (present(lbound2)) lbound2_=lbound2 
    range2_=(/1,size(array_out, 2)/)
    if (present(range2)) range2_=range2 
    !
    d2_start = range2_(1) -lbound2_ +1
    d2_end   = range2_(2) -lbound2_ +1
    d2_size  = range2_(2) -range2_(1) + 1
    d2_ld    = size(array_out, 2)
    lbound3_=1
    if (present(lbound3)) lbound3_=lbound3 
    range3_=(/1,size(array_out, 3)/)
    if (present(range3)) range3_=range3 
    !
    d3_start = range3_(1) -lbound3_ +1
    d3_end   = range3_(2) -lbound3_ +1
    d3_size  = range3_(2) -range3_(1) + 1
    d3_ld    = size(array_out, 3)
    lbound4_=1
    if (present(lbound4)) lbound4_=lbound4 
    range4_=(/1,size(array_out, 4)/)
    if (present(range4)) range4_=range4 
    !
    d4_start = range4_(1) -lbound4_ +1
    d4_end   = range4_(2) -lbound4_ +1
    d4_size  = range4_(2) -range4_(1) + 1
    d4_ld    = size(array_out, 4)
    !
#if defined(__CUDA)
    !call dxlib_errore('cu_memsync_','4D arrays not implemented yet',1)
    array_out(d1_start:d1_end,d2_start:d2_end,d3_start:d3_end,d4_start:d4_end) = &
              array_in(d1_start:d1_end,d2_start:d2_end,d3_start:d3_end,d4_start:d4_end)
#else
    array_out(d1_start:d1_end,d2_start:d2_end,d3_start:d3_end,d4_start:d4_end) = &
              array_in(d1_start:d1_end,d2_start:d2_end,d3_start:d3_end,d4_start:d4_end)
#endif
    !
end subroutine sp_h2d_memsync_c4d
!
!
subroutine dp_d2h_memsync_r1d(array_out, array_in, &
                                             range1, lbound1  )
#if defined(__CUDA)
    use cudafor
#endif
    implicit none
    !
    integer, parameter :: PRCSN = selected_real_kind(14,200)
    real(PRCSN), intent(inout) :: array_out(:)
    real(PRCSN), intent(in)    :: array_in(:)
    integer, optional, intent(in) ::  range1(2)
    integer, optional, intent(in) ::  lbound1
#if defined(__CUDA)
    attributes(device) :: array_in
    integer :: ierr
#endif
    !
    integer :: d1_start, d1_end, d1_size, d1_ld
    integer :: lbound1_, range1_(2)
    !
    lbound1_=1
    if (present(lbound1)) lbound1_=lbound1 
    range1_=(/1,size(array_out, 1)/)
    if (present(range1)) range1_=range1 
    !
    d1_start = range1_(1) -lbound1_ +1
    d1_end   = range1_(2) -lbound1_ +1
    d1_size  = range1_(2) -range1_(1) + 1
    d1_ld    = size(array_out, 1)
    !
#if defined(__CUDA)
    ierr = cudaMemcpy( array_out(d1_start), array_in(d1_start), d1_size, cudaMemcpyDeviceToHost )
#else
    array_out(d1_start:d1_end) = &
              array_in(d1_start:d1_end)
#endif
    !
end subroutine dp_d2h_memsync_r1d
!
subroutine dp_d2h_memsync_r2d(array_out, array_in, &
                                             range1, lbound1 , &
                                             range2, lbound2  )
#if defined(__CUDA)
    use cudafor
#endif
    implicit none
    !
    integer, parameter :: PRCSN = selected_real_kind(14,200)
    real(PRCSN), intent(inout) :: array_out(:,:)
    real(PRCSN), intent(in)    :: array_in(:,:)
    integer, optional, intent(in) ::  range1(2), range2(2)
    integer, optional, intent(in) ::  lbound1, lbound2
#if defined(__CUDA)
    attributes(device) :: array_in
    integer :: ierr
#endif
    !
    integer :: d1_start, d1_end, d1_size, d1_ld
    integer :: lbound1_, range1_(2)
    integer :: d2_start, d2_end, d2_size, d2_ld
    integer :: lbound2_, range2_(2)
    !
    lbound1_=1
    if (present(lbound1)) lbound1_=lbound1 
    range1_=(/1,size(array_out, 1)/)
    if (present(range1)) range1_=range1 
    !
    d1_start = range1_(1) -lbound1_ +1
    d1_end   = range1_(2) -lbound1_ +1
    d1_size  = range1_(2) -range1_(1) + 1
    d1_ld    = size(array_out, 1)
    lbound2_=1
    if (present(lbound2)) lbound2_=lbound2 
    range2_=(/1,size(array_out, 2)/)
    if (present(range2)) range2_=range2 
    !
    d2_start = range2_(1) -lbound2_ +1
    d2_end   = range2_(2) -lbound2_ +1
    d2_size  = range2_(2) -range2_(1) + 1
    d2_ld    = size(array_out, 2)
    !
#if defined(__CUDA)
    ierr = cudaMemcpy2D( array_out(d1_start, d2_start) , d1_ld, array_in(d1_start, d2_start), d2_ld, d1_size, d2_size )
      !CALL dev_memsync( pinned_buffer, sc_d, (/ nb1, notcnv /), 1, (/ n_start, n_end-n_start+1 /), 1 )
      CALL dev_memsync( sc_d, pinned_buffer, (/ nb1, nb1 + notcnv /), 1, (/ n_start, n_end /), 1 )
           cudaMemcpy2D( pinned_buffer(nb1, n_start) , nvecx, sc_d( nb1, n_start ), nvecx, notcnv, n_end-n_start+1 )
#else
    array_out(d1_start:d1_end,d2_start:d2_end) = &
              array_in(d1_start:d1_end,d2_start:d2_end)
#endif
    !
end subroutine dp_d2h_memsync_r2d
!
subroutine dp_d2h_memsync_r3d(array_out, array_in, &
                                             range1, lbound1 , &
                                             range2, lbound2 , &
                                             range3, lbound3  )
#if defined(__CUDA)
    use cudafor
#endif
    implicit none
    !
    integer, parameter :: PRCSN = selected_real_kind(14,200)
    real(PRCSN), intent(inout) :: array_out(:,:,:)
    real(PRCSN), intent(in)    :: array_in(:,:,:)
    integer, optional, intent(in) ::  range1(2), range2(2), range3(2)
    integer, optional, intent(in) ::  lbound1, lbound2, lbound3
#if defined(__CUDA)
    attributes(device) :: array_in
    integer :: ierr
#endif
    !
    integer :: d1_start, d1_end, d1_size, d1_ld
    integer :: lbound1_, range1_(2)
    integer :: d2_start, d2_end, d2_size, d2_ld
    integer :: lbound2_, range2_(2)
    integer :: d3_start, d3_end, d3_size, d3_ld
    integer :: lbound3_, range3_(2)
    !
    lbound1_=1
    if (present(lbound1)) lbound1_=lbound1 
    range1_=(/1,size(array_out, 1)/)
    if (present(range1)) range1_=range1 
    !
    d1_start = range1_(1) -lbound1_ +1
    d1_end   = range1_(2) -lbound1_ +1
    d1_size  = range1_(2) -range1_(1) + 1
    d1_ld    = size(array_out, 1)
    lbound2_=1
    if (present(lbound2)) lbound2_=lbound2 
    range2_=(/1,size(array_out, 2)/)
    if (present(range2)) range2_=range2 
    !
    d2_start = range2_(1) -lbound2_ +1
    d2_end   = range2_(2) -lbound2_ +1
    d2_size  = range2_(2) -range2_(1) + 1
    d2_ld    = size(array_out, 2)
    lbound3_=1
    if (present(lbound3)) lbound3_=lbound3 
    range3_=(/1,size(array_out, 3)/)
    if (present(range3)) range3_=range3 
    !
    d3_start = range3_(1) -lbound3_ +1
    d3_end   = range3_(2) -lbound3_ +1
    d3_size  = range3_(2) -range3_(1) + 1
    d3_ld    = size(array_out, 3)
    !
#if defined(__CUDA)
    array_out(d1_start:d1_end,d2_start:d2_end,d3_start:d3_end) = &
              array_in(d1_start:d1_end,d2_start:d2_end,d3_start:d3_end)
#else
    array_out(d1_start:d1_end,d2_start:d2_end,d3_start:d3_end) = &
              array_in(d1_start:d1_end,d2_start:d2_end,d3_start:d3_end)
#endif
    !
end subroutine dp_d2h_memsync_r3d
!
subroutine dp_d2h_memsync_r4d(array_out, array_in, &
                                             range1, lbound1 , &
                                             range2, lbound2 , &
                                             range3, lbound3 , &
                                             range4, lbound4  )
#if defined(__CUDA)
    use cudafor
#endif
    implicit none
    !
    integer, parameter :: PRCSN = selected_real_kind(14,200)
    real(PRCSN), intent(inout) :: array_out(:,:,:,:)
    real(PRCSN), intent(in)    :: array_in(:,:,:,:)
    integer, optional, intent(in) ::  range1(2), range2(2), range3(2), range4(2)
    integer, optional, intent(in) ::  lbound1, lbound2, lbound3, lbound4
#if defined(__CUDA)
    attributes(device) :: array_in
    integer :: ierr
#endif
    !
    integer :: d1_start, d1_end, d1_size, d1_ld
    integer :: lbound1_, range1_(2)
    integer :: d2_start, d2_end, d2_size, d2_ld
    integer :: lbound2_, range2_(2)
    integer :: d3_start, d3_end, d3_size, d3_ld
    integer :: lbound3_, range3_(2)
    integer :: d4_start, d4_end, d4_size, d4_ld
    integer :: lbound4_, range4_(2)
    !
    lbound1_=1
    if (present(lbound1)) lbound1_=lbound1 
    range1_=(/1,size(array_out, 1)/)
    if (present(range1)) range1_=range1 
    !
    d1_start = range1_(1) -lbound1_ +1
    d1_end   = range1_(2) -lbound1_ +1
    d1_size  = range1_(2) -range1_(1) + 1
    d1_ld    = size(array_out, 1)
    lbound2_=1
    if (present(lbound2)) lbound2_=lbound2 
    range2_=(/1,size(array_out, 2)/)
    if (present(range2)) range2_=range2 
    !
    d2_start = range2_(1) -lbound2_ +1
    d2_end   = range2_(2) -lbound2_ +1
    d2_size  = range2_(2) -range2_(1) + 1
    d2_ld    = size(array_out, 2)
    lbound3_=1
    if (present(lbound3)) lbound3_=lbound3 
    range3_=(/1,size(array_out, 3)/)
    if (present(range3)) range3_=range3 
    !
    d3_start = range3_(1) -lbound3_ +1
    d3_end   = range3_(2) -lbound3_ +1
    d3_size  = range3_(2) -range3_(1) + 1
    d3_ld    = size(array_out, 3)
    lbound4_=1
    if (present(lbound4)) lbound4_=lbound4 
    range4_=(/1,size(array_out, 4)/)
    if (present(range4)) range4_=range4 
    !
    d4_start = range4_(1) -lbound4_ +1
    d4_end   = range4_(2) -lbound4_ +1
    d4_size  = range4_(2) -range4_(1) + 1
    d4_ld    = size(array_out, 4)
    !
#if defined(__CUDA)
    array_out(d1_start:d1_end,d2_start:d2_end,d3_start:d3_end,d4_start:d4_end) = &
              array_in(d1_start:d1_end,d2_start:d2_end,d3_start:d3_end,d4_start:d4_end)
#else
    array_out(d1_start:d1_end,d2_start:d2_end,d3_start:d3_end,d4_start:d4_end) = &
              array_in(d1_start:d1_end,d2_start:d2_end,d3_start:d3_end,d4_start:d4_end)
#endif
    !
end subroutine dp_d2h_memsync_r4d
!
subroutine dp_d2h_memsync_c1d(array_out, array_in, &
                                             range1, lbound1  )
#if defined(__CUDA)
    use cudafor
#endif
    implicit none
    !
    integer, parameter :: PRCSN = selected_real_kind(14,200)
    complex(PRCSN), intent(inout) :: array_out(:)
    complex(PRCSN), intent(in)    :: array_in(:)
    integer, optional, intent(in) ::  range1(2)
    integer, optional, intent(in) ::  lbound1
#if defined(__CUDA)
    attributes(device) :: array_in
    integer :: ierr
#endif
    !
    integer :: d1_start, d1_end, d1_size, d1_ld
    integer :: lbound1_, range1_(2)
    !
    lbound1_=1
    if (present(lbound1)) lbound1_=lbound1 
    range1_=(/1,size(array_out, 1)/)
    if (present(range1)) range1_=range1 
    !
    d1_start = range1_(1) -lbound1_ +1
    d1_end   = range1_(2) -lbound1_ +1
    d1_size  = range1_(2) -range1_(1) + 1
    d1_ld    = size(array_out, 1)
    !
#if defined(__CUDA)
    ierr = cudaMemcpy( array_out(d1_start), array_in(d1_start), d1_size, cudaMemcpyDeviceToHost )
#else
    array_out(d1_start:d1_end) = &
              array_in(d1_start:d1_end)
#endif
    !
end subroutine dp_d2h_memsync_c1d
!
subroutine dp_d2h_memsync_c2d(array_out, array_in, &
                                             range1, lbound1 , &
                                             range2, lbound2  )
#if defined(__CUDA)
    use cudafor
#endif
    implicit none
    !
    integer, parameter :: PRCSN = selected_real_kind(14,200)
    complex(PRCSN), intent(inout) :: array_out(:,:)
    complex(PRCSN), intent(in)    :: array_in(:,:)
    integer, optional, intent(in) ::  range1(2), range2(2)
    integer, optional, intent(in) ::  lbound1, lbound2
#if defined(__CUDA)
    attributes(device) :: array_in
    integer :: ierr
#endif
    !
    integer :: d1_start, d1_end, d1_size, d1_ld
    integer :: lbound1_, range1_(2)
    integer :: d2_start, d2_end, d2_size, d2_ld
    integer :: lbound2_, range2_(2)
    !
    lbound1_=1
    if (present(lbound1)) lbound1_=lbound1 
    range1_=(/1,size(array_out, 1)/)
    if (present(range1)) range1_=range1 
    !
    d1_start = range1_(1) -lbound1_ +1
    d1_end   = range1_(2) -lbound1_ +1
    d1_size  = range1_(2) -range1_(1) + 1
    d1_ld    = size(array_out, 1)
    lbound2_=1
    if (present(lbound2)) lbound2_=lbound2 
    range2_=(/1,size(array_out, 2)/)
    if (present(range2)) range2_=range2 
    !
    d2_start = range2_(1) -lbound2_ +1
    d2_end   = range2_(2) -lbound2_ +1
    d2_size  = range2_(2) -range2_(1) + 1
    d2_ld    = size(array_out, 2)
    !
#if defined(__CUDA)
    ierr = cudaMemcpy2D( array_out(d1_start, d2_start) , d1_ld, array_in(d1_start, d2_start), d2_ld, d1_size, d2_size )
#else
    array_out(d1_start:d1_end,d2_start:d2_end) = &
              array_in(d1_start:d1_end,d2_start:d2_end)
#endif
    !
end subroutine dp_d2h_memsync_c2d
!
subroutine dp_d2h_memsync_c3d(array_out, array_in, &
                                             range1, lbound1 , &
                                             range2, lbound2 , &
                                             range3, lbound3  )
#if defined(__CUDA)
    use cudafor
#endif
    implicit none
    !
    integer, parameter :: PRCSN = selected_real_kind(14,200)
    complex(PRCSN), intent(inout) :: array_out(:,:,:)
    complex(PRCSN), intent(in)    :: array_in(:,:,:)
    integer, optional, intent(in) ::  range1(2), range2(2), range3(2)
    integer, optional, intent(in) ::  lbound1, lbound2, lbound3
#if defined(__CUDA)
    attributes(device) :: array_in
    integer :: ierr
#endif
    !
    integer :: d1_start, d1_end, d1_size, d1_ld
    integer :: lbound1_, range1_(2)
    integer :: d2_start, d2_end, d2_size, d2_ld
    integer :: lbound2_, range2_(2)
    integer :: d3_start, d3_end, d3_size, d3_ld
    integer :: lbound3_, range3_(2)
    !
    lbound1_=1
    if (present(lbound1)) lbound1_=lbound1 
    range1_=(/1,size(array_out, 1)/)
    if (present(range1)) range1_=range1 
    !
    d1_start = range1_(1) -lbound1_ +1
    d1_end   = range1_(2) -lbound1_ +1
    d1_size  = range1_(2) -range1_(1) + 1
    d1_ld    = size(array_out, 1)
    lbound2_=1
    if (present(lbound2)) lbound2_=lbound2 
    range2_=(/1,size(array_out, 2)/)
    if (present(range2)) range2_=range2 
    !
    d2_start = range2_(1) -lbound2_ +1
    d2_end   = range2_(2) -lbound2_ +1
    d2_size  = range2_(2) -range2_(1) + 1
    d2_ld    = size(array_out, 2)
    lbound3_=1
    if (present(lbound3)) lbound3_=lbound3 
    range3_=(/1,size(array_out, 3)/)
    if (present(range3)) range3_=range3 
    !
    d3_start = range3_(1) -lbound3_ +1
    d3_end   = range3_(2) -lbound3_ +1
    d3_size  = range3_(2) -range3_(1) + 1
    d3_ld    = size(array_out, 3)
    !
#if defined(__CUDA)
    array_out(d1_start:d1_end,d2_start:d2_end,d3_start:d3_end) = &
              array_in(d1_start:d1_end,d2_start:d2_end,d3_start:d3_end)
#else
    array_out(d1_start:d1_end,d2_start:d2_end,d3_start:d3_end) = &
              array_in(d1_start:d1_end,d2_start:d2_end,d3_start:d3_end)
#endif
    !
end subroutine dp_d2h_memsync_c3d
!
subroutine dp_d2h_memsync_c4d(array_out, array_in, &
                                             range1, lbound1 , &
                                             range2, lbound2 , &
                                             range3, lbound3 , &
                                             range4, lbound4  )
#if defined(__CUDA)
    use cudafor
#endif
    implicit none
    !
    integer, parameter :: PRCSN = selected_real_kind(14,200)
    complex(PRCSN), intent(inout) :: array_out(:,:,:,:)
    complex(PRCSN), intent(in)    :: array_in(:,:,:,:)
    integer, optional, intent(in) ::  range1(2), range2(2), range3(2), range4(2)
    integer, optional, intent(in) ::  lbound1, lbound2, lbound3, lbound4
#if defined(__CUDA)
    attributes(device) :: array_in
    integer :: ierr
#endif
    !
    integer :: d1_start, d1_end, d1_size, d1_ld
    integer :: lbound1_, range1_(2)
    integer :: d2_start, d2_end, d2_size, d2_ld
    integer :: lbound2_, range2_(2)
    integer :: d3_start, d3_end, d3_size, d3_ld
    integer :: lbound3_, range3_(2)
    integer :: d4_start, d4_end, d4_size, d4_ld
    integer :: lbound4_, range4_(2)
    !
    lbound1_=1
    if (present(lbound1)) lbound1_=lbound1 
    range1_=(/1,size(array_out, 1)/)
    if (present(range1)) range1_=range1 
    !
    d1_start = range1_(1) -lbound1_ +1
    d1_end   = range1_(2) -lbound1_ +1
    d1_size  = range1_(2) -range1_(1) + 1
    d1_ld    = size(array_out, 1)
    lbound2_=1
    if (present(lbound2)) lbound2_=lbound2 
    range2_=(/1,size(array_out, 2)/)
    if (present(range2)) range2_=range2 
    !
    d2_start = range2_(1) -lbound2_ +1
    d2_end   = range2_(2) -lbound2_ +1
    d2_size  = range2_(2) -range2_(1) + 1
    d2_ld    = size(array_out, 2)
    lbound3_=1
    if (present(lbound3)) lbound3_=lbound3 
    range3_=(/1,size(array_out, 3)/)
    if (present(range3)) range3_=range3 
    !
    d3_start = range3_(1) -lbound3_ +1
    d3_end   = range3_(2) -lbound3_ +1
    d3_size  = range3_(2) -range3_(1) + 1
    d3_ld    = size(array_out, 3)
    lbound4_=1
    if (present(lbound4)) lbound4_=lbound4 
    range4_=(/1,size(array_out, 4)/)
    if (present(range4)) range4_=range4 
    !
    d4_start = range4_(1) -lbound4_ +1
    d4_end   = range4_(2) -lbound4_ +1
    d4_size  = range4_(2) -range4_(1) + 1
    d4_ld    = size(array_out, 4)
    !
#if defined(__CUDA)
    array_out(d1_start:d1_end,d2_start:d2_end,d3_start:d3_end,d4_start:d4_end) = &
              array_in(d1_start:d1_end,d2_start:d2_end,d3_start:d3_end,d4_start:d4_end)
#else
    array_out(d1_start:d1_end,d2_start:d2_end,d3_start:d3_end,d4_start:d4_end) = &
              array_in(d1_start:d1_end,d2_start:d2_end,d3_start:d3_end,d4_start:d4_end)
#endif
    !
end subroutine dp_d2h_memsync_c4d
!
subroutine sp_d2h_memsync_r1d(array_out, array_in, &
                                             range1, lbound1  )
#if defined(__CUDA)
    use cudafor
#endif
    implicit none
    !
    integer, parameter :: PRCSN = selected_real_kind(6, 37)
    real(PRCSN), intent(inout) :: array_out(:)
    real(PRCSN), intent(in)    :: array_in(:)
    integer, optional, intent(in) ::  range1(2)
    integer, optional, intent(in) ::  lbound1
#if defined(__CUDA)
    attributes(device) :: array_in
    integer :: ierr
#endif
    !
    integer :: d1_start, d1_end, d1_size, d1_ld
    integer :: lbound1_, range1_(2)
    !
    lbound1_=1
    if (present(lbound1)) lbound1_=lbound1 
    range1_=(/1,size(array_out, 1)/)
    if (present(range1)) range1_=range1 
    !
    d1_start = range1_(1) -lbound1_ +1
    d1_end   = range1_(2) -lbound1_ +1
    d1_size  = range1_(2) -range1_(1) + 1
    d1_ld    = size(array_out, 1)
    !
#if defined(__CUDA)
    ierr = cudaMemcpy( array_out(d1_start), array_in(d1_start), d1_size, cudaMemcpyDeviceToHost )
#else
    array_out(d1_start:d1_end) = &
              array_in(d1_start:d1_end)
#endif
    !
end subroutine sp_d2h_memsync_r1d
!
subroutine sp_d2h_memsync_r2d(array_out, array_in, &
                                             range1, lbound1 , &
                                             range2, lbound2  )
#if defined(__CUDA)
    use cudafor
#endif
    implicit none
    !
    integer, parameter :: PRCSN = selected_real_kind(6, 37)
    real(PRCSN), intent(inout) :: array_out(:,:)
    real(PRCSN), intent(in)    :: array_in(:,:)
    integer, optional, intent(in) ::  range1(2), range2(2)
    integer, optional, intent(in) ::  lbound1, lbound2
#if defined(__CUDA)
    attributes(device) :: array_in
    integer :: ierr
#endif
    !
    integer :: d1_start, d1_end, d1_size, d1_ld
    integer :: lbound1_, range1_(2)
    integer :: d2_start, d2_end, d2_size, d2_ld
    integer :: lbound2_, range2_(2)
    !
    lbound1_=1
    if (present(lbound1)) lbound1_=lbound1 
    range1_=(/1,size(array_out, 1)/)
    if (present(range1)) range1_=range1 
    !
    d1_start = range1_(1) -lbound1_ +1
    d1_end   = range1_(2) -lbound1_ +1
    d1_size  = range1_(2) -range1_(1) + 1
    d1_ld    = size(array_out, 1)
    lbound2_=1
    if (present(lbound2)) lbound2_=lbound2 
    range2_=(/1,size(array_out, 2)/)
    if (present(range2)) range2_=range2 
    !
    d2_start = range2_(1) -lbound2_ +1
    d2_end   = range2_(2) -lbound2_ +1
    d2_size  = range2_(2) -range2_(1) + 1
    d2_ld    = size(array_out, 2)
    !
#if defined(__CUDA)
    ierr = cudaMemcpy2D( array_out(d1_start, d2_start) , d1_ld, array_in(d1_start, d2_start), d2_ld, d1_size, d2_size )
#else
    array_out(d1_start:d1_end,d2_start:d2_end) = &
              array_in(d1_start:d1_end,d2_start:d2_end)
#endif
    !
end subroutine sp_d2h_memsync_r2d
!
subroutine sp_d2h_memsync_r3d(array_out, array_in, &
                                             range1, lbound1 , &
                                             range2, lbound2 , &
                                             range3, lbound3  )
#if defined(__CUDA)
    use cudafor
#endif
    implicit none
    !
    integer, parameter :: PRCSN = selected_real_kind(6, 37)
    real(PRCSN), intent(inout) :: array_out(:,:,:)
    real(PRCSN), intent(in)    :: array_in(:,:,:)
    integer, optional, intent(in) ::  range1(2), range2(2), range3(2)
    integer, optional, intent(in) ::  lbound1, lbound2, lbound3
#if defined(__CUDA)
    attributes(device) :: array_in
    integer :: ierr
#endif
    !
    integer :: d1_start, d1_end, d1_size, d1_ld
    integer :: lbound1_, range1_(2)
    integer :: d2_start, d2_end, d2_size, d2_ld
    integer :: lbound2_, range2_(2)
    integer :: d3_start, d3_end, d3_size, d3_ld
    integer :: lbound3_, range3_(2)
    !
    lbound1_=1
    if (present(lbound1)) lbound1_=lbound1 
    range1_=(/1,size(array_out, 1)/)
    if (present(range1)) range1_=range1 
    !
    d1_start = range1_(1) -lbound1_ +1
    d1_end   = range1_(2) -lbound1_ +1
    d1_size  = range1_(2) -range1_(1) + 1
    d1_ld    = size(array_out, 1)
    lbound2_=1
    if (present(lbound2)) lbound2_=lbound2 
    range2_=(/1,size(array_out, 2)/)
    if (present(range2)) range2_=range2 
    !
    d2_start = range2_(1) -lbound2_ +1
    d2_end   = range2_(2) -lbound2_ +1
    d2_size  = range2_(2) -range2_(1) + 1
    d2_ld    = size(array_out, 2)
    lbound3_=1
    if (present(lbound3)) lbound3_=lbound3 
    range3_=(/1,size(array_out, 3)/)
    if (present(range3)) range3_=range3 
    !
    d3_start = range3_(1) -lbound3_ +1
    d3_end   = range3_(2) -lbound3_ +1
    d3_size  = range3_(2) -range3_(1) + 1
    d3_ld    = size(array_out, 3)
    !
#if defined(__CUDA)
    array_out(d1_start:d1_end,d2_start:d2_end,d3_start:d3_end) = &
              array_in(d1_start:d1_end,d2_start:d2_end,d3_start:d3_end)
#else
    array_out(d1_start:d1_end,d2_start:d2_end,d3_start:d3_end) = &
              array_in(d1_start:d1_end,d2_start:d2_end,d3_start:d3_end)
#endif
    !
end subroutine sp_d2h_memsync_r3d
!
subroutine sp_d2h_memsync_r4d(array_out, array_in, &
                                             range1, lbound1 , &
                                             range2, lbound2 , &
                                             range3, lbound3 , &
                                             range4, lbound4  )
#if defined(__CUDA)
    use cudafor
#endif
    implicit none
    !
    integer, parameter :: PRCSN = selected_real_kind(6, 37)
    real(PRCSN), intent(inout) :: array_out(:,:,:,:)
    real(PRCSN), intent(in)    :: array_in(:,:,:,:)
    integer, optional, intent(in) ::  range1(2), range2(2), range3(2), range4(2)
    integer, optional, intent(in) ::  lbound1, lbound2, lbound3, lbound4
#if defined(__CUDA)
    attributes(device) :: array_in
    integer :: ierr
#endif
    !
    integer :: d1_start, d1_end, d1_size, d1_ld
    integer :: lbound1_, range1_(2)
    integer :: d2_start, d2_end, d2_size, d2_ld
    integer :: lbound2_, range2_(2)
    integer :: d3_start, d3_end, d3_size, d3_ld
    integer :: lbound3_, range3_(2)
    integer :: d4_start, d4_end, d4_size, d4_ld
    integer :: lbound4_, range4_(2)
    !
    lbound1_=1
    if (present(lbound1)) lbound1_=lbound1 
    range1_=(/1,size(array_out, 1)/)
    if (present(range1)) range1_=range1 
    !
    d1_start = range1_(1) -lbound1_ +1
    d1_end   = range1_(2) -lbound1_ +1
    d1_size  = range1_(2) -range1_(1) + 1
    d1_ld    = size(array_out, 1)
    lbound2_=1
    if (present(lbound2)) lbound2_=lbound2 
    range2_=(/1,size(array_out, 2)/)
    if (present(range2)) range2_=range2 
    !
    d2_start = range2_(1) -lbound2_ +1
    d2_end   = range2_(2) -lbound2_ +1
    d2_size  = range2_(2) -range2_(1) + 1
    d2_ld    = size(array_out, 2)
    lbound3_=1
    if (present(lbound3)) lbound3_=lbound3 
    range3_=(/1,size(array_out, 3)/)
    if (present(range3)) range3_=range3 
    !
    d3_start = range3_(1) -lbound3_ +1
    d3_end   = range3_(2) -lbound3_ +1
    d3_size  = range3_(2) -range3_(1) + 1
    d3_ld    = size(array_out, 3)
    lbound4_=1
    if (present(lbound4)) lbound4_=lbound4 
    range4_=(/1,size(array_out, 4)/)
    if (present(range4)) range4_=range4 
    !
    d4_start = range4_(1) -lbound4_ +1
    d4_end   = range4_(2) -lbound4_ +1
    d4_size  = range4_(2) -range4_(1) + 1
    d4_ld    = size(array_out, 4)
    !
#if defined(__CUDA)
    array_out(d1_start:d1_end,d2_start:d2_end,d3_start:d3_end,d4_start:d4_end) = &
              array_in(d1_start:d1_end,d2_start:d2_end,d3_start:d3_end,d4_start:d4_end)
#else
    array_out(d1_start:d1_end,d2_start:d2_end,d3_start:d3_end,d4_start:d4_end) = &
              array_in(d1_start:d1_end,d2_start:d2_end,d3_start:d3_end,d4_start:d4_end)
#endif
    !
end subroutine sp_d2h_memsync_r4d
!
subroutine sp_d2h_memsync_c1d(array_out, array_in, &
                                             range1, lbound1  )
#if defined(__CUDA)
    use cudafor
#endif
    implicit none
    !
    integer, parameter :: PRCSN = selected_real_kind(6, 37)
    complex(PRCSN), intent(inout) :: array_out(:)
    complex(PRCSN), intent(in)    :: array_in(:)
    integer, optional, intent(in) ::  range1(2)
    integer, optional, intent(in) ::  lbound1
#if defined(__CUDA)
    attributes(device) :: array_in
    integer :: ierr
#endif
    !
    integer :: d1_start, d1_end, d1_size, d1_ld
    integer :: lbound1_, range1_(2)
    !
    lbound1_=1
    if (present(lbound1)) lbound1_=lbound1 
    range1_=(/1,size(array_out, 1)/)
    if (present(range1)) range1_=range1 
    !
    d1_start = range1_(1) -lbound1_ +1
    d1_end   = range1_(2) -lbound1_ +1
    d1_size  = range1_(2) -range1_(1) + 1
    d1_ld    = size(array_out, 1)
    !
#if defined(__CUDA)
    ierr = cudaMemcpy( array_out(d1_start), array_in(d1_start), d1_size, cudaMemcpyDeviceToHost )
#else
    array_out(d1_start:d1_end) = &
              array_in(d1_start:d1_end)
#endif
    !
end subroutine sp_d2h_memsync_c1d
!
subroutine sp_d2h_memsync_c2d(array_out, array_in, &
                                             range1, lbound1 , &
                                             range2, lbound2  )
#if defined(__CUDA)
    use cudafor
#endif
    implicit none
    !
    integer, parameter :: PRCSN = selected_real_kind(6, 37)
    complex(PRCSN), intent(inout) :: array_out(:,:)
    complex(PRCSN), intent(in)    :: array_in(:,:)
    integer, optional, intent(in) ::  range1(2), range2(2)
    integer, optional, intent(in) ::  lbound1, lbound2
#if defined(__CUDA)
    attributes(device) :: array_in
    integer :: ierr
#endif
    !
    integer :: d1_start, d1_end, d1_size, d1_ld
    integer :: lbound1_, range1_(2)
    integer :: d2_start, d2_end, d2_size, d2_ld
    integer :: lbound2_, range2_(2)
    !
    lbound1_=1
    if (present(lbound1)) lbound1_=lbound1 
    range1_=(/1,size(array_out, 1)/)
    if (present(range1)) range1_=range1 
    !
    d1_start = range1_(1) -lbound1_ +1
    d1_end   = range1_(2) -lbound1_ +1
    d1_size  = range1_(2) -range1_(1) + 1
    d1_ld    = size(array_out, 1)
    lbound2_=1
    if (present(lbound2)) lbound2_=lbound2 
    range2_=(/1,size(array_out, 2)/)
    if (present(range2)) range2_=range2 
    !
    d2_start = range2_(1) -lbound2_ +1
    d2_end   = range2_(2) -lbound2_ +1
    d2_size  = range2_(2) -range2_(1) + 1
    d2_ld    = size(array_out, 2)
    !
#if defined(__CUDA)
    ierr = cudaMemcpy2D( array_out(d1_start, d2_start) , d1_ld, array_in(d1_start, d2_start), d2_ld, d1_size, d2_size )
#else
    array_out(d1_start:d1_end,d2_start:d2_end) = &
              array_in(d1_start:d1_end,d2_start:d2_end)
#endif
    !
end subroutine sp_d2h_memsync_c2d
!
subroutine sp_d2h_memsync_c3d(array_out, array_in, &
                                             range1, lbound1 , &
                                             range2, lbound2 , &
                                             range3, lbound3  )
#if defined(__CUDA)
    use cudafor
#endif
    implicit none
    !
    integer, parameter :: PRCSN = selected_real_kind(6, 37)
    complex(PRCSN), intent(inout) :: array_out(:,:,:)
    complex(PRCSN), intent(in)    :: array_in(:,:,:)
    integer, optional, intent(in) ::  range1(2), range2(2), range3(2)
    integer, optional, intent(in) ::  lbound1, lbound2, lbound3
#if defined(__CUDA)
    attributes(device) :: array_in
    integer :: ierr
#endif
    !
    integer :: d1_start, d1_end, d1_size, d1_ld
    integer :: lbound1_, range1_(2)
    integer :: d2_start, d2_end, d2_size, d2_ld
    integer :: lbound2_, range2_(2)
    integer :: d3_start, d3_end, d3_size, d3_ld
    integer :: lbound3_, range3_(2)
    !
    lbound1_=1
    if (present(lbound1)) lbound1_=lbound1 
    range1_=(/1,size(array_out, 1)/)
    if (present(range1)) range1_=range1 
    !
    d1_start = range1_(1) -lbound1_ +1
    d1_end   = range1_(2) -lbound1_ +1
    d1_size  = range1_(2) -range1_(1) + 1
    d1_ld    = size(array_out, 1)
    lbound2_=1
    if (present(lbound2)) lbound2_=lbound2 
    range2_=(/1,size(array_out, 2)/)
    if (present(range2)) range2_=range2 
    !
    d2_start = range2_(1) -lbound2_ +1
    d2_end   = range2_(2) -lbound2_ +1
    d2_size  = range2_(2) -range2_(1) + 1
    d2_ld    = size(array_out, 2)
    lbound3_=1
    if (present(lbound3)) lbound3_=lbound3 
    range3_=(/1,size(array_out, 3)/)
    if (present(range3)) range3_=range3 
    !
    d3_start = range3_(1) -lbound3_ +1
    d3_end   = range3_(2) -lbound3_ +1
    d3_size  = range3_(2) -range3_(1) + 1
    d3_ld    = size(array_out, 3)
    !
#if defined(__CUDA)
    array_out(d1_start:d1_end,d2_start:d2_end,d3_start:d3_end) = &
              array_in(d1_start:d1_end,d2_start:d2_end,d3_start:d3_end)
#else
    array_out(d1_start:d1_end,d2_start:d2_end,d3_start:d3_end) = &
              array_in(d1_start:d1_end,d2_start:d2_end,d3_start:d3_end)
#endif
    !
end subroutine sp_d2h_memsync_c3d
!
subroutine sp_d2h_memsync_c4d(array_out, array_in, &
                                             range1, lbound1 , &
                                             range2, lbound2 , &
                                             range3, lbound3 , &
                                             range4, lbound4  )
#if defined(__CUDA)
    use cudafor
#endif
    implicit none
    !
    integer, parameter :: PRCSN = selected_real_kind(6, 37)
    complex(PRCSN), intent(inout) :: array_out(:,:,:,:)
    complex(PRCSN), intent(in)    :: array_in(:,:,:,:)
    integer, optional, intent(in) ::  range1(2), range2(2), range3(2), range4(2)
    integer, optional, intent(in) ::  lbound1, lbound2, lbound3, lbound4
#if defined(__CUDA)
    attributes(device) :: array_in
    integer :: ierr
#endif
    !
    integer :: d1_start, d1_end, d1_size, d1_ld
    integer :: lbound1_, range1_(2)
    integer :: d2_start, d2_end, d2_size, d2_ld
    integer :: lbound2_, range2_(2)
    integer :: d3_start, d3_end, d3_size, d3_ld
    integer :: lbound3_, range3_(2)
    integer :: d4_start, d4_end, d4_size, d4_ld
    integer :: lbound4_, range4_(2)
    !
    lbound1_=1
    if (present(lbound1)) lbound1_=lbound1 
    range1_=(/1,size(array_out, 1)/)
    if (present(range1)) range1_=range1 
    !
    d1_start = range1_(1) -lbound1_ +1
    d1_end   = range1_(2) -lbound1_ +1
    d1_size  = range1_(2) -range1_(1) + 1
    d1_ld    = size(array_out, 1)
    lbound2_=1
    if (present(lbound2)) lbound2_=lbound2 
    range2_=(/1,size(array_out, 2)/)
    if (present(range2)) range2_=range2 
    !
    d2_start = range2_(1) -lbound2_ +1
    d2_end   = range2_(2) -lbound2_ +1
    d2_size  = range2_(2) -range2_(1) + 1
    d2_ld    = size(array_out, 2)
    lbound3_=1
    if (present(lbound3)) lbound3_=lbound3 
    range3_=(/1,size(array_out, 3)/)
    if (present(range3)) range3_=range3 
    !
    d3_start = range3_(1) -lbound3_ +1
    d3_end   = range3_(2) -lbound3_ +1
    d3_size  = range3_(2) -range3_(1) + 1
    d3_ld    = size(array_out, 3)
    lbound4_=1
    if (present(lbound4)) lbound4_=lbound4 
    range4_=(/1,size(array_out, 4)/)
    if (present(range4)) range4_=range4 
    !
    d4_start = range4_(1) -lbound4_ +1
    d4_end   = range4_(2) -lbound4_ +1
    d4_size  = range4_(2) -range4_(1) + 1
    d4_ld    = size(array_out, 4)
    !
#if defined(__CUDA)
    array_out(d1_start:d1_end,d2_start:d2_end,d3_start:d3_end,d4_start:d4_end) = &
              array_in(d1_start:d1_end,d2_start:d2_end,d3_start:d3_end,d4_start:d4_end)
#else
    array_out(d1_start:d1_end,d2_start:d2_end,d3_start:d3_end,d4_start:d4_end) = &
              array_in(d1_start:d1_end,d2_start:d2_end,d3_start:d3_end,d4_start:d4_end)
#endif
    !
end subroutine sp_d2h_memsync_c4d
!


subroutine dp_dev_conjg_c1d(array_inout, &
                                    
                                   range1, lbound1 )
    implicit none
    !
    integer, parameter :: PRCSN = selected_real_kind(14,200)
    complex(PRCSN), intent(inout) :: array_inout(:)
    integer, optional, intent(in) ::  range1(2)
    integer, optional, intent(in) ::  lbound1
#if defined(__CUDA)
    attributes(device) :: array_inout
#endif
    !
    integer :: i1, d1s, d1e
    integer :: lbound1_, range1_(2)
    !
    lbound1_=1
    if (present(lbound1)) lbound1_=lbound1 
    range1_=(/1,size(array_inout, 1)/)
    if (present(range1)) range1_=range1 
    !
    d1s = range1_(1) -lbound1_ +1
    d1e = range1_(2) -lbound1_ +1
    !
    !
    !$cuf kernel do(1)
    do i1 = d1s, d1e
        array_inout(i1 ) = &
            conjg (  array_inout (i1 ) )
    enddo
    !
end subroutine dp_dev_conjg_c1d
!
subroutine dp_dev_conjg_c2d(array_inout, &
                                    
                                   range1, lbound1, & 
                                   range2, lbound2 )
    implicit none
    !
    integer, parameter :: PRCSN = selected_real_kind(14,200)
    complex(PRCSN), intent(inout) :: array_inout(:,:)
    integer, optional, intent(in) ::  range1(2), range2(2)
    integer, optional, intent(in) ::  lbound1, lbound2
#if defined(__CUDA)
    attributes(device) :: array_inout
#endif
    !
    integer :: i1, d1s, d1e
    integer :: lbound1_, range1_(2)
    integer :: i2, d2s, d2e
    integer :: lbound2_, range2_(2)
    !
    lbound1_=1
    if (present(lbound1)) lbound1_=lbound1 
    range1_=(/1,size(array_inout, 1)/)
    if (present(range1)) range1_=range1 
    !
    d1s = range1_(1) -lbound1_ +1
    d1e = range1_(2) -lbound1_ +1
    !
    lbound2_=1
    if (present(lbound2)) lbound2_=lbound2 
    range2_=(/1,size(array_inout, 2)/)
    if (present(range2)) range2_=range2 
    !
    d2s = range2_(1) -lbound2_ +1
    d2e = range2_(2) -lbound2_ +1
    !
    !
    !$cuf kernel do(2)
    do i2 = d2s, d2e
    do i1 = d1s, d1e
        array_inout(i1,i2 ) = &
            conjg (  array_inout (i1,i2 ) )
    enddo
    enddo
    !
end subroutine dp_dev_conjg_c2d
!
subroutine dp_dev_conjg_c3d(array_inout, &
                                    
                                   range1, lbound1, & 
                                   range2, lbound2, & 
                                   range3, lbound3 )
    implicit none
    !
    integer, parameter :: PRCSN = selected_real_kind(14,200)
    complex(PRCSN), intent(inout) :: array_inout(:,:,:)
    integer, optional, intent(in) ::  range1(2), range2(2), range3(2)
    integer, optional, intent(in) ::  lbound1, lbound2, lbound3
#if defined(__CUDA)
    attributes(device) :: array_inout
#endif
    !
    integer :: i1, d1s, d1e
    integer :: lbound1_, range1_(2)
    integer :: i2, d2s, d2e
    integer :: lbound2_, range2_(2)
    integer :: i3, d3s, d3e
    integer :: lbound3_, range3_(2)
    !
    lbound1_=1
    if (present(lbound1)) lbound1_=lbound1 
    range1_=(/1,size(array_inout, 1)/)
    if (present(range1)) range1_=range1 
    !
    d1s = range1_(1) -lbound1_ +1
    d1e = range1_(2) -lbound1_ +1
    !
    lbound2_=1
    if (present(lbound2)) lbound2_=lbound2 
    range2_=(/1,size(array_inout, 2)/)
    if (present(range2)) range2_=range2 
    !
    d2s = range2_(1) -lbound2_ +1
    d2e = range2_(2) -lbound2_ +1
    !
    lbound3_=1
    if (present(lbound3)) lbound3_=lbound3 
    range3_=(/1,size(array_inout, 3)/)
    if (present(range3)) range3_=range3 
    !
    d3s = range3_(1) -lbound3_ +1
    d3e = range3_(2) -lbound3_ +1
    !
    !
    !$cuf kernel do(3)
    do i3 = d3s, d3e
    do i2 = d2s, d2e
    do i1 = d1s, d1e
        array_inout(i1,i2,i3 ) = &
            conjg (  array_inout (i1,i2,i3 ) )
    enddo
    enddo
    enddo
    !
end subroutine dp_dev_conjg_c3d
!
subroutine dp_dev_conjg_c4d(array_inout, &
                                    
                                   range1, lbound1, & 
                                   range2, lbound2, & 
                                   range3, lbound3, & 
                                   range4, lbound4 )
    implicit none
    !
    integer, parameter :: PRCSN = selected_real_kind(14,200)
    complex(PRCSN), intent(inout) :: array_inout(:,:,:,:)
    integer, optional, intent(in) ::  range1(2), range2(2), range3(2), range4(2)
    integer, optional, intent(in) ::  lbound1, lbound2, lbound3, lbound4
#if defined(__CUDA)
    attributes(device) :: array_inout
#endif
    !
    integer :: i1, d1s, d1e
    integer :: lbound1_, range1_(2)
    integer :: i2, d2s, d2e
    integer :: lbound2_, range2_(2)
    integer :: i3, d3s, d3e
    integer :: lbound3_, range3_(2)
    integer :: i4, d4s, d4e
    integer :: lbound4_, range4_(2)
    !
    lbound1_=1
    if (present(lbound1)) lbound1_=lbound1 
    range1_=(/1,size(array_inout, 1)/)
    if (present(range1)) range1_=range1 
    !
    d1s = range1_(1) -lbound1_ +1
    d1e = range1_(2) -lbound1_ +1
    !
    lbound2_=1
    if (present(lbound2)) lbound2_=lbound2 
    range2_=(/1,size(array_inout, 2)/)
    if (present(range2)) range2_=range2 
    !
    d2s = range2_(1) -lbound2_ +1
    d2e = range2_(2) -lbound2_ +1
    !
    lbound3_=1
    if (present(lbound3)) lbound3_=lbound3 
    range3_=(/1,size(array_inout, 3)/)
    if (present(range3)) range3_=range3 
    !
    d3s = range3_(1) -lbound3_ +1
    d3e = range3_(2) -lbound3_ +1
    !
    lbound4_=1
    if (present(lbound4)) lbound4_=lbound4 
    range4_=(/1,size(array_inout, 4)/)
    if (present(range4)) range4_=range4 
    !
    d4s = range4_(1) -lbound4_ +1
    d4e = range4_(2) -lbound4_ +1
    !
    !
    !$cuf kernel do(4)
    do i4 = d4s, d4e
    do i3 = d3s, d3e
    do i2 = d2s, d2e
    do i1 = d1s, d1e
        array_inout(i1,i2,i3,i4 ) = &
            conjg (  array_inout (i1,i2,i3,i4 ) )
    enddo
    enddo
    enddo
    enddo
    !
end subroutine dp_dev_conjg_c4d
!
subroutine sp_dev_conjg_c1d(array_inout, &
                                    
                                   range1, lbound1 )
    implicit none
    !
    integer, parameter :: PRCSN = selected_real_kind(6, 37)
    complex(PRCSN), intent(inout) :: array_inout(:)
    integer, optional, intent(in) ::  range1(2)
    integer, optional, intent(in) ::  lbound1
#if defined(__CUDA)
    attributes(device) :: array_inout
#endif
    !
    integer :: i1, d1s, d1e
    integer :: lbound1_, range1_(2)
    !
    lbound1_=1
    if (present(lbound1)) lbound1_=lbound1 
    range1_=(/1,size(array_inout, 1)/)
    if (present(range1)) range1_=range1 
    !
    d1s = range1_(1) -lbound1_ +1
    d1e = range1_(2) -lbound1_ +1
    !
    !
    !$cuf kernel do(1)
    do i1 = d1s, d1e
        array_inout(i1 ) = &
            conjg (  array_inout (i1 ) )
    enddo
    !
end subroutine sp_dev_conjg_c1d
!
subroutine sp_dev_conjg_c2d(array_inout, &
                                    
                                   range1, lbound1, & 
                                   range2, lbound2 )
    implicit none
    !
    integer, parameter :: PRCSN = selected_real_kind(6, 37)
    complex(PRCSN), intent(inout) :: array_inout(:,:)
    integer, optional, intent(in) ::  range1(2), range2(2)
    integer, optional, intent(in) ::  lbound1, lbound2
#if defined(__CUDA)
    attributes(device) :: array_inout
#endif
    !
    integer :: i1, d1s, d1e
    integer :: lbound1_, range1_(2)
    integer :: i2, d2s, d2e
    integer :: lbound2_, range2_(2)
    !
    lbound1_=1
    if (present(lbound1)) lbound1_=lbound1 
    range1_=(/1,size(array_inout, 1)/)
    if (present(range1)) range1_=range1 
    !
    d1s = range1_(1) -lbound1_ +1
    d1e = range1_(2) -lbound1_ +1
    !
    lbound2_=1
    if (present(lbound2)) lbound2_=lbound2 
    range2_=(/1,size(array_inout, 2)/)
    if (present(range2)) range2_=range2 
    !
    d2s = range2_(1) -lbound2_ +1
    d2e = range2_(2) -lbound2_ +1
    !
    !
    !$cuf kernel do(2)
    do i2 = d2s, d2e
    do i1 = d1s, d1e
        array_inout(i1,i2 ) = &
            conjg (  array_inout (i1,i2 ) )
    enddo
    enddo
    !
end subroutine sp_dev_conjg_c2d
!
subroutine sp_dev_conjg_c3d(array_inout, &
                                    
                                   range1, lbound1, & 
                                   range2, lbound2, & 
                                   range3, lbound3 )
    implicit none
    !
    integer, parameter :: PRCSN = selected_real_kind(6, 37)
    complex(PRCSN), intent(inout) :: array_inout(:,:,:)
    integer, optional, intent(in) ::  range1(2), range2(2), range3(2)
    integer, optional, intent(in) ::  lbound1, lbound2, lbound3
#if defined(__CUDA)
    attributes(device) :: array_inout
#endif
    !
    integer :: i1, d1s, d1e
    integer :: lbound1_, range1_(2)
    integer :: i2, d2s, d2e
    integer :: lbound2_, range2_(2)
    integer :: i3, d3s, d3e
    integer :: lbound3_, range3_(2)
    !
    lbound1_=1
    if (present(lbound1)) lbound1_=lbound1 
    range1_=(/1,size(array_inout, 1)/)
    if (present(range1)) range1_=range1 
    !
    d1s = range1_(1) -lbound1_ +1
    d1e = range1_(2) -lbound1_ +1
    !
    lbound2_=1
    if (present(lbound2)) lbound2_=lbound2 
    range2_=(/1,size(array_inout, 2)/)
    if (present(range2)) range2_=range2 
    !
    d2s = range2_(1) -lbound2_ +1
    d2e = range2_(2) -lbound2_ +1
    !
    lbound3_=1
    if (present(lbound3)) lbound3_=lbound3 
    range3_=(/1,size(array_inout, 3)/)
    if (present(range3)) range3_=range3 
    !
    d3s = range3_(1) -lbound3_ +1
    d3e = range3_(2) -lbound3_ +1
    !
    !
    !$cuf kernel do(3)
    do i3 = d3s, d3e
    do i2 = d2s, d2e
    do i1 = d1s, d1e
        array_inout(i1,i2,i3 ) = &
            conjg (  array_inout (i1,i2,i3 ) )
    enddo
    enddo
    enddo
    !
end subroutine sp_dev_conjg_c3d
!
subroutine sp_dev_conjg_c4d(array_inout, &
                                    
                                   range1, lbound1, & 
                                   range2, lbound2, & 
                                   range3, lbound3, & 
                                   range4, lbound4 )
    implicit none
    !
    integer, parameter :: PRCSN = selected_real_kind(6, 37)
    complex(PRCSN), intent(inout) :: array_inout(:,:,:,:)
    integer, optional, intent(in) ::  range1(2), range2(2), range3(2), range4(2)
    integer, optional, intent(in) ::  lbound1, lbound2, lbound3, lbound4
#if defined(__CUDA)
    attributes(device) :: array_inout
#endif
    !
    integer :: i1, d1s, d1e
    integer :: lbound1_, range1_(2)
    integer :: i2, d2s, d2e
    integer :: lbound2_, range2_(2)
    integer :: i3, d3s, d3e
    integer :: lbound3_, range3_(2)
    integer :: i4, d4s, d4e
    integer :: lbound4_, range4_(2)
    !
    lbound1_=1
    if (present(lbound1)) lbound1_=lbound1 
    range1_=(/1,size(array_inout, 1)/)
    if (present(range1)) range1_=range1 
    !
    d1s = range1_(1) -lbound1_ +1
    d1e = range1_(2) -lbound1_ +1
    !
    lbound2_=1
    if (present(lbound2)) lbound2_=lbound2 
    range2_=(/1,size(array_inout, 2)/)
    if (present(range2)) range2_=range2 
    !
    d2s = range2_(1) -lbound2_ +1
    d2e = range2_(2) -lbound2_ +1
    !
    lbound3_=1
    if (present(lbound3)) lbound3_=lbound3 
    range3_=(/1,size(array_inout, 3)/)
    if (present(range3)) range3_=range3 
    !
    d3s = range3_(1) -lbound3_ +1
    d3e = range3_(2) -lbound3_ +1
    !
    lbound4_=1
    if (present(lbound4)) lbound4_=lbound4 
    range4_=(/1,size(array_inout, 4)/)
    if (present(range4)) range4_=range4 
    !
    d4s = range4_(1) -lbound4_ +1
    d4e = range4_(2) -lbound4_ +1
    !
    !
    !$cuf kernel do(4)
    do i4 = d4s, d4e
    do i3 = d3s, d3e
    do i2 = d2s, d2e
    do i1 = d1s, d1e
        array_inout(i1,i2,i3,i4 ) = &
            conjg (  array_inout (i1,i2,i3,i4 ) )
    enddo
    enddo
    enddo
    enddo
    !
end subroutine sp_dev_conjg_c4d
!


subroutine dxlib_errore( calling_routine, message, ierr )
    implicit none
    character(len=*), intent(in) :: calling_routine, message
    ! the name of the calling calling_routine
    ! the output message
    integer,          intent(in) :: ierr
    !
    print *, trim(calling_routine), trim(message), ierr
    !
end subroutine dxlib_errore
