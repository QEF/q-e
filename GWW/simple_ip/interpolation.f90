! This subroutine uses the trilinear interpolation to interpolate the non-local part of the pseudopotential.
! Algorithm from Wikipedia


subroutine  trilinear_parallel_nc(kptns,ik,sh,c)
  USE kinds, ONLY : DP 
  USE constants, ONLY : pi   ! DEBUG
  USE io_global, ONLY : stdout ! DEBUG
  USE simple_ip_objects
  !
  IMPLICIT NONE
  !
  TYPE(kpoints) :: kptns
  TYPE(shirley) :: sh
  INTEGER, INTENT(in) :: ik ! dense local k-grid index
  COMPLEX(kind=DP), INTENT(out) :: c(sh%nkb,sh%npol,sh%ntot_e)  ! Trilinear interpolated projector (non-collinear)
  COMPLEX(kind=DP), DIMENSION(:,:,:), ALLOCATABLE :: c000, c100, c010, c110, c001, c101, c011, c111 
  COMPLEX(kind=DP), DIMENSION(:,:,:), ALLOCATABLE :: c00, c01, c10, c11, c0, c1
  REAL(kind=DP) ::  delta(3)

  allocate(c000(sh%nkb,sh%npol,sh%ntot_e),c100(sh%nkb,sh%npol,sh%ntot_e),&
  & c010(sh%nkb,sh%npol,sh%ntot_e),c110(sh%nkb,sh%npol,sh%ntot_e), &
  & c001(sh%nkb,sh%npol,sh%ntot_e), c101(sh%nkb,sh%npol,sh%ntot_e), &
  & c011(sh%nkb,sh%npol,sh%ntot_e), c111(sh%nkb,sh%npol,sh%ntot_e) )

  allocate(c00(sh%nkb,sh%npol,sh%ntot_e),c01(sh%nkb,sh%npol,sh%ntot_e),&
  & c10(sh%nkb,sh%npol,sh%ntot_e),c11(sh%nkb,sh%npol,sh%ntot_e), &
  & c0(sh%nkb,sh%npol,sh%ntot_e), c1(sh%nkb,sh%npol,sh%ntot_e))
  
  c000(1:sh%nkb,1:sh%npol,1:sh%ntot_e) = sh%beck_nc(1:sh%nkb,1:sh%npol,1:sh%ntot_e, kptns%coord_cube(1,ik))
  c100(1:sh%nkb,1:sh%npol,1:sh%ntot_e) = sh%beck_nc(1:sh%nkb,1:sh%npol,1:sh%ntot_e, kptns%coord_cube(2,ik))
  c010(1:sh%nkb,1:sh%npol,1:sh%ntot_e) = sh%beck_nc(1:sh%nkb,1:sh%npol,1:sh%ntot_e, kptns%coord_cube(3,ik))
  c110(1:sh%nkb,1:sh%npol,1:sh%ntot_e) = sh%beck_nc(1:sh%nkb,1:sh%npol,1:sh%ntot_e, kptns%coord_cube(4,ik))
  c001(1:sh%nkb,1:sh%npol,1:sh%ntot_e) = sh%beck_nc(1:sh%nkb,1:sh%npol,1:sh%ntot_e, kptns%coord_cube(5,ik))
  c101(1:sh%nkb,1:sh%npol,1:sh%ntot_e) = sh%beck_nc(1:sh%nkb,1:sh%npol,1:sh%ntot_e, kptns%coord_cube(6,ik))
  c011(1:sh%nkb,1:sh%npol,1:sh%ntot_e) = sh%beck_nc(1:sh%nkb,1:sh%npol,1:sh%ntot_e, kptns%coord_cube(7,ik))
  c111(1:sh%nkb,1:sh%npol,1:sh%ntot_e) = sh%beck_nc(1:sh%nkb,1:sh%npol,1:sh%ntot_e, kptns%coord_cube(8,ik))
     
  delta(1:3) = kptns%pos_cube(1:3,ik)
  
  c00 = c000*(1.0 - delta(1)) + c100*delta(1)
  c01 = c001*(1.0 - delta(1)) + c101*delta(1)
  c10 = c010*(1.0 - delta(1)) + c110*delta(1)
  c11 = c011*(1.0 - delta(1)) + c111*delta(1)

  c0 = c00*(1.0 - delta(2)) + c10*delta(2)
  c1 = c01*(1.0 - delta(2)) + c11*delta(2)
  
  c = c0*(1.0 - delta(3)) + c1*delta(3)

  deallocate(c000, c100, c010, c110, c001, c101, c011, c111, c00, c01, c10, c11, c0, c1)

end subroutine trilinear_parallel_nc

subroutine  trilinear_parallelc(kptns,ik,sh,c)
  USE kinds, ONLY : DP 
  USE simple_ip_objects

  implicit none

  TYPE(kpoints) :: kptns
  TYPE(shirley) :: sh
  INTEGER, INTENT(in) :: ik ! dense local k-grid index
  COMPLEX(kind=DP), INTENT(out) :: c(sh%nkb,sh%ntot_e)  ! Trilinear interpolated projector (collinear)
  COMPLEX(kind=DP), DIMENSION(:,:), ALLOCATABLE :: c000, c100, c010, c110, c001, c101, c011, c111 
  COMPLEX(kind=DP), DIMENSION(:,:), ALLOCATABLE :: c00, c01, c10, c11, c0, c1
  REAL(kind=DP) ::  delta(3)

  allocate(c000(sh%nkb,sh%ntot_e),c100(sh%nkb,sh%ntot_e),&
  & c010(sh%nkb,sh%ntot_e),c110(sh%nkb,sh%ntot_e), &
  & c001(sh%nkb,sh%ntot_e), c101(sh%nkb,sh%ntot_e), &
  & c011(sh%nkb,sh%ntot_e), c111(sh%nkb,sh%ntot_e) )

  allocate(c00(sh%nkb,sh%ntot_e),c01(sh%nkb,sh%ntot_e),&
  & c10(sh%nkb,sh%ntot_e),c11(sh%nkb,sh%ntot_e), &
  & c0(sh%nkb,sh%ntot_e), c1(sh%nkb,sh%ntot_e))
 
  c000(1:sh%nkb,1:sh%ntot_e) = sh%beckc(1:sh%nkb,1:sh%ntot_e, kptns%coord_cube(1,ik))
  c100(1:sh%nkb,1:sh%ntot_e) = sh%beckc(1:sh%nkb,1:sh%ntot_e, kptns%coord_cube(2,ik))
  c010(1:sh%nkb,1:sh%ntot_e) = sh%beckc(1:sh%nkb,1:sh%ntot_e, kptns%coord_cube(3,ik))
  c110(1:sh%nkb,1:sh%ntot_e) = sh%beckc(1:sh%nkb,1:sh%ntot_e, kptns%coord_cube(4,ik))
  c001(1:sh%nkb,1:sh%ntot_e) = sh%beckc(1:sh%nkb,1:sh%ntot_e, kptns%coord_cube(5,ik))
  c101(1:sh%nkb,1:sh%ntot_e) = sh%beckc(1:sh%nkb,1:sh%ntot_e, kptns%coord_cube(6,ik))
  c011(1:sh%nkb,1:sh%ntot_e) = sh%beckc(1:sh%nkb,1:sh%ntot_e, kptns%coord_cube(7,ik))
  c111(1:sh%nkb,1:sh%ntot_e) = sh%beckc(1:sh%nkb,1:sh%ntot_e, kptns%coord_cube(8,ik))
     
  delta(1:3) = kptns%pos_cube(1:3,ik)
  
  c00 = c000*(1.0 - delta(1)) + c100*delta(1)
  c01 = c001*(1.0 - delta(1)) + c101*delta(1)
  c10 = c010*(1.0 - delta(1)) + c110*delta(1)
  c11 = c011*(1.0 - delta(1)) + c111*delta(1)

  c0 = c00*(1.0 - delta(2)) + c10*delta(2)
  c1 = c01*(1.0 - delta(2)) + c11*delta(2)
  
  c = c0*(1.0 - delta(3)) + c1*delta(3)

  deallocate(c000, c100, c010, c110, c001, c101, c011, c111, c00, c01, c10, c11, c0, c1)

end subroutine trilinear_parallelc

subroutine  trilinear_parallel_commut(kptns,ik,sh,c)
  USE kinds, ONLY : DP 
  USE simple_ip_objects

  implicit none

  TYPE(kpoints) :: kptns
  TYPE(shirley) :: sh
  INTEGER, INTENT(in) :: ik ! dense local k-grid index
  COMPLEX(kind=DP), INTENT(out) :: c(sh%ntot_e,sh%ntot_e,3)  ! Trilinear interpolated commutator
  COMPLEX(kind=DP), DIMENSION(:,:,:), ALLOCATABLE :: c000, c100, c010, c110, c001, c101, c011, c111
  COMPLEX(kind=DP), DIMENSION(:,:,:), ALLOCATABLE :: c00, c01, c10, c11, c0, c1
  REAL(kind=DP) ::  delta(3)

  allocate(c000(sh%ntot_e,sh%ntot_e,3),c100(sh%ntot_e,sh%ntot_e,3),&
  & c010(sh%ntot_e,sh%ntot_e,3),c110(sh%ntot_e,sh%ntot_e,3), &
  & c001(sh%ntot_e,sh%ntot_e,3), c101(sh%ntot_e,sh%ntot_e,3), &
  & c011(sh%ntot_e,sh%ntot_e,3), c111(sh%ntot_e,sh%ntot_e,3) )

  allocate(c00(sh%ntot_e,sh%ntot_e,3),c01(sh%ntot_e,sh%ntot_e,3),&
  & c10(sh%ntot_e,sh%ntot_e,3),c11(sh%ntot_e,sh%ntot_e,3), &
  & c0(sh%ntot_e,sh%ntot_e,3), c1(sh%ntot_e,sh%ntot_e,3))
 
  c000(1:sh%ntot_e,1:sh%ntot_e,1:3) = sh%commut(1:sh%ntot_e,1:sh%ntot_e,1:3, kptns%coord_cube(1,ik))
  c100(1:sh%ntot_e,1:sh%ntot_e,1:3) = sh%commut(1:sh%ntot_e,1:sh%ntot_e,1:3, kptns%coord_cube(2,ik))
  c010(1:sh%ntot_e,1:sh%ntot_e,1:3) = sh%commut(1:sh%ntot_e,1:sh%ntot_e,1:3, kptns%coord_cube(3,ik))
  c110(1:sh%ntot_e,1:sh%ntot_e,1:3) = sh%commut(1:sh%ntot_e,1:sh%ntot_e,1:3, kptns%coord_cube(4,ik))
  c001(1:sh%ntot_e,1:sh%ntot_e,1:3) = sh%commut(1:sh%ntot_e,1:sh%ntot_e,1:3, kptns%coord_cube(5,ik))
  c101(1:sh%ntot_e,1:sh%ntot_e,1:3) = sh%commut(1:sh%ntot_e,1:sh%ntot_e,1:3, kptns%coord_cube(6,ik))
  c011(1:sh%ntot_e,1:sh%ntot_e,1:3) = sh%commut(1:sh%ntot_e,1:sh%ntot_e,1:3, kptns%coord_cube(7,ik))
  c111(1:sh%ntot_e,1:sh%ntot_e,1:3) = sh%commut(1:sh%ntot_e,1:sh%ntot_e,1:3, kptns%coord_cube(8,ik))
     
  delta(1:3) = kptns%pos_cube(1:3,ik)
  
  c00 = c000*(1.0 - delta(1)) + c100*delta(1)
  c01 = c001*(1.0 - delta(1)) + c101*delta(1)
  c10 = c010*(1.0 - delta(1)) + c110*delta(1)
  c11 = c011*(1.0 - delta(1)) + c111*delta(1)

  c0 = c00*(1.0 - delta(2)) + c10*delta(2)
  c1 = c01*(1.0 - delta(2)) + c11*delta(2)
  
  c = c0*(1.0 - delta(3)) + c1*delta(3)

  deallocate(c000, c100, c010, c110, c001, c101, c011, c111, c00, c01, c10, c11, c0, c1)

end subroutine trilinear_parallel_commut



! OLD ROUTINE (to be deleted)
subroutine  trilinear(q,nkpoints,bg,at,ndim,beck,c)   
  USE kinds, ONLY : DP 
  USE constants, ONLY : pi   ! DEBUG
  USE io_global, ONLY : stdout ! DEBUG

  implicit none

  REAL(kind=DP), INTENT(in) :: bg(3,3) ! reciprocal-space basis vectors
  REAL(kind=DP), INTENT(in) :: at(3,3) ! real-space basis vectors
  REAL(kind=DP), INTENT(in) :: q(3)  ! k-point
  INTEGER, INTENT(in) :: nkpoints(3) ! interpolation grid
  INTEGER, INTENT(in) :: ndim
  INTEGER :: n(3), np(3), i, j
  COMPLEX(kind=DP), INTENT(in) :: beck(ndim,nkpoints(1)*nkpoints(2)*nkpoints(3)) ! Projector 
  COMPLEX(kind=DP), DIMENSION(ndim), INTENT(out) :: c  ! Trilinear interpolated projector

  COMPLEX(kind=DP), DIMENSION(ndim) :: c000, c100, c010, c110, c001, c101, c011, c111 
  COMPLEX(kind=DP), DIMENSION(ndim) :: c00, c01, c10, c11, c0, c1
  REAL(kind=DP) :: qproj(3), delta(3)


  ! Project q on the bg basis (it is obtained doing the scalar product of q with the direct lattice vectors at)
  qproj(1:3) = 0
  do i=1,3
     do j=1,3
        qproj(i) = qproj(i) + q(j)*at(j,i) 
     enddo
  enddo
  
  n(1:3) = 0
  do i=1,3
     n(i) = int(qproj(i)*nkpoints(i))
     np(i) = n(i) + 1
     if ( np(i) >= nkpoints(i) )  np(i) = 0
  enddo   
  
  c000(1:ndim) = beck(1:ndim, n(1)*nkpoints(2)*nkpoints(3) + n(2)*nkpoints(3) + n(3) + 1)
  c100(1:ndim) = beck(1:ndim, np(1)*nkpoints(2)*nkpoints(3) + n(2)*nkpoints(3) + n(3) + 1)
  c010(1:ndim) = beck(1:ndim, n(1)*nkpoints(2)*nkpoints(3) + np(2)*nkpoints(3) + n(3) + 1)
  c110(1:ndim) = beck(1:ndim, np(1)*nkpoints(2)*nkpoints(3) + np(2)*nkpoints(3) + n(3) + 1)
  c001(1:ndim) = beck(1:ndim, n(1)*nkpoints(2)*nkpoints(3) + n(2)*nkpoints(3) + np(3) + 1)
  c101(1:ndim) = beck(1:ndim, np(1)*nkpoints(2)*nkpoints(3) + n(2)*nkpoints(3) + np(3) + 1)
  c011(1:ndim) = beck(1:ndim, n(1)*nkpoints(2)*nkpoints(3) + np(2)*nkpoints(3) + np(3) + 1)
  c111(1:ndim) = beck(1:ndim, np(1)*nkpoints(2)*nkpoints(3) + np(2)*nkpoints(3) + np(3) + 1)
  
  do i=1,3
     delta(i) = ( qproj(i) - dble(n(i))/dble(nkpoints(i)) ) * dble(nkpoints(i)) !  qproj(i) --> q(i)
  enddo

  c00 = c000*(1.0 - delta(1)) + c100*delta(1)
  c01 = c001*(1.0 - delta(1)) + c101*delta(1)
  c10 = c010*(1.0 - delta(1)) + c110*delta(1)
  c11 = c011*(1.0 - delta(1)) + c111*delta(1)

  c0 = c00*(1.0 - delta(2)) + c10*delta(2)
  c1 = c01*(1.0 - delta(2)) + c11*delta(2)
  
  c = c0*(1.0 - delta(3)) + c1*delta(3)

  ! DEBUG
!  write(stdout, *) '**********************INTERPOLATION*************************'
!  write(stdout, *) 'qproj', qproj(1:3)
!
!  write(stdout, *) 'n', n(1:3)
!  write(stdout, *) 'np', np(1:3)
!  write(stdout, *) n(1)*nkpoints(2)*nkpoints(3) + n(2)*nkpoints(3) + n(3) + 1
!  write(stdout, *) np(1)*nkpoints(2)*nkpoints(3) + n(2)*nkpoints(3) + n(3) + 1
!  write(stdout, *) n(1)*nkpoints(2)*nkpoints(3) + np(2)*nkpoints(3) + n(3) + 1
!  write(stdout, *) np(1)*nkpoints(2)*nkpoints(3) + np(2)*nkpoints(3) + n(3) + 1
!  write(stdout, *) n(1)*nkpoints(2)*nkpoints(3) + n(2)*nkpoints(3) + np(3) + 1
!  write(stdout, *) np(1)*nkpoints(2)*nkpoints(3) + n(2)*nkpoints(3) + np(3) + 1
!  write(stdout, *) n(1)*nkpoints(2)*nkpoints(3) + np(2)*nkpoints(3) + np(3) + 1
!  write(stdout, *) np(1)*nkpoints(2)*nkpoints(3) + np(2)*nkpoints(3) + np(3) + 1
!
!  write(stdout, *) 'delta', delta(1:3)
!
!  write(stdout, *) 'c000', c000(1:1)
!  write(stdout, *) 'c100', c100(1:1)
!  write(stdout, *) 'c010', c010(1:1)
!  write(stdout, *) 'c110', c110(1:1)
!  write(stdout, *) 'c001', c001(1:1)
!  write(stdout, *) 'c101', c101(1:1)
!  write(stdout, *) 'c011', c011(1:1)
!  write(stdout, *) 'c111', c111(1:1)
!
!  write(stdout, *) 'c00', c00(1:1)
!  write(stdout, *) 'c01', c01(1:1)
!  write(stdout, *) 'c10', c10(1:1)
!  write(stdout, *) 'c11', c11(1:1)
!  write(stdout, *) 'c0', c0(1:1)
!  write(stdout, *) 'c1', c1(1:1)
!  write(stdout, *) 'c', c(1:1)
  ! DEBUG


  return
end subroutine trilinear



