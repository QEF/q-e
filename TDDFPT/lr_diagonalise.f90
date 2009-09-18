!-----------------------------------------------------------------------
subroutine lr_diagonalise(iter)
  !---------------------------------------------------------------------
  ! Brent Walker, ICTP, 2004
  !---------------------------------------------------------------------
  ! ... diagonalises the coefficient matrix from the Lanczos using
  ! ... LAPACK/BLAS routines
  !---------------------------------------------------------------------
  !
#include "f_defs.h"
  !
  use io_global,            only : stdout
  use kinds,                only : dp
  USE lr_variables,   ONLY : lr_verbosity
  !
  implicit none
  !
  integer,intent(in) :: iter
  !
  integer :: dimen
  real(kind=dp),allocatable :: coeff_mat(:,:) 
  !
  dimen=2*iter
  allocate(coeff_mat(dimen,dimen))
  coeff_mat(:,:)=0.0d0
  call lr_build_matrix_spectrum(coeff_mat,iter)
  !
  call lr_diagonalise_matrix(coeff_mat,dimen)
  !
  deallocate(coeff_mat)
  !
  return
end subroutine lr_diagonalise
!-----------------------------------------------------------------------
subroutine lr_diagonalise_matrix(coeff_mat,dimen)
  !
#include "f_defs.h"
  !
  use io_global,            only : stdout
  use kinds,                only : dp
  use lr_variables,         only : vl,vr,eval1,eval2
   USE io_global,      ONLY : stdout
 !
  implicit none
  !
  ! input variables
  integer,intent(in) :: dimen
  real(kind=dp),intent(in) :: coeff_mat(dimen,dimen)
  !
  ! local variables
  real(kind=dp),allocatable :: work(:)
  integer :: info
  integer :: i,j
  real(kind=dp) :: temp1,temp2
  real(kind=dp),allocatable :: vl_temp(:),vr_temp(:)
  !
  If (lr_verbosity > 5) THEN
     WRITE(stdout,'("<lr_diagonalise>")')
  endif
  if(allocated(eval1)) deallocate(eval1)
  if(allocated(eval2)) deallocate(eval2)
  if(allocated(vl)) deallocate(vl)
  if(allocated(vr)) deallocate(vr)
  !
  allocate(eval1(dimen))
  allocate(eval2(dimen))
  allocate(vl(dimen,dimen))
  allocate(vr(dimen,dimen))
  eval1(:)=0.0d0
  eval2(:)=0.0d0
  vl(:,:)=0.0d0
  vr(:,:)=0.0d0
  !
  allocate(work(8*dimen))
  allocate(vl_temp(dimen),vr_temp(dimen))
  work(:)=0.0d0
  vl_temp(:)=0.0d0
  vr_temp(:)=0.0d0
  !
  info=0
  !
  call dgeev('v','v',dimen,coeff_mat,dimen,eval1,eval2, &
       vl,dimen,vr,dimen,work,8*dimen,info)
  !
  if (info/=0) then
     call errore(' lr_main ', 'Diagonalisation of coefficient ' // &
          & 'matrix unsuccessful',1)
  endif
  !
  if(.true.) then
     ! sort the eigenvalues (inefficient)
     do i=1,dimen
        do j=1,i
           if(eval1(i)<eval1(j)) then
              temp1=eval1(i)
              temp2=eval2(i)
              vl_temp(:)=vl(:,i)
              vr_temp(:)=vr(:,i)
              eval1(i)=eval1(j)
              eval2(i)=eval2(j)
              vl(:,i)=vl(:,j)
              vr(:,i)=vr(:,j)
              eval1(j)=temp1
              eval2(j)=temp2
              vl(:,j)=vl_temp(:)
              vr(:,j)=vr_temp(:)
           endif
        enddo
     enddo
  endif
  !
  write(stdout,'(5X,"# < start of eigenvalue listing >")')
  do i=1,dimen
     write(stdout,'(5X,i5,2(1X,f20.12))') i,eval1(i),eval2(i)
  enddo
  write(stdout,'(5X,"# < end of eigenvalue listing >")')
  !
  deallocate(work)
  deallocate(vl_temp)
  deallocate(vr_temp)
  !
  return
  !
end subroutine lr_diagonalise_matrix
!-----------------------------------------------------------------------
subroutine lr_build_matrix_spectrum(coeff_mat,iter)
  !---------------------------------------------------------------------
  ! ... version for non-hermitian lanczos scheme
  !---------------------------------------------------------------------
  !
#include "f_defs.h"
  !
  use io_global,            only : stdout
  use kinds,                only : dp
  use lr_variables,         only : a,b
  !
  implicit none
  !
  integer,intent(in) :: iter
  real(kind=dp) :: coeff_mat(2*iter,2*iter)
  !
  integer :: i,j
  !
  do i=1,iter
     coeff_mat(i,i)=1.0d0
  enddo
  coeff_mat(1,1)=0.0d0
  coeff_mat(2,2)=0.0d0
  coeff_mat(1,2)=a(1,1)
  coeff_mat(2,1)=a(2,1)
  !
  trirecursion: do i = 2, iter
     coeff_mat(2*i-1,2*i-1)=0.0d0
     coeff_mat(2*i+0,2*i+0)=0.0d0
     coeff_mat(2*i-1,2*i+0)=a(1,i)
     coeff_mat(2*i+0,2*i-1)=a(2,i)
     !    write(*,*) "A and B", i,"-th iteration"
     !    write(*,"(4f12.4)") a(1,i),a(2,i),b(1,i),b(2,i)
     !    write(*,*) ""
     coeff_mat(2*i-3,2*i-1)=b(1,i)
     coeff_mat(2*i-2,2*i+0)=b(2,i)
     coeff_mat(2*i-1,2*i-3)=abs(b(2,i))
     coeff_mat(2*i+0,2*i-2)=abs(b(1,i))
  end do trirecursion
  !
  if(.false.) then
     do i=1,2*iter
        do j=1,2*iter
           write(stdout,'(f10.5,$)') coeff_mat(i,j)
        enddo
        write(stdout,*)
     enddo
  endif
  !
  return
end subroutine lr_build_matrix_spectrum
!-----------------------------------------------------------------------
