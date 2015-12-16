!
! Copyright (C) 2001-2015 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
SUBROUTINE lr_diagonalise(iter)
  !---------------------------------------------------------------------
  ! Brent Walker, ICTP, 2004
  !---------------------------------------------------------------------
  ! ... diagonalises the coefficient matrix from the Lanczos using
  ! ... LAPACK/BLAS routines
  !---------------------------------------------------------------------
  !
  ! I. Timrov: This routine is not used any longer in turboTDDFPT...
  !
  USE io_global,            ONLY : stdout
  USE kinds,                ONLY : dp
  USE lr_variables,   ONLY : lr_verbosity
  !
  IMPLICIT NONE
  !
  INTEGER,INTENT(in) :: iter
  !
  INTEGER :: dimen
  real(kind=dp),ALLOCATABLE :: coeff_mat(:,:)
  !
  dimen=2*iter
  ALLOCATE(coeff_mat(dimen,dimen))
  coeff_mat(:,:)=0.0d0
  CALL lr_build_matrix_spectrum(coeff_mat,iter)
  !
  CALL lr_diagonalise_matrix(coeff_mat,dimen)
  !
  DEALLOCATE(coeff_mat)
  !
  RETURN
END SUBROUTINE lr_diagonalise
!-----------------------------------------------------------------------
SUBROUTINE lr_diagonalise_matrix(coeff_mat,dimen)
  !
  !
  USE io_global,            ONLY : stdout
  USE kinds,                ONLY : dp
  USE lr_variables,         ONLY : vl,vr,eval1,eval2
   USE io_global,      ONLY : stdout
 !
  IMPLICIT NONE
  !
  ! input variables
  INTEGER,INTENT(in) :: dimen
  real(kind=dp),INTENT(in) :: coeff_mat(dimen,dimen)
  !
  ! local variables
  real(kind=dp),ALLOCATABLE :: work(:)
  INTEGER :: info
  INTEGER :: i,j
  real(kind=dp) :: temp1,temp2
  real(kind=dp),ALLOCATABLE :: vl_temp(:),vr_temp(:)
  !
  IF (lr_verbosity > 5) THEN
     WRITE(stdout,'("<lr_diagonalise>")')
  ENDIF
  IF(allocated(eval1)) DEALLOCATE(eval1)
  IF(allocated(eval2)) DEALLOCATE(eval2)
  IF(allocated(vl)) DEALLOCATE(vl)
  IF(allocated(vr)) DEALLOCATE(vr)
  !
  ALLOCATE(eval1(dimen))
  ALLOCATE(eval2(dimen))
  ALLOCATE(vl(dimen,dimen))
  ALLOCATE(vr(dimen,dimen))
  eval1(:)=0.0d0
  eval2(:)=0.0d0
  vl(:,:)=0.0d0
  vr(:,:)=0.0d0
  !
  ALLOCATE(work(8*dimen))
  ALLOCATE(vl_temp(dimen),vr_temp(dimen))
  work(:)=0.0d0
  vl_temp(:)=0.0d0
  vr_temp(:)=0.0d0
  !
  info=0
  !
  CALL dgeev('v','v',dimen,coeff_mat,dimen,eval1,eval2, &
       vl,dimen,vr,dimen,work,8*dimen,info)
  !
  IF (info/=0) THEN
     CALL errore(' lr_main ', 'Diagonalisation of coefficient ' // &
          & 'matrix unsuccessful',1)
  ENDIF
  !
  IF(.true.) THEN
     ! sort the eigenvalues (inefficient)
     DO i=1,dimen
        DO j=1,i
           IF(eval1(i)<eval1(j)) THEN
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
           ENDIF
        ENDDO
     ENDDO
  ENDIF
  !
  WRITE(stdout,'(5X,"# < start of eigenvalue listing >")')
  DO i=1,dimen
     WRITE(stdout,'(5X,i5,2(1X,f20.12))') i,eval1(i),eval2(i)
  ENDDO
  WRITE(stdout,'(5X,"# < end of eigenvalue listing >")')
  !
  DEALLOCATE(work)
  DEALLOCATE(vl_temp)
  DEALLOCATE(vr_temp)
  !
  RETURN
  !
END SUBROUTINE lr_diagonalise_matrix
!-----------------------------------------------------------------------
SUBROUTINE lr_build_matrix_spectrum(coeff_mat,iter)
  !---------------------------------------------------------------------
  ! ... version for non-hermitian lanczos scheme
  !---------------------------------------------------------------------
  !
  !
  USE io_global,            ONLY : stdout
  USE kinds,                ONLY : dp
  USE lr_variables,         ONLY : a,b
  !
  IMPLICIT NONE
  !
  INTEGER,INTENT(in) :: iter
  real(kind=dp) :: coeff_mat(2*iter,2*iter)
  !
  INTEGER :: i,j
  !
  DO i=1,iter
     coeff_mat(i,i)=1.0d0
  ENDDO
  coeff_mat(1,1)=0.0d0
  coeff_mat(2,2)=0.0d0
  coeff_mat(1,2)=a(1,1)
  coeff_mat(2,1)=a(2,1)
  !
  trirecursion: DO i = 2, iter
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
  ENDDO trirecursion
  !
  IF(.false.) THEN
     DO i=1,2*iter
        DO j=1,2*iter
           WRITE(stdout,'(f10.5)',ADVANCE='no') coeff_mat(i,j)
        ENDDO
        WRITE(stdout,*)
     ENDDO
  ENDIF
  !
  RETURN
END SUBROUTINE lr_build_matrix_spectrum
!-----------------------------------------------------------------------
