! Copyright (C) 2008 Dmitry Korotin dmitry@korotin.name, Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#define ZERO (0.d0,0.d0)
#define ONE (1.d0,0.d0)


SUBROUTINE ortho_wfc(lda,ldb,wfc,ierr)
!This subroutine orthogonalizes wfcs.

  USE kinds,      ONLY : DP
  USE io_global,  ONLY : stdout
  USE mp_bands,   ONLY : intra_bgrp_comm
  USE mp,         ONLY : mp_sum
  USE noncollin_module, ONLY : noncolin, npol
  
  implicit none
  INTEGER, intent(in) :: lda,ldb
  INTEGER, intent(out) :: ierr
  COMPLEX(DP), intent(inout) :: wfc(lda,ldb)
  
  INTEGER :: i,j,k
  COMPLEX(DP), allocatable :: overlap(:,:),work(:,:), wfc_ortho(:,:)
  REAL(DP) , ALLOCATABLE :: e (:)

  ierr = 0

  ALLOCATE (overlap( lda , lda))    
  ALLOCATE (work   ( lda , lda))    
  ALLOCATE (e      ( lda))    
  ALLOCATE (wfc_ortho( lda , ldb))    
  !
  ! calculate overlap matrix
  !
  overlap = ZERO
  work = ZERO
  e = 0.d0

  CALL ZGEMM ('n', 'c', lda, lda, ldb, (1.d0, 0.d0), &
       wfc, lda, wfc, lda, (0.d0, 0.d0), overlap, lda)
    
  CALL mp_sum(  overlap, intra_bgrp_comm )

  ! find O^-.5
  !
  CALL cdiagh (lda, overlap, lda, e, work)
  DO i = 1, lda
     IF(ABS(e(i)).lt.1.d-10) THEN
        ierr = 1
        RETURN
     ELSE
        e (i) = 1.d0/dsqrt(e(i))
     END IF
  ENDDO
  
  overlap = ZERO
  
  DO i = 1, lda
     DO j = 1, lda
        overlap (i, j) = ZERO
        DO k = 1, lda
           overlap (i, j) = overlap (i, j) + e(k)*work(i, k)*DCONJG(work (j, k) )
        ENDDO
     ENDDO
  ENDDO
  !
  ! trasform wfs O^-.5 psi
  !
  wfc_ortho(:,:) = ZERO
  call ZGEMM('N', 'N', lda, ldb, lda, ONE, overlap, lda, &
       wfc, lda, ZERO, wfc_ortho, lda)
 
  wfc(:,:) = wfc_ortho(:,:)
         
  DEALLOCATE (overlap)
  DEALLOCATE (work)
  DEALLOCATE (e)
  DEALLOCATE (wfc_ortho)
  
  RETURN

END SUBROUTINE

SUBROUTINE check_ortho(lda,ldb,wfc)
!This subroutine checks orthogonality of wfs. Created for debug purposes.

  USE kinds,      ONLY : DP
  USE io_global,  ONLY : stdout
  USE noncollin_module, ONLY : noncolin, npol
  
  implicit none
  INTEGER, intent(in) :: lda,ldb
  COMPLEX(DP), intent(in) :: wfc(lda,ldb)
  
  INTEGER :: i,j,k
  COMPLEX(DP), allocatable :: overlap(:,:)
  
  
  ALLOCATE (overlap( lda , lda))    
  overlap = ZERO
  
  !
  ! calculate overlap matrix
  !
  CALL ZGEMM ('n', 'c', lda, lda, ldb, ONE, &
       wfc, lda, wfc, lda, ZERO, overlap, lda)

  write(stdout,'(5x,a45,2i5)') 'check_ortho for wavefunction with dimentions ', lda,ldb
  do i=1,lda
     write(stdout,'(5x,8f8.4)') (dreal(overlap(i,j)),j=1,lda)
  end do
  write(stdout,'(5x,a18)') 'end of check_ortho'

  DEALLOCATE (overlap)
  
  RETURN

END SUBROUTINE
