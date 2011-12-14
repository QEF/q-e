! Copyright (C) 2006-2008 Dmitry Korotin dmitry@korotin.name
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#define ZERO (0.d0,0.d0)
#define ONE (1.d0,0.d0)

!----------------------------------------------------------------------
subroutine wannier_occupancies(occ)
  !----------------------------------------------------------------------
  !    
  ! ... This routine computes occupation of each wannier. It is assumed that WF generated  already and stored if the buffer.
  !
  use kinds, only: DP
  use wannier_new, only: nwan, pp
  use io_global, only : stdout 
  use wvfct, only: nbnd, et, wg
  use klist, only: nks
  use lsda_mod, only: current_spin, lsda, nspin, isk
  USE io_files
  USE buffers
  
  implicit none
  real(DP), intent(out) :: occ(nwan,nwan,nspin)

  integer :: i,j,k,ik

  occ = ZERO
  current_spin = 1
 
  DO ik=1, nks
     IF (lsda) current_spin  = isk(ik)
     CALL get_buffer( pp, nwordwpp, iunwpp, ik)
     DO i=1, nwan
        DO j=1,nwan
           DO k=1, nbnd
              occ(i,j,current_spin) = occ(i,j,current_spin) + pp(i,k)*conjg(pp(j,k))*wg(k,ik)
           END DO
        END DO
     END DO
  END DO

  IF(nspin.eq.1) occ=occ*0.5D0

  return
end subroutine wannier_occupancies
