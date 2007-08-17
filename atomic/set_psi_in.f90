SUBROUTINE set_psi_in(ik,l,j,e,psi_out)
!
!  This subroutine calculates the all electron wavefunction psi at the 
!  input energy e
!
use kinds, only : dp
use radial_grids, only : ndmx
USE ld1inc, only : grid, rel, zed, vpot
IMPLICIT NONE
INTEGER :: l, ik    ! input: angular momentum and index of the cut-off radius
REAL(DP) :: e, j    ! input: energy and total angular momentum
REAL(DP) :: psi_out(ndmx) ! output: the function psi.
REAL(DP) :: psi_dir(ndmx,2) ! auxiliary function.
REAL(DP) :: ze2, jnor
integer  :: n

IF (rel == 1) THEN
   CALL lschps(3,zed,grid,grid%mesh,grid%mesh,1,l,e,psi_out,vpot)
ELSEIF (rel == 2) THEN
   CALL dir_outward(ndmx,grid%mesh,l,j,e,grid%dx,psi_dir,grid%r,grid%rab,vpot)
   psi_out(:)=psi_dir(:,1)
ELSE
   ze2=-zed*2.0_dp
   CALL intref(l,e,grid%mesh,grid,vpot,ze2,psi_out)
ENDIF
!
!    fix arbitrarily the norm at the cut-off radius equal to 0.5
!
jnor=psi_out(ik)
DO n=1,grid%mesh
   psi_out(n)=psi_out(n)*0.5_dp/jnor
ENDDO

RETURN
END SUBROUTINE set_psi_in

