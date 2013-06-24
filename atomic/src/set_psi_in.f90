!
! Copyright (C) 2010 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

SUBROUTINE set_psi_in(ik,l,j,e,psi_out,psi_out_rel)
  !
  !  This subroutine calculates the all electron wavefunction psi at the
  !  input energy e
  !
  USE kinds, ONLY : dp
  USE radial_grids, ONLY : ndmx
  USE ld1inc, ONLY : grid, rel, zed, vpot
  IMPLICIT NONE
  INTEGER :: l, ik    ! input: angular momentum and index of the cut-off radius
  REAL(DP) :: e, j    ! input: energy and total angular momentum
  REAL(DP) :: psi_out(ndmx) ! output: the function psi.
  REAL(DP) :: psi_out_rel(ndmx) ! output: the function psi (small component).
  REAL(DP) :: psi_dir(ndmx,2) ! auxiliary function.
  REAL(DP) :: ze2, jnor, thrdum=0.0_dp
  INTEGER  :: n, i, nstop

  psi_out_rel=0.0_DP
  IF (rel == 1) THEN
     n =  grid%mesh
     CALL lschps (3, zed, thrdum, grid, n, 1, l, e, vpot, psi_out, nstop )
  ELSEIF (rel == 2) THEN
     CALL dir_outward (ndmx, grid%mesh, l, j, e, grid%dx, psi_dir, &
          grid%r, grid%rab, vpot )
     psi_out(:)=psi_dir(:,1)
     psi_out_rel(:)=psi_dir(:,2)
  ELSE
     ze2=-zed*2.0_dp
     CALL intref(l,e,grid%mesh,grid,vpot,ze2,psi_out)
  ENDIF
  !
  !    fix arbitrarily the norm at the cut-off radius equal to (about) 0.5**2
  !
  jnor=0.0_dp
  DO n=1,ik
     jnor = jnor + grid%dx*grid%r(n)*psi_out(n)**2
  ENDDO
  jnor = sqrt(jnor)
  DO n=1,grid%mesh
     psi_out(n)=psi_out(n)*0.5_dp/jnor
  ENDDO
  IF (rel==2) THEN
     DO n=1,grid%mesh
        psi_out_rel(n)=psi_out_rel(n)*0.5_dp/jnor
     ENDDO
  ENDIF
  !
  !  Set the wavefunction to zero if it diverges too much. In any case
  !  this part of the scattering wavefunction is not used and generates 
  !  noise in the code.
  !
  Do n=ik+1,grid%mesh
     IF (ABS(psi_out(n))>1.d9) THEN
        DO i=n,grid%mesh
           psi_out(i)=0.0_DP
           IF (rel==2) psi_out_rel(i)=0.0_DP
        ENDDO
     ENDIF
  ENDDO 

  RETURN
END SUBROUTINE set_psi_in

