SUBROUTINE set_psi_in(ik,l,j,e,psi_out)
!
!  This subroutine calculates the all electron wavefunction psi at the 
!  input energy e
!
USE ld1inc
IMPLICIT NONE
INTEGER :: l, ik    ! input: angular momentum and index of the cut-off radius
REAL(DP) :: e, j    ! input: energy and total angular momentum
REAL(DP) :: psi_out(ndm) ! output: the function psi.
REAL(DP) :: psi_dir(ndm,2) ! auxiliary function.
REAL(DP) :: ze2, jnor
integer  :: n

IF (rel == 1) THEN
   CALL lschps(3,zed,exp(dx),dx,mesh,mesh,mesh,1,l,e,psi_out,r,vpot)
ELSEIF (rel == 2) THEN
   CALL dir_outward(ndm,mesh,l,j,e,dx,psi_dir,r,rab,vpot)
   psi_out(:)=psi_dir(:,1)
ELSE
   ze2=-zed*2.0_dp
   CALL intref(l,e,mesh,dx,r,r2,sqr,vpot,ze2,psi_out)
ENDIF
!
!    fix arbitrarily the norm at the cut-off radius equal to 0.5
!
jnor=psi_out(ik)
DO n=1,mesh
   psi_out(n)=psi_out(n)*0.5_dp/jnor
ENDDO

RETURN
END SUBROUTINE set_psi_in

