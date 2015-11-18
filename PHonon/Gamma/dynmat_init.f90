!
! Copyright (C) 2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
SUBROUTINE dynmat_init
  !-----------------------------------------------------------------------
  !
  !  Calculate part of the terms appearing in the dynamical matrix
  !
  USE ions_base, ONLY : ntyp => nsp, nat, ityp, zv, tau
  USE cell_base, ONLY : at, bg, omega, alat
  USE gvect,     ONLY : ngm, g, gg
  USE cgcom
  IMPLICIT NONE
  real(DP), ALLOCATABLE:: dyn0(:,:),dyn1(:,:), dyncc(:,:)
  INTEGER :: i,j, na,nb
  !
  CALL start_clock('dynmat_init')
  !
  ALLOCATE  ( dyn0 ( 3*nat, nmodes))
  ALLOCATE  ( dyn1 ( 3*nat, nmodes))
  ALLOCATE  ( dyncc( 3*nat, nmodes))
  !
  !  first electronic contribution arising from the term  <psi|d2v|psi>
  !
  CALL rhod2vkb(dyn0)
  !
  !  ionic contribution
  !
  CALL d2ion (nat,ntyp,ityp,zv,tau,alat,omega,                      &
       at,bg,g,gg,ngm,nmodes,u,has_equivalent,dyn1)
  !
  !  core-correction contribution
  !
  CALL dynmatcc(dyncc)
  !
  DO j=1,nmodes
     DO i=1,3*nat
        dyn(i,j)=dyn0(i,j)+dyn1(i,j)+dyncc(i,j)
     ENDDO
  ENDDO
  !
  DEALLOCATE(dyncc)
  DEALLOCATE(dyn1 )
  DEALLOCATE(dyn0 )
  !
  CALL stop_clock('dynmat_init')
  !
  RETURN
END SUBROUTINE dynmat_init
