!
! Copyright (C) 2004-2007 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!=----------------------------------------------------------------------------=!
  MODULE upf_to_internal
!=----------------------------------------------------------------------------=!
  !! Contains two subroutines:
  !! add_upf_grid  generates the radial grid from data contained in upf
  !! set_upf_q     builds the Q(r) functions if not explicitly present
  !!
  USE pseudo_types
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: add_upf_grid, set_upf_q
  SAVE

!=----------------------------------------------------------------------------=!
  CONTAINS
!=----------------------------------------------------------------------------=!
!
!---------------------------------------------------------------------
SUBROUTINE add_upf_grid (upf, grid)
  !---------------------------------------------------------------------
  !
  USE radial_grids, ONLY: radial_grid_type, allocate_radial_grid
  !
  IMPLICIT NONE
  !
  TYPE (pseudo_upf), INTENT(in) :: upf
  TYPE (radial_grid_type), INTENT(out) :: grid
  !
  CALL allocate_radial_grid(grid,upf%mesh)
  !
  grid%dx   = upf%dx
  grid%xmin = upf%xmin
  grid%zmesh= upf%zmesh
  grid%mesh = upf%mesh
  grid%r  (1:upf%mesh) = upf%r  (1:upf%mesh)
  grid%rab(1:upf%mesh) = upf%rab(1:upf%mesh)
  !
  grid%r2 = upf%r**2
  grid%sqr= sqrt(upf%r)
  ! Prevent FP error if r(1) = 0
  IF ( upf%r(1) > 1.0D-16 ) THEN
     grid%rm1 = upf%r**(-1)
     grid%rm2 = upf%r**(-2)
     grid%rm3 = upf%r**(-3)
  ELSE
     grid%rm1(1) =0.0_dp
     grid%rm2(1) =0.0_dp
     grid%rm3(1) =0.0_dp
     grid%rm1(2:)= upf%r(2:)**(-1)
     grid%rm2(2:)= upf%r(2:)**(-2)
     grid%rm3(2:)= upf%r(2:)**(-3)
  END IF
  !
END SUBROUTINE add_upf_grid
!
!---------------------------------------------------------------------
SUBROUTINE set_upf_q (upf)
  !---------------------------------------------------------------------
  !
  ! For USPP we set the augmentation charge as an l-dependent array in all
  ! cases. This is already the case when upf%tpawp or upf%q_with_l are .true.
  ! For vanderbilt US pseudos, where nqf and rinner are non zero, we do here
  ! what otherwise would be done multiple times in many parts of the code
  ! whenever the q_l(r) were to be constructed. 
  ! For simple rrkj3 pseudos we duplicate the information contained in q(r)
  ! for all q_l(r).
  !
  ! This requires a little extra memory but unifies the treatment of q_l(r)
  ! and allows further weaking with the augmentation charge.
  !
  IMPLICIT NONE
  !
  TYPE (pseudo_upf) :: upf
  !
  !     Local variables
  !
  INTEGER :: nb, mb, ijv, ir, ilast, l, l1, l2
  !
  IF ( upf%tvanp .and. .not.upf%q_with_l ) THEN
     ALLOCATE( upf%qfuncl ( upf%mesh, upf%nbeta*(upf%nbeta+1)/2, 0:upf%nqlc-1 ) )
     upf%qfuncl  = 0.0_DP
     
     DO nb = 1, upf%nbeta
        DO mb = nb, upf%nbeta
           ! ijv is the combined (nb,mb) index
           ijv = mb * (mb-1) / 2 + nb
           l1=upf%lll(nb) ; l2=upf%lll(mb)
           ! copy q(r) to the l-dependent grid 
           DO l=abs(l1-l2),l1+l2,2
              upf%qfuncl(1:upf%mesh,ijv,l) = upf%qfunc(1:upf%mesh,ijv)
           END DO
! adjust the inner values on the l-dependent grid if nqf and rinner are defined
           IF ( upf%nqf > 0 ) THEN
              DO l = abs(l1-l2),l1+l2, 2
                 IF ( upf%rinner (l+1) > 0.0_dp) THEN
                    DO ir = 1, upf%kkbeta
                       if (upf%r(ir) <upf%rinner (l+1) ) ilast = ir
                    END DO
                    CALL setqfnew( upf%nqf,upf%qfcoef(:,l+1,nb,mb), ilast, &
                         upf%r, l, 2, upf%qfuncl(:,ijv,l) )
                 END IF
              END DO
           END IF
        END DO
     END DO
  END IF

END SUBROUTINE set_upf_q
!------------------------------------------------------------------------
SUBROUTINE setqfnew( nqf, qfcoef, mesh, r, l, n, rho )
  !-----------------------------------------------------------------------
  !
  ! ... Computes the Q function from its polynomial expansion (r < rinner)
  ! ... On input: nqf = number of polynomial coefficients
  ! ...    qfcoef(nqf)= the coefficients defining Q
  ! ...          mesh = number of mesh point
  ! ...        r(mesh)= the radial mesh
  ! ...             l = angular momentum
  ! ...             n = additional exponent, result is multiplied by r^n
  ! ... On output:
  ! ...      rho(mesh)= r^n * Q(r)
  !
  USE upf_kinds, ONLY: dp
  !
  IMPLICIT NONE
  !
  INTEGER,  INTENT(in):: nqf, l, mesh, n
  REAL(dp), INTENT(in) :: r(mesh), qfcoef(nqf)
  REAL(dp), INTENT(out) :: rho(mesh)
  !
  INTEGER  :: ir, i
  REAL(dp) :: rr
  !
  DO ir = 1, mesh
     rr = r(ir)**2
     rho(ir) = qfcoef(1)
     DO i = 2, nqf
        rho(ir) = rho(ir) + qfcoef(i)*rr**(i-1)
     ENDDO
     rho(ir) = rho(ir)*r(ir)**(l+n)
  ENDDO
  !
  RETURN
  !
END SUBROUTINE setqfnew
!

!=----------------------------------------------------------------------------=!
  END MODULE upf_to_internal
!=----------------------------------------------------------------------------=!
