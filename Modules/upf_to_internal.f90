!
! Copyright (C) 2004-2007 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! This module is USEd, for the time being, as an interface
! between the UPF pseudo type and the pseudo variables internal representation

!=----------------------------------------------------------------------------=!
  MODULE upf_to_internal
!=----------------------------------------------------------------------------=!

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
  !   Complete pseudopotential "upf" read from old-style PP files
  !   by reconstructing the radial grid and the Q(r) functions 
  !   Obsolescent, to be used with old formats only
  !
  USE radial_grids, ONLY: radial_grid_type, allocate_radial_grid
  !
  IMPLICIT NONE
  !
  TYPE (pseudo_upf) :: upf
  TYPE (radial_grid_type), target :: grid
  !
  CALL allocate_radial_grid(grid,upf%mesh)
  grid%dx   = upf%dx
  grid%xmin = upf%xmin
  grid%zmesh= upf%zmesh
  grid%mesh = upf%mesh
  !
  grid%r  (1:upf%mesh) = upf%r  (1:upf%mesh)
  grid%rab(1:upf%mesh) = upf%rab(1:upf%mesh)
  upf%grid => grid
  !
  CALL set_upf_q (upf)
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
  ! (such as in init_us_1, addusforce_r, bp_calc_btq, compute_qdipol)
  ! whenever the q_l(r) were to be constructed. 
  ! For simple rrkj3 pseudos we duplicate the infomration contained in q(r)
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
                    CALL setqfnew( upf%nqf,upf%qfcoef(1,l+1,nb,mb), ilast, upf%r, l, 2, upf%qfuncl(1,ijv,l) )
                 END IF
              END DO
           END IF
        END DO
     END DO
  END IF

END SUBROUTINE set_upf_q
!=----------------------------------------------------------------------------=!
  END MODULE upf_to_internal
!=----------------------------------------------------------------------------=!
