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

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: set_pseudo_upf
  SAVE

!=----------------------------------------------------------------------------=!
  CONTAINS
!=----------------------------------------------------------------------------=!
!
!---------------------------------------------------------------------
subroutine set_pseudo_upf (is, upf, grid)
  !---------------------------------------------------------------------
  !
  !   set "is"-th pseudopotential using the Unified Pseudopotential Format
  !   "upf" - convert and copy to internal variables
  !   If "grid" is present, reconstruct radial grid.
  !   Obsolescent - for old-style PP formats only.
  !
  USE funct, ONLY: set_dft_from_name, set_dft_from_indices
  !
  USE pseudo_types
  USE radial_grids, ONLY: radial_grid_type, allocate_radial_grid
  !
  implicit none
  !
  INTEGER :: is
  TYPE (pseudo_upf) :: upf
  TYPE (radial_grid_type), target, optional :: grid
  !
  !     Local variables
  !
  integer :: iexch,icorr,igcx,igcc
  INTEGER :: nb, mb, ijv, ir, ilast, l, l1, l2
  !
  ! old formats never contain "1/r" pseudopotentials
  !
  upf%tcoulombp = .false.
  !
  ! workaround for rrkj format - it contains the indices, not the name
  !
  if ( upf%dft(1:6)=='INDEX:') then
     read( upf%dft(7:10), '(4i1)') iexch,icorr,igcx,igcc
     call set_dft_from_indices(iexch,icorr,igcx,igcc, 0) !Cannot read nonloc in this format
  else
     call set_dft_from_name( upf%dft )
  end if
  !
  if(present(grid)) then
    call allocate_radial_grid(grid,upf%mesh)
    grid%dx   = upf%dx
    grid%xmin = upf%xmin
    grid%zmesh= upf%zmesh
    grid%mesh = upf%mesh
    !
    grid%r  (1:upf%mesh) = upf%r  (1:upf%mesh)
    grid%rab(1:upf%mesh) = upf%rab(1:upf%mesh)
    upf%grid => grid
  endif
  !
  ! For USPP we set the augmentation charge as an l-dependent array in all cases.
  ! This is already the case when upf%tpawp or upf%q_with_l are .true. .
  ! For vanderbilt US pseudos, where nqf and rinner are non zero, we do here what otherwise
  ! would be done multiple times in many parts of the code (such as in init_us_1, addusforce_r, 
  ! bp_calc_btq, compute_qdipol) whenever the q_l(r) were to be constructed. 
  ! For simple rrkj3 pseudos we duplicate the infomation contained in q(r) for all q_l(r).
  !
  ! This requires a little extra memory but unifies the treatment of q_l(r) and allows further 
  ! tweaking with the augmentation charge.
  !
  if ( upf%tvanp .and. .not.upf%q_with_l ) then
     ALLOCATE( upf%qfuncl ( upf%mesh, upf%nbeta*(upf%nbeta+1)/2, 0:upf%nqlc-1 ) )
     upf%qfuncl  = 0.0_DP

     do nb = 1, upf%nbeta
        do mb = nb, upf%nbeta
           ! ijv is the combined (nb,mb) index
           ijv = mb * (mb-1) / 2 + nb
           l1=upf%lll(nb) ; l2=upf%lll(mb)
! copy q(r) to the l-dependent grid 
           DO l=abs(l1-l2),l1+l2,2
              upf%qfuncl(1:upf%mesh,ijv,l) = upf%qfunc(1:upf%mesh,ijv)
           END DO
! adjust the inner values on the l-dependent grid if nqf and rinner are defined
           if ( upf%nqf > 0 ) then
              do l = abs(l1-l2),l1+l2, 2
                 if ( upf%rinner (l+1) > 0.0_dp) then
                    do ir = 1, upf%kkbeta
                       if (upf%r(ir) <upf%rinner (l+1) ) ilast = ir
                    enddo
                    call setqfnew( upf%nqf,upf%qfcoef(1,l+1,nb,mb), ilast, upf%r, l, 2, upf%qfuncl(1,ijv,l) )
                 end if
              end do
           end if
        enddo
     enddo
  end if

end subroutine set_pseudo_upf

!=----------------------------------------------------------------------------=!
  END MODULE upf_to_internal
!=----------------------------------------------------------------------------=!
