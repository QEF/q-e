
! This module is used, for the time being, as an interface
! between the UPF pseudo type and the pseudo variables internal representation

!=----------------------------------------------------------------------------=!
  MODULE upf_to_internal
!=----------------------------------------------------------------------------=!

  IMPLICIT NONE
  SAVE

!=----------------------------------------------------------------------------=!
  CONTAINS
!=----------------------------------------------------------------------------=!

!
!---------------------------------------------------------------------
subroutine set_pseudo (is, upf, ierr)
  !---------------------------------------------------------------------
  !
  !   set "is"-th pseudopotential using the Unified Pseudopotential Format
  !   dummy argument ( upf ) - convert and copy to internal PWscf variables
  !   return error code in "ierr" (success: ierr=0)
  !
  ! PWSCF modules
  !
  use pwcom
  use funct
  !
  use pseudo_types
  use read_pseudo_module
  !
  implicit none
  !
  real(kind=DP), parameter :: rcut = 10.d0
  integer :: is, ierr, ir
  !
  !     Local variables
  !
  integer :: nb
  TYPE (pseudo_upf) :: upf
  !
  !
  if (ierr .ne. 0) return
  !
  zp(is)  = upf%zp
  psd (is)= upf%psd
  tvanp(is)=upf%tvanp
  nlcc(is) = upf%nlcc
  dft = upf%dft
  call which_dft (upf%dft, iexch, icorr, igcx, igcc)
  mesh(is) = upf%mesh
  !
  nchi(is) = upf%nwfc
  lchi(1:upf%nwfc, is) = upf%lchi(1:upf%nwfc)
  oc(1:upf%nwfc, is) = upf%oc(1:upf%nwfc)
  chi(1:upf%mesh, 1:upf%nwfc, is) = upf%chi(1:upf%mesh, 1:upf%nwfc)
  !
  nbeta(is)= upf%nbeta
  kkbeta(is)=0
  do nb=1,upf%nbeta
     kkbeta(is)=max(upf%kkbeta(nb),kkbeta(is))
  end do
  betar(1:upf%mesh, 1:upf%nbeta, is) = upf%beta(1:upf%mesh, 1:upf%nbeta)
  dion(1:upf%nbeta, 1:upf%nbeta, is) = upf%dion(1:upf%nbeta, 1:upf%nbeta)
  !
  lmax(is) = upf%lmax
  nqlc(is) = upf%nqlc
  nqf (is) = upf%nqf
  lll(1:upf%nbeta,is) = upf%lll(1:upf%nbeta)
  rinner(1:upf%nqlc,is) = upf%rinner(1:upf%nqlc)
  qqq(1:upf%nbeta,1:upf%nbeta,is) = upf%qqq(1:upf%nbeta,1:upf%nbeta)
  qfunc (1:upf%mesh, 1:upf%nbeta, 1:upf%nbeta, is) = &
       upf%qfunc(1:upf%mesh,1:upf%nbeta,1:upf%nbeta)
  qfcoef(1:upf%nqf, 1:upf%nqlc, 1:upf%nbeta, 1:upf%nbeta, is ) = &
       upf%qfcoef( 1:upf%nqf, 1:upf%nqlc, 1:upf%nbeta, 1:upf%nbeta )
  !
  r  (1:upf%mesh, is) = upf%r  (1:upf%mesh)
  rab(1:upf%mesh, is) = upf%rab(1:upf%mesh)
  !
  if ( upf%nlcc) then
     rho_atc(1:upf%mesh, is) = upf%rho_atc(1:upf%mesh)
  else
     rho_atc(:,is) = 0.d0
  end if
  rho_at (1:upf%mesh, is) = upf%rho_at (1:upf%mesh)
  !!! TEMP
  lloc(is) = 0
  vnl(1:upf%mesh,lloc(is),is) = upf%vloc(1:upf%mesh)
  !!!

  do ir = 1, mesh (is)
    if (r (ir, is) .gt.rcut) then
        msh (is) = ir
        goto 5
    endif
  enddo
  msh (is) = mesh (is)
  !
  ! force msh to be odd for simpson integration
  !
5 msh (is) = 2 * ( (msh (is) + 1) / 2) - 1

  zv(is) = zp(is)

end subroutine set_pseudo

!=----------------------------------------------------------------------------=!
  END MODULE upf_to_internal
!=----------------------------------------------------------------------------=!
