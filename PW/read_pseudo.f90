!
!---------------------------------------------------------------------
subroutine read_pseudo (is, iunps, ierr)
  !---------------------------------------------------------------------
  !
  !   read "is"-th pseudopotential in the Unified Pseudopotential Format
  !   from unit "iunps" - convert and copy to internal PWscf variables
  !   return error code in "ierr" (success: ierr=0)
  !
  ! PWSCF modules
  !
  use atom,  only: zmesh, mesh,dx, r, rab, vloc_at, chi, oc, nchi, lchi, &
       rho_at, rho_atc
  use char,  only: psd
  use pseud, only: zp, lmax, lloc
  use nl_c_c,only: nlcc
  use us,    only: dion, betar, qqq, qfcoef, qfunc, nqf, nqlc, rinner, &
       nh, nbeta, kkbeta, lll, tvanp
  use funct, only: dft, which_dft
  !
  use pseudo_types
  use read_pseudo_module
  !
  implicit none
  !
  integer :: is, iunps, ierr
  !
  !     Local variables
  !
  integer :: nb
  TYPE (pseudo_upf) :: upf
  !
  !
  call read_pseudo_upf(iunps, upf, ierr)
  !
  if (ierr .ne. 0) return
  !
  zp(is)  = upf%zp
  psd (is)= upf%psd
  tvanp(is)=upf%tvanp
  nlcc(is) = upf%nlcc
  dft = upf%dft
  call which_dft (upf%dft)
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
  !!!
  vloc_at(1:upf%mesh,is) = upf%vloc(1:upf%mesh)
  CALL deallocate_pseudo_upf( upf )

end subroutine read_pseudo
