!
! Copyright (C) 2004 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

!
!---------------------------------------------------------------------
subroutine read_pseudoupf 
  !---------------------------------------------------------------------
  !
  !   read "is"-th pseudopotential in the Unified Pseudopotential Format
  !   from unit "iunps" - convert and copy to internal PWscf variables
  !   return error code in "ierr" (success: ierr=0)
  !
  ! PWSCF modules
  !
  use ld1inc
  use funct, only: set_dft_from_name
  !
  use pseudo_types
  use read_upf_module
  !
  implicit none
  !
  integer :: iunps, ierr
  !
  !     Local variables
  !
  integer :: nb, ios
  real(DP) :: fpi
  TYPE (pseudo_upf) :: upf
  !
  !
  iunps=2
  open(unit=iunps,file=file_pseudo,status='old',form='formatted', &
       err=100, iostat=ios)
100   call errore('read_pseudoupf','open error on file '//file_pseudo,ios)

  call read_pseudo_upf(iunps, upf, ierr)
  !
  if (ierr .ne. 0) &
     call errore('read_pseudoupf','reading pseudo upf',abs(ierr))
  !
  zval  = upf%zp
  nlcc = upf%nlcc
  call set_dft_from_name (upf%dft)

  if (upf%typ.eq.'NC') then
     pseudotype=2
  else
     pseudotype=3
  endif
  etots=upf%etotps
  lmax = upf%lmax
  mesh = upf%mesh
  r  (1:mesh) = upf%r  (1:upf%mesh)
  rab(1:mesh) = upf%rab(1:upf%mesh)
  r2 (1:mesh) = r(1:mesh)**2
  sqr(1:mesh) = sqrt(r(1:mesh))
  if (.not.upf%has_so) then
     if (r(1) > 0.0_dp) then
        !
        ! r(i+1) = exp(xmin)/zmesh * exp(i*dx)
        !
        dx=log(r(2)/r(1))
        rmax=r(mesh)
        xmin=log(zed*r(1))
        zmesh=zed
     else
        !
        ! r(i+1) = exp(xmin)/zmesh * ( exp(i*dx) - 1 )
        !
        dx=log( (r(3)-r(2)) / r(2) )
        rmax=r(mesh)
        zmesh=zed
        xmin=log(zed*r(2)**2/(r(3)-2.0_dp*r(2)))
     end if
  else
     dx=upf%dx
     xmin=upf%xmin
     zmesh=upf%zmesh
     rmax=exp(xmin+(mesh-1)*dx)/zmesh
  endif
  if (abs(exp(xmin+(mesh-1)*dx)/zed-rmax).gt.1.e-6_dp) &
       call errore('read_pseudoup','mesh not supported',1)

  nwfs = upf%nwfc

  nbeta= upf%nbeta
  lls(1:nbeta)=upf%lll(1:nbeta)

  if (upf%has_so) then
     jjs(1:nbeta)=upf%jjj(1:nbeta)
  else
     jjs=0.0_dp
  endif
  !
  !
  do nb=1,nbeta
     ikk(nb)=upf%kkbeta(nb)
     els(nb)=upf%els_beta(nb)
     rcut(nb)=upf%rcut(nb)
     rcutus(nb)=upf%rcutus(nb)
  end do
  betas(1:mesh, 1:nbeta) = upf%beta(1:upf%mesh, 1:upf%nbeta)
  bmat(1:nbeta, 1:nbeta) = upf%dion(1:upf%nbeta, 1:upf%nbeta)
  !
  if (pseudotype.eq.3) then
     qq(1:nbeta,1:nbeta) = upf%qqq(1:upf%nbeta,1:upf%nbeta)
     qvan (1:mesh, 1:nbeta, 1:nbeta) = &
          upf%qfunc(1:upf%mesh,1:upf%nbeta,1:upf%nbeta)
  else
     qq=0.0_dp
     qvan=0.0_dp
  endif
  !
  !
  if (upf%nlcc) then
     fpi=16.0_dp*atan(1.0_dp)
     rhoc(1:mesh) = upf%rho_atc(1:upf%mesh)*fpi*r2(1:upf%mesh)
  else
     rhoc(:) = 0.0_dp
  end if
  rhos=0.0_dp
  rhos (1:mesh,1) = upf%rho_at (1:upf%mesh)
  phis(1:mesh,1:nwfs)=upf%chi(1:mesh,1:nwfs)
  !!! TEMP
  lloc = -1
  vpsloc(1:mesh) = upf%vloc(1:upf%mesh)
  !!!

  CALL deallocate_pseudo_upf( upf )

end subroutine read_pseudoupf
