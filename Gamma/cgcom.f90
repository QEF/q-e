!
! Copyright (C) 2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
! These are PHONON-specific modules (Conjugate Gradient version)
!
module phunits
  character(len=20) :: fildyn
  character(len=75) :: title_ph
  integer iuwfc, iubar, iudwf, iuscf, iuvkb, lrwfc, lrbar, lrdwf, lrscf
end module phunits

module flags
  logical :: trans, epsil, raman, equil, nlcc_any, asr
end module flags

module dielectric
  use parameters, only: DP
  real(kind=DP) :: epsilon0(3,3)
  real(kind=DP), allocatable :: zstar(:,:,:)
end module dielectric
module modes1
  use parameters, only: DP
  integer :: nmodes
  real(kind=DP), allocatable::  dyn(:,:), u(:,:)
end module modes1

module cgconv
  use parameters, only: DP
  integer :: niter_ph
  real(kind=DP) :: tr2_ph
end module cgconv

module AA
  use parameters, only: DP
  complex(kind=DP), allocatable, target :: aux2(:), aux3(:)
  real(kind=DP), allocatable, target :: auxr(:)
end module AA

module dmu
  use parameters, only: DP
  real(kind=DP), allocatable:: &
       dmuxc(:),   &! d V_xc / d rho
       grho(:,:,:), &! gradient of the unperturbed density
       dvxc_rr(:,:,:), &!
       dvxc_sr(:,:,:), &! derivatives of the E_xc functional w.r.t.
       dvxc_ss(:,:,:), &! r=rho and s=|grad(rho)|
       dvxc_s (:,:,:)
end module dmu

module phon
  use parameters, only: DP
  complex(kind=DP), allocatable:: dvpsi(:,:), dpsi(:,:)
end module phon

module symmetry
  integer :: n_diff_sites, nasr
  integer, allocatable::  equiv_atoms(:,:), n_equiv_atoms(:)
  integer, allocatable:: has_equivalent(:)
end module symmetry

module diffs
  use parameters, only: DP
  integer :: nderiv, first, last
  real(kind=DP) :: deltatau
end module diffs

module cgcom
  use cgconv
  use phunits
  use flags
  use modes1
  use AA
  use phon
  use diffs
  use dmu
  use symmetry
  use dielectric
end module cgcom
