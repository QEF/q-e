!
! Copyright (C) 2001-2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!---------------------------------------------------------------------
!    Common variables for the phonon program
!    The variables needed to describe the modes and the small group of q
!
module modes
  use parameters, only : DP
  integer :: irgq(48), nsymq, irotmq, nirr, nmodes
  ! selects the operations of the small group
  ! the number of symmetry of the small group
  ! selects the symmetry sending q <-> -q+G
  ! the number of irreducible representation
  !    contained in the dynamical matrix
  ! the number of modes
  integer, allocatable, target :: npert (:) !3 * nat )
  ! the number of perturbations per IR
  integer :: npertx, invs(48)
  ! max number of perturbations per IR
  ! the inver of each matrix

  real (kind=DP), allocatable    :: rtau (:,:,:) !3, 48, nat)
  ! coordinates of direct translations
  real (kind=DP) :: gi(3,48), gimq(3)
  ! the possible G associated to each symmetry
  ! the G associated to the symmetry q<->-q+G

  integer , parameter :: max_irr_dim = 4    ! maximal allowed dimension for
                                           ! irreducible representattions

  complex (kind=DP), pointer :: u (:,:), & ! 3 * nat, 3 * nat),
       ubar(:),&                           ! 3 * nat), &
       t (:,:,:,:),&                       ! max_irr_dim, max_irr_dim, 48,3 * nat),
       tmq (:,:,:)                         ! max_irr_dim, max_irr_dim, 3 * nat)
  ! the transformation modes patterns
  ! the mode for deltarho
  ! the symmetry in the base of the pattern
  ! the symmetry q<->-q in the base of the pa

  logical :: minus_q
  ! if true there is the symmetry sending q<->-q
end module modes

!
!   The dynamical matrix
!
module dynmat
  use parameters, only : DP
  complex (kind=DP), allocatable :: dyn00 (:,:),& ! 3 * nat, 3 * nat),
       dyn (:,:)                              ! 3 * nat, 3 * nat)
  ! the initial dynamical matrix
  ! the dynamical matrix

  real (kind=DP), allocatable    :: w2 (:)        ! 3 * nat)
  ! omega^2
end module dynmat

!
!   The q point
!
module qpoint
  use parameters, only : DP
  integer, pointer :: igkq (:)     ! npwx)
  ! correspondence k+q+G <-> G
  integer :: nksq, npwq
  ! the real number of k points
  ! the number of plane waves for q

  real (kind=DP) :: xq(3)
  ! the coordinates of the q point

  complex (kind=DP), allocatable :: eigqts (:) ! nat)
  ! the phases associated to the q
end module qpoint
!
!    The wavefunctions at point k+q
!
module eqv
  use parameters, only : DP
  complex (kind=DP), pointer :: evq (:,:)
!
!    The variable describing the linear response problem
!

  complex (kind=DP), allocatable :: dvpsi (:,:), dpsi (:,:)
  ! the product of dV psi
  ! the change of the wavefunctions

  real (kind=DP), allocatable :: dmuxc (:,:,:) ! nrxx, nspin, nspin),
  real (kind=DP), allocatable, target :: vlocq (:,:) ! ngm, ntyp)
  ! the derivative of the xc potential
  ! the local potential at q+G
end module eqv

!
!     the variables for the electric field perturbation
!
module efield
  use parameters, only : DP
  real (kind=DP) :: epsilon (3, 3)
  real (kind=DP), allocatable ::  zstareu (:,:,:),& ! 3, 3, nat),
       zstarue (:,:,:)                          ! 3, nat, 3)
  ! the dielectric constant
  ! the effective charges Z(E,Us) (E=scf,Us=bare)
  ! the effective charges Z(Us,E) (Us=scf,E=bare)

  complex (kind=DP), allocatable :: zstareu0 (:,:),& ! 3, 3 * nat),
       zstarue0 (:,:)                            ! 3 * nat, 3)
  ! the effective charges
  ! the effective charges
end module efield

!
!    The variables needed for non-linear core correction
!
module nlcc_ph
  use parameters, only : DP
  complex (kind=DP), allocatable, target :: drc (:,:) ! ngm, ntyp)
  ! contain the rhoc (without structure fac) for all atomic types
  logical           :: nlcc_any
  ! .t. if any atom-type has nlcc
end module nlcc_ph

!
!    The variables needed for gradient corrected calculations
!
module gc_ph
  use parameters, only : DP
  real (kind=DP), allocatable :: grho (:,:,:),& ! 3, nrxx, nspin),
       dvxc_rr (:,:,:),&           ! nrxx, nspin, nspin), &
       dvxc_sr (:,:,:),&           ! nrxx, nspin, nspin),
       dvxc_ss (:,:,:),&           ! nrxx, nspin, nspin), &
       dvxc_s (:,:,:)              ! nrxx, nspin, nspin)
  ! gradient of the unpert. density
  !
  ! derivatives of the E_xc functiona
  ! r=rho and s=|grad(rho)|
  !
end module gc_ph

!
!   These are additional variables needed for the linear response
!   program with the US pseudopotentials
!
module phus
  use parameters, only : DP
  real (kind=DP), allocatable  :: alphasum (:,:,:,:),&! nhm*(nhm+1)/2,3,nat,nspin)
                                                  ! used to compute modes
       dpqq(:,:,:,:)                              ! dipole moment of each Q
  complex (kind=DP), allocatable :: &
       int1 (:,:,:,:,:),&  ! nhm, nhm, 3, nat, nspin),&
       int2 (:,:,:,:,:),&  ! nhm, nhm, 3,nat, nat),&
       int3 (:,:,:,:,:),&  ! nhm, nhm, 3, nat, nspin),&
       int4 (:,:,:,:,:),&  ! nhm*(nhm+1)/2, 3, 3, nat, nspin),&
       int5 (:,:,:,:,:)    ! nhm*(nhm+1)/2, 3, 3, nat, nat),&
  complex (kind=DP), allocatable, target :: &
       becp1 (:,:,:),   &  ! nkbtot, nbnd, nksq),&
       alphap (:,:,:,:)    ! nkbtot, nbnd, 3, nksq)
  ! integrals of dQ and V_eff
  ! integrals of dQ and V_loc
  ! integrals of Q and dV_Hxc
  ! integrals of d^2Q and V
  ! integrals of dQ and dV_lo
  ! the becq used in ch_psi
  ! the derivative of the bec
end module phus
!
!   the variables needed for partial computation of dynamical matrix
!
module partial
  use parameters, only : DP
  integer, allocatable :: comp_irr (:),& ! 3 * nat ),
       ifat (:),&                    ! nat),
       done_irr (:),&                ! 3 * nat), &
       list (:),&                    ! 3 * nat),
       atomo (:)                     ! nat)
  ! if 1 this representation has to be computed
  ! if 1 this matrix element is computed
  ! if 1 this representation has been done
  ! a list of representations
  ! which atom
  integer ::   nat_todo, nrapp
  ! number of atoms to compute
  ! The representation to do
  logical :: all_comp
  ! if true all representation have been computed
end module partial

!
!   the variable controlling the phonon run
!
module control_ph
  use parameters, only : DP, npk
  integer, parameter :: maxter = 100
  ! maximum number of iterations
  integer :: niter_ph, nmix_ph, nbnd_occ(npk), irr0, iter0, maxirr
  ! maximum number of iterations (read from input)
  ! mixing type
  ! occupated bands in metals
  ! starting representation
  ! starting iteration
  ! maximum number of representation

  real (kind=DP) :: tr2_ph, alpha_mix(maxter), time_now, alpha_pv
  ! convergence threshold
  ! the mixing parameter
  ! CPU time up to now
  ! the alpha value for shifting the bands

  logical :: lgamma, convt, epsil, trans, elph, zue, recover
  ! if true this is a q=0 computation
  ! if true the phonon has converged
  ! if true computes dielec. const and eff. c
  ! if true computes phonons
  ! if true computes electron-phonon interact
  ! if true computes eff.cha. with ph
  ! if true the run restart

end module control_ph

!
! a character common for phonon
!
module char_ph
  character(len=75) :: title_ph! * 75
  ! title of the phonon run
end module char_ph
!
!   the units of the files and the record lengths
!
module units_ph

  integer :: iuwfc, lrwfc, iuvkb, iubar, lrbar, iudwf, iupsir, &
       lrdwf, iudrhous, lrdrhous, iudyn, iupdyn, iunrec, iudvscf, iudrho, &
       lrdrho
  ! iunit with the wavefunctions
  ! the length of wavefunction record
  ! unit with vkb
  ! unit with the part DV_{bare}
  ! length of the DV_{bare}
  ! unit with D psi
  ! unit with evc in real space
  ! length of D psi record
  ! the unit with the products
  ! the lenght of the products
  ! the unit for the dynamical matrix
  ! the unit for the partial dynamical matrix
  ! the unit with the recover data
  ! the unit where the delta Vscf is written
  ! the unit where the delta rho is written
  ! the length of the deltarho files

end module units_ph
!
!    the name of the files
!
module output
  character (len=14) :: fildyn, filelph, fildvscf, fildrho
  ! output file for the dynamical matrix
  ! output file for electron-phonon coefficie
  ! output file for deltavscf
  ! output file for deltarho
end module output

module phcom
  use  modes
  use dynmat
  use qpoint
  use eqv
  use efield
  use nlcc_ph
  use gc_ph
  use phus
  use partial
  use control_ph
  use char_ph
  use units_ph
  use output
end module phcom
