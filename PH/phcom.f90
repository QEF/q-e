!
! Copyright (C) 2001-2004 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
!
! ... Common variables for the phonon program
!  
MODULE modes
  USE kinds,  ONLY : DP
  !
  ! ... The variables needed to describe the modes and the small group of q
  !
  SAVE
  !
  INTEGER :: irgq(48), nsymq, irotmq, nirr, nmodes
  ! selects the operations of the small group
  ! the number of symmetry of the small group
  ! selects the symmetry sending q <-> -q+G
  ! the number of irreducible representation
  !    contained in the dynamical matrix
  ! the number of modes
  INTEGER, ALLOCATABLE, TARGET :: npert(:) !3 * nat )
  ! the number of perturbations per IR
  INTEGER :: npertx, invs(48)
  ! max number of perturbations per IR
  ! the inver of each matrix
  REAL (DP), ALLOCATABLE :: rtau(:,:,:) !3, 48, nat)
  ! coordinates of direct translations
  REAL (DP) :: gi(3,48), gimq(3)
  ! the possible G associated to each symmetry
  ! the G associated to the symmetry q<->-q+G
  INTEGER, PARAMETER :: max_irr_dim = 4    ! maximal allowed dimension for
                                           ! irreducible representattions
  COMPLEX (DP), POINTER :: &
       u(:,:),                     &!  3 * nat, 3 * nat),
       ubar(:),                    &!  3 * nat), &
       t(:,:,:,:),                 &! max_irr_dim, max_irr_dim, 48,3 * nat),
       tmq(:,:,:)                   ! max_irr_dim, max_irr_dim, 3 * nat)
  ! the transformation modes patterns
  ! the mode for deltarho
  ! the symmetry in the base of the pattern
  ! the symmetry q<->-q in the base of the pa
  LOGICAL :: &
       minus_q       !  if .TRUE. there is the symmetry sending q<->-q
  !     
END MODULE modes
!
!
MODULE dynmat
  USE kinds, ONLY :  DP
  !
  ! ... The dynamical matrix 
  !
  SAVE
  !
  COMPLEX (DP), ALLOCATABLE :: &
       dyn00(:,:),           &! 3 * nat, 3 * nat),
       dyn(:,:)               ! 3 * nat, 3 * nat)
  ! the initial dynamical matrix
  ! the dynamical matrix
  REAL (DP), ALLOCATABLE :: &
       w2(:)                  ! 3 * nat)
  ! omega^2
  !
END MODULE dynmat
!
!
MODULE qpoint
  USE kinds, ONLY :  DP
  !
  ! ... The q point
  !
  SAVE
  !
  INTEGER, POINTER :: igkq(:)     ! npwx)
  ! correspondence k+q+G <-> G
  INTEGER :: nksq, npwq
  ! the real number of k points
  ! the number of plane waves for q
  REAL (DP) :: xq(3)
  ! the coordinates of the q point
  COMPLEX (DP), ALLOCATABLE :: eigqts(:) ! nat)
  ! the phases associated to the q
  !
END MODULE qpoint
!
!
MODULE eqv
  USE kinds, ONLY :  DP
  !
  ! ... The wavefunctions at point k+q 
  !
  SAVE
  !
  COMPLEX (DP), POINTER :: evq(:,:)
  !
  ! ... The variable describing the linear response problem 
  !
  COMPLEX (DP), ALLOCATABLE :: dvpsi(:,:), dpsi(:,:)
  ! the product of dV psi
  ! the change of the wavefunctions
  REAL (DP), ALLOCATABLE :: dmuxc(:,:,:)        ! nrxx, nspin, nspin),
  REAL (DP), ALLOCATABLE, TARGET :: vlocq(:,:)  ! ngm, ntyp)
  ! the derivative of the xc potential
  ! the local potential at q+G
  !
END MODULE eqv
!
!
MODULE efield_mod
  USE kinds, ONLY :  DP
  !
  ! ... the variables for the electric field perturbation
  !  
  SAVE
  !
  REAL (DP) :: epsilon (3, 3)
  REAL (DP), ALLOCATABLE :: &
       zstareu(:,:,:),       &! 3, 3, nat),
       zstarue(:,:,:)         ! 3, nat, 3)
  ! the dielectric constant
  ! the effective charges Z(E,Us) (E=scf,Us=bare)
  ! the effective charges Z(Us,E) (Us=scf,E=bare)
  COMPLEX (DP), ALLOCATABLE :: &
       zstareu0(:,:),        &! 3, 3 * nat),
       zstarue0(:,:)          ! 3 * nat, 3)
  ! the effective charges
  !
END MODULE efield_mod
!
!
MODULE nlcc_ph
  USE kinds, ONLY :  DP
  !
  ! ... The variables needed for non-linear core correction
  !
  SAVE
  !
  COMPLEX (DP), ALLOCATABLE, TARGET :: drc(:,:) ! ngm, ntyp)
  ! contain the rhoc (without structure fac) for all atomic types
  LOGICAL :: nlcc_any
  ! .T. if any atom-type has nlcc
  !
END MODULE nlcc_ph
!
!
MODULE gc_ph
  USE kinds, ONLY :  DP
  !
  ! ... The variables needed for gradient corrected calculations
  !
  SAVE
  !
  REAL (DP), ALLOCATABLE :: &
       grho(:,:,:),              &! 3, nrxx, nspin),
       gmag(:,:,:),              &! 3, nrxx, nspin),
       vsgga(:),                 &! nrxx
       segni(:),                 &! nrxx
       dvxc_rr(:,:,:),           &! nrxx, nspin, nspin), &
       dvxc_sr(:,:,:),           &! nrxx, nspin, nspin),
       dvxc_ss(:,:,:),           &! nrxx, nspin, nspin), &
       dvxc_s(:,:,:)              ! nrxx, nspin, nspin)
  !
  ! in the noncollinear case gmag contains the gradient of the magnetization
  ! grho the gradient of rho+ and of rho-, the eigenvalues of the spin density
  ! vsgga= 0.5* (V_up-V_down) to be used in the calculation of the change
  ! of the exchange and correlation magnetic field.
  ! gradient of the unpert. density
  !
  ! derivatives of the E_xc functiona
  ! r=rho and s=|grad(rho)|
  !
END MODULE gc_ph
!
!
MODULE phus
  USE kinds, ONLY :  DP
  !
  ! ... These are additional variables needed for the linear response
  ! ... program with the US pseudopotentials
  !
  SAVE
  !
  REAL (DP), ALLOCATABLE :: &
       alphasum(:,:,:,:),   &! nhm*(nhm+1)/2,3,nat,nspin)
                             ! used to compute modes
       dpqq(:,:,:,:)         ! dipole moment of each Q
  COMPLEX (DP), ALLOCATABLE :: &
       int1(:,:,:,:,:),     &! nhm, nhm, 3, nat, nspin),&
       int2(:,:,:,:,:),     &! nhm, nhm, 3,nat, nat),&
       int3(:,:,:,:,:),     &! nhm, nhm, max_irr_dim, nat, nspin),&
       int4(:,:,:,:,:),     &! nhm*(nhm+1)/2, 3, 3, nat, nspin),&
       int5(:,:,:,:,:),     &! nhm*(nhm+1)/2, 3, 3, nat, nat),&
       int1_nc(:,:,:,:,:),     &! nhm, nhm, 3, nat, nspin),&
       int2_so(:,:,:,:,:,:),   &! nhm, nhm, 3,nat,nat,nspin),&
       int3_nc(:,:,:,:,:),     &! nhm, nhm, max_irr_dim, nat, nspin),&
       int4_nc(:,:,:,:,:,:),   &! nhm, nhm, 3, 3, nat, nspin),&
       int5_so(:,:,:,:,:,:,:), &! nhm*(nhm+1)/2, 3, 3, nat, nat, nspin),&
       becsum_nc(:,:,:,:),     &! nhm*(nhm+1)/2,nat,npol,npol)
       alphasum_nc(:,:,:,:,:), &! nhm*(nhm+1)/2,3,nat,npol,npol)
       dpqq_so(:,:,:,:,:)       ! dipole moment of each Q and the fcoef factors

  COMPLEX (DP), ALLOCATABLE, TARGET :: &
       becp1(:,:,:),        &! nkbtot, nbnd, nksq),&
       becp1_nc(:,:,:,:),   &! nkbtot, npol, nbnd, nksq),&
       alphap(:,:,:,:),     &! nkbtot, nbnd, 3, nksq)
       alphap_nc(:,:,:,:,:)  ! nkbtot, npol, nbnd, 3, nksq)
  ! integrals of dQ and V_eff
  ! integrals of dQ and V_loc
  ! integrals of Q and dV_Hxc
  ! integrals of d^2Q and V
  ! integrals of dQ and dV_lo
  ! the becq used in ch_psi
  ! the derivative of the bec
  !
END MODULE phus
!
!
MODULE partial
  USE kinds, ONLY :  DP
  !
  ! ... the variables needed for partial computation of dynamical matrix
  !
  SAVE
  !  
  INTEGER, ALLOCATABLE :: &
       comp_irr(:),           &! 3 * nat ),
       ifat(:),               &! nat),
       done_irr(:),           &! 3 * nat), &
       list(:),               &! 3 * nat),
       atomo(:)                ! nat)
  ! if 1 this representation has to be computed
  ! if 1 this matrix element is computed
  ! if 1 this representation has been done
  ! a list of representations
  ! which atom
  INTEGER :: nat_todo, nrapp
  ! number of atoms to compute
  ! The representation to do
  LOGICAL :: all_comp
  ! if .TRUE. all representation have been computed
  !
END MODULE partial
!
MODULE gamma_gamma
  INTEGER, ALLOCATABLE :: &
           has_equivalent(:),  &  ! 0 if the atom has to be calculated
           n_equiv_atoms(:),   &  ! number of equivalent atoms
           equiv_atoms(:,:)       ! which atoms are equivalent

  INTEGER :: n_diff_sites,    &   ! Number of different sites
             nasr                 ! atom calculated with asr
                                  !
  LOGICAL :: asr                  ! if true apply the asr

END MODULE gamma_gamma
!
MODULE control_ph
  USE kinds, ONLY :  DP
  USE parameters, ONLY: npk
  !
  ! ... the variable controlling the phonon run
  !
  SAVE
  !
  INTEGER, PARAMETER :: maxter = 100
  ! maximum number of iterations
  INTEGER :: niter_ph, nmix_ph, nbnd_occ(npk), irr0, maxirr
  ! maximum number of iterations (read from input)
  ! mixing type
  ! occupated bands in metals
  ! starting representation
  ! maximum number of representation
  real(DP) :: tr2_ph
  ! threshold for phonon calculation
  REAL (DP) :: alpha_mix(maxter), time_now, alpha_pv
  ! the mixing parameter
  ! CPU time up to now
  ! the alpha value for shifting the bands
  LOGICAL :: lgamma,      &! if .TRUE. this is a q=0 computation
             lgamma_gamma,&! if .TRUE. this is a q=0 computation with k=0 only 
             convt,       &! if .TRUE. the phonon has converged
             epsil,       &! if .TRUE. computes dielec. const and eff. charges
             trans,       &! if .TRUE. computes phonons
             elph,        &! if .TRUE. computes electron-ph interaction coeffs
             zue,         &! if .TRUE. computes eff. charges as induced polarization
             recover,     &! if .TRUE. the run restarts
             lrpa,        &! if .TRUE. calculates the RPA dielectric constant
             lnoloc,      &! if .TRUE. calculates the dielectric constant
                           ! neglecting local field effects
             search_sym,  &! if .TRUE. search the mode symmetry
             lnscf,       &! if .TRUE. the run makes first a nscf calculation
             ldisp,       &! if .TRUE. the run calculates full phonon dispersion
             reduce_io     ! if .TRUE. reduces needed I/O
  !
END MODULE control_ph
!
!
MODULE freq_ph
  !
  USE kinds,   ONLY : DP
  !
  SAVE
  !
  ! ... the variables for computing frequency dependent dielectric constant
  !
  LOGICAL :: fpol ! if .TRUE. dynamic dielectric constant is computed
  !
  INTEGER, PARAMETER :: nfsmax=50  ! # of maximum frequencies
  INTEGER :: nfs                   ! # of frequencies
  !
  REAL (KIND=DP) :: fiu(nfsmax)    ! values  of frequency
  !
END MODULE freq_ph
!
!
MODULE units_ph
  !
  ! ... the units of the files and the record lengths
  !
  SAVE
  !
  INTEGER :: &
       iuwfc, lrwfc, iuvkb, iubar, lrbar, iuebar, lrebar, iudwf, iupsir, &
       lrdwf, iudrhous, lrdrhous, iudyn, iupdyn, iunrec, iudvscf, iudrho, &
       lrdrho, iucom, lrcom, iudvkb3, lrdvkb3
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
  ! the unit of the bare commutator in US case
  ! the length  of the bare commutator in US case
  logical, ALLOCATABLE :: this_dvkb3_is_on_file(:), &
                          this_pcxpsi_is_on_file(:,:)
  !
END MODULE units_ph
!
!
MODULE output
  !
  ! ... the name of the files
  !
  SAVE
  !
  CHARACTER (LEN=256) :: fildyn, fildvscf, fildrho
  ! output file for the dynamical matrix
  ! output file for deltavscf
  ! output file for deltarho
  !
END MODULE output
!
!
MODULE disp
  !
  USE kinds, ONLY: DP
  !
  SAVE
  !
  INTEGER, PARAMETER :: nqmax = 1000
  !
  INTEGER :: nq1, nq2, nq3
    ! number of q-points in each direction
  INTEGER :: iq1, iq2, iq3
    ! specific q point from the regular grid
    ! (i.e., iq1/nq1,iq2/nq2,iq3/nq3)
  INTEGER :: nqs
    ! number of q points to be calculated 
  REAL (DP), ALLOCATABLE :: x_q(:,:)
    ! coordinates of the q points
  !
END MODULE disp
!
!
MODULE phcom
  USE modes
  USE dynmat
  USE qpoint
  USE eqv
  USE efield_mod
  USE nlcc_ph
  USE gc_ph
  USE phus
  USE partial
  USE control_ph
  USE freq_ph
  USE units_ph
  USE output
  USE gamma_gamma
  USE disp 
END MODULE phcom
