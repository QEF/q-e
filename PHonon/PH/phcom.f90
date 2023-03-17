!
! Copyright (C) 2001-2018 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
!
MODULE modes
  !
  !! Common variables for the phonon program needed to describe the modes 
  !! and the small group of q.
  !
  USE kinds,  ONLY : DP
  !
  SAVE
  !
  INTEGER :: nirr
  !! number of irreducible representations contained in the dynamical matrix
  INTEGER :: nmodes
  !! number of modes
  INTEGER, ALLOCATABLE, TARGET :: npert(:) !(3 * nat )
  !! the number of perturbations per IR
  INTEGER :: npertx
  !! max number of perturbations per IR
  COMPLEX (DP), ALLOCATABLE :: u(:,:)     !(3 * nat, 3 * nat)
  !! the transformation modes patterns
  COMPLEX (DP), ALLOCATABLE :: t(:,:,:,:) !(npertx, npertx, 48,3 * nat),
  !! the mode for deltarho
  COMPLEX (DP), ALLOCATABLE :: tmq(:,:,:) !(npertx, npertx, 3 * nat)
  !! the symmetry in the base of the pattern
  ! the symmetry q<->-q in the base of the pa
  !
  CHARACTER(15), ALLOCATABLE :: name_rap_mode(:)
  !! symmetry type of each mode
  INTEGER, ALLOCATABLE :: num_rap_mode(:)
  !! number of the representation for each mode
  !
END MODULE modes
!
MODULE cryst_ph
   !  
   !! This module contains the variables that describe properties of the
   !! crystal that are needed by the phonon program and are not in pw data
   !! probably these variables should be in the common of pw.
   !! These variables are sets immediately after reading the pw variables.
   !
   USE kinds,  ONLY : DP
   !
   SAVE
   !
   LOGICAL :: magnetic_sym
   !! TRUE in the non-collinear magnetic case
   !
END MODULE cryst_ph
!
MODULE dynmat
  !
  !! The dynamical matrix
  !
  USE kinds, ONLY :  DP
  !
  SAVE
  !
  COMPLEX(DP), ALLOCATABLE :: dyn00(:,:)    ! (3 * nat, 3 * nat)
  !! the initial dynamical matrix
  COMPLEX(DP), ALLOCATABLE :: dyn(:,:)      ! (3 * nat, 3 * nat)
  !! the dynamical matrix
  COMPLEX(DP), ALLOCATABLE :: dyn_rec(:,:)  ! (3 * nat, 3 * nat)
  !! the contribution of each representation to the dynamical matrix
  REAL (DP), ALLOCATABLE :: w2(:)           ! (3 * nat)
  !! omega^2
  !
  ! DFPT+U
  COMPLEX(DP), ALLOCATABLE :: dyn_hub_bare(:,:) ! (3*nat,*3nat)
  !! the bare part of the Hubbard dynamical matrix
  COMPLEX(DP), ALLOCATABLE :: dyn_hub_scf(:,:)  ! (3*nat,*3nat)
  !! the scf part of the  Hubbard dynamical matrix
  !
END MODULE dynmat
!
!
MODULE efield_mod
  !
  !! The variables for the electric field perturbation.
  !
  USE kinds, ONLY :  DP
  !
  SAVE
  !
  REAL (DP) :: epsilon (3, 3)
  !! the dielectric constant
  REAL (DP), ALLOCATABLE :: zstareu(:,:,:)     !(3, 3, nat),
  !! the effective charges Z(E,Us) (E=scf,Us=bare)
  REAL (DP), ALLOCATABLE :: zstarue(:,:,:)     !(3, nat, 3)
  !! the effective charges Z(Us,E) (Us=scf,E=bare)
  COMPLEX (DP), ALLOCATABLE :: &
       zstareu0(:,:),        &! 3, 3 * nat),
       zstarue0(:,:),        &! 3 * nat, 3)
       zstarue0_rec(:,:)      ! 3 * nat, 3)
  ! the effective charges
  !
END MODULE efield_mod
!
!
MODULE nlcc_ph
  !
  !! The variables needed for non-linear core correction.
  !
  USE kinds, ONLY :  DP
  !
  SAVE
  !
  COMPLEX (DP), ALLOCATABLE, TARGET :: drc(:,:) ! (ngm, ntyp)
  !! contain the rhoc (without structure fac) for all atomic types.
  !
END MODULE nlcc_ph
!
!
MODULE phus
  !
  !! These are additional variables needed for the linear response
  !! program with the US pseudopotentials.
  !
  USE kinds, ONLY :  DP
  USE becmod, ONLY : bec_type
#if defined(__CUDA)
  USE becmod_gpum, ONLY : bec_type_d
#endif
  !
  SAVE
  !
  REAL (DP), ALLOCATABLE :: alphasum(:,:,:,:) ! (nhm*(nhm+1)/2,3,nat,nspin)
  !! used to compute modes. It contains \(\sum_i \langle \psi_i| d/du
  !! (|\beta_n><beta_m|) | \psi_i\rangle + (m-n)\)
  ! dipole moment of each Q.
  ! These variables contains the five integrals defined in PRB 64, 35118 (2001)
  COMPLEX (DP), ALLOCATABLE :: int1(:,:,:,:,:)     ! nhm, nhm, 3, nat, nspin),&
  !! int1 -> \int V_eff d/du (Q) d^3r
  COMPLEX (DP), ALLOCATABLE :: int2(:,:,:,:,:)     ! nhm, nhm, 3,nat, nat),&
  !! int2 -> \int d/du (V_loc) Q d^3r
  COMPLEX (DP), ALLOCATABLE :: int4(:,:,:,:,:)     ! nhm*(nhm+1)/2, 3, 3, nat, nspin),&
  !! int4 -> \int V_eff d^2/dudu (Q) d^3r
  COMPLEX (DP), ALLOCATABLE :: int5(:,:,:,:,:)     ! nhm*(nhm+1)/2, 3, 3, nat, nat),&
  !! int5 -> \int d/du (V_loc) d/du (Q) d^3r
  COMPLEX (DP), ALLOCATABLE :: int1_nc(:,:,:,:,:)     ! nhm, nhm, 3, nat, nspin),&
  !! int1 - noncollinear
  COMPLEX (DP), ALLOCATABLE :: int2_so(:,:,:,:,:,:)   ! nhm, nhm, 3, nat,nat,nspin),&
  !! int2 - spin-orbit
  COMPLEX (DP), ALLOCATABLE :: int4_nc(:,:,:,:,:,:)   ! nhm, nhm, 3, 3, nat, nspin),&
  !! int4 - noncollinear
  COMPLEX (DP), ALLOCATABLE :: int5_so(:,:,:,:,:,:,:) ! nhm*(nhm+1)/2, 3, 3, nat, nat, nspin),&
  !! int5 - spin-orbit
  !
  !  int3 -> \int d\du (V_Hxc) Q d^3r .... generalized to Delta V_Hxc and move to lr_us in LR_Modules
  !
  COMPLEX (DP), ALLOCATABLE :: becsum_nc(:,:,:,:)     ! nhm*(nhm+1)/2,nat,npol,npol)
  !! it contains \(\sum_i \langle\psi_i | \beta_n\rangle\langle\beta_m| \psi_i \rangle + (m-n)\)
  COMPLEX (DP), ALLOCATABLE :: becsumort(:,:,:,:)     ! nhm*(nhm+1)/2,nat,nspin,3*nat)
  !! it contains \(\text{alphasum}+\sum_i \langle\psi_i | \beta_n\rangle\langle\beta_m| \delta \psi_i \rangle\)
  COMPLEX (DP), ALLOCATABLE :: alphasum_nc(:,:,:,:,:)   ! nhm*(nhm+1)/2,3,nat,npol,npol)
  !
  !
  type(bec_type),  ALLOCATABLE, TARGET :: alphap(:,:)  ! nkbtot, nbnd, 3, nksq)
#if defined(__CUDA)
  type(bec_type_d),  ALLOCATABLE, TARGET :: alphap_d(:,:)
#endif
  !! contains \( \langle d\du (\beta_n) | \psi_i \rangle\)
  !
END MODULE phus
!
!
MODULE partial
  !
  !! The variables needed for partial computation of dynamical matrix.
  !
  USE kinds, ONLY :  DP
  !
  SAVE
  !
  INTEGER, ALLOCATABLE :: atomo(:)
  !! (nat) : list of the atoms that moves
  INTEGER :: nat_todo
  !! number of atoms to compute
  INTEGER :: nat_todo_input
  !! nat_todo given in input
  LOGICAL, ALLOCATABLE :: comp_irr(:) !(3*nat)
  !! TRUE if this irr.rep. has to be computed
  LOGICAL, ALLOCATABLE :: done_irr(:) ! (3*nat)
  !! TRUE if this irr.rep. has been done
  LOGICAL :: all_comp
  !! if TRUE all representation have been computed
  INTERFACE 
    SUBROUTINE set_local_atomo(nat, nat_todo_, atomo_, nsym, irt,  nat_l, atomo_l) 
      IMPLICIT NONE 
      INTEGER,INTENT(IN)               :: nat, nat_todo_, nsym, atomo_(nat_todo_), irt(48,nat)
      !! :nat: total number of atoms
      !! :nat_todo: number of atoms effectively displaced 
      !! :nsym: number of symmetries in the system
      !! :atomo: list of atoms to be displaced before symmetrization 
      !! :irt: atoms corresponding atom for each sym operation and atom
      INTEGER,INTENT(OUT)              :: nat_l 
      !! actual number of atoms to be displaced considering symmetries
      INTEGER,ALLOCATABLE,INTENT(OUT)  :: atomo_l(:)
      !! list with the indeces of all the atoms to be displaced 
    END SUBROUTINE set_local_atomo
  END INTERFACE
END MODULE partial
!
MODULE gamma_gamma
  !
  INTEGER, ALLOCATABLE :: has_equivalent(:)
  !! 0 if the atom has to be calculated
  INTEGER, ALLOCATABLE :: with_symmetry(:)
  !! calculated by symmetry
  INTEGER, ALLOCATABLE :: n_equiv_atoms(:)
  !! number of equivalent atoms
  INTEGER, ALLOCATABLE :: equiv_atoms(:,:)
  !! which atoms are equivalent
  INTEGER :: n_diff_sites
  !! Number of different sites
  INTEGER :: nasr
  !! atom calculated with asr
  LOGICAL :: asr
  !! if TRUE apply the asr
  !
END MODULE gamma_gamma
!
MODULE control_ph
  !
  !! The variables controlling the phonon run.
  !
  USE kinds, ONLY :  DP
  USE parameters, ONLY: npk
  !
  SAVE
  !
  INTEGER, PARAMETER :: maxter = 150
  !! maximum number of iterations
  INTEGER :: niter_ph
  !! maximum number of iterations (read from input)
  INTEGER :: nmix_ph
  !! mixing type
  INTEGER :: start_irr
  !! initial representation
  INTEGER :: last_irr
  !! last representation of this run
  INTEGER :: current_iq
  !! current q point
  INTEGER :: start_q
  !! initial q in the list
  INTEGER :: last_q
  !! last_q in the list
  !
  REAL(DP) :: tr2_ph
  !! threshold for phonon calculation
  REAL(DP) :: alpha_mix(maxter)
  !! the mixing parameter
  CHARACTER(LEN=10) :: where_rec='no_recover'
  !! where the ph run recovered
  CHARACTER(LEN=12) :: electron_phonon
  CHARACTER(LEN=256) :: flmixdpot, tmp_dir_ph, tmp_dir_phq
  INTEGER :: rec_code=-1000
  !! code for recover
  INTEGER :: rec_code_read=-1000
  !! code for recover. Not changed during the run
  LOGICAL :: lgamma_gamma
  !! if TRUE this is a q=0 computation with k=0 only
  LOGICAL :: convt
  !! if TRUE the phonon has converged
  LOGICAL :: epsil
  !! if TRUE computes dielec. const and eff. charges
  LOGICAL :: done_epsil=.FALSE.
  !! TRUE when diel. constant is available
  LOGICAL :: trans
  !! if TRUE computes phonons
  LOGICAL :: zue
  !! if TRUE computes eff. charges as induced polarization
  LOGICAL :: done_zue=.FALSE.
  !! TRUE when the eff. charges are available
  LOGICAL :: zeu
  !! if TRUE computes eff. charges as induced forces
  LOGICAL :: done_zeu=.FALSE.
  !! TRUE when the eff. charges are available
  LOGICAL :: done_start_zstar=.FALSE.
  !
  LOGICAL :: only_wfc=.FALSE.
  !! if TRUE computes only bands
  LOGICAL :: only_init=.FALSE.
  !! if TRUE computes only initial stuff
  LOGICAL :: with_ext_images=.FALSE.
  !! if TRUE use an external driver to decide what each image does.
  LOGICAL :: always_run=.FALSE.
  !! if TRUE the code do not stop after doing partial representations
  !! always_run=.TRUE., only for testing purposes
  LOGICAL :: recover
  !! if TRUE the run restarts
  LOGICAL :: low_directory_check=.FALSE.
  !! if TRUE search on the phsave directory only the representations
  !! requested in input.
  LOGICAL :: ext_restart
  !! if TRUE there is a restart file
  LOGICAL :: ext_recover
  !! if TRUE there is a recover file
  LOGICAL :: lnoloc
  !! if TRUE calculates the dielectric constant neglecting local field effects
  LOGICAL :: search_sym=.TRUE.
  !! if TRUE search the mode symmetry
  LOGICAL :: search_sym_save=.TRUE.
  !! save search symmetry 
  LOGICAL :: lnscf
  !! if TRUE the run makes first a nscf calculation
  LOGICAL :: ldisp
  !! if TRUE the run calculates full phonon dispersion
  LOGICAL :: reduce_io
  !! if TRUE reduces needed I/O
  LOGICAL :: done_bands
  !! if TRUE the bands have been calculated
  LOGICAL :: bands_computed=.FALSE.
  !! if TRUE the bands were computed in this run
  LOGICAL :: nogg
  !! if TRUE gamma_gamma tricks are disabled
  LOGICAL :: u_from_file=.FALSE.
  !! if TRUE the u are on file
  LOGICAL :: recover_read=.FALSE.
  !! if TRUE the recover data have been read
  LOGICAL :: ldiag=.FALSE.
  !! if TRUE force the diagonalization
  LOGICAL :: lqdir=.FALSE.
  !! if TRUE each q writes in its directory
  LOGICAL :: qplot=.FALSE.
  !! if TRUE the q are read from input
  LOGICAL :: xmldyn=.FALSE.
  !! if TRUE the dynamical matrix is in xml form
  LOGICAL :: all_done
  !! if TRUE all representations have been done
  !
  LOGICAL :: newgrid=.FALSE.
  !! if TRUE use new k-point grid nk1,nk2,nk3
  INTEGER :: nk1,nk2,nk3, k1,k2,k3
  !! new Monkhorst-Pack k-point grid
  !
  CHARACTER(LEN=256) :: dftd3_hess 
  ! file from where the dftd3 hessian is read
  !
END MODULE control_ph
!
!
MODULE freq_ph
  !
  !! The variables for computing frequency dependent dielectric constant.
  !
  USE kinds,   ONLY : DP
  !
  SAVE
  !
  LOGICAL :: fpol
  !! if TRUE dynamic dielectric constant is computed
  LOGICAL :: done_fpol
  !! if TRUE all dynamic dielectric constant is computed
  !
  INTEGER :: nfs
  !! number of frequencies
  !
  INTEGER :: current_iu
  !! the current frequency 
  !
  REAL (KIND=DP), ALLOCATABLE :: fiu(:)
  !! values  of frequency
  !
  REAL (KIND=DP), ALLOCATABLE :: polar(:,:,:)
  !! values  of frequency
  LOGICAL, ALLOCATABLE :: comp_iu(:)
  !! values  of frequency to calculate in this ru
  !
  LOGICAL, ALLOCATABLE :: done_iu(:)
  !! values of frequency already calculated

  !
END MODULE freq_ph
!
!
MODULE units_ph
  !
  !! The units of the files and the record lengths.
  !
  SAVE
  !
  INTEGER :: iuvkb
  !! unit with vkb
  INTEGER :: iubar
  !! unit with the part DV_bare
  INTEGER :: lrbar
  !! length of the DV_bare
  INTEGER :: iuebar
  !! unit with the part DV_bare for the electric field
  INTEGER :: lrebar
  !! length of the DV_bare fro the electric field
  INTEGER :: iupsir
  !! unit with evc in real space
  INTEGER :: iudrhous, lrdrhous
  INTEGER :: iudyn
  !! the unit for the dynamical matrix
  INTEGER :: iupdyn
  !! the unit for the partial dynamical matrix
  INTEGER :: iunrec
  !! the unit with the recover data
  INTEGER :: iudvscf
  !! the unit where the delta Vscf is written
  INTEGER :: iudrho
  !! the unit where the delta rho is written
  INTEGER :: lrdrho
  !! the length of the deltarho files
  INTEGER :: iucom
  !! the unit of the bare commutator in US case
  INTEGER :: lrcom
  !! the length  of the bare commutator in US case
  INTEGER :: iudvkb3, lrdvkb3
  INTEGER :: iuint3paw
  !! the unit of the int3_paw coefficients
  INTEGER :: lint3paw
  !! the length of the int3_paw coefficients
  INTEGER :: iundnsscf
  !! the unit of dnsscf, for DFPT+U
  INTEGER :: iudvpsi
  !! unit of DV_SCF * psi
  INTEGER :: lrdvpsi
  !! length of DV_SCF * psi
  INTEGER :: iugauge
  !! Unit for reading and writing gauge information in ahc.f90
  !
  LOGICAL, ALLOCATABLE :: this_dvkb3_is_on_file(:), &
                          this_pcxpsi_is_on_file(:,:)
  !
END MODULE units_ph
!
!
MODULE output
  !
  !! The name of the files.
  !
  SAVE
  !
  CHARACTER (LEN=256) :: fildyn
  !! output file for the dynamical matrix
  CHARACTER (LEN=256) :: fildvscf
  !! output file for deltavscf
  CHARACTER (LEN=256) :: fildrho
  !! output file for deltarho
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
  INTEGER :: nq1, nq2, nq3
  !! number of q-points in each direction
  INTEGER :: nqs
  !! number of q points to be calculated
  REAL(DP), ALLOCATABLE :: x_q(:,:)
  !! coordinates of the q points
  REAL(DP), ALLOCATABLE :: wq(:)
  !! for plot
  REAL(DP), ALLOCATABLE :: omega_disp(:,:)
  LOGICAL, ALLOCATABLE :: lgamma_iq(:)
  !! if TRUE this q is gamma.
  LOGICAL, ALLOCATABLE :: done_iq(:)
  !! if TRUE this q point has been already calculated
  LOGICAL, ALLOCATABLE :: comp_iq(:)
  !! if TRUE this q point has to be calculated
  !
END MODULE disp

MODULE grid_irr_iq
   !
   INTEGER, ALLOCATABLE :: npert_irr_iq(:,:) 
   !! for each q and irr: the number of perturbations
   INTEGER, ALLOCATABLE :: irr_iq(:)
   !! number of irreducible representation per q point
   INTEGER, ALLOCATABLE :: nsymq_iq(:)
   !! dimension of the small group of q for each q
   !
   LOGICAL, ALLOCATABLE :: comp_irr_iq(:,:)
   !! for each q and irr: if TRUE this representation has to be calculated
   LOGICAL, ALLOCATABLE :: done_irr_iq(:,:)
   !! for each q and irr: if TRUE this representation has been already 
   !! calculated
   LOGICAL, ALLOCATABLE :: done_elph_iq(:,:)
   !! for each q and irr: if TRUE the elph of this representation has 
   !! been already calculated
   LOGICAL, ALLOCATABLE :: done_bands(:)
   !! nqs, if TRUE the bands of this q have been calculated
   !
END MODULE grid_irr_iq

MODULE ldaU_ph
  !
  USE kinds,      ONLY : DP
  USE parameters, ONLY : ntypx
  !
  SAVE
  ! ... atomic wfc's at k
  COMPLEX(DP), ALLOCATABLE, TARGET :: wfcatomk(:,:)
  !! atomic wfc at k
  COMPLEX(DP), ALLOCATABLE, TARGET :: dwfcatomk(:,:,:)
  !! derivative of atomic wfc at k
  COMPLEX(DP), ALLOCATABLE, TARGET :: sdwfcatomk(:,:)
  !! S * derivative of atomic wfc at k
  !
  ! ... atomic wfc's at k+q
  COMPLEX(DP), POINTER :: wfcatomkpq(:,:)
  !! atomic wfc at k+q
  COMPLEX(DP), POINTER :: dwfcatomkpq(:,:,:)
  !! derivative of atomic wfc at k+q
  COMPLEX(DP), POINTER :: sdwfcatomkpq(:,:)
  !! S * derivative of atomic wfc at k+q
  ! 
  COMPLEX(DP), ALLOCATABLE, TARGET :: dvkb(:,:,:)
  !! derivative of beta funtions at k  
  COMPLEX(DP), POINTER :: vkbkpq(:,:)
  !! beta funtions at k+q
  COMPLEX(DP), POINTER :: dvkbkpq(:,:,:)
  !! derivative of beta funtions at k+q
  !
  ! Various arrays for the response occupation matrix
  COMPLEX(DP), ALLOCATABLE :: dnsbare(:,:,:,:,:,:)
  !! bare derivative of ns
  COMPLEX(DP), ALLOCATABLE :: dnsbare_all_modes(:,:,:,:,:)
  !! bare derivative of ns for all modes
  COMPLEX(DP), ALLOCATABLE :: dnsscf_all_modes(:,:,:,:,:)
  !! SCF  derivative of ns for all modes
  COMPLEX(DP), ALLOCATABLE :: dnsorth(:,:,:,:,:)
  !! valence component of dns
  COMPLEX(DP), ALLOCATABLE :: dnsorth_cart(:,:,:,:,:,:)
  !! same as above, but in cart. coordinates
  !
  COMPLEX (DP), ALLOCATABLE :: proj1(:,:),    &
                               proj2(:,:),    &
                               projpb(:,:),   &
                               projpdb(:,:,:)
  ! Arrays to store scalar products between vectors
  ! projpb  = <psi|beta>
  ! projpdb = <psi|dbeta>
  !
  !
  LOGICAL  :: read_dns_bare
  !! if TRUE read the first bare derivative of ns from file
  CHARACTER(LEN=4) :: d2ns_type
  !! type of approximation to compute the second bare derivative 
  !! of atomic occupation matrix ns
  !
END MODULE ldaU_ph

MODULE nc_mag_aux
  USE kinds,      ONLY : DP
  SAVE
  
  COMPLEX (DP), ALLOCATABLE ::  &
                               deeq_nc_save(:,:,:,:,:), &
                               int1_nc_save(:,:,:,:,:,:), &
                               int3_save(:, :, :, :, :, :)
END MODULE nc_mag_aux

!MODULE qpoint_aux
!  USE kinds,      ONLY : DP
!  USE becmod,     ONLY : bec_type
!  SAVE
  
!  INTEGER, ALLOCATABLE :: ikmks(:)    ! index of -k for magnetic calculations

!  INTEGER, ALLOCATABLE :: ikmkmqs(:)  ! index of -k-q for magnetic calculations

!  TYPE(bec_type), ALLOCATABLE :: becpt(:), alphapt(:,:)

!END MODULE qpoint_aux

MODULE phcom
  USE dynmat
  USE eqv
  USE efield_mod
  USE nlcc_ph
  USE phus
  USE partial
  USE control_ph
  USE freq_ph
  USE units_ph
  USE output
  USE gamma_gamma
  USE disp
  USE grid_irr_iq
  USE ldaU_ph
  USE nc_mag_aux
!  USE qpoint_aux
END MODULE phcom
