!
! Copyright (C) 2001-2016 Quantum ESPRESSO group
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
 INTEGER :: nirr, nmodes
  ! number of irreducible representations contained in the dynamical matrix
  ! number of modes
  INTEGER, ALLOCATABLE, TARGET :: npert(:) !3 * nat )
  ! the number of perturbations per IR
  INTEGER :: npertx
  ! max number of perturbations per IR
  COMPLEX (DP), POINTER :: &
       u(:,:),                     &!  3 * nat, 3 * nat),
       t(:,:,:,:),                 &! npertx, npertx, 48,3 * nat),
       tmq(:,:,:)                   ! npertx, npertx, 3 * nat)
  ! the transformation modes patterns
  ! the mode for deltarho
  ! the symmetry in the base of the pattern
  ! the symmetry q<->-q in the base of the pa

  CHARACTER(15), ALLOCATABLE :: name_rap_mode(:) ! symmetry type of each mode
  INTEGER, ALLOCATABLE :: num_rap_mode(:)  ! number of the representation for
                                           ! each mode
  !
END MODULE modes
!
MODULE cryst_ph
   !  
   USE kinds,  ONLY : DP
   !
   SAVE
   !
   ! This modeule contains the variables that describe properties of the
   ! crystal that are needed by the phonon program and are not in pw data
   ! probably these variables should be in the common of pw
   ! These variables are sets immediately after reading the pw variables
   !
   LOGICAL :: magnetic_sym   ! true in the non-collinear magnetic case

END MODULE cryst_ph
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
       dyn(:,:),             &! 3 * nat, 3 * nat)
       dyn_rec(:,:)           ! 3 * nat, 3 * nat)
  ! the initial dynamical matrix
  ! the dynamical matrix
  ! the contribution of each representation to the dynamical matrix
  REAL (DP), ALLOCATABLE :: &
       w2(:)                  ! 3 * nat)
  ! omega^2
  !
END MODULE dynmat
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
       zstarue0(:,:),        &! 3 * nat, 3)
       zstarue0_rec(:,:)      ! 3 * nat, 3)
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
  !
END MODULE nlcc_ph
!
!
MODULE phus
  USE kinds, ONLY :  DP
  USE becmod, ONLY : bec_type
  !
  ! ... These are additional variables needed for the linear response
  ! ... program with the US pseudopotentials
  !
  SAVE
  !
  REAL (DP), ALLOCATABLE :: &
       alphasum(:,:,:,:)     ! nhm*(nhm+1)/2,3,nat,nspin)
                             ! used to compute modes
  ! alphasum contains \sum_i <psi_i| d/du (|\beta_n><beta_m|) | psi_i> + (m-n)
  ! dipole moment of each Q
  COMPLEX (DP), ALLOCATABLE :: &
       int1(:,:,:,:,:),     &! nhm, nhm, 3, nat, nspin),&
       int2(:,:,:,:,:),     &! nhm, nhm, 3,nat, nat),&
       int4(:,:,:,:,:),     &! nhm*(nhm+1)/2, 3, 3, nat, nspin),&
       int5(:,:,:,:,:),     &! nhm*(nhm+1)/2, 3, 3, nat, nat),&
       int1_nc(:,:,:,:,:),     &! nhm, nhm, 3, nat, nspin),&
       int2_so(:,:,:,:,:,:),   &! nhm, nhm, 3, nat,nat,nspin),&
       int4_nc(:,:,:,:,:,:),   &! nhm, nhm, 3, 3, nat, nspin),&
       int5_so(:,:,:,:,:,:,:), &! nhm*(nhm+1)/2, 3, 3, nat, nat, nspin),&
!
!  These variables contains the five integrals defined in PRB 64, 35118 (2001)
!  int1 -> \int V_eff d/du (Q) d^3r
!  int2 -> \int d/du (V_loc) Q d^3r
!  int3 -> \int d\du (V_Hxc) Q d^3r .... generalized to Delta V_Hxc and move to lr_us in LR_Modules
!  int4 -> \int V_eff d^2/dudu (Q) d^3r
!  int5 -> \int d/du (V_loc) d/du (Q) d^3r
!
       becsum_nc(:,:,:,:),     &! nhm*(nhm+1)/2,nat,npol,npol)
       becsumort(:,:,:,:),     &! nhm*(nhm+1)/2,nat,nspin,3*nat)
       alphasum_nc(:,:,:,:,:)   ! nhm*(nhm+1)/2,3,nat,npol,npol)
!
!  becsum contains \sum_i <\psi_i | \beta_n><\beta_m| \psi_i > + (m-n)
!  besumort contains alphasum+\sum_i <\psi_i | \beta_n><\beta_m| \delta \psi_i >
!
  type (bec_type),  ALLOCATABLE, TARGET :: &
       alphap(:,:)           ! nkbtot, nbnd, 3, nksq)
  !
  ! alphap contains < d\du (\beta_n) | psi_i>
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
       atomo(:)         ! (nat) : list of the atoms that moves
  INTEGER :: nat_todo,    & ! number of atoms to compute
             nat_todo_input ! nat_todo given in input
  LOGICAL, ALLOCATABLE :: &
       comp_irr(:),    &! (3*nat) : .true. if this irr.rep. has to be computed
       done_irr(:)      ! (3*nat) : .true. if this irr.rep. has been done
  LOGICAL :: all_comp   ! if .TRUE. all representation have been computed
  !
END MODULE partial
!
MODULE gamma_gamma
  INTEGER, ALLOCATABLE :: &
           has_equivalent(:),  &  ! 0 if the atom has to be calculated
           with_symmetry(:),   &  ! calculated by symmetry
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
  INTEGER, PARAMETER :: maxter = 100 ! maximum number of iterations
  INTEGER :: niter_ph,      & ! maximum number of iterations (read from input)
             nmix_ph,       & ! mixing type
             start_irr,     & ! initial representation
             last_irr,      & ! last representation of this run
             current_iq,    & ! current q point
             start_q, last_q  ! initial q in the list, last_q in the list
  REAL(DP) :: tr2_ph  ! threshold for phonon calculation
  REAL(DP) :: alpha_mix(maxter), & ! the mixing parameter
              time_now             ! CPU time up to now
  CHARACTER(LEN=10)  :: where_rec='no_recover'! where the ph run recovered
  CHARACTER(LEN=12) :: electron_phonon
  CHARACTER(LEN=256) :: flmixdpot, tmp_dir_ph, tmp_dir_phq
  INTEGER :: rec_code=-1000,    &! code for recover
             rec_code_read=-1000 ! code for recover. Not changed during the run
  LOGICAL :: lgamma_gamma,&! if .TRUE. this is a q=0 computation with k=0 only
             convt,       &! if .TRUE. the phonon has converged
             epsil,       &! if .TRUE. computes dielec. const and eff. charges
             done_epsil=.FALSE.,  &! .TRUE. when diel. constant is available
             trans,       &! if .TRUE. computes phonons
             zue,         &! if .TRUE. computes eff. charges as induced polarization
             done_zue=.FALSE., &! .TRUE. when the eff. charges are available
             zeu,         &! if .TRUE. computes eff. charges as induced forces
             done_zeu=.FALSE., &! .TRUE. when the eff. charges are available
             done_start_zstar=.FALSE., &!
             only_wfc=.FALSE.,  &! if .TRUE. computes only bands
             only_init=.FALSE.,  &! if .TRUE. computes only initial stuff
             with_ext_images=.FALSE., & ! if .TRUE. use an external driver
                                        ! to decide what each image does.
             always_run=.FALSE., & ! if .TRUE. the code do not stop after
                                   ! doing partial representations
             recover,     &! if .TRUE. the run restarts
             low_directory_check=.FALSE., & ! if .TRUE. search on the phsave 
                           ! directory only the representations requested 
                           ! in input.
             ext_restart, &! if .TRUE. there is a restart file
             ext_recover, &! if .TRUE. there is a recover file
             lnoloc,      &! if .TRUE. calculates the dielectric constant
                           ! neglecting local field effects
             search_sym=.TRUE.,  &! if .TRUE. search the mode symmetry
             search_sym_save=.TRUE.,  &! save search symmetry 
             lnscf,       &! if .TRUE. the run makes first a nscf calculation
             ldisp,       &! if .TRUE. the run calculates full phonon dispersion
             reduce_io,   &! if .TRUE. reduces needed I/O
             done_bands,  &! if .TRUE. the bands have been calculated
             bands_computed=.FALSE., & ! if .TRUE. the bands were computed
                                       ! in this run
             nogg,        &! if .TRUE. gamma_gamma tricks are disabled
             u_from_file=.FALSE.,  & ! if true the u are on file
             recover_read=.FALSE., & ! if true the recover data have been read
             ldiag=.FALSE.,        & ! if true force the diagonalization
             lqdir=.FALSE.,        & ! if true each q writes in its directory
             qplot=.FALSE.,        & ! if true the q are read from input
             xmldyn=.FALSE.,   & ! if true the dynamical matrix is in xml form
             all_done, &      ! if .TRUE. all representations have been done
             newgrid=.FALSE.  ! if .TRUE. use new k-point grid nk1,nk2,nk3
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
  LOGICAL :: fpol, & ! if .TRUE. dynamic dielectric constant is computed
             done_fpol ! if .TRUE. all dynamic dielectric constant is computed
  !
  INTEGER :: nfs                   ! # of frequencies
  !
  INTEGER :: current_iu            ! the current frequency 
  !
  REAL (KIND=DP), ALLOCATABLE :: fiu(:)    ! values  of frequency
  !
  REAL (KIND=DP), ALLOCATABLE :: polar(:,:,:)    ! values  of frequency

  LOGICAL, ALLOCATABLE :: comp_iu(:) ! values  of frequency to calculate in this ru
  !
  LOGICAL, ALLOCATABLE :: done_iu(:)    ! values of frequency already calculated

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
       iuwfc,     & ! iunit with the wavefunctions
       lrwfc,     & ! the length of wavefunction record
       iuvkb,     & ! unit with vkb
       iubar,     & ! unit with the part DV_{bare}
       lrbar,     & ! length of the DV_{bare}
       iuebar,    & ! unit with the part DV_{bare} for the electric field
       lrebar,    & ! length of the DV_{bare} fro the electric field
       iudwf,     & ! unit with D psi
       iupsir,    & ! unit with evc in real space
       lrdwf,     & ! length of D psi record
       iudrhous, lrdrhous, &
       iudyn,     & ! the unit for the dynamical matrix
       iupdyn,    & ! the unit for the partial dynamical matrix
       iunrec,    & ! the unit with the recover data
       iudvscf,   & ! the unit where the delta Vscf is written
       iudrho,    & ! the unit where the delta rho is written
       lrdrho,    & ! the length of the deltarho files
       iucom,     & ! the unit of the bare commutator in US case
       lrcom,     & ! the length  of the bare commutator in US case
       iudvkb3, lrdvkb3, &
       iuint3paw, & ! the unit of the int3_paw coefficients
       lint3paw     ! the length of the int3_paw coefficients
  ! the unit with the products
  ! the length of the products

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
  INTEGER :: nq1, nq2, nq3  ! number of q-points in each direction
  INTEGER :: nqs            ! number of q points to be calculated
  REAL(DP), ALLOCATABLE :: x_q(:,:), & ! coordinates of the q points
                        wq(:) ! for plot

  REAL(DP), ALLOCATABLE :: omega_disp(:,:)

  LOGICAL, ALLOCATABLE :: &
       lgamma_iq(:),    &! if .true. this q is gamma.
       done_iq(:),      &! if .true. this q point has been already calculated
       comp_iq(:)        ! if .true. this q point has to be calculated
  !
END MODULE disp

MODULE grid_irr_iq

   INTEGER, ALLOCATABLE ::  &
       npert_irr_iq(:,:),&! for each q and irr: the number of perturbations
       irr_iq(:),        &! number of irreducible representation per q point
       nsymq_iq(:)        ! dimension of the small group of q for each q

  LOGICAL, ALLOCATABLE ::  &
       comp_irr_iq(:,:),   & ! for each q and irr: if .TRUE. this 
                             ! representation has  to be calculated
       done_irr_iq(:,:),   & ! for each q and irr: if .TRUE. this 
                             ! representation has been already calculated
       done_elph_iq(:,:),   & ! for each q and irr: if .TRUE. the elph of this 
                             ! representation has been already calculated
       done_bands(:)         ! nqs, if .TRUE. the bands of this q have been 
                             ! calculated
END MODULE grid_irr_iq

!
!
MODULE phcom
  USE modes
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
END MODULE phcom
