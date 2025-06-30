!
! Copyright (C) 2001-2025 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
!
! ... Common variables for LR_Modules routines
!
MODULE qpoint
  !
  USE kinds,      ONLY : DP
  !
  ! ... The variables needed to specify various indices,
  ! ... number of plane waves and k points and their coordiantes.
  !
  SAVE
  !
  INTEGER, POINTER :: igkq(:)     ! npwx)
  ! correspondence k+q+G <-> G
  INTEGER :: nksq, npwq, nksqtot
  ! the real number of k points
  ! the number of plane waves for q
  ! the total number of q points
  INTEGER, ALLOCATABLE :: ikks(:), ikqs(:)
  ! the index of k point in the list of k
  ! the index of k+q point in the list of k
  REAL (DP) :: xq(3)
  ! the coordinates of the q point
  COMPLEX (DP), ALLOCATABLE :: eigqts(:) ! nat)
  ! the phases associated to the q
  REAL (DP), ALLOCATABLE :: xk_col(:,:)
  !
END MODULE qpoint
!
!
!
MODULE qpoint_aux
  USE kinds,      ONLY : DP
  USE becmod,     ONLY : bec_type
 
  SAVE
  
  INTEGER, ALLOCATABLE :: ikmks(:)    ! index of -k for magnetic calculations

  INTEGER, ALLOCATABLE :: ikmkmqs(:)  ! index of -k-q for magnetic calculations

  TYPE(bec_type), ALLOCATABLE :: becpt(:), alphapt(:,:)

END MODULE qpoint_aux
!
!
MODULE control_lr
  !
  USE kinds,      ONLY : DP
  !
  ! ... The variables controlling the run of linear response codes
  !
  SAVE
  !
  INTEGER, ALLOCATABLE :: nbnd_occ(:)  ! occupied bands in metals
  REAL(DP) :: alpha_pv       ! the alpha value for shifting the bands
  LOGICAL  :: lgamma         ! if .TRUE. this is a q=0 computation
  LOGICAL  :: lrpa           ! if .TRUE. uses the Random Phace Approximation
  REAL(DP) :: ethr_nscf      ! convergence threshol for KS eigenvalues in the
                             ! NSCF calculation
  ! Sternheimer case 
  LOGICAL :: lgamma_gamma
  !! if TRUE this is a q=0 computation with k=0 only
  LOGICAL :: ext_recover, &! if .TRUE. there is a recover file
             lnoloc        ! if .TRUE. calculates the dielectric constant
                           ! neglecting local field effects
  !
  ! Variables for recover
  !
  CHARACTER(LEN=10) :: where_rec = 'no_recover'
  !! where the ph run recovered
  INTEGER :: rec_code = -1000
  !! code for recover
  INTEGER :: rec_code_read = -1000  ! code for recover. Not changed during the run
  !
  INTEGER :: nbnd_occx        ! maximun value of nbnd_occ(:)
  LOGICAL :: reduce_io
  !! if TRUE reduces needed I/O
  !
  ! Parameters controlling DFPT self-consistent iteration
  !
  INTEGER, PARAMETER :: maxter = 150
  !! maximum number of iterations
  LOGICAL :: convt
  !! if TRUE the DFPT has converged
  CHARACTER(LEN=256) :: flmixdpot
  !! File for writing history for potential mixing
  INTEGER :: niter_ph
  !! maximum number of iterations (read from input)
  INTEGER :: nmix_ph
  !! mixing type
  REAL(DP) :: tr2_ph
  !! threshold for DFPT calculation
  REAL(DP) :: alpha_mix(maxter)
  !! the mixing parameter
  !
  LOGICAL :: lmultipole=.FALSE.   
  !! if TRUE macroscopic density response to q-potential perturbation is written as output
  LOGICAL :: lnolr=.FALSE.   
  !! if TRUE G=0 component of the Hartree term is not added in dv_of_drho
  !
END MODULE control_lr
!
MODULE eqv
  !
  USE kinds,  ONLY : DP
  !
  ! ... The variables describing the linear response problem
  !
  SAVE
  !
  COMPLEX (DP), POINTER :: evq(:,:)
  ! the wavefunctions at point k+q
  COMPLEX (DP), ALLOCATABLE :: dvpsi(:,:), dpsi(:,:)
  ! the product of dV psi
  ! the change of the wavefunctions
  COMPLEX (DP), ALLOCATABLE :: drhos(:,:,:)
  !! the change of the density (smooth part only, dffts)
  REAL (DP), ALLOCATABLE :: dmuxc(:,:,:)        ! nrxx, nspin, nspin)
  ! the derivative of the xc potential
  REAL (DP), ALLOCATABLE, TARGET :: vlocq(:,:)  ! ngm, ntyp)
  ! the local potential at q+G
  !
END MODULE eqv
!
MODULE gc_lr
  !
  USE kinds, ONLY : DP
  !
  ! ... The variables needed for gradient corrected calculations
  !
  SAVE
  !
  REAL (DP), ALLOCATABLE :: &
       grho(:,:,:),    &! gradient of the unperturbed density  (3,nrxx,nspin)
       gmag(:,:,:),    &! 3, nrxx, nspin)
       vsgga(:),       &! nrxx)
       segni(:),       &! nrxx)
       dvxc_rr(:,:,:), &! derivatives of the E_xc functional w.r.t. r and s  
       dvxc_sr(:,:,:), &! r=rho and s=|grad(rho)|
       dvxc_ss(:,:,:), &! dimensions: (nrxx, nspin, nspin)
       dvxc_s(:,:,:)
  !
  ! in the noncollinear case gmag contains the gradient of the magnetization
  ! grho the gradient of rho+ and of rho-, the eigenvalues of the spin density
  ! vsgga= 0.5* (V_up-V_down) to be used in the calculation of the change
  ! of the exchange and correlation magnetic field.
  !
END MODULE gc_lr
!
MODULE lr_symm_base
  !
  USE kinds,  ONLY : DP
  !
  ! ... The variables needed to describe the modes and the small group of q
  !
  SAVE
  !
  INTEGER :: irgq(48), nsymq=0, irotmq
  ! selects the operations of the small group
  ! the number of symmetry of the small group
  ! selects the symmetry sending q <-> -q+G
  REAL (DP), ALLOCATABLE :: rtau(:,:,:) !3, 48, nat)
  ! coordinates of direct translations
  REAL (DP) :: gi(3,48), gimq(3)
  ! the possible G associated to each symmetry
  ! the G associated to the symmetry q<->-q+G
  LOGICAL :: minus_q, & ! if .TRUE. there is the symmetry sending q<->-q
             invsymq    ! if .TRUE. the small group of q has inversion
  !
  ! Symmetry representation of the perturbations
  !
  INTEGER :: lr_npert
  !! Number of perturbations considered at the same time.
  !! e.g., for phonons: dimension of the irreducible representation
  !! e.g., for electric fields: 3
  COMPLEX(DP), ALLOCATABLE :: upert(:, :, :)
  !! Representation of the symmetry in the perturbation basis. Size (lr_npert, lr_npert, nsymq)
  !! e.g., for phonons: transformation matrix of the patterns
  !! e.g., for electric fields: transformation matrix of Cartesian vectors
  COMPLEX(DP), ALLOCATABLE :: upert_mq(:, :)
  !! Representation of the symmetry that transforms q to -q. Size (lr_npert, lr_npert)
  !
END MODULE lr_symm_base
!
MODULE lrus
  !
  USE kinds,  ONLY : DP
  USE becmod, ONLY : bec_type
  !
  ! ... These are additional variables needed for the linear response
  ! ... with US pseudopotentials and a generic perturbation Delta Vscf
  !
  SAVE
  !
  COMPLEX (DP), ALLOCATABLE :: &
       int3(:,:,:,:,:),     &! nhm, nhm, nat, nspin, npert)
       int3_paw(:,:,:,:,:), &! nhm, nhm, nat, nspin, npert)
       int3_nc(:,:,:,:,:),  &! nhm, nhm, nat, nspin, npert)
       intq(:,:,:),         &! nhm, nhm, nat)
       intq_nc(:,:,:,:)      ! nhm, nhm, nat, nspin)
  ! int3 -> \int (Delta V_Hxc) Q d^3r
  ! similarly for int_nc while
  ! int3_paw contains Delta (D^1-\tilde D^1)
  ! intq integral of e^iqr Q
  ! intq_nc integral of e^iqr Q in the noncollinear case
  !
  REAL (DP), ALLOCATABLE ::    dpqq(:,:,:,:)       ! nhm, nhm, 3, ntyp)
  COMPLEX (DP), ALLOCATABLE :: dpqq_so(:,:,:,:,:)  ! nhm, nhm, nspin, 3, ntyp)
  ! dpqq and dpqq_so: dipole moment of each Q multiplied by the fcoef factors 
  !
  type (bec_type), ALLOCATABLE, TARGET :: becp1(:) ! nksq)
  ! becp1 contains < beta_n | psi_i >
  !
  REAL (DP),    ALLOCATABLE :: bbg(:,:)      ! nkb, nkb)
  ! for gamma_only     
  COMPLEX (DP), ALLOCATABLE :: bbk(:,:,:)    ! nkb, nkb, nks)
  ! for k points
  COMPLEX (DP), ALLOCATABLE :: bbnc(:,:,:) ! nkb*npol, nkb*npol, nks)
  ! for the noncollinear case
  ! bbg = < beta^N_i | beta^P_j > 
  ! bbg/bbk/bbnc are the scalar products of beta functions 
  ! localized on atoms N and P.
  !
END MODULE lrus
!
MODULE units_lr
  !
  USE kinds,  ONLY : DP
  !
  ! ... These are the units used in the linear response calculations
  !
  SAVE
  !
  INTEGER :: iuwfc,   & ! unit for wavefunctions
             lrwfc,   & ! the length of wavefunction record
             iuatwfc, & ! unit for atomic wavefunctions
             iuatswfc,& ! unit for atomic wavefunctions * S
             iudwf,   & ! unit with D psi
             lrdwf      ! length of D psi record
  !
END MODULE units_lr

MODULE ldaU_lr
  !
  USE kinds,      ONLY : DP
  USE parameters, ONLY : ntypx
  !
  REAL(DP) :: effU(ntypx)
  ! effective Hubbard parameter: effU = Hubbard_U - Hubbard_J0
  ! TODO: Can be moved to PW/ldaU
  !
  COMPLEX(DP), ALLOCATABLE :: dnsscf(:,:,:,:,:)
  !! SCF derivative of ns
  !
  COMPLEX(DP), ALLOCATABLE :: vh_u_save(:,:,:,:,:)
  COMPLEX(DP), ALLOCATABLE :: vh_uv_save(:,:,:,:,:,:)
  ! to save the two unperturbed Hubbard U and V potentials;
  ! one normal, and the other with the time-reversed m_hubb
  !
  COMPLEX(DP), ALLOCATABLE, TARGET :: swfcatomk(:,:)
  !! S * atomic wfc at k
  COMPLEX(DP), POINTER :: swfcatomkpq(:,:)
  !! S * atomic wfc at k+q
  !
  LOGICAL :: lr_has_dnsorth = .FALSE.
  !! If true, add lr_dnsorth to dnsscf.
  COMPLEX(DP), ALLOCATABLE :: lr_dnsorth(:, :, :, :, :)
  !! Fixed term to be added to dnsscf. Size (ldim, ldim, nspin, nat, 3*nat)
  !
END MODULE ldaU_lr
