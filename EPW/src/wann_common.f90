  !
  ! Copyright (C) 2016-2023 EPW-Collaboration
  ! Copyright (C) 2010-2016 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino
  !
  ! Copyright (C) 2003 PWSCF group
  ! This file is distributed under the terms of the
  ! GNU General Public License. See the file `License'
  ! in the root directory of the present distribution,
  ! or http://www.gnu.org/copyleft/gpl.txt .
  !
  ! Adapted from QE.
  !
  !-----------------------------------------------------------------------------------------
  MODULE wann_common
  !-----------------------------------------------------------------------------------------
  !!
  !! Interface between Wannier and EPW
  !!
  USE kinds, ONLY : DP
  !
  CHARACTER(LEN = 15) :: wan_mode
  !! running mode
  CHARACTER(LEN = 256) :: seedname2
  !! prepended to file names in wannier90. For implementation of wannier_lib
  CHARACTER(LEN = 3), ALLOCATABLE :: atsym(:)
  !! atomic symbols. atsym(nat)
  LOGICAL :: logwann
  !!
  LOGICAL :: write_unk
  !! Set to .TRUE. to write the periodic part of the Bloch functions. Default is .FALSE.
  LOGICAL :: wvfn_formatted
  !! Set to .TRUE. to write formatted wavefunctions. Default is .FALSE. (only relevant if write_unk=.TRUE.)
  LOGICAL :: write_amn
  !! write A_mn(k) matrices to file (not used in library mode)
  LOGICAL :: write_mmn
  !! write M_mn(k,b) matrices to file
  LOGICAL :: write_spn
  !! write S matrices between Bloch states (non-collinear spin calculation only)
  LOGICAL, ALLOCATABLE :: zerophase(:, :)
  !! Phase
  LOGICAL, ALLOCATABLE :: lwindow(:, :)
  !! lwindow(iband,ik) is .TRUE. if the band ibnd lies within the outer enery window ak k-point ik
  !! lwindow(num_bands,iknum)
  LOGICAL, ALLOCATABLE :: excluded_band(:)
  !! list of bands to exclude from the calculation of WFs
  INTEGER :: nnb
  !! number of b-vectors
  INTEGER :: iun_nnkp
  !! unit to write nnkp
  INTEGER :: iun_mmn
  !! unit to write mmn matrices
  INTEGER :: iun_amn
  !! unit to write amn matrices
  INTEGER :: iun_band
  !! unit to write eigenvalues
  INTEGER :: iun_spn
  !! unit to write spinors for noncollinear calculations
  INTEGER :: iun_plot
  !! unit to write wave functions
  INTEGER :: nnbx
  !! max number of b-vectors
  INTEGER :: nexband
  !! number of excluded bands
  INTEGER :: mp_grid(3)
  !! dimensions of MP k-point grid
  INTEGER :: n_wannier
  !! number of WFs
  INTEGER :: n_proj
  !! number of projections (=#WF unless spinors then =#WF/2)
  INTEGER :: ispinw
  !! ispinw = 1, ikstart = 1, ikstop=nkstot/2 for spin-up
  INTEGER :: ikstart
  !! ispinw=2, ikstart=nkstot/2+1, ikstop=nkstot for spin-down
  INTEGER :: ikstop
  !! ispinw=0, ikstart = 1, ikstop=nkstot for unpolarized and non-collinear
  INTEGER :: iknum
  !! number of k-points, iknum = nkstot/2 for spin-polarized case, iknum = nkstot for unpolarized and non-collinear
  INTEGER :: num_bands
  !! number of bands left after exclusions
  INTEGER :: num_nnmax = 12
  !!
  INTEGER, ALLOCATABLE :: kpb(:, :)
  !! list of nearest neighbours for eack k-point k+b(ik,ib)
  INTEGER, ALLOCATABLE :: g_kpb(:, :, :)
  !! G_k+b(ipol,ik,ib) vector that brings the nn-th nearest neighbour of each
  !! k-point to its periodic image (crystal coords.)
  INTEGER, ALLOCATABLE :: ig_(:, :)
  !! index of G_k+b(ipol,ik,ib)
  INTEGER, ALLOCATABLE :: lw(:, :)
  !! l of wannier (16, n_wannier)
  INTEGER, ALLOCATABLE :: mw(:, :)
  !! m of wannier (16, n_wannier)
  INTEGER, ALLOCATABLE :: num_sph(:)
  !! num. func. in lin. comb., num_sph(n_wannier)
  INTEGER, ALLOCATABLE :: spin_eig(:)
  !! "1" or "-1" to denote projection onto up or down spin states
  INTEGER, ALLOCATABLE :: l_w(:)
  !! angular part l and mr of wannier (n_wannier) as from table 3.1,3.2 of spec.
  INTEGER, ALLOCATABLE :: mr_w(:)
  !!! angular part l and mr of wannier (n_wannier) as from table 3.1,3.2 of spec.
  INTEGER, ALLOCATABLE :: r_w(:)
  !! radial part of wannier (n_wannier) as from table 3.3 of spec.
  REAL(KIND = DP) :: rlatt(3, 3)
  !! real lattices (Cartesian coords., units of Angstrom)
  REAL(KIND = DP) :: glatt(3, 3)
  !! recip. lattices (Cartesian coords., units of reciproval Angstrom)
  REAL(KIND = DP) :: spreads(3)
  !! values of \Omega, \Omega_I, and \tilde{\Omega}
  REAL(KIND = DP), ALLOCATABLE :: center_w(:, :)
  !! projection function center (crystal coords.) center_w(3,n_wannier)
  REAL(KIND = DP), ALLOCATABLE :: spin_qaxis(:, :)
  !! spin quantisation axis (Cartesian coords.)
  REAL(KIND = DP), ALLOCATABLE :: zaxis(:, :)
  !! defines the axis from which the polar angle \theta in spherical coords. i is measured. zaxis(3,n_wannier)
  REAL(KIND = DP), ALLOCATABLE :: xaxis(:, :)
  !! defines the axis from which the azimuthal angle \phi in spherical coords is measured. xaxis(3,n_wannier)
  REAL(KIND = DP), ALLOCATABLE :: alpha_w(:)
  !! alpha_w(n_wannier) (called proj_zona in wannier spec)
  REAL(KIND = DP), ALLOCATABLE :: csph(:, :)
  !! expansion coefficients of gf on QE ylm function (16,n_wannier)
  REAL(KIND = DP), ALLOCATABLE :: kpt_latt(:, :)
  !! k-points in crystal coords. kpt_latt(3,iknum)
  REAL(KIND = DP), ALLOCATABLE :: atcart(:, :)
  !! atom centres (Cartesian coords., units of Angstrom). atcart(3,nat)
  REAL(KIND = DP), ALLOCATABLE :: wann_centers(:, :)
  !! centers of WFs (Cartesian coords., units of Angstrom) wann_centers(3,n_wannier)
  REAL(KIND = DP), ALLOCATABLE :: wann_spreads(:)
  !! spread of WFs in Angstrom^2. wann_spreads(n_wannier)
  REAL(KIND = DP), ALLOCATABLE :: eigval(:, :)
  !! eigenvalues \epsilon_nk corresponding to the \psi_nk eigenstates
  COMPLEX(KIND = DP), ALLOCATABLE :: gf(:, :)
  !! guding_function(npwx,n_proj)
  COMPLEX(KIND = DP), ALLOCATABLE :: gf_spinor(:, :)
  !!
  COMPLEX(KIND = DP), ALLOCATABLE :: sgf_spinor(:, :)
  !!
  COMPLEX(KIND = DP), ALLOCATABLE :: m_mat(:, :, :, :)
  !! overlap matrices between neighbouring periodic parts of Bloch eigenstates at each k-point M_mn(k,b)=<u_mk|u_nk+b>
  !! m_mat(num_bands,num_bands,nnb,iknum)
  COMPLEX(KIND = DP), ALLOCATABLE :: a_mat(:, :, :)
  !! matrices describing the projection of n_wannier trial orbitals
  !! on num_bands Bloch states at each k-point A_mn(k)=<psi_mk|g_n> a_mat(num_bands,n_wannier,iknum)
  COMPLEX(KIND = DP), ALLOCATABLE :: u_mat(:, :, :)
  !! unitary matrix at each k-point u_mat(n_wannier,n_wannier,iknum)
  COMPLEX(KIND = DP), ALLOCATABLE :: u_mat_opt(:, :, :)
  !! unitary matrix for the optimal sub-space at each k-point u_mat_opt(num_bands,n_wannier,iknum)
  !-----------------------------------------------------------------------------------------
  END MODULE wann_common
  !-----------------------------------------------------------------------------------------
