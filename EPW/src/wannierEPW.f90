  !
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
MODULE  wannierEPW
   USE kinds, only : DP
   !INTEGER, ALLOCATABLE :: nnb(:)       ! #b  (ik)
   INTEGER              :: nnb          ! number of b-vectors
   INTEGER, ALLOCATABLE :: kpb(:,:)     ! list of nearest neighbours for eack k-point k+b(ik,ib)
   INTEGER, ALLOCATABLE :: g_kpb(:,:,:) ! G_k+b(ipol,ik,ib)
                                        ! vector that brings the nn-th nearest neighbour of each
                                        ! k-point to its periodic image (crystal coords.)
   INTEGER, ALLOCATABLE :: ig_(:,:)     ! index of G_k+b(ipol,ik,ib)
   INTEGER, ALLOCATABLE :: lw(:,:), mw(:,:) ! l and m of wannier (16,n_wannier)
   INTEGER, ALLOCATABLE :: num_sph(:)   ! num. func. in lin. comb., num_sph(n_wannier)
   LOGICAL, ALLOCATABLE :: excluded_band(:) ! list of bands to exclude from the calculation of WFs
   INTEGER  :: iun_nnkp, iun_mmn, iun_amn, iun_band, iun_spn, iun_plot, nnbx, nexband
   INTEGER  :: n_wannier ! number of WFs
   INTEGER  :: n_proj    ! number of projections (=#WF unless spinors then =#WF/2)
   COMPLEX(DP), ALLOCATABLE :: gf(:,:)  ! guding_function(npwx,n_proj)
   COMPLEX(DP), ALLOCATABLE :: gf_spinor(:,:)
   COMPLEX(DP), ALLOCATABLE :: sgf_spinor(:,:)
   INTEGER  :: ispinw, ikstart, ikstop ! ispinw=1, ikstart=1, ikstop=nkstot/2 for spin-up
                                       ! ispinw=2, ikstart=nkstot/2+1, ikstop=nkstot for spin-down
                                       ! ispinw=0, ikstart=1, ikstop=nkstot for unpolarized and non-collinear
   INTEGER  :: iknum  ! number of k-points
                      ! iknum = nkstot/2 for spin-polarized case
                      ! iknum = nkstot for unpolarized and non-collinear
   CHARACTER(LEN=15)  :: wan_mode    ! running mode
   LOGICAL            :: logwann
   LOGICAL            :: write_unk   ! Set to .true. to write the periodic part of the Bloch functions. Default is .false.
   LOGICAL            :: reduce_unk  ! Set to .true. to reduce file-size (and resolution) of Bloch functions by a factor of 8.
                                     ! Default is .false. (only relevant if write_unk=.true.)
   LOGICAL            :: wvfn_formatted ! Set to .true. to write formatted wavefunctions.
                                        ! Default is .false. (only relevant if write_unk=.true.)
   LOGICAL            :: write_amn      ! write A_mn(k) matrices to file (not used in library mode)
   LOGICAL            :: write_mmn      ! write M_mn(k,b) matrices to file
   LOGICAL            :: write_spn      ! write S matrices between Bloch states (non-collinear spin calculation only)
   ! input data from nnkp file
   REAL(DP), ALLOCATABLE :: center_w(:,:)   ! projection function center (crystal coords.) center_w(3,n_wannier)
   INTEGER,  ALLOCATABLE :: spin_eig(:)     ! "1" or "-1" to denote projection onto up or down spin states
   REAL(DP), ALLOCATABLE :: spin_qaxis(:,:) ! spin quantisation axis (Cartesian coords.)
   INTEGER, ALLOCATABLE  :: l_w(:), mr_w(:) ! angular part l and mr of wannier (n_wannier) as from table 3.1,3.2 of spec.
   INTEGER, ALLOCATABLE  :: r_w(:)          ! radial part of wannier (n_wannier) as from table 3.3 of spec.
   REAL(DP), ALLOCATABLE :: zaxis(:,:) ! defines the axis from which the polar angle \theta in spherical coords. i
                                       ! is measured. zaxis(3,n_wannier)
   REAL(DP), ALLOCATABLE :: xaxis(:,:) ! defines the axis from which the azimuthal angle \phi in spherical coords.
                                       ! is measured. xaxis(3,n_wannier)
   REAL(DP), ALLOCATABLE :: alpha_w(:) ! alpha_w(n_wannier) ( called proj_zona in wannier spec)
   !
   REAL(DP), ALLOCATABLE :: csph(:,:)  ! expansion coefficients of gf on QE ylm function (16,n_wannier)
   ! SP: We had to make it an array in order to be able to allocate & pass the
   ! info to the wannier_lib
   CHARACTER(len=256) :: seedname2 ! prepended to file names in wannier90
   ! For implementation of wannier_lib
   INTEGER               :: mp_grid(3)     ! dimensions of MP k-point grid
   REAL(DP)              :: rlatt(3,3)     ! real lattices (Cartesian coords., units of Angstrom)
   REAL(DP)              :: glatt(3,3)     ! recip. lattices (Cartesian coords., units of reciproval Angstrom)
   REAL(DP), ALLOCATABLE :: kpt_latt(:,:)  ! k-points in crystal coords. kpt_latt(3,iknum)
   REAL(DP), ALLOCATABLE :: atcart(:,:)    ! atom centres (Cartesian coords., units of Angstrom). atcart(3,nat)
   INTEGER               :: num_bands      ! number of bands left after exclusions
   CHARACTER(len=3), ALLOCATABLE :: atsym(:) ! atomic symbols. atsym(nat)
   INTEGER               :: num_nnmax=12
   COMPLEX(DP), ALLOCATABLE :: m_mat(:,:,:,:) ! overlap matrices between neighbouring periodic parts of Bloch
                                              ! eigenstates at each k-point M_mn(k,b)=<u_mk|u_nk+b>
                                              ! m_mat(num_bands,num_bands,nnb,iknum)
   COMPLEX(DP), ALLOCATABLE :: a_mat(:,:,:)   ! matrices describing the projection of n_wannier trial orbitals
                                              ! on num_bands Bloch states at each k-point A_mn(k)=<psi_mk|g_n>
                                              ! a_mat(num_bands,n_wannier,iknum)
   COMPLEX(DP), ALLOCATABLE :: u_mat(:,:,:)   ! unitary matrix at each k-point u_mat(n_wannier,n_wannier,iknum)
   COMPLEX(DP), ALLOCATABLE :: u_mat_opt(:,:,:) ! unitary matrix for the optimal sub-space at each k-point
                                                ! u_mat_opt(num_bands,n_wannier,iknum)
   LOGICAL, ALLOCATABLE     :: lwindow(:,:) ! lwindow(iband,ik) is .true. if the band ibnd lies within the
                                            ! outer enery window ak k-point ik
                                            ! lwindow(num_bands,iknum)
   REAL(DP), ALLOCATABLE    :: wann_centers(:,:) ! centers of WFs (Cartesian coords., units of Angstrom)
                                                 ! wann_centers(3,n_wannier)
   REAL(DP), ALLOCATABLE    :: wann_spreads(:)   ! spread of WFs in Angstrom^2. wann_spreads(n_wannier)
   REAL(DP)                 :: spreads(3)        ! values of \Omega, \Omega_I, and \tilde{\Omega}
   REAL(DP), ALLOCATABLE    :: eigval(:,:)       ! eigenvalues \epsilon_nk corresponding to the \psi_nk eigenstates
                                                 ! eigval(num_bands,iknum)
END MODULE  wannierEPW
