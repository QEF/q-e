  !
  ! Copyright (C) 2016-2023 EPW-Collaboration
  ! Copyright (C) 2010-2016 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino
  ! Copyright (C) 2007-2009 Jesse Noffsinger, Brad Malone, Feliciano Giustino
  !
  ! This file is distributed under the terms of the GNU General Public
  ! License. See the file `LICENSE' in the root directory of the
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .
  !
  !--------------------------------------------------------------------------
  MODULE tdbe_eph_common
  !--------------------------------------------------------------------------
  !!
  !! Global variables for time-dependent Boltzmann simulations
  !! Authored by Yiming Pan (pan@physik.uni-kiel.de) 
  !! and Fabio Caruso (caruso@physik.uni-kiel.de) (Apr. 2025)
  !!
  USE kinds, ONLY : DP
  !
  IMPLICIT NONE
  !
  SAVE
  !  
  REAL(KIND = DP), ALLOCATABLE :: enk_all(:, :)
  !! Electron energies 
  REAL(KIND = DP), ALLOCATABLE :: wnuq_all(:, :)
  !! Phonon freqencies
  REAL(KIND = DP), ALLOCATABLE :: wkf_all(:) 
  !! weights of k points
  REAL(KIND = DP), ALLOCATABLE :: wqf_all(:) 
  !! weights of q points
  REAL(KIND = DP), ALLOCATABLE :: xkf_all(:,:)
  !! crystal coordinates of k points
  REAL(KIND = DP), ALLOCATABLE :: xqf_all(:,:)
  !! crystal coordinates of q points
  REAL(KIND = DP), ALLOCATABLE :: vnk_all(:,:,:)
  !! Electron velcity
  REAL(KIND = DP), ALLOCATABLE :: vnuq_all(:,:,:)
  !! phonon velcity
  REAL(KIND = DP) :: ef0
  !! Fermi energy
  COMPLEX(KIND = DP), ALLOCATABLE :: unk_all(:, :, :)
  !! Electron eigenvectors
  COMPLEX(KIND = DP), ALLOCATABLE :: Amnk_all(:, :, :, :)
  !! Berry connections
  COMPLEX(KIND = DP), ALLOCATABLE :: uf_all(:,:, :)
  !! Phonon eigenvectors
  REAL(KIND = DP), ALLOCATABLE :: epsil_ipa(:, :)
  !! Fermi energy
  INTEGER:: totq
  !! Total number of q-points within the fsthick window.
  INTEGER, ALLOCATABLE :: selecq(:)
  !! Array of selected q-points for e-ph scattering
  INTEGER, ALLOCATABLE :: bztoibz_dg(:)
  !! Map from k point to its irr-k point
  INTEGER, ALLOCATABLE :: s_bztoibz_dg(:)
  !! The symmetry operation mapping k and irr-k
  INTEGER :: nkirr_dg
  !! Number of irrducible k points on double grid
  REAL(KIND = DP), ALLOCATABLE :: ie_ph(:,:,:)
  !! Collision integral for e-ph interaction
  REAL(KIND = DP), ALLOCATABLE :: iph_e(:,:,:)
  !! Collision integral for e-ph interaction
  REAL(KIND = DP), ALLOCATABLE ::  iph_ph(:,:,:) 
  !! Collision integral for ph-ph interaction
  REAL(KIND = DP), ALLOCATABLE :: ie_e(:,:,:)
  !! Collision integral for el-el interaction
  COMPLEX(KIND = DP), ALLOCATABLE :: rho_e(:,:,:)
  !! Density matrix for electrons rho_e(nbndfst,nbndfst,nktotf)
  ! REAL (KIND = DP), ALLOCATABLE :: Qnu(:), Pnu(:)
  ! !! Cohrent phonon amplitude
  !! For future development of semiconductor Bloch equation
  REAL(KIND = DP), ALLOCATABLE :: felec(:,:)
  !! felec(nbndfst,nkfs) is the electron occupation at time t
  REAL(KIND = DP), ALLOCATABLE :: nphon(:,:) 
  !! nphon(nmodes,totq) is the phnonon occupation at time t
  REAL(KIND = DP), ALLOCATABLE :: nphon_pre(:,:) 
  !! Phonon distribution at the previous ph-ph time step
  INTEGER, ALLOCATABLE :: iq2nqfs(:, :)
  !! map of q-point from the full grid to the q-point for which
  !! both k and k\pm q falls into fermi shell
  INTEGER, ALLOCATABLE :: indx_mapk(:)
  !! map of k-point from nkfd to nkf
  INTEGER, ALLOCATABLE :: indx_mapq(:)
  !! map of q point from nqfd to nqf
  INTEGER :: nkfsden
  !! number of k-points in the fermi shell
  INTEGER, ALLOCATABLE :: ikfsdx(:)
  !! map of k index from the Fermi shell to the full grid.
  REAL(KIND = DP), PARAMETER :: fs2sec   = 1.0E-15_DP
  !! femtosecond to second
  REAL(KIND = DP), PARAMETER :: ps2sec   = 1.0E-12_DP
  !! picosecond to second
  INTEGER :: nstg
  !! Number of stages for Runge-Kutta mathods
  REAL(KIND = DP) :: b_tab(6)
  !! Butcher Tableaux parameters b_tab
  REAL(KIND = DP) :: c_tab(6)
  !! Butcher Tableaux parameters c_tab
  REAL(KIND = DP) :: a_tab(6,6)
  !! Butcher Tableaux parameters a_tab
  REAL(KIND = DP) :: tstart
  !! The starting time of simulation
  !--------------------------------------------------------------------------
  END MODULE tdbe_eph_common
  !--------------------------------------------------------------------------
  !
  !--------------------------------------------------------------------------
  MODULE tdbe_phph_common
  !--------------------------------------------------------------------------
  !!
  !! Global variables for time-dependent Boltzmann simulations
  !! and Semiconductor Bloch equations which will be implemented
  !! in the future.
  !!
  USE kinds, ONLY : DP
  !
  IMPLICIT NONE
  !
  SAVE
  !
  INTEGER :: nifc
  !! Number of anharmonic interatomic force constants
  INTEGER, ALLOCATABLE :: ind_atms(:, :)
  !! indexes of three atoms
  REAL(KIND = DP), ALLOCATABLE :: ifc3(:, :, :, :)
  !! 3rd force constants
  REAL(KIND = DP), ALLOCATABLE :: rvecs(:, :, :)
  !! vectors of  unit cells in crystal coordinates
  REAL(kind = DP),ALLOCATABLE :: psi2_p(:)
  !! phonon-phonon matrix elements multiplied by the delta function
  REAL(kind = DP),ALLOCATABLE :: psi2_m(:)
  !! phonon-phonon matrix elements multiplied by the delta function
  INTEGER, ALLOCATABLE :: ind_p(:,:)
  !! indexes of the ph-ph matrix elements
  INTEGER, ALLOCATABLE :: ind_m(:,:)
  !! indexes of the ph-ph matrix elements
  REAL(KIND = DP) :: dt_in_ps
  !! time step in picosecond
  INTEGER(KIND = 8) :: nind_p
  !! Number of ph-ph matrix elements
  INTEGER(KIND = 8) :: nind_m
  !! Number of ph-ph matrix elements
  REAL(KIND = DP),allocatable :: ph_lw(:,:)
  !! phonon linewidth 
  INTEGER, parameter :: ntph = 50000
  !! number of temperatures
  REAL(KIND = DP) :: phtemp(ntph) 
  !! array of temperatures
  REAL(KIND = DP) :: e_latt(ntph)
  !! array of lattice energy at temperatures in phtemp
  integer :: istart
  !! Index of temperature, used for deciding the average temperature of the 
  !! lattice
  !--------------------------------------------------------------------------
  END MODULE tdbe_phph_common
  !--------------------------------------------------------------------------
  !
  !
  !--------------------------------------------------------------------------
  MODULE tdbe_common
  !--------------------------------------------------------------------------
  !
  USE tdbe_eph_common
  USE tdbe_phph_common
  !
  !--------------------------------------------------------------------------
  END MODULE tdbe_common
  !--------------------------------------------------------------------------
