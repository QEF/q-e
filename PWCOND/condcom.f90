!
! Copyright (C) 2003 A. Smogunov 
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!    Common variables for conductance calculation  
!
module filnam_cond
  use parameters
  !
  !  Names of the files     
  !
  character(len=80) :: &
        band_file,       &
        tran_file,       &
        fil_loc        ! file with 2D eigenvectors and eigenvalues
end module filnam_cond

module zdir_cond
  use parameters
  !
  !  description of z direction
  !
  integer     :: &
     nrz,        &      ! total number of the slabs in z direction
     nrzp,       &      ! number of slabs per CPU
     dnslab,     &      ! additional slabs (nrz=nr3s+dnslab) 
     nz1                ! number of subslabs in the slab
  integer, allocatable :: &
     nkofz(:)           ! local CPU order of the slabs   
  real(kind=DP) :: &
     bdl1, bdl2,   &    ! boundaries for the left tip
     bdr1, bdr2,   &    !   --- || ---      right tip 
     zl                 ! the lenght of the unit cell in z direction
  real(kind=DP), allocatable ::  &
     z(:)               ! the unit cell division in z direction
end module zdir_cond      

module perp_cond
  use parameters
  !
  !  description of the direction perpendicular to z (XY)
  !
  integer      :: &
     nrx,        &      ! the number of mesh points in x direction
     nry,        &      !            --             in y direction
     ngper,       &     ! number of perpendicular G vectors
     ngpsh,       &     ! number of shells for G
     nkpts,       &
     n2d                ! dimension of reduced vector space in XY
  integer, allocatable :: &
     ninsh(:)           ! number of G in shell 
  real(kind=DP)    :: &
     sarea              ! the cross section  
  real(kind=DP), allocatable :: &
     xyk(:,:),        & ! coordinates of perpendicular k point
     wkpt(:),         &
     gper(:,:),       & ! coordinates of perpendicular G
     gnsh(:)            ! the norm of the shell
end module perp_cond       

module eigen_cond
  use parameters
  !
  !   Eigenvalue equation for local potential
  !
  complex(kind=DP), allocatable  :: &
     vppot(:,:),      & ! Fourier comp. of local potential in each slab
     psiper(:,:,:),   & ! eigenvectors in each slab
     newbg(:,:),      & ! basis of reduced set --> exp(G)
     zk(:,:)            ! the k for each eigenvalue (computed through zkr)
  real(kind=DP), allocatable ::  &
     zkr(:,:)           ! eigenvalues for first energy point
end module eigen_cond     

module control_cond
  use parameters
  !
  !   control of the run
  !
  integer      :: &
     ikind,       &      ! the kind of calculation
     nenergy,     &      ! number of energies computed             
     iofspin             ! spin index for calculation
  real(kind=DP)  :: &
     energy0,       &  ! initial energy
     denergy,       &  ! delta of energy
     eryd,          &  ! current energy
     ecut2d,        &  ! 2D cutoff 
     ewind,         &  ! the window above energy for 2D computation
     delgep,        &  ! infinitesimal for GEP 
     epsproj,       &  ! accuracy of n2d reduction
     cutplot           ! cutoff of Im(k) for CB plotting
  real(kind=DP), allocatable  :: &
     earr(:),         & ! energy array
     tran_tot(:)        ! transmission array
  logical        :: &
     lwrite_loc,    &   ! if .t. save eigenproblem result on fil_loc
     lread_loc,     &   ! if .t. read eigenproblem result from fil_loc
     llapack,       &   ! if .t. use LAPACK routine for GEP
     llocal             ! if .t. the local implementation
end module control_cond
         
module cross_cond
  use parameters
  !   
  !  description of nonlocal orbitals 
  !
  integer     ::     &
     nocrosl,        &  ! number of crossing orbitals for left tip
     nocrosr,        &  !           -- || --         for right tip 
     noinsl,         &  ! number of interior orbitals for left tip
     noinss,         &  !           -- || --         for scat. reg. 
     noinsr,         &  !           -- || --         for right tip 
     norb,           &  ! total number of orbitals
     norbf,          &  ! max number of needed orbitals
     norbs              ! tot. number of orbitals in scat. region
  integer, allocatable :: & 
     itnew(:),       &  ! the type of each orbital
     nbnew(:),       &  ! the nb of each orbital 
     ls(:),          &  ! the l of each orbital 
     mnew(:),        &  ! the m of each orbital
     cross(:,:),     &  ! 1 if the orbital crosses the slab
     natih(:,:)         ! na and ih of the orbital
  real(kind=DP), allocatable :: &
     rsph(:,:),      &  ! the radius of nonlocal sphere
     taunew(:,:),    &  ! center of each orbital
     zpseu(:,:,:)       ! coefficients of nonlocal pseudopotential
end module cross_cond

module scattnl_cond
  use parameters
  !
  !   The variables computed by scatter_forw
  !
  complex(kind=DP), allocatable :: &
     fun0(:,:),      &  ! local fun. on         left boundary
     fun1(:,:),      &  !      --              right boundary
     fund0(:,:),     &  ! local fun.' on        left boundary
     fund1(:,:),     &  !      --              right boundary
     funl0(:,:),     &  ! nonloc. fun. on       left boundary
     funl1(:,:),     &  !      --              right boundary
     fundl0(:,:),    &  ! nonlocal fun.' on     left boundary
     fundl1(:,:),    &  !      --               right boundary
     intw1(:,:),     &  ! integrals with beta-fun. of loc. fun.
     intw2(:,:)         !      --                   nonloc fun. 
end module scattnl_cond                        

module cb_cond
  use parameters
  !
  !   Some variables of CBS for the tips needed for matching
  !
  integer     ::     &
     nchanl,         &  ! number of prop. channels in the left tip
     nchanr             !           -- || --             right tip     
  complex(kind=DP), allocatable  :: &
     kvall(:),       &  ! k           for left tip 
     kfunl(:,:),     &  ! phi_k(z=d)  for left tip
     kfundl(:,:),    &  ! phi_k'(z=d) for left tip
     kintl(:,:),     &  ! integral of phi_k with beta-fun.
     kcoefl(:,:),    &  ! coeff. of phi_k over nonloc. fun.
     kvalr(:),       &  ! k           for right tip
     kfunr(:,:),     &  ! phi_k(z=0)  for right tip
     kfundr(:,:),    &  ! phi_k'(z=0) for right tip
     kintr(:,:),     &  ! integral of phi_k with beta-fun.
     kcoefr(:,:)        ! coeff. of phi_k over nonloc. fun.
end module cb_cond

module cond
    use filnam_cond
    use zdir_cond
    use perp_cond
    use eigen_cond
    use control_cond
    use cross_cond
    use scattnl_cond
    use cb_cond
end module cond               

