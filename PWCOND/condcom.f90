
!
! Copyright (C) 2003 A. Smogunov 
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
!
! ... Common variables for conductance calculation  
!
!
!
MODULE zdir_cond
  USE kinds
  !
  ! ... description of z direction
  !
  SAVE
  !
  INTEGER :: &
      nrz,             &!  total number of the slabs in z direction
      nrzp,            &!  number of slabs per CPU
      dnslab,          &!  additional slabs (nrz=nr3s+dnslab) 
      nz1               !  number of subslabs in the slab
  INTEGER, ALLOCATABLE :: &
      nkofz(:)          !  local CPU order of the slabs   
  REAL(KIND=DP) :: &
      bdl1, bdl2,      &!  boundaries for the left tip
      bdr1, bdr2,      &!    --- || ---      right tip 
      zl                !  the lenght of the unit cell in z direction
  REAL(KIND=DP), ALLOCATABLE ::  &
      z(:)              !  the unit cell division in z direction
  !
END MODULE zdir_cond      
!
!
MODULE perp_cond
  USE kinds
  !
  ! ... description of the direction perpendicular to z (XY)
  !
  SAVE
  !
  INTEGER :: &
       nrx,            &!  the number of mesh points in x direction
       nry,            &!             --             in y direction
       ngper,          &!  number of perpendicular G vectors
       ngpsh,          &!  number of shells for G
       nkpts,          &!  number of kpts in the perpendicular direction
       n2d              !  dimension of reduced vector space in XY
  INTEGER, ALLOCATABLE :: &
       ninsh(:)         !  number of G in shell 
  REAL(KIND=DP) :: &
       sarea            !  the cross section  
  REAL(KIND=DP), ALLOCATABLE :: &
       xyk(:,:),       &!  coordinates of perpendicular k point
       wkpt(:),        &
       gper(:,:),      &!  coordinates of perpendicular G
       gnsh(:)          !  the norm of the shell
  !
END MODULE perp_cond       
!
!
MODULE eigen_cond
  USE kinds
  !
  ! ... Eigenvalue equation for local potential
  !
  SAVE
  !
  COMPLEX(KIND=DP), ALLOCATABLE :: &
       vppot(:,:,:,:), &!  Fourier comp. of local potential in each slab
       psiper(:,:,:),  &!  eigenvectors in each slab
       newbg(:,:),     &!  basis of reduced set --> exp(G)
       zk(:,:)          !  the k for each eigenvalue (computed through zkr)
  REAL(KIND=DP), ALLOCATABLE ::  &
     zkr(:,:)           !  eigenvalues for first energy point
  !
END MODULE eigen_cond     
!
!
MODULE control_cond
  USE kinds
  !
  ! ... control of the run
  !
  SAVE
  !
  INTEGER :: &
      ikind,           &!  the kind of calculation
      nenergy,         &!  number of energies computed             
      iofspin           !  spin index for calculation
  REAL(KIND=DP) :: &
      energy0,         &!  initial energy
      denergy,         &!  delta of energy
      eryd,            &!  current energy
      ecut2d,          &!  2D cutoff 
      ewind,           &!  the window above energy for 2D computation
      delgep,          &!  infinitesimal for GEP 
      epsproj,         &!  accuracy of n2d reduction
      cutplot           !  cutoff of Im(k) for CB plotting
  REAL(KIND=DP), ALLOCATABLE :: &
      earr(:),         &!  energy array
      tran_tot(:)       !  transmission array
  LOGICAL :: &
      lwrite_loc,      &!  if .T. save eigenproblem result on fil_loc
      lread_loc,       &!  if .T. read eigenproblem result from fil_loc
      llapack,         &!  if .T. use LAPACK routine for GEP
      llocal            !  if .T. the local implementation
  !
END MODULE control_cond
!
!         
MODULE cross_cond
  USE kinds
  !   
  ! ... description of nonlocal orbitals 
  !
  SAVE
  !
  INTEGER :: &
      nocrosl,         &!  number of crossing orbitals for left tip
      nocrosr,         &!            -- || --         for right tip 
      noinsl,          &!  number of interior orbitals for left tip
      noinss,          &!            -- || --         for scat. reg. 
      noinsr,          &!            -- || --         for right tip 
      norb,            &!  total number of orbitals
      norbf,           &!  max number of needed orbitals
      norbs             !  tot. number of orbitals in scat. region
  INTEGER, ALLOCATABLE :: & 
      itnew(:),        &!  the type of each orbital
      nbnew(:),        &!  the nb of each orbital 
      ls(:),           &!  the l of each orbital 
      mnew(:),         &!  the m of each orbital
      cross(:,:),      &!  1 if the orbital crosses the slab
      natih(:,:)        !  na and ih of the orbital
  REAL(KIND=DP), ALLOCATABLE :: &
      rsph(:,:),       &!  the radius of nonlocal sphere
      taunew(:,:),     &!  center of each orbital
      zpseu(:,:,:)      !  coefficients of nonlocal pseudopotential
  COMPLEX(KIND=DP), ALLOCATABLE :: &
      zpseu_nc(:,:,:)   !  coefficients of nonlocal pseudopotential (nc case)
  !    
END MODULE cross_cond
!
!
MODULE scattnl_cond
  USE kinds
  !
  ! ... The variables computed by scatter_forw
  !
  SAVE
  !
  COMPLEX(KIND=DP), ALLOCATABLE :: &
       fun0(:,:),      &!  local fun. on         left boundary
       fun1(:,:),      &!       --              right boundary
       fund0(:,:),     &!  local fun.' on        left boundary
       fund1(:,:),     &!       --              right boundary
       funl0(:,:),     &!  nonloc. fun. on       left boundary
       funl1(:,:),     &!       --              right boundary
       fundl0(:,:),    &!  nonlocal fun.' on     left boundary
       fundl1(:,:),    &!       --              right boundary
       intw1(:,:),     &!  integrals with beta-fun. of loc. fun.
       intw2(:,:)       !       --                   nonloc fun. 
  !
END MODULE scattnl_cond                        
!
!
MODULE cb_cond
  USE kinds
  !
  ! ... Some variables of CBS for the tips needed for matching
  !
  SAVE
  !
  INTEGER :: &
       nchanl,         &!  number of prop. channels in the left tip
       nchanr           !            -- || --             right tip     
  COMPLEX(KIND=DP), ALLOCATABLE :: &
       kvall(:),       &!  k           for left tip 
       kfunl(:,:),     &!  phi_k(z=d)  for left tip
       kfundl(:,:),    &!  phi_k'(z=d) for left tip
       kintl(:,:),     &!  integral of phi_k with beta-fun.
       kcoefl(:,:),    &!  coeff. of phi_k over nonloc. fun.
       kvalr(:),       &!  k           for right tip
       kfunr(:,:),     &!  phi_k(z=0)  for right tip
       kfundr(:,:),    &!  phi_k'(z=0) for right tip
       kintr(:,:),     &!  integral of phi_k with beta-fun.
       kcoefr(:,:)      !  coeff. of phi_k over nonloc. fun.
  !
END MODULE cb_cond
!
!
MODULE cond
  USE zdir_cond
  USE perp_cond
  USE eigen_cond
  USE control_cond
  USE cross_cond
  USE scattnl_cond
  USE cb_cond
END MODULE cond               

