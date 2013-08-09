!
! Copyright (C) 2013 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!=----------------------------------------------------------------------------=!
!
MODULE environ_input
!
!=----------------------------------------------------------------------------=!
!  this module contains all variables in namelist &ENVIRON and related routines
!  performing initialization and broadcast - Written by Oliviero Andreussi
!=----------------------------------------------------------------------------=!
  !
  USE kinds,      ONLY : DP
  USE parameters, ONLY : nsx
  !
  IMPLICIT NONE
  !
  SAVE
  !
!=----------------------------------------------------------------------------=!
!  ENVIRON Namelist Input Parameters
!=----------------------------------------------------------------------------=!
!
! Global parameters
!
        INTEGER  :: verbose = 0
        ! verbosity  0: only prints summary of polarization charge calculation; 
        !    1: prints an extra file with details of iterative convergence;
        !    2: prints 3D cube files of physical properties
        REAL(DP) :: environ_thr = 1.d-1
        ! how early in scf should the corrective pot start being calculated
        CHARACTER( LEN = 256 ) :: environ_type = 'input'
        ! keyword to set up all the environment parameters at once to a specific set
        ! vacuum = all the flags are off (perm=1.d0, surf=0.0, pres=0.0)
        ! water = parameters optimized for water solutions in Andreussi et al. 
        !         J. Chem. Phys. 136, 064102 (perm=78, surf=50, pres=-0.35)
        ! input = do not use any predefined set, use paramters from input
!
! Switching function parameters
!
        REAL(DP) :: stype = 1
        ! type of switching functions used in the solvation models
        !    0: original Fattebert-Gygi
        !    1: ultrasoft dielectric function as defined in Andreussi et al.
        REAL(DP) :: rhomax = 0.005
        ! first parameter of the sw function, roughly corresponding 
        ! to the density threshold of the solvation model
        REAL(DP) :: rhomin = 0.0001
        ! second parameter of the sw function when stype=1
        REAL(DP) :: tbeta = 4.8
        ! second parameter of the sw function when stype=0
!
! Dielectric solvent parameters
!
        REAL(DP) :: env_static_permittivity = 1.D0
        ! static dielectric permittivity of the solvation model. If set equal
        ! to one (=vacuum) no dielectric effects
        CHARACTER( LEN = 256 ) :: eps_mode = 'electronic'
        !  eps_mode method for calculating the density that sets 
        !  the dielectric constant
        !  electronic = dielectric depends self-consist. on electronic density
        !  ionic = dielectric defined on a fictitious ionic density, generated
        !          as the sum of exponential functions centered on atomic 
        !          positions of width specified in input by solvationrad(ityp)
        !  full  = similar to electronic, but an extra density is added to 
        !          represent the core electrons and the nuclei. This extra 
        !          density is defined as the sum of gaussian functions centered
        !          on atomic positions of width equal to atomicspread(ityp)
        REAL(DP) :: solvationrad(nsx) = 3.D0
        ! solvationrad radius of the solvation shell for each species when the
        ! ionic dielectric function is adopted
        REAL(DP) :: atomicspread(nsx) = 0.5D0
        ! gaussian spreads of the atomic density of charge
        LOGICAL :: add_jellium = .false.
        ! depending on periodic boundary corrections, one may need to explicitly
        ! polarize the compensatinig jellium background
!
! Numerical differentiators paramters
!
        INTEGER  :: ifdtype = 1 
        ! type of numerical differentiator: 1=central differences, 
        ! 2=low-noise lanczos (m=2), 3=low-noise lanczos (m=4), 
        ! 4=smooth noise-robust (n=2), 5=smooth noise-robust (n=4)
        INTEGER  :: nfdpoint = 1
        ! number of points used in the numerical differentiator 
        ! N = 2*nfdpoint+1
!
! Iterative solver parameters
!
        CHARACTER( LEN=256 ) :: mixtype = 'linear'
        ! mixing method for iterative calculation of polarization charges
        ! 'linear', 'anderson', 'diis', 'broyden'
        REAL(DP) :: mixrhopol = 0.5D0
        ! mixing param to be used in iter calculation of polarization charges
        REAL(DP) :: tolrhopol = 1.D-10
        ! convergence threshold for polarization charges in iterative procedure
        INTEGER :: ndiis=1
        ! order of DIIS interpolation of iterative polarization charge
!
! Cavitation energy parameters
!
        REAL(DP) :: env_surface_tension = 0.D0
        ! solvent surface tension, if equal to zero no cavitation term 
        REAL(DP) :: delta = 0.00001D0
        ! finite difference parameter to compute the molecular surface
!
! PV energy parameters
!
        REAL(DP) :: env_pressure = 0.D0
        ! external pressure for PV energy, if equal to zero no pressure term 
!
! Ionic countercharge parameters
!
        REAL(DP) :: env_ioncc_concentration = 0.D0
        ! molar concentration of ionic countercharge (M=mol/L)
        REAL(DP) :: zion = 1.D0
        ! valence of ionic countercharge
        REAL(DP) :: rhopb = 0.0001D0
        ! density threshold for the onset of ionic countercharge
        REAL(DP) :: solvent_temperature = 300.D0
        ! temperature of the solution

        NAMELIST / environ /                                           &
             verbose, environ_thr, environ_type,                       &
             stype, rhomax, rhomin, tbeta,                             &
             env_static_permittivity, eps_mode,                        &
             solvationrad, atomicspread, add_jellium,                  &
             ifdtype, nfdpoint,                                        &
             mixtype, ndiis, mixrhopol, tolrhopol,                     &
             env_surface_tension, delta,                               &
             env_pressure,                                             &
             env_ioncc_concentration, zion, rhopb,                     &
             solvent_temperature

  CONTAINS
     !
     !=----------------------------------------------------------------------=!
     !
     !  Variables initialization for Namelist ENVIRON
     !
     !=----------------------------------------------------------------------=!
     !
     !-----------------------------------------------------------------------
     SUBROUTINE environ_defaults( prog )
       !-----------------------------------------------------------------------
       !
       IMPLICIT NONE
       !
       CHARACTER(LEN=2) :: prog   ! ... specify the calling program
       !
       !
       verbose      = 0
       environ_thr  = 1.D-1
       environ_type = 'input'
       !
       stype   = 1
       rhomax  = 0.005
       rhomin  = 0.0001
       tbeta   = 4.8
       !
       env_static_permittivity = 1.D0
       eps_mode        = 'electronic'
       solvationrad(:) = 3.D0
       atomicspread(:) = 0.5D0
       add_jellium = .false.
       !
       ifdtype  = 1
       nfdpoint = 2
       !
       mixtype   = 'linear'
       ndiis     = 1
       mixrhopol = 0.5
       tolrhopol = 1.D-10
       !
       env_surface_tension = 0.D0
       delta = 0.00001D0
       !
       env_pressure = 0.D0
       !
       env_ioncc_concentration = 0.0D0
       zion = 1.0D0
       rhopb = 0.0001D0
       solvent_temperature = 300.0D0
       !
       RETURN
       !
     END SUBROUTINE
     !
     !=----------------------------------------------------------------------=!
     !
     !  Broadcast variables values for Namelist ENVIRON
     !
     !=----------------------------------------------------------------------=!
     !
     !-----------------------------------------------------------------------
     SUBROUTINE environ_bcast()
       !-----------------------------------------------------------------------
       !
       USE io_global, ONLY : ionode_id
       USE mp,        ONLY : mp_bcast
       USE mp_images, ONLY : intra_image_comm
       !
       IMPLICIT NONE
       !
       CALL mp_bcast( verbose,                    ionode_id, intra_image_comm )
       CALL mp_bcast( environ_thr,                ionode_id, intra_image_comm )
       CALL mp_bcast( environ_type,               ionode_id, intra_image_comm )
       !
       CALL mp_bcast( stype,                      ionode_id, intra_image_comm )
       CALL mp_bcast( rhomax,                     ionode_id, intra_image_comm )
       CALL mp_bcast( rhomin,                     ionode_id, intra_image_comm )
       CALL mp_bcast( tbeta,                      ionode_id, intra_image_comm )
       !
       CALL mp_bcast( env_static_permittivity,    ionode_id, intra_image_comm )
       CALL mp_bcast( eps_mode,                   ionode_id, intra_image_comm )
       CALL mp_bcast( solvationrad,               ionode_id, intra_image_comm )
       CALL mp_bcast( atomicspread,               ionode_id, intra_image_comm )
       CALL mp_bcast( add_jellium,                ionode_id, intra_image_comm )
       !
       CALL mp_bcast( ifdtype,                    ionode_id, intra_image_comm )
       CALL mp_bcast( nfdpoint,                   ionode_id, intra_image_comm )
       !
       CALL mp_bcast( mixtype,                    ionode_id, intra_image_comm )
       CALL mp_bcast( ndiis,                      ionode_id, intra_image_comm )
       CALL mp_bcast( mixrhopol,                  ionode_id, intra_image_comm )
       CALL mp_bcast( tolrhopol,                  ionode_id, intra_image_comm )
       !
       CALL mp_bcast( env_surface_tension,        ionode_id, intra_image_comm )
       CALL mp_bcast( delta,                      ionode_id, intra_image_comm )
       !
       CALL mp_bcast( env_pressure,               ionode_id, intra_image_comm )
       !
       CALL mp_bcast( env_ioncc_concentration,    ionode_id, intra_image_comm )
       CALL mp_bcast( zion,                       ionode_id, intra_image_comm )
       CALL mp_bcast( rhopb,                      ionode_id, intra_image_comm )
       CALL mp_bcast( solvent_temperature,        ionode_id, intra_image_comm )
       !
      RETURN
       !
     END SUBROUTINE
     !
!=----------------------------------------------------------------------------=!
!
END MODULE environ_input
!
!=----------------------------------------------------------------------------=!
