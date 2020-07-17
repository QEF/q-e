! Copyright (C) 2001-2009 Quantum ESPRESSO group
! Copyright (C) 2019 Brian Kolb, Timo Thonhauser
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
! ----------------------------------------------------------------------


MODULE vdW_DF


! This module calculates the non-local correlation contribution to the
! energy and potential according to
!
!    M. Dion, H. Rydberg, E. Schroeder, D.C. Langreth, and
!    B.I. Lundqvist, Phys. Rev. Lett. 92, 246401 (2004).
!
! henceforth referred to as DION. Further information about the
! functional and its corresponding potential can be found in:
!
!    T. Thonhauser, V.R. Cooper, S. Li, A. Puzder, P. Hyldgaard,
!    and D.C. Langreth, Phys. Rev. B 76, 125112 (2007).
!
! The proper spin extension of vdW-DF, i.e. svdW-DF, is derived in
!
!    T. Thonhauser, S. Zuluaga, C.A. Arter, K. Berland, E. Schroder,
!    and P. Hyldgaard, Phys. Rev. Lett. 115, 136402 (2015).
!
! henceforth referred to as THONHAUSER.
!
!
! Two review articles show many of the vdW-DF applications:
!
!    D.C. Langreth et al., J. Phys.: Condens. Matter 21, 084203 (2009).
!
!    K. Berland et al., Rep. Prog. Phys. 78, 066501 (2015).
!
!
! The method implemented is based on the method of G. Roman-Perez and
! J.M. Soler described in:
!
!    G. Roman-Perez and J.M. Soler, Phys. Rev. Lett. 103, 096102 (2009).
!
! henceforth referred to as SOLER.
!
!
! xc_vdW_DF and xc_vdW_DF_spin are the driver routines for vdW-DF
! calculations and are called from Modules/funct.f90. The routines here
! set up the parallel run (if any) and carry out the calls necessary to
! calculate the non-local correlation contributions to the energy and
! potential.


USE kinds,             ONLY : dp
USE constants,         ONLY : pi, fpi, e2
USE mp,                ONLY : mp_sum, mp_barrier, mp_get, mp_size, mp_rank, mp_bcast
USE mp_images,         ONLY : intra_image_comm
USE mp_bands,          ONLY : intra_bgrp_comm
USE io_global,         ONLY : stdout, ionode
USE fft_base,          ONLY : dfftp
USE fft_interfaces,    ONLY : fwfft, invfft
USE control_flags,     ONLY : iverbosity, gamma_only
USE corr_lda,          ONLY : pw, pw_spin


! ----------------------------------------------------------------------
! No implicit variables

IMPLICIT NONE


! ----------------------------------------------------------------------
! By default everything is private

PRIVATE


! ----------------------------------------------------------------------
! Save all objects in this module

SAVE


! ----------------------------------------------------------------------
! Public functions

PUBLIC  :: xc_vdW_DF, xc_vdW_DF_spin, vdW_DF_stress,                   &
           vdW_DF_energy, vdW_DF_potential,                            &
           generate_kernel, interpolate_kernel,                        &
           initialize_spline_interpolation, spline_interpolation


! ----------------------------------------------------------------------
! Public variables

PUBLIC  :: inlc, vdW_DF_analysis, Nr_points, r_max, q_min, q_cut, Nqs, q_mesh


! ----------------------------------------------------------------------
! General variables

INTEGER                  :: inlc            = 1
! The non-local correlation

INTEGER                  :: vdW_DF_analysis = 0
! vdW-DF analysis tool as described in PRB 97, 085115 (2018)

REAL(DP), PARAMETER      :: epsr            = 1.0D-12
! A small number to cut off densities

INTEGER                  :: idx
! Indexing variable


! ----------------------------------------------------------------------
! Kernel specific parameters and variables

INTEGER, PARAMETER       :: Nr_points = 1024
! The number of radial points (also the number of k points) used in the
! formation of the kernel functions for each pair of q values.
! Increasing this value will help in case you get a run-time error
! saying that you are trying to use a k value that is larger than the
! largest tabulated k point since the largest k point will be 2*pi/r_max
! * Nr_points. Memory usage of the vdW_DF piece of PWSCF will increase
! roughly linearly with this variable.

REAL(DP), PARAMETER      :: r_max     = 100.0D0
! The value of the maximum radius to use for the real-space kernel
! functions for each pair of q values. The larger this value is the
! smaller the smallest k value will be since the smallest k point value
! is 2*pi/r_max. Be careful though, since this will also decrease the
! maximum k point value and the vdW_DF code will crash if it encounters
! a g-vector with a magnitude greater than 2*pi/r_max *Nr_points.

REAL(DP), PARAMETER      :: dr = r_max/Nr_points, dk = 2.0D0*pi/r_max
! Real space and k-space spacing of grid points.

REAL(DP), PARAMETER      :: q_min = 1.0D-5, q_cut = 5.0D0
! The maximum and minimum values of q. During a vdW run, values of q0
! found larger than q_cut will be saturated (SOLER equation 5) to q_cut.

INTEGER,  PARAMETER                 :: Nqs    = 20
REAL(DP), PARAMETER, DIMENSION(Nqs) :: q_mesh = (/  &
   q_min              , 0.0449420825586261D0, 0.0975593700991365D0, 0.159162633466142D0, &
   0.231286496836006D0, 0.315727667369529D0 , 0.414589693721418D0 , 0.530335368404141D0, &
   0.665848079422965D0, 0.824503639537924D0 , 1.010254382520950D0 , 1.227727621364570D0, &
   1.482340921174910D0, 1.780437058359530D0 , 2.129442028133640D0 , 2.538050036534580D0, &
   3.016440085356680D0, 3.576529545442460D0 , 4.232271035198720D0 , q_cut /)

! The above two parameters define the q mesh to be used in the vdW_DF
! code. These are perhaps the most important to have set correctly.
! Increasing the number of q points will DRAMATICALLY increase the
! memory usage of the vdW_DF code because the memory consumption depends
! quadratically on the number of q points in the mesh. Increasing the
! number of q points may increase accuracy of the vdW_DF code, although,
! in testing it was found to have little effect. The largest value of
! the q mesh is q_cut. All values of q0 (DION equation 11) larger than
! this value during a run will be saturated to this value using equation
! 5 of SOLER. In testing, increasing the value of q_cut was found to
! have little impact on the results, although it is possible that in
! some systems it may be more important. Always make sure that the
! variable Nqs is consistent with the number of q points that are
! actually in the variable q_mesh. Also, do not set any q value to 0.
! This will cause an infinity in the Fourier transform.

INTEGER,  PARAMETER      :: Nintegration_points = 256
! Number of integration points for real-space kernel generation (see
! DION equation 14). This is how many a's and b's there will be.

REAL(DP), PARAMETER      :: a_min = 0.0D0, a_max = 64.0D0
! Min/max values for the a and b integration in DION equation 14.

REAL(DP) :: kernel( 0:Nr_points, Nqs, Nqs ), d2phi_dk2( 0:Nr_points, Nqs, Nqs )
! Matrices holding the Fourier transformed kernel function  and its
! second derivative for each pair of q values. The ordering is
! kernel(k_point, q1_value, q2_value).

REAL(DP) :: W_ab( Nintegration_points, Nintegration_points )
! Defined in DION equation 16.

REAL(DP) :: a_points( Nintegration_points ), a_points2( Nintegration_points )
! The values of the "a" points (DION equation 14) and their squares.

CONTAINS








  ! ####################################################################
  !                           |             |
  !                           |  functions  |
  !                           |_____________|
  !
  ! Functions to be used in get_q0_on_grid, get_q0_on_grid_spin, and
  ! phi_value.

  FUNCTION Fs(s)

     IMPLICIT NONE
     REAL(DP) :: s, Fs, Z_ab = 0.0D0

     IF ( inlc == 1 .OR. inlc == 3 ) THEN
        Z_ab = -0.8491D0
     ELSE IF ( inlc == 2 .OR. inlc == 4 .OR. inlc == 5 ) THEN
        Z_ab = -1.887D0
     END IF

     Fs = 1.0D0 - Z_ab * s * s / 9.0D0

  END FUNCTION Fs




  FUNCTION dFs_ds(s)

     IMPLICIT NONE
     REAL(DP)             :: s, dFs_ds, Z_ab = 0.0D0
     REAL(DP), PARAMETER  :: prefac = -2.0D0/9.0D0

     IF ( inlc == 1 .OR. inlc == 3 ) THEN
        Z_ab = -0.8491D0
     ELSE IF ( inlc == 2 .OR. inlc == 4 .OR. inlc == 5 ) THEN
        Z_ab = -1.887D0
     END IF

     dFs_ds =  prefac * s * Z_ab

  END FUNCTION dFs_ds




  FUNCTION kF(rho)

     IMPLICIT NONE
     REAL(DP)             :: rho, kF
     REAL(DP), PARAMETER  :: ex = 1.0D0/3.0D0

     kF = ( 3.0D0 * pi * pi * rho )**ex

  END FUNCTION kF




  FUNCTION dkF_drho(rho)

     IMPLICIT NONE
     REAL(DP)             :: rho, dkF_drho
     REAL(DP), PARAMETER  :: prefac = 1.0D0/3.0D0

     dkF_drho = prefac * kF(rho) / rho

  END FUNCTION dkF_drho




  FUNCTION ds_drho(rho, s)

     IMPLICIT NONE
     REAL(DP) :: rho, s, ds_drho

     ds_drho = -s * ( dkF_drho(rho) / kF(rho) + 1.0D0 / rho )

  END FUNCTION ds_drho




  FUNCTION ds_dgradrho(rho)

     IMPLICIT NONE
     REAL(DP) :: rho, ds_dgradrho

     ds_dgradrho = 0.5D0 / (kF(rho) * rho)

  END FUNCTION ds_dgradrho




  FUNCTION dqx_drho(rho, s)

     IMPLICIT NONE
     REAL(DP) :: rho, s, dqx_drho

     dqx_drho = dkF_drho(rho) * Fs(s) + kF(rho) * dFs_ds(s) * ds_drho(rho, s)

  END FUNCTION dqx_drho




  FUNCTION h_function(y)

     IMPLICIT NONE
     REAL(DP)             :: y, y2, y4, h_function
     REAL(DP), PARAMETER  :: g1 = fpi/9.0D0                                     ! vdW-DF1/2
     REAL(DP), PARAMETER  :: a3 = 0.94950D0, g3 = 1.12D0, g32 = g3*g3           ! vdW-DF3-opt1
     REAL(DP), PARAMETER  :: a4 = 0.28248D0, g4 = 1.29D0, g42 = g4*g4           ! vdW-DF3-opt2
     REAL(DP), PARAMETER  :: a5 = 2.01059D0, b5 = 8.17471D0, g5 = 1.84981D0, &  ! vdW-DF-C6
                             AA = ( b5 + a5*(a5/2.0D0-g5) ) / ( 1.0D0+g5-a5 )   !


     y2 = y*y

     IF ( inlc == 1 .OR. inlc == 2 ) THEN

        h_function = 1.0D0 - EXP( -g1*y2 )
     
     ELSE IF ( inlc == 3 ) THEN

        y4 = y2*y2
        h_function = 1.0D0 - 1.0D0 / ( 1.0D0 + g3*y2 + g32*y4 + a3*y4*y4 )

     ELSE IF ( inlc == 4 ) THEN 

        y4 = y2*y2
        h_function = 1.0D0 - 1.0D0 / ( 1.0D0 + g4*y2 + g42*y4 + a4*y4*y4 )

     ELSE IF ( inlc == 5 ) THEN
     
        y4 = y2*y2
        h_function = 1.0D0 - ( 1.0D0 + ( (a5-g5)*y2 + AA*y4 ) / ( 1.0D0+AA*y2 ) ) * EXP( -a5*y2 ) 

     END IF

  END FUNCTION








  ! ####################################################################
  !                           |             |
  !                           |  XC_VDW_DF  |
  !                           |_____________|

  SUBROUTINE xc_vdW_DF (rho_valence, rho_core, etxc, vtxc, v)

  USE gvect,                 ONLY : ngm, g
  USE cell_base,             ONLY : omega, tpiba

  IMPLICIT NONE

  ! --------------------------------------------------------------------
  ! Local variables
  !                                               _
  REAL(DP), INTENT(IN)    :: rho_valence(:,:)    !
  REAL(DP), INTENT(IN)    :: rho_core(:)         !  PWSCF input variables
  REAL(DP), INTENT(INOUT) :: etxc, vtxc, v(:,:)  !_

  INTEGER :: i_grid, theta_i, i_proc             ! Indexing variables over grid points,
                                                 ! theta functions, and processors.

  REAL(DP) :: grid_cell_volume                   ! The volume of the unit cell per G-grid point.

  REAL(DP), ALLOCATABLE ::  q0(:)                ! The saturated value of q (equations 11 and 12
                                                 ! of DION). This saturation is that of
                                                 ! equation 5 in SOLER.

  REAL(DP), ALLOCATABLE :: grad_rho(:,:)         ! The gradient of the charge density. The
                                                 ! format is as follows:
                                                 ! grad_rho(cartesian_component,grid_point).

  REAL(DP), ALLOCATABLE :: potential(:)          ! The vdW contribution to the potential.

  REAL(DP), ALLOCATABLE :: dq0_drho(:)           ! The derivative of the saturated q0
                                                 ! (equation 5 of SOLER) with respect
                                                 ! to the charge density (see
                                                 ! get_q0_on_grid subroutine for details).

  REAL(DP), ALLOCATABLE :: dq0_dgradrho(:)       ! The derivative of the saturated q0
                                                 ! (equation 5 of SOLER) with respect
                                                 ! to the gradient of the charge density
                                                 ! (again, see get_q0_on_grid subroutine).

  COMPLEX(DP), ALLOCATABLE :: thetas(:,:)        ! These are the functions of equation 8 of
                                                 ! SOLER. They will be forward Fourier transformed
                                                 ! in place to get theta(k) and worked on in
                                                 ! place to get the u_alpha(r) of equation 11
                                                 ! in SOLER. They are formatted as follows:
                                                 ! thetas(grid_point, theta_i).

  REAL(DP) :: Ec_nl                              ! The non-local vdW contribution to the energy.

  REAL(DP), ALLOCATABLE :: total_rho(:)          ! This is the sum of the valence and core
                                                 ! charge. This just holds the piece assigned
                                                 ! to this processor.

  LOGICAL, SAVE :: first_iteration = .TRUE.      ! Whether this is the first time this
                                                 ! routine has been called.




  ! --------------------------------------------------------------------
  ! Write out the vdW-DF information and initialize the calculation.

  IF ( first_iteration ) THEN
     IF ( inlc > 5 ) CALL errore( 'xc_vdW_DF', 'inlc not implemented', 1 )
     CALL generate_kernel
     IF ( ionode ) CALL vdW_info
     first_iteration = .FALSE.
  END IF


  ! --------------------------------------------------------------------
  ! Allocate arrays. nnr is a PWSCF variable that holds the number of
  ! points assigned to a given processor.

  allocate( total_rho(dfftp%nnr), grad_rho(3,dfftp%nnr),                &
            potential(dfftp%nnr), thetas(dfftp%nnr, Nqs),               &
            q0(dfftp%nnr), dq0_drho(dfftp%nnr), dq0_dgradrho(dfftp%nnr) )


  ! --------------------------------------------------------------------
  ! Add together the valence and core charge densities to get the total
  ! charge density. Note that rho_core is not the true core density and
  ! it is only non-zero for pseudopotentials with non-local core
  ! corrections.

  total_rho = rho_valence(:,1) + rho_core(:)


  ! --------------------------------------------------------------------
  ! Here we calculate the gradient in reciprocal space using FFT.

  CALL fft_gradient_r2r (dfftp, total_rho, g, grad_rho)


  ! --------------------------------------------------------------------
  ! Find the value of q0 for all assigned grid points. q is defined in
  ! equations 11 and 12 of DION and q0 is the saturated version of q
  ! defined in equation 5 of SOLER. This routine also returns the
  ! derivatives of the q0s with respect to the charge-density and the
  ! gradient of the charge-density. These are needed for the potential
  ! calculated below. This routine also calculates the thetas.

  CALL get_q0_on_grid (total_rho, grad_rho, q0, dq0_drho, dq0_dgradrho, thetas)


  ! --------------------------------------------------------------------
  ! Carry out the integration in equation 8 of SOLER. This also turns
  ! the theta arrays into the precursor to the u_i(k) array which is
  ! inverse fourier transformed to get the u_i(r) functions of SOLER
  ! equation 11. Add the energy we find to the output variable etxc.

  CALL vdW_DF_energy (thetas, Ec_nl)
  etxc = etxc + Ec_nl

  IF ( iverbosity > 0 ) THEN
     CALL mp_sum(Ec_nl, intra_bgrp_comm)
     IF ( ionode ) THEN
        WRITE(stdout,'(/ / A)')       "     -----------------------------------------------"
        WRITE(stdout,'(A, F15.8, A)') "     Non-local corr. energy    =  ", Ec_nl, " Ry"
        WRITE(stdout,'(A /)')         "     -----------------------------------------------"
     END IF
  END IF


  ! --------------------------------------------------------------------
  ! Here we calculate the potential. This is calculated via equation 10
  ! of SOLER, using the u_i(r) calculated from quations 11 and 12 of
  ! SOLER. Each processor allocates the array to be the size of the full
  ! grid because, as can be seen in SOLER equation 10, processors need
  ! to access grid points outside their allocated regions. Begin by
  ! FFTing the u_i(k) to get the u_i(r) of SOLER equation 11.

  DO theta_i = 1, Nqs
     CALL invfft('Rho', thetas(:,theta_i), dfftp)
  END DO

  CALL vdW_DF_potential (q0, dq0_drho, dq0_dgradrho, grad_rho, thetas, potential)
  v(:,1) = v(:,1) + e2 * potential(:)


  ! --------------------------------------------------------------------
  ! The integral of rho(r)*potential(r) for the vtxc output variable.

  grid_cell_volume = omega/(dfftp%nr1*dfftp%nr2*dfftp%nr3)

  DO i_grid = 1, dfftp%nnr
     vtxc = vtxc + e2 * grid_cell_volume * rho_valence(i_grid,1) * potential(i_grid)
  END DO

  DEALLOCATE ( total_rho, grad_rho, potential, thetas, q0, dq0_drho, dq0_dgradrho )

  END SUBROUTINE xc_vdW_DF








  ! ####################################################################
  !                          |                  |
  !                          |  XC_VDW_DF_spin  |
  !                          |__________________|
  !
  ! This subroutine is as similar to xc_vdW_DF as possible, but handles
  ! the collinear nspin=2 case.

  SUBROUTINE xc_vdW_DF_spin (rho_valence, rho_core, etxc, vtxc, v)

  USE gvect,                 ONLY : ngm, g
  USE cell_base,             ONLY : omega, tpiba

  IMPLICIT NONE

  ! --------------------------------------------------------------------
  ! Local variables
  !                                              _
  REAL(DP), INTENT(IN) :: rho_valence(:,:)      !
  REAL(DP), INTENT(IN) :: rho_core(:)           ! PWSCF input variables.
  REAL(DP), INTENT(INOUT) :: etxc, vtxc, v(:,:) !_


  INTEGER :: i_grid, theta_i, i_proc            ! Indexing variables over grid points,
                                                ! theta functions, and processors, and a
                                                ! generic index.

  REAL(DP) :: grid_cell_volume                  ! The volume of the unit cell per G-grid point.

  REAL(DP), ALLOCATABLE :: q0(:)                ! The saturated value of q (equations 11 and 12
                                                ! of DION). This saturation is that of
                                                ! equation 5 in SOLER.

  REAL(DP), ALLOCATABLE :: grad_rho(:,:)        ! The gradient of the charge density. The
                                                ! format is as follows:
                                                ! grad_rho(cartesian_component,grid_point).
  REAL(DP), ALLOCATABLE :: grad_rho_up(:,:)     ! The gradient of the up charge density.
                                                ! Same format as grad_rho.
  REAL(DP), ALLOCATABLE :: grad_rho_down(:,:)   ! The gradient of the down charge density.
                                                ! Same format as grad_rho.

  REAL(DP), ALLOCATABLE :: potential_up(:)      ! The vdW contribution to the potential.
  REAL(DP), ALLOCATABLE :: potential_down(:)    ! The vdW contribution to the potential.

  REAL(DP), ALLOCATABLE :: dq0_drho_up(:)       ! The derivative of the saturated q0
  REAL(DP), ALLOCATABLE :: dq0_drho_down(:)     ! (equation 5 of SOLER) with respect
                                                ! to the charge density (see
                                                ! get_q0_on_grid subroutine for details).

  REAL(DP), ALLOCATABLE :: dq0_dgradrho_up(:)   ! The derivative of the saturated q0
  REAL(DP), ALLOCATABLE :: dq0_dgradrho_down(:) ! (equation 5 of SOLER) with respect
                                                ! to the gradient of the charge density
                                                ! (again, see get_q0_on_grid subroutine).

  COMPLEX(DP), ALLOCATABLE :: thetas(:,:)       ! These are the functions of equation 8 of
                                                ! SOLER. They will be forward Fourier transformed
                                                ! in place to get theta(k) and worked on in
                                                ! place to get the u_alpha(r) of equation 11
                                                ! in SOLER. They are formatted as follows:
                                                ! thetas(grid_point, theta_i).

  REAL(DP) :: Ec_nl                             ! The non-local vdW contribution to the energy.

  REAL(DP), ALLOCATABLE :: total_rho(:)         ! This is the sum of the valence (up and down)
                                                ! and core charge. This just holds the piece
                                                ! assigned to this processor.
  REAL(DP), ALLOCATABLE :: rho_up(:)            ! This is the just the up valence charge.
                                                ! This just holds the piece assigned
                                                ! to this processor.
  REAL(DP), ALLOCATABLE :: rho_down(:)          ! This is the just the down valence charge.
                                                ! This just holds the piece assigned
                                                ! to this processor.

  LOGICAL, SAVE :: first_iteration = .TRUE.     ! Whether this is the first time this
                                                ! routine has been called.




  ! --------------------------------------------------------------------
  ! Write out the vdW-DF information and initialize the calculation.

  IF ( first_iteration ) THEN
     IF ( inlc > 5 ) CALL errore( 'xc_vdW_DF_spin', 'inlc not implemented', 1 )
     CALL generate_kernel
     IF ( ionode ) CALL vdW_info
     first_iteration = .FALSE.
  END IF


  ! --------------------------------------------------------------------
  ! Allocate arrays. nnr is a PWSCF variable that holds the number of
  ! points assigned to a given processor.

  ALLOCATE( total_rho(dfftp%nnr), rho_up(dfftp%nnr), rho_down(dfftp%nnr),         &
     grad_rho(3,dfftp%nnr), grad_rho_up(3,dfftp%nnr), grad_rho_down(3,dfftp%nnr), &
     potential_up(dfftp%nnr), potential_down(dfftp%nnr), thetas(dfftp%nnr, Nqs),  &
     q0(dfftp%nnr), dq0_drho_up(dfftp%nnr), dq0_dgradrho_up(dfftp%nnr),           &
     dq0_drho_down(dfftp%nnr), dq0_dgradrho_down(dfftp%nnr) )


  ! --------------------------------------------------------------------
  ! Add together the valence and core charge densities to get the total
  ! charge density. Note that rho_core is not the true core density and
  ! it is only non-zero for pseudopotentials with non-local core
  ! corrections.

  rho_up    = ( rho_valence(:,1) + rho_valence(:,2) + rho_core(:) )*0.5D0
  rho_down  = ( rho_valence(:,1) - rho_valence(:,2) + rho_core(:) )*0.5D0
  total_rho = rho_up + rho_down

#if defined (__SPIN_BALANCED)
  rho_up   = total_rho*0.5D0
  rho_down = rho_up
  WRITE(stdout,'(/,/,"     Performing spin-balanced Ecnl calculation!")')
#endif


  ! --------------------------------------------------------------------
  ! Here we calculate the gradient in reciprocal space using FFT.

  CALL fft_gradient_r2r (dfftp, total_rho, g, grad_rho)
  CALL fft_gradient_r2r (dfftp, rho_up,    g, grad_rho_up)
  CALL fft_gradient_r2r (dfftp, rho_down,  g, grad_rho_down)


  ! --------------------------------------------------------------------
  ! Find the value of q0 for all assigned grid points. q is defined in
  ! equations 11 and 12 of DION and q0 is the saturated version of q
  ! defined in equation 5 of SOLER. In the spin case, q0 is defined by
  ! equation 8 (and text above that equation) of THONHAUSER. This
  ! routine also returns the derivatives of the q0s with respect to the
  ! charge-density and the gradient of the charge-density. These are
  ! needed for the potential calculated below.

  CALL get_q0_on_grid_spin (total_rho, rho_up, rho_down, grad_rho, &
       grad_rho_up, grad_rho_down, q0, dq0_drho_up, dq0_drho_down, &
       dq0_dgradrho_up, dq0_dgradrho_down, thetas)


  ! --------------------------------------------------------------------
  ! Carry out the integration in equation 8 of SOLER. This also turns
  ! the thetas array into the precursor to the u_i(k) array which is
  ! inverse fourier transformed to get the u_i(r) functions of SOLER
  ! equation 11. Add the energy we find to the output variable etxc.

  CALL vdW_DF_energy(thetas, Ec_nl)
  etxc = etxc + Ec_nl

  IF ( iverbosity > 0 ) THEN
     CALL mp_sum(Ec_nl, intra_bgrp_comm)
     IF (ionode) THEN
        WRITE(stdout,'(/ / A)')       "     -----------------------------------------------"
        WRITE(stdout,'(A, F15.8, A)') "     Non-local corr. energy    =  ", Ec_nl, " Ry"
        WRITE(stdout,'(A /)')         "     -----------------------------------------------"
     END IF
  END IF


  ! --------------------------------------------------------------------
  ! Here we calculate the potential. This is calculated via equation 10
  ! of SOLER, using the u_i(r) calculated from quations 11 and 12 of
  ! SOLER. Each processor allocates the array to be the size of the full
  ! grid because, as can be seen in SOLER equation 10, processors need
  ! to access grid points outside their allocated regions. Begin by
  ! FFTing the u_i(k) to get the u_i(r) of SOLER equation 11.

  DO theta_i = 1, Nqs
     CALL invfft('Rho', thetas(:,theta_i), dfftp)
  END DO

  CALL vdW_DF_potential (q0, dq0_drho_up  , dq0_dgradrho_up  , grad_rho_up  , thetas, potential_up  )
  CALL vdW_DF_potential (q0, dq0_drho_down, dq0_dgradrho_down, grad_rho_down, thetas, potential_down)

  v(:,1) = v(:,1) + e2 * potential_up  (:)
  v(:,2) = v(:,2) + e2 * potential_down(:)


  ! --------------------------------------------------------------------
  ! The integral of rho(r)*potential(r) for the vtxc output variable

  grid_cell_volume = omega/(dfftp%nr1*dfftp%nr2*dfftp%nr3)

  DO i_grid = 1, dfftp%nnr
     vtxc = vtxc + e2 * grid_cell_volume * (rho_valence(i_grid,1) +   &
            rho_valence(i_grid,2)) * 0.5_dp * potential_up  (i_grid)  &
                 + e2 * grid_cell_volume * (rho_valence(i_grid,1) -   &
            rho_valence(i_grid,2)) * 0.5_dp * potential_down(i_grid)
  END DO

  DEALLOCATE( total_rho, rho_up, rho_down, grad_rho, grad_rho_up, grad_rho_down, &
              potential_up, potential_down, thetas,                              &
              q0, dq0_drho_up, dq0_dgradrho_up, dq0_drho_down, dq0_dgradrho_down )

  END SUBROUTINE xc_vdW_DF_spin








  ! ####################################################################
  !                       |                  |
  !                       |  GET_Q0_ON_GRID  |
  !                       |__________________|
  !
  ! This routine first calculates the q value defined in (DION equations
  ! 11 and 12), then saturates it according to (SOLER equation 5). More
  ! specifically it calculates the following:
  !
  !     q0(ir) = q0 as defined above
  !     dq0_drho(ir) = total_rho * d q0 /d rho
  !     dq0_dgradrho = total_rho / |grad_rho| * d q0 / d |grad_rho|

  SUBROUTINE get_q0_on_grid (total_rho, grad_rho, q0, dq0_drho, dq0_dgradrho, thetas)

  IMPLICIT NONE

  REAL(DP), INTENT(IN)      :: total_rho(:), grad_rho(:,:)         ! Input variables needed.

  REAL(DP), INTENT(OUT)     :: q0(:), dq0_drho(:), dq0_dgradrho(:) ! Output variables that have been allocated
                                                                   ! outside this routine but will be set here.
  COMPLEX(DP), INTENT(INOUT):: thetas(:,:)                         ! The thetas from SOLER.

  INTEGER, PARAMETER        :: m_cut = 12                          ! How many terms to include in the sum
                                                                   ! of SOLER equation 5.

  REAL(DP)                  :: rho                                 ! Local variable for the density.
  REAL(DP)                  :: r_s                                 ! Wigner-Seitz radius.
  REAL(DP)                  :: s                                   ! Reduced gradient.
  REAL(DP)                  :: q
  REAL(DP)                  :: ec
  REAL(DP)                  :: dq0_dq                              ! The derivative of the saturated
                                                                   ! q0 with respect to q.

  INTEGER                   :: i_grid                              ! Indexing variable.




  ! --------------------------------------------------------------------
  ! Initialize q0-related arrays.

  q0(:)           = q_cut
  dq0_drho(:)     = 0.0D0
  dq0_dgradrho(:) = 0.0D0


  DO i_grid = 1, dfftp%nnr

     rho = total_rho(i_grid)


     ! -----------------------------------------------------------------
     ! This prevents numerical problems. If the charge density is
     ! negative (an unphysical situation), we simply treat it as very
     ! small. In that case, q0 will be very large and will be saturated.
     ! For a saturated q0 the derivative dq0_dq will be 0 so we set q0 =
     ! q_cut and dq0_drho = dq0_dgradrho = 0 and go on to the next
     ! point.

     IF ( rho < epsr ) CYCLE


     ! -----------------------------------------------------------------
     ! Calculate some intermediate values needed to find q.

     r_s = ( 3.0D0 / (4.0D0*pi*rho) )**(1.0D0/3.0D0)

     s   = SQRT( grad_rho(1,i_grid)**2 + grad_rho(2,i_grid)**2 + grad_rho(3,i_grid)**2 ) / &
           (2.0D0 * kF(rho) * rho )


     ! -----------------------------------------------------------------
     ! This is the q value defined in equations 11 and 12 of DION.
     ! Use pw() from flib/functionals.f90 to get qc = kf/eps_x * eps_c.
     !
     CALL pw(r_s, 1, ec, dq0_drho(i_grid))
     q = -4.0D0*pi/3.0D0 * ec + kF(rho) * Fs(s)


     ! -----------------------------------------------------------------
     ! Bring q into its proper bounds.

     CALL saturate_q ( q, q_cut, q0(i_grid), dq0_dq )
     IF (q0(i_grid) < q_min) q0(i_grid) = q_min


     ! -----------------------------------------------------------------
     ! Here we find derivatives. These are actually the density times
     ! the derivative of q0 with respect to rho and grad_rho. The
     ! density factor comes in since we are really differentiating
     ! theta = (rho)*P(q0) with respect to density (or its gradient)
     ! which will be
     !
     !    dtheta_drho = P(q0) + dP_dq0 * [rho * dq0_dq * dq_drho]
     !
     ! and
     !
     !    dtheta_dgrad_rho = dP_dq0 * [rho * dq0_dq * dq_dgrad_rho]
     !
     ! The parts in square brackets are what is calculated here. The
     ! dP_dq0 term will be interpolated later.

     dq0_drho(i_grid)     = dq0_dq * rho * ( -4.0D0*pi/3.0D0 * &
                            (dq0_drho(i_grid) - ec)/rho + dqx_drho(rho, s) )
     dq0_dgradrho(i_grid) = dq0_dq * rho * kF(rho) * dFs_ds(s) * ds_dgradrho(rho)

  END DO


  ! --------------------------------------------------------------------
  ! Here we calculate the theta functions of SOLER equation 8. These are
  ! defined as
  !
  !    rho * P_i(q0(rho, grad_rho))
  !
  ! where P_i is a polynomial that interpolates a Kroneker delta
  ! function at the point q_i (taken from the q_mesh) and q0 is the
  ! saturated version of q. q is defined in equations 11 and 12 of DION
  ! and the saturation proceedure is defined in equation 5 of SOLER.
  ! This is the biggest memory consumer in the method since the thetas
  ! array is (total # of FFT points)*Nqs complex numbers. In a parallel
  ! run, each processor will hold the values of all the theta functions
  ! on just the points assigned to it. thetas are stored in reciprocal
  ! space as theta_i(k) because this is the way they are used later for
  ! the convolution (equation 8 of SOLER). Start by interpolating the
  ! P_i polynomials defined in equation 3 in SOLER for the particular q0
  ! values we have.

  CALL spline_interpolation (q_mesh, q0, thetas)

  DO i_grid = 1, dfftp%nnr
     thetas(i_grid,:) = thetas(i_grid,:) * total_rho(i_grid)
  END DO

  DO idx = 1, Nqs
     CALL fwfft ('Rho', thetas(:,idx), dfftp)
  END DO

  END SUBROUTINE get_q0_on_grid








  ! ####################################################################
  !                       |                       |
  !                       |  GET_Q0_ON_GRID_spin  |
  !                       |_______________________|

  SUBROUTINE get_q0_on_grid_spin (total_rho, rho_up, rho_down, grad_rho, &
             grad_rho_up, grad_rho_down, q0, dq0_drho_up, dq0_drho_down, &
             dq0_dgradrho_up, dq0_dgradrho_down, thetas)

  IMPLICIT NONE

  REAL(DP),  INTENT(IN)      :: total_rho(:), grad_rho(:,:)              ! Input variables.
  REAL(DP),  INTENT(IN)      :: rho_up(:), grad_rho_up(:,:)              ! Input variables.
  REAL(DP),  INTENT(IN)      :: rho_down(:), grad_rho_down(:,:)          ! Input variables.

  REAL(DP),  INTENT(OUT)     :: q0(:), dq0_drho_up(:), dq0_drho_down(:)  ! Output variables.
  REAL(DP),  INTENT(OUT)     :: dq0_dgradrho_up(:), dq0_dgradrho_down(:) ! Output variables.
  COMPLEX(DP), INTENT(INOUT) :: thetas(:,:)                              ! The thetas from SOLER.

  REAL(DP)                   :: rho, up, down                            ! Local copy of densities.
  REAL(DP)                   :: zeta                                     ! Spin polarization.
  REAL(DP)                   :: r_s                                      ! Wigner-Seitz radius.
  REAL(DP)                   :: q, qc, qx, qx_up, qx_down                ! q for exchange and correlation.
  REAL(DP)                   :: q0x_up, q0x_down                         ! Saturated q values.
  REAL(DP)                   :: fac
  REAL(DP)                   :: ec, vc(2)
  REAL(DP)                   :: dq0_dq, dq0x_up_dq, dq0x_down_dq         ! Derivative of q0 w.r.t q.
  REAL(DP)                   :: dqc_drho_up, dqc_drho_down               ! Intermediate values.
  REAL(DP)                   :: dqx_drho_up, dqx_drho_down               ! Intermediate values.
  REAL(DP)                   :: s_up, s_down                             ! Reduced gradients.
  INTEGER                    :: i_grid                                   ! Indexing variable.
  LOGICAL                    :: calc_qx_up, calc_qx_down


  fac = 2.0D0**(-1.0D0/3.0D0)


  ! --------------------------------------------------------------------
  ! Initialize q0-related arrays.

  q0(:)                = q_cut
  dq0_drho_up(:)       = 0.0D0
  dq0_drho_down(:)     = 0.0D0
  dq0_dgradrho_up(:)   = 0.0D0
  dq0_dgradrho_down(:) = 0.0D0


  DO i_grid = 1, dfftp%nnr

     rho  = total_rho(i_grid)
     up   = rho_up(i_grid)
     down = rho_down(i_grid)


     ! -----------------------------------------------------------------
     ! This prevents numerical problems. If the charge density is
     ! negative (an unphysical situation), we simply treat it as very
     ! small. In that case, q0 will be very large and will be saturated.
     ! For a saturated q0 the derivative dq0_dq will be 0 so we set q0 =
     ! q_cut and dq0_drho = dq0_dgradrho = 0 and go on to the next
     ! point.

     IF ( rho < epsr ) CYCLE

     calc_qx_up   = .TRUE.
     calc_qx_down = .TRUE.

     IF ( up   < epsr/2.0D0 ) calc_qx_up   = .FALSE.
     IF ( down < epsr/2.0D0 ) calc_qx_down = .FALSE.


     ! -----------------------------------------------------------------
     ! The spin case is numerically even more tricky and we have to
     ! saturate each spin channel separately. Note that we are
     ! saturating at a higher value here, so that very large q values
     ! get saturated to exactly q_cut in the second, overall saturation.

     q0x_up        = 0.0D0
     q0x_down      = 0.0D0
     dqx_drho_up   = 0.0D0
     dqx_drho_down = 0.0D0


     IF ( calc_qx_up ) THEN
        s_up    = SQRT( grad_rho_up(1,i_grid)**2 + grad_rho_up(2,i_grid)**2 + &
                  grad_rho_up(3,i_grid)**2 ) / (2.0D0 * kF(up) * up)
        qx_up   = kF(2.0D0*up) * Fs(fac*s_up)
        CALL saturate_q (qx_up, 4.0D0*q_cut, q0x_up, dq0x_up_dq)
     END IF

     IF ( calc_qx_down ) THEN
        s_down  = SQRT( grad_rho_down(1,i_grid)**2 + grad_rho_down(2,i_grid)**2 + &
                  grad_rho_down(3,i_grid)**2) / (2.0D0 * kF(down) * down)
        qx_down = kF(2.0D0*down) * Fs(fac*s_down)
        CALL saturate_q (qx_down, 4.0D0*q_cut, q0x_down, dq0x_down_dq)
     END IF


     ! -----------------------------------------------------------------
     ! This is the q value defined in equations 11 and 12 of DION and
     ! equation 8 of THONHAUSER (also see text above that equation).

     r_s  = ( 3.0D0 / (4.0D0*pi*rho) )**(1.0D0/3.0D0)
     zeta = (up - down) / rho
     IF ( ABS(zeta) > 1.0D0 ) zeta = SIGN(1.0D0, zeta)
     call pw_spin( r_s, zeta, ec, vc(1), vc(2) )
     dqc_drho_up   = vc(1)
     dqc_drho_down = vc(2)

     qx = ( up * q0x_up + down * q0x_down ) / rho
     qc = -4.0D0*pi/3.0D0 * ec
     q  = qx + qc


     ! -----------------------------------------------------------------
     ! Bring q into its proper bounds.

     CALL saturate_q (q, q_cut, q0(i_grid), dq0_dq)
     IF (q0(i_grid) < q_min) q0(i_grid) = q_min


     ! -----------------------------------------------------------------
     ! Here we find derivatives. These are actually the density times
     ! the derivative of q0 with respect to rho and grad_rho. The
     ! density factor comes in since we are really differentiating
     ! theta = (rho)*P(q0) with respect to density (or its gradient)
     ! which will be
     !
     !    dtheta_drho = P(q0) + dP_dq0 * [rho * dq0_dq * dq_drho]
     !
     ! and
     !
     !    dtheta_dgrad_rho = dP_dq0 * [rho * dq0_dq * dq_dgrad_rho]
     !
     ! The parts in square brackets are what is calculated here. The
     ! dP_dq0 term will be interpolated later.

     IF ( calc_qx_up ) THEN
        dqx_drho_up   = 2.0D0*dq0x_up_dq*up*dqx_drho(2.0D0*up, fac*s_up) + q0x_up*down/rho
        dq0_dgradrho_up (i_grid) = 2.0D0 * dq0_dq * dq0x_up_dq * up * kF(2.0D0*up) * &
                        dFs_ds(fac*s_up) * ds_dgradrho(2.0D0*up)
     END IF

     IF ( calc_qx_down ) THEN
        dqx_drho_down = 2.0D0*dq0x_down_dq*down*dqx_drho(2.0D0*down, fac*s_down) + q0x_down*up/rho
        dq0_dgradrho_down(i_grid) = 2.0D0 * dq0_dq * dq0x_down_dq * down * kF(2.0D0*down) * &
                        dFs_ds(fac*s_down) * ds_dgradrho(2.0D0*down)
     END IF

     IF ( calc_qx_down ) dqx_drho_up   = dqx_drho_up   - q0x_down*down/rho
     IF ( calc_qx_up )   dqx_drho_down = dqx_drho_down - q0x_up  *up  /rho

     dqc_drho_up   = -4.0D0*pi/3.0D0 * (dqc_drho_up   - ec)
     dqc_drho_down = -4.0D0*pi/3.0D0 * (dqc_drho_down - ec)

     dq0_drho_up  (i_grid) = dq0_dq * (dqc_drho_up   + dqx_drho_up  )
     dq0_drho_down(i_grid) = dq0_dq * (dqc_drho_down + dqx_drho_down)

  END DO


  ! --------------------------------------------------------------------
  ! Here we calculate the theta functions of SOLER equation 8. These are
  ! defined as
  !
  !    rho * P_i(q0(rho, grad_rho))
  !
  ! where P_i is a polynomial that interpolates a Kroneker delta
  ! function at the point q_i (taken from the q_mesh) and q0 is the
  ! saturated version of q. q is defined in equations 11 and 12 of DION
  ! and the saturation proceedure is defined in equation 5 of SOLER.
  ! This is the biggest memory consumer in the method since the thetas
  ! array is (total # of FFT points)*Nqs complex numbers. In a parallel
  ! run, each processor will hold the values of all the theta functions
  ! on just the points assigned to it. thetas are stored in reciprocal
  ! space as theta_i(k) because this is the way they are used later for
  ! the convolution (equation 8 of SOLER). Start by interpolating the
  ! P_i polynomials defined in equation 3 in SOLER for the particular q0
  ! values we have.

  CALL spline_interpolation (q_mesh, q0, thetas)

  DO i_grid = 1, dfftp%nnr
     thetas(i_grid,:) = thetas(i_grid,:) * total_rho(i_grid)
  END DO

  DO idx = 1, Nqs
     CALL fwfft ('Rho', thetas(:,idx), dfftp)
  END Do

  END SUBROUTINE get_q0_on_grid_spin








  ! ####################################################################
  !                            |              |
  !                            |  saturate_q  |
  !                            |______________|

  SUBROUTINE saturate_q (q, q_cutoff, q0, dq0_dq)

  IMPLICIT NONE

  REAL(DP),  INTENT(IN)      :: q             ! Input q.
  REAL(DP),  INTENT(IN)      :: q_cutoff      ! Cutoff q.
  REAL(DP),  INTENT(OUT)     :: q0            ! Output saturated q.
  REAL(DP),  INTENT(OUT)     :: dq0_dq        ! Derivative of dq0/dq.

  REAL(DP)                   :: e_exp         ! Exponent.
  INTEGER,   PARAMETER       :: m_cut = 12    ! How many terms to include in
                                              ! the sum of SOLER equation 5.




  ! --------------------------------------------------------------------
  ! Here, we calculate q0 by saturating q according to equation 5 of
  ! SOLER. Also, we find the derivative dq0_dq needed for the
  ! derivatives dq0_drho and dq0_dgradrh0 discussed below.

  e_exp  = 0.0D0
  dq0_dq = 0.0D0

  DO idx = 1, m_cut
     e_exp  = e_exp + (q/q_cutoff)**idx/idx
     dq0_dq = dq0_dq + (q/q_cutoff)**(idx-1)
  END Do

  q0     = q_cutoff*(1.0D0 - EXP(-e_exp))
  dq0_dq = dq0_dq * EXP(-e_exp)

  END SUBROUTINE saturate_q








  ! ####################################################################
  !                          |               |
  !                          | vdW_DF_energy |
  !                          |_______________|
  !
  ! This routine carries out the integration of equation 8 of SOLER.  It
  ! returns the non-local exchange-correlation energy and the u_alpha(k)
  ! arrays used to find the u_alpha(r) arrays via equations 11 and 12 in
  ! SOLER.

  SUBROUTINE vdW_DF_energy (thetas, vdW_xc_energy)

  USE gvect,           ONLY : gg, ngm, igtongl, gl, ngl, gstart
  USE cell_base,       ONLY : tpiba, omega

  IMPLICIT NONE

  COMPLEX(DP), INTENT(INOUT) :: thetas(:,:)            ! On input this variable holds the theta
                                                       ! functions (equation 8, SOLER) in the format
                                                       ! thetas(grid_point, theta_i). On output
                                                       ! this array holds u_alpha(k) =
                                                       ! Sum_j[theta_beta(k)phi_alpha_beta(k)].

  REAL(DP), INTENT(OUT) :: vdW_xc_energy               ! The non-local correlation energy.

  REAL(DP), ALLOCATABLE :: kernel_of_k(:,:)            ! This array will hold the interpolated kernel
                                                       ! values for each pair of q values in the q_mesh.

  REAL(DP)    :: g                                     ! The magnitude of the current g vector.
  INTEGER     :: last_g                                ! The shell number of the last g vector.


  INTEGER     :: g_i, q1_i, q2_i, i_grid               ! Index variables.

  COMPLEX(DP) :: theta(Nqs), thetam(Nqs), theta_g(Nqs) ! Temporary storage arrays used since we
                                                       ! are overwriting the thetas array here.
  REAL(DP)    :: G0_term, G_multiplier

  COMPLEX(DP), ALLOCATABLE :: u_vdw(:,:)               ! Temporary array holding u_alpha(k).




  ALLOCATE ( u_vdW(dfftp%nnr,Nqs), kernel_of_k(Nqs, Nqs) )
  vdW_xc_energy = 0.0D0
  u_vdW(:,:)    = CMPLX(0.0_DP, 0.0_DP, kind=dp)


  ! --------------------------------------------------------------------
  ! Loop over PWSCF's array of magnitude-sorted g-vector shells. For
  ! each shell, interpolate the kernel at this magnitude of g, then find
  ! all points on the shell and carry out the integration over those
  ! points. The PWSCF variables used here are ngm = number of g-vectors
  ! on this processor, nl = an array that gives the indices into the FFT
  ! grid for a particular g vector, igtongl = an array that gives the
  ! index of which shell a particular g vector is in, gl = an array that
  ! gives the magnitude of the g vectors for each shell. In essence, we
  ! are forming the reciprocal-space u(k) functions of SOLER equation
  ! 11. These are kept in thetas array. Here we should use gstart,ngm
  ! but all the cases are handled by conditionals inside the loop

  G_multiplier = 1.0D0
  IF ( gamma_only ) G_multiplier = 2.0D0

  last_g = -1

  DO g_i = 1, ngm

     IF ( igtongl(g_i) .NE. last_g) THEN
        g = SQRT(gl(igtongl(g_i))) * tpiba
        CALL interpolate_kernel(g, kernel_of_k)
        last_g = igtongl(g_i)
     END IF

     theta = thetas(dfftp%nl(g_i),:)

     DO q2_i = 1, Nqs
        DO q1_i = 1, Nqs
           u_vdW(dfftp%nl(g_i),q2_i)  = u_vdW(dfftp%nl(g_i),q2_i) + kernel_of_k(q2_i,q1_i)*theta(q1_i)
        END DO
        vdW_xc_energy = vdW_xc_energy + G_multiplier * (u_vdW(dfftp%nl(g_i),q2_i)*CONJG(theta(q2_i)))
     END DO

     IF (g_i < gstart ) vdW_xc_energy = vdW_xc_energy / G_multiplier

  END DO

  IF ( gamma_only ) u_vdW(dfftp%nlm(:),:) = CONJG( u_vdW(dfftp%nl(:),:) )


  ! --------------------------------------------------------------------
  ! Apply scaling factors. The e2 comes from PWSCF's choice of units.
  ! This should be 0.5 * e2 * vdW_xc_energy * (2pi)^3/omega * (omega)^2,
  ! with the (2pi)^3/omega being the volume element for the integral
  ! (the volume of the reciprocal unit cell) and the 2 factors of omega
  ! being used to cancel the factor of 1/omega PWSCF puts on forward
  ! FFTs of the 2 theta factors. 1 omega cancels and the (2pi)^3
  ! cancels because there should be a factor of 1/(2pi)^3 on the radial
  ! Fourier transform of phi that was left out to cancel with this
  ! factor.

  vdW_xc_energy = 0.5D0 * e2 * omega * vdW_xc_energy

  thetas(:,:) = u_vdW(:,:)
  DEALLOCATE ( u_vdW, kernel_of_k )

  END SUBROUTINE vdW_DF_energy








  ! ####################################################################
  !                        |                   |
  !                        |  vdW_DF_potential |
  !                        |___________________|
  !
  ! This routine finds the non-local correlation contribution to the
  ! potential (i.e. the derivative of the non-local piece of the energy
  ! with respect to density) given in SOLER equation 10. The u_alpha(k)
  ! functions were found while calculating the energy. They are passed
  ! in as the matrix u_vdW. Most of the required derivatives were
  ! calculated in the "get_q0_on_grid" routine, but the derivative of
  ! the interpolation polynomials, P_alpha(q), (SOLER equation 3) with
  ! respect to q is interpolated here, along with the polynomials
  ! themselves.

  SUBROUTINE vdW_DF_potential (q0, dq0_drho, dq0_dgradrho, grad_rho, u_vdW, potential)

  USE gvect,               ONLY : g
  USE cell_base,           ONLY : alat, tpiba

  IMPLICIT NONE

  REAL(DP), INTENT(IN) ::  q0(:), grad_rho(:,:)       ! Input arrays holding the value of q0 for
                                                      ! all points assigned to this processor and
                                                      ! the gradient of the charge density for
                                                      ! points assigned to this processor.

  REAL(DP), INTENT(IN) :: dq0_drho(:), dq0_dgradrho(:)! The derivative of q0 with respect to the
                                                      ! charge density and gradient of the charge
                                                      ! density (almost). See comments in the
                                                      ! get_q0_on_grid subroutine above.

  COMPLEX(DP), INTENT(IN) :: u_vdW(:,:)               ! The functions u_alpha(r) obtained by
                                                      ! inverse transforming the functions
                                                      ! u_alph(k). See equations 11 and 12 in SOLER.

  REAL(DP), INTENT(INOUT) :: potential(:)             ! The non-local correlation potential for
                                                      ! points on the grid over the whole cell (not
                                                      ! just those assigned to this processor).

  REAL(DP), ALLOCATABLE, SAVE :: d2y_dx2(:,:)         ! Second derivatives of P_alpha polynomials
                                                      ! for interpolation.

  INTEGER :: i_grid, P_i,icar                         ! Index variables.

  INTEGER :: q_low, q_hi, q                           ! Variables to find the bin in the q_mesh that
                                                      ! a particular q0 belongs to (for interpolation).
  REAL(DP) :: dq, a, b, c, d, e, f                    ! Intermediate variables used in the
                                                      ! interpolation of the polynomials.

  REAL(DP) :: y(Nqs), dP_dq0, P                       ! The y values for a given polynomial (all 0
                                                      ! exept for element i of P_i) The derivative
                                                      ! of P at a given q0 and the value of P at a
                                                      ! given q0. Both of these are interpolated
                                                      ! below.

  REAL(DP) :: gradient2                               ! Squared gradient.

  REAL(DP)   , ALLOCATABLE ::h_prefactor(:)
  COMPLEX(DP), ALLOCATABLE ::h(:)




  ALLOCATE ( h_prefactor(dfftp%nnr), h(dfftp%nnr) )

  potential     = 0.0D0
  h_prefactor   = 0.0D0


  ! --------------------------------------------------------------------
  ! Get the second derivatives of the P_i functions for interpolation.
  ! We have already calculated this once but it is very fast and it's
  ! just as easy to calculate it again.

  IF (.NOT. ALLOCATED( d2y_dx2) ) THEN

     ALLOCATE( d2y_dx2(Nqs, Nqs) )
     CALL initialize_spline_interpolation ( q_mesh, d2y_dx2(:,:) )

  end if


  DO i_grid = 1, dfftp%nnr

     q_low = 1
     q_hi = Nqs


     ! -----------------------------------------------------------------
     ! Figure out which bin our value of q0 is in in the q_mesh.

     DO WHILE ( (q_hi - q_low) > 1)

        q = INT((q_hi + q_low)/2)

        IF (q_mesh(q) > q0(i_grid)) THEN
           q_hi = q
        ELSE
           q_low = q
        END IF

     END DO

     IF ( q_hi == q_low ) CALL errore('vdW_DF_potential','qhi == qlow',1)

     dq = q_mesh(q_hi) - q_mesh(q_low)

     a = (q_mesh(q_hi) - q0(i_grid))/dq
     b = (q0(i_grid) - q_mesh(q_low))/dq
     c = (a**3 - a)*dq**2/6.0D0
     d = (b**3 - b)*dq**2/6.0D0
     e = (3.0D0*a**2 - 1.0D0)*dq/6.0D0
     f = (3.0D0*b**2 - 1.0D0)*dq/6.0D0

     DO P_i = 1, Nqs
        y = 0.0D0
        y(P_i) = 1.0D0

        P      = a*y(q_low) + b*y(q_hi)  + c*d2y_dx2(P_i,q_low) + d*d2y_dx2(P_i,q_hi)
        dP_dq0 = (y(q_hi) - y(q_low))/dq - e*d2y_dx2(P_i,q_low) + f*d2y_dx2(P_i,q_hi)


        ! --------------------------------------------------------------
        ! The first term in equation 10 of SOLER.
        potential(i_grid) = potential(i_grid) + u_vdW(i_grid,P_i)* (P + dP_dq0 * dq0_drho(i_grid))
        IF (q0(i_grid) .NE. q_mesh(Nqs)) THEN
           h_prefactor(i_grid) = h_prefactor(i_grid) + u_vdW(i_grid,P_i)*dP_dq0*dq0_dgradrho(i_grid)
        END IF

     END DO

  END DO

  DO icar = 1,3

     h(:) = CMPLX( h_prefactor(:) * grad_rho(icar,:), 0.0_DP, kind=dp )

     DO i_grid = 1, dfftp%nnr
        gradient2 = grad_rho(1,i_grid)**2 + grad_rho(2,i_grid)**2 + grad_rho(3,i_grid)**2
        IF ( gradient2 > 0.0D0 ) h(i_grid) = h(i_grid) / SQRT( gradient2 )
     END DO

     CALL fwfft ('Rho', h, dfftp)
     h(dfftp%nl(:)) = CMPLX(0.0_DP,1.0_DP, kind=dp) * tpiba * g(icar,:) * h(dfftp%nl(:))
     IF (gamma_only) h(dfftp%nlm(:)) = CONJG(h(dfftp%nl(:)))
     CALL invfft ('Rho', h, dfftp)
     potential(:) = potential(:) - REAL(h(:))

  END Do

  DEALLOCATE ( h_prefactor, h )

  END SUBROUTINE vdW_DF_potential








  ! ####################################################################
  !                       |                        |
  !                       |  SPLINE_INTERPOLATION  |
  !                       |________________________|
  !
  ! This routine is modeled after an algorithm from NUMERICAL_RECIPES It
  ! was adapted for Fortran, of course and for the problem at hand, in
  ! that it finds the bin a particular x value is in and then loops over
  ! all the P_i functions so we only have to find the bin once.

  SUBROUTINE spline_interpolation (x, evaluation_points, values)

  IMPLICIT NONE

  REAL(DP), INTENT(IN) :: x(:), evaluation_points(:)     ! Input variables. The x values used to
                                                         ! form the interpolation (q_mesh in this
                                                         ! case) and the values of q0 for which we
                                                         ! are interpolating the function.

  COMPLEX(DP), INTENT(INOUT) :: values(:,:)              ! An output array (allocated outside this
                                                         ! routine) that stores the interpolated
                                                         ! values of the P_i (SOLER equation 3)
                                                         ! polynomials. The format is
                                                         ! values(grid_point, P_i).

  INTEGER :: Ngrid_points, Nx                            ! Total number of grid points to evaluate
                                                         ! and input x points.

  REAL(DP), ALLOCATABLE, SAVE :: d2y_dx2(:,:)            ! The second derivatives required to do
                                                         ! the interpolation.

  INTEGER :: i_grid, lower_bound, upper_bound, P_i       ! Some indexing variables.

  REAL(DP), ALLOCATABLE :: y(:)                          ! Temporary variables needed for the
  REAL(DP) :: a, b, c, d, dx                             ! interpolation.




  Nx = size(x)
  Ngrid_points = size(evaluation_points)


  ! --------------------------------------------------------------------
  ! Allocate the temporary array.

  ALLOCATE( y(Nx) )


  ! --------------------------------------------------------------------
  ! If this is the first time this routine has been called we need to
  ! get the second derivatives (d2y_dx2) required to perform the
  ! interpolations. So we allocate the array and call
  ! initialize_spline_interpolation to get d2y_dx2.

  IF (.NOT. ALLOCATED(d2y_dx2) ) THEN

     ALLOCATE( d2y_dx2(Nx,Nx) )
     CALL initialize_spline_interpolation(x, d2y_dx2)

  END IF

  DO i_grid=1, Ngrid_points

     lower_bound = 1
     upper_bound = Nx

     DO WHILE ( (upper_bound - lower_bound) > 1 )

        idx = (upper_bound+lower_bound) / 2

        IF ( evaluation_points(i_grid) > x(idx) ) THEN
           lower_bound = idx
        ELSE
           upper_bound = idx
        END IF

     END DO

     dx = x(upper_bound)-x(lower_bound)

     a = (x(upper_bound) - evaluation_points(i_grid))/dx
     b = (evaluation_points(i_grid) - x(lower_bound))/dx
     c = ((a**3-a)*dx**2)/6.0D0
     d = ((b**3-b)*dx**2)/6.0D0

     DO P_i = 1, Nx

        y = 0
        y(P_i) = 1

        values(i_grid, P_i) = a*y(lower_bound) + b*y(upper_bound) &
             + (c*d2y_dx2(P_i,lower_bound) + d*d2y_dx2(P_i, upper_bound))
     END DO

  END DO

  DEALLOCATE( y )

  END SUBROUTINE spline_interpolation








  ! ####################################################################
  !                  |                                   |
  !                  |  INITIALIZE_SPLINE_INTERPOLATION  |
  !                  |___________________________________|
  !
  ! This routine is modeled after an algorithm from NUMERICAL_RECIPES It
  ! was adapted for Fortran and for the problem at hand.

  SUBROUTINE initialize_spline_interpolation (x, d2y_dx2)

  IMPLICIT NONE

  REAL(DP), INTENT(IN)    :: x(:)               ! The input abscissa values.
  REAL(DP), INTENT(INOUT) :: d2y_dx2(:,:)       ! The output array (allocated outside this routine)
                                                ! that holds the second derivatives required for
                                                ! interpolating the function.

  INTEGER :: Nx, P_i                            ! The total number of x points and some indexing
                                                ! variables.

  REAL(DP), ALLOCATABLE :: temp_array(:), y(:)  ! Some temporary arrays required. y is the array
                                                ! that holds the funcion values (all either 0 or
                                                ! 1 here).

  REAL(DP) :: temp1, temp2                      ! Some temporary variables required.




  Nx = SIZE(x)

  ALLOCATE( temp_array(Nx), y(Nx) )

  DO P_i=1, Nx

     ! -----------------------------------------------------------------
     ! In the Soler method, the polynomicals that are interpolated are Kroneker
     ! delta funcions at a particular q point. So, we set all y values to 0
     ! except the one corresponding to the particular function P_i.

     y = 0.0D0
     y(P_i) = 1.0D0

     d2y_dx2(P_i,1) = 0.0D0
     temp_array(1) = 0.0D0

     DO idx = 2, Nx-1

        temp1 = (x(idx)-x(idx-1))/(x(idx+1)-x(idx-1))
        temp2 = temp1 * d2y_dx2(P_i,idx-1) + 2.0D0
        d2y_dx2(P_i,idx) = (temp1-1.0D0)/temp2

        temp_array(idx) = (y(idx+1)-y(idx))/(x(idx+1)-x(idx)) &
             - (y(idx)-y(idx-1))/(x(idx)-x(idx-1))
        temp_array(idx) = (6.0D0*temp_array(idx)/(x(idx+1)-x(idx-1)) &
             - temp1*temp_array(idx-1))/temp2

     END DO

     d2y_dx2(P_i,Nx) = 0.0D0

     DO idx=Nx-1, 1, -1

        d2y_dx2(P_i,idx) = d2y_dx2(P_i,idx) * d2y_dx2(P_i,idx+1) + temp_array(idx)

     END DO

  END DO

  DEALLOCATE( temp_array, y )

  END SUBROUTINE initialize_spline_interpolation








  ! ####################################################################
  !                          |                    |
  !                          | INTERPOLATE_KERNEL |
  !                          |____________________|
  !
  ! This routine is modeled after an algorithm from NUMERICAL_RECIPES
  ! Adapted for Fortran and the problem at hand. This function is used
  ! to find the Phi_alpha_beta needed for equations
  ! 8 and 11 of SOLER.

  SUBROUTINE interpolate_kernel (k, kernel_of_k)

  IMPLICIT NONE

  REAL(DP), INTENT(IN)    :: k                 ! Input value, the magnitude of the g-vector
                                               ! for the current point.

  REAL(DP), INTENT(INOUT) :: kernel_of_k(:,:)  ! An output array (allocated outside this routine)
                                               ! that holds the interpolated value of the kernel
                                               ! for each pair of q points (i.e. the phi_alpha_beta
                                               ! of the Soler method.

  INTEGER  :: q1_i, q2_i, k_i                  ! Indexing variables.

  REAL(DP) :: A, B, C, D                       ! Intermediate values for the interpolation.




  ! --------------------------------------------------------------------
  ! Check to make sure that the kernel table we have is capable of
  ! dealing with this value of k. If k is larger than
  ! Nr_points*2*pi/r_max then we can't perform the interpolation. In
  ! that case, a kernel file should be generated with a larger number of
  ! radial points.

  IF ( k >= Nr_points*dk ) THEN

     WRITE(*,'(A,F10.5,A,F10.5)') "k =  ", k, "     k_max =  ", Nr_points*dk
     CALL errore('interpolate kernel', 'k value requested is out of range',1)

  END IF

  kernel_of_k = 0.0D0


  ! --------------------------------------------------------------------
  ! This integer division figures out which bin k is in since the kernel
  ! is set on a uniform grid.

  k_i = INT(k/dk)


  ! --------------------------------------------------------------------
  ! Test to see if we are trying to interpolate a k that is one of the
  ! actual function points we have. The value is just the value of the
  ! function in that case.

  IF ( MOD(k,dk) == 0 ) THEN

     DO q1_i = 1, Nqs
        DO q2_i = 1, q1_i

           kernel_of_k(q1_i, q2_i) = kernel(k_i,q1_i, q2_i)
           kernel_of_k(q2_i, q1_i) = kernel(k_i,q2_i, q1_i)

        END DO
     END DO

     RETURN

  END IF


  ! --------------------------------------------------------------------
  ! If we are not on a function point then we carry out the
  ! interpolation.

  A = (dk*(k_i+1.0D0) - k)/dk
  B = (k - dk*k_i)/dk
  C = (A**3-A)*dk**2/6.0D0
  D = (B**3-B)*dk**2/6.0D0

  DO q1_i = 1, Nqs
     DO q2_i = 1, q1_i

        kernel_of_k(q1_i, q2_i) = A*kernel(k_i, q1_i, q2_i) + B*kernel(k_i+1, q1_i, q2_i) &
                        + (C*d2phi_dk2(k_i, q1_i, q2_i) + D*d2phi_dk2(k_i+1, q1_i, q2_i))

        kernel_of_k(q2_i, q1_i) = kernel_of_k(q1_i, q2_i)

     END DO
  END DO

  END SUBROUTINE interpolate_kernel








  ! ####################################################################
  !                         |                 |
  !                         |  VDW_DF_STRESS  |
  !                         |_________________|

  SUBROUTINE vdW_DF_stress (rho_valence, rho_core, nspin, sigma)

  use gvect,           ONLY : ngm, g
  USE cell_base,       ONLY : tpiba

  IMPLICIT NONE

  REAL(DP), INTENT(IN)     :: rho_valence(:)         !
  REAL(dp), INTENT(IN)     :: rho_core(:)            ! Input variables.
  INTEGER,  INTENT(IN)     :: nspin                  !
  REAL(dp), INTENT(INOUT)  :: sigma(3,3)             !

  REAL(DP), ALLOCATABLE    :: grad_rho(:,:)          !
  REAL(DP), ALLOCATABLE    :: total_rho(:)           ! Rho values.

  REAL(DP), ALLOCATABLE    :: q0(:)                  !
  REAL(DP), ALLOCATABLE    :: dq0_drho(:)            ! q-values.
  REAL(DP), ALLOCATABLE    :: dq0_dgradrho(:)        !

  COMPLEX(DP), ALLOCATABLE :: thetas(:,:)            ! Thetas.
  INTEGER                  :: i_proc, theta_i, l, m

  REAL(DP)                 :: sigma_grad(3,3)
  REAL(DP)                 :: sigma_ker(3,3)




  ! --------------------------------------------------------------------
  ! Tests
 
  IF ( inlc > 5 ) CALL errore( 'xc_vdW_DF', 'inlc not implemented', 1 )

#if defined (__SPIN_BALANCED)
  IF ( nspin==2 ) THEN
     WRITE(stdout,'(/,/ "     Performing spin-balanced Ecnl stress calculation!")')
  ELSE IF ( nspin > 2 ) THEN
     CALL errore ('vdW_DF_stress', 'noncollinear vdW stress not implemented', 1)
  END IF
#else
  IF ( nspin>=2 ) THEN
     CALL errore ('vdW_DF_stress', 'vdW stress not implemented for nspin > 1', 1)
  END IF
#endif

  sigma(:,:)      = 0.0_DP
  sigma_grad(:,:) = 0.0_DP
  sigma_ker(:,:)  = 0.0_DP


  ! --------------------------------------------------------------------
  ! Allocations

  ALLOCATE( total_rho(dfftp%nnr), grad_rho(3,dfftp%nnr), thetas(dfftp%nnr, Nqs), &
            q0(dfftp%nnr), dq0_drho(dfftp%nnr), dq0_dgradrho(dfftp%nnr) )


  ! --------------------------------------------------------------------
  ! Charge

  total_rho = rho_valence(:) + rho_core(:)


  ! --------------------------------------------------------------------
  ! Here we calculate the gradient in reciprocal space using FFT.

  CALL fft_gradient_r2r (dfftp, total_rho,  g, grad_rho)


  ! --------------------------------------------------------------------
  ! Get q0.

  CALL get_q0_on_grid (total_rho, grad_rho, q0, dq0_drho, dq0_dgradrho, thetas)


  ! --------------------------------------------------------------------
  ! Stress

  CALL vdW_DF_stress_gradient (total_rho, grad_rho, q0, dq0_drho, dq0_dgradrho, thetas, sigma_grad)
  CALL vdW_DF_stress_kernel   (total_rho, q0, thetas, sigma_ker)

  sigma = - (sigma_grad + sigma_ker)

  DO l = 1, 3
     DO m = 1, l - 1
        sigma (m, l) = sigma (l, m)
     END DO
  END DO

  DEALLOCATE( total_rho, grad_rho, thetas, q0, dq0_drho, dq0_dgradrho )

  END SUBROUTINE vdW_DF_stress








  ! ####################################################################
  !                     |                          |
  !                     |  VDW_DF_STRESS_GRADIENT  |
  !                     |__________________________|

  SUBROUTINE vdW_DF_stress_gradient (total_rho, grad_rho, q0, &
             dq0_drho, dq0_dgradrho, thetas, sigma)

  USE gvect,                 ONLY : ngm, g, gg, igtongl, gl, ngl, gstart
  USE cell_base,             ONLY : omega, tpiba, alat, at, tpiba2

  IMPLICIT NONE

  REAL(DP), INTENT(IN)     :: total_rho(:)           !
  REAL(DP), INTENT(IN)     :: grad_rho(:, :)         !
  REAL(DP), INTENT(INOUT)  :: sigma(:,:)             !
  REAL(DP), INTENT(IN)     :: q0(:)                  ! Input variables.
  REAL(DP), INTENT(IN)     :: dq0_drho(:)            !
  REAL(DP), INTENT(IN)     :: dq0_dgradrho(:)        !
  COMPLEX(DP), INTENT(IN)  :: thetas(:,:)            !

  COMPLEX(DP), ALLOCATABLE :: u_vdW(:,:)             !

  REAL(DP), ALLOCATABLE    :: d2y_dx2(:,:)           !
  REAL(DP) :: y(Nqs), dP_dq0, P, a, b, c, d, e, f    ! Interpolation.
  REAL(DP) :: dq                                     !

  INTEGER  :: q_low, q_hi, q, q1_i, q2_i , g_i       ! Loop and q-points.

  INTEGER  :: l, m
  REAL(DP) :: prefactor                              ! Final summation of sigma
  REAL(DP) :: grad2                                  ! magnitude of density gradient.

  INTEGER  :: i_proc, theta_i, i_grid, q_i, &        ! Iterators.
              ix, iy, iz                             !

  CHARACTER(LEN=1) :: intvar




  ALLOCATE( d2y_dx2(Nqs, Nqs) )
  ALLOCATE( u_vdW(dfftp%nnr, Nqs) )

  sigma(:,:) = 0.0_DP
  prefactor  = 0.0_DP


  ! --------------------------------------------------------------------
  ! Get u in k-space.

  CALL thetas_to_uk(thetas, u_vdW)


  ! --------------------------------------------------------------------
  ! Get u in real space.

  DO theta_i = 1, Nqs
     CALL invfft('Rho', u_vdW(:,theta_i), dfftp)
  END DO


  ! --------------------------------------------------------------------
  ! Get the second derivatives for interpolating the P_i.

  CALL initialize_spline_interpolation(q_mesh, d2y_dx2(:,:))


  ! --------------------------------------------------------------------
  ! Do the real space integration to obtain the stress component.

  DO i_grid = 1, dfftp%nnr

     IF ( total_rho(i_grid) < epsr ) CYCLE

     q_low = 1
     q_hi  = Nqs
     grad2 = sqrt( grad_rho(1,i_grid)**2 + grad_rho(2,i_grid)**2 + grad_rho(3,i_grid)**2 )

     IF ( grad2 == 0.0_dp ) CYCLE

     ! -----------------------------------------------------------------
     ! Figure out which bin our value of q0 is in the q_mesh.
     DO WHILE ( (q_hi - q_low) > 1)

        q = INT((q_hi + q_low)/2)

        IF (q_mesh(q) > q0(i_grid)) THEN
            q_hi = q
        ELSE
            q_low = q
        END IF

     END DO

     IF (q_hi == q_low) call errore('stress_vdW_gradient','qhi == qlow', 1)

     dq = q_mesh(q_hi) - q_mesh(q_low)
     a  = (q_mesh(q_hi) - q0(i_grid))/dq
     b  = (q0(i_grid) - q_mesh(q_low))/dq
     c  = (a**3 - a)*dq**2/6.0D0
     d  = (b**3 - b)*dq**2/6.0D0
     e  = (3.0D0*a**2 - 1.0D0)*dq/6.0D0
     f  = (3.0D0*b**2 - 1.0D0)*dq/6.0D0

     DO q_i = 1, Nqs

        y(:)   = 0.0D0
        y(q_i) = 1.0D0

        dP_dq0 = (y(q_hi) - y(q_low))/dq - e*d2y_dx2(q_i,q_low) + f*d2y_dx2(q_i,q_hi)

        prefactor = u_vdW(i_grid,q_i) * dP_dq0 * dq0_dgradrho(i_grid) / grad2

        DO l = 1, 3
        DO m = 1, l

            sigma (l, m) = sigma (l, m) -  e2 * prefactor * &
                           (grad_rho(l,i_grid) * grad_rho(m,i_grid))
        END DO
        END DO

     END DO

  END DO

  CALL mp_sum(  sigma, intra_bgrp_comm )

  CALL dscal (9, 1.0D0 / (dfftp%nr1 * dfftp%nr2 * dfftp%nr3), sigma, 1)

  DEALLOCATE( d2y_dx2, u_vdW )

  END SUBROUTINE vdW_DF_stress_gradient








  ! ####################################################################
  !                      |                        |
  !                      |  VDW_DF_STRESS_KERNEL  |
  !                      |________________________|

  SUBROUTINE vdW_DF_stress_kernel (total_rho, q0, thetas, sigma)

  USE gvect,                 ONLY : ngm, g, gg, igtongl, gl, ngl, gstart
  USE cell_base,             ONLY : omega, tpiba, tpiba2

  IMPLICIT NONE

  REAL(DP), INTENT(IN)    :: q0(:)
  REAL(DP), INTENT(IN)    :: total_rho(:)
  REAL(DP), INTENT(INOUT) :: sigma(3,3)
  COMPLEX(DP), INTENT(IN) :: thetas(:,:)

  REAL(DP), ALLOCATABLE   :: dkernel_of_dk(:,:)

  INTEGER                 :: l, m, q1_i, q2_i , g_i
  INTEGER                 :: last_g, theta_i
  REAL(DP)                :: g2, ngmod2, g_kernel, G_multiplier




  ALLOCATE( dkernel_of_dk(Nqs, Nqs) )

  sigma(:,:) = 0.0_DP

  ! --------------------------------------------------------------------
  ! Integration in g-space.

  last_g = -1

  G_multiplier = 1.0D0

  IF ( gamma_only ) G_multiplier = 2.0D0

  DO g_i = gstart, ngm

     g2 = gg (g_i) * tpiba2
     g_kernel = SQRT(g2)

     IF ( igtongl(g_i) .NE. last_g) THEN

        CALL interpolate_Dkernel_Dk(g_kernel, dkernel_of_dk)  ! Gets the derivatives.
        last_g = igtongl(g_i)

     END IF

     DO q2_i = 1, Nqs
     DO q1_i = 1, Nqs
        DO l = 1, 3
        DO m = 1, l

           sigma (l, m) = sigma (l, m) - G_multiplier * 0.5 * e2 * thetas(dfftp%nl(g_i),q1_i) * &
                          dkernel_of_dk(q1_i,q2_i)*conjg(thetas(dfftp%nl(g_i),q2_i))* &
                          (g (l, g_i) * g (m, g_i) * tpiba2) / g_kernel
        END DO
        END DO
     END DO
     END DO

     IF ( g_i < gstart ) sigma(:,:) = sigma(:,:) / G_multiplier

  END DO

  CALL mp_sum( sigma, intra_bgrp_comm )

  DEALLOCATE( dkernel_of_dk )

  END SUBROUTINE vdW_DF_stress_kernel








  ! ####################################################################
  !                        |                        |
  !                        | INTERPOLATE_DKERNEL_DK |
  !                        |________________________|

  SUBROUTINE interpolate_Dkernel_Dk (k, dkernel_of_dk)

  IMPLICIT NONE

  REAL(DP), INTENT(IN)    :: k                      ! Input value, the magnitude of the g-vector
                                                    ! for the current point.

  REAL(DP), INTENT(INOUT) :: dkernel_of_dk(Nqs,Nqs) ! An output array (allocated outside this
                                                    ! routine) that holds the interpolated value of
                                                    ! the kernel for each pair of q points (i.e. the
                                                    ! phi_alpha_beta of the Soler method.

  INTEGER :: q1_i, q2_i, k_i                        ! Indexing variables.

  REAL(DP) :: A, B, dAdk, dBdk, dCdk, dDdk          ! Intermediate values for the interpolation.




  IF ( k >= Nr_points*dk ) THEN

     WRITE(*,'(A,F10.5,A,F10.5)') "k =  ", k, "     k_max =  ", Nr_points*dk
     CALL errore('interpolate kernel', 'k value requested is out of range',1)

  END IF

  dkernel_of_dk = 0.0D0

  k_i  = INT(k/dk)

  A    = (dk*(k_i+1.0D0) - k)/dk
  B    = (k - dk*k_i)/dk

  dAdk = -1.0D0/dk
  dBdk = 1.0D0/dk
  dCdk = -((3*A**2 -1.0D0)/6.0D0)*dk
  dDdk = ((3*B**2 -1.0D0)/6.0D0)*dk

  DO q1_i = 1, Nqs
     DO q2_i = 1, q1_i

        dkernel_of_dk(q1_i, q2_i) = dAdk*kernel(k_i, q1_i, q2_i) + dBdk*kernel(k_i+1, q1_i, q2_i) &
                            + dCdk*d2phi_dk2(k_i, q1_i, q2_i) + dDdk*d2phi_dk2(k_i+1, q1_i, q2_i)

        dkernel_of_dk(q2_i, q1_i) = dkernel_of_dk(q1_i, q2_i)

     END DO
  END DO

  END SUBROUTINE interpolate_Dkernel_Dk








  ! ####################################################################
  !                          |              |
  !                          | thetas_to_uk |
  !                          |______________|


  SUBROUTINE thetas_to_uk (thetas, u_vdW)

  USE gvect,           ONLY : gg, ngm, igtongl, gl, ngl, gstart
  USE cell_base,       ONLY : tpiba, omega

  IMPLICIT NONE

  COMPLEX(DP), INTENT(IN)  :: thetas(:,:)       ! On input this variable holds the theta functions
                                                ! (equation 8, SOLER) in the format
                                                ! thetas(grid_point, theta_i).
  COMPLEX(DP), INTENT(OUT) :: u_vdW(:,:)        ! On output this array holds u_alpha(k) =
                                                ! Sum_j[theta_beta(k)phi_alpha_beta(k)].

  REAL(DP), ALLOCATABLE    :: kernel_of_k(:,:)  ! This array will hold the interpolated kernel
                                                ! values for each pair of q values in the q_mesh.

  REAL(DP) :: g
  INTEGER  :: last_g, g_i, q1_i, q2_i, i_grid   ! Index variables.

  COMPLEX(DP) :: theta(Nqs)                     ! Temporary storage vector used since we are
                                                ! overwriting the thetas array here.




  ALLOCATE( kernel_of_k(Nqs, Nqs) )

  u_vdW(:,:) = CMPLX(0.0_DP, 0.0_DP, kind=dp)

  last_g = -1

  DO g_i = 1, ngm

     IF ( igtongl(g_i) .ne. last_g) THEN

        g = SQRT(gl(igtongl(g_i))) * tpiba
        CALL interpolate_kernel(g, kernel_of_k)
        last_g = igtongl(g_i)

     END IF

     theta = thetas(dfftp%nl(g_i),:)

     DO q2_i = 1, Nqs
        DO q1_i = 1, Nqs
           u_vdW(dfftp%nl(g_i),q2_i) = u_vdW(dfftp%nl(g_i),q2_i) + kernel_of_k(q2_i,q1_i)*theta(q1_i)
        END DO
     END DO

  END Do

  IF ( gamma_only ) u_vdW(dfftp%nlm(:),:) = CONJG(u_vdW(dfftp%nl(:),:))

  DEALLOCATE( kernel_of_k )

  END SUBROUTINE thetas_to_uk








  ! ####################################################################
  !                           |                 |
  !                           | GENERATE_KERNEL |
  !                           |_________________|
  !
  ! The original definition of the kernel function is given in DION
  ! equations 14-16. The Soler method makes the kernel function a
  ! function of only 1 variable (r) by first putting it in the form
  ! phi(q1*r, q2*r). Then, the q-dependence is removed by expanding the
  ! function in a special way (see SOLER equation 3). This yields a
  ! separate function for each pair of q points that is a function of r
  ! alone. There are (Nqs^2+Nqs)/2 unique functions, where Nqs is the
  ! number of q points used. In the Soler method, the kernel is first
  ! made in the form phi(d1, d2) but this is not done here. It was found
  ! that, with q's chosen judiciously ahead of time, the kernel and the
  ! second derivatives required for interpolation could be tabulated
  ! ahead of time for faster use of the vdW-DF functional. Through
  ! testing we found no need to soften the kernel and correct for this
  ! later (see SOLER eqations 6-7).
  !
  ! The algorithm employed here is "embarrassingly parallel," meaning
  ! that it parallelizes very well up to (Nqs^2+Nqs)/2 processors,
  ! where, again, Nqs is the number of q points chosen. However,
  ! parallelization on this scale is unnecessary. In testing the code
  ! runs in under a minute on 16 Intel Xeon processors.
  !
  ! IMPORTANT NOTICE: Results are very sensitive to compilation details.
  ! In particular, the usage of FMA (Fused Multiply-and-Add)
  ! instructions used by modern CPUs such as AMD Interlagos (Bulldozer)
  ! and Intel Ivy Bridge may affect quite heavily some components of the
  ! kernel (communication by Ake Sandberg, Umea University). In practice
  ! this should not be a problem, since most affected elements are the
  ! less relevant ones.
  !
  ! Some of the algorithms here are somewhat modified versions of those
  ! found in:
  !
  !    Numerical Recipes in C; William H. Press, Brian P. Flannery, Saul
  !    A. Teukolsky, and William T. Vetterling. Cambridge University
  !    Press (1988).
  !
  ! hereafter referred to as NUMERICAL_RECIPES. The routines were
  ! translated to Fortran, of course and variable names are generally
  ! different.
  !
  ! For the calculation of the kernel we have benefited from access to
  ! earlier vdW-DF implementation into PWscf and ABINIT, written by Timo
  ! Thonhauser, Valentino Cooper, and David Langreth. These codes, in
  ! turn, benefited from earlier codes written by Maxime Dion and Henrik
  ! Rydberg.

  SUBROUTINE generate_kernel

  IMPLICIT NONE

  INTEGER  :: a_i, b_i, q1_i, q2_i, r_i
  ! Indexing variables.

  REAL(DP) :: weights( Nintegration_points )
  ! Array to hold dx values for the Gaussian-Legendre integration of the kernel.

  REAL(DP) :: sin_a( Nintegration_points ), cos_a( Nintegration_points )
  ! Sine and cosine values of the aforementioned points a.

  REAL(DP) :: d1, d2, d, integral
  ! Intermediate values.

  ! --------------------------------------------------------------------
  ! The following variables control the parallel environment.

  INTEGER :: my_start_q, my_end_q, Ntotal
  ! Starting and ending q value for each  processor, also the total
  ! number of calculations to do, i.e. (Nqs^2 + Nqs)/2.

  REAL(DP), ALLOCATABLE :: phi(:,:), phi_deriv(:,:)
  ! Arrays to store the kernel functions and their second derivatives.
  ! They are stored as phi(radial_point, idx).

  INTEGER, ALLOCATABLE  :: indices(:,:), proc_indices(:,:)
  ! Indices holds the values of q1 and q2 as partitioned out to the
  ! processors. It is an Ntotal x 2 array stored as indices(index of
  ! point number, q1:q2). Proc_indices holds the section of the indices
  ! array that is assigned to each processor. This is a Nproc x 2
  ! array, stored as proc_indices(processor_number,
  ! starting_index:ending_index)

  INTEGER :: Nper, Nextra, start_q, end_q
  ! Baseline number of jobs per processor, number of processors that
  ! get an extra job in case the number of jobs doesn't split evenly
  ! over the number of processors, starting index into the indices
  ! array, ending index into the indices array.

  INTEGER :: nproc, mpime
  ! Number or procs, rank of current processor.

  INTEGER :: proc_i, my_Nqs




  ! --------------------------------------------------------------------
  ! Start the timer.

  CALL start_clock ( 'vdW_kernel' )


  ! --------------------------------------------------------------------
  ! The total number of phi_alpha_beta functions that have to be
  ! calculated.

  Ntotal = (Nqs**2 + Nqs)/2
  ALLOCATE ( indices(Ntotal, 2) )


  ! --------------------------------------------------------------------
  ! This part fills in the indices array. It just loops through the q1
  ! and q2 values and stores them. Sections of this array will be
  ! assigned to each of the processors later.

  idx = 1

  DO q1_i = 1, Nqs
     DO q2_i = 1, q1_i
        indices(idx, 1) = q1_i
        indices(idx, 2) = q2_i
        idx = idx + 1
     END DO
  END DO


  ! --------------------------------------------------------------------
  ! Figure out the baseline number of functions to be calculated by each
  ! processor and how many processors get one extra job.

  nproc  = mp_size( intra_image_comm )
  mpime  = mp_rank( intra_image_comm )
  Nper   = Ntotal/nproc
  Nextra = MOD(Ntotal, nproc)

  ALLOCATE( proc_indices(nproc, 2) )

  start_q = 0
  end_q   = 0


  ! --------------------------------------------------------------------
  ! Loop over all the processors and figure out which section of the
  ! indices array each processor should do. All processors figure this
  ! out for every processor so there is no need to communicate results.

  DO proc_i = 1, nproc

     start_q = end_q + 1
     end_q   = start_q + (Nper - 1)
     IF (proc_i <= Nextra) end_q = end_q + 1

     ! This is to prevent trouble if number of processors exceeds Ntotal.
     IF ( proc_i > Ntotal ) THEN
        start_q    = Ntotal
        end_q      = Ntotal
     END IF

     IF ( proc_i == (mpime+1) ) THEN
        my_start_q = start_q
        my_end_q   = end_q
     END IF

     proc_indices(proc_i, 1) = start_q
     proc_indices(proc_i, 2) = end_q

  END DO


  ! --------------------------------------------------------------------
  ! Store how many jobs are assigned to me.

  my_Nqs    = my_end_q - my_start_q + 1
  ALLOCATE( phi( 0:Nr_points, my_Nqs ), phi_deriv( 0:Nr_points, my_Nqs ) )

  phi       = 0.0D0
  phi_deriv = 0.0D0
  kernel    = 0.0D0
  d2phi_dk2 = 0.0D0


  ! --------------------------------------------------------------------
  ! Find the integration points we are going to use in the
  ! Gaussian-Legendre integration.

  CALL prep_gaussian_quadrature( weights )


  ! --------------------------------------------------------------------
  ! Get a, a^2, sin(a), cos(a) and the weights for the Gaussian-Legendre
  ! integration.

  DO a_i=1, Nintegration_points
     a_points (a_i) = TAN( a_points(a_i) )
     a_points2(a_i) = a_points(a_i)**2
     weights(a_i)   = weights(a_i)*( 1 + a_points2(a_i) )
     cos_a(a_i)     = COS( a_points(a_i) )
     sin_a(a_i)     = SIN( a_points(a_i) )
  END DO


  ! --------------------------------------------------------------------
  ! Calculate the value of the W function defined in DION equation 16
  ! for each value of a and b.

  DO a_i = 1, Nintegration_points
  DO b_i = 1, Nintegration_points
     W_ab(a_i, b_i) = 2.0D0 * weights(a_i)*weights(b_i) * (           &
        (3.0D0-a_points2(a_i))*a_points(b_i) *sin_a(a_i)*cos_a(b_i) + &
        (3.0D0-a_points2(b_i))*a_points(a_i) *cos_a(a_i)*sin_a(b_i) + &
        (a_points2(a_i)+a_points2(b_i)-3.0D0)*sin_a(a_i)*sin_a(b_i) - &
        3.0D0*a_points(a_i)*a_points(b_i)*cos_a(a_i)*cos_a(b_i) )   / &
        (a_points(a_i)*a_points(b_i) )
  END DO
  END DO


  ! --------------------------------------------------------------------
  ! vdW-DF analysis tool as described in PRB 97, 085115 (2018).

  IF      ( vdW_DF_analysis == 1 ) THEN

     DO a_i = 1, Nintegration_points
     DO b_i = 1, Nintegration_points
        W_ab(a_i, b_i) = weights(a_i)*weights(b_i) *                  &
           a_points(a_i)*a_points(b_i)*sin_a(a_i)*sin_a(b_i)
     END DO
     END DO

  ELSE IF ( vdW_DF_analysis == 2 ) THEN

     DO a_i = 1, Nintegration_points
     DO b_i = 1, Nintegration_points
        W_ab(a_i, b_i) = W_ab(a_i, b_i) - weights(a_i)*weights(b_i) *  &
           a_points(a_i)*a_points(b_i)*sin_a(a_i)*sin_a(b_i)
     END DO
     END DO

  END IF


  ! --------------------------------------------------------------------
  ! Now, we loop over all the pairs (q1,q2) that are assigned to us and
  ! perform our calculations.

  DO idx = 1, my_Nqs

     ! -----------------------------------------------------------------
     ! First, get the value of phi(q1*r, q2*r) for each r and the
     ! particular values of q1 and q2 we are using.

     DO r_i = 1, Nr_points
        d1  = q_mesh( indices(idx+my_start_q-1, 1) ) * dr * r_i
        d2  = q_mesh( indices(idx+my_start_q-1, 2) ) * dr * r_i
        phi(r_i, idx) = phi_value(d1, d2)
     END DO


     ! -----------------------------------------------------------------
     ! Now, perform a radial FFT to turn our phi_alpha_beta(r) into
     ! phi_alpha_beta(k) needed for SOLER equation 8.

     CALL radial_fft( phi(:,idx) )


     ! -----------------------------------------------------------------
     ! Determine the spline interpolation coefficients for the Fourier
     ! transformed kernel function.

     CALL set_up_splines( phi(:, idx), phi_deriv(:, idx) )

  END DO


  ! --------------------------------------------------------------------
  ! Finally, we collect the results after letting everybody catch up.

  CALL mp_barrier( intra_image_comm )

  DO proc_i = 0, nproc-1

     IF ( proc_i >= Ntotal ) EXIT

     CALL mp_get ( phi      , phi      , mpime, 0, proc_i, 0, intra_image_comm )
     CALL mp_get ( phi_deriv, phi_deriv, mpime, 0, proc_i, 0, intra_image_comm )

     IF ( mpime == 0 ) THEN

        DO idx = proc_indices(proc_i+1,1), proc_indices(proc_i+1,2)
           q1_i = indices(idx, 1)
           q2_i = indices(idx, 2)
           kernel    (:, q1_i, q2_i) = phi       (:, idx - proc_indices(proc_i+1,1) + 1)
           d2phi_dk2 (:, q1_i, q2_i) = phi_deriv (:, idx - proc_indices(proc_i+1,1) + 1)
           kernel    (:, q2_i, q1_i) = kernel    (:, q1_i, q2_i)
           d2phi_dk2 (:, q2_i, q1_i) = d2phi_dk2 (:, q1_i, q2_i)
        END DO

     END IF

  END DO

  CALL mp_bcast ( kernel   , 0, intra_image_comm )
  CALL mp_bcast ( d2phi_dk2, 0, intra_image_comm )


  ! --------------------------------------------------------------------
  ! Keep the lines below for testing and combatibility with the old
  ! kernel file reading/writing method.
  !
  ! Writing the calculated kernel.
  !
  ! IF ( ionode ) THEN
  !    WRITE(stdout,'(/ / A)') "     vdW-DF kernel table calculated and written to file."
  !    OPEN(UNIT=21, FILE='kernel_table', STATUS='replace', FORM='formatted', ACTION='write')
  !    WRITE(21, '(2i5,f13.8)') Nqs, Nr_points
  !    WRITE(21, '(1p4e23.14)') r_max
  !    WRITE(21, '(1p4e23.14)') q_mesh
  !    DO q1_i = 1, Nqs
  !       DO q2_i = 1, q1_i
  !          WRITE(21, '(1p4e23.14)') kernel(:, q1_i, q2_i)
  !       END DO
  !    END DO
  !    DO q1_i = 1, Nqs
  !       DO q2_i = 1, q1_i
  !          WRITE(21, '(1p4e23.14)') d2phi_dk2(:, q1_i, q2_i)
  !       END DO
  !    END DO
  !    CLOSE (21)
  ! END IF
  !
  !
  ! Reading the kernel from an old kernel file.
  !
  ! IF (ionode) WRITE(stdout,'(/ / A)') "     vdW-DF kernel read from file."
  ! OPEN(UNIT=21, FILE='vdW_kernel_table', STATUS='old', FORM='formatted', ACTION='read')
  ! read(21, '(/ / / / / /)')
  ! DO q1_i = 1, Nqs
  !    DO q2_i = 1, q1_i
  !       READ(21, '(1p4e23.14)') kernel(:, q1_i, q2_i)
  !       kernel(:, q2_i, q1_i) = kernel(:, q1_i, q2_i)
  !    END DO
  ! END DO
  ! DO q1_i = 1, Nqs
  !    DO q2_i = 1, q1_i
  !       READ(21, '(1p4e23.14)')    d2phi_dk2(:, q1_i, q2_i)
  !       d2phi_dk2(:, q2_i, q1_i) = d2phi_dk2(:, q1_i, q2_i)
  !    END DO
  ! END DO
  ! CLOSE (21)


  DEALLOCATE( indices, proc_indices, phi, phi_deriv )


  ! --------------------------------------------------------------------
  ! Stop the timer.

  CALL stop_clock ( 'vdW_kernel' )

  END SUBROUTINE generate_kernel








  ! ####################################################################
  !                    |                            |
  !                    |  PREP_GAUSSIAN_QUADRATURE  |
  !                    |____________________________|
  !
  ! Routine to calculate the points and weights for the
  ! Gaussian-Legendre integration. This routine is modeled after the
  ! routine GAULEG from NUMERICAL_RECIPES.

  SUBROUTINE prep_gaussian_quadrature( weights )

  REAL(DP), INTENT(INOUT) :: weights(:)
  ! The points and weights for the Gaussian-Legendre integration.

  INTEGER  :: Npoints
  ! The number of points we actually have to calculate. The rest will
  ! be obtained from symmetry.

  REAL(DP) :: poly_1, poly_2, poly_3
  ! Temporary storage for Legendre polynomials.

  INTEGER  :: i_point, i_poly
  ! Indexing variables.

  REAL(DP) :: root, dp_dx, last_root
  ! The value of the root of a given Legendre polynomial, the derivative
  ! of the polynomial at that root and the value of the root in the last
  ! iteration (to check for convergence of Newton's method).

  real(dp) :: midpoint, length
  ! The middle of the x-range and the length to that point.




  Npoints  = (Nintegration_points + 1)/2
  midpoint = 0.5D0 * ( ATAN(a_min) + ATAN(a_max) )
  length   = 0.5D0 * ( ATAN(a_max) - ATAN(a_min) )

  DO i_point = 1, Npoints
     ! -----------------------------------------------------------------
     ! Make an initial guess for the root.

     root = COS(DBLE(pi*(i_point - 0.25D0)/(Nintegration_points + 0.5D0)))

     DO
        ! --------------------------------------------------------------
        ! Use the recurrence relations to find the desired polynomial,
        ! evaluated at the approximate root. See NUMERICAL_RECIPES.

        poly_1 = 1.0D0
        poly_2 = 0.0D0

        DO i_poly = 1, Nintegration_points

           poly_3 = poly_2
           poly_2 = poly_1
           poly_1 = ((2.0D0 * i_poly - 1.0D0)*root*poly_2 - (i_poly-1.0D0)*poly_3)/i_poly

        END DO


        ! --------------------------------------------------------------
        ! Use the recurrence relations to find the desired polynomial.
        ! Find the derivative of the polynomial and use it in Newton's
        ! method to refine our guess for the root.

        dp_dx = Nintegration_points * (root*poly_1 - poly_2)/(root**2 - 1.0D0)

        last_root = root
        root      = last_root - poly_1/dp_dx


        ! --------------------------------------------------------------
        ! Check for convergence.

        IF (abs(root - last_root) <= 1.0D-14) EXIT

     END DO


     ! -----------------------------------------------------------------
     ! Fill in the array of evaluation points.

     a_points(i_point) = midpoint - length*root
     a_points(Nintegration_points + 1 - i_point) = midpoint + length*root


     ! -----------------------------------------------------------------
     ! Fill in the array of weights.

     weights(i_point) = 2.0D0 * length/((1.0D0 - root**2)*dp_dx**2)
     weights(Nintegration_points + 1 - i_point) = weights(i_point)

  END DO

  END SUBROUTINE prep_gaussian_quadrature








  ! ####################################################################
  !                            |             |
  !                            |  PHI_VALUE  |
  !                            |_____________|
  !
  ! This function returns the value of the kernel calculated via DION
  ! equation 14.

  REAL(DP) FUNCTION phi_value(d1, d2)

  REAL(DP), INTENT(IN) :: d1, d2
  ! The point at which to evaluate the kernel. d1 = q1*r and d2 = q2*r.

  REAL(DP) :: w, x, y, z, T
  ! Intermediate values.

  REAL(DP) :: nu(Nintegration_points), nu1(Nintegration_points)
  ! Defined in the discussio below equation 16 of DION.

  INTEGER  :: a_i, b_i
  ! Indexing variables.




  ! --------------------------------------------------------------------
  ! Loop over all integration points and calculate the value of the nu
  ! functions defined in the discussion below equation 16 in DION.

  DO a_i = 1, Nintegration_points
     nu(a_i)  = a_points2(a_i)/( 2.0D0 * h_function( a_points(a_i)/d1 ))
     nu1(a_i) = a_points2(a_i)/( 2.0D0 * h_function( a_points(a_i)/d2 ))
  END DO


  ! --------------------------------------------------------------------
  ! Carry out the integration of DION equation 13.

  phi_value = 0.0D0

  DO a_i = 1, Nintegration_points
     w = nu(a_i)
     y = nu1(a_i)
     DO b_i = 1, Nintegration_points
        x = nu(b_i)
        z = nu1(b_i)
        T = (1.0D0/(w+x) + 1.0D0/(y+z))*(1.0D0/((w+y)*(x+z)) + 1.0D0/((w+z)*(y+x)))
        phi_value = phi_value + T * W_ab(a_i, b_i)
     END DO
  END DO

  phi_value = 1.0D0/pi**2*phi_value

  END FUNCTION phi_value








  ! ####################################################################
  !                            |              |
  !                            |  RADIAL_FFT  |
  !                            |______________|
  !
  ! This subroutine performs a radial Fourier transform on the
  ! real-space kernel functions. Basically, this is just
  ! int(4*pi*r^2*phi*sin(k*r)/(k*r))dr integrated from 0 to r_max. That
  ! is, it is the kernel function phi integrated with the 0^th spherical
  ! Bessel function radially, with a 4*pi assumed from angular
  ! integration since we have spherical symmetry. The spherical symmetry
  ! comes in because the kernel function depends only on the magnitude
  ! of the vector between two points. The integration is done using the
  ! trapezoid rule.

  SUBROUTINE radial_fft(phi)

  REAL(DP), INTENT(INOUT) :: phi(0:Nr_points)
  ! On input holds the real-space function phi_q1_q2(r).
  ! On output hold the reciprocal-space function phi_q1_q2(k).

  REAL(DP) :: phi_k(0:Nr_points)
  ! Temporary storage for phi_q1_q2(k).

  INTEGER  :: k_i, r_i
  ! Indexing variables.

  REAL(DP) :: r, k
  ! The real and reciprocal space points.




  phi_k = 0.0D0

  ! --------------------------------------------------------------------
  ! Handle the k=0 point separately.

  DO r_i = 1, Nr_points
     r        = r_i * dr
     phi_k(0) = phi_k(0) + phi(r_i)*r**2
  END DO


  ! --------------------------------------------------------------------
  ! Subtract half of the last value off because of the trapezoid rule.

  phi_k(0) = phi_k(0) - 0.5D0 * (Nr_points*dr)**2 * phi(Nr_points)


  ! --------------------------------------------------------------------
  ! Integration for the rest of the k-points.

  DO k_i = 1, Nr_points
     k = k_i * dk
     DO r_i = 1, Nr_points
        r          = r_i * dr
        phi_k(k_i) = phi_k(k_i) + phi(r_i) * r * SIN(k*r) / k
     END DO
     phi_k(k_i) = phi_k(k_i) - 0.5D0 * phi(Nr_points) * r * SIN(k*r) / k
  END DO


  ! --------------------------------------------------------------------
  ! Add in the 4*pi and the dr factor for the integration.

  phi = 4.0D0 * pi * phi_k * dr

  END SUBROUTINE radial_fft








  ! ####################################################################
  !                          |                  |
  !                          |  SET UP SPLINES  |
  !                          |__________________|
  !
  ! This subroutine accepts a function (phi) and finds at each point the
  ! second derivative (D2) for use with spline interpolation. This
  ! function assumes we are using the expansion described in SOLER
  ! equation 3. That is, the derivatives are those needed to interpolate
  ! Kronecker delta functions at each of the q values. Other than some
  ! special modification to speed up the algorithm in our particular
  ! case, this algorithm is taken directly from NUMERICAL_RECIPES.

  SUBROUTINE set_up_splines(phi, D2)

  REAL(DP), INTENT(IN)    :: phi(0:Nr_points)
  ! The k-space kernel function for a particular q1 and q2.

  REAL(DP), INTENT(INOUT) :: D2(0:Nr_points)
  ! The second derivatives to be used in the interpolation expansion
  ! (SOLER equation 3).

  REAL(DP), ALLOCATABLE   :: temp_array(:)         ! Temporary storage.
  REAL(DP)                :: temp_1, temp_2

  INTEGER  :: r_i
  ! Indexing variable.




  ALLOCATE( temp_array(0:Nr_points) )

  D2         = 0
  temp_array = 0

  DO r_i = 1, Nr_points - 1
     temp_1  = DBLE(r_i - (r_i - 1))/DBLE( (r_i + 1) - (r_i - 1) )
     temp_2  = temp_1 * D2(r_i-1) + 2.0D0
     D2(r_i) = (temp_1 - 1.0D0)/temp_2
     temp_array(r_i) = ( phi(r_i+1) - phi(r_i))/DBLE( dk*((r_i+1) - r_i) ) - &
          ( phi(r_i) - phi(r_i-1))/DBLE( dk*(r_i - (r_i-1)) )
     temp_array(r_i) = (6.0D0*temp_array(r_i)/DBLE( dk*((r_i+1) - (r_i-1)) )-&
          temp_1*temp_array(r_i-1))/temp_2
  END DO

  D2(Nr_points) = 0.0D0
  DO  r_i = Nr_points-1, 0, -1
     D2(r_i) = D2(r_i)*D2(r_i+1) + temp_array(r_i)
  END DO

  DEALLOCATE( temp_array )

  END SUBROUTINE set_up_splines








  ! ####################################################################
  !                          |            |
  !                          |  VDW_INFO  |
  !                          |____________|

  SUBROUTINE vdW_info

  IMPLICIT NONE




  WRITE(stdout,'(/)')
  WRITE(stdout,'(5x,"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")')
  WRITE(stdout,'(5x,"%                                                                      %")')
  WRITE(stdout,'(5x,"% You are using vdW-DF, which was implemented by the Thonhauser group. %")')
  WRITE(stdout,'(5x,"% Please cite the following two papers that made this development      %")')
  WRITE(stdout,'(5x,"% possible and the two reviews that describe the various versions:     %")')
  WRITE(stdout,'(5x,"%                                                                      %")')
  WRITE(stdout,'(5x,"%   T. Thonhauser et al., PRL 115, 136402 (2015).                      %")')
  WRITE(stdout,'(5x,"%   T. Thonhauser et al., PRB 76, 125112 (2007).                       %")')
  WRITE(stdout,'(5x,"%   K. Berland et al., Rep. Prog. Phys. 78, 066501 (2015).             %")')
  WRITE(stdout,'(5x,"%   D.C. Langreth et al., J. Phys.: Condens. Matter 21, 084203 (2009). %")')
  WRITE(stdout,'(5x,"%                                                                      %")')
  WRITE(stdout,'(5x,"%                                                                      %")')
  WRITE(stdout,'(5x,"% If you are calculating the stress with vdW-DF, please also cite:     %")')
  WRITE(stdout,'(5x,"%                                                                      %")')
  WRITE(stdout,'(5x,"%   R. Sabatini et al., J. Phys.: Condens. Matter 24, 424209 (2012).   %")')
  WRITE(stdout,'(5x,"%                                                                      %")')
  WRITE(stdout,'(5x,"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")')
  WRITE(stdout,'(/)')

  IF ( iverbosity > 0 ) THEN
     WRITE(stdout,'(5x,"Carrying out vdW-DF run using the following parameters:")')
     WRITE(stdout,'(5X,A,I3,A,I5,A,F8.3)' ) "Nqs    = ", Nqs, "  Npoints = ", Nr_points, &
                  "  r_max = ", r_max
     WRITE(stdout,'(5X,"q_mesh =",4F12.8)') (q_mesh(idx), idx=1, 4)
     WRITE(stdout,'(13X,4F12.8)') (q_mesh(idx), idx=5, Nqs)
  END IF

  END SUBROUTINE




END MODULE vdW_DF
