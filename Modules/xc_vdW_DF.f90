! Copyright (C) 2001-2009 Quantum ESPRESSO group
! Copyright (C) 2015 Brian Kolb, Timo Thonhauser
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
! ----------------------------------------------------------------------


MODULE vdW_DF


! This module calculates the non-local correlation contribution to the
! energy and potential according to
!
!    M. Dion, H. Rydberg, E. Schroeder, D. C. Langreth, and
!    B. I. Lundqvist, Phys. Rev. Lett. 92, 246401 (2004).
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
! Two review article show many of the vdW-DF applications:
!
!    D. C. Langreth et al., J. Phys.: Condens. Matter 21, 084203 (2009).
!
!    K. Berland et al, Rep. Prog. Phys. 78, 066501 (2015).
!
!
! The method implemented is based on the method of G. Roman-Perez and
! J. M. Soler described in:
!
!    G. Roman-Perez and J. M. Soler, PRL 103, 096102 (2009).
!
! henceforth referred to as SOLER.
!
!
! There are a number of subroutines in this file. All are used only by
! other subroutines here except for the xc_vdW_DF subroutine, which is
! the driver routine for the vdW-DF calculations and is called from
! v_of_rho. This routine handles setting up the parallel run (if any)
! and carries out the calls necessary to calculate the non-local
! correlation contributions to the energy and potential.


USE kinds,             ONLY : dp
USE constants,         ONLY : pi, e2
USE kernel_table,      ONLY : q_mesh, Nr_points, Nqs, r_max, q_cut, q_min, kernel, d2phi_dk2, dk
USE mp,                ONLY : mp_bcast, mp_sum, mp_barrier
USE mp_bands,          ONLY : intra_bgrp_comm
USE io_global,         ONLY : stdout, ionode
USE fft_base,          ONLY : dfftp
USE fft_interfaces,    ONLY : fwfft, invfft
USE control_flags,     ONLY : iverbosity, gamma_only

implicit none

REAL(DP), PARAMETER :: epsr =1.d-12  ! a small number to cut off densities
integer :: vdw_type = 1

private
public  :: xc_vdW_DF, xc_vdW_DF_spin, stress_vdW_DF, interpolate_kernel, &
           vdw_type, numerical_gradient, initialize_spline_interpolation

CONTAINS








  ! ####################################################################
  !                           |             |
  !                           |  functions  |
  !                           |_____________|
  !
  ! Functions to be used in get_q0_on_grid and get_q0_on_grid_spin().

  function Fs(s)

     implicit none
     real(dp) :: s, Fs, Z_ab
     real(dp) :: fa=0.1234D0, fb=17.33D0, fc=0.163D0  ! Reparameterized values
                                                      ! from JCTC 5, 2745 (2009).

     if(vdw_type == 4) then
         Fs = ( 1 + 15.0D0*fa*s**2 + fb*s**4 + fc*s**6 )**(1.0D0/15.0D0)
     else
         ! ------------------------------------------------------------
         ! Original functional choice for Fs, as definded in DION
         if (vdw_type == 1) Z_ab = -0.8491D0
         if (vdw_type == 2) Z_ab = -1.887D0
         Fs = 1.0D0 - Z_ab * s**2 / 9.0D0
     end if

  end function Fs




  function dFs_ds(s)

     implicit none
     real(dp) :: s, dFs_ds, Z_ab
     real(dp) :: fa=0.1234D0, fb=17.33D0, fc=0.163D0  ! Reparameterized values
                                                      ! from JCTC 5, 2745 (2009).

     if(vdw_type == 4) then
         dFs_ds = ( 30.0D0*fa*s + 4.0D0*fb*s**3 + 6.0D0*fc*s**5 ) &
                / ( 15.0D0*( 1.0D0 + 15.0D0*fa*s**2 + fb*s**4 + fc*s**6 )**(14.0D0/15.0D0) )
     else
         ! ------------------------------------------------------------
         ! Original functional choice for Fs, as definded in DION
         if (vdw_type == 1) Z_ab = -0.8491D0
         if (vdw_type == 2) Z_ab = -1.887D0
         dFs_ds =  -2.0D0 * s * Z_ab / 9.0D0
     end if

  end function dFs_ds




  function kF(rho)

     implicit none
     real(dp) :: rho, kF

     kF = ( 3.0D0 * pi**2 * rho )**(1.0D0/3.0D0)

  end function kF




  function dkF_drho(rho)

     implicit none
     real(dp) :: rho, dkF_drho

     dkF_drho = (1.0D0/3.0D0) * kF(rho) / rho

  end function dkF_drho




  function ds_drho(rho, s)

     implicit none
     real(dp) :: rho, s, ds_drho

     ds_drho = -s * ( dkF_drho(rho) / kF(rho) + 1.0D0 / rho )

  end function ds_drho




  function ds_dgradrho(rho)

     implicit none
     real(dp) :: rho, ds_dgradrho

     ds_dgradrho = 1.0D0 / (2.0D0 * kF(rho) * rho)

  end function ds_dgradrho




  function dqx_drho(rho, s)

     implicit none
     real(dp) :: rho, s, dqx_drho

     dqx_drho = dkF_drho(rho) * Fs(s) + kF(rho) * dFs_ds(s) * ds_drho(rho, s)

  end function dqx_drho








  ! ####################################################################
  !                           |             |
  !                           |  XC_VDW_DF  |
  !                           |_____________|

  SUBROUTINE xc_vdW_DF (rho_valence, rho_core, etxc, vtxc, v)

  USE gvect,                 ONLY : ngm, nl, g, nlm
  USE cell_base,             ONLY : omega, tpiba

  implicit none

  ! --------------------------------------------------------------------
  ! Local variables
  !                                               _
  real(dp), intent(IN) :: rho_valence(:,:)       !
  real(dp), intent(IN) :: rho_core(:)            !  PWSCF input variables
  real(dp), intent(inout) :: etxc, vtxc, v(:,:)  !_


  integer :: i_grid, theta_i, i_proc             ! Indexing variables over grid points,
                                                 ! theta functions, and processors.

  real(dp) :: grid_cell_volume                   ! The volume of the unit cell per G-grid point.

  real(dp), allocatable ::  q0(:)                ! The saturated value of q (equations 11 and 12
                                                 ! of DION). This saturation is that of
                                                 ! equation 5 in SOLER.

  real(dp), allocatable :: grad_rho(:,:)         ! The gradient of the charge density. The
                                                 ! format is as follows:
                                                 ! grad_rho(grid_point, cartesian_component).

  real(dp), allocatable :: potential(:)          ! The vdW contribution to the potential

  real(dp), allocatable :: dq0_drho(:)           ! The derivative of the saturated q0
                                                 ! (equation 5 of SOLER) with respect
                                                 ! to the charge density (see
                                                 ! get_q0_on_grid subroutine for details).

  real(dp), allocatable :: dq0_dgradrho(:)       ! The derivative of the saturated q0
                                                 ! (equation 5 of SOLER) with respect
                                                 ! to the gradient of the charge density
                                                 ! (again, see get_q0_on_grid subroutine).

  complex(dp), allocatable :: thetas(:,:)        ! These are the functions of equation 8 of
                                                 ! SOLER. They will be forward Fourier transformed
                                                 ! in place to get theta(k) and worked on in
                                                 ! place to get the u_alpha(r) of equation 11
                                                 ! in SOLER. They are formatted as follows:
                                                 ! thetas(grid_point, theta_i).

  real(dp) :: Ec_nl                              ! The non-local vdW contribution to the energy.

  real(dp), allocatable :: total_rho(:)          ! This is the sum of the valence and core
                                                 ! charge. This just holds the piece assigned
                                                 ! to this processor.

  logical, save :: first_iteration = .TRUE.      ! Whether this is the first time this
                                                 ! routine has been called.




  ! --------------------------------------------------------------------
  ! Check that the requested non-local functional is implemented.

  if ( vdW_type /= 1 .AND. vdW_type /= 2) call errore('xc_vdW_DF','E^nl_c not implemented',1)


  ! --------------------------------------------------------------------
  ! Write out the vdW-DF imformation.

  if ( ionode .AND. first_iteration ) call vdW_info
  first_iteration = .FALSE.


  ! --------------------------------------------------------------------
  ! Allocate arrays. nnr is a PWSCF variable that holds the number of
  ! points assigned to a given processor.

  allocate( q0(dfftp%nnr), dq0_drho(dfftp%nnr), dq0_dgradrho(dfftp%nnr), grad_rho(dfftp%nnr, 3) )
  allocate( total_rho(dfftp%nnr), potential(dfftp%nnr), thetas(dfftp%nnr, Nqs) )


  ! --------------------------------------------------------------------
  ! Add together the valence and core charge densities to get the total
  ! charge density. Note that rho_core is not the true core density and
  ! it is only non-zero for pseudopotentials with non-local core
  ! corrections.

  total_rho = rho_valence(:,1) + rho_core(:)


  ! --------------------------------------------------------------------
  ! Here we calculate the gradient in reciprocal space using FFT.

  call numerical_gradient (total_rho, grad_rho)


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
  ! the thetas array into the precursor to the u_i(k) array which is
  ! inverse fourier transformed to get the u_i(r) functions of SOLER
  ! equation 11. Add the energy we find to the output variable etxc.

  call vdW_energy (thetas, Ec_nl)
  etxc = etxc + Ec_nl

  if (iverbosity > 0) then
     call mp_sum(Ec_nl, intra_bgrp_comm)
     if (ionode) then
        write(stdout,'(/ / A)')       "     -----------------------------------------------"
        write(stdout,'(A, F15.8, A)') "     Non-local corr. energy    =  ", Ec_nl, " Ry"
        write(stdout,'(A /)')         "     -----------------------------------------------"
     end if
  end if


  ! --------------------------------------------------------------------
  ! Here we calculate the potential. This is calculated via equation 10
  ! of SOLER, using the u_i(r) calculated from quations 11 and 12 of
  ! SOLER. Each processor allocates the array to be the size of the full
  ! grid  because, as can be seen in SOLER equation 10, processors need
  ! to access grid points outside their allocated regions. Begin by
  ! FFTing the u_i(k) to get the u_i(r) of SOLER equation 11.

  do theta_i = 1, Nqs
     CALL invfft('Dense', thetas(:,theta_i), dfftp)
  end do

  call get_potential (q0, dq0_drho, dq0_dgradrho, grad_rho, thetas, potential)
  v(:,1) = v(:,1) + e2 * potential(:)


  ! --------------------------------------------------------------------
  ! The integral of rho(r)*potential(r) for the vtxc output variable.

  grid_cell_volume = omega/(dfftp%nr1*dfftp%nr2*dfftp%nr3)

  do i_grid = 1, dfftp%nnr
     vtxc = vtxc + e2 * grid_cell_volume * rho_valence(i_grid,1) * potential(i_grid)
  end do

  deallocate ( potential, q0, grad_rho, dq0_drho, dq0_dgradrho, total_rho, thetas )

  END SUBROUTINE xc_vdW_DF








  ! ####################################################################
  !                          |                  |
  !                          |  XC_VDW_DF_spin  |
  !                          |__________________|
  !
  ! This subroutine is as similar to xc_vdW_DF() as possible, but
  ! handles the collinear nspin=2 case.

  SUBROUTINE xc_vdW_DF_spin (rho_valence, rho_core, etxc, vtxc, v)

  USE gvect,                 ONLY : ngm, nl, g, nlm
  USE cell_base,             ONLY : omega, tpiba

  implicit none

  ! --------------------------------------------------------------------
  ! Local variables
  !                                              _
  real(dp), intent(IN) :: rho_valence(:,:)      !
  real(dp), intent(IN) :: rho_core(:)           !  PWSCF input variables
  real(dp), intent(inout) :: etxc, vtxc, v(:,:) !_


  integer :: i_grid, theta_i, i_proc            ! Indexing variables over grid points,
                                                ! theta functions, and processors, and a
                                                ! generic index.

  real(dp) :: grid_cell_volume                  ! The volume of the unit cell per G-grid point.

  real(dp), allocatable :: q0(:)                ! The saturated value of q (equations 11 and 12
                                                ! of DION). This saturation is that of
                                                ! equation 5 in SOLER.

  real(dp), allocatable :: grad_rho(:,:)        ! The gradient of the charge density. The
                                                ! format is as follows:
                                                ! grad_rho(grid_point, cartesian_component)
  real(dp), allocatable :: grad_rho_up(:,:)     ! The gradient of the up charge density. The
                                                ! format is as follows:
                                                ! grad_rho(grid_point, cartesian_component)
  real(dp), allocatable :: grad_rho_down(:,:)   ! The gradient of the down charge density. The
                                                ! format is as follows:
                                                ! grad_rho(grid_point, cartesian_component)

  real(dp), allocatable :: potential_up(:)      ! The vdW contribution to the potential.
  real(dp), allocatable :: potential_down(:)    ! The vdW contribution to the potential.

  real(dp), allocatable :: dq0_drho_up(:)       ! The derivative of the saturated q0
  real(dp), allocatable :: dq0_drho_down(:)     ! (equation 5 of SOLER) with respect
                                                ! to the charge density (see
                                                ! get_q0_on_grid subroutine for details).

  real(dp), allocatable :: dq0_dgradrho_up(:)   ! The derivative of the saturated q0
  real(dp), allocatable :: dq0_dgradrho_down(:) ! (equation 5 of SOLER) with respect
                                                ! to the gradient of the charge density
                                                ! (again, see get_q0_on_grid subroutine).

  complex(dp), allocatable :: thetas(:,:)       ! These are the functions of equation 8 of
                                                ! SOLER. They will be forward Fourier transformed
                                                ! in place to get theta(k) and worked on in
                                                ! place to get the u_alpha(r) of equation 11
                                                ! in SOLER. They are formatted as follows:
                                                ! thetas(grid_point, theta_i).

  real(dp) :: Ec_nl                             ! The non-local vdW contribution to the energy.

  real(dp), allocatable :: total_rho(:)         ! This is the sum of the valence (up and down)
                                                ! and core charge. This just holds the piece
                                                ! assigned to this processor.
  real(dp), allocatable :: rho_up(:)            ! This is the just the up valence charge.
                                                ! This just holds the piece assigned
                                                ! to this processor.
  real(dp), allocatable :: rho_down(:)          ! This is the just the down valence charge.
                                                ! This just holds the piece assigned
                                                ! to this processor.

  logical, save :: first_iteration = .TRUE.     ! Whether this is the first time this
                                                ! routine has been called.




  ! --------------------------------------------------------------------
  ! Check that the requested non-local functional is implemented.

  if ( vdW_type /= 1 .AND. vdW_type /= 2) call errore('xc_vdW_DF','E^nl_c not implemented',1)


  ! --------------------------------------------------------------------
  ! Write out the vdW-DF imformation.

  if ( ionode .AND. first_iteration ) call vdW_info
  first_iteration = .FALSE.


  ! --------------------------------------------------------------------
  ! Allocate arrays. nnr is a PWSCF variable that holds the number of
  ! points assigned to a given processor.

  allocate( q0(dfftp%nnr), total_rho(dfftp%nnr), grad_rho(dfftp%nnr, 3) )
  allocate( rho_up(dfftp%nnr), rho_down(dfftp%nnr) )
  allocate( dq0_drho_up (dfftp%nnr), dq0_dgradrho_up  (dfftp%nnr) )
  allocate( dq0_drho_down(dfftp%nnr), dq0_dgradrho_down(dfftp%nnr) )
  allocate( grad_rho_up(dfftp%nnr, 3), grad_rho_down(dfftp%nnr, 3) )
  allocate( potential_up(dfftp%nnr), potential_down(dfftp%nnr) )
  allocate( thetas(dfftp%nnr, Nqs) )


  ! --------------------------------------------------------------------
  ! Add together the valence and core charge densities to get the total
  ! charge density. Note that rho_core is not the true core density and
  ! it is only non-zero for pseudopotentials with non-local core
  ! corrections.

  rho_up    = rho_valence(:,1) + 0.5D0*rho_core(:)
  rho_down  = rho_valence(:,2) + 0.5D0*rho_core(:)
  total_rho = rho_up + rho_down


  ! --------------------------------------------------------------------
  ! Here we calculate the gradient in reciprocal space using FFT.

  call numerical_gradient ( total_rho, grad_rho      )
  call numerical_gradient ( rho_up,    grad_rho_up   )
  call numerical_gradient ( rho_down,  grad_rho_down )


  ! --------------------------------------------------------------------
  ! Find the value of q0 for all assigned grid points. q is defined in
  ! equations 11 and 12 of DION and q0 is the saturated version of q
  ! defined in equation 5 of SOLER. In the spin case, q0 is defined by
  ! equation 8 (and text above that equation) of THONHAUSER. This
  ! routine also returns the derivatives of the q0s with respect to the
  ! charge-density and the gradient of the charge-density. These are
  ! needed for the potential calculated below.

  CALL get_q0_on_grid_spin (total_rho, rho_up, rho_down, grad_rho, grad_rho_up, grad_rho_down, &
       q0, dq0_drho_up, dq0_drho_down, dq0_dgradrho_up, dq0_dgradrho_down, thetas)


  ! --------------------------------------------------------------------
  ! Carry out the integration in equation 8 of SOLER. This also turns
  ! the thetas array into the precursor to the u_i(k) array which is
  ! inverse fourier transformed to get the u_i(r) functions of SOLER
  ! equation 11. Add the energy we find to the output variable etxc.

  call vdW_energy(thetas, Ec_nl)
  etxc = etxc + Ec_nl

  if (iverbosity > 0) then
     call mp_sum(Ec_nl, intra_bgrp_comm)
     if (ionode) then
        write(stdout,'(/ / A)')       "     -----------------------------------------------"
        write(stdout,'(A, F15.8, A)') "     Non-local corr. energy    =  ", Ec_nl, " Ry"
        write(stdout,'(A /)')         "     -----------------------------------------------"
     end if
  end if


  ! --------------------------------------------------------------------
  ! Here we calculate the potential. This is calculated via equation 10
  ! of SOLER, using the u_i(r) calculated from quations 11 and 12 of
  ! SOLER. Each processor allocates the array to be the size of the full
  ! grid because, as can be seen in SOLER equation 10, processors need
  ! to access grid points outside their allocated regions. Begin by
  ! FFTing the u_i(k) to get the u_i(r) of SOLER equation 11.

  do theta_i = 1, Nqs
     CALL invfft('Dense', thetas(:,theta_i), dfftp)
  end do

  call get_potential (q0, dq0_drho_up  , dq0_dgradrho_up  , grad_rho_up  , thetas, potential_up  )
  call get_potential (q0, dq0_drho_down, dq0_dgradrho_down, grad_rho_down, thetas, potential_down)
  v(:,1) = v(:,1) + e2 * potential_up  (:)
  v(:,2) = v(:,2) + e2 * potential_down(:)


  ! --------------------------------------------------------------------
  ! The integral of rho(r)*potential(r) for the vtxc output variable

  grid_cell_volume = omega/(dfftp%nr1*dfftp%nr2*dfftp%nr3)

  do i_grid = 1, dfftp%nnr
     vtxc = vtxc + e2 * grid_cell_volume * rho_valence(i_grid,1) * potential_up  (i_grid) &
                 + e2 * grid_cell_volume * rho_valence(i_grid,2) * potential_down(i_grid)
  end do


  deallocate( potential_up, potential_down, q0, grad_rho, grad_rho_up, &
              grad_rho_down, dq0_drho_up, dq0_drho_down, thetas, &
              dq0_dgradrho_up, dq0_dgradrho_down, total_rho, rho_up, rho_down )

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

  implicit none

  real(dp),  intent(IN)     :: total_rho(:), grad_rho(:,:)         ! Input variables needed

  real(dp),  intent(OUT)    :: q0(:), dq0_drho(:), dq0_dgradrho(:) ! Output variables that have been allocated
                                                                   ! outside this routine but will be set here.
  complex(dp), intent(inout):: thetas(:,:)                         ! The thetas from SOLER.

  integer,   parameter      :: m_cut = 12                          ! How many terms to include in the sum
                                                                   ! of SOLER equation 5.

  real(dp)                  :: rho                                 ! Local variable for the density.
  real(dp)                  :: r_s                                 ! Wigner–Seitz radius.
  real(dp)                  :: s                                   ! Reduced gradient.
  real(dp)                  :: q, ec
  real(dp)                  :: dq0_dq                              ! The derivative of the saturated
                                                                   ! q0 with respect to q.

  integer                   :: i_grid, idx                         ! Indexing variables.




  ! --------------------------------------------------------------------
  ! Initialize q0-related arrays.

  q0(:)           = q_cut
  dq0_drho(:)     = 0.0D0
  dq0_dgradrho(:) = 0.0D0

  do i_grid = 1, dfftp%nnr

     rho = total_rho(i_grid)


     ! -----------------------------------------------------------------
     ! This prevents numerical problems. If the charge density is
     ! negative (an unphysical situation), we simply treat it as very
     ! small. In that case, q0 will be very large and will be saturated.
     ! For a saturated q0 the derivative dq0_dq will be 0 so we set q0 =
     ! q_cut and dq0_drho = dq0_dgradrho = 0 and go on to the next
     ! point.

     if ( rho < epsr ) cycle


     ! -----------------------------------------------------------------
     ! Calculate some intermediate values needed to find q.

     r_s = ( 3.0D0 / (4.0D0*pi*rho) )**(1.0D0/3.0D0)

     s   = sqrt( grad_rho(i_grid,1)**2 + grad_rho(i_grid,2)**2 + grad_rho(i_grid,3)**2 ) &
         / (2.0D0 * kF(rho) * rho )


     ! -----------------------------------------------------------------
     ! This is the q value defined in equations 11 and 12 of DION.
     ! Use pw() from flib/functionals.f90 to get qc = kf/eps_x * eps_c.

     call pw(r_s, 1, ec, dq0_drho(i_grid))
     q = -4.0D0*pi/3.0D0 * ec + kF(rho) * Fs(s)


     ! -----------------------------------------------------------------
     ! Bring q into its proper bounds.

     CALL saturate_q ( q, q_cut, q0(i_grid), dq0_dq )
     if (q0(i_grid) < q_min) q0(i_grid) = q_min


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

  end do


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

  do i_grid = 1, dfftp%nnr
     thetas(i_grid,:) = thetas(i_grid,:) * total_rho(i_grid)
  end do

  do idx = 1, Nqs
     CALL fwfft ('Dense', thetas(:,idx), dfftp)
  end do

  END SUBROUTINE get_q0_on_grid








  ! ####################################################################
  !                       |                       |
  !                       |  GET_Q0_ON_GRID_spin  |
  !                       |_______________________|

  SUBROUTINE get_q0_on_grid_spin (total_rho, rho_up, rho_down, grad_rho, &
             grad_rho_up, grad_rho_down, q0, dq0_drho_up, dq0_drho_down, &
             dq0_dgradrho_up, dq0_dgradrho_down, thetas)

  implicit none

  real(dp),  intent(IN)      :: total_rho(:), grad_rho(:,:)              ! Input variables.
  real(dp),  intent(IN)      :: rho_up(:), grad_rho_up(:,:)              ! Input variables.
  real(dp),  intent(IN)      :: rho_down(:), grad_rho_down(:,:)          ! Input variables.

  real(dp),  intent(OUT)     :: q0(:), dq0_drho_up(:), dq0_drho_down(:)  ! Output variables.
  real(dp),  intent(OUT)     :: dq0_dgradrho_up(:), dq0_dgradrho_down(:) ! Output variables.
  complex(dp), intent(inout) :: thetas(:,:)                              ! The thetas from SOLER.

  real(dp)                   :: rho, up, down                            ! Local copy of densities.
  real(dp)                   :: zeta                                     ! Spin polarization.
  real(dp)                   :: r_s                                      ! Wigner-Seitz radius.
  real(dp)                   :: q, qc, qx, qx_up, qx_down                ! q for exchange and correlation.
  real(dp)                   :: q0x_up, q0x_down                         ! Saturated q values.
  real(dp)                   :: ec, fac
  real(dp)                   :: dq0_dq, dq0x_up_dq, dq0x_down_dq         ! Derivative of q0 w.r.t q.
  real(dp)                   :: dqc_drho_up, dqc_drho_down               ! Intermediate values.
  real(dp)                   :: dqx_drho_up, dqx_drho_down               ! Intermediate values.
  real(dp)                   :: s_up, s_down                             ! Reduced gradients.
  integer                    :: i_grid, idx                              ! Indexing variables
  logical                    :: calc_qx_up, calc_qx_down




  fac = 2.0D0**(-1.0D0/3.0D0)


  ! --------------------------------------------------------------------
  ! Initialize q0-related arrays.

  q0(:)                = q_cut
  dq0_drho_up(:)       = 0.0D0
  dq0_drho_down(:)     = 0.0D0
  dq0_dgradrho_up(:)   = 0.0D0
  dq0_dgradrho_down(:) = 0.0D0


  do i_grid = 1, dfftp%nnr

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

     if ( rho < epsr ) cycle

     calc_qx_up   = .TRUE.
     calc_qx_down = .TRUE.

     if ( up   < epsr/2.0D0 ) calc_qx_up   = .FALSE.
     if ( down < epsr/2.0D0 ) calc_qx_down = .FALSE.


     ! -----------------------------------------------------------------
     ! The spin case is numerically even more tricky and we have to
     ! saturate each spin channel separately. Note that we are
     ! saturating at a higher value here, so that very large q values
     ! get saturated to exactly q_cut in the second, overall saturation.

     q0x_up        = 0.0D0
     q0x_down      = 0.0D0
     dqx_drho_up   = 0.0D0
     dqx_drho_down = 0.0D0


     if (calc_qx_up) then
        s_up    = sqrt( grad_rho_up(i_grid,1)**2 + grad_rho_up(i_grid,2)**2 + &
                  grad_rho_up(i_grid,3)**2 ) / (2.0D0 * kF(up) * up)
        qx_up   = kF(2.0D0*up) * Fs(fac*s_up)
        CALL saturate_q (qx_up, 4.0D0*q_cut, q0x_up, dq0x_up_dq)
     end if

     if (calc_qx_down) then
        s_down  = sqrt( grad_rho_down(i_grid,1)**2 + grad_rho_down(i_grid,2)**2 + &
                  grad_rho_down(i_grid,3)**2) / (2.0D0 * kF(down) * down)
        qx_down = kF(2.0D0*down) * Fs(fac*s_down)
        CALL saturate_q (qx_down, 4.0D0*q_cut, q0x_down, dq0x_down_dq)
     end if


     ! -----------------------------------------------------------------
     ! This is the q value defined in equations 11 and 12 of DION and
     ! equation 8 of THONHAUSER (also see text above that equation).

     r_s  = ( 3.0D0 / (4.0D0*pi*rho) )**(1.0D0/3.0D0)
     zeta = (up - down) / rho
     IF ( ABS(zeta) > 1.0D0 ) zeta = SIGN(1.0D0, zeta)
     call pw_spin(r_s, zeta, ec, dqc_drho_up, dqc_drho_down)

     qx = ( up * q0x_up + down * q0x_down ) / rho
     qc = -4.0D0*pi/3.0D0 * ec
     q  = qx + qc


     ! -----------------------------------------------------------------
     ! Bring q into its proper bounds.

     CALL saturate_q (q, q_cut, q0(i_grid), dq0_dq)
     if (q0(i_grid) < q_min) q0(i_grid) = q_min


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

     if (calc_qx_up) then
        dqx_drho_up   = 2.0D0*dq0x_up_dq*up*dqx_drho(2.0D0*up, fac*s_up) + q0x_up*down/rho
        dq0_dgradrho_up (i_grid) = 2.0D0 * dq0_dq * dq0x_up_dq * up * kF(2.0D0*up) * &
                        dFs_ds(fac*s_up) * ds_dgradrho(2.0D0*up)
     end if

     if (calc_qx_down) then
        dqx_drho_down = 2.0D0*dq0x_down_dq*down*dqx_drho(2.0D0*down, fac*s_down) + q0x_down*up/rho
        dq0_dgradrho_down(i_grid) = 2.0D0 * dq0_dq * dq0x_down_dq * down * kF(2.0D0*down) * &
                        dFs_ds(fac*s_down) * ds_dgradrho(2.0D0*down)
     end if

     if (calc_qx_down) dqx_drho_up   = dqx_drho_up   - q0x_down*down/rho
     if (calc_qx_up)   dqx_drho_down = dqx_drho_down - q0x_up  *up  /rho

     dqc_drho_up   = -4.0D0*pi/3.0D0 * (dqc_drho_up   - ec)
     dqc_drho_down = -4.0D0*pi/3.0D0 * (dqc_drho_down - ec)

     dq0_drho_up  (i_grid) = dq0_dq * (dqc_drho_up   + dqx_drho_up  )
     dq0_drho_down(i_grid) = dq0_dq * (dqc_drho_down + dqx_drho_down)

  end do


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

  do i_grid = 1, dfftp%nnr
     thetas(i_grid,:) = thetas(i_grid,:) * total_rho(i_grid)
  end do

  do idx = 1, Nqs
     CALL fwfft ('Dense', thetas(:,idx), dfftp)
  end do

  END SUBROUTINE get_q0_on_grid_spin








  ! ####################################################################
  !                            |              |
  !                            |  saturate_q  |
  !                            |______________|

  SUBROUTINE saturate_q (q, q_cut, q0, dq0_dq)

  implicit none

  real(dp),  intent(IN)      :: q             ! Input q.
  real(dp),  intent(IN)      :: q_cut         ! Cutoff q.
  real(dp),  intent(OUT)     :: q0            ! Output saturated q.
  real(dp),  intent(OUT)     :: dq0_dq        ! Derivative of dq0/dq.

  integer,    parameter      :: m_cut = 12    ! How many terms to include
                                              ! in the sum of SOLER equation 5.

  real(dp)                   :: e             ! Exponent.
  integer                    :: idx           ! Indexing variable.




  ! --------------------------------------------------------------------
  ! Here, we calculate q0 by saturating q according to equation 5 of
  ! SOLER. Also, we find the derivative dq0_dq needed for the
  ! derivatives dq0_drho and dq0_dgradrh0 discussed below.

  e      = 0.0D0
  dq0_dq = 0.0D0

  do idx = 1, m_cut
     e      = e + (q/q_cut)**idx/idx
     dq0_dq = dq0_dq + (q/q_cut)**(idx-1)
  end do

  q0     = q_cut*(1.0D0 - exp(-e))
  dq0_dq = dq0_dq * exp(-e)

  END SUBROUTINE saturate_q








  ! ####################################################################
  !                            |             |
  !                            | VDW_ENERGY  |
  !                            |_____________|
  !
  ! This routine carries out the integration of equation 8 of SOLER.  It
  ! returns the non-local exchange-correlation energy and the u_alpha(k)
  ! arrays used to find the u_alpha(r) arrays via equations 11 and 12 in
  ! SOLER.

  SUBROUTINE vdW_energy (thetas, vdW_xc_energy)

  USE gvect,           ONLY : nl, nlm, gg, ngm, igtongl, gl, ngl, gstart
  USE cell_base,       ONLY : tpiba, omega

  implicit none

  complex(dp), intent(inout) :: thetas(:,:)            ! On input this variable holds the theta
                                                       ! functions (equation 8, SOLER) in the format
                                                       ! thetas(grid_point, theta_i). On output
                                                       ! this array holds u_alpha(k) =
                                                       ! Sum_j[theta_beta(k)phi_alpha_beta(k)]

  real(dp), intent(out) :: vdW_xc_energy               ! The non-local correlation energy.

  real(dp), allocatable :: kernel_of_k(:,:)            ! This array will hold the interpolated kernel
                                                       ! values for each pair of q values in the q_mesh.

  real(dp)    :: g                                     ! The magnitude of the current g vector.
  integer     :: last_g                                ! The shell number of the last g vector.


  integer     :: g_i, q1_i, q2_i, i_grid               ! Index variables.

  complex(dp) :: theta(Nqs), thetam(Nqs), theta_g(Nqs) ! Temporary storage arrays used since we
                                                       ! are overwriting the thetas array here.
  real(dp)    :: G0_term, G_multiplier

  complex(dp), allocatable :: u_vdw(:,:)               ! Temporary array holding u_alpha(k).




  vdW_xc_energy = 0.0D0
  allocate (u_vdW(dfftp%nnr,Nqs), kernel_of_k(Nqs, Nqs))
  u_vdW(:,:) = CMPLX(0.0_DP,0.0_DP)


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
  if (gamma_only) G_multiplier = 2.0D0

  last_g = -1

  do g_i = 1, ngm

     if ( igtongl(g_i) .ne. last_g) then

        g = sqrt(gl(igtongl(g_i))) * tpiba
        call interpolate_kernel(g, kernel_of_k)
        last_g = igtongl(g_i)

     end if

     theta = thetas(nl(g_i),:)

     do q2_i = 1, Nqs
        do q1_i = 1, Nqs
           u_vdW(nl(g_i),q2_i)  = u_vdW(nl(g_i),q2_i) + kernel_of_k(q2_i,q1_i)*theta(q1_i)
        end do
        vdW_xc_energy = vdW_xc_energy + G_multiplier * (u_vdW(nl(g_i),q2_i)*conjg(theta(q2_i)))
     end do

     if (g_i < gstart ) vdW_xc_energy = vdW_xc_energy / G_multiplier

  end do

  if (gamma_only) u_vdW(nlm(:),:) = CONJG(u_vdW(nl(:),:))


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
  deallocate (u_vdW, kernel_of_k)

  END SUBROUTINE vdW_energy








  ! ####################################################################
  !                          |                 |
  !                          |  GET_POTENTIAL  |
  !                          |_________________|
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

  SUBROUTINE get_potential (q0, dq0_drho, dq0_dgradrho, grad_rho, u_vdW, potential)

  USE gvect,               ONLY : nl, g, nlm
  USE cell_base,           ONLY : alat, tpiba

  implicit none

  real(dp), intent(in) ::  q0(:), grad_rho(:,:)       ! Input arrays holding the value of q0 for
                                                      ! all points assigned to this processor and
                                                      ! the gradient of the charge density for
                                                      ! points assigned to this processor.

  real(dp), intent(in) :: dq0_drho(:), dq0_dgradrho(:)! The derivative of q0 with respect to the
                                                      ! charge density and gradient of the charge
                                                      ! density (almost). See comments in the
                                                      ! get_q0_on_grid subroutine above.

  complex(dp), intent(in) :: u_vdW(:,:)               ! The functions u_alpha(r) obtained by
                                                      ! inverse transforming the functions
                                                      ! u_alph(k). See equations 11 and 12 in SOLER

  real(dp), intent(inout) :: potential(:)             ! The non-local correlation potential for
                                                      ! points on the grid over the whole cell (not
                                                      ! just those assigned to this processor).

  real(dp), allocatable, save :: d2y_dx2(:,:)         ! Second derivatives of P_alpha polynomials
                                                      ! for interpolation.

  integer :: i_grid, P_i,icar                         ! Index variables.

  integer :: q_low, q_hi, q                           ! Variables to find the bin in the q_mesh that
                                                      ! a particular q0 belongs to (for interpolation).
  real(dp) :: dq, a, b, c, d, e, f                    ! Intermediate variables used in the
                                                      ! interpolation of the polynomials.

  real(dp) :: y(Nqs), dP_dq0, P                       ! The y values for a given polynomial (all 0
                                                      ! exept for element i of P_i) The derivative
                                                      ! of P at a given q0 and the value of P at a
                                                      ! given q0. Both of these are interpolated
                                                      ! below.

  real(dp) :: gradient2                               ! Squared gradient.

  real(dp)   , allocatable ::h_prefactor(:)
  complex(dp), allocatable ::h(:)




  allocate (h_prefactor(dfftp%nnr), h(dfftp%nnr))

  potential     = 0.0D0
  h_prefactor   = 0.0D0


  ! --------------------------------------------------------------------
  ! Get the second derivatives of the P_i functions for interpolation.
  ! We have already calculated this once but it is very fast and it's
  ! just as easy to calculate it again.

  if (.not. allocated( d2y_dx2) ) then

     allocate( d2y_dx2(Nqs, Nqs) )
     call initialize_spline_interpolation (q_mesh, d2y_dx2(:,:))

  end if


  do i_grid = 1, dfftp%nnr

     q_low = 1
     q_hi = Nqs


     ! -----------------------------------------------------------------
     ! Figure out which bin our value of q0 is in in the q_mesh.

     do while ( (q_hi - q_low) > 1)

        q = int((q_hi + q_low)/2)

        if (q_mesh(q) > q0(i_grid)) then
           q_hi = q
        else
           q_low = q
        end if

     end do

     if (q_hi == q_low) call errore('get_potential','qhi == qlow',1)

     dq = q_mesh(q_hi) - q_mesh(q_low)

     a = (q_mesh(q_hi) - q0(i_grid))/dq
     b = (q0(i_grid) - q_mesh(q_low))/dq
     c = (a**3 - a)*dq**2/6.0D0
     d = (b**3 - b)*dq**2/6.0D0
     e = (3.0D0*a**2 - 1.0D0)*dq/6.0D0
     f = (3.0D0*b**2 - 1.0D0)*dq/6.0D0

     do P_i = 1, Nqs
        y = 0.0D0
        y(P_i) = 1.0D0

        P      = a*y(q_low) + b*y(q_hi)  + c*d2y_dx2(P_i,q_low) + d*d2y_dx2(P_i,q_hi)
        dP_dq0 = (y(q_hi) - y(q_low))/dq - e*d2y_dx2(P_i,q_low) + f*d2y_dx2(P_i,q_hi)


        ! --------------------------------------------------------------
        ! The first term in equation 10 of SOLER.

        potential(i_grid) = potential(i_grid) + u_vdW(i_grid,P_i)* (P + dP_dq0 * dq0_drho(i_grid))
        if (q0(i_grid) .ne. q_mesh(Nqs)) then
           h_prefactor(i_grid) = h_prefactor(i_grid) + u_vdW(i_grid,P_i)*dP_dq0*dq0_dgradrho(i_grid)
        end if
     end do
  end do

  do icar = 1,3

     h(:) = CMPLX( h_prefactor(:) * grad_rho(:,icar), 0.0_DP )

     do i_grid = 1, dfftp%nnr
        gradient2 = grad_rho(i_grid,1)**2 + grad_rho(i_grid,2)**2 + grad_rho(i_grid,3)**2
        if ( gradient2 > 0.0D0 ) h(i_grid) = h(i_grid) / SQRT( gradient2 )
     end do

     CALL fwfft ('Dense', h, dfftp)
     h(nl(:)) = CMPLX(0.0_DP,1.0_DP) * tpiba * g(icar,:) * h(nl(:))
     if (gamma_only) h(nlm(:)) = CONJG(h(nl(:)))
     CALL invfft ('Dense', h, dfftp)
     potential(:) = potential(:) - REAL(h(:))

  end do

  deallocate (h_prefactor, h)

  END SUBROUTINE get_potential








  ! ####################################################################
  !                       |                        |
  !                       |  SPLINE_INTERPOLATION  |
  !                       |________________________|
  !
  ! This routine is modeled after an algorithm from "Numerical Recipes
  ! in C" by Cambridge University press, page 97.  It was adapted for
  ! Fortran, of course and for the problem at hand, in that it finds the
  ! bin a particular x value is in and then loops over all the P_i
  ! functions so we only have to find the bin once.

  SUBROUTINE spline_interpolation (x, evaluation_points, values)

  implicit none

  real(dp), intent(in) :: x(:), evaluation_points(:)     ! Input variables. The x values used to
                                                         ! form the interpolation (q_mesh in this
                                                         ! case) and the values of q0 for which we
                                                         ! are interpolating the function.

  complex(dp), intent(inout) :: values(:,:)              ! An output array (allocated outside this
                                                         ! routine) that stores the interpolated
                                                         ! values of the P_i (SOLER equation 3)
                                                         ! polynomials. The format is
                                                         ! values(grid_point, P_i).

  integer :: Ngrid_points, Nx                            ! Total number of grid points to evaluate
                                                         ! and input x points.

  real(dp), allocatable, save :: d2y_dx2(:,:)            ! The second derivatives required to do
                                                         ! the interpolation.

  integer :: i_grid, lower_bound, upper_bound, idx, P_i  ! Some indexing variables.

  real(dp), allocatable :: y(:)                          ! Temporary variables needed for the
  real(dp) :: a, b, c, d, dx                             ! interpolation.




  Nx = size(x)
  Ngrid_points = size(evaluation_points)


  ! --------------------------------------------------------------------
  ! Allocate the temporary array.

  allocate( y(Nx) )


  ! --------------------------------------------------------------------
  ! If this is the first time this routine has been called we need to
  ! get the second derivatives (d2y_dx2) required to perform the
  ! interpolations. So we allocate the array and call
  ! initialize_spline_interpolation to get d2y_dx2.

  if (.not. allocated(d2y_dx2) ) then

     allocate( d2y_dx2(Nx,Nx) )
     call initialize_spline_interpolation(x, d2y_dx2)

  end if

  do i_grid=1, Ngrid_points

     lower_bound = 1
     upper_bound = Nx

     do while ( (upper_bound - lower_bound) > 1 )

        idx = (upper_bound+lower_bound) / 2

        if ( evaluation_points(i_grid) > x(idx) ) then
           lower_bound = idx
        else
           upper_bound = idx
        end if

     end do

     dx = x(upper_bound)-x(lower_bound)

     a = (x(upper_bound) - evaluation_points(i_grid))/dx
     b = (evaluation_points(i_grid) - x(lower_bound))/dx
     c = ((a**3-a)*dx**2)/6.0D0
     d = ((b**3-b)*dx**2)/6.0D0

     do P_i = 1, Nx

        y = 0
        y(P_i) = 1

        values(i_grid, P_i) = a*y(lower_bound) + b*y(upper_bound) &
             + (c*d2y_dx2(P_i,lower_bound) + d*d2y_dx2(P_i, upper_bound))
     end do

  end do

  deallocate( y )

  END SUBROUTINE spline_interpolation








  ! ####################################################################
  !                  |                                   |
  !                  |  INITIALIZE_SPLINE_INTERPOLATION  |
  !                  |___________________________________|
  !
  ! This routine is modeled after an algorithm from "Numerical Recipes
  ! in C" by Cambridge University Press, pages 96-97. It was adapted
  ! for Fortran and for the problem at hand.

  SUBROUTINE initialize_spline_interpolation (x, d2y_dx2)

  implicit none

  real(dp), intent(in)  :: x(:)                 ! The input abscissa values.
  real(dp), intent(inout) :: d2y_dx2(:,:)       ! The output array (allocated outside this routine)
                                                ! that holds the second derivatives required for
                                                ! interpolating the function.

  integer :: Nx, P_i, idx                       ! The total number of x points and some indexing
                                                ! variables.

  real(dp), allocatable :: temp_array(:), y(:)  ! Some temporary arrays required. y is the array
                                                ! that holds the funcion values (all either 0 or
                                                ! 1 here).

  real(dp) :: temp1, temp2                      ! Some temporary variables required.




  Nx = size(x)

  allocate( temp_array(Nx), y(Nx) )

  do P_i=1, Nx

     ! -----------------------------------------------------------------
     ! In the Soler method, the polynomicals that are interpolated are Kroneker
     ! delta funcions at a particular q point. So, we set all y values to 0
     ! except the one corresponding to the particular function P_i.

     y = 0.0D0
     y(P_i) = 1.0D0

     d2y_dx2(P_i,1) = 0.0D0
     temp_array(1) = 0.0D0

     do idx = 2, Nx-1

        temp1 = (x(idx)-x(idx-1))/(x(idx+1)-x(idx-1))
        temp2 = temp1 * d2y_dx2(P_i,idx-1) + 2.0D0
        d2y_dx2(P_i,idx) = (temp1-1.0D0)/temp2

        temp_array(idx) = (y(idx+1)-y(idx))/(x(idx+1)-x(idx)) &
             - (y(idx)-y(idx-1))/(x(idx)-x(idx-1))
        temp_array(idx) = (6.0D0*temp_array(idx)/(x(idx+1)-x(idx-1)) &
             - temp1*temp_array(idx-1))/temp2

     end do

     d2y_dx2(P_i,Nx) = 0.0D0

     do idx=Nx-1, 1, -1

        d2y_dx2(P_i,idx) = d2y_dx2(P_i,idx) * d2y_dx2(P_i,idx+1) + temp_array(idx)

     end do

  end do

  deallocate( temp_array, y)

  END SUBROUTINE initialize_spline_interpolation








  ! ####################################################################
  !                          |                    |
  !                          | INTERPOLATE_KERNEL |
  !                          |____________________|
  !
  ! This routine is modeled after an algorithm from "Numerical Recipes in C" by
  ! Cambridge University Press, page 97.  Adapted for Fortran and the problem at
  ! hand.  This function is used to find the Phi_alpha_beta needed for equations
  ! 8 and 11 of SOLER.

  SUBROUTINE interpolate_kernel (k, kernel_of_k)

  implicit none

  real(dp), intent(in) :: k                    ! Input value, the magnitude of the g-vector
                                               ! for the current point.

  real(dp), intent(inout) :: kernel_of_k(:,:)  ! An output array (allocated outside this routine)
                                               ! that holds the interpolated value of the kernel
                                               ! for each pair of q points (i.e. the phi_alpha_beta
                                               ! of the Soler method.

  integer :: q1_i, q2_i, k_i                   ! Indexing variables.

  real(dp) :: A, B, C, D                       ! Intermediate values for the interpolation.




  ! --------------------------------------------------------------------
  ! Check to make sure that the kernel table we have is capable of
  ! dealing with this value of k. If k is larger than
  ! Nr_points*2*pi/r_max then we can't perform the interpolation. In
  ! that case, a kernel file should be generated with a larger number of
  ! radial points.

  if ( k >= Nr_points*dk ) then

     write(*,'(A,F10.5,A,F10.5)') "k =  ", k, "     k_max =  ",Nr_points*dk
     call errore('interpolate kernel', 'k value requested is out of range',1)

  end if

  kernel_of_k = 0.0D0


  ! --------------------------------------------------------------------
  ! This integer division figures out which bin k is in since the kernel
  ! is set on a uniform grid.

  k_i = int(k/dk)


  ! --------------------------------------------------------------------
  ! Test to see if we are trying to interpolate a k that is one of the
  ! actual function points we have.  The value is just the value of the
  ! function in that case.

  if (mod(k,dk) == 0) then

     do q1_i = 1, Nqs
        do q2_i = 1, q1_i

           kernel_of_k(q1_i, q2_i) = kernel(k_i,q1_i, q2_i)
           kernel_of_k(q2_i, q1_i) = kernel(k_i,q2_i, q1_i)

        end do
     end do

     return

  end if


  ! --------------------------------------------------------------------
  ! If we are not on a function point then we carry out the
  ! interpolation.

  A = (dk*(k_i+1.0D0) - k)/dk
  B = (k - dk*k_i)/dk
  C = (A**3-A)*dk**2/6.0D0
  D = (B**3-B)*dk**2/6.0D0

  do q1_i = 1, Nqs
     do q2_i = 1, q1_i

        kernel_of_k(q1_i, q2_i) = A*kernel(k_i, q1_i, q2_i) + B*kernel(k_i+1, q1_i, q2_i) &
             +(C*d2phi_dk2(k_i, q1_i, q2_i) + D*d2phi_dk2(k_i+1, q1_i, q2_i))

        kernel_of_k(q2_i, q1_i) = kernel_of_k(q1_i, q2_i)

     end do
  end do

  END SUBROUTINE interpolate_kernel








  ! ####################################################################
  !                         |                 |
  !                         |  STRESS_VDW_DF  |
  !                         |_________________|

  SUBROUTINE stress_vdW_DF (rho_valence, rho_core, nspin, sigma)

  use gvect,           ONLY : ngm, nl, g, nlm
  USE cell_base,       ONLY : tpiba

  implicit none

  real(dp), intent(IN)     :: rho_valence(:,:)       !
  real(dp), intent(IN)     :: rho_core(:)            ! Input variables
  integer,  intent(IN)     :: nspin                  !
  real(dp), intent(inout)  :: sigma(3,3)             !

  real(dp), allocatable    :: grad_rho(:,:)          !
  real(dp), allocatable    :: total_rho(:)           ! Rho values

  real(dp), allocatable    :: q0(:)                  !
  real(dp), allocatable    :: dq0_drho(:)            ! Q-values
  real(dp), allocatable    :: dq0_dgradrho(:)        !

  complex(dp), allocatable :: thetas(:,:)            ! Thetas
  integer                  :: i_proc, theta_i, l, m

  real(dp)                 :: sigma_grad(3,3)
  real(dp)                 :: sigma_ker(3,3)




  ! --------------------------------------------------------------------
  ! Tests

  if ( nspin>1 ) call errore('stress_vdW_DF','spin polarized vdW stress not implemented', nspin)

  sigma(:,:)      = 0.0_DP
  sigma_grad(:,:) = 0.0_DP
  sigma_ker(:,:)  = 0.0_DP


  ! --------------------------------------------------------------------
  ! Allocations

  allocate( grad_rho(dfftp%nnr, 3) )
  allocate( total_rho(dfftp%nnr) )
  allocate( q0(dfftp%nnr) )
  allocate( dq0_drho(dfftp%nnr), dq0_dgradrho(dfftp%nnr) )
  allocate( thetas(dfftp%nnr, Nqs) )


  ! --------------------------------------------------------------------
  ! Charge

  total_rho = rho_valence(:,1) + rho_core(:)


  ! --------------------------------------------------------------------
  ! Here we calculate the gradient in reciprocal space using FFT.

  call numerical_gradient (total_rho, grad_rho)


  ! --------------------------------------------------------------------
  ! Get q0.

  CALL get_q0_on_grid (total_rho, grad_rho, q0, dq0_drho, dq0_dgradrho, thetas)


  ! --------------------------------------------------------------------
  ! Stress

  CALL stress_vdW_DF_gradient (total_rho, grad_rho, q0, dq0_drho, dq0_dgradrho, thetas, sigma_grad)
  CALL stress_vdW_DF_kernel   (total_rho, q0, thetas, sigma_ker)

  sigma = - (sigma_grad + sigma_ker)

  do l = 1, 3
     do m = 1, l - 1
        sigma (m, l) = sigma (l, m)
     enddo
  enddo

  deallocate( grad_rho, total_rho, q0, dq0_drho, dq0_dgradrho, thetas )

  END SUBROUTINE stress_vdW_DF








  ! ####################################################################
  !                     |                          |
  !                     |  STRESS_VDW_DF_GRADIENT  |
  !                     |__________________________|

  SUBROUTINE stress_vdW_DF_gradient (total_rho, grad_rho, q0, dq0_drho, &
                                     dq0_dgradrho, thetas, sigma)

  USE gvect,                 ONLY : ngm, nl, g, nlm, nl, gg, igtongl, &
                                    gl, ngl, gstart
  USE cell_base,             ONLY : omega, tpiba, alat, at, tpiba2

  implicit none

  real(dp), intent(IN)     :: total_rho(:)           !
  real(dp), intent(IN)     :: grad_rho(:, :)         ! Input variables.
  real(dp), intent(inout)  :: sigma(:,:)             !
  real(dp), intent(IN)     :: q0(:)                  !
  real(dp), intent(IN)     :: dq0_drho(:)            !
  real(dp), intent(IN)     :: dq0_dgradrho(:)        !
  complex(dp), intent(IN)  :: thetas(:,:)            !

  complex(dp), allocatable :: u_vdW(:,:)             !

  real(dp), allocatable    :: d2y_dx2(:,:)           !
  real(dp) :: y(Nqs), dP_dq0, P, a, b, c, d, e, f    ! Interpolation.
  real(dp) :: dq                                     !

  integer  :: q_low, q_hi, q, q1_i, q2_i , g_i       ! Loop and q-points.

  integer  :: l, m
  real(dp) :: prefactor                              ! Final summation of sigma.
  real(dp) :: grad2                                  ! Magnitude of density gradient.

  integer  :: i_proc, theta_i, i_grid, q_i, &        !
              ix, iy, iz                             ! Iterators.

  character(LEN=1) :: intvar




  allocate( d2y_dx2(Nqs, Nqs) )
  allocate( u_vdW(dfftp%nnr, Nqs) )

  sigma(:,:) = 0.0_DP
  prefactor  = 0.0_DP


  ! --------------------------------------------------------------------
  ! Get u in k-space.

  call thetas_to_uk(thetas, u_vdW)


  ! --------------------------------------------------------------------
  ! Get u in real space.

  do theta_i = 1, Nqs
     CALL invfft('Dense', u_vdW(:,theta_i), dfftp)
  end do


  ! --------------------------------------------------------------------
  ! Get the second derivatives for interpolating the P_i.

  call initialize_spline_interpolation(q_mesh, d2y_dx2(:,:))


  ! --------------------------------------------------------------------
  ! Do the real space integration to obtain the stress component.

  do i_grid = 1, dfftp%nnr

     if ( total_rho(i_grid) < epsr ) cycle

     q_low = 1
     q_hi  = Nqs
     grad2 = sqrt( grad_rho(i_grid,1)**2 + grad_rho(i_grid,2)**2 + grad_rho(i_grid,3)**2 )

     if ( grad2 == 0.0_dp ) cycle
     ! -----------------------------------------------------------------
     ! Figure out which bin our value of q0 is in the q_mesh.

     do while ( (q_hi - q_low) > 1)

        q = int((q_hi + q_low)/2)

        if (q_mesh(q) > q0(i_grid)) then
            q_hi = q
        else
            q_low = q
        end if

     end do

     if (q_hi == q_low) call errore('stress_vdW_gradient','qhi == qlow',1)

     dq = q_mesh(q_hi) - q_mesh(q_low)

     a = (q_mesh(q_hi) - q0(i_grid))/dq
     b = (q0(i_grid) - q_mesh(q_low))/dq
     c = (a**3 - a)*dq**2/6.0D0
     d = (b**3 - b)*dq**2/6.0D0
     e = (3.0D0*a**2 - 1.0D0)*dq/6.0D0
     f = (3.0D0*b**2 - 1.0D0)*dq/6.0D0

     do q_i = 1, Nqs

        y(:)   = 0.0D0
        y(q_i) = 1.0D0

        dP_dq0 = (y(q_hi) - y(q_low))/dq - e*d2y_dx2(q_i,q_low) + f*d2y_dx2(q_i,q_hi)

        prefactor = u_vdW(i_grid,q_i) * dP_dq0 * dq0_dgradrho(i_grid) / grad2

        do l = 1, 3
            do m = 1, l

                sigma (l, m) = sigma (l, m) -  e2 * prefactor * &
                               (grad_rho(i_grid,l) * grad_rho(i_grid,m))
            end do
        end do

     end do

  end do

#if defined(__MPI)
  call mp_sum(  sigma, intra_bgrp_comm )
#endif

  call dscal (9, 1.d0 / (dfftp%nr1 * dfftp%nr2 * dfftp%nr3), sigma, 1)

  deallocate( d2y_dx2, u_vdW )

  END SUBROUTINE stress_vdW_DF_gradient








  ! ####################################################################
  !                      |                          |
  !                      |  STRESS_VDW_DF_KERNEL    |
  !                      |__________________________|

  SUBROUTINE stress_vdW_DF_kernel (total_rho, q0, thetas, sigma)

  USE gvect,                 ONLY : ngm, nl, g, nl, gg, igtongl, gl, ngl, gstart
  USE cell_base,             ONLY : omega, tpiba, tpiba2

  implicit none

  real(dp), intent(IN)    :: q0(:)
  real(dp), intent(IN)    :: total_rho(:)
  real(dp), intent(inout) :: sigma(3,3)
  complex(dp), intent(IN) :: thetas(:,:)

  real(dp), allocatable   :: dkernel_of_dk(:,:)

  integer                 :: l, m, q1_i, q2_i , g_i
  real(dp)                :: g2, ngmod2, g_kernel, G_multiplier
  integer                 :: last_g, theta_i




  allocate( dkernel_of_dk(Nqs, Nqs) )

  sigma(:,:) = 0.0_DP

  ! --------------------------------------------------------------------
  ! Integration in g-space.

  last_g = -1

  G_multiplier = 1.0D0

  if (gamma_only) G_multiplier = 2.0D0

  do g_i = gstart, ngm

     g2 = gg (g_i) * tpiba2
     g_kernel = sqrt(g2)

     if ( igtongl(g_i) .ne. last_g) then

        call interpolate_Dkernel_Dk(g_kernel, dkernel_of_dk)  ! Gets the derivatives.
        last_g = igtongl(g_i)

     end if

     do q2_i = 1, Nqs
        do q1_i = 1, Nqs
           do l = 1, 3
              do m = 1, l

              sigma (l, m) = sigma (l, m) - G_multiplier * 0.5 * e2 * thetas(nl(g_i),q1_i) * &
                             dkernel_of_dk(q1_i,q2_i)*conjg(thetas(nl(g_i),q2_i))* &
                             (g (l, g_i) * g (m, g_i) * tpiba2) / g_kernel
              end do
           end do
        enddo
     end do

     if (g_i < gstart ) sigma(:,:) = sigma(:,:) / G_multiplier

  enddo

#if defined(__MPI)
  call mp_sum(  sigma, intra_bgrp_comm )
#endif

  deallocate( dkernel_of_dk )

  END SUBROUTINE stress_vdW_DF_kernel








  ! ####################################################################
  !                        |                        |
  !                        | INTERPOLATE_DKERNEL_DK |
  !                        |________________________|


  SUBROUTINE interpolate_Dkernel_Dk (k, dkernel_of_dk)

  implicit none

  real(dp), intent(in) :: k                         ! Input value, the magnitude of the g-vector
                                                    ! for the current point.

  real(dp), intent(inout) :: dkernel_of_dk(Nqs,Nqs) ! An output array (allocated outside this
                                                    ! routine) that holds the interpolated value of
                                                    ! the kernel for each pair of q points (i.e. the
                                                    ! phi_alpha_beta of the Soler method.

  integer :: q1_i, q2_i, k_i                        ! Indexing variables.

  real(dp) :: A, B, dAdk, dBdk, dCdk, dDdk          ! Intermediate values for the interpolation.




  if ( k >= Nr_points*dk ) then

     write(*,'(A,F10.5,A,F10.5)') "k =  ", k, "     k_max =  ",Nr_points*dk
     call errore('interpolate kernel', 'k value requested is out of range',1)

  end if

  dkernel_of_dk = 0.0D0

  k_i = int(k/dk)

  A = (dk*(k_i+1.0D0) - k)/dk
  B = (k - dk*k_i)/dk

  dAdk = -1.0D0/dk
  dBdk = 1.0D0/dk
  dCdk = -((3*A**2 -1.0D0)/6.0D0)*dk
  dDdk = ((3*B**2 -1.0D0)/6.0D0)*dk

  do q1_i = 1, Nqs
     do q2_i = 1, q1_i

        dkernel_of_dk(q1_i, q2_i) = dAdk*kernel(k_i, q1_i, q2_i) + dBdk*kernel(k_i+1, q1_i, q2_i) &
             + dCdk*d2phi_dk2(k_i, q1_i, q2_i) + dDdk*d2phi_dk2(k_i+1, q1_i, q2_i)

        dkernel_of_dk(q2_i, q1_i) = dkernel_of_dk(q1_i, q2_i)

     end do
  end do

  END SUBROUTINE interpolate_Dkernel_Dk








  ! ####################################################################
  !                       |                       |
  !                       |   NUMERICAL_GRADIENT  |
  !                       |_______________________|
  !
  ! Calculates the gradient of the charge density numerically on the
  ! grid.  We use the PWSCF gradient style.

  SUBROUTINE numerical_gradient (total_rho, grad_rho)

  USE gvect,             ONLY : ngm, nl, g, nlm
  USE cell_base,         ONLY : tpiba

  implicit none

  real(dp), intent(in) :: total_rho(:)    ! Input array holding total charge density.

  real(dp), intent(out) :: grad_rho(:,:)  ! Output array that will holds the gradient
                                          ! of the charge density.

  integer :: icar                         ! Counter on cartesian components.
  complex(dp), allocatable :: c_rho(:)    ! auxiliary complex array for rho.
  complex(dp), allocatable :: c_grho(:)   ! auxiliary complex array for grad rho.




  ! --------------------------------------------------------------------
  ! Rho in G space.

  allocate ( c_rho(dfftp%nnr), c_grho(dfftp%nnr) )
  c_rho(1:dfftp%nnr) = CMPLX(total_rho(1:dfftp%nnr),0.0_DP)
  CALL fwfft ('Dense', c_rho, dfftp)

  do icar=1,3

     ! -----------------------------------------------------------------
     ! Compute gradient in G space.

     c_grho(:) =CMPLX(0.0_DP,0.0_DP)
     c_grho(nl(:)) = CMPLX (0.0_DP,1.0_DP) * tpiba * g(icar,:) * c_rho(nl(:))
     if (gamma_only) c_grho( nlm(:) ) = CONJG( c_grho( nl(:) ) )


     ! -----------------------------------------------------------------
     ! Back in real space.

     CALL invfft ('Dense', c_grho, dfftp)
     grad_rho(:,icar) = REAL( c_grho(:) )

  end do

  deallocate ( c_rho, c_grho )

  END SUBROUTINE numerical_gradient








  ! ####################################################################
  !                          |              |
  !                          | thetas_to_uk |
  !                          |______________|


  SUBROUTINE thetas_to_uk (thetas, u_vdW)

  USE gvect,           ONLY : nl, nlm, gg, ngm, igtongl, gl, ngl, gstart
  USE cell_base,       ONLY : tpiba, omega

  implicit none

  complex(dp), intent(in)  :: thetas(:,:)       ! On input this variable holds the theta functions
                                                ! (equation 8, SOLER) in the format
                                                ! thetas(grid_point, theta_i).
  complex(dp), intent(out) :: u_vdW(:,:)        ! On output this array holds u_alpha(k) =
                                                ! Sum_j[theta_beta(k)phi_alpha_beta(k)].

  real(dp), allocatable    :: kernel_of_k(:,:)  ! This array will hold the interpolated kernel
                                                ! values for each pair of q values in the q_mesh.

  real(dp) :: g
  integer :: last_g, g_i, q1_i, q2_i, i_grid    ! Index variables.

  complex(dp) :: theta(Nqs)                     ! Temporary storage vector used since we are
                                                ! overwriting the thetas array here.




  allocate( kernel_of_k(Nqs, Nqs) )

  u_vdW(:,:) = CMPLX(0.0_DP,0.0_DP)

  last_g = -1

  do g_i = 1, ngm

     if ( igtongl(g_i) .ne. last_g) then

        g = sqrt(gl(igtongl(g_i))) * tpiba
        call interpolate_kernel(g, kernel_of_k)
        last_g = igtongl(g_i)

     end if

     theta = thetas(nl(g_i),:)

     do q2_i = 1, Nqs
        do q1_i = 1, Nqs
           u_vdW(nl(g_i),q2_i) = u_vdW(nl(g_i),q2_i) + kernel_of_k(q2_i,q1_i)*theta(q1_i)
        end do
     end do

  end do

  if (gamma_only) u_vdW(nlm(:),:) = CONJG(u_vdW(nl(:),:))

  deallocate( kernel_of_k )

  END SUBROUTINE thetas_to_uk








  ! ####################################################################
  !                          |            |
  !                          |  VDW_INFO  |
  !                          |____________|

  SUBROUTINE vdW_info

  implicit none

  integer :: I




  ! --------------------------------------------------------------------
  ! Here we output some of the parameters being used in the run. This is
  ! important because these parameters are read from the
  ! vdW_kernel_table file. The user should ensure that these are the
  ! parameters they were intending to use on each run.

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

  if (iverbosity > 0) then
  WRITE(stdout,'(5x,"Carrying out vdW-DF run using the following parameters:")')
  WRITE(stdout,'(5X,A,I3,A,I5,A,F8.3)' ) "Nqs    = ",Nqs,"  Npoints = ",Nr_points,"  r_max = ",r_max
  WRITE(stdout,'(5X,"q_mesh =",4F12.8)') (q_mesh(I), I=1, 4)
  WRITE(stdout,'(13X,4F12.8)') (q_mesh(I), I=5, Nqs)
  end if

  END SUBROUTINE


END MODULE vdW_DF
