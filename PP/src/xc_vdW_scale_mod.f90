!-----------------------------------------------------------------------
! Program written by Yang Jiao,  Oct 2016, GPL, No warranties.
! 
! This module is adapted from MODULE vdW_DF
! ----------------------------------------------------------------------
MODULE vdW_DF_scale
!-----------------------------------------------------------------------
!
! ... This module calculates the coupling constant dependancy of
!     non-local correlation 
!     E_{c,\lambda}^{nl},    \lambda \in [0:1]
! ... Eq. (48) of PRB 97, 085115 (2018).
!
!
!----------------------------------------
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


CONTAINS
  ! ####################################################################
  !                           |                 |
  !                           |  XC_VDW_DF_ncc  |
  !                           |_________________|

  SUBROUTINE xc_vdW_DF_ncc(cc,lecnl_qx,etcnlccc) 

  USE gvect,                 ONLY : ngm, g
  USE cell_base,             ONLY : omega, tpiba
  USE scf,                   ONLY : rho, rho_core
  USE vdW_DF,                ONLY : vdw_type,vdW_energy

  implicit none

  ! --------------------------------------------------------------------
  ! Local variables
  !                                               _
  real(dp), intent(IN) :: cc                     !
  logical, intent(IN)  :: lecnl_qx               !  PPACF input variables
  real(dp), intent(inout) :: etcnlccc            !_


  integer :: i_grid, theta_i, i_proc             ! Indexing variables over grid points,
                                                 ! theta functions, and processors.

  real(dp) :: grid_cell_volume                   ! The volume of the unit cell per G-grid point.

  real(dp), allocatable ::  q0(:)                ! The saturated value of q (equations 11 and 12
                                                 ! of DION). This saturation is that of
                                                 ! equation 5 in SOLER.

  real(dp), allocatable :: grad_rho(:,:)         ! The gradient of the charge density. The
                                                 ! format is as follows:
                                                 ! grad_rho(cartesian_component,grid_point).

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




  ! --------------------------------------------------------------------
  ! Check that the requested non-local functional is implemented.

  if ( vdW_type /= 1 .AND. vdW_type /= 2) call errore('xc_vdW_DF','E^nl_c not implemented',1)


  ! --------------------------------------------------------------------
  ! Allocate arrays. nnr is a PWSCF variable that holds the number of
  ! points assigned to a given processor.

  allocate( q0(dfftp%nnr), grad_rho(3,dfftp%nnr) )
  allocate( total_rho(dfftp%nnr), thetas(dfftp%nnr, Nqs) )


  ! --------------------------------------------------------------------
  ! Add together the valence and core charge densities to get the total
  ! charge density. Note that rho_core is not the true core density and
  ! it is only non-zero for pseudopotentials with non-local core
  ! corrections.

  total_rho = rho%of_r(:,1) + rho_core(:)


  ! --------------------------------------------------------------------
  ! Here we calculate the gradient in reciprocal space using FFT.

  call fft_gradient_r2r (dfftp, total_rho, g, grad_rho)


  ! --------------------------------------------------------------------
  ! Find the value of q0 for all assigned grid points. q is defined in
  ! equations 11 and 12 of DION and q0 is the saturated version of q
  ! defined in equation 5 of SOLER. This routine also returns the
  ! derivatives of the q0s with respect to the charge-density and the
  ! gradient of the charge-density. These are needed for the potential
  ! calculated below. This routine also calculates the thetas.

  CALL get_q0cc_on_grid (cc,lecnl_qx,total_rho, grad_rho, q0, thetas)


  ! --------------------------------------------------------------------
  ! Carry out the integration in equation 8 of SOLER. This also turns
  ! the thetas array into the precursor to the u_i(k) array which is
  ! inverse fourier transformed to get the u_i(r) functions of SOLER
  ! equation 11. Add the energy we find to the output variable etxc.

  call vdW_energy (thetas, Ec_nl)
  call mp_sum(Ec_nl, intra_bgrp_comm)
  etcnlccc = Ec_nl

  deallocate ( q0, grad_rho, total_rho, thetas )

  END SUBROUTINE xc_vdW_DF_ncc








  ! ####################################################################
  !                          |                      |
  !                          |  XC_VDW_DF_spin_ncc  |
  !                          |______________________|
  !
  ! This subroutine is as similar to xc_vdW_DF_ncc() as possible, but
  ! handles the collinear nspin=2 case.

  SUBROUTINE xc_vdW_DF_spin_ncc (cc,lecnl_qx,etcnlccc)

  USE gvect,                 ONLY : ngm, g
  USE cell_base,             ONLY : omega, tpiba
  USE scf,                   ONLY : rho, rho_core
  USE vdW_DF,                ONLY : vdw_type,vdW_energy

  implicit none

  ! --------------------------------------------------------------------
  ! Local variables
  !                                              _
  real(dp), intent(IN) :: cc                    !
  logical, intent(IN)  :: lecnl_qx              !  PPACF input variables
  real(dp), intent(inout) :: etcnlccc           !_


  integer :: i_grid, theta_i, i_proc            ! Indexing variables over grid points,
                                                ! theta functions, and processors, and a
                                                ! generic index.

  real(dp) :: grid_cell_volume                  ! The volume of the unit cell per G-grid point.

  real(dp), allocatable :: q0(:)                ! The saturated value of q (equations 11 and 12
                                                ! of DION). This saturation is that of
                                                ! equation 5 in SOLER.

  real(dp), allocatable :: grad_rho(:,:)        ! The gradient of the charge density. The
                                                ! format is as follows:
                                                ! grad_rho(cartesian_component,grid_point)
  real(dp), allocatable :: grad_rho_up(:,:)     ! The gradient of the up charge density.
                                                ! Same format as grad_rho
  real(dp), allocatable :: grad_rho_down(:,:)   ! The gradient of the down charge density.
                                                ! Same format as grad_rho

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




  ! --------------------------------------------------------------------
  ! Check that the requested non-local functional is implemented.

  if ( vdW_type /= 1 .AND. vdW_type /= 2) call errore('xc_vdW_DF','E^nl_c not implemented',1)



  ! --------------------------------------------------------------------
  ! Allocate arrays. nnr is a PWSCF variable that holds the number of
  ! points assigned to a given processor.

  allocate( q0(dfftp%nnr), total_rho(dfftp%nnr), grad_rho(3,dfftp%nnr) )
  allocate( rho_up(dfftp%nnr), rho_down(dfftp%nnr) )
  allocate( grad_rho_up(3,dfftp%nnr), grad_rho_down(3,dfftp%nnr) )
  allocate( thetas(dfftp%nnr, Nqs) )


  ! --------------------------------------------------------------------
  ! Add together the valence and core charge densities to get the total
  ! charge density. Note that rho_core is not the true core density and
  ! it is only non-zero for pseudopotentials with non-local core
  ! corrections.

  rho_up    = ( rho%of_r(:,1) + rho%of_r(:,2) + rho_core(:) ) * 0.5_dp
  rho_down  = ( rho%of_r(:,1) - rho%of_r(:,2) + rho_core(:) ) * 0.5_dp
  total_rho = rho_up + rho_down

#if defined (__SPIN_BALANCED)
     rho_up   = total_rho*0.5D0
     rho_down = rho_up
     write(stdout,'(/,/,"     Performing spin-balanced Ecnl calculation!")')
#endif

  ! --------------------------------------------------------------------
  ! Here we calculate the gradient in reciprocal space using FFT.

  call fft_gradient_r2r (dfftp, total_rho, g, grad_rho)
  call fft_gradient_r2r (dfftp, rho_up,    g, grad_rho_up)
  call fft_gradient_r2r (dfftp, rho_down,  g, grad_rho_down)

  ! --------------------------------------------------------------------
  ! Find the value of q0 for all assigned grid points. q is defined in
  ! equations 11 and 12 of DION and q0 is the saturated version of q
  ! defined in equation 5 of SOLER. In the spin case, q0 is defined by
  ! equation 8 (and text above that equation) of THONHAUSER. This
  ! routine also returns the derivatives of the q0s with respect to the
  ! charge-density and the gradient of the charge-density. These are
  ! needed for the potential calculated below.

  CALL get_q0cc_on_grid_spin (cc,lecnl_qx,total_rho, rho_up, rho_down, &
       grad_rho, grad_rho_up, grad_rho_down, q0, thetas)


  ! --------------------------------------------------------------------
  ! Carry out the integration in equation 8 of SOLER. This also turns
  ! the thetas array into the precursor to the u_i(k) array which is
  ! inverse fourier transformed to get the u_i(r) functions of SOLER
  ! equation 11. Add the energy we find to the output variable etxc.

  call vdW_energy(thetas, Ec_nl)
  call mp_sum(Ec_nl, intra_bgrp_comm)
  etcnlccc = Ec_nl

  deallocate( q0, grad_rho, grad_rho_up, grad_rho_down, thetas, &
              total_rho, rho_up, rho_down )

  END SUBROUTINE xc_vdW_DF_spin_ncc








  ! ####################################################################
  !                       |                    |
  !                       |  GET_Q0CC_ON_GRID  |
  !                       |____________________|
  !
  ! This routine calculates the scaling of the inverse length scale of 
  ! q value defined in (DION equations 11 and 12), then saturates it 
  ! according to (SOLER equation 5). More specifically it calculates 
  ! Eq. (45) of PRB 97, 085115 (2018).
  !

  SUBROUTINE get_q0cc_on_grid (cc,lecnl_qx,total_rho, grad_rho, q0, thetas)

  USE vdW_DF,               ONLY : saturate_q, Fs, kF
  USE vdW_DF,               ONLY : spline_interpolation

  implicit none

  real(dp),  intent(IN)     :: cc
  logical,   intent(IN)     :: lecnl_qx
  real(dp),  intent(IN)     :: total_rho(:), grad_rho(:,:)         ! Input variables needed

  real(dp),  intent(OUT)    :: q0(:)                               ! Output variables that have been allocated
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

  real(dp)                  :: dqc_drho
  integer                   :: i_grid, idx                         ! Indexing variables.




  ! --------------------------------------------------------------------
  ! Initialize q0-related arrays.

  q0(:)           = q_cut

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

     s   = sqrt( grad_rho(1,i_grid)**2 + grad_rho(2,i_grid)**2 + grad_rho(3,i_grid)**2 ) &
         / (2.0D0 * kF(rho) * rho )


     ! -----------------------------------------------------------------
     ! This is the q value defined in equations 11 and 12 of DION.
     ! Use pw() from flib/functionals.f90 to get qc = kf/eps_x * eps_c.

     call pw(cc*r_s, 1, ec, dqc_drho)
     q = -4.0D0*pi/3.0D0 * ec + kF(rho) * Fs(s)/cc
     if(lecnl_qx) then
       q=kF(rho) * Fs(s)/cc
     endif


     ! -----------------------------------------------------------------
     ! Bring q into its proper bounds.

     CALL saturate_q ( q*cc, q_cut, q0(i_grid), dq0_dq )
     if (q0(i_grid) < q_min) q0(i_grid) = q_min

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
     CALL fwfft ('Rho', thetas(:,idx), dfftp)
  end do

  END SUBROUTINE get_q0cc_on_grid








  ! ####################################################################
  !                       |                         |
  !                       |  GET_Q0CC_ON_GRID_spin  |
  !                       |_________________________|

  SUBROUTINE get_q0cc_on_grid_spin (cc,lecnl_qx,total_rho, rho_up, rho_down, grad_rho, &
             grad_rho_up, grad_rho_down, q0, thetas)

  USE vdW_DF,               ONLY : saturate_q, Fs, kF
  USE vdW_DF,               ONLY : spline_interpolation

  implicit none

  real(dp),  intent(IN)      :: cc                                       ! Input variables.
  logical,   intent(IN)      :: lecnl_qx                                 ! Input variables.
  real(dp),  intent(IN)      :: total_rho(:), grad_rho(:,:)              ! Input variables.
  real(dp),  intent(IN)      :: rho_up(:), grad_rho_up(:,:)              ! Input variables.
  real(dp),  intent(IN)      :: rho_down(:), grad_rho_down(:,:)          ! Input variables.

  real(dp),  intent(OUT)     :: q0(:)                                    ! Output variables.
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
        s_up    = sqrt( grad_rho_up(1,i_grid)**2 + grad_rho_up(2,i_grid)**2 + &
                  grad_rho_up(3,i_grid)**2 ) / (2.0D0 * kF(up) * up)
        qx_up   = kF(2.0D0*up) * Fs(fac*s_up)
        CALL saturate_q (qx_up, 4.0D0*q_cut, q0x_up, dq0x_up_dq)
     end if

     if (calc_qx_down) then
        s_down  = sqrt( grad_rho_down(1,i_grid)**2 + grad_rho_down(2,i_grid)**2 + &
                  grad_rho_down(3,i_grid)**2) / (2.0D0 * kF(down) * down)
        qx_down = kF(2.0D0*down) * Fs(fac*s_down)
        CALL saturate_q (qx_down, 4.0D0*q_cut, q0x_down, dq0x_down_dq)
     end if


     ! -----------------------------------------------------------------
     ! This is the q value defined in equations 11 and 12 of DION and
     ! equation 8 of THONHAUSER (also see text above that equation).

     r_s  = ( 3.0D0 / (4.0D0*pi*rho) )**(1.0D0/3.0D0)
     zeta = (up - down) / rho
     IF ( ABS(zeta) > 1.0D0 ) zeta = SIGN(1.0D0, zeta)
     call pw_spin(cc*r_s, zeta, ec, dqc_drho_up, dqc_drho_down)

     qx = ( up * q0x_up + down * q0x_down ) / rho / cc
     qc = -4.0D0*pi/3.0D0 * ec
     q  = qx + qc
     if(lecnl_qx) then
        q=qx
     endif


     ! -----------------------------------------------------------------
     ! Bring q into its proper bounds.

     CALL saturate_q (q*cc, q_cut, q0(i_grid), dq0_dq)
     if (q0(i_grid) < q_min) q0(i_grid) = q_min


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
     CALL fwfft ('Rho', thetas(:,idx), dfftp)
  end do

  END SUBROUTINE get_q0cc_on_grid_spin




END MODULE vdW_DF_scale
