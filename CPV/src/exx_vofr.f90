module multipole_expansion
#ifdef __CUDA
  use cudafor
#endif
  use kinds, only : dp
  implicit none
contains

#ifdef __CUDA
  attributes(device, host) &
#endif
    function get_plm(x, y, l, m) result (plm)
    !
    ! this subroutine calculates all of the associated Legendre polynomials up 
    ! to a supplied lpole (maximum lpole is 9). 
    !
    IMPLICIT NONE
    !
    ! INPUT VARIABLES
    ! cos(theta),sin(theta), for a given grid point
    real(dp), value :: x,y
    ! order of multipole expansion
    integer, value :: l, m
    real(dp) :: plm
    !
    ! WORK VARIABLES:
    ! powers of x, y: xn=x^n, yn=y^n
    REAL(dp) x2, x3, x4, x5, x6, x7, x8, x9
    REAL(dp) y2, y3, y4, y5, y6, y7, y8, y9
    !
    select case (l)
    case (0)
      plm = 1.d0
    case (1)
      if (m.eq.0) then
        plm = x
      else
        plm = -y
      end if ! m.eq.0
    case (2)
      x2 = x*x
      y2 = y*y
      select case (m)
      case (0)
        plm =  1.5d0*x2 - 0.5d0
      case (1)
        plm = -3.0d0*x*y
      case (2)
        plm =  3.0d0*y2
      end select ! m
    case (3)
      x2 = x*x
      y2 = y*y
      x3 = x2*x
      y3 = y2*y
      select case (m)
      case (0)
        plm =  2.5d0*x3 - 1.5d0*x
      case (1)
        plm = (-7.5d0*x2 + 1.5d0)*y
      case (2)
        plm =  15.0d0*x*y2
      case (3)
        plm = -15.0d0*y3
      end select ! m
    case (4)
      x2 = x*x
      y2 = y*y
      x3 = x2*x
      y3 = y2*y
      x4 = x2*x2
      y4 = y2*y2
      select case (m)
      case (0)
        plm =  4.375d0*x4 - 3.75d0*x2 + 0.375d0
      case (1)
        plm = (-17.5d0*x3 + 7.5d0*x)*y
      case (2)
        plm = ( 52.5d0*x2 - 7.5d0  )*y2
      case (3)
        plm = -105.0d0*x*y3
      case (4)
        plm =  105.0d0*y4
      end select ! m
    case (5)
      x2 = x*x
      y2 = y*y
      x3 = x2*x
      y3 = y2*y
      x4 = x2*x2
      y4 = y2*y2
      x5 = x3*x2
      y5 = y3*y2
      select case (m)
      case (0)
        plm =  7.875d0*x5 - 8.75d0*x3 + 1.875d0*x
      case (1)
        plm = (-39.375d0*x4 + 26.25d0*x2 - 1.875d0)*y
      case (2)
        plm = ( 157.5d0*x3 - 52.5d0*x)*y2
      case (3)
        plm = (-472.5d0*x2 + 52.5d0)  *y3
      case (4)
        plm =  945.0d0*x*y4
      case (5)
        plm = -945.0d0*y5
      end select ! m
    case (6)
      x2 = x*x
      y2 = y*y
      x3 = x2*x
      y3 = y2*y
      x4 = x2*x2
      y4 = y2*y2
      x5 = x3*x2
      y5 = y3*y2
      x6 = x3*x3
      y6 = y3*y3
      select case (m)
      case (0)
        plm = 14.4375d0*x6 - 19.6875d0*x4 + 6.5625d0*x2 - 0.3125d0
      case (1)
        plm = (-86.625d0*x5 + 78.75d0*x3 - 13.125d0*x)*y
      case (2)
        plm = ( 433.125d0*x4 - 236.25d0*x2 + 13.125d0)*y2
      case (3)
        plm = (-1732.5d0*x3 + 472.5d0*x            )*y3
      case (4)
        plm = ( 5197.5d0*x2 - 472.5d0              )*y4
      case (5)
        plm = -10395.0d0*x*y5
      case (6)
        plm =  10395.0d0*y6
      end select ! m
    case (7)
      x2 = x*x
      y2 = y*y
      x3 = x2*x
      y3 = y2*y
      x4 = x2*x2
      y4 = y2*y2
      x5 = x3*x2
      y5 = y3*y2
      x6 = x3*x3
      y6 = y3*y3
      x7 = x4*x3
      y7 = y4*y3
      select case (m)
      case (0)
        plm = 26.8125d0*x7 - 43.3125d0*x5 + 19.6875d0*x3 - 2.1875d0*x
      case (1)
        plm = -187.6875d0*x6*y + 216.5625d0*x4*y - 59.0625d0*x2*y + 2.1875d0*y
      case (2)
        plm = 1126.125d0*x5*y2 - 866.25d0*x3*y2 + 118.125d0*x*y2
      case (3)
        plm = -5630.625d0*x4*y3 + 2598.75d0*x2*y3 - 118.125d0*y3
      case (4)
        plm = 22522.5d0*x3*y4 - 5197.5d0*x*y4
      case (5)
        plm = -67567.5d0*x2*y5 + 5197.5d0*y5
      case (6)
        plm = 135135.0d0*x*y6
      case (7)
        plm = -135135.0d0*y7
      end select ! m
    case (8)
      x2 = x*x
      y2 = y*y
      x3 = x2*x
      y3 = y2*y
      x4 = x2*x2
      y4 = y2*y2
      x5 = x3*x2
      y5 = y3*y2
      x6 = x3*x3
      y6 = y3*y3
      x7 = x4*x3
      y7 = y4*y3
      x8 = x4*x4
      y8 = y4*y4
      select case (m)
      case (0)
        plm = 50.2734375d0*x8 - 93.84375d0*x6 + 54.140625d0*x4 - 9.84375d0*x2 + 0.2734375d0
      case (1)
        plm = -402.1875d0*x7*y + 563.0625d0*x5*y - 216.5625d0*x3*y + 19.6875d0*x*y
      case (2)
        plm = 2815.3125d0*x6*y2 - 2815.3125d0*x4*y2 + 649.6875d0*x2*y2 - 19.6875d0*y2
      case (3)
        plm = -16891.875d0*x5*y3 + 11261.25d0*x3*y3 - 1299.375d0*x*y3
      case (4)
        plm = 84459.375d0*x4*y4 - 33783.75d0*x2*y4 + 1299.375d0*y4
      case (5)
        plm = -337837.5d0*x3*y5 + 67567.5d0*x*y5
      case (6)
        plm = 1013512.5d0*x2*y6 - 67567.5d0*y6
      case (7)
        plm = -2027025.0d0*x*y7
      case (8)
        plm = 2027025.0d0*y8
      end select ! m
    case (9)
      x2 = x*x
      y2 = y*y
      x3 = x2*x
      y3 = y2*y
      x4 = x2*x2
      y4 = y2*y2
      x5 = x3*x2
      y5 = y3*y2
      x6 = x3*x3
      y6 = y3*y3
      x7 = x4*x3
      y7 = y4*y3
      x8 = x4*x4
      y8 = y4*y4
      y9 = y5*y4
      x9 = x5*x4
      select case (m)
      case (0)
        plm = 94.9609375d0*x9 - 201.09375d0*x7 + 140.765625d0*x5 - 36.09375d0*x3 + 2.4609375d0*x
      case (1)
        plm = -854.6484375d0*x8*y + 1407.65625d0*x6*y - 703.828125d0*x4*y + 108.28125d0*x2*y - 2.4609375d0*y
      case (2)
        plm = 6837.1875d0*x7*y2 - 8445.9375d0*x5*y2 + 2815.3125d0*x3*y2 - 216.5625d0*x*y2
      case (3)
        plm = -47860.3125d0*x6*y3 + 42229.6875d0*x4*y3 - 8445.9375d0*x2*y3 + 216.5625d0*y3
      case (4)
        plm = 287161.875d0*x5*y4 - 168918.75d0*x3*y4 + 16891.875d0*x*y4
      case (5)
        plm = -1435809.375d0*x4*y5 + 506756.25d0*x2*y5 - 16891.875d0*y5
      case (6)
        plm = 5743237.5d0*x3*y6 - 1013512.5d0*x*y6
      case (7)
        plm = -17229712.5d0*x2*y7 + 1013512.5d0*y7
      case (8)
        plm = 34459425.0d0*x*y8
      case (9)
        plm = -34459425.0d0*y9
      end select ! m
    end select ! l
    return
  end function get_plm 

end module multipole_expansion

!===========================================================================================
SUBROUTINE getvofr(me_r, ps_r, n_me, n_ps, hcub, rhops, potme, guess_state, psgsn, rhops_old, potps_old, nstep)
    !=======================================================================================
    ! Code Version 1.0 (Princeton University, September 2014)
    !=======================================================================================
    ! Given charge density, get the potential using multipole expansion
    ! for boundary region and solving poisson equation for inside box
    ! Adapted from PARSEC by Lingzhu Kong,  http://parsec.ices.utexas.edu/
    !=======================================================================================
    !
    USE kinds,                   ONLY  :  DP
    USE io_global,               ONLY  :  stdout
    USE fft_scalar,              ONLY  :  cfft3d
    USE wannier_base,            ONLY  :  poisson_eps
#ifdef __CUDA
    USE exx_module,              ONLY  : coemicf_d, coeke_d  !MCA/HK : dirty hack for std CG
#endif
    USE exx_module,              ONLY  :  coemicf, coeke
    USE exx_module,              ONLY  :  fbsscale
    USE exx_module,              ONLY  :  nord2
    USE exx_module,              ONLY  :  lmax, lm_mx
    USE exx_module,              ONLY  :  n_exx
    USE funct,                   ONLY  :  get_screening_parameter
    USE mp_global,               ONLY  :  me_image
    USE parallel_include
    USE exx_module,              ONLY  :  pot_ps, rho_ps
    !
    IMPLICIT NONE
    !
    !-----------------------------------------------------------------------
    ! --- call variable ---
    !-----------------------------------------------------------------------
    INTEGER                :: me_r(6)
    INTEGER                :: ps_r(6)
    INTEGER                :: n_me
    INTEGER                :: n_ps
    REAL(DP)               :: hcub
    REAL(DP)               :: rhops(n_ps)
    REAL(DP)               :: potme(n_me)
    INTEGER                :: guess_state
    INTEGER                :: psgsn
    REAL(DP)               :: rhops_old(n_ps, psgsn)
    REAL(DP)               :: potps_old(n_ps, psgsn)
    INTEGER                :: nstep
    !-----------------------------------------------------------------------

    !-----------------------------------------------------------------------
    ! --- local variable ---
    !-----------------------------------------------------------------------
    REAL(DP)                  :: eps
!!    REAL(DP), ALLOCATABLE     :: rho_ps(:)
!!    REAL(DP), ALLOCATABLE     :: pot_ps(:)
    !REAL(DP), ALLOCATABLE     :: rho_me(:)
    !REAL(DP), ALLOCATABLE     :: plr_me(:)
    !REAL(DP), ALLOCATABLE     :: klr_me(:)
    !COMPLEX(DP), ALLOCATABLE  :: rho_ce(:)
    !COMPLEX(DP), ALLOCATABLE  :: plr_ce(:)
    !COMPLEX(DP), ALLOCATABLE  :: klr_ce(:)
    !REAL(DP), ALLOCATABLE     :: GS_rho(:,:)
    !REAL(DP), ALLOCATABLE     :: GS_pot(:,:)
    !REAL(DP), ALLOCATABLE     :: GS_coe(:)
    !REAL(DP), ALLOCATABLE     :: GS_prj(:)
    REAL(DP)                  :: normfactor
    REAL(DP)                  :: rhosum
    !COMPLEX(DP)               :: qlm(0:lmax, 0:lmax)
    COMPLEX(DP), ALLOCATABLE  :: qlm(:, :)
    COMPLEX(DP)               :: qlm_1d(lm_mx)
    INTEGER                   :: gsn, ngsn
    INTEGER                   :: ncb(3)
    INTEGER                   :: itr,jtr
    !-----------------------------------------------------------------------
    REAL(DP)                  :: omega
    !-----------------------------------------------------------------------
    LOGICAL                   :: anti_alising = .FALSE.
    REAL(DP)                  :: sigma = 0.0D0
    !-----------------------------------------------------------------------
#ifdef __CUDA
    !REAL(DP), ALLOCATABLE, DEVICE     :: rhops_d(:)
    !REAL(DP), ALLOCATABLE, DEVICE     :: rho_ps_d(:)
    !REAL(DP), ALLOCATABLE, DEVICE     :: pot_ps_d(:)
    !REAL(DP), ALLOCATABLE, DEVICE     :: potme_d(:) ! don't do implicitly
    attributes(device) :: qlm
    attributes(device) :: potme, rhops
#endif
    

    !-----------------------------------------------------------------------
    ! --- external functions ---
    !-----------------------------------------------------------------------
    REAL(DP), EXTERNAL     :: dnrm2
    !-----------------------------------------------------------------------
    !write(*,*) 'MCA: Starting getvofr'
    !-----------------------------------------------------------------------
    ! --- initialize ---
    !-----------------------------------------------------------------------
    if (.not.allocated(rho_ps)) allocate( rho_ps(n_ps))
    if (.not.allocated(pot_ps)) then
      allocate( pot_ps(n_ps)        ); pot_ps  = 0.0d0
    end if
    !ALLOCATE( rho_me(n_me)        ); rho_me  = 0.0D0
    !ALLOCATE( plr_me(n_me)        ); plr_me  = 0.0D0
    !ALLOCATE( klr_me(n_me)        ); klr_me  = 0.0D0
    !ALLOCATE( rho_ce(n_me)        ); rho_ce  = CMPLX(0.0D0, 0.0D0, kind=DP)
    !ALLOCATE( plr_ce(n_me)        ); plr_ce  = CMPLX(0.0D0, 0.0D0, kind=DP)
    !ALLOCATE( klr_ce(n_me)        ); klr_ce  = CMPLX(0.0D0, 0.0D0, kind=DP)
    !ALLOCATE( GS_rho(n_ps, psgsn) ); GS_rho  = 0.0D0
    !ALLOCATE( GS_pot(n_ps, psgsn) ); GS_pot  = 0.0D0
    !ALLOCATE( GS_coe(psgsn)       ); GS_coe  = 0.0D0
    !ALLOCATE( GS_prj(psgsn)       ); GS_prj  = 0.0D0
    !-----------------------------------------------------------------------
    ncb(1)  = ps_r(4)-ps_r(1)+1
    ncb(2)  = ps_r(5)-ps_r(2)+1
    ncb(3)  = ps_r(6)-ps_r(3)+1
    !-----------------------------------------------------------------------
    gsn     = MIN(guess_state, psgsn)
    !-----------------------------------------------------------------------
    potme   = 0.0D0
    !-----------------------------------------------------------------------
    rho_ps = rhops ! TODO : loop (device)
    !-----------------------------------------------------------------------
    omega = get_screening_parameter()
    !-----------------------------------------------------------------------
!#ifdef __CUDA
!    if (.not.allocated(potme_d)) allocate(potme_d,   mold=potme)
!    potme_d = potme
!#endif
    ALLOCATE(qlm(0:lmax, 0:lmax))

!   !========================================================================
!   ! First, we solve the diffusion equation to get the diffused density
!   !------------------------------------------------------------------------
!   CALL start_clock('getvofr_derk')
!   !------------------------------------------------------------------------
!   rhosum = SUM(rho_ps)
!   !------------------------------------------------------------------------
!   CALL DERK(ncb, 1.0/(omega**2.0), rho_ps, coeke)
!   !------------------------------------------------------------------------
!   IF (me_image .EQ. 0) WRITE(*,"(I3, E15.7, E15.7)") ncb(1), rhosum, SUM(rho_ps)
!   !------------------------------------------------------------------------
!   CALL stop_clock('getvofr_derk')
!   !------------------------------------------------------------------------

    !---------------------------------------------------------------------------
    ! Then, we calculate the potential outside the inner sphere, which will
    ! also give the boundary values for potential inside the sphere.
    !---------------------------------------------------------------------------
    CALL start_clock('getvofr_qlm')
!#ifdef __CUDA
!    !ALLOCATE(rhops_d,  source=rhops)
!    CALL getqlm(ps_r, hcub, rho_ps, qlm)
!#else
!    CALL getqlm(ps_r, hcub, rhops, qlm)
!#endif
    !write(*,*) 'MCA: getqlm'
    CALL getqlm(ps_r, hcub, rho_ps, qlm)
    CALL stop_clock('getvofr_qlm')
    !-----------------------------------------------------------------------
    CALL start_clock('getvofr_bound')
!#ifdef __CUDA
!    CALL exx_boundaryv(me_r, ps_r, potme_d, qlm)
!    potme = potme_d
!#else
    !write(*,*) 'MCA: exx_boundary'
    CALL exx_boundaryv(me_r, ps_r, potme, qlm)
!#endif
    !CALL exx_boundaryv_cuda(me_r, ps_r, potme, qlm_1d)
    CALL stop_clock('getvofr_bound')
    !-----------------------------------------------------------------------
    CALL start_clock('getvofr_geterho')
!#ifdef __CUDA
!    CALL geterho(me_r, ps_r, potme_d, rho_ps)
!#else
    !write(*,*) 'MCA: geterho'
    CALL geterho(me_r, ps_r, potme, rho_ps)
!#endif
    CALL stop_clock('getvofr_geterho')
    !========================================================================
    !TODO : copy out rho_ps before the next step


   ! ! HYK: TODO rethink the GS
   ! !---------------------------------------------------------------------------
   ! CALL start_clock('getvofr_gs')
   ! !-------------------------------------------------------------------------------
   ! !             NEXT, USE GS TO CALCULATE ORTHOGONAL RHO-V PAIRS
   ! !-------------------------------------------------------------------------------
   ! DO itr = 1, gsn
   !   !--------------------------------------
   !   GS_rho(:,itr) = rhops_old(:,itr)
   !   GS_pot(:,itr) = potps_old(:,itr)
   !   !--------------------------------------
   !   DO jtr = 1, itr-1
   !     GS_prj(jtr) = DOT_PRODUCT(GS_rho(:,itr), GS_rho(:,jtr))
   !   END DO
   !   !--------------------------------------
   !   DO jtr = 1, itr-1
   !     GS_rho(:,itr) = GS_rho(:,itr) - GS_prj(jtr) * GS_rho(:,jtr)
   !     GS_pot(:,itr) = GS_pot(:,itr) - GS_prj(jtr) * GS_pot(:,jtr)
   !   END DO
   !   !--------------------------------------
   !   normfactor = 1.0D0/dnrm2(n_ps, GS_rho(:,itr), 1)
   !   !--------------------------------------
   !   GS_rho(:,itr) = GS_rho(:,itr) * normfactor
   !   GS_pot(:,itr) = GS_pot(:,itr) * normfactor
   !   !--------------------------------------
   ! END DO
   ! !---------------------------------------------------------------------------
   ! !                              CALCULATE GUESS
   ! !---------------------------------------------------------------------------
   ! DO itr = 1, gsn
   !   GS_coe(itr) = DOT_PRODUCT(rho_ps, GS_rho(:,itr))
   ! END DO
   ! !---------------------------------------------------------------------------
   ! DO itr = 1, gsn
   !   pot_ps(:) = pot_ps(:) + GS_coe(itr) * GS_pot(:,itr)
   ! END DO
   ! !---------------------------------------------------------------------------
   ! CALL stop_clock('getvofr_gs')
   ! !---------------------------------------------------------------------------

    !---------------------------------------------------------------------------
    ! --- Poisson solver ---
    !---------------------------------------------------------------------------
    CALL start_clock('getvofr_solver')
    !---------------------------------------------------------------------------

!#ifdef __CUDA
!    allocate(pot_ps_d, source=pot_ps)
!#endif
    !write(*,*) 'MCA: cg_solver'
    call cg_solver_stdcg
!#ifdef __CUDA
    !pot_ps = pot_ps_d
!#endif

    !---------------------------------------------------------------------------
    CALL stop_clock('getvofr_solver')
    !---------------------------------------------------------------------------

    ! CALL write_rho_pot(ps_r, rhops, rho_ps, pot_ps)

    !!---------------------------------------------------------------------------
    !ngsn = MIN(gsn+1,psgsn)
    !!---------------------------------------------------------------------------
    !DO itr = ngsn, 2, -1
    !  rhops_old(:,itr) = rhops_old(:,itr-1)
    !  potps_old(:,itr) = potps_old(:,itr-1)
    !END DO
    !!---------------------------------------------------------------------------
    !rhops_old(:,1) = rho_ps(:)
    !potps_old(:,1) = pot_ps(:)
    !!---------------------------------------------------------------------------

    !---------------------------------------------------------------------------
    !write(*,*) 'MCA: ps2me'
    CALL ps2me(me_r, ps_r, potme, pot_ps)
    !write(*,*) 'MCA: after ps2me'
    !---------------------------------------------------------------------------

    ! TODO: potme = potme_d

    !---------------------------------------------------------------------------
    ! up to this point, the Exx potential is computed and stored in me_r next,
    ! we are going to compute the long range potential with convolution
    ! and subtract it from the total Exx potential
    !---------------------------------------------------------------------------
    !IF (omega .NE. 0.0_DP) THEN
    !    !-----------------------------------------------------------------------
    !    IF (anti_alising) sigma = 0.5*hcub**0.3333333333333333D0
    !    !-----------------------------------------------------------------------
    !    CALL ps2me(me_r, ps_r, rho_me, rhops)               ! put rho to ME cube
    !    !-----------------------------------------------------------------------
    !    ! write(*,*) "rho_me: "
    !    ! write(*,'(100000000E20.10)')  rho_me
    !    !-----------------------------------------------------------------------
    !    CALL extend_and_subsample(me_r, rho_me, anti_alising, sigma)
    !    !-----------------------------------------------------------------------
    !    IF (anti_alising) THEN
    !        CALL kernel_lr(me_r, klr_me, (1/omega**2.0D0 - 2.0D0*sigma**2.0D0)**(-0.5D0))
    !    ELSE
    !        CALL kernel_lr(me_r, klr_me, omega)
    !    END IF
    !    !-----------------------------------------------------------------------

    !    !-----------------------------------------------------------------------
    !    !$omp parallel do private(itr)
    !    !-----------------------------------------------------------------------
    !    !DIR$ SIMD
    !    !-----------------------------------------------------------------------
    !    DO itr=1,n_me
    !        rho_ce(itr) = CMPLX(rho_me(itr),0.0D0, kind=DP)
    !        klr_ce(itr) = CMPLX(klr_me(itr),0.0D0, kind=DP)
    !    END DO
    !    !-----------------------------------------------------------------------

    !    CALL start_clock('getvofr_rs')
    !    !-----------------------------------------------------------------------
    !    CALL cfft3d(rho_ce, me_r(4)-me_r(1)+1, me_r(5)-me_r(2)+1, me_r(6)-me_r(3)+1, &
    !                        me_r(4)-me_r(1)+1, me_r(5)-me_r(2)+1, me_r(6)-me_r(3)+1, -1)
    !    CALL cfft3d(klr_ce, me_r(4)-me_r(1)+1, me_r(5)-me_r(2)+1, me_r(6)-me_r(3)+1, &
    !                        me_r(4)-me_r(1)+1, me_r(5)-me_r(2)+1, me_r(6)-me_r(3)+1, -1)
    !    !-----------------------------------------------------------------------
    !    CALL stop_clock('getvofr_rs')

    !    !-----------------------------------------------------------------------
    !    !$omp parallel do private(itr)
    !    !-----------------------------------------------------------------------
    !    !DIR$ SIMD
    !    !-----------------------------------------------------------------------
    !    DO itr=1,n_me
    !        plr_ce(itr) = rho_ce(itr) * klr_ce(itr) * hcub*8.0D0 * n_me
    !    END DO
    !    !-----------------------------------------------------------------------
    !    
    !    !-----------------------------------------------------------------------
    !    CALL cfft3d(plr_ce, me_r(4)-me_r(1)+1, me_r(5)-me_r(2)+1, me_r(6)-me_r(3)+1, &
    !                        me_r(4)-me_r(1)+1, me_r(5)-me_r(2)+1, me_r(6)-me_r(3)+1, 1)
    !    !-----------------------------------------------------------------------

    !    !-----------------------------------------------------------------------
    !    !$omp parallel do private(itr)
    !    !-----------------------------------------------------------------------
    !    !DIR$ SIMD
    !    !-----------------------------------------------------------------------
    !    DO itr=1,n_me
    !        plr_me(itr) = DBLE(plr_ce(itr))
    !    END DO
    !    !-----------------------------------------------------------------------
    !    
    !    !-----------------------------------------------------------------------
    !    CALL shrink_and_interpolate(me_r, plr_me)
    !    !-----------------------------------------------------------------------

    !    !-----------------------------------------------------------------------
    !    ! write(*,*) "plr_me: "
    !    ! write(*,'(100000000E20.10)')  plr_me
    !    !-----------------------------------------------------------------------

    !    !-----------------------------------------------------------------------
    !    !$omp parallel do private(itr)
    !    !-----------------------------------------------------------------------
    !    !DIR$ SIMD
    !    !-----------------------------------------------------------------------
    !    DO itr=1,n_me
    !        potme(itr) = potme(itr) - plr_me(itr)
    !    END DO
    !    !-----------------------------------------------------------------------
    !END IF
    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------

    !---------------------------------------------------------------------------
    !DEALLOCATE( rho_ps )
    !DEALLOCATE( pot_ps )
    !DEALLOCATE( rho_me )
    !DEALLOCATE( plr_me )
    !DEALLOCATE( klr_me )
    !DEALLOCATE( rho_ce )
    !DEALLOCATE( plr_ce )
    !DEALLOCATE( klr_ce )
    !DEALLOCATE( GS_rho )
    !DEALLOCATE( GS_pot )
    !DEALLOCATE( GS_coe )
    !DEALLOCATE( GS_prj )
    DEALLOCATE( qlm    )
!#ifdef __CUDA
!    !deallocate(rhops_d)
!    !deallocate(rho_ps_d)
!    !deallocate(potme_d)
!    !deallocate(pot_ps_d)
!#endif
    !---------------------------------------------------------------------------

    !---------------------------------------------------------------------------
    RETURN
    !---------------------------------------------------------------------------
  contains

    !subroutine  cg_solver_cpu_pcg()
    !  implicit none
    !  CALL CGMIC(nstep, ncb, poisson_eps, fbsscale, coemicf, coeke, rho_ps, pot_ps)
    !  return
    !end subroutine cg_solver_cpu_pcg

    subroutine  cg_solver_stdcg()
      implicit none
#ifdef __CUDA
      !ALLOCATE(rho_ps_d,  source=rho_ps)
      !ALLOCATE(pot_ps_d,  source=pot_ps)
      CALL CGMIC_STDCG(nstep, ncb, poisson_eps, fbsscale, coemicf_d, coeke_d, rho_ps, pot_ps)
      !pot_ps = pot_ps_d
      !DEALLOCATE(rho_ps_d)
      !DEALLOCATE(pot_ps_d)
#else
      CALL CGMIC_STDCG(nstep, ncb, poisson_eps, fbsscale, coemicf, coeke, rho_ps, pot_ps)
#endif
      !
      return
    end subroutine cg_solver_stdcg

    !subroutine  cg_solver_gpu_stdcg()
    !  implicit none
    !  CALL CGMIC_CUDA(nstep, ncb, poisson_eps, fbsscale, coemicf, coeke, rho_ps, pot_ps)
    !  return
    !end subroutine cg_solver_gpu_stdcg

!    subroutine  cg_solver_gpu()
!!			use laplacian_natan_kronik, only : eps10
!!			use laplacian_natan_kronik, only : output_prefactor
!!			use laplacian_natan_kronik, only : output_stencil_coefficients
!!			use laplacian_natan_kronik, only : nord2
!!			use laplacian_natan_kronik, only : eta
!!			use laplacian_natan_kronik, only : drct_aa, drct_bb, drct_cc
!!			use laplacian_natan_kronik, only : drct_ab, drct_ac, drct_bc
!      implicit none
!      integer, device :: ncb_d(3)
!      real(8), device :: coeke_d(-nord2:nord2, 3,3)
!      real(8), device :: coemicf_d(-nord2:nord2, 3,3)
!      real(8), device, allocatable :: rho_ps_d(:,:,:)
!      real(8), device, allocatable :: pot_ps_d(:,:,:)
!			!real(8), device :: prefac_d(drct_aa:drct_bc)   
!			!real(8), device :: stencil_coe_d(-nord2:nord2,drct_aa:drct_bc)
!      !
!      !real(8) :: prefac(drct_aa:drct_bc)
!      !real(8) :: stencil_coe(-nord2:nord2,drct_aa:drct_bc)
!      type(dim3) :: grid, tBlock
!      integer :: istat
!      !
!      !call output_prefactor(prefac)
!			!call output_stencil_coefficients(stencil_coe)
!      !
!      ! alloc GPU vars
!      allocate(rho_ps_d(ncb(1),ncb(2),ncb(3)))
!      allocate(pot_ps_d(ncb(1),ncb(2),ncb(3)))
!      ! copy to GPU
!      ncb_d = ncb
!      !prefac_d = prefac
!      !stencil_coe_d = stencil_coe
!      coeke_d = coeke
!      coemicf_d = coemicf
!      istat = cudaMemcpy(rho_ps_d,rho_ps, product(ncb)) !more efficient memcp
!      tBlock = dim3(8,8,8) ! TODO: arch specific settins...
!      grid = dim3(ceiling(real(ncb(1))/tBlock%x), &
!        &         ceiling(real(ncb(2))/tBlock%y), &
!        &         ceiling(real(ncb(3))/tBlock%z)) 
!      ! GPU kernel
!      call cgmic_cuda<<<grid, tBlock>>>(nstep, ncb_d, poisson_eps, fbsscale, coemicf_d, coeke_d, rho_ps_d, pot_ps_d)
!                                       ! out , in ,    in (value), in(value), 
!      istat = cudaGetLastError() 
!      istat = cudaMemcpy(pot_ps, pot_ps_d, product(ncb)) !more efficient memcp
!      !
!      deallocate(rho_ps_d)
!      deallocate(pot_ps_d)
!      return
!    end subroutine cg_solver_gpu

END SUBROUTINE getvofr
!===============================================================================

!===============================================================================
SUBROUTINE getqlm(ps_r, hcub, rho, qlm)
    !
    USE kinds,            ONLY  :  DP
    USE exx_module,       ONLY  :  lpole=>lmax
    USE exx_module,       ONLY  :  clm
#ifdef __CUDA
    USE exx_module,       ONLY  :  me_cs => me_cs_d
    USE exx_module,       ONLY  :  me_rc => me_rc_d
    USE exx_module,       ONLY  :  me_ri => me_ri_d
    USE exx_module,       ONLY  :  me_rs => me_rs_d
#else
    USE exx_module,       ONLY  :  me_cs
    USE exx_module,       ONLY  :  me_rc
    USE exx_module,       ONLY  :  me_ri
    USE exx_module,       ONLY  :  me_rs
#endif
    USE multipole_expansion, ONLY  :  get_plm
    !
    IMPLICIT NONE
    !------------------------------------------------------------------------------
    ! --- pass in variable ---
    !------------------------------------------------------------------------------
    INTEGER     :: ps_r(6)
    REAL(DP)    :: hcub
    REAL(DP)    :: rho(ps_r(1):ps_r(4),ps_r(2):ps_r(5),ps_r(3):ps_r(6))
    COMPLEX(DP) :: qlm(0:lpole, 0:lpole)
    COMPLEX(DP) :: qlm_tmp
    real(dp) :: coef_l
#ifdef __CUDA
    attributes(device) :: rho, qlm
#endif
    !------------------------------------------------------------------------------
    !
    !------------------------------------------------------------------------------
    REAL(DP)    :: x, y, z, costheta, sintheta, sqrxy
    INTEGER     :: i,j,k,l,m
    !------------------------------------------------------------------------------
    !
    !------------------------------------------------------------------------------
    !------------------------------------------------------------------------------
    DO l = 0, lpole
      DO m = 0, l
        if (m.eq.0) then
          coef_l = clm(l,0)*hcub
        else
          coef_l = 2.d0*clm(l,m)*hcub
        end if ! m.eq.0
        qlm_tmp = (0.d0, 0.d0)
#ifdef __CUDA
        !$cuf kernel do (3)
#else
        !$omp parallel do collapse(3) private(i,j,k) reduction(+:qlm_tmp)
#endif
        DO k=ps_r(3),ps_r(6)
          DO j=ps_r(2),ps_r(5)
            DO i=ps_r(1),ps_r(4)
              !----------------------------
              x = me_cs(1,i,j,k)
              y = me_cs(2,i,j,k)
              z = me_cs(3,i,j,k)
              !------------------------------------------------------------------------------
              sqrxy = DSQRT(x*x+y*y)
              !------------------------------------------------------------------------------
              costheta =     z*me_ri(1,i,j,k)
              sintheta = sqrxy*me_ri(1,i,j,k)
              !------------------------------------------------------------------------------
              qlm_tmp = qlm_tmp + rho(i,j,k)*me_rs(l,i,j,k)*get_plm(costheta,sintheta,l,m)*me_rc(m,i,j,k)*coef_l
            END DO ! i
          END DO ! j
        END DO ! k
#ifndef __CUDA
        !$omp end parallel do
#endif
        qlm(l,m) = qlm_tmp
      END DO ! m
    END DO ! l

    !!------------------------------------------------------------------------------
    !!cuf kernel do (1)
    !DO l = 0, lpole
    !  qlm(l,0) = qlm(l,0)*clm(l,0)*hcub
    !END DO
    !!------------------------------------------------------------------------------
    !!
    !!------------------------------------------------------------------------------
    !DO l = 1, lpole
    !  DO m = 1, l, 1
    !    qlm(l,m) = qlm(l,m)*clm(l,m)*hcub*2.0_DP
    !  END DO
    !END DO
    !!------------------------------------------------------------------------------
    !
    !---------------------------------------------------------------------------
    RETURN
    !---------------------------------------------------------------------------
END SUBROUTINE getqlm
!===============================================================================

!subroutine  compute_multipole_moments(ps_r, hcub, rho, qlm_1d)
!  use kinds,            only : dp
!  use exx_module,       only : lmax, lm_mx
!  use exx_module,       only : clm
!  use exx_module,       only : me_cs, me_rs, me_ri
!  implicit none
!  !
!  integer,  intent(in) :: ps_r(6)
!  real(dp), intent(in) :: hcub
!  real(dp), intent(in) :: rho(ps_r(1):ps_r(4),ps_r(2):ps_r(5),ps_r(3):ps_r(6))
!  !complex(dp), intent(out) :: qlm(0:lmax, 0:lmax) ! HK: todo: should be real in this case
!  complex(dp), intent(out) :: qlm_1d(lm_mx)
!  integer     :: i,j,k,l,m
!  real(dp)    :: plm_1d(lm_mx)
!  real(dp)    :: r, ri, r_l, costheta, sintheta
!  complex(dp) :: exp_niphi !  exp (-i*phi)
!  complex(dp) :: exp_iphi  !  exp (i*phi)
!  complex(dp) :: exp_iphi_m ! exp (i*m*phi)
!  real(dp)    :: buffer(ps_r(1):ps_r(4),ps_r(2):ps_r(5),ps_r(3):ps_r(6),lm_mx)
!  !
!  qlm_1d(:,:) = (0.0d0, 0.0d0)
!  do k=ps_r(3),ps_r(6)
!    do j=ps_r(2),ps_r(5)
!      do i=ps_r(1),ps_r(4)
!        call compute_polar_coordinates(r, ri, costheta, sintheta, exp_niphi)
!        exp_iphi = conjg(exp_niphi)
!        !call setplm(costheta,sintheta,lmax,plm)
!        call setplm_1d(costheta,sintheta,lmax,plm_1d)
!        r_l = 1.d0
!        do l=0, lmax
!          if (l.gt.0) r_l = r_l * ri
!          exp_iphi_m = (1.d0, 0.d0)
!          do m=0, l
!            exp_iphi_m = exp_iphi_m * exp_iphi
!            buffer(i,j,k,lm_1d(l,m)) = rho(i,j,k)*r_l*plm_1d(lm_1d(l,m))*dble(exp_iphi_m)
!          end do
!        end do
!      end do
!    end do
!  end do
!  !
!  do l=0, lmax
!    do m=0, l
!      tmp = 0.d0
!      !$cuf kernel do <<<*,*>>>
!      do k=ps_r(3),ps_r(6)
!        do j=ps_r(2),ps_r(5)
!          do i=ps_r(1),ps_r(4)
!            tmp = tmp + buffer(i,j,k,lm_1d(l,m))
!          end do
!        end do
!      end do
!      if (m.eq.0) then
!        qlm_1d(lm_1d(l,m)) = tmp*clm(l,m)*hcub
!      else
!        qlm_1d(lm_1d(l,m)) = tmp*clm(l,m)*hcub*2.0_dp
!      end if
!    end do
!  end do
!  !
!  return
!end subroutine compute_multipole_moments


!!===============================================================================
!SUBROUTINE qlmc_r(qlm,rho_r,i,j,k,lpole)
!  ! HK: does not seem to be used ... TODO : elliminate
!    !
!    USE kinds,            ONLY  :  DP
!    USE exx_module,       ONLY  :  me_cs, me_rs, me_ri
!    !
!    IMPLICIT NONE
!    !
!    !------------------------------------------------------------------------------
!    ! --- pass in variable ---
!    !------------------------------------------------------------------------------
!    COMPLEX(DP) :: qlm(0:lpole, 0:lpole)
!    REAL(DP)    :: rho_r
!    INTEGER     :: i,j,k
!    INTEGER     :: lpole
!    !------------------------------------------------------------------------------
!    !
!    !------------------------------------------------------------------------------
!    ! --- local variable ---
!    !------------------------------------------------------------------------------
!    INTEGER     :: l,m
!    REAL(DP)    :: x,y,z
!    REAL(DP)    :: sqrxy
!    REAL(DP)    :: costheta,sintheta
!    !------------------------------------------------------------------------------
!    REAL(DP)    :: eps=1.0E-10            ! threshold
!    REAL(DP)    :: plm(0:lpole, 0:lpole)  ! temporary storage of the associated Legendre polynom
!    COMPLEX(DP) :: cxy(1:lpole)           ! coefficient array: e^{i m \phi_j} = cos(phi_j) + i*sin(phi_j)
!    !------------------------------------------------------------------------------
!    !
!
!    !------------------------------------------------------------------------------
!    x = me_cs(1,i,j,k)
!    y = me_cs(2,i,j,k)
!    z = me_cs(3,i,j,k)
!    !------------------------------------------------------------------------------
!    ! WRITE(*,"(I4,I4,I4,F15.11,F15.11,F15.11,F15.11,F15.11,F15.11)") i,j,k,x,y,z,me_rs(1,i,j,k),me_ri(1,i,j,k),rho_r
!    !------------------------------------------------------------------------------
!    sqrxy = DSQRT(x*x+y*y)
!    !------------------------------------------------------------------------------
!    costheta =     z*me_ri(1,i,j,k)
!    sintheta = sqrxy*me_ri(1,i,j,k)
!    !------------------------------------------------------------------------------
!    CALL setplm(costheta,sintheta,lpole,plm)
!    !------------------------------------------------------------------------------
!    qlm(0,0) = qlm(0,0)+rho_r
!    !------------------------------------------------------------------------------
!    DO l = 1,lpole
!      qlm(l,0) = qlm(l,0)+rho_r*me_rs(l,i,j,k)*plm(l,0)                ! qlm(l,m=0)
!    END DO
!    !------------------------------------------------------------------------------
!    IF (sqrxy .GT. eps) THEN
!      !----------------------------------------------------------------------------
!      cxy(1) = CMPLX(x,y,kind=DP)/sqrxy
!      !----------------------------------------------------------------------------
!      DO m = 2,lpole
!        cxy(m) = cxy(m-1)*cxy(1)                          ! cxy(m) = exp(i*m*phi_j)
!      END DO                         
!      !----------------------------------------------------------------------------
!      DO l = 1,lpole
!        DO m = 1,l
!          qlm(l,m) = qlm(l,m) + rho_r*me_rs(l,i,j,k)*plm(l,m)*cxy(m)
!        END DO
!      END DO
!      !----------------------------------------------------------------------------
!    END IF
!    !------------------------------------------------------------------------------
!    !
!    !---------------------------------------------------------------------------
!    RETURN
!    !---------------------------------------------------------------------------
!END SUBROUTINE qlmc_r
!!===============================================================================


!===============================================================================
SUBROUTINE qlmc_rr(qlm,rho_r,i,j,k,lpole)
    !
    USE kinds,            ONLY  :  DP
    USE exx_module,       ONLY  :  me_cs, me_rs, me_ri, me_rc
    !
    IMPLICIT NONE
    !
    !------------------------------------------------------------------------------
    ! --- pass in variable ---
    !------------------------------------------------------------------------------
    COMPLEX(DP) :: qlm(0:lpole, 0:lpole)
    REAL(DP)    :: rho_r
    INTEGER     :: i,j,k
    INTEGER     :: lpole
    !------------------------------------------------------------------------------
    !
    !------------------------------------------------------------------------------
    ! --- local variable ---
    !------------------------------------------------------------------------------
    INTEGER     :: l,m
    REAL(DP)    :: x,y,z
    REAL(DP)    :: sqrxy
    REAL(DP)    :: costheta,sintheta
    !------------------------------------------------------------------------------
    REAL(DP)    :: eps=1.0E-10            ! threshold
    REAL(DP)    :: plm(0:lpole, 0:lpole)  ! temporary storage of the associated Legendre polynom
    COMPLEX(DP) :: cxy(1:lpole)           ! coefficient array: e^{i m \phi_j} = cos(phi_j) + i*sin(phi_j)
    !------------------------------------------------------------------------------

    !------------------------------------------------------------------------------
    x = me_cs(1,i,j,k)
    y = me_cs(2,i,j,k)
    z = me_cs(3,i,j,k)
    !------------------------------------------------------------------------------
    ! WRITE(*,"(I4,I4,I4,F15.11,F15.11,F15.11,F15.11,F15.11,F15.11)") i,j,k,x,y,z,me_rs(1,i,j,k),me_ri(1,i,j,k),rho_r
    !------------------------------------------------------------------------------
    sqrxy = DSQRT(x*x+y*y)
    !------------------------------------------------------------------------------
    costheta =     z*me_ri(1,i,j,k)
    sintheta = sqrxy*me_ri(1,i,j,k)
    !------------------------------------------------------------------------------
    CALL setplm(costheta,sintheta,lpole,plm)
    !------------------------------------------------------------------------------
    DO l=0,lpole
      DO m=0,l
        qlm(l,m) = qlm(l,m) + rho_r*me_rs(l,i,j,k)*plm(l,m)*me_rc(m,i,j,k)
      END DO
    END DO
    !------------------------------------------------------------------------------
     
    !---------------------------------------------------------------------------
    RETURN
    !---------------------------------------------------------------------------
END SUBROUTINE qlmc_rr
!===============================================================================


!===============================================================================
SUBROUTINE exx_boundaryv(me_r, ps_r, potme, qlm)
  !HK: TODO-GPU: swap l loop outside and use grid loop for Cuda acceleration
    !
    USE kinds,            ONLY  :  DP
    USE exx_module,       ONLY  :  lpole=>lmax 
    USE exx_module,       ONLY  :  clm
#ifdef __CUDA
    USE exx_module,       ONLY  :  me_cs => me_cs_d
    USE exx_module,       ONLY  :  me_rc => me_rc_d
    USE exx_module,       ONLY  :  me_ri => me_ri_d
    USE exx_module,       ONLY  :  me_rs => me_rs_d
#else
    USE exx_module,       ONLY  :  me_cs
    USE exx_module,       ONLY  :  me_rc
    USE exx_module,       ONLY  :  me_ri
    USE exx_module,       ONLY  :  me_rs
#endif
    USE multipole_expansion, ONLY  : get_plm
    !
    IMPLICIT  NONE
    !--------------------------------------------------------------------
    ! pass in variables
    !--------------------------------------------------------------------
    INTEGER       :: me_r(6), ps_r(6)
    REAL(DP)      :: potme(me_r(1):me_r(4),me_r(2):me_r(5),me_r(3):me_r(6))
    COMPLEX(DP)   :: qlm(0:lpole, 0:lpole)
#ifdef __CUDA
    attributes(device) :: potme, qlm
#endif
    !--------------------------------------------------------------------
    !
    !--------------------------------------------------------------------
    ! local variables
    !--------------------------------------------------------------------
    INTEGER       :: i,j,k
    !------------------------------------------------------------------------------
    INTEGER       :: l, m
    INTEGER       :: ps_r1, ps_r2, ps_r3, ps_r4, ps_r5, ps_r6
    complex(dp)   :: cpot_r
    REAL(DP)      :: costheta, sintheta, x, y, z, sqrxy
    !------------------------------------------------------------------------------

    ! WRITE(*,*) "multipole expansion pot"
    ! make ps_r as integer to improve cuf kernel perf
    ps_r1 = ps_r(1); ps_r2 = ps_r(2); ps_r3 = ps_r(3)
    ps_r4 = ps_r(4); ps_r5 = ps_r(5); ps_r6 = ps_r(6)

    !------------------------------------------------------------------------------
    ! JJ: to continue
    !------------------------------------------------------------------------------
#ifdef __CUDA
    !$cuf kernel do (3)
#else
    !$omp parallel do private(i,j,k) schedule(guided)
#endif
    DO k=me_r(3),me_r(6)
      DO j=me_r(2),me_r(5)
        DO i=me_r(1),me_r(4)
          !------------------------------------------------------------------------
          IF(.NOT.( (i.GE.ps_r1).AND.(i.LE.ps_r4).AND. &
                    (j.GE.ps_r2).AND.(j.LE.ps_r5).AND. &
                    (k.GE.ps_r3).AND.(k.LE.ps_r6) )) THEN
                  cpot_r = (0.0_DP,0.0_DP)
                  !------------------------------------------------------------------------------
                  x = me_cs(1,i,j,k) ! HK: TODO: dev var
                  y = me_cs(2,i,j,k)
                  z = me_cs(3,i,j,k)
                  !------------------------------------------------------------------------------
                  sqrxy = DSQRT(x*x+y*y)
                  !------------------------------------------------------------------------------
                  costheta =     z*me_ri(1,i,j,k) ! HK: TODO: GPU method should compute me_ri directly as 1/r
                  sintheta = sqrxy*me_ri(1,i,j,k)
                  !------------------------------------------------------------------------------
                  DO l=0,lpole
                    DO m=0,l
                      cpot_r = cpot_r + qlm(l,m)*me_ri(l+1,i,j,k)*get_plm(costheta,sintheta,l,m)*DCONJG(me_rc(m,i,j,k))
                      ! TODO: get_plm device function (need to be in mod)
                    END DO
                  END DO
                  potme(i,j,k) = dble(cpot_r)
          END IF
          !------------------------------------------------------------------------
        END DO
      END DO
    END DO
#ifndef __CUDA
    !$omp end parallel do
#endif
    !------------------------------------------------------------------------------
    !------------------------------------------------------------------------------
    !
    !------------------------------------------------------------------------------
    RETURN
    !---------------------------------------------------------------------------
END SUBROUTINE exx_boundaryv
!===========================================================================

!subroutine  exx_boundaryv_cuda(me_r, ps_r, potme, qlm_1d)
!  ! rename to compute_farfield_potential
!  ! HK: do not store large array, but rather compute via direct algorithm
!  !
!  use kinds,            only  :  dp
!  use exx_module,       only  :  lmax, lm_mx
!  use exx_module,       only  :  clm
!  use exx_module,       only  :  me_cs, me_rs, me_ri
!  implicit none
!  !
!  integer,     intent(in)    :: me_r(6), ps_r(6)
!  real(dp),    intent(inout) :: potme(me_r(1):me_r(4),me_r(2):me_r(5),me_r(3):me_r(6))
!  complex(dp), intent(in)    :: qlm_1d(lm_mx)
!  real(dp)                   :: plm_1d(lm_mx)  ! temporary storage of the associated legendre polynomials
!  integer       :: l,m,i,j,k
!  real(dp)      :: r, ri, ri_lplus1, costheta, sintheta
!  complex(dp)   :: exp_niphi !  exp (-i*phi)
!  complex(dp)   :: exp_niphi_m ! exp (-i*m*phi)
!  do k = me_r(3), me_r(6)
!    do j = me_r(2), me_r(5)
!      do i = me_r(1), me_r(4)
!        if (is_grid_outside_pe(i,j,k,ps_r)) then
!          potme(i,j,k) = 0.d0
!          !call compute_polar_coordinates ! yields ri, costheta, sintheta, exp_niphi
!          call compute_polar_coordinates(r,ri, costheta, sintheta, exp_niphi)
!          call setplm_1d(costheta, sintheta, plm_1d) ! HK: TODO: only needed for variable cell simulations use condition with thdyn
!          !                                          !   : also may use recursion relation to reduce computational cost
!          ri_lplus1 = 1.d0
!          do l = 0, lmax
!            ri_lplus1 = ri_lplus1 * ri
!            exp_niphi_m = (1.d0, 0.d0)
!            do m = 0, l
!              exp_niphi_m = exp_niphi_m * exp_niphi
!              potme(i,j,k) = potme(i,j,k) &
!                &            + qlm_1d(lm_1d(l,m))*ri_lplus1*plm_1d(lm_1d(l,m))*dble(exp_niphi_m)
!            end do ! m
!          end do ! l
!        end if ! is_grid_outside_pe(i,j,k,ps_r)
!      end do ! i
!    end do ! j
!  end do ! k
!  return
!end subroutine exx_boundaryv_cuda
!
!logical function is_grid_outside_pe(i,j,k,pe_r)
!  implicit none
!  integer, intent(in) :: i, j, k, ps_r(6)
!  is_grid_outside_pe = .not.((i.ge.ps_r(1)).and.(i.le.ps_r(4)).and. &
!    & (j.ge.ps_r(2)).and.(j.le.ps_r(5)).and.(k.ge.ps_r(3)).and.(k.le.ps_r(6)))
!end function is_grid_outside_pe
!
!subroutine  compute_polar_coordinates(r, ri, costheta, sintheta, exp_niphi)
!  ! HK: TODO: decide on using direct/indirect algorithm (see also compute_polar_coordinates_original)
!  use cell_base,   only : h
!  use fft_base,    only : dfftp
!  implicit none
!  real(dp) :: x, y, z, sqrxy, s(3) ! for GPU, used as input argument to prevent allocation of memory
!  real(8), intent(out) :: r, ri, costheta, sintheta, exp_niphi
!  !
!  s(1) = (DBLE(i)/DBLE(dfftp%nr1)) - DBLE(INT(dfftp%nr1/2))/DBLE(dfftp%nr1)
!  s(2) = (DBLE(j)/DBLE(dfftp%nr2)) - DBLE(INT(dfftp%nr2/2))/DBLE(dfftp%nr2)
!  s(3) = (DBLE(k)/DBLE(dfftp%nr3)) - DBLE(INT(dfftp%nr3/2))/DBLE(dfftp%nr3)
!  !
!  x = h(1,1)*s(1)+h(1,2)*s(2)+h(1,3)*s(3)   !r_i = h s_i
!  y = h(2,1)*s(1)+h(2,2)*s(2)+h(2,3)*s(3)   !r_i = h s_i
!  z = h(3,1)*s(1)+h(3,2)*s(2)+h(3,3)*s(3)   !r_i = h s_i
!  !
!  sqrxy = x*x+y*y
!  r = sqrt(sqrxy+z*z)
!  sqrxy = sqrt(sqrxy)
!  ri = 1.d0/r
!  ! 
!  costheta = z     * ri
!  sintheta = sqrxy * ri
!  !
!  if (sqrxy.gt.1.d-10) then
!    exp_niphi = cmplx(x,-y)/sqrxy
!  else
!    exp_niphi = (1.d0,0.d0)
!  end if
!  !
!  return
!end subroutine compute_polar_coordinates
!
!subroutine  compute_polar_coordinates_tmp(r, ri, costheta, sintheta, exp_niphi)
!  ! HK: TODO: decide on using direct/indirect algorithm (see also compute_polar_coordinates_original)
!  use cell_base,   only : h
!  use fft_base,    only : dfftp
!  implicit none
!  real(dp) :: x, y, z, sqrxy, s(3) ! for GPU, used as input argument to prevent allocation of memory
!  real(8), intent(out) :: r, ri, costheta, sintheta, exp_niphi
!  !
!  s(1) = (DBLE(i)/DBLE(dfftp%nr1)) - DBLE(INT(dfftp%nr1/2))/DBLE(dfftp%nr1)
!  s(2) = (DBLE(j)/DBLE(dfftp%nr2)) - DBLE(INT(dfftp%nr2/2))/DBLE(dfftp%nr2)
!  s(3) = (DBLE(k)/DBLE(dfftp%nr3)) - DBLE(INT(dfftp%nr3/2))/DBLE(dfftp%nr3)
!  !
!  s = matmul(h,s)
!  !
!  r = sqrt(s(1)*s(1)+s(2)*s(2)+s(3)*s(3))
!  !sqrxy = s(1)*s(1)+s(2)*s(2)
!  !r = sqrt(sqrxy+z*z)
!  !sqrxy = sqrt(sqrxy)
!  !ri = 1.d0/r
!  ! 
!  costheta = s(3) / r
!  sintheta = sqrt(s(1)*s(1)+s(2)*s(2)) / r
!  !
!  if (sqrxy.gt.1.d-10) then
!    exp_niphi = cmplx(s(1),-s(2))/sqrt(s(1)*s(1)+s(2)*s(2))
!  else
!    exp_niphi = (1.d0,0.d0)
!  end if
!  !
!  return
!end subroutine compute_polar_coordinates_tmp

!!===========================================================================
!FUNCTION qlme_r(qlm,i,j,k,lpole) RESULT(pot_r)
!    !
!    USE kinds,            ONLY  :  DP
!    USE exx_module,       ONLY  :  me_cs, me_rs, me_ri
!    !
!    IMPLICIT NONE
!    !
!    !------------------------------------------------------------------------------
!    ! --- pass in variable ---
!    !------------------------------------------------------------------------------
!    COMPLEX(DP) :: qlm(0:lpole, 0:lpole) ! HK: TODO: flattern to 1D using starting point s(l)=lmax^2+1 for each l subspace
!    INTEGER     :: i,j,k
!    INTEGER     :: lpole
!    !------------------------------------------------------------------------------
!    REAL(DP)    :: pot_r
!    !------------------------------------------------------------------------------
!    !
!    !------------------------------------------------------------------------------
!    ! --- local variable ---
!    !------------------------------------------------------------------------------
!    INTEGER     :: l,m
!    REAL(DP)    :: rho_r
!    REAL(DP)    :: x,y,z
!    REAL(DP)    :: sqrxy
!    REAL(DP)    :: costheta,sintheta
!    !------------------------------------------------------------------------------
!    REAL(DP)    :: eps=1.0E-10            ! threshold
!    REAL(DP)    :: plm(0:lpole, 0:lpole)  ! temporary storage of the associated Legendre polynom
!    !                                     ! HK: TODO: flattern to 1D using starting point s(l)=lmax^2+1 for each l subspace
!    COMPLEX(DP) :: cxy(1:lpole)           ! coefficient array: e^{i m \phi_j} = cos(phi_j) + i*sin(phi_j)
!    COMPLEX(DP) :: cpot_r
!    !------------------------------------------------------------------------------
!
!    !------------------------------------------------------------------------------
!    cpot_r = 0.0_DP
!    !------------------------------------------------------------------------------
!    x = me_cs(1,i,j,k)
!    y = me_cs(2,i,j,k)
!    z = me_cs(3,i,j,k)
!    !------------------------------------------------------------------------------
!    sqrxy = DSQRT(x*x+y*y)
!    !------------------------------------------------------------------------------
!    costheta =     z*me_ri(1,i,j,k)
!    sintheta = sqrxy*me_ri(1,i,j,k)
!    !------------------------------------------------------------------------------
!    CALL setplm(costheta,sintheta,lpole,plm)
!    !------------------------------------------------------------------------------
!    DO l = 0,lpole
!      cpot_r = cpot_r + qlm(l,0)*me_ri(l+1,i,j,k)*plm(l,0)              ! m=0 terms
!    END DO
!    !------------------------------------------------------------------------------
!    IF (sqrxy .GT. eps) THEN
!      !----------------------------------------------------------------------------
!      cxy(1) = CMPLX(x,-y,kind=DP)/sqrxy
!      !----------------------------------------------------------------------------
!      DO m = 2,lpole
!        cxy(m) = cxy(m-1)*cxy(1)                         ! cxy(m) = exp(-i*m*phi_j)
!      END DO                         
!      !----------------------------------------------------------------------------
!      DO l = 1,lpole
!        DO m = 1,l
!          cpot_r = cpot_r + qlm(l,m)*me_ri(l+1,i,j,k)*plm(l,m)*cxy(m)
!        END DO
!      END DO
!      !----------------------------------------------------------------------------
!    END IF
!    !------------------------------------------------------------------------------
!    !
!    !------------------------------------------------------------------------------
!    pot_r = REAL(cpot_r)
!    !------------------------------------------------------------------------------
!    !
!    ! WRITE(*,"(I3,I3,I3,F14.7,F14.7,F14.7,F14.7,F14.7,ES17.7)") i,j,k,x,y,z,me_rs(1,i,j,k),me_ri(1,i,j,k),pot_r
!    !
!    !---------------------------------------------------------------------------
!    RETURN
!    !---------------------------------------------------------------------------
!END FUNCTION qlme_r
!!===========================================================================

!===========================================================================
FUNCTION qlme_rr(qlm,i,j,k,lpole) RESULT(pot_r)
    !
    USE kinds,            ONLY  :  DP
    USE exx_module,       ONLY  :  me_cs, me_rs, me_ri, me_rc
    !
    IMPLICIT NONE
    !
    !------------------------------------------------------------------------------
    ! --- pass in variable ---
    !------------------------------------------------------------------------------
    COMPLEX(DP) :: qlm(0:lpole, 0:lpole)
    INTEGER     :: i,j,k
    INTEGER     :: lpole
    !------------------------------------------------------------------------------
    REAL(DP)    :: pot_r
    !------------------------------------------------------------------------------
    !
    !------------------------------------------------------------------------------
    ! --- local variable ---
    !------------------------------------------------------------------------------
    INTEGER     :: l,m
    REAL(DP)    :: rho_r
    REAL(DP)    :: x,y,z
    REAL(DP)    :: sqrxy
    REAL(DP)    :: costheta,sintheta
    !------------------------------------------------------------------------------
    REAL(DP)    :: eps=1.0E-10            ! threshold
    REAL(DP)    :: plm(0:lpole, 0:lpole)  ! temporary storage of the associated Legendre polynom
    COMPLEX(DP) :: cxy(1:lpole)           ! coefficient array: e^{i m \phi_j} = cos(phi_j) + i*sin(phi_j)
    COMPLEX(DP) :: cpot_r
    !------------------------------------------------------------------------------

    !------------------------------------------------------------------------------
    cpot_r = 0.0_DP
    !------------------------------------------------------------------------------
    x = me_cs(1,i,j,k) ! HK: TODO: GPU method should compute x,y,z directly
    y = me_cs(2,i,j,k)
    z = me_cs(3,i,j,k)
    !------------------------------------------------------------------------------
    sqrxy = DSQRT(x*x+y*y)
    !------------------------------------------------------------------------------
    costheta =     z*me_ri(1,i,j,k) ! HK: TODO: GPU method should compute me_ri directly as 1/r
    sintheta = sqrxy*me_ri(1,i,j,k)
    !------------------------------------------------------------------------------
    CALL setplm(costheta,sintheta,lpole,plm) ! HK: TODO: only needed for variable cell simulations use condition with thdyn
    !------------------------------------------------------------------------------
!   DO l=0,lpole
!     cpot_r = cpot_r + qlm(l,0)*me_ri(l+1,i,j,k)*plm(l,m)*CONJG(me_rc(0,i,j,k))
!     write(*,*) cpot_r
!   END DO
!   !------------------------------------------------------------------------------
    ! HK: TODO: consider how to improve the GPU support of this integration
    DO l=0,lpole
      DO m=0,l
        cpot_r = cpot_r + qlm(l,m)*me_ri(l+1,i,j,k)*plm(l,m)*DCONJG(me_rc(m,i,j,k))
      END DO
    END DO
    !------------------------------------------------------------------------------

    !------------------------------------------------------------------------------
    pot_r = DBLE(cpot_r)
    !------------------------------------------------------------------------------
    !
    ! WRITE(*,"(I3,I3,I3,F14.7,F14.7,F14.7,F14.7,F14.7,ES17.7)") i,j,k,x,y,z,me_rs(1,i,j,k),me_ri(1,i,j,k),pot_r
    !
    !---------------------------------------------------------------------------
    RETURN
    !---------------------------------------------------------------------------
END FUNCTION qlme_rr
!===========================================================================


!===========================================================================
SUBROUTINE geterho(me_r, ps_r, potme, rhops)
    !
    USE kinds,            ONLY  :  DP
    USE exx_module,       ONLY  :  nord2
#ifdef __CUDA
    USE exx_module,       ONLY  :  coeke => coeke_d
#else
    USE exx_module,       ONLY  :  coeke
#endif
    USE constants,        ONLY  :  eps12
    !
    IMPLICIT NONE
    !--------------------------------------------------------------------
    ! pass in variables
    !--------------------------------------------------------------------
    INTEGER       :: me_r(6), ps_r(6)
    REAL(DP)      :: potme(me_r(1):me_r(4),me_r(2):me_r(5),me_r(3):me_r(6))
    REAL(DP)      :: rhops(ps_r(1):ps_r(4),ps_r(2):ps_r(5),ps_r(3):ps_r(6))
#ifdef __CUDA
    attributes(device) :: potme, rhops
#endif
    !--------------------------------------------------------------------
    !
    !--------------------------------------------------------------------
    ! local variables
    !--------------------------------------------------------------------
    INTEGER       :: i,j,k,ish
    INTEGER       :: ps_r1, ps_r2, ps_r3, ps_r4, ps_r5, ps_r6
    !------------------------------------------------------------------------------

    ! WRITE(*,*) "wrapped rho"
    !TODO: to GPU : rhops, potme
    ps_r1 = ps_r(1); ps_r2 = ps_r(2); ps_r3 = ps_r(3)
    ps_r4 = ps_r(4); ps_r5 = ps_r(5); ps_r6 = ps_r(6)

    !------------------------------------------------------------------------------
#ifdef __CUDA
    !$cuf kernel do(3)
#else
    !$omp parallel do private(i,j,k,ish) schedule(guided)
#endif
    !------------------------------------------------------------------------------
    DO k=ps_r3,ps_r6
      DO j=ps_r2,ps_r5
        DO i=ps_r1,ps_r4
          !------------------------------------------------------------------------
          IF(.NOT.( (i.GE.ps_r1+nord2).AND.(i.LE.ps_r4-nord2).AND. &
                    (j.GE.ps_r2+nord2).AND.(j.LE.ps_r5-nord2).AND. &
                    (k.GE.ps_r3+nord2).AND.(k.LE.ps_r6-nord2) )) THEN
            !----------------------------------------------------------------------
            DO ish=1,nord2
              rhops(i,j,k) = rhops(i,j,k) - coeke(ish,1,1)*potme(i+ish, j,     k    ) &
                                          - coeke(ish,1,1)*potme(i-ish, j,     k    ) &
                                          - coeke(ish,2,2)*potme(i,     j+ish, k    ) &
                                          - coeke(ish,2,2)*potme(i,     j-ish, k    ) &
                                          - coeke(ish,3,3)*potme(i,     j,     k+ish) &
                                          - coeke(ish,3,3)*potme(i,     j,     k-ish) &
                                          - coeke(ish,1,2)*potme(i+ish, j+ish, k    ) &
                                          + coeke(ish,1,2)*potme(i+ish, j-ish, k    ) &
                                          + coeke(ish,1,2)*potme(i-ish, j+ish, k    ) &
                                          - coeke(ish,1,2)*potme(i-ish, j-ish, k    ) &
                                          - coeke(ish,1,3)*potme(i+ish, j,     k+ish) &
                                          + coeke(ish,1,3)*potme(i+ish, j,     k-ish) &
                                          + coeke(ish,1,3)*potme(i-ish, j,     k+ish) &
                                          - coeke(ish,1,3)*potme(i-ish, j,     k-ish) &
                                          - coeke(ish,2,3)*potme(i,     j+ish, k+ish) &
                                          + coeke(ish,2,3)*potme(i,     j+ish, k-ish) &
                                          + coeke(ish,2,3)*potme(i,     j-ish, k+ish) &
                                          - coeke(ish,2,3)*potme(i,     j-ish, k-ish)
            END DO
            !----------------------------------------------------------------------
          END IF
          !------------------------------------------------------------------------
          ! WRITE(*,"(I4,I4,I4,F15.11)") i,j,k, rhops(i,j,k)
          !------------------------------------------------------------------------
        END DO
      END DO
    END DO
#ifndef __CUDA
    !$omp end parallel do 
#endif
    !-----------------------------------------------------------------------
    RETURN
    !-----------------------------------------------------------------------
END SUBROUTINE geterho
!==================================================================

!==================================================================
SUBROUTINE setplm(x, y, lpole, plm)
    !
    ! this subroutine calculates all of the associated Legendre polynomials up 
    ! to a supplied lpole (maximum lpole is 9). 
    !
    IMPLICIT NONE
    !
    ! INPUT VARIABLES
    ! cos(theta),sin(theta), for a given grid point
    REAL*8 x,y
    ! order of multipole expansion
    INTEGER lpole
    ! OUTPUT VARIABLES
    ! array containing P_{lm}, the associated Legendre polynomials
    REAL*8 plm(0:lpole, 0:lpole)
    ! WORK VARIABLES:
    ! powers of x, y: xn=x^n, yn=y^n
    REAL*8 x2, x3, x4, x5, x6, x7, x8, x9
    REAL*8 y2, y3, y4, y5, y6, y7, y8, y9
    !
    plm(0,0) = 1.0
    IF (lpole .GE. 1) THEN
      plm(1,0) = x
      plm(1,1) = -y
    END IF
    IF (lpole .GE. 2) THEN
      x2 = x*x
      y2 = y*y
      plm(2,0) =  1.5d0*x2 - 0.5d0
      plm(2,1) = -3.0d0*x*y
      plm(2,2) =  3.0d0*y2
    END IF
    IF (lpole .GE. 3) THEN
      x3 = x2*x
      y3 = y2*y
      plm(3,0) =   2.5d0*x3 - 1.5d0*x
      plm(3,1) = (-7.5d0*x2 + 1.5d0)*y
      plm(3,2) =  15.0d0*x*y2
      plm(3,3) = -15.0d0*y3
    END IF
    IF (lpole .GE. 4) THEN
      x4 = x2*x2
      y4 = y2*y2
      plm(4,0) =  4.375d0*x4 - 3.75d0*x2 + 0.375d0
      plm(4,1) = (-17.5d0*x3 + 7.5d0*x)*y
      plm(4,2) = ( 52.5d0*x2 - 7.5d0  )*y2
      plm(4,3) = -105.0d0*x*y3
      plm(4,4) =  105.0d0*y4
    END IF
    IF (lpole .GE. 5) THEN
      x5 = x3*x2
      y5 = y3*y2
      plm(5,0) =  7.875d0*x5 - 8.75d0*x3 + 1.875d0*x
      plm(5,1) = (-39.375d0*x4 + 26.25d0*x2 - 1.875d0)*y
      plm(5,2) = ( 157.5d0*x3 - 52.5d0*x)*y2
      plm(5,3) = (-472.5d0*x2 + 52.5d0)  *y3
      plm(5,4) =  945.0d0*x*y4
      plm(5,5) = -945.0d0*y5
    END IF
    IF (lpole .GE. 6) THEN
      x6 = x3*x3
      y6 = y3*y3
      plm(6,0) = 14.4375d0*x6 - 19.6875d0*x4 + 6.5625d0*x2 - 0.3125d0
      plm(6,1) = (-86.625d0*x5 + 78.75d0*x3 - 13.125d0*x)*y
      plm(6,2) = ( 433.125d0*x4 - 236.25d0*x2 + 13.125d0)*y2
      plm(6,3) = (-1732.5d0*x3 + 472.5d0*x            )*y3
      plm(6,4) = ( 5197.5d0*x2 - 472.5d0              )*y4
      plm(6,5) = -10395.0d0*x*y5
      plm(6,6) =  10395.0d0*y6
    END IF
    IF (lpole .GE. 7) THEN
      x7 = x4*x3
      y7 = y4*y3
      plm(7,0) = 26.8125d0*x7 - 43.3125d0*x5 + 19.6875d0*x3 - 2.1875d0*x
      plm(7,1) = -187.6875d0*x6*y + 216.5625d0*x4*y - 59.0625d0*x2*y + 2.1875d0*y
      plm(7,2) = 1126.125d0*x5*y2 - 866.25d0*x3*y2 + 118.125d0*x*y2
      plm(7,3) = -5630.625d0*x4*y3 + 2598.75d0*x2*y3 - 118.125d0*y3
      plm(7,4) = 22522.5d0*x3*y4 - 5197.5d0*x*y4
      plm(7,5) = -67567.5d0*x2*y5 + 5197.5d0*y5
      plm(7,6) = 135135.0d0*x*y6
      plm(7,7) = -135135.0d0*y7
    END IF
    IF (lpole .GE. 8) THEN
      x8 = x4*x4
      y8 = y4*y4
      plm(8,0) = 50.2734375d0*x8 - 93.84375d0*x6 + 54.140625d0*x4 - 9.84375d0*x2 + 0.2734375d0
      plm(8,1) = -402.1875d0*x7*y + 563.0625d0*x5*y - 216.5625d0*x3*y + 19.6875d0*x*y
      plm(8,2) = 2815.3125d0*x6*y2 - 2815.3125d0*x4*y2 + 649.6875d0*x2*y2 - 19.6875d0*y2
      plm(8,3) = -16891.875d0*x5*y3 + 11261.25d0*x3*y3 - 1299.375d0*x*y3
      plm(8,4) = 84459.375d0*x4*y4 - 33783.75d0*x2*y4 + 1299.375d0*y4
      plm(8,5) = -337837.5d0*x3*y5 + 67567.5d0*x*y5
      plm(8,6) = 1013512.5d0*x2*y6 - 67567.5d0*y6
      plm(8,7) = -2027025.0d0*x*y7
      plm(8,8) = 2027025.0d0*y8
    END IF
    IF (lpole .GE. 9) THEN
      y9 = y5*y4
      x9 = x5*x4
      plm(9,0) = 94.9609375d0*x9 - 201.09375d0*x7 + 140.765625d0*x5 - 36.09375d0*x3 + 2.4609375d0*x
      plm(9,1) = -854.6484375d0*x8*y + 1407.65625d0*x6*y - 703.828125d0*x4*y + 108.28125d0*x2*y - 2.4609375d0*y
      plm(9,2) = 6837.1875d0*x7*y2 - 8445.9375d0*x5*y2 + 2815.3125d0*x3*y2 - 216.5625d0*x*y2
      plm(9,3) = -47860.3125d0*x6*y3 + 42229.6875d0*x4*y3 - 8445.9375d0*x2*y3 + 216.5625d0*y3
      plm(9,4) = 287161.875d0*x5*y4 - 168918.75d0*x3*y4 + 16891.875d0*x*y4
      plm(9,5) = -1435809.375d0*x4*y5 + 506756.25d0*x2*y5 - 16891.875d0*y5
      plm(9,6) = 5743237.5d0*x3*y6 - 1013512.5d0*x*y6
      plm(9,7) = -17229712.5d0*x2*y7 + 1013512.5d0*y7
      plm(9,8) = 34459425.0d0*x*y8
      plm(9,9) = -34459425.0d0*y9
    END IF
    !
    !--------------------------------------------------------------------
    RETURN
    !--------------------------------------------------------------------
END SUBROUTINE setplm 
!==================================================================

!==================================================================
SUBROUTINE write_rho_pot(ps_r, rhops, rho_ps, pot_ps)
    USE kinds,            ONLY  :  DP
    !
    IMPLICIT  NONE
    !--------------------------------------------------------------------
    INTEGER       :: ps_r(6)
    REAL(DP)      :: rhops(ps_r(1):ps_r(4),ps_r(2):ps_r(5),ps_r(3):ps_r(6))
    REAL(DP)      :: rho_ps(ps_r(1):ps_r(4),ps_r(2):ps_r(5),ps_r(3):ps_r(6))
    REAL(DP)      :: pot_ps(ps_r(1):ps_r(4),ps_r(2):ps_r(5),ps_r(3):ps_r(6))
    !--------------------------------------------------------------------
    INTEGER       :: i,j,k
    !--------------------------------------------------------------------
    DO k=ps_r(3),ps_r(6)
      DO j=ps_r(2),ps_r(5)
        DO i=ps_r(1),ps_r(4)
        !----------------------------------------------------------------
        WRITE(*,"(I4,I4,I4,F15.11,F15.11,F15.11)") i,j,k, rhops(i,j,k),rho_ps(i,j,k),pot_ps(i,j,k)
        !----------------------------------------------------------------
        END DO
      END DO
    END DO
    !--------------------------------------------------------------------
    RETURN
    !--------------------------------------------------------------------
END SUBROUTINE write_rho_pot
!==================================================================

!==================================================================
SUBROUTINE ps2me(me_r, ps_r, potme, potps)
    USE kinds,            ONLY  :  DP
    !
    IMPLICIT  NONE
    !--------------------------------------------------------------------
    INTEGER       :: me_r(6), ps_r(6)
    REAL(DP)      :: potme(me_r(1):me_r(4),me_r(2):me_r(5),me_r(3):me_r(6))
    REAL(DP)      :: potps(ps_r(1):ps_r(4),ps_r(2):ps_r(5),ps_r(3):ps_r(6))
#ifdef __CUDA
    attributes(device) :: potps, potme
#endif
    !--------------------------------------------------------------------
    ! local variables
    !--------------------------------------------------------------------
    INTEGER       :: i,j,k
    INTEGER       :: ps_r1, ps_r2, ps_r3, ps_r4, ps_r5, ps_r6
    !------------------------------------------------------------------------------
    !
    ps_r1 = ps_r(1); ps_r2 = ps_r(2); ps_r3 = ps_r(3)
    ps_r4 = ps_r(4); ps_r5 = ps_r(5); ps_r6 = ps_r(6)
    !
#ifdef __CUDA
    !$cuf kernel do (3)
#else
    !$omp parallel do private(i,j,k)
#endif
    DO k=ps_r3,ps_r6
      DO j=ps_r2,ps_r5
        DO i=ps_r1,ps_r4
          potme(i,j,k) = potps(i,j,k)
        END DO
      END DO
    END DO
#ifndef __CUDA
    !$omp end parallel do 
#endif
    !--------------------------------------------------------------------
    !potme(ps_r(1):ps_r(4),ps_r(2):ps_r(5),ps_r(3):ps_r(6)) = potps(:,:,:)
    !--------------------------------------------------------------------
    RETURN
    !--------------------------------------------------------------------
END SUBROUTINE ps2me
!==================================================================


!!==================================================================
!SUBROUTINE extend_and_subsample(me_r, rho_me, anti_alising, sigma)
!    USE kinds,            ONLY  :  DP
!    USE exx_module,       ONLY  :  me_rs
!    USE exx_module,       ONLY  :  nrgr
!    !
!    IMPLICIT  NONE
!    !--------------------------------------------------------------------
!    INTEGER                :: me_r(6)
!    REAL(DP)               :: rho_me(me_r(1):me_r(4),me_r(2):me_r(5),me_r(3):me_r(6))
!    LOGICAL                :: anti_alising
!    REAL(DP)               :: sigma
!    !--------------------------------------------------------------------
!    INTEGER                :: npd(3)
!    INTEGER                :: mc_r(6)
!    REAL(DP), ALLOCATABLE  :: rho_mc(:,:,:)
!    !--------------------------------------------------------------------
!    mc_r(1)=1
!    mc_r(2)=1
!    mc_r(3)=1
!    mc_r(4)=(me_r(4)-me_r(1)+2)/2-1
!    mc_r(5)=(me_r(5)-me_r(2)+2)/2-1
!    mc_r(6)=(me_r(6)-me_r(3)+2)/2-1
!    !--------------------------------------------------------------------
!    npd(1) = ((me_r(4)-me_r(1)+1) - (mc_r(4)-mc_r(1)+1)) / 2
!    npd(2) = ((me_r(5)-me_r(2)+1) - (mc_r(5)-mc_r(2)+1)) / 2
!    npd(3) = ((me_r(6)-me_r(3)+1) - (mc_r(6)-mc_r(3)+1)) / 2
!    !--------------------------------------------------------------------
!    ALLOCATE(rho_mc(mc_r(1):mc_r(4), mc_r(2):mc_r(5), mc_r(3):mc_r(6))); rho_mc=0.0D0
!    !--------------------------------------------------------------------
!    IF (anti_alising) THEN
!        CALL interp_f2c_gaussian(me_r,me_r,rho_me,mc_r,mc_r,rho_mc, me_rs(1,nrgr(1)-1:nrgr(1)+1, &
!                                                                            nrgr(2)-1:nrgr(2)+1, &
!                                                                            nrgr(3)-1:nrgr(3)+1)/sigma)
!    ELSE
!        CALL interp_f2c_linear(me_r,me_r,rho_me,mc_r,mc_r,rho_mc)
!    END IF
!    !--------------------------------------------------------------------
!    rho_me = 0.0D0
!    !--------------------------------------------------------------------
!    rho_me(me_r(1)+npd(1):me_r(4)-npd(1), me_r(2)+npd(2):me_r(5)-npd(2), me_r(3)+npd(3):me_r(6)-npd(3)) = rho_mc(:,:,:)
!    !--------------------------------------------------------------------
!    DEALLOCATE(rho_mc)
!    !--------------------------------------------------------------------
!END SUBROUTINE extend_and_subsample
!!==================================================================


!==================================================================
SUBROUTINE kernel_lr(me_r, klr_me, omega)
    USE kinds,            ONLY  :  DP
    USE exx_module,       ONLY  :  me_rs
    !
    IMPLICIT  NONE
    !--------------------------------------------------------------------
    INTEGER                :: me_r(6)
    REAL(DP)               :: klr_me(me_r(1):me_r(4),me_r(2):me_r(5),me_r(3):me_r(6))
    REAL(DP)               :: omega
    !--------------------------------------------------------------------
    REAL(DP)               :: r
    INTEGER                :: nb(3)
    INTEGER                :: i,j,k
    INTEGER                :: is,js,ks
    !--------------------------------------------------------------------

    nb(1) = (me_r(4)-me_r(1)+1)
    nb(2) = (me_r(5)-me_r(2)+1)
    nb(3) = (me_r(6)-me_r(3)+1)

    !--------------------------------------------------------------------
    !$omp parallel do private(i,j,k,r) schedule(guided)
    !--------------------------------------------------------------------
    DO k=me_r(3),me_r(6)
      DO j=me_r(2),me_r(5)
        DO i=me_r(1),me_r(4)
        !----------------------------------------------------------------
        is = MOD(i-me_r(1)+nb(1)/2, nb(1))+me_r(1)
        js = MOD(j-me_r(2)+nb(2)/2, nb(2))+me_r(2)
        ks = MOD(k-me_r(3)+nb(3)/2, nb(3))+me_r(3)
        !----------------------------------------------------------------
        r = me_rs(1,is,js,ks) * 2.0D0
        !----------------------------------------------------------------
        klr_me(i,j,k) = erf(omega*r)/r
        !----------------------------------------------------------------
        END DO
      END DO
    END DO
    !--------------------------------------------------------------------
    klr_me(me_r(1), me_r(2), me_r(3)) = (2*omega)/(3.1415926535897932D0**0.5D0)
    !--------------------------------------------------------------------
END SUBROUTINE kernel_lr
!==================================================================


!==================================================================
SUBROUTINE gaussian(ps_r, rho_ps, p)
    USE kinds,            ONLY  :  DP
    USE exx_module,       ONLY  :  me_rs
    !
    IMPLICIT  NONE
    !--------------------------------------------------------------------
    INTEGER                :: ps_r(6)
    REAL(DP)               :: rho_ps(ps_r(1):ps_r(4),ps_r(2):ps_r(5),ps_r(3):ps_r(6))
    REAL(DP)               :: p
    !--------------------------------------------------------------------
    REAL(DP)               :: r
    INTEGER                :: i,j,k
    !--------------------------------------------------------------------

    !--------------------------------------------------------------------
    !$omp parallel do private(i,j,k,r) schedule(guided)
    !--------------------------------------------------------------------
    DO k=ps_r(3),ps_r(6)
      DO j=ps_r(2),ps_r(5)
        DO i=ps_r(1),ps_r(4)
        !----------------------------------------------------------------
        r = me_rs(1,i,j,k)
        !----------------------------------------------------------------
        rho_ps(i,j,k) = (p/3.1415926535897932D0)**1.5D0 * exp(-p * r**2.0D0)
        !----------------------------------------------------------------
        END DO
      END DO
    END DO
    !--------------------------------------------------------------------
END SUBROUTINE gaussian
!==================================================================


!!==================================================================
!SUBROUTINE shrink_and_interpolate(me_r, pot_me)
!    USE kinds,            ONLY  :  DP
!    !
!    IMPLICIT  NONE
!    !--------------------------------------------------------------------
!    INTEGER                :: me_r(6)
!    REAL(DP)               :: pot_me(me_r(1):me_r(4),me_r(2):me_r(5),me_r(3):me_r(6))
!    !--------------------------------------------------------------------
!    INTEGER                :: npd(3)
!    INTEGER                :: mc_r(6)
!    INTEGER                :: mf_r(6)
!    REAL(DP), ALLOCATABLE  :: pot_mc(:,:,:)
!    REAL(DP), ALLOCATABLE  :: pot_mf(:,:,:)
!    INTEGER                :: i,j,k
!    !--------------------------------------------------------------------
!    mc_r(1)=1
!    mc_r(2)=1
!    mc_r(3)=1
!    mc_r(4)=(me_r(4)-me_r(1)+2)/2-1
!    mc_r(5)=(me_r(5)-me_r(2)+2)/2-1
!    mc_r(6)=(me_r(6)-me_r(3)+2)/2-1
!    !--------------------------------------------------------------------
!    mf_r(1)=me_r(1)-1
!    mf_r(2)=me_r(2)-1
!    mf_r(3)=me_r(3)-1
!    mf_r(4)=me_r(4)+1
!    mf_r(5)=me_r(5)+1
!    mf_r(6)=me_r(6)+1
!    !--------------------------------------------------------------------
!    npd(1) = ((me_r(4)-me_r(1)+1) - (mc_r(4)-mc_r(1)+1)) / 2
!    npd(2) = ((me_r(5)-me_r(2)+1) - (mc_r(5)-mc_r(2)+1)) / 2
!    npd(3) = ((me_r(6)-me_r(3)+1) - (mc_r(6)-mc_r(3)+1)) / 2
!    !--------------------------------------------------------------------
!
!    !--------------------------------------------------------------------
!    ALLOCATE(pot_mc(mc_r(1):mc_r(4), mc_r(2):mc_r(5), mc_r(3):mc_r(6))); pot_mc=0.0D0
!    ALLOCATE(pot_mf(mf_r(1):mf_r(4), mf_r(2):mf_r(5), mf_r(3):mf_r(6))); pot_mf=0.0D0
!    !--------------------------------------------------------------------
!
!    !--------------------------------------------------------------------
!    pot_mc(:,:,:) = pot_me(me_r(1)+npd(1):me_r(4)-npd(1), me_r(2)+npd(2):me_r(5)-npd(2), me_r(3)+npd(3):me_r(6)-npd(3))
!    !--------------------------------------------------------------------
!
!    !--------------------------------------------------------------------
!    CALL interp_c2f_spline(mc_r,mc_r,pot_mc,me_r,mf_r,pot_mf)
!    !--------------------------------------------------------------------
!
!    !--------------------------------------------------------------------
!    pot_me(:,:,:) = pot_mf(me_r(1):me_r(4), me_r(2):me_r(5), me_r(3):me_r(6))
!    !--------------------------------------------------------------------
!
!!   !--------------------------------------------------------------------
!!   DO k=me_r(3),me_r(6)
!!     DO j=me_r(2),me_r(5)
!!       DO i=me_r(1),me_r(4)
!!       !----------------------------------------------------------------
!!       WRITE(*,"(I4,I4,I4,F15.11)") i,j,k,pot_me(i,j,k)
!!       !----------------------------------------------------------------
!!       END DO
!!     END DO
!!   END DO
!!   !--------------------------------------------------------------------
!    DEALLOCATE(pot_mc)
!    DEALLOCATE(pot_mf)
!    !--------------------------------------------------------------------
!END SUBROUTINE shrink_and_interpolate
!==================================================================

!!==================================================================
!subroutine setplm_1d(x, y, plm_1d)
!  use exx_module,  only : lmax, lm_mx
!  !
!  ! This subroutine calculates all of the associated Legendre polynomials up 
!  ! to a supplied lmax (maximum lmax is 9). 
!  !
!  ! HK: TODO: add evaluation for general plm (recursion); rename subroutine to compute_legendre_polynomial_plm 
!  implicit none
!  !
!  ! INPUT VARIABLES
!  ! cos(theta),sin(theta), for a given grid point
!  real(8) :: x,y
!  ! order of multipole expansion
!  ! OUTPUT VARIABLES
!  ! array containing P_{lm}, the associated Legendre polynomials
!  real(8) :: plm_1d(lm_mx)
!  ! WORK VARIABLES:
!  ! powers of x, y: xn=x^n, yn=y^n
!  real(8) :: x2, x3, x4, x5, x6, x7, x8, x9
!  real(8) :: y2, y3, y4, y5, y6, y7, y8, y9
!  !
!  ! HK: TODO: instead of calling lm_1d many times we should simply loop through l, m as these routine is already ordered along lm_1d
!  plm(lm_1d(0,0)) = 1.0
!  if (lmax .ge. 1) then
!    plm(lm_1d(1,0)) = x
!    plm(lm_1d(1,1)) = -y
!  end if
!  if (lmax .ge. 2) then
!    x2 = x*x
!    y2 = y*y
!    plm(lm_1d(2,0)) =  1.5d0*x2 - 0.5d0
!    plm(lm_1d(2,1)) = -3.0d0*x*y
!    plm(lm_1d(2,2)) =  3.0d0*y2
!  end if
!  if (lmax .ge. 3) then
!    x3 = x2*x
!    y3 = y2*y
!    plm(lm_1d(3,0)) =   2.5d0*x3 - 1.5d0*x
!    plm(lm_1d(3,1)) = (-7.5d0*x2 + 1.5d0)*y
!    plm(lm_1d(3,2)) =  15.0d0*x*y2
!    plm(lm_1d(3,3)) = -15.0d0*y3
!  end if
!  if (lmax .ge. 4) then
!    x4 = x2*x2
!    y4 = y2*y2
!    plm(lm_1d(4,0)) =  4.375d0*x4 - 3.75d0*x2 + 0.375d0
!    plm(lm_1d(4,1)) = (-17.5d0*x3 + 7.5d0*x)*y
!    plm(lm_1d(4,2)) = ( 52.5d0*x2 - 7.5d0  )*y2
!    plm(lm_1d(4,3)) = -105.0d0*x*y3
!    plm(lm_1d(4,4)) =  105.0d0*y4
!  end if
!  if (lmax .ge. 5) then
!    x5 = x3*x2
!    y5 = y3*y2
!    plm(lm_1d(5,0)) =  7.875d0*x5 - 8.75d0*x3 + 1.875d0*x
!    plm(lm_1d(5,1)) = (-39.375d0*x4 + 26.25d0*x2 - 1.875d0)*y
!    plm(lm_1d(5,2)) = ( 157.5d0*x3 - 52.5d0*x)*y2
!    plm(lm_1d(5,3)) = (-472.5d0*x2 + 52.5d0)  *y3
!    plm(lm_1d(5,4)) =  945.0d0*x*y4
!    plm(lm_1d(5,5)) = -945.0d0*y5
!  end if
!  if (lmax .ge. 6) then
!    x6 = x3*x3
!    y6 = y3*y3
!    plm(lm_1d(6,0)) = 14.4375d0*x6 - 19.6875d0*x4 + 6.5625d0*x2 - 0.3125d0
!    plm(lm_1d(6,1)) = (-86.625d0*x5 + 78.75d0*x3 - 13.125d0*x)*y
!    plm(lm_1d(6,2)) = ( 433.125d0*x4 - 236.25d0*x2 + 13.125d0)*y2
!    plm(lm_1d(6,3)) = (-1732.5d0*x3 + 472.5d0*x            )*y3
!    plm(lm_1d(6,4)) = ( 5197.5d0*x2 - 472.5d0              )*y4
!    plm(lm_1d(6,5)) = -10395.0d0*x*y5
!    plm(lm_1d(6,6)) =  10395.0d0*y6
!  end if
!  if (lmax .ge. 7) then
!    x7 = x4*x3
!    y7 = y4*y3
!    plm(lm_1d(7,0)) = 26.8125d0*x7 - 43.3125d0*x5 + 19.6875d0*x3 - 2.1875d0*x
!    plm(lm_1d(7,1)) = -187.6875d0*x6*y + 216.5625d0*x4*y - 59.0625d0*x2*y + 2.1875d0*y
!    plm(lm_1d(7,2)) = 1126.125d0*x5*y2 - 866.25d0*x3*y2 + 118.125d0*x*y2
!    plm(lm_1d(7,3)) = -5630.625d0*x4*y3 + 2598.75d0*x2*y3 - 118.125d0*y3
!    plm(lm_1d(7,4)) = 22522.5d0*x3*y4 - 5197.5d0*x*y4
!    plm(lm_1d(7,5)) = -67567.5d0*x2*y5 + 5197.5d0*y5
!    plm(lm_1d(7,6)) = 135135.0d0*x*y6
!    plm(lm_1d(7,7)) = -135135.0d0*y7
!  end if
!  if (lmax .ge. 8) then
!    x8 = x4*x4
!    y8 = y4*y4
!    plm(lm_1d(8,0)) = 50.2734375d0*x8 - 93.84375d0*x6 + 54.140625d0*x4 - 9.84375d0*x2 + 0.2734375d0
!    plm(lm_1d(8,1)) = -402.1875d0*x7*y + 563.0625d0*x5*y - 216.5625d0*x3*y + 19.6875d0*x*y
!    plm(lm_1d(8,2)) = 2815.3125d0*x6*y2 - 2815.3125d0*x4*y2 + 649.6875d0*x2*y2 - 19.6875d0*y2
!    plm(lm_1d(8,3)) = -16891.875d0*x5*y3 + 11261.25d0*x3*y3 - 1299.375d0*x*y3
!    plm(lm_1d(8,4)) = 84459.375d0*x4*y4 - 33783.75d0*x2*y4 + 1299.375d0*y4
!    plm(lm_1d(8,5)) = -337837.5d0*x3*y5 + 67567.5d0*x*y5
!    plm(lm_1d(8,6)) = 1013512.5d0*x2*y6 - 67567.5d0*y6
!    plm(lm_1d(8,7)) = -2027025.0d0*x*y7
!    plm(lm_1d(8,8)) = 2027025.0d0*y8
!  end if
!  if (lmax .ge. 9) then
!    y9 = y5*y4
!    x9 = x5*x4
!    plm(lm_1d(9,0)) = 94.9609375d0*x9 - 201.09375d0*x7 + 140.765625d0*x5 - 36.09375d0*x3 + 2.4609375d0*x
!    plm(lm_1d(9,1)) = -854.6484375d0*x8*y + 1407.65625d0*x6*y - 703.828125d0*x4*y + 108.28125d0*x2*y - 2.4609375d0*y
!    plm(lm_1d(9,2)) = 6837.1875d0*x7*y2 - 8445.9375d0*x5*y2 + 2815.3125d0*x3*y2 - 216.5625d0*x*y2
!    plm(lm_1d(9,3)) = -47860.3125d0*x6*y3 + 42229.6875d0*x4*y3 - 8445.9375d0*x2*y3 + 216.5625d0*y3
!    plm(lm_1d(9,4)) = 287161.875d0*x5*y4 - 168918.75d0*x3*y4 + 16891.875d0*x*y4
!    plm(lm_1d(9,5)) = -1435809.375d0*x4*y5 + 506756.25d0*x2*y5 - 16891.875d0*y5
!    plm(lm_1d(9,6)) = 5743237.5d0*x3*y6 - 1013512.5d0*x*y6
!    plm(lm_1d(9,7)) = -17229712.5d0*x2*y7 + 1013512.5d0*y7
!    plm(lm_1d(9,8)) = 34459425.0d0*x*y8
!    plm(lm_1d(9,9)) = -34459425.0d0*y9
!  end if
!  !
!  !--------------------------------------------------------------------
!  return
!  !--------------------------------------------------------------------
!end subroutine setplm_1d
!!==================================================================
!
!
!! HK: For better efficiency, we could also tabulate lm_1d <-> l,m relationship upto lmax (initialize stage): TODO
!integer function lm_1d(l,m)
!  implicit none
!  integer, intent(in) :: l, m
!  lm_1d = (l*(l+1))/2 + 1 + m
!end function lm_1d
!!TODO: add a note...
!
!function lm_1d_to_l_and_m(lm_1d) result(lm_v)
!  implicit none
!  integer, intent(in) :: lm_1d
!  integer :: l,m
!  integer :: lm_v (2)
!  logical :: is_right_l
!  is_right_l = .false.
!  do l = floor(sqrt(2.0*lm_1d))+1, ceiling(sqrt(2.0*lm_1d)+0.5)
!    m = lm_1d - lm_1d(l,0)
!    is_right_l = (m.ge.0.and.m.le.l)
!    if (is_right_l) exit
!  end do ! l
!  if (.not.is_right_l) then
!    write(*,*) 'Error: cannot find lm_1d_to_l_and_m relation. STOP!'
!    stop
!  else
!    lm_v = (/l, m/)
!  end if
!end function lm_1d_to_l_and_m
