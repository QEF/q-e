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
    REAL(DP)                  :: normfactor
    REAL(DP)                  :: rhosum
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
    attributes(device) :: qlm
    attributes(device) :: potme, rhops
#endif
    

    !-----------------------------------------------------------------------
    ! --- external functions ---
    !-----------------------------------------------------------------------
    REAL(DP), EXTERNAL     :: dnrm2
    !-----------------------------------------------------------------------
    write(*,*) 'MCA: Starting getvofr'
    !-----------------------------------------------------------------------
    ! --- initialize ---
    !-----------------------------------------------------------------------
    if (.not.allocated(rho_ps)) allocate( rho_ps(n_ps))
    if (.not.allocated(pot_ps)) then
      allocate( pot_ps(n_ps)        ); pot_ps  = 0.0d0
    end if
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
    ALLOCATE(qlm(0:lmax, 0:lmax))
    !---------------------------------------------------------------------------
    ! Then, we calculate the potential outside the inner sphere, which will
    ! also give the boundary values for potential inside the sphere.
    !---------------------------------------------------------------------------
    CALL start_clock('getvofr_qlm')
    write(*,*) 'MCA: getqlm'
    CALL getqlm(ps_r, hcub, rho_ps, qlm)
    CALL stop_clock('getvofr_qlm')
    !-----------------------------------------------------------------------
    CALL start_clock('getvofr_bound')
    write(*,*) 'MCA: exx_boundary'
    CALL exx_boundaryv(me_r, ps_r, potme, qlm)
    rhosum = sum(potme)
    write (*,*), 'sum potme: ', rhosum
    CALL stop_clock('getvofr_bound')
    !-----------------------------------------------------------------------
    CALL start_clock('getvofr_geterho')
    write(*,*) 'MCA: geterho'
    CALL geterho(me_r, ps_r, potme, rho_ps)
    rhosum = sum(rho_ps)
    write (*,*), 'sum rho_ps: ', rhosum
    CALL stop_clock('getvofr_geterho')
    !========================================================================
    !TODO : copy out rho_ps before the next step
    !---------------------------------------------------------------------------
    ! --- Poisson solver ---
    !---------------------------------------------------------------------------
    CALL start_clock('getvofr_solver')
    !---------------------------------------------------------------------------

    write(*,*) 'MCA: cg_solver'
    call cg_solver_stdcg

    !---------------------------------------------------------------------------
    CALL stop_clock('getvofr_solver')
    !---------------------------------------------------------------------------

    !---------------------------------------------------------------------------
    write(*,*) 'MCA: ps2me'
    CALL ps2me(me_r, ps_r, potme, pot_ps)
    write(*,*) 'MCA: after ps2me'
    !---------------------------------------------------------------------------
    DEALLOCATE( qlm    )
    !---------------------------------------------------------------------------
    RETURN
    !---------------------------------------------------------------------------
  contains

    subroutine  cg_solver_stdcg()
      implicit none
#ifdef __CUDA
      CALL CGMIC_STDCG(nstep, ncb, poisson_eps, fbsscale, coemicf_d, coeke_d, rho_ps, pot_ps)
#else
      CALL CGMIC_STDCG(nstep, ncb, poisson_eps, fbsscale, coemicf, coeke, rho_ps, pot_ps)
#endif
      !
      return
    end subroutine cg_solver_stdcg

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
        !$omp parallel do collapse(3) private(i,j,k,x,y,z,sqrxy,costheta,sintheta) reduction(+:qlm_tmp)
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
    !
    !---------------------------------------------------------------------------
    RETURN
    !---------------------------------------------------------------------------
END SUBROUTINE getqlm
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
    !$omp parallel do private(x,y,z,sqrxy,costheta,sintheta,cpot_r) schedule(guided)
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

