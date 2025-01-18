!
! Copyright (C) 2024 Quantum ESPRESSO Foundation
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
MODULE beta_mod
  !
  !! Variables and routines for nonlocal beta functions in numerical form
  !! Contains generation of interpolation tables in reciprocal space,
  !! interpolation routines and other utility routines
  !! Code moved to upflib and restructured by Paolo Giannozzi, 2024
  !
  USE upf_kinds,    ONLY : dp
  USE upf_const,    ONLY : fpi, e2
  !
  IMPLICIT NONE
  !
  PRIVATE
  PUBLIC :: init_tab_beta
  PUBLIC :: deallocate_tab_beta
  PUBLIC :: scale_tab_beta
  PUBLIC :: interp_beta
  PUBLIC :: interp_dbeta
  !
  SAVE
  !
  INTEGER :: nqx = 0
  !! size of interpolation table
  REAL(DP), PARAMETER:: dq = 0.01_dp
  !! grid step for interpolation table
  REAL(DP) :: qmax = 0.0_dp 
  !! max q covered by the interpolation table
  REAL(DP), ALLOCATABLE :: tab_beta(:,:,:)
  !! interpolation table for numerical beta functions in reciprocal space
  !
CONTAINS
!
!----------------------------------------------------------------------
SUBROUTINE init_tab_beta ( qmax_, omega, comm, ierr ) 
  !----------------------------------------------------------------------
  !
  ! Compute interpolation table for beta(G) radial functions
  !
  USE upf_kinds,    ONLY : dp
  USE atom,         ONLY : rgrid
  USE uspp_param,   ONLY : upf, lmaxq, nbetam, nsp
  USE mp,           ONLY : mp_sum
  USE m_gth,        ONLY : mk_ffnl_gth
  !
  IMPLICIT NONE
  !
  REAL(dp), INTENT(IN) :: qmax_
  !! Interpolate q up to qmax_ (sqrt(Ry), q^2 is an energy)
  REAL(dp), INTENT(IN) :: omega
  !! Unit-cell volume
  INTEGER, INTENT(IN)  :: comm
  !! MPI communicator, to split the workload
  INTEGER, INTENT(OUT) :: ierr
  !! error code: ierr = 0 if interpolation table (IT) was allocated
  !!             ierr =-1 if IT had insufficient size and was re-allocated
  !!             ierr =-2 if IT was already present and nothing is done
  !
  INTEGER :: ndm, startq, lastq, nt, l, nb, iq, ir
  REAL(dp) :: qi
  ! q-point grid for interpolation
  REAL(dp) :: pref
  ! the prefactor of the Q functions
  real(DP) ::  vqint, d1
  !
  REAL(dp), allocatable :: aux (:)
  ! work space
  REAL(dp), allocatable :: besr(:)
  ! work space
  !
  IF ( .NOT. ALLOCATED(tab_beta) ) THEN
     !! table not yet allocated
     qmax = qmax_
     ierr = 0
  ELSE IF ( qmax_ > qmax ) THEN
     !! table ìs allocated but dimension insufficient: re-allocate
     !! (with some margin so that this does not happen too often)
     !$acc exit data delete(tab_beta)
     DEALLOCATE ( tab_beta )
     qmax = qmax_ + MAX(dq*100,qmax_-qmax)
     ierr =-1
  ELSE
     !! table already computed: exit
     ierr =-2
     RETURN
  END IF
  nqx = INT( qmax/dq + 4)
  allocate(tab_beta(nqx,nbetam,nsp))
  !$acc enter data create(tab_beta)
  ndm = MAXVAL ( upf(:)%kkbeta )
  allocate( aux (ndm) )
  allocate (besr( ndm))
  pref = fpi / sqrt (omega)
  call divide (comm, nqx, startq, lastq)
  tab_beta (:,:,:) = 0.d0
  do nt = 1, nsp
     do nb = 1, upf(nt)%nbeta
        l = upf(nt)%lll (nb)
        do iq = startq, lastq
           qi = (iq - 1) * dq
           if ( upf(nt)%is_gth ) then
              CALL mk_ffnl_gth( nt, nb, 1, omega, [ qi ] , tab_beta(iq,nb,nt) )
           else
              call sph_bes (upf(nt)%kkbeta, rgrid(nt)%r, qi, l, besr)
              do ir = 1, upf(nt)%kkbeta
                 aux (ir) = upf(nt)%beta (ir, nb) * besr (ir) * rgrid(nt)%r(ir)
              enddo
              call simpson (upf(nt)%kkbeta, aux, rgrid(nt)%rab, vqint)
              tab_beta (iq, nb, nt) = vqint * pref
           end if
        enddo
     enddo
  enddo
  deallocate (besr)
  deallocate (aux)
  !
  call mp_sum(  tab_beta, comm )
!$acc update device (tab_beta)
  !
END SUBROUTINE init_tab_beta
!
!----------------------------------------------------------------------
SUBROUTINE interp_beta( nt, npw_, qg, vq )
  !----------------------------------------------------------------------
  !
  USE upf_kinds,  ONLY : dp
  USE uspp_param, ONLY : upf, nbetam
  !
  implicit none
  integer, intent(in) :: nt, npw_
  real(dp), intent(in ) :: qg(npw_)
  real(dp), intent(out) :: vq(npw_,nbetam)
  !
  integer :: i0, i1, i2, i3, nbnt, nb, ig
  real(dp):: qgr, px, ux, vx, wx
  !
  nbnt = upf(nt)%nbeta
  !$acc data present (tab_beta) present_or_copyin (qg) present_or_copyout (vq)
  !$acc parallel loop collapse(2)
  do nb = 1, nbnt
     DO ig = 1, npw_
        qgr = qg(ig)
        px = qgr / dq - DBLE(INT(qgr/dq))
        ux = 1.0_dp - px
        vx = 2.0_dp - px
        wx = 3.0_dp - px
        i0 = INT(qgr/dq) + 1
        i1 = i0 + 1
        i2 = i0 + 2
        i3 = i0 + 3
        if ( i3 <= nqx ) then
           vq(ig,nb) = &
             tab_beta(i0,nb,nt) * ux * vx * wx / 6.0_dp + &
             tab_beta(i1,nb,nt) * px * vx * wx / 2.0_dp - &
             tab_beta(i2,nb,nt) * px * ux * wx / 2.0_dp + &
             tab_beta(i3,nb,nt) * px * ux * vx / 6.0_dp
        else
           !! This case should never happen if tab_beta is properly allocated
           !! (setting q_max to be large enough) - for compatibility with GWW
           vq(ig,nb) = 0.0_dp
        end if
     END DO
  END DO
  !$acc end data
  !----------------------------------------------------------------------
END SUBROUTINE interp_beta
!----------------------------------------------------------------------
!
!----------------------------------------------------------------------
SUBROUTINE interp_dbeta( nt, npw, qg, vq )
  !----------------------------------------------------------------------
  !
  USE upf_kinds,  ONLY : dp
  USE uspp_param, ONLY : upf, nbetam
  !
  implicit none
  integer, intent(in) :: nt, npw
  real(dp), intent(in ) :: qg(npw)
  real(dp), intent(out) :: vq(npw,nbetam)
  !
  integer :: i0, i1, i2, i3, nbnt, nb, ig
  real(dp):: qgr, px, ux, vx, wx
  !
  nbnt = upf(nt)%nbeta
  !$acc data present (tab_beta) present_or_copyin (qg) present_or_copyout (vq)
  !$acc parallel loop collapse(2)
  DO nb = 1, nbnt
     DO ig = 1, npw
        qgr = qg(ig)
        px = qgr / dq - INT(qgr/dq)
        ux = 1.0_dp - px
        vx = 2.0_dp - px
        wx = 3.0_dp - px
        i0 = qgr / dq + 1
        i1 = i0 + 1
        i2 = i0 + 2
        i3 = i0 + 3
        IF ( i3 <= nqx ) THEN
            vq(ig,nb) = ( tab_beta(i0,nb,nt) * (-vx*wx-ux*wx-ux*vx)/6.0_dp + &
                          tab_beta(i1,nb,nt) * (+vx*wx-px*wx-px*vx)/2.0_dp - &
                          tab_beta(i2,nb,nt) * (+ux*wx-px*wx-px*ux)/2.0_dp + &
                          tab_beta(i3,nb,nt) * (+ux*vx-px*vx-px*ux)/6.0_dp )/dq
        ELSE
            vq(ig,nb) = 0.0_dp 
        END IF
     ENDDO
  END DO
  !$acc end data
  !----------------------------------------------------------------------
END SUBROUTINE interp_dbeta
!
  subroutine deallocate_tab_beta ()
    implicit none
    !
    !$acc exit data delete(tab_beta)
    if( allocated( tab_beta ) )  deallocate( tab_beta )
    !
  end subroutine deallocate_tab_beta
  !
  subroutine scale_tab_beta( vol_ratio_m1 )
    ! vol_ratio_m1 = omega_old / omega
    implicit none
    real(DP), intent(in) :: vol_ratio_m1
    !
    tab_beta(:,:,:) = tab_beta(:,:,:) * SQRT(vol_ratio_m1)
!$acc update device ( tab_beta)
  end subroutine scale_tab_beta
   !
END MODULE beta_mod
