!
! Copyright (C) 2024 Quantum ESPRESSO Foundation
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
MODULE atwfc_mod
  !
  !! Variables and routines atomic wavefunctions in numerical form
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
  PUBLIC :: init_tab_atwfc
  PUBLIC :: deallocate_tab_atwfc
  PUBLIC :: scale_tab_atwfc
  PUBLIC :: interp_atwfc
  PUBLIC :: interp_atdwfc
  !
  SAVE
  !
  INTEGER :: nqx = 0
  !! size of interpolation table
  REAL(DP), PARAMETER:: dq = 0.01_dp
  !! grid step for interpolation table
  REAL(DP) :: qmax = 0.0_dp 
  !! max q covered by the interpolation table
  REAL(DP), ALLOCATABLE :: tab_atwfc(:,:,:)
  !! interpolation table for numerical beta functions in reciprocal space
  !
CONTAINS
  !
  !-----------------------------------------------------------------------
  SUBROUTINE init_tab_atwfc( qmax_, omega, comm, ierr)
  !-----------------------------------------------------------------------
  !! This routine computes a table with the radial Fourier transform 
  !! of the atomic wavefunctions.
  !
  USE upf_kinds,    ONLY : DP
  USE atom,         ONLY : rgrid, msh
  USE upf_const,    ONLY : fpi
  USE uspp_param,   ONLY : nsp, upf, nwfcm, nsp
  USE mp,           ONLY : mp_sum
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
  INTEGER :: nt, nb, iq, ir, l, startq, lastq, ndm
  !
  REAL(DP), ALLOCATABLE :: aux(:), vchi(:)
  REAL(DP) :: vqint, pref, q
  !
  IF ( .NOT. ALLOCATED(tab_atwfc) ) THEN
     !! table not yet allocated
     qmax = qmax_
     ierr = 0
  ELSE IF ( qmax_ > qmax ) THEN
     !! table ìs allocated but dimension insufficient: re-allocate
     !! (with some margin so that this does not happen too often)
     !$acc exit data delete(tab_atwfc)
     DEALLOCATE ( tab_atwfc )
     qmax = qmax_ + MAX(dq*100,qmax_-qmax)
     ierr =-1
  ELSE
     !! table already computed: exit
     ierr =-2
     RETURN
  END IF
  nqx = NINT( qmax/dq + 4)
  allocate(tab_atwfc(nqx,nwfcm,nsp))
  !$acc enter data create(tab_atwfc)
  !
  ndm = MAXVAL(msh(1:nsp))
  ALLOCATE( aux(ndm), vchi(ndm) )
  !
  ! chiq = radial fourier transform of atomic orbitals chi
  !
  pref = fpi / SQRT(omega)
  ! needed to normalize atomic wfcs (not a bad idea in general and 
  ! necessary to compute correctly lda+U projections)
  CALL divide( comm, nqx, startq, lastq )
  !
  tab_atwfc(:,:,:) = 0.0_DP
  !
  DO nt = 1, nsp
     DO nb = 1, upf(nt)%nwfc
        !
        IF (upf(nt)%oc(nb) >= 0.0_DP) THEN
           l = upf(nt)%lchi (nb)
           !
           DO iq = startq, lastq
              q = dq * (iq - 1)
              CALL sph_bes( msh(nt), rgrid(nt)%r, q, l, aux )
              DO ir = 1, msh(nt)
                 vchi(ir) = upf(nt)%chi(ir,nb) * aux(ir) * rgrid(nt)%r(ir)
              ENDDO
              CALL simpson( msh(nt), vchi, rgrid(nt)%rab, vqint )
              tab_atwfc( iq, nb, nt ) = vqint * pref
           ENDDO
           !
        ENDIF
        !
     ENDDO
  ENDDO
  !
  CALL mp_sum( tab_atwfc, comm )
  !
  !$acc update device(tab_atwfc)
  !
  DEALLOCATE( aux, vchi )
  !
  RETURN
  !
END SUBROUTINE init_tab_atwfc
  !
  !-----------------------------------------------------------------------
  SUBROUTINE interp_atwfc ( npw, qg, nwfcm, chiq )
  !-----------------------------------------------------------------------
  !
  ! computes chiq: radial fourier transform of atomic orbitals chi
  !
  USE uspp_param, ONLY : upf, nsp
  !
  INTEGER, INTENT(IN)  :: npw
  INTEGER, INTENT(IN)  :: nwfcm
  REAL(dp), INTENT(IN) :: qg(npw)
  REAL(dp), INTENT(OUT):: chiq(npw,nwfcm,nsp)
  !
  INTEGER :: nt, nb, ig
  INTEGER :: i0, i1, i2, i3
  REAL(dp):: qgr, px, ux, vx, wx
  !
  !$acc data present_or_copyin(qg) present_or_copyout(chiq) present(tab_atwfc)
  DO nt = 1, nsp
     DO nb = 1, upf(nt)%nwfc
        IF ( upf(nt)%oc(nb) >= 0.d0 ) THEN
           !$acc parallel loop
           DO ig = 1, npw
              qgr = qg(ig)
              px = qgr / dq - DBLE(INT(qgr/dq))
              ux = 1.d0 - px
              vx = 2.d0 - px
              wx = 3.d0 - px
              i0 = INT(qgr/dq) + 1
              i1 = i0 + 1
              i2 = i0 + 2
              i3 = i0 + 3
              chiq(ig,nb,nt) = &
                     tab_atwfc(i0,nb,nt) * ux * vx * wx / 6.d0 + &
                     tab_atwfc(i1,nb,nt) * px * vx * wx / 2.d0 - &
                     tab_atwfc(i2,nb,nt) * px * ux * wx / 2.d0 + &
                     tab_atwfc(i3,nb,nt) * px * ux * vx / 6.d0
           END DO
           !
        END IF
     END DO
  END DO
  !$acc end data
END SUBROUTINE interp_atwfc
!
!-----------------------------------------------------------------------
SUBROUTINE interp_atdwfc ( npw, qg, nwfcm, dchiq )
  !-----------------------------------------------------------------------
  !
  ! computes dchi/dq
  !
  USE upf_kinds,  ONLY : dp
  USE uspp_param, ONLY : upf, nsp
  !
  INTEGER, INTENT(IN)  :: npw
  INTEGER, INTENT(IN)  :: nwfcm
  REAL(dp), INTENT(IN) :: qg(npw)
  REAL(dp), INTENT(OUT):: dchiq(npw,nwfcm,nsp)
  !
  INTEGER :: nt, nb, ig
  INTEGER :: i0, i1, i2, i3
  REAL(dp):: px, ux, vx, wx
  !
  !$acc data present_or_copyin(qg) present_or_copyout(dchiq) present(tab_atwfc)
  DO nt=1,nsp
     DO nb=1,upf(nt)%nwfc
        IF (upf(nt)%oc(nb) >= 0.d0) THEN
           !$acc parallel loop
           DO ig = 1, npw
              px = qg(ig) / dq - INT(qg(ig)/dq)
              ux = 1.d0 - px
              vx = 2.d0 - px
              wx = 3.d0 - px
              i0 = qg(ig) / dq + 1
              i1 = i0 + 1
              i2 = i0 + 2
              i3 = i0 + 3
              dchiq(ig,nb,nt) = &
                 ( tab_atwfc(i0, nb, nt) * (-vx*wx-ux*wx-ux*vx)/6.d0 + &
                   tab_atwfc(i1, nb, nt) * (+vx*wx-px*wx-px*vx)/2.d0 - &
                   tab_atwfc(i2, nb, nt) * (+ux*wx-px*wx-px*ux)/2.d0 + &
                   tab_atwfc(i3, nb, nt) * (+ux*vx-px*vx-px*ux)/6.d0 )/dq
           ENDDO
        ENDIF
     ENDDO
  ENDDO
  !$acc end data
  END SUBROUTINE interp_atdwfc
  !
  subroutine allocate_tab_atwfc ( nqx_, nwfcm, nsp)
    implicit none
    INTEGER, INTENT(IN) :: nqx_, nwfcm, nsp
    !
    nqx = nqx_
    allocate(tab_atwfc(nqx_,nwfcm,nsp))
    !$acc enter data create(tab_atwfc)
    !
  end subroutine allocate_tab_atwfc
  !
  subroutine deallocate_tab_atwfc ()
    implicit none
    !
    !$acc exit data delete(tab_atwfc)
    if( allocated( tab_atwfc ) )  deallocate( tab_atwfc )
    !
  end subroutine deallocate_tab_atwfc
  !
  subroutine scale_tab_atwfc( vol_ratio_m1 )
    ! vol_ratio_m1 = omega_old / omega
    implicit none
    real(DP), intent(in) :: vol_ratio_m1
    !
    if ( allocated(tab_atwfc) ) tab_atwfc(:,:,:) = tab_atwfc(:,:,:) * SQRT(vol_ratio_m1)
    !$acc update device ( tab_atwfc)
  end subroutine scale_tab_atwfc
   !
END MODULE atwfc_mod
