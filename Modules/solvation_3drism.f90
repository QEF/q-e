!
! Copyright (C) 2015-2016 Satomichi Nishihara
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!---------------------------------------------------------------------------
SUBROUTINE solvation_3drism(rismt, ierr)
  !---------------------------------------------------------------------------
  !
  ! ... setup charge density, solvation potential and solvation energy for DFT,
  ! ... which are derived from 3D-RISM
  !
  USE cell_base,      ONLY : tpiba2, omega
  USE constants,      ONLY : fpi, e2
  USE err_rism,       ONLY : IERR_RISM_NULL, IERR_RISM_INCORRECT_DATA_TYPE
  USE fft_interfaces, ONLY : fwfft
  USE io_global,      ONLY : stdout
  USE kinds,          ONLY : DP
  USE mp,             ONLY : mp_sum
  USE rism,           ONLY : rism_type, ITYPE_3DRISM
  USE solvmol,        ONLY : get_nuniq_in_solVs, solVs, iuniq_to_nsite, &
                           & iuniq_to_isite, isite_to_isolV, isite_to_iatom
  !
  IMPLICIT NONE
  !
  TYPE(rism_type), INTENT(INOUT) :: rismt
  INTEGER,         INTENT(OUT)   :: ierr
  !
  INTEGER                  :: nq
  INTEGER                  :: iq
  INTEGER                  :: iiq
  INTEGER                  :: iv
  INTEGER                  :: nv
  INTEGER                  :: isolV
  INTEGER                  :: iatom
  INTEGER                  :: ir
  INTEGER                  :: ig
  REAL(DP)                 :: rhov
  REAL(DP)                 :: qv
  REAL(DP)                 :: ntmp
  REAL(DP)                 :: gg0
  REAL(DP)                 :: domega
  REAL(DP)                 :: fac
  REAL(DP)                 :: rhotot
  REAL(DP),    ALLOCATABLE :: rhor(:)
  COMPLEX(DP), ALLOCATABLE :: aux(:)
  !
  ! ... number of sites in solvents
  nq = get_nuniq_in_solVs()
  !
  ! ... check data type
  IF (rismt%itype /= ITYPE_3DRISM) THEN
    ierr = IERR_RISM_INCORRECT_DATA_TYPE
    RETURN
  END IF
  !
  IF (rismt%mp_site%nsite < nq) THEN
    ierr = IERR_RISM_INCORRECT_DATA_TYPE
    RETURN
  END IF
  !
  IF (rismt%nr < rismt%dfft%nnr) THEN
    ierr = IERR_RISM_INCORRECT_DATA_TYPE
    RETURN
  END IF
  !
  IF (rismt%ng < rismt%gvec%ngm) THEN
    ierr = IERR_RISM_INCORRECT_DATA_TYPE
    RETURN
  END IF
  !
  ! ... allocate working memory
  IF (rismt%dfft%nnr > 0) THEN
    ALLOCATE(rhor(rismt%dfft%nnr))
    ALLOCATE(aux(rismt%dfft%nnr))
  END IF
  !
  ! ... set variables
  domega = omega / DBLE(rismt%dfft%nr1) &
               & / DBLE(rismt%dfft%nr2) &
               & / DBLE(rismt%dfft%nr3)
  !
  ! ... make nsol, qsol
  DO iq = rismt%mp_site%isite_start, rismt%mp_site%isite_end
    iiq   = iq - rismt%mp_site%isite_start + 1
    iv    = iuniq_to_isite(1, iq)
    nv    = iuniq_to_nsite(iq)
    isolV = isite_to_isolV(iv)
    iatom = isite_to_iatom(iv)
    rhov  = DBLE(nv) * solVs(isolV)%density
    qv    = solVs(isolV)%charge(iatom)
    fac   = rhov * domega
    !
    ntmp = 0.0_DP
!$omp parallel do default(shared) private(ir) reduction(+:ntmp)
    DO ir = 1, rismt%dfft%nnr
      ntmp = ntmp + fac * rismt%gr(ir, iiq)
    END DO
!$omp end parallel do
    rismt%nsol(iiq) = ntmp
    rismt%qsol(iiq) = qv * ntmp
  END DO
  !
  IF (rismt%nsite > 0) THEN
    CALL mp_sum(rismt%nsol, rismt%mp_site%intra_sitg_comm)
    CALL mp_sum(rismt%qsol, rismt%mp_site%intra_sitg_comm)
  END IF
  !
  ! ... make qtot
  rismt%qtot = 0.0_DP
  DO iq = rismt%mp_site%isite_start, rismt%mp_site%isite_end
    iiq = iq - rismt%mp_site%isite_start + 1
    rismt%qtot = rismt%qtot + rismt%qsol(iiq)
  END DO
  !
  CALL mp_sum(rismt%qtot, rismt%mp_site%inter_sitg_comm)
  !
  ! ... make rhotot, rhor
  rhotot = 0.0_DP
  IF (rismt%dfft%nnr > 0) THEN
    rhor = 0.0_DP
  END IF
  !
  DO iq = rismt%mp_site%isite_start, rismt%mp_site%isite_end
    iiq    = iq - rismt%mp_site%isite_start + 1
    iv     = iuniq_to_isite(1, iq)
    nv     = iuniq_to_nsite(iq)
    isolV  = isite_to_isolV(iv)
    iatom  = isite_to_iatom(iv)
    rhov   = DBLE(nv) * solVs(isolV)%density
    qv     = solVs(isolV)%charge(iatom)
    rhotot = rhotot + qv * rhov
    IF (rismt%dfft%nnr > 0) THEN
      rhor(1:rismt%dfft%nnr) = rhor(1:rismt%dfft%nnr) &
      & + qv * rhov * rismt%gr(1:rismt%dfft%nnr, iiq)
    END IF
  END DO
  !
  CALL mp_sum(rhotot, rismt%mp_site%inter_sitg_comm)
  IF (rismt%dfft%nnr > 0) THEN
    CALL mp_sum(rhor, rismt%mp_site%inter_sitg_comm)
  END IF
  !
  ! ... make rhog
!$omp parallel do default(shared) private(ir)
  DO ir = 1, rismt%dfft%nnr
    aux(ir) = CMPLX(rhor(ir), 0.0_DP, kind=DP)
  END DO
!$omp end parallel do
  !
  IF (rismt%dfft%nnr > 0) THEN
    CALL fwfft('Rho', aux, rismt%dfft)
  END IF
  !
!$omp parallel do default(shared) private(ig)
  DO ig = 1, rismt%gvec%ngm
    rismt%rhog(ig) = aux(rismt%dfft%nl(ig))
  END DO
!$omp end parallel do
  !
  ! ... renormalize rhog
  IF (rismt%gvec%gstart > 1) THEN
    rismt%rhog(1) = CMPLX(rhotot, 0.0_DP, kind=DP)
  END IF
  !
  WRITE(stdout, '(/,5X,"solvent charge ",F10.5, &
                  & ", renormalised to ",F10.5)') rismt%qtot, rhotot * omega
  !
  ! ... make vpot
  fac = e2 * fpi / tpiba2
  IF (rismt%gvec%ngm > 0) THEN
    rismt%vpot = CMPLX(0.0_DP, 0.0_DP, kind=DP)
  END IF
  !
!$omp parallel do default(shared) private(ig, gg0)
  DO ig = rismt%gvec%gstart, rismt%gvec%ngm
    gg0 = rismt%gvec%gg(ig)
    rismt%vpot(ig) = fac * rismt%rhog(ig) / gg0
  END DO
!$omp end parallel do
  !
  ! ... make esol
  rismt%esol = 0.0_DP
  !
  DO iq = rismt%mp_site%isite_start, rismt%mp_site%isite_end
    iiq = iq - rismt%mp_site%isite_start + 1
    rismt%esol = rismt%esol + rismt%usol(iiq)
  END DO
  !
  CALL mp_sum(rismt%esol, rismt%mp_site%inter_sitg_comm)
  !
  ! ... deallocate working memory
  IF (rismt%dfft%nnr > 0) THEN
    DEALLOCATE(rhor)
    DEALLOCATE(aux)
  END IF
  !
  ! ... normally done
  ierr = IERR_RISM_NULL
  !
END SUBROUTINE solvation_3drism
