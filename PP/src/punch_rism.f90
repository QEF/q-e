!
! Copyright (C) 2016 National Institute of Advanced Industrial Science and Technology (AIST)
! [ This code is written by Satomichi Nishihara. ]
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------------------------------------------
SUBROUTINE punch_rism(filplot)
  !--------------------------------------------------------------------------
  !
  ! ... This subroutine writes data of 3D-RISM or Laue-RISM on file filplot.
  ! ... The data are following:
  ! ...   1) charge of solvent
  ! ...   2) potential of solvent
  ! ...   3) potential of solute
  ! ...   4) total potential
  ! ...   5) solvent-atomic densities
  !
  USE cell_base,      ONLY : at, celldm, ibrav, alat, omega
  USE control_flags,  ONLY : gamma_only
  USE fft_interfaces, ONLY : invfft
  USE io_global,      ONLY : stdout, ionode
  USE ions_base,      ONLY : nat, ntyp => nsp, ityp, tau, zv, atm
  USE kinds,          ONLY : DP
  USE lauefft,        ONLY : gather_lauefft_to_real
  USE mp,             ONLY : mp_rank, mp_barrier, mp_get, mp_sum
  USE rism,           ONLY : ITYPE_3DRISM, ITYPE_LAUERISM
  USE rism3d_facade,  ONLY : lrism3d, rism3t, rism3d_is_laue
  USE run_info,       ONLY : title
  USE scatter_mod,    ONLY : gather_grid
  USE solvmol,        ONLY : solVs, get_nuniq_in_solVs, &
                           & iuniq_to_isite, iuniq_to_nsite, isite_to_isolV, isite_to_iatom
  USE parallel_include
  !
  IMPLICIT NONE
  !
  CHARACTER(LEN=*) :: filplot
  !
  LOGICAL                        :: laue
  INTEGER                        :: nq
  INTEGER                        :: iq, iiq
  INTEGER                        :: iv, nv
  INTEGER                        :: isolV
  INTEGER                        :: iatom
  INTEGER                        :: n1, nx1
  INTEGER                        :: n2, nx2
  INTEGER                        :: n3, nx3
  INTEGER                        :: nz, nxz
  INTEGER                        :: nfft
  INTEGER                        :: ndata
  REAL(DP),          ALLOCATABLE :: raux(:,:)
  CHARACTER(LEN=12), ALLOCATABLE :: sdata(:)
  !
  ! ... check conditions
  IF (LEN_TRIM(filplot) < 1) THEN
    CALL errore('punch_rism', 'file name is empty', 1)
  END IF
  !
  IF (.NOT. lrism3d) THEN
    CALL errore('punch_rism', 'there are not data of 3D-RISM', 1)
  END IF
  !
  ! ... set Laue
  laue = rism3d_is_laue()
  !
  ! ... set number of data
  nq    = get_nuniq_in_solVs()
  ndata = 4 + nq
  !
  ! ... set FFT dimensions
  IF (.NOT. laue) THEN
    n1    = rism3t%dfft%nr1
    n2    = rism3t%dfft%nr2
    n3    = rism3t%dfft%nr3
    nx1   = rism3t%dfft%nr1x
    nx2   = rism3t%dfft%nr2x
    nx3   = rism3t%dfft%nr3x
    nz    = n3
    nxz   = nx3
    nfft  = nx1 * nx2 * nx3
  ELSE
    n1    = rism3t%lfft%dfft%nr1
    n2    = rism3t%lfft%dfft%nr2
    n3    = rism3t%lfft%dfft%nr3
    nx1   = rism3t%lfft%dfft%nr1x
    nx2   = rism3t%lfft%dfft%nr2x
    nx3   = rism3t%lfft%dfft%nr3x
    nz    = rism3t%lfft%nrz
    nxz   = rism3t%lfft%nrzx
    nfft  = nx1 * nx2 * nxz
  END IF
  !
  ! ... allocate working memory
  IF (ionode) THEN
    ALLOCATE(raux(nfft, ndata))
    ALLOCATE(sdata(ndata))
    raux  = 0.0_DP
    sdata = ''
  END IF
  !
  ! ... set name of data
  IF (ionode) THEN
    sdata(1) = 'chg'
    sdata(2) = 'v_solv'
    sdata(3) = 'v_solu'
    sdata(4) = 'v_tot'
    DO iq = 1, nq
      iv    = iuniq_to_isite(1, iq)
      isolV = isite_to_isolV(iv)
      iatom = isite_to_iatom(iv)
      sdata(4 + iq) = TRIM(ADJUSTL(solVs(isolV)%aname(iatom)))
    END DO
  END IF
  !
  ! ... set data to plot
  IF (.NOT. laue) THEN
    ! ... for 3D-RISM
    WRITE(stdout, '(/5x,"Punching data of 3D-RISM")')
    CALL gather_3drism()
    !
  ELSE
    ! ... for Laue-RISM
    WRITE(stdout, '(/5x,"Punching data of Laue-RISM")')
    CALL gather_lauerism()
    CALL define_expandcell()
    !
  END IF
  !
  ! ... plot data
  IF (ionode) THEN
    CALL plot_rism(TRIM(ADJUSTL(filplot)), title, ndata, nx1, nx2, nxz, n1, n2, nz, &
                 & nat, ntyp, ibrav, celldm, at, rism3t%gvec%gcutm, laue, &
                 & sdata, atm, ityp, zv, tau, raux, +1)
  END IF
  !
  ! ... deallocate working memory
  IF (ionode) THEN
    DEALLOCATE(raux)
    DEALLOCATE(sdata)
  END IF
  !
CONTAINS
  !
  !--------------------------------------------------------------------------
  SUBROUTINE gather_3drism()
    !--------------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    INTEGER                  :: ig
    INTEGER                  :: ir
    INTEGER                  :: fft_rank
    INTEGER                  :: my_group_id
    INTEGER                  :: io_group_id
    INTEGER                  :: owner_group_id
    REAL(DP)                 :: rhov
    REAL(DP),    ALLOCATABLE :: rhor(:)
    REAL(DP),    ALLOCATABLE :: vpor(:)
    REAL(DP),    ALLOCATABLE :: rtmp(:)
    COMPLEX(DP), ALLOCATABLE :: aux(:)
    !
    ! ... check condition
    IF (rism3t%itype /= ITYPE_3DRISM) THEN
      CALL errore('punch_rism', 'incorrect type of RISM', 1)
    END IF
    !
    IF (rism3t%mp_site%nsite < nq) THEN
      CALL errore('punch_rism', 'nsite < nq', 1)
    END IF
    !
    IF (rism3t%nr < rism3t%dfft%nnr) THEN
      CALL errore('punch_rism', 'nr < nnr', 1)
    END IF
    !
    IF (rism3t%ng < rism3t%gvec%ngm) THEN
      CALL errore('punch_rism', 'ng < ngmt', 1)
    END IF
    !
    ! ... allocate working memory
    IF (rism3t%dfft%nnr > 0) THEN
      ALLOCATE(aux(rism3t%dfft%nnr))
      ALLOCATE(rhor(rism3t%dfft%nnr))
      ALLOCATE(vpor(rism3t%dfft%nnr))
    END IF
    ALLOCATE(rtmp(nfft))
    !
    ! ... set MPI-ranks
    my_group_id = mp_rank(rism3t%mp_site%inter_sitg_comm)
    !
    io_group_id = 0
    IF (ionode) THEN
      io_group_id = my_group_id
    END IF
    CALL mp_sum(io_group_id, rism3t%mp_site%intra_sitg_comm)
    CALL mp_sum(io_group_id, rism3t%mp_site%inter_sitg_comm)
    !
    fft_rank = mp_rank(rism3t%dfft%comm)
    !
    ! ... (1) gather charge of solvent
    IF (my_group_id == io_group_id) THEN
      IF (rism3t%dfft%nnr > 0) THEN
        aux = CMPLX(0.0_DP, 0.0_DP, kind=DP)
      END IF
      DO ig = 1, rism3t%gvec%ngm
        aux(rism3t%dfft%nl(ig)) = rism3t%rhog(ig)
      END DO
      IF (gamma_only) THEN
        DO ig = rism3t%gvec%gstart, rism3t%gvec%ngm
          aux(rism3t%dfft%nlm(ig)) = CONJG(aux(rism3t%dfft%nl(ig)))
        END DO
      END IF
      !
      IF (rism3t%dfft%nnr > 0) THEN
        CALL invfft('Rho', aux, rism3t%dfft)
      END IF
      !
      rhor = 0.0_DP
      DO ir = 1, rism3t%dfft%nnr
        rhor(ir) = DBLE(aux(ir))
      END DO
      !
#if defined(__MPI)
      rtmp = 0.0_DP
      CALL gather_grid(rism3t%dfft, rhor, rtmp)
#else
      rtmp = rhor
#endif
      IF (ionode) THEN
        raux(:, 1) = rtmp(:)
      END IF
    END IF
    !
    ! ... (2) gather potential of solvent
    IF (my_group_id == io_group_id) THEN
      IF (rism3t%dfft%nnr > 0) THEN
        aux = CMPLX(0.0_DP, 0.0_DP, kind=DP)
      END IF
      DO ig = 1, rism3t%gvec%ngm
        aux(rism3t%dfft%nl(ig)) = rism3t%vpot(ig)
      END DO
      IF (gamma_only) THEN
        DO ig = rism3t%gvec%gstart, rism3t%gvec%ngm
          aux(rism3t%dfft%nlm(ig)) = CONJG(aux(rism3t%dfft%nl(ig)))
        END DO
      END IF
      !
      IF (rism3t%dfft%nnr > 0) THEN
        CALL invfft('Rho', aux, rism3t%dfft)
      END IF
      !
      vpor = 0.0_DP
      DO ir = 1, rism3t%dfft%nnr
        vpor(ir) = -DBLE(aux(ir))
      END DO
      !
#if defined(__MPI)
      rtmp = 0.0_DP
      CALL gather_grid(rism3t%dfft, vpor, rtmp)
#else
      rtmp = vpor
#endif
      IF (ionode) THEN
        raux(:, 2) = rtmp(:)
      END IF
    END IF
    !
    ! ... (3) gather potential of solute
    IF (my_group_id == io_group_id) THEN
      vpor = 0.0_DP
      DO ir = 1, rism3t%dfft%nnr
        vpor(ir) = -(rism3t%vsr(ir) + rism3t%vlr(ir))
      END DO
      !
#if defined(__MPI)
      rtmp = 0.0_DP
      CALL gather_grid(rism3t%dfft, vpor, rtmp)
#else
      rtmp = vpor
#endif
      IF (ionode) THEN
        raux(:, 3) = rtmp(:)
      END IF
    END IF
    !
    ! ... (4) gather total potential
    IF (ionode) THEN
      raux(:, 4) = raux(:, 2) + raux(:, 3)
    END IF
    !
    ! ... (5) gather solvent-atomic densities
    DO iq = 1, nq
      !
      IF (rism3t%mp_site%isite_start <= iq .AND. iq <= rism3t%mp_site%isite_end) THEN
        owner_group_id = my_group_id
        iiq = iq - rism3t%mp_site%isite_start + 1
#if defined(__MPI)
        rtmp = 0.0_DP
        CALL gather_grid(rism3t%dfft, rism3t%gr(:, iiq), rtmp)
#else
        rtmp = rism3t%gr(1:rism3t%dfft%nnr, iiq)
#endif
      ELSE
        owner_group_id = 0
        IF (fft_rank == rism3t%dfft%root) THEN
          rtmp = 0.0_DP
        END IF
      END IF
      !
      CALL mp_sum(owner_group_id, rism3t%mp_site%inter_sitg_comm)
      !
      IF (fft_rank == rism3t%dfft%root) THEN
        CALL mp_get(rtmp, rtmp, my_group_id, io_group_id, &
                  & owner_group_id, iq, rism3t%mp_site%inter_sitg_comm)
      END IF
      !
      IF (ionode) THEN
        iv    = iuniq_to_isite(1, iq)
        nv    = iuniq_to_nsite(iq)
        isolV = isite_to_isolV(iv)
        rhov  = DBLE(nv) * solVs(isolV)%density
        raux(:, 4 + iq) = rhov * rtmp(:)
      END IF
      !
      CALL mp_barrier(rism3t%intra_comm)
      !
    END DO
    !
    ! ... deallocate working memory
    IF (rism3t%dfft%nnr > 0) THEN
      DEALLOCATE(aux)
      DEALLOCATE(rhor)
      DEALLOCATE(vpor)
    END IF
    DEALLOCATE(rtmp)
    !
  END SUBROUTINE gather_3drism
  !
  !--------------------------------------------------------------------------
  SUBROUTINE gather_lauerism()
    !--------------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    INTEGER                  :: ista
    INTEGER                  :: iend
    INTEGER                  :: fft_rank
    INTEGER                  :: my_group_id
    INTEGER                  :: io_group_id
    INTEGER                  :: owner_group_id
    REAL(DP)                 :: rhovr
    REAL(DP)                 :: rhovl
    REAL(DP),    ALLOCATABLE :: rhor(:)
    REAL(DP),    ALLOCATABLE :: vpor(:)
    REAL(DP),    ALLOCATABLE :: rtmp(:)
    REAL(DP),    ALLOCATABLE :: stmp(:)
    COMPLEX(DP), ALLOCATABLE :: ltmp(:)
    !
    ! ... check condition
    IF (rism3t%itype /= ITYPE_LAUERISM) THEN
      CALL errore('punch_rism', 'incorrect type of RISM', 1)
    END IF
    !
    IF (rism3t%mp_site%nsite < nq) THEN
      CALL errore('punch_rism', 'nsite < nq', 1)
    END IF
    !
    IF (rism3t%nr < rism3t%lfft%dfft%nnr) THEN
      CALL errore('punch_rism', 'nr < nnr', 1)
    END IF
    !
    IF (rism3t%nrzl < rism3t%lfft%nrz) THEN
      CALL errore('punch_rism', 'nrzl < nrz', 1)
    END IF
    !
    ! ... allocate working memory
    IF (fft_rank == rism3t%lfft%dfft%root) THEN
      ALLOCATE(rtmp(nfft))
    ELSE
      ALLOCATE(rtmp(1))
    END IF
    ALLOCATE(stmp(nx1 * nx2 * nx3))
    ALLOCATE(ltmp(rism3t%nrzl * rism3t%ngxy))
    !
    ! ... set MPI-ranks
    my_group_id = mp_rank(rism3t%mp_site%inter_sitg_comm)
    !
    io_group_id = 0
    IF (ionode) THEN
      io_group_id = my_group_id
    END IF
    CALL mp_sum(io_group_id, rism3t%mp_site%intra_sitg_comm)
    CALL mp_sum(io_group_id, rism3t%mp_site%inter_sitg_comm)
    !
    fft_rank = mp_rank(rism3t%lfft%dfft%comm)
    !
    ! ... (1) gather charge of solvent
    IF (my_group_id == io_group_id) THEN
      rtmp = 0.0_DP
      CALL gather_lauefft_to_real(rism3t%lfft, rism3t%rhog, rism3t%nrzl, rtmp)
      IF (ionode) THEN
        raux(:, 1) = rtmp(:)
      END IF
    END IF
    !
    ! ... (2) gather potential of solvent
    IF (my_group_id == io_group_id) THEN
      rtmp = 0.0_DP
      CALL gather_lauefft_to_real(rism3t%lfft, rism3t%vpot, rism3t%nrzl, rtmp)
      IF (ionode) THEN
        raux(:, 2) = -rtmp(:)
      END IF
    END IF
    !
    ! ... (3) gather potential of solute
    IF (my_group_id == io_group_id) THEN
#if defined(__MPI)
      stmp = 0.0_DP
      CALL gather_grid(rism3t%lfft%dfft, rism3t%vsr, stmp)
#else
      stmp = rism3t%vsr(1:rism3t%lfft%dfft%nnr)
#endif
      rtmp = 0.0_DP
      CALL gather_lauefft_to_real(rism3t%lfft, rism3t%vlgz, rism3t%nrzl, rtmp)
      IF (fft_rank == rism3t%lfft%dfft%root) THEN
        CALL set_on_unit_cell(stmp, rtmp, .TRUE.)
      END IF
      IF (ionode) THEN
        raux(:, 3) = -rtmp(:)
      END IF
    END IF
    !
    ! ... (4) gather total potential
    IF (ionode) THEN
      raux(:, 4) = raux(:, 2) + raux(:, 3)
    END IF
    !
    ! ... (5) gather solvent-atomic densities
    DO iq = 1, nq
      !
      IF (rism3t%mp_site%isite_start <= iq .AND. iq <= rism3t%mp_site%isite_end) THEN
        owner_group_id = my_group_id
        iiq = iq - rism3t%mp_site%isite_start + 1
#if defined(__MPI)
        stmp = 0.0_DP
        CALL gather_grid(rism3t%lfft%dfft, rism3t%gr(:, iiq), stmp)
#else
        stmp = rism3t%gr(1:rism3t%lfft%dfft%nnr, iiq)
#endif
        rtmp = 0.0_DP
        ltmp = rism3t%hsgz(:, iiq) + rism3t%hlgz(:, iiq)
        CALL gather_lauefft_to_real(rism3t%lfft, ltmp, rism3t%nrzl, rtmp)
        IF (fft_rank == rism3t%lfft%dfft%root) THEN
          rtmp = rtmp + 1.0_DP
          CALL set_on_unit_cell(stmp, rtmp, .FALSE.)
        END IF
      ELSE
        owner_group_id = 0
        IF (fft_rank == rism3t%lfft%dfft%root) THEN
          rtmp = 0.0_DP
        END IF
      END IF
      !
      CALL mp_sum(owner_group_id, rism3t%mp_site%inter_sitg_comm)
      !
      IF (fft_rank == rism3t%lfft%dfft%root) THEN
        CALL mp_get(rtmp, rtmp, my_group_id, io_group_id, &
                  & owner_group_id, iq, rism3t%mp_site%inter_sitg_comm)
      END IF
      !
      IF (ionode) THEN
        iv    = iuniq_to_isite(1, iq)
        nv    = iuniq_to_nsite(iq)
        isolV = isite_to_isolV(iv)
        rhovr = DBLE(nv) * solVs(isolV)%density
        rhovl = DBLE(nv) * solVs(isolV)%subdensity
        !
        ista = 1
        iend = nx1 * nx2 * rism3t%lfft%izleft_end
        rtmp(ista:iend) = rhovl * rtmp(ista:iend)
        !
        ista = nx1 * nx2 * (rism3t%lfft%izright_start - 1) + 1
        iend = nx1 * nx2 * nz
        rtmp(ista:iend) = rhovr * rtmp(ista:iend)
        !
        raux(:, 4 + iq) = rtmp(:)
      END IF
      !
      CALL mp_barrier(rism3t%intra_comm)
      !
    END DO
    !
    ! ... deallocate working memory
    DEALLOCATE(rtmp)
    DEALLOCATE(stmp)
    DEALLOCATE(ltmp)
    !
  END SUBROUTINE gather_lauerism
  !
  !--------------------------------------------------------------------------
  SUBROUTINE set_on_unit_cell(rsrc, rdst, ladd)
    !--------------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    REAL(DP), INTENT(IN)    :: rsrc(:)  !dimension(nx1*nx2*nx3)
    REAL(DP), INTENT(INOUT) :: rdst(:)  !dimension(nx1*nx2*nxz)
    LOGICAL,  INTENT(IN)    :: ladd
    !
    INTEGER :: ir, iir
    INTEGER :: i1, i2, i3
    INTEGER :: iz, iz0
    !
    iz0 = rism3t%lfft%izcell_start - 1
    !
    DO i3 = 1, n3
      !
      IF (i3 > (n3 - (n3 / 2))) THEN
        iz = i3 - n3
      ELSE
        iz = i3
      END IF
      iz = iz + (n3 / 2)
      iz = iz + iz0
      !
      DO i2 = 1, n2
        IF (ladd) THEN
          DO i1 = 1, n1
            ir  = i1 + (i2 - 1) * nx1 + (i3 - 1) * nx1 * nx2
            iir = i1 + (i2 - 1) * nx1 + (iz - 1) * nx1 * nx2
            rdst(iir) = rdst(iir) + rsrc(ir)
          END DO
        ELSE
          DO i1 = 1, n1
            ir  = i1 + (i2 - 1) * nx1 + (i3 - 1) * nx1 * nx2
            iir = i1 + (i2 - 1) * nx1 + (iz - 1) * nx1 * nx2
            rdst(iir) = rsrc(ir)
          END DO
        END IF
      END DO
      !
    END DO
    !
  END SUBROUTINE set_on_unit_cell
  !
  !--------------------------------------------------------------------------
  SUBROUTINE define_expandcell()
    !--------------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    INTEGER  :: ia
    REAL(DP) :: zleft
    REAL(DP) :: zright
    !
    zleft  = rism3t%lfft%zleft
    zright = rism3t%lfft%zright
    !
    ! ... shift atomic positions
    DO ia = 1, nat
      tau(3, ia) = tau(3, ia) - zleft
    END DO
    !
    ! ... re-define cell
    ibrav    = 0
    celldm   = 0.0_DP
    at(3, 3) = zright - zleft
    at       = alat * at
    !
    CALL latgen(ibrav, celldm, at(1, 1), at(1, 2), at(1, 3), omega)
    alat = celldm(1)
    at   = at / alat
    !
  END SUBROUTINE define_expandcell
  !
END SUBROUTINE punch_rism
