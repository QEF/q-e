  !
  ! Copyright (C) 2016-2023 EPW-Collaboration
  ! Copyright (C) 2010-2016 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino
  ! Copyright (C) 2007-2009 Roxana Margine
  !
  ! This file is distributed under the terms of the GNU General Public
  ! License. See the file `LICENSE' in the root directory of the
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .
  !
  !----------------------------------------------------------------------
  MODULE supercond_aniso
  !----------------------------------------------------------------------
  !!
  !! This module contains all the subroutines linked with superconductivity using
  !! the isotropic or anisotropic Eliashberg formalism.
  !!
  IMPLICIT NONE
  !
  CONTAINS
    !
    !-----------------------------------------------------------------------
    SUBROUTINE eliashberg_aniso_iaxis
    !-----------------------------------------------------------------------
    !!
    !! This routine is the driver of the self-consistent cycle for the anisotropic
    !! Eliashberg equations on the imaginary-axis.
    !!
    !! SH: Modified to allow for fbw calculations (Nov 2021).
    !! HM: Modified to allow for fbw calculation with sparse-ir sampling (May 2024)
    !!
    USE kinds,         ONLY : DP
    USE control_flags, ONLY : iverbosity
    USE global_var,    ONLY : gtemp
    USE input,         ONLY : nqstep, nsiter, nstemp, broyden_beta, broyden_ndim, &
                              limag, lpade, lacon, fsthick, imag_read, npade, &
                              fbw, positive_matsu, icoulomb, filirobj, gridsamp, &
                              positive_matsu
    USE supercond_common, ONLY : nsw, nsiw, lacon_fly, adelta, adeltap, adeltai, &
                              adeltaip, aznormi, aznormip, ashifti, ashiftip, &
                              nkfs, nqfs, nbndfs, ekfs, ef0, agap, nbnd_cl, &
                              adeltai_cl, adeltaip_cl, nkstot_cl, siz_ir, &
                              siz_ir_cl
    USE supercond,     ONLY : dos_quasiparticle, gen_freqgrid_iaxis, &
                              eliashberg_grid
    USE ep_constants,  ONLY : kelvin2eV, ci, zero, czero
    USE io_global,     ONLY : stdout, ionode_id
    USE mp_world,      ONLY : mpime
    USE mp_global,     ONLY : inter_pool_comm
    USE mp,            ONLY : mp_barrier, mp_sum
    USE io_supercond,  ONLY : eliashberg_write_iaxis, eliashberg_read_aniso_iaxis, &
                              eliashberg_write_raxis
    USE utilities,     ONLY : mix_wrap
    USE low_lvl,       ONLY : mem_size_eliashberg
    USE printing,      ONLY : prtheader_supercond
    USE parallelism,   ONLY : fkbounds, para_bounds
    USE io_var,        ONLY : iunirobj
    USE sparse_ir,     ONLY : IR, finalize_ir, set_beta
    USE io_sparse_ir,  ONLY : read_ir_epw
    !
    IMPLICIT NONE
    !
    ! Local variables
    CHARACTER(LEN = 256) :: cname
    !! character in file name
    LOGICAL :: conv
    !! True if calculation is converged
    INTEGER :: itemp
    !! Counter on temperature index
    INTEGER :: iter
    !! Counter on iteration steps
    INTEGER :: N
    !! Maximum nr. frequency points in Pade approx
    INTEGER :: ik
    !! Counter on k-points
    INTEGER :: ibnd
    !! Counter on bands
    INTEGER :: lower_bnd, upper_bnd
    !! Lower/upper bound index after k paral
    INTEGER :: nks
    !! Number of k points per pool
    INTEGER :: istate
    !! Counter on states, used when icoulomb > 0
    INTEGER :: is_start, is_stop
    !! Lower/upper bound index after paral for states
    INTEGER :: num_states
    !! Number of states per pool: num_states = is_stop - is_start + 1
    INTEGER :: narray(2)
    !! integers to determine the size of arrays of FFT
    TYPE(IR) :: ir_obj
    !! which contains ir objects such as basis functions
    INTEGER(8) :: imelt
    !! Counter memory
    INTEGER :: ierr
    !! Error status
    !
    REAL(KIND = DP) :: tcpu
    !! cpu time
    REAL(KIND = DP) :: beta
    !! Invese temperature beta = 1 / gtemp.
    REAL(KIND = DP), EXTERNAL :: get_clock
    !! get the time spent
    REAL(KIND = DP), ALLOCATABLE :: rdeltain(:), rdeltaout(:)
    !! Temporary variables for mix_broyden in analytic continuation
    REAL(KIND = DP), ALLOCATABLE :: cdeltain(:), cdeltaout(:)
    !! Temporary variables for mix_broyden in analytic continuation
    REAL(KIND = DP), ALLOCATABLE :: dv1(:, :, :, :), df1(:, :, :, :)
    !! Temporary variables for mix_broyden
    REAL(KIND = DP), ALLOCATABLE :: dv2(:, :, :, :), df2(:, :, :, :)
    !! Temporary variables for mix_broyden
    REAL(KIND = DP), ALLOCATABLE :: dv3(:, :, :, :), df3(:, :, :, :)
    !! Temporary variables for mix_broyden
    REAL(KIND = DP), ALLOCATABLE :: dv4(:, :, :), df4(:, :, :)
    !! Temporary variables for mix_broyden
    CHARACTER(LEN=100) filename
    !! the file name of ir object file
    !
    CALL start_clock('aniso_iaxis')
    !
    CALL eliashberg_grid()
    !
    CALL fkbounds(nkfs, lower_bnd, upper_bnd)
    !
    nks = upper_bnd - lower_bnd + 1
    !
    IF (fbw .AND. (gridsamp == 2) .AND. (icoulomb > 0)) THEN
      !
      ! HM: Because nkstot_cl can be quite small,
      !     we do a combined parallelization on band index and k point index.
      CALL para_bounds(is_start, is_stop, nbnd_cl * nkstot_cl)
      num_states = is_stop - is_start + 1
      !
    ENDIF
    !
    DO itemp = 1, nstemp ! loop over temperature
      !
      ! HM: change the order of subroutines
      !     because nsiw(itemp) is not defined until calling gen_freqgrid_iaxis
      !     if using sparse-ir method.
      CALL start_clock('iaxis_imag')
      !
      IF (fbw .AND. (gridsamp == 2)) THEN
        !
        beta = 1.0E0_DP / gtemp(itemp)
        !
        IF (itemp == 1) THEN
          !
          WRITE(stdout, '(/5x, a/)') 'Start reading ir object file'
          IF (mpime == ionode_id) THEN
            filename = TRIM(filirobj)
            OPEN(UNIT = iunirobj, FILE=filename, STATUS='old',IOSTAT=ierr)
            IF (ierr /= 0) CALL errore('eliashberg_aniso_iaxis', 'Error while opening the ir object file', 1)
          ENDIF
          ir_obj = read_ir_epw(iunirobj, beta, positive_matsu)
          IF (mpime == ionode_id) THEN
            CLOSE(iunirobj)
          ENDIF
          WRITE(stdout, '(/5x, a/)') 'Finish reading ir object file'
        ELSE
          !
          ! HM: In the current implementation, the same IR objects is used for all temperatures.
          !     If you want to use different IR object files depending on the temperature,
          !     please declare 'imatches' and make changes in the following IF statement.
          !
          !IF (.NOT. imatches(TRIM(filirobj), TRIM(filirobj_old))) THEN
          IF (1 == 0) THEN
            WRITE(stdout, '(/5x, a/)') 'Start reading ir object file'
            IF (mpime == ionode_id) THEN
              filename = TRIM(filirobj)
              OPEN(UNIT = iunirobj, FILE=filename, STATUS='old',IOSTAT=ierr)
              IF (ierr /= 0) CALL errore('eliashberg_aniso_iaxis', 'Error while opening the ir object file', 1)
            ENDIF
            ir_obj = read_ir_epw(iunirobj, beta)
            IF (mpime == ionode_id) THEN
              CLOSE(iunirobj)
            ENDIF
          ELSE
            ! HM: If the same IR objects are used,
            !     update the temperature using set_beta.
            CALL set_beta(ir_obj, beta)
          ENDIF !
        ENDIF ! itemp == 1
        !
        CALL gen_freqgrid_iaxis(itemp, ir_obj)
        !
      ELSE
        CALL gen_freqgrid_iaxis(itemp)
      ENDIF
      !
      CALL prtheader_supercond(itemp, 1)
      !
      IF (itemp == 1) THEN
        ALLOCATE(agap(nbndfs, nkfs), STAT = ierr)
        IF (ierr /= 0) CALL errore('eliashberg_aniso_iaxis', 'Error allocating agap', 1)
        agap(:, :) = zero
      ENDIF
      !
      IF ((limag .AND. .NOT. imag_read) .OR. (limag .AND. imag_read .AND. itemp /= 1)) THEN
        !
        ! get the size of required memory for df1, dv1
        imelt = 2 * nbndfs * nks * broyden_ndim * nsiw(itemp)
        CALL mem_size_eliashberg(2, imelt)
        !
        ! RM - df1, dv1 are defined per k-points per pool
        ALLOCATE(df1(nbndfs, lower_bnd:upper_bnd, nsiw(itemp), broyden_ndim), STAT = ierr)
        IF (ierr /= 0) CALL errore('eliashberg_aniso_iaxis', 'Error allocating df1', 1)
        ALLOCATE(dv1(nbndfs, lower_bnd:upper_bnd, nsiw(itemp), broyden_ndim), STAT = ierr)
        IF (ierr /= 0) CALL errore('eliashberg_aniso_iaxis', 'Error allocating dv1', 1)
        df1(:, :, :, :) = zero
        dv1(:, :, :, :) = zero
        !
        IF (fbw) THEN
          !
          ! get the size of required memory for df2, dv2, df3, dv3
          imelt = 4 * nbndfs * nks * broyden_ndim * nsiw(itemp)
          CALL mem_size_eliashberg(2, imelt)
          !
          ! RM - df2, dv2, df3, dv3 are defined per k-points per pool
          ALLOCATE(df2(nbndfs, lower_bnd:upper_bnd, nsiw(itemp), broyden_ndim), STAT = ierr)
          IF (ierr /= 0) CALL errore('eliashberg_aniso_iaxis', 'Error allocating df2', 1)
          ALLOCATE(dv2(nbndfs, lower_bnd:upper_bnd, nsiw(itemp), broyden_ndim), STAT = ierr)
          IF (ierr /= 0) CALL errore('eliashberg_aniso_iaxis', 'Error allocating dv2', 1)
          ALLOCATE(df3(nbndfs, lower_bnd:upper_bnd, nsiw(itemp), broyden_ndim), STAT = ierr)
          IF (ierr /= 0) CALL errore('eliashberg_aniso_iaxis', 'Error allocating df3', 1)
          ALLOCATE(dv3(nbndfs, lower_bnd:upper_bnd, nsiw(itemp), broyden_ndim), STAT = ierr)
          IF (ierr /= 0) CALL errore('eliashberg_aniso_iaxis', 'Error allocating dv3', 1)
          df2(:, :, :, :) = zero
          dv2(:, :, :, :) = zero
          df3(:, :, :, :) = zero
          dv3(:, :, :, :) = zero
          !
          IF ((gridsamp == 2) .AND. (icoulomb > 0)) THEN
            ! get the size of required memory for df4, dv4
            imelt = 2 * num_states * broyden_ndim * nsiw(itemp)
            CALL mem_size_eliashberg(2, imelt)
            !
            ALLOCATE(df4(is_start:is_stop, nsiw(itemp), broyden_ndim), STAT = ierr)
            IF (ierr /= 0) CALL errore('eliashberg_aniso_iaxis', 'Error allocating df4', 1)
            ALLOCATE(dv4(is_start:is_stop, nsiw(itemp), broyden_ndim), STAT = ierr)
            IF (ierr /= 0) CALL errore('eliashberg_aniso_iaxis', 'Error allocating dv4', 1)
            df4(:, :, :) = zero
            dv4(:, :, :) = zero
          ENDIF
          !
        ENDIF ! fbw
        !
        iter = 1
        conv = .FALSE.
        DO WHILE (.NOT. conv .AND. iter <= nsiter)
          IF (fbw .AND. (gridsamp == 2)) THEN
            CALL sum_eliashberg_aniso_iaxis_wrapper(itemp, iter, conv, ir_obj)
          ELSE
            CALL sum_eliashberg_aniso_iaxis_wrapper(itemp, iter, conv)
          ENDIF
          IF (fbw) THEN
            DO ik = lower_bnd, upper_bnd
              DO ibnd = 1, nbndfs
                IF (ABS(ekfs(ibnd, ik) - ef0) < fsthick) THEN
                  CALL mix_wrap(nsiw(itemp), adeltai(:, ibnd, ik), adeltaip(:, ibnd, ik), &
                                broyden_beta, iter, broyden_ndim, conv, &
                                df1(ibnd, ik, :, :), dv1(ibnd, ik, :, :))
                  CALL mix_wrap(nsiw(itemp), aznormi(:, ibnd, ik), aznormip(:, ibnd, ik), &
                                broyden_beta, iter, broyden_ndim, &
                                conv, df2(ibnd, ik, :, :), dv2(ibnd, ik, :, :))
                  CALL mix_wrap(nsiw(itemp), ashifti(:, ibnd, ik), ashiftip(:, ibnd, ik), &
                                broyden_beta, iter, broyden_ndim, conv, &
                                df3(ibnd, ik, :, :), dv3(ibnd, ik, :, :))
                ENDIF
              ENDDO
            ENDDO
            ! Make everything 0 except the range of k-points we are working on
            aznormip(:, :, 1:lower_bnd - 1) = zero
            ashiftip(:, :, 1:lower_bnd - 1) = zero
            adeltaip(:, :, 1:lower_bnd - 1) = zero
            aznormip(:, :, lower_bnd + nks:nkfs) = zero
            ashiftip(:, :, lower_bnd + nks:nkfs) = zero
            adeltaip(:, :, lower_bnd + nks:nkfs) = zero
            !
            ! collect contributions from all pools
            CALL mp_sum(aznormip, inter_pool_comm)
            CALL mp_sum(ashiftip, inter_pool_comm)
            CALL mp_sum(adeltaip, inter_pool_comm)
            IF ((gridsamp == 2) .AND. (icoulomb > 0)) THEN
              IF (num_states > 0) THEN
                DO istate = is_start, is_stop
                  ik = (istate - 1) / nbnd_cl + 1
                  ibnd = MOD(istate - 1, nbnd_cl) + 1
                  CALL mix_wrap(nsiw(itemp), adeltai_cl(:, istate), adeltaip_cl(:, ibnd, ik), &
                                broyden_beta, iter, broyden_ndim, conv, &
                                df4(istate, :, :), dv4(istate, :, :))
                ENDDO
              ENDIF
              ! Make everything 0 except the range of states we are working on
              IF (num_states < 1) THEN
                adeltaip_cl(:, : ,:) = zero
              ELSE
                DO istate = 1, is_start - 1
                  ik = (istate - 1) / nbnd_cl + 1
                  ibnd = MOD(istate - 1, nbnd_cl) + 1
                  !
                  adeltaip_cl(:, ibnd, ik) = zero
                ENDDO
                DO istate = is_stop + 1, nbnd_cl * nkstot_cl
                  ik = (istate - 1) / nbnd_cl + 1
                  ibnd = MOD(istate - 1, nbnd_cl) + 1
                  !
                  adeltaip_cl(:, ibnd, ik) = zero
                ENDDO
              ENDIF
              !
              ! collect contributions from all pools
              CALL mp_sum(adeltaip_cl, inter_pool_comm)
              !
            ENDIF
            CALL mp_barrier(inter_pool_comm)
            !
          ELSE ! not fbw
            DO ik = lower_bnd, upper_bnd
              DO ibnd = 1, nbndfs
                IF (ABS(ekfs(ibnd, ik) - ef0) < fsthick) THEN
                  CALL mix_wrap(nsiw(itemp), adeltai(:, ibnd, ik), adeltaip(:, ibnd, ik), &
                                broyden_beta, iter, broyden_ndim, conv, &
                                df1(ibnd, ik, :, :), dv1(ibnd, ik, :, :))
                ENDIF
              ENDDO
            ENDDO
            ! Make everything 0 except the range of k-points we are working on
            adeltaip(:, :, 1:lower_bnd - 1) = zero
            adeltaip(:, :, lower_bnd + nks:nkfs) = zero
            !
            ! collect contributions from all pools
            CALL mp_sum(adeltaip, inter_pool_comm)
            CALL mp_barrier(inter_pool_comm)
          ENDIF ! fbw
          iter = iter + 1
        ENDDO ! iter
        !
        DEALLOCATE(df1, STAT = ierr)
        IF (ierr /= 0) CALL errore('eliashberg_aniso_iaxis', 'Error deallocating df1', 1)
        DEALLOCATE(dv1, STAT = ierr)
        IF (ierr /= 0) CALL errore('eliashberg_aniso_iaxis', 'Error deallocating dv1', 1)
        !
        ! remove memory allocated for df1, dv1
        imelt = 2 * nbndfs * nks * broyden_ndim * nsiw(itemp)
        CALL mem_size_eliashberg(2, -imelt)
        !
        IF (fbw) THEN
          DEALLOCATE(df2, STAT = ierr)
          IF (ierr /= 0) CALL errore('eliashberg_aniso_iaxis', 'Error deallocating df2', 1)
          DEALLOCATE(dv2, STAT = ierr)
          IF (ierr /= 0) CALL errore('eliashberg_aniso_iaxis', 'Error deallocating dv2', 1)
          DEALLOCATE(df3, STAT = ierr)
          IF (ierr /= 0) CALL errore('eliashberg_aniso_iaxis', 'Error deallocating df3', 1)
          DEALLOCATE(dv3, STAT = ierr)
          IF (ierr /= 0) CALL errore('eliashberg_aniso_iaxis', 'Error deallocating dv3', 1)
          !
          ! remove memory allocated for df2, dv2, df3, dv3
          imelt = 4 * nbndfs * nks * broyden_ndim * nsiw(itemp)
          CALL mem_size_eliashberg(2, -imelt)
          !
          IF ((gridsamp == 2) .AND. (icoulomb > 0)) THEN
            DEALLOCATE(df4, STAT = ierr)
            IF (ierr /= 0) CALL errore('eliashberg_aniso_iaxis', 'Error allocating df4', 1)
            DEALLOCATE(dv4, STAT = ierr)
            IF (ierr /= 0) CALL errore('eliashberg_aniso_iaxis', 'Error allocating dv4', 1)
            ! remove memory allocated for df4, dv4
            imelt = 2 * num_states * broyden_ndim * nsiw(itemp)
            CALL mem_size_eliashberg(2, -imelt)
          ENDIF
          !
        ENDIF ! fbw
        !
        IF (conv) THEN
          CALL deallocate_aniso_iaxis(1)
          !
          ! remove memory allocated for deltai, znormi, adeltaip
          imelt = (2 + nbndfs * nkfs) * nsiw(itemp)
          CALL mem_size_eliashberg(2, -imelt)
          !
          IF (fbw) THEN
            ! remove memory allocated for shifti, aznormip, ashiftip
            imelt = (1 + 2 * nbndfs * nkfs) * nsiw(itemp)
            CALL mem_size_eliashberg(2, -imelt)
          ENDIF
          !
          CALL eliashberg_write_iaxis(itemp)
          !
          CALL stop_clock('iaxis_imag')
          CALL print_clock('iaxis_imag')
          WRITE(stdout, '(a)') ' '
        ELSEIF (.NOT. conv .AND. (iter - 1) == nsiter) THEN
          CALL deallocate_aniso_iaxis(5)
          CALL deallocate_aniso()
          CALL stop_clock('iaxis_imag')
          CALL print_clock('iaxis_imag')
          WRITE(stdout, '(a)') ' '
          RETURN
        ENDIF
      ELSEIF (limag .AND. imag_read .AND. itemp == 1) THEN
        CALL eliashberg_read_aniso_iaxis(itemp)
      ENDIF ! limag
      !
      IF (lpade) THEN
        CALL prtheader_supercond(itemp, 2)
        CALL start_clock('raxis_pade')
        IF (fbw .AND. (.NOT. positive_matsu)) THEN
          N = npade * (nsiw(itemp) / 2) / 100
          IF (mod(N, 2) /= 0 ) N = N + 1
          IF (N > (nsiw(itemp) / 2)) N = N - 2
        ELSE
          N = npade * nsiw(itemp) / 100
          IF (mod(N, 2) /= 0 ) N = N + 1
          IF (N > nsiw(itemp)) N = N - 2
        ENDIF
        CALL pade_cont_aniso(itemp, N)
        !
        cname = 'pade'
        CALL eliashberg_write_raxis(itemp, cname)
        !
        CALL dos_quasiparticle(itemp)
        !
        CALL stop_clock('raxis_pade')
        CALL print_clock('raxis_pade')
        WRITE(stdout, '(a)') ' '
      ENDIF ! lpade
      !
      ! HP: lacon is not implemented for fbw runs
      IF (.NOT. fbw .AND. lacon) THEN
        CALL prtheader_supercond(itemp, 3)
        CALL start_clock('raxis_acon')
        !
        iter = 1
        conv = .FALSE.
        !
        ! get the size of required memory for rdeltain, cdeltain, cdeltaout, rdeltaout,
        ! df1, dv1, df2, dv2
        imelt = 4 * (1 + nbndfs * nks * broyden_ndim) * nsw
        CALL mem_size_eliashberg(2, imelt)
        !
        ALLOCATE(rdeltain(nsw), STAT = ierr)
        IF (ierr /= 0) CALL errore('eliashberg_aniso_iaxis', 'Error allocating rdeltain', 1)
        ALLOCATE(cdeltain(nsw), STAT = ierr)
        IF (ierr /= 0) CALL errore('eliashberg_aniso_iaxis', 'Error allocating cdeltain', 1)
        ALLOCATE(rdeltaout(nsw), STAT = ierr)
        IF (ierr /= 0) CALL errore('eliashberg_aniso_iaxis', 'Error allocating rdeltaout', 1)
        ALLOCATE(cdeltaout(nsw), STAT = ierr)
        IF (ierr /= 0) CALL errore('eliashberg_aniso_iaxis', 'Error allocating cdeltaout', 1)
        ! RM - df1, dv1, df2, dv2 are defined per k-points per pool
        ALLOCATE(df1(nbndfs, lower_bnd:upper_bnd, nsw, broyden_ndim), STAT = ierr)
        IF (ierr /= 0) CALL errore('eliashberg_aniso_iaxis', 'Error allocating df1', 1)
        ALLOCATE(dv1(nbndfs, lower_bnd:upper_bnd, nsw, broyden_ndim), STAT = ierr)
        IF (ierr /= 0) CALL errore('eliashberg_aniso_iaxis', 'Error allocating dv1', 1)
        ALLOCATE(df2(nbndfs, lower_bnd:upper_bnd, nsw, broyden_ndim), STAT = ierr)
        IF (ierr /= 0) CALL errore('eliashberg_aniso_iaxis', 'Error allocating df2', 1)
        ALLOCATE(dv2(nbndfs, lower_bnd:upper_bnd, nsw, broyden_ndim), STAT = ierr)
        IF (ierr /= 0) CALL errore('eliashberg_aniso_iaxis', 'Error allocating dv2', 1)
        rdeltain(:)  = zero
        cdeltain(:)  = zero
        rdeltaout(:) = zero
        cdeltaout(:) = zero
        df1(:, :, :, :) = zero
        dv1(:, :, :, :) = zero
        df2(:, :, :, :) = zero
        dv2(:, :, :, :) = zero
        !
        DO WHILE (.NOT. conv .AND. iter <= nsiter)
          CALL analytic_cont_aniso(itemp, iter, conv)
          DO ik = lower_bnd, upper_bnd
            DO ibnd = 1, nbndfs
              IF (ABS(ekfs(ibnd, ik) - ef0) < fsthick) THEN
                rdeltain(:)  = REAL(adeltap(:, ibnd, ik))
                cdeltain(:)  = AIMAG(adeltap(:, ibnd, ik))
                rdeltaout(:) = REAL(adelta(:, ibnd, ik))
                cdeltaout(:) = AIMAG(adelta(:, ibnd, ik))
                CALL mix_wrap(nsw, rdeltaout, rdeltain, broyden_beta, iter, broyden_ndim, &
                              conv, df1(ibnd, ik, :, :), dv1(ibnd, ik, :, :))
                CALL mix_wrap(nsw, cdeltaout, cdeltain, broyden_beta, iter, broyden_ndim, &
                              conv, df2(ibnd, ik, :, :), dv2(ibnd, ik, :, :))
                adeltap(:, ibnd, ik) = rdeltain(:) + ci * cdeltain(:)
              ENDIF
            ENDDO
          ENDDO
          ! Make everything 0 except the range of k-points we are working on
          adeltap(:, :, 1:lower_bnd - 1) = czero
          adeltap(:, :, lower_bnd + nks:nkfs) = czero
          !
          ! collect contributions from all pools
          CALL mp_sum(adeltap, inter_pool_comm)
          CALL mp_barrier(inter_pool_comm)
          iter = iter + 1
        ENDDO ! iter
        !
        DEALLOCATE(rdeltain, STAT = ierr)
        IF (ierr /= 0) CALL errore('eliashberg_aniso_iaxis', 'Error deallocating rdeltain', 1)
        DEALLOCATE(cdeltain, STAT = ierr)
        IF (ierr /= 0) CALL errore('eliashberg_aniso_iaxis', 'Error deallocating cdeltain', 1)
        DEALLOCATE(rdeltaout, STAT = ierr)
        IF (ierr /= 0) CALL errore('eliashberg_aniso_iaxis', 'Error deallocating rdeltaout', 1)
        DEALLOCATE(cdeltaout, STAT = ierr)
        IF (ierr /= 0) CALL errore('eliashberg_aniso_iaxis', 'Error deallocating cdeltaout', 1)
        DEALLOCATE(df1, STAT = ierr)
        IF (ierr /= 0) CALL errore('eliashberg_aniso_iaxis', 'Error deallocating df1', 1)
        DEALLOCATE(dv1, STAT = ierr)
        IF (ierr /= 0) CALL errore('eliashberg_aniso_iaxis', 'Error deallocating dv1', 1)
        DEALLOCATE(df2, STAT = ierr)
        IF (ierr /= 0) CALL errore('eliashberg_aniso_iaxis', 'Error deallocating df2', 1)
        DEALLOCATE(dv2, STAT = ierr)
        IF (ierr /= 0) CALL errore('eliashberg_aniso_iaxis', 'Error deallocating dv2', 1)
        !
        ! remove memory allocated for rdeltain, cdeltain, cdeltaout, rdeltaout, df1, dv1, df2, dv2
        imelt = 4 * (1 + nbndfs * nks * broyden_ndim) * nsw
        CALL mem_size_eliashberg(2, -imelt)
        !
        IF (conv) THEN
          IF (limag .AND. imag_read .AND. itemp == 1) THEN
            CALL deallocate_aniso_iaxis(2)
          ELSE
            CALL deallocate_aniso_iaxis(4)
          ENDIF
          !
          ! remove memory allocated for wsi, wsn, adeltai, aznormi
          imelt = (2 + 2 * nbndfs * nks) * nsiw(itemp)
          CALL mem_size_eliashberg(2, -imelt)
          !
          CALL deallocate_aniso_raxis(1)
          !
          ! remove memory allocated for adeltap, aznormp
          imelt = 2 * nbndfs * nkfs * nsw
          CALL mem_size_eliashberg(2, -imelt)
          !
          ! remove memory allocated for gp, gm, adsumi, azsumi
          imelt = 2 * (nqstep + nks * nbndfs) * nsw
          CALL mem_size_eliashberg(2, -imelt)
          !
          IF (.NOT. lacon_fly) THEN
            ! remove memory allocated for a2fij
            imelt = nks * MAXVAL(nqfs(:)) * nbndfs**2 * nqstep
            CALL mem_size_eliashberg(2, -imelt)
          ENDIF
          !
          cname = 'acon'
          CALL eliashberg_write_raxis(itemp, cname)
          !
          CALL dos_quasiparticle(itemp) ! - change to work on pool
          !
          CALL stop_clock('raxis_acon')
          CALL print_clock('raxis_acon')
          WRITE(stdout, '(a)') ' '
        ELSEIF (.NOT. conv .AND. (iter - 1) == nsiter) THEN
          IF (limag .AND. imag_read .AND. itemp == 1) THEN
            CALL deallocate_aniso_iaxis(2)
          ELSE
            CALL deallocate_aniso_iaxis(4)
          ENDIF
          !
          CALL deallocate_aniso_raxis(3)
          CALL deallocate_aniso()
          CALL stop_clock('raxis_acon')
          CALL print_clock('raxis_acon')
          WRITE(stdout, '(a)') ' '
          RETURN
        ENDIF
      ENDIF ! lacon
      !
      IF (.NOT. lacon) THEN
        IF (limag .AND. imag_read .AND. itemp == 1) THEN
          CALL deallocate_aniso_iaxis(2)
        ELSE
          CALL deallocate_aniso_iaxis(4)
        ENDIF
        !
        IF (fbw) THEN
          ! remove memory allocated for wsi, wsn, adeltai, aznormi, ashifti
          imelt = (2 + 3 * nbndfs * nks) * nsiw(itemp)
        ELSE
          ! remove memory allocated for wsi, adeltai, aznormi
          imelt = (1 + 2 * nbndfs * nks) * nsiw(itemp)
        ENDIF
        CALL mem_size_eliashberg(2, -imelt)
        IF (fbw) THEN
          IF (gridsamp == 2) THEN
            ! If using sparse-ir
            ! remove memory allocated for arrays related to IR (complex)
            imelt = 6 * ir_obj%nfreq_f + 1 * ir_obj%nfreq_b
            imelt = imelt * siz_ir
            CALL mem_size_eliashberg(4, -imelt)
            IF (positive_matsu) THEN
              ! remove memory allocated for arrays related to IR (real)
              imelt = 7 * (ir_obj%size + ir_obj%ntau)
              imelt = imelt * siz_ir
              CALL mem_size_eliashberg(2, -imelt)
            ELSE
              ! remove memory allocated for arrays related to IR (complex)
              imelt = 7 * (ir_obj%size + ir_obj%ntau)
              imelt = imelt * siz_ir
              CALL mem_size_eliashberg(4, -imelt)
            ENDIF
            ! remove memory allocated for weight_q
            imelt = siz_ir
            CALL mem_size_eliashberg(2, -imelt)
            ! remove memory allocated for num_js1
            imelt = nbndfs * nks
            CALL mem_size_eliashberg(2, -imelt)
            IF (icoulomb > 0) THEN
              ! remove memory allocated for adeltaip_cl, w_stat
              imelt = (nkstot_cl * nsiw(itemp) + nbnd_cl) * nbnd_cl
              ! remove memory allocated for adeltai_cl
              imelt = imelt + num_states * nsiw(itemp)
              CALL mem_size_eliashberg(2, -imelt)
              ! remove memory allocated for arrays related to IR (complex)
              imelt = ir_obj%nfreq_f
              imelt = imelt * siz_ir_cl
              CALL mem_size_eliashberg(4, -imelt)
              IF (positive_matsu) THEN
                ! remove memory allocated for arrays related to IR (real)
                imelt = ir_obj%size + ir_obj%ntau
                imelt = imelt * siz_ir_cl
                CALL mem_size_eliashberg(2, -imelt)
              ELSE
                ! remove memory allocated for arrays related to IR (complex)
                imelt = ir_obj%size + ir_obj%ntau
                imelt = imelt * siz_ir_cl
                CALL mem_size_eliashberg(4, -imelt)
              ENDIF
              ! remove memory allocated for weight_cl
              imelt = siz_ir_cl
              CALL mem_size_eliashberg(2, -imelt)
              ! remove memory allocated for num_js2, num_js3
              imelt = nbndfs * nks + num_states
              CALL mem_size_eliashberg(2, -imelt)
            ENDIF
            IF (iverbosity == 4) THEN
              ! remove memory allocated for gl_abs, fl_abs, and knll_abs
              imelt = 3 * nbndfs * nkfs * ir_obj%size
              CALL mem_size_eliashberg(2, -imelt)
            ENDIF
          ELSEIF (gridsamp == 3) THEN
            ! If using FFT
            !
            IF (positive_matsu) THEN
              n = 2 * nsiw(itemp)
            ELSE
              n = nsiw(itemp)
            ENDIF
            !
            narray(1) = 6
            narray(2) = 4
            ! remove memory allocated for fft_in1, fft_out1, fft_in2, fft_out2
            imelt = 2 * (narray(1) + narray(2)) * n
            ! "2" indicates that arrays are allocated for both input and output
            CALL mem_size_eliashberg(2, -imelt)
          ENDIF
        ELSE
          ! not fbw
          IF (gridsamp == 3) THEN
            ! If using FFT
            !
            IF (positive_matsu) THEN
              n = 2 * nsiw(itemp)
            ELSE
              n = nsiw(itemp)
            ENDIF
            !
            narray(1) = 5
            narray(2) = 3
            ! remove memory allocated for fft_in1, fft_out1, fft_in2, fft_out2
            imelt = 2 * (narray(1) + narray(2)) * n
            ! "2" indicates that arrays are allocated for both input and output
            CALL mem_size_eliashberg(2, -imelt)
          ENDIF
        ENDIF
      ENDIF
      !
      tcpu = get_clock('aniso_iaxis')
      WRITE(stdout, '(5x, a, i3, a, f18.2, a)') 'itemp = ', itemp, '   total cpu time :', tcpu, ' secs'
      WRITE(stdout, '(a)') ' '
      !
      IF (lpade .OR. lacon) THEN
        CALL deallocate_aniso_raxis(2)
        IF (fbw) THEN
          ! remove memory allocated for ws, delta, znorm, shift, adelta, aznorm, ashift
          imelt = nsw + 2 * (3 + 3 * nbndfs * nks) * nsw
        ELSE
          ! remove memory allocated for ws, delta, znorm, adelta, aznorm
          imelt = nsw + 2 * (2 + 2 * nbndfs * nks) * nsw
        ENDIF
        CALL mem_size_eliashberg(2, -imelt)
      ENDIF
    ENDDO ! itemp
    !
    IF (fbw .AND. (gridsamp == 2)) THEN
      ! remove memory allocated for ir objects
      CALL finalize_ir(ir_obj)
    ENDIF
    !
    CALL deallocate_aniso()
    !
    CALL stop_clock('aniso_iaxis')
    !
    RETURN
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE eliashberg_aniso_iaxis
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    SUBROUTINE sum_eliashberg_aniso_iaxis_wrapper(itemp, iter, conv, ir_obj)
    !-----------------------------------------------------------------------
    !!
    !! This routine solves the anisotropic Eliashberg equations on the imaginary-axis
    !!
    !! SH: Modified to allow for fbw calculations; and
    !!       re-ordered "deltai, ..." arrays' indices for efficiency (Nov 2021).
    !! HM: Modified to allow for fbw calculations with sparse-ir sampling (May 2024).
    !!
    USE kinds,             ONLY : DP
    USE control_flags,     ONLY : iverbosity
    USE global_var,        ONLY : gtemp
    USE input,             ONLY : nsiter, conv_thr_iaxis, fsthick, fbw, muchem, &
                                  imag_read, a_gap0, ngaussw, degaussw, gridsamp, &
                                  positive_matsu, icoulomb
    USE supercond_common,  ONLY : nqfs, wkfs, w0g, ekfs, nkfs, nbndfs, dosef, ef0, &
                                  nsiw, wsn, wsi, wsphmax, gap0, akeri, limag_fly, &
                                  deltai, znormi, shifti, adeltai, adeltaip, aznormi, aznormip, &
                                  naznormi, ashifti, ashiftip, muintr, fft_in1, fft_out1, &
                                  fft_in2, fft_out2, ir_giw, ir_gl, ir_gtau, ir_knliw, &
                                  ir_knll, ir_knltau, ir_cvliw, ir_cvll, ir_cvltau, &
                                  nbnd_offset, nbnd_cl, nbndfs_all, ik_cl_to_fs, nkstot_cl, &
                                  ibnd_kfs_all_to_kfs, adeltai_cl, adeltaip_cl, w_stat, &
                                  ir_giw_cl, ir_gl_cl, ir_gtau_cl, siz_ir, siz_ir_cl, weight_q, &
                                  weight_cl, num_js1, num_js2, num_js3, gl_abs, fl_abs, knll_abs, &
                                  ir_gl_d, ir_gtau_d, ir_knll_d, ir_knltau_d, ir_cvll_d, &
                                  ir_cvltau_d, ir_gl_cl_d, ir_gtau_cl_d
    USE ep_constants,      ONLY : kelvin2eV, pi, zero, one, eps8, czero, cone, two, &
                                  ryd2ev, eps16
    USE io_global,         ONLY : stdout, ionode_id
    USE mp_global,         ONLY : inter_pool_comm
    USE mp_world,          ONLY : mpime
    USE mp,                ONLY : mp_bcast, mp_barrier, mp_sum, mp_max, mp_min
    USE parallelism,       ONLY : fkbounds, para_bounds
    USE low_lvl,           ONLY : mem_size_eliashberg, memlt_eliashberg
    USE fft_scalar,        ONLY : cft_1z
    USE utilities,         ONLY : dos_ef_seq
    USE sparse_ir,         ONLY : IR
    !
    IMPLICIT NONE
    !
    LOGICAL, INTENT(inout) :: conv
    !! True if the calculation is converged
    INTEGER, INTENT(in) :: itemp
    !! Counter on temperature index
    INTEGER, INTENT(in) :: iter
    !! Counter on iteration steps
    TYPE(IR), INTENT(IN), OPTIONAL :: ir_obj
    !! which contains ir objects such as basis functions
    !
    ! Local variables
    LOGICAL :: linsidei
    !! used to determine if the state (ibnd_cl, ik_cl) is outside the window.
    LOGICAL :: linsidei2
    !! used to determine if the state (ibnd_cl, ik_cl) is outside the window.
    LOGICAL :: linsidei3
    !! used to determine if the state (ibnd_cl, ik_cl) is outside the window.
    INTEGER :: iw
    !! Counter on frequency imag-axis
    INTEGER :: narray(2)
    !! integers to determine the size of arrays of FFT
    INTEGER :: ik
    !! Counter on k-points
    INTEGER :: lower_bnd, upper_bnd
    !! Lower/upper bound index after k paral
    INTEGER :: nks
    !! Number of k points per pool
    INTEGER :: ibnd
    !! Counter on bands at k
    INTEGER :: ik_cl
    !! Counter on k-points
    INTEGER :: ibnd_cl
    !! Counter on bands at k
    INTEGER :: ibndfs
    !! Counter on bands in the fsthick window
    INTEGER :: n
    !! (wsn(nsiw(itemp)) + 1) + 1
    INTEGER, SAVE :: ns
    !! mu_inter parameters
    INTEGER :: is_start, is_stop
    !! Lower/upper bound index after paral for states
    INTEGER :: num_states
    !! Number of states per pool: num_states = is_stop - is_start + 1
    INTEGER(8) :: imelt
    !! Counter memory
    INTEGER :: ierr
    !! Error status
    !
    REAL(KIND = DP) :: absdelta, reldelta, errdelta
    !! Errors in supercond. gap
    REAL(KIND = DP) :: weight
    !! Factor in supercond. equations
    REAL(KIND = DP), SAVE :: nel, nstate
    !! mu_inter parameters
    REAL(KIND = DP) :: inv_dos
    !! Invese dos inv_dos = 1/dosef. Defined for efficiency reason
    REAL(KIND = DP) :: dos_muintr
    !! Density of states at the updated chemical potential
    REAL(KIND = DP) :: omega
    !! sqrt{w^2+\delta^2}
    REAL(KIND = DP) :: dFE
    !! free energy difference between supercond and normal states
    REAL(KIND = DP), EXTERNAL :: wgauss
    !! Compute the approximate theta function. Used to calculate nel
    REAL(KIND = DP), ALLOCATABLE, SAVE :: inv_wsi(:)
    !! Invese imaginary freq. inv_wsi = 1/wsi. Defined for efficiency reason
    REAL(KIND = DP), ALLOCATABLE, SAVE :: deltaold(:)
    !! supercond. gap from previous iteration
    !
    IF (fbw .AND. (gridsamp == 2)) THEN
      IF (.NOT. PRESENT(ir_obj)) THEN
        CALL errore('sum_eliashberg_aniso_iaxis_wrapper', 'Error: ir_obj is not given while gridsamp = 2', 1)
      ENDIF
    ENDIF
    !
    inv_dos = one / dosef
    !
    CALL fkbounds(nkfs, lower_bnd, upper_bnd)
    !
    nks = upper_bnd - lower_bnd + 1
    !
    IF (fbw .AND. (gridsamp == 2) .AND. (icoulomb > 0)) THEN
      !
      ! HM: Because nkstot_cl can be quite small,
      !     we do a combined parallelization on band index and k point index.
      CALL para_bounds(is_start, is_stop, nbnd_cl * nkstot_cl)
      num_states = is_stop - is_start + 1
      !
    ENDIF
    !
    IF (iter == 1) THEN
      IF (itemp == 1 .OR. (itemp == 2 .AND. imag_read)) THEN
        ! SH: calculate the input parameters for mu_inter
        nel = zero
        nstate = zero
        ns = 0
        DO ik = lower_bnd, upper_bnd
          DO ibnd = 1, nbndfs
            IF (ABS(ekfs(ibnd, ik) - ef0) < fsthick) THEN
              ns = ns + 1
              ! HP: The initial guess is based on the FD dist. at 0 K.
              nstate = nstate + 0.5d0 * wkfs(ik)
              nel = nel + wkfs(ik) * wgauss( (ef0 - ekfs(ibnd, ik)) / degaussw, ngaussw )
            ENDIF
          ENDDO
        ENDDO
        CALL mp_sum(ns, inter_pool_comm)
        CALL mp_sum(nstate, inter_pool_comm)
        CALL mp_sum(nel, inter_pool_comm)
        CALL mp_barrier(inter_pool_comm)
        !
      ENDIF
      !
      IF (fbw .AND. (gridsamp == 2)) THEN
        !
        ! get the size of required memory for num_js1
        imelt = nbndfs * nks
        !
        CALL mem_size_eliashberg(2, imelt)
        !
        ALLOCATE(num_js1(nbndfs, lower_bnd:upper_bnd), STAT = ierr)
        IF (ierr /= 0) CALL errore('sum_eliashberg_aniso_iaxis_wrapper', 'Error allocating num_js1', 1)
        !
        num_js1(:, :) = 0
        !
        IF (icoulomb > 0) THEN
          !
          ! get the size of required memory for num_js2, num_js3
          imelt = nbndfs * nks + num_states
          !
          CALL mem_size_eliashberg(2, imelt)
          !
          ALLOCATE(num_js2(nbndfs, lower_bnd:upper_bnd), STAT = ierr)
          IF (ierr /= 0) CALL errore('sum_eliashberg_aniso_iaxis_wrapper', 'Error allocating num_js2', 1)
          ALLOCATE(num_js3(is_start:is_stop), STAT = ierr)
          IF (ierr /= 0) CALL errore('sum_eliashberg_aniso_iaxis_wrapper', 'Error allocating num_js3', 1)
          !
          num_js2(:, :) = 0
          num_js3(:) = 0
          !
        ENDIF
        !
      ENDIF ! fbw
      !
      IF (fbw .AND. (gridsamp == 2)) THEN
        CALL count_states()
      ENDIF
      !
      ! RM - adeltai, aznormi, naznormi, are defined per k-points per pool
      !
      ! get the size of required memory for inv_wsi, deltai, znormi,
      ! adeltai, adeltaip, aznormi, naznormi
      imelt = (3 + 3 * nbndfs * nks + nbndfs * nkfs) * nsiw(itemp)
      CALL mem_size_eliashberg(2, imelt)
      !
      ALLOCATE(inv_wsi(nsiw(itemp)), STAT = ierr)
      IF (ierr /= 0) CALL errore('sum_eliashberg_aniso_iaxis_wrapper', 'Error allocating inv_wsi', 1)
      !
      ALLOCATE(deltai(nsiw(itemp)), STAT = ierr)
      IF (ierr /= 0) CALL errore('sum_eliashberg_aniso_iaxis_wrapper', 'Error allocating deltai', 1)
      ALLOCATE(znormi(nsiw(itemp)), STAT = ierr)
      IF (ierr /= 0) CALL errore('sum_eliashberg_aniso_iaxis_wrapper', 'Error allocating znormi', 1)
      !
      ALLOCATE(adeltai(nsiw(itemp), nbndfs, lower_bnd:upper_bnd), STAT = ierr)
      IF (ierr /= 0) CALL errore('sum_eliashberg_aniso_iaxis_wrapper', 'Error allocating adeltai', 1)
      ALLOCATE(aznormi(nsiw(itemp), nbndfs, lower_bnd:upper_bnd), STAT = ierr)
      IF (ierr /= 0) CALL errore('sum_eliashberg_aniso_iaxis_wrapper', 'Error allocating aznormi', 1)
      ALLOCATE(naznormi(nsiw(itemp), nbndfs, lower_bnd:upper_bnd), STAT = ierr)
      IF (ierr /= 0) CALL errore('sum_eliashberg_aniso_iaxis_wrapper', 'Error allocating naznormi', 1)
      !
      ALLOCATE(adeltaip(nsiw(itemp), nbndfs, nkfs), STAT = ierr)
      IF (ierr /= 0) CALL errore('sum_eliashberg_aniso_iaxis_wrapper', 'Error allocating adeltaip', 1)
      !
      adeltaip(:, :, :) = zero
      IF (fbw) THEN
        ! get the size of required memory for shifti, aznormip, ashifti, ashiftip
        imelt = (1 + 1 * nbndfs * nks + 2 * nbndfs * nkfs) * nsiw(itemp)
        CALL mem_size_eliashberg(2, imelt)
        !
        ! SH: to allocate and initiate the fbw run variables
        ALLOCATE(shifti(nsiw(itemp)), STAT = ierr)
        IF (ierr /= 0) CALL errore('sum_eliashberg_aniso_iaxis_fbw', 'Error allocating shifti', 1)
        ALLOCATE(ashifti(nsiw(itemp), nbndfs, lower_bnd:upper_bnd), STAT = ierr)
        IF (ierr /= 0) CALL errore('sum_eliashberg_aniso_iaxis_fbw', 'Error allocating ashifti', 1)
        ALLOCATE(aznormip(nsiw(itemp), nbndfs, nkfs), STAT = ierr)
        IF (ierr /= 0) CALL errore('sum_eliashberg_aniso_iaxis_fbw', 'Error allocating aznormip', 1)
        ALLOCATE(ashiftip(nsiw(itemp), nbndfs, nkfs), STAT = ierr)
        IF (ierr /= 0) CALL errore('sum_eliashberg_aniso_iaxis_fbw', 'Error allocating ashiftip', 1)
        !
        aznormip(:, :, :) = one
        ashiftip(:, :, :) = zero
        !
        IF (gridsamp == 2) THEN
          !
          ! get the size of required memory for arrays related to IR (complex)
          imelt = 6 * ir_obj%nfreq_f + 1 * ir_obj%nfreq_b
          imelt = imelt * siz_ir
          !
          CALL mem_size_eliashberg(4, imelt)
          !
          ALLOCATE(ir_giw(3 * siz_ir, ir_obj%nfreq_f), STAT = ierr)
          IF (ierr /= 0) CALL errore('sum_eliashberg_aniso_iaxis_wrapper', 'Error allocating ir_giw', 1)
          ALLOCATE(ir_knliw(siz_ir, ir_obj%nfreq_b), STAT = ierr)
          IF (ierr /= 0) CALL errore('sum_eliashberg_aniso_iaxis_wrapper', 'Error allocating ir_knliw', 1)
          ALLOCATE(ir_cvliw(3 * siz_ir, ir_obj%nfreq_f), STAT = ierr)
          IF (ierr /= 0) CALL errore('sum_eliashberg_aniso_iaxis_wrapper', 'Error allocating ir_cvliw', 1)
          ir_giw(:, :) = czero
          ir_knliw(:, :) = czero
          ir_cvliw(:, :) = czero
          !
          ! get the size of required memory for weight_q
          imelt = siz_ir
          CALL mem_size_eliashberg(2, imelt)
          ALLOCATE(weight_q(siz_ir), STAT = ierr)
          IF (ierr /= 0) CALL errore('sum_eliashberg_aniso_iaxis_wrapper', 'Error allocating weight_q', 1)
          !
          IF (positive_matsu) THEN
            ! get the size of required memory for arrays related to IR (real)
            imelt = 7 * (ir_obj%size + ir_obj%ntau)
            imelt = imelt * siz_ir
            !
            CALL mem_size_eliashberg(2, imelt)
            ALLOCATE(ir_gl_d(3 * siz_ir, ir_obj%size), STAT = ierr)
            IF (ierr /= 0) CALL errore('sum_eliashberg_aniso_iaxis_wrapper', 'Error allocating ir_gl_d', 1)
            ALLOCATE(ir_gtau_d(3 * siz_ir, ir_obj%ntau), STAT = ierr)
            IF (ierr /= 0) CALL errore('sum_eliashberg_aniso_iaxis_wrapper', 'Error allocating ir_gtau_d', 1)
            ALLOCATE(ir_knll_d(siz_ir, ir_obj%size), STAT = ierr)
            IF (ierr /= 0) CALL errore('sum_eliashberg_aniso_iaxis_wrapper', 'Error allocating ir_knll_d', 1)
            ALLOCATE(ir_knltau_d(siz_ir, ir_obj%ntau), STAT = ierr)
            IF (ierr /= 0) CALL errore('sum_eliashberg_aniso_iaxis_wrapper', 'Error allocating ir_knltau_d', 1)
            ALLOCATE(ir_cvll_d(3 * siz_ir, ir_obj%size), STAT = ierr)
            IF (ierr /= 0) CALL errore('sum_eliashberg_aniso_iaxis_wrapper', 'Error allocating ir_cvll_d', 1)
            ALLOCATE(ir_cvltau_d(3 * siz_ir, ir_obj%ntau), STAT = ierr)
            IF (ierr /= 0) CALL errore('sum_eliashberg_aniso_iaxis_wrapper', 'Error allocating ir_cvltau_d', 1)
            ir_gl_d(:, :) = zero
            ir_gtau_d(:, :) = zero
            ir_knll_d(:, :) = zero
            ir_knltau_d(:, :) = zero
            ir_cvll_d(:, :) = zero
            ir_cvltau_d(:, :) = zero
          ELSE
            ! get the size of required memory for arrays related to IR (complex)
            imelt = 7 * (ir_obj%size + ir_obj%ntau)
            imelt = imelt * siz_ir
            !
            CALL mem_size_eliashberg(4, imelt)
            ALLOCATE(ir_gl(3 * siz_ir, ir_obj%size), STAT = ierr)
            IF (ierr /= 0) CALL errore('sum_eliashberg_aniso_iaxis_wrapper', 'Error allocating ir_gl', 1)
            ALLOCATE(ir_gtau(3 * siz_ir, ir_obj%ntau), STAT = ierr)
            IF (ierr /= 0) CALL errore('sum_eliashberg_aniso_iaxis_wrapper', 'Error allocating ir_gtau', 1)
            ALLOCATE(ir_knll(siz_ir, ir_obj%size), STAT = ierr)
            IF (ierr /= 0) CALL errore('sum_eliashberg_aniso_iaxis_wrapper', 'Error allocating ir_knll', 1)
            ALLOCATE(ir_knltau(siz_ir, ir_obj%ntau), STAT = ierr)
            IF (ierr /= 0) CALL errore('sum_eliashberg_aniso_iaxis_wrapper', 'Error allocating ir_knltau', 1)
            ALLOCATE(ir_cvll(3 * siz_ir, ir_obj%size), STAT = ierr)
            IF (ierr /= 0) CALL errore('sum_eliashberg_aniso_iaxis_wrapper', 'Error allocating ir_cvll', 1)
            ALLOCATE(ir_cvltau(3 * siz_ir, ir_obj%ntau), STAT = ierr)
            IF (ierr /= 0) CALL errore('sum_eliashberg_aniso_iaxis_wrapper', 'Error allocating ir_cvltau', 1)
            ir_gl(:, :) = czero
            ir_gtau(:, :) = czero
            ir_knll(:, :) = czero
            ir_knltau(:, :) = czero
            ir_cvll(:, :) = czero
            ir_cvltau(:, :) = czero
          ENDIF
          !
          IF (icoulomb > 0) THEN
            !
            ! get the size of required memory for adeltaip_cl, w_stat
            imelt = (nkstot_cl * nsiw(itemp) + nbnd_cl) * nbnd_cl
            ! get the size of required memory for adeltai_cl
            imelt = imelt + num_states * nsiw(itemp)
            CALL mem_size_eliashberg(2, imelt)
            !
            ALLOCATE(adeltai_cl(nsiw(itemp), is_start:is_stop), STAT = ierr)
            IF (ierr /= 0) CALL errore('sum_eliashberg_aniso_iaxis_wrapper', 'Error allocating adeltai_cl', 1)
            ALLOCATE(adeltaip_cl(nsiw(itemp), nbnd_cl, nkstot_cl), STAT = ierr)
            IF (ierr /= 0) CALL errore('sum_eliashberg_aniso_iaxis_wrapper', 'Error allocating adeltaip_cl', 1)
            ALLOCATE(w_stat(nbnd_cl, nbnd_cl), STAT = ierr)
            IF (ierr /= 0) CALL errore('sum_eliashberg_aniso_iaxis_wrapper', 'Error allocating w_stat', 1)
            adeltaip(:, :, :) = zero
            adeltaip_cl(:, :, :) = zero
            !
            ! get the size of required memory for arrays related to IR (complex)
            imelt = ir_obj%nfreq_f
            imelt = imelt * siz_ir_cl
            CALL mem_size_eliashberg(4, imelt)
            !
            ALLOCATE(ir_giw_cl(siz_ir_cl, ir_obj%nfreq_f), STAT = ierr)
            IF (ierr /= 0) CALL errore('sum_eliashberg_aniso_iaxis_wrapper', 'Error allocating ir_giw_cl', 1)
            ir_giw_cl(:, :) = czero
            !
            IF (positive_matsu) THEN
              ! get the size of required memory for arrays related to IR (real)
              imelt = ir_obj%size + ir_obj%ntau
              imelt = imelt * siz_ir_cl
              CALL mem_size_eliashberg(2, imelt)
              !
              ALLOCATE(ir_gl_cl_d(siz_ir_cl, ir_obj%size), STAT = ierr)
              IF (ierr /= 0) CALL errore('sum_eliashberg_aniso_iaxis_wrapper', 'Error allocating ir_gl_cl_d', 1)
              ALLOCATE(ir_gtau_cl_d(siz_ir_cl, ir_obj%ntau), STAT = ierr)
              IF (ierr /= 0) CALL errore('sum_eliashberg_aniso_iaxis_wrapper', 'Error allocating ir_gtau_cl_d', 1)
              ir_gl_cl_d(:, :) = zero
              ir_gtau_cl_d(:, :) = zero
              !
            ELSE
              ! get the size of required memory for arrays related to IR (complex)
              imelt = ir_obj%size + ir_obj%ntau
              imelt = imelt * siz_ir_cl
              CALL mem_size_eliashberg(4, imelt)
              !
              ALLOCATE(ir_gl_cl(siz_ir_cl, ir_obj%size), STAT = ierr)
              IF (ierr /= 0) CALL errore('sum_eliashberg_aniso_iaxis_wrapper', 'Error allocating ir_gl_cl', 1)
              ALLOCATE(ir_gtau_cl(siz_ir_cl, ir_obj%ntau), STAT = ierr)
              IF (ierr /= 0) CALL errore('sum_eliashberg_aniso_iaxis_wrapper', 'Error allocating ir_gtau_cl', 1)
              ir_gl_cl(:, :) = czero
              ir_gtau_cl(:, :) = czero
              !
            ENDIF
            ! get the size of required memory for weight_cl
            imelt = siz_ir_cl
            CALL mem_size_eliashberg(2, imelt)
            ALLOCATE(weight_cl(siz_ir_cl), STAT = ierr)
            IF (ierr /= 0) CALL errore('sum_eliashberg_aniso_iaxis_wrapper', 'Error allocating weight_cl', 1)
            !
          ENDIF
          !
          IF (iverbosity == 4) THEN
            !
            ! get the size of required memory for gl_abs, fl_abs, and knll_abs
            imelt = 3 * nbndfs * nkfs * ir_obj%size
            CALL mem_size_eliashberg(2, imelt)
            !
            ALLOCATE(gl_abs(ir_obj%size, nbndfs, nkfs), STAT = ierr)
            IF (ierr /= 0) CALL errore('sum_eliashberg_aniso_iaxis_wrapper', 'Error allocating gl_abs', 1)
            ALLOCATE(fl_abs(ir_obj%size, nbndfs, nkfs), STAT = ierr)
            IF (ierr /= 0) CALL errore('sum_eliashberg_aniso_iaxis_wrapper', 'Error allocating fl_abs', 1)
            ALLOCATE(knll_abs(ir_obj%size, nbndfs, nkfs), STAT = ierr)
            IF (ierr /= 0) CALL errore('sum_eliashberg_aniso_iaxis_wrapper', 'Error allocating knll_abs', 1)
            gl_abs(:, :, :) = zero
            fl_abs(:, :, :) = zero
            knll_abs(:, :, :) = zero
            !
          ENDIF
          !
        ELSEIF (gridsamp == 3) THEN
          !
          IF (positive_matsu) THEN
            n = 2 * nsiw(itemp)
          ELSE
            n = nsiw(itemp)
          ENDIF
          !
          narray(1) = 6
          narray(2) = 4
          ! get the size of required memory for arrays in FFT
          imelt = 2 * (narray(1) + narray(2)) * n
          ! "2" indicates that arrays are allocated for both input and output
          CALL mem_size_eliashberg(4, imelt)
          !
          ALLOCATE(fft_in1(n * narray(1)), STAT = ierr)
          IF (ierr /= 0) CALL errore('sum_eliashberg_aniso_iaxis_wrapper', 'Error allocating fft_in1', 1)
          ALLOCATE(fft_out1(n * narray(1)), STAT = ierr)
          IF (ierr /= 0) CALL errore('sum_eliashberg_aniso_iaxis_wrapper', 'Error allocating fft_out1', 1)
          ALLOCATE(fft_in2(n * narray(2)), STAT = ierr)
          IF (ierr /= 0) CALL errore('sum_eliashberg_aniso_iaxis_wrapper', 'Error allocating fft_in2', 1)
          ALLOCATE(fft_out2(n * narray(2)), STAT = ierr)
          IF (ierr /= 0) CALL errore('sum_eliashberg_aniso_iaxis_wrapper', 'Error allocating fft_out2', 1)
          fft_in1(:) = czero
          fft_out1(:) = czero
          fft_in2(:) = czero
          fft_out2(:) = czero
          !
        ENDIF ! positive_matsu, gridsamp
        !
      ELSE
        ! not fbw
        IF (gridsamp == 3) THEN
          !
          n = 2 * nsiw(itemp)
          !
          narray(1) = 5
          narray(2) = 3
          ! get the size of required memory for arrays in FFT
          imelt = 2 * (narray(1) + narray(2)) * n
          ! "2" indicates that arrays are allocated for both input and output
          CALL mem_size_eliashberg(4, imelt)
          !
          ALLOCATE(fft_in1(n * narray(1)), STAT = ierr)
          IF (ierr /= 0) CALL errore('sum_eliashberg_aniso_iaxis_wrapper', 'Error allocating fft_in1', 1)
          ALLOCATE(fft_out1(n * narray(1)), STAT = ierr)
          IF (ierr /= 0) CALL errore('sum_eliashberg_aniso_iaxis_wrapper', 'Error allocating fft_out1', 1)
          ALLOCATE(fft_in2(n * narray(2)), STAT = ierr)
          IF (ierr /= 0) CALL errore('sum_eliashberg_aniso_iaxis_wrapper', 'Error allocating fft_in2', 1)
          ALLOCATE(fft_out2(n * narray(2)), STAT = ierr)
          IF (ierr /= 0) CALL errore('sum_eliashberg_aniso_iaxis_wrapper', 'Error allocating fft_out2', 1)
          fft_in1(:) = czero
          fft_out1(:) = czero
          fft_in2(:) = czero
          fft_out2(:) = czero
          !
        ENDIF ! gridsamp
      ENDIF
      !
      ! SH: set the initial value of superconducting chemical potential;
      !     will be used only for the fbw calculations!
      muintr = ef0
      !
      DO iw = 1, nsiw(itemp)
        inv_wsi(iw) = one / wsi(iw)
      ENDDO
      !
      DO ik = 1, nkfs
        DO ibnd = 1, nbndfs
          IF (ABS(ekfs(ibnd, ik) - ef0) < fsthick) THEN
            DO iw = 1, nsiw(itemp)
              ! a_gap0 is set to 1.0d0 as default
              IF (a_gap0 <= -eps8) THEN
                IF (ABS(wsi(iw)) < 2.d0 * wsphmax) THEN
                  adeltaip(iw, ibnd, ik) = gap0
                ELSE
                  adeltaip(iw, ibnd, ik) = zero
                ENDIF
              ELSEIF (ABS(a_gap0) < eps8) THEN
                adeltaip(iw, ibnd, ik) = gap0
              ELSE
                adeltaip(iw, ibnd, ik) = gap0 / (one + a_gap0 * (wsi(iw) / wsphmax)**2)
              ENDIF
            ENDDO
          ENDIF
        ENDDO ! ibnd
      ENDDO ! ik
      !
      IF (fbw .AND. (gridsamp == 2) .AND. (icoulomb > 0)) THEN
        linsidei = .TRUE.
        linsidei2 = .TRUE.
        linsidei3 = .TRUE.
        DO ik_cl = 1, nkstot_cl
          !
          ! initilize linsidei
          linsidei = .TRUE.
          IF (ik_cl_to_fs(ik_cl) == 0) linsidei = .FALSE.
          DO ibnd_cl = 1, nbnd_cl
            ! initilize linsidei2
            linsidei2 = .TRUE.
            IF (((ibnd_cl - nbnd_offset) < 1) .OR. ((ibnd_cl - nbnd_offset) > nbndfs_all)) linsidei2 = .FALSE.
            !
            ! initilize linsidei3
            linsidei3 = .TRUE.
            IF (linsidei .AND. linsidei2) THEN
              ibndfs = ibnd_kfs_all_to_kfs(ibnd_cl - nbnd_offset, ik_cl_to_fs(ik_cl))
              IF ((ibndfs < 1) .OR. (ibndfs > nbndfs)) THEN
                linsidei3 = .FALSE.
              ELSE
                ! Always use ekfs instead of ek_cl to determine
                ! if the state (ibnd_cl, ik_cl) is outside the window.
                IF (ABS(ekfs(ibndfs, ik_cl_to_fs(ik_cl)) - ef0) >= fsthick) linsidei3 = .FALSE.
              ENDIF
            ENDIF
            !
            IF (linsidei .AND. linsidei2 .AND. linsidei3) THEN
              ! When the state (ibnd_cl, ik_cl) is inside of the fsthick window,
              ! no need to calculate adeltai_cl(:, ibnd_cl, ik_cl)
              adeltaip_cl(:, ibnd_cl, ik_cl) = zero
            ELSE
              DO iw = 1, nsiw(itemp)
                ! a_gap0 is set to 1.0d0 as default
                IF (a_gap0 <= -eps8) THEN
                  IF (ABS(wsi(iw)) < 2.d0 * wsphmax) THEN
                    adeltaip_cl(iw, ibnd_cl, ik_cl) = gap0
                  ELSE
                    adeltaip_cl(iw, ibnd_cl, ik_cl) = zero
                  ENDIF
                ELSEIF (ABS(a_gap0) < eps8) THEN
                  adeltaip_cl(iw, ibnd_cl, ik_cl) = gap0
                ELSE
                  adeltaip_cl(iw, ibnd_cl, ik_cl) = gap0 / (one + a_gap0 * (wsi(iw) / wsphmax)**2)
                ENDIF
              ENDDO ! iw
            ENDIF ! linsidei .AND. linsidei2 .AND. linsidei3
          ENDDO ! ibnd_cl
        ENDDO ! ik_cl
      ENDIF
      !
      CALL memlt_eliashberg(itemp, 'imag')
      IF (.NOT. limag_fly) CALL kernel_aniso_iaxis(itemp)
      !
    ENDIF ! iter
    !
    ! SH: for the case of fbw runs
    IF (fbw .AND. ((gridsamp <= 0) .OR. (positive_matsu .AND. (gridsamp == 1)))) THEN
      CALL sum_eliashberg_aniso_iaxis_fbw_simple(itemp, nel, nstate)
    ELSEIF (fbw .AND. (gridsamp == 2) .AND. (icoulomb == 0)) THEN
      CALL sum_eliashberg_aniso_iaxis_fbw_ir(itemp, iter, ns, nel, nstate, ir_obj)
    ELSEIF (fbw .AND. (gridsamp == 2) .AND. (icoulomb > 0)) THEN
      CALL sum_eliash_aniso_iaxis_fbw_ir_coul(itemp, iter, &
                                                     ns, nel, nstate, ir_obj)
    ELSEIF (fbw .AND. (gridsamp == 3)) THEN
      CALL sum_eliashberg_aniso_iaxis_fbw_fft(itemp, nel, nstate, narray)
    ELSEIF ((.NOT. fbw) .AND. (positive_matsu .AND. (gridsamp <= 1))) THEN
      CALL sum_eliashberg_aniso_iaxis_fsr_simple(itemp)
    ELSEIF ((.NOT. fbw) .AND. (positive_matsu .AND. (gridsamp == 3))) THEN
      CALL sum_eliashberg_aniso_iaxis_fsr_fft(itemp, narray)
    ELSE
      CALL errore('sum_eliashberg_aniso_iaxis_wrapper', 'Unexpected combination of input values.', 1)
    ENDIF ! fbw
    !
    IF (fbw) THEN
      deltai(:) = zero
      znormi(:) = zero
      shifti(:) = zero
      DO iw = 1, nsiw(itemp) ! loop over omega
        DO ik = lower_bnd, upper_bnd
          DO ibnd = 1, nbndfs
            IF (ABS(ekfs(ibnd, ik) - ef0) < fsthick) THEN
              weight     = 0.5d0 * wkfs(ik) * w0g(ibnd, ik) * inv_dos
              znormi(iw) = znormi(iw) + weight * aznormi(iw, ibnd, ik)
              deltai(iw) = deltai(iw) + weight * adeltai(iw, ibnd, ik)
              shifti(iw) = shifti(iw) + weight * ashifti(iw, ibnd, ik)
              naznormi(iw, ibnd, ik) = 1.d0 + gtemp(itemp) * naznormi(iw, ibnd, ik) * inv_wsi(iw)
              ! Eqs.(21)-(22) in Margine and Giustino, PRB 87, 024505 (2013)
              aznormi(iw, ibnd, ik) = 1.d0 + gtemp(itemp) * aznormi(iw, ibnd, ik) * inv_wsi(iw)
              adeltai(iw, ibnd, ik) = gtemp(itemp) * adeltai(iw, ibnd, ik) / aznormi(iw, ibnd, ik)
              ashifti(iw, ibnd, ik) = - gtemp(itemp) * ashifti(iw, ibnd, ik)
            ENDIF
          ENDDO ! ibnd
        ENDDO ! ik
      ENDDO ! iw
      !
      ! collect contributions from all pools
      CALL mp_sum(deltai, inter_pool_comm)
      CALL mp_sum(znormi, inter_pool_comm)
      CALL mp_sum(shifti, inter_pool_comm)
      CALL mp_barrier(inter_pool_comm)
    ELSE ! not fbw
      deltai(:) = zero
      znormi(:) = zero
      DO iw = 1, nsiw(itemp) ! loop over omega
        DO ik = lower_bnd, upper_bnd
          DO ibnd = 1, nbndfs
            IF (ABS(ekfs(ibnd, ik) - ef0) < fsthick) THEN
              weight     = 0.5d0 * wkfs(ik) * w0g(ibnd, ik) * inv_dos
              znormi(iw) = znormi(iw) + weight * aznormi(iw, ibnd, ik)
              deltai(iw) = deltai(iw) + weight * adeltai(iw, ibnd, ik)
              naznormi(iw, ibnd, ik) = 1.d0 + pi * gtemp(itemp) * naznormi(iw, ibnd, ik) * inv_wsi(iw)
              ! Eqs.(21)-(22) in Margine and Giustino, PRB 87, 024505 (2013)
              aznormi(iw, ibnd, ik) = 1.d0 + pi * gtemp(itemp) * aznormi(iw, ibnd, ik) * inv_wsi(iw)
              adeltai(iw, ibnd, ik) = pi * gtemp(itemp) * adeltai(iw, ibnd, ik) / aznormi(iw, ibnd, ik)
            ENDIF
          ENDDO ! ibnd
        ENDDO ! ik
      ENDDO ! iw
      !
      ! collect contributions from all pools
      CALL mp_sum(deltai, inter_pool_comm)
      CALL mp_sum(znormi, inter_pool_comm)
      CALL mp_barrier(inter_pool_comm)
    ENDIF
    !
    IF (mpime == ionode_id) THEN
      !
      IF (iter == 1) THEN
        ALLOCATE(deltaold(nsiw(itemp)), STAT = ierr)
        IF (ierr /= 0) CALL errore('sum_eliashberg_aniso_iaxis_wrapper', 'Error allocating deltaold', 1)
        deltaold(:) = gap0
      ENDIF
      !
      absdelta = zero
      reldelta = zero
      DO iw = 1, nsiw(itemp) ! loop over omega
        IF (fbw) THEN
          znormi(iw) = 1.d0 + gtemp(itemp) * znormi(iw) * inv_wsi(iw)
          deltai(iw) = gtemp(itemp) * deltai(iw) / znormi(iw)
          shifti(iw) = - gtemp(itemp) * shifti(iw)
        ELSE
          znormi(iw) = 1.d0 + pi * gtemp(itemp) * znormi(iw) * inv_wsi(iw)
          deltai(iw) = pi * gtemp(itemp) * deltai(iw) / znormi(iw)
        ENDIF
        !
        reldelta   = reldelta + ABS(deltai(iw) - deltaold(iw))
        absdelta   = absdelta + ABS(deltai(iw))
      ENDDO ! iw
      !
      errdelta = reldelta / absdelta
      deltaold(:) = deltai(:)
      !
      IF (iter == 1 .AND. fbw)       WRITE(stdout, '(5x, a)') &
        '   iter      ethr          znormi      deltai [meV]   shifti [meV]       mu [eV]'
      IF (iter == 1 .AND. .NOT. fbw) WRITE(stdout, '(5x, a)') &
        '   iter      ethr          znormi      deltai [meV]'
      IF (fbw) THEN
        IF (.NOT. positive_matsu) THEN
          ! if positive_matsu =.FALSE. the index of the lowest frequency is not 1
          n = nsiw(itemp)/2 + 1
        ELSE
          n = 1
        ENDIF
        WRITE(stdout, '(5x, i6, 5ES15.6)') &
          iter, errdelta, znormi(n), deltai(n) * 1000.d0, shifti(n) * 1000.d0, muintr
      ELSEIF (.NOT. fbw) THEN
        WRITE(stdout, '(5x, i6, 3ES15.6)') &
          iter, errdelta, znormi(1), deltai(1) * 1000.d0
      ENDIF
      !      WRITE(stdout, '(5x, a, i6, a, ES15.6, a, ES15.6, a, ES15.6)') 'iter = ', iter, &
      !                    '   ethr = ', errdelta, '   znormi(1) = ', znormi(1), &
      !                    '   deltai(1) = ', deltai(1)
      !
      IF (errdelta < conv_thr_iaxis) conv = .TRUE.
      IF (conv .OR. iter == nsiter) THEN
        IF (fbw) THEN
          IF (.NOT. positive_matsu) THEN
            ! if positive_matsu =.FALSE. the index of the lowest frequency is not 1
            n = nsiw(itemp) / 2 + 1
            gap0 = deltai(n)
          ELSE
            gap0 = deltai(1)
          ENDIF
        ELSEIF (.NOT. fbw) THEN
          gap0 = deltai(1)
        ENDIF
      ENDIF
      !
      IF (conv .OR. iter == nsiter) THEN
        DEALLOCATE(deltaold, STAT = ierr)
        IF (ierr /= 0) CALL errore('sum_eliashberg_aniso_iaxis_wrapper', 'Error deallocating deltaold', 1)
      ENDIF
      !
      IF (conv) THEN
        WRITE(stdout, '(5x, a, i6)') 'Convergence was reached in nsiter = ', iter
        WRITE(stdout, '(a)') ' '
      ELSEIF (.NOT. conv .AND. iter == nsiter) THEN
        WRITE(stdout, '(5x, a, i6)') 'Convergence was not reached in nsiter = ', iter
        WRITE(stdout, '(5x, a)') 'Increase nsiter or reduce conv_thr_iaxis'
        WRITE(stdout, '(a)') ' '
      ENDIF
    ENDIF
    CALL mp_bcast(gap0, ionode_id, inter_pool_comm)
    CALL mp_bcast(conv, ionode_id, inter_pool_comm)
    CALL mp_barrier(inter_pool_comm)
    !
    IF (conv .OR. iter == nsiter) THEN
      !
      IF (fbw) THEN
        ! SH: write the chemical potential for fbw runs
        WRITE(stdout, '(5x, a, i3, a, ES20.10, a)') &
          'Chemical potential (itemp = ', itemp, ') = ', muintr, ' eV'
        IF (muchem) THEN
          ! degaussw is already converted to be in the unit of eV in read_eigenvalues
          dos_muintr = dos_ef_seq(ngaussw, degaussw, muintr, ekfs, wkfs, nkfs, nbndfs)
          dos_muintr = dos_muintr / two
          ! Because we already have DOS in the unit of 1/eV, we do not have to divide it by ryd2ev.
          WRITE(stdout, '(5x, a, i3, a)') &
            'DOS at the chemical potential (itemp = ', itemp, ')'
          WRITE(stdout, '(11x, a, ES20.10)') &
            '(states/spin/eV/Unit Cell) = ', dos_muintr
          ! degaussw is already converted to be in the unit of eV in read_eigenvalues
          dos_muintr = dos_ef_seq(ngaussw, degaussw, ef0, ekfs, wkfs, nkfs, nbndfs)
          dos_muintr = dos_muintr / two
          ! Because we already have DOS in the unit of 1/eV, we do not have to divide it by ryd2ev.
          WRITE(stdout, '(5x, a)') &
            'DOS at the non-interacting Fermi energy'
          WRITE(stdout, '(11x, a, ES20.10)') &
            '(states/spin/eV/Unit Cell) = ', dos_muintr
        ENDIF
        WRITE(stdout, '(a)') ' '
      ENDIF
      !
      ! Compute the free energy difference between the superconducting and normal states
      dFE = zero
      DO ik = lower_bnd, upper_bnd
        DO ibnd = 1, nbndfs
          IF (ABS(ekfs(ibnd, ik) - ef0) < fsthick) THEN
            DO iw = 1, nsiw(itemp)
              weight = 0.5d0 * wkfs(ik) * w0g(ibnd,ik)
              omega = DSQRT(wsi(iw) * wsi(iw) + adeltai(iw, ibnd, ik) * adeltai(iw, ibnd, ik))
              dFE = dFE - weight * (omega - ABS(wsi(iw))) &
                * (aznormi(iw, ibnd, ik) - naznormi(iw, ibnd, ik) * ABS(wsi(iw)) / omega)
            ENDDO
          ENDIF
        ENDDO
      ENDDO
      IF (positive_matsu) THEN
        ! HM: We multiply it by two to account for the contribution
        ! from the negative Matsubara frequencies.
        dFE = dFE * 2.0d0
      ENDIF
      ! collect contributions from all pools
      CALL mp_sum(dFE, inter_pool_comm)
      CALL mp_barrier(inter_pool_comm)
      !
      dFE = dFE * pi * gtemp(itemp)
      !
      WRITE(stdout, '(5x, a, i3, a, f8.3, a, a, f12.6, a)') &
              'Temp (itemp = ', itemp, ') = ', gtemp(itemp) / kelvin2eV, ' K', &
              '  Free energy = ', dFE * 1d3, ' meV'
      WRITE(stdout, '(a)') ' '
      !
      DEALLOCATE(inv_wsi, STAT = ierr)
      IF (ierr /= 0) CALL errore('sum_eliashberg_aniso_iaxis_wrapper', 'Error deallocating inv_wsi', 1)
      DEALLOCATE(naznormi, STAT = ierr)
      IF (ierr /= 0) CALL errore('sum_eliashberg_aniso_iaxis_wrapper', 'Error deallocating naznormi', 1)
      !
      ! remove memory allocated for inv_wsi, deltaold, naznormi
      imelt = (2 + nbndfs * nks) * nsiw(itemp)
      CALL mem_size_eliashberg(2, -imelt)
      !
      IF (.NOT. limag_fly) THEN
        !
        DEALLOCATE(akeri, STAT = ierr)
        IF (ierr /= 0) CALL errore('sum_eliashberg_aniso_iaxis_wrapper', 'Error deallocating akeri', 1)
        !
        ! remove memory allocated for akeri (SH: this is adjusted for sparse sampling)
        imelt = nks * MAXVAL(nqfs(:)) * nbndfs**2 * 2 * (wsn(nsiw(itemp)) + 1)
        CALL mem_size_eliashberg(2, -imelt)
        !
      ENDIF
      !
    ENDIF
    !
    RETURN
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE sum_eliashberg_aniso_iaxis_wrapper
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    SUBROUTINE count_states()
    !-----------------------------------------------------------------------
    !!
    !! This routine counts the number of states (n',k') for each state (n, k).
    !!
    !! num_js1: the count of states (n',k') within fsthick window for each state (n, k) within fsthick window
    !! num_js2: the count of states (n',k') outside fsthick window for each state (n, k) within fsthick window
    !! num_js3: the count of states (n',k') outside fsthick window for each state (n, k) outside fsthick window
    !!
    USE kinds,             ONLY : DP
    USE input,             ONLY : fsthick, icoulomb, emax_coulomb, emin_coulomb
    USE supercond_common,  ONLY : ixkqf, ixqfs, nkfs, nqfs, ekfs, nbndfs, &
                                  ef0, num_js1, num_js2, num_js3, nk1_cl, nk2_cl, &
                                  nk3_cl, ik_bz_to_ibz_cl, nbnd_offset, nbnd_cl, &
                                  nbndfs_all, ik_cl_to_fs, nkstot_cl, ek_cl, &
                                  ibnd_kfs_all_to_kfs
    USE parallelism,       ONLY : fkbounds, para_bounds
    !
    ! Local variables
    LOGICAL :: linsidei
    !! used to determine if the state (ibnd_cl, ik_cl) is outside the window.
    LOGICAL :: linsidei2
    !! used to determine if the state (ibnd_cl, ik_cl) is outside the window.
    LOGICAL :: linsidei3
    !! used to determine if the state (ibnd_cl, ik_cl) is outside the window.
    LOGICAL :: linsidej
    !! used to determine if the state (jbnd_cl, jk_cl) is outside the window.
    LOGICAL :: linsidej2
    !! used to determine if the state (jbnd_cl, jk_cl) is outside the window.
    LOGICAL :: linsidej3
    !! used to determine if the state (jbnd_cl, jk_cl) is outside the window.
    INTEGER :: ik
    !! Counter on k-points
    INTEGER :: lower_bnd, upper_bnd
    !! Lower/upper bound index after k paral
    INTEGER :: iq
    !! Counter on q-points for which k+sign*q is within the Fermi shell
    INTEGER :: iq0
    !! Index of iq on full q-mesh
    INTEGER :: ibnd
    !! Counter on bands at k
    INTEGER :: jbnd
    !! Counter on bands at k+q
    INTEGER :: ik_cl
    !! Counter on k-points
    INTEGER :: jk
    !! Counter on k-points
    INTEGER :: jk_cl
    !! Counter on k-points: jk_cl = ik_bz_to_ibz_cl(jk)
    INTEGER :: ibnd_cl
    !! Counter on bands at k
    INTEGER :: jbnd_cl
    !! Counter on bands at k'
    INTEGER :: ibndfs
    !! Counter on bands in the fsthick window
    INTEGER :: jbndfs
    !! Counter on bands in the fsthick window
    INTEGER :: istate
    !! Counter on states, used when icoulomb > 0
    INTEGER :: is_start, is_stop
    !! Lower/upper bound index after paral for states
    INTEGER :: num_states
    !! Number of states per pool: num_states = is_stop - is_start + 1
    !
    CALL fkbounds(nkfs, lower_bnd, upper_bnd)
    !
    ! HM: Because nkstot_cl can be quite small,
    !     we do a combined parallelization on band index and k point index.
    CALL para_bounds(is_start, is_stop, nbnd_cl * nkstot_cl)
    num_states = is_stop - is_start + 1
    !
    ! HM: If sparse-ir sampling is employed,
    !     count the number of states (n',k') inside
    !     the fsthick window for each state (n, k).
    num_js1(:, :) = 0
    IF (icoulomb > 0) num_js2(:, :) = 0
    DO ik = lower_bnd, upper_bnd
      DO ibnd = 1, nbndfs
        IF (ABS(ekfs(ibnd, ik) - ef0) < fsthick) THEN
          DO iq = 1, nqfs(ik)
            ! iq0 - index of q-point on the full q-mesh
            iq0 = ixqfs(ik, iq)
            DO jbnd = 1, nbndfs
              IF (ABS(ekfs(jbnd, ixkqf(ik, iq0)) - ef0) < fsthick) THEN
                !
                num_js1(ibnd, ik) = num_js1(ibnd, ik) + 1
                !
              ENDIF
            ENDDO
          ENDDO
          IF (icoulomb > 0) THEN
            linsidej = .TRUE.
            linsidej2 = .TRUE.
            linsidej3 = .TRUE.
            !nn = 0 ! DEBUG
            DO jk = 1, (nk1_cl * nk2_cl * nk3_cl)
              ! find the index of the irreducible k point corresponding to the current k point
              jk_cl = ik_bz_to_ibz_cl(jk)
              !
              ! initilize linsidej
              linsidej = .TRUE.
              IF (ik_cl_to_fs(jk_cl) == 0) linsidej = .FALSE.
              !IF (.NOT.(linsidej)) CYCLE ! DEBUG
              DO jbnd_cl = 1, nbnd_cl
                ! initilize linsidej2
                linsidej2 = .TRUE.
                IF (((jbnd_cl - nbnd_offset) < 1) .OR. ((jbnd_cl - nbnd_offset) > nbndfs_all)) &
                  linsidej2 = .FALSE.
                !
                ! initilize linsidej3
                linsidej3 = .TRUE.
                IF (linsidej .AND. linsidej2) THEN
                  jbndfs = ibnd_kfs_all_to_kfs(jbnd_cl - nbnd_offset, ik_cl_to_fs(jk_cl))
                  IF ((jbndfs < 1) .OR. (jbndfs > nbndfs)) THEN
                    linsidej3 = .FALSE.
                  ELSE
                    ! Always use ekfs instead of ek_cl to determine
                    ! if the state (jbnd_cl, jk_cl) is outside the window.
                    IF (ABS(ekfs(jbndfs, ik_cl_to_fs(jk_cl)) - ef0) >= fsthick) &
                      linsidej3 = .FALSE.
                  ENDIF
                ENDIF
                !
                !IF (.NOT.(linsidej .AND. linsidej2 .AND. linsidej3)) CYCLE ! DEBUG
                !nn = nn + 1 ! DEBUG
                IF (.NOT.(linsidej .AND. linsidej2 .AND. linsidej3)) THEN
                  IF (((ek_cl(jbnd_cl, jk_cl) - ef0) >= emax_coulomb) .OR. &
                      ((ek_cl(jbnd_cl, jk_cl) - ef0) <= emin_coulomb)) THEN
                    CYCLE
                  ENDIF
                ENDIF
                !
                num_js2(ibnd, ik) = num_js2(ibnd, ik) + 1
                !
              ENDDO
            ENDDO
          ENDIF
        ENDIF
      ENDDO
    ENDDO
    !
    IF (icoulomb > 0) THEN
      ! HM: Count the number of states (n',k') between
      !     emin_coulomb and emax_coulomb for each istate corresponding to (n, k).
      !     Note that if a state of istate lies in the fsthick window, no need to
      !     count the number because adeltai_cl is will not calculated for
      !     the corresponding istate.
      !     In other words, what we have to count here is the number of final states
      !     inside of [emin_coulomb: emax_coulomb] for each initial state being
      !     outside of the fsthick window and inside of [emin_coulomb: emax_coulomb].
      !
      IF (num_states > 0) THEN
        num_js3(:) = 0
        DO istate = is_start, is_stop
          ik_cl = (istate - 1) / nbnd_cl + 1
          ibnd_cl = MOD(istate - 1, nbnd_cl) + 1
          !
          ! initilize linsidei
          linsidei = .TRUE.
          IF (ik_cl_to_fs(ik_cl) == 0) linsidei = .FALSE.
          ! initilize linsidei2
          linsidei2 = .TRUE.
          IF (((ibnd_cl - nbnd_offset) < 1) .OR. ((ibnd_cl - nbnd_offset) > nbndfs_all)) &
            linsidei2 = .FALSE.
          !
          ! initilize linsidei3
          linsidei3 = .TRUE.
          IF (linsidei .AND. linsidei2) THEN
            ibndfs = ibnd_kfs_all_to_kfs(ibnd_cl - nbnd_offset, ik_cl_to_fs(ik_cl))
            IF ((ibndfs < 1) .OR. (ibndfs > nbndfs)) THEN
              linsidei3 = .FALSE.
            ELSE
              ! Always use ekfs instead of ek_cl to determine
              ! if the state (ibnd_cl, ik_cl) is outside the window.
              IF (ABS(ekfs(ibndfs, ik_cl_to_fs(ik_cl)) - ef0) >= fsthick) &
                linsidei3 = .FALSE.
            ENDIF
          ENDIF
          !
          IF (linsidei .AND. linsidei2 .AND. linsidei3) THEN
            CYCLE
          ELSEIF ((ek_cl(ibnd_cl, ik_cl) - ef0) >= emax_coulomb .OR. &
                  (ek_cl(ibnd_cl, ik_cl) - ef0) <= emin_coulomb) THEN
            CYCLE
          ELSE
            DO jk = 1, (nk1_cl * nk2_cl * nk3_cl)
              ! find the index of the irreducible k point corresponding to the current k point
              jk_cl = ik_bz_to_ibz_cl(jk)
              !
              ! initilize linsidej
              linsidej = .TRUE.
              IF (ik_cl_to_fs(jk_cl) == 0) linsidej = .FALSE.
              DO jbnd_cl = 1, nbnd_cl
                ! initilize linsidej2
                linsidej2 = .TRUE.
                IF (((jbnd_cl - nbnd_offset) < 1) .OR. ((jbnd_cl - nbnd_offset) > nbndfs_all)) &
                  linsidej2 = .FALSE.
                !
                ! initilize linsidej3
                linsidej3 = .TRUE.
                IF (linsidej .AND. linsidej2) THEN
                  jbndfs = ibnd_kfs_all_to_kfs(jbnd_cl - nbnd_offset, ik_cl_to_fs(jk_cl))
                  IF ((jbndfs < 1) .OR. (jbndfs > nbndfs)) THEN
                    linsidej3 = .FALSE.
                  ELSE
                    ! Always use ekfs instead of ek_cl to determine
                    ! if the state (jbnd_cl, jk_cl) is outside the window.
                    IF (ABS(ekfs(jbndfs, ik_cl_to_fs(jk_cl)) - ef0) >= fsthick) &
                      linsidej3 = .FALSE.
                  ENDIF
                ENDIF
                !
                IF (.NOT.(linsidej .AND. linsidej2 .AND. linsidej3)) THEN
                  IF (((ek_cl(jbnd_cl, jk_cl) - ef0) >= emax_coulomb) .OR. &
                      ((ek_cl(jbnd_cl, jk_cl) - ef0) <= emin_coulomb)) THEN
                    CYCLE
                  ENDIF
                ENDIF
                !
                num_js3(istate) = num_js3(istate) + 1
              ENDDO
            ENDDO
          ENDIF ! linsidei .AND. linsidei2
        ENDDO
      ENDIF
    ENDIF
    !
    RETURN
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE count_states
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    SUBROUTINE sum_eliashberg_aniso_iaxis_fsr_simple(itemp)
    !-----------------------------------------------------------------------
    !!
    !! This routine solves the anisotropic FSR Eliashberg equations on the imaginary-axis
    !!
    USE kinds,             ONLY : DP
    USE global_var,        ONLY : wqf
    USE input,             ONLY : muc, fsthick
    USE supercond_common,  ONLY : ixkqf, ixqfs, nqfs, w0g, ekfs, nkfs, nbndfs, spin_fac, &
                                  dosef, ef0, nsiw, wsn, wsi, akeri, limag_fly, adeltai, &
                                  adeltaip, aznormi, naznormi
    USE parallelism,       ONLY : fkbounds
    USE low_lvl,           ONLY : mem_size_eliashberg
    USE ep_constants,      ONLY : pi, zero, one
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: itemp
    !! Counter on temperature index
    !
    ! Local variables
    INTEGER :: iw, iwp
    !! Counter on frequency imag-axis
    INTEGER :: ik
    !! Counter on k-points
    INTEGER :: iq
    !! Counter on q-points for which k+sign*q is within the Fermi shell
    INTEGER :: iq0
    !! Index of iq on full q-mesh
    INTEGER :: ibnd
    !! Counter on bands at k
    INTEGER :: jbnd
    !! Counter on bands at k+q
    INTEGER :: lower_bnd, upper_bnd
    !! Lower/upper bound index after k paral
    INTEGER(8) :: imelt
    !! Counter memory
    INTEGER :: ierr
    !! Error status
    !
    REAL(KIND = DP) :: lambdam
    !! K_{-}(n,n',T))
    REAL(KIND = DP) :: lambdap
    !! K_{+}(n,n',T)
    REAL(KIND = DP) :: kernelm
    !! kernelm = lambdam - lambdap
    REAL(KIND = DP) :: kernelp
    !! kernelp = lambdam + lambdap
    REAL(KIND = DP) :: weight
    !! Factor in supercond. equations
    REAL(KIND = DP) :: esqrt, zesqrt, desqrt
    !! Temporary variables
    REAL(KIND = DP) :: inv_dos
    !! Invese dos inv_dos = 1/dosef. Defined for efficiency reason
    REAL(KIND = DP), ALLOCATABLE :: inv_wsi(:)
    !! Invese imaginary freq. inv_wsi = 1/wsi. Defined for efficiency reason
    !
    inv_dos = one / dosef
    !
    CALL fkbounds(nkfs, lower_bnd, upper_bnd)
    !
    ! get the size of required memory for inv_wsi
    imelt = nsiw(itemp)
    CALL mem_size_eliashberg(2, imelt)
    ALLOCATE(inv_wsi(nsiw(itemp)), STAT = ierr)
    IF (ierr /= 0) CALL errore('sum_eliashberg_aniso_iaxis_fsr_simple', 'Error allocating inv_wsi', 1)
    !
    DO iw = 1, nsiw(itemp)
      inv_wsi(iw) = one / wsi(iw)
    ENDDO
    !
    naznormi(:, :, :) = zero
    adeltai(:, :, :)  = zero
    aznormi(:, :, :)  = zero
    DO ik = lower_bnd, upper_bnd
      DO ibnd = 1, nbndfs
        IF (ABS(ekfs(ibnd, ik) - ef0) < fsthick) THEN
          DO iq = 1, nqfs(ik)
            ! iq0 - index of q-point on the full q-mesh
            iq0 = ixqfs(ik, iq)
            DO jbnd = 1, nbndfs
              IF (ABS(ekfs(jbnd, ixkqf(ik, iq0)) - ef0) < fsthick) THEN
                weight = wqf(iq) * w0g(jbnd, ixkqf(ik, iq0)) * inv_dos
                DO iwp = 1, nsiw(itemp) ! loop over omega_prime
                  esqrt  = weight / DSQRT(wsi(iwp)**2.d0 + adeltaip(iwp, jbnd, ixkqf(ik, iq0))**2.d0)
                  zesqrt = esqrt * wsi(iwp)
                  desqrt = esqrt * adeltaip(iwp, jbnd, ixkqf(ik, iq0))
                  !
                  DO iw = 1, nsiw(itemp) ! loop over omega
                    IF (limag_fly) THEN
                      CALL lambdar_aniso_ver1(ik, iq, ibnd, jbnd, wsi(iw) - wsi(iwp), lambdam)
                      CALL lambdar_aniso_ver1(ik, iq, ibnd, jbnd, wsi(iw) + wsi(iwp), lambdap)
                    ELSE
                      ! SH: For general case (including sparse sampling)
                      !       "actual" Matsubara indices n1/n2 are needed instead of iw/iwp
                      lambdam = akeri(ABS(wsn(iw) - wsn(iwp)) + 1,     jbnd, iq, ibnd, ik)
                      lambdap = akeri(ABS(wsn(iw) + wsn(iwp) + 1) + 1, jbnd, iq, ibnd, ik)
                    ENDIF
                    ! Eq. (4.4) in Picket, PRB 26, 1186 (1982)
                    kernelm = lambdam - lambdap
                    kernelp = lambdam + lambdap
                    naznormi(iw, ibnd, ik) = naznormi(iw, ibnd, ik) + weight * kernelm
                    ! Eqs.(21)-(22) in Margine and Giustino, PRB 87, 024505 (2013)
                    ! using kernelm and kernelp the sum over |wp| < wscut in Eqs. (21)-(22)
                    ! is rewritten as a sum over iwp = 1, nsiw(itemp)
                    aznormi(iw, ibnd, ik) = aznormi(iw, ibnd, ik) + zesqrt * kernelm
                    adeltai(iw, ibnd, ik) = adeltai(iw, ibnd, ik) + desqrt * (kernelp - 2.d0 * muc * spin_fac)
                  ENDDO ! iw
                ENDDO ! iwp
              ENDIF
            ENDDO ! jbnd
          ENDDO ! iq
        ENDIF
      ENDDO ! ibnd
    ENDDO ! ik
    !
    DEALLOCATE(inv_wsi, STAT = ierr)
    IF (ierr /= 0) CALL errore('sum_eliashberg_aniso_iaxis_fsr_simple', 'Error deallocating inv_wsi', 1)
    ! remove memory allocated for inv_wsi
    imelt = nsiw(itemp)
    CALL mem_size_eliashberg(2, -imelt)
    !
    RETURN
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE sum_eliashberg_aniso_iaxis_fsr_simple
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    SUBROUTINE sum_eliashberg_aniso_iaxis_fsr_fft(itemp, narray)
    !-----------------------------------------------------------------------
    !!
    !! This routine solves the anisotropic FSR Eliashberg equations on the imaginary-axis
    !! using FFT.
    !! This is used for uniform sampling (gridsamp == 3).
    !!
    USE kinds,             ONLY : DP
    USE global_var,        ONLY : wqf, gtemp
    USE input,             ONLY : muc, fsthick, muchem, positive_matsu
    USE supercond_common,  ONLY : ixkqf, ixqfs, nqfs, w0g, ekfs, nkfs, nbndfs, spin_fac, &
                                  dosef, ef0, nsiw, wsi, akeri, limag_fly, adeltai, &
                                  adeltaip, aznormi, aznormip, naznormi, fft_in1, &
                                  fft_out1, fft_in2, fft_out2
    USE parallelism,       ONLY : fkbounds
    USE low_lvl,           ONLY : mem_size_eliashberg
    USE ep_constants,      ONLY : pi, zero, one, czero, cone
    USE fft_scalar,        ONLY : cft_1z
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: itemp
    !! Counter on temperature index
    INTEGER, INTENT(in) :: narray(2)
    !! integers to determine the size of arrays of FFT
    !
    ! Local variables
    INTEGER :: iw
    !! Counter on frequency imag-axis
    INTEGER :: iw_f
    !! Counter for Fermionic variables
    INTEGER :: iw_b
    !! Counter for akeri
    INTEGER :: wsnx
    !! integer part of Matsubara freq.
    INTEGER :: ik
    !! Counter on k-points
    INTEGER :: iq
    !! Counter on q-points for which k+sign*q is within the Fermi shell
    INTEGER :: iq0
    !! Index of iq on full q-mesh
    INTEGER :: ibnd
    !! Counter on bands at k
    INTEGER :: jbnd
    !! Counter on bands at k+q
    INTEGER :: lower_bnd, upper_bnd
    !! Lower/upper bound index after k paral
    INTEGER :: n
    !! n = nsiw(itemp)
    INTEGER(8) :: imelt
    !! Counter memory
    INTEGER :: ierr
    !! Error status
    !
    REAL(KIND = DP) :: kernel
    !! K_{-}(n,n',T))
    REAL(KIND = DP) :: weight
    !! Factor in supercond. equations
    REAL(KIND = DP) :: esqrt
    !! Temporary variables
    REAL(KIND = DP) :: esqrt0
    !! esqrt and zesqrt calculated with adeltaip = 0.0, which is used for naznormi
    REAL(KIND = DP) :: inv_dos
    !! Invese dos inv_dos = 1/dosef. Defined for efficiency reason
    REAL(KIND = DP) :: omega
    !! Matsubara frequency
    REAL(KIND = DP), ALLOCATABLE :: inv_wsi(:)
    !! Invese imaginary freq. inv_wsi = 1/wsi. Defined for efficiency reason
    !
    inv_dos = one / dosef
    !
    CALL fkbounds(nkfs, lower_bnd, upper_bnd)
    !
    ! get the size of required memory for inv_wsi
    imelt = nsiw(itemp)
    CALL mem_size_eliashberg(2, imelt)
    ALLOCATE(inv_wsi(nsiw(itemp)), STAT = ierr)
    IF (ierr /= 0) CALL errore('sum_eliashberg_aniso_iaxis_fsr_fft', 'Error allocating inv_wsi', 1)
    !
    DO iw = 1, nsiw(itemp)
      inv_wsi(iw) = one / wsi(iw)
    ENDDO
    !
    n = 2 * nsiw(itemp)
    !
    naznormi(:, :, :) = zero
    adeltai(:, :, :)  = zero
    aznormi(:, :, :)  = zero
    DO ik = lower_bnd, upper_bnd
      DO ibnd = 1, nbndfs
        IF (ABS(ekfs(ibnd, ik) - ef0) < fsthick) THEN
          DO iq = 1, nqfs(ik)
            ! iq0 - index of q-point on the full q-mesh
            iq0 = ixqfs(ik, iq)
            DO jbnd = 1, nbndfs
              IF (ABS(ekfs(jbnd, ixkqf(ik, iq0)) - ef0) < fsthick) THEN
                !
                fft_in1(:) = czero
                fft_out1(:) = czero
                fft_in2(:) = czero
                fft_out2(:) = czero
                !
                DO iw = 1, n
                  !
                  ! FOR Fermionic functions
                  !
                  !!   iw =   1,   2,   3,   4,   5,...
                  !! wsnx =   1,   3,   5,   7,   9,...
                  !! iw_f =   1,   2,   3,   4,   5,...
                  !!
                  !!   iw = ...,  N-1,    N,   N+1,   N+2,   N+3,...
                  !! wsnx = ..., 2N-3, 2N-1, -2N+1, -2N+3, -2N+5,...
                  !! iw_f = ...,  N-1,    N,     N,   N-1,   N-2,...
                  !!
                  !!   iw = ..., 2N-3, 2N-2, 2N-1,   2N
                  !! wsnx = ...,   -7,   -5,   -3,   -1
                  !! iw_f = ...,    4,    3,    2,    1
                  !!
                  !! where N = nsiw(itemp) = (wsn(nsiw(itemp)) + 1)
                  !! NOTE: 2N = n = 2 * nsiw(itemp)
                  !!
                  ! wsnx: integer part of Matsubara freq.
                  ! omega: Matsubara frequencies in the unit of eV
                  !
                  IF (iw .LE. n/2) THEN
                    wsnx = 2 * iw - 1
                    iw_f = iw
                  ELSEIF (iw .GT. n/2) THEN
                    wsnx = (2 * (iw - n/2) - 1) - n ! (2 * (iw - N) - 1) - 2N
                    iw_f = n - iw + 1 ! 2N - iw + 1
                  ENDIF
                  !
                  omega = DBLE(wsnx) * pi * gtemp(itemp)
                  !
                  weight = w0g(jbnd, ixkqf(ik, iq0))
                  esqrt  = weight / DSQRT(omega**2.d0 + adeltaip(iw_f, jbnd, ixkqf(ik, iq0))**2.d0)
                  fft_in1(0 * n + iw) = cone * esqrt * omega
                  fft_in1(1 * n + iw) = cone * esqrt * adeltaip(iw_f, jbnd, ixkqf(ik, iq0))
                  esqrt0  = weight / DSQRT(omega**2.d0)
                  fft_in1(2 * n + iw) = cone * esqrt0 * omega
                ENDDO
                !
                DO iw = 1, n
                  !
                  ! FOR Bosonic functions
                  !
                  !!   iw = 1, 2, 3, 4, 5,...
                  !! wsnx = 0, 2, 4, 6, 8,...
                  !! iw_b = 1, 2, 3, 4, 5,...
                  !!
                  !!   iw = ...,  N-1,    N, N+1,   N+2,   N+3,...
                  !! wsnx = ..., 2N-4, 2N-2, -2N, -2N+2, -2N+4,...
                  !! iw_b = ...,  N-1,    N, N+1,     N,   N-1,...
                  !!
                  !!   iw = ..., 2N-3, 2N-2, 2N-1,   2N
                  !! wsnx = ...,   -8,   -6,   -4,   -2
                  !! iw_b = ...,    5,    4,    3,    2
                  !!
                  ! In this case, only akeri(1 : n/2 + 1) are used.
                  !
                  ! wsnx: integer part of Matsubara freq.
                  ! omega: Matsubara frequencies in the unit of eV
                  !
                  IF (iw .LE. n/2) THEN
                    wsnx = 2 * (iw - 1)
                    iw_b = iw
                  ELSEIF (iw .GT. n/2) THEN
                    wsnx = 2 * (iw - n - 1) ! 2 * (iw - 2N - 1)
                    iw_b = ABS(iw - n - 1) + 1 ! ABS(iw - 2N - 1) + 1
                  ENDIF
                  !
                  omega = DBLE(wsnx) * pi * gtemp(itemp)
                  !
                  IF (limag_fly) THEN
                    CALL lambdar_aniso_ver1(ik, iq, ibnd, jbnd, omega, kernel)
                  ELSE
                    kernel = akeri(iw_b, jbnd, iq, ibnd, ik)
                  ENDIF
                  fft_in1(3 * n + iw) = cone * kernel * inv_dos
                  fft_in1(4 * n + iw) = cone * (kernel - muc * spin_fac) * inv_dos
                ENDDO
                !
                CALL cft_1z(fft_in1, narray(1), n, n, -1, fft_out1)
                ! divided by the factor n in cft_1z
                !
                fft_in2(0 * n + 1 : 1 * n) = &
                fft_out1(0 * n + 1 : 1 * n) * fft_out1(3 * n + 1 : 4 * n) * REAL(n, KIND=DP)
                fft_in2(1 * n + 1 : 2 * n) = &
                fft_out1(1 * n + 1 : 2 * n) * fft_out1(4 * n + 1 : 5 * n) * REAL(n, KIND=DP)
                fft_in2(2 * n + 1 : 3 * n) = &
                fft_out1(2 * n + 1 : 3 * n) * fft_out1(3 * n + 1 : 4 * n) * REAL(n, KIND=DP)
                !
                CALL cft_1z(fft_in2, narray(2), n, n, +1, fft_out2)
                !
                DO iw = 1, n/2
                  !
                  iw_f = iw
                  !
                  weight = wqf(iq)
                  aznormi(iw_f, ibnd, ik) = &
                  aznormi(iw_f, ibnd, ik) + weight * REAL(fft_out2(0 * n + iw), KIND=DP)
                  adeltai(iw_f, ibnd, ik) = &
                  adeltai(iw_f, ibnd, ik) + weight * REAL(fft_out2(1 * n + iw), KIND=DP)
                  naznormi(iw_f, ibnd, ik) = &
                  naznormi(iw_f, ibnd, ik) + weight * REAL(fft_out2(2 * n + iw), KIND=DP)
                ENDDO
                !
              ENDIF
            ENDDO
          ENDDO
        ENDIF
      ENDDO
    ENDDO
    !
    DEALLOCATE(inv_wsi, STAT = ierr)
    IF (ierr /= 0) CALL errore('sum_eliashberg_aniso_iaxis_fsr_fft', 'Error deallocating inv_wsi', 1)
    ! remove memory allocated for inv_wsi
    imelt = nsiw(itemp)
    CALL mem_size_eliashberg(2, -imelt)
    !
    RETURN
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE sum_eliashberg_aniso_iaxis_fsr_fft
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    SUBROUTINE sum_eliashberg_aniso_iaxis_fbw_simple(itemp, nel, nstate)
    !-----------------------------------------------------------------------
    !!
    !! This routine solves the anisotropic FBW Eliashberg equations on the imaginary-axis
    !! using simple summation.
    !! This is used for uniform sampling and sparse sampling (gridsamp <= 1).
    !!
    USE kinds,             ONLY : DP
    USE global_var,        ONLY : wqf
    USE input,             ONLY : muc, fsthick, muchem, positive_matsu
    USE supercond_common,  ONLY : ixkqf, ixqfs, nqfs, ekfs, nkfs, nbndfs, spin_fac, &
                                  dosef, ef0, nsiw, wsn, wsi, akeri, limag_fly, adeltai, &
                                  adeltaip, aznormi, aznormip, naznormi, ashifti, &
                                  ashiftip, muintr
    USE parallelism,       ONLY : fkbounds
    USE low_lvl,           ONLY : mem_size_eliashberg
    USE ep_constants,      ONLY : pi, zero, one
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: itemp
    !! Counter on temperature index
    REAL(KIND = DP), INTENT(in) :: nel
    !! mu_inter parameters
    REAL(KIND = DP), INTENT(in) :: nstate
    !! mu_inter parameters
    !
    ! Local variables
    INTEGER :: iw, iwp
    !! Counter on frequency imag-axis
    INTEGER :: ik
    !! Counter on k-points
    INTEGER :: iq
    !! Counter on q-points for which k+sign*q is within the Fermi shell
    INTEGER :: iq0
    !! Index of iq on full q-mesh
    INTEGER :: ibnd
    !! Counter on bands at k
    INTEGER :: jbnd
    !! Counter on bands at k+q
    INTEGER :: lower_bnd, upper_bnd
    !! Lower/upper bound index after k paral
    INTEGER(8) :: imelt
    !! Counter memory
    INTEGER :: ierr
    !! Error status
    !
    REAL(KIND = DP) :: lambdam
    !! K_{-}(n,n',T))
    REAL(KIND = DP) :: lambdap
    !! K_{+}(n,n',T)
    REAL(KIND = DP) :: kernelm
    !! kernelm = lambdam - lambdap
    REAL(KIND = DP) :: kernelp
    !! kernelp = lambdam + lambdap
    REAL(KIND = DP) :: kernel
    !! K_{-}(n,n',T))
    REAL(KIND = DP) :: weight
    !! Factor in supercond. equations
    REAL(KIND = DP) :: esqrt, zesqrt, desqrt, sesqrt
    !! Temporary variables
    REAL(KIND = DP) :: esqrt0, zesqrt0
    !! esqrt and zesqrt calculated with adeltaip = 0.0, which is used for naznormi
    REAL(KIND = DP) :: inv_dos
    !! Invese dos inv_dos = 1/dosef. Defined for efficiency reason
    REAL(KIND = DP), ALLOCATABLE :: inv_wsi(:)
    !! Invese imaginary freq. inv_wsi = 1/wsi. Defined for efficiency reason
    !
    inv_dos = one / dosef
    !
    CALL fkbounds(nkfs, lower_bnd, upper_bnd)
    !
    ! get the size of required memory for inv_wsi
    imelt = nsiw(itemp)
    CALL mem_size_eliashberg(2, imelt)
    ALLOCATE(inv_wsi(nsiw(itemp)), STAT = ierr)
    IF (ierr /= 0) CALL errore('sum_eliashberg_aniso_iaxis_fbw_simple', 'Error allocating inv_wsi', 1)
    !
    DO iw = 1, nsiw(itemp)
      inv_wsi(iw) = one / wsi(iw)
    ENDDO
    !
    ! SH: update the chemical potential from the inital guess
    IF (muchem) CALL mu_inter_aniso(itemp, muintr, nel, nstate)
    !
    naznormi(:, :, :) = zero
    adeltai(:, :, :)  = zero
    aznormi(:, :, :)  = zero
    ashifti(:, :, :)  = zero
    DO ik = lower_bnd, upper_bnd
      DO ibnd = 1, nbndfs
        IF (ABS(ekfs(ibnd, ik) - ef0) < fsthick) THEN
          DO iq = 1, nqfs(ik)
            ! iq0 - index of q-point on the full q-mesh
            iq0 = ixqfs(ik, iq)
            DO jbnd = 1, nbndfs
              IF (ABS(ekfs(jbnd, ixkqf(ik, iq0)) - ef0) < fsthick) THEN
                !! this is for FBW case
                weight = wqf(iq) * inv_dos
                DO iwp = 1, nsiw(itemp) ! loop over omega_prime
                  esqrt = weight / ((wsi(iwp) * aznormip(iwp, jbnd, ixkqf(ik, iq0)))**2.d0 &
                            + (ekfs(jbnd, ixkqf(ik, iq0)) - muintr + ashiftip(iwp, jbnd, ixkqf(ik, iq0)))**2.d0 &
                            + (adeltaip(iwp, jbnd, ixkqf(ik, iq0)) * aznormip(iwp, jbnd, ixkqf(ik, iq0)))**2.d0)
                  zesqrt = esqrt * wsi(iwp) * aznormip(iwp, jbnd, ixkqf(ik, iq0))
                  desqrt = esqrt * adeltaip(iwp, jbnd, ixkqf(ik, iq0)) * aznormip(iwp, jbnd, ixkqf(ik, iq0))
                  sesqrt = esqrt * (ekfs(jbnd, ixkqf(ik, iq0)) - muintr + ashiftip(iwp, jbnd, ixkqf(ik, iq0)))
                  esqrt0 = weight / ((wsi(iwp) * aznormip(iwp, jbnd, ixkqf(ik, iq0)))**2.d0 &
                            + (ekfs(jbnd, ixkqf(ik, iq0)) - muintr + ashiftip(iwp, jbnd, ixkqf(ik, iq0)))**2.d0)
                  zesqrt0 = esqrt0 * wsi(iwp) * aznormip(iwp, jbnd, ixkqf(ik, iq0))
                  DO iw = 1, nsiw(itemp) ! loop over omega
                    IF (positive_matsu) THEN
                      IF (limag_fly) THEN
                        CALL lambdar_aniso_ver1(ik, iq, ibnd, jbnd, wsi(iw) - wsi(iwp), lambdam)
                        CALL lambdar_aniso_ver1(ik, iq, ibnd, jbnd, wsi(iw) + wsi(iwp), lambdap)
                      ELSE
                        ! SH: For general case (including sparse sampling)
                        !       "actual" Matsubara indices n1/n2 are needed instead of iw/iwp
                        lambdam = akeri(ABS(wsn(iw) - wsn(iwp)) + 1,     jbnd, iq, ibnd, ik)
                        lambdap = akeri(ABS(wsn(iw) + wsn(iwp) + 1) + 1, jbnd, iq, ibnd, ik)
                      ENDIF
                      ! Eq. (4.4) in Picket, PRB 26, 1186 (1982)
                      kernelm = lambdam - lambdap
                      kernelp = lambdam + lambdap
                      naznormi(iw, ibnd, ik) = naznormi(iw, ibnd, ik) + zesqrt0 * kernelm
                      ! Eqs.(15-17) in Margine and Giustino, PRB 87, 024505 (2013)
                      ! using kernelm and kernelp the sum over |wp| < wscut
                      ! is rewritten as a sum over iwp = 1, nsiw(itemp)
                      aznormi(iw, ibnd, ik) = aznormi(iw, ibnd, ik) + zesqrt * kernelm
                      adeltai(iw, ibnd, ik) = adeltai(iw, ibnd, ik) + desqrt * (kernelp - 2.d0 * muc * spin_fac)
                      ashifti(iw, ibnd, ik) = ashifti(iw, ibnd, ik) + sesqrt * kernelp
                    ELSE
                      IF (limag_fly) THEN
                        CALL lambdar_aniso_ver1(ik, iq, ibnd, jbnd, wsi(iw) - wsi(iwp), kernel)
                      ELSE
                        ! SH: For general case (including sparse sampling)
                        !       "actual" Matsubara indices n1/n2 are needed instead of iw/iwp
                        kernel = akeri(ABS(wsn(iw) - wsn(iwp)) + 1,     jbnd, iq, ibnd, ik)
                      ENDIF
                      !
                      naznormi(iw, ibnd, ik) = naznormi(iw, ibnd, ik) + zesqrt0 * kernel
                      ! Eqs.(15-17) in Margine and Giustino, PRB 87, 024505 (2013)
                      aznormi(iw, ibnd, ik) = aznormi(iw, ibnd, ik) + zesqrt * kernel
                      adeltai(iw, ibnd, ik) = adeltai(iw, ibnd, ik) + desqrt * (kernel - muc * spin_fac)
                      ashifti(iw, ibnd, ik) = ashifti(iw, ibnd, ik) + sesqrt * kernel
                    ENDIF
                  ENDDO ! iw
                ENDDO ! iwp
              ENDIF
            ENDDO ! jbnd
          ENDDO ! iq
        ENDIF
      ENDDO ! ibnd
    ENDDO ! ik
    !
    DEALLOCATE(inv_wsi, STAT = ierr)
    IF (ierr /= 0) CALL errore('sum_eliashberg_aniso_iaxis_fbw_simple', 'Error deallocating inv_wsi', 1)
    ! remove memory allocated for inv_wsi
    imelt = nsiw(itemp)
    CALL mem_size_eliashberg(2, -imelt)
    !
    RETURN
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE sum_eliashberg_aniso_iaxis_fbw_simple
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    SUBROUTINE sum_eliashberg_aniso_iaxis_fbw_fft(itemp, nel, nstate, narray)
    !-----------------------------------------------------------------------
    !!
    !! This routine solves the anisotropic FBW Eliashberg equations on the imaginary-axis
    !! using FFT.
    !! This is used for uniform sampling (gridsamp == 3).
    !!
    USE kinds,             ONLY : DP
    USE global_var,        ONLY : wqf, gtemp
    USE input,             ONLY : muc, fsthick, muchem, positive_matsu
    USE supercond_common,  ONLY : ixkqf, ixqfs, nqfs, ekfs, nkfs, nbndfs, spin_fac, &
                                  dosef, ef0, nsiw, wsi, akeri, limag_fly, adeltai, &
                                  adeltaip, aznormi, aznormip, naznormi, ashifti, &
                                  ashiftip, muintr, fft_in1, fft_out1, &
                                  fft_in2, fft_out2
    USE parallelism,       ONLY : fkbounds
    USE low_lvl,           ONLY : mem_size_eliashberg
    USE ep_constants,      ONLY : pi, zero, one, czero, cone
    USE fft_scalar,        ONLY : cft_1z
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: itemp
    !! Counter on temperature index
    INTEGER, INTENT(in) :: narray(2)
    !! integers to determine the size of arrays of FFT
    REAL(KIND = DP), INTENT(in) :: nel
    !! mu_inter parameters
    REAL(KIND = DP), INTENT(in) :: nstate
    !! mu_inter parameters
    !
    ! Local variables
    INTEGER :: iw
    !! Counter on frequency imag-axis
    INTEGER :: iw_f
    !! Counter for Fermionic variables
    INTEGER :: iw_b
    !! Counter for akeri
    INTEGER :: wsnx
    !! integer part of Matsubara freq.
    INTEGER :: ik
    !! Counter on k-points
    INTEGER :: iq
    !! Counter on q-points for which k+sign*q is within the Fermi shell
    INTEGER :: iq0
    !! Index of iq on full q-mesh
    INTEGER :: ibnd
    !! Counter on bands at k
    INTEGER :: jbnd
    !! Counter on bands at k+q
    INTEGER :: lower_bnd, upper_bnd
    !! Lower/upper bound index after k paral
    INTEGER :: n
    !! n = nsiw(itemp)
    INTEGER(8) :: imelt
    !! Counter memory
    INTEGER :: ierr
    !! Error status
    !
    REAL(KIND = DP) :: kernel
    !! K_{-}(n,n',T))
    REAL(KIND = DP) :: weight
    !! Factor in supercond. equations
    REAL(KIND = DP) :: esqrt
    !! Temporary variables
    REAL(KIND = DP) :: esqrt0
    !! esqrt and zesqrt calculated with adeltaip = 0.0, which is used for naznormi
    REAL(KIND = DP) :: inv_dos
    !! Invese dos inv_dos = 1/dosef. Defined for efficiency reason
    REAL(KIND = DP) :: omega
    !! Matsubara frequency
    REAL(KIND = DP), ALLOCATABLE :: inv_wsi(:)
    !! Invese imaginary freq. inv_wsi = 1/wsi. Defined for efficiency reason
    !
    inv_dos = one / dosef
    !
    CALL fkbounds(nkfs, lower_bnd, upper_bnd)
    !
    ! get the size of required memory for inv_wsi
    imelt = nsiw(itemp)
    CALL mem_size_eliashberg(2, imelt)
    ALLOCATE(inv_wsi(nsiw(itemp)), STAT = ierr)
    IF (ierr /= 0) CALL errore('sum_eliashberg_aniso_iaxis_fbw_fft', 'Error allocating inv_wsi', 1)
    !
    DO iw = 1, nsiw(itemp)
      inv_wsi(iw) = one / wsi(iw)
    ENDDO
    !
    IF (positive_matsu) THEN
      n = 2 * nsiw(itemp)
    ELSE
      n = nsiw(itemp)
    ENDIF
    ! SH: update the chemical potential from the inital guess
    IF (muchem) CALL mu_inter_aniso(itemp, muintr, nel, nstate)
    !
    naznormi(:, :, :) = zero
    adeltai(:, :, :)  = zero
    aznormi(:, :, :)  = zero
    ashifti(:, :, :)  = zero
    !
    DO ik = lower_bnd, upper_bnd
      DO ibnd = 1, nbndfs
        IF (ABS(ekfs(ibnd, ik) - ef0) < fsthick) THEN
          DO iq = 1, nqfs(ik)
            ! iq0 - index of q-point on the full q-mesh
            iq0 = ixqfs(ik, iq)
            DO jbnd = 1, nbndfs
              IF (ABS(ekfs(jbnd, ixkqf(ik, iq0)) - ef0) < fsthick) THEN
                !
                fft_in1(:) = czero
                fft_out1(:) = czero
                fft_in2(:) = czero
                fft_out2(:) = czero
                !
                DO iw = 1, n
                  !
                  IF (positive_matsu) THEN
                    ! FOR Fermionic functions
                    !
                    !!   iw =   1,   2,   3,   4,   5,...
                    !! wsnx =   1,   3,   5,   7,   9,...
                    !! iw_f =   1,   2,   3,   4,   5,...
                    !!
                    !!   iw = ...,  N-1,    N,   N+1,   N+2,   N+3,...
                    !! wsnx = ..., 2N-3, 2N-1, -2N+1, -2N+3, -2N+5,...
                    !! iw_f = ...,  N-1,    N,     N,   N-1,   N-2,...
                    !!
                    !!   iw = ..., 2N-3, 2N-2, 2N-1,   2N
                    !! wsnx = ...,   -7,   -5,   -3,   -1
                    !! iw_f = ...,    4,    3,    2,    1
                    !!
                    !! where N = nsiw(itemp) = (wsn(nsiw(itemp)) + 1)
                    !! NOTE: 2N = n = 2 * nsiw(itemp)
                    !!
                    ! wsnx: integer part of Matsubara freq.
                    ! omega: Matsubara frequencies in the unit of eV
                    !
                    IF (iw .LE. n/2) THEN
                      wsnx = 2 * iw - 1
                      iw_f = iw
                    ELSEIF (iw .GT. n/2) THEN
                      wsnx = (2 * (iw - n/2) - 1) - n ! (2 * (iw - N) - 1) - 2N
                      iw_f = n - iw + 1 ! 2N - iw + 1
                    ENDIF
                  ELSE
                    ! FOR Fermionic functions
                    !
                    !!   iw =   1,   2,   3,   4,   5,...
                    !! wsnx =   1,   3,   5,   7,   9,...
                    !! iw_f = N+1, N+2, N+3, N+4, N+5,...
                    !!
                    !!   iw = ...,  N-1,    N,   N+1,   N+2,   N+3,...
                    !! wsnx = ..., 2N-3, 2N-1, -2N+1, -2N+3, -2N+5,...
                    !! iw_f = ..., 2N-1,   2N,     1,     2,     3,...
                    !!
                    !!   iw = ..., 2N-3, 2N-2, 2N-1,   2N
                    !! wsnx = ...,   -7,   -5,   -3,   -1
                    !! iw_f = ...,  N-3,  N-2,  N-1,    N
                    !!
                    !! where N = nsiw(itemp)/2 = (wsn(nsiw(itemp)) + 1)
                    !! NOTE: 2N = n = nsiw(itemp)
                    !!
                    ! wsnx: integer part of Matsubara freq.
                    ! omega: Matsubara frequencies in the unit of eV
                    !
                    IF (iw .LE. n/2) THEN
                      wsnx = 2 * iw - 1
                      iw_f = iw + n/2 ! iw + N
                    ELSEIF (iw .GT. n/2) THEN
                      wsnx = 2 * (iw - n) - 1 ! 2 * (iw - 2N) - 1
                      iw_f = iw - n/2 ! iw - N
                    ENDIF
                  ENDIF
                  !
                  omega = DBLE(wsnx) * pi * gtemp(itemp)
                  !
                  esqrt = one / ((omega * aznormip(iw_f, jbnd, ixkqf(ik, iq0)))**2.d0 &
                                  + (ekfs(jbnd, ixkqf(ik, iq0)) - muintr + ashiftip(iw_f, jbnd, ixkqf(ik, iq0)))**2.d0 &
                                  + (adeltaip(iw_f, jbnd, ixkqf(ik, iq0)) * aznormip(iw_f, jbnd, ixkqf(ik, iq0)))**2.d0)
                  fft_in1(0 * n + iw) = &
                  cone * esqrt * omega * aznormip(iw_f, jbnd, ixkqf(ik, iq0))
                  fft_in1(1 * n + iw) = &
                  cone * esqrt * (ekfs(jbnd, ixkqf(ik, iq0)) - muintr + ashiftip(iw_f, jbnd, ixkqf(ik, iq0)))
                  fft_in1(2 * n + iw) = &
                  cone * esqrt * (adeltaip(iw_f, jbnd, ixkqf(ik, iq0)) * aznormip(iw_f, jbnd, ixkqf(ik, iq0)))
                  esqrt0 = one / ((omega * aznormip(iw_f, jbnd, ixkqf(ik, iq0)))**2.d0 &
                                  + (ekfs(jbnd, ixkqf(ik, iq0)) - muintr + ashiftip(iw_f, jbnd, ixkqf(ik, iq0)))**2.d0)
                  fft_in1(3 * n + iw) = &
                  cone * esqrt0 * omega * aznormip(iw_f, jbnd, ixkqf(ik, iq0))
                ENDDO
                !
                DO iw = 1, n
                  !
                  ! FOR Bosonic functions
                  !
                  !!   iw = 1, 2, 3, 4, 5,...
                  !! wsnx = 0, 2, 4, 6, 8,...
                  !! iw_b = 1, 2, 3, 4, 5,...
                  !!
                  !!   iw = ...,  N-1,    N, N+1,   N+2,   N+3,...
                  !! wsnx = ..., 2N-4, 2N-2, -2N, -2N+2, -2N+4,...
                  !! iw_b = ...,  N-1,    N, N+1,     N,   N-1,...
                  !!
                  !!   iw = ..., 2N-3, 2N-2, 2N-1,   2N
                  !! wsnx = ...,   -8,   -6,   -4,   -2
                  !! iw_b = ...,    5,    4,    3,    2
                  !!
                  ! In this case, only akeri(1 : n/2 + 1) are used.
                  !
                  ! wsnx: integer part of Matsubara freq.
                  ! omega: Matsubara frequencies in the unit of eV
                  !
                  IF (iw .LE. n/2) THEN
                    wsnx = 2 * (iw - 1)
                    iw_b = iw
                  ELSEIF (iw .GT. n/2) THEN
                    wsnx = 2 * (iw - n - 1) ! 2 * (iw - 2N - 1)
                    iw_b = ABS(iw - n - 1) + 1 ! ABS(iw - 2N - 1) + 1
                  ENDIF
                  !
                  omega = DBLE(wsnx) * pi * gtemp(itemp)
                  !
                  IF (limag_fly) THEN
                    CALL lambdar_aniso_ver1(ik, iq, ibnd, jbnd, omega, kernel)
                  ELSE
                    kernel = akeri(iw_b, jbnd, iq, ibnd, ik)
                  ENDIF
                  fft_in1(4 * n + iw) = cone * kernel * inv_dos
                  fft_in1(5 * n + iw) = cone * (kernel - muc * spin_fac) * inv_dos
                ENDDO
                !
                CALL cft_1z(fft_in1, narray(1), n, n, -1, fft_out1)
                ! divided by the factor n in cft_1z
                !
                fft_in2(0 * n + 1 : 1 * n) = &
                fft_out1(0 * n + 1 : 1 * n) * fft_out1(4 * n + 1 : 5 * n) * REAL(n, KIND=DP)
                fft_in2(1 * n + 1 : 2 * n) = &
                fft_out1(1 * n + 1 : 2 * n) * fft_out1(4 * n + 1 : 5 * n) * REAL(n, KIND=DP)
                fft_in2(2 * n + 1 : 3 * n) = &
                fft_out1(2 * n + 1 : 3 * n) * fft_out1(5 * n + 1 : 6 * n) * REAL(n, KIND=DP)
                fft_in2(3 * n + 1 : 4 * n) = &
                fft_out1(3 * n + 1 : 4 * n) * fft_out1(4 * n + 1 : 5 * n) * REAL(n, KIND=DP)
                !
                CALL cft_1z(fft_in2, narray(2), n, n, +1, fft_out2)
                !
                IF (positive_matsu) THEN
                  DO iw = 1, n/2
                    !
                    iw_f = iw
                    !
                    weight = wqf(iq)
                    aznormi(iw_f, ibnd, ik) = &
                    aznormi(iw_f, ibnd, ik) + weight * REAL(fft_out2(0 * n + iw), KIND=DP)
                    ashifti(iw_f, ibnd, ik) = &
                    ashifti(iw_f, ibnd, ik) + weight * REAL(fft_out2(1 * n + iw), KIND=DP)
                    adeltai(iw_f, ibnd, ik) = &
                    adeltai(iw_f, ibnd, ik) + weight * REAL(fft_out2(2 * n + iw), KIND=DP)
                    naznormi(iw_f, ibnd, ik) = &
                    naznormi(iw_f, ibnd, ik) + weight * REAL(fft_out2(3 * n + iw), KIND=DP)
                  ENDDO
                ELSE
                  DO iw = 1, n
                    !
                    IF (iw .LE. n/2) THEN
                      iw_f = iw + n/2 ! iw + N
                    ELSEIF (iw .GT. n/2) THEN
                      iw_f = iw - n/2 ! iw - N
                    ENDIF
                    !
                    weight = wqf(iq)
                    aznormi(iw_f, ibnd, ik) = &
                    aznormi(iw_f, ibnd, ik) + weight * REAL(fft_out2(0 * n + iw), KIND=DP)
                    ashifti(iw_f, ibnd, ik) = &
                    ashifti(iw_f, ibnd, ik) + weight * REAL(fft_out2(1 * n + iw), KIND=DP)
                    adeltai(iw_f, ibnd, ik) = &
                    adeltai(iw_f, ibnd, ik) + weight * REAL(fft_out2(2 * n + iw), KIND=DP)
                    naznormi(iw_f, ibnd, ik) = &
                    naznormi(iw_f, ibnd, ik) + weight * REAL(fft_out2(3 * n + iw), KIND=DP)
                  ENDDO
                ENDIF
                !
              ENDIF
            ENDDO
          ENDDO
        ENDIF
      ENDDO
    ENDDO
    !
    DEALLOCATE(inv_wsi, STAT = ierr)
    IF (ierr /= 0) CALL errore('sum_eliashberg_aniso_iaxis_fbw_fft', 'Error deallocating inv_wsi', 1)
    ! remove memory allocated for inv_wsi
    imelt = nsiw(itemp)
    CALL mem_size_eliashberg(2, -imelt)
    !
    RETURN
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE sum_eliashberg_aniso_iaxis_fbw_fft
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    SUBROUTINE sum_eliashberg_aniso_iaxis_fbw_ir(itemp, iter, ns, nel, nstate, ir_obj)
    !-----------------------------------------------------------------------
    !!
    !! This routine solves the anisotropic FBW Eliashberg equations on the imaginary-axis
    !! using sparse-ir sampling.
    !! This is used for uniform sampling and sparse sampling (gridsamp == 2).
    !!
    !! When considering only positive Matsubara frequencies,
    !! the data size of the Matsubara Green's function G(iw) is halved.
    !! Consequently, the imaginary-time Green's function G(tau)
    !! and the IR coefficients G_l can be defined as real-valued functions.
    !! The same applies to the kernel.
    !!
    USE kinds,             ONLY : DP
    USE control_flags,     ONLY : iverbosity
    USE global_var,        ONLY : wqf, gtemp
    USE input,             ONLY : muc, fsthick, muchem, eps_cut_ir, positive_matsu
    USE supercond_common,  ONLY : ixkqf, ixqfs, nqfs, ekfs, nkfs, nbndfs, spin_fac, &
                                  dosef, ef0, nsiw, wsi, adeltai, &
                                  adeltaip, aznormi, aznormip, naznormi, ashifti, &
                                  ashiftip, muintr, ir_giw, ir_gl, ir_gtau, ir_knliw, &
                                  ir_knll, ir_knltau, ir_cvliw, ir_cvll, ir_cvltau, &
                                  ir_gl_d, ir_gtau_d, ir_knll_d, ir_knltau_d, &
                                  ir_cvll_d, ir_cvltau_d, weight_q, num_js1, &
                                  gl_abs, fl_abs, knll_abs, siz_ir
    USE parallelism,       ONLY : fkbounds
    USE low_lvl,           ONLY : mem_size_eliashberg
    USE ep_constants,      ONLY : pi, zero, one, eps8, czero, cone, ci, eps16
    USE io_supercond,      ONLY : print_gl, print_kernell
    USE mp_global,         ONLY : inter_pool_comm
    USE mp,                ONLY : mp_sum
    USE sparse_ir,         ONLY : IR, fit_matsubara_b, fit_matsubara_f, fit_tau, &
                                  evaluate_matsubara_f, evaluate_tau
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: itemp
    !! Counter on temperature index
    INTEGER, INTENT(in) :: iter
    !! Counter on iteration steps
    INTEGER, INTENT(in) :: ns
    !! mu_inter parameters
    REAL(KIND = DP), INTENT(in) :: nel
    !! mu_inter parameters
    REAL(KIND = DP), INTENT(in) :: nstate
    !! mu_inter parameters
    TYPE(IR), INTENT(IN) :: ir_obj
    !! which contains ir objects such as basis functions
    !
    ! Local variables
    INTEGER :: iw
    !! Counter on frequency imag-axis
    INTEGER :: it
    !! Counter on imaginary time
    INTEGER :: il
    !! Counter on index of IR basis function
    INTEGER :: ik
    !! Counter on k-points
    INTEGER :: iq
    !! Counter on q-points for which k+sign*q is within the Fermi shell
    INTEGER :: iq0
    !! Index of iq on full q-mesh
    INTEGER :: ibnd
    !! Counter on bands at k
    INTEGER :: jbnd
    !! Counter on bands at k+q
    INTEGER :: lower_bnd, upper_bnd
    !! Lower/upper bound index after k paral
    INTEGER :: jstate
    !! Counter on states
    INTEGER :: jx
    !! Counter on states
    INTEGER :: jx2
    !! Counter on states
    INTEGER :: offset
    !! offset used when gridsamp == 2
    INTEGER(8) :: imelt
    !! Counter memory
    INTEGER :: ierr
    !! Error status
    !
    REAL(KIND = DP) :: kernel
    !! K_{-}(n,n',T))
    REAL(KIND = DP) :: esqrt
    !! Temporary variables
    REAL(KIND = DP) :: esqrt0
    !! esqrt and zesqrt calculated with adeltaip = 0.0, which is used for naznormi
    REAL(KIND = DP) :: inv_dos
    !! Invese dos inv_dos = 1/dosef. Defined for efficiency reason
    REAL(KIND = DP) :: omega
    !! Matsubara frequency
    REAL(KIND = DP) :: gl_ratio
    !! dummy for G(lmax)/G(l=1)
    REAL(KIND = DP) :: knll_ratio
    !! dummy for G(lmax)/G(l=1)
    REAL(KIND = DP) :: knll_ratio_max
    !! Maximum of  G(lmax)/G(l=1)
    REAL(KIND = DP) :: gl_max
    !! dummy for MAXVAL(ABS(ir_gl(jx2, :)))
    REAL(KIND = DP) :: knll_max
    !! dummy for MAXVAL(ABS(ir_knll(jx2, :)))
    REAL(KIND = DP), ALLOCATABLE :: inv_wsi(:)
    !! Invese imaginary freq. inv_wsi = 1/wsi. Defined for efficiency reason
    !
    inv_dos = one / dosef
    !
    CALL fkbounds(nkfs, lower_bnd, upper_bnd)
    !
    ! get the size of required memory for inv_wsi
    imelt = nsiw(itemp)
    CALL mem_size_eliashberg(2, imelt)
    ALLOCATE(inv_wsi(nsiw(itemp)), STAT = ierr)
    IF (ierr /= 0) CALL errore('sum_eliashberg_aniso_iaxis_fbw_ir', 'Error allocating inv_wsi', 1)
    !
    DO iw = 1, nsiw(itemp)
      inv_wsi(iw) = one / wsi(iw)
    ENDDO
    !
    ! HM: update the chemical potential from the previous one
    IF (muchem) CALL mu_inter_aniso(itemp, muintr, nel, nstate, ns, ir_obj)
    !
    naznormi(:, :, :) = zero
    adeltai(:, :, :)  = zero
    aznormi(:, :, :)  = zero
    ashifti(:, :, :)  = zero
    !
    ir_giw(:, :) = czero
    ir_knliw(:, :) = czero
    ir_cvliw(:, :) = czero
    weight_q(:) = zero
    IF (positive_matsu) THEN
      ir_gl_d(:, :) = zero
      ir_gtau_d(:, :) = zero
      ir_knll_d(:, :) = zero
      ir_knltau_d(:, :) = zero
      ir_cvll_d(:, :) = zero
      ir_cvltau_d(:, :) = zero
    ELSE
      ir_gl(:, :) = czero
      ir_gtau(:, :) = czero
      ir_knll(:, :) = czero
      ir_knltau(:, :) = czero
      ir_cvll(:, :) = czero
      ir_cvltau(:, :) = czero
    ENDIF
    !
    IF ((iverbosity == 4) .AND. (iter <= 5 .OR. iter == 15 .OR. MOD(iter, 10) == 0)) THEN
      gl_abs(:, :, :) = zero
      fl_abs(:, :, :) = zero
      IF (iter == 1) knll_abs(:, :, :) = zero
    ENDIF
    !
    DO ik = lower_bnd, upper_bnd
      DO ibnd = 1, nbndfs
        IF (ABS(ekfs(ibnd, ik) - ef0) < fsthick) THEN
          ! The below IF statement is just for outputting IR coefficients
          IF ((iverbosity == 4) .AND. (iter <= 5 .OR. iter == 15 .OR. MOD(iter, 10) == 0)) THEN
            IF (iter == 1) knll_ratio_max = zero
            DO iw = 1, ir_obj%nfreq_f
              ! ir_obj%nfreq_f == nsiw(itemp)
              !
              ! FOR Fermionic functions
              !
              ! wsi(iw) == DBLE(ir_obj%freq_f(iw)) * pi * gtemp(itemp)
              !
              esqrt = - one / ((wsi(iw) * aznormip(iw, ibnd, ik))**2.d0 &
                        + (ekfs(ibnd, ik) - muintr + ashiftip(iw, ibnd, ik))**2.d0 &
                        + (adeltaip(iw, ibnd, ik) * aznormip(iw, ibnd, ik))**2.d0)
              !
              ! imaginary part of normal Green's function
              ir_giw(1, iw) = ci * esqrt * wsi(iw) * aznormip(iw, ibnd, ik)
              !
              ! real part of normal Green's function
              ir_giw(1, iw) = ir_giw(1, iw) + &
                              cone * esqrt * (ekfs(ibnd, ik) - muintr + ashiftip(iw, ibnd, ik))
              !
              ! anomalous Green's function
              ir_giw(2, iw) = cone * esqrt * adeltaip(iw, ibnd, ik) * aznormip(iw, ibnd, ik)
              !
            ENDDO
            IF (positive_matsu) THEN
              CALL fit_matsubara_f (ir_obj, ir_giw, ir_gl_d)
              !
              gl_abs(1:ir_obj%size, ibnd, ik) = ABS(ir_gl_d(1, 1:ir_obj%size))
              fl_abs(1:ir_obj%size, ibnd, ik) = ABS(ir_gl_d(2, 1:ir_obj%size))
            ELSE
              CALL fit_matsubara_f (ir_obj, ir_giw, ir_gl)
              !
              gl_abs(1:ir_obj%size, ibnd, ik) = ABS(ir_gl(1, 1:ir_obj%size))
              fl_abs(1:ir_obj%size, ibnd, ik) = ABS(ir_gl(2, 1:ir_obj%size))
            ENDIF
          ENDIF
          jstate = 0
          ! HM: jstate is the index combining iq and jbnd.
          !     Each time jstate is divisible by siz_ir,
          !     calculate the convolutions and add the partial
          !     summation for each array.
          DO iq = 1, nqfs(ik)
            ! iq0 - index of q-point on the full q-mesh
            iq0 = ixqfs(ik, iq)
            DO jbnd = 1, nbndfs
              IF (ABS(ekfs(jbnd, ixkqf(ik, iq0)) - ef0) < fsthick) THEN
                !
                jstate = jstate + 1
                jx = MOD(jstate - 1, siz_ir) + 1
                !
                DO iw = 1, ir_obj%nfreq_f
                  ! ir_obj%nfreq_f == nsiw(itemp)
                  !
                  ! FOR Fermionic functions
                  !
                  ! wsi(iw) == DBLE(ir_obj%freq_f(iw)) * pi * gtemp(itemp)
                  !
                  esqrt = - one / ((wsi(iw) * aznormip(iw, jbnd, ixkqf(ik, iq0)))**2.d0 &
                            + (ekfs(jbnd, ixkqf(ik, iq0)) - muintr + ashiftip(iw, jbnd, ixkqf(ik, iq0)))**2.d0 &
                            + (adeltaip(iw, jbnd, ixkqf(ik, iq0)) * aznormip(iw, jbnd, ixkqf(ik, iq0)))**2.d0)
                  !
                  ! imaginary part of normal Green's function
                  offset = 0
                  ir_giw(offset + jx, iw) = ci * esqrt * wsi(iw) * aznormip(iw, jbnd, ixkqf(ik, iq0))
                  !
                  ! real part of normal Green's function
                  ir_giw(offset + jx, iw) = ir_giw(offset + jx, iw) + &
                  cone * esqrt * (ekfs(jbnd, ixkqf(ik, iq0)) - muintr + ashiftip(iw, jbnd, ixkqf(ik, iq0)))
                  !
                  ! anomalous Green's function
                  offset = siz_ir
                  ir_giw(offset + jx, iw) = &
                  cone * esqrt * adeltaip(iw, jbnd, ixkqf(ik, iq0)) * aznormip(iw, jbnd, ixkqf(ik, iq0))
                  !
                  esqrt0 = - one / ((wsi(iw) * aznormip(iw, jbnd, ixkqf(ik, iq0)))**2.d0 &
                            + (ekfs(jbnd, ixkqf(ik, iq0)) - muintr + ashiftip(iw, jbnd, ixkqf(ik, iq0)))**2.d0)
                  !
                  ! imaginary part of normal Green's function for normal state
                  offset = 2 * siz_ir
                  ir_giw(offset + jx, iw) = ci * esqrt0 * wsi(iw) * aznormip(iw, jbnd, ixkqf(ik, iq0))
                  !
                ENDDO
                !
                DO iw = 1, ir_obj%nfreq_b
                  !
                  ! FOR Bosonic functions
                  !
                  omega = DBLE(ir_obj%freq_b(iw)) * pi * gtemp(itemp)
                  !
                  CALL lambdar_aniso_ver1(ik, iq, ibnd, jbnd, omega, kernel)
                  !
                  ir_knliw(jx, iw) = cone * kernel * inv_dos
                  !
                ENDDO
                !
                weight_q(jx) = wqf(iq)
                !
                IF ((jx == siz_ir) .OR. (jstate == num_js1(ibnd, ik))) THEN
                  ! To obtain expansion coefficients of Green's function
                  ! and kernel in the IR basis
                  IF (positive_matsu) THEN
                    CALL fit_matsubara_f (ir_obj, ir_giw, ir_gl_d)
                    CALL fit_matsubara_b (ir_obj, ir_knliw, ir_knll_d)
                  ELSE
                    CALL fit_matsubara_f (ir_obj, ir_giw, ir_gl)
                    CALL fit_matsubara_b (ir_obj, ir_knliw, ir_knll)
                  ENDIF
                  !
                  ! The below IF statement is just for outputting IR coefficients
                  IF (positive_matsu) THEN
                    IF ((iverbosity == 4) .AND. (iter == 1)) THEN
                      DO jx2 = 1, jx
                        IF (ABS(ir_knll_d(jx2, 1)) > eps16) THEN
                          knll_ratio = ABS(ir_knll_d(jx2, ir_obj%size)) / MAXVAL(ABS(ir_knll_d(jx2, :)))
                          IF (knll_ratio > knll_ratio_max) THEN
                            knll_abs(1:ir_obj%size, ibnd, ik) = ABS(ir_knll_d(jx2, 1:ir_obj%size))
                            knll_ratio_max = knll_ratio
                          ENDIF
                        ENDIF
                      ENDDO
                    ENDIF
                  ELSE
                    IF ((iverbosity == 4) .AND. (iter == 1)) THEN
                      DO jx2 = 1, jx
                        IF (ABS(ir_knll(jx2, 1)) > eps16) THEN
                          knll_ratio = ABS(ir_knll(jx2, ir_obj%size)) / MAXVAL(ABS(ir_knll(jx2, :)))
                          IF (knll_ratio > knll_ratio_max) THEN
                            knll_abs(1:ir_obj%size, ibnd, ik) = ABS(ir_knll(jx2, 1:ir_obj%size))
                            knll_ratio_max = knll_ratio
                          ENDIF
                        ENDIF
                      ENDDO
                    ENDIF
                  ENDIF
                  !
                  IF (eps_cut_ir > zero) THEN
                    IF (positive_matsu) THEN
                      DO jx2 = 1, jx
                        gl_max = MAXVAL(ABS(ir_gl_d(jx2, :)))
                        knll_max = MAXVAL(ABS(ir_knll_d(jx2, :)))
                        DO il = 1, ir_obj%size ! loop over the index of IR basis functions
                          gl_ratio = ABS(ir_gl_d(jx2, il)) / gl_max
                          knll_ratio = ABS(ir_knll_d(jx2, il)) / knll_max
                          IF (gl_ratio < eps_cut_ir) THEN
                            ir_gl_d(jx2, il) = czero
                          ENDIF
                          IF (knll_ratio < eps_cut_ir) THEN
                            ir_knll_d(jx2, il) = czero
                          ENDIF
                        ENDDO
                      ENDDO
                    ELSE
                      DO jx2 = 1, jx
                        gl_max = MAXVAL(ABS(ir_gl(jx2, :)))
                        knll_max = MAXVAL(ABS(ir_knll(jx2, :)))
                        DO il = 1, ir_obj%size ! loop over the index of IR basis functions
                          gl_ratio = ABS(ir_gl(jx2, il)) / gl_max
                          knll_ratio = ABS(ir_knll(jx2, il)) / knll_max
                          IF (gl_ratio < eps_cut_ir) THEN
                            ir_gl(jx2, il) = czero
                          ENDIF
                          IF (knll_ratio < eps_cut_ir) THEN
                            ir_knll(jx2, il) = czero
                          ENDIF
                        ENDDO
                      ENDDO
                    ENDIF
                  ENDIF
                  !
                  IF (positive_matsu) THEN
                    ! To obtain Green's function and kernel of imaginary time
                    ! from the corresponding expansion coefficients
                    CALL evaluate_tau (ir_obj, ir_gl_d, ir_gtau_d)
                    CALL evaluate_tau (ir_obj, ir_knll_d, ir_knltau_d)
                    !
                    DO it = 1, ir_obj%ntau
                      offset = 0
                      ir_cvltau_d(offset + 1: offset + jx, it) = &
                      ir_gtau_d(offset + 1: offset + jx, it) * ir_knltau_d(1:jx, it) / gtemp(itemp)
                      offset = siz_ir
                      ir_cvltau_d(offset + 1: offset + jx, it) = &
                      ir_gtau_d(offset + 1: offset + jx, it) * ir_knltau_d(1:jx, it) / gtemp(itemp)
                      offset = 2 * siz_ir
                      ir_cvltau_d(offset + 1: offset + jx, it) = &
                      ir_gtau_d(offset + 1: offset + jx, it) * ir_knltau_d(1:jx, it) / gtemp(itemp)
                    ENDDO
                    !
                    ! To obtain expansion coefficients of the convolutions
                    ! in the IR basis
                    CALL fit_tau (ir_obj, ir_cvltau_d, ir_cvll_d)
                    !
                    ! To reconstruct the convolutions of Matsubara freq.
                    ! from the coefficients in the IR basis
                    CALL evaluate_matsubara_f (ir_obj, ir_cvll_d, ir_cvliw)
                    !
                  ELSE
                    ! To obtain green's function and kernel of imaginary time
                    ! from the corresponding expansion coefficients
                    CALL evaluate_tau (ir_obj, ir_gl, ir_gtau)
                    CALL evaluate_tau (ir_obj, ir_knll, ir_knltau)
                    !
                    DO it = 1, ir_obj%ntau
                      offset = 0
                      ir_cvltau(offset + 1: offset + jx, it) = &
                      ir_gtau(offset + 1: offset + jx, it) * ir_knltau(1:jx, it) / gtemp(itemp)
                      offset = siz_ir
                      ir_cvltau(offset + 1: offset + jx, it) = &
                      ir_gtau(offset + 1: offset + jx, it) * ir_knltau(1:jx, it) / gtemp(itemp)
                      offset = 2 * siz_ir
                      ir_cvltau(offset + 1: offset + jx, it) = &
                      ir_gtau(offset + 1: offset + jx, it) * ir_knltau(1:jx, it) / gtemp(itemp)
                    ENDDO
                    !
                    ! To obtain expansion coefficients of the convolutions
                    ! in the IR basis
                    CALL fit_tau (ir_obj, ir_cvltau, ir_cvll)
                    !
                    ! To reconstruct the convolutions of Matsubara freq.
                    ! from the coefficients in the IR basis
                    CALL evaluate_matsubara_f (ir_obj, ir_cvll, ir_cvliw)
                    !
                  ENDIF
                  !
                  DO iw = 1, ir_obj%nfreq_f
                    !
                    offset = 0
                    aznormi(iw, ibnd, ik) = &
                    aznormi(iw, ibnd, ik) - SUM(weight_q(1:jx) * AIMAG(ir_cvliw(offset + 1: offset + jx, iw)))
                    ashifti(iw, ibnd, ik) = &
                    ashifti(iw, ibnd, ik) - SUM(weight_q(1:jx) * REAL(ir_cvliw(offset + 1: offset + jx, iw), KIND=DP))
                    offset = siz_ir
                    adeltai(iw, ibnd, ik) = &
                    adeltai(iw, ibnd, ik) - SUM(weight_q(1:jx) * REAL(ir_cvliw(offset + 1: offset + jx, iw), KIND=DP))
                    offset = 2 * siz_ir
                    naznormi(iw, ibnd, ik) = &
                    naznormi(iw, ibnd, ik) - SUM(weight_q(1:jx) * AIMAG(ir_cvliw(offset + 1: offset + jx, iw)))
                    !
                    IF (ABS(muc) .GT. eps8) THEN
                      offset = siz_ir
                      ! HM: IR cannot reproduce a constant function on frequencies,
                      !     so we have to treat the term multiplied by muc separately.
                      !     Be careful so as not to mistake the sign!
                      IF (positive_matsu) THEN
                        !     ir_gtau_d(siz_ir + 1: siz_ir + jx, 1) is an anomalous Green's function of imaginary time at tau = 0+
                        adeltai(iw, ibnd, ik) = &
                        adeltai(iw, ibnd, ik) + &
                        muc * spin_fac * SUM(weight_q(1:jx) * REAL(ir_gtau_d(offset + 1: offset + jx, 1), KIND=DP)) &
                        * inv_dos / gtemp(itemp)
                      ELSE
                        !     ir_gtau(siz_ir + 1: siz_ir + jx, 1) is an anomalous Green's function of imaginary time at tau = 0+
                        adeltai(iw, ibnd, ik) = &
                        adeltai(iw, ibnd, ik) + &
                        muc * spin_fac * SUM(weight_q(1:jx) * REAL(ir_gtau(offset + 1: offset + jx, 1), KIND=DP)) &
                        * inv_dos / gtemp(itemp)
                      ENDIF
                    ENDIF
                    !
                  ENDDO
                  !
                  ir_giw(:, :) = czero
                  ir_knliw(:, :) = czero
                  ir_cvliw(:, :) = czero
                  IF (positive_matsu) THEN
                    ir_gl_d(:, :) = zero
                    ir_gtau_d(:, :) = zero
                    ir_knll_d(:, :) = zero
                    ir_knltau_d(:, :) = zero
                    ir_cvll_d(:, :) = zero
                    ir_cvltau_d(:, :) = zero
                  ELSE
                    ir_gl(:, :) = czero
                    ir_gtau(:, :) = czero
                    ir_knll(:, :) = czero
                    ir_knltau(:, :) = czero
                    ir_cvll(:, :) = czero
                    ir_cvltau(:, :) = czero
                  ENDIF
                ENDIF ! (jx == siz_ir) .OR. (jstate == num_js1(ibnd, ik))
                !
              ENDIF
            ENDDO
          ENDDO
        ENDIF
      ENDDO
    ENDDO
    IF ((iverbosity == 4) .AND. (iter == 1)) THEN
      ! collect contributions from all pools
      CALL mp_sum(knll_abs, inter_pool_comm)
      CALL print_kernell(itemp, ir_obj%size, knll_abs)
    ENDIF
    IF ((iverbosity == 4) .AND. (iter <= 5 .OR. iter == 15 .OR. MOD(iter, 10) == 0)) THEN
      ! collect contributions from all pools
      CALL mp_sum(gl_abs,   inter_pool_comm)
      CALL mp_sum(fl_abs,   inter_pool_comm)
      CALL print_gl(itemp, iter, ir_obj%size, gl_abs, fl_abs)
    ENDIF
    !
    DEALLOCATE(inv_wsi, STAT = ierr)
    IF (ierr /= 0) CALL errore('sum_eliashberg_aniso_iaxis_fbw_ir', 'Error deallocating inv_wsi', 1)
    ! remove memory allocated for inv_wsi
    imelt = nsiw(itemp)
    CALL mem_size_eliashberg(2, -imelt)
    !
    RETURN
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE sum_eliashberg_aniso_iaxis_fbw_ir
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    SUBROUTINE sum_eliash_aniso_iaxis_fbw_ir_coul(itemp, iter, ns, nel, nstate, ir_obj)
    !-----------------------------------------------------------------------
    !!
    !! This routine solves the anisotropic FBW Eliashberg equations,
    !! which include the Coulomb contribution arising from outside fsthick window,
    !! on the imaginary-axis using sparse-ir sampling (icoulomb > 0).
    !! This is used for uniform sampling and sparse sampling (gridsamp == 2).
    !!
    !! When considering only positive Matsubara frequencies,
    !! the data size of the Matsubara Green's function G(iw) is halved.
    !! Consequently, the imaginary-time Green's function G(tau)
    !! and the IR coefficients G_l can be defined as real-valued functions.
    !! The same applies to the kernel.
    !!
    USE kinds,             ONLY : DP
    USE control_flags,     ONLY : iverbosity
    USE global_var,        ONLY : wqf, gtemp
    USE input,             ONLY : muc, fsthick, muchem, eps_cut_ir, positive_matsu, &
                                  emax_coulomb, emin_coulomb
    USE supercond_common,  ONLY : ixkqf, ixqfs, nqfs, ekfs, nkfs, nbndfs, spin_fac, &
                                  dosef, ef0, nsiw, wsi, adeltai, &
                                  adeltaip, aznormi, aznormip, naznormi, ashifti, &
                                  ashiftip, muintr, ir_giw, ir_gl, ir_gtau, ir_knliw, &
                                  ir_knll, ir_knltau, ir_cvliw, ir_cvll, ir_cvltau, &
                                  ir_gl_d, ir_gtau_d, ir_knll_d, ir_knltau_d, &
                                  ir_cvll_d, ir_cvltau_d, weight_q, weight_cl, num_js1, &
                                  num_js2, num_js3, ir_giw_cl, ir_gl_cl, ir_gtau_cl, &
                                  ir_gl_cl_d, ir_gtau_cl_d, nk1_cl, nk2_cl, nk3_cl, &
                                  ik_bz_to_ibz_cl, nbnd_offset, ibnd_kfs_to_kfs_all, &
                                  nbnd_cl, nbndfs_all, ik_cl_to_fs, nkstot_cl, ek_cl, &
                                  ibnd_kfs_all_to_kfs, adeltai_cl, adeltaip_cl, w_stat, &
                                  gl_abs, fl_abs, knll_abs, siz_ir, siz_ir_cl
    USE parallelism,       ONLY : fkbounds, para_bounds
    USE low_lvl,           ONLY : mem_size_eliashberg
    USE ep_constants,      ONLY : pi, zero, one, eps8, czero, cone, ci, eps16
    USE io_supercond,      ONLY : print_gl, print_kernell
    USE mp_global,         ONLY : inter_pool_comm
    USE mp,                ONLY : mp_sum
    USE sparse_ir,         ONLY : IR, fit_matsubara_b, fit_matsubara_f, fit_tau, &
                                  evaluate_matsubara_f, evaluate_tau
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: itemp
    !! Counter on temperature index
    INTEGER, INTENT(in) :: iter
    !! Counter on iteration steps
    INTEGER, INTENT(in) :: ns
    !! mu_inter parameters
    REAL(KIND = DP), INTENT(in) :: nel
    !! mu_inter parameters
    REAL(KIND = DP), INTENT(in) :: nstate
    !! mu_inter parameters
    TYPE(IR), INTENT(IN) :: ir_obj
    !! which contains ir objects such as basis functions
    !
    ! Local variables
    LOGICAL :: linsidei
    !! used to determine if the state (ibnd_cl, ik_cl) is outside the window.
    LOGICAL :: linsidei2
    !! used to determine if the state (ibnd_cl, ik_cl) is outside the window.
    LOGICAL :: linsidei3
    !! used to determine if the state (ibnd_cl, ik_cl) is outside the window.
    LOGICAL :: linsidej
    !! used to determine if the state (jbnd_cl, jk_cl) is outside the window.
    LOGICAL :: linsidej2
    !! used to determine if the state (jbnd_cl, jk_cl) is outside the window.
    LOGICAL :: linsidej3
    !! used to determine if the state (jbnd_cl, jk_cl) is outside the window.
    INTEGER :: iw
    !! Counter on frequency imag-axis
    INTEGER :: it
    !! Counter on imaginary time
    INTEGER :: il
    !! Counter on index of IR basis function
    INTEGER :: ik
    !! Counter on k-points
    INTEGER :: iq
    !! Counter on q-points for which k+sign*q is within the Fermi shell
    INTEGER :: iq0
    !! Index of iq on full q-mesh
    INTEGER :: ibnd
    !! Counter on bands at k
    INTEGER :: jbnd
    !! Counter on bands at k+q
    INTEGER :: ik_cl
    !! Counter on k-points
    INTEGER :: jk
    !! Counter on k-points
    INTEGER :: jk_cl
    !! Counter on k-points: jk_cl = ik_bz_to_ibz_cl(jk)
    INTEGER :: ibnd_cl
    !! Counter on bands at k
    INTEGER :: jbnd_cl
    !! Counter on bands at k'
    INTEGER :: ibndfs
    !! Counter on bands in the fsthick window
    INTEGER :: jbndfs
    !! Counter on bands in the fsthick window
    INTEGER :: lower_bnd, upper_bnd
    !! Lower/upper bound index after k paral
    INTEGER :: jstate
    !! Counter on states
    INTEGER :: jx
    !! Counter on states
    INTEGER :: jx2
    !! Counter on states
    INTEGER :: istate
    !! Counter on states, used when icoulomb > 0
    INTEGER :: is_start, is_stop
    !! Lower/upper bound index after paral for states
    INTEGER :: num_states
    !! Number of states per pool: num_states = is_stop - is_start + 1
    INTEGER :: offset
    !! offset used when gridsamp == 2
    INTEGER(8) :: imelt
    !! Counter memory
    INTEGER :: ierr
    !! Error status
    !
    REAL(KIND = DP) :: kernel
    !! K_{-}(n,n',T))
    REAL(KIND = DP) :: esqrt
    !! Temporary variables
    REAL(KIND = DP) :: esqrt0
    !! esqrt and zesqrt calculated with adeltaip = 0.0, which is used for naznormi
    REAL(KIND = DP) :: inv_dos
    !! Invese dos inv_dos = 1/dosef. Defined for efficiency reason
    REAL(KIND = DP) :: omega
    !! Matsubara frequency
    REAL(KIND = DP) :: gl_ratio
    !! dummy for G(lmax)/G(l=1)
    REAL(KIND = DP) :: knll_ratio
    !! dummy for G(lmax)/G(l=1)
    REAL(KIND = DP) :: knll_ratio_max
    !! Maximum of  G(lmax)/G(l=1)
    REAL(KIND = DP) :: gl_max
    !! dummy for MAXVAL(ABS(ir_gl(jx2, :)))
    REAL(KIND = DP) :: knll_max
    !! dummy for MAXVAL(ABS(ir_knll(jx2, :)))
    REAL(KIND = DP), ALLOCATABLE :: inv_wsi(:)
    !! Invese imaginary freq. inv_wsi = 1/wsi. Defined for efficiency reason
    !
    inv_dos = one / dosef
    !
    CALL fkbounds(nkfs, lower_bnd, upper_bnd)
    !
    ! HM: Because nkstot_cl can be quite small,
    !     we do a combined parallelization on band index and k point index.
    CALL para_bounds(is_start, is_stop, nbnd_cl * nkstot_cl)
    num_states = is_stop - is_start + 1
    !
    ! get the size of required memory for inv_wsi
    imelt = nsiw(itemp)
    CALL mem_size_eliashberg(2, imelt)
    ALLOCATE(inv_wsi(nsiw(itemp)), STAT = ierr)
    IF (ierr /= 0) CALL errore('sum_eliash_aniso_iaxis_fbw_ir_coul', 'Error allocating inv_wsi', 1)
    !
    DO iw = 1, nsiw(itemp)
      inv_wsi(iw) = one / wsi(iw)
    ENDDO
    !
    ! coulomb interaction given by the following form will be implemented in the future
    ! w = w_stat(jbnd_cl, jk_cl, ibnd_cl, ik_cl) + w_dyn(iw_b, jbnd_cl, jk_cl, ibnd_cl, ik_cl)
    w_stat(:,:) = muc * spin_fac * inv_dos
    !
    ! HM: update the chemical potential from the previous one
    IF (muchem) CALL mu_inter_aniso(itemp, muintr, nel, nstate, ns, ir_obj)
    !
    naznormi(:, :, :) = zero
    adeltai(:, :, :)  = zero
    aznormi(:, :, :)  = zero
    ashifti(:, :, :)  = zero
    !
    ir_giw(:, :) = czero
    ir_knliw(:, :) = czero
    ir_cvliw(:, :) = czero
    IF (positive_matsu) THEN
      ir_gl_d(:, :) = zero
      ir_gtau_d(:, :) = zero
      ir_knll_d(:, :) = zero
      ir_knltau_d(:, :) = zero
      ir_cvll_d(:, :) = zero
      ir_cvltau_d(:, :) = zero
      !
      ir_giw_cl(:, :) = czero
      ir_gl_cl_d(:, :) = zero
      ir_gtau_cl_d(:, :) = zero
    ELSE
      ir_gl(:, :) = czero
      ir_gtau(:, :) = czero
      ir_knll(:, :) = czero
      ir_knltau(:, :) = czero
      ir_cvll(:, :) = czero
      ir_cvltau(:, :) = czero
      !
      ir_giw_cl(:, :) = czero
      ir_gl_cl(:, :) = czero
      ir_gtau_cl(:, :) = czero
    ENDIF
    !
    IF ((iverbosity == 4) .AND. (iter <= 5 .OR. iter == 15 .OR. MOD(iter, 10) == 0)) THEN
      gl_abs(:, :, :) = zero
      fl_abs(:, :, :) = zero
      IF (iter == 1) knll_abs(:, :, :) = zero
    ENDIF
    !
    DO ik = lower_bnd, upper_bnd
      DO ibnd = 1, nbndfs
        IF (ABS(ekfs(ibnd, ik) - ef0) < fsthick) THEN
          ! The below IF statement is just for outputting IR coefficients
          IF ((iverbosity == 4) .AND. (iter <= 5 .OR. iter == 15 .OR. MOD(iter, 10) == 0)) THEN
            IF (iter == 1) knll_ratio_max = zero
            DO iw = 1, ir_obj%nfreq_f
              ! ir_obj%nfreq_f == nsiw(itemp)
              !
              ! FOR Fermionic functions
              !
              ! wsi(iw) == DBLE(ir_obj%freq_f(iw)) * pi * gtemp(itemp)
              !
              esqrt = - one / ((wsi(iw) * aznormip(iw, ibnd, ik))**2.d0 &
                        + (ekfs(ibnd, ik) - muintr + ashiftip(iw, ibnd, ik))**2.d0 &
                        + (adeltaip(iw, ibnd, ik) * aznormip(iw, ibnd, ik))**2.d0)
              !
              ! imaginary part of normal Green's function
              ir_giw(1, iw) = ci * esqrt * wsi(iw) * aznormip(iw, ibnd, ik)
              !
              ! real part of normal Green's function
              ir_giw(1, iw) = ir_giw(1, iw) + &
                              cone * esqrt * (ekfs(ibnd, ik) - muintr + ashiftip(iw, ibnd, ik))
              !
              ! anomalous Green's function
              ir_giw(2, iw) = cone * esqrt * adeltaip(iw, ibnd, ik) * aznormip(iw, ibnd, ik)
              !
            ENDDO
            IF (positive_matsu) THEN
              CALL fit_matsubara_f (ir_obj, ir_giw, ir_gl_d)
              !
              gl_abs(1:ir_obj%size, ibnd, ik) = ABS(ir_gl_d(1, 1:ir_obj%size))
              fl_abs(1:ir_obj%size, ibnd, ik) = ABS(ir_gl_d(2, 1:ir_obj%size))
            ELSE
              CALL fit_matsubara_f (ir_obj, ir_giw, ir_gl)
              !
              gl_abs(1:ir_obj%size, ibnd, ik) = ABS(ir_gl(1, 1:ir_obj%size))
              fl_abs(1:ir_obj%size, ibnd, ik) = ABS(ir_gl(2, 1:ir_obj%size))
            ENDIF
          ENDIF
          !
          ibnd_cl = ibnd_kfs_to_kfs_all(ibnd, ik) + nbnd_offset
          !nn = 0 ! DEBUG
          jstate = 0
          ! HM: jstate is the index combining iq and jbnd.
          !     Each time jstate is divisible by siz_ir,
          !     calculate the convolutions and add the partial
          !     summation for each array.
          DO iq = 1, nqfs(ik)
            ! iq0 - index of q-point on the full q-mesh
            iq0 = ixqfs(ik, iq)
            DO jbnd = 1, nbndfs
              IF (ABS(ekfs(jbnd, ixkqf(ik, iq0)) - ef0) < fsthick) THEN
                !
                jstate = jstate + 1
                jx = MOD(jstate - 1, siz_ir) + 1
                !nn = nn + 1 ! DEBUG
                DO iw = 1, ir_obj%nfreq_f
                  ! ir_obj%nfreq_f == nsiw(itemp)
                  !
                  ! FOR Fermionic functions
                  !
                  ! wsi(iw) == DBLE(ir_obj%freq_f(iw)) * pi * gtemp(itemp)
                  !
                  esqrt = - one / ((wsi(iw) * aznormip(iw, jbnd, ixkqf(ik, iq0)))**2.d0 &
                            + (ekfs(jbnd, ixkqf(ik, iq0)) - muintr + ashiftip(iw, jbnd, ixkqf(ik, iq0)))**2.d0 &
                            + (adeltaip(iw, jbnd, ixkqf(ik, iq0)) * aznormip(iw, jbnd, ixkqf(ik, iq0)))**2.d0)
                  !
                  ! imaginary part of normal Green's function
                  offset = 0
                  ir_giw(offset + jx, iw) = ci * esqrt * wsi(iw) * aznormip(iw, jbnd, ixkqf(ik, iq0))
                  !
                  ! real part of normal Green's function
                  ir_giw(offset + jx, iw) = ir_giw(offset + jx, iw) + &
                  cone * esqrt * (ekfs(jbnd, ixkqf(ik, iq0)) - muintr + ashiftip(iw, jbnd, ixkqf(ik, iq0)))
                  !
                  ! anomalous Green's function
                  offset = siz_ir
                  ir_giw(offset + jx, iw) = &
                  cone * esqrt * adeltaip(iw, jbnd, ixkqf(ik, iq0)) * aznormip(iw, jbnd, ixkqf(ik, iq0))
                  !
                  esqrt0 = - one / ((wsi(iw) * aznormip(iw, jbnd, ixkqf(ik, iq0)))**2.d0 &
                            + (ekfs(jbnd, ixkqf(ik, iq0)) - muintr + ashiftip(iw, jbnd, ixkqf(ik, iq0)))**2.d0)
                  !
                  ! imaginary part of normal Green's function for normal state
                  offset = 2 * siz_ir
                  ir_giw(offset + jx, iw) = ci * esqrt0 * wsi(iw) * aznormip(iw, jbnd, ixkqf(ik, iq0))
                  !
                ENDDO
                !
                DO iw = 1, ir_obj%nfreq_b
                  !
                  ! FOR Bosonic functions
                  !
                  omega = DBLE(ir_obj%freq_b(iw)) * pi * gtemp(itemp)
                  !
                  CALL lambdar_aniso_ver1(ik, iq, ibnd, jbnd, omega, kernel)
                  !
                  ir_knliw(jx, iw) = cone * kernel * inv_dos
                  !
                ENDDO
                !
                weight_q(jx) = wqf(iq)
                !
                IF ((jx == siz_ir) .OR. (jstate == num_js1(ibnd, ik))) THEN
                  ! To obtain expansion coefficients of green's function
                  ! and kernel in the IR basis
                  IF (positive_matsu) THEN
                    CALL fit_matsubara_f (ir_obj, ir_giw, ir_gl_d)
                    CALL fit_matsubara_b (ir_obj, ir_knliw, ir_knll_d)
                  ELSE
                    CALL fit_matsubara_f (ir_obj, ir_giw, ir_gl)
                    CALL fit_matsubara_b (ir_obj, ir_knliw, ir_knll)
                  ENDIF
                  !
                  ! The below IF statement is just for outputting IR coefficients
                  IF (positive_matsu) THEN
                    IF ((iverbosity == 4) .AND. (iter == 1)) THEN
                      DO jx2 = 1, jx
                        IF (ABS(ir_knll_d(jx2, 1)) > eps16) THEN
                          knll_ratio = ABS(ir_knll_d(jx2, ir_obj%size)) / MAXVAL(ABS(ir_knll_d(jx2, :)))
                          IF (knll_ratio > knll_ratio_max) THEN
                            knll_abs(1:ir_obj%size, ibnd, ik) = ABS(ir_knll_d(jx2, 1:ir_obj%size))
                            knll_ratio_max = knll_ratio
                          ENDIF
                        ENDIF
                      ENDDO
                    ENDIF
                  ELSE
                    IF ((iverbosity == 4) .AND. (iter == 1)) THEN
                      DO jx2 = 1, jx
                        IF (ABS(ir_knll(jx2, 1)) > eps16) THEN
                          knll_ratio = ABS(ir_knll(jx2, ir_obj%size)) / MAXVAL(ABS(ir_knll(jx2, :)))
                          IF (knll_ratio > knll_ratio_max) THEN
                            knll_abs(1:ir_obj%size, ibnd, ik) = ABS(ir_knll(jx2, 1:ir_obj%size))
                            knll_ratio_max = knll_ratio
                          ENDIF
                        ENDIF
                      ENDDO
                    ENDIF
                  ENDIF
                  !
                  IF (eps_cut_ir > zero) THEN
                    IF (positive_matsu) THEN
                      DO jx2 = 1, jx
                        gl_max = MAXVAL(ABS(ir_gl_d(jx2, :)))
                        knll_max = MAXVAL(ABS(ir_knll_d(jx2, :)))
                        DO il = 1, ir_obj%size ! loop over the index of IR basis functions
                          gl_ratio = ABS(ir_gl_d(jx2, il)) / gl_max
                          knll_ratio = ABS(ir_knll_d(jx2, il)) / knll_max
                          IF (gl_ratio < eps_cut_ir) THEN
                            ir_gl_d(jx2, il) = czero
                          ENDIF
                          IF (knll_ratio < eps_cut_ir) THEN
                            ir_knll_d(jx2, il) = czero
                          ENDIF
                        ENDDO
                      ENDDO
                    ELSE
                      DO jx2 = 1, jx
                        gl_max = MAXVAL(ABS(ir_gl(jx2, :)))
                        knll_max = MAXVAL(ABS(ir_knll(jx2, :)))
                        DO il = 1, ir_obj%size ! loop over the index of IR basis functions
                          gl_ratio = ABS(ir_gl(jx2, il)) / gl_max
                          knll_ratio = ABS(ir_knll(jx2, il)) / knll_max
                          IF (gl_ratio < eps_cut_ir) THEN
                            ir_gl(jx2, il) = czero
                          ENDIF
                          IF (knll_ratio < eps_cut_ir) THEN
                            ir_knll(jx2, il) = czero
                          ENDIF
                        ENDDO
                      ENDDO
                    ENDIF
                  ENDIF
                  !
                  IF (positive_matsu) THEN
                    ! To obtain Green's function and kernel of imaginary time
                    ! from the corresponding expansion coefficients
                    CALL evaluate_tau (ir_obj, ir_gl_d, ir_gtau_d)
                    CALL evaluate_tau (ir_obj, ir_knll_d, ir_knltau_d)
                    !
                    DO it = 1, ir_obj%ntau
                      offset = 0
                      ir_cvltau_d(offset + 1: offset + jx, it) = &
                      ir_gtau_d(offset + 1: offset + jx, it) * ir_knltau_d(1:jx, it) / gtemp(itemp)
                      offset = siz_ir
                      ir_cvltau_d(offset + 1: offset + jx, it) = &
                      ir_gtau_d(offset + 1: offset + jx, it) * ir_knltau_d(1:jx, it) / gtemp(itemp)
                      offset = 2 * siz_ir
                      ir_cvltau_d(offset + 1: offset + jx, it) = &
                      ir_gtau_d(offset + 1: offset + jx, it) * ir_knltau_d(1:jx, it) / gtemp(itemp)
                    ENDDO
                    !
                    ! To obtain expansion coefficients of the convolutions
                    ! in the IR basis
                    CALL fit_tau (ir_obj, ir_cvltau_d, ir_cvll_d)
                    !
                    ! To reconstruct the convolutions of Matsubara freq.
                    ! from the coefficients in the IR basis
                    CALL evaluate_matsubara_f (ir_obj, ir_cvll_d, ir_cvliw)
                    !
                  ELSE
                    ! To obtain Green's function and kernel of imaginary time
                    ! from the corresponding expansion coefficients
                    CALL evaluate_tau (ir_obj, ir_gl, ir_gtau)
                    CALL evaluate_tau (ir_obj, ir_knll, ir_knltau)
                    !
                    DO it = 1, ir_obj%ntau
                      offset = 0
                      ir_cvltau(offset + 1: offset + jx, it) = &
                      ir_gtau(offset + 1: offset + jx, it) * ir_knltau(1:jx, it) / gtemp(itemp)
                      offset = siz_ir
                      ir_cvltau(offset + 1: offset + jx, it) = &
                      ir_gtau(offset + 1: offset + jx, it) * ir_knltau(1:jx, it) / gtemp(itemp)
                      offset = 2 * siz_ir
                      ir_cvltau(offset + 1: offset + jx, it) = &
                      ir_gtau(offset + 1: offset + jx, it) * ir_knltau(1:jx, it) / gtemp(itemp)
                    ENDDO
                    !
                    ! To obtain expansion coefficients of the convolutions
                    ! in the IR basis
                    CALL fit_tau (ir_obj, ir_cvltau, ir_cvll)
                    !
                    ! To reconstruct the convolutions of Matsubara freq.
                    ! from the coefficients in the IR basis
                    CALL evaluate_matsubara_f (ir_obj, ir_cvll, ir_cvliw)
                    !
                  ENDIF
                  !
                  DO iw = 1, ir_obj%nfreq_f
                    !
                    !IF ((ik == lower_bnd) .AND. (ibnd == 1) .AND. (nn == 1)) WRITE(stdout, '(5x, a, ES20.10)')    'weight = ', weight ! DEBUG
                    offset = 0
                    aznormi(iw, ibnd, ik) = &
                    aznormi(iw, ibnd, ik) - SUM(weight_q(1:jx) * AIMAG(ir_cvliw(offset + 1: offset + jx, iw)))
                    ashifti(iw, ibnd, ik) = &
                    ashifti(iw, ibnd, ik) - SUM(weight_q(1:jx) * REAL(ir_cvliw(offset + 1: offset + jx, iw), KIND=DP))
                    offset = siz_ir
                    adeltai(iw, ibnd, ik) = &
                    adeltai(iw, ibnd, ik) - SUM(weight_q(1:jx) * REAL(ir_cvliw(offset + 1: offset + jx, iw), KIND=DP))
                    offset = 2 * siz_ir
                    naznormi(iw, ibnd, ik) = &
                    naznormi(iw, ibnd, ik) - SUM(weight_q(1:jx) * AIMAG(ir_cvliw(offset + 1: offset + jx, iw)))
                    !
                  ENDDO
                  !
                  ir_giw(:, :) = czero
                  ir_knliw(:, :) = czero
                  ir_cvliw(:, :) = czero
                  IF (positive_matsu) THEN
                    ir_gl_d(:, :) = zero
                    ir_gtau_d(:, :) = zero
                    ir_knll_d(:, :) = zero
                    ir_knltau_d(:, :) = zero
                    ir_cvll_d(:, :) = zero
                    ir_cvltau_d(:, :) = zero
                  ELSE
                    ir_gl(:, :) = czero
                    ir_gtau(:, :) = czero
                    ir_knll(:, :) = czero
                    ir_knltau(:, :) = czero
                    ir_cvll(:, :) = czero
                    ir_cvltau(:, :) = czero
                  ENDIF
                ENDIF
                !
              ENDIF
            ENDDO
          ENDDO
          !IF ((ik == upper_bnd) .AND. (ibnd == 1)) WRITE(stdout, '(5x, a, I4)')    'nn = ', nn ! DEBUG
          !IF ((ik == upper_bnd) .AND. (ibnd == 1)) WRITE(1000+mpime, '(5x, a, I4)')    'nn = ', nn ! DEBUG
          !
          jstate = 0
          ! HM: jstate is the index combining jk_cl and jbnd_cl.
          !     Each time jstate is divisible by siz_ir_cl,
          !     calculate the convolutions and add the partial
          !     summation for each array.
          linsidej = .TRUE.
          linsidej2 = .TRUE.
          linsidej3 = .TRUE.
          !nn = 0 ! DEBUG
          DO jk = 1, (nk1_cl * nk2_cl * nk3_cl)
            ! find the index of the irreducible k point corresponding to the current k point
            jk_cl = ik_bz_to_ibz_cl(jk)
            !
            ! initilize linsidej
            linsidej = .TRUE.
            IF (ik_cl_to_fs(jk_cl) == 0) linsidej = .FALSE.
            !IF (.NOT.(linsidej)) CYCLE ! DEBUG
            DO jbnd_cl = 1, nbnd_cl
              ! initilize linsidej2
              linsidej2 = .TRUE.
              IF (((jbnd_cl - nbnd_offset) < 1) .OR. ((jbnd_cl - nbnd_offset) > nbndfs_all)) &
                linsidej2 = .FALSE.
              !
              ! initilize linsidej3
              linsidej3 = .TRUE.
              IF (linsidej .AND. linsidej2) THEN
                jbndfs = ibnd_kfs_all_to_kfs(jbnd_cl - nbnd_offset, ik_cl_to_fs(jk_cl))
                IF ((jbndfs < 1) .OR. (jbndfs > nbndfs)) THEN
                  linsidej3 = .FALSE.
                ELSE
                  ! Always use ekfs instead of ek_cl to determine
                  ! if the state (jbnd_cl, jk_cl) is outside the window.
                  IF (ABS(ekfs(jbndfs, ik_cl_to_fs(jk_cl)) - ef0) >= fsthick) &
                    linsidej3 = .FALSE.
                ENDIF
              ENDIF
              !
              !IF (.NOT.(linsidej .AND. linsidej2 .AND. linsidej3)) CYCLE ! DEBUG
              !nn = nn + 1 ! DEBUG
              IF (.NOT.(linsidej .AND. linsidej2 .AND. linsidej3)) THEN
                IF (((ek_cl(jbnd_cl, jk_cl) - ef0) >= emax_coulomb) .OR. &
                    ((ek_cl(jbnd_cl, jk_cl) - ef0) <= emin_coulomb)) THEN
                  CYCLE
                ENDIF
              ENDIF
              !
              jstate = jstate + 1
              jx = MOD(jstate - 1, siz_ir_cl) + 1
              DO iw = 1, ir_obj%nfreq_f
                ! ir_obj%nfreq_f == nsiw(itemp)
                !
                ! FOR Fermionic functions
                !
                ! wsi(iw) == DBLE(ir_obj%freq_f(iw)) * pi * gtemp(itemp)
                !
                IF (linsidej .AND. linsidej2 .AND. linsidej3) THEN
                  esqrt = - one / ((wsi(iw) * aznormip(iw, jbndfs, ik_cl_to_fs(jk_cl)))**2.d0 &
                  + (ekfs(jbndfs, ik_cl_to_fs(jk_cl)) - muintr + ashiftip(iw, jbndfs, ik_cl_to_fs(jk_cl)))**2.d0 &
                  + (adeltaip(iw, jbndfs, ik_cl_to_fs(jk_cl)) * aznormip(iw, jbndfs, ik_cl_to_fs(jk_cl)))**2.d0)
                  !
                  ! anomalous Green's function
                  ir_giw_cl(jx, iw) = cone * esqrt * adeltaip(iw, jbndfs, ik_cl_to_fs(jk_cl)) &
                                          * aznormip(iw, jbndfs, ik_cl_to_fs(jk_cl))
                ELSE
                  esqrt = - one / (wsi(iw)**2.d0 &
                  + (ek_cl(jbnd_cl, jk_cl) - muintr)**2.d0 &
                  + (adeltaip_cl(iw, jbnd_cl, jk_cl))**2.d0)
                  !
                  ! anomalous Green's function
                  ir_giw_cl(jx, iw) = cone * esqrt * adeltaip_cl(iw, jbnd_cl, jk_cl)
                ENDIF
                !
              ENDDO
              !
              !DO iw = 1, ir_obj%nfreq_b
                !
                !ir_knliw_cl(jx, iw) = w_dyn(iw_b, jbnd_cl, jk_cl, ibnd_cl, ik?)
                !
              !ENDDO
              !
              ! always use the weight of nscf k-point.
              weight_cl(jx) = 1.0_DP / REAL((nk1_cl * nk2_cl * nk3_cl), KIND=DP)
              !IF ((ik == lower_bnd) .AND. (ibnd == 1) .AND. (nn == 1))&
              !  WRITE(stdout, '(5x, a, ES20.10)')    'weight = ', weight_cl(jx) ! DEBUG
              IF (ABS(w_stat(jbnd_cl, ibnd_cl)) .GT. eps8 * inv_dos) THEN
                weight_cl(jx) = weight_cl(jx) * w_stat(jbnd_cl, ibnd_cl)
              ELSE
                weight_cl(jx) = zero
              ENDIF
              !
              IF ((jx == siz_ir_cl) .OR. (jstate == num_js2(ibnd, ik))) THEN
                ! To obtain expansion coefficients of green's function
                ! and kernel in the IR basis
                IF (positive_matsu) THEN
                  CALL fit_matsubara_f (ir_obj, ir_giw_cl, ir_gl_cl_d)
                  !CALL fit_matsubara_b (ir_obj, ir_knliw_cl, ir_knll_cl_d)
                ELSE
                  CALL fit_matsubara_f (ir_obj, ir_giw_cl, ir_gl_cl)
                  !CALL fit_matsubara_b (ir_obj, ir_knliw_cl, ir_knll_cl)
                ENDIF
                !
                IF (eps_cut_ir > zero) THEN
                  IF (positive_matsu) THEN
                    DO jx2 = 1, jx
                      gl_max = MAXVAL(ABS(ir_gl_cl_d(jx2, :)))
                      !knll_max = MAXVAL(ABS(ir_knll_cl_d(jx2, :)))
                      DO il = 1, ir_obj%size ! loop over the index of IR basis functions
                        gl_ratio = ABS(ir_gl_cl_d(jx2, il)) / gl_max
                        !knll_ratio = ABS(ir_knll_cl_d(jx2, il)) / knll_max
                        IF (gl_ratio < eps_cut_ir) THEN
                          ir_gl_cl_d(jx2, il) = czero
                        ENDIF
                        !IF (knll_ratio < eps_cut_ir) THEN
                        !  ir_knll_cl_d(jx2, il) = czero
                        !ENDIF
                      ENDDO
                    ENDDO
                  ELSE
                    DO jx2 = 1, jx
                      gl_max = MAXVAL(ABS(ir_gl_cl(jx2, :)))
                      !knll_max = MAXVAL(ABS(ir_knll_cl(jx2, :)))
                      DO il = 1, ir_obj%size ! loop over the index of IR basis functions
                        gl_ratio = ABS(ir_gl_cl(jx2, il)) / gl_max
                        !knll_ratio = ABS(ir_knll_cl(jx2, il)) / knll_max
                        IF (gl_ratio < eps_cut_ir) THEN
                          ir_gl_cl(jx2, il) = czero
                        ENDIF
                        !IF (knll_ratio < eps_cut_ir) THEN
                        !  ir_knll_cl(jx2, il) = czero
                        !ENDIF
                      ENDDO
                    ENDDO
                  ENDIF
                ENDIF
                !
                IF (positive_matsu) THEN
                  ! To obtain Green's function and kernel of imaginary time
                  ! from the corresponding expansion coefficients
                  CALL evaluate_tau (ir_obj, ir_gl_cl_d, ir_gtau_cl_d)
                  !!!!! Frequency dependence of Coulomb kernel !!!!!
                  ! HM: We leave the following comments for when
                  !     we incorporate the dynamic part of Coulomb kernel.
                  !CALL evaluate_tau (ir_obj, ir_knll_cl_d, ir_knltau_cl_d)
                  !
                  !DO it = 1, ir_obj%ntau
                  !  ir_cvltau_cl_d(1, it) = ir_gtau_cl_d(1, it) * ir_knltau_cl_d(1, it) / gtemp(itemp)
                  !ENDDO
                  !
                  ! To obtain expansion coefficients of the convolutions
                  ! in the IR basis
                  !CALL fit_tau (ir_obj, ir_cvltau_cl_d, ir_cvll_cl_d)
                  !
                  ! To reconstruct the convolutions of Matsubara freq.
                  ! from the coefficients in the IR basis
                  !CALL evaluate_matsubara_f (ir_obj, ir_cvll_cl_d, ir_cvliw_cl)
                  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                  !
                  DO iw = 1, ir_obj%nfreq_f
                    !
                    ! HM: IR cannot reproduce a constant function on frequencies,
                    !     so we have to treat the term multiplied by muc separately.
                    !     Be careful so as not to mistake the sign!
                    !     ir_gtau_cl_d(1:jx, 1) is an anomalous Green's function of imaginary time at tau = 0+
                    adeltai(iw, ibnd, ik) = adeltai(iw, ibnd, ik) &
                                            + SUM(weight_cl(1:jx) * REAL(ir_gtau_cl_d(1:jx, 1), KIND=DP)) &
                                            / gtemp(itemp)
                    !
                  ENDDO
                  !
                  ir_giw_cl(:, :) = czero
                  ir_gl_cl_d(:, :) = zero
                  ir_gtau_cl_d(:, :) = zero
                ELSE
                  ! To obtain Green's function and kernel of imaginary time
                  ! from the corresponding expansion coefficients
                  CALL evaluate_tau (ir_obj, ir_gl_cl, ir_gtau_cl)
                  !!!!! Frequency dependence of Coulomb kernel !!!!!
                  ! HM: We leave the following comments for when
                  !     we incorporate the dynamic part of Coulomb kernel.
                  !CALL evaluate_tau (ir_obj, ir_knll_cl, ir_knltau_cl)
                  !
                  !DO it = 1, ir_obj%ntau
                  !  ir_cvltau_cl(1, it) = ir_gtau_cl(1, it) * ir_knltau_cl(1, it) / gtemp(itemp)
                  !ENDDO
                  !
                  ! To obtain expansion coefficients of the convolutions
                  ! in the IR basis
                  !CALL fit_tau (ir_obj, ir_cvltau_cl, ir_cvll_cl)
                  !
                  ! To reconstruct the convolutions of Matsubara freq.
                  ! from the coefficients in the IR basis
                  !CALL evaluate_matsubara_f (ir_obj, ir_cvll_cl, ir_cvliw_cl)
                  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                  !
                  DO iw = 1, ir_obj%nfreq_f
                    !
                    ! HM: IR cannot reproduce a constant function on frequencies,
                    !     so we have to treat the term multiplied by muc separately.
                    !     Be careful so as not to mistake the sign!
                    !     ir_gtau_cl(1:jx, 1) is an anomalous Green's function of imaginary time at tau = 0+
                    adeltai(iw, ibnd, ik) = adeltai(iw, ibnd, ik) &
                                            + SUM(weight_cl(1:jx) * REAL(ir_gtau_cl(1:jx, 1), KIND=DP)) &
                                            / gtemp(itemp)
                    !
                  ENDDO
                  !
                  ir_giw_cl(:, :) = czero
                  ir_gl_cl(:, :) = czero
                  ir_gtau_cl(:, :) = czero
                ENDIF
              ENDIF
              !
            ENDDO
          ENDDO
        ENDIF ! ABS(ekfs(ibnd, ik) - ef0) < fsthick
      ENDDO
    ENDDO
    !
    IF ((iverbosity == 4) .AND. (iter == 1)) THEN
      ! collect contributions from all pools
      CALL mp_sum(knll_abs, inter_pool_comm)
      CALL print_kernell(itemp, ir_obj%size, knll_abs)
    ENDIF
    IF ((iverbosity == 4) .AND. (iter <= 5 .OR. iter == 15 .OR. MOD(iter, 10) == 0)) THEN
      ! collect contributions from all pools
      CALL mp_sum(gl_abs,   inter_pool_comm)
      CALL mp_sum(fl_abs,   inter_pool_comm)
      CALL print_gl(itemp, iter, ir_obj%size, gl_abs, fl_abs)
    ENDIF
    !WRITE(stdout, '(5x, a, I4)')    'nn = ', nn ! DEBUG
    !WRITE(1000+mpime, '(5x, a, I4)')    'nn = ', nn ! DEBUG
    !
    adeltai_cl(:, :) = zero
    !
    ir_giw_cl(:, :) = czero
    IF (positive_matsu) THEN
      ir_gl_cl_d(:, :) = zero
      ir_gtau_cl_d(:, :) = zero
    ELSE
      ir_gl_cl(:, :) = czero
      ir_gtau_cl(:, :) = czero
    ENDIF
    !
    linsidei = .TRUE.
    linsidei2 = .TRUE.
    linsidei3 = .TRUE.
    linsidej = .TRUE.
    linsidej2 = .TRUE.
    linsidej3 = .TRUE.
    !
    IF (num_states > 0) THEN
      DO istate = is_start, is_stop
        ik_cl = (istate - 1) / nbnd_cl + 1
        ibnd_cl = MOD(istate - 1, nbnd_cl) + 1
        !
        ! initilize linsidei
        linsidei = .TRUE.
        IF (ik_cl_to_fs(ik_cl) == 0) linsidei = .FALSE.
        ! initilize linsidei2
        linsidei2 = .TRUE.
        IF (((ibnd_cl - nbnd_offset) < 1) .OR. ((ibnd_cl - nbnd_offset) > nbndfs_all)) &
          linsidei2 = .FALSE.
        !
        ! initilize linsidei3
        linsidei3 = .TRUE.
        IF (linsidei .AND. linsidei2) THEN
          ibndfs = ibnd_kfs_all_to_kfs(ibnd_cl - nbnd_offset, ik_cl_to_fs(ik_cl))
          IF ((ibndfs < 1) .OR. (ibndfs > nbndfs)) THEN
            linsidei3 = .FALSE.
          ELSE
            ! Always use ekfs instead of ek_cl to determine
            ! if the state (ibnd_cl, ik_cl) is outside the window.
            IF (ABS(ekfs(ibndfs, ik_cl_to_fs(ik_cl)) - ef0) >= fsthick) &
              linsidei3 = .FALSE.
          ENDIF
        ENDIF
        !
        IF (linsidei .AND. linsidei2 .AND. linsidei3) THEN
          ! When the state (ibnd_cl, ik_cl) is inside of the fsthick window,
          ! no need to calculate adeltai_cl(:, istate)
          adeltai_cl(:, istate) = zero
        ELSEIF ((ek_cl(ibnd_cl, ik_cl) - ef0) >= emax_coulomb .OR. &
                (ek_cl(ibnd_cl, ik_cl) - ef0) <= emin_coulomb) THEN
          ! When the state (ibnd_cl, ik_cl) is outside of [emin_coulomb+ef0: emax_coulomb+ef0],
          ! no need to calculate adeltai_cl(:, istate)
          adeltai_cl(:, istate) = zero
        ELSE
          jstate = 0
          ! HM: jstate is the index combining jk_cl and jbnd_cl.
          !     Each time jstate is divisible by siz_ir_cl,
          !     calculate the convolutions and add the partial
          !     summation for each array.
          DO jk = 1, (nk1_cl * nk2_cl * nk3_cl)
            ! find the index of the irreducible k point corresponding to the current k point
            jk_cl = ik_bz_to_ibz_cl(jk)
            !
            ! initilize linsidej
            linsidej = .TRUE.
            IF (ik_cl_to_fs(jk_cl) == 0) linsidej = .FALSE.
            DO jbnd_cl = 1, nbnd_cl
              ! initilize linsidej2
              linsidej2 = .TRUE.
              IF (((jbnd_cl - nbnd_offset) < 1) .OR. ((jbnd_cl - nbnd_offset) > nbndfs_all)) &
                linsidej2 = .FALSE.
              !
              ! initilize linsidej3
              linsidej3 = .TRUE.
              IF (linsidej .AND. linsidej2) THEN
                jbndfs = ibnd_kfs_all_to_kfs(jbnd_cl - nbnd_offset, ik_cl_to_fs(jk_cl))
                IF ((jbndfs < 1) .OR. (jbndfs > nbndfs)) THEN
                  linsidej3 = .FALSE.
                ELSE
                  ! Always use ekfs instead of ek_cl to determine
                  ! if the state (jbnd_cl, jk_cl) is outside the window.
                  IF (ABS(ekfs(jbndfs, ik_cl_to_fs(jk_cl)) - ef0) >= fsthick) &
                    linsidej3 = .FALSE.
                ENDIF
              ENDIF
              !
              IF (.NOT.(linsidej .AND. linsidej2 .AND. linsidej3)) THEN
                IF (((ek_cl(jbnd_cl, jk_cl) - ef0) >= emax_coulomb) .OR. &
                    ((ek_cl(jbnd_cl, jk_cl) - ef0) <= emin_coulomb)) THEN
                  CYCLE
                ENDIF
              ENDIF
              !
              jstate = jstate + 1
              jx = MOD(jstate - 1, siz_ir_cl) + 1
              DO iw = 1, ir_obj%nfreq_f
                ! ir_obj%nfreq_f == nsiw(itemp)
                !
                ! FOR Fermionic functions
                !
                ! wsi(iw) == DBLE(ir_obj%freq_f(iw)) * pi * gtemp(itemp)
                !
                IF (linsidej .AND. linsidej2 .AND. linsidej3) THEN
                  esqrt = - one / ((wsi(iw) * aznormip(iw, jbndfs, ik_cl_to_fs(jk_cl)))**2.d0 &
                  + (ekfs(jbndfs, ik_cl_to_fs(jk_cl)) - muintr + ashiftip(iw, jbndfs, ik_cl_to_fs(jk_cl)))**2.d0 &
                  + (adeltaip(iw, jbndfs, ik_cl_to_fs(jk_cl)) * aznormip(iw, jbndfs, ik_cl_to_fs(jk_cl)))**2.d0)
                  !
                  ! anomalous Green's function
                  ir_giw_cl(jx, iw) = cone * esqrt * adeltaip(iw, jbndfs, ik_cl_to_fs(jk_cl)) &
                                          * aznormip(iw, jbndfs, ik_cl_to_fs(jk_cl))
                ELSE
                  esqrt = - one / (wsi(iw)**2.d0 &
                  + (ek_cl(jbnd_cl, jk_cl) - muintr)**2.d0 &
                  + (adeltaip_cl(iw, jbnd_cl, jk_cl))**2.d0)
                  !
                  ! anomalous Green's function
                  ir_giw_cl(jx, iw) = cone * esqrt * adeltaip_cl(iw, jbnd_cl, jk_cl)
                ENDIF
                !
              ENDDO
              !
              !DO iw = 1, ir_obj%nfreq_b
                !
                !ir_knliw_cl(1, iw) = w_dyn(iw_b, jbnd_cl, jk_cl, ibnd_cl, ik?)
                !
              !ENDDO
              !
              ! always use the weight of nscf k-point.
              weight_cl(jx) = 1.0_DP / REAL((nk1_cl * nk2_cl * nk3_cl), KIND=DP)
              !IF ((ik == lower_bnd) .AND. (ibnd == 1) .AND. (nn == 1))&
              !  WRITE(stdout, '(5x, a, ES20.10)')    'weight = ', weight_cl(jx) ! DEBUG
              IF (ABS(w_stat(jbnd_cl, ibnd_cl)) .GT. eps8 * inv_dos) THEN
                weight_cl(jx) = weight_cl(jx) * w_stat(jbnd_cl, ibnd_cl)
              ELSE
                weight_cl(jx) = zero
              ENDIF
              !
              IF ((jx == siz_ir_cl) .OR. (jstate == num_js3(istate))) THEN
                ! To obtain expansion coefficients of green's function
                ! and kernel in the IR basis
                IF (positive_matsu) THEN
                  CALL fit_matsubara_f (ir_obj, ir_giw_cl, ir_gl_cl_d)
                  !CALL fit_matsubara_b (ir_obj, ir_knliw_cl, ir_knll_cl_d)
                ELSE
                  CALL fit_matsubara_f (ir_obj, ir_giw_cl, ir_gl_cl)
                  !CALL fit_matsubara_b (ir_obj, ir_knliw_cl, ir_knll_cl)
                ENDIF
                !
                IF (eps_cut_ir > zero) THEN
                  IF (positive_matsu) THEN
                    DO jx2 = 1, jx
                      gl_max = MAXVAL(ABS(ir_gl_cl_d(jx2, :)))
                      !knll_max = MAXVAL(ABS(ir_knll_cl_d(jx2, :)))
                      DO il = 1, ir_obj%size ! loop over the index of IR basis functions
                        gl_ratio = ABS(ir_gl_cl_d(jx2, il)) / gl_max
                        !knll_ratio = ABS(ir_knll_cl_d(jx2, il)) / knll_max
                        IF (gl_ratio < eps_cut_ir) THEN
                          ir_gl_cl_d(jx2, il) = czero
                        ENDIF
                        !IF (knll_ratio < eps_cut_ir) THEN
                        !  ir_knll_cl_d(jx2, il) = czero
                        !ENDIF
                      ENDDO
                    ENDDO
                  ELSE
                    DO jx2 = 1, jx
                      gl_max = MAXVAL(ABS(ir_gl_cl(jx2, :)))
                      !knll_max = MAXVAL(ABS(ir_knll_cl(jx2, :)))
                      DO il = 1, ir_obj%size ! loop over the index of IR basis functions
                        gl_ratio = ABS(ir_gl_cl(jx2, il)) / gl_max
                        !knll_ratio = ABS(ir_knll_cl(jx2, il)) / knll_max
                        IF (gl_ratio < eps_cut_ir) THEN
                          ir_gl_cl(jx2, il) = czero
                        ENDIF
                        !IF (knll_ratio < eps_cut_ir) THEN
                        !  ir_knll_cl(jx2, il) = czero
                        !ENDIF
                      ENDDO
                    ENDDO
                  ENDIF
                ENDIF
                !
                IF (positive_matsu) THEN
                  ! To obtain green's function and kernel of imaginary time
                  ! from the corresponding expansion coefficients
                  CALL evaluate_tau (ir_obj, ir_gl_cl_d, ir_gtau_cl_d)
                  !!!!! Frequency dependence of Coulomb kernel !!!!!
                  ! HM: We leave the following comments for when
                  !     we incorporate the dynamic part of Coulomb kernel.
                  !CALL evaluate_tau (ir_obj, ir_knll_cl_d, ir_knltau_cl_d)
                  !
                  !DO it = 1, ir_obj%ntau
                  !  ir_cvltau_cl_d(1, it) = ir_gtau_cl_d(1, it) * ir_knltau_cl_d(1, it) / gtemp(itemp)
                  !ENDDO
                  !
                  ! To obtain expansion coefficients of the convolutions
                  ! in the IR basis
                  !CALL fit_tau (ir_obj, ir_cvltau_cl_d, ir_cvll_cl_d)
                  !
                  ! To reconstruct the convolutions of Matsubara freq.
                  ! from the coefficients in the IR basis
                  !CALL evaluate_matsubara_f (ir_obj, ir_cvll_cl_d, ir_cvliw_cl)
                  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                  !
                  DO iw = 1, ir_obj%nfreq_f
                    !
                    ! HM: IR cannot reproduce a constant function on frequencies,
                    !     so we have to treat the term multiplied by muc separately.
                    !     Be careful so as not to mistake the sign!
                    !     ir_gtau_cl_d(1:jx, 1) is an anomalous Green's function of imaginary time at tau = 0+
                    adeltai_cl(iw, istate) = adeltai_cl(iw, istate) &
                                          + SUM(weight_cl(1:jx) * REAL(ir_gtau_cl_d(1:jx, 1), KIND=DP)) &
                                          / gtemp(itemp)
                    !
                  ENDDO
                  !
                  ir_giw_cl(:, :) = czero
                  ir_gl_cl_d(:, :) = zero
                  ir_gtau_cl_d(:, :) = zero
                ELSE
                  ! To obtain green's function and kernel of imaginary time
                  ! from the corresponding expansion coefficients
                  CALL evaluate_tau (ir_obj, ir_gl_cl, ir_gtau_cl)
                  !!!!! Frequency dependence of Coulomb kernel !!!!!
                  ! HM: We leave the following comments for when
                  !     we incorporate the dynamic part of Coulomb kernel.
                  !CALL evaluate_tau (ir_obj, ir_knll_cl, ir_knltau_cl)
                  !
                  !DO it = 1, ir_obj%ntau
                  !  ir_cvltau_cl(1, it) = ir_gtau_cl(1, it) * ir_knltau_cl(1, it) / gtemp(itemp)
                  !ENDDO
                  !
                  ! To obtain expansion coefficients of the convolutions
                  ! in the IR basis
                  !CALL fit_tau (ir_obj, ir_cvltau_cl, ir_cvll_cl)
                  !
                  ! To reconstruct the convolutions of Matsubara freq.
                  ! from the coefficients in the IR basis
                  !CALL evaluate_matsubara_f (ir_obj, ir_cvll_cl, ir_cvliw_cl)
                  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                  !
                  DO iw = 1, ir_obj%nfreq_f
                    !
                    ! HM: IR cannot reproduce a constant function on frequencies,
                    !     so we have to treat the term multiplied by muc separately.
                    !     Be careful so as not to mistake the sign!
                    !     ir_gtau_cl(1:jx, 1) is an anomalous Green's function of imaginary time at tau = 0+
                    adeltai_cl(iw, istate) = adeltai_cl(iw, istate) &
                                          + SUM(weight_cl(1:jx) * REAL(ir_gtau_cl(1:jx, 1), KIND=DP)) &
                                          / gtemp(itemp)
                    !
                  ENDDO
                  !
                  ir_giw_cl(:, :) = czero
                  ir_gl_cl(:, :) = czero
                  ir_gtau_cl(:, :) = czero
                ENDIF
              ENDIF
            ENDDO
          ENDDO
        ENDIF ! linsidei .AND. linsidei2
      ENDDO
    ENDIF
    !
    DO iw = 1, nsiw(itemp) ! loop over omega
      IF (num_states > 0) THEN
        DO istate = is_start, is_stop
          ik_cl = (istate - 1) / nbnd_cl + 1
          ibnd_cl = MOD(istate - 1, nbnd_cl) + 1
          ! assume aznormi is 1.0 outside the fsthick window
          adeltai_cl(iw, istate) = gtemp(itemp) * adeltai_cl(iw, istate) ! / 1.0_DP
        ENDDO ! istate
      ENDIF
    ENDDO ! iw
    !
    DEALLOCATE(inv_wsi, STAT = ierr)
    IF (ierr /= 0) CALL errore('sum_eliash_aniso_iaxis_fbw_ir_coul', 'Error deallocating inv_wsi', 1)
    ! remove memory allocated for inv_wsi
    imelt = nsiw(itemp)
    CALL mem_size_eliashberg(2, -imelt)
    !
    RETURN
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE sum_eliash_aniso_iaxis_fbw_ir_coul
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    SUBROUTINE analytic_cont_aniso(itemp, iter, conv)
    !-----------------------------------------------------------------------
    !
    ! This routine does the analytic continuation of the anisotropic Eliashberg equations
    ! from the imaginary-axis to the real axis
    ! reference F. Marsiglio, M. Schossmann, and J. Carbotte, Phys. Rev. B 37, 4965 (1988)
    !
    USE kinds,         ONLY : DP
    USE modes,         ONLY : nmodes
    USE global_var,    ONLY : wqf, wf, gtemp
    USE input,         ONLY : nqstep, degaussq, nsiter, conv_thr_racon, fsthick, &
                              lpade, eps_acoustic
    USE supercond_common,     ONLY : ixkqf, ixqfs, nqfs, wkfs, w0g, ekfs, nkfs, nbndfs, dosef, ef0, &
                              nsw, dwsph, ws, wsph, gap0, agap, gp, gm, adsumi, azsumi, a2fij, &
                              g2, lacon_fly, delta, znorm, adelta, adeltap, aznorm, aznormp
    USE supercond,     ONLY : gamma_acont
    USE ep_constants,  ONLY : ci, zero, one, czero, cone
    USE ep_constants,  ONLY : pi
    USE io_global,     ONLY : stdout, ionode_id
    USE mp_global,     ONLY : inter_pool_comm
    USE mp_world,      ONLY : mpime
    USE mp,            ONLY : mp_bcast, mp_barrier, mp_sum
    USE parallelism,   ONLY : fkbounds
    USE low_lvl,       ONLY : mem_size_eliashberg, memlt_eliashberg
    !
    IMPLICIT NONE
    !
    LOGICAL, INTENT(inout) :: conv
    !! True if the calculation is converged
    INTEGER, INTENT(in) :: itemp
    !! Counter on temperature index
    INTEGER, INTENT(in) :: iter
    !! Counter on the iteration number
    !
    ! Local variables
    INTEGER :: i, iw, iwp
    !! Counter on frequency real-axis
    INTEGER :: ik
    !! Counter on k-points
    INTEGER :: iq
    !! Counter on q-points for which k+sign*q is within the Fermi shell
    INTEGER :: iq0
    !! Index of iq on full q-mesh
    INTEGER :: lower_bnd, upper_bnd
    !! Lower/upper bound index after k paral
    INTEGER :: nks
    !! Number of k points per pool
    INTEGER :: ibnd
    !! Counter on bands at k
    INTEGER :: jbnd
    !! Counter on bands at k+q
    INTEGER :: iwph
    !! Counter on frequency in a2f
    INTEGER :: imode
    !! Counter on phonon modes
    INTEGER(8) :: imelt
    !! Counter memory
    INTEGER :: ierr
    !! Error status
    !
    REAL(KIND = DP) :: rgammap
    !! - bose_einstein(w') - fermi_dirac(w + w')
    REAL(KIND = DP) :: rgammam
    !!   bose_einstein(w') + fermi_dirac(-w + w')
    REAL(KIND = DP) :: absdelta, reldelta, errdelta
    !! Errors in supercond. gap
    REAL(KIND = DP) :: weight
    !! Factor in supercond. equations
    REAL(KIND = DP) :: inv_degaussq
    !! 1.0/degaussq. Defined for efficiency reasons
    REAL(KIND = DP) :: inv_dos
    !! Invese dos inv_dos = 1/dosef. Defined for efficiency reason
    REAL(KIND = DP) :: inv_ws
    !! Invese real frequency inv_ws = 1/ws. Defined for efficiency reason
    REAL(KIND = DP), EXTERNAL :: w0gauss
    !! The derivative of wgauss:  an approximation to the delta function
    !!RM - updated
    REAL(KIND = DP) :: a2f_tmp(nqstep)
    !! Temporary variable for Eliashberg spectral function
    !
    REAL(KIND = DP) :: root_im
    !! Temporary variable
    COMPLEX(KIND = DP) :: az2, ad2, esqrt, root
    !! Temporary variables
    COMPLEX(KIND = DP), ALLOCATABLE, SAVE :: deltaold(:)
    !! supercond. gap from previous iteration
    !
    inv_degaussq = one / degaussq
    !
    inv_dos = one / dosef
    !
    CALL fkbounds(nkfs, lower_bnd, upper_bnd)
    !
    nks = upper_bnd - lower_bnd + 1
    !
    IF (iter == 1) THEN
      !
      ! RM - adelta, aznorm are defined per k-points per pool
      !
      ! get the size of required allocated memory for
      ! delta, znorm, deltaold, adelta, aznorm, adeltap, aznormp, gp, gm
      IF (lpade) THEN
        imelt = 2 * (1 + 2 * nbndfs * nkfs) * nsw + 2 * nqstep * nsw
      ELSE
        imelt = 2 * (3 + 2 * nbndfs * nks + 2 * nbndfs * nkfs) * nsw + 2 * nqstep * nsw
      ENDIF
      CALL mem_size_eliashberg(2, imelt)
      !
      IF (.NOT. lpade) THEN
        ALLOCATE(delta(nsw), STAT = ierr)
        IF (ierr /= 0) CALL errore('analytic_cont_aniso', 'Error allocating delta', 1)
        ALLOCATE(znorm(nsw), STAT = ierr)
        IF (ierr /= 0) CALL errore('analytic_cont_aniso', 'Error allocating znorm', 1)
        ALLOCATE(adelta(nsw, nbndfs, lower_bnd:upper_bnd), STAT = ierr)
        IF (ierr /= 0) CALL errore('analytic_cont_aniso', 'Error allocating adelta', 1)
        ALLOCATE(aznorm(nsw, nbndfs, lower_bnd:upper_bnd), STAT = ierr)
        IF (ierr /= 0) CALL errore('analytic_cont_aniso', 'Error allocating aznorm', 1)
      ENDIF
      ALLOCATE(adeltap(nsw, nbndfs, nkfs), STAT = ierr)
      IF (ierr /= 0) CALL errore('analytic_cont_aniso', 'Error allocating adeltap', 1)
      ALLOCATE(aznormp(nsw, nbndfs, nkfs), STAT = ierr)
      IF (ierr /= 0) CALL errore('analytic_cont_aniso', 'Error allocating aznormp', 1)
      ALLOCATE(deltaold(nsw), STAT = ierr)
      IF (ierr /= 0) CALL errore('analytic_cont_aniso', 'Error allocating deltaold', 1)
      ALLOCATE(gp(nsw, nqstep), STAT = ierr)
      IF (ierr /= 0) CALL errore('analytic_cont_aniso', 'Error allocating gp', 1)
      ALLOCATE(gm(nsw, nqstep), STAT = ierr)
      IF (ierr /= 0) CALL errore('analytic_cont_aniso', 'Error allocating gm', 1)
      !
      adeltap(:, :, :) = czero
      aznormp(:, :, :) = czero
      deltaold(:) = czero
      IF (lpade) THEN
        adeltap(:, :, lower_bnd:upper_bnd) = adelta(:, :, lower_bnd:upper_bnd)
        deltaold(:) = delta(:)
      ELSE
        DO iw = 1, nsw
          adeltap(iw, :, lower_bnd:upper_bnd) = agap(:, lower_bnd:upper_bnd)
        ENDDO
        deltaold(:) = gap0
      ENDIF
      ! collect contributions from all pools
      ! (only at iter=1 since for iter>1 it is done in eliashberg_aniso_iaxis)
      CALL mp_sum(adeltap, inter_pool_comm)
      CALL mp_barrier(inter_pool_comm)
      !
      aznormp(:, :, lower_bnd:upper_bnd) = cone
      !
      !! Eq.(28) in Margine and Giustino, PRB 87, 024505 (2013)
      DO iwp = 1, nqstep ! loop over omega_prime
        DO iw = 1, nsw ! loop over omega
          CALL gamma_acont(ws(iw), ws(iwp), gtemp(itemp), rgammap, rgammam)
          gp(iw, iwp) = rgammap
          gm(iw, iwp) = rgammam
        ENDDO
      ENDDO
      CALL kernel_aniso_analytic_cont(itemp)
      CALL memlt_eliashberg(itemp, 'acon')
      IF (.NOT. lacon_fly) CALL evaluate_a2fij()
    ENDIF ! iter
    !
    ! RM - calculate esqrt (esqrt is stored as aznormp to avoid allocation of a new large array)
    ! Z(ik, ibnd, w \pm w') / sqrt{[Z(ik, ibnd, w \pm w')]^2 * [(w \pm w')^2] - [D(ik, ibnd, w \pm w')]^2]}
    DO ik = lower_bnd, upper_bnd
      DO ibnd = 1, nbndfs
        IF (ABS(ekfs(ibnd, ik) - ef0) < fsthick) THEN
          !
          ! nvfortran bug workaround https://gitlab.com/QEF/q-e/-/issues/593
#if defined(__NVCOMPILER) && ( __NVCOMPILER_MAJOR__ > 23 || ( __NVCOMPILER_MAJOR__ == 23 &&  __NVCOMPILER_MINOR__ >= 3) )
          !pgi$l novector
#endif
          DO iw = 1, nsw ! loop over omega
            az2 = aznormp(iw, ibnd, ik) * aznormp(iw, ibnd, ik)
            ad2 = adeltap(iw, ibnd, ik) * adeltap(iw, ibnd, ik)
            root = SQRT(az2 * (ws(iw) * ws(iw) - ad2))
            root_im = AIMAG(root)
            IF (root_im < zero) &
              root = CONJG(root)
            esqrt = aznormp(iw, ibnd, ik) / root
            aznormp(iw, ibnd, ik) = esqrt
          ENDDO
        ENDIF ! fsthick
      ENDDO ! ibnd
    ENDDO ! ik
    !
    ! collect contributions from all pools
    CALL mp_sum(aznormp, inter_pool_comm)
    CALL mp_barrier(inter_pool_comm)
    !
    adelta(:, :, :) = czero
    aznorm(:, :, :) = czero
    a2f_tmp(:) = zero
    DO ik = lower_bnd, upper_bnd
      DO ibnd = 1, nbndfs
        IF (ABS(ekfs(ibnd, ik) - ef0) < fsthick) THEN
          DO iq = 1, nqfs(ik)
            ! iq0 - index of q-point on the full q-mesh
            iq0 = ixqfs(ik, iq)
            DO jbnd = 1, nbndfs
              IF (ABS(ekfs(jbnd, ixkqf(ik, iq0)) - ef0) < fsthick) THEN
                !
                !!RM - updated
                IF (lacon_fly) THEN ! evaluate a2fij on the fly
                  a2f_tmp(:) = zero
                  DO imode = 1, nmodes
                    IF (wf(imode, iq0) > eps_acoustic) THEN
                      DO iwph = 1, nqstep
                        weight = w0gauss((wsph(iwph) - wf(imode, iq0)) * inv_degaussq, 0) * inv_degaussq
                        a2f_tmp(iwph) =  a2f_tmp(iwph) + weight * dosef * g2(ik, iq, ibnd, jbnd, imode)
                      ENDDO ! iwph
                    ENDIF ! wf
                  ENDDO ! imode
                ELSE
                  a2f_tmp(:) = a2fij(:, jbnd, iq, ibnd, ik)
                ENDIF ! lacon_fly
                !
                weight = wqf(iq) * w0g(jbnd, ixkqf(ik, iq0)) * inv_dos
                DO iw = 1, nsw ! loop over omega
                  DO iwp = 1, nqstep ! loop over omega_prime
                    !
                    i = iw + iwp
                    IF (i <= nsw) THEN
                      esqrt = weight * gp(iw, iwp) * a2f_tmp(iwp) * aznormp(i, jbnd, ixkqf(ik, iq0))
                      !
                      aznorm(iw, ibnd, ik) = aznorm(iw, ibnd, ik) - esqrt * ws(i)
                      adelta(iw, ibnd, ik) = adelta(iw, ibnd, ik) - esqrt * adeltap(i, jbnd, ixkqf(ik, iq0))
                    ENDIF
                    !
                    i = ABS(iw - iwp)
                    IF (i > 0) THEN
                      esqrt = weight * gm(iw, iwp) * a2f_tmp(iwp) * aznormp(i, jbnd, ixkqf(ik, iq0))
                      !
                      aznorm(iw, ibnd, ik) = aznorm(iw, ibnd, ik) + esqrt * ws(i) * SIGN(1, iw - iwp)
                      adelta(iw, ibnd, ik) = adelta(iw, ibnd, ik) + esqrt * adeltap(i, jbnd, ixkqf(ik, iq0))
                    ENDIF
                  ENDDO ! iwp
                ENDDO ! iw
              ENDIF ! fsthick
            ENDDO ! jbnd
          ENDDO ! iq
          ! Eqs.(26)-(27) in Margine and Giustino, PRB 87, 024505 (2013)
          DO iw = 1, nsw ! loop over omega
            aznorm(iw, ibnd, ik) = - gtemp(itemp) * azsumi(iw, ibnd, ik) + ci * aznorm(iw, ibnd, ik) * dwsph
            adelta(iw, ibnd, ik) =   gtemp(itemp) * adsumi(iw, ibnd, ik) + ci * adelta(iw, ibnd, ik) * dwsph
          ENDDO ! iw
        ENDIF ! fsthick
      ENDDO ! ibnd
    ENDDO ! ik
    !
    delta(:) = czero
    znorm(:) = czero
    DO iw = 1, nsw ! loop over omega
      inv_ws = one / ws(iw)
      DO ik = lower_bnd, upper_bnd
        DO ibnd = 1, nbndfs
          IF (ABS(ekfs(ibnd, ik) - ef0) < fsthick) THEN
            weight = 0.5d0 * wkfs(ik) * w0g(ibnd, ik) * inv_dos
            znorm(iw) = znorm(iw) + weight * aznorm(iw, ibnd, ik)
            delta(iw) = delta(iw) + weight * adelta(iw, ibnd, ik)
            ! Eqs.(26)-(27) in Margine and Giustino, PRB 87, 024505 (2013)
            aznorm(iw, ibnd, ik) = 1.d0 + pi * aznorm(iw, ibnd, ik) * inv_ws
            adelta(iw, ibnd, ik) = pi * adelta(iw, ibnd, ik) / aznorm(iw, ibnd, ik)
          ENDIF
        ENDDO ! ibnd
      ENDDO ! ik
    ENDDO ! iw
    !
    ! collect contributions from all pools
    CALL mp_sum(delta, inter_pool_comm)
    CALL mp_sum(znorm, inter_pool_comm)
    CALL mp_barrier(inter_pool_comm)
    !
    ! RM - update aznormp for next iteration
    ! Make everything 0 except the range of k-points we are working on
    ! In principle aznormp can be mixed as done for adeltap
    aznormp(:, :, :) = czero
    aznormp(:, :, lower_bnd:upper_bnd) = aznorm(:, :, lower_bnd:upper_bnd)
    !
    IF (mpime == ionode_id) THEN
      !
      absdelta = zero
      reldelta = zero
      DO iw = 1, nsw
        inv_ws = one / ws(iw)
        znorm(iw) = 1.0d0 + pi * znorm(iw) * inv_ws
        delta(iw) = pi * delta(iw) / znorm(iw)
        reldelta = reldelta + ABS(delta(iw) - deltaold(iw))
        absdelta = absdelta + ABS(delta(iw))
      ENDDO
      errdelta = reldelta / absdelta
      deltaold(:) = delta(:)
      !
      IF (iter == 1) &
        WRITE(stdout, '(5x, a)') '   iter      ethr         Re[znorm]   Re[delta] [meV]'
        WRITE(stdout, '(5x, i6, 3ES15.6)') iter, errdelta, REAL(znorm(1)), REAL(delta(1)) * 1000.d0
!      WRITE(stdout, '(5x, a, i6, a, ES15.6, a, ES15.6, a, ES15.6)') 'iter = ', iter, &
!                    '   ethr = ', errdelta, '   Re[znorm(1)] = ', REAL(znorm(1)), &
!                    '   Re[delta(1)] = ', REAL(delta(1))
      !
      IF (errdelta < conv_thr_racon) conv = .TRUE.
      !
      IF (conv) THEN
        WRITE(stdout, '(5x, a, i6)') 'Convergence was reached in nsiter = ', iter
        WRITE(stdout, '(a)') ' '
      ELSEIF (.NOT. conv .AND. iter == nsiter) THEN
        WRITE(stdout, '(5x, a, i6)') 'Convergence was not reached in nsiter = ', iter
        WRITE(stdout, '(5x, a)') 'Increase nsiter or reduce conv_thr_racon'
        WRITE(stdout, '(a)') ' '
      ENDIF
    ENDIF !mpime
    CALL mp_bcast(conv, ionode_id, inter_pool_comm)
    CALL mp_barrier(inter_pool_comm)
    !
    IF (conv .OR. iter == nsiter) THEN
      DEALLOCATE(deltaold, STAT = ierr)
      IF (ierr /= 0) CALL errore('analytic_cont_aniso', 'Error deallocating deltaold', 1)
      !
      ! remove memory allocated for deltaold
      imelt = nsw
      CALL mem_size_eliashberg(2, -imelt)
      !
    ENDIF
    !
    RETURN
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE analytic_cont_aniso
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    SUBROUTINE pade_cont_aniso(itemp, N)
    !-----------------------------------------------------------------------
    !!
    !! This routine uses pade approximants to continue the anisotropic Eliashberg equations
    !! from the imaginary-axis to the real-axis
    !!
    !! SH: Adopted for the full-bandwidth calculations (Nov 2021).
    !!
    !
    USE kinds,         ONLY : DP
    USE input,         ONLY : fsthick, fbw, positive_matsu
    USE supercond_common,     ONLY : nsw, ws, nsiw, wsi, delta, &
                              znorm, adelta, aznorm, adeltai, aznormi, &
                              wkfs, dosef, w0g, nkfs, nbndfs, ef0, ekfs, &
                              shift, ashift, ashifti
    USE utilities,     ONLY : pade_coeff, pade_eval
    USE ep_constants,  ONLY : cone, ci, zero, czero, one
    USE io_global,     ONLY : stdout, ionode_id
    USE mp_global,     ONLY : inter_pool_comm
    USE mp_world,      ONLY : mpime
    USE mp,            ONLY : mp_barrier, mp_sum
    USE parallelism,   ONLY : fkbounds
    USE low_lvl,       ONLY : mem_size_eliashberg
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: itemp
    !! Counter on temperature index
    INTEGER, INTENT(in) :: N
    !! Nr. of frequency points in the Pade approx
    !
    ! Local variable
    !
    INTEGER :: iw
    !! Counter on frequency imag- and real-axis
    INTEGER :: ik
    !! Counter on k-points
    INTEGER :: lower_bnd, upper_bnd
    !! Lower/upper bound index after k paral
    INTEGER :: nks
    !! Number of k points per pool
    INTEGER :: ibnd
    !! Counter on bands
    INTEGER(8) :: imelt
    !! Counter memory
    INTEGER :: ierr
    !! Error status
    !
    REAL(KIND = DP) :: weight
    !! Temporary variable
    REAL(KIND = DP) :: inv_dos
    !! Invese dos inv_dos = 1/dosef. Defined for efficiency reason
    !
    ! arrays used in pade_coeff and pade_eval
    COMPLEX(KIND = DP) :: omega
    !! frequency real-axis
    COMPLEX(KIND = DP) :: padapp
    !! znorm or delta on real-axis after pade_eval
    COMPLEX(KIND = DP), ALLOCATABLE :: a(:)
    !! a - pade coeff for deltai
    COMPLEX(KIND = DP), ALLOCATABLE :: b(:)
    !! b - pade coeff for znormi
    COMPLEX(KIND = DP), ALLOCATABLE :: c(:)
    !! c - pade coeff for shifti
    COMPLEX(KIND = DP), ALLOCATABLE :: z(:)
    !! z - frequency imag-axis
    COMPLEX(KIND = DP), ALLOCATABLE :: u(:)
    !! u - deltai
    COMPLEX(KIND = DP), ALLOCATABLE :: v(:)
    !! v - znormi
    COMPLEX(KIND = DP), ALLOCATABLE :: w(:)
    !! w - shifti
    !
    inv_dos = one / dosef
    !
    CALL fkbounds(nkfs, lower_bnd, upper_bnd)
    !
    nks = upper_bnd - lower_bnd + 1
    !
    ! RM - adelta, aznorm are defined per k-points per pool
    ! get the size of required allocated memory for
    ! a, b, z, u, v, delta, znorm, adelta, aznorm
    imelt = 2 * 5 * N + 2 * (2 + 2 * nbndfs * nks) * nsw
    CALL mem_size_eliashberg(2, imelt)
    !
    ALLOCATE(delta(nsw), STAT = ierr)
    IF (ierr /= 0) CALL errore('pade_cont_aniso', 'Error allocating delta', 1)
    ALLOCATE(znorm(nsw), STAT = ierr)
    IF (ierr /= 0) CALL errore('pade_cont_aniso', 'Error allocating znorm', 1)
    ALLOCATE(adelta(nsw, nbndfs, lower_bnd:upper_bnd), STAT = ierr)
    IF (ierr /= 0) CALL errore('pade_cont_aniso', 'Error allocating adelta', 1)
    ALLOCATE(aznorm(nsw, nbndfs, lower_bnd:upper_bnd), STAT = ierr)
    IF (ierr /= 0) CALL errore('pade_cont_aniso', 'Error allocating aznorm', 1)
    ALLOCATE(a(N), STAT = ierr)
    IF (ierr /= 0) CALL errore('pade_cont_aniso', 'Error allocating a', 1)
    ALLOCATE(b(N), STAT = ierr)
    IF (ierr /= 0) CALL errore('pade_cont_aniso', 'Error allocating b', 1)
    ALLOCATE(z(N), STAT = ierr)
    IF (ierr /= 0) CALL errore('pade_cont_aniso', 'Error allocating z', 1)
    ALLOCATE(u(N), STAT = ierr)
    IF (ierr /= 0) CALL errore('pade_cont_aniso', 'Error allocating u', 1)
    ALLOCATE(v(N), STAT = ierr)
    IF (ierr /= 0) CALL errore('pade_cont_aniso', 'Error allocating v', 1)
    delta(:) = czero
    znorm(:) = czero
    adelta(:, :, :) = czero
    aznorm(:, :, :) = czero
    a(:) = czero
    b(:) = czero
    z(:) = czero
    u(:) = czero
    v(:) = czero
    !
    IF (fbw) THEN
      ! get the size of required allocated memory for
      ! c, w, shift, ashift
      imelt = 2 * 2 * N + 2 * (1 + 1 * nbndfs * nks) * nsw
      CALL mem_size_eliashberg(2, imelt)
      !
      ALLOCATE(shift(nsw), STAT = ierr)
      IF (ierr /= 0) CALL errore('pade_cont_aniso_fbw', 'Error allocating shift', 1)
      ALLOCATE(ashift(nsw, nbndfs, lower_bnd:upper_bnd), STAT = ierr)
      IF (ierr /= 0) CALL errore('pade_cont_aniso_fbw', 'Error allocating ashift', 1)
      ALLOCATE(c(N), STAT = ierr)
      IF (ierr /= 0) CALL errore('pade_cont_aniso_fbw', 'Error allocating c', 1)
      ALLOCATE(w(N), STAT = ierr)
      IF (ierr /= 0) CALL errore('pade_cont_aniso_fbw', 'Error allocating w', 1)
      shift(:) = czero
      ashift(:, :, :) = czero
      c(:) = czero
      w(:) = czero
    ENDIF
      !
    IF (fbw) THEN
      DO ik = lower_bnd, upper_bnd
        DO ibnd = 1, nbndfs
          IF (ABS(ekfs(ibnd, ik) - ef0) < fsthick) THEN
            weight = 0.5d0 * wkfs(ik) * w0g(ibnd, ik) * inv_dos
            IF (positive_matsu) THEN
              DO iw = 1, N
                z(iw) = ci * wsi(iw)
                u(iw) = cone * adeltai(iw, ibnd, ik)
                v(iw) = cone * aznormi(iw, ibnd, ik)
                w(iw) = cone * ashifti(iw, ibnd, ik)
              ENDDO
            ELSE
              DO iw = 1, N
                !
                !! if positive_matsu =.FALSE.
                !! use the data only on the positive Matsubara frequencies.
                !
                z(iw) = ci * wsi(nsiw(itemp)/2 + iw)
                u(iw) = cone * adeltai(nsiw(itemp)/2 + iw, ibnd, ik)
                v(iw) = cone * aznormi(nsiw(itemp)/2 + iw, ibnd, ik)
                w(iw) = cone * ashifti(nsiw(itemp)/2 + iw, ibnd, ik)
              ENDDO
            ENDIF ! positive_matsu
            CALL pade_coeff(N, z, u, a)
            CALL pade_coeff(N, z, v, b)
            CALL pade_coeff(N, z, w, c)
            DO iw = 1, nsw
              omega = cone * ws(iw)
              CALL pade_eval(N, z, a, omega, padapp)
              adelta(iw, ibnd, ik) = padapp
              CALL pade_eval(N, z, b, omega, padapp)
              aznorm(iw, ibnd, ik) = padapp
              CALL pade_eval(N, z, c, omega, padapp)
              ashift(iw, ibnd, ik) = padapp
              !
              znorm(iw) = znorm(iw) + weight * aznorm(iw, ibnd, ik)
              delta(iw) = delta(iw) + weight * adelta(iw, ibnd, ik)
              shift(iw) = shift(iw) + weight * ashift(iw, ibnd, ik)
            ENDDO
          ENDIF
        ENDDO ! ibnd
      ENDDO ! ik
      ! collect contributions from all pools
      CALL mp_sum(znorm,  inter_pool_comm)
      CALL mp_sum(delta,  inter_pool_comm)
      CALL mp_sum(shift,  inter_pool_comm)
      CALL mp_barrier(inter_pool_comm)
    ELSE
      DO ik = lower_bnd, upper_bnd
        DO ibnd = 1, nbndfs
          IF (ABS(ekfs(ibnd, ik) - ef0) < fsthick) THEN
            weight = 0.5d0 * wkfs(ik) * w0g(ibnd, ik) * inv_dos
            DO iw = 1, N
              z(iw) = ci * wsi(iw)
              u(iw) = cone * adeltai(iw, ibnd, ik)
              v(iw) = cone * aznormi(iw, ibnd, ik)
            ENDDO
            CALL pade_coeff(N, z, u, a)
            CALL pade_coeff(N, z, v, b)
            DO iw = 1, nsw
              omega = cone * ws(iw)
              CALL pade_eval(N, z, a, omega, padapp)
              adelta(iw, ibnd, ik) = padapp
              CALL pade_eval(N, z, b, omega, padapp)
              aznorm(iw, ibnd, ik) = padapp
              !
              znorm(iw) = znorm(iw) + weight * aznorm(iw, ibnd, ik)
              delta(iw) = delta(iw) + weight * adelta(iw, ibnd, ik)
            ENDDO
          ENDIF
        ENDDO ! ibnd
      ENDDO ! ik
      ! collect contributions from all pools
      CALL mp_sum(znorm,  inter_pool_comm)
      CALL mp_sum(delta,  inter_pool_comm)
      CALL mp_barrier(inter_pool_comm)
    ENDIF ! fbw
    !
    IF (mpime == ionode_id) THEN
      IF (fbw) THEN
        WRITE(stdout, '(5x, a)') '   pade    Re[znorm]   Re[delta] [meV]   Re[shift] [meV]'
        WRITE(stdout, '(5x, i6, 3ES15.6)') N, REAL(znorm(1)), REAL(delta(1)) * 1000.d0, REAL(shift(1)) * 1000.d0
      ELSE
        WRITE(stdout, '(5x, a)') '   pade    Re[znorm]   Re[delta] [meV]'
        WRITE(stdout, '(5x, i6, 2ES15.6)') N, REAL(znorm(1)), REAL(delta(1)) * 1000.d0
        !
        !      WRITE(stdout, '(5x, a, i6, a, ES15.6, a, ES15.6)') 'pade = ', N, &
        !                    '   Re[znorm(1)] = ', REAL(znorm(1)), &
        !                    '   Re[delta(1)] = ', REAL(delta(1))
      ENDIF
      !
      WRITE(stdout, '(a)') ' '
      WRITE(stdout, '(5x, a, i6, a)') 'Convergence was reached for N = ', N, ' Pade approximants'
      WRITE(stdout, '(a)') ' '
    ENDIF
    !
    DEALLOCATE(a, STAT = ierr)
    IF (ierr /= 0) CALL errore('pade_cont_aniso', 'Error deallocating a', 1)
    DEALLOCATE(b, STAT = ierr)
    IF (ierr /= 0) CALL errore('pade_cont_aniso', 'Error deallocating b', 1)
    DEALLOCATE(z, STAT = ierr)
    IF (ierr /= 0) CALL errore('pade_cont_aniso', 'Error deallocating z', 1)
    DEALLOCATE(u, STAT = ierr)
    IF (ierr /= 0) CALL errore('pade_cont_aniso', 'Error deallocating u', 1)
    DEALLOCATE(v, STAT = ierr)
    IF (ierr /= 0) CALL errore('pade_cont_aniso', 'Error deallocating v', 1)
    ! remove memory allocated for a, b, z, u, v
    imelt = 2 * 5 * N
    CALL mem_size_eliashberg(2, -imelt)
    !
    IF (fbw) THEN
      DEALLOCATE(c, STAT = ierr)
      IF (ierr /= 0) CALL errore('pade_cont_aniso', 'Error deallocating c', 1)
      DEALLOCATE(w, STAT = ierr)
      IF (ierr /= 0) CALL errore('pade_cont_aniso', 'Error deallocating w', 1)
      ! remove memory allocated for c, w
      imelt = 2 * 2 * N
      CALL mem_size_eliashberg(2, -imelt)
    ENDIF
    !
    RETURN
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE pade_cont_aniso
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    SUBROUTINE kernel_aniso_iaxis(itemp)
    !-----------------------------------------------------------------------
    !!
    !! Compute kernels K_{+}(ik, iq, ibnd, jbnd; n, n', T) and
    !! K_{-}(ik, iq, ibnd, jbnd; n, n', T) and store them in memory
    !!
    !! SH: Modified to allow for sparse sampling of Matsubara freq. (Nov 2021).
    !! HM: It is assumed that this subroutine is not called when using sparse-ir sampling
    !!
    USE kinds,         ONLY : DP
    USE input,         ONLY : fsthick
    USE global_var,    ONLY : gtemp
    USE supercond_common,     ONLY : nkfs, nbndfs, nsiw, akeri, ekfs, ef0, ixkqf, ixqfs, nqfs, wsn
    USE ep_constants,  ONLY : zero
    USE ep_constants,  ONLY : pi
    USE parallelism,   ONLY : fkbounds
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: itemp
    !! Counter on temperature
    !
    ! Local variables
    INTEGER :: iw
    !! Counter on frequency imag-axis
    INTEGER :: n
    !! Nr. of Matsubara frequencies
    INTEGER :: ik
    !! Counter on k-points
    INTEGER :: iq
    !! Counter on q-points for which k+sign*q is within the Fermi shell
    INTEGER :: iq0
    !! Index of iq on full q-mesh
    INTEGER :: lower_bnd, upper_bnd
    !! Lower/upper bound index after k paral
    INTEGER :: ibnd
    !! Counter on bandst k
    INTEGER :: jbnd
    !! Counter on bands at k+q
    INTEGER :: ierr
    !! Error status
    !
    REAL(KIND = DP) :: omega
    !! frequency imag-axis
    REAL(KIND = DP) :: lambda_eph
    !! electron-phonon coupling
    !
    ! SH: nsiw is replaced with "1+largest Matsubara index"
    !       for compatibility with the general case of sparse sampling
    ! RM: n = nsiw(itemp) for uniform samplin
    !
    n = wsn(nsiw(itemp)) + 1
    !
    CALL fkbounds(nkfs, lower_bnd, upper_bnd)
    !
    ALLOCATE(akeri(2 * n, nbndfs, MAXVAL(nqfs(:)), nbndfs, lower_bnd:upper_bnd), STAT = ierr)
    IF (ierr /= 0) CALL errore('kernel_aniso_iaxis', 'Error allocating akeri', 1)
    akeri(:, :, :, :, :) = zero
    !
    ! RM - if lambdar_aniso_ver2 is used then one needs to CALL evaluate_a2fij
    !
    DO ik = lower_bnd, upper_bnd
      DO ibnd = 1, nbndfs
        IF (ABS(ekfs(ibnd, ik) - ef0) < fsthick) THEN
          DO iq = 1, nqfs(ik)
            ! iq0 - index of q-point on the full q-mesh
            iq0 = ixqfs(ik, iq)
            DO jbnd = 1, nbndfs
              IF (ABS(ekfs(jbnd, ixkqf(ik, iq0)) - ef0) < fsthick) THEN
                DO iw = 1, 2 * n
                  omega = DBLE(2 * (iw - 1)) * pi * gtemp(itemp)
                  CALL lambdar_aniso_ver1(ik, iq, ibnd, jbnd, omega, lambda_eph)
                  !CALL lambdar_aniso_ver2(ik, iq, ibnd, jbnd, omega, lambda_eph)
                  akeri(iw, jbnd, iq, ibnd, ik) = lambda_eph
                ENDDO ! iw
              ENDIF ! fsthick
            ENDDO ! jbnd
          ENDDO ! iq
        ENDIF ! fsthick
      ENDDO ! ibnd
    ENDDO ! ik
    !
    RETURN
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE kernel_aniso_iaxis
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    SUBROUTINE lambdar_aniso_ver1(ik, iq, ibnd, jbnd, omega, lambda_eph)
    !-----------------------------------------------------------------------
    !
    ! computes lambda(ik, iq, ibnd, jbnd; n-n')
    ! reference H. Choi et. al, Physica C 385, 66 (2003)
    !
    USE kinds,         ONLY : DP
    USE modes,         ONLY : nmodes
    USE global_var,    ONLY : wf
    USE input,         ONLY : eps_acoustic
    USE supercond_common,     ONLY : ixqfs, g2, dosef
    USE ep_constants,  ONLY : zero
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: ik
    !! Counter on k-points
    INTEGER, INTENT(in) :: iq
    !! Counter on q-points for which k+sign*q is within the Fermi shell
    INTEGER, INTENT(in) :: ibnd
    !! Counter on bands at k
    INTEGER, INTENT(in) :: jbnd
    !! Counter on bands k+q
    !
    REAL(KIND = DP), INTENT(in) :: omega
    !! frequency on imag-axis
    REAL(KIND = DP), INTENT(out) :: lambda_eph
    !! electron-phonon coupling lambda_ij(k, k+q; n-n')
    !
    ! Local variables
    INTEGER :: iq0
    !! Index of iq on full q-mesh
    INTEGER :: imode
    !! Counter on phonon modes
    !
    ! iq0 - index of q-point on the full q-mesh
    iq0 = ixqfs(ik, iq)
    lambda_eph = zero
    DO imode = 1, nmodes  ! loop over frequency modes
      IF (wf(imode, iq0) > eps_acoustic) THEN
        lambda_eph = lambda_eph + g2(ik, iq, ibnd, jbnd, imode) * wf(imode, iq0) &
                   / (wf(imode, iq0) * wf(imode, iq0) + omega * omega)
      ENDIF
    ENDDO
    lambda_eph = 2.d0 * lambda_eph * dosef
    !
    RETURN
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE lambdar_aniso_ver1
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    SUBROUTINE lambdar_aniso_ver2(ik, iq, ibnd, jbnd, omega, lambda_eph)
    !-----------------------------------------------------------------------
    !
    ! computes lambda(ik, iq, ibnd, jbnd; n-n')
    ! reference H. Choi et. al, Physica C 385, 66 (2003)
    !
    USE kinds,         ONLY : DP
    USE input,         ONLY : nqstep
    USE supercond_common,     ONLY : a2fij, wsph, dwsph
    USE ep_constants,  ONLY : zero
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: ik
    !! Counter on k-points
    INTEGER, INTENT(in) :: iq
    !! Counter on q-points for which k+sign*q is within the Fermi shell
    INTEGER, INTENT(in) :: ibnd
    !! Counter on bands at k
    INTEGER, INTENT(in) :: jbnd
    !! Counter on bands at k+q
    REAL(KIND = DP), INTENT(in) :: omega
    !! frequency on imag-axis
    REAL(KIND = DP), INTENT(out) :: lambda_eph
    !! electron-phonon coupling lambda_ij(k, k+q; n-n')
    !
    ! Local variables
    INTEGER :: iwph
    !! Counter on frequency
    !
    ! Eq.(18) in Margine and Giustino, PRB 87, 024505 (2013)
    lambda_eph = zero
    DO iwph = 1, nqstep
      lambda_eph = lambda_eph + wsph(iwph) * a2fij(iwph, jbnd, iq, ibnd, ik) &
                 / (wsph(iwph) * wsph(iwph) + omega * omega)
    ENDDO
    lambda_eph = 2.d0 * lambda_eph * dwsph
    !
    RETURN
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE lambdar_aniso_ver2
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    SUBROUTINE kernel_aniso_analytic_cont(itemp)
    !-----------------------------------------------------------------------
    !!
    !! computes kernels K_{+}(w, iw_n, T) and K_{-}(w, iw_n, T)
    !! reference F. Masiglio, M. Schossmann, and J. Carbotte, PRB 37, 4965 (1988)
    !!
    USE kinds,         ONLY : DP
    USE global_var,    ONLY : wqf
    USE input,         ONLY : muc, fsthick
    USE supercond_common,     ONLY : nsw, nsiw, ws, wsi, adeltai, nkfs, nbndfs, &
                              spin_fac, dosef, ixkqf, ixqfs, nqfs, &
                              w0g, ekfs, ef0, adsumi, azsumi
    USE ep_constants,  ONLY : zero, one
    USE mp,            ONLY : mp_barrier, mp_sum
    USE mp_global,     ONLY : inter_pool_comm
    USE parallelism,   ONLY : fkbounds
    USE low_lvl,       ONLY : mem_size_eliashberg
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: itemp
    !! Counter on temperature index
    !
    ! Local variables
    INTEGER :: iw
    !! Counter on frequency real-axis
    INTEGER :: iwp
    !! Counter on frequency imag-axis
    INTEGER :: ik
    !! Counter on k-points
    INTEGER :: iq
    !! Counter on q-points for which k+sign*q is within the Fermi shell
    INTEGER :: iq0
    !! Index of iq on full q-mesh
    INTEGER :: lower_bnd, upper_bnd
    !! Lower/upper bound index after k paral
    INTEGER :: ibnd
    !! Counter on bands
    INTEGER :: jbnd
    !! Counter on bands
    INTEGER(8) :: imelt
    !! Counter memory
    INTEGER :: ierr
    !! Error status
    !
    REAL(KIND = DP) :: esqrt
    !! 1 / sqrt{w^2+\delta^2}
    REAL(KIND = DP) :: zesqrt
    !! w / sqrt{w^2+\delta^2}
    REAL(KIND = DP) :: desqrt
    !! \delta / sqrt{w^2+\delta^2}
    REAL(KIND = DP) :: kernelr
    !! 2 * Re[lambda(w - iw_n)]
    REAL(KIND = DP) :: kerneli
    !! 2 * Im[lambda(w - iw_n)]
    REAL(KIND = DP) :: weight
    !! factor in supercond. gap equations
    REAL(KIND = DP) :: inv_dos
    !! Invese dos inv_dos = 1/dosef. Defined for efficiency reason
    REAL(KIND = DP), ALLOCATABLE :: adeltai_tmp(:, :, :)
    !! Temporary array to collect adeltai from all pools
    !
    COMPLEX(KIND = DP) :: lambda_eph
    !! electron-phonon coupling lambda_ij(k, k+q; w-iw_n)
    !
    CALL fkbounds(nkfs, lower_bnd, upper_bnd)
    !
    ! get memory size required for adsumi, azsumi
    imelt = 2 * (upper_bnd - lower_bnd + 1) * nbndfs * nsw
    CALL mem_size_eliashberg(2, imelt)
    !
    ALLOCATE(adeltai_tmp(nsiw(itemp), nbndfs, nkfs), STAT = ierr)
    IF (ierr /= 0) CALL errore('kernel_aniso_analytic_cont', 'Error allocating adeltai_tmp', 1)
    adeltai_tmp(:, :, :)  = zero
    adeltai_tmp(:, :, lower_bnd:upper_bnd) = adeltai(:, :, lower_bnd:upper_bnd)
    ! collect contributions from all pools
    CALL mp_sum(adeltai_tmp, inter_pool_comm)
    CALL mp_barrier(inter_pool_comm)
    !
    ALLOCATE(adsumi(nsw, nbndfs, lower_bnd:upper_bnd), STAT = ierr)
    IF (ierr /= 0) CALL errore('kernel_aniso_analytic_cont', 'Error allocating adsumi', 1)
    ALLOCATE(azsumi(nsw, nbndfs, lower_bnd:upper_bnd), STAT = ierr)
    IF (ierr /= 0) CALL errore('kernel_aniso_analytic_cont', 'Error allocating azsumi', 1)
    adsumi(:, :, :) = zero
    azsumi(:, :, :) = zero
    !
    inv_dos = one / dosef
    !
    ! RM - if lambdai_aniso_ver2 is used then one needs to CALL evaluate_a2fij
    !
    DO ik = lower_bnd, upper_bnd
      DO ibnd = 1, nbndfs
        IF (ABS(ekfs(ibnd, ik) - ef0) < fsthick) THEN
          DO iq = 1, nqfs(ik)
            ! iq0 - index of q-point on the full q-mesh
            iq0 = ixqfs(ik, iq)
            DO jbnd = 1, nbndfs
              IF (ABS(ekfs(jbnd, ixkqf(ik, iq0)) - ef0) < fsthick) THEN
                weight = wqf(iq) * w0g(jbnd, ixkqf(ik, iq0)) * inv_dos
                DO iwp = 1, nsiw(itemp) ! loop over iw_n
                  esqrt = weight / DSQRT(wsi(iwp) * wsi(iwp) + &
!                          adeltai(iwp, jbnd, ixkqf(ik, iq0)) * adeltai(iwp, jbnd, ixkqf(ik, iq0)))
                          adeltai_tmp(iwp, jbnd, ixkqf(ik, iq0)) * adeltai_tmp(iwp, jbnd, ixkqf(ik, iq0)))
                  zesqrt =  esqrt * wsi(iwp)
!                  desqrt =  esqrt * adeltai(iwp, jbnd, ixkqf(ik, iq0))
                  desqrt =  esqrt * adeltai_tmp(iwp, jbnd, ixkqf(ik, iq0))
                  DO iw = 1, nsw ! loop over omega
                    CALL lambdai_aniso_ver1(ik, iq, ibnd, jbnd, ws(iw), wsi(iwp), lambda_eph)
                    !CALL lambdai_aniso_ver2(ik, iq, ibnd, jbnd, ws(iw), wsi(iwp), lambda_eph)
                    kernelr = 2.d0 * REAL(lambda_eph)
                    kerneli = 2.d0 * AIMAG(lambda_eph)
                    azsumi(iw, ibnd, ik) = azsumi(iw, ibnd, ik) + zesqrt * kerneli
                    adsumi(iw, ibnd, ik) = adsumi(iw, ibnd, ik) + desqrt * (kernelr - 2.d0 * muc * spin_fac)
                  ENDDO ! iw
                ENDDO ! iwp
              ENDIF
            ENDDO ! jbnd
          ENDDO ! iq
        ENDIF
      ENDDO ! ibnd
    ENDDO ! ik
    !
    DEALLOCATE(adeltai_tmp, STAT = ierr)
    IF (ierr /= 0) CALL errore('kernel_aniso_analytic_cont', 'Error deallocating adeltai_tmp', 1)
    !
    RETURN
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE kernel_aniso_analytic_cont
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    SUBROUTINE lambdai_aniso_ver1(ik, iq, ibnd, jbnd, omega, omegap, lambda_eph)
    !-----------------------------------------------------------------------
    !!
    !! computes lambda_ij(k, k+q; w-iw_n)
    !! reference F. Masiglio, M. Schossmann, and J. Carbotte, PRB 37, 4965 (1988)
    !!
    USE kinds,         ONLY : DP
    USE modes,         ONLY : nmodes
    USE global_var,    ONLY : wf
    USE input,         ONLY : eps_acoustic
    USE supercond_common,     ONLY : ixqfs, g2, dosef
    USE ep_constants,  ONLY : ci, czero
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: ik
    !! Counter on k-points
    INTEGER, INTENT(in) :: iq
    !! Counter on q-points for which k+sign*q is within the Fermi shell
    INTEGER, INTENT(in) :: ibnd
    !! Counter on bands at k
    INTEGER, INTENT(in) :: jbnd
    !! Counter on bands at k+q
    REAL(KIND = DP), INTENT(in) :: omega
    !! frequency w at point iw on real-axis
    REAL(KIND = DP), INTENT(in) :: omegap
    !! frequency w_n at point iwp on imag-axis
    COMPLEX(KIND = DP), INTENT(out) :: lambda_eph
    !! electron-phonon coupling lambda_ij(k, k+q; w-iw_n)
    !
    ! Local variables
    INTEGER :: iq0
    !! Index of iq on full q-mesh
    INTEGER :: imode
    !! Counter on phonon modes
    !
    ! iq0 - index of q-point on the full q-mesh
    iq0 = ixqfs(ik, iq)
    lambda_eph = czero
    DO imode = 1, nmodes  ! loop over frequency modes
      IF (wf(imode, iq0) > eps_acoustic) THEN
        lambda_eph = lambda_eph +  g2(ik, iq, ibnd, jbnd, imode) * wf(imode, iq0) &
                   / (wf(imode, iq0) * wf(imode, iq0) - (omega - ci * omegap) * (omega - ci * omegap))
      ENDIF
    ENDDO ! iwph
    lambda_eph = 2.d0 * lambda_eph * dosef
    !
    RETURN
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE lambdai_aniso_ver1
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    SUBROUTINE lambdai_aniso_ver2(ik, iq, ibnd, jbnd, omega, omegap, lambda_eph)
    !-----------------------------------------------------------------------
    !!
    !! computes lambda(w-iw_n)
    !! reference F. Masiglio, M. Schossmann, and J. Carbotte, PRB 37, 4965 (1988)
    !!
    USE kinds,         ONLY : DP
    USE input,         ONLY : nqstep
    USE supercond_common,     ONLY : a2fij, dwsph, wsph
    USE ep_constants,  ONLY : ci, czero
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: ik
    !! Counter on k-points
    INTEGER, INTENT(in) :: iq
    !! Counter on q-points for which k+sign*q is within the Fermi shell
    INTEGER, INTENT(in) :: ibnd
    !! Counter on bands at k
    INTEGER, INTENT(in) :: jbnd
    !! Counter on bands at k+q
    REAL(KIND = DP), INTENT(in) :: omega
    !! frequency w at point iw on real-axis
    REAL(KIND = DP), INTENT(in) :: omegap
    !! frequency w_n at point iwp on imag-axis
    COMPLEX(KIND = DP), INTENT(out) :: lambda_eph
    !! electron-phonon coupling lambda_ij(k, k+q; w-iw_n)
    !
    ! Local variables
    INTEGER :: iwph
    !! Counter on frequency
    !
    lambda_eph = czero
    DO iwph = 1, nqstep
      lambda_eph = lambda_eph + wsph(iwph) * a2fij(iwph, jbnd, iq, ibnd, ik) &
                 / (wsph(iwph) * wsph(iwph) - (omega - ci * omegap) * (omega - ci * omegap))
    ENDDO ! iwph
    lambda_eph = 2.d0 * lambda_eph * dwsph
    !
    RETURN
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE lambdai_aniso_ver2
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    SUBROUTINE evaluate_a2fij
    !-----------------------------------------------------------------------
    !!
    !! computes the anisotropic spectral function a2F(k, k', w)
    !!
    USE kinds,         ONLY : DP
    USE modes,         ONLY : nmodes
    USE global_var,    ONLY : wf
    USE input,         ONLY : fsthick, eps_acoustic, nqstep, degaussq
    USE supercond_common,     ONLY : nkfs, nbndfs, g2, a2fij, ixkqf, ixqfs, nqfs, ekfs, ef0, &
                              dosef, wsph
    USE ep_constants,  ONLY : zero, one
    USE parallelism,   ONLY : fkbounds
    !
    IMPLICIT NONE
    !
    ! Local variables
    INTEGER :: ik
    !! Counter on k-points
    INTEGER :: iq
    !! Counter on q-points for which k+sign*q is within the Fermi shell
    INTEGER :: ibnd
    !! Counter on bands at k
    INTEGER :: jbnd
    !! Counter on bands at k+q
    INTEGER :: iq0
    !! Index of iq on full q-mesh
    INTEGER :: lower_bnd, upper_bnd
    !! Lower/upper bound index after k paral
    INTEGER :: imode
    !! Counter on phonon modes
    INTEGER :: iwph
    !! Counter on frequency
    REAL(KIND = DP) :: inv_degaussq
    !! 1.0/degaussq. Defined for efficiency reasons
    INTEGER :: ierr
    !! Error status
    !
    REAL(KIND = DP) :: weight
    !! Factor in a2fij
    REAL(KIND = DP), EXTERNAL :: w0gauss
    !! The derivative of wgauss: an approximation to the delta function
    !
    CALL fkbounds(nkfs, lower_bnd, upper_bnd)
    ALLOCATE(a2fij(nqstep, nbndfs, MAXVAL(nqfs(:)), nbndfs, lower_bnd:upper_bnd), STAT = ierr)
    IF (ierr /= 0) CALL errore('evaluate_a2fij', 'Error allocating a2fij', 1)
    a2fij(:, :, :, :, :) = zero
    !
    inv_degaussq = one / degaussq
    !
    DO ik = lower_bnd, upper_bnd
      DO ibnd = 1, nbndfs
        IF (ABS(ekfs(ibnd, ik) - ef0) < fsthick) THEN
          DO iq = 1, nqfs(ik)
            ! iq0 - index of q-point on the full q-mesh
            iq0 = ixqfs(ik, iq)
            DO jbnd = 1, nbndfs
              IF (ABS(ekfs(jbnd, ixkqf(ik, iq0)) - ef0) < fsthick) THEN
                DO imode = 1, nmodes
                  IF (wf(imode, iq0) > eps_acoustic) THEN
                    DO iwph = 1, nqstep
                      weight  = w0gauss((wsph(iwph) - wf(imode, iq0)) * inv_degaussq, 0) * inv_degaussq
                      a2fij(iwph, jbnd, iq, ibnd, ik) = a2fij(iwph, jbnd, iq, ibnd, ik) &
                                             + weight * dosef * g2(ik, iq, ibnd, jbnd, imode)
                    ENDDO ! iwph
                  ENDIF ! wf
                ENDDO ! imode
              ENDIF
            ENDDO ! jbnd
          ENDDO ! iq
        ENDIF
      ENDDO ! ibnd
    ENDDO ! ik
    !
    RETURN
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE evaluate_a2fij
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    SUBROUTINE mu_inter_aniso_newton(itemp, muintr, nel, nstate)
    !-----------------------------------------------------------------------
    !!
    !! SH: To find the superconducting state chemical potential
    !!       using the Newton-Raphson method (Nov 2021).
    !!
    USE kinds,          ONLY : DP
    USE supercond_common,      ONLY : ekfs, ef0, nkfs, nbndfs, wsi, nsiw, &
                               adeltaip, aznormip, ashiftip, wkfs
    USE global_var,     ONLY : gtemp
    USE input,          ONLY : fsthick, nsiter, broyden_beta, positive_matsu
    USE ep_constants,   ONLY : zero, eps6
    !
    IMPLICIT NONE
    !
    INTEGER         :: itemp
    !! Temperature
    REAL(KIND = DP) :: muintr
    !! Interacting chemical potential: initial value
    REAL(KIND = DP) :: nel
    !! Without SOC: Nr. of electrons within Fermi window
    !! With SOC: 0.5 * Nr. of electrons within Fermi window
    REAL(KIND = DP) :: nstate
    !! Nr. of states within Fermi window
    !
    ! Local variables
    LOGICAL :: conv
    !! Convergence parameter
    INTEGER :: ik
    !! Index for momentum state
    INTEGER :: ibnd
    !! Index for band
    INTEGER :: iw
    !! Index for Matsubara frequencies
    INTEGER :: iter
    !! Counter for iterations
    REAL(KIND = DP) :: delta
    !! Temporary variable to store energy difference
    REAL(KIND = DP) :: theta
    !! Theta in Eliashberg equations
    REAL(KIND = DP) :: muout, muin
    !! Variables for Newton's minimization
    REAL(KIND = DP) :: fmu
    !! To store f(mu) value
    REAL(KIND = DP) :: dmu
    !! To store d_f(mu)/d_mu value
    !
    muin = muintr
    !
    iter = 1
    conv = .FALSE.
    !
    DO WHILE (.NOT. conv .AND. iter <= nsiter)
      !
      fmu = zero ! f(mu)
      dmu = zero ! d_f(mu)/d_mu
      !
      ! SH: Here, we have an explicit summation over full range of frequencies
      !      for even functions, and multiply that by two (loop over iw).
      !
      DO ik = 1, nkfs
        DO ibnd = 1, nbndfs
          IF (ABS(ekfs(ibnd, ik) - ef0) < fsthick) THEN
            DO iw = 1, nsiw(itemp)
              delta = ekfs(ibnd, ik) - muin + ashiftip(iw, ibnd, ik)
              theta = (wsi(iw) * aznormip(iw, ibnd, ik))**2.d0 + &
                      (ekfs(ibnd, ik) - muin + ashiftip(iw, ibnd, ik))**2.d0 + &
                      (aznormip(iw, ibnd, ik) * adeltaip(iw, ibnd, ik))**2.d0
              IF (positive_matsu) THEN
                fmu = fmu + 2.d0 * wkfs(ik) * delta / theta
                dmu = dmu + 2.d0 * wkfs(ik) * (2.d0 * delta**2.d0 - theta) / (theta**2.d0)
              ELSE
                fmu = fmu + wkfs(ik) * delta / theta
                dmu = dmu + wkfs(ik) * (2.d0 * delta**2.d0 - theta) / (theta**2.d0)
              ENDIF
            ENDDO ! iw
          ENDIF
        ENDDO ! ibnd
      ENDDO ! ik
      !
      ! HP: factor 2 is already included above in wkfs(ik)
      fmu = nstate - gtemp(itemp) * fmu - nel
      dmu = - gtemp(itemp) * dmu
      !
      muout = muin - fmu / dmu
      !
      ! HP: linear mixing
      muout = (1.d0 - ABS(broyden_beta)) * muin + ABS(broyden_beta) * muout
      !
      IF (ABS((muout - muin) / muin) <= eps6) conv = .TRUE.
      !
      muin = muout
      iter = iter + 1
      !
    END DO
    !
    IF (.NOT. conv .AND. (iter - 1) == nsiter) &
      CALL errore('mu_inter_aniso_newton', 'Error failed to find the mu_inter_aniso_newton value',1)
    !
    muintr = muout
    !
    RETURN
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE mu_inter_aniso_newton
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    SUBROUTINE mu_inter_aniso(itemp, muintr, nel, nstate, ns, ir_obj)
    !-----------------------------------------------------------------------
    !! HM: To find the superconducting state chemical potential
    !!     using Brent's method, which is similar to the bisection
    !!     method but more efficient for converging fastly.
    !
    USE kinds,         ONLY : DP
    USE input,         ONLY : fsthick, nsiter, fbw, gridsamp
    USE ep_constants,  ONLY : eps4, eps6, zero, one
    USE sparse_ir,     ONLY : IR
    !
    IMPLICIT NONE
    !
    INTEGER         :: itemp
    !! Temperature
    REAL(KIND = DP) :: muintr
    !! Interacting chemical potential: initial value
    REAL(KIND = DP) :: nel
    !! Without SOC: Nr. of electrons within Fermi window
    !! With SOC: 0.5 * Nr. of electrons within Fermi window
    REAL(KIND = DP) :: nstate
    !! Nr. of states within Fermi window
    INTEGER, INTENT(IN), OPTIONAL :: ns
    !! Total number states within the fsthick window
    TYPE(IR), INTENT(IN), OPTIONAL :: ir_obj
    !! which contains ir objects such as basis functions
    !
    ! Local variables
    LOGICAL :: conv
    !! Convergence parameter
    INTEGER :: iter
    !! Counter for iterations
    REAL(KIND = DP) :: a
    !! point determined so that f(a) has a different sign from f(b).
    REAL(KIND = DP) :: b
    !! current iterate
    !! |f(b)| is always smaller than |f(a)|; b is assumed to be closer to solution than a.
    REAL(KIND = DP) :: c
    !! previous iterate
    REAL(KIND = DP) :: xm
    !! 0.5_DP * (c - b)
    REAL(KIND = DP) :: d, e
    !! dummys
    REAL(KIND = DP) :: fmu_a
    !! To store f(mu = a) value
    REAL(KIND = DP) :: fmu_b
    !! To store f(mu = b) value
    REAL(KIND = DP) :: fmu_c
    !! To store f(mu = c) value
    REAL(KIND = DP) :: p, q, r, s
    !! dummys
    REAL(KIND = DP) :: tol
    !! eps6 * ABS(muintr)
    REAL(KIND = DP) :: tol1
    !! specific numerical tolerance
    REAL(KIND = DP), PARAMETER :: eps = EPSILON(one)
    !! a positive model number that is almost negligible compared to unity
    !
    IF (gridsamp == 2) THEN
      IF (.NOT. PRESENT(ns)) THEN
        CALL errore('mu_inter_aniso', 'Error: ns is not given while gridsamp = 2', 1)
      ENDIF
      IF (.NOT. PRESENT(ir_obj)) THEN
        CALL errore('mu_inter_aniso', 'Error: ir_obj is not given while gridsamp = 2', 1)
      ENDIF
    ENDIF
    !
    tol = eps6 * ABS(muintr)
    tol1 = one + eps
    iter = 1
    conv = .FALSE.
    !
    a = muintr - 1.D-1 * fsthick
    b = muintr + 1.D-1 * fsthick
    !c = muintr
    IF (gridsamp .NE. 2) THEN
      CALL calc_fmu_seq(itemp, nel, nstate, a, fmu_a)
      CALL calc_fmu_seq(itemp, nel, nstate, b, fmu_b)
      !CALL calc_fmu_seq(itemp, nel, nstate, c, fmu_c)
    ELSE
      CALL calc_fmu_seq_ir(itemp, ns, nel, nstate, a, fmu_a, ir_obj)
      CALL calc_fmu_seq_ir(itemp, ns, nel, nstate, b, fmu_b, ir_obj)
      !CALL calc_fmu_seq_ir(itemp, ns, nel, nstate, c, fmu_c, ir_obj)
    ENDIF
    !
    DO WHILE ((iter .LE. nsiter) .AND. (fmu_a*fmu_b .GE. zero))
      IF (fmu_a .GE. zero) THEN
        a = a - 1.D-1 * fsthick
        IF (gridsamp .NE. 2) THEN
          CALL calc_fmu_seq(itemp, nel, nstate, a, fmu_a)
        ELSE
          CALL calc_fmu_seq_ir(itemp, ns, nel, nstate, a, fmu_a, ir_obj)
        ENDIF
      ELSEIF (fmu_b .LE. zero) THEN
        b = b + 1.D-1 * fsthick
        IF (gridsamp .NE. 2) THEN
          CALL calc_fmu_seq(itemp, nel, nstate, b, fmu_b)
        ELSE
          CALL calc_fmu_seq_ir(itemp, ns, nel, nstate, b, fmu_b, ir_obj)
        ENDIF
      ENDIF
      iter = iter + 1
    ENDDO
    !
    IF (fmu_a*fmu_b > zero) THEN
      CALL errore('mu_inter_aniso', 'Error: initial guess is quite far from the solution, &
                  &or wscut is too small.',1)
    ENDIF
    !
    IF (fmu_a == zero) THEN
      muintr = a
      conv = .TRUE.
      RETURN
    ELSEIF (fmu_b == zero) THEN
      muintr = b
      conv = .TRUE.
      RETURN
    ENDIF
    !
    iter = 1
    !
    DO WHILE (.NOT. conv .AND. iter .LE. nsiter)
      !
      IF ((iter == 1) .OR. &
          ((fmu_b * (fmu_c / ABS(fmu_c))) > zero)) THEN
        c = a
        fmu_c = fmu_a
        d = b - a
        e = d
      ENDIF
      !
      IF (ABS(fmu_c) .LT. ABS(fmu_b)) THEN
        a = b
        b = c
        c = a
        fmu_a = fmu_b
        fmu_b = fmu_c
        fmu_c = fmu_a
      ENDIF
      !
      tol1 = 2.0_DP * eps * ABS(b) + 0.5_DP * tol
      xm = 0.5_DP * (c - b)
      IF ((ABS(xm) .LE. tol1) .OR. (fmu_b == zero)) THEN
        muintr = b
        conv = .TRUE.
        RETURN
      ENDIF
      !
      ! see if a bisection is forced
      IF ((ABS(e) .GE. tol1) .AND. (ABS(fmu_a) .GT. ABS(fmu_b))) THEN
        s = fmu_b / fmu_a
        IF (a .NE. c) THEN
          ! inverse quadratic interpolation
          q = fmu_a / fmu_c
          r = fmu_b / fmu_c
          p = s * (2.0_DP * xm * q * (q - r) - (b - a) * (r - one))
          q = (q - one) * (r - one) * (s - one)
        ELSE
          ! linear interpolation
          p = 2.0_DP * xm * s
          q = one - s
        ENDIF
        IF (p .LE. zero) THEN
          p = -p
        ELSE
          q = -q
        ENDIF
        s = e
        e = d
        IF (((2.0_DP * p) .GE. (3.0_DP * xm * q - ABS(tol1 * q))) .OR. &
            (p .GE. ABS(0.5_DP * s * q))) THEN
          d = xm
          e = d
        ELSE
          d = p / q
        ENDIF
      ELSE
        d = xm
        e = d
      ENDIF
      !
      a = b
      fmu_a = fmu_b
      IF (ABS(d) .LE. tol1) THEN
        IF (xm .LE. zero) THEN
          b = b - tol1
        ELSE
          b = b + tol1
        ENDIF
      ELSE
        b = b + d
      ENDIF
      IF (gridsamp .NE. 2) THEN
        CALL calc_fmu_seq(itemp, nel, nstate, b, fmu_b)
      ELSE
        CALL calc_fmu_seq_ir(itemp, ns, nel, nstate, b, fmu_b, ir_obj)
      ENDIF
      !
      iter = iter + 1
      !
    ENDDO
    !
    IF (.NOT. conv .AND. (iter - 1) == nsiter) &
      CALL errore('mu_inter_aniso', 'Error failed to find the mu_inter_aniso value',1)
    !
    muintr = b
    !
    RETURN
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE mu_inter_aniso
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    SUBROUTINE calc_fmu_seq(itemp, nel, nstate, mu, fmu)
    !-----------------------------------------------------------------------
    !! This subroutine calculates fmu
    !!
    !!    fmu = nstate + SUM_{ik, ibnd}( wkfs * ( ekfs - mu + ashiftip ) ) - nel
    !!
    !!  where fmu is the difference between the target nel and the calculated nel
    !
    USE kinds,          ONLY : DP
    USE supercond_common,      ONLY : ekfs, ef0, nkfs, nbndfs, wsi, nsiw, &
                               adeltaip, aznormip, ashiftip, wkfs
    USE global_var,     ONLY : gtemp
    USE input,          ONLY : fsthick, broyden_beta, positive_matsu
    USE ep_constants,   ONLY : zero, one
    !
    IMPLICIT NONE
    !
    INTEGER         :: itemp
    !! Temperature
    REAL(KIND = DP) :: nel
    !! Without SOC: Nr. of electrons within Fermi window
    !! With SOC: 0.5 * Nr. of electrons within Fermi window
    REAL(KIND = DP) :: nstate
    !! Nr. of states within Fermi window
    REAL(KIND = DP) :: mu
    !! trial chemical potential
    REAL(KIND = DP) :: fmu
    !! fmu, which is the difference between the target nel and the calculated nel
    !
    ! Local variables
    INTEGER :: ik
    !! Index for momentum state
    INTEGER :: ibnd
    !! Index for band
    INTEGER :: iw
    !! Index for Matsubara frequencies
    INTEGER :: ierr
    !! Error status
    REAL(KIND = DP) :: delta
    !! Temporary variable to store energy difference
    REAL(KIND = DP) :: inv_theta
    !! inverse Theta in Eliashberg equations
    !
    fmu = zero
    !
    DO ik = 1, nkfs
      DO ibnd = 1, nbndfs
        IF (ABS(ekfs(ibnd, ik) - ef0) < fsthick) THEN
          DO iw = 1, nsiw(itemp)
            delta = ekfs(ibnd, ik) - mu + ashiftip(iw, ibnd, ik)
            inv_theta = one / ((wsi(iw) * aznormip(iw, ibnd, ik))**2.d0 + &
                    (ekfs(ibnd, ik) - mu + ashiftip(iw, ibnd, ik))**2.d0 + &
                    (aznormip(iw, ibnd, ik) * adeltaip(iw, ibnd, ik))**2.d0)
            IF (positive_matsu) THEN
              fmu = fmu + 2.d0 * wkfs(ik) * delta * inv_theta
            ELSE
              fmu = fmu + wkfs(ik) * delta * inv_theta
            ENDIF
          ENDDO ! iw
        ENDIF
      ENDDO ! ibnd
    ENDDO ! ik
    !
    ! HP: factor 2 is already included above in wkfs(ik)
    fmu = nstate - gtemp(itemp) * fmu - nel
    !
    RETURN
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE calc_fmu_seq
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    SUBROUTINE calc_fmu_seq_ir(itemp, ns, nel, nstate, mu, fmu, ir_obj)
    !-----------------------------------------------------------------------
    !! This subroutine calculates fmu
    !!
    !!    F(mu) = SUM_{ik, ibnd, ispin}(1 + G(tau = 0+)) - Nel
    !!
    !! but in practical,
    !!
    !!    fmu = 2.0D0 * nstate + SUM_{ik, ibnd}( wkfs * G(tau = 0+) ) - nel
    !!
    !!  where fmu is the difference between the target nel and the calculated nel
    !
    USE kinds,          ONLY : DP
    USE supercond_common,  ONLY : ekfs, ef0, nkfs, nbndfs, wsi, nsiw, &
                               adeltaip, aznormip, ashiftip, wkfs
    USE input,          ONLY : fsthick, positive_matsu
    USE ep_constants,   ONLY : one, zero, eps6, ci, cone, czero
    USE sparse_ir,      ONLY : IR, fit_matsubara_f
    !
    IMPLICIT NONE
    !
    INTEGER         :: itemp
    !! Temperature
    INTEGER         :: ns
    !! Total number states within the fsthick window
    REAL(KIND = DP) :: nel
    !! Without SOC: Nr. of electrons within Fermi window
    !! With SOC: 0.5 * Nr. of electrons within Fermi window
    REAL(KIND = DP) :: nstate
    !! Nr. of states within Fermi window
    REAL(KIND = DP) :: mu
    !! trial chemical potential
    REAL(KIND = DP) :: fmu
    !! fmu, which is the difference between the target nel and the calculated nel
    TYPE(IR), INTENT(IN) :: ir_obj
    !! which contains ir objects such as basis functions
    !
    ! Local variables
    INTEGER :: n
    !! Counter for states within the fsthick window
    INTEGER :: ik
    !! Index for momentum state
    INTEGER :: ibnd
    !! Index for band
    INTEGER :: iw
    !! Index for Matsubara frequencies
    INTEGER :: ierr
    !! Error status
    REAL(KIND = DP) :: inv_theta
    !! inverse Theta in Eliashberg equations
    REAL(KIND = DP), ALLOCATABLE :: gl_d(:, :)
    !! Green's functions in IR
    REAL(KIND = DP), ALLOCATABLE :: u0p_d(:)
    !! IR basis functions
    REAL(KIND = DP), ALLOCATABLE :: gtau0p_d(:)
    !! Green's functions of imaginary time at tau = 0+
    COMPLEX(KIND = DP), ALLOCATABLE :: gl(:, :)
    !! Green's functions in IR
    COMPLEX(KIND = DP), ALLOCATABLE :: u0p(:)
    !! IR basis functions
    COMPLEX(KIND = DP), ALLOCATABLE :: gtau0p(:)
    !! Green's functions of imaginary time at tau = 0+
    COMPLEX(KIND = DP) :: giw(ns, ir_obj%nfreq_f)
    !! Green's functions of Matsubara freq.: G = -(iw*Z + chi) / theta
    !COMPLEX(KIND = DP) :: gtau(ns, ir_obj%ntau)
    !!
    !
    fmu = zero
    !
    n = 0
    DO ik = 1, nkfs
      DO ibnd = 1, nbndfs
        IF (ABS(ekfs(ibnd, ik) - ef0) < fsthick) THEN
          n = n + 1
          DO iw = 1, ir_obj%nfreq_f
            !
            inv_theta = one / ((wsi(iw) * aznormip(iw, ibnd, ik))**2.d0 + &
                    (ekfs(ibnd, ik) - mu + ashiftip(iw, ibnd, ik))**2.d0 + &
                    (aznormip(iw, ibnd, ik) * adeltaip(iw, ibnd, ik))**2.d0)
            giw(n, iw) = ci * wsi(iw) * aznormip(iw, ibnd, ik) + &
                         cone * (ekfs(ibnd, ik) - mu + ashiftip(iw, ibnd, ik))
            giw(n, iw) = - giw(n, iw) * inv_theta
            !
          ENDDO
        ENDIF
      ENDDO
    ENDDO
    !
    IF (positive_matsu) THEN
      !
      ALLOCATE(gl_d(ns, ir_obj%size), STAT = ierr)
      IF (ierr /= 0) CALL errore('calc_fmu_seq', 'Error allocating gl', 1)
      ALLOCATE(u0p_d(ir_obj%size), STAT = ierr)
      IF (ierr /= 0) CALL errore('calc_fmu_seq', 'Error allocating u0p', 1)
      ALLOCATE(gtau0p_d(ns), STAT = ierr)
      IF (ierr /= 0) CALL errore('calc_fmu_seq', 'Error allocating gtau0p', 1)
      !
      ! To obtain expansion coefficients of Green's function
      ! and kernel in the IR basis
      CALL fit_matsubara_f(ir_obj, giw, gl_d)
      !
      !CALL evaluate_tau(ir_obj, gl_d, gtau_d)
      !! execute
      !!   gtau = MATMUL(gl, transpose(ir_obj%u%a_real))
      !! in evaluate_tau
      !
      !gtau0p_d(:) = gtau_d(:, 1)
      !
      ! To obtain Green's function of imaginary time
      ! from the corresponding expansion coefficients
      ! Here, we want it only for tau = 0+,
      ! so using DGEMV instead of evaluate_tau.
      u0p_d(:) = ir_obj%u%a_real(1, :)
      CALL DGEMV ('n', ns, ir_obj%size, one, gl_d, ns, u0p_d, 1, zero, gtau0p_d, 1)
      !
      n = 0
      DO ik = 1, nkfs
        DO ibnd = 1, nbndfs
          IF (ABS(ekfs(ibnd, ik) - ef0) < fsthick) THEN
            n = n + 1
            !
            fmu = fmu + wkfs(ik) * gtau0p_d(n)
          ENDIF
        ENDDO
      ENDDO
      !
      DEALLOCATE(gl_d, STAT = ierr)
      IF (ierr /= 0) CALL errore('calc_fmu_seq', 'Error deallocating gl', 1)
      DEALLOCATE(u0p_d, STAT = ierr)
      IF (ierr /= 0) CALL errore('calc_fmu_seq', 'Error deallocating u0p', 1)
      DEALLOCATE(gtau0p_d, STAT = ierr)
      IF (ierr /= 0) CALL errore('calc_fmu_seq', 'Error deallocating gtau0p', 1)
      !
    ELSE
      !
      ALLOCATE(gl(ns, ir_obj%size), STAT = ierr)
      IF (ierr /= 0) CALL errore('calc_fmu_seq', 'Error allocating gl', 1)
      ALLOCATE(u0p(ir_obj%size), STAT = ierr)
      IF (ierr /= 0) CALL errore('calc_fmu_seq', 'Error allocating u0p', 1)
      ALLOCATE(gtau0p(ns), STAT = ierr)
      IF (ierr /= 0) CALL errore('calc_fmu_seq', 'Error allocating gtau0p', 1)
      !
      ! To obtain expansion coefficients of green's function
      ! and kernel in the IR basis
      CALL fit_matsubara_f (ir_obj, giw, gl)
      !
      !CALL evaluate_tau (ir_obj, gl, gtau)
      !! execute
      !!   gtau = MATMUL(gl, transpose(ir_obj%u%a))
      !! in evaluate_tau
      !
      !gtau0p(:) = gtau(:, 1)
      !
      ! To obtain green's function of imaginary time
      ! from the corresponding expansion coefficients
      ! Here, we want it only for tau = 0+,
      ! so using ZGEMV instead of evaluate_tau.
      u0p(:) = ir_obj%u%a(1, :)
      CALL ZGEMV ('n', ns, ir_obj%size, cone, gl, ns, u0p, 1, czero, gtau0p, 1)
      !
      n = 0
      DO ik = 1, nkfs
        DO ibnd = 1, nbndfs
          IF (ABS(ekfs(ibnd, ik) - ef0) < fsthick) THEN
            n = n + 1
            !
            fmu = fmu + wkfs(ik) * REAL(gtau0p(n), KIND=DP)
          ENDIF
        ENDDO
      ENDDO
      !
      DEALLOCATE(gl, STAT = ierr)
      IF (ierr /= 0) CALL errore('calc_fmu_seq', 'Error deallocating gl', 1)
      DEALLOCATE(u0p, STAT = ierr)
      IF (ierr /= 0) CALL errore('calc_fmu_seq', 'Error deallocating u0p', 1)
      DEALLOCATE(gtau0p, STAT = ierr)
      IF (ierr /= 0) CALL errore('calc_fmu_seq', 'Error deallocating gtau0p', 1)
      !
    ENDIF
    !
    fmu = 2.0D0 * nstate + fmu - nel
    !
    RETURN
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE calc_fmu_seq_ir
    !-----------------------------------------------------------------------
    !
    !----------------------------------------------------------------------
    SUBROUTINE deallocate_aniso_iaxis(iset)
    !----------------------------------------------------------------------
    !!
    !!  deallocates the variables allocated for imag-axis solutions
    !!
    !!  iset = 1 is for deallocating arrays except for the arrays used in Pade.
    !!  iset = 2 is for deallocating the remaining arrays that were not deallocated when iset = 1.
    !!           It should be used to avoid deallocating arrays not allocated when limag .AND. imag_read .AND. itemp == 1.
    !!  iset = 3 is for deallocating arrays that are deallocated when iset = 1 and 2 at once.
    !!  iset = 4 is for deallocating the remaining arrays that were not deallocated when iset = 1,
    !!           which include the arrays skipped if iset = 2.
    !!  iset = 5 is for deallocating arrays that are deallocated when iset = 1 and 4 at once.
    !!
    !!  Note that iset = 3 is not used currently.
    !
    USE control_flags, ONLY : iverbosity
    USE input,         ONLY : fbw, icoulomb, gridsamp, positive_matsu
    USE supercond_common,     ONLY : wsi, wsn, deltai, znormi, adeltai, adeltaip, &
                              aznormi, aznormip, shifti, ashifti, ashiftip, &
                              fft_in1, fft_out1, fft_in2, fft_out2, &
                              ir_giw, ir_gl, ir_gtau, ir_knliw, ir_knll, &
                              ir_knltau, ir_cvliw, ir_cvll, ir_cvltau, &
                              adeltai_cl, adeltaip_cl, w_stat, &
                              ir_giw_cl, ir_gl_cl, ir_gtau_cl, ir_gl_d, &
                              ir_gtau_d, ir_knll_d, ir_knltau_d, ir_cvll_d, &
                              ir_cvltau_d, ir_gl_cl_d, ir_gtau_cl_d, &
                              weight_q, weight_cl, num_js1, num_js2, num_js3, &
                              gl_abs, fl_abs, knll_abs
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: iset
    !! Define what set of variables are deallocated
    !
    INTEGER :: ierr
    !! Error status
    !
    IF (iset == 1 .OR. iset == 3 .OR. iset == 5) THEN
      ! sum_eliashberg_aniso_iaxis
      DEALLOCATE(deltai, STAT = ierr)
      IF (ierr /= 0) CALL errore('deallocate_aniso_iaxis', 'Error deallocating deltai', 1)
      DEALLOCATE(znormi, STAT = ierr)
      IF (ierr /= 0) CALL errore('deallocate_aniso_iaxis', 'Error deallocating znormi', 1)
      DEALLOCATE(adeltaip, STAT = ierr)
      IF (ierr /= 0) CALL errore('deallocate_aniso_iaxis', 'Error deallocating adeltaip', 1)
      !
      IF (fbw) THEN
        ! SH: deallocate fbw-related arrays in sum_eliashberg_aniso_iaxis
        DEALLOCATE(shifti, STAT = ierr)
        IF (ierr /= 0) CALL errore('deallocate_aniso_iaxis', 'Error deallocating shifti', 1)
        DEALLOCATE(aznormip, STAT = ierr)
        IF (ierr /= 0) CALL errore('deallocate_aniso_iaxis', 'Error deallocating aznormip', 1)
        DEALLOCATE(ashiftip, STAT = ierr)
        IF (ierr /= 0) CALL errore('deallocate_aniso_iaxis', 'Error deallocating ashiftip', 1)
      ENDIF
    ENDIF
    !
    IF (iset == 2 .OR. iset == 3 .OR. iset == 4 .OR. iset == 5) THEN
      ! gen_freqgrid_iaxis
      DEALLOCATE(wsi, STAT = ierr)
      IF (ierr /= 0) CALL errore('deallocate_aniso_iaxis', 'Error deallocating wsi', 1)
      DEALLOCATE(wsn, STAT = ierr)
      IF (ierr /= 0) CALL errore('deallocate_aniso_iaxis', 'Error deallocating wsn', 1)
      !
      ! sum_eliashberg_aniso_iaxis
      DEALLOCATE(adeltai, STAT = ierr)
      IF (ierr /= 0) CALL errore('deallocate_aniso_iaxis', 'Error deallocating adeltai', 1)
      DEALLOCATE(aznormi, STAT = ierr)
      IF (ierr /= 0) CALL errore('deallocate_aniso_iaxis', 'Error deallocating aznormi', 1)
      !
      IF (fbw) THEN
        ! SH: deallocate fbw-related arrays in sum_eliashberg_aniso_iaxis
        DEALLOCATE(ashifti, STAT = ierr)
        IF (ierr /= 0) CALL errore('deallocate_aniso_iaxis', 'Error deallocating ashifti', 1)
      ENDIF
    ENDIF
    !
    IF  (iset == 4 .OR. iset == 5) THEN
      ! If limag .AND. imag_read .AND. itemp == 1,
      ! the following arrays is not allocated where this subroutine is called.
      IF (gridsamp == 2) THEN
        ! HM: deallocate ir-related arrays in sum_eliashberg_aniso_iaxis
        DEALLOCATE(ir_giw, STAT = ierr)
        IF (ierr /= 0) CALL errore('deallocate_aniso_iaxis', 'Error deallocating ir_giw', 1)
        DEALLOCATE(ir_knliw, STAT = ierr)
        IF (ierr /= 0) CALL errore('deallocate_aniso_iaxis', 'Error deallocating ir_knliw', 1)
        DEALLOCATE(ir_cvliw, STAT = ierr)
        IF (ierr /= 0) CALL errore('deallocate_aniso_iaxis', 'Error deallocating ir_cvliw', 1)
        IF (positive_matsu) THEN
          DEALLOCATE(ir_gl_d, STAT = ierr)
          IF (ierr /= 0) CALL errore('deallocate_aniso_iaxis', 'Error deallocating ir_gl_d', 1)
          DEALLOCATE(ir_gtau_d, STAT = ierr)
          IF (ierr /= 0) CALL errore('deallocate_aniso_iaxis', 'Error deallocating ir_gtau_d', 1)
          DEALLOCATE(ir_knll_d, STAT = ierr)
          IF (ierr /= 0) CALL errore('deallocate_aniso_iaxis', 'Error deallocating ir_knll_d', 1)
          DEALLOCATE(ir_knltau_d, STAT = ierr)
          IF (ierr /= 0) CALL errore('deallocate_aniso_iaxis', 'Error deallocating ir_knltau_d', 1)
          DEALLOCATE(ir_cvll_d, STAT = ierr)
          IF (ierr /= 0) CALL errore('deallocate_aniso_iaxis', 'Error deallocating ir_cvll_d', 1)
          DEALLOCATE(ir_cvltau_d, STAT = ierr)
          IF (ierr /= 0) CALL errore('deallocate_aniso_iaxis', 'Error deallocating ir_cvltau_d', 1)
        ELSE
          DEALLOCATE(ir_gl, STAT = ierr)
          IF (ierr /= 0) CALL errore('deallocate_aniso_iaxis', 'Error deallocating ir_gl', 1)
          DEALLOCATE(ir_gtau, STAT = ierr)
          IF (ierr /= 0) CALL errore('deallocate_aniso_iaxis', 'Error deallocating ir_gtau', 1)
          DEALLOCATE(ir_knll, STAT = ierr)
          IF (ierr /= 0) CALL errore('deallocate_aniso_iaxis', 'Error deallocating ir_knll', 1)
          DEALLOCATE(ir_knltau, STAT = ierr)
          IF (ierr /= 0) CALL errore('deallocate_aniso_iaxis', 'Error deallocating ir_knltau', 1)
          DEALLOCATE(ir_cvll, STAT = ierr)
          IF (ierr /= 0) CALL errore('deallocate_aniso_iaxis', 'Error deallocating ir_cvll', 1)
          DEALLOCATE(ir_cvltau, STAT = ierr)
          IF (ierr /= 0) CALL errore('deallocate_aniso_iaxis', 'Error deallocating ir_cvltau', 1)
        ENDIF
        DEALLOCATE(weight_q, STAT = ierr)
        IF (ierr /= 0) CALL errore('deallocate_aniso_iaxis', 'Error deallocating weight_q', 1)
        DEALLOCATE(num_js1, STAT = ierr)
        IF (ierr /= 0) CALL errore('deallocate_aniso_iaxis', 'Error deallocating num_js1', 1)
        !
        IF (icoulomb > 0) THEN
          ! HM: deallocate arrays for outer bands in sum_eliashberg_aniso_iaxis
          DEALLOCATE(adeltai_cl, STAT = ierr)
          IF (ierr /= 0) CALL errore('deallocate_aniso_iaxis', 'Error deallocating adeltai_cl', 1)
          DEALLOCATE(adeltaip_cl, STAT = ierr)
          IF (ierr /= 0) CALL errore('deallocate_aniso_iaxis', 'Error deallocating adeltaip_cl', 1)
          DEALLOCATE(w_stat, STAT = ierr)
          IF (ierr /= 0) CALL errore('deallocate_aniso_iaxis', 'Error deallocating w_stat', 1)
          DEALLOCATE(weight_cl, STAT = ierr)
          IF (ierr /= 0) CALL errore('deallocate_aniso_iaxis', 'Error deallocating weight_cl', 1)
          !!
          ! HM: deallocate ir-related arrays for outer bands in sum_eliashberg_aniso_iaxis
          DEALLOCATE(ir_giw_cl, STAT = ierr)
          IF (ierr /= 0) CALL errore('deallocate_aniso_iaxis', 'Error deallocating ir_giw_cl', 1)
          IF (positive_matsu) THEN
            DEALLOCATE(ir_gl_cl_d, STAT = ierr)
            IF (ierr /= 0) CALL errore('deallocate_aniso_iaxis', 'Error deallocating ir_gl_cl_d', 1)
            DEALLOCATE(ir_gtau_cl_d, STAT = ierr)
            IF (ierr /= 0) CALL errore('deallocate_aniso_iaxis', 'Error deallocating ir_gtau_cl_d', 1)
          ELSE
            DEALLOCATE(ir_gl_cl, STAT = ierr)
            IF (ierr /= 0) CALL errore('deallocate_aniso_iaxis', 'Error deallocating ir_gl_cl', 1)
            DEALLOCATE(ir_gtau_cl, STAT = ierr)
            IF (ierr /= 0) CALL errore('deallocate_aniso_iaxis', 'Error deallocating ir_gtau_cl', 1)
          ENDIF
          DEALLOCATE(num_js2, STAT = ierr)
          IF (ierr /= 0) CALL errore('deallocate_aniso_iaxis', 'Error deallocating num_js2', 1)
          DEALLOCATE(num_js3, STAT = ierr)
          IF (ierr /= 0) CALL errore('deallocate_aniso_iaxis', 'Error deallocating num_js3', 1)
        ENDIF
        !
        IF (iverbosity == 4) THEN
          ! HM: deallocate gl_abs, fl_abs, and knll_abs
          DEALLOCATE(gl_abs, STAT = ierr)
          IF (ierr /= 0) CALL errore('deallocate_aniso_iaxis', 'Error allocating gl_abs', 1)
          DEALLOCATE(fl_abs, STAT = ierr)
          IF (ierr /= 0) CALL errore('deallocate_aniso_iaxis', 'Error allocating fl_abs', 1)
          DEALLOCATE(knll_abs, STAT = ierr)
          IF (ierr /= 0) CALL errore('deallocate_aniso_iaxis', 'Error allocating knll_abs', 1)
        ENDIF
      ELSEIF (gridsamp == 3) THEN
        ! HM: deallocate FFT-related arrays in sum_eliashberg_aniso_iaxis
        DEALLOCATE(fft_in1, STAT = ierr)
        IF (ierr /= 0) CALL errore('deallocate_aniso_iaxis', 'Error deallocating fft_in1', 1)
        DEALLOCATE(fft_out1, STAT = ierr)
        IF (ierr /= 0) CALL errore('deallocate_aniso_iaxis', 'Error deallocating fft_out1', 1)
        DEALLOCATE(fft_in2, STAT = ierr)
        IF (ierr /= 0) CALL errore('deallocate_aniso_iaxis', 'Error deallocating fft_in2', 1)
        DEALLOCATE(fft_out2, STAT = ierr)
        IF (ierr /= 0) CALL errore('deallocate_aniso_iaxis', 'Error deallocating fft_out2', 1)
      ENDIF
    ENDIF
    !
    RETURN
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE deallocate_aniso_iaxis
    !-----------------------------------------------------------------------
    !
    !----------------------------------------------------------------------
    SUBROUTINE deallocate_aniso_raxis(iset)
    !----------------------------------------------------------------------
    !!
    !!  deallocates the variables allocated for real-axis solutions
    !!
    USE input,         ONLY : fbw, lacon
    USE supercond_common,     ONLY : ws, delta, znorm, adelta, adeltap, &
                              aznorm, aznormp, gp, gm, adsumi, &
                              azsumi, a2fij, lacon_fly, shift, ashift
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: iset
    !! Define what set of variables are deallocated
    !
    INTEGER :: ierr
    !! Error status
    !
    IF (iset == 1 .OR. iset == 3) THEN
      IF (lacon) THEN
        ! analytic_cont_aniso
        DEALLOCATE(adeltap, STAT = ierr)
        IF (ierr /= 0) CALL errore('deallocate_aniso_raxis', 'Error deallocating adeltap', 1)
        DEALLOCATE(aznormp, STAT = ierr)
        IF (ierr /= 0) CALL errore('deallocate_aniso_raxis', 'Error deallocating aznormp', 1)
        DEALLOCATE(gp, STAT = ierr)
        IF (ierr /= 0) CALL errore('deallocate_aniso_raxis', 'Error deallocating gp', 1)
        DEALLOCATE(gm, STAT = ierr)
        IF (ierr /= 0) CALL errore('deallocate_aniso_raxis', 'Error deallocating gm', 1)
        ! kernel_aniso_analytic_cont
        DEALLOCATE(adsumi, STAT = ierr)
        IF (ierr /= 0) CALL errore('deallocate_aniso_raxis', 'Error deallocating adsumi', 1)
        DEALLOCATE(azsumi, STAT = ierr)
        IF (ierr /= 0) CALL errore('deallocate_aniso_raxis', 'Error deallocating azsumi', 1)
        ! evaluate_a2fij
        IF (.NOT. lacon_fly) THEN
          DEALLOCATE(a2fij, STAT = ierr)
          IF (ierr /= 0) CALL errore('deallocate_aniso_raxis', 'Error deallocating a2fij', 1)
        ENDIF
      ENDIF
    ENDIF
    !
    IF (iset == 2 .OR. iset == 3) THEN
      ! gen_freqgrid_iaxis
      DEALLOCATE(ws, STAT = ierr)
      IF (ierr /= 0) CALL errore('deallocate_aniso_raxis', 'Error deallocating ws', 1)
      !
      ! pade_cont_aniso or analytic_cont_aniso
      DEALLOCATE(delta, STAT = ierr)
      IF (ierr /= 0) CALL errore('deallocate_aniso_raxis', 'Error deallocating delta', 1)
      DEALLOCATE(adelta, STAT = ierr)
      IF (ierr /= 0) CALL errore('deallocate_aniso_raxis', 'Error deallocating adelta', 1)
      DEALLOCATE(znorm, STAT = ierr)
      IF (ierr /= 0) CALL errore('deallocate_aniso_raxis', 'Error deallocating znorm', 1)
      DEALLOCATE(aznorm, STAT = ierr)
      IF (ierr /= 0) CALL errore('deallocate_aniso_raxis', 'Error deallocating aznorm', 1)
      !
      IF (fbw) THEN
        ! pade_cont_aniso fbw-related arrays
        DEALLOCATE(shift, STAT = ierr)
        IF (ierr /= 0) CALL errore('deallocate_aniso_raxis', 'Error deallocating shift', 1)
        DEALLOCATE(ashift, STAT = ierr)
        IF (ierr /= 0) CALL errore('deallocate_aniso_raxis', 'Error deallocating ashift', 1)
      ENDIF
      !
    ENDIF
    !
    RETURN
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE deallocate_aniso_raxis
    !-----------------------------------------------------------------------
    !
    !----------------------------------------------------------------------
    SUBROUTINE deallocate_aniso()
    !----------------------------------------------------------------------
    !!
    !!  deallocates the variables allocated by eliashberg_init,
    !!  eliashberg_grid, read_frequencies, read_eigenvalues, read_kqmap,
    !!  read_ephmat and evaluate_a2f_lambda
    !!
    USE global_var,       ONLY : wf, wqf, xqf, gtemp, bztoibz
    USE input,            ONLY : fbw, icoulomb, gridsamp
    USE supercond_common, ONLY : ekfs, xkfs, wkfs, g2, w0g, a2f_tmp, &
                                 ixkff, ixkqf, ixqfs, nqfs, memlt_pool, &
                                 wsph, nsiw, agap, ibnd_kfs_all_to_kfs, &
                                 ibnd_kfs_to_kfs_all, ekfs_all, &
                                 xkfs_all, wkfs_all
    USE supercond_coul,   ONLY : deallocate_coulomb
    !
    IMPLICIT NONE
    !
    INTEGER :: ierr
    !! Error status
    !
    IF (fbw .AND. (gridsamp == 2) .AND. (icoulomb > 0)) THEN
      CALL deallocate_coulomb()
    ENDIF
    ! eliashberg_init
    DEALLOCATE(gtemp, STAT = ierr)
    IF (ierr /= 0) CALL errore('deallocate_aniso', 'Error deallocating gtemp', 1)
    ! eliashberg_grid
    DEALLOCATE(nsiw, STAT = ierr)
    IF (ierr /= 0) CALL errore('deallocate_aniso', 'Error deallocating nsiw', 1)
    ! read_frequencies
    DEALLOCATE(wf, STAT = ierr)
    IF (ierr /= 0) CALL errore('deallocate_aniso', 'Error deallocating wf', 1)
    DEALLOCATE(wqf, STAT = ierr)
    IF (ierr /= 0) CALL errore('deallocate_aniso', 'Error deallocating wqf', 1)
    DEALLOCATE(xqf, STAT = ierr)
    IF (ierr /= 0) CALL errore('deallocate_aniso', 'Error deallocating xqf', 1)
    ! read_eigenvalues
    DEALLOCATE(ibnd_kfs_all_to_kfs, STAT = ierr)
    IF (ierr /= 0) CALL errore('deallocate_aniso', 'Error deallocating ibnd_kfs_all_to_kfs', 1)
    DEALLOCATE(ibnd_kfs_to_kfs_all, STAT = ierr)
    IF (ierr /= 0) CALL errore('deallocate_aniso', 'Error deallocating ibnd_kfs_to_kfs_all', 1)
    DEALLOCATE(ekfs, STAT = ierr)
    IF (ierr /= 0) CALL errore('deallocate_aniso', 'Error deallocating ekfs', 1)
    DEALLOCATE(xkfs, STAT = ierr)
    IF (ierr /= 0) CALL errore('deallocate_aniso', 'Error deallocating xkfs', 1)
    DEALLOCATE(wkfs, STAT = ierr)
    IF (ierr /= 0) CALL errore('deallocate_aniso', 'Error deallocating wkfs', 1)
    DEALLOCATE(ekfs_all, STAT = ierr)
    IF (ierr /= 0) CALL errore('deallocate_aniso', 'Error deallocating ekfs_all', 1)
    DEALLOCATE(xkfs_all, STAT = ierr)
    IF (ierr /= 0) CALL errore('deallocate_aniso', 'Error deallocating xkfs_all', 1)
    DEALLOCATE(wkfs_all, STAT = ierr)
    IF (ierr /= 0) CALL errore('deallocate_aniso', 'Error deallocating wkfs_all', 1)
    DEALLOCATE(w0g, STAT = ierr)
    IF (ierr /= 0) CALL errore('deallocate_aniso', 'Error deallocating w0g', 1)
    ! read_kqmap
    DEALLOCATE(ixkff, STAT = ierr)
    IF (ierr /= 0) CALL errore('deallocate_aniso', 'Error deallocating ixkff', 1)
    DEALLOCATE(bztoibz, STAT = ierr)
    IF (ierr /= 0) CALL errore('deallocate_aniso', 'Error deallocating bztoibz', 1)
    DEALLOCATE(ixkqf, STAT = ierr)
    IF (ierr /= 0) CALL errore('deallocate_aniso', 'Error deallocating ixkqf', 1)
    DEALLOCATE(ixqfs, STAT = ierr)
    IF (ierr /= 0) CALL errore('deallocate_aniso', 'Error deallocating ixqfs', 1)
    DEALLOCATE(nqfs, STAT = ierr)
    IF (ierr /= 0) CALL errore('deallocate_aniso', 'Error deallocating nqfs', 1)
    DEALLOCATE(memlt_pool, STAT = ierr)
    IF (ierr /= 0) CALL errore('deallocate_aniso', 'Error deallocating memlt_pool', 1)
    ! read_ephmat
    DEALLOCATE(g2, STAT = ierr)
    IF (ierr /= 0) CALL errore('deallocate_aniso', 'Error deallocating g2', 1)
    ! read_a2f
    DEALLOCATE(wsph, STAT = ierr)
    IF (ierr /= 0) CALL errore('deallocate_aniso', 'Error deallocating wsph', 1)
    DEALLOCATE(a2f_tmp, STAT = ierr)
    IF (ierr /= 0) CALL errore('deallocate_aniso', 'Error deallocating a2f_tmp', 1)
    ! sum_eliashberg_aniso_iaxis
    DEALLOCATE(agap, STAT = ierr)
    IF (ierr /= 0) CALL errore('deallocate_aniso', 'Error deallocating agap', 1)
    !
    RETURN
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE deallocate_aniso
    !-----------------------------------------------------------------------
    !
  !-----------------------------------------------------------------------
  END MODULE supercond_aniso
  !-----------------------------------------------------------------------
