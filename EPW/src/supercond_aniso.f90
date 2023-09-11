  !
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
    !!
    USE kinds,         ONLY : DP
    USE control_flags, ONLY : iverbosity
    USE epwcom,        ONLY : nsiter, nstemp, broyden_beta, broyden_ndim, &
                              limag, lpade, lacon, fsthick, imag_read, npade, &
                              fbw 
    USE elph2,         ONLY : gtemp
    USE eliashbergcom, ONLY : nsw, nsiw, adelta, adeltap, adeltai, adeltaip, &
                              aznormi, aznormip, ashifti, ashiftip, &
                              nkfs, nbndfs, ekfs, ef0, agap
    USE supercond,     ONLY : dos_quasiparticle, gen_freqgrid_iaxis, &
                              eliashberg_grid
    USE constants_epw, ONLY : kelvin2eV, ci, zero, czero
    USE io_global,     ONLY : stdout
    USE mp_global,     ONLY : inter_pool_comm, my_pool_id
    USE mp,            ONLY : mp_barrier, mp_sum
    USE io_eliashberg, ONLY : eliashberg_write_iaxis, eliashberg_read_aniso_iaxis, &
                              eliashberg_write_raxis
    USE utilities,     ONLY : mix_wrap
    USE low_lvl,       ONLY : mem_size_eliashberg
    USE printing,      ONLY : prtheader_supercond
    USE division,      ONLY : fkbounds
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
    INTEGER(8) :: imelt
    !! Counter memory
    INTEGER :: ierr
    !! Error status
    !
    REAL(KIND = DP) :: tcpu
    !! cpu time
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
    !
    CALL start_clock('aniso_iaxis')
    !
    CALL eliashberg_grid()
    !
    CALL fkbounds(nkfs, lower_bnd, upper_bnd)
    !
    nks = upper_bnd - lower_bnd + 1
    !
    DO itemp = 1, nstemp ! loop over temperature
      !
      CALL prtheader_supercond(itemp, 1)
      CALL start_clock('iaxis_imag')
      CALL gen_freqgrid_iaxis(itemp)
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
        ENDIF ! fbw
        !
        iter = 1
        conv = .FALSE.
        DO WHILE (.NOT. conv .AND. iter <= nsiter)
          CALL sum_eliashberg_aniso_iaxis(itemp, iter, conv)
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
          CALL deallocate_aniso_iaxis(3)
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
        N = npade * nsiw(itemp) / 100
        IF (mod(N, 2) /= 0 ) N = N + 1
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
          CALL deallocate_aniso_iaxis(2)
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
          cname = 'acon'
          CALL eliashberg_write_raxis(itemp, cname) 
          !
          CALL dos_quasiparticle(itemp) ! - change to work on pool
          !
          CALL stop_clock('raxis_acon')
          CALL print_clock('raxis_acon')
          WRITE(stdout, '(a)') ' '
        ELSEIF (.NOT. conv .AND. (iter - 1) == nsiter) THEN
          CALL deallocate_aniso_iaxis(2)      
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
        CALL deallocate_aniso_iaxis(2)
        !
        IF (fbw) THEN
          ! remove memory allocated for wsi, wsn, adeltai, aznormi, ashifti
          imelt = (2 + 3 * nbndfs * nks) * nsiw(itemp)
        ELSE
          ! remove memory allocated for wsi, adeltai, aznormi
          imelt = (1 + 2 * nbndfs * nks) * nsiw(itemp)
        ENDIF
        CALL mem_size_eliashberg(2, -imelt)
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
    SUBROUTINE sum_eliashberg_aniso_iaxis(itemp, iter, conv)
    !-----------------------------------------------------------------------
    !!
    !! This routine solves the anisotropic Eliashberg equations on the imaginary-axis
    !!
    !! SH: Modified to allow for fbw calculations; and 
    !!       re-ordered "deltai, ..." arrays' indices for efficiency (Nov 2021).
    !!
    USE kinds,             ONLY : DP
    USE elph2,             ONLY : wqf, gtemp
    USE epwcom,            ONLY : nsiter, nstemp, muc, conv_thr_iaxis, fsthick, fbw, muchem, imag_read
    USE eliashbergcom,     ONLY : ixkqf, ixqfs, nqfs, wkfs, w0g, ekfs, nkfs, nbndfs, dosef, ef0, &
                                  nsiw, wsn, wsi, wsphmax, gap0, agap, akeri, limag_fly, &
                                  deltai, znormi, shifti, adeltai, adeltaip, aznormi, aznormip, &
                                  naznormi, ashifti, ashiftip, muintr
    USE constants_epw,     ONLY : kelvin2eV, zero, one
    USE constants,         ONLY : pi
    USE io_global,         ONLY : stdout, ionode_id
    USE mp_global,         ONLY : inter_pool_comm
    USE mp_world,          ONLY : mpime
    USE mp,                ONLY : mp_bcast, mp_barrier, mp_sum
    USE division,          ONLY : fkbounds
    USE low_lvl,           ONLY : mem_size_eliashberg, memlt_eliashberg
    !
    IMPLICIT NONE
    !
    LOGICAL, INTENT(inout) :: conv
    !! True if the calculation is converged
    INTEGER, INTENT(in) :: itemp
    !! Counter on temperature index
    INTEGER, INTENT(in) :: iter
    !! Counter on iteration steps
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
    INTEGER :: lower_bnd, upper_bnd
    !! Lower/upper bound index after k paral
    INTEGER :: nks
    !! Number of k points per pool
    INTEGER :: ibnd
    !! Counter on bands at k
    INTEGER :: jbnd
    !! Counter on bands at k+q
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
    REAL(KIND = DP) :: absdelta, reldelta, errdelta
    !! Errors in supercond. gap
    REAL(KIND = DP) :: weight
    !! Factor in supercond. equations
    REAL(KIND = DP) :: esqrt, zesqrt, desqrt, sesqrt
    !! Temporary variables
    REAL(KIND = DP), SAVE :: nel, nstate
    !! mu_inter parameters
    REAL(KIND = DP) :: inv_dos
    !! Invese dos inv_dos = 1/dosef. Defined for efficiency reason
    REAL(KIND = DP) :: omega
    !! sqrt{w^2+\delta^2}
    REAL(KIND = DP) :: dFE
    !! free energy difference between supercond and normal states
    REAL(KIND = DP), ALLOCATABLE, SAVE :: inv_wsi(:)
    !! Invese imaginary freq. inv_wsi = 1/wsi. Defined for efficiency reason
    REAL(KIND = DP), ALLOCATABLE, SAVE :: deltaold(:)
    !! supercond. gap from previous iteration
    !
    inv_dos = one / dosef
    !
    CALL fkbounds(nkfs, lower_bnd, upper_bnd)
    !
    nks = upper_bnd - lower_bnd + 1
    !
    IF (iter == 1) THEN
      IF (itemp == 1 .OR. (itemp == 2 .AND. imag_read)) THEN
        ! SH: calculate the input parameters for mu_inter
        nel = zero
        nstate = zero
        DO ik = lower_bnd, upper_bnd
          DO ibnd = 1, nbndfs
            IF (ABS(ekfs(ibnd, ik) - ef0) < fsthick) THEN
              ! HP: The initial guess is based on the FD dist. at 0 K.
              nstate = nstate + 0.5d0 * wkfs(ik)
              IF ((ekfs(ibnd, ik) - ef0) < zero) nel = nel + wkfs(ik)
            ENDIF
          ENDDO
        ENDDO
        CALL mp_sum(nstate, inter_pool_comm)
        CALL mp_sum(nel, inter_pool_comm)
        CALL mp_barrier(inter_pool_comm)
        !
      ENDIF
      !
      ! RM - adeltai, aznormi, naznormi, are defined per k-points per pool
      !
      ! get the size of required memory for inv_wsi, deltaold, deltai, znormi, 
      ! adeltai, adeltaip, aznormi, naznormi
      imelt = (4 + 3 * nbndfs * nks + nbndfs * nkfs) * nsiw(itemp)
      CALL mem_size_eliashberg(2, imelt)
      !
      ALLOCATE(inv_wsi(nsiw(itemp)), STAT = ierr)
      IF (ierr /= 0) CALL errore('sum_eliashberg_aniso_iaxis', 'Error allocating inv_wsi', 1)
      !      
      ALLOCATE(deltai(nsiw(itemp)), STAT = ierr)
      IF (ierr /= 0) CALL errore('sum_eliashberg_aniso_iaxis', 'Error allocating deltai', 1)
      ALLOCATE(znormi(nsiw(itemp)), STAT = ierr)
      IF (ierr /= 0) CALL errore('sum_eliashberg_aniso_iaxis', 'Error allocating znormi', 1)
      !
      ALLOCATE(adeltai(nsiw(itemp), nbndfs, lower_bnd:upper_bnd), STAT = ierr)
      IF (ierr /= 0) CALL errore('sum_eliashberg_aniso_iaxis', 'Error allocating adeltai', 1)
      ALLOCATE(aznormi(nsiw(itemp), nbndfs, lower_bnd:upper_bnd), STAT = ierr)
      IF (ierr /= 0) CALL errore('sum_eliashberg_aniso_iaxis', 'Error allocating aznormi', 1)
      ALLOCATE(naznormi(nsiw(itemp), nbndfs, lower_bnd:upper_bnd), STAT = ierr)
      IF (ierr /= 0) CALL errore('sum_eliashberg_aniso_iaxis', 'Error allocating naznormi', 1)
      !
      ALLOCATE(adeltaip(nsiw(itemp), nbndfs, nkfs), STAT = ierr)
      IF (ierr /= 0) CALL errore('sum_eliashberg_aniso_iaxis', 'Error allocating adeltaip', 1)
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
      ENDIF ! fbw
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
              IF (wsi(iw) < 2.d0 * wsphmax) THEN
                adeltaip(iw, ibnd, ik) = gap0
              ELSE
                adeltaip(iw, ibnd, ik) = zero
              ENDIF
            ENDDO
          ENDIF
        ENDDO ! ibnd
      ENDDO ! ik
      !
      CALL memlt_eliashberg(itemp, 'imag')
      IF (.NOT. limag_fly) CALL kernel_aniso_iaxis(itemp)
      !
    ENDIF ! iter
    !
    ! SH: for the case of fbw runs
    IF (fbw) THEN
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
                      ! Eqs.(15-17) in Margine and Giustino, PRB 87, 024505 (2013)
                      ! using kernelm and kernelp the sum over |wp| < wscut
                      ! is rewritten as a sum over iwp = 1, nsiw(itemp)
                      aznormi(iw, ibnd, ik) = aznormi(iw, ibnd, ik) + zesqrt * kernelm
                      adeltai(iw, ibnd, ik) = adeltai(iw, ibnd, ik) + desqrt * (kernelp - 2.d0 * muc)
                      ashifti(iw, ibnd, ik) = ashifti(iw, ibnd, ik) + sesqrt * kernelp
                    ENDDO ! iw
                  ENDDO ! iwp
                ENDIF
              ENDDO ! jbnd
            ENDDO ! iq
          ENDIF
        ENDDO ! ibnd
      ENDDO ! ik
      !
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
      CALL mp_sum(deltai,   inter_pool_comm)
      CALL mp_sum(znormi,   inter_pool_comm)
      CALL mp_sum(shifti,   inter_pool_comm)
      CALL mp_barrier(inter_pool_comm)
      !
    ELSE ! not fbw
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
                      adeltai(iw, ibnd, ik) = adeltai(iw, ibnd, ik) + desqrt * (kernelp - 2.d0 * muc)
                    ENDDO ! iw
                  ENDDO ! iwp
                ENDIF
              ENDDO ! jbnd
            ENDDO ! iq
          ENDIF
        ENDDO ! ibnd
      ENDDO ! ik
      !
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
    ENDIF ! fbw
    !
    IF (mpime == ionode_id) THEN
      !
      IF (iter == 1) THEN
        ALLOCATE(deltaold(nsiw(itemp)), STAT = ierr)
        IF (ierr /= 0) CALL errore('sum_eliashberg_aniso_iaxis', 'Error allocating deltaold', 1)
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
      IF (fbw)                       WRITE(stdout, '(5x, i6, 5ES15.6)') &
        iter, errdelta, znormi(1), deltai(1) * 1000.d0, shifti(1) * 1000.d0, muintr
      IF (.NOT. fbw)                 WRITE(stdout, '(5x, i6, 3ES15.6)') &
        iter, errdelta, znormi(1), deltai(1) * 1000.d0
      !      WRITE(stdout, '(5x, a, i6, a, ES15.6, a, ES15.6, a, ES15.6)') 'iter = ', iter, &
      !                    '   ethr = ', errdelta, '   znormi(1) = ', znormi(1), &
      !                    '   deltai(1) = ', deltai(1)
      !
      IF (errdelta < conv_thr_iaxis) conv = .TRUE.
      IF (conv .OR. iter == nsiter) THEN
        gap0 = deltai(1)
      ENDIF
      !
      IF (conv .OR. iter == nsiter) THEN
        DEALLOCATE(deltaold, STAT = ierr)
        IF (ierr /= 0) CALL errore('sum_eliashberg_aniso_iaxis', 'Error deallocating deltaold', 1)
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
      ! SH: write the chemical potential for fbw runs
      WRITE(stdout, '(5x, a, i3, a, ES20.10, a)') &
        'Chemical potential (itemp = ', itemp, ') = ', muintr, ' eV'
      WRITE(stdout, '(a)') ' '
      !
      ! Compute the free energy difference between the superconducting and normal states      
      dFE = zero
      DO ik = lower_bnd, upper_bnd
        DO ibnd = 1, nbndfs
          IF (ABS(ekfs(ibnd, ik) - ef0) < fsthick) THEN
            DO iw = 1, nsiw(itemp)
              weight = 0.5d0 * wkfs(ik) * w0g(ibnd,ik)
              omega = DSQRT(wsi(iw) * wsi(iw) + adeltai(iw, ibnd, ik) * adeltai(iw, ibnd, ik))
              dFE = dFE - weight * (omega - wsi(iw)) &
                * (aznormi(iw, ibnd, ik) - naznormi(iw, ibnd, ik) * wsi(iw) / omega)
            ENDDO
          ENDIF
        ENDDO
      ENDDO
      ! collect contributions from all pools
      CALL mp_sum(dFE, inter_pool_comm)
      CALL mp_barrier(inter_pool_comm)
      !
      dFE = dFE * pi * gtemp(itemp)
      !
      WRITE(stdout, '(5x, a, i3, a, f8.3, a, a, f12.6, a)') &
              'Temp (itemp = ', itemp, ') = ', gtemp(itemp) / kelvin2eV, ' K', &
              '  Free energy = ', dFE * 1000.0, ' meV'
      WRITE(stdout, '(a)') ' '
      !
      DEALLOCATE(inv_wsi, STAT = ierr)
      IF (ierr /= 0) CALL errore('sum_eliashberg_aniso_iaxis', 'Error deallocating inv_wsi', 1)
      DEALLOCATE(naznormi, STAT = ierr)
      IF (ierr /= 0) CALL errore('sum_eliashberg_aniso_iaxis', 'Error deallocating naznormi', 1)      
      !      
      ! remove memory allocated for inv_wsi, deltaold, naznormi
      imelt = (2 + nbndfs * nks) * nsiw(itemp) 
      CALL mem_size_eliashberg(2, -imelt)
      !
      IF (.NOT. limag_fly) THEN
        !
        DEALLOCATE(akeri, STAT = ierr)
        IF (ierr /= 0) CALL errore('sum_eliashberg_aniso_iaxis', 'Error deallocating akeri', 1)
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
    END SUBROUTINE sum_eliashberg_aniso_iaxis
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
    USE elph2,         ONLY : wqf, wf, gtemp
    USE epwcom,        ONLY : nqstep, nstemp, degaussq, nsiter, conv_thr_racon, fsthick, &
                              lpade, eps_acustic
    USE eliashbergcom, ONLY : ixkqf, ixqfs, nqfs, wkfs, w0g, ekfs, nkfs, nbndfs, dosef, ef0, &
                              nsw, dwsph, ws, wsph, gap0, agap, gp, gm, adsumi, azsumi, a2fij, &
                              g2, lacon_fly, delta, znorm, adelta, adeltap, aznorm, aznormp
    USE supercond,     ONLY : gamma_acont
    USE constants_epw, ONLY : ci, zero, one, czero, cone
    USE constants,     ONLY : pi
    USE io_global,     ONLY : stdout, ionode_id
    USE mp_global,     ONLY : inter_pool_comm, my_pool_id
    USE mp_world,      ONLY : mpime
    USE mp,            ONLY : mp_bcast, mp_barrier, mp_sum
    USE division,      ONLY : fkbounds
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
          DO iw = 1, nsw ! loop over omega
            az2 = aznormp(iw, ibnd, ik) * aznormp(iw, ibnd, ik)
            ad2 = adeltap(iw, ibnd, ik) * adeltap(iw, ibnd, ik)
            root = SQRT(az2 * (ws(iw) * ws(iw) - ad2))
            IF (AIMAG(root) < zero) & 
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
                    IF (wf(imode, iq0) > eps_acustic) THEN
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
      ! remove memory allocated for deltaold, gp, gm, adsumi, azsumi
      imelt = 2 * (1 + nqstep + nks * nbndfs) * nsw
      CALL mem_size_eliashberg(2, -imelt)
      !
      IF (.NOT. lacon_fly) THEN
        ! remove memory allocated for a2fij
        imelt = nks * MAXVAL(nqfs(:)) * nbndfs**2 * nqstep
        CALL mem_size_eliashberg(2, -imelt)
      ENDIF
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
    USE epwcom,        ONLY : fsthick, fbw
    USE eliashbergcom, ONLY : nsw, ws, nsiw, wsi, delta, znorm, &
                              adelta, aznorm, adeltai, aznormi, &
                              wkfs, dosef, w0g, nkfs, nbndfs, ef0, ekfs, &
                              shift, ashift, ashifti
    USE utilities,     ONLY : pade_coeff, pade_eval
    USE constants_epw, ONLY : cone, ci, zero, czero, one
    USE io_global,     ONLY : stdout, ionode_id
    USE mp_global,     ONLY : inter_pool_comm
    USE mp_world,      ONLY : mpime
    USE mp,            ONLY : mp_barrier, mp_sum
    USE division,      ONLY : fkbounds
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
            DO iw = 1, N
              z(iw) = ci * wsi(iw)
              u(iw) = cone * adeltai(iw, ibnd, ik)
              v(iw) = cone * aznormi(iw, ibnd, ik)
              w(iw) = cone * ashifti(iw, ibnd, ik)
            ENDDO
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
    !!
    USE kinds,         ONLY : DP
    USE epwcom,        ONLY : fsthick
    USE elph2,         ONLY : gtemp
    USE eliashbergcom, ONLY : nkfs, nbndfs, nsiw, akeri, ekfs, ef0, ixkqf, ixqfs, nqfs, wsn
    USE constants_epw, ONLY : zero
    USE constants,     ONLY : pi
    USE division,      ONLY : fkbounds
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
    USE elph2,         ONLY : wf
    USE epwcom,        ONLY : eps_acustic
    USE eliashbergcom, ONLY : ixqfs, g2, dosef
    USE constants_epw, ONLY : zero
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
      IF (wf(imode, iq0) > eps_acustic) THEN
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
    USE epwcom,        ONLY : nqstep
    USE eliashbergcom, ONLY : a2fij, wsph, dwsph
    USE constants_epw, ONLY : zero
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
    USE elph2,         ONLY : wqf
    USE epwcom,        ONLY : muc, fsthick
    USE eliashbergcom, ONLY : nsw, nsiw, ws, wsi, adeltai, nkfs, nbndfs, dosef, ixkqf, ixqfs, nqfs, &
                              w0g, ekfs, ef0, adsumi, azsumi
    USE constants_epw, ONLY : zero, one
    USE mp,            ONLY : mp_barrier, mp_sum
    USE mp_global,     ONLY : inter_pool_comm
    USE division,      ONLY : fkbounds
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
                    adsumi(iw, ibnd, ik) = adsumi(iw, ibnd, ik) + desqrt * (kernelr - 2.d0 * muc)
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
    USE elph2,         ONLY : wf
    USE epwcom,        ONLY : eps_acustic
    USE eliashbergcom, ONLY : ixqfs, g2, dosef
    USE constants_epw, ONLY : ci, czero
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
      IF (wf(imode, iq0) > eps_acustic) THEN
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
    USE epwcom,        ONLY : nqstep
    USE eliashbergcom, ONLY : a2fij, dwsph, wsph
    USE constants_epw, ONLY : ci, czero
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
    USE elph2,         ONLY : wf
    USE epwcom,        ONLY : fsthick, eps_acustic, nqstep, degaussq
    USE eliashbergcom, ONLY : nkfs, nbndfs, g2, a2fij, ixkqf, ixqfs, nqfs, ekfs, ef0, &
                              dosef, wsph
    USE constants_epw, ONLY : zero, one
    USE division,      ONLY : fkbounds
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
                  IF (wf(imode, iq0) > eps_acustic) THEN
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
    SUBROUTINE mu_inter_aniso(itemp, muintr, nel, nstate)
    !-----------------------------------------------------------------------
    !!
    !! SH: To find the superconducting state chemical potential
    !!       using the Newton-Raphson method (Nov 2021).
    !!
    USE kinds,          ONLY : DP
    USE eliashbergcom,  ONLY : ekfs, ef0, nkfs, nbndfs, wsi, nsiw, &
                               adeltaip, aznormip, ashiftip, wkfs
    USE elph2,          ONLY : gtemp
    USE epwcom,         ONLY : fsthick, nsiter, broyden_beta
    USE constants_epw,  ONLY : zero, eps6
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
    INTEGER :: lower_bnd, upper_bnd
    !! Temporary variables
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
              fmu = fmu + 2.d0 * wkfs(ik) * delta / theta
              dmu = dmu + 2.d0 * wkfs(ik) * (2.d0 * delta**2.d0 - theta) / (theta**2.d0)
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
      CALL errore('mu_inter_aniso', 'Error failed to find the mu_inter_aniso value',1)
    !
    muintr = muout
    !
    RETURN
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE mu_inter_aniso
    !-----------------------------------------------------------------------
    !
    !----------------------------------------------------------------------
    SUBROUTINE deallocate_aniso_iaxis(iset)
    !----------------------------------------------------------------------
    !!
    !!  deallocates the variables allocated for imag-axis solutions
    !!
    !
    USE epwcom,        ONLY : fbw
    USE eliashbergcom, ONLY : wsi, wsn, deltai, znormi, adeltai, adeltaip, & 
                              aznormi, aznormip, shifti, ashifti, ashiftip
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
        IF (ierr /= 0) CALL errore('deallocate_aniso_iaxis_fbw', 'Error deallocating shifti', 1)
        DEALLOCATE(aznormip, STAT = ierr)
        IF (ierr /= 0) CALL errore('deallocate_aniso_iaxis_fbw', 'Error deallocating aznormip', 1)
        DEALLOCATE(ashiftip, STAT = ierr)
        IF (ierr /= 0) CALL errore('deallocate_aniso_iaxis_fbw', 'Error deallocating ashiftip', 1)
      ENDIF
    ENDIF
    !
    IF (iset == 2 .OR. iset == 3) THEN
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
        IF (ierr /= 0) CALL errore('deallocate_aniso_iaxis_fbw', 'Error deallocating ashifti', 1)
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
    USE epwcom,        ONLY : fbw, lacon
    USE eliashbergcom, ONLY : ws, delta, znorm, adelta, adeltap, & 
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
    USE elph2,         ONLY : wf, wqf, xqf, gtemp
    USE eliashbergcom, ONLY : ekfs, xkfs, wkfs, g2, w0g, a2f_iso, &
                              ixkff, ixkqf, ixqfs, nqfs, memlt_pool, &
                              wsph, nsiw, agap
    !
    IMPLICIT NONE
    !
    INTEGER :: ierr
    !! Error status
    !
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
    DEALLOCATE(ekfs, STAT = ierr)
    IF (ierr /= 0) CALL errore('deallocate_aniso', 'Error deallocating ekfs', 1)
    DEALLOCATE(xkfs, STAT = ierr)
    IF (ierr /= 0) CALL errore('deallocate_aniso', 'Error deallocating xkfs', 1)
    DEALLOCATE(wkfs, STAT = ierr)
    IF (ierr /= 0) CALL errore('deallocate_aniso', 'Error deallocating wkfs', 1)
    DEALLOCATE(w0g, STAT = ierr)
    IF (ierr /= 0) CALL errore('deallocate_aniso', 'Error deallocating w0g', 1)
    ! read_kqmap
    DEALLOCATE(ixkff, STAT = ierr)
    IF (ierr /= 0) CALL errore('deallocate_aniso', 'Error deallocating ixkff', 1)
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
    DEALLOCATE(a2f_iso, STAT = ierr)
    IF (ierr /= 0) CALL errore('deallocate_aniso', 'Error deallocating a2f_iso', 1)
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
