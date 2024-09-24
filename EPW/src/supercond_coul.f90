  !
  ! Copyright (C) 2024 EPW-Collaboration
  !
  ! This file is distributed under the terms of the GNU General Public
  ! License. See the file `LICENSE' in the root directory of the
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .
  !
  !----------------------------------------------------------------------
  MODULE supercond_coul
  !----------------------------------------------------------------------
  !!
  !! This module contains all the coulomb parts of the superconductivity function of EPW
  !!
  IMPLICIT NONE
  !
  CONTAINS
    !-----------------------------------------------------------------------
    SUBROUTINE read_eigen_cl()
    !-----------------------------------------------------------------------
    !!
    !! Read the eigenvalues from the file created by PP/bands2epw.x
    !! k points written in the file are assumed to be in cartesian coord.
    !! After reading them, this routine converts them to crystal coordinates.
    !!
    USE kinds,         ONLY : DP
    USE cell_base,     ONLY : at, alat
    USE symm_base,     ONLY : nrot, nsym, t_rev, s
    USE control_flags, ONLY : noinv
    USE noncollin_module, ONLY : noncolin, lspinorb
    USE lsda_mod,      ONLY : lsda!, nspin
    USE supercond_common, ONLY : ef0, nk1_cl, nk2_cl, nk3_cl, nbnd_cl, nkstot_cl, &
                              wk_cl, xk_cl, ek_cl, xk_bz_cl
    USE io_global,     ONLY : stdout, ionode_id
    USE io_var,        ONLY : iufilnscf
    USE io_files,      ONLY : prefix, tmp_dir
    USE input,         ONLY : nkf1, nkf2, nkf3, emax_coulomb, emin_coulomb, &
                              filnscf_coul
    USE ep_constants,  ONLY : ryd2ev, zero, eps6, eps4
    USE mp,            ONLY : mp_bcast
    USE mp_global,     ONLY : inter_pool_comm
    USE mp_world,      ONLY : mpime
    USE low_lvl,       ONLY : mem_size_eliashberg
    !
    IMPLICIT NONE
    !
    CHARACTER(LEN = 256) :: filband
    !! File name
    !CHARACTER(LEN = 256) :: dirname
    !!! Name of the directory where ikmap/egnv/freq/ephmat files are saved
    !
    LOGICAL :: lsda_read
    !! Local variable for lsda
    LOGICAL :: noinv_read
    !! Local variable for noinv
    LOGICAL :: time_reversal_read
    !! Local variable for time_reversal
    LOGICAL :: noncolin_read
    !! Local variable for noncolin
    LOGICAL :: lspinorb_read
    !! Local variable for lspinorb
    LOGICAL :: ldum
    !! false if at given by the file is not consistent
    !
    INTEGER(8) :: imelt
    !! Required allocation of memory
    INTEGER :: ios
    !! IO error message
    INTEGER :: ierr
    !! Error status
    INTEGER :: nsym_read
    !! Local variable for nsym
    INTEGER :: nrot_read
    !! Local variable for nrot
    INTEGER :: nspin_read
    !! Local variable for nspin
    INTEGER :: t_rev_read(48)
    !! Local variable for t_rev
    INTEGER :: s_read(3,3,48)
    !! Local variable for s
    INTEGER :: k1, k2, k3
    !! the grid offsets given by the file
    INTEGER :: i, j, k
    !! Counters
    INTEGER :: ik
    !! Counter on k-point
    INTEGER :: ibnd
    !! Counter on bands
    INTEGER :: ibndmin
    !! Lower band bound
    INTEGER :: ibndmax
    !! Upper band bound
    INTEGER :: isym
    !! Counter on symmetry
    INTEGER :: ik_dum
    !! Dummy used for counter read from the file
    INTEGER :: nbnd_dum
    !! Dummy used for number of bands
    !
    REAL(KIND = DP) :: alat_read
    !! Local variable for alat
    REAL(KIND = DP) :: at_read(3, 3)
    !! Local variable for at
    REAL(KIND = DP) :: at_diff(3, 3)
    !! Difference between at and at_read
    REAL(KIND = DP) :: ebnd
    !! Eigenvalue at (ibnd, ik)
    REAL(KIND = DP), ALLOCATABLE :: ek_dum(:, :)
    !! eigenvalues at E_i(k), ek_dum(nbnd_dum, nkstot_cl)
    !
    WRITE(stdout, '(/5x, a/)') 'Start reading nscf file for Coulomb'
    !
    IF (filnscf_coul == ' ') THEN
      filband = TRIM(tmp_dir) // TRIM(prefix) // '.bands.out'
    ELSE
      filband = TRIM(filnscf_coul)
    ENDIF
    !
    IF (mpime == ionode_id) THEN
      OPEN(UNIT = iufilnscf, FILE = filband, STATUS = 'unknown', FORM ='formatted', IOSTAT = ios)
      IF (ios /= 0) CALL errore('read_eigen_cl', 'error opening file ' // filband, iufilnscf)
      READ(iufilnscf, '(ES23.15)') alat_read
      READ(iufilnscf, '(ES23.15,2(1x,ES23.15))') at_read(:,1)
      READ(iufilnscf, '(ES23.15,2(1x,ES23.15))') at_read(:,2)
      READ(iufilnscf, '(ES23.15,2(1x,ES23.15))') at_read(:,3)
      READ(iufilnscf, '(I5,1x,I5)') nrot_read, nsym_read
      READ(iufilnscf, '(L1)') time_reversal_read
      READ(iufilnscf, '(I2,15(1x,I2))') t_rev_read(1:16)
      READ(iufilnscf, '(I2,15(1x,I2))') t_rev_read(17:32)
      READ(iufilnscf, '(I2,15(1x,I2))') t_rev_read(33:48)
      DO isym = 1, 48
        READ(iufilnscf, '(I3,8(1x,I3))') ((s_read(i, j, isym), i=1,3), j=1,3)
      ENDDO
      READ(iufilnscf, '(L1,1x,I3)') lsda_read, nspin_read
      READ(iufilnscf, '(L1)') noinv_read
      READ(iufilnscf, '(L1,1x,L1)') noncolin_read, lspinorb_read
      READ(iufilnscf, '(I6,2(1x,I6))') nk1_cl, nk2_cl, nk3_cl
      READ(iufilnscf, '(I6,2(1x,I6))') k1, k2, k3
      READ(iufilnscf, '(I6,1x,I10)') nbnd_dum, nkstot_cl
    ENDIF
    !
    !CALL mp_bcast(alat_read, ionode_id, inter_pool_comm)
    !CALL mp_bcast(at_read, ionode_id, inter_pool_comm)
    !CALL mp_bcast(nrot_read, ionode_id, inter_pool_comm)
    !CALL mp_bcast(nsym_read, ionode_id, inter_pool_comm)
    !CALL mp_bcast(time_reversal_read, ionode_id, inter_pool_comm)
    !CALL mp_bcast(t_rev_read, ionode_id, inter_pool_comm)
    !CALL mp_bcast(s_read, ionode_id, inter_pool_comm)
    !CALL mp_bcast(nsym_read, ionode_id, inter_pool_comm)
    !CALL mp_bcast(lsda_read, ionode_id, inter_pool_comm)
    !CALL mp_bcast(nspin_read, ionode_id, inter_pool_comm)
    !CALL mp_bcast(noinv_read, ionode_id, inter_pool_comm)
    !CALL mp_bcast(noncolin_read, ionode_id, inter_pool_comm)
    !CALL mp_bcast(lspinorb_read, ionode_id, inter_pool_comm)
    !CALL mp_bcast(k1, ionode_id, inter_pool_comm)
    !CALL mp_bcast(k2, ionode_id, inter_pool_comm)
    !CALL mp_bcast(k3, ionode_id, inter_pool_comm)
    CALL mp_bcast(nk1_cl, ionode_id, inter_pool_comm)
    CALL mp_bcast(nk2_cl, ionode_id, inter_pool_comm)
    CALL mp_bcast(nk3_cl, ionode_id, inter_pool_comm)
    CALL mp_bcast(nbnd_dum, ionode_id, inter_pool_comm)
    CALL mp_bcast(nkstot_cl, ionode_id, inter_pool_comm)
    !
    IF (mpime == ionode_id) THEN
      IF (ABS(alat_read - alat) > eps6) THEN
        !WRITE(stdout, '(/5x, a, F15.9/)') 'alat_read = ', alat_read
        !WRITE(stdout, '(/5x, a, F15.9/)') 'alat      = ', alat
        WRITE(stdout, '(/5x, a/)') 'WARNING: alat from nscf file for Coulomb is not consistent.'
      ENDIF
      at_diff = at_read - at
      ldum = .TRUE.
      DO i = 1, 3
        DO j = 1, 3
          IF (ABS(at_diff(i, j)) > eps6) THEN
            !WRITE(stdout, '(/5x, a, F15.9/)') 'at_read = ', at_read(i, j)
            !WRITE(stdout, '(/5x, a, F15.9/)') 'at      = ', at(i, j)
            ldum = .FALSE.
          ENDIF
        ENDDO
      ENDDO
      !
      IF (.NOT. ldum) THEN
        WRITE(stdout, '(/5x, a/)') 'WARNING: at from nscf file for Coulomb is not consistent.'
      ENDIF
      !
      ldum = .TRUE.
      DO isym = 1, nrot
        IF (t_rev_read(isym) .NE. t_rev(isym)) THEN
          ldum = .FALSE.
        ENDIF
      ENDDO
      !
      IF (.NOT. ldum) THEN
        WRITE(stdout, '(/5x, a/)') 'WARNING: t_rev from nscf file for Coulomb is not consistent.'
      ENDIF
      !
      ldum = .TRUE.
      DO isym = 1, nrot
        DO j = 1, 3
          DO i = 1, 3
            IF (s_read(i, j, isym) .NE. s(i, j, isym)) THEN
              ldum = .FALSE.
            ENDIF
          ENDDO
        ENDDO
      ENDDO
      !
      IF (.NOT. ldum) THEN
        WRITE(stdout, '(/5x, a/)') 'WARNING: s from nscf file for Coulomb is not consistent.'
      ENDIF
      !
      IF (nsym_read .NE. nsym) THEN
        WRITE(stdout, '(/5x, a/)') 'WARNING: nsym from nscf file for Coulomb is not consistent.'
      ENDIF
      IF (nrot_read .NE. nrot) THEN
        WRITE(stdout, '(/5x, a/)') 'WARNING: nrot from nscf file for Coulomb is not consistent.'
      ENDIF
      IF (lsda_read .NEQV. lsda) THEN
        WRITE(stdout, '(/5x, a/)') 'WARNING: lsda from nscf file for Coulomb is not consistent.'
      ENDIF
      !IF (nspin_read .NE. nspin) THEN
      !  WRITE(stdout, '(/5x, a/)') 'WARNING: nspin from nscf file for Coulomb is not consistent.'
      !ENDIF
      IF (noinv_read .NEQV. noinv) THEN
        WRITE(stdout, '(/5x, a/)') 'WARNING: noinv from nscf file for Coulomb is not consistent.'
      ENDIF
      IF (noncolin_read .NEQV. noncolin) THEN
        WRITE(stdout, '(/5x, a/)') 'WARNING: noncolin from nscf file for Coulomb is not consistent.'
      ENDIF
      IF (lspinorb_read .NEQV. lspinorb) THEN
        WRITE(stdout, '(/5x, a/)') 'WARNING: lspinorb from nscf file for Coulomb is not consistent.'
      ENDIF
      IF (.NOT.((k1 == 0).AND.(k2 == 0).AND.(k3 == 0))) THEN
        CALL errore ('read_eigen_cl', 'k grid from nscf must be zero-centered.', 1)
      ENDIF
    ENDIF ! mpime == ionode_id
    !
    IF (.NOT.((MOD(nkf1, nk1_cl) == 0) .AND. &
              (MOD(nkf2, nk2_cl) == 0) .AND. &
              (MOD(nkf3, nk3_cl) == 0))) THEN
      CALL errore ('read_eigen_cl', &
                   'eliashberg requires nk1, nk2, nk3 to be divisors of nkf1, nkf2, nkf3', 1)
      !
    ENDIF
    !
    ! get the size of the weights and the k vectors that need to be stored in each pool
    imelt = (1 + 3) * nkstot_cl
    CALL mem_size_eliashberg(2, imelt)
    ALLOCATE(wk_cl(nkstot_cl), STAT = ierr)
    IF (ierr /= 0) CALL errore('read_eigen_cl', 'Error allocating wk_cl', 1)
    ALLOCATE(xk_cl(3, nkstot_cl), STAT = ierr)
    IF (ierr /= 0) CALL errore('read_eigen_cl', 'Error allocating xk_cl', 1)
    ! get the size of ek_dum
    imelt = nbnd_dum * nkstot_cl
    CALL mem_size_eliashberg(2, imelt)
    ALLOCATE(ek_dum(nbnd_dum, nkstot_cl), STAT = ierr)
    IF (ierr /= 0) CALL errore('read_eigen_cl', 'Error allocating ek_dum', 1)
    wk_cl(:) = zero
    xk_cl(:, :) = zero
    ek_dum(:, :) = zero
    !
    IF (mpime == ionode_id) THEN
      DO ik = 1, nkstot_cl
        READ(iufilnscf, '(I10,1x,ES23.15)') ik_dum, wk_cl(ik)
        IF (ik_dum .NE. ik) CALL errore ('read_eigen_cl', 'Wrong ik', ik_dum)
        READ(iufilnscf, '(ES23.15,2(1x,ES23.15))') xk_cl(:, ik)
        READ(iufilnscf, '(4(ES24.15))') ek_dum(:, ik)
      ENDDO
      !
      ! go from Ryd to eV
      ek_dum(:, :) = ek_dum(:, :) * ryd2ev
      !  bring all the k points from cartesian to crystal coordinates using at
      CALL cryst_to_cart(nkstot_cl, xk_cl, at, -1)
      !
      IF (noncolin_read .OR. lsda_read) THEN
        IF (ABS(SUM(wk_cl) - 1.d0) > eps4) &
          WRITE(stdout,'(5x,"WARNING: k-point weigths do not add up to 1 [read_eigen_cl]")')
      ELSE
        IF (ABS(SUM(wk_cl) - 2.d0) > eps4) &
          WRITE(stdout,'(5x,"WARNING: k-point weigths do not add up to 2 [read_eigen_cl]")')
      ENDIF
      !
    ENDIF
    !
    CALL mp_bcast(wk_cl, ionode_id, inter_pool_comm)
    CALL mp_bcast(xk_cl, ionode_id, inter_pool_comm)
    !
    WRITE(stdout, '(/5x, a/)') 'Finish reading nscf file for Coulomb'
    !
    IF (mpime == ionode_id) THEN
      CLOSE(iufilnscf)
    ENDIF
    !
    IF (mpime == ionode_id) THEN
      ! Note: Only IO node has the data of ek_dum.
      WRITE(stdout, '(/5x, 3a, 3I4)') 'k-grid read from ', TRIM(filband), ' : ', nk1_cl, nk2_cl, nk3_cl
      WRITE(stdout, '(5x, 3a, I9)') 'Nr irreducible k-points read from ', TRIM(filband), ' : ', nkstot_cl
      WRITE(stdout, '(5x, a, ES20.10)') 'Minimum eigenvalue of bands taken from the file (eV) = ', &
                                        MINVAL(ek_dum(:, :))
      WRITE(stdout, '(5x, a, ES20.10)') 'Maximum eigenvalue of bands taken from the file (eV) = ', &
                                        MAXVAL(ek_dum(:, :))
      WRITE(stdout, '(5x, a, ES20.10)') 'emin_coulomb + "Fermi level" (eV)                    = ', &
                                        emin_coulomb + ef0
      WRITE(stdout, '(5x, a, ES20.10)') 'emax_coulomb + "Fermi level" (eV)                    = ', &
                                        emax_coulomb + ef0
      IF((MAXVAL(ek_dum(:, :)) >= (emin_coulomb + ef0)) .OR. &
         (MINVAL(ek_dum(:, :)) <= (emin_coulomb + ef0))) THEN
        WRITE(stdout, '(5x,a,f14.6,a,f14.6,a)' ) 'Only states taken from nscf file between ', &
                       (emin_coulomb + ef0), ' eV and ', &
                       (emax_coulomb + ef0), ' eV' 
        WRITE(stdout, '(5x,a)' ) 'will be included in the Eliashberg calculations.'
      ENDIF
      !
      ibndmin = 100000
      ibndmax = 0
      !
      DO ik = 1, nkstot_cl
        DO ibnd = 1, nbnd_dum
          ebnd = ek_dum(ibnd, ik)
          !
          IF (((ebnd - ef0) < emax_coulomb) .AND. &
              ((ebnd - ef0) > emin_coulomb)) THEN
            ibndmin = MIN(ibnd, ibndmin)
            ibndmax = MAX(ibnd, ibndmax)
          ENDIF
          !
        ENDDO
      ENDDO
      !
      nbnd_cl = ibndmax - ibndmin + 1
    ENDIF
    CALL mp_bcast(nbnd_cl, ionode_id, inter_pool_comm)
    !
    WRITE(stdout,'(5x, i7, a23, a60/)') nbnd_cl, ' bands in the interval ', &
                 '[emin_coulomb + "Fermi level", emax_coulomb + "Fermi level"]'
    !
    ! get the size of ek_cl
    imelt = nbnd_cl * nkstot_cl
    CALL mem_size_eliashberg(2, imelt)
    ALLOCATE(ek_cl(nbnd_cl, nkstot_cl), STAT = ierr)
    IF (ierr /= 0) CALL errore('read_eigen_cl', 'Error allocating ek_cl', 1)
    !
    IF (mpime == ionode_id) THEN
      ! Copy the values of ek_dum into ek_cl
      ! The values are already in the unit of eV.
      DO ik = 1, nkstot_cl
        ek_cl(1:nbnd_cl, ik) = ek_dum(ibndmin:ibndmax, ik)
      ENDDO
    ENDIF
    !
    CALL mp_bcast(ek_cl, ionode_id, inter_pool_comm)
    !
    !
    ALLOCATE(xk_bz_cl(3, nk1_cl * nk2_cl * nk3_cl), STAT = ierr)
    IF (ierr /= 0) CALL errore('read_eigen_cl', 'Error allocating xk_bz_cl', 1)
    IF (mpime == ionode_id) THEN
      DO i = 1, nk1_cl
        DO j = 1, nk2_cl
          DO k = 1, nk3_cl
            ik = (i - 1) * nk2_cl * nk3_cl + (j - 1) * nk3_cl + k
            xk_bz_cl(1, ik) = DBLE(i - 1) / DBLE(nk1_cl)
            xk_bz_cl(2, ik) = DBLE(j - 1) / DBLE(nk2_cl)
            xk_bz_cl(3, ik) = DBLE(k - 1) / DBLE(nk3_cl)
            xk_bz_cl(1, ik) = xk_bz_cl(1, ik) - NINT(xk_bz_cl(1, ik))
            xk_bz_cl(2, ik) = xk_bz_cl(2, ik) - NINT(xk_bz_cl(2, ik))
            xk_bz_cl(3, ik) = xk_bz_cl(3, ik) - NINT(xk_bz_cl(3, ik))
          ENDDO
        ENDDO
      ENDDO
    ENDIF
    CALL mp_bcast(xk_bz_cl, ionode_id, inter_pool_comm)
    !
    RETURN
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE read_eigen_cl
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    SUBROUTINE find_indices_ik_cl()
    !-----------------------------------------------------------------------
    !!
    !! this routine finds k points, which are given by nscf file for Coulomb 
    !! output by bands2epw.x on the fine k-mesh and then obtains conversion 
    !! function for band index, "ik_cl_to_fs" and "ik_fs_to_cl".
    !!
    !! ik_cl_to_fs(i) returns 0 if xk_cl(i) does not exist within 
    !! the Fermi shell, and it returns the corresponding index value otherwise.
    !! That means
    !! xk_cl(:, i) = xkfs(:, ik_cl_to_fs(i))
    !! if ik_cl_to_fs(i) is not zero.
    !!
    !! ik_fs_to_cl(j) returns 0 if xkfs(j) does not exist on 
    !! the coase k grid given by the file "prefix.bands.out", and it returns 
    !! the corresponding index value otherwise.
    !! That means:
    !! xkfs(:, j) = xk_cl(:, ik_fs_to_cl(j))
    !! if ik_fs_to_cl(j) is not zero.
    !!
    USE kinds,         ONLY : DP
    USE io_global,     ONLY : stdout
    USE input,         ONLY : nkf1, nkf2, nkf3
    USE supercond_common, ONLY : nkfs, xkfs, ik_cl_to_fs, nkstot_cl, xk_cl, &
                              ik_fs_to_cl, nk1_cl, nk2_cl, nk3_cl, &
                              xk_bz_cl, ik_bz_to_ibz_cl, ik_ibz_to_bz_cl
    USE ep_constants,  ONLY : eps6, eps5, eps8
    USE mp,            ONLY : mp_sum
    USE low_lvl,       ONLY : mem_size_eliashberg
    USE kfold,         ONLY : backtoBZ
    USE symm_base,     ONLY : nsym, t_rev, s
    !
    IMPLICIT NONE
    !
    ! Local variables
    LOGICAL :: in_the_list
    !! Is the point in the list
    INTEGER :: ik, ipol
    !! Counter on k points
    INTEGER :: ik_cl
    !! Counter on k points given by the file
    INTEGER :: i, j, k
    !! Counters
    INTEGER :: n
    !! Current index of k point
    INTEGER :: ns
    !! Counter on symmetry
    !
    INTEGER(8) :: imelt
    !! Required allocation of memory
    INTEGER :: ierr
    !! Error status
    REAL(KIND = DP):: xk_dum(3)
    !! Dummy
    !
    REAL(KIND = DP) :: xx1, yy1, zz1
    !! Temporary variables
    REAL(KIND = DP) :: xx2, yy2, zz2
    !! Temporary variables
    !
    imelt = nkstot_cl
    CALL mem_size_eliashberg(1, imelt)
    ALLOCATE(ik_cl_to_fs(nkstot_cl), STAT = ierr)
    IF (ierr /= 0) CALL errore('find_indices_ik_cl', 'Error allocating ik_cl_to_fs', 1)
    ik_cl_to_fs(:) = 0
    !
    imelt = nkfs
    CALL mem_size_eliashberg(1, imelt)
    ALLOCATE(ik_fs_to_cl(nkfs), STAT = ierr)
    IF (ierr /= 0) CALL errore('find_indices_ik_cl', 'Error allocating ik_fs_to_cl', 1)
    ik_fs_to_cl(:) = 0
    !
    ! xk_cl and xkf are assumed to be in crys coord
    !
    !DO ik_cl = 1, nkstot_cl
    !  xx1 = xk_cl(1, ik_cl) * nkf1
    !  yy1 = xk_cl(2, ik_cl) * nkf2
    !  zz1 = xk_cl(3, ik_cl) * nkf3
    !  in_the_list = ABS(xx1 - NINT(xx1)) <= eps6 .AND. &
    !                ABS(yy1 - NINT(yy1)) <= eps6 .AND. &
    !                ABS(zz1 - NINT(zz1)) <= eps6
    !  IF (.NOT. in_the_list) &
    !    CALL errore('find_indices_ik_cl', 'nscf k point does not fall on fine k-grid', 1)
    !  !
    !  CALL backtoBZ(xx1, yy1, zz1, nkf1, nkf2, nkf3)
    !  !
    !  DO ik = 1, nkfs
    !    xx2 = xkfs(1, ik) * nkf1
    !    yy2 = xkfs(2, ik) * nkf2
    !    zz2 = xkfs(3, ik) * nkf3
    !    ! no need to check whether xkfs fall on fine k-grid
    !    ! NOTE: All pools have the same values for xkfs and nkfs.
    !    !
    !    CALL backtoBZ(xx2, yy2, zz2, nkf1, nkf2, nkf3)
    !    !
    !    IF ((NINT(xx1) == NINT(xx2)) .AND. &
    !        (NINT(yy1) == NINT(yy2)) .AND. &
    !        (NINT(zz1) == NINT(zz2))) THEN
    !      !
    !      IF (ik_cl_to_fs(ik_cl) .NE. 0) & 
    !        CALL errore('find_indices_ik_cl', 'ik_cl_to_fs already have a non-zero value', 1)
    !      !
    !      ik_cl_to_fs(ik_cl) = ik
    !      !
    !              !
    !      IF (ik_fs_to_cl(ik) .NE. 0) & 
    !        CALL errore('find_indices_ik_cl', 'ik_fs_to_cl already have a non-zero value', 2)
    !      !
    !      ik_fs_to_cl(ik) = ik_cl 
    !      !
    !    ENDIF
    !    !
    !    IF ((ABS(xx1 - xx2) < eps5) .AND. &
    !        (ABS(yy1 - yy2) < eps5) .AND. &
    !        (ABS(zz1 - zz2) < eps5)) EXIT
    !    !
    !  ENDDO
    !  !
    !ENDDO
    DO ik = 1, nkfs
      xx1 = xkfs(1, ik) * nkf1
      yy1 = xkfs(2, ik) * nkf2
      zz1 = xkfs(3, ik) * nkf3
      ! no need to check whether xkfs fall on fine k-grid
      ! NOTE: All pools have the same values for xkfs and nkfs.
      !
      CALL backtoBZ(xx1, yy1, zz1, nkf1, nkf2, nkf3)
      !
      DO ik_cl = 1, nkstot_cl
        xx2 = xk_cl(1, ik_cl) * nkf1
        yy2 = xk_cl(2, ik_cl) * nkf2
        zz2 = xk_cl(3, ik_cl) * nkf3
        in_the_list = ABS(xx2 - NINT(xx2)) <= eps6 .AND. &
                      ABS(yy2 - NINT(yy2)) <= eps6 .AND. &
                      ABS(zz2 - NINT(zz2)) <= eps6
        IF (.NOT. in_the_list) &
          CALL errore('find_indices_ik_cl', 'nscf k point does not fall on fine k-grid', 1)
        !
        CALL backtoBZ(xx2, yy2, zz2, nkf1, nkf2, nkf3)
        !
        !IF ((ABS(xx1 - xx2) < eps5) .AND. &
        !    (ABS(yy1 - yy2) < eps5) .AND. &
        !    (ABS(zz1 - zz2) < eps5)) THEN
        IF ((NINT(xx1) == NINT(xx2)) .AND. &
            (NINT(yy1) == NINT(yy2)) .AND. &
            (NINT(zz1) == NINT(zz2))) THEN
          !
          IF (ik_cl_to_fs(ik_cl) .NE. 0) & 
            CALL errore('find_indices_ik_cl', 'ik_cl_to_fs already have a non-zero value', 1)
          !
          ik_cl_to_fs(ik_cl) = ik
          !
                  !
          IF (ik_fs_to_cl(ik) .NE. 0) & 
            CALL errore('find_indices_ik_cl', 'ik_fs_to_cl already have a non-zero value', 2)
          !
          ik_fs_to_cl(ik) = ik_cl 
          !
        ENDIF
        !
        IF ((ABS(xx1 - xx2) < eps5) .AND. &
            (ABS(yy1 - yy2) < eps5) .AND. &
            (ABS(zz1 - zz2) < eps5)) EXIT
        !
      ENDDO
      !
    ENDDO
    !
    imelt = nk1_cl * nk2_cl * nk3_cl
    CALL mem_size_eliashberg(1, imelt)
    ALLOCATE(ik_bz_to_ibz_cl(nk1_cl * nk2_cl * nk3_cl), STAT = ierr)
    IF (ierr /= 0) CALL errore('find_indices_ik_cl', 'Error allocating ik_bz_to_ibz_cl', 1)
    ik_bz_to_ibz_cl(:) = 0
    !
    imelt = nkstot_cl
    CALL mem_size_eliashberg(1, imelt)
    ALLOCATE(ik_ibz_to_bz_cl(nkstot_cl), STAT = ierr)
    IF (ierr /= 0) CALL errore('find_indices_ik_cl', 'Error allocating ik_ibz_to_bz_cl', 1)
    ik_ibz_to_bz_cl(:) = 0
    !
    DO ik_cl = 1, nkstot_cl
      DO ns = 1, nsym
        DO i = 1, 3
          xk_dum(i) = s(i, 1, ns) * xk_cl(1, ik_cl) &
                    + s(i, 2, ns) * xk_cl(2, ik_cl) &
                    + s(i, 3, ns) * xk_cl(3, ik_cl)
          xk_dum(i) = xk_dum(i) - NINT(xk_dum(i) + eps8) 
        ENDDO
        IF(t_rev(ns) == 1) xk_dum = -xk_dum
        xx2 = xk_dum(1) * nk1_cl
        yy2 = xk_dum(2) * nk2_cl
        zz2 = xk_dum(3) * nk3_cl
        !
        i = MOD(NINT(xk_dum(1) * nk1_cl + 2 * nk1_cl), nk1_cl) + 1
        j = MOD(NINT(xk_dum(2) * nk2_cl + 2 * nk2_cl), nk2_cl) + 1
        k = MOD(NINT(xk_dum(3) * nk3_cl + 2 * nk3_cl), nk3_cl) + 1
        n = (k - 1) + (j - 1) * nk3_cl + (i - 1) * nk2_cl * nk3_cl + 1
        !
        IF ((ABS(xk_dum(1) - xk_bz_cl(1, n)) > eps5) .OR. &
            (ABS(xk_dum(2) - xk_bz_cl(2, n)) > eps5) .OR. &
            (ABS(xk_dum(3) - xk_bz_cl(3, n)) > eps5)) THEN
              WRITE(stdout, '(8x,"xk_cl (",i5,") = (",3f12.7,")")') ik_cl, &
                   (xk_cl(ipol, ik_cl) , ipol = 1, 3)
              WRITE(stdout, '(8x,"S * xk_cl     = (",3f12.7,")")') &
                   (xk_dum(ipol) , ipol = 1, 3)
              WRITE(stdout, '(8x,"xk_bz (",i5,") = (",3f12.7,")")') ik, &
                   (xk_bz_cl(ipol, n) , ipol = 1, 3)
          CALL errore('find_indices_ik_cl', 'xk_cl does not match.', 1)
        ENDIF
        !
        ik_bz_to_ibz_cl(n) = ik_cl
        IF (ik_ibz_to_bz_cl(ik_cl) == 0) &
          ik_ibz_to_bz_cl(ik_cl) = n
        !
      ENDDO ! ns
      !
    ENDDO
    !
    DO ik = 1, nk1_cl * nk2_cl * nk3_cl
      IF (ik_bz_to_ibz_cl(ik) == 0) &
        CALL errore('find_indices_ik_cl', 'There is no k point equivalent to some xk_cl in the 1st BZ.', ik)
    ENDDO
    !
    DO ik_cl = 1, nkstot_cl
      IF (ik_ibz_to_bz_cl(ik_cl) == 0) &
      CALL errore('find_indices_ik_cl', 'There is no k point equivalent to some xk_cl in the 1st BZ (2).', ik)
    ENDDO
    !
    !!!!!!!! FOR DEBUG !!!!!!!!
    !WRITE(stdout, '(/)')
    !DO ik = 1, nk1_cl * nk2_cl * nk3_cl
    !  WRITE(stdout, '(8x,"xk_bz (",i5,") = (",3f12.7,")")') ik, &
    !  (xk_bz_cl(ipol, ik) , ipol = 1, 3)
    !ENDDO
    !WRITE(stdout, '(/)')
    !DO ik_cl = 1, nkstot_cl
    !  WRITE(stdout, '(8x,"xk_bz (",i5,") = (",3f12.7,")")') ik_ibz_to_bz_cl(ik_cl), &
    !  (xk_bz_cl(ipol, ik_ibz_to_bz_cl(ik_cl)) , ipol = 1, 3)
    !ENDDO
    !WRITE(stdout, '(/)')
    !DO ik_cl = 1, nkstot_cl
    !  WRITE(stdout, '(8x,"xk_cl (",i5,") = (",3f12.7,")")') ik_cl, &
    !  (xk_cl(ipol, ik_cl) , ipol = 1, 3)
    !ENDDO
    !WRITE(stdout, '(/)')
    !ik = 0
    !DO ik_cl = 1, nkstot_cl
    !  IF (ik_cl_to_fs(ik_cl) .NE. 0) THEN
    !    ik = ik + 1
    !    WRITE(stdout, '(8x,"xk_cl(",i5,") = (",3f12.7,")")') ik_cl, &
    !         (xk_cl(ipol, ik_cl) , ipol = 1, 3)
    !    WRITE(stdout, '(8x,"xkfs (",i5,") = (",3f12.7,")")') ik_cl_to_fs(ik_cl), &
    !         (xkfs(ipol, ik_cl_to_fs(ik_cl)) , ipol = 1, 3)
    !  ENDIF
    !ENDDO
    !WRITE(stdout, '(/5x, a, I4/)')    'nkstot_cl_tmp = ', ik
    !ik_cl = 0
    !DO ik = 1, nkfs
    !  IF (ik_fs_to_cl(ik) .NE. 0) THEN
    !    ik_cl = ik_cl + 1
    !    WRITE(stdout, '(8x,"xk_cl(",i5,") = (",3f12.7,")")') ik_fs_to_cl(ik), &
    !         (xk_cl(ipol, ik_fs_to_cl(ik)) , ipol = 1, 3)
    !    WRITE(stdout, '(8x,"xkfs (",i5,") = (",3f12.7,")")') ik, &
    !         (xkfs(ipol, ik) , ipol = 1, 3)
    !  ENDIF
    !ENDDO
    !WRITE(stdout, '(/5x, a, I4/)')    'nkstot_cl_tmp = ', ik_cl
    !!
    !!WRITE(1000 + mpime, '(/)')
    !!DO ik = 1, nk1_cl * nk2_cl * nk3_cl
    !!  WRITE(1000 + mpime, '(8x,"xk_bz (",i5,") = (",3f12.7,")")') ik, &
    !!  (xk_bz_cl(ipol, ik) , ipol = 1, 3)
    !!ENDDO
    !!WRITE(1000 + mpime, '(/)')
    !!DO ik_cl = 1, nkstot_cl
    !!  WRITE(1000 + mpime, '(8x,"xk_bz (",i5,") = (",3f12.7,")")') ik_ibz_to_bz_cl(ik_cl), &
    !!  (xk_bz_cl(ipol, ik_ibz_to_bz_cl(ik_cl)) , ipol = 1, 3)
    !!ENDDO
    !!WRITE(1000 + mpime, '(/)')
    !!DO ik_cl = 1, nkstot_cl
    !!  WRITE(1000 + mpime, '(8x,"xk_cl (",i5,") = (",3f12.7,")")') ik_cl, &
    !!  (xk_cl(ipol, ik_cl) , ipol = 1, 3)
    !!ENDDO
    !!WRITE(1000 + mpime, '(/)')
    !!ik = 0
    !!DO ik_cl = 1, nkstot_cl
    !!  IF (ik_cl_to_fs(ik_cl) .NE. 0) THEN
    !!    ik = ik + 1
    !!    WRITE(1000 + mpime, '(8x,"xk_cl(",i5,") = (",3f12.7,")")') ik_cl, &
    !!         (xk_cl(ipol, ik_cl) , ipol = 1, 3)
    !!    WRITE(1000 + mpime, '(8x,"xkfs (",i5,") = (",3f12.7,")")') ik_cl_to_fs(ik_cl), &
    !!         (xkfs(ipol, ik_cl_to_fs(ik_cl)) , ipol = 1, 3)
    !!  ENDIF
    !!ENDDO
    !!WRITE(1000 + mpime, '(/5x, a, I4/)')    'nkstot_cl_tmp = ', ik
    !!ik_cl = 0
    !!DO ik = 1, nkfs
    !!  IF (ik_fs_to_cl(ik) .NE. 0) THEN
    !!    ik_cl = ik_cl + 1
    !!    WRITE(1000 + mpime, '(8x,"xk_cl(",i5,") = (",3f12.7,")")') ik_fs_to_cl(ik), &
    !!         (xk_cl(ipol, ik_fs_to_cl(ik)) , ipol = 1, 3)
    !!    WRITE(1000 + mpime, '(8x,"xkfs (",i5,") = (",3f12.7,")")') ik, &
    !!         (xkfs(ipol, ik) , ipol = 1, 3)
    !!  ENDIF
    !!ENDDO
    !!WRITE(1000 + mpime, '(/5x, a, I4/)')    'nkstot_cl_tmp = ', ik_cl
    !!!!!!!! FOR DEBUG !!!!!!!!
    !
    RETURN
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE find_indices_ik_cl
    !-----------------------------------------------------------------------
    !
    !-----------------------------------------------------------------------
    SUBROUTINE find_nbnd_offset()
    !-----------------------------------------------------------------------
    !!
    !! This routine obtains the conversion function for band index
    !! from ek_cl to ekf, and vice versa, "nbnd_offset".
    !! Using "nbnd_offset", we can convert the band index from ek_cl 
    !! to ekfs, and vice versa. That is,
    !!   ek_cl(ibnd_cl, ik_fs_to_cl(ik)) ~= ekfs(ibnd, ik)
    !!   ekfs(ibndfs, ik_cl_to_fs(ik_cl)) ~= ek_cl(ibnd, ik_cl)
    !! where
    !!   ibnd_cl = ibnd_kfs_to_kfs_all(ibnd, ik) + nbnd_offset
    !!   ibndfs = ibnd_kfs_all_to_kfs(ibnd - nbnd_offset, ik_cl_to_fs(ik_cl))
    !! within the Fermi shell.
    !!
    !! NOTE:
    !! Because the data stored in ekfs is transferred from etk in a
    !! complicated way, it is not straightforward to suppose the original
    !! band index of etk from ekfs.
    !! See also read_eigenvalues in io_eliashberg.
    !!
    !! What this routine do are:
    !! 1. Find the index numbers of the ek_cl and ekf with the closest 
    !!    energy to ef0, "nbndmin_cl" and "nbndmin_fs". 
    !!    Because ekf is already unallocated here, ekfs is used instead.
    !!    Please note that 
    !!    ekfs(ibnd, ik) = ekf(ibnd2, ik2)
    !!    with 
    !!    ibnd2 = ibnd_kfs_to_kfs_all(ibnd, ik),
    !!    ik2 = ikfs_to_kf(ik).
    !! 2. Find the function "nbnd_offset", which converts index of ekf to 
    !!    index of ek_cl, by assuming that nbnd_offset is around 
    !!    "nbndmin_cl - nbndmin_fs".
    !! 3. Check whether all the values of ek_cl coincide with those of
    !!    ekfs using "nbnd_offset".
    !!
    USE kinds,         ONLY : DP
    USE io_global,     ONLY : stdout
    USE input,        ONLY : fsthick
    USE supercond_common, ONLY : ekfs, nkfs, nbndfs, ef0, &
                              ek_cl, nbnd_cl, nkstot_cl, &
                              ik_cl_to_fs, ik_fs_to_cl, &
                              ibnd_kfs_to_kfs_all, ibnd_kfs_all_to_kfs, &
                              nbnd_offset, nbndfs_all
    USE ep_constants,  ONLY : zero
    USE mp,            ONLY : mp_sum
    USE kfold,         ONLY : backtoBZ
    !
    IMPLICIT NONE
    !
    ! Local variables
    INTEGER :: nbndfs_cl
    !! The maximum number of bands within the Fermi shell
    INTEGER :: nbndfs_dum
    !! The maximum number of bands within the Fermi shell
    INTEGER :: n1
    !! Counter on the number of bands within the Fermi shell for each k point
    INTEGER :: n2
    !! Counter on the number of bands within the Fermi shell for each k point
    INTEGER :: nbndmin_cl
    !! The band index of ek_cl which is the closest to ef0
    INTEGER :: nbndmin_fs
    !! The band index of ekfs which is the closest to ef0
    INTEGER :: ik_cl
    !! Counter on k points for ek_cl
    INTEGER :: ik
    !! Counter on k points
    INTEGER :: ibnd
    !! Counter on bands
    INTEGER :: ipol
    !! Counter for ediff
    INTEGER :: ibnd_cl
    !! Counter on bands for ek_cl
    INTEGER :: ibndfs
    !! Counter on bands for ekfs
    REAL(KIND = DP) :: mindiff_cl
    !! To store the minimum value of ABS(ek_cl - ef0)
    REAL(KIND = DP) :: mindiff_fs
    !! To store the minimum value of ABS(ekfs - ef0)
    REAL(KIND = DP) :: ebnd
    !! Local variable for energy
    REAL(KIND = DP) :: ediff(-5:5)
    !! To store the sum of absolute differences, ABS(ekfs - ek_cl)
    REAL(KIND = DP) :: ediffmin
    !! used to find the minimum value of ediff
    REAL(KIND = DP) :: ebndmax1
    !! temporarily store the maximum value of ABS(ekfs-ek_cl)
    REAL(KIND = DP) :: ebndmax2
    !! temporarily store the maximum value of ABS(ekfs-ek_cl)
    REAL(KIND = DP) :: ebndmax3
    !! temporarily store the maximum value of ABS(ekfs-ek_cl)
    REAL(KIND = DP) :: ebndmax4
    !! temporarily store the maximum value of ABS(ekfs-ek_cl)
    !
    !
    !! 1. Find the index numbers of the ek_cl and ekf with the closest 
    !!    energy to ef0, "nbndmin_cl" and "nbndmin_fs". 
    !!    Because ekf is already unallocated here, ekfs is used instead.
    !!    Please note that 
    !!    ekfs(ibnd, ik) = ekf(ibnd2, ik2)
    !!    with 
    !!    ibnd2 = ibnd_kfs_to_kfs_all(ibnd, ik),
    !!    ik2 = ikfs_to_kf(ik).
    !
    nbndfs_cl = 0
    n1 = 0
    mindiff_cl = 1D5
    nbndmin_cl = 0
    nbndfs_dum = 0
    n2 = 0
    mindiff_fs = 1D5
    nbndmin_fs = 0
    !
    !!!!!!!! FOR DEBUG !!!!!!!!
    !DO ibnd = 1, nbndfs
    !  WRITE(stdout, '(/5x, a, f15.7/)')    'ekfs  = ', ekfs(ibnd, 2)
    !ENDDO
    !DO ibnd = 1, nbnd_cl
    !  WRITE(stdout, '(/5x, a, f15.7/)')    'ek_cl = ', ek_cl(ibnd, 2)
    !ENDDO
    !!!!!!!! FOR DEBUG !!!!!!!!
    DO ik_cl = 1, nkstot_cl
      IF (ik_cl_to_fs(ik_cl) == 0) CYCLE
      n1 = 0
      n2 = 0
      DO ibnd = 1, nbnd_cl
        ebnd = ek_cl(ibnd, ik_cl)
        IF (ABS(ebnd - ef0) < fsthick) THEN
          n1 = n1 + 1
          IF (nbndfs_cl < n1) nbndfs_cl = n1
          IF (ABS(ebnd - ef0) < mindiff_cl) THEN
            mindiff_cl = ABS(ebnd - ef0)
            nbndmin_cl = ibnd
          ENDIF 
        ENDIF
      ENDDO
      DO ibnd = 1, nbndfs
        ebnd = ekfs(ibnd, ik_cl_to_fs(ik_cl))
        IF (ABS(ebnd - ef0) < fsthick) THEN
          n2 = n2 + 1
          IF (nbndfs_dum < n2) nbndfs_dum = n2
          IF (ABS(ebnd - ef0) < mindiff_fs) THEN
            mindiff_fs = ABS(ebnd - ef0)
            nbndmin_fs = ibnd_kfs_to_kfs_all(ibnd, ik_cl_to_fs(ik_cl))
          ENDIF 
        ENDIF
      ENDDO
      !
      !IF (n1 .NE. n2) CALL errore('find_ibnd_cl', &
      !'The number of eigenvalues of ek_cl does not coincide with that of ekfs within the Fermi shell.', 1)
      !
      IF (nbndmin_fs == 0) CALL errore('find_ibnd_cl', &
      'nbndmin_fs must be larger than zero.', 1)
      !
      ! nbndfs_cl and nbndfs need not match.
      !
    ENDDO
    !
    !! 2. Find the function "nbnd_offset", which converts index of ekf to 
    !!    index of ek_cl, by assuming that nbnd_offset is around 
    !!    "nbndmin_cl - nbndmin_fs".
    !
    ! 
    ediff(:) = zero
    ediffmin = 1d5
    !
    ! Here we put the trial nbnd_offset: 
    ! nbnd_offset = ipol + nbndmin_cl - nbndmin_fs.
    ! While changing ipol, take a summation of the absolute values of 
    ! the energy differences, ediff(ipol), over the states in the fsthick window.
    DO ipol = -5, 5
      nbnd_offset = ipol + nbndmin_cl - nbndmin_fs
      IF (ediff(ipol) > 9.99d3) THEN
        CYCLE
      ENDIF
      DO ik_cl = 1, nkstot_cl
        IF (ik_cl_to_fs(ik_cl) == 0) CYCLE
        !DO ibnd = 1, nbndfs_all
        !  IF (ibnd_kfs_all_to_kfs(ibnd, ik_cl_to_fs(ik_cl)) == 0) CYCLE
        !  ibnd_cl = ibnd + nbnd_offset
        !  IF ((ibnd_cl < 1) .OR. (ibnd_cl > nbnd_cl)) THEN
        !    ediff(ipol) = 1d4
        !    CYCLE
        !  ENDIF
        !  ibndfs = ibnd_kfs_all_to_kfs(ibnd, ik_cl_to_fs(ik_cl))
        !  ediff = ediff + &
        !  ABS(ekfs(ibndfs, ik_cl_to_fs(ik_cl)) - ek_cl(ibnd_cl, ik_cl))
        !ENDDO
        DO ibnd = 1, nbndfs
          ibnd_cl = ibnd_kfs_to_kfs_all(ibnd, ik_cl_to_fs(ik_cl)) + nbnd_offset
          IF ((ibnd_cl < 1) .OR. (ibnd_cl > nbnd_cl)) THEN
            !
            ! no need to calculate ediff if ibnd_cl is outside of the bound of ek_cl
            !
            ediff(ipol) = 1d4
            CYCLE
          ENDIF
          IF (ABS(ekfs(ibnd, ik_cl_to_fs(ik_cl)) - ef0) < fsthick) THEN
            ediff(ipol) = ediff(ipol) + &
            ABS(ekfs(ibnd, ik_cl_to_fs(ik_cl)) - ek_cl(ibnd_cl, ik_cl))
          ENDIF
        ENDDO
      ENDDO
    ENDDO
    !
    ! Next, find one ipol making ediff the smallest and determine nbnd_offset from the ipol.
    DO ipol = -5, 5
      !!!!!!!! FOR DEBUG !!!!!!!!
      !WRITE(stdout, '(/5x, a, I4/)')    'ipol            = ', ipol
      !WRITE(stdout, '(/5x, a, F15.9/)') 'ediff(ipol)     = ', ediff(ipol)
      !!!!!!!! FOR DEBUG !!!!!!!!
      IF(ediff(ipol) < ediffmin) THEN
        ediffmin = ediff(ipol)
        nbnd_offset = ipol + nbndmin_cl - nbndmin_fs
      ENDIF
    ENDDO
    !
    !!!!!!!! FOR DEBUG !!!!!!!!
    !WRITE(stdout, '(/5x, a, I4/)')    'nbnd_offset = ', nbnd_offset
    !!!!!!!! FOR DEBUG !!!!!!!!
    !
    !! 3. Check whether all the values of ek_cl coincide with those of
    !!    ekfs using "nbnd_offset".
    !
    ebndmax1 = zero
    ebndmax2 = zero
    ebndmax3 = zero
    ebndmax4 = zero
    !
    DO ik_cl = 1, nkstot_cl
      IF (ik_cl_to_fs(ik_cl) == 0) CYCLE
      DO ibnd = 1, nbndfs
        IF (ABS(ekfs(ibnd, ik_cl_to_fs(ik_cl)) - ef0) < fsthick) THEN
          ibnd_cl = ibnd_kfs_to_kfs_all(ibnd, ik_cl_to_fs(ik_cl)) + nbnd_offset
          IF ((ibnd_cl < 1) .OR. (ibnd_cl > nbnd_cl)) THEN
            CALL errore('find_ibnd_cl', &
                       'ibnd_cl must be in [1, nbnd_cl].', 1)
          ENDIF
          ebnd = ekfs(ibnd, ik_cl_to_fs(ik_cl)) - ek_cl(ibnd_cl, ik_cl)
          ebnd = ABS(ebnd)
          ebndmax1 = MAX(ebnd, ebndmax1)
          !!!!!!!! FOR DEBUG !!!!!!!!
          !WRITE(stdout, '(5x, f15.7, 1x, f15.7)') &
          !ekfs(ibnd, ik_cl_to_fs(ik_cl)), ek_cl(ibnd_cl, ik_cl)
          !!!!!!!! FOR DEBUG !!!!!!!!
        ENDIF
      ENDDO
    ENDDO
    !
    !!!!!!!! FOR DEBUG !!!!!!!!
    !WRITE(stdout, '(/)')
    DO ik_cl = 1, nkstot_cl
      IF (ik_cl_to_fs(ik_cl) == 0) CYCLE
      DO ibnd = 1, nbnd_cl
        IF (((ibnd - nbnd_offset) < 1) .OR. ((ibnd - nbnd_offset) > nbndfs_all)) CYCLE
        ibndfs = ibnd_kfs_all_to_kfs(ibnd - nbnd_offset, ik_cl_to_fs(ik_cl))
        IF ((ibndfs < 1) .OR. (ibndfs > nbndfs)) CYCLE
        IF (ABS(ekfs(ibndfs, ik_cl_to_fs(ik_cl)) - ef0) < fsthick) THEN
          ebnd = ekfs(ibndfs, ik_cl_to_fs(ik_cl)) - ek_cl(ibnd, ik_cl)
          ebnd = ABS(ebnd)
          ebndmax2 = MAX(ebnd, ebndmax2)
          !!!!!!!! FOR DEBUG !!!!!!!!
          !WRITE(stdout, '(5x, f15.7, 1x, f15.7)') &
          !ekfs(ibndfs, ik_cl_to_fs(ik_cl)), ek_cl(ibnd, ik_cl)
          !!!!!!!! FOR DEBUG !!!!!!!!
        ENDIF
      ENDDO
    ENDDO
    !
    !WRITE(stdout, '(/)')
    DO ik = 1, nkfs
      IF (ik_fs_to_cl(ik) == 0) CYCLE
      DO ibnd = 1, nbndfs
        IF (ABS(ekfs(ibnd, ik) - ef0) < fsthick) THEN
          ibnd_cl = ibnd_kfs_to_kfs_all(ibnd, ik) + nbnd_offset
          IF ((ibnd_cl < 1) .OR. (ibnd_cl > nbnd_cl)) THEN
            CALL errore('find_ibnd_cl', &
                       'ibnd_cl must be in [1, nbnd_cl].', 1)
          ENDIF
          ebnd = ekfs(ibnd, ik) - ek_cl(ibnd_cl, ik_fs_to_cl(ik))
          ebnd = ABS(ebnd)
          ebndmax3 = MAX(ebnd, ebndmax3)
          !!!!!!!! FOR DEBUG !!!!!!!!
          !WRITE(stdout, '(5x, f15.7, 1x, f15.7)') &
          !ekfs(ibnd, ik), ek_cl(ibnd_cl, ik_fs_to_cl(ik))
          !!!!!!!! FOR DEBUG !!!!!!!!
        ENDIF
      ENDDO
    ENDDO
    !
    !WRITE(stdout, '(/)')
    DO ik = 1, nkfs
      IF (ik_fs_to_cl(ik) == 0) CYCLE
      DO ibnd = 1, nbnd_cl
        IF (((ibnd - nbnd_offset) < 1) .OR. ((ibnd - nbnd_offset) > nbndfs_all)) CYCLE
        ibndfs = ibnd_kfs_all_to_kfs(ibnd - nbnd_offset, ik)
        IF ((ibndfs < 1) .OR. (ibndfs > nbndfs)) CYCLE
        IF (ABS(ekfs(ibndfs, ik) - ef0) < fsthick) THEN
          ebnd = ekfs(ibndfs, ik) - ek_cl(ibnd, ik_fs_to_cl(ik))
          ebnd = ABS(ebnd)
          ebndmax4 = MAX(ebnd, ebndmax4)
          !!!!!!!! FOR DEBUG !!!!!!!!
          !WRITE(stdout, '(5x, f15.7, 1x, f15.7)') &
          !ekfs(ibndfs, ik), ek_cl(ibnd, ik_fs_to_cl(ik))
          !!!!!!!! FOR DEBUG !!!!!!!!
        ENDIF
      ENDDO
    ENDDO
    !
    WRITE(stdout, '(/5x, a, ES10.2, a/)') &
      'ekfs is in agreement with ek_cl within a difference of ', ebndmax1, ' eV.'
    IF (ABS(ebndmax2 - ebndmax1) > ebndmax1 * 1d-2) THEN
      WRITE(stdout, '(/5x, a, ES10.2, a/)') &
        'ekfs is in agreement with ek_cl within a difference(2) of ', ebndmax2, ' eV.'
    ELSEIF (ABS(ebndmax3 - ebndmax1) > ebndmax1 * 1d-2) THEN
      WRITE(stdout, '(/5x, a, ES10.2, a/)') &
        'ekfs is in agreement with ek_cl within a difference(3) of ', ebndmax3, ' eV.'
    ELSEIF (ABS(ebndmax4 - ebndmax1) > ebndmax1 * 1d-2) THEN
      WRITE(stdout, '(/5x, a, ES10.2, a/)') &
        'ekfs is in agreement with ek_cl within a difference(4) of ', ebndmax4, ' eV.'
    ENDIF
    IF (ebndmax1 > 0.2_DP) THEN
      WRITE(stdout, '(/5x, a/)') &
      'WARNING: The difference between ekfs and ek_cl is large on some of k points.'
    ENDIF
    !
    !WRITE(1000+mpime, '(/5x, a, ES10.2, a/)') &
    !  'ekfs is in agreement with ek_cl within a difference of ', ebndmax1, ' eV.'
    !IF (ebndmax1 > 0.2_DP) THEN
    !  WRITE(1000+mpime, '(/5x, a/)') &
    !  'WARNING: The difference between ekfs and ek_cl is large on some of k points.'
    !ENDIF
    !
    RETURN
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE find_nbnd_offset
    !-----------------------------------------------------------------------
    !
    !----------------------------------------------------------------------
    SUBROUTINE deallocate_coulomb()
    !----------------------------------------------------------------------
    !!
    !!  deallocates the variables allocated for Coulomb
    !!
    !
    USE supercond_common, ONLY : nkfs, nbnd_cl, nkstot_cl, &
                              wk_cl, xk_cl, ek_cl, &
                              ik_cl_to_fs, ik_fs_to_cl, &
                              xk_bz_cl, ik_bz_to_ibz_cl, &
                              ik_ibz_to_bz_cl, &
                              nk1_cl, nk2_cl, nk3_cl
    USE low_lvl,           ONLY : mem_size_eliashberg
    !
    IMPLICIT NONE
    !
    INTEGER :: ierr
    !! Error status
    INTEGER(8) :: imelt
    !! Counter memory
    !
    DEALLOCATE(wk_cl, STAT = ierr)
    IF (ierr /= 0) CALL errore('deallocate_coulomb', 'Error deallocating wk_cl', 1)
    DEALLOCATE(xk_cl, STAT = ierr)
    IF (ierr /= 0) CALL errore('deallocate_coulomb', 'Error deallocating xk_cl', 1)
    DEALLOCATE(ek_cl, STAT = ierr)
    IF (ierr /= 0) CALL errore('deallocate_coulomb', 'Error deallocating ek_cl', 1)
    DEALLOCATE(xk_bz_cl, STAT = ierr)
    IF (ierr /= 0) CALL errore('deallocate_coulomb', 'Error deallocating xk_bz_cl', 1)
    !
    ! remove memory allocated for wk_cl, xk_cl, ek_cl, xk_bz_cl
    imelt = (1 + 3 + nbnd_cl) * nkstot_cl + 3 * nk1_cl * nk2_cl * nk3_cl
    CALL mem_size_eliashberg(2, -imelt)
    !
    DEALLOCATE(ik_cl_to_fs, STAT = ierr)
    IF (ierr /= 0) CALL errore('deallocate_coulomb', 'Error deallocating ik_cl_to_fs', 1)
    DEALLOCATE(ik_fs_to_cl, STAT = ierr)
    IF (ierr /= 0) CALL errore('deallocate_coulomb', 'Error deallocating ik_fs_to_cl', 1)
    ! remove memory allocated for ik_cl_to_fs, ik_fs_to_cl
    imelt = nkstot_cl + nkfs
    CALL mem_size_eliashberg(1, -imelt)
    !
    DEALLOCATE(ik_bz_to_ibz_cl, STAT = ierr)
    IF (ierr /= 0) CALL errore('deallocate_coulomb', 'Error deallocating ik_bz_to_ibz_cl', 1)
    imelt = nk1_cl * nk2_cl * nk3_cl
    CALL mem_size_eliashberg(1, -imelt)
    DEALLOCATE(ik_ibz_to_bz_cl, STAT = ierr)
    IF (ierr /= 0) CALL errore('deallocate_coulomb', 'Error deallocating ik_ibz_to_bz_cl', 1)
    imelt = nkstot_cl
    CALL mem_size_eliashberg(1, -imelt)
    !
    RETURN
    !
    !-----------------------------------------------------------------------
    END SUBROUTINE deallocate_coulomb
    !-----------------------------------------------------------------------
    !
  END MODULE supercond_coul