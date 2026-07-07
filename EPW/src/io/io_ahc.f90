!
! Copyright (C) 2023-2026 EPW-Collaboration
! Copyright (C) 2010-2019 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino
!
! This file is distributed under the terms of the GNU General Public
! License. See the file `LICENSE' in the root directory of the
! present distribution, or http://www.gnu.org/copyleft.gpl.txt .
!
!----------------------------------------------------------------------------
MODULE io_ahc
!----------------------------------------------------------------------------
!!
!! This module contains routines for reading and processing files written
!! by ph.x with electron-phonon = 'ahc'.
!!
!! The read data are used for calculating the real part of the electron
!! self-energy within the Allen-Heine-Cardona (AHC) theory.
!!
IMPLICIT NONE
!
INTEGER :: nbnd_from_ahc_file = -1
!! nbnd read from ahc_info.dat. Set by read_ahc_info(), checked by check_ahc_bands().
!! This should equal nbnd from pwcom (the NSCF total number of bands).
!
CONTAINS
  !
  SUBROUTINE read_sthmat(iq_irr, lwin, lwinq)
  !----------------------------------------------------------------------------
  !!
  !! Read and calculate Sternheimer matrix from ahc_dir and save at sthmatq.
  !! See Eqs. (8, 32) of Lihm and Park, PRX 11, 041053 (2021).
  !!
  !! We want the following matrix element:
  !! sthmatq(ib, jb, ik, imode, jmode)
  !! = <d_{q, imode} psi_ib(k) | Q_outer(k+q) dV_{q, jmode} | psi_jb(k)>
  !!
  !! Q_outer(k+q): projection to the states at k+q above the upper limit of the
  !!               outer energy window (dis_win_max).
  !!
  !! However, in ph.x, a different projection operator is used:
  !! sth_mat_ph = <d_{q, imode} psi_ib(k) | Q_band(k+q) dV_{q, jmode} | psi_jb(k)>
  !!
  !! Q_band(k+q): projection to the states at k+q not included in the NSCF
  !!              calculation (bands with index nbnd + 1 and above).
  !!
  !! So, we add the missing contribution.
  !! sthmatq(m, n, ik, imode, jmode)
  !! = sth_mat_ph(m, n, ik, imode, jmode)
  !! + sum_{p for lwinq(p) == .FALSE.} (g_pm(k, q, imode))* * g_pn(k, q, jmode)
  !!                                   / (e_m(k) - e_p(k+q))
  !!
  !! where g_mn(k, q, imode) = <psi_m(k+q) | V_{q, imode} | psi_n(k)>
  !! computed in ph.x and read from file here.
  !!
  !! (imode and jmode are in the Cartesian basis)
  !!
  !! We assume the following constraints:
  !! 1) The bands included in the Wannierization should pass the
  !!    degeneracy_check.py test.
  !!    One should use bands_skipped input to exclude other bands.
  !! 2) The bands computed in ahc.f90 should be identical to the bands included
  !!    in the Wannierization (i.e. the bands that are not excluded by
  !!    bands_skipped).
  !! 3) One should not exclude intermediate bands by bands_skipped.
  !!    ex) bands_skipped = 'exclude_bands = 1-5, 8, 10-12' is not allowed.
  !! Then, one can safely assume that the band indices from ph.x ahc output
  !! is identical to that of EPW, from 1 to nbndep.
  !! (nbndep = ahc_nbnd)
  !!
  !! TODO: band (image) parallelization
  !!
  !----------------------------------------------------------------------------
  !
  USE kinds,            ONLY : DP
  USE ep_constants,     ONLY : zero, one, czero, eps4, ryd2ev
  USE pwcom,            ONLY : nbnd, nks
  USE klist,            ONLY : nkstot
  USE modes,            ONLY : nmodes
  USE input,            ONLY : ahc_nbnd, ahc_nbndskip, dvscf_dir, &
                               ahc_rest_fan_broadening
  USE global_var,       ONLY : nbndep, sthmatq
  USE kfold,            ONLY : ktokpmq
  USE parallelism,      ONLY : fkbounds
  USE io_var,           ONLY : iuahcsth, iuahcgkk, iuahcet
  !
  IMPLICIT NONE
  !
  LOGICAL, INTENT(in) :: lwin(nbndep, nks)
  !! Outer windows at k+q
  LOGICAL, INTENT(in) :: lwinq(nbndep, nks)
  !! Outer windows at k+q
  INTEGER, INTENT(in) :: iq_irr
  !! Current q-point index in the list of irreducible q-points.
  !
  CHARACTER(LEN=256) :: filesth
  !! Filename of the ph.x upper Fan data
  CHARACTER(LEN=256) :: fileepph
  !! Filename of the ph.x e-ph matrix data
  CHARACTER(LEN=256) :: fileetk
  !! Filename of the ph.x energy at k data
  CHARACTER(LEN=256) :: fileetq
  !! Filename of the ph.x energy at k+q data
  INTEGER :: ik
  !! Counter on k-points in the pool
  INTEGER :: ibnd
  !! Counter on bands
  INTEGER :: jbnd
  !! Counter on bands
  INTEGER :: kbnd
  !! Counter on bands
  INTEGER :: imode
  !! Counter on modes
  INTEGER :: jmode
  !! Counter on modes
  INTEGER :: lower_bnd
  !! Lower bounds index after k paral
  INTEGER :: upper_bnd
  !! Upper bounds index after k paral
  INTEGER :: ib_win_fst
  !! First band inside the outer window
  INTEGER :: ib_win_lst
  !! First band inside the outer window
  INTEGER :: ik_global
  !! Global index of k point (for pool parallelization)
  INTEGER :: recl
  !! Record length of the ph.x data
  INTEGER :: ios
  !! Integer variable for I/O control
  INTEGER :: ierr
  !! Error status
  REAL(KIND = DP) :: delta_e
  !! Energy difference
  REAL(KIND = DP), ALLOCATABLE :: inv_delta_e(:, :)
  !! 1 / delta_e
  REAL(KIND = DP), ALLOCATABLE :: etk_ph(:, :)
  !! Electron energy at k read from ph.x output
  REAL(KIND = DP), ALLOCATABLE :: etq_ph(:, :)
  !! Electron energy at k+q read from ph.x output
  COMPLEX(KIND = DP), ALLOCATABLE :: ep_mat_ph(:, :, :)
  !! E-ph matrix read from ph.x output
  COMPLEX(KIND = DP), ALLOCATABLE :: sth_mat(:, :, :, :)
  !! Upper Fan matrix read from file.
  COMPLEX(KIND = DP), ALLOCATABLE :: sth_mat_temp(:, :)
  !! Temporary storage for upper Fan matrix.
  !
  CHARACTER(LEN=6), EXTERNAL :: int_to_char
  !
  CALL start_clock('read_sthmat')
  !
  IF (nbndep /= ahc_nbnd) CALL errore('read_sthmat', 'nbndep /= ahc_nbnd', 1)
  !
  ALLOCATE(sth_mat_temp(nbndep, nbndep), STAT = ierr)
  IF (ierr /= 0) CALL errore('read_sthmat', 'Error allocating sth_mat_temp', 1)
  ALLOCATE(sth_mat(ahc_nbnd, ahc_nbnd, nmodes, nmodes), STAT = ierr)
  IF (ierr /= 0) CALL errore('read_sthmat', 'Error allocating sth_mat', 1)
  ALLOCATE(ep_mat_ph(nbnd, ahc_nbnd, nmodes), STAT = ierr)
  IF (ierr /= 0) CALL errore('read_sthmat', 'Error allocating ep_mat_ph', 1)
  ALLOCATE(etk_ph(nbnd, nkstot), STAT = ierr)
  IF (ierr /= 0) CALL errore('read_sthmat', 'Error allocating etk_ph', 1)
  ALLOCATE(etq_ph(nbnd, nkstot), STAT = ierr)
  IF (ierr /= 0) CALL errore('read_sthmat', 'Error allocating etq_ph', 1)
  ALLOCATE(inv_delta_e(nbnd, ahc_nbnd), STAT = ierr)
  IF (ierr /= 0) CALL errore('read_sthmat', 'Error allocating inv_delta_e', 1)
  !
  ! find the bounds of k-dependent arrays in the parallel case in each pool
  CALL fkbounds(nkstot, lower_bnd, upper_bnd)
  !
  ! Open ph.x data file
  !
  filesth = TRIM(dvscf_dir) // 'ahc_dir/ahc_upfan_iq' // TRIM(int_to_char(iq_irr)) // '.bin'
  fileepph = TRIM(dvscf_dir) // 'ahc_dir/ahc_gkk_iq' // TRIM(int_to_char(iq_irr)) // '.bin'
  fileetk = TRIM(dvscf_dir) // 'ahc_dir/ahc_etk_iq' // TRIM(int_to_char(iq_irr)) // '.bin'
  fileetq = TRIM(dvscf_dir) // 'ahc_dir/ahc_etq_iq' // TRIM(int_to_char(iq_irr)) // '.bin'
  !
  INQUIRE(IOLENGTH=recl) sth_mat
  OPEN(UNIT=iuahcsth, FILE=TRIM(filesth), ACTION='read', FORM='unformatted', &
       ACCESS='direct', STATUS='old', RECL=recl, IOSTAT=ios)
  IF (ios /= 0) CALL errore('read_sthmat', 'error opening '//TRIM(filesth), 1)
  !
  INQUIRE(IOLENGTH=recl) ep_mat_ph
  OPEN(UNIT=iuahcgkk, FILE=TRIM(fileepph), ACTION='read', FORM='unformatted', &
       ACCESS='direct', STATUS='old', RECL=recl, IOSTAT=ios)
  IF (ios /= 0) CALL errore('read_sthmat', 'error opening '//TRIM(fileepph), 1)
  !
  INQUIRE(IOLENGTH=recl) etk_ph
  OPEN(UNIT=iuahcet, FILE=TRIM(fileetk), ACTION='read', FORM='unformatted', &
       ACCESS='direct', STATUS='old', RECL=recl, IOSTAT=ios)
  IF (ios /= 0) CALL errore('read_sthmat', 'error opening '//TRIM(fileetk), 1)
  READ(iuahcet, REC=1) etk_ph
  CLOSE(iuahcet)
  !
  INQUIRE(IOLENGTH=recl) etq_ph
  OPEN(UNIT=iuahcet, FILE=TRIM(fileetq), ACTION='read', FORM='unformatted', &
       ACCESS='direct', STATUS='old', RECL=recl, IOSTAT=ios)
  IF (ios /= 0) CALL errore('read_sthmat', 'error opening '//TRIM(fileetq), 1)
  READ(iuahcet, REC=1) etq_ph
  CLOSE(iuahcet)
  !
  sthmatq(:, :, :, :, :) = czero
  !
  DO ik = 1, nks
    !
    ik_global = ik + lower_bnd - 1
    !
    ! Read ahc_dir files.
    !
    READ(iuahcsth, REC=ik_global) sth_mat
    READ(iuahcgkk, REC=ik_global) ep_mat_ph
    !
    ! Add additional term coming from states outside the outer window.
    !
    ! inv_delta_e(kbnd, ibnd) = 1 / ( e_ibnd(k) - e_kbnd(k+q) )
    ! inv_delta_e(kbnd, ibnd) = Re 1 / ( e_ibnd(k) - e_kbnd(k+q) + im * 100 meV)
    !
    DO ibnd = 1, ahc_nbnd
      DO kbnd = 1, nbnd
        delta_e = etk_ph(ibnd + ahc_nbndskip, ik_global) - etq_ph(kbnd, ik_global)
        inv_delta_e(kbnd, ibnd) = REAL(one / CMPLX(delta_e, ahc_rest_fan_broadening / ryd2ev, KIND=DP), KIND=DP)
      ENDDO ! kbnd
    ENDDO ! ibnd
    !
    ! Compute Sternheimer matrix elements
    !
    DO jmode = 1, nmodes
      DO imode = 1, nmodes
        DO jbnd = 1, ahc_nbnd
          !
          ! Skip jbnd outside the outer window
          IF (.NOT. lwin(jbnd, ik)) CYCLE
          !
          DO ibnd = 1, ahc_nbnd
            !
            ! Skip ibnd outside the outer window
            IF (.NOT. lwin(ibnd, ik)) CYCLE
            !
            DO kbnd = 1, nbnd
              !
              ! Only use psi_kbnd(k+q) outside the outer window
              !
              IF ( (kbnd - ahc_nbndskip >= 1) .AND. &
                   (kbnd - ahc_nbndskip <= nbndep) ) THEN
                IF (lwinq(kbnd - ahc_nbndskip, ik)) CYCLE
              ENDIF
              !
              sth_mat(ibnd, jbnd, imode, jmode) = sth_mat(ibnd, jbnd, imode, jmode) &
              + CONJG(ep_mat_ph(kbnd, ibnd, imode)) * ep_mat_ph(kbnd, jbnd, jmode) &
              * inv_delta_e(kbnd, ibnd)
              !
            ENDDO ! kbnd
          ENDDO ! ibnd
        ENDDO ! jbnd
      ENDDO ! imode
    ENDDO ! jmode
    !
    sthmatq(:, :, ik, :, :) = sth_mat
    !
  ENDDO
  !
  ! Slim down sthmatq to the states inside the outer window. Matrix elements
  ! for bands outside the outer window are removed.
  !
  DO ik = 1, nks
    !
    ! Find index of first and last band inside the outer window
    !
    ib_win_fst = -1
    ib_win_lst = -1
    !
    DO ibnd = 1, nbndep
      IF (lwin(ibnd, ik)) THEN
        ib_win_fst = ibnd
        EXIT
      ENDIF
    ENDDO
    !
    DO ibnd = nbndep, 1, -1
      IF (lwin(ibnd, ik)) THEN
        ib_win_lst = ibnd
        EXIT
      ENDIF
    ENDDO
    !
    IF (ib_win_fst == -1) CALL errore('read_sthmat', 'ib_win_fst not found', 1)
    IF (ib_win_lst == -1) CALL errore('read_sthmat', 'ib_win_lst not found', 1)
    !
    DO jmode = 1, nmodes
      DO imode = 1, nmodes
        !
        sth_mat_temp = czero
        !
        DO jbnd = ib_win_fst, ib_win_lst
          DO ibnd = ib_win_fst, ib_win_lst
            sth_mat_temp(ibnd - ib_win_fst + 1, jbnd - ib_win_fst + 1) &
            = sthmatq(ibnd, jbnd, ik, imode, jmode)
          ENDDO
        ENDDO
        !
        sthmatq(:, :, ik, imode, jmode) = sth_mat_temp
        !
      ENDDO ! imode
    ENDDO ! jmode
  ENDDO ! ik
  !
  CLOSE(iuahcsth)
  CLOSE(iuahcgkk)
  !
  DEALLOCATE(sth_mat, STAT = ierr)
  IF (ierr /= 0) CALL errore('read_sthmat', 'Error deallocating sth_mat_smooth', 1)
  DEALLOCATE(ep_mat_ph, STAT = ierr)
  IF (ierr /= 0) CALL errore('read_sthmat', 'Error deallocating ep_mat_ph', 1)
  DEALLOCATE(etk_ph, STAT = ierr)
  IF (ierr /= 0) CALL errore('read_sthmat', 'Error deallocating etk_ph', 1)
  DEALLOCATE(etq_ph, STAT = ierr)
  IF (ierr /= 0) CALL errore('read_sthmat', 'Error deallocating etq_ph', 1)
  DEALLOCATE(inv_delta_e, STAT = ierr)
  IF (ierr /= 0) CALL errore('read_sthmat', 'Error deallocating inv_delta_e', 1)
  !
  CALL stop_clock('read_sthmat')
  !
  !----------------------------------------------------------------------------
  END SUBROUTINE read_sthmat
  !----------------------------------------------------------------------------
  !
  !----------------------------------------------------------------------------
  SUBROUTINE read_dwmat(lwin)
  !----------------------------------------------------------------------------
  !!
  !! Read Debye-Waller matrix from output of ph.x and store it to dw_mat.
  !! The file is located in dvscf_dir/ahc_dir.
  !!
  !! TODO: band (image) parallelization
  !!
  !----------------------------------------------------------------------------
  !
  USE kinds,            ONLY : DP
  USE ep_constants,     ONLY : czero
  USE pwcom,            ONLY : nks
  USE klist,            ONLY : nkstot
  USE modes,            ONLY : nmodes
  USE input,            ONLY : ahc_nbnd, dvscf_dir
  USE global_var,       ONLY : dw_mat, nbndep
  USE parallelism,      ONLY : fkbounds
  USE io_var,           ONLY : iuahcdw
  !
  IMPLICIT NONE
  !
  LOGICAL, INTENT(in) :: lwin(nbndep, nks)
  !! Bands at k within outer energy window
  !
  CHARACTER(LEN=256) :: filedw
  !! Filename of the ph.x Debye-Waller data
  INTEGER :: lower_bnd
  !! Lower bounds index after k paral
  INTEGER :: upper_bnd
  !! Upper bounds index after k paral
  INTEGER :: ik
  !! Counter on k-points in the pool
  INTEGER :: ik_global
  !! Global index of k point (for pool parallelization)
  INTEGER :: ibnd
  !! COunter on bands
  INTEGER :: ib_win_fst
  !! Index of first band inside the outer window
  INTEGER :: ib_win_lst
  !! Index of last band inside the outer window
  INTEGER :: ndimwin_ik
  !! Number of bands inside the window. Identical to ndimwin(ik).
  INTEGER :: imode
  !! Counter on modes
  INTEGER :: idir
  !! Counter on directions
  INTEGER :: recl
  !! Record length of the ph.x data
  INTEGER :: ios
  !! Integer variable for I/O control
  INTEGER :: ierr
  !! Error status
  COMPLEX(KIND = DP), ALLOCATABLE :: dw_mat_temp(:, :, :, :)
  !! Temporary containor for Debye-Waller matrix element
  !
  CALL start_clock('read_dwmat')
  !
  IF (nbndep /= ahc_nbnd) CALL errore('read_dwmat', 'nbndep /= ahc_nbnd', 1)
  !
  ALLOCATE(dw_mat_temp(ahc_nbnd, ahc_nbnd, nmodes, 3), STAT = ierr)
  IF (ierr /= 0) CALL errore('read_dwmat', 'Error allocating dw_mat_temp', 1)
  !
  ! find the bounds of k-dependent arrays in the parallel case in each pool
  CALL fkbounds(nkstot, lower_bnd, upper_bnd)
  !
  dw_mat = czero
  !
  ! Open ph.x data file
  !
  filedw = TRIM(dvscf_dir) // 'ahc_dir/ahc_dw.bin'
  !
  INQUIRE(IOLENGTH=recl) dw_mat_temp
  OPEN(UNIT=iuahcdw, FILE=TRIM(filedw), ACTION='read', FORM='unformatted', &
       ACCESS='direct', STATUS='old', RECL=recl, IOSTAT=ios)
  IF (ios /= 0) CALL errore('read_dwmat', 'error opening '//TRIM(filedw), 1)
  !
  ! Read ahc_dir output
  !
  DO ik = 1, nks
    ik_global = ik + lower_bnd - 1
    READ(iuahcdw, REC=ik_global) dw_mat_temp
    !
    ! Srim down to states only inside the outer window
    ! (See subroutine rotate_epmat)
    !
    ib_win_fst = -1
    ib_win_lst = -1
    !
    DO ibnd = 1, nbndep
      IF (lwin(ibnd, ik)) THEN
        ib_win_fst = ibnd
        EXIT
      ENDIF
    ENDDO
    !
    DO ibnd = nbndep, 1, -1
      IF (lwin(ibnd, ik)) THEN
        ib_win_lst = ibnd
        EXIT
      ENDIF
    ENDDO
    !
    IF (ib_win_fst == -1) CALL errore('read_dwmat', 'ib_win_fst not found', 1)
    IF (ib_win_lst == -1) CALL errore('read_dwmat', 'ib_win_lst not found', 1)
    !
    ndimwin_ik = ib_win_lst - ib_win_fst + 1
    !
    ! Reshape:
    ! dw_mat_temp (from ph.x) has shape (ahc_nbnd, ahc_nbnd, nmodes, 3)
    ! dw_mat (for epw)        has shape (ahc_nbnd, ahc_nbnd, nks, 3, nmodes)
    !
    DO imode = 1, nmodes
      DO idir = 1, 3
        dw_mat(1:ndimwin_ik, 1:ndimwin_ik, ik, idir, imode) &
        = dw_mat_temp(ib_win_fst:ib_win_lst, ib_win_fst:ib_win_lst, imode, idir)
      ENDDO
    ENDDO
    !
  ENDDO
  !
  CLOSE(iuahcdw)
  !
  DEALLOCATE(dw_mat_temp, STAT = ierr)
  IF (ierr /= 0) CALL errore('read_dwmat', 'Error deallocating dw_mat_temp', 1)
  !
  CALL stop_clock('read_dwmat')
  !
  !--------------------------------------------------------------------------
  END SUBROUTINE read_dwmat
  !--------------------------------------------------------------------------
  !
  !--------------------------------------------------------------------------
  SUBROUTINE read_ahc_info()
  !--------------------------------------------------------------------------
  !!
  !! Read ahc_info.dat from ahc_dir to obtain ahc_nbnd, ahc_nbndskip, and
  !! nbnd used in the ph.x AHC calculation.
  !! If the user set ahc_nbnd or ahc_nbndskip in the EPW input, print a
  !! deprecation warning and override with the values from the file.
  !!
  !--------------------------------------------------------------------------
  !
  USE io_global,        ONLY : stdout, ionode, ionode_id
  USE input,            ONLY : ahc_nbnd, ahc_nbndskip, dvscf_dir
  USE mp,               ONLY : mp_bcast
  USE mp_global,        ONLY : world_comm
  !
  IMPLICIT NONE
  !
  CHARACTER(LEN=256) :: filename
  !! Path to ahc_info.dat
  CHARACTER(LEN=256) :: label
  !! Label read from file
  INTEGER :: iun
  !! File unit
  INTEGER :: ios
  !! I/O status
  INTEGER :: ahc_nbnd_file
  !! ahc_nbnd read from file
  INTEGER :: ahc_nbndskip_file
  !! ahc_nbndskip read from file
  INTEGER :: nbnd_file
  !! nbnd read from file
  LOGICAL :: file_exists
  !! Whether ahc_info.dat exists
  !
  filename = TRIM(dvscf_dir) // 'ahc_dir/ahc_info.dat'
  !
  IF (ionode) THEN
    !
    INQUIRE(FILE=TRIM(filename), EXIST=file_exists)
    !
    IF (.NOT. file_exists) THEN
      WRITE(stdout, '(5x,a)') 'ERROR: File not found: ' // TRIM(filename)
      WRITE(stdout, '(5x,a)') 'Please rerun ph.x with electron_phonon = "ahc" to generate this file,'
      WRITE(stdout, '(5x,a)') 'or manually create ahc_dir/ahc_info.dat (see PHonon/PH/ahc.f90 for format).'
    ENDIF
    !
  ENDIF
  !
  CALL mp_bcast(file_exists, ionode_id, world_comm)
  IF (.NOT. file_exists) CALL errore('read_ahc_info', &
      'ahc_info.dat not found in ahc_dir. Rerun ph.x AHC or create it manually.', 1)
  !
  IF (ionode) THEN
    !
    OPEN(NEWUNIT=iun, FILE=TRIM(filename), FORM='formatted', STATUS='old', IOSTAT=ios)
    IF (ios /= 0) CALL errore('read_ahc_info', 'Error opening ' // TRIM(filename), 1)
    !
    READ(iun, '(A14, I8)', IOSTAT=ios) label, ahc_nbnd_file
    IF (ios /= 0) CALL errore('read_ahc_info', 'Error reading ahc_nbnd from ahc_info.dat', 1)
    READ(iun, '(A14, I8)', IOSTAT=ios) label, ahc_nbndskip_file
    IF (ios /= 0) CALL errore('read_ahc_info', 'Error reading ahc_nbndskip from ahc_info.dat', 1)
    READ(iun, '(A14, I8)', IOSTAT=ios) label, nbnd_file
    IF (ios /= 0) CALL errore('read_ahc_info', 'Error reading nbnd from ahc_info.dat', 1)
    !
    CLOSE(iun)
    !
    WRITE(stdout, '(5x,a,I5)') 'AHC info: ahc_nbnd     = ', ahc_nbnd_file
    WRITE(stdout, '(5x,a,I5)') 'AHC info: ahc_nbndskip = ', ahc_nbndskip_file
    WRITE(stdout, '(5x,a,I5)') 'AHC info: nbnd         = ', nbnd_file
    !
    ! Override with values from file
    !
    ahc_nbnd = ahc_nbnd_file
    ahc_nbndskip = ahc_nbndskip_file
    nbnd_from_ahc_file = nbnd_file
    !
  ENDIF
  !
  CALL mp_bcast(ahc_nbnd, ionode_id, world_comm)
  CALL mp_bcast(ahc_nbndskip, ionode_id, world_comm)
  CALL mp_bcast(nbnd_from_ahc_file, ionode_id, world_comm)
  !
  !--------------------------------------------------------------------------
  END SUBROUTINE read_ahc_info
  !--------------------------------------------------------------------------
  !
  !--------------------------------------------------------------------------
  SUBROUTINE check_ahc_bands(excluded_band)
  !--------------------------------------------------------------------------
  !!
  !! Check consistency between ahc_nbnd, ahc_nbndskip, nbndep, nbnd, and
  !! the excluded_band array from bands_skipped.
  !! Must be called after read_ahc_info() and setup_nnkp().
  !!
  !--------------------------------------------------------------------------
  !
  USE io_global,        ONLY : stdout
  USE pwcom,            ONLY : nbnd
  USE input,            ONLY : ahc_nbnd, ahc_nbndskip
  USE global_var,       ONLY : nbndep
  !
  IMPLICIT NONE
  !
  LOGICAL, INTENT(IN) :: excluded_band(:)
  !! Boolean array of excluded bands (size nbnd)
  !
  INTEGER :: ibnd
  !! Band counter
  LOGICAL :: has_error
  !! Whether any error was found
  CHARACTER(LEN=256) :: correct_exclude
  !! The correct exclude_bands string
  CHARACTER(LEN=6), EXTERNAL :: int_to_char
  !
  ! Check nbnd consistency (nbnd from pwcom is set by this point)
  !
  IF (nbnd_from_ahc_file > 0 .AND. nbnd_from_ahc_file /= nbnd) THEN
    WRITE(stdout, '(5x,a,I5)') 'ERROR: nbnd from ahc_info.dat = ', nbnd_from_ahc_file
    WRITE(stdout, '(5x,a,I5)') '       nbnd from NSCF         = ', nbnd
    CALL errore('check_ahc_bands', &
        'nbnd mismatch between AHC calculation and current NSCF. Rerun AHC with the same nbnd.', 1)
  ENDIF
  !
  has_error = .FALSE.
  !
  ! Check low bands (1:ahc_nbndskip) are all excluded
  !
  DO ibnd = 1, ahc_nbndskip
    IF (.NOT. excluded_band(ibnd)) THEN
      WRITE(stdout, '(5x,a,I5,a)') 'ERROR: Band ', ibnd, &
          ' is not excluded but must be (below ahc_nbndskip).'
      has_error = .TRUE.
    ENDIF
  ENDDO
  !
  ! Check high bands (ahc_nbndskip+ahc_nbnd+1:nbnd) are all excluded
  !
  DO ibnd = ahc_nbndskip + ahc_nbnd + 1, nbnd
    IF (.NOT. excluded_band(ibnd)) THEN
      WRITE(stdout, '(5x,a,I5,a)') 'ERROR: Band ', ibnd, &
          ' is not excluded but must be (above ahc_nbndskip + ahc_nbnd).'
      has_error = .TRUE.
    ENDIF
  ENDDO
  !
  ! Check middle bands (ahc_nbndskip+1:ahc_nbndskip+ahc_nbnd) are NOT excluded
  !
  DO ibnd = ahc_nbndskip + 1, ahc_nbndskip + ahc_nbnd
    IF (excluded_band(ibnd)) THEN
      WRITE(stdout, '(5x,a,I5,a)') 'ERROR: Band ', ibnd, &
          ' is excluded but must not be (inside AHC window).'
      has_error = .TRUE.
    ENDIF
  ENDDO
  !
  IF (has_error) THEN
    !
    ! Build the correct exclude_bands string
    !
    IF (ahc_nbndskip > 0 .AND. ahc_nbndskip + ahc_nbnd < nbnd) THEN
      correct_exclude = '1:' // TRIM(ADJUSTL(int_to_char(ahc_nbndskip))) &
          // ',' // TRIM(ADJUSTL(int_to_char(ahc_nbndskip + ahc_nbnd + 1))) &
          // ':' // TRIM(ADJUSTL(int_to_char(nbnd)))
    ELSEIF (ahc_nbndskip > 0) THEN
      correct_exclude = '1:' // TRIM(ADJUSTL(int_to_char(ahc_nbndskip)))
    ELSEIF (ahc_nbndskip + ahc_nbnd < nbnd) THEN
      correct_exclude = TRIM(ADJUSTL(int_to_char(ahc_nbndskip + ahc_nbnd + 1))) &
          // ':' // TRIM(ADJUSTL(int_to_char(nbnd)))
    ELSE
      correct_exclude = '(none)'
    ENDIF
    !
    WRITE(stdout, '(5x,a)')
    WRITE(stdout, '(5x,a,I5)') 'ahc_nbndskip (from ahc_info.dat) = ', ahc_nbndskip
    WRITE(stdout, '(5x,a,I5)') 'ahc_nbnd     (from ahc_info.dat) = ', ahc_nbnd
    WRITE(stdout, '(5x,a,I5)') 'nbnd         (from ahc_info.dat) = ', nbnd
    WRITE(stdout, '(5x,a)')
    WRITE(stdout, '(5x,a)') 'You must set bands_skipped = ''exclude_bands = ' &
        // TRIM(correct_exclude) // ''''
    WRITE(stdout, '(5x,a)')
    !
    CALL errore('check_ahc_bands', 'Inconsistent bands_skipped for WFPT/AHC', 1)
    !
  ENDIF
  !
  ! Also check nbndep == ahc_nbnd (should be guaranteed by the above checks,
  ! but verify explicitly)
  !
  IF (nbndep /= ahc_nbnd) THEN
    WRITE(stdout, '(5x,a,I5)') 'ERROR: nbndep = ', nbndep
    WRITE(stdout, '(5x,a,I5)') '       ahc_nbnd = ', ahc_nbnd
    CALL errore('check_ahc_bands', 'nbndep /= ahc_nbnd (internal consistency error)', 1)
  ENDIF
  !
  WRITE(stdout, '(5x,a)') 'AHC band consistency checks passed.'
  !
  !--------------------------------------------------------------------------
  END SUBROUTINE check_ahc_bands
  !--------------------------------------------------------------------------
  !
!------------------------------------------------------------------------------
END MODULE io_ahc
!------------------------------------------------------------------------------
