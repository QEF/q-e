  !
  ! Copyright (C) 2024 EPW-Collaboration
  !
  ! This file is distributed under the terms of the GNU General Public
  ! License. See the file `LICENSE' in the root directory of the
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .
  !
  !-----------------------------------------------------------------------
  PROGRAM nscf2supercond
  !-----------------------------------------------------------------------
  !!
  !! This program extracts the data of the nscf calculation.
  !!
  USE io_files,  ONLY : prefix, tmp_dir
  USE mp_global, ONLY : mp_startup
  USE mp_pools,  ONLY : npool
  USE environment,   ONLY : environment_start, environment_end
  USE io_global, ONLY : meta_ionode, meta_ionode_id, stdout
  USE mp_global, ONLY : world_comm
  USE mp,        ONLY : mp_bcast, mp_barrier
  !
  IMPLICIT NONE
  !
  CHARACTER(LEN=256), EXTERNAL :: trimcheck
  !
  CHARACTER (len=256) :: filband, outdir
  INTEGER :: ios
  !! IO error message
  !
  NAMELIST / bands / outdir, prefix, filband
  !
  ! initialise environment
  !
#if defined(__MPI)
  CALL mp_startup ( )
#endif
  CALL environment_start ( 'BANDS' )
  !
  !   set default values for variables in namelist
  !
  prefix = 'pwscf'
  CALL get_environment_variable( 'ESPRESSO_TMPDIR', outdir )
  IF ( trim( outdir ) == ' ' ) outdir = './'
  filband = 'bands.out'
  !
  ios = 0
  !
  IF ( meta_ionode )  THEN
     !
     CALL input_from_file ( )
     !
     READ (5, bands, iostat = ios)
     !
     tmp_dir = trimcheck (outdir)
     !
  ENDIF
  !
  !
  CALL mp_bcast( ios, meta_ionode_id, world_comm)
  IF (ios /= 0) CALL errore ('bands2epw', 'reading bands namelist', abs(ios) )
  !
  ! ... Broadcast variables
  !
  CALL mp_bcast( tmp_dir, meta_ionode_id, world_comm)
  CALL mp_bcast( prefix, meta_ionode_id, world_comm)
  CALL mp_bcast( filband, meta_ionode_id, world_comm)
  !
  !   Now allocate space for pwscf variables, read and check them.
  !
  CALL read_file()
  !
  CALL write_eigen_band (filband)
  !
  CALL environment_end ( 'BANDS' )
  !
  CALL stop_pp()
  STOP
  !
  !-----------------------------------------------------------------------
  END PROGRAM nscf2supercond
  !-----------------------------------------------------------------------
  !
  !-----------------------------------------------------------------------
  SUBROUTINE write_eigen_band (filband)
  !-----------------------------------------------------------------------
  !!
  !! This routine writes the band energies on a file for epw.
  !!
  !
  USE kinds,                ONLY : dp
  USE cell_base,            ONLY : alat, at
  USE symm_base,            ONLY : nrot, nsym, s, t_rev, time_reversal
  USE constants,            ONLY : rytoev
  USE klist,                ONLY : xk, nks, nkstot, wk
  USE start_k,              ONLY : nk1, nk2, nk3, k1, k2, k3
  USE io_files,             ONLY : iunpun
  USE wvfct,                ONLY : nbnd, et
  USE noncollin_module,     ONLY : noncolin, lspinorb
  USE control_flags,        ONLY : noinv
  USE lsda_mod,             ONLY : nspin, lsda
  USE io_global,            ONLY : meta_ionode, meta_ionode_id, stdout
  USE mp,                   ONLY : mp_bcast
  !
  IMPLICIT NONE
  !
  CHARACTER (len=*) :: filband
  !
  INTEGER :: t_rev_tmp(48)
  !! it is used to copy t_rev(1:nrot)
  INTEGER :: s_tmp(3, 3, 48)
  !! it is used to copy s(:,:,1:nrot)
  INTEGER :: ik
  !! Counter on k points
  INTEGER :: isym
  !! Counter on symmetry
  INTEGER :: i
  !! Counter
  INTEGER :: j
  !! Counter
  INTEGER :: ios
  !! IO error message
  !
  ! HM: Since s is not properly initialized,
  !     s(:,:,nrot+1:48) may contain uninitialized values, depending on the compiler.
  s_tmp(:, :, :) = 0
  s_tmp(:, :, 1:nrot) = s(:, :, 1:nrot)
  t_rev_tmp(:) = 0
  t_rev_tmp(1:nrot) = t_rev(1:nrot)
  !
  IF ( meta_ionode ) THEN
    OPEN (unit = iunpun, file = filband, status = 'unknown', form = &
    'formatted', position = 'rewind', iostat = ios)
    WRITE(iunpun, '(ES23.15)') alat
    WRITE(iunpun, '(ES23.15,2(1x,ES23.15))') at(:,1)
    WRITE(iunpun, '(ES23.15,2(1x,ES23.15))') at(:,2)
    WRITE(iunpun, '(ES23.15,2(1x,ES23.15))') at(:,3)
    WRITE(iunpun, '(I5,1x,I5)') nrot, nsym
    WRITE(iunpun, '(L1)') time_reversal
    WRITE(iunpun, '(I2,15(1x,I2))') t_rev_tmp(1:16)
    WRITE(iunpun, '(I2,15(1x,I2))') t_rev_tmp(17:32)
    WRITE(iunpun, '(I2,15(1x,I2))') t_rev_tmp(33:48)
    DO isym = 1, 48
      WRITE(iunpun, '(I3,8(1x,I3))') ((s_tmp(i, j, isym), i=1,3), j=1,3)
    ENDDO
    WRITE(iunpun, '(L1,1x,I3)') lsda, nspin
    WRITE(iunpun, '(L1)') noinv
    WRITE(iunpun, '(L1,1x,L1)') noncolin, lspinorb
    WRITE(iunpun, '(I6,2(1x,I6))') nk1, nk2, nk3
    WRITE(iunpun, '(I6,2(1x,I6))') k1, k2, k3
    WRITE(iunpun, '(I6,1x,I10)') nbnd, nkstot
    !
    DO ik = 1, nkstot
      WRITE(iunpun, '(I10,1x,ES23.15)') ik, wk(ik)
      WRITE(iunpun, '(ES23.15,2(1x,ES23.15))') xk(:, ik)
      WRITE(iunpun, '(4(ES24.15))') et(:, ik)
    ENDDO
    CLOSE(iunpun)
  ENDIF
  !
  RETURN
  !
  !-----------------------------------------------------------------------
  END SUBROUTINE write_eigen_band
  !-----------------------------------------------------------------------
  
  
