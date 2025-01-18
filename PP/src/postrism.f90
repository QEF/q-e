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
PROGRAM pprism
  !--------------------------------------------------------------------------
  !
  ! ... Program to plot solvent distributions
  ! ... calculated by 3D-RISM or Laue-RISM
  !
  USE io_global,   ONLY : ionode
  USE mp_global,   ONLY : mp_startup
  USE environment, ONLY : environment_start, environment_end
  !
  IMPLICIT NONE
  !
  CHARACTER(LEN=256) :: filplot
  LOGICAL            :: lpunch
  !
#if defined(__MPI)
  CALL mp_startup()
#endif
  CALL environment_start('POST-RISM')
  !
  IF (ionode) CALL input_from_file()
  !
  CALL extract_rism(filplot, lpunch)
  !
  CALL solvdens(filplot, lpunch)
  !
  CALL environment_end('POST-RISM')
  !
  CALL stop_pp()
  !
END PROGRAM pprism
!
!--------------------------------------------------------------------------
SUBROUTINE extract_rism(filplot, lpunch)
  !--------------------------------------------------------------------------
  !
  USE io_files,      ONLY : tmp_dir, prefix
  USE io_global,     ONLY : stdin, ionode, ionode_id
  USE mp,            ONLY : mp_bcast
  USE mp_images,     ONLY : intra_image_comm
  USE rism3d_facade, ONLY : lrism3d
  !
  IMPLICIT NONE
  !
  CHARACTER(LEN=256), INTENT(OUT) :: filplot
  LOGICAL,            INTENT(OUT) :: lpunch
  !
  LOGICAL                      :: needwf
  INTEGER                      :: ios
  CHARACTER(LEN=256)           :: outdir
  CHARACTER(LEN=256), EXTERNAL :: trimcheck
  !
  NAMELIST / inputpp / outdir, prefix, filplot, lpunch
  !
  ! ... set default values for variables in namelist
  prefix  = 'pwscf'
  filplot = ''
  lpunch  = .FALSE.
  !
  CALL get_environment_variable('ESPRESSO_TMPDIR', outdir)
  IF (LEN_TRIM(outdir) < 1) THEN
    outdir = './'
  END IF
  !
  ! ... read the namelist inputpp
  IF (ionode)  THEN
    READ(stdin, nml=inputpp, iostat=ios)
    !
    tmp_dir = trimcheck(outdir)
    !
    IF (LEN_TRIM(filplot) < 1) THEN
      filplot = TRIM(ADJUSTL(prefix)) // '.pprism'
    END IF
  ENDIF
  !
  CALL mp_bcast(ios, ionode_id, intra_image_comm)
  !
  IF (ios /= 0) THEN
    CALL errore('postrism', 'reading inputpp namelist', ABS(ios))
  END IF
  !
  ! ... broadcast variables
  CALL mp_bcast(tmp_dir, ionode_id, intra_image_comm)
  CALL mp_bcast(prefix,  ionode_id, intra_image_comm)
  CALL mp_bcast(filplot, ionode_id, intra_image_comm)
  CALL mp_bcast(lpunch,  ionode_id, intra_image_comm)
  !
  IF (.NOT. lpunch) THEN
    RETURN
  END IF
  !
  ! ... read data of pw.x
  needwf = .FALSE.
  CALL read_file_new(needwf)
  !
  ! ... check condition
  IF (.NOT. lrism3d) THEN
    CALL errore('postrism', 'no data about 3D-RISM', 1)
  END IF
  !
  ! ... punch data
  CALL punch_rism(filplot)
  !
END SUBROUTINE extract_rism
