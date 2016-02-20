!
! Copyright (C) 2001-2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
SUBROUTINE openfild3
  !-----------------------------------------------------------------------
  !
  !     This subroutine opens all the files necessary for the
  !     third derivative calculation.
  !
  USE pwcom
  USE phcom
  USE d3com
  USE fft_base,        ONLY : dfftp
  USE control_flags,   ONLY : twfcollect
  USE io_files,        ONLY : iunigk, prefix, tmp_dir, diropn, seqopn
  USE io_global,       ONLY : ionode
  USE mp_global,       ONLY : kunit, me_pool, root_pool
  USE uspp,            ONLY : nlcc_any

  USE control_lr, ONLY : lgamma
  !
  IMPLICIT NONE
  !
  INTEGER :: ios
  ! integer variable for I/O control
  CHARACTER (len=256) :: filint, tmp_dir_save
  ! the name of the file
  LOGICAL :: exst
  ! logical variable to check file existe
  INTEGER       :: ndr, kunittmp, ierr
  REAL(DP) :: edum(1,1), wdum(1,1)

  twfcollect=.FALSE.

  IF (LEN_TRIM(prefix) == 0) CALL errore ('openfild3', 'wrong prefix', 1)
  !
  !     The file with the wavefunctions
  !
  iuwfc = 20

  lrwfc = 2 * nbnd * npwx
  CALL diropn (iuwfc, 'wfc', lrwfc, exst)
  IF (.NOT.exst) THEN
     CALL errore ('openfild3', 'file ' // TRIM(prefix) //'.wfc not found', 1)
  END IF
  !
  !    The file with deltaV_{bare} * psi
  !
  iubar = 21
  lrbar = 2 * nbnd * npwx
  CALL diropn (iubar, 'bar', lrbar, exst)
  IF (recover.AND..NOT.exst) &
     CALL errore ('openfild3', 'file ' // TRIM(prefix) //'.bar not found', 1)
  !
  !    The file with the solution delta psi
  !
  iudwf = 22
  lrdwf = 2 * nbnd * npwx
  CALL diropn (iudwf, 'dwf', lrdwf, exst)
  IF (recover.AND..NOT.exst) &
     CALL errore ('openfild3', 'file ' // TRIM(prefix) //'.dwf not found', 1)
  !
  !   Here the sequential files
  !
  !   The igk at a given k (and k+q if q!=0)
  !
  iunigk = 24
  CALL seqopn (iunigk, 'igk', 'unformatted', exst)
  !
  !   a formatted file which contains the dynamical matrix in cartesian
  !   coordinates is opened in the current directory
  !
  ! ... by the first node only, other nodes write on unit 6 (i.e. /dev/null)
  !
  IF ( ionode ) THEN
     !
     iudyn = 26
     OPEN (unit=iudyn, file=fildyn, status='unknown', err=110, iostat=ios)
110  CALL errore ('openfild3', 'opening file'//fildyn, ABS (ios) )
     REWIND (iudyn)
     !
  ELSE
     !
     iudyn = 6
     !
  END IF
  !cccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !
  ! Variation of the charge density with respect to a perturbation
  ! with a generic q
  !
  iudrho = 25
  iud0rho = 33
  IF (lgamma) iud0rho = iudrho
  lrdrho = 2 * dfftp%nr1x * dfftp%nr2x * dfftp%nr3x * nspin
  !
  !   is opened only by the first task of each pool
  !
  IF ( me_pool == root_pool ) THEN
     !
     filint = TRIM(fildrho) !//".u"
     ! FIXME: workaround for filename mess
     tmp_dir_save=tmp_dir
     if ( lgamma) tmp_dir=TRIM(tmp_dir)//'_ph0/'
     !
     CALL diropn (iudrho, filint, lrdrho, exst)
     IF(nlcc_any) CALL diropn (iudrho+1000, trim(filint)//"_cc", lrdrho, exst)
     !
     tmp_dir=tmp_dir_save
     ! FIXME END
     !
     ! Variation of the charge density with respect to a perturbation with q=
     ! Not needed if q=0
     !
     IF (.NOT.lgamma) THEN
        filint = TRIM(fild0rho) !//".u"
        CALL diropn (iud0rho, filint, lrdrho, exst)
        IF(nlcc_any) CALL diropn (iud0rho+1000, trim(filint)//"_cc", lrdrho, exst)
     ENDIF
     !
  END IF
  !
  ! If q=0, we need only one file with the variation of the wavefunctions
  !
  iud0qwf = iudwf
  iudqwf = iudwf
  IF (.NOT.lgamma) THEN
     !
     !    Open the file with the solution q=0 delta psi
     !
     iud0qwf = 34
     CALL diropn (iud0qwf, 'd0wf', lrdwf, exst)
     !
     !    Open the file with the solution q=0 delta psi
     !
     iudqwf = 35
     CALL diropn (iudqwf, 'dqwf', lrdwf, exst)
  ENDIF
  !
  !    The file with   <psi| dqV |psi>
  !
  iupdqvp = 36
  lrpdqvp = 2 * nbnd * nbnd
  CALL diropn (iupdqvp, 'pdp' , lrpdqvp, exst)
  !
  !    The file with   <psi| d0V |psi>
  !
  iupd0vp = iupdqvp
  IF (.NOT.lgamma) THEN
     iupd0vp = 37
     CALL diropn (iupd0vp, 'p0p', lrpdqvp, exst)

  ENDIF
  IF (degauss.NE.0.d0) THEN
     !
     !    The file with   <dqpsi| dqV |psi> (only in the metallic case)
     !
     iudpdvp_1 = 38
     lrdpdvp = 2 * nbnd * nbnd
     CALL diropn (iudpdvp_1, 'pv1' , lrdpdvp, exst)
     !
     !    The file with   <dqpsi| d0V |psi>
     !
     iudpdvp_2 = iudpdvp_1
     iudpdvp_3 = iudpdvp_1
     IF (.NOT.lgamma) THEN
        iudpdvp_2 = 39
        CALL diropn (iudpdvp_2, 'pv2' , lrdpdvp, exst)
        !
        !    The file with   <d0psi| dqV |psi>
        !
        iudpdvp_3 = 40
        CALL diropn (iudpdvp_3, 'pv3', lrdpdvp, exst)
     ENDIF
     !
     ! The file containing the variation of the FermiEnergy ef_sh
     !
     ! opened only by the first task of the first pool
     !

     IF ( ionode ) THEN
        !
        iuef = 41
        CALL seqopn (iuef, 'efs', 'unformatted', exst)
        !
     END IF
     !
  ENDIF
  RETURN
END SUBROUTINE openfild3
