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
  USE control_flags,   ONLY : twfcollect
  USE io_files,        ONLY : iunigk, prefix
  USE restart_module,  ONLY : readfile_new
  USE io_global,       ONLY : ionode
  USE mp_global,       ONLY : kunit, me_pool, root_pool
  !
  IMPLICIT NONE
  !
  INTEGER :: ios
  ! integer variable for I/O control
  CHARACTER (len=256) :: filint
  ! the name of the file
  LOGICAL :: exst
  ! logical variable to check file existe
  INTEGER       :: ndr, kunittmp, ierr
  REAL(KIND=DP) :: edum(1,1), wdum(1,1)

  twfcollect=.FALSE.

  IF (LEN_TRIM(prefix) == 0) CALL errore ('openfild3', 'wrong prefix', 1)
  !
  !     The file with the wavefunctions
  !
  iuwfc = 20

  lrwfc = 2 * nbnd * npwx
  filint = TRIM(prefix) //'.wfc'
  CALL diropn (iuwfc, filint, lrwfc, exst)
  IF (.NOT.exst) THEN 
     ndr      = 4
     kunittmp = 1
     kunittmp = kunit
     CALL readfile_new( 'wave', ndr, edum, wdum, kunittmp, lrwfc, &
                        iuwfc, ierr )
     IF ( ierr > 0 ) &
        CALL errore ('openfild3', 'file '//filint//' not found', 1)
     twfcollect=.NOT.exst
  END IF
  !
  !    The file with deltaV_{bare} * psi
  !
  iubar = 21
  lrbar = 2 * nbnd * npwx
  filint = TRIM(prefix) //'.bar'
  CALL diropn (iubar, filint, lrbar, exst)
  IF (recover.AND..NOT.exst) CALL errore ('openfild3', 'file bar not &
       &found', 1)
  !
  !    The file with the solution delta psi
  !
  iudwf = 22
  lrdwf = 2 * nbnd * npwx
  filint = TRIM(prefix) //'.dwf'
  CALL diropn (iudwf, filint, lrdwf, exst)
  IF (recover.AND..NOT.exst) CALL errore ('openfild3', 'file dwf not &
       &found', 1)
  !
  !   Here the sequential files
  !
  !   The igk at a given k (and k+q if q!=0)
  !
  iunigk = 24
  filint = TRIM(prefix) //'.igk'
  CALL seqopn (iunigk, filint, 'unformatted', exst)
  !
  !   a formatted file which contains the dynamical matrix in cartesian
  !   coordinates is opened in the current directory
  !
  ! ... by the first node only, other nodes write on unit 6 (i.e. /dev/null)
  !
  IF ( ionode ) THEN
     !
     iudyn = 26
     OPEN (unit = iudyn, file = fildyn, status = 'unknown', err = 110, &
          iostat = ios)
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
  lrdrho = 2 * nrx1 * nrx2 * nrx3 * nspin
  !
  !   is opened only by the first task of each pool
  !
  IF ( me_pool == root_pool ) THEN
     !
     filint = TRIM(fildrho)//".u"
     CALL diropn (iudrho, filint, lrdrho, exst)
     !
     ! Variation of the charge density with respect to a perturbation with q=
     ! Not needed if q=0
     !
     IF (.NOT.lgamma) THEN
        filint = TRIM(fild0rho)//".u"
        CALL diropn (iud0rho, filint, lrdrho, exst)
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
     filint = TRIM(prefix) //'.d0wf'
     CALL diropn (iud0qwf, filint, lrdwf, exst)
     !
     !    Open the file with the solution q=0 delta psi
     !
     iudqwf = 35
     filint = TRIM(prefix) //'.dqwf'
     CALL diropn (iudqwf, filint, lrdwf, exst)
  ENDIF
  !
  !    The file with   <psi| dqV |psi>
  !
  iupdqvp = 36
  lrpdqvp = 2 * nbnd * nbnd
  filint = TRIM(prefix) //'.pdp'
  CALL diropn (iupdqvp, filint, lrpdqvp, exst)
  !
  !    The file with   <psi| d0V |psi>
  !
  iupd0vp = iupdqvp
  IF (.NOT.lgamma) THEN
     iupd0vp = 37
     filint = TRIM(prefix) //'.p0p'
     CALL diropn (iupd0vp, filint, lrpdqvp, exst)

  ENDIF
  IF (degauss.NE.0.d0) THEN
     !
     !    The file with   <dqpsi| dqV |psi> (only in the metallic case)
     !
     iudpdvp_1 = 38
     lrdpdvp = 2 * nbnd * nbnd
     filint = TRIM(prefix) //'.pv1'
     CALL diropn (iudpdvp_1, filint, lrdpdvp, exst)
     !
     !    The file with   <dqpsi| d0V |psi>
     !
     iudpdvp_2 = iudpdvp_1
     iudpdvp_3 = iudpdvp_1
     IF (.NOT.lgamma) THEN
        iudpdvp_2 = 39
        filint = TRIM(prefix) //'.pv2'
        CALL diropn (iudpdvp_2, filint, lrdpdvp, exst)
        !
        !    The file with   <d0psi| dqV |psi>
        !
        iudpdvp_3 = 40
        filint = TRIM(prefix) //'.pv3'
        CALL diropn (iudpdvp_3, filint, lrdpdvp, exst)
     ENDIF
     !
     ! The file containing the variation of the FermiEnergy ef_sh
     !
     ! opened only by the first task of the first pool
     !

     IF ( ionode ) THEN
        !
        iuef = 41
        filint = TRIM(prefix) //'.efs'
        CALL seqopn (iuef, filint, 'unformatted', exst)
        !
     END IF
     !
  ENDIF
  RETURN
END SUBROUTINE openfild3
