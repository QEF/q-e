!
! Copyright (C) 2001-2004 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
SUBROUTINE openfilq()
  !----------------------------------------------------------------------------
  !
  ! ... This subroutine opens all the files necessary for the phononq
  ! ... calculation.
  !
  USE kinds,          ONLY : DP
  USE units_ph,       ONLY : iuwfc, iudwf, iubar, iucom, iudvkb3, &
                             iudrhous, iuebar, iudrho, iudyn, iudvscf, &
                             lrwfc, lrdwf, lrbar, lrcom, lrdvkb3, &
                             lrdrhous, lrebar, lrdrho
  USE control_ph,     ONLY : epsil, zue, recover, trans, elph
  USE output,         ONLY : fildrho, fildyn, fildvscf
  USE us,             ONLY : okvan
  USE wvfct,          ONLY : nbnd, npwx
  USE gvect,          ONLY : nrx1, nrx2, nrx3, nrxx
  USE lsda_mod,       ONLY : nspin
  USE uspp,           ONLY : nkb
  USE io_files,       ONLY : prefix, iunigk
  USE control_flags,  ONLY : twfcollect
  USE restart_module, ONLY : readfile_new
  USE mp_global,      ONLY : me_pool, kunit
  USE io_global,      ONLY : ionode
  !
  IMPLICIT NONE
  !
  INTEGER :: ios
  ! integer variable for I/O control
  CHARACTER (len=256) :: filint
  ! the name of the file
  LOGICAL :: exst
  ! logical variable to check file existe
  !
  REAL(kind=DP) :: edum(1,1), wdum(1,1)
  INTEGER :: ndr, ierr, kunittmp

  IF (LEN_TRIM(prefix) == 0) CALL errore ('openfilq', 'wrong prefix', 1)

  !
  !     There are six direct access files to be opened in the tmp area
  !
  !     The file with the wavefunctions
  !
  iuwfc = 20

  lrwfc = 2 * nbnd * npwx
  filint = TRIM(prefix)//'.wfc'
  CALL diropn (iuwfc, filint, lrwfc, exst)

  IF (.NOT.exst) THEN
    ndr      = 4
    kunittmp = 1

    kunittmp = kunit

    CALL readfile_new( 'wave', ndr, edum, wdum, kunittmp, lrwfc, iuwfc, ierr )
    IF( ierr > 0 ) &
      CALL errore ('openfilq', 'file '//filint//' not found', 1)
    twfcollect=.NOT.exst
  END IF
  !
  !    The file with deltaV_{bare} * psi
  !
  iubar = 21
  lrbar = 2 * nbnd * npwx
  filint = TRIM(prefix) //'.bar'
  CALL diropn (iubar, filint, lrbar, exst)
  IF (recover.AND..NOT.exst) CALL errore ('openfilq','file bar not found', 1)
  !
  !    The file with the solution delta psi
  !
  iudwf = 22
  lrdwf = 2 * nbnd * npwx
  filint = TRIM(prefix) //'.dwf'
  CALL diropn (iudwf, filint, lrdwf, exst)
  IF (recover.AND..NOT.exst) CALL errore ('openfilq','file dwf not found', 1)
  !
  !   open a file with the static change of the charge
  !
  IF (okvan) THEN
     iudrhous = 25
     lrdrhous = 2 * nrxx * nspin
     filint = TRIM(prefix) //'.prd'
     CALL diropn (iudrhous, filint, lrdrhous, exst)
     IF (recover.AND..NOT.exst) CALL errore ('openfilq','file prod not found', 1)
  ENDIF
  !
  !   An optional file for testing purposes containing the deltarho
  !
  lrdrho = 2 * nrx1 * nrx2 * nrx3 * nspin
  !
  IF (fildrho.NE.' ') THEN
     iudrho = 23
     IF ( me_pool == 0 ) THEN
        !
        IF(epsil) THEN
           filint = TRIM(fildrho)//".E"
        ELSE
           filint = TRIM(fildrho)//".u"
        END IF
        !
        CALL diropn (iudrho, filint, lrdrho, exst)
        !
     END IF

  ENDIF
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

  !   ... by the first node only, other nodes write on unit 6 (i.e./dev/nu
  !   exception: electron-phonon calculation from saved data
  !  (iudyn is read, not written, by all nodes)
  !
  IF ( ( .NOT. ionode ) .AND. (.NOT.elph.OR.trans) ) THEN
     iudyn = 6
     GOTO 400
  ENDIF

  IF (trans.OR.elph) THEN
     iudyn = 26
     OPEN (unit = iudyn, file = fildyn, status = 'unknown', err = &
          100, iostat = ios)
100  CALL errore ('openfilq', 'opening file'//fildyn, ABS (ios) )
     REWIND (iudyn)
  ENDIF
  !
  !   An optional file for electron-phonon calculations containing deltaVs
  !
400 IF (fildvscf.NE.' ') THEN
     iudvscf = 27
     !!! CALL seqopn (iudvscf, fildvscf, 'unformatted', exst)
     !!! REWIND (iudvscf)
     IF ( me_pool == 0 ) THEN
        CALL diropn (iudvscf, fildvscf, lrdrho, exst)
     END IF
  END IF
  !
  !    In the USPP case we need two files for the Commutator, the first is
  !    given by filbar and a second which just contains P_c x |psi>,
  !    which is required for the calculation of the Born effective carges
  !
  IF (okvan .AND. (epsil .OR. zue)) THEN
     iucom = 28
     lrcom = 2 * nbnd * npwx
     filint = TRIM(prefix) //'.com'
     CALL diropn (iucom, filint, lrcom, exst)
     IF (recover.AND..NOT.exst) &
         CALL errore ('openfilq', 'file com not found', 1)
  !
  !    In the USPP case we also need a file in  order to store derivatives 
  !    ok kb projectors
  !  
     iudvkb3 = 29
     lrdvkb3 = 2 * npwx * nkb * 3
     filint = TRIM(prefix) //'.dvkb3'
     CALL diropn (iudvkb3, filint, lrdvkb3, exst)
     IF (recover.AND..NOT.exst) &
         CALL errore ('openfilq', 'file dvkb3 not found', 1)
  ENDIF
  IF (epsil .OR. zue) THEN
     iuebar = 30
     lrebar = 2 * nbnd * npwx
     filint = TRIM(prefix) //'.ebar'
     CALL diropn (iuebar, filint, lrebar, exst)
     IF (recover.AND..NOT.exst) &
        CALL errore ('openfilq','file ebar not found', 1)
  ENDIF
  !
  RETURN
  !
END SUBROUTINE openfilq
