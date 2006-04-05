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
  USE output,         ONLY : fildyn, fildvscf
  USE wvfct,          ONLY : nbnd, npwx
  USE gvect,          ONLY : nrx1, nrx2, nrx3, nrxx
  USE lsda_mod,       ONLY : nspin
  USE uspp,           ONLY : nkb, okvan
  USE io_files,       ONLY : prefix, iunigk
  USE control_flags,  ONLY : twfcollect
  USE mp_global,      ONLY : me_pool, kunit
  USE io_global,      ONLY : ionode
  USE ramanm, ONLY: lraman, elop, iuchf, iud2w, iuba2, lrchf, lrd2w, lrba2
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
  REAL(DP) :: edum(1,1), wdum(1,1)
  INTEGER :: ndr, ierr, kunittmp
  !
  !
  IF (LEN_TRIM(prefix) == 0) CALL errore ('openfilq', 'wrong prefix', 1)
  !
  !     There are six direct access files to be opened in the tmp area
  !
  !     The file with the wavefunctions
  !
  iuwfc = 20
  lrwfc = 2 * nbnd * npwx
  CALL diropn (iuwfc, 'wfc', lrwfc, exst)
  IF (.NOT.exst) THEN
     CALL errore ('openfilq', 'file '//trim(prefix)//'.wfc not found', 1)
  END IF
  !
  !    The file with deltaV_{bare} * psi
  !
  iubar = 21
  lrbar = 2 * nbnd * npwx
  CALL diropn (iubar, 'bar', lrbar, exst)
  IF (recover.AND..NOT.exst) &
     CALL errore ('openfilq','file '//trim(prefix)//'.bar not found', 1)
  !
  !    The file with the solution delta psi
  !
  iudwf = 22
  lrdwf = 2 * nbnd * npwx
  CALL diropn (iudwf, 'dwf', lrdwf, exst)
  IF (recover.AND..NOT.exst) &
     CALL errore ('openfilq','file '//trim(prefix)//'.dwf not found', 1)
  !
  !   open a file with the static change of the charge
  !
  IF (okvan) THEN
     iudrhous = 25
     lrdrhous = 2 * nrxx * nspin
     CALL diropn (iudrhous, 'prd', lrdrhous, exst)
     IF (recover.AND..NOT.exst) &
        CALL errore ('openfilq','file '//trim(prefix)//'.prd not found', 1)
  ENDIF
  !
  !  Optional file(s) containing Delta\rho (opened and written in solve_e
  !  and solve_linter). Used for third-order calculations.
  !
  iudrho = 23
  lrdrho = 2 * nrx1 * nrx2 * nrx3 * nspin
  !
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

  !   ... by the first node only, other nodes write on unit 6 (i.e./dev/null
  !   exception: electron-phonon calculation from saved data
  !  (iudyn is read, not written, by all nodes)
  !
  IF ( ( .NOT. ionode ) .AND. (.NOT.elph.OR.trans) ) THEN
     iudyn = 6
     GOTO 400
  ENDIF

  IF (trans.OR.elph) THEN
     iudyn = 26
     OPEN (unit=iudyn, file=fildyn, status='unknown', err=100, iostat=ios)
100  CALL errore ('openfilq', 'opening file'//fildyn, ABS (ios) )
     REWIND (iudyn)
  ENDIF
  !
  !   An optional file for electron-phonon calculations containing deltaVscf
  !
400 IF (fildvscf.NE.' ') THEN
     iudvscf = 27
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
     CALL diropn (iucom, 'com', lrcom, exst)
     IF (recover.AND..NOT.exst) &
         CALL errore ('openfilq', 'file '//trim(prefix)//'.com not found', 1)
  !
  !    In the USPP case we also need a file in  order to store derivatives 
  !    of kb projectors
  !  
     iudvkb3 = 29
     lrdvkb3 = 2 * npwx * nkb * 3
     CALL diropn (iudvkb3, 'dvkb3', lrdvkb3, exst)
     IF (recover.AND..NOT.exst) &
         CALL errore ('openfilq', 'file '//trim(prefix)//'.dvkb3 not found', 1)
  ENDIF
  IF (epsil .OR. zue) THEN
     iuebar = 30
     lrebar = 2 * nbnd * npwx
     CALL diropn (iuebar, 'ebar', lrebar, exst)
     IF (recover.AND..NOT.exst) &
        CALL errore ('openfilq','file '//trim(prefix)//'.ebar not found', 1)
  ENDIF
  !
  !    files used by raman calculation
  !
  IF (lraman .OR.elop) THEN
     iuchf = 31
     lrchf = 2 * nbnd * npwx
     CALL diropn (iuchf, 'cwf', lrchf, exst)

     iud2w = 32
     lrd2w = 2 * nbnd * npwx
     CALL diropn (iud2w, 'd2w', lrd2w, exst)

     iuba2 = 33
     lrba2 = 2 * nbnd * npwx
     CALL diropn(iuba2, 'ba2', lrba2, exst)
  ENDIF

  RETURN
  !
END SUBROUTINE openfilq
