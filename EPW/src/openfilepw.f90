  !                                                                            
  ! Copyright (C) 2010-2016 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino 
  ! Copyright (C) 2007-2009 Jesse Noffsinger, Brad Malone, Feliciano Giustino  
  !                                                                            
  ! This file is distributed under the terms of the GNU General Public         
  ! License. See the file `LICENSE' in the root directory of the               
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .             
  !                                                                            
  ! Adapted from the code PH/openfilq - Quantum-ESPRESSO group                
  !-----------------------------------------------------------------------
  subroutine openfilepw
  !-----------------------------------------------------------------------
  !
  !     This subroutine opens all the files necessary for the EPW
  !     calculation.
  !
  ! RM - Nov/Dec 2014
  ! Imported the noncolinear case implemented by xlzhang
  !
  !-----------------------------------------------------------------------
  use mp_global,        ONLY : me_pool
  use io_files,         ONLY : prefix, tmp_dir, &
                               diropn, seqopn
  use units_ph,         ONLY : iudrhous, lrdrhous, iudvkb3, iuwfc
  USE uspp,             ONLY : nkb, okvan
  ! nwordwfc is the record length for the direct-access file containing
  ! wavefunctions
  USE wvfct,            ONLY : nbnd, npwx
  USE noncollin_module, ONLY : npol,nspin_mag
  use phcom,            ONLY : lrwfc, fildvscf, iudvscf, lrdvkb3, lrdrho
  use epwcom,           ONLY : elinterp, &
                               iuncuf, nbndsub, lrcuf
  USE fft_base,         ONLY : dfftp
  !USE epwcom,           ONLY : iudvscf0, tphases, fildvscf0
  !
  implicit none
  INTEGER, EXTERNAL :: find_free_unit
  integer :: spot
  ! integer variable for I/O control
  ! used for extracting the fildvscf0 directory
  character (len=256) :: filint,tmp_dir_save
  ! the name of the file
  ! safe place for the tmp_dir variable  
  logical :: exst
  ! logical variable to check file existe
  ! logical variable to check file exists in memory
  !
  IF (len_trim(prefix) == 0) call errore ('openfilepw', 'wrong prefix', 1)
  !
  !     The file with the wavefunctions
  !
  iuwfc = 20 
  lrwfc = 2 * nbnd * npwx * npol 
  CALL diropn(iuwfc,'wfc',lrwfc,exst) 
  IF (.not. exst) CALL errore ('openfilepw','file '//TRIM( prefix )//'.wfc'//' not found',1)
  !
  !
#ifdef __PARA
      IF (me_pool /= 0) goto 300 
#endif
#ifdef __PARA
300  continue
#endif
  !
  !   file for setting unitary gauges of eigenstates
  !
  !
  ! RM - nspin corresponds to nspin_mag according to QE5.0.3
  !    - this will have to change when we move to QE5.0.3 
  !
  !lrdrho = 2 * dfftp%nr1x *dfftp%nr2x *dfftp%nr3x * nspin_mag
  !IF (fildvscf0 .eq. fildvscf) THEN
  !   iudvscf0 = iudvscf
  !ELSE
  !   iudvscf0 = find_free_unit()
  !   IF ( me_pool == 0 .and. tphases) THEN
  !      tmp_dir_save=tmp_dir 
  !      spot=INDEX(fildvscf0,'/',.true.)
  !      tmp_dir=fildvscf0(1:spot) 
  !      CALL diropn (iudvscf0, 'dvscf', lrdrho, exst)
  !      tmp_dir=tmp_dir_save
  !   END IF
  !ENDIF
  !
  !
  IF (elinterp) then
    !
    !  open the file for writing the rotation matrix and the
    !  electron-phonon matrix on the fine mesh (too big to stay
    !  in memory for BC53). Use direct access because the nested
    !  loops on modes and k points. @ FG
    !
    !  currently not used.  Should probably add this back in for big systems
    !  The size is about nbnd^2
    !
    iuncuf    = find_free_unit()
    IF (nbndsub .ne. 0) THEN
       lrcuf  = 2 * nbndsub * nbndsub
    ELSE
       lrcuf  = 2 * nbnd * nbnd
    ENDIF
    filint    = trim(prefix)//'.cuf'
    CALL diropn (iuncuf, 'cuf', lrcuf, exst)  
    !
  ENDIF
  !
  !
  !    In the USPP case we also need a file in  order to store derivatives 
  !    of kb projectors
  !  
  IF (okvan) THEN
     iudvkb3 = find_free_unit()
     lrdvkb3 = 2 * npwx * nkb * 3
     CALL diropn (iudvkb3, 'dvkb3', lrdvkb3, exst)
     !
     iudrhous = find_free_unit()
     ! RM - nspin corresponds to nspin_mag according to QE5.0.3
     !    - this will have to change when we move to QE5.0.3 
     lrdrhous = 2 * dfftp%nnr * nspin_mag
     CALL diropn (iudrhous, 'prd', lrdrhous, exst)
  ENDIF
  !
  RETURN
  !
  end subroutine openfilepw
