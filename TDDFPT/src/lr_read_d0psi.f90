!
! Copyright (C) 2001-2015 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
SUBROUTINE lr_read_d0psi()
  !---------------------------------------------------------------------
  !
  ! This subroutine reads in and stores vectors necessary to
  ! restart the Lanczos recursion.
  !
  ! Modified by Osman Baris Malcioglu (2009)
  !
  USE klist,                ONLY : nks,degauss
  USE io_files,             ONLY : prefix, diropn, tmp_dir, wfc_dir
  USE lr_variables,         ONLY : d0psi, n_ipol, LR_polarization, lr_verbosity, &
                                 & nwordd0psi, iund0psi, eels
  USE wvfct,                ONLY : nbnd, npwx, et
  USE io_global,            ONLY : stdout
  !
  IMPLICIT NONE
  !
  ! local variables
  !
  INTEGER :: ip
  CHARACTER(len=6), EXTERNAL :: int_to_char
  LOGICAL :: exst
  CHARACTER(len=256) :: tmp_dir_saved
  !
  If (lr_verbosity > 5) THEN
    WRITE(stdout,'("<lr_read_d0psi>")')
  endif
  !
  ! This is a parallel read, done in wfc_dir
  !
  tmp_dir_saved = tmp_dir
  !
  IF ( wfc_dir /= 'undefined' ) tmp_dir = wfc_dir
  !
  DO ip=1,n_ipol
     !
     IF (n_ipol==1) THEN
        !
        CALL diropn ( iund0psi, 'd0psi.'//trim(int_to_char(LR_polarization)), nwordd0psi, exst)
        !
        IF (.not.exst .and. wfc_dir /= 'undefined') THEN
           !
           WRITE( stdout, '(/5x,"Attempting to read d0psi from outdir instead of wfcdir")' )
           CLOSE( UNIT = iund0psi)
           tmp_dir = tmp_dir_saved 
           !
           CALL diropn ( iund0psi, 'd0psi.'//trim(int_to_char(LR_polarization)), nwordd0psi, exst)
           !
           IF (.not.exst) CALL errore('lr_read_d0psi', &
                  & trim( prefix )//'.d0psi.'//trim(int_to_char(LR_polarization))//' not found',1)
           !
        ENDIF
        !
     ENDIF
     !
     IF (n_ipol==3 .and. .not.eels) THEN
        !
        CALL diropn ( iund0psi, 'd0psi.'//trim(int_to_char(ip)), nwordd0psi, exst)
        !
        IF (.not.exst .and. wfc_dir /= 'undefined') THEN
           ! 
           WRITE( stdout, '(/5x,"Attempting to read d0psi from outdir instead of wfcdir")' )
           CLOSE( UNIT = iund0psi)
           tmp_dir = tmp_dir_saved
           !
           CALL diropn ( iund0psi, 'd0psi.'//trim(int_to_char(LR_polarization)), nwordd0psi, exst)
           !
           IF (.not.exst) CALL errore('lr_read_d0psi', &
                  & trim( prefix )//'.d0psi.'//trim(int_to_char(ip))//' not found',1)
           !
        ENDIF
        !
     ENDIF
     !
     CALL davcio(d0psi(1,1,1,ip),nwordd0psi,iund0psi,1,-1)
     !
     CLOSE( UNIT = iund0psi)
     !
  ENDDO
  !
  ! End of file i/o
  !
  tmp_dir = tmp_dir_saved
  !
END SUBROUTINE lr_read_d0psi
!-----------------------------------------------------------------------
