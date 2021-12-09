!
! Copyright (C) 2001-2021 Quantum ESPRESSO group
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
  ! Modified by Iurii Timrov (2015)
  !
  USE klist,                ONLY : nks,degauss
  USE io_files,             ONLY : prefix, diropn, tmp_dir, wfc_dir
  USE lr_variables,         ONLY : d0psi, d0psi2, n_ipol, LR_polarization, &
                                   & lr_verbosity, nwordd0psi, iund0psi, eels, &
                                   & magnons, V0psi, O_psi, ipol, n_op
  USE wvfct,                ONLY : nbnd, npwx, et
  USE io_global,            ONLY : stdout
  USE qpoint,               ONLY : nksq
  USE noncollin_module,     ONLY : npol
  USE control_lr,           ONLY : nbnd_occx
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
  nwordd0psi = 2 * nbnd * npwx * npol * nksq
  !  
  ! This is a parallel read, done in wfc_dir
  !
  tmp_dir_saved = tmp_dir
  !
  IF ( wfc_dir /= 'undefined' ) tmp_dir = wfc_dir
  !
  IF (magnons) THEN
     !
     nwordd0psi = 4 * nbnd_occx * npwx * npol * nksq
     !
     ! read V0psi
     !
     V0psi = (0.0d0, 0.0d0)
     !
     DO ip=1,n_ipol
        !        
        IF (n_ipol==1) THEN
           !
           CALL diropn ( iund0psi, 'V0psi.'//trim(int_to_char(ipol)), nwordd0psi, exst)
           !
           IF (.not.exst .and. wfc_dir /= 'undefined') THEN
              !
              WRITE( stdout, '(/5x,"Attempting to read V0psi from outdir instead of wfcdir")' )
              CLOSE( UNIT = iund0psi)
              tmp_dir = tmp_dir_saved 
              !
              CALL diropn ( iund0psi, 'V0psi.'//trim(int_to_char(ipol)), nwordd0psi, exst)
              !
              IF (.not.exst) CALL errore('lr_read_d0psi', &
                     & trim( prefix )//'.V0psi.'//trim(int_to_char(ipol))//' not found',1)
              !
           ENDIF
           !
        ENDIF
        !
        IF (n_ipol==3 ) THEN
           !
           CALL diropn ( iund0psi, 'V0psi.'//trim(int_to_char(ip)), nwordd0psi, exst)
           !
           IF (.not.exst .and. wfc_dir /= 'undefined') THEN
              ! 
              WRITE( stdout, '(/5x,"Attempting to read V0psi from outdir instead of wfcdir")' )
              CLOSE( UNIT = iund0psi)
              tmp_dir = tmp_dir_saved
              !
              CALL diropn ( iund0psi, 'V0psi.'//trim(int_to_char(ip)), nwordd0psi, exst)
              !
              IF (.not.exst) CALL errore('lr_read_d0psi', &
                     & trim( prefix )//'.V0psi.'//trim(int_to_char(ip))//' not found',1)
              !
           ENDIF
           !
        ENDIF
        !
        CALL davcio(V0psi(:,:,:,:,ip),nwordd0psi,iund0psi,1,-1)
        !
        CLOSE( UNIT = iund0psi)
        !
     ENDDO
     !
     ! read O_psi
     !
     O_psi = (0.0d0, 0.0d0)
     !
     DO ip=1,n_op
        !
        IF (n_op==3) THEN
           !
           CALL diropn ( iund0psi, 'O_psi.'//trim(int_to_char(ip)), nwordd0psi, exst)
           !
           IF (.not.exst .and. wfc_dir /= 'undefined') THEN
              ! 
              WRITE( stdout, '(/5x,"Attempting to read O_psi from outdir instead of wfcdir")' )
              CLOSE( UNIT = iund0psi)
              tmp_dir = tmp_dir_saved
              !
              CALL diropn ( iund0psi, 'O_psi.'//trim(int_to_char(ip)), nwordd0psi, exst)
              !
              IF (.not.exst) CALL errore('lr_read_d0psi', &
                     & trim( prefix )//'.O_psi.'//trim(int_to_char(ip))//' not found',1)
              !
           ENDIF
           !
        ENDIF
        !
        CALL davcio(O_psi(:,:,:,:,ip),nwordd0psi,iund0psi,1,-1)
        !
        CLOSE( UNIT = iund0psi)
        !
     ENDDO
  ELSE 
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
  ENDIF
  !
  ! EELS: Reading of d0psi2 (n_ipol=1)
  !
  IF (eels) THEN
     !
     CALL diropn ( iund0psi, 'd0psi2.'//trim(int_to_char(LR_polarization)), nwordd0psi, exst)
     CALL davcio(d0psi2(1,1,1,1),nwordd0psi,iund0psi,1,-1)
     CLOSE( UNIT = iund0psi)
     !
  ENDIF
  !
  ! End of file i/o
  !
  tmp_dir = tmp_dir_saved
  !
  RETURN
  !
END SUBROUTINE lr_read_d0psi
!-----------------------------------------------------------------------
