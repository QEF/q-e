!
! Copyright (C) 2003-2021 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
SUBROUTINE read_hr()
  !---------------------------------------------------------------------
  !
  ! ...  This routine calculates KC matrix H(R) - if not calculated
  ! ...  already - and it prints it into a formatted file
  !
  USE control_kcw,          ONLY : Hamlt_R, num_wann, irvect, nkstot_eff
  USE klist,                ONLY : nkstot
  USE lsda_mod,             ONLY : nspin
  USE interpolation,        ONLY : real_ham
  USE io_files,             ONLY : prefix
  USE io_global,            ONLY : ionode
  USE constants,            ONLY : rytoev
  USE mp_global,            ONLY : intra_image_comm
  USE mp,                   ONLY : mp_bcast
  USE io_global,            ONLY : ionode, ionode_id
  !
  !
  IMPLICIT NONE
  !
  CHARACTER(LEN=33) :: header
  INTEGER :: ir
  INTEGER :: iband, jband, num_r, idum, ipol
  INTEGER :: ierr
  LOGICAL :: exst
  !
  !
  IF (nspin == 4) THEN
    nkstot_eff = nkstot
  ELSE
    nkstot_eff = nkstot/nspin
  ENDIF
  !
  !
  !
  IF ( ionode ) THEN
    !
    INQUIRE( file=TRIM(prefix)//'.kcw_hr.dat', exist=exst )
    !WRITE(*,*) "HERE", exst
    IF ( .NOT. exst ) call errore('read_hr', 'File H(R) not FOUND', 1 )
    
    OPEN( 100, file=TRIM(prefix)//'.kcw_hr.dat', status='unknown', IOSTAT=ierr )
    IF (ierr /= 0 ) call errore('read_hr', 'Error while opening H(R) file ', abs (ierr) )
    !
    READ( 100, *) header
    READ( 100, '(10x,i5,16x,i5)' ) num_r, num_wann
    !WRITE(*,*) num_r, nkstot
    IF (num_r /= nkstot_eff) & 
      CALL errore('read_hr', 'Number of R/k point DOES not MATCH',num_r)
  ENDIF
  !
  CALL mp_bcast( num_r, ionode_id, intra_image_comm )
  CALL mp_bcast( num_wann, ionode_id, intra_image_comm )
  ALLOCATE( Hamlt_R(num_r,num_wann,num_wann), irvect(3,num_r)) 
  !
  IF (ionode ) THEN
    !
    DO ir = 1, num_r
      DO iband = 1, num_wann
        DO jband = 1, num_wann
          !
          READ( 100, '(5i5,2f12.6)' ) (irvect(ipol,ir), ipol=1,3), idum, idum, Hamlt_R(ir,jband,iband)
          Hamlt_R(ir,jband,iband) = Hamlt_R(ir,jband,iband)/rytoev
          !
        ENDDO
      ENDDO
    ENDDO
    !
    CLOSE( 100 )
    !
  ENDIF
  !
  CALL mp_bcast( Hamlt_R, ionode_id, intra_image_comm )
  !
  !
END SUBROUTINE read_hr
