!
! Copyright (C) 2013 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
SUBROUTINE restart_in_electrons (iter, dr2, ethr, et)
  !-----------------------------------------------------------------------
  USE kinds,         ONLY: dp
  USE io_global,     ONLY: stdout
  USE io_files,      ONLY: iunres, seqopn
  USE klist,         ONLY: nks
  USE wvfct,         ONLY: nbnd
  USE add_dmft_occ,  ONLY: dmft
#if defined(__OSCDFT)
  USE plugin_flags,     ONLY : use_oscdft
  USE oscdft_base,      ONLY : oscdft_ctx
  USE oscdft_functions, ONLY : oscdft_electrons_restart
#endif
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT (inout) :: iter
  REAL(dp), INTENT(inout) :: dr2, ethr, et(nbnd,nks)
  !
  REAL(dp), ALLOCATABLE :: et_(:,:)
  REAL(dp):: dr2_, ethr_
  INTEGER :: ios
  LOGICAL :: exst
  !
  CALL seqopn (iunres, 'restart_scf', 'formatted', exst)
  IF ( exst ) THEN
     ios = 0
     READ (iunres, *, iostat=ios) iter, dr2_, ethr_
     IF ( ios /= 0 ) THEN
        iter = 0
     ELSE IF ( iter < 1 .AND. .NOT. dmft) THEN
        iter = 0
     ELSE
        ALLOCATE (et_(nbnd,nks))
        READ (iunres, *, iostat=ios) et_
        IF ( ios /= 0 ) THEN
           iter = 0
        ELSE
           IF (dmft) THEN
               WRITE( stdout, &
               '(5x,"Calculation restarted from a previous run, expecting DMFT hdf5 archive")' )
           ELSE
               WRITE( stdout, &
               '(5x,"Calculation restarted from scf iteration #",i6)' ) iter + 1
           ENDIF
           dr2 = dr2_
           ethr= ethr_
           et (:,:) = et_(:,:)
        END IF
        DEALLOCATE (et_)
     END IF
#if defined(__OSCDFT)
     IF (use_oscdft .AND. (oscdft_ctx%inp%oscdft_type==1)) THEN
        CALL oscdft_electrons_restart(oscdft_ctx)
     END IF
#endif
  ELSE
     iter = 0
  END IF
  CLOSE ( unit=iunres, status='delete')
  !
END SUBROUTINE restart_in_electrons
