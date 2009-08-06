!
! Copyright (C) 2009 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
SUBROUTINE new_occ()
  !-----------------------------------------------------------------------
  !
  ! This routine computes the occupations of the bands according to
  ! their projections on the initial atomic wavefunctios. It is
  ! used in isolated atoms.
  !
  USE io_global,            ONLY : stdout
  USE kinds,                ONLY : DP
  USE basis,                ONLY : natomwfc
  USE klist,                ONLY : nks, ngk
  USE ldaU,                 ONLY : swfcatom
  USE lsda_mod,             ONLY : lsda, current_spin, nspin, isk
  USE wvfct,                ONLY : nbnd, npw, npwx, igk, wg
  USE control_flags,        ONLY : gamma_only
  USE wavefunctions_module, ONLY : evc
  USE noncollin_module,     ONLY : noncolin
  USE gvect,                ONLY : gstart
  USE io_files,             ONLY : iunigk, nwordwfc, iunwfc, nwordatwfc, iunsat
  USE buffers,              ONLY : get_buffer
  USE fixed_occ,            ONLY : one_atom_occupations, f_inp
  USE mp_global,            ONLY : intra_pool_comm
  USE mp,                   ONLY : mp_sum


  IMPLICIT NONE
  !
  ! I/O variables
  !
  INTEGER :: ik, ibnd, i
  ! ounter on k points
  !    "    "  bands

  REAL(DP), EXTERNAL :: ddot
  COMPLEX(DP) :: zdotc
  COMPLEX(DP) , ALLOCATABLE :: proj(:,:)

  REAL(DP) :: max_value, save_value
  INTEGER :: select_ibnd, iatwfc, first_available_band

  ALLOCATE( proj(natomwfc,nbnd))  
  !
  !    we start a loop on k points
  !

  IF (nks > 1) REWIND (iunigk)

  DO ik = 1, nks
     IF (lsda) current_spin = isk(ik)
     npw = ngk (ik)
     IF (nks > 1) THEN
        READ (iunigk) igk
        CALL get_buffer  (evc, nwordwfc, iunwfc, ik)
     END IF
     CALL davcio (swfcatom, nwordatwfc, iunsat, ik, - 1)
     !
     ! make the projection
     !
     DO ibnd = 1, nbnd
        DO i = 1, natomwfc
           IF ( gamma_only ) THEN 
              proj (i, ibnd) = 2.d0 * &
                    ddot(2*npw, swfcatom (1, i), 1, evc (1, ibnd), 1) 
              IF (gstart.EQ.2) proj (i, ibnd) = proj (i, ibnd) - &
                    swfcatom (1, i) * evc (1, ibnd)
           ELSE 
              proj (i, ibnd) = zdotc (npw, swfcatom (1, i), 1, evc (1, ibnd), 1)
              IF (noncolin) &
                 proj (i, ibnd) = proj(i, ibnd) + &
                   zdotc (npw, swfcatom(npwx+1,i), 1, evc(npwx+1,ibnd), 1)
           ENDIF
        ENDDO
     ENDDO

#ifdef __PARA
     CALL mp_sum ( proj, intra_pool_comm )
#endif

     IF (one_atom_occupations) THEN
        IF (nbnd > natomwfc) CALL errore('new_occ','too many bands',1)
!        DO ibnd=1,nbnd
!           write(6,*) 'bands ', ibnd
!           write(6,'(8f9.4)') (ABS(proj(i,ibnd)), i=1,nbnd)
!        END DO
        wg(:,ik)=0.0_DP
        first_available_band=1
        DO iatwfc=1,nbnd
           IF (f_inp(iatwfc,ik)>0.0_DP) THEN
              max_value=ABS(proj(iatwfc,first_available_band))
              select_ibnd=first_available_band
              DO ibnd=nbnd,first_available_band+1,-1
!
!   Search the band with maximum projection
!
                 IF (wg(ibnd,ik)==0.0_DP) THEN
                    save_value=ABS(proj(iatwfc,ibnd))
                    IF (save_value> max_value) THEN
                        max_value=save_value
                        select_ibnd=ibnd
                    END IF
                 ENDIF
              END DO
              wg(select_ibnd,ik)=f_inp(iatwfc,ik)
              DO ibnd=1, nbnd
                 IF (wg(ibnd,ik)==0.0_DP) THEN
                    first_available_band=ibnd
                    EXIT
                 ENDIF
              ENDDO
           END IF      
        END DO
        IF (nks==1) THEN
           WRITE( stdout, '(5x, "Occupations" )' ) 
        ELSE
           IF (ik==1) WRITE( stdout, '(5x, "Occupations for spin-up" )' ) 
           IF (ik==2) WRITE( stdout, '(5x, "Occupations for spin-down" )' ) 
        ENDIF
        WRITE(6,9030) (wg(ibnd,ik), ibnd=1,nbnd)
     END IF
  ENDDO
  DEALLOCATE ( proj )
9030 FORMAT( '  ',8F9.4 )

  RETURN

END SUBROUTINE new_occ
