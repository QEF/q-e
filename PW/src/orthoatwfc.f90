!
! Copyright (C) 2001-2007 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
SUBROUTINE orthoatwfc
  !-----------------------------------------------------------------------
  !
  ! This routine is meant to orthogonalize all the atomic wfcs. This is
  ! useful when we want to compute the occupation of the atomic orbitals
  ! in order to make lda+U calculations
  !
  USE kinds,      ONLY : DP
  USE io_global,  ONLY : stdout
  USE io_files,   ONLY : iunat, iunsat, nwordatwfc, iunigk
  USE ions_base,  ONLY : nat
  USE basis,      ONLY : natomwfc
  USE klist,      ONLY : nks, xk, ngk
  USE ldaU,       ONLY : swfcatom, U_projection
  USE wvfct,      ONLY : npwx, npw, igk
  USE uspp,       ONLY : nkb, vkb
  USE becmod,     ONLY : allocate_bec_type, deallocate_bec_type, &
                         bec_type, becp, calbec
  USE mp_global,  ONLY : intra_pool_comm, intra_bgrp_comm
  USE mp,         ONLY : mp_sum
  USE control_flags,    ONLY : gamma_only
  USE noncollin_module, ONLY : noncolin, npol
  ! 
  IMPLICIT NONE
  !
  !
  INTEGER :: ik, ibnd, info, i, j, k, na, nb, nt, isym, n, ntemp, m, &
       l, lm, ltot, ntot, ipol
  ! the k point under consideration
  ! counter on bands
  REAL(DP) :: t0, scnds
  ! cpu time spent
  LOGICAL :: orthogonalize_wfc
     
  COMPLEX(DP) :: temp, t (5)
  COMPLEX(DP) , ALLOCATABLE :: wfcatom (:,:), work (:,:), overlap (:,:)
  REAL(DP) , ALLOCATABLE :: e (:)

  t0 = scnds ()
  
  IF (noncolin) THEN
     ALLOCATE (wfcatom( npwx*npol, natomwfc))    
  ELSE
     ALLOCATE (wfcatom( npwx, natomwfc))    
  END IF
  ALLOCATE (overlap( natomwfc , natomwfc))    
  ALLOCATE (work   ( natomwfc , natomwfc))    
  ALLOCATE (e      ( natomwfc))    

  IF (U_projection=="file") THEN
     WRITE( stdout,*) 'LDA+U Projector read from file '
     RETURN
  END IF

  IF (U_projection=="atomic") THEN
     orthogonalize_wfc = .FALSE.
     WRITE( stdout,*) 'Atomic wfc used for LDA+U Projector are NOT orthogonalized'
  ELSE IF (U_projection=="ortho-atomic") THEN
     orthogonalize_wfc = .TRUE.
     WRITE( stdout,*) 'Atomic wfc used for LDA+U Projector are orthogonalized'
     IF (gamma_only) THEN
        WRITE( stdout,*) 'Gamma-only calculation for this case not implemented'
        STOP
     END IF
  ELSE IF (U_projection=="norm-atomic") THEN
     orthogonalize_wfc = .TRUE.
     WRITE( stdout,*) 'Atomic wfc used for LDA+U Projector are normalized but NOT orthogonalized'
     IF (gamma_only) THEN
        WRITE( stdout,*) 'Gamma-only calculation for this case not implemented'
        STOP
     END IF
  ELSE
     WRITE( stdout,*) "U_projection_type =", U_projection
     CALL errore ("orthoatwfc"," this U_projection_type is not valid",1)
  END IF

  ! Allocate the array becp = <beta|wfcatom>
  CALL allocate_bec_type (nkb,natomwfc, becp) 
  
  IF (nks > 1) REWIND (iunigk)
  
  DO ik = 1, nks
     
     npw = ngk (ik)
     IF (nks > 1) READ (iunigk) igk
     
     overlap(:,:) = (0.d0,0.d0)
     work(:,:) = (0.d0,0.d0)
     
     IF (noncolin) THEN
       CALL atomic_wfc_nc_updown (ik, wfcatom)
     ELSE
       CALL atomic_wfc (ik, wfcatom)
     ENDIF

     !
     ! write atomic wfc on unit iunat
     !
     CALL davcio (wfcatom, nwordatwfc, iunat, ik, 1)
     CALL init_us_2 (npw, igk, xk (1, ik), vkb)
     CALL calbec (npw, vkb, wfcatom, becp) 
     CALL s_psi (npwx, npw, natomwfc, wfcatom, swfcatom)

     IF (orthogonalize_wfc) THEN
        !
        ! calculate overlap matrix
        !
        IF (noncolin) THEN
           CALL zgemm ('c', 'n', natomwfc, natomwfc, npwx*npol, (1.d0, 0.d0), &
             wfcatom, npwx*npol, swfcatom, npwx*npol, (0.d0,0.d0), overlap, &
                                                                   natomwfc)
        ELSE
           CALL zgemm ('c', 'n', natomwfc, natomwfc, npw, (1.d0, 0.d0), &
             wfcatom, npwx, swfcatom, npwx, (0.d0, 0.d0), overlap, natomwfc)
        END IF
        !
        CALL mp_sum(  overlap, intra_bgrp_comm )
        !
        IF (U_projection=="norm-atomic") THEN
           DO i = 1, natomwfc
              DO j = i+1, natomwfc
                 overlap(i,j) = CMPLX(0.d0,0.d0, kind=dp)
                 overlap(j,i) = CMPLX(0.d0,0.d0, kind=dp)
              ENDDO
           ENDDO
        END IF
        !
        ! find O^-.5
        !
        CALL cdiagh (natomwfc, overlap, natomwfc, e, work)
        DO i = 1, natomwfc
           e (i) = 1.d0 / dsqrt (e (i) )
        ENDDO
        DO i = 1, natomwfc
           DO j = i, natomwfc
              temp = (0.d0, 0.d0)
              DO k = 1, natomwfc
                 temp = temp + e (k) * work (j, k) * CONJG (work (i, k) )
              ENDDO
              overlap (i, j) = temp
              IF (j.NE.i) overlap (j, i) = CONJG (temp)
           ENDDO
        ENDDO
        !
        ! trasform atomic orbitals O^-.5 psi
        !
        DO i = 1, npw
           work(:,1) = (0.d0,0.d0)
           IF (noncolin) THEN
              DO ipol=1,npol
                 j = i + (ipol-1)*npwx
                 CALL zgemv ('n',natomwfc,natomwfc,(1.d0,0.d0),overlap, &
                      natomwfc,swfcatom(j,1),npwx*npol, (0.d0,0.d0),work,1)
                 CALL zcopy (natomwfc,work,1,swfcatom(j,1),npwx*npol)
              END DO
           ELSE
              CALL zgemv ('n', natomwfc, natomwfc, (1.d0, 0.d0) , overlap, &
                   natomwfc, swfcatom (i, 1) , npwx, (0.d0, 0.d0) , work, 1)
              CALL zcopy (natomwfc, work, 1, swfcatom (i, 1), npwx)
           END IF
        ENDDO
        
     END IF ! orthogonalize_wfc
     !
     ! write S * atomic wfc to unit iunsat
     !
     CALL davcio (swfcatom, nwordatwfc, iunsat, ik, 1)
     
  ENDDO
  DEALLOCATE (overlap)
  DEALLOCATE (work)
  DEALLOCATE (e)
  DEALLOCATE (wfcatom)
  CALL deallocate_bec_type ( becp )
  !
  RETURN
     
END SUBROUTINE orthoatwfc
