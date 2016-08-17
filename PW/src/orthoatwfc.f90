!
! Copyright (C) 2001-2013 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
SUBROUTINE orthoUwfc
  !-----------------------------------------------------------------------
  !
  ! This routine saves to buffer "iunhub" atomic wavefunctions having an
  ! associated Hubbard U term, for DFT+U calculations. Atomic wavefunctions
  ! are orthogonalized if desired, depending upon the value of "U_projection"
  ! "swfcatom" must NOT be allocated on input.
  !
  USE kinds,      ONLY : DP
  USE buffers,    ONLY : get_buffer, save_buffer
  USE io_global,  ONLY : stdout
  USE io_files,   ONLY : iunhub, nwordwfcU
  USE ions_base,  ONLY : nat
  USE basis,      ONLY : natomwfc, swfcatom
  USE klist,      ONLY : nks, xk, ngk, igk_k
  USE ldaU,       ONLY : U_projection, wfcU, nwfcU, copy_U_wfc
  USE wvfct,      ONLY : npwx
  USE uspp,       ONLY : nkb, vkb
  USE becmod,     ONLY : allocate_bec_type, deallocate_bec_type, &
                         bec_type, becp, calbec
  USE control_flags,    ONLY : gamma_only
  USE noncollin_module, ONLY : noncolin, npol
  ! 
  IMPLICIT NONE
  !
  !
  INTEGER :: ik, ibnd, info, i, j, k, na, nb, nt, isym, n, ntemp, m, &
       l, lm, ltot, ntot, ipol, npw
  ! ik: the k point under consideration
  ! ibnd: counter on bands
  LOGICAL :: orthogonalize_wfc, normalize_only
  COMPLEX(DP) , ALLOCATABLE :: wfcatom (:,:)

  IF ( U_projection == "pseudo" ) THEN
     WRITE( stdout,*) 'Beta functions used for LDA+U Projector'
     RETURN
  ELSE IF (U_projection=="file") THEN
     !
     ! Read atomic wavefunctions from file (produced by pmw.x). In this case,
     ! U-specific atomic wavefunctions wfcU coincide with atomic wavefunctions 
     !
     WRITE( stdout,*) 'LDA+U Projector read from file '
     DO ik = 1, nks
        CALL get_buffer (wfcU, nwordwfcU, iunhub, ik)
     END DO
     RETURN
  ELSE IF (U_projection=="atomic") THEN
     orthogonalize_wfc = .FALSE.
     normalize_only = .FALSE.
     WRITE( stdout,*) 'Atomic wfc used for LDA+U Projector are NOT orthogonalized'
  ELSE IF (U_projection=="ortho-atomic") THEN
     orthogonalize_wfc = .TRUE.
     normalize_only = .FALSE.    
     WRITE( stdout,*) 'Atomic wfc used for LDA+U Projector are orthogonalized'
     IF (gamma_only) CALL errore('orthoatwfc', &
          'Gamma-only calculation for this case not implemented', 1 )
  ELSE IF (U_projection=="norm-atomic") THEN
     orthogonalize_wfc = .TRUE.
     normalize_only = .TRUE.
     WRITE( stdout,*) 'Atomic wfc used for LDA+U Projector are normalized but NOT orthogonalized'
     IF (gamma_only) CALL errore('orthoatwfc', &
          'Gamma-only calculation for this case not implemented', 1 )
  ELSE
     WRITE( stdout,*) "U_projection_type =", U_projection
     CALL errore ("orthoatwfc"," this U_projection_type is not valid",1)
  END IF

  ALLOCATE ( wfcatom(npwx*npol, natomwfc), swfcatom(npwx*npol, natomwfc) )

  ! Allocate the array becp = <beta|wfcatom>
  CALL allocate_bec_type (nkb,natomwfc, becp) 
  
  DO ik = 1, nks
     
     IF (noncolin) THEN
       CALL atomic_wfc_nc_updown (ik, wfcatom)
     ELSE
       CALL atomic_wfc (ik, wfcatom)
     ENDIF
     npw = ngk (ik)
     CALL init_us_2 (npw, igk_k(1,ik), xk (1, ik), vkb)
     CALL calbec (npw, vkb, wfcatom, becp) 
     CALL s_psi (npwx, npw, natomwfc, wfcatom, swfcatom)

     IF (orthogonalize_wfc) &
        CALL ortho_swfc ( normalize_only, natomwfc, wfcatom, swfcatom )
     !
     ! copy atomic wavefunctions with Hubbard U term only in wfcU
     ! save to unit iunhub
     !
     CALL copy_U_wfc (swfcatom, noncolin)
     IF ( nks > 1 ) &
          CALL save_buffer (wfcU, nwordwfcU, iunhub, ik)
     !
  ENDDO
  DEALLOCATE (wfcatom, swfcatom)
  CALL deallocate_bec_type ( becp )
  !
  RETURN
     
END SUBROUTINE orthoUwfc
!
!-----------------------------------------------------------------------
SUBROUTINE orthoatwfc (orthogonalize_wfc)
  !-----------------------------------------------------------------------
  !
  ! This routine calculates atomic wavefunctions, orthogonalizes them
  ! if "orthogonalzie_wfc" is .true., saves them into buffer "iunsat".
  ! "swfcatom" must be allocated on input.
  ! Useful for options "wannier" and "one_atom_occupations"
  !
  USE kinds,      ONLY : DP
  USE buffers,    ONLY : save_buffer
  USE io_global,  ONLY : stdout
  USE io_files,   ONLY : iunsat, nwordatwfc
  USE ions_base,  ONLY : nat
  USE basis,      ONLY : natomwfc, swfcatom
  USE klist,      ONLY : nks, xk, ngk, igk_k
  USE wvfct,      ONLY : npwx
  USE uspp,       ONLY : nkb, vkb
  USE becmod,     ONLY : allocate_bec_type, deallocate_bec_type, &
                         bec_type, becp, calbec
  USE control_flags,    ONLY : gamma_only
  USE noncollin_module, ONLY : noncolin, npol
  ! 
  IMPLICIT NONE
  !
  LOGICAL, INTENT(in) :: orthogonalize_wfc
  !
  INTEGER :: ik, ibnd, info, i, j, k, na, nb, nt, isym, n, ntemp, m, &
       l, lm, ltot, ntot, ipol, npw
  ! ik: the k point under consideration
  ! ibnd: counter on bands
  LOGICAL :: normalize_only = .FALSE.
  COMPLEX(DP) , ALLOCATABLE :: wfcatom (:,:)

  normalize_only=.FALSE.
  ALLOCATE (wfcatom( npwx*npol, natomwfc))

  ! Allocate the array becp = <beta|wfcatom>
  CALL allocate_bec_type (nkb,natomwfc, becp) 
  
  DO ik = 1, nks
     
     IF (noncolin) THEN
       CALL atomic_wfc_nc_updown (ik, wfcatom)
     ELSE
       CALL atomic_wfc (ik, wfcatom)
     ENDIF
     npw = ngk (ik)
     CALL init_us_2 (npw, igk_k(1,ik), xk (1, ik), vkb)
     CALL calbec (npw, vkb, wfcatom, becp) 
     CALL s_psi (npwx, npw, natomwfc, wfcatom, swfcatom)

     IF (orthogonalize_wfc) &
        CALL ortho_swfc ( normalize_only, natomwfc, wfcatom, swfcatom )
     !
     ! write S * atomic wfc to unit iunsat
     !
     CALL save_buffer (swfcatom, nwordatwfc, iunsat, ik)
     !
  ENDDO
  DEALLOCATE (wfcatom)
  CALL deallocate_bec_type ( becp )
  !
  RETURN
     
END SUBROUTINE orthoatwfc
!
!-----------------------------------------------------------------------
SUBROUTINE ortho_swfc ( normalize_only, m, wfc, swfc )
  !-----------------------------------------------------------------------
  !
  ! On input : wfc (npwx*npol,m) =  \psi = a set of "m" (atomic) wavefcts
  !            swfc(npwx*npol,m) = S\psi 
  !            normalize_only    = only normalize, do not orthonormalize
  ! On output: swfc = S^{-1/2}\psi = orthonormalized wavefunctions
  !                                 (i.e. <swfc_i|S|swfc_j> = \delta_{ij})
  !             wfc = currently unchanged
  !
  USE kinds,      ONLY : DP
  USE wvfct,      ONLY : npwx, npw
  USE mp_bands,   ONLY : intra_bgrp_comm
  USE mp,         ONLY : mp_sum
  USE noncollin_module, ONLY : noncolin, npol
  IMPLICIT NONE
  !
  INTEGER, INTENT(in) :: m
  LOGICAL, INTENT(in) :: normalize_only
  COMPLEX(dp), INTENT(IN   ) :: wfc (npwx*npol,m)
  COMPLEX(dp), INTENT(INOUT) :: swfc(npwx*npol,m)

  COMPLEX(DP) :: temp 
  COMPLEX(DP) , ALLOCATABLE ::  work (:,:), overlap (:,:)
  REAL(DP) , ALLOCATABLE :: e (:)
  INTEGER :: i, j, k, ipol

  ALLOCATE (overlap( m , m))    
  ALLOCATE (work   ( m , m))    
  ALLOCATE (e      ( m))    
  
  overlap(:,:) = (0.d0,0.d0)
  work(:,:) = (0.d0,0.d0)
  !
  ! calculate overlap matrix
  !
  IF (noncolin) THEN
     CALL zgemm ('c', 'n', m, m, npwx*npol, (1.d0, 0.d0), wfc, &
          npwx*npol, swfc, npwx*npol, (0.d0,0.d0), overlap, m)
  ELSE
     CALL zgemm ('c', 'n', m, m, npw, (1.d0, 0.d0), wfc, &
          npwx, swfc, npwx, (0.d0, 0.d0), overlap, m)
  END IF
  !
  CALL mp_sum(  overlap, intra_bgrp_comm )
  !
  IF ( normalize_only ) THEN
     DO i = 1, m
        DO j = i+1, m
           overlap(i,j) = CMPLX(0.d0,0.d0, kind=dp)
           overlap(j,i) = CMPLX(0.d0,0.d0, kind=dp)
        ENDDO
     ENDDO
  END IF
  !
  ! find O^-.5
  !
  CALL cdiagh (m, overlap, m, e, work)
  DO i = 1, m
     e (i) = 1.d0 / SQRT (e (i) )
  ENDDO
  DO i = 1, m
     DO j = i, m
        temp = (0.d0, 0.d0)
        DO k = 1, m
           temp = temp + e (k) * work (j, k) * CONJG (work (i, k) )
        ENDDO
        overlap (i, j) = temp
        IF (j.NE.i) overlap (j, i) = CONJG (temp)
     ENDDO
  ENDDO
  !
  ! trasform atomic orbitals O^-.5 psi
  ! FIXME: can be done in a faster way by using wfc as work space 
  !
  DO i = 1, npw
     work(:,1) = (0.d0,0.d0)
     IF (noncolin) THEN
        DO ipol=1,npol
           j = i + (ipol-1)*npwx
           CALL zgemv ('n',m,m,(1.d0,0.d0),overlap, &
                m, swfc(j,1),npwx*npol, (0.d0,0.d0),work,1)
           CALL zcopy (m,work,1,swfc(j,1),npwx*npol)
        END DO
     ELSE
        CALL zgemv ('n', m, m, (1.d0, 0.d0) , overlap, &
             m, swfc (i, 1) , npwx, (0.d0, 0.d0) , work, 1)
        CALL zcopy (m, work, 1, swfc (i, 1), npwx)
     END IF
  ENDDO
  
  DEALLOCATE (overlap)
  DEALLOCATE (work)
  DEALLOCATE (e)
        
END SUBROUTINE ortho_swfc
