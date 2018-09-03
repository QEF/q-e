!
! Copyright (C) 2001-2018 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
SUBROUTINE hp_sphi
  !-----------------------------------------------------------------------
  !
  ! This routine computes and writes on file the atomic functions
  ! S(k)*phi(k) and S(k+q)*phi(k+q). In the NCPP case S(k)=1 and S(k+q)=1. 
  ! phi(k) is defined by Eq. (A6) in Ref. [1]
  ! [1] Phys. Rev. B 98, 085127 (2018)
  !
  USE kinds,            ONLY : DP
  USE io_files,         ONLY : nwordatwfc
  USE basis,            ONLY : natomwfc
  USE klist,            ONLY : xk, ngk, igk_k
  USE wvfct,            ONLY : npwx
  USE uspp,             ONLY : vkb, nkb, okvan
  USE becmod,           ONLY : allocate_bec_type, deallocate_bec_type, becp
  USE buffers,          ONLY : save_buffer
  USE io_global,        ONLY : stdout
  USE mp,               ONLY : mp_sum
  USE mp_global,        ONLY : intra_pool_comm
  USE noncollin_module, ONLY : noncolin
  USE eqv,              ONLY : swfcatomk, swfcatomkpq
  USE qpoint,           ONLY : nksq, ikks, ikqs
  USE control_lr,       ONLY : lgamma
  USE units_lr,         ONLY : iuatwfc
  USE ldaU,             ONLY : U_projection
  ! 
  IMPLICIT NONE
  INTEGER :: ik, ikk, ikq
  COMPLEX(DP), ALLOCATABLE :: wfcatom(:,:)  ! atomic wfc
  !
  CALL start_clock ('hp_sphi')
  !
  IF (okvan) CALL allocate_bec_type (nkb, natomwfc, becp)
  !
  ALLOCATE (wfcatom(npwx, natomwfc))
  !
  IF (U_projection=="atomic") THEN
     WRITE( stdout, '(/5x,"Atomic wfc used for the DFT+U projector are NOT orthogonalized")')
  ELSEIF (U_projection=="ortho-atomic") THEN
     WRITE( stdout, '(/5x,"Atomic wfc used for the DFT+U projector are orthogonalized")')
  ELSEIF (U_projection=="norm-atomic") THEN
     WRITE( stdout, '(/5x,"Atomic wfc used for the DFT+U projector are normalized but NOT orthogonalized")')
  ELSE
     WRITE(stdout,*) "U_projection_type =", U_projection
     CALL errore ("hp_sphi"," This U_projection_type is not valid",1)
  ENDIF
  !
  DO ik = 1, nksq
     !
     ikk = ikks(ik)
     ikq = ikqs(ik)
     !
     ! Determine the atomic orbital at k : phi(k)
     !
     CALL atomic_wfc (ikk, wfcatom)
     !
     ! Compute S(k)*phi(k) (phi means the atomic orbital)
     !
     CALL s_phi (ikk, wfcatom, swfcatomk)     
     !
     ! Orthogonalize the atomic orbitals (if needed) 
     !
     IF ( U_projection=="ortho-atomic" .OR. U_projection=="norm-atomic") &
          CALL hp_orthogonalize_or_normalize_wfc (ikk, wfcatom, swfcatomk)
     !
     ! Write S(k)*phi(k) on file
     !
     CALL save_buffer (swfcatomk, nwordatwfc, iuatwfc, ikk)
     !
     IF (.NOT.lgamma) THEN
        !
        ! Determine the atomic orbital at k+q : phi(k+q)
        ! 
        CALL atomic_wfc (ikq, wfcatom)
        !
        ! Compute S(k+q)*phi(k+q) 
        !
        CALL s_phi (ikq, wfcatom, swfcatomkpq)
        !
        ! Orthogonalize the atomic orbitals (if needed)
        !
        IF ( U_projection=="ortho-atomic" .OR. U_projection=="norm-atomic") &
             CALL hp_orthogonalize_or_normalize_wfc (ikq, wfcatom, swfcatomkpq)
        !
        ! Write S(k+q)*phi(k+q) on file
        !
        CALL save_buffer (swfcatomkpq, nwordatwfc, iuatwfc, ikq)
        !
     ENDIF
     !
  ENDDO  
  !
  DEALLOCATE(wfcatom)
  ! 
  IF (okvan) CALL deallocate_bec_type ( becp )
  !
  CALL stop_clock ('hp_sphi')
  !
  RETURN
  !
CONTAINS
  !
SUBROUTINE s_phi (ik_, wfc, swfc)
  !-----------------------------------------------------------------------
  !
  ! NCPP: swfc = wfc
  ! USPP: swfc = S * wfc
  !
  USE kinds,          ONLY : DP
  USE becmod,         ONLY : calbec
  !
  IMPLICIT NONE
  INTEGER,     INTENT(IN)  :: ik_
  COMPLEX(DP), INTENT(IN)  :: wfc  (npwx, natomwfc)
  COMPLEX(DP), INTENT(OUT) :: swfc (npwx, natomwfc)
  INTEGER :: npw_
  !
  ! NCPP case
  !
  IF ( nkb == 0 .OR. .NOT. okvan ) THEN
     swfc = wfc
     RETURN
  ENDIF
  !
  ! USPP case
  !
  ! Number of plane waves at ik_
  !
  npw_ = ngk(ik_)
  !
  ! Compute beta functions vkb at ik_
  !
  CALL init_us_2 (npw_, igk_k(1,ik_), xk(1,ik_), vkb)
  !
  ! Compute the product of beta functions vkb
  ! with the functions wfc : becp = <vkb|wfc>
  !
  CALL calbec (npw_, vkb, wfc, becp)
  !
  ! Calculate S*|wfc> = |wfc> + \sum qq * |vkb> * becp 
  !  
  CALL s_psi (npwx, npw_, natomwfc, wfc, swfc)
  !
  RETURN
  !
END SUBROUTINE s_phi
  !
SUBROUTINE hp_orthogonalize_or_normalize_wfc (ik_, wfc, swfc)
  !
  ! This routine:
  ! - orthogonalizes atomic orbitals if U_projection=="ortho-atomic" 
  !   or
  ! - normalizes atomic orbitals if U_projection=="norm-atomic" 
  !
  ! Inspired by PW/src/orthoatwfc.f90
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: ik_
  COMPLEX(DP), INTENT(IN)    :: wfc  (npwx, natomwfc)  
  COMPLEX(DP), INTENT(INOUT) :: swfc (npwx, natomwfc)  
  !
  ! Local variables
  !
  INTEGER :: npw_
  COMPLEX(DP), ALLOCATABLE :: work(:,:), overlap(:,:)
  REAL(DP),    ALLOCATABLE :: e(:)
  COMPLEX(DP) :: temp
  INTEGER :: i, j, k
  !
  IF (noncolin) CALL errore('hp_orthogonalize_or_normalize_wfc', &
                     & 'Noncolliner case is not implemented',1)
  !
  ! Number of plane waves at ik_
  !
  npw_ = ngk(ik_)
  !
  ALLOCATE (overlap(natomwfc, natomwfc))
  ALLOCATE (work   (natomwfc, natomwfc))
  ALLOCATE (e      (natomwfc))
  !
  overlap(:,:) = (0.d0,0.d0)
  work(:,:)    = (0.d0,0.d0)
  !
  ! Calculate the overlap matrix
  !
  CALL ZGEMM ('C', 'N', natomwfc, natomwfc, npw_, (1.d0, 0.d0), &
             wfc, npwx, swfc, npwx, (0.d0, 0.d0), overlap, natomwfc)  
  !
  CALL mp_sum(overlap, intra_pool_comm)
  !
  IF (U_projection=="norm-atomic") THEN
     ! Set to zero the off-diagonal matrix elements
     DO i = 1, natomwfc
        DO j = i+1, natomwfc
           overlap(i,j) = (0.d0, 0.d0)
           overlap(j,i) = (0.d0, 0.d0)
        ENDDO
     ENDDO
  ENDIF 
  !
  ! Find O^-.5
  !
  CALL cdiagh (natomwfc, overlap, natomwfc, e, work)
  !  
  DO i = 1, natomwfc
     e(i) = 1.d0 / DSQRT(e(i))
  ENDDO
  DO i = 1, natomwfc
     DO j = i, natomwfc
        temp = (0.d0, 0.d0)
        DO k = 1, natomwfc
           temp = temp + e(k) * work(j,k) * CONJG(work(i,k))
        ENDDO
        overlap(i, j) = temp
        IF (j.NE.i) overlap(j, i) = CONJG(temp)
     ENDDO
  ENDDO
  !
  ! Trasform atomic orbitals O^-.5 psi
  !
  DO i = 1, npw_
     work(:,1) = (0.d0, 0.d0)
     CALL ZGEMV ('N', natomwfc, natomwfc, (1.d0, 0.d0), overlap, &
                 natomwfc, swfc(i,1), npwx, (0.d0, 0.d0), work, 1)
     CALL ZCOPY (natomwfc, work, 1, swfc(i,1), npwx)
  ENDDO
  !
  RETURN
  !
END SUBROUTINE hp_orthogonalize_or_normalize_wfc
  !
END SUBROUTINE hp_sphi
