!
! Copyright (C) 2001-2020 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
#if !defined(__CUDA)
#define cublasZgemm zgemm
#endif
!
!-----------------------------------------------------------------------
SUBROUTINE orthoUwfc_gpu
  !-----------------------------------------------------------------------
  !
  ! This routine saves to buffer "iunhub" atomic wavefunctions having an
  ! associated Hubbard U term * S, for DFT+U(+V) calculations. Same for 
  ! "iunhub2" but without S (this is then used to computed Hubbard forces 
  ! and stresses). Atomic wavefunctions
  ! are orthogonalized if desired, depending upon the value of "U_projection"
  ! "swfcatom" must NOT be allocated on input.
  !
  USE kinds,      ONLY : DP
  USE buffers,    ONLY : get_buffer, save_buffer
  USE io_global,  ONLY : stdout
  USE io_files,   ONLY : iunhub, iunhub2, nwordwfcU
  USE ions_base,  ONLY : nat
  USE basis,      ONLY : natomwfc, swfcatom
  USE klist,      ONLY : nks, xk, ngk, igk_k, igk_k_d
  USE ldaU,       ONLY : U_projection, wfcU, nwfcU, copy_U_wfc
  USE wvfct,      ONLY : npwx
  USE uspp,       ONLY : nkb, vkb
  USE becmod,     ONLY : allocate_bec_type, deallocate_bec_type, &
                         bec_type, becp, calbec
  USE control_flags,    ONLY : gamma_only
  USE noncollin_module, ONLY : noncolin, npol
  USE mp_bands,         ONLY : use_bgrp_in_hpsi
  !
  USE uspp_gpum,        ONLY : using_vkb, using_vkb_d, vkb_d
  USE becmod_gpum,      ONLY : becp_d
  USE becmod_subs_gpum, ONLY : using_becp_auto, using_becp_d_auto, &
                               calbec_gpu
  ! 
  IMPLICIT NONE
  !
  !
  INTEGER :: ik, ibnd, info, i, j, k, na, nb, nt, isym, n, ntemp, m, &
       l, lm, ltot, ntot, ipol, npw
  ! ik: the k point under consideration
  ! ibnd: counter on bands
  LOGICAL :: orthogonalize_wfc, normalize_only, save_flag
  COMPLEX(DP) , ALLOCATABLE :: wfcatom (:,:)
  !
  COMPLEX(DP) , ALLOCATABLE :: wfcatom_d(:,:)
  COMPLEX(DP) , ALLOCATABLE :: swfcatom_d(:,:)
  !
#if defined(__CUDA)
  attributes(DEVICE) :: wfcatom_d, swfcatom_d
#endif  
  !
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
        CALL get_buffer(wfcU, nwordwfcU, iunhub, ik)
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
  !
  ALLOCATE( wfcatom(npwx*npol,natomwfc), swfcatom(npwx*npol,natomwfc) )
  ALLOCATE( wfcatom_d(npwx*npol,natomwfc), swfcatom_d(npwx*npol,natomwfc) )               !
  save_flag = use_bgrp_in_hpsi ; use_bgrp_in_hpsi=.false.
  !
  ! Allocate the array becp = <beta|wfcatom>
  CALL allocate_bec_type( nkb, natomwfc, becp ) 
  CALL using_becp_auto(2)
  !
  !
  DO ik = 1, nks
     !
     IF (noncolin) THEN
       CALL atomic_wfc_nc_updown( ik, wfcatom )
       wfcatom_d = wfcatom
     ELSE
       CALL atomic_wfc_gpu( ik, wfcatom_d )
     ENDIF
     !
     npw = ngk(ik)
     !CALL using_vkb(1)
     CALL using_vkb_d(2)
     CALL init_us_2_gpu( npw, igk_k_d(1,ik), xk(1,ik), vkb_d )
     !
     CALL using_becp_d_auto(2)
     CALL calbec_gpu( npw, vkb_d, wfcatom_d, becp_d )
     !
     CALL s_psi_gpu( npwx, npw, natomwfc, wfcatom_d, swfcatom_d )
     !
     IF (orthogonalize_wfc) CALL ortho_swfc_gpu( npw, normalize_only, &
                                 natomwfc, wfcatom_d, swfcatom_d, .TRUE. )
     !
     ! copy atomic wavefunctions with Hubbard U term only in wfcU
     ! (this is then used to compute Hubbard forces and stresses)
     ! save to unit iunhub2
     !
     wfcatom  = wfcatom_d
     swfcatom = swfcatom_d
     !
     CALL copy_U_wfc( wfcatom, noncolin )
     CALL save_buffer( wfcU, nwordwfcU, iunhub2, ik )
     !
     ! copy S * atomic wavefunctions with Hubbard U term only in wfcU
     ! (this is used during the self-consistent solution of Kohn-Sham equations)
     ! save to unit iunhub
     !
     CALL copy_U_wfc( swfcatom, noncolin )
     IF ( nks > 1 ) &
          CALL save_buffer( wfcU, nwordwfcU, iunhub, ik )
     !
  ENDDO
  !
  DEALLOCATE( wfcatom, swfcatom )
  DEALLOCATE( wfcatom_d, swfcatom_d )
  !
  CALL deallocate_bec_type( becp )
  CALL using_becp_auto(2)
  !
  use_bgrp_in_hpsi = save_flag
  !
  RETURN
     
END SUBROUTINE orthoUwfc_gpu
!
!-----------------------------------------------------------------------
SUBROUTINE orthoatwfc_gpu( orthogonalize_wfc )
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
  USE klist,      ONLY : nks, xk, ngk, igk_k, igk_k_d
  USE wvfct,      ONLY : npwx
  USE uspp,       ONLY : nkb, vkb
  USE becmod,     ONLY : allocate_bec_type, deallocate_bec_type, &
                         bec_type, becp, calbec
  USE control_flags,    ONLY : gamma_only
  USE noncollin_module, ONLY : noncolin, npol
  !
  USE uspp_gpum,        ONLY : using_vkb, using_vkb_d, vkb_d
  USE becmod_gpum,      ONLY : becp_d
  USE becmod_subs_gpum, ONLY : using_becp_auto, using_becp_d_auto, &
                               calbec_gpu
  ! 
  IMPLICIT NONE
  !
  LOGICAL, INTENT(IN) :: orthogonalize_wfc
  !
  INTEGER :: ik, ibnd, info, i, j, k, na, nb, nt, isym, n, ntemp, m, &
             l, lm, ltot, ntot, ipol, npw
  ! ik: the k point under consideration
  ! ibnd: counter on bands
  LOGICAL :: normalize_only = .FALSE.
  COMPLEX(DP), ALLOCATABLE :: wfcatom(:,:)
  COMPLEX(DP), ALLOCATABLE :: wfcatom_d(:,:), swfcatom_d(:,:)
  !
#if defined(__CUDA)
  attributes(DEVICE) :: wfcatom_d, swfcatom_d
#endif
  !
  normalize_only=.FALSE.
  ALLOCATE( wfcatom(npwx*npol,natomwfc) )
  ALLOCATE( wfcatom_d(npwx*npol,natomwfc), swfcatom_d(npwx*npol,natomwfc) )
  !
  ! Allocate the array becp = <beta|wfcatom>
  CALL allocate_bec_type( nkb, natomwfc, becp )
  !
  
  DO ik = 1, nks
     !
     IF (noncolin) THEN
       CALL atomic_wfc_nc_updown( ik, wfcatom )
       wfcatom_d = wfcatom
     ELSE
       CALL atomic_wfc_gpu( ik, wfcatom_d )
     ENDIF
     !
     npw = ngk(ik)
     !CALL using_vkb(1)
     CALL using_vkb_d(2)
     CALL init_us_2_gpu( npw, igk_k_d(1,ik), xk(1,ik), vkb_d )
     !
     CALL using_becp_auto(2)
     CALL using_becp_d_auto(2)
     CALL calbec_gpu( npw, vkb_d, wfcatom_d, becp_d )
     !
     CALL s_psi_gpu( npwx, npw, natomwfc, wfcatom_d, swfcatom_d )
     !
     IF (orthogonalize_wfc) CALL ortho_swfc_gpu( npw, normalize_only, &
                                 natomwfc, wfcatom_d, swfcatom_d, .FALSE. )
     !
     ! write S * atomic wfc to unit iunsat
     !
     swfcatom = swfcatom_d
     !
     CALL save_buffer( swfcatom, nwordatwfc, iunsat, ik )
     !
  ENDDO
  !
  DEALLOCATE( wfcatom )
  DEALLOCATE( wfcatom_d, swfcatom_d )
  CALL deallocate_bec_type( becp )
  !
  RETURN
  !
END SUBROUTINE orthoatwfc_gpu
!
!-----------------------------------------------------------------------
SUBROUTINE ortho_swfc_gpu( npw, normalize_only, m, wfc_d, swfc_d, lflag )
  !-----------------------------------------------------------------------
  !
  ! On input : wfc (npwx*npol,m) =  \psi = a set of "m" (atomic) wavefcts
  !            swfc(npwx*npol,m) = S\psi 
  !            normalize_only    = only normalize, do not orthonormalize
  !
  ! This routine will compute the overlap matrix O: 
  ! O_ij = <wfc_i|S|wfc_j> = <wfc_i|swfc_j>
  !
  ! On output: swfc = O^{-1/2} S\psi, i.e. S * orthonormalized wavefunctions
  ! If lflag=.FALSE. : wfc are unchanged on output (not orthonormalized), i.e.
  !                    wfc = \psi
  ! If lflag=.TRUE.  : wfc are orthonormalized on output, i.e.
  !                    wfc = O^{-1/2} \psi, <wfc_i|S|wfc_j> = \delta_{ij}
  !
#if defined(__CUDA)
  USE cublas
#endif
  !
  USE kinds,            ONLY : DP
  USE wvfct,            ONLY : npwx
  USE mp_bands,         ONLY : intra_bgrp_comm, me_bgrp, root_bgrp
  USE mp,               ONLY : mp_sum
  USE noncollin_module, ONLY : noncolin, npol
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: m, npw
  LOGICAL, INTENT(IN) :: normalize_only
  COMPLEX(DP), INTENT(INOUT) :: wfc_d(npwx*npol,m)
  COMPLEX(DP), INTENT(INOUT) :: swfc_d(npwx*npol,m)
  LOGICAL, INTENT(IN) :: lflag
  !
  ! ... local variables
  !
  INTEGER :: i, j, k, ipol
  COMPLEX(DP) :: temp 
  !
  COMPLEX(DP), ALLOCATABLE :: overlap_d(:,:)
  COMPLEX(DP), ALLOCATABLE :: work_d(:,:), s_d(:,:)
  REAL(DP) , ALLOCATABLE :: e_d(:)
  !
#if defined(__CUDA)
  attributes(DEVICE) :: wfc_d, swfc_d, overlap_d, work_d, s_d, e_d
#endif  
  !
  ALLOCATE( overlap_d(m,m) )
  ALLOCATE( work_d(m,m), s_d(m,m) )
  ALLOCATE( e_d(m) )
  !
  overlap_d(:,:) = (0.d0,0.d0)
  work_d(:,:) = (0.d0,0.d0)
  !
  ! calculate overlap matrix
  
  !
  IF (noncolin) THEN
     CALL cublasZgemm( 'C', 'N', m, m, npwx*npol, (1.d0,0.d0), wfc_d, &
                      npwx*npol, swfc_d, npwx*npol, (0.d0,0.d0), overlap_d, m )
  ELSE
     CALL cublasZgemm( 'C', 'N', m, m, npw, (1.d0,0.d0), wfc_d, &
                      npwx, swfc_d, npwx, (0.d0,0.d0), overlap_d, m )
  END IF
  !
  CALL mp_sum( overlap_d, intra_bgrp_comm )
  !
  IF ( normalize_only ) THEN
     !$cuf kernel do (1) <<<*,*>>>
     DO i = 1, m
        DO j = i+1, m
           overlap_d(i,j) = CMPLX(0.d0,0.d0, kind=dp)
           overlap_d(j,i) = CMPLX(0.d0,0.d0, kind=dp)
        ENDDO
     ENDDO
  END IF
  !
  ! find O^(-1/2)
  !
  s_d = CMPLX(0.d0,0.d0, kind=dp)
  !$cuf kernel do (1)
  DO i = 1, m
    s_d(i,i) = CMPLX(1.d0,0.d0, kind=dp)
  ENDDO
  !
  CALL laxlib_cdiaghg_gpu( m, m, overlap_d, s_d, m, e_d, work_d, me_bgrp, &
                           root_bgrp, intra_bgrp_comm )
  !
  !$cuf kernel do (1) <<<*,*>>>
  DO i = 1, m
     e_d(i) = 1.d0 / SQRT(e_d(i))
  ENDDO
  !$cuf kernel do (1) <<<*,*>>>
  DO i = 1, m
     DO j = i, m
        temp = (0.d0, 0.d0)
        DO k = 1, m
           temp = temp + e_d(k) * work_d(j,k) * CONJG(work_d(i,k))
        ENDDO
        overlap_d(i,j) = temp
        IF (j /= i) overlap_d(j,i) = CONJG(temp)
     ENDDO
  ENDDO
  !
  ! transform atomic orbitals O^(-1/2) S\psi
  ! FIXME: can be done in a faster way by using wfc as work space 
  !
  DO i = 1, npw
     work_d(:,1) = (0.d0,0.d0)
     IF (noncolin) THEN
        DO ipol=1,npol
           j = i + (ipol-1)*npwx
           CALL cublasZgemv( 'N', m, m, (1.d0,0.d0), overlap_d, &
                             m, swfc_d(j,1), npwx*npol, (0.d0,0.d0), work_d, 1 )
           CALL cublasZcopy( m, work_d, 1, swfc_d(j,1), npwx*npol )
        END DO
     ELSE
        CALL cublasZgemv( 'N', m, m, (1.d0,0.d0), overlap_d, &
                          m, swfc_d(i,1), npwx, (0.d0,0.d0), work_d, 1 )
        CALL cublasZcopy( m, work_d, 1, swfc_d(i,1), npwx )
     END IF
  ENDDO
  !
  ! If lflag=.TRUE. transform atomic orbitals without
  ! the ultrasoft S operator O^(-1/2) \psi
  !
  IF (lflag) THEN
   DO i = 1, npw
     work_d(:,1) = (0.d0,0.d0)
     IF (noncolin) THEN
        DO ipol = 1, npol
           j = i + (ipol-1)*npwx
           CALL cublasZgemv( 'N', m, m, (1.d0,0.d0), overlap_d, &
                             m, wfc_d(j,1), npwx*npol, (0.d0,0.d0), work_d, 1 )
           CALL cublasZcopy( m, work_d, 1, wfc_d(j,1), npwx*npol )
        END DO
     ELSE
        CALL cublasZgemv( 'N', m, m, (1.d0,0.d0), overlap_d, &
                          m, wfc_d(i,1), npwx, (0.d0,0.d0), work_d, 1 )
        CALL cublasZcopy( m, work_d, 1, wfc_d(i,1), npwx )
     END IF
   ENDDO
  ENDIF
  !
  DEALLOCATE( overlap_d, s_d )
  DEALLOCATE( work_d, e_d )
  !
  RETURN
  !      
END SUBROUTINE ortho_swfc_gpu
