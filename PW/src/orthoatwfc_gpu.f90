!
! Copyright (C) 2001-2022 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
#if !defined(__CUDA)
#define cublasZgemm zgemm
#define cublasZgemv zgemv
#define cublasZcopy zcopy
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
  ! are orthogonalized if desired, depending upon the value of "Hubbard_projectors"
  ! "swfcatom" must NOT be allocated on input.
  !
  USE kinds,      ONLY : DP
  USE buffers,    ONLY : get_buffer, save_buffer
  USE io_global,  ONLY : stdout
  USE io_files,   ONLY : iunhub, nwordwfcU
  USE ions_base,  ONLY : nat
  USE basis,      ONLY : natomwfc, swfcatom
  USE klist,      ONLY : nks, xk, ngk, igk_k
  USE ldaU,       ONLY : Hubbard_projectors, wfcU, nwfcU, copy_U_wfc
  USE wvfct,      ONLY : npwx
  USE uspp,       ONLY : nkb, vkb
  USE becmod,     ONLY : allocate_bec_type, deallocate_bec_type, &
                         bec_type, becp, calbec
  USE control_flags,    ONLY : gamma_only
  USE noncollin_module, ONLY : noncolin, npol
  USE mp_bands,         ONLY : use_bgrp_in_hpsi
  !
  USE becmod_gpum,      ONLY : becp_d
  USE becmod_subs_gpum, ONLY : using_becp_auto, using_becp_d_auto, &
                               calbec_gpu
  USE uspp_init,        ONLY : init_us_2
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
  IF ( Hubbard_projectors == "pseudo" ) THEN
     WRITE( stdout,*) 'Beta functions used for Hubbard projectors'
     RETURN
  ELSE IF (Hubbard_projectors=="wf") THEN
     !
     ! Read Wannier functions from file (produced by pmw.x).
     !
     WRITE( stdout,*) 'Hubbard projectors are read from file produced by pmw.x'
     DO ik = 1, nks
        CALL get_buffer (wfcU, nwordwfcU, iunhub, ik)
     END DO
     RETURN
     !
  ELSE IF (Hubbard_projectors=="atomic") THEN
     orthogonalize_wfc = .FALSE.
     normalize_only = .FALSE.
     WRITE( stdout,'(/5x,a,/)') 'Atomic wfc used for Hubbard projectors are NOT orthogonalized'
  ELSE IF (Hubbard_projectors=="ortho-atomic") THEN
     orthogonalize_wfc = .TRUE.
     normalize_only = .FALSE.    
     WRITE( stdout,'(/5x,a,/)') 'Atomic wfc used for Hubbard projectors are orthogonalized'
     IF (gamma_only) CALL errore('orthoUwfc', &
          'Gamma-only calculation for this case not implemented', 1 )
  ELSE IF (Hubbard_projectors=="norm-atomic") THEN
     orthogonalize_wfc = .TRUE.
     normalize_only = .TRUE.
     WRITE( stdout,'(/5x,a,/)') 'Atomic wfc used for Hubbard projectors are normalized but NOT orthogonalized'
     IF (gamma_only) CALL errore('orthoUwfc', &
          'Gamma-only calculation for this case not implemented', 1 )
  ELSE
     WRITE(stdout,'(/5x,"Hubbard_projectors = ",a)') Hubbard_projectors
     CALL errore ("orthoUwfc"," This type of Hubbard projectors is not valid",1)
  END IF
  !
  ALLOCATE( wfcatom(npwx*npol,natomwfc), swfcatom(npwx*npol,natomwfc) )
  ALLOCATE( wfcatom_d(npwx*npol,natomwfc), swfcatom_d(npwx*npol,natomwfc) )               !
  save_flag = use_bgrp_in_hpsi ; use_bgrp_in_hpsi=.false.

  ! Allocate the array becp = <beta|wfcatom>
  CALL allocate_bec_type (nkb,natomwfc, becp) 
  CALL using_becp_auto(2)
  
  DO ik = 1, nks
     
     IF (noncolin) THEN
       CALL atomic_wfc_nc_updown( ik, wfcatom )
       wfcatom_d = wfcatom
     ELSE
       CALL atomic_wfc_gpu( ik, wfcatom_d )
     ENDIF
     !
     npw = ngk(ik)
     !
     CALL init_us_2( npw, igk_k(1,ik), xk(1,ik), vkb, .true. )
     !
     CALL using_becp_d_auto(2)
!$acc data present(vkb(:,:))
!$acc host_data use_device(vkb)
     CALL calbec_gpu( npw, vkb, wfcatom_d, becp_d )
!$acc end host_data
!$acc end data
     !
     CALL s_psi_gpu( npwx, npw, natomwfc, wfcatom_d, swfcatom_d )
     !
     IF (orthogonalize_wfc) &
        CALL ortho_swfc_gpu( npw, normalize_only, natomwfc, wfcatom_d, swfcatom_d, .FALSE. )
     !
     ! copy S * atomic wavefunctions with Hubbard U term only in wfcU
     ! (this is used during the self-consistent solution of Kohn-Sham equations)
     ! save to unit iunhub
     !
     swfcatom = swfcatom_d
     CALL copy_U_wfc (swfcatom, noncolin)
     IF ( nks > 1 ) &
          CALL save_buffer (wfcU, nwordwfcU, iunhub, ik)
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
  ! if "orthogonalize_wfc" is .true., saves them into buffer "iunsat".
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
  USE becmod_gpum,      ONLY : becp_d
  USE becmod_subs_gpum, ONLY : using_becp_auto, using_becp_d_auto, &
                               calbec_gpu
  USE uspp_init,        ONLY : init_us_2
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
     !
     CALL init_us_2( npw, igk_k(1,ik), xk(1,ik), vkb, .true. )
     !
     CALL using_becp_auto(2)
     CALL using_becp_d_auto(2)
!$acc data present(vkb(:,:))
!$acc host_data use_device(vkb)
     CALL calbec_gpu( npw, vkb, wfcatom_d, becp_d )
!$acc end host_data
!$acc end data
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
  ! On input : 
  ! wfc (npwx*npol,m) =  \phi = a set of "m" (atomic) wavefcts
  ! swfc(npwx*npol,m) = S\phi 
  ! normalize_only    = only normalize, do not orthonormalize
  !
  ! On output this routine will compute the overlap matrix O: 
  ! O_ij = <wfc_i|S|wfc_j> = <wfc_i|swfc_j>
  ! If lflag=.FALSE. : wfc are unchanged,   swfc = O^{-1/2} S\phi.
  ! If lflag=.TRUE.  : wfc = O^{-1/2} \phi, swfc are unchanged.
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
  USE force_mod,        ONLY : eigenval, eigenvect, overlap_inv
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
  ! s_d = CMPLX(0.d0,0.d0, kind=dp)  ! fused below
  !$cuf kernel do (2)
  DO i = 1, m
     DO j = 1, m
        s_d(j,i) = CMPLX(0.d0,0.d0, kind=dp)
        IF (i == j) s_d(j,i) = CMPLX(1.d0,0.d0, kind=dp)
     ENDDO
  ENDDO
  ! THIS SHOULD BE A SIMPLE CDIAGH (NOT GENERALIZED!) DRIVER NEEDED IN LAXLIB
  CALL laxlib_cdiaghg_gpu( m, m, overlap_d, s_d, m, e_d, work_d, me_bgrp, &
                           root_bgrp, intra_bgrp_comm )
  !
  !$cuf kernel do (2) <<<*,*>>>
  DO i = 1, m
     DO j = 1, m
        IF ( j < i ) CYCLE
        temp = (0.d0, 0.d0)
        DO k = 1, m
           temp = temp + work_d(j,k) * 1.d0/SQRT(e_d(k)) * CONJG(work_d(i,k))
        ENDDO
        overlap_d(i,j) = temp
        IF (j /= i) overlap_d(j,i) = CONJG(temp)
     ENDDO
  ENDDO
  !
  IF (lflag) THEN
     !
     ! Save quantities which are needed for 
     ! calculations of Hubbard forces and stress
     eigenval(:) = e_d(:)
     eigenvect(:,:) = work_d(:,:)
     overlap_inv(:,:) = overlap_d(:,:)
     !
     ! Transform atomic orbitals WITHOUT the ultrasoft S operator 
     ! O^(-1/2) \psi (note the transposition):
     ! \phi_I = \sum_J O^{-1/2}_JI \phi_J
     !
     DO i = 1, npw
        work_d(:,1) = (0.d0,0.d0)
        IF (noncolin) THEN
           DO ipol=1,npol
              j = i + (ipol-1)*npwx
              CALL cublasZgemv ('n',m,m,(1.d0,0.d0),overlap_d, &
                   m, wfc_d(j,1), npwx*npol, (0.d0,0.d0),work_d,1)
              CALL zcopy (m,work_d,1,wfc_d(j,1),npwx*npol)
           END DO
        ELSE
           CALL cublasZgemv ('n', m, m, (1.d0, 0.d0) , overlap_d, &
                m, wfc_d (i, 1) , npwx, (0.d0, 0.d0) , work_d, 1)
           CALL zcopy (m, work_d, 1, wfc_d (i, 1), npwx)
        END IF
     ENDDO
     !
  ELSE
     !
     ! Transform atomic orbitals WITH the ultrasoft S operator 
     ! O^(-1/2) \Spsi (note the transposition):
     ! \Sphi_I = \sum_J O^{-1/2}_JI \Sphi_J
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
  ENDIF
  !
  DEALLOCATE( overlap_d, s_d )
  DEALLOCATE( work_d, e_d )
  !
  RETURN
  !      
END SUBROUTINE ortho_swfc_gpu

!
!-----------------------------------------------------------------------
SUBROUTINE calculate_doverlap_inv_gpu (m, e, work, doverlap, doverlap_inv)
  !---------------------------------------------------------------------
  !! This routine computes the derivative of O^{-1/2}, i.e.
  !! [d((O^{-1/2}))]_IJ, where O_IJ is the overlap matrix. 
  !! Note, on the input this routine requires dO (not transposed).
  !! The solution is written in a closed form by solving the Lyapunov
  !! equation (a particular case of the Sylvester equation).
  !! See Eq. (32) in PRB 102, 235159 (2020).
  !! Written by I. Timrov (June 2020)
  !
#if defined(__CUDA)  
  USE cublas
#endif
  USE kinds,       ONLY : DP
  !
  IMPLICIT NONE
  !
  ! I/O variables
  !
  INTEGER, INTENT(IN)      :: m
  !! The total number of atomic functions 
  REAL(DP), INTENT(IN)     :: e(m)
  !! The eigenvalues of the overlap matrix
  COMPLEX(DP), INTENT(IN)  :: work(m,m)
  !! The eigenvectors of the overlap matrix
  COMPLEX(DP), INTENT(IN)  :: doverlap(m,m)
  !! The derivative of the overlap matrix O_IJ (not transposed)  
  COMPLEX(DP), INTENT(OUT) :: doverlap_inv(m,m)
  !! The derivative of transposed O^{-1/2}
  !
  ! Local variables
  INTEGER :: m1, m2, m3, m4
  COMPLEX(DP), ALLOCATABLE :: aux(:,:)
  !! eigenvectors of the overlap matrix
  !! auxiliary array
  !
#if defined(__CUDA)  
  attributes(DEVICE) :: e, work, doverlap, doverlap_inv, aux
#endif
  ALLOCATE (aux(m,m))
  !
  ! Compute (work^H) * doverlap * work 
  ! and put the result back in doverlap
  !
  ! Compute aux = doverlap * work
  CALL ZGEMM('N','N', m, m, m, (1.d0,0.d0), doverlap, &
              m, work, m, (0.d0,0.d0), aux, m)
  ! Compute (work^H) * aux
  CALL ZGEMM('C','N', m, m, m, (1.d0,0.d0), work, &
              m, aux, m, (0.d0,0.d0), doverlap, m)
  !
  !$cuf kernel do(2)
  DO m1 = 1, m
     DO m2 = 1, m
        aux(m1,m2) = doverlap(m1,m2) / &
                    (e(m1)*DSQRT(e(m2))+e(m2)*DSQRT(e(m1)))
     ENDDO
  ENDDO
  !
  ! Compute work * aux * (work^H)
  !
  ! Compute doverlap = aux * (work^H)
  CALL ZGEMM('N','C', m, m, m, (1.d0,0.d0), aux, &
              m, work, m, (0.d0,0.d0), doverlap, m)
  ! Compute doverlap_inv = work * doverlap
  CALL ZGEMM('N','N', m, m, m, (-1.d0,0.d0), work, &
              m, doverlap, m, (0.d0,0.d0), doverlap_inv, m)
  !
  DEALLOCATE (aux)
  !
  RETURN
  !
END SUBROUTINE calculate_doverlap_inv_gpu
