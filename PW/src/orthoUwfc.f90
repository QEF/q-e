!
! Copyright (C) 2001-2024 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
SUBROUTINE orthoUwfc(save_wfcatom)
  !-----------------------------------------------------------------------
  !
  ! This routine saves to buffer "iunhub" atomic wavefunctions having an
  ! associated Hubbard U term * S, for DFT+U(+V) calculations. Same for 
  ! "iunhub_noS" but without S (this is then used for plotting Hubbard projector 
  ! functions or other post-processing operations). Atomic wavefunctions
  ! are orthogonalized if desired, depending upon the value of "Hubbard_projectors"
  ! If save_wfcatom == .TRUE., also write atomic wavefunctions before
  ! applying S to buffer.
  !
  USE kinds,      ONLY : DP
  USE buffers,    ONLY : get_buffer, save_buffer
  USE io_global,  ONLY : stdout
  USE io_files,   ONLY : iunhub, iunhub_noS, nwordwfcU
  USE ions_base,  ONLY : nat
  USE basis,      ONLY : natomwfc
  USE klist,      ONLY : nks, xk, ngk, igk_k
  USE ldaU,       ONLY : Hubbard_projectors, wfcU, nwfcU, copy_U_wfc
  USE wvfct,      ONLY : npwx
  USE uspp,       ONLY : nkb, vkb
  USE becmod,     ONLY : allocate_bec_type_acc, deallocate_bec_type_acc, &
                         bec_type, becp, calbec
  USE control_flags,    ONLY : gamma_only, use_gpu, offload_type
  USE noncollin_module, ONLY : noncolin, npol
  USE mp_bands,         ONLY : use_bgrp_in_hpsi
  USE uspp_init,        ONLY : init_us_2
  IMPLICIT NONE
  !
  LOGICAL, INTENT(IN) :: save_wfcatom
  !! If .TRUE., write atomic wavefunction before applying S to buffer
  !
  INTEGER :: ik, ibnd, info, i, j, k, na, nb, nt, isym, n, ntemp, m, &
       l, lm, ltot, ntot, ipol, npw
  ! ik: the k point under consideration
  ! ibnd: counter on bands
  LOGICAL :: orthogonalize_wfc, normalize_only, save_flag
  COMPLEX(DP) , ALLOCATABLE :: wfcatom (:,:), swfcatom(:,:), wfcUaux (:,:)
  !
  IF ( Hubbard_projectors == "pseudo" ) THEN
     WRITE( stdout,'(/5x,a,/)') 'Beta functions used for Hubbard projectors'
     RETURN
  ELSE IF (Hubbard_projectors=="wf") THEN
     ! Read Wannier functions from file (produced by Wannier90 or pmw.x).
     WRITE( stdout,'(/5x,a,/)') 'Wannier functions used for Hubbard projectors'
     CALL read_wf_projectors (save_wfcatom)
     RETURN
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
  ALLOCATE ( wfcatom(npwx*npol, natomwfc), swfcatom(npwx*npol, natomwfc) )
  !$acc enter data create(wfcatom, swfcatom)
  !
  save_flag = use_bgrp_in_hpsi ; use_bgrp_in_hpsi=.false.
  !
  ! Allocate the array becp = <beta|wfcatom>
  CALL allocate_bec_type_acc (nkb,natomwfc, becp)
  !
  DO ik = 1, nks
     !
     IF (noncolin) THEN
       CALL atomic_wfc_nc_updown (ik, wfcatom)
     ELSE
       CALL atomic_wfc (ik, wfcatom)
     ENDIF
     npw = ngk (ik)
     CALL init_us_2 (npw, igk_k(1,ik), xk (1, ik), vkb, use_gpu)
     CALL calbec (offload_type, npw, vkb, wfcatom, becp)
     CALL s_psi_acc (npwx, npw, natomwfc, wfcatom, swfcatom)
     !
     IF (orthogonalize_wfc) CALL ortho_swfc ( npw, normalize_only, natomwfc, wfcatom, swfcatom, .FALSE. )
     !
     ! copy S * atomic wavefunctions with Hubbard U term only in wfcU
     ! (this is used during the self-consistent solution of Kohn-Sham equations)
     ! save to unit iunhub
     !
     CALL copy_U_wfc (swfcatom, noncolin)
     !$acc update host(wfcU)
     IF ( nks > 1 ) CALL save_buffer (wfcU, nwordwfcU, iunhub, ik)
     !
     ! If save_wfcatom=.TRUE. copy the orthonormalized wfcatom to wfcU and save
     ! to unit iunhub_noS
     !
     IF (save_wfcatom) THEN
        ! Calculate swfcatom = S * \phi
        CALL s_psi_acc (npwx, npw, natomwfc, wfcatom, swfcatom)
        IF (orthogonalize_wfc) CALL ortho_swfc ( npw, normalize_only, natomwfc, wfcatom, swfcatom, .TRUE. )
        ! If nks=1, wfcU is kept in memory, while we want to use wfcU as a 
        ! workspace. Hence we store it temporarily in the auxiliary array wfcUaux
        IF (nks==1) THEN
           ALLOCATE (wfcUaux(npwx*npol,nwfcU))
           wfcUaux = wfcU
        ENDIF
        CALL copy_U_wfc (wfcatom, noncolin)
        !$acc update host(wfcU)
        ! Write wfcU = O^{-1/2} \phi (no ultrasoft S)
        CALL save_buffer (wfcU, nwordwfcU, iunhub_noS, ik)
        IF (nks==1) THEN
           wfcU = wfcUaux
           !$acc update device(wfcU)
           DEALLOCATE (wfcUaux)
        ENDIF
     ENDIF
     !
  ENDDO
  !$acc exit data delete(wfcatom, swfcatom)
  DEALLOCATE (wfcatom, swfcatom)
  CALL deallocate_bec_type_acc ( becp )
  !
  use_bgrp_in_hpsi = save_flag
  !
  RETURN
  !
END SUBROUTINE orthoUwfc
!
!-----------------------------------------------------------------------
SUBROUTINE orthoUwfc_k (ik, lflag)
  !-----------------------------------------------------------------------
  !
  ! For a given k point "ik", this routine computes (ortho-)atomic wavefunctions 
  ! having an associated Hubbard U term * S, for DFT+U(+V) calculations. 
  ! Also without S (this is then used to computed Hubbard forces and stresses). 
  ! wfcatom and swfcatom must be allocated on input.
  ! Beta functions vkb must be already computed before.
  !
  ! lflag=.TRUE.  : wfcU = O^{-1/2}  \phi (w/o ultrasoft S)
  ! lflag=.FALSE. : wfcU = O^{-1/2} S\phi (w/  ultrasoft S)
  !
  USE kinds,            ONLY : DP
  USE io_global,        ONLY : stdout
  USE ions_base,        ONLY : nat
  USE basis,            ONLY : natomwfc, wfcatom, swfcatom
  USE klist,            ONLY : nks, xk, ngk, igk_k
  USE ldaU,             ONLY : Hubbard_projectors, wfcU, nwfcU, copy_U_wfc
  USE wvfct,            ONLY : npwx
  USE uspp,             ONLY : nkb, vkb
  USE becmod,           ONLY : allocate_bec_type_acc, deallocate_bec_type_acc, &
                               bec_type, becp, calbec
  USE control_flags,    ONLY : gamma_only, offload_type
  USE noncollin_module, ONLY : noncolin, npol
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: ik ! the k point under consideration
  LOGICAL, INTENT(IN) :: lflag
  !
  INTEGER :: ibnd, info, i, j, k, na, nb, nt, isym, n, ntemp, m, &
             l, lm, ltot, ntot, ipol, npw
  LOGICAL :: orthogonalize_wfc, normalize_only, save_flag
  COMPLEX(DP), ALLOCATABLE :: aux(:,:)
  !$acc declare device_resident(aux)

  IF ( Hubbard_projectors == "pseudo" ) THEN
     CALL errore ("orthoUwfc_k","Hubbard_projectors=pseudo is not supported",1)
  ELSE IF (Hubbard_projectors=="wf") THEN
     CALL errore ("orthoUwfc_k","Hubbard_projectors=wf is not supported",1)
  ELSE IF (Hubbard_projectors=="atomic") THEN
     orthogonalize_wfc = .FALSE.
     normalize_only = .FALSE.
  ELSE IF (Hubbard_projectors=="ortho-atomic") THEN
     orthogonalize_wfc = .TRUE.
     normalize_only = .FALSE.    
     IF (gamma_only) CALL errore('orthoUwfc_k', &
          'Gamma-only calculation for this case not implemented', 1 )
  ELSE IF (Hubbard_projectors=="norm-atomic") THEN
     CALL errore ("orthoUwfc_k","Hubbard_projectors=norm-atomic is not supported",1)
  ELSE
     WRITE(stdout,'(/5x,"Hubbard_projectors = ",a)') Hubbard_projectors
     CALL errore ("orthoUwfc_k"," this Hubbard_projectors type is not valid",1)
  END IF
  !
  !$acc data present(wfcatom, swfcatom)
  ! Compute atomic wfc at this k (phi)
  IF (noncolin) THEN
     CALL atomic_wfc_nc_updown (ik, wfcatom)
  ELSE
     CALL atomic_wfc (ik, wfcatom)
  ENDIF
  !
  IF (Hubbard_projectors=="ortho-atomic") THEN
     ALLOCATE(aux(npwx*npol,natomwfc))
     ! Copy atomic wfcs (phi)
     !$acc kernels
     aux(:,:) = wfcatom(:,:)
     !$acc end kernels
  ENDIF
  !
  ! Number of plane waves at this k point
  npw = ngk(ik)
  !
  IF (orthogonalize_wfc .OR. .NOT.lflag) THEN
     ! Allocate the array becp = <beta|wfcatom>
     CALL allocate_bec_type_acc (nkb,natomwfc, becp)
     CALL calbec (offload_type, npw, vkb, wfcatom, becp)
     ! Calculate swfcatom = S * phi
     CALL s_psi_acc (npwx, npw, natomwfc, wfcatom, swfcatom)
     CALL deallocate_bec_type_acc (becp)
  ENDIF
  !
  ! Compute the overlap matrix
  ! lflag=.FALSE. : On the output wfcatom are unchanged, swfcatom = O^{-1/2} S\phi.
  ! lflag=.TRUE.  : On the output wfcatom = O^{-1/2} \phi (no ultrasoft S), swfcatom are unchanged.
  IF (orthogonalize_wfc) THEN
     CALL ortho_swfc ( npw, normalize_only, natomwfc, wfcatom, swfcatom, lflag )
  END IF
  !
  IF (lflag) THEN
     ! Copy (ortho-)atomic wavefunctions with Hubbard U term only
     ! in wfcU (no ultrasoft S): wfcatom = O^{-1/2} \phi.
     CALL copy_U_wfc (wfcatom, noncolin)
  ELSE
     ! Copy (ortho-)atomic wavefunctions with Hubbard U term only
     ! in wfcU (with ultrasoft S): swfcatom = O^{-1/2} S\phi.
     CALL copy_U_wfc (swfcatom, noncolin)
  ENDIF
  !$acc update host(wfcU)
  !
  IF (Hubbard_projectors=="ortho-atomic") THEN
     ! Copy atomic wfcs
     !$acc kernels
     wfcatom(:,:) = aux(:,:)
     !$acc end kernels
     DEALLOCATE(aux)
  ENDIF
  !$acc end data
  !
  RETURN
  !   
END SUBROUTINE orthoUwfc_k
!
!-----------------------------------------------------------------------
SUBROUTINE ortho_swfc ( npw, normalize_only, m, wfc, swfc, lflag )
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
  USE kinds,            ONLY : DP
  USE wvfct,            ONLY : npwx
  USE mp_bands,         ONLY : intra_bgrp_comm, me_bgrp, root_bgrp
  USE mp,               ONLY : mp_sum
  USE noncollin_module, ONLY : noncolin, npol
  USE force_mod,        ONLY : eigenval, eigenvect, overlap_inv
  USE control_flags,    ONLY : use_gpu
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: m, npw
  LOGICAL, INTENT(IN) :: normalize_only
  COMPLEX(dp), INTENT(INOUT) :: wfc (npwx*npol,m)
  COMPLEX(dp), INTENT(INOUT) :: swfc(npwx*npol,m)
  LOGICAL, INTENT(IN) :: lflag
  !
  ! ... local variables
  !
  COMPLEX(DP) :: temp 
  COMPLEX(DP) , ALLOCATABLE ::  work (:,:), overlap (:,:)
  REAL(DP) , ALLOCATABLE :: e (:)
  COMPLEX(DP) , ALLOCATABLE ::  s(:,:)
  !$acc declare device_resident(work, overlap, e, s)
  INTEGER :: i, j, k, ipol
  !
  ALLOCATE (overlap(m,m), work(m,m), e(m), s(m,m))
  ! 
  !$acc kernels
  overlap(:,:) = (0.d0,0.d0)
  work(:,:) = (0.d0,0.d0)
  !$acc end kernels
  !
  ! calculate overlap matrix
  !
  IF (noncolin) THEN
     !$acc host_data use_device(wfc, swfc, overlap)
     CALL MYZGEMM ('c', 'n', m, m, npwx*npol, (1.d0, 0.d0), wfc, &
          npwx*npol, swfc, npwx*npol, (0.d0,0.d0), overlap, m)
     !$acc end host_data
  ELSE
     !$acc host_data use_device(wfc, swfc, overlap)
     CALL MYZGEMM ('c', 'n', m, m, npw, (1.d0, 0.d0), wfc, &
          npwx, swfc, npwx, (0.d0, 0.d0), overlap, m)
     !$acc end host_data
  END IF
  !
  !$acc host_data use_device(overlap)
  CALL mp_sum(  overlap, intra_bgrp_comm )
  !$acc end host_data
  !
  IF ( normalize_only ) THEN
     !$acc parallel
     !$acc loop gang
     DO i = 1, m
        !$acc loop vector
        DO j = i+1, m
           overlap(i,j) = CMPLX(0.d0,0.d0, kind=dp)
           overlap(j,i) = CMPLX(0.d0,0.d0, kind=dp)
        ENDDO
     ENDDO
     !$acc end parallel
  END IF
  !
  ! find O^(-1/2) (actually, its transpose)
  !
!civn: ZHEEV not available in cuBLAS/cuSOLVER?
  IF(use_gpu) THEN
    !
    ! s_d = CMPLX(0.d0,0.d0, kind=dp)  ! fused below
    !$acc kernels
    s(:,:) = CMPLX(0.d0,0.d0, kind=dp) 
    DO i = 1, m
       s(i,i) = CMPLX(1.d0,0.d0, kind=dp)
    ENDDO
    !$acc end kernels
    ! THIS SHOULD BE A SIMPLE CDIAGH (NOT GENERALIZED!) DRIVER NEEDED IN LAXLIB
    !$acc host_data use_device(overlap, s, e, work)
    CALL laxlib_cdiaghg_gpu( m, m, overlap, s, m, e, work, me_bgrp, &
                             root_bgrp, intra_bgrp_comm )
    !$acc end host_data
    !
  ELSE
    CALL cdiagh (m, overlap, m, e, work)
  END IF 
  !
  !$acc parallel loop collapse(2) private(temp)
  DO i = 1, m
     DO j = 1, m
        IF ( j < i ) CYCLE
        temp = (0.d0, 0.d0)
        !$acc loop seq
        DO k = 1, m
           temp = temp + work(j,k) * 1.d0/SQRT(e(k)) * CONJG(work(i,k))
        ENDDO
        overlap(i,j) = temp
        IF (j /= i) overlap(j,i) = CONJG(temp)
     ENDDO
  ENDDO
  !
  IF (lflag) THEN
     !
     ! Save quantities which are needed for 
     ! calculations of Hubbard forces and stress
     IF (allocated(eigenval)) THEN
        !$acc kernels
        eigenval(:) = e(:)
        !$acc end kernels
     ENDIF
     IF (allocated(eigenvect)) THEN
        !$acc kernels
        eigenvect(:,:) = work(:,:)
        !$acc end kernels
     ENDIF
     IF (allocated (overlap_inv)) THEN
        !$acc kernels 
        overlap_inv(:,:) = overlap(:,:)
        !$acc end kernels
     ENDIF
     !
  END IF 
  !
  DEALLOCATE( work )
  !
  ALLOCATE( work(m, npwx*npol ) )
  !
  !$acc kernels
  work(:,:) = (0.d0,0.d0)
  !$acc end kernels
  !
  IF (lflag) THEN
     !
     ! Transform atomic orbitals WITHOUT the ultrasoft S operator 
     ! O^(-1/2) \psi (note the transposition):
     ! \phi_I = \sum_J O^{-1/2}_JI \phi_J
     !
     IF(noncolin) THEN 
       !$acc host_data use_device(overlap, wfc, work)
       CALL MYZGEMM('n', 't', m, npwx*npol, m, (1.d0,0.d0), overlap, m, wfc, npwx*npol, (0.d0,0.d0), work, m )
       !$acc end host_data
       !$acc parallel loop collapse(2) 
       DO i = 1, npwx*npol
         DO j = 1, m
           wfc(i,j) = work(j,i)
         END DO 
       END DO
     ELSE
       !$acc host_data use_device(overlap, wfc, work)
       CALL MYZGEMM('n', 't', m, npw, m, (1.d0,0.d0), overlap, m, wfc, npwx*npol, (0.d0,0.d0), work, m )
       !$acc end host_data
       !$acc parallel loop collapse(2)
       DO i = 1, npw
         DO j = 1, m
           wfc(i,j) = work(j,i)
         END DO 
       END DO
     END IF
    
     !
     !
  ELSE
     !
     ! Transform atomic orbitals WITH the ultrasoft S operator 
     ! O^(-1/2) \Spsi (note the transposition):
     ! \Sphi_I = \sum_J O^{-1/2}_JI \Sphi_J
     ! FIXME: can be done in a faster way by using wfc as work space 
     !
     IF(noncolin) THEN 
       !$acc host_data use_device(overlap, swfc, work)
       CALL MYZGEMM('n', 't', m, npwx*npol, m, (1.d0,0.d0), overlap, m, swfc, npwx*npol, (0.d0,0.d0), work, m )
       !$acc end host_data 
       !$acc parallel loop collapse(2)
       DO i = 1, npwx*npol
         DO j = 1, m
           swfc(i,j) = work(j,i)
         END DO 
       END DO
     ELSE
       !$acc host_data use_device(overlap, swfc, work)
       CALL MYZGEMM('n', 't', m, npw, m, (1.d0,0.d0), overlap, m, swfc, npwx*npol, (0.d0,0.d0), work, m )
       !$acc end host_data
       !$acc parallel loop collapse(2)
       DO i = 1, npw
         DO j = 1, m
           swfc(i,j) = work(j,i)
         END DO 
       END DO
     END IF
     !
  ENDIF
  !
  DEALLOCATE (overlap, work, e, s)
  !
  RETURN
  !      
END SUBROUTINE ortho_swfc
!
!-----------------------------------------------------------------------
SUBROUTINE calculate_doverlap_inv (m, e, work, doverlap, doverlap_inv)
  !---------------------------------------------------------------------
  !! This routine computes the derivative of O^{-1/2}, i.e.
  !! [d((O^{-1/2}))]_IJ, where O_IJ is the overlap matrix. 
  !! Note, on the input this routine requires dO (not transposed).
  !! The solution is written in a closed form by solving the Lyapunov
  !! equation (a particular case of the Sylvester equation).
  !! See Eq. (32) in PRB 105, 199901(E) (2022).
  !! See Eq. (32) in PRB 102, 235159 (2020).
  !! Written by I. Timrov (June 2020)
  !
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
  !$acc declare device_resident(aux)
  !! eigenvectors of the overlap matrix
  !! auxiliary array
  !
  ALLOCATE (aux(m,m))
  !
  ! Compute (work^H) * doverlap * work 
  ! and put the result back in doverlap
  !
  ! Compute aux = doverlap * work
  !$acc host_data use_device(doverlap, work, aux)
  CALL MYZGEMM('N','N', m, m, m, (1.d0,0.d0), doverlap, &
              m, work, m, (0.d0,0.d0), aux, m)
  !$acc end host_data
  ! Compute (work^H) * aux
  !$acc host_data use_device(work, aux, doverlap)
  CALL MYZGEMM('C','N', m, m, m, (1.d0,0.d0), work, &
              m, aux, m, (0.d0,0.d0), doverlap, m)
  !$acc end host_data
  !
  !$acc parallel loop collapse(2)
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
  !$acc host_data use_device(aux, work, doverlap)
  CALL MYZGEMM('N','C', m, m, m, (1.d0,0.d0), aux, &
              m, work, m, (0.d0,0.d0), doverlap, m)
  !$acc end host_data
  ! Compute doverlap_inv = work * doverlap
  !$acc host_data use_device(work, doverlap, doverlap_inv)
  CALL MYZGEMM('N','N', m, m, m, (-1.d0,0.d0), work, &
              m, doverlap, m, (0.d0,0.d0), doverlap_inv, m)
  !$acc end host_data
  !
  DEALLOCATE (aux)
  !
  RETURN
  !
END SUBROUTINE calculate_doverlap_inv

SUBROUTINE read_wf_projectors (save_wfcatom)
  !--------------------------------------------------------------------------
  !
  !! This routine reads Hubbard projectors (Wannier functions)
  !! from file produced by Wannier90 or pmw.x.
  !! The routine reads collected wfcU (no S), and it writes 
  !! S*wfcU in a distributed format
  !!
  !! Written by I. Timrov (October 2024)
  !
  USE kinds,            ONLY : DP
  USE io_global,        ONLY : stdout
  USE klist,            ONLY : nkstot, nks, xk, ngk, igk_k
  USE io_files,         ONLY : nwordwfcU, iunhub, iunhub_noS
  USE wvfct,            ONLY : npwx
  USE buffers,          ONLY : get_buffer, save_buffer
  USE ldaU,             ONLY : wfcU, nwfcU
  USE uspp,             ONLY : nkb, vkb, okvan
  USE uspp_init,        ONLY : init_us_2
  USE becmod,           ONLY : becp, calbec, allocate_bec_type, deallocate_bec_type
  USE noncollin_module, ONLY : npol
  !
  IMPLICIT NONE
  INTEGER :: ik,  & ! dummy index to number k points
             npw, & ! number of plane waves
             ierr
  LOGICAL, INTENT(IN) :: save_wfcatom
  COMPLEX(DP), ALLOCATABLE :: swfcU(:,:)
  !
  WRITE(stdout, '(5x,A)') 'Reading WFs to construct the Hubbard projectors...'
  !
  IF (okvan) THEN
     ! Allocate the array containing <beta|wfcU>
     CALL allocate_bec_type (nkb, nwfcU, becp)
     ALLOCATE (swfcU(npwx*npol,nwfcU))
  ENDIF
  !
  DO ik = 1, nks
     !
     ! Read wfcU (no S) from iunhub
     CALL get_buffer (wfcU, nwordwfcU, iunhub, ik) 
     !
     ! Write wfcU (no S) to iunhub_noS for further post-processing
     IF (save_wfcatom) CALL save_buffer (wfcU, nwordwfcU, iunhub_noS, ik)
     !
     IF (okvan) THEN
        npw = ngk(ik)
        ! USPP: Compute swfcU = S*wfcU
        CALL init_us_2 (npw, igk_k(1,ik), xk(1,ik), vkb)
        CALL calbec (npw, vkb, wfcU, becp)
        CALL s_psi (npwx, npw, nwfcU, wfcU, swfcU)
        ! Write swfcU
        CALL save_buffer (swfcU, nwordwfcU, iunhub, ik)
     ENDIF
     !
  ENDDO
  !
  IF (okvan) THEN
     CALL deallocate_bec_type (becp)
     DEALLOCATE (swfcU)
  ENDIF
  !
  RETURN
  !
END SUBROUTINE read_wf_projectors
