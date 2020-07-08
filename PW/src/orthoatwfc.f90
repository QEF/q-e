!
! Copyright (C) 2001-2020 Quantum ESPRESSO group
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
  ! associated Hubbard U term * S, for DFT+U(+V) calculations. Same for 
  ! Atomic wavefunctions are orthogonalized if desired, depending upon 
  ! the value of "U_projection". "swfcatom" must NOT be allocated on input.
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
  USE mp_bands,         ONLY : use_bgrp_in_hpsi
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
     IF (gamma_only) CALL errore('orthoUwfc', &
          'Gamma-only calculation for this case not implemented', 1 )
  ELSE IF (U_projection=="norm-atomic") THEN
     orthogonalize_wfc = .TRUE.
     normalize_only = .TRUE.
     WRITE( stdout,*) 'Atomic wfc used for LDA+U Projector are normalized but NOT orthogonalized'
     IF (gamma_only) CALL errore('orthoUwfc', &
          'Gamma-only calculation for this case not implemented', 1 )
  ELSE
     WRITE( stdout,*) "U_projection_type =", U_projection
     CALL errore ("orthoUwfc"," this U_projection_type is not valid",1)
  END IF

  ALLOCATE ( wfcatom(npwx*npol, natomwfc), swfcatom(npwx*npol, natomwfc) )
  
  save_flag = use_bgrp_in_hpsi ; use_bgrp_in_hpsi=.false.

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
        CALL ortho_swfc ( npw, normalize_only, natomwfc, wfcatom, swfcatom, .FALSE. )
     !
     ! copy S * atomic wavefunctions with Hubbard U term only in wfcU
     ! (this is used during the self-consistent solution of Kohn-Sham equations)
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
  use_bgrp_in_hpsi = save_flag
  !
  RETURN
  !     
END SUBROUTINE orthoUwfc
!
!-----------------------------------------------------------------------
SUBROUTINE orthoUwfc2 (ik)
  !-----------------------------------------------------------------------
  !
  ! For a given k point "ik", this routine computes (ortho-)atomic wavefunctions 
  ! having an associated Hubbard U term * S, for DFT+U(+V) calculations. 
  ! Also without S (this is then used to computed Hubbard forces 
  ! and stresses). 
  ! wfcatom and swfcatom must be allocated on input.
  ! Beta functions vkb must be already computed before.
  !
  USE kinds,            ONLY : DP
  USE io_global,        ONLY : stdout
  USE io_files,         ONLY : iunhub, nwordwfcU
  USE ions_base,        ONLY : nat
  USE basis,            ONLY : natomwfc, wfcatom, swfcatom
  USE klist,            ONLY : nks, xk, ngk, igk_k
  USE ldaU,             ONLY : U_projection, wfcU, nwfcU, copy_U_wfc
  USE wvfct,            ONLY : npwx
  USE uspp,             ONLY : nkb, vkb
  USE becmod,           ONLY : allocate_bec_type, deallocate_bec_type, &
                               bec_type, becp, calbec
  USE control_flags,    ONLY : gamma_only
  USE noncollin_module, ONLY : noncolin 
  ! 
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: ik ! the k point under consideration
  !
  INTEGER :: ibnd, info, i, j, k, na, nb, nt, isym, n, ntemp, m, &
             l, lm, ltot, ntot, ipol, npw
  LOGICAL :: orthogonalize_wfc, normalize_only, save_flag
  COMPLEX(DP), ALLOCATABLE :: aux(:,:)

  IF ( U_projection == "pseudo" ) THEN
     CALL errore ("orthoUwfc2","U_projection=pseudo is not supported",1)
  ELSE IF (U_projection=="file") THEN
     CALL errore ("orthoUwfc2","U_projection=file is not supported",1)
  ELSE IF (U_projection=="atomic") THEN
     orthogonalize_wfc = .FALSE.
     normalize_only = .FALSE.
  ELSE IF (U_projection=="ortho-atomic") THEN
     orthogonalize_wfc = .TRUE.
     normalize_only = .FALSE.    
     IF (gamma_only) CALL errore('orthoUwfc2', &
          'Gamma-only calculation for this case not implemented', 1 )
  ELSE IF (U_projection=="norm-atomic") THEN
     CALL errore ("orthoUwfc2","U_projection=norm-atomic is not supported",1)
  ELSE
     WRITE( stdout,*) "U_projection_type =", U_projection
     CALL errore ("orthoUwfc2"," this U_projection_type is not valid",1)
  END IF
  !
  ! Compute atomic wfc at this k (phi)
  IF (noncolin) THEN
     CALL atomic_wfc_nc_updown (ik, wfcatom)
  ELSE
     CALL atomic_wfc (ik, wfcatom)
  ENDIF
  !
  IF (U_projection=="ortho-atomic") THEN
     ALLOCATE(aux(npwx,natomwfc))
     ! Copy atomic wfcs (phi)
     aux(:,:) = wfcatom(:,:)
  ENDIF
  !
  IF (orthogonalize_wfc) THEN
     !
     ! Number of plane waves at this k point
     npw = ngk(ik)
     !
     ! Allocate the array becp = <beta|wfcatom>
     CALL allocate_bec_type (nkb,natomwfc, becp)
     CALL calbec (npw, vkb, wfcatom, becp)
     ! Calculate swfcatom = S * phi
     CALL s_psi (npwx, npw, natomwfc, wfcatom, swfcatom)
     CALL deallocate_bec_type (becp)
     !  
     ! Compute the overlap matrix
     ! On the output: wfcatom = O^{-1/2} \phi (no ultrasoft S)
     CALL ortho_swfc ( npw, normalize_only, natomwfc, wfcatom, swfcatom, .TRUE. )
     !
  ENDIF
  !
  ! Copy (ortho-)atomic wavefunctions with Hubbard U term only 
  ! in wfcU (no ultrasoft S).
  CALL copy_U_wfc (wfcatom, noncolin)
  !
  IF (U_projection=="ortho-atomic") THEN
     ! Copy atomic wfcs
     wfcatom(:,:) = aux(:,:)
     DEALLOCATE(aux)
  ENDIF
  !
  RETURN
  !   
END SUBROUTINE orthoUwfc2
!
!-----------------------------------------------------------------------
SUBROUTINE orthoatwfc (orthogonalize_wfc)
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
        CALL ortho_swfc ( npw, normalize_only, natomwfc, wfcatom, swfcatom, .FALSE. )
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
  USE mp_bands,         ONLY : intra_bgrp_comm
  USE mp,               ONLY : mp_sum
  USE noncollin_module, ONLY : noncolin, npol
  USE force_mod,        ONLY : eigenval, eigenvect, overlap_inv
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: m, npw
  LOGICAL, INTENT(IN) :: normalize_only
  COMPLEX(dp), INTENT(INOUT) :: wfc (npwx*npol,m)
  COMPLEX(dp), INTENT(INOUT) :: swfc(npwx*npol,m)
  LOGICAL, INTENT(IN) :: lflag

  COMPLEX(DP) :: temp 
  COMPLEX(DP) , ALLOCATABLE ::  work (:,:), overlap (:,:)
  REAL(DP) , ALLOCATABLE :: e (:)
  INTEGER :: i, j, k, ipol

  ALLOCATE (overlap( m , m))    
  ALLOCATE (work   ( m , m))    
  ALLOCATE (e      ( m))    
  ! 
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
  ! find O^(-1/2) (actually, its transpose)
  !
  CALL cdiagh (m, overlap, m, e, work)
  !
  DO i = 1, m
     DO j = i, m
        temp = (0.d0, 0.d0)
        DO k = 1, m
           temp = temp + work (j, k) * (1.d0/SQRT(e(k))) * CONJG (work (i, k) )
        ENDDO
        overlap (i, j) = temp
        IF (j.NE.i) overlap (j, i) = CONJG (temp)
     ENDDO
  ENDDO
  !
  IF (lflag) THEN
     !
     ! Save quantities which are needed for 
     ! calculations of Hubbard forces and stress
     eigenval(:) = e(:)
     eigenvect(:,:) = work(:,:)
     overlap_inv(:,:) = overlap(:,:)
     !
     ! Transform atomic orbitals WITHOUT the ultrasoft S operator 
     ! O^(-1/2) \psi (note the transposition):
     ! \phi_I = \sum_J O^{-1/2}_JI \phi_J
     !
     DO i = 1, npw
        work(:,1) = (0.d0,0.d0)
        IF (noncolin) THEN
           DO ipol=1,npol
              j = i + (ipol-1)*npwx
              CALL zgemv ('n',m,m,(1.d0,0.d0),overlap, &
                   m, wfc(j,1), npwx*npol, (0.d0,0.d0),work,1)
              CALL zcopy (m,work,1,wfc(j,1),npwx*npol)
           END DO
        ELSE
           CALL zgemv ('n', m, m, (1.d0, 0.d0) , overlap, &
                m, wfc (i, 1) , npwx, (0.d0, 0.d0) , work, 1)
           CALL zcopy (m, work, 1, wfc (i, 1), npwx)
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
        work(:,1) = (0.d0,0.d0)
        IF (noncolin) THEN
           DO ipol=1,npol
              j = i + (ipol-1)*npwx
              CALL zgemv ('n',m,m,(1.d0,0.d0),overlap, &
                   m, swfc(j,1), npwx*npol, (0.d0,0.d0),work,1)
              CALL zcopy (m,work,1,swfc(j,1),npwx*npol)
           END DO
        ELSE
           CALL zgemv ('n', m, m, (1.d0, 0.d0) , overlap, &
                m, swfc (i, 1) , npwx, (0.d0, 0.d0) , work, 1)
           CALL zcopy (m, work, 1, swfc (i, 1), npwx)
        END IF
     ENDDO
     !
  ENDIF
  !
  DEALLOCATE (overlap)
  DEALLOCATE (work)
  DEALLOCATE (e)
  !
  RETURN
  !      
END SUBROUTINE ortho_swfc
!
!-----------------------------------------------------------------------
SUBROUTINE calculate_doverlap_inv (m, e, work, doverlap, doverlap_inv)
  !---------------------------------------------------------------------
  !! This routine computes the derivative of transposed O^{-1/2}, i.e.
  !! [d((O^{-1/2})^T)]_IJ, where O_IJ is the overlap matrix. 
  !! Note, on the input this routine requires dO not transposed.
  !! The solution is written in a closed form by solving the Lyapunov
  !! equation (a particular case of the Sylvester equation).
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
  !! eigenvectors of the overlap matrix
  !! auxiliary array
  !
  ALLOCATE (aux(m,m))
  !
  ! Compute (work^T) * (doverlap^T) * ((work^H)^T) 
  ! and put the result back in doverlap
  !
  ! Compute aux = (work^H)*doverlap
  CALL ZGEMM('C','N', m, m, m, (1.d0,0.d0), work, &
              m, doverlap, m, (0.d0,0.d0), aux, m)
  ! Compute (work^T) * (aux^T)
  CALL ZGEMM('T','T', m, m, m, (1.d0,0.d0), work, &
              m, aux, m, (0.d0,0.d0), doverlap, m)
  !
  DO m1 = 1, m
     DO m2 = 1, m
        aux(m1,m2) = doverlap(m1,m2) / &
                    (e(m1)*DSQRT(e(m2))+e(m2)*DSQRT(e(m1)))
     ENDDO
  ENDDO
  !
  ! Compute ((work^H)^T) * aux * (work^T)
  !
  ! Compute doverlap = (aux^T)*(work^H)
  CALL ZGEMM('T','C', m, m, m, (1.d0,0.d0), aux, &
              m, work, m, (0.d0,0.d0), doverlap, m)
  ! Compute doverlap_inv = (doverlap^T)*(work^T)
  CALL ZGEMM('T','T', m, m, m, (-1.d0,0.d0), doverlap, &
              m, work, m, (0.d0,0.d0), doverlap_inv, m)
  !
  DEALLOCATE (aux)
  !
  RETURN
  !
END SUBROUTINE calculate_doverlap_inv
