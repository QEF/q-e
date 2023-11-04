!
! Copyright (C) 2019 National Institute of Advanced Industrial Science and Technology (AIST)
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
SUBROUTINE rotate_xpsi_driver_cuf ( h_psi_hptr, s_psi_hptr, h_psi_dptr, s_psi_dptr, &
              npwx, npw, nstart, nbnd, psi_d, npol, overlap, evc_d, hevc_d, sevc_d, e_d, use_para_diag, gamma_only )
  !----------------------------------------------------------------------------
  !
  !! Driver routine for Hamiltonian diagonalization in the subspace 
  !! spanned by nstart states psi ( atomic or random wavefunctions ).
  !! Produces on output nbnd eigenvectors ( nbnd <= nstart ) in evc.
  !! Calls h_psi_ptr, s_psi_ptr to calculate H|psi> and S|psi>,
  !! which are saved in hevc and sevc.
  !
  USE util_param,         ONLY : DP
  !
  IMPLICIT NONE
  !
  ! ... I/O variables
  !
  INTEGER, INTENT(IN) :: npw, npwx, nstart, nbnd, npol
  !! dimension of the matrix to be diagonalized
  !! leading dimension of matrix psi, as declared in the calling pgm unit
  !! input number of states
  !! output number of states
  !! number of spin polarizations
  LOGICAL, INTENT(IN) :: overlap
  !! if .FALSE. : S|psi> not needed
  COMPLEX(DP), INTENT(INOUT) :: psi_d(npwx*npol,nstart)
  !! vectors spannign the subspace 
  COMPLEX(DP), INTENT(INOUT)   :: evc_d(npwx*npol,nbnd)
  !! input and output eigenvectors (may overlap)
  COMPLEX(DP), INTENT(OUT)   :: hevc_d(npwx*npol,nbnd), sevc_d(npwx*npol,nbnd)
  !! H|psi> and S|psi>
  REAL(DP),  INTENT(OUT) :: e_d(nbnd)
  !! eigenvalues
  LOGICAL, INTENT(IN) :: use_para_diag 
  !! if true, use parallel diagonalization 
  LOGICAL, INTENT(IN) :: gamma_only 
  !! set to true if H matrix is real 
#if defined(__CUDA)
  attributes(DEVICE)       :: psi_d, evc_d, hevc_d, sevc_d, e_d
#endif

  COMPLEX(DP), ALLOCATABLE         :: psi_h(:,:)
  COMPLEX(DP), ALLOCATABLE, TARGET :: evc_h(:,:)
  COMPLEX(DP), ALLOCATABLE         :: hevc_h(:,:) 
  COMPLEX(DP), POINTER             :: sevc_h(:,:)
  REAL(DP), ALLOCATABLE            :: e_h(:)
  !
  EXTERNAL :: h_psi_hptr, h_psi_dptr, &  ! host pointers
              s_psi_hptr, s_psi_dptr     ! device pointers
    ! h_psi_... (npwx,npw,nbnd,psi,hpsi)
    !     calculates H|psi>
    ! s_psi_... (npwx,npw,nbnd,spsi)
    !     calculates S|psi> (if needed)
    !     Vectors psi,hpsi,spsi are dimensioned (npwx,npol,nbnd)
  !
  CALL start_clock_gpu( 'wfcrot' ); !write (*,*) 'start wfcrot' ; FLUSH(6)
  !write (*,*) 'gamma_only' , gamma_only; FLUSH(6)
  !
  IF( use_para_diag ) THEN
     !
     !Allocate arrays to workaround parallel case
     !
     ALLOCATE(psi_h(npwx*npol,nstart), evc_h(npwx*npol,nbnd), hevc_h(npwx*npol,nbnd), &
              e_h(nbnd))
     IF(overlap) THEN 
        ALLOCATE(sevc_h(npwx*npol,nbnd))
     ELSE
        sevc_h => evc_h
     END IF
     !
     psi_h(1:npwx*npol,1:nstart) = psi_d(1:npwx*npol,1:nstart)
     evc_h(1:npwx*npol,1:nbnd)   = evc_d(1:npwx*npol,1:nbnd)
     !
     ! use data distributed subroutine
     !
     IF ( gamma_only ) THEN
  !write (*,*) 'inside para gamma'; FLUSH(6)
        !
        call protate_xpsi_gamma ( h_psi_hptr, s_psi_hptr, overlap, &
                                  npwx, npw, nstart, nbnd, psi_h, evc_h, hevc_h, sevc_h, e_h )
        !
     ELSE
  !write (*,*) 'inside para k'; FLUSH(6)
        !
        call protate_xpsi_k ( h_psi_hptr, s_psi_hptr, overlap, &
                              npwx, npw, nstart, nbnd, npol, psi_h, evc_h, hevc_h, sevc_h, e_h )
        !
     END IF
     psi_d(1:npwx*npol,1:nstart) = psi_h(1:npwx*npol,1:nstart)
     evc_d(1:npwx*npol,1:nbnd)   = evc_h(1:npwx*npol,1:nbnd)
     hevc_d(1:npwx*npol,1:nbnd)  = hevc_h(1:npwx*npol,1:nbnd)
     e_d(1:nbnd)                 = e_h(1:nbnd)
     !
     DEALLOCATE(psi_h, evc_h, hevc_h, e_h)
     IF(overlap) THEN 
       sevc_d(1:npwx*npol,1:nbnd)  = sevc_h(1:npwx*npol,1:nbnd)
       DEALLOCATE(sevc_h)
     ELSE
        NULLIFY(sevc_h)
     END IF
     !
  ELSE
     !
     ! use serial subroutines
     !
     IF ( gamma_only ) THEN
  !write (*,*) 'inside serial gamma'; FLUSH(6)
        !
        CALL rotate_xpsi_gamma_gpu ( h_psi_dptr, s_psi_dptr, overlap, &
                                 npwx, npw, nstart, nbnd, psi_d, evc_d, hevc_d, sevc_d, e_d )
        !
     ELSE
  !write (*,*) 'inside serial k'; FLUSH(6)
        !
        CALL rotate_xpsi_k_gpu ( h_psi_dptr, s_psi_dptr, overlap, &
                             npwx, npw, nstart, nbnd, npol, psi_d, evc_d, hevc_d, sevc_d, e_d )
        !
     END IF
     !
  END IF
  !
  CALL stop_clock_gpu( 'wfcrot' )
  !
END SUBROUTINE rotate_xpsi_driver_cuf
