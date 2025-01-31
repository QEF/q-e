!
! Copyright (C) 2001-2007 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------------
SUBROUTINE rotate_wfc_gpu &
            ( npwx, npw, nstart, gstart, nbnd, psi, npol, overlap, evc, e )
  !----------------------------------------------------------------------------
  !
  ! ... Driver routine (maybe it should be an interface) for
  ! ... Hamiltonian diagonalization in the subspace spanned
  ! ... by nstart states psi ( atomic or random wavefunctions ).
  ! ... Produces on output nbnd eigenvectors ( nbnd <= nstart ) in evc.
  ! ... Calls h_psi, s_psi to calculate H|psi> ans S|psi>
  ! ... It only uses an auxiliary array of the same size as psi.
  !
  USE kinds,         ONLY : DP
  USE control_flags, ONLY : use_para_diag, gamma_only
  !
  IMPLICIT NONE
  !
  ! ... I/O variables
  !
  INTEGER, INTENT(IN) :: npw, npwx, nstart, nbnd, gstart, npol
    ! dimension of the matrix to be diagonalized
    ! leading dimension of matrix psi, as declared in the calling pgm unit
    ! input number of states
    ! output number of states
    ! first G with nonzero norm
    ! number of spin polarizations
  LOGICAL, INTENT(IN) :: overlap
    ! if .FALSE. : S|psi> not needed
  COMPLEX(DP), INTENT(INOUT) :: psi(npwx*npol,nstart), evc(npwx*npol,nbnd)
    ! input and output eigenvectors (may overlap)
  REAL(DP), INTENT(OUT) :: e(nbnd)
    ! eigenvalues
  COMPLEX(DP), ALLOCATABLE :: psi_h(:,:), evc_h(:,:)
  REAL(DP), ALLOCATABLE    :: e_h(:)
  EXTERNAL h_psi, s_psi, h_psi_gpu, s_psi_acc
  !
    ! h_psi(npwx,npw,nvec,psi,hpsi)
    !     calculates H|psi>
    ! s_psi(npwx,npw,nvec,spsi)
    !     calculates S|psi> (if needed)
    !     Vectors psi,hpsi,spsi are dimensioned (npwx,npol,nvec)
  !
  CALL start_clock_gpu( 'wfcrot' ); !write (*,*) 'start wfcrot' ; FLUSH(6)
  !write (*,*) 'gamma_only' , gamma_only; FLUSH(6)

  IF( use_para_diag ) THEN
     !
     ! Allocate arrays to workaround parallel case
     !
     ALLOCATE(psi_h(npwx*npol,nstart), evc_h(npwx*npol,nbnd), e_h(nbnd))
     !$acc update self(psi)
     psi_h(1:npwx*npol,1:nstart) = psi(1:npwx*npol,1:nstart)
     evc_h(1:npwx*npol,1:nbnd) = evc(1:npwx*npol,1:nbnd)
     !
     ! use data distributed subroutine
     !
     IF ( gamma_only ) THEN
  !write (*,*) 'inside para gamma'; FLUSH(6)
        !
        CALL protate_wfc_gamma ( h_psi, s_psi, overlap, &
                                 npwx, npw, nstart, nbnd, psi_h, evc_h, e_h )
        !
     ELSE
  !write (*,*) 'inside para k'; FLUSH(6)
        !
        CALL protate_wfc_k ( h_psi, s_psi, overlap, &
                             npwx, npw, nstart, nbnd, npol, psi_h, evc_h, e_h )
        !
     END IF
     psi(1:npwx*npol,1:nstart) = psi_h(1:npwx*npol,1:nstart)
     !$acc update device(psi)
     evc(1:npwx*npol,1:nbnd)   = evc_h(1:npwx*npol,1:nbnd)
     e(1:nbnd)                 = e_h(1:nbnd)
     DEALLOCATE(psi_h, evc_h, e_h)
     !
  ELSE
     !
     ! use serial subroutines
     !
     IF ( gamma_only ) THEN
  !write (*,*) 'inside serial gamma'; FLUSH(6)
        !
        CALL rotate_wfc_gamma ( h_psi_gpu, s_psi_acc, overlap, &
                                    npwx, npw, nstart, nbnd, psi, evc, e )
        !
     ELSE
  !write (*,*) 'inside serial k'; FLUSH(6)
        !
        CALL rotate_wfc_k ( h_psi_gpu, s_psi_acc, overlap, &
                                npwx, npw, nstart, nbnd, npol, psi, evc, e )
        !
     END IF
     !
  END IF
  !
  CALL stop_clock_gpu( 'wfcrot' )!; write (*,*) 'stop wfcrot' ; FLUSH(6)
  !
END SUBROUTINE rotate_wfc_gpu
