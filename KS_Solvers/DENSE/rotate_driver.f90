!
! Copyright (C) 2019 National Institute of Advanced Industrial Science and Technology (AIST)
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
SUBROUTINE rotate_xpsi_driver ( h_psi_hptr,  s_psi_hptr, h_psi_dptr,  s_psi_dptr, &
              npwx, npw, nstart, nbnd, psi, npol, overlap, evc, hevc, sevc, e, use_para_diag, gamma_only )
  !----------------------------------------------------------------------------
  !
  ! ... Driver routine (maybe it should be an interface) for
  ! ... Hamiltonian diagonalization in the subspace spanned
  ! ... by nstart states psi ( atomic or random wavefunctions ).
  ! ... Produces on output nbnd eigenvectors ( nbnd <= nstart ) in evc.
  ! ... Calls h_psi_ptr, s_psi_ptr to calculate H|psi> and S|psi>,
  ! ... which are saved in hevc and sevc.
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
  COMPLEX(DP), INTENT(INOUT) :: psi(npwx*npol,nstart)
  !! vectors spanning the subspace 
  COMPLEX(DP), INTENT(OUT)   :: evc(npwx*npol,nbnd)
  !! input and output eigenvectors (may overlap)
  COMPLEX(DP), INTENT(OUT) :: hevc(npwx*npol,nbnd), sevc(npwx*npol,nbnd)
  !! H|psi> and S|psi>
  REAL(DP), INTENT(OUT) :: e(nbnd)
  !! eigenvalues
  LOGICAL,INTENT(IN) :: use_para_diag 
  !! if true parallel diagonalization will be used 
  LOGICAL,INTENT(IN) :: gamma_only
  !! set to true when H  is real 

  !
  EXTERNAL :: h_psi_hptr, s_psi_hptr, h_psi_dptr, s_psi_dptr
    ! [hptr, dptr] --> [host, device] pointers (without GPU dptr = hptr)
    ! h_psi_ptr(npwx,npw,nbnd,psi,hpsi)
    !     calculates H|psi>
    ! s_psi_ptr(npwx,npw,nbnd,spsi)
    !     calculates S|psi> (if needed)
    !     Vectors psi,hpsi,spsi are dimensioned (npwx,npol,nbnd)
  !
  CALL start_clock( 'wfcrot' )
  !
  IF( use_para_diag ) THEN
     !
     ! use data distributed subroutine
     !
     !$acc update host(psi, evc, sevc)
     !$civn: the following subroutines work on CPU
     IF ( gamma_only ) THEN
        !
        CALL protate_xpsi_gamma ( h_psi_hptr, s_psi_hptr, overlap, &
                                  npwx, npw, nstart, nbnd, psi, evc, hevc, sevc, e )
        !
     ELSE
        !
        CALL protate_xpsi_k ( h_psi_hptr, s_psi_hptr, overlap, &
                              npwx, npw, nstart, nbnd, npol, psi, evc, hevc, sevc, e )
        !
     END IF
     !$acc update device(psi, evc, hevc, sevc, e)
     !
  ELSE
     !
     ! use serial subroutines
     !
     !$civn: the following subroutines work on GPU 
     IF ( gamma_only ) THEN
        !
        CALL rotate_xpsi_gamma ( h_psi_dptr, s_psi_dptr, overlap, &
                                 npwx, npw, nstart, nbnd, psi, evc, hevc, sevc, e )
        !
     ELSE
        !
        CALL rotate_xpsi_k ( h_psi_dptr, s_psi_dptr, overlap, &
                             npwx, npw, nstart, nbnd, npol, psi, evc, hevc, sevc, e )
        !
     END IF
     !
  END IF
  !
  CALL stop_clock( 'wfcrot' )
  !
END SUBROUTINE rotate_xpsi_driver

