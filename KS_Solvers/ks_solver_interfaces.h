! Copyright (C) 2013-2021 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

INTERFACE rotate_xpsi
  SUBROUTINE rotate_xpsi_driver ( h_psi_ptr, s_psi_ptr, &
      npwx, npw, nstart, nbnd, psi, npol, overlap, evc, hevc, sevc, e, use_para_diag, gamma_only )
    !
    IMPORT :: DP  
    IMPLICIT NONE
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
    EXTERNAL :: h_psi_ptr, s_psi_ptr
  END SUBROUTINE rotate_xpsi_driver 
#if defined (__CUDA) 
  SUBROUTINE rotate_xpsi_driver_cuf ( h_psi_hptr, s_psi_hptr, h_psi_dptr, s_psi_dptr, &
      npwx, npw, nstart, nbnd, psi, npol, overlap, evc, hevc, sevc, e_d, use_para_diag, gamma_only )
  !! Driver routine for Hamiltonian diagonalization in the subspace 
  !! spanned by nstart states psi ( atomic or random wavefunctions ).
  !! Interface for the CUDA-Fortran case. 
  !! Produces on output nbnd eigenvectors ( nbnd <= nstart ) in evc.
  !! Calls h_psi, s_psi to calculate H|psi> and S|psi>,
  !! which are saved in hevc and sevc.
  IMPORT :: DP
  IMPLICIT NONE 
  INTEGER, INTENT(IN) :: npw, npwx, nstart, nbnd, npol
  !! dimension of the matrix to be diagonalized
  !! leading dimension of matrix psi, as declared in the calling pgm unit
  !! input number of states
  !! output number of states
  !! number of spin polarizations
  LOGICAL, INTENT(IN) :: overlap
  !! if .FALSE. : S|psi> not needed
  COMPLEX(DP), INTENT(INOUT) :: psi(npwx*npol,nstart)
  !! vectors spannign the subspace 
  COMPLEX(DP), INTENT(INOUT)   :: evc(npwx*npol,nbnd)
  !! input and output eigenvectors (may overlap)
  COMPLEX(DP), INTENT(OUT)   :: hevc(npwx*npol,nbnd), sevc(npwx*npol,nbnd)
  !! H|psi> and S|psi>
  REAL(DP),  INTENT(OUT) :: e_d(nbnd)
  !! eigenvalues
  LOGICAL, INTENT(IN) :: use_para_diag 
  !! if true, use parallel diagonalization 
  LOGICAL, INTENT(IN) :: gamma_only 
  !! set to true if H matrix is real 
  attributes(DEVICE)       :: e_d
  EXTERNAL :: h_psi_hptr, h_psi_dptr, s_psi_hptr, s_psi_dptr 
END SUBROUTINE rotate_xpsi_driver_cuf
#endif 
END INTERFACE 

  
  
