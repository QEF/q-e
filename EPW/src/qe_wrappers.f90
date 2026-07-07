  !
  ! Copyright (C) 2023-2026 EPW-Collaboration
  ! Copyright (C) 2024 EPW-Collaboration
  !
  ! This file is distributed under the terms of the GNU General Public
  ! License. See the file `LICENSE' in the root directory of the
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .
  !
  !----------------------------------------------------------------------
  MODULE qe_wrappers
  !----------------------------------------------------------------------
  !!
  !! EPW wrappers to use Quantum Espresso subroutines
  !!
  IMPLICIT NONE
  !
  CONTAINS
    !-----------------------------------------------------------------------
    SUBROUTINE adddvepsi_us_epw(becp1, becp2, ipol, ik, dvpsi)
    !-----------------------------------------------------------------------  
    ! AACA: LSDA implementation requires a wrapper on the QE surbroutine adddvepsi_us
    ! found in LR_Modules in the file with the same name beacuase of the different scheemes 
    ! of parallelization in QE and in wannier90. Mainly ngk dimensions are not correct in 
    ! the LSDA case instead they are changed with the correct dimensions and values and then 
    ! the QE subroutine is called
    !
    USE kinds,           ONLY : DP
    USE pwcom,           ONLY : npwx, nbnd
    USE input,           ONLY : ngk_loc
    USE klist,           ONLY : ngk
    USE becmod,          ONLY : bec_type
    USE noncollin_module,ONLY : npol
    USE control_lr,      ONLY : nbnd_occ
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: ipol
    !! Number of polarizations
    INTEGER, INTENT(in) :: ik
    !! k point index
    TYPE(bec_type), INTENT(in) :: becp1
    !!  eigenvector to apply  dimensions ( nkb, nbnd )
    TYPE(bec_type), INTENT(in) :: becp2
    !! vector 2 dimensions ( nkb, nbnd )
    COMPLEX(KIND=DP), INTENT(inout) :: dvpsi(npwx*npol,nbnd)
    !! Derivative of the KS potential
    !
    ! Local variables
    !
    INTEGER(KIND = DP) :: ngk_store
    !! Store QE native ngk at poit ik
    INTEGER(KIND = DP) :: nbnd_store
    !! Store QE native nbnd_occ at point ik
    !
    ! Make a copy of the original values 
    ngk_store = ngk(ik)
    nbnd_store = nbnd_occ(ik)
    !
    ! Copy the values of the EPW native arrays into QE native arrays
    ngk(ik) = ngk_loc(ik)
    nbnd_occ(ik) = nbnd
    !
    CALL adddvepsi_us(becp1, becp2, ipol, ik, dvpsi)
    !
    ! Resotre the correct values
    ngk(ik) = ngk_store
    nbnd_occ(ik) = nbnd_store
    !
    RETURN
    !-----------------------------------------------------------------------
    END SUBROUTINE adddvepsi_us_epw
    !-----------------------------------------------------------------------
    !
    !------------------------------------------------------------------------
    SUBROUTINE commutator_Hx_psi_epw(ik, nbnd_calc, vpol, becp1, becp2, dpsi)
    !------------------------------------------------------------------------
    ! AACA: LSDA implementation requires a wrapper on the QE surbroutine commutator_Hx_psi
    ! foind in PW/src in the file with the same name beacuase of the different scheemes 
    ! of parallelization in QE and in wannier90. Mainly ngk, xk, igk_k and et dimensions are 
    ! not correct in  the LSDA case instead they are changed with the correct dimensions and  
    ! values and then the QE subroutine is called
    !
    USE kinds,           ONLY : DP
    USE pwcom,           ONLY : npwx
    USE wvfct,           ONLY : nbnd, et
    USE noncollin_module,ONLY : npol
    USE becmod,          ONLY : bec_type
    USE klist,           ONLY : ngk, igk_k, xk
    USE input,           ONLY : ngk_loc, igk_k_loc, et_loc, xk_loc, isk_loc
    USE lsda_mod,        ONLY : isk
  
    IMPLICIT NONE
    !
    COMPLEX(DP), INTENT(inout)    :: dpsi(npwx*npol, nbnd_calc)
    !! Derivative of the wave function
    TYPE(bec_type), INTENT(in)  :: becp1 
    !! eigenvector to apply  dimensions ( nkb, nbnd )
    TYPE(bec_type), INTENT(inout) :: becp2 
    !! vector 2 dimensions ( nkb, nbnd )
    INTEGER, INTENT(in) :: ik
    !! k point index
    INTEGER, INTENT(in) :: nbnd_calc
    !! number of bands
    REAL(DP), INTENT(in) :: vpol(3)
    !! polarization vector in Cartesian coordinates
    !
    ! Local variables
    !
    INTEGER :: ierr
    !! Error
    INTEGER(KIND = DP) :: ngk_store
    !! Store ngk QE native variable
    INTEGER(KIND = DP) :: isk_store
    !! Store isk QE native variable
    INTEGER(KIND = DP), ALLOCATABLE :: igk_k_store(:)
    !! Store igk_k QE native variable
    REAL(KIND = DP) xk_store(3)
    !! Store xk QE native variable
    REAL(KIND = DP), ALLOCATABLE :: et_store(:)
    !! Store et QE native variable
    !
    ! Create the arrays to save the original values
    ALLOCATE(igk_k_store(npwx), STAT = ierr)
    IF (ierr /= 0) CALL errore('commutator_Hx_psi_epw', 'Error allocating igk_k_store', 1)
    ALLOCATE(et_store(nbnd), STAT = ierr)
    IF (ierr /= 0) CALL errore('commutator_Hx_psi_epw', 'Error allocating et_store', 1)
    !
    ngk_store = ngk(ik)
    isk_store = isk(ik)
    igk_k_store(:) = igk_k(:, ik)
    xk_store(:) = xk(:, ik)
    et_store(:) = et(:, ik)
    !
    ! Change the content of the QE native variables
    !
    ngk(ik) = ngk_loc(ik)
    isk(ik) = isk_loc(ik)
    igk_k(:, ik) = igk_k_loc(:, ik)
    xk(:, ik) = xk_loc(:, ik)
    et(:, ik) = et_loc(:, ik)
    !
    !$acc update device(igk_k(:,ik),et(:,ik))
    CALL commutator_Hx_psi(ik, nbnd_calc, vpol, becp1, becp2, dpsi)
    !
    ! Restore the original values to the QE native variables  
    !
    ngk(ik) = ngk_store
    isk(ik) = isk_store
    igk_k(:, ik) = igk_k_store(:)
    xk(:, ik) = xk_store(:)
    et(:, ik) = et_store(:)
    !
    ! Deallocate the store arrays.
    DEALLOCATE(igk_k_store, STAT = ierr)
    IF (ierr /= 0) CALL errore('commutator_Hx_psi_epw', 'Error deallocating igk_k_store', 1)
    DEALLOCATE(et_store, STAT = ierr)
    IF (ierr /= 0) CALL errore('commutator_Hx_psi_epw', 'Error deallocating et_store', 1)
    !
    RETURN
    END SUBROUTINE commutator_Hx_psi_epw
    !-----------------------------------------------------------------------------
  !----------------------------------------------------------------------
  END MODULE qe_wrappers
  !----------------------------------------------------------------------
