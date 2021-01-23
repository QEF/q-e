!
! Copyright (C) 2009-2010 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-------------------------------------------------------------------------
SUBROUTINE compute_deff_gpu( deff_d, et )
  !-----------------------------------------------------------------------
  !! This routine computes the effective value of the D-eS coefficients
  !! which appear often in many expressions in the US or PAW case. 
  !! This routine is for the collinear case.
  !
  USE kinds,       ONLY: DP
  USE ions_base,   ONLY: nat
  USE uspp,        ONLY: okvan
  USE uspp_param,  ONLY: nhm
  USE lsda_mod,    ONLY: current_spin
  !
  USE uspp_gpum,   ONLY: using_deeq_d, deeq_d, &
                         using_qq_at_d, qq_at_d
  !
  IMPLICIT NONE
  !
  REAL(DP), INTENT(IN) :: et
  !! The eigenvalues of the hamiltonian
  COMPLEX(DP), INTENT(OUT) :: deff_d(nhm,nhm,nat)
  !! Effective values of the D-eS coefficients
  !
  ! ... local variables
  !
  INTEGER :: nt, na, is, i, j
  !
#if defined(__CUDA)
  attributes(DEVICE) ::  deff_d
#endif  
  !
  CALL using_deeq_d(0)
  CALL using_qq_at_d(0)
  !
  IF (.NOT. okvan) THEN
     !
     !$cuf kernel do (3)  <<<*,*>>>
     DO na = 1, nat
       DO i = 1, nhm
         DO j = 1,nhm
           deff_d(i,j,na) = CMPLX(deeq_d(i,j,na,current_spin))
         ENDDO
       ENDDO
     ENDDO
     !
  ELSE
     !
     !$cuf kernel do (3)  <<<*,*>>>
     DO na = 1, nat
       DO i = 1, nhm
         DO j = 1,nhm
           deff_d(i,j,na) = CMPLX(deeq_d(i,j,na,current_spin) - et*qq_at_d(i,j,na))
         ENDDO
       ENDDO
     ENDDO
     !
  ENDIF
  !
  !
  RETURN
  !
END SUBROUTINE compute_deff_gpu
!
!
!---------------------------------------------------------------------------
SUBROUTINE compute_deff_nc_gpu( deff_d, et )
  !-------------------------------------------------------------------------
  !! This routine computes the effective value of the D-eS coefficients
  !! which appears often in many expressions. This routine is for the
  !! noncollinear case.
  !
  USE kinds,            ONLY: DP
  USE ions_base,        ONLY: nsp, nat, ityp
  USE spin_orb,         ONLY: lspinorb
  USE noncollin_module, ONLY: noncolin, npol
  USE uspp,             ONLY: okvan
  USE uspp_param,       ONLY: nhm
  USE lsda_mod,         ONLY: nspin
  !
  USE uspp_gpum,        ONLY: using_deeq_nc_d, using_qq_at_d,    &
                              using_qq_so_d, deeq_nc_d, qq_so_d, &
                              qq_at_d
  USE device_memcpy_m,    ONLY: dev_memcpy
  !
  IMPLICIT NONE
  !
  REAL(DP), INTENT(IN) :: et
  !! The eigenvalues of the hamiltonian
  COMPLEX(DP), INTENT(OUT) :: deff_d(nhm,nhm,nat,nspin) 
  !! Effective values of the D-eS coefficients
  !
  ! ... local variables
  !
  INTEGER :: nt, na, is, js, ijs, i, j, ias
  INTEGER :: na_v(nat), nt_v(nat)
  INTEGER, ALLOCATABLE :: na_d(:), nt_d(:)
  !
#if defined(__CUDA)
  attributes(DEVICE) ::  deff_d, na_d, nt_d
#endif  
  !
  CALL using_deeq_nc_d(0)
  IF (.NOT. lspinorb) CALL using_qq_at_d(0)
  IF (lspinorb) CALL using_qq_so_d(0)
  !
  CALL dev_memcpy( deff_d, deeq_nc_d )
  !
  IF ( okvan ) THEN
    !
    ALLOCATE( nt_d(nat), na_d(nat) )
    !
    i = 0
    DO nt = 1, nsp
      DO na = 1, nat
        IF ( ityp(na)/=nt ) CYCLE
        i = i + 1
        nt_v(i) = nt
        na_v(i) = na
      ENDDO
    ENDDO
    !
    nt_d = nt_v ; na_d = na_v
    !
    IF (lspinorb) THEN
      !
      !$cuf kernel do (3) <<<*,*>>>
      DO ias = 1, nat
        DO i = 1, nhm
          DO j = 1, nhm
            na = na_d(ias)
            nt = nt_d(ias)
            deff_d(i,j,na,:) = deeq_nc_d(i,j,na,:) - et*qq_so_d(i,j,:,nt)
          ENDDO
        ENDDO
      ENDDO
      !
    ELSE  
      !
      !$cuf kernel do (3) <<<*,*>>>
      DO ias = 1, nat
        DO i = 1, nhm
          DO j = 1, nhm
            na = na_d(ias)
            nt = nt_d(ias)
            DO is = 1, npol
              ijs = (is-1)*npol+is
              deff_d(i,j,na,ijs) = deeq_nc_d(i,j,na,ijs) - &
                                                 et*qq_at_d(i,j,na)
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      !
    ENDIF
    !
    DEALLOCATE( nt_d, na_d )
    !
  ENDIF  
  !
  !
  RETURN
  !
END SUBROUTINE compute_deff_nc_gpu
