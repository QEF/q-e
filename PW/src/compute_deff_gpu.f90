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
  USE ions_base,   ONLY: nat !, nsp, ityp
  USE uspp,        ONLY: okvan
  USE uspp_param,  ONLY: nhm
  USE lsda_mod,    ONLY: current_spin
  !
  USE uspp_gpum, ONLY : using_deeq_d, deeq_d, &
                        using_qq_at_d, qq_at_d
  !
  !USE device_util_m,  ONLY : dev_memcpy
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
  INTEGER :: nt, na, is, i ,j
  !
#if defined(__CUDA)
  attributes(DEVICE) ::  deff_d
#endif  
  !
  !call using_deeq(0) ; call using_qq_at(0)
  CALL using_deeq_d(0)
  CALL using_qq_at_d(0)
  !
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
! SUBROUTINE compute_deff_nc( deff, et )
!   !-------------------------------------------------------------------------
!   !! This routine computes the effective value of the D-eS coefficients
!   !! which appears often in many expressions. This routine is for the
!   !! noncollinear case.
!   !
!   USE kinds,            ONLY: DP
!   USE ions_base,        ONLY: nsp, nat, ityp
!   USE spin_orb,         ONLY: lspinorb
!   USE noncollin_module, ONLY: noncolin, npol
!   USE uspp,             ONLY: deeq_nc, qq_at, qq_so, okvan
!   USE uspp_param,       ONLY: nhm
!   USE lsda_mod,         ONLY: nspin
!   !
!   USE uspp_gpum, ONLY : using_deeq_nc, using_qq_at, using_qq_so
!   !
!   IMPLICIT NONE
!   !
!   REAL(DP), INTENT(IN) :: et
!   !! The eigenvalues of the hamiltonian
!   COMPLEX(DP), INTENT(OUT) :: deff(nhm,nhm,nat,nspin) 
!   !! Effective values of the D-eS coefficients
!   !
!   ! ... local variables
!   !
!   INTEGER :: nt, na, is, js, ijs
!   !
!   CALL using_deeq_nc(0)
!   IF (.not. lspinorb) CALL using_qq_at(0)
!   IF (lspinorb) CALL using_qq_so(0)
!   !
!   deff=deeq_nc
!   IF (okvan) THEN
!      !
!      DO nt = 1, nsp
!         DO na = 1, nat
!            !
!            IF ( ityp(na) == nt ) THEN
!               IF (lspinorb) THEN
!                  deff(:,:,na,:) = deff(:,:,na,:) - et * qq_so(:,:,:,nt)
!               ELSE
!                  ijs=0
!                  !
!                  DO is=1,npol
!                     DO js=1,npol
!                        !
!                        ijs=ijs+1
!                        IF (is==js) deff(:,:,na,ijs)=deff(:,:,na,ijs)-et*qq_at(:,:,na)
!                        !
!                     ENDDO
!                  ENDDO
!                  !
!               ENDIF
!            ENDIF
!            !
!         ENDDO
!      ENDDO
!      !
!   ENDIF
!   !
!   !
!   RETURN
!   !
! END SUBROUTINE compute_deff_nc
