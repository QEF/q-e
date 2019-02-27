!
! Copyright (C) 2001-2015 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-------------------------------------------------------------------------
SUBROUTINE lr_sym_eels (dvtosym)
  !-----------------------------------------------------------------------
  !
  ! This subroutine symmetrizes the response charge-density.
  ! Inspired by PH/symdvscf.f90
  !
  ! TODO: Try to use the rotuine PH/symdvscf.f90 directly.
  !
  ! Created by Iurii Timrov (2013)
  !
  USE kinds,            only : DP
  USE constants,        ONLY : tpi
  USE fft_base,         ONLY : dfftp
  USE cell_base,        ONLY : at
  USE symm_base,        ONLY : s, ftau
  USE noncollin_module, ONLY : nspin_lsda, nspin_mag

  USE lr_symm_base, ONLY : minus_q, nsymq, irotmq, gi, gimq
  !  
  IMPLICIT NONE
  !
  COMPLEX(DP) :: dvtosym(dfftp%nr1x, dfftp%nr2x, dfftp%nr3x, nspin_mag)
  ! the charge density response to be symmetrized
  INTEGER :: is, ri, rj, rk, i, j, k, ipol, isym, irot
  ! counters
  REAL(DP) :: gf(3), n(3)
  ! temp variables
  COMPLEX(DP), ALLOCATABLE :: dvsym(:,:,:)
  ! the symmetrized charge density response
  COMPLEX(DP) ::  aux2, term(3, 48), phase(48)
  ! auxiliary space
  ! the multiplication factor
  ! the phase factor
  !
  IF (nsymq==1) RETURN
  !if (nsymq==1 .and. (.not.minus_q) ) return
  !
  CALL start_clock ('lr_sym_eels')
  !
  ALLOCATE(dvsym(dfftp%nr1x, dfftp%nr2x, dfftp%nr3x))
  !
  n(1) = tpi / DBLE (dfftp%nr1)
  n(2) = tpi / DBLE (dfftp%nr2)
  n(3) = tpi / DBLE (dfftp%nr3)
  !
  !------------------------------------------------------------------------!
  !  If necessary, symmetrize with respect to the sym.op.  S*q = -q + G    !
  !------------------------------------------------------------------------!
  !
  !if (minus_q) then
  !   !
  !   gf(:) =  gimq(1) * at(1,:) * n(:) + &
  !            gimq(2) * at(2,:) * n(:) + &
  !            gimq(3) * at(3,:) * n(:)
  !   !
  !   term(:,1) = CMPLX(cos(gf(:)), sin(gf(:)), kind=DP)
  !   !
  !   do is = 1, nspin_lsda
  !      !
  !      phase(1) = (1.d0, 0.d0)
  !      !
  !      do k = 1, dfftp%nr3
  !         do j = 1, dfftp%nr2
  !            do i = 1, dfftp%nr1
  !               !
  !               ! Rotation and fractional translation: S^-1 * r - ftau
  !               !
  !               call ruotaijk (s(1,1,irotmq), ftau(1,irotmq), i, j, k, &
  !               dfftp%nr1, dfftp%nr2, dfftp%nr3, ri, rj, rk)
  !               !
  !               ! drho(S^-1 * r - ftau) * exp(i G*a) 
  !               !
  !               aux2 = (0.d0, 0.d0)
  !               !
  !               aux2 = aux2 + dvtosym(ri,rj,rk,is) * phase(1)
  !               !
  !               ! 1/2 * [ drho(r) + conjg{ drho(S^-1 * r - ftau)*exp(i G*a) } ]
  !               !
  !               dvsym(i,j,k) = ( dvtosym(i,j,k,is) + CONJG(aux2) ) * 0.5d0
  !               !
  !               phase (1) = phase (1) * term (1, 1)
  !               !
  !            enddo
  !            phase (1) = phase (1) * term (2, 1)
  !         enddo
  !         phase (1) = phase (1) * term (3, 1)
  !      enddo
  !      !
  !      dvtosym(:,:,:,is) = dvsym(:,:,:)
  !      !
  !   enddo
  !   !
  !endif
  !
  !----------------------------------------------------------------!
  ! Symmetrize with respect to the small group of q : S*q = q + G  !
  !----------------------------------------------------------------!
  !
  ! Calculation of the phase exp(i G*r).
  !
  DO isym = 1, nsymq
     !
     gf(:) =  gi(1,isym) * at(1, :) * n(:) + &
              gi(2,isym) * at(2, :) * n(:) + &
              gi(3,isym) * at(3, :) * n(:)
     !
     term(:,isym) = CMPLX(cos(gf(:)), sin(gf(:)), kind=DP)
     !
  ENDDO
  !
  DO is = 1, nspin_lsda
     !
     dvsym(:,:,:) = (0.d0, 0.d0)
     !
     DO isym = 1, nsymq
        phase(isym) = (1.d0, 0.d0)
     ENDDO
     !
     DO k = 1, dfftp%nr3
        DO j = 1, dfftp%nr2
           DO i = 1, dfftp%nr1
              !
              ! Loop on the symmetry operations of the small group of q.
              !
              DO isym = 1, nsymq
                 !
                 ! Rotation and fractional translation: S^-1 * r - ftau
                 !
                 CALL ruotaijk (s(1,1,isym), ftau(1,isym), i, j, k, &
                 dfftp%nr1, dfftp%nr2, dfftp%nr3, ri, rj, rk)
                 !
                 ! Calculate drho(S^-1 * r - ftau) * exp(i G*r)
                 !
                 dvsym(i,j,k) = dvsym(i,j,k) + dvtosym(ri,rj,rk,is) * phase(isym)
                 !
              ENDDO
              !
              DO isym = 1, nsymq
                 phase (isym) = phase (isym) * term (1, isym)
              ENDDO
              !
           ENDDO
           !
           DO isym = 1, nsymq
              phase (isym) = phase (isym) * term (2, isym)
           ENDDO
           !
        ENDDO
        !
        DO isym = 1, nsymq
           phase (isym) = phase (isym) * term (3, isym)
        ENDDO
        !
     ENDDO
     !
     ! Final normalization on the number of the symmetry operations
     ! in the small group of q.
     !
     dvtosym(:,:,:,is) = dvsym(:,:,:) / DBLE (nsymq)
     !
  ENDDO
  !
  DEALLOCATE(dvsym)
  !
  CALL stop_clock ('lr_sym_eels')
  !
  RETURN
  !
END SUBROUTINE lr_sym_eels
