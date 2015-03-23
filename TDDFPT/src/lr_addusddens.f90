!
! Copyright (C) 2001-2015 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------------
SUBROUTINE lr_addusddens (drhoscf, dbecsum)
  !---------------------------------------------------------------------------
  !
  ! Calculate the additional charge in reciprocal space due to US PP's
  ! See Eq.(36) in B. Walker and R. Gebauer, J. Chem. Phys. 127, 164106 (2007)
  ! Then sum up the normal and ultrasoft charges.
  ! It assumes that the array dbecsum has already been computed.
  ! Inspired by PH/addusddens.f90
  !
  ! Created by Iurii Timrov (2013)
  !
  USE kinds,                ONLY : DP
  USE ions_base,            ONLY : nat, ityp, tau, ntyp => nsp
  USE fft_base,             ONLY : dfftp
  USE fft_interfaces,       ONLY : invfft
  USE gvect,                ONLY : gg, ngm, nl, g, eigts1, eigts2, eigts3, mill
  USE uspp,                 ONLY : okvan
  USE wavefunctions_module, ONLY : psic
  USE uspp_param,           ONLY : upf, lmaxq, nh, nhm
  USE paw_variables,        ONLY : okpaw
  USE qpoint,               ONLY : xq, eigqts
  USE noncollin_module,     ONLY : nspin_mag
  !
  IMPLICIT NONE
  !
  COMPLEX(DP), INTENT(inout) :: drhoscf(dfftp%nnr, nspin_mag)
  ! input/output : change of the charge density
  COMPLEX(DP), INTENT(in)    :: dbecsum(nhm*(nhm+1)/2, nat, nspin_mag)
  ! input : the ultrasoft term
  !
  ! the local variables
  !
  INTEGER :: ig, na, nt, ih, jh, is, ijh, ir
  ! counter on G vectors
  ! counter on atoms
  ! counter on atomic type
  ! counter on beta functions
  ! counter on beta functions
  ! counter on r vectors
  ! counter on spin
  ! counter on combined beta functions
  !
  REAL(DP), ALLOCATABLE :: qmod(:), qpg(:,:), ylmk0(:,:)
  ! the modulus of q+G
  ! the values of q+G
  ! the spherical harmonics
  !
  COMPLEX(DP), ALLOCATABLE :: sk(:), qgm(:), aux(:,:)
  ! the structure factor
  ! q_lm(G)
  ! auxiliary variable for drho(G)
  !
  IF (.NOT.okvan) RETURN
  !
  CALL start_clock ('lr_addusddens')
  !
  ALLOCATE (aux(ngm,nspin_mag))
  ALLOCATE (sk(ngm))
  ALLOCATE (ylmk0(ngm,lmaxq * lmaxq))
  ALLOCATE (qgm(ngm))
  ALLOCATE (qmod(ngm))
  ALLOCATE (qpg(3,ngm))
  !
  aux(:,:) = (0.d0, 0.d0)
  !
  ! Calculate the q+G vector, its modulus, and the spherical harmonics.
  !
  CALL setqmod (ngm, xq, g, qmod, qpg)
  !
  CALL ylmr2 (lmaxq * lmaxq, ngm, qpg, qmod, ylmk0)
  !
  DO ig = 1, ngm
     qmod(ig) = sqrt(qmod(ig))
  ENDDO
  !
  DO nt = 1, ntyp
     IF (upf(nt)%tvanp) THEN
        ijh = 0
        DO ih = 1, nh (nt)
           DO jh = ih, nh (nt)
              !
              ! Calculate the Fourier transform of the Q functions,
              ! and put the result in qgm.
              !
              CALL qvan2 (ngm, ih, jh, nt, qmod, qgm, ylmk0)
              !
              ijh = ijh + 1
              DO na = 1, nat
                 IF (ityp (na) .eq.nt) THEN
                    !
                    ! Calculate the second term in Eq.(36) of the ultrasoft paper.
                    !
                    DO is = 1, nspin_mag
                       DO ig = 1, ngm
                          !
                          ! Calculate the structure factor
                          !
                          sk(ig) = eigts1(mill(1,ig),na) * &
                                   eigts2(mill(2,ig),na) * &
                                   eigts3(mill(3,ig),na) * &
                                   eigqts(na) 
                          !
                          aux(ig,is) = aux(ig,is) + 2.0d0 * qgm(ig) * sk(ig) * dbecsum(ijh,na,is)
                          !
                       ENDDO
                    ENDDO
                    !
                 ENDIF
              ENDDO
           ENDDO
        ENDDO
     ENDIF
  ENDDO
  !
  ! Convert aux to real space, and add to the charge density.
  !
  DO is = 1, nspin_mag
      !
      psic(:) = (0.d0, 0.d0)
      !
      DO ig = 1, ngm
         psic(nl(ig)) = aux(ig,is)
      ENDDO
      !
      CALL invfft ('Dense', psic, dfftp)
      !
      DO ir = 1, dfftp%nnr
         drhoscf(ir,is) = drhoscf(ir,is) + psic(ir) 
      ENDDO
      !
  ENDDO
  !
  DEALLOCATE (qpg)
  DEALLOCATE (qmod)
  DEALLOCATE (qgm)
  DEALLOCATE (ylmk0)
  DEALLOCATE (sk)
  DEALLOCATE (aux)
  !
  CALL stop_clock ('lr_addusddens')
  !
  RETURN
  !
END SUBROUTINE lr_addusddens
