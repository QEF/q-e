  !
  ! Copyright (C) 2010-2016 Samuel Ponce', Roxana Margine, Carla Verdi,
  ! Feliciano Giustino
  !
  ! Copyright (C) 2001-2003 PWSCF group
  ! This file is distributed under the terms of the
  ! GNU General Public License. See the file `License'
  ! in the root directory of the present distribution,
  ! or http://www.gnu.org/copyleft/gpl.txt .
  !
  ! Adapted from LR_Modules/newdq.f90 (QE)
  !
  !----------------------------------------------------------------------
  SUBROUTINE newdq2(dvscf, npe, xq0, timerev)
  !----------------------------------------------------------------------
  !!
  !! This routine computes the contribution of the selfconsistent
  !! change of the potential to the known part of the linear
  !! system and adds it to dvpsi.
  !!
  !! Roxana Margine - Jan 2019: Updated based on QE 6.3
  !!
  USE kinds,                ONLY : DP
  USE ions_base,            ONLY : nat, ityp, ntyp => nsp
  USE noncollin_module,     ONLY : noncolin, nspin_mag
  USE cell_base,            ONLY : omega
  USE fft_base,             ONLY : dfftp
  USE fft_interfaces,       ONLY : fwfft
  USE gvect,                ONLY : g, gg, ngm, mill, eigts1, eigts2, eigts3
  USE uspp,                 ONLY : okvan
  USE uspp_param,           ONLY : upf, lmaxq, nh, nhm
  USE mp_global,            ONLY : intra_pool_comm
  USE mp,                   ONLY : mp_sum
  USE lrus,                 ONLY : int3
  USE qpoint,               ONLY : eigqts
  USE constants_epw,        ONLY : czero, zero
  !
  IMPLICIT NONE
  !
  LOGICAL, INTENT(in) :: timerev
  !!  true if we are using time reversal
  INTEGER, INTENT(in) :: npe
  !! Number of perturbations for this irr representation
  REAL(KIND = DP), INTENT(in) :: xq0(3)
  !! The first q-point in the star (cartesian coords.)
  COMPLEX(KIND = DP), INTENT(in) :: dvscf(dfftp%nnr, nspin_mag, npe)
  !! Change of the selfconsistent potential
  !
  ! Local variables
  INTEGER :: na
  !! counter on atoms
  INTEGER :: ig
  !! counter on G vectors
  INTEGER :: nt
  !! counter on atomic types
  INTEGER ::  ir
  !! counter on real mesh
  INTEGER :: ipert
  !! counter on change of Vscf due to perturbations
  INTEGER :: is
  !! counter on spin
  INTEGER :: ih
  !! Counter on beta functions
  INTEGER :: jh
  !! Counter on beta functions
  INTEGER :: ierr
  !! Error status
  REAL(KIND = DP), ALLOCATABLE :: qmod(:)
  !! the modulus of q+G
  REAL(KIND = DP), ALLOCATABLE :: qg(:, :)
  !! the values of q+G
  REAL(KIND = DP), ALLOCATABLE :: ylmk0(:, :)
  !! the spherical harmonics at q+G
  COMPLEX(KIND = DP), EXTERNAL :: ZDOTC
  !! the scalar product function
  COMPLEX(KIND = DP), ALLOCATABLE :: aux1(:), aux2(:, :)
  !! Auxillary variable
  COMPLEX(KIND = DP), ALLOCATABLE :: qgm(:)
  !! the augmentation function at q+G
  COMPLEX(KIND = DP), ALLOCATABLE :: veff(:)
  !! effective potential
  !
  IF (.NOT. okvan) RETURN
  !
  CALL start_clock('newdq2')
  !
  int3(:, :, :, :, :) = czero
  ALLOCATE(aux1(ngm), STAT = ierr)
  IF (ierr /= 0) CALL errore('newdq2', 'Error allocating aux1', 1)
  ALLOCATE(aux2(ngm, nspin_mag), STAT = ierr)
  IF (ierr /= 0) CALL errore('newdq2', 'Error allocating aux2', 1)
  ALLOCATE(veff(dfftp%nnr), STAT = ierr)
  IF (ierr /= 0) CALL errore('newdq2', 'Error allocating veff', 1)
  ALLOCATE(ylmk0(ngm, lmaxq * lmaxq), STAT = ierr)
  IF (ierr /= 0) CALL errore('newdq2', 'Error allocating ylmk0', 1)
  ALLOCATE(qgm(ngm), STAT = ierr)
  IF (ierr /= 0) CALL errore('newdq2', 'Error allocating qgm', 1)
  ALLOCATE(qmod(ngm), STAT = ierr)
  IF (ierr /= 0) CALL errore('newdq2', 'Error allocating qmod', 1)
  ALLOCATE(qg(3, ngm), STAT = ierr)
  IF (ierr /= 0) CALL errore('newdq2', 'Error allocating qg', 1)
  aux1(:)     = czero
  aux2(:, :)  = czero
  veff(:)     = czero
  ylmk0(:, :) = zero
  qgm(:)      = czero
  qmod(:)     = zero
  qg(:, :)    = zero
  !
  ! first compute the spherical harmonics
  !
  CALL setqmod(ngm, xq0, g, qmod, qg)
  CALL ylmr2(lmaxq * lmaxq, ngm, qg, qmod, ylmk0)
  !
  DO ig = 1, ngm
    qmod(ig) = SQRT(qmod(ig))
  ENDDO
  !
  ! and for each perturbation of this irreducible representation
  ! integrate the change of the self consistent potential and
  ! the Q functions
  !
  DO ipert = 1, npe
    DO is = 1, nspin_mag
      DO ir = 1, dfftp%nnr
        IF (timerev) THEN
          veff(ir) = CONJG(dvscf(ir, is, ipert))
        ELSE
          veff(ir) = dvscf(ir, is, ipert)
        ENDIF
      ENDDO
      CALL fwfft('Rho', veff, dfftp)
      DO ig = 1, ngm
        aux2(ig, is) = veff(dfftp%nl(ig))
      ENDDO
    ENDDO
    !
    DO nt = 1, ntyp
      IF (upf(nt)%tvanp) THEN
        !
        DO ih = 1, nh(nt)
          DO jh = ih, nh(nt)
            !
            CALL qvan2(ngm, ih, jh, nt, qmod, qgm, ylmk0)
            !
            DO na = 1, nat
              IF (ityp(na) == nt) THEN
                DO ig = 1, ngm
                  aux1(ig) = qgm(ig) * eigts1(mill(1, ig), na) * &
                                       eigts2(mill(2, ig), na) * &
                                       eigts3(mill(3, ig), na) * &
                                       eigqts(na)
                ENDDO
                DO is = 1, nspin_mag
                  int3(ih, jh, na, is, ipert) = omega * ZDOTC(ngm, aux1, 1, aux2(1, is), 1)
                ENDDO
              ENDIF
            ENDDO
          ENDDO ! jh
        ENDDO ! ih
        ! 
        DO na = 1, nat
          IF (ityp(na) == nt) THEN
            !
            ! We use the symmetry properties of the ps factor
            !
            DO ih = 1, nh(nt)
              DO jh = ih, nh(nt)
                DO is = 1, nspin_mag
                  int3(jh, ih, na, is, ipert) = int3(ih, jh, na, is, ipert)
                ENDDO
              ENDDO
            ENDDO
          ENDIF ! ityp
        ENDDO ! na
      ENDIF ! upf
    ENDDO ! nt
  ENDDO ! ipert
  !
  CALL mp_sum(int3, intra_pool_comm)
  !
  IF (noncolin) CALL set_int3_nc(npe)
  !
  DEALLOCATE(aux1, STAT = ierr)
  IF (ierr /= 0) CALL errore('newdq2', 'Error deallocating aux1', 1)
  DEALLOCATE(aux2, STAT = ierr)
  IF (ierr /= 0) CALL errore('newdq2', 'Error deallocating aux2', 1)
  DEALLOCATE(veff, STAT = ierr)
  IF (ierr /= 0) CALL errore('newdq2', 'Error deallocating veff', 1)
  DEALLOCATE(ylmk0, STAT = ierr)
  IF (ierr /= 0) CALL errore('newdq2', 'Error deallocating ylmk0', 1)
  DEALLOCATE(qgm, STAT = ierr)
  IF (ierr /= 0) CALL errore('newdq2', 'Error deallocating qgm', 1)
  DEALLOCATE(qmod, STAT = ierr)
  IF (ierr /= 0) CALL errore('newdq2', 'Error deallocating qmod', 1)
  DEALLOCATE(qg, STAT = ierr)
  IF (ierr /= 0) CALL errore('newdq2', 'Error deallocating qg', 1)
  !
  CALL stop_clock('newdq2')
  !
  RETURN
  !
  !----------------------------------------------------------------------
  END SUBROUTINE newdq2
  !----------------------------------------------------------------------
