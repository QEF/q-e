!
! Copyright (C) 2001-2015 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------------
SUBROUTINE lr_addusddens (npert, dbecsum, drhop)
  !---------------------------------------------------------------------------
  !
  ! Self-consistent ultrasoft-augmentation contribution to drho.
  ! Given dbecsum (built from the Sternheimer dpsi), adds to drhop the
  ! corresponding augmentation density: Eq. (36) of B. Walker and
  ! R. Gebauer, J. Chem. Phys. 127, 164106 (2007).
  !
  ! This is the SCF augmentation piece that applies to every DFPT
  ! perturbation type (phonon, electric field, EELS, Hubbard, ...). For
  ! phonons there is an additional non-self-consistent Pulay-like
  ! contribution coming from the displacement of the projector centers;
  ! that piece is computed by PHonon/PH/addusddens_pulay.f90 and stored in
  ! dfpt_data%drhop_pulay. It is added in response_kernels.f90, not here.
  !
  ! It assumes that the array dbecsum has already been computed.
  !
  ! Created by Iurii Timrov (2013)
  !
  USE kinds,                ONLY : DP
  USE ions_base,            ONLY : nat, ityp, ntyp => nsp
  USE cell_base,            ONLY : tpiba
  USE fft_base,             ONLY : dfftp
  USE fft_interfaces,       ONLY : invfft
  USE gvect,                ONLY : ngm, g, eigts1, eigts2, eigts3, mill
  USE noncollin_module,     ONLY : nspin_mag
  USE uspp,                 ONLY : okvan
  USE uspp_param,           ONLY : upf, lmaxq, nh, nhm
  USE qpoint,               ONLY : xq, eigqts
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(in) :: npert
  ! input : number of perturbations
  COMPLEX(DP), INTENT(in)    :: dbecsum(nhm*(nhm+1)/2, nat, nspin_mag, npert)
  ! input : the ultrasoft term
  COMPLEX(DP), INTENT(inout) :: drhop(dfftp%nnr, nspin_mag, npert)
  ! input/output : change of the charge density
  !
  ! the local variables
  !
  COMPLEX(DP) :: dbec
  ! temporary scalar to hold dbecsum value
  !
  INTEGER :: ig, na, nt, ih, jh, is, ijh, ir, nab, nb
  ! counter on G vectors
  ! counter on atoms
  ! counter on atomic type
  ! counter on beta functions
  ! counter on beta functions
  ! counter on spin
  ! counter on combined beta functions
  ! counter on r vectors
  ! max number of atoms of type nt
  ! counter on atoms of the same type
  !
  INTEGER :: ipert
  !! counter on perturbations
  !
  REAL(DP), ALLOCATABLE :: qmod(:), qpg(:,:), ylmk0(:,:)
  ! the modulus of q+G
  ! the values of q+G
  ! the spherical harmonics
  !
  COMPLEX(DP), ALLOCATABLE :: sk(:,:), qgm(:), aux(:, :, :), aux_r(:)
  ! the structure factor
  ! q_lm(G)
  ! auxiliary variable for drho(G)
  ! auxiliary variable for drho(r)
  !
  IF (.NOT.okvan) RETURN
  !
  CALL start_clock ('lr_addusddens')
  !
  ALLOCATE (aux(ngm, nspin_mag, npert))
  ALLOCATE (aux_r(dfftp%nnr))
  ALLOCATE (ylmk0(ngm,lmaxq * lmaxq))
  ALLOCATE (qgm(ngm))
  ALLOCATE (qmod(ngm))
  ALLOCATE (qpg(3,ngm))
  !
  ! Calculate the q+G vector, its modulus, and the spherical harmonics.
  !
  CALL setqmod (ngm, xq, g, qmod, qpg)
  !
  !$acc data create(ylmk0, qgm) copyin(qmod, eigqts) copyout(aux)
  !
  ! ylmr2 has copyin(qpg, qmod) and copyout(ylmk0)
  CALL ylmr2 (lmaxq * lmaxq, ngm, qpg, qmod, ylmk0)
  DEALLOCATE (qpg)
  !
  !$acc parallel loop
  DO ig = 1, ngm
     qmod(ig) = sqrt(qmod(ig)) * tpiba
  ENDDO
  !
  !$acc kernels
  aux(:, :, :) = (0.d0, 0.d0)
  !$acc end kernels
  !
  DO nt = 1, ntyp
     IF (upf(nt)%tvanp) THEN
        !
        ! count the number of atoms of type nt and allocate sk accordingly
        !
        nab = 0
        DO na = 1, nat
           IF ( ityp(na) == nt ) nab = nab + 1
        ENDDO
        !
        ALLOCATE( sk(ngm,nab) )
        !$acc data create(sk)
        !
        nb = 0
        DO na = 1, nat
           IF ( ityp(na) == nt ) THEN
              nb = nb + 1
              !
              ! Calculate the structure factor for all atoms of type nt
              !
              !$acc parallel loop present(eigts1,eigts2,eigts3,mill,eigqts)
              DO ig = 1, ngm
                 sk(ig,nb) = eigts1(mill(1,ig),na) * &
                             eigts2(mill(2,ig),na) * &
                             eigts3(mill(3,ig),na) * &
                             eigqts(na)
              ENDDO
           ENDIF
        ENDDO
        !
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
              nb = 0
              DO na = 1, nat
                 IF (ityp (na) .eq.nt) THEN
                    nb = nb + 1
                    !
                    ! Calculate the second term in Eq.(36) of the ultrasoft paper.
                    !
                    DO ipert = 1, npert
                       DO is = 1, nspin_mag
                          dbec = dbecsum(ijh, na, is, ipert)
                          !$acc parallel loop present(aux, qgm, sk)
                          DO ig = 1, ngm
                             !
                             aux(ig, is, ipert) = aux(ig, is, ipert) &
                                + 2.0d0 * qgm(ig) * sk(ig, nb) * dbec
                             !
                          ENDDO
                       ENDDO
                    ENDDO ! ipert
                    !
                 ENDIF
              ENDDO
           ENDDO
        ENDDO
        !
        !$acc end data ! affects only the device-resident sk array
        DEALLOCATE( sk )
        !
     ENDIF
  ENDDO
  !$acc end data
  ! aux is copied back to host due to copyout(aux)
  !
  ! Convert aux to real space, and add to the charge density.
  !
  DO ipert = 1, npert
     DO is = 1, nspin_mag
         !
         aux_r(:) = (0.d0, 0.d0)
         !
         DO ig = 1, ngm
            aux_r(dfftp%nl(ig)) = aux(ig, is, ipert)
         ENDDO
         !
         CALL invfft('Rho', aux_r, dfftp)
         !
         DO ir = 1, dfftp%nnr
            drhop(ir, is, ipert) = drhop(ir, is, ipert) + aux_r(ir)
         ENDDO
         !
     ENDDO
  ENDDO ! ipert
  !
  DEALLOCATE (qmod)
  DEALLOCATE (qgm)
  DEALLOCATE (ylmk0)
  DEALLOCATE (aux)
  DEALLOCATE (aux_r)
  !
  CALL stop_clock ('lr_addusddens')
  !
  RETURN
  !
END SUBROUTINE lr_addusddens
