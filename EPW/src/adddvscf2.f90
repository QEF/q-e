  !
  ! Copyright (C) 2010-2016 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino
  !
  ! Copyright (C) 2001-2003 PWSCF group
  ! This file is distributed under the terms of the
  ! GNU General Public License. See the file `License'
  ! in the root directory of the present distribution,
  ! or http://www.gnu.org/copyleft/gpl.txt .
  !
  ! adapted from LR_Modules/adddvscf.f90 (QE)
  !
  !----------------------------------------------------------------------
  SUBROUTINE adddvscf2( ipert, ik )
  !----------------------------------------------------------------------
  !!
  !! This routine computes the contribution of the selfconsistent
  !! change of the potential to the known part of the linear
  !! system and adds it to dvpsi.
  !! It implements the second term in Eq. B30 of PRB 64, 235118 (2001).
  !! Only used in the case of USPP.
  !!
  !! Roxana Margine - Jan 2019: Updated based on QE 6.3
  !! SP - Jan 2019: Clean
  !!
  USE kinds,      ONLY : DP
  USE uspp_param, ONLY : upf, nh
  USE uspp,       ONLY : vkb, okvan
  USE lsda_mod,   ONLY : lsda, current_spin, isk
  USE ions_base,  ONLY : ntyp => nsp, nat, ityp
  USE wvfct,      ONLY : npwx
  USE lrus,       ONLY : int3, int3_nc, becp1
  USE qpoint,     ONLY : npwq
  USE eqv,        ONLY : dvpsi
  USE elph2,      ONLY : lower_band, upper_band
  USE noncollin_module, ONLY : noncolin, npol
  USE constants_epw,    ONLY : czero
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(in) :: ik
  !! Counter on k-point
  INTEGER, INTENT(in) :: ipert
  !! Counter on Vscf perturbations
  !
  !   Local variables
  !
  INTEGER :: na
  !! Counter on atoms
  INTEGER :: nt
  !! Counter on atomic types
  INTEGER :: ibnd
  !! Counter on bands
  INTEGER :: ih
  !! Counter on beta functions
  INTEGER :: jh
  !! Counter on beta functions
  INTEGER :: ijkb0
  !! Auxiliary variable for counting
  INTEGER :: ikb
  !! Counter on becp functions
  INTEGER :: jkb
  !! Counter on becp functions
  INTEGER :: is
  !! Counter on polarization
  INTEGER :: js
  !! Counter on polarization
  INTEGER ::  ijs
  !! Counter on combined is and js polarization
  !
  COMPLEX(kind=DP) :: sum_k
  !! auxiliary sum variable
  COMPLEX(kind=DP) :: sum_nc(npol)
  !! auxiliary sum variable non-collinear case
  !
  IF (.not.okvan) RETURN
  !
  CALL start_clock('adddvscf2')
  !
  IF (lsda) current_spin = isk(ik)
  !
  ijkb0 = 0
  DO nt = 1, ntyp
    IF ( upf(nt)%tvanp ) THEN
      DO na = 1, nat
        IF (ityp(na) .eq. nt) THEN
          !
          !   we multiply the integral for the becp term and the beta_n
          !
          DO ibnd = lower_band, upper_band
            DO ih = 1, nh(nt)
               ikb = ijkb0 + ih
               IF (noncolin) THEN
                 sum_nc = czero
               ELSE
                 sum_k = czero
               ENDIF
               DO jh = 1, nh(nt)
                 jkb = ijkb0 + jh
                 IF (noncolin) THEN
                   ijs = 0
                   DO is = 1, npol
                     DO js = 1, npol
                       ijs = ijs + 1
                       sum_nc(is) = sum_nc(is) + int3_nc(ih,jh,na,ijs,ipert) * &
                             becp1(ik)%nc(jkb,js,ibnd)
                     ENDDO
                   ENDDO
                 ELSE
                   sum_k = sum_k + int3(ih,jh,na,current_spin,ipert) * &
                               becp1(ik)%k(jkb,ibnd)
                 ENDIF
               ENDDO
               IF (noncolin) THEN
                  CALL zaxpy( npwq, sum_nc(1), vkb(1,ikb), 1, dvpsi(1,ibnd), 1 )
                  CALL zaxpy( npwq, sum_nc(2), vkb(1,ikb), 1, dvpsi(1+npwx,ibnd), 1 )
               ELSE
                  CALL zaxpy( npwq, sum_k, vkb(1,ikb), 1, dvpsi(1,ibnd), 1 )
               ENDIF
            ENDDO
          ENDDO
          ijkb0 = ijkb0 + nh(nt)
        ENDIF
      ENDDO
    ELSE
       DO na = 1, nat
          IF (ityp(na) .eq. nt) ijkb0 = ijkb0 + nh(nt)
       ENDDO
    ENDIF
  ENDDO
  !
  CALL stop_clock('adddvscf2')
  !
  RETURN
  !
  END SUBROUTINE adddvscf2
