!
! Copyright (C) 2001-2008 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------
SUBROUTINE compute_intq
  !----------------------------------------------------------------------
  !
  !     This routine computes the contribution of the selfconsistent
  !     change of the potential to the known part of the linear
  !     system and adds it to dvpsi.
  !
  USE kinds,                ONLY : DP
  USE ions_base,            ONLY : nat, ityp, ntyp => nsp
  USE noncollin_module,     ONLY : noncolin
  USE cell_base,            ONLY : omega, tpiba
  USE uspp,                 ONLY : okvan
  USE uspp_param,           ONLY : upf, lmaxq, nh, nhm

  USE lrus,                 ONLY : intq
  USE qpoint,               ONLY : xq, eigqts

  IMPLICIT NONE

  INTEGER :: na, ig, nt, ir, ih, jh
  ! countera

  REAL(DP), ALLOCATABLE ::  ylmk0 (:,:)
  ! the modulus of q+G
  ! the values of q+G
  ! the spherical harmonics


  ! work space
  COMPLEX(DP) :: qgm(1), aux1
  REAL(DP)    :: qmod(1), zero(3,1), qg(3,1)

  IF (.NOT.okvan) RETURN
  CALL start_clock ('compute_intq')

  intq (:,:,:) = (0.D0, 0.0D0)
  ALLOCATE (ylmk0(1 , lmaxq * lmaxq))
  !
  !    first compute the spherical harmonics
  !
  zero=0.0_DP
  CALL setqmod (1, xq, zero, qmod, qg)
  CALL ylmr2 (lmaxq * lmaxq, 1, qg, qmod, ylmk0)
  qmod(1) = SQRT (qmod(1))*tpiba

  DO nt = 1, ntyp
     IF (upf(nt)%tvanp ) THEN
        DO ih = 1, nh (nt)
           DO jh = ih, nh (nt)
              CALL qvan2 (1, ih, jh, nt, qmod, qgm, ylmk0)
              DO na = 1, nat
                 IF (ityp (na) == nt) THEN
                    aux1 = qgm(1) * eigqts(na)
                    intq(ih,jh,na) = omega * CONJG(aux1)
                 ENDIF
              ENDDO
           ENDDO
        ENDDO
        DO na = 1, nat
           IF (ityp(na) == nt) THEN
              !
              !    We use the symmetry properties of the ps factor
              !
              DO ih = 1, nh (nt)
                 DO jh = ih, nh (nt)
                    intq(jh,ih,na) = intq(ih,jh,na)
                 ENDDO
              ENDDO
           ENDIF
        ENDDO
     ENDIF
  ENDDO

  IF (noncolin) CALL set_intq_nc()

  DEALLOCATE (ylmk0)

  CALL stop_clock ('compute_intq')
  RETURN
END SUBROUTINE compute_intq
