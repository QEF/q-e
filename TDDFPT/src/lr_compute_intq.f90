!
! Copyright (C) 2001-2019 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------
SUBROUTINE lr_compute_intq
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
!  USE lr_variables,         ONLY : intq
  USE qpoint,               ONLY : xq, eigqts
  USE lrus,                 ONLY : intq

  IMPLICIT NONE

  INTEGER :: na, ig, nt, ir, ih, jh
  ! counters
  COMPLEX(DP) :: qq_nt(nhm,nhm,ntyp)
  ! work space

  IF (.NOT.okvan) RETURN
  CALL start_clock ('lr_compute_intq')
  !
  CALL compute_qqc ( tpiba, xq, omega, qq_nt )
  !
  intq (:,:,:) = (0.D0, 0.0D0)
  DO nt = 1, ntyp
     IF (upf(nt)%tvanp ) THEN
        DO ih = 1, nh (nt)
           DO jh = ih, nh (nt)
              DO na = 1, nat
                 IF (ityp (na) == nt) THEN
                    intq(ih,jh,na) = CONJG( qq_nt(ih,jh,nt) * eigqts(na) )
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

  IF (noncolin) CALL lr_set_intq_nc()

  CALL stop_clock ('lr_compute_intq')
  RETURN
END SUBROUTINE lr_compute_intq
