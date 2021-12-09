
!
! Copyright (C) 2006 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------
SUBROUTINE compute_qdipol_so( dpqq, dpqq_so )
  !----------------------------------------------------------------------
  !! This routine multiplies the dpqq coefficients with spin-orbit 
  !! "fcoef" coefficients - to be called in spin-orbit case only.
  !
  USE kinds,      ONLY : DP
  USE ions_base,  ONLY : ntyp => nsp
  USE lsda_mod,   ONLY : nspin
  USE uspp_param, ONLY : upf, nh, nhm
  USE upf_spinorb,ONLY : fcoef
  !
  IMPLICIT NONE
  !
  REAL(DP) :: dpqq(nhm,nhm,3,ntyp)
  !! dpqq coefficients
  COMPLEX(DP) :: dpqq_so(nhm,nhm,nspin,3,ntyp)
  !! spin orbit fcoef coefficients
  !
  !  ... local variables
  !
  INTEGER :: ipol
  INTEGER :: nt, ih, jh, kh, lh, ijs, is1, is2, is
  !
  dpqq_so = (0.d0,0.d0)
  !
  DO ipol = 1, 3
    DO nt = 1, ntyp
      !
      IF ( upf(nt)%tvanp ) THEN
        IF (upf(nt)%has_so) THEN
          DO ih = 1, nh(nt)
            DO jh = 1, nh(nt)
              DO kh = 1, nh(nt)
                DO lh = 1, nh(nt)
                  ijs = 0
                  !
                  DO is1 = 1, 2
                    DO is2 = 1, 2
                      ijs = ijs+1
                      DO is = 1, 2
                        dpqq_so(kh,lh,ijs,ipol,nt)=dpqq_so(kh,lh,ijs,ipol,nt)&
                              +dpqq(ih,jh,ipol,nt)*fcoef(kh,ih,is1,is,nt)    &
                                                  *fcoef(jh,lh,is,is2,nt)
                      ENDDO
                    ENDDO
                  ENDDO
                  !
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ELSE
          DO ih = 1, nh (nt)
             DO jh = ih, nh (nt)
                dpqq_so(ih,jh,1,ipol,nt) = dpqq(ih,jh,ipol,nt)
                dpqq_so(jh,ih,1,ipol,nt) = dpqq_so(ih,jh,1,ipol,nt)
                dpqq_so(ih,jh,4,ipol,nt) = dpqq_so(ih,jh,1,ipol,nt)
                dpqq_so(jh,ih,4,ipol,nt) = dpqq_so(ih,jh,4,ipol,nt)
             ENDDO
          ENDDO
        ENDIF
      ENDIF
      !
    ENDDO
  ENDDO
  !
  RETURN
  !
END SUBROUTINE compute_qdipol_so
