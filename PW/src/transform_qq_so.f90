
! Copyright (C) 2012 Quantum-Espresso group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------
SUBROUTINE transform_qq_so(qq,qq_so)
  !----------------------------------------------------------------------
  !
  !
  USE kinds,        ONLY : DP
  USE ions_base,    ONLY : ntyp => nsp
  USE uspp_param,   ONLY : upf, nhm, nh
  USE spin_orb,     ONLY : lspinorb, fcoef
  !
  implicit none
  !
  !     here a few local variables
  !
  integer :: nt, ih, jh, kh, lh, ijs, is1, is2, is
  complex(DP) :: qq(nhm,nhm,ntyp), qq_so(nhm,nhm,4,ntyp)

  qq_so=(0.0_DP, 0.0_DP)
  DO nt = 1, ntyp
    IF ( upf(nt)%tvanp ) THEN
      IF (upf(nt)%has_so) THEN
        DO ih=1,nh(nt)
          DO jh=1,nh(nt)
            DO kh=1,nh(nt)
              DO lh=1,nh(nt)
                ijs=0
                DO is1=1,2
                  DO is2=1,2
                    ijs=ijs+1
                    DO is=1,2
                      qq_so(kh,lh,ijs,nt) = qq_so(kh,lh,ijs,nt)       &
                          + qq(ih,jh,nt)*fcoef(kh,ih,is1,is,nt)&
                                               *fcoef(jh,lh,is,is2,nt)
                    ENDDO
                  ENDDO
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ELSE
        DO ih = 1, nh (nt)
          DO jh = ih, nh (nt)
             IF (lspinorb) THEN
                 qq_so (ih, jh, 1, nt) = qq (ih, jh, nt) 
                 qq_so (jh, ih, 1, nt) = qq_so (ih, jh, 1, nt)
                 qq_so (ih, jh, 4, nt) = qq_so (ih, jh, 1, nt)
                 qq_so (jh, ih, 4, nt) = qq_so (ih, jh, 4, nt)
             ENDIF
          ENDDO
        ENDDO
      ENDIF
    ENDIF
  ENDDO

  RETURN
END SUBROUTINE transform_qq_so
