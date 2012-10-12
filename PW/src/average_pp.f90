!
! Copyright (C) 2005-2006 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
SUBROUTINE average_pp ( ntyp ) 
  !----------------------------------------------------------------------------
  !
  USE kinds,            ONLY : DP
  USE atom,             ONLY : rgrid
  USE uspp_param,       ONLY : upf
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: ntyp
  !
  INTEGER  :: nt, nb, nbe, ind, ind1, l
  REAL(DP) :: vionl
  !
  !
  DO nt = 1, ntyp
     !
     IF ( upf(nt)%has_so ) THEN
        !
        IF ( upf(nt)%tvanp ) &
             CALL errore( 'average_pp', 'FR-PP please use lspinorb=.true.', 1 )
        !
        nbe = 0
        !
        DO nb = 1, upf(nt)%nbeta
           !
           nbe = nbe + 1
           !
           IF ( upf(nt)%lll(nb) /= 0 .AND. &
                ABS( upf(nt)%jjj(nb) - upf(nt)%lll(nb) - 0.5D0 ) < 1.D-7 ) &
              nbe = nbe - 1
        END DO
        !
        upf(nt)%nbeta = nbe
        !
        nbe = 0
        !
        DO nb = 1, upf(nt)%nbeta
           !
           nbe = nbe + 1
           !
           l = upf(nt)%lll(nbe)
           !
           IF ( l /= 0 ) THEN
              !
              IF (ABS(upf(nt)%jjj(nbe)-upf(nt)%lll(nbe)+0.5d0) < 1.d-7) THEN
                 IF ( ABS( upf(nt)%jjj(nbe+1)-upf(nt)%lll(nbe+1)-0.5d0 ) &
                      > 1.d-7 ) call errore('average_pp','wrong beta functions',1)
                 ind=nbe+1
                 ind1=nbe
              ELSE
                 IF (ABS(upf(nt)%jjj(nbe+1)-upf(nt)%lll(nbe+1)+0.5d0) > 1.d-7) &
                      call errore('average_pp','wrong beta functions',2)
                 ind=nbe
                 ind1=nbe+1
              ENDIF
              !
              vionl = ( ( l + 1.D0 ) * upf(nt)%dion(ind,ind) + &
                   l * upf(nt)%dion(ind1,ind1) ) / ( 2.D0 * l + 1.D0 )
              !
              upf(nt)%beta(1:rgrid(nt)%mesh,nb) = 1.D0 / ( 2.D0 * l + 1.D0 ) * &
                   ( ( l + 1.D0 ) * SQRT( upf(nt)%dion(ind,ind) / vionl ) * &
                   upf(nt)%beta(1:rgrid(nt)%mesh,ind) + &
                   l * SQRT( upf(nt)%dion(ind1,ind1) / vionl ) * &
                   upf(nt)%beta(1:rgrid(nt)%mesh,ind1) )
              !
              upf(nt)%dion(nb,nb) = vionl
              !
              nbe = nbe + 1
                 !
           ELSE
              !
              upf(nt)%beta(1:rgrid(nt)%mesh,nb) = &
                  upf(nt)%beta(1:rgrid(nt)%mesh,nbe)
              !
              upf(nt)%dion(nb,nb) = upf(nt)%dion(nbe,nbe)
              !
           END IF
           !
           upf(nt)%lll(nb)=upf(nt)%lll(nbe)
           !
        END DO
        !
        nbe = 0
        !
        DO nb = 1, upf(nt)%nwfc
           !
           nbe = nbe + 1
           !
           IF ( upf(nt)%lchi(nb) /= 0 .AND. &
                ABS(upf(nt)%jchi(nb)-upf(nt)%lchi(nb)-0.5D0 ) < 1.D-7 ) &
              nbe = nbe - 1
           !
        END DO
        !
        upf(nt)%nwfc = nbe
        ! 
        nbe = 0
        !
        do nb = 1, upf(nt)%nwfc
           !
           nbe = nbe + 1
           !
           l = upf(nt)%lchi(nbe)
           !
           IF ( l /= 0 ) THEN
              !
              IF (ABS(upf(nt)%jchi(nbe)-upf(nt)%lchi(nbe)+0.5d0) < 1.d-7) THEN
                 IF ( ABS(upf(nt)%jchi(nbe+1)-upf(nt)%lchi(nbe+1)-0.5d0) > &
                      1.d-7) call errore('average_pp','wrong chi functions',3)
                 ind=nbe+1
                 ind1=nbe
              ELSE
                 IF ( ABS(upf(nt)%jchi(nbe+1)-upf(nt)%lchi(nbe+1)+0.5d0) > &
                      1.d-7) call errore('average_pp','wrong chi functions',4)
                 ind=nbe
                 ind1=nbe+1
              END IF
              !
              upf(nt)%chi(1:rgrid(nt)%mesh,nb) = &
                 ((l+1.D0) * upf(nt)%chi(1:rgrid(nt)%mesh,ind)+ &
                   l * upf(nt)%chi(1:rgrid(nt)%mesh,ind1)) / ( 2.D0 * l + 1.D0 )
              !   
              nbe = nbe + 1
              !
           ELSE
              !
              upf(nt)%chi(1:rgrid(nt)%mesh,nb) = upf(nt)%chi(1:rgrid(nt)%mesh,nbe)
              !
           END IF
           !
           upf(nt)%lchi(nb)= upf(nt)%lchi(nbe)
           !
        END DO
        !
     END IF
     !
     upf(nt)%has_so = .FALSE.
     !
  END DO
  !
END SUBROUTINE average_pp
