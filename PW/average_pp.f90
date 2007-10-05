!
! Copyright (C) 2005-2006 Quantum-ESPRESSO group
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
  USE atom,             ONLY : chi, nchi, lchi, jchi, rgrid
  USE spin_orb,         ONLY : so
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
     IF ( so(nt) ) THEN
        !
        IF ( upf(nt)%tvanp ) &
             CALL errore( 'setup', 'US j-average not yet implemented', 1 )
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
                      > 1.d-7 ) call errore('setup','wrong beta functions',1)
                 ind=nbe+1
                 ind1=nbe
              ELSE
                 IF (ABS(upf(nt)%jjj(nbe+1)-upf(nt)%lll(nbe+1)+0.5d0) > 1.d-7) &
                      call errore('setup','wrong beta functions',1)
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
        DO nb = 1, nchi(nt)
           !
           nbe = nbe + 1
           !
           IF ( lchi(nb,nt) /= 0 .AND. &
                ABS(jchi(nb,nt)-lchi(nb,nt)-0.5D0 ) < 1.D-7 ) &
              nbe = nbe - 1
           !
        END DO
        !
        nchi(nt) = nbe
        ! 
        nbe = 0
        !
        do nb = 1, nchi(nt)
           !
           nbe = nbe + 1
           !
           l = lchi(nbe,nt)
           !
           IF ( l /= 0 ) THEN
              !
              IF (ABS(jchi(nbe,nt)-lchi(nbe,nt)+0.5d0) < 1.d-7) THEN
                 IF ( ABS(jchi(nbe+1,nt)-lchi(nbe+1,nt)-0.5d0) > &
                      1.d-7) call errore('setup','wrong chi functions',1)
                 ind=nbe+1
                 ind1=nbe
              ELSE
                 IF ( ABS(jchi(nbe+1,nt)-lchi(nbe+1,nt)+0.5d0) > &
                      1.d-7) call errore('setup','wrong chi functions',1)
                 ind=nbe
                 ind1=nbe+1
              END IF
              !
              chi(1:rgrid(nt)%mesh,nb,nt) = &
                 ((l+1.D0) * chi(1:rgrid(nt)%mesh,ind,nt)+ &
                   l * chi(1:rgrid(nt)%mesh,ind1,nt)) / ( 2.D0 * l + 1.D0 )
              !   
              nbe = nbe + 1
              !
           ELSE
              !
              chi(1:rgrid(nt)%mesh,nb,nt) = chi(1:rgrid(nt)%mesh,nbe,nt)
              !
           END IF
           !
           lchi(nb,nt)= lchi(nbe,nt)
           !
        END DO
        !
     END IF
     !
     so(nt) = .FALSE.
     !
  END DO
  !
END SUBROUTINE average_pp
