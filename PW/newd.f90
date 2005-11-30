!
! Copyright (C) 2001-2004 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "f_defs.h"
!
SUBROUTINE newd()
  use realus, ONLY: tqr,  newdrealsub
  if (tqr) then
     call newdrealsub()
  else
     call newd_g()
  end if
  return
END SUBROUTINE newd
!----------------------------------------------------------------------------
SUBROUTINE newd_g()
  !----------------------------------------------------------------------------
  !
  ! ... This routine computes the integral of the effective potential with
  ! ... the Q function and adds it to the bare ionic D term which is used
  ! ... to compute the non-local term in the US scheme.
  !
  USE kinds,                ONLY : DP
  USE ions_base,            ONLY : nat, ntyp => nsp, ityp
  USE cell_base,            ONLY : omega
  USE gvect,                ONLY : nr1, nr2, nr3, nrx1, nrx2, nrx3, &
                                   g, gg, ngm, gstart, ig1, ig2, ig3, &
                                   eigts1, eigts2, eigts3, nl
  USE lsda_mod,             ONLY : nspin
  USE scf,                  ONLY : vr, vltot
  USE uspp,                 ONLY : deeq, dvan, deeq_nc, dvan_so, okvan
  USE uspp_param,           ONLY : lmaxq, nh, nhm, tvanp
  USE wvfct,                ONLY : gamma_only
  USE wavefunctions_module, ONLY : psic
  USE spin_orb,             ONLY : lspinorb
  USE noncollin_module,     ONLY : noncolin
  !
  IMPLICIT NONE
  !
  INTEGER :: ig, nt, ih, jh, na, is, nht
    ! counters on g vectors, atom type, beta functions x 2, atoms, spin
  COMPLEX(DP), ALLOCATABLE :: aux(:,:), qgm(:), qgm_na(:)
    ! work space
  REAL(DP), ALLOCATABLE :: ylmk0(:,:), qmod(:)
    ! spherical harmonics, modulus of G
  REAL(DP) :: fact, DDOT
  !
  !
  IF ( .NOT. okvan ) THEN
     !
     ! ... no ultrasoft potentials: use bare coefficients for projectors
     !
     DO na = 1, nat
        !
        nt  = ityp(na)
        nht = nh(nt)
        !
        IF ( lspinorb ) THEN
           !
           deeq_nc(1:nht,1:nht,na,1:nspin) = dvan_so(1:nht,1:nht,1:nspin,nt)
           !
        ELSE IF ( noncolin ) THEN
           !
           deeq_nc(1:nht,1:nht,na,1) = dvan(1:nht,1:nht,nt)
           deeq_nc(1:nht,1:nht,na,2) = ( 0.D0, 0.D0 )
           deeq_nc(1:nht,1:nht,na,3) = ( 0.D0, 0.D0 )
           deeq_nc(1:nht,1:nht,na,4) = dvan(1:nht,1:nht,nt)
           !
        ELSE
           !
           DO is = 1, nspin
              !
              deeq(1:nht,1:nht,na,is) = dvan(1:nht,1:nht,nt)
              !
           END DO
           !
        END IF
        !
     END DO
     !
     ! ... early return
     !
     RETURN
     !
  END IF
  !
  IF ( gamma_only ) THEN
     !
     fact = 2.D0
     !
  ELSE
     !
     fact = 1.D0
     !
  END IF
  !
  CALL start_clock( 'newd' )
  !
  ALLOCATE( aux( ngm, nspin ), qgm_na( ngm ), &
            qgm( ngm ), qmod( ngm ), ylmk0( ngm, lmaxq*lmaxq ) )
  !
  deeq(:,:,:,:) = 0.D0
  !
  CALL ylmr2( lmaxq * lmaxq, ngm, g, gg, ylmk0 )
  !
  qmod(1:ngm) = SQRT( gg(1:ngm) )
  !
  ! ... fourier transform of the total effective potential
  !
  DO is = 1, nspin
     !
     IF ( nspin == 4 .AND. is /= 1 ) THEN 
        !
        psic(:) = vr(:,is)
        !
     ELSE
        !
        psic(:) = vltot(:) + vr(:,is)
        !
     END IF
     !
     CALL cft3( psic, nr1, nr2, nr3, nrx1, nrx2, nrx3, - 1 )
     !
     aux(1:ngm,is) = psic( nl(1:ngm) )
     !
  END DO
  !
  ! ... here we compute the integral Q*V for each atom,
  ! ...       I = sum_G exp(-iR.G) Q_nm v^*
  !
  DO nt = 1, ntyp
     !
     IF ( tvanp(nt) ) THEN
        !
        DO ih = 1, nh(nt)
           !
           DO jh = ih, nh(nt)
              !
              ! ... The Q(r) for this atomic species without structure factor
              !
              CALL qvan2( ngm, ih, jh, nt, qmod, qgm, ylmk0 )
              !
              DO na = 1, nat
                 !
                 IF ( ityp(na) == nt ) THEN
                    !
                    ! ... The Q(r) for this specific atom
                    !
                    qgm_na(1:ngm) = qgm(1:ngm) * eigts1(ig1(1:ngm),na) &
                                               * eigts2(ig2(1:ngm),na) &
                                               * eigts3(ig3(1:ngm),na)
                    !
                    ! ... and the product with the Q functions
                    !
                    DO is = 1, nspin
                       !
                       deeq(ih,jh,na,is) = fact * omega * &
                                        DDOT( 2 * ngm, aux(1,is), 1, qgm_na, 1 )
                       !
                       IF ( gamma_only .AND. gstart == 2 ) &
                           deeq(ih,jh,na,is) = deeq(ih,jh,na,is) - &
                                           omega * DBLE( aux(1,is) * qgm_na(1) )
                       !
                       deeq(jh,ih,na,is) = deeq(ih,jh,na,is)
                       !
                    END DO
                    !
                 END IF
                 !
              END DO
              !
           END DO
           !
        END DO
        !
     END IF
     !
  END DO
  !
  CALL reduce( nhm * nhm * nat * nspin, deeq )
  !
  IF ( lspinorb ) THEN
     !
     CALL newd_so()
     !
  ELSE IF ( noncolin ) THEN 
     !
     CALL newd_nc()
     !
  ELSE
     !
     DO na = 1, nat
        !
        nt = ityp(na)
        !
        DO is = 1, nspin
           !
           DO ih = 1, nh(nt)
              !
              DO jh = ih, nh(nt)
                 !
                 deeq(ih,jh,na,is) = deeq(ih,jh,na,is) + dvan(ih,jh,nt)
                 deeq(jh,ih,na,is) = deeq(ih,jh,na,is)
                 !
              END DO
              !
           END DO
           !
        END DO
        !
     END DO
     !
  END IF
  !
  DEALLOCATE( aux, qgm_na, qgm, qmod, ylmk0 )
  !
  CALL stop_clock( 'newd' )
  !
  RETURN
  !
  CONTAINS
    !
    !------------------------------------------------------------------------
    SUBROUTINE newd_so()
      !------------------------------------------------------------------------
      !
      USE spin_orb, ONLY : fcoef
      !
      IMPLICIT NONE
      !
      INTEGER :: ijs, is1, is2, kh, lh
      !
      !
      DO na = 1, nat
         !
         nt  = ityp(na)
         ijs = 0
         !
         DO is1 = 1, 2
            !
            DO is2 =1, 2
               !
               ijs = ijs + 1
               !
               DO ih = 1, nh(nt)
                  !
                  DO jh = 1, nh(nt)
                     !
                     deeq_nc(ih,jh,na,ijs) = dvan_so(ih,jh,ijs,nt)
                     !
                     DO kh = 1, nh(nt)
                        !
                        DO lh = 1, nh(nt)
                           !
                           deeq_nc(ih,jh,na,ijs) = deeq_nc(ih,jh,na,ijs) +   &
                                deeq (kh,lh,na,1)*            &
                                (fcoef(ih,kh,is1,1,nt)*fcoef(lh,jh,1,is2,nt)  + &
                                fcoef(ih,kh,is1,2,nt)*fcoef(lh,jh,2,is2,nt)) + &
                                deeq (kh,lh,na,2)*            &
                                (fcoef(ih,kh,is1,1,nt)*fcoef(lh,jh,2,is2,nt)  + &
                                fcoef(ih,kh,is1,2,nt)*fcoef(lh,jh,1,is2,nt)) + &
                                (0.D0,-1.D0)*deeq (kh,lh,na,3)*            &
                                (fcoef(ih,kh,is1,1,nt)*fcoef(lh,jh,2,is2,nt)  - &
                                fcoef(ih,kh,is1,2,nt)*fcoef(lh,jh,1,is2,nt)) + &
                                deeq (kh,lh,na,4)*            &
                                (fcoef(ih,kh,is1,1,nt)*fcoef(lh,jh,1,is2,nt)  - &
                                fcoef(ih,kh,is1,2,nt)*fcoef(lh,jh,2,is2,nt))   
                           !
                        END DO
                        !
                     END DO
                     !
                  END DO
                  !
               END DO
               !
            END DO
            !
         END DO
         !
      END DO
      !
      RETURN
      !
    END SUBROUTINE newd_so
    !
    !------------------------------------------------------------------------
    SUBROUTINE newd_nc()
      !------------------------------------------------------------------------
      !
      IMPLICIT NONE
      !
      !
      DO na = 1, nat
         !
         nt = ityp(na)
         !
         DO is = 1, nspin
            !
            DO ih = 1, nh(nt)
               !
               DO jh = 1, nh(nt)
                  !
                  deeq_nc(ih,jh,na,1) = dvan(ih,jh,nt) + &
                                        deeq(ih,jh,na,1) + deeq(ih,jh,na,4)
                  !                      
                  deeq_nc(ih,jh,na,2) = deeq(ih,jh,na,2) - &
                                        ( 0.D0, 1.D0 ) * deeq(ih,jh,na,3)
                  !                      
                  deeq_nc(ih,jh,na,3) = deeq(ih,jh,na,2) + &
                                        ( 0.D0, 1.D0 ) * deeq(ih,jh,na,3)
                  !                      
                  deeq_nc(ih,jh,na,4) = dvan(ih,jh,nt) + &
                                        deeq(ih,jh,na,1) - deeq(ih,jh,na,4)
                  !
               END DO
               !
            END DO
            !
         END DO
         !
      END DO
      !
    END SUBROUTINE newd_nc
    !
END SUBROUTINE newd_g
