!
! Copyright (C) 2001-2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "f_defs.h"
!
!----------------------------------------------------------------------------
SUBROUTINE add_shift_us( shift_nl )
  !----------------------------------------------------------------------------
  !
  ! ... nonlocal potential contribution to forces
  ! ... wrapper
  !
  USE kinds,                ONLY : DP
  USE wvfct,                ONLY : g2kin
  USE control_flags,        ONLY : gamma_only
  USE cell_base,            ONLY : at, bg, tpiba, tpiba2
  USE ions_base,            ONLY : nat, ntyp => nsp , ityp
  USE klist,                ONLY : nks, xk
  USE gvect,                ONLY : g, ngm, ecutwfc
  USE uspp,                 ONLY : nkb, vkb, qq, deeq
  USE uspp_param,           ONLY : upf, nh, newpseudo
  USE wvfct,                ONLY : nbnd, npw, npwx, igk, wg, et
  USE lsda_mod,             ONLY : lsda, current_spin, isk
  USE symme,                ONLY : irt, s, nsym
  USE wavefunctions_module, ONLY : evc
  USE io_files,             ONLY : iunwfc, nwordwfc
  USE becmod,               ONLY : calbec
  !
  IMPLICIT NONE
  !
  ! ... the dummy variable
  !
  REAL(DP) :: shift_nl(nat)
  ! output: the nonlocal contribution
  !
  !
  IF ( gamma_only ) THEN
     !
     CALL add_shift_us_gamma()
     !
  ELSE
     !
     CALL add_shift_us_k()
     !
  END IF  
  !
  RETURN
  !
  CONTAINS
     !
     !-----------------------------------------------------------------------
     SUBROUTINE add_shift_us_gamma()
       !-----------------------------------------------------------------------
       !
       ! ... calculation at gamma
       !
       IMPLICIT NONE
       !
       REAL(DP), ALLOCATABLE    :: rbecp(:,:), shift_(:)
       ! auxiliary variables contain <beta|psi> 
       REAL(DP) :: ps
       INTEGER       :: ik, ibnd, ig, ih, jh, na, nt, ikb, jkb, ijkb0
       ! counters
       !
       !
       !
       ALLOCATE( rbecp( nkb, nbnd ), shift_(nat) )    
       !   
       shift_(:) = 0.d0
       !
       ! ... the forces are a sum over the K points and the bands
       !
       DO ik = 1, nks
          IF ( lsda ) current_spin = isk(ik)
          !
          call gk_sort (xk(1,ik), ngm, g, ecutwfc / tpiba2, npw, igk, g2kin)
          IF ( nks > 1 ) THEN
             CALL davcio( evc, nwordwfc, iunwfc, ik, -1 )
             IF ( nkb > 0 ) CALL init_us_2( npw, igk, xk(1,ik), vkb )
          END IF
          !
          CALL calbec ( npw, vkb, evc, rbecp )
          !
          ijkb0 = 0
          DO nt = 1, ntyp
             DO na = 1, nat
                IF ( ityp(na) == nt ) THEN
                   DO ih = 1, nh(nt)
                      ikb = ijkb0 + ih
                      DO ibnd = 1, nbnd
                         ps = deeq(ih,ih,na,current_spin) - &
                              et(ibnd,ik) * qq(ih,ih,nt)
                         shift_(na) = shift_(na) + ps * wg(ibnd,ik) * &
                                      rbecp(ikb,ibnd) * rbecp(ikb,ibnd)
                      END DO
                      !
                      IF ( upf(nt)%tvanp .OR. newpseudo(nt) ) THEN
                         !
                         ! ... in US case there is a contribution for jh<>ih. 
                         ! ... We use here the symmetry in the interchange 
                         ! ... of ih and jh
                         !
                         DO jh = ( ih + 1 ), nh(nt)
                            jkb = ijkb0 + jh
                            DO ibnd = 1, nbnd
                               ps = deeq(ih,jh,na,current_spin) - &
                                    et(ibnd,ik) * qq(ih,jh,nt)
                               shift_(na) = shift_(na) + ps * wg(ibnd,ik) * &
                                     2.d0 *rbecp(ikb,ibnd) *rbecp(jkb,ibnd) 
                            END DO
                         END DO
                      END IF
                   END DO
                   ijkb0 = ijkb0 + nh(nt)
                END IF
             END DO
          END DO
       END DO
       !
#ifdef __PARA
       !
       ! ... collect contributions across pools
       !
       CALL poolreduce( nat, shift_ )
#endif
       !
       ! ... Since our summation over k points was only on the irreducible 
       ! ... BZ we have to symmetrize the shifts. 
       !
       CALL symscalar( nat, shift_, nsym, s, irt )
       !
       shift_nl(:) = shift_nl(:) + shift_(:)
       !
       DEALLOCATE( rbecp, shift_ ) 
       !
       RETURN
       !
     END SUBROUTINE add_shift_us_gamma
     !     
     !-----------------------------------------------------------------------
     SUBROUTINE add_shift_us_k()
       !-----------------------------------------------------------------------
       !  
       IMPLICIT NONE
       !
       REAL(DP), ALLOCATABLE :: shift_(:)
       ! auxiliary variable 
       COMPLEX(DP), ALLOCATABLE :: becp(:,:) 
       !  contains products of wavefunctions and beta

       REAL(DP) :: ps
       INTEGER       :: ik, ibnd, ig, ih, jh, na, nt, ikb, jkb, ijkb0
       ! counters
       !
       ALLOCATE( becp(nkb,nbnd), shift_( nat ) )
       shift_(:) = 0.D0
       ! 
       ! ... the shifts are a sum over the K points and the bands
       !
       DO ik = 1, nks
          IF ( lsda ) current_spin = isk(ik)
          !
          call gk_sort (xk(1,ik), ngm, g, ecutwfc / tpiba2, npw, igk, g2kin)
          IF ( nks > 1 ) THEN
             CALL davcio( evc, nwordwfc, iunwfc, ik, -1 )
             IF ( nkb > 0 ) CALL init_us_2( npw, igk, xk(1,ik), vkb )
          END IF
          !
          CALL calbec( npw, vkb, evc, becp )
          !
          !
          ijkb0 = 0
          DO nt = 1, ntyp
             DO na = 1, nat
                IF ( ityp(na) == nt ) THEN
                   DO ih = 1, nh(nt)
                      ikb = ijkb0 + ih
                      DO ibnd = 1, nbnd
                         ps = deeq(ih,ih,na,current_spin) - &
                              et(ibnd,ik) * qq(ih,ih,nt)
                         shift_(na) = shift_(na) + ps * wg(ibnd,ik) * &
                                      DBLE( CONJG( becp(ikb,ibnd) ) * &
                                                   becp(ikb,ibnd) )
                      END DO
                      !
                      IF ( upf(nt)%tvanp .OR. newpseudo(nt) ) THEN
                         !
                         ! ... in US case there is a contribution for jh<>ih. 
                         ! ... We use here the symmetry in the interchange 
                         ! ... of ih and jh
                         !
                         DO jh = ( ih + 1 ), nh(nt)
                            jkb = ijkb0 + jh
                            DO ibnd = 1, nbnd
                               ps = deeq(ih,jh,na,current_spin) - &
                                    et(ibnd,ik) * qq (ih,jh,nt)
                               shift_(na) = shift_ (na) + ps * wg(ibnd,ik) * &
                                      2.d0 * DBLE( CONJG( becp(ikb,ibnd) ) * &
                                                          becp(jkb,ibnd) )
                            END DO
                         END DO
                      END IF
                   END DO
                   ijkb0 = ijkb0 + nh(nt)
                END IF
             END DO
          END DO
       END DO
       !
#ifdef __PARA
       !
       ! ... collect contributions across pools
       !
       CALL poolreduce( nat, shift_ )
#endif
       !
       ! ... Since our summation over k points was only on the irreducible 
       ! ... BZ we have to symmetrize the forces. 
       !
       CALL symscalar( nat, shift_, nsym, s, irt )
       !
       shift_nl(:) = shift_nl(:) + shift_(:)

       DEALLOCATE( shift_ , becp)
       !
       RETURN
       !
     END SUBROUTINE add_shift_us_k
     !     
END SUBROUTINE add_shift_us
