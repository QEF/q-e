!
! Copyright (C) 2001-2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------------
SUBROUTINE add_shift_us( shift_nl )
  !----------------------------------------------------------------------------
  !
  ! ... nonlocal potential contribution to forces
  ! ... wrapper
  !
  USE kinds,                ONLY : DP
  USE control_flags,        ONLY : gamma_only
  USE cell_base,            ONLY : at, bg
  USE ions_base,            ONLY : nat, ntyp => nsp , ityp
  USE klist,                ONLY : nks, xk, ngk, igk_k
  USE gvect,                ONLY : g, ngm
  USE uspp,                 ONLY : nkb, vkb, qq, deeq
  USE uspp_param,           ONLY : upf, nh, newpseudo
  USE wvfct,                ONLY : nbnd, wg, et
  USE lsda_mod,             ONLY : lsda, isk
  USE symme,                ONLY : symscalar
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
  ENDIF
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
       USE mp_global,            ONLY: inter_pool_comm, intra_pool_comm
       USE mp,                   ONLY: mp_sum

       IMPLICIT NONE
       !
       REAL(DP), ALLOCATABLE    :: rbecp(:,:), shift_(:)
       ! auxiliary variables contain <beta|psi>
       REAL(DP) :: ps
       INTEGER  :: npw, ik, is, ibnd, ih, jh, na, nt, ikb, jkb, ijkb0
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
          !
          is = isk(ik)
          npw = ngk(ik)
          IF ( nks > 1 ) THEN
             CALL davcio( evc, 2*nwordwfc, iunwfc, ik, -1 )
             IF ( nkb > 0 ) CALL init_us_2( npw, igk_k(1,ik), xk(1,ik), vkb )
          ENDIF
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
                         ps = deeq(ih,ih,na,is) - &
                              et(ibnd,ik) * qq(ih,ih,nt)
                         shift_(na) = shift_(na) + ps * wg(ibnd,ik) * &
                                      rbecp(ikb,ibnd) * rbecp(ikb,ibnd)
                      ENDDO
                      !
                      IF ( upf(nt)%tvanp .or. newpseudo(nt) ) THEN
                         !
                         ! ... in US case there is a contribution for jh<>ih.
                         ! ... We use here the symmetry in the interchange
                         ! ... of ih and jh
                         !
                         DO jh = ( ih + 1 ), nh(nt)
                            jkb = ijkb0 + jh
                            DO ibnd = 1, nbnd
                               ps = deeq(ih,jh,na,is) - &
                                    et(ibnd,ik) * qq(ih,jh,nt)
                               shift_(na) = shift_(na) + ps * wg(ibnd,ik) * &
                                     2.d0 *rbecp(ikb,ibnd) *rbecp(jkb,ibnd)
                            ENDDO
                         ENDDO
                      ENDIF
                   ENDDO
                   ijkb0 = ijkb0 + nh(nt)
                ENDIF
             ENDDO
          ENDDO
       ENDDO
       !
#if defined(__MPI)
       !
       ! ... collect contributions across pools
       !
       CALL mp_sum( shift_, inter_pool_comm )
#endif
       !
       ! ... Since our summation over k points was only on the irreducible
       ! ... BZ we have to symmetrize the shifts.
       !
       CALL symscalar( nat, shift_ )
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
       USE mp_global,            ONLY: inter_pool_comm, intra_pool_comm
       USE mp,                   ONLY: mp_sum

       IMPLICIT NONE
       !
       REAL(DP), ALLOCATABLE :: shift_(:)
       ! auxiliary variable
       COMPLEX(DP), ALLOCATABLE :: becp(:,:)
       !  contains products of wavefunctions and beta

       REAL(DP) :: ps
       INTEGER  :: npw, ik, is, ibnd, ih, jh, na, nt, ikb, jkb, ijkb0
       ! counters
       !
       ALLOCATE( becp(nkb,nbnd), shift_( nat ) )
       shift_(:) = 0.D0
       !
       ! ... the shifts are a sum over the K points and the bands
       !
       DO ik = 1, nks
          !
          is = isk(ik)
          npw = ngk(ik)
          IF ( nks > 1 ) THEN
             CALL davcio( evc, 2*nwordwfc, iunwfc, ik, -1 )
             IF ( nkb > 0 ) CALL init_us_2( npw, igk_k(1,ik), xk(1,ik), vkb )
          ENDIF
          !
          CALL calbec( npw, vkb, evc, becp )
          !
          ijkb0 = 0
          DO nt = 1, ntyp
             DO na = 1, nat
                IF ( ityp(na) == nt ) THEN
                   DO ih = 1, nh(nt)
                      ikb = ijkb0 + ih
                      DO ibnd = 1, nbnd
                         ps = deeq(ih,ih,na,is) - &
                              et(ibnd,ik) * qq(ih,ih,nt)
                         shift_(na) = shift_(na) + ps * wg(ibnd,ik) * &
                                      dble( conjg( becp(ikb,ibnd) ) * &
                                                   becp(ikb,ibnd) )
                      ENDDO
                      !
                      IF ( upf(nt)%tvanp .or. newpseudo(nt) ) THEN
                         !
                         ! ... in US case there is a contribution for jh<>ih.
                         ! ... We use here the symmetry in the interchange
                         ! ... of ih and jh
                         !
                         DO jh = ( ih + 1 ), nh(nt)
                            jkb = ijkb0 + jh
                            DO ibnd = 1, nbnd
                               ps = deeq(ih,jh,na,is) - &
                                    et(ibnd,ik) * qq (ih,jh,nt)
                               shift_(na) = shift_ (na) + ps * wg(ibnd,ik) * &
                                      2.d0 * dble( conjg( becp(ikb,ibnd) ) * &
                                                          becp(jkb,ibnd) )
                            ENDDO
                         ENDDO
                      ENDIF
                   ENDDO
                   ijkb0 = ijkb0 + nh(nt)
                ENDIF
             ENDDO
          ENDDO
       ENDDO
       !
#if defined(__MPI)
       !
       ! ... collect contributions across pools
       !
       CALL mp_sum( shift_, inter_pool_comm )
#endif
       !
       ! ... Since our summation over k points was only on the irreducible
       ! ... BZ we have to symmetrize the forces.
       !
       CALL symscalar( nat, shift_ )
       !
       shift_nl(:) = shift_nl(:) + shift_(:)

       DEALLOCATE( shift_ , becp)
       !
       RETURN
       !
     END SUBROUTINE add_shift_us_k
     !
END SUBROUTINE add_shift_us
