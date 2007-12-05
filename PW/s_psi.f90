!
! Copyright (C) 2001-2007 Quantum-Espresso group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "f_defs.h"
!
!----------------------------------------------------------------------------
SUBROUTINE s_psi( lda, n, m, psi, spsi )
  !----------------------------------------------------------------------------
  !
  ! ... This routine applies the S matrix to m wavefunctions psi
  ! ... and puts the results in spsi.
  ! ... Requires the products of psi with all beta functions
  ! ... in array becp(nkb,m) (calculated in h_psi or by calbec)
  !
  ! ... input:
  !
  ! ...    lda   leading dimension of arrays psi, spsi
  ! ...    n     true dimension of psi, spsi
  ! ...    m     number of states psi
  ! ...    psi
  !
  ! ... output:
  !
  ! ...    spsi  S*psi
  !
  USE kinds,      ONLY : DP
  USE uspp,       ONLY : vkb, nkb, qq, okvan
  USE uspp_param, ONLY : upf, nh 
  USE wvfct,      ONLY : igk, g2kin
  USE gsmooth,    ONLY : nls, nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, nrxxs
  USE ldaU,       ONLY : lda_plus_u
  USE ions_base,  ONLY : nat, nsp, ityp
  USE control_flags,    ONLY: gamma_only 
  USE noncollin_module, ONLY: npol, noncolin
  !
  IMPLICIT NONE
  !
  ! ... First the dummy variables
  !
  INTEGER          :: lda, n, m
  COMPLEX(DP) :: psi(lda*npol,m), spsi(lda*npol,m)
  !
  CALL start_clock( 's_psi' )  
  !
  IF ( gamma_only ) THEN
     !
     CALL s_psi_gamma()
     !
  ELSE IF ( noncolin ) THEN
     !
     ! ... Unlike s_psi_gamma and s_psi_k, this is an external routine
     !
     CALL s_psi_nc( lda, n, m, psi, spsi )
     !
  ELSE 
     !
     CALL s_psi_k()
     !
  END IF    
  !
  CALL stop_clock( 's_psi' )
  !
  RETURN
  !
  CONTAINS
     !
     !-----------------------------------------------------------------------
     SUBROUTINE s_psi_gamma()
       !-----------------------------------------------------------------------
       ! 
       ! ... gamma version
       !
       USE becmod, ONLY : rbecp
       !
       IMPLICIT NONE  
       !
       ! ... here the local variables
       !
       INTEGER :: ikb, jkb, ih, jh, na, nt, ijkb0, ibnd
         ! counters
       REAL(DP), ALLOCATABLE :: ps(:,:)
         ! the product vkb and psi
       !
       !
       ! ... initialize  spsi
       !
       spsi = psi
       !
       ! ... The product with the beta functions
       !
       IF ( nkb == 0 .OR. .NOT. okvan ) RETURN
       !
       ALLOCATE( ps( nkb, m ) )
       !    
       ps(:,:) = 0.D0
       !
       ijkb0 = 0
       DO nt = 1, nsp
          IF ( upf(nt)%tvanp ) THEN
             DO na = 1, nat
                IF ( ityp(na) == nt ) THEN
                   DO ibnd = 1, m
                      DO jh = 1, nh(nt)
                         jkb = ijkb0 + jh
                         DO ih = 1, nh(nt)
                            ikb = ijkb0 + ih
                            ps(ikb,ibnd) = ps(ikb,ibnd) + &
                                           qq(ih,jh,nt) * rbecp(jkb,ibnd)
                         END DO
                      END DO
                   END DO
                   ijkb0 = ijkb0 + nh(nt)
                END IF
             END DO
          ELSE
             DO na = 1, nat
                IF ( ityp(na) == nt ) ijkb0 = ijkb0 + nh(nt)
             END DO
          END IF
       END DO
       !
       IF ( m == 1 ) THEN
          !
          CALL DGEMV( 'N', 2 * n, nkb, 1.D0, vkb, &
                      2 * lda, ps, 1, 1.D0, spsi, 1 )
          !
       ELSE
          !
          CALL DGEMM( 'N', 'N', 2 * n, m, nkb, 1.D0, vkb, &
                      2 * lda, ps, nkb, 1.D0, spsi, 2 * lda )
          !
       END IF
       !
       DEALLOCATE( ps ) 
       !
       RETURN
       !
     END SUBROUTINE s_psi_gamma
     !
     !-----------------------------------------------------------------------
     SUBROUTINE s_psi_k()
       !-----------------------------------------------------------------------
       !
       ! ... k-points version
       !
       USE becmod,  ONLY : becp
       !
       IMPLICIT NONE
       !
       ! ... local variables
       !
       INTEGER :: ikb, jkb, ih, jh, na, nt, ijkb0, ibnd
         ! counters
       COMPLEX(DP), ALLOCATABLE :: ps(:,:)
         ! the product vkb and psi
       !
       !
       ! ... initialize  spsi
       !
       spsi = psi
       !
       ! ... The product with the beta functions
       !
       IF ( nkb == 0 .OR. .NOT. okvan ) RETURN
       !
       ALLOCATE( ps( nkb, m ) )    
       !
       ps(:,:) = ( 0.D0, 0.D0 )
       !
       ijkb0 = 0
       DO nt = 1, nsp
          IF ( upf(nt)%tvanp ) THEN
             DO na = 1, nat
                IF ( ityp(na) == nt ) THEN
                   DO ibnd = 1, m
                      DO jh = 1, nh(nt)
                         jkb = ijkb0 + jh
                         DO ih = 1, nh(nt)
                            ikb = ijkb0 + ih
                            ps(ikb,ibnd) = ps(ikb,ibnd) + &
                                           qq(ih,jh,nt) * becp(jkb,ibnd)
                         END DO
                      END DO
                   END DO
                   ijkb0 = ijkb0 + nh(nt)
                END IF
             END DO
          ELSE
             DO na = 1, nat
                IF ( ityp(na) == nt ) ijkb0 = ijkb0 + nh(nt)
             END DO
          END IF
       END DO
       !
       IF ( m == 1 ) THEN
          !
          CALL ZGEMV( 'N', n, nkb, ( 1.D0, 0.D0 ), vkb, &
                      lda, ps, 1, ( 1.D0, 0.D0 ), spsi, 1 )
          !
       ELSE
          !
          CALL ZGEMM( 'N', 'N', n, m, nkb, ( 1.D0, 0.D0 ), vkb, &
                      lda, ps, nkb, ( 1.D0, 0.D0 ), spsi, lda )
          !
       END IF
       !
       DEALLOCATE( ps )
       !
       RETURN
       !
     END SUBROUTINE s_psi_k     
     !
END subroutine s_psi
!
!-----------------------------------------------------------------------
subroutine s_psi_nc (lda, n, m, psi, spsi )
  !-----------------------------------------------------------------------
  !
  !    This routine applies the S matrix to m wavefunctions psi
  !    and puts the results in spsi.
  !    Requires the products of psi with all beta functions
  !    in array becp_nc(nkb,m) (calculated by calbec)
  ! input:
  !     lda   leading dimension of arrays psi, spsi
  !     n     true dimension of psi, spsi
  !     m     number of states psi
  !     psi
  ! output:
  !     spsi  S*psi
  !
  USE kinds,     ONLY : DP
  USE ions_base, ONLY: nat, ityp, ntyp => nsp
  USE uspp_param, ONLY: upf, nh
  USE uspp, ONLY: nkb, vkb, qq, qq_so, okvan
  use wvfct, ONLY: igk, g2kin
  use gsmooth, ONLY: nls, nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, nrxxs
  use ldaU, ONLY: lda_plus_u
  use becmod
  use noncollin_module, ONLY: npol
  use spin_orb, ONLY: lspinorb
  implicit none
  !
  !     First the dummy variables
  !
  integer :: lda, n, m
  complex(DP) :: psi (lda, npol, m), spsi (lda,npol, m)
  complex(DP), external ::  ZDOTU
  !
  !    here the local variables
  !
  integer :: ikb, jkb, ih, jh, na, nt, ijkb0, ibnd, ipol
  ! counters
  complex(DP), allocatable :: ps (:,:,:)
  ! the product vkb and psi
  !
  !   initialize  spsi
  !
  call ZCOPY (lda * m *npol, psi(1,1,1), 1, spsi(1,1,1), 1)
  !
  !  The product with the beta functions
  !
  if (nkb.eq.0.or..not.okvan) return
  !
  allocate (ps(nkb,npol,m))    
  ps(:,:,:) = (0.D0,0.D0)
  !
  ijkb0 = 0
  do nt = 1, ntyp
     if ( upf(nt)%tvanp ) then
        do na = 1, nat
           if (ityp (na) .eq.nt) then
              do ih = 1,nh(nt)
                 ikb = ijkb0 + ih
                 do ibnd = 1, m
                    do jh = 1, nh (nt)
                       jkb = ijkb0 + jh
                       if (lspinorb) then
                          ps(ikb,1,ibnd)=ps(ikb,1,ibnd) + &
                                 qq_so(ih,jh,1,nt)*becp_nc(jkb,1,ibnd)+ &
                                 qq_so(ih,jh,2,nt)*becp_nc(jkb,2,ibnd)
                          ps(ikb,2,ibnd)=ps(ikb,2,ibnd) + &
                                 qq_so(ih,jh,3,nt)*becp_nc(jkb,1,ibnd)+ &
                                 qq_so(ih,jh,4,nt)*becp_nc(jkb,2,ibnd)
                       else
                          do ipol=1,npol
                             ps(ikb,ipol,ibnd)=ps(ikb,ipol,ibnd) + &
                                        qq(ih,jh,nt)*becp_nc(jkb,ipol,ibnd)
                          enddo
                       endif
                    enddo
                 enddo
              enddo
              ijkb0 = ijkb0 + nh (nt)
           endif
        enddo
     else
        do na = 1, nat
           if (ityp (na) .eq.nt) ijkb0 = ijkb0 + nh (nt)
        enddo
     endif

  enddo
  call ZGEMM ('N', 'N', n, m*npol, nkb, (1.d0, 0.d0) , vkb, &
       lda, ps, nkb, (1.d0, 0.d0) , spsi(1,1,1), lda)

  DEALLOCATE(ps)

  return
end subroutine s_psi_nc
