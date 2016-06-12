!
! Copyright (C) 2001-2016 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------
subroutine addusdbec (ik, wgt, dpsi, dbecsum)
  !----------------------------------------------------------------------
  !
  !  This routine adds to dbecsum the contribution of this
  !  k point. It implements Eq. B15 of PRB 64, 235118 (2001).
  !
  USE kinds,      ONLY : DP
  USE ions_base,  ONLY : nat, ityp, ntyp => nsp
  USE becmod,     ONLY : calbec
  USE wvfct,      ONLY : npwx, nbnd
  USE uspp,       ONLY : nkb, vkb, okvan, ijtoh
  USE uspp_param, ONLY : upf, nh, nhm
  USE mp_bands,   ONLY : intra_bgrp_comm
  USE klist,      ONLY : ngk
  USE lrus,       ONLY : becp1
  USE qpoint,     ONLY : ikks, ikqs
  USE control_lr, ONLY : nbnd_occ
  !
  IMPLICIT NONE
  !
  !   the dummy variables
  !
  COMPLEX(DP) :: dbecsum (nhm*(nhm+1)/2, nat), dpsi(npwx,nbnd)
  ! inp/out: the sum kv of bec *
  ! input  : contains delta psi
  INTEGER :: ik
  ! input: the k point
  REAL(DP) :: wgt
  ! input: the weight of this k point
  !
  !     here the local variables
  !
  INTEGER :: na, nt, ih, jh, ibnd, ikb, jkb, ijh, startb, &
       lastb, ijkb0
  ! counter on atoms
  ! counter on atomic type
  ! counter on solid beta functions
  ! counter on solid beta functions
  ! counter on the bands
  ! the real k point
  ! counter on solid becp
  ! counter on solid becp
  ! composite index for dbecsum
  ! divide among processors the sum
  ! auxiliary variable for counting
  INTEGER :: ikk, ikq, npwq
  ! index of the point k
  ! index of the point k+q
  ! number of the plane-waves at point k+q

  COMPLEX(DP), ALLOCATABLE :: dbecq (:,:)
  ! the change of becq

  IF (.NOT.okvan) RETURN
  !
  CALL start_clock ('addusdbec')
  !
  ALLOCATE (dbecq( nkb, nbnd))
  !
  ikk  = ikks(ik)
  ikq  = ikqs(ik)
  npwq = ngk(ikq)
  !
  ! First compute the product of dpsi and vkb
  !
  CALL calbec (npwq, vkb, dpsi, dbecq)
  !
  !  And then we add the product to becsum
  !
  !  Band parallelization: each processor takes care of its slice of bands
  !
  CALL divide (intra_bgrp_comm, nbnd_occ (ikk), startb, lastb)
  !
  ijkb0 = 0
  do nt = 1, ntyp
     if (upf(nt)%tvanp ) then
        do na = 1, nat
           if (ityp (na) .eq.nt) then
              !
              !  And qgmq and becp and dbecq
              !
              do ih = 1, nh(nt)
                 ikb = ijkb0 + ih
                 ijh=ijtoh(ih,ih,nt)
                 do ibnd = startb, lastb
                    dbecsum (ijh, na) = dbecsum (ijh, na) + &
                         wgt * ( CONJG(becp1(ik)%k(ikb,ibnd)) * dbecq(ikb,ibnd) )
                 enddo
                 do jh = ih + 1, nh (nt)
                    ijh=ijtoh(ih,jh,nt)
                    jkb = ijkb0 + jh
                    do ibnd = startb, lastb
                       dbecsum (ijh, na) = dbecsum (ijh, na) + &
                         wgt*( CONJG(becp1(ik)%k(ikb,ibnd))*dbecq(jkb,ibnd) + &
                               CONJG(becp1(ik)%k(jkb,ibnd))*dbecq(ikb,ibnd) )
                    enddo
                    ijh = ijh + 1
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
  !
  DEALLOCATE (dbecq)
  !
  CALL stop_clock ('addusdbec')
  !
  RETURN
  !
end subroutine addusdbec
