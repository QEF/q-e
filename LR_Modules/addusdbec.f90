!
! Copyright (C) 2001-2008 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------
subroutine addusdbec (ik, wgt, psi, dbecsum)
  !----------------------------------------------------------------------
  !
  !  This routine adds to dbecsum the contribution of this
  !  k point. It implements Eq. B15 of PRB 64, 235118 (2001).
  !
  USE kinds, only : DP
  USE ions_base, ONLY : nat, ityp, ntyp => nsp
  USE becmod, ONLY : calbec
  USE wvfct, only: npw, npwx, nbnd
  USE uspp, only: nkb, vkb, okvan, ijtoh
  USE uspp_param, only: upf, nh, nhm

  USE lrus,   ONLY : becp1
  USE qpoint, ONLY : npwq, ikks
  USE control_lr, ONLY : nbnd_occ
  !
  USE mp_bands, ONLY : intra_bgrp_comm
  !
  implicit none
  !
  !   the dummy variables
  !
  complex(DP) :: dbecsum (nhm*(nhm+1)/2, nat), psi(npwx,nbnd)
  ! inp/out: the sum kv of bec *
  ! input  : contains delta psi
  integer :: ik
  ! input: the k point
  real(DP) :: wgt
  ! input: the weight of this k point
  !
  !     here the local variables
  !
  integer :: na, nt, ih, jh, ibnd, ikk, ikb, jkb, ijh, startb, &
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

  complex(DP), allocatable :: dbecq (:,:)
  ! the change of becq

  if (.not.okvan) return

  call start_clock ('addusdbec')

  allocate (dbecq( nkb, nbnd))
  ikk = ikks(ik)
  !
  !     First compute the product of psi and vkb
  !
  call calbec (npwq, vkb, psi, dbecq)
  !
  !  And then we add the product to becsum
  !
  !  Band parallelization: each processor takes care of its slice of bands
  !
  call divide (intra_bgrp_comm, nbnd_occ (ikk), startb, lastb)
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
  deallocate (dbecq)

  call stop_clock ('addusdbec')
  return
end subroutine addusdbec
