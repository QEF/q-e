!
! Copyright (C) 2001-2016 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------
subroutine addusdbec_nc (ik, wgt, dpsi, dbecsum_nc)
  !----------------------------------------------------------------------
  !
  !  This routine adds to the dbecsum the term which correspond to this
  !  k point. After the accumulation the additional part of the charge
  !  is computed in addusddens.
  !
  USE kinds,            ONLY : DP
  USE lsda_mod,         ONLY : nspin
  USE ions_base,        ONLY : nat, ityp, ntyp => nsp
  USE becmod,           ONLY : calbec
  USE wvfct,            ONLY : npwx, nbnd
  USE uspp,             ONLY : nkb, vkb, okvan
  USE noncollin_module, ONLY : noncolin, npol
  USE uspp_param,       ONLY : upf, nh, nhm
  USE mp_bands,         ONLY : intra_bgrp_comm
  USE klist,            ONLY : ngk
  USE lrus,             ONLY : becp1
  USE qpoint,           ONLY : ikks, ikqs
  USE control_lr,       ONLY : nbnd_occ
  !
  IMPLICIT NONE
  !
  !   the dummy variables
  !
  COMPLEX(DP) :: dbecsum_nc (nhm,nhm,nat,nspin), dpsi(npwx*npol,nbnd)
  ! inp/out: the sum kv of bec *
  ! input  : contains delta psi
  INTEGER :: ik
  ! input: the k point
  REAL(DP) :: wgt
  ! input: the weight of this k point
  !
  !     here the local variables
  !
  INTEGER :: na, nt, ih, jh, ibnd, ikb, jkb, startb, &
       lastb, ijkb0, is1, is2, ijs
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

  COMPLEX(DP), ALLOCATABLE :: dbecq_nc(:,:,:)
  ! the change of becq

  IF (.NOT.okvan) RETURN
  !
  CALL start_clock ('addusdbec_nc')
  !
  ALLOCATE (dbecq_nc( nkb,npol, nbnd))
  !
  ikk  = ikks(ik)
  ikq  = ikqs(ik)
  npwq = ngk(ikq)
  !
  ! First compute the product of dpsi and vkb
  !
  CALL calbec (npwq, vkb, dpsi, dbecq_nc)
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
              do ih = 1, nh (nt)
                 ikb = ijkb0 + ih
                 do jh = 1, nh (nt)
                    jkb = ijkb0 + jh
                    DO ibnd = startb, lastb
                       ijs=0
                       DO is1=1,npol
                          DO is2=1,npol
                             ijs=ijs+1
                             dbecsum_nc(ih,jh,na,ijs)=dbecsum_nc(ih,jh,na,ijs)+&
                                wgt*CONJG(becp1(ik)%nc(ikb,is1,ibnd))          &
                                        *dbecq_nc(jkb,is2,ibnd)
                          ENDDO
                       ENDDO
                    ENDDO
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
  DEALLOCATE (dbecq_nc)
  !
  CALL stop_clock ('addusdbec_nc')
  !
  RETURN
  !
end subroutine addusdbec_nc
