!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
subroutine compute_becsum_ph
  !-----------------------------------------------------------------------
  !
  !   This routine computes the becsum term which is used to compute the
  !   change of the charge due to the displacement of the augmentation
  !   term.
  !   It implements Eq. B16 of Ref.[1].
  !   [1] PRB 64, 235118 (2001).
  !
  !
  !
  USE kinds, only : DP
  USE ions_base, ONLY : nat, ityp, ntyp => nsp
  USE lsda_mod,   ONLY : current_spin, isk, lsda
  USE wvfct,      ONLY : nbnd, wg
  USE noncollin_module, ONLY : noncolin, npol
  USE uspp, ONLY: okvan, becsum
  USE uspp_param, ONLY: upf, nh
  USE paw_variables, ONLY : okpaw

  USE phus,       ONLY : alphasum, alphasum_nc, becsum_nc

  USE lrus,       ONLY : becp1
  USE qpoint,     ONLY : nksq, ikks, ikqs
  USE control_ph, ONLY : rec_code_read
  USE control_lr, ONLY : nbnd_occ

  implicit none

  integer :: ik, ikk, ikq, ijkb0, ijh, ikb, jkb, ih, jh, na, nt, ibnd
  ! counter on k points, beta functions, atoms and bands
  integer :: ijs, is1, is2
  real(DP) :: wgg1 ! auxiliary weight

  IF (rec_code_read >= -20.and..not.okpaw) return
  if (.not.okvan) return
  IF (noncolin) becsum_nc = (0.d0,0.d0)
  becsum = 0.d0
  do ik = 1, nksq
     ikk = ikks(ik)
     ikq = ikqs(ik)
     if (lsda) current_spin = isk (ikk)
     ijkb0 = 0
     do nt = 1, ntyp
        if (upf(nt)%tvanp ) then
           do na = 1, nat
              if (ityp (na) == nt) then
                 ijh = 0
                 do ih = 1, nh (nt)
                    ikb = ijkb0 + ih
                    ijh = ijh + 1
                    do ibnd = 1, nbnd_occ (ikk)
                       wgg1 = wg (ibnd, ikk)
                       IF (noncolin) THEN
                          DO is1=1,npol
                             DO is2=1,npol
                                becsum_nc(ijh,na,is1,is2) = &
                                   becsum_nc(ijh,na,is1,is2) + wgg1* &
                                   CONJG(becp1(ik)%nc(ikb,is1,ibnd))* &
                                         becp1(ik)%nc(ikb,is2,ibnd)
                             END DO
                          END DO
                       ELSE
                          becsum(ijh,na,current_spin) = &
                            becsum(ijh,na,current_spin) + wgg1 * &
                         DBLE ( CONJG(becp1(ik)%k(ikb,ibnd)) *   &
                                      becp1(ik)%k(ikb,ibnd) )
                       END IF
                    enddo
                    do jh = ih+1, nh (nt)
                       jkb = ijkb0 + jh
                       ijh = ijh + 1
                       do ibnd = 1, nbnd
                          wgg1 = wg (ibnd, ikk)
                          IF (noncolin) THEN
                             DO is1=1,npol
                                DO is2=1,npol
                                   becsum_nc(ijh,na,is1,is2) = &
                                      becsum_nc(ijh,na,is1,is2)+ wgg1 * &
                                      (CONJG(becp1(ik)%nc(ikb,is1,ibnd)) * &
                                             becp1(ik)%nc(jkb,is2,ibnd)  )
                                END DO
                             END DO
                          ELSE
                             becsum(ijh,na,current_spin) = &
                                 becsum(ijh,na,current_spin)+wgg1 * 2.d0 * &
                                DBLE ( CONJG(becp1(ik)%k(ikb,ibnd)) * &
                                             becp1(ik)%k(jkb,ibnd) )
                          END IF
                       enddo
                    enddo
                 enddo
                 ijkb0 = ijkb0 + nh (nt)
              endif
           enddo
        else
           do na = 1, nat
              if (ityp(na) == nt) ijkb0 = ijkb0 + nh (nt)
           enddo
        endif
     enddo
  enddo

  IF (noncolin.and.okvan) THEN
     DO nt = 1, ntyp
        IF ( upf(nt)%tvanp ) THEN
           DO na = 1, nat
              IF (ityp(na)==nt) THEN
                 IF (upf(nt)%has_so) THEN
                    CALL transform_becsum_so(becsum_nc,becsum,na)
                 ELSE
                    CALL transform_becsum_nc(becsum_nc,becsum,na)
                 END IF
              END IF
           END DO
        END IF
     END DO
  END IF

  !      do na=1,nat
  !         nt=ityp(na)
  !         do ijh=1,nh(nt)*(nh(nt)+1)/2
  !            WRITE( stdout,'(2i5,f20.10)') na, ijh, becsum(ijh,na,1)
  !         enddo
  !      enddo
  !      call stop_ph(.true.)
  return
end subroutine compute_becsum_ph
