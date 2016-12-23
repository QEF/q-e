!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------
subroutine adddvscf (ipert, ik)
  !----------------------------------------------------------------------
  !
  !     This routine computes the contribution of the selfconsistent
  !     change of the potential to the known part of the linear
  !     system and adds it to dvpsi.
  !     It implements the second term in Eq. B30 of PRB 64, 235118 (2001).
  !

  USE kinds,      ONLY : DP
  USE uspp_param, ONLY : upf, nh
  USE uspp,       ONLY : vkb, okvan
! modules from pwcom
  USE lsda_mod,   ONLY : lsda, current_spin, isk
  USE ions_base,  ONLY : ntyp => nsp, nat, ityp
  USE wvfct,      ONLY : nbnd, npwx
  USE noncollin_module, ONLY : noncolin, npol
! modules from phcom
  USE phus,       ONLY : becp1

  USE lrus,       ONLY : int3, int3_nc
  USE qpoint,     ONLY : npwq, ikks
  USE eqv,        ONLY : dvpsi
  USE el_phon,    ONLY : elph_nbnd_min, elph_nbnd_max
  USE elph_tetra_mod, ONLY : elph_tetra

  implicit none
  !
  !   The dummy variables
  !
  integer :: ik, ipert
  ! input: the k point
  ! input: the perturbation
  !
  !   And the local variables
  !
  integer :: na, nt, ibnd, ih, jh, ijkb0, ikk, ikb, jkb, is, js, ijs
  integer :: ibnd_fst, ibnd_lst
  ! counter on atoms
  ! counter on atomic types
  ! counter on bands
  ! counter on beta functions
  ! counter on beta functions
  ! auxiliary variable for indexing
  ! counter on the k points
  ! counter on vkb
  ! counter on vkb
  complex(DP) :: sum, sum_nc(npol)
  ! auxiliary variable

  if(elph_tetra == 0) then
     ibnd_fst = 1
     ibnd_lst = nbnd
  else
     ibnd_fst = elph_nbnd_min
     ibnd_lst = elph_nbnd_max
  end if
  
  if (.not.okvan) return
  call start_clock ('adddvscf')
  ikk = ikks(ik)
  if (lsda) current_spin = isk (ikk)
  ijkb0 = 0
  do nt = 1, ntyp
     if (upf(nt)%tvanp  ) then
        do na = 1, nat
           if (ityp (na) .eq.nt) then
              !
              !   we multiply the integral for the becp term and the beta_n
              !
              do ibnd = ibnd_fst, ibnd_lst
                 do ih = 1, nh (nt)
                    ikb = ijkb0 + ih
                    IF (noncolin) THEN
                       sum_nc = (0.d0, 0.d0)
                    ELSE
                       sum = (0.d0, 0.d0)
                    END IF
                    do jh = 1, nh (nt)
                       jkb = ijkb0 + jh
                       IF (noncolin) THEN
                          ijs=0
                          do is=1,npol
                             do js=1,npol
                                ijs=ijs+1
                                sum_nc(is)=sum_nc(is)+               &
                                     int3_nc(ih,jh,na,ijs,ipert)*    &
                                     becp1(ik)%nc(jkb, js, ibnd)
                             enddo
                          enddo
                       ELSE
                          sum = sum + int3 (ih, jh, na, current_spin, ipert)*&
                                   becp1(ik)%k(jkb, ibnd)
                       END IF
                    enddo
                    IF (noncolin) THEN
                       call zaxpy(npwq,sum_nc(1),vkb(1,ikb),1,dvpsi(1,ibnd),1)
                       call zaxpy(npwq,sum_nc(2),vkb(1,ikb),1, &
                                                 dvpsi(1+npwx,ibnd),1)
                    ELSE
                       call zaxpy(npwq,sum,vkb(1,ikb),1,dvpsi(1,ibnd),1)
                    END IF
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

  call stop_clock ('adddvscf')
  return
end subroutine adddvscf
