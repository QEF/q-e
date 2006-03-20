
!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!--------------------------------------------------------------------
subroutine local_dos1d (ik, kband, plan)
  !--------------------------------------------------------------------
  !
  !     calculates |psi|^2 for band kband at point ik
  !
#include "f_defs.h"
  USE cell_base, ONLY: omega
  USE ions_base, ONLY: nat, ntyp=>nsp, ityp
  USE gvect
  USE gsmooth
  USE lsda_mod, ONLY: current_spin
  USE uspp, ONLY: becsum, indv, nhtol, nhtoj
  USE uspp_param, ONLY: nh, tvanp, nhm
  USE wvfct, ONLY: npw, wg, igk
  USE noncollin_module, ONLY: noncolin, npol
  USE spin_orb, ONLY: lspinorb, so, fcoef
  USE wavefunctions_module,  ONLY: evc, psic, evc_nc, psic_nc
  USE becmod, ONLY: becp, becp_nc
  implicit none
  !
  ! input variables
  !
  integer :: ik, kband
  ! input: the k point
  ! input: the band

  real(DP) :: plan (nr3)
  ! output: the planar average of this state
  !
  !    Additional local variables for Ultrasoft PP's
  !

  integer :: ikb, jkb, ijkb0, ih, jh, na, ijh, ipol, np
  ! counter on beta functions
  ! counter on beta functions
  ! auxiliary variable for ijkb0
  ! counter on solid beta functions
  ! counter on solid beta functions
  ! counter on atoms
  ! counter on composite beta functions
  ! the pseudopotential
  !
  !    And here the local variables
  !
  integer :: ir, is, ig, ibnd, is1, is2, kkb, kh
  ! counter on 3D r points
  ! counter on spin polarizations
  ! counter on g vectors
  ! counter on bands

  real(DP) :: w, w1
  ! the weight of one k point
  real(DP), allocatable :: aux (:)
  ! auxiliary for rho

  complex(DP), allocatable :: prho (:), be1(:,:), be2(:,:)
  ! complex charge for fft

  allocate (prho(nrxx))    
  allocate (aux(nrxx))    
  if (lspinorb) then
     allocate(be1(nhm,2))
     allocate(be2(nhm,2))
  endif

  aux(:) = 0.d0
  becsum(:,:,:) = 0.d0

  wg (kband, ik) = 1.d0
  !
  !
  !     First compute the square modulus of the state kband,ik on the smooth
  !     mesh
  !
  if (noncolin) then
     psic_nc = (0.d0,0.d0)
     do ipol=1,npol
        do ig = 1, npw
           psic_nc (nls (igk (ig) ), ipol ) = evc_nc (ig, ipol, kband)
        enddo
        call cft3s (psic_nc(1,ipol), nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, 2)
     enddo

     w1 = wg (kband, ik) / omega
     do ipol=1,npol
        do ir = 1, nrxxs
           aux(ir) = aux(ir) + w1 * (DBLE(psic_nc(ir,ipol))**2 + &
                                     AIMAG(psic_nc(ir,ipol))**2)
        enddo
     enddo
  else
     psic(1:nrxxs) = (0.d0,0.d0)
     do ig = 1, npw
        psic (nls (igk (ig) ) ) = evc (ig, kband)
     enddo
     call cft3s (psic, nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, 2)

     w1 = wg (kband, ik) / omega
     do ir = 1, nrxxs
        aux(ir) = aux(ir) + w1 * (DBLE(psic(ir))**2 + AIMAG(psic(ir))**2)
     enddo
  endif

  !
  !    If we have a US pseudopotential we compute here the becsum term
  !
  ibnd = kband

  w1 = wg (ibnd, ik)
  ijkb0 = 0
  do np = 1, ntyp
     if (tvanp (np) ) then
        do na = 1, nat
           if (ityp (na) == np) then
              if (noncolin) then
                 if (so(np)) then
                    be1=(0.d0,0.d0)
                    be2=(0.d0,0.d0)
                    do ih = 1, nh(np)
                       ikb = ijkb0 + ih
                       do kh = 1, nh(np)
                          if ((nhtol(kh,np).eq.nhtol(ih,np)).and. &
                              (nhtoj(kh,np).eq.nhtoj(ih,np)).and. &
                              (indv(kh,np).eq.indv(ih,np))) then
                             kkb=ijkb0 + kh
                             do is1=1,2
                                do is2=1,2
                                   be1(ih,is1)=be1(ih,is1)+ &
                                        fcoef(ih,kh,is1,is2,np)* &
                                        becp_nc(kkb,is2,ibnd)
                                   be2(ih,is1)=be2(ih,is1)+ &
                                        fcoef(kh,ih,is2,is1,np)* &
                                     CONJG(becp_nc(kkb,is2,ibnd))
                                enddo
                             enddo
                          endif
                       enddo
                    enddo
                 endif
              endif
              ijh = 1
              do ih = 1, nh (np)
                 ikb = ijkb0 + ih
                 if (noncolin) then
                    if (so(np)) then
                        becsum(ijh,na,1)=becsum(ijh,na,1)+ w1*    &
                            (be1(ih,1)*be2(ih,1)+be1(ih,2)*be2(ih,2))
                    else
                       do ipol=1,npol
                          becsum(ijh,na,current_spin) = &
                             becsum(ijh,na,current_spin) + w1 * &
                               DBLE( CONJG(becp_nc(ikb,ipol,ibnd)) * &
                                           becp_nc(ikb,ipol,ibnd) )
                       enddo
                    endif
                 else
                    becsum(ijh,na,current_spin) = &
                        becsum(ijh,na,current_spin) + w1 * &
                        DBLE( CONJG(becp(ikb,ibnd)) * becp(ikb,ibnd) )
                 endif
                 ijh = ijh + 1
                 do jh = ih + 1, nh (np)
                    jkb = ijkb0 + jh
                    if (noncolin) then
                       if (so(np)) then
                          becsum(ijh,na,1)=becsum(ijh,na,1) &
                              + w1*((be1(jh,1)*be2(ih,1)+   &
                                     be1(jh,2)*be2(ih,2))+  &
                                    (be1(ih,1)*be2(jh,1)+   &
                                     be1(ih,2)*be2(jh,2)) )
                       else
                          do ipol=1,npol
                             becsum(ijh,na,current_spin) = &
                                becsum(ijh,na,current_spin) + w1 * 2.d0 * &
                                DBLE( CONJG(becp_nc(ikb,ipol,ibnd))  &
                                        * becp_nc(jkb,ipol,ibnd) )
                          enddo
                       endif
                    else
                       becsum(ijh,na,current_spin) = &
                           becsum(ijh,na,current_spin) + w1 * 2.d0 * &
                           DBLE( CONJG(becp(ikb,ibnd)) * becp(jkb,ibnd) )
                    endif
                    ijh = ijh + 1
                 enddo
              enddo
              ijkb0 = ijkb0 + nh (np)
           endif
        enddo
     else
        do na = 1, nat
           if (ityp (na) .eq.np) ijkb0 = ijkb0 + nh (np)
        enddo
     endif
  enddo
  !
  !    Interpolate on the thick mesh and pass to reciprocal space
  !
  if (doublegrid) then
     call interpolate (aux, aux, 1)
  endif
  do ir = 1, nrxx
     prho (ir) = CMPLX (aux (ir), 0.d0)
  enddo
  call cft3 (prho, nr1, nr2, nr3, nrx1, nrx2, nrx3, - 1)
  !
  !    Here we add the US contribution to the charge for the atoms which n
  !    it. Or compute the planar average in the NC case.
  !
  call addusdens1d (plan, prho)
  !
  deallocate (aux)
  deallocate (prho)
  if (lspinorb) then
     deallocate(be1)
     deallocate(be2)
  endif
  !
  return
end subroutine local_dos1d
