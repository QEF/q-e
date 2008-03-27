!
! Copyright (C) 2001-2007 Quantum-Espresso group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "f_defs.h"
!
!--------------------------------------------------------------------
subroutine local_dos (iflag, lsign, kpoint, kband, spin_component, &
                      emin, emax, dos)
  !--------------------------------------------------------------------
  !
  !     iflag=0: calculates |psi|^2 for band "kband" at point "kpoint"
  !     iflag=1: calculates the local density of state at e_fermi
  !              (only for metals)
  !     iflag=2: calculates the local density of  electronic entropy
  !              (only for metals with fermi spreading)
  !     iflag=3: calculates the integral of local dos from "emin" to "emax"
  !              (emin, emax in Ry)
  !
  !     lsign:   if true and k=gamma and iflag=0, write |psi|^2 * sign(psi)
  !     spin_component: for iflag=3 and LSDA calculations only
  !                     0 for up+down dos,  1 for up dos, 2 for down dos
  !
  USE kinds,                ONLY : DP
  USE cell_base,            ONLY : omega, tpiba2
  USE ions_base,            ONLY : nat, ntyp => nsp, ityp
  USE ener,                 ONLY : ef
  USE gvect,                ONLY : nr1, nr2, nr3, nrx1, nrx2, nrx3, nrxx, &
                                   ngm, g, ecutwfc
  USE gsmooth,              ONLY : nls, nlsm, nr1s, nr2s, nr3s, &
                                   nrx1s, nrx2s, nrx3s, nrxxs, doublegrid
  USE klist,                ONLY : lgauss, degauss, ngauss, nks, wk, xk
  USE lsda_mod,             ONLY : lsda, nspin, current_spin, isk
  USE scf,                  ONLY : rho
  USE symme,                ONLY : nsym, s, ftau
  USE uspp,                 ONLY : nkb, vkb, becsum, nhtol, nhtoj, indv
  USE uspp_param,           ONLY : upf, nh, nhm
  USE wavefunctions_module, ONLY : evc, psic, psic_nc
  USE wvfct,                ONLY : nbnd, npwx, npw, igk, wg, et, g2kin
  USE control_flags,        ONLY : gamma_only
  USE noncollin_module,     ONLY : noncolin, npol
  USE spin_orb,             ONLY : lspinorb, so, fcoef
  USE io_files,             ONLY : iunwfc, nwordwfc
  USE mp_global,            ONLY : me_pool, nproc_pool
  USE mp,                   ONLY : mp_bcast, mp_sum
  USE mp_global,            ONLY : inter_pool_comm, intra_pool_comm
  USE becmod,               ONLY : calbec
  implicit none
  !
  ! input variables
  !
  integer, intent(in) :: iflag, kpoint, kband, spin_component
  logical, intent(in) :: lsign
  real(DP), intent(in) :: emin, emax
  !
  real(DP), intent(out) :: dos (nrxx)
  !
  !    local variables
  !
  integer :: ikb, jkb, ijkb0, ih, jh, kh, na, ijh, np
  ! counters for US PPs
  integer :: ir, is, ig, ibnd, ik, irm, isup, isdw, ipol, kkb, is1, is2
  ! counters
  real(DP) :: w, w1, modulus
  real(DP), allocatable :: rbecp(:,:), segno(:), maxmod(:)
  complex(DP), allocatable :: becp(:,:),  &
                                   becp_nc(:,:,:), be1(:,:), be2(:,:)
  integer :: who_calculate, iproc
  complex(DP) :: phase 
  real(DP), external :: w0gauss, w1gauss
  !
  ! input checks
  !
  if (noncolin.and. lsign) call errore('local_dos','not available yet',1)
  if (noncolin.and. gamma_only) call errore('local_dos','not available yet',1)
  !
  if ( (iflag == 0) .and. ( kband < 1 .or. kband > nbnd ) ) &
       call errore ('local_dos', 'wrong band specified', 1)
  if ( (iflag == 0) .and. ( kpoint < 1 .or. kpoint > nks ) ) &
       call errore ('local_dos', 'wrong kpoint specified', 1)
  if (lsign) then
     if (iflag /= 0) call errore ('local_dos', 'inconsistent flags', 1)
     if (sqrt(xk(1,kpoint)**2+xk(2,kpoint)**2+xk(3,kpoint)**2) > 1d-9 )  &
        call errore ('local_dos', 'k must be zero', 1)
  end if
  !
  if (gamma_only) then 
     allocate (rbecp(nkb,nbnd))
  else 
     if (noncolin) then
        allocate (becp_nc(nkb,npol,nbnd))
        if (lspinorb) then
          allocate(be1(nhm,2))
          allocate(be2(nhm,2))
        endif
     else
        allocate (becp(nkb,nbnd))
     endif
  endif
  rho%of_r(:,:) = 0.d0
  dos(:) = 0.d0
  becsum(:,:,:) = 0.d0
  if (lsign) allocate(segno(nrxx))
  !
  !   calculate the correct weights
  !
  if (iflag /= 0 .and. .not.lgauss) call errore ('local_dos', &
       'gaussian broadening needed', 1)
  if (iflag == 2 .and. ngauss /= -99) call errore ('local_dos', &
       ' beware: not using Fermi-Dirac function ',  - ngauss)
  do ik = 1, nks
     do ibnd = 1, nbnd
        if (iflag == 0) then
           wg (ibnd, ik) = 0.d0
        elseif (iflag == 1) then
           wg (ibnd, ik) = wk (ik) * w0gauss ( (ef - et (ibnd, ik) ) &
                / degauss, ngauss) / degauss
        elseif (iflag == 2) then
           wg (ibnd, ik) = - wk (ik) * w1gauss ( (ef - et (ibnd, ik) ) &
                / degauss, ngauss)
        elseif (iflag == 3) then
           if (et (ibnd, ik) <=  emax .and. et (ibnd, ik) >= emin) then
              wg (ibnd, ik) = wk (ik)
           else
              wg (ibnd, ik) = 0.d0
           endif
        else
           call errore ('local_dos', ' iflag not allowed', abs (iflag) )
        endif
     enddo
  enddo
  if (iflag == 0) wg (kband, kpoint) = 1.d0
  !
  !     here we sum for each k point the contribution
  !     of the wavefunctions to the density of states
  !
  do ik = 1, nks
     if (ik == kpoint .or. iflag /= 0) then
        if (lsda) current_spin = isk (ik)
        call gk_sort (xk (1, ik), ngm, g, ecutwfc / tpiba2, npw, igk, g2kin)
        call davcio (evc, nwordwfc, iunwfc, ik, - 1)
        call init_us_2 (npw, igk, xk (1, ik), vkb)

        if (gamma_only) then
           call calbec ( npw, vkb, evc, rbecp )
        else if (noncolin) then
           call calbec ( npw, vkb, evc, becp_nc )
        else
           call calbec ( npw, vkb, evc, becp )
        end if
     !
     !     here we compute the density of states
     !
        do ibnd = 1, nbnd
           if (ibnd == kband .or. iflag /= 0) then
              if (noncolin) then
                 psic_nc = (0.d0,0.d0)
                 do ig = 1, npw
                    psic_nc(nls(igk(ig)),1)=evc(ig     ,ibnd)
                    psic_nc(nls(igk(ig)),2)=evc(ig+npwx,ibnd)
                 enddo
                 do ipol=1,npol
                    call cft3s (psic_nc(1,ipol),nr1s,nr2s,nr3s, &
                                                nrx1s,nrx2s,nrx3s,2)
                 enddo
              else
                 psic(1:nrxxs) = (0.d0,0.d0)
                 do ig = 1, npw
                    psic (nls (igk (ig) ) ) = evc (ig, ibnd)
                 enddo
                 if (gamma_only) then
                    do ig = 1, npw
                       psic (nlsm(igk (ig) ) ) = CONJG(evc (ig, ibnd))
                    enddo
                 end if
                 call cft3s (psic, nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, 2)
              endif
              w1 = wg (ibnd, ik) / omega
!
!  Compute and save the sign of the wavefunction at the gamma point
!
              if (lsign) then
                 if (gamma_only) then
                    !  psi(r) is real by construction
                    segno(1:nrxxs) = DBLE(psic(1:nrxxs))
                 else
                    !  determine the phase factor that makes psi(r) real.
                    allocate(maxmod(nproc_pool))
                    maxmod(me_pool+1)=0.0_DP
                    do ir = 1, nrxxs
                       modulus=abs(psic(ir))
                       if (modulus > maxmod(me_pool+1)) then
                          irm=ir
                          maxmod(me_pool+1)=modulus
                       endif
                    enddo
                    who_calculate=1
#ifdef __PARA
                    call reduce(nproc_pool,maxmod)
                    do iproc=2,nproc_pool
                       if (maxmod(iproc)>maxmod(who_calculate)) &
                          who_calculate=iproc
                    enddo
#endif
                    if (maxmod(who_calculate) < 1.d-10) &
                       call errore('local_dos','zero wavefunction',1)
                    IF (me_pool+1==who_calculate) &
                          phase = psic(irm)/maxmod(who_calculate)
                    deallocate(maxmod)
#ifdef __PARA
                    call mp_bcast(phase,who_calculate-1,intra_pool_comm)
#endif
                    segno(1:nrxxs) = DBLE( psic(1:nrxxs)*CONJG(phase) )
                 endif
                 if (doublegrid) call interpolate (segno, segno, 1)
                 segno(:) = sign( 1.d0, segno(:) )
              endif
              !
              if (noncolin) then
                 do ipol=1,npol
                    do ir=1,nrxxs
                       rho%of_r(ir,current_spin)=rho%of_r(ir,current_spin)+&
                          w1*(DBLE(psic_nc(ir,ipol))**2+ &
                             AIMAG(psic_nc(ir,ipol))**2)
                    enddo
                 enddo
              else
                 do ir=1,nrxxs
                    rho%of_r(ir,current_spin)=rho%of_r(ir,current_spin) + &
                      w1 * (DBLE( psic (ir) ) **2 + AIMAG (psic (ir) ) **2)
                 enddo
              endif
        !
        !    If we have a US pseudopotential we compute here the becsum term
        !
              w1 = wg (ibnd, ik)
              ijkb0 = 0
              do np = 1, ntyp
                if (upf(np)%tvanp  ) then
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
                        ijh = 1
                        do ih = 1, nh (np)
                          ikb = ijkb0 + ih
                          if (so(np)) then
                            becsum(ijh,na,1)=becsum(ijh,na,1)+ w1*    &
                               (be1(ih,1)*be2(ih,1)+be1(ih,2)*be2(ih,2))
                          else
                            becsum(ijh,na,1) = becsum(ijh,na,1)+  &
                             w1*(CONJG(becp_nc(ikb,1,ibnd))*      &
                                       becp_nc(ikb,1,ibnd)+       &
                                 CONJG(becp_nc(ikb,2,ibnd))*      &
                                       becp_nc(ikb,2,ibnd))
                          endif
                          ijh = ijh + 1
                          do jh = ih + 1, nh (np)
                            jkb = ijkb0 + jh
                            if (so(np)) then 
                              becsum(ijh,na,1)=becsum(ijh,na,1) &
                                 + w1*((be1(jh,1)*be2(ih,1)+   &
                                        be1(jh,2)*be2(ih,2))+  &
                                       (be1(ih,1)*be2(jh,1)+   &
                                        be1(ih,2)*be2(jh,2)) )
                            else
                              becsum(ijh,na,1)= becsum(ijh,na,1)+ &
                                   w1*2.d0*DBLE(CONJG(becp_nc(ikb,1,ibnd)) &
                                     *becp_nc(jkb,1,ibnd) + &
                                CONJG(becp_nc(ikb,2,ibnd)) &
                                     *becp_nc(jkb,2,ibnd) )
                            endif
                            ijh = ijh + 1
                          enddo
                        enddo
                      else
                        ijh = 1
                        do ih = 1, nh (np)
                          ikb = ijkb0 + ih
                          if (gamma_only) then
                              becsum(ijh,na,current_spin) = &
                                    becsum(ijh,na,current_spin) + w1 * &
                                    rbecp(ikb,ibnd)*rbecp(ikb,ibnd)
                          else
                              becsum(ijh,na,current_spin) = &
                                   becsum(ijh,na,current_spin) + w1 * &
                               DBLE(CONJG(becp(ikb,ibnd))*becp(ikb,ibnd))
                          end if
                          ijh = ijh + 1
                          do jh = ih + 1, nh (np)
                             jkb = ijkb0 + jh
                             if (gamma_only) then
                                becsum(ijh,na,current_spin) = &
                                   becsum(ijh,na,current_spin) + 2.d0*w1 * &
                                   rbecp(ikb,ibnd)*rbecp(jkb,ibnd)
                             else
                                becsum(ijh,na,current_spin) = &
                                  becsum(ijh,na,current_spin) + 2.d0*w1 * &
                                  DBLE(CONJG(becp(ikb,ibnd))*becp(jkb,ibnd))
                             endif
                             ijh = ijh + 1
                          enddo
                        enddo
                      endif
                      ijkb0 = ijkb0 + nh (np)
                    endif
                  enddo
                else
                  do na = 1, nat
                    if (ityp (na) == np) ijkb0 = ijkb0 + nh (np)
                  enddo
                endif
              enddo
           endif
        enddo
     endif
  enddo
  if (gamma_only) then
     deallocate(rbecp)
  else
     if (noncolin) then
        if (lspinorb) then
           deallocate(be1)
           deallocate(be2)
        endif
        deallocate(becp_nc)
     else
        deallocate(becp)
     endif
  endif
  if (doublegrid) then
     if (noncolin) then
       call interpolate(rho%of_r, rho%of_r, 1)
     else
       do is = 1, nspin
         call interpolate(rho%of_r(1, is), rho%of_r(1, is), 1)
       enddo
     endif
  endif
  !
  !    Here we add the US contribution to the charge
  !
  call addusdens
  !
  if (nspin == 1 .or. nspin==4) then
     is = 1 
     dos(:) = rho%of_r (:, is)
  else
     IF ( iflag==3 .AND. (spin_component==1 .OR. spin_component==2 ) ) THEN
        dos(:) = rho%of_r (:, spin_component)
     ELSE
        isup = 1
        isdw = 2
        dos(:) = rho%of_r (:, isup) + rho%of_r (:, isdw)
     END IF
  end if
  if (lsign) then
     dos(:) = dos(:) * segno(:)
     deallocate(segno)
  endif
  if (iflag == 0) return
  !
  !    symmetrization of the local dos
  !
#ifdef __PARA
  !
  ! reduce charge density across pools
  !
  call mp_sum( dos, inter_pool_comm )
  call psymrho(dos, nrx1, nrx2, nrx3, nr1, nr2, nr3, nsym, s, ftau)
#else
  call symrho (dos, nrx1, nrx2, nrx3, nr1, nr2, nr3, nsym, s, ftau)
#endif
  return

end subroutine local_dos
