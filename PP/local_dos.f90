!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!--------------------------------------------------------------------
subroutine local_dos (iflag, lsign, kpoint, kband, emin, emax, dos)
  !--------------------------------------------------------------------
  !
  !     iflag=0: calculates |psi|^2 for band kband at point kpoint
  !     iflag=1: calculates the local density of state at e_fermi
  !              (only for metals)
  !     iflag=2: calculates the local density of  electronic entropy
  !              (only for metals with fermi spreading)
  !     iflag=3: calculates the integral of local dos from emin to emax
  !
#include "machine.h"
  USE kinds, only: DP
  USE cell_base,            ONLY : omega, tpiba2
  USE basis,                ONLY : nat, ntyp, ityp
  USE ener,                 ONLY : ef
  USE gvect,                ONLY : nr1, nr2, nr3, nrx1, nrx2, nrx3, nrxx, &
                                   ngm, g, ecutwfc
  USE gsmooth,              ONLY : nls, nlsm, nr1s, nr2s, nr3s, &
                                   nrx1s, nrx2s, nrx3s, nrxxs, doublegrid
  USE klist,                ONLY : lgauss, degauss, ngauss, nks, wk, xk
  USE lsda_mod,             ONLY : lsda, nspin, current_spin, isk
  USE scf,                  ONLY : rho
  USE symme,                ONLY : nsym, s, ftau
  USE us,                   ONLY : becsum, nh, nkb, vkb, tvanp
  USE wavefunctions_module, ONLY : evc, psic
  USE wvfct,                ONLY : nbnd, npwx, npw, igk, wg, et, g2kin, &
       gamma_only
  USE io_files,             ONLY:  iunwfc, nwordwfc
#ifdef __PARA
  use mp, only: mp_bcast
#endif
  implicit none
  !
  ! input variables
  !
  integer :: iflag, kpoint, kband
  !
  real(kind=DP) :: emin, emax
  ! output as determined by iflag
  real(kind=DP) :: dos (nrxx)

  logical :: lsign    ! if true and k=gamma and iflag=0 
                      ! write |psi|^2 * sign(psi)
                      
  !
  !    local variables
  !
  integer :: ikb, jkb, ijkb0, ih, jh, na, ijh, np
  ! counters for US PPs
  integer :: ir, is, ig, ibnd, ik, irm, isup, isdw
  ! counters
  real(kind=DP) :: w, w1, modulus, maxmod
  real(kind=DP), allocatable :: rbecp(:,:), segno(:)
  complex(kind=DP), allocatable :: becp(:,:)
  complex(kind=DP) :: phase 
  real(kind=DP), external :: w0gauss, w1gauss
  !
  ! input checks
  !
  if (lsign .and.  (iflag /= 0)) call errore ('local_dos', &
       'inconsistent flags', 1)
  if (lsign .and. (kpoint < 1) .or. (kpoint > nks)) &
       call errore ('local_dos', 'wrong k-point specified', 1)
  if (lsign .and. ( sqrt(xk(1,kpoint)**2 + xk(2,kpoint)**2 + &
                         xk(3,kpoint)**2) > 1d-9 ) ) &
       call errore ('local_dos', 'k must be zero', 1)
  !
  if ( (iflag == 0) .and. ( ( kband < 1 .or. kband > nbnd ) .or. &
        (kpoint < 1) .or. (kpoint > nks)) ) &
       call errore ('local_dos', 'wrong kpoint or band specified', 1)
  !
  if (gamma_only) then 
     allocate (rbecp(nkb,nbnd))
  else 
     allocate (becp(nkb,nbnd))
  end if
  rho(:,:) = 0.d0
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
           call pw_gemm( 'Y', nkb, nbnd, npw, vkb, npwx, evc, npwx, rbecp, nkb)
        else
           call ccalbec (nkb, npwx, npw, nbnd, becp, vkb, evc)
        end if
     !
     !     here we compute the density of states
     !
        do ibnd = 1, nbnd
           if (ibnd == kband .or. iflag /= 0) then
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
              w1 = wg (ibnd, ik) / omega
!
!  Compute and save the sign of the wavefunction at the gamma point
!
              if (lsign) then
                 if (gamma_only) then
                    !  psi(r) is real by construction
                    segno(1:nrxxs) = real(psic(1:nrxxs))
                 else
                    !  determine the phase factor that makes psi(r) real.
                    maxmod=0.d0
                    do ir = 1, nrxxs
                       modulus=abs(psic(ir))
                       if (modulus > maxmod) then
                          irm=ir
                          maxmod=modulus
                       endif
                    enddo
                    if (maxmod > 1.d-10) then
                       phase = psic(irm)/maxmod
                    else
                       call errore('local_dos','zero wavefunction',1)
                    endif
#ifdef __PARA
                    call mp_bcast(phase,0)
#endif
                    segno(1:nrxxs) = real( psic(1:nrxxs)*CONJG(phase) )
                 endif
                 if (doublegrid) call interpolate (segno, segno, 1)
                 segno(:) = sign( 1.d0, segno(:) )
              endif
              !
              do ir=1,nrxxs
                 rho (ir, current_spin) = rho (ir, current_spin) + &
                   w1 * (real ( psic (ir) ) **2 + DIMAG (psic (ir) ) **2)
              enddo
        !
        !    If we have a US pseudopotential we compute here the sumbec term
        !

              w1 = wg (ibnd, ik)
              ijkb0 = 0
              do np = 1, ntyp
                 if (tvanp (np) ) then
                    do na = 1, nat
                       if (ityp (na) == np) then
                          ijh = 1
                          do ih = 1, nh (np)
                             ikb = ijkb0 + ih
                             if(gamma_only) then
                                becsum(ijh,na,current_spin) = &
                                     becsum(ijh,na,current_spin) + w1 * &
                                     rbecp(ikb,ibnd)*rbecp(ikb,ibnd)
                             else
                                becsum(ijh,na,current_spin) = &
                                     becsum(ijh,na,current_spin) + w1 * &
                                     real(conjg(becp(ikb,ibnd))*becp(ikb,ibnd))
                             end if
                             ijh = ijh + 1
                             do jh = ih + 1, nh (np)
                                jkb = ijkb0 + jh
                                if(gamma_only) then
                                   becsum(ijh,na,current_spin) = &
                                     becsum(ijh,na,current_spin) + 2.d0*w1 * &
                                     rbecp(ikb,ibnd)*rbecp(jkb,ibnd)
                                else
                                   becsum(ijh,na,current_spin) = &
                                     becsum(ijh,na,current_spin) + 2.d0*w1 * &
                                     real(conjg(becp(ikb,ibnd))*becp(jkb,ibnd))
                                end if
                                ijh = ijh + 1
                             enddo
                          enddo
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
     deallocate(becp)
  end if
  if (doublegrid) then
     do is = 1, nspin
        call interpolate (rho (1, is), rho (1, is), 1)
     enddo
  endif
  !
  !    Here we add the US contribution to the charge
  !
  call addusdens
  !
  if (nspin == 1) then
     is = 1 
     dos(:) = rho (:, is)
  else
     isup = 1
     isdw = 2
     dos(:) = rho (:, isup) + rho (:, isdw)
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
  call poolreduce (nrxx, dos)
  call psymrho(dos, nrx1, nrx2, nrx3, nr1, nr2, nr3, nsym, s, ftau)
#else
  call symrho (dos, nrx1, nrx2, nrx3, nr1, nr2, nr3, nsym, s, ftau)
#endif
  return

end subroutine local_dos
