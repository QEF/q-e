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
  use parameters, only: DP
  use pwcom
  USE io_files, ONLY: iunwfc, nwordwfc
  USE wavefunctions_module,  ONLY: evc, psic
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
  integer :: ir, is, ig, ibnd, ik, irm
  ! counters
  real(kind=DP) :: w0gauss, w1gauss, w, w1
  real(kind=DP) :: seno, coseno, modulus, maxmod
  complex(kind=DP), allocatable :: segno(:), becp(:,:)
  logical :: lgamma
  external w0gauss, w1gauss
  !
  allocate (becp(nkb,nbnd))
  rho(:,:) = 0.d0
  dos(:) = 0.d0
  becsum(:,:,:) = 0.d0
  !
  !   calculate the correct weights
  !
  if (iflag.ne.0.and..not.lgauss) call errore ('local_dos', &
       'gaussian broadening needed', 1)
  if (iflag.eq.2.and.ngauss.ne. - 99) call errore ('local_dos', &
       ' beware: not using Fermi-Dirac function ',  - ngauss)
  do ik = 1, nks
     do ibnd = 1, nbnd
        if (iflag.eq.0) then
           wg (ibnd, ik) = 0.d0
        elseif (iflag.eq.1) then
           wg (ibnd, ik) = wk (ik) * w0gauss ( (ef - et (ibnd, ik) ) &
                / degauss, ngauss) / degauss
        elseif (iflag.eq.2) then
           wg (ibnd, ik) = - wk (ik) * w1gauss ( (ef - et (ibnd, ik) ) &
                / degauss, ngauss)
        elseif (iflag.eq.3) then
           if (et (ibnd, ik) .le. emax .and. et (ibnd, ik) .ge.emin) then
              wg (ibnd, ik) = wk (ik)
           else
              wg (ibnd, ik) = 0.d0
           endif
        else
           call errore ('local_dos', ' iflag not allowed', abs (iflag) )
        endif
     enddo

  enddo

  if (iflag.eq.0) wg (kband, kpoint) = 1.d0
  !
  !     here we sum for each k point the contribution
  !     of the wavefunctions to the density of states
  !
  lgamma=.false.
  if (iflag.eq.0) &
     lgamma=(sqrt(xk(1,kpoint)**2+xk(2,kpoint)**2+xk(3,kpoint)**2).lt.1d-9)
  do ik = 1, nks
     if (ik.eq.kpoint.or.iflag.ne.0) then
        if (lsda) current_spin = isk (ik)
        call gk_sort (xk (1, ik), ngm, g, ecutwfc / tpiba2, npw, igk, g2kin)
        call davcio (evc, nwordwfc, iunwfc, ik, - 1)
        call init_us_2 (npw, igk, xk (1, ik), vkb)

        call ccalbec (nkb, npwx, npw, nbnd, becp, vkb, evc)
     !
     !     here we compute the density of states
     !
        do ibnd = 1, nbnd
           if (ibnd.eq.kband.or.iflag.ne.0) then
              psic(1:nrxxs) = (0.d0,0.d0)
              do ig = 1, npw
                 psic (nls (igk (ig) ) ) = evc (ig, ibnd)
              enddo
              call cft3s (psic, nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, 2)
              w1 = wg (ibnd, ik) / omega
!
!  Compute and save the sign of the wavefunction at the gamma point
!  after multiplication for a phase factor that makes it real.
!
              if (iflag.eq.0.and.lgamma.and.lsign) then
                 if (ik.eq.kpoint.and.ibnd.eq.kband) then
                    allocate(segno(nrxx))
                    if (doublegrid) then
                       call cinterpolate (segno, psic, 1)
                    else
                       call ZCOPY(nrxx,psic,1,segno,1)
                    endif   
                    maxmod=0.d0
                    do ir = 1, nrxx
                       modulus=abs(segno(ir))
                       if (modulus.gt.maxmod) then
                          irm=ir
                          maxmod=modulus
                       endif
                    enddo
                    if (maxmod.gt.1.d-10) then
                       coseno=real(segno(irm))/maxmod
                       seno= DIMAG(segno(irm))/maxmod
                    else   
                       call errore('local_dos','zero wavefuntion',1)
                    endif
#ifdef __PARA
                    call mp_bcast(coseno,0)
                    call mp_bcast(seno,0)
#endif
                    do ir=1,nrxx
                       segno(ir)=segno(ir)*cmplx(coseno,-seno)
                       segno(ir)=cmplx(sign(1.d0,real(segno(ir))),0.d0)
                    enddo
                 endif
              endif

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
                       if (ityp (na) .eq.np) then
                          ijh = 1
                          do ih = 1, nh (np)
                             ikb = ijkb0 + ih
                             becsum(ijh,na,current_spin)=becsum(ijh,na, &
                               current_spin)+w1*real(conjg(becp(ikb,ibnd)) &
                                            * becp (ikb, ibnd) )
                             ijh = ijh + 1
                             do jh = ih + 1, nh (np)
                                jkb = ijkb0 + jh
                                becsum(ijh,na,current_spin)=becsum(ijh,na, &
                                current_spin)+w1*2.d0*real(conjg(becp(ikb,ibnd)) &
                               * becp (jkb, ibnd) )
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
           endif
        enddo
     endif
  enddo
  deallocate(becp)
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
  call DCOPY (nrxx, rho (1, 1), 1, dos, 1)
  do is = 2, nspin
     call DAXPY (nrxx, 1.d0, rho (1, is), 1, dos, 1)
  enddo
  if (iflag.eq.0.and.lgamma.and.lsign) then
     do ir=1,nrxx
        dos(ir)=dos(ir)*segno(ir)
     enddo
     deallocate(segno)
  endif
  if (iflag.eq.0) return
  !
  !    symmetrization of the local dos
  !
#ifdef __PARA
  !
  ! reduce charge density across pools
  !
  call poolreduce (nrxx, dos)
  call psymrho (dos, nrx1, nrx2, nrx3, nr1, nr2, nr3, nsym, s, ftau)
#else
  call symrho (dos, nrx1, nrx2, nrx3, nr1, nr2, nr3, nsym, s, ftau)
#endif
  return

end subroutine local_dos
