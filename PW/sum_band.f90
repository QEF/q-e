!
! Copyright (C) 2001-2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!--------------------------------------------------------------------
subroutine sum_band
  !--------------------------------------------------------------------
  !
  !     calculates the symmetrized charge density and sum of occupied
  !     eigenvalues.
  !     this version works also for metals (gaussian spreading technique)
  !
  !
#include "machine.h"
  use pwcom
  use io_files, only: iunwfc, nwordwfc, iunigk
  USE wavefunctions_module,    ONLY : evc, psic
  use becmod
#ifdef __PARA
  use para
#endif
  implicit none
  !
  !    local variables
  !
  integer :: ikb, jkb, ijkb0, ih, jh, ijh, na, np
  ! counters on beta functions, atoms, pseudopotentials
  integer :: ir, is, ig, ibnd, ik
  ! counter on 3D r points
  ! counter on spin polarizations
  ! counter on g vectors
  ! counter on bands
  ! counter on k points
  real(kind=DP) :: w1
  ! weight
  !
  call start_clock ('sum_band')
  becsum(:,:,:) = 0.d0
  rho(:,:) = 0.d0
  eband = 0.d0
  demet = 0.d0
  !
  !   calculate weights for the insulator case
  !
  if (.not.lgauss.and..not.ltetra.and..not.tfixed_occ) then
     call iweights (nks, wk, nbnd, nelec, wg)
     !
     !    calculate weights for the metallic case
     !
  elseif (ltetra) then
#ifdef __PARA
     call poolrecover (et, nbnd, nkstot, nks)
     if (me.eq.1.and.mypool.eq.1) then
#endif
        call tweights (nkstot, nspin, nbnd, nelec, ntetra, tetra, et, ef, wg)
#ifdef __PARA
     endif
     call poolscatter (nbnd, nkstot, wg, nks, wg)
     if (me.eq.1) call poolbcast (1, ef)
     call broadcast (1, ef)
#endif
  elseif (lgauss) then
     call gweights (nks, wk, nbnd, nelec, degauss, ngauss, et, ef, demet, wg)
  elseif (tfixed_occ) then
     do is=1,nspin
        do ibnd=1,nbnd
           wg(ibnd,is)=f_inp(ibnd,is)
        enddo
     enddo
  endif
  !
  ! Needed for LDA+U
  !
  if (lda_plus_u) call new_ns
  !
  !     here we sum for each k point the contribution
  !     of the wavefunctions to the charge
  !
  if (nks.gt.1) rewind (iunigk)
  do ik = 1, nks
     if (lsda) current_spin = isk (ik)
     if (nks.gt.1) then
        read (iunigk) npw, igk
        call davcio (evc, nwordwfc, iunwfc, ik, - 1)
     endif
     call init_us_2 (npw, igk, xk (1, ik), vkb)
     !
     !     here we compute the band energy: the sum of the eigenvalues
     !
     do ibnd = 1, nbnd
        eband = eband+et (ibnd, ik) * wg (ibnd, ik)
        !            WRITE( stdout,'(4x, " ibnd = ", i5, " ik= ",i5,
        !     +                       f15.5)') ibnd, ik,
        !     +                           wg(ibnd,ik)
        !
        !  the sum of eband and demet is the integral for e < ef of e n(e)
        !  which reduces for degauss=0 to the sum of the eigenvalues.
        !  the factors two is for spin degeneracy
        !
        psic(:) = (0.d0,0.d0)
        do ig = 1, npw
           psic (nls (igk (ig) ) ) = evc (ig, ibnd)
        enddo
        call cft3s (psic, nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, 2)
        w1 = wg (ibnd, ik) / omega
        !
        ! increment the charge density ...
        !
        do ir = 1, nrxxs
           rho (ir, current_spin) = rho (ir, current_spin) + &
                w1 * (DREAL (psic (ir) ) **2 + DIMAG (psic (ir) ) **2)
        enddo
     end do
     !
     !    If we have a US pseudopotential we compute here the sumbec term
     !
     if (.not.okvan) go to 10
     !
     call ccalbec (nkb, npwx, npw, nbnd, becp, vkb, evc)
     !
     call start_clock ('sumbec')
     do ibnd = 1, nbnd
        w1 = wg (ibnd, ik)
        ijkb0 = 0
        do np = 1, ntyp
           if (tvanp (np) ) then
              do na = 1, nat
                 if (ityp (na) .eq.np) then
                    ijh = 1
                    do ih = 1, nh (np)
                       ikb = ijkb0 + ih
                       becsum (ijh, na, current_spin) = &
                            becsum (ijh, na, current_spin) + &
                            w1 * DREAL(conjg(becp (ikb, ibnd))*becp(ikb,ibnd))
                       ijh = ijh + 1
                       do jh = ih + 1, nh (np)
                          jkb = ijkb0 + jh
                          becsum (ijh, na, current_spin) = &
                               becsum (ijh, na, current_spin) + w1 * 2.d0 * &
                               DREAL(conjg(becp(ikb,ibnd))*becp(jkb,ibnd) )
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
     enddo
     call stop_clock ('sumbec')
10   continue
  enddo
  !
  ! If a double grid is used, interpolate onto the fine grid
  !
  if (doublegrid) then
     do is = 1, nspin
        call interpolate (rho (1, is), rho (1, is), 1)
     enddo
  endif
  !
  !    Here we add the Ultrasoft contribution to the charge
  !
  if (okvan) call addusdens
#ifdef __PARA
  call poolreduce (1, eband)
  call poolreduce (1, demet)
#endif
  !
  !    symmetrization of the charge density (and local magnetization)
  !
#ifdef __PARA
  !
  ! reduce charge density across pools
  !
  call poolreduce (nspin * nrxx, rho)
  do is = 1, nspin
     call psymrho (rho (1, is), nrx1, nrx2, nrx3, nr1, nr2, nr3, nsym, &
          s, ftau)
  enddo
#else
  do is = 1, nspin
     call symrho (rho (1, is), nrx1, nrx2, nrx3, nr1, nr2, nr3, nsym, &
          s, ftau)
  enddo
#endif

  call stop_clock ('sum_band')
  return
end subroutine sum_band

