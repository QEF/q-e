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
#include "machine.h"
  USE cell_base, ONLY: omega
  USE ions_base, ONLY: nat, ntyp=>nsp, ityp
  USE gvect
  USE gsmooth
  USE lsda_mod, ONLY: current_spin
  USE us, ONLY: becsum, nh, tvanp
  USE wvfct, ONLY: npw, wg, igk
  USE wavefunctions_module,  ONLY: evc, psic
  USE becmod, ONLY: becp
  implicit none
  !
  ! input variables
  !
  integer :: ik, kband
  ! input: the k point
  ! input: the band

  real(kind=DP) :: plan (nr3)
  ! output: the planar average of this state
  !
  !    Additional local variables for Ultrasoft PP's
  !

  integer :: ikb, jkb, ijkb0, ih, jh, na, ijh, np
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
  integer :: ir, is, ig, ibnd
  ! counter on 3D r points
  ! counter on spin polarizations
  ! counter on g vectors
  ! counter on bands

  real(kind=DP) :: w, w1
  ! the weight of one k point
  real(kind=DP), allocatable :: aux (:)
  ! auxiliary for rho

  complex(kind=DP), allocatable :: prho (:)
  ! complex charge for fft

  allocate (prho(nrxx))    
  allocate (aux(nrxx))    

  aux(:) = 0.d0
  becsum(:,:,:) = 0.d0

  wg (kband, ik) = 1.d0
  !
  !
  !     First compute the square modulus of the state kband,ik on the smoo
  !     mesh
  !

  psic(1:nrxxs) = (0.d0,0.d0)
  do ig = 1, npw
     psic (nls (igk (ig) ) ) = evc (ig, kband)
  enddo
  call cft3s (psic, nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, 2)

  w1 = wg (kband, ik) / omega
  do ir = 1, nrxxs
     aux(ir) = aux(ir) + w1 * (real(psic(ir))**2 + DIMAG(psic(ir))**2)
  enddo
  !
  !    If we have a US pseudopotential we compute here the sumbec term
  !
  ibnd = kband

  w1 = wg (ibnd, ik)
  ijkb0 = 0
  do np = 1, ntyp
     if (tvanp (np) ) then
        do na = 1, nat
           if (ityp (na) == np) then
              ijh = 1
              do ih = 1, nh (np)
                 ikb = ijkb0 + ih
                 becsum(ijh,na,current_spin) = &
                        becsum(ijh,na,current_spin) + w1 * &
                        real ( conjg(becp(ikb,ibnd)) * becp(ikb,ibnd) )
                 ijh = ijh + 1
                 do jh = ih + 1, nh (np)
                    jkb = ijkb0 + jh
                    becsum(ijh,na,current_spin) = &
                           becsum(ijh,na,current_spin) + w1 * 2.d0 * &
                           real( conjg(becp(ikb,ibnd)) * becp(jkb,ibnd) )
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
     prho (ir) = cmplx (aux (ir), 0.d0)
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
  !
  return
end subroutine local_dos1d
