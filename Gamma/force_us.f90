!
! Copyright (C) 2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------
subroutine force_us (forcenl)
  !----------------------------------------------------------------------
  !
  ! nonlocal potential contribution to forces
  !
#include "machine.h"

  use pwcom
  USE wavefunctions,  ONLY: evc
  implicit none
  !
  real(kind=DP) :: forcenl (3, nat)
  ! output: the nonlocal contribution

  real(kind=DP), allocatable :: becp(:,:), dbecp (:,:,:)
  ! auxiliary variables contain <beta|psi> and <dbeta|psi>
  complex(kind=DP), allocatable :: vkb1 (:,:)
  ! auxiliary variable contains g*|beta>
  real(kind=DP) :: ps
  integer :: ik, ipol, ibnd, ig, ih, jh, na, nt, ikb, jkb, ijkb0
  ! counters
  !
  forcenl(:,:) = 0.d0
  allocate (becp(nkb, nbnd), dbecp(nkb, nbnd, 3))    
  allocate (vkb1(  npwx, nkb))    
  if (nks.gt.1) rewind iunigk
  !
  !   the forces are a sum over the K points and the bands
  !
  do ik = 1, nks
     if (lsda) current_spin = isk (ik)
     if (nks.gt.1) then
        read (iunigk) npw, igk
        call davcio (evc, nwordwfc, iunwfc, ik, - 1)
        call init_us_2 (npw, igk, xk (1, ik), vkb)
     endif
     !
     call pw_gemm ('Y', nkb, nbnd, npw, vkb, npwx, evc, npwx, becp, nkb)
     !
     do ipol = 1, 3
        do jkb = 1,nkb
           do ig = 1, npw
              vkb1 (ig, jkb) = vkb(ig,jkb) *(0.d0,-1.d0)*g(ipol,igk(ig) )
           enddo
        enddo
        !
        call pw_gemm ('Y', nkb, nbnd, npw, vkb1, npwx, evc, npwx, &
             dbecp(1,1,ipol), nkb)
        !
     enddo
     !
     ijkb0 = 0
     do nt = 1, ntyp
        do na = 1, nat
           if (ityp (na) .eq.nt) then
              do ih = 1, nh (nt)
                 ikb = ijkb0 + ih
                 do ibnd = 1, nbnd
                    ps = deeq (ih, ih, na, current_spin) - &
                         et (ibnd, ik) * qq (ih,ih, nt)
                    do ipol = 1, 3
                       forcenl (ipol, na) = forcenl (ipol, na) - &
                            ps * wg (ibnd, ik) * 2.d0 * tpiba * &
                            dbecp (ikb, ibnd, ipol) * becp (ikb, ibnd)
                    enddo
                 enddo

                 if (tvanp (nt) .or.newpseudo (nt) ) then
                    !
                    ! in US case there is a contribution for jh<>ih. We use
                    ! here the symmetry in the interchange of ih and jh
                    !
                    do jh = ih + 1, nh (nt)
                       jkb = ijkb0 + jh
                       do ibnd = 1, nbnd
                          ps = deeq (ih, jh, na, current_spin) - &
                               et (ibnd, ik) * qq (ih, jh, nt)
                          do ipol = 1, 3
                             forcenl (ipol, na) = forcenl (ipol, na) - &
                                  ps * wg (ibnd, ik) * 2.d0 *  tpiba * &
                                  ( dbecp(ikb, ibnd, ipol)*becp(jkb, ibnd) + &
                                    dbecp(jkb, ibnd, ipol)*becp(ikb, ibnd) )
                          enddo
                       enddo
                    enddo
                 endif
              enddo
              ijkb0 = ijkb0 + nh (nt)
           endif
        enddo
     enddo
  enddo
  !
  ! The total D matrix depends on the ionic position via the augmentation
  ! part \int V_eff Q dr, the term deriving from the derivative of Q
  ! is added in the routine addusforce
  !
  call addusforce (forcenl)
#ifdef __PARA
  ! collect contributions across pools
  call poolreduce (3 * nat, forcenl)
#endif
  !
  ! Since our summation over k points was only on the irreducible BZ we
  ! have to symmetrize the forces. The symmetry matrices are in the
  ! crystal basis so...
  ! Transform to crystal axis...
  !
  do na = 1, nat
     call trnvect (forcenl (1, na), at, bg, - 1)
  enddo
  !
  ! ...symmetrize...
  !
  call symvect (nat, forcenl, nsym, s, irt)
  !
  ! ... and transform back to cartesian axis
  !
  do na = 1, nat
     call trnvect (forcenl (1, na), at, bg, 1)
  enddo
  deallocate (vkb1)
  deallocate (becp, dbecp)
  return
end subroutine force_us

