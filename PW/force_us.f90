!
! Copyright (C) 2001 PWSCF group
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
  use becmod
  implicit none
  !
  !   the dummy variable
  !
  real(kind=DP) :: forcenl (3, nat)
  ! output: the nonlocal contribution

  complex(kind=DP), allocatable :: dbecp (:,:,:), work1 (:,:), aux1 (:)
  ! auxiliary variable contains <dbeta|psi>
  ! auxiliary variable contains g*psi
  ! auxiliary variable contains beta(k+g)
  complex(kind=DP) ::  fact, ZDOTC
  ! a multiplicative factor
  ! the scalar product function

  real(kind=DP) :: fac, ps

  integer :: ik, ipol, ibnd, ig, ih, jh, na, nt, ikb, jkb, jkb2, &
       ijkb0
  ! counter on k points
  ! counter on polarization
  ! counter on bands
  ! counter on G vectors
  ! counters on atomic beta functions
  ! counter on atoms
  ! counter on type of atoms
  ! counters on beta functions
  !
  forcenl(:,:) = 0.d0
  allocate (dbecp(  nkb,  3, nbnd))    
  allocate (work1(  npwx,  3))    
  allocate (aux1(  npwx))    
  if (nks.gt.1) rewind iunigk
  !
  !   the forces are a sum over the K points and the bands
  !
  fact = DCMPLX (0.d0, tpiba)
  do ik = 1, nks
     if (lsda) current_spin = isk (ik)
     if (nks.gt.1) then
        read (iunigk) npw, igk
        call davcio (evc, nwordwfc, iunwfc, ik, - 1)
        call init_us_2 (npw, igk, xk (1, ik), vkb)
     endif
     call ccalbec (nkb, npwx, npw, nbnd, becp, vkb, evc)
     jkb2 = 0
     do nt = 1, ntyp
        do na = 1, nat
           if (ityp (na) .eq.nt) then
              do ih = 1, nh (nt)
                 jkb2 = jkb2 + 1
                 do ig = 1, npw
                    do ipol = 1, 3
                       work1 (ig, ipol) = vkb(ig,jkb2)* g(ipol,igk(ig) )
                    enddo
                 enddo
                 do ibnd = 1, nbnd
                    do ipol = 1, 3
                       dbecp (jkb2, ipol, ibnd) = fact * &
                            ZDOTC (npw,work1(1, ipol),1,evc(1,ibnd),1)
                    enddo
                 enddo
              enddo
           endif
        enddo

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
                    fac = wg (ibnd, ik)
                    do ipol = 1, 3
                       forcenl (ipol, na) = forcenl (ipol, na) - &
                            ps * fac * 2.d0 * &
                            DREAL (conjg (dbecp (ikb, ipol, ibnd) ) &
                                   * becp (ikb, ibnd) )
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
                          fac = wg (ibnd, ik)
                          ps = deeq (ih, jh, na, current_spin) - &
                               et (ibnd, ik) * qq (ih, jh, nt)
                          do ipol = 1, 3
                             forcenl (ipol, na) = forcenl (ipol, na) - &
                                  ps * fac * 2.d0 * &
                                  DREAL (conjg (dbecp (ikb, ipol, ibnd) ) * &
                                         becp (jkb, ibnd) + &
                                         dbecp (jkb, ipol, ibnd) * &
                                         conjg (becp (ikb, ibnd) ) )
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
#ifdef __PARA
  call reduce (3 * nat, forcenl)
#endif
  !
  ! The total D matrix depends on the ionic position due to the augmentati
  ! part \int V_eff Q dr, the term deriving from the derivative of Q
  ! is added in the routine addusforce
  !
  call addusforce (forcenl)
#ifdef __PARA
  ! collect contributions across poo
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
  deallocate (aux1)
  deallocate (work1)
  deallocate (dbecp)
  return
end subroutine force_us

