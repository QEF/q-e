!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
subroutine drhodvnl (ik, ikk, nper, nu_i0, wdyn, dbecq, dalpq)
  !-----------------------------------------------------------------------
  !
  !  This routine compute the term of the dynamical matrix due to
  !  the orthogonality constraint. Only the part which is due to
  !  the nonlocal terms is computed here
  !
#include "machine.h"


  use pwcom
  use parameters, only : DP
  use phcom
  implicit none
  integer :: ik, ikk, nper, nu_i0
  ! input: the current k point
  ! input: the number of perturbations
  ! input: the initial mode

  complex(kind=DP) :: dbecq (nkb, nbnd, nper), dalpq (nkb, nbnd,3, nper),&
          wdyn (3 * nat, 3 * nat)
  ! input: the becp with psi_{k+
  ! input: the alphap with psi_{
  ! output: the term of the dynamical

  complex(kind=DP) :: ps, dynwrk (3 * nat, 3 * nat)
  ! dynamical matrix
  complex(kind=DP) , allocatable :: ps1 (:,:), ps2 (:,:,:)

  integer :: ibnd, ijkb0, ijkb0b, ih, jh, ikb, jkb, ipol, jpol, &
       startb, lastb, iper, na, nb, nt, ntb, mu, nu
  ! counter on bands
  ! counter on beta functions
  ! counter on beta functions
  ! counter on polarization
  ! divide band among processors
  ! counter on perturbations
  ! counter on atoms and types
  ! counter on the weights

  allocate (ps1 (  nkb , nbnd))    
  allocate (ps2 (  nkb , nbnd , 3))    

  call setv (2 * 3 * nat * 3 * nat, 0.d0, dynwrk, 1)

  call divide (nbnd, startb, lastb)
  call setv (2 * nkb * nbnd, 0.d0, ps1, 1)
  call setv (2 * nkb * nbnd * 3, 0.d0, ps2, 1)
  !
  !   Here we prepare the two terms
  !
  ijkb0 = 0
  do nt = 1, ntyp
     do na = 1, nat
        if (ityp (na) .eq.nt) then
           do ih = 1, nh (nt)
              ikb = ijkb0 + ih
              do jh = 1, nh (nt)
                 jkb = ijkb0 + jh
                 do ibnd = startb, lastb
                    ps1 (ikb, ibnd) = ps1 (ikb, ibnd) + &
                         (deeq (ih, jh, na,current_spin) - &
                          et (ibnd, ikk) * qq (ih, jh, nt) )*becp1(jkb,ibnd,ik)
                    do ipol = 1, 3
                       ps2 (ikb, ibnd, ipol) = ps2 (ikb, ibnd, ipol) + &
                            (deeq (ih, jh,na, current_spin) - &
                             et (ibnd, ikk) * qq (ih, jh, nt) ) * &
                            alphap (jkb, ibnd, ipol, ik)
                       if (okvan) ps2 (ikb, ibnd, ipol) = &
                            ps2 (ikb, ibnd, ipol) + &
                               int1 (ih, jh, ipol, na, current_spin) * &
                               becp1 (jkb, ibnd, ik)
                    enddo
                 enddo
              enddo
           enddo
           ijkb0 = ijkb0 + nh (nt)
        endif
     enddo
  enddo
  !
  !     Here starts the loop on the atoms (rows)
  !
  ijkb0 = 0
  do nt = 1, ntyp
     do na = 1, nat
        if (ityp (na) .eq.nt) then
           do ipol = 1, 3
              mu = 3 * (na - 1) + ipol
              do ibnd = startb, lastb
                 do ih = 1, nh (nt)
                    ikb = ijkb0 + ih
                    do iper = 1, nper
                       nu = nu_i0 + iper
                       dynwrk (nu, mu) = dynwrk (nu, mu) + &
                            2.d0 * wk (ikk) * (ps2 (ikb, ibnd, ipol) * &
                            conjg (dbecq (ikb, ibnd, iper) ) + &
                            ps1(ikb,ibnd) * conjg (dalpq(ikb,ibnd,ipol,iper)) )
                    enddo
                 enddo
                 if (okvan) then
                    ijkb0b = 0
                    do ntb = 1, ntyp
                       do nb = 1, nat
                          if (ityp (nb) .eq.ntb) then
                             do ih = 1, nh (ntb)
                                ikb = ijkb0b + ih
                                ps = (0.d0, 0.d0)
                                do jh = 1, nh (ntb)
                                   jkb = ijkb0b + jh
                                   ps = ps + int2 (ih, jh, ipol, na, nb) * &
                                        becp1 (jkb, ibnd,ik)
                                enddo
                                do iper = 1, nper
                                   nu = nu_i0 + iper
                                   dynwrk (nu, mu) = dynwrk (nu, mu) + &
                                        2.d0 * wk (ikk) * ps * &
                                        conjg (dbecq (ikb, ibnd, iper) )
                                enddo
                             enddo
                             ijkb0b = ijkb0b + nh (ntb)
                          endif
                       enddo
                    enddo
                 endif
              enddo
           enddo
           ijkb0 = ijkb0 + nh (nt)
        endif
     enddo
  enddo
#ifdef __PARA
  call reduce (2 * 3 * nat * 3 * nat, dynwrk)
#endif
  call DAXPY (18 * nat * nat, 1.d0, dynwrk, 1, wdyn, 1)

  deallocate (ps2)
  deallocate (ps1)
  return
end subroutine drhodvnl
