!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
function cgracsc (nkb, bec1, bec2, nhm, ntyp, nh, qq, nat, ityp, &
     npw, psi1, psi2, tvanp)
  !-----------------------------------------------------------------------
  !
  !     This function computes the scalar product between two wavefunction
  !     and the S matrix of the US pseudopotential: <psi1 | S | psi2 >.
  !     It assumes that the product of psi1 with all the beta functions
  !     is in bec1, and the product of psi2 is in bec2.
  !
  !
#include "machine.h"
  USE kinds
  implicit none
  !
  !     here the dummy variables
  !

  integer :: nkb, npw, nhm, ntyp, nat, ityp (nat), nh (ntyp)
  ! input: the number of beta functions
  ! input: the number of plane waves
  ! input: the maximum number of solid be
  ! input: the number of types of atoms
  ! input: the number of atoms
  ! input: the type of each atom
  ! input: the number of beta for each ty

  complex(kind=DP) :: bec1 (nkb), bec2 (nkb), psi1 (npw), psi2 (npw), &
       cgracsc
  ! input: the product of beta and psi1
  ! input: the product of beta and psi2
  ! input: the first wavefunction
  ! input: the second wavefunction
  ! output: the value of the scalar produ

  real(kind=DP) :: qq (nhm, nhm, ntyp)
  ! input: the q values defining S
  logical :: tvanp (ntyp)
  ! input: if true the pseudo is vanderb
  !
  !    Here the local variables
  !

  integer :: ikb, jkb, na, np, ijkb0, ih, jh
  ! counter on total beta functions
  ! counter on total beta functions
  ! counter on atoms
  ! the pseudopotential of each atom
  ! auxiliary variable to compute ikb and jkb
  ! counter on solid beta functions
  ! counter on solid beta functions

  complex(kind=DP) :: scal, ZDOTC
  !
  scal = ZDOTC (npw, psi1, 1, psi2, 1)
#ifdef __PARA
  call reduce (2, scal)
#endif
  ijkb0 = 0
  do np = 1, ntyp
     if (tvanp (np) ) then
        do na = 1, nat
           if (ityp (na) .eq.np) then
              do ih = 1, nh (np)
                 ikb = ijkb0 + ih
                 do jh = 1, nh (np)
                    jkb = ijkb0 + jh
                    scal = scal + qq (ih,jh,np)*conjg(bec1(ikb))*bec2(jkb)
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

  cgracsc = scal
  return
end function cgracsc

