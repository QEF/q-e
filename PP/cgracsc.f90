!
! Copyright (C) 2001-2005 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
function cgracsc (nkb, bec1, bec2, nhm, ntyp, nh, qq, nat, ityp, &
     npw, psi1, psi2, upf)
  !-----------------------------------------------------------------------
  !
  !     This function computes the scalar product between two wavefunction
  !     and the S matrix of the US pseudopotential: <psi1 | S | psi2 >.
  !     It assumes that the product of psi1 with all the beta functions
  !     is in bec1, and the product of psi2 is in bec2.
  !
  !
#include "f_defs.h"
  USE kinds
  USE pseudo_types, ONLY : pseudo_upf 
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

  complex(DP) :: bec1 (nkb), bec2 (nkb), psi1 (npw), psi2 (npw), &
       cgracsc
  ! input: the product of beta and psi1
  ! input: the product of beta and psi2
  ! input: the first wavefunction
  ! input: the second wavefunction
  ! output: the value of the scalar produ

  real(DP) :: qq (nhm, nhm, ntyp)
  ! input: the q values defining S
  type(pseudo_upf) :: upf (ntyp)
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

  complex(DP) :: scal, ZDOTC
  !
  scal = ZDOTC (npw, psi1, 1, psi2, 1)
#ifdef __PARA
  call reduce (2, scal)
#endif
  ijkb0 = 0
  do np = 1, ntyp
     if (upf(np)%tvanp ) then
        do na = 1, nat
           if (ityp (na) .eq.np) then
              do ih = 1, nh (np)
                 ikb = ijkb0 + ih
                 do jh = 1, nh (np)
                    jkb = ijkb0 + jh
                    scal = scal + qq (ih,jh,np)*CONJG(bec1(ikb))*bec2(jkb)
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

!
!-----------------------------------------------------------------------
function cgracsc_nc (nkb, bec1, bec2, nhm, ntyp, nh, nat, ityp, &
     npw, npol, psi1, psi2, upf)
  !-----------------------------------------------------------------------
  !
  !     This function computes the scalar product between two wavefunction
  !     and the S matrix of the US pseudopotential: <psi1 | S | psi2 >.
  !     It assumes that the product of psi1 with all the beta functions
  !     is in bec1, and the product of psi2 is in bec2.
  !
  !
#include "f_defs.h"
  USE kinds
  USE uspp, ONLY: qq, qq_so
  USE spin_orb, ONLY: lspinorb
  USE pseudo_types, ONLY : pseudo_upf
  implicit none
  !
  !     here the dummy variables
  !

  integer :: nkb, npw, npol, nhm, ntyp, nat, ityp (nat), nh (ntyp)
  ! input: the number of beta functions
  ! input: the number of plane waves
  ! input: the maximum number of solid be
  ! input: the number of types of atoms
  ! input: the number of atoms
  ! input: the type of each atom
  ! input: the number of beta for each ty

  complex(DP) :: bec1 (nkb,npol), bec2 (nkb,npol), &
                      psi1 (npw,npol), psi2 (npw,npol), cgracsc_nc
  ! input: the product of beta and psi1
  ! input: the product of beta and psi2
  ! input: the first wavefunction
  ! input: the second wavefunction
  ! output: the value of the scalar produ

  type(pseudo_upf) :: upf (ntyp)
  ! input: if true the pseudo is vanderb
  !
  !    Here the local variables
  !

  integer :: ikb, jkb, na, np, ijkb0, ih, jh, ipol, jpol, ijh
  ! counter on total beta functions
  ! counter on total beta functions
  ! counter on atoms
  ! the pseudopotential of each atom
  ! auxiliary variable to compute ikb and jkb
  ! counter on solid beta functions
  ! counter on solid beta functions

  complex(DP) :: scal, ZDOTC
  !
  scal = ZDOTC (npw*npol, psi1, 1, psi2, 1)
#ifdef __PARA
  call reduce (2, scal)
#endif
  ijkb0 = 0
  do np = 1, ntyp
     if (upf(np)%tvanp ) then
        do na = 1, nat
           if (ityp (na) .eq.np) then
              do ih = 1, nh (np)
                 ikb = ijkb0 + ih
                 do jh = 1, nh (np)
                    jkb = ijkb0 + jh
                    if (lspinorb) then
                       ijh=0
                       do ipol=1,npol
                          do jpol=1,npol
                            ijh=ijh+1
                            scal=scal+qq_so(ih,jh,ijh,np)* &
                                      CONJG(bec1(ikb,ipol))*bec2(jkb,jpol)
                          end do
                       end do
                    else
                       do ipol=1,npol
                          scal=scal+qq(ih,jh,np)* &
                               CONJG(bec1(ikb,ipol))*bec2(jkb,ipol)
                       end do
                    end if
                 end do
              end do
              ijkb0 = ijkb0 + nh (np)
           end if
        end do
     else
        do na = 1, nat
           if (ityp (na) .eq.np) ijkb0 = ijkb0 + nh (np)
        enddo
     endif
  enddo

  cgracsc_nc = scal
  return
end function cgracsc_nc
