!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

!-----------------------------------------------------------------------
subroutine addusdynmat (dynwrk)
  !-----------------------------------------------------------------------
  !
  !     This routine computes the additional terms which are contained in
  !     <psi|V"|psi> part of the dynamical matrix and which are due
  !     to the change of the self consistent D term in the pseudopotential
  !     There are four additional terms which we compute here.
  !
#include "machine.h"


  use pwcom
  USE kinds, only : DP
  USE uspp_param, only: tvanp, nh
  use phcom
  implicit none

  complex(kind=DP) :: dynwrk (3 * nat, 3 * nat)
  ! inp/out: the dynamical matrix

  integer :: ipol, jpol, np, na, nb, nu_i, nu_j, ih, jh, ijh, dim, &
       is
  ! counter on polarizations
  ! counter on pseudopotentials
  ! counter on atoms
  ! counter on modes
  ! counter on solid beta functions
  ! composed dimension of the beta
  ! counter on spin

  complex(kind=DP) :: term (3, 3), dyn1 (3 * nat, 3 * nat)
  ! auxiliary space
  ! auxiliary dynamical matrix


  if (.not.okvan) return
  call start_clock ('addusdynmat')

  dyn1 (:,:) = (0.d0, 0.d0)
  !
  !  We compute the four terms required
  !
  do na = 1, nat
     np = ityp (na)
     if (tvanp (np) ) then
        dim = (nh (np) * (nh (np) + 1) ) / 2
        do ipol = 1, 3
           nu_i = 3 * (na - 1) + ipol
           do jpol = 1, 3
              nu_j = 3 * (na - 1) + jpol
              do is = 1, nspin
                 do ijh = 1, dim
                    dynwrk(nu_i, nu_j)=dynwrk(nu_i, nu_j)+ &
                           int4(ijh,ipol,jpol,na,is) * becsum(ijh,na,is)
                 enddo
              enddo
           enddo
        enddo
        !
        !   The second term requires an exchange of the components.
        !
        term (:,:) = (0.d0, 0.d0)
        do ipol = 1, 3
           do jpol = 1, 3
              ijh = 0
              do ih = 1, nh (np)
                 do jh = ih, nh (np)
                    ijh = ijh + 1
                    do is = 1, nspin
                       term(ipol,jpol) = term(ipol,jpol) + &
                       conjg(int1(ih,jh,ipol,na,is))*alphasum(ijh,jpol,na,is)
                    enddo
                 enddo
              enddo
           enddo
        enddo
        !
        !  And then we add the appropriate terms to the dynamical matrix
        !
        do ipol = 1, 3
           nu_i = 3 * (na - 1) + ipol
           do jpol = 1, 3
              nu_j = 3 * (na - 1) + jpol
              dynwrk (nu_i, nu_j) = dynwrk (nu_i, nu_j) + &
                                    term (ipol, jpol) + term (jpol, ipol)
           enddo
        enddo
        !
        !   the other two terms do not contain a delta ss'
        !
        do nb = 1, nat
           do ipol = 1, 3
              nu_i = 3 * (nb - 1) + ipol
              do jpol = 1, 3
                 nu_j = 3 * (na - 1) + jpol
                 ijh = 0
                 do ih = 1, nh (np)
                    do jh = ih, nh (np)
                       ijh = ijh + 1
                       do is = 1, nspin
                          dyn1(nu_i,nu_j)=dyn1(nu_i,nu_j) + &
                                          conjg(int2(ih,jh,ipol,nb,na)) * &
                                               alphasum(ijh,jpol,na,is) + &
                                          int5(ijh,ipol,jpol,nb,na) *     &
                                                  becsum(ijh,na,is)
                       enddo
                    enddo
                 enddo
              enddo
           enddo
        enddo
     endif

  enddo
  do nu_i = 1, nmodes
     do nu_j = 1, nmodes
        dynwrk (nu_i, nu_j) = dynwrk (nu_i, nu_j) + &
                              dyn1 (nu_i, nu_j) + conjg (dyn1 (nu_j, nu_i) )
     enddo

  enddo
  deallocate (int4)
  deallocate (int5)

  call stop_clock ('addusdynmat')
  return
end subroutine addusdynmat
