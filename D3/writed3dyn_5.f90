!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
subroutine writed3dyn_5 (d3dyn_x, filename, isw)
  !-----------------------------------------------------------------------
  !
  !     writes in a file the third derivative of dynamical matrix
  !     isw = +1  :  d3dyn_x is in cartesian axis
  !     isw = -1  :  rotates d3dyn_x from the basis of pattern to
  !                      cartesian axis
#include "machine.h"
  USE kinds, only : DP
  use pwcom
  use phcom
  use d3com
#ifdef __PARA
  use para
#endif
  implicit none
  integer :: isw, iud3dyn, n_d3, na, nb, icart, jcart, kcart, na_i, &
       na_j, na_k
  ! input: switch
  ! index on cartesian coordinates
  ! index on cartesian coordinates
  ! index on cartesian coordinates
  ! index on modes
  ! index on modes
  ! index on modes

  complex (kind = dp) :: d3dyn_x (3 * nat, 3 * nat, 3 * nat), work
  ! input: the third derivative of the dynamical matrix
  complex (kind = dp), allocatable :: aux (:,:,:)
  ! auxiliary space

  character (len=*) :: filename
  ! input: the name of the file

#ifdef __PARA
  if (me.ne.1.or.mypool.ne.1) return
#endif

  allocate  (aux( 3 * nat, 3 * nat, 3 * nat))    
  if (isw.eq. + 1) then
     call ZCOPY (27 * nat * nat * nat, d3dyn_x, 1, aux, 1)
  elseif (isw.eq. - 1) then
     !
     !   Rotates third derivative of the dynamical basis from the basis
     !   of modes to cartesisn axis
     !
     do kcart = 1, 3 * nat
        do icart = 1, 3 * nat
           do jcart = 1, 3 * nat
              work = (0.d0, 0.d0)
              do na_k = 1, 3 * nat
                 do na_i = 1, 3 * nat
                    do na_j = 1, 3 * nat
                       work = work + conjg (ug0 (kcart, na_k) ) * u (icart, na_i) &
                            * d3dyn_x (na_k, na_i, na_j) * conjg (u (jcart, na_j) )
                    enddo
                 enddo
              enddo
              aux (kcart, icart, jcart) = work
           enddo
        enddo

     enddo

  endif
  iud3dyn = 57

  open (unit = iud3dyn, file = trim(filename), status = 'unknown')
  do n_d3 = 1, 3 * nat
     write (iud3dyn, * )
     write (iud3dyn,  * ) '               modo:', n_d3
     write (iud3dyn, * )
     do na = 1, nat
        do nb = 1, nat
           write (iud3dyn, '(2i3)') na, nb
           do icart = 1, 3
              write (iud3dyn, '(3E24.12)') (aux (n_d3, icart + 3 * (na - 1) , &
                   jcart + 3 * (nb - 1) ) , jcart = 1, 3)
           enddo
        enddo
     enddo

  enddo

  close (iud3dyn)

  deallocate (aux)
  return
end subroutine writed3dyn_5
