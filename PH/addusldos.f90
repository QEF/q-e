!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------
subroutine addusldos (ldos, becsum1)
  !----------------------------------------------------------------------
  !
  !  This routine adds to the local DOS the part which is due to
  !  the US augmentation.
  !
#include "machine.h"
  use pwcom
  implicit none
  complex(kind=DP) :: ldos (nrxx, nspin)
  ! local density of states

  real(kind=DP) :: becsum1 ( (nhm * (nhm + 1) ) / 2, nat, nspin)
  ! input: the becsum1 ter
  !
  !     here the local variables
  !

  integer :: ig, na, nt, ih, jh, ijh, is, ir
  ! counter on G vectors
  ! counter on atoms
  ! counter on atomic types
  ! counter on beta functions
  ! counter on beta functions
  ! composite index ih jh
  ! counter on spin
  ! counter on mesh points

  real(kind=DP),pointer :: ylmk0 (:,:), qmod (:)
  ! the spherical harmonics
  ! the modulus of G

  complex(kind=DP), allocatable :: qg (:), aux (:,:)
  ! auxiliary variable for FFT
  ! auxiliary variable for rho(G)

  allocate (aux ( ngm , nspin))    
  allocate (qg  ( nrxx))    
  allocate (qmod( ngm))    
  allocate (ylmk0 ( ngm , lqx * lqx))    

  call setv (2 * ngm * nspin, 0.d0, aux, 1)
  call ylmr2 (lqx * lqx, ngm, g, gg, ylmk0)
  do ig = 1, ngm
     qmod (ig) = sqrt (gg (ig) )

  enddo
  do nt = 1, ntyp
     if (tvanp (nt) ) then
        ijh = 0
        do ih = 1, nh (nt)
           do jh = ih, nh (nt)
              call qvan2 (ngm, ih, jh, nt, qmod, qgm, ylmk0)
              ijh = ijh + 1
              do na = 1, nat
                 if (ityp (na) .eq.nt) then
                    !
                    !  Multiply becsum and qg with the correct structure factor
                    !
                    do is = 1, nspin
                       do ig = 1, ngm
                          aux (ig, is) = aux (ig, is) + qgm (ig) * becsum1 (ijh, na, &
                               is) * (eigts1 (ig1 (ig), na) * eigts2 (ig2 (ig), na) &
                               * eigts3 (ig3 (ig), na) )
                       enddo
                    enddo
                 endif
              enddo
           enddo
        enddo
     endif
  enddo
  !
  !     convert aux to real space and adds to the charge density
  !
  if (okvan) then
     do is = 1, nspin
        call setv (2 * nrxx, 0.d0, qg, 1)
        do ig = 1, ngm
           qg (nl (ig) ) = aux (ig, is)

        enddo

        call cft3 (qg, nr1, nr2, nr3, nrx1, nrx2, nrx3, 1)
        do ir = 1, nrxx
           ldos (ir, is) = ldos (ir, is) + DREAL (qg (ir) )
        enddo

     enddo

  endif
  deallocate (ylmk0)
  deallocate (qmod)
  deallocate (qg)
  deallocate (aux)
  return
end subroutine addusldos
