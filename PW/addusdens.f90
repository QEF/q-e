!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------
subroutine addusdens
  !----------------------------------------------------------------------
  !
  !  This routine adds to the charge density the part which is due to
  !  the US augmentation.
  !
#include "machine.h"
  use pwcom
  implicit none
  !
  !     here the local variables
  !

  integer :: ig, na, nt, ih, jh, ijh, ir, is
  ! counter on G vectors
  ! counter on atoms
  ! the atom type
  ! counter on beta functions
  ! counter on beta functions
  ! composite index ih jh
  ! counter on mesh points
  ! counter on spin polarization

  real(kind=DP), allocatable :: qmod (:), ylmk0 (:,:)
  ! the modulus of G
  ! the spherical harmonics

  complex(kind=DP) :: skk
  complex(kind=DP), allocatable ::  qg (:), aux (:,:)
  ! work space for FFT
  ! work space for rho(G,nspin)


  call start_clock ('addusdens')

  allocate (aux ( ngm, nspin))    
  allocate (qg( nrxx))    
  allocate (qmod( ngm))    
  allocate (ylmk0( ngm, lqx * lqx))    

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
                          skk = eigts1 (ig1 (ig), na) * &
                                eigts2 (ig2 (ig), na) * &
                                eigts3 (ig3 (ig), na)
                          aux(ig,is)=aux(ig,is) + qgm(ig)*skk*becsum(ijh,na,is)
                       enddo
                    enddo
                 endif
              enddo
           enddo
        enddo
     endif
  enddo
  !
  !     convert aux to real space and add to the charge density
  !
  if (okvan) then
     do is = 1, nspin
        call setv (2 * nrxx, 0.d0, qg, 1)
        do ig = 1, ngm
           qg (nl (ig) ) = aux (ig, is)
        enddo

        call cft3 (qg, nr1, nr2, nr3, nrx1, nrx2, nrx3, 1)
        do ir = 1, nrxx
           rho (ir, is) = rho (ir, is) + DREAL (qg (ir) )
        enddo
     enddo
  endif

  deallocate (ylmk0)
  deallocate (qmod)
  deallocate (qg)
  deallocate (aux)

  call stop_clock ('addusdens')
  return
end subroutine addusdens

