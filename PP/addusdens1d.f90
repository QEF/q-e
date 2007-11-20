!
! Copyright (C) 2001-2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------
subroutine addusdens1d (plan, prho)
  !----------------------------------------------------------------------
  !
  !  This routine adds to the charge density the part which is due to
  !  the US augmentation. This is done only along the G_z direction in
  !  reciprocal space. The output of the routine is the planar average
  !  of the charge density.
  !
#include "f_defs.h"
  USE kinds, only: DP
  USE cell_base, ONLY: alat, omega, celldm
  USE ions_base, ONLY: nat, ntyp => nsp, ityp
  USE gvect, ONLY: nr3, nrx3, nrxx, nl, eigts1, eigts2, eigts3, ig1,ig2,ig3
  USE lsda_mod, ONLY: current_spin
  USE uspp, ONLY: becsum
  USE uspp_param, ONLY: upf, lmaxq, nh
  !
  !     here the local variables
  !
  implicit none
  integer :: ig, na, nt, ih, jh, ijh, ngm1d, ig1dto3d (nr3), &
       igtongl1d (nr3), nl1d (nr3)
  ! counter on G vectors
  ! counter on atoms
  ! counter on atomic types
  ! counter on beta functions
  ! counter on beta functions
  ! composite index ih jh
  ! the number of 1D G vectors on this processor
  ! correspondence 1D with 3D G vectors
  ! the correspondence 1D with the 3D shells
  ! correspondence 1D FFT mesh G with array G

  real(DP) :: plan (nr3), dimz, g1d (3, nr3), gg1d (nr3), qmod (nr3), &
       qgr (nr3), qgi (nr3), ylmk0 (nr3, lmaxq * lmaxq)
  !  the planar average
  !  dimension along z
  !  ngm1d 3D vectors with the 1D G of this proc
  !  ngm1d scalars with the modulus of 1D G
  ! the modulus of G
  ! real and
  ! imaginary part of qg
  ! the spherical harmonics

  complex(DP) :: skk, prho (nrxx), qg (nrx3)
  ! auxiliary variable
  ! auxiliary space for the charge
  ! auxiliary variable for FFT
  ! auxiliary variable for rho(G,nspin)
  complex(DP), allocatable :: qgm(:), aux (:)


  call ggen1d (ngm1d, g1d, gg1d, ig1dto3d, nl1d, igtongl1d)
  allocate (qgm(ngm1d), aux(ngm1d))
  do ig = 1, ngm1d
     qmod (ig) = sqrt (gg1d (ig) )
  enddo
  aux(:) = (0.d0, 0.d0)

  if (ngm1d > 0) then
     call ylmr2 (lmaxq * lmaxq, ngm1d, g1d, gg1d, ylmk0)
     do nt = 1, ntyp
        if (upf(nt)%tvanp  ) then
           ijh = 0
           do ih = 1, nh (nt)
              do jh = ih, nh (nt)
                 call qvan2 (ngm1d, ih, jh, nt, qmod, qgm, ylmk0)
                 ijh = ijh + 1
                 do na = 1, nat

                    if (ityp (na) == nt) then
                       !
                       !  Multiply becsum and qg with the correct structure factor
                       !
                       do ig = 1, ngm1d

                          skk = eigts1 (ig1 (ig1dto3d (ig) ), na) * &
                                eigts2 (ig2 (ig1dto3d (ig) ), na) * &
                                eigts3 (ig3 (ig1dto3d (ig) ), na)
                          aux (ig) = aux (ig) + qgm (ig) * skk * &
                               becsum (ijh, na, current_spin)
                       enddo
                    endif
                 enddo
              enddo
           enddo
        endif
     enddo
     !
     !     adds to the charge density and converts to real space
     !
     qg(:) = (0.d0, 0.d0)
     do ig = 1, ngm1d
        qg (nl1d (ig) ) = aux (ig) + prho (nl (ig1dto3d (ig) ) )
     enddo
  else
     qg(:) = (0.d0, 0.d0)
  endif
#ifdef __PARA
  call reduce (2 * nrx3, qg)
#endif
  dimz = alat * celldm (3)
  do ig = 1, nr3
     qgr (ig) =  DBLE (qg (ig) )
     qgi (ig) = AIMAG (qg (ig) )
  enddo
  call cft (qgr, qgi, nr3, nr3, nr3, 1)
  do ig = 1, nr3
     plan (ig) = qgr (ig) * omega / dimz
  enddo
  deallocate (aux, qgm)

  return
end subroutine addusdens1d
