!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------


subroutine do_elf (elf)
  !-----------------------------------------------------------------------
  !  calculate electron localization function
  !
  !  elf = 1/(1+d**2)
  !
  !          where
  !
  !  d = ( t(r) - t_von_Weitzacker(r) ) / t_Thomas-Fermi(r)
  !
  !          and
  !
  !  t (r) = (hbar**2/2m) * \sum_{k,i} |grad psi_{k,i}|**2
  !
  !  t_von_Weitzaecker(r) = (hbar**2/2m) * 0.25 * |grad rho(r)|**2/rho
  !  t_von_Weitzaecker(r) == t_noninteracting-boson
  !
  !  t_Thomas-Fermi (r) =  (hbar**2/2m) * 3/5 * (3*pi**2)**(2/3) * rho**(5
  !
  !
#include "machine.h"
  use parameters
  use pwcom
  !
  ! I/O variables
  !
  implicit none
  real(kind=DP) :: elf (nrxx)
  !
  ! local variables
  !
  integer :: i, j, k, ng, ibnd, ik, is
  real(kind=DP) :: gv(3), w1, d, arho, fac
  real(kind=DP), allocatable :: kkin (:), tbos (:)
  complex(kind=DP), allocatable :: aux (:), aux2 (:)

  allocate (kkin( nrxx))    
  allocate (tbos( nrxx))    
  allocate (aux ( nrxx))    
  allocate (aux2( nrxx))    
  call errore ('do_elf', 'elf + US not fully implemented', - 1)
  !
  ! put the total (up+down) charge density in rho(*,1)
  !
  do is = 2, nspin
     call DAXPY (nrxx, 1.d0, rho (1, is), 1, rho (1, 1), 1)
  enddo
  call setv (2 * nrxx, 0d0, aux, 1)

  call setv (nrxx, 0d0, kkin, 1)
  !
  ! Calculates local kinetic energy, stored in kkin
  !
  do ik = 1, nks
     !
     !    prepare the indices of this k point
     !

     call gk_sort (xk (1, ik), ngm, g, ecutwfc / tpiba2, npw, igk, g2kin)
     !
     !   reads the eigenfunctions
     !

     call davcio (evc, nwordwfc, iunwfc, ik, - 1)
     do ibnd = 1, nbnd

        do j = 1, 3
           call setv (2 * nrxx, 0d0, aux, 1)

           w1 = wg (ibnd, ik) / omega
           do i = 1, npw
              gv (j) = (xk (j, ik) + g (j, igk (i) ) ) * tpiba
              aux (nl (igk (i) ) ) = cmplx (0d0, gv (j) ) * evc (i, ibnd)

           enddo

           call cft3 (aux, nr1, nr2, nr3, nrx1, nrx2, nrx3, 1)
           do i = 1, nrxx
              kkin (i) = kkin (i) + w1 * (real (aux (i) ) **2 + DIMAG (aux (i) ) &
                   **2)

           enddo
           ! j
        enddo
        ! ibnd

     enddo
     ! ik

  enddo
#ifdef __PARA
  !
  ! reduce local kinetic energy across pools
  !
  call poolreduce (nrxx, kkin)
  call psymrho (kkin, nrx1, nrx2, nrx3, nr1, nr2, nr3, nsym, s, &
       ftau)
#else
  call symrho (kkin, nrx1, nrx2, nrx3, nr1, nr2, nr3, nsym, s, ftau)
#endif
  !
  ! Calculates the bosonic kinetic density, stored in tbos
  !          aux --> charge density in Fourier space
  !         aux2 --> iG * rho(G)
  !
  call setv (nrxx, 0d0, tbos, 1)
  call setv (2 * nrxx, 0d0, aux, 1)

  call DCOPY (nrxx, rho, 1, aux, 2)

  call cft3 (aux, nr1, nr2, nr3, nrx1, nrx2, nrx3, - 1)
  do j = 1, 3
     call setv (2 * nrxx, 0d0, aux2, 1)
     do i = 1, ngm
        aux2 (nl (i) ) = aux (nl (i) ) * cmplx (0.0d0, g (j, i) * tpiba)
     enddo
     call cft3 (aux2, nr1, nr2, nr3, nrx1, nrx2, nrx3, 1)
     do i = 1, nrxx
        tbos (i) = tbos (i) + real (aux2 (i) ) **2
     enddo

  enddo
  !
  ! Calculates ELF
  !

  fac = 5.d0 / (3.d0 * (3.d0 * pi**2) ** (2.d0 / 3.d0) )
  call setv (nrxx, 0d0, elf, 1)
  do i = 1, nrxx
     arho = abs (rho (i, 1) )
     if (arho.gt.1.d-30) then
        d = fac / (rho (i, 1) ** (5d0 / 3d0) ) * (kkin (i) - 0.25d0 * &
             tbos (i) / rho (i, 1) )
        elf (i) = 1.0d0 / (1.0d0 + d**2)
     endif

  enddo
  deallocate (aux)
  deallocate (aux2)
  deallocate (tbos)
  deallocate (kkin)
  return
end subroutine do_elf
