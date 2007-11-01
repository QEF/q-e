!
! Copyright (C) 2001-2005 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
subroutine do_elf (elf)
  !-----------------------------------------------------------------------
  !
  !  calculatation of the electron localization function
  !
  !  elf = 1/(1+d**2)
  !
  !          where
  !
  !  d = ( t(r) - t_von_Weizacker(r) ) / t_Thomas-Fermi(r)
  !
  !          and
  !
  !  t (r) = (hbar**2/2m) * \sum_{k,i} |grad psi_{k,i}|**2
  !
  !  t_von_Weizaecker(r) = (hbar**2/2m) * 0.25 * |grad rho(r)|**2/rho
  !  t_von_Weizaecker(r) == t_noninteracting-boson
  !
  !  t_Thomas-Fermi (r) = (hbar**2/2m) * 3/5 * (3*pi**2)**(2/3) * rho**(5/3)
  !
  !
#include "f_defs.h"
  USE kinds, ONLY: DP
  USE constants, ONLY: pi
  USE cell_base, ONLY: omega, tpiba, tpiba2
  USE gvect, ONLY: nr1,nr2,nr3, nrx1,nrx2,nrx3, nrxx, gcutm, ecutwfc, &
       dual, g, ngm, nl, nlm
  USE gsmooth, ONLY : nls, nlsm, nr1s, nr2s, nr3s, &
                      nrx1s, nrx2s, nrx3s, nrxxs, doublegrid
  USE io_files, ONLY: iunwfc, nwordwfc
  USE klist, ONLY: nks, xk
  USE lsda_mod, ONLY: nspin
  USE scf, ONLY: rho
  USE symme, ONLY: nsym, ftau, s
  USE wvfct, ONLY: npw, igk, g2kin, nbnd, wg, gamma_only
  USE wavefunctions_module,  ONLY: evc
  !
  ! I/O variables
  !
  implicit none
  real(DP) :: elf (nrxx)
  !
  ! local variables
  !
  integer :: i, j, k, ng, ibnd, ik, is
  real(DP) :: gv(3), w1, d, fac
  real(DP), allocatable :: kkin (:), tbos (:)
  complex(DP), allocatable :: aux (:), aux2 (:)
  !
  call infomsg ('do_elf', 'elf + US not fully implemented')
  !
  allocate (kkin( nrxx))    
  allocate (aux ( nrxxs))
  aux(:) = (0.d0,0.d0)
  kkin(:) = 0.d0
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
     !
     do ibnd = 1, nbnd
        do j = 1, 3
           aux(:) = (0.d0,0.d0)
           w1 = wg (ibnd, ik) / omega
           do i = 1, npw
              gv (j) = (xk (j, ik) + g (j, igk (i) ) ) * tpiba
              aux (nls(igk (i) ) ) = CMPLX (0d0, gv (j) ) * evc (i, ibnd)
              IF (gamma_only) THEN
                 aux (nlsm(igk (i) ) ) = CMPLX (0d0, -gv (j) ) * &
                      CONJG ( evc (i, ibnd) )
              END IF
           enddo
           call cft3s (aux, nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, 2)
           do i = 1, nrxxs
              kkin(i) = kkin(i) + w1 * (DBLE(aux(i))**2 + AIMAG(aux(i))**2)
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
  call poolreduce (nrxxs, kkin)
#endif
  !
  ! interpolate the local kinetic energy to the dense grid
  ! Note that for US PP this term is incomplete: it contains
  ! only the contribution from the smooth part of the wavefunction
  !
  if (doublegrid) call interpolate (kkin, kkin, 1)
#ifdef __PARA
  call psymrho (kkin, nrx1, nrx2, nrx3, nr1, nr2, nr3, nsym, s, ftau)
#else
  call symrho  (kkin, nrx1, nrx2, nrx3, nr1, nr2, nr3, nsym, s, ftau)
#endif
  !
  ! Calculates the bosonic kinetic density, stored in tbos
  !          aux --> charge density in Fourier space
  !         aux2 --> iG * rho(G)
  !
  deallocate (aux)
  allocate ( tbos(nrxx), aux2(nrxx), aux(nrxx) )
  tbos(:) = 0.d0
  !
  ! put the total (up+down) charge density in rho%of_r(*,1)
  !
  do is = 2, nspin
     rho%of_r (:, 1) =  rho%of_r (:, 1) + rho%of_r (:, is)
  enddo
  !
  aux(:) = CMPLX ( rho%of_r(:, 1), 0.d0 )
  call cft3 (aux, nr1, nr2, nr3, nrx1, nrx2, nrx3, - 1)
  !
  do j = 1, 3
     aux2(:) = (0.d0,0.d0)
     do i = 1, ngm
        aux2(nl(i)) = aux(nl(i)) * CMPLX (0.0d0, g(j,i)*tpiba)
     enddo
     IF (gamma_only) THEN
        do i = 1, ngm
           aux2(nlm(i)) = aux(nlm(i)) * CMPLX (0.0d0,-g(j,i)*tpiba)
        enddo
     END IF

     call cft3 (aux2, nr1, nr2, nr3, nrx1, nrx2, nrx3, 1)
     do i = 1, nrxx
        tbos (i) = tbos (i) + DBLE(aux2(i))**2
     enddo
  enddo
  !
  ! Calculates ELF
  !
  fac = 5.d0 / (3.d0 * (3.d0 * pi**2) ** (2.d0 / 3.d0) )
  elf(:) = 0.d0
  do i = 1, nrxx
     if (rho%of_r (i,1) > 1.d-30) then
        d = fac / rho%of_r(i,1)**(5d0/3d0) * (kkin(i)-0.25d0*tbos(i)/rho%of_r(i,1))
        elf (i) = 1.0d0 / (1.0d0 + d**2)
     endif
  enddo
  deallocate (aux, aux2, tbos, kkin)
  return
end subroutine do_elf
