!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
subroutine dvscf (nu_i, dvloc, xq_x)
!-----------------------------------------------------------------------
!
!   It reads the variation of the charge density from a file and
!   calculates the variation of the local part of the variation of the
!   K-S potential.
!
#include "machine.h"
  use pwcom
  use phcom
  use d3com
  !
  implicit none
  integer :: nu_i
  ! input: mode under consideration
  real (8) :: xq_x (3)
  ! input: coordinates of the q point
  complex (8) :: dvloc (nrxx)
  ! output: local part of the variation
  !         of the K_S potential
  !
  ! Local variables
  !

  integer :: iudrho_x, ig, ir, mu, na, nt
  ! unit containing the charge variation
  ! counter on G points
  ! counter on real mesh
  ! counter on modes
  ! counter on atoms
  ! the type of atom

  real (kind=DP) :: qg2, gtau
  ! the modulus of (q+G)^2
  ! auxiliary variable: g*tau

  complex (kind=DP) ::  guexp
  ! auxiliary variable: g*u*exp(gtau)

  real (kind=DP), pointer :: vloc_x (:,:)
  ! the local potential at G+q
  complex (kind=DP), pointer :: u_x(:,:), drc_x (:,:)
  complex (kind=DP), allocatable :: aux1 (:), aux2 (:)
  ! the transformation modes patterns
  ! contain drho_core for all atomic types
  logical :: q_eq_zero
  ! true if xq equal zero

  allocate  (aux1( nrxx))    
  allocate  (aux2( nrxx))    

  q_eq_zero = xq_x(1).eq.0.d0 .and. xq_x(2).eq.0.d0 .and. xq_x(3).eq.0.d0
  if (q_eq_zero) then
     u_x   => ug0
     drc_x => d0rc
     vloc_x => vlocg0
     iudrho_x = iud0rho
  else
     u_x   => u
     drc_x => drc
     vloc_x => vlocq
     iudrho_x = iudrho
  endif

  call davcio_drho (aux2, lrdrho, iudrho_x, nu_i, - 1)
  do ir = 1, nrxx
     dvloc (ir) = aux2(ir) * dmuxc(ir,1,1)
  enddo

  call cft3 (aux2, nr1, nr2, nr3, nrx1, nrx2, nrx3, - 1)
  call setv (2 * nrxx, 0.d0, aux1, 1)
  do ig = 1, ngm
     qg2 = (g(1,ig)+xq_x(1))**2 + (g(2,ig)+xq_x(2))**2 + (g(3,ig)+xq_x(3))**2
     if (qg2.gt.1.d-8) then
        aux1(nl(ig)) = e2 * fpi * aux2(nl(ig)) / (tpiba2 * qg2)
     endif

  enddo

  if (nlcc_any) call setv (2 * nrxx, 0.d0, aux2, 1)
  do na = 1, nat
     mu = 3 * (na - 1)
     if (abs(u_x(mu+1,nu_i)) + abs(u_x(mu+2,nu_i)) + &
         abs(u_x(mu+3,nu_i)) .gt. 1.0d-12) then
        nt = ityp (na)
        do ig = 1, ngm

           gtau = tpi * ( (g(1,ig) + xq_x(1)) * tau(1,na) + &
                          (g(2,ig) + xq_x(2)) * tau(2,na) + &
                          (g(3,ig) + xq_x(3)) * tau(3,na) )

           guexp = tpiba * ( (g(1,ig) + xq_x(1)) * u_x(mu+1,nu_i) + &
                             (g(2,ig) + xq_x(2)) * u_x(mu+2,nu_i) + &
                             (g(3,ig) + xq_x(3)) * u_x(mu+3,nu_i) ) * &
                   DCMPLX(0.d0,-1.d0) * DCMPLX(cos(gtau),-sin(gtau))
           aux1 (nl(ig)) = aux1 (nl(ig)) + vloc_x (ig,nt) * guexp
           if (nlcc(nt)) aux2 (nl(ig)) = aux2 (nl(ig)) + drc_x(ig,nt) * guexp
        enddo
     endif

  enddo
  call cft3 (aux1, nr1, nr2, nr3, nrx1, nrx2, nrx3, + 1)

  call DAXPY (2 * nrxx, 1.d0, aux1, 1, dvloc, 1)
  if (nlcc_any) then
     call cft3 (aux2, nr1, nr2, nr3, nrx1, nrx2, nrx3, + 1)
     do ir = 1, nrxx
        aux2 (ir) = aux2 (ir) * dmuxc(ir,1,1)
     enddo
     call DAXPY (2 * nrxx, 1.d0, aux2, 1, dvloc, 1)

  endif

  if (doublegrid) call cinterpolate (dvloc, dvloc, - 1)
  deallocate (aux1)
  deallocate (aux2)
  return
end subroutine dvscf
