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
  USE ions_base, ONLY : nat, ityp, tau
  USE kinds,     ONLY : DP
  USE fft_base,  ONLY : dfftp
  USE fft_interfaces, ONLY : fwfft, invfft
  use pwcom
  USE uspp_param, ONLY: upf
  USE uspp,       ONLY: nlcc_any
  use phcom
  use d3com
  !
  implicit none
  integer :: nu_i
  ! input: mode under consideration
  real (DP) :: xq_x (3)
  ! input: coordinates of the q point
  complex (DP) :: dvloc (dfftp%nnr)
  ! output: local part of the variation
  !         of the K_S potential
  !
  ! Local variables
  !

  integer :: iudrho_x, ig, ir, mu, na, nt
  ! unit containing the charge variation
  ! countes

  real (DP) :: qg2, gtau
  ! the modulus of (q+G)^2
  ! auxiliary variable: g*tau

  complex (DP) ::  guexp
  ! auxiliary variable: g*u*exp(gtau)

  real (DP), pointer :: vloc_x (:,:)
  ! the local potential at G+q
  complex (DP), pointer :: u_x(:,:), drc_x (:,:)
  complex (DP), allocatable :: aux1 (:), aux2 (:)
  ! the transformation modes patterns
  ! contain drho_core for all atomic types
  logical :: q_eq_zero
  ! true if xq equal zero

  allocate  (aux1( dfftp%nnr))
  allocate  (aux2( dfftp%nnr))

  q_eq_zero = xq_x(1) == 0.d0 .and. xq_x(2) == 0.d0 .and. xq_x(3) == 0.d0
  if (q_eq_zero) then
     u_x   => ug0
     if (nlcc_any) drc_x => d0rc
     vloc_x => vlocg0
     iudrho_x = iud0rho
  else
     u_x   => u
     if (nlcc_any) drc_x => drc
     vloc_x => vlocq
     iudrho_x = iudrho
  endif

  call davcio_drho (aux2, lrdrho, iudrho_x, nu_i, - 1)

! IT: Warning, if you uncomment the following line,
! you have to precompute the response core charge density
! and pass it as the input to dv_of_drho.
! 
!  call dv_of_drho (aux2(1), .true., ???)

!  dvloc = aux2(:)

!  deallocate (aux1, aux2)

!  return

  dvloc (:) = aux2(:) * dmuxc(:,1,1)
  CALL fwfft ('Dense', aux2, dfftp)

  aux1 (:) = (0.d0, 0.d0)
  do ig = 1, ngm
     qg2 = (g(1,ig)+xq_x(1))**2 + (g(2,ig)+xq_x(2))**2 + (g(3,ig)+xq_x(3))**2
     if (qg2 > 1.d-8) then
        aux1(nl(ig)) = e2 * fpi * aux2(nl(ig)) / (tpiba2 * qg2)
     endif
  enddo

  if (nlcc_any) aux2 (:) = (0.d0, 0.d0)
  do na = 1, nat
     mu = 3 * (na - 1)
     if (abs(u_x(mu+1,nu_i)) + abs(u_x(mu+2,nu_i)) + &
         abs(u_x(mu+3,nu_i)) > 1.0d-12) then
        nt = ityp (na)
        do ig = 1, ngm

           gtau = tpi * ( (g(1,ig) + xq_x(1)) * tau(1,na) + &
                          (g(2,ig) + xq_x(2)) * tau(2,na) + &
                          (g(3,ig) + xq_x(3)) * tau(3,na) )

           guexp = tpiba * ( (g(1,ig) + xq_x(1)) * u_x(mu+1,nu_i) + &
                             (g(2,ig) + xq_x(2)) * u_x(mu+2,nu_i) + &
                             (g(3,ig) + xq_x(3)) * u_x(mu+3,nu_i) ) * &
                   (0.d0,-1.d0) * CMPLX(cos(gtau),-sin(gtau),kind=DP)
           aux1 (nl(ig)) = aux1 (nl(ig)) + vloc_x (ig,nt) * guexp
           if (upf(nt)%nlcc) then
              aux2 (nl(ig)) = aux2 (nl(ig)) + drc_x(ig,nt) * guexp
           end if
        enddo
     endif

  enddo
  CALL invfft ('Dense', aux1, dfftp)
  dvloc (:) = dvloc(:) + aux1 (:)
  if (nlcc_any) then
     CALL invfft ('Dense', aux2, dfftp)
     dvloc (:) = dvloc(:) + aux2 (:) * dmuxc(:,1,1)
  endif

  if (doublegrid) call cinterpolate (dvloc, dvloc, - 1)
  deallocate (aux1)
  deallocate (aux2)
  return
end subroutine dvscf
