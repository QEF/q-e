! Copyright (C) 2001-2005 Quantum-ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "f_defs.h"
!-----------------------------------------------------------------------
SUBROUTINE apply_vel(psi, vel_psi, ik, ipol, q)
  !-----------------------------------------------------------------------
  !
  ! ... Apply the velocity operator
  ! ...   v = p + dV^{NL}_{k+q,k}/dk
  ! ...
  ! ... Here we use Hartree atomic units, so that:
  ! ...   V^{NL} => V^{NL} * ryd_to_hartree
  !-----------------------------------------------------------------------
  USE kinds,                ONLY : DP
  USE ions_base,            ONLY : nat, ityp, ntyp => nsp
  USE klist,                ONLY : xk
  USE wvfct,                ONLY : nbnd, npwx, npw, igk  
  USE becmod,               ONLY : becp, calbec, allocate_bec, deallocate_bec
  USE pwcom,                ONLY : nkb, vkb, tpiba
  USE gipaw_module,         ONLY : q_gipaw, nbnd_occ

  !-- paramters ----------------------------------------------------------
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: ipol       ! cartesian direction (1..3)
  INTEGER, INTENT(IN) :: ik         ! k-point
  COMPLEX(DP), INTENT(IN) :: psi(npwx,nbnd)
  COMPLEX(DP), INTENT(OUT) :: vel_psi(npwx,nbnd)
  REAL(DP), INTENT(IN) :: q(3)

  !-- local variables ----------------------------------------------------
  real(dp), parameter :: ryd_to_hartree = 0.5d0
  complex(dp), allocatable :: aux(:,:), vkb_save(:,:)
  real(dp) :: dk, dxk(3)
  integer :: isign
  logical :: q_is_zero


  ! first apply p
  call apply_p(psi, vel_psi, ik, ipol, q)

  call start_clock('apply_vel')

  ! set dk
  dk = q_gipaw/2.d0
  
  ! if no projectors, return
  if (nkb == 0) return

  ! check if |q| is zero
  q_is_zero = .false.
  if (sqrt(q(1)*q(1)+q(2)*q(2)+q(3)*q(3)) < 1d-8) q_is_zero = .true.

  ! allocate temporary arrays, save old NL-potential
  allocate(aux(npwx,nbnd), vkb_save(npwx,nkb))
  vkb_save = vkb

#if 1
  !====================================================================
  ! compute (1/2|dk|) ( V^{NL}_{k+dk+q,k+dk} |psi> - 
  !                     V^{NL}_{k-dk+q,k-dk} |psi> )
  !====================================================================
  call allocate_bec(nkb,nbnd)
  do isign = -1,1,2
    dxk(:) = xk(:,ik)
    dxk(ipol) = dxk(ipol) + isign * dk     ! k + dk

    ! compute <\beta(k \pm dk)| and project on |psi>
    call init_us_2_no_phase(npw, igk, dxk, vkb)
    !it was: call ccalbec (nkb, npwx, npw, nbnd_occ(ik), becp, vkb, psi)
    call calbec (npw, vkb, psi, becp, nbnd_occ(ik))

    ! |q|!=0 => compute |\beta(k \pm dk + q)>
    if (.not. q_is_zero) then
      dxk(:) = dxk(:) + q(:)
      call init_us_2_no_phase(npw, igk, dxk, vkb)
    endif

    ! apply |\beta(k \pm dk+q)>D<\beta(k \pm dk)| to |psi>
    aux = (0.d0,0.d0)
    call add_vuspsi(npwx, npw, nbnd_occ(ik), psi, aux)
    vel_psi = vel_psi + dble(isign) * ryd_to_hartree * aux/(2.d0*dk*tpiba)
  enddo
  call deallocate_bec
#else

  do isign = -1,1,2
    ! compute <\beta(k)| and project on |psi>
    call init_us_2_no_phase(npw, igk, xk(1,ik), vkb)
    !it was: call ccalbec (nkb, npwx, npw, nbnd_occ(ik), becp, vkb, psi)
    call calbec (npw, vkb, psi, becp, nbnd_occ(ik))

    dxk(ipol) = xk(ipol,ik) + isign * dk     ! k + dk
    call init_us_2_no_phase(npw, igk, dxk, vkb)

    ! apply |\beta(k \pm dk+q)>D<\beta(k)| to |psi>
    aux = (0.d0,0.d0)
    call add_vuspsi(npwx, npw, nbnd_occ(ik), psi, aux)
    vel_psi = vel_psi + 0.5d0*dble(isign) * ryd_to_hartree * aux/(2.d0*dk*tpiba)


    ! compute <\beta(k \pm dk)| and project on |psi>
    dxk(:) = xk(:,ik)
    dxk(ipol) = dxk(ipol) + isign * dk     ! k + dk
    call init_us_2_no_phase(npw, igk, xk(1,ik), vkb)
    !it was: call ccalbec (nkb, npwx, npw, nbnd_occ(ik), becp, vkb, psi)
    call calbec (npw, vkb, psi, becp, nbnd_occ(ik))

    dxk(:) = xk(:,ik)
    call init_us_2_no_phase(npw, igk, dxk, vkb)

    ! apply |\beta(k+q)>D<\beta(k \pm dk)| to |psi>
    aux = (0.d0,0.d0)
    call add_vuspsi(npwx, npw, nbnd_occ(ik), psi, aux)
    vel_psi = vel_psi + 0.5d0*dble(isign) * ryd_to_hartree * aux/(2.d0*dk*tpiba)
  enddo
#endif

  ! restore NL-potential at k
  vkb = vkb_save
  
  ! free memory
  deallocate(aux, vkb_save)

  call stop_clock('apply_vel')

END SUBROUTINE apply_vel
