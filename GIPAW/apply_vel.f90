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
  USE becmod,               ONLY : becp
  USE pwcom
  USE gipaw_module

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
  do isign = -1,1,2
    dxk(:) = xk(:,ik)
    dxk(ipol) = dxk(ipol) + isign * dk     ! k + dk

    ! compute <\beta(k \pm dk)| and project on |psi>
    call init_us_2_no_phase(npw, igk, dxk, vkb)
    call ccalbec (nkb, npwx, npw, nbnd_occ(ik), becp, vkb, psi)

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
#else

  do isign = -1,1,2
    ! compute <\beta(k)| and project on |psi>
    call init_us_2_no_phase(npw, igk, xk(1,ik), vkb)
    call ccalbec (nkb, npwx, npw, nbnd_occ(ik), becp, vkb, psi)

    dxk(ipol) = xk(ipol,ik) + isign * dk     ! k + dk
    call init_us_2_no_phase(npw, igk, dxk, vkb)

    ! apply |\beta(k \pm dk+q)>D<\beta(k)| to |psi>
    aux = (0.d0,0.d0)
    call add_vuspsi(npwx, npw, nbnd_occ(ik), psi, aux)
    vel_psi = vel_psi + 0.5*dble(isign) * ryd_to_hartree * aux/(2.d0*dk*tpiba)


    ! compute <\beta(k \pm dk)| and project on |psi>
    dxk(:) = xk(:,ik)
    dxk(ipol) = dxk(ipol) + isign * dk     ! k + dk
    call init_us_2_no_phase(npw, igk, xk(1,ik), vkb)
    call ccalbec (nkb, npwx, npw, nbnd_occ(ik), becp, vkb, psi)

    dxk(:) = xk(:,ik)
    call init_us_2_no_phase(npw, igk, dxk, vkb)

    ! apply |\beta(k+q)>D<\beta(k \pm dk)| to |psi>
    aux = (0.d0,0.d0)
    call add_vuspsi(npwx, npw, nbnd_occ(ik), psi, aux)
    vel_psi = vel_psi + 0.5*dble(isign) * ryd_to_hartree * aux/(2.d0*dk*tpiba)
  enddo
#endif

  ! restore NL-potential at k
  vkb = vkb_save
  
  ! free memory
  deallocate(aux, vkb_save)

  call stop_clock('apply_vel')

END SUBROUTINE apply_vel




! Here are some previous and testing version of the subroutine
! it looks like that the analytic derivative has some troubles.
!
! The numerical one (that computes just |dbeta/dk>) is just fine
! but is is more involved that the routine above.

#if 0
  real(DP), allocatable  :: gk (:,:)
  complex(DP), allocatable :: ps2(:,:,:), dvkb (:,:), dvkb1 (:,:), &
       work (:,:), becp2(:,:), becp1(:,:), aux(:,:), vkb_save(:,:)
  complex(DP), external :: ZDOTC
  integer :: ig, na, ibnd, ikb, jkb, nt, ih, jh, ijkb0
  integer :: isign
  real(DP) :: axis(3)



  do ibnd = 1, nbnd_occ(ik) 
  print*, ibnd
  print '(2(F12.6,2X))', vel_psi(1:5,ibnd)
  print*
  enddo

  !!!vl_psi = 0.d0

  ! this is the old routine, applying |beta><dbeta/dk| + |dbeta/dk><beta|
  ! compute (k+G)/|k+G|
  allocate (gk(3,npwx))
  do ig = 1, npw
     gk (1, ig) = (xk (1, ik) + g (1, igk (ig) ) ) * tpiba
     gk (2, ig) = (xk (2, ik) + g (2, igk (ig) ) ) * tpiba
     gk (3, ig) = (xk (3, ik) + g (3, igk (ig) ) ) * tpiba
     g2kin (ig) = gk (1, ig) **2 + gk (2, ig) **2 + gk (3, ig) **2

     if (g2kin (ig) < 1.0d-10) then
        gk (1, ig) = 0.d0
        gk (2, ig) = 0.d0
        gk (3, ig) = 0.d0
     else
        gk (1, ig) = gk (1, ig) / sqrt (g2kin (ig) )
        gk (2, ig) = gk (2, ig) / sqrt (g2kin (ig) )
        gk (3, ig) = gk (3, ig) / sqrt (g2kin (ig) )
     endif
  enddo

  ! allocate variable for the derivative of the beta's
  allocate(dvkb(npwx,nkb), dvkb1(npwx,nkb))
  dvkb = (0.d0,0.d0)
  dvkb1 = (0.d0,0.d0)
  work = (0.d0,0.d0)

  ! and this is the contribution from nonlocal pseudopotentials
  axis = 0.d0; axis(ipol) = 1.d0
  !!!call gen_us_dy (ik, at (1, ipol), dvkb1)
  call gen_us_dy (ik, axis(1), dvkb1)
  call gen_us_dj (ik, dvkb)
  print*, dvkb1(2,:)

  jkb = 0
  do nt = 1, ntyp
     do na = 1, nat
        if (nt == ityp (na)) then
           do ikb = 1, nh (nt)
              jkb = jkb + 1
              do ig = 1, npw
                 !!!work (ig,jkb) = dvkb1 (ig, jkb) + dvkb (ig, jkb) * &
                 !!!     (at (1, ipol) * gk (1, ig) + &
                 !!!     at (2, ipol) * gk (2, ig) + &
                 !!!     at (3, ipol) * gk (3, ig) )
                 work (ig,jkb) = dvkb1 (ig, jkb) + dvkb (ig, jkb) * gk(ipol,ig)
              enddo
           enddo
        endif
     enddo
  enddo
  work(1,3) = (0.d0,-0.133684d0)
  jkb = 3
  write(98,'(2F12.6,2X))') work(:,jkb)
  print*

  ! this was the previous numerical derivative
  !allocate(dvkb(npwx,nkb), dvkb1(npwx,nkb))
  dvkb = (0.d0,0.d0)
  dvkb1 = (0.d0,0.d0)
  work = (0.d0,0.d0)
  vel_psi = 0.d0

  ! compute it by finite differences
  work = (0.d0,0.d0)
  dxk(:) = xk(:,ik); dxk(ipol) = dxk(ipol) + dk
  call init_us_2(npw,igk,dxk,work)

  dxk(:) = xk(:,ik); dxk(ipol) = dxk(ipol) - dk
  call init_us_2(npw,igk,dxk,dvkb)

  work(:,:) = work(:,:) - dvkb(:,:)
  work = 0.5d0 * work / (dk * tpiba)
  write(99,'(2F12.6,2X))') work(:,jkb)
  print*
 
 
  allocate( becp2(nkb,nbnd), becp1(nkb,nbnd))
  becp2(:,:) = (0.d0,0.d0) 
  becp1(:,:) = (0.d0,0.d0) 
  call ccalbec (nkb, npwx, npw, nbnd, becp2, work, psi)
  call init_us_2(npw,igk,xk(1,ik),vkb)
  call ccalbec (nkb, npwx, npw, nbnd, becp1, vkb, psi)

  ijkb0 = 0
  allocate (ps2 ( nkb, nbnd, 2))
  ps2(:,:,:)=(0.d0,0.d0)
  do nt = 1, ntyp
     do na = 1, nat
        if (nt == ityp (na)) then
           do ih = 1, nh (nt)
              ikb = ijkb0 + ih
              do jh = 1, nh (nt)
                 jkb = ijkb0 + jh
                 do ibnd = 1, nbnd_occ(ik)
                    !! CHECK factor (-1.d0,0.d0)
                    !!!ps2(ikb,ibnd,1) = ps2(ikb,ibnd,1)+ becp2(jkb,ibnd)* &
                    !!!     (-1.d0,0.d0) * (deeq(ih,jh,na,current_spin) &
                    !!!     -et(ibnd,ik)*qq(ih,jh,nt))
                    !!!ps2(ikb,ibnd,2) = ps2(ikb,ibnd,2) +becp1(jkb,ibnd,ik) * &
                    !!!     (-1.d0,0.d0) * (deeq(ih,jh,na,current_spin)&
                    !!!     -et(ibnd,ik)*qq(ih,jh,nt))
                    ps2(ikb,ibnd,1) = ps2(ikb,ibnd,1) + becp2(jkb,ibnd)* &
                         (1.d0,0.d0) * deeq(ih,jh,na,current_spin)
                    ps2(ikb,ibnd,2) = ps2(ikb,ibnd,2) + becp1(jkb,ibnd) * &
                         (1.d0,0.d0) * deeq(ih,jh,na,current_spin)
                 enddo
              enddo
           enddo
           ijkb0=ijkb0+nh(nt)
        end if
     end do
  end do
  if (ikb /= nkb .OR. jkb /= nkb) call errore ('apply_vel', 'unexpected error',1)
  !print*, 'ps2(:,:,1)='
  !print*, ps2(:,:,1)
  !print*, 'ps2(:,:,2)='
  !print*, ps2(:,:,2)

  CALL ZGEMM( 'N', 'N', npw, nbnd_occ(ik), nkb, &
       (1.d0,0.d0), vkb(1,1), npwx, ps2(1,1,1), nkb, (1.d0,0.d0), &
       vel_psi(1,1), npwx )

  CALL ZGEMM( 'N', 'N', npw, nbnd_occ(ik), nkb, &
       (1.d0,0.d0), work(1,1), npwx, ps2(1,1,2), nkb, (1.d0,0.d0), &
       vel_psi(1,1), npwx )

  do ibnd = 1, nbnd_occ(ik) 
  print*, ibnd
  print '(2(F12.6,2X))', vel_psi(1:5,ibnd)
  print*
  enddo
  STOP

  deallocate (ps2)
  deallocate (gk)
#endif



