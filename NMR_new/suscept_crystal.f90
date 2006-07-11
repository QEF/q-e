!
! Copyright (C) 2001-2005 Quantum-ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "f_defs.h"
!-----------------------------------------------------------------------
SUBROUTINE suscept_crystal
  !-----------------------------------------------------------------------
  !
  ! ... Compute the "bare" susceptibility as in Eq.(64-65) of
  ! ... PRB 63, 245101 (2001)
  ! ... add more comments
  !
  USE kinds,                       ONLY : DP
  USE io_global,                   ONLY : stdout
  USE io_files,                    ONLY : nwordwfc, iunwfc
  USE cell_base,                   ONLY : at, bg, omega, tpiba, tpiba2
  USE wavefunctions_module,        ONLY : evc
  USE klist,                       ONLY : nks, nkstot, wk, xk, nelec
  USE wvfct,                       ONLY : nbnd, npwx, npw, igk, wg, g2kin, &
                                          current_k
  USE lsda_mod,                    ONLY : current_spin, lsda, isk
  USE becmod,                      ONLY : becp  
  USE symme,                       ONLY : nsym, s, ftau
  USE mp_global,                   ONLY : my_pool_id, npool, me_pool, root_pool
  USE pwcom
  USE nmr_module

  !-- local variables ----------------------------------------------------
  IMPLICIT NONE
  complex(dp), allocatable, dimension(:,:,:) :: p_evc, vel_evc, g_vel_evc
  complex(dp), allocatable :: aux(:,:)

  ! Q tensor of eq. (65) (pGv => HH in Paratec, vGv => VV in Paratec)
  real(dp) :: q_pGv(3,3,-1:1), q_vGv(3,3,-1:1)

  ! F tensor of eq. (64)
  real(dp) :: f_pGv(3,3,-1:1), f_vGv(3,3,-1:1)

  ! chi_bare tensor of eq. (64)
  real(dp) :: chi_bare_pGv(3,3), chi_bare_vGv(3,3)

  ! f-sum rule
  real(dp) :: f_sum(3,3)

  real(dp) :: tmp(3,3), q(3), braket
  integer :: ia, ib, ik, ipol, jpol, i, ibnd, isign
  complex(dp), external :: ZDOTC
  !-----------------------------------------------------------------------

  ! allocate memory
  allocate ( p_evc(npwx,nbnd,3), vel_evc(npwx,nbnd,3), &
             aux(npwx,nbnd), g_vel_evc(npwx,nbnd,3) )

  ! zero the Q tensors
  q_pGv(:,:,:) = 0.d0
  q_vGv(:,:,:) = 0.d0

  ! zero the current and the field
  j_bare(:,:,:,:) = (0.d0,0.d0)
  b_ind(:,:,:) = (0.d0,0.d0)

  write(stdout, '(5X,''Computing the magnetic susceptibility'')')
  !====================================================================
  ! loop over k-points
  !====================================================================
  do ik = 1, nks
#ifdef __PARA
    if (me_pool == root_pool) &
    write(*, '(5X,''k-point #'',I5,'' of '',I5,6X,''pool #'',I3)') &
      ik, nks, my_pool_id+1
#else
    write(stdout, '(5X,''k-point #'',I5,'' of '',I5)') ik, nks
#endif
    current_k = ik
    current_spin = isk(ik)

    ! initialize at k-point k 
    call gk_sort(xk(1,ik), ngm, g, ecutwfc/tpiba2, npw, igk, g2kin)
    g2kin(:) = g2kin(:) * tpiba2
    call init_us_2(npw,igk,xk(1,ik),vkb)

    ! read wfcs from file and compute becp
    call davcio (evc, nwordwfc, iunwfc, ik, -1)
    call ccalbec (nkb, npwx, npw, nbnd, becp, vkb, evc)

    ! this is the case q = 0 (like the case of the f-sum rule)
    q(:) = 0.d0
    !!!write(*,'(''q='',3(F12.4))') q
    call compute_u_kq(ik, q)
    evc = evq

    ! compute p_k|evc>, v_k|evc> and G_k v_{k,k}|evc>
    call apply_operators

    !------------------------------------------------------------------
    ! f-sum rule (pGv term only) 
    !------------------------------------------------------------------
    do ia = 1, 3 
      do ib = 1, 3
        do ibnd = 1, nbnd_occ(ik)
          braket = 2.d0*real(ZDOTC(npw, p_evc(1,ibnd,ia), 1, &
                                        g_vel_evc(1,ibnd,ib), 1), DP)
          f_sum(ia,ib) = f_sum(ia,ib) + wg(ibnd,ik) * braket
        enddo
      enddo
    enddo

    !------------------------------------------------------------------
    ! pGv and vGv contribution to chi_{bare}
    !------------------------------------------------------------------
    do i = 1, 3
      call add_to_tensor(q_pGv(:,:,0), p_evc, g_vel_evc)
      call add_to_tensor(q_vGv(:,:,0), vel_evc, g_vel_evc)
    enddo

    !------------------------------------------------------------------
    ! loop over -q and +q
    !------------------------------------------------------------------
    do isign = -1, 1, 2

      ! loop over cartesian directions
      do i = 1, 3
        ! set the q vector
        q(:) = 0.d0
        q(i) = dble(isign) * q_nmr
        !!!write(*,'(''q='',3(F12.4))') q
                       
        ! compute the wfcs at k+q
        call compute_u_kq(ik, q)

        ! compute p_k|evc>, v_k|evc> and G_{k+q} v_{k+q,k}|evc>
        call apply_operators

        ! pGv and vGv contribution to chi_bare
        call add_to_tensor(q_pGv(:,:,isign), p_evc, g_vel_evc)
        call add_to_tensor(q_vGv(:,:,isign), vel_evc, g_vel_evc)

        ! now the j_bare term  
        call add_to_current(j_bare(:,:,:,current_spin), evc, g_vel_evc)
      enddo  ! i

    enddo  ! isign
  enddo  ! ik

#ifdef __PARA
  ! reduce over G-vectors
  call reduce(9, f_sum)
  call reduce(3*9, q_pGv)
  call reduce(3*9, q_vGv)
#endif

#ifdef __PARA
  ! reduce over k-points
  call poolreduce(9, f_sum)
  call poolreduce(3*9, q_pGv)
  call poolreduce(3*9, q_vGv)
  call poolreduce(nrxxs*nspin*9, j_bare)
#endif

  !====================================================================
  ! print out results
  !====================================================================
  write(stdout,'(5X,''End of magnetic susceptibility calculation'')')
  write(stdout,*)

  ! f-sum rule
  if (iverbosity > 0) then
    write(stdout, '(5X,''f-sum rule:'')')
    write(stdout, tens_fmt) f_sum
  endif
  call sym_cart_tensor(f_sum)
  write(stdout, '(5X,''f-sum rule (symmetrized):'')')
  write(stdout, tens_fmt) f_sum

  ! F_{ij} = (2 - \delta_{ij}) Q_{ij}
  do ipol = 1, 3
    do jpol = 1, 3
      f_pGv(ipol,jpol,:) = 2.d0*q_pGv(ipol,jpol,:)
      if (ipol == jpol) f_pGv(ipol,jpol,:) = q_pGv(ipol,jpol,:)

      f_vGv(ipol,jpol,:) = 2.d0*q_vGv(ipol,jpol,:)
      if (ipol == jpol) f_vGv(ipol,jpol,:) = q_vGv(ipol,jpol,:)
    enddo
  enddo

  ! compute chi_bare both pGv and vGv terms
  chi_bare_pGv(:,:) = f_pGv(:,:,1) - 2.d0*f_pGv(:,:,0) + f_pGv(:,:,-1)
  chi_bare_pGv(:,:) = -0.5d0 * chi_bare_pGv(:,:) / (c * q_nmr * tpiba)**2
  if (iverbosity > 0) then
    write(stdout, '(5X,''chi_bare pGv (HH) in paratec units:'')')
    write(stdout, '(3(5X,3(F12.6,2X)/))') chi_bare_pGv(:,:) * c**2
  endif
  call sym_cart_tensor(chi_bare_pGv)
  if (iverbosity > 0) then
    write(stdout, '(3(5X,3(F12.6,2X)/))') chi_bare_pGv(:,:) * c**2
  endif

  chi_bare_vGv(:,:) = f_vGv(:,:,1) - 2.d0*f_vGv(:,:,0) + f_vGv(:,:,-1)
  chi_bare_vGv(:,:) = -0.5d0 * chi_bare_vGv(:,:) / (c * q_nmr * tpiba)**2
  if (iverbosity > 0) then
    write(stdout, '(5X,''chi_bare vGv (VV) in paratec units:'')')
    write(stdout, '(3(5X,3(F12.6,2X)/))') chi_bare_vGv(:,:) * c**2
  endif
  call sym_cart_tensor(chi_bare_vGv)
  if (iverbosity > 0) then
    write(stdout, '(3(5X,3(F12.6,2X)/))') chi_bare_vGv(:,:) * c**2
  endif

  ! convert from atomic units to 10^{-6} cm^3 / mol
  tmp(:,:) = chi_bare_pGv(:,:) * 1d6 * a0_to_cm**3.d0 * avogadro
  write(stdout, '(5X,''chi_bare pGv (HH) in 10^{-6} cm^3/mol:'')')
  write(stdout, tens_fmt) tmp(:,:)

  tmp(:,:) = chi_bare_vGv(:,:) * 1d6 * a0_to_cm**3.d0 * avogadro
  write(stdout, '(5X,''chi_bare vGv (VV) in 10^{-6} cm^3/mol:'')')
  write(stdout, tens_fmt) tmp(:,:)

  !--------------------------------------------------------------------
  ! now get the current, induced field and chemical shifts
  !--------------------------------------------------------------------
  chi_bare_pGv(:,:) = chi_bare_pGv(:,:) / omega
  j_bare(:,:,:,:) = j_bare(:,:,:,:) / (2.d0 * q_nmr * tpiba * c * omega)

  !nsym = 1
  ! either you symmetrize the current ...
  do i = 1, nspin
#ifdef __PARA
    call psymmetrize_field(j_bare(:,:,:,i),1)
#else
    call symmetrize_field(j_bare(:,:,:,i),1)
#endif
  enddo

  ! compute induced field
  do ipol = 1, 3
    call biot_savart(ipol)
  enddo

  do i = 1, nspin
    if (trim(filcurr) /= '') &
      call write_tensor_field(filcurr, i, j_bare(1,1,1,i))
  enddo
  if (trim(filfield) /= '') &
    call write_tensor_field(filfield, 0, b_ind_r)


  ! ... or you symmetrize the induced field
  !call symmetrize_field(b_ind_r,0)
  !call field_to_reciprocal_space
  
  ! compute chemical shifts
  call compute_sigma_bare(chi_bare_pGv)

  deallocate(p_evc, vel_evc, aux, g_vel_evc, j_bare, b_ind)

  CONTAINS


  !====================================================================
  ! compute p_k|evc>, v_k|evc> and G_k v_{k+q,k}|evc>
  !====================================================================
  SUBROUTINE apply_operators
    implicit none
    integer ipol

    do ipol = 1, 3
      call apply_p(evc, p_evc(1,1,ipol), ik, ipol, q)
      call apply_vel(evc, vel_evc(1,1,ipol), ik, ipol, q)
      ! necessary because aux is overwritten by subroutine greenfunction
      aux(:,:) = vel_evc(:,:,ipol)
      call greenfunction(ik, aux, g_vel_evc(1,1,ipol), q)
    enddo
  END SUBROUTINE apply_operators



  !====================================================================
  ! add contribution the Q tensors
  ! Q_{\alpha,\beta} += <(e_i \times ul)_\alpha | (e_i \times ur)_\beta>
  !====================================================================
  SUBROUTINE add_to_tensor(qt, ul, ur)
    implicit none
    real(dp), intent(inout) :: qt(3,3)
    complex(dp), intent(in) :: ul(npwx,nbnd,3), ur(npwx,nbnd,3)
    real(dp) :: braket
    integer :: ibnd, ia, ib, comp_ia, comp_ib, ind(3,3), mult(3,3)

    ! index for the cross product
    ind(:,1) = (/ 1, 3, 2/);  mult(:,1) = (/ 0,-1, 1 /)
    ind(:,2) = (/ 3, 2, 1/);  mult(:,2) = (/ 1, 0,-1 /)
    ind(:,3) = (/ 2, 1, 3/);  mult(:,3) = (/-1, 1, 0 /)

    do ia = 1, 3    ! ia = alpha
      comp_ia = ind(ia,i)
      if (mult(ia,i) == 0) cycle

      do ib = 1, 3    ! ib = beta
        comp_ib = ind(ib,i)
        if (mult(ib,i) == 0) cycle

        do ibnd = 1, nbnd_occ(ik)
          braket = real(ZDOTC(npw, ul(1,ibnd,comp_ia), 1, &
                                   ur(1,ibnd,comp_ib), 1), DP)
          qt(ia,ib) = qt(ia,ib) + wg(ibnd,ik) * &
                      braket * mult(ia,i) * mult(ib,i)
        enddo  ! ibnd

      enddo  ! ib
    enddo  ! ia
!#ifdef __PARA
!    call reduce(9, qt)
!#endif
  END SUBROUTINE add_to_tensor



  !====================================================================
  ! add contribution the the current
  ! j(r)_{\alpha,\beta} += <ul|J(r)|(B\times e_i \cdot ur)>
  !====================================================================
  SUBROUTINE add_to_current(j, ul, ur)
    implicit none
    real(dp), intent(inout) :: j(nrxxs,3,3)
    complex(dp), intent(in) :: ul(npwx,nbnd), ur(npwx,nbnd,3)
    real(dp) :: braket, fact
    integer :: ibdir, icomp, ind(3,3), mult(3,3)

    ! index for the cross product
    ind(:,1) = (/ 1, 3, 2/);  mult(:,1) = (/ 0,-1, 1 /)
    ind(:,2) = (/ 3, 2, 1/);  mult(:,2) = (/ 1, 0,-1 /)
    ind(:,3) = (/ 2, 1, 3/);  mult(:,3) = (/-1, 1, 0 /)

    ! loop over B direction
    do ibdir = 1, 3
      if (i == ibdir) cycle
      icomp = ind(ibdir, i)
      fact = real(mult(ibdir,i)*isign)
      call j_para(fact, ul(1,1), ur(1,1,icomp), ik, q, j(1,1,ibdir))
    enddo
  END SUBROUTINE add_to_current
  
END SUBROUTINE suscept_crystal
