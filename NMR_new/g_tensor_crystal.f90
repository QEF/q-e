!
! Copyright (C) 2001-2005 Quantum-ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "f_defs.h"
!-----------------------------------------------------------------------
SUBROUTINE g_tensor_crystal
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
  USE lsda_mod,                    ONLY : current_spin, lsda, isk, nspin
  USE becmod,                      ONLY : becp  
  USE symme,                       ONLY : nsym, s, ftau
  USE scf,                         ONLY : vr, vltot, rho
  USE gvect,                       ONLY : ngm, nr1, nr2, nr3, nrx1, nrx2, &
                                          nrx3, nrxx, nlm, g
  USE mp_global,                   ONLY : my_pool_id
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

  real(dp) :: q(3), braket, delta_rmc
  integer :: ia, ib, ik, ipol, jpol, i, ibnd, isign, ispin
  complex(dp), external :: ZDOTC

  integer :: s_maj, s_min
  real(dp) :: e_hartree, charge, s_weight, rho_diff, d_omega
  real(dp), allocatable :: grad_vr(:,:), v_local(:,:)
  real(dp), allocatable :: grad_vh(:,:), vh(:)
  real(dp), dimension ( 3, 3 ) :: delta_g_rmc, delta_g_bare, delta_g_soo, &
       delta_g_soo_2, delta_g_paramagn, delta_g_diamagn, delta_g_total
  real(dp) :: g_e = 2.0023192778_DP, g_prime, units_Ry2Ha = 0.5_DP
  
  logical :: tevaluate_chi, tcalculate_correct_delta_g_soo
  
  !-----------------------------------------------------------------------
  
  !<apsi> TMPTMPTMP Until the reconstruction has been implemented
  delta_g_paramagn = 0.0_DP
  delta_g_diamagn = 0.0_DP
  !</apsi>
  
  g_prime = 2 * ( g_e - 1 )
  
  !<apsi> TMPTMPTMP Move to input or make it default?
  tcalculate_correct_delta_g_soo = .true.
  !</apsi>
  if ( tcalculate_correct_delta_g_soo ) then
     tevaluate_chi = .true.
  else
     tevaluate_chi = .false.
  end if
  
  ! Select majority and minority spin components
  rho_diff = SUM ( rho ( :, 1 ) - rho ( :, nspin ) )
  if ( rho_diff > +1.0e-3 ) then
     s_maj = 1
     s_min = nspin
  else if ( rho_diff < -1.0e-3 ) then
     s_maj = nspin
     s_min = 1
  else
     write ( stdout, * ) "WARNING: rho_diff zero!"
  end if
  
  ! allocate memory
  allocate ( p_evc(npwx,nbnd,3), vel_evc(npwx,nbnd,3), &
       aux(npwx,nbnd), g_vel_evc(npwx,nbnd,3) )
  
  ! zero the Q tensors
  q_pGv(:,:,:) = 0.d0
  q_vGv(:,:,:) = 0.d0
  
  ! zero the current and the field
  j_bare(:,:,:,:) = (0.d0,0.d0)
  b_ind(:,:,:) = (0.d0,0.d0)
  
  delta_rmc = 0.0_DP
  
  write(stdout, '(5X,''Computing the magnetic susceptibility'')')
  
  !====================================================================
  ! loop over k-points
  !====================================================================
  do ik = 1, nks
#ifdef __PARA
    if (my_pool_id == 0) &
    write(*, '(5X,''k-point #'',I5,'' of '',I5,'' (my_pool_id='',I3)') &
      ik, nks, my_pool_id
#else
    write(*, '(5X,''k-point #'',I5,'' of '',I5)') ik, nks
#endif
    current_k = ik
    current_spin = isk(ik)
    
    if ( current_spin == s_maj ) then
       s_weight = +1
    else
       s_weight = -1
    end if
    
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
    
    do ibnd = 1, nbnd_occ(ik)
       delta_rmc = delta_rmc - s_weight / c ** 2 * g_e &
            * wg(ibnd,ik) * SUM ( g2kin(:) * conjg(evc(:,ibnd)) * evc(:,ibnd) )
    end do
    
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
#ifdef __PARA
    call reduce(9, f_sum)
#endif

    !------------------------------------------------------------------
    ! pGv and vGv contribution to chi_{bare}
    !------------------------------------------------------------------
    if ( tevaluate_chi ) then
      do i = 1, 3
        call add_to_tensor(q_pGv(:,:,0), p_evc, g_vel_evc)
        call add_to_tensor(q_vGv(:,:,0), vel_evc, g_vel_evc)
      enddo
    end if
    
    !------------------------------------------------------------------
    ! loop over -q and +q
    !------------------------------------------------------------------
    do isign = -1, 1, 2

      ! loop over cartesian directions
      do i = 1, 3
        if (iverbosity > 10) then
          write(stdout,*) "  QQQ: ", i * isign
        end if
        ! set the q vector
        q(:) = 0.d0
        q(i) = dble(isign) * q_nmr
        !!!write(*,'(''q='',3(F12.4))') q
                       
        ! compute the wfcs at k+q
        call compute_u_kq(ik, q)

        ! compute p_k|evc>, v_k|evc> and G_{k+q} v_{k+q,k}|evc>
        call apply_operators

        if ( tevaluate_chi ) then
          ! pGv and vGv contribution to chi_bare
          call add_to_tensor(q_pGv(:,:,isign), p_evc, g_vel_evc)
          call add_to_tensor(q_vGv(:,:,isign), vel_evc, g_vel_evc)
        end if

        ! now the j_bare term  
        call add_to_current(j_bare(:,:,:,current_spin), evc, g_vel_evc)
      enddo  ! i

    enddo  ! isign
  enddo  ! ik
  ! TODO: put a lot of poolreduce here
#ifdef __PARA
  call poolreduce(9, f_sum)
  call poolreduce(9, q_pGv)
  call poolreduce(9, q_vGv)
  ! TODO: non working yet in parallel!!!
#endif
  
  !====================================================================
  ! print out results
  !====================================================================
  write(stdout,'(5X,''End of magnetic susceptibility calculation'',/)')
  
  ! f-sum rule
  if (iverbosity > 0) then
    write(stdout, '(5X,''f-sum rule:'')')
    write(stdout, '(3(5X,3(F12.6,2X)/))') f_sum
  endif
  call sym_cart_tensor(f_sum)
  write(stdout, '(5X,''f-sum rule (symmetrized):'')')
  write(stdout, '(3(5X,3(F12.6,2X)/))') f_sum
  
  !--------------------------------------------------------------------
  ! now get the current, induced field and delta_g
  !--------------------------------------------------------------------
  chi_bare_pGv(:,:) = chi_bare_pGv(:,:) / omega
  j_bare(:,:,:,:) = j_bare(:,:,:,:) / (2.d0 * q_nmr * tpiba * c * omega)
  
  ! either you symmetrize the current ...
  do ispin = 1, nspin
    call symmetrize_field(j_bare(:,:,:,ispin),1)
  end do
  
  !
  ! calculate the susceptibility
  !
  ! F_{ij} = (2 - \delta_{ij}) Q_{ij}
  if ( tevaluate_chi ) then
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
  end if
  
  ! compute induced field
  do ipol = 1, 3
    call biot_savart(ipol)
  enddo
  
  ! compute chemical shifts
  !call compute_sigma_bare(chi_bare_pGv)
  
  
  d_omega = omega / REAL ( nr1 * nr2 * nr3, DP )
  
  !
  ! ***************** spin-orbit-bare *******************
  !
  
  ! <apsi> TMPTMPTMP PLEASE CHECK FOR VANDERBILT/HARD GRIDS
  allocate ( grad_vr ( 3, nrxx ), v_local ( nrxx, nspin ) )
  do ispin = 1, nspin
     v_local(:,ispin) = vltot(:) + vr(:,ispin)
  end do
  call gradient ( nrx1, nrx2, nrx3, nr1, nr2, nr3, nrxx, v_local, &
       ngm, g, nl, alat, grad_vr )
  grad_vr = grad_vr * units_Ry2Ha
  deallocate ( v_local )
  ! </apsi>
  
  do ipol = 1, 3
     delta_g_bare ( ipol, 1 ) = SUM ( &
          ( j_bare(:,2,ipol,s_maj)-j_bare(:,2,ipol,s_min)) * grad_vr ( 3, : ) &
          -(j_bare(:,3,ipol,s_maj)-j_bare(:,3,ipol,s_min)) * grad_vr ( 2, : ) )
     delta_g_bare ( ipol, 2 ) = SUM ( &
          ( j_bare(:,3,ipol,s_maj)-j_bare(:,3,ipol,s_min)) * grad_vr ( 1, : ) &
          -(j_bare(:,1,ipol,s_maj)-j_bare(:,1,ipol,s_min)) * grad_vr ( 3, : ) )
     delta_g_bare ( ipol, 3 ) = SUM ( &
          ( j_bare(:,1,ipol,s_maj)-j_bare(:,1,ipol,s_min)) * grad_vr ( 2, : ) &
          -(j_bare(:,2,ipol,s_maj)-j_bare(:,2,ipol,s_min)) * grad_vr ( 1, : ) )
  end do
  deallocate ( grad_vr )
  
#ifdef __PARA
  call reduce(9, delta_g_bare)
#endif
  
  delta_g_bare = delta_g_bare * d_omega
  delta_g_bare = - delta_g_bare * g_prime / 2.d0 / c * 1.0e6
  
  !
  ! ***************** spin-other-orbit *******************
  !
  
  ! This is the form used in the implementation in 'paratec'
  
  ! calculate the spin-other-orbit term a'la paratec:
  !   int_r j_up(r) x v_h[n_unpaired] d^3r
  allocate ( grad_vh ( 3, nrxx ), vh ( nrxx ) )
  call v_h( rho(:,s_maj)-rho(:,s_min), nr1, nr2, nr3, nrx1, nrx2, nrx3, nrxx, &
       nl, ngm, gg, gstart, 1, alat, omega, e_hartree, charge, vh )
  call gradient ( nrx1, nrx2, nrx3, nr1, nr2, nr3, nrxx, vh, &
       ngm, g, nl, alat, grad_vh )
  grad_vh = grad_vh * units_Ry2Ha
  deallocate ( vh )
  
  ! -2*j_bare_dw(r) because of the self-interaction correction:
  !   j_bare(r) - [j_bare_up(r)-j_bare_dw(r)] = -2*j_bare_dw(r)
  do ipol = 1, 3
     delta_g_soo ( ipol, 1 ) = 2 * SUM &
          ( j_bare(:,2,ipol,s_min) * grad_vh ( 3, : ) &
          - j_bare(:,3,ipol,s_min) * grad_vh ( 2, : ) )
     delta_g_soo ( ipol, 2 ) = 2 * SUM &
          ( j_bare(:,3,ipol,s_min) * grad_vh ( 1, : ) &
          - j_bare(:,1,ipol,s_min) * grad_vh ( 3, : ) )
     delta_g_soo ( ipol, 3 ) = 2 * SUM &
          ( j_bare(:,1,ipol,s_min) * grad_vh ( 2, : ) &
          - j_bare(:,2,ipol,s_min) * grad_vh ( 1, : ) )
  end do
  deallocate ( grad_vh )
  
#ifdef __PARA
  call reduce ( 9, delta_g_soo )
#endif
  
  delta_g_soo = delta_g_soo * d_omega
  delta_g_soo = delta_g_soo * 2 / c * 1.0e6
  
  ! This is obtained using the equation (7) of Pickard et Mauri, PRL 88/086043
  if ( tcalculate_correct_delta_g_soo ) then
     !<apsi> The G=0 term is still not yet in
     !</apsi>
     do jpol = 1, 3
        do ipol = 1, 3
           delta_g_soo_2 ( ipol, jpol ) = SUM ( b_ind_r(:,ipol,jpol) &
                * ( rho(:,s_maj)-rho(:,s_min) ) )
        end do
     end do
  end if
  
  delta_g_soo_2 = delta_g_soo_2 * d_omega
  delta_g_soo_2 = delta_g_soo_2 * 2 * 1.0e6
  
  !
  ! ***************** relativistic-mass-correction *******************
  !
  
  delta_rmc = delta_rmc * units_Ry2Ha * 1.0e6
  delta_g_rmc = 0.0_DP
  do i = 1, 3
     delta_g_rmc(i,i) = delta_rmc
  end do
  
  !
  ! ***************** total delta_g *******************
  !
  
  delta_g_total = delta_g_rmc + delta_g_bare + delta_g_diamagn &
       + delta_g_paramagn + delta_g_soo 
  
  write (stdout,*)
  write (stdout,*) '**********************************************'
  write (stdout,*)
  write (stdout,*) 'Delta g - relativistic-mass-correction'
  write (stdout,*) 
  write ( stdout, '(3(5X,3(F12.6,2X)/))' ) delta_g_rmc(:,:)
  write (stdout,*) '**********************************************'
  write (stdout,*)
  write (stdout,*) 'Delta g - spin-orbit-bare'
  write (stdout,*) 
  write ( stdout, '(3(5X,3(F12.6,2X)/))' ) delta_g_bare(:,:)
  write (stdout,*) '**********************************************'
  write (stdout,*)
  write (stdout,*) 'Delta g - spin-orbit diamagnetic correction (GIPAW)'
  write (stdout,*) 
  write ( stdout, '(3(5X,3(F12.6,2X)/))' ) delta_g_diamagn(:,:)
  write (stdout,*) '**********************************************'
  write (stdout,*)
  write (stdout,*) 'Delta g - spin-orbit paramagnetic correction (GIPAW)'
  write (stdout,*) 
  write ( stdout, '(3(5X,3(F12.6,2X)/))' ) delta_g_paramagn(:,:)
  write (stdout,*) '**********************************************'
  write (stdout,*)
  write (stdout,*) 'Delta g - spin-other-orbit'
  write (stdout,*) 
  write ( stdout, '(3(5X,3(F12.6,2X)/))' ) delta_g_soo(:,:)
  write (stdout,*) '**********************************************'
  write (stdout,*)
  write (stdout,*) 'Delta g - spin-other-orbit, version 2'
  write (stdout,*) 
  write ( stdout, '(3(5X,3(F12.6,2X)/))' ) delta_g_soo_2(:,:)
  write (stdout,*) '**********************************************'
  write (stdout,*)
  write (stdout,*) 'Delta g - total'
  write (stdout,*) 
  write ( stdout, '(3(5X,3(F12.6,2X)/))' ) delta_g_total(:,:)
  write (stdout,*) '**********************************************'
  write (stdout,*)
  
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
#ifdef __PARA
    call reduce(9, qt)
#endif
  END SUBROUTINE add_to_tensor
  
  
  !====================================================================
  ! add contribution the the current
  ! j(r)_{\alpha,\beta} += <ul|J(r)|(B\times e_i \cdot ur)>
  !====================================================================
  SUBROUTINE add_to_current(j, ul, ur)
    implicit none
    real(dp), intent(inout) :: j(nrxxs,3,3)
    complex(dp), intent(in) :: ul(npwx,nbnd), ur(npwx,nbnd,3)
    real(dp) :: fact
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
  
END SUBROUTINE g_tensor_crystal
