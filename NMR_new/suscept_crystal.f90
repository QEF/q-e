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
  
  !<apsi>
  USE paw,                         ONLY: paw_vkb, paw_becp, paw_nkb, aephi, &
                                         psphi, paw_nh, paw_nhtol, &
                                         paw_nhtom, paw_indv, paw_nbeta
  USE ions_base, ONLY : nat
  !</apsi>

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

  integer :: ia, ib, ik, ipol, jpol, i, ibnd, isign
  real(dp) :: tmp(3,3), q(3), braket, sigma_bare(3,3,nat)
  real(dp) :: diamagnetic_corr_tensor(3,3,nat)
  real(dp) :: paramagnetic_corr_tensor(3,3,nat)
  real(dp) :: sigma_diamagnetic(3,3,nat)
  real(dp) :: sigma_paramagnetic(3,3,nat)
  complex(dp), external :: ZDOTC
  !-----------------------------------------------------------------------

  ! allocate memory
  allocate ( p_evc(npwx,nbnd,3), vel_evc(npwx,nbnd,3), &
             aux(npwx,nbnd), g_vel_evc(npwx,nbnd,3) )

  ! zero the Q tensors
  q_pGv(:,:,:) = 0.0_dp
  q_vGv(:,:,:) = 0.0_dp

  ! zero the current and the field
  j_bare(:,:,:,:) = (0.0_dp,0.0_dp)
  b_ind(:,:,:) = (0.0_dp,0.0_dp)
  
  sigma_diamagnetic = 0.0_dp
  sigma_paramagnetic = 0.0_dp
  
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
    
    !<apsi>
    call init_paw_2_no_phase (npw, igk, xk (1, ik), paw_vkb)
    call ccalbec (paw_nkb, npwx, npw, nbnd, paw_becp, paw_vkb, evc)
    diamagnetic_corr_tensor = 0.0
    call diamagnetic_correction ( diamagnetic_corr_tensor )
    sigma_diamagnetic = sigma_diamagnetic + diamagnetic_corr_tensor
    !</apsi>
    
    ! this is the case q = 0 (like the case of the f-sum rule)
    q(:) = 0.0_dp
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
          braket = 2.0_dp*real(ZDOTC(npw, p_evc(1,ibnd,ia), 1, &
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
        q(:) = 0.0_dp
        q(i) = dble(isign) * q_nmr
        !!!write(*,'(''q='',3(F12.4))') q
        
        ! compute the wfcs at k+q
        call compute_u_kq(ik, q)
        
        ! compute p_k|evc>, v_k|evc> and G_{k+q} v_{k+q,k}|evc>
        call apply_operators
        
        !<apsi>
        call init_paw_2_no_phase (npw, igk, xk(:,ik)+q(:), paw_vkb)
        call paramagnetic_correction ( paramagnetic_corr_tensor )
        call add_to_sigma_para ( paramagnetic_corr_tensor, sigma_paramagnetic )
        !</apsi>
        
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
  call poolreduce(9*nat, sigma_diamagnetic)
  call poolreduce(9*nat, sigma_paramagnetic)
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
      f_pGv(ipol,jpol,:) = 2.0_dp*q_pGv(ipol,jpol,:)
      if (ipol == jpol) f_pGv(ipol,jpol,:) = q_pGv(ipol,jpol,:)

      f_vGv(ipol,jpol,:) = 2.0_dp*q_vGv(ipol,jpol,:)
      if (ipol == jpol) f_vGv(ipol,jpol,:) = q_vGv(ipol,jpol,:)
    enddo
  enddo
  
  ! compute chi_bare both pGv and vGv terms
  chi_bare_pGv(:,:) = f_pGv(:,:,1) - 2.0_dp*f_pGv(:,:,0) + f_pGv(:,:,-1)
  chi_bare_pGv(:,:) = -0.5_dp * chi_bare_pGv(:,:) / (c * q_nmr * tpiba)**2
  if (iverbosity > 0) then
    write(stdout, '(5X,''chi_bare pGv (HH) in paratec units:'')')
    write(stdout, '(3(5X,3(F12.6,2X)/))') chi_bare_pGv(:,:) * c**2
  endif
  call sym_cart_tensor(chi_bare_pGv)
  if (iverbosity > 0) then
    write(stdout, '(3(5X,3(F12.6,2X)/))') chi_bare_pGv(:,:) * c**2
  endif
  
  chi_bare_vGv(:,:) = f_vGv(:,:,1) - 2.0_dp*f_vGv(:,:,0) + f_vGv(:,:,-1)
  chi_bare_vGv(:,:) = -0.5_dp * chi_bare_vGv(:,:) / (c * q_nmr * tpiba)**2
  if (iverbosity > 0) then
    write(stdout, '(5X,''chi_bare vGv (VV) in paratec units:'')')
    write(stdout, '(3(5X,3(F12.6,2X)/))') chi_bare_vGv(:,:) * c**2
  endif
  call sym_cart_tensor(chi_bare_vGv)
  if (iverbosity > 0) then
    write(stdout, '(3(5X,3(F12.6,2X)/))') chi_bare_vGv(:,:) * c**2
  endif

  ! convert from atomic units to 10^{-6} cm^3 / mol
  tmp(:,:) = chi_bare_pGv(:,:) * 1e6_dp * a0_to_cm**3.0_dp * avogadro
  write(stdout, '(5X,''chi_bare pGv (HH) in 10^{-6} cm^3/mol:'')')
  write(stdout, tens_fmt) tmp(:,:)

  tmp(:,:) = chi_bare_vGv(:,:) * 1e6_dp * a0_to_cm**3.0_dp * avogadro
  write(stdout, '(5X,''chi_bare vGv (VV) in 10^{-6} cm^3/mol:'')')
  write(stdout, tens_fmt) tmp(:,:)

  !--------------------------------------------------------------------
  ! now get the current, induced field and chemical shifts
  !--------------------------------------------------------------------
  chi_bare_pGv(:,:) = chi_bare_pGv(:,:) / omega
  j_bare(:,:,:,:) = j_bare(:,:,:,:) / (2.0_dp * q_nmr * tpiba * c * omega)
  
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
  call compute_sigma_bare(chi_bare_pGv, sigma_bare)
  
  call compute_sigma_diamagnetic( sigma_diamagnetic )
  
  call compute_sigma_paramagnetic( sigma_paramagnetic )
  
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
  

  !====================================================================
  ! ...
  !====================================================================
  SUBROUTINE diamagnetic_correction ( diamagnetic_tensor )
    
    USE atom,       ONLY : r, rab
    USE ions_base,  ONLY : nat, ityp, ntyp => nsp
    USE nmr_module, ONLY : c
    
    implicit none
    
    ! Arguments
    real(dp), intent(inout):: diamagnetic_tensor(3,3,nat)

    integer :: l1, m1, lm1, l2, m2, lm2, ih, ikb, nbs1, jh, jkb, nbs2
    integer :: nt, ibnd, na, lm, nrc, ijkb0
    complex(dp) , allocatable :: efg_corr(:,:)
    complex(dp) :: bec_product
    
    allocate (efg_corr(lmaxx**2,nat))
    efg_corr = 0.0_dp
    
    !  rc=1.60_dp
    !  nrc=count(r(1:msh(1),1).le.rc)
    !
    ! calculate radial integration on atom site 
    ! <aephi|1/r^3|aephi>-<psphi|1/r^3|psphi>
    !
    
    !
    !  calculation of the reconstruction part
    !
    
    do ibnd = 1, nbnd
       ijkb0 = 0
       do nt = 1, ntyp
          do na = 1, nat
             if (ityp (na) .eq.nt) then
                do ih = 1, paw_nh (nt)
                   ikb = ijkb0 + ih
                   nbs1=paw_indv(ih,nt)
                   l1=paw_nhtol(ih,nt)
                   m1=paw_nhtom(ih,nt)
                   lm1=m1+l1**2
                   do jh = 1, paw_nh (nt) 
                      jkb = ijkb0 + jh
                      nbs2=paw_indv(jh,nt)
                      l2=paw_nhtol(jh,nt)
                      m2=paw_nhtom(jh,nt)
                      lm2=m2+l2**2 
                      
                      bec_product = paw_becp(jkb,ibnd) &
                           * CONJG(paw_becp(ikb,ibnd))
                      
                      !<apsi> s/non-trace-zero component
                      ! 2/3 to separate the non-trace vanishing component
                      ! 1/(2c^2) from the equation (59) in PM-PRB
                      IF ( l1 == l2 .AND. m1 == m2 ) THEN
                         diamagnetic_tensor(1,1,na) &
                              = diamagnetic_tensor(1,1,na) &
                              + 2.0_dp / 3.0_dp * bec_product &
                              * radial_integral_diamagnetic(nbs1,nbs2,nt) &
                              * wg(ibnd,ik) / (2.0_dp*c**2)
                         diamagnetic_tensor(2,2,na) &
                              = diamagnetic_tensor(2,2,na) &
                              + 2.0_dp / 3.0_dp * bec_product &
                              * radial_integral_diamagnetic(nbs1,nbs2,nt) &
                              * wg(ibnd,ik) / (2.0_dp*c**2)
                         diamagnetic_tensor(3,3,na) &
                              = diamagnetic_tensor(3,3,na) &
                              + 2.0_dp / 3.0_dp * bec_product &
                              * radial_integral_diamagnetic(nbs1,nbs2,nt) &
                              * wg(ibnd,ik) / (2.0_dp*c**2)
                      END IF
                      
                      ! 2/3 to separate the non-trace vanishing component
                      do lm = 5, 9
                         efg_corr(lm,na) =  efg_corr(lm,na) &
                              + bec_product / 3.0_dp &
                              * radial_integral_diamagnetic(nbs1,nbs2,nt) &
                              * ap(lm,lm1,lm2) * wg(ibnd,ik) / (2.0_dp*c**2)
                      enddo
                   enddo
                enddo
                ijkb0 = ijkb0 + paw_nh (nt)
             endif
          enddo
       enddo
    enddo
    
!    write(6,'("DDD1",5F14.8)') efg_corr(5:9,1)
!    write(6,'("DDD2",5F14.8)') efg_corr(5:9,2)
    
    !
    !  transform in cartesian coordinates
    !
    
    efg_corr(5:9,:nat) = - sqrt(4.0_dp * pi/5.0_dp) * efg_corr(5:9,:nat)
    
    diamagnetic_tensor(1,1,:) = diamagnetic_tensor(1,1,:) &
         + sqrt(3.0_dp) * efg_corr(8,:) - efg_corr(5,:)
    diamagnetic_tensor(2,2,:) = diamagnetic_tensor(2,2,:) &
         - sqrt(3.0_dp) * efg_corr(8,:) - efg_corr(5,:)
    diamagnetic_tensor(3,3,:) = diamagnetic_tensor(3,3,:) &
         + efg_corr(5,:) * 2.0_dp
    diamagnetic_tensor(1,2,:) = diamagnetic_tensor(1,2,:) &
         +  efg_corr(9,:) * sqrt(3.0_dp)
    diamagnetic_tensor(2,1,:) = diamagnetic_tensor(1,2,:)
    diamagnetic_tensor(1,3,:) = diamagnetic_tensor(1,3,:) &
         - efg_corr(6,:) * sqrt(3.0_dp)
    diamagnetic_tensor(3,1,:) = diamagnetic_tensor(1,3,:)
    diamagnetic_tensor(2,3,:) = diamagnetic_tensor(2,3,:) &
         - efg_corr(7,:) * sqrt(3.0_dp)
    diamagnetic_tensor(3,2,:) = diamagnetic_tensor(2,3,:)
    
    ! efg_corr(5,:) = 3z^2-1
    ! efg_corr(6,:) = -xz
    ! efg_corr(7,:) = -yz
    ! efg_corr(8,:) = x^2-y^2
    ! efg_corr(9,:) = xy
    
    deallocate ( efg_corr )
    
  END SUBROUTINE diamagnetic_correction
  

  !====================================================================
  ! ...
  !====================================================================
  SUBROUTINE paramagnetic_correction ( paramagnetic_tensor )
    
    USE ions_base,  ONLY : nat, ityp, ntyp => nsp
    USE nmr_module, ONLY : c
    
    implicit none
    
    ! Arguments
    real(dp), intent(inout):: paramagnetic_tensor(3,3,nat)
    
    integer :: l1, m1, lm1, l2, m2, lm2, ih, ikb, nbs1, jh, jkb, nbs2
    integer :: nt, ibnd, na, lm, j, ijkb0, ipol
    complex(dp) :: bec_product
    complex(dp) , allocatable :: efg_corr(:,:)
    
    integer, parameter :: ng_ = 27, lmax2_ = 16
    integer :: mg, i1, i2, i3
    real(DP) :: g_ (3, ng_), gg_ (ng_)
    real(DP) :: ylm_ (ng_,lmax2_)
    
    !--------------------------------------------------------------------------
    
    allocate (efg_corr(3,nat))
    
    !
    !  calculation of the reconstruction part
    !
    
    do ipol = 1, 3 
       
       if ( ipol == i ) cycle !TESTTESTTEST
       
       call ccalbec (paw_nkb, npwx, npw, nbnd, paw_becp2, paw_vkb, &
            g_vel_evc(1,1,ipol))
       
       efg_corr = 0.0_dp
       
       do ibnd = 1, nbnd
          ijkb0 = 0
          do nt = 1, ntyp
             do na = 1, nat
                
                if (ityp (na) .eq.nt) then
                   do ih = 1, paw_nh (nt)
                      ikb = ijkb0 + ih
                      nbs1=paw_indv(ih,nt)
                      l1=paw_nhtol(ih,nt)
                      m1=paw_nhtom(ih,nt)
                      lm1=m1+l1**2
                      
                      do jh = 1, paw_nh (nt) 
                         jkb = ijkb0 + jh
                         nbs2=paw_indv(jh,nt)
                         l2=paw_nhtol(jh,nt)
                         m2=paw_nhtom(jh,nt)
                         lm2=m2+l2**2
                         
                         if ( l1 /= l2 ) cycle
                         
                         bec_product = CONJG(paw_becp(ikb,ibnd)) &
                              * paw_becp2(jkb,ibnd)
                         
                         efg_corr(1,na) = efg_corr(1,na) &
                              + bec_product &
                              * radial_integral_paramagnetic(nbs1,nbs2,nt) &
                              * lx ( lm1, lm2 ) * wg(ibnd,ik) / c ** 2
                         efg_corr(2,na) = efg_corr(2,na) &
                              + bec_product &
                              * radial_integral_paramagnetic(nbs1,nbs2,nt) &
                              * ly ( lm1, lm2 ) * wg(ibnd,ik) / c ** 2
                         efg_corr(3,na) = efg_corr(3,na) &
                              + bec_product &
                              * radial_integral_paramagnetic(nbs1,nbs2,nt) &
                              * lz ( lm1, lm2 ) * wg(ibnd,ik) / c ** 2
                         
                      enddo
                   enddo
                   ijkb0 = ijkb0 + paw_nh (nt)
                endif
             enddo
          enddo
       enddo
       
       paramagnetic_tensor ( :, ipol, : ) = REAL ( efg_corr, dp )
       
       write(6,'("DDD1",2I3,3(F16.7,2X)') &
            ipol, i*isign, REAL ( efg_corr(1:3,1) ) * 1e6
       !write(6,'("DDD2",2I3,3(F16.7,2X)') &
       !     ipol, i*isign, REAL ( efg_corr(1:3,2) ) * 1e6
       
    end do
    
    deallocate(efg_corr)
    
  END SUBROUTINE paramagnetic_correction
  
  !====================================================================
  ! ...
  !====================================================================
  SUBROUTINE add_to_sigma_para( paramagnetic_correction, sigma_paramagnetic )
    implicit none
    real(dp), intent(in) :: paramagnetic_correction(3,3,nat)
    real(dp), intent(inout) :: sigma_paramagnetic(3,3,nat)
    real(dp) :: fact
    integer :: ibdir, icomp, ipol, ind(3,3), mult(3,3)
    
    ! index for the cross product
    ind(:,1) = (/ 1, 3, 2/);  mult(:,1) = (/ 0,-1, 1 /)
    ind(:,2) = (/ 3, 2, 1/);  mult(:,2) = (/ 1, 0,-1 /)
    ind(:,3) = (/ 2, 1, 3/);  mult(:,3) = (/-1, 1, 0 /)
    
    ! loop over B direction
    do ibdir = 1, 3
      if (i == ibdir) cycle
      icomp = ind(ibdir, i)
      fact = real(mult(ibdir,i)*isign)
      
      do ipol = 1, 3
         sigma_paramagnetic ( ipol, icomp, : ) &
              = sigma_paramagnetic ( ipol, icomp, : ) &
              + fact * paramagnetic_correction ( ipol, ibdir, : ) &
              / ( 2 * q_nmr * tpiba )
         
!   Do iq=1, 3    ! loop over all q-points
!      ...
!      Do iperm0 =1,-1,-2
!         b0 = Mod(iq+iperm0+2,3) +1 ! gives different b0 for iperm0
!         p0 = Mod(iq-iperm0+2,3) +1 ! gives p0 = abs(q x b0)
!         ...
!         Call take_nonloc_deriv_kq_k(p0,k_gspace, rk(1),u_k(1,i),&
!         ...
!         para_corr_shift_tot(b0,:,:,:) = para_corr_shift_tot(b0,:,:,:)-&
!              Real(iperm0,dp)*Real(iqsign,dp)*para_corr_shift*&
!              kpoints%w(irk)/qmag/dtwo/Real(crys%nspin,dp)
         
      end do
    enddo
  END SUBROUTINE add_to_sigma_para
  
END SUBROUTINE suscept_crystal
