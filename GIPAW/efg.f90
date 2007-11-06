!
! Copyright (C) 2001-2005 Quantum-ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "f_defs.h"
!-----------------------------------------------------------------------

!******************************************************************************

subroutine efg
  
  USE io_files,     ONLY : nd_nmbr
  USE io_global,    ONLY : stdout
  USE kinds,        ONLY : dp 
  USE parameters,   ONLY : ntypx
  USE constants,    ONLY : pi, tpi, fpi, angstrom_au, rytoev, electronvolt_si
  USE scf,          ONLY : rho
  USE gvect,        ONLY : nr1,nr2,nr3,nrx1,nrx2,nrx3,nrxx,&
                         g,gg,nl,gstart,ngm
                        ! gvectors and parameters for the FFT
  USE cell_base,    ONLY : at,bg         !parameters of the cell
  USE ions_base,    ONLY : nat, atm, tau, ityp, zv
  USE symme,        ONLY : nsym, s, irt
  USE gipaw_module, ONLY : q_efg, job
  
  implicit none
  
  integer :: alpha, beta, ig, na, i
  real(DP) :: eta, Cq
  real(DP) :: fac, trace, arg, e2
  complex(DP), allocatable:: aux(:)
  complex(DP), allocatable:: efgg_el(:,:,:),efgr_el(:,:,:)
  complex(DP), allocatable:: efg_io(:,:,:)
  real(DP), allocatable:: zion(:), efg_corr_tens(:,:,:), efg_total(:,:,:)
  real(DP):: efg_eig(3), v(3)
  complex(DP) :: work(3,3), efg_vect(3,3)
  
  allocate(aux(nrxx))
  allocate(efgg_el(nrxx,3,3))
  allocate(efgr_el(nat,3,3))
  allocate(efg_io(nat,3,3))
  allocate(zion(nat))
  allocate(efg_corr_tens(3,3,nat))
  allocate(efg_total(3,3,nat))
  
  ! e2 = 2.0_dp ! rydberg
  e2 = 1.0_dp  ! hartree
  fac= fpi * e2
  aux(:) = rho%of_r(:,1)
  efgg_el(:,:,:) = (0.0_dp,0.0_dp)
  
  call cft3(aux,nr1,nr2,nr3,nrx1,nrx2,nrx3,-1)
  
  !
  ! calculation of the electric field gradient in the G-space
  !
  do ig= gstart, ngm
     trace = 1.0_dp/3.0_dp * gg(ig)
     
     do alpha=1,3
        efgg_el(ig,alpha,alpha)= -trace
        
        do beta=1,3
           efgg_el(ig,alpha,beta) = ( efgg_el(ig,alpha,beta) &
                + g(alpha,ig) * g(beta,ig)) * fac * aux(nl(ig)) / gg(ig)
        end do
     end do
  end do
  
  ! 
  ! fourier transform on the atomic position
  !
  
  efgr_el=(0.0_dp,0.0_dp)
  do alpha=1,3
     do beta=1,3
        do na=1,nat
           do ig= gstart, ngm
              arg=(tau(1,na)*g(1,ig)+tau(2,na)*g(2,ig)+tau(3,na)*g(3,ig))*tpi
              efgr_el(na,alpha,beta)=efgr_el(na,alpha,beta) &
                   + efgg_el(ig,alpha,beta) * CMPLX(cos(arg),sin(arg))
           end do
        end do
     end do
  end do
  
#ifdef __PARA
  call reduce (2*3*3*nat, efgr_el) !2*, efgr_el is a complex array
#endif
  
  write ( stdout, '( / )' )
  
  do na=1,nat
     do beta=1,3
        write ( stdout, 1000 ) atm(ityp(na)), na, "efgr_el", &
             ( REAL(efgr_el(na,alpha,beta),dp), alpha = 1, 3 )
     end do
     write ( stdout, '( / )' )
  end do
  
1000 FORMAT(1x,a,i3,2x,a,3(1x,f9.6))
  
  !
  ! Ionic contribution
  !
  
  call ewald_dipole (efg_io, zv)
  
  do na=1,nat
     do beta=1,3
        write (stdout,1000) atm(ityp(na)), na, "efg_ion", &
             ( REAL(efg_io(na,alpha,beta),dp), alpha = 1, 3 )
     end do
     write ( stdout, '( / )' )
  end do
  
  call efg_correction ( efg_corr_tens )
  
  !
  ! symmetrise efg_tensor
  !
  
  do na = 1,nat
     call trntns (efg_corr_tens(:,:,na),at, bg, -1)
  end do
  call symz(efg_corr_tens, nsym, s, nat, irt)
  do na = 1,nat
     call trntns (efg_corr_tens(:,:,na),at, bg, 1)
  end do
  
  !
  ! print results
  !
  
  do na=1,nat
     do beta=1,3
        write (stdout,1000) atm(ityp(na)),na,"efg_corr", &
             (2*REAL(efg_corr_tens(alpha,beta,na),dp), alpha = 1, 3 )
     end do
     write ( stdout, '( / )' )
  end do
  
  do na=1,nat
     efg_total(:,:,na) = REAL ( 2 * efg_corr_tens(:,:,na) &
          + efgr_el(na,:,:) + efg_io(na,:,:), dp )
     
     do beta=1,3
        write (stdout,1000) atm(ityp(na)),na,"efg",&
             (efg_total(alpha,beta,na),alpha=1,3)
     end do
     write ( stdout, '( / )' )
     
     do alpha=1,3
        do beta=1,3
           work(beta,alpha) = CMPLX ( efg_total(alpha,beta,na), 0.0_dp )
        end do
     end do
     
     !
     ! diagonalise the tensor to extract the quadrupolar parameters Cq and eta
     !
     
     call cdiagh(3,work,3,efg_eig,efg_vect)
     
     v(2) = efg_eig(2)
     if ( abs(efg_eig(1)) > abs(efg_eig(3)) ) then
        v(1)=efg_eig(1)
        v(3)=efg_eig(3)
     else
        v(1)=efg_eig(3)
        v(3)=efg_eig(1)
     end if
     
     if (abs(v(1))<1e-5) then
        eta = 0.0_dp
     else
        eta = (v(2)-v(3))/v(1)
     end if
     
     Cq = v(1) * q_efg(ityp(na)) * rytoev*2.0_dp * angstrom_au ** 2 &
          * electronvolt_si * 1.d18 / 6.62620_dp
     
     write ( stdout, 1200 ) atm(ityp(na)), na, q_efg(ityp(na)), Cq, eta
     write ( stdout, '( / )' )
  end do
  
1200 FORMAT ( 1x, a, 1x, i3, 5x, 'Q= ', f5.2, ' 10e-30 m^2', &
          5x, ' Cq=', f9.4, ' MHz', 5x, 'eta= ', f8.5 )
  
end subroutine efg

!******************************************************************************

subroutine hyperfine
  
  USE io_files,     ONLY : nd_nmbr
  USE io_global,    ONLY : stdout
  USE kinds,        ONLY : dp 
  USE parameters,   ONLY : ntypx
  USE constants,    ONLY : pi, tpi, fpi, angstrom_au, rytoev, electronvolt_si
  USE scf,          ONLY : rho
  USE gvect,        ONLY : nr1,nr2,nr3,nrx1,nrx2,nrx3,nrxx,&
                         g,gg,nl,gstart,ngm
                        ! gvectors and parameters for the FFT
  USE cell_base,    ONLY : at,bg         ! parameters of the cell
  USE ions_base,    ONLY : nat, atm, tau, ityp, zv
  USE symme,        ONLY : nsym, s, irt
  USE lsda_mod,     ONLY : current_spin, nspin
  USE wvfct,        ONLY : current_k
  USE gipaw_module, ONLY : hfi_nuclear_g_factor, hfi_output_unit, &
                         hfi_isotope, job
  
  implicit none
  
  integer :: alpha, beta, ig, na, i
  INTEGER :: s_weight, spin, s_min, s_maj
  REAL ( dp ) :: common_factor, output_factor, rho_diff
  real(DP) :: trace, arg, e2, fact
  real(DP) :: efg_eig(3), v(3), delta_t(ngm,ntypx)
  complex(DP) :: work(3,3), efg_vect(3,3)
  
  complex(DP), allocatable :: aux(:)
  complex(DP), allocatable :: efgg_el(:,:,:),efgr_el(:,:,:)
  complex(DP), allocatable :: efg_io(:,:,:)
  complex(DP), allocatable :: efgr_fc(:)
  real(DP), allocatable :: efg_corr_tens(:,:,:), efg_total(:,:,:)
  real(DP), allocatable :: fc_recon(:), fc_recon_zora(:)
  complex(DP), allocatable :: efgg_zora(:,:)
  complex(DP), allocatable :: efgr_zora(:)
  
  REAL ( dp ) :: mu0_by_fpi = 1e-7
  REAL ( dp ) :: mu_n = 5.051e-27
  REAL ( dp ) :: Bohr_radius = 0.529e-10
  REAL ( dp ) :: gamma_e = 28.0e6
  REAL ( dp ) :: lambda = 2.997
  
  
  common_factor = mu0_by_fpi * mu_n / Bohr_radius ** 3
  
  SELECT CASE ( TRIM ( hfi_output_unit ) )
  CASE ( 'MHz' )
     output_factor = 1e-3 * common_factor * gamma_e
  CASE ( 'mT' )
     output_factor = 1e3  * common_factor
  CASE ( 'G', 'Gauss' )
     output_factor = 1e4  * common_factor
  CASE ( '10e-4cm^-1' )
     output_factor = 1e-3 * common_factor * gamma_e &
          / lambda
  CASE DEFAULT
     CALL errore ( "hyperfine", "unknown units for output", 1 )
  END SELECT
  
  allocate(aux(nrxx))
  allocate(efgg_el(nrxx,3,3))
  allocate(efgg_zora(ngm,nat))
  allocate(efgr_el(nat,3,3))
  allocate(efg_io(nat,3,3))
  allocate(efgr_fc(nat))
  allocate(efg_corr_tens(3,3,nat))
  allocate(efg_total(3,3,nat))
  allocate(fc_recon(nat))
  allocate(fc_recon_zora(nat))
  allocate(efgr_zora(nat))
  
  ! Select majority and minority spin components
  rho_diff = SUM ( rho%of_r( :, 1 ) - rho%of_r( :, nspin ) )
  if ( rho_diff > +1.0d-3 ) then
     s_maj = 1
     s_min = nspin
  else if ( rho_diff < -1.0d-3 ) then
     s_maj = nspin
     s_min = 1
  else
     write ( stdout, * ) "WARNING: rho_diff zero!"
  end if
  
  aux(:) = rho%of_r(:,s_maj) - rho%of_r(:,s_min)
  efgg_el(:,:,:) = (0.0_dp,0.0_dp)
  
  call cft3(aux,nr1,nr2,nr3,nrx1,nrx2,nrx3,-1)
  
  !
  ! calculation of the electric field gradient in the G-space
  !
  do ig= gstart, ngm
     trace = 1.0_dp/3.0_dp * gg(ig)
     
     do alpha=1, 3
        efgg_el(ig,alpha,alpha) = - trace
        
        do beta=1, 3
           efgg_el(ig,alpha,beta) = ( efgg_el(ig,alpha,beta) &
                + g(alpha,ig) * g(beta,ig) ) * fpi * aux(nl(ig)) / gg(ig)
        end do
     end do
  end do
  
  ! 
  ! fourier transform on the atomic position
  !
  
  efgr_el = (0.0_dp,0.0_dp)
  
  do alpha=1,3
     do beta=1,3
        do na=1,nat
           do ig= gstart, ngm
              arg=(tau(1,na)*g(1,ig)+tau(2,na)*g(2,ig)+tau(3,na)*g(3,ig))*tpi
              efgr_el(na,alpha,beta) = efgr_el(na,alpha,beta) &
                   - efgg_el(ig,alpha,beta) * CMPLX(cos(arg),sin(arg))
           end do
        end do
     end do
  end do
  
  !
  ! Fermi contact term
  ! 
  efgr_fc = (0.0_dp,0.0_dp)
  
  do na=1,nat
     do ig = gstart, ngm
        arg = (tau(1,na)*g(1,ig)+tau(2,na)*g(2,ig)+tau(3,na)*g(3,ig))*tpi
        efgr_fc(na) = efgr_fc(na) + aux(nl(ig)) * CMPLX(cos(arg),sin(arg))
     end do
  end do
  
#ifdef __PARA
  call reduce (2*3*3*nat, efgr_el) !2*, efgr_el is a complex array
  call reduce (2*nat, efgr_fc)     !2*, efgr_el is a complex array
#endif
  
  write ( stdout, '( / )' )
  
  do na=1,nat
     do beta=1,3
        fact = hfi_nuclear_g_factor(na) * output_factor
        write ( stdout, 1000 ) atm(ityp(na)), na, "hfi_bare ", &
             ( fact * REAL(efgr_el(na,alpha,beta),dp), alpha = 1, 3 )
     end do
     write ( stdout, '( / )' )
  end do
  
1000 FORMAT(1x,a,i3,2x,a,3(1x,f14.6))
  
  ! Dipole-dipole interaction
  
  CALL efg_correction ( efg_corr_tens )
  
  !
  ! PAW reconstruction for Fermi contact term
  !
  CALL fermi_contact_reconstruction ( fc_recon, fc_recon_zora )
  
  ! 
  ! Fourier transform of Thomson's delta function
  !
  
  CALL delta_Thomson_radial_ft ( delta_t )
  
  ! 
  ! Fourier transform on the atomic position
  !
  
  efgr_zora = (0.0_dp,0.0_dp)
  
  do na=1, nat
     do ig= gstart, ngm
        arg = (tau(1,na)*g(1,ig)+tau(2,na)*g(2,ig)+tau(3,na)*g(3,ig))*tpi
        efgr_zora(na) = efgr_zora(na) &
             + delta_t(ig,ityp(na)) * aux(nl(ig)) * CMPLX(cos(arg),sin(arg))
     end do
  end do
  
  !
  ! symmetrise efg_tensor
  !
  
  do na = 1,nat
     call trntns (efg_corr_tens(:,:,na),at, bg, -1)
  end do
  call symz(efg_corr_tens, nsym, s, nat, irt)
  do na = 1,nat
     call trntns (efg_corr_tens(:,:,na),at, bg, 1)
  end do
  
  !
  ! print results
  !
  
  do na=1,nat
     do beta=1,3
        fact = hfi_nuclear_g_factor(na) * output_factor
        write (stdout,1000) atm(ityp(na)),na,"hfi_recon", &
             ( fact * REAL(efg_corr_tens(alpha,beta,na),dp), alpha = 1, 3 )
     end do
     write ( stdout, '( / )' )
  end do
  
  write ( stdout, '( A, 2/ )' ) &
       " ******************************** HFI ********************************"
  
  do na=1,nat
     
     efg_total(:,:,na) = REAL ( efgr_el(na,:,:) + efg_corr_tens(:,:,na), dp )
     
     do beta=1,3
        fact = hfi_nuclear_g_factor(na) * output_factor
        write (stdout,1000) atm(ityp(na)),na,"hfi",&
             (fact * efg_total(alpha,beta,na),alpha=1,3)
     end do
     write ( stdout, '( / )' )
     
     do alpha=1,3
        do beta=1,3
           work(beta,alpha) = CMPLX ( efg_total(alpha,beta,na), 0.0_dp )
        end do
     end do
     
     !
     ! diagonalise the tensor to extract the parametres
     !
     
     call cdiagh(3,work,3,efg_eig,efg_vect)
     
     v(2) = efg_eig(2)
     if ( abs(efg_eig(1)) < abs(efg_eig(3)) ) then
        v(1) = efg_eig(1)
        v(3) = efg_eig(3)
     else
        v(1) = efg_eig(3)
        v(3) = efg_eig(1)
     end if
     
     
     write ( stdout, '( A, 1x, i3 )' ) atm(ityp(na)), na
     
     write ( stdout, '( A )' ) &
          "Fermi contact term:   bare        reconstruction  total"
     
     common_factor = 8*pi/3 * output_factor * hfi_nuclear_g_factor(ityp(na))
     write ( stdout, '( 14x, 3F16.6, / )' ) &
          common_factor * REAL ( efgr_fc(na) ), &
          common_factor * fc_recon(na), &
          common_factor * ( REAL ( efgr_fc(na) ) + fc_recon(na) )
     
     write ( stdout, '( A )' ) &
          " **** ZORA ****"
     
     common_factor = 8*pi/3 * output_factor * hfi_nuclear_g_factor(ityp(na))
     write ( stdout, '( 14x, 3F16.6, / )' ) &
          common_factor * REAL ( efgr_zora(na) ), &
          common_factor * fc_recon_zora(na), &
          common_factor * ( REAL ( efgr_zora(na) ) + fc_recon_zora(na) )
     
     write ( stdout, 1200 ) atm(ityp(na)), na, &
          hfi_nuclear_g_factor(ityp(na)), hfi_output_unit, &
          v(1:3) * output_factor * hfi_nuclear_g_factor(ityp(na))
     write ( stdout, '( / )' )
  end do
  
1200 FORMAT ( a, 1x, i3, ': g_N', f11.6, 1x, a, 2x, 3f14.6 )
 
end subroutine hyperfine

!******************************************************************************

subroutine efg_correction ( efg_corr_tens )
  
  USE io_files,              ONLY : nwordwfc, iunwfc
  USE kinds,                 ONLY : dp
  USE uspp,                  ONLY : ap
  USE parameters,            ONLY : lmaxx, ntypx
  USE atom,                  ONLY : rgrid
  USE gvect,                 ONLY : g,ngm,ecutwfc
  USE klist,                 ONLY : nks, xk, wk
  USE cell_base,             ONLY : tpiba2
  USE ions_base,             ONLY : nat, ityp, ntyp => nsp
  USE wvfct,                 ONLY : npwx, nbnd, npw, igk, g2kin
  USE wavefunctions_module,  ONLY : evc
  USE paw_gipaw,             ONLY : paw_recon, paw_nkb, paw_vkb, paw_becp
  USE constants,             ONLY : pi
  USE buffers
  USE scf,                   ONLY : rho
  USE lsda_mod,              ONLY : current_spin, nspin, isk
  USE wvfct,                 ONLY : current_k, wg
  USE io_global,             ONLY : stdout
  USE gipaw_module,          ONLY : job, nbnd_occ
  
  implicit none
  
  ! Argument
  real(DP), intent(out) :: efg_corr_tens(3,3,nat)
  
  ! Local
  integer :: j, ill, nt, ibnd, il1, il2, ik, iat, nbs1, nbs2, kkpsi
  integer :: lm,l,m,m1,m2,lm1,lm2, l1, l2,m3,m4,n,n1,nrc
  integer :: ijkb0,ijkb,ih,jh,na,np,ikb,jkb
  INTEGER :: s_min, s_maj, s_weight
  REAL ( dp ) :: rc, rho_diff, sum_occ(nspin)
  complex(DP) :: bec_product
  
  real(DP), allocatable :: at_efg(:,:,:), work(:)
  complex(DP), allocatable :: efg_corr(:,:)
  
  !----------------------------------------------------------------------------
  
  allocate ( efg_corr(lmaxx**2,nat) )
  
  efg_corr = 0.0_dp
  
  ! Select majority and minority spin components
  rho_diff = SUM ( rho%of_r( :, 1 ) - rho%of_r( :, nspin ) )
  if ( rho_diff > +1.0d-3 ) then
     s_maj = 1
     s_min = nspin
  else if ( rho_diff < -1.0d-3 ) then
     s_maj = nspin
     s_min = 1
  else
     write ( stdout, * ) "WARNING: rho_diff zero!"
  end if
  
  allocate ( at_efg ( paw_nkb, paw_nkb, ntypx) ) 
  
  !
  ! calculate radial integration on atom site 
  ! <aephi|1/r^3|aephi>-<psphi|1/r^3|psphi>
  !
  
  at_efg = 0.0_dp
  do nt = 1, ntyp
     
     kkpsi = paw_recon(nt)%aephi(1)%kkpsi
     allocate ( work(kkpsi) )
     
     do il1 = 1, paw_recon(nt)%paw_nbeta
        nrc = paw_recon(nt)%psphi(il1)%label%nrc
        do il2 = 1, paw_recon(nt)%paw_nbeta
           work = 0.0_dp
           do j = 2,nrc
              work(j) = &
                   ( paw_recon(nt)%aephi(il1)%psi(j) &
                   * paw_recon(nt)%aephi(il2)%psi(j) &
                   - paw_recon(nt)%psphi(il1)%psi(j) &
                   * paw_recon(nt)%psphi(il2)%psi(j) ) &
                   / rgrid(nt)%r(j) ** 3
           end do
           call simpson(nrc,work,rgrid(nt)%rab,at_efg(il1,il2,nt))
           !!!print*, nt, il1, il2, at_efg(il1,il2,nt)
        end do
     end do
     
     deallocate ( work )
  end do
  
  !
  !  calculation of the reconstruction part
  !
  
  do ik = 1, nks
     
     current_k = ik
     current_spin = isk(ik)
     
     ! Different sign for spins only in "hyperfine", not "efg"
     if ( current_spin == s_min .AND. job == "hyperfine" ) then
        s_weight = -1
     else
        s_weight = +1
     end if
     
     call gk_sort ( xk(1,ik), ngm, g, ecutwfc / tpiba2, npw, igk, g2kin )
     g2kin(:) = g2kin(:) * tpiba2
     call get_buffer ( evc, nwordwfc, iunwfc, ik)
     
     call init_gipaw_2 ( npw, igk, xk(1,ik), paw_vkb )
     call ccalbec ( paw_nkb, npwx, npw, nbnd, paw_becp, paw_vkb, evc )
     
     do ibnd = 1, nbnd_occ(ik)
        ijkb0 = 0
        do nt = 1, ntyp
           do na = 1, nat
              
              if ( ityp(na) == nt ) then
                 do ih = 1, paw_recon(nt)%paw_nh
                    ikb = ijkb0 + ih
                    nbs1 = paw_recon(nt)%paw_indv(ih)
                    l1 = paw_recon(nt)%paw_nhtol(ih)
                    m1 = paw_recon(nt)%paw_nhtom(ih)
                    lm1 = m1 + l1**2
                    do jh = 1, paw_recon(nt)%paw_nh
                       jkb = ijkb0 + jh
                       nbs2 = paw_recon(nt)%paw_indv(jh)
                       l2 = paw_recon(nt)%paw_nhtol(jh)
                       m2 = paw_recon(nt)%paw_nhtom(jh)
                       lm2 = m2 + l2**2 
                       
                       bec_product = paw_becp(jkb,ibnd) &
                            * CONJG( paw_becp(ikb,ibnd) )

                       do lm = 5, 9
                          efg_corr(lm,na) = efg_corr(lm,na) &
                               + s_weight * bec_product &
                               * at_efg(nbs1,nbs2,nt) &
                               * ap(lm,lm1,lm2) * wg(ibnd,ik)
                       end do
                    end do
                 end do
                 
                 ijkb0 = ijkb0 + paw_recon(nt)%paw_nh
              end if
              
           end do
        end do
     end do
     
  end do
  
  !<apsi> test only
  sum_occ = 0
  DO ik = 1, nks
     current_k = ik
     current_spin = isk(ik)
     sum_occ(current_spin) = sum_occ(current_spin) + SUM(wg(:,ik))
  END DO
  write(6,*) "TMPTMPTMP: ", sum_occ(:nspin)
  !</apsi>
  
  ! For hyper-fine interaction the function used is
  !   (3r_ir_j/r^2 - delta_i,j) / r^3, for efg the other way round;
  ! the former is implemented above, so...
  IF ( job == "hyperfine" ) THEN
     efg_corr = - efg_corr
  END IF
  
  !
  !  transform in cartesian coordinates
  !
  
  efg_corr_tens(1,1,:) =  sqrt(3.0_dp) * efg_corr(8,:) - efg_corr(5,:)
  efg_corr_tens(2,2,:) = -sqrt(3.0_dp) * efg_corr(8,:) - efg_corr(5,:)
  efg_corr_tens(3,3,:) = 2.0_dp * efg_corr(5,:)
  efg_corr_tens(1,2,:) = sqrt(3.0_dp) * efg_corr(9,:)
  efg_corr_tens(2,1,:) = efg_corr_tens(1,2,:)
  efg_corr_tens(1,3,:) = -efg_corr(6,:) * sqrt(3.0_dp)
  efg_corr_tens(3,1,:) = efg_corr_tens(1,3,:)
  efg_corr_tens(2,3,:) = -efg_corr(7,:) * sqrt(3.0_dp)
  efg_corr_tens(3,2,:) = efg_corr_tens(2,3,:)
  
  efg_corr_tens = - sqrt(4.0_dp*pi/5.0_dp)*efg_corr_tens
  
  ! dia_corr(5,:) = 3z^2-1
  ! dia_corr(6,:) = -xz
  ! dia_corr(7,:) = -yz
  ! dia_corr(8,:) = x^2-y^2
  ! dia_corr(9,:) = xy
  
  deallocate ( efg_corr )
  deallocate ( at_efg )
  
end subroutine efg_correction

!******************************************************************************

subroutine fermi_contact_reconstruction ( fc_recon, fc_recon_zora )
  
  USE io_files,              ONLY : nwordwfc, iunwfc
  USE kinds,                 ONLY : dp
  USE uspp,                  ONLY : ap
  USE parameters,            ONLY : lmaxx, ntypx
  USE atom,                  ONLY : rgrid
  USE gvect,                 ONLY : g,ngm,ecutwfc, gg
  USE klist,                 ONLY : nks, xk, wk
  USE cell_base,             ONLY : tpiba2
  USE ions_base,             ONLY : nat, ityp, ntyp => nsp, atm
  USE wvfct,                 ONLY : npwx, nbnd, npw, igk, g2kin
  USE wavefunctions_module,  ONLY : evc
  USE paw_gipaw,             ONLY : paw_recon, paw_nkb, paw_vkb, paw_becp
  USE constants,             ONLY : pi
  USE buffers
  USE scf,                   ONLY : rho
  USE lsda_mod,              ONLY : current_spin, nspin, isk
  USE wvfct,                 ONLY : current_k, wg
  USE io_global,             ONLY : stdout
  USE gipaw_module,          ONLY : job, nbnd_occ, alpha, iverbosity
  USE constants,             ONLY : fpi
  
  implicit none
  
  ! Argument
  real(DP), intent(out) :: fc_recon(nat)
  real(DP), intent(out) :: fc_recon_zora(nat)
  
  ! Local
  integer :: j, ill, nt, ibnd, il1, il2, ik, iat, nbs1, nbs2, kkpsi
  integer :: lm,l,m,m1,m2,lm1,lm2, l1, l2,m3,m4,n,n1,nrc
  integer :: ijkb0,ijkb,ih,jh,na,np,ikb,jkb
  INTEGER :: s_min, s_maj, s_weight, r_first, gv
  REAL ( dp ) :: rc, rho_diff, arg, r_Thomson
  complex(DP) :: bec_product
  
  real(DP), allocatable :: at_efg(:,:,:), at_efg_zora(:,:,:), work(:)
  complex(DP), allocatable :: efg_corr(:,:)
  
  INTEGER, EXTERNAL :: atomic_number
  
  !----------------------------------------------------------------------------
  
  allocate ( efg_corr(lmaxx**2,nat) )
  
  efg_corr = 0.0_dp
  
  ! Select majority and minority spin components
  rho_diff = SUM ( rho%of_r( :, 1 ) - rho%of_r( :, nspin ) )
  if ( rho_diff > +1.0d-3 ) then
     s_maj = 1
     s_min = nspin
  else if ( rho_diff < -1.0d-3 ) then
     s_maj = nspin
     s_min = 1
  else
     write ( stdout, * ) "WARNING: rho_diff zero!"
  end if
  
  allocate ( at_efg(paw_nkb,paw_nkb,ntyp) ) 
  allocate ( at_efg_zora(paw_nkb,paw_nkb,ntyp) ) 
  
  !
  ! calculate radial integration on atom site 
  ! <aephi|1/r^3|aephi>-<psphi|1/r^3|psphi>
  !
  
  at_efg = 0.0_dp
  at_efg_zora = 0.0_dp
  do nt = 1, ntyp
     
     kkpsi = paw_recon(nt)%aephi(1)%kkpsi
     allocate ( work(kkpsi) )
     
     IF ( ABS ( rgrid(nt)%r(1) ) < 1e-8 ) THEN
        r_first = 2
     ELSE
        r_first = 1
     END IF
     
     do il1 = 1, paw_recon(nt)%paw_nbeta
        nrc = paw_recon(nt)%psphi(il1)%label%nrc
        l1 = paw_recon(nt)%psphi(il1)%label%l
        IF ( l1 /= 0 ) CYCLE
        
        do il2 = 1, paw_recon(nt)%paw_nbeta
           
           l2 = paw_recon(nt)%psphi(il2)%label%l
           IF ( l2 /= 0 ) CYCLE
           
           ! Dirac delta function
           work = 0.0_dp
           do j = r_first, nrc
              work(j) = &
                   ( paw_recon(nt)%aephi(il1)%psi(j) &
                   * paw_recon(nt)%aephi(il2)%psi(j) &
                   - paw_recon(nt)%psphi(il1)%psi(j) &
                   * paw_recon(nt)%psphi(il2)%psi(j) ) &
                   / rgrid(nt)%r(j) ** 2 / fpi
           end do
           
           at_efg(il1,il2,nt) = work(2)
           
           IF ( iverbosity > 100 ) THEN
              do j = r_first, nrc
                 write(1010+nt,*) rgrid(nt)%r(j), work(j) !* fpi
              end do
              write(1010+nt,*) ""
           END IF
           
           r_Thomson = atomic_number ( atm(nt) ) * alpha ** 2
           
           ! Dirac delta function
           work = 0.0_dp
           do j = r_first, nrc
              work(j) = fpi * &
                   ( paw_recon(nt)%aephi(il1)%psi(r_first) &
                   * paw_recon(nt)%aephi(il2)%psi(r_first) &
                   - paw_recon(nt)%psphi(il1)%psi(r_first) &
                   * paw_recon(nt)%psphi(il2)%psi(r_first) ) &
                   / rgrid(nt)%r(r_first) ** 2 / fpi &
                   * 2 / ( fpi * rgrid(nt)%r(j) ** 2 * r_Thomson &
                   * ( 1 + 2 * rgrid(nt)%r(j) / r_Thomson ) ** 2 ) &
                   * rgrid(nt)%r(j) ** 2
           end do

!           ! Dirac delta function
!           work = 0.0_dp
!           do j = r_first, nrc
!              work(j) = fpi * &
!                   ( paw_recon(nt)%aephi(il1)%psi(j) &
!                   * paw_recon(nt)%aephi(il2)%psi(j) &
!                   - paw_recon(nt)%psphi(il1)%psi(j) &
!                   * paw_recon(nt)%psphi(il2)%psi(j) ) &
!                   / r(j,nt) ** 2 / fpi &
!                   * 2 / ( fpi * r(j,nt) ** 2 * r_Thomson &
!                   * ( 1 + 2 * r(j,nt) / r_Thomson ) ** 2 ) &
!                   * r(j,nt) ** 2
!           end do
           
!           r_Thomson = 2
!           do j = r_first, nrc
!              work(j) = fpi * &
!                   2.0 / ( fpi * r(j,nt) ** 2 * r_Thomson &
!                   * ( 1.0 + 2.0 * r(j,nt) / r_Thomson ) ** 2 ) &
!                   * r(j,nt) ** 2
!              work(j) = fpi * &
!                   1/(sqrt(2*pi)*r_Thomson) ** 3 &
!                   * exp(-r(j,nt)**2/(2*r_Thomson**2)) * r(j,nt) ** 2
!           end do
           
           IF ( iverbosity > 100 ) THEN
              do j = r_first, nrc
                 write(1000+nt,*) rgrid(nt)%r(j), work(j) !* fpi
              end do
              write(1000+nt,*) ""
           END IF
           
           CALL simpson(nrc,work,rgrid(nt)%rab(:nrc),at_efg_zora(il1,il2,nt))
           
           IF ( iverbosity > 100 ) THEN
              write(6,'(A,2i6,2F10.5)') "DDD: ", &
                   il1, il2, at_efg_zora(il1,il2,nt), &
                   at_efg(il1,il2,nt)
           END IF
           
        end do
     end do
     
     deallocate ( work )
  end do
  
  !
  !  calculation of the reconstruction part
  !
  
  fc_recon = 0.0
  fc_recon_zora = 0.0
  
  do ik = 1, nks
     
     current_k = ik
     current_spin = isk(ik)
     
     ! Different sign for spins only in "hyperfine", not "efg"
     if ( current_spin == s_min .AND. job == "hyperfine" ) then
        s_weight = -1
     else
        s_weight = +1
     end if
     
     call gk_sort ( xk(1,ik), ngm, g, ecutwfc / tpiba2, npw, igk, g2kin )
     call get_buffer ( evc, nwordwfc, iunwfc, ik)
     
     call init_gipaw_2 ( npw, igk, xk(1,ik), paw_vkb )
     call ccalbec ( paw_nkb, npwx, npw, nbnd, paw_becp, paw_vkb, evc )
     
     do ibnd = 1, nbnd_occ(ik)
        ijkb0 = 0
        do nt = 1, ntyp
           do na = 1, nat
              
              if ( ityp(na) == nt ) then
                 do ih = 1, paw_recon(nt)%paw_nh
                    ikb = ijkb0 + ih
                    nbs1 = paw_recon(nt)%paw_indv(ih)
                    l1 = paw_recon(nt)%paw_nhtol(ih)
                    m1 = paw_recon(nt)%paw_nhtom(ih)
                    lm1 = m1 + l1**2
                    IF ( l1 /= 0 ) CYCLE
                    
                    do jh = 1, paw_recon(nt)%paw_nh
                       jkb = ijkb0 + jh
                       nbs2 = paw_recon(nt)%paw_indv(jh)
                       l2 = paw_recon(nt)%paw_nhtol(jh)
                       m2 = paw_recon(nt)%paw_nhtom(jh)
                       lm2 = m2 + l2**2 
                       IF ( l2 /= 0 ) CYCLE
                       
                       bec_product = paw_becp(jkb,ibnd) &
                            * CONJG( paw_becp(ikb,ibnd) )
                       
                       fc_recon(na) = fc_recon(na) &
                            + s_weight * at_efg(nbs1,nbs2,nt) &
                            * bec_product * wg(ibnd,ik)
                       
                       fc_recon_zora(na) = fc_recon_zora(na) &
                            + s_weight * at_efg_zora(nbs1,nbs2,nt) &
                            * bec_product * wg(ibnd,ik)
                       
                    end do
                 end do
                 
                 ijkb0 = ijkb0 + paw_recon(nt)%paw_nh
              end if
              
           end do
        end do
     end do
     
  end do
  
  !
  !  transform in cartesian coordinates
  !
  
  deallocate ( efg_corr )
  deallocate ( at_efg )
  deallocate ( at_efg_zora )
  
end subroutine fermi_contact_reconstruction

!******************************************************************************

subroutine delta_Thomson_radial_ft ( delta_t )
  
  USE kinds,                 ONLY : dp
  USE atom,                  ONLY : rgrid
  USE gvect,                 ONLY : ngm, gg
  USE ions_base,             ONLY : ntyp => nsp, atm
  USE constants,             ONLY : pi, fpi
  USE cell_base,             ONLY : tpiba
  USE gipaw_module,          ONLY : alpha, iverbosity
  
  implicit none
  
  ! Argument
  real(DP), intent(out) :: delta_t(ngm,ntyp)
  
  ! Local
  INTEGER :: gv, j, nt, kkpsi, nrc
  REAL ( dp ) :: gr, r_Thomson
  
  REAL(dp), allocatable :: f_radial(:), work(:)
  
  INTEGER, EXTERNAL :: atomic_number
  
  !----------------------------------------------------------------------------
  
  do nt = 1, ntyp
     
     allocate ( work(rgrid(nt)%mesh), f_radial(rgrid(nt)%mesh) )
     
     ! Thomson's delta function
     
     r_Thomson = atomic_number ( atm(nt) ) * alpha ** 2
     
     ! Terms r(j,nt) ** 2 from the definition of delta_Thomson
     !    and the radial volume element r^2 in integral cancel each other
     DO j = 1, rgrid(nt)%mesh
        f_radial(j) = 2 / ( fpi * r_Thomson &
             * ( 1 + 2 * rgrid(nt)%r(j) / r_Thomson ) ** 2 )
     END DO
     
     DO gv = 1, ngm
        
        ! Thomson delta function
        work = 0.0_dp
        do j = 1, rgrid(nt)%mesh
           gr = SQRT(gg(gv)) * tpiba * rgrid(nt)%r(j)
           IF ( gr < 1.0e-8 ) THEN
              work(j) = f_radial(j) * fpi
           ELSE
              work(j) = f_radial(j) * fpi * SIN ( gr ) / gr
           END IF
        end do
        
        CALL simpson(rgrid(nt)%mesh,work,rgrid(nt)%rab,delta_t(gv,nt))
        
        IF ( iverbosity > 100 ) THEN
           write(1020+nt,*) SQRT(gg(gv))*tpiba, delta_t(gv,nt)
        END IF
        
     END DO
     
     deallocate ( work, f_radial )
  end do
  
end subroutine delta_Thomson_radial_ft
