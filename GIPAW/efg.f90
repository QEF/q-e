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
  complex(DP) :: workc(3,3), efg_vect(3,3)
  
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
             ( 2*sqrt(6.0_dp)*efg_corr_tens(alpha,beta,na), alpha = 1, 3 )
     end do
     write ( stdout, '( / )' )
  end do
  
  write ( stdout, '( "Total EFG calculation:", / )' )
  
  do na=1,nat
     efg_total(:,:,na) = efg_corr_tens(:,:,na) &
          + REAL ( efgr_el(na,:,:) + efg_io(na,:,:), dp )
     
     do beta=1,3
        write (stdout,1000) atm(ityp(na)),na,"efg",&
             (efg_total(alpha,beta,na),alpha=1,3)
     end do
     write ( stdout, '( / )' )
     
     do alpha=1,3
        do beta=1,3
           workc(beta,alpha) = CMPLX ( efg_total(alpha,beta,na), 0.0_dp )
        end do
     end do
     
     !
     ! diagonalise the tensor to extract the quadrupolar parameters Cq and eta
     !
     
     call cdiagh(3,workc,3,efg_eig,efg_vect)
     
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
     
     write ( stdout, '( 1x, a, 1x, i3, 5x, "eig= ", 3(F10.6,3X) )' ) &
          atm(ityp(na)), na, v(1:3)
     
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
                         hfi_isotope, job, hfi_via_reconstruction_only, &
                         iverbosity, radial_integral_splines
  
  implicit none
  
  integer :: alpha, beta, ig, na, i
  INTEGER :: spin, s_min, s_maj
  REAL ( dp ) :: common_factor, output_factor, rho_diff
  real(DP) :: trace, arg, e2, fact
  real(DP) :: efg_eig(3), v(3), delta_Th(ngm,ntypx)
  complex(DP) :: workc(3,3), efg_vect(3,3)
  
  complex(DP), allocatable :: aux(:)
  complex(DP), allocatable :: efgg_el(:,:,:),efgr_el(:,:,:)
  complex(DP), allocatable :: efg_io(:,:,:)
  complex(DP), allocatable :: efgr_fc_bare(:)
  real(DP), allocatable :: efg_corr_tens(:,:,:), efg_total(:,:,:)
  real(DP), allocatable :: fc_recon(:)
  real(DP), allocatable :: fc_recon_extrapolated(:)
  real(DP), allocatable :: fc_recon_zora(:)
  real(DP), allocatable :: fc_recon_zora_extrapolated(:)
  complex(DP), allocatable :: efgr_fc_bare_zora(:)
  
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
  allocate(efgr_el(nat,3,3))
  allocate(efg_io(nat,3,3))
  allocate(efgr_fc_bare(nat))
  allocate(efg_corr_tens(3,3,nat))
  allocate(efg_total(3,3,nat))
  allocate(fc_recon(nat))
  allocate(fc_recon_extrapolated(nat))
  allocate(fc_recon_zora(nat))
  allocate(fc_recon_zora_extrapolated(nat))
  allocate(efgr_fc_bare_zora(nat))
  
  ! Select majority and minority spin components
  rho_diff = SUM ( rho%of_r( :, 1 ) - rho%of_r( :, nspin ) )
#ifdef __PARA
  call reduce(1, rho_diff)
#endif
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
#ifdef __PARA
  call reduce (2*3*3*nat, efgr_el) !2*, efgr_el is a complex array
#endif
  
  write ( stdout, '( / )' )
  
  do na=1,nat
     do beta=1,3
        fact = hfi_nuclear_g_factor(ityp(na)) * output_factor
        write ( stdout, 1000 ) atm(ityp(na)), na, "hfi_bare ", &
             ( fact * REAL(efgr_el(na,alpha,beta),dp), alpha = 1, 3 )
     end do
     write ( stdout, '( / )' )
  end do
  
1000 FORMAT(1x,a,i3,2x,a,3(1x,f14.6))
  
  !
  ! Fermi contact term - non-relativistic bare part
  ! 
  efgr_fc_bare = (0.0_dp,0.0_dp)
  
  IF ( .NOT. hfi_via_reconstruction_only ) THEN
     do na=1,nat
        do ig = gstart, ngm
           arg = (tau(1,na)*g(1,ig)+tau(2,na)*g(2,ig)+tau(3,na)*g(3,ig))*tpi
           efgr_fc_bare(na) = efgr_fc_bare(na) &
                + aux(nl(ig)) * CMPLX(cos(arg),sin(arg))
        end do
     end do
#ifdef __PARA
     call reduce (2*nat, efgr_fc_bare) !2*, efgr_fc_bare is a complex array
#endif
  END IF
  
  !
  ! Dipole-dipole interaction
  !  
  CALL efg_correction ( efg_corr_tens )
  
  !
  ! PAW reconstruction for Fermi contact term
  !
  CALL fermi_contact_reconstruction ( fc_recon, fc_recon_extrapolated, &
       fc_recon_zora, fc_recon_zora_extrapolated )
  
  efgr_fc_bare_zora = (0.0_dp,0.0_dp)
  
  IF ( .NOT. hfi_via_reconstruction_only ) THEN
     
     ! 
     ! Fourier transform of Thomson's delta function
     !
     CALL delta_Thomson_radial_ft ( delta_Th )
     
     ! 
     ! Fourier transform on the atomic position
     !
     do na=1, nat
        do ig= gstart, ngm
           arg = (tau(1,na)*g(1,ig)+tau(2,na)*g(2,ig)+tau(3,na)*g(3,ig))*tpi
           efgr_fc_bare_zora(na) = efgr_fc_bare_zora(na) &
                + delta_Th(ig,ityp(na)) * aux(nl(ig)) &
                * CMPLX ( cos(arg), sin(arg) )
        end do
     end do
#ifdef __PARA
     call reduce (2*nat, efgr_fc_bare_zora) !2*, efgr_fc_bare_zora is a complex array
#endif
  END IF
  
  !
  ! symmetrise efg_tensor (dipole-dipole interaction)
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
        fact = hfi_nuclear_g_factor(ityp(na)) * output_factor
        write (stdout,1000) atm(ityp(na)),na,"hfi_recon", &
             ( fact * efg_corr_tens(alpha,beta,na), alpha = 1, 3 )
     end do
     write ( stdout, '( / )' )
  end do
  
  write ( stdout, '( A, 2/ )' ) &
       " ******************************** HFI ********************************"
  
  do na=1,nat
     
     efg_total(:,:,na) = REAL ( efgr_el(na,:,:), dp ) + efg_corr_tens(:,:,na)
     
     write ( stdout, 1200 ) atm(ityp(na)), na, &
          hfi_nuclear_g_factor(ityp(na)), hfi_output_unit
     do beta=1,3
        fact = hfi_nuclear_g_factor(ityp(na)) * output_factor
        write (stdout,1000) atm(ityp(na)),na,"hfi_dipole",&
             (fact * efg_total(alpha,beta,na),alpha=1,3)
     end do
     
     do alpha=1,3
        do beta=1,3
           workc(beta,alpha) = CMPLX ( efg_total(alpha,beta,na), 0.0_dp )
        end do
     end do
     
     !
     ! diagonalise the tensor to extract the parametres
     !
     
     call cdiagh(3,workc,3,efg_eig,efg_vect)
     
     v(2) = efg_eig(2)
     if ( abs(efg_eig(1)) < abs(efg_eig(3)) ) then
        v(1) = efg_eig(1)
        v(3) = efg_eig(3)
     else
        v(1) = efg_eig(3)
        v(3) = efg_eig(1)
     end if
     
     write ( stdout, '( )' )
     
     write (stdout,1000) atm(ityp(na)),na,"hfi_dipole",&
          v(1:3) * output_factor * hfi_nuclear_g_factor(ityp(na))
     
     write ( stdout, '( )' )
     
     !write ( stdout, '( A, 1x, i3 )' ) atm(ityp(na)), na
     
     write ( stdout, '( A )' ) &
          "Fermi contact term:        bare     reconstruction      total"
     
     IF ( .NOT. radial_integral_splines ) THEN
        WRITE ( stdout, * ) "No extrapolation with Simpson"
     END IF
     
     common_factor = 8*pi/3 * output_factor * hfi_nuclear_g_factor(ityp(na))
     IF ( iverbosity > 5 ) THEN
        write ( stdout, '( 18x, 3(1x,F14.6) )' ) &
             common_factor * REAL ( efgr_fc_bare(na) ), &
             common_factor * fc_recon(na), &
             common_factor * ( REAL ( efgr_fc_bare(na) ) + fc_recon(na) )
     END IF
     
     write ( stdout, '( 18x, 3(1x,F14.6) )' ) &
          common_factor * REAL ( efgr_fc_bare(na) ), &
          common_factor * fc_recon_extrapolated(na), &
          common_factor &
          * ( REAL ( efgr_fc_bare(na) ) + fc_recon_extrapolated(na) )
     
     write ( stdout, '( A )' ) " **** ZORA ****"
     
     common_factor = 8*pi/3 * output_factor * hfi_nuclear_g_factor(ityp(na))
     IF ( iverbosity > 5 ) THEN
        write ( stdout, '( 18x, 3(1x,F14.6) )' ) &
             common_factor * REAL ( efgr_fc_bare_zora(na) ), &
             common_factor * fc_recon_zora(na), &
             common_factor &
             * ( REAL ( efgr_fc_bare_zora(na) ) + fc_recon_zora(na) )
     END IF
     
     write ( stdout, '( 18x, 3(1x,F14.6), / )' ) &
          common_factor * REAL ( efgr_fc_bare_zora(na) ), &
          common_factor * fc_recon_zora_extrapolated(na), &
          common_factor &
          * ( REAL ( efgr_fc_bare_zora(na) ) + fc_recon_zora_extrapolated(na) )
     
     write ( stdout, '( / )' )
     
  end do
  
1200 FORMAT ( 1x, a, i3, ': g_n =', f10.6, 18x, a, 2x, 3f14.6 )
 
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
  USE becmod,                ONLY : calbec
  USE constants,             ONLY : pi, fpi
  USE buffers
  USE scf,                   ONLY : rho
  USE lsda_mod,              ONLY : current_spin, nspin, isk
  USE wvfct,                 ONLY : current_k, wg
  USE io_global,             ONLY : stdout
  USE gipaw_module,          ONLY : job, nbnd_occ, spline_integration, &
                                   radial_integral_splines, &
                                   radial_integral_diamagnetic
  
  implicit none
  
  ! Argument
  real(DP), intent(out) :: efg_corr_tens(3,3,nat)
  
  ! Local
  integer :: j, ill, nt, ibnd, il1, il2, ik, iat, nbs1, nbs2, kkpsi
  integer :: lm,l,m,m1,m2,lm1,lm2, l1, l2,m3,m4,n,n1,nrc
  integer :: ijkb0,ijkb,ih,jh,na,np,ikb,jkb, r_first
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
#ifdef __PARA
    call reduce(1, rho_diff)
#endif
  if ( rho_diff > +1.0d-3 ) then
     s_maj = 1
     s_min = nspin
  else if ( rho_diff < -1.0d-3 ) then
     s_maj = nspin
     s_min = 1
  else
     IF ( nspin > 1 ) THEN
        write ( stdout, * ) "WARNING: rho_diff zero!"
     END IF
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
     
     IF ( ABS ( rgrid(nt)%r(1) ) < 1e-8 ) THEN
        r_first = 2
     ELSE
        r_first = 1
     END IF
     
     do il1 = 1, paw_recon(nt)%paw_nbeta
        nrc = paw_recon(nt)%psphi(il1)%label%nrc
        do il2 = 1, paw_recon(nt)%paw_nbeta
           work = 0.0_dp
           do j = r_first, nrc
              work(j) = &
                   ( paw_recon(nt)%aephi(il1)%psi(j) &
                   * paw_recon(nt)%aephi(il2)%psi(j) &
                   - paw_recon(nt)%psphi(il1)%psi(j) &
                   * paw_recon(nt)%psphi(il2)%psi(j) ) &
                   / rgrid(nt)%r(j) ** 3
           end do
           
           IF ( radial_integral_splines ) THEN
              at_efg(il1,il2,nt) = &
                   spline_integration ( rgrid(nt)%r(:nrc), work(:nrc) )
           ELSE
              call simpson(nrc,work,rgrid(nt)%rab,at_efg(il1,il2,nt))
           END IF
           
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
!     g2kin(:) = g2kin(:) * tpiba2
     call get_buffer ( evc, nwordwfc, iunwfc, ik)
     
     call init_gipaw_2 ( npw, igk, xk(1,ik), paw_vkb )
     !it was: call ccalbec ( paw_nkb, npwx, npw, nbnd, paw_becp, paw_vkb, evc )
     call calbec ( npw, paw_vkb, evc, paw_becp )
     
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

subroutine fermi_contact_reconstruction ( fc_recon, fc_recon_extrapolated, &
     fc_recon_zora, fc_recon_zora_extrapolated )
  
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
  USE becmod,                ONLY : calbec
  USE constants,             ONLY : pi, fpi
  USE buffers
  USE scf,                   ONLY : rho
  USE lsda_mod,              ONLY : current_spin, nspin, isk
  USE wvfct,                 ONLY : current_k, wg
  USE io_global,             ONLY : stdout
  USE gipaw_module,          ONLY : job, nbnd_occ, alpha, iverbosity, &
                                    spline_integration, &
                                    spline_integration_mirror, &
                                    radial_integral_splines, &
                                    radial_extrapolation, &
                                    spline_mirror_extrapolate, &
                                    hfi_extrapolation_npoints, &
                                    hfi_via_reconstruction_only
  
  implicit none
  
  ! Argument
  real(DP), intent(out) :: fc_recon(nat)
  real(DP), intent(out) :: fc_recon_extrapolated(nat)
  real(DP), intent(out) :: fc_recon_zora(nat)
  real(DP), intent(out) :: fc_recon_zora_extrapolated(nat)
  
  ! Local
  integer :: j, ill, nt, ibnd, il1, il2, ik, iat, nbs1, nbs2, kkpsi
  integer :: lm,l,m,m1,m2,lm1,lm2, l1, l2,m3,m4,n,n1,nrc
  integer :: ijkb0,ijkb,ih,jh,na,np,ikb,jkb
  INTEGER :: s_min, s_maj, s_weight, r_first, gv
  REAL ( dp ) :: rc, rho_diff, arg, r_Thomson
  complex(DP) :: bec_product
  
  real(DP), allocatable :: at_efg(:,:,:), at_efg_zora(:,:,:)
  real(DP), allocatable :: at_efg_extrapolated(:,:,:)
  real(DP), allocatable :: at_efg_zora_extrapolated(:,:,:)
  real(DP), allocatable :: work(:)
  REAL ( dp ), ALLOCATABLE :: x_extrapolate(:), y_extrapolate(:)
  
  !<apsi>
  REAL ( dp ) :: startu, startd, x
  real(DP), allocatable :: r_symmetric(:), work_symmetric(:), d2y(:)
  !</apsi>
  
  INTEGER :: norder_extrapolate = 3
  
  INTEGER, EXTERNAL :: atomic_number
  
  !----------------------------------------------------------------------------
  
  allocate ( at_efg(paw_nkb,paw_nkb,ntyp) )
  allocate ( at_efg_extrapolated(paw_nkb,paw_nkb,ntyp) )
  allocate ( at_efg_zora(paw_nkb,paw_nkb,ntyp) )
  allocate ( at_efg_zora_extrapolated(paw_nkb,paw_nkb,ntyp) )
  
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
     
     r_Thomson = atomic_number ( atm(nt) ) * alpha ** 2
     
     ALLOCATE ( x_extrapolate(hfi_extrapolation_npoints) )
     ALLOCATE ( y_extrapolate(hfi_extrapolation_npoints) )
     
     DO j = 1, hfi_extrapolation_npoints
        IF ( radial_integral_splines ) THEN
           x_extrapolate(j) = j / REAL ( hfi_extrapolation_npoints, dp ) &
                * rgrid(nt)%r(r_first)
        ELSE
           x_extrapolate(j) = j / REAL ( hfi_extrapolation_npoints + 1, dp ) &
                * rgrid(nt)%r(r_first)
        END IF
     END DO
     
     do il1 = 1, paw_recon(nt)%paw_nbeta
        nrc = paw_recon(nt)%psphi(il1)%label%nrc
        l1 = paw_recon(nt)%psphi(il1)%label%l
        IF ( l1 /= 0 ) CYCLE
        
        do il2 = 1, paw_recon(nt)%paw_nbeta
           
           l2 = paw_recon(nt)%psphi(il2)%label%l
           IF ( l2 /= 0 ) CYCLE
           
           work = 0.0_dp
           
           IF ( hfi_via_reconstruction_only ) THEN
              do j = r_first, nrc
                 work(j) = &
                      ( paw_recon(nt)%aephi(il1)%psi(j) &
                      * paw_recon(nt)%aephi(il2)%psi(j) ) &
                      / rgrid(nt)%r(j) ** 2 / fpi
              end do
           ELSE
              do j = r_first, nrc
                 work(j) = &
                      ( paw_recon(nt)%aephi(il1)%psi(j) &
                      * paw_recon(nt)%aephi(il2)%psi(j) &
                      - paw_recon(nt)%psphi(il1)%psi(j) &
                      * paw_recon(nt)%psphi(il2)%psi(j) ) &
                      / rgrid(nt)%r(j) ** 2 / fpi
              end do
           END IF
           
           at_efg(il1,il2,nt) = work(r_first)
           
           ! Extrapolation a'la paratec
           x = 0.0
           at_efg_extrapolated(il1,il2,nt) &
                = spline_mirror_extrapolate ( rgrid(nt)%r(:nrc), work(:nrc), x )
           
           CALL radial_extrapolation ( rgrid(nt)%r(:), work(:nrc), &
                x_extrapolate, y_extrapolate, norder_extrapolate )
           
           ! Value at the first extrapolated radial point
           !at_efg_extrapolated(il1,il2,nt) = y_extrapolate(1)
           !at_efg_extrapolated(il1,il2,nt) = y_extrapolate(1)
           
           ! Multiply with Thomson's delta function
           do j = r_first, nrc
              work(j) = work(j) &
                   * 2 / ( r_Thomson &
                   * ( 1 + 2 * rgrid(nt)%r(j) / r_Thomson ) ** 2 )
           end do
           
           IF ( radial_integral_splines ) THEN
              at_efg_zora(il1,il2,nt) &
                   = spline_integration ( rgrid(nt)%r(:nrc), work(:nrc) )
           ELSE
              CALL simpson(nrc,work,rgrid(nt)%rab(:),at_efg_zora(il1,il2,nt))
           END IF
           
           ! ... jetzt wird's extrapoliert...
           
           do j = 1, hfi_extrapolation_npoints
              y_extrapolate(j) = y_extrapolate(j) &
                   * 2 / ( r_Thomson &
                   * ( 1 + 2 * x_extrapolate(j) / r_Thomson ) ** 2 )
           end do
           
           IF ( radial_integral_splines ) THEN
              at_efg_zora_extrapolated(il1,il2,nt) &
                   = at_efg_zora(il1,il2,nt) &
                   + spline_integration_mirror ( x_extrapolate, y_extrapolate )
           ELSE
              at_efg_zora_extrapolated(il1,il2,nt) &
                   = at_efg_zora(il1,il2,nt)
           END IF
           
        end do
     end do
     
     deallocate ( x_extrapolate, y_extrapolate )
     
     deallocate ( work )
  end do
  
  !
  !  calculation of the reconstruction part
  !
  
  ! Select majority and minority spin components
  rho_diff = SUM ( rho%of_r( :, 1 ) - rho%of_r( :, nspin ) )
#ifdef __PARA
    call reduce(1, rho_diff)
#endif
  if ( rho_diff > +1.0d-3 ) then
     s_maj = 1
     s_min = nspin
  else if ( rho_diff < -1.0d-3 ) then
     s_maj = nspin
     s_min = 1
  else
     write ( stdout, * ) "WARNING: rho_diff zero!"
  end if
  
  fc_recon = 0.0
  fc_recon_extrapolated = 0.0
  fc_recon_zora = 0.0
  fc_recon_zora_extrapolated = 0.0
  
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
     !it was: call ccalbec ( paw_nkb, npwx, npw, nbnd, paw_becp, paw_vkb, evc )
     call calbec ( npw, paw_vkb, evc, paw_becp )
     
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
                       
                       fc_recon_extrapolated(na) = fc_recon_extrapolated(na) &
                            + s_weight * at_efg_extrapolated(nbs1,nbs2,nt) &
                            * bec_product * wg(ibnd,ik)
                       
                       fc_recon_zora(na) = fc_recon_zora(na) &
                            + s_weight * at_efg_zora(nbs1,nbs2,nt) &
                            * bec_product * wg(ibnd,ik)
                       
                       fc_recon_zora_extrapolated(na) &
                            = fc_recon_zora_extrapolated(na) &
                            + s_weight &
                            * at_efg_zora_extrapolated(nbs1,nbs2,nt) &
                            * bec_product * wg(ibnd,ik)
                       
                    end do
                 end do
                 
                 ijkb0 = ijkb0 + paw_recon(nt)%paw_nh
              end if
              
           end do
        end do
     end do
     
  end do
  
  deallocate ( at_efg )
  deallocate ( at_efg_extrapolated )
  deallocate ( at_efg_zora )
  deallocate ( at_efg_zora_extrapolated )
  
end subroutine fermi_contact_reconstruction

!******************************************************************************

subroutine delta_Thomson_radial_ft ( delta_Th )
  
  USE kinds,                 ONLY : dp
  USE atom,                  ONLY : rgrid
  USE gvect,                 ONLY : ngm, gg
  USE ions_base,             ONLY : ntyp => nsp, atm
  USE constants,             ONLY : pi, fpi
  USE cell_base,             ONLY : tpiba
  USE gipaw_module,          ONLY : alpha, iverbosity, spline_integration, &
                                    radial_integral_splines
  
  implicit none
  
  ! Argument
  real(DP), intent(out) :: delta_Th(ngm,ntyp)
  
  ! Local
  INTEGER :: gv, j, nt, r_first
  REAL ( dp ) :: gr, r_Thomson, delta_Th_correction
  
  REAL(dp), allocatable :: f_radial(:), work(:)
  
  INTEGER, EXTERNAL :: atomic_number
  
  !----------------------------------------------------------------------------
  
  do nt = 1, ntyp
     
     allocate ( work(rgrid(nt)%mesh), f_radial(rgrid(nt)%mesh) )
     
     ! Thomson's delta function
     
     r_Thomson = atomic_number ( atm(nt) ) * alpha ** 2
     
     ! Terms rgrid(nt)%r(j) ** 2 from the definition of delta_Thomson
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
        
        IF ( radial_integral_splines ) THEN
           delta_Th(gv,nt) &
                = spline_integration ( rgrid(nt)%r(:rgrid(nt)%mesh), &
                work(:rgrid(nt)%mesh) )
        ELSE
           CALL simpson(rgrid(nt)%mesh,work,rgrid(nt)%rab(:),delta_Th(gv,nt))
        END IF
        
        IF ( iverbosity > 100 ) THEN
           write(1020+nt,*) SQRT(gg(gv))*tpiba, delta_Th(gv,nt)
        END IF
        
     END DO
     
     ! Neglect the dependence on sin(gr)/(gr) - r assumed to be small enough
     IF ( ABS ( rgrid(nt)%r(1) ) < 1e-8 ) THEN
        r_first = 2
     ELSE
        r_first = 1
     END IF
     
     delta_Th(:ngm,nt) = delta_Th(:ngm,nt) + delta_Th_correction
     
     deallocate ( work, f_radial )
     
  end do
  
end subroutine delta_Thomson_radial_ft
