!
! Copyright (C) 2001-2005 Quantum-ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "f_defs.h"
!-----------------------------------------------------------------------

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
  USE gipaw_module, ONLY : q_efg
  
  implicit none
  
  real(DP) :: eta, Cq
  real(DP) :: fac, trace, arg, e2
  integer :: alpha, beta, ig, na, i
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
  aux(:) = rho(:,1)
  efgg_el(:,:,:) = (0.0_dp,0.0_dp)
  
  call cft3(aux,nr1,nr2,nr3,nrx1,nrx2,nrx3,-1)
  
  !
  ! calculation of the electic field gradient in the G-space
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
              efgr_el(na,alpha,beta)=efgr_el(na,alpha,beta)+ &
                   efgg_el(ig,alpha,beta) * CMPLX(cos(arg),sin(arg))
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

subroutine efg_correction ( efg_corr_tens )
  
  USE io_files ,ONLY : nwordwfc, iunwfc
  USE kinds , ONLY : dp
  USE uspp ,ONLY : ap
  USE parameters, ONLY : lmaxx, ntypx
  USE atom , ONLY : r,rab,msh
  USE gvect, ONLY : g,ngm,ecutwfc
  USE klist, ONLY : nks, xk, wk
  USE cell_base, ONLY : tpiba2
  USE ions_base, ONLY : nat, ityp, ntyp => nsp
  USE wvfct, ONLY :npwx, nbnd, npw, igk, g2kin
  USE wavefunctions_module, ONLY : evc
  USE paw, ONLY : paw_vkb, paw_becp, paw_nkb, aephi, psphi, &
       paw_nh, paw_nhtol, paw_nhtom, paw_indv, paw_nbeta
  USE constants, ONLY : pi
  
  implicit none
  
  ! Argument
  real(DP), intent(out) :: efg_corr_tens(3,3,nat)
  
  ! Local
  integer :: j, ill, nt, ibnd, il1, il2, ik, iat, nbs1, nbs2, kkpsi
  integer :: lm,l,m,m1,m2,lm1,lm2, l1, l2,m3,m4,n,n1,nrc
  integer :: ijkb0,ijkb,ih,jh,na,np,ikb,jkb
  real(DP), allocatable :: at_efg(:,:,:), work(:)
  complex(DP) :: bec_product
  complex(DP), allocatable :: efg_corr(:,:)
  real(DP) :: rc
  
  !----------------------------------------------------------------------------
  
  allocate ( efg_corr(lmaxx**2,nat) )
  
  efg_corr = 0.0_dp
  
  call init_paw_1
  
  allocate ( at_efg ( paw_nkb, paw_nkb, ntypx) ) 
  allocate ( paw_vkb ( npwx,  paw_nkb ) )
  allocate ( paw_becp ( paw_nkb, nbnd ) )
  
  !
  ! calculate radial integration on atom site 
  ! <aephi|1/r^3|aephi>-<psphi|1/r^3|psphi>
  !
  
  at_efg = 0.0_dp
  do nt = 1, ntyp
     
     kkpsi = aephi(nt,1)%kkpsi
     allocate (work(kkpsi))
     
     do il1 = 1, paw_nbeta(nt)
        nrc = psphi(nt,il1)%label%nrc
        do il2 = 1, paw_nbeta(nt)
           work = 0.0_dp
           do j = 2,nrc
              work(j)=(aephi(nt,il1)%psi(j)*aephi(nt,il2)%psi(j)-&
                   psphi(nt,il1)%psi(j)*psphi(nt,il2)%psi(j))/r(j,nt)**3
           end do
           call simpson(nrc,work,rab(:,nt),at_efg(il1,il2,nt))
        end do
     end do
     
     deallocate ( work )
  end do
  
  !
  !  calculation of the reconstruction part
  !
  
  do ik = 1, nks
     
     call gk_sort ( xk(1,ik), ngm, g, ecutwfc / tpiba2, npw, igk, g2kin )
     call davcio ( evc, nwordwfc, iunwfc, ik, -1 )
     
     call init_paw_2 ( npw, igk, xk(1,ik), paw_vkb )
     call ccalbec ( paw_nkb, npwx, npw, nbnd, paw_becp, paw_vkb, evc )
     
     do ibnd = 1, nbnd
        ijkb0 = 0
        do nt = 1, ntyp
           do na = 1, nat
              
              if ( ityp(na) == nt ) then
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
                            * CONJG( paw_becp(ikb,ibnd) )
                       
                       do lm = 5, 9
                          efg_corr(lm,na) = efg_corr(lm,na) &
                               + bec_product &
                               * at_efg(nbs1,nbs2,nt) &
                               * ap(lm,lm1,lm2) * wk(ik) / 2.0_dp
                       end do
                    end do
                 end do
                 
                 ijkb0 = ijkb0 + paw_nh (nt)
              end if
              
           end do
        end do
     end do
     
  end do
  
  !
  !  transforme in cartesian coordinates
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
  
  efg_corr_tens = -sqrt(4.0_dp*pi/5.0_dp)*efg_corr_tens
  
  ! dia_corr(5,:) = 3z^2-1
  ! dia_corr(6,:) = -xz
  ! dia_corr(7,:) = -yz
  ! dia_corr(8,:) = x^2-y^2
  ! dia_corr(9,:) = xy
  
  deallocate ( efg_corr )
  deallocate ( at_efg )
  deallocate ( paw_vkb )
  
end subroutine efg_correction
