! FOR GWW
!
! Copyright (C) 2001-2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! Author: P. Umari
! Modified by G. Stenuit
!
!----------------------------------------------------------------------------
SUBROUTINE energies_xc( lda, n, m, e_xc, e_h )
  !----------------------------------------------------------------------------
  !
  ! computes the expectation values of the exchange and correlation potential
  ! and of hartree potential
  ! ... input:
  ! ...    lda   leading dimension of arrays psi, spsi, hpsi
  ! ...    n     true dimension of psi, spsi, hpsi
  ! ...    m     number of states psi
  ! ...    psi
  !
  ! ... output:
  !       e_xc
  !       e_h
  USE kinds,            ONLY : DP
  USE control_flags,    ONLY : gamma_only
  USE uspp,             ONLY : vkb, nkb
  USE wvfct,            ONLY : igk, g2kin, ecutwfc
  USE fft_base,         ONLY : dffts, dfftp
  USE fft_interfaces,   ONLY : fwfft, invfft
  USE gvecs,          ONLY : nls, doublegrid
  USE gvect,            ONLY : ngm, gstart, nl, nlm, g, gg, gcutm
  USE cell_base,        ONLY : alat, omega
  USE lsda_mod,         ONLY : nspin
  USE ldaU,             ONLY : lda_plus_u
  USE lsda_mod,         ONLY : current_spin
  USE gvect,            ONLY : gstart
  USE io_global,        ONLY : stdout
  USE scf,              ONLY : scf_type, rho, v, vrs, rho_core, rhog_core, vltot
  USE constants,        ONLY : rytoev
  USE io_files,         ONLY : find_free_unit, diropn
  USE ener,             ONLY : etxc, vtxc, ehart

#ifdef EXX
  USE exx,      ONLY : vexx !Suriano
  USE funct,    ONLY : exx_is_active
#endif

  !
  IMPLICIT NONE
  !
  ! ... input/output arguments
  !
  INTEGER          :: lda, n, m
!  COMPLEX(DP) :: psi(lda,m)
  REAL(kind=DP) :: e_xc(m), e_h(m)
  !
  !
  CALL start_clock( 'h_psi' )
  !
  IF ( gamma_only ) THEN
     !
     write(stdout,*) 'BEFORE energies_xc_gamma'
     call flush_unit(stdout)
     CALL energies_xc_gamma()
     write(stdout,*) 'AFTER energies_xc_gamma'
     call flush_unit(stdout)
     !
  ELSE
     !
     CALL energies_xc_k( )
     !
  END IF
  !
  CALL stop_clock( 'h_psi' )
  !
  RETURN
  !
  CONTAINS
     !
     !-----------------------------------------------------------------------
     SUBROUTINE energies_xc_k( )
       !-----------------------------------------------------------------------
       !
       ! ... k-points version
       !
       USE wavefunctions_module, ONLY : psic, evc
       USE becmod,  ONLY : bec_type, becp
       !
       IMPLICIT NONE
       !

       INTEGER :: ibnd, j,is, ig
       !!! REAL(dp) :: etxc,vtxc
       ! counters
       !
       !
       !
       ! ... Here we apply the kinetic energy (k+G)^2 psi
       !
       !
       !
       ! ... Here we add the Hubbard potential times psi
       !
       !
       ! ... the local potential V_Loc psi. First the psi in real space
!set exchange and correlation potential
      !if(.not.allocated(rho%of_r)) write(stdout,*) 'rho not allocated'
      if(.not.allocated(psic)) write(stdout,*) 'psic not allocated'
      !if(.not.allocated(v%of_r)) write(stdout,*) 'v not allocated'
      !
      call v_xc (rho, rho_core, rhog_core, etxc, vtxc, v%of_r)
      !!!! CALL v_xc(rho,rho_core,nr1,nr2,nr3,nr1x,nr2x,nr3x,&
      !!!!&nrxx, nl,ngm,g,nspin,alat,omega,etxc,vtxc,v%of_r)

       do is=1,nspin
          vrs(:,is)=v%of_r(:,is)
          if(doublegrid) call interpolate(vrs(1,is),vrs(1,is),-1)
       enddo
       !
       DO ibnd = 1, m
          !
          CALL start_clock( 'firstfft' )
          !
          psic(1:dffts%nnr) = ( 0.D0, 0.D0 )
          !
          psic(nls(igk(1:n))) = evc(1:n,ibnd)
          !
          CALL invfft ('Wave', psic, dffts)
          !
          CALL stop_clock( 'firstfft' )
          !
          ! ... product with the potential vrs = (vltot+vr) on the smooth grid
          !
          psic(1:dffts%nnr) = psic(1:dffts%nnr) * vrs(1:dffts%nnr,current_spin)
          !
          ! ... back to reciprocal space
          !
          CALL start_clock( 'secondfft' )
          !
          CALL fwfft ('Wave', psic, dffts)
          !
          ! ... addition to the total product
          !
          e_xc(ibnd)=0.d0
          do ig=1,n
             e_xc(ibnd)=e_xc(ibnd)+real(conjg(evc(ig,ibnd))*psic(nls(igk(ig))))
          enddo
          write(stdout,*) 'Routine energies_xc :', ibnd, e_xc(ibnd)*rytoev
!
          CALL stop_clock( 'secondfft' )
          !
       END DO
       !
       !
       !
       RETURN
       !
     END SUBROUTINE energies_xc_k
!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!
     SUBROUTINE energies_xc_gamma

       USE uspp,                 ONLY : okvan
       USE wannier_gw,           ONLY : becp_gw, restart_gww
       USE realus,               ONLY : adduspos_gamma_r
       USE wvfct,                ONLY : npwx, npw, nbnd, et, g2kin
       USE wavefunctions_module, ONLY : evc
       USE klist,                ONLY : xk
       USE mp,                   ONLY : mp_sum
       USE gvect,                ONLY : gstart,g
       USE constants,            ONLY : rytoev
       USE becmod,               ONLY : bec_type, becp, allocate_bec_type, calbec
       USE cell_base,            ONLY : tpiba2


       implicit none

       INTEGER :: ibnd,ir,ig
       INTEGER :: iunwfcreal
       REAL(kind=DP) :: charge
       REAL(kind=DP), ALLOCATABLE :: psi_r(:),psi_rs(:)
       LOGICAL :: exst
       REAL(kind=DP), ALLOCATABLE :: rho_fake_core(:)
       COMPLEX(kind=DP), ALLOCATABLE :: hpsi(:,:), psi(:,:)
       REAL(kind=DP) :: norm, diff
       REAL(kind=DP) :: tresh=1.d-5

!if required calculates also the KS energies
      if(restart_gww==-1) then
         write(stdout,*) 'ATTENZIONE1 and doublegrid=', doublegrid
         write(stdout,*) 'nkb=', nkb, 'nbnd=', nbnd
         call flush_unit(stdout)

         call allocate_bec_type (nkb, nbnd, becp)

         IF ( nkb > 0 )  CALL init_us_2( npw, igk, xk(1,1), vkb )
         call calbec(npw, vkb, evc, becp, nbnd) !! N.B. nbnd is optional !

         allocate(hpsi(npwx,nbnd))
         allocate(psi(npwx,nbnd))
         psi(:,:)=0.D0
         hpsi(:,:)=0.D0

         psi(:,1:nbnd)=evc(:,1:nbnd)

!         g2kin(1:npw) = ( ( xk(1,1) + g(1,igk(1:npw)) )**2 + &
!              ( xk(2,1) + g(2,igk(1:npw)) )**2 + &
!              ( xk(3,1) + g(3,igk(1:npw)) )**2 ) * tpiba2

         call g2_kin(1)

         call h_psi( npwx, npw, nbnd, psi, hpsi )
         et(:,1)=0.d0

         !!! for info:
          write(stdout,*) 'gstart=', gstart, ' and nbnd=', nbnd, &
                          ' and npw=', npw,  ' and npwx=', npwx
         call flush_unit(stdout)
         do ibnd=1,nbnd

           !   call dgemm('T','N',1,1,2*npw,2.d0,evc(:,ibnd),2*npwx,hpsi(:,ibnd),2*npwx,&
           !        &0.d0,et(ibnd,1),1)
           do ig=1,npw
              et(ibnd,1)=et(ibnd,1)+2.d0*dble(conjg(evc(ig,ibnd))*hpsi(ig,ibnd))
           enddo
           if(gstart==2) then
              et(ibnd,1)=et(ibnd,1)-dble(conjg(evc(1,ibnd))*hpsi(1,ibnd))
           endif
         enddo

         call mp_sum(et(:,1))

         do ibnd=1,nbnd
            write(stdout,*) 'KS energy:', ibnd, et(ibnd,1)*rytoev
         enddo
         call flush_unit(stdout)
         deallocate(hpsi,psi,becp%r)
       endif
       !
       allocate(psi_r(dfftp%nnr),psi_rs(dffts%nnr))
       !
       iunwfcreal=find_free_unit()
       CALL diropn( iunwfcreal, 'real_whole', dffts%nnr, exst )
       !
!calculate xc potential on fine grid
       !
       !if(.not.allocated(rho%of_r)) write(stdout,*) 'rho not allocated'
       !if(.not.allocated(v%of_r)) write(stdout,*) 'v not allocated'
       allocate(rho_fake_core(dfftp%nnr))
       rho_fake_core(:)=0.d0
       CALL v_xc( rho, rho_core, rhog_core, etxc, vtxc, v%of_r )
       !!CALL v_xc(rho,rho_core,nr1,nr2,nr3,nr1x,nr2x,nr3x,&
       !!     &nrxx, nl,ngm,g,nspin,alat,omega,etxc,vtxc,v%of_r)
       !
       !!! to test non-linear core correction :
       !! replace rho_core by rho_fake_core in
       !! CALL v_xc( rho, rho_core, rhog_core, etxc, vtxc, v%of_r )
       !
       !  CALL v_xc(rho,rho_fake_core,nr1,nr2,nr3,nr1x,nr2x,nr3x,&
       !    &nrxx, nl,ngm,g,nspin,alat,omega,etxc,vtxc,vr)
       deallocate(rho_fake_core)
       !
       !
       do ibnd=1,m!loop on states
         !read from disk wfc on coarse grid
         !
         CALL davcio( psi_rs, dffts%nnr, iunwfcreal, ibnd, -1)
         !
         !
         if(doublegrid) then
           call interpolate(psi_r,psi_rs,1)
         else
           psi_r(:)=psi_rs(:)
         endif
         !
         !
         do ir=1,dfftp%nnr
            psi_r(ir)=psi_r(ir)**2.d0
         enddo
         !
         !
         if(okvan) call adduspos_gamma_r(ibnd,ibnd,psi_r,1,becp_gw(:,ibnd),becp_gw(:,ibnd))
         !
         !
         ! CHECK if the norm if equal to 1.0
         !
         norm=0.0
         do ir=1,dfftp%nnr
            norm=norm+psi_r(ir)
         enddo
         norm=norm/dble(dfftp%nr1*dfftp%nr2*dfftp%nr3)
         call  mp_sum(norm)
         !
         diff=abs( norm - 1.0d0 )
         !
         ! Write a warning if the norm is not equal to one (within a given threshold)
         if(diff .gt. tresh ) then
            write(stdout,*) 'WARNING: in e_xc part : for the band i=', ibnd
            write(stdout,*) 'WARNING: norm after US contribution is not exactly equal to 1.0 ', norm
            call flush_unit(stdout)
         endif
         !
         e_xc(ibnd)=0.d0
         do ir=1,dfftp%nnr
           e_xc(ibnd)=e_xc(ibnd)+psi_r(ir)*v%of_r(ir,1)!the 1 is for the spin NOT IMPLEMENTED YET
         enddo
         !
         e_xc(ibnd)=e_xc(ibnd)/dble(dfftp%nr1*dfftp%nr2*dfftp%nr3)
         !
         call mp_sum(e_xc(ibnd))
         !
         write(stdout,*) 'Routine energies_xc :', ibnd, e_xc(ibnd)*rytoev
         call flush_unit(stdout)
         !
       enddo
       !
       v%of_r(:,:)=0.d0
       !
       CALL v_h( rho%of_g, ehart, charge, v%of_r )
       !!CALL v_h( rho, nr1, nr2, nr3, nr1x, nr2x, nr3x, nrxx, nl, &
       !!         ngm, gg, gstart, nspin, alat, omega, ehart, charge, v%of_r )
       !
       do ibnd=1,m!loop on states
           !read from disk wfc on coarse grid
           CALL davcio( psi_rs,dffts%nnr,iunwfcreal,ibnd,-1)
           if(doublegrid) then
              call interpolate(psi_r,psi_rs,1)
           else
              psi_r(:)=psi_rs(:)
           endif
           !
           do ir=1,dfftp%nnr
              psi_r(ir)=psi_r(ir)**2.d0
           enddo
           !
           if(okvan) call adduspos_gamma_r(ibnd,ibnd,psi_r,1,becp_gw(:,ibnd),becp_gw(:,ibnd))
           !
           ! CHECK the norm
           norm=0.0
           do ir=1,dfftp%nnr
              norm=norm+psi_r(ir)
           enddo
           norm=norm/dble(dfftp%nr1*dfftp%nr2*dfftp%nr3)
           call  mp_sum(norm)
           !
           diff=abs( norm - 1.0d0 )
           !
           ! Write a warning if the norm is not equal to one (within a given threshold)
           if(diff .gt. tresh ) then
              write(stdout,*) 'WARNING: in e_h part : for the band i=', ibnd
              write(stdout,*) 'WARNING: norm after US contribution is not exactly equal to 1.0 ', norm
              call flush_unit(stdout)
           endif
           !
           e_h(ibnd)=0.d0
           do ir=1,dfftp%nnr
              e_h(ibnd)=e_h(ibnd)+psi_r(ir)*v%of_r(ir,1)!the 1 is for the spin NOT IMPLEMENTED YET
           enddo
           e_h(ibnd)=e_h(ibnd)/dble(dfftp%nr1*dfftp%nr2*dfftp%nr3)
           !
           call mp_sum(e_h(ibnd))
           write(stdout,*) 'Routine energies_h :', ibnd, e_h(ibnd)*rytoev
           !
       enddo
       !
       deallocate(psi_r,psi_rs)
       !
       close(iunwfcreal)
       !
       return
       !
  END SUBROUTINE energies_xc_gamma
     !
     !
END SUBROUTINE energies_xc


SUBROUTINE write_energies_xc(e_xc)
  !
  !
  USE kinds, ONLY : DP
  USE wannier_gw, ONLY : num_nbnds, nbnd_normal
  USE io_files, ONLY : find_free_unit, prefix
  USE io_global, ONLY : ionode
  USE wvfct,    ONLY : nbnd

  implicit none

  REAL(kind=DP), intent(in) :: e_xc(nbnd)!exchange and correlation energies

  INTEGER :: iunu, iw


  if(ionode) then
     iunu = find_free_unit()

     open(unit=iunu,file=trim(prefix)//'.dft_xc',status='unknown',form='unformatted')

     write(iunu) nbnd_normal


     do iw=1,nbnd_normal
        write(iunu) e_xc(iw)
     enddo

     close(iunu)
  endif
  !
  !
  return

END SUBROUTINE write_energies_xc
