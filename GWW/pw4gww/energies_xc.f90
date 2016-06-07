!
! Copyright (C) 2001-2013 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!

!
!----------------------------------------------------------------------------
SUBROUTINE energies_xc( lda, n, m, psi, e_xc, e_h,ispin )
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
  USE kinds,    ONLY : DP
  USE uspp,     ONLY : vkb, nkb
  USE gvecs,  ONLY : nls, doublegrid
  USE gvect,                ONLY : ngm, gstart, nl, nlm, g, gg, gcutm
  USE cell_base,            ONLY :  alat, omega
  USE lsda_mod,             ONLY : nspin
  USE ldaU,     ONLY : lda_plus_u
  USE lsda_mod, ONLY : current_spin
  USE gvect,    ONLY : gstart
  USE io_global, ONLY :stdout
  USE scf,       ONLY : rho, vltot, vrs, rho_core,rhog_core, scf_type
  USE constants,  ONLY :rytoev
  USE io_files, ONLY : diropn
  USE mp, ONLY : mp_sum, mp_barrier
  USE mp_world, ONLY : world_comm
  USE control_flags,        ONLY : gamma_only
  USE funct,            ONLY : dft_is_meta
  USE fft_base,             ONLY : dfftp, dffts
  USE fft_interfaces,       ONLY : fwfft, invfft

  USE exx,      ONLY : vexx !Suriano
  USE funct,    ONLY : exx_is_active,dft_is_hybrid
  USE klist, ONLY : igk_k

  !
  IMPLICIT NONE
  INTEGER, EXTERNAL :: find_free_unit
  !
  ! ... input/output arguments
  !
  INTEGER          :: lda, n, m
  COMPLEX(DP) :: psi(lda,m)
  REAL(kind=DP) :: e_xc(m), e_h(m)
  INTEGER, INTENT(in) :: ispin !spin 1,2

  REAL(kind=DP), ALLOCATABLE :: vr(:,:)
  !
  !
  CALL start_clock( 'h_psi' )
  allocate(vr(dfftp%nnr,nspin))
  !
  IF ( gamma_only ) THEN
     !
     CALL energies_xc_gamma()
     !
  ELSE
     !
     CALL energies_xc_k( )
     !
  END IF

  !
  CALL stop_clock( 'h_psi' )
  deallocate(vr)
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
       USE wavefunctions_module, ONLY : psic
       USE becmod,  ONLY : becp
       !
       IMPLICIT NONE
       !

       INTEGER :: ibnd, j,is, ig
       REAL(dp) :: etxc,vtxc
       REAL(kind=DP) :: ehart, charge
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
          if(.not.allocated(psic)) write(stdout,*) 'psic not allocated'
       if (dft_is_meta()) then
!         call v_xc_meta( rho, rho_core, rhog_core, etxc, vtxc, v%of_r, v%kin_r )
      else
         CALL v_xc( rho, rho_core, rhog_core, etxc, vtxc, vr )
      endif


       do is=1,nspin
          vrs(:,is)=vr(:,is)
          if(doublegrid) call interpolate(vrs(1,is),vrs(1,is),-1)
       enddo
       !
       DO ibnd = 1, m
          !
          CALL start_clock( 'firstfft' )
          !
          psic(1:dffts%nnr) = ( 0.D0, 0.D0 )
          !
          psic(nls(igk_k(1:n,1))) = psi(1:n,ibnd)
          !
          CALL invfft ('Wave', psic, dffts)

          !
          CALL stop_clock( 'firstfft' )
          !
          ! ... product with the potential vrs = (vltot+vr) on the smooth grid
          !
          psic(1:dffts%nnr) = psic(1:dffts%nnr) * vrs(1:dffts%nnr,1)
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
             e_xc(ibnd)=e_xc(ibnd)+real(conjg(psi(ig,ibnd))*psic(nls(igk_k(ig,1))))
          enddo
          call mp_sum(e_xc(ibnd),world_comm)
          write(stdout,*) 'energies_xc :', ibnd, e_xc(ibnd)*rytoev
!
          CALL stop_clock( 'secondfft' )
          !
       END DO
       vr(:,:)=0.d0
       call  v_h(rho%of_g , ehart, charge, vr )
       do is=1,nspin
          vrs(:,is)=vr(:,is)
          if(doublegrid) call interpolate(vrs(1,is),vrs(1,is),-1)
       enddo


       DO ibnd = 1, m

          CALL start_clock( 'firstfft' )
          psic(1:dffts%nnr) = ( 0.D0, 0.D0 )
          psic(nls(igk_k(1:n,1))) = psi(1:n,ibnd)

          CALL invfft ('Wave', psic, dffts)


          CALL stop_clock( 'firstfft' )


          psic(1:dffts%nnr) = psic(1:dffts%nnr) * vrs(1:dffts%nnr,1)

          CALL start_clock( 'secondfft' )
          CALL fwfft ('Wave', psic, dffts)
          e_h(ibnd)=0.d0
          do ig=1,n
             e_h(ibnd)=e_h(ibnd)+real(conjg(psi(ig,ibnd))*psic(nls(igk_k(ig,1))))
          enddo
          call mp_sum(e_h(ibnd),world_comm)
          write(stdout,*) 'energies_h :', ibnd, e_h(ibnd)*rytoev

          CALL stop_clock( 'secondfft' )

       enddo!

       !
       !
       !
       RETURN
       !
     END SUBROUTINE energies_xc_k
     SUBROUTINE energies_xc_gamma

       USE uspp, ONLY : okvan
       USE wannier_gw, ONLY : becp_gw, restart_gww,l_whole_s,l_verbose,&
                               &l_scissor,scissor,num_nbndv,num_nbnds
      ! USE realus,  ONLY : adduspos_gamma_r
       USE wvfct,    ONLY : npwx,npw,nbnd, et,g2kin
       USE wavefunctions_module, ONLY : evc
       USE klist,                ONLY : xk
       USE mp, ONLY : mp_sum
       USE mp_world, ONLY : world_comm
       USE gvect,  ONLY : gstart,g
       USE constants, ONLY : rytoev
       USE becmod,           ONLY : becp, calbec,allocate_bec_type,deallocate_bec_type
       USE cell_base,            ONLY : tpiba2
       USE io_global, ONLY : ionode
       USE io_files, ONLY :prefix,tmp_dir
     USE exx, ONLY : exxalfa

       implicit none

       INTEGER, EXTERNAL :: find_free_unit

       INTEGER :: ibnd,jbnd,ir,ig
       INTEGER :: iunwfcreal,iunu
       REAL(kind=DP) :: etxc,vtxc,ehart, charge
       REAL(kind=DP), ALLOCATABLE :: psi_r(:),psi_rs(:)
       LOGICAL :: exst
       REAL(kind=DP), ALLOCATABLE :: rho_fake_core(:)
       COMPLEX(kind=DP), ALLOCATABLE :: hpsi(:,:)
       REAL(kind=DP), ALLOCATABLE :: exact_x(:)
       REAL(kind=DP), ALLOCATABLE :: e_hub(:)!Hubbard contribution to KS energies
       REAL(kind=DP), ALLOCATABLE :: et_off(:,:)!off-diagonal energies


       allocate(exact_x(nbnd))
       allocate(e_hub(nbnd))
       if(l_whole_s) then
          allocate (et_off(nbnd,nbnd))
       endif
!if required calculates also the KS energies
 !     if(restart_gww==-1) then
         if(l_verbose) write(stdout,*) 'ATTENZIONE1'
         FLUSH(stdout)
         !allocate( becp%r( nkb, nbnd ) )
         call allocate_bec_type ( nkb, nbnd, becp)
         if(l_verbose) write(stdout,*) 'ATTENZIONE2'
         FLUSH(stdout)

         IF ( nkb > 0 )  CALL init_us_2( npw, igk_k(1,1), xk(1,1), vkb )
         !call ccalbec( nkb, npwx, npw, nbnd, becp%r, vkb, evc )
        !if(nkb> 0)CALL calbec ( npw, vkb, psi, becp, nbnd)

        if(l_verbose)write(stdout,*) 'ATTENZIONE3'
         FLUSH(stdout)

         allocate(hpsi(npwx,nbnd))
         if(l_verbose)write(stdout,*) 'ATTENZIONE4'
         FLUSH(stdout)

         g2kin(1:npw) = ( ( xk(1,1) + g(1,igk_k(1:npw,1)) )**2 + &
              ( xk(2,1) + g(2,igk_k(1:npw,1)) )**2 + &
              ( xk(3,1) + g(3,igk_k(1:npw,1)) )**2 ) * tpiba2


         if(l_verbose)write(stdout,*) 'ATTENZIONE5'
          FLUSH(stdout)


 !         exxalfa=0.d0!ATTENZIONE
         call h_psi( npwx, npw, nbnd, psi, hpsi )
         et(:,ispin)=0.d0
         if(l_verbose)write(stdout,*) 'ATTENZIONE6'
         if(l_verbose)write(stdout,*) 'EXXALFA', exxalfa
         FLUSH(stdout)

          do ibnd=1,nbnd

           !   call dgemm('T','N',1,1,2*npw,2.d0,evc(:,ibnd),2*npwx,hpsi(:,ibnd),2*npwx,&
           !        &0.d0,et(ibnd,1),1)
              do ig=1,npw
                 et(ibnd,ispin)=et(ibnd,ispin)+2.d0*dble(conjg(evc(ig,ibnd))*hpsi(ig,ibnd))
              enddo
              if(gstart==2) then
                et(ibnd,ispin)=et(ibnd,ispin)-dble(conjg(evc(1,ibnd))*hpsi(1,ibnd))
             endif
          enddo
          call mp_sum(et(:,ispin),world_comm)
          if(l_scissor) then
             et(1:num_nbndv(ispin),ispin)=et(1:num_nbndv(ispin),ispin)+scissor(1)/rytoev
             et(num_nbndv(ispin)+1:num_nbnds,ispin)=et(num_nbndv(ispin)+1:num_nbnds,ispin)+scissor(2)/rytoev
          endif

          if(l_verbose)write(stdout,*) 'ATTENZIONE7'
          FLUSH(stdout)
!if required calculate Hubbard U contribution to eigen-energies
         e_hub(:)=0.d0
         if ( lda_plus_u ) then
            hpsi(:,:)=(0.d0,0.d0)
            CALL vhpsi( npwx, npw, nbnd, psi, hpsi )

            do ibnd=1,nbnd

               do ig=1,npw
                  e_hub(ibnd)=e_hub(ibnd)+2.d0*dble(conjg(psi(ig,ibnd))*hpsi(ig,ibnd))
               enddo
               if(gstart==2) then
                  e_hub(ibnd)=e_hub(ibnd)-dble(conjg(psi(1,ibnd))*hpsi(1,ibnd))
               endif
            enddo
            call mp_sum(e_hub(:),world_comm)
            do ibnd=1,nbnd
               write(stdout,*) 'Hubbard U energy:',ibnd,e_hub(ibnd)*rytoev
            enddo
            FLUSH(stdout)

         endif
          do ibnd=1,nbnd
             write(stdout,*) 'KS energy:',ibnd,et(ibnd,ispin)*rytoev
          enddo
          FLUSH(stdout)

!in case of hybrid functionals and HF we have to calculated also the exact exchange part

          if(dft_is_hybrid()) then
!NOT_TO_BE_INCLUDED_START
             hpsi(:,:)=(0.d0,0.d0)
             call vexx( npwx, npw, nbnd, psi, hpsi )
             do ibnd=1,nbnd
                call dgemm('T','N',1,1,2*npw,2.d0,evc(:,ibnd),2*npwx,hpsi(:,ibnd),2*npwx,&
                     &0.d0,exact_x(ibnd),1)
                if(gstart==2) then
                   exact_x(ibnd)=exact_x(ibnd)-dble(conjg(evc(1,ibnd))*hpsi(1,ibnd))
                endif
                call mp_sum(exact_x(ibnd),world_comm)
                write(stdout,*) 'Exact exchange :',ibnd, exact_x(ibnd)
             enddo
!NOT_TO_BE_INCLUDED_END
          end if

 !         deallocate(hpsi,becp%r)
          call deallocate_bec_type ( becp)
!       endif



       allocate(psi_r(dfftp%nnr),psi_rs(dfftp%nnr))

       iunwfcreal=find_free_unit()
       CALL diropn( iunwfcreal, 'real_whole', dffts%nnr, exst )


!calculate xc potential on fine grid


       if(.not.allocated(vr)) write(stdout,*) 'vr not allocated'
       allocate(rho_fake_core(dfftp%nnr))
       rho_fake_core(:)=0.d0

       if (dft_is_meta()) then
      !    call v_xc_meta( rho, rho_core, rhog_core, etxc, vtxc, v%of_r, v%kin_r )
       else
          CALL v_xc( rho, rho_core, rhog_core, etxc, vtxc, vr )
       endif


     deallocate(rho_fake_core)


       if(l_whole_s) then
!NOT_TO_BE_INCLUDED_START
          allocate(hpsi(npwx,nbnd))
          hpsi(:,:)=(0.d0,0.d0)
          CALL vloc_psi_gamma ( npwx, npw, nbnd, evc, vr(1,ispin), hpsi )
          call dgemm('T','N',nbnd,nbnd,2*npw,2.d0,evc,2*npwx,hpsi,2*npwx,&
               &0.d0,et_off,nbnd)
          if(gstart==2) then
             do ibnd=1,nbnd
                do jbnd=1,nbnd
                   et_off(ibnd,jbnd)=et_off(ibnd,jbnd)-dble(conjg(evc(1,ibnd))*hpsi(1,jbnd))
                enddo
             enddo
          endif
          deallocate(hpsi)
          call mp_sum(et_off,world_comm)
!write on file
          if(ionode) then
             iunu = find_free_unit()
             if(ispin==1) then
                open(unit=iunu,file=trim(tmp_dir)//trim(prefix)//'.exc_off',status='unknown',form='unformatted')
             else
                open(unit=iunu,file=trim(tmp_dir)//trim(prefix)//'.exc_off2',status='unknown',form='unformatted')
             endif
             write(iunu) nbnd
             do ibnd=1,nbnd
                write(iunu) et_off(1:nbnd,ibnd)
             enddo
             close(iunu)
          endif
!NOT_TO_BE_INCLUDED_END
       endif




       do ibnd=1,m!loop on states
!read from disk wfc on coarse grid
         CALL davcio( psi_rs,dffts%nnr,iunwfcreal,ibnd+(ispin-1)*nbnd,-1)
         if(doublegrid) then
           call interpolate(psi_r,psi_rs,1)
         else
           psi_r(:)=psi_rs(:)
         endif

         do ir=1,dfftp%nnr
            psi_r(ir)=psi_r(ir)**2.d0
         enddo

         !if(okvan) call adduspos_gamma_r(ibnd,ibnd,psi_r,1,becp_gw(:,ibnd),becp_gw(:,ibnd))

         e_xc(ibnd)=0.d0
         do ir=1,dfftp%nnr
           e_xc(ibnd)=e_xc(ibnd)+psi_r(ir)*vr(ir,ispin)!the 1 is for the spin NOT IMPLEMENTED YET
         enddo
         e_xc(ibnd)=e_xc(ibnd)/dble(dfftp%nr1*dfftp%nr2*dfftp%nr3)

         call mp_sum(e_xc(ibnd),world_comm)

!ifrequired add the contribution from exact exchange for hybrids and HF
         if(dft_is_hybrid()) then
!NOT_TO_BE_INCLUDED_START
            e_xc(ibnd)=e_xc(ibnd)+exact_x(ibnd)
!NOT_TO_BE_INCLUDED_END
         endif

         write(stdout,*) 'Routine energies_xc :', ibnd, e_xc(ibnd)*rytoev

!now hartree term


      enddo

!if required add to e_xc Hubbard U terms

      if(lda_plus_u) then

         e_xc(1:nbnd)=e_xc(1:nbnd)+e_hub(1:nbnd)
      endif


       vr(:,:)=0.d0

       call  v_h(rho%of_g , ehart, charge, vr )


        do ibnd=1,m!loop on states
!read from disk wfc on coarse grid
           CALL davcio( psi_rs,dffts%nnr,iunwfcreal,ibnd+(ispin-1)*nbnd,-1)
           if(doublegrid) then
              call interpolate(psi_r,psi_rs,1)
           else
              psi_r(:)=psi_rs(:)
           endif

           do ir=1,dfftp%nnr
              psi_r(ir)=psi_r(ir)**2.d0
           enddo

           !if(okvan) call adduspos_gamma_r(ibnd,ibnd,psi_r,1,becp_gw(:,ibnd),becp_gw(:,ibnd))

           e_h(ibnd)=0.d0
           do ir=1,dfftp%nnr
              e_h(ibnd)=e_h(ibnd)+psi_r(ir)*vr(ir,ispin)
           enddo
           e_h(ibnd)=e_h(ibnd)/dble(dfftp%nr1*dfftp%nr2*dfftp%nr3)

           call mp_sum(e_h(ibnd),world_comm)
           write(stdout,*) 'Routine energies_h :', ibnd, e_h(ibnd)*rytoev

!now hartree term


        enddo





       deallocate(psi_r,psi_rs)
       deallocate(exact_x)
      close(iunwfcreal)
      deallocate(e_hub)
      if(l_whole_s) then
!NOT_TO_BE_INCLUDED_START
         deallocate(et_off)
!NOT_TO_BE_INCLUDED_END
      endif

      return

     END SUBROUTINE energies_xc_gamma
     !
END SUBROUTINE energies_xc


SUBROUTINE write_energies_xc(e_xc)

  USE kinds, ONLY : DP
  USE wannier_gw, ONLY : num_nbnds, l_verbose
  USE io_files, ONLY : prefix,tmp_dir
  USE io_global, ONLY : ionode
  USE wvfct,    ONLY : nbnd
  USE lsda_mod, ONLY : nspin

  implicit none

  INTEGER, EXTERNAL :: find_free_unit

  REAL(kind=DP) :: e_xc(nbnd,nspin)!exchange and correlation energies

  INTEGER :: iunu, iw


  if(ionode) then
     iunu = find_free_unit()

     open(unit=iunu,file=trim(tmp_dir)//trim(prefix)//'.dft_xc',status='unknown',form='unformatted')

     write(iunu) nbnd


     do iw=1,nbnd
        write(iunu) e_xc(iw,1)
        if(l_verbose) WRITE(*,*) 'SCRITTO e_XC 1', e_xc(iw,1)
     enddo
     if(nspin==2) then
        do iw=1,nbnd
           write(iunu) e_xc(iw,2)
           if(l_verbose) WRITE(*,*) 'SCRITTO e_XC 2', e_xc(iw,2)
        enddo
     endif
     close(iunu)
  endif

  return

END SUBROUTINE write_energies_xc
