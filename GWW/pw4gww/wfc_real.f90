! FOR GWW
!
! Copyright (C) 2001-2005 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! Author: P. Umari
! Modified by G. Stenuit
!
!----------------------------------------------------------------------------
SUBROUTINE wfc_gamma_real(itask)

!this subroutine writes the wfcs on real space - coarse grid
!on disk
!it works only at GAMMA
!it supports US pseudopotentials
!it also sets uo the array bec_gw


  USE kinds,                ONLY : DP
  USE control_flags,        ONLY : gamma_only
  USE fft_base,             ONLY : dffts
  USE fft_interfaces,       ONLY : invfft
  USE gvecs,              ONLY : nls, nlsm, doublegrid
  USE io_files,             ONLY : iunwfc, nwordwfc, iunigk
  USE io_files,             ONLY : find_free_unit, diropn
  USE wvfct,                ONLY : nbnd, npwx, npw, igk, wg, et
  USE mp,                   ONLY : mp_bcast
  USE io_global,            ONLY : stdout
  USE klist,                ONLY : lgauss, degauss, ngauss, nks, &
                                   nkstot, wk, xk, nelec, nelup, neldw, &
                                   two_fermi_energies, ngk
  USE lsda_mod,             ONLY : lsda, nspin, current_spin, isk
  USE wavefunctions_module, ONLY : evc, psic
  USE wannier_gw,           ONLY : becp_gw, becp_gw_c
  USE uspp
  USE becmod,               ONLY : calbec

  IMPLICIT NONE

  INTEGER, INTENT(in) :: itask !if ==1 consider subspace{c'}

  !
  INTEGER :: ikb, jkb, ijkb0, ih, jh, ijh, na, np
  INTEGER :: ir, is, ig, ibnd, ik
  INTEGER :: iunwfcreal
  LOGICAL :: exst
  REAL(kind=DP), ALLOCATABLE :: tmpreal(:)

  integer :: ipw, iigk

  write(stdout,*) 'FUNCTION WFC_REAL' !ATTENZIONE
    call flush_unit(stdout)

    allocate(tmpreal(dffts%nnr))
  !!! modified
  IF(.not.gamma_only) THEN
       write(stdout,*) ' wfc_gamma_real only for GAMMA'
       !!! this has been removed since in the pw4gww
       !!! file, we don't impose any more gamma point only
       stop
       !write(stdout,*) 'But the calculation will be done by taking only the WF at ik=1 (Gamma point)'
       !write(stdout,*) 'So, nks is set to 1'
       !call flush_unit(stdout)
       !nks=1
  ENDIF

  iunwfcreal=find_free_unit()
  CALL diropn( iunwfcreal, 'real_whole', dffts%nnr, exst )

       !
  IF ( nks > 1 ) REWIND( iunigk )
  !
  !!! for me, here, we have to remove the do loop
  !!! since we would like to keep only the wfc at gamma point
  k_loop: DO ik = 1, nks
     !
     IF ( lsda ) current_spin = isk(ik)
     !
     IF ( nks > 1 ) THEN
             !
        !!!!READ( iunigk ) npw, igk
        npw = ngk(ik)
!        READ( iunigk ) igk
!        CALL davcio(evc, nwordwfc,iunwfc, ik, -1 )
        !!!!CALL get_buffer ( evc, nwordwfc, iunwfc, ik)
        !
     END IF

     if  ( nkb > 0 .and. okvan) then
        CALL init_us_2( npw, igk, xk(1,ik), vkb )
        if(itask/=1) then
           call calbec(npw, vkb, evc, becp_gw, nbnd)
        else
           call calbec(npw, vkb, evc, becp_gw_c, nbnd)
        endif
     endif
          !

          !
          ! ... here we compute the band energy: the sum of the eigenvalues
          !
! --------------------------------------------------------
! pourquoi??? nei altri posti dove si fa fft non e cosi!!!
!     if(gstart==2) then
!       do ibnd=1,nbnd
!          evc(1,ibnd)=dble(evc(1,ibnd))
!       enddo
!     endif
! ----------------------------------------------------------
     DO ibnd = 1, nbnd, 2
        write(stdout,*) 'IBND:',ibnd
        call flush_unit(stdout)
             !
!  write(stdout,*) "nrxx=",nrxx
!  write(stdout,*) "nls=", ubound(nls(:))
   write(stdout,*) 'lbound and ubound of psic: ', lbound(psic), ubound(psic)
!  write(stdout,*)"npw, lbound evc", ubound(evc(:,:),1), ubound(evc(:,:),2)
!  write(stdout,*) "ubound(igk)", ubound(igk(:))
!  write(stdout,*) "nbnd", nbnd
!  write(stdout,*) "EVC 1"
!  write(stdout,*) "igk(1)=", igk(1)
!  write(stdout,*) nls(1)
!  write(stdout,*) "---------------"
        call flush_unit(stdout)
        psic(:) = ( 0.D0, 0.D0 )
        !
        IF ( ibnd < nbnd ) THEN
                !
                ! ... two ffts at the same time
                !
           psic(nls(igk(1:npw)))  = evc(1:npw,ibnd) + &
                ( 0.D0, 1.D0 ) * evc(1:npw,ibnd+1)
           psic(nlsm(igk(1:npw))) = CONJG( evc(1:npw,ibnd) - &
                ( 0.D0, 1.D0 ) * evc(1:npw,ibnd+1) )
                !
        ELSE
                !
           psic(nls(igk(1:npw)))  = evc(1:npw,ibnd)
           psic(nlsm(igk(1:npw))) = CONJG( evc(1:npw,ibnd) )
           !
        END IF
             !
        write(stdout,*) 'before'
        call flush_unit(stdout)

        CALL invfft ('Wave', psic, dffts)
             !
        !
        ! ... increment the charge density ...
        !
             !

        write(stdout,*) 'after'
        call flush_unit(stdout)

        tmpreal(:)= DBLE(psic(:))
        CALL davcio( tmpreal,dffts%nnr,iunwfcreal,ibnd,1)
        if(ibnd+1 <= nbnd) then
           tmpreal(:)=dimag(psic(:))
           CALL davcio( tmpreal,dffts%nnr,iunwfcreal,ibnd+1,1)
        endif

             !
     END DO ! over the bands
          !

          !
  END DO k_loop
  CLOSE(iunwfcreal)
  deallocate(tmpreal)
       !
  END SUBROUTINE
  !
  !
SUBROUTINE wfc_gamma_real_after_rot(itask)
  !
!this subroutine writes the wfcs on real space - coarse grid
!on disk
!it works only at GAMMA
!it supports US pseudopotentials
!it also sets uo the array bec_gw
! HAS TO BE USED AFTER rotate_wannier_gamma in order to read the evc1 from files wfc_w
  USE kinds,                ONLY : DP
  USE control_flags,        ONLY : gamma_only
  USE fft_base,             ONLY : dffts
  USE fft_interfaces,       ONLY : invfft
  USE gvecs,              ONLY : nls, nlsm,  doublegrid
  USE io_files,             ONLY : iunwfc, nwordwfc, iunigk
  USE io_files,             ONLY : find_free_unit, diropn
  USE wvfct,                ONLY : nbnd, npwx, npw, igk, wg, et
  USE mp,                   ONLY : mp_bcast
  USE io_global,            ONLY : stdout
  USE klist,                ONLY : lgauss, degauss, ngauss, nks, &
                                   nkstot, wk, xk, nelec, nelup, neldw, &
                                   two_fermi_energies, ngk
  USE lsda_mod,             ONLY : lsda, nspin, current_spin, isk
  USE wavefunctions_module, ONLY : psic
  USE wannier_gw,           ONLY : becp_gw, becp_gw_c, nbnd_normal
  USE uspp
  USE becmod,               ONLY : calbec
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(in) :: itask !if ==1 consider subspace{c'}
  !
  INTEGER :: ikb, jkb, ijkb0, ih, jh, ijh, na, np
  INTEGER :: ir, is, ig, ibnd, ik
  INTEGER :: iunwfcreal, iun_wannier
  LOGICAL :: exst
  REAL(kind=DP), ALLOCATABLE :: tmpreal(:)
  integer :: ipw, iigk
  !
  ! declaration of temporary wannierized evc --> evc1
  COMPLEX(kind=DP), ALLOCATABLE :: evc1(:,:)
  !
  write(stdout,*) 'FUNCTION WFC_REAL' !ATTENZIONE
  call flush_unit(stdout)
  !
  allocate(tmpreal(dffts%nnr))
  !
  IF(.not.gamma_only) THEN
       write(stdout,*) ' wfc_gamma_real only for GAMMA'
       stop
  ENDIF
  !
  iunwfcreal=find_free_unit()
  CALL diropn( iunwfcreal, 'real_whole', dffts%nnr, exst )
  !
  allocate(evc1(npw,nbnd_normal))
  iun_wannier = find_free_unit()
  call diropn(iun_wannier,"wfc_w",2*nwordwfc,exst)
  call davcio(evc1,nwordwfc,iun_wannier,1,-1)
  close(iun_wannier)
  !
  IF ( nks > 1 ) REWIND( iunigk )
  !
  !!! for me, here, we have to remove the do loop
  !!! since we would like to keep only the wfc at gamma point
  k_loop: DO ik = 1, nks
     !
     IF ( lsda ) current_spin = isk(ik)
     !
     IF ( nks > 1 ) THEN
        !
        npw = ngk(ik)
        !
     END IF
     !
     if  ( nkb > 0 .and. okvan) then
        CALL init_us_2( npw, igk, xk(1,ik), vkb )
        if(itask/=1) then
           call calbec(npw, vkb, evc1, becp_gw, nbnd)
        else
           call calbec(npw, vkb, evc1, becp_gw_c, nbnd)
        endif
     endif
     !
     DO ibnd = 1, nbnd, 2
        write(stdout,*) 'IBND:',ibnd
        call flush_unit(stdout)
        !
        write(stdout,*) 'lbound and ubound of psic: ', lbound(psic), ubound(psic)
        call flush_unit(stdout)
        psic(:) = ( 0.D0, 0.D0 )
        !
        IF ( ibnd < nbnd ) THEN
                !
                ! ... two ffts at the same time
                !
           psic(nls(igk(1:npw)))  = evc1(1:npw,ibnd) + &
                ( 0.D0, 1.D0 ) * evc1(1:npw,ibnd+1)
           psic(nlsm(igk(1:npw))) = CONJG( evc1(1:npw,ibnd) - &
                ( 0.D0, 1.D0 ) * evc1(1:npw,ibnd+1) )
                !
        ELSE
                !
           psic(nls(igk(1:npw)))  = evc1(1:npw,ibnd)
           psic(nlsm(igk(1:npw))) = CONJG( evc1(1:npw,ibnd) )
           !
        END IF
        !
        CALL invfft ('Wave', psic, dffts)
        write(stdout,*) 'after'
        call flush_unit(stdout)

        tmpreal(:)= DBLE(psic(:))
        CALL davcio( tmpreal,dffts%nnr,iunwfcreal,ibnd,1)
        if(ibnd+1 <= nbnd) then
           tmpreal(:)=dimag(psic(:))
           CALL davcio( tmpreal,dffts%nnr,iunwfcreal,ibnd+1,1)
        endif
        !
     END DO ! over the bands
     !
  END DO k_loop
  CLOSE(iunwfcreal)
  deallocate(tmpreal)
  !
END SUBROUTINE wfc_gamma_real_after_rot
  !
  !
SUBROUTINE write_wfc_plot(itask)
!save wannier functions on disk for plotting
     USE io_files,             ONLY : nwordwfc
     USE wavefunctions_module, ONLY : evc
     USE io_files,             ONLY : find_free_unit, diropn

    implicit none
    INTEGER, INTENT(in) :: itask!0 save MLWF, 1 save ULWF

    INTEGER :: iunplot
    LOGICAL :: exst

     iunplot=find_free_unit()
     if(itask==0) then
        CALL diropn( iunplot, 'wfc_mlwf', nwordwfc, exst )
     else
        CALL diropn( iunplot, 'wfc_ulwf', nwordwfc, exst )
     endif
     CALL davcio(evc,nwordwfc,iunplot,1,1)
     close(iunplot)
    return
  END SUBROUTINE write_wfc_plot
