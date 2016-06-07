!
! Copyright (C) 2001-2013 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------------
SUBROUTINE wfc_gamma_real(itask,ispin)

!this subroutine writes the wfcs on real space - coarse grid
!on disk
!it works only at GAMMA
!it supports US pseudopotentials
!it also sets uo the array bec_gw


  USE kinds,                ONLY : DP
  USE gvect,                ONLY : gstart
  USE gvecs,              ONLY : nls, nlsm, doublegrid
   USE io_files,             ONLY : iunwfc, nwordwfc, diropn
  USE wvfct,                ONLY : nbnd, npwx, npw, wg, et
  USE mp,                   ONLY : mp_bcast
  USE io_global,            ONLY : stdout
  USE klist,                ONLY : lgauss, degauss, ngauss, nks, &
                                   nkstot, wk, xk, nelec, nelup, neldw, &
                                   two_fermi_energies, igk_k
  USE lsda_mod,             ONLY : lsda, nspin, current_spin, isk
  USE wavefunctions_module, ONLY : evc, psic
  USE io_files,             ONLY : diropn
  USE wannier_gw,           ONLY : becp_gw, becp_gw_c, l_verbose
  USE uspp
  USE control_flags,    ONLY : gamma_only
  USE fft_base,             ONLY : dfftp, dffts
  USE fft_interfaces,       ONLY : fwfft, invfft
 
  IMPLICIT NONE

  INTEGER, EXTERNAL :: find_free_unit

  INTEGER, INTENT(in) :: itask !if ==1 consider subspace{c'}
  INTEGER, INTENT(in) :: ispin!spin variable 1,2

  !
  INTEGER :: ikb, jkb, ijkb0, ih, jh, ijh, na, np
  INTEGER :: ir, is, ig, ibnd, ik
  INTEGER :: iunwfcreal
  LOGICAL :: exst
  REAL(kind=DP), ALLOCATABLE :: tmpreal(:)

  if(l_verbose) write(stdout,*) 'FUNCTION WFC_REAL' !ATTENZIONE
    FLUSH(stdout)

    allocate(tmpreal(dffts%nnr))

  IF(.not.gamma_only) THEN
       write(stdout,*) ' wfc_gamma_real only for GAMMA'
       stop
  ENDIF

  iunwfcreal=find_free_unit()
  CALL diropn( iunwfcreal, 'real_whole', dffts%nnr, exst )

       !
  
  !
  
  if  ( nkb > 0 .and. okvan) then
     CALL init_us_2( npw, igk_k(1,1), xk(1,1), vkb )
     if(itask/=1) then
        !CALL ccalbec( nkb, npwx, npw, nbnd, becp_gw, vkb, evc )
     else
        !CALL ccalbec( nkb, npwx, npw, nbnd, becp_gw_c, vkb, evc )
     endif
  endif
     !
  
          !
          ! ... here we compute the band energy: the sum of the eigenvalues
          !
 
  if(gstart==2) then
     do ibnd=1,nbnd
        evc(1,ibnd)=dble(evc(1,ibnd))
     enddo
  endif

  DO ibnd = 1, nbnd, 2
      if(l_verbose) write(stdout,*) 'IBND:',ibnd
     FLUSH(stdout)
     !
     psic(:) = ( 0.D0, 0.D0 )
        !
     IF ( ibnd < nbnd ) THEN
                !
                ! ... two ffts at the same time
                !
        psic(nls(igk_k(1:npw,1)))  = evc(1:npw,ibnd) + &
             ( 0.D0, 1.D0 ) * evc(1:npw,ibnd+1)
        psic(nlsm(igk_k(1:npw,1))) = CONJG( evc(1:npw,ibnd) - &
             ( 0.D0, 1.D0 ) * evc(1:npw,ibnd+1) )
                !
     ELSE
        !
        psic(nls(igk_k(1:npw,1)))  = evc(1:npw,ibnd)
        psic(nlsm(igk_k(1:npw,1))) = CONJG( evc(1:npw,ibnd) )
           !
     END IF
             !
      if(l_verbose) write(stdout,*) 'before'
     FLUSH(stdout)

     CALL invfft ('Wave', psic, dffts)
                    !
        !
        ! ... increment the charge density ...
        !
             !
     
     if(l_verbose) write(stdout,*) 'after'
     FLUSH(stdout)

     tmpreal(:)= DBLE(psic(:))
     CALL davcio( tmpreal,dffts%nnr,iunwfcreal,ibnd+(ispin-1)*nbnd,1)
     if(ibnd+1 <= nbnd) then
        tmpreal(:)=dimag(psic(:))
        CALL davcio( tmpreal,dffts%nnr,iunwfcreal,ibnd+1+(ispin-1)*nbnd,1)
     endif

             !
  END DO
          !
         
  close(iunwfcreal)
  deallocate(tmpreal)
       !
  END SUBROUTINE


  SUBROUTINE write_wfc_plot(itask)
!save wannier functions on disk for plotting
     USE io_files,             ONLY : nwordwfc, diropn
     USE wavefunctions_module, ONLY : evc
   
    implicit none

    INTEGER, EXTERNAL :: find_free_unit

    INTEGER, INTENT(in) :: itask!0 save MLWF, 1 save ULWF

    INTEGER :: iunplot
    LOGICAL :: exst

     iunplot=find_free_unit()
     if(itask==0) then
        CALL diropn( iunplot, 'wfc_mlwf', nwordwfc, exst )
     else
        CALL diropn( iunplot, 'wfc_ulwf', nwordwfc, exst )
     endif
     CALL davcio(evc,2*nwordwfc,iunplot,1,1)
     close(iunplot)
    return
  END SUBROUTINE write_wfc_plot
