!
! Copyright (C) 2001-2009 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!--------------------------------------------------------------------
SUBROUTINE stm (sample_bias, stmdos, istates)
  !--------------------------------------------------------------------
  !
  !     This routine calculates an stm image defined as the local density
  !     of states at the fermi energy.
  !     The bias of the sample is decided by sample_bias, states between
  !     ef and ef + sample_bias are taken into account.
  !     On output istates is the number of states used to compute the image.
  !     The slab must be oriented with the main axis along celldm(3).
  !     It may not properly work if the slab has two symmetric surfaces.
  !
  USE kinds, ONLY: DP
  USE constants, ONLY: tpi, rytoev
  USE io_global, ONLY : stdout
  USE cell_base, ONLY: omega, at
  USE fft_base,  ONLY: dfftp
  USE scatter_mod,  ONLY: gather_grid
  USE fft_interfaces, ONLY : fwfft, invfft
  USE gvect, ONLY: ngm, g, nl, nlm
  USE klist, ONLY: xk, lgauss, degauss, ngauss, wk, nks, nelec, ngk, igk_k
  USE ener, ONLY: ef
  USE symme, ONLY : sym_rho, sym_rho_init
  USE scf, ONLY: rho
  USE wvfct, ONLY: npwx, nbnd, wg, et
  USE control_flags, ONLY : gamma_only
  USE wavefunctions_module,  ONLY : evc, psic
  USE io_files, ONLY: iunwfc, nwordwfc
  USE constants,      ONLY : degspin
  USE mp,        ONLY : mp_max, mp_min, mp_sum
  USE mp_global, ONLY : inter_pool_comm
  !
  IMPLICIT NONE
  !
  REAL(DP), INTENT(IN) :: sample_bias
  REAL(DP), INTENT(OUT):: stmdos (dfftp%nr1x*dfftp%nr2x*dfftp%nr3x)
  ! the stm density of states
  INTEGER, INTENT(OUT):: istates
  ! the number of states to compute the image
  !
  !    And here the local variables
  !
  INTEGER :: ir, ig, ibnd, ik, nbnd_ocp, first_band, last_band, npw
  ! counters on 3D r points
  ! counter on g vectors
  ! counter on bands
  ! counter on k points
  ! number of occupied bands
  ! first band close enough to the specified energy range [down1:up1]
  ! last  band close enough to the specified energy range [down1:up1]

  real(DP) :: emin, emax, x, y, &
       w1, w2, up, up1, down, down1, t0, scnds
  COMPLEX(DP), PARAMETER :: i= (0.d0, 1.d0)

  real(DP), ALLOCATABLE :: gs (:,:)
  COMPLEX(DP), ALLOCATABLE :: psi (:,:)
  ! plane stm wfc

  real(DP), EXTERNAL :: w0gauss

  t0 = scnds ()
  ALLOCATE (gs( 2, npwx))
  ALLOCATE (psi(dfftp%nr1x, dfftp%nr2x))
  !
  stmdos(:) = 0.d0
  rho%of_r(:,:) = 0.d0
  WRITE( stdout, '(5x,"Use the true wfcs")')
  WRITE( stdout, '(5x,"Sample bias          =",f8.4, &
          &       " eV")') sample_bias * rytoev
  !
  IF (.not.lgauss) THEN
     !
     !  for semiconductors, add small broadening
     !
     nbnd_ocp = nint (nelec) / degspin

     IF (nbnd<=nbnd_ocp + 1) CALL errore ('stm', 'not enough bands', 1)
     emin = et (nbnd_ocp + 1, 1)
     DO ik = 2, nks
        emin = min (emin, et (nbnd_ocp + 1, ik) )
     ENDDO
#if defined(__MPI)
     ! find the minimum across pools
     CALL mp_min( emin, inter_pool_comm )
#endif
     emax = et (nbnd_ocp, 1)
     DO ik = 2, nks
        emax = max (emax, et (nbnd_ocp, ik) )
     ENDDO
#if defined(__MPI)
     ! find the maximum across pools
     CALL mp_max( emax, inter_pool_comm )
#endif
     ef = (emin + emax) * 0.5d0
     degauss = 0.00001d0

     ngauss = 0
     WRITE( stdout, '(/5x,"Occupied bands: ",i6)') nbnd_ocp
     WRITE( stdout, '(/5x,"  Fermi energy: ",f10.2," eV")') ef * rytoev
     WRITE( stdout, '(/5x,"    Gap energy: ",f10.2," eV")')  (emax - emin)  * rytoev
  ENDIF
  !
  !     take only the states in the energy window above or below the fermi
  !     energy as determined by the bias of the sample
  !
  IF (sample_bias>0) THEN
     up = ef + sample_bias
     down = ef
  ELSE
     up = ef
     down = ef + sample_bias
  ENDIF
  up1   = up   + 3.d0 * degauss
  down1 = down - 3.d0 * degauss

  DO ik = 1, nks
     DO ibnd = 1, nbnd
        IF (et (ibnd, ik) > down .and. et (ibnd, ik) < up) THEN
           wg (ibnd, ik) = wk (ik)
        ELSEIF (et (ibnd, ik) < down) THEN
           wg (ibnd, ik) = wk (ik) * w0gauss ( (down - et (ibnd, ik) ) &
                / degauss, ngauss)
        ELSEIF (et (ibnd, ik) > up) THEN
           wg (ibnd, ik) = wk (ik) * w0gauss ( (up - et (ibnd, ik) ) &
                / degauss, ngauss)
        ENDIF
     ENDDO
  ENDDO
  !
  istates = 0
  !
  !     here we sum for each k point the contribution
  !     of the wavefunctions to the stm dos
  !
  DO ik = 1, nks
     DO ibnd = 1, nbnd
        IF (et(ibnd,ik) < down1) first_band= ibnd+1
        IF (et(ibnd,ik) < up1)   last_band = ibnd
     ENDDO
     istates = istates +  (last_band - first_band + 1)

     npw = ngk(ik)
     CALL davcio (evc, 2*nwordwfc, iunwfc, ik, - 1)
     !
     IF (gamma_only) THEN
        !
        !     gamma only version of STM.
        !     Two bands computed in a single FT as in the main (PW) code
        !
        DO ibnd = first_band, last_band, 2
           w1 = wg (ibnd, ik) / omega
           !!! WRITE( stdout, * ) w1, ibnd, ik

           IF ( ibnd < last_band ) THEN
              w2 = wg (ibnd+1, ik) / omega
              !!! WRITE( stdout, * ) w2, ibnd+1, ik
           ELSE
              w2= 0.d0
           ENDIF
           !
           !     Compute the contribution of these states only if needed
           !
           psic(:) = (0.d0, 0.d0)
           IF ( ibnd < last_band ) THEN
              DO ig = 1, npw
                 psic(nl(igk_k(ig,ik)))  = &
                             evc(ig,ibnd) + (0.D0,1.D0) * evc(ig,ibnd+1)
                 psic(nlm(igk_k(ig,ik))) = &
                      conjg( evc(ig,ibnd) - (0.D0,1.D0) * evc(ig,ibnd+1) )
              ENDDO
           ELSE
              DO ig = 1, npw
                 psic(nl (igk_k(ig,ik))) =        evc(ig,ibnd)
                 psic(nlm(igk_k(ig,ik))) = conjg( evc(ig,ibnd) )
              ENDDO
           ENDIF

           CALL invfft ('Dense', psic, dfftp)
           DO ir = 1, dfftp%nnr
              rho%of_r (ir, 1) = rho%of_r (ir, 1) + w1* dble( psic(ir) )**2 + &
                                                    w2*aimag( psic(ir) )**2
           ENDDO
        ENDDO
     ELSE
        !
        !     k-point version of STM.
        !
        DO ibnd = first_band, last_band

           w1 = wg (ibnd, ik) / omega
           !!! WRITE( stdout, * ) w1, ibnd, ik
           !
           !     Compute the contribution of this state only if needed
           !
           psic(:) = (0.d0, 0.d0)
           DO ig = 1, npw
              psic(nl(igk_k(ig,ik)))  = evc(ig,ibnd)
           ENDDO

           CALL invfft ('Dense', psic, dfftp)
           DO ir = 1, dfftp%nnr
              rho%of_r (ir, 1) = rho%of_r (ir, 1) + w1 * &
                                ( dble(psic (ir) ) **2 + aimag(psic (ir) ) **2)
           ENDDO
        ENDDO
     ENDIF
  ENDDO
#if defined(__MPI)
  CALL mp_sum( rho%of_r, inter_pool_comm )
#endif
  !
  !     symmetrization of the stm dos
  !
  IF ( .not. gamma_only) THEN
     !
     CALL sym_rho_init (gamma_only)
     !
     psic(:) = cmplx ( rho%of_r(:,1), 0.0_dp, kind=dp)
     CALL fwfft ('Dense', psic, dfftp)
     rho%of_g(:,1) = psic(nl(:))
     CALL sym_rho (1, rho%of_g)
     psic(:) = (0.0_dp, 0.0_dp)
     psic(nl(:)) = rho%of_g(:,1)
     CALL invfft ('Dense', psic, dfftp)
     rho%of_r(:,1) = dble(psic(:))
  ENDIF
#if defined(__MPI)
  CALL gather_grid (dfftp, rho%of_r(:,1), stmdos)
#else
  stmdos(:) = rho%of_r(:,1)
#endif
  DEALLOCATE(psi)
  DEALLOCATE(gs)
  WRITE( stdout, '(/5x,"STM:",f10.2,"s cpu time")') scnds ()-t0
  !
#if defined(__MPI)
  CALL mp_sum( istates, inter_pool_comm )
#endif

  RETURN
END SUBROUTINE stm
