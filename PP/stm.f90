!
! Copyright (C) 2001-2009 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!--------------------------------------------------------------------
subroutine stm (wf, sample_bias, z, dz, stmdos)
  !--------------------------------------------------------------------
  !
  !     This routine calculates an stm image defined as the local density
  !     of states at the fermi energy.
  !     The bias of the sample is decided by sample_bias, states between
  !     ef and ef + sample_bias are taken into account.
  !     It needs the workfunction wf. On output wf contains the number of
  !     states used to compute the image.
  !     The slab must be oriented with the main axis along celldm(3).
  !     It may not properly work if the slab has two symmetric surfaces.
  !
  USE kinds, ONLY: DP
  USE constants, ONLY: tpi, rytoev
  USE io_global, ONLY : stdout
  USE cell_base, ONLY: tpiba2, tpiba, omega, at, alat
  USE gvect, ONLY: nrx1, nrx2, nrx3, nr1, nr2, nr3, ngm, g, ecutwfc, &
       nl, nlm, nrxx
  USE klist, ONLY: xk, lgauss, degauss, ngauss, wk, nks, nelec
  USE ener, ONLY: ef
  USE symme, ONLY : sym_rho, sym_rho_init
  USE scf, ONLY: rho
  USE wvfct, ONLY: npwx, npw, nbnd, wg, et, g2kin, igk
  USE control_flags, ONLY : gamma_only
  USE wavefunctions_module,  ONLY : evc, psic
  USE io_files, ONLY: iunwfc, nwordwfc
  USE constants,      ONLY : degspin
  USE mp,        ONLY : mp_max, mp_min, mp_sum
  USE mp_global, ONLY : inter_pool_comm, intra_pool_comm 
  USE fft_base,  ONLY : grid_gather
!
  implicit none
  real(DP) :: sample_bias, z, dz, stmdos (nrx1 * nrx2 * nrx3)
  ! the stm density of states
  !
  !    And here the local variables
  !

  logical :: uguale

  integer :: istates, igs, npws, ir, irx, iry, irz, ig, ibnd, &
       ik, nbnd_ocp, first_band, last_band
  ! the number of states to compute the image
  ! counter on surface g vectors
  ! number of surface g vectors
  ! counters on 3D r points
  ! counter on g vectors
  ! counter on bands
  ! counter on k points
  ! number of occupied bands
  ! first band close enough to the specified energy range [down1:up1]
  ! last  band close enough to the specified energy range [down1:up1]

  real(DP) :: emin, emax, fac, wf, wf1, x, y, zz, &
       w1, w2, up, up1, down, down1, t0, scnds
  complex(DP), parameter :: i= (0.d0, 1.d0)

  real(DP), allocatable :: gs (:,:)
  complex(DP), allocatable :: a (:), psi (:,:)
  ! the coefficients of the matching wfc
  ! plane stm wfc

  real(DP), external :: w0gauss

  t0 = scnds ()
  allocate (gs( 2, npwx))    
  allocate (a ( npwx))    
  allocate (psi(nrx1, nrx2))    
  !
  stmdos(:) = 0.d0
  rho%of_r(:,:) = 0.d0
  WRITE( stdout, '(5x,"Use the true wfcs")')
  WRITE( stdout, '(5x,"Sample bias          =",f8.4, &
          &       " eV")') sample_bias * rytoev
  !
  if (.not.lgauss) then
     !
     !  for semiconductors, add small broadening
     !
     nbnd_ocp = nint (nelec) / degspin

     if (nbnd.le.nbnd_ocp + 1) call errore ('stm', 'not enough bands', 1)
     emin = et (nbnd_ocp + 1, 1)
     do ik = 2, nks
        emin = min (emin, et (nbnd_ocp + 1, ik) )
     enddo
#ifdef __PARA
     ! find the minimum across pools
     call mp_min( emin, inter_pool_comm )
#endif
     emax = et (nbnd_ocp, 1)
     do ik = 2, nks
        emax = max (emax, et (nbnd_ocp, ik) )
     enddo
#ifdef __PARA
     ! find the maximum across pools
     call mp_max( emax, inter_pool_comm )
#endif
     ef = (emin + emax) * 0.5d0
     degauss = 0.00001d0

     ngauss = 0
     WRITE( stdout, '(/5x,"Occupied bands: ",i6)') nbnd_ocp
     WRITE( stdout, '(/5x,"  Fermi energy: ",f10.2," eV")') ef * rytoev
     WRITE( stdout, '(/5x,"    Gap energy: ",f10.2," eV")')  (emax - emin)  * rytoev
  endif
  !
  !     take only the states in the energy window above or below the fermi
  !     energy as determined by the bias of the sample
  !
  if (sample_bias.gt.0) then
     up = ef + sample_bias
     down = ef
  else
     up = ef
     down = ef + sample_bias
  endif
  up1   = up   + 3.d0 * degauss
  down1 = down - 3.d0 * degauss

  do ik = 1, nks
     do ibnd = 1, nbnd
        if (et (ibnd, ik) > down .and. et (ibnd, ik) < up) then
           wg (ibnd, ik) = wk (ik)
        elseif (et (ibnd, ik) < down) then
           wg (ibnd, ik) = wk (ik) * w0gauss ( (down - et (ibnd, ik) ) &
                / degauss, ngauss)
        elseif (et (ibnd, ik) > up) then
           wg (ibnd, ik) = wk (ik) * w0gauss ( (up - et (ibnd, ik) ) &
                / degauss, ngauss)
        endif
     enddo
  enddo
  !
  istates = 0
  !
  !     here we sum for each k point the contribution
  !     of the wavefunctions to the stm dos
  !
  do ik = 1, nks
     DO ibnd = 1, nbnd
        if (et(ibnd,ik) < down1) first_band= ibnd+1
        if (et(ibnd,ik) < up1)   last_band = ibnd
     END DO
     istates = istates +  (last_band - first_band + 1)

     call gk_sort (xk (1, ik), ngm, g, ecutwfc / tpiba2, npw, igk, g2kin)
     call davcio (evc, nwordwfc, iunwfc, ik, - 1)
     !
     if (gamma_only) then
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
           END IF
           !
           !     Compute the contribution of these states only if needed
           !
           psic(:) = (0.d0, 0.d0)
           IF ( ibnd < last_band ) THEN
              do ig = 1, npw
                 psic(nl(igk(ig)))  = &
                             evc(ig,ibnd) + (0.D0,1.D0) * evc(ig,ibnd+1)
                 psic(nlm(igk(ig))) = &
                      CONJG( evc(ig,ibnd) - (0.D0,1.D0) * evc(ig,ibnd+1) )
              enddo
           ELSE
              do ig = 1, npw
                 psic(nl (igk(ig))) =        evc(ig,ibnd)
                 psic(nlm(igk(ig))) = CONJG( evc(ig,ibnd) )
              end do
           END IF

           call cft3 (psic, nr1, nr2, nr3, nrx1, nrx2, nrx3, 1)
           do ir = 1, nrxx
              rho%of_r (ir, 1) = rho%of_r (ir, 1) + w1* DBLE( psic(ir) )**2 + &
                                                    w2*AIMAG( psic(ir) )**2
           enddo
        END DO
     else
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
           do ig = 1, npw
              psic(nl(igk(ig)))  = evc(ig,ibnd)
           end do

           call cft3 (psic, nr1, nr2, nr3, nrx1, nrx2, nrx3, 1)
           do ir = 1, nrxx
              rho%of_r (ir, 1) = rho%of_r (ir, 1) + w1 * &
                                ( DBLE(psic (ir) ) **2 + AIMAG(psic (ir) ) **2)
           enddo
        END DO
     end if
  enddo
#ifdef __PARA
  call mp_sum( rho%of_r, inter_pool_comm )
#endif
  !
  !     symmetrization of the stm dos
  !
  IF ( .NOT. gamma_only) THEN
     !
     CALL sym_rho_init (gamma_only) 
     !
     psic(:) = CMPLX ( rho%of_r(:,1), 0.0_dp, KIND=dp)
     CALL cft3s (psic, nr1, nr2, nr3, nrx1, nrx2, nrx3, -1)
     rho%of_g(:,1) = psic(nl(:))
     CALL sym_rho (1, rho%of_g)
     psic(:) = (0.0_dp, 0.0_dp)
     psic(nl(:)) = rho%of_g(:,1)
     CALL cft3s (psic, nr1, nr2, nr3, nrx1, nrx2, nrx3, 1)
     rho%of_r(:,1) = DBLE(psic(:))
  END IF
#ifdef __PARA
  call grid_gather (rho%of_r(:,1), stmdos)
#else
  stmdos(:) = rho%of_r(:,1)
#endif
  deallocate(psi)
  deallocate(a)
  deallocate(gs)
  WRITE( stdout, '(/5x,"STM:",f10.2,"s cpu time")') scnds ()-t0
  !
  !     use wf to store istates
  !
  wf = istates
#ifdef __PARA
  call mp_sum( wf, inter_pool_comm )
#endif
  z = z / alat

  dz = dz / alat
  return
end subroutine stm
