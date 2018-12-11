  !
  ! Copyright (C) 2010-2016 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino  
  ! Copyright (C) 2007-2009 Roxana Margine
  ! 
  ! This file is distributed under the terms of the GNU General Public         
  ! License. See the file `LICENSE' in the root directory of the               
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .
  !
  !----------------------------------------------------------------------
  MODULE superconductivity_aniso
  !----------------------------------------------------------------------
  !! 
  !! This module contains all the subroutines linked with superconductivity using  
  !! the isotropic or anisotropic Eliashberg formalism. 
  !! 
  IMPLICIT NONE
  ! 
  CONTAINS
    !                                                                            
    !-----------------------------------------------------------------------
    SUBROUTINE eliashberg_aniso_iaxis
    !-----------------------------------------------------------------------
    !
    ! This routine is the driver of the self-consistent cycle for the anisotropic 
    ! Eliashberg equations on the imaginary-axis.  
    !
    USE kinds,         ONLY : DP
    USE io_global,     ONLY : stdout
    USE control_flags, ONLY : iverbosity
    USE epwcom,        ONLY : nsiter, nstemp, broyden_beta, broyden_ndim, & 
                              limag, lpade, lacon, fsthick, imag_read, wscut
    USE eliashbergcom, ONLY : nsw, nsiw, ADelta, ADeltap, ADeltai, ADeltaip, &
                              estemp, nkfs, nbndfs, ekfs, ef0
    USE superconductivity, ONLY : free_energy, dos_quasiparticle, gen_freqgrid_iaxis, &
                                  mem_size_eliashberg, deallocate_eliashberg_aniso_iaxis, & 
                                  deallocate_eliashberg_aniso_raxis, deallocate_eliashberg
    USE constants_epw, ONLY : kelvin2eV, ci, pi
    USE io_global,     ONLY : ionode_id
    USE mp_global,     ONLY : inter_pool_comm
    USE mp,            ONLY : mp_bcast, mp_barrier
    USE mp_world,      ONLY : mpime
    USE io_eliashberg, ONLY : eliashberg_read_aniso_iaxis
    ! 
    IMPLICIT NONE
    !
    ! Local variables
    INTEGER :: itemp
    !! Counter on temperature index
    INTEGER :: iter
    !! Counter on iteration steps
    INTEGER :: N
    !! Maximum nr. frequency points in Pade approx
    INTEGER :: ik
    !! Counter on k-points
    INTEGER :: ibnd
    !! Counter on bands
    INTEGER :: imelt
    !! Counter memory
    REAL(DP) :: rdeltaout(nsw), rdeltain(nsw), cdeltaout(nsw), cdeltain(nsw)
    !! Temporary variables for mix_broyden in analytic continuation
    REAL(DP) :: tcpu
    REAL(DP), EXTERNAL :: get_clock
    LOGICAL :: conv
    !! True if calculation is converged
    !
    CALL start_clock( 'aniso_iaxis' )
    !
    DO itemp = 1, nstemp ! loop over temperature
       !
       WRITE(stdout,'(a)') '  '
       WRITE(stdout,'(5x,a,i3,a,f8.4,a,a,i3,a)') 'temp(', itemp, ') = ', estemp(itemp)/kelvin2eV, ' K '
       WRITE(stdout,'(a)') '  '
       IF ( limag .AND. .not. imag_read ) THEN 
          WRITE(stdout,'(5x,a)') 'Solve anisotropic Eliashberg equations on imaginary-axis ' 
       ELSEIF ( limag .AND. imag_read ) THEN
          WRITE(stdout,'(5x,a)') 'Read from file Delta and Znorm on imaginary-axis '
       ENDIF
       WRITE(stdout,'(a)') '  '
       WRITE(stdout,'(5x,a,i6,a,i6)') 'Total number of frequency points nsiw ( ', itemp, ' ) = ', nsiw(itemp)
       WRITE(stdout,'(5x,a,f10.4)') 'Cutoff frequency wscut = ', (2.d0*nsiw(itemp)+1)*pi*estemp(itemp)
       WRITE(stdout,'(a)') '  '
       CALL start_clock( 'iaxis_imag' )
       CALL gen_freqgrid_iaxis( itemp )
       !
       IF ( ( limag .AND. .not. imag_read ) .OR. ( limag .AND. imag_read .AND. itemp .ne. 1 ) ) THEN
          iter = 1
          conv = .false.
          DO WHILE ( .not. conv .AND. iter .le. nsiter ) 
             CALL sum_eliashberg_aniso_iaxis( itemp, iter, conv )
             IF (mpime .eq. ionode_id) THEN
                DO ik = 1, nkfs
                   DO ibnd = 1, nbndfs
                      IF ( abs( ekfs(ibnd,ik) - ef0 ) .lt. fsthick ) THEN
                         CALL mix_broyden_aniso( ik, ibnd, nsiw(itemp), & 
                              ADeltai(ibnd,ik,:), ADeltaip(ibnd,ik,:), broyden_beta, iter, broyden_ndim, conv )
                      ENDIF
                   ENDDO
                ENDDO
             ENDIF
             CALL mp_bcast( ADeltai, ionode_id, inter_pool_comm )
             CALL mp_bcast( ADeltaip, ionode_id, inter_pool_comm )
             CALL mp_barrier(inter_pool_comm)
             iter = iter + 1
          ENDDO ! iter
          !
          IF ( conv ) THEN
             IF ( ALLOCATED(ADeltaip) ) DEALLOCATE(ADeltaip)
             !
             ! SP : Only print the Free energy if the user want it
             !
             IF ( iverbosity .eq. 2 ) THEN
                IF (mpime .eq. ionode_id) THEN
                   CALL free_energy( itemp )
                ENDIF
                CALL mp_barrier(inter_pool_comm)
             ENDIF
             !
             WRITE(stdout,'(a)') '  '
             CALL stop_clock( 'iaxis_imag' )
             CALL print_clock( 'iaxis_imag' )
          ELSEIF ( .not. conv .AND. (iter-1) .eq. nsiter ) THEN
             CALL deallocate_eliashberg
             WRITE(stdout,'(a)') 'not converged  '
             CALL stop_clock( 'iaxis_imag' )
             CALL print_clock( 'iaxis_imag' )
             CALL errore('sum_eliashberg_aniso_iaxis','convergence was not reached',1)
             RETURN
          ENDIF
       ELSEIF ( limag .AND. imag_read .AND. itemp .eq. 1 ) THEN
          CALL eliashberg_read_aniso_iaxis( itemp )
       ENDIF
       !
       IF ( lpade ) THEN 
          WRITE(stdout,'(a)') '  '  
          WRITE(stdout,'(5x,a)') 'Pade approximant of anisotropic Eliashberg equations from imaginary-axis to real-axis'
          WRITE(stdout,'(5x,a,f10.4)') 'Cutoff frequency wscut = ', wscut
          WRITE(stdout,'(a)') '  '
          CALL start_clock( 'raxis_pade' )
          conv = .false.
          N = 90 * nsiw(itemp) / 100
          IF ( mod(N,2) .ne. 0 ) N = N + 1
          CALL pade_cont_aniso_iaxis_to_raxis( itemp, N, conv )
          !
          IF ( conv ) THEN
             IF (mpime .eq. ionode_id) THEN
                CALL dos_quasiparticle( itemp )
             ENDIF
             CALL mp_barrier(inter_pool_comm)
             CALL stop_clock( 'raxis_pade' )
             CALL print_clock( 'raxis_pade' )
             WRITE(stdout,'(a)') '  '
          ELSEIF ( .not. conv  .AND. (iter-1) .eq. nsiter ) THEN
             CALL deallocate_eliashberg
             WRITE(stdout,'(a)') '  '
             CALL stop_clock( 'raxis_pade' )
             CALL print_clock( 'raxis_pade' )
             CALL errore('pade_cont_iso_iaxis_to_raxis','converged was not reached',1)
             RETURN
          ENDIF
       ENDIF
       !
       IF ( lacon ) THEN 
          WRITE(stdout,'(a)') '  '
          WRITE(stdout,'(5x,a)') 'Analytic continuation of anisotropic Eliashberg equations from imaginary-axis to real-axis'
          WRITE(stdout,'(a)') '  '
          WRITE(stdout,'(5x,a,i6)') 'Total number of frequency points nsw = ', nsw
          WRITE(stdout,'(5x,a,f10.4)') 'Cutoff frequency wscut = ', wscut
          WRITE(stdout,'(a)') '    '
          CALL start_clock( 'raxis_acon' )
          !
          iter = 1
          conv = .false.
          DO WHILE ( .not. conv .AND. iter .le. nsiter )
            CALL analytic_cont_aniso_iaxis_to_raxis( itemp, iter, conv )
            IF (mpime .eq. ionode_id) THEN
              DO ik = 1, nkfs
                DO ibnd = 1, nbndfs
                  IF ( abs( ekfs(ibnd,ik) - ef0 ) .lt. fsthick ) THEN
                    rdeltain(:)  = real(ADeltap(ibnd,ik,:))
                    cdeltain(:)  = aimag(ADeltap(ibnd,ik,:))
                    rdeltaout(:) = real(ADelta(ibnd,ik,:))
                    cdeltaout(:) = aimag(ADelta(ibnd,ik,:))
                    CALL mix_broyden_aniso ( ik, ibnd, nsw, rdeltaout, rdeltain, broyden_beta, iter, broyden_ndim, conv )
                    CALL mix_broyden2_aniso( ik, ibnd, nsw, cdeltaout, cdeltain, broyden_beta, iter, broyden_ndim, conv )
                    ADeltap(ik,ibnd,:) = rdeltain(:) + ci * cdeltain(:)
                  ENDIF
                ENDDO
              ENDDO
            ENDIF
            CALL mp_bcast( ADelta, ionode_id, inter_pool_comm )
            CALL mp_bcast( ADeltap, ionode_id, inter_pool_comm )
            CALL mp_barrier(inter_pool_comm)
            iter = iter + 1
          ENDDO ! iter
          !
          IF ( conv ) THEN 
            IF (mpime .eq. ionode_id) THEN
               CALL dos_quasiparticle( itemp )
            ENDIF
            CALL mp_barrier(inter_pool_comm)
            WRITE(stdout,'(a)') '  '
            CALL stop_clock( 'raxis_acon' )
            CALL print_clock( 'raxis_acon' )
          ELSEIF ( .not. conv .AND. (iter-1) .eq. nsiter ) THEN
            CALL deallocate_eliashberg
            WRITE(stdout,'(a)') '  '
            CALL stop_clock( 'raxis_acon' )
            CALL print_clock( 'raxis_acon' )
            CALL errore('analytic_cont_aniso_iaxis_to_raxis','convergence was not reached',1)
            RETURN
          ENDIF
          !
       ENDIF
       !
       CALL deallocate_eliashberg_aniso_iaxis
       ! remove memory allocated for wsi, Deltai, Znormi, NZnormi, ADeltai, AZnormi, NAZnormi
       imelt = ( 4 + 3 * nbndfs * nkfs ) * nsiw(itemp) 
       CALL mem_size_eliashberg( -imelt )
       !
       CALL deallocate_eliashberg_aniso_raxis
       IF ( lpade ) THEN
         ! remove memory allocated for ws, Delta, Znorm, ADelta, AZnorm
         imelt = nsw + 2 * ( 2 + 2 * nbndfs * nkfs ) * nsw
         CALL mem_size_eliashberg( -imelt )
       ELSEIF ( lacon ) THEN 
         ! remove memory allocated for ws, Delta, Znorm, ADelta, ADeltap, AZnorm, AZnormp
         imelt = nsw + 2 * ( 2 + 4 * nbndfs * nkfs ) * nsw
         CALL mem_size_eliashberg( -imelt )
       ENDIF
       ! 
       tcpu = get_clock('aniso_iaxis')
       WRITE(stdout,'(5x,a,i3,a,f18.2,a)') 'itemp = ', itemp, '   total cpu time :', tcpu, ' secs'
       !
    ENDDO ! itemp
    !
    CALL stop_clock( 'aniso_iaxis' )
    !
    RETURN
    !
    END SUBROUTINE eliashberg_aniso_iaxis
    !
    !-----------------------------------------------------------------------
    SUBROUTINE sum_eliashberg_aniso_iaxis( itemp, iter, conv ) 
    !-----------------------------------------------------------------------
    !
    ! This routine solves the anisotropic Eliashberg equations on the imaginary-axis
    !
    USE kinds,         ONLY : DP
    USE io_global,     ONLY : stdout
    USE elph2,         ONLY : wqf
    USE epwcom,        ONLY : nsiter, nstemp, muc, conv_thr_iaxis, fsthick
    USE eliashbergcom, ONLY : nsiw, estemp, gap0, gap, Agap, wsi, AKeri, limag_fly, & 
                              NAZnormi, AZnormi, ADeltai, ADeltaip, NZnormi, Znormi, & 
                              Deltai, wsphmax, nkfs, nbndfs, dosef, ef0, ixkqf, ixqfs, & 
                              nqfs, wkfs, w0g, ekfs
    USE superconductivity, ONLY : mem_size_eliashberg, eliashberg_memlt_aniso_iaxis
    USE constants_epw, ONLY : pi, zero, czero 
    USE io_global,     ONLY : ionode_id
    USE mp_global,     ONLY : inter_pool_comm
    USE mp_world,      ONLY : mpime
    USE mp,            ONLY : mp_bcast, mp_barrier, mp_sum
    USE io_eliashberg, ONLY : eliashberg_write_iaxis
    USE division,      ONLY : fkbounds
    ! 
    IMPLICIT NONE
    !
    INTEGER, INTENT (in) :: itemp
    !! Counter on temperature index
    INTEGER, INTENT(in) :: iter
    !! Counter on iteration steps
    LOGICAL, INTENT(inout) :: conv
    !! True if the calculation is converged
    !
    ! Local variables
    INTEGER :: iw, iwp
    !! Counter on frequency imag-axis
    INTEGER :: ik
    !! Counter on k-points
    INTEGER :: iq
    !! Counter on q-points for which k+sign*q is within the Fermi shell
    INTEGER :: iq0
    !! Index of iq on full q-mesh
    INTEGER :: lower_bnd
    !! Lower bounds index after k paral
    INTEGER :: upper_bnd
    !! Upper bounds index after k paral
    INTEGER :: ibnd
    !! Counter on bands
    INTEGER :: jbnd
    !! Counter on bands
    INTEGER :: imelt
    !! Counter memory
    REAL(DP) :: kernelp, kernelm, lambdap, lambdam
    !! lambdam - K_{-}(n,n'T))
    !! lambdap - K_{+}(n,n',T)
    !! kernelm = lambdam - lambdap
    !! kernelp = lambdam + lambdap
    REAL(DP) :: absdelta, reldelta, errdelta
    !! Errors in Delta
    REAL(DP) :: esqrt, weight
    REAL(DP), ALLOCATABLE :: wesqrt(:,:,:), desqrt(:,:,:)
    !! Temporary variables
    REAL(DP), ALLOCATABLE, SAVE :: Deltaold(:)
    !! gap
    !
    IF ( .not. ALLOCATED(wesqrt) ) ALLOCATE( wesqrt(nbndfs,nkfs,nsiw(itemp)) )
    IF ( .not. ALLOCATED(desqrt) ) ALLOCATE( desqrt(nbndfs,nkfs,nsiw(itemp)) )
    wesqrt(:,:,:) = zero
    desqrt(:,:,:) = zero
    !
    IF ( iter .eq. 1 ) THEN
       !
       IF ( itemp .eq. 1 ) THEN 
          ! get the size of required memory for  gap, Agap
          imelt = ( 1 + nbndfs * nkfs ) * nstemp 
          CALL mem_size_eliashberg( imelt )
       ENDIF
       !
       ! get the size of required memory for  
       ! wesqrt, desqrt, Deltai, Znormi, NZnormi, ADeltai, ADeltaip, AZnormi, NAZnormi, Deltaold
       imelt = ( 4 + 6 * nbndfs * nkfs ) * nsiw(itemp)
       CALL mem_size_eliashberg( imelt )
       !
       IF ( .not. ALLOCATED(gap) )       ALLOCATE( gap(nstemp) )
       IF ( .not. ALLOCATED(Agap) )      ALLOCATE( Agap(nbndfs,nkfs,nstemp) )
       IF ( .not. ALLOCATED(Deltai) )    ALLOCATE( Deltai(nsiw(itemp)) )
       IF ( .not. ALLOCATED(Znormi) )    ALLOCATE( Znormi(nsiw(itemp)) )
       IF ( .not. ALLOCATED(NZnormi) )   ALLOCATE( NZnormi(nsiw(itemp)) )
       IF ( .not. ALLOCATED(ADeltai) )   ALLOCATE( ADeltai(nbndfs,nkfs,nsiw(itemp)) )
       IF ( .not. ALLOCATED(ADeltaip) )  ALLOCATE( ADeltaip(nbndfs,nkfs,nsiw(itemp)) )
       IF ( .not. ALLOCATED(AZnormi) )   ALLOCATE( AZnormi(nbndfs,nkfs,nsiw(itemp)) )
       IF ( .not. ALLOCATED(NAZnormi) )  ALLOCATE( NAZnormi(nbndfs,nkfs,nsiw(itemp)) )
       gap(itemp) = zero
       Agap(:,:,itemp) = zero
       ADeltaip(:,:,:) = zero
       !
       DO ik = 1, nkfs
          DO ibnd = 1, nbndfs
             IF ( abs( ekfs(ibnd,ik) - ef0 ) .lt. fsthick ) THEN
                DO iw = 1, nsiw(itemp)
                   IF ( wsi(iw) .lt. 2.d0*wsphmax ) THEN
                      ADeltaip(ibnd,ik,iw) = gap0
                   ELSE
                      ADeltaip(ibnd,ik,iw) = zero
                   ENDIF
                ENDDO
             ENDIF
          ENDDO ! ibnd
       ENDDO ! ik
       !
       CALL eliashberg_memlt_aniso_iaxis( itemp )
       IF ( .not. limag_fly ) CALL kernel_aniso_iaxis( itemp )
       !
    ENDIF 
    Deltai(:) = zero
    Znormi(:) = zero
    NZnormi(:) = zero
    ADeltai(:,:,:) = zero
    AZnormi(:,:,:) = zero
    NAZnormi(:,:,:) = zero
    !
    CALL fkbounds( nkfs, lower_bnd, upper_bnd )
    !
    DO ik = lower_bnd, upper_bnd
       DO ibnd = 1, nbndfs
          IF ( ( abs( ekfs(ibnd,ik) - ef0 ) .lt. fsthick ) ) THEN
             DO iq = 1, nqfs(ik)
                ! iq0 - index of q-point on the full q-mesh
                iq0 = ixqfs(ik,iq)
                DO jbnd = 1, nbndfs
                   IF ( ( abs( ekfs(jbnd,ixkqf(ik,iq0)) - ef0 ) .lt. fsthick ) ) THEN
                      weight = wqf(iq) * w0g(jbnd,ixkqf(ik,iq0)) / dosef
                      DO iw = 1, nsiw(itemp) ! loop over omega
                         DO iwp = 1, nsiw(itemp) ! loop over omega_prime
                            !
                            ! this step is performed at each iter step only for iw=1 
                            IF ( iw .eq. 1 ) THEN
                               esqrt = 1.d0 / sqrt( wsi(iwp)**2.d0 + ADeltaip(jbnd,ixkqf(ik,iq0),iwp)**2.d0 )
                               wesqrt(jbnd,ixkqf(ik,iq0),iwp) = wsi(iwp) * esqrt 
                               desqrt(jbnd,ixkqf(ik,iq0),iwp) = ADeltaip(jbnd,ixkqf(ik,iq0),iwp) * esqrt 
                            ENDIF
                            IF ( limag_fly ) THEN 
                               CALL lambdar_aniso_ver1( ik, iq, ibnd, jbnd, wsi(iw) - wsi(iwp), lambdam )
                               CALL lambdar_aniso_ver1( ik, iq, ibnd, jbnd, wsi(iw) + wsi(iwp), lambdap )
                            ELSE
                               lambdam = AKeri( ik, iq, ibnd, jbnd, abs(iw-iwp)+1 )
                               lambdap = AKeri( ik, iq, ibnd, jbnd, abs(iw+iwp) )
                            ENDIF
                            kernelm = lambdam - lambdap
                            kernelp = lambdam + lambdap
                            NAZnormi(ibnd,ik,iw) = NAZnormi(ibnd,ik,iw) + weight * kernelm
                            AZnormi(ibnd,ik,iw) = AZnormi(ibnd,ik,iw) + weight * wesqrt(jbnd,ixkqf(ik,iq0),iwp) &
                                                * kernelm
                            ADeltai(ibnd,ik,iw) = ADeltai(ibnd,ik,iw) + weight * desqrt(jbnd,ixkqf(ik,iq0),iwp) &
                                                * ( kernelp - 2.d0 * muc )
                         ENDDO ! iwp
                      ENDDO ! iw
                   ENDIF
                ENDDO ! jbnd
             ENDDO ! iq
          ENDIF
       ENDDO ! ibnd
    ENDDO ! ik
    !
    IF( ALLOCATED(wesqrt) ) DEALLOCATE(wesqrt)
    IF( ALLOCATED(desqrt) ) DEALLOCATE(desqrt)
    !
    ! collect contributions from all pools 
    CALL mp_sum( AZnormi, inter_pool_comm )
    CALL mp_sum( NAZnormi, inter_pool_comm )
    CALL mp_sum( ADeltai, inter_pool_comm )
    CALL mp_barrier(inter_pool_comm)
    !
    IF (mpime .eq. ionode_id) THEN
      IF ( iter .eq. 1 ) THEN
         IF ( .not. ALLOCATED(Deltaold) ) ALLOCATE( Deltaold(nsiw(itemp)) )
         Deltaold(:) = gap0
      ENDIF
      absdelta = zero
      reldelta = zero
      DO iw = 1, nsiw(itemp) ! loop over omega
         DO ik = 1, nkfs
            DO ibnd = 1, nbndfs
               IF ( abs( ekfs(ibnd,ik) - ef0 ) .lt. fsthick ) THEN
                  weight = 0.5d0 * wkfs(ik) * w0g(ibnd,ik) / dosef
                  Znormi(iw) = Znormi(iw) + weight * AZnormi(ibnd,ik,iw)
                  Deltai(iw) = Deltai(iw) + weight * ADeltai(ibnd,ik,iw)
                  NAZnormi(ibnd,ik,iw) = 1.d0 + pi * estemp(itemp) * NAZnormi(ibnd,ik,iw) / wsi(iw)
                  AZnormi(ibnd,ik,iw) = 1.d0 + pi * estemp(itemp) * AZnormi(ibnd,ik,iw) / wsi(iw)
                  ADeltai(ibnd,ik,iw) = pi * estemp(itemp) * ADeltai(ibnd,ik,iw) / AZnormi(ibnd,ik,iw)
               ENDIF
            ENDDO ! ibnd
         ENDDO ! ik
         NZnormi(iw) = 1.d0 + pi * estemp(itemp) * NZnormi(iw) / wsi(iw)
         Znormi(iw) = 1.d0 + pi * estemp(itemp) * Znormi(iw) / wsi(iw)
         Deltai(iw) = pi * estemp(itemp) * Deltai(iw) / Znormi(iw)
         reldelta = reldelta + abs( Deltai(iw) - Deltaold(iw) )
         absdelta = absdelta + abs( Deltai(iw) )
      ENDDO ! iw
      errdelta = reldelta / absdelta
      Deltaold(:) = Deltai(:)
      !
      WRITE(stdout,'(5x,a,i6,a,ES20.10,a,ES20.10,a,ES20.10,a,ES20.10)') 'iter = ', iter, & 
                   '   relerr = ', errdelta, '   abserr = ', reldelta / dble(nsiw(itemp)), &
                   '   Znormi(1) = ', Znormi(1), '   Deltai(1) = ', Deltai(1)
      !
      IF ( errdelta .lt. conv_thr_iaxis) conv = .true.
      IF ( errdelta .lt. conv_thr_iaxis .OR. iter .eq. nsiter ) THEN
         gap(itemp) = Deltai(1)
         gap0 = gap(itemp)
         !
         CALL eliashberg_write_iaxis( itemp )
         !
      ENDIF
      !
      IF ( conv .OR. iter .eq. nsiter ) THEN
         IF( ALLOCATED(Deltaold) ) DEALLOCATE(Deltaold)
         WRITE(stdout,'(5x,a,i6)') 'Convergence was reached in nsiter = ', iter
      ENDIF
      IF ( .not. conv .AND. iter .eq. nsiter ) THEN
         WRITE(stdout,'(a)') ' '
         WRITE(stdout,'(5x,a,i6)') 'Convergence was not reached in nsiter = ', iter
         CALL errore('sum_eliashberg_aniso_iaxis','increase nsiter or reduce conv_thr_iaxis',1)
      ENDIF
    ENDIF
    CALL mp_bcast( Deltai, ionode_id, inter_pool_comm )
    CALL mp_bcast( Znormi, ionode_id, inter_pool_comm )
    CALL mp_bcast( NZnormi, ionode_id, inter_pool_comm )
    CALL mp_bcast( AZnormi, ionode_id, inter_pool_comm )
    CALL mp_bcast( NAZnormi, ionode_id, inter_pool_comm )
    CALL mp_bcast( gap0, ionode_id, inter_pool_comm )
    CALL mp_bcast( gap, ionode_id, inter_pool_comm )
    CALL mp_bcast( Agap, ionode_id, inter_pool_comm )
    CALL mp_bcast( conv, ionode_id, inter_pool_comm )
    CALL mp_barrier(inter_pool_comm)
    !
    IF ( conv .OR. iter .eq. nsiter ) THEN 
       !
       ! remove memory allocated for wesqrt, desqrt, ADeltaip, Deltaold
       imelt = ( 1 + 3 * nbndfs * nkfs ) * nsiw(itemp)
       CALL mem_size_eliashberg( -imelt )
       !
       IF ( .not. limag_fly ) THEN
          !
          IF ( ALLOCATED(AKeri) ) DEALLOCATE(AKeri)
          !
          ! remove memory allocated for AKeri 
          imelt = ( upper_bnd - lower_bnd + 1 ) * maxval(nqfs(:)) * nbndfs**2 * ( 2 * nsiw(itemp) )
          CALL mem_size_eliashberg( -imelt )
          !
       ENDIF
       !
    ENDIF
    !
    RETURN
    !
    END SUBROUTINE sum_eliashberg_aniso_iaxis
    !
    !-----------------------------------------------------------------------
    SUBROUTINE analytic_cont_aniso_iaxis_to_raxis( itemp, iter, conv ) 
    !-----------------------------------------------------------------------
    !
    ! This routine does the analytic continuation of the anisotropic Eliashberg equations 
    ! from the imaginary-axis to the real axis
    ! reference F. Marsiglio, M. Schossmann, and J. Carbotte, Phys. Rev. B 37, 4965 (1988)
    !
    USE kinds,         ONLY : DP
    USE io_global,     ONLY : stdout
    USE phcom,         ONLY : nmodes
    USE elph2,         ONLY : wqf, wf
    USE epwcom,        ONLY : nqstep, degaussq, nsiter, conv_thr_racon, fsthick, & 
                              lpade, eps_acustic
    USE eliashbergcom, ONLY : nsw, estemp, dwsph, ws, wsph, gap, Agap, Gp, Gm, ADsumi, AZsumi, &                           
                              Delta, Znorm, ADelta, ADeltap, AZnorm, AZnormp, g2, lacon_fly, & 
                              a2fij, wkfs, dosef, ixkqf, ixqfs, nqfs, w0g, nkfs, nbndfs, ef0, ekfs
    USE superconductivity_iso, ONLY : gamma_acont
    USE superconductivity, ONLY : mem_size_eliashberg, eliashberg_memlt_aniso_acon
    USE constants_epw, ONLY : pi, ci, zero, czero, cone
    USE io_global,     ONLY : ionode_id
    USE mp_global,     ONLY : inter_pool_comm
    USE mp_world,      ONLY : mpime
    USE mp,            ONLY : mp_bcast, mp_barrier, mp_sum
    USE io_eliashberg, ONLY : eliashberg_write_raxis
    USE division,      ONLY : fkbounds
    ! 
    IMPLICIT NONE
    !
    INTEGER, INTENT (in) :: itemp
    !! Counter on temperature index
    INTEGER, INTENT (in) :: iter
    !! Counter on the iteration number
    LOGICAL, INTENT (inout) :: conv
    !! True if the calculation is converged
    !
    ! Local variables
    INTEGER :: i, iw, iwp
    !! Counter on frequency real-axis
    INTEGER :: ik
    !! Counter on k-points
    INTEGER :: iq
    !! Counter on q-points for which k+sign*q is within the Fermi shell
    INTEGER :: iq0
    !! Index of iq on full q-mesh
    INTEGER :: lower_bnd
    !! Lower bounds index after k paral
    INTEGER :: upper_bnd
    !! Upper bounds index after k paral
    INTEGER :: ibnd
    !! Counter on bands at k
    INTEGER :: jbnd
    !! Counter on bands at k+q
    INTEGER :: iwph
    !! Counter on frequency in a2f
    INTEGER :: imode
    !! Counter on phonon modes
    INTEGER :: imelt
    !! Counter memory

    REAL(kind=DP) :: rgammap
    !! -bose_einstein( w' ) - fermi_dirac(  w + w' )
    REAL(kind=DP) :: rgammam
    !! bose_einstein( w' ) + fermi_dirac( -w + w' )
    REAL(DP) :: absdelta, reldelta, errdelta
    !! Errors in Delta
    REAL(kind=DP) :: weight, a2f_
    !! Temporary variables
    REAL(kind=DP), EXTERNAL :: w0gauss
    !!
    COMPLEX(kind=DP) :: esqrt, root
    !! Temporary variables
    COMPLEX(DP), ALLOCATABLE, SAVE :: Deltaold(:)
    !! gap
    CHARACTER (len=256) :: cname
    !! character in file name
    !
    ! SP: Need initialization
    a2f_ = zero
    !
    IF ( iter .eq. 1 ) THEN
       !
       ! get the size of required allocated memory for 
       ! Delta, Znorm, Deltaold, ADelta, ADeltap, AZnorm, AZnormp, Gp, Gm 
       IF ( lpade ) THEN                 
          imelt = 2 * ( 1 + 2 * nbndfs * nkfs ) * nsw + 2 * nqstep * nsw
       ELSE
          imelt = 2 * ( 3 + 4 * nbndfs * nkfs ) * nsw + 2 * nqstep * nsw
       ENDIF
       CALL mem_size_eliashberg( imelt )
       !
       IF ( .not. ALLOCATED(Delta) )    ALLOCATE( Delta(nsw) )
       IF ( .not. ALLOCATED(Znorm) )    ALLOCATE( Znorm(nsw) )
       IF ( .not. ALLOCATED(ADelta) )   ALLOCATE( ADelta(nbndfs,nkfs,nsw) )
       IF ( .not. ALLOCATED(ADeltap) )  ALLOCATE( ADeltap(nbndfs,nkfs,nsw) )
       IF ( .not. ALLOCATED(AZnorm) )   ALLOCATE( AZnorm(nbndfs,nkfs,nsw) )
       IF ( .not. ALLOCATED(AZnormp) )  ALLOCATE( AZnormp(nbndfs,nkfs,nsw) )
       IF ( .not. ALLOCATED(Deltaold) ) ALLOCATE( Deltaold(nsw) )
       ADeltap(:,:,:) = czero
       AZnormp(:,:,:) = cone
       Deltaold(:) = czero
       IF ( lpade ) THEN 
          ADeltap(:,:,:) = ADelta(:,:,:)
          Deltaold(:) = Delta(:)
       ELSE
          DO ik = 1, nkfs
             DO ibnd = 1, nbndfs
                IF ( abs( ekfs(ibnd,ik) - ef0 ) .lt. fsthick ) THEN
                   ADeltap(ibnd,ik,:) = Agap(ibnd,ik,itemp) 
                ENDIF
             ENDDO ! ibnd
          ENDDO ! ik
          Deltaold(:) = gap(itemp)
       ENDIF
       !
       IF ( .not. ALLOCATED(Gp) ) ALLOCATE( Gp(nsw,nqstep) )
       IF ( .not. ALLOCATED(Gm) ) ALLOCATE( Gm(nsw,nqstep) )
       DO iw = 1, nsw ! loop over omega
          DO iwp = 1, nqstep ! loop over omega_prime
             CALL gamma_acont( ws(iw), ws(iwp), estemp(itemp), rgammap, rgammam )
             Gp(iw,iwp) = rgammap
             Gm(iw,iwp) = rgammam
          ENDDO
       ENDDO
       CALL kernel_aniso_iaxis_analytic_cont( itemp )
       CALL eliashberg_memlt_aniso_acon
       IF ( .not. lacon_fly ) CALL evaluate_a2fij
    ENDIF
    Delta(:) = czero
    Znorm(:) = czero
    ADelta(:,:,:) = czero
    AZnorm(:,:,:) = czero
    !
    CALL fkbounds( nkfs, lower_bnd, upper_bnd )
    !
    DO ik = lower_bnd, upper_bnd
       DO ibnd = 1, nbndfs
          IF ( abs( ekfs(ibnd,ik) - ef0 ) .lt. fsthick ) THEN
             DO iq = 1, nqfs(ik)
                ! iq0 - index of q-point on the full q-mesh
                iq0 = ixqfs(ik,iq)
                DO jbnd = 1, nbndfs
                   IF ( abs( ekfs(jbnd,ixkqf(ik,iq0)) - ef0 ) .lt. fsthick ) THEN
                      !
                      IF ( lacon_fly ) THEN ! evaluate a2fij on the fly
                         DO imode = 1, nmodes
                            IF ( wf(imode,iq0) .gt. eps_acustic ) THEN
                               DO iwph = 1, nqstep
                                  weight  = w0gauss( ( wsph(iwph) - wf(imode,iq0) ) / degaussq, 0 ) / degaussq
                                  a2f_ = weight * dosef * g2(ik,iq,ibnd,jbnd,imode)
                               ENDDO ! iwph
                            ENDIF ! wf
                         ENDDO ! imode
                      ENDIF
                      !
                      weight = wqf(iq) * w0g(jbnd,ixkqf(ik,iq0)) / dosef
                      DO iw = 1, nsw ! loop over omega
                         DO iwp = 1, nqstep ! loop over omega_prime
                            !
                            i = iw + iwp - 1
                            IF ( i .le. nsw ) THEN
                               root = sqrt(   AZnormp(jbnd,ixkqf(ik,iq0),i)**2.d0 & 
                                            * ( ws(i)**2.d0 - ADeltap(jbnd,ixkqf(ik,iq0),i)**2.d0 ) )
                               IF ( aimag(root) .lt. zero ) THEN 
                                  esqrt = AZnormp(jbnd,ixkqf(ik,iq0),i) / conjg(root)
                               ELSE  
                                  esqrt = AZnormp(jbnd,ixkqf(ik,iq0),i) / root
                               ENDIF
                               IF ( lacon_fly ) THEN 
                                  esqrt = esqrt * weight * Gp(iw,iwp) * a2f_
                               ELSE
                                  esqrt = esqrt * weight * Gp(iw,iwp) * a2fij(ik,iq,ibnd,jbnd,iwp) 
                               ENDIF
                               AZnorm(ibnd,ik,iw) = AZnorm(ibnd,ik,iw) - ws(i) * esqrt 
                               ADelta(ibnd,ik,iw) = ADelta(ibnd,ik,iw) - ADeltap(jbnd,ixkqf(ik,iq0),i) * esqrt
                            ENDIF
                            ! 
                            i = abs(iw - iwp) + 1
                            root = sqrt(   AZnormp(jbnd,ixkqf(ik,iq0),i)**2.d0 & 
                                         * ( ws(i)**2.d0 - ADeltap(jbnd,ixkqf(ik,iq0),i)**2.d0 ) )
                            IF ( aimag(root) .lt. zero ) THEN 
                               esqrt = AZnormp(jbnd,ixkqf(ik,iq0),i) / conjg(root)
                            ELSE  
                               esqrt = AZnormp(jbnd,ixkqf(ik,iq0),i) / root
                            ENDIF
                            esqrt = esqrt * weight * Gm(iw,iwp) * a2fij(ik,iq,ibnd,jbnd,iwp)
                            IF ( iw .lt. iwp ) THEN 
                               AZnorm(ibnd,ik,iw) = AZnorm(ibnd,ik,iw) - ws(i) * esqrt 
                            ELSE
                               AZnorm(ibnd,ik,iw) = AZnorm(ibnd,ik,iw) + ws(i) * esqrt 
                            ENDIF
                            ADelta(ibnd,ik,iw) = ADelta(ibnd,ik,iw) + ADeltap(jbnd,ixkqf(ik,iq0),i) * esqrt
                         ENDDO ! iwp
                      ENDDO ! iw
                   ENDIF
                ENDDO ! jbnd
             ENDDO ! iq
             DO iw = 1, nsw ! loop over omega
                AZnorm(ibnd,ik,iw) = - estemp(itemp) * AZsumi(ibnd,ik,iw) + ci * AZnorm(ibnd,ik,iw) * dwsph
                ADelta(ibnd,ik,iw) =   estemp(itemp) * ADsumi(ibnd,ik,iw) + ci * ADelta(ibnd,ik,iw) * dwsph
             ENDDO ! iw
          ENDIF
       ENDDO ! ibnd
    ENDDO ! ik
    !
    ! collect contributions from all pools 
    CALL mp_sum( AZnorm, inter_pool_comm )
    CALL mp_sum( ADelta, inter_pool_comm )
    CALL mp_barrier(inter_pool_comm)
    !
    IF (mpime .eq. ionode_id) THEN
      absdelta = zero
      reldelta = zero
      DO iw = 1, nsw ! loop over omega
         DO ik = 1, nkfs
            DO ibnd = 1, nbndfs
               IF ( abs( ekfs(ibnd,ik) - ef0 ) .lt. fsthick ) THEN
                  weight = 0.5d0 * wkfs(ik) * w0g(ibnd,ik) / dosef
                  Znorm(iw) = Znorm(iw) + weight * AZnorm(ibnd,ik,iw)
                  Delta(iw) = Delta(iw) + weight * ADelta(ibnd,ik,iw)
                  AZnorm(ibnd,ik,iw) = 1.d0 + pi * AZnorm(ibnd,ik,iw) / ws(iw)
                  ADelta(ibnd,ik,iw) = pi * ADelta(ibnd,ik,iw) / AZnorm(ibnd,ik,iw)
               ENDIF
            ENDDO ! ibnd                   
         ENDDO ! ik
         Znorm(iw) = 1.0d0 + pi * Znorm(iw) / ws(iw)
         Delta(iw) = pi * Delta(iw) / Znorm(iw)
         reldelta = reldelta + abs( Delta(iw) - Deltaold(iw) ) 
         absdelta = absdelta + abs( Delta(iw) ) 
      ENDDO ! iw
      errdelta = reldelta / absdelta
      Deltaold(:) = Delta(:)
      !
      WRITE(stdout,'(5x,a,i6,a,ES20.10,a,ES20.10,a,ES20.10)') 'iter = ', iter, & 
                   '   error = ', errdelta, '   Re[Znorm(1)] = ', real(Znorm(1)), & 
                   '   Re[Delta(1)] = ', real(Delta(1))
      !
      IF ( errdelta .lt. conv_thr_racon ) conv = .true.
      IF ( errdelta .lt. conv_thr_racon .OR. iter .eq. nsiter ) THEN
         cname = 'acon'
         CALL eliashberg_write_raxis( itemp, cname )
      ENDIF
      !
      IF ( conv .OR. iter .eq. nsiter ) THEN
         WRITE(stdout,'(5x,a,i6)') 'Convergence was reached in nsiter = ', iter
      ENDIF
      IF ( .not. conv .AND. iter .eq. nsiter ) THEN
         WRITE(stdout,'(5x,a,i6)') 'Convergence was not reached in nsiter = ', iter
         CALL errore('analytic_cont_aniso_iaxis_to_raxis','increase nsiter or reduce conv_thr_racon',1)
      ENDIF
    ENDIF
    CALL mp_bcast( Delta, ionode_id, inter_pool_comm )
    CALL mp_bcast( Znorm, ionode_id, inter_pool_comm )
    CALL mp_bcast( AZnorm, ionode_id, inter_pool_comm )
    CALL mp_bcast( gap, ionode_id, inter_pool_comm )
    CALL mp_bcast( Agap, ionode_id, inter_pool_comm )
    CALL mp_bcast( conv, ionode_id, inter_pool_comm )
    CALL mp_barrier(inter_pool_comm)
    !
    IF ( conv .OR. iter .eq. nsiter ) THEN
       !
       IF( ALLOCATED(Deltaold) ) DEALLOCATE(Deltaold)
       IF( ALLOCATED(Gp) )       DEALLOCATE(Gp)
       IF( ALLOCATED(Gm) )       DEALLOCATE(Gm)
       IF( ALLOCATED(ADsumi) )   DEALLOCATE(ADsumi)
       IF( ALLOCATED(AZsumi) )   DEALLOCATE(AZsumi)
       !
       ! remove memory allocated for Deltaold, Gp, Gm, ADsumi, AZsumi
       imelt = 2 * nsw + 2 * nqstep * nsw + 2 * ( upper_bnd - lower_bnd + 1 ) * nbndfs * nsw
       CALL mem_size_eliashberg( -imelt )
       !
       IF ( .not. lacon_fly ) THEN
          !
          IF ( ALLOCATED(a2fij) ) DEALLOCATE(a2fij)
          !
          ! remove memory allocated for a2fij
          imelt = ( upper_bnd - lower_bnd + 1 ) * maxval(nqfs(:)) * nbndfs**2 * nqstep
          CALL mem_size_eliashberg( -imelt )
          !
       ENDIF
       !
    ENDIF
    !
    RETURN
    !
    END SUBROUTINE analytic_cont_aniso_iaxis_to_raxis
    !
    !-----------------------------------------------------------------------
    SUBROUTINE pade_cont_aniso_iaxis_to_raxis( itemp, N, conv )
    !-----------------------------------------------------------------------
    !
    ! This routine uses pade approximants to continue the anisotropic Eliashberg equations 
    ! from the imaginary-axis to the real-axis
    !
    USE kinds,         ONLY : DP
    USE io_global,     ONLY : stdout
    USE epwcom,        ONLY : fsthick
    USE eliashbergcom, ONLY : nsw, ws, wsi, gap, Agap, Delta, Znorm, & 
                              ADelta, AZnorm, ADeltai, AZnormi, &              
                              wkfs, dosef, w0g, nkfs, nbndfs, ef0, ekfs
    USE superconductivity, ONLY : pade_coeff, pade_eval, mem_size_eliashberg
    USE constants_epw, ONLY : cone, ci, zero, czero
    USE io_global,     ONLY : ionode_id
    USE mp_global,     ONLY : inter_pool_comm
    USE mp_world,      ONLY : mpime
    USE mp,            ONLY : mp_bcast, mp_barrier, mp_sum
    USE io_eliashberg, ONLY : eliashberg_write_raxis
    USE division,      ONLY : fkbounds
    ! 
    IMPLICIT NONE
    !
    INTEGER, INTENT (in) :: itemp
    !! Counter on temperature index
    INTEGER, INTENT (in) :: N
    !! Maximum number of frequency
    LOGICAL, INTENT (inout) :: conv
    !! True if the calculation is converged
    !
    ! Local variable
    INTEGER :: iw
    !! Counter on frequency imag- and real-axis
    INTEGER :: ik
    !! Counter on k-points
    INTEGER :: lower_bnd
    !! Lower bounds index after k paral
    INTEGER :: upper_bnd
    !! Upper bounds index after k paral
    INTEGER :: ibnd
    !! Counter on bands
    INTEGER :: imelt
    !! Counter memory
    !
    REAL(DP) :: absdelta, reldelta, errdelta
    !! Errors in Delta
    REAL(DP) :: weight
    !! Temporary variable
    ! 
    ! arrays used in pade_coeff and pade_eval
    COMPLEX(DP), ALLOCATABLE :: a(:)
    !! a - pade coeff for Deltai
    COMPLEX(DP), ALLOCATABLE :: b(:)
    !! b - pade coeff for Znormi
    COMPLEX(DP), ALLOCATABLE :: z(:)
    !! z - frequency imag-axis
    COMPLEX(DP), ALLOCATABLE :: u(:)
    !! u - Deltai
    COMPLEX(DP), ALLOCATABLE :: v(:)
    !! v - Znormi
    COMPLEX(DP) :: omega
    !! frequency real-axis
    COMPLEX(DP) :: padapp
    !! Znorm or Delta on real-axis after pade_eval
    COMPLEX(DP), ALLOCATABLE :: Deltaold(:)
    !! gap
    CHARACTER (len=256) :: cname
    !! character in file name
    !
    ! get the size of required allocated memory for 
    ! a, b, z, u, v, Delta, Znorm, Deltaold, ADelta, AZnorm
    imelt = 2 * 5 * N + 2 * ( 3 + 2 * nbndfs * nkfs ) * nsw
    CALL mem_size_eliashberg( imelt )
    !
    IF ( .not. ALLOCATED(Delta) )    ALLOCATE( Delta(nsw) )
    IF ( .not. ALLOCATED(Znorm) )    ALLOCATE( Znorm(nsw) )
    IF ( .not. ALLOCATED(ADelta) )   ALLOCATE( ADelta(nbndfs,nkfs,nsw) )
    IF ( .not. ALLOCATED(AZnorm) )   ALLOCATE( AZnorm(nbndfs,nkfs,nsw) )
    IF ( .not. ALLOCATED(Deltaold) ) ALLOCATE( Deltaold(nsw) )
    IF ( .not. ALLOCATED(a) )        ALLOCATE( a(N) )
    IF ( .not. ALLOCATED(b) )        ALLOCATE( b(N) )
    IF ( .not. ALLOCATED(z) )        ALLOCATE( z(N) )
    IF ( .not. ALLOCATED(u) )        ALLOCATE( u(N) )
    IF ( .not. ALLOCATED(v) )        ALLOCATE( v(N) )
    Delta(:) = czero
    Znorm(:) = czero
    ADelta(:,:,:) = czero
    AZnorm(:,:,:) = czero
    Deltaold(:) = gap(itemp)
    absdelta = zero
    reldelta = zero
    a(:) = czero
    b(:) = czero
    z(:) = czero
    u(:) = czero
    v(:) = czero
    !
    CALL fkbounds( nkfs, lower_bnd, upper_bnd )
    !
    DO ik = lower_bnd, upper_bnd
       DO ibnd = 1, nbndfs
          IF ( abs( ekfs(ibnd,ik) - ef0 ) .lt. fsthick ) THEN
             DO iw = 1, N
                z(iw) = ci * wsi(iw)
                u(iw) = cone * ADeltai(ibnd,ik,iw) 
                v(iw) = cone * AZnormi(ibnd,ik,iw)
             ENDDO
             CALL pade_coeff( N, z, u, a )
             CALL pade_coeff( N, z, v, b )
             DO iw = 1, nsw
                omega = cone * ws(iw)
                CALL pade_eval( N, z, a, omega, padapp )
                ADelta(ibnd,ik,iw) = padapp
                CALL pade_eval( N, z, b, omega, padapp )
                AZnorm(ibnd,ik,iw) = padapp
             ENDDO
          ENDIF
       ENDDO ! ibnd
    ENDDO ! ik
    !
    ! collect contributions from all pools 
    CALL mp_sum( AZnorm, inter_pool_comm )
    CALL mp_sum( ADelta, inter_pool_comm )
    CALL mp_barrier(inter_pool_comm)
    !
    IF (mpime .eq. ionode_id) THEN
      DO iw = 1, nsw ! loop over omega
         DO ik = 1, nkfs
            DO ibnd = 1, nbndfs
               IF ( abs( ekfs(ibnd,ik) - ef0 ) .lt. fsthick ) THEN
                  weight = 0.5d0 * wkfs(ik) * w0g(ibnd,ik) / dosef
                  Znorm(iw) = Znorm(iw) + weight * AZnorm(ibnd,ik,iw)
                  Delta(iw) = Delta(iw) + weight * ADelta(ibnd,ik,iw)
               ENDIF
            ENDDO ! ibnd                   
         ENDDO ! ik
         reldelta = reldelta + abs( Delta(iw) - Deltaold(iw) )
         absdelta = absdelta + abs( Delta(iw) )
      ENDDO ! iw
      errdelta = reldelta / absdelta
      !
      IF ( errdelta .gt. zero ) THEN
         conv = .true.
         WRITE(stdout,'(5x,a,i6,a,ES20.10,a,ES20.10,a,ES20.10)') 'pade = ', N, & 
                '   error = ', errdelta, '   Re[Znorm(1)] = ', real(Znorm(1)), & 
                '   Re[Delta(1)] = ', real(Delta(1))
         cname = 'pade'
         CALL eliashberg_write_raxis( itemp, cname )
      ENDIF
    ENDIF
    CALL mp_bcast( Delta, ionode_id, inter_pool_comm )
    CALL mp_bcast( Znorm, ionode_id, inter_pool_comm )
    CALL mp_bcast( gap, ionode_id, inter_pool_comm )
    CALL mp_bcast( Agap, ionode_id, inter_pool_comm )
    CALL mp_bcast( conv, ionode_id, inter_pool_comm )
    CALL mp_barrier(inter_pool_comm)
    !
    IF( ALLOCATED(Deltaold) ) DEALLOCATE(Deltaold)
    IF( ALLOCATED(a) )        DEALLOCATE(a)
    IF( ALLOCATED(b) )        DEALLOCATE(b)
    IF( ALLOCATED(z) )        DEALLOCATE(z)
    IF( ALLOCATED(u) )        DEALLOCATE(u)
    IF( ALLOCATED(v) )        DEALLOCATE(v)
    !
    ! remove memory allocated for Deltaold, a, b, z, u, v
    imelt = 2 * ( nsw + 5 * N )
    CALL mem_size_eliashberg( -imelt )
    !
    RETURN
    !
    END SUBROUTINE pade_cont_aniso_iaxis_to_raxis
    !
    !-----------------------------------------------------------------------
    SUBROUTINE kernel_aniso_iaxis( itemp )
    !-----------------------------------------------------------------------
    !!  
    !! Compute kernels K_{+}(ik,iq,ibnd,jbnd;n,n',T) and K_{-}(ik,iq,ibnd,jbnd;n,n',T)
    !! and store them in memory
    !!
    USE kinds,         ONLY : DP
    USE epwcom,        ONLY : fsthick
    USE eliashbergcom, ONLY : nkfs, nbndfs, nsiw, estemp, AKeri, ekfs, ef0, ixkqf, ixqfs, nqfs
    USE constants_epw, ONLY : pi, zero 
    USE division,      ONLY : fkbounds
    ! 
    IMPLICIT NONE
    !
    INTEGER, INTENT (in) :: itemp
    !! Counter on temperature
    !
    ! Local variables
    INTEGER :: iw, n
    !! Counter on frequency imag-axis
    INTEGER :: ik
    !! Counter on k-points
    INTEGER :: iq
    !! Counter on q-points for which k+sign*q is within the Fermi shell
    INTEGER :: iq0
    !! Index of iq on full q-mesh
    INTEGER :: lower_bnd
    !! Lower bounds index after k paral
    INTEGER :: upper_bnd
    !! Upper bounds index after k paral
    INTEGER :: ibnd
    !! Counter on bands
    INTEGER :: jbnd
    !! Counter on bands
    !
    REAL(DP) :: omega
    !! frequency imag-axis
    REAL(DP) :: lambda_eph
    !! electron-phonon coupling
    !
    CALL fkbounds( nkfs, lower_bnd, upper_bnd )
    IF ( .not. ALLOCATED(AKeri) ) ALLOCATE( AKeri(lower_bnd:upper_bnd,maxval(nqfs(:)),nbndfs,nbndfs,2*nsiw(itemp)) )
    AKeri(:,:,:,:,:) = zero
    !
    ! RM - if lambdar_aniso_ver2 is used then one needs to CALL evaluate_a2fij
    !
    DO ik = lower_bnd, upper_bnd
       DO ibnd = 1, nbndfs
          IF ( abs( ekfs(ibnd,ik) - ef0 ) .lt. fsthick ) THEN
             DO iq = 1, nqfs(ik)
                ! iq0 - index of q-point on the full q-mesh
                iq0 = ixqfs(ik,iq)
                DO jbnd = 1, nbndfs
                   IF ( abs( ekfs(jbnd,ixkqf(ik,iq0)) - ef0 ) .lt. fsthick ) THEN
                      DO iw = 1, 2*nsiw(itemp)
                         n = iw - 1
                         omega = dble(2*n) * pi * estemp(itemp)
                         CALL lambdar_aniso_ver1( ik, iq, ibnd, jbnd, omega, lambda_eph )
                         !CALL lambdar_aniso_ver2( ik, iq, ibnd, jbnd, omega, lambda_eph )
                         AKeri(ik,iq,ibnd,jbnd,iw) = lambda_eph
                      ENDDO ! iw
                   ENDIF
                ENDDO ! jbnd
             ENDDO ! iq
          ENDIF
       ENDDO ! ibnd
    ENDDO ! ik
    !
    RETURN
    !
    END SUBROUTINE kernel_aniso_iaxis
    !
    !-----------------------------------------------------------------------
    SUBROUTINE lambdar_aniso_ver1( ik, iq, ibnd, jbnd, omega, lambda_eph )
    !-----------------------------------------------------------------------
    !
    ! computes lambda(ik,iq,ibnd,jbnd;n-n')   
    ! reference H. Choi et. al, Physica C 385, 66 (2003)
    !
    USE kinds,         ONLY : DP
    USE phcom,         ONLY : nmodes
    USE elph2,         ONLY : wf
    USE epwcom,        ONLY : eps_acustic
    USE eliashbergcom, ONLY : ixqfs, g2, dosef
    USE constants_epw, ONLY : zero
    ! 
    IMPLICIT NONE
    !
    INTEGER, INTENT (in) :: ik
    !! Counter on k-points
    INTEGER, INTENT (in) :: iq
    !! Counter on q-points for which k+sign*q is within the Fermi shell
    INTEGER, INTENT (in) :: ibnd
    !! Counter on bands
    INTEGER, INTENT (in) :: jbnd
    !! Counter on bands
    !
    REAL(DP), INTENT (in) :: omega
    !! frequency on imag-axis
    REAL(DP), INTENT (out) :: lambda_eph
    !! electron-phonon coupling lambda_ij(k,k+q;n-n')
    !
    ! Local variables
    INTEGER :: iq0
    !! Index of iq on full q-mesh
    INTEGER :: imode
    !! Counter on phonon modes
    !
    ! iq0 - index of q-point on the full q-mesh
    iq0 = ixqfs(ik,iq)
    lambda_eph = zero
    DO imode = 1, nmodes  ! loop over frequency modes
       IF ( wf(imode,iq0) .gt. eps_acustic ) THEN 
          lambda_eph = lambda_eph + g2(ik,iq,ibnd,jbnd,imode) * wf(imode,iq0) & 
                     / ( wf(imode,iq0)**2.d0 + omega**2.d0 )
       ENDIF
    ENDDO 
    lambda_eph = 2.d0 * lambda_eph * dosef
    !
    RETURN
    !
    END SUBROUTINE lambdar_aniso_ver1
    !
    !-----------------------------------------------------------------------
    SUBROUTINE lambdar_aniso_ver2( ik, iq, ibnd, jbnd, omega, lambda_eph )
    !-----------------------------------------------------------------------
    !
    ! computes lambda(ik,iq,ibnd,jbnd;n-n')   
    ! reference H. Choi et. al, Physica C 385, 66 (2003)
    !
    USE kinds,         ONLY : DP
    USE epwcom,        ONLY : nqstep
    USE eliashbergcom, ONLY : a2fij, wsph, dwsph
    USE constants_epw, ONLY : zero
    !                 
    IMPLICIT NONE        
    !
    INTEGER, INTENT (in) :: ik
    !! Counter on k-points
    INTEGER, INTENT (in) :: iq
    !! Counter on q-points for which k+sign*q is within the Fermi shell
    INTEGER, INTENT (in) :: ibnd
    !! Counter on bands
    INTEGER, INTENT (in) :: jbnd
    !! Counter on bands
    REAL(DP), INTENT (in) :: omega
    !! frequency on imag-axis
    REAL(DP), INTENT (out) :: lambda_eph
    !! electron-phonon coupling lambda_ij(k,k+q;n-n')
    !
    ! Local variables
    INTEGER :: iwph
    !! Counter on frequency
    !                          
    lambda_eph = zero
    DO iwph = 1, nqstep
       lambda_eph = lambda_eph + wsph(iwph) * a2fij(ik,iq,ibnd,jbnd,iwph) & 
                  / ( wsph(iwph)**2.d0 + omega**2.d0 )
    ENDDO 
    lambda_eph = 2.d0 * lambda_eph * dwsph
    !
    RETURN
    !
    END SUBROUTINE lambdar_aniso_ver2
    !
    !-----------------------------------------------------------------------
    SUBROUTINE kernel_aniso_iaxis_analytic_cont( itemp )
    !-----------------------------------------------------------------------
    !!  
    !! computes kernels K_{+}(w,iw_n,T) and K_{-}(w,iw_n,T)
    !! reference F. Masiglio, M. Schossmann, and J. Carbotte, PRB 37, 4965 (1988)
    !!
    USE kinds,         ONLY : DP
    USE elph2,         ONLY : wqf
    USE epwcom,        ONLY : muc, fsthick
    USE eliashbergcom, ONLY : nsw, nsiw, ws, wsi, ADeltai, nkfs, nbndfs, dosef, ixkqf, ixqfs, nqfs, & 
                              w0g, ekfs, ef0, ADsumi, AZsumi
    USE superconductivity, ONLY : mem_size_eliashberg
    USE constants_epw, ONLY : zero
    USE mp,            ONLY : mp_bcast, mp_barrier, mp_sum
    USE division,      ONLY : fkbounds
    !      
    IMPLICIT NONE
    !
    INTEGER, INTENT (in) :: itemp
    !! Counter on temperature index
    !
    ! Local variables
    INTEGER :: iw
    !! Counter on frequency real-axis
    INTEGER :: iwp
    !! Counter on frequency imag-axis
    INTEGER :: ik
    !! Counter on k-points
    INTEGER :: iq
    !! Counter on q-points for which k+sign*q is within the Fermi shell
    INTEGER :: iq0
    !! Index of iq on full q-mesh
    INTEGER :: lower_bnd
    !! Lower bounds index after k paral
    INTEGER :: upper_bnd
    !! Upper bounds index after k paral
    INTEGER :: ibnd
    !! Counter on bands
    INTEGER :: jbnd
    !! Counter on bands
    INTEGER :: imelt
    !! Counter memory
    ! 
    REAL(DP) :: esqrt, kernelp, kernelm, weight
    REAL(DP), ALLOCATABLE :: wesqrt(:,:,:), desqrt(:,:,:)
    !! Temporary working variables
    !
    COMPLEX(DP) :: lambda_eph
    !! electron-phonon coupling lambda_ij(k,k+q;w-iw_n)
    !
    CALL fkbounds( nkfs, lower_bnd, upper_bnd )
    !
    ! get memory size required for wesqrt, desqrt, ADsumi, AZsumi
    imelt = 2 * nbndfs * nkfs * nsiw(itemp) + 2 * ( upper_bnd - lower_bnd + 1 ) * nbndfs * nsw 
    CALL mem_size_eliashberg( imelt )
    !
    IF ( .not. ALLOCATED(wesqrt) ) ALLOCATE( wesqrt(nbndfs,nkfs,nsiw(itemp)) )
    IF ( .not. ALLOCATED(desqrt) ) ALLOCATE( desqrt(nbndfs,nkfs,nsiw(itemp)) )
    !
    DO ik = lower_bnd, upper_bnd
       IF ( .not. ALLOCATED(ADsumi) ) ALLOCATE( ADsumi(nbndfs,lower_bnd:upper_bnd,nsw) )
       IF ( .not. ALLOCATED(AZsumi) ) ALLOCATE( AZsumi(nbndfs,lower_bnd:upper_bnd,nsw) )
    ENDDO
    ADsumi(:,:,:) = zero
    AZsumi(:,:,:) = zero
    !
    ! RM - if lambdai_aniso_ver2 is used then one needs to CALL evaluate_a2fij
    !
    DO ik = lower_bnd, upper_bnd
       DO ibnd = 1, nbndfs
          IF ( abs( ekfs(ibnd,ik) - ef0 ) .lt. fsthick ) THEN
             DO iq = 1, nqfs(ik)
                ! iq0 - index of q-point on the full q-mesh
                iq0 = ixqfs(ik,iq)
                DO jbnd = 1, nbndfs
                   IF ( abs( ekfs(jbnd,ixkqf(ik,iq0)) - ef0 ) .lt. fsthick ) THEN
                      weight = wqf(iq) * w0g(jbnd,ixkqf(ik,iq0)) / dosef
                      DO iw = 1, nsw ! loop over omega 
                         DO iwp = 1, nsiw(itemp) ! loop over iw_n
                            CALL lambdai_aniso_ver1( ik, iq, ibnd, jbnd, ws(iw), wsi(iwp), lambda_eph )
                            !CALL lambdai_aniso_ver2( ik, iq, ibnd, jbnd, ws(iw), wsi(iwp), lambda_eph )
                            kernelp = 2.d0 * real(lambda_eph)
                            kernelm = 2.d0 * aimag(lambda_eph)
                            IF ( iw .eq. 1 ) THEN
                               esqrt = 1.d0 / sqrt( wsi(iwp)**2.d0 + ADeltai(jbnd,ixkqf(ik,iq0),iwp)**2.d0 )
                               wesqrt(jbnd,ixkqf(ik,iq0),iwp) =  wsi(iwp) * esqrt
                               desqrt(jbnd,ixkqf(ik,iq0),iwp) =  ADeltai(jbnd,ixkqf(ik,iq0),iwp) * esqrt
                            ENDIF
                            AZsumi(ibnd,ik,iw) = AZsumi(ibnd,ik,iw) & 
                                               + weight * wesqrt(jbnd,ixkqf(ik,iq0),iwp) * kernelm 
                            ADsumi(ibnd,ik,iw) = ADsumi(ibnd,ik,iw) & 
                                               + weight * desqrt(jbnd,ixkqf(ik,iq0),iwp) * ( kernelp - 2.d0 * muc ) 
                         ENDDO ! iwp
                      ENDDO ! iw
                   ENDIF
                ENDDO ! jbnd
             ENDDO ! iq
          ENDIF
       ENDDO ! ibnd
    ENDDO ! ik
    !
    IF( ALLOCATED(wesqrt) ) DEALLOCATE(wesqrt)
    IF( ALLOCATED(desqrt) ) DEALLOCATE(desqrt)
    !  
    ! remove memory allocated for wesqrt, desqrt
    imelt = 2 * nbndfs * nkfs * nsiw(itemp) 
    CALL mem_size_eliashberg ( -imelt )
    !
    RETURN
    !
    END SUBROUTINE kernel_aniso_iaxis_analytic_cont
    !                                                
    !-----------------------------------------------------------------------
    SUBROUTINE lambdai_aniso_ver1( ik, iq, ibnd, jbnd, omega, omegap, lambda_eph )
    !-----------------------------------------------------------------------
    !!
    !! computes lambda_ij(k,k+q;w-iw_n) 
    !! reference F. Masiglio, M. Schossmann, and J. Carbotte, PRB 37, 4965 (1988)
    !!
    USE kinds,         ONLY : DP
    USE phcom,         ONLY : nmodes
    USE elph2,         ONLY : wf
    USE epwcom,        ONLY : eps_acustic
    USE eliashbergcom, ONLY : ixqfs, g2, dosef
    USE constants_epw, ONLY : ci, czero
    !     
    IMPLICIT NONE
    !
    INTEGER, INTENT (in) :: ik
    !! Counter on k-points
    INTEGER, INTENT (in) :: iq
    !! Counter on q-points for which k+sign*q is within the Fermi shell
    INTEGER, INTENT (in) :: ibnd
    !! Counter on bands at k
    INTEGER, INTENT (in) :: jbnd
    !! Counter on bands at k+q
    REAL(DP), INTENT (in) :: omega
    !! frequency w at point iw on real-axis
    REAL(DP), INTENT (in) :: omegap
    !! frequency w_n at point iwp on imag-axis
    COMPLEX(DP), INTENT (out) :: lambda_eph
    !! electron-phonon coupling lambda_ij(k,k+q;w-iw_n)
    !
    ! Local variables
    INTEGER :: iq0
    !! Index of iq on full q-mesh
    INTEGER :: imode
    !! Counter on phonon modes
    !
    ! iq0 - index of q-point on the full q-mesh
    iq0 = ixqfs(ik,iq)
    lambda_eph = czero
    DO imode = 1, nmodes  ! loop over frequency modes
       IF ( wf(imode,iq0) .gt. eps_acustic ) THEN 
          lambda_eph = lambda_eph +  g2(ik,iq,ibnd,jbnd,imode) * wf(imode,iq0) & 
                     / ( wf(imode,iq0)**2.d0 - (omega - ci*omegap)**2.d0 )
       ENDIF
    ENDDO ! iwph
    lambda_eph = 2.d0 * lambda_eph * dosef
    !
    RETURN
    !
    END SUBROUTINE lambdai_aniso_ver1
    !
    !-----------------------------------------------------------------------               
    SUBROUTINE lambdai_aniso_ver2( ik, iq, ibnd, jbnd, omega, omegap, lambda_eph )
    !-----------------------------------------------------------------------
    !!
    !! computes lambda(w-iw_n)   
    !! reference F. Masiglio, M. Schossmann, and J. Carbotte, PRB 37, 4965 (1988)
    !!
    USE kinds,         ONLY : DP
    USE epwcom,        ONLY : nqstep
    USE eliashbergcom, ONLY : a2fij, dwsph, wsph
    USE constants_epw, ONLY : ci, czero
    !     
    IMPLICIT NONE
    ! 
    INTEGER, INTENT (in) :: ik
    !! Counter on k-points
    INTEGER, INTENT (in) :: iq
    !! Counter on q-points for which k+sign*q is within the Fermi shell
    INTEGER, INTENT (in) :: ibnd
    !! Counter on bands at k
    INTEGER, INTENT (in) :: jbnd
    !! Counter on bands at k+q
    REAL(DP), INTENT (in) :: omega
    !! frequency w at point iw on real-axis
    REAL(DP), INTENT (in) :: omegap
    !! frequency w_n at point iwp on imag-axis
    COMPLEX(DP), INTENT (out) :: lambda_eph
    !! electron-phonon coupling lambda_ij(k,k+q;w-iw_n)
    !
    ! Local variables
    INTEGER :: iwph
    !! Counter on frequency
    !
    lambda_eph = czero
    DO iwph = 1, nqstep
       lambda_eph = lambda_eph + wsph(iwph) * a2fij(ik,iq,ibnd,jbnd,iwph) &
                  / ( wsph(iwph)**2.d0 - (omega - ci*omegap)**2.d0 )
    ENDDO ! iwph
    lambda_eph = 2.d0 * lambda_eph * dwsph
    !
    RETURN
    !
    END SUBROUTINE lambdai_aniso_ver2
    !
    !-----------------------------------------------------------------------
    SUBROUTINE evaluate_a2fij
    !-----------------------------------------------------------------------
    !!
    !! computes the anisotropic spectral function a2F(k,k',w) 
    !!
    USE kinds,         ONLY : DP
    USE phcom,         ONLY : nmodes
    USE elph2,         ONLY : wf
    USE epwcom,        ONLY : fsthick, eps_acustic, nqstep, degaussq
    USE eliashbergcom, ONLY : nkfs, nbndfs, g2, a2fij, ixkqf, ixqfs, nqfs, ekfs, ef0, & 
                              dosef, wsph
    USE constants_epw, ONLY : zero
    USE division,      ONLY : fkbounds
    ! 
    IMPLICIT NONE
    !
    ! Local variables
    INTEGER :: ik
    !! Counter on k-points
    INTEGER :: iq
    !! Counter on q-points for which k+sign*q is within the Fermi shell
    INTEGER :: ibnd
    !! Counter on bands at k
    INTEGER :: jbnd
    !! Counter on bands at k+q
    INTEGER :: iq0
    !! Index of iq on full q-mesh
    INTEGER :: lower_bnd
    !! Lower bounds index after k paral
    INTEGER :: upper_bnd
    !! Upper bounds index after k paral
    INTEGER :: imode
    !! Counter on phonon modes
    INTEGER :: iwph
    !! Counter on frequency
    !
    REAL(DP) :: weight
    !! Temporary variable 
    REAL(DP), EXTERNAL :: w0gauss
    !
    CALL fkbounds( nkfs, lower_bnd, upper_bnd )
    IF ( .not. ALLOCATED(a2fij) ) ALLOCATE(a2fij(lower_bnd:upper_bnd,maxval(nqfs(:)),nbndfs,nbndfs,nqstep))
    a2fij(:,:,:,:,:) = zero
    !
    DO ik = lower_bnd, upper_bnd 
       DO ibnd = 1, nbndfs
          IF ( abs( ekfs(ibnd,ik) - ef0 ) .lt. fsthick ) THEN
             DO iq = 1, nqfs(ik)
                ! iq0 - index of q-point on the full q-mesh
                iq0 = ixqfs(ik,iq)
                DO jbnd = 1, nbndfs
                   IF ( abs( ekfs(jbnd,ixkqf(ik,iq0)) - ef0 ) .lt. fsthick ) THEN
                      DO imode = 1, nmodes
                         IF ( wf(imode,iq0) .gt. eps_acustic ) THEN 
                            DO iwph = 1, nqstep
                               weight  = w0gauss( ( wsph(iwph) - wf(imode,iq0) ) / degaussq, 0 ) / degaussq
                               a2fij(ik,iq,ibnd,jbnd,iwph) = a2fij(ik,iq,ibnd,jbnd,iwph) &
                                                           + weight * dosef * g2(ik,iq,ibnd,jbnd,imode)
                            ENDDO ! iwph
                         ENDIF ! wf
                      ENDDO ! imode
                   ENDIF
                ENDDO ! jbnd
             ENDDO ! iq
          ENDIF
       ENDDO ! ibnd
    ENDDO ! ik
    !
    RETURN
    !
    END SUBROUTINE evaluate_a2fij
    ! 
    !-----------------------------------------------------------------------
    SUBROUTINE evaluate_a2f_lambda
    !-----------------------------------------------------------------------
    !
    ! computes the isotropic spectral function a2F(w), total lambda, and 
    ! distribution of lambda
    !
    USE kinds,         ONLY : DP
    USE io_global,     ONLY : stdout
    USE io_epw,        ONLY : iua2ffil, iudosfil, iufillambda, iufillambdaFS
    USE io_files,      ONLY : prefix
    USE phcom,         ONLY : nmodes
    USE cell_base,     ONLY : bg
    USE control_flags, ONLY : iverbosity
    USE elph2,         ONLY : nqtotf, wqf, wf
    USE epwcom,        ONLY : fsthick, eps_acustic, nqstep, degaussq, delta_qsmear, nqsmear, & 
                              degaussw, nkf1, nkf2, nkf3
    USE eliashbergcom, ONLY : nkfs, nbndfs, g2, ixkqf, ixqfs, nqfs, w0g, ekfs, ef0, dosef, wsph, &
                              wkfs, dwsph, a2f_iso, ixkff
    USE constants_epw, ONLY : ryd2ev
    USE io_global,     ONLY : ionode_id
    USE mp_global,     ONLY : inter_pool_comm, my_pool_id, npool
    USE mp_world,      ONLY : mpime
    USE mp,            ONLY : mp_bcast, mp_barrier, mp_sum
    USE division,      ONLY : fkbounds
    ! 
    IMPLICIT NONE
    !
    INTEGER :: ik, iq, iq0, iwph, ibnd, jbnd, imode, lower_bnd, upper_bnd, &
         ismear, ibin, nbin, nbink, i, j, k
    REAL(DP) :: weight, weightq, l_sum, lambda_eph, lambda_max(npool), & 
         sigma, dbin, dbink, x1, x2, x3
    REAL(DP), ALLOCATABLE :: a2f(:,:), phdos(:,:), l_a2f(:), lambda_k(:,:), &
         lambda_k_bin(:), lambda_pairs(:), a2f_modeproj(:,:), phdos_modeproj(:,:)
    REAL(DP), EXTERNAL :: w0gauss
    CHARACTER (len=256) :: name1
    ! 
    ! This is only a quick fix since the subroutine was written for parallel execution - FG June 2014
#if ! defined(__MPI)
    npool = 1
    my_pool_id = 0
#endif
    !
    ! degaussq is read from the input file in meV and converted to Ryd in epw_readin.f90
    ! go from Ryd to eV
    degaussq = degaussq * ryd2ev
    delta_qsmear = delta_qsmear * ryd2ev
    !
    CALL fkbounds( nkfs, lower_bnd, upper_bnd )
    !
    IF ( .not. ALLOCATED(a2f_iso) ) ALLOCATE( a2f_iso(nqstep) )
    IF ( .not. ALLOCATED(a2f) )     ALLOCATE( a2f(nqstep,nqsmear) )
    IF ( .not. ALLOCATED(a2f_modeproj) ) ALLOCATE( a2f_modeproj(nmodes,nqstep) )
    a2f_iso(:) = 0.d0
    a2f(:,:) = 0.d0
    a2f_modeproj(:,:) = 0.d0
    !
    ! RM - the 0 index in k is required when printing out values of lambda_k 
    ! When the k-point is outside the Fermi shell, ixkff(ik)=0
    IF ( .not. ALLOCATED(lambda_k) ) ALLOCATE(lambda_k(0:nkfs,nbndfs))
    lambda_k(:,:) = 0.d0
    !
    l_sum = 0.d0
    lambda_max(:) = 0.d0
    DO ismear = 1, nqsmear
       sigma = degaussq + (ismear-1) * delta_qsmear
       DO ik = lower_bnd, upper_bnd
          DO ibnd = 1, nbndfs
             IF ( abs( ekfs(ibnd,ik) - ef0 ) .lt. fsthick ) THEN 
                DO iq = 1, nqfs(ik)
                   ! iq0 - index of q-point on the full q-mesh
                   iq0 = ixqfs(ik,iq)
                   DO jbnd = 1, nbndfs
                      IF ( abs( ekfs(jbnd,ixkqf(ik,iq0)) - ef0 ) .lt. fsthick ) THEN
                         weight = wkfs(ik) * wqf(iq) * w0g(ibnd,ik) * w0g(jbnd,ixkqf(ik,iq0))
                         lambda_eph = 0.d0
                         DO imode = 1, nmodes
                            IF ( wf(imode,iq0) .gt. eps_acustic ) THEN
                               IF ( ismear .eq. 1 ) THEN 
                                  lambda_eph = lambda_eph + g2(ik,iq,ibnd,jbnd,imode) / wf(imode,iq0)
                               ENDIF
                               DO iwph = 1, nqstep
                                  weightq  = w0gauss( ( wsph(iwph) - wf(imode,iq0) ) / sigma, 0 ) / sigma
                                  a2f(iwph,ismear) = a2f(iwph,ismear) + weight * weightq * g2(ik,iq,ibnd,jbnd,imode)
                                  IF ( ismear .eq. 1 ) THEN
                                     a2f_modeproj(imode,iwph) = a2f_modeproj(imode,iwph) +&
                                         weight * weightq * g2(ik,iq,ibnd,jbnd,imode)
                                  ENDIF
                               ENDDO ! iwph
                            ENDIF ! wf
                         ENDDO ! imode
                         IF ( ismear .eq. 1 .AND. lambda_eph .gt. 0.d0 ) THEN
                            l_sum = l_sum + weight * lambda_eph
                            weight = wqf(iq) * w0g(jbnd,ixkqf(ik,iq0)) 
                            lambda_k(ik,ibnd) = lambda_k(ik,ibnd) + weight * lambda_eph
                            IF ( lambda_eph .gt. lambda_max(my_pool_id+1) ) THEN
                               lambda_max(my_pool_id+1) = lambda_eph
                            ENDIF
                         ENDIF
                      ENDIF ! ekq
                   ENDDO ! jbnd
                ENDDO ! iq
             ENDIF ! ekk
          ENDDO ! ibnd
       ENDDO ! ik
    ENDDO ! ismear
    !
    a2f(:,:) = 0.5d0 * a2f(:,:) / dosef
    a2f_modeproj(:,:) = 0.5d0 * a2f_modeproj(:,:) / dosef
    l_sum = l_sum / dosef
    lambda_k(:,:) = 2.d0 * lambda_k(:,:)
    lambda_max(:) = 2.d0 * dosef * lambda_max(:)
    !
    ! collect contributions from all pools (sum over k-points)
    CALL mp_sum( l_sum, inter_pool_comm )
    CALL mp_sum( a2f, inter_pool_comm )
    CALL mp_sum( a2f_modeproj, inter_pool_comm )
    CALL mp_sum( lambda_max, inter_pool_comm )
    CALL mp_sum( lambda_k, inter_pool_comm )
    CALL mp_barrier(inter_pool_comm)
    !
    IF ( mpime .eq. ionode_id ) THEN
      !
      OPEN( unit = iua2ffil, file = TRIM(prefix)//".a2f", form = 'formatted')
      OPEN( unit = iudosfil, file = TRIM(prefix)//".phdos", form = 'formatted')
      !
      IF ( .not. ALLOCATED(phdos) )     ALLOCATE( phdos(nqstep,nqsmear) )
      IF ( .not. ALLOCATED(phdos_modeproj) ) ALLOCATE( phdos_modeproj(nmodes,nqstep) )
      phdos(:,:) = 0.d0
      phdos_modeproj(:,:) = 0.d0
      !
      DO ismear = 1, nqsmear
         sigma = degaussq + (ismear-1) * delta_qsmear
         DO iq = 1, nqtotf
            DO imode = 1, nmodes
               IF ( wf(imode,iq) .gt. eps_acustic ) THEN
                  DO iwph = 1, nqstep
                     weightq  = w0gauss( ( wsph(iwph) - wf(imode,iq)) / sigma, 0 ) / sigma
                     phdos(iwph,ismear) = phdos(iwph,ismear) + wqf(iq) * weightq
                     IF ( ismear .eq. 1 ) THEN
                        phdos_modeproj(imode,iwph) = phdos_modeproj(imode,iwph) + wqf(iq) * weightq
                     ENDIF
                  ENDDO ! iwph
               ENDIF ! wf
            ENDDO ! imode
         ENDDO ! iq
      ENDDO ! ismear
      !
      IF ( .not. ALLOCATED(l_a2f) )   ALLOCATE( l_a2f(nqsmear) )
      l_a2f(:) = 0.d0
      !
      DO ismear = 1, nqsmear
         DO iwph = 1, nqstep
            l_a2f(ismear) = l_a2f(ismear) + a2f(iwph,ismear) / wsph(iwph)
            ! wsph in meV (from eV) and phdos in states/meV (from states/eV)
            IF (ismear .eq. nqsmear) WRITE (iua2ffil,'(f12.7,15f12.7)') wsph(iwph)*1000.d0, a2f(iwph,:)
            IF (ismear .eq. nqsmear) WRITE (iudosfil,'(f12.7,15f15.7)') wsph(iwph)*1000.d0, phdos(iwph,:)/1000.d0
         ENDDO
         l_a2f(ismear) = 2.d0 * l_a2f(ismear) * dwsph
      ENDDO
      !
      WRITE(iua2ffil,*) "Integrated el-ph coupling"
      WRITE(iua2ffil,'("  #         ", 15f12.7)') l_a2f(:)
      WRITE(iua2ffil,*) "Phonon smearing (meV)" 
      WRITE(iua2ffil,'("  #         ", 15f12.7)') ( (degaussq+(ismear-1)*delta_qsmear)*1000.d0,ismear=1,nqsmear )
      WRITE(iua2ffil,'(" Electron smearing (eV)", f12.7)') degaussw
      WRITE(iua2ffil,'(" Fermi window (eV)", f12.7)') fsthick
      WRITE(iua2ffil,'(" Summed el-ph coupling ", f12.7)') l_sum
      CLOSE(iua2ffil)
      CLOSE(iudosfil)
      !
      a2f_iso(:) = a2f(:,1)
      OPEN( unit = iua2ffil, file = TRIM(prefix)//".a2f_iso", form = 'formatted')
      OPEN( unit = iudosfil, file = TRIM(prefix)//".phdos_proj", form = 'formatted')
      DO iwph = 1, nqstep
         ! wsph in meV (from eV) and phdos in states/meV (from states/eV)
         WRITE(iua2ffil,'(f12.7,100f12.7)') wsph(iwph)*1000.d0, a2f_iso(iwph), a2f_modeproj(:,iwph)
         WRITE(iudosfil,'(f12.7,100f15.7)') wsph(iwph)*1000.d0, phdos(iwph,1)/1000.d0, phdos_modeproj(:,iwph)/1000.d0
      ENDDO
      WRITE(iua2ffil,'(a,f18.7,a,f18.7)') 'lambda_int = ', l_a2f(1), '   lambda_sum = ',l_sum
      CLOSE(iua2ffil)
      CLOSE(iudosfil)
      !
      IF ( ALLOCATED(phdos) )          DEALLOCATE( phdos )
      IF ( ALLOCATED(phdos_modeproj) ) DEALLOCATE( phdos_modeproj )
      IF ( ALLOCATED(l_a2f) )          DEALLOCATE( l_a2f )
      !
    ENDIF
    !
    CALL mp_bcast( a2f_iso, ionode_id, inter_pool_comm )
    CALL mp_barrier(inter_pool_comm)
    !
    IF ( ALLOCATED(a2f) )            DEALLOCATE( a2f )
    IF ( ALLOCATED(a2f_modeproj) )   DEALLOCATE( a2f_modeproj )
    !
    nbink = int( 1.25d0 * maxval(lambda_k(:,:)) / 0.005d0 )
    dbink = 1.25d0 * maxval(lambda_k(:,:)) / dble(nbink)
    IF ( .not. ALLOCATED(lambda_k_bin) ) ALLOCATE ( lambda_k_bin(nbink) )
    lambda_k_bin(:) = 0.d0
    !
    !SP : Should be initialized
    nbin = 0
    dbin = 0.0_DP
    !
    IF ( iverbosity == 2 ) THEN
       nbin = int( 1.25d0 * maxval(lambda_max(:)) / 0.005d0 )
       dbin = 1.25d0 * maxval(lambda_max(:)) / dble(nbin)
       IF ( .not. ALLOCATED(lambda_pairs) ) ALLOCATE ( lambda_pairs(nbin) )
       lambda_pairs(:) = 0.d0
    ENDIF
    ! 
    WRITE(stdout,'(5x,a13,f21.7,a18,f21.7)') 'lambda_max = ', maxval(lambda_max(:)), '   lambda_k_max = ', maxval(lambda_k(:,:))
    WRITE(stdout,'(a)') ' '
    !
    lambda_k(:,:) = 0.d0
    DO ik = lower_bnd, upper_bnd
       DO ibnd = 1, nbndfs
          IF ( abs( ekfs(ibnd,ik) - ef0 ) .lt. fsthick ) THEN
             DO iq = 1, nqfs(ik)
                ! iq0 - index of q-point on the full q-mesh
                iq0 = ixqfs(ik,iq)
                DO jbnd = 1, nbndfs
                   IF ( abs( ekfs(jbnd,ixkqf(ik,iq0)) - ef0 ) .lt. fsthick ) THEN
                      weight = wqf(iq) * w0g(jbnd,ixkqf(ik,iq0)) / dosef
                      CALL lambdar_aniso_ver1( ik, iq, ibnd, jbnd, 0.d0, lambda_eph )
                      lambda_k(ik,ibnd) = lambda_k(ik,ibnd) +  weight * lambda_eph
                      IF ( iverbosity == 2 ) THEN
                         DO ibin = 1, nbin
                            sigma = 1.d0 * dbin
                            weight = w0gauss( ( lambda_eph - dble(ibin) * dbin ) / sigma, 0 ) / sigma
                            lambda_pairs(ibin) = lambda_pairs(ibin) + weight
                         ENDDO
                      ENDIF
                   ENDIF
                ENDDO ! jbnd
             ENDDO ! iq
             DO ibin = 1, nbink
                sigma = 1.d0 * dbink
                weight = w0gauss( ( lambda_k(ik,ibnd) - dble(ibin) * dbink ) / sigma, 0 ) / sigma
                lambda_k_bin(ibin) = lambda_k_bin(ibin) + weight
             ENDDO
          ENDIF
       ENDDO ! ibnd
    ENDDO ! ik
    !
    ! collect contributions from all pools 
    CALL mp_sum( lambda_k, inter_pool_comm )
    IF ( iverbosity .eq. 2 ) THEN  
      CALL mp_sum( lambda_pairs, inter_pool_comm )
    ENDIF
    CALL mp_sum( lambda_k_bin, inter_pool_comm )
    CALL mp_barrier(inter_pool_comm)
    !
    IF ( mpime .eq. ionode_id ) THEN
      !
      ! SP: Produced if user really wants it 
      IF ( iverbosity .eq. 2 ) THEN
        OPEN(unit = iufillambda, file = TRIM(prefix)//".lambda_aniso", form = 'formatted')
        WRITE(iufillambda,'(2a12,2a7)') '# enk-e0[eV]','  lambda_nk','# kpt','# band'
        DO ik = 1, nkfs
           DO ibnd = 1, nbndfs
              IF ( abs( ekfs(ibnd,ik) - ef0 ) .lt. fsthick ) THEN
                 WRITE(iufillambda,'(2f12.7,2i7)') ekfs(ibnd,ik) - ef0, lambda_k(ik,ibnd), ik, ibnd
              ENDIF
           ENDDO
        ENDDO
        CLOSE(iufillambda)
      ENDIF
      !
      OPEN(unit = iufillambda, file = TRIM(prefix)//".lambda_k_pairs", form = 'formatted')
      WRITE(iufillambda,'(a12,a30)') '# lambda_nk','  \rho(lambda_nk) scaled to 1'
      DO ibin = 1, nbink
        WRITE(iufillambda,'(2f21.7)') dbink*dble(ibin), lambda_k_bin(ibin)/maxval(lambda_k_bin(:))
      ENDDO
      CLOSE(iufillambda)
      !
      ! SP: Produced if user really wants it 
      IF ( iverbosity == 2 ) THEN  
        OPEN( unit = iufillambda, file = TRIM(prefix)//".lambda_pairs", form = 'formatted')
      WRITE(iufillambda,'(a12,a30)') "# lambda_nk,n'k'", "  \rho(lambda_nk,n'k') scaled to 1"
        DO ibin = 1, nbin
          WRITE(iufillambda,'(2f21.7)') dbin*dble(ibin), lambda_pairs(ibin)/maxval(lambda_pairs(:))
        ENDDO
        CLOSE(iufillambda)
      ENDIF
      !
      ! SP & RM: .cube file for VESTA plotting (only if iverbosity = 2)
      !
      ! RM - If the k-point is outside the Fermi shell,
      ! ixkff(ik)=0 and lambda_k(0,ibnd) = 0.0
      !
      IF ( iverbosity .eq. 2 ) THEN
        !
        DO ibnd = 1, nbndfs
          !
          IF ( ibnd < 10 ) THEN
            WRITE(name1,'(a,a8,i1,a5)') TRIM(prefix),'.lambda_', ibnd, '.cube'
          ELSEIF ( ibnd < 100 ) THEN
            WRITE(name1,'(a,a8,i2,a5)') TRIM(prefix),'.lambda_', ibnd, '.cube'
          ELSEIF( ibnd < 1000 ) THEN
            WRITE(name1,'(a,a8,i3,a5)') TRIM(prefix),'.lambda_', ibnd, '.cube'
          ELSE 
            CALL errore( 'eliashberg_setup', 'Too many bands ',1)  
          ENDIF  
          !  
          OPEN(iufillambdaFS, file=name1, form='formatted')
          WRITE(iufillambdaFS,*) 'Cubfile created from EPW calculation'
          WRITE(iufillambdaFS,*) 'lambda'
          WRITE(iufillambdaFS,'(i5,3f12.6)') 1, 0.0d0, 0.0d0, 0.0d0
          WRITE(iufillambdaFS,'(i5,3f12.6)') nkf1, (bg(i,1)/DBLE(nkf1),i=1,3)
          WRITE(iufillambdaFS,'(i5,3f12.6)') nkf2, (bg(i,2)/DBLE(nkf2),i=1,3)
          WRITE(iufillambdaFS,'(i5,3f12.6)') nkf3, (bg(i,3)/DBLE(nkf3),i=1,3)
          WRITE(iufillambdaFS,'(i5,4f12.6)') 1, 1.0d0, 0.0d0, 0.0d0, 0.0d0
          WRITE(iufillambdaFS,'(6f12.6)') ( lambda_k(ixkff(ik),ibnd), ik=1,nkf1*nkf2*nkf3 )
          CLOSE(iufillambdaFS)
        ENDDO
        !
      ENDIF
      !
      ! SP & RM : Write on file the lambda close to the Fermi surface along with 
      ! Cartesian coordinate, band index, energy distance from Fermi level
      ! and lambda value.
      !
      OPEN(unit = iufillambdaFS, file = TRIM(prefix)//".lambda_FS", form='formatted')
      WRITE(iufillambdaFS,'(a75)') '#               k-point                  Band Enk-Ef [eV]            lambda'
      DO i = 1, nkf1
         DO j = 1, nkf2
            DO k = 1, nkf3
               ik = k + (j-1)*nkf3 + (i-1)*nkf2*nkf3
               IF ( ixkff(ik) .gt. 0 ) THEN
                  DO ibnd = 1, nbndfs
                     ! SP: Here take a 0.2 eV interval around the FS.
                     IF ( abs( ekfs(ibnd,ixkff(ik)) - ef0 ) .lt. fsthick ) THEN
                     !IF ( abs( ekfs(ibnd,ixkff(ik)) - ef0 ) .lt. 0.2 ) THEN
                        x1 = bg(1,1)*(i-1)/nkf1+bg(1,2)*(j-1)/nkf2+bg(1,3)*(k-1)/nkf3
                        x2 = bg(2,1)*(i-1)/nkf1+bg(2,2)*(j-1)/nkf2+bg(2,3)*(k-1)/nkf3
                        x3 = bg(3,1)*(i-1)/nkf1+bg(3,2)*(j-1)/nkf2+bg(3,3)*(k-1)/nkf3
                        WRITE(iufillambdaFS,'(3f12.6,i8,f12.6,f24.15)') x1, x2, x3, ibnd, &
                                         ekfs(ibnd,ixkff(ik))-ef0, lambda_k(ixkff(ik),ibnd)
                     ENDIF
                  ENDDO ! ibnd
               ENDIF
            ENDDO  ! k
         ENDDO ! j
      ENDDO ! i
      CLOSE(iufillambdaFS)
    ENDIF
    CALL mp_barrier(inter_pool_comm)
    !
    IF ( ALLOCATED(lambda_k) )     DEALLOCATE(lambda_k)
    IF ( ALLOCATED(lambda_pairs) ) DEALLOCATE(lambda_pairs)
    IF ( ALLOCATED(lambda_k_bin) ) DEALLOCATE(lambda_k_bin)
    !
    RETURN
    !
    END SUBROUTINE evaluate_a2f_lambda
    ! 




    !
  END MODULE superconductivity_aniso
