  !                                                                            
  ! Copyright (C) 2010-2016 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino   
  ! Copyright (C) 2007-2009 Roxana Margine
  !                                                                            
  ! This file is distributed under the terms of the GNU General Public         
  ! License. See the file `LICENSE' in the root directory of the               
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .             
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
  USE constants_epw, ONLY : kelvin2eV, ci, pi
  USE io_global,     ONLY : ionode_id
  USE mp_global,     ONLY : inter_pool_comm
  USE mp,            ONLY : mp_bcast, mp_barrier
  USE mp_world,      ONLY : mpime
  ! 
  IMPLICIT NONE
  !
  INTEGER :: itemp, iter, N, ik, ibnd, imelt
  REAL(DP) :: tcpu, rdeltaout(nsw), rdeltain(nsw), cdeltaout(nsw), cdeltain(nsw)
  REAL(DP), EXTERNAL :: get_clock
  LOGICAL :: conv 
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
  ! input
  !
  ! itemp  - temperature point
  ! iter   - iteration number
  ! conv   - convergence flag 
  !
  ! output 
  !
  ! conv   - convergence flag 
  !
  USE kinds,         ONLY : DP
  USE io_global,     ONLY : stdout
  USE elph2,         ONLY : wqf
  USE epwcom,        ONLY : nsiter, nstemp, muc, conv_thr_iaxis, fsthick
  USE eliashbergcom, ONLY : nsiw, estemp, gap0, gap, Agap, wsi, AKeri, limag_fly, & 
                            NAZnormi, AZnormi, ADeltai, ADeltaip, NZnormi, Znormi, & 
                            Deltai, wsphmax, nkfs, nbndfs, dosef, ef0, ixkqf, ixqfs, & 
                            nqfs, wkfs, w0g, ekfs
  USE constants_epw, ONLY : pi  
  USE io_global,     ONLY : ionode_id
  USE mp_global,     ONLY : inter_pool_comm
  USE mp_world,      ONLY : mpime
  USE mp,            ONLY : mp_bcast, mp_barrier, mp_sum
  ! 
  IMPLICIT NONE
  !
  INTEGER  :: iw, iwp, itemp, iter, ik, iq, iq0, ibnd, jbnd, & 
             lower_bnd, upper_bnd, imelt
  REAL(DP) :: esqrt, absdelta, reldelta, errdelta, weight
  REAL(DP) :: kernelp, kernelm, lambdap, lambdam
  REAL(DP), ALLOCATABLE :: wesqrt(:,:,:), desqrt(:,:,:)
  REAL(DP), ALLOCATABLE, SAVE :: Deltaold(:)
  LOGICAL  :: conv
  !
  IF ( .not. ALLOCATED(wesqrt) ) ALLOCATE( wesqrt(nbndfs,nkfs,nsiw(itemp)) )
  IF ( .not. ALLOCATED(desqrt) ) ALLOCATE( desqrt(nbndfs,nkfs,nsiw(itemp)) )
  wesqrt(:,:,:) = 0.d0
  desqrt(:,:,:) = 0.d0
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
     gap(itemp) = 0.d0
     Agap(:,:,itemp) = 0.d0
     ADeltaip(:,:,:) = 0.d0
     !
     DO ik = 1, nkfs
        DO ibnd = 1, nbndfs
           IF ( abs( ekfs(ibnd,ik) - ef0 ) .lt. fsthick ) THEN
              DO iw = 1, nsiw(itemp)
                 IF ( wsi(iw) .lt. 2.d0*wsphmax ) THEN
                    ADeltaip(ibnd,ik,iw) = gap0
                 ELSE
                    ADeltaip(ibnd,ik,iw) = 0.d0
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
  Deltai(:) = 0.d0
  Znormi(:) = 0.d0
  NZnormi(:) = 0.d0
  ADeltai(:,:,:) = 0.d0
  AZnormi(:,:,:) = 0.d0
  NAZnormi(:,:,:) = 0.d0
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
    absdelta = 0.d0
    reldelta = 0.d0
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
  SUBROUTINE eliashberg_read_aniso_iaxis( itemp )
  !-----------------------------------------------------------------------
  !!  
  !! This routine reads from file the anisotropic Delta and Znorm on the imaginary-axis
  !! 
  !! input
  !!
  !! itemp  - temperature point
  !!
  !---------------------------------------------------------------------- 
  USE kinds,         ONLY : DP
  USE io_epw,        ONLY : iufilgap
  USE io_files,      ONLY : prefix
  USE control_flags, ONLY : iverbosity
  USE epwcom,        ONLY : nstemp, fsthick
  USE eliashbergcom, ONLY : nsiw, estemp, gap0, gap, Agap, wsi, NZnormi, Znormi, Deltai, & 
                            AZnormi, NAZnormi, ADeltai, nkfs, nbndfs, ef0, ekfs, &
                            dosef, wkfs, w0g
  USE constants_epw, ONLY : kelvin2eV
  USE io_global, ONLY : ionode_id
  USE mp_global, ONLY : inter_pool_comm
  USE mp_world,  ONLY : mpime
  USE mp,        ONLY : mp_bcast, mp_barrier, mp_sum
  ! 
  IMPLICIT NONE
  !
  INTEGER, INTENT(in) :: itemp
  !
  ! Local variables
  INTEGER :: iw
  !! Counter on frequency
  INTEGER :: ik
  !! Counter on k-poin
  INTEGER :: ibnd
  !! Counter on band
  INTEGER :: imelt
  !! Required allocation of memory
  INTEGER :: ios
  !! Status variables when reading a file
  REAL(DP) :: temp, eband, omega, weight
  REAL(DP) :: eps=1.0d-6
  CHARACTER (len=256) :: name1, word
  !
  ! get the size of required allocated memory 
  imelt = ( 1 + nbndfs * nkfs ) * nstemp + ( 3 + 3 * nbndfs * nkfs ) * nsiw(itemp)
  CALL mem_size_eliashberg( imelt )
  !
  IF ( .not. ALLOCATED(gap) )      ALLOCATE( gap(nstemp) )
  IF ( .not. ALLOCATED(Agap) )     ALLOCATE( Agap(nbndfs,nkfs,nstemp) )
  IF ( .not. ALLOCATED(Deltai) )   ALLOCATE( Deltai(nsiw(itemp)) )
  IF ( .not. ALLOCATED(Znormi) )   ALLOCATE( Znormi(nsiw(itemp)) )
  IF ( .not. ALLOCATED(NZnormi) )  ALLOCATE( NZnormi(nsiw(itemp)) )
  IF ( .not. ALLOCATED(ADeltai) )  ALLOCATE( ADeltai(nbndfs,nkfs,nsiw(itemp)) )
  IF ( .not. ALLOCATED(AZnormi) )  ALLOCATE( AZnormi(nbndfs,nkfs,nsiw(itemp)) )
  IF ( .not. ALLOCATED(NAZnormi) ) ALLOCATE( NAZnormi(nbndfs,nkfs,nsiw(itemp)) )
  gap(:) = 0.d0
  Agap(:,:,:) = 0.d0
  Deltai(:) = 0.d0
  Znormi(:) = 0.d0
  NZnormi(:) = 0.d0
  ADeltai(:,:,:) = 0.d0
  AZnormi(:,:,:) = 0.d0
  NAZnormi(:,:,:) = 0.d0
  !
  IF (mpime .eq. ionode_id) THEN     
    !   
    temp = estemp(itemp) / kelvin2eV
    ! anisotropic case
    IF ( temp .lt. 10.d0 ) THEN
       WRITE(name1,'(a,a13,f4.2)') TRIM(prefix),'.imag_aniso_0', temp
    ELSEIF ( temp .ge. 10.d0 ) THEN
       WRITE(name1,'(a,a12,f5.2)') TRIM(prefix),'.imag_aniso_', temp
    ENDIF 
    OPEN(iufilgap, file=name1, form='formatted', err=100, iostat=ios)
100 CALL errore('eliashberg_read_aniso_iaxis','opening file '//name1,abs(ios))
    READ(iufilgap,'(a)') word
    DO iw = 1, nsiw(itemp) ! loop over omega
       DO ik = 1, nkfs
          DO ibnd = 1, nbndfs
             IF ( abs( ekfs(ibnd,ik) - ef0 ) .lt. fsthick ) THEN
                READ(iufilgap,'(5ES20.10)') omega, eband, AZnormi(ibnd,ik,iw), ADeltai(ibnd,ik,iw), NAZnormi(ibnd,ik,iw)
                IF ( iw .eq. 1 ) & 
                   Agap(ibnd,ik,itemp) = ADeltai(ibnd,ik,1)
             ENDIF
          ENDDO ! ibnd
       ENDDO ! ik             
       IF ( abs(wsi(iw)-omega) .gt. eps ) &
          CALL errore('eliashberg_read_aniso_iaxis','temperature not the same with the input',1)
    ENDDO ! iw
    CLOSE(iufilgap)
    !
    DO iw = 1, nsiw(itemp) ! loop over omega
      DO ik = 1, nkfs
         DO ibnd = 1, nbndfs
            IF ( abs( ekfs(ibnd,ik) - ef0 ) .lt. fsthick ) THEN
               weight = 0.5d0 * wkfs(ik) * w0g(ibnd,ik) / dosef
               Znormi(iw) = Znormi(iw) + weight * AZnormi(ibnd,ik,iw)
               Deltai(iw) = Deltai(iw) + weight * ADeltai(ibnd,ik,iw)
               NZnormi(iw) = NZnormi(iw) + weight * NAZnormi(ibnd,ik,iw)
            ENDIF
         ENDDO ! ibnd
      ENDDO ! ik
    ENDDO ! iw
    gap(itemp) = Deltai(1)
    gap0 = gap(itemp)
    !
    CALL gap_FS( itemp )
    !
    IF ( iverbosity .eq. 2 ) &
       CALL free_energy( itemp )
    !
  ENDIF
  CALL mp_bcast( Deltai, ionode_id, inter_pool_comm )
  CALL mp_bcast( Znormi, ionode_id, inter_pool_comm )
  CALL mp_bcast( NZnormi, ionode_id, inter_pool_comm )
  CALL mp_bcast( ADeltai, ionode_id, inter_pool_comm )
  CALL mp_bcast( AZnormi, ionode_id, inter_pool_comm )
  CALL mp_bcast( NAZnormi, ionode_id, inter_pool_comm )
  CALL mp_bcast( gap0, ionode_id, inter_pool_comm )
  CALL mp_bcast( gap, ionode_id, inter_pool_comm )
  CALL mp_bcast( Agap, ionode_id, inter_pool_comm )
  CALL mp_barrier(inter_pool_comm)
  !
  RETURN
  !
  END SUBROUTINE eliashberg_read_aniso_iaxis
  !
  !-----------------------------------------------------------------------
