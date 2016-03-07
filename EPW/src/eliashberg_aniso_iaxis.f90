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
#include "f_defs.h"
  !
  USE kinds,         ONLY : DP
  USE io_global,     ONLY : stdout
  USE io_files,      ONLY : prefix
  USE control_flags, ONLY : iverbosity
  USE epwcom,        ONLY : nsiter, nstemp, broyden_beta, broyden_ndim, & 
                            limag, lpade, lacon, fsthick, imag_read, pade_read
  USE eliashbergcom, ONLY : nsw, nsiw, ADelta, ADeltap, ADeltai, ADeltaip, Agap, gap, &
                            Delta, Deltai, estemp, nkfs, nbndfs, ekfs, ef0
  USE constants_epw, ONLY : kelvin2eV, ci
#ifdef __PARA
  USE io_global,     ONLY : ionode_id
  USE mp_global,     ONLY : inter_pool_comm, my_pool_id, npool
  USE mp,            ONLY : mp_bcast, mp_barrier, mp_sum
  USE mp_world,      ONLY : mpime
#endif
  ! 
  IMPLICIT NONE
  !
  INTEGER :: itemp, iter, N, ik, ibnd, imelt
  REAL(DP) :: dFE, tcpu, rdeltaout(nsw), rdeltain(nsw), cdeltaout(nsw), cdeltain(nsw)
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
    WRITE(stdout,'(a)') '  '
    CALL start_clock( 'iaxis_imag' )
    CALL gen_freqgrid_iaxis( itemp )
    !
    IF ( ( limag .AND. .not. imag_read ) .OR. ( limag .AND. imag_read .AND. itemp .ne. 1 ) ) THEN
      iter = 1
      conv = .false.
      DO WHILE ( .not. conv .AND. iter .le. nsiter ) 
        CALL sum_eliashberg_aniso_iaxis( itemp, iter, conv )
#ifdef __PARA 
        IF (mpime .eq. ionode_id) THEN
#endif  
          DO ik = 1, nkfs
            DO ibnd = 1, nbndfs
              IF ( abs( ekfs(ibnd,ik) - ef0 ) .lt. fsthick ) THEN
                CALL mix_broyden_aniso( ik, ibnd, nsiw(itemp), & 
                     ADeltai(ibnd,ik,:), ADeltaip(ibnd,ik,:), broyden_beta, iter, broyden_ndim, conv )
              ENDIF
            ENDDO
          ENDDO
#ifdef __PARA
        ENDIF
        CALL mp_bcast( ADeltai, ionode_id, inter_pool_comm )
        CALL mp_bcast( ADeltaip, ionode_id, inter_pool_comm )
        CALL mp_barrier(inter_pool_comm)
#endif
        iter = iter + 1
      ENDDO ! iter
      !
      IF ( conv ) THEN
        IF ( ALLOCATED(ADeltaip) ) DEALLOCATE(ADeltaip)
          !
          ! SP : Only print the Free energy if the user want it
          !
          IF ( iverbosity .eq. 0 ) THEN
#ifdef __PARA 
            IF (mpime .eq. ionode_id) THEN
#endif 
            CALL free_energy( itemp, dFE )
#ifdef __PARA
            ENDIF
            CALL mp_barrier(inter_pool_comm)
#endif
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
        CALL read_eliashberg_aniso_iaxis( itemp )
      ENDIF
      !
      If ( ( lpade .AND. .not. pade_read ) .OR. ( lpade .AND. pade_read .AND. itemp .ne. 1 ) ) THEN 
        WRITE(stdout,'(a)') '  '  
        WRITE(stdout,'(5x,a)') 'Pade approximant of anisotropic Eliashberg equations from imaginary-axis to real-axis'
        WRITE(stdout,'(a)') '  '
        CALL start_clock( 'raxis_pade' )
        conv = .false.
        N = 90 * nsiw(itemp) / 100
        IF ( mod(N,2) .ne. 0 ) N = N + 1
        CALL pade_cont_aniso_iaxis_to_raxis( itemp, N, conv )
        !
        IF ( conv ) THEN
#ifdef __PARA 
          IF (mpime .eq. ionode_id) THEN
#endif
            CALL dos_quasiparticle( itemp )
#ifdef __PARA
          ENDIF
          CALL mp_barrier(inter_pool_comm)
#endif     
          WRITE(stdout,'(a)') '  '
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
      ELSEIF ( lpade .AND. pade_read .AND. itemp .eq. 1 ) THEN
        WRITE(stdout,'(5x,a)') 'Read from file Delta and Znorm on real-axis '
        CALL read_eliashberg_aniso_raxis_pade( itemp )
      ENDIF
      !
      IF ( lacon ) THEN 
        WRITE(stdout,'(a)') '  '
        WRITE(stdout,'(5x,a)') 'Analytic continuation of anisotropic Eliashberg equations from imaginary-axis to real-axis'
        WRITE(stdout,'(a)') '  '
        WRITE(stdout,'(5x,a,i6)') 'Total number of frequency points nsw = ', nsw
        WRITE(stdout,'(a)') '    '
        CALL start_clock( 'raxis_acon' )
        !
        iter = 1
        conv = .false.
        DO WHILE ( .not. conv .AND. iter .le. nsiter )
          CALL analytic_cont_aniso_iaxis_to_raxis( itemp, iter, conv )
#ifdef __PARA 
          IF (mpime .eq. ionode_id) THEN
#endif
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
#ifdef __PARA
          ENDIF
          CALL mp_bcast( ADelta, ionode_id, inter_pool_comm )
          CALL mp_bcast( ADeltap, ionode_id, inter_pool_comm )
          CALL mp_barrier(inter_pool_comm)
#endif
          iter = iter + 1
        ENDDO ! iter
        !
        IF ( conv ) THEN 
#ifdef __PARA 
          IF (mpime .eq. ionode_id) THEN
#endif
            CALL dos_quasiparticle( itemp )
#ifdef __PARA
        ENDIF
        CALL mp_barrier(inter_pool_comm)
#endif
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
    ELSE
      ! remove memory allocated for ws
      imelt = nsw 
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
#include "f_defs.h"
  !
  USE kinds,         ONLY : DP
  USE io_global,     ONLY : stdout
  USE io_epw,        ONLY : iufilgap, iufilgapFS
  USE io_files,      ONLY : prefix
  USE elph2,         ONLY : wqf
  USE cell_base,     ONLY : bg
  USE control_flags, ONLY : iverbosity
  USE epwcom,        ONLY : nsiter, nstemp, muc, conv_thr_iaxis, fsthick, nkf1, &
                            nkf2, nkf3
  USE eliashbergcom, ONLY : nsiw, estemp, gap0, gap, Agap, wsi, AKeri, limag_fly, & 
                            NAZnormi, AZnormi, ADeltai, ADeltaip, NZnormi, Znormi, & 
                            Deltai, wsphmax, nkfs, nbndfs, dosef, ef0, ixkqf, ixqfs, & 
                            nqfs, wkfs, w0g, ekfs, ixkff
  USE constants_epw, ONLY : pi, kelvin2eV 
#ifdef __PARA
  USE io_global,     ONLY : ionode_id
  USE mp_global,     ONLY : inter_pool_comm, my_pool_id, npool
  USE mp_world,      ONLY : mpime
  USE mp,            ONLY : mp_bcast, mp_barrier, mp_sum
#endif
  ! 
  IMPLICIT NONE
  !
  INTEGER  :: i, j, k, iw, iwp, itemp, iter, ik, iq, iq0, ibnd, jbnd, & 
             lower_bnd, upper_bnd, ibin, nbin, imelt
  REAL(DP) :: esqrt, absdelta, reldelta, errdelta, weight, temp, delta_max, dbin, sigma
  REAL(DP) :: kernelp, kernelm, lambdap, lambdam, x1, x2, x3
  REAL(DP), ALLOCATABLE :: wesqrt(:,:,:), desqrt(:,:,:), delta_k_bin(:), Agap_tmp(:,:)
  REAL(DP), ALLOCATABLE, SAVE :: Deltaold(:)
  REAL(DP), EXTERNAL :: w0gauss
  LOGICAL  :: conv
  CHARACTER (len=256) :: name1, name2
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
     DO ik = 1, nkfs
        DO ibnd = 1, nbndfs
           IF ( abs( ekfs(ibnd,ik) - ef0 ) .lt. fsthick ) THEN
              IF ( itemp .eq. 1 ) THEN
                 DO iw = 1, nsiw(itemp)
                    IF ( wsi(iw) .lt. 2.d0*wsphmax ) THEN 
                       ADeltaip(ibnd,ik,iw) = gap0 
                    ELSE
                       ADeltaip(ibnd,ik,iw) = 0.d0
                    ENDIF
                 ENDDO
              ELSEIF ( itemp .ne. 1 ) THEN
                 DO iw = 1, nsiw(itemp)
                    IF ( wsi(iw) .lt. 2.d0*wsphmax ) THEN
                       ADeltaip(ibnd,ik,iw) = gap(itemp-1)
                    ELSE
                       ADeltaip(ibnd,ik,iw) = 0.d0
                    ENDIF
                 ENDDO
              ENDIF
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
#ifdef __PARA
  ! collect contributions from all pools 
  CALL mp_sum( AZnormi, inter_pool_comm )
  CALL mp_sum( NAZnormi, inter_pool_comm )
  CALL mp_sum( ADeltai, inter_pool_comm )
  CALL mp_barrier(inter_pool_comm)
#endif
  !
#ifdef __PARA 
  IF (mpime .eq. ionode_id) THEN
#endif 
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
  WRITE(stdout,'(5x,a,i6,a,d18.9,a,d18.9,a,d18.9,a,d18.9)') 'iter = ', iter, '   relerr = ', errdelta, &
                                      '   abserr = ', reldelta / dble(nsiw(itemp)), &
                                      '   Znormi(1) = ', Znormi(1), '   Deltai(1) = ', Deltai(1)
  !
  IF ( errdelta .lt. conv_thr_iaxis) conv = .true.
  IF ( errdelta .lt. conv_thr_iaxis .OR. iter .eq. nsiter ) THEN
     temp = estemp(itemp) / kelvin2eV
     IF ( temp .lt. 10.d0 ) THEN
        WRITE(name1,'(a,a13,f4.2)') TRIM(prefix),'.imag_aniso_0', temp
     ELSEIF ( temp .ge. 10.d0 ) THEN
        WRITE(name1,'(a,a12,f5.2)') TRIM(prefix),'.imag_aniso_', temp
     ENDIF
     OPEN(iufilgap, file=name1, form='formatted')
     WRITE(iufilgap,'(5a20)') '#        w [eV]', 'Enk-Ef [eV]', 'Znorm(w) [eV]', 'Delta(w) [eV]', 'NZnorm(w) [eV]'
     DO iw = 1, nsiw(itemp) ! loop over omega
        DO ik = 1, nkfs
           DO ibnd = 1, nbndfs
              IF ( abs( ekfs(ibnd,ik) - ef0 ) .lt. fsthick ) THEN
                 WRITE(iufilgap,'(5ES20.10)') wsi(iw), ekfs(ibnd,ik)-ef0,&
                       AZnormi(ibnd,ik,iw), ADeltai(ibnd,ik,iw), NAZnormi(ibnd,ik,iw)
                 IF ( iw .eq. 1 ) Agap(ibnd,ik,itemp) = ADeltai(ibnd,ik,iw)
              ENDIF
           ENDDO ! ibnd                   
        ENDDO ! ik
     ENDDO ! iw
     CLOSE(iufilgap)
     gap(itemp) = Deltai(1)
     !
     delta_max = 1.25d0 * maxval(Agap(:,:,itemp)) 
     nbin = int(delta_max/(0.005d0/1000.d0))
     dbin = delta_max / dble(nbin)
     IF ( .not. ALLOCATED(delta_k_bin) ) ALLOCATE( delta_k_bin(nbin) )
     delta_k_bin(:) = 0.d0
     !
     DO ik = 1, nkfs
        DO ibnd = 1, nbndfs
           IF ( abs( ekfs(ibnd,ik) - ef0 ) .lt. fsthick ) THEN
              DO ibin = 1, nbin
                 sigma = 1.d0 * dbin
                 weight = w0gauss( ( Agap(ibnd,ik,itemp) - dble(ibin) * dbin) / sigma, 0 ) / sigma
                 delta_k_bin(ibin) = delta_k_bin(ibin) + weight
              ENDDO
           ENDIF
        ENDDO
     ENDDO
     !
     IF ( temp .lt. 10.d0 ) THEN
        WRITE(name1,'(a,a18,f4.2)') TRIM(prefix),'.imag_aniso_gap0_0', temp
     ELSEIF ( temp .ge. 10.d0 ) THEN
        WRITE(name1,'(a,a17,f5.2)') TRIM(prefix),'.imag_aniso_gap0_', temp
     ENDIF
     OPEN(iufilgap, file=name1, form='formatted')
     DO ibin = 1, nbin
        WRITE(iufilgap,'(2ES20.10)') temp + delta_k_bin(ibin)/maxval(delta_k_bin(:)), dbin*dble(ibin)
     ENDDO
     CLOSE(iufilgap)
     !
     IF ( ALLOCATED(delta_k_bin) ) DEALLOCATE(delta_k_bin) 
     !
     ! RM - If the k-point is outside the Fermi shell, 
     ! ixkff(ik)=0 and Agap_tmp(:,0) = 0.0
     !
     IF ( .not. ALLOCATED(Agap_tmp) ) ALLOCATE(Agap_tmp(nbndfs,0:nkfs))
     Agap_tmp(:,1:nkfs) = Agap(:,1:nkfs,itemp)
     Agap_tmp(:,0) = 0.0d0
     !
     ! SP & RM: .cube file for VESTA plotting (only if iverbosity = 2)
     !
     IF ( iverbosity .eq. 2 ) THEN
       ! 
       DO ibnd = 1, nbndfs
          IF ( temp .lt. 10.d0 ) THEN
             WRITE(name1,'(a,a18,f4.2,a1,i1,a5)')TRIM(prefix),'_imag_aniso_gap0_0', temp, '_', ibnd, '.cube'
          ELSEIF ( temp .ge. 10.d0 ) THEN
             WRITE(name1,'(a,a17,f5.2,a1,i1,a5)')TRIM(prefix),'.imag_aniso_gap0_', temp, '_', ibnd, '.cube'
          ENDIF
          OPEN(iufilgap, file=name1, form='formatted')
          WRITE(iufilgap,*) 'Cubfile created from EPW calculation'
          WRITE(iufilgap,*) 'gap'
          WRITE(iufilgap,'(i5,3f12.6)') 1, 0.0d0, 0.0d0, 0.0d0
          WRITE(iufilgap,'(i5,3f12.6)') nkf1, (bg(i,1)/dble(nkf1),i=1,3)
          WRITE(iufilgap,'(i5,3f12.6)') nkf2, (bg(i,2)/dble(nkf2),i=1,3)
          WRITE(iufilgap,'(i5,3f12.6)') nkf3, (bg(i,3)/dble(nkf3),i=1,3)
          WRITE(iufilgap,'(i5,4f12.6)') 1, 1.0d0, 0.0d0, 0.0d0, 0.0d0
          WRITE(iufilgap,'(6f12.6)') ( Agap_tmp(ibnd,ixkff(ik)),ik=1,nkf1*nkf2*nkf3 )
          CLOSE(iufilgap)
       ENDDO
       ! 
     ENDIF
     ! 
     ! SP & RM : Write on file the superconducting gap close to the Fermi surface along with 
     !     Cartesian coordinate, band index, energy distance from Fermi level and gap value.
     ! 
     IF ( temp .lt. 10.d0 ) THEN
        WRITE(name2,'(a,a20,f4.2)') TRIM(prefix),'.imag_aniso_gap_FS_0', temp
     ELSEIF ( temp .ge. 10.d0 ) THEN
        WRITE(name2,'(a,a19,f5.2)') TRIM(prefix),'.imag_aniso_gap_FS_', temp
     ENDIF
     OPEN(iufilgapFS, file=name2, form='formatted')
     WRITE(iufilgapFS,'(a78)') '#               k-point                  Band Enk-Ef [eV]        Delta(0) [eV]'
     DO i = 1, nkf1
       DO j = 1, nkf2
         DO k = 1, nkf3
           ik = k + (j-1)*nkf3 + (i-1)*nkf2*nkf3
           IF ( ixkff(ik) .gt. 0 ) THEN
             DO ibnd = 1, nbndfs
               ! RM: Everything is in eV here. 
               ! SP: Here take a 0.2 eV interval around the FS. 
               !IF ( abs( ekfs(ibnd,ixkff(ik)) - ef0 ) .lt. fsthick ) THEN
               IF ( abs( ekfs(ibnd,ixkff(ik)) - ef0 ) .lt. 0.2 ) THEN
                 x1 = bg(1,1)*(i-1)/nkf1+bg(1,2)*(j-1)/nkf2+bg(1,3)*(k-1)/nkf3
                 x2 = bg(2,1)*(i-1)/nkf1+bg(2,2)*(j-1)/nkf2+bg(2,3)*(k-1)/nkf3
                 x3 = bg(3,1)*(i-1)/nkf1+bg(3,2)*(j-1)/nkf2+bg(3,3)*(k-1)/nkf3
                 WRITE(iufilgapFS,'(3f12.6,i8,f12.6,f24.15)') x1, x2, x3, ibnd, & 
                                  ekfs(ibnd,ixkff(ik))-ef0, Agap_tmp(ibnd,ixkff(ik))
               ENDIF
             ENDDO ! ibnd    
           ENDIF  
         ENDDO  ! k           
       ENDDO ! j
     ENDDO ! i
     CLOSE(iufilgapFS)
     !
     ! isotropic case
     ! SP: Only write isotropic if user really want that
     IF ( iverbosity .eq. 2 ) THEN
       IF ( temp .lt. 10.d0 ) THEN
          WRITE(name1,'(a,a11,f4.2)') TRIM(prefix),'.imag_iso_0', temp
       ELSEIF ( temp .ge. 10.d0 ) THEN
          WRITE(name1,'(a,a10,f5.2)') TRIM(prefix),'.imag_iso_', temp
       ENDIF
       OPEN(iufilgap, file=name1, form='formatted')
       WRITE(iufilgap,'(4a24)') 'w', 'Znorm(w)', 'Delta(w)', 'NZnorm(w)'
       DO iw = 1, nsiw(itemp) ! loop over omega
          WRITE(iufilgap,'(4ES20.10)') wsi(iw), Znormi(iw), Deltai(iw), NZnormi(iw)
       ENDDO
       CLOSE(iufilgap)
     ENDIF 
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
#ifdef __PARA
  ENDIF
  CALL mp_bcast( Deltai, ionode_id, inter_pool_comm )
  CALL mp_bcast( Znormi, ionode_id, inter_pool_comm )
  CALL mp_bcast( NZnormi, ionode_id, inter_pool_comm )
  CALL mp_bcast( AZnormi, ionode_id, inter_pool_comm )
  CALL mp_bcast( NAZnormi, ionode_id, inter_pool_comm )
  CALL mp_bcast( gap, ionode_id, inter_pool_comm )
  CALL mp_bcast( Agap, ionode_id, inter_pool_comm )
  CALL mp_bcast( conv, ionode_id, inter_pool_comm )
  CALL mp_barrier(inter_pool_comm)
#endif
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
  ! This routine does the analytic continuation of the isotropic Eliashberg equations 
  ! from the imaginary-axis to the real axis
  ! reference F. Marsiglio, M. Schossmann, and J. Carbotte, Phys. Rev. B 37, 4965 (1988)
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
#include "f_defs.h"
  !
  USE kinds,         ONLY : DP
  USE io_global,     ONLY : stdout
  USE io_epw,        ONLY : iufilgap
  USE io_files,      ONLY : prefix
  USE control_flags, ONLY : iverbosity
  USE phcom,         ONLY : nmodes
  USE elph2,         ONLY : wqf, wf
  USE epwcom,        ONLY : nqstep, degaussq, nsiter, conv_thr_racon, fsthick, & 
                            lpade, eps_acustic
  USE eliashbergcom, ONLY : nsw, estemp, dwsph, ws, wsph, gap, Agap, Gp, Gm, ADsumi, AZsumi, &                           
                            Delta, Znorm, ADelta, ADeltap, AZnorm, AZnormp, g2, lacon_fly, & 
                            a2fij, wkfs, dosef, ixkqf, ixqfs, nqfs, w0g, nkfs, nbndfs, ef0, ekfs
  USE constants_epw, ONLY : pi, kelvin2eV, ci
#ifdef __PARA
  USE io_global, ONLY : ionode_id
  USE mp_global, ONLY : inter_pool_comm, my_pool_id, npool
  USE mp_world,  ONLY : mpime
  USE mp,        ONLY : mp_bcast, mp_barrier, mp_sum
#endif
  ! 
  IMPLICIT NONE
  !
  INTEGER :: i, iw, iwp, iwph, itemp, iter, ik, iq, iq0, ibnd, jbnd, imode, & 
             lower_bnd, upper_bnd, ibin, nbin, imelt
  REAL(DP) :: rgammap, rgammam, absdelta, reldelta, errdelta, weight, temp, delta_max, dbin, & 
              a2f_, sigma
  REAL(DP), ALLOCATABLE :: delta_k_bin(:)
  REAL(DP), EXTERNAL :: w0gauss
  COMPLEX(DP) :: esqrt, root
  COMPLEX(DP), ALLOCATABLE, SAVE :: Deltaold(:)
  LOGICAL :: conv, lgap
  CHARACTER (len=256) :: name1
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
     ADeltap(:,:,:) = (0.d0, 0.d0)
     AZnormp(:,:,:) = (1.d0, 0.d0)
     Deltaold(:) = (0.d0, 0.d0)
     IF ( lpade ) THEN 
        ADeltap(:,:,:) = ADelta(:,:,:)
        Deltaold(:) = Delta(:)
     ELSE
        DO ik = 1, nkfs
           DO ibnd = 1, nbndfs
              IF ( abs( ekfs(ibnd,ik) - ef0 ) .lt. fsthick ) THEN
                 ADeltap(ibnd,ik,:) = Agap(ibnd,ik,itemp) ! gap(itemp)
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
  Delta(:) = (0.d0, 0.d0)
  Znorm(:) = (0.d0, 0.d0)
  ADelta(:,:,:) = (0.d0, 0.d0)
  AZnorm(:,:,:) = (0.d0, 0.d0)
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
                             IF ( aimag(root) .lt. 0.d0 ) THEN 
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
                          IF ( aimag(root) .lt. 0.d0 ) THEN 
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
#ifdef __PARA
  ! collect contributions from all pools 
  CALL mp_sum( AZnorm, inter_pool_comm )
  CALL mp_sum( ADelta, inter_pool_comm )
  CALL mp_barrier(inter_pool_comm)
#endif
  !
#ifdef __PARA 
  IF (mpime .eq. ionode_id) THEN
#endif 
  absdelta = 0.d0
  reldelta = 0.d0
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
  WRITE(stdout,'(5x,a,i6,a,d18.9,a,d18.9,a,d18.9)') 'iter = ', iter, '   error = ', errdelta, &
                         '   Re[Znorm(1)] = ', real(Znorm(1)), '   Re[Delta(1)] = ', real(Delta(1))
  !
  IF ( errdelta .lt. conv_thr_racon ) conv = .true.
  IF ( errdelta .lt. conv_thr_racon .OR. iter .eq. nsiter ) THEN
     temp = estemp(itemp) / kelvin2eV
     IF ( iverbosity .eq. 2 ) THEN 
        IF ( temp .lt. 10.d0 ) THEN
           WRITE(name1,'(a,a13,f4.2)') TRIM(prefix),'.acon_aniso_0', temp
        ELSEIF ( temp .ge. 10.d0 ) THEN
           WRITE(name1,'(a,a12,f5.2)') TRIM(prefix),'.acon_aniso_', temp
        ENDIF
        OPEN(iufilgap, file=name1, form='formatted')
        WRITE(iufilgap,'(6a24)') 'w', 'Enk-Ef', 'Re[Znorm(w)]', 'Im[Znorm(w)]', 'Re[Delta(w)]', 'Im[Delta(w)]'
     ENDIF
     
     DO ik = 1, nkfs
        DO ibnd = 1, nbndfs
           IF ( abs( ekfs(ibnd,ik) - ef0 ) .lt. fsthick ) THEN
              lgap = .true.
                 ! DO iw = 1, nsw
                 DO iw = 1, nsw-1 ! FG: this change is prevent segfault in ADelta(*,*,iw+1) and  ws(iw+1)
                 IF ( lgap .AND. iw .lt. nqstep .AND. real(ADelta(ibnd,ik,iw)) .gt. 0.d0 & 
                      .AND. real(ADelta(ibnd,ik,iw+1)) .gt. 0.d0 & 
                      .AND. ( ws(iw) - real(ADelta(ibnd,ik,iw)) )*( ws(iw+1) - real(ADelta(ibnd,ik,iw+1)) ) .lt. 0.d0 ) THEN 
                    Agap(ibnd,ik,itemp) = (   ( real(ADelta(ibnd,ik,iw))   - ws(iw)   ) * ws(iw+1) & 
                                            - ( real(ADelta(ibnd,ik,iw+1)) - ws(iw+1) ) * ws(iw) ) &
                                        / ( ( real(ADelta(ibnd,ik,iw)) - ws(iw) ) - ( real(ADelta(ibnd,ik,iw+1)) - ws(iw+1) ) )  
                    !WRITE(stdout,'(5x,a,i6,4ES20.10)') 'iw = ', iw, ws(iw), real(ADelta(ibnd,ik,iw)), & 
                    !                                    aimag(ADelta(ibnd,ik,iw)), Agap(ibnd,ik,itemp)
                    !WRITE(stdout,'(5x,a,i6,4ES20.10)') 'iw = ', iw+1, ws(iw+1), real(ADelta(ibnd,ik,iw+1)), & 
                    !                                    aimag(ADelta(ibnd,ik,iw+1)), Agap(ibnd,ik,itemp)
                    !WRITE(stdout,'(a)') '   '
                    lgap = .false. 
                 ENDIF
                 IF ( iverbosity .eq. 2 ) THEN
                    WRITE(iufilgap,'(6ES20.10)') ws(iw), ekfs(ibnd,ik)-ef0, &
                                                real(AZnorm(ibnd,ik,iw)), aimag(AZnorm(ibnd,ik,iw)), &
                                                real(ADelta(ibnd,ik,iw)), aimag(ADelta(ibnd,ik,iw))
                 ENDIF
              ENDDO ! iw
           ENDIF
        ENDDO ! ibnd                   
     ENDDO ! ik
     IF ( iverbosity .eq. 2 ) CLOSE(iufilgap)
     !
     delta_max = 1.25d0 * maxval(Agap(:,:,itemp))
     nbin = int(delta_max/(0.005d0/1000.d0))
     dbin = delta_max / dble(nbin)
     IF ( .not. ALLOCATED(delta_k_bin) ) ALLOCATE ( delta_k_bin(nbin) )
     delta_k_bin(:) = 0.d0 
     !
     DO ik = 1, nkfs
        DO ibnd = 1, nbndfs
           IF ( abs( ekfs(ibnd,ik) - ef0 ) .lt. fsthick ) THEN
              DO ibin = 1, nbin
                 sigma = 1.d0 * dbin
                 weight = w0gauss( ( Agap(ibnd,ik,itemp) - dble(ibin) * dbin) / sigma, 0 ) / sigma
                 delta_k_bin(ibin) = delta_k_bin(ibin) + weight
              ENDDO
           ENDIF
        ENDDO
     ENDDO
     ! 
     IF ( temp .lt. 10.d0 ) THEN
        WRITE(name1,'(a,a18,f4.2)') TRIM(prefix),'.acon_aniso_gap0_0', temp
     ELSEIF ( temp .ge. 10.d0 ) THEN
        WRITE(name1,'(a,a17,f5.2)') TRIM(prefix),'.acon_aniso_gap0_', temp
     ENDIF
     OPEN(iufilgap, file=name1, form='formatted') 
     DO ibin = 1, nbin
        WRITE(iufilgap,'(2ES20.10)') temp + delta_k_bin(ibin)/maxval(delta_k_bin(:)), dbin*dble(ibin)
     ENDDO
     CLOSE(iufilgap)
     !
     IF ( ALLOCATED(delta_k_bin) ) DEALLOCATE(delta_k_bin)
     !
     ! isotropic case
     IF ( temp .lt. 10.d0 ) THEN
        WRITE(name1,'(a,a11,f4.2)') TRIM(prefix),'.acon_iso_0', temp
     ELSEIF ( temp .ge. 10.d0 ) THEN
        WRITE(name1,'(a,a10,f5.2)') TRIM(prefix),'.acon_iso_', temp
     ENDIF
     OPEN(iufilgap, file=name1, form='formatted')
     WRITE(iufilgap,'(5a18)') 'w', 'Re[Znorm(w)]', 'Im[Znorm(w)]', 'Re[Delta(w)]', 'Im[Delta(w)]'
     lgap = .true.
     ! DO iw = 1, nsw
     DO iw = 1, nsw-1   ! FG: this change is prevent segfault in Delta(iw+1) and ws(iw+1)
        IF ( lgap .AND. iw .lt. nqstep .AND. real(Delta(iw)) .gt. 0.d0 .AND. real(Delta(iw+1)) .gt. 0.d0 .AND. &
             ( ws(iw) - real(Delta(iw)) )*( ws(iw+1) - real(Delta(iw+1)) ) .lt. 0.d0 ) THEN
           gap(itemp) = ( ( real(Delta(iw)) - ws(iw) ) * ws(iw+1) - ( real(Delta(iw+1)) - ws(iw+1) ) * ws(iw) ) &
                      / ( ( real(Delta(iw)) - ws(iw) ) - ( real(Delta(iw+1)) - ws(iw+1) ) )
           !WRITE(stdout,'(a)') '   '
           !WRITE(stdout,'(5x,a)') 'gap_edge'
           !WRITE(stdout,'(5x,a,i6,4ES20.10)') 'iw = ', iw, ws(iw), real(Delta(iw)), aimag(Delta(iw)), &
           !                                 gap(itemp)
           !WRITE(stdout,'(5x,a,i6,4ES20.10)') 'iw = ', iw+1, ws(iw+1), real(Delta(iw+1)), aimag(Delta(iw+1)), &
           !                                 gap(itemp)
           !WRITE(stdout,'(a)') '   '
           lgap = .false. 
        ENDIF
        WRITE(iufilgap,'(5ES20.10)') ws(iw), real(Znorm(iw)), aimag(Znorm(iw)), &
                                    real(Delta(iw)), aimag(Delta(iw))
     ENDDO ! iw
     CLOSE(iufilgap)
  ENDIF
  !
  IF ( conv .OR. iter .eq. nsiter ) THEN
     WRITE(stdout,'(5x,a,i6)') 'Convergence was reached in nsiter = ', iter
  ENDIF
  IF ( .not. conv .AND. iter .eq. nsiter ) THEN
     WRITE(stdout,'(5x,a,i6)') 'Convergence was not reached in nsiter = ', iter
     CALL errore('analytic_cont_aniso_iaxis_to_raxis','increase nsiter or reduce conv_thr_racon',1)
  ENDIF
#ifdef __PARA
  ENDIF
  CALL mp_bcast( Delta, ionode_id, inter_pool_comm )
  CALL mp_bcast( Znorm, ionode_id, inter_pool_comm )
  CALL mp_bcast( AZnorm, ionode_id, inter_pool_comm )
  CALL mp_bcast( gap, ionode_id, inter_pool_comm )
  CALL mp_bcast( Agap, ionode_id, inter_pool_comm )
  CALL mp_bcast( conv, ionode_id, inter_pool_comm )
  CALL mp_barrier(inter_pool_comm)
#endif
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
  ! input
  !
  ! itemp  - temperature point
  ! N      - number Pade approximants
  ! conv   - convergence flag 
  !
  ! output 
  !
  ! conv   - convergence flag 
  !
#include "f_defs.h"
  !
  USE kinds,         ONLY : DP
  USE io_global,     ONLY : stdout
  USE io_epw,        ONLY : iufilgap
  USE io_files,      ONLY : prefix
  USE control_flags, ONLY : iverbosity
  USE epwcom,        ONLY : nqstep, fsthick
  USE eliashbergcom, ONLY : nsw, estemp, ws, wsi, gap, Agap, &                                            
                            Delta, Znorm, ADelta, AZnorm, ADeltai, AZnormi, &              
                            wkfs, dosef, w0g, nkfs, nbndfs, ef0, ekfs
  USE constants_epw, ONLY : kelvin2eV, cone, ci
#ifdef __PARA
  USE io_global, ONLY : ionode_id
  USE mp_global, ONLY : inter_pool_comm, my_pool_id, npool
  USE mp_world,  ONLY : mpime
  USE mp,        ONLY : mp_bcast, mp_barrier, mp_sum
#endif
  ! 
  IMPLICIT NONE
  !
  INTEGER :: iw, itemp, N, ik, ibnd, lower_bnd, upper_bnd, ibin, nbin, imelt
  REAL(DP) :: absdelta, reldelta, errdelta, weight, temp, delta_max, dbin, sigma
  REAL(DP), ALLOCATABLE :: delta_k_bin(:)
  REAL(DP), EXTERNAL :: w0gauss
  COMPLEX(DP) :: omega, padapp
  COMPLEX(DP), ALLOCATABLE :: a(:), b(:), z(:), u(:), v(:), Deltaold(:)
  LOGICAL :: conv, lgap
  CHARACTER(len=256) :: name1
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
  Delta(:) = (0.d0, 0.d0)
  Znorm(:) = (0.d0, 0.d0)
  ADelta(:,:,:) = (0.d0, 0.d0)
  AZnorm(:,:,:) = (0.d0, 0.d0)
  Deltaold(:) = gap(itemp)
  absdelta = 0.d0
  reldelta = 0.d0
  a(:) = (0.d0, 0.d0)
  b(:) = (0.d0, 0.d0)
  z(:) = (0.d0, 0.d0)
  u(:) = (0.d0, 0.d0)
  v(:) = (0.d0, 0.d0)
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
#ifdef __PARA
  ! collect contributions from all pools 
  CALL mp_sum( AZnorm, inter_pool_comm )
  CALL mp_sum( ADelta, inter_pool_comm )
  CALL mp_barrier(inter_pool_comm)
#endif
  !
#ifdef __PARA 
  IF (mpime .eq. ionode_id) THEN
#endif 
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
  IF ( errdelta .gt. 0.d0 ) THEN
     conv = .true.
     WRITE(stdout,'(5x,a,i6,a,d18.9,a,d18.9,a,d18.9)') 'pade = ', N, '   error = ', errdelta, &
                  '   Re[Znorm(1)] = ', real(Znorm(1)), '   Re[Delta(1)] = ', real(Delta(1))
     temp = estemp(itemp) / kelvin2eV
     IF ( iverbosity .eq. 2 ) THEN 
        IF ( temp .lt. 10.d0 ) THEN
           WRITE(name1,'(a,a13,f4.2)') TRIM(prefix),'.pade_aniso_0', temp
        ELSEIF ( temp .ge. 10.d0 ) THEN
           WRITE(name1,'(a,a12,f5.2)') TRIM(prefix),'.pade_aniso_', temp
        ENDIF
        OPEN(iufilgap, file=name1, form='formatted')
        WRITE(iufilgap,'(6a24)') 'w', 'Enk-Ef', 'Re[Znorm(w)]', 'Im[Znorm(w)]', 'Re[Delta(w)]', 'Im[Delta(w)]'
     ENDIF
     !
     DO ik = 1, nkfs
        DO ibnd = 1, nbndfs
           IF ( abs( ekfs(ibnd,ik) - ef0 ) .lt. fsthick ) THEN
              lgap = .true.
              ! DO iw = 1, nsw
              DO iw = 1, nsw-1   ! FG: this change is to prevent segfault in ws(iw+1) and ADelta(*,*,iw+1)
                 IF ( lgap .AND. iw .lt. nqstep .AND. real(ADelta(ibnd,ik,iw)) .gt. 0.d0 & 
                      .AND. real(ADelta(ibnd,ik,iw+1)) .gt. 0.d0 & 
                      .AND. ( ws(iw) - real(ADelta(ibnd,ik,iw)) )*( ws(iw+1) - real(ADelta(ibnd,ik,iw+1)) ) .lt. 0.d0 ) THEN 
                    Agap(ibnd,ik,itemp) = (   ( real(ADelta(ibnd,ik,iw))   - ws(iw)   ) * ws(iw+1) & 
                                            - ( real(ADelta(ibnd,ik,iw+1)) - ws(iw+1) ) * ws(iw) ) &
                                        / ( ( real(ADelta(ibnd,ik,iw)) - ws(iw) ) - ( real(ADelta(ibnd,ik,iw+1)) - ws(iw+1) ) )  
                    !WRITE(stdout,'(5x,a,i6,4ES20.10)') 'iw = ', iw, ws(iw), real(ADelta(ibnd,ik,iw)), & 
                    !                                 aimag(ADelta(ibnd,ik,iw)), Agap(ibnd,ik,itemp)
                    !WRITE(stdout,'(5x,a,i6,4ES20.10)') 'iw = ', iw+1, ws(iw+1), real(ADelta(ibnd,ik,iw+1)), & 
                    !                                 aimag(ADelta(ibnd,ik,iw+1)), Agap(ibnd,ik,itemp)
                    !WRITE(stdout,'(a)') '   '
                    lgap = .false. 
                 ENDIF
                 IF ( iverbosity .eq. 2 ) THEN 
                    WRITE(iufilgap,'(6ES20.10)') ws(iw), ekfs(ibnd,ik)-ef0, &
                                                real(AZnorm(ibnd,ik,iw)), aimag(AZnorm(ibnd,ik,iw)), &
                                                real(ADelta(ibnd,ik,iw)), aimag(ADelta(ibnd,ik,iw))
                 ENDIF
              ENDDO ! iw
           ENDIF
        ENDDO ! ibnd                   
     ENDDO ! ik
     IF ( iverbosity .eq. 2 ) CLOSE(iufilgap)
     !
     delta_max = 1.25d0 * maxval(Agap(:,:,itemp))
     nbin = int(delta_max/(0.005d0/1000.d0))
     dbin = delta_max / dble(nbin)
     IF ( .not. ALLOCATED(delta_k_bin) ) ALLOCATE ( delta_k_bin(nbin) )
     delta_k_bin(:) = 0.d0 
     !
     DO ik = 1, nkfs
        DO ibnd = 1, nbndfs
           IF ( abs( ekfs(ibnd,ik) - ef0 ) .lt. fsthick ) THEN
              DO ibin = 1, nbin
                 sigma = 1.d0 * dbin
                 weight = w0gauss( ( Agap(ibnd,ik,itemp) - dble(ibin) * dbin) / sigma, 0 ) / sigma
                 delta_k_bin(ibin) = delta_k_bin(ibin) + weight
              ENDDO
           ENDIF
        ENDDO
     ENDDO
     ! 
     IF ( temp .lt. 10.d0 ) THEN
        WRITE(name1,'(a,a18,f4.2)') TRIM(prefix),'.pade_aniso_gap0_0', temp
     ELSEIF ( temp .ge. 10.d0 ) THEN
        WRITE(name1,'(a,a17,f5.2)') TRIM(prefix),'.pade_aniso_gap0_', temp
     ENDIF
     OPEN(iufilgap, file=name1, form='formatted') 
     DO ibin = 1, nbin
        WRITE(iufilgap,'(2ES20.10)') temp + delta_k_bin(ibin)/maxval(delta_k_bin(:)), dbin*dble(ibin)
     ENDDO
     CLOSE(iufilgap)
     !
     IF ( ALLOCATED(delta_k_bin) ) DEALLOCATE(delta_k_bin)
     !
     ! isotropic case
     ! SP: Only write isotropic if user really want that
     IF ( iverbosity .eq. 2 ) THEN     
       IF ( temp .lt. 10.d0 ) THEN
          WRITE(name1,'(a,a11,f4.2)') TRIM(prefix),'.pade_iso_0', temp
       ELSEIF ( temp .ge. 10.d0 ) THEN
          WRITE(name1,'(a,a10,f5.2)') TRIM(prefix),'.pade_iso_', temp
       ENDIF
       OPEN(iufilgap, file=name1, form='formatted')
       WRITE(iufilgap,'(5a18)') 'w', 'Re[Znorm(w)]', 'Im[Znorm(w)]', 'Re[Delta(w)]', 'Im[Delta(w)]'
       lgap = .true.
       !    DO iw = 1, nsw
       DO iw = 1, nsw-1   ! this change is to prevent segfault in Delta(iw+1) and ws(iw+1)
          IF ( lgap .AND. iw .lt. nqstep .AND. real(Delta(iw)) .gt. 0.d0 .AND. real(Delta(iw+1)) .gt. 0.d0 .AND. &
               ( ws(iw) - real(Delta(iw)) )*( ws(iw+1) - real(Delta(iw+1)) ) .lt. 0.d0 ) THEN
             gap(itemp) = ( ( real(Delta(iw)) - ws(iw) ) * ws(iw+1) - ( real(Delta(iw+1)) - ws(iw+1) ) * ws(iw) ) &
                        / ( ( real(Delta(iw)) - ws(iw) ) - ( real(Delta(iw+1)) - ws(iw+1) ) )
             !WRITE(stdout,'(a)') '   '
             !WRITE(stdout,'(5x,a)') 'gap_edge'
             !WRITE(stdout,'(5x,a,i6,4d18.9)') 'iw = ', iw, ws(iw), real(Delta(iw)), aimag(Delta(iw)), &
             !                                 gap(itemp)
             !WRITE(stdout,'(5x,a,i6,4d18.9)') 'iw = ', iw+1, ws(iw+1), real(Delta(iw+1)), aimag(Delta(iw+1)), &
             !                                 gap(itemp)
             !WRITE(stdout,'(a)') '   '
             lgap = .false. 
          ENDIF
          WRITE(iufilgap,'(5ES20.10)') ws(iw), real(Znorm(iw)), aimag(Znorm(iw)), &
                                     real(Delta(iw)), aimag(Delta(iw))
       ENDDO ! iw
       CLOSE(iufilgap)
     ENDIF
  ENDIF
#ifdef __PARA
  ENDIF
  CALL mp_bcast( Delta, ionode_id, inter_pool_comm )
  CALL mp_bcast( Znorm, ionode_id, inter_pool_comm )
  CALL mp_bcast( gap, ionode_id, inter_pool_comm )
  CALL mp_bcast( Agap, ionode_id, inter_pool_comm )
  CALL mp_bcast( conv, ionode_id, inter_pool_comm )
  CALL mp_barrier(inter_pool_comm)
#endif
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
  SUBROUTINE read_eliashberg_aniso_iaxis( itemp )
  !-----------------------------------------------------------------------
  !  
  ! This routine reads from file the anisotropic Delta and Znorm on the imaginary-axis
  ! 
  ! input
  !
  ! itemp  - temperature point
  ! 
#include "f_defs.h"
  !     
  USE kinds,         ONLY : DP
  USE io_epw,        ONLY : iufilgap
  USE io_files,      ONLY : prefix
  USE epwcom,        ONLY : nstemp, fsthick
  USE eliashbergcom, ONLY : nsiw, estemp, gap, Agap, wsi, Znormi, Deltai, & 
                            AZnormi, NAZnormi, ADeltai, nkfs, nbndfs, ef0, ekfs
  USE constants_epw, ONLY : kelvin2eV
#ifdef __PARA
  USE io_global, ONLY : ionode_id
  USE mp_global, ONLY : inter_pool_comm, my_pool_id, npool
  USE mp_world,  ONLY : mpime
  USE mp,        ONLY : mp_bcast, mp_barrier, mp_sum
#endif
  ! 
  IMPLICIT NONE
  !
  INTEGER :: iw, itemp, ik, ibnd, imelt, ios
  REAL(DP) :: temp, eband, omega
  CHARACTER (len=256) :: name1, word
  !
  ! get the size of required allocated memory 
  imelt = ( 1 + nbndfs * nkfs ) * nstemp + ( 2 + 3 * nbndfs * nkfs ) * nsiw(itemp)
  CALL mem_size_eliashberg( imelt )
  !
  IF ( .not. ALLOCATED(gap) )      ALLOCATE( gap(nstemp) )
  IF ( .not. ALLOCATED(Agap) )     ALLOCATE( Agap(nbndfs,nkfs,nstemp) )
  IF ( .not. ALLOCATED(Deltai) )   ALLOCATE( Deltai(nsiw(itemp)) )
  IF ( .not. ALLOCATED(Znormi) )   ALLOCATE( Znormi(nsiw(itemp)) )
  IF ( .not. ALLOCATED(ADeltai) )  ALLOCATE( ADeltai(nbndfs,nkfs,nsiw(itemp)) )
  IF ( .not. ALLOCATED(AZnormi) )  ALLOCATE( AZnormi(nbndfs,nkfs,nsiw(itemp)) )
  IF ( .not. ALLOCATED(NAZnormi) ) ALLOCATE( NAZnormi(nbndfs,nkfs,nsiw(itemp)) )
  gap(:) = 0.d0
  Agap(:,:,:) = 0.d0
  Deltai(:) = 0.d0
  Znormi(:) = 0.d0
  ADeltai(:,:,:) = 0.d0
  AZnormi(:,:,:) = 0.d0
  NAZnormi(:,:,:) = 0.d0
  !
#ifdef __PARA                            
  IF (mpime .eq. ionode_id) THEN     
#endif 
  !   
  temp = estemp(itemp) / kelvin2eV
  ! anisotropic case
  IF ( temp .lt. 10.d0 ) THEN
     WRITE(name1,'(a,a13,f4.2)') TRIM(prefix),'.imag_aniso_0', temp
  ELSEIF ( temp .ge. 10.d0 ) THEN
     WRITE(name1,'(a,a12,f5.2)') TRIM(prefix),'.imag_aniso_', temp
  ENDIF 
  OPEN(iufilgap, file=name1, form='formatted', err=100, iostat=ios)
100 CALL errore('read_eliashberg_aniso_iaxis','opening file '//name1,abs(ios))
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
  ENDDO ! iw
  CLOSE(iufilgap)
  !
  ! isotropic case
  IF ( temp .lt. 10.d0 ) THEN
     WRITE(name1,'(a,a11,f4.2)') TRIM(prefix),'.imag_iso_0', temp
  ELSEIF ( temp .ge. 10.d0 ) THEN
     WRITE(name1,'(a,a10,f5.2)') TRIM(prefix),'.imag_iso_', temp
  ENDIF
  OPEN(iufilgap, file=name1, form='formatted', err=200, iostat=ios)
200 CALL errore('read_eliashberg_aniso_iaxis','opening file '//name1,abs(ios))
  READ(iufilgap,'(a)') word
  DO iw = 1, nsiw(itemp) ! loop over omega
     READ(iufilgap,'(3ES20.10)') omega, Znormi(iw), Deltai(iw)
     IF ( abs(wsi(iw)-omega) .gt. 10.d-8 ) &
        CALL errore('read_eliashberg_aniso_iaxis','temperature not the same with the input',1)
  ENDDO
  gap(itemp) = Deltai(1)
  CLOSE(iufilgap)
  !
#ifdef __PARA
  ENDIF
  CALL mp_bcast( Deltai, ionode_id, inter_pool_comm )
  CALL mp_bcast( Znormi, ionode_id, inter_pool_comm )
  CALL mp_bcast( ADeltai, ionode_id, inter_pool_comm )
  CALL mp_bcast( AZnormi, ionode_id, inter_pool_comm )
  CALL mp_bcast( NAZnormi, ionode_id, inter_pool_comm )
  CALL mp_bcast( gap, ionode_id, inter_pool_comm )
  CALL mp_bcast( Agap, ionode_id, inter_pool_comm )
  CALL mp_barrier(inter_pool_comm)
#endif
  !
  RETURN
  !
  END SUBROUTINE read_eliashberg_aniso_iaxis
  !
  !-----------------------------------------------------------------------
  SUBROUTINE read_eliashberg_aniso_raxis_analytic_cont( itemp )
  !-----------------------------------------------------------------------
  !  
  ! This routine reads from file the anisotropic Delta and Znorm on the real-axis
  ! 
  ! input
  !
  ! itemp  - temperature point
  ! 
#include "f_defs.h"
  !     
  USE kinds,         ONLY : DP
  USE io_epw,        ONLY : iufilgap
  USE io_files,      ONLY : prefix
  USE epwcom,        ONLY : nstemp, fsthick
  USE eliashbergcom, ONLY : nsw, estemp, gap, Agap, ws, Znorm, Delta, AZnorm, ADelta, & 
                            nkfs, nbndfs, ef0, ekfs
  USE constants_epw, ONLY : kelvin2eV, ci 
#ifdef __PARA
  USE io_global, ONLY : ionode_id
  USE mp_global, ONLY : inter_pool_comm, my_pool_id, npool
  USE mp_world,  ONLY : mpime
  USE mp,        ONLY : mp_bcast, mp_barrier, mp_sum
#endif
  ! 
  IMPLICIT NONE
  !
  INTEGER :: iw, itemp, ik, ibnd, imelt, ios
  REAL(DP) :: a, b, c, d, temp, eband, omega
  CHARACTER (len=256) :: name1, word
  !
  ! get the size of required allocated memory 
  imelt = (1 + nbndfs * nkfs ) * nstemp + 2 * ( 2 + 2 * nbndfs * nkfs ) * nsw
  CALL mem_size_eliashberg( imelt )
  !
  IF ( .not. ALLOCATED(gap) )    ALLOCATE( gap(nstemp) )
  IF ( .not. ALLOCATED(Agap) )   ALLOCATE( Agap(nbndfs,nkfs,nstemp) )
  IF ( .not. ALLOCATED(Delta) )  ALLOCATE( Delta(nsw) )
  IF ( .not. ALLOCATED(Znorm) )  ALLOCATE( Znorm(nsw) )
  IF ( .not. ALLOCATED(ADelta) ) ALLOCATE( ADelta(nbndfs,nkfs,nsw) )
  IF ( .not. ALLOCATED(AZnorm) ) ALLOCATE( AZnorm(nbndfs,nkfs,nsw) )
  gap(:) = 0.d0
  Agap(:,:,:) = 0.d0
  Delta(:) = (0.d0, 0.d0)
  Znorm(:) = (0.d0, 0.d0)
  ADelta(:,:,:) = (0.d0, 0.d0)
  AZnorm(:,:,:) = (0.d0, 0.d0)
  !
#ifdef __PARA                            
  IF (mpime .eq. ionode_id) THEN
#endif 
  !   
  temp = estemp(itemp) / kelvin2eV
  ! anisotropic case
  IF ( temp .lt. 10.d0 ) THEN
     WRITE(name1,'(a,a13,f4.2)') TRIM(prefix),'.acon_aniso_0', temp
  ELSEIF ( temp .ge. 10.d0 ) THEN
     WRITE(name1,'(a,a12,f5.2)') TRIM(prefix),'.acon_aniso_', temp
  ENDIF
  OPEN(iufilgap, file=name1, form='formatted', err=100, iostat=ios)
100 CALL errore('read_eliashberg_aniso_raxis_analytic_cont','opening file '//name1,abs(ios))
  READ(iufilgap,'(a)') word
  DO ik = 1, nkfs
     DO ibnd = 1, nbndfs
        IF ( abs( ekfs(ibnd,ik) - ef0 ) .lt. fsthick ) THEN
           DO iw = 1, nsw ! loop over omega
              READ(iufilgap,'(6ES20.10)') omega, eband, a, b, c, d 
              AZnorm(ibnd,ik,iw) = a + ci*b
              ADelta(ibnd,ik,iw) = c + ci*d
              IF ( iw .eq. 1 ) & 
                 Agap(ibnd,ik,itemp) = real(ADelta(ibnd,ik,1))
           ENDDO ! iw
        ENDIF
     ENDDO ! ibnd             
  ENDDO ! ik
  CLOSE(iufilgap)
  !
  ! isotropic case
  IF ( temp .lt. 10.d0 ) THEN
     WRITE(name1,'(a,a11,f4.2)') TRIM(prefix),'.acon_iso_0', temp
  ELSEIF ( temp .ge. 10.d0 ) THEN
     WRITE(name1,'(a,a10,f5.2)') TRIM(prefix),'.acon_iso_', temp
  ENDIF
  OPEN(iufilgap, file=name1, form='formatted', err=200, iostat=ios)
200 CALL errore('read_eliashberg_aniso_raxis_analytic_cont','opening file '//name1,abs(ios))
  READ(iufilgap,'(a)') word 
  DO iw = 1, nsw ! loop over omega
     READ(iufilgap,'(5ES20.10)') omega, a, b, c, d
           Znorm(iw) = a + ci*b
           Delta(iw) = c + ci*d
     IF ( abs(ws(iw)-omega) .gt. 10.d-8 ) &
        CALL errore('read_eliashberg_aniso_raxis_analytic_cont','temperature not the same with the input',1)
  ENDDO
  gap(itemp) = REAL(Delta(1))
  CLOSE(iufilgap)
  !
#ifdef __PARA
  ENDIF
  CALL mp_bcast( Delta, ionode_id, inter_pool_comm )
  CALL mp_bcast( Znorm, ionode_id, inter_pool_comm )
  CALL mp_bcast( ADelta, ionode_id, inter_pool_comm )
  CALL mp_bcast( AZnorm, ionode_id, inter_pool_comm )
  CALL mp_bcast( gap, ionode_id, inter_pool_comm )
  CALL mp_bcast( Agap, ionode_id, inter_pool_comm )
  CALL mp_barrier(inter_pool_comm)
#endif
  !
  RETURN
  !
  END SUBROUTINE read_eliashberg_aniso_raxis_analytic_cont
  !
  !-----------------------------------------------------------------------
  SUBROUTINE read_eliashberg_aniso_raxis_pade( itemp )
  !-----------------------------------------------------------------------
  !  
  ! This routine reads from file the anisotropic Delta and Znorm on the real-axis
  ! 
  ! input
  !
  ! itemp  - temperature point
  ! 
#include "f_defs.h"
  !     
  USE kinds,         ONLY : DP
  USE io_epw,        ONLY : iufilgap
  USE io_files,      ONLY : prefix
  USE epwcom,        ONLY : nstemp, fsthick
  USE eliashbergcom, ONLY : nsw, estemp, gap, Agap, ws, Znorm, Delta, AZnorm, ADelta, & 
                            nkfs, nbndfs, ef0, ekfs
  USE constants_epw, ONLY : kelvin2eV, ci 
#ifdef __PARA
  USE io_global, ONLY : ionode_id
  USE mp_global, ONLY : inter_pool_comm, my_pool_id, npool
  USE mp_world,  ONLY : mpime
  USE mp,        ONLY : mp_bcast, mp_barrier, mp_sum
#endif
  ! 
  IMPLICIT NONE
  !
  INTEGER :: iw, itemp, ik, ibnd, imelt, ios
  REAL(DP) :: a, b, c, d, temp, eband, omega
  CHARACTER (len=256) :: name1, word
  !
  ! get the size of required allocated memory 
  imelt = (1 + nbndfs * nkfs ) * nstemp + 2 * ( 2 + 2 * nbndfs * nkfs ) * nsw
  CALL mem_size_eliashberg( imelt )
  !
  IF ( .not. ALLOCATED(gap) )    ALLOCATE( gap(nstemp) )
  IF ( .not. ALLOCATED(Agap) )   ALLOCATE( Agap(nbndfs,nkfs,nstemp) )
  IF ( .not. ALLOCATED(Delta) )  ALLOCATE( Delta(nsw) )
  IF ( .not. ALLOCATED(Znorm) )  ALLOCATE( Znorm(nsw) )
  IF ( .not. ALLOCATED(ADelta) ) ALLOCATE( ADelta(nbndfs,nkfs,nsw) )
  IF ( .not. ALLOCATED(AZnorm) ) ALLOCATE( AZnorm(nbndfs,nkfs,nsw) )
  gap(:) = 0.d0
  Agap(:,:,:) = 0.d0
  Delta(:) = (0.d0, 0.d0)
  Znorm(:) = (0.d0, 0.d0)
  ADelta(:,:,:) = (0.d0, 0.d0)
  AZnorm(:,:,:) = (0.d0, 0.d0)
  !        
#ifdef __PARA                            
  IF (mpime .eq. ionode_id) THEN 
#endif  
  !   
  temp = estemp(itemp) / kelvin2eV
  ! anisotropic case
  IF ( temp .lt. 10.d0 ) THEN
     WRITE(name1,'(a,a13,f4.2)') TRIM(prefix),'.pade_aniso_0', temp
  ELSEIF ( temp .ge. 10.d0 ) THEN
     WRITE(name1,'(a,a12,f5.2)') TRIM(prefix),'.pade_aniso_', temp
  ENDIF
  OPEN(iufilgap, file=name1, form='formatted', err=100, iostat=ios)
100 CALL errore('read_eliashberg_aniso_raxis_pade','opening file '//name1,abs(ios))
  READ(iufilgap,'(a)') word
  DO ik = 1, nkfs
     DO ibnd = 1, nbndfs
        IF ( abs( ekfs(ibnd,ik) - ef0 ) .lt. fsthick ) THEN
           DO iw = 1, nsw ! loop over omega
              READ(iufilgap,'(6ES20.10)') omega, eband, a, b, c, d
              AZnorm(ibnd,ik,iw) = a + ci*b
              ADelta(ibnd,ik,iw) = c + ci*d
              IF ( iw .eq. 1 ) & 
                 Agap(ibnd,ik,itemp) = real(ADelta(ibnd,ik,1))
           ENDDO ! iw
        ENDIF
     ENDDO ! ibnd             
  ENDDO ! ik 
  CLOSE(iufilgap)
  !
  ! isotropic case
  IF ( temp .lt. 10.d0 ) THEN
     WRITE(name1,'(a,a11,f4.2)') TRIM(prefix),'.pade_iso_0', temp
  ELSEIF ( temp .ge. 10.d0 ) THEN
     WRITE(name1,'(a,a10,f5.2)') TRIM(prefix),'.pade_iso_', temp
  ENDIF
  OPEN(iufilgap, file=name1, form='formatted', err=200, iostat=ios)
200 CALL errore('read_eliashberg_aniso_raxis_pade','opening file '//name1,abs(ios))
  READ(iufilgap,'(a)') word
  DO iw = 1, nsw ! loop over omega
     READ(iufilgap,'(5ES20.10)') omega, a, b, c, d
           Znorm(iw) = a + ci*b
           Delta(iw) = c + ci*d
     IF ( abs(ws(iw)-omega) .gt. 10.d-8 ) &
        CALL errore('read_eliashberg_aniso_raxis_pade','temperature not the same with the input',1)
  ENDDO
  gap(itemp) = REAL(Delta(1)) 
  CLOSE(iufilgap)
  !
#ifdef __PARA
  ENDIF
  CALL mp_bcast( Delta, ionode_id, inter_pool_comm )
  CALL mp_bcast( Znorm, ionode_id, inter_pool_comm )
  CALL mp_bcast( ADelta, ionode_id, inter_pool_comm )
  CALL mp_bcast( AZnorm, ionode_id, inter_pool_comm )
  CALL mp_bcast( gap, ionode_id, inter_pool_comm )
  CALL mp_bcast( Agap, ionode_id, inter_pool_comm )
  CALL mp_barrier(inter_pool_comm)
#endif
  !
  RETURN
  !
  END SUBROUTINE read_eliashberg_aniso_raxis_pade
  !
  !-----------------------------------------------------------------------

