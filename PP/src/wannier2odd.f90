!
! Copyright (C) 2003-2013 Quantum ESPRESSO and Wannier90 groups
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
! Written by Riccardo De Gennaro, EPFL (Sept 2020).
!
!
!----------------------------------------------------------------------------
MODULE wannier2odd
  !--------------------------------------------------------------------------
  !
  USE kinds,               ONLY : DP
  !
  !
  IMPLICIT NONE
  !
  PRIVATE
  !
  PUBLIC :: wan2odd
  !
  COMPLEX(DP), ALLOCATABLE :: evcx(:,:)
  COMPLEX(DP), ALLOCATABLE :: evcx_dis(:,:)
  COMPLEX(DP), ALLOCATABLE :: evcw(:)
  COMPLEX(DP), ALLOCATABLE :: ewan(:,:)
  COMPLEX(DP), ALLOCATABLE :: psicx(:)
  COMPLEX(DP), ALLOCATABLE :: rhor(:), rhow(:)         ! real space supercell density
  COMPLEX(DP), ALLOCATABLE :: rhog(:,:), rhowg(:,:)    ! G-space supercell density
  !
  !
  CONTAINS
  !
  !--------------------------------------------------------------------------
  SUBROUTINE wan2odd( seedname, ikstart, plot, split_evc_file, gamma_trick )
    !------------------------------------------------------------------------
    !
    ! ...  This routine:
    !
    ! ...  1) reads the KS states u_nk(G) from PW and the rotation 
    ! ...     matrices U(k) from Wannier90
    !
    ! ...  2) Fourier transforms u_nk to real-space and extends them 
    ! ...     to the supercell defined by the k-points sampling
    !
    ! ...  3) applies the matrices U(k) and realizes the Wannier 
    ! ...     functions in G-space
    !
    ! ...  4) Wannier functions are finally written in a CP-readable
    ! ...     file
    !
    USE io_global,           ONLY : stdout
    USE io_files,            ONLY : nwordwfc, iunwfc, restart_dir
    USE io_base,             ONLY : write_rhog
    USE mp_pools,            ONLY : my_pool_id
    USE mp_bands,            ONLY : my_bgrp_id, root_bgrp_id, &
                                    root_bgrp, intra_bgrp_comm
    USE wavefunctions,       ONLY : evc, psic
    USE fft_base,            ONLY : dffts
    USE fft_interfaces,      ONLY : invfft, fwfft
    USE buffers,             ONLY : open_buffer, close_buffer, &
                                    save_buffer, get_buffer
    USE lsda_mod,            ONLY : nspin
    USE klist,               ONLY : xk, ngk, igk_k, nelec
    USE gvect,               ONLY : ngm
    USE wvfct,               ONLY : wg, npwx, nbnd
    USE cell_base,           ONLY : tpiba, omega, at
    USE control_flags,       ONLY : gamma_only
    USE constants,           ONLY : tpi
    USE noncollin_module,    ONLY : npol
    USE scell_wfc,           ONLY : extend_wfc
    USE fft_supercell,       ONLY : dfftcp, setup_scell_fft, bg_cp, at_cp, omega_cp, &
                                    gamma_only_x, npwxcp, ngmcp, mill_cp, ig_l2g_cp, &
                                    iunwann, nwordwann, check_fft
    USE cp_files,            ONLY : write_wannier_cp
    USE read_wannier
    !
    !
    IMPLICIT NONE
    !
    CHARACTER(LEN=256), INTENT(IN) :: seedname
    INTEGER, INTENT(IN) :: ikstart
    LOGICAL, INTENT(IN) :: plot
    LOGICAL, INTENT(IN) :: split_evc_file
    LOGICAL, INTENT(IN) :: gamma_trick
    !
    CHARACTER(LEN=256) :: dirname
    INTEGER :: ik, ikevc, ibnd, iw, ip
    INTEGER :: i, j, k, ir, n, counter
    INTEGER :: num_inc
    INTEGER :: npw
    INTEGER :: nwordwfcx                            ! record length for supercell wfcs
    INTEGER :: iunwfcx = 24                         ! unit for supercell wfc file
    INTEGER :: io_level = 1
    LOGICAL :: exst, opnd
    LOGICAL :: calc_rho=.true.
    REAL(DP) :: kvec(3), rvec(3)
    REAL(DP) :: dot_prod
    COMPLEX(DP) :: phase
    LOGICAL :: wf_is_cmplx
    REAL(DP) :: ratio
    !
    !
    CALL start_clock( 'wannier2odd' )
    !
    !
    ! ... read Wannier90 output and initialize the supercell
    !
    CALL read_wannier_chk( seedname )
    gamma_only_x = gamma_trick
    IF ( gamma_trick ) WRITE( stdout, 10 ) 
    CALL setup_scell_fft( )
    !
    !
    ! ... if all the occupied states are in, we calculate also
    ! ... the total density
    !
    IF ( ANY(excluded_band(1:nelec/2)) ) calc_rho = .false.
    !
    !
    nwordwfcx = num_bands*npwxcp*npol
    nwordwann = num_wann*npwxcp*npol
    !
    CALL dealloc_w2odd( )
    ALLOCATE( evcx(npwxcp*npol,num_bands) )
    ALLOCATE( evcw(npwxcp*npol) )
    ALLOCATE( ewan(npwxcp*npol,num_wann) )
    ALLOCATE( psicx(dfftcp%nnr) )
    !
    IF ( calc_rho ) THEN
      !
      ALLOCATE( rhor(dfftcp%nnr) )
      ALLOCATE( rhog(ngmcp,nspin) )
      rhor(:) = ( 0.D0, 0.D0 )
      !
      ALLOCATE( rhow(dfftcp%nnr) )
      ALLOCATE( rhowg(ngmcp,nspin) )
      rhow(:) = ( 0.D0, 0.D0 )
      !
    ENDIF
    !
    IF ( have_disentangled ) ALLOCATE( evcx_dis(npwxcp*npol,MAXVAL(ndimwin)) )
    !
    !
    ! ... open buffer for direct-access to the extended wavefunctions
    !
    CALL open_buffer( iunwfcx, 'wfcx', nwordwfcx, io_level, exst )
    !
    ! ... loop to read the primitive cell wavefunctions and 
    ! ... extend them to the supercell
    !
    DO ik = 1, num_kpts
      !
      ikevc = ik + ikstart - 1
      CALL davcio( evc, 2*nwordwfc, iunwfc, ikevc, -1 )
      npw = ngk(ik)
      kvec(:) = xk(:,ik)
      !
      counter = 0
      !
      DO ibnd = 1, nbnd
        !
        IF ( excluded_band(ibnd) ) CYCLE
        !
        counter = counter + 1
        !
        psic(:) = ( 0.D0, 0.D0 )
        psicx(:) = ( 0.D0, 0.D0 )
        psic( dffts%nl(igk_k(1:npw,ik)) ) = evc(1:npw,ibnd)
        IF( gamma_only ) psic( dffts%nlm(igk_k(1:npw,ik)) ) = CONJG(evc(1:npw,ibnd))
        CALL invfft( 'Wave', psic, dffts )
        !
        ! ... here we extend the wfc to the whole supercell
        !
        ! ... NB: the routine extend_wfc applies also the phase factor
        ! ...     e^(ikr) so the output wfc (psicx) is a proper Bloch
        ! ...     function and not just its periodic part
        !
        CALL extend_wfc( psic, psicx, dfftcp, kvec )
        !
        ! ... calculate the total density in the supercell
        ! NB: for some reason the weights wg sum up to 2 (and not to 1!),
        !     maybe it includes the spin (?). The result is that, although
        !     here we are considering only the psicx of one spin channel,
        !     the charge corresponding to this density is total and so it 
        !     matches directly nelec (and not simply nelup)
        !
        IF ( calc_rho ) &
          rhor(:) = rhor(:) + ( DBLE( psicx(:) )**2 + &
                               AIMAG( psicx(:) )**2 ) * wg(ibnd,ik) / omega
        !
        CALL fwfft( 'Wave', psicx, dfftcp )
        evcx(1:npwxcp,counter) = psicx( dfftcp%nl(1:npwxcp) )
        !
      ENDDO ! ibnd
      !
      IF ( counter .ne. num_bands ) &
        CALL errore( 'wan2odd', 'wrong number of included bands', counter )
      !
      ! ... save the extended wavefunctions into the buffer
      !
      CALL save_buffer( evcx, nwordwfcx, iunwfcx, ik )
      !
    ENDDO ! ik
    !
    !
    IF ( calc_rho ) THEN
      !
      rhog(:,:) = ( 0.D0, 0.D0 )
      CALL fwfft( 'Rho', rhor, dfftcp )
      rhog(1:ngmcp,1) = rhor( dfftcp%nl(1:ngmcp) )
      !
      CALL check_rho( rhog )
      !
    ENDIF
    !
    !
    ! ... here the Wannier functions are realized
    ! w_Rn(G) = sum_k e^(-ikR) sum_m U_mn(k)*psi_km(G) / Nk^(1/2)
    !
    CALL open_buffer( iunwann, 'wann', nwordwann, io_level, exst )
    !
    ir = 0
    !
    DO i = 1, kgrid(1)
      DO j = 1, kgrid(2)
        DO k = 1, kgrid(3)
          !
          ir = ir + 1
          ewan(:,:) = ( 0.D0, 0.D0 )
          !
          rvec(:) = (/ i-1, j-1, k-1 /)
          CALL cryst_to_cart( 1, rvec, at, 1 )
          !
          DO iw = 1, num_wann
            DO ik = 1, num_kpts
              !
              ! ... phase factor e^(-ikR)
              !
              kvec(:) = xk(:,ik)
              dot_prod = tpi * SUM( kvec(:) * rvec(:) )
              phase = CMPLX( COS(dot_prod), -SIN(dot_prod), KIND=DP )
              !
              ! ... read the supercell-extended Bloch functions
              !
              evcx(:,:) = ( 0.D0, 0.D0 )
              CALL get_buffer( evcx, nwordwfcx, iunwfcx, ik )
              !
              ! ... selecting disentangled bands
              !
              IF ( have_disentangled ) THEN
                !
                num_inc = ndimwin(ik)
                counter = 0
                !
                DO n = 1, num_bands
                  IF ( lwindow(n,ik) ) THEN
                    counter = counter + 1
                    evcx_dis(:,counter) = evcx(:,n)
                  ENDIF
                ENDDO
                !
                IF ( counter .ne. num_inc ) &
                  CALL errore( 'wan2odd', 'Wrong number of included bands &
                                           in disentanglement', counter )
                !
              ENDIF
              !
              ! ... calculate the Wannier function (ir,iw)
              !
              DO ip = 1, num_wann
                !
                ! ... applies disentanglement optimal matrix
                !
                evcw(:) = ( 0.D0, 0.D0 )
                IF ( have_disentangled ) THEN
                  DO n = 1, num_inc
                    evcw(:) = evcw(:) + u_mat_opt(n,ip,ik) * evcx_dis(:,n)
                  ENDDO
                ELSE
                  evcw(:) = evcx(:,ip)
                ENDIF
                !
                ewan(:,iw) = ewan(:,iw) + phase * u_mat(ip,iw,ik) * evcw(:) / SQRT(DBLE(num_kpts)) 
                !
              ENDDO
              !
            ENDDO ! ik
            !
!            CALL check_complex_wfc( ewan(:,iw), wf_is_cmplx, ratio )
!            !
!            ! ... if gamma_only_x=.true. the Wannier functions must be real;
!            ! ... if one of the realized WFs is found to be complex then the
!            ! ... code will restart without using the gamma-trick (complex wfc)
!            ! 
!            IF ( gamma_only_x .and. wf_is_cmplx ) THEN
!              !
!              WRITE( stdout, 20 ) ir, iw
!              WRITE( stdout, 21 ) ratio
!              CALL errore( 'wan2odd', 'complex Wannier functions are incompatible with gamma_trick=.true.', 1 )
!              !
!            ENDIF
            !
            ! ... recalculate the total density from the WFs
            !
            IF ( calc_rho ) THEN
              !
              psicx(:) = ( 0.D0, 0.D0 )
              psicx( dfftcp%nl(1:npwxcp) ) = ewan(1:npwxcp,iw)
              IF( gamma_only_x ) psicx( dfftcp%nlm(1:npwxcp) ) = CONJG(ewan(1:npwxcp,iw))
              CALL invfft( 'Wave', psicx, dfftcp )
              rhow(:) = rhow(:) + DBLE( psicx(:) )**2 + AIMAG( psicx(:) )**2
              !
            ENDIF
            !
          ENDDO ! iw
          !
          CALL save_buffer( ewan, nwordwann, iunwann, ir )
          !
        ENDDO
      ENDDO
    ENDDO
    !
    !
    CALL close_buffer( iunwfcx, 'delete' )
    !
    ! ... checks the consistency between the Wannier and Bloch
    ! ... densities and write the G-space density to file
    !
    IF ( calc_rho ) THEN
      !
      rhowg(:,:) = ( 0.D0, 0.D0 )
      CALL fwfft( 'Rho', rhow, dfftcp )
      rhowg(1:ngmcp,1) = rhow( dfftcp%nl(1:ngmcp) ) / omega_cp
      IF ( nspin == 1 ) rhowg(:,1) = rhowg(:,1) * 2    !!! factor 2 for the other spin-component
      CALL check_rho( rhowg, rhog )
      !
      dirname = './'
      IF ( my_pool_id == 0 .AND. my_bgrp_id == root_bgrp_id ) &
           CALL write_rhog( TRIM(dirname) // "charge-density-x", &
           root_bgrp, intra_bgrp_comm, &
           bg_cp(:,1)*tpiba, bg_cp(:,2)*tpiba, bg_cp(:,3)*tpiba, &
           gamma_only_x, mill_cp, ig_l2g_cp, rhowg(:,:) )
      !
    ENDIF
    !
    !
    ! ... write the WFs to a CP-Koopmans-readable file
    !
    CALL write_wannier_cp( iunwann, nwordwann, npwxcp, num_wann, &
                           num_kpts, ig_l2g_cp, split_evc_file )
    !
    IF ( .not. plot ) CALL close_buffer( iunwann, 'delete' )
    !
    CALL stop_clock( 'wannier2odd' )
    !
    !
    10 FORMAT( /, 2X, 'WARNING: gamma_trick=.true. forces the Wannier functions to be real.', /, &
                 11X, 'DO THIS WITH CAUTION !', / )
    20 FORMAT( 5X, 'The Wannier function ( ', I4, ', ', I4, ' ) is found to be complex' )
    21 FORMAT( 5X, 'Maximum Im/Re ratio = ', E15.8 )
    !
    !
  END SUBROUTINE wan2odd
  !
  !
  !---------------------------------------------------------------------
  SUBROUTINE check_rho( rhog, rhogref )
    !-------------------------------------------------------------------
    !
    ! ...  this routine performs some checks on the supercell total density:
    ! ...  1) the total charge
    ! ...  2) (if rhogref is present) it checks that rhog matches with rhogref
    !
    USE mp,                  ONLY : mp_sum
    USE mp_bands,            ONLY : intra_bgrp_comm
    USE klist,               ONLY : nelec
    USE scf,                 ONLY : rho
    USE constants,           ONLY : eps6
    USE fft_supercell,       ONLY : gstart_cp, omega_cp, check_fft, ngmcp
    USE read_wannier,        ONLY : num_kpts
    !
    !
    IMPLICIT NONE
    !
    COMPLEX(DP), INTENT(IN) :: rhog(:,:) 
    COMPLEX(DP), INTENT(IN), OPTIONAL :: rhogref(:,:)
    !
    REAL(DP) :: nelec_, charge
    INTEGER :: ik
    !
    !
    ! ... check the total charge
    !
    charge = 0.D0
    IF ( gstart_cp == 2 ) THEN
      charge = rhog(1,1) * omega_cp
    ENDIF
    !
    CALL mp_sum( charge, intra_bgrp_comm )
    !
    nelec_ = nelec * num_kpts
    IF ( check_fft ) nelec_ = nelec
    IF ( ABS( charge - nelec_ ) > 1.D-3 * charge ) &
         CALL errore( 'wan2odd', 'wrong total charge', 1 )
    !
    !
    ! ... check rho(G) when dfftcp is taken equal to dffts
    !
    IF ( check_fft ) THEN
      DO ik = 1, ngmcp
        IF ( ABS( DBLE(rhog(ik,1) - rho%of_g(ik,1)) ) .ge. eps6 .or. &
             ABS( AIMAG(rhog(ik,1) - rho%of_g(ik,1)) ) .ge. eps6 ) THEN
          CALL errore( 'wan2odd', 'rhog and rho%of_g differ', ik )
        ENDIF
      ENDDO
    ENDIF
    !
    !
    ! ... when present, rhogref is compared to rhog
    !
    IF ( PRESENT(rhogref) ) THEN
      DO ik = 1, ngmcp
        IF ( ABS( DBLE(rhog(ik,1) - rhogref(ik,1)) ) .ge. eps6 .or. &
             ABS( AIMAG(rhog(ik,1) - rhogref(ik,1)) ) .ge. eps6 ) THEN
          CALL errore( 'wan2odd', 'rhog and rhogref differ', ik )
        ENDIF
      ENDDO 
    ENDIF
    !
    !
  END SUBROUTINE check_rho
  !
  !
  !---------------------------------------------------------------------
  SUBROUTINE check_complex_wfc( evc, wf_is_cmplx, ratio )
    !-------------------------------------------------------------------
    !
    ! ...  this routine checks whether the input wavefunction (in PWs)
    ! ...  is real or complex: if the maximum value of the Im/Re ratio
    ! ...  of the wavefunction on the real grid is bigger the chosen
    ! ...  threshold, evc is considered complex
    !
    ! ... WIP: NOT ABLE TO REPRODUCE WANNIER90  maximum Im/Re ratio yet !!!
    !
    USE mp_bands,            ONLY : intra_bgrp_comm
    USE mp,                  ONLY : mp_max
    USE fft_interfaces,      ONLY : invfft, fwfft
    USE fft_supercell,       ONLY : dfftcp, npwxcp
    !
    !
    IMPLICIT NONE
    !
    COMPLEX(DP), INTENT(IN) :: evc(:) 
    LOGICAL, INTENT(OUT) :: wf_is_cmplx
    REAL(DP), INTENT(OUT) :: ratio
    !
    COMPLEX(DP), ALLOCATABLE :: psic(:)
    REAL(DP) :: thr=1.D-3
    !
    !
    wf_is_cmplx = .false.
    !
    ALLOCATE( psic(dfftcp%nnr) )
    psic(:) = ( 0.D0, 0.D0 )
    psic( dfftcp%nl(1:npwxcp) ) = evc(1:npwxcp)
    CALL invfft( 'Wave', psic, dfftcp )
    !
    ratio = MAXVAL( ABS( AIMAG(psic(:)) / DBLE(psic(:)) ) )
    CALL mp_max( ratio, intra_bgrp_comm )
    !
    IF ( ratio .gt. thr ) wf_is_cmplx = .true.
    !
    !
  END SUBROUTINE check_complex_wfc
  !
  !
  !---------------------------------------------------------------------
  SUBROUTINE dealloc_w2odd
    !-------------------------------------------------------------------
    !
    !
    IMPLICIT NONE
    !
    !
    IF ( ALLOCATED(evcx) ) DEALLOCATE( evcx )
    IF ( ALLOCATED(evcw) ) DEALLOCATE( evcw )
    IF ( ALLOCATED(ewan) ) DEALLOCATE( ewan )
    IF ( ALLOCATED(psicx) ) DEALLOCATE( psicx )
    IF ( ALLOCATED(rhor) ) DEALLOCATE( rhor )
    IF ( ALLOCATED(rhog) ) DEALLOCATE( rhog )
    IF ( ALLOCATED(rhow) ) DEALLOCATE( rhow )
    IF ( ALLOCATED(rhowg) ) DEALLOCATE( rhowg )
    IF ( ALLOCATED(evcx_dis) ) DEALLOCATE( evcx_dis )
    !
    !
  END SUBROUTINE dealloc_w2odd
  !
  !
END MODULE wannier2odd
