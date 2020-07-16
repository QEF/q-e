! Copyrigh(C) 2005-2018 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------------
MODULE exx
  !-----------------------------------------------------------------------------
  !! Variables and subroutines for calculation of exact-exchange contribution.  
  !! Implements ACE: Lin Lin, J. Chem. Theory Comput. 2016, 12, 2242.  
  !! Contains code for band parallelization over pairs of bands: see T. Barnes,
  !! T. Kurth, P. Carrier, N. Wichmann, D. Prendergast, P.R.C. Kent, J. Deslippe
  !! Computer Physics Communications 2017, dx.doi.org/10.1016/j.cpc.2017.01.008.
  !
  USE kinds,                ONLY : DP
  USE noncollin_module,     ONLY : noncolin, npol
  USE io_global,            ONLY : ionode, stdout
  !
  USE control_flags,        ONLY : gamma_only, tqr
  USE fft_types,            ONLY : fft_type_descriptor
  USE stick_base,           ONLY : sticks_map, sticks_map_deallocate
  !
  IMPLICIT NONE
  !
  SAVE
  !
  ! ... general purpose vars
  !
  REAL(DP):: exxalfa=0._DP
  !! the parameter multiplying the exact-exchange part
  REAL(DP), ALLOCATABLE :: x_occupation(:,:)
  !! the weights of auxiliary functions in the density matrix
  INTEGER :: x_nbnd_occ
  !! number of bands of auxiliary functions with at least 
  !! some x_occupation > eps_occ
  REAL(DP), PARAMETER :: eps_occ = 1.d-8
  !! occupation threshold
  COMPLEX(DP), ALLOCATABLE :: exxbuff(:,:,:)
  !! Buffers: temporary (complex) buffer for wfc storage
  REAL(DP), ALLOCATABLE :: locbuff(:,:,:)
  !! temporary (real) buffer for wfc storage
  REAL(DP), ALLOCATABLE :: locmat(:,:,:)
  !! buffer for matrix of localization integrals
  REAL(DP), ALLOCATABLE :: exxmat(:,:,:,:)
  !! buffer for matrix of localization integrals (K)
  !
#if defined(__USE_INTEL_HBM_DIRECTIVES)
!DIR$ ATTRIBUTES FASTMEM :: exxbuff
#elif defined(__USE_CRAY_HBM_DIRECTIVES)
!DIR$ memory(bandwidth) exxbuff
#endif
  !
  LOGICAL :: use_ace 
  !! if .TRUE. use Lin Lin's ACE method, if .FALSE. do not use ACE, 
  !! use old algorithm instead
  COMPLEX(DP), ALLOCATABLE :: xi(:,:,:)
  !! ACE projectors
  COMPLEX(DP), ALLOCATABLE :: evc0(:,:,:)
  !! old wfc (G-space) needed to compute fock3
  INTEGER :: nbndproj
  !! 
  LOGICAL :: domat
  !! 
  REAL(DP)::  local_thr 
  !! threshold for Lin Lin's SCDM localized orbitals: discard 
  !! contribution to V_x if overlap between localized orbitals
  !! is smaller than "local_thr".
  !
  ! ... energy related variables
  !
  REAL(DP) :: fock0 = 0.0_DP
  !! sum <old|Vx(old)|old>
  REAL(DP) :: fock1 = 0.0_DP
  !! sum <new|vx(old)|new>
  REAL(DP) :: fock2 = 0.0_DP
  !! sum <new|vx(new)|new>
  REAL(DP) :: fock3 = 0.0_DP
  !! sum <old|vx(new)|old>
  REAL(DP) :: dexx  = 0.0_DP
  !! fock1 - 0.5*(fock2+fock0)
  !
  ! ... custom fft grid and related G-vectors
  !
  TYPE(fft_type_descriptor) :: dfftt 
  LOGICAL :: exx_fft_initialized = .FALSE.
  REAL(kind=DP), DIMENSION(:), POINTER :: ggt => null()
  !! G-vectors in custom gri
  REAL(kind=DP), DIMENSION(:,:),POINTER :: gt => null()
  !! G-vectors in custom grid
  INTEGER :: gstart_t
  !! gstart_t=2 if ggt(1)=0, =1 otherwise
  INTEGER :: npwt
  !! number of plane waves in custom grid (Gamma-only)
  INTEGER :: ngmt_g
  !! Total number of G-vectors in custom grid
  REAL(DP)  :: ecutfock
  !! energy cutoff for custom grid
  INTEGER :: ibnd_start = 0
  !! starting band index used in bgrp parallelization
  INTEGER :: ibnd_end = 0
  !! ending band index used in bgrp parallelization
  INTEGER :: ibnd_buff_start
  !! starting buffer index used in bgrp parallelization
  INTEGER :: ibnd_buff_end
  !! ending buffer index used in bgrp parallelization
  !
 CONTAINS
#define _CX(A)  CMPLX(A,0._dp,kind=DP)
#define _CY(A)  CMPLX(0._dp,-A,kind=DP)
  !
  !------------------------------------------------------------------------
  SUBROUTINE exx_fft_create()
    !------------------------------------------------------------------------
    !! Initialise the custom grid that allows to put the wavefunction
    !! onto the new (smaller) grid for \rho=\psi_{k+q}\psi^*_k and vice versa.  
    !! Set up fft descriptors, including parallel stuff: sticks, planes, etc.
    !
    USE gvecw,          ONLY : ecutwfc
    USE gvect,          ONLY : ecutrho, ngm, g, gg, gstart, mill
    USE cell_base,      ONLY : at, bg, tpiba2
    USE recvec_subs,    ONLY : ggen, ggens
    USE fft_base,       ONLY : smap
    USE fft_types,      ONLY : fft_type_init
    USE symm_base,      ONLY : fft_fact
    USE mp_exx,         ONLY : nproc_egrp, negrp, intra_egrp_comm
    USE mp_bands,       ONLY : nproc_bgrp, intra_bgrp_comm, nyfft
    !
    USE klist,          ONLY : nks, xk
    USE mp_pools,       ONLY : inter_pool_comm
    USE mp,             ONLY : mp_max, mp_sum
    !
    USE control_flags,  ONLY : tqr
    USE realus,         ONLY : qpointlist, tabxx, tabp
    USE exx_band,       ONLY : smap_exx
    USE command_line_options, ONLY : nmany_
    !
    IMPLICIT NONE
    !
    ! ... local variables
    !
    INTEGER :: ik, ngmt
    INTEGER, ALLOCATABLE :: ig_l2gt(:), millt(:,:)
    INTEGER, EXTERNAL :: n_plane_waves
    REAL(DP) :: gkcut, gcutmt
    LOGICAL :: lpara
    !
    IF ( exx_fft_initialized ) RETURN
    !
    ! Initialise the custom grid that allows us to put the wavefunction
    ! onto the new (smaller) grid for \rho=\psi_{k+q}\psi^*_k and vice versa
    !
    ! gkcut is such that all |k+G|^2 < gkcut (in units of (2pi/a)^2)
    ! Note that with k-points, gkcut > ecutwfc/(2pi/a)^2
    ! gcutmt is such that |q+G|^2 < gcutmt
    !
    IF ( gamma_only ) THEN
       gkcut = ecutwfc / tpiba2
       gcutmt = ecutfock / tpiba2
    ELSE
       gkcut = 0.0_DP
       DO ik = 1, nks
          gkcut = MAX(gkcut, SQRT(SUM(xk(:,ik)**2)))
       ENDDO
       CALL mp_max( gkcut, inter_pool_comm )
       ! Alternatively, variable "qnorm" earlier computed in "exx_grid_init"
       ! could be used as follows:
       ! gkcut = ( SQRT(ecutwfc/tpiba2) + qnorm )**2
       gkcut = ( SQRT(ecutwfc/tpiba2) + gkcut )**2
       ! 
       ! ... the following instruction may be needed if ecutfock \simeq ecutwfc
       ! and guarantees that all k+G are included
       !
       gcutmt = MAX(ecutfock/tpiba2, gkcut)
    ENDIF
    !
    ! ... set up fft descriptors, including parallel stuff: sticks, planes, etc.
    !
    IF (negrp == 1) THEN
       !
       ! ... no band parallelization: exx grid is a subgrid of general grid
       !
       lpara = ( nproc_bgrp > 1 )
       CALL fft_type_init( dfftt, smap, "rho", gamma_only, lpara,         &
                           intra_bgrp_comm, at, bg, gcutmt, gcutmt/gkcut, &
                           fft_fact=fft_fact, nyfft=nyfft, nmany=nmany_ )
       CALL ggens( dfftt, gamma_only, at, g, gg, mill, gcutmt, ngmt, gt, ggt )
       gstart_t = gstart
       npwt = n_plane_waves(ecutwfc/tpiba2, nks, xk, gt, ngmt)
       ngmt_g = ngmt
       CALL mp_sum( ngmt_g, intra_bgrp_comm )
       !
    ELSE
       !
       WRITE( 6, "(5X,'Exchange parallelized over bands (',i4,' band groups)')" ) &
              negrp
       lpara = ( nproc_egrp > 1 )
       CALL fft_type_init( dfftt, smap_exx, "rho", gamma_only, lpara,     &
                           intra_egrp_comm, at, bg, gcutmt, gcutmt/gkcut, &
                           fft_fact=fft_fact, nyfft=nyfft, nmany=nmany_ )
       ngmt = dfftt%ngm
       ngmt_g = ngmt
       CALL mp_sum( ngmt_g, intra_egrp_comm )
       ALLOCATE( gt(3,dfftt%ngm) )
       ALLOCATE( ggt(dfftt%ngm)  )
       ALLOCATE( millt(3,dfftt%ngm) )
       ALLOCATE( ig_l2gt(dfftt%ngm) )
       !
       CALL ggen( dfftt, gamma_only, at, bg, gcutmt, ngmt_g, ngmt, &
                  gt, ggt, millt, ig_l2gt, gstart_t )
       !
       DEALLOCATE( ig_l2gt )
       DEALLOCATE( millt )
       npwt = n_plane_waves( ecutwfc/tpiba2, nks, xk, gt, ngmt )
       !
    ENDIF
    ! define clock labels (this enables the corresponding fft too)
    dfftt%rho_clock_label = 'fftc' ; dfftt%wave_clock_label = 'fftcw' 
    !
    WRITE( stdout, '(/5x,"EXX grid: ",i8," G-vectors", 5x,       &
         &   "FFT dimensions: (",i4,",",i4,",",i4,")")') ngmt_g, &
         &   dfftt%nr1, dfftt%nr2, dfftt%nr3
    !
    exx_fft_initialized = .TRUE.
    !
    IF (tqr) THEN
       IF (ecutfock == ecutrho) THEN
          WRITE( stdout, '(5x,"Real-space augmentation: EXX grid -> DENSE grid")' )
          tabxx => tabp
       ELSE
          WRITE( stdout, '(5x,"Real-space augmentation: initializing EXX grid")' )
          CALL qpointlist( dfftt, tabxx )
       ENDIF
    ENDIF
    !
    RETURN
    !
  END SUBROUTINE exx_fft_create
  !
  !
  !------------------------------------------------------------------------
  SUBROUTINE exx_gvec_reinit( at_old )
    !----------------------------------------------------------------------
    !! Re-initialize g-vectors after rescaling.
    !
    USE cell_base,  ONLY : bg
    !
    IMPLICIT NONE
    !
    REAL(DP), INTENT(IN) :: at_old(3,3)
    !! the lattice vectors at the previous step
    !
    ! ... local variables
    !
    INTEGER :: ig
    REAL(DP) :: gx, gy, gz
    !
    ! ... rescale g-vectors
    !
    CALL cryst_to_cart( dfftt%ngm, gt, at_old, -1 )
    CALL cryst_to_cart( dfftt%ngm, gt, bg,     +1 )
    !
    DO ig = 1, dfftt%ngm
       gx = gt(1,ig)
       gy = gt(2,ig)
       gz = gt(3,ig)
       ggt(ig) = gx*gx + gy*gy + gz*gz
    ENDDO
    !
  END SUBROUTINE exx_gvec_reinit
  !
  !
  !------------------------------------------------------------------------
  SUBROUTINE deallocate_exx()
    !------------------------------------------------------------------------
    !! Deallocates exx objects.
    !
    USE becmod,    ONLY : deallocate_bec_type, is_allocated_bec_type, bec_type
    USE us_exx,    ONLY : becxx
    USE exx_base,  ONLY : xkq_collect, index_xkq, index_xk, index_sym, rir, &
                          working_pool, exx_grid_initialized
    !
    IMPLICIT NONE
    !
    INTEGER :: ikq
    !
    exx_grid_initialized = .FALSE.
    !
    IF ( ALLOCATED(index_xkq) )    DEALLOCATE( index_xkq )
    IF ( ALLOCATED(index_xk ) )    DEALLOCATE( index_xk  )
    IF ( ALLOCATED(index_sym) )    DEALLOCATE( index_sym )
    IF ( ALLOCATED(rir) )          DEALLOCATE( rir )
    IF ( ALLOCATED(x_occupation) ) DEALLOCATE( x_occupation )
    IF ( ALLOCATED(xkq_collect ) ) DEALLOCATE( xkq_collect  )
    IF ( ALLOCATED(exxbuff) )      DEALLOCATE( exxbuff )
    IF ( ALLOCATED(locbuff) )      DEALLOCATE( locbuff )
    IF ( ALLOCATED(locmat) )       DEALLOCATE( locmat )
    IF ( ALLOCATED(exxmat) )       DEALLOCATE( exxmat )
    IF ( ALLOCATED(xi)   )         DEALLOCATE( xi   )
    IF ( ALLOCATED(evc0) )         DEALLOCATE( evc0 )
    !
    IF ( ALLOCATED(becxx) ) THEN
      DO ikq = 1, SIZE(becxx)
        IF (is_allocated_bec_type(becxx(ikq))) CALL deallocate_bec_type( becxx(ikq) )
      ENDDO
      !
      DEALLOCATE( becxx )
    ENDIF
    !
    IF ( ALLOCATED(working_pool) )  DEALLOCATE( working_pool )
    !
    exx_fft_initialized = .FALSE.
    IF ( ASSOCIATED(gt)  )  DEALLOCATE( gt  )
    IF ( ASSOCIATED(ggt) )  DEALLOCATE( ggt )
    !
  END SUBROUTINE deallocate_exx
  !
  !
  !------------------------------------------------------------------------
  SUBROUTINE exxinit( DoLoc )
    !------------------------------------------------------------------------
    !! This subroutine is run before the first H_psi() of each iteration. 
    !! It saves the wavefunctions for the right density matrix, in real space.
    !
    USE wavefunctions,        ONLY : evc
    USE io_files,             ONLY : nwordwfc, iunwfc_exx
    USE buffers,              ONLY : get_buffer
    USE wvfct,                ONLY : nbnd, npwx, wg, current_k
    USE klist,                ONLY : ngk, nks, nkstot, xk, wk, igk_k
    USE symm_base,            ONLY : nsym, s, sr
    USE mp_pools,             ONLY : npool, nproc_pool, me_pool, inter_pool_comm
    USE mp_exx,               ONLY : me_egrp, negrp, init_index_over_band,  &
                                     my_egrp_id, inter_egrp_comm,           &
                                     intra_egrp_comm, iexx_start, iexx_end, &
                                     all_start, all_end
    USE mp,                   ONLY : mp_sum, mp_bcast
    USE funct,                ONLY : get_exx_fraction, start_exx,exx_is_active, &
                                     get_screening_parameter, get_gau_parameter
    USE scatter_mod,          ONLY : gather_grid, scatter_grid
    USE fft_interfaces,       ONLY : invfft
    USE uspp,                 ONLY : nkb, vkb, okvan
    USE us_exx,               ONLY : rotate_becxx
    USE paw_variables,        ONLY : okpaw
    USE paw_exx,              ONLY : PAW_init_fock_kernel
    USE mp_orthopools,        ONLY : intra_orthopool_comm
    USE exx_base,             ONLY : nkqs, xkq_collect, index_xk, index_sym,  &
                                     exx_set_symm, rir, working_pool, exxdiv, &
                                     erfc_scrlen, gau_scrlen, exx_divergence
    USE exx_band,             ONLY : change_data_structure, nwordwfc_exx, &
                                     transform_evc_to_exx, igk_exx, evc_exx
    !
    IMPLICIT NONE
    !
    LOGICAL :: DoLoc
    !! TRUE:  Real Array locbuff(ir, nbnd, nkqs);  
    !! FALSE: Complex Array exxbuff(ir, nbnd/2, nkqs).
    !
    ! ... local variables
    !
    INTEGER :: ik, ibnd, i, j, k, ir, isym, ikq, ig
    INTEGER :: ibnd_loop_start
    INTEGER :: ipol, jpol
    REAL(DP), ALLOCATABLE :: occ(:,:)
    COMPLEX(DP),ALLOCATABLE :: temppsic(:)
#if defined(__USE_INTEL_HBM_DIRECTIVES)
!DIR$ ATTRIBUTES FASTMEM :: temppsic
#elif defined(__USE_CRAY_HBM_DIRECTIVES)
!DIR$ memory(bandwidth) temppsic
#endif
    COMPLEX(DP),ALLOCATABLE :: temppsic_nc(:,:), psic_nc(:,:)
    COMPLEX(DP),ALLOCATABLE :: psic_exx(:)
    INTEGER :: nxxs, nrxxs
#if defined(__MPI)
    COMPLEX(DP),ALLOCATABLE  :: temppsic_all(:), psic_all(:)
    COMPLEX(DP), ALLOCATABLE :: temppsic_all_nc(:,:), psic_all_nc(:,:)
#endif
    COMPLEX(DP) :: d_spin(2,2,48)
    INTEGER :: npw, current_ik
    INTEGER, EXTERNAL :: global_kpoint_index
    INTEGER :: ibnd_start_new, ibnd_end_new, max_buff_bands_per_egrp
    INTEGER :: ibnd_exx, evc_offset
    !
    CALL start_clock ('exxinit')
    IF ( Doloc ) THEN
        WRITE(stdout,'(/,5X,"Using localization algorithm with threshold: ",&
                & D10.2)') local_thr
        ! IF (.NOT.gamma_only) CALL errore('exxinit','SCDM with K-points NYI',1)
        IF (okvan .OR. okpaw) CALL errore( 'exxinit','SCDM with USPP/PAW not &
                                           &implemented', 1 )
    ENDIF 
    IF ( use_ace ) &
        WRITE(stdout,'(/,5X,"Using ACE for calculation of exact exchange")') 
    !
    CALL transform_evc_to_exx( 2 )
    !
    ! ... prepare the symmetry matrices for the spin part
    !
    IF (noncolin) THEN
       DO isym = 1, nsym
          CALL find_u( sr(:,:,isym), d_spin(:,:,isym) )
       ENDDO
    ENDIF
    !
    CALL exx_fft_create()
    !
    ! Note that nxxs is not the same as nrxxs in parallel case
    nxxs = dfftt%nr1x * dfftt%nr2x * dfftt%nr3x
    nrxxs = dfftt%nnr
#if defined(__MPI)
    IF (noncolin) THEN
       ALLOCATE( psic_all_nc(nxxs,npol), temppsic_all_nc(nxxs,npol) )
    ELSEIF ( .NOT. gamma_only ) THEN
       ALLOCATE( psic_all(nxxs), temppsic_all(nxxs) )
    ENDIF
#endif
    IF (noncolin) THEN
       ALLOCATE( temppsic_nc(nrxxs, npol), psic_nc(nrxxs, npol) )
    ELSEIF ( .NOT. gamma_only ) THEN
       ALLOCATE( temppsic(nrxxs) )
    ENDIF
    !
    ALLOCATE( psic_exx(nrxxs) )
    !
    IF (.NOT.exx_is_active()) THEN
       !
       erfc_scrlen = get_screening_parameter()
       gau_scrlen = get_gau_parameter()
       exxdiv  = exx_divergence()
       exxalfa = get_exx_fraction()
       !
       CALL start_exx()
    ENDIF
    !
    IF (.NOT. gamma_only) CALL exx_set_symm( dfftt%nr1,  dfftt%nr2,  dfftt%nr3, &
                                             dfftt%nr1x, dfftt%nr2x, dfftt%nr3x )
    ! set occupations of wavefunctions used in the calculation of exchange term
    IF (.NOT. ALLOCATED(x_occupation)) ALLOCATE( x_occupation(nbnd,nkstot) )
    ALLOCATE( occ(nbnd,nks) )
    !
    DO ik = 1, nks
       IF (ABS(wk(ik)) > eps_occ) THEN
          occ(1:nbnd,ik) = wg(1:nbnd,ik) / wk(ik)
       ELSE
          occ(1:nbnd,ik) = 0._DP
       ENDIF
    ENDDO
    !
    CALL poolcollect( nbnd, nks, occ, nkstot, x_occupation )
    !
    DEALLOCATE( occ )
    !
    ! ... find an upper bound to the number of bands with non zero occupation.
    ! Useful to distribute bands among band groups
    !
    x_nbnd_occ = 0
    DO ik = 1, nkstot
       DO ibnd = MAX(1,x_nbnd_occ), nbnd
          IF (ABS(x_occupation(ibnd,ik)) > eps_occ) x_nbnd_occ = ibnd
       ENDDO
    ENDDO
    !
    IF (nbndproj == 0) nbndproj = nbnd
    !
    CALL divide( inter_egrp_comm, x_nbnd_occ, ibnd_start, ibnd_end )
    CALL init_index_over_band( inter_egrp_comm, nbnd, nbnd )
    !
    ! ... this will cause exxbuff to be calculated for every band
    ibnd_start_new = iexx_start
    ibnd_end_new = iexx_end
    !
    IF ( gamma_only ) THEN
        ibnd_buff_start = (ibnd_start_new+1)/2
        ibnd_buff_end   = (ibnd_end_new+1)/2
        max_buff_bands_per_egrp = MAXVAL((all_end(:)+1)/2-(all_start(:)+1)/2)+1
    ELSE
        ibnd_buff_start = ibnd_start_new
        ibnd_buff_end   = ibnd_end_new
        max_buff_bands_per_egrp = MAXVAL(all_end(:)-all_start(:))+1
    ENDIF
    !
    IF (DoLoc) THEN
      !
      IF (gamma_only) THEN
        IF (.NOT. ALLOCATED(locbuff)) ALLOCATE( locbuff(nrxxs*npol,nbnd,nks) )
        IF (.NOT. ALLOCATED(locmat))  ALLOCATE( locmat(nbnd,nbnd,nks) )
        locbuff = 0.0d0
        locmat = 0.0d0
      ELSE 
        IF (.NOT. ALLOCATED(exxbuff)) ALLOCATE( exxbuff(nrxxs*npol,nbnd,nkqs) )
        IF (.NOT. ALLOCATED(exxmat) ) ALLOCATE( exxmat(nbnd,nkqs,nbnd,nks) )
        exxbuff = (0.0d0, 0.0d0)
        exxmat = 0.0d0
      ENDIF
      !
      IF (.NOT. ALLOCATED(evc0)) then 
        ALLOCATE( evc0(npwx*npol,nbndproj,nks) )
        evc0 = (0.0d0,0.0d0)
      ENDIF
      !
    ELSE
      !
      IF (.NOT. ALLOCATED(exxbuff)) THEN
         IF (gamma_only) THEN
            ALLOCATE( exxbuff(nrxxs*npol,ibnd_buff_start:ibnd_buff_start + &
                                          max_buff_bands_per_egrp-1,nkqs) ) ! THIS WORKS as for k
         ELSE
            ALLOCATE( exxbuff(nrxxs*npol,ibnd_buff_start:ibnd_buff_start + &
                                          max_buff_bands_per_egrp-1,nkqs) )
         ENDIF
      ENDIF
    ENDIF
    !
    !assign buffer
    IF(DoLoc) THEN
      IF(gamma_only) THEN
!$omp parallel do collapse(3) default(shared) firstprivate(npol,nrxxs,nkqs, &
!$omp                ibnd_buff_start,ibnd_buff_end) private(ir,ibnd,ikq,ipol)
        DO ikq=1,SIZE(locbuff,3) 
          DO ibnd=1, x_nbnd_occ 
            DO ir=1,nrxxs*npol
              locbuff(ir,ibnd,ikq)=0.0_DP
            ENDDO
          ENDDO
        ENDDO
      ENDIF
    ELSE
!$omp parallel do collapse(3) default(shared) firstprivate(npol,nrxxs,nkqs, &
!$omp                ibnd_buff_start,ibnd_buff_end) private(ir,ibnd,ikq,ipol)
      DO ikq = 1, SIZE(exxbuff,3) 
         DO ibnd = ibnd_buff_start, ibnd_buff_end
            DO ir = 1, nrxxs*npol
               exxbuff(ir,ibnd,ikq) = (0.0_DP,0.0_DP)
            ENDDO
         ENDDO
      ENDDO
      ! the above loops will replaced with the following line soon
      !CALL threaded_memset(exxbuff, 0.0_DP, nrxxs*npol*SIZE(exxbuff,2)*nkqs*2)
    ENDIF
    !
    ! ... This is parallelized over pools. Each pool computes only its k-points
    !
    KPOINTS_LOOP : &
    DO ik = 1, nks
       !
       IF ( nks > 1 ) CALL get_buffer( evc_exx, nwordwfc_exx, iunwfc_exx, ik )
       !
       ! ik         = index of k-point in this pool
       ! current_ik = index of k-point over all pools
       !
       current_ik = global_kpoint_index( nkstot, ik )
       !
       IF_GAMMA_ONLY : &
       IF (gamma_only) THEN
          !
          IF (MOD(iexx_start,2) == 0) THEN
             ibnd_loop_start = iexx_start-1
          ELSE
             ibnd_loop_start = iexx_start
          ENDIF
          !
          evc_offset = 0
          DO ibnd = ibnd_loop_start, iexx_end, 2
             !
             psic_exx(:) = ( 0._DP, 0._DP )
             !
             IF ( ibnd < iexx_end ) THEN
                IF ( ibnd == ibnd_loop_start .AND. MOD(iexx_start,2) == 0 ) THEN
                   DO ig = 1, npwt
                      psic_exx(dfftt%nl(ig))  = ( 0._DP, 1._DP )*evc_exx(ig,1)
                      psic_exx(dfftt%nlm(ig)) = ( 0._DP, 1._DP )*CONJG(evc_exx(ig,1))
                   ENDDO
                   evc_offset = -1
                ELSE
                   DO ig = 1, npwt
                      psic_exx(dfftt%nl(ig))  = evc_exx(ig,ibnd-ibnd_loop_start+evc_offset+1) &
                           + ( 0._DP, 1._DP ) * evc_exx(ig,ibnd-ibnd_loop_start+evc_offset+2)
                      psic_exx(dfftt%nlm(ig)) = CONJG( evc_exx(ig,ibnd-ibnd_loop_start+evc_offset+1) ) &
                           + ( 0._DP, 1._DP ) * CONJG( evc_exx(ig,ibnd-ibnd_loop_start+evc_offset+2) )
                   ENDDO
                ENDIF
             ELSE
                DO ig=1,npwt
                   psic_exx(dfftt%nl (ig)) = evc_exx(ig,ibnd-ibnd_loop_start+evc_offset+1)
                   psic_exx(dfftt%nlm(ig)) = CONJG( evc_exx(ig,ibnd-ibnd_loop_start+evc_offset+1) )
                ENDDO
             ENDIF
             !
             CALL invfft( 'Wave', psic_exx, dfftt )
             !
             IF (DoLoc) THEN
               locbuff(1:nrxxs,ibnd-ibnd_loop_start+evc_offset+1,ik) = DBLE(  psic_exx(1:nrxxs) )
               IF (ibnd-ibnd_loop_start+evc_offset+2 <= nbnd) &
                  locbuff(1:nrxxs,ibnd-ibnd_loop_start+evc_offset+2,ik) = AIMAG( psic_exx(1:nrxxs) )
             ELSE
               exxbuff(1:nrxxs,(ibnd+1)/2,current_ik)=psic_exx(1:nrxxs) 
             ENDIF
             !
          ENDDO
          !
       ELSE IF_GAMMA_ONLY
          !
          npw = ngk (ik)
          IBND_LOOP_K : &
          DO ibnd = iexx_start, iexx_end
             !
             ibnd_exx = ibnd
             IF (noncolin) THEN
!$omp parallel do default(shared) private(ir) firstprivate(nrxxs)
                DO ir = 1, nrxxs
                   temppsic_nc(ir,1) = ( 0._DP, 0._DP )
                   temppsic_nc(ir,2) = ( 0._DP, 0._DP )
                ENDDO
!$omp parallel do default(shared) private(ig) firstprivate(npw,ik,ibnd_exx)
                DO ig = 1, npw
                   temppsic_nc(dfftt%nl(igk_exx(ig,ik)),1) = evc_exx(ig,ibnd-iexx_start+1)
                ENDDO
!$omp end parallel do
                CALL invfft( 'Wave', temppsic_nc(:,1), dfftt )
!$omp parallel do default(shared) private(ig) firstprivate(npw,ik,ibnd_exx,npwx)
                DO ig = 1, npw
                   temppsic_nc(dfftt%nl(igk_exx(ig,ik)),2) = evc_exx(ig+npwx,ibnd-iexx_start+1)
                ENDDO
!$omp end parallel do
                CALL invfft( 'Wave', temppsic_nc(:,2), dfftt )
             ELSE
!$omp parallel do default(shared) private(ir) firstprivate(nrxxs)
                DO ir = 1, nrxxs
                   temppsic(ir) = ( 0._DP, 0._DP )
                ENDDO
!$omp parallel do default(shared) private(ig) firstprivate(npw,ik,ibnd_exx)
                DO ig = 1, npw
                   temppsic(dfftt%nl(igk_exx(ig,ik))) = evc_exx(ig,ibnd-iexx_start+1)
                ENDDO
!$omp end parallel do
                CALL invfft( 'Wave', temppsic, dfftt )
             ENDIF
             !
             DO ikq = 1, nkqs
                !
                IF (index_xk(ikq) /= current_ik) CYCLE
                isym = ABS(index_sym(ikq) )
                !
                IF (noncolin) THEN ! noncolinear
#if defined(__MPI)
                   DO ipol = 1, npol
                      CALL gather_grid( dfftt, temppsic_nc(:,ipol), temppsic_all_nc(:,ipol) )
                   ENDDO
                   !
                   IF ( me_egrp == 0 ) THEN
!$omp parallel do default(shared) private(ir) firstprivate(npol,nxxs)
                      DO ir = 1, nxxs
                         !DIR$ UNROLL_AND_JAM (2)
                         DO ipol = 1, npol
                            psic_all_nc(ir,ipol) = (0.0_DP, 0.0_DP)
                         ENDDO
                      ENDDO
!$omp end parallel do
!$omp parallel do default(shared) private(ir) firstprivate(npol,isym,nxxs) reduction(+:psic_all_nc)
                      DO ir = 1, nxxs
                         !DIR$ UNROLL_AND_JAM (4)
                         DO ipol = 1, npol
                            DO jpol = 1, npol
                               psic_all_nc(ir,ipol) = psic_all_nc(ir,ipol) + &
                                             CONJG(d_spin(jpol,ipol,isym)) * &
                                             temppsic_all_nc(rir(ir,isym),jpol)
                            ENDDO
                         ENDDO
                      ENDDO
!$omp end parallel do
                   ENDIF
                   !
                   DO ipol = 1, npol
                      CALL scatter_grid( dfftt, psic_all_nc(:,ipol), psic_nc(:,ipol) )
                   ENDDO
#else
!$omp parallel do default(shared) private(ir) firstprivate(npol,nxxs)
                   DO ir = 1, nxxs
                      !DIR$ UNROLL_AND_JAM (2)
                      DO ipol = 1, npol
                         psic_nc(ir,ipol) = (0._DP,0._DP)
                      ENDDO
                   ENDDO
!$omp end parallel do
!$omp parallel do default(shared) private(ipol,jpol,ir) firstprivate(npol,isym,nxxs) reduction(+:psic_nc)
                   DO ir = 1, nxxs
                      !DIR$ UNROLL_AND_JAM (4)
                      DO ipol = 1, npol
                         DO jpol = 1, npol
                            psic_nc(ir,ipol) = psic_nc(ir,ipol) + CONJG(d_spin(jpol,ipol,isym))* &
                                               temppsic_nc(rir(ir,isym),jpol)
                         ENDDO
                      ENDDO
                   ENDDO
!$omp end parallel do
#endif
                   IF (index_sym(ikq) > 0 ) THEN
                      ! sym. op. without time reversal: normal case
!$omp parallel do default(shared) private(ir) firstprivate(ibnd,isym,ikq)
                      DO ir = 1, nrxxs
                         exxbuff(ir,ibnd,ikq) = psic_nc(ir,1)
                         exxbuff(ir+nrxxs,ibnd,ikq) = psic_nc(ir,2)
                      ENDDO
!$omp end parallel do
                   ELSE
                      ! sym. op. with time reversal: spin 1->2*, 2->-1*
!$omp parallel do default(shared) private(ir) firstprivate(ibnd,isym,ikq)
                      DO ir = 1, nrxxs
                         exxbuff(ir,ibnd,ikq) = CONJG(psic_nc(ir,2))
                         exxbuff(ir+nrxxs,ibnd,ikq) = -CONJG(psic_nc(ir,1))
                      ENDDO
!$omp end parallel do
                   ENDIF
                ELSE ! noncolinear
#if defined(__MPI)
                   CALL gather_grid( dfftt, temppsic, temppsic_all )
                   IF ( me_egrp == 0 ) THEN
!$omp parallel do default(shared) private(ir) firstprivate(isym)
                      DO ir = 1, nxxs
                         psic_all(ir) = temppsic_all(rir(ir,isym))
                      ENDDO
!$omp end parallel do
                   ENDIF
                   CALL scatter_grid( dfftt, psic_all, psic_exx )
#else
!$omp parallel do default(shared) private(ir) firstprivate(isym)
                   DO ir = 1, nrxxs
                      psic_exx(ir) = temppsic(rir(ir,isym))
                   ENDDO
!$omp end parallel do
#endif
!$omp parallel do default(shared) private(ir) firstprivate(isym,ibnd,ikq)
                   DO ir = 1, nrxxs
                      IF (index_sym(ikq) < 0 ) THEN
                         psic_exx(ir) = CONJG(psic_exx(ir))
                      ENDIF
                      exxbuff(ir,ibnd,ikq) = psic_exx(ir)
                   ENDDO
!$omp end parallel do
                   !
                ENDIF ! noncolinear
                !
             ENDDO
             !
          ENDDO&
          IBND_LOOP_K
          !
       ENDIF&
       IF_GAMMA_ONLY
    ENDDO&
    KPOINTS_LOOP
    !
    DEALLOCATE( psic_exx )
    IF (noncolin) THEN
       DEALLOCATE( temppsic_nc, psic_nc )
#if defined(__MPI)
       DEALLOCATE( temppsic_all_nc, psic_all_nc )
#endif
    ELSE IF ( .NOT. gamma_only ) THEN
       DEALLOCATE( temppsic )
#if defined(__MPI)
       DEALLOCATE( temppsic_all, psic_all )
#endif
    ENDIF
    !
    ! Each wavefunction in exxbuff is computed by a single pool, collect among 
    ! pools in a smart way (i.e. without doing all-to-all sum and bcast)
    ! See also the initialization of working_pool in exx_mp_init
    ! Note that in Gamma-only LSDA can be parallelized over two pools, and there
    ! is no need to communicate anything: each pools deals with its own spin
    !
    IF ( .NOT. gamma_only ) THEN
       DO ikq = 1, nkqs
         CALL mp_bcast( exxbuff(:,:,ikq), working_pool(ikq), intra_orthopool_comm ) 
       ENDDO
    ENDIF
    !
    ! For US/PAW only: compute <beta_I|psi_j,k+q> for the entire 
    ! de-symmetrized k+q grid by rotating the ones from the irreducible wedge
    !
    IF (okvan) CALL rotate_becxx( nkqs, index_xk, index_sym, xkq_collect )
    !
    ! Initialize 4-wavefunctions one-center Fock integrals
    !    \int \psi_a(r)\phi_a(r)\phi_b(r')\psi_b(r')/|r-r'|
    !
    IF (okpaw) CALL PAW_init_fock_kernel()
    !
    CALL change_data_structure( .FALSE. )
    !
    CALL stop_clock( 'exxinit' )
    !
  END SUBROUTINE exxinit
  !
  !
  !-----------------------------------------------------------------------
  SUBROUTINE vexx( lda, n, m, psi, hpsi, becpsi )
    !-----------------------------------------------------------------------
    !! Wrapper routine computing V_x\psi, V_x = exchange potential. 
    !! Calls generic version vexx_k or Gamma-specific one vexx_gamma.
    !
    USE becmod,         ONLY : bec_type
    USE uspp,           ONLY : okvan
    USE paw_variables,  ONLY : okpaw
    USE us_exx,         ONLY : becxx
    USE mp_exx,         ONLY : negrp, inter_egrp_comm, init_index_over_band
    USE wvfct,          ONLY : nbnd
    USE exx_band,       ONLY : transform_psi_to_exx, transform_hpsi_to_local, &
                               psi_exx, hpsi_exx, igk_exx
    !
    IMPLICIT NONE
    !
    INTEGER :: lda
    !! input: leading dimension of arrays psi and hpsi
    INTEGER :: n
    !! input: true dimension of psi and hpsi
    INTEGER :: m
    !! input: number of states psi
    COMPLEX(DP) :: psi(lda*npol,m)
    !! input: m wavefunctions
    COMPLEX(DP) :: hpsi(lda*npol,m)
    !! output: V_x*psi
    TYPE(bec_type), OPTIONAL :: becpsi
    !! input: <beta|psi>, optional but needed for US and PAW case
    !
    INTEGER :: i
    !
    IF ((okvan.OR.okpaw) .AND. .NOT. PRESENT(becpsi)) &
       CALL errore( 'vexx','becpsi needed for US/PAW case', 1 )
    CALL start_clock( 'vexx' )
    !
    IF (negrp > 1) THEN
       CALL init_index_over_band( inter_egrp_comm, nbnd, m )
       !
       ! ... transform psi to the EXX data structure
       CALL transform_psi_to_exx( lda, n, m, psi )
    ENDIF
    !
    ! ... calculate the EXX contribution to hpsi
    !
    IF ( gamma_only ) THEN
       IF (negrp == 1)THEN
          CALL vexx_gamma( lda, n, m, psi, hpsi, becpsi )
       ELSE
          CALL vexx_gamma( lda, n, m, psi_exx, hpsi_exx, becpsi )
       ENDIF
    ELSE
       IF (negrp.eq.1)THEN
          CALL vexx_k( lda, n, m, psi, hpsi, becpsi )
       ELSE
          CALL vexx_k( lda, n, m, psi_exx, hpsi_exx, becpsi )
       ENDIF
    ENDIF
    !
    IF (negrp > 1) THEN
       !
       ! ... transform hpsi to the local data structure
       !
       CALL transform_hpsi_to_local(lda,n,m,hpsi)
       !
    ENDIF
    !
    CALL stop_clock( 'vexx' )
    !
  END SUBROUTINE vexx
  !
  !
  !-----------------------------------------------------------------------
  SUBROUTINE vexx_gamma( lda, n, m, psi, hpsi, becpsi )
    !-----------------------------------------------------------------------
    !! Gamma-specific version of vexx.
    !
    USE constants,      ONLY : fpi, e2, pi
    USE cell_base,      ONLY : omega
    USE gvect,          ONLY : ngm, g
    USE wvfct,          ONLY : npwx, current_k, nbnd
    USE klist,          ONLY : xk, nks, nkstot, igk_k
    USE fft_interfaces, ONLY : fwfft, invfft
    USE becmod,         ONLY : bec_type
    USE mp_exx,         ONLY : inter_egrp_comm, my_egrp_id, &
                               intra_egrp_comm, me_egrp, &
                               negrp, max_pairs, egrp_pairs, ibands, nibands, &
                               iexx_istart, iexx_iend, &
                               all_start, all_end, iexx_start, jblock
    USE mp,             ONLY : mp_sum, mp_barrier, mp_circular_shift_left
    USE uspp,           ONLY : nkb, okvan
    USE paw_variables,  ONLY : okpaw
    USE us_exx,         ONLY : bexg_merge, becxx, addusxx_g, addusxx_r, &
                               newdxx_g, newdxx_r, add_nlxx_pot, &
                               qvan_init, qvan_clean
    USE paw_exx,        ONLY : PAW_newdxx
    USE exx_base,       ONLY : nqs, index_xkq, index_xk, xkq_collect, &
                               coulomb_fac, g2_convolution_all
    USE exx_band,       ONLY : result_sum, igk_exx
    !
    IMPLICIT NONE
    !
    INTEGER :: lda
    !! input: leading dimension of arrays psi and hpsi
    INTEGER :: n
    !! input: true dimension of psi and hpsi
    INTEGER :: m
    !! input: number of states psi
    COMPLEX(DP) :: psi(lda*npol,m)
    !! input: m wavefunctions
    COMPLEX(DP) :: hpsi(lda*npol,m)
    !! output: V_x*psi
    TYPE(bec_type), OPTIONAL :: becpsi ! or call a calbec(...psi) instead
    !! input: <beta|psi>, optional but needed for US and PAW case
    !
    ! ... local variables
    !
    COMPLEX(DP), ALLOCATABLE :: result(:,:)
    REAL(DP), ALLOCATABLE :: temppsic_dble (:)
    REAL(DP), ALLOCATABLE :: temppsic_aimag(:)
    !
    COMPLEX(DP), ALLOCATABLE :: vc(:,:), deexx(:,:)
    REAL(DP), ALLOCATABLE :: fac(:)
    INTEGER :: ibnd, ik, im , ikq, iq, ipol
    INTEGER :: ir, ig
    INTEGER :: current_ik
    INTEGER :: ibnd_loop_start
    INTEGER :: nrxxs
    REAL(DP) :: x1, x2, xkp(3)
    REAL(DP) :: xkq(3)
    INTEGER, EXTERNAL :: global_kpoint_index
    INTEGER :: ialloc
    COMPLEX(DP), ALLOCATABLE :: big_result(:,:)
    INTEGER :: iproc, nproc_egrp, ii, ipair
    INTEGER :: jbnd, jstart, jend
    ! ... scratch space for fft of psi and rho
    COMPLEX(DP), ALLOCATABLE :: psi_rhoc_work(:)
    INTEGER :: jblock_start, jblock_end
    INTEGER :: iegrp, wegrp
    INTEGER :: exxbuff_index
    INTEGER :: ending_im
    !
    ialloc = nibands(my_egrp_id+1)
    !
    ALLOCATE( fac(dfftt%ngm) )
    !
    nrxxs = dfftt%nnr
    !
    !ALLOCATE( result(nrxxs), temppsic_DBLE(nrxxs), temppsic_aimag(nrxxs) )
    ALLOCATE( result(nrxxs,ialloc), temppsic_DBLE(nrxxs) )
    ALLOCATE( temppsic_aimag(nrxxs) )
    ALLOCATE( psi_rhoc_work(nrxxs) )
    !
    ALLOCATE( vc(nrxxs,ialloc) )
    IF (okvan) ALLOCATE( deexx(nkb,ialloc) )
    !
    current_ik = global_kpoint_index( nkstot, current_k )
    xkp = xk(:,current_k)
    !
    ALLOCATE( big_result(n,m) )
    big_result = 0.0_DP
    result = 0.0_DP
    !
    DO ii = 1, nibands(my_egrp_id+1)
       IF (okvan) deexx(:,ii) = 0.0_DP
    ENDDO
    !
    ! Here the loops start
    !
    INTERNAL_LOOP_ON_Q : &
    DO iq = 1, nqs
       !
       ikq  = index_xkq(current_ik,iq)
       ik   = index_xk(ikq)
       xkq  = xkq_collect(:,ikq)
       !
       ! calculate the 1/|r-r'| (actually, k+q+g) factor and place it in fac
       CALL g2_convolution_all( dfftt%ngm, gt, xkp, xkq, iq, current_k )
       IF ( okvan .AND..NOT.tqr ) CALL qvan_init( dfftt%ngm, xkq, xkp )
       !
       DO iegrp = 1, negrp
          !
          ! compute the id of group whose data is currently worked on
          wegrp = MOD(iegrp+my_egrp_id-1, negrp)+1
          !
          jblock_start = all_start(wegrp)
          jblock_end   = all_end(wegrp)
          !
          LOOP_ON_PSI_BANDS : &
          DO ii = 1,  nibands(my_egrp_id+1)
             !
             ibnd = ibands(ii,my_egrp_id+1)
             !
             IF (ibnd==0 .OR. ibnd>m) CYCLE
             !
             IF (MOD(ii,2) == 1) THEN
                !
                psi_rhoc_work = (0._DP,0._DP)
                !
                IF ((ii+1) <= MIN(m,nibands(my_egrp_id+1))) THEN
                   ! deal with double bands
!$omp parallel do  default(shared), private(ig)
                   DO ig = 1, npwt
                      psi_rhoc_work( dfftt%nl(ig) )  =       psi(ig, ii) + (0._DP,1._DP) * psi(ig, ii+1)
                      psi_rhoc_work( dfftt%nlm(ig) ) = CONJG(psi(ig, ii) - (0._DP,1._DP) * psi(ig, ii+1))
                   ENDDO
!$omp end parallel do
                ENDIF
                !
                IF ( ii == MIN(m,nibands(my_egrp_id+1)) ) THEN
                   ! deal with a single last band
!$omp parallel do  default(shared), private(ig)
                   DO ig = 1, npwt
                      psi_rhoc_work( dfftt%nl(ig) )  =       psi(ig,ii)
                      psi_rhoc_work( dfftt%nlm(ig) ) = CONJG(psi(ig,ii))
                   ENDDO
!$omp end parallel do
                ENDIF
                !
                CALL invfft( 'Wave', psi_rhoc_work, dfftt )
!$omp parallel do default(shared), private(ir)
                DO ir = 1, nrxxs
                   temppsic_DBLE(ir)  = DBLE( psi_rhoc_work(ir) )
                   temppsic_aimag(ir) = AIMAG( psi_rhoc_work(ir) )
                ENDDO
!$omp end parallel do
                !
             ENDIF
             !
             !
             !determine which j-bands to calculate
             jstart = 0
             jend = 0
             DO ipair = 1, max_pairs
                IF (egrp_pairs(1,ipair,my_egrp_id+1) == ibnd) THEN
                   IF (jstart == 0)THEN
                      jstart = egrp_pairs(2,ipair,my_egrp_id+1)
                      jend = jstart
                   ELSE
                      jend = egrp_pairs(2,ipair,my_egrp_id+1)
                   ENDIF
                ENDIF
             ENDDO
             !
             jstart = MAX(jstart,jblock_start)
             jend = MIN(jend,jblock_end)
             !
             IF (MOD(jstart,2) ==0 ) THEN
                ibnd_loop_start = jstart-1
             ELSE
                ibnd_loop_start = jstart
             ENDIF
             !
             IBND_LOOP_GAM : &
             DO jbnd = ibnd_loop_start, jend, 2 !for each band of psi
                !
                exxbuff_index = (jbnd+1)/2-(all_start(wegrp)+1)/2+(iexx_start+1)/2
                !
                IF ( jbnd < jstart ) THEN
                   x1 = 0.0_DP
                ELSE
                   x1 = x_occupation(jbnd,ik)
                ENDIF
                IF ( jbnd == jend) THEN
                   x2 = 0.0_DP
                ELSE
                   x2 = x_occupation(jbnd+1,ik)
                ENDIF
                IF ( ABS(x1) < eps_occ .AND. ABS(x2) < eps_occ ) CYCLE
                !
                ! ... calculate rho in real space. Gamma tricks are used.
                ! temppsic is real; tempphic contains one band in the real part,
                ! another one in the imaginary part; the same applies to rhoc
                !
                IF ( MOD(ii,2) == 0 ) THEN
!$omp parallel do default(shared), private(ir)
                   DO ir = 1, nrxxs
                      psi_rhoc_work(ir) = exxbuff(ir,exxbuff_index,ikq) * temppsic_aimag(ir) / omega
                   ENDDO
!$omp end parallel do
                ELSE
!$omp parallel do default(shared), private(ir)
                   DO ir = 1, nrxxs
                      psi_rhoc_work(ir) = exxbuff(ir,exxbuff_index,ikq) * temppsic_DBLE(ir) / omega
                   ENDDO
!$omp end parallel do
                ENDIF
                !
                ! ... bring rho to G-space
                !
                !   >>> add augmentation in REAL SPACE here
                IF (okvan .AND. tqr) THEN
                   IF (jbnd >= jstart) &
                        CALL addusxx_r( psi_rhoc_work, &
                       _CX(becxx(ikq)%r(:,jbnd)), _CX(becpsi%r(:,ibnd)))
                   IF (jbnd < jend) &
                        CALL addusxx_r( psi_rhoc_work, &
                       _CY(becxx(ikq)%r(:,jbnd+1)),_CX(becpsi%r(:,ibnd)))
                ENDIF
                !
                CALL fwfft( 'Rho', psi_rhoc_work, dfftt )
                !   >>>> add augmentation in G SPACE here
                IF (okvan .AND. .NOT. tqr) THEN
                   ! ... contribution from one band added to real (in real space) part of rhoc
                   IF (jbnd >= jstart) &
                        CALL addusxx_g( dfftt, psi_rhoc_work, xkq,  xkp, 'r', &
                        becphi_r=becxx(ikq)%r(:,jbnd), becpsi_r=becpsi%r(:,ibnd) )
                   ! ... contribution from following band added to imaginary (in real space) part of rhoc
                   IF (jbnd < jend) &
                        CALL addusxx_g( dfftt, psi_rhoc_work, xkq,  xkp, 'i', &
                        becphi_r=becxx(ikq)%r(:,jbnd+1), becpsi_r=becpsi%r(:,ibnd) )
                ENDIF
                !   >>>> charge density done
                !
                vc(:,ii) = 0._DP
                !
!$omp parallel do default(shared), private(ig)
                DO ig = 1, dfftt%ngm
                   !
                   vc(dfftt%nl(ig),ii)  = coulomb_fac(ig,iq,current_k) * psi_rhoc_work(dfftt%nl(ig))
                   vc(dfftt%nlm(ig),ii) = coulomb_fac(ig,iq,current_k) * psi_rhoc_work(dfftt%nlm(ig))
                   !
                ENDDO
!$omp end parallel do
                !
                !   >>>>  compute <psi|H_fock G SPACE here
                IF (okvan .AND. .NOT. tqr) THEN
                   IF (jbnd >= jstart) &
                        CALL newdxx_g( dfftt, vc(:,ii), xkq, xkp, 'r', deexx(:,ii), &
                           becphi_r=x1*becxx(ikq)%r(:,jbnd) )
                   IF (jbnd<jend) &
                        CALL newdxx_g( dfftt, vc(:,ii), xkq, xkp, 'i', deexx(:,ii), &
                            becphi_r=x2*becxx(ikq)%r(:,jbnd+1) )
                ENDIF
                !
                !brings back v in real space
                CALL invfft( 'Rho', vc(:,ii), dfftt )
                !
                !   >>>>  compute <psi|H_fock REAL SPACE here
                IF (okvan .AND. tqr) THEN
                   IF (jbnd >= jstart) &
                        CALL newdxx_r( dfftt,vc(:,ii), _CX(x1*becxx(ikq)%r(:,jbnd)), deexx(:,ii) )
                   IF (jbnd < jend) &
                        CALL newdxx_r( dfftt,vc(:,ii), _CY(x2*becxx(ikq)%r(:,jbnd+1)), deexx(:,ii) )
                ENDIF
                !
                IF (okpaw) THEN
                   IF (jbnd >= jstart) &
                        CALL PAW_newdxx( x1/nqs, _CX(becxx(ikq)%r(:,jbnd)),&
                                                _CX(becpsi%r(:,ibnd)), deexx(:,ii) )
                   IF (jbnd < jend) &
                        CALL PAW_newdxx( x2/nqs, _CX(becxx(ikq)%r(:,jbnd+1)),&
                                                _CX(becpsi%r(:,ibnd)), deexx(:,ii) )
                ENDIF
                !
                ! ... accumulates over bands and k points
                !
!$omp parallel do default(shared), private(ir)
                DO ir = 1, nrxxs
                   result(ir,ii) = result(ir,ii) &
                                 + x1* DBLE(vc(ir,ii))* DBLE(exxbuff(ir,exxbuff_index,ikq)) &
                                 + x2*AIMAG(vc(ir,ii))*AIMAG(exxbuff(ir,exxbuff_index,ikq))
                ENDDO
!$omp end parallel do
                !
             ENDDO &
             IBND_LOOP_GAM
             !
          ENDDO &
          LOOP_ON_PSI_BANDS
          !
          ! get the next nbnd/negrp data
          IF (negrp > 1) CALL mp_circular_shift_left( exxbuff(:,:,ikq), me_egrp, inter_egrp_comm )
          !
       ENDDO ! iegrp
       IF (okvan .AND. .NOT.tqr) CALL qvan_clean()
    ENDDO &
    INTERNAL_LOOP_ON_Q
    !
    DO ii = 1, nibands(my_egrp_id+1)
       !
       ibnd = ibands(ii,my_egrp_id+1)
       !
       IF (ibnd == 0 .OR. ibnd>m) CYCLE
       !
       IF (okvan) THEN
          CALL mp_sum( deexx(:,ii), intra_egrp_comm )
       ENDIF
       !
       ! ... brings back result in G-space
       !
       CALL fwfft( 'Wave' , result(:,ii), dfftt )
       !
       ! ... communicate result
       DO ig = 1, n
          big_result(ig,ibnd) = big_result(ig,ibnd) - exxalfa*result(dfftt%nl(igk_exx(ig,current_k)),ii)
       ENDDO
       !
       ! ... add non-local \sum_I |beta_I> \alpha_Ii (the sum on i is outside)
       IF (okvan) CALL add_nlxx_pot( lda, big_result(:,ibnd), xkp, n, &
                                     igk_exx(1,current_k), deexx(:,ii), eps_occ, exxalfa )
    ENDDO
    !
    CALL result_sum( n*npol, m, big_result )
    !
    IF (iexx_istart(my_egrp_id+1) > 0) THEN
       IF (negrp == 1) THEN
          ending_im = m
       ELSE
          ending_im = iexx_iend(my_egrp_id+1) - iexx_istart(my_egrp_id+1) + 1
       END IF
       DO im = 1, ending_im
!$omp parallel do default(shared), private(ig) firstprivate(im,n)
           DO ig = 1, n
              hpsi(ig,im) = hpsi(ig,im) + big_result(ig,im+iexx_istart(my_egrp_id+1)-1)
           ENDDO
!$omp end parallel do
       ENDDO
    ENDIF
    !
    DEALLOCATE( big_result )
    DEALLOCATE( result, temppsic_dble, temppsic_aimag )
    DEALLOCATE( psi_rhoc_work )
    DEALLOCATE( vc, fac )
    IF (okvan) DEALLOCATE( deexx )
    !
  END SUBROUTINE vexx_gamma
  !
  !
  !-----------------------------------------------------------------------
  SUBROUTINE vexx_k( lda, n, m, psi, hpsi, becpsi )
    !-----------------------------------------------------------------------
    !! Generic, k-point version of vexx.
    !
    USE constants,      ONLY : fpi, e2, pi
    USE cell_base,      ONLY : omega
    USE gvect,          ONLY : ngm, g
    USE wvfct,          ONLY : npwx, current_k, nbnd
    USE klist,          ONLY : xk, nks, nkstot
    USE fft_interfaces, ONLY : fwfft, invfft
    USE becmod,         ONLY : bec_type
    USE mp_exx,         ONLY : inter_egrp_comm, my_egrp_id, negrp, &
                               intra_egrp_comm, me_egrp, &
                               max_pairs, egrp_pairs, ibands, nibands, &
                               max_ibands, iexx_istart, iexx_iend, &
                               all_start, all_end, iexx_start, jblock
    USE mp,             ONLY : mp_sum, mp_barrier, mp_circular_shift_left
    USE uspp,           ONLY : nkb, okvan
    USE paw_variables,  ONLY : okpaw
    USE us_exx,         ONLY : bexg_merge, becxx, addusxx_g, addusxx_r, &
                               newdxx_g, newdxx_r, add_nlxx_pot, &
                               qvan_init, qvan_clean
    USE paw_exx,        ONLY : PAW_newdxx
    USE exx_base,       ONLY : nqs, xkq_collect, index_xkq, index_xk, &
                               coulomb_fac, g2_convolution_all
    USE exx_band,       ONLY : result_sum, igk_exx
    USE io_global,      ONLY : stdout
    !
    !
    IMPLICIT NONE
    !
    INTEGER :: lda
    !! input: leading dimension of arrays psi and hpsi
    INTEGER :: n
    !! input: true dimension of psi and hpsi
    INTEGER :: m
    !! input: number of states psi
    COMPLEX(DP) :: psi(lda*npol,max_ibands)
    !! input: m wavefunctions
    COMPLEX(DP) :: hpsi(lda*npol,max_ibands)
    !! output: V_x*psi
    TYPE(bec_type), OPTIONAL :: becpsi  ! or call a calbec(...psi) instead
    !! input: <beta|psi>, optional but needed for US and PAW case
    !
    ! ... local variables
    !
    COMPLEX(DP),ALLOCATABLE :: temppsic(:,:), result(:,:)
#if defined(__USE_INTEL_HBM_DIRECTIVES)
!DIR$ ATTRIBUTES FASTMEM :: result
#elif defined(__USE_CRAY_HBM_DIRECTIVES)
!DIR$ memory(bandwidth) result
#endif
    COMPLEX(DP),ALLOCATABLE :: temppsic_nc(:,:,:),result_nc(:,:,:)
    INTEGER :: request_send, request_recv
    !
    COMPLEX(DP),ALLOCATABLE :: deexx(:,:)
    COMPLEX(DP),ALLOCATABLE,TARGET :: rhoc(:,:), vc(:,:)
#if defined(__USE_MANY_FFT)
    COMPLEX(DP),POINTER :: prhoc(:), pvc(:)
#endif
#if defined(__USE_INTEL_HBM_DIRECTIVES)
!DIR$ ATTRIBUTES FASTMEM :: rhoc, vc
#elif defined(__USE_CRAY_HBM_DIRECTIVES)
!DIR$ memory(bandwidth) rhoc, vc
#endif
    REAL(DP),   ALLOCATABLE :: fac(:), facb(:)
    INTEGER :: ibnd, ik, im , ikq, iq, ipol
    INTEGER :: ir, ig, ir_start, ir_end
    INTEGER :: irt, nrt, nblock
    INTEGER :: current_ik
    INTEGER :: ibnd_loop_start
    INTEGER :: nrxxs
    REAL(DP) :: x1, x2, xkp(3), omega_inv, nqs_inv
    REAL(DP) :: xkq(3)
    INTEGER, EXTERNAL :: global_kpoint_index
    DOUBLE PRECISION :: max, tempx
    COMPLEX(DP), ALLOCATABLE :: big_result(:,:)
    INTEGER :: ir_out, ipair, jbnd
    INTEGER :: ii, jstart, jend, jcount, jind
    INTEGER :: ialloc, ending_im
    INTEGER :: ijt, njt, jblock_start, jblock_end
    INTEGER :: iegrp, wegrp
    !
    ialloc = nibands(my_egrp_id+1)
    !
    ALLOCATE( fac(dfftt%ngm) )
    nrxxs= dfftt%nnr
    ALLOCATE( facb(nrxxs) )
    !
    IF (noncolin) THEN
       ALLOCATE( temppsic_nc(nrxxs,npol,ialloc), result_nc(nrxxs,npol,ialloc) )
    ELSE
       ALLOCATE( temppsic(nrxxs,ialloc), result(nrxxs,ialloc) )
    ENDIF
    !
    IF (okvan) ALLOCATE( deexx(nkb,ialloc) )
    !
    current_ik = global_kpoint_index( nkstot, current_k )
    xkp = xk(:,current_k)
    !
    ALLOCATE( big_result(n*npol,m) )
    big_result = 0.0_DP
    !
    ! allocate arrays for rhoc and vc
    ALLOCATE( rhoc(nrxxs,jblock), vc(nrxxs,jblock) )
#if defined(__USE_MANY_FFT)
    prhoc(1:nrxxs*jblock) => rhoc(:,:)
    pvc(1:nrxxs*jblock) => vc(:,:)
#endif
    !
    DO ii = 1, nibands(my_egrp_id+1)
       !
       ibnd = ibands(ii,my_egrp_id+1)
       !
       IF (ibnd==0 .OR. ibnd>m) CYCLE
       IF (okvan) deexx(:,ii) = 0._DP
       !
       IF (noncolin) THEN
          temppsic_nc(:,:,ii) = 0._DP
       ELSE
!$omp parallel do  default(shared), private(ir) firstprivate(nrxxs)
          DO ir = 1, nrxxs
             temppsic(ir,ii) = 0._DP
          ENDDO
       ENDIF
       !
       IF (noncolin) THEN
          !
!$omp parallel do  default(shared), private(ig)
          DO ig = 1, n
             temppsic_nc(dfftt%nl(igk_exx(ig,current_k)),1,ii) = psi(ig,ii)
             temppsic_nc(dfftt%nl(igk_exx(ig,current_k)),2,ii) = psi(npwx+ig,ii)
          ENDDO
!$omp end parallel do
          !
          CALL invfft( 'Wave', temppsic_nc(:,1,ii), dfftt )
          CALL invfft( 'Wave', temppsic_nc(:,2,ii), dfftt )
          !
       ELSE
          !
!$omp parallel do  default(shared), private(ig)
          DO ig = 1, n
             temppsic( dfftt%nl(igk_exx(ig,current_k)), ii ) = psi(ig,ii)
          ENDDO
!$omp end parallel do
          !
          CALL invfft( 'Wave', temppsic(:,ii), dfftt )
          !
       ENDIF
       !
       IF (noncolin) THEN
!$omp parallel do default(shared) firstprivate(nrxxs) private(ir)
          DO ir = 1, nrxxs
             result_nc(ir,1,ii) = 0.0_DP
             result_nc(ir,2,ii) = 0.0_DP
          ENDDO
       ELSE
!$omp parallel do default(shared) firstprivate(nrxxs) private(ir)
          DO ir = 1, nrxxs
             result(ir,ii) = 0.0_DP
          ENDDO
       ENDIF
       !
    ENDDO
    !
    !precompute these guys
    omega_inv = 1.0 / omega
    nqs_inv = 1.0 / nqs
    !
    !------------------------------------------------------------------------!
    ! Beginning of main loop
    !------------------------------------------------------------------------!
    DO iq = 1, nqs
       !
       ikq = index_xkq(current_ik,iq)
       ik  = index_xk(ikq)
       xkq = xkq_collect(:,ikq)
       !
       ! calculate the 1/|r-r'| (actually, k+q+g) factor and place it in fac
       CALL g2_convolution_all( dfftt%ngm, gt, xkp, xkq, iq, current_k )
       !
! JRD - below not threaded
       facb = 0D0
       DO ig = 1, dfftt%ngm
          facb(dfftt%nl(ig)) = coulomb_fac(ig,iq,current_k)
       ENDDO
       !
       IF ( okvan .AND..NOT.tqr ) CALL qvan_init( dfftt%ngm, xkq, xkp )
       !
       DO iegrp = 1, negrp
          !
          ! compute the id of group whose data is currently worked on
          wegrp = MOD(iegrp+my_egrp_id-1, negrp)+1
          njt = (all_end(wegrp)-all_start(wegrp)+jblock)/jblock
          !
          DO ijt = 1, njt
             !
             jblock_start = (ijt - 1) * jblock + all_start(wegrp)
             jblock_end = MIN(jblock_start+jblock-1,all_end(wegrp))
             !
             DO ii = 1, nibands(my_egrp_id+1)
                !
                ibnd = ibands(ii,my_egrp_id+1)
                !
                IF (ibnd==0 .OR. ibnd>m) CYCLE
                !
                !determine which j-bands to calculate
                jstart = 0
                jend = 0
                !
                DO ipair = 1, max_pairs
                   IF (egrp_pairs(1,ipair,my_egrp_id+1) == ibnd) THEN
                      IF (jstart == 0) THEN
                         jstart = egrp_pairs(2,ipair,my_egrp_id+1)
                         jend = jstart
                      ELSE
                         jend = egrp_pairs(2,ipair,my_egrp_id+1)
                      ENDIF
                   ENDIF
                ENDDO
                !
                jstart = MAX( jstart, jblock_start )
                jend = MIN( jend, jblock_end )
                !
                !how many iters
                jcount = jend-jstart+1
                IF (jcount <= 0) CYCLE
                !
                !----------------------------------------------------------------------!
                !INNER LOOP START
                !----------------------------------------------------------------------!
                !
                nblock = 2048
                nrt = (nrxxs+nblock-1)/nblock
                !
!$omp parallel do collapse(2) private(ir_start,ir_end)
                DO irt = 1, nrt
                   DO jbnd = jstart, jend
                      ir_start = (irt - 1) * nblock + 1
                      ir_end = MIN(ir_start+nblock-1, nrxxs)
                      IF (noncolin) THEN
                         DO ir = ir_start, ir_end
                            rhoc(ir,jbnd-jstart+1) = ( CONJG(exxbuff(ir,jbnd-all_start(wegrp)+ &
                                                       iexx_start,ikq))*temppsic_nc(ir,1,ii) + &
                                                       CONJG(exxbuff(nrxxs+ir,jbnd-all_start(wegrp)+ &
                                                       iexx_start,ikq))*temppsic_nc(ir,2,ii) )/omega
                         ENDDO
                      ELSE
!DIR$ vector nontemporal (rhoc)
                         DO ir = ir_start, ir_end
                            rhoc(ir,jbnd-jstart+1) = CONJG(exxbuff(ir,jbnd-all_start(wegrp)+ &
                                                     iexx_start,ikq))*temppsic(ir,ii)*omega_inv
                         ENDDO
                      ENDIF
                   ENDDO
                ENDDO
!$omp end parallel do
                !
                !   ... add augmentation in REAL space HERE
                IF (okvan .AND. tqr) THEN ! augment the "charge" in real space
                   DO jbnd = jstart, jend
                      CALL addusxx_r(rhoc(:,jbnd-jstart+1), becxx(ikq)%k(:,jbnd), becpsi%k(:,ibnd))
                   ENDDO
                ENDIF
                !
                !   ... brings it to G-space
#if defined(__USE_MANY_FFT)
                CALL fwfft( 'Rho', prhoc, dfftt, howmany=jcount )
#else
                DO jbnd=jstart, jend
                   CALL fwfft( 'Rho', rhoc(:,jbnd-jstart+1), dfftt )
                ENDDO
#endif
                !
                !   ... add augmentation in G space HERE
                IF (okvan .AND. .NOT. tqr) THEN
                   DO jbnd = jstart, jend
                      CALL addusxx_g( dfftt, rhoc(:,jbnd-jstart+1), xkq, xkp, &
                      'c', becphi_c=becxx(ikq)%k(:,jbnd),becpsi_c=becpsi%k(:,ibnd) )
                   ENDDO
                ENDIF
                !   ... charge done
                !
!call start_collection()
!$omp parallel do collapse(2) private(ir_start,ir_end)
                DO irt = 1, nrt
                   DO jbnd = jstart, jend
                      ir_start = (irt - 1) * nblock + 1
                      ir_end = MIN(ir_start+nblock-1,nrxxs)
!DIR$ vector nontemporal (vc)
                      DO ir = ir_start, ir_end
                         vc(ir,jbnd-jstart+1) = facb(ir) * rhoc(ir,jbnd-jstart+1)*&
                                                x_occupation(jbnd,ik) * nqs_inv
                      ENDDO
                   ENDDO
                ENDDO
!$omp end parallel do
!call stop_collection()
                !
                ! Add ultrasoft contribution (RECIPROCAL SPACE)
                ! compute alpha_I,j,k+q = \sum_J \int <beta_J|phi_j,k+q> V_i,j,k,q Q_I,J(r) d3r
                IF (okvan .AND. .NOT. tqr) THEN
                   DO jbnd=jstart, jend
                      CALL newdxx_g( dfftt, vc(:,jbnd-jstart+1), xkq, xkp, 'c',&
                                     deexx(:,ii), becphi_c=becxx(ikq)%k(:,jbnd) )
                   ENDDO
                ENDIF
                !
                !brings back v in real space
#if defined(__USE_MANY_FFT)
                !fft many
                CALL invfft( 'Rho', pvc, dfftt, howmany=jcount )
#else
                DO jbnd = jstart, jend
                   CALL invfft( 'Rho', vc(:,jbnd-jstart+1), dfftt )
                ENDDO
#endif
                !
                ! Add ultrasoft contribution (REAL SPACE)
                IF (okvan .AND. tqr) THEN
                   DO jbnd = jstart, jend
                      CALL newdxx_r( dfftt, vc(:,jbnd-jstart+1), becxx(ikq)%k(:,jbnd),deexx(:,ii) )
                   ENDDO
                ENDIF
                !
                ! ... Add PAW one-center contribution
                IF (okpaw) THEN
                   DO jbnd = jstart, jend
                      CALL PAW_newdxx( x_occupation(jbnd,ik)/nqs, becxx(ikq)%k(:,jbnd), &
                                       becpsi%k(:,ibnd), deexx(:,ii) )
                   ENDDO
                ENDIF
                !
                ! ... accumulates over bands and k points
                !
!call start_collection()
!$omp parallel do private(ir_start,ir_end)
                DO irt = 1, nrt
                   DO jbnd = jstart, jend
                      ir_start = (irt - 1) * nblock + 1
                      ir_end = MIN(ir_start+nblock-1, nrxxs)
                      IF (noncolin) THEN
                         DO ir = ir_start, ir_end
                            result_nc(ir,1,ii) = result_nc(ir,1,ii) + vc(ir,jbnd-jstart+1) * &
                                                 exxbuff(ir,jbnd-all_start(wegrp)+iexx_start,ikq)
                            result_nc(ir,2,ii) = result_nc(ir,2,ii) + vc(ir,jbnd-jstart+1) * &
                                                 exxbuff(ir+nrxxs,jbnd-all_start(wegrp)+iexx_start,ikq)
                         ENDDO
                      ELSE
!!dir$ vector nontemporal (result)
                         DO ir = ir_start, ir_end
                            result(ir,ii) = result(ir,ii) + vc(ir,jbnd-jstart+1)* &
                                            exxbuff(ir,jbnd-all_start(wegrp)+iexx_start,ikq)
                         ENDDO
                      ENDIF
                   ENDDO
                ENDDO
!$omp end parallel do
!call stop_collection()
                !
                !----------------------------------------------------------------------!
                !INNER LOOP END
                !----------------------------------------------------------------------!
                !
             ENDDO !I-LOOP
          ENDDO !IJT
          !
          ! get the next nbnd/negrp data
          IF (negrp > 1) CALL mp_circular_shift_left( exxbuff(:,:,ikq), me_egrp, inter_egrp_comm )
          !
       ENDDO !iegrp
       !
       IF ( okvan .AND..NOT.tqr ) CALL qvan_clean()
       !
    ENDDO
    !
    !
    !
    DO ii = 1, nibands(my_egrp_id+1)
       !
       ibnd = ibands(ii,my_egrp_id+1)
       !
       IF (ibnd==0 .OR. ibnd>m) CYCLE
       !
       IF (okvan) THEN
          CALL mp_sum( deexx(:,ii), intra_egrp_comm )
       ENDIF
       !
       IF (noncolin) THEN
          !brings back result in G-space
          CALL fwfft( 'Wave', result_nc(:,1,ii), dfftt )
          CALL fwfft( 'Wave', result_nc(:,2,ii), dfftt )
          !
          DO ig = 1, n
             big_result(ig,ibnd) = big_result(ig,ibnd) - exxalfa* &
                                   result_nc(dfftt%nl(igk_exx(ig,current_k)),1,ii)
             big_result(n+ig,ibnd) = big_result(n+ig,ibnd) - exxalfa* &
                                     result_nc(dfftt%nl(igk_exx(ig,current_k)),2,ii)
          ENDDO
       ELSE
          !
          CALL fwfft( 'Wave', result(:,ii), dfftt )
          !
          DO ig = 1, n
             big_result(ig,ibnd) = big_result(ig,ibnd) - exxalfa* &
                                   result(dfftt%nl(igk_exx(ig,current_k)),ii)
          ENDDO
       ENDIF
       !
       ! add non-local \sum_I |beta_I> \alpha_Ii (the sum on i is outside)
       IF (okvan) CALL add_nlxx_pot( lda, big_result(:,ibnd), xkp, n, igk_exx(:,current_k), &
                                    deexx(:,ii), eps_occ, exxalfa )
       !
    ENDDO
    !
    !deallocate temporary arrays
    DEALLOCATE( rhoc, vc )
    !
    !sum result
    CALL result_sum( n*npol, m, big_result )
    !
    IF (iexx_istart(my_egrp_id+1) > 0) THEN
       !
       IF (negrp == 1) THEN
          ending_im = m
       ELSE
          ending_im = iexx_iend(my_egrp_id+1) - iexx_istart(my_egrp_id+1) + 1
       ENDIF
       !
       IF (noncolin) THEN
          DO im = 1, ending_im
!$omp parallel do default(shared), private(ig) firstprivate(im,n)
             DO ig = 1, n
                hpsi(ig,im) = hpsi(ig,im) + big_result(ig,im+iexx_istart(my_egrp_id+1)-1)
             ENDDO
!$omp end parallel do
!$omp parallel do default(shared), private(ig) firstprivate(im,n)
             DO ig = 1, n
                hpsi(lda+ig,im) = hpsi(lda+ig,im) + big_result(n+ig,im+iexx_istart(my_egrp_id+1)-1)
             ENDDO
!$omp end parallel do
          ENDDO
       ELSE
          DO im = 1, ending_im
!$omp parallel do default(shared), private(ig) firstprivate(im,n)
             DO ig = 1, n
                hpsi(ig,im) = hpsi(ig,im) + big_result(ig,im+iexx_istart(my_egrp_id+1)-1)
             ENDDO
!$omp end parallel do
          ENDDO
       ENDIF
    ENDIF
    !
    IF (noncolin) THEN
       DEALLOCATE( temppsic_nc, result_nc )
    ELSE
       DEALLOCATE( temppsic, result )
    ENDIF
    !
    DEALLOCATE( big_result )
    DEALLOCATE( fac, facb )
    IF (okvan) DEALLOCATE( deexx )
    !
  END SUBROUTINE vexx_k
  !
  !
  !-----------------------------------------------------------------------
  FUNCTION exxenergy()
    !-----------------------------------------------------------------------
    !! NB: This function is meant to give the SAME RESULT as exxenergy2.
    !! It is worth keeping it in the repository because in spite of being
    !! slower it is a simple driver using vexx potential routine so it is
    !! good, from time to time, to replace exxenergy2 with it to check that
    !! everything is ok and energy and potential are consistent as they should.
    !
    USE io_files,               ONLY : iunwfc_exx, nwordwfc
    USE buffers,                ONLY : get_buffer
    USE wvfct,                  ONLY : nbnd, npwx, wg, current_k
    USE gvect,                  ONLY : gstart
    USE wavefunctions,          ONLY : evc
    USE lsda_mod,               ONLY : lsda, current_spin, isk
    USE klist,                  ONLY : ngk, nks, xk
    USE mp_pools,               ONLY : inter_pool_comm
    USE mp_exx,                 ONLY : intra_egrp_comm, intra_egrp_comm, &
                                       negrp
    USE mp,                     ONLY : mp_sum
    USE becmod,                 ONLY : bec_type, allocate_bec_type, &
                                       deallocate_bec_type, calbec
    USE uspp,                   ONLY : okvan,nkb,vkb
    USE exx_band,               ONLY : nwordwfc_exx, igk_exx
    !
    IMPLICIT NONE
    !
    TYPE(bec_type) :: becpsi
    REAL(DP) :: exxenergy,  energy
    INTEGER :: npw, ibnd, ik
    COMPLEX(DP) :: vxpsi(npwx*npol,nbnd), psi(npwx*npol,nbnd)
    !
    exxenergy = 0._DP
    !
    CALL start_clock( 'exxenergy' )
    !
    IF (okvan) CALL allocate_bec_type( nkb, nbnd, becpsi )
    energy = 0._dp
    !
    DO ik = 1, nks
       npw = ngk(ik)
       ! setup variables for usage by vexx (same logic as for H_psi)
       current_k = ik
       IF ( lsda ) current_spin = isk(ik)
       ! end setup
       IF ( nks > 1 ) THEN
          CALL get_buffer( psi, nwordwfc_exx, iunwfc_exx, ik )
       ELSE
          psi(1:npwx*npol,1:nbnd) = evc(1:npwx*npol,1:nbnd)
       ENDIF
       !
       IF (okvan) THEN
          ! prepare the |beta> function at k+q
          CALL init_us_2( npw, igk_exx(1,ik), xk(:,ik), vkb )
          ! compute <beta_I|psi_j> at this k+q point, for all band and all projectors
          CALL calbec( npw, vkb, psi, becpsi, nbnd )
       ENDIF
       !
       vxpsi(:,:) = (0._dp, 0._dp)
       CALL vexx( npwx, npw, nbnd, psi, vxpsi, becpsi )
       !
       DO ibnd = 1, nbnd
          energy = energy + DBLE(wg(ibnd,ik) * dot_product(psi(1:npw,ibnd),vxpsi(1:npw,ibnd)))
          IF (noncolin) energy = energy + &
                  DBLE(wg(ibnd,ik) * dot_product(psi(npwx+1:npwx+npw,ibnd),vxpsi(npwx+1:npwx+npw,ibnd)))
          !
       ENDDO
       IF (gamma_only .AND. gstart == 2) THEN
           DO ibnd = 1, nbnd
              energy = energy - &
                       DBLE(0.5_dp * wg(ibnd,ik) * CONJG(psi(1,ibnd)) * vxpsi(1,ibnd))
           ENDDO
       ENDIF
    ENDDO
    !
    IF (gamma_only) energy = 2 * energy
    !
    CALL mp_sum( energy, intra_egrp_comm )
    CALL mp_sum( energy, inter_pool_comm )
    IF (okvan)  CALL deallocate_bec_type( becpsi )
    !
    exxenergy = energy
    !
    CALL stop_clock( 'exxenergy' )
    !
  END FUNCTION exxenergy
  !
  !
  !-----------------------------------------------------------------------
  FUNCTION exxenergy2()
    !-----------------------------------------------------------------------
    !! Wrapper to \(\texttt{exxenergy2_gamma}\) and \(\texttt{exxenergy2_k}\).
    !
    IMPLICIT NONE
    !
    REAL(DP) :: exxenergy2
    !
    CALL start_clock( 'exxenergy' )
    !
    IF ( gamma_only ) THEN
       exxenergy2 = exxenergy2_gamma()
    ELSE
       exxenergy2 = exxenergy2_k()
    ENDIF
    !
    CALL stop_clock( 'exxenergy' )
    !
  END FUNCTION  exxenergy2
  !
  !
  !-----------------------------------------------------------------------
  FUNCTION exxenergy2_gamma()
    !-----------------------------------------------------------------------
    !
    USE constants,               ONLY : fpi, e2, pi
    USE io_files,                ONLY : iunwfc_exx, nwordwfc
    USE buffers,                 ONLY : get_buffer
    USE cell_base,               ONLY : alat, omega, bg, at, tpiba
    USE symm_base,               ONLY : nsym, s
    USE gvect,                   ONLY : ngm, gstart, g
    USE wvfct,                   ONLY : nbnd, npwx, wg
    USE wavefunctions,           ONLY : evc
    USE klist,                   ONLY : xk, ngk, nks, nkstot
    USE lsda_mod,                ONLY : lsda, current_spin, isk
    USE mp_pools,                ONLY : inter_pool_comm
    USE mp_bands,                ONLY : intra_bgrp_comm
    USE mp_exx,                  ONLY : inter_egrp_comm, my_egrp_id, negrp, &
                                        intra_egrp_comm, me_egrp, &
                                        max_pairs, egrp_pairs, ibands, nibands, &
                                        iexx_istart, iexx_iend, &
                                        all_start, all_end, iexx_start, &
                                        init_index_over_band, jblock
    USE mp,                      ONLY : mp_sum, mp_circular_shift_left
    USE fft_interfaces,          ONLY : fwfft, invfft
    USE gvect,                   ONLY : ecutrho
    USE klist,                   ONLY : wk
    USE uspp,                    ONLY : okvan,nkb,vkb
    USE becmod,                  ONLY : bec_type, allocate_bec_type, &
                                        deallocate_bec_type, calbec
    USE paw_variables,           ONLY : okpaw
    USE paw_exx,                 ONLY : PAW_xx_energy
    USE us_exx,                  ONLY : bexg_merge, becxx, addusxx_g, &
                                        addusxx_r, qvan_init, qvan_clean
    USE exx_base,                ONLY : nqs, xkq_collect, index_xkq, index_xk, &
                                        coulomb_fac, g2_convolution_all
    USE exx_band,                ONLY : igk_exx, change_data_structure, &
                                        transform_evc_to_exx, nwordwfc_exx, &
                                        evc_exx
    !
    IMPLICIT NONE
    !
    REAL(DP)   :: exxenergy2_gamma
    !
    ! ... local variables
    !
    REAL(DP) :: energy
    COMPLEX(DP), ALLOCATABLE :: temppsic(:)
    COMPLEX(DP), ALLOCATABLE :: rhoc(:)
    REAL(DP),    ALLOCATABLE :: fac(:)
    COMPLEX(DP), ALLOCATABLE :: vkb_exx(:,:)
    INTEGER  :: jbnd, ibnd, ik, ikk, ig, ikq, iq, ir
    INTEGER  :: nrxxs, current_ik, ibnd_loop_start
    REAL(DP) :: x1, x2
    REAL(DP) :: xkq(3), xkp(3), vc
    INTEGER, EXTERNAL :: global_kpoint_index
    !
    TYPE(bec_type) :: becpsi
    COMPLEX(DP), ALLOCATABLE :: psi_t(:), prod_tot(:)
    REAL(DP),ALLOCATABLE :: temppsic_dble (:)
    REAL(DP),ALLOCATABLE :: temppsic_aimag(:)
    INTEGER :: npw
    INTEGER :: istart, iend, ipair, ii, ialloc
    INTEGER :: ijt, njt, jblock_start, jblock_end
    INTEGER :: exxbuff_index
    INTEGER :: calbec_start, calbec_end
    INTEGER :: intra_bgrp_comm_
    INTEGER :: iegrp, wegrp
    !
    CALL init_index_over_band( inter_egrp_comm, nbnd, nbnd )
    !
    CALL transform_evc_to_exx( 0 )
    !
    ialloc = nibands(my_egrp_id+1)
    !
    nrxxs = dfftt%nnr
    ALLOCATE( fac(dfftt%ngm) )
    !
    ALLOCATE( temppsic(nrxxs), temppsic_DBLE(nrxxs), temppsic_aimag(nrxxs) )
    ALLOCATE( rhoc(nrxxs) )
    ALLOCATE( vkb_exx(npwx,nkb) )
    !
    energy = 0.0_DP
    !
    CALL allocate_bec_type( nkb, nbnd, becpsi )
    !
    IKK_LOOP : &
    DO ikk = 1, nks
       current_ik = global_kpoint_index( nkstot, ikk )
       xkp = xk(:,ikk)
       !
       IF ( lsda ) current_spin = isk(ikk)
       npw = ngk (ikk)
       IF ( nks > 1 ) CALL get_buffer( evc_exx, nwordwfc_exx, iunwfc_exx, ikk )
       !
       ! ... prepare the |beta> function at k+q
       CALL init_us_2( npw, igk_exx(:,ikk), xkp, vkb_exx )
       !
       ! ... compute <beta_I|psi_j> at this k+q point, for all band and all projectors
       calbec_start = ibands(1,my_egrp_id+1)
       calbec_end = ibands(nibands(my_egrp_id+1),my_egrp_id+1)
       !
       intra_bgrp_comm_ = intra_bgrp_comm
       intra_bgrp_comm = intra_egrp_comm
       !
       CALL calbec( npw, vkb_exx, evc_exx, becpsi, nibands(my_egrp_id+1) )
       !
       intra_bgrp_comm = intra_bgrp_comm_
       !
       IQ_LOOP : &
       DO iq = 1,nqs
          !
          ikq  = index_xkq(current_ik,iq)
          ik   = index_xk(ikq)
          !
          xkq = xkq_collect(:,ikq)
          !
          CALL g2_convolution_all( dfftt%ngm, gt, xkp, xkq, iq, current_ik )
          !
          fac = coulomb_fac(:,iq,current_ik)
          fac(gstart_t:) = 2 * coulomb_fac(gstart_t:,iq,current_ik)
          !
          IF ( okvan .AND. .NOT.tqr ) CALL qvan_init( dfftt%ngm, xkq, xkp )
          !
          DO iegrp = 1, negrp
             !
             ! ... compute the id of group whose data is currently worked on
             wegrp = MOD(iegrp+my_egrp_id-1, negrp)+1
             !
             jblock_start = all_start(wegrp)
             jblock_end   = all_end(wegrp)
             !
             JBND_LOOP : &
             DO ii = 1, nibands(my_egrp_id+1)
                !
                jbnd = ibands(ii,my_egrp_id+1)
                !
                IF (jbnd==0 .OR. jbnd>nbnd) CYCLE
                !
                IF ( MOD(ii,2)==1 ) THEN
                   !
                   temppsic = (0._DP,0._DP)
                   !
                   IF ( (ii+1) <= nibands(my_egrp_id+1) ) THEN
                      ! deal with double bands
!$omp parallel do  default(shared), private(ig)
                      DO ig = 1, npwt
                         temppsic( dfftt%nl(ig) )  = &
                              evc_exx(ig,ii) + (0._DP,1._DP) * evc_exx(ig,ii+1)
                         temppsic( dfftt%nlm(ig) ) = &
                              CONJG(evc_exx(ig,ii) - (0._DP,1._DP) * evc_exx(ig,ii+1))
                      ENDDO
!$omp end parallel do
                   ENDIF
                   !
                   IF (ii == nibands(my_egrp_id+1)) THEN
                      ! deal with a single last band
!$omp parallel do  default(shared), private(ig)
                      DO ig = 1, npwt
                         temppsic( dfftt%nl(ig) ) = evc_exx(ig,ii)
                         temppsic( dfftt%nlm(ig) ) = CONJG(evc_exx(ig,ii))
                      ENDDO
!$omp end parallel do
                   ENDIF
                   !
                   CALL invfft( 'Wave', temppsic, dfftt )
!$omp parallel do default(shared), private(ir)
                   DO ir = 1, nrxxs
                      temppsic_DBLE(ir) = DBLE( temppsic(ir) )
                      temppsic_aimag(ir) = AIMAG( temppsic(ir) )
                   ENDDO
!$omp end parallel do
                   !
                ENDIF
                !
                !determine which j-bands to calculate
                istart = 0
                iend = 0
                !
                DO ipair = 1, max_pairs
                   IF (egrp_pairs(1,ipair,my_egrp_id+1) == jbnd) THEN
                      IF (istart == 0) THEN
                         istart = egrp_pairs(2,ipair,my_egrp_id+1)
                         iend = istart
                      ELSE
                         iend = egrp_pairs(2,ipair,my_egrp_id+1)
                      ENDIF
                   ENDIF
                ENDDO
                !
                istart = MAX(istart,jblock_start)
                iend = MIN(iend,jblock_end)
                !
                IF (MOD(istart,2) == 0) THEN
                   ibnd_loop_start = istart-1
                ELSE
                   ibnd_loop_start = istart
                ENDIF
                !
                IBND_LOOP_GAM : &
                DO ibnd = ibnd_loop_start, iend, 2       !for each band of psi
                   !
                   exxbuff_index = (ibnd+1)/2-(all_start(wegrp)+1)/2+(iexx_start+1)/2
                   !
                   IF ( ibnd < istart ) THEN
                      x1 = 0.0_DP
                   ELSE
                      x1 = x_occupation(ibnd,ik)
                   ENDIF
                   !
                   IF ( ibnd < iend ) THEN
                      x2 = x_occupation(ibnd+1,ik)
                   ELSE
                      x2 = 0.0_DP
                   ENDIF
                   ! calculate rho in real space. Gamma tricks are used.
                   ! temppsic is real; tempphic contains band 1 in the real part,
                   ! band 2 in the imaginary part; the same applies to rhoc
                   !
                   IF ( MOD(ii,2) == 0 ) THEN
                      rhoc = 0.0_DP
!$omp parallel do default(shared), private(ir)
                      DO ir = 1, nrxxs
                         rhoc(ir) = exxbuff(ir,exxbuff_index,ikq) * temppsic_aimag(ir) / omega
                      ENDDO
!$omp end parallel do
                   ELSE
!$omp parallel do default(shared), private(ir)
                      DO ir = 1, nrxxs
                         rhoc(ir) = exxbuff(ir,exxbuff_index,ikq) * temppsic_DBLE(ir) / omega
                      ENDDO
!$omp end parallel do
                   ENDIF
                   !
                   IF (okvan .AND. tqr) THEN
                      IF (ibnd >= istart) &
                           CALL addusxx_r( rhoc,_CX(becxx(ikq)%r(:,ibnd)), &
                                          _CX(becpsi%r(:,jbnd)) )
                      IF (ibnd<iend) &
                           CALL addusxx_r(rhoc,_CY(becxx(ikq)%r(:,ibnd+1)), &
                                          _CX(becpsi%r(:,jbnd)))
                   ENDIF
                   !
                   ! bring rhoc to G-space
                   CALL fwfft( 'Rho', rhoc, dfftt )
                   !
                   IF (okvan .AND. .NOT.tqr) THEN
                      IF (ibnd >= istart) &
                           CALL addusxx_g( dfftt, rhoc, xkq, xkp, 'r', &
                           becphi_r=becxx(ikq)%r(:,ibnd), becpsi_r=becpsi%r(:,jbnd-calbec_start+1) )
                      IF (ibnd < iend) &
                           CALL addusxx_g( dfftt, rhoc, xkq, xkp, 'i', &
                           becphi_r=becxx(ikq)%r(:,ibnd+1), becpsi_r=becpsi%r(:,jbnd-calbec_start+1) )
                   ENDIF
                   !
                   vc = 0.0_DP
!$omp parallel do  default(shared), private(ig),  reduction(+:vc)
                   DO ig = 1, dfftt%ngm
                      !
                      ! The real part of rhoc contains the contribution from band ibnd
                      ! The imaginary part    contains the contribution from band ibnd+1
                      !
                      vc = vc + fac(ig) * ( x1 * &
                           ABS( rhoc(dfftt%nl(ig)) + CONJG(rhoc(dfftt%nlm(ig))) )**2 &
                                 +x2 * &
                           ABS( rhoc(dfftt%nl(ig)) - CONJG(rhoc(dfftt%nlm(ig))) )**2 )
                   ENDDO
!$omp end parallel do
                   !
                   vc = vc * omega * 0.25_DP / nqs
                   energy = energy - exxalfa * vc * wg(jbnd,ikk)
                   !
                   IF (okpaw) THEN
                      IF (ibnd >= ibnd_start) &
                           energy = energy + exxalfa*wg(jbnd,ikk)*&
                           x1 * PAW_xx_energy(_CX(becxx(ikq)%r(:,ibnd)),_CX(becpsi%r(:,jbnd)) )
                      IF (ibnd < ibnd_end) &
                           energy = energy + exxalfa*wg(jbnd,ikk)*&
                           x2 * PAW_xx_energy(_CX(becxx(ikq)%r(:,ibnd+1)), _CX(becpsi%r(:,jbnd)) )
                   ENDIF
                   !
                ENDDO &
                IBND_LOOP_GAM
             ENDDO &
             JBND_LOOP
             !
             ! get the next nbnd/negrp data
             IF (negrp > 1) CALL mp_circular_shift_left( exxbuff(:,:,ikq), me_egrp, inter_egrp_comm )
             !
          ENDDO ! iegrp
          IF ( okvan .AND. .NOT.tqr ) CALL qvan_clean( )
          !
       ENDDO &
       IQ_LOOP
    ENDDO &
    IKK_LOOP
    !
    DEALLOCATE( temppsic, temppsic_dble, temppsic_aimag )
    !
    DEALLOCATE( rhoc, fac )
    CALL deallocate_bec_type( becpsi )
    !
    CALL mp_sum( energy, inter_egrp_comm )
    CALL mp_sum( energy, intra_egrp_comm )
    CALL mp_sum( energy, inter_pool_comm )
    !
    exxenergy2_gamma = energy
    !
    CALL change_data_structure( .FALSE. )
    !
  END FUNCTION  exxenergy2_gamma
  !
  !
  !-----------------------------------------------------------------------
  FUNCTION exxenergy2_k()
    !-----------------------------------------------------------------------
    !
    USE constants,               ONLY : fpi, e2, pi
    USE io_files,                ONLY : iunwfc_exx, nwordwfc
    USE buffers,                 ONLY : get_buffer
    USE cell_base,               ONLY : alat, omega, bg, at, tpiba
    USE symm_base,               ONLY : nsym, s
    USE gvect,                   ONLY : ngm, gstart, g
    USE wvfct,                   ONLY : nbnd, npwx, wg
    USE wavefunctions,           ONLY : evc
    USE klist,                   ONLY : xk, ngk, nks, nkstot
    USE lsda_mod,                ONLY : lsda, current_spin, isk
    USE mp_pools,                ONLY : inter_pool_comm
    USE mp_exx,                  ONLY : inter_egrp_comm, my_egrp_id, negrp, &
                                        intra_egrp_comm, me_egrp, &
                                        max_pairs, egrp_pairs, ibands, nibands, &
                                        iexx_istart, iexx_iend, &
                                        all_start, all_end, iexx_start, &
                                        init_index_over_band, jblock
    USE mp_bands,                ONLY : intra_bgrp_comm
    USE mp,                      ONLY : mp_sum, mp_circular_shift_left
    USE fft_interfaces,          ONLY : fwfft, invfft
    USE gvect,                   ONLY : ecutrho
    USE klist,                   ONLY : wk
    USE uspp,                    ONLY : okvan,nkb,vkb
    USE becmod,                  ONLY : bec_type, allocate_bec_type, &
                                        deallocate_bec_type, calbec
    USE paw_variables,           ONLY : okpaw
    USE paw_exx,                 ONLY : PAW_xx_energy
    USE us_exx,                  ONLY : bexg_merge, becxx, addusxx_g, &
                                        addusxx_r, qvan_init, qvan_clean
    USE exx_base,                ONLY : nqs, xkq_collect, index_xkq, index_xk, &
                                        coulomb_fac, g2_convolution_all
    USE exx_band,                ONLY : change_data_structure, &
                                        transform_evc_to_exx, nwordwfc_exx, &
                                        igk_exx, evc_exx
    !
    IMPLICIT NONE
    !
    REAL(DP)   :: exxenergy2_k
    !
    ! ... local variables
    !
    REAL(DP) :: energy
    COMPLEX(DP), ALLOCATABLE :: temppsic(:,:)
    COMPLEX(DP), ALLOCATABLE :: temppsic_nc(:,:,:)
    COMPLEX(DP), ALLOCATABLE,TARGET :: rhoc(:,:)
#if defined(__USE_MANY_FFT)
    COMPLEX(DP), POINTER :: prhoc(:)
#endif
    REAL(DP),    ALLOCATABLE :: fac(:)
    INTEGER  :: npw, jbnd, ibnd, ibnd_inner_start, ibnd_inner_end, ibnd_inner_count, &
                ik, ikk, ig, ikq, iq, ir
    INTEGER  :: h_ibnd, nrxxs, current_ik, ibnd_loop_start, nblock, nrt, irt, &
                ir_start, ir_end
    REAL(DP) :: x1, x2
    REAL(DP) :: xkq(3), xkp(3), vc, omega_inv
    INTEGER, EXTERNAL :: global_kpoint_index
    !
    TYPE(bec_type) :: becpsi
    COMPLEX(DP), ALLOCATABLE :: psi_t(:), prod_tot(:)
    INTEGER :: intra_bgrp_comm_
    INTEGER :: ii, ialloc, jstart, jend, ipair
    INTEGER :: ijt, njt, jblock_start, jblock_end
    INTEGER :: iegrp, wegrp
    !
    CALL init_index_over_band( inter_egrp_comm, nbnd, nbnd )
    !
    CALL transform_evc_to_exx( 0 )
    !
    ialloc = nibands(my_egrp_id+1)
    !
    nrxxs = dfftt%nnr
    ALLOCATE( fac(dfftt%ngm) )
    !
    IF (noncolin) THEN
       ALLOCATE( temppsic_nc(nrxxs,npol,ialloc) )
    ELSE
       ALLOCATE( temppsic(nrxxs,ialloc) )
    ENDIF
    !
    energy = 0.0_DP
    !
    CALL allocate_bec_type( nkb, nbnd, becpsi )
    !
    !precompute that stuff
    omega_inv = 1.0/omega
    !
    IKK_LOOP : &
    DO ikk = 1, nks
       !
       current_ik = global_kpoint_index ( nkstot, ikk )
       xkp = xk(:,ikk)
       !
       IF ( lsda ) current_spin = isk(ikk)
       npw = ngk(ikk)
       IF ( nks > 1 ) CALL get_buffer( evc_exx, nwordwfc_exx, iunwfc_exx, ikk )
       !
       ! compute <beta_I|psi_j> at this k+q point, for all band and all projectors
       intra_bgrp_comm_ = intra_bgrp_comm
       intra_bgrp_comm = intra_egrp_comm
       !
       IF (okvan .OR. okpaw) THEN
          CALL compute_becpsi( npw, igk_exx(:,ikk), xkp, evc_exx, &
                               becpsi%k(:,ibands(1,my_egrp_id+1)) )
       ENDIF
       !
       intra_bgrp_comm = intra_bgrp_comm_
       !
       ! ... precompute temppsic
       !
       IF (noncolin) THEN
          temppsic_nc = 0.0_DP
       ELSE
          temppsic = 0.0_DP
       ENDIF
       !
       DO ii = 1, nibands(my_egrp_id+1)
          !
          jbnd = ibands(ii,my_egrp_id+1)
          !
          IF (jbnd == 0 .OR. jbnd > nbnd) CYCLE
          !
          !IF ( abs(wg(jbnd,ikk)) < eps_occ) CYCLE
          !
          IF (noncolin) THEN
             !
!$omp parallel do default(shared), private(ig)
             DO ig = 1, npw
                temppsic_nc(dfftt%nl(igk_exx(ig,ikk)),1,ii) = evc_exx(ig,ii)
                temppsic_nc(dfftt%nl(igk_exx(ig,ikk)),2,ii) = evc_exx(npwx+ig,ii)
             ENDDO
!$omp end parallel do
             !
             CALL invfft( 'Wave', temppsic_nc(:,1,ii), dfftt )
             CALL invfft( 'Wave', temppsic_nc(:,2,ii), dfftt )
             !
          ELSE
!$omp parallel do default(shared), private(ig)
             DO ig = 1, npw
                temppsic(dfftt%nl(igk_exx(ig,ikk)),ii) = evc_exx(ig,ii)
             ENDDO
!$omp end parallel do
             !
             CALL invfft( 'Wave', temppsic(:,ii), dfftt )
             !
          ENDIF
       ENDDO
       !
       IQ_LOOP : &
       DO iq = 1,nqs
          !
          ikq  = index_xkq(current_ik,iq)
          ik   = index_xk(ikq)
          !
          xkq = xkq_collect(:,ikq)
          !
          CALL g2_convolution_all( dfftt%ngm, gt, xkp, xkq, iq, ikk )
          IF ( okvan .AND..NOT.tqr ) CALL qvan_init( dfftt%ngm, xkq, xkp )
          !
          DO iegrp = 1, negrp
             !
             ! ... compute the id of group whose data is currently worked on
             wegrp = MOD(iegrp+my_egrp_id-1, negrp)+1
             njt = (all_end(wegrp)-all_start(wegrp)+jblock)/jblock
             !
             IJT_LOOP : &
             DO ijt = 1, njt
                !
                jblock_start = (ijt - 1) * jblock + all_start(wegrp)
                jblock_end = MIN(jblock_start+jblock-1,all_end(wegrp))
                !
                JBND_LOOP : &
                DO ii = 1, nibands(my_egrp_id+1)
                   !
                   jbnd = ibands(ii,my_egrp_id+1)
                   !
                   IF (jbnd==0 .OR. jbnd>nbnd) CYCLE
                   !
                   !determine which j-bands to calculate
                   jstart = 0
                   jend = 0
                   !
                   DO ipair = 1, max_pairs
                      IF (egrp_pairs(1,ipair,my_egrp_id+1) == jbnd) THEN
                         IF (jstart == 0) THEN
                            jstart = egrp_pairs(2,ipair,my_egrp_id+1)
                            jend = jstart
                         ELSE
                            jend = egrp_pairs(2,ipair,my_egrp_id+1)
                         ENDIF
                      ENDIF
                   ENDDO
                   !
                   !these variables prepare for inner band parallelism
                   jstart = MAX(jstart,jblock_start)
                   jend = MIN(jend,jblock_end)
                   ibnd_inner_start = jstart
                   ibnd_inner_end = jend
                   ibnd_inner_count = jend-jstart+1
                   !
                   !allocate arrays
                   ALLOCATE( rhoc(nrxxs,ibnd_inner_count) )
#if defined(__USE_MANY_FFT)
                   prhoc(1:nrxxs*ibnd_inner_count) => rhoc
#endif 
                   !calculate rho in real space
                   nblock = 2048
                   nrt = (nrxxs+nblock-1) / nblock
!$omp parallel do collapse(2) private(ir_start,ir_end)
                   DO irt = 1, nrt
                      DO ibnd = ibnd_inner_start, ibnd_inner_end
                         ir_start = (irt - 1) * nblock + 1
                         ir_end = MIN(ir_start+nblock-1,nrxxs)
                         IF (noncolin) THEN
                            DO ir = ir_start, ir_end
                               rhoc(ir,ibnd-ibnd_inner_start+1) = &
                                 ( CONJG(exxbuff(ir,ibnd-all_start(wegrp)+iexx_start,ikq)) * &
                                 temppsic_nc(ir,1,ii) + &
                                 CONJG(exxbuff(nrxxs+ir,ibnd-all_start(wegrp)+iexx_start,ikq)) * &
                                 temppsic_nc(ir,2,ii) ) * omega_inv
                            ENDDO
                         ELSE
                            DO ir = ir_start, ir_end
                               rhoc(ir,ibnd-ibnd_inner_start+1) = omega_inv * &
                                 CONJG(exxbuff(ir,ibnd-all_start(wegrp)+iexx_start,ikq)) * &
                                 temppsic(ir,ii)
                            ENDDO
                         ENDIF
                      ENDDO
                   ENDDO
!$omp end parallel do
                   !
                   ! augment the "charge" in real space
                   IF (okvan .AND. tqr) THEN
!$omp parallel do default(shared) private(ibnd) firstprivate(ibnd_inner_start,ibnd_inner_end)
                      DO ibnd = ibnd_inner_start, ibnd_inner_end
                         CALL addusxx_r( rhoc(:,ibnd-ibnd_inner_start+1), &
                                        becxx(ikq)%k(:,ibnd), becpsi%k(:,jbnd))
                      ENDDO
!$omp end parallel do
                   ENDIF
                   !
                   ! bring rhoc to G-space
#if defined(__USE_MANY_FFT)
                   CALL fwfft ('Rho', prhoc, dfftt, howmany=ibnd_inner_count)
#else
                   DO ibnd = ibnd_inner_start, ibnd_inner_end
                      CALL fwfft('Rho', rhoc(:,ibnd-ibnd_inner_start+1), dfftt)
                   ENDDO
#endif
                   ! augment the "charge" in G space
                   IF (okvan .AND. .NOT. tqr) THEN
                      DO ibnd = ibnd_inner_start, ibnd_inner_end
                         CALL addusxx_g(dfftt, rhoc(:,ibnd-ibnd_inner_start+1), &
                              xkq, xkp, 'c', becphi_c=becxx(ikq)%k(:,ibnd),     &
                              becpsi_c=becpsi%k(:,jbnd))
                      ENDDO
                   ENDIF
                   !
!$omp parallel do reduction(+:energy) private(vc)
                   DO ibnd = ibnd_inner_start, ibnd_inner_end
                      vc=0.0_DP
                      DO ig=1,dfftt%ngm
                         vc = vc + coulomb_fac(ig,iq,ikk) * &
                             DBLE(rhoc(dfftt%nl(ig),ibnd-ibnd_inner_start+1) *&
                             CONJG(rhoc(dfftt%nl(ig),ibnd-ibnd_inner_start+1)))
                      ENDDO
                      vc = vc * omega * x_occupation(ibnd,ik) / nqs
                      energy = energy - exxalfa * vc * wg(jbnd,ikk)
                      !
                      IF (okpaw) THEN
                         energy = energy +exxalfa*x_occupation(ibnd,ik)/nqs*wg(jbnd,ikk) &
                              *PAW_xx_energy(becxx(ikq)%k(:,ibnd), becpsi%k(:,jbnd))
                      ENDIF
                   ENDDO
!$omp end parallel do
                   !
                   !deallocate memory
                   DEALLOCATE( rhoc )
                ENDDO &
                JBND_LOOP
                !
             ENDDO&
             IJT_LOOP
             ! get the next nbnd/negrp data
             IF (negrp > 1) call mp_circular_shift_left( exxbuff(:,:,ikq), me_egrp, inter_egrp_comm )
             !
          END DO !iegrp
          !
          IF ( okvan .AND. .NOT.tqr ) CALL qvan_clean()
       ENDDO &
       IQ_LOOP
       !
    ENDDO &
    IKK_LOOP
    !
    IF (noncolin) THEN
       DEALLOCATE( temppsic_nc )
    ELSE
       DEALLOCATE( temppsic )
    ENDIF
    !
    DEALLOCATE( fac )
    !
    CALL deallocate_bec_type( becpsi )
    !
    CALL mp_sum( energy, inter_egrp_comm )
    CALL mp_sum( energy, intra_egrp_comm )
    CALL mp_sum( energy, inter_pool_comm )
    !
    exxenergy2_k = energy
    CALL change_data_structure( .FALSE. )
    !
  END FUNCTION  exxenergy2_k
  !
  !
  !-----------------------------------------------------------------------
  FUNCTION exx_stress()
    !-----------------------------------------------------------------------
    !! This is Eq.(10) of PRB 73, 125120 (2006).
    !
    USE constants,            ONLY : fpi, e2, pi, tpi
    USE io_files,             ONLY : iunwfc_exx, nwordwfc
    USE buffers,              ONLY : get_buffer
    USE cell_base,            ONLY : alat, omega, bg, at, tpiba
    USE symm_base,            ONLY : nsym, s
    USE wvfct,                ONLY : nbnd, npwx, wg, current_k
    USE wavefunctions,        ONLY : evc
    USE klist,                ONLY : xk, ngk, nks
    USE lsda_mod,             ONLY : lsda, current_spin, isk
    USE gvect,                ONLY : g
    USE mp_pools,             ONLY : npool, inter_pool_comm
    USE mp_exx,               ONLY : inter_egrp_comm, intra_egrp_comm, &
                                     ibands, nibands, my_egrp_id, jblock, &
                                     egrp_pairs, max_pairs, negrp, me_egrp, &
                                     all_start, all_end, iexx_start
    USE mp,                   ONLY : mp_sum, mp_circular_shift_left
    USE fft_base,             ONLY : dffts
    USE fft_interfaces,       ONLY : fwfft, invfft
    USE uspp,                 ONLY : okvan
    !
    USE exx_base,             ONLY : nq1, nq2, nq3, nqs, eps, exxdiv,       &
                                     x_gamma_extrapolation, on_double_grid, &
                                     grid_factor, yukawa, erfc_scrlen,      &
                                     use_coulomb_vcut_ws, use_coulomb_vcut_spheric, &
                                     gau_scrlen, vcut, index_xkq, index_xk, index_sym
    USE exx_band,             ONLY : change_data_structure, transform_evc_to_exx, &
                                     g_exx, igk_exx, nwordwfc_exx, evc_exx
    USE coulomb_vcut_module,  ONLY : vcut_get,  vcut_spheric_get
    !
    IMPLICIT NONE
    !
    ! ... local variables
    !
    REAL(DP) :: exx_stress(3,3), exx_stress_(3,3)
    !
    COMPLEX(DP),ALLOCATABLE :: tempphic(:), temppsic(:), result(:)
    COMPLEX(DP),ALLOCATABLE :: tempphic_nc(:,:), temppsic_nc(:,:), &
                               result_nc(:,:)
    COMPLEX(DP),ALLOCATABLE :: rhoc(:)
    REAL(DP),   ALLOCATABLE :: fac(:), fac_tens(:,:,:), fac_stress(:)
    INTEGER  :: npw, jbnd, ibnd, ik, ikk, ig, ir, ikq, iq, isym
    INTEGER  :: nqi, iqi, beta, nrxxs, ngm
    INTEGER  :: ibnd_loop_start
    REAL(DP) :: x1, x2
    REAL(DP) :: qq, xk_cryst(3), sxk(3), xkq(3), vc(3,3), x, q(3)
    REAL(DP) :: delta(3,3)
    INTEGER :: jstart, jend, ii, ipair, jblock_start, jblock_end
    INTEGER :: iegrp, wegrp
    INTEGER :: exxbuff_index
    !
    CALL start_clock( 'exx_stress' )
    !
    CALL transform_evc_to_exx( 0 )
    !
    IF (npool>1) CALL errore( 'exx_stress', 'stress not available with pools', 1 )
    IF (noncolin) CALL errore( 'exx_stress', 'noncolinear stress not implemented', 1 )
    IF (okvan) CALL infomsg( 'exx_stress', 'USPP stress not tested' )
    !
    nrxxs = dfftt%nnr
    ngm   = dfftt%ngm
    delta = RESHAPE( (/1._dp,0._dp,0._dp, 0._dp,1._dp,0._dp, 0._dp,0._dp,1._dp/), (/3,3/))
    exx_stress_ = 0._dp
    !
    ALLOCATE( tempphic(nrxxs), temppsic(nrxxs), rhoc(nrxxs), fac(ngm) )
    ALLOCATE( fac_tens(3,3,ngm), fac_stress(ngm) )
    !
    nqi = nqs
    !
    ! ... loop over k-points
    DO ikk = 1, nks
       current_k = ikk
       IF (lsda) current_spin = isk(ikk)
       npw = ngk(ikk)
       !
       IF (nks > 1) CALL get_buffer( evc_exx, nwordwfc_exx, iunwfc_exx, ikk )
       !
       DO iqi = 1, nqi
          !
          iq = iqi
          !
          ikq  = index_xkq(current_k,iq)
          ik   = index_xk(ikq)
          isym = ABS(index_sym(ikq))      
          !      

          ! FIXME: use cryst_to_cart and company as above..      
          xk_cryst(:) = at(1,:)*xk(1,ik)+at(2,:)*xk(2,ik)+at(3,:)*xk(3,ik)      
          IF (index_sym(ikq) < 0) xk_cryst = -xk_cryst      
          sxk(:) = s(:,1,isym)*xk_cryst(1) + &      
                   s(:,2,isym)*xk_cryst(2) + &      
                   s(:,3,isym)*xk_cryst(3)      
          xkq(:) = bg(:,1)*sxk(1) + bg(:,2)*sxk(2) + bg(:,3)*sxk(3)      
          !      
          !CALL start_clock ('exxen2_ngmloop')      
          !      
!$omp parallel do default(shared), private(ig, beta, q, qq, on_double_grid, x)
          DO ig = 1, ngm      
             IF (negrp == 1) THEN      
                q(1) = xk(1,current_k) - xkq(1) + g(1,ig)      
                q(2) = xk(2,current_k) - xkq(2) + g(2,ig)      
                q(3) = xk(3,current_k) - xkq(3) + g(3,ig)      
             ELSE      
                q(1) = xk(1,current_k) - xkq(1) + g_exx(1,ig)      
                q(2) = xk(2,current_k) - xkq(2) + g_exx(2,ig)      
                q(3) = xk(3,current_k) - xkq(3) + g_exx(3,ig)      
             ENDIF      
             !      
             q = q * tpiba      
             qq = ( q(1)*q(1) + q(2)*q(2) + q(3)*q(3) )      
             !      
             DO beta = 1, 3      
                fac_tens(1:3,beta,ig) = q(1:3)*q(beta)      
             ENDDO      
             !      
             IF (x_gamma_extrapolation) THEN      
                on_double_grid = .TRUE.      
                x= 0.5d0/tpiba*(q(1)*at(1,1)+q(2)*at(2,1)+q(3)*at(3,1))*nq1      
                on_double_grid = on_double_grid .AND. (ABS(x-NINT(x))<eps)      
                x= 0.5d0/tpiba*(q(1)*at(1,2)+q(2)*at(2,2)+q(3)*at(3,2))*nq2      
                on_double_grid = on_double_grid .AND. (ABS(x-NINT(x))<eps)      
                x= 0.5d0/tpiba*(q(1)*at(1,3)+q(2)*at(2,3)+q(3)*at(3,3))*nq3      
                on_double_grid = on_double_grid .AND. (ABS(x-NINT(x))<eps)      
             ELSE      
                on_double_grid = .FALSE.      
             ENDIF      
             !      
             IF (use_coulomb_vcut_ws) THEN      
                fac(ig) = vcut_get(vcut, q)      
                fac_stress(ig) = 0._dp   ! not implemented      
                IF (gamma_only .AND. qq > 1.d-8) fac(ig) = 2.d0 * fac(ig)      
                !      
             ELSEIF ( use_coulomb_vcut_spheric ) THEN      
                fac(ig) = vcut_spheric_get(vcut, q)      
                fac_stress(ig) = 0._dp   ! not implemented      
                IF (gamma_only .AND. qq > 1.d-8) fac(ig) = 2.d0 * fac(ig)      
                !      
             ELSEIF (gau_scrlen > 0) THEN      
                fac(ig) = e2*((pi/gau_scrlen)**(1.5d0))* &       
                          EXP(-qq/4.d0/gau_scrlen) * grid_factor       
                fac_stress(ig) =  e2*2.d0/4.d0/gau_scrlen  * &       
                                  EXP(-qq/4.d0/gau_scrlen) * &
                                  ((pi/gau_scrlen)**(1.5d0))*grid_factor       
                IF (gamma_only) fac(ig) = 2.d0 * fac(ig)       
                IF (gamma_only) fac_stress(ig) = 2.d0 * fac_stress(ig)       
                IF (on_double_grid) fac(ig) = 0._dp       
                IF (on_double_grid) fac_stress(ig) = 0._dp       
                !       
             ELSEIF (qq > 1.d-8) THEN      
                IF ( erfc_scrlen > 0 ) THEN       
                  fac(ig)=e2*fpi/qq*(1._dp-EXP(-qq/4.d0/erfc_scrlen**2)) * grid_factor       
                  fac_stress(ig) = -e2*fpi * 2.d0/qq**2 * ( &       
                      (1._dp+qq/4.d0/erfc_scrlen**2)*EXP(-qq/4.d0/erfc_scrlen**2) - 1._dp) * &       
                      grid_factor       
                ELSE       
                  fac(ig)=e2*fpi/( qq + yukawa ) * grid_factor       
                  fac_stress(ig) = 2.d0 * e2*fpi/(qq+yukawa)**2 * grid_factor       
                ENDIF       
                !       
                IF (gamma_only) fac(ig) = 2.d0 * fac(ig)       
                IF (gamma_only) fac_stress(ig) = 2.d0 * fac_stress(ig)       
                IF (on_double_grid) fac(ig) = 0._dp       
                IF (on_double_grid) fac_stress(ig) = 0._dp       
                !      
             ELSE 
                ! 
                fac(ig) = -exxdiv ! or rather something else (see f.gygi)       
                fac_stress(ig) = 0._dp  ! or -exxdiv_stress (not yet implemented)       
                IF ( yukawa> 0._dp .AND. .NOT. x_gamma_extrapolation) THEN       
                  fac(ig) = fac(ig) + e2*fpi/( qq + yukawa )       
                  fac_stress(ig) = 2.d0 * e2*fpi/(qq+yukawa)**2       
                ENDIF       
                IF (erfc_scrlen > 0._dp .AND. .NOT. x_gamma_extrapolation) THEN       
                  fac(ig) = e2*fpi / (4.d0*erfc_scrlen**2)       
                  fac_stress(ig) = e2*fpi / (8.d0*erfc_scrlen**4)       
                ENDIF 
                !
             ENDIF
             !
          ENDDO      
!$omp end parallel do
          !CALL stop_clock ('exxen2_ngmloop')      
          DO iegrp = 1, negrp      
             !      
             ! compute the id of group whose data is currently worked on      
             wegrp = MOD(iegrp+my_egrp_id-1, negrp)+1      
             !      
             jblock_start = all_start(wegrp)      
             jblock_end = all_end(wegrp)      
             !      
             ! loop over bands      
             DO ii = 1, nibands(my_egrp_id+1)      
                !      
                jbnd = ibands(ii,my_egrp_id+1)      
                !      
                IF (jbnd==0 .OR. jbnd>nbnd) CYCLE      
                !      
                !determine which j-bands to calculate      
                jstart = 0      
                jend = 0      
                DO ipair=1, max_pairs      
                   IF (egrp_pairs(1,ipair,my_egrp_id+1).eq.jbnd)THEN      
                      IF (jstart == 0)THEN      
                         jstart = egrp_pairs(2,ipair,my_egrp_id+1)      
                         jend = jstart      
                      ELSE      
                         jend = egrp_pairs(2,ipair,my_egrp_id+1)      
                      ENDIF      
                   ENDIF      
                ENDDO      
                !      
                jstart = MAX(jstart,jblock_start)      
                jend = MIN(jend,jblock_end)      
                !      
                temppsic(:) = ( 0._dp, 0._dp )      
!$omp parallel do default(shared), private(ig)
                DO ig = 1, npw      
                   temppsic(dfftt%nl(igk_exx(ig,ikk))) = evc_exx(ig,ii)      
                ENDDO      
!$omp end parallel do
                !      
                IF (gamma_only) THEN      
!$omp parallel do default(shared), private(ig)
                   DO ig = 1, npw      
                      temppsic(dfftt%nlm(igk_exx(ig,ikk))) = CONJG(evc_exx(ig,ii))      
                   ENDDO      
!$omp end parallel do
                ENDIF      
                !      
                CALL invfft( 'Wave', temppsic, dfftt )      
                !      
                IF (gamma_only) THEN      
                   !      
                   IF (MOD(jstart,2) == 0) THEN      
                      ibnd_loop_start = jstart-1      
                   ELSE      
                      ibnd_loop_start = jstart      
                   ENDIF      
                   !      
                   DO ibnd = ibnd_loop_start, jend, 2     !for each band of psi      
                      !      
                      exxbuff_index = (ibnd+1)/2-(all_start(wegrp)+1)/2+(iexx_start+1)/2      
                      !      
                      IF ( ibnd < jstart ) THEN      
                         x1 = 0._dp      
                      ELSE      
                         x1 = x_occupation(ibnd,ik)      
                      ENDIF      
                      !      
                      IF ( ibnd == jend) THEN      
                         x2 = 0._dp      
                      ELSE      
                         x2 = x_occupation(ibnd+1,ik)      
                      ENDIF      
                      !      
                      IF ( ABS(x1) < eps_occ .AND. ABS(x2) < eps_occ ) CYCLE      
                      !      
                      ! calculate rho in real space      
!$omp parallel do default(shared), private(ir)
                      DO ir = 1, nrxxs      
                         tempphic(ir) = exxbuff(ir,exxbuff_index,ikq)      
                         rhoc(ir) = CONJG(tempphic(ir))*temppsic(ir) / omega      
                      ENDDO      
!$omp end parallel do
                      ! bring it to G-space      
                      CALL fwfft( 'Rho', rhoc, dfftt )      
                      !      
                      vc = 0._dp      
!$omp parallel do default(shared), private(ig), reduction(+:vc)
                      DO ig = 1, ngm      
                         !      
                         vc(:,:) = vc(:,:) + x1 * 0.25_dp * &      
                                   ABS( rhoc(dfftt%nl(ig)) + &      
                                   CONJG(rhoc(dfftt%nlm(ig))))**2 * &      
                                   (fac_tens(:,:,ig)*fac_stress(ig)/2.d0 - delta(:,:)*fac(ig))      
                         vc(:,:) = vc(:,:) + x2 * 0.25_dp * &      
                                   ABS( rhoc(dfftt%nl(ig)) - &      
                                   CONJG(rhoc(dfftt%nlm(ig))))**2 * &      
                                   (fac_tens(:,:,ig)*fac_stress(ig)/2.d0 - delta(:,:)*fac(ig))      
                      ENDDO      
!$omp end parallel do
                      vc = vc / nqs / 4.d0
                      exx_stress_ = exx_stress_ + exxalfa * vc * wg(jbnd,ikk)
                   ENDDO
                   !
                ELSE
                   !
                   DO ibnd = jstart, jend    !for each band of psi
                      !
                      IF ( ABS(x_occupation(ibnd,ik)) < 1.d-6) CYCLE
                      !
                      ! calculate rho in real space
!$omp parallel do default(shared), private(ir)
                      DO ir = 1, nrxxs
                         tempphic(ir) = exxbuff(ir,ibnd-all_start(wegrp)+iexx_start,ikq)
                         rhoc(ir) = CONJG(tempphic(ir))*temppsic(ir) / omega
                      ENDDO
!$omp end parallel do
                      !
                      ! bring it to G-space
                      CALL fwfft( 'Rho', rhoc, dfftt )
                      !
                      vc = 0._dp
!$omp parallel do default(shared), private(ig), reduction(+:vc)
                      DO ig = 1, ngm
                         vc(:,:) = vc(:,:) + rhoc(dfftt%nl(ig))  * &
                                   CONJG(rhoc(dfftt%nl(ig))) *     &
                                   (fac_tens(:,:,ig)*fac_stress(ig)/2.d0 - &
                                   delta(:,:)*fac(ig))
                      ENDDO
!$omp end parallel do
                      !
                      vc = vc * x_occupation(ibnd,ik) / nqs / 4.d0
                      exx_stress_ = exx_stress_ + exxalfa * vc * wg(jbnd,ikk)
                      !
                   ENDDO
                   !
                ENDIF ! gamma or k-points
                !
             ENDDO ! jbnd
             !
             ! get the next nbnd/negrp data
             IF (negrp > 1) CALL mp_circular_shift_left( exxbuff(:,:,ikq), me_egrp, &
                                                         inter_egrp_comm )
             !
          ENDDO ! iegrp
          !
       ENDDO ! iqi
       !
    ENDDO ! ikk
    !
    DEALLOCATE( tempphic, temppsic, rhoc, fac, fac_tens, fac_stress )
    !
    CALL mp_sum( exx_stress_, intra_egrp_comm )
    CALL mp_sum( exx_stress_, inter_egrp_comm )
    CALL mp_sum( exx_stress_, inter_pool_comm )
    !
    exx_stress = exx_stress_
    !
    CALL change_data_structure( .FALSE. )
    !
    CALL stop_clock( 'exx_stress' )
    !
  END FUNCTION exx_stress
  !
  !
  !----------------------------------------------------------------------
  SUBROUTINE compute_becpsi( npw_, igk_, q_, evc_exx, becpsi_k )
  !----------------------------------------------------------------------
  !! Calculates beta functions (Kleinman-Bylander projectors), with
  !! structure factor, for all atoms, in reciprocal space.
  !! FIXME: why so much replicated code?  
  !
  USE kinds,         ONLY : DP
  USE ions_base,     ONLY : nat, ntyp => nsp, ityp, tau
  USE cell_base,     ONLY : tpiba, omega
  USE constants,     ONLY : tpi
  USE gvect,         ONLY : eigts1, eigts2, eigts3, mill, g
  USE wvfct,         ONLY : npwx, nbnd
  USE us,            ONLY : nqx, dq, tab, tab_d2y, spline_ps
  USE m_gth,         ONLY : mk_ffnl_gth
  USE splinelib
  USE uspp,          ONLY : nkb, nhtol, nhtolm, indv
  USE uspp_param,    ONLY : upf, lmaxkb, nhm, nh
  USE becmod,        ONLY : calbec
  USE mp_exx,        ONLY : ibands, nibands, my_egrp_id
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: npw_
  !! number of PWs
  INTEGER, INTENT(IN) :: igk_(npw_)
  !! indices of G in the list of q+G vectors
  REAL(DP), INTENT(IN) :: q_(3)
  !! q vector (2pi/a units)
  COMPLEX(DP), INTENT(IN) :: evc_exx(npwx,nibands(my_egrp_id+1))
  !! wavefunctions from the PW set to exx
  COMPLEX(DP), INTENT(OUT) :: becpsi_k(nkb,nibands(my_egrp_id+1))
  !! <beta|psi> for k points
  !
  ! ... local variables
  !
  COMPLEX(DP) :: vkb_(npwx,1) !beta functions (npw_ <= npwx)
  !
  INTEGER :: i0, i1, i2, i3, ig, lm, na, nt, nb, ih, jkb
  !
  REAL(DP) :: px, ux, vx, wx, arg
  REAL(DP), ALLOCATABLE :: gk(:,:), qg(:), vq(:), ylm(:,:), vkb1(:,:)
  !
  COMPLEX(DP) :: phase, pref
  COMPLEX(DP), ALLOCATABLE :: sk(:)
  !
  REAL(DP), ALLOCATABLE :: xdata(:)
  INTEGER :: iq
  INTEGER :: istart, iend
  !
  istart = ibands(1,my_egrp_id+1)
  iend = ibands(nibands(my_egrp_id+1),my_egrp_id+1)
  !
  IF (lmaxkb < 0) RETURN
  !
  ALLOCATE( vkb1(npw_,nhm) )
  ALLOCATE( sk(npw_) )    
  ALLOCATE( qg(npw_) )    
  ALLOCATE( vq(npw_) )    
  ALLOCATE( ylm(npw_,(lmaxkb+1)**2) )    
  ALLOCATE( gk(3,npw_) )
  !
  ! write(*,'(3i4,i5,3f10.5)') size(tab,1), size(tab,2), size(tab,3), size(vq), q_
  !
  DO ig = 1, npw_
     gk(1,ig) = q_(1) + g(1,igk_(ig))
     gk(2,ig) = q_(2) + g(2,igk_(ig))
     gk(3,ig) = q_(3) + g(3,igk_(ig))
     qg(ig) = gk(1,ig)**2 + gk(2,ig)**2 + gk(3,ig)**2
  ENDDO
  !
  CALL ylmr2( (lmaxkb+1)**2, npw_, gk, qg, ylm )
  !
  ! ... set now qg=|q+G| in atomic units
  !
  DO ig = 1, npw_
     qg(ig) = SQRT(qg(ig))*tpiba
  ENDDO
  !
  IF (spline_ps) THEN
     ALLOCATE( xdata(nqx) )
     DO iq = 1, nqx
       xdata(iq) = (iq - 1) * dq
     ENDDO
  ENDIF
  ! |beta_lm(q)> = (4pi/omega).Y_lm(q).f_l(q).(i^l).S(q)
  jkb = 0
  !
  DO nt = 1, ntyp
     ! ... calculate beta in G-space using an interpolation table:
     ! f_l(q)=\int _0 ^\infty dr r^2 f_l(r) j_l(q.r)
     DO nb = 1, upf(nt)%nbeta
        IF ( upf(nt)%is_gth ) THEN
           CALL mk_ffnl_gth( nt, nb, npw_, omega, qg, vq )
        ELSE
           DO ig = 1, npw_
              IF (spline_ps) THEN
                vq(ig) = splint(xdata, tab(:,nb,nt), tab_d2y(:,nb,nt), qg(ig))
              ELSE
                px = qg (ig) / dq - INT(qg (ig) / dq)
                ux = 1.d0 - px
                vx = 2.d0 - px
                wx = 3.d0 - px
                i0 = INT( qg (ig) / dq ) + 1
                i1 = i0 + 1
                i2 = i0 + 2
                i3 = i0 + 3
                vq (ig) = tab (i0, nb, nt) * ux * vx * wx / 6.d0 + &
                          tab (i1, nb, nt) * px * vx * wx / 2.d0 - &
                          tab (i2, nb, nt) * px * ux * wx / 2.d0 + &
                          tab (i3, nb, nt) * px * ux * vx / 6.d0
              ENDIF
           ENDDO
        ENDIF
        !
        ! ... add spherical harmonic part  (Y_lm(q)*f_l(q)) 
        DO ih = 1, nh(nt)
           IF (nb == indv(ih,nt)) THEN
              !l = nhtol (ih,nt)
              lm = nhtolm(ih,nt)
              DO ig = 1, npw_
                 vkb1(ig,ih) = ylm(ig,lm) * vq(ig)
              ENDDO
           ENDIF
        ENDDO
        !
     ENDDO
     !
     ! ... vkb1 contains all betas including angular part for type nt
     ! now add the structure factor and factor (-i)^l
     !
     DO na = 1, nat
        ! ordering: first all betas for atoms of type 1
        !           then  all betas for atoms of type 2  and so on
        IF (ityp(na) == nt) THEN
           arg = (q_(1) * tau(1,na) + &
                  q_(2) * tau(2,na) + &
                  q_(3) * tau(3,na) ) * tpi
           phase = CMPLX(COS(arg), - SIN(arg), KIND=DP)
           DO ig = 1, npw_
              sk (ig) = eigts1(mill(1,igk_(ig)), na) * &
                        eigts2(mill(2,igk_(ig)), na) * &
                        eigts3(mill(3,igk_(ig)), na)
           ENDDO
           !
           DO ih = 1, nh (nt)
              jkb = jkb + 1
              pref = (0.d0, -1.d0)**nhtol(ih,nt) * phase
              DO ig = 1, npw_
                 vkb_(ig,1) = vkb1(ig,ih) * sk(ig) * pref
              ENDDO
              !
              DO ig = npw_+1, npwx
                 vkb_(ig, 1) = (0.0_DP, 0.0_DP)
              ENDDO
              !
              CALL calbec( npw_, vkb_, evc_exx, becpsi_k(jkb:jkb,:), &
                           nibands(my_egrp_id+1) )
           ENDDO
           !
        ENDIF
        !
     ENDDO
  ENDDO
  !
  DEALLOCATE( gk )
  DEALLOCATE( ylm )
  DEALLOCATE( vq )
  DEALLOCATE( qg )
  DEALLOCATE( sk )
  DEALLOCATE( vkb1 )
  !
  RETURN
  !
  END SUBROUTINE compute_becpsi
  !
  !
  !-----------------------------------------------------------------------------
  SUBROUTINE aceinit( DoLoc, exex )
    !----------------------------------------------------------------------------
    !! ACE Initialization
    !
    USE wvfct,            ONLY : nbnd, npwx, current_k
    USE klist,            ONLY : nks, xk, ngk, igk_k
    USE uspp,             ONLY : nkb, vkb, okvan
    USE becmod,           ONLY : allocate_bec_type, deallocate_bec_type, &
                                 bec_type, calbec
    USE lsda_mod,         ONLY : current_spin, lsda, isk
    USE io_files,         ONLY : nwordwfc, iunwfc
    USE buffers,          ONLY : get_buffer
    USE mp_pools,         ONLY : inter_pool_comm
    USE mp_bands,         ONLY : intra_bgrp_comm
    USE mp,               ONLY : mp_sum
    USE wavefunctions,    ONLY : evc
    !
    IMPLICIT NONE
    !
    LOGICAL, INTENT(IN) :: DoLoc
    !! if TRUE calculates exact exchange with SCDM orbitals
    REAL(DP), OPTIONAL, INTENT(OUT) :: exex
    !! ACE energy
    !
    ! ... local variables
    !
    REAL(DP) :: ee, eexx
    INTEGER :: ik, npw
    TYPE(bec_type) :: becpsi
    !
    IF (nbndproj < x_nbnd_occ .OR. nbndproj > nbnd) THEN 
       WRITE( stdout, '(3(A,I4))' ) ' occ = ', x_nbnd_occ, ' proj = ', nbndproj, &
                                    ' tot = ', nbnd
       CALL errore( 'aceinit', 'n_proj must be between occ and tot.', 1 )
    ENDIF
    !
    IF (.NOT. ALLOCATED(xi)) ALLOCATE( xi(npwx*npol,nbndproj,nks) )
    IF ( okvan ) CALL allocate_bec_type( nkb, nbnd, becpsi )
    !
    eexx = 0.0d0
    xi = (0.0d0,0.0d0)
    !
    DO ik = 1, nks
       npw = ngk(ik)
       current_k = ik
       IF ( lsda ) current_spin = isk(ik)
       IF ( nks > 1 ) CALL get_buffer( evc, nwordwfc, iunwfc, ik )
       IF ( okvan ) THEN
          CALL init_us_2( npw, igk_k(1,ik), xk(:,ik), vkb )
          CALL calbec( npw, vkb, evc, becpsi, nbnd )
       ENDIF
       IF (gamma_only) THEN
          CALL aceinit_gamma( DoLoc, npw, nbnd, evc, xi(1,1,ik), becpsi, ee )
       ELSE
          CALL aceinit_k( DoLoc, npw, nbnd, evc, xi(1,1,ik), becpsi, ee )
       ENDIF
       eexx = eexx + ee
    ENDDO
    !
    CALL mp_sum( eexx, inter_pool_comm )
    ! WRITE(stdout,'(/,5X,"ACE energy",f15.8)') eexx
    !
    IF (PRESENT(exex)) exex = eexx
    IF ( okvan ) CALL deallocate_bec_type( becpsi )
    !
    domat = .FALSE.
    !
  END SUBROUTINE aceinit
  !
  !
  !---------------------------------------------------------------------------------
  SUBROUTINE aceinit_gamma( DoLoc, nnpw, nbnd, phi, xitmp, becpsi, exxe )
    !-------------------------------------------------------------------------------
    !! Compute xi(npw,nbndproj) for the ACE method.
    !
    USE becmod,         ONLY : bec_type
    USE lsda_mod,       ONLY : current_spin
    USE mp,             ONLY : mp_stop
    !
    IMPLICIT NONE
    !
    LOGICAL, INTENT(IN) :: DoLoc
    !! if TRUE calculates exact exchange with SCDM orbitals
    INTEGER :: nnpw
    !! number of pw
    INTEGER :: nbnd
    !! number of bands
    COMPLEX(DP) :: phi(nnpw,nbnd)
    !! wavefunction
    COMPLEX(DP) :: xitmp(nnpw,nbndproj)
    !! xi(npw,nbndproj)
    TYPE(bec_type), INTENT(IN) :: becpsi
    !! <beta|psi>
    REAL(DP) :: exxe
    !! exx energy
    !
    ! ... local variables
    !
    INTEGER :: nrxxs
    REAL(DP), ALLOCATABLE :: mexx(:,:)
    REAL(DP), PARAMETER :: Zero=0.0d0, One=1.0d0, Two=2.0d0, Pt5=0.50d0  
    LOGICAL :: domat0  
    !
    CALL start_clock( 'aceinit' )  
    !
    nrxxs = dfftt%nnr * npol  
    !
    ALLOCATE( mexx(nbndproj,nbndproj) )  
    xitmp = (Zero,Zero)  
    mexx = Zero  
    !  
    IF ( DoLoc ) then    
      CALL vexx_loc( nnpw, nbndproj, xitmp, mexx )
      CALL MatSymm( 'S', 'L', mexx,nbndproj )
    ELSE  
    ! |xi> = Vx[phi]|phi>
    CALL vexx( nnpw, nnpw, nbndproj, phi, xitmp, becpsi )
    ! mexx = <phi|Vx[phi]|phi>
    CALL matcalc( 'exact', .TRUE., 0, nnpw, nbndproj, nbndproj, phi, xitmp, mexx, exxe )
    ! |xi> = -One * Vx[phi]|phi> * rmexx^T
    ENDIF  
    !
    CALL aceupdate( nbndproj, nnpw, xitmp, mexx )
    !
    DEALLOCATE( mexx )  
    !
    IF ( local_thr > 0.0d0 ) THEN
      domat0 = domat
      domat = .TRUE.  
      CALL vexxace_gamma( nnpw, nbndproj, evc0(1,1,current_spin), exxe )  
      evc0(:,:,current_spin) = phi(:,:)  
      domat = domat0  
    ENDIF
    !
    CALL stop_clock( 'aceinit' )  
    !
  END SUBROUTINE aceinit_gamma
  !
  !
  !----------------------------------------------------------------------------------
  SUBROUTINE vexxace_gamma( nnpw, nbnd, phi, exxe, vphi )
    !-------------------------------------------------------------------------------
    !! Do the ACE potential and (optional) print the ACE matrix representation.
    !
    USE wvfct,        ONLY : current_k, wg
    USE lsda_mod,     ONLY : current_spin
    !
    IMPLICIT NONE
    !
    INTEGER :: nnpw
    !! number of plane waves
    INTEGER :: nbnd
    !! number of bands
    COMPLEX(DP) :: phi(nnpw,nbnd) 
    !! wave function
    REAL(DP) :: exxe
    !! exx energy
    COMPLEX(DP), OPTIONAL :: vphi(nnpw,nbnd)
    !! v times phi
    !
    ! ... local variables
    !
    INTEGER :: i, ik
    REAL*8, ALLOCATABLE :: rmexx(:,:)  
    COMPLEX(DP),ALLOCATABLE :: cmexx(:,:), vv(:,:)  
    REAL*8, PARAMETER :: Zero=0.0d0, One=1.0d0, Two=2.0d0, Pt5=0.50d0  
    !
    CALL start_clock( 'vexxace' )
    !
    ALLOCATE( vv(nnpw,nbnd) )
    !
    IF (PRESENT(vphi)) THEN  
      vv = vphi  
    ELSE  
      vv = (Zero,Zero)  
    ENDIF  
    !
    ! do the ACE potential
    ALLOCATE( rmexx(nbndproj,nbnd), cmexx(nbndproj,nbnd) )
    !
    rmexx = Zero  
    cmexx = (Zero,Zero)  
    ! <xi|phi>
    CALL matcalc( '<xi|phi>', .FALSE. , 0, nnpw, nbndproj, nbnd, xi(1,1,current_k), &
                  phi, rmexx, exxe )
    ! |vv> = |vphi> + (-One) * |xi> * <xi|phi>
    cmexx = (One,Zero)*rmexx
    !
    CALL ZGEMM( 'N', 'N', nnpw, nbnd, nbndproj, -(One,Zero), xi(1,1,current_k), &
                nnpw, cmexx, nbndproj, (One,Zero), vv, nnpw )
    !
    DEALLOCATE( cmexx, rmexx )  
    !
    IF (domat) THEN
      ALLOCATE( rmexx(nbnd,nbnd) )
      CALL matcalc( 'ACE', .TRUE., 0, nnpw, nbnd, nbnd, phi, vv, rmexx, exxe )
      DEALLOCATE( rmexx )
#if defined(__DEBUG)
        WRITE(stdout,'(4(A,I3),A,I9,A,f12.6)') 'vexxace: nbnd=', nbnd, ' nbndproj=',nbndproj, &
                                               ' k=',current_k,' spin=',current_spin,' npw=', &
                                               nnpw, ' E=',exxe
      ELSE
        WRITE(stdout,'(4(A,I3),A,I9)')         'vexxace: nbnd=', nbnd, ' nbndproj=',nbndproj, &
                                               ' k=',current_k,' spin=',current_spin,' npw=', &
                                               nnpw
#endif
      ENDIF
      !
      IF (PRESENT(vphi)) vphi = vv
      DEALLOCATE( vv )
      !
      CALL stop_clock( 'vexxace' )
      !
  END SUBROUTINE vexxace_gamma
  !
  !
  !-------------------------------------------------------------------------------------------
  SUBROUTINE aceupdate( nbndproj, nnpw, xitmp, rmexx )
    !----------------------------------------------------------------------------------------
    !! Build the ACE operator from the potential amd matrix (rmexx is assumed symmetric
    !! and only the Lower Triangular part is considered).
    !
    IMPLICIT NONE
    !
    INTEGER :: nbndproj
    !! number of bands
    INTEGER :: nnpw
    !! number of PW
    COMPLEX(DP) :: xitmp(nnpw,nbndproj)
    !! xi(nnpw,nbndproj)
    REAL(DP) :: rmexx(nbndproj,nbndproj)
    !! |xi> = -One * Vx[phi]|phi> * rmexx^T
    !
    ! ... local variables
    !
    COMPLEX(DP), ALLOCATABLE :: cmexx(:,:)
    REAL(DP), PARAMETER :: Zero=0.0d0, One=1.0d0, Two=2.0d0, Pt5=0.50d0
    !
    CALL start_clock( 'aceupdate' )
    !
    ! rmexx = -(Cholesky(rmexx))^-1
    rmexx = -rmexx
    ! CALL invchol( nbndproj, rmexx )
    CALL MatChol( nbndproj, rmexx )
    CALL MatInv( 'L', nbndproj, rmexx )
    !
    ! |xi> = -One * Vx[phi]|phi> * rmexx^T
    ALLOCATE( cmexx(nbndproj,nbndproj) )
    cmexx = (One,Zero)*rmexx
    CALL ZTRMM( 'R', 'L', 'C', 'N', nnpw, nbndproj, (One,Zero), cmexx, nbndproj, xitmp, nnpw )
    !
    DEALLOCATE( cmexx )
    !
    CALL stop_clock( 'aceupdate' )
    !
  END SUBROUTINE
  !
  !
  !---------------------------------------------------------------------------------------------
  SUBROUTINE aceinit_k( DoLoc, nnpw, nbnd, phi, xitmp, becpsi, exxe )
    !-----------------------------------------------------------------------------------------
    !! Compute xi(npw,nbndproj) for the ACE method.
    !
    USE becmod,               ONLY : bec_type
    USE wvfct,                ONLY : current_k, npwx
    USE klist,                ONLY : wk
    USE noncollin_module,     ONLY : npol
    !
    IMPLICIT NONE
    !
    LOGICAL, INTENT(IN) :: DoLoc
    !! if TRUE calculates exact exchange with SCDM orbitals
    INTEGER :: nnpw
    !! number of PW
    INTEGER :: nbnd
    !! number of bands
    COMPLEX(DP) :: phi(npwx*npol,nbnd)
    !! wave function
    COMPLEX(DP) :: xitmp(npwx*npol,nbndproj)
    !! xi(nnpw,nbndproj)
    TYPE(bec_type), INTENT(IN) :: becpsi
    !! <beta|psi>
    REAL(DP) :: exxe
    !! exx energy
    !
    ! ... local variables
    !
    COMPLEX(DP), ALLOCATABLE :: mexx(:,:)
    REAL(DP) :: exxe0
    REAL(DP), PARAMETER :: Zero=0.0d0, One=1.0d0, Two=2.0d0, Pt5=0.50d0
    INTEGER :: i
    LOGICAL :: domat0
    !
    CALL start_clock( 'aceinit' )
    !
    IF (nbndproj>nbnd) CALL errore( 'aceinit_k', 'nbndproj greater than nbnd.', 1 )
    IF (nbndproj<=0)   CALL errore( 'aceinit_k', 'nbndproj le 0.', 1 )
    !
    ALLOCATE( mexx(nbndproj,nbndproj) )
    xitmp = (Zero,Zero)
    mexx  = (Zero,Zero)
    IF ( DoLoc ) THEN
      CALL vexx_loc_k( nnpw, nbndproj, xitmp, mexx, exxe )
      CALL MatSymm_k( 'S', 'L', mexx, nbndproj )
    ELSE
      ! |xi> = Vx[phi]|phi>
      CALL vexx( npwx, nnpw, nbndproj, phi, xitmp, becpsi )
      ! mexx = <phi|Vx[phi]|phi>
      CALL matcalc_k( 'exact', .TRUE., 0, current_k, npwx*npol, nbndproj, nbndproj, &
                      phi, xitmp, mexx, exxe )
    ENDIF
#if defined(__DEBUG)
      WRITE( stdout,'(3(A,I3),A,I9,A,f12.6)') 'aceinit_k: nbnd=', nbnd, ' nbndproj=',nbndproj, &
                                              ' k=',current_k,' npw=',nnpw,' Ex(k)=',exxe
#endif
    ! Skip k-points that have exactly zero weight
    IF(wk(current_k)/=0._dp)THEN
      ! |xi> = -One * Vx[phi]|phi> * rmexx^T
      CALL aceupdate_k( nbndproj, nnpw, xitmp, mexx )
    ENDIF
    !
    DEALLOCATE( mexx )
    !
    IF ( DoLoc ) THEN
       domat0 = domat
       domat = .TRUE.
       CALL vexxace_k( nnpw, nbnd, evc0(1,1,current_k), exxe )
       evc0(:,:,current_k) = phi(:,:)
       domat = domat0
    ENDIF 
    !
    CALL stop_clock( 'aceinit' )
    !
  END SUBROUTINE aceinit_k
  !
  !
  !------------------------------------------------------------------------------
  SUBROUTINE aceupdate_k( nbndproj, nnpw, xitmp, mexx )
    !----------------------------------------------------------------------------
    !! Updates xi(npw,nbndproj) for the ACE method.
    !
    USE wvfct,                ONLY : npwx
    USE noncollin_module,     ONLY : noncolin, npol
    !
    IMPLICIT NONE
    !
    INTEGER :: nbndproj
    !! number of bands
    INTEGER :: nnpw
    !! number of PW
    COMPLEX(DP) :: mexx(nbndproj,nbndproj)
    !! mexx = -(Cholesky(mexx))^-1
    COMPLEX(DP) :: xitmp(npwx*npol,nbndproj)
    !! |xi> = -One * Vx[phi]|phi> * mexx^T
    !
    CALL start_clock( 'aceupdate' )
    !
    ! mexx = -(Cholesky(mexx))^-1
    mexx = -mexx
    CALL invchol_k( nbndproj, mexx )
    !
    ! |xi> = -One * Vx[phi]|phi> * mexx^T
    CALL ZTRMM( 'R', 'L', 'C', 'N', npwx*npol, nbndproj, (1.0_dp,0.0_dp), mexx,nbndproj, &
                xitmp, npwx*npol )
    !
    CALL stop_clock( 'aceupdate' )
    !
  END SUBROUTINE aceupdate_k
  !
  !
  !--------------------------------------------------------------------------------------
  SUBROUTINE vexxace_k( nnpw, nbnd, phi, exxe, vphi )
    !-----------------------------------------------------------------------------------
    !! Do the ACE potential and (optional) print the ACE matrix representation.
    !
    USE becmod,               ONLY : calbec
    USE wvfct,                ONLY : current_k, npwx
    USE noncollin_module,     ONLY : npol
    !
    IMPLICIT NONE
    !
    REAL(DP) :: exxe
    !! exx energy
    INTEGER :: nnpw
    !! number of PW
    INTEGER :: nbnd
    !! number of bands
    COMPLEX(DP) :: phi(npwx*npol,nbnd)
    !! wave function
    COMPLEX(DP), OPTIONAL :: vphi(npwx*npol,nbnd)
    !! ACE potential
    !
    ! ... local variables
    !
    INTEGER :: i
    COMPLEX(DP), ALLOCATABLE :: cmexx(:,:), vv(:,:)
    REAL*8, PARAMETER :: Zero=0.0d0, One=1.0d0, Two=2.0d0, Pt5=0.50d0
    !
    CALL start_clock( 'vexxace' )
    !
    ALLOCATE( vv(npwx*npol,nbnd) )  
    IF (PRESENT(vphi)) THEN  
      vv = vphi  
    ELSE  
      vv = (Zero, Zero)
    ENDIF  
    !
    ! do the ACE potential! 
    ALLOCATE( cmexx(nbndproj,nbnd) )  
    cmexx = (Zero,Zero)  
    ! <xi|phi>
    CALL matcalc_k( '<xi|phi>', .FALSE., 0, current_k, npwx*npol, nbndproj, nbnd, &
                    xi(1,1,current_k), phi, cmexx, exxe )
    !
    ! |vv> = |vphi> + (-One) * |xi> * <xi|phi>! 
    CALL ZGEMM( 'N', 'N', npwx*npol, nbnd, nbndproj, -(One,Zero), xi(1,1,current_k), &
                npwx*npol, cmexx, nbndproj, (One,Zero), vv, npwx*npol )
    !
    IF (domat) THEN
       !
       IF ( nbndproj /= nbnd) THEN
          DEALLOCATE( cmexx )
          ALLOCATE( cmexx(nbnd,nbnd) )
       ENDIF
       !
       CALL matcalc_k( 'ACE', .TRUE., 0, current_k, npwx*npol, nbnd, nbnd, phi, vv, cmexx, exxe )
       !
#if defined(__DEBUG)
       WRITE(stdout,'(3(A,I3),A,I9,A,f12.6)') 'vexxace_k: nbnd=', nbnd, ' nbndproj=',nbndproj, &
                    ' k=',current_k,' npw=',nnpw, ' Ex(k)=',exxe
    ELSE
       WRITE(stdout,'(3(A,I3),A,I9)') 'vexxace_k: nbnd=', nbnd, ' nbndproj=',nbndproj, &
                       ' k=',current_k,' npw=',nnpw
#endif
    ENDIF
    !  
    IF (PRESENT(vphi)) vphi = vv
    DEALLOCATE( vv, cmexx )
    !
    CALL stop_clock( 'vexxace' )
    !
  END SUBROUTINE vexxace_k
  !
  !
  !---------------------------------------------------------------------------------
  SUBROUTINE vexx_loc( npw, nbnd, hpsi, mexx )
    !---------------------------------------------------------------------------------
    !! Exact exchange with SCDM orbitals.  
    !! Vx|phi> =  Vx|psi> <psi|Vx|psi>^(-1) <psi|Vx|phi>.  
    !! locmat contains localization integrals.
    !
    USE noncollin_module,  ONLY : npol
    USE cell_base,         ONLY : omega, alat
    USE wvfct,             ONLY : current_k
    USE klist,             ONLY : xk, nks, nkstot
    USE fft_interfaces,    ONLY : fwfft, invfft
    USE mp,                ONLY : mp_stop, mp_barrier, mp_sum
    USE mp_bands,          ONLY : intra_bgrp_comm, me_bgrp, nproc_bgrp
    USE exx_base,          ONLY : nqs, xkq_collect, index_xkq, index_xk, &
                                  g2_convolution
    !
    IMPLICIT NONE
    !
    INTEGER :: npw
    !! number of PW
    INTEGER :: nbnd
    !! number of bands
    COMPLEX(DP) :: hpsi(npw,nbnd)
    !! hpsi
    REAL(DP) :: mexx(nbnd,nbnd)
    !! mexx contains in output the exchange matrix
    !
    ! ... local variables
    !
    INTEGER :: nrxxs, npairs, ntot, NBands   
    INTEGER :: ig, ir, ik, ikq, iq, ibnd, jbnd, kbnd, NQR  
    INTEGER :: current_ik  
    REAL(DP) :: exxe  
    COMPLEX(DP), ALLOCATABLE :: rhoc(:), vc(:), RESULT(:,:)   
    REAL(DP), ALLOCATABLE :: fac(:)  
    REAL(DP) :: xkp(3), xkq(3)  
    INTEGER, EXTERNAL  :: global_kpoint_index  
    !
    WRITE( stdout, '(5X,A)' ) ' '   
    WRITE( stdout, '(5X,A)' ) 'Exact-exchange with localized orbitals'  
    !
    CALL start_clock( 'vexxloc' )
    !
    WRITE( stdout,'(7X,A,f24.12)' ) 'local_thr =', local_thr  
    nrxxs = dfftt%nnr  
    NQR = nrxxs*npol  
    !
    ! ... exchange projected onto localized orbitals 
    WRITE( stdout,'(A)' ) 'Allocating exx quantities...'
    ALLOCATE( fac(dfftt%ngm) )
    ALLOCATE( rhoc(nrxxs), vc(NQR) )
    ALLOCATE( RESULT(nrxxs,nbnd) ) 
    WRITE( stdout,'(A)' ) 'Allocations done.'
    !
    current_ik = global_kpoint_index( nkstot, current_k )
    xkp = xk(:,current_k)
    !
    vc = (0.0d0, 0.0d0)
    npairs = 0 
    !
    DO iq = 1, nqs
       ikq  = index_xkq(current_ik,iq)  
       ik   = index_xk(ikq)  
       xkq  = xkq_collect(:,ikq)  
       !  
       CALL g2_convolution( dfftt%ngm, gt, xkp, xkq, fac )  
       !  
       RESULT = (0.0d0, 0.0d0)  
       !  
       DO ibnd = 1, nbnd  
         !
         IF (x_occupation(ibnd,ikq) > 0.0d0) THEN
           !
           DO ir = 1, NQR   
             rhoc(ir) = locbuff(ir,ibnd,ikq) * locbuff(ir,ibnd,ikq) / omega  
           ENDDO
           !
           CALL fwfft( 'Rho', rhoc, dfftt )
           !
           vc = (0.0d0, 0.0d0)  
           DO ig = 1, dfftt%ngm  
               vc(dfftt%nl(ig))  = fac(ig) * rhoc(dfftt%nl(ig))   
               vc(dfftt%nlm(ig)) = fac(ig) * rhoc(dfftt%nlm(ig))  
           ENDDO  
           !
           CALL invfft( 'Rho', vc, dfftt )
           !
           DO ir = 1, NQR   
             RESULT(ir,ibnd) = RESULT(ir,ibnd) + locbuff(ir,ibnd,ikq) * vc(ir)   
           ENDDO  
           !
         ENDIF   
         !
         DO kbnd = 1, ibnd-1  
           IF ( (locmat(ibnd,kbnd,ikq) > local_thr) .AND. &  
                ( (x_occupation(ibnd,ikq) > 0.0d0) .OR.   &
                  (x_occupation(kbnd,ikq) > 0.0d0) ) ) THEN
             !
             !write(stdout,'(3I4,3f12.6,A)') ikq, ibnd, kbnd, x_occupation(ibnd,ikq), &
             !                    x_occupation(kbnd,ikq), locmat(ibnd,kbnd,ikq), ' IN '
             !
             DO ir = 1, NQR   
               rhoc(ir) = locbuff(ir,ibnd,ikq) * locbuff(ir,kbnd,ikq) / omega  
             ENDDO
             !
             npairs = npairs + 1  
             !
             CALL fwfft( 'Rho', rhoc, dfftt )
             !
             vc = (0.0d0, 0.0d0)
             !
             DO ig = 1, dfftt%ngm  
                 vc(dfftt%nl(ig))  = fac(ig) * rhoc(dfftt%nl(ig))   
                 vc(dfftt%nlm(ig)) = fac(ig) * rhoc(dfftt%nlm(ig))  
             ENDDO
             !
             CALL invfft( 'Rho', vc, dfftt )
             !
             DO ir = 1, NQR   
               RESULT(ir,kbnd) = RESULT(ir,kbnd) + x_occupation(ibnd,ikq) * locbuff(ir,ibnd,ikq) * vc(ir)   
             ENDDO
             !
             DO ir = 1, NQR   
               RESULT(ir,ibnd) = RESULT(ir,ibnd) + x_occupation(kbnd,ikq) * locbuff(ir,kbnd,ikq) * vc(ir)   
             ENDDO
             ! ELSE   
             !   write(stdout,'(3I4,3f12.6,A)') ikq, ibnd, kbnd, x_occupation(ibnd,ikq), &
             !               x_occupation(kbnd,ikq), locmat(ibnd,kbnd,ikq), '      OUT '  
           ENDIF
           !
         ENDDO
         !
       ENDDO   
       !
       DO jbnd = 1, nbnd  
         !
         CALL fwfft( 'Wave', RESULT(:,jbnd), dfftt )
         !
         DO ig = 1, npw  
            hpsi(ig,jbnd) = hpsi(ig,jbnd) - exxalfa*RESULT(dfftt%nl(ig),jbnd)   
         ENDDO
         !
       ENDDO
       !
    ENDDO
    !
    DEALLOCATE( fac, vc )
    DEALLOCATE( RESULT )
    !
    ! ... localized functions to G-space and exchange matrix onto localized functions
    ALLOCATE( RESULT(npw,nbnd) )
    RESULT = (0.0d0,0.0d0)
    !
    DO jbnd = 1, nbnd
      rhoc(:) = DBLE(locbuff(:,jbnd,ikq)) + (0.0d0,1.0d0)*0.0d0
      CALL fwfft( 'Wave' , rhoc, dfftt )
      DO ig = 1, npw
        RESULT(ig,jbnd) = rhoc(dfftt%nl(ig))
      ENDDO
    ENDDO
    !
    DEALLOCATE( rhoc )
    !
    CALL matcalc( 'M1-', .TRUE., 0, npw, nbnd, nbnd, RESULT, hpsi, mexx, exxe )
    !
    DEALLOCATE( RESULT )
    !
    NBands = INT(SUM(x_occupation(:,ikq)))
    ntot = NBands * (NBands-1)/2 + NBands * (nbnd-NBands)
    WRITE( stdout,'(7X,2(A,I12),A,f12.2)') '  Pairs(full): ',      ntot, &
                  '   Pairs(included): ', npairs, &
                  '   Pairs(%): ', DBLE(npairs)/DBLE(ntot)*100.0d0
    !
    CALL stop_clock( 'vexxloc' )
    !
  END SUBROUTINE vexx_loc
  !
  !
  !----------------------------------------------------------------------------------------------
  SUBROUTINE compute_density( DoPrint, Shift, CenterPBC, SpreadPBC, Overlap, PsiI, PsiJ, NQR, &
                              ibnd, jbnd )
    !-------------------------------------------------------------------------------------------
    !! Manipulate density: get pair density, center, spread, absolute overlap.  
    !! Shift:  
    !! .FALSE. refer the centers to the cell -L/2 ... +L/2;  
    !! .TRUE.  shift the centers to the cell 0 ... L (presumably the one given in input).
    !
    USE constants,        ONLY : pi, bohr_radius_angs 
    USE cell_base,        ONLY : alat, omega
    USE mp,               ONLY : mp_sum
    USE mp_bands,         ONLY : intra_bgrp_comm
    USE fft_types,        ONLY : fft_index_to_3d
    !
    IMPLICIT NONE
    !
    LOGICAL, INTENT(IN) :: DoPrint
    !! whether to print or not the quantities
    LOGICAL,  INTENT(IN) :: Shift
    !! .FALSE. Centers with respect to the minimum image cell convention;  
    !! .TRUE.  Centers shifted to the input cell.
    REAL(DP), INTENT(OUT) :: CenterPBC(3)
    REAL(DP), INTENT(OUT) :: SpreadPBC(3)
    REAL(DP), INTENT(OUT) :: Overlap
    INTEGER,  INTENT(IN) :: NQR
    REAL(DP), INTENT(IN) :: PsiI(NQR)
    REAL(DP), INTENT(IN) :: PsiJ(NQR) 
    INTEGER, INTENT(IN) :: ibnd
    INTEGER, INTENT(IN) :: jbnd
    !
    ! ... local variables
    !
    REAL(DP) :: vol, rbuff, TotSpread
    INTEGER :: ir, i, j, k
    LOGICAL :: offrange
    COMPLEX(DP) :: cbuff(3)
    REAL(DP), PARAMETER :: Zero=0.0d0, One=1.0d0, Two=2.0d0 
    !
    vol = omega / DBLE(dfftt%nr1 * dfftt%nr2 * dfftt%nr3)
    !
    CenterPBC = Zero
    SpreadPBC = Zero
    Overlap = Zero
    cbuff = (Zero,Zero)
    rbuff = Zero
    !
    DO ir = 1, dfftt%nr1x*dfftt%my_nr2p*dfftt%my_nr3p
       !
       ! ... three dimensional indexes
       !
       CALL fft_index_to_3d (ir, dfftt, i,j,k, offrange)
       IF ( offrange ) CYCLE
       !
       rbuff = PsiI(ir) * PsiJ(ir) / omega
       Overlap = Overlap + ABS(rbuff)*vol
       cbuff(1) = cbuff(1) + rbuff*EXP((Zero,One)*Two*pi*DBLE(i)/DBLE(dfftt%nr1))*vol 
       cbuff(2) = cbuff(2) + rbuff*EXP((Zero,One)*Two*pi*DBLE(j)/DBLE(dfftt%nr2))*vol
       cbuff(3) = cbuff(3) + rbuff*EXP((Zero,One)*Two*pi*DBLE(k)/DBLE(dfftt%nr3))*vol
    ENDDO
    !  
    CALL mp_sum( cbuff, intra_bgrp_comm )  
    CALL mp_sum( Overlap, intra_bgrp_comm )  
    !  
    CenterPBC(1) =  alat/Two/pi*AIMAG( LOG(cbuff(1)) )  
    CenterPBC(2) =  alat/Two/pi*AIMAG( LOG(cbuff(2)) )  
    CenterPBC(3) =  alat/Two/pi*AIMAG( LOG(cbuff(3)) )  
    !  
    IF (Shift) THEN  
      IF (CenterPBC(1) < Zero) CenterPBC(1) = CenterPBC(1) + alat  
      IF (CenterPBC(2) < Zero) CenterPBC(2) = CenterPBC(2) + alat  
      IF (CenterPBC(3) < Zero) CenterPBC(3) = CenterPBC(3) + alat  
    ENDIF  
    !  
    rbuff = DBLE(cbuff(1))**2 + AIMAG(cbuff(1))**2   
    SpreadPBC(1) = -(alat/Two/pi)**2 * DLOG(rbuff)   
    rbuff = DBLE(cbuff(2))**2 + AIMAG(cbuff(2))**2   
    SpreadPBC(2) = -(alat/Two/pi)**2 * DLOG(rbuff)   
    rbuff = DBLE(cbuff(3))**2 + AIMAG(cbuff(3))**2   
    SpreadPBC(3) = -(alat/Two/pi)**2 * DLOG(rbuff)   
    TotSpread = (SpreadPBC(1) + SpreadPBC(2) + SpreadPBC(3))*bohr_radius_angs**2    
    !  
    IF (DoPrint) THEN
      WRITE(stdout,'(A,2I4)')     'MOs:                  ', ibnd, jbnd
      WRITE(stdout,'(A,10f12.6)') 'Absolute Overlap:     ', Overlap
      WRITE(stdout,'(A,10f12.6)') 'Center(PBC)[A]:       ', CenterPBC(1)*bohr_radius_angs, &
              CenterPBC(2)*bohr_radius_angs, CenterPBC(3)*bohr_radius_angs
      WRITE(stdout,'(A,10f12.6)') 'Spread [A**2]:        ', SpreadPBC(1)*bohr_radius_angs**2, &
              SpreadPBC(2)*bohr_radius_angs**2, SpreadPBC(3)*bohr_radius_angs**2
      WRITE(stdout,'(A,10f12.6)') 'Total Spread [A**2]:  ', TotSpread
    ENDIF  
    !  
    IF (TotSpread < Zero) CALL errore( 'compute_density', 'Negative spread found', 1 )  
    !  
  END SUBROUTINE compute_density 
  !
  !
  !------------------------------------------------------------------------------------
  SUBROUTINE compute_density_k( DoPrint, Shift, CenterPBC, SpreadPBC, Overlap, PsiI, &
                                PsiJ, NQR, ibnd, jbnd )
    !----------------------------------------------------------------------------------
    !! Manipulate density: get pair density, center, spread, absolute overlap.  
    !! Shift:  
    !! .FALSE. refer the centers to the cell -L/2 ... +L/2   
    !! .TRUE.  shift the centers to the cell 0 ... L (presumably the 
    !!         one given in input )
    !
    USE constants,        ONLY : pi, bohr_radius_angs 
    USE cell_base,        ONLY : alat, omega
    USE mp,               ONLY : mp_sum
    USE mp_bands,         ONLY : intra_bgrp_comm
    USE fft_types,        ONLY : fft_index_to_3d
    !
    IMPLICIT NONE
    !
    LOGICAL, INTENT(IN) :: DoPrint
    !! whether to print or not the quantities 
    LOGICAL, INTENT(IN) :: Shift
    !! .FALSE. Centers with respect to the minimum image cell convention
    !! .TRUE.  Centers shifted to the input cell
    REAL(DP), INTENT(OUT) :: CenterPBC(3)
    !! Coordinates of the center
    REAL(DP), INTENT(OUT) :: SpreadPBC(3)
    REAL(DP), INTENT(OUT) :: Overlap
    INTEGER,  INTENT(IN) :: NQR
    COMPLEX(DP), INTENT(IN) :: PsiI(NQR)
    COMPLEX(DP), INTENT(IN) :: PsiJ(NQR) 
    INTEGER,  INTENT(IN) :: ibnd
    INTEGER,  INTENT(IN) :: jbnd
    !
    ! ... local variables
    !
    REAL(DP) :: vol, TotSpread, rbuff
    INTEGER :: ir, i, j, k
    LOGICAL :: offrange
    COMPLEX(DP) :: cbuff(3)
    REAL(DP), PARAMETER :: Zero=0.0d0, One=1.0d0, Two=2.0d0 
    !
    vol = omega / DBLE(dfftt%nr1 * dfftt%nr2 * dfftt%nr3)
    !
    CenterPBC = Zero 
    SpreadPBC = Zero 
    Overlap = Zero
    cbuff = (Zero, Zero) 
    rbuff = Zero
    !
    DO ir = 1, dfftt%nr1x*dfftt%my_nr2p*dfftt%my_nr3p
       !
       ! ... three dimensional indexes
       !
       CALL fft_index_to_3d (ir, dfftt, i,j,k, offrange)
       IF ( offrange ) CYCLE
       !
       rbuff = ABS(PsiI(ir) * CONJG(PsiJ(ir)) / omega )
       Overlap = Overlap + rbuff*vol
       cbuff(1) = cbuff(1) + rbuff*EXP((Zero,One)*Two*pi*DBLE(i)/DBLE(dfftt%nr1))*vol 
       cbuff(2) = cbuff(2) + rbuff*EXP((Zero,One)*Two*pi*DBLE(j)/DBLE(dfftt%nr2))*vol
       cbuff(3) = cbuff(3) + rbuff*EXP((Zero,One)*Two*pi*DBLE(k)/DBLE(dfftt%nr3))*vol
    ENDDO
    !
    CALL mp_sum( cbuff, intra_bgrp_comm )
    CALL mp_sum( Overlap, intra_bgrp_comm )
    !
    CenterPBC(1) =  alat/Two/pi*AIMAG(LOG(cbuff(1)))
    CenterPBC(2) =  alat/Two/pi*AIMAG(LOG(cbuff(2)))
    CenterPBC(3) =  alat/Two/pi*AIMAG(LOG(cbuff(3)))
    !
    IF (Shift) THEN
      IF (CenterPBC(1) < Zero) CenterPBC(1) = CenterPBC(1) + alat
      IF (CenterPBC(2) < Zero) CenterPBC(2) = CenterPBC(2) + alat
      IF (CenterPBC(3) < Zero) CenterPBC(3) = CenterPBC(3) + alat
    ENDIF
    !
    rbuff = DBLE(cbuff(1))**2 + AIMAG(cbuff(1))**2 
    SpreadPBC(1) = -(alat/Two/pi)**2 * LOG(rbuff) 
    rbuff = DBLE(cbuff(2))**2 + AIMAG(cbuff(2))**2 
    SpreadPBC(2) = -(alat/Two/pi)**2 * LOG(rbuff) 
    rbuff = DBLE(cbuff(3))**2 + AIMAG(cbuff(3))**2 
    SpreadPBC(3) = -(alat/Two/pi)**2 * LOG(rbuff) 
    TotSpread = (SpreadPBC(1) + SpreadPBC(2) + SpreadPBC(3))*bohr_radius_angs**2  
    !
    IF (DoPrint) THEN
      WRITE(stdout,'(A,2I4)')     'MOs:                  ', ibnd, jbnd
      WRITE(stdout,'(A,10f12.6)') 'Absolute Overlap:     ', Overlap 
      WRITE(stdout,'(A,10f12.6)') 'Center(PBC)[A]:       ', CenterPBC(1)*bohr_radius_angs, &
              CenterPBC(2)*bohr_radius_angs, CenterPBC(3)*bohr_radius_angs
      WRITE(stdout,'(A,10f12.6)') 'Spread [A**2]:        ', SpreadPBC(1)*bohr_radius_angs**2, &
              SpreadPBC(2)*bohr_radius_angs**2, SpreadPBC(3)*bohr_radius_angs**2
      WRITE(stdout,'(A,10f12.6)') 'Total Spread [A**2]:  ', TotSpread
    ENDIF 
    !
    IF (TotSpread < Zero) CALL errore( 'compute_density_k', 'Negative spread found', 1 )
    !
  END SUBROUTINE compute_density_k
  !
  !
  !------------------------------------------------------------------------
  SUBROUTINE vexx_loc_k( npw, NBands, hpsi, mexx, exxe )
    !-----------------------------------------------------------------------
    !! Generic, k-point version of \(\texttt{vexx}\).
    !
    USE cell_base,       ONLY : omega
    USE gvect,           ONLY : ngm, g
    USE wvfct,           ONLY : current_k, npwx
    USE klist,           ONLY : xk, nks, nkstot
    USE fft_interfaces,  ONLY : fwfft, invfft
    USE exx_base,        ONLY : index_xkq, nqs, index_xk, xkq_collect, &
                                g2_convolution
    USE exx_band,        ONLY : igk_exx 
    !
    IMPLICIT NONE
    !
    INTEGER :: npw
    !! number of PW
    INTEGER :: NBands 
    !! number of bands
    COMPLEX(DP) :: hpsi(npwx*npol,NBands)
    !! h psi
    COMPLEX(DP) :: mexx(NBands,NBands)
    !! mexx contains in output the exchange matrix
    REAL(DP) :: exxe
    !! exx energy
    !
    ! ... local variables
    !
    COMPLEX(DP), ALLOCATABLE :: RESULT(:), RESULT2(:,:)
    COMPLEX(DP), ALLOCATABLE :: rhoc(:), vc(:)
    REAL(DP), ALLOCATABLE :: fac(:)
    INTEGER :: ibnd, jbnd, ik, ikq, iq
    INTEGER :: ir, ig, NBin, NBtot
    INTEGER :: current_ik, current_jk
    INTEGER :: nrxxs
    REAL(DP) :: xkp(3)
    REAL(DP) :: xkq(3)
    !
    INTEGER, EXTERNAL :: global_kpoint_index
    !
    CALL start_clock( 'vexxloc' )
    !
    ALLOCATE( fac(dfftt%ngm) )
    !
    nrxxs = dfftt%nnr
    !
    ALLOCATE( RESULT(nrxxs) )
    ALLOCATE( rhoc(nrxxs), vc(nrxxs) )
    !
    current_ik = global_kpoint_index ( nkstot, current_k )
    current_jk = index_xkq(current_ik,1)
    !
    NBin = 0  
    NBtot  = 0
    xkp = xk(:,current_k)
    DO jbnd = 1, NBands 
       RESULT = (0.0_DP, 0.0_DP) 
       DO iq = 1, nqs
          ikq = index_xkq(current_ik,iq)
          ik  = index_xk(ikq)
          xkq = xkq_collect(:,ikq)
          CALL g2_convolution( dfftt%ngm, gt, xkp, xkq, fac )
          DO ibnd = 1, NBands 
             ! IF ( abs(x_occupation(ibnd,ik)) < eps_occ) CYCLE 
             ! 
             NBtot = NBtot + 1 
             IF ((exxmat(ibnd,ikq,jbnd,current_k) > local_thr).AND. &
                ((x_occupation(ibnd,ik) > eps_occ))) THEN 
                  NBin = NBin + 1
               !
               ! write(stdout,'(4I4,f12.6,A)') ibnd, ikq, jbnd, current_k, exxmat(ibnd,ikq,jbnd,current_k), ' IN '
!$omp parallel do default(shared), private(ir)
               DO ir = 1, nrxxs
                  rhoc(ir)=CONJG(exxbuff(ir,ibnd,ikq))*exxbuff(ir,jbnd,current_jk) / omega
               ENDDO
!$omp end parallel do
               CALL fwfft( 'Rho', rhoc, dfftt )
               vc = (0._DP, 0._DP)
!$omp parallel do default(shared), private(ig)
               DO ig = 1, dfftt%ngm  
                  vc(dfftt%nl(ig)) = & 
                        fac(ig) * rhoc(dfftt%nl(ig)) * x_occupation(ibnd,ik) / nqs
               ENDDO
!$omp end parallel do
               CALL invfft( 'Rho', vc, dfftt )
!$omp parallel do default(shared), private(ir)
               DO ir = 1, nrxxs
                  RESULT(ir) = RESULT(ir) + vc(ir)*exxbuff(ir,ibnd,ikq)
               ENDDO
!$omp end parallel do
!            ELSE
!              write(stdout,'(4I4,f12.6,A)') ibnd, ikq, jbnd, current_k, exxmat(ibnd,ikq,jbnd,current_k), ' OUT'
             ENDIF 
         ENDDO
       ENDDO 
       !
       CALL fwfft( 'Wave', RESULT, dfftt )
!$omp parallel do default(shared), private(ig)
       DO ig = 1, npw
          hpsi(ig,jbnd) = hpsi(ig,jbnd) - exxalfa*RESULT(dfftt%nl(igk_exx(ig,current_k)))
       ENDDO
!$omp end parallel do
    ENDDO 
    !
    DEALLOCATE( RESULT )
    DEALLOCATE( vc, fac )
    !
    ! ... Localized functions to G-space and exchange matrix onto localized functions
    ALLOCATE( RESULT2(npwx,NBands) )
    RESULT2 = (0.0d0,0.0d0)
    !
    DO jbnd = 1, NBands
      rhoc(:) = exxbuff(:,jbnd,current_jk)
      CALL fwfft( 'Wave' , rhoc, dfftt )
      DO ig = 1, npw
        RESULT2(ig,jbnd) = rhoc(dfftt%nl(igk_exx(ig,current_k)))
      ENDDO
    ENDDO
    !
    DEALLOCATE( rhoc )
    CALL matcalc_k( 'M1-', .TRUE., 0, current_k, npwx*npol, NBands, NBands, RESULT2, hpsi, mexx, exxe )
    DEALLOCATE( RESULT2 )
    !
    WRITE(stdout,'(7X,2(A,I12),A,f12.2)') '  Pairs(full): ',  NBtot, &
            '   Pairs(included): ', NBin, &
            '   Pairs(%): ', DBLE(NBin)/DBLE(NBtot)*100.0d0
    !
    CALL stop_clock( 'vexxloc' )
    !
  END SUBROUTINE vexx_loc_k
  !
  !
END MODULE exx
