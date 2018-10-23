! Copyrigh(C) 2005-2018 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------
MODULE exx
  !--------------------------------------
  !
  ! Variables and subroutines for calculation of exact-exchange contribution
  ! Implements ACE: Lin Lin, J. Chem. Theory Comput. 2016, 12, 2242
  ! Contains code for band parallelization over pairs of bands: see T. Barnes,
  ! T. Kurth, P. Carrier, N. Wichmann, D. Prendergast, P.R.C. Kent, J. Deslippe
  ! Computer Physics Communications 2017, dx.doi.org/10.1016/j.cpc.2017.01.008
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
  SAVE
  !
  ! general purpose vars
  !
  ! exxalfa is the parameter multiplying the exact-exchange part
  REAL(DP):: exxalfa=0._dp
  ! x_occupation(nbnd,nkstot) the weight of
  ! auxiliary functions in the density matrix
  REAL(DP),    ALLOCATABLE :: x_occupation(:,:)
#if defined(__CUDA)
  REAL(DP),    ALLOCATABLE, DEVICE :: x_occupation_d(:,:)
#endif
  ! number of bands of auxiliary functions with
  ! at least some x_occupation > eps_occ
  INTEGER :: x_nbnd_occ
  REAL(DP), PARAMETER :: eps_occ  = 1.d-8
  ! 
  ! Buffers: temporary (complex) buffer for wfc storage
  COMPLEX(DP), ALLOCATABLE :: exxbuff(:,:,:)
#if defined(__CUDA)
  COMPLEX(DP), ALLOCATABLE, DEVICE :: exxbuff_d(:,:,:)
#endif
  ! temporary (real) buffer for wfc storage
  REAL(DP), ALLOCATABLE    :: locbuff(:,:,:)
  ! buffer for matrix of localization integrals
  REAL(DP), ALLOCATABLE    :: locmat(:,:,:)
  !

  LOGICAL :: use_ace        !  true: Use Lin Lin's ACE method
                            !  false: do not use ACE, use old algorithm instead
  COMPLEX(DP), ALLOCATABLE :: xi(:,:,:)  ! ACE projectors
  COMPLEX(DP), ALLOCATABLE :: evc0(:,:,:)! old wfc (G-space) needed to compute fock3
  INTEGER :: nbndproj
  LOGICAL :: domat
  REAL(DP)::  local_thr        ! threshold for Lin Lin's SCDM localized orbitals:
                               ! discard contribution to V_x if overlap between
                               ! localized orbitals is smaller than "local_thr"
  !
  ! energy related variables
  !
  REAL(DP) :: fock0 = 0.0_DP, & !   sum <old|Vx(old)|old>
              fock1 = 0.0_DP, & !   sum <new|vx(old)|new>
              fock2 = 0.0_DP, & !   sum <new|vx(new)|new>
              fock3 = 0.0_DP, & !   sum <old|vx(new)|old>
              dexx  = 0.0_DP    !   fock1  - 0.5*(fock2+fock0)
  !
  ! custom fft grid and related G-vectors
  !
  TYPE ( fft_type_descriptor ) :: dfftt
  LOGICAL :: exx_fft_initialized = .FALSE.
  ! G^2 in custom grid
  REAL(kind=DP), DIMENSION(:), POINTER :: ggt
  ! G-vectors in custom grid
  REAL(kind=DP), DIMENSION(:,:),POINTER :: gt
  ! gstart_t=2 if ggt(1)=0, =1 otherwise
  INTEGER :: gstart_t
  ! number of plane waves in custom grid (Gamma-only)
  INTEGER :: npwt
  ! Total number of G-vectors in custom grid
  INTEGER :: ngmt_g
  ! energy cutoff for custom grid
  REAL(DP)  :: ecutfock
  !
  ! starting and ending band index used in bgrp parallelization
  INTEGER :: ibnd_start = 0, ibnd_end = 0
  ! starting and ending buffer index used in bgrp parallelization
  integer :: ibnd_buff_start, ibnd_buff_end
  !
 CONTAINS
#define _CX(A)  CMPLX(A,0._dp,kind=DP)
#define _CY(A)  CMPLX(0._dp,-A,kind=DP)
  !
  !------------------------------------------------------------------------
  SUBROUTINE exx_fft_create ()

    USE gvecw,        ONLY : ecutwfc
    USE gvect,        ONLY : ecutrho, ngm, g, gg, gstart, mill
    USE cell_base,    ONLY : at, bg, tpiba2
    USE recvec_subs,  ONLY : ggen, ggens
    USE fft_base,     ONLY : smap
    USE fft_types,    ONLY : fft_type_init
    USE symm_base,    ONLY : fft_fact
    USE mp_exx,       ONLY : nproc_egrp, negrp, intra_egrp_comm
    USE mp_bands,     ONLY : nproc_bgrp, intra_bgrp_comm, nyfft
    !
    USE klist,        ONLY : nks, xk
    USE mp_pools,     ONLY : inter_pool_comm
    USE mp,           ONLY : mp_max, mp_sum
    !
    USE control_flags,ONLY : tqr
    USE realus,       ONLY : qpointlist, tabxx, tabp
    USE exx_band,     ONLY : smap_exx

    IMPLICIT NONE
    INTEGER :: ik, ngmt
    INTEGER, ALLOCATABLE :: ig_l2gt(:), millt(:,:)
    INTEGER, EXTERNAL :: n_plane_waves
    REAL(dp) :: gkcut, gcutmt
    LOGICAL :: lpara

    IF( exx_fft_initialized) RETURN
    !
    ! Initialise the custom grid that allows us to put the wavefunction
    ! onto the new (smaller) grid for \rho=\psi_{k+q}\psi^*_k and vice versa
    !
    ! gkcut is such that all |k+G|^2 < gkcut (in units of (2pi/a)^2)
    ! Note that with k-points, gkcut > ecutwfc/(2pi/a)^2
    ! gcutmt is such that |q+G|^2 < gcutmt
    !
    IF ( gamma_only ) THEN       
       gkcut = ecutwfc/tpiba2
       gcutmt = ecutfock / tpiba2
    ELSE
       !
       gkcut = 0.0_dp
       DO ik = 1,nks
          gkcut = MAX ( gkcut, sqrt( sum(xk(:,ik)**2) ) )
       ENDDO
       CALL mp_max( gkcut, inter_pool_comm )
       ! Alternatively, variable "qnorm" earlier computed in "exx_grid_init"
       ! could be used as follows:
       ! gkcut = ( sqrt(ecutwfc/tpiba2) + qnorm )**2
       gkcut = ( sqrt(ecutwfc/tpiba2) + gkcut )**2
       ! 
       ! The following instruction may be needed if ecutfock \simeq ecutwfc
       ! and guarantees that all k+G are included
       !
       gcutmt = max(ecutfock/tpiba2,gkcut)
       !
    ENDIF
    !
    ! ... set up fft descriptors, including parallel stuff: sticks, planes, etc.
    !
    IF( negrp == 1 ) THEN
       !
       ! ... no band parallelization: exx grid is a subgrid of general grid
       !
       lpara = ( nproc_bgrp > 1 )
       CALL fft_type_init( dfftt, smap, "rho", gamma_only, lpara, &
            intra_bgrp_comm, at, bg, gcutmt, gcutmt/gkcut, &
            fft_fact=fft_fact, nyfft=nyfft )
       CALL ggens( dfftt, gamma_only, at, g, gg, mill, gcutmt, ngmt, gt, ggt )
       gstart_t = gstart
       npwt = n_plane_waves (ecutwfc/tpiba2, nks, xk, gt, ngmt)
       ngmt_g = ngmt
       CALL mp_sum (ngmt_g, intra_bgrp_comm )
       !
    ELSE
       !
       WRITE(6,"(5X,'Exchange parallelized over bands (',i4,' band groups)')")&
            negrp
       lpara = ( nproc_egrp > 1 )
       CALL fft_type_init( dfftt, smap_exx, "rho", gamma_only, lpara, &
            intra_egrp_comm, at, bg, gcutmt, gcutmt/gkcut, &
            fft_fact=fft_fact, nyfft=nyfft )
       ngmt = dfftt%ngm
       ngmt_g = ngmt
       CALL mp_sum( ngmt_g, intra_egrp_comm )
       ALLOCATE ( gt(3,dfftt%ngm) )
       ALLOCATE ( ggt(dfftt%ngm) )
       ALLOCATE ( millt(3,dfftt%ngm) )
       ALLOCATE ( ig_l2gt(dfftt%ngm) )
       CALL ggen( dfftt, gamma_only, at, bg, gcutmt, ngmt_g, ngmt, &
            gt, ggt, millt, ig_l2gt, gstart_t )
       DEALLOCATE ( ig_l2gt )
       DEALLOCATE ( millt )
       npwt = n_plane_waves (ecutwfc/tpiba2, nks, xk, gt, ngmt)
       !
    END IF
    ! define clock labels (this enables the corresponding fft too)
    dfftt%rho_clock_label = 'fftc' ; dfftt%wave_clock_label = 'fftcw' 
    !
    WRITE( stdout, '(/5x,"EXX grid: ",i8," G-vectors", 5x, &
         &   "FFT dimensions: (",i4,",",i4,",",i4,")")') ngmt_g, &
         &   dfftt%nr1, dfftt%nr2, dfftt%nr3
    exx_fft_initialized = .true.
    !
    IF(tqr) THEN
       IF(ecutfock==ecutrho) THEN
          WRITE(stdout,'(5x,"Real-space augmentation: EXX grid -> DENSE grid")')
          tabxx => tabp
       ELSE
          WRITE(stdout,'(5x,"Real-space augmentation: initializing EXX grid")')
          CALL qpointlist(dfftt, tabxx)
       ENDIF
    ENDIF

    RETURN
    !------------------------------------------------------------------------
  END SUBROUTINE exx_fft_create
  !------------------------------------------------------------------------
  !
  SUBROUTINE exx_gvec_reinit( at_old )
    !
    USE cell_base,  ONLY : bg
    IMPLICIT NONE
    REAL(dp), INTENT(in) :: at_old(3,3)
    INTEGER :: ig
    REAL(dp) :: gx, gy, gz
    !
    ! ... rescale g-vectors
    !
    CALL cryst_to_cart(dfftt%ngm, gt, at_old, -1)
    CALL cryst_to_cart(dfftt%ngm, gt, bg,     +1)
    !
    DO ig = 1, dfftt%ngm
       gx = gt(1, ig)
       gy = gt(2, ig)
       gz = gt(3, ig)
       ggt(ig) = gx * gx + gy * gy + gz * gz
    END DO
    !
  END SUBROUTINE exx_gvec_reinit
  !------------------------------------------------------------------------
  SUBROUTINE deallocate_exx ()
    !------------------------------------------------------------------------
    !
    USE becmod, ONLY : deallocate_bec_type, is_allocated_bec_type, bec_type
    USE us_exx, ONLY : becxx
    USE exx_base, ONLY : xkq_collect, index_xkq, index_xk, index_sym, rir, &
         working_pool, exx_grid_initialized
    !
    IMPLICIT NONE
    INTEGER :: ikq
    !
    exx_grid_initialized = .false.
    IF ( allocated(index_xkq) ) DEALLOCATE(index_xkq)
    IF ( allocated(index_xk ) ) DEALLOCATE(index_xk )
    IF ( allocated(index_sym) ) DEALLOCATE(index_sym)
    IF ( allocated(rir)       ) DEALLOCATE(rir)
    IF ( allocated(x_occupation) ) DEALLOCATE(x_occupation)
#if defined (__CUDA)
    IF ( allocated(x_occupation_d) ) DEALLOCATE(x_occupation_d)
#endif
    IF ( allocated(xkq_collect ) ) DEALLOCATE(xkq_collect)
    IF ( allocated(exxbuff) ) DEALLOCATE(exxbuff)
#if defined(__CUDA)
    IF ( allocated(exxbuff_d) ) DEALLOCATE(exxbuff_d)
#endif
    IF ( allocated(locbuff) ) DEALLOCATE(locbuff)
    IF ( allocated(locmat) )  DEALLOCATE(locmat)
    IF ( allocated(xi) )      DEALLOCATE(xi)
    IF ( allocated(evc0) )    DEALLOCATE(evc0)
    !
    IF(allocated(becxx)) THEN
      DO ikq = 1, SIZE(becxx)
        IF(is_allocated_bec_type(becxx(ikq))) CALL deallocate_bec_type(becxx(ikq))
      ENDDO
      DEALLOCATE(becxx)
    ENDIF
    !
    IF ( allocated(working_pool) )  DEALLOCATE(working_pool)
    !
    exx_fft_initialized = .false.
    IF ( ASSOCIATED (gt)  )  DEALLOCATE(gt)
    IF ( ASSOCIATED (ggt) )  DEALLOCATE(ggt)
    !
    !------------------------------------------------------------------------
  END SUBROUTINE deallocate_exx
  !------------------------------------------------------------------------
  !
  !
  !------------------------------------------------------------------------
  SUBROUTINE exxinit(DoLoc)
  !------------------------------------------------------------------------

    ! This SUBROUTINE is run before the first H_psi() of each iteration.
    ! It saves the wavefunctions for the right density matrix, in real space
    !
    ! DoLoc = .true.  ...    Real Array locbuff(ir, nbnd, nkqs)
    !         .false. ... Complex Array exxbuff(ir, nbnd/2, nkqs)
    !
    USE wavefunctions, ONLY : evc
    USE io_files,             ONLY : nwordwfc, iunwfc_exx
    USE buffers,              ONLY : get_buffer
    USE wvfct,                ONLY : nbnd, npwx, wg, current_k
    USE klist,                ONLY : ngk, nks, nkstot, xk, wk, igk_k
    USE symm_base,            ONLY : nsym, s, sr
    USE mp_pools,             ONLY : npool, nproc_pool, me_pool, inter_pool_comm
    USE mp_exx,               ONLY : me_egrp, negrp, &
                                     init_index_over_band, my_egrp_id,  &
                                     inter_egrp_comm, intra_egrp_comm, &
                                     iexx_start, iexx_end, all_start, all_end
    USE mp,                   ONLY : mp_sum, mp_bcast
    USE funct,                ONLY : get_exx_fraction, start_exx,exx_is_active,&
                                     get_screening_parameter, get_gau_parameter
    USE scatter_mod,          ONLY : gather_grid, scatter_grid
    USE fft_interfaces,       ONLY : invfft
    USE uspp,                 ONLY : nkb, vkb, okvan
    USE us_exx,               ONLY : rotate_becxx
    USE paw_variables,        ONLY : okpaw
    USE paw_exx,              ONLY : PAW_init_fock_kernel
    USE mp_orthopools,        ONLY : intra_orthopool_comm
    USE exx_base,             ONLY : nkqs, xkq_collect, index_xk, index_sym, &
         exx_set_symm, rir, working_pool, exxdiv, erfc_scrlen, gau_scrlen, &
         exx_divergence
    USE exx_band,             ONLY : change_data_structure, nwordwfc_exx, &
         transform_evc_to_exx, igk_exx, evc_exx

    USE wavefunctions_gpum, ONLY : using_evc
    USE uspp_gpum,                 ONLY : using_vkb ! is this needed?
    !
    IMPLICIT NONE
    INTEGER :: ik,ibnd, i, j, k, ir, isym, ikq, ig
    INTEGER :: ibnd_loop_start
    INTEGER :: ipol, jpol
    REAL(dp), ALLOCATABLE   :: occ(:,:)
    COMPLEX(DP),ALLOCATABLE :: temppsic(:)
    COMPLEX(DP),ALLOCATABLE :: temppsic_nc(:,:), psic_nc(:,:)
    COMPLEX(DP),ALLOCATABLE :: psic_exx(:)
    INTEGER :: nxxs, nrxxs
#if defined(__MPI)
    COMPLEX(DP),ALLOCATABLE  :: temppsic_all(:),      psic_all(:)
    COMPLEX(DP), ALLOCATABLE :: temppsic_all_nc(:,:), psic_all_nc(:,:)
#endif
    COMPLEX(DP) :: d_spin(2,2,48)
    INTEGER :: npw, current_ik
    INTEGER, EXTERNAL :: global_kpoint_index
    INTEGER :: ibnd_start_new, ibnd_end_new, max_buff_bands_per_egrp
    INTEGER :: ibnd_exx, evc_offset
    LOGICAL :: DoLoc 
    !
    CALL start_clock ('exxinit')
    IF ( Doloc ) THEN
        WRITE(stdout,'(/,5X,"Using localization algorithm with threshold: ",&
                & D10.2)') local_thr
        IF(.NOT.gamma_only) CALL errore('exxinit','SCDM with K-points NYI',1)
        IF(okvan.OR.okpaw) CALL errore('exxinit','SCDM with USPP/PAW not implemented',1)
    END IF 
    IF ( use_ace ) &
        WRITE(stdout,'(/,5X,"Using ACE for calculation of exact exchange")') 
    !
    CALL transform_evc_to_exx(2)
    !
    !  prepare the symmetry matrices for the spin part
    !
    IF (noncolin) THEN
       DO isym=1,nsym
          CALL find_u(sr(:,:,isym), d_spin(:,:,isym))
       ENDDO
    ENDIF

    CALL exx_fft_create()

    ! Note that nxxs is not the same as nrxxs in parallel case
    nxxs = dfftt%nr1x *dfftt%nr2x *dfftt%nr3x
    nrxxs= dfftt%nnr
#if defined(__MPI)
    IF (noncolin) THEN
       ALLOCATE(psic_all_nc(nxxs,npol), temppsic_all_nc(nxxs,npol) )
    ELSEIF ( .not. gamma_only ) THEN
       ALLOCATE(psic_all(nxxs), temppsic_all(nxxs) )
    ENDIF
#endif
    IF (noncolin) THEN
       ALLOCATE(temppsic_nc(nrxxs, npol), psic_nc(nrxxs, npol))
    ELSEIF ( .not. gamma_only ) THEN
       ALLOCATE(temppsic(nrxxs))
    ENDIF
    ALLOCATE(psic_exx(nrxxs))
    !
    IF (.not.exx_is_active()) THEN
       !
       erfc_scrlen = get_screening_parameter()
       gau_scrlen = get_gau_parameter()
       exxdiv  = exx_divergence()
       exxalfa = get_exx_fraction()
       !
       CALL start_exx()
    ENDIF

    IF ( .not.gamma_only ) CALL exx_set_symm (dfftt%nr1, dfftt%nr2, dfftt%nr3, &
                                              dfftt%nr1x,dfftt%nr2x,dfftt%nr3x )
    ! set occupations of wavefunctions used in the calculation of exchange term
    IF( .NOT. ALLOCATED(x_occupation)) ALLOCATE( x_occupation(nbnd,nkstot) )
#if defined (__CUDA)
    IF( .NOT. ALLOCATED(x_occupation_d)) ALLOCATE( x_occupation_d(nbnd,nkstot) )
#endif
    ALLOCATE ( occ(nbnd,nks) )
    DO ik =1,nks
       IF(abs(wk(ik)) > eps_occ ) THEN
          occ(1:nbnd,ik) = wg (1:nbnd, ik) / wk(ik)
       ELSE
          occ(1:nbnd,ik) = 0._dp
       ENDIF
    ENDDO
    CALL poolcollect(nbnd, nks, occ, nkstot, x_occupation)
#if defined (__CUDA)
    x_occupation_d = x_occupation
#endif

    DEALLOCATE ( occ )

    ! find an upper bound to the number of bands with non zero occupation.
    ! Useful to distribute bands among band groups

    x_nbnd_occ = 0
    DO ik =1,nkstot
       DO ibnd = max(1,x_nbnd_occ), nbnd
          IF (abs(x_occupation(ibnd,ik)) > eps_occ ) x_nbnd_occ = ibnd
       ENDDO
    ENDDO

    ! FIXME: IF(nbndproj.eq.0) nbndproj = x_nbnd_occ doesn't work 
    IF(nbndproj.eq.0) nbndproj = nbnd

    CALL divide ( inter_egrp_comm, x_nbnd_occ, ibnd_start, ibnd_end )
    CALL init_index_over_band(inter_egrp_comm,nbnd,nbnd)

    !this will cause exxbuff to be calculated for every band
    ibnd_start_new = iexx_start
    ibnd_end_new = iexx_end

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
    IF( DoLoc) then 
      IF (.not. allocated(locbuff)) ALLOCATE( locbuff(nrxxs*npol, nbnd, nks))
      IF (.not. allocated(locmat))  ALLOCATE( locmat(nbnd, nbnd, nks))
      IF (.not. allocated(evc0)) then 
        ALLOCATE( evc0(npwx*npol,nbndproj,nks) )
        evc0 = (0.0d0,0.0d0)
      END IF
    ELSE
      IF (.not. allocated(exxbuff)) THEN
         IF (gamma_only) THEN
            ALLOCATE( exxbuff(nrxxs*npol, ibnd_buff_start:ibnd_buff_start+max_buff_bands_per_egrp-1, nks))
         ELSE
            ALLOCATE( exxbuff(nrxxs*npol, ibnd_buff_start:ibnd_buff_start+max_buff_bands_per_egrp-1, nkqs))
#if defined(__CUDA)
            ALLOCATE( exxbuff_d(nrxxs*npol, ibnd_buff_start:ibnd_buff_start+max_buff_bands_per_egrp-1, nkqs))
#endif
         END IF
      END IF
    END IF

    !assign buffer
    IF(DoLoc) then
!$omp parallel do collapse(3) default(shared) firstprivate(npol,nrxxs,nkqs,ibnd_buff_start,ibnd_buff_end) private(ir,ibnd,ikq,ipol)
      DO ikq=1,SIZE(locbuff,3) 
         DO ibnd=1, x_nbnd_occ 
            DO ir=1,nrxxs*npol
               locbuff(ir,ibnd,ikq)=0.0_DP
            ENDDO
         ENDDO
      ENDDO
    ELSE
!$omp parallel do collapse(3) default(shared) firstprivate(npol,nrxxs,nkqs,ibnd_buff_start,ibnd_buff_end) private(ir,ibnd,ikq,ipol)
      DO ikq=1,SIZE(exxbuff,3) 
         DO ibnd=ibnd_buff_start,ibnd_buff_end
            DO ir=1,nrxxs*npol
               exxbuff(ir,ibnd,ikq)=(0.0_DP,0.0_DP)
            ENDDO
         ENDDO
      ENDDO
      ! the above loops will replaced with the following line soon
      !CALL threaded_memset(exxbuff, 0.0_DP, nrxxs*npol*SIZE(exxbuff,2)*nkqs*2)
    END IF
    !
    !   This is parallelized over pools. Each pool computes only its k-points
    !
    KPOINTS_LOOP : &
    DO ik = 1, nks

       IF ( nks > 1 ) &
          CALL get_buffer(evc_exx, nwordwfc_exx, iunwfc_exx, ik)
       !
       ! ik         = index of k-point in this pool
       ! current_ik = index of k-point over all pools
       !
       current_ik = global_kpoint_index ( nkstot, ik )
       !
       IF_GAMMA_ONLY : &
       IF (gamma_only) THEN
          !
          IF(mod(iexx_start,2)==0) THEN
             ibnd_loop_start=iexx_start-1
          ELSE
             ibnd_loop_start=iexx_start
          ENDIF

          evc_offset = 0
          DO ibnd = ibnd_loop_start, iexx_end, 2
             !
             psic_exx(:) = ( 0._dp, 0._dp )
             !
             IF ( ibnd < iexx_end ) THEN
                IF ( ibnd == ibnd_loop_start .and. MOD(iexx_start,2) == 0 ) THEN
                   DO ig=1,npwt
                      psic_exx(dfftt%nl(ig))  = ( 0._dp, 1._dp )*evc_exx(ig,1)
                      psic_exx(dfftt%nlm(ig)) = ( 0._dp, 1._dp )*conjg(evc_exx(ig,1))
                   ENDDO
                   evc_offset = -1
                ELSE
                   DO ig=1,npwt
                      psic_exx(dfftt%nl(ig))  = evc_exx(ig,ibnd-ibnd_loop_start+evc_offset+1)  &
                           + ( 0._dp, 1._dp ) * evc_exx(ig,ibnd-ibnd_loop_start+evc_offset+2)
                      psic_exx(dfftt%nlm(ig)) = conjg( evc_exx(ig,ibnd-ibnd_loop_start+evc_offset+1) ) &
                           + ( 0._dp, 1._dp ) * conjg( evc_exx(ig,ibnd-ibnd_loop_start+evc_offset+2) )
                   ENDDO
                END IF
             ELSE
                DO ig=1,npwt
                   psic_exx(dfftt%nl (ig)) = evc_exx(ig,ibnd-ibnd_loop_start+evc_offset+1)
                   psic_exx(dfftt%nlm(ig)) = conjg( evc_exx(ig,ibnd-ibnd_loop_start+evc_offset+1) )
                ENDDO
             ENDIF

             CALL invfft ('Wave', psic_exx, dfftt)
             IF(DoLoc) then
               locbuff(1:nrxxs,ibnd-ibnd_loop_start+evc_offset+1,ik)=Dble(  psic_exx(1:nrxxs) )
               IF(ibnd-ibnd_loop_start+evc_offset+2.le.nbnd) &
                  locbuff(1:nrxxs,ibnd-ibnd_loop_start+evc_offset+2,ik)=Aimag( psic_exx(1:nrxxs) )
             ELSE
               exxbuff(1:nrxxs,(ibnd+1)/2,ik)=psic_exx(1:nrxxs)
             END IF
             
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
                DO ir=1,nrxxs
                   temppsic_nc(ir,1) = ( 0._dp, 0._dp )
                   temppsic_nc(ir,2) = ( 0._dp, 0._dp )
                ENDDO
!$omp parallel do default(shared) private(ig) firstprivate(npw,ik,ibnd_exx)
                DO ig=1,npw
                   temppsic_nc(dfftt%nl(igk_exx(ig,ik)),1) = evc_exx(ig,ibnd-iexx_start+1)
                ENDDO
!$omp end parallel do
                CALL invfft ('Wave', temppsic_nc(:,1), dfftt)
!$omp parallel do default(shared) private(ig) firstprivate(npw,ik,ibnd_exx,npwx)
                DO ig=1,npw
                   temppsic_nc(dfftt%nl(igk_exx(ig,ik)),2) = evc_exx(ig+npwx,ibnd-iexx_start+1)
                ENDDO
!$omp end parallel do
                CALL invfft ('Wave', temppsic_nc(:,2), dfftt)
             ELSE
!$omp parallel do default(shared) private(ir) firstprivate(nrxxs)
                DO ir=1,nrxxs
                   temppsic(ir) = ( 0._dp, 0._dp )
                ENDDO
!$omp parallel do default(shared) private(ig) firstprivate(npw,ik,ibnd_exx)
                DO ig=1,npw
                   temppsic(dfftt%nl(igk_exx(ig,ik))) = evc_exx(ig,ibnd-iexx_start+1)
                ENDDO
!$omp end parallel do
                CALL invfft ('Wave', temppsic, dfftt)
             ENDIF
             !
             DO ikq=1,nkqs
                !
                IF (index_xk(ikq) /= current_ik) CYCLE
                isym = abs(index_sym(ikq) )
                !
                IF (noncolin) THEN ! noncolinear
#if defined(__MPI)
                   DO ipol=1,npol
                      CALL gather_grid(dfftt, temppsic_nc(:,ipol), temppsic_all_nc(:,ipol))
                   ENDDO
                   IF ( me_egrp == 0 ) THEN
!$omp parallel do default(shared) private(ir) firstprivate(npol,nxxs)
                      DO ir=1,nxxs
                         !DIR$ UNROLL_AND_JAM (2)
                         DO ipol=1,npol
                            psic_all_nc(ir,ipol) = (0.0_DP, 0.0_DP)
                         ENDDO
                      ENDDO
!$omp end parallel do
!$omp parallel do default(shared) private(ir) firstprivate(npol,isym,nxxs) reduction(+:psic_all_nc)
                      DO ir=1,nxxs
                         !DIR$ UNROLL_AND_JAM (4)
                         DO ipol=1,npol
                            DO jpol=1,npol
                               psic_all_nc(ir,ipol) = psic_all_nc(ir,ipol) + &
                                  conjg(d_spin(jpol,ipol,isym)) * &
                                  temppsic_all_nc(rir(ir,isym),jpol)
                            ENDDO
                         ENDDO
                      ENDDO
!$omp end parallel do
                   ENDIF
                   DO ipol=1,npol
                      CALL scatter_grid(dfftt,psic_all_nc(:,ipol), psic_nc(:,ipol))
                   ENDDO
#else
!$omp parallel do default(shared) private(ir) firstprivate(npol,nxxs)
                   DO ir=1,nxxs
                      !DIR$ UNROLL_AND_JAM (2)
                      DO ipol=1,npol
                         psic_nc(ir,ipol) = (0._dp, 0._dp)
                      ENDDO
                   ENDDO
!$omp end parallel do
!$omp parallel do default(shared) private(ipol,jpol,ir) firstprivate(npol,isym,nxxs) reduction(+:psic_nc)
                   DO ir=1,nxxs
                      !DIR$ UNROLL_AND_JAM (4)
                      DO ipol=1,npol
                         DO jpol=1,npol
                            psic_nc(ir,ipol) = psic_nc(ir,ipol) + conjg(d_spin(jpol,ipol,isym))* temppsic_nc(rir(ir,isym),jpol)
                         ENDDO
                      ENDDO
                   ENDDO
!$omp end parallel do
#endif
                   IF (index_sym(ikq) > 0 ) THEN
                      ! sym. op. without time reversal: normal case
!$omp parallel do default(shared) private(ir) firstprivate(ibnd,isym,ikq)
                      DO ir=1,nrxxs
                         exxbuff(ir,ibnd,ikq)=psic_nc(ir,1)
                         exxbuff(ir+nrxxs,ibnd,ikq)=psic_nc(ir,2)
                      ENDDO
!$omp end parallel do
                   ELSE
                      ! sym. op. with time reversal: spin 1->2*, 2->-1*
!$omp parallel do default(shared) private(ir) firstprivate(ibnd,isym,ikq)
                      DO ir=1,nrxxs
                         exxbuff(ir,ibnd,ikq)=CONJG(psic_nc(ir,2))
                         exxbuff(ir+nrxxs,ibnd,ikq)=-CONJG(psic_nc(ir,1))
                      ENDDO
!$omp end parallel do
                   ENDIF
                ELSE ! noncolinear
#if defined(__MPI)
                   CALL gather_grid(dfftt,temppsic,temppsic_all)
                   IF ( me_egrp == 0 ) THEN
!$omp parallel do default(shared) private(ir) firstprivate(isym)
                      DO ir=1,nxxs
                         psic_all(ir) = temppsic_all(rir(ir,isym))
                      ENDDO
!$omp end parallel do
                   ENDIF
                   CALL scatter_grid(dfftt,psic_all,psic_exx)
#else
!$omp parallel do default(shared) private(ir) firstprivate(isym)
                   DO ir=1,nrxxs
                      psic_exx(ir) = temppsic(rir(ir,isym))
                   ENDDO
!$omp end parallel do
#endif
!$omp parallel do default(shared) private(ir) firstprivate(isym,ibnd,ikq)
                   DO ir=1,nrxxs
                      IF (index_sym(ikq) < 0 ) THEN
                         psic_exx(ir) = conjg(psic_exx(ir))
                      ENDIF
                      exxbuff(ir,ibnd,ikq)=psic_exx(ir)
                   ENDDO
!$omp end parallel do
                   !
                ENDIF ! noncolinear
                
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
    DEALLOCATE (psic_exx)
    IF (noncolin) THEN
       DEALLOCATE(temppsic_nc, psic_nc)
#if defined(__MPI)
       DEALLOCATE(temppsic_all_nc, psic_all_nc)
#endif
    ELSE IF ( .not. gamma_only ) THEN
       DEALLOCATE(temppsic)
#if defined(__MPI)
       DEALLOCATE(temppsic_all, psic_all)
#endif
    ENDIF
    !
    ! Each wavefunction in exxbuff is computed by a single pool, collect among 
    ! pools in a smart way (i.e. without doing all-to-all sum and bcast)
    ! See also the initialization of working_pool in exx_mp_init
    ! Note that in Gamma-only LSDA can be parallelized over two pools, and there
    ! is no need to communicate anything: each pools deals with its own spin
    !
    IF ( .NOT.gamma_only ) THEN
       DO ikq = 1, nkqs
         CALL mp_bcast(exxbuff(:,:,ikq), working_pool(ikq), intra_orthopool_comm)
       ENDDO
    END IF
    !
    ! For US/PAW only: compute <beta_I|psi_j,k+q> for the entire 
    ! de-symmetrized k+q grid by rotating the ones from the irreducible wedge
    !
    IF(okvan) CALL rotate_becxx(nkqs, index_xk, index_sym, xkq_collect)
    !
    ! Initialize 4-wavefunctions one-center Fock integrals
    !    \int \psi_a(r)\phi_a(r)\phi_b(r')\psi_b(r')/|r-r'|
    !
    IF(okpaw) CALL PAW_init_fock_kernel()
    !
    CALL change_data_structure(.FALSE.)
    !
#if defined(__CUDA)
    exxbuff_d=exxbuff
#endif
    CALL stop_clock ('exxinit')
    !
    !-----------------------------------------------------------------------
  END SUBROUTINE exxinit
  !-----------------------------------------------------------------------
  !
  !-----------------------------------------------------------------------
  SUBROUTINE vexx(lda, n, m, psi, hpsi, becpsi)
  !-----------------------------------------------------------------------
    !
    ! ... Wrapper routine computing V_x\psi, V_x = exchange potential
    ! ... Calls generic version vexx_k or Gamma-specific one vexx_gamma
    !
    ! ... input:
    ! ...    lda   leading dimension of arrays psi and hpsi
    ! ...    n     true dimension of psi and hpsi
    ! ...    m     number of states psi
    ! ...    psi   m wavefunctions
    ! ..     becpsi <beta|psi>, optional but needed for US and PAW case
    !
    ! ... output:
    ! ...    hpsi  V_x*psi
    !
    USE becmod,         ONLY : bec_type
    USE uspp,           ONLY : okvan
    USE paw_variables,  ONLY : okpaw
    USE us_exx,         ONLY : becxx
    USE mp_exx,         ONLY : negrp, inter_egrp_comm, init_index_over_band
    USE wvfct,          ONLY : nbnd
    USE exx_band,       ONLY : transform_psi_to_exx, transform_hpsi_to_local,&
                               psi_exx, hpsi_exx, igk_exx
    !
    IMPLICIT NONE
    !
    INTEGER                  :: lda, n, m
    COMPLEX(DP)              :: psi(lda*npol,m)
    COMPLEX(DP)              :: hpsi(lda*npol,m)
    TYPE(bec_type), OPTIONAL :: becpsi
    INTEGER :: i
    !
    IF ( (okvan.or.okpaw) .and. .not. present(becpsi)) &
       CALL errore('vexx','becpsi needed for US/PAW case',1)
    CALL start_clock ('vexx')
    !
    IF(negrp.gt.1)THEN
       CALL init_index_over_band(inter_egrp_comm,nbnd,m)
       !
       ! transform psi to the EXX data structure
       !
       CALL transform_psi_to_exx(lda,n,m,psi)
       !
    END IF
    !
    ! calculate the EXX contribution to hpsi
    !
    IF(gamma_only) THEN
       IF(negrp.eq.1)THEN
          CALL vexx_gamma(lda, n, m, psi, hpsi, becpsi)
       ELSE
          CALL vexx_gamma(lda, n, m, psi_exx, hpsi_exx, becpsi)
       END IF
    ELSE
       IF(negrp.eq.1)THEN
          CALL vexx_k(lda, n, m, psi, hpsi, becpsi)
       ELSE
          CALL vexx_k(lda, n, m, psi_exx, hpsi_exx, becpsi)
       END IF
    ENDIF
    !
    IF(negrp.gt.1)THEN
       !
       ! transform hpsi to the local data structure
       !
       CALL transform_hpsi_to_local(lda,n,m,hpsi)
       !
    END IF
    !
    CALL stop_clock ('vexx')
    !
    !-----------------------------------------------------------------------
  END SUBROUTINE vexx
  !-----------------------------------------------------------------------
  !
  !-----------------------------------------------------------------------
  SUBROUTINE vexx_gamma(lda, n, m, psi, hpsi, becpsi)
  !-----------------------------------------------------------------------
    !
    ! ... Gamma-specific version of vexx
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
    !
    IMPLICIT NONE
    !
    INTEGER                  :: lda, n, m
    COMPLEX(DP)              :: psi(lda*npol,m)
    COMPLEX(DP)              :: hpsi(lda*npol,m)
    TYPE(bec_type), OPTIONAL :: becpsi ! or call a calbec(...psi) instead
    !
    ! local variables
    COMPLEX(DP),ALLOCATABLE :: result(:,:)
    REAL(DP),ALLOCATABLE :: temppsic_dble (:)
    REAL(DP),ALLOCATABLE :: temppsic_aimag(:)
    !
    COMPLEX(DP),ALLOCATABLE :: vc(:,:), deexx(:,:)
    REAL(DP),   ALLOCATABLE :: fac(:)
    INTEGER          :: ibnd, ik, im , ikq, iq, ipol
    INTEGER          :: ir, ig
    INTEGER          :: current_ik
    INTEGER          :: ibnd_loop_start
    INTEGER          :: nrxxs
    REAL(DP) :: x1, x2, xkp(3)
    REAL(DP) :: xkq(3)
    INTEGER, EXTERNAL :: global_kpoint_index
    INTEGER :: ialloc
    COMPLEX(DP), ALLOCATABLE :: big_result(:,:)
    INTEGER :: iproc, nproc_egrp, ii, ipair
    INTEGER :: jbnd, jstart, jend
    ! scratch space for fft of psi and rho
    COMPLEX(DP), ALLOCATABLE :: psi_rhoc_work(:)
    INTEGER :: jblock_start, jblock_end
    INTEGER :: iegrp, wegrp
    INTEGER :: exxbuff_index
    INTEGER :: ending_im
    !
    ialloc = nibands(my_egrp_id+1)
    !
    ALLOCATE( fac(dfftt%ngm) )
    nrxxs= dfftt%nnr
    !
    !ALLOCATE( result(nrxxs), temppsic_dble(nrxxs), temppsic_aimag(nrxxs) )
    ALLOCATE( result(nrxxs,ialloc), temppsic_dble(nrxxs) )
    ALLOCATE( temppsic_aimag(nrxxs) )
    ALLOCATE( psi_rhoc_work(nrxxs) )
    !
    ALLOCATE( vc(nrxxs,ialloc))
    IF(okvan) ALLOCATE(deexx(nkb,ialloc))
    !
    current_ik = global_kpoint_index ( nkstot, current_k )
    xkp = xk(:,current_k)
    !
    allocate(big_result(n,m))
    big_result = 0.0_DP
    result = 0.0_DP
    !
    DO ii=1, nibands(my_egrp_id+1)
       IF(okvan) deexx(:,ii) = 0.0_DP
    END DO
    !
    ! Here the loops start
    !
    INTERNAL_LOOP_ON_Q : &
    DO iq=1,nqs
       !
       ikq  = index_xkq(current_ik,iq)
       ik   = index_xk(ikq)
       xkq  = xkq_collect(:,ikq)
       !
       ! calculate the 1/|r-r'| (actually, k+q+g) factor and place it in fac
       CALL g2_convolution_all(dfftt%ngm, gt, xkp, xkq, iq, current_k)
       IF ( okvan .and..not.tqr ) CALL qvan_init (dfftt%ngm, xkq, xkp)
       !
       DO iegrp=1, negrp
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
             IF (ibnd.eq.0.or.ibnd.gt.m) CYCLE
             !
             IF ( mod(ii,2)==1 ) THEN
                !
                psi_rhoc_work = (0._DP,0._DP)
                !
                IF ( (ii+1)<=min(m,nibands(my_egrp_id+1)) ) THEN
                   ! deal with double bands
!$omp parallel do  default(shared), private(ig)
                   DO ig = 1, npwt
                      psi_rhoc_work( dfftt%nl(ig) )  =       psi(ig, ii) + (0._DP,1._DP) * psi(ig, ii+1)
                      psi_rhoc_work( dfftt%nlm(ig) ) = conjg(psi(ig, ii) - (0._DP,1._DP) * psi(ig, ii+1))
                   ENDDO
!$omp end parallel do
                ENDIF
                !
                IF ( ii==min(m,nibands(my_egrp_id+1)) ) THEN
                   ! deal with a single last band
!$omp parallel do  default(shared), private(ig)
                   DO ig = 1, npwt
                      psi_rhoc_work( dfftt%nl(ig) )  =       psi(ig,ii)
                      psi_rhoc_work( dfftt%nlm(ig) ) = conjg(psi(ig,ii))
                   ENDDO
!$omp end parallel do
                ENDIF
                !
                CALL invfft ('Wave', psi_rhoc_work, dfftt)
!$omp parallel do default(shared), private(ir)
                DO ir = 1, nrxxs
                   temppsic_dble(ir)  = dble ( psi_rhoc_work(ir) )
                   temppsic_aimag(ir) = aimag( psi_rhoc_work(ir) )
                ENDDO
!$omp end parallel do
                !
             ENDIF
             !
             !
             !determine which j-bands to calculate
             jstart = 0
             jend = 0
             DO ipair=1, max_pairs
                IF(egrp_pairs(1,ipair,my_egrp_id+1).eq.ibnd)THEN
                   IF(jstart.eq.0)THEN
                      jstart = egrp_pairs(2,ipair,my_egrp_id+1)
                      jend = jstart
                   ELSE
                      jend = egrp_pairs(2,ipair,my_egrp_id+1)
                   END IF
                END IF
             END DO
             !
             jstart = max(jstart,jblock_start)
             jend = min(jend,jblock_end)
             !
             IF(mod(jstart,2)==0) THEN
                ibnd_loop_start=jstart-1
             ELSE
                ibnd_loop_start=jstart
             ENDIF
             !
             IBND_LOOP_GAM : &
             DO jbnd=ibnd_loop_start,jend, 2 !for each band of psi
                !
                exxbuff_index = (jbnd+1)/2-(all_start(wegrp)+1)/2+(iexx_start+1)/2
                !
                IF( jbnd < jstart ) THEN
                   x1 = 0.0_DP
                ELSE
                   x1 = x_occupation(jbnd,  ik)
                ENDIF
                IF( jbnd == jend) THEN
                   x2 = 0.0_DP
                ELSE
                   x2 = x_occupation(jbnd+1,  ik)
                ENDIF
                IF ( abs(x1) < eps_occ .and. abs(x2) < eps_occ ) CYCLE
                !
                ! calculate rho in real space. Gamma tricks are used.
                ! temppsic is real; tempphic contains one band in the real part,
                ! another one in the imaginary part; the same applies to rhoc
                !
                IF( mod(ii,2) == 0 ) THEN
!$omp parallel do default(shared), private(ir)
                   DO ir = 1, nrxxs
                      psi_rhoc_work(ir) = exxbuff(ir,exxbuff_index,ikq) * temppsic_aimag(ir) / omega
                   ENDDO
!$omp end parallel do
                ELSE
!$omp parallel do default(shared), private(ir)
                   DO ir = 1, nrxxs
                      psi_rhoc_work(ir) = exxbuff(ir,exxbuff_index,ikq) * temppsic_dble(ir) / omega
                   ENDDO
!$omp end parallel do
                ENDIF
                !
                ! bring rho to G-space
                !
                !   >>>> add augmentation in REAL SPACE here
                IF(okvan .and. tqr) THEN
                   IF(jbnd>=jstart) &
                        CALL addusxx_r(psi_rhoc_work, &
                       _CX(becxx(ikq)%r(:,jbnd)), _CX(becpsi%r(:,ibnd)))
                   IF(jbnd<jend) &
                        CALL addusxx_r(psi_rhoc_work, &
                       _CY(becxx(ikq)%r(:,jbnd+1)),_CX(becpsi%r(:,ibnd)))
                ENDIF
                !
                CALL fwfft ('Rho', psi_rhoc_work, dfftt)
                !   >>>> add augmentation in G SPACE here
                IF(okvan .and. .not. tqr) THEN
                   ! contribution from one band added to real (in real space) part of rhoc
                   IF(jbnd>=jstart) &
                        CALL addusxx_g(dfftt, psi_rhoc_work, xkq,  xkp, 'r', &
                        becphi_r=becxx(ikq)%r(:,jbnd), becpsi_r=becpsi%r(:,ibnd) )
                   ! contribution from following band added to imaginary (in real space) part of rhoc
                   IF(jbnd<jend) &
                        CALL addusxx_g(dfftt, psi_rhoc_work, xkq,  xkp, 'i', &
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
                IF(okvan .and. .not. tqr) THEN
                   IF(jbnd>=jstart) &
                        CALL newdxx_g(dfftt, vc(:,ii), xkq, xkp, 'r', deexx(:,ii), &
                           becphi_r=x1*becxx(ikq)%r(:,jbnd))
                   IF(jbnd<jend) &
                        CALL newdxx_g(dfftt, vc(:,ii), xkq, xkp, 'i', deexx(:,ii), &
                            becphi_r=x2*becxx(ikq)%r(:,jbnd+1))
                ENDIF
                !
                !brings back v in real space
                CALL invfft ('Rho', vc(:,ii), dfftt)
                !
                !   >>>>  compute <psi|H_fock REAL SPACE here
                IF(okvan .and. tqr) THEN
                   IF(jbnd>=jstart) &
                        CALL newdxx_r(dfftt,vc(:,ii), _CX(x1*becxx(ikq)%r(:,jbnd)), deexx(:,ii))
                   IF(jbnd<jend) &
                        CALL newdxx_r(dfftt,vc(:,ii), _CY(x2*becxx(ikq)%r(:,jbnd+1)), deexx(:,ii))
                ENDIF
                !
                IF(okpaw) THEN
                   IF(jbnd>=jstart) &
                        CALL PAW_newdxx(x1/nqs, _CX(becxx(ikq)%r(:,jbnd)),&
                                                _CX(becpsi%r(:,ibnd)), deexx(:,ii))
                   IF(jbnd<jend) &
                        CALL PAW_newdxx(x2/nqs, _CX(becxx(ikq)%r(:,jbnd+1)),&
                                                _CX(becpsi%r(:,ibnd)), deexx(:,ii))
                ENDIF
                !
                ! accumulates over bands and k points
                !
!$omp parallel do default(shared), private(ir)
                DO ir = 1, nrxxs
                   result(ir,ii) = result(ir,ii) &
                                 + x1* dble(vc(ir,ii))* dble(exxbuff(ir,exxbuff_index,ikq)) &
                                 + x2*aimag(vc(ir,ii))*aimag(exxbuff(ir,exxbuff_index,ikq))
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
          IF (negrp>1) call mp_circular_shift_left( exxbuff(:,:,ikq), me_egrp, inter_egrp_comm )
          !
       ENDDO ! iegrp
       IF ( okvan .and..not.tqr ) CALL qvan_clean ()
    ENDDO &
    INTERNAL_LOOP_ON_Q
    !
    DO ii=1, nibands(my_egrp_id+1)
       !
       ibnd = ibands(ii,my_egrp_id+1)
       !
       IF (ibnd.eq.0.or.ibnd.gt.m) CYCLE
       !
       IF(okvan) THEN
          CALL mp_sum(deexx(:,ii),intra_egrp_comm)
       ENDIF
       !
       !
       ! brings back result in G-space
       !
       CALL fwfft( 'Wave' , result(:,ii), dfftt )
       !communicate result
       DO ig = 1, n
          big_result(ig,ibnd) = big_result(ig,ibnd) - exxalfa*result(dfftt%nl(igk_exx(ig,current_k)),ii)
       END DO
       !
       ! add non-local \sum_I |beta_I> \alpha_Ii (the sum on i is outside)
       IF(okvan) CALL add_nlxx_pot (lda, big_result(:,ibnd), xkp, n, &
            igk_exx(1,current_k), deexx(:,ii), eps_occ, exxalfa)
    END DO
    !
    CALL result_sum(n*npol, m, big_result)
    IF (iexx_istart(my_egrp_id+1).gt.0) THEN
       IF (negrp == 1) then
          ending_im = m
       ELSE
          ending_im = iexx_iend(my_egrp_id+1) - iexx_istart(my_egrp_id+1) + 1
       END IF
       DO im=1, ending_im
!$omp parallel do default(shared), private(ig) firstprivate(im,n)
           DO ig = 1, n
              hpsi(ig,im)=hpsi(ig,im) + big_result(ig,im+iexx_istart(my_egrp_id+1)-1)
           ENDDO
!$omp end parallel do
       END DO
    END IF
    !
    DEALLOCATE(big_result)
    DEALLOCATE( result, temppsic_dble, temppsic_aimag)
    DEALLOCATE( vc, fac )
    IF(okvan) DEALLOCATE( deexx )
    !
    !-----------------------------------------------------------------------
  END SUBROUTINE vexx_gamma
  !-----------------------------------------------------------------------
  !
  !
  !-----------------------------------------------------------------------
  SUBROUTINE vexx_k(lda, n, m, psi, hpsi, becpsi)
  !-----------------------------------------------------------------------
    !
    ! ... generic, k-point version of vexx
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
#if defined(__CUDA)
    USE exx_band,       ONLY : igk_exx_d
#endif
    USE io_global,      ONLY : stdout
    !
    !
    IMPLICIT NONE

#if defined(__CUDA)
#define DEV_ATTRIBUTES , DEVICE
#else
#define DEV_ATTRIBUTES 
#endif

    !
    INTEGER                  :: lda, n, m
    COMPLEX(DP)              :: psi(:,:)
#if defined(__CUDA)
    COMPLEX(DP),ALLOCATABLE,DEVICE :: psi_d(:,:)
#endif
    COMPLEX(DP)              :: hpsi(lda*npol,max_ibands)
    TYPE(bec_type), OPTIONAL :: becpsi ! or call a calbec(...psi) instead
    !

    ! local variables
    COMPLEX(DP),ALLOCATABLE DEV_ATTRIBUTES :: temppsic(:,:)
    COMPLEX(DP),ALLOCATABLE DEV_ATTRIBUTES :: temppsic_nc(:,:,:)

    COMPLEX(DP),ALLOCATABLE :: result(:,:), result_nc(:,:,:)
#if defined(__CUDA)
    COMPLEX(DP),ALLOCATABLE,DEVICE :: result_d(:,:), result_nc_d(:,:,:)
#endif

    INTEGER          :: request_send, request_recv
    !
    COMPLEX(DP),ALLOCATABLE :: deexx(:,:)

#if defined(__CUDA)
    COMPLEX(DP),ALLOCATABLE,TARGET,DEVICE :: rhoc_d(:,:), vc_d(:,:)
#endif
    COMPLEX(DP),ALLOCATABLE,TARGET :: rhoc(:,:), vc(:,:)

!#if defined(__USE_MANY_FFT)
!    COMPLEX(DP),POINTER :: prhoc(:), pvc(:)
!#endif


    REAL(DP),   ALLOCATABLE :: fac(:), facb(:)
#if defined(__CUDA)
    REAL(DP),ALLOCATABLE,DEVICE :: facb_d(:)
#endif
    INTEGER  :: ibnd, ik, im , ikq, iq, ipol
    INTEGER  :: ir, ig, ir_start, ir_end
    INTEGER  :: irt, nrt, nblock
    INTEGER  :: current_ik
    INTEGER  :: ibnd_loop_start
    INTEGER  :: nrxxs
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
    INTEGER :: all_start_tmp
    !hack around PGI bug
#if defined(__CUDA)
    INTEGER, POINTER, DEVICE :: dfftt__nl(:)
    dfftt__nl=>dfftt%nl_d
#endif
    CALL start_clock( 'vexx_k_setup' )


    ialloc = nibands(my_egrp_id+1)
    !
    ALLOCATE( fac(dfftt%ngm) )
    nrxxs= dfftt%nnr
    ALLOCATE( facb(nrxxs) )
#if defined(__CUDA)
    ALLOCATE( psi_d, source=psi )
    ALLOCATE( facb_d(nrxxs) )
#endif
    !

    IF (noncolin) THEN
#if defined(__CUDA)
       ALLOCATE( result_nc_d(nrxxs,npol,ialloc) )
#endif
       ALLOCATE( result_nc(nrxxs,npol,ialloc) )

       !temppsic knows where it is
       ALLOCATE( temppsic_nc(nrxxs,npol,ialloc) )
    ELSE
#if defined(__CUDA)
       ALLOCATE( result_d(nrxxs,ialloc) )
#endif
       ALLOCATE( result(nrxxs,ialloc) )
       
       !temppsic knows
       ALLOCATE( temppsic(nrxxs,ialloc) )
    ENDIF
    !
    IF(okvan) ALLOCATE(deexx(nkb,ialloc))
    !
    current_ik = global_kpoint_index ( nkstot, current_k )
    xkp = xk(:,current_k)
    !
    allocate(big_result(n*npol,m))
    big_result = 0.0_DP
    !
    !allocate arrays for rhoc and vc
#if defined(__CUDA)
    ALLOCATE(rhoc_d(nrxxs,jblock), vc_d(nrxxs,jblock))
#endif
    ALLOCATE(rhoc(nrxxs,jblock), vc(nrxxs,jblock))


!#if defined(__USE_MANY_FFT)
!    prhoc(1:nrxxs*jblock) => rhoc(:,:)
!    pvc(1:nrxxs*jblock) => vc(:,:)
!#endif
    !

    DO ii=1, nibands(my_egrp_id+1)
       !
       ibnd = ibands(ii,my_egrp_id+1)
       !
       IF (ibnd.eq.0.or.ibnd.gt.m) CYCLE
       !
       IF(okvan) deexx(:,ii) = 0._DP
       !
       IF (noncolin) THEN
          temppsic_nc(:,:,ii) = 0._DP
       ELSE
#if defined(__CUDA)
          temppsic(:,ii) = 0._DP
#else
!$omp parallel do  default(shared), private(ir) firstprivate(nrxxs)
          DO ir = 1, nrxxs
             temppsic(ir,ii) = 0._DP
          ENDDO
#endif
       END IF

#if defined(__CUDA)
       !
       IF (noncolin) THEN
          !$cuf kernel do (1)
          DO ig = 1, n
             temppsic_nc(dfftt__nl(igk_exx_d(ig,current_k)),1,ii) = psi_d(ig,ii)
             temppsic_nc(dfftt__nl(igk_exx_d(ig,current_k)),2,ii) = psi_d(npwx+ig,ii)
          ENDDO
          CALL invfft ('Wave', temppsic_nc(:,1,ii), dfftt)
          CALL invfft ('Wave', temppsic_nc(:,2,ii), dfftt)
       ELSE
          !$cuf kernel do (1)
          DO ig = 1, n
             temppsic( dfftt__nl(igk_exx_d(ig,current_k)), ii ) = psi_d(ig,ii)
          ENDDO
          CALL invfft ('Wave', temppsic(:,ii), dfftt)
       END IF
#else
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
          CALL invfft ('Wave', temppsic_nc(:,1,ii), dfftt)
          CALL invfft ('Wave', temppsic_nc(:,2,ii), dfftt)
          !
       ELSE
          !
!$omp parallel do  default(shared), private(ig)
          DO ig = 1, n
             temppsic( dfftt%nl(igk_exx(ig,current_k)), ii ) = psi(ig,ii)
          ENDDO
!$omp end parallel do
          !
          CALL invfft ('Wave', temppsic(:,ii), dfftt)
          !
       END IF
#endif

       !
#if defined (__CUDA)
       IF (noncolin) THEN
          result_nc_d(:,:,ii) = 0.0_DP
       ELSE
          result_d(:,ii) = 0.0_DP
       ENDIF
#else
       IF (noncolin) THEN
!$omp parallel do default(shared) firstprivate(nrxxs) private(ir)
          DO ir=1,nrxxs
             result_nc(ir,1,ii) = 0.0_DP
             result_nc(ir,2,ii) = 0.0_DP
          ENDDO
       ELSE
!$omp parallel do default(shared) firstprivate(nrxxs) private(ir)
          DO ir=1,nrxxs
             result(ir,ii) = 0.0_DP
          ENDDO
       END IF
       !
#endif
    END DO

#if defined (__CUDA)
    ! no longer need psi_d
    DEALLOCATE(psi_d)
#endif

    !
    !precompute these guys
    omega_inv = 1.0 / omega
    nqs_inv = 1.0 / nqs
    !
    CALL stop_clock( 'vexx_k_setup' )
    CALL start_clock( 'vexx_k_main' )
    !------------------------------------------------------------------------!
    ! Beginning of main loop
    !------------------------------------------------------------------------!
    vexxmain: DO iq=1, nqs
       !
       ikq  = index_xkq(current_ik,iq)
       ik   = index_xk(ikq)
       xkq  = xkq_collect(:,ikq)
       !
       ! calculate the 1/|r-r'| (actually, k+q+g) factor and place it in fac
       CALL g2_convolution_all(dfftt%ngm, gt, xkp, xkq, iq, current_k)
       !
! JRD - below not threaded
       facb = 0D0
       DO ig = 1, dfftt%ngm
          facb(dfftt%nl(ig)) = coulomb_fac(ig,iq,current_k)
       ENDDO
#if defined(__CUDA)
       facb_d = facb
#endif
       !
       IF ( okvan .and..not.tqr ) CALL qvan_init (dfftt%ngm, xkq, xkp)
       !
       DO iegrp=1, negrp
          !
          ! compute the id of group whose data is currently worked on
          wegrp = MOD(iegrp+my_egrp_id-1, negrp)+1
          njt = (all_end(wegrp)-all_start(wegrp)+jblock)/jblock
          !
          DO ijt=1, njt
             !
             jblock_start = (ijt - 1) * jblock + all_start(wegrp)
             jblock_end = min(jblock_start+jblock-1,all_end(wegrp))
             !
             DO ii=1, nibands(my_egrp_id+1)
                !
                ibnd = ibands(ii,my_egrp_id+1)
                !
                IF (ibnd.eq.0.or.ibnd.gt.m) CYCLE
                !
                !determine which j-bands to calculate
                jstart = 0
                jend = 0
                DO ipair=1, max_pairs
                   IF(egrp_pairs(1,ipair,my_egrp_id+1).eq.ibnd)THEN
                      IF(jstart.eq.0)THEN
                         jstart = egrp_pairs(2,ipair,my_egrp_id+1)
                         jend = jstart
                      ELSE
                         jend = egrp_pairs(2,ipair,my_egrp_id+1)
                      END IF
                   END IF
                END DO
                !
                jstart = max(jstart,jblock_start)
                jend = min(jend,jblock_end)
                !
                !how many iters
                jcount=jend-jstart+1
                if(jcount<=0) cycle
                !
                !----------------------------------------------------------------------!
                !INNER LOOP START
                !----------------------------------------------------------------------!
                !
                nblock=2048
                nrt = (nrxxs+nblock-1)/nblock
                !
#if defined(__CUDA)
associate(rhoc=>rhoc_d, exxbuff=>exxbuff_d)
                all_start_tmp=all_start(wegrp)
                !$cuf kernel do (2)
                DO jbnd=jstart, jend
                   DO ir = 1, nrxxs

                     IF (noncolin) THEN
                       rhoc(ir,jbnd-jstart+1) = &
                       (conjg(exxbuff(ir,jbnd-all_start_tmp+iexx_start,ikq))*temppsic_nc(ir,1,ii) +&
                       conjg(exxbuff(nrxxs+ir,jbnd-all_start_tmp+iexx_start,ikq))*temppsic_nc(ir,2,ii)) * omega_inv
                     ELSE

                       rhoc(ir,jbnd-jstart+1) = &
                       conjg(exxbuff(ir,jbnd-all_start_tmp+iexx_start,ikq))*temppsic(ir,ii)* omega_inv
                     ENDIF

                   ENDDO
                ENDDO
end associate
#else

!$omp parallel do collapse(2) private(ir_start,ir_end)
                DO irt = 1, nrt
                   DO jbnd=jstart, jend
                      ir_start = (irt - 1) * nblock + 1
                      ir_end = min(ir_start+nblock-1,nrxxs)
                      IF (noncolin) THEN
                         DO ir = ir_start, ir_end
                            rhoc(ir,jbnd-jstart+1) = &
                               ( conjg(exxbuff(ir,jbnd-all_start(wegrp)+iexx_start,ikq))*temppsic_nc(ir,1,ii) +&
                               conjg(exxbuff(nrxxs+ir,jbnd-all_start(wegrp)+iexx_start,ikq))*temppsic_nc(ir,2,ii) )/omega
                         ENDDO
                      ELSE
!DIR$ vector nontemporal (rhoc)
                         DO ir = ir_start, ir_end
                            rhoc(ir,jbnd-jstart+1) = &
                               conjg(exxbuff(ir,jbnd-all_start(wegrp)+iexx_start,ikq))*temppsic(ir,ii) * omega_inv
                         ENDDO
                      ENDIF
                   ENDDO
                ENDDO
!$omp end parallel do
#endif
                !
                !   >>>> add augmentation in REAL space HERE
                IF(okvan .and. tqr) THEN ! augment the "charge" in real space
                   DO jbnd=jstart, jend
                      CALL addusxx_r(rhoc(:,jbnd-jstart+1), becxx(ikq)%k(:,jbnd), becpsi%k(:,ibnd))
                   ENDDO
                ENDIF
                !
                !   >>>> brings it to G-space
!#if defined(__USE_MANY_FFT)
!                CALL fwfft ('Rho', prhoc, dfftt, howmany=jcount)
!#else
#if defined (__CUDA)
                !rhoc_d = rhoc
                DO jbnd=jstart, jend
                   CALL fwfft('Rho', rhoc_d(:,jbnd-jstart+1), dfftt)
                ENDDO
                !rhoc = rhoc_d
#else
                DO jbnd=jstart, jend
                   CALL fwfft('Rho', rhoc(:,jbnd-jstart+1), dfftt)
                ENDDO
#endif
!#endif
                !
                !   >>>> add augmentation in G space HERE
                IF(okvan .and. .not. tqr) THEN
                   print*, "In bad branch!"
                   DO jbnd=jstart, jend
                      CALL addusxx_g(dfftt, rhoc(:,jbnd-jstart+1), xkq, xkp, &
                      'c', becphi_c=becxx(ikq)%k(:,jbnd),becpsi_c=becpsi%k(:,ibnd))
                   ENDDO
                ENDIF
                !   >>>> charge done
                !
#if defined(__CUDA)
associate(vc=>vc_d, facb=>facb_d, rhoc=>rhoc_d, x_occupation=>x_occupation_d)
                !$cuf kernel do (2)
                DO jbnd=jstart, jend
                   DO ir = 1, nrxxs
                         vc(ir,jbnd-jstart+1) = facb(ir) * rhoc(ir,jbnd-jstart+1)*&
                                                x_occupation(jbnd,ik) * nqs_inv
                   ENDDO
                ENDDO
end associate
#else
!call start_collection()
!$omp parallel do collapse(2) private(ir_start,ir_end)
                DO irt = 1, nrt
                   DO jbnd=jstart, jend
                      ir_start = (irt - 1) * nblock + 1
                      ir_end = min(ir_start+nblock-1,nrxxs)
!DIR$ vector nontemporal (vc)
                      DO ir = ir_start, ir_end
                         vc(ir,jbnd-jstart+1) = facb(ir) * rhoc(ir,jbnd-jstart+1)*&
                                                x_occupation(jbnd,ik) * nqs_inv
                      ENDDO
                   ENDDO
                ENDDO
!$omp end parallel do
!call stop_collection()
#endif
                !
                ! Add ultrasoft contribution (RECIPROCAL SPACE)
                ! compute alpha_I,j,k+q = \sum_J \int <beta_J|phi_j,k+q> V_i,j,k,q Q_I,J(r) d3r
                IF(okvan .and. .not. tqr) THEN
#if defined (__CUDA)
                   vc = vc_d
#endif
                   DO jbnd=jstart, jend
                      CALL newdxx_g(dfftt, vc(:,jbnd-jstart+1), xkq, xkp, 'c',&
                                    deexx(:,ii), becphi_c=becxx(ikq)%k(:,jbnd))
                   ENDDO
#if defined (__CUDA)
                   vc_d = vc
#endif
                ENDIF
                !
                !brings back v in real space
!#if defined(__USE_MANY_FFT)
!                !fft many
!                CALL invfft ('Rho', pvc, dfftt, howmany=jcount)
!#else
#if defined (__CUDA)
                DO jbnd=jstart, jend
                   CALL invfft('Rho', vc_d(:,jbnd-jstart+1), dfftt)
                ENDDO
#else
                DO jbnd=jstart, jend
                   CALL invfft('Rho', vc(:,jbnd-jstart+1), dfftt)
                ENDDO
#endif
!#endif
                !
                ! Add ultrasoft contribution (REAL SPACE)
                IF(okvan .and. tqr) THEN
#if defined(__CUDA)
                   vc = vc_d
#endif
                   DO jbnd=jstart, jend
                      CALL newdxx_r(dfftt, vc(:,jbnd-jstart+1), becxx(ikq)%k(:,jbnd),deexx(:,ii))
                   ENDDO
#if defined(__CUDA)
                   vc_d = vc
#endif
                ENDIF
                !
                ! Add PAW one-center contribution
                IF(okpaw) THEN
#if defined(__CUDA)
                   vc = vc_d
#endif
                   DO jbnd=jstart, jend
                      CALL PAW_newdxx(x_occupation(jbnd,ik)/nqs, becxx(ikq)%k(:,jbnd), becpsi%k(:,ibnd), deexx(:,ii))
                   ENDDO
#if defined(__CUDA)
                   vc_d = vc
#endif
                ENDIF
                !
                !accumulates over bands and k points
                !

#if defined(__CUDA)    
associate(exxbuff=>exxbuff_d, result_nc=>result_nc_d, result=>result_d, vc=>vc_d)
                all_start_tmp=all_start(wegrp)
                DO jbnd=jstart, jend
                   !$cuf kernel do (1)
                   DO ir = 1, nrxxs
                      IF (noncolin) THEN
                         result_nc(ir,1,ii) = result_nc(ir,1,ii) &
                              + vc(ir,jbnd-jstart+1) * exxbuff(ir,jbnd-all_start_tmp+iexx_start,ikq)
                         result_nc(ir,2,ii) = result_nc(ir,2,ii) &
                              + vc(ir,jbnd-jstart+1) * exxbuff(ir+nrxxs,jbnd-all_start_tmp+iexx_start,ikq)
                      ELSE
                         result(ir,ii) = result(ir,ii) &
                              + vc(ir,jbnd-jstart+1)*exxbuff(ir,jbnd-all_start_tmp+iexx_start,ikq)
                      ENDIF
                   ENDDO
                ENDDO
end associate
#else
!call start_collection()
!$omp parallel do private(ir_start,ir_end)
                DO irt = 1, nrt
                   DO jbnd=jstart, jend
                      ir_start = (irt - 1) * nblock + 1
                      ir_end = min(ir_start+nblock-1,nrxxs)
                      IF (noncolin) THEN
                         DO ir = ir_start, ir_end
                            result_nc(ir,1,ii) = result_nc(ir,1,ii) &
                                               + vc(ir,jbnd-jstart+1) * exxbuff(ir,jbnd-all_start(wegrp)+iexx_start,ikq)
                            result_nc(ir,2,ii) = result_nc(ir,2,ii) &
                                               + vc(ir,jbnd-jstart+1) * exxbuff(ir+nrxxs,jbnd-all_start(wegrp)+iexx_start,ikq)
                         ENDDO
                      ELSE
!!dir$ vector nontemporal (result)
                         DO ir = ir_start, ir_end
                            result(ir,ii) = result(ir,ii) &
                                          + vc(ir,jbnd-jstart+1)*exxbuff(ir,jbnd-all_start(wegrp)+iexx_start,ikq)
                         ENDDO
                      ENDIF
                   ENDDO
                ENDDO
!$omp end parallel do
!call stop_collection()
#endif

                !
                !----------------------------------------------------------------------!
                !INNER LOOP END
                !----------------------------------------------------------------------!
                !
             END DO !I-LOOP
          END DO !IJT
          !
          ! get the next nbnd/negrp data
          IF (negrp>1) call mp_circular_shift_left( exxbuff(:,:,ikq), me_egrp, inter_egrp_comm )
          !
       END DO !iegrp
       !
       IF ( okvan .and..not.tqr ) CALL qvan_clean ()
    END DO vexxmain

#if defined(__CUDA)
    IF (noncolin) THEN
       result_nc = result_nc_d
    ELSE
       result = result_d
    ENDIF
#endif

    CALL stop_clock( 'vexx_k_main' )
    CALL start_clock( 'vexx_k_fin' )
    !
    !
    !
    DO ii=1, nibands(my_egrp_id+1)
       !
       ibnd = ibands(ii,my_egrp_id+1)
       !
       IF (ibnd.eq.0.or.ibnd.gt.m) CYCLE
       !
       IF(okvan) THEN
          CALL mp_sum(deexx(:,ii),intra_egrp_comm)
       ENDIF
       !
       IF (noncolin) THEN
          !brings back result in G-space
          CALL fwfft ('Wave', result_nc(:,1,ii), dfftt)
          CALL fwfft ('Wave', result_nc(:,2,ii), dfftt)
          DO ig = 1, n
             big_result(ig,ibnd) = big_result(ig,ibnd) - exxalfa*result_nc(dfftt%nl(igk_exx(ig,current_k)),1,ii)
             big_result(n+ig,ibnd) = big_result(n+ig,ibnd) - exxalfa*result_nc(dfftt%nl(igk_exx(ig,current_k)),2,ii)
          ENDDO
       ELSE
          !
          CALL fwfft ('Wave', result(:,ii), dfftt)
          DO ig = 1, n
             big_result(ig,ibnd) = big_result(ig,ibnd) - exxalfa*result(dfftt%nl(igk_exx(ig,current_k)),ii)
          ENDDO
       ENDIF
       !
       ! add non-local \sum_I |beta_I> \alpha_Ii (the sum on i is outside)
       IF(okvan) CALL add_nlxx_pot (lda, big_result(:,ibnd), xkp, n, igk_exx(:,current_k),&
            deexx(:,ii), eps_occ, exxalfa)
       !
    END DO
    !
    !deallocate temporary arrays
    DEALLOCATE(rhoc, vc)
#if defined(__CUDA)
    DEALLOCATE(rhoc_d, vc_d)
#endif
    !
    !sum result
    CALL result_sum(n*npol, m, big_result)
    IF (iexx_istart(my_egrp_id+1).gt.0) THEN
       IF (negrp == 1) then
          ending_im = m
       ELSE
          ending_im = iexx_iend(my_egrp_id+1) - iexx_istart(my_egrp_id+1) + 1
       END IF
       IF(noncolin) THEN
          DO im=1, ending_im
!$omp parallel do default(shared), private(ig) firstprivate(im,n)
             DO ig = 1, n
                hpsi(ig,im)=hpsi(ig,im) + big_result(ig,im+iexx_istart(my_egrp_id+1)-1)
             ENDDO
!$omp end parallel do
!$omp parallel do default(shared), private(ig) firstprivate(im,n)
             DO ig = 1, n
                hpsi(lda+ig,im)=hpsi(lda+ig,im) + big_result(n+ig,im+iexx_istart(my_egrp_id+1)-1)
             ENDDO
!$omp end parallel do
          END DO
       ELSE
          DO im=1, ending_im
!$omp parallel do default(shared), private(ig) firstprivate(im,n)
             DO ig = 1, n
                hpsi(ig,im)=hpsi(ig,im) + big_result(ig,im+iexx_istart(my_egrp_id+1)-1)
             ENDDO
!$omp end parallel do
          ENDDO
       END IF
    END IF
    !

    !these need to be deallocated anyhow
    DEALLOCATE(big_result)
    !
#if defined (__CUDA)
    DEALLOCATE(facb_d)
#endif
    DEALLOCATE(fac, facb )

    IF (noncolin) THEN
       DEALLOCATE(temppsic_nc)
    ELSE
       DEALLOCATE(temppsic)
    ENDIF

#if defined (__CUDA)
    IF (noncolin) THEN
       DEALLOCATE(result_nc_d)
    ELSE
       DEALLOCATE(result_d)
    ENDIF
#endif
    IF (noncolin) THEN
       DEALLOCATE( result_nc )
    ELSE
       DEALLOCATE( result )
    ENDIF
    !
    IF(okvan) DEALLOCATE( deexx)
    CALL stop_clock( 'vexx_k_fin' )
    !
    !------------------------------------------------------------------------
  END SUBROUTINE vexx_k
  !-----------------------------------------------------------------------
  !
  !-----------------------------------------------------------------------
  FUNCTION exxenergy ()
    !-----------------------------------------------------------------------
    !
    ! NB: This function is meant to give the SAME RESULT as exxenergy2.
    !     It is worth keeping it in the repository because in spite of being
    !     slower it is a simple driver using vexx potential routine so it is
    !     good, from time to time, to replace exxenergy2 with it to check that
    !     everything is ok and energy and potential are consistent as they should.
    !
    USE io_files,               ONLY : iunwfc_exx, nwordwfc
    USE buffers,                ONLY : get_buffer
    USE wvfct,                  ONLY : nbnd, npwx, wg, current_k
    USE gvect,                  ONLY : gstart
    USE wavefunctions,   ONLY : evc
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

    USE wavefunctions_gpum, ONLY : using_evc
    USE uspp_gpum,                 ONLY : using_vkb 

    IMPLICIT NONE

    TYPE(bec_type) :: becpsi
    REAL(DP)       :: exxenergy,  energy
    INTEGER        :: npw, ibnd, ik
    COMPLEX(DP)    :: vxpsi ( npwx*npol, nbnd ), psi(npwx*npol,nbnd)
    COMPLEX(DP),EXTERNAL :: zdotc
    !
    exxenergy=0._dp

    CALL start_clock ('exxenergy')

    IF(okvan) CALL allocate_bec_type( nkb, nbnd, becpsi)
    energy = 0._dp

    CALL using_evc(0)

    DO ik=1,nks
       npw = ngk (ik)
       ! setup variables for usage by vexx (same logic as for H_psi)
       current_k = ik
       IF ( lsda ) current_spin = isk(ik)
       ! end setup
       IF ( nks > 1 ) THEN
          CALL get_buffer(psi, nwordwfc_exx, iunwfc_exx, ik)
       ELSE
          psi(1:npwx*npol,1:nbnd) = evc(1:npwx*npol,1:nbnd)
       ENDIF
       !
       IF(okvan)THEN
          CALL using_vkb(1)
          ! prepare the |beta> function at k+q
          CALL init_us_2(npw, igk_exx(1,ik), xk(:,ik), vkb)
          ! compute <beta_I|psi_j> at this k+q point, for all band and all projectors
          CALL calbec(npw, vkb, psi, becpsi, nbnd)
       ENDIF
       !
       vxpsi(:,:) = (0._dp, 0._dp)
       CALL vexx(npwx,npw,nbnd,psi,vxpsi,becpsi)
       !
       DO ibnd=1,nbnd
          energy = energy + dble(wg(ibnd,ik) * zdotc(npw,psi(1,ibnd),1,vxpsi(1,ibnd),1))
          IF (noncolin) energy = energy + &
                            dble(wg(ibnd,ik) * zdotc(npw,psi(npwx+1,ibnd),1,vxpsi(npwx+1,ibnd),1))
          !
       ENDDO
       IF (gamma_only .and. gstart == 2) THEN
           DO ibnd=1,nbnd
              energy = energy - &
                       dble(0.5_dp * wg(ibnd,ik) * conjg(psi(1,ibnd)) * vxpsi(1,ibnd))
           ENDDO
       ENDIF
    ENDDO
    !
    IF (gamma_only) energy = 2 * energy

    CALL mp_sum( energy, intra_egrp_comm)
    CALL mp_sum( energy, inter_pool_comm )
    IF(okvan)  CALL deallocate_bec_type(becpsi)
    !
    exxenergy = energy
    !
    CALL stop_clock ('exxenergy')
    !-----------------------------------------------------------------------
  END FUNCTION exxenergy
  !-----------------------------------------------------------------------
  !
  !-----------------------------------------------------------------------
  FUNCTION exxenergy2()
    !-----------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    REAL(DP)   :: exxenergy2
    !
    CALL start_clock ('exxenergy')
    !
    IF( gamma_only ) THEN
       exxenergy2 = exxenergy2_gamma()
    ELSE
       exxenergy2 = exxenergy2_k()
    ENDIF
    !
    CALL stop_clock ('exxenergy')
    !
    !-----------------------------------------------------------------------
  END FUNCTION  exxenergy2
  !-----------------------------------------------------------------------
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
    USE wavefunctions,    ONLY : evc
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
    USE exx_band,                ONLY : igk_exx, &
         change_data_structure, transform_evc_to_exx, nwordwfc_exx,  &
         evc_exx
    !
    IMPLICIT NONE
    !
    REAL(DP)   :: exxenergy2_gamma
    !
    ! local variables
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
    CALL init_index_over_band(inter_egrp_comm,nbnd,nbnd)
    !
    CALL transform_evc_to_exx(0)
    !
    ialloc = nibands(my_egrp_id+1)
    !
    nrxxs= dfftt%nnr
    ALLOCATE( fac(dfftt%ngm) )
    !
    ALLOCATE(temppsic(nrxxs), temppsic_dble(nrxxs),temppsic_aimag(nrxxs))
    ALLOCATE( rhoc(nrxxs) )
    ALLOCATE( vkb_exx(npwx,nkb) )
    !
    energy=0.0_DP
    !
    CALL allocate_bec_type( nkb, nbnd, becpsi)
    !
    IKK_LOOP : &
    DO ikk=1,nks
       current_ik = global_kpoint_index ( nkstot, ikk )
       xkp = xk(:,ikk)
       !
       IF ( lsda ) current_spin = isk(ikk)
       npw = ngk (ikk)
       IF ( nks > 1 ) &
          CALL get_buffer (evc_exx, nwordwfc_exx, iunwfc_exx, ikk)
       !
       ! prepare the |beta> function at k+q
       CALL init_us_2(npw, igk_exx(:,ikk), xkp, vkb_exx)
       ! compute <beta_I|psi_j> at this k+q point, for all band and all projectors
       calbec_start = ibands(1,my_egrp_id+1)
       calbec_end = ibands(nibands(my_egrp_id+1),my_egrp_id+1)
       !
       intra_bgrp_comm_ = intra_bgrp_comm
       intra_bgrp_comm = intra_egrp_comm
       CALL calbec(npw, vkb_exx, evc_exx, &
            becpsi, nibands(my_egrp_id+1) )
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
          CALL g2_convolution_all(dfftt%ngm, gt, xkp, xkq, iq, &
               current_ik)
          fac = coulomb_fac(:,iq,current_ik)
          fac(gstart_t:) = 2 * coulomb_fac(gstart_t:,iq,current_ik)
          IF ( okvan .and..not.tqr ) CALL qvan_init (dfftt%ngm, xkq, xkp)
          !
          DO iegrp=1, negrp
             !
             ! compute the id of group whose data is currently worked on
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
                IF (jbnd.eq.0.or.jbnd.gt.nbnd) CYCLE
                !
                IF ( mod(ii,2)==1 ) THEN
                   !
                   temppsic = (0._DP,0._DP)
                   !
                   IF ( (ii+1)<=nibands(my_egrp_id+1) ) THEN
                      ! deal with double bands
!$omp parallel do  default(shared), private(ig)
                      DO ig = 1, npwt
                         temppsic( dfftt%nl(ig) )  = &
                              evc_exx(ig,ii) + (0._DP,1._DP) * evc_exx(ig,ii+1)
                         temppsic( dfftt%nlm(ig) ) = &
                              conjg(evc_exx(ig,ii) - (0._DP,1._DP) * evc_exx(ig,ii+1))
                      ENDDO
!$omp end parallel do
                   ENDIF
                   !
                   IF ( ii==nibands(my_egrp_id+1) ) THEN
                      ! deal with a single last band
!$omp parallel do  default(shared), private(ig)
                      DO ig = 1, npwt
                         temppsic( dfftt%nl(ig) )  =       evc_exx(ig,ii)
                         temppsic( dfftt%nlm(ig) ) = conjg(evc_exx(ig,ii))
                      ENDDO
!$omp end parallel do
                   ENDIF
                   !
                   CALL invfft ('Wave', temppsic, dfftt)
!$omp parallel do default(shared), private(ir)
                   DO ir = 1, nrxxs
                      temppsic_dble(ir)  = dble ( temppsic(ir) )
                      temppsic_aimag(ir) = aimag( temppsic(ir) )
                   ENDDO
!$omp end parallel do
                   !
                ENDIF
                !
                !determine which j-bands to calculate
                istart = 0
                iend = 0
                DO ipair=1, max_pairs
                   IF(egrp_pairs(1,ipair,my_egrp_id+1).eq.jbnd)THEN
                      IF(istart.eq.0)THEN
                         istart = egrp_pairs(2,ipair,my_egrp_id+1)
                         iend = istart
                      ELSE
                         iend = egrp_pairs(2,ipair,my_egrp_id+1)
                      END IF
                   END IF
                END DO
                !
                istart = max(istart,jblock_start)
                iend = min(iend,jblock_end)
                !
                IF(mod(istart,2)==0) THEN
                   ibnd_loop_start=istart-1
                ELSE
                   ibnd_loop_start=istart
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
                   IF( mod(ii,2) == 0 ) THEN
                      rhoc = 0.0_DP
!$omp parallel do default(shared), private(ir)
                      DO ir = 1, nrxxs
                         rhoc(ir) = exxbuff(ir,exxbuff_index,ikq) * temppsic_aimag(ir) / omega
                      ENDDO
!$omp end parallel do
                   ELSE
!$omp parallel do default(shared), private(ir)
                      DO ir = 1, nrxxs
                         rhoc(ir) = exxbuff(ir,exxbuff_index,ikq) * temppsic_dble(ir) / omega
                      ENDDO
!$omp end parallel do
                   ENDIF
                   !
                   IF(okvan .and.tqr) THEN
                      IF(ibnd>=istart) &
                           CALL addusxx_r(rhoc,_CX(becxx(ikq)%r(:,ibnd)),&
                                          _CX(becpsi%r(:,jbnd)))
                      IF(ibnd<iend) &
                           CALL addusxx_r(rhoc,_CY(becxx(ikq)%r(:,ibnd+1)),&
                                          _CX(becpsi%r(:,jbnd)))
                   ENDIF
                   !
                   ! bring rhoc to G-space
                   CALL fwfft ('Rho', rhoc, dfftt)
                   !
                   IF(okvan .and..not.tqr) THEN
                      IF(ibnd>=istart ) &
                           CALL addusxx_g( dfftt, rhoc, xkq, xkp, 'r', &
                           becphi_r=becxx(ikq)%r(:,ibnd), becpsi_r=becpsi%r(:,jbnd-calbec_start+1) )
                      IF(ibnd<iend) &
                           CALL addusxx_g( dfftt, rhoc, xkq, xkp, 'i', &
                           becphi_r=becxx(ikq)%r(:,ibnd+1), becpsi_r=becpsi%r(:,jbnd-calbec_start+1) )
                   ENDIF
                   !
                   vc = 0.0_DP
!$omp parallel do  default(shared), private(ig),  reduction(+:vc)
                   DO ig = 1,dfftt%ngm
                      !
                      ! The real part of rhoc contains the contribution from band ibnd
                      ! The imaginary part    contains the contribution from band ibnd+1
                      !
                      vc = vc + fac(ig) * ( x1 * &
                           abs( rhoc(dfftt%nl(ig)) + conjg(rhoc(dfftt%nlm(ig))) )**2 &
                                 +x2 * &
                           abs( rhoc(dfftt%nl(ig)) - conjg(rhoc(dfftt%nlm(ig))) )**2 )
                   ENDDO
!$omp end parallel do
                   !
                   vc = vc * omega * 0.25_DP / nqs
                   energy = energy - exxalfa * vc * wg(jbnd,ikk)
                   !
                   IF(okpaw) THEN
                      IF(ibnd>=ibnd_start) &
                           energy = energy +exxalfa*wg(jbnd,ikk)*&
                           x1 * PAW_xx_energy(_CX(becxx(ikq)%r(:,ibnd)),_CX(becpsi%r(:,jbnd)) )
                      IF(ibnd<ibnd_end) &
                           energy = energy +exxalfa*wg(jbnd,ikk)*&
                           x2 * PAW_xx_energy(_CX(becxx(ikq)%r(:,ibnd+1)), _CX(becpsi%r(:,jbnd)) )
                   ENDIF
                   !
                ENDDO &
                IBND_LOOP_GAM
             ENDDO &
             JBND_LOOP
             !
             ! get the next nbnd/negrp data
             IF (negrp>1) call mp_circular_shift_left( exxbuff(:,:,ikq), me_egrp, inter_egrp_comm )
             !
          ENDDO ! iegrp
          IF ( okvan .and..not.tqr ) CALL qvan_clean ( )
          !
       ENDDO &
       IQ_LOOP
    ENDDO &
    IKK_LOOP
    !
    DEALLOCATE(temppsic,temppsic_dble,temppsic_aimag)
    !
    DEALLOCATE(rhoc, fac )
    CALL deallocate_bec_type(becpsi)
    !
    CALL mp_sum( energy, inter_egrp_comm )
    CALL mp_sum( energy, intra_egrp_comm )
    CALL mp_sum( energy, inter_pool_comm )
    !
    exxenergy2_gamma = energy
    CALL change_data_structure(.FALSE.)
    !
    !-----------------------------------------------------------------------
  END FUNCTION  exxenergy2_gamma
  !-----------------------------------------------------------------------
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
    USE wavefunctions,    ONLY : evc
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
         transform_evc_to_exx, nwordwfc_exx, igk_exx, evc_exx
    !
    IMPLICIT NONE
    !
    REAL(DP)   :: exxenergy2_k
    !
    ! local variables
    REAL(DP) :: energy
    COMPLEX(DP), ALLOCATABLE :: temppsic(:,:)
    COMPLEX(DP), ALLOCATABLE :: temppsic_nc(:,:,:)
    COMPLEX(DP), ALLOCATABLE,TARGET :: rhoc(:,:)
#if defined(__USE_MANY_FFT)
    COMPLEX(DP), POINTER :: prhoc(:)
#endif
    REAL(DP),    ALLOCATABLE :: fac(:)
    INTEGER  :: npw, jbnd, ibnd, ibnd_inner_start, ibnd_inner_end, ibnd_inner_count, ik, ikk, ig, ikq, iq, ir
    INTEGER  :: h_ibnd, nrxxs, current_ik, ibnd_loop_start, nblock, nrt, irt, ir_start, ir_end
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
    CALL init_index_over_band(inter_egrp_comm,nbnd,nbnd)
    !
    CALL transform_evc_to_exx(0)
    !
    ialloc = nibands(my_egrp_id+1)
    !
    nrxxs = dfftt%nnr
    ALLOCATE( fac(dfftt%ngm) )
    !
    IF (noncolin) THEN
       ALLOCATE(temppsic_nc(nrxxs,npol,ialloc))
    ELSE
       ALLOCATE(temppsic(nrxxs,ialloc))
    ENDIF
    !
    energy=0.0_DP
    !
    CALL allocate_bec_type( nkb, nbnd, becpsi)
    !
    !precompute that stuff
    omega_inv=1.0/omega
    !
    IKK_LOOP : &
    DO ikk=1,nks
       !
       current_ik = global_kpoint_index ( nkstot, ikk )
       xkp = xk(:,ikk)
       !
       IF ( lsda ) current_spin = isk(ikk)
       npw = ngk (ikk)
       IF ( nks > 1 ) &
          CALL get_buffer (evc_exx, nwordwfc_exx, iunwfc_exx, ikk)
       !
       ! compute <beta_I|psi_j> at this k+q point, for all band and all projectors
       intra_bgrp_comm_ = intra_bgrp_comm
       intra_bgrp_comm = intra_egrp_comm
       IF (okvan.or.okpaw) THEN
          CALL compute_becpsi (npw, igk_exx(:,ikk), xkp, evc_exx, &
               becpsi%k(:,ibands(1,my_egrp_id+1)) )
       END IF
       intra_bgrp_comm = intra_bgrp_comm_
       !
       ! precompute temppsic
       !
       IF (noncolin) THEN
          temppsic_nc = 0.0_DP
       ELSE
          temppsic    = 0.0_DP
       ENDIF
       DO ii=1, nibands(my_egrp_id+1)
          !
          jbnd = ibands(ii,my_egrp_id+1)
          !
          IF (jbnd.eq.0.or.jbnd.gt.nbnd) CYCLE
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
             CALL invfft ('Wave', temppsic_nc(:,1,ii), dfftt)
             CALL invfft ('Wave', temppsic_nc(:,2,ii), dfftt)
             !
          ELSE
!$omp parallel do default(shared), private(ig)
             DO ig = 1, npw
                temppsic(dfftt%nl(igk_exx(ig,ikk)),ii) = evc_exx(ig,ii)
             ENDDO
!$omp end parallel do
             !
             CALL invfft ('Wave', temppsic(:,ii), dfftt)
             !
          ENDIF
       END DO
       !
       IQ_LOOP : &
       DO iq = 1,nqs
          !
          ikq  = index_xkq(current_ik,iq)
          ik   = index_xk(ikq)
          !
          xkq = xkq_collect(:,ikq)
          !
          CALL g2_convolution_all(dfftt%ngm, gt, xkp, xkq, iq, ikk)
          IF ( okvan .and..not.tqr ) CALL qvan_init (dfftt%ngm, xkq, xkp)
          !
          DO iegrp=1, negrp
             !
             ! compute the id of group whose data is currently worked on
             wegrp = MOD(iegrp+my_egrp_id-1, negrp)+1
             njt = (all_end(wegrp)-all_start(wegrp)+jblock)/jblock
             !
             IJT_LOOP : &
             DO ijt=1, njt
                !
                jblock_start = (ijt - 1) * jblock + all_start(wegrp)
                jblock_end = min(jblock_start+jblock-1,all_end(wegrp))
                !
                JBND_LOOP : &
                DO ii=1, nibands(my_egrp_id+1)
                   !
                   jbnd = ibands(ii,my_egrp_id+1)
                   !
                   IF (jbnd.eq.0.or.jbnd.gt.nbnd) CYCLE
                   !
                   !determine which j-bands to calculate
                   jstart = 0
                   jend = 0
                   DO ipair=1, max_pairs
                      IF(egrp_pairs(1,ipair,my_egrp_id+1).eq.jbnd)THEN
                         IF(jstart.eq.0)THEN
                            jstart = egrp_pairs(2,ipair,my_egrp_id+1)
                            jend = jstart
                         ELSE
                            jend = egrp_pairs(2,ipair,my_egrp_id+1)
                         END IF
                      END IF
                   END DO
                   !
                   !these variables prepare for inner band parallelism
                   jstart = max(jstart,jblock_start)
                   jend = min(jend,jblock_end)
                   ibnd_inner_start=jstart
                   ibnd_inner_end=jend
                   ibnd_inner_count=jend-jstart+1
                   !
                   !allocate arrays
                   ALLOCATE( rhoc(nrxxs,ibnd_inner_count) )
#if defined(__USE_MANY_FFT)
                   prhoc(1:nrxxs*ibnd_inner_count) => rhoc
#endif 
                   !calculate rho in real space
                   nblock=2048
                   nrt = (nrxxs+nblock-1) / nblock
!$omp parallel do collapse(2) private(ir_start,ir_end)
                   DO irt = 1, nrt
                      DO ibnd=ibnd_inner_start, ibnd_inner_end
                         ir_start = (irt - 1) * nblock + 1
                         ir_end = min(ir_start+nblock-1,nrxxs)
                         IF (noncolin) THEN
                            DO ir = ir_start, ir_end
                               rhoc(ir,ibnd-ibnd_inner_start+1) = &
                                 ( conjg(exxbuff(ir,ibnd-all_start(wegrp)+iexx_start,ikq)) * temppsic_nc(ir,1,ii) + &
                                 conjg(exxbuff(nrxxs+ir,ibnd-all_start(wegrp)+iexx_start,ikq))*temppsic_nc(ir,2,ii) ) * omega_inv
                            ENDDO
                         ELSE
                            DO ir = ir_start, ir_end
                               rhoc(ir,ibnd-ibnd_inner_start+1) = omega_inv * &
                                 conjg(exxbuff(ir,ibnd-all_start(wegrp)+iexx_start,ikq)) * &
                                 temppsic(ir,ii)
                            ENDDO
                         ENDIF
                      ENDDO
                   ENDDO
!$omp end parallel do
                   !
                   ! augment the "charge" in real space
                   IF(okvan .and. tqr) THEN
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
                   DO ibnd=ibnd_inner_start, ibnd_inner_end
                      CALL fwfft('Rho', rhoc(:,ibnd-ibnd_inner_start+1), dfftt)
                   ENDDO
#endif
                   ! augment the "charge" in G space
                   IF(okvan .and. .not. tqr) THEN
                      DO ibnd = ibnd_inner_start, ibnd_inner_end
                         CALL addusxx_g(dfftt, rhoc(:,ibnd-ibnd_inner_start+1), &
                              xkq, xkp, 'c', becphi_c=becxx(ikq)%k(:,ibnd), &
                              becpsi_c=becpsi%k(:,jbnd))
                      ENDDO
                   ENDIF
                   !
!$omp parallel do reduction(+:energy) private(vc)
                   DO ibnd = ibnd_inner_start, ibnd_inner_end
                      vc=0.0_DP
                      DO ig=1,dfftt%ngm
                         vc = vc + coulomb_fac(ig,iq,ikk) * &
                             dble(rhoc(dfftt%nl(ig),ibnd-ibnd_inner_start+1) *&
                             conjg(rhoc(dfftt%nl(ig),ibnd-ibnd_inner_start+1)))
                      ENDDO
                      vc = vc * omega * x_occupation(ibnd,ik) / nqs
                      energy = energy - exxalfa * vc * wg(jbnd,ikk)
                      
                      IF(okpaw) THEN
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
             IF (negrp>1) call mp_circular_shift_left( exxbuff(:,:,ikq), me_egrp, inter_egrp_comm )
             !
          END DO !iegrp
          !
          IF ( okvan .and..not.tqr ) CALL qvan_clean ( )
       ENDDO &
       IQ_LOOP
       !
    ENDDO &
    IKK_LOOP
    !
    IF (noncolin) THEN
       DEALLOCATE(temppsic_nc)
    ELSE
       DEALLOCATE(temppsic)
    ENDIF
    !
    DEALLOCATE(fac)
    CALL deallocate_bec_type(becpsi)
    !
    CALL mp_sum( energy, inter_egrp_comm )
    CALL mp_sum( energy, intra_egrp_comm )
    CALL mp_sum( energy, inter_pool_comm )
    !
    exxenergy2_k = energy
    CALL change_data_structure(.FALSE.)
    !
    !-----------------------------------------------------------------------
  END FUNCTION  exxenergy2_k
  !-----------------------------------------------------------------------
  !
  !-----------------------------------------------------------------------
  FUNCTION exx_stress()
    !-----------------------------------------------------------------------
    !
    ! This is Eq.(10) of PRB 73, 125120 (2006).
    !
    USE constants,            ONLY : fpi, e2, pi, tpi
    USE io_files,             ONLY : iunwfc_exx, nwordwfc
    USE buffers,              ONLY : get_buffer
    USE cell_base,            ONLY : alat, omega, bg, at, tpiba
    USE symm_base,            ONLY : nsym, s
    USE wvfct,                ONLY : nbnd, npwx, wg, current_k
    USE wavefunctions, ONLY : evc
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
    USE exx_base,             ONLY : nq1, nq2, nq3, nqs, eps, exxdiv, &
         x_gamma_extrapolation, on_double_grid, grid_factor, yukawa,  &
         erfc_scrlen, use_coulomb_vcut_ws, use_coulomb_vcut_spheric,  &
         gau_scrlen, vcut, index_xkq, index_xk, index_sym
    USE exx_band,             ONLY : change_data_structure, &
         transform_evc_to_exx, g_exx, igk_exx,  &
         nwordwfc_exx, evc_exx
    USE coulomb_vcut_module,  ONLY : vcut_get,  vcut_spheric_get
    !
    ! ---- local variables -------------------------------------------------
    !
    IMPLICIT NONE
    !
    ! local variables
    REAL(DP)   :: exx_stress(3,3), exx_stress_(3,3)
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

    CALL start_clock ('exx_stress')

    CALL transform_evc_to_exx(0)

    IF (npool>1) CALL errore('exx_stress','stress not available with pools',1)
    IF (noncolin) CALL errore('exx_stress','noncolinear stress not implemented',1)
    IF (okvan) CALL infomsg('exx_stress','USPP stress not tested')

    nrxxs = dfftt%nnr
    ngm   = dfftt%ngm
    delta = reshape( (/1._dp,0._dp,0._dp, 0._dp,1._dp,0._dp, 0._dp,0._dp,1._dp/), (/3,3/))
    exx_stress_ = 0._dp
    ALLOCATE( tempphic(nrxxs), temppsic(nrxxs), rhoc(nrxxs), fac(ngm) )
    ALLOCATE( fac_tens(3,3,ngm), fac_stress(ngm) )
    !
    nqi=nqs
    !
    ! loop over k-points
    DO ikk = 1, nks
        current_k = ikk
        IF (lsda) current_spin = isk(ikk)
        npw = ngk(ikk)

        IF (nks > 1) &
            CALL get_buffer(evc_exx, nwordwfc_exx, iunwfc_exx, ikk)

            DO iqi = 1, nqi
                !
                iq=iqi
                !
                ikq  = index_xkq(current_k,iq)
                ik   = index_xk(ikq)
                isym = abs(index_sym(ikq))
                !

                ! FIXME: use cryst_to_cart and company as above..
                xk_cryst(:)=at(1,:)*xk(1,ik)+at(2,:)*xk(2,ik)+at(3,:)*xk(3,ik)
                IF (index_sym(ikq) < 0) xk_cryst = -xk_cryst
                sxk(:) = s(:,1,isym)*xk_cryst(1) + &
                         s(:,2,isym)*xk_cryst(2) + &
                         s(:,3,isym)*xk_cryst(3)
                xkq(:) = bg(:,1)*sxk(1) + bg(:,2)*sxk(2) + bg(:,3)*sxk(3)

                !CALL start_clock ('exxen2_ngmloop')

!$omp parallel do default(shared), private(ig, beta, q, qq, on_double_grid, x)
                DO ig = 1, ngm
                   IF(negrp.eq.1) THEN
                      q(1)= xk(1,current_k) - xkq(1) + g(1,ig)
                      q(2)= xk(2,current_k) - xkq(2) + g(2,ig)
                      q(3)= xk(3,current_k) - xkq(3) + g(3,ig)
                   ELSE
                      q(1)= xk(1,current_k) - xkq(1) + g_exx(1,ig)
                      q(2)= xk(2,current_k) - xkq(2) + g_exx(2,ig)
                      q(3)= xk(3,current_k) - xkq(3) + g_exx(3,ig)
                   END IF

                  q = q * tpiba
                  qq = ( q(1)*q(1) + q(2)*q(2) + q(3)*q(3) )

                  DO beta = 1, 3
                      fac_tens(1:3,beta,ig) = q(1:3)*q(beta)
                  ENDDO

                  IF (x_gamma_extrapolation) THEN
                      on_double_grid = .true.
                      x= 0.5d0/tpiba*(q(1)*at(1,1)+q(2)*at(2,1)+q(3)*at(3,1))*nq1
                      on_double_grid = on_double_grid .and. (abs(x-nint(x))<eps)
                      x= 0.5d0/tpiba*(q(1)*at(1,2)+q(2)*at(2,2)+q(3)*at(3,2))*nq2
                      on_double_grid = on_double_grid .and. (abs(x-nint(x))<eps)
                      x= 0.5d0/tpiba*(q(1)*at(1,3)+q(2)*at(2,3)+q(3)*at(3,3))*nq3
                      on_double_grid = on_double_grid .and. (abs(x-nint(x))<eps)
                  ELSE
                      on_double_grid = .false.
                  ENDIF

                  IF (use_coulomb_vcut_ws) THEN
                      fac(ig) = vcut_get(vcut, q)
                      fac_stress(ig) = 0._dp   ! not implemented
                      IF (gamma_only .and. qq > 1.d-8) fac(ig) = 2.d0 * fac(ig)

                  ELSEIF ( use_coulomb_vcut_spheric ) THEN
                      fac(ig) = vcut_spheric_get(vcut, q)
                      fac_stress(ig) = 0._dp   ! not implemented
                      IF (gamma_only .and. qq > 1.d-8) fac(ig) = 2.d0 * fac(ig)

                  ELSE IF (gau_scrlen > 0) THEN
                      fac(ig)=e2*((pi/gau_scrlen)**(1.5d0))* &
                            exp(-qq/4.d0/gau_scrlen) * grid_factor
                      fac_stress(ig) =  e2*2.d0/4.d0/gau_scrlen * &
                            exp(-qq/4.d0/gau_scrlen) *((pi/gau_scrlen)**(1.5d0))* &
                                                                     grid_factor
                      IF (gamma_only) fac(ig) = 2.d0 * fac(ig)
                      IF (gamma_only) fac_stress(ig) = 2.d0 * fac_stress(ig)
                      IF (on_double_grid) fac(ig) = 0._dp
                      IF (on_double_grid) fac_stress(ig) = 0._dp

                  ELSEIF (qq > 1.d-8) THEN
                      IF ( erfc_scrlen > 0 ) THEN
                        fac(ig)=e2*fpi/qq*(1._dp-exp(-qq/4.d0/erfc_scrlen**2)) * grid_factor
                        fac_stress(ig) = -e2*fpi * 2.d0/qq**2 * ( &
                            (1._dp+qq/4.d0/erfc_scrlen**2)*exp(-qq/4.d0/erfc_scrlen**2) - 1._dp) * &
                            grid_factor
                      ELSE
                        fac(ig)=e2*fpi/( qq + yukawa ) * grid_factor
                        fac_stress(ig) = 2.d0 * e2*fpi/(qq+yukawa)**2 * grid_factor
                      ENDIF

                      IF (gamma_only) fac(ig) = 2.d0 * fac(ig)
                      IF (gamma_only) fac_stress(ig) = 2.d0 * fac_stress(ig)
                      IF (on_double_grid) fac(ig) = 0._dp
                      IF (on_double_grid) fac_stress(ig) = 0._dp

                  ELSE
                      fac(ig)= -exxdiv ! or rather something else (see f.gygi)
                      fac_stress(ig) = 0._dp  ! or -exxdiv_stress (not yet implemented)
                      IF ( yukawa> 0._dp .and. .not. x_gamma_extrapolation) THEN
                        fac(ig) = fac(ig) + e2*fpi/( qq + yukawa )
                        fac_stress(ig) = 2.d0 * e2*fpi/(qq+yukawa)**2
                      ENDIF
                      IF (erfc_scrlen > 0._dp .and. .not. x_gamma_extrapolation) THEN
                        fac(ig) = e2*fpi / (4.d0*erfc_scrlen**2)
                        fac_stress(ig) = e2*fpi / (8.d0*erfc_scrlen**4)
                      ENDIF
                  ENDIF
                ENDDO
!$omp end parallel do
                !CALL stop_clock ('exxen2_ngmloop')
        DO iegrp=1, negrp
           !
           ! compute the id of group whose data is currently worked on
           wegrp = MOD(iegrp+my_egrp_id-1, negrp)+1
           !
           jblock_start = all_start(wegrp)
           jblock_end   = all_end(wegrp)
           !
        ! loop over bands
        DO ii = 1, nibands(my_egrp_id+1)
            !
            jbnd = ibands(ii,my_egrp_id+1)
            !
            IF (jbnd.eq.0.or.jbnd.gt.nbnd) CYCLE
            !
            !determine which j-bands to calculate
            jstart = 0
            jend = 0
            DO ipair=1, max_pairs
               IF(egrp_pairs(1,ipair,my_egrp_id+1).eq.jbnd)THEN
                  IF(jstart.eq.0)THEN
                     jstart = egrp_pairs(2,ipair,my_egrp_id+1)
                     jend = jstart
                  ELSE
                     jend = egrp_pairs(2,ipair,my_egrp_id+1)
                  END IF
               END IF
            ENDDO
            !
            jstart = max(jstart,jblock_start)
            jend = min(jend,jblock_end)
            !
            temppsic(:) = ( 0._dp, 0._dp )
!$omp parallel do default(shared), private(ig)
            DO ig = 1, npw
                temppsic(dfftt%nl(igk_exx(ig,ikk))) = evc_exx(ig,ii)
            ENDDO
!$omp end parallel do
            !
            IF(gamma_only) THEN
!$omp parallel do default(shared), private(ig)
                DO ig = 1, npw
                    temppsic(dfftt%nlm(igk_exx(ig,ikk))) = &
                         conjg(evc_exx(ig,ii))
                ENDDO
!$omp end parallel do
            ENDIF

            CALL invfft ('Wave', temppsic, dfftt)

                IF (gamma_only) THEN
                    !
                    IF(mod(jstart,2)==0) THEN
                      ibnd_loop_start=jstart-1
                    ELSE
                      ibnd_loop_start=jstart
                    ENDIF
                    !
                    DO ibnd = ibnd_loop_start, jend, 2     !for each band of psi
                        !
                        exxbuff_index = (ibnd+1)/2-(all_start(wegrp)+1)/2+(iexx_start+1)/2
                        !
                        IF( ibnd < jstart ) THEN
                            x1 = 0._dp
                        ELSE
                            x1 = x_occupation(ibnd,  ik)
                        ENDIF

                        IF( ibnd == jend) THEN
                            x2 = 0._dp
                        ELSE
                            x2 = x_occupation(ibnd+1,  ik)
                        ENDIF
                        IF ( abs(x1) < eps_occ .and. abs(x2) < eps_occ ) CYCLE
                        !
                        ! calculate rho in real space
!$omp parallel do default(shared), private(ir)
                        DO ir = 1, nrxxs
                            tempphic(ir) = exxbuff(ir,exxbuff_index,ikq)
                            rhoc(ir)     = conjg(tempphic(ir))*temppsic(ir) / omega
                        ENDDO
!$omp end parallel do
                        ! bring it to G-space
                        CALL fwfft ('Rho', rhoc, dfftt)

                        vc = 0._dp
!$omp parallel do default(shared), private(ig), reduction(+:vc)
                        DO ig = 1, ngm
                            !
                            vc(:,:) = vc(:,:) + x1 * 0.25_dp * &
                                      abs( rhoc(dfftt%nl(ig)) + &
                                      conjg(rhoc(dfftt%nlm(ig))))**2 * &
                                      (fac_tens(:,:,ig)*fac_stress(ig)/2.d0 - delta(:,:)*fac(ig))
                            vc(:,:) = vc(:,:) + x2 * 0.25_dp * &
                                      abs( rhoc(dfftt%nl(ig)) - &
                                      conjg(rhoc(dfftt%nlm(ig))))**2 * &
                                      (fac_tens(:,:,ig)*fac_stress(ig)/2.d0 - delta(:,:)*fac(ig))
                        ENDDO
!$omp end parallel do
                        vc = vc / nqs / 4.d0
                        exx_stress_ = exx_stress_ + exxalfa * vc * wg(jbnd,ikk)
                    ENDDO

                ELSE

                    DO ibnd = jstart, jend    !for each band of psi
                      !
                      IF ( abs(x_occupation(ibnd,ik)) < 1.d-6) CYCLE
                      !
                      ! calculate rho in real space
!$omp parallel do default(shared), private(ir)
                      DO ir = 1, nrxxs
                          tempphic(ir) = exxbuff(ir,ibnd-all_start(wegrp)+iexx_start,ikq)
                          rhoc(ir)     = conjg(tempphic(ir))*temppsic(ir) / omega
                      ENDDO
!$omp end parallel do

                      ! bring it to G-space
                      CALL fwfft ('Rho', rhoc, dfftt)

                      vc = 0._dp
!$omp parallel do default(shared), private(ig), reduction(+:vc)
                      DO ig = 1, ngm
                          vc(:,:) = vc(:,:) + rhoc(dfftt%nl(ig))  * &
                                        conjg(rhoc(dfftt%nl(ig)))* &
                                    (fac_tens(:,:,ig)*fac_stress(ig)/2.d0 - delta(:,:)*fac(ig))
                      ENDDO
!$omp end parallel do
                      vc = vc * x_occupation(ibnd,ik) / nqs / 4.d0
                      exx_stress_ = exx_stress_ + exxalfa * vc * wg(jbnd,ikk)

                    ENDDO

                ENDIF ! gamma or k-points

        ENDDO ! jbnd
           !
           ! get the next nbnd/negrp data
           IF (negrp>1) call mp_circular_shift_left( exxbuff(:,:,ikq), me_egrp, inter_egrp_comm )
           !
        ENDDO ! iegrp
            ENDDO ! iqi
    ENDDO ! ikk

    DEALLOCATE(tempphic, temppsic, rhoc, fac, fac_tens, fac_stress )
    !
    CALL mp_sum( exx_stress_, intra_egrp_comm )
    CALL mp_sum( exx_stress_, inter_egrp_comm )
    CALL mp_sum( exx_stress_, inter_pool_comm )
    exx_stress = exx_stress_

    CALL change_data_structure(.FALSE.)

    CALL stop_clock ('exx_stress')
    !-----------------------------------------------------------------------
  END FUNCTION exx_stress
  !-----------------------------------------------------------------------
  !
!----------------------------------------------------------------------
SUBROUTINE compute_becpsi (npw_, igk_, q_, evc_exx, becpsi_k)
  !----------------------------------------------------------------------
  !
  !   Calculates beta functions (Kleinman-Bylander projectors), with
  !   structure factor, for all atoms, in reciprocal space. On input:
  !      npw_       : number of PWs 
  !      igk_(npw_) : indices of G in the list of q+G vectors
  !      q_(3)      : q vector (2pi/a units)
  !  On output:
  !      vkb_(npwx,nkb) : beta functions (npw_ <= npwx)
  !
  USE kinds,      ONLY : DP
  USE ions_base,  ONLY : nat, ntyp => nsp, ityp, tau
  USE cell_base,  ONLY : tpiba
  USE constants,  ONLY : tpi
  USE gvect,      ONLY : eigts1, eigts2, eigts3, mill, g
  USE wvfct,      ONLY : npwx, nbnd
  USE us,         ONLY : nqx, dq, tab, tab_d2y, spline_ps
  USE m_gth,      ONLY : mk_ffnl_gth
  USE splinelib
  USE uspp,       ONLY : nkb, nhtol, nhtolm, indv
  USE uspp_param, ONLY : upf, lmaxkb, nhm, nh
  USE becmod,     ONLY : calbec
  USE mp_exx,     ONLY : ibands, nibands, my_egrp_id
  !
  USE us_gpum,    ONLY : using_tab, using_tab_d2y
  !
  implicit none
  !
  INTEGER, INTENT (IN) :: npw_, igk_ (npw_)
  REAL(dp), INTENT(IN) :: q_(3)
  COMPLEX(dp), INTENT(IN) :: evc_exx(npwx,nibands(my_egrp_id+1))
  COMPLEX(dp), INTENT(OUT) :: becpsi_k(nkb,nibands(my_egrp_id+1))
  COMPLEX(dp) :: vkb_ (npwx, 1)
  !
  !     Local variables
  !
  integer :: i0,i1,i2,i3, ig, lm, na, nt, nb, ih, jkb

  real(DP) :: px, ux, vx, wx, arg
  real(DP), allocatable :: gk (:,:), qg (:), vq (:), ylm (:,:), vkb1(:,:)

  complex(DP) :: phase, pref
  complex(DP), allocatable :: sk(:)

  real(DP), allocatable :: xdata(:)
  integer :: iq
  integer :: istart, iend

  istart = ibands(1,my_egrp_id+1)
  iend = ibands(nibands(my_egrp_id+1),my_egrp_id+1)
  !
  !
  if (lmaxkb.lt.0) return
  !
  call using_tab(0)
  if (spline_ps) call using_tab_d2y(0)
  !
  allocate (vkb1( npw_,nhm))    
  allocate (  sk( npw_))    
  allocate (  qg( npw_))    
  allocate (  vq( npw_))    
  allocate ( ylm( npw_, (lmaxkb + 1) **2))    
  allocate (  gk( 3, npw_))    
  !
!   write(*,'(3i4,i5,3f10.5)') size(tab,1), size(tab,2), size(tab,3), size(vq), q_

  do ig = 1, npw_
     gk (1,ig) = q_(1) + g(1, igk_(ig) )
     gk (2,ig) = q_(2) + g(2, igk_(ig) )
     gk (3,ig) = q_(3) + g(3, igk_(ig) )
     qg (ig) = gk(1, ig)**2 +  gk(2, ig)**2 + gk(3, ig)**2
  enddo
  !
  call ylmr2 ((lmaxkb+1)**2, npw_, gk, qg, ylm)
  !
  ! set now qg=|q+G| in atomic units
  !
  do ig = 1, npw_
     qg(ig) = sqrt(qg(ig))*tpiba
  enddo

  if (spline_ps) then
    allocate(xdata(nqx))
    do iq = 1, nqx
      xdata(iq) = (iq - 1) * dq
    enddo
  endif
  ! |beta_lm(q)> = (4pi/omega).Y_lm(q).f_l(q).(i^l).S(q)
  jkb = 0
  do nt = 1, ntyp
     ! calculate beta in G-space using an interpolation table f_l(q)=\int _0 ^\infty dr r^2 f_l(r) j_l(q.r)
     do nb = 1, upf(nt)%nbeta
        if ( upf(nt)%is_gth ) then
           call mk_ffnl_gth( nt, nb, npw_, qg, vq )
        else
           do ig = 1, npw_
              if (spline_ps) then
                vq(ig) = splint(xdata, tab(:,nb,nt), tab_d2y(:,nb,nt), qg(ig))
              else
                px = qg (ig) / dq - int (qg (ig) / dq)
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
              endif
           enddo
        endif
        ! add spherical harmonic part  (Y_lm(q)*f_l(q)) 
        do ih = 1, nh (nt)
           if (nb.eq.indv (ih, nt) ) then
              !l = nhtol (ih, nt)
              lm =nhtolm (ih, nt)
              do ig = 1, npw_
                 vkb1 (ig,ih) = ylm (ig, lm) * vq (ig)
              enddo
           endif
        enddo
     enddo
     !
     ! vkb1 contains all betas including angular part for type nt
     ! now add the structure factor and factor (-i)^l
     !
     do na = 1, nat
        ! ordering: first all betas for atoms of type 1
        !           then  all betas for atoms of type 2  and so on
        if (ityp (na) .eq.nt) then
           arg = (q_(1) * tau (1, na) + &
                  q_(2) * tau (2, na) + &
                  q_(3) * tau (3, na) ) * tpi
           phase = CMPLX(cos (arg), - sin (arg) ,kind=DP)
           do ig = 1, npw_
              sk (ig) = eigts1 (mill(1,igk_(ig)), na) * &
                        eigts2 (mill(2,igk_(ig)), na) * &
                        eigts3 (mill(3,igk_(ig)), na)
           enddo
           do ih = 1, nh (nt)
              jkb = jkb + 1
              pref = (0.d0, -1.d0) **nhtol (ih, nt) * phase
              do ig = 1, npw_
                 vkb_(ig, 1) = vkb1 (ig,ih) * sk (ig) * pref
              enddo
              do ig = npw_+1, npwx
                 vkb_(ig, 1) = (0.0_dp, 0.0_dp)
              enddo
              !
              CALL calbec(npw_, vkb_, evc_exx, becpsi_k(jkb:jkb,:), &
                   nibands(my_egrp_id+1))
              !
           enddo
        endif
     enddo
  enddo
  deallocate (gk)
  deallocate (ylm)
  deallocate (vq)
  deallocate (qg)
  deallocate (sk)
  deallocate (vkb1)

  return
END SUBROUTINE compute_becpsi
  !-----------------------------------------------------------------------

SUBROUTINE aceinit( exex )
  !
  USE wvfct,      ONLY : nbnd, npwx, current_k
  USE klist,      ONLY : nks, xk, ngk, igk_k
  USE uspp,       ONLY : nkb, vkb, okvan
  USE becmod,     ONLY : allocate_bec_type, deallocate_bec_type, &
                         bec_type, calbec
  USE lsda_mod,   ONLY : current_spin, lsda, isk
  USE io_files,   ONLY : nwordwfc, iunwfc
  USE buffers,    ONLY : get_buffer
  USE mp_pools,   ONLY : inter_pool_comm
  USE mp_bands,   ONLY : intra_bgrp_comm
  USE mp,         ONLY : mp_sum
  USE wavefunctions, ONLY : evc
  !
  USE wavefunctions_gpum, ONLY : using_evc
  USE uspp_gpum,                 ONLY : using_vkb
  !
  IMPLICIT NONE
  !
  REAL (DP) :: ee, eexx
  INTEGER :: ik, npw
  TYPE(bec_type) :: becpsi
  REAL (DP), OPTIONAL, INTENT(OUT) :: exex
  !
  IF(nbndproj<x_nbnd_occ.or.nbndproj>nbnd) THEN 
    write(stdout,'(3(A,I4))') ' occ = ', x_nbnd_occ, ' proj = ', nbndproj, ' tot = ', nbnd
    CALL errore('aceinit','n_proj must be between occ and tot.',1)
  END IF 

  CALL using_evc(0)

  IF (.not. allocated(xi)) ALLOCATE( xi(npwx*npol,nbndproj,nks) )
  IF ( okvan ) CALL allocate_bec_type( nkb, nbnd, becpsi)
  eexx = 0.0d0
  xi = (0.0d0,0.0d0)
  DO ik = 1, nks
     npw = ngk (ik)
     current_k = ik
     IF ( lsda ) current_spin = isk(ik)
     IF ( nks > 1 ) CALL get_buffer(evc, nwordwfc, iunwfc, ik)
     IF ( nks > 1 ) CALL using_evc(2)
     IF ( okvan ) THEN
        CALL using_vkb(1)
        CALL init_us_2(npw, igk_k(1,ik), xk(:,ik), vkb)
        CALL calbec ( npw, vkb, evc, becpsi, nbnd )
     ENDIF
     IF (gamma_only) THEN
        CALL aceinit_gamma(npw,nbnd,evc,xi(1,1,ik),becpsi,ee)
     ELSE
        CALL aceinit_k(npw,nbnd,evc,xi(1,1,ik),becpsi,ee)
     ENDIF
     eexx = eexx + ee
  ENDDO
  CALL mp_sum( eexx, inter_pool_comm)
  ! WRITE(stdout,'(/,5X,"ACE energy",f15.8)') eexx
  IF(present(exex)) exex = eexx
  IF ( okvan ) CALL deallocate_bec_type(becpsi)
  domat = .false.
END SUBROUTINE aceinit
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE aceinit_gamma(nnpw,nbnd,phi,xitmp,becpsi,exxe)
USE becmod,         ONLY : bec_type
USE lsda_mod,       ONLY : current_spin
USE mp,             ONLY : mp_stop
!
! compute xi(npw,nbndproj) for the ACE method
!
IMPLICIT NONE
  INTEGER :: nnpw,nbnd, nrxxs
  COMPLEX(DP) :: phi(nnpw,nbnd)
  real(DP), ALLOCATABLE :: mexx(:,:)
  COMPLEX(DP) :: xitmp(nnpw,nbndproj)
  real(DP) :: exxe
  real(DP), PARAMETER :: Zero=0.0d0, One=1.0d0, Two=2.0d0, Pt5=0.50d0
  TYPE(bec_type), INTENT(in) :: becpsi
  logical :: domat0

  CALL start_clock( 'aceinit' )

  nrxxs= dfftt%nnr * npol

  ALLOCATE( mexx(nbndproj,nbndproj) )
  xitmp = (Zero,Zero)
  mexx = Zero
  !
  IF( local_thr.gt.0.0d0 ) then  
    CALL vexx_loc(nnpw, nbndproj, xitmp, mexx)
    Call MatSymm('S','L',mexx,nbndproj)
  ELSE
!   |xi> = Vx[phi]|phi>
    CALL vexx(nnpw, nnpw, nbndproj, phi, xitmp, becpsi)
!   mexx = <phi|Vx[phi]|phi>
    CALL matcalc('exact',.true.,0,nnpw,nbndproj,nbndproj,phi,xitmp,mexx,exxe)
!   |xi> = -One * Vx[phi]|phi> * rmexx^T
  END IF

  CALL aceupdate(nbndproj,nnpw,xitmp,mexx)

  DEALLOCATE( mexx )

  IF( local_thr.gt.0.0d0 ) then  
    domat0 = domat
    domat = .true.
    CALL vexxace_gamma(nnpw,nbndproj,evc0(1,1,current_spin),exxe)
    evc0(:,:,current_spin) = phi(:,:)
    domat = domat0
  END IF 

  CALL stop_clock( 'aceinit' )

END SUBROUTINE aceinit_gamma
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE vexxace_gamma(nnpw,nbnd,phi,exxe,vphi)
USE wvfct,    ONLY : current_k, wg
USE lsda_mod, ONLY : current_spin
!
! do the ACE potential and
! (optional) print the ACE matrix representation
!
IMPLICIT NONE
  real(DP) :: exxe
  INTEGER :: nnpw,nbnd,i,ik
  COMPLEX(DP) :: phi(nnpw,nbnd)
  COMPLEX(DP),OPTIONAL :: vphi(nnpw,nbnd)
  real*8,ALLOCATABLE :: rmexx(:,:)
  COMPLEX(DP),ALLOCATABLE :: cmexx(:,:), vv(:,:)
  real*8, PARAMETER :: Zero=0.0d0, One=1.0d0, Two=2.0d0, Pt5=0.50d0

  CALL start_clock('vexxace')

  ALLOCATE( vv(nnpw,nbnd) )
  IF(present(vphi)) THEN
    vv = vphi
  ELSE
    vv = (Zero, Zero)
  ENDIF

! do the ACE potential
  ALLOCATE( rmexx(nbndproj,nbnd),cmexx(nbndproj,nbnd) )
  rmexx = Zero
  cmexx = (Zero,Zero)
! <xi|phi>
  CALL matcalc('<xi|phi>',.false.,0,nnpw,nbndproj,nbnd,xi(1,1,current_k),phi,rmexx,exxe)
! |vv> = |vphi> + (-One) * |xi> * <xi|phi>
  cmexx = (One,Zero)*rmexx
  CALL ZGEMM ('N','N',nnpw,nbnd,nbndproj,-(One,Zero),xi(1,1,current_k), &
                      nnpw,cmexx,nbndproj,(One,Zero),vv,nnpw)
  DEALLOCATE( cmexx,rmexx )

  IF(domat) THEN
    ALLOCATE( rmexx(nbnd,nbnd) )
    CALL matcalc('ACE',.true.,0,nnpw,nbnd,nbnd,phi,vv,rmexx,exxe)
    DEALLOCATE( rmexx )
#if defined(__DEBUG)
    WRITE(stdout,'(4(A,I3),A,I9,A,f12.6)') 'vexxace: nbnd=', nbnd, ' nbndproj=',nbndproj, &
                                              ' k=',current_k,' spin=',current_spin,' npw=',nnpw, ' E=',exxe
  ELSE
    WRITE(stdout,'(4(A,I3),A,I9)')         'vexxace: nbnd=', nbnd, ' nbndproj=',nbndproj, &
                                              ' k=',current_k,' spin=',current_spin,' npw=',nnpw
#endif
  ENDIF

  IF(present(vphi)) vphi = vv
  DEALLOCATE( vv )

  CALL stop_clock('vexxace')

END SUBROUTINE vexxace_gamma
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE aceupdate(nbndproj,nnpw,xitmp,rmexx)
!
!  Build the ACE operator from the potential amd matrix
!  (rmexx is assumed symmetric and only the Lower Triangular part is considered)
!
IMPLICIT NONE
  INTEGER :: nbndproj,nnpw
  REAL(DP) :: rmexx(nbndproj,nbndproj)
  COMPLEX(DP),ALLOCATABLE :: cmexx(:,:)
  COMPLEX(DP) ::  xitmp(nnpw,nbndproj)
  REAL(DP), PARAMETER :: Zero=0.0d0, One=1.0d0, Two=2.0d0, Pt5=0.50d0

  CALL start_clock('aceupdate')

! rmexx = -(Cholesky(rmexx))^-1
  rmexx = -rmexx
! CALL invchol( nbndproj, rmexx )
  CALL MatChol( nbndproj, rmexx )
  CALL MatInv( 'L',nbndproj, rmexx )

! |xi> = -One * Vx[phi]|phi> * rmexx^T
  ALLOCATE( cmexx(nbndproj,nbndproj) )
  cmexx = (One,Zero)*rmexx
  CALL ZTRMM('R','L','C','N',nnpw,nbndproj,(One,Zero),cmexx,nbndproj,xitmp,nnpw)
  DEALLOCATE( cmexx )

  CALL stop_clock('aceupdate')

END SUBROUTINE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE aceinit_k(nnpw,nbnd,phi,xitmp,becpsi,exxe)
USE becmod,               ONLY : bec_type
USE wvfct,                ONLY : current_k, npwx
USE noncollin_module,     ONLY : npol
!
! compute xi(npw,nbndproj) for the ACE method
!
IMPLICIT NONE
INTEGER :: nnpw,nbnd,i
COMPLEX(DP) :: phi(npwx*npol,nbnd),xitmp(npwx*npol,nbndproj)
COMPLEX(DP), ALLOCATABLE :: mexx(:,:), mexx0(:,:)
real(DP) :: exxe, exxe0
real(DP), PARAMETER :: Zero=0.0d0, One=1.0d0, Two=2.0d0, Pt5=0.50d0
TYPE(bec_type), INTENT(in) :: becpsi

  CALL start_clock( 'aceinit' )

  IF(nbndproj>nbnd) CALL errore('aceinit_k','nbndproj greater than nbnd.',1)
  IF(nbndproj<=0) CALL errore('aceinit_k','nbndproj le 0.',1)

  ALLOCATE( mexx(nbndproj,nbndproj), mexx0(nbndproj,nbndproj) )
  xitmp = (Zero,Zero)
  mexx  = (Zero,Zero)
  mexx0 = (Zero,Zero)
! |xi> = Vx[phi]|phi>
  CALL vexx(npwx, nnpw, nbndproj, phi, xitmp, becpsi)
! mexx = <phi|Vx[phi]|phi>
  CALL matcalc_k('exact',.true.,0,current_k,npwx*npol,nbndproj,nbndproj,phi,xitmp,mexx,exxe)
#if defined(__DEBUG)
  WRITE(stdout,'(3(A,I3),A,I9,A,f12.6)') 'aceinit_k: nbnd=', nbnd, ' nbndproj=',nbndproj, &
                                    ' k=',current_k,' npw=',nnpw,' Ex(k)=',exxe
#endif
! |xi> = -One * Vx[phi]|phi> * rmexx^T
  CALL aceupdate_k(nbndproj,nnpw,xitmp,mexx)

  DEALLOCATE( mexx )

  CALL stop_clock( 'aceinit' )

END SUBROUTINE aceinit_k
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE aceupdate_k(nbndproj,nnpw,xitmp,mexx)
USE wvfct ,               ONLY : npwx
USE noncollin_module,     ONLY : noncolin, npol
IMPLICIT NONE
  INTEGER :: nbndproj,nnpw
  COMPLEX(DP) :: mexx(nbndproj,nbndproj), xitmp(npwx*npol,nbndproj)

  CALL start_clock('aceupdate')

! mexx = -(Cholesky(mexx))^-1
  mexx = -mexx
  CALL invchol_k( nbndproj, mexx )

! |xi> = -One * Vx[phi]|phi> * mexx^T
  CALL ZTRMM('R','L','C','N',npwx*npol,nbndproj,(1.0_dp,0.0_dp),mexx,nbndproj,&
          xitmp,npwx*npol)

  CALL stop_clock('aceupdate')

END SUBROUTINE aceupdate_k
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE vexxace_k(nnpw,nbnd,phi,exxe,vphi)
USE becmod,               ONLY : calbec
USE wvfct,                ONLY : current_k, npwx
USE noncollin_module,     ONLY : npol
!
! do the ACE potential and
! (optional) print the ACE matrix representation
!
IMPLICIT NONE
  real(DP) :: exxe
  INTEGER :: nnpw,nbnd,i
  COMPLEX(DP) :: phi(npwx*npol,nbnd)
  COMPLEX(DP),OPTIONAL :: vphi(npwx*npol,nbnd)
  COMPLEX(DP),ALLOCATABLE :: cmexx(:,:), vv(:,:)
  real*8, PARAMETER :: Zero=0.0d0, One=1.0d0, Two=2.0d0, Pt5=0.50d0

  CALL start_clock('vexxace')

  ALLOCATE( vv(npwx*npol,nbnd) )
  IF(present(vphi)) THEN
    vv = vphi
  ELSE
    vv = (Zero, Zero)
  ENDIF

! do the ACE potential
  ALLOCATE( cmexx(nbndproj,nbnd) )
  cmexx = (Zero,Zero)
! <xi|phi>
  CALL matcalc_k('<xi|phi>',.false.,0,current_k,npwx*npol,nbndproj,nbnd,xi(1,1,current_k),phi,cmexx,exxe)

! |vv> = |vphi> + (-One) * |xi> * <xi|phi>
  CALL ZGEMM ('N','N',npwx*npol,nbnd,nbndproj,-(One,Zero),xi(1,1,current_k),npwx*npol,cmexx,nbndproj,(One,Zero),vv,npwx*npol)

  IF(domat) THEN
     IF ( nbndproj /= nbnd) THEN
        DEALLOCATE( cmexx )
        ALLOCATE( cmexx(nbnd,nbnd) )
     END IF
     CALL matcalc_k('ACE',.true.,0,current_k,npwx*npol,nbnd,nbnd,phi,vv,cmexx,exxe)
#if defined(__DEBUG)
    WRITE(stdout,'(3(A,I3),A,I9,A,f12.6)') 'vexxace_k: nbnd=', nbnd, ' nbndproj=',nbndproj, &
                   ' k=',current_k,' npw=',nnpw, ' Ex(k)=',exxe
  ELSE
    WRITE(stdout,'(3(A,I3),A,I9)') 'vexxace_k: nbnd=', nbnd, ' nbndproj=',nbndproj, &
                   ' k=',current_k,' npw=',nnpw
#endif
  ENDIF

  IF(present(vphi)) vphi = vv
  DEALLOCATE( vv,cmexx )

  CALL stop_clock('vexxace')

END SUBROUTINE vexxace_k
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE vexx_loc(npw, nbnd, hpsi, mexx)
USE noncollin_module,  ONLY : npol
USE cell_base,         ONLY : omega, alat
USE wvfct,             ONLY : current_k
USE klist,             ONLY : xk, nks, nkstot
USE fft_interfaces,    ONLY : fwfft, invfft
USE mp,                ONLY : mp_stop, mp_barrier, mp_sum
USE mp_bands,          ONLY : intra_bgrp_comm, me_bgrp, nproc_bgrp
USE exx_base,          ONLY : nqs, xkq_collect, index_xkq, index_xk, &
     g2_convolution
implicit none
!
!  Exact exchange with SCDM orbitals. 
!    Vx|phi> =  Vx|psi> <psi|Vx|psi>^(-1) <psi|Vx|phi>
!  locmat contains localization integrals
!  mexx contains in output the exchange matrix
!   
  integer :: npw, nbnd, nrxxs, npairs, ntot, NBands 
  integer :: ig, ir, ik, ikq, iq, ibnd, jbnd, kbnd, NQR
  INTEGER :: current_ik
  real(DP) :: ovpairs(2), mexx(nbnd,nbnd), exxe
  complex(DP) :: hpsi(npw,nbnd)
  COMPLEX(DP),ALLOCATABLE :: rhoc(:), vc(:), RESULT(:,:) 
  REAL(DP),ALLOCATABLE :: fac(:)
  REAL(DP) :: xkp(3), xkq(3)
  INTEGER, EXTERNAL  :: global_kpoint_index


  write(stdout,'(5X,A)') ' ' 
  write(stdout,'(5X,A)') 'Exact-exchange with localized orbitals'

  CALL start_clock('vexxloc')

  write(stdout,'(7X,A,f24.12)') 'local_thr =', local_thr
  nrxxs= dfftt%nnr
  NQR = nrxxs*npol

! exchange projected onto localized orbitals 
  ALLOCATE( fac(dfftt%ngm) )
  ALLOCATE( rhoc(nrxxs), vc(NQR) )
  ALLOCATE( RESULT(nrxxs,nbnd) ) 

  current_ik = global_kpoint_index ( nkstot, current_k )
  xkp = xk(:,current_k)

  vc=(0.0d0, 0.0d0)
  npairs = 0 
  ovpairs = 0.0d0

  DO iq=1,nqs
     ikq  = index_xkq(current_ik,iq)
     ik   = index_xk(ikq)
     xkq  = xkq_collect(:,ikq)
     CALL g2_convolution(dfftt%ngm, gt, xkp, xkq, fac)
     RESULT = (0.0d0, 0.0d0)
     DO ibnd = 1, nbnd
       IF(x_occupation(ibnd,ikq).gt.0.0d0) THEN 
         DO ir = 1, NQR 
           rhoc(ir) = locbuff(ir,ibnd,ikq) * locbuff(ir,ibnd,ikq) / omega
         ENDDO
         CALL fwfft ('Rho', rhoc, dfftt)
         vc=(0.0d0, 0.0d0)
         DO ig = 1, dfftt%ngm
             vc(dfftt%nl(ig))  = fac(ig) * rhoc(dfftt%nl(ig)) 
             vc(dfftt%nlm(ig)) = fac(ig) * rhoc(dfftt%nlm(ig))
         ENDDO
         CALL invfft ('Rho', vc, dfftt)
         DO ir = 1, NQR 
           RESULT(ir,ibnd) = RESULT(ir,ibnd) + locbuff(ir,ibnd,ikq) * vc(ir) 
         ENDDO
       END IF 

       DO kbnd = 1, ibnd-1
         IF((locmat(ibnd,kbnd,ikq).gt.local_thr).and. &
            ((x_occupation(ibnd,ikq).gt.0.0d0).or.(x_occupation(kbnd,ikq).gt.0.0d0))) then 
           IF((x_occupation(ibnd,ikq).gt.0.0d0).or. &
              (x_occupation(kbnd,ikq).gt.0.0d0) ) &
              ovpairs(2) = ovpairs(2) + locmat(ibnd,kbnd,ikq) 
!          write(stdout,'(3I4,3f12.6,A)') ikq, ibnd, kbnd, x_occupation(ibnd,ikq), x_occupation(kbnd,ikq), locmat(ibnd,kbnd,ikq), ' IN '
           ovpairs(1) = ovpairs(1) + locmat(ibnd,kbnd,ikq) 
           DO ir = 1, NQR 
             rhoc(ir) = locbuff(ir,ibnd,ikq) * locbuff(ir,kbnd,ikq) / omega
           ENDDO
           npairs = npairs + 1
           CALL fwfft ('Rho', rhoc, dfftt)
           vc=(0.0d0, 0.0d0)
           DO ig = 1, dfftt%ngm
               vc(dfftt%nl(ig))  = fac(ig) * rhoc(dfftt%nl(ig)) 
               vc(dfftt%nlm(ig)) = fac(ig) * rhoc(dfftt%nlm(ig))
           ENDDO
           CALL invfft ('Rho', vc, dfftt)
           DO ir = 1, NQR 
             RESULT(ir,kbnd) = RESULT(ir,kbnd) + x_occupation(ibnd,ikq) * locbuff(ir,ibnd,ikq) * vc(ir) 
           ENDDO
           DO ir = 1, NQR 
             RESULT(ir,ibnd) = RESULT(ir,ibnd) + x_occupation(kbnd,ikq) * locbuff(ir,kbnd,ikq) * vc(ir) 
           ENDDO
!        ELSE 
!          write(stdout,'(3I4,3f12.6,A)') ikq, ibnd, kbnd, x_occupation(ibnd,ikq), x_occupation(kbnd,ikq), locmat(ibnd,kbnd,ikq), '      OUT '
         END IF 
       ENDDO 
     ENDDO 

     DO jbnd = 1, nbnd
       CALL fwfft( 'Wave' , RESULT(:,jbnd), dfftt )
       DO ig = 1, npw
          hpsi(ig,jbnd) = hpsi(ig,jbnd) - exxalfa*RESULT(dfftt%nl(ig),jbnd) 
       ENDDO
     ENDDO

   ENDDO 
   DEALLOCATE( fac, vc )
   DEALLOCATE( RESULT )

!  Localized functions to G-space and exchange matrix onto localized functions
   allocate( RESULT(npw,nbnd) )
   RESULT = (0.0d0,0.0d0)
   DO jbnd = 1, nbnd
     rhoc(:) = dble(locbuff(:,jbnd,ikq)) + (0.0d0,1.0d0)*0.0d0
     CALL fwfft( 'Wave' , rhoc, dfftt )
     DO ig = 1, npw
       RESULT(ig,jbnd) = rhoc(dfftt%nl(ig))
     ENDDO
   ENDDO
   deallocate ( rhoc )
   CALL matcalc('M1-',.true.,0,npw,nbnd,nbnd,RESULT,hpsi,mexx,exxe)
   deallocate( RESULT )

   NBands = int(sum(x_occupation(:,ikq)))
   ntot = NBands * (NBands-1)/2 + NBands * (nbnd-NBands)
   write(stdout,'(7X,2(A,I12),A,f12.2)') '  Pairs(full): ',      ntot, &
           '   Pairs(included): ', npairs, &
           '   Pairs(%): ', dble(npairs)/dble(ntot)*100.0d0
   write(stdout,'(7X,3(A,f12.6))')       'OvPairs(full): ', ovpairs(2), &
           ' OvPairs(included): ', ovpairs(1), &
           ' OvPairs(%): ', ovpairs(1)/ovpairs(2)*100.0d0

  CALL stop_clock('vexxloc')

END SUBROUTINE vexx_loc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE compute_density(DoPrint,Shift,CenterPBC,SpreadPBC,Overlap,PsiI,PsiJ,NQR,ibnd,jbnd)
USE constants,            ONLY : pi,  bohr_radius_angs 
USE cell_base,            ONLY : alat, omega
USE mp,                   ONLY : mp_sum
USE mp_bands,             ONLY : intra_bgrp_comm
! 
! Manipulate density: get pair density, center, spread, absolute overlap 
!    DoPrint : ... whether to print or not the quantities
!    Shift   : ... .false. refer the centers to the cell -L/2 ... +L/2
!                  .true.  shift the centers to the cell 0 ... L (presumably the one given in input) 
!
IMPLICIT NONE
  INTEGER,  INTENT(IN) :: NQR, ibnd, jbnd
  REAL(DP), INTENT(IN) :: PsiI(NQR), PsiJ(NQR) 
  LOGICAL,  INTENT(IN) :: DoPrint
  LOGICAL,  INTENT(IN) :: Shift ! .false. Centers with respect to the minimum image cell convention
                                ! .true.  Centers shifted to the input cell
  REAL(DP),    INTENT(OUT) :: Overlap, CenterPBC(3), SpreadPBC(3)

  REAL(DP) :: vol, rbuff, TotSpread
  INTEGER :: ir, i, j, k , idx, j0, k0
  COMPLEX(DP) :: cbuff(3)
  REAL(DP), PARAMETER :: Zero=0.0d0, One=1.0d0, Two=2.0d0 

  vol = omega / dble(dfftt%nr1 * dfftt%nr2 * dfftt%nr3)

  CenterPBC = Zero 
  SpreadPBC = Zero 
  Overlap = Zero
  cbuff = (Zero, Zero) 
  rbuff = Zero  

  j0 = dfftt%my_i0r2p ; k0 = dfftt%my_i0r3p
  DO ir = 1, dfftt%nr1x*dfftt%my_nr2p*dfftt%my_nr3p
     !
     ! ... three dimensional indexes
     !
     idx = ir -1
     k   = idx / (dfftt%nr1x*dfftt%my_nr2p)
     idx = idx - (dfftt%nr1x*dfftt%my_nr2p)*k
     k   = k + k0
     IF ( k .GE. dfftt%nr3 ) CYCLE
     j   = idx / dfftt%nr1x
     idx = idx - dfftt%nr1x * j
     j   = j + j0
     IF ( j .GE. dfftt%nr2 ) CYCLE
     i   = idx
     IF ( i .GE. dfftt%nr1 ) CYCLE
     !
     rbuff = PsiI(ir) * PsiJ(ir) / omega
     Overlap = Overlap + abs(rbuff)*vol
     cbuff(1) = cbuff(1) + rbuff*exp((Zero,One)*Two*pi*DBLE(i)/DBLE(dfftt%nr1))*vol 
     cbuff(2) = cbuff(2) + rbuff*exp((Zero,One)*Two*pi*DBLE(j)/DBLE(dfftt%nr2))*vol
     cbuff(3) = cbuff(3) + rbuff*exp((Zero,One)*Two*pi*DBLE(k)/DBLE(dfftt%nr3))*vol
  ENDDO

  call mp_sum(cbuff,intra_bgrp_comm)
  call mp_sum(Overlap,intra_bgrp_comm)

  CenterPBC(1) =  alat/Two/pi*aimag(log(cbuff(1)))
  CenterPBC(2) =  alat/Two/pi*aimag(log(cbuff(2)))
  CenterPBC(3) =  alat/Two/pi*aimag(log(cbuff(3)))

  IF(Shift) then 
    if(CenterPBC(1).lt.Zero) CenterPBC(1) = CenterPBC(1) + alat
    if(CenterPBC(2).lt.Zero) CenterPBC(2) = CenterPBC(2) + alat
    if(CenterPBC(3).lt.Zero) CenterPBC(3) = CenterPBC(3) + alat
  END IF 

  rbuff = dble(cbuff(1))**2 + aimag(cbuff(1))**2 
  SpreadPBC(1) = -(alat/Two/pi)**2 * dlog(rbuff) 
  rbuff = dble(cbuff(2))**2 + aimag(cbuff(2))**2 
  SpreadPBC(2) = -(alat/Two/pi)**2 * dlog(rbuff) 
  rbuff = dble(cbuff(3))**2 + aimag(cbuff(3))**2 
  SpreadPBC(3) = -(alat/Two/pi)**2 * dlog(rbuff) 
  TotSpread = (SpreadPBC(1) + SpreadPBC(2) + SpreadPBC(3))*bohr_radius_angs**2  

  IF(DoPrint) then 
    write(stdout,'(A,2I4)')     'MOs:                  ', ibnd, jbnd
    write(stdout,'(A,10f12.6)') 'Absolute Overlap:     ', Overlap 
    write(stdout,'(A,10f12.6)') 'Center(PBC)[A]:       ',  CenterPBC(1)*bohr_radius_angs, &
            CenterPBC(2)*bohr_radius_angs, CenterPBC(3)*bohr_radius_angs
    write(stdout,'(A,10f12.6)') 'Spread [A**2]:        ', SpreadPBC(1)*bohr_radius_angs**2,&
            SpreadPBC(2)*bohr_radius_angs**2, SpreadPBC(3)*bohr_radius_angs**2
    write(stdout,'(A,10f12.6)') 'Total Spread [A**2]:  ', TotSpread
  END IF 

  IF(TotSpread.lt.Zero) Call errore('compute_density','Negative spread found',1)

END SUBROUTINE compute_density 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END MODULE exx
!-----------------------------------------------------------------------
