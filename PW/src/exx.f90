! Copyright (C) 2005-2015 Quantum ESPRESSO group
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
  USE coulomb_vcut_module,  ONLY : vcut_init, vcut_type, vcut_info, &
                                   vcut_get,  vcut_spheric_get
  USE noncollin_module,     ONLY : noncolin, npol
  USE io_global,            ONLY : ionode
  USE fft_custom,           ONLY : fft_cus
  !
  USE control_flags,        ONLY : gamma_only, tqr
  USE fft_types,            ONLY : fft_type_descriptor
  USE stick_base,           ONLY : sticks_map, sticks_map_deallocate
  !
  USE io_global,            ONLY : stdout

  IMPLICIT NONE
  SAVE
  COMPLEX(DP), ALLOCATABLE :: psi_exx(:,:), hpsi_exx(:,:)
  COMPLEX(DP), ALLOCATABLE :: evc_exx(:,:), psic_exx(:)
  INTEGER :: lda_original, n_original
  !
  ! general purpose vars
  !
  REAL(DP):: exxalfa=0._dp                ! 1 if exx, 0 elsewhere
  !
  ! variables defining the auxiliary k-point grid
  ! used in X BZ integration
  !
  INTEGER :: nq1=1, nq2=1, nq3=1         ! integers defining the X integration mesh
  INTEGER :: nqs=1                       ! number of points in the q-grid
  INTEGER :: nkqs                        ! total number of different k+q
  !
  REAL(DP),    ALLOCATABLE :: xkq_collect(:,:)  ! xkq(3,nkqs) the auxiliary k+q set
  REAL(DP),    ALLOCATABLE :: x_occupation(:,:)
                                         ! x_occupation(nbnd,nkstot) the weight of
                                         ! auxiliary functions in the density matrix

  INTEGER :: x_nbnd_occ                  ! number of bands of auxiliary functions with
                                         ! at least some x_occupation > eps_occ

  INTEGER :: ibnd_start = 0              ! starting band index used in bgrp parallelization
  INTEGER :: ibnd_end = 0                ! ending band index used in bgrp parallelization

  COMPLEX(DP), ALLOCATABLE :: exxbuff(:,:,:)
                                         ! temporary buffer for wfc storage
  !
  LOGICAL :: use_ace=.true.  !  true: Use Lin Lin's ACE method
                             !  false: do not use ACE, use old algorithm instead
  COMPLEX(DP), ALLOCATABLE :: xi(:,:,:)  ! ACE projectors
  INTEGER :: nbndproj
  LOGICAL :: domat
  !
  LOGICAL :: use_scdm=.false. ! if .true. enable Lin Lin's SCDM localization
                              ! currently implemented only within ACE formalism
  REAL(DP):: local_thr        ! threshold for Lin Lin's SCDM localized orbitals:
                              ! discard contribution to V_x if overlap between
                              ! localized orbitals is smaller than "local_thr"
  !
#if defined(__USE_INTEL_HBM_DIRECTIVES)
!DIR$ ATTRIBUTES FASTMEM :: exxbuff
#elif defined(__USE_CRAY_HBM_DIRECTIVES)
!DIR$ memory(bandwidth) exxbuff
#endif
  !
  !
  ! let xk(:,ik) + xq(:,iq) = xkq(:,ikq) = S(isym)*xk(ik') + G
  !
  !     index_xkq(ik,iq) = ikq
  !     index_xk(ikq)    = ik'
  !     index_sym(ikq)   = isym
  !
  INTEGER, ALLOCATABLE :: index_xkq(:,:) ! index_xkq(nks,nqs)
  INTEGER, ALLOCATABLE :: index_xk(:)    ! index_xk(nkqs)
  INTEGER, ALLOCATABLE :: index_sym(:)   ! index_sym(nkqs)
  INTEGER, ALLOCATABLE :: rir(:,:)       ! rotations to take k to q
  INTEGER, ALLOCATABLE :: working_pool(:)! idx of pool that has the k+q-point before rotation
  !
  ! Internal:
  LOGICAL :: exx_grid_initialized = .false.
  !
  ! variables to deal with Coulomb divergence
  ! and related issues
  !
  REAL(DP)         :: eps  = 1.d-6
  REAL(DP)         :: eps_qdiv = 1.d-8 ! |q| > eps_qdiv
  REAL(DP), PARAMETER :: eps_occ  = 1.d-8 ! skip band with occupation < eps_occ
  REAL(DP)         :: exxdiv = 0._dp
  CHARACTER(32)    :: exxdiv_treatment  = ' '
  !
  ! x_gamma_extrapolation
  LOGICAL           :: x_gamma_extrapolation =.true.
  LOGICAl           :: on_double_grid =.false.
  REAL(DP)          :: grid_factor = 1.d0 !8.d0/7.d0
  !
  ! Gygi-Baldereschi
  LOGICAL           :: use_regularization = .true.
  !
  ! yukawa method
  REAL(DP)          :: yukawa = 0._dp
  !
  ! erfc screening
  REAL(DP)          :: erfc_scrlen = 0._dp
  !
  ! erf screening
  REAL(DP)          :: erf_scrlen = 0._dp
  !
  ! gau-pbe screening
  REAL (DP)         :: gau_scrlen = 0.d0
  !
  ! cutoff techniques
  LOGICAL           :: use_coulomb_vcut_ws = .false.
  LOGICAL           :: use_coulomb_vcut_spheric = .false.
  REAL(DP)          :: ecutvcut
  TYPE(vcut_type)   :: vcut

  !
  ! energy related variables
  !
  REAL(DP) :: fock0 = 0.0_DP, & !   sum <phi|Vx(phi)|phi>
              fock1 = 0.0_DP, & !   sum <psi|vx(phi)|psi>
              fock2 = 0.0_DP, & !   sum <psi|vx(psi)|psi>
              dexx  = 0.0_DP    !   fock1  - 0.5*(fock2+fock0)
  !
  ! custom fft grids
  !
  TYPE(fft_cus) exx_fft         ! Custom grid for Vx*psi calculation
  REAL(DP)  :: ecutfock         ! energy cutoff for custom grid
  !
  ! mapping for the data structure conversion
  !
  TYPE comm_packet
     INTEGER :: size
     INTEGER, ALLOCATABLE :: indices(:)
     COMPLEX(DP), ALLOCATABLE :: msg(:,:,:)
  END TYPE comm_packet
  TYPE(comm_packet), ALLOCATABLE :: comm_recv(:,:), comm_send(:,:)
  TYPE(comm_packet), ALLOCATABLE :: comm_recv_reverse(:,:)
  TYPE(comm_packet), ALLOCATABLE :: comm_send_reverse(:,:,:)
  INTEGER, ALLOCATABLE :: lda_local(:,:)
  INTEGER, ALLOCATABLE :: lda_exx(:,:)
  INTEGER, ALLOCATABLE :: ngk_local(:), ngk_exx(:)
  INTEGER, ALLOCATABLE :: igk_exx(:,:)
  INTEGER :: npwx_local = 0
  INTEGER :: npwx_exx = 0
  INTEGER :: npw_local = 0
  INTEGER :: npw_exx = 0
  INTEGER :: n_local = 0
  INTEGER :: nwordwfc_exx
  LOGICAL :: first_data_structure_change = .TRUE.

  INTEGER :: ngm_loc, ngm_g_loc, gstart_loc
  INTEGER, ALLOCATABLE :: ig_l2g_loc(:)
  REAL(DP), ALLOCATABLE :: g_loc(:,:), gg_loc(:)
  INTEGER, ALLOCATABLE :: mill_loc(:,:), nl_loc(:)
  INTEGER :: ngms_loc, ngms_g_loc
  INTEGER, ALLOCATABLE :: nls_loc(:)
  INTEGER, ALLOCATABLE :: nlm_loc(:)
  INTEGER, ALLOCATABLE :: nlsm_loc(:)

  !gcutm
  INTEGER :: ngm_exx, ngm_g_exx, gstart_exx
  INTEGER, ALLOCATABLE :: ig_l2g_exx(:)
  REAL(DP), ALLOCATABLE :: g_exx(:,:), gg_exx(:)
  INTEGER, ALLOCATABLE :: mill_exx(:,:), nl_exx(:)
  INTEGER :: ngms_exx, ngms_g_exx
  INTEGER, ALLOCATABLE :: nls_exx(:)
  INTEGER, ALLOCATABLE :: nlm_exx(:)
  INTEGER, ALLOCATABLE :: nlsm_exx(:)
  !gcutms

  !the coulomb factor is reused between iterations
  REAL(DP), ALLOCATABLE :: coulomb_fac(:,:,:)
  !list of which coulomb factors have been calculated already
  LOGICAL, ALLOCATABLE :: coulomb_done(:,:)

  TYPE(fft_type_descriptor) :: dfftp_loc, dffts_loc
  TYPE(fft_type_descriptor) :: dfftp_exx, dffts_exx
  TYPE (sticks_map) :: smap_exx ! Stick map descriptor
  INTEGER :: ngw_loc, ngs_loc
  INTEGER :: ngw_exx, ngs_exx

 CONTAINS
#define _CX(A)  CMPLX(A,0._dp,kind=DP)
#define _CY(A)  CMPLX(0._dp,-A,kind=DP)
  !
  !------------------------------------------------------------------------
  SUBROUTINE exx_fft_create ()
    USE gvecw,        ONLY : ecutwfc
    USE gvect,        ONLY : ecutrho, ig_l2g
    USE klist,        ONLY : qnorm
    USE cell_base,    ONLY : at, bg, tpiba2
    USE fft_custom,   ONLY : set_custom_grid, ggent
    USE control_flags,ONLY : tqr
    USE realus,       ONLY : qpointlist, tabxx, tabp!, tabs

    IMPLICIT NONE

    IF( exx_fft%initialized) RETURN

    ! Initialise the custom grid that allows us to put the wavefunction
    ! onto the new (smaller) grid for rho (and vice versa)
    !
    exx_fft%ecutt=ecutwfc
    ! with k-points the following instructions guarantees that the sphere in
    ! G space contains k+G points - needed if ecutfock \simeq ecutwfc
    IF ( gamma_only ) THEN
       exx_fft%dual_t = ecutfock/ecutwfc
    ELSE
       exx_fft%dual_t = max(ecutfock,(sqrt(ecutwfc)+qnorm)**2)/ecutwfc
    ENDIF
    !
    exx_fft%gcutmt = exx_fft%dual_t*exx_fft%ecutt / tpiba2
    CALL data_structure_custom(exx_fft, smap_exx, gamma_only)
    CALL ggent(exx_fft)
    exx_fft%initialized = .true.
    !
    IF(tqr)THEN
      WRITE(stdout, '(5x,a)') "Initializing real-space augmentation for EXX grid"
      IF(ecutfock==ecutrho)THEN
        WRITE(stdout, '(7x,a)') " EXX grid -> DENSE grid"
        tabxx => tabp
!       ELSEIF(ecutfock==ecutwfc)THEN
!         WRITE(stdout, '(7x,a)') " EXX grid -> SMOOTH grid"
!         tabxx => tabs
      ELSE
        CALL qpointlist(exx_fft%dfftt, tabxx)
      ENDIF
    ENDIF

    RETURN
    !------------------------------------------------------------------------
  END SUBROUTINE exx_fft_create
  !------------------------------------------------------------------------
  !
  !------------------------------------------------------------------------
  SUBROUTINE deallocate_exx ()
    !------------------------------------------------------------------------
    !
    USE becmod, ONLY : deallocate_bec_type, is_allocated_bec_type, bec_type
    USE us_exx, ONLY : becxx
    USE fft_custom,  ONLY : deallocate_fft_custom
    !
    IMPLICIT NONE
    INTEGER :: ikq
    !
    IF ( allocated(index_xkq) ) DEALLOCATE(index_xkq)
    IF ( allocated(index_xk ) ) DEALLOCATE(index_xk )
    IF ( allocated(index_sym) ) DEALLOCATE(index_sym)
    IF ( allocated(rir)       ) DEALLOCATE(rir)
    IF ( allocated(x_occupation) ) DEALLOCATE(x_occupation)
    IF ( allocated(xkq_collect) )  DEALLOCATE(xkq_collect)
    IF ( allocated(exxbuff) )      DEALLOCATE(exxbuff)
!civn
    IF ( allocated(xi) )      DEALLOCATE(xi)
    !
    IF(allocated(becxx)) THEN
      DO ikq = 1, nkqs
        IF(is_allocated_bec_type(becxx(ikq))) CALL deallocate_bec_type(becxx(ikq))
      ENDDO
      DEALLOCATE(becxx)
    ENDIF
    !
    CALL deallocate_fft_custom(exx_fft)
    !
    !------------------------------------------------------------------------
  END SUBROUTINE deallocate_exx
  !------------------------------------------------------------------------
  !
  SUBROUTINE exx_grid_reinit()
    IMPLICIT NONE
    DEALLOCATE(xkq_collect,index_xk,index_sym)
    exx_grid_initialized = .false.
    nkqs = 0
    CALL exx_grid_init()
    !
    DEALLOCATE(working_pool)
    CALL exx_mp_init()
  END SUBROUTINE exx_grid_reinit
  
  !------------------------------------------------------------------------
  SUBROUTINE exx_mp_init()
    !------------------------------------------------------------------------
    ! 1. setup the orthopool communicators, which include the (n-1)th CPU of
    !    each pool (i.e. orthopool 0, contains CPU 0 of each pool)
    ! 2. setup global variable "working_pool" which contains the index of the
    !    pool which has the local copy prior to rotation of the ikq-th pool
    ! 3. working_pool(ikq) is also index of the CPU in each orthopool which has
    !    to broadcast the wavefunction at k-point ikq
    USE mp_images,      ONLY : intra_image_comm
    USE mp_pools,       ONLY : my_pool_id
    USE mp_orthopools,  ONLY : mp_start_orthopools, intra_orthopool_comm
    USE mp,             ONLY : mp_sum
    USE klist,          ONLY : nkstot
    !
    IMPLICIT NONE
    INTEGER :: ikq, current_ik, ik
    INTEGER, EXTERNAL :: local_kpoint_index
    !
    IF(nkqs<=0) CALL errore("exx_mp_init","exx_grid_init must be called first!",1)
    !
    ! mp_start_orthopools contains a check to avoid double initialisation
    CALL mp_start_orthopools ( intra_image_comm )
    !
    ALLOCATE(working_pool(nkqs))
    !
    DO ikq = 1, nkqs
      DO current_ik = 1, nkstot
        ! find the index of point k+q from the global list
        IF ( index_xk(ikq) == current_ik) EXIT
      ENDDO
      IF(current_ik>nkstot) CALL errore("exx_mp_init", "could not find the local index", ikq)
      ! Find the pool-index of the global point, returns -1 if the point
      ! was not in my pool
      ik = local_kpoint_index(nkstot, current_ik)
      !
      ! Set a temporary variable to the index of the pool who did the work
      IF( ik>0) THEN
        working_pool(ikq) = my_pool_id
      ELSE
        working_pool(ikq) = 0
      ENDIF
      !
    ENDDO    
    !
      ! Collect the variable among all the nth CPUs of each pool,
      ! i.e. inside each orto-pool group
    CALL mp_sum(working_pool, intra_orthopool_comm)
    !
  END SUBROUTINE exx_mp_init
  !
  !------------------------------------------------------------------------
  SUBROUTINE exx_grid_init()
    !------------------------------------------------------------------------
    !
    USE symm_base,  ONLY : nsym, s
    USE cell_base,  ONLY : bg, at
    USE spin_orb,   ONLY : domag
    USE noncollin_module, ONLY : nspin_lsda
    USE klist,      ONLY : xk, wk, nkstot, nks, qnorm
    USE wvfct,      ONLY : nbnd
    USE start_k,    ONLY : nk1,nk2,nk3
    USE control_flags, ONLY : iverbosity
    !
    IMPLICIT NONE
    !
    CHARACTER(13) :: sub_name='exx_grid_init'
    INTEGER       :: iq1, iq2, iq3, isym, ik, ikq, iq, max_nk, temp_nkqs
    INTEGER, ALLOCATABLE :: temp_index_xk(:), temp_index_sym(:)
    INTEGER, ALLOCATABLE :: temp_index_ikq(:), new_ikq(:)
    REAL(DP),ALLOCATABLE :: temp_xkq(:,:), xk_collect(:,:)
    LOGICAL      :: xk_not_found
    REAL(DP)     :: sxk(3), dxk(3), xk_cryst(3)
    REAL(DP)     :: dq1, dq2, dq3
    CHARACTER (len=6), EXTERNAL :: int_to_char
    !
    CALL start_clock ('exx_grid')
    !
    IF(nq1<=0) nq1 = nk1
    IF(nq2<=0) nq2 = nk2
    IF(nq3<=0) nq3 = nk3
    IF(nkstot==nspin_lsda) THEN
      nq1=1; nq2=1; nq3=1
    ENDIF

    IF(any((/nq1,nq2,nq3/)<=0)) CALL errore('exx_grid_init',"wrong EXX q grid", 1)
    !
    IF(exx_grid_initialized) CALL errore('exx_grid_init', "grid already initialized",1)
    exx_grid_initialized = .true.
    !
    ! definitions and checks
    !
    grid_factor = 1._dp
    IF (x_gamma_extrapolation) &
        grid_factor = 8.d0/7.d0
    !
    nqs = nq1 * nq2 * nq3
    !
    ! all processors on all pools need to have access to all k+q points
    !
    ALLOCATE(xk_collect(3,nkstot))
    CALL poolcollect(3, nks, xk, nkstot, xk_collect)
    !
    ! set a safe limit as the maximum number of auxiliary points we may need
    ! and allocate auxiliary arrays
    max_nk = nkstot * min(48, 2 * nsym)
    ALLOCATE( temp_index_xk(max_nk), temp_index_sym(max_nk) )
    ALLOCATE( temp_index_ikq(max_nk), new_ikq(max_nk) )
    ALLOCATE( temp_xkq(3,max_nk) )
    !
    ! find all k-points equivalent by symmetry to the points in the k-list
    !
    temp_nkqs = 0
    DO isym=1,nsym
      DO ik =1, nkstot
          ! go to crystalline coordinates
          xk_cryst(:) = xk_collect(:,ik)
          CALL cryst_to_cart(1, xk_cryst, at, -1)
          ! rotate with this sym.op.
          sxk(:) = s(:,1,isym)*xk_cryst(1) + &
                   s(:,2,isym)*xk_cryst(2) + &
                   s(:,3,isym)*xk_cryst(3)
          ! add sxk to the auxiliary list IF it is not already present
          xk_not_found = .true.
          ! *** do-loop skipped the first time because temp_nkqs == 0
          DO ikq=1, temp_nkqs
            IF (xk_not_found ) THEN
                dxk(:) = sxk(:)-temp_xkq(:,ikq) - nint(sxk(:)-temp_xkq(:,ikq))
                IF ( abs(dxk(1))<=eps .and. &
                     abs(dxk(2))<=eps .and. &
                     abs(dxk(3))<=eps ) xk_not_found = .false.
            ENDIF
          ENDDO
          IF (xk_not_found) THEN
            temp_nkqs                 = temp_nkqs + 1
            temp_xkq(:,temp_nkqs)     = sxk(:)
            temp_index_xk(temp_nkqs)  = ik
            temp_index_sym(temp_nkqs) = isym
          ENDIF

          sxk(:) = - sxk(:)
          xk_not_found = .true.
          DO ikq=1, temp_nkqs
            IF (xk_not_found ) THEN
                dxk(:) = sxk(:) - temp_xkq(:,ikq) - nint(sxk(:) - temp_xkq(:,ikq))
                IF ( abs(dxk(1))<=eps .and. &
                     abs(dxk(2))<=eps .and. &
                     abs(dxk(3))<=eps ) xk_not_found = .false.
            ENDIF
          ENDDO
          IF (xk_not_found .and. .not. (noncolin.and.domag) ) THEN
            temp_nkqs                 = temp_nkqs + 1
            temp_xkq(:,temp_nkqs)     = sxk(:)
            temp_index_xk(temp_nkqs)  = ik
            temp_index_sym(temp_nkqs) =-isym
          ENDIF

      ENDDO
    ENDDO

    !
    ! define the q-mesh step-sizes
    !
    dq1= 1._dp/dble(nq1)
    dq2= 1._dp/dble(nq2)
    dq3= 1._dp/dble(nq3)
    !
    ! allocate and fill the array index_xkq(nkstot,nqs)
    !
    IF(.not.allocated(index_xkq))    ALLOCATE( index_xkq(nkstot,nqs) )
    IF(.not.allocated(x_occupation)) ALLOCATE( x_occupation(nbnd,nkstot) )
    nkqs = 0
    new_ikq(:) = 0
    DO ik=1,nkstot
      ! go to crystalline coordinates
      xk_cryst(:) = xk_collect(:,ik)
      CALL cryst_to_cart(1, xk_cryst, at, -1)
      !
      iq = 0
      DO iq1=1, nq1
        sxk(1) = xk_cryst(1) + (iq1-1) * dq1
        DO iq2 =1, nq2
          sxk(2) = xk_cryst(2) + (iq2-1) * dq2
          DO iq3 =1, nq3
              sxk(3) = xk_cryst(3) + (iq3-1) * dq3
              iq = iq + 1

              xk_not_found = .true.
              !
              DO ikq=1, temp_nkqs
                IF ( xk_not_found ) THEN
                    dxk(:) = sxk(:)-temp_xkq(:,ikq) - nint(sxk(:)-temp_xkq(:,ikq))
                    IF ( all(abs(dxk) < eps ) ) THEN
                        xk_not_found = .false.
                        IF ( new_ikq(ikq) == 0) THEN
                            nkqs = nkqs + 1
                            temp_index_ikq(nkqs) = ikq
                            new_ikq(ikq) = nkqs
                        ENDIF
                        index_xkq(ik,iq) = new_ikq(ikq)
                    ENDIF
                ENDIF
              ENDDO ! ikq
              !
              IF (xk_not_found) THEN
                WRITE (*,*) ik, iq, temp_nkqs
                WRITE (*,*) sxk(:)
                CALL errore(sub_name, ' k + q is not an S*k ', (ik-1) * nqs + iq )
              ENDIF

          ENDDO
        ENDDO
      ENDDO

    ENDDO
    !
    ! allocate and fill the arrays xkq(3,nkqs), index_xk(nkqs) and index_sym(nkqs)
    ! NOTE: nkqs will be redefined as nspin_lsda*nkqs later 
    !
    ALLOCATE( xkq_collect(3,nspin_lsda*nkqs), index_xk(nspin_lsda*nkqs),  &
              index_sym(nspin_lsda*nkqs) )

    DO ik =1, nkqs
      ikq               = temp_index_ikq(ik)
      xkq_collect(:,ik) = temp_xkq(:,ikq)
      index_xk(ik)      = temp_index_xk(ikq)
      index_sym(ik)     = temp_index_sym(ikq)
    ENDDO
    CALL cryst_to_cart(nkqs, xkq_collect, bg, +1)

    IF( nkqs > 1) THEN
      WRITE(stdout, '(5x,3a)') "EXX: setup a grid of "//trim(int_to_char(nkqs))&
                           //" q-points centered on each k-point"
      IF ( nkqs < 100 .OR. iverbosity > 0 ) THEN
          WRITE( stdout, '(5x,a)' ) '(k+q)-points:'
          DO ik = 1, nkqs
            WRITE( stdout, '(3f12.7,5x,2i5)') (xkq_collect(ikq,ik), ikq=1,3), &
                 index_xk(ik), index_sym(ik)
          ENDDO
      ELSE
          WRITE( stdout, '(5x,a)' ) "(set verbosity='high' to see the list)"
      END IF
    ELSE
      WRITE(stdout, '("EXX: grid of k+q points same as grid of k-points")')
    ENDIF

    ! if nspin == 2, the kpoints are repeated in couples (spin up, spin down)
    IF (nspin_lsda == 2) THEN
      DO ik = 1, nkstot/2
          DO iq =1, nqs
            index_xkq(nkstot/2+ik,iq) = index_xkq(ik,iq) + nkqs
          ENDDO
      ENDDO
      DO ikq=1,nkqs
          xkq_collect(:,ikq + nkqs) = xkq_collect(:,ikq)
          index_xk(ikq + nkqs)  = index_xk(ikq) + nkstot/2
          index_sym(ikq + nkqs) = index_sym(ikq)
      ENDDO
      nkqs = 2 * nkqs
    ENDIF
    !
    ! clean up
    DEALLOCATE(temp_index_xk, temp_index_sym, temp_index_ikq, new_ikq, temp_xkq)
    !
    ! check that everything is what it should be
    CALL exx_grid_check ( xk_collect(:,:) )
    DEALLOCATE( xk_collect )
    !
    ! qnorm = max |k+q|, useful for reduced-cutoff calculations with k-points
    !
    qnorm = 0.0_dp
    DO iq = 1,nkqs
       DO ik = 1,nks
          qnorm = max(qnorm, sqrt( sum((xk(:,ik)-xkq_collect(:,iq))**2) ))
       ENDDO
    ENDDO
    !
    CALL stop_clock ('exx_grid')
    !
    RETURN
    !------------------------------------------------------------------------
  END SUBROUTINE exx_grid_init
  !------------------------------------------------------------------------
  !
  !------------------------------------------------------------------------
  SUBROUTINE exx_div_check()
    !------------------------------------------------------------------------
    !
    USE cell_base,  ONLY : at, alat
    USE funct,      ONLY : get_screening_parameter
    !
    IMPLICIT NONE
    !
    REAL(DP)     :: atws(3,3)
    CHARACTER(13) :: sub_name='exx_div_check'

    !
    ! EXX singularity treatment
    !
    SELECT CASE ( trim(exxdiv_treatment) )
    CASE ( "gygi-baldereschi", "gygi-bald", "g-b", "gb" )
      !
      use_regularization = .true.
      !
    CASE ( "vcut_ws" )
      !
      use_regularization = .true.
      use_coulomb_vcut_ws = .true.
      IF ( x_gamma_extrapolation ) &
            CALL errore(sub_name,'cannot USE x_gamm_extrap and vcut_ws', 1)
      !
    CASE ( "vcut_spherical" )
      !
      use_regularization = .true.
      use_coulomb_vcut_spheric = .true.
      IF ( x_gamma_extrapolation ) &
            CALL errore(sub_name,'cannot USE x_gamm_extrap and vcut_spherical', 1)
      !
    CASE ( "none" )
      use_regularization = .false.
      !
    CASE DEFAULT
      CALL errore(sub_name,'invalid exxdiv_treatment: '//trim(exxdiv_treatment), 1)
    END SELECT
    !
    ! Set variables for Coulomb vcut
    ! NOTE: some memory is allocated inside this routine (in the var vcut)
    !       and should be deallocated somewehre, at the end of the run
    !
    IF ( use_coulomb_vcut_ws .or. use_coulomb_vcut_spheric ) THEN
        !
        ! build the superperiodicity direct lattice
        !
        atws = alat * at
        !
        atws(:,1) = atws(:,1) * nq1
        atws(:,2) = atws(:,2) * nq2
        atws(:,3) = atws(:,3) * nq3
        !
        CALL vcut_init( vcut, atws, ecutvcut )
        !
        IF ( ionode ) CALL vcut_info( stdout, vcut )
        !
    ENDIF
    RETURN
  !------------------------------------------------------------------------
  END SUBROUTINE exx_div_check
  !------------------------------------------------------------------------


  !------------------------------------------------------------------------
  SUBROUTINE exx_grid_check ( xk_collect )
    !------------------------------------------------------------------------
    USE symm_base,  ONLY : s
    USE cell_base,  ONLY : at
    USE klist,      ONLY : nkstot, xk
    IMPLICIT NONE
    REAL(dp), INTENT(in) :: xk_collect(:,:) 
    !
    REAL(DP) :: sxk(3), dxk(3), xk_cryst(3), xkk_cryst(3)
    INTEGER :: iq1, iq2, iq3, isym, ik, ikk, ikq, iq
    REAL(DP) :: dq1, dq2, dq3
    dq1= 1._dp/dble(nq1)
    dq2= 1._dp/dble(nq2)
    dq3= 1._dp/dble(nq3)

    DO ik =1, nkstot
      xk_cryst(:) = xk_collect(:,ik)
      CALL cryst_to_cart(1, xk_cryst, at, -1)
      !
      iq = 0
      DO iq1=1, nq1
        sxk(1) = xk_cryst(1) + (iq1-1) * dq1
        DO iq2 =1, nq2
          sxk(2) = xk_cryst(2) + (iq2-1) * dq2
          DO iq3 =1, nq3
              sxk(3) = xk_cryst(3) + (iq3-1) * dq3
              iq = iq + 1

              ikq  = index_xkq(ik,iq)
              ikk  = index_xk(ikq)
              isym = index_sym(ikq)

              xkk_cryst(:) = at(1,:)*xk_collect(1,ikk) + &
                             at(2,:)*xk_collect(2,ikk) + &
                             at(3,:)*xk_collect(3,ikk)
              IF (isym < 0 ) xkk_cryst(:) = - xkk_cryst(:)
              isym = abs (isym)
              dxk(:) = s(:,1,isym)*xkk_cryst(1) + &
                       s(:,2,isym)*xkk_cryst(2) + &
                       s(:,3,isym)*xkk_cryst(3) - sxk(:)
              dxk(:) = dxk(:) - nint(dxk(:))
              IF ( .not. ( abs(dxk(1))<=eps .and. &
                           abs(dxk(2))<=eps .and. &
                           abs(dxk(3))<=eps )   ) THEN
                  WRITE(*,*) ik,iq
                  WRITE(*,*) ikq,ikk,isym
                  WRITE(*,*) dxk(:)
                  CALL errore('exx_grid_check', 'something wrong', 1 )
              ENDIF

          ENDDO
        ENDDO
      ENDDO
    ENDDO
    !
    RETURN

    !------------------------------------------------------------------------
  END SUBROUTINE exx_grid_check
  !------------------------------------------------------------------------
  !
  !------------------------------------------------------------------------
  SUBROUTINE exx_restart( set_ace )
     !------------------------------------------------------------------------
     !This SUBROUTINE is called when restarting an exx calculation
     USE funct,                ONLY : get_exx_fraction, start_exx, &
                                      exx_is_active, get_screening_parameter

     IMPLICIT NONE
     LOGICAL, INTENT(in) :: set_ace
     !
     use_ace = set_ace
     erfc_scrlen = get_screening_parameter()
     exxdiv = exx_divergence()
     exxalfa = get_exx_fraction()
     CALL start_exx()
     CALL weights()
     CALL exxinit()
     fock0 = exxenergy2()
     RETURN
     !------------------------------------------------------------------------
  END SUBROUTINE exx_restart
  !------------------------------------------------------------------------
  !
  !------------------------------------------------------------------------
  SUBROUTINE exxinit()
  !------------------------------------------------------------------------

    ! This SUBROUTINE is run before the first H_psi() of each iteration.
    ! It saves the wavefunctions for the right density matrix, in real space
    !
    USE wavefunctions_module, ONLY : evc, psic
    USE io_files,             ONLY : nwordwfc, iunwfc_exx
    USE buffers,              ONLY : get_buffer
    USE wvfct,                ONLY : nbnd, npwx, wg, current_k
    USE klist,                ONLY : ngk, nks, nkstot, xk, wk, igk_k
    USE symm_base,            ONLY : nsym, s, sr
    USE mp_pools,             ONLY : npool, nproc_pool, me_pool, inter_pool_comm
    USE mp_exx,               ONLY : me_egrp, set_egrp_indices, negrp, &
                                     init_index_over_band, my_egrp_id,  &
                                     inter_egrp_comm, intra_egrp_comm, &
                                     iexx_start, iexx_end
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

    IMPLICIT NONE
    INTEGER :: ik,ibnd, i, j, k, ir, isym, ikq, ig
    INTEGER :: h_ibnd
    INTEGER :: ibnd_loop_start, ibnd_buff_start, ibnd_buff_end
    INTEGER :: ipol, jpol
    REAL(dp), ALLOCATABLE   :: occ(:,:)
    COMPLEX(DP),ALLOCATABLE :: temppsic(:)
#if defined(__USE_INTEL_HBM_DIRECTIVES)
!DIR$ ATTRIBUTES FASTMEM :: temppsic
#elif defined(__USE_CRAY_HBM_DIRECTIVES)
!DIR$ memory(bandwidth) temppsic
#endif
    COMPLEX(DP),ALLOCATABLE :: temppsic_nc(:,:), psic_nc(:,:)
    INTEGER :: nxxs, nrxxs
#if defined(__MPI)
    COMPLEX(DP),ALLOCATABLE  :: temppsic_all(:),      psic_all(:)
    COMPLEX(DP), ALLOCATABLE :: temppsic_all_nc(:,:), psic_all_nc(:,:)
#endif
    COMPLEX(DP) :: d_spin(2,2,48)
    INTEGER :: npw, current_ik
    INTEGER, EXTERNAL :: global_kpoint_index
    INTEGER :: ibnd_start_new, ibnd_end_new
    INTEGER :: ibnd_exx, evc_offset
    !
    CALL start_clock ('exxinit')
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
    nxxs = exx_fft%dfftt%nr1x *exx_fft%dfftt%nr2x *exx_fft%dfftt%nr3x
    nrxxs= exx_fft%dfftt%nnr

    !allocate psic_exx
    IF(.not.allocated(psic_exx))THEN
       ALLOCATE(psic_exx(nrxxs))
       psic_exx = 0.0_DP
    END IF
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

    IF ( .not. gamma_only ) CALL exx_set_symm ( )

    ! set occupations of wavefunctions used in the calculation of exchange term

    ALLOCATE ( occ(nbnd,nks) )
    DO ik =1,nks
       IF(abs(wk(ik)) > eps_occ ) THEN
          occ(1:nbnd,ik) = wg (1:nbnd, ik) / wk(ik)
       ELSE
          occ(1:nbnd,ik) = 0._dp
       ENDIF
    ENDDO
    CALL poolcollect(nbnd, nks, occ, nkstot, x_occupation)
    DEALLOCATE ( occ )

    ! find an upper bound to the number of bands with non zero occupation.
    ! Useful to distribute bands among band groups

    x_nbnd_occ = 0
    DO ik =1,nkstot
       DO ibnd = max(1,x_nbnd_occ), nbnd
          IF (abs(x_occupation(ibnd,ik)) > eps_occ ) x_nbnd_occ = ibnd
       ENDDO
    ENDDO

    CALL set_egrp_indices(x_nbnd_occ,ibnd_start,ibnd_end)
    CALL init_index_over_band(inter_egrp_comm,nbnd,nbnd)

    !this will cause exxbuff to be calculated for every band
    ibnd_start_new = iexx_start
    ibnd_end_new = iexx_end

    IF ( gamma_only ) THEN
        ibnd_buff_start = ibnd_start_new/2
        IF(mod(iexx_start,2)==1) ibnd_buff_start = ibnd_buff_start +1
        !
        ibnd_buff_end = ibnd_end_new/2
        IF(mod(iexx_end,2)==1) ibnd_buff_end = ibnd_buff_end +1
    ELSE
        ibnd_buff_start = ibnd_start_new
        ibnd_buff_end   = ibnd_end_new
    ENDIF
    !
    IF (.not. allocated(exxbuff)) THEN
       IF (gamma_only) THEN
          ALLOCATE( exxbuff(nrxxs*npol, ibnd_buff_start:ibnd_buff_end, nks))
       ELSE
          ALLOCATE( exxbuff(nrxxs*npol, ibnd_buff_start:ibnd_buff_end, nkqs))
       END IF
    END IF

    !assign buffer
!$omp parallel do collapse(3) default(shared) firstprivate(npol,nrxxs,nkqs,ibnd_buff_start,ibnd_buff_end) private(ir,ibnd,ikq,ipol)
    DO ikq=1,SIZE(exxbuff,3)
       DO ibnd=ibnd_buff_start,ibnd_buff_end
          DO ir=1,nrxxs*npol
             exxbuff(ir,ibnd,ikq)=(0.0_DP,0.0_DP)
          ENDDO
       ENDDO
    ENDDO
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
          h_ibnd = iexx_start/2
          !
          IF(mod(iexx_start,2)==0) THEN
             h_ibnd=h_ibnd-1
             ibnd_loop_start=iexx_start-1
          ELSE
             ibnd_loop_start=iexx_start
          ENDIF

          evc_offset = 0
          DO ibnd = ibnd_loop_start, iexx_end, 2
             ibnd_exx = ibnd
             h_ibnd = h_ibnd + 1
             !
             psic_exx(:) = ( 0._dp, 0._dp )
             !
             IF ( ibnd < iexx_end ) THEN
                IF ( ibnd == ibnd_loop_start .and. MOD(iexx_start,2) == 0 ) THEN
                   DO ig=1,exx_fft%npwt
                      psic_exx(exx_fft%nlt(ig))  = ( 0._dp, 1._dp )*evc_exx(ig,1)
                      psic_exx(exx_fft%nltm(ig)) = ( 0._dp, 1._dp )*conjg(evc_exx(ig,1))
                   ENDDO
                   evc_offset = -1
                ELSE
                   DO ig=1,exx_fft%npwt
                      psic_exx(exx_fft%nlt(ig))  = evc_exx(ig,ibnd-ibnd_loop_start+evc_offset+1)  &
                           + ( 0._dp, 1._dp ) * evc_exx(ig,ibnd-ibnd_loop_start+evc_offset+2)
                      psic_exx(exx_fft%nltm(ig)) = conjg( evc_exx(ig,ibnd-ibnd_loop_start+evc_offset+1) ) &
                           + ( 0._dp, 1._dp ) * conjg( evc_exx(ig,ibnd-ibnd_loop_start+evc_offset+2) )
                   ENDDO
                END IF
             ELSE
                DO ig=1,exx_fft%npwt
                   psic_exx(exx_fft%nlt (ig)) = evc_exx(ig,ibnd-ibnd_loop_start+evc_offset+1)
                   psic_exx(exx_fft%nltm(ig)) = conjg( evc_exx(ig,ibnd-ibnd_loop_start+evc_offset+1) )
                ENDDO
             ENDIF

             CALL invfft ('CustomWave', psic_exx, exx_fft%dfftt)

             exxbuff(1:nrxxs,h_ibnd,ik)=psic_exx(1:nrxxs)
             
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
                   temppsic_nc(exx_fft%nlt(igk_exx(ig,ik)),1) = evc_exx(ig,ibnd-iexx_start+1)
                ENDDO
!$omp end parallel do
                CALL invfft ('CustomWave', temppsic_nc(:,1), exx_fft%dfftt)
!$omp parallel do default(shared) private(ig) firstprivate(npw,ik,ibnd_exx,npwx)
                DO ig=1,npw
                   temppsic_nc(exx_fft%nlt(igk_exx(ig,ik)),2) = evc_exx(ig+npwx,ibnd-iexx_start+1)
                ENDDO
!$omp end parallel do
                CALL invfft ('CustomWave', temppsic_nc(:,2), exx_fft%dfftt)
             ELSE
!$omp parallel do default(shared) private(ir) firstprivate(nrxxs)
                DO ir=1,nrxxs
                   temppsic(ir) = ( 0._dp, 0._dp )
                ENDDO
!$omp parallel do default(shared) private(ig) firstprivate(npw,ik,ibnd_exx)
                DO ig=1,npw
                   temppsic(exx_fft%nlt(igk_exx(ig,ik))) = evc_exx(ig,ibnd-iexx_start+1)
                ENDDO
!$omp end parallel do
                CALL invfft ('CustomWave', temppsic, exx_fft%dfftt)
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
                      CALL gather_grid(exx_fft%dfftt, temppsic_nc(:,ipol), temppsic_all_nc(:,ipol))
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
                      CALL scatter_grid(exx_fft%dfftt,psic_all_nc(:,ipol), psic_nc(:,ipol))
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
                   CALL gather_grid(exx_fft%dfftt,temppsic,temppsic_all)
                   IF ( me_egrp == 0 ) THEN
!$omp parallel do default(shared) private(ir) firstprivate(isym)
                      DO ir=1,nxxs
                         psic_all(ir) = temppsic_all(rir(ir,isym))
                      ENDDO
!$omp end parallel do
                   ENDIF
                   CALL scatter_grid(exx_fft%dfftt,psic_all,psic_exx)
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
    IF ( use_ace) CALL aceinit ( )
    !
    CALL stop_clock ('exxinit')
    !
    !-----------------------------------------------------------------------
  END SUBROUTINE exxinit
  !-----------------------------------------------------------------------
  !
  !-----------------------------------------------------------------------
  SUBROUTINE exx_set_symm ( )
    !-----------------------------------------------------------------------
    !
    ! Uses nkqs and index_sym from module exx, computes rir
    !
    USE fft_custom,           ONLY : fft_cus
    USE symm_base,            ONLY : nsym, s, sr, ft
    !
    IMPLICIT NONE
    !
    INTEGER :: nxxs, nr1,nr2,nr3, nr1x,nr2x,nr3x
    INTEGER :: ftau(3), ikq, isym, i,j,k, ri,rj,rk, ir
    LOGICAL :: ispresent(nsym)
    REAL (dp) :: ft_(3), eps2 = 1.0d-5
    !
    nr1 = exx_fft%dfftt%nr1
    nr2 = exx_fft%dfftt%nr2
    nr3 = exx_fft%dfftt%nr3
    nr1x= exx_fft%dfftt%nr1x
    nr2x= exx_fft%dfftt%nr2x
    nr3x= exx_fft%dfftt%nr3x
    nxxs = nr1x*nr2x*nr3x
    IF(.not. allocated(rir)) ALLOCATE(rir(nxxs,nsym))
    rir = 0
    ispresent(1:nsym) = .false.

    DO ikq =1,nkqs
       isym = abs(index_sym(ikq))
       IF (.not. ispresent(isym) ) THEN
          ispresent(isym) = .true.
          IF ( mod(s(2, 1, isym) * nr1, nr2) /= 0 .or. &
               mod(s(3, 1, isym) * nr1, nr3) /= 0 .or. &
               mod(s(1, 2, isym) * nr2, nr1) /= 0 .or. &
               mod(s(3, 2, isym) * nr2, nr3) /= 0 .or. &
               mod(s(1, 3, isym) * nr3, nr1) /= 0 .or. &
               mod(s(2, 3, isym) * nr3, nr2) /= 0 ) THEN
             CALL errore ('exx_set_symm',' EXX smooth grid is not compatible &
                          & with symmetry: change ecutfock',isym)
          ENDIF
          DO ir=1, nxxs
             rir(ir,isym) = ir
          ENDDO
          ! fractional translation in FFT grid coordinates
          ft_(1) = ft(1,isym)*nr1
          ft_(2) = ft(2,isym)*nr2
          ft_(3) = ft(3,isym)*nr3
          ftau(:) = NINT(ft_(:))
          IF ( abs (ft_(1) - ftau(1) ) / nr1 > eps2 .or. &
               abs (ft_(2) - ftau(2) ) / nr2 > eps2 .or. &
               abs (ft_(3) - ftau(3) ) / nr3 > eps2 ) THEN
             CALL infomsg ('exx_set_symm',' EXX smooth grid is not compatible &
                  & with fractional translation: change ecutfock')
          ENDIF
          DO k = 1, nr3
             DO j = 1, nr2
                DO i = 1, nr1
                   CALL ruotaijk (s(1,1,isym), ftau, i,j,k, nr1,nr2,nr3, ri,rj,rk)
                   ir =   i + ( j-1)*nr1x + ( k-1)*nr1x*nr2x
                   rir(ir,isym) = ri + (rj-1)*nr1x + (rk-1)*nr1x*nr2x
                ENDDO
             ENDDO
          ENDDO

       ENDIF
    ENDDO
  END SUBROUTINE exx_set_symm
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
    USE mp_exx,       ONLY : inter_egrp_comm, intra_egrp_comm, my_egrp_id, negrp, &
                             negrp, max_pairs, egrp_pairs, ibands, nibands, &
                             max_ibands, iexx_istart, iexx_iend, jblock
    USE mp,             ONLY : mp_sum, mp_barrier, mp_bcast
    USE uspp,           ONLY : nkb, okvan
    USE paw_variables,  ONLY : okpaw
    USE us_exx,         ONLY : bexg_merge, becxx, addusxx_g, addusxx_r, &
                               newdxx_g, newdxx_r, add_nlxx_pot, &
                               qvan_init, qvan_clean
    USE paw_exx,        ONLY : PAW_newdxx
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
    COMPLEX(DP),ALLOCATABLE :: result(:,:), result_g(:)
    REAL(DP),ALLOCATABLE :: temppsic_dble (:)
    REAL(DP),ALLOCATABLE :: temppsic_aimag(:)
    !
    COMPLEX(DP),ALLOCATABLE :: rhoc(:,:), vc(:,:), deexx(:,:)
    REAL(DP),   ALLOCATABLE :: fac(:)
    INTEGER          :: ibnd, ik, im , ikq, iq, ipol
    INTEGER          :: ir, ig
    INTEGER          :: current_ik
    INTEGER          :: ibnd_loop_start
    INTEGER          :: h_ibnd, nrxxs
    REAL(DP) :: x1, x2, xkp(3)
    REAL(DP) :: xkq(3)
    ! <LMS> temp array for vcut_spheric
    INTEGER, EXTERNAL :: global_kpoint_index
    LOGICAL :: l_fft_doubleband
    LOGICAL :: l_fft_singleband
    INTEGER :: ialloc
    COMPLEX(DP), ALLOCATABLE :: big_result(:,:)
    INTEGER :: iegrp, iproc, nproc_egrp, ii, ipair
    INTEGER :: jbnd, jstart, jend
    COMPLEX(DP), ALLOCATABLE :: exxtemp(:,:), psiwork(:)
    INTEGER :: ijt, njt, jblock_start, jblock_end
    INTEGER :: index_start, index_end, exxtemp_index
    INTEGER :: ending_im
    !
    ialloc = nibands(my_egrp_id+1)
    !
    ALLOCATE( fac(exx_fft%ngmt) )
    nrxxs= exx_fft%dfftt%nnr
    !
    !ALLOCATE( result(nrxxs), temppsic_dble(nrxxs), temppsic_aimag(nrxxs) )
    ALLOCATE( result(nrxxs,ialloc), temppsic_dble(nrxxs) )
    ALLOCATE( temppsic_aimag(nrxxs) )
    ALLOCATE( result_g(n) )
    ALLOCATE( psiwork(nrxxs) )
    !
    ALLOCATE( exxtemp(nrxxs*npol, jblock) )
    !
    ALLOCATE(rhoc(nrxxs,nbnd), vc(nrxxs,nbnd))
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
       CALL g2_convolution_all(exx_fft%ngmt, exx_fft%gt, xkp, xkq, iq, current_k)
       IF ( okvan .and..not.tqr ) CALL qvan_init (exx_fft%ngmt, xkq, xkp)
       !
       njt = nbnd / (2*jblock)
       if (mod(nbnd, (2*jblock)) .ne. 0) njt = njt + 1
       !
       IJT_LOOP : &
       DO ijt=1, njt
          !
          jblock_start = (ijt - 1) * (2*jblock) + 1
          jblock_end = min(jblock_start+(2*jblock)-1,nbnd)
          index_start = (ijt - 1) * jblock + 1
          index_end = min(index_start+jblock-1,nbnd/2)
          !
          !gather exxbuff for jblock_start:jblock_end
          call exxbuff_comm_gamma(exxtemp,current_k,nrxxs*npol,jblock_start,jblock_end,jblock)
          !
          LOOP_ON_PSI_BANDS : &
          DO ii = 1,  nibands(my_egrp_id+1)
             !
             ibnd = ibands(ii,my_egrp_id+1)
             !
             IF (ibnd.eq.0.or.ibnd.gt.m) CYCLE
             !
             psiwork = 0.0_DP
             !
             l_fft_doubleband = .false.
             l_fft_singleband = .false.
             !
             IF ( mod(ii,2)==1 .and. (ii+1)<=min(m,nibands(my_egrp_id+1)) ) l_fft_doubleband = .true.
             IF ( mod(ii,2)==1 .and. ii==min(m,nibands(my_egrp_id+1)) )     l_fft_singleband = .true.
             !
             IF( l_fft_doubleband ) THEN
!$omp parallel do  default(shared), private(ig)
                DO ig = 1, exx_fft%npwt
                   psiwork( exx_fft%nlt(ig) )  =       psi(ig, ii) + (0._DP,1._DP) * psi(ig, ii+1)
                   psiwork( exx_fft%nltm(ig) ) = conjg(psi(ig, ii) - (0._DP,1._DP) * psi(ig, ii+1))
                ENDDO
!$omp end parallel do
             ENDIF
             !
             IF( l_fft_singleband ) THEN
!$omp parallel do  default(shared), private(ig)
                DO ig = 1, exx_fft%npwt
                   psiwork( exx_fft%nlt(ig) )  =       psi(ig,ii) 
                   psiwork( exx_fft%nltm(ig) ) = conjg(psi(ig,ii))
                ENDDO
!$omp end parallel do
             ENDIF
             !
             IF( l_fft_doubleband.or.l_fft_singleband) THEN
                CALL invfft ('CustomWave', psiwork, exx_fft%dfftt)
!$omp parallel do default(shared), private(ir)
                DO ir = 1, nrxxs
                   temppsic_dble(ir)  = dble ( psiwork(ir) )
                   temppsic_aimag(ir) = aimag( psiwork(ir) )
                ENDDO
!$omp end parallel do
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
             h_ibnd = jstart/2
             IF(mod(jstart,2)==0) THEN
                h_ibnd=h_ibnd-1
                ibnd_loop_start=jstart-1
             ELSE
                ibnd_loop_start=jstart
             ENDIF
             !
             exxtemp_index = max(0, (jstart-jblock_start)/2 )
             !
             IBND_LOOP_GAM : &
             DO jbnd=ibnd_loop_start,jend, 2 !for each band of psi
                !
                h_ibnd = h_ibnd + 1
                !
                exxtemp_index = exxtemp_index + 1
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
                      rhoc(ir,ii) = exxtemp(ir,exxtemp_index) * temppsic_aimag(ir) / omega
                   ENDDO
!$omp end parallel do
                ELSE
!$omp parallel do default(shared), private(ir)
                   DO ir = 1, nrxxs
                      rhoc(ir,ii) = exxtemp(ir,exxtemp_index) * temppsic_dble(ir) / omega
                   ENDDO
!$omp end parallel do
                ENDIF
                !
                ! bring rho to G-space
                !
                !   >>>> add augmentation in REAL SPACE here
                IF(okvan .and. tqr) THEN
                   IF(jbnd>=jstart) &
                        CALL addusxx_r(exx_fft,rhoc(:,ii), &
                       _CX(becxx(ikq)%r(:,jbnd)), _CX(becpsi%r(:,ibnd)))
                   IF(jbnd<jend) &
                        CALL addusxx_r(exx_fft,rhoc(:,ii), &
                       _CY(becxx(ikq)%r(:,jbnd+1)),_CX(becpsi%r(:,ibnd)))
                ENDIF
                !
                CALL fwfft ('Custom', rhoc(:,ii), exx_fft%dfftt)
                !   >>>> add augmentation in G SPACE here
                IF(okvan .and. .not. tqr) THEN
                   ! contribution from one band added to real (in real space) part of rhoc
                   IF(jbnd>=jstart) &
                        CALL addusxx_g(exx_fft, rhoc(:,ii), xkq,  xkp, 'r', &
                        becphi_r=becxx(ikq)%r(:,jbnd), becpsi_r=becpsi%r(:,ibnd) )
                   ! contribution from following band added to imaginary (in real space) part of rhoc
                   IF(jbnd<jend) &
                        CALL addusxx_g(exx_fft, rhoc(:,ii), xkq,  xkp, 'i', &
                        becphi_r=becxx(ikq)%r(:,jbnd+1), becpsi_r=becpsi%r(:,ibnd) )
                ENDIF
                !   >>>> charge density done
                !
                vc(:,ii) = 0._DP
                !
!$omp parallel do default(shared), private(ig)
                DO ig = 1, exx_fft%ngmt
                   !
                   vc(exx_fft%nlt(ig),ii)  = coulomb_fac(ig,iq,current_k) * rhoc(exx_fft%nlt(ig),ii)
                   vc(exx_fft%nltm(ig),ii) = coulomb_fac(ig,iq,current_k) * rhoc(exx_fft%nltm(ig),ii)
                   !
                ENDDO
!$omp end parallel do
                !
                !   >>>>  compute <psi|H_fock G SPACE here
                IF(okvan .and. .not. tqr) THEN
                   IF(jbnd>=jstart) &
                        CALL newdxx_g(exx_fft, vc(:,ii), xkq, xkp, 'r', deexx(:,ii), &
                           becphi_r=x1*becxx(ikq)%r(:,jbnd))
                   IF(jbnd<jend) &
                        CALL newdxx_g(exx_fft, vc(:,ii), xkq, xkp, 'i', deexx(:,ii), &
                            becphi_r=x2*becxx(ikq)%r(:,jbnd+1))
                ENDIF
                !
                !brings back v in real space
                CALL invfft ('Custom', vc(:,ii), exx_fft%dfftt)
                !
                !   >>>>  compute <psi|H_fock REAL SPACE here
                IF(okvan .and. tqr) THEN
                   IF(jbnd>=jstart) &
                        CALL newdxx_r(exx_fft,vc(:,ii), _CX(x1*becxx(ikq)%r(:,jbnd)), deexx(:,ii))
                   IF(jbnd<jend) &
                        CALL newdxx_r(exx_fft,vc(:,ii), _CY(x2*becxx(ikq)%r(:,jbnd+1)), deexx(:,ii))
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
                   result(ir,ii) = result(ir,ii)+x1* dble(vc(ir,ii))* dble(exxtemp(ir,exxtemp_index))&
                                                +x2*aimag(vc(ir,ii))*aimag(exxtemp(ir,exxtemp_index))
                ENDDO
!$omp end parallel do
                !
             ENDDO &
             IBND_LOOP_GAM
             !
          ENDDO &
          LOOP_ON_PSI_BANDS
       ENDDO &
       IJT_LOOP
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
       CALL fwfft( 'CustomWave' , result(:,ii), exx_fft%dfftt )
       !communicate result
       DO ig = 1, n
          big_result(ig,ibnd) = big_result(ig,ibnd) - exxalfa*result(exx_fft%nlt(igk_exx(ig,current_k)),ii)
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
    !
    DEALLOCATE( result, temppsic_dble, temppsic_aimag)
    !
    DEALLOCATE(rhoc, vc, fac )
    !
    DEALLOCATE(exxtemp)
    !
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
    USE mp_exx,       ONLY : inter_egrp_comm, intra_egrp_comm, my_egrp_id, &
                             negrp, max_pairs, egrp_pairs, ibands, nibands, &
                             max_ibands, iexx_istart, iexx_iend, jblock
    USE mp,             ONLY : mp_sum, mp_barrier, mp_bcast
    USE uspp,           ONLY : nkb, okvan
    USE paw_variables,  ONLY : okpaw
    USE us_exx,         ONLY : bexg_merge, becxx, addusxx_g, addusxx_r, &
                               newdxx_g, newdxx_r, add_nlxx_pot, &
                               qvan_init, qvan_clean
    USE paw_exx,        ONLY : PAW_newdxx
    USE io_global,      ONLY : stdout
    !
    !
    IMPLICIT NONE
    !
    INTEGER                  :: lda, n, m
    COMPLEX(DP)              :: psi(lda*npol,max_ibands)
    COMPLEX(DP)              :: hpsi(lda*npol,max_ibands)
    TYPE(bec_type), OPTIONAL :: becpsi ! or call a calbec(...psi) instead
    !
    ! local variables
    COMPLEX(DP),ALLOCATABLE :: temppsic(:,:), result(:,:)
#if defined(__USE_INTEL_HBM_DIRECTIVES)
!DIR$ ATTRIBUTES FASTMEM :: result
#elif defined(__USE_CRAY_HBM_DIRECTIVES)
!DIR$ memory(bandwidth) result
#endif
    COMPLEX(DP),ALLOCATABLE :: temppsic_nc(:,:,:),result_nc(:,:,:)
    INTEGER          :: request_send, request_recv
    !
	COMPLEX(DP),ALLOCATABLE :: deexx(:,:)
    COMPLEX(DP),ALLOCATABLE,TARGET :: rhoc(:,:), vc(:,:)
#if defined(__USE_3D_FFT) & defined(__USE_MANY_FFT)
    COMPLEX(DP),POINTER :: prhoc(:), pvc(:)
#endif
#if defined(__USE_INTEL_HBM_DIRECTIVES)
!DIR$ ATTRIBUTES FASTMEM :: rhoc, vc
#elif defined(__USE_CRAY_HBM_DIRECTIVES)
!DIR$ memory(bandwidth) rhoc, vc
#endif
    REAL(DP),   ALLOCATABLE :: fac(:), facb(:)
    INTEGER          :: ibnd, ik, im , ikq, iq, ipol
    INTEGER          :: ir, ig, ir_start, ir_end
	INTEGER			 :: irt, nrt, nblock
    INTEGER          :: current_ik
    INTEGER          :: ibnd_loop_start
    INTEGER          :: h_ibnd, nrxxs
    REAL(DP) :: x1, x2, xkp(3), omega_inv, nqs_inv
    REAL(DP) :: xkq(3)
    ! <LMS> temp array for vcut_spheric
    INTEGER, EXTERNAL :: global_kpoint_index
    DOUBLE PRECISION :: max, tempx
    COMPLEX(DP), ALLOCATABLE :: big_result(:,:)
    INTEGER :: ir_out, ipair, jbnd
    INTEGER :: ii, jstart, jend, jcount, jind
    INTEGER :: ialloc, ending_im
    INTEGER :: ijt, njt, jblock_start, jblock_end
    COMPLEX(DP), ALLOCATABLE :: exxtemp(:,:)
    !
    ialloc = nibands(my_egrp_id+1)
    !
    ALLOCATE( fac(exx_fft%ngmt) )
    nrxxs= exx_fft%dfftt%nnr
    ALLOCATE( facb(nrxxs) )
    !
    IF (noncolin) THEN
       ALLOCATE( temppsic_nc(nrxxs,npol,ialloc), result_nc(nrxxs,npol,ialloc) )
    ELSE
       ALLOCATE( temppsic(nrxxs,ialloc), result(nrxxs,ialloc) )
    ENDIF
    !
    ALLOCATE( exxtemp(nrxxs*npol, jblock) )
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
    ALLOCATE(rhoc(nrxxs,jblock), vc(nrxxs,jblock))
#if defined(__USE_3D_FFT) & defined(__USE_MANY_FFT)
    prhoc(1:nrxxs*jblock) => rhoc(:,:)
    pvc(1:nrxxs*jblock) => vc(:,:)
#endif
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
!$omp parallel do  default(shared), private(ir) firstprivate(nrxxs)
          DO ir = 1, nrxxs
             temppsic(ir,ii) = 0._DP
          ENDDO
       END IF
       !
       IF (noncolin) THEN
          !
!$omp parallel do  default(shared), private(ig)
          DO ig = 1, n
             temppsic_nc(exx_fft%nlt(igk_exx(ig,current_k)),1,ii) = psi(ig,ii)
             temppsic_nc(exx_fft%nlt(igk_exx(ig,current_k)),2,ii) = psi(npwx+ig,ii)
          ENDDO
!$omp end parallel do
          !
          CALL invfft ('CustomWave', temppsic_nc(:,1,ii), exx_fft%dfftt)
          CALL invfft ('CustomWave', temppsic_nc(:,2,ii), exx_fft%dfftt)
          !
       ELSE
          !
!$omp parallel do  default(shared), private(ig)
          DO ig = 1, n
             temppsic( exx_fft%nlt(igk_exx(ig,current_k)), ii ) = psi(ig,ii)
          ENDDO
!$omp end parallel do
          !
          CALL invfft ('CustomWave', temppsic(:,ii), exx_fft%dfftt)
          !
       END IF
       !
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
    END DO
    !
    !
    !
    !
    !
    !
    !
    !precompute these guys
    omega_inv = 1.0 / omega
    nqs_inv = 1.0 / nqs
    !
    !
    !
    !------------------------------------------------------------------------!
    ! Beginning of main loop
    !------------------------------------------------------------------------!
    DO iq=1, nqs
       !
       ikq  = index_xkq(current_ik,iq)
       ik   = index_xk(ikq)
       xkq  = xkq_collect(:,ikq)
       !
       ! calculate the 1/|r-r'| (actually, k+q+g) factor and place it in fac
       CALL g2_convolution_all(exx_fft%ngmt, exx_fft%gt, xkp, xkq, iq, current_k)
       !
! JRD - below not threaded
       facb = 0D0
       DO ig = 1, exx_fft%ngmt
          facb(exx_fft%nlt(ig)) = coulomb_fac(ig,iq,current_k)
       ENDDO
       !
       IF ( okvan .and..not.tqr ) CALL qvan_init (exx_fft%ngmt, xkq, xkp)
       !
       njt = nbnd / jblock
       if (mod(nbnd, jblock) .ne. 0) njt = njt + 1
       !
       DO ijt=1, njt
          !
          jblock_start = (ijt - 1) * jblock + 1
          jblock_end = min(jblock_start+jblock-1,nbnd)
          !
          !gather exxbuff for jblock_start:jblock_end
          call exxbuff_comm(exxtemp,ikq,nrxxs*npol,jblock_start,jblock_end)
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
             !loads the phi from file
             !
             IF (noncolin) THEN
!$omp parallel do collapse(2) default(shared) firstprivate(jstart,jend) private(jbnd,ir)
                DO jbnd=jstart, jend
                   DO ir = 1, nrxxs
                      rhoc(ir,jbnd-jstart+1) = ( conjg(exxtemp(ir,jbnd-jblock_start+1))*temppsic_nc(ir,1,ii) +&
                      conjg(exxtemp(nrxxs+ir,jbnd-jblock_start+1))*temppsic_nc(ir,2,ii) )/omega
                   END DO
                END DO
!$omp end parallel do
             ELSE
                nblock=2048
                nrt = nrxxs / nblock
                if (mod(nrxxs, nblock) .ne. 0) nrt = nrt + 1
!$omp parallel do collapse(2) default(shared) firstprivate(jstart,jend,nblock,nrxxs,omega_inv) private(irt,ir,ir_start,ir_end,jbnd)
                DO irt = 1, nrt
                   DO jbnd=jstart, jend
                      ir_start = (irt - 1) * nblock + 1
                      ir_end = min(ir_start+nblock-1,nrxxs)
!DIR$ vector nontemporal (rhoc)
                      DO ir = ir_start, ir_end
                         rhoc(ir,jbnd-jstart+1)=conjg(exxtemp(ir,jbnd-jblock_start+1))*temppsic(ir,ii) * omega_inv
                      ENDDO
                   ENDDO
                ENDDO
!$omp end parallel do
             ENDIF
             !
             !   >>>> add augmentation in REAL space HERE
             IF(okvan .and. tqr) THEN ! augment the "charge" in real space
                DO jbnd=jstart, jend
                   CALL addusxx_r(exx_fft, rhoc(:,jbnd-jstart+1), becxx(ikq)%k(:,jbnd), becpsi%k(:,ibnd))
                ENDDO
             ENDIF
             !
             !   >>>> brings it to G-space
#if defined(__USE_3D_FFT) & defined(__USE_MANY_FFT)
             CALL fwfft ('Custom', prhoc, exx_fft%dfftt, howmany=jcount)
#else
             DO jbnd=jstart, jend
                CALL fwfft('Custom', rhoc(:,jbnd-jstart+1), exx_fft%dfftt)
             ENDDO
#endif
             !
             !   >>>> add augmentation in G space HERE
             IF(okvan .and. .not. tqr) THEN
                DO jbnd=jstart, jend
                   CALL addusxx_g(exx_fft, rhoc(:,jbnd-jstart+1), xkq, xkp, &
                   'c', becphi_c=becxx(ikq)%k(:,jbnd),becpsi_c=becpsi%k(:,ibnd))
                ENDDO
             ENDIF
             !   >>>> charge done
             !
!call start_collection()
             nblock=2048
             nrt = nrxxs / nblock
             if (mod(nrxxs, nblock) .ne. 0) nrt = nrt + 1
!$omp parallel do collapse(2) default(shared) firstprivate(jstart,jend,nblock,nrxxs,nqs_inv) private(irt,ir,ir_start,ir_end,jbnd)
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
             !
             ! Add ultrasoft contribution (RECIPROCAL SPACE)
             ! compute alpha_I,j,k+q = \sum_J \int <beta_J|phi_j,k+q> V_i,j,k,q Q_I,J(r) d3r
             IF(okvan .and. .not. tqr) THEN
                DO jbnd=jstart, jend
                   CALL newdxx_g(exx_fft, vc(:,jbnd-jstart+1), xkq, xkp, 'c',&
                                 deexx(:,ii), becphi_c=becxx(ikq)%k(:,jbnd))
                ENDDO
             ENDIF
             !
             !brings back v in real space
#if defined(__USE_3D_FFT) & defined(__USE_MANY_FFT)
             !fft many
             CALL invfft ('Custom', pvc, exx_fft%dfftt, howmany=jcount)
#else
             DO jbnd=jstart, jend
                CALL invfft('Custom', vc(:,jbnd-jstart+1), exx_fft%dfftt)
             ENDDO
#endif
             !
             ! Add ultrasoft contribution (REAL SPACE)
             IF(okvan .and. tqr) THEN
                DO jbnd=jstart, jend
                   CALL newdxx_r(exx_fft, vc(:,jbnd-jstart+1), becxx(ikq)%k(:,jbnd),deexx(:,ii))
                ENDDO
             ENDIF
             !
             ! Add PAW one-center contribution
             IF(okpaw) THEN
                DO jbnd=jstart, jend
                   CALL PAW_newdxx(x_occupation(jbnd,ik)/nqs, becxx(ikq)%k(:,jbnd), becpsi%k(:,ibnd), deexx(:,ii))
                ENDDO
             ENDIF
             !
             !accumulates over bands and k points
             !
             IF (noncolin) THEN
!$omp parallel do default(shared) firstprivate(jstart,jend) private(ir,jbnd)
                DO ir = 1, nrxxs
                   DO jbnd=jstart, jend
                      result_nc(ir,1,ii)= result_nc(ir,1,ii) + vc(ir,jbnd-jstart+1) * exxtemp(ir,jbnd-jblock_start+1)
                      result_nc(ir,2,ii)= result_nc(ir,2,ii) + vc(ir,jbnd-jstart+1) * exxtemp(ir+nrxxs,jbnd-jblock_start+1)
                   ENDDO
                ENDDO
!$omp end parallel do
             ELSE
!call start_collection()
                nblock=2048
                nrt = nrxxs / nblock
                if (mod(nrxxs, nblock) .ne. 0) nrt = nrt + 1
!$omp parallel do default(shared) firstprivate(jstart,jend,nblock,nrxxs) private(ir,ir_start,ir_end,jbnd)
                DO irt = 1, nrt
                   DO jbnd=jstart, jend
                      ir_start = (irt - 1) * nblock + 1
                      ir_end = min(ir_start+nblock-1,nrxxs)
!!dir$ vector nontemporal (result)
                      DO ir = ir_start, ir_end
                         result(ir,ii) = result(ir,ii) + vc(ir,jbnd-jstart+1)*exxtemp(ir,jbnd-jblock_start+1)
                      ENDDO
                   ENDDO
                ENDDO
!$omp end parallel do
!call stop_collection()
             ENDIF
             !
             !----------------------------------------------------------------------!
             !INNER LOOP END
             !----------------------------------------------------------------------!
             !
          END DO !I-LOOP
       END DO !IJT
       !
       IF ( okvan .and..not.tqr ) CALL qvan_clean ()
    END DO
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
          CALL fwfft ('CustomWave', result_nc(:,1,ii), exx_fft%dfftt)
          CALL fwfft ('CustomWave', result_nc(:,2,ii), exx_fft%dfftt)
          DO ig = 1, n
             big_result(ig,ibnd) = big_result(ig,ibnd) - exxalfa*result_nc(exx_fft%nlt(igk_exx(ig,current_k)),1,ii)
             big_result(n+ig,ibnd) = big_result(n+ig,ibnd) - exxalfa*result_nc(exx_fft%nlt(igk_exx(ig,current_k)),2,ii)
          ENDDO
       ELSE
          !
          CALL fwfft ('CustomWave', result(:,ii), exx_fft%dfftt)
          DO ig = 1, n
             big_result(ig,ibnd) = big_result(ig,ibnd) - exxalfa*result(exx_fft%nlt(igk_exx(ig,current_k)),ii)
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
    IF (noncolin) THEN
       DEALLOCATE(temppsic_nc, result_nc)
    ELSE
       DEALLOCATE(temppsic, result)
    ENDIF
    DEALLOCATE(big_result)
    !
    DEALLOCATE(fac, facb )
    !
    DEALLOCATE(exxtemp)
    !
    IF(okvan) DEALLOCATE( deexx)
    !
    !------------------------------------------------------------------------
  END SUBROUTINE vexx_k
  !-----------------------------------------------------------------------
  !
  !-----------------------------------------------------------------------
  SUBROUTINE g2_convolution_all(ngm, g, xk, xkq, iq, current_k)
    !-----------------------------------------------------------------------
    ! Wrapper for g2_convolution
    USE kinds,     ONLY : DP
    USE klist,     ONLY : nks
    !
    IMPLICIT NONE
    !
    INTEGER,  INTENT(in)    :: ngm    ! Number of G vectors
    REAL(DP), INTENT(in)    :: g(3,ngm) ! Cartesian components of G vectors
    REAL(DP), INTENT(in)    :: xk(3)  ! current k vector
    REAL(DP), INTENT(in)    :: xkq(3) ! current q vector
    INTEGER, INTENT(in) :: current_k, iq
    !
    ! Check if coulomb_fac has been allocated
    !
    IF( .NOT.ALLOCATED( coulomb_fac ) ) ALLOCATE( coulomb_fac(ngm,nqs,nks) )
    !
    ! Check if coulomb_done has been allocated
    !
    IF( .NOT.ALLOCATED( coulomb_done) ) THEN
       ALLOCATE( coulomb_done(nqs,nks) )
       coulomb_done = .FALSE.
    END IF
    !
    ! return if this k and k' already computed, otherwise compute it
    !
    IF ( coulomb_done(iq,current_k) ) RETURN
    !
    CALL g2_convolution(ngm, g, xk, xkq, coulomb_fac(:,iq,current_k))
    coulomb_done(iq,current_k) = .TRUE.
    !
  END SUBROUTINE g2_convolution_all
  !
  !-----------------------------------------------------------------------
  SUBROUTINE g2_convolution(ngm, g, xk, xkq, fac)
  !-----------------------------------------------------------------------
    ! This routine calculates the 1/|r-r'| part of the exact exchange
    ! expression in reciprocal space (the G^-2 factor).
    ! It then regularizes it according to the specified recipe
    USE kinds,     ONLY : DP
    USE cell_base, ONLY : tpiba, at, tpiba2
    USE constants, ONLY : fpi, e2, pi
    !
    IMPLICIT NONE
    !
    INTEGER,  INTENT(in)    :: ngm   ! Number of G vectors
    REAL(DP), INTENT(in)    :: g(3,ngm) ! Cartesian components of G vectors
    REAL(DP), INTENT(in)    :: xk(3) ! current k vector
    REAL(DP), INTENT(in)    :: xkq(3) ! current q vector
    !
    REAL(DP), INTENT(inout) :: fac(ngm) ! Calculated convolution
    !
    !Local variables
    INTEGER :: ig !Counters
    REAL(DP) :: q(3), qq, x
    REAL(DP) :: grid_factor_track(ngm), qq_track(ngm)
    REAL(DP) :: nqhalf_dble(3)
    LOGICAL :: odg(3)
    !
    ! First the types of Coulomb potential that need q(3) and an external call
    !
    IF( use_coulomb_vcut_ws ) THEN
       DO ig = 1, ngm
          q(:)= ( xk(:) - xkq(:) + g(:,ig) ) * tpiba
          fac(ig) = vcut_get(vcut,q)
       ENDDO
       RETURN
    ENDIF
    !
    IF ( use_coulomb_vcut_spheric ) THEN
       DO ig = 1, ngm
          q(:)= ( xk(:) - xkq(:) + g(:,ig) ) * tpiba
          fac(ig) = vcut_spheric_get(vcut,q)
       ENDDO
       RETURN
    ENDIF
    !
    ! Now the Coulomb potential that are computed on the fly
    !
    nqhalf_dble(1:3) = (/ dble(nq1)*0.5_DP, dble(nq2)*0.5_DP, dble(nq3)*0.5_DP /)
    !
    ! Set the grid_factor_track and qq_track
    !
    IF( x_gamma_extrapolation ) THEN
!$omp parallel do default(shared), private(ig,q,x,odg)
       DO ig = 1, ngm
          q(:)= xk(:) - xkq(:) + g(:,ig)
          qq_track(ig) = sum(q(:)**2) * tpiba2
          x = (q(1)*at(1,1)+q(2)*at(2,1)+q(3)*at(3,1))*nqhalf_dble(1)
          odg(1) = abs(x-nint(x))<eps
          x = (q(1)*at(1,2)+q(2)*at(2,2)+q(3)*at(3,2))*nqhalf_dble(2)
          odg(2) = abs(x-nint(x))<eps
          x = (q(1)*at(1,3)+q(2)*at(2,3)+q(3)*at(3,3))*nqhalf_dble(3)
          odg(3) = abs(x-nint(x))<eps
          IF( all ( odg(:) ) ) THEN
             grid_factor_track(ig) = 0._DP ! on double grid
          ELSE
             grid_factor_track(ig) = grid_factor ! not on double grid
          ENDIF
       ENDDO
!$omp end parallel do
    ELSE
!$omp parallel do default(shared), private(ig,q)
       DO ig = 1, ngm
          q(:)= xk(:) - xkq(:) + g(:,ig)
          qq_track(ig) = sum(q(:)**2) * tpiba2
       ENDDO
!$omp end parallel do
       grid_factor_track = 1._DP
    ENDIF
    !
    ! The big loop
    !
!$omp parallel do default(shared), private(ig,qq)
    DO ig=1,ngm
      !
      qq = qq_track(ig)
      !
      IF(gau_scrlen > 0) THEN
         fac(ig)=e2*((pi/gau_scrlen)**(1.5_DP))*exp(-qq/4._DP/gau_scrlen) * grid_factor_track(ig)
         !
      ELSEIF (qq > eps_qdiv) THEN
         !
         IF ( erfc_scrlen > 0  ) THEN
            fac(ig)=e2*fpi/qq*(1._DP-exp(-qq/4._DP/erfc_scrlen**2)) * grid_factor_track(ig)
         ELSEIF( erf_scrlen > 0 ) THEN
            fac(ig)=e2*fpi/qq*(exp(-qq/4._DP/erf_scrlen**2)) * grid_factor_track(ig)
         ELSE
            fac(ig)=e2*fpi/( qq + yukawa ) * grid_factor_track(ig) ! as HARTREE
         ENDIF
         !
      ELSE
         !
         fac(ig)= - exxdiv ! or rather something ELSE (see F.Gygi)
         !
         IF ( yukawa > 0._DP.and. .not. x_gamma_extrapolation ) fac(ig) = fac(ig) + e2*fpi/( qq + yukawa )
         IF( erfc_scrlen > 0._DP.and. .not. x_gamma_extrapolation ) fac(ig) = fac(ig) + e2*pi/(erfc_scrlen**2)
         !
      ENDIF
      !
    ENDDO
!$omp end parallel do
  END SUBROUTINE g2_convolution
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
    USE wavefunctions_module,   ONLY : evc
    USE lsda_mod,               ONLY : lsda, current_spin, isk
    USE klist,                  ONLY : ngk, nks, xk
    USE mp_pools,               ONLY : inter_pool_comm
    USE mp_exx,                 ONLY : intra_egrp_comm, intra_egrp_comm, &
                                       negrp
    USE mp,                     ONLY : mp_sum
    USE becmod,                 ONLY : bec_type, allocate_bec_type, &
                                       deallocate_bec_type, calbec
    USE uspp,                   ONLY : okvan,nkb,vkb

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
    USE gvect,                   ONLY : ngm, gstart, g, nl
    USE wvfct,                   ONLY : nbnd, npwx, wg
    USE wavefunctions_module,    ONLY : evc
    USE klist,                   ONLY : xk, ngk, nks, nkstot
    USE lsda_mod,                ONLY : lsda, current_spin, isk
    USE mp_pools,                ONLY : inter_pool_comm
    USE mp_bands,                ONLY : intra_bgrp_comm
    USE mp_exx,                  ONLY : inter_egrp_comm, intra_egrp_comm, negrp, &
                                        init_index_over_band, max_pairs, egrp_pairs, &
                                        my_egrp_id, nibands, ibands, jblock
    USE mp,                      ONLY : mp_sum
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
    INTEGER  :: h_ibnd, nrxxs, current_ik, ibnd_loop_start
    REAL(DP) :: x1, x2
    REAL(DP) :: xkq(3), xkp(3), vc
    ! temp array for vcut_spheric
    INTEGER, EXTERNAL :: global_kpoint_index
    !
    TYPE(bec_type) :: becpsi
    COMPLEX(DP), ALLOCATABLE :: psi_t(:), prod_tot(:)
    REAL(DP),ALLOCATABLE :: temppsic_dble (:)
    REAL(DP),ALLOCATABLE :: temppsic_aimag(:)
    LOGICAL :: l_fft_doubleband
    LOGICAL :: l_fft_singleband
    INTEGER :: jmax, npw
    INTEGER :: istart, iend, ipair, ii, ialloc
    COMPLEX(DP), ALLOCATABLE :: exxtemp(:,:)
    INTEGER :: ijt, njt, jblock_start, jblock_end
    INTEGER :: index_start, index_end, exxtemp_index
    INTEGER :: calbec_start, calbec_end
    INTEGER :: intra_bgrp_comm_
    !
    CALL init_index_over_band(inter_egrp_comm,nbnd,nbnd)
    !
    CALL transform_evc_to_exx(0)
    !
    ialloc = nibands(my_egrp_id+1)
    !
    nrxxs= exx_fft%dfftt%nnr
    ALLOCATE( fac(exx_fft%ngmt) )
    !
    ALLOCATE(temppsic(nrxxs), temppsic_dble(nrxxs),temppsic_aimag(nrxxs))
    ALLOCATE( rhoc(nrxxs) )
    ALLOCATE( vkb_exx(npwx,nkb) )
    !
    ALLOCATE( exxtemp(nrxxs*npol, jblock) )
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
          CALL g2_convolution_all(exx_fft%ngmt, exx_fft%gt, xkp, xkq, iq, &
               current_ik)
          fac = coulomb_fac(:,iq,current_ik)
          fac(exx_fft%gstart_t:) = 2 * coulomb_fac(exx_fft%gstart_t:,iq,current_ik)
          IF ( okvan .and..not.tqr ) CALL qvan_init (exx_fft%ngmt, xkq, xkp)
          !
          jmax = nbnd
          DO jbnd = nbnd,1, -1
             IF ( abs(wg(jbnd,ikk)) < eps_occ) CYCLE
             jmax = jbnd
             exit
          ENDDO
          !
          njt = nbnd / (2*jblock)
          if (mod(nbnd, (2*jblock)) .ne. 0) njt = njt + 1
          !
          IJT_LOOP : &
          DO ijt=1, njt
             !
             jblock_start = (ijt - 1) * (2*jblock) + 1
             jblock_end = min(jblock_start+(2*jblock)-1,nbnd)
             index_start = (ijt - 1) * jblock + 1
             index_end = min(index_start+jblock-1,nbnd/2)
             !
             !gather exxbuff for jblock_start:jblock_end
             call exxbuff_comm_gamma(exxtemp,ikk,nrxxs*npol,jblock_start,jblock_end,jblock)
             !
             JBND_LOOP : &
             DO ii = 1, nibands(my_egrp_id+1)
                !
                jbnd = ibands(ii,my_egrp_id+1)
                !
                IF (jbnd.eq.0.or.jbnd.gt.nbnd) CYCLE
                !
                temppsic = 0._DP
                !
                l_fft_doubleband = .false.
                l_fft_singleband = .false.
                !
                IF ( mod(ii,2)==1 .and. (ii+1)<=nibands(my_egrp_id+1) ) l_fft_doubleband = .true.
                IF ( mod(ii,2)==1 .and. ii==nibands(my_egrp_id+1) )     l_fft_singleband = .true.
                !
                IF( l_fft_doubleband ) THEN
!$omp parallel do  default(shared), private(ig)
                   DO ig = 1, exx_fft%npwt
                      temppsic( exx_fft%nlt(ig) )  = &
                           evc_exx(ig,ii) + (0._DP,1._DP) * evc_exx(ig,ii+1)
                      temppsic( exx_fft%nltm(ig) ) = &
                           conjg(evc_exx(ig,ii) - (0._DP,1._DP) * evc_exx(ig,ii+1))
                   ENDDO
!$omp end parallel do
                ENDIF
                !
                IF( l_fft_singleband ) THEN
!$omp parallel do  default(shared), private(ig)
                   DO ig = 1, exx_fft%npwt
                      temppsic( exx_fft%nlt(ig) )  =       evc_exx(ig,ii)
                      temppsic( exx_fft%nltm(ig) ) = conjg(evc_exx(ig,ii))
                   ENDDO
!$omp end parallel do
                ENDIF
                !
                IF( l_fft_doubleband.or.l_fft_singleband) THEN
                   CALL invfft ('CustomWave', temppsic, exx_fft%dfftt)
!$omp parallel do default(shared), private(ir)
                   DO ir = 1, nrxxs
                      temppsic_dble(ir)  = dble ( temppsic(ir) )
                      temppsic_aimag(ir) = aimag( temppsic(ir) )
                   ENDDO
!$omp end parallel do
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
                h_ibnd = istart/2
                IF(mod(istart,2)==0) THEN
                   h_ibnd=h_ibnd-1
                   ibnd_loop_start=istart-1
                ELSE
                   ibnd_loop_start=istart
                ENDIF
                !
                exxtemp_index = max(0, (istart-jblock_start)/2 )
                !
                IBND_LOOP_GAM : &
                DO ibnd = ibnd_loop_start, iend, 2       !for each band of psi
                   !
                   h_ibnd = h_ibnd + 1
                   !
                   exxtemp_index = exxtemp_index + 1
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
                         rhoc(ir) = exxtemp(ir,exxtemp_index) * temppsic_aimag(ir) / omega
                      ENDDO
!$omp end parallel do
                   ELSE
!$omp parallel do default(shared), private(ir)
                      DO ir = 1, nrxxs
                         rhoc(ir) = exxtemp(ir,exxtemp_index) * temppsic_dble(ir) / omega
                      ENDDO
!$omp end parallel do
                   ENDIF
                   !
                   IF(okvan .and.tqr) THEN
                      IF(ibnd>=istart) &
                           CALL addusxx_r(exx_fft,rhoc,_CX(becxx(ikq)%r(:,ibnd)),&
                                          _CX(becpsi%r(:,jbnd)))
                      IF(ibnd<iend) &
                           CALL addusxx_r(exx_fft,rhoc,_CY(becxx(ikq)%r(:,ibnd+1)),&
                                          _CX(becpsi%r(:,jbnd)))
                   ENDIF
                   !
                   ! bring rhoc to G-space
                   CALL fwfft ('Custom', rhoc, exx_fft%dfftt)
                   !
                   IF(okvan .and..not.tqr) THEN
                      IF(ibnd>=istart ) &
                           CALL addusxx_g( exx_fft, rhoc, xkq, xkp, 'r', &
                           becphi_r=becxx(ikq)%r(:,ibnd), becpsi_r=becpsi%r(:,jbnd-calbec_start+1) )
                      IF(ibnd<iend) &
                           CALL addusxx_g( exx_fft, rhoc, xkq, xkp, 'i', &
                           becphi_r=becxx(ikq)%r(:,ibnd+1), becpsi_r=becpsi%r(:,jbnd-calbec_start+1) )
                   ENDIF
                   !
                   vc = 0.0_DP
!$omp parallel do  default(shared), private(ig),  reduction(+:vc)
                   DO ig = 1,exx_fft%ngmt
                      !
                      ! The real part of rhoc contains the contribution from band ibnd
                      ! The imaginary part    contains the contribution from band ibnd+1
                      !
                      vc = vc + fac(ig) * ( x1 * &
                           abs( rhoc(exx_fft%nlt(ig)) + conjg(rhoc(exx_fft%nltm(ig))) )**2 &
                                 +x2 * &
                           abs( rhoc(exx_fft%nlt(ig)) - conjg(rhoc(exx_fft%nltm(ig))) )**2 )
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
          ENDDO &
          IJT_LOOP
          IF ( okvan .and..not.tqr ) CALL qvan_clean ( )
          !
       ENDDO &
       IQ_LOOP
    ENDDO &
    IKK_LOOP
    !
    DEALLOCATE(temppsic,temppsic_dble,temppsic_aimag)
    !
    DEALLOCATE(exxtemp)
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
    USE gvect,                   ONLY : ngm, gstart, g, nl
    USE wvfct,                   ONLY : nbnd, npwx, wg
    USE wavefunctions_module,    ONLY : evc
    USE klist,                   ONLY : xk, ngk, nks, nkstot
    USE lsda_mod,                ONLY : lsda, current_spin, isk
    USE mp_pools,                ONLY : inter_pool_comm
    USE mp_exx,                  ONLY : inter_egrp_comm, intra_egrp_comm, negrp, &
                                        ibands, nibands, my_egrp_id, max_pairs, &
                                        egrp_pairs, init_index_over_band, &
                                        jblock
    USE mp_bands,                ONLY : intra_bgrp_comm
    USE mp,                      ONLY : mp_sum
    USE fft_interfaces,          ONLY : fwfft, invfft
!#if defined(__USE_3D_FFT) & defined(__USE_MANY_FFT)
!    USE fft_interfaces,          ONLY : fwfftm, invfftm
!#endif
    USE gvect,                   ONLY : ecutrho
    USE klist,                   ONLY : wk
    USE uspp,                    ONLY : okvan,nkb,vkb
    USE becmod,                  ONLY : bec_type, allocate_bec_type, &
                                        deallocate_bec_type, calbec
    USE paw_variables,           ONLY : okpaw
    USE paw_exx,                 ONLY : PAW_xx_energy
    USE us_exx,                  ONLY : bexg_merge, becxx, addusxx_g, &
                                        addusxx_r, qvan_init, qvan_clean
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
#if defined(__USE_3D_FFT) & defined(__USE_MANY_FFT)
    COMPLEX(DP), POINTER :: prhoc(:)
#endif
    REAL(DP),    ALLOCATABLE :: fac(:)
    INTEGER  :: npw, jbnd, ibnd, ibnd_inner_start, ibnd_inner_end, ibnd_inner_count, ik, ikk, ig, ikq, iq, ir
    INTEGER  :: h_ibnd, nrxxs, current_ik, ibnd_loop_start, nblock, nrt, irt, ir_start, ir_end
    REAL(DP) :: x1, x2
    REAL(DP) :: xkq(3), xkp(3), vc, omega_inv
	!REAL(DP), ALLOCATEABLE :: vcarr(:)
    ! temp array for vcut_spheric
    INTEGER, EXTERNAL :: global_kpoint_index
    !
    TYPE(bec_type) :: becpsi
    COMPLEX(DP), ALLOCATABLE :: psi_t(:), prod_tot(:)
    INTEGER :: intra_bgrp_comm_
    INTEGER :: ii, ialloc, jstart, jend, ipair
    INTEGER :: ijt, njt, jblock_start, jblock_end
    COMPLEX(DP), ALLOCATABLE :: exxtemp(:,:)
    !
    CALL init_index_over_band(inter_egrp_comm,nbnd,nbnd)
    !
    CALL transform_evc_to_exx(0)
    !
    ialloc = nibands(my_egrp_id+1)
    !
    nrxxs = exx_fft%dfftt%nnr
    ALLOCATE( fac(exx_fft%ngmt) )
    !
    IF (noncolin) THEN
       ALLOCATE(temppsic_nc(nrxxs,npol,ialloc))
    ELSE
       ALLOCATE(temppsic(nrxxs,ialloc))
    ENDIF
    !
    ALLOCATE( exxtemp(nrxxs*npol, jblock) )
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
                temppsic_nc(exx_fft%nlt(igk_exx(ig,ikk)),1,ii) = evc_exx(ig,ii)
                temppsic_nc(exx_fft%nlt(igk_exx(ig,ikk)),2,ii) = evc_exx(npwx+ig,ii)
             ENDDO
!$omp end parallel do
             !
             CALL invfft ('CustomWave', temppsic_nc(:,1,ii), exx_fft%dfftt)
             CALL invfft ('CustomWave', temppsic_nc(:,2,ii), exx_fft%dfftt)
             !
          ELSE
!$omp parallel do default(shared), private(ig)
             DO ig = 1, npw
                temppsic(exx_fft%nlt(igk_exx(ig,ikk)),ii) = evc_exx(ig,ii)
             ENDDO
!$omp end parallel do
             !
             CALL invfft ('CustomWave', temppsic(:,ii), exx_fft%dfftt)
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
          CALL g2_convolution_all(exx_fft%ngmt, exx_fft%gt, xkp, xkq, iq, ikk)
          IF ( okvan .and..not.tqr ) CALL qvan_init (exx_fft%ngmt, xkq, xkp)
          !
          njt = nbnd / jblock
          if (mod(nbnd, jblock) .ne. 0) njt = njt + 1
          !
          IJT_LOOP : &
          DO ijt=1, njt
             !
             jblock_start = (ijt - 1) * jblock + 1
             jblock_end = min(jblock_start+jblock-1,nbnd)
             !
             !gather exxbuff for jblock_start:jblock_end
             call exxbuff_comm(exxtemp,ikq,nrxxs*npol,jblock_start,jblock_end)
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
#if defined(__USE_3D_FFT) & defined(__USE_MANY_FFT)
                prhoc(1:nrxxs*ibnd_inner_count) => rhoc
#endif 
                !
                ! load the phi at this k+q and band
                IF (noncolin) THEN
!$omp parallel do collapse(2) default(shared) private(ir,ibnd) firstprivate(ibnd_inner_start,ibnd_inner_end)
                   DO ibnd = ibnd_inner_start, ibnd_inner_end
                      DO ir = 1, nrxxs
                         rhoc(ir,ibnd-ibnd_inner_start+1)=(conjg(exxtemp(ir,ibnd-jblock_start+1))*temppsic_nc(ir,1,ii) + &
                              conjg(exxtemp(ir+nrxxs,ibnd-jblock_start+1))*temppsic_nc(ir,2,ii) ) * omega_inv
                      ENDDO
                   ENDDO
!$omp end parallel do
                ELSE
					
                   !calculate rho in real space
                   nblock=2048
                   nrt = nrxxs / nblock
                   if (mod(nrxxs, nblock) .ne. 0) nrt = nrt + 1
!$omp parallel do collapse(2) default(shared) &
!$omp& firstprivate(ibnd_inner_start,ibnd_inner_end,nblock,nrxxs,omega_inv) &
!$omp& private(ir,irt,ir_start,ir_end,ibnd)
                   DO irt = 1, nrt
                      DO ibnd=ibnd_inner_start, ibnd_inner_end
                         ir_start = (irt - 1) * nblock + 1
                         ir_end = min(ir_start+nblock-1,nrxxs)
!DIR$ vector nontemporal (rhoc)
                         DO ir = ir_start, ir_end
                            rhoc(ir,ibnd-ibnd_inner_start+1) = omega_inv * &
                              conjg(exxtemp(ir,ibnd-jblock_start+1)) * &
                              temppsic(ir,ii)
                         ENDDO
                      ENDDO
                   ENDDO
!$omp end parallel do

                ENDIF
                
                ! augment the "charge" in real space
                IF(okvan .and. tqr) THEN
!$omp parallel do default(shared) private(ibnd) firstprivate(ibnd_inner_start,ibnd_inner_end)
                   DO ibnd = ibnd_inner_start, ibnd_inner_end
                      CALL addusxx_r(exx_fft, rhoc(:,ibnd-ibnd_inner_start+1), &
                                     becxx(ikq)%k(:,ibnd), becpsi%k(:,jbnd))
                   ENDDO
!$omp end parallel do
                ENDIF
                !
                ! bring rhoc to G-space
#if defined(__USE_3D_FFT) & defined(__USE_MANY_FFT)
                CALL fwfft ('Custom', prhoc, exx_fft%dfftt, howmany=ibnd_inner_count)
#else
                DO ibnd=ibnd_inner_start, ibnd_inner_end
                   CALL fwfft('Custom', rhoc(:,ibnd-ibnd_inner_start+1), exx_fft%dfftt)
                ENDDO
#endif
                
                ! augment the "charge" in G space
                IF(okvan .and. .not. tqr) THEN
                   DO ibnd = ibnd_inner_start, ibnd_inner_end
                      CALL addusxx_g(exx_fft, rhoc(:,ibnd-ibnd_inner_start+1), &
                           xkq, xkp, 'c', becphi_c=becxx(ikq)%k(:,ibnd), &
                           becpsi_c=becpsi%k(:,jbnd))
                   ENDDO
                ENDIF
                !
                
!$omp parallel do default(shared) reduction(+:energy) private(ig,ibnd,vc) firstprivate(ibnd_inner_start,ibnd_inner_end)
                DO ibnd = ibnd_inner_start, ibnd_inner_end
                   vc=0.0_DP
                   DO ig=1,exx_fft%ngmt
                      vc = vc + coulomb_fac(ig,iq,ikk) * &
                          dble(rhoc(exx_fft%nlt(ig),ibnd-ibnd_inner_start+1) *&
                          conjg(rhoc(exx_fft%nlt(ig),ibnd-ibnd_inner_start+1)))
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
    DEALLOCATE(exxtemp)
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
  FUNCTION exx_divergence ()
    !-----------------------------------------------------------------------
     USE constants,      ONLY : fpi, e2, pi
     USE cell_base,      ONLY : bg, at, alat, omega
     USE gvect,          ONLY : ngm, g
     USE gvecw,          ONLY : gcutw
     USE mp_exx,         ONLY : intra_egrp_comm
     USE mp,             ONLY : mp_sum

     IMPLICIT NONE
     REAL(DP) :: exx_divergence

     ! local variables
     INTEGER :: iq1,iq2,iq3, ig
     REAL(DP) :: div, dq1, dq2, dq3, xq(3), q_, qq, tpiba2, alpha, x, q(3)

     INTEGER :: nqq, iq
     REAL(DP) :: aa, dq

     CALL start_clock ('exx_div')

     tpiba2 = (fpi / 2.d0 / alat) **2

     alpha  = 10._dp / gcutw

     IF ( .not. use_regularization ) THEN
        exx_divergence = 0._dp
        RETURN
     ENDIF

     dq1= 1._dp/dble(nq1)
     dq2= 1._dp/dble(nq2) 
     dq3= 1._dp/dble(nq3) 

     div = 0._dp
     DO iq1=1,nq1
        DO iq2=1,nq2
           DO iq3=1,nq3
              xq(:) = bg(:,1) * (iq1-1) * dq1 + &
                      bg(:,2) * (iq2-1) * dq2 + &
                      bg(:,3) * (iq3-1) * dq3
              DO ig=1,ngm
                 q(1)= xq(1) + g(1,ig)
                 q(2)= xq(2) + g(2,ig)
                 q(3)= xq(3) + g(3,ig)
                 qq = ( q(1)**2 + q(2)**2 + q(3)**2 )
                 IF (x_gamma_extrapolation) THEN
                    on_double_grid = .true.
                    x= 0.5d0*(q(1)*at(1,1)+q(2)*at(2,1)+q(3)*at(3,1))*nq1
                    on_double_grid = on_double_grid .and. (abs(x-nint(x))<eps)
                    x= 0.5d0*(q(1)*at(1,2)+q(2)*at(2,2)+q(3)*at(3,2))*nq2
                    on_double_grid = on_double_grid .and. (abs(x-nint(x))<eps)
                    x= 0.5d0*(q(1)*at(1,3)+q(2)*at(2,3)+q(3)*at(3,3))*nq3
                    on_double_grid = on_double_grid .and. (abs(x-nint(x))<eps)
                 ENDIF
                 IF (.not.on_double_grid) THEN
                    IF ( qq > 1.d-8 ) THEN
                       IF ( erfc_scrlen > 0 ) THEN
                          div = div + exp( -alpha * qq) / qq * &
                                (1._dp-exp(-qq*tpiba2/4.d0/erfc_scrlen**2)) * grid_factor
                       ELSEIF ( erf_scrlen >0 ) THEN
                          div = div + exp( -alpha * qq) / qq * &
                                (exp(-qq*tpiba2/4.d0/erf_scrlen**2)) * grid_factor
                       ELSE

                          div = div + exp( -alpha * qq) / (qq + yukawa/tpiba2) &
                                                     * grid_factor
                       ENDIF
                    ENDIF
                 ENDIF
              ENDDO
           ENDDO
        ENDDO
     ENDDO
     CALL mp_sum(  div, intra_egrp_comm )
     IF (gamma_only) THEN
        div = 2.d0 * div
     ENDIF
     IF ( .not. x_gamma_extrapolation ) THEN
        IF ( yukawa > 0._dp) THEN
           div = div + tpiba2/yukawa
        ELSEIF( erfc_scrlen > 0._dp ) THEN
           div = div + tpiba2/4.d0/erfc_scrlen**2
        ELSE
           div = div - alpha
        ENDIF
     ENDIF

     div = div * e2 * fpi / tpiba2 / nqs

     alpha = alpha / tpiba2

     nqq = 100000
     dq = 5.0d0 / sqrt(alpha) /nqq
     aa = 0._dp
     DO iq=0,  nqq
        q_ = dq * (iq+0.5d0)
        qq = q_ * q_
        IF ( erfc_scrlen > 0 ) THEN
           aa = aa  -exp( -alpha * qq) * exp(-qq/4.d0/erfc_scrlen**2) * dq
        ELSEIF ( erf_scrlen > 0 ) THEN
           aa = 0._dp
        ELSE
           aa = aa - exp( -alpha * qq) * yukawa / (qq + yukawa) * dq
        ENDIF
     ENDDO
     aa = aa * 8.d0 /fpi
     aa = aa + 1._dp/sqrt(alpha*0.25d0*fpi)
     IF( erf_scrlen > 0) aa = 1._dp/sqrt((alpha+1._dp/4.d0/erf_scrlen**2)*0.25d0*fpi)
     div = div - e2*omega * aa

     exx_divergence = div * nqs
     CALL stop_clock ('exx_div')

     RETURN
    !-----------------------------------------------------------------------
  END FUNCTION exx_divergence
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
    USE wavefunctions_module, ONLY : evc
    USE klist,                ONLY : xk, ngk, nks
    USE lsda_mod,             ONLY : lsda, current_spin, isk
    USE gvect,                ONLY : g, nl
    USE mp_pools,             ONLY : npool, inter_pool_comm
    USE mp_exx,               ONLY : inter_egrp_comm, intra_egrp_comm, &
                                     ibands, nibands, my_egrp_id, jblock, &
                                     egrp_pairs, max_pairs, negrp
    USE mp,                   ONLY : mp_sum
    USE fft_base,             ONLY : dffts
    USE fft_interfaces,       ONLY : fwfft, invfft
    USE uspp,                 ONLY : okvan
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
    INTEGER  :: h_ibnd, nqi, iqi, beta, nrxxs, ngm
    INTEGER  :: ibnd_loop_start
    REAL(DP) :: x1, x2
    REAL(DP) :: qq, xk_cryst(3), sxk(3), xkq(3), vc(3,3), x, q(3)
    ! temp array for vcut_spheric
    REAL(DP) :: delta(3,3)
    COMPLEX(DP), ALLOCATABLE :: exxtemp(:,:)
    INTEGER :: jstart, jend, ipair, jblock_start, jblock_end
    INTEGER :: ijt, njt, ii, jcount, exxtemp_index

    CALL start_clock ('exx_stress')

    CALL transform_evc_to_exx(0)

    IF (npool>1) CALL errore('exx_stress','stress not available with pools',1)
    IF (noncolin) CALL errore('exx_stress','noncolinear stress not implemented',1)
    IF (okvan) CALL infomsg('exx_stress','USPP stress not tested')

    nrxxs = exx_fft%dfftt%nnr
    ngm   = exx_fft%ngmt
    delta = reshape( (/1._dp,0._dp,0._dp, 0._dp,1._dp,0._dp, 0._dp,0._dp,1._dp/), (/3,3/))
    exx_stress_ = 0._dp
    ALLOCATE( tempphic(nrxxs), temppsic(nrxxs), rhoc(nrxxs), fac(ngm) )
    ALLOCATE( fac_tens(3,3,ngm), fac_stress(ngm) )
    !
    ALLOCATE( exxtemp(nrxxs*npol, nbnd) )
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
        ! loop over bands
        njt = nbnd / jblock
        if (mod(nbnd, jblock) .ne. 0) njt = njt + 1
        !
        DO ijt=1, njt
           !
           !gather exxbuff for jblock_start:jblock_end
           IF (gamma_only) THEN
              !
              jblock_start = (ijt - 1) * (2*jblock) + 1
              jblock_end = min(jblock_start+(2*jblock)-1,nbnd)
              !
              call exxbuff_comm_gamma(exxtemp,ikq,nrxxs*npol,jblock_start,&
                   jblock_end,jblock)
           ELSE
              !
              jblock_start = (ijt - 1) * jblock + 1
              jblock_end = min(jblock_start+jblock-1,nbnd)
              !
              call exxbuff_comm(exxtemp,ikq,nrxxs*npol,jblock_start,jblock_end)
           END IF
           !
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
            !how many iters
            jcount=jend-jstart+1
            if(jcount<=0) cycle
            !
            temppsic(:) = ( 0._dp, 0._dp )
!$omp parallel do default(shared), private(ig)
            DO ig = 1, npw
                temppsic(exx_fft%nlt(igk_exx(ig,ikk))) = evc_exx(ig,ii)
            ENDDO
!$omp end parallel do
            !
            IF(gamma_only) THEN
!$omp parallel do default(shared), private(ig)
                DO ig = 1, npw
                    temppsic(exx_fft%nltm(igk_exx(ig,ikk))) = &
                         conjg(evc_exx(ig,ii))
                ENDDO
!$omp end parallel do
            ENDIF

            CALL invfft ('CustomWave', temppsic, exx_fft%dfftt)

                IF (gamma_only) THEN
                    !
                    h_ibnd = jstart/2
                    !
                    IF(mod(jstart,2)==0) THEN
                      h_ibnd=h_ibnd-1
                      ibnd_loop_start=jstart-1
                    ELSE
                      ibnd_loop_start=jstart
                    ENDIF
                    !
                    exxtemp_index = max(0, (jstart-jblock_start)/2 )
                    !
                    DO ibnd = ibnd_loop_start, jend, 2     !for each band of psi
                        !
                        h_ibnd = h_ibnd + 1
                        !
                        exxtemp_index = exxtemp_index + 1
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
                            tempphic(ir) = exxtemp(ir,exxtemp_index)
                            rhoc(ir)     = conjg(tempphic(ir))*temppsic(ir) / omega
                        ENDDO
!$omp end parallel do
                        ! bring it to G-space
                        CALL fwfft ('Custom', rhoc, exx_fft%dfftt)

                        vc = 0._dp
!$omp parallel do default(shared), private(ig), reduction(+:vc)
                        DO ig = 1, ngm
                            !
                            vc(:,:) = vc(:,:) + fac(ig) * x1 * &
                                      abs( rhoc(exx_fft%nlt(ig)) + &
                                      conjg(rhoc(exx_fft%nltm(ig))))**2 * &
                                      (fac_tens(:,:,ig)*fac_stress(ig)/2.d0 - delta(:,:)*fac(ig))
                            vc(:,:) = vc(:,:) + fac(ig) * x2 * &
                                      abs( rhoc(exx_fft%nlt(ig)) - &
                                      conjg(rhoc(exx_fft%nltm(ig))))**2 * &
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
                          tempphic(ir) = exxtemp(ir,ibnd-jblock_start+1)
                          rhoc(ir)     = conjg(tempphic(ir))*temppsic(ir) / omega
                      ENDDO
!$omp end parallel do

                      ! bring it to G-space
                      CALL fwfft ('Custom', rhoc, exx_fft%dfftt)

                      vc = 0._dp
!$omp parallel do default(shared), private(ig), reduction(+:vc)
                      DO ig = 1, ngm
                          vc(:,:) = vc(:,:) + rhoc(exx_fft%nlt(ig))  * &
                                        conjg(rhoc(exx_fft%nlt(ig)))* &
                                    (fac_tens(:,:,ig)*fac_stress(ig)/2.d0 - delta(:,:)*fac(ig))
                      ENDDO
!$omp end parallel do
                      vc = vc * x_occupation(ibnd,ik) / nqs / 4.d0
                      exx_stress_ = exx_stress_ + exxalfa * vc * wg(jbnd,ikk)

                    ENDDO

                ENDIF ! gamma or k-points

        ENDDO ! jbnd
        ENDDO !ijt
            ENDDO ! iqi
    ENDDO ! ikk

    DEALLOCATE(tempphic, temppsic, rhoc, fac, fac_tens, fac_stress )
    DEALLOCATE(exxtemp)
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
  !
  !-----------------------------------------------------------------------
  SUBROUTINE transform_evc_to_exx(type)
  !-----------------------------------------------------------------------
    USE mp_exx,               ONLY : negrp
    USE wvfct,                ONLY : npwx, nbnd
    USE io_files,             ONLY : nwordwfc, iunwfc, iunwfc_exx
    USE klist,                ONLY : nks, ngk, igk_k
    USE uspp,                 ONLY : nkb, vkb
    USE wavefunctions_module, ONLY : evc
    USE control_flags,        ONLY : io_level
    USE buffers,              ONLY : open_buffer, get_buffer, save_buffer
    USE mp_exx,               ONLY : max_ibands
    !
    !
    IMPLICIT NONE
    !
    INTEGER, intent(in) :: type
    INTEGER :: lda, n, ik
    LOGICAL :: exst_mem, exst_file
    !
    IF (negrp.eq.1) THEN
       !
       ! no change in data structure is necessary
       ! just copy all of the required data
       !
       ! get evc_exx
       !
       IF(.not.allocated(evc_exx))THEN
          ALLOCATE(evc_exx(npwx*npol,nbnd))
       END IF
       evc_exx = evc
       !
       ! get igk_exx
       !
       IF(.not.allocated(igk_exx)) THEN
          ALLOCATE( igk_exx( npwx, nks ) )
          igk_exx = igk_k
       END IF
       !
       ! get the wfc buffer is used
       !
       iunwfc_exx = iunwfc
       nwordwfc_exx = nwordwfc
       !
       RETURN
       !
    END IF
    !
    ! change the data structure of evc and igk
    !
    lda = npwx
    n = npwx 
    npwx_local = npwx
    IF( .not.allocated(ngk_local) ) allocate(ngk_local(nks))
    ngk_local = ngk
    !
    IF ( .not.allocated(comm_recv) ) THEN
       !
       ! initialize all of the conversion maps and change the data structure
       !
       CALL initialize_local_to_exact_map(lda, nbnd)
    ELSE
       !
       ! change the data structure
       !
       CALL change_data_structure(.TRUE.)
    END IF
    !
    lda = npwx
    n = npwx 
    npwx_exx = npwx
    IF( .not.allocated(ngk_exx) ) allocate(ngk_exx(nks))
    ngk_exx = ngk
    !
    ! get evc_exx
    !
    IF(.not.allocated(evc_exx))THEN
       ALLOCATE(evc_exx(lda*npol,max_ibands+2))
       !
       ! ... open files/buffer for wavefunctions (nwordwfc set in openfil)
       ! ... io_level > 1 : open file, otherwise: open buffer
       !
       nwordwfc_exx  = size(evc_exx)
       CALL open_buffer( iunwfc_exx, 'wfc_exx', nwordwfc_exx, io_level, &
            exst_mem, exst_file )
    END IF
    !
    DO ik=1, nks
       !
       ! read evc for the local data structure
       !
       IF ( nks > 1 ) CALL get_buffer(evc, nwordwfc, iunwfc, ik)
       !
       ! transform evc to the EXX data structure
       !
       CALL transform_to_exx(lda, n, nbnd, nbnd, ik, evc, evc_exx, type)
       !
       ! save evc to a buffer
       !
       IF ( nks > 1 ) CALL save_buffer ( evc_exx, nwordwfc_exx, iunwfc_exx, ik )
    END DO
    !
    !-----------------------------------------------------------------------
  END SUBROUTINE transform_evc_to_exx
  !-----------------------------------------------------------------------
  !
  !-----------------------------------------------------------------------
  SUBROUTINE transform_psi_to_exx(lda, n, m, psi)
  !-----------------------------------------------------------------------
    USE wvfct,        ONLY : current_k, npwx, nbnd
    USE mp_exx,       ONLY : negrp, nibands, my_egrp_id, max_ibands
    !
    !
    IMPLICIT NONE
    !
    Integer, INTENT(in) :: lda
    INTEGER, INTENT(in) :: m
    INTEGER, INTENT(inout) :: n
    COMPLEX(DP), INTENT(in) :: psi(lda*npol,m) 
    !
    ! change to the EXX data strucutre
    !
    npwx_local = npwx
    n_local = n
    IF ( .not.allocated(comm_recv) ) THEN
       !
       ! initialize all of the conversion maps and change the data structure
       !
       CALL initialize_local_to_exact_map(lda, m)
    ELSE
       !
       ! change the data structure
       !
       CALL change_data_structure(.TRUE.)
    END IF
    npwx_exx = npwx
    n = ngk_exx(current_k)
    !
    ! get igk
    !
    CALL update_igk(.TRUE.)
    !
    ! get psi_exx
    !
    CALL transform_to_exx(lda, n, m, max_ibands, &
         current_k, psi, psi_exx, 0)
    !
    ! zero hpsi_exx
    !
    hpsi_exx = 0.d0
    !
    !-----------------------------------------------------------------------
  END SUBROUTINE transform_psi_to_exx
  !-----------------------------------------------------------------------
  !
  !-----------------------------------------------------------------------
  SUBROUTINE transform_hpsi_to_local(lda, n, m, hpsi)
  !-----------------------------------------------------------------------
    USE mp_exx,         ONLY : iexx_istart, iexx_iend, my_egrp_id
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(in) :: lda
    INTEGER, INTENT(in) :: m
    INTEGER, INTENT(inout) :: n
    COMPLEX(DP), INTENT(out) :: hpsi(lda_original*npol,m)
    INTEGER :: m_exx
    !
    ! change to the local data structure
    !
    CALL change_data_structure(.FALSE.)
    !
    ! get igk
    !
    CALL update_igk(.FALSE.)
    n = n_local
    !
    ! transform hpsi_exx to the local data structure
    !
    m_exx = iexx_iend(my_egrp_id+1) - iexx_istart(my_egrp_id+1) + 1
    CALL transform_to_local(m,m_exx,hpsi_exx,hpsi)
    !
    !-----------------------------------------------------------------------
  END SUBROUTINE transform_hpsi_to_local
  !-----------------------------------------------------------------------
  !
  !-----------------------------------------------------------------------
  SUBROUTINE initialize_local_to_exact_map(lda, m)
  !-----------------------------------------------------------------------
    !
    ! determine the mapping between the local and EXX data structures
    !
    USE wvfct,          ONLY : npwx, nbnd
    USE klist,          ONLY : nks, igk_k
    USE mp_exx,         ONLY : nproc_egrp, negrp, my_egrp_id, me_egrp, &
                               intra_egrp_comm, inter_egrp_comm, &
                               ibands, nibands, init_index_over_band, &
                               iexx_istart, iexx_iend, max_ibands
    USE mp_pools,       ONLY : nproc_pool, me_pool, intra_pool_comm
    USE mp,             ONLY : mp_sum
    USE gvect,          ONLY : ig_l2g
    USE uspp,           ONLY : nkb
    !
    !
    IMPLICIT NONE
    !
    Integer :: lda
    INTEGER :: n, m
    INTEGER, ALLOCATABLE :: local_map(:,:), exx_map(:,:)
    INTEGER, ALLOCATABLE :: l2e_map(:,:), e2l_map(:,:,:)
    INTEGER, ALLOCATABLE :: psi_source(:), psi_source_exx(:,:)
    INTEGER :: current_index
    INTEGER :: i, j, k, ik, im, ig, count, iproc, prev, iegrp
    INTEGER :: total_lda(nks), prev_lda(nks)
    INTEGER :: total_lda_exx(nks), prev_lda_exx(nks)
    INTEGER :: n_local
    INTEGER :: lda_max_local, lda_max_exx
    INTEGER :: max_lda_egrp
    INTEGER :: request_send(nproc_egrp), request_recv(nproc_egrp)
    INTEGER :: ierr
    INTEGER :: egrp_base, total_lda_egrp(nks), prev_lda_egrp(nks)
    INTEGER :: igk_loc(npwx)
    !
    ! initialize the pair assignments
    !
    CALL init_index_over_band(inter_egrp_comm,nbnd,m)
    !
    ! allocate bookeeping arrays
    !
    IF ( .not.allocated(comm_recv) ) THEN
       ALLOCATE(comm_recv(nproc_egrp,nks),comm_send(nproc_egrp,nks))
    END IF
    IF ( .not.allocated(lda_local) ) THEN
       ALLOCATE(lda_local(nproc_pool,nks))
       ALLOCATE(lda_exx(nproc_egrp,nks))
    END IF
    !
    ! store the original values of lda and n
    !
    lda_original = lda
    n_original = n
    !
    ! construct the local map
    !
    lda_local = 0
    DO ik = 1, nks
       igk_loc = igk_k(:,ik)
       n = 0
       DO i = 1, size(igk_loc)
          IF(igk_loc(i).gt.0)n = n + 1
       END DO
       lda_local(me_pool+1,ik) = n
       CALL mp_sum(lda_local(:,ik),intra_pool_comm)
       total_lda(ik) = sum(lda_local(:,ik))
       prev_lda(ik) = sum(lda_local(1:me_pool,ik))
    END DO
    ALLOCATE(local_map(maxval(total_lda),nks))
    local_map = 0
    DO ik = 1, nks
       local_map(prev_lda(ik)+1:prev_lda(ik)+lda_local(me_pool+1,ik),ik) = &
            ig_l2g(igk_k(1:lda_local(me_pool+1,ik),ik))
    END DO
    CALL mp_sum(local_map,intra_pool_comm)
    !
    !-----------------------------------------!
    ! Switch to the exx data structure        !
    !-----------------------------------------!
    !
    CALL change_data_structure(.TRUE.)
    !
    ! construct the exx map
    !
    lda_exx = 0
    DO ik = 1, nks
       n = 0
       DO i = 1, size(igk_exx(:,ik))
          IF(igk_exx(i,ik).gt.0)n = n + 1
       END DO
       lda_exx(me_egrp+1,ik) = n
       CALL mp_sum(lda_exx(:,ik),intra_egrp_comm)
       total_lda_exx(ik) = sum(lda_exx(:,ik))
       prev_lda_exx(ik) = sum(lda_exx(1:me_egrp,ik))
    END DO
    ALLOCATE(exx_map(maxval(total_lda_exx),nks))
    exx_map = 0
    DO ik = 1, nks
       exx_map(prev_lda_exx(ik)+1:prev_lda_exx(ik)+lda_exx(me_egrp+1,ik),ik) = &
            ig_l2g(igk_exx(1:lda_exx(me_egrp+1,ik),ik))    
    END DO
    CALL mp_sum(exx_map,intra_egrp_comm)
    !
    ! construct the l2e_map
    !
    allocate( l2e_map(maxval(total_lda_exx),nks) )
    l2e_map = 0
    DO ik = 1, nks
       DO ig = 1, lda_exx(me_egrp+1,ik)
          DO j=1, total_lda(ik)
             IF( local_map(j,ik).EQ.exx_map(ig+prev_lda_exx(ik),ik) ) exit
          END DO
          l2e_map(ig+prev_lda_exx(ik),ik) = j
       END DO
    END DO
    CALL mp_sum(l2e_map,intra_egrp_comm)
    !
    ! plan communication for the data structure change
    !
    lda_max_local = maxval(lda_local)
    lda_max_exx = maxval(lda_exx)
    allocate(psi_source(maxval(total_lda_exx)))
    !
    DO ik = 1, nks
       !
       ! determine which task each value will come from
       !
       psi_source = 0
       DO ig = 1, lda_exx(me_egrp+1,ik)
          j = 1
          DO i = 1, nproc_pool
             j = j + lda_local(i,ik)
             IF( j.gt.l2e_map(ig+prev_lda_exx(ik),ik) ) exit
          END DO
          psi_source(ig+prev_lda_exx(ik)) = i-1
       END DO
       CALL mp_sum(psi_source,intra_egrp_comm)
       !
       ! allocate communication packets to receive psi and hpsi
       !
       DO iproc=0, nproc_egrp-1
          !
          ! determine how many values need to come from iproc
          !
          count = 0
          DO ig=1, lda_exx(me_egrp+1,ik)
             IF ( MODULO(psi_source(ig+prev_lda_exx(ik)),nproc_egrp)&
                  .eq.iproc ) THEN
                count = count + 1
             END IF
          END DO
          !
          ! allocate the communication packet
          !
          comm_recv(iproc+1,ik)%size = count
          IF (count.gt.0) THEN
             IF (.not.ALLOCATED(comm_recv(iproc+1,ik)%msg)) THEN
                ALLOCATE(comm_recv(iproc+1,ik)%indices(count))
                ALLOCATE(comm_recv(iproc+1,ik)%msg(count,npol,max_ibands+2))
             END IF
          END IF
          !
          ! determine which values need to come from iproc
          !
          count = 0
          DO ig=1, lda_exx(me_egrp+1,ik)
             IF ( MODULO(psi_source(ig+prev_lda_exx(ik)),nproc_egrp)&
                  .eq.iproc ) THEN
                count = count + 1
                comm_recv(iproc+1,ik)%indices(count) = ig
             END IF
          END DO
          !
       END DO
       !
       ! allocate communication packets to send psi and hpsi
       !
       prev = 0
       DO iproc=0, nproc_egrp-1
          !
          ! determine how many values need to be sent to iproc
          !
          count = 0
          DO ig=1, lda_exx(iproc+1,ik)
             IF ( MODULO(psi_source(ig+prev),nproc_egrp).eq.me_egrp ) THEN
                count = count + 1
             END IF
          END DO
          !
          ! allocate the communication packet
          !
          comm_send(iproc+1,ik)%size = count
          IF (count.gt.0) THEN
             IF (.not.ALLOCATED(comm_send(iproc+1,ik)%msg)) THEN
                ALLOCATE(comm_send(iproc+1,ik)%indices(count))
                ALLOCATE(comm_send(iproc+1,ik)%msg(count,npol,max_ibands+2))
             END IF
          END IF
          !
          ! determine which values need to be sent to iproc
          !
          count = 0
          DO ig=1, lda_exx(iproc+1,ik)
             IF ( MODULO(psi_source(ig+prev),nproc_egrp).eq.me_egrp ) THEN
                count = count + 1
                comm_send(iproc+1,ik)%indices(count) = l2e_map(ig+prev,ik)
             END IF
          END DO
          !
          prev = prev + lda_exx(iproc+1,ik)
          !
       END DO
       !
    END DO
    !
    ! allocate psi_exx and hpsi_exx
    !
    IF(allocated(psi_exx))DEALLOCATE(psi_exx)
    ALLOCATE(psi_exx(npwx*npol, max_ibands ))
    IF(allocated(hpsi_exx))DEALLOCATE(hpsi_exx)
    ALLOCATE(hpsi_exx(npwx*npol, max_ibands ))
    !
    ! allocate communication arrays for the exx to local transformation
    !
    IF ( .not.allocated(comm_recv_reverse) ) THEN
       ALLOCATE(comm_recv_reverse(nproc_egrp,nks))
       ALLOCATE(comm_send_reverse(nproc_egrp,negrp,nks))
    END IF
    !
    egrp_base = my_egrp_id*nproc_egrp
    DO ik = 1, nks
       total_lda_egrp(ik) = &
            sum( lda_local(egrp_base+1:(egrp_base+nproc_egrp),ik) )
       prev_lda_egrp(ik) = &
            sum( lda_local(egrp_base+1:(egrp_base+me_egrp),ik) )
    END DO
    !
    max_lda_egrp = 0
    DO j = 1, negrp
       DO ik = 1, nks
          max_lda_egrp = max(max_lda_egrp, &
               sum( lda_local((j-1)*nproc_egrp+1:j*nproc_egrp,ik) ) )
       END DO
    END DO
    !
    ! determine the e2l_map
    !
    allocate( e2l_map(max_lda_egrp,nks,negrp) )
    e2l_map = 0
    DO ik = 1, nks
       DO ig = 1, lda_local(me_pool+1,ik)
          DO j=1, total_lda_exx(ik)
             IF( local_map(ig+prev_lda(ik),ik).EQ.exx_map(j,ik) ) exit
          END DO
          e2l_map(ig+prev_lda_egrp(ik),ik,my_egrp_id+1) = j
       END DO
    END DO
    CALL mp_sum(e2l_map(:,:,my_egrp_id+1),intra_egrp_comm)
    CALL mp_sum(e2l_map,inter_egrp_comm)
    !
    ! plan communication for the local to EXX data structure transformation
    !
    allocate(psi_source_exx( max_lda_egrp, negrp ))
    DO ik = 1, nks
       !
       ! determine where each value is coming from
       !
       psi_source_exx = 0
       DO ig = 1, lda_local(me_pool+1,ik)
          j = 1
          DO i = 1, nproc_egrp
             j = j + lda_exx(i,ik)
             IF( j.gt.e2l_map(ig+prev_lda_egrp(ik),ik,my_egrp_id+1) ) exit
          END DO
          psi_source_exx(ig+prev_lda_egrp(ik),my_egrp_id+1) = i-1
       END DO
       CALL mp_sum(psi_source_exx(:,my_egrp_id+1),intra_egrp_comm)
       CALL mp_sum(psi_source_exx,inter_egrp_comm)
       !
       ! allocate communication packets to receive psi and hpsi (reverse)
       !
       DO iegrp=my_egrp_id+1, my_egrp_id+1
          DO iproc=0, nproc_egrp-1
             !
             ! determine how many values need to come from iproc
             !
             count = 0
             DO ig=1, lda_local(me_pool+1,ik)
                IF ( psi_source_exx(ig+prev_lda_egrp(ik),iegrp).eq.iproc ) THEN
                   count = count + 1
                END IF
             END DO
             !
             ! allocate the communication packet
             !
             comm_recv_reverse(iproc+1,ik)%size = count
             IF (count.gt.0) THEN
                IF (.not.ALLOCATED(comm_recv_reverse(iproc+1,ik)%msg)) THEN
                   ALLOCATE(comm_recv_reverse(iproc+1,ik)%indices(count))
                   ALLOCATE(comm_recv_reverse(iproc+1,ik)%msg(count,npol,m+2))
                END IF
             END IF
             !
             ! determine which values need to come from iproc
             !
             count = 0
             DO ig=1, lda_local(me_pool+1,ik)
                IF ( psi_source_exx(ig+prev_lda_egrp(ik),iegrp).eq.iproc ) THEN
                   count = count + 1
                   comm_recv_reverse(iproc+1,ik)%indices(count) = ig
                END IF
             END DO
             
          END DO
       END DO
       !
       ! allocate communication packets to send psi and hpsi
       !
       DO iegrp=1, negrp
          prev = 0
          DO iproc=0, nproc_egrp-1
             !
             ! determine how many values need to be sent to iproc
             !
             count = 0
             DO ig=1, lda_local(iproc+(iegrp-1)*nproc_egrp+1,ik)
                IF ( psi_source_exx(ig+prev,iegrp).eq.me_egrp ) THEN
                   count = count + 1
                END IF
             END DO
             !
             ! allocate the communication packet
             !
             comm_send_reverse(iproc+1,iegrp,ik)%size = count
             IF (count.gt.0) THEN
                IF (.not.ALLOCATED(comm_send_reverse(iproc+1,iegrp,ik)%msg))THEN
                   ALLOCATE(comm_send_reverse(iproc+1,iegrp,ik)%indices(count))
                   ALLOCATE(comm_send_reverse(iproc+1,iegrp,ik)%msg(count,npol,&
                        iexx_iend(my_egrp_id+1)-iexx_istart(my_egrp_id+1)+3))
                END IF
             END IF
             !
             ! determine which values need to be sent to iproc
             !
             count = 0
             DO ig=1, lda_local(iproc+(iegrp-1)*nproc_egrp+1,ik)
                IF ( psi_source_exx(ig+prev,iegrp).eq.me_egrp ) THEN
                   count = count + 1
                   comm_send_reverse(iproc+1,iegrp,ik)%indices(count) = &
                        e2l_map(ig+prev,ik,iegrp)
                END IF
             END DO
             !
             prev = prev + lda_local( (iegrp-1)*nproc_egrp+iproc+1, ik )
             !
          END DO
       END DO
       !
    END DO
    !
    ! deallocate arrays
    !
    DEALLOCATE( local_map, exx_map )
    DEALLOCATE( l2e_map, e2l_map )
    DEALLOCATE( psi_source, psi_source_exx )
    !
    !-----------------------------------------------------------------------
  END SUBROUTINE initialize_local_to_exact_map
  !-----------------------------------------------------------------------
  !
  !-----------------------------------------------------------------------
  SUBROUTINE transform_to_exx(lda, n, m, m_out, ik, psi, psi_out, type)
  !-----------------------------------------------------------------------
    !
    ! transform psi into the EXX data structure, and place the result in psi_out
    !
    USE wvfct,        ONLY : nbnd
    USE mp,           ONLY : mp_sum
    USE mp_pools,     ONLY : nproc_pool, me_pool
    USE mp_exx,       ONLY : intra_egrp_comm, inter_egrp_comm, &
         nproc_egrp, me_egrp, negrp, my_egrp_id, nibands, ibands, &
         max_ibands, all_start, all_end
    USE parallel_include
    USE klist,        ONLY : xk, wk, nkstot, nks, qnorm
    !
    !
    IMPLICIT NONE
    !
    Integer :: lda
    INTEGER :: n, m, m_out
    COMPLEX(DP) :: psi(npwx_local*npol,m) 
    COMPLEX(DP) :: psi_out(npwx_exx*npol,m_out)
    INTEGER, INTENT(in) :: type

    COMPLEX(DP), ALLOCATABLE :: psi_work(:,:,:,:), psi_gather(:,:)
    INTEGER :: i, j, im, iproc, ig, ik, iegrp
    INTEGER :: prev, lda_max_local

    INTEGER :: request_send(nproc_egrp), request_recv(nproc_egrp)
    INTEGER :: ierr
    INTEGER :: current_ik
    !INTEGER, EXTERNAL :: find_current_k
    INTEGER :: recvcount(negrp)
    INTEGER :: displs(negrp)
#if defined(__MPI)
    INTEGER :: istatus(MPI_STATUS_SIZE)
#endif
    INTEGER :: requests(max_ibands+2,negrp)
    !
    INTEGER :: ipol, my_lda, lda_offset, count
    lda_max_local = maxval(lda_local)
    current_ik = ik
    !
    !-------------------------------------------------------!
    !Communication Part 1
    !-------------------------------------------------------!
    !
    allocate(psi_work(lda_max_local,npol,m,negrp))
    allocate(psi_gather(lda_max_local*npol,m))
    DO im=1, m
       my_lda = lda_local(me_pool+1,current_ik)
       DO ipol=1, npol
          lda_offset = lda_max_local*(ipol-1)
          DO ig=1, my_lda
             psi_gather(lda_offset+ig,im) = psi(npwx_local*(ipol-1)+ig,im)
          END DO
       END DO
    END DO
    recvcount = lda_max_local*npol
    count = lda_max_local*npol
#if defined (__MPI_NONBLOCKING)
    IF ( type.eq.0 ) THEN
       DO iegrp=1, negrp
          displs(iegrp) = (iegrp-1)*(count*m)
       END DO
       DO iegrp=1, negrp
          DO im=1, nibands(iegrp)
             IF ( my_egrp_id.eq.(iegrp-1) ) THEN
                DO j=1, negrp
                   displs(j) = (j-1)*(count*m) + count*(ibands(im,iegrp)-1)
                END DO
             END IF
#if defined(__MPI)
             CALL MPI_IGATHERV( psi_gather(:, ibands(im,iegrp) ), &
                  count, MPI_DOUBLE_COMPLEX, &
                  psi_work, &
                  recvcount, displs, MPI_DOUBLE_COMPLEX, &
                  iegrp-1, &
                  inter_egrp_comm, requests(im,iegrp), ierr )
#endif
          END DO
       END DO

    ELSE IF(type.eq.1) THEN
#elif defined(__MPI)
       CALL MPI_ALLGATHER( psi_gather, &
            count*m, MPI_DOUBLE_COMPLEX, &
            psi_work, &
            count*m, MPI_DOUBLE_COMPLEX, &
            inter_egrp_comm, ierr )
#endif

#if defined (__MPI_NONBLOCKING)
    ELSE IF(type.eq.2) THEN !evc2

       DO iegrp=1, negrp
          displs(iegrp) = (iegrp-1)*(count*m)
       END DO
       DO iegrp=1, negrp
          DO im=1, all_end(iegrp) - all_start(iegrp) + 1
             !
             IF(all_start(iegrp).eq.0) CYCLE
             !
             IF ( my_egrp_id.eq.(iegrp-1) ) THEN
                DO j=1, negrp
                   displs(j) = (j-1)*(count*m) + count*(im+all_start(iegrp)-2)
                END DO
             END IF
#if defined(__MPI)
             CALL MPI_IGATHERV( psi_gather(:, im+all_start(iegrp)-1 ), &
                  count, MPI_DOUBLE_COMPLEX, &
                  psi_work, &
                  recvcount, displs, MPI_DOUBLE_COMPLEX, &
                  iegrp-1, &
                  inter_egrp_comm, requests(im,iegrp), ierr )
#endif
          END DO
       END DO
          
    END IF
    !
    IF(type.eq.0)THEN
       DO iegrp=1, negrp
          DO im=1, nibands(iegrp)
             CALL MPI_WAIT(requests(im,iegrp), istatus, ierr)
          END DO
       END DO
    ELSEIF(type.eq.2)THEN
       DO iegrp=1, negrp
          DO im=1, all_end(iegrp) - all_start(iegrp) + 1
             IF(all_start(iegrp).eq.0) CYCLE
             CALL MPI_WAIT(requests(im,iegrp), istatus, ierr)
          END DO
       END DO
    END IF
#endif
    !
    !-------------------------------------------------------!
    !Communication Part 2
    !-------------------------------------------------------!
    !
    !send communication packets
    !
    DO iproc=0, nproc_egrp-1
       IF ( comm_send(iproc+1,current_ik)%size.gt.0) THEN
          DO i=1, comm_send(iproc+1,current_ik)%size
             ig = comm_send(iproc+1,current_ik)%indices(i)
             !
             !determine which egrp this corresponds to
             !
             prev = 0
             DO j=1, nproc_pool
                IF ((prev+lda_local(j,current_ik)).ge.ig) THEN 
                   ig = ig - prev
                   exit
                END IF
                prev = prev + lda_local(j,current_ik)
             END DO
             !
             ! prepare the message
             !
             DO ipol=1, npol
             IF ( type.eq.0 ) THEN !psi or hpsi
                DO im=1, nibands(my_egrp_id+1)
                   comm_send(iproc+1,current_ik)%msg(i,ipol,im) = &
                        psi_work(ig,ipol,ibands(im,my_egrp_id+1),1+(j-1)/nproc_egrp)
                END DO
             ELSE IF (type.eq.1) THEN !evc
                DO im=1, m
                   comm_send(iproc+1,current_ik)%msg(i,ipol,im) = &
                        psi_work(ig,ipol,im,1+(j-1)/nproc_egrp)
                END DO
             ELSE IF ( type.eq.2 ) THEN !evc2
                IF(all_start(my_egrp_id+1).gt.0) THEN
                   DO im=1, all_end(my_egrp_id+1) - all_start(my_egrp_id+1) + 1
                      comm_send(iproc+1,current_ik)%msg(i,ipol,im) = &
                           psi_work(ig,ipol,im+all_start(my_egrp_id+1)-1,1+(j-1)/nproc_egrp)
                   END DO
                END IF
             END IF
             END DO
             !
          END DO
          !
          ! send the message
          !
#if defined(__MPI)
          IF ( type.eq.0 ) THEN !psi or hpsi
             CALL MPI_ISEND( comm_send(iproc+1,current_ik)%msg, &
                  comm_send(iproc+1,current_ik)%size*npol*nibands(my_egrp_id+1), &
                  MPI_DOUBLE_COMPLEX, &
                  iproc, 100+iproc*nproc_egrp+me_egrp, &
                  intra_egrp_comm, request_send(iproc+1), ierr )
          ELSE IF (type.eq.1) THEN !evc
             CALL MPI_ISEND( comm_send(iproc+1,current_ik)%msg, &
                  comm_send(iproc+1,current_ik)%size*npol*m, MPI_DOUBLE_COMPLEX, &
                  iproc, 100+iproc*nproc_egrp+me_egrp, &
                  intra_egrp_comm, request_send(iproc+1), ierr )
          ELSE IF (type.eq.2) THEN !evc2
             CALL MPI_ISEND( comm_send(iproc+1,current_ik)%msg, &
                  comm_send(iproc+1,current_ik)%size*npol*(all_end(my_egrp_id+1)-all_start(my_egrp_id+1)+1), &
                  MPI_DOUBLE_COMPLEX, &
                  iproc, 100+iproc*nproc_egrp+me_egrp, &
                  intra_egrp_comm, request_send(iproc+1), ierr )
          END IF
#endif
          !
       END IF
    END DO
    !
    ! begin receiving the messages
    !
    DO iproc=0, nproc_egrp-1
       IF ( comm_recv(iproc+1,current_ik)%size.gt.0) THEN
          !
          ! receive the message
          !
#if defined(__MPI)
          IF (type.eq.0) THEN !psi or hpsi
             CALL MPI_IRECV( comm_recv(iproc+1,current_ik)%msg, &
                  comm_recv(iproc+1,current_ik)%size*npol*nibands(my_egrp_id+1), &
                  MPI_DOUBLE_COMPLEX, &
                  iproc, 100+me_egrp*nproc_egrp+iproc, &
                  intra_egrp_comm, request_recv(iproc+1), ierr )
          ELSE IF (type.eq.1) THEN !evc
             CALL MPI_IRECV( comm_recv(iproc+1,current_ik)%msg, &
                  comm_recv(iproc+1,current_ik)%size*npol*m, MPI_DOUBLE_COMPLEX, &
                  iproc, 100+me_egrp*nproc_egrp+iproc, &
                  intra_egrp_comm, request_recv(iproc+1), ierr )
          ELSE IF (type.eq.2) THEN !evc2
             CALL MPI_IRECV( comm_recv(iproc+1,current_ik)%msg, &
                  comm_recv(iproc+1,current_ik)%size*npol*(all_end(my_egrp_id+1)-all_start(my_egrp_id+1)+1), &
                  MPI_DOUBLE_COMPLEX, &
                  iproc, 100+me_egrp*nproc_egrp+iproc, &
                  intra_egrp_comm, request_recv(iproc+1), ierr )
          END IF
#endif
          !
       END IF
    END DO
    !
    ! assign psi_out
    !
    DO iproc=0, nproc_egrp-1
       IF ( comm_recv(iproc+1,current_ik)%size.gt.0 ) THEN
          !
          ! wait for the message to be received
          !
#if defined(__MPI)
          CALL MPI_WAIT(request_recv(iproc+1), istatus, ierr)
#endif
          !
          DO i=1, comm_recv(iproc+1,current_ik)%size
             ig = comm_recv(iproc+1,current_ik)%indices(i)
             !
             ! place the message into the correct elements of psi_out
             !
             IF (type.eq.0) THEN !psi or hpsi
                DO im=1, nibands(my_egrp_id+1)
                   DO ipol=1, npol
                      psi_out(ig+npwx_exx*(ipol-1),im) = &
                           comm_recv(iproc+1,current_ik)%msg(i,ipol,im)
                   END DO
                END DO
             ELSE IF (type.eq.1) THEN !evc
                DO im=1, m
                   DO ipol=1, npol
                      psi_out(ig+npwx_exx*(ipol-1),im) = &
                           comm_recv(iproc+1,current_ik)%msg(i,ipol,im)
                   END DO
                END DO
             ELSE IF (type.eq.2) THEN !evc2
                DO im=1, all_end(my_egrp_id+1) - all_start(my_egrp_id+1) + 1
                   DO ipol=1, npol
                      psi_out(ig+npwx_exx*(ipol-1),im) = comm_recv(iproc+1,current_ik)%msg(i,ipol,im)
                   END DO
                END DO
             END IF
             !
          END DO
          !
       END IF
    END DO
    !
    ! wait for everything to finish sending
    !
#if defined(__MPI)
    DO iproc=0, nproc_egrp-1
       IF ( comm_send(iproc+1,current_ik)%size.gt.0 ) THEN
          CALL MPI_WAIT(request_send(iproc+1), istatus, ierr)
       END IF
    END DO
#endif
    !
    ! deallocate arrays
    !
    DEALLOCATE( psi_work, psi_gather )
    !
    !-----------------------------------------------------------------------
  END SUBROUTINE transform_to_exx
  !-----------------------------------------------------------------------
  !
  !-----------------------------------------------------------------------
  SUBROUTINE change_data_structure(is_exx)
  !-----------------------------------------------------------------------
    !
    ! change between the local and EXX data structures
    ! is_exx = .TRUE. - change to the EXX data structure
    ! is_exx = .FALSE. - change to the local data strucutre
    !
    USE cell_base,      ONLY : at, bg, tpiba2
    USE cellmd,         ONLY : lmovecell
    USE wvfct,          ONLY : npwx
    USE gvect,          ONLY : gcutm, ig_l2g, g, gg, nl, nlm, ngm, ngm_g, mill, &
                               gstart, gvect_init, deallocate_gvect_exx
    USE gvecs,          ONLY : gcutms, ngms, ngms_g, nls, nlsm, gvecs_init, &
                               deallocate_gvecs
    USE gvecw,          ONLY : gkcut, ecutwfc, gcutw
    USE klist,          ONLY : xk, nks, ngk
    USE mp_bands,       ONLY : intra_bgrp_comm, ntask_groups
    USE mp_exx,         ONLY : intra_egrp_comm, me_egrp, exx_mode, nproc_egrp, &
                               negrp, root_egrp
    USE io_global,      ONLY : stdout
    USE fft_base,       ONLY : dfftp, dffts, dtgs, smap, fft_base_info
    USE fft_types,      ONLY : fft_type_init
    USE recvec_subs,    ONLY : ggen
    USE task_groups,    ONLY : task_groups_init
    !
    !
    IMPLICIT NONE
    !
    LOGICAL, intent(in) :: is_exx
    INTEGER, EXTERNAL  :: n_plane_waves
    COMPLEX(DP), ALLOCATABLE :: work_space(:)
    INTEGER :: ik, i
    INTEGER :: ngm_, ngs_, ngw_
    LOGICAL exst
#if defined (__MPI)
  LOGICAL :: lpara = .true.
#else
  LOGICAL :: lpara = .false.
#endif
    !
    IF (negrp.eq.1) RETURN
    !
    IF (first_data_structure_change) THEN
       allocate( ig_l2g_loc(ngm), g_loc(3,ngm), gg_loc(ngm) )
       allocate( mill_loc(3,ngm), nl_loc(ngm) )
       allocate( nls_loc(size(nls)) )
       allocate( nlm_loc(size(nlm)) )
       allocate( nlsm_loc(size(nlsm)) )
       ig_l2g_loc = ig_l2g
       g_loc = g
       gg_loc = gg
       mill_loc = mill
       nl_loc = nl
       nls_loc = nls
       nlm_loc = nlm
       nlsm_loc = nlsm
       ngm_loc = ngm
       ngm_g_loc = ngm_g
       gstart_loc = gstart
       ngms_loc = ngms
       ngms_g_loc = ngms_g
    END IF
    !
    ! generate the gvectors for the new data structure
    !
    IF (is_exx) THEN
       exx_mode = 1
       IF(first_data_structure_change)THEN
          dfftp_loc = dfftp
          dffts_loc = dffts

          CALL fft_type_init( dffts_exx, smap_exx, "wave", gamma_only, &
               lpara, intra_egrp_comm, at, bg, gkcut, gcutms/gkcut, &
               ntask_groups )
          CALL fft_type_init( dfftp_exx, smap_exx, "rho", gamma_only, &
               lpara, intra_egrp_comm, at, bg,  gcutm )
          CALL fft_base_info( ionode, stdout )
          ngs_ = dffts_exx%ngl( dffts_exx%mype + 1 )
          ngm_ = dfftp_exx%ngl( dfftp_exx%mype + 1 )
          IF( gamma_only ) THEN
             ngs_ = (ngs_ + 1)/2
             ngm_ = (ngm_ + 1)/2
          END IF
          dfftp = dfftp_exx
          dffts = dffts_exx
          ngm = ngm_
          ngms = ngs_
       ELSE
          dfftp = dfftp_exx
          dffts = dffts_exx
          ngm = ngm_exx
          ngms = ngms_exx
       END IF
       call deallocate_gvect_exx()
       call deallocate_gvecs()
       call gvect_init( ngm , intra_egrp_comm )
       call gvecs_init( ngms , intra_egrp_comm )
    ELSE
       exx_mode = 2
       dfftp = dfftp_loc
       dffts = dffts_loc
       ngm = ngm_loc
       ngms = ngms_loc
       call deallocate_gvect_exx()
       call deallocate_gvecs()
       call gvect_init( ngm , intra_bgrp_comm )
       call gvecs_init( ngms , intra_bgrp_comm )
       exx_mode = 0
    END IF
    !
    IF (first_data_structure_change) THEN
       CALL ggen( gamma_only, at, bg, intra_egrp_comm, no_global_sort = .FALSE. )
       allocate( ig_l2g_exx(ngm), g_exx(3,ngm), gg_exx(ngm) )
       allocate( mill_exx(3,ngm), nl_exx(ngm) )
       allocate( nls_exx(size(nls)) )
       allocate( nlm_exx(size(nlm) ) )
       allocate( nlsm_exx(size(nlsm) ) )
       ig_l2g_exx = ig_l2g
       g_exx = g
       gg_exx = gg
       mill_exx = mill
       nl_exx = nl
       nls_exx = nls
       nlm_exx = nlm
       nlsm_exx = nlsm
       ngm_exx = ngm
       ngm_g_exx = ngm_g
       gstart_exx = gstart
       ngms_exx = ngms
       ngms_g_exx = ngms_g
    ELSE IF ( is_exx ) THEN
       ig_l2g = ig_l2g_exx
       g = g_exx
       gg = gg_exx
       mill = mill_exx
       nl = nl_exx
       nls = nls_exx
       nlm = nlm_exx
       nlsm = nlsm_exx
       ngm = ngm_exx
       ngm_g = ngm_g_exx
       gstart = gstart_exx
       ngms = ngms_exx
       ngms_g = ngms_g_exx
    ELSE ! not is_exx
       ig_l2g = ig_l2g_loc
       g = g_loc
       gg = gg_loc
       mill = mill_loc
       nl = nl_loc
       nls = nls_loc
       nlm = nlm_loc
       nlsm = nlsm_loc
       ngm = ngm_loc
       ngm_g = ngm_g_loc
       gstart = gstart_loc
       ngms = ngms_loc
       ngms_g = ngms_g_loc
    END IF
    !
    ! get npwx and ngk
    !
    IF ( is_exx.and.npwx_exx.gt.0 ) THEN
       npwx = npwx_exx
       ngk = ngk_exx
    ELSE IF ( .not.is_exx.and.npwx_local.gt.0 ) THEN
       npwx = npwx_local
       ngk = ngk_local
    ELSE
       npwx = n_plane_waves (gcutw, nks, xk, g, ngm)
    END IF
    !
    ! get igk
    !
    IF( first_data_structure_change ) THEN
       allocate(igk_exx(npwx,nks),work_space(npwx))
       first_data_structure_change = .FALSE.
       IF ( nks.eq.1 ) THEN
          CALL gk_sort( xk, ngm, g, ecutwfc / tpiba2, ngk, igk_exx, work_space )
       END IF
       IF ( nks > 1 ) THEN
          !
          DO ik = 1, nks
             CALL gk_sort( xk(1,ik), ngm, g, ecutwfc / tpiba2, ngk(ik), &
                  igk_exx(1,ik), work_space )
          END DO
          !
       END IF
       DEALLOCATE( work_space )
    END IF
    !
    ! generate ngl and igtongl
    !
    CALL gshells( lmovecell )
    !
    !-----------------------------------------------------------------------
  END SUBROUTINE change_data_structure
  !-----------------------------------------------------------------------
  !
  !-----------------------------------------------------------------------
  SUBROUTINE update_igk(is_exx)
  !-----------------------------------------------------------------------
    USE cell_base,      ONLY : tpiba2
    USE gvect,          ONLY : ngm, g
    USE gvecw,          ONLY : ecutwfc
    USE wvfct,          ONLY : npwx, npw, current_k
    USE klist,          ONLY : xk, igk_k
    USE mp_exx,         ONLY : negrp
    !
    !
    IMPLICIT NONE
    !
    LOGICAL, intent(in) :: is_exx
    COMPLEX(DP), ALLOCATABLE :: work_space(:)
    INTEGER :: comm
    INTEGER :: ik, i
    LOGICAL exst
    !
    IF (negrp.eq.1) RETURN
    !
    ! get igk
    !
    allocate(work_space(npwx))
    ik = current_k
    IF(is_exx) THEN
       CALL gk_sort( xk(1,ik), ngm, g, ecutwfc / tpiba2, npw, igk_exx(1,ik), &
            work_space )
    ELSE
       CALL gk_sort( xk(1,ik), ngm, g, ecutwfc / tpiba2, npw, igk_k(1,ik), &
            work_space )
    END IF
    !
    DEALLOCATE( work_space )
    !
    !-----------------------------------------------------------------------
  END SUBROUTINE update_igk
  !-----------------------------------------------------------------------
  !
  !-----------------------------------------------------------------------
  SUBROUTINE communicate_exxbuff (ipair, request_send, request_recv)
  !-----------------------------------------------------------------------
    USE mp_exx,       ONLY : iexx_start, iexx_end, inter_egrp_comm, &
                               intra_egrp_comm, my_egrp_id, negrp, &
                               max_pairs, egrp_pairs
    USE parallel_include
    USE io_global,      ONLY : stdout
    INTEGER, intent(in)      :: ipair
    INTEGER                  :: nrxxs
    INTEGER                  :: request_send, request_recv
    INTEGER                  :: dest, sender, ierr, jnext, jnext_dest
#if defined(__MPI)
    INTEGER :: istatus(MPI_STATUS_SIZE)
#endif
    !
    nrxxs= exx_fft%dfftt%nnr
    !
    IF (ipair.lt.max_pairs) THEN
       !
       IF (ipair.gt.1) THEN
#if defined(__MPI)
          CALL MPI_WAIT(request_send, istatus, ierr)
          CALL MPI_WAIT(request_recv, istatus, ierr)
#endif
       END IF
       !
       sender = my_egrp_id + 1
       IF (sender.ge.negrp) sender = 0
       jnext = egrp_pairs(2,ipair+1,sender+1)
       !
       dest = my_egrp_id - 1
       IF (dest.lt.0) dest = negrp - 1
       jnext_dest = egrp_pairs(2,ipair+1,dest+1)
       !
#if defined(__MPI)
       CALL MPI_ISEND( exxbuff(:,:,jnext_dest), nrxxs*npol*nqs, &
            MPI_DOUBLE_COMPLEX, dest, 101, inter_egrp_comm, request_send, ierr )
#endif
       !
#if defined(__MPI)
       CALL MPI_IRECV( exxbuff(:,:,jnext), nrxxs*npol*nqs, &
            MPI_DOUBLE_COMPLEX, sender, 101, inter_egrp_comm, &
            request_recv, ierr )
#endif
       !
    END IF
    !
    !-----------------------------------------------------------------------
  END SUBROUTINE communicate_exxbuff
  !-----------------------------------------------------------------------
  !
  !-----------------------------------------------------------------------
  SUBROUTINE result_sum (n, m, data)
  !-----------------------------------------------------------------------
    USE parallel_include
    USE mp_exx,       ONLY : iexx_start, iexx_end, inter_egrp_comm, &
                               intra_egrp_comm, my_egrp_id, negrp, &
                               max_pairs, egrp_pairs, max_contributors, &
                               contributed_bands, all_end, &
                               iexx_istart, iexx_iend, band_roots
    USE mp,           ONLY : mp_sum, mp_bcast

    INTEGER, INTENT(in) :: n, m
    COMPLEX(DP), INTENT(inout) :: data(n,m)
#if defined(__MPI)
    INTEGER :: istatus(MPI_STATUS_SIZE)
    COMPLEX(DP), ALLOCATABLE :: recvbuf(:,:)
    COMPLEX(DP) :: data_sum(n,m), test(negrp)
    INTEGER :: im, iegrp, ibuf, i, j, nsending(m)
    INTEGER :: ncontributing(m)
    INTEGER :: contrib_this(negrp,m), displs(negrp,m)

    INTEGER sendcount, sendtype, ierr, root, request(m)
    INTEGER sendc(negrp), sendd(negrp)
#endif
    !
    IF (negrp.eq.1) RETURN
    !
#if defined(__MPI)
    ! gather data onto the correct nodes
    !
    ALLOCATE( recvbuf( n*max_contributors, max(1,iexx_end-iexx_start+1) ) )
    displs = 0
    ibuf = 0
    nsending = 0
    contrib_this = 0
    !
    DO im=1, m
       !
       IF(contributed_bands(im,my_egrp_id+1)) THEN
          sendcount = n
       ELSE
          sendcount = 0
       END IF
       !
       root = band_roots(im)
       !
       IF(my_egrp_id.eq.root) THEN
          !
          ! determine the number of sending processors
          !
          ibuf = ibuf + 1
          ncontributing(im) = 0
          DO iegrp=1, negrp
             IF(contributed_bands(im,iegrp)) THEN
                ncontributing(im) = ncontributing(im) + 1
                contrib_this(iegrp,im) = n
                IF(iegrp.lt.negrp) displs(iegrp+1,im) = displs(iegrp,im) + n
                nsending(im) = nsending(im) + 1
             ELSE
                contrib_this(iegrp,im) = 0
                IF(iegrp.lt.negrp) displs(iegrp+1,im) = displs(iegrp,im)
             END IF
          END DO
       END IF
       !
#if defined(__MPI_NONBLOCKING)
       CALL MPI_IGATHERV(data(:,im), sendcount, MPI_DOUBLE_COMPLEX, &
            recvbuf(:,max(1,ibuf)), contrib_this(:,im), &
            displs(:,im), MPI_DOUBLE_COMPLEX, &
            root, inter_egrp_comm, request(im), ierr)
#else
       CALL MPI_GATHERV(data(:,im), sendcount, MPI_DOUBLE_COMPLEX, &
            recvbuf(:,max(1,ibuf)), contrib_this(:,im), &
            displs(:,im), MPI_DOUBLE_COMPLEX, &
            root, inter_egrp_comm, ierr)
#endif
       !
    END DO
    !
#if defined(__MPI_NONBLOCKING)
    DO im=1, m
       CALL MPI_WAIT(request(im), istatus, ierr)
    END DO
#endif
    !
    ! perform the sum
    !
    DO im=iexx_istart(my_egrp_id+1), iexx_iend(my_egrp_id+1)
       IF(im.eq.0)exit
       data(:,im) = 0._dp
       ibuf = im - iexx_istart(my_egrp_id+1) + 1
       DO j=1, nsending(im)
          DO i=1, n
             data(i,im) = data(i,im) + recvbuf(i+n*(j-1),ibuf)
          END DO
       END DO
    END DO
#endif
    !
    !-----------------------------------------------------------------------
  END SUBROUTINE result_sum
  !-----------------------------------------------------------------------
  !
  !-----------------------------------------------------------------------
  SUBROUTINE transform_to_local(m, m_exx, psi, psi_out)
  !-----------------------------------------------------------------------
    USE mp,           ONLY : mp_sum
    USE mp_pools,     ONLY : nproc_pool, me_pool, intra_pool_comm
    USE mp_exx,       ONLY : intra_egrp_comm, inter_egrp_comm, &
         nproc_egrp, me_egrp, negrp, my_egrp_id, iexx_istart, iexx_iend
    USE parallel_include
    USE klist,        ONLY : xk, wk, nkstot, nks, qnorm
    USE wvfct,        ONLY : current_k
    !
    !
    IMPLICIT NONE
    !
    INTEGER :: m, m_exx
    COMPLEX(DP) :: psi(npwx_exx*npol,m_exx)
    COMPLEX(DP) :: psi_out(npwx_local*npol,m)
    !
    INTEGER :: i, j, im, iproc, ig, ik, current_ik, iegrp
    INTEGER :: prev, lda_max_local, prev_lda_exx
    INTEGER :: my_bands, recv_bands, tag
    !
    INTEGER :: request_send(nproc_egrp,negrp), request_recv(nproc_egrp,negrp)
    INTEGER :: ierr
    INTEGER :: ipol
#if defined(__MPI)
    INTEGER :: istatus(MPI_STATUS_SIZE)
#endif
    !INTEGER, EXTERNAL :: find_current_k
    !
    current_ik = current_k
    prev_lda_exx = sum( lda_exx(1:me_egrp,current_ik) )
    !
    my_bands = iexx_iend(my_egrp_id+1) - iexx_istart(my_egrp_id+1) + 1
    !
    ! send communication packets
    !
    IF ( iexx_istart(my_egrp_id+1).gt.0 ) THEN
       DO iegrp=1, negrp
          DO iproc=0, nproc_egrp-1
             IF ( comm_send_reverse(iproc+1,iegrp,current_ik)%size.gt.0) THEN
                DO i=1, comm_send_reverse(iproc+1,iegrp,current_ik)%size
                   ig = comm_send_reverse(iproc+1,iegrp,current_ik)%indices(i)
                   ig = ig - prev_lda_exx
                   !
                   DO im=1, my_bands
                      DO ipol=1, npol
                         comm_send_reverse(iproc+1,iegrp,current_ik)%msg(i,ipol,im) = &
                              psi(ig+npwx_exx*(ipol-1),im)
                      END DO
                   END DO
                   !
                END DO
                !
                ! send the message
                !
                tag = 0
#if defined(__MPI)
                CALL MPI_ISEND( comm_send_reverse(iproc+1,iegrp,current_ik)%msg, &
                     comm_send_reverse(iproc+1,iegrp,current_ik)%size*npol*my_bands, &
                     MPI_DOUBLE_COMPLEX, &
                     iproc+(iegrp-1)*nproc_egrp, &
                     tag, &
                     intra_pool_comm, request_send(iproc+1,iegrp), ierr )
#endif
                !
             END IF
          END DO
       END DO
    END IF
    !
    ! begin receiving the communication packets
    !
    DO iegrp=1, negrp
       !
       IF ( iexx_istart(iegrp).le.0 ) CYCLE
       !
       recv_bands = iexx_iend(iegrp) - iexx_istart(iegrp) + 1
       !
       DO iproc=0, nproc_egrp-1
          IF ( comm_recv_reverse(iproc+1,current_ik)%size.gt.0) THEN
             !
             !receive the message
             !
             tag = 0
#if defined(__MPI)
             CALL MPI_IRECV( comm_recv_reverse(iproc+1,current_ik)%msg(:,:,iexx_istart(iegrp)), &
                  comm_recv_reverse(iproc+1,current_ik)%size*npol*recv_bands, &
                  MPI_DOUBLE_COMPLEX, &
                  iproc+(iegrp-1)*nproc_egrp, &
                  tag, &
                  intra_pool_comm, request_recv(iproc+1,iegrp), ierr )
#endif
             !
          END IF
       END DO
    END DO
    !
    ! assign psi
    !
    DO iproc=0, nproc_egrp-1
       IF ( comm_recv_reverse(iproc+1,current_ik)%size.gt.0 ) THEN
          !
#if defined(__MPI)
          DO iegrp=1, negrp
             IF ( iexx_istart(iegrp).le.0 ) CYCLE
             CALL MPI_WAIT(request_recv(iproc+1,iegrp), istatus, ierr)
          END DO
#endif
          !
          DO i=1, comm_recv_reverse(iproc+1,current_ik)%size
             ig = comm_recv_reverse(iproc+1,current_ik)%indices(i)
             !
             ! set psi_out
             !
             DO im=1, m
                DO ipol=1, npol
                   psi_out(ig+npwx_local*(ipol-1),im) = psi_out(ig+npwx_local*(ipol-1),im) + &
                        comm_recv_reverse(iproc+1,current_ik)%msg(i,ipol,im)
                END DO
             END DO
             !
          END DO
          !
       END IF
    END DO
    !
    ! wait for everything to finish sending
    !
#if defined(__MPI)
    IF ( iexx_istart(my_egrp_id+1).gt.0 ) THEN
       DO iproc=0, nproc_egrp-1
          DO iegrp=1, negrp
             IF ( comm_send_reverse(iproc+1,iegrp,current_ik)%size.gt.0 ) THEN
                CALL MPI_WAIT(request_send(iproc+1,iegrp), istatus, ierr)
             END IF
          END DO
       END DO
    END IF
#endif
    !
    !-----------------------------------------------------------------------
  END SUBROUTINE transform_to_local
  !-----------------------------------------------------------------------
  !
  !-----------------------------------------------------------------------
  SUBROUTINE exxbuff_comm(exxtemp,ikq,lda,jstart,jend)
  !-----------------------------------------------------------------------
    USE mp_exx,               ONLY : my_egrp_id, inter_egrp_comm, jblock, &
                                     all_end, negrp
    USE mp,             ONLY : mp_bcast
    USE parallel_include
    COMPLEX(DP), intent(inout) :: exxtemp(lda,jend-jstart+1)
    INTEGER, intent(in) :: ikq, lda, jstart, jend
    INTEGER :: jbnd, iegrp, ierr, request_exxbuff(jend-jstart+1)
#if defined(__MPI)
    INTEGER :: istatus(MPI_STATUS_SIZE)
#endif

    DO jbnd=jstart, jend
       !determine which band group has this value of exxbuff
       DO iegrp=1, negrp
          IF(all_end(iegrp) >= jbnd) exit
       END DO

       IF (iegrp == my_egrp_id+1) THEN
          exxtemp(:,jbnd-jstart+1) = exxbuff(:,jbnd,ikq)
       END IF

#if defined(__MPI_NONBLOCKING)
       CALL MPI_IBCAST(exxtemp(:,jbnd-jstart+1), &
            lda, &
            MPI_DOUBLE_COMPLEX, &
            iegrp-1, &
            inter_egrp_comm, &
            request_exxbuff(jbnd-jstart+1), ierr)
#elif defined(__MPI)
       CALL mp_bcast(exxtemp(:,jbnd-jstart+1),iegrp-1,inter_egrp_comm)
#endif
    END DO

    DO jbnd=jstart, jend
#if defined(__MPI_NONBLOCKING)
       CALL MPI_WAIT(request_exxbuff(jbnd-jstart+1), istatus, ierr)
#endif
    END DO

  !
  !-----------------------------------------------------------------------
  END SUBROUTINE exxbuff_comm
  !-----------------------------------------------------------------------
  !
  !-----------------------------------------------------------------------
  SUBROUTINE exxbuff_comm_gamma(exxtemp,ikq,lda,jstart,jend,jlength)
  !-----------------------------------------------------------------------
    USE mp_exx,               ONLY : my_egrp_id, inter_egrp_comm, jblock, &
                                     all_end, negrp
    USE mp,             ONLY : mp_bcast
    USE parallel_include
    COMPLEX(DP), intent(inout) :: exxtemp(lda*npol,jlength)
    INTEGER, intent(in) :: ikq, lda, jstart, jend, jlength
    INTEGER :: jbnd, iegrp, ierr, request_exxbuff(jend-jstart+1), ir
#if defined(__MPI)
    INTEGER :: istatus(MPI_STATUS_SIZE)
#endif
    REAL(DP), ALLOCATABLE :: work(:,:)

    CALL start_clock ('comm_buff')

    ALLOCATE ( work(lda*npol,jend-jstart+1) )
    DO jbnd=jstart, jend
       !determine which band group has this value of exxbuff
       DO iegrp=1, negrp
          IF(all_end(iegrp) >= jbnd) exit
       END DO

       IF (iegrp == my_egrp_id+1) THEN
!          exxtemp(:,jbnd-jstart+1) = exxbuff(:,jbnd,ikq)
          IF( MOD(jbnd,2) == 1 ) THEN
             work(:,jbnd-jstart+1) = REAL(exxbuff(:,jbnd/2+1,ikq))
          ELSE
             work(:,jbnd-jstart+1) = AIMAG(exxbuff(:,jbnd/2,ikq))
          END IF
       END IF

#if defined(__MPI_NONBLOCKING)
       CALL MPI_IBCAST(work(:,jbnd-jstart+1), &
            lda*npol, &
            MPI_DOUBLE_PRECISION, &
            iegrp-1, &
            inter_egrp_comm, &
            request_exxbuff(jbnd-jstart+1), ierr)
#elif defined(__MPI)
       CALL mp_bcast(work(:,jbnd-jstart+1),iegrp-1,inter_egrp_comm)
#endif
    END DO

    DO jbnd=jstart, jend
#if defined(__MPI_NONBLOCKING)
       CALL MPI_WAIT(request_exxbuff(jbnd-jstart+1), istatus, ierr)
#endif
    END DO
    DO jbnd=jstart, jend, 2
       IF (jbnd < jend) THEN
          ! two real bands can be coupled into a single complex one
          DO ir=1, lda*npol
             exxtemp(ir,1+(jbnd-jstart+1)/2) = CMPLX( work(ir,jbnd-jstart+1),&
                                                      work(ir,jbnd-jstart+2),&
                                                      KIND=dp )
          END DO
       ELSE
          ! case of lone last band
          DO ir=1, lda*npol
             exxtemp(ir,1+(jbnd-jstart+1)/2) = CMPLX( work(ir,jbnd-jstart+1),&
                                                      0.0_dp, KIND=dp )
          END DO
       ENDIF
    END DO
    DEALLOCATE ( work )

    CALL stop_clock ('comm_buff')
  !
  !-----------------------------------------------------------------------
  END SUBROUTINE exxbuff_comm_gamma
  !-----------------------------------------------------------------------






!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE matprt(label,n,m,A)
IMPLICIT NONE
  INTEGER :: n,m,i
  real*8 :: A(n,m)
  CHARACTER(len=50) :: frmt
  CHARACTER(len=*) :: label

  WRITE(stdout,'(A)') label
  frmt = ' '
  WRITE(frmt,'(A,I4,A)') '(',m,'f16.10)'
  DO i = 1,n
    WRITE(stdout,frmt) A(i,:)
  ENDDO
END SUBROUTINE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE errinfo(routine,message,INFO)
IMPLICIT NONE
  INTEGER :: INFO
  CHARACTER(len=*) :: routine,message

  IF(INFO/=0) THEN
    WRITE(*,*) routine,' exited with INFO= ',INFO
    CALL errore(routine,message,1)
  ENDIF

END SUBROUTINE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE aceinit( )
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
  USE wavefunctions_module, ONLY : evc
  !
  IMPLICIT NONE
  !
  REAL (DP) :: ee, eexx
  INTEGER :: ik, npw
  TYPE(bec_type) :: becpsi
  !
  nbndproj = nbnd
  IF (.not. allocated(xi)) ALLOCATE( xi(npwx*npol,nbndproj,nks) )
  IF ( okvan ) CALL allocate_bec_type( nkb, nbnd, becpsi)
  eexx = 0.0d0
  xi = (0.0d0,0.0d0)
  DO ik = 1, nks
     npw = ngk (ik)
     current_k = ik
     IF ( lsda ) current_spin = isk(ik)
     IF ( nks > 1 ) CALL get_buffer(evc, nwordwfc, iunwfc, ik)
     IF ( okvan ) THEN
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
  WRITE(stdout,'(/,5X,"ACE energy",f15.8)') eexx
  IF ( okvan ) CALL deallocate_bec_type(becpsi)
  domat = .false.
END SUBROUTINE aceinit
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE aceinit_gamma(nnpw,nbnd,phi,xitmp,becpsi,exxe)
USE becmod,               ONLY : bec_type
USE wvfct,                ONLY : current_k
!
! compute xi(npw,nbndproj) for the ACE method
!
IMPLICIT NONE
  INTEGER :: nnpw,nbnd
  COMPLEX(DP) :: phi(nnpw,nbnd)
  real(DP), ALLOCATABLE :: mexx(:,:)
  COMPLEX(DP) :: xitmp(nnpw,nbndproj)
  INTEGER :: i
  real(DP) :: exxe
  real(DP), PARAMETER :: Zero=0.0d0, One=1.0d0, Two=2.0d0, Pt5=0.50d0
  TYPE(bec_type), INTENT(in) :: becpsi

  CALL start_clock( 'aceinit' )

  IF(nbndproj>nbnd) CALL errore('aceinit','nbndproj greater than nbnd.',1)
  IF(nbndproj<=0) CALL errore('aceinit','nbndproj le 0.',1)

  ALLOCATE( mexx(nbndproj,nbndproj) )
  xitmp = (Zero,Zero)
  mexx = Zero
! |xi> = Vx[phi]|phi>
  CALL vexx(nnpw, nnpw, nbndproj, phi, xitmp, becpsi)
! mexx = <phi|Vx[phi]|phi>
  CALL matcalc('exact',.true.,.false.,nnpw,nbndproj,nbndproj,phi,xitmp,mexx,exxe)
! |xi> = -One * Vx[phi]|phi> * rmexx^T
  CALL aceupdate(nbndproj,nnpw,xitmp,mexx)
  DEALLOCATE( mexx )

  CALL stop_clock( 'aceinit' )

END SUBROUTINE
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
  CALL matcalc('<xi|phi>',.false.,.false.,nnpw,nbndproj,nbnd,xi(1,1,current_k),phi,rmexx,exxe)
! |vv> = |vphi> + (-One) * |xi> * <xi|phi>
  cmexx = (One,Zero)*rmexx
  CALL ZGEMM ('N','N',nnpw,nbnd,nbndproj,-(One,Zero),xi(1,1,current_k), &
                      nnpw,cmexx,nbndproj,(One,Zero),vv,nnpw)
  DEALLOCATE( cmexx,rmexx )

  IF(domat) THEN
    ALLOCATE( rmexx(nbnd,nbnd) )
    CALL matcalc('ACE',.true.,.false.,nnpw,nbnd,nbnd,phi,vv,rmexx,exxe)
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
SUBROUTINE matcalc(label,DoE,PrtMat,ninner,n,m,U,V,mat,ee)
USE becmod,   ONLY : calbec
USE wvfct,    ONLY : current_k, wg
IMPLICIT NONE
!
! compute the (n,n) matrix representation <U|V>
! and energy from V (m,n) and U(m,n)
!
  INTEGER :: ninner,n,m,i
  real(DP) :: ee
  COMPLEX(DP) :: U(ninner,n), V(ninner,m)
  real(DP) :: mat(n,m)
  real(DP), PARAMETER :: Zero=0.0d0, One=1.0d0, Two=2.0d0, Pt5=0.50d0
  CHARACTER(len=*) :: label
  CHARACTER(len=2) :: string
  LOGICAL :: DoE,PrtMat

  CALL start_clock('matcalc')

  string = 'M-'
  mat = Zero
  CALL calbec(ninner, U, V, mat, m)

  IF(DoE) THEN
    IF(n/=m) CALL errore('matcalc','no trace for rectangular matrix.',1)
    IF(PrtMat) CALL matprt(string//label,n,m,mat)
    string = 'E-'
    ee = Zero
    DO i = 1,n
     ee = ee + wg(i,current_k)*mat(i,i)
    ENDDO
    !WRITE(stdout,'(A,f16.8,A)') string//label, ee, ' Ry'
  ENDIF

  CALL stop_clock('matcalc')

END SUBROUTINE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE aceupdate(nbndproj,nnpw,xitmp,rmexx)
IMPLICIT NONE
  INTEGER :: INFO,nbndproj,nnpw
  real(DP) :: rmexx(nbndproj,nbndproj)
  COMPLEX(DP),ALLOCATABLE :: cmexx(:,:)
  COMPLEX(DP) ::  xitmp(nnpw,nbndproj)
  real(DP), PARAMETER :: Zero=0.0d0, One=1.0d0, Two=2.0d0, Pt5=0.50d0

  CALL start_clock('aceupdate')

! rmexx = -(Cholesky(rmexx))^-1
  INFO = -1
  rmexx = -rmexx
  CALL DPOTRF( 'L', nbndproj, rmexx, nbndproj, INFO )
  CALL errinfo('DPOTRF','Cholesky failed in aceupdate.',INFO)
  INFO = -1
  CALL DTRTRI( 'L', 'N', nbndproj, rmexx, nbndproj, INFO )
  CALL errinfo('DTRTRI','inversion failed in aceupdate.',INFO)

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
  CALL matcalc_k('exact',.true.,.false.,current_k,npwx*npol,nbndproj,nbndproj,phi,xitmp,mexx,exxe)
#if defined(__DEBUG)
  WRITE(stdout,'(3(A,I3),A,I9,A,f12.6)') 'aceinit_k: nbnd=', nbnd, ' nbndproj=',nbndproj, &
                                    ' k=',current_k,' npw=',nnpw,' Ex(k)=',exxe
#endif
! |xi> = -One * Vx[phi]|phi> * rmexx^T
  CALL aceupdate_k(nbndproj,nnpw,xitmp,mexx)

  DEALLOCATE( mexx )

  CALL stop_clock( 'aceinit' )

END SUBROUTINE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE matcalc_k(label,DoE,PrtMat,ik,ninner,n,m,U,V,mat,ee)
USE wvfct,                ONLY : wg, npwx
USE becmod,               ONLY : calbec
USE noncollin_module,     ONLY : noncolin, npol
IMPLICIT NONE
!
! compute the (n,n) matrix representation <U|V>
! and energy from V (m,n) and U(m,n)
!
  INTEGER :: ninner,n,m,i,ik, neff
  real(DP) :: ee
  COMPLEX(DP) :: U(ninner,n), V(ninner,m), mat(n,m)
  real(DP), PARAMETER :: Zero=0.0d0, One=1.0d0, Two=2.0d0, Pt5=0.50d0
  CHARACTER(len=*) :: label
  CHARACTER(len=2) :: string
  LOGICAL :: DoE,PrtMat

  CALL start_clock('matcalc')

  string = 'M-'
  mat = (Zero,Zero)
  IF(noncolin) THEN
    noncolin = .false.
    CALL calbec(ninner, U, V, mat, m)
    noncolin = .true.
  ELSE
    CALL calbec(ninner, U, V, mat, m)
  ENDIF

  IF(DoE) THEN
    IF(n/=m) CALL errore('matcalc','no trace for rectangular matrix.',1)
    IF(PrtMat) CALL matprt_k(string//label,n,m,mat)
    string = 'E-'
    ee = Zero
    DO i = 1,n
      ee = ee + wg(i,ik)*DBLE(mat(i,i))
    ENDDO
    !write(stdout,'(A,f16.8,A)') string//label, ee, ' Ry'
  ENDIF

  CALL stop_clock('matcalc')

END SUBROUTINE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE matprt_k(label,n,m,A)
IMPLICIT NONE
  INTEGER :: n,m,i
  COMPLEX(DP) :: A(n,m)
  CHARACTER(len=50) :: frmt
  CHARACTER(len=*) :: label

  WRITE(stdout,'(A)') label//'(real)'
  frmt = ' '
  WRITE(frmt,'(A,I4,A)') '(',m,'f12.6)'
  DO i = 1,n
    WRITE(stdout,frmt) dreal(A(i,:))
  ENDDO

  WRITE(stdout,'(A)') label//'(imag)'
  frmt = ' '
  WRITE(frmt,'(A,I4,A)') '(',m,'f12.6)'
  DO i = 1,n
    WRITE(stdout,frmt) aimag(A(i,:))
  ENDDO
END SUBROUTINE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE aceupdate_k(nbndproj,nnpw,xitmp,mexx)
USE wvfct ,               ONLY : npwx
USE noncollin_module,     ONLY : noncolin, npol
IMPLICIT NONE
  INTEGER :: INFO,nbndproj,nnpw
  COMPLEX(DP) :: mexx(nbndproj,nbndproj), xitmp(npwx*npol,nbndproj)
  real(DP), PARAMETER :: Zero=0.0d0, One=1.0d0, Two=2.0d0, Pt5=0.50d0

  CALL start_clock('aceupdate')

! mexx = -(Cholesky(mexx))^-1
  INFO = -1
  mexx = -mexx
  CALL ZPOTRF( 'L', nbndproj, mexx, nbndproj, INFO )
  CALL errinfo('DPOTRF','Cholesky failed in aceupdate.',INFO)
  INFO = -1
  CALL ZTRTRI( 'L', 'N', nbndproj, mexx, nbndproj, INFO )
  CALL errinfo('DTRTRI','inversion failed in aceupdate.',INFO)
! |xi> = -One * Vx[phi]|phi> * mexx^T
  CALL ZTRMM('R','L','C','N',npwx*npol,nbndproj,(One,Zero),mexx,nbndproj,xitmp,npwx*npol)

  CALL stop_clock('aceupdate')

END SUBROUTINE



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
  CALL matcalc_k('<xi|phi>',.false.,.false.,current_k,npwx*npol,nbndproj,nbnd,xi(1,1,current_k),phi,cmexx,exxe)

! |vv> = |vphi> + (-One) * |xi> * <xi|phi>
  CALL ZGEMM ('N','N',npwx*npol,nbnd,nbndproj,-(One,Zero),xi(1,1,current_k),npwx*npol,cmexx,nbndproj,(One,Zero),vv,npwx*npol)

  IF(domat) THEN
     CALL matcalc_k('ACE',.true.,.false.,current_k,npwx*npol,nbnd,nbnd,phi,vv,cmexx,exxe)
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

END SUBROUTINE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-----------------------------------------------------------------------
END MODULE exx
!-----------------------------------------------------------------------
