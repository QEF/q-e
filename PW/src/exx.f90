
! Copyright (C) 2005-2015 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------
MODULE exx
  !--------------------------------------
!   !
  USE kinds,                ONLY : DP
  USE coulomb_vcut_module,  ONLY : vcut_init, vcut_type, vcut_info, &
                                   vcut_get,  vcut_spheric_get
  USE noncollin_module,     ONLY : noncolin, npol
  USE io_global,            ONLY : ionode
  USE fft_custom,           ONLY : fft_cus
  !
  USE control_flags,        ONLY : gamma_only, tqr

  IMPLICIT NONE
  SAVE
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
!civn
  COMPLEX(DP), ALLOCATABLE :: xi(:,:,:)
  INTEGER :: nbndproj
  LOGICAL :: domat
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
  LOGICAL           :: on_double_grid =.false.
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
    USE mp_bands,     ONLY : intra_bgrp_comm
    USE control_flags,ONLY : tqr
    USE realus,       ONLY : qpointlist, tabxx, tabp!, tabs
    USE io_global,    ONLY : stdout

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
    CALL data_structure_custom(exx_fft, gamma_only)
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
    !USE exx,            ONLY : nkqs, working_pool
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
    USE io_global,  ONLY : stdout
    USE start_k,    ONLY : nk1,nk2,nk3
    USE mp_pools,   ONLY : npool
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
    USE io_global,  ONLY : stdout
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
      !
    CASE ( "vcut_ws" )
      !
      use_coulomb_vcut_ws = .true.
      IF ( x_gamma_extrapolation ) &
            CALL errore(sub_name,'cannot USE x_gamm_extrap and vcut_ws', 1)
      !
    CASE ( "vcut_spherical" )
      !
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
    USE mp_pools,   ONLY : npool
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
  SUBROUTINE exx_restart(l_exx_was_active)
     !------------------------------------------------------------------------
     !This SUBROUTINE is called when restarting an exx calculation
     USE funct,                ONLY : get_exx_fraction, start_exx, &
                                      exx_is_active, get_screening_parameter

     IMPLICIT NONE
     LOGICAL, INTENT(in) :: l_exx_was_active

     IF (.not. l_exx_was_active ) RETURN ! nothing had happened yet
     !
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
    USE io_files,             ONLY : nwordwfc, iunwfc
    USE buffers,              ONLY : get_buffer
    USE wvfct,                ONLY : nbnd, npwx, wg, current_k
    USE klist,                ONLY : ngk, nks, nkstot, xk, wk, igk_k
    USE symm_base,            ONLY : nsym, s, sr, ftau
    USE mp_pools,             ONLY : npool, nproc_pool, me_pool, inter_pool_comm
    USE mp_bands,             ONLY : me_bgrp, set_bgrp_indices, nbgrp
    USE mp,                   ONLY : mp_sum, mp_bcast
    USE funct,                ONLY : get_exx_fraction, start_exx,exx_is_active,&
                                     get_screening_parameter, get_gau_parameter
    USE scatter_mod,          ONLY : gather_grid, scatter_grid
    USE fft_interfaces,       ONLY : invfft
    USE uspp,                 ONLY : nkb, vkb, okvan
    USE us_exx,               ONLY : rotate_becxx
    USE paw_variables,        ONLY : okpaw
    USE paw_exx,              ONLY : PAW_init_keeq
    USE mp_pools,             ONLY : me_pool, my_pool_id, root_pool, nproc_pool, &
                                     inter_pool_comm, my_pool_id, intra_pool_comm
    USE mp_orthopools,        ONLY : intra_orthopool_comm

    IMPLICIT NONE
    INTEGER :: ik,ibnd, i, j, k, ir, isym, ikq, ig
    INTEGER :: h_ibnd
    INTEGER :: ibnd_loop_start, ibnd_buff_start, ibnd_buff_end
    INTEGER :: ipol, jpol
    REAL(dp), ALLOCATABLE   :: occ(:,:)
    COMPLEX(DP),ALLOCATABLE :: temppsic(:)
    COMPLEX(DP),ALLOCATABLE :: temppsic_nc(:,:), psic_nc(:,:)
    INTEGER :: nxxs, nrxxs
#if defined(__MPI)
    COMPLEX(DP),ALLOCATABLE  :: temppsic_all(:),      psic_all(:)
    COMPLEX(DP), ALLOCATABLE :: temppsic_all_nc(:,:), psic_all_nc(:,:)
#endif
    COMPLEX(DP) :: d_spin(2,2,48)
    INTEGER :: npw, current_ik
    INTEGER,EXTERNAL :: global_kpoint_index
    !
    CALL start_clock ('exxinit')
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

    CALL set_bgrp_indices(x_nbnd_occ,ibnd_start,ibnd_end)

    IF ( gamma_only ) THEN
        ibnd_buff_start = ibnd_start/2
        IF(mod(ibnd_start,2)==1) ibnd_buff_start = ibnd_buff_start +1
        !
        ibnd_buff_end = ibnd_end/2
        IF(mod(ibnd_end,2)==1) ibnd_buff_end = ibnd_buff_end +1
    ELSE
        ibnd_buff_start = ibnd_start
        ibnd_buff_end   = ibnd_end
    ENDIF
    !
    IF (.not. allocated(exxbuff)) &
        ALLOCATE( exxbuff(nrxxs*npol, ibnd_buff_start:ibnd_buff_end, nkqs))
    exxbuff=(0.0_DP,0.0_DP)
    !
    !   This is parallelized over pools. Each pool computes only its k-points
    !
    KPOINTS_LOOP : &
    DO ik = 1, nks

       IF ( nks > 1 ) &
          CALL get_buffer(evc, nwordwfc, iunwfc, ik)
       !
       ! ik         = index of k-point in this pool
       ! current_ik = index of k-point over all pools
       !
       current_ik = global_kpoint_index ( nkstot, ik )
       !
       IF_GAMMA_ONLY : &
       IF (gamma_only) THEN
          !
          h_ibnd = ibnd_start/2
          !
          IF(mod(ibnd_start,2)==0) THEN
             h_ibnd=h_ibnd-1
             ibnd_loop_start=ibnd_start-1
          ELSE
             ibnd_loop_start=ibnd_start
          ENDIF

          DO ibnd = ibnd_loop_start, ibnd_end, 2
             h_ibnd = h_ibnd + 1
             !
             psic(:) = ( 0._dp, 0._dp )
             !
             IF ( ibnd < ibnd_end ) THEN
                DO ig=1,exx_fft%npwt
                   psic(exx_fft%nlt(ig))  = evc(ig,ibnd)  &
                        + ( 0._dp, 1._dp ) * evc(ig,ibnd+1)
                   psic(exx_fft%nltm(ig)) = conjg( evc(ig,ibnd) ) &
                        + ( 0._dp, 1._dp ) * conjg( evc(ig,ibnd+1) )
                ENDDO
             ELSE
                DO ig=1,exx_fft%npwt
                   psic(exx_fft%nlt (ig)) = evc(ig,ibnd)
                   psic(exx_fft%nltm(ig)) = conjg( evc(ig,ibnd) )
                ENDDO
             ENDIF

             CALL invfft ('CustomWave', psic, exx_fft%dfftt)

             exxbuff(1:nrxxs,h_ibnd,ik)=psic(1:nrxxs)

          ENDDO
          !
       ELSE IF_GAMMA_ONLY
          !
          npw = ngk (ik)
          IBND_LOOP_K : &
          DO ibnd = ibnd_start, ibnd_end
             !
             IF (noncolin) THEN
                temppsic_nc(:,:) = ( 0._dp, 0._dp )
                temppsic_nc(exx_fft%nlt(igk_k(1:npw,ik)),1) = evc(1:npw,ibnd)
                CALL invfft ('CustomWave', temppsic_nc(:,1), exx_fft%dfftt)
                temppsic_nc(exx_fft%nlt(igk_k(1:npw,ik)),2) = evc(npwx+1:npwx+npw,ibnd)
                CALL invfft ('CustomWave', temppsic_nc(:,2), exx_fft%dfftt)
             ELSE
                temppsic(:) = ( 0._dp, 0._dp )
                temppsic(exx_fft%nlt(igk_k(1:npw,ik))) = evc(1:npw,ibnd)
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
                   IF ( me_bgrp == 0 ) THEN
                      psic_all_nc(:,:) = (0.0_DP, 0.0_DP)
                      DO ipol=1,npol
                         DO jpol=1,npol
                            psic_all_nc(:,ipol)=psic_all_nc(:,ipol) &
                              +  conjg(d_spin(jpol,ipol,isym))* &
                                 temppsic_all_nc(rir(:,isym),jpol)
                         ENDDO
                      ENDDO
                   ENDIF
                   DO ipol=1,npol
                      CALL scatter_grid(exx_fft%dfftt,psic_all_nc(:,ipol), psic_nc(:,ipol))
                   ENDDO
#else
                   psic_nc(:,:) = (0._dp, 0._dp)
                   DO ipol=1,npol
                      DO jpol=1,npol
                         psic_nc(:,ipol) = psic_nc(:,ipol) + &
                              conjg(d_spin(jpol,ipol,isym))* &
                                        temppsic_nc(rir(:,isym),jpol)
                      ENDDO
                   ENDDO
#endif
                   exxbuff(      1:  nrxxs,ibnd,ikq)=psic_nc(:,1)
                   exxbuff(nrxxs+1:2*nrxxs,ibnd,ikq)=psic_nc(:,2)
                ELSE ! noncolinear
#if defined(__MPI)
                  CALL gather_grid(exx_fft%dfftt,temppsic,temppsic_all)
                  IF ( me_bgrp == 0 ) &
                    psic_all(1:nxxs) = temppsic_all(rir(1:nxxs,isym))
                  CALL scatter_grid(exx_fft%dfftt,psic_all,psic)
#else
                  psic(1:nrxxs) = temppsic(rir(1:nrxxs,isym))
#endif
                  IF (index_sym(ikq) < 0 ) psic(1:nrxxs) = conjg(psic(1:nrxxs))
                  exxbuff(1:nrxxs,ibnd,ikq)=psic(1:nrxxs)
                  !
                ENDIF ! noncolinear
             ENDDO
             !
          ENDDO &
          IBND_LOOP_K
          !
       ENDIF &
       IF_GAMMA_ONLY
    ENDDO &
    KPOINTS_LOOP
    !
    IF (noncolin) THEN
       DEALLOCATE(temppsic_nc, psic_nc)
#if defined(__MPI)
       DEALLOCATE(temppsic_all_nc, psic_all_nc)
#endif
    ELSEIF ( .not. gamma_only ) THEN
       DEALLOCATE(temppsic)
#if defined(__MPI)
       DEALLOCATE(temppsic_all, psic_all)
#endif
    ENDIF
    !
    ! Each wavefunction in exxbuff is computed by a single pool, collect among 
    ! pools in a smart way (i.e. without doing all-to-all sum and bcast)
    ! See also the initialization of working_pool in exx_mp_init
!      IF (npool>1 ) CALL mp_sum(exxbuff, inter_pool_comm)
    DO ikq = 1, nkqs
      CALL mp_bcast(exxbuff(:,:,ikq), working_pool(ikq), intra_orthopool_comm)
    ENDDO
    !
    ! For US/PAW only: compute <beta_I|psi_j,k+q> for the entire 
    ! de-symmetrized k+q grid by rotating the ones fro mthe irreducible wedge
    IF(okvan) CALL rotate_becxx(nkqs, index_xk, index_sym, xkq_collect)
    !
    ! Initialize 4-wavefunctions one-center Fock integrals
    !    \int \psi_a(r)\phi_a(r)\phi_b(r')\psi_b(r')/|r-r'|
    IF(okpaw) CALL PAW_init_keeq()
    !
#if defined(__EXX_ACE)
    CALL aceinit ( )
#endif
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
    USE symm_base,            ONLY : nsym, s, sr, ftau
    !
    IMPLICIT NONE
    !
    INTEGER :: nxxs, nr1,nr2,nr3, nr1x,nr2x,nr3x
    INTEGER :: ikq, isym, i,j,k, ri,rj,rk, ir
    LOGICAL :: ispresent(nsym)
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
             CALL errore ('exxinit',' EXX smooth grid is not compatible with &
                                    & symmetry: change ecutfock',isym)
          ENDIF
          DO ir=1, nxxs
             rir(ir,isym) = ir
          ENDDO
          DO k = 1, nr3
             DO j = 1, nr2
                DO i = 1, nr1
                   CALL ruotaijk (s(1,1,isym), ftau(1,isym), i,j,k, nr1,nr2,nr3, ri,rj,rk)
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
    !
    IMPLICIT NONE
    !
    INTEGER                  :: lda, n, m
    COMPLEX(DP)              :: psi(lda*npol,m)
    COMPLEX(DP)              :: hpsi(lda*npol,m)
    TYPE(bec_type), OPTIONAL :: becpsi
    !
    IF ( (okvan.or.okpaw) .and. .not. present(becpsi)) &
       CALL errore('vexx','becpsi needed for US/PAW case',1)
    CALL start_clock ('vexx')
    !
    IF(gamma_only) THEN
       CALL vexx_gamma(lda, n, m, psi, hpsi, becpsi)
    ELSE
       CALL vexx_k(lda, n, m, psi, hpsi, becpsi)
    ENDIF
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
    USE wvfct,          ONLY : npwx, current_k
    USE klist,          ONLY : xk, nks, nkstot, igk_k
    USE fft_interfaces, ONLY : fwfft, invfft
    USE becmod,         ONLY : bec_type
    USE mp_bands,       ONLY : inter_bgrp_comm, intra_bgrp_comm, my_bgrp_id, nbgrp
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
    COMPLEX(DP),ALLOCATABLE :: RESULT(:)
    REAL(DP),ALLOCATABLE :: temppsic_dble (:)
    REAL(DP),ALLOCATABLE :: temppsic_aimag(:)
    !
    COMPLEX(DP),ALLOCATABLE :: rhoc(:), vc(:), deexx(:)
    REAL(DP),   ALLOCATABLE :: fac(:)
    INTEGER          :: ibnd, ik, im , ikq, iq, ipol
    INTEGER          :: ir, ig
    INTEGER          :: current_ik
    INTEGER          :: ibnd_loop_start
    INTEGER          :: h_ibnd, nrxxs
    REAL(DP) :: x1, x2, xkp(3)
    REAL(DP) :: xkq(3)
    ! <LMS> temp array for vcut_spheric
    INTEGER, EXTERNAL  :: global_kpoint_index
    LOGICAL :: l_fft_doubleband
    LOGICAL :: l_fft_singleband
    !
    ALLOCATE( fac(exx_fft%ngmt) )
    nrxxs= exx_fft%dfftt%nnr
    !
    ALLOCATE( RESULT(nrxxs), temppsic_dble(nrxxs), temppsic_aimag(nrxxs) )
    !
    ALLOCATE(rhoc(nrxxs), vc(nrxxs))
    IF(okvan) ALLOCATE(deexx(nkb))
    !
    current_ik = global_kpoint_index ( nkstot, current_k )
    xkp = xk(:,current_k)
    !
    ! This is to stop numerical inconsistencies creeping in through the band parallelization.
    !
    IF(my_bgrp_id>0) THEN
       hpsi=0.0_DP
       psi =0.0_DP
    ENDIF
    IF (nbgrp>1) THEN
       CALL mp_bcast(hpsi,0,inter_bgrp_comm)
       CALL mp_bcast(psi,0,inter_bgrp_comm)
    ENDIF
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
       CALL g2_convolution(exx_fft%ngmt, exx_fft%gt, xkp, xkq, fac)
       IF ( okvan .and..not.tqr ) CALL qvan_init (exx_fft%ngmt, xkq, xkp)
       !
       LOOP_ON_PSI_BANDS : &
       DO im = 1,m !for each band of psi (the k cycle is outside band)
          IF(okvan) deexx(:) = 0.0_DP
          !
          RESULT = 0.0_DP
          !
          l_fft_doubleband = .false.
          l_fft_singleband = .false.
          !
          IF ( mod(im,2)==1 .and. (im+1)<=m ) l_fft_doubleband = .true.
          IF ( mod(im,2)==1 .and. im==m )     l_fft_singleband = .true.
          !
          IF( l_fft_doubleband ) THEN
!$omp parallel do  default(shared), private(ig)
             DO ig = 1, exx_fft%npwt
                RESULT( exx_fft%nlt(ig) )  =       psi(ig, im) + (0._DP,1._DP) * psi(ig, im+1)
                RESULT( exx_fft%nltm(ig) ) = conjg(psi(ig, im) - (0._DP,1._DP) * psi(ig, im+1))
             ENDDO
!$omp end parallel do
          ENDIF
          !
          IF( l_fft_singleband ) THEN
!$omp parallel do  default(shared), private(ig)
             DO ig = 1, exx_fft%npwt
                RESULT( exx_fft%nlt(ig) )  =       psi(ig,im)
                RESULT( exx_fft%nltm(ig) ) = conjg(psi(ig,im))
             ENDDO
!$omp end parallel do
          ENDIF
          !
          IF( l_fft_doubleband.or.l_fft_singleband) THEN
             CALL invfft ('CustomWave', RESULT, exx_fft%dfftt)
!$omp parallel do default(shared), private(ir)
             DO ir = 1, nrxxs
                temppsic_dble(ir)  = dble ( RESULT(ir) )
                temppsic_aimag(ir) = aimag( RESULT(ir) )
             ENDDO
!$omp end parallel do
          ENDIF
          !
          RESULT = 0.0_DP
          !
          h_ibnd = ibnd_start/2
          IF(mod(ibnd_start,2)==0) THEN
             h_ibnd=h_ibnd-1
             ibnd_loop_start=ibnd_start-1
          ELSE
             ibnd_loop_start=ibnd_start
          ENDIF
          !
          IBND_LOOP_GAM : &
          DO ibnd=ibnd_loop_start,ibnd_end, 2 !for each band of psi
             !
             h_ibnd = h_ibnd + 1
             IF( ibnd < ibnd_start ) THEN
                x1 = 0.0_DP
             ELSE
                x1 = x_occupation(ibnd,  ik)
             ENDIF
             IF( ibnd == ibnd_end) THEN
                x2 = 0.0_DP
             ELSE
                x2 = x_occupation(ibnd+1,  ik)
             ENDIF
             IF ( abs(x1) < eps_occ .and. abs(x2) < eps_occ ) CYCLE
             !
             ! calculate rho in real space. Gamma tricks are used.
             ! temppsic is real; tempphic contains one band in the real part,
             ! another one in the imaginary part; the same applies to rhoc
             !
             IF( mod(im,2) == 0 ) THEN
!$omp parallel do default(shared), private(ir)
                DO ir = 1, nrxxs
                   rhoc(ir) = exxbuff(ir,h_ibnd,ikq) * temppsic_aimag(ir) / omega
                ENDDO
!$omp end parallel do
             ELSE
!$omp parallel do default(shared), private(ir)
                DO ir = 1, nrxxs
                   rhoc(ir) = exxbuff(ir,h_ibnd,ikq) * temppsic_dble(ir) / omega
                ENDDO
!$omp end parallel do
             ENDIF
             !
             ! bring rho to G-space
             !
             !   >>>> add augmentation in REAL SPACE here
             IF(okvan .and. tqr) THEN
                IF(ibnd>=ibnd_start) &
                CALL addusxx_r(exx_fft,rhoc,_CX(becxx(ikq)%r(:,ibnd)), _CX(becpsi%r(:,im)))
                IF(ibnd<ibnd_end) &
                CALL addusxx_r(exx_fft,rhoc,_CY(becxx(ikq)%r(:,ibnd+1)),_CX(becpsi%r(:,im)))
             ENDIF
             !
             CALL fwfft ('Custom', rhoc, exx_fft%dfftt)
             !   >>>> add augmentation in G SPACE here
             IF(okvan .and. .not. tqr) THEN
                ! contribution from one band added to real (in real space) part of rhoc
                IF(ibnd>=ibnd_start) &
                   CALL addusxx_g(exx_fft, rhoc, xkq,  xkp, 'r', &
                   becphi_r=becxx(ikq)%r(:,ibnd), becpsi_r=becpsi%r(:,im) )
                ! contribution from following band added to imaginary (in real space) part of rhoc
                IF(ibnd<ibnd_end) &
                   CALL addusxx_g(exx_fft, rhoc, xkq,  xkp, 'i', &
                   becphi_r=becxx(ikq)%r(:,ibnd+1), becpsi_r=becpsi%r(:,im) )
             ENDIF
             !   >>>> charge density done
             !
             vc = 0._DP
             !
!$omp parallel do default(shared), private(ig)
             DO ig = 1, exx_fft%ngmt
                !
                vc(exx_fft%nlt(ig))  = fac(ig) * rhoc(exx_fft%nlt(ig))
                vc(exx_fft%nltm(ig)) = fac(ig) * rhoc(exx_fft%nltm(ig))
                !
             ENDDO
!$omp end parallel do
             !
             !   >>>>  compute <psi|H_fock G SPACE here
             IF(okvan .and. .not. tqr) THEN
                IF(ibnd>=ibnd_start) &
                CALL newdxx_g(exx_fft, vc, xkq, xkp, 'r', deexx, &
                              becphi_r=x1*becxx(ikq)%r(:,ibnd))
                IF(ibnd<ibnd_end) &
                CALL newdxx_g(exx_fft, vc, xkq, xkp, 'i', deexx, &
                              becphi_r=x2*becxx(ikq)%r(:,ibnd+1))
             ENDIF
             !
             !brings back v in real space
             CALL invfft ('Custom', vc, exx_fft%dfftt)
             !
             !   >>>>  compute <psi|H_fock REAL SPACE here
             IF(okvan .and. tqr) THEN
                IF(ibnd>=ibnd_start) &
                CALL newdxx_r(exx_fft,vc, _CX(x1*becxx(ikq)%r(:,ibnd)), deexx)
                IF(ibnd<ibnd_end) &
                CALL newdxx_r(exx_fft,vc, _CY(x2*becxx(ikq)%r(:,ibnd+1)), deexx)
             ENDIF
             !
             IF(okpaw) THEN
                IF(ibnd>=ibnd_start) &
                CALL PAW_newdxx(x1/nqs, _CX(becxx(ikq)%r(:,ibnd)),   _CX(becpsi%r(:,im)), deexx)
                IF(ibnd<ibnd_end) &
                CALL PAW_newdxx(x2/nqs, _CX(becxx(ikq)%r(:,ibnd+1)), _CX(becpsi%r(:,im)), deexx)
             ENDIF
             !
             ! accumulates over bands and k points
             !
!$omp parallel do default(shared), private(ir)
             DO ir = 1, nrxxs
                RESULT(ir) = RESULT(ir)+x1* dble(vc(ir))* dble(exxbuff(ir,h_ibnd,ikq))&
                                       +x2*aimag(vc(ir))*aimag(exxbuff(ir,h_ibnd,ikq))
             ENDDO
!$omp end parallel do
             !
          ENDDO &
          IBND_LOOP_GAM
          !
          !
          IF(okvan) THEN
             CALL mp_sum(deexx,intra_bgrp_comm)
             CALL mp_sum(deexx,inter_bgrp_comm)
          ENDIF
          !
          CALL mp_sum( RESULT(1:nrxxs), inter_bgrp_comm)
          !
          ! brings back result in G-space
          !
          CALL fwfft( 'CustomWave' , RESULT, exx_fft%dfftt )
          !
!$omp parallel do default(shared), private(ig)
          DO ig = 1, n
             hpsi(ig,im)=hpsi(ig,im) - exxalfa*RESULT(exx_fft%nlt(ig))
          ENDDO
!$omp end parallel do
          ! add non-local \sum_I |beta_I> \alpha_Ii (the sum on i is outside)
          IF(okvan) CALL add_nlxx_pot (lda, hpsi(:,im), xkp, n, &
                           igk_k(1,current_k), deexx, eps_occ, exxalfa)
       ENDDO &
       LOOP_ON_PSI_BANDS
       IF ( okvan .and..not.tqr ) CALL qvan_clean ()
       !
    ENDDO &
    INTERNAL_LOOP_ON_Q
    !
    DEALLOCATE( RESULT, temppsic_dble, temppsic_aimag)
    !
    DEALLOCATE(rhoc, vc, fac )
    !
    IF(okvan) DEALLOCATE( deexx )
    !
    !-----------------------------------------------------------------------
  END SUBROUTINE vexx_gamma
  !-----------------------------------------------------------------------
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
    USE wvfct,          ONLY : npwx, current_k
    USE klist,          ONLY : xk, nks, nkstot, igk_k
    USE fft_interfaces, ONLY : fwfft, invfft
    USE becmod,         ONLY : bec_type
    USE mp_bands,       ONLY : inter_bgrp_comm, intra_bgrp_comm, my_bgrp_id, nbgrp
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
    COMPLEX(DP),ALLOCATABLE :: temppsic(:), RESULT(:)
    COMPLEX(DP),ALLOCATABLE :: temppsic_nc(:,:),result_nc(:,:)
    COMPLEX(DP),ALLOCATABLE :: result_g(:), result_nc_g(:,:)
    !
    COMPLEX(DP),ALLOCATABLE :: rhoc(:), vc(:), deexx(:)
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
    !
    ALLOCATE( fac(exx_fft%ngmt) )
    nrxxs= exx_fft%dfftt%nnr
    !
    IF (noncolin) THEN
       ALLOCATE( temppsic_nc(nrxxs,npol), result_nc(nrxxs,npol) )
       ALLOCATE( result_nc_g(n,npol) )
    ELSE
       ALLOCATE( temppsic(nrxxs), RESULT(nrxxs) )
       ALLOCATE( result_g(n) )
    ENDIF
    !
    ALLOCATE(rhoc(nrxxs), vc(nrxxs))
    IF(okvan) ALLOCATE(deexx(nkb))
    !
    current_ik = global_kpoint_index ( nkstot, current_k )
    xkp = xk(:,current_k)
    !
    ! This is to stop numerical inconsistencies creeping in through the band parallelization.
    !
    IF(my_bgrp_id>0) THEN
       hpsi=0.0_DP
       psi =0.0_DP
    ENDIF
    IF (nbgrp>1) THEN
       CALL mp_bcast(hpsi,0,inter_bgrp_comm)
       CALL mp_bcast(psi,0,inter_bgrp_comm)
    ENDIF
    !
    LOOP_ON_PSI_BANDS : &
    DO im = 1,m !for each band of psi (the k cycle is outside band)
       IF(okvan) deexx = 0.0_DP
       !
       IF (noncolin) THEN
          temppsic_nc = 0._DP
       ELSE
          temppsic    = 0.0_DP
       ENDIF
       !
       IF (noncolin) THEN
          !
!$omp parallel do  default(shared), private(ig)
          DO ig = 1, n
             temppsic_nc(exx_fft%nlt(igk_k(ig,current_k)),1) = psi(ig,im)
          ENDDO
!$omp end parallel do
!$omp parallel do  default(shared), private(ig)
          DO ig = 1, n
             temppsic_nc(exx_fft%nlt(igk_k(ig,current_k)),2) = psi(npwx+ig,im)
          ENDDO
!$omp end parallel do
          !
          CALL invfft ('CustomWave', temppsic_nc(:,1), exx_fft%dfftt)
          CALL invfft ('CustomWave', temppsic_nc(:,2), exx_fft%dfftt)
          !
       ELSE
          !
!$omp parallel do  default(shared), private(ig)
          DO ig = 1, n
             temppsic( exx_fft%nlt(igk_k(ig,current_k)) ) = psi(ig,im)
          ENDDO
!$omp end parallel do
          CALL invfft ('CustomWave', temppsic, exx_fft%dfftt)
          !
       ENDIF
       !
       IF (noncolin) THEN
          result_nc = 0.0_DP
       ELSE
          RESULT    = 0.0_DP
       ENDIF
       !
       INTERNAL_LOOP_ON_Q : &
       DO iq=1,nqs
          !
          ikq  = index_xkq(current_ik,iq)
          ik   = index_xk(ikq)
          xkq  = xkq_collect(:,ikq)
          !
          ! calculate the 1/|r-r'| (actually, k+q+g) factor and place it in fac
          CALL g2_convolution(exx_fft%ngmt, exx_fft%gt, xkp, xkq, fac)
          IF ( okvan .and..not.tqr ) CALL qvan_init (exx_fft%ngmt, xkq, xkp)
          !
          IBND_LOOP_K : &
          DO ibnd=ibnd_start,ibnd_end !for each band of psi
             !
             IF ( abs(x_occupation(ibnd,ik)) < eps_occ) CYCLE IBND_LOOP_K
             !
             !loads the phi from file
             !
             !   >>>> calculate rho in real space
             IF (noncolin) THEN
!$omp parallel do default(shared), private(ir)
                DO ir = 1, nrxxs
                   rhoc(ir) = ( conjg(exxbuff(ir,ibnd,ikq))*temppsic_nc(ir,1) + &
                                conjg(exxbuff(nrxxs+ir,ibnd,ikq))*temppsic_nc(ir,2) )/omega
                ENDDO
!$omp end parallel do
             ELSE
!$omp parallel do default(shared), private(ir)
                DO ir = 1, nrxxs
                   rhoc(ir)=conjg(exxbuff(ir,ibnd,ikq))*temppsic(ir) / omega
                ENDDO
!$omp end parallel do
             ENDIF
             !   >>>> add augmentation in REAL space HERE
             IF(okvan .and. tqr) THEN ! augment the "charge" in real space
                CALL addusxx_r(exx_fft,rhoc, becxx(ikq)%k(:,ibnd), becpsi%k(:,im))
             ENDIF
             !
             !   >>>> brings it to G-space
             CALL fwfft('Custom', rhoc, exx_fft%dfftt)
             !
             !   >>>> add augmentation in G space HERE
             IF(okvan .and. .not. tqr) THEN
                CALL addusxx_g(exx_fft, rhoc, xkq, xkp, 'c', &
                   becphi_c=becxx(ikq)%k(:,ibnd),becpsi_c=becpsi%k(:,im))
             ENDIF
             !   >>>> charge done
             !
             vc = 0._DP
             !
!$omp parallel do default(shared), private(ig)
             DO ig = 1, exx_fft%ngmt
                vc(exx_fft%nlt(ig)) = fac(ig) * rhoc(exx_fft%nlt(ig)) * &
                                             x_occupation(ibnd,ik) / nqs
             ENDDO
!$omp end parallel do
             !
             ! Add ultrasoft contribution (RECIPROCAL SPACE)
             ! compute alpha_I,j,k+q = \sum_J \int <beta_J|phi_j,k+q> V_i,j,k,q Q_I,J(r) d3r
             IF(okvan .and. .not. tqr) THEN
                CALL newdxx_g(exx_fft, vc, xkq, xkp, 'c', deexx, &
                              becphi_c=becxx(ikq)%k(:,ibnd))
             ENDIF
             !
             !brings back v in real space
             CALL invfft ('Custom', vc, exx_fft%dfftt)
             !
             ! Add ultrasoft contribution (REAL SPACE)
             IF(okvan .and. tqr) CALL newdxx_r(exx_fft,vc, becxx(ikq)%k(:,ibnd),deexx)
             !
             ! Add PAW one-center contribution
             IF(okpaw) THEN
                CALL PAW_newdxx(x_occupation(ibnd,ik)/nqs, becxx(ikq)%k(:,ibnd), becpsi%k(:,im), deexx)
             ENDIF
             !
             !accumulates over bands and k points
             !
             IF (noncolin) THEN
!$omp parallel do default(shared), private(ir)
                DO ir = 1, nrxxs
                   result_nc(ir,1)= result_nc(ir,1) + vc(ir) * exxbuff(ir,ibnd,ikq)
                ENDDO
!$omp end parallel do
!$omp parallel do default(shared), private(ir)
                DO ir = 1, nrxxs
                   result_nc(ir,2)= result_nc(ir,2) + vc(ir) * exxbuff(ir+nrxxs,ibnd,ikq)
                ENDDO
!$omp end parallel do
             ELSE
!$omp parallel do default(shared), private(ir)
                DO ir = 1, nrxxs
                   RESULT(ir) = RESULT(ir) + vc(ir)*exxbuff(ir,ibnd,ikq)
                ENDDO
!$omp end parallel do
             ENDIF
             !
          ENDDO &
          IBND_LOOP_K
          IF ( okvan .and..not.tqr ) CALL qvan_clean ()
          !
       ENDDO &
       INTERNAL_LOOP_ON_Q
       !
       IF(okvan) THEN
         CALL mp_sum(deexx,intra_bgrp_comm)
         CALL mp_sum(deexx,inter_bgrp_comm)
       ENDIF
       !
       ! bring result back to G-space
       !
       IF (noncolin) THEN
          !
          CALL fwfft ('CustomWave', result_nc(:,1), exx_fft%dfftt)
          CALL fwfft ('CustomWave', result_nc(:,2), exx_fft%dfftt)
          !
          !communicate result
          DO ig = 1, n
             result_nc_g(ig,1:npol) = result_nc(exx_fft%nlt(igk_k(ig,current_k)),1:npol)
          ENDDO
          CALL mp_sum( result_nc_g(1:n,1:npol), inter_bgrp_comm)
          !
          !adds it to hpsi
!$omp parallel do default(shared), private(ig)
          DO ig = 1, n
             hpsi(ig,im)    = hpsi(ig,im)     - exxalfa*result_nc_g(ig,1)
          ENDDO
!$omp end parallel do
!$omp parallel do default(shared), private(ig)
          DO ig = 1, n
             hpsi(lda+ig,im)= hpsi(lda+ig,im) - exxalfa*result_nc_g(ig,2)
          ENDDO
!$omp end parallel do
          !
       ELSE
          !
          CALL fwfft ('CustomWave', RESULT, exx_fft%dfftt)
          !
          !communicate result
          DO ig = 1, n
             result_g(ig) = RESULT(exx_fft%nlt(igk_k(ig,current_k)))
          ENDDO
          CALL mp_sum( result_g(1:n), inter_bgrp_comm)
          !
          !adds it to hpsi
!$omp parallel do default(shared), private(ig)
          DO ig = 1, n
             hpsi(ig,im)=hpsi(ig,im) - exxalfa*result_g(ig)
          ENDDO
!$omp end parallel do
       ENDIF
       !
       ! add non-local \sum_I |beta_I> \alpha_Ii (the sum on i is outside)
       IF(okvan) CALL add_nlxx_pot(lda, hpsi(:,im), xkp, n, igk_k(1,current_k),&
                                       deexx, eps_occ, exxalfa)
       !
    ENDDO &
    LOOP_ON_PSI_BANDS
    !
    IF (noncolin) THEN
       DEALLOCATE(temppsic_nc, result_nc, result_nc_g )
    ELSE
       DEALLOCATE(temppsic, RESULT, result_g )
    ENDIF
    !
    DEALLOCATE(rhoc, vc, fac )
    !
    IF(okvan) DEALLOCATE( deexx)
    !
    !-----------------------------------------------------------------------
  END SUBROUTINE vexx_k
  !-----------------------------------------------------------------------
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
    USE io_files,               ONLY : iunwfc, nwordwfc
    USE buffers,                ONLY : get_buffer
    USE wvfct,                  ONLY : nbnd, npwx, wg, current_k
    USE gvect,                  ONLY : gstart
    USE wavefunctions_module,   ONLY : evc
    USE lsda_mod,               ONLY : lsda, current_spin, isk
    USE klist,                  ONLY : ngk, nks, xk, igk_k
    USE mp_pools,               ONLY : inter_pool_comm
    USE mp_bands,               ONLY : intra_bgrp_comm, intra_bgrp_comm, nbgrp
    USE mp,                     ONLY : mp_sum
    USE becmod,                 ONLY : bec_type, allocate_bec_type, deallocate_bec_type, calbec
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
          CALL get_buffer(psi, nwordwfc, iunwfc, ik)
       ELSE
          psi(1:npwx*npol,1:nbnd) = evc(1:npwx*npol,1:nbnd)
       ENDIF
       !
       IF(okvan)THEN
          ! prepare the |beta> function at k+q
          CALL init_us_2(npw, igk_k(1,ik), xk(:,ik), vkb)
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

    CALL mp_sum( energy, intra_bgrp_comm)
    CALL mp_sum( energy, inter_pool_comm )
    IF(okvan)  CALL deallocate_bec_type(becpsi)
    !
    exxenergy = energy
    !
    CALL stop_clock ('exxenergy')
    !-----------------------------------------------------------------------
  END FUNCTION exxenergy
  !-----------------------------------------------------------------------


  !-----------------------------------------------------------------------
  FUNCTION exxenergy2()
    !-----------------------------------------------------------------------
    !
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
    USE io_files,                ONLY : iunwfc, nwordwfc
    USE buffers,                 ONLY : get_buffer
    USE cell_base,               ONLY : alat, omega, bg, at, tpiba
    USE symm_base,               ONLY : nsym, s
    USE gvect,                   ONLY : ngm, gstart, g, nl
    USE wvfct,                   ONLY : nbnd, npwx, wg
    USE wavefunctions_module,    ONLY : evc
    USE klist,                   ONLY : xk, ngk, nks, nkstot, igk_k
    USE lsda_mod,                ONLY : lsda, current_spin, isk
    USE mp_pools,                ONLY : inter_pool_comm
    USE mp_bands,                ONLY : inter_bgrp_comm, intra_bgrp_comm, nbgrp
    USE mp,                      ONLY : mp_sum
    USE fft_interfaces,          ONLY : fwfft, invfft
    USE gvect,                   ONLY : ecutrho
    USE klist,                   ONLY : wk
    USE uspp,                    ONLY : okvan,nkb,vkb
    USE becmod,                  ONLY : bec_type, allocate_bec_type, deallocate_bec_type, calbec
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
    !
    nrxxs= exx_fft%dfftt%nnr
    ALLOCATE( fac(exx_fft%ngmt) )
    !
    ALLOCATE(temppsic(nrxxs), temppsic_dble(nrxxs),temppsic_aimag(nrxxs))
    ALLOCATE( rhoc(nrxxs) )
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
          CALL get_buffer (evc, nwordwfc, iunwfc, ikk)
       !
       ! prepare the |beta> function at k+q
       CALL init_us_2(npw, igk_k(:,ikk), xkp, vkb)
       ! compute <beta_I|psi_j> at this k+q point, for all band and all projectors
       CALL calbec(npw, vkb, evc, becpsi, nbnd)
       !
       IQ_LOOP : &
       DO iq = 1,nqs
          !
          ikq  = index_xkq(current_ik,iq)
          ik   = index_xk(ikq)
          !
          xkq = xkq_collect(:,ikq)
          !
          CALL g2_convolution(exx_fft%ngmt, exx_fft%gt, xkp, xkq, fac)
          fac(exx_fft%gstart_t:) = 2 * fac(exx_fft%gstart_t:)
          IF ( okvan .and..not.tqr ) CALL qvan_init (exx_fft%ngmt, xkq, xkp)
          !
          jmax = nbnd
          DO jbnd = nbnd,1, -1
             IF ( abs(wg(jbnd,ikk)) < eps_occ) CYCLE
             jmax = jbnd
             exit
          ENDDO
          !
          JBND_LOOP : &
          DO jbnd = 1, jmax     !for each band of psi (the k cycle is outside band)
             !
             temppsic = 0._DP
             !
             l_fft_doubleband = .false.
             l_fft_singleband = .false.
             !
             IF ( mod(jbnd,2)==1 .and. (jbnd+1)<=jmax ) l_fft_doubleband = .true.
             IF ( mod(jbnd,2)==1 .and. jbnd==jmax )     l_fft_singleband = .true.
             !
             IF( l_fft_doubleband ) THEN
!$omp parallel do  default(shared), private(ig)
                DO ig = 1, exx_fft%npwt
                   temppsic( exx_fft%nlt(ig) )  =       evc(ig,jbnd) + (0._DP,1._DP) * evc(ig,jbnd+1)
                   temppsic( exx_fft%nltm(ig) ) = conjg(evc(ig,jbnd) - (0._DP,1._DP) * evc(ig,jbnd+1))
                ENDDO
!$omp end parallel do
             ENDIF
             !
             IF( l_fft_singleband ) THEN
!$omp parallel do  default(shared), private(ig)
                DO ig = 1, exx_fft%npwt
                   temppsic( exx_fft%nlt(ig) )  =       evc(ig,jbnd)
                   temppsic( exx_fft%nltm(ig) ) = conjg(evc(ig,jbnd))
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
             h_ibnd = ibnd_start/2
             IF(mod(ibnd_start,2)==0) THEN
                h_ibnd=h_ibnd-1
                ibnd_loop_start=ibnd_start-1
             ELSE
                ibnd_loop_start=ibnd_start
             ENDIF
             !
             IBND_LOOP_GAM : &
             DO ibnd = ibnd_loop_start, ibnd_end, 2       !for each band of psi
                !
                h_ibnd = h_ibnd + 1
                !
                IF ( ibnd < ibnd_start ) THEN
                   x1 = 0.0_DP
                ELSE
                   x1 = x_occupation(ibnd,ik)
                ENDIF
                !
                IF ( ibnd < ibnd_end ) THEN
                   x2 = x_occupation(ibnd+1,ik)
                ELSE
                   x2 = 0.0_DP
                ENDIF
                IF ( abs(x1) < eps_occ .and. abs(x2) < eps_occ ) CYCLE IBND_LOOP_GAM
                ! calculate rho in real space. Gamma tricks are used.
                ! temppsic is real; tempphic contains band 1 in the real part,
                ! band 2 in the imaginary part; the same applies to rhoc
                !
                IF( mod(jbnd,2) == 0 ) THEN
!$omp parallel do default(shared), private(ir)
                   DO ir = 1, nrxxs
                      rhoc(ir) = exxbuff(ir,h_ibnd,ikq) * temppsic_aimag(ir) / omega
                   ENDDO
!$omp end parallel do
                ELSE
!$omp parallel do default(shared), private(ir)
                   DO ir = 1, nrxxs
                      rhoc(ir) = exxbuff(ir,h_ibnd,ikq) * temppsic_dble(ir) / omega
                   ENDDO
!$omp end parallel do
                ENDIF
                !
                IF(okvan .and.tqr) THEN
                   IF(ibnd>=ibnd_start) &
                   CALL addusxx_r(exx_fft,rhoc, _CX(becxx(ikq)%r(:,ibnd)), _CX(becpsi%r(:,jbnd)))
                   IF(ibnd<ibnd_end) &
                   CALL addusxx_r(exx_fft,rhoc,_CY(becxx(ikq)%r(:,ibnd+1)),_CX(becpsi%r(:,jbnd)))
                ENDIF
                !
                ! bring rhoc to G-space
                CALL fwfft ('Custom', rhoc, exx_fft%dfftt)
                !
                IF(okvan .and..not.tqr) THEN
                   IF(ibnd>=ibnd_start ) &
                      CALL addusxx_g( exx_fft, rhoc, xkq, xkp, 'r', &
                      becphi_r=becxx(ikq)%r(:,ibnd), becpsi_r=becpsi%r(:,jbnd) )
                   IF(ibnd<ibnd_end) &
                      CALL addusxx_g( exx_fft, rhoc, xkq, xkp, 'i', &
                      becphi_r=becxx(ikq)%r(:,ibnd+1), becpsi_r=becpsi%r(:,jbnd) )
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
    CALL mp_sum( energy, inter_bgrp_comm )
    CALL mp_sum( energy, intra_bgrp_comm )
    CALL mp_sum( energy, inter_pool_comm )
    !
    exxenergy2_gamma = energy
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
    USE io_files,                ONLY : iunwfc, nwordwfc
    USE buffers,                 ONLY : get_buffer
    USE cell_base,               ONLY : alat, omega, bg, at, tpiba
    USE symm_base,               ONLY : nsym, s
    USE gvect,                   ONLY : ngm, gstart, g, nl
    USE wvfct,                   ONLY : nbnd, npwx, wg
    USE wavefunctions_module,    ONLY : evc
    USE klist,                   ONLY : xk, ngk, nks, nkstot, igk_k
    USE lsda_mod,                ONLY : lsda, current_spin, isk
    USE mp_pools,                ONLY : inter_pool_comm
    USE mp_bands,                ONLY : inter_bgrp_comm, intra_bgrp_comm, nbgrp
    USE mp,                      ONLY : mp_sum
    USE fft_interfaces,          ONLY : fwfft, invfft
    USE gvect,                   ONLY : ecutrho
    USE klist,                   ONLY : wk
    USE uspp,                    ONLY : okvan,nkb,vkb
    USE becmod,                  ONLY : bec_type, allocate_bec_type, deallocate_bec_type, calbec
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
    COMPLEX(DP), ALLOCATABLE :: temppsic(:)
    COMPLEX(DP), ALLOCATABLE :: temppsic_nc(:,:)
    COMPLEX(DP), ALLOCATABLE :: rhoc(:)
    REAL(DP),    ALLOCATABLE :: fac(:)
    INTEGER  :: npw, jbnd, ibnd, ik, ikk, ig, ikq, iq, ir
    INTEGER  :: h_ibnd, nrxxs, current_ik, ibnd_loop_start
    REAL(DP) :: x1, x2
    REAL(DP) :: xkq(3), xkp(3), vc
    ! temp array for vcut_spheric
    INTEGER, EXTERNAL :: global_kpoint_index
    !
    TYPE(bec_type) :: becpsi
    COMPLEX(DP), ALLOCATABLE :: psi_t(:), prod_tot(:)
    !
    nrxxs = exx_fft%dfftt%nnr
    ALLOCATE( fac(exx_fft%ngmt) )
    !
    IF (noncolin) THEN
       ALLOCATE(temppsic_nc(nrxxs,npol))
    ELSE
       ALLOCATE(temppsic(nrxxs))
    ENDIF
    ALLOCATE( rhoc(nrxxs) )
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
          CALL get_buffer (evc, nwordwfc, iunwfc, ikk)
       !
       ! prepare the |beta> function at k+q
       CALL init_us_2(npw, igk_k(:,ikk), xkp, vkb)
       ! compute <beta_I|psi_j> at this k+q point, for all band and all projectors
       CALL calbec(npw, vkb, evc, becpsi, nbnd)
       !
       JBND_LOOP : &
       DO jbnd = 1, nbnd     !for each band of psi (the k cycle is outside band)
          !
          IF ( abs(wg(jbnd,ikk)) < eps_occ) CYCLE
          !
          IF (noncolin) THEN
             temppsic_nc = 0.0_DP
          ELSE
             temppsic    = 0.0_DP
          ENDIF
          !
          IF (noncolin) THEN
             !
!$omp parallel do default(shared), private(ig)
             DO ig = 1, npw
                temppsic_nc(exx_fft%nlt(igk_k(ig,ikk)),1) = evc(ig,jbnd)
             ENDDO
!$omp end parallel do
!$omp parallel do default(shared), private(ig)
             DO ig = 1, npw
                temppsic_nc(exx_fft%nlt(igk_k(ig,ikk)),2) = evc(npwx+ig,jbnd)
             ENDDO
!$omp end parallel do
             !
             CALL invfft ('CustomWave', temppsic_nc(:,1), exx_fft%dfftt)
             CALL invfft ('CustomWave', temppsic_nc(:,2), exx_fft%dfftt)
             !
          ELSE
!$omp parallel do default(shared), private(ig)
             DO ig = 1, npw
                temppsic(exx_fft%nlt(igk_k(ig,ikk))) = evc(ig,jbnd)
             ENDDO
!$omp end parallel do
             !
             CALL invfft ('CustomWave', temppsic, exx_fft%dfftt)
             !
          ENDIF
          !
          IQ_LOOP : &
          DO iq = 1,nqs
             !
             ikq  = index_xkq(current_ik,iq)
             ik   = index_xk(ikq)
             !
             xkq = xkq_collect(:,ikq)
             !
             CALL g2_convolution(exx_fft%ngmt, exx_fft%gt, xkp, xkq, fac)
             IF ( okvan .and..not.tqr ) CALL qvan_init (exx_fft%ngmt, xkq, xkp)
             !
             IBND_LOOP_K : &
             DO ibnd = ibnd_start, ibnd_end
                !
                IF ( abs(x_occupation(ibnd,ik)) < eps_occ) CYCLE
                !
                ! load the phi at this k+q and band
                IF (noncolin) THEN
                   !
!$omp parallel do  default(shared), private(ir)
                   DO ir = 1, nrxxs
                      rhoc(ir)=(conjg(exxbuff(ir      ,ibnd,ikq))*temppsic_nc(ir,1) + &
                                conjg(exxbuff(ir+nrxxs,ibnd,ikq))*temppsic_nc(ir,2) )/omega
                   ENDDO
!$omp end parallel do
                ELSE
                   !calculate rho in real space
!$omp parallel do  default(shared), private(ir)
                   DO ir = 1, nrxxs
                      rhoc(ir)=conjg(exxbuff(ir,ibnd,ikq))*temppsic(ir) / omega
                   ENDDO
!$omp end parallel do
                ENDIF
                ! augment the "charge" in real space
                IF(okvan .and. tqr) CALL addusxx_r(exx_fft,rhoc, becxx(ikq)%k(:,ibnd), becpsi%k(:,jbnd))
                !
                ! bring rhoc to G-space
                CALL fwfft ('Custom', rhoc, exx_fft%dfftt)
                ! augment the "charge" in G space
                IF(okvan .and. .not. tqr) &
                   CALL addusxx_g(exx_fft, rhoc, xkq, xkp, 'c', &
                   becphi_c=becxx(ikq)%k(:,ibnd),becpsi_c=becpsi%k(:,jbnd))
                !
                vc = 0.0_DP
!$omp parallel do  default(shared), private(ig), reduction(+:vc)
                DO ig=1,exx_fft%ngmt
                   vc = vc + fac(ig) * dble(rhoc(exx_fft%nlt(ig)) * &
                                      conjg(rhoc(exx_fft%nlt(ig))))
                ENDDO
!$omp end parallel do
                vc = vc * omega * x_occupation(ibnd,ik) / nqs
                !
                energy = energy - exxalfa * vc * wg(jbnd,ikk)
                !
                IF(okpaw) THEN
                   energy = energy +exxalfa*x_occupation(ibnd,ik)/nqs*wg(jbnd,ikk) &
                              *PAW_xx_energy(becxx(ikq)%k(:,ibnd), becpsi%k(:,jbnd))
                ENDIF
                !
             ENDDO &
             IBND_LOOP_K
             IF ( okvan .and..not.tqr ) CALL qvan_clean ( )
          ENDDO &
          IQ_LOOP
       ENDDO &
       JBND_LOOP
    ENDDO &
    IKK_LOOP
    !
    IF (noncolin) THEN
       DEALLOCATE(temppsic_nc)
    ELSE
       DEALLOCATE(temppsic)
    ENDIF
    !
    DEALLOCATE(rhoc, fac )
    CALL deallocate_bec_type(becpsi)
    !
    CALL mp_sum( energy, inter_bgrp_comm )
    CALL mp_sum( energy, intra_bgrp_comm )
    CALL mp_sum( energy, inter_pool_comm )
    !
    exxenergy2_k = energy
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
     USE mp_bands,       ONLY : intra_bgrp_comm
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
     CALL mp_sum(  div, intra_bgrp_comm )
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
    USE io_files,             ONLY : iunwfc, nwordwfc
    USE buffers,              ONLY : get_buffer
    USE cell_base,            ONLY : alat, omega, bg, at, tpiba
    USE symm_base,            ONLY : nsym, s
    USE wvfct,                ONLY : nbnd, npwx, wg, current_k
    USE wavefunctions_module, ONLY : evc
    USE klist,                ONLY : xk, ngk, nks, igk_k
    USE lsda_mod,             ONLY : lsda, current_spin, isk
    USE gvect,                ONLY : g, nl
    USE mp_pools,             ONLY : npool, inter_pool_comm
    USE mp_bands,             ONLY : inter_bgrp_comm, intra_bgrp_comm
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
    COMPLEX(DP),ALLOCATABLE :: tempphic(:), temppsic(:), RESULT(:)
    COMPLEX(DP),ALLOCATABLE :: tempphic_nc(:,:), temppsic_nc(:,:), &
                               result_nc(:,:)
    COMPLEX(DP),ALLOCATABLE :: rhoc(:)
    REAL(DP),    ALLOCATABLE :: fac(:), fac_tens(:,:,:), fac_stress(:)
    INTEGER  :: npw, jbnd, ibnd, ik, ikk, ig, ir, ikq, iq, isym
    INTEGER  :: h_ibnd, nqi, iqi, beta, nrxxs, ngm
    INTEGER  :: ibnd_loop_start
    REAL(DP) :: x1, x2
    REAL(DP) :: qq, xk_cryst(3), sxk(3), xkq(3), vc(3,3), x, q(3)
    ! temp array for vcut_spheric
    REAL(DP) :: delta(3,3)

    CALL start_clock ('exx_stress')

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
    nqi=nqs
    !
    ! loop over k-points
    DO ikk = 1, nks
        current_k = ikk
        IF (lsda) current_spin = isk(ikk)
        npw = ngk(ikk)

        IF (nks > 1) &
            CALL get_buffer(evc, nwordwfc, iunwfc, ikk)

        ! loop over bands
        DO jbnd = 1, nbnd
            !
            temppsic(:) = ( 0._dp, 0._dp )
!$omp parallel do default(shared), private(ig)
            DO ig = 1, npw
                temppsic(exx_fft%nlt(igk_k(ig,ikk))) = evc(ig,jbnd)
            ENDDO
!$omp end parallel do
            !
            IF(gamma_only) THEN
!$omp parallel do default(shared), private(ig)
                DO ig = 1, npw
                    temppsic(exx_fft%nltm(igk_k(ig,ikk))) = conjg(evc(ig,jbnd))
                ENDDO
!$omp end parallel do
            ENDIF

            CALL invfft ('CustomWave', temppsic, exx_fft%dfftt)

            DO iqi = 1, nqi
                !
                iq=iqi
                !
                ikq  = index_xkq(current_k,iq)
                ik   = index_xk(ikq)
                isym = abs(index_sym(ikq))

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
                  q(1)= xk(1,current_k) - xkq(1) + g(1,ig)
                  q(2)= xk(2,current_k) - xkq(2) + g(2,ig)
                  q(3)= xk(3,current_k) - xkq(3) + g(3,ig)

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

                  ELSEIF (gau_scrlen > 0) THEN
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

                IF (gamma_only) THEN
                    !
                    h_ibnd = ibnd_start/2
                    !
                    IF(mod(ibnd_start,2)==0) THEN
                      h_ibnd=h_ibnd-1
                      ibnd_loop_start=ibnd_start-1
                    ELSE
                      ibnd_loop_start=ibnd_start
                    ENDIF
                    !
                    DO ibnd = ibnd_loop_start, ibnd_end, 2     !for each band of psi
                        !
                        h_ibnd = h_ibnd + 1
                        !
                        IF( ibnd < ibnd_start ) THEN
                            x1 = 0._dp
                        ELSE
                            x1 = x_occupation(ibnd,  ik)
                        ENDIF

                        IF( ibnd == ibnd_end) THEN
                            x2 = 0._dp
                        ELSE
                            x2 = x_occupation(ibnd+1,  ik)
                        ENDIF
                        IF ( abs(x1) < eps_occ .and. abs(x2) < eps_occ ) CYCLE
                        !
                        ! calculate rho in real space
!$omp parallel do default(shared), private(ir)
                        DO ir = 1, nrxxs
                            tempphic(ir) = exxbuff(ir,h_ibnd,ikq)
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

                    DO ibnd = ibnd_start, ibnd_end    !for each band of psi
                      !
                      IF ( abs(x_occupation(ibnd,ik)) < 1.d-6) CYCLE
                      !
                      ! calculate rho in real space
!$omp parallel do default(shared), private(ir)
                      DO ir = 1, nrxxs
                          tempphic(ir) = exxbuff(ir,ibnd,ikq)
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

            ENDDO ! iqi
        ENDDO ! jbnd
    ENDDO ! ikk

    DEALLOCATE(tempphic, temppsic, rhoc, fac, fac_tens, fac_stress )
    !
    CALL mp_sum( exx_stress_, intra_bgrp_comm )
    CALL mp_sum( exx_stress_, inter_bgrp_comm )
    CALL mp_sum( exx_stress_, inter_pool_comm )
    exx_stress = exx_stress_

    CALL stop_clock ('exx_stress')
    !-----------------------------------------------------------------------
  END FUNCTION exx_stress
  !-----------------------------------------------------------------------
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE matprt(label,n,m,A)
IMPLICIT NONE
  INTEGER :: n,m,i
  real*8 :: A(n,m)
  CHARACTER(len=50) :: frmt
  CHARACTER(len=*) :: label

  WRITE(*,'(A)') label
  frmt = ' '
  WRITE(frmt,'(A,I4,A)') '(',m,'f16.10)'
  DO i = 1,n
    WRITE(*,frmt) A(i,:)
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
  USE io_global,  ONLY : stdout
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
    WRITE(*,'(4(A,I3),A,I9,A,f12.6)') 'vexxace: nbnd=', nbnd, ' nbndproj=',nbndproj, &
                                              ' k=',current_k,' spin=',current_spin,' npw=',nnpw, ' E=',exxe
  ELSE
    WRITE(*,'(4(A,I3),A,I9)')         'vexxace: nbnd=', nbnd, ' nbndproj=',nbndproj, &
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
    WRITE(*,'(A,f16.8,A)') string//label, ee, ' Ry'
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
  WRITE(*,'(3(A,I3),A,I9,A,f12.6)') 'aceinit_k: nbnd=', nbnd, ' nbndproj=',nbndproj, &
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
    !write(*,'(A,f16.8,A)') string//label, ee, ' Ry'
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

  WRITE(*,'(A)') label//'(real)'
  frmt = ' '
  WRITE(frmt,'(A,I4,A)') '(',m,'f12.6)'
  DO i = 1,n
    WRITE(*,frmt) dreal(A(i,:))
  ENDDO

  WRITE(*,'(A)') label//'(imag)'
  frmt = ' '
  WRITE(frmt,'(A,I4,A)') '(',m,'f12.6)'
  DO i = 1,n
    WRITE(*,frmt) aimag(A(i,:))
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
    WRITE(*,'(3(A,I3),A,I9,A,f12.6)') 'vexxace_k: nbnd=', nbnd, ' nbndproj=',nbndproj, &
                   ' k=',current_k,' npw=',nnpw, ' Ex(k)=',exxe
  ELSE
    WRITE(*,'(3(A,I3),A,I9)') 'vexxace_k: nbnd=', nbnd, ' nbndproj=',nbndproj, &
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
