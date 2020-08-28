!
! Copyright (C) 2005-2018 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------
MODULE exx_base
  !--------------------------------------
  !! Basic variables and subroutines for calculation of exact exchange (EXX).
  !
  USE kinds,                ONLY : DP
  USE coulomb_vcut_module,  ONLY : vcut_init, vcut_type, vcut_info, &
                                   vcut_get,  vcut_spheric_get
  USE noncollin_module,     ONLY : noncolin
  USE io_global,            ONLY : ionode
  !
  USE control_flags,        ONLY : gamma_only
  USE fft_types,            ONLY : fft_type_descriptor
  USE stick_base,           ONLY : sticks_map, sticks_map_deallocate
  !
  USE io_global,            ONLY : stdout
  !
  IMPLICIT NONE
  !
  SAVE
  !
  ! ... variables defining the auxiliary k-point grid
  !     used in X BZ integration
  !
  INTEGER :: nq1=1, nq2=1, nq3=1
  !! integers defining the X integration mesh
  INTEGER :: nqs=1
  !! number of points in the q-grid
  INTEGER :: nkqs
  !! total number of different k+q
  REAL(DP), ALLOCATABLE :: xkq_collect(:,:)
  !! xkq(3,nkqs) the auxiliary k+q set
  !
  ! ... let xk(:,ik) + xq(:,iq) = xkq(:,ikq) = S(isym)*xk(ik') + G
  !
  !     index_xkq(ik,iq) = ikq
  !     index_xk(ikq)    = ik'
  !     index_sym(ikq)   = isym
  !
  INTEGER, ALLOCATABLE :: index_xkq(:,:)
  !! index_xkq(nks,nqs)
  INTEGER, ALLOCATABLE :: index_xk(:)
  !! index_xk(nkqs)
  INTEGER, ALLOCATABLE :: index_sym(:)
  !! index_sym(nkqs)
  INTEGER, ALLOCATABLE :: rir(:,:)
  !! rotations to take k to q
  INTEGER, ALLOCATABLE :: working_pool(:)
  !! idx of pool that has the k+q-point before rotation
  !
  ! ... Internal:
  LOGICAL :: exx_grid_initialized = .FALSE.
  !
  ! ... variables to deal with Coulomb divergence
  !     and related issues
  REAL(DP) :: eps  = 1.d-6
  REAL(DP) :: eps_qdiv = 1.d-8 ! |q| > eps_qdiv
  REAL(DP) :: exxdiv = 0._dp
  CHARACTER(32) :: exxdiv_treatment  = ' '
  !
  ! ... x_gamma_extrapolation
  LOGICAL :: x_gamma_extrapolation =.TRUE.
  LOGICAl :: on_double_grid =.FALSE.
  REAL(DP) :: grid_factor = 1.d0 !8.d0/7.d0
  !
  LOGICAL :: use_regularization = .TRUE.
  !! Gygi-Baldereschi
  !
  REAL(DP) :: yukawa = 0._dp
  !! Yukawa method
  !
  REAL(DP) :: erfc_scrlen = 0._dp
  !! erfc screening
  !
  REAL(DP) :: erf_scrlen = 0._dp
  !! erf screening
  !
  REAL(DP) :: gau_scrlen = 0.d0
  !! gau-pbe screening
  !
  ! ... cutoff techniques
  LOGICAL :: use_coulomb_vcut_ws = .FALSE.
  LOGICAL :: use_coulomb_vcut_spheric = .FALSE.
  REAL(DP) :: ecutvcut
  TYPE(vcut_type) :: vcut
  !
  REAL(DP), ALLOCATABLE :: coulomb_fac(:,:,:)
  !! the Coulomb factor is reused between iterations
  !
  LOGICAL, ALLOCATABLE :: coulomb_done(:,:)
  !! list of which Coulomb factors have been calculated already
  !
 CONTAINS
  !
  !------------------------------------------------------------------------
  SUBROUTINE exx_mp_init()
    !------------------------------------------------------------------------
    !! 1) Setup the orthopool communicators, which include the (n-1)-th CPU of
    !!    each pool (i.e. orthopool 0, contains CPU 0 of each pool).
    !
    !! 2) Setup global variable \(\text{working_pool}\) which contains the index of the
    !!    pool which has the local copy prior to rotation of the \(\text{ikq}\)-th
    !!    pool.  
    !
    !! 3) \(\text{working_pool(ikq)}\) is also index of the CPU in each orthopool
    !!    which has to broadcast the wavefunction at k-point \(\text{ikq}\).
    !
    USE mp_images,      ONLY : intra_image_comm
    USE mp_pools,       ONLY : my_pool_id
    USE mp_orthopools,  ONLY : mp_start_orthopools, intra_orthopool_comm
    USE mp,             ONLY : mp_sum
    USE klist,          ONLY : nkstot
    !
    IMPLICIT NONE
    !
    INTEGER :: ikq, current_ik, ik
    INTEGER, EXTERNAL :: local_kpoint_index
    !
    IF (nkqs <= 0) CALL errore( "exx_mp_init", "exx_grid_init must be called first!", 1 )
    !
    ! ... mp_start_orthopools contains a check to avoid double initialisation
    CALL mp_start_orthopools( intra_image_comm )
    !
    IF (ALLOCATED(working_pool)) DEALLOCATE( working_pool )
    ALLOCATE( working_pool(nkqs) )
    !
    DO ikq = 1, nkqs
      DO current_ik = 1, nkstot
        ! find the index of point k+q from the global list
        IF (index_xk(ikq) == current_ik) EXIT
      ENDDO
      IF (current_ik>nkstot) CALL errore( "exx_mp_init", "could not find the local index", ikq )
      ! Find the pool-index of the global point, returns -1 if the point
      ! was not in my pool
      ik = local_kpoint_index(nkstot, current_ik)
      !
      ! Set a temporary variable to the index of the pool who did the work
      IF( ik > 0 ) THEN
        working_pool(ikq) = my_pool_id
      ELSE
        working_pool(ikq) = 0
      ENDIF
      !
    ENDDO    
      !
      ! Collect the variable among all the nth CPUs of each pool,
      ! i.e. inside each orto-pool group
    CALL mp_sum( working_pool, intra_orthopool_comm )
    !
  END SUBROUTINE exx_mp_init
  !
  !------------------------------------------------------------------------
  SUBROUTINE exx_grid_init( reinit )
    !------------------------------------------------------------------------
    !!
    !
    USE symm_base,         ONLY : nsym, s
    USE cell_base,         ONLY : bg, at, tpiba
    USE spin_orb,          ONLY : domag
    USE noncollin_module,  ONLY : nspin_lsda
    USE klist,             ONLY : xk, wk, nkstot, nks, qnorm
    USE wvfct,             ONLY : nbnd
    USE start_k,           ONLY : nk1,nk2,nk3
    USE control_flags,     ONLY : iverbosity
    USE mp_pools,          ONLY : inter_pool_comm
    !
    IMPLICIT NONE
    !
    LOGICAL, OPTIONAL :: reinit
    !! reinitialize exx if .TRUE.
    !
    ! ... local variables
    !
    CHARACTER(13) :: sub_name='exx_grid_init'
    INTEGER :: iq1, iq2, iq3, isym, ik, ikq, iq, max_nk, temp_nkqs, idx, sign_
    INTEGER :: nqx(3)
    INTEGER, ALLOCATABLE :: temp_index_xk(:), temp_index_sym(:)
    INTEGER, ALLOCATABLE :: temp_index_ikq(:)
    REAL(DP), ALLOCATABLE :: temp_xkq(:,:), xk_collect(:,:)
    LOGICAL :: xk_not_found
    REAL(DP) :: sxk(3), dxk(3), xk_cryst(3)
    REAL(DP) :: dq1, dq2, dq3
    CHARACTER(LEN=6), EXTERNAL :: int_to_char
    !
    CALL start_clock ( 'exx_grid' )
    !
    IF ( PRESENT (reinit) ) THEN
       IF ( reinit ) THEN
          IF (ALLOCATED(xkq_collect))  DEALLOCATE( xkq_collect )
          IF (ALLOCATED(index_xk)   )  DEALLOCATE( index_xk    )
          IF (ALLOCATED(index_sym)  )  DEALLOCATE( index_sym   )
          exx_grid_initialized = .FALSE.
          nkqs = 0
       ENDIF
    ENDIF
    !
    IF (nq1 <= 0) nq1 = nk1
    IF (nq2 <= 0) nq2 = nk2
    IF (nq3 <= 0) nq3 = nk3
    !
    IF (nkstot == nspin_lsda) THEN
      nq1=1; nq2=1; nq3=1
    ENDIF
    !
    IF (ANY( (/nq1,nq2,nq3/) <=0 )) CALL errore( 'exx_grid_init', "wrong EXX q grid", 1 )
    !
    IF (exx_grid_initialized) CALL errore( 'exx_grid_init', "grid already initialized", 1 )
    exx_grid_initialized = .TRUE.
    !
    ! definitions and checks
    !
    grid_factor = 1._dp
    IF (x_gamma_extrapolation) grid_factor = 8.d0/7.d0
    !
    nqs = nq1 * nq2 * nq3
    !
    ! all processors on all pools need to have access to all k+q points
    !
    ALLOCATE( xk_collect(3,nkstot) )
    !
    CALL poolcollect( 3, nks, xk, nkstot, xk_collect )
    !
    ! set a safe limit as the maximum number of auxiliary points we may need
    ! and allocate auxiliary arrays
    max_nk = nkstot * MIN(48, 2 * nsym)
    ALLOCATE( temp_index_xk(max_nk), temp_index_sym(max_nk) )
    ALLOCATE( temp_index_ikq(max_nk) )
    ALLOCATE( temp_xkq(3,max_nk) )
    !
    ! find all k-points equivalent by symmetry to the points in the k-list
    !
    temp_nkqs = 0
    !
    DO isym = 1, nsym
      DO ik = 1, nkstot
          !
          ! go to crystalline coordinates
          xk_cryst(:) = xk_collect(:,ik)
          CALL cryst_to_cart( 1, xk_cryst, at, -1 )
          ! rotate with this sym.op.
          sxk(:) = s(:,1,isym)*xk_cryst(1) + &
                   s(:,2,isym)*xk_cryst(2) + &
                   s(:,3,isym)*xk_cryst(3)
          ! add sxk to the auxiliary list IF it is not already present
          xk_not_found = .TRUE.
          ! *** do-loop skipped the first time because temp_nkqs == 0
          DO ikq = 1, temp_nkqs
            IF (xk_not_found ) THEN
                dxk(:) = sxk(:)-temp_xkq(:,ikq) - NINT(sxk(:)-temp_xkq(:,ikq))
                IF ( ABS(dxk(1))<=eps .AND. &
                     ABS(dxk(2))<=eps .AND. &
                     ABS(dxk(3))<=eps ) xk_not_found = .FALSE.
            ENDIF
          ENDDO
          !
          IF (xk_not_found) THEN
            temp_nkqs                 = temp_nkqs + 1
            temp_xkq(:,temp_nkqs)     = sxk(:)
            temp_index_xk(temp_nkqs)  = ik
            temp_index_sym(temp_nkqs) = isym
          ENDIF
          !
          sxk(:) = - sxk(:)
          xk_not_found = .TRUE.
          DO ikq = 1, temp_nkqs
            IF (xk_not_found ) THEN
                dxk(:) = sxk(:) - temp_xkq(:,ikq) - NINT(sxk(:) - temp_xkq(:,ikq))
                IF ( ABS(dxk(1))<=eps .AND. &
                     ABS(dxk(2))<=eps .AND. &
                     ABS(dxk(3))<=eps ) xk_not_found = .FALSE.
            ENDIF
          ENDDO
          !
          IF (xk_not_found .AND. .NOT. (noncolin.AND.domag) ) THEN
            temp_nkqs                 = temp_nkqs + 1
            temp_xkq(:,temp_nkqs)     = sxk(:)
            temp_index_xk(temp_nkqs)  = ik
            temp_index_sym(temp_nkqs) =-isym
          ENDIF
          !
      ENDDO
    ENDDO
    !
    ! Find good q-point grid. Decrease the nqX until a good grid is found or
    ! until it is 1 x 1 x 1 (always good)
    idx = 1
    sign_ = -1
    nqx = (/nq1, nq2, nq3/)
    DO WHILE (.TRUE.)
      CALL exx_qgrid_init(temp_nkqs, xk_collect, temp_xkq, &
                          nkqs, temp_index_ikq, dxk)

      ! Good q-point mesh
      IF (ALL(ABS(dxk) < eps ) ) THEN
        !
        IF (idx > 1) &
          WRITE(stdout, '(5x,a)') "EXX: WARNING: q-point mesh has been updated!"
        !
        WRITE(stdout, '(5x,a,3i5)') "EXX: q-point mesh: ", nq1, nq2, nq3
        EXIT ! DO WHILE
      ENDIF
      !
      ! Try q-points around the input mesh, prioritizing smaller mesh
      !
      nq1 = nqx(1) + idx * sign_
      nq2 = nqx(2) + idx * sign_
      nq3 = nqx(3) + idx * sign_
      !
      ! Ensure no values smaller than 1
      IF (nq1 < 1) nq1 = 1
      IF (nq2 < 1) nq2 = 1
      IF (nq3 < 1) nq3 = 1
      !
      ! Enforce nqX <= nkX. This is important for surfaces to keep the
      ! Z q-point 1.
      !
      IF (nq1 > nk1) nq1 = nk1
      IF (nq2 > nk2) nq2 = nk2
      IF (nq3 > nk3) nq3 = nk3
      !
      nqs = nq1 * nq2 * nq3
      !
      sign_ = -1 * sign_
      !
      ! Increase idx every other time sign is changed
      IF (sign_ < 0) idx = idx + 1
      !
    ENDDO
    !
    ! allocate and fill the arrays xkq(3,nkqs), index_xk(nkqs) and index_sym(nkqs)
    ! NOTE: nkqs will be redefined as nspin_lsda*nkqs later 
    !
    ALLOCATE( xkq_collect(3,nspin_lsda*nkqs), index_xk(nspin_lsda*nkqs),  &
              index_sym(nspin_lsda*nkqs) )
    !
    DO ik = 1, nkqs
      ikq               = temp_index_ikq(ik)
      xkq_collect(:,ik) = temp_xkq(:,ikq)
      index_xk(ik)      = temp_index_xk(ikq)
      index_sym(ik)     = temp_index_sym(ikq)
    ENDDO
    !
    CALL cryst_to_cart( nkqs, xkq_collect, bg, +1 )
    !
    IF (nkqs > 1) THEN
      WRITE(stdout, '(5x,3a)') "EXX: setup a grid of "//TRIM(int_to_char(nkqs))&
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
      WRITE(stdout, '(5X,"EXX: grid of k+q points same as grid of k-points")')
    ENDIF
    !
    ! if nspin == 2, the kpoints are repeated in couples (spin up, spin down)
    IF (nspin_lsda == 2) THEN
      DO ik = 1, nkstot/2
          DO iq = 1, nqs
            index_xkq(nkstot/2+ik,iq) = index_xkq(ik,iq) + nkqs
          ENDDO
      ENDDO
      !
      DO ikq = 1, nkqs
          xkq_collect(:,ikq+nkqs) = xkq_collect(:,ikq)
          index_xk(ikq + nkqs)  = index_xk(ikq) + nkstot/2
          index_sym(ikq + nkqs) = index_sym(ikq)
      ENDDO
      nkqs = 2 * nkqs
    ENDIF
    !
    ! clean up
    DEALLOCATE( temp_index_xk, temp_index_sym, temp_index_ikq, temp_xkq )
    !
    ! check that everything is what it should be
    CALL exx_grid_check( xk_collect(:,:) )
    DEALLOCATE( xk_collect )
    !
    ! qnorm = max |q|, used in allocate_nlpot to compute the maximum size
    !         of some arrays (e.g. qrad) - beware: needed for US/PAW+EXX
    !
    qnorm = 0.0_dp
    DO iq = 1, nkqs
       DO ik = 1, nks
          qnorm = MAX(qnorm, SQRT( SUM((xk(:,ik)-xkq_collect(:,iq))**2) ))
       ENDDO
    ENDDO
    qnorm = qnorm * tpiba
    !
    CALL stop_clock( 'exx_grid' )
    !
    RETURN
    !
  END SUBROUTINE exx_grid_init
  !
  !
  !------------------------------------------------------------------------
  SUBROUTINE exx_qgrid_init(temp_nkqs, xk_collect, temp_xkq, nkqs, temp_index_ikq, dxk)
    !------------------------------------------------------------------------
    !! Generate q-point mesh compatible with the k-point mesh
    !
    USE klist,     ONLY : nkstot, xk
    USE cell_base, ONLY : at
    USE symm_base, ONLY : nsym
    !
    IMPLICIT NONE
    INTEGER, INTENT (IN) :: temp_nkqs
    REAL(DP), INTENT(IN) :: xk_collect(:,:), temp_xkq(:,:)
    !
    REAL(DP), INTENT (OUT) :: dxk(:)
    INTEGER, INTENT (OUT) :: temp_index_ikq(:), nkqs
    !
    INTEGER :: ik, ikq, iq, iq1, iq2, iq3, j, max_nk
    INTEGER, ALLOCATABLE :: new_ikq(:)
    REAL(DP) :: sxk(3), xk_cryst(3), dq1, dq2, dq3
    LOGICAL :: xk_not_found
    !
    max_nk = nkstot * MIN(48, 2 * nsym)
    ALLOCATE( new_ikq(max_nk) )

    IF ( ALLOCATED(index_xkq) ) DEALLOCATE( index_xkq )
    ALLOCATE( index_xkq(nkstot, nqs) )
    !
    nkqs = 0
    new_ikq(:) = 0
    !
    ! define the q-mesh step-sizes
    !
    dq1 = 1._dp / DBLE(nq1)
    dq2 = 1._dp / DBLE(nq2)
    dq3 = 1._dp / DBLE(nq3)
    !
    DO ik = 1, nkstot
      ! go to crystalline coordinates
      xk_cryst(:) = xk_collect(:,ik)
      CALL cryst_to_cart( 1, xk_cryst, at, -1 )
      !
      iq = 0
      !
      DO iq1 = 1, nq1
        sxk(1) = xk_cryst(1) + (iq1-1) * dq1
        DO iq2 = 1, nq2
          sxk(2) = xk_cryst(2) + (iq2-1) * dq2
          DO iq3 = 1, nq3
              sxk(3) = xk_cryst(3) + (iq3-1) * dq3
              iq = iq + 1
              xk_not_found = .TRUE.
              !
              DO ikq = 1, temp_nkqs
                IF ( xk_not_found ) THEN
                    dxk(:) = sxk(:)-temp_xkq(:,ikq) - NINT(sxk(:)-temp_xkq(:,ikq))
                    IF ( ALL(ABS(dxk) < eps ) ) THEN
                        xk_not_found = .FALSE.
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
                DEALLOCATE( new_ikq )
                RETURN
              ENDIF
              !
          ENDDO
        ENDDO
      ENDDO
      !
    ENDDO
    !
    DEALLOCATE( new_ikq )
  END SUBROUTINE exx_qgrid_init
  !
  !
  !------------------------------------------------------------------------
  SUBROUTINE exx_div_check()
    !------------------------------------------------------------------------
    !! EXX singularity treatment.
    !
    USE cell_base,  ONLY : at, alat
    USE funct,      ONLY : get_screening_parameter
    !
    IMPLICIT NONE
    !
    REAL(DP) :: atws(3,3)
    CHARACTER(13) :: sub_name='exx_div_check'
    !
    SELECT CASE ( TRIM(exxdiv_treatment) )
    CASE ( "gygi-baldereschi", "gygi-bald", "g-b", "gb" )
      !
      use_regularization = .TRUE.
      !
    CASE ( "vcut_ws" )
      !
      use_regularization = .TRUE.
      use_coulomb_vcut_ws = .TRUE.
      IF ( x_gamma_extrapolation ) &
            CALL errore( sub_name, 'cannot USE x_gamm_extrap and vcut_ws', 1 )
      !
    CASE ( "vcut_spherical" )
      !
      use_regularization = .TRUE.
      use_coulomb_vcut_spheric = .TRUE.
      IF ( x_gamma_extrapolation ) &
            CALL errore( sub_name, 'cannot USE x_gamm_extrap and vcut_spherical', 1 )
      !
    CASE ( "none" )
      use_regularization = .FALSE.
      !
    CASE DEFAULT
      CALL errore( sub_name, 'invalid exxdiv_treatment: '//TRIM(exxdiv_treatment), 1 )
      !
    END SELECT
    !
    ! Set variables for Coulomb vcut
    ! NOTE: some memory is allocated inside this routine (in variable vcut)
    !       and should be deallocated at the end of the run
    !
    IF ( use_coulomb_vcut_ws .OR. use_coulomb_vcut_spheric ) THEN
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
       IF (ionode) CALL vcut_info( stdout, vcut )
       !
    ENDIF
    !
    RETURN
    !
  END SUBROUTINE exx_div_check
  !
  !
  !------------------------------------------------------------------------
  SUBROUTINE exx_grid_check( xk_collect )
    !------------------------------------------------------------------------
    !
    USE symm_base,  ONLY : s
    USE cell_base,  ONLY : at
    USE klist,      ONLY : nkstot, xk
    !
    IMPLICIT NONE
    !
    REAL(DP), INTENT(IN) :: xk_collect(:,:) 
    !!
    !
    ! ... local variables
    !
    REAL(DP) :: sxk(3), dxk(3), xk_cryst(3), xkk_cryst(3)
    INTEGER :: iq1, iq2, iq3, isym, ik, ikk, ikq, iq
    REAL(DP) :: dq1, dq2, dq3
    !
    dq1 = 1._dp/DBLE(nq1)
    dq2 = 1._dp/DBLE(nq2)
    dq3 = 1._dp/DBLE(nq3)
    !
    DO ik = 1, nkstot
      xk_cryst(:) = xk_collect(:,ik)
      CALL cryst_to_cart( 1, xk_cryst, at, -1 )
      !
      iq = 0
      DO iq1=1, nq1
        sxk(1) = xk_cryst(1) + (iq1-1) * dq1
        DO iq2 =1, nq2
          sxk(2) = xk_cryst(2) + (iq2-1) * dq2
          DO iq3 =1, nq3
              sxk(3) = xk_cryst(3) + (iq3-1) * dq3
              iq = iq + 1
              !
              ikq  = index_xkq(ik,iq)
              ikk  = index_xk(ikq)
              isym = index_sym(ikq)
              !
              xkk_cryst(:) = at(1,:)*xk_collect(1,ikk) + &
                             at(2,:)*xk_collect(2,ikk) + &
                             at(3,:)*xk_collect(3,ikk)
              IF (isym < 0 ) xkk_cryst(:) = - xkk_cryst(:)
              isym = ABS(isym)
              dxk(:) = s(:,1,isym)*xkk_cryst(1) + &
                       s(:,2,isym)*xkk_cryst(2) + &
                       s(:,3,isym)*xkk_cryst(3) - sxk(:)
              dxk(:) = dxk(:) - NINT(dxk(:))
              IF ( .NOT. ( ABS(dxk(1)) <= eps .AND. &
                           ABS(dxk(2)) <= eps .AND. &
                           ABS(dxk(3)) <= eps )   ) THEN
                  WRITE(*,*) ik,iq
                  WRITE(*,*) ikq,ikk,isym
                  WRITE(*,*) dxk(:)
                  CALL errore( 'exx_grid_check', 'something wrong', 1 )
              ENDIF
              !
          ENDDO
        ENDDO
      ENDDO
    ENDDO
    !
    RETURN
    !
    !
  END SUBROUTINE exx_grid_check
  !
  !
  !-----------------------------------------------------------------------
  SUBROUTINE exx_set_symm( nr1, nr2, nr3, nr1x, nr2x, nr3x )
    !-----------------------------------------------------------------------
    !! Uses \(\text{nkqs}\) and \(\text{index_sym}\) from module \(\texttt{exx}\),
    !! computes \(\text{rir}\).
    !
    USE symm_base,  ONLY : nsym, s, ft
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN) :: nr1, nr2, nr3, nr1x, nr2x, nr3x 
    !
    ! ... local variables
    !
    INTEGER :: ikq, isym, i,j,k, ri,rj,rk, ir, nxxs
    INTEGER, allocatable :: ftau(:,:), s_scaled(:,:,:)
    !
    nxxs = nr1x*nr2x*nr3x
    !
    IF (.NOT. ALLOCATED(rir)) THEN   
        ALLOCATE( rir(nxxs,nsym) )
    ELSEIF ((SIZE(rir,1) /= nxxs) ) THEN 
        DEALLOCATE( rir )
        ALLOCATE( rir(nxxs,nsym) )
    ENDIF
    !
    rir = 0
    ALLOCATE ( ftau(3,nsym), s_scaled(3,3,nsym) )
    CALL scale_sym_ops (nsym, s, ft, nr1, nr2, nr3, s_scaled, ftau)
    DO isym = 1, nsym
       DO k = 1, nr3
          DO j = 1, nr2
             DO i = 1, nr1
                CALL rotate_grid_point( s_scaled(1,1,isym), ftau(1,isym), &
                     i, j, k, nr1, nr2, nr3, ri, rj, rk )
                ir = i + (j-1)*nr1x + (k-1)*nr1x*nr2x
                rir(ir,isym) = ri + (rj-1)*nr1x + (rk-1)*nr1x*nr2x
             ENDDO
          ENDDO
       ENDDO
    ENDDO
    !
    DEALLOCATE ( s_scaled, ftau )
    !
  END SUBROUTINE exx_set_symm
  !
  !
  !-----------------------------------------------------------------------
  SUBROUTINE g2_convolution_all( ngm, g, xk, xkq, iq, current_k )
    !-----------------------------------------------------------------------
    !! Wrapper for g2_convolution.
    !
    USE kinds,     ONLY : DP
    USE klist,     ONLY : nks
    !
    IMPLICIT NONE
    !
    INTEGER,  INTENT(IN) :: ngm
    !! Number of G vectors
    REAL(DP), INTENT(IN) :: g(3,ngm)
    !! Cartesian components of G vectors
    REAL(DP), INTENT(IN) :: xk(3)
    !! current k vector
    REAL(DP), INTENT(IN) :: xkq(3)
    !! current q vector
    INTEGER, INTENT(IN) :: current_k
    !! current k-point index
    INTEGER, INTENT(IN) :: iq
    !! q-grid point index
    !
    ! ... Check if coulomb_fac has been allocated
    IF( .NOT. ALLOCATED( coulomb_fac ) ) ALLOCATE( coulomb_fac(ngm,nqs,nks) )
    !
    ! ... Check if coulomb_done has been allocated
    IF( .NOT. ALLOCATED( coulomb_done) ) THEN
       ALLOCATE( coulomb_done(nqs,nks) )
       coulomb_done = .FALSE.
    ENDIF
    !
    ! ... return if this k and k' already computed, otherwise compute it
    IF ( coulomb_done(iq,current_k) ) RETURN
    !
    CALL g2_convolution( ngm, g, xk, xkq, coulomb_fac(:,iq,current_k) )
    !
    coulomb_done(iq,current_k) = .TRUE.
    !
  END SUBROUTINE g2_convolution_all
  !
  !
  !-----------------------------------------------------------------------
  SUBROUTINE g2_convolution( ngm, g, xk, xkq, fac )
    !-----------------------------------------------------------------------
    !! This routine calculates the 1/|r-r'| part of the exact exchange
    !! expression in reciprocal space (the G^-2 factor).
    !! It then regularizes it according to the specified recipe.
    !
    USE kinds,      ONLY : DP
    USE cell_base,  ONLY : tpiba, at, tpiba2
    USE constants,  ONLY : fpi, e2, pi
    !
    IMPLICIT NONE
    !
    INTEGER,  INTENT(IN) :: ngm
    !! Number of G vectors
    REAL(DP), INTENT(IN) :: g(3,ngm)
    !! Cartesian components of G vectors
    REAL(DP), INTENT(IN) :: xk(3)
    !! current k vector
    REAL(DP), INTENT(IN) :: xkq(3)
    !! current q vector
    REAL(DP), INTENT(INOUT) :: fac(ngm)
    !! Calculated convolution
    !
    ! ... local variables
    !
    INTEGER :: ig !Counters
    REAL(DP) :: q(3), qq, x
    REAL(DP) :: grid_factor_track(ngm), qq_track(ngm)
    REAL(DP) :: nqhalf_dble(3)
    LOGICAL :: odg(3)
    !
    ! ... First the types of Coulomb potential that need q(3) and an external call
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
    nqhalf_dble(1:3) = (/ DBLE(nq1)*0.5_DP, DBLE(nq2)*0.5_DP, DBLE(nq3)*0.5_DP /)
    !
    ! Set the grid_factor_track and qq_track
    !
    IF ( x_gamma_extrapolation ) THEN
!$omp parallel do default(shared), private(ig,q,x,odg)
       DO ig = 1, ngm
          q(:)= xk(:) - xkq(:) + g(:,ig)
          qq_track(ig) = SUM(q(:)**2) * tpiba2
          x = (q(1)*at(1,1)+q(2)*at(2,1)+q(3)*at(3,1))*nqhalf_dble(1)
          odg(1) = ABS(x-NINT(x)) < eps
          x = (q(1)*at(1,2)+q(2)*at(2,2)+q(3)*at(3,2))*nqhalf_dble(2)
          odg(2) = ABS(x-NINT(x)) < eps
          x = (q(1)*at(1,3)+q(2)*at(2,3)+q(3)*at(3,3))*nqhalf_dble(3)
          odg(3) = ABS(x-NINT(x)) < eps
          IF( ALL( odg(:) ) ) THEN
             grid_factor_track(ig) = 0._DP ! on double grid
          ELSE
             grid_factor_track(ig) = grid_factor ! not on double grid
          ENDIF
       ENDDO
!$omp end parallel do
    ELSE
!$omp parallel do default(shared), private(ig,q)
       DO ig = 1, ngm
          q(:) = xk(:) - xkq(:) + g(:,ig)
          qq_track(ig) = SUM(q(:)**2) * tpiba2
       ENDDO
!$omp end parallel do
       grid_factor_track = 1._DP
    ENDIF
    !
    ! ... The big loop
    !
!$omp parallel do default(shared), private(ig,qq)
    DO ig = 1, ngm
      !
      qq = qq_track(ig)
      !
      IF(gau_scrlen > 0) THEN
         fac(ig) = e2*((pi/gau_scrlen)**(1.5_DP))*EXP(-qq/4._DP/gau_scrlen) * grid_factor_track(ig)
         !
      ELSEIF (qq > eps_qdiv) THEN
         !
         IF ( erfc_scrlen > 0  ) THEN
            fac(ig) = e2*fpi/qq*(1._DP-EXP(-qq/4._DP/erfc_scrlen**2)) * grid_factor_track(ig)
         ELSEIF( erf_scrlen > 0 ) THEN
            fac(ig) = e2*fpi/qq*(EXP(-qq/4._DP/erf_scrlen**2)) * grid_factor_track(ig)
         ELSE
            fac(ig) = e2*fpi/( qq + yukawa ) * grid_factor_track(ig) ! as HARTREE
         ENDIF
         !
      ELSE
         !
         fac(ig) = - exxdiv ! or rather something ELSE (see F.Gygi)
         !
         IF (yukawa>0._DP .AND. .NOT.x_gamma_extrapolation) fac(ig) = fac(ig) + &
                                                            e2*fpi/( qq + yukawa )
         !
         IF (erfc_scrlen>0._DP .AND. .NOT.x_gamma_extrapolation) fac(ig) = fac(ig) + &
                                                                  e2*pi/(erfc_scrlen**2)
         !
      ENDIF
      !
    ENDDO
!$omp end parallel do
    !
  END SUBROUTINE g2_convolution
  !
  !
  !-----------------------------------------------------------------------
  FUNCTION exx_divergence()
     !-----------------------------------------------------------------------
     !
     USE constants,      ONLY : fpi, e2, pi
     USE cell_base,      ONLY : bg, at, alat, omega
     USE gvect,          ONLY : ngm, g
     USE gvecw,          ONLY : gcutw
     USE mp_exx,         ONLY : intra_egrp_comm
     USE mp,             ONLY : mp_sum
     !
     IMPLICIT NONE
     !
     REAL(DP) :: exx_divergence
     !
     ! ... local variables
     !
     INTEGER :: iq1,iq2,iq3, ig
     REAL(DP) :: div, dq1, dq2, dq3, xq(3), q_, qq, &
                 tpiba2, alpha, x, q(3)
     INTEGER :: nqq, iq
     REAL(DP) :: aa, dq
     !
     CALL start_clock( 'exx_div' )
     !
     tpiba2 = (fpi / 2.d0 / alat)**2
     !
     alpha  = 10._dp / gcutw
     !
     IF ( .NOT. use_regularization ) THEN
        exx_divergence = 0._dp
        RETURN
     ENDIF
     !
     dq1 = 1._dp / DBLE(nq1)
     dq2 = 1._dp / DBLE(nq2) 
     dq3 = 1._dp / DBLE(nq3) 
     !
     div = 0._dp
     !
     DO iq1 = 1, nq1
        DO iq2 = 1, nq2
           DO iq3 = 1, nq3
              xq(:) = bg(:,1) * (iq1-1) * dq1 + &
                      bg(:,2) * (iq2-1) * dq2 + &
                      bg(:,3) * (iq3-1) * dq3
              !
              DO ig = 1, ngm
                 !
                 q(1) = xq(1) + g(1,ig)
                 q(2) = xq(2) + g(2,ig)
                 q(3) = xq(3) + g(3,ig)
                 qq = ( q(1)**2 + q(2)**2 + q(3)**2 )
                 !
                 IF (x_gamma_extrapolation) THEN
                    on_double_grid = .TRUE.
                    x = 0.5d0*(q(1)*at(1,1)+q(2)*at(2,1)+q(3)*at(3,1))*nq1
                    on_double_grid = on_double_grid .AND. (ABS(x-NINT(x))<eps)
                    x = 0.5d0*(q(1)*at(1,2)+q(2)*at(2,2)+q(3)*at(3,2))*nq2
                    on_double_grid = on_double_grid .AND. (ABS(x-NINT(x))<eps)
                    x = 0.5d0*(q(1)*at(1,3)+q(2)*at(2,3)+q(3)*at(3,3))*nq3
                    on_double_grid = on_double_grid .AND. (ABS(x-NINT(x))<eps)
                 ENDIF
                 !
                 IF (.NOT.on_double_grid) THEN
                    IF ( qq > 1.d-8 ) THEN
                       IF ( erfc_scrlen > 0 ) THEN
                          div = div + EXP( -alpha * qq) / qq * &
                                (1._dp-EXP(-qq*tpiba2/4.d0/erfc_scrlen**2)) * grid_factor
                       ELSEIF ( erf_scrlen >0 ) THEN
                          div = div + EXP( -alpha * qq) / qq * &
                                (EXP(-qq*tpiba2/4.d0/erf_scrlen**2)) * grid_factor
                       ELSE

                          div = div + EXP( -alpha * qq) / (qq + yukawa/tpiba2) &
                                                     * grid_factor
                       ENDIF
                    ENDIF
                 ENDIF
                 !
              ENDDO
              !
           ENDDO
        ENDDO
     ENDDO
     !
     CALL mp_sum( div, intra_egrp_comm )
     !
     IF (gamma_only) THEN
        div = 2.d0 * div
     ENDIF
     !
     IF ( .NOT. x_gamma_extrapolation ) THEN
        IF ( yukawa > 0._dp) THEN
           div = div + tpiba2/yukawa
        ELSEIF( erfc_scrlen > 0._dp ) THEN
           div = div + tpiba2/4.d0/erfc_scrlen**2
        ELSE
           div = div - alpha
        ENDIF
     ENDIF
     !
     div = div * e2 * fpi / tpiba2 / nqs
     !
     alpha = alpha / tpiba2
     !
     nqq = 100000
     dq = 5.0d0 / SQRT(alpha) / nqq
     aa = 0._dp
     !
     DO iq = 0, nqq
        q_ = dq * (iq+0.5d0)
        qq = q_ * q_
        IF ( erfc_scrlen > 0 ) THEN
           aa = aa  -EXP( -alpha * qq) * EXP(-qq/4.d0/erfc_scrlen**2)*dq
        ELSEIF ( erf_scrlen > 0 ) THEN
           aa = 0._dp
        ELSE
           aa = aa - EXP( -alpha * qq) * yukawa / (qq + yukawa)*dq
        ENDIF
     ENDDO
     !
     aa = aa * 8.d0/fpi
     aa = aa + 1._dp/SQRT(alpha*0.25d0*fpi)
     IF ( erf_scrlen > 0) aa = 1._dp/SQRT((alpha+1._dp/4.d0/erf_scrlen**2)*0.25d0*fpi)
     div = div - e2*omega * aa
     !
     exx_divergence = div * nqs
     !
     CALL stop_clock( 'exx_div' )
     !
     RETURN
     !
  END FUNCTION exx_divergence
  !
  !
END MODULE exx_base
