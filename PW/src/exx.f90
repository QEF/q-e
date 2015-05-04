!
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
  USE kinds,                ONLY : DP
  USE coulomb_vcut_module,  ONLY : vcut_init, vcut_type, vcut_info, &
                                   vcut_get,  vcut_spheric_get
  USE noncollin_module,     ONLY : noncolin, npol
  USE io_global,            ONLY : ionode
  USE fft_custom,           ONLY : fft_cus
  !
  USE control_flags, ONLY : tqr

  IMPLICIT NONE
  SAVE
  !
  ! general purpose vars
  !
  REAL(DP):: exxalfa=0._dp                ! 1 if exx, 0 elsewhere
  INTEGER :: exx_nwordwfc, ji
  CHARACTER(len=1) :: exx_augmented = 'x' ! r -> real space
                                          ! k -> reciprocal space
                                          ! x -> do not augment
  !
  ! variables defining the auxiliary k-point grid 
  ! used in X BZ integration
  !
  INTEGER :: nq1=1, nq2=1, nq3=1         ! integers defining the X integration mesh
  INTEGER :: nqs=1                       ! number of points in the q-gridd
  INTEGER :: nkqs                        ! total number of different k+q
  !
  REAL(DP),    ALLOCATABLE :: xkq_collect(:,:)  ! xkq(3,nkqs) the auxiliary k+q set
  REAL(DP),    ALLOCATABLE :: x_occupation(:,:)           
                                         ! x_occupation(nbnd,nks) the weight of 
                                         ! auxiliary functions in the density matrix
  COMPLEX(DP), ALLOCATABLE :: exxbuff(:,:,:)
                                         ! temporary buffer for wfc storage
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
!
!  Used for k points pool parallelization. All pools need these quantities.
!  They are allocated only IF needed.
!
  REAL(DP),    ALLOCATABLE :: xk_collect(:,:)
  REAL(DP),    ALLOCATABLE :: wk_collect(:)
  REAL(DP),    ALLOCATABLE :: wg_collect(:,:)
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
  CHARACTER(32)    :: exxdiv_treatment  = ''
  !
  ! x_gamma_extrapolation
  LOGICAL           :: x_gamma_extrapolation =.TRUE.
  LOGICAl           :: on_double_grid =.FALSE.
  REAL(DP)          :: grid_factor = 1.d0 !8.d0/7.d0 
  !
  ! Gygi-Baldereschi 
  LOGICAL           :: use_regularization = .TRUE.
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
  LOGICAL           :: use_coulomb_vcut_ws = .FALSE.
  LOGICAL           :: use_coulomb_vcut_spheric = .FALSE.
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
  TYPE(fft_cus) exx_fft_g2r     ! Grid for wfcs -> real space
  TYPE(fft_cus) exx_fft_r2g     ! Grid for real space -> restricted G space
  REAL(DP)  :: ecutfock         ! energy cutoff for custom grid
  REAL(DP)  :: exx_dual = 4.0_DP! dual for the custom grid
 CONTAINS
#define _CX(A)  CMPLX(A,0._dp,kind=DP)
#define _CY(A)  CMPLX(0._dp,-A,kind=DP)
  !------------------------------------------------------------------------
  SUBROUTINE exx_grid_convert( psi, npw, fft, psi_t, sign, igkt )
    !------------------------------------------------------------------------
    ! 
    ! This routine reorders the gvectors of the wavefunction psi and
    ! puts the result in psi_t. This reordering is needed when going
    ! between two different fft grids.
    !
    ! sign > 0 goes from the smooth grid to the grid defined in fft
    ! sign < 0 goes from the grid defined in fft to the smooth grid 
    !

    USE mp_bands,   ONLY : me_bgrp, nproc_bgrp, intra_bgrp_comm
    USE fft_custom, ONLY : reorderwfp_col
    USE gvect,      ONLY : ig_l2g

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: npw
    COMPLEX(kind=DP), INTENT(IN) :: psi(npw)
    COMPLEX(kind=DP), INTENT(INOUT) :: psi_t(:)
    INTEGER, OPTIONAL, INTENT(INOUT) :: igkt(:)
    INTEGER, INTENT(IN) :: sign
    TYPE(fft_cus), INTENT(IN) :: fft
    
    INTEGER :: ig
    
    CALL start_clock('exx_grid_convert')
    
    IF(sign > 0 .AND. PRESENT(igkt) ) THEN
       DO ig=1, fft%ngmt
          igkt(ig)=ig
       ENDDO
    ENDIF

    
    IF( fft%dual_t==4.d0) THEN
       psi_t(1:fft%npwt)=psi(1:fft%npwt)
    ELSE
       IF (sign > 0 ) THEN
          CALL reorderwfp_col ( 1, npw, fft%npwt, psi, psi_t, npw, fft%npwt,&
               & ig_l2g, fft%ig_l2gt, fft%ngmt_g, me_bgrp, nproc_bgrp,&
               & intra_bgrp_comm )
       ELSE
          CALL reorderwfp_col ( 1, fft%npwt, npw, psi, psi_t, fft%npwt, npw,&
               & fft%ig_l2gt, ig_l2g, fft%ngmt_g, me_bgrp, nproc_bgrp,&
               & intra_bgrp_comm )
       ENDIF
    ENDIF

    CALL stop_clock('exx_grid_convert')

    RETURN
    
  END SUBROUTINE exx_grid_convert
  !------------------------------------------------------------------------
  !
  !------------------------------------------------------------------------
  SUBROUTINE exx_fft_create ()
    USE wvfct,        ONLY : ecutwfc, npw
    USE gvect,        ONLY : ecutrho, ig_l2g
    USE uspp,         ONLY : okvan
    USE paw_variables,ONLY : okpaw
    USE control_flags,ONLY : gamma_only

    IMPLICIT NONE

    ! Initalise the g2r grid that allows us to put the wavefunction
    ! onto the new (smaller) grid for rho.
    exx_fft_g2r%ecutt=ecutwfc
    exx_fft_g2r%dual_t=ecutfock/ecutwfc
    CALL allocate_fft_custom(exx_fft_g2r)

    IF (MAXVAL( ABS(ig_l2g(1:npw)-exx_fft_g2r%ig_l2gt(1:npw))) /= 0) THEN
       CALL errore('exx_fft_create', ' exx fft grid not compatible with &
            &the smooth fft grid. ', 1 )
       
    ENDIF

    ! Initalise the r2g grid that we then use when applying the Fock
    ! operator in our new restricted space.
    exx_fft_r2g%ecutt=ecutfock/exx_dual
    exx_fft_r2g%dual_t=exx_dual
    CALL allocate_fft_custom(exx_fft_r2g)
    !------------------------------------------------------------------------
  END SUBROUTINE exx_fft_create
  !------------------------------------------------------------------------
  !
  !------------------------------------------------------------------------
  SUBROUTINE exx_fft_destroy ()
  !------------------------------------------------------------------------
    USE fft_custom,  ONLY : deallocate_fft_custom

    IMPLICIT NONE

    CALL deallocate_fft_custom(exx_fft_g2r)
    CALL deallocate_fft_custom(exx_fft_r2g)
    !------------------------------------------------------------------------
  END SUBROUTINE exx_fft_destroy
  !------------------------------------------------------------------------
  !
  !------------------------------------------------------------------------
  SUBROUTINE deallocate_exx ()
    !------------------------------------------------------------------------
    !
    USE becmod, ONLY : deallocate_bec_type, is_allocated_bec_type, bec_type
    USE us_exx, ONLY : becxx
    IMPLICIT NONE
    INTEGER :: ikq
    !
    IF ( allocated(index_xkq) ) DEALLOCATE(index_xkq)
    IF ( allocated(index_xk ) ) DEALLOCATE(index_xk )
    IF ( allocated(index_sym) ) DEALLOCATE(index_sym)
    IF ( ALLOCATED (rir)       ) DEALLOCATE (rir)
    IF ( allocated(x_occupation) ) DEALLOCATE(x_occupation)
    IF ( allocated(xkq_collect) )  DEALLOCATE(xkq_collect)
    IF ( allocated(exxbuff) )      DEALLOCATE(exxbuff)
    !
    IF(ALLOCATED(becxx)) THEN
      DO ikq = 1, nkqs
        IF(is_allocated_bec_type(becxx(ikq))) CALL deallocate_bec_type(becxx(ikq))
      ENDDO
      DEALLOCATE(becxx)
   ENDIF
    !
    CALL exx_fft_destroy()
    !
    !  Pool variables deallocation
    !
    IF ( allocated (xk_collect) )  DEALLOCATE( xk_collect )
    IF ( allocated (wk_collect) )  DEALLOCATE( wk_collect )
    IF ( allocated (wg_collect) )  DEALLOCATE( wg_collect )
    !
    !
    !------------------------------------------------------------------------
  END SUBROUTINE deallocate_exx
  !------------------------------------------------------------------------
  !
  SUBROUTINE exx_grid_reinit()
    IMPLICIT NONE
    DEALLOCATE(xkq_collect,index_xk,index_sym)
    exx_grid_initialized = .false.
    CALL exx_grid_init()
  END SUBROUTINE exx_grid_reinit
  !------------------------------------------------------------------------
  SUBROUTINE exx_grid_init()
    !------------------------------------------------------------------------
    !
    USE symm_base,  ONLY : nsym, s
    USE cell_base,  ONLY : bg, at
    USE spin_orb,   ONLY : domag
    USE noncollin_module, ONLY : nspin_lsda
    USE klist,      ONLY : xk, wk, nkstot, nks
    USE wvfct,      ONLY : nbnd
    USE io_global,  ONLY : stdout
    USE start_k,    ONLY : nk1,nk2,nk3
    USE mp_pools,   ONLY : npool
    !
    IMPLICIT NONE
    !
    CHARACTER(13) :: sub_name='exx_grid_init'
    INTEGER       :: iq1, iq2, iq3, isym, ik, ikq, iq, max_nk, temp_nkqs
    INTEGER, allocatable :: temp_index_xk(:), temp_index_sym(:)
    INTEGER, allocatable :: temp_index_ikq(:), new_ikq(:)
    REAL(DP),allocatable :: temp_xkq(:,:)
    LOGICAL      :: xk_not_found
    REAL(DP)     :: sxk(3), dxk(3), xk_cryst(3)
    REAL(DP)     :: dq1, dq2, dq3
    CHARACTER (LEN=6), EXTERNAL :: int_to_char
    !
    CALL start_clock ('exx_grid')
    !
    IF(nq1<=0) nq1 = nk1
    IF(nq2<=0) nq2 = nk2
    IF(nq3<=0) nq3 = nk3
    IF(nkstot==nspin_lsda) THEN
      nq1=1; nq2=1; nq3=1
    ENDIF
     
    IF(ANY((/nq1,nq2,nq3/)<=0)) CALL errore('exx_grid_init',"wrong EXX q grid", 1)
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
    ! all processors need to have access to all k+q points
    !
    IF ( .NOT.allocated (xk_collect) )  ALLOCATE(xk_collect(3,nkstot))
    IF ( .NOT.allocated (wk_collect) )  ALLOCATE(wk_collect(nkstot))
    ! the next if/then if probably not necessary, as xk_wk collect can
    ! deal with npool==1, leaving it for clarity.
    IF ( npool > 1 ) THEN
      CALL xk_wk_collect(xk_collect, wk_collect, xk, wk, nkstot, nks)
    ELSE
      xk_collect(:,1:nks) = xk(:,1:nks)
      wk_collect(1:nks) = wk(1:nks)
    ENDIF
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
          ! *** do-loop skipped the first time because temp_nksq == 0
          DO ikq=1, temp_nkqs
            IF (xk_not_found ) THEN
                dxk(:) = sxk(:)-temp_xkq(:,ikq) - nint(sxk(:)-temp_xkq(:,ikq))
                IF ( abs(dxk(1)).le.eps .and. &
                     abs(dxk(2)).le.eps .and. &
                     abs(dxk(3)).le.eps ) xk_not_found = .false.
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
                IF ( abs(dxk(1)).le.eps .and. &
                     abs(dxk(2)).le.eps .and. &
                     abs(dxk(3)).le.eps ) xk_not_found = .false.
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
    dq1= 1._dp/DBLE(nq1)
    dq2= 1._dp/DBLE(nq2)
    dq3= 1._dp/DBLE(nq3)
    !
    ! allocate and fill the array index_xkq(nkstot,nqs)
    !
    if(.not.allocated(index_xkq))    ALLOCATE( index_xkq(nkstot,nqs) )
    if(.not.allocated(x_occupation)) ALLOCATE( x_occupation(nbnd,nkstot) )
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
                    IF ( ALL(abs(dxk) < eps ) ) THEN
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
                write (*,*) ik, iq, temp_nkqs
                write (*,*) sxk(:)
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
      WRITE(stdout, '(5x,3a)') "EXX: setup a grid of "//TRIM(int_to_char(nkqs))&
                           //" q-points centered on each k-point"
      WRITE( stdout, '(5x,a)' ) '(k+q)-points:'
      do ik = 1, nkqs
          WRITE( stdout, '(3f12.7,5x,i2,i5)') (xkq_collect (ikq, ik) , ikq = 1, 3) , &
                 index_xk(ik), index_sym(ik)
      enddo
    ELSE
      WRITE(stdout, '("EXX: grid of k+q points same as grid of k-points")')
    ENDIF
    
    ! if nspin == 2, the kpoints are repeated in couples (spin up, spin down)
    IF (nspin_lsda == 2) THEN
      DO ik = 1, nkstot/2
          DO iq =1, nqs
            index_xkq(nkstot/2+ik,iq) = index_xkq(ik,iq) + nkqs
          END DO
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
    CALL exx_grid_check () 
    !
    CALL exx_set_qnorm(nkqs, xkq_collect)
    !
    CALL stop_clock ('exx_grid')
    !
    RETURN
    !------------------------------------------------------------------------
  END SUBROUTINE exx_grid_init
  !------------------------------------------------------------------------
  !
  !-----------------------------------------------------------------------
  SUBROUTINE exx_n_plane_waves(ecutwfc, tpiba2, g, ngm, npwx)
    !-----------------------------------------------------------------------
    !
    ! Find number of plane waves for each k-point, keeping in mind that EXX uses
    ! a larger non-reduced grid of k and k+q points
    USE kinds, ONLY : DP
    !USE exx,              ONLY : nkqs, xkq_collect,exx_grid_initialized
    USE funct, ONLY : dft_is_hybrid
    USE uspp,  ONLY : okvan
    IMPLICIT NONE
    !
    integer, intent(in) :: ngm
    real(DP),intent(in) :: ecutwfc, tpiba2, g (3, ngm)
    integer, intent(out):: npwx
    integer,allocatable :: ngkq(:)
    IF(.not. okvan) RETURN
    IF(.not.dft_is_hybrid()) RETURN
    IF(.not.exx_grid_initialized) &
      CALL errore("exx_n_plane_waves","you must initialize the grid first",1)
    ALLOCATE(ngkq(nkqs))
    CALL n_plane_waves (ecutwfc, tpiba2, nkqs, xkq_collect, g, ngm, npwx, ngkq)
    DEALLOCATE(ngkq)
    RETURN
    !------------------------------------------------------------------------
  END SUBROUTINE exx_n_plane_waves
  !------------------------------------------------------------------------
  !
  !------------------------------------------------------------------------
  SUBROUTINE exx_set_qnorm(nq, xkq)
    !------------------------------------------------------------------------
    !
    ! WARNING: setting qnorm increases the amount of space allocated in allocate_nlpot
    !           DOING IT HERE AND TO USE QNORM FOR THIS IS WRONG because:
    ! 1. phonon may overwrite it later with something smaller, 
    !    the actual qnorm should include the rotated(k+q_phonon) in it
    ! 2. the value we use is sufficient, but may be too large, not clear to me
    USE kinds,      ONLY : DP
    USE klist,      ONLY : qnorm
    USE uspp,       ONLY : okvan
    USE klist,      ONLY : xk, nks
    !
    IMPLICIT NONE
    INTEGER,INTENT(in)  :: nq
    REAL(DP),INTENT(in) :: xkq(3,nq)
    INTEGER :: i,j
    !
    IF(.not.okvan) RETURN
!     DO i = 1,nq
!       qnorm = MAX(qnorm,  SQRT(SUM(xkq(:,i)**2))  )
!     ENDDO

    DO i = 1,nq
    DO j = 1,nks
      qnorm = MAX(qnorm, SQRT( SUM((xk(:,j)-xkq(:,i))**2) ))
    ENDDO
    ENDDO
  
    RETURN
    !------------------------------------------------------------------------
  END SUBROUTINE exx_set_qnorm
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
    SELECT CASE ( TRIM(exxdiv_treatment) ) 
    CASE ( "gygi-baldereschi", "gygi-bald", "g-b", "gb" )
      !
      use_regularization = .TRUE.
      !
      !
    CASE ( "vcut_ws" )
      !
      use_coulomb_vcut_ws = .TRUE.
      IF ( x_gamma_extrapolation ) &
            CALL errore(sub_name,'cannot USE x_gamm_extrap and vcut_ws', 1)
      !
    CASE ( "vcut_spherical" ) 
      !
      use_coulomb_vcut_spheric = .TRUE.
      IF ( x_gamma_extrapolation ) &
            CALL errore(sub_name,'cannot USE x_gamm_extrap and vcut_spherical', 1)
      !
    CASE ( "none" )
      use_regularization = .FALSE.
      !
    CASE DEFAULT
      CALL errore(sub_name,'invalid exxdiv_treatment: '//TRIM(exxdiv_treatment), 1)
    END SELECT
    !
    ! Set variables for Coulomb vcut
    ! NOTE: some memory is allocated inside this routine (in the var vcut)
    !       and should be deallocated somewehre, at the end of the run
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
        IF ( ionode ) CALL vcut_info( stdout, vcut )
        !          
    ENDIF
    RETURN
  !------------------------------------------------------------------------
  END SUBROUTINE exx_div_check 
  !------------------------------------------------------------------------


  !------------------------------------------------------------------------
  SUBROUTINE exx_grid_check ( )
    !------------------------------------------------------------------------
    USE symm_base,  ONLY : s
    USE cell_base,  ONLY : at
    USE klist,      ONLY : nkstot, xk
    USE mp_pools,   ONLY : npool
    IMPLICIT NONE
    REAL(DP) :: sxk(3), dxk(3), xk_cryst(3), xkk_cryst(3)
    INTEGER :: iq1, iq2, iq3, isym, ik, ikk, ikq, iq
    REAL(DP) :: dq1, dq2, dq3
    dq1= 1._dp/DBLE(nq1)
    dq2= 1._dp/DBLE(nq2)
    dq3= 1._dp/DBLE(nq3)

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

              IF (npool>1) THEN
                xkk_cryst(:) = at(1,:)*xk_collect(1,ikk)+at(2,:)*xk_collect(2,ikk)+at(3,:)*xk_collect(3,ikk)
              ELSE
                xkk_cryst(:) = at(1,:)*xk(1,ikk)+at(2,:)*xk(2,ikk)+at(3,:)*xk(3,ikk)
              ENDIF

              IF (isym < 0 ) xkk_cryst(:) = - xkk_cryst(:)
              isym = abs (isym)
              dxk(:) = s(:,1,isym)*xkk_cryst(1) + &
                       s(:,2,isym)*xkk_cryst(2) + &
                       s(:,3,isym)*xkk_cryst(3) - sxk(:)
              dxk(:) = dxk(:) - nint(dxk(:))
              IF ( .not. ( abs(dxk(1)).le.eps .and. &
                           abs(dxk(2)).le.eps .and. &
                           abs(dxk(3)).le.eps )   ) THEN
                  write(*,*) ik,iq
                  write(*,*) ikq,ikk,isym
                  write(*,*) dxk(:)
                  CALL errore('exx_grid_check', 'something wrong', 1 )
              ENDIF

          ENDDO
        ENDDO
      ENDDO
    ENDDO
    !
    return

    !------------------------------------------------------------------------
  END SUBROUTINE exx_grid_check
  !------------------------------------------------------------------------
  !
  !------------------------------------------------------------------------
  SUBROUTINE exx_restart(l_exx_was_active)
     !------------------------------------------------------------------------
     !This SUBROUTINE is called when restarting an exx calculation
     USE funct,                ONLY : get_exx_fraction, start_exx, exx_is_active, &
                                     get_screening_parameter
     USE fft_base,             ONLY : dffts

     IMPLICIT NONE
     LOGICAL, INTENT(IN) :: l_exx_was_active

     IF (.not. l_exx_was_active ) return ! nothing had happpened yet
     !
     exx_nwordwfc=2*dffts%nnr
     erfc_scrlen = get_screening_parameter()
     exxdiv = exx_divergence() 
     exxalfa = get_exx_fraction()
     CALL start_exx
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
    USE wavefunctions_module, ONLY : evc  
    USE io_files,             ONLY : nwordwfc, iunwfc, iunigk
    USE buffers,              ONLY : get_buffer
    USE gvecs,                ONLY : nls
    USE wvfct,                ONLY : nbnd, npwx, npw, igk, wg
    USE control_flags,        ONLY : gamma_only
    USE klist,                ONLY : ngk, nks, nkstot
    USE symm_base,            ONLY : nsym, s, sr, ftau
    USE mp_pools,             ONLY : npool, nproc_pool, me_pool, inter_pool_comm
    USE mp_bands,             ONLY : nproc_bgrp, me_bgrp, init_index_over_band,&
                                     inter_bgrp_comm, ibnd_start, ibnd_end,nbgrp
    USE mp,                   ONLY : mp_sum
    USE funct,                ONLY : get_exx_fraction, start_exx,exx_is_active,&
                                     get_screening_parameter, get_gau_parameter
    USE fft_base,             ONLY : cgather_smooth, cscatter_smooth,&
                                     dffts, cgather_custom, cscatter_custom
    USE fft_interfaces,       ONLY : invfft
    USE becmod,               ONLY : allocate_bec_type, bec_type
    USE uspp,                 ONLY : nkb, okvan
    USE us_exx,               ONLY : becxx
    USE paw_variables,        ONLY : okpaw
    USE paw_exx,              ONLY : PAW_init_keeq

    IMPLICIT NONE
    INTEGER :: ik,ibnd, i, j, k, ir, ri, rj, rk, isym, ikq, ig
    INTEGER :: h_ibnd
    INTEGER :: ibnd_loop_start, ibnd_buff_start, ibnd_buff_end
    INTEGER :: ipol, jpol
    COMPLEX(DP),ALLOCATABLE :: temppsic(:),      psic(:), tempevc(:,:)
    COMPLEX(DP),ALLOCATABLE :: temppsic_nc(:,:), psic_nc(:,:)
    INTEGER :: nxxs, nrxxs, nr1x,nr2x,nr3x,nr1,nr2,nr3
#ifdef __MPI
    COMPLEX(DP),allocatable  :: temppsic_all(:),      psic_all(:)
    COMPLEX(DP), ALLOCATABLE :: temppsic_all_nc(:,:), psic_all_nc(:,:)
#endif
    COMPLEX(DP) :: d_spin(2,2,48)
    INTEGER :: current_ik
    logical, allocatable :: ispresent(:)
    integer       :: find_current_k

    CALL start_clock ('exxinit')
    !
    !  prepare the symmetry matrices for the spin part
    !
    IF (noncolin) THEN
       DO isym=1,nsym
          CALL find_u(sr(:,:,isym), d_spin(:,:,isym))
       ENDDO
    ENDIF

    ! Beware: not the same as nrxxs in parallel case
    IF(gamma_only) THEN
       CALL exx_fft_create()
       nxxs =exx_fft_g2r%dfftt%nr1x *exx_fft_g2r%dfftt%nr2x *exx_fft_g2r%dfftt%nr3x 
       nrxxs= exx_fft_g2r%dfftt%nnr
       nr1  = exx_fft_g2r%dfftt%nr1
       nr2  = exx_fft_g2r%dfftt%nr2
       nr3  = exx_fft_g2r%dfftt%nr3
       nr1x = exx_fft_g2r%dfftt%nr1x
       nr2x = exx_fft_g2r%dfftt%nr2x
       nr3x = exx_fft_g2r%dfftt%nr3x
    ELSE
       nxxs = dffts%nr1x * dffts%nr2x * dffts%nr3x
       nrxxs= dffts%nnr
       nr1  = dffts%nr1
       nr2  = dffts%nr2
       nr3  = dffts%nr3
       nr1x = dffts%nr1x
       nr2x = dffts%nr2x
       nr3x = dffts%nr3x
    ENDIF
#ifdef __MPI
    IF (noncolin) THEN
       ALLOCATE(psic_all_nc(nxxs,npol), temppsic_all_nc(nxxs,npol) )
    ELSE
       ALLOCATE(psic_all(nxxs), temppsic_all(nxxs) )
    ENDIF
#endif
    CALL init_index_over_band(inter_bgrp_comm,nbnd)
    IF (noncolin) THEN
       ALLOCATE(temppsic_nc(nrxxs, npol), psic_nc(nrxxs, npol))
    ELSE
       ALLOCATE(temppsic(nrxxs), psic(nrxxs))
    ENDIF
    !
    ! prepare space to keep the <beta_I|phi_j> scalar products (for ultrasoft/paw only)
    IF(.not. allocated(becxx) .and. okvan) THEN 
        ALLOCATE(becxx(nkqs))
        DO ikq = 1,nkqs
            CALL allocate_bec_type( nkb, nbnd, becxx(ikq))
        ENDDO
    ENDIF
    !
    IF ( gamma_only ) THEN
        ibnd_buff_start = ibnd_start/2
        IF(MOD(ibnd_start,2)==0) ibnd_buff_start = ibnd_buff_start -1
        !
        ibnd_buff_end = ibnd_end/2
        IF(MOD(ibnd_end,2)==1) ibnd_buff_end = ibnd_buff_end +1
    ELSE
        ibnd_buff_start = ibnd_start
        ibnd_buff_end   = ibnd_end
    ENDIF
    !
    IF (.NOT. allocated(exxbuff)) &
        ALLOCATE( exxbuff(nrxxs*npol, ibnd_buff_start:ibnd_buff_end, nkqs))
    !
    ALLOCATE(tempevc( npwx*npol, nbnd ))
    ALLOCATE(ispresent(nsym))
    IF(.NOT. ALLOCATED(rir)) ALLOCATE(rir(nxxs,nsym))
    rir = 0
    exx_nwordwfc=2*nrxxs
    IF (.not.exx_is_active()) THEN 
       !
       erfc_scrlen = get_screening_parameter()
       gau_scrlen = get_gau_parameter()
       exxdiv  = exx_divergence() 
       exxalfa = get_exx_fraction()
       !
       CALL start_exx()
    ENDIF

    IF (.NOT.allocated (wg_collect)) ALLOCATE(wg_collect(nbnd,nkstot))
    IF (npool>1) THEN
      CALL wg_all(wg_collect, wg, nkstot, nks)
    ELSE
      wg_collect = wg
    ENDIF

    IF ( nks > 1 ) REWIND( iunigk )

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
             CALL errore ('exxinit',' EXX smooth grid is not compatible with symmetry: &
                                    & change ecutfock',isym)
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

    exxbuff=(0.0_DP,0.0_DP)
    ! set appropriately the x_occupation
    DO ik =1,nkstot
       IF(ABS(wk_collect(ik)) > eps_occ ) THEN
          x_occupation(1:nbnd,ik) = wg_collect (1:nbnd, ik) / wk_collect(ik)
       ELSE
          x_occupation(1:nbnd,ik) = 0._dp
       ENDIF
    ENDDO
    !
    !   This is parallelized over pool. Each pool computes only its k-points
    !
    KPOINTS_LOOP : &
    DO ik = 1, nks
       npw = ngk (ik)
       IF ( nks > 1 ) THEN
          READ( iunigk ) igk
          CALL get_buffer(tempevc, nwordwfc, iunwfc, ik)
       ELSE
          tempevc(1:npwx*npol,1:nbnd) = evc(1:npwx*npol,1:nbnd)
       ENDIF
       !
       ! only useful for npool>1, but always work
       current_ik=find_current_k(ik, nkstot, nks)
       !
       IF_GAMMA_ONLY : & 
       IF (gamma_only) THEN
          !
          h_ibnd = ibnd_start/2
          !
          IF(MOD(ibnd_start,2)==0) THEN
             h_ibnd=h_ibnd-1
             ibnd_loop_start=ibnd_start-1
          ELSE
             ibnd_loop_start=ibnd_start
          ENDIF

          DO ibnd = ibnd_loop_start, ibnd_end, 2
             h_ibnd = h_ibnd + 1
             !
             temppsic(:) = ( 0._dp, 0._dp )
             !
             if ( ibnd < ibnd_end ) then
                DO ig=1,exx_fft_g2r%npwt
                   temppsic(exx_fft_g2r%nlt(ig))  = tempevc(ig,ibnd)  &
                        + ( 0._dp, 1._dp ) * tempevc(ig,ibnd+1)
                   temppsic(exx_fft_g2r%nltm(ig)) = CONJG( tempevc(ig,ibnd) ) &
                        + ( 0._dp, 1._dp ) * CONJG( tempevc(ig,ibnd+1) )
                END DO
             else
                DO ig=1,exx_fft_g2r%npwt
                   temppsic(exx_fft_g2r%nlt (ig)) = tempevc(ig,ibnd) 
                   temppsic(exx_fft_g2r%nltm(ig)) = CONJG( tempevc(ig,ibnd) ) 
                END DO
             end if

             CALL invfft ('CustomWave', temppsic, exx_fft_g2r%dfftt)

             DO ikq=1,nkqs
                IF (index_xk(ikq) .ne. current_ik) cycle
                isym = abs(index_sym(ikq) )
#ifdef __MPI
                CALL cgather_custom(temppsic,temppsic_all, exx_fft_g2r%dfftt)
                IF ( me_bgrp == 0 ) &
                psic_all(1:nxxs) = temppsic_all(rir(1:nxxs,isym))
                CALL cscatter_custom(psic_all,psic, exx_fft_g2r%dfftt)
#else
                psic(1:nrxxs) = temppsic(rir(1:nrxxs,isym))
#endif
                IF (index_sym(ikq) < 0 ) &
                   CALL errore('exxinit','index_sym < 0 with gamma_only (!?)',1)

                exxbuff(1:nrxxs,h_ibnd,ikq)=psic(1:nrxxs)
             ENDDO
          END DO
          !
       ELSE IF_GAMMA_ONLY 
          !
          IBND_LOOP_K : &
          DO ibnd = ibnd_start, ibnd_end
             !
             IF (noncolin) THEN
                temppsic_nc(:,:) = ( 0._dp, 0._dp )
                temppsic_nc(nls(igk(1:npw)),1) = tempevc(1:npw,ibnd)
                CALL invfft ('Wave', temppsic_nc(:,1), dffts)
                temppsic_nc(nls(igk(1:npw)),2) = tempevc(npwx+1:npwx+npw,ibnd)
                CALL invfft ('Wave', temppsic_nc(:,2), dffts)
             ELSE
                temppsic(:) = ( 0._dp, 0._dp )
                temppsic(nls(igk(1:npw))) = tempevc(1:npw,ibnd)
                CALL invfft ('Wave', temppsic, dffts)
             ENDIF
             !
             DO ikq=1,nkqs
                !
                IF (index_xk(ikq) /= current_ik) CYCLE
                isym = abs(index_sym(ikq) )
                !
                IF (noncolin) THEN ! noncolinear
#ifdef __MPI
                   DO ipol=1,npol
                      CALL cgather_smooth(temppsic_nc(:,ipol), temppsic_all_nc(:,ipol))
                   ENDDO
                   IF ( me_bgrp == 0 ) THEN
                      psic_all_nc(:,:) = (0.0_DP, 0.0_DP)
                      DO ipol=1,npol
                         DO jpol=1,npol
                            psic_all_nc(:,ipol)=psic_all_nc(:,ipol) &
                              +  CONJG(d_spin(jpol,ipol,isym))* &
                                 temppsic_all_nc(rir(:,isym),jpol)
                         ENDDO
                      ENDDO
                   ENDIF
                   DO ipol=1,npol
                      CALL cscatter_smooth(psic_all_nc(:,ipol), psic_nc(:,ipol))
                   ENDDO
#else
                   psic_nc(:,:) = (0._dp, 0._dp)
                   DO ipol=1,npol
                      DO jpol=1,npol
                         psic_nc(:,ipol) = psic_nc(:,ipol) + &
                              CONJG(d_spin(jpol,ipol,isym))* &
                                        temppsic_nc(rir(:,isym),jpol)
                      END DO
                   END DO
#endif
                   exxbuff(      1:  nrxxs,ibnd,ikq)=psic_nc(:,1)
                   exxbuff(nrxxs+1:2*nrxxs,ibnd,ikq)=psic_nc(:,2)
                ELSE ! noncolinear
#ifdef __MPI
                  CALL cgather_smooth(temppsic,temppsic_all)
                  IF ( me_bgrp == 0 ) &
                    psic_all(1:nxxs) = temppsic_all(rir(1:nxxs,isym))
                  CALL cscatter_smooth(psic_all,psic)
#else
                  psic(1:nrxxs) = temppsic(rir(1:nrxxs,isym))
#endif
                  IF (index_sym(ikq) < 0 ) psic(1:nrxxs) = CONJG(psic(1:nrxxs))
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
    !   All pools must have the complete set of wavefunctions (i.e. from every kpoint)
    IF (npool>1) CALL mp_sum(exxbuff, inter_pool_comm)
    !
    ! compute <beta_I|psi_j,k+q> for the entire de-symmetrized k+q grid
    !
    CALL compute_becxx()
    !
    ! CHECKME: probably it's enough that each pool computes its own bec
    !          and then I sum them like exxbuff, but check it. In this case this
    !          call should only act when index_xk(ikq) = current_ik
    !
    ! Initialize 4-wavefunctions one-center Fock integral \int \psi_a(r)\phi_a(r)\phi_b(r')\psi_b(r')/|r-r'|
    IF(okpaw) CALL PAW_init_keeq()
    !
    DEALLOCATE(tempevc)
    DEALLOCATE(ispresent)
    IF (noncolin) THEN
       DEALLOCATE(temppsic_nc, psic_nc)
#ifdef __MPI
       DEALLOCATE(temppsic_all_nc, psic_all_nc)
#endif 
    ELSE
       DEALLOCATE(temppsic, psic)
#ifdef __MPI
       DEALLOCATE(temppsic_all, psic_all)
#endif 
    ENDIF

    CALL stop_clock ('exxinit')  
    !
    !-----------------------------------------------------------------------
  END SUBROUTINE exxinit
  !-----------------------------------------------------------------------
  !
  !-----------------------------------------------------------------------
  SUBROUTINE compute_becxx ( )
    !-----------------------------------------------------------------------
    !
    ! prepare the necessary quantities, then call calbec to compute <beta_I|phi_j,k+q>
    ! and store it becxx(ikq). This must be called AFTER exxbuff and xkq_collected are done
    ! (i.e. at the end of exxinit)
    !
    USE kinds,                ONLY : DP
    USE wvfct,                ONLY : g2kin, npwx, ecutwfc, nbnd
    USE gvect,                ONLY : g, ngm
    USE gvecs,                ONLY : nls, nlsm
    USE cell_base,            ONLY : tpiba2
    USE uspp,                 ONLY : nkb, okvan
    USE becmod,               ONLY : calbec
    USE fft_base,             ONLY : dffts
    USE fft_interfaces,       ONLY : fwfft
    USE control_flags,        ONLY : gamma_only
    USE us_exx,               ONLY : becxx
    USE mp_bands,             ONLY : ibnd_start, ibnd_end

    IMPLICIT NONE
    !
    INTEGER  :: npwq, ibnd, i, ikq, j, h_ibnd, ibnd_loop_start
    REAL(DP) :: gcutwfc
    INTEGER,ALLOCATABLE     :: igkq(:)   ! order of wavefunctions at k+q[+G]
    COMPLEX(DP),ALLOCATABLE :: vkbq(:,:) ! |beta_I> 
    COMPLEX(DP),ALLOCATABLE :: evcq(:,:) ! |psi_j,k> in g-space
    COMPLEX(DP),ALLOCATABLE :: phi(:)    ! aux space for fwfft
    COMPLEX(DP) :: fp, fm
    !
    ! NOTE: I do not want to use vkb from uspp, as you never know if it is going to be used again or not,
    !       this way we are wasting some memory, but the fault is with uspp that should not use global
    !       variables for temporary data (lp-2012-10-03)
    !
    IF(.not. okvan) RETURN
    !
    CALL start_clock('becxx')
    !
    gcutwfc = ecutwfc / tpiba2
    ALLOCATE(igkq(npwx))
    ALLOCATE(vkbq(npwx,nkb))
    ALLOCATE(phi(dffts%nnr))
    ALLOCATE(evcq(npwx,nbnd))
    !
    DO ikq = 1,nkqs
      ! each pool only does its own k-points, then it calls mp_sum (to be tested)
      ! bands count is reset at each k-point
      !
      ! prepare the g-vectors mapping
      CALL gk_sort(xkq_collect(:, ikq), ngm, g, gcutwfc, npwq, igkq, g2kin )
      ! prepare the |beta> function at k+q
      CALL init_us_2(npwq, igkq, xkq_collect(:, ikq), vkbq)
      !
      ! take rotated phi to G space
      IF (gamma_only) THEN
         !
         h_ibnd=ibnd_start/2
         !
         IF(MOD(ibnd_start,2)==0) THEN
            h_ibnd=h_ibnd-1
            ibnd_loop_start=ibnd_start-1
         ELSE
            ibnd_loop_start=ibnd_start
         ENDIF

         DO ibnd = ibnd_loop_start,ibnd_end,2
            h_ibnd = h_ibnd + 1
            phi(:) = exxbuff(:,h_ibnd,ikq)
            CALL fwfft ('Wave', phi, dffts)
            IF (ibnd < ibnd_end) THEN
               ! two ffts at the same time
               DO j = 1, npwq
                  fp = (phi (nls(igkq(j))) + phi (nlsm(igkq(j))))*0.5d0
                  fm = (phi (nls(igkq(j))) - phi (nlsm(igkq(j))))*0.5d0
                  evcq( j, ibnd)   = CMPLX( DBLE(fp), AIMAG(fm),kind=DP)
                  evcq( j, ibnd+1) = CMPLX(AIMAG(fp),- DBLE(fm),kind=DP)
               ENDDO
            ELSE
               DO j = 1, npwq
                  evcq(j, ibnd)   =  phi(nls(igkq(j)))
               ENDDO
            ENDIF
         ENDDO
      ELSE
         DO ibnd = ibnd_start,ibnd_end
            phi(:) = exxbuff(:,ibnd,ikq)
            CALL fwfft ('Wave', phi, dffts)
            FORALL(i=1:npwq) evcq(i,ibnd) = phi(nls(igkq(i)))
         ENDDO
      ENDIF
      !
      ! compute <beta_I|psi_j> at this k+q point, for all band and all projectors
      CALL calbec(npwq, vkbq, evcq, becxx(ikq), nbnd)
      !
    ENDDO
    !
    ! only work for k (only to be called once...):
    ! CALL mp_sum(becxx%k, inter_pool_comm)
    !
    DEALLOCATE(igkq, vkbq, phi, evcq)
    !
    CALL stop_clock('becxx')
    !-----------------------------------------------------------------------
  END SUBROUTINE compute_becxx
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
    USE control_flags,  ONLY : gamma_only
    USE uspp,           ONLY : okvan
    USE paw_variables,  ONLY : okpaw
    !
    IMPLICIT NONE
    !
    INTEGER                  :: lda, n, m
    COMPLEX(DP)              :: psi(lda*npol,m) 
    COMPLEX(DP)              :: hpsi(lda*npol,m)
    TYPE(bec_type), OPTIONAL :: becpsi
    !
    IF ( (okvan.OR.okpaw) .AND. .NOT. PRESENT(becpsi)) &
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
    USE gvecs,          ONLY : nls, ngms
    USE wvfct,          ONLY : npwx, npw, igk, current_k, ecutwfc
    USE control_flags,  ONLY : gamma_only
    USE klist,          ONLY : xk, nks, nkstot
    USE fft_base,       ONLY : dffts
    USE fft_interfaces, ONLY : fwfft, invfft
    USE becmod,         ONLY : bec_type
    USE mp_bands,       ONLY : ibnd_start, ibnd_end, inter_bgrp_comm, &
                               intra_bgrp_comm, my_bgrp_id, nbgrp
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
    COMPLEX(DP),ALLOCATABLE :: result(:)
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
    INTEGER  :: find_current_k
    LOGICAL :: l_fft_doubleband
    LOGICAL :: l_fft_singleband
    !
    ALLOCATE( fac(exx_fft_r2g%ngmt) )
    nrxxs= exx_fft_g2r%dfftt%nnr
    !
    ALLOCATE( result(nrxxs), temppsic_dble(nrxxs), temppsic_aimag(nrxxs) )
    !
    ALLOCATE(rhoc(nrxxs), vc(nrxxs))
    IF(okvan) ALLOCATE(deexx(nkb))
    !
    current_ik=find_current_k(current_k,nkstot,nks)
    xkp = xk_collect(:,current_ik)
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
       CALL g2_convolution(exx_fft_r2g%ngmt, exx_fft_r2g%gt, xk(:,current_k), xkq, fac) 
       IF ( okvan .AND..NOT.tqr ) CALL qvan_init (xkq, xkp)
       !
       LOOP_ON_PSI_BANDS : &
       DO im = 1,m !for each band of psi (the k cycle is outside band)
          IF(okvan) deexx(:) = 0.0_DP
          !
          result = 0.0_DP
          !
          l_fft_doubleband = .FALSE.
          l_fft_singleband = .FALSE.
          !
          IF ( MOD(im,2)==1 .AND. (im+1)<=m ) l_fft_doubleband = .TRUE.
          IF ( MOD(im,2)==1 .AND. im==m )     l_fft_singleband = .TRUE.
          !
          IF( l_fft_doubleband ) THEN 
!$omp parallel do  default(shared), private(ig)
             DO ig = 1, exx_fft_g2r%npwt
                result( exx_fft_g2r%nlt(ig) )  =       psi(ig, im) + (0._DP,1._DP) * psi(ig, im+1)
                result( exx_fft_g2r%nltm(ig) ) = CONJG(psi(ig, im) - (0._DP,1._DP) * psi(ig, im+1))
             ENDDO
!$omp end parallel do
          ENDIF
          !
          IF( l_fft_singleband ) THEN 
!$omp parallel do  default(shared), private(ig)
             DO ig = 1, exx_fft_g2r%npwt
                result( exx_fft_g2r%nlt(ig) )  =       psi(ig,im) 
                result( exx_fft_g2r%nltm(ig) ) = CONJG(psi(ig,im))
             ENDDO
!$omp end parallel do
          ENDIF
          !
          IF( l_fft_doubleband.OR.l_fft_singleband) THEN
             CALL invfft ('CustomWave', result, exx_fft_g2r%dfftt)
!$omp parallel do default(shared), private(ir)
             DO ir = 1, nrxxs
                temppsic_dble(ir)  = DBLE ( result(ir) )
                temppsic_aimag(ir) = AIMAG( result(ir) )
             ENDDO
!$omp end parallel do
          ENDIF
          !
          result = 0.0_DP
          !
          h_ibnd = ibnd_start/2
          IF(MOD(ibnd_start,2)==0) THEN
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
             IF ( ABS(x1) < eps_occ .AND. ABS(x2) < eps_occ ) CYCLE
             !
             ! calculate rho in real space. Gamma tricks are used. 
             ! temppsic is real; tempphic contains one band in the real part, 
             ! another one in the imaginary part; the same applies to rhoc
             !
             IF( MOD(im,2) == 0 ) THEN 
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
             IF(okvan .AND. tqr) THEN
                IF(ibnd>=ibnd_start) &
                CALL addusxx_r(rhoc, _CX(becxx(ikq)%r(:,ibnd)), _CX(becpsi%r(:,im)))
                IF(ibnd<ibnd_end) &
                CALL addusxx_r(rhoc,_CY(becxx(ikq)%r(:,ibnd+1)),_CX(becpsi%r(:,im)))
             ENDIF
             !
             CALL fwfft ('Custom', rhoc, exx_fft_r2g%dfftt)
             !   >>>> add augmentation in G SPACE here
             IF(okvan .AND. .NOT. TQR) THEN
                ! contribution from one band added to real (in real space) part of rhoc
                IF(ibnd>=ibnd_start) &
                   CALL addusxx_g(rhoc, xkq,  xkp, 'r', &
                   becphi_r=becxx(ikq)%r(:,ibnd), becpsi_r=becpsi%r(:,im) )
                ! contribution from following band added to imaginary (in real space) part of rhoc
                IF(ibnd<ibnd_end) &
                   CALL addusxx_g(rhoc, xkq,  xkp, 'i', &
                   becphi_r=becxx(ikq)%r(:,ibnd+1), becpsi_r=becpsi%r(:,im) )
             ENDIF
             !   >>>> charge density done
             !
             vc = 0._DP
             !
!$omp parallel do default(shared), private(ig)
             DO ig = 1, exx_fft_r2g%ngmt
                ! 
                vc(exx_fft_r2g%nlt(ig))  = fac(ig) * rhoc(exx_fft_r2g%nlt(ig)) 
                vc(exx_fft_r2g%nltm(ig)) = fac(ig) * rhoc(exx_fft_r2g%nltm(ig)) 
                !                 
             ENDDO
!$omp end parallel do
             !
             !   >>>>  compute <psi|H_fock G SPACE here
             IF(okvan .and. .not. TQR) THEN
                IF(ibnd>=ibnd_start) &
                CALL newdxx_g(vc, xkq, xkp, 'r', deexx, becphi_r=x1*becxx(ikq)%r(:,ibnd))
                IF(ibnd<ibnd_end) &
                CALL newdxx_g(vc, xkq, xkp, 'i', deexx,becphi_r=x2*becxx(ikq)%r(:,ibnd+1))
             ENDIF
             !
             !brings back v in real space
             CALL invfft ('Custom', vc, exx_fft_r2g%dfftt) 
             !
             !   >>>>  compute <psi|H_fock REAL SPACE here
             IF(okvan .and. tqr) THEN
                IF(ibnd>=ibnd_start) &
                CALL newdxx_r(vc, CMPLX(x1*becxx(ikq)%r(:,ibnd), 0.0_DP, KIND=DP), deexx)
                IF(ibnd<ibnd_end) &
                CALL newdxx_r(vc, CMPLX(0.0_DP,-x2*becxx(ikq)%r(:,ibnd+1), KIND=DP), deexx)
             ENDIF
             !
             IF(okpaw) THEN
                IF(ibnd>=ibnd_start) &
                CALL PAW_newdxx(x1/nqs, _CX(becxx(ikq)%r(:,ibnd)), _CX(becpsi%r(:,im)), deexx)
                IF(ibnd<ibnd_end) &
                CALL PAW_newdxx(x2/nqs, _CX(becxx(ikq)%r(:,ibnd+1)), _CX(becpsi%r(:,im)), deexx)
             ENDIF
             !
             ! accumulates over bands and k points
             !
!$omp parallel do default(shared), private(ir)
             DO ir = 1, nrxxs
                result(ir) = result(ir)+x1* DBLE(vc(ir))* DBLE(exxbuff(ir,h_ibnd,ikq))&
                                       +x2*AIMAG(vc(ir))*AIMAG(exxbuff(ir,h_ibnd,ikq))
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
          CALL mp_sum( result(1:nrxxs), inter_bgrp_comm)
          !
          ! brings back result in G-space
          !
          CALL fwfft( 'CustomWave' , result, exx_fft_g2r%dfftt )
          !
!$omp parallel do default(shared), private(ig)
          DO ig = 1, n
             hpsi(ig,im)=hpsi(ig,im) - exxalfa*result(exx_fft_g2r%nlt(ig))
          ENDDO
!$omp end parallel do
          ! add non-local \sum_I |beta_I> \alpha_Ii (the sum on i is outside)
          IF(okvan) CALL add_nlxx_pot (lda, hpsi(:,im), xkp, npw, igk, &
                                       deexx, eps_occ, exxalfa)
       ENDDO &
       LOOP_ON_PSI_BANDS
       IF ( okvan .AND..NOT.tqr ) CALL qvan_clean ()
       !
    ENDDO &
    INTERNAL_LOOP_ON_Q
    !  
    DEALLOCATE( result, temppsic_dble, temppsic_aimag) 
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
    USE gvecs,          ONLY : nls, ngms
    USE wvfct,          ONLY : npwx, npw, igk, current_k, ecutwfc
    USE control_flags,  ONLY : gamma_only
    USE klist,          ONLY : xk, nks, nkstot
    USE fft_base,       ONLY : dffts
    USE fft_interfaces, ONLY : fwfft, invfft
    USE becmod,         ONLY : bec_type
    USE mp_bands,       ONLY : ibnd_start, ibnd_end, inter_bgrp_comm, &
                               intra_bgrp_comm, my_bgrp_id, nbgrp
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
    COMPLEX(DP),ALLOCATABLE :: temppsic(:), result(:)
    COMPLEX(DP),ALLOCATABLE :: temppsic_nc(:,:),result_nc(:,:)
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
    INTEGER  :: find_current_k
    !
    ALLOCATE( fac(ngm) )
    nrxxs= dffts%nnr
    !
    IF (noncolin) THEN
       ALLOCATE( temppsic_nc(nrxxs,npol), result_nc(nrxxs,npol) )
    ELSE
       ALLOCATE( temppsic(nrxxs), result(nrxxs) )
    ENDIF
    !
    ALLOCATE(rhoc(nrxxs), vc(nrxxs))
    IF(okvan) ALLOCATE(deexx(nkb))
    !
    current_ik=find_current_k(current_k,nkstot,nks)
    xkp = xk_collect(:,current_ik)
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
             temppsic_nc(nls(igk(ig)),1) = psi(ig,im)
          ENDDO
!$omp end parallel do
!$omp parallel do  default(shared), private(ig)
          DO ig = 1, n
             temppsic_nc(nls(igk(ig)),2) = psi(npwx+ig,im)
          ENDDO
!$omp end parallel do
          !
          CALL invfft ('Wave', temppsic_nc(:,1), dffts)
          CALL invfft ('Wave', temppsic_nc(:,2), dffts)
          !
       ELSE
          !
!$omp parallel do  default(shared), private(ig)
          DO ig = 1, n
             temppsic( nls(igk(ig)) ) = psi(ig,im)
          ENDDO
!$omp end parallel do
          CALL invfft ('Wave', temppsic, dffts)
          !
       ENDIF
       !
       IF (noncolin) THEN
          result_nc = 0.0_DP
       ELSE
          result    = 0.0_DP
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
          CALL g2_convolution(ngms, g, xk(:,current_k), xkq, fac)
          IF ( okvan .AND..NOT.tqr ) CALL qvan_init (xkq, xkp)
          !
          IBND_LOOP_K : &
          DO ibnd=ibnd_start,ibnd_end !for each band of psi
             !
             IF ( ABS(x_occupation(ibnd,ik)) < eps_occ) CYCLE IBND_LOOP_K
             !
             !loads the phi from file
             !
             !   >>>> calculate rho in real space
             IF (noncolin) THEN
!$omp parallel do default(shared), private(ir)
                DO ir = 1, nrxxs
                   rhoc(ir) = ( CONJG(exxbuff(ir,ibnd,ikq))*temppsic_nc(ir,1) + &
                                 CONJG(exxbuff(nrxxs+ir,ibnd,ikq))*temppsic_nc(ir,2) )/omega
                ENDDO
!$omp end parallel do
             ELSE
!$omp parallel do default(shared), private(ir)
                DO ir = 1, nrxxs
                   rhoc(ir)=CONJG(exxbuff(ir,ibnd,ikq))*temppsic(ir) / omega
                ENDDO
!$omp end parallel do
             ENDIF
             !   >>>> add augmentation in REAL space HERE
             IF(okvan .AND. tqr) THEN ! augment the "charge" in real space
                CALL addusxx_r(rhoc, becxx(ikq)%k(:,ibnd), becpsi%k(:,im))
             ENDIF
             !
             !   >>>> brings it to G-space
             CALL fwfft('Smooth', rhoc, dffts)
             !
             !   >>>> add augmentation in G space HERE
             IF(okvan .AND. .NOT. tqr) THEN
                CALL addusxx_g(rhoc, xkq, xkp, 'c', &
                   becphi_c=becxx(ikq)%k(:,ibnd),becpsi_c=becpsi%k(:,im))
             ENDIF
             !   >>>> charge done
             !
             vc = 0._DP
             !
!$omp parallel do default(shared), private(ig)
             DO ig = 1, ngms
                vc(nls(ig)) = fac(ig) * rhoc(nls(ig)) * x_occupation(ibnd,ik) / nqs
             ENDDO
!$omp end parallel do
             !
             ! Add ultrasoft contribution (RECIPROCAL SPACE)
             ! compute alpha_I,j,k+q = \sum_J \int <beta_J|phi_j,k+q> V_i,j,k,q Q_I,J(r) d3r
             IF(okvan .AND. .NOT. tqr) THEN
                CALL newdxx_g(vc, xkq, xkp, 'c', deexx, becphi_c=becxx(ikq)%k(:,ibnd))
             ENDIF
             !
             !brings back v in real space
             CALL invfft ('Smooth', vc, dffts)
             !
             ! Add ultrasoft contribution (REAL SPACE)
             IF(okvan .AND. TQR) CALL newdxx_r(vc, becxx(ikq)%k(:,ibnd),deexx)
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
                   result(ir) = result(ir) + vc(ir)*exxbuff(ir,ibnd,ikq)
                ENDDO
!$omp end parallel do
             ENDIF
             !
          ENDDO &
          IBND_LOOP_K
          IF ( okvan .AND..NOT.tqr ) CALL qvan_clean ()
          !
       ENDDO &
       INTERNAL_LOOP_ON_Q
       !
       IF(okvan) THEN
         CALL mp_sum(deexx,intra_bgrp_comm)
         CALL mp_sum(deexx,inter_bgrp_comm)
       ENDIF
       !
       IF (noncolin) THEN
          CALL mp_sum( result_nc(1:nrxxs,1:npol), inter_bgrp_comm)
       ELSE
          CALL mp_sum( result(1:nrxxs), inter_bgrp_comm)
       ENDIF
       !
       !brings back result in G-space
       !
       IF (noncolin) THEN
          !brings back result in G-space
          CALL fwfft ('Wave', result_nc(:,1), dffts)
          CALL fwfft ('Wave', result_nc(:,2), dffts)
          !
          !adds it to hpsi
!$omp parallel do default(shared), private(ig)
          DO ig = 1, n
             hpsi(ig,im)     = hpsi(ig,im)     - exxalfa*result_nc(nls(igk(ig)),1)
          ENDDO
!$omp end parallel do
!$omp parallel do default(shared), private(ig)
          DO ig = 1, n
             hpsi(lda+ig,im) = hpsi(lda+ig,im) - exxalfa*result_nc(nls(igk(ig)),2)
          ENDDO
!$omp end parallel do
          !
       ELSE
          !
          CALL fwfft ('Wave', result, dffts)
          !
          !adds it to hpsi
!$omp parallel do default(shared), private(ig)
          DO ig = 1, n
             hpsi(ig,im)=hpsi(ig,im) - exxalfa*result(nls(igk(ig)))
          ENDDO
!$omp end parallel do
       ENDIF
       !
       ! add non-local \sum_I |beta_I> \alpha_Ii (the sum on i is outside)
       IF(okvan) CALL add_nlxx_pot (lda, hpsi(:,im), xkp, npw, igk, &
                                       deexx, eps_occ, exxalfa)
       !
    ENDDO &
    LOOP_ON_PSI_BANDS
    !
    IF (noncolin) THEN
       DEALLOCATE(temppsic_nc, result_nc) 
    ELSE
       DEALLOCATE(temppsic, result) 
    END IF
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
    INTEGER,  INTENT(IN)    :: ngm   ! Number of G vectors
    REAL(DP), INTENT(IN)    :: g(3,ngm) ! Cartesian components of G vectors
    REAL(DP), INTENT(IN)    :: xk(3) ! current k vector
    REAL(DP), INTENT(IN)    :: xkq(3) ! current q vector
    !
    REAL(DP), INTENT(INOUT) :: fac(ngm) ! Calculated convolution
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
    nqhalf_dble(1:3) = (/ DBLE(nq1)*0.5_DP, DBLE(nq2)*0.5_DP, DBLE(nq3)*0.5_DP /) 
    !
    ! Set the grid_factor_track and qq_track
    !
    IF( x_gamma_extrapolation ) THEN 
!$omp parallel do default(shared), private(ig,q,x,odg)
       DO ig = 1, ngm 
          q(:)= xk(:) - xkq(:) + g(:,ig) 
          qq_track(ig) = SUM(q(:)**2) * tpiba2
          x = (q(1)*at(1,1)+q(2)*at(2,1)+q(3)*at(3,1))*nqhalf_dble(1)
          odg(1) = ABS(x-NINT(x))<eps
          x = (q(1)*at(1,2)+q(2)*at(2,2)+q(3)*at(3,2))*nqhalf_dble(2)
          odg(2) = ABS(x-NINT(x))<eps
          x = (q(1)*at(1,3)+q(2)*at(2,3)+q(3)*at(3,3))*nqhalf_dble(3)
          odg(3) = ABS(x-NINT(x))<eps
          IF( ALL ( odg(:) ) ) THEN
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
          qq_track(ig) = SUM(q(:)**2) * tpiba2
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
         fac(ig)=e2*((pi/gau_scrlen)**(1.5_DP))*EXP(-qq/4._DP/gau_scrlen) * grid_factor_track(ig)
         !
      ELSE IF (qq > eps_qdiv) THEN
         !
         IF ( erfc_scrlen > 0  ) THEN
            fac(ig)=e2*fpi/qq*(1._DP-EXP(-qq/4._DP/erfc_scrlen**2)) * grid_factor_track(ig)
         ELSEIF( erf_scrlen > 0 ) THEN
            fac(ig)=e2*fpi/qq*(EXP(-qq/4._DP/erf_scrlen**2)) * grid_factor_track(ig)
         ELSE
            fac(ig)=e2*fpi/( qq + yukawa ) * grid_factor_track(ig) ! as HARTREE
         ENDIF
         !
      ELSE
         !
         fac(ig)= - exxdiv ! or rather something ELSE (see F.Gygi)
         !
         IF ( yukawa > 0._DP.AND. .NOT. x_gamma_extrapolation ) fac(ig) = fac(ig) + e2*fpi/( qq + yukawa )
         IF( erfc_scrlen > 0._DP.AND. .NOT. x_gamma_extrapolation ) fac(ig) = fac(ig) + e2*pi/(erfc_scrlen**2)
         !
      ENDIF
      !
    ENDDO
!$omp end parallel do
  END SUBROUTINE g2_convolution
  !-----------------------------------------------------------------------
  !
  !-----------------------------------------------------------------------
  FUNCTION exxenergy2()
    !-----------------------------------------------------------------------
    !
    USE control_flags,           ONLY : gamma_only
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
    USE io_files,                ONLY : iunigk,iunwfc, nwordwfc
    USE buffers,                 ONLY : get_buffer
    USE cell_base,               ONLY : alat, omega, bg, at, tpiba
    USE symm_base,               ONLY : nsym, s
    USE gvect,                   ONLY : ngm, gstart, g, nl
    USE gvecs,                   ONLY : ngms, nls, nlsm, doublegrid
    USE wvfct,                   ONLY : nbnd, npwx, npw, igk, wg, ecutwfc
    USE control_flags,           ONLY : gamma_only
    USE wavefunctions_module,    ONLY : evc
    USE klist,                   ONLY : xk, ngk, nks, nkstot
    USE lsda_mod,                ONLY : lsda, current_spin, isk
    USE mp_pools,                ONLY : inter_pool_comm
    USE mp_bands,                ONLY : inter_bgrp_comm, intra_bgrp_comm, &
                                        nbgrp, ibnd_start, ibnd_end
    USE mp,                      ONLY : mp_sum
    USE fft_base,                ONLY : dffts
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
    INTEGER,        EXTERNAL :: find_current_k
    !
    TYPE(bec_type) :: becpsi
    COMPLEX(DP), ALLOCATABLE :: psi_t(:), prod_tot(:)
    INTEGER,     ALLOCATABLE :: igkt(:)
    REAL(DP),ALLOCATABLE :: temppsic_dble (:)
    REAL(DP),ALLOCATABLE :: temppsic_aimag(:)
    LOGICAL :: l_fft_doubleband
    LOGICAL :: l_fft_singleband
    INTEGER :: jmax
    !
    nrxxs= exx_fft_g2r%dfftt%nnr
    ALLOCATE( fac(exx_fft_r2g%ngmt) )
    !
    ALLOCATE(temppsic(nrxxs), temppsic_dble(nrxxs),temppsic_aimag(nrxxs)) 
    ALLOCATE( rhoc(nrxxs) )
    !
    energy=0.0_DP
    !
    CALL allocate_bec_type( nkb, nbnd, becpsi)
    !
    IF ( nks > 1 ) REWIND( iunigk )
    !
    IKK_LOOP : &
    DO ikk=1,nks
       current_ik=find_current_k(ikk,nkstot,nks)
       xkp = xk_collect(:,current_ik)
       !
       IF ( lsda ) current_spin = isk(ikk)
       npw = ngk (ikk)
       IF ( nks > 1 ) THEN
          READ( iunigk ) igk
          CALL get_buffer (evc, nwordwfc, iunwfc, ikk)
       END IF
       !
       ! prepare the |beta> function at k+q
       CALL init_us_2(npw, igk, xk(:,ikk), vkb)
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
          CALL g2_convolution(exx_fft_r2g%ngmt, exx_fft_r2g%gt, xk(:,current_ik), xkq, fac) 
          fac(exx_fft_r2g%gstart_t:) = 2 * fac(exx_fft_r2g%gstart_t:)
          IF ( okvan .AND..NOT.tqr ) CALL qvan_init (xkq, xkp)
          !
          jmax = nbnd 
          DO jbnd = nbnd,1, -1
             IF ( ABS(wg(jbnd,ikk)) < eps_occ) CYCLE
             jmax = jbnd 
             EXIT
          ENDDO
          !
          JBND_LOOP : &
          DO jbnd = 1, jmax     !for each band of psi (the k cycle is outside band)
             !
             temppsic = 0._DP
             !
             l_fft_doubleband = .FALSE.
             l_fft_singleband = .FALSE.
             !
             IF ( MOD(jbnd,2)==1 .AND. (jbnd+1)<=jmax ) l_fft_doubleband = .TRUE.
             IF ( MOD(jbnd,2)==1 .AND. jbnd==jmax )     l_fft_singleband = .TRUE.
             !
             IF( l_fft_doubleband ) THEN 
!$omp parallel do  default(shared), private(ig)
                DO ig = 1, exx_fft_g2r%npwt
                   temppsic( exx_fft_g2r%nlt(ig) )  =       evc(ig,jbnd) + (0._DP,1._DP) * evc(ig,jbnd+1)
                   temppsic( exx_fft_g2r%nltm(ig) ) = CONJG(evc(ig,jbnd) - (0._DP,1._DP) * evc(ig,jbnd+1))
                ENDDO
!$omp end parallel do
             ENDIF
             !
             IF( l_fft_singleband ) THEN 
!$omp parallel do  default(shared), private(ig)
                DO ig = 1, exx_fft_g2r%npwt
                   temppsic( exx_fft_g2r%nlt(ig) )  =       evc(ig,jbnd) 
                   temppsic( exx_fft_g2r%nltm(ig) ) = CONJG(evc(ig,jbnd))
                ENDDO
!$omp end parallel do
             ENDIF
             !
             IF( l_fft_doubleband.OR.l_fft_singleband) THEN
                CALL invfft ('CustomWave', temppsic, exx_fft_g2r%dfftt)
!$omp parallel do default(shared), private(ir)
                DO ir = 1, nrxxs
                   temppsic_dble(ir)  = DBLE ( temppsic(ir) )
                   temppsic_aimag(ir) = AIMAG( temppsic(ir) )
                ENDDO
!$omp end parallel do
             ENDIF
             !
             h_ibnd = ibnd_start/2
             IF(MOD(ibnd_start,2)==0) THEN
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
                IF( MOD(jbnd,2) == 0 ) THEN
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
                   CALL addusxx_r(rhoc, _CX(becxx(ikq)%r(:,ibnd)), _CX(becpsi%r(:,jbnd)))
                   IF(ibnd<ibnd_end) &
                   CALL addusxx_r(rhoc,_CY(becxx(ikq)%r(:,ibnd+1)),_CX(becpsi%r(:,jbnd)))
                ENDIF
                !
                ! bring rhoc to G-space
                CALL fwfft ('Custom', rhoc, exx_fft_r2g%dfftt)
                !
                IF(okvan .and..not.tqr) THEN
                   IF(ibnd>=ibnd_start ) &
                      CALL addusxx_g( rhoc, xkq, xkp, 'r', &
                      becphi_r=becxx(ikq)%r(:,ibnd), becpsi_r=becpsi%r(:,jbnd) )
                   IF(ibnd<ibnd_end) &
                      CALL addusxx_g( rhoc, xkq, xkp, 'i', &
                      becphi_r=becxx(ikq)%r(:,ibnd+1), becpsi_r=becpsi%r(:,jbnd) )
                ENDIF
                !
                vc = 0.0_DP
!$omp parallel do  default(shared), private(ig),  reduction(+:vc)
                DO ig = 1,exx_fft_r2g%ngmt
                   !
                   ! The real part of rhoc contains the contribution from band ibnd
                   ! The imaginary part    contains the contribution from band ibnd+1
                   !
                   vc = vc + fac(ig) * ( x1 * &
                        ABS( rhoc(exx_fft_r2g%nlt(ig)) + CONJG(rhoc(exx_fft_r2g%nltm(ig))) )**2 &
                                        +x2 * &
                        ABS( rhoc(exx_fft_r2g%nlt(ig)) - CONJG(rhoc(exx_fft_r2g%nltm(ig))) )**2 )
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
          IF ( okvan .AND..NOT.tqr ) CALL qvan_clean ( )
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
    USE io_files,                ONLY : iunigk,iunwfc, nwordwfc
    USE buffers,                 ONLY : get_buffer
    USE cell_base,               ONLY : alat, omega, bg, at, tpiba
    USE symm_base,               ONLY : nsym, s
    USE gvect,                   ONLY : ngm, gstart, g, nl
    USE gvecs,                   ONLY : ngms, nls, nlsm, doublegrid
    USE wvfct,                   ONLY : nbnd, npwx, npw, igk, wg, ecutwfc
    USE control_flags,           ONLY : gamma_only
    USE wavefunctions_module,    ONLY : evc
    USE klist,                   ONLY : xk, ngk, nks, nkstot
    USE lsda_mod,                ONLY : lsda, current_spin, isk
    USE mp_pools,                ONLY : inter_pool_comm
    USE mp_bands,                ONLY : inter_bgrp_comm, intra_bgrp_comm, &
                                        nbgrp, ibnd_start, ibnd_end
    USE mp,                      ONLY : mp_sum
    USE fft_base,                ONLY : dffts
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
    INTEGER  :: jbnd, ibnd, ik, ikk, ig, ikq, iq, ir
    INTEGER  :: h_ibnd, nrxxs, current_ik, ibnd_loop_start
    REAL(DP) :: x1, x2
    REAL(DP) :: xkq(3), xkp(3), vc
    ! temp array for vcut_spheric
    INTEGER,        EXTERNAL :: find_current_k
    !
    TYPE(bec_type) :: becpsi
    COMPLEX(DP), ALLOCATABLE :: psi_t(:), prod_tot(:)
    INTEGER,     ALLOCATABLE :: igkt(:)
    !
    nrxxs = dffts%nnr
    ALLOCATE( fac(ngms) )
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
    IF ( nks > 1 ) REWIND( iunigk )
    !
    IKK_LOOP : &
    DO ikk=1,nks
       current_ik=find_current_k(ikk,nkstot,nks)
       xkp = xk_collect(:,current_ik)
       !
       IF ( lsda ) current_spin = isk(ikk)
       npw = ngk (ikk)
       IF ( nks > 1 ) THEN
          READ( iunigk ) igk
          CALL get_buffer (evc, nwordwfc, iunwfc, ikk)
       END IF
       !
       ! prepare the |beta> function at k+q
       CALL init_us_2(npw, igk, xk(:,ikk), vkb)
       ! compute <beta_I|psi_j> at this k+q point, for all band and all projectors
       CALL calbec(npw, vkb, evc, becpsi, nbnd)
       !
       JBND_LOOP : &
       DO jbnd = 1, nbnd     !for each band of psi (the k cycle is outside band)
          !
          IF ( ABS(wg(jbnd,ikk)) < eps_occ) CYCLE
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
                temppsic_nc(nls(igk(ig)),1) = evc(ig,jbnd)
             ENDDO
!$omp end parallel do
!$omp parallel do default(shared), private(ig)
             DO ig = 1, npw
                temppsic_nc(nls(igk(ig)),2) = evc(npwx+ig,jbnd)
             ENDDO
!$omp end parallel do
             !
             CALL invfft ('Wave', temppsic_nc(:,1), dffts)
             CALL invfft ('Wave', temppsic_nc(:,2), dffts)
             !
          ELSE
!$omp parallel do default(shared), private(ig)
             DO ig = 1, npw
                temppsic(nls(igk(ig))) = evc(ig,jbnd)
             ENDDO
!$omp end parallel do
             !
             CALL invfft ('Wave', temppsic, dffts)
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
             CALL g2_convolution(ngms, g, xk(:,current_ik), xkq, fac)
             IF ( okvan .AND..NOT.tqr ) CALL qvan_init (xkq, xkp)
             !
             IBND_LOOP_K : &
             DO ibnd = ibnd_start, ibnd_end
                !
                IF ( ABS(x_occupation(ibnd,ik)) < eps_occ) CYCLE
                !
                ! load the phi at this k+q and band
                IF (noncolin) THEN
                   !
!$omp parallel do  default(shared), private(ir) 
                   DO ir = 1, nrxxs
                      rhoc(ir)=(CONJG(exxbuff(ir      ,ibnd,ikq))*temppsic_nc(ir,1) + &
                                CONJG(exxbuff(ir+nrxxs,ibnd,ikq))*temppsic_nc(ir,2) )/omega
                   ENDDO
!$omp end parallel do
                ELSE
                   !calculate rho in real space
!$omp parallel do  default(shared), private(ir)
                   DO ir = 1, nrxxs
                      rhoc(ir)=CONJG(exxbuff(ir,ibnd,ikq))*temppsic(ir) / omega
                   ENDDO
!$omp end parallel do
                ENDIF
                ! augment the "charge" in real space
                IF(okvan .AND. tqr) CALL addusxx_r(rhoc, becxx(ikq)%k(:,ibnd), becpsi%k(:,jbnd))
                !
                ! bring rhoc to G-space
                CALL fwfft ('Smooth', rhoc, dffts)
                ! augment the "charge" in G space
                IF(okvan .AND. .NOT. tqr) & 
                   CALL addusxx_g(rhoc, xkq, xkp, 'c', &
                   becphi_c=becxx(ikq)%k(:,ibnd),becpsi_c=becpsi%k(:,jbnd))
                !
                vc = 0.0_DP
!$omp parallel do  default(shared), private(ig), reduction(+:vc)
                DO ig=1,ngms
                   vc = vc + fac(ig) * DBLE(rhoc(nls(ig))*CONJG(rhoc(nls(ig))))
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
             IF ( okvan .AND..NOT.tqr ) CALL qvan_clean ( )
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
     USE wvfct,          ONLY : ecutwfc
     USE io_global,      ONLY : stdout
     USE control_flags,  ONLY : gamma_only
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

     alpha  = 10._dp * tpiba2 / ecutwfc

     IF ( .NOT. use_regularization ) THEN
        exx_divergence = 0._dp
        RETURN
     END IF

     dq1= 1._dp/DBLE(nq1)
     dq2= 1._dp/DBLE(nq2) 
     dq3= 1._dp/DBLE(nq3) 

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
     if( erf_scrlen > 0) aa = 1._dp/sqrt((alpha+1._dp/4.d0/erf_scrlen**2)*0.25d0*fpi)
     div = div - e2*omega * aa

     exx_divergence = div * nqs
     CALL stop_clock ('exx_div')

     return
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
    USE io_files,             ONLY : iunigk,iunwfc, nwordwfc
    USE buffers,              ONLY : get_buffer
    USE cell_base,            ONLY : alat, omega, bg, at, tpiba
    USE symm_base,            ONLY : nsym, s
    USE gvect,                ONLY : ngm
    USE gvecs,                ONLY : ngms, nls, nlsm, doublegrid
    USE wvfct,                ONLY : nbnd, npwx, npw, igk, wg, current_k
    USE control_flags,        ONLY : gamma_only
    USE wavefunctions_module, ONLY : evc
    USE klist,                ONLY : xk, ngk, nks
    USE lsda_mod,             ONLY : lsda, current_spin, isk
    USE gvect,                ONLY : g, nl
    USE mp_pools,             ONLY : npool, inter_pool_comm
    USE mp_bands,             ONLY : inter_bgrp_comm, intra_bgrp_comm, &
                                     ibnd_start, ibnd_end
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
    REAL(DP),    allocatable :: fac(:), fac_tens(:,:,:), fac_stress(:)
    INTEGER  :: jbnd, ibnd, ik, ikk, ig, ir, ikq, iq, isym
    INTEGER  :: h_ibnd, nqi, iqi, beta, nrxxs
    INTEGER  :: ibnd_loop_start
    REAL(DP) :: x1, x2
    REAL(DP) :: qq, xk_cryst(3), sxk(3), xkq(3), vc(3,3), x, q(3)
    ! temp array for vcut_spheric
    REAL(DP) :: delta(3,3)
  
    CALL start_clock ('exx_stress')

    IF (npool>1) CALL errore('exx_stress','stress not available with pools',1)
    IF (noncolin) CALL errore('exx_stress','noncolinear stress not implemented',1)
    IF (okvan) CALL infomsg('exx_stress','USPP stress not tested')

    nrxxs = dffts%nnr
    delta = reshape( (/1._dp,0._dp,0._dp, 0._dp,1._dp,0._dp, 0._dp,0._dp,1._dp/), (/3,3/))
    exx_stress_ = 0._dp
    allocate( tempphic(nrxxs), temppsic(nrxxs), rhoc(nrxxs), fac(ngm) )
    allocate( fac_tens(3,3,ngm), fac_stress(ngm) )

    IF ( nks > 1 ) rewind( iunigk )
    !
    nqi=nqs
    !
    ! loop over k-points
    DO ikk = 1, nks
        current_k = ikk
        IF (lsda) current_spin = isk(ikk)
        npw = ngk(ikk)

        IF (nks > 1) THEN
            read(iunigk) igk
            CALL get_buffer(evc, nwordwfc, iunwfc, ikk)
        ENDIF

        ! loop over bands
        DO jbnd = 1, nbnd
            !
            temppsic(:) = ( 0._dp, 0._dp )
!$omp parallel do default(shared), private(ig)
            DO ig = 1, npw
                temppsic(nls(igk(ig))) = evc(ig,jbnd)
            ENDDO
!$omp end parallel do
            !
            IF(gamma_only) THEN
!$omp parallel do default(shared), private(ig)
                DO ig = 1, npw
                    temppsic(nlsm(igk(ig))) = conjg(evc(ig,jbnd))
                ENDDO
!$omp end parallel do
            ENDIF

            CALL invfft ('Wave', temppsic, dffts)       

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
                      on_double_grid = .FALSE.
                  ENDIF

                  IF (use_coulomb_vcut_ws) THEN
                      fac(ig) = vcut_get(vcut, q)
                      fac_stress(ig) = 0._dp   ! not implemented
                      IF (gamma_only .and. qq > 1.d-8) fac(ig) = 2.d0 * fac(ig)

                  ELSE IF ( use_coulomb_vcut_spheric ) THEN
                      fac(ig) = vcut_spheric_get(vcut, q)
                      fac_stress(ig) = 0._dp   ! not implemented
                      IF (gamma_only .and. qq > 1.d-8) fac(ig) = 2.d0 * fac(ig) 

                  ELSE IF (gau_scrlen > 0) then
                      fac(ig)=e2*((pi/gau_scrlen)**(1.5d0))* &
                            exp(-qq/4.d0/gau_scrlen) * grid_factor
                      fac_stress(ig) =  e2*2.d0/4.d0/gau_scrlen * &
                            exp(-qq/4.d0/gau_scrlen) *((pi/gau_scrlen)**(1.5d0))* &
                                                                     grid_factor
                      IF (gamma_only) fac(ig) = 2.d0 * fac(ig)
                      IF (gamma_only) fac_stress(ig) = 2.d0 * fac_stress(ig)
                      IF (on_double_grid) fac(ig) = 0._dp
                      IF (on_double_grid) fac_stress(ig) = 0._dp

                  ELSE IF (qq > 1.d-8) THEN
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
                    IF(MOD(ibnd_start,2)==0) THEN
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
                            rhoc(ir)     = CONJG(tempphic(ir))*temppsic(ir) / omega
                        ENDDO
!$omp end parallel do
                        ! bring it to G-space
                        CALL fwfft ('Smooth', rhoc, dffts)
    
                        vc = 0._dp
!$omp parallel do default(shared), private(ig), reduction(+:vc)
                        DO ig = 1, ngms
                            !
                            vc(:,:) = vc(:,:) + fac(ig) * x1 * &
                                      abs( rhoc(nls(ig))+CONJG(rhoc(nlsm(ig))))**2 * &
                                      (fac_tens(:,:,ig)*fac_stress(ig)/2.d0 - delta(:,:)*fac(ig))
                            vc(:,:) = vc(:,:) + fac(ig) * x2 * &
                                      abs( rhoc(nls(ig))-CONJG(rhoc(nlsm(ig))))**2 * &
                                      (fac_tens(:,:,ig)*fac_stress(ig)/2.d0 - delta(:,:)*fac(ig))
                        enddo
!$omp end parallel do
                        vc = vc / nqs / 4.d0
                        exx_stress_ = exx_stress_ + exxalfa * vc * wg(jbnd,ikk)
                    ENDDO

                ELSE

                    DO ibnd = ibnd_start, ibnd_end    !for each band of psi
                      !
                      IF ( abs(x_occupation(ibnd,ik)) < 1.d-6) cycle
                      !
                      ! calculate rho in real space
!$omp parallel do default(shared), private(ir)
                      DO ir = 1, nrxxs
                          tempphic(ir) = exxbuff(ir,ibnd,ikq)
                          rhoc(ir)     = CONJG(tempphic(ir))*temppsic(ir) / omega
                      ENDDO
!$omp end parallel do

                      ! bring it to G-space
                      CALL fwfft ('Smooth', rhoc, dffts)

                      vc = 0._dp
!$omp parallel do default(shared), private(ig), reduction(+:vc)
                      DO ig = 1, ngms
                          vc(:,:) = vc(:,:) + rhoc(nls(ig))*CONJG(rhoc(nls(ig))) * &
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

    deALLOCATE(tempphic, temppsic, rhoc, fac, fac_tens, fac_stress )
    !
    CALL mp_sum( exx_stress_, intra_bgrp_comm )
    CALL mp_sum( exx_stress_, inter_bgrp_comm )
    CALL mp_sum( exx_stress_, inter_pool_comm )
    exx_stress = exx_stress_

    CALL stop_clock ('exx_stress')
    !-----------------------------------------------------------------------
  END FUNCTION exx_stress
  !-----------------------------------------------------------------------
!-----------------------------------------------------------------------
END MODULE exx
!-----------------------------------------------------------------------
