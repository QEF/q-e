!
! Copyright (C) 2005-2012 Quantum ESPRESSO group
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
  USE io_global,            ONLY : ionode, stdout
  USE fft_custom,           ONLY : fft_cus
  !
  ! FIXME: move down when ready
  USE mp_global,  ONLY : npool

  IMPLICIT NONE
  SAVE

  !
  ! general purpose vars
  !
  REAL(DP):: exxalfa=0._dp                ! 1 if exx, 0 elsewhere
  INTEGER :: exx_nwordwfc, ji

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
                                         ! temporay buffer to store wfc 
  COMPLEX(DP), ALLOCATABLE :: exxbuff_nc(:,:,:,:)
                                         ! temporay buffer to store wfc in the
                                         ! noncollinear case

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
  ! variables to deal with Coulomb divergence
  ! and related issues
  !
  REAL(DP)         :: eps  = 1.d-6
  REAL(DP)         :: eps_qdiv = 1.d-8 ! |q| > eps_qdiv
  REAL(DP)         :: eps_occ  = 1.d-6 ! skip band where occupation is less than this
  REAL(DP)         :: exxdiv = 0._dp
  CHARACTER(32)    :: exxdiv_treatment 
  !
  ! x_gamma_extrapolation
  LOGICAL           :: x_gamma_extrapolation =.TRUE.
  LOGICAl           :: on_double_grid =.FALSE.
  REAL(DP)          :: grid_factor = 8.d0/7.d0 
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

    USE mp_global,  ONLY : me_bgrp, nproc_bgrp, intra_bgrp_comm
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

    IMPLICIT NONE

    IF(ecutfock <= 0.0_DP) ecutfock = ecutrho
    IF(ecutfock < ecutwfc) CALL errore('exx_fft_create', &
            'ecutfock can not be smaller than ecutwfc!', 1) 

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
    IMPLICIT NONE
    !
    IF ( allocated(index_xkq) ) DEALLOCATE(index_xkq)
    IF ( allocated(index_xk ) ) DEALLOCATE(index_xk )
    IF ( allocated(index_sym) ) DEALLOCATE(index_sym)
    IF ( ALLOCATED (rir)       ) DEALLOCATE (rir)
    IF ( allocated(x_occupation) ) DEALLOCATE(x_occupation)
    IF ( allocated(xkq_collect) )  DEALLOCATE(xkq_collect)
    IF ( allocated(exxbuff) )      DEALLOCATE(exxbuff)
    IF ( allocated(exxbuff_nc) )   DEALLOCATE(exxbuff_nc)
    !
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
  !------------------------------------------------------------------------
  SUBROUTINE exx_grid_init()
    !------------------------------------------------------------------------
    !
    USE symm_base,  ONLY : nsym, s
    USE cell_base,  ONLY : bg, at, alat
    USE lsda_mod,   ONLY : nspin
    USE spin_orb,   ONLY : domag
    USE noncollin_module, ONLY : nspin_lsda
    USE klist,      ONLY : xk, wk, nkstot, nks
    USE wvfct,      ONLY : nbnd
    USE io_global,  ONLY : stdout
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

    CALL start_clock ('exx_grid')

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
    IF (npool>1) THEN
      CALL xk_wk_collect(xk_collect, wk_collect, xk, wk, nkstot, nks)
    ELSE
      xk_collect = xk
      wk_collect = wk
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
          ! *** do-loop skipped the first time becasue temp_nksq == 0
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
    WRITE(stdout, '(5x,a,i10)') "EXX: grid of k+q point setup nkqs = ", nkqs

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
    USE cell_base,  ONLY : bg, at, alat
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
    CASE ( "gygi-baldereschi", "gygi-bald", "g-b" )
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
    ! <AF>
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
        !CALL start_clock ('exx_vcut_init')
        CALL vcut_init( vcut, atws, ecutvcut )
        !CALL stop_clock ('exx_vcut_init')
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
    USE symm_base, ONLY : s
    USE cell_base, ONLY : bg, at
    USE lsda_mod,  ONLY : nspin
    USE io_global, ONLY : stdout
    USE klist,     ONLY : nkstot, xk
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
                  CALL errore('exx_grid_check', &
                              'something wrong', 1 )
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
    USE io_global,            ONLY : stdout

    IMPLICIT NONE
    LOGICAL, INTENT(IN) :: l_exx_was_active

    IF (.not. l_exx_was_active ) return ! nothing had happpened yet
    !!
    exx_nwordwfc=2*dffts%nnr
    !iunexx = find_free_unit()
    !CALL diropn(iunexx,'exx', exx_nwordwfc, exst) 
    erfc_scrlen = get_screening_parameter()
    exxdiv = exx_divergence() 
    exxalfa = get_exx_fraction()
    CALL start_exx
    CALL weights()
    CALL exxinit()
    fock0 = exxenergy2()
 
    return
    !------------------------------------------------------------------------
  END SUBROUTINE exx_restart
  !------------------------------------------------------------------------
  !------------------------------------------------------------------------
  SUBROUTINE exxinit()
  !------------------------------------------------------------------------

    !This SUBROUTINE is run before the first H_psi() of each iteration.
    !It saves the wavefunctions for the right density matrix. in real space
    !It saves all the wavefunctions in a single file called prefix.exx
    !
    USE wavefunctions_module, ONLY : evc  
    USE io_files,             ONLY : nwordwfc, iunwfc, iunigk, &
                                     tmp_dir, prefix
    USE io_global,            ONLY : stdout
    USE buffers,              ONLY : get_buffer
    USE gvecs,                ONLY : nls, nlsm, ngms, doublegrid
    USE wvfct,                ONLY : nbnd, npwx, npw, igk, wg, et
    USE control_flags,        ONLY : gamma_only
    USE klist,                ONLY : wk, ngk, nks, nkstot
    USE symm_base,            ONLY : nsym, s, sr, ftau

    USE mp_global,            ONLY : nproc_pool, me_pool, nproc_bgrp, me_bgrp, &
                                     init_index_over_band, inter_bgrp_comm, &
                                     inter_pool_comm
    USE mp,                   ONLY : mp_sum
    USE funct,                ONLY : get_exx_fraction, start_exx, exx_is_active, &
                                     get_screening_parameter 
    USE fft_base,             ONLY : cgather_smooth, cscatter_smooth,&
                                     dffts, cgather_custom, cscatter_custom
    USE fft_interfaces,       ONLY : invfft
    USE uspp,                 ONLY : nkb, okvan

    IMPLICIT NONE
    INTEGER :: ik,ibnd, i, j, k, ir, ri, rj, rk, isym, ikq
    INTEGER :: h_ibnd, half_nbnd
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
       nxxs=exx_fft_g2r%dfftt%nr1x *exx_fft_g2r%dfftt%nr2x *exx_fft_g2r%dfftt%nr3x 
       nrxxs= exx_fft_g2r%dfftt%nnr
       nr1 = exx_fft_g2r%dfftt%nr1
       nr2 = exx_fft_g2r%dfftt%nr2
       nr3 = exx_fft_g2r%dfftt%nr3
       nr1x = exx_fft_g2r%dfftt%nr1x
       nr2x = exx_fft_g2r%dfftt%nr2x
       nr3x = exx_fft_g2r%dfftt%nr3x
    ELSE
       nxxs = dffts%nr1x * dffts%nr2x * dffts%nr3x
       nrxxs= dffts%nnr
       nr1 = dffts%nr1
       nr2 = dffts%nr2
       nr3 = dffts%nr3
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
       IF (.NOT. allocated(exxbuff_nc)) ALLOCATE( exxbuff_nc(nrxxs,npol,nbnd,nkqs))
    ELSE
       ALLOCATE(temppsic(nrxxs), psic(nrxxs))
       if( .not. allocated( exxbuff ) ) ALLOCATE( exxbuff(nrxxs,nbnd,nkqs) )
    ENDIF
    !
    ALLOCATE(ispresent(nsym))
    ALLOCATE(tempevc( npwx*npol, nbnd ))

    IF(.NOT. ALLOCATED(rir)) ALLOCATE(rir(nxxs,nsym))
    rir = 0
    exx_nwordwfc=2*nrxxs
    IF (.not.exx_is_active()) THEN 
       !iunexx = find_free_unit()
       !CALL diropn(iunexx,'exx', exx_nwordwfc, exst) 
       erfc_scrlen = get_screening_parameter()
       exxdiv = exx_divergence() 
       exxalfa = get_exx_fraction()
       !
       CALL start_exx
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
             CALL errore ('exxinit',' EXX smooth grid is not compatible with symmetry: change ecutfock',isym)
          ENDIF
          DO ir=1, nxxs
             rir(ir,isym) = ir
          ENDDO
          DO k = 1, nr3
             DO j = 1, nr2
                DO i = 1, nr1
                   CALL ruotaijk (s(1,1,isym), ftau(1,isym), i, j, k, nr1,nr2,nr3, ri, rj, rk )
                   ir =   i + ( j-1)*nr1x + ( k-1)*nr1x*nr2x
                   rir(ir,isym) = ri + (rj-1)*nr1x + (rk-1)*nr1x*nr2x
                ENDDO
             ENDDO
          ENDDO

       ENDIF
    ENDDO

    IF (noncolin) THEN
       exxbuff_nc=(0.0_DP,0.0_DP)
    ELSE
       exxbuff=(0.0_DP,0.0_DP)
    ENDIF
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
       GAMMA_OR_NOT : & 
       IF (gamma_only) THEN
          half_nbnd = ( nbnd + 1 )/2
          h_ibnd = 0

          do ibnd =1, nbnd, 2     
             h_ibnd = h_ibnd + 1
             !
             temppsic(:) = ( 0._dp, 0._dp )
             !
             if (ibnd < nbnd) then
                temppsic(exx_fft_g2r%nlt(1:exx_fft_g2r%npwt)) = tempevc(1:exx_fft_g2r%npwt,ibnd)  &
                          + ( 0.D0, 1.D0 ) * tempevc(1:exx_fft_g2r%npwt,ibnd+1)
                temppsic(exx_fft_g2r%nltm(1:exx_fft_g2r%npwt)) = CONJG( tempevc(1:exx_fft_g2r%npwt,ibnd) ) &
                          + ( 0.D0, 1.D0 ) * CONJG( tempevc(1:exx_fft_g2r%npwt,ibnd+1) )
             else
                temppsic(exx_fft_g2r%nlt (1:exx_fft_g2r%npwt)) = tempevc(1:exx_fft_g2r%npwt,ibnd) 
                temppsic(exx_fft_g2r%nltm(1:exx_fft_g2r%npwt)) = CONJG( tempevc(1:exx_fft_g2r%npwt,ibnd) ) 
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
                !CALL davcio(psic,exx_nwordwfc,iunexx,(ikq-1)*half_nbnd+h_ibnd,1)
             ENDDO
          END DO
          !
       ELSE GAMMA_OR_NOT 
          !
          IBND_LOOP_K : &
          DO ibnd =1, nbnd     
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
                   exxbuff_nc(:,:,ibnd,ikq)=psic_nc(:,:)
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
                  !CALL davcio(psic,exx_nwordwfc,iunexx,(ikq-1)*nbnd+ibnd,1)
                ENDIF ! noncolinear
             ENDDO
             !
          ENDDO &
          IBND_LOOP_K 
          !
       ENDIF & 
       GAMMA_OR_NOT
    ENDDO &
    KPOINTS_LOOP
    !
    !   All pools must have the complete set of wavefunctions (i.e. from every kpoint)
    IF (npool>1) THEN
       IF (noncolin) THEN
          CALL mp_sum(exxbuff_nc, inter_pool_comm)
       ELSE
          CALL mp_sum(exxbuff, inter_pool_comm)
       END IF
    END IF
    !
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
  SUBROUTINE vexx(lda, n, m, psi, hpsi)
  !-----------------------------------------------------------------------
    !This routine calculates V_xx \Psi
    
    ! ... This routine computes the product of the Hamiltonian
    ! ... matrix with m wavefunctions contained in psi
    !
    ! ... input:
    ! ...    lda   leading dimension of arrays psi, spsi, hpsi
    ! ...    n     true dimension of psi, spsi, hpsi
    ! ...    m     number of states psi
    ! ...    psi
    !
    ! ... output:
    ! ...    hpsi  Vexx*psi
    !
    USE constants,      ONLY : fpi, e2, pi
    USE cell_base,      ONLY : alat, omega, bg, at, tpiba
    USE symm_base,      ONLY : nsym, s
    USE gvect,          ONLY : ngm
    USE gvecs,          ONLY : nls, nlsm, ngms, doublegrid
    USE wvfct,          ONLY : nbnd, npwx, npw, igk, current_k, wg
    USE control_flags,  ONLY : gamma_only
    USE klist,          ONLY : xk, wk, nks, nkstot
    USE lsda_mod,       ONLY : lsda, current_spin, isk
    USE gvect,          ONLY : g, nl
    USE fft_base,       ONLY : dffts
    USE fft_interfaces, ONLY : fwfft, invfft

    USE mp_global,      ONLY : ibnd_start, ibnd_end, inter_bgrp_comm, &
                               intra_bgrp_comm, my_bgrp_id, nbgrp
    USE mp,             ONLY : mp_sum, mp_barrier
    USE gvect,          ONLY : ecutrho
    USE wavefunctions_module, ONLY : psic


    IMPLICIT NONE

    INTEGER                  :: lda, n, m
    COMPLEX(DP)              :: psi(lda*npol,m) 
    COMPLEX(DP)              :: hpsi(lda*npol,m)
    !
    ! local variables
    COMPLEX(DP),ALLOCATABLE :: tempphic(:), temppsic(:), result(:)
    COMPLEX(DP),ALLOCATABLE :: tempphic_nc(:,:), temppsic_nc(:,:), &
                               result_nc(:,:)
    !
    COMPLEX(DP),ALLOCATABLE :: rhoc(:), vc(:), deexx(:)!,:,:)
    REAL(DP),   ALLOCATABLE :: fac(:)
    INTEGER          :: ibnd, ik, im , ikq, iq, isym, ipol
    INTEGER          :: h_ibnd, half_nbnd, nrxxs
    INTEGER          :: current_ik
    INTEGER          :: ibnd_loop_start
    REAL(DP) :: x1, x2, xkp(3)
    REAL(DP) :: xk_cryst(3), sxk(3), xkq(3)
    ! <LMS> temp array for vcut_spheric
    INTEGER  :: find_current_k

    CALL start_clock ('vexx')

    IF(gamma_only) THEN
       ALLOCATE( fac(exx_fft_r2g%ngmt) )
       nrxxs= exx_fft_g2r%dfftt%nnr
    ELSE
       ALLOCATE( fac(ngm) )
       nrxxs = dffts%nnr
    ENDIF

    IF (noncolin) THEN
       ALLOCATE(tempphic_nc(nrxxs,npol), temppsic_nc(nrxxs,npol), result_nc(nrxxs,npol))
    ELSE
       ALLOCATE( tempphic(nrxxs), temppsic(nrxxs), result(nrxxs) )
    ENDIF

    ALLOCATE(rhoc(nrxxs), vc(nrxxs))
    !
    current_ik=find_current_k(current_k,nkstot,nks)
    xkp = xk_collect(:,current_ik)
    !
    ! This is to stop numerical inconsistencies creeping in through the band parallelization.
    !
    IF(my_bgrp_id>0) THEN
       hpsi=(0.0_DP,0.0_DP)
       psi=(0.0_DP,0.0_DP)
    ENDIF
    IF (nbgrp>1) THEN
       CALL mp_sum(hpsi,inter_bgrp_comm)
       CALL mp_sum(psi,inter_bgrp_comm)
    ENDIF
    !
    LOOP_ON_PSI_BANDS : &
    DO im=1,m !for each band of psi (the k cycle is outside band)
      
       IF (noncolin) THEN
          temppsic_nc = ( 0.D0, 0.D0 )
       ELSE
          temppsic(:) = ( 0.D0, 0.D0 )
       ENDIF

       IF(gamma_only) THEN
          !
          temppsic(exx_fft_g2r%nlt(1:exx_fft_g2r%npwt)) =&
               & psi(1:exx_fft_g2r%npwt, im) 
          temppsic(exx_fft_g2r%nltm(1:exx_fft_g2r%npwt)) =&
               & CONJG(psi(1:exx_fft_g2r%npwt,im))
          !
          CALL invfft ('CustomWave', temppsic, exx_fft_g2r%dfftt)
          !
       ELSE
          IF (noncolin) THEN
             temppsic_nc(nls(igk(1:npw)),1) = psi(1:npw,im)
             CALL invfft ('Wave', temppsic_nc(:,1), dffts)
             temppsic_nc(nls(igk(1:npw)),2) = psi(npwx+1:npwx+npw,im)
             CALL invfft ('Wave', temppsic_nc(:,2), dffts)
           ELSE
             temppsic(nls(igk(1:npw))) = psi(1:npw,im)
             CALL invfft ('Wave', temppsic, dffts)
           ENDIF
       ENDIF

      IF (noncolin) THEN
          result_nc(:,:) = (0.0_DP,0.0_DP)
      ELSE
          result(:)   = (0._dp,0._dp)
      ENDIF
        
      INTERNAL_LOOP_ON_Q : &
      DO iq=1,nqs
        !
        ikq  = index_xkq(current_ik,iq)
        ik   = index_xk(ikq)
        isym = ABS(index_sym(ikq))
        xkq = xkq_collect(:,ikq)
        !
        ! calculate the 1/|r-r'| (actually, k+q+g) factor and place it in fac
        IF(gamma_only) THEN
            CALL g2_convolution(exx_fft_r2g%ngmt, exx_fft_r2g%gt, xk(:,current_k), xkq, fac) 
        ELSE
            CALL g2_convolution(ngms, g, xk(:,current_k), xkq, fac)
        ENDIF
        !
        GAMMA_OR_NOT : &
        IF (gamma_only) THEN
            half_nbnd = ( nbnd + 1 ) / 2
            h_ibnd = ibnd_start/2
            IF(MOD(ibnd_start,2)==0) THEN
              h_ibnd=h_ibnd-1
              ibnd_loop_start=ibnd_start-1
            ELSE
              ibnd_loop_start=ibnd_start
            ENDIF

            IBND_LOOP_GAM : &
            DO ibnd=ibnd_loop_start,ibnd_end, 2 !for each band of psi
              h_ibnd = h_ibnd + 1
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
              IF ( ABS(x1) < 1.d-6 .AND.  ABS(x2) < 1.d-6 ) CYCLE
              !
              !loads the phi from file
              tempphic = exxbuff(:,h_ibnd,ikq)
              !CALL davcio ( tempphic, exx_nwordwfc, iunexx, &
              !                    (ikq-1)*half_nbnd+h_ibnd, -1 )
              !calculate rho in real space
              rhoc(:)=(0._dp,0._dp)
              rhoc(1:nrxxs)=CONJG(tempphic(1:nrxxs))*temppsic(1:nrxxs) / omega
              !brings it to G-space
              IF (ecutfock == ecutrho) THEN
                 CALL fwfft ('Custom', rhoc, exx_fft_r2g%dfftt)
                 vc(:) = ( 0.D0, 0.D0 )
                 vc(exx_fft_r2g%nlt(1:exx_fft_r2g%ngmt))  =&
                      & fac(1:exx_fft_r2g%ngmt) * rhoc(exx_fft_r2g&
                      &%nlt(1:exx_fft_r2g%ngmt)) 
                 vc(exx_fft_r2g%nltm(1:exx_fft_r2g%ngmt)) =&
                      & fac(1:exx_fft_r2g%ngmt) * rhoc(exx_fft_r2g&
                      &%nltm(1:exx_fft_r2g%ngmt)) 
                 !brings back v in real space
                 CALL invfft ('Custom', vc, exx_fft_r2g%dfftt) 
              ELSE
                 CALL fwfft ('CustomWave', rhoc, exx_fft_r2g%dfftt)
                 vc(:) = ( 0.D0, 0.D0 )
                 vc(exx_fft_r2g%nlt(1:exx_fft_r2g%npwt))  =&
                      & fac(1:exx_fft_r2g%npwt) * rhoc(exx_fft_r2g&
                      &%nlt(1:exx_fft_r2g%npwt)) 
                 vc(exx_fft_r2g%nltm(1:exx_fft_r2g%npwt)) =&
                      & fac(1:exx_fft_r2g%npwt) * rhoc(exx_fft_r2g&
                      &%nltm(1:exx_fft_r2g%npwt)) 
                 !brings back v in real space
                 CALL invfft ('CustomWave', vc, exx_fft_r2g%dfftt) 
              ENDIF
   
              vc = CMPLX( x1 * DBLE (vc), x2 * AIMAG(vc) ,kind=DP)/ nqs

              !accumulates over bands and k points
              result(1:nrxxs) = result(1:nrxxs) + DBLE( vc(1:nrxxs) &
                   &* tempphic(1:nrxxs))
          END DO &
          IBND_LOOP_GAM
          !
      ELSE GAMMA_OR_NOT
          !
          IBND_LOOP_K : &
          DO ibnd=ibnd_start,ibnd_end !for each band of psi
              IF ( ABS(x_occupation(ibnd,ik)) < eps_occ) CYCLE IBND_LOOP_K
            !
            !loads the phi from file
            !
            IF (noncolin) THEN
                tempphic_nc(:,:)=exxbuff_nc(:,:,ibnd,ikq)
            ELSE
                tempphic(:)=exxbuff(:,ibnd,ikq)
            ENDIF
            !calculate rho in real space
            IF (noncolin) THEN
                rhoc(:) = ( CONJG(tempphic_nc(:,1))*temppsic_nc(:,1) + &
                            CONJG(tempphic_nc(:,2))*temppsic_nc(:,2) )/omega
            ELSE
                rhoc(:)=CONJG(tempphic(:))*temppsic(:) / omega
            ENDIF

            !brings it to G-space
            CALL fwfft('Smooth', rhoc, dffts)

            vc(:) = ( 0._dp, 0._dp )
            vc(nls(1:ngms)) = fac(1:ngms) * rhoc(nls(1:ngms))
            vc = vc * x_occupation(ibnd,ik) / nqs
            !
            !brings back v in real space
            CALL invfft ('Smooth', vc, dffts)
              
            !accumulates over bands and k points
            IF (noncolin) THEN
                DO ipol=1,npol
                  result_nc(:,ipol)= result_nc(:,ipol) &
                                    +vc(:) * tempphic_nc(:,ipol)
                ENDDO
            ELSE
                result(1:nrxxs)=result(1:nrxxs)+vc(1:nrxxs)*tempphic(1:nrxxs)
            END IF
          END DO &
          IBND_LOOP_K
        END IF &
        GAMMA_OR_NOT
        !
      END DO &
      INTERNAL_LOOP_ON_Q
      !
      IF (noncolin) THEN
         CALL mp_sum( result_nc(1:nrxxs,1:npol), inter_bgrp_comm)
      ELSE
         CALL mp_sum( result(1:nrxxs), inter_bgrp_comm)
      END IF
      !
      !brings back result in G-space
      IF( gamma_only) THEN
         !
         CALL fwfft( 'CustomWave' , result, exx_fft_g2r%dfftt )
         !
         hpsi(1:npw,im)=hpsi(1:npw,im) - exxalfa*result(exx_fft_g2r%nlt(1:npw))
         !
      ELSE
          IF (noncolin) THEN
            !brings back result in G-space
            CALL fwfft ('Wave', result_nc(:,1), dffts)
            CALL fwfft ('Wave', result_nc(:,2), dffts)
            !adds it to hpsi
            hpsi(1:n,im)         = hpsi(1:n,im)         - exxalfa*result_nc(nls(igk(1:n)),1)
            hpsi(lda+1:lda+n,im) = hpsi(lda+1:lda+n,im) - exxalfa*result_nc(nls(igk(1:n)),2)
          ELSE
            CALL fwfft ('Wave', result, dffts)
            !adds it to hpsi
            hpsi(1:npw,im)=hpsi(1:npw,im) - exxalfa*result(nls(igk(1:npw)))
          ENDIF
      ENDIF
      !
    END DO &
    LOOP_ON_PSI_BANDS
      
    IF (noncolin) THEN
      DEALLOCATE(tempphic_nc, temppsic_nc, result_nc) 
    ELSE
      DEALLOCATE(tempphic, temppsic, result) 
    END IF

    DEALLOCATE(rhoc, vc, fac )

    !
    CALL stop_clock ('vexx')
    !
    !-----------------------------------------------------------------------
  END SUBROUTINE vexx
  !-----------------------------------------------------------------------
  !
  !-----------------------------------------------------------------------
  SUBROUTINE g2_convolution(ngm, g, xk, xkq, fac)
  !-----------------------------------------------------------------------
    ! This routine calculates the 1/|r-r'| part of the exact exchange 
    ! expression in reciprocal space (the G^-2 factor).
    ! It then regularizes it according to the specified recipe
    USE kinds,     ONLY : DP
    USE cell_base, ONLY : tpiba, at
    USE constants, ONLY : fpi, e2, pi
    
    IMPLICIT NONE
    
    INTEGER,  INTENT(IN)    :: ngm   ! Number of G vectors
    REAL(DP), INTENT(IN)    :: g(3,ngm) ! Cartesian components of G vectors
    REAL(DP), INTENT(IN)    :: xk(3) ! current k vector
    REAL(DP), INTENT(IN)    :: xkq(3) ! current q vector
    
    REAL(DP), INTENT(INOUT) :: fac(ngm) ! Calculated convolution
    
    
    !Local variables
    INTEGER :: ig !Counters 
    REAL(DP) :: q(3), qq, x
    
    !CALL start_clock ('vexx_ngmloop')
    DO ig=1,ngm
      !
      q(:)= xk(:) - xkq(:) + g(:,ig)
      !
      q = q * tpiba
      !
      qq = SUM(q(:)**2) 
      !
      IF (x_gamma_extrapolation) THEN
          on_double_grid = .TRUE.
          x= 0.5d0/tpiba*(q(1)*at(1,1)+q(2)*at(2,1)+q(3)*at(3,1))*nq1
          on_double_grid = on_double_grid .AND. (ABS(x-NINT(x))<eps)
          x= 0.5d0/tpiba*(q(1)*at(1,2)+q(2)*at(2,2)+q(3)*at(3,2))*nq2
          on_double_grid = on_double_grid .AND. (ABS(x-NINT(x))<eps)
          x= 0.5d0/tpiba*(q(1)*at(1,3)+q(2)*at(2,3)+q(3)*at(3,3))*nq3
          on_double_grid = on_double_grid .AND. (ABS(x-NINT(x))<eps)
      ENDIF
      
      IF ( use_coulomb_vcut_ws ) THEN
          !
          fac(ig) = vcut_get(vcut,q)
          !
      ELSE IF ( use_coulomb_vcut_spheric ) THEN
          !
          fac(ig) = vcut_spheric_get(vcut,q)
          !
      ELSE IF (qq > eps_qdiv) THEN
          !
          IF ( erfc_scrlen > 0  ) THEN
            fac(ig)=e2*fpi/qq*(1._dp-EXP(-qq/4.d0/erfc_scrlen**2)) * grid_factor
          ELSEIF( erf_scrlen > 0 ) THEN
            fac(ig)=e2*fpi/qq*(EXP(-qq/4.d0/erf_scrlen**2)) * grid_factor
          ELSE
            fac(ig)=e2*fpi/( qq + yukawa ) * grid_factor
          END IF
          IF (on_double_grid) fac(ig) = 0._dp
          !
      ELSE
          !
          fac(ig)= - exxdiv ! or rather something ELSE (see F.Gygi)
          !
          IF ( yukawa > 0._dp.AND. .NOT. x_gamma_extrapolation ) &
              fac(ig) = fac(ig) + e2*fpi/( qq + yukawa )
          IF( erfc_scrlen > 0._dp.AND. .NOT. x_gamma_extrapolation ) fac(ig) = fac(ig) + e2*pi/(erfc_scrlen**2)
          !
      ENDIF
      !
    ENDDO
    !CALL stop_clock ('vexx_ngmloop')
  END SUBROUTINE g2_convolution
  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  FUNCTION exxenergy ()
    !-----------------------------------------------------------------------
    ! This function is called to correct the deband value and have 
    ! the correct energy 
    USE io_files,   ONLY : iunigk,iunwfc, nwordwfc
    USE buffers,    ONLY : get_buffer
    USE wvfct,      ONLY : nbnd, npwx, npw, igk, wg, current_k
    USE control_flags, ONLY : gamma_only
    USE gvect,      ONLY : gstart
    USE wavefunctions_module, ONLY : evc
    USE lsda_mod,   ONLY : lsda, current_spin, isk
    USE klist,      ONLY : ngk, nks, xk
    USE mp_global,  ONLY : inter_pool_comm, inter_bgrp_comm, intra_bgrp_comm, nbgrp
    USE mp,         ONLY : mp_sum

    IMPLICIT NONE

    REAL(DP)       :: exxenergy,  energy
    INTEGER        :: ibnd, ik
    COMPLEX(DP)    :: vxpsi ( npwx*npol, nbnd ), psi(npwx*npol,nbnd)
    COMPLEX(DP),EXTERNAL :: ZDOTC

    CALL start_clock ('exxenergy')

    energy=0._dp

    IF ( nks > 1 ) REWIND( iunigk )
    DO ik=1,nks
       current_k = ik
       IF ( lsda ) current_spin = isk(ik)
       npw = ngk (ik)
       IF ( nks > 1 ) THEN
          READ( iunigk ) igk
          CALL get_buffer(psi, nwordwfc, iunwfc, ik)
       ELSE
          psi(1:npwx*npol,1:nbnd) = evc(1:npwx*npol,1:nbnd)
       END IF
       !
       vxpsi(:,:) = (0._dp, 0._dp)

       CALL vexx(npwx,npw,nbnd,psi,vxpsi)

       DO ibnd=1,nbnd
          energy = energy + wg(ibnd,ik) * ZDOTC(npw,psi(1,ibnd),1,vxpsi(1,ibnd),1)
          IF (noncolin) energy = energy + &
                            wg(ibnd,ik) * ZDOTC(npw,psi(npwx+1,ibnd),1,vxpsi(npwx+1,ibnd),1)
       ENDDO
       IF (gamma_only .and. gstart == 2) THEN
           DO ibnd=1,nbnd
              energy = energy - &
                       0.5_dp * wg(ibnd,ik) * CONJG(psi(1,ibnd)) * vxpsi(1,ibnd)
           ENDDO
       ENDIF
    END DO
    !
    IF (gamma_only) energy = 2 * energy

    CALL mp_sum( energy, intra_bgrp_comm)
    CALL mp_sum( energy, inter_pool_comm )
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
    USE constants, ONLY : fpi, e2, pi
    USE io_files,  ONLY : iunigk,iunwfc, nwordwfc
    USE buffers,   ONLY : get_buffer
    USE cell_base, ONLY : alat, omega, bg, at, tpiba
    USE symm_base, ONLY : nsym, s
    USE gvect,     ONLY : ngm, gstart, g, nl
    USE gvecs,     ONLY : ngms, nls, nlsm, doublegrid
    USE wvfct,     ONLY : nbnd, npwx, npw, igk, wg
    USE control_flags,        ONLY : gamma_only
    USE wavefunctions_module, ONLY : evc
    USE klist,     ONLY : xk, ngk, nks, nkstot
    USE lsda_mod,  ONLY : lsda, current_spin, isk
    USE mp_global, ONLY : inter_pool_comm, inter_image_comm, inter_bgrp_comm, intra_bgrp_comm, nbgrp
    USE mp_global, ONLY : my_image_id, nimage, ibnd_start, ibnd_end
    USE mp,        ONLY : mp_sum
    USE fft_base,  ONLY : dffts
    USE fft_interfaces,       ONLY : fwfft, invfft
    USE gvect,     ONLY : ecutrho
    USE klist,     ONLY : wk

    IMPLICIT NONE
    !
    REAL(DP)   :: exxenergy2,  energy
    !
    ! local variables
    COMPLEX(DP), allocatable :: tempphic(:), temppsic(:)
    COMPLEX(DP), ALLOCATABLE :: tempphic_nc(:,:), temppsic_nc(:,:)
    COMPLEX(DP), ALLOCATABLE :: rhoc(:)
    REAL(DP),   ALLOCATABLE :: fac(:)
    INTEGER          :: jbnd, ibnd, ik, ikk, ig, ikq, iq, isym
    INTEGER          :: half_nbnd, h_ibnd, nrxxs, current_ik
    REAL(DP)    :: x1, x2
    REAL(DP) :: xk_cryst(3), sxk(3), xkq(3), vc
    ! temp array for vcut_spheric
    INTEGER,EXTERNAL :: find_current_k

    COMPLEX(kind=DP), ALLOCATABLE :: psi_t(:), prod_tot(:)
    INTEGER, ALLOCATABLE :: igkt(:)
    !
    CALL start_clock ('exxen2')

    IF(gamma_only) THEN
       nrxxs= exx_fft_g2r%dfftt%nnr
       ALLOCATE( fac(exx_fft_r2g%ngmt) )
    ELSE
       nrxxs = dffts%nnr
       ALLOCATE( fac(ngms) )
    ENDIF

    IF (noncolin) THEN
       ALLOCATE(tempphic_nc(nrxxs,npol), temppsic_nc(nrxxs,npol))
    ELSE
       ALLOCATE(tempphic(nrxxs), temppsic(nrxxs)) 
    ENDIF
    ALLOCATE( rhoc(nrxxs) )

    energy=0._dp

    IF ( nks > 1 ) REWIND( iunigk )

    IKK_LOOP : &
    DO ikk=1,nks
       current_ik=find_current_k(ikk,nkstot,nks)
       !
       IF ( lsda ) current_spin = isk(ikk)
       npw = ngk (ikk)
       IF ( nks > 1 ) THEN
          READ( iunigk ) igk
          CALL get_buffer (evc, nwordwfc, iunwfc, ikk)
       END IF
       !
       JBND_LOOP : &
       DO jbnd = ibnd_start, ibnd_end !for each band of psi (the k cycle is outside band)
          IF (noncolin) THEN
              temppsic_nc = ( 0._dp, 0._dp )
          ELSE
              temppsic = ( 0._dp, 0._dp )
          ENDIF
          IF(gamma_only) THEN
             !
             temppsic(exx_fft_g2r%nlt(1:exx_fft_g2r%npwt)) =&
                  & evc(1:exx_fft_g2r%npwt,jbnd) 
             temppsic(exx_fft_g2r%nltm(1:exx_fft_g2r%npwt)) =&
                  & CONJG(evc(1:exx_fft_g2r%npwt,jbnd))
             CALL invfft ('CustomWave', temppsic, exx_fft_g2r%dfftt)
             !
          ELSE
              IF (noncolin) THEN
                temppsic_nc(nls(igk(1:npw)),1) = evc(1:npw,jbnd)
                CALL invfft ('Wave', temppsic_nc(:,1), dffts)
                temppsic_nc(nls(igk(1:npw)),2) = evc(npwx+1:npwx+npw,jbnd)
                CALL invfft ('Wave', temppsic_nc(:,2), dffts)
              ELSE
                temppsic(nls(igk(1:npw))) = evc(1:npw,jbnd)
                CALL invfft ('Wave', temppsic, dffts)
              ENDIF
          ENDIF
          
          IQ_LOOP : &
          DO iq = 1,nqs
            !
            ikq  = index_xkq(current_ik,iq)
            ik   = index_xk(ikq)
            isym = abs(index_sym(ikq))
            !
            xkq = xkq_collect(:,ikq)
            !
            IF(gamma_only) THEN
              CALL g2_convolution(exx_fft_r2g%ngmt, exx_fft_r2g%gt, xk(:,current_ik), xkq, fac) 
              fac(exx_fft_r2g%gstart_t:) = 2.d0 * fac(exx_fft_r2g%gstart_t:)
            ELSE
              CALL g2_convolution(ngms, g, xk(:,current_ik), xkq, fac)
            ENDIF

            GAMMA_OR_NOT : &
            IF (gamma_only) THEN
              half_nbnd = ( nbnd + 1) / 2
              h_ibnd = 0
              !
              IBND_LOOP_GAM : &
              DO ibnd=1,nbnd,2 !for each band of psi
                  h_ibnd = h_ibnd + 1
                  x1 = x_occupation(ibnd,ik)
                  IF ( ibnd < nbnd ) THEN
                    x2 = x_occupation(ibnd+1,ik)
                  ELSE
                    x2 = 0._dp
                  ENDIF
                  IF ( abs(x1) < 1.d-6 .and. abs(x2) < 1.d-6 ) cycle
                  !
                  !loads the phi from file
                  !
                  tempphic(1:nrxxs)=exxbuff(1:nrxxs,h_ibnd,ikq)
                  !
                  !CALL davcio (tempphic, exx_nwordwfc, iunexx, &
                  !                       (ikq-1)*half_nbnd+h_ibnd, -1 )
                  !calculate rho in real space
                  rhoc(:)=(0._dp, 0._dp)
                  rhoc(1:nrxxs)=CONJG(tempphic(1:nrxxs))*temppsic(1:nrxxs) / omega
                  IF_ECUTFOCK : &
                  IF (ecutfock == ecutrho) THEN
                    !brings it to G-space
                    CALL fwfft ('Custom', rhoc, exx_fft_r2g%dfftt)
                    vc = 0._dp
                    DO ig=1,exx_fft_r2g%ngmt
                       vc = vc + fac(ig) * x1 * &
                            ABS( rhoc(exx_fft_r2g%nlt(ig)) + CONJG(rhoc(exx_fft_r2g%nltm(ig))) )**2
                       vc = vc + fac(ig) * x2 * &
                            ABS( rhoc(exx_fft_r2g%nlt(ig)) - CONJG(rhoc(exx_fft_r2g%nltm(ig))) )**2
                    END DO
                    !
                  ELSE IF_ECUTFOCK
                    !
                    !brings it to G-space
                    CALL fwfft ('CustomWave', rhoc, exx_fft_r2g%dfftt)
                    vc = 0._dp
                    DO ig=1,exx_fft_r2g%npwt
                       vc = vc + fac(ig) * x1 * &
                            ABS( rhoc(exx_fft_r2g%nlt(ig)) + CONJG(rhoc(exx_fft_r2g%nltm(ig))) )**2
                       vc = vc + fac(ig) * x2 * &
                            ABS( rhoc(exx_fft_r2g%nlt(ig)) - CONJG(rhoc(exx_fft_r2g%nltm(ig))) )**2
                    END DO
                  ENDIF&
                  IF_ECUTFOCK
                  !
                  vc = vc * omega * 0.25d0 / nqs
                  energy = energy - exxalfa * vc * wg(jbnd,ikk)
                  !
                END DO &
                IBND_LOOP_GAM
                !
             ELSE GAMMA_OR_NOT
                !
                IBND_LOOP_K : &
                DO ibnd=1,nbnd !for each band of psi
                   IF ( abs(x_occupation(ibnd,ik)) < 1.d-6) cycle
                   !
                   !loads the phi from file
                   !
                   IF (noncolin) THEN
                      tempphic_nc(:,:)=exxbuff_nc(:,:,ibnd,ikq)
                      rhoc(:)=(CONJG(tempphic_nc(:,1))*temppsic_nc(:,1) + &
                               CONJG(tempphic_nc(:,2))*temppsic_nc(:,2) )/omega
                   ELSE
                      tempphic(:)=exxbuff(:,ibnd,ikq)

                      !CALL davcio (tempphic, exx_nwordwfc, iunexx, &
                      !                       (ikq-1)*nbnd+ibnd, -1 )
                      !calculate rho in real space
                      rhoc(:)=CONJG(tempphic(:))*temppsic(:) / omega
                   ENDIF
                   !brings it to G-space
                   CALL fwfft ('Smooth', rhoc, dffts)
                   !
                   vc = 0._dp
                   DO ig=1,ngms
                      vc = vc + fac(ig) * rhoc(nls(ig)) * CONJG(rhoc(nls(ig)))
                   ENDDO
                   vc = vc * omega * x_occupation(ibnd,ik) / nqs
 
                   energy = energy - exxalfa * vc * wg(jbnd,ikk)
              END DO &
              IBND_LOOP_K 
            END IF  &
            GAMMA_OR_NOT
            !
          END DO IQ_LOOP
       END DO JBND_LOOP
    END DO IKK_LOOP

    IF (noncolin) THEN
       DEALLOCATE(tempphic_nc, temppsic_nc) 
    ELSE
       DEALLOCATE(tempphic, temppsic) 
    ENDIF

    DEALLOCATE(rhoc, fac )

!
! Was used for image parallelization
!    CALL mp_sum( energy, inter_image_comm )
!
    CALL mp_sum( energy, inter_bgrp_comm )
    CALL mp_sum( energy, intra_bgrp_comm )
    CALL mp_sum( energy, inter_pool_comm )
    !
    exxenergy2 = energy
    !
    CALL stop_clock ('exxen2')

    !-----------------------------------------------------------------------
  END FUNCTION  exxenergy2
  !-----------------------------------------------------------------------
  !
  !-----------------------------------------------------------------------
  FUNCTION exx_divergence ()
    !-----------------------------------------------------------------------
     USE constants, ONLY : fpi, e2, pi
     USE cell_base, ONLY : bg, at, alat, omega
     USE gvect,     ONLY : ngm, g
     USE wvfct,     ONLY : ecutwfc
     USE io_global, ONLY : stdout
     USE control_flags, ONLY : gamma_only
     USE mp_global, ONLY : intra_bgrp_comm
     USE mp,        ONLY : mp_sum

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
        return
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

!     div = div - e2*omega/sqrt(alpha*0.25d0*fpi)

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
    USE constants, ONLY : fpi, e2, pi, tpi
    USE io_files,  ONLY : iunigk,iunwfc, nwordwfc
    USE buffers,   ONLY : get_buffer
    USE cell_base, ONLY : alat, omega, bg, at, tpiba
    USE symm_base,ONLY : nsym, s
    USE gvect,     ONLY : ngm
    USE gvecs,   ONLY : nls, nlsm, doublegrid
    USE wvfct,     ONLY : nbnd, npwx, npw, igk, wg, current_k
    USE control_flags, ONLY : gamma_only
    USE wavefunctions_module, ONLY : evc
    USE klist,     ONLY : xk, ngk, nks
    USE lsda_mod,  ONLY : lsda, current_spin, isk
    USE gvect,     ONLY : g, nl
    USE mp_global, ONLY : inter_pool_comm, inter_bgrp_comm, intra_bgrp_comm
    USE mp_global, ONLY : my_image_id, nimage
    USE mp,        ONLY : mp_sum 
    USE fft_base,  ONLY : dffts
    USE fft_interfaces, ONLY : fwfft, invfft
    ! ---- local variables -------------------------------------------------
    IMPLICIT NONE
    REAL(DP)   :: exx_stress(3,3), exx_stress_(3,3)
    complex(dp), allocatable :: tempphic(:), temppsic(:)
    complex(dp), allocatable :: rhoc(:)
    REAL(DP), allocatable :: fac(:), fac_tens(:,:,:), fac_stress(:)
    INTEGER :: jbnd, ibnd, ik, ikk, ig, ikq, iq, isym
    INTEGER :: half_nbnd, h_ibnd, nqi, iqi, beta, nrxxs
    REAL(DP) :: x1, x2
    REAL(DP) :: qq, xk_cryst(3), sxk(3), xkq(3), vc(3,3), x, q(3)
    ! temp array for vcut_spheric
    REAL(DP) :: delta(3,3)
  
    CALL start_clock ('exx_stress')

    IF (npool>1) CALL errore('exx_stress','stress not available with pools',1)
    IF (noncolin) CALL errore('exx_stress','stress not available with noncolin',1)

    nrxxs = dffts%nnr
    delta = reshape( (/1._dp,0._dp,0._dp, 0._dp,1._dp,0._dp, 0._dp,0._dp,1._dp/), (/3,3/))
    exx_stress_ = 0._dp
    allocate( tempphic(nrxxs), temppsic(nrxxs), rhoc(nrxxs), fac(ngm) )
    allocate( fac_tens(3,3,ngm), fac_stress(ngm) )

    IF ( nks > 1 ) rewind( iunigk )
  !
  ! Was used for image parallelization
  !
    nqi=nqs
  !  nqi = nqs/nimage
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
            temppsic(:) = ( 0._dp, 0._dp )
            temppsic(nls(igk(1:npw))) = evc(1:npw,jbnd)
            if(gamma_only) temppsic(nlsm(igk(1:npw))) = conjg(evc(1:npw,jbnd))
            CALL invfft ('Wave', temppsic, dffts)       

            DO iqi = 1, nqi
  !
  ! Was used for image parallelization
  !
  !              iq = iqi + nqi*my_image_id
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
                DO ig = 1, ngm
                  q(1)= xk(1,current_k) - xkq(1) + g(1,ig)
                  q(2)= xk(2,current_k) - xkq(2) + g(2,ig)
                  q(3)= xk(3,current_k) - xkq(3) + g(3,ig)

                  q = q * tpiba
                  qq = ( q(1)*q(1) + q(2)*q(2) + q(3)*q(3) )

                  DO beta = 1, 3
                      fac_tens(1:3,beta,ig) = q(1:3)*q(beta)
                  enddo

                  IF (x_gamma_extrapolation) THEN
                      on_double_grid = .true.
                      x= 0.5d0/tpiba*(q(1)*at(1,1)+q(2)*at(2,1)+q(3)*at(3,1))*nq1
                      on_double_grid = on_double_grid .and. (abs(x-nint(x))<eps)
                      x= 0.5d0/tpiba*(q(1)*at(1,2)+q(2)*at(2,2)+q(3)*at(3,2))*nq2
                      on_double_grid = on_double_grid .and. (abs(x-nint(x))<eps)
                      x= 0.5d0/tpiba*(q(1)*at(1,3)+q(2)*at(2,3)+q(3)*at(3,3))*nq3
                      on_double_grid = on_double_grid .and. (abs(x-nint(x))<eps)
                  ENDIF

                  IF (use_coulomb_vcut_ws) THEN
                      fac(ig) = vcut_get(vcut, q)
                      fac_stress(ig) = 0._dp   ! not implemented
                      IF (gamma_only .and. qq > 1.d-8) fac(ig) = 2.d0 * fac(ig)

                  ELSE IF ( use_coulomb_vcut_spheric ) THEN
                      fac(ig) = vcut_spheric_get(vcut, q)
                      fac_stress(ig) = 0._dp   ! not implemented
                      IF (gamma_only .and. qq > 1.d-8) fac(ig) = 2.d0 * fac(ig) 

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
                enddo
                !CALL stop_clock ('exxen2_ngmloop')

                IF (gamma_only) THEN
                    half_nbnd = (nbnd + 1) / 2
                    h_ibnd = 0
                    DO ibnd=1,nbnd, 2 !for each band of psi
                        h_ibnd = h_ibnd + 1
                        x1 = x_occupation(ibnd,ik)
                        IF ( ibnd < nbnd ) THEN
                          x2 = x_occupation(ibnd+1,ik)
                        ELSE
                          x2 = 0._dp
                        ENDIF
                        IF ( abs(x1) < 1.d-6 .and. abs(x2) < 1.d-6 ) cycle

                        ! loads the phi from file
                        tempphic(1:nrxxs)=exxbuff(1:nrxxs,h_ibnd,ikq)
                  
                        ! calculate rho in real space
                        rhoc(:)=CONJG(tempphic(:))*temppsic(:) / omega
                        ! brings it to G-space
                        CALL fwfft ('Smooth', rhoc, dffts)
    
                        vc = 0._dp
                        DO ig=1,ngm
                          vc(:,:) = vc(:,:) + fac(ig) * x1 * &
                                    abs( rhoc(nls(ig))+CONJG(rhoc(nlsm(ig))))**2 * &
                                    (fac_tens(:,:,ig)*fac_stress(ig)/2.d0 - delta(:,:)*fac(ig))
                          vc(:,:) = vc(:,:) + fac(ig) * x2 * &
                                    abs( rhoc(nls(ig))-CONJG(rhoc(nlsm(ig))))**2 * &
                                    (fac_tens(:,:,ig)*fac_stress(ig)/2.d0 - delta(:,:)*fac(ig))
                        enddo
                        vc = vc / nqs / 4.d0
                        exx_stress_ = exx_stress_ + exxalfa * vc * wg(jbnd,ikk)
                    enddo
                ELSE
                    DO ibnd=1,nbnd !for each band of psi
                      IF ( abs(x_occupation(ibnd,ik)) < 1.d-6) cycle

                      ! loads the phi from file
                      tempphic(1:nrxxs)=exxbuff(1:nrxxs,ibnd,ikq)
                      !
                      ! calculate rho in real space
                      rhoc(:)=CONJG(tempphic(:))*temppsic(:) / omega

                      ! brings it to G-space
                      CALL fwfft ('Smooth', rhoc, dffts)

                      vc = 0._dp
                      DO ig=1,ngm
                        vc(:,:) = vc(:,:) + rhoc(nls(ig))*CONJG(rhoc(nls(ig))) * &
                                  (fac_tens(:,:,ig)*fac_stress(ig)/2.d0 - delta(:,:)*fac(ig))
                      ENDDO
                      vc = vc * x_occupation(ibnd,ik) / nqs / 4.d0
                      exx_stress_ = exx_stress_ + exxalfa * vc * wg(jbnd,ikk)
                    ENDDO
                ENDIF ! gamma or k-points

            enddo ! iqi
        enddo ! jbnd
    enddo ! ikk

    DEALLOCATE(tempphic, temppsic, rhoc, fac, fac_tens, fac_stress )
  !
  ! Was used for image parallelization
  !  CALL mp_sum( exx_stress_, inter_image_comm )
  !
    CALL mp_sum( exx_stress_, intra_bgrp_comm )
    CALL mp_sum( exx_stress_, inter_pool_comm )
    exx_stress = exx_stress_

    CALL stop_clock ('exx_stress')
    !-----------------------------------------------------------------------
  END FUNCTION exx_stress
  !-----------------------------------------------------------------------
  !
!-----------------------------------------------------------------------
END MODULE exx
!-----------------------------------------------------------------------
