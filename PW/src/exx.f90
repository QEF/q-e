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
  USE io_global,            ONLY : ionode, stdout
  USE fft_custom,           ONLY : fft_cus
  !
  IMPLICIT NONE
  SAVE

  !
  ! general purpose vars
  !
  REAL(DP):: exxalfa=0.d0                ! 1 if exx, 0 elsewhere
  INTEGER :: exx_nwordwfc, ji

  !
  ! variables defining the auxiliary k-point grid 
  ! used in X BZ integration
  !
  INTEGER :: nq1=1, nq2=1, nq3=1         ! integers defining the X integration mesh
  INTEGER :: nqs=1                       ! number of points in the q-gridd
  INTEGER :: nkqs                        ! total number of different k+q
  !
  REAL(DP),    ALLOCATABLE :: xkq(:,:)   ! xkq(3,nkqs) the auxiliary k+q set
  REAL(DP),    ALLOCATABLE :: x_occupation(:,:)           
                                         ! x_occupation(nbnd,nks) the weight of 
                                         ! auxiliary functions in the density matrix
  COMPLEX(DP), ALLOCATABLE :: exxbuff(:,:,:)
                                         ! temporay buffer to store wfc 

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
!
!  Used for k points pool parallelization. All pools needs these quantities.
!  They are allocated only if needed.
!
  REAL(DP),    ALLOCATABLE :: xk_collect(:,:)
  REAL(DP),    ALLOCATABLE :: wk_collect(:)
  REAL(DP),    ALLOCATABLE :: wg_collect(:,:)
  LOGICAL :: pool_para=.FALSE.
  !
  ! variables to deal with Coulomb divergence
  ! and related issues
  !
  REAL (DP)         :: eps =1.d-6
  REAL (DP)         :: exxdiv = 0.d0
  CHARACTER(80)     :: exxdiv_treatment 
  !
  ! x_gamma_extrapolation
  LOGICAL           :: x_gamma_extrapolation =.TRUE.
  LOGICAl           :: on_double_grid =.FALSE.
  REAL (DP)         :: grid_factor = 8.d0/7.d0 
  !
  ! Gygi-Baldereschi 
  LOGICAL           :: use_regularization = .TRUE.
  !
  ! yukawa method
  REAL (DP)         :: yukawa = 0.d0
  !
  ! erfc screening
  REAL (DP)         :: erfc_scrlen = 0.d0
  !
  ! erf screening
  REAL (DP)         :: erf_scrlen = 0.d0
  ! cutoff techniques
  LOGICAL           :: use_coulomb_vcut_ws = .FALSE.
  LOGICAL           :: use_coulomb_vcut_spheric = .FALSE.
  REAL (DP)         :: ecutvcut
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

    USE io_global,  ONLY : ionode_id
    USE mp_global,  ONLY : mpime, nproc, intra_bgrp_comm
    USE mp_wave,    ONLY : mergewf, splitwf
    USE gvect, ONLY : ig_l2g

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: npw
    COMPLEX(kind=DP), INTENT(IN) :: psi(npw)
    COMPLEX(kind=DP), INTENT(INOUT) :: psi_t(:)
    INTEGER, OPTIONAL, INTENT(INOUT) :: igkt(:)
    INTEGER, INTENT(IN) :: sign
    TYPE(fft_cus), INTENT(IN) :: fft
    
    COMPLEX(kind=DP), ALLOCATABLE :: evc_g(:)
    INTEGER :: ig
    
    ALLOCATE( evc_g( fft%ngmt_g ) )
    
    IF(sign > 0 .AND. PRESENT(igkt) ) THEN
       DO ig=1, fft%ngmt
          igkt(ig)=ig
       ENDDO
    ENDIF
    
    IF( fft%dual_t==4.d0) THEN
       psi_t(1:fft%npwt)=psi(1:fft%npwt)
    ELSE
       IF (sign > 0 ) THEN
          CALL mergewf(psi, evc_g, npw, ig_l2g, mpime, nproc,&
               & ionode_id, intra_bgrp_comm)  
          CALL splitwf(psi_t(:), evc_g, fft%npwt, fft%ig_l2gt, mpime,&
               & nproc, ionode_id, intra_bgrp_comm)  
       ELSE
          CALL mergewf(psi, evc_g, fft%npwt, fft%ig_l2gt, mpime,&
               & nproc, ionode_id, intra_bgrp_comm)  
          CALL splitwf(psi_t, evc_g, npw, ig_l2g, mpime, nproc,&
               & ionode_id, intra_bgrp_comm)  
       ENDIF
    ENDIF

    DEALLOCATE( evc_g ) 

  END SUBROUTINE exx_grid_convert
  !------------------------------------------------------------------------
  SUBROUTINE exx_fft_create ()
  !------------------------------------------------------------------------
    
    USE wvfct,        ONLY : ecutwfc
    
    IMPLICIT NONE

    IF(ecutfock < 0.0_DP) ecutfock = ecutwfc

    exx_fft_g2r%ecutt=ecutwfc
    exx_fft_g2r%dual_t=exx_dual
    CALL allocate_fft_custom(exx_fft_g2r)

    exx_fft_r2g%ecutt=ecutfock
    exx_fft_r2g%dual_t=ecutwfc*exx_dual/ecutfock
    CALL allocate_fft_custom(exx_fft_r2g)

  END SUBROUTINE exx_fft_create
  !------------------------------------------------------------------------
  SUBROUTINE exx_fft_destroy ()
  !------------------------------------------------------------------------
    USE fft_custom,  ONLY : deallocate_fft_custom

    IMPLICIT NONE

    CALL deallocate_fft_custom(exx_fft_g2r)
    CALL deallocate_fft_custom(exx_fft_r2g)

  END SUBROUTINE exx_fft_destroy
  !------------------------------------------------------------------------
  SUBROUTINE deallocate_exx ()
  !------------------------------------------------------------------------
  !
  IF ( ALLOCATED (index_xkq) ) DEALLOCATE (index_xkq)
  IF ( ALLOCATED (index_xk ) ) DEALLOCATE (index_xk )
  IF ( ALLOCATED (index_sym) ) DEALLOCATE (index_sym)
  IF ( ALLOCATED (x_occupation) ) DEALLOCATE (x_occupation)
  IF ( ALLOCATED (xkq) ) DEALLOCATE (xkq)
  IF ( ALLOCATED (exxbuff) ) DEALLOCATE (exxbuff)
  !
  CALL exx_fft_destroy()
  !
  !  Pool variables deallocation
  !
  IF ( ALLOCATED (xk_collect) )  DEALLOCATE ( xk_collect )
  IF ( ALLOCATED (wk_collect) )  DEALLOCATE ( wk_collect )
  IF ( ALLOCATED (wg_collect) )  DEALLOCATE ( wg_collect )
  !
  !
  END SUBROUTINE deallocate_exx
  !------------------------------------------------------------------------
  subroutine exx_grid_init()
  !------------------------------------------------------------------------
  !
  USE symm_base,  ONLY : nsym, s
  USE cell_base,  ONLY : bg, at, alat
  USE lsda_mod,   ONLY : nspin
  USE noncollin_module, ONLY : nspin_lsda
  USE klist,      ONLY : xk, wk, nkstot, nks
  USE wvfct,      ONLY : nbnd
  USE io_global,  ONLY : stdout
  !
  USE mp_global,  ONLY : nproc, npool, nimage
  !
  IMPLICIT NONE
  !
  CHARACTER(13) :: sub_name='exx_grid_init'
  integer       :: iq1, iq2, iq3, isym, ik, ikq, iq, max_nk, temp_nkqs
  integer,   allocatable :: temp_index_xk(:), temp_index_sym(:)
  integer,   allocatable :: temp_index_ikq(:), new_ikq(:)
  real (DP), allocatable :: temp_xkq(:,:)
  logical       :: xk_not_found
  real (DP)     :: sxk(3), dxk(3), xk_cryst(3)
  real (DP)     :: dq1, dq2, dq3
  logical       :: no_pool_para
  integer       :: find_current_k

  CALL start_clock ('exx_grid')

  !
  ! definitions and checks
  !
#ifdef EXXDEBUG
  IF (ionode) WRITE(stdout,'(/,2x,a,3i4)') "EXX : q-grid dimensions are ", nq1,nq2,nq3
#endif
  !
  grid_factor = 1.d0
  !
  IF (x_gamma_extrapolation) THEN
      !
#ifdef EXXDEBUG
      IF (ionode) WRITE (stdout,'(2x,a)') "EXX : q->0 dealt with 8/7 -1/7 trick"
#endif
      grid_factor = 8.d0/7.d0
      !
  ENDIF
  

  !
  nqs = nq1 * nq2 * nq3
  !
  ! all processors need to have access to all k+q points
  !
  pool_para =  npool>1
  if (pool_para ) then
     IF ( .NOT.ALLOCATED (xk_collect) )  ALLOCATE (xk_collect(3,nkstot))
     IF ( .NOT.ALLOCATED (wk_collect) )  ALLOCATE (wk_collect(nkstot))
     CALL xk_wk_collect(xk_collect, wk_collect, xk, wk, nkstot, nks)
  end if
  !
  ! set a safe limit as the maximum number of auxiliary points we may need
  ! and allocate auxiliary arrays
  !
  max_nk = nkstot * min(48, 2 * nsym)
  allocate ( temp_index_xk(max_nk), temp_index_sym(max_nk) )
  allocate ( temp_index_ikq(max_nk), new_ikq(max_nk) )
  allocate ( temp_xkq(3,max_nk) )
  !
  ! find all k-points equivalent by symmetry to the points in the k-list
  !
  temp_nkqs = 0
  do isym=1,nsym
     do ik =1, nkstot
        IF (pool_para) THEN
           xk_cryst(:) = at(1,:)*xk_collect(1,ik) + at(2,:)*xk_collect(2,ik)&
                       + at(3,:)*xk_collect(3,ik)
        ELSE
           xk_cryst(:) = at(1,:)*xk(1,ik) + at(2,:)*xk(2,ik) &
                       + at(3,:)*xk(3,ik)
        ENDIF
        sxk(:) = s(:,1,isym)*xk_cryst(1) + &
                 s(:,2,isym)*xk_cryst(2) + &
                 s(:,3,isym)*xk_cryst(3)
        ! add sxk to the auxiliary list if it is not already present
        xk_not_found = .true.
        do ikq=1, temp_nkqs
           if (xk_not_found ) then
              dxk(:) = sxk(:)-temp_xkq(:,ikq) - nint(sxk(:)-temp_xkq(:,ikq))
              if ( abs(dxk(1)).le.eps .and. &
                   abs(dxk(2)).le.eps .and. &
                   abs(dxk(3)).le.eps ) xk_not_found = .false.
           end if
        end do
        if (xk_not_found) then
           temp_nkqs                 = temp_nkqs + 1
           temp_xkq(:,temp_nkqs)     = sxk(:)
           temp_index_xk(temp_nkqs)  = ik
           temp_index_sym(temp_nkqs) = isym 
        end if

        sxk(:) = - sxk(:)
        xk_not_found = .true.
        do ikq=1, temp_nkqs
           if (xk_not_found ) then
              dxk(:) = sxk(:) - temp_xkq(:,ikq) - nint(sxk(:) - temp_xkq(:,ikq))
              if ( abs(dxk(1)).le.eps .and. &
                   abs(dxk(2)).le.eps .and. &
                   abs(dxk(3)).le.eps ) xk_not_found = .false.
           end if
        end do
        if (xk_not_found) then
           temp_nkqs                 = temp_nkqs + 1
           temp_xkq(:,temp_nkqs)     = sxk(:)
           temp_index_xk(temp_nkqs)  = ik
           temp_index_sym(temp_nkqs) =-isym 
        end if

     end do
  end do

  !
  ! define the q-mesh step-sizes
  !
  dq1= 1.d0/DBLE(nq1)
  dq2= 1.d0/DBLE(nq2)
  dq3= 1.d0/DBLE(nq3)
  !
  ! allocate and fill the array index_xkq(nkstot,nqs)
  !
  if(.not.ALLOCATED(index_xkq)) allocate ( index_xkq(nkstot,nqs) )
  if(.not.ALLOCATED(x_occupation)) allocate ( x_occupation(nbnd,nkstot) )
  nkqs = 0
  new_ikq(:) = 0
  do ik=1,nkstot 
     IF (pool_para) THEN
        xk_cryst(:) = at(1,:)*xk_collect(1,ik) + at(2,:)*xk_collect(2,ik)&
                    + at(3,:)*xk_collect(3,ik)
     ELSE
        xk_cryst(:) = at(1,:)*xk(1,ik) + at(2,:)*xk(2,ik) + at(3,:)*xk(3,ik)
     ENDIF

     iq = 0
     do iq1=1, nq1
       sxk(1) = xk_cryst(1) + (iq1-1) * dq1
       do iq2 =1, nq2
         sxk(2) = xk_cryst(2) + (iq2-1) * dq2
         do iq3 =1, nq3
            sxk(3) = xk_cryst(3) + (iq3-1) * dq3
            iq = iq + 1

            xk_not_found = .true.
            do ikq=1, temp_nkqs
               if ( xk_not_found ) then
                  dxk(:) = sxk(:)-temp_xkq(:,ikq) - nint(sxk(:)-temp_xkq(:,ikq))
                  if ( abs(dxk(1)).le.eps .and. &
                       abs(dxk(2)).le.eps .and. &
                       abs(dxk(3)).le.eps ) then

                       xk_not_found = .false.

                       if( new_ikq(ikq) == 0) then
                           nkqs = nkqs + 1
                           temp_index_ikq(nkqs) = ikq
                           new_ikq(ikq) = nkqs
                       end if
                       index_xkq(ik,iq) = new_ikq(ikq)

                  end if
               end if
            end do
            if (xk_not_found) then
               write (*,*) ik, iq, temp_nkqs
               write (*,*) sxk(:)
               call errore('exx_grid_init', &
                           ' k + q is not an S*k ', (ik-1) * nqs + iq )
            end if

         end do
       end do
     end do

  end do
  !
  ! allocate and fill the arrays xkq(3,nkqs), index_xk(nkqs) and index_sym(nkqs)
  !
  allocate ( xkq(3,nspin_lsda*nkqs), index_xk(nspin_lsda*nkqs),  &
             index_sym(nspin_lsda*nkqs) )

  do ik =1, nkqs
     ikq = temp_index_ikq(ik)
     xkq(:,ik) = bg(:,1)*temp_xkq(1,ikq) + &
                 bg(:,2)*temp_xkq(2,ikq) + &
                 bg(:,3)*temp_xkq(3,ikq)
     index_xk(ik)  = temp_index_xk(ikq)
     index_sym(ik) = temp_index_sym(ikq)
  end do

  IF (nspin == 2) THEN
     DO ik = 1, nkstot/2
        DO iq =1, nqs
           index_xkq(nkstot/2+ik,iq) = index_xkq(ik,iq) + nkqs
        END DO
     ENDDO
     do ikq=1,nkqs
        xkq(:,ikq + nkqs)     = xkq(:,ikq)
        index_xk(ikq + nkqs)  = index_xk(ikq) + nkstot/2
        index_sym(ikq + nkqs) = index_sym(ikq)
     end do
     nkqs = 2 * nkqs
  ENDIF
  !
  ! clean up
  !
  deallocate (temp_index_xk, temp_index_sym, temp_index_ikq, new_ikq, temp_xkq)
  !
  ! check that everything is what it should be
  !
  call exx_grid_check () 

  CALL stop_clock ('exx_grid')
  !
  RETURN
  END SUBROUTINE exx_grid_init

  !------------------------------------------------------------------------
  subroutine exx_div_check()
  !------------------------------------------------------------------------
  !
  USE cell_base,  ONLY : bg, at, alat
  USE io_global,  ONLY : stdout
  USE funct,      ONLY : get_screening_parameter
  !
  IMPLICIT NONE
  !
  REAL (DP)     :: atws(3,3)
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
          CALL errore(sub_name,'cannot use x_gamm_extrap and vcut_ws', 1)
     !
  CASE ( "vcut_spherical" ) 
     !
     use_coulomb_vcut_spheric = .TRUE.
     IF ( x_gamma_extrapolation ) &
          CALL errore(sub_name,'cannot use x_gamm_extrap and vcut_spherical', 1)
     !
  CASE ( "none" )
     use_regularization = .FALSE.
     !
  CASE DEFAULT
     CALL errore(sub_name,'invalid exxdiv_treatment: '//TRIM(exxdiv_treatment), 1)
  END SELECT
  !
#ifdef EXXDEBUG
  IF ( ionode ) WRITE (stdout,'(2x,"EXX : q->0 dealt with ",a, " trick" )') &
                TRIM(exxdiv_treatment)
#endif
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
#ifdef EXXDEBUG
  write (stdout,"(2x,'EXX : exx div treatment check successful')")
#endif 
  RETURN
  END SUBROUTINE exx_div_check 


  !------------------------------------------------------------------------
  SUBROUTINE exx_grid_check ( )
  !------------------------------------------------------------------------
  USE symm_base, ONLY : s
  USE cell_base, ONLY : bg, at
  USE lsda_mod,  ONLY : nspin
  USE io_global, ONLY : stdout
  USE klist,     ONLY : nkstot, xk
  implicit none
  real (DP) :: sxk(3), dxk(3), xk_cryst(3), xkk_cryst(3)
  integer :: iq1, iq2, iq3, isym, ik, ikk, ikq, iq
  real (DP) :: eps, dq1, dq2, dq3
  eps = 1.d-6
  dq1= 1.d0/DBLE(nq1)
  dq2= 1.d0/DBLE(nq2)
  dq3= 1.d0/DBLE(nq3)

  do ik =1, nkstot
     IF (pool_para) THEN
        xk_cryst(:) = at(1,:)*xk_collect(1,ik) + at(2,:)*xk_collect(2,ik) + &
                      at(3,:)*xk_collect(3,ik)
     ELSE
        xk_cryst(:) = at(1,:)*xk(1,ik) + at(2,:)*xk(2,ik) + at(3,:)*xk(3,ik)
     ENDIF

     iq = 0
     do iq1=1, nq1
       sxk(1) = xk_cryst(1) + (iq1-1) * dq1
       do iq2 =1, nq2
         sxk(2) = xk_cryst(2) + (iq2-1) * dq2
         do iq3 =1, nq3
            sxk(3) = xk_cryst(3) + (iq3-1) * dq3
            iq = iq + 1
            
            ikq  = index_xkq(ik,iq) 
            ikk  = index_xk(ikq)
            isym = index_sym(ikq)

            IF (pool_para) THEN
               xkk_cryst(:) = at(1,:)*xk_collect(1,ikk)+at(2,:)* &
                                xk_collect(2,ikk)+at(3,:)*xk_collect(3,ikk)
            ELSE
               xkk_cryst(:) = at(1,:)*xk(1,ikk)+at(2,:)*xk(2,ikk) &
                            + at(3,:)*xk(3,ikk)
            ENDIF

            if (isym < 0 ) xkk_cryst(:) = - xkk_cryst(:)
            isym = abs (isym)
            dxk(:) = s(:,1,isym)*xkk_cryst(1) + &
                     s(:,2,isym)*xkk_cryst(2) + &
                     s(:,3,isym)*xkk_cryst(3) - sxk(:)
            dxk(:) = dxk(:) - nint(dxk(:))
            if ( .not. ( abs(dxk(1)).le.eps .and. &
                         abs(dxk(2)).le.eps .and. &
                         abs(dxk(3)).le.eps )   ) then
                 write(*,*) ik,iq
                 write(*,*) ikq,ikk,isym
                 write(*,*) dxk(:)
                 call errore('exx_grid_check', &
                             'something wrong', 1 )
            end if

         end do
       end do
     end do
  end do
#ifdef EXXDEBUG
  write (stdout,"(2x,'EXX : grid check successful')")
#endif
  return

  end subroutine exx_grid_check

  !------------------------------------------------------------------------
  subroutine exx_restart(l_exx_was_active)
  !------------------------------------------------------------------------
    !This subroutine is called when restarting an exx calculation
    use funct,                ONLY : get_exx_fraction, start_exx, exx_is_active, &
                                     get_screening_parameter
    USE fft_base,             ONLY : dffts
    USE io_global,            ONLY : stdout

    implicit none
    logical, intent(in) :: l_exx_was_active
    logical :: exst

    if (.not. l_exx_was_active ) return ! nothing had happpened yet
    !!
    exx_nwordwfc=2*dffts%nnr
    !iunexx = find_free_unit()
    !call diropn(iunexx,'exx', exx_nwordwfc, exst) 
    erfc_scrlen = get_screening_parameter()
    exxdiv = exx_divergence() 
    exxalfa = get_exx_fraction()
#ifdef EXXDEBUG
    write (stdout,*) " ! EXXALFA SET TO ", exxalfa
#endif
    call start_exx
    call weights()
    call exxinit()
    fock0 = exxenergy2()
 
    return
  end subroutine exx_restart
  !------------------------------------------------------------------------
  subroutine exxinit()
  !------------------------------------------------------------------------

    !This subroutine is run before the first H_psi() of each iteration.
    !It saves the wavefunctions for the right density matrix. in real space
    !It saves all the wavefunctions in a single file called prefix.exx
    !
    USE wavefunctions_module, ONLY : evc  
    USE io_files,             ONLY : nwordwfc, iunwfc, iunigk, &
                                     tmp_dir, prefix
    USE io_global,            ONLY : stdout
    USE buffers,              ONLY : get_buffer
    USE gvecs,              ONLY : nls, nlsm, doublegrid
    USE wvfct,                ONLY : nbnd, npwx, npw, igk, wg, et
    USE control_flags,        ONLY : gamma_only
    USE klist,                ONLY : wk, ngk, nks, nkstot
    USE symm_base,            ONLY : nsym, s, ftau

    use mp_global,            ONLY : nproc_pool, me_pool, nproc_bgrp, me_bgrp, &
                                     init_index_over_band, inter_bgrp_comm, &
                                     mpime, inter_pool_comm
    use mp,                   ONLY : mp_sum
    use funct,                ONLY : get_exx_fraction, start_exx, exx_is_active, &
                                     get_screening_parameter 
    USE fft_base,             ONLY : cgather_smooth, cscatter_smooth,&
         & dffts, cgather_custom, cscatter_custom
    use fft_interfaces,       ONLY : invfft

    implicit none
    integer :: ik,ibnd, i, j, k, ir, ri, rj, rk, isym, ikq
    integer :: h_ibnd, half_nbnd
    COMPLEX(DP),allocatable :: temppsic(:), psic(:), tempevc(:,:)
    INTEGER :: nxxs, nrxxs, nr1x,nr2x,nr3x,nr1,nr2,nr3
#ifdef __MPI
    COMPLEX(DP),allocatable :: temppsic_all(:), psic_all(:)
#endif
    INTEGER :: current_ik
    logical, allocatable :: present(:)
    logical :: exst
    INTEGER, ALLOCATABLE :: rir(:,:)
    
    COMPLEX(kind=DP), ALLOCATABLE :: state_fc_t(:,:),evc_g(:)
    integer       :: find_current_k


    call start_clock ('exxinit')

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
    ALLOCATE(psic_all(nxxs), temppsic_all(nxxs) )
#endif
    CALL init_index_over_band(inter_bgrp_comm,nbnd)
    ALLOCATE(temppsic(nrxxs))
    allocate(present(nsym),rir(nxxs,nsym))
    allocate( psic(nrxxs),tempevc( npwx, nbnd ))

    if( .not. allocated( exxbuff ) ) allocate( exxbuff( nrxxs, nkqs, nbnd ) )


    !write(*,*) 'debug: exxbuff size=',size(exxbuff)
    exx_nwordwfc=2*nrxxs
    if (.not.exx_is_active()) then 
       !iunexx = find_free_unit()
       !call diropn(iunexx,'exx', exx_nwordwfc, exst) 
       erfc_scrlen = get_screening_parameter()
       exxdiv = exx_divergence() 
       exxalfa = get_exx_fraction()
#ifdef EXXDEBUG
       write (stdout,*) " ! EXXALFA SET TO ", exxalfa
write(stdout,*) "exxinit, erfc_scrlen set to: ", erfc_scrlen
write(stdout,*) "exxinit, yukawa set to: ", yukawa
#endif
     !
       call start_exx
    endif

#ifdef __MPI
    IF (pool_para) THEN
       IF ( .NOT.ALLOCATED (wg_collect) ) ALLOCATE(wg_collect(nbnd,nkstot))
       CALL wg_all(wg_collect, wg, nkstot, nks)
    ENDIF
#endif

    IF ( nks > 1 ) REWIND( iunigk )

    present(1:nsym) = .false.
    do ikq =1,nkqs
       isym = abs(index_sym(ikq))
       if (.not. present(isym) ) then
          present(isym) = .true.
          if ( mod (s (2, 1, isym) * nr1, nr2) .ne.0 .or. &
               mod (s (3, 1, isym) * nr1, nr3) .ne.0 .or. &
               mod (s (1, 2, isym) * nr2, nr1) .ne.0 .or. &
               mod (s (3, 2, isym) * nr2, nr3) .ne.0 .or. &
               mod (s (1, 3, isym) * nr3, nr1) .ne.0 .or. &
               mod (s (2, 3, isym) * nr3, nr2) .ne.0 ) then
             call errore ('exxinit',' EXX + smooth grid is not working',isym)
          end if
          do ir=1, nxxs
             rir(ir,isym) = ir
          end do
          do k = 1, nr3
             do j = 1, nr2
                do i = 1, nr1
                   call ruotaijk (s(1,1,isym), ftau(1,isym), i, j, k, &
                                  nr1,nr2,nr3, ri, rj , rk )
                   ir =   i + ( j-1)*nr1x + ( k-1)*nr1x*nr2x
                   rir(ir,isym) = ri + (rj-1)*nr1x + (rk-1)*nr1x*nr2x
                end do
             end do
          end do

       end if
    end do

    exxbuff=(0.0_DP,0.0_DP)
    ! set appropriately the x_occupation
    do ik =1,nkstot
       IF (pool_para) THEN
          x_occupation(1:nbnd,ik) = wg_collect (1:nbnd, ik) / wk_collect(ik)
       ELSE
          x_occupation(1:nbnd,ik) = wg(1:nbnd, ik) / wk(ik)
       ENDIF
    end do
!
!   This is parallelized over pool. Each pool computes only its k-points
!
    DO ik = 1, nks
       npw = ngk (ik)
       IF ( nks > 1 ) THEN
          READ( iunigk ) igk
          CALL get_buffer (tempevc, nwordwfc, iunwfc, ik)
       ELSE
          tempevc(1:npwx,1:nbnd) = evc(1:npwx,1:nbnd)
       ENDIF
       IF (pool_para) THEN
          current_ik=find_current_k(ik, nkstot, nks)
       ELSE
          current_ik=ik
       ENDIF

       if (gamma_only) then
          half_nbnd = ( nbnd + 1 )/2
          h_ibnd = 0

          ALLOCATE(state_fc_t(exx_fft_g2r%npwt,nbnd))
          ALLOCATE( evc_g( exx_fft_g2r%ngmt_g ) )
          state_fc_t=( 0.D0, 0.D0 )
          DO ibnd=1,nbnd !exx_cus%nbndv
             CALL exx_grid_convert(tempevc(:,ibnd), npw, exx_fft_g2r,&
                  & state_fc_t(:,ibnd), 1 )  
          ENDDO
          do ibnd =1, nbnd, 2     
             h_ibnd = h_ibnd + 1
             !
             temppsic(:) = ( 0.D0, 0.D0 )
             !
             if (ibnd < nbnd) then
                temppsic(exx_fft_g2r%nlt(1:exx_fft_g2r%npwt)) = state_fc_t(1:exx_fft_g2r%npwt,ibnd)  &
                          + ( 0.D0, 1.D0 ) * state_fc_t(1:exx_fft_g2r%npwt,ibnd+1)
                temppsic(exx_fft_g2r%nltm(1:exx_fft_g2r%npwt)) = CONJG( state_fc_t(1:exx_fft_g2r%npwt,ibnd) ) &
                          + ( 0.D0, 1.D0 ) * CONJG( state_fc_t(1:exx_fft_g2r%npwt,ibnd+1) )
             else
                temppsic(exx_fft_g2r%nlt (1:exx_fft_g2r%npwt)) = state_fc_t(1:exx_fft_g2r%npwt,ibnd) 
                temppsic(exx_fft_g2r%nltm(1:exx_fft_g2r%npwt)) = CONJG( state_fc_t(1:exx_fft_g2r%npwt,ibnd) ) 
             end if

             CALL invfft ('CustomWave', temppsic, exx_fft_g2r%dfftt)

             do ikq=1,nkqs
                if (index_xk(ikq) .ne. current_ik) cycle
                isym = abs(index_sym(ikq) )
#ifdef __MPI
                CALL cgather_custom(temppsic,temppsic_all, exx_fft_g2r%dfftt)
                IF ( me_bgrp == 0 ) &
                psic_all(1:nxxs) = temppsic_all(rir(1:nxxs,isym))
                CALL cscatter_custom(psic_all,psic, exx_fft_g2r%dfftt)
#else
                psic(1:nrxxs) = temppsic(rir(1:nrxxs,isym))
#endif
                if (index_sym(ikq) < 0 ) &
                   call errore('exxinit','index_sym < 0 with gamma_only (!?)',1)

                do ji=1, nrxxs
                   exxbuff(ji,ikq,h_ibnd)=psic(ji)
                enddo
                !CALL davcio(psic,exx_nwordwfc,iunexx,(ikq-1)*half_nbnd+h_ibnd,1)
             end do
          END DO
          DEALLOCATE(state_fc_t, evc_g)
       else
          do ibnd =1, nbnd     
             temppsic(:) = ( 0.D0, 0.D0 )
             temppsic(nls(igk(1:npw))) = tempevc(1:npw,ibnd)
             CALL invfft ('Wave', temppsic, dffts)

             do ikq=1,nkqs
                if (index_xk(ikq) .ne. current_ik) cycle

                isym = abs(index_sym(ikq) )
#ifdef __MPI
                call cgather_smooth(temppsic,temppsic_all)
                IF ( me_bgrp == 0 ) &
                    psic_all(1:nxxs) = temppsic_all(rir(1:nxxs,isym))
                call cscatter_smooth(psic_all,psic)
#else
                psic(1:nrxxs) = temppsic(rir(1:nrxxs,isym))
#endif
                if (index_sym(ikq) < 0 ) psic(1:nrxxs) = CONJG(psic(1:nrxxs))
                do ji=1, nrxxs
                   exxbuff(ji,ikq,ibnd)=psic(ji)
                   
                enddo
                !CALL davcio(psic,exx_nwordwfc,iunexx,(ikq-1)*nbnd+ibnd,1)
             end do
          end do
       end if
    end do

    IF (pool_para) CALL mp_sum(exxbuff, inter_pool_comm)

    deallocate(temppsic, psic,tempevc)
    deallocate(present,rir)
#ifdef __MPI
    deallocate(temppsic_all, psic_all)
#endif 

    call stop_clock ('exxinit')  

  end subroutine exxinit
  
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
    USE constants, ONLY : fpi, e2, pi
    USE cell_base, ONLY : alat, omega, bg, at, tpiba
    USE symm_base, ONLY : nsym, s
    USE gvect,     ONLY : ngm
    USE gvecs,   ONLY : nls, nlsm, doublegrid
    USE wvfct,     ONLY : nbnd, npwx, npw, igk, current_k
    USE control_flags, ONLY : gamma_only
    USE klist,     ONLY : xk, nks, nkstot
    USE lsda_mod,  ONLY : lsda, current_spin, isk
    USE gvect,     ONLY : g, nl
    use fft_base,  ONLY : dffts
    use fft_interfaces, ONLY : fwfft, invfft

    USE parallel_include  
    USE mp_global, ONLY : nproc, inter_image_comm, me_image, my_image_id,&
         & nimage, nproc_image, ibnd_start, ibnd_end, mpime, inter_bgrp_comm, intra_bgrp_comm,&
         & my_bgrp_id, nbgrp
    USE mp,        ONLY : mp_sum, mp_barrier
    USE wvfct,        ONLY : ecutwfc
    USE wavefunctions_module, ONLY : psic

    IMPLICIT NONE

    INTEGER                     :: lda, n, m
    COMPLEX(DP)                 :: psi(lda,m) 
    COMPLEX(DP)                 :: hpsi(lda,m)

    INTEGER          :: nqi, myrank, mysize

    ! local variables
    COMPLEX(DP), allocatable :: tempphic(:), temppsic(:), result(:)
    COMPLEX(DP), allocatable :: rhoc(:), vc(:)
    REAL (DP),   ALLOCATABLE :: fac(:)
    INTEGER          :: ibnd, ik, im , ig, ikq, iq, isym, iqi
    INTEGER          :: h_ibnd, half_nbnd, ierr, nrxxs
    INTEGER          :: current_ik
    REAL(DP) :: x1, x2
    REAL(DP) :: qq, xk_cryst(3), sxk(3), xkq(3), x, q(3)
    ! <LMS> temp array for vcut_spheric
    REAL(DP) :: atws(3,3)
    integer       :: find_current_k

    COMPLEX(kind=DP), ALLOCATABLE :: psi_t(:), prod_tot(:)
    INTEGER, ALLOCATABLE :: igkt(:)

    CALL start_clock ('vexx')

    IF(gamma_only) THEN
       ALLOCATE( fac(exx_fft_r2g%ngmt) )
       nrxxs= exx_fft_g2r%dfftt%nnr
       ALLOCATE( psi_t( exx_fft_g2r%npwt ) )
       ALLOCATE( prod_tot(nrxxs) )
       ALLOCATE( igkt( exx_fft_g2r%ngmt ) )
    ELSE
       ALLOCATE( fac(ngm) )
       nrxxs = dffts%nnr
    ENDIF

    ALLOCATE (tempphic(nrxxs), temppsic(nrxxs), result(nrxxs), &
              rhoc(nrxxs), vc(nrxxs))

    ! write (*,*) exx_nwordwfc,lda,n,m, lda*n
!
! Was used for parallelization on images
!    nqi=nqs/nimage
    nqi=nqs
!
    IF (pool_para) THEN
       current_ik=find_current_k(current_k,nkstot,nks)
    ELSE
       current_ik = current_k
    ENDIF

    if(my_bgrp_id>0) then
      hpsi=(0.0_DP,0.0_DP)
      psi=(0.0_DP,0.0_DP)
    endif
    if (nbgrp>1) then
       call mp_sum(hpsi,inter_bgrp_comm)
       call mp_sum(psi,inter_bgrp_comm)
    endif

    DO im=1,m !for each band of psi (the k cycle is outside band)
       temppsic(:) = ( 0.D0, 0.D0 )

       IF(gamma_only) THEN
          prod_tot(:) = (0.d0,0.d0)
          CALL exx_grid_convert( psi(:,im), npw, exx_fft_g2r, psi_t,&
            & 1, igkt )
          temppsic(exx_fft_g2r%nlt(1:exx_fft_g2r%npwt)) =&
               & psi_t(1:exx_fft_g2r%npwt) 
          temppsic(exx_fft_g2r%nltm(1:exx_fft_g2r%npwt)) =&
               & CONJG(psi_t(1:exx_fft_g2r%npwt))
          CALL invfft ('CustomWave', temppsic, exx_fft_g2r%dfftt)

       ELSE
          temppsic(nls(igk(1:npw))) = psi(1:npw,im)
          CALL invfft ('Wave', temppsic, dffts)
       ENDIF


       result(:)   = (0.d0,0.d0)

       DO iqi=1,nqi
!
! Was used for parallelization on images
!          iq=iqi+nqi*my_image_id
          iq=iqi
!
          ikq  = index_xkq(current_ik,iq)
          ik   = index_xk(ikq)
          isym = ABS(index_sym(ikq))

          IF (pool_para) THEN
             xk_cryst(:) = at(1,:)*xk_collect(1,ik) + at(2,:)*xk_collect(2,ik)&
                         + at(3,:)*xk_collect(3,ik)
          ELSE
             xk_cryst(:) = at(1,:)*xk(1,ik) + at(2,:)*xk(2,ik) &
                         + at(3,:)*xk(3,ik)
          ENDIF
          IF (index_sym(ikq) < 0 ) xk_cryst = - xk_cryst
          sxk(:) = s(:,1,isym)*xk_cryst(1) + &
               s(:,2,isym)*xk_cryst(2) + &
               s(:,3,isym)*xk_cryst(3) 
          xkq(:) = bg(:,1)*sxk(1) + bg(:,2)*sxk(2) + bg(:,3)*sxk(3)
          !
          ! calculate the 1/|r-r'| factor and place it in fac
          !
          IF(gamma_only) THEN
             CALL g2_convolution(exx_fft_r2g%ngmt, exx_fft_r2g%gt, xk(:&
                  &,current_k), xkq, fac) 
          ELSE
             CALL g2_convolution(ngm, g, xk(:,current_k), xkq, fac)
          ENDIF
          !
          IF (gamma_only) THEN
             half_nbnd = ( nbnd + 1 ) / 2
             h_ibnd = ibnd_start/2
             DO ibnd=ibnd_start,ibnd_end, 2 !for each band of psi
                h_ibnd = h_ibnd + 1
                x1 = x_occupation(ibnd,  ik)
                IF (ibnd < ibnd_end) THEN
                   x2 = x_occupation(ibnd + 1,ik)
                ELSE
                   x2 = 0.d0
                END IF
                IF ( ABS(x1) < 1.d-6 .AND.  ABS(x2) < 1.d-6 ) CYCLE
                !
                !loads the phi from file
                !
#ifdef EXXDEBUG
WRITE(stdout,*) "from vexx, nrxxs: ", nrxxs, ikq, h_ibnd 
call flush_unit(stdout)
#endif
              DO ji=1, nrxxs
                 tempphic(ji)=exxbuff(ji,ikq,h_ibnd)
                   
              ENDDO
              !CALL davcio ( tempphic, exx_nwordwfc, iunexx, &
              !                    (ikq-1)*half_nbnd+h_ibnd, -1 )
              !calculate rho in real space
              rhoc(:)=(0.d0,0.d0)
              rhoc(1:nrxxs)=CONJG(tempphic(1:nrxxs))*temppsic(1:nrxxs) / omega
              !brings it to G-space
              IF (ecutfock == ecutwfc) THEN
                 CALL fwfft ('Custom', rhoc, exx_fft_r2g%dfftt)
                 vc(:) = ( 0.D0, 0.D0 )
                 vc(exx_fft_r2g%nlt(igkt(1:exx_fft_r2g%ngmt)))  =&
                      & fac(1:exx_fft_r2g%ngmt) * rhoc(exx_fft_r2g&
                      &%nlt(igkt(1:exx_fft_r2g%ngmt))) 
                 vc(exx_fft_r2g%nltm(igkt(1:exx_fft_r2g%ngmt))) =&
                      & fac(1:exx_fft_r2g%ngmt) * rhoc(exx_fft_r2g&
                      &%nltm(igkt(1:exx_fft_r2g%ngmt))) 
                 !brings back v in real space
                 CALL invfft ('Custom', vc, exx_fft_r2g%dfftt) 
              ELSE
                 CALL fwfft ('CustomWave', rhoc, exx_fft_r2g%dfftt)
                 vc(:) = ( 0.D0, 0.D0 )
                 vc(exx_fft_r2g%nlt(igkt(1:exx_fft_r2g%npwt)))  =&
                      & fac(1:exx_fft_r2g%npwt) * rhoc(exx_fft_r2g&
                      &%nlt(igkt(1:exx_fft_r2g%npwt))) 
                 vc(exx_fft_r2g%nltm(igkt(1:exx_fft_r2g%npwt))) =&
                      & fac(1:exx_fft_r2g%npwt) * rhoc(exx_fft_r2g&
                      &%nltm(igkt(1:exx_fft_r2g%npwt))) 
                 !brings back v in real space
                 CALL invfft ('CustomWave', vc, exx_fft_r2g%dfftt) 
              ENDIF
   

              vc = CMPLX( x1 * DBLE (vc), x2 * AIMAG(vc) ,kind=DP)/ nqs

              !accumulates over bands and k points
!              result(1:nrxxs) = result(1:nrxxs) + &
!                   DBLE( vc(1:nrxxs) * tempphic(1:nrxxs) )
              prod_tot(:) = prod_tot(:) + DBLE( vc(:) * tempphic(:))
           END DO

        ELSE
           DO ibnd=ibnd_start,ibnd_end !for each band of psi
              IF ( ABS(x_occupation(ibnd,ik)) < 1.d-6) CYCLE
              !
              !loads the phi from file
              !
#ifdef EXXDEBUG
WRITE(stdout,*) "from vexx, nrxxs: ", nrxxs, ikq, ibnd
if(allocated(tempphic)) write(stdout,*) "tempphic allcated"
if(allocated(exxbuff)) write(stdout,*) "exxbuff allocated"
call flush_unit(stdout)
#endif
              DO ji=1, nrxxs
                 tempphic(ji)=exxbuff(ji,ikq,ibnd)                
              ENDDO
                
              !CALL davcio ( tempphic, exx_nwordwfc, iunexx, &
              !                        (ikq-1)*nbnd+ibnd, -1 )
              !calculate rho in real space
              rhoc(:)=CONJG(tempphic(:))*temppsic(:) / omega
              !brings it to G-space
              CALL fwfft ('Smooth', rhoc, dffts)
   
              vc(:) = ( 0.D0, 0.D0 )
              vc(nls(1:ngm)) = fac(1:ngm) * rhoc(nls(1:ngm))

              vc = vc * x_occupation(ibnd,ik) / nqs

              !brings back v in real space
              CALL invfft ('Smooth', vc, dffts) 
                
              !accumulates over bands and k points
              result(1:nrxxs)=result(1:nrxxs)+vc(1:nrxxs)*tempphic(1:nrxxs)
           END DO
        END IF
     END DO


     IF(gamma_only) THEN
        CALL mp_sum( prod_tot(1:nrxxs), inter_bgrp_comm)
     ELSE
        CALL mp_sum( result(1:nrxxs), inter_bgrp_comm)
     ENDIF

     !
     ! Was used for parallelization on images
     !       CALL mp_sum( result(1:nrxxs), inter_image_comm )
     !brings back result in G-space
        IF( gamma_only) THEN
           CALL fwfft( 'CustomWave' , prod_tot, exx_fft_g2r%dfftt )
           psic(1:exx_fft_g2r%npwt) = prod_tot(exx_fft_g2r&
                &%nlt(igkt(1:exx_fft_g2r%npwt)))

           result(:)=(0.d0,0.d0)
           CALL exx_grid_convert(psic,npw, exx_fft_g2r, result, -1)
           hpsi(1:npw,im)=hpsi(1:npw,im) - exxalfa*result(1:npw)
        ELSE
           CALL fwfft ('Wave', result, dffts)
           !adds it to hpsi
           hpsi(1:npw,im)=hpsi(1:npw,im) - exxalfa*result(nls(igk(1:npw)))
        ENDIF
  END DO
    
  DEALLOCATE (tempphic,temppsic, result, rhoc, vc, fac )

  IF(gamma_only) DEALLOCATE( psi_t, prod_tot, igkt )

  CALL stop_clock ('vexx')

END SUBROUTINE vexx
!-----------------------------------------------------------------------
SUBROUTINE g2_convolution(ngm, g, xk, xkq, fac)
!-----------------------------------------------------------------------
  ! This routine calculates the 1/|r-r'| part of the exact exchange 
  ! expression in reciprocal space (the G^-2 factor).
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
     qq = SUM(q(:)**2.0_DP) 
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
     ELSE IF (qq.GT.1.d-8) THEN
        !
        IF ( erfc_scrlen > 0  ) THEN
           fac(ig)=e2*fpi/qq*(1.d0-EXP(-qq/4.d0/erfc_scrlen**2)) * grid_factor
        ELSEIF( erf_scrlen > 0 ) THEN
           fac(ig)=e2*fpi/qq*(EXP(-qq/4.d0/erf_scrlen**2)) * grid_factor
        ELSE
           fac(ig)=e2*fpi/( qq + yukawa ) * grid_factor
        END IF
        IF (on_double_grid) fac(ig) = 0.d0
        !
     ELSE
        !
        fac(ig)= - exxdiv ! or rather something else (see F.Gygi)
        !
        IF ( yukawa > 0.d0.AND. .NOT. x_gamma_extrapolation ) &
             fac(ig) = fac(ig) + e2*fpi/( qq + yukawa )
        IF( erfc_scrlen > 0.d0.AND. .NOT. x_gamma_extrapolation ) fac(ig) = fac(ig) + e2*pi/(erfc_scrlen**2)
        !
     ENDIF
     !
  ENDDO
  !CALL stop_clock ('vexx_ngmloop')
END SUBROUTINE g2_convolution
!-----------------------------------------------------------------------

  function exxenergy ()
    ! This function is called to correct the deband value and have 
    ! the correct energy 
    USE io_files,   ONLY : iunigk,iunwfc, nwordwfc
    USE buffers,    ONLY : get_buffer
    USE wvfct,      ONLY : nbnd, npwx, npw, igk, wg, current_k
    USE control_flags, ONLY : gamma_only
    USE gvect,      ONLY : gstart
    USE wavefunctions_module, ONLY : evc
    USE lsda_mod,   ONLY : lsda, current_spin, isk
    USE klist,      ONLY : ngk, nks
    USE mp_global,  ONLY : inter_pool_comm, inter_bgrp_comm, intra_bgrp_comm, nbgrp
    USE mp,         ONLY : mp_sum

    implicit none
    REAL (DP)   :: exxenergy,  energy
    INTEGER          :: ibnd, ik
    COMPLEX(DP) :: vxpsi ( npwx, nbnd ), psi(npwx,nbnd)
    COMPLEX(DP) :: zdotc

    call start_clock ('exxenergy')

    energy=0.d0
    IF ( nks > 1 ) REWIND( iunigk )
    do ik=1,nks
       current_k = ik
       IF ( lsda ) current_spin = isk(ik)
       npw = ngk (ik)
       IF ( nks > 1 ) THEN
          READ( iunigk ) igk
          call get_buffer  (psi, nwordwfc, iunwfc, ik)
       ELSE
          psi(1:npwx,1:nbnd) = evc(1:npwx,1:nbnd)
       END IF
       vxpsi(:,:) = (0.d0, 0.d0)
       call vexx(npwx,npw,nbnd,psi,vxpsi)
       do ibnd=1,nbnd
          energy = energy + &
                   wg(ibnd,ik) * zdotc(npw,psi(1,ibnd),1,vxpsi(1,ibnd),1)
       end do
       if (gamma_only .and. gstart == 2) then
           do ibnd=1,nbnd
              energy = energy - &
                       0.5d0 * wg(ibnd,ik) * CONJG(psi(1,ibnd)) * vxpsi(1,ibnd)
           end do
       end if
    end do

    if (gamma_only) energy = 2.d0 * energy

    call mp_sum( energy, intra_bgrp_comm)
    call mp_sum( energy, inter_pool_comm )
    ! 
    exxenergy = energy
    !
    call stop_clock ('exxenergy')
  end function exxenergy

  !-----------------------------------------------------------------------
  function exxenergy2()
  !-----------------------------------------------------------------------
    !
    USE constants, ONLY : fpi, e2, pi
    USE io_files,  ONLY : iunigk,iunwfc, nwordwfc
    USE buffers,   ONLY : get_buffer
    USE cell_base, ONLY : alat, omega, bg, at, tpiba
    USE symm_base,ONLY : nsym, s
    USE gvect,     ONLY : ngm, gstart
    USE gvecs,   ONLY : nls, nlsm, doublegrid
    USE wvfct,     ONLY : nbnd, npwx, npw, igk, wg, current_k
    USE control_flags, ONLY : gamma_only
    USE wavefunctions_module, ONLY : evc
    USE klist,     ONLY : xk, ngk, nks, nkstot
    USE lsda_mod,  ONLY : lsda, current_spin, isk
    USE gvect,     ONLY : g, nl
    USE mp_global, ONLY : inter_pool_comm, inter_image_comm, inter_bgrp_comm, intra_bgrp_comm, nbgrp
    USE mp_global, ONLY : my_image_id, nimage, ibnd_start, ibnd_end
    USE mp,        ONLY : mp_sum
    use fft_base,  ONLY : dffts
    use fft_interfaces, ONLY : fwfft, invfft
    USE wvfct,     ONLY : ecutwfc

    IMPLICIT NONE
    REAL (DP)   :: exxenergy2,  energy

    ! local variables
    COMPLEX(DP), allocatable :: tempphic(:), temppsic(:)
    COMPLEX(DP), ALLOCATABLE :: rhoc(:)
    REAL (DP),   ALLOCATABLE :: fac(:)
    integer          :: jbnd, ibnd, ik, ikk, ig, ikq, iq, isym
    integer          :: half_nbnd, h_ibnd, nqi, iqi, nrxxs
    real(DP)    :: x1, x2
    REAL(DP) :: qq, xk_cryst(3), sxk(3), xkq(3), vc, x, q(3)
    ! temp array for vcut_spheric
    real(DP) :: atws(3,3) 
    integer       :: find_current_k

    COMPLEX(kind=DP), ALLOCATABLE :: psi_t(:), prod_tot(:)
    INTEGER, ALLOCATABLE :: igkt(:)
    INTEGER :: ngm_fft

    call start_clock ('exxen2')

    IF(gamma_only) THEN
       nrxxs= exx_fft_g2r%dfftt%nnr
!       nrxxs = exx_fft_g2r%nrxxt
       ALLOCATE( fac(exx_fft_r2g%ngmt) )
       ALLOCATE( psi_t( exx_fft_g2r%npwt ) )
       ALLOCATE( prod_tot(exx_fft_g2r%nrxxt) )
       ALLOCATE( igkt( exx_fft_g2r%ngmt ) )
    ELSE
       nrxxs = dffts%nnr
       ALLOCATE( fac(ngm) )
    ENDIF

    ALLOCATE (tempphic(nrxxs), temppsic(nrxxs), rhoc(nrxxs) )

    energy=0.d0

!
! Was used for parallelization on images    
!    nqi=nqs/nimage
    nqi=nqs
!
    IF ( nks > 1 ) REWIND( iunigk )
    do ikk=1,nks
       IF (pool_para) THEN
          current_k=find_current_k(ikk,nkstot,nks)
       ELSE
          current_k = ikk
       ENDIF
       IF ( lsda ) current_spin = isk(ikk)
       npw = ngk (ikk)
       IF ( nks > 1 ) THEN
          READ( iunigk ) igk
          call get_buffer (evc, nwordwfc, iunwfc, ikk)
       END IF
       do jbnd=ibnd_start, ibnd_end !for each band of psi (the k cycle is outside band)
          temppsic(:) = ( 0.D0, 0.D0 )
       IF(gamma_only) THEN
          CALL exx_grid_convert( evc(:,jbnd), npw, exx_fft_g2r, psi_t,&
            & 1, igkt )
          temppsic(exx_fft_g2r%nlt(1:exx_fft_g2r%npwt)) =&
               & psi_t(1:exx_fft_g2r%npwt) 
          temppsic(exx_fft_g2r%nltm(1:exx_fft_g2r%npwt)) =&
               & CONJG(psi_t(1:exx_fft_g2r%npwt))
          CALL invfft ('CustomWave', temppsic, exx_fft_g2r%dfftt)
       ELSE
          temppsic(nls(igk(1:npw))) = evc(1:npw,jbnd)
          CALL invfft ('Wave', temppsic, dffts)
       ENDIF
       
          do iqi=1,nqi
!
! Was used for parallelization on images
!             iq=iqi+nqi*my_image_id
             iq=iqi
!
             ikq  = index_xkq(current_k,iq)
             ik   = index_xk(ikq)
             isym = abs(index_sym(ikq))

             IF (pool_para) THEN
                xk_cryst(:) = at(1,:)*xk_collect(1,ik) + &
                              at(2,:)*xk_collect(2,ik) + &
                              at(3,:)*xk_collect(3,ik)
             ELSE
                xk_cryst(:) = at(1,:)*xk(1,ik) + at(2,:)*xk(2,ik) + &
                              at(3,:)*xk(3,ik)
             ENDIF

             if (index_sym(ikq) < 0 ) xk_cryst = - xk_cryst
             sxk(:) = s(:,1,isym)*xk_cryst(1) + &
                      s(:,2,isym)*xk_cryst(2) + &
                      s(:,3,isym)*xk_cryst(3) 
             xkq(:) = bg(:,1)*sxk(1) + bg(:,2)*sxk(2) + bg(:,3)*sxk(3)

          IF(gamma_only) THEN
             CALL g2_convolution(exx_fft_r2g%ngmt, exx_fft_r2g%gt, xk(:&
                  &,current_k), xkq, fac) 
             fac(exx_fft_r2g%gstart_t:) = 2.d0 * fac(exx_fft_r2g%gstart_t:)
          ELSE
             CALL g2_convolution(ngm, g, xk(:,current_k), xkq, fac)
          ENDIF

             if (gamma_only) then
                half_nbnd = ( nbnd + 1) / 2
                h_ibnd = 0
                do ibnd=1,nbnd, 2 !for each band of psi
                   h_ibnd = h_ibnd + 1
                   x1 = x_occupation(ibnd,ik)
                   if ( ibnd < nbnd ) then
                      x2 = x_occupation(ibnd+1,ik)
                   else
                      x2 = 0.d0
                   end if
                   if ( abs(x1) < 1.d-6 .and. abs(x2) < 1.d-6 ) cycle
                   !
                   !loads the phi from file
                   !
                   do ji=1, nrxxs
                      tempphic(ji)=exxbuff(ji,ikq,h_ibnd)
                   enddo
                
                   !CALL davcio (tempphic, exx_nwordwfc, iunexx, &
                   !                       (ikq-1)*half_nbnd+h_ibnd, -1 )
                   !calculate rho in real space
                   rhoc(:)=(0.d0, 0.d0)
                   rhoc(1:nrxxs)=CONJG(tempphic(1:nrxxs))*temppsic(1:nrxxs) / omega
                   IF (ecutfock == ecutwfc) THEN
                      !brings it to G-space
                      CALL fwfft ('Custom', rhoc, exx_fft_r2g%dfftt)
                      vc = 0.D0
                      DO ig=1,exx_fft_r2g%ngmt
                         vc = vc + fac(ig) * x1 * &
                              ABS( rhoc(exx_fft_r2g%nlt(igkt(ig))) + CONJG(rhoc(exx_fft_r2g%nltm(igkt(ig)))) )**2
                         vc = vc + fac(ig) * x2 * &
                              ABS( rhoc(exx_fft_r2g%nlt(igkt(ig))) - CONJG(rhoc(exx_fft_r2g%nltm(igkt(ig)))) )**2
                      END DO
                   ELSE
                      !brings it to G-space
                      CALL fwfft ('CustomWave', rhoc, exx_fft_r2g%dfftt)
                      vc = 0.D0
                      DO ig=1,exx_fft_r2g%npwt
                         vc = vc + fac(ig) * x1 * &
                              ABS( rhoc(exx_fft_r2g%nlt(igkt(ig))) + CONJG(rhoc(exx_fft_r2g%nltm(igkt(ig)))) )**2
                         vc = vc + fac(ig) * x2 * &
                              ABS( rhoc(exx_fft_r2g%nlt(igkt(ig))) - CONJG(rhoc(exx_fft_r2g%nltm(igkt(ig)))) )**2
                      END DO
                   ENDIF
                   vc = vc * omega * 0.25d0 / nqs

                   energy = energy - exxalfa * vc * wg(jbnd,ikk)
                end do
             else
                do ibnd=1,nbnd !for each band of psi
                   if ( abs(x_occupation(ibnd,ik)) < 1.d-6) cycle
                   !
                   !loads the phi from file
                   !
                   do ji=1, nrxxs
                      tempphic(ji)=exxbuff(ji,ikq,ibnd)

                   enddo

                   !CALL davcio (tempphic, exx_nwordwfc, iunexx, &
                   !                       (ikq-1)*nbnd+ibnd, -1 )
                   !calculate rho in real space
                   rhoc(:)=CONJG(tempphic(:))*temppsic(:) / omega
                   !brings it to G-space
                   CALL fwfft ('Smooth', rhoc, dffts)
   
                   vc = 0.D0
                   do ig=1,ngm
                      vc = vc + fac(ig) * rhoc(nls(ig)) * CONJG(rhoc(nls(ig)))
                   end do
                   vc = vc * omega * x_occupation(ibnd,ik) / nqs

                   energy = energy - exxalfa * vc * wg(jbnd,ikk)
                end do
             end if
          end do
       end do
    end do

    DEALLOCATE (tempphic, temppsic, rhoc, fac )
    IF(gamma_only) DEALLOCATE( psi_t, prod_tot, igkt )
!
! Was used for image parallelization
!    call mp_sum( energy, inter_image_comm )
!
    call mp_sum( energy, inter_bgrp_comm )
    call mp_sum( energy, intra_bgrp_comm )
    call mp_sum( energy, inter_pool_comm )
    !
    exxenergy2 = energy
    !
    call stop_clock ('exxen2')

  end function  exxenergy2

  function exx_divergence ()

     USE constants, ONLY : fpi, e2, pi
     USE cell_base, ONLY : bg, at, alat, omega
     USE gvect,     ONLY : ngm, g
     USE wvfct,     ONLY : ecutwfc
     USE io_global, ONLY : stdout
     USE control_flags, ONLY : gamma_only
     USE mp_global, ONLY : intra_bgrp_comm
     USE mp,        ONLY : mp_sum

     implicit none
     real(DP) :: exx_divergence

     ! local variables
     integer :: iq1,iq2,iq3, ig
     real(DP) :: div, dq1, dq2, dq3, xq(3), q_, qq, tpiba2, alpha, x, q(3)

     integer :: nqq, iq
     real(DP) :: aa, dq

     call start_clock ('exx_div')

     tpiba2 = (fpi / 2.d0 / alat) **2

     alpha  = 10.d0 * tpiba2 / ecutwfc

     IF ( .NOT. use_regularization ) THEN
        exx_divergence = 0.d0
        return
     END IF

     dq1= 1.d0/DBLE(nq1)
     dq2= 1.d0/DBLE(nq2) 
     dq3= 1.d0/DBLE(nq3) 

     div = 0.d0
     do iq1=1,nq1
        do iq2=1,nq2
           do iq3=1,nq3
              xq(:) = bg(:,1) * (iq1-1) * dq1 + &
                      bg(:,2) * (iq2-1) * dq2 + &
                      bg(:,3) * (iq3-1) * dq3 
              do ig=1,ngm
                 q(1)= xq(1) + g(1,ig)
                 q(2)= xq(2) + g(2,ig)
                 q(3)= xq(3) + g(3,ig)
                 qq = ( q(1)*q(1) + q(2)*q(2) + q(3)*q(3) ) 
                 if (x_gamma_extrapolation) then
                    on_double_grid = .true.
                    x= 0.5d0*(q(1)*at(1,1)+q(2)*at(2,1)+q(3)*at(3,1))*nq1
                    on_double_grid = on_double_grid .and. (abs(x-nint(x))<eps)
                    x= 0.5d0*(q(1)*at(1,2)+q(2)*at(2,2)+q(3)*at(3,2))*nq2
                    on_double_grid = on_double_grid .and. (abs(x-nint(x))<eps)
                    x= 0.5d0*(q(1)*at(1,3)+q(2)*at(2,3)+q(3)*at(3,3))*nq3
                    on_double_grid = on_double_grid .and. (abs(x-nint(x))<eps)
                 end if
                 if (.not.on_double_grid) then
                    if ( qq > 1.d-8 ) then
                       if ( erfc_scrlen > 0 ) then
                          div = div + exp( -alpha * qq) / qq * &
                                (1.d0-exp(-qq*tpiba2/4.d0/erfc_scrlen**2)) * grid_factor
                       elseif ( erf_scrlen >0 ) then
                          div = div + exp( -alpha * qq) / qq * &
                                (exp(-qq*tpiba2/4.d0/erf_scrlen**2)) * grid_factor
                       else

                          div = div + exp( -alpha * qq) / (qq + yukawa/tpiba2) &
                                                     * grid_factor
                       endif
                    end if
                 end if
              end do
           end do
        end do
     end do
     call mp_sum(  div, intra_bgrp_comm )
     if (gamma_only) then
        div = 2.d0 * div
     endif
     if ( .not. x_gamma_extrapolation ) then
        if ( yukawa > 0.d0) then
           div = div + tpiba2/yukawa
        elseif( erfc_scrlen > 0.d0 ) then
           div = div + tpiba2/4.d0/erfc_scrlen**2
        else
           div = div - alpha
        end if
     end if

     div = div * e2 * fpi / tpiba2 / nqs

     alpha = alpha / tpiba2

     nqq = 100000
     dq = 5.0d0 / sqrt(alpha) /nqq
     aa = 0.d0
     do iq=0,  nqq
        q_ = dq * (iq+0.5d0)
        qq = q_ * q_
        if ( erfc_scrlen > 0 ) then
           aa = aa  -exp( -alpha * qq) * exp(-qq/4.d0/erfc_scrlen**2) * dq
        elseif ( erf_scrlen > 0 ) then
           aa = 0.d0
        else
           aa = aa - exp( -alpha * qq) * yukawa / (qq + yukawa) * dq
        end if
     end do
     aa = aa * 8.d0 /fpi
     aa = aa + 1.d0/sqrt(alpha*0.25d0*fpi) 
     if( erf_scrlen > 0) aa = 1.d0/sqrt((alpha+1.d0/4.d0/erf_scrlen**2)*0.25d0*fpi)
#ifdef EXXDEBUG
     write (stdout,*) aa, 1.d0/sqrt(alpha*0.25d0*fpi)
#endif   
     div = div - e2*omega * aa

!     div = div - e2*omega/sqrt(alpha*0.25d0*fpi)

     exx_divergence = div * nqs
#ifdef EXXDEBUG
     write (stdout,'(a,i4,a,3f12.4)') 'EXX divergence (',nq1,')= ', &
                                  div, alpha
#endif
     call stop_clock ('exx_div')

     return
  end function exx_divergence 
  



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
  use fft_base,  ONLY : dffts
  use fft_interfaces, ONLY : fwfft, invfft
  ! ---- local variables -------------------------------------------------
  IMPLICIT NONE
  real (dp)   :: exx_stress(3,3), exx_stress_(3,3)
  complex(dp), allocatable :: tempphic(:), temppsic(:)
  complex(dp), allocatable :: rhoc(:)
  real(dp), allocatable :: fac(:), fac_tens(:,:,:), fac_stress(:)
  integer :: jbnd, ibnd, ik, ikk, ig, ikq, iq, isym
  integer :: half_nbnd, h_ibnd, nqi, iqi, beta, nrxxs
  real(dp) :: x1, x2
  real(dp) :: qq, xk_cryst(3), sxk(3), xkq(3), vc(3,3), x, q(3)
  ! temp array for vcut_spheric
  real(dp) :: atws(3,3) 
  real(dp) :: delta(3,3)
 
  call start_clock ('exx_stress')

  IF (pool_para) call errore('exx_stress','stress not available with pools',1)
  nrxxs = dffts%nnr
  delta = reshape( (/1.d0,0.d0,0.d0, 0.d0,1.d0,0.d0, 0.d0,0.d0,1.d0/), (/3,3/))
  exx_stress_ = 0.d0
  allocate( tempphic(nrxxs), temppsic(nrxxs), rhoc(nrxxs), fac(ngm) )
  allocate( fac_tens(3,3,ngm), fac_stress(ngm) )

  if ( nks > 1 ) rewind( iunigk )
!
! Was used for image parallelization
!
  nqi=nqs
!  nqi = nqs/nimage
!
  ! loop over k-points
  do ikk = 1, nks
      current_k = ikk
      if (lsda) current_spin = isk(ikk)
      npw = ngk(ikk)

      if (nks > 1) then
          read(iunigk) igk
          call get_buffer(evc, nwordwfc, iunwfc, ikk)
      end if

      ! loop over bands
      do jbnd = 1, nbnd
          temppsic(:) = ( 0.d0, 0.d0 )
          temppsic(nls(igk(1:npw))) = evc(1:npw,jbnd)
          if(gamma_only) temppsic(nlsm(igk(1:npw))) = conjg(evc(1:npw,jbnd))
          CALL invfft ('Wave', temppsic, dffts)       

          do iqi = 1, nqi
!
! Was used for image parallelization
!
!              iq = iqi + nqi*my_image_id
              iq=iqi
!
              ikq  = index_xkq(current_k,iq)
              ik   = index_xk(ikq)
              isym = abs(index_sym(ikq))

              xk_cryst(:)=at(1,:)*xk(1,ik)+at(2,:)*xk(2,ik)+at(3,:)*xk(3,ik)
              if (index_sym(ikq) < 0) xk_cryst = -xk_cryst
              sxk(:) = s(:,1,isym)*xk_cryst(1) + &
                       s(:,2,isym)*xk_cryst(2) + &
                       s(:,3,isym)*xk_cryst(3) 
              xkq(:) = bg(:,1)*sxk(1) + bg(:,2)*sxk(2) + bg(:,3)*sxk(3)

              !CALL start_clock ('exxen2_ngmloop')
              do ig = 1, ngm
                 q(1)= xk(1,current_k) - xkq(1) + g(1,ig)
                 q(2)= xk(2,current_k) - xkq(2) + g(2,ig)
                 q(3)= xk(3,current_k) - xkq(3) + g(3,ig)

                 q = q * tpiba
                 qq = ( q(1)*q(1) + q(2)*q(2) + q(3)*q(3) )

                 do beta = 1, 3
                     fac_tens(1:3,beta,ig) = q(1:3)*q(beta)
                 enddo

                 if (x_gamma_extrapolation) then
                    on_double_grid = .true.
                    x= 0.5d0/tpiba*(q(1)*at(1,1)+q(2)*at(2,1)+q(3)*at(3,1))*nq1
                    on_double_grid = on_double_grid .and. (abs(x-nint(x))<eps)
                    x= 0.5d0/tpiba*(q(1)*at(1,2)+q(2)*at(2,2)+q(3)*at(3,2))*nq2
                    on_double_grid = on_double_grid .and. (abs(x-nint(x))<eps)
                    x= 0.5d0/tpiba*(q(1)*at(1,3)+q(2)*at(2,3)+q(3)*at(3,3))*nq3
                    on_double_grid = on_double_grid .and. (abs(x-nint(x))<eps)
                 endif

                 if (use_coulomb_vcut_ws) then
                    fac(ig) = vcut_get(vcut, q)
                    fac_stress(ig) = 0.d0   ! not implemented
                    if (gamma_only .and. qq > 1.d-8) fac(ig) = 2.d0 * fac(ig)

                 else if ( use_coulomb_vcut_spheric ) then
                    fac(ig) = vcut_spheric_get(vcut, q)
                    fac_stress(ig) = 0.d0   ! not implemented
                    if (gamma_only .and. qq > 1.d-8) fac(ig) = 2.d0 * fac(ig) 

                 else if (qq > 1.d-8) then
                    if ( erfc_scrlen > 0 ) then
                       fac(ig)=e2*fpi/qq*(1.d0-exp(-qq/4.d0/erfc_scrlen**2)) * grid_factor
                       fac_stress(ig) = -e2*fpi * 2.d0/qq**2 * ( &
                           (1.d0+qq/4.d0/erfc_scrlen**2)*exp(-qq/4.d0/erfc_scrlen**2) - 1.d0) * &
                           grid_factor
                    else
                       fac(ig)=e2*fpi/( qq + yukawa ) * grid_factor
                       fac_stress(ig) = 2.d0 * e2*fpi/(qq+yukawa)**2 * grid_factor
                    end if

                    if (gamma_only) fac(ig) = 2.d0 * fac(ig)
                    if (gamma_only) fac_stress(ig) = 2.d0 * fac_stress(ig)
                    if (on_double_grid) fac(ig) = 0.d0
                    if (on_double_grid) fac_stress(ig) = 0.d0

                 else
                    fac(ig)= -exxdiv ! or rather something else (see f.gygi)
                    fac_stress(ig) = 0.d0  ! or -exxdiv_stress (not yet implemented)
                    if ( yukawa> 0.d0 .and. .not. x_gamma_extrapolation) then
                       fac(ig) = fac(ig) + e2*fpi/( qq + yukawa )
                       fac_stress(ig) = 2.d0 * e2*fpi/(qq+yukawa)**2
                    endif
                    if (erfc_scrlen > 0.d0 .and. .not. x_gamma_extrapolation) then
                       fac(ig) = e2*fpi / (4.d0*erfc_scrlen**2)
                       fac_stress(ig) = e2*fpi / (8.d0*erfc_scrlen**4)
                    endif
                 endif
              enddo
              !CALL stop_clock ('exxen2_ngmloop')

              if (gamma_only) then
                  half_nbnd = (nbnd + 1) / 2
                  h_ibnd = 0
                  do ibnd=1,nbnd, 2 !for each band of psi
                      h_ibnd = h_ibnd + 1
                      x1 = x_occupation(ibnd,ik)
                      if ( ibnd < nbnd ) then
                         x2 = x_occupation(ibnd+1,ik)
                      else
                         x2 = 0.d0
                      end if
                      if ( abs(x1) < 1.d-6 .and. abs(x2) < 1.d-6 ) cycle

                      ! loads the phi from file
                      do ji=1, nrxxs
                        tempphic(ji)=exxbuff(ji,ikq,h_ibnd)
                      enddo
                
                      ! calculate rho in real space
                      rhoc(:)=CONJG(tempphic(:))*temppsic(:) / omega
                      ! brings it to G-space
                      CALL fwfft ('Smooth', rhoc, dffts)
   
                      vc = 0.d0
                      do ig=1,ngm
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
              else
                  do ibnd=1,nbnd !for each band of psi
                     if ( abs(x_occupation(ibnd,ik)) < 1.d-6) cycle

                     ! loads the phi from file
                     do ji=1, nrxxs
                       tempphic(ji)=exxbuff(ji,ikq,ibnd)
                     enddo

                     ! calculate rho in real space
                     rhoc(:)=CONJG(tempphic(:))*temppsic(:) / omega

                     ! brings it to G-space
                     CALL fwfft ('Smooth', rhoc, dffts)

                     vc = 0.d0
                     do ig=1,ngm
                       vc(:,:) = vc(:,:) + rhoc(nls(ig))*CONJG(rhoc(nls(ig))) * &
                                 (fac_tens(:,:,ig)*fac_stress(ig)/2.d0 - delta(:,:)*fac(ig))
                     end do
                     vc = vc * x_occupation(ibnd,ik) / nqs / 4.d0
                     exx_stress_ = exx_stress_ + exxalfa * vc * wg(jbnd,ikk)
                  end do
              endif ! gamma or k-points

          enddo ! iqi
      enddo ! jbnd
  enddo ! ikk

  deallocate (tempphic, temppsic, rhoc, fac, fac_tens, fac_stress )
!
! Was used for image parallelization
!  call mp_sum( exx_stress_, inter_image_comm )
!
  call mp_sum( exx_stress_, intra_bgrp_comm )
  call mp_sum( exx_stress_, inter_pool_comm )
  exx_stress = exx_stress_

  call stop_clock ('exx_stress')

  END FUNCTION exx_stress



end module exx
