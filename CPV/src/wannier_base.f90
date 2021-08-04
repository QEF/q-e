!
! Copyright (C) 2002-2005 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
MODULE wannier_base
  !----------------------------------------------------------------------------
  !
  USE kinds, ONLY : DP
  !
  IMPLICIT NONE
  !
  ! ... input variables
  !
  LOGICAL              :: wf_efield
  LOGICAL              :: wf_switch
  INTEGER              :: sw_len
  REAL(DP)             :: efx0, efy0, efz0
  REAL(DP)             :: efx1, efy1, efz1
  INTEGER              :: wfsd
  REAL(DP)             :: wfdt
!==============================================================
!exx_wf related
  INTEGER              :: neigh
  INTEGER              :: vnbsp 
  REAL(DP)             :: poisson_eps
  REAL(DP)             :: dis_cutoff
  REAL(DP)             :: exx_ps_rcut_s
  REAL(DP)             :: exx_ps_rcut_p
  REAL(DP)             :: exx_me_rcut_s
  REAL(DP)             :: exx_me_rcut_p
  LOGICAL              :: texx_cube
  !!
!==============================================================
  REAL(DP)             :: maxwfdt
  REAL(DP)             :: wf_q
  REAL(DP)             :: wf_friction
  INTEGER              :: nit
  INTEGER              :: nsd
  INTEGER              :: nsteps
  REAL(DP)             :: tolw
  LOGICAL              :: adapt
  INTEGER              :: calwf
  INTEGER              :: nwf
  INTEGER              :: wffort
  LOGICAL              :: writev
  INTEGER, ALLOCATABLE :: iplot(:)
  !
  ! ... other internal variables
  !
  INTEGER                  :: nw, nwrwf, iwf, jwf
  INTEGER,     ALLOCATABLE :: wfg1(:), wfg(:,:)
  INTEGER,     ALLOCATABLE :: indexplus(:,:), indexminus(:,:)
  INTEGER,     ALLOCATABLE :: indexplusz(:), indexminusz(:)
  INTEGER,     ALLOCATABLE :: tag(:,:), tagp(:,:)
  REAL(DP),    ALLOCATABLE :: weight(:)            ! weights of G vectors
  REAL(DP),    ALLOCATABLE :: gnx(:,:)
  INTEGER ,    ALLOCATABLE :: gnn(:,:)
  COMPLEX(DP), ALLOCATABLE :: expo(:,:)
  !
  ! ... reusable scratch buffers to prevent memory fragmentation
  !
  REAL(DP),    ALLOCATABLE :: becwf(:,:)
  COMPLEX(DP), ALLOCATABLE :: cwf(:,:), bec2(:), bec3(:), bec2up(:)
  COMPLEX(DP), ALLOCATABLE :: bec2dw(:), bec3up(:), bec3dw(:)
  COMPLEX(DP), ALLOCATABLE :: c_m(:,:), c_p(:,:), c_psp(:,:)
  COMPLEX(DP), ALLOCATABLE :: c_msp(:,:)
  INTEGER,     ALLOCATABLE :: tagz(:)
  REAL(DP),    ALLOCATABLE :: Uspin(:,:)
  COMPLEX(DP), ALLOCATABLE :: X(:,:), Xsp(:,:), X2(:,:), X3(:,:)
  COMPLEX(DP), ALLOCATABLE :: O(:,:,:), Ospin(:,:,:), Oa(:,:,:)
  COMPLEX(DP), ALLOCATABLE :: qv(:), fg1(:)
  REAL(DP),    ALLOCATABLE :: gr(:,:), mt(:), mt0(:), wr(:), W(:,:), EW(:,:)
  INTEGER,     ALLOCATABLE :: f3(:), f4(:)
  COMPLEX(DP), ALLOCATABLE :: U2(:,:)
#if defined (__MPI)
  COMPLEX(DP), ALLOCATABLE :: psitot(:,:), psitot_pl(:,:)
  COMPLEX(DP), ALLOCATABLE :: psitot_mi(:,:)
  INTEGER,     ALLOCATABLE :: ns(:)
#endif
  !
  CONTAINS
    !
    !------------------------------------------------------------------------
!=============================================================================
!exx_wf related
    SUBROUTINE wannier_init( wf_efield_, wf_switch_, sw_len_, efx0_, efy0_, &
                             efz0_, efx1_, efy1_, efz1_, wfsd_, wfdt_,      &
                             neigh_, poisson_eps_ ,dis_cutoff_, exx_ps_rcut_s_, exx_me_rcut_s_,&
                             exx_ps_rcut_p_, exx_me_rcut_p_, texx_cube_,    &
                             vnbsp_, maxwfdt_, wf_q_, wf_friction_, nit_, nsd_,     &
                             nsteps_, tolw_, adapt_, calwf_, nwf_, wffort_, &
                             writev_, iplot_, restart_mode_ )
!=============================================================================
      !------------------------------------------------------------------------
      !
      IMPLICIT NONE
      !
      LOGICAL,          INTENT(IN) :: wf_efield_
      LOGICAL,          INTENT(IN) :: wf_switch_
      INTEGER,          INTENT(IN) :: sw_len_
      REAL(DP),         INTENT(IN) :: efx0_, efy0_, efz0_
      REAL(DP),         INTENT(IN) :: efx1_, efy1_, efz1_
      INTEGER,          INTENT(IN) :: wfsd_
      REAL(DP),         INTENT(IN) :: wfdt_
!=============================================================================
!exx_wf related
      INTEGER,          INTENT(IN) :: neigh_
      INTEGER,          INTENT(IN) :: vnbsp_
      REAL(DP),         INTENT(IN) :: poisson_eps_
      REAL(DP),         INTENT(IN) :: dis_cutoff_
      REAL(DP),         INTENT(IN) :: exx_ps_rcut_s_
      REAL(DP),         INTENT(IN) :: exx_me_rcut_s_
      REAL(DP),         INTENT(IN) :: exx_ps_rcut_p_
      REAL(DP),         INTENT(IN) :: exx_me_rcut_p_
      LOGICAL,          INTENT(IN) :: texx_cube_
      !!
!=============================================================================
      REAL(DP),         INTENT(IN) :: maxwfdt_
      REAL(DP),         INTENT(IN) :: wf_q_
      REAL(DP),         INTENT(IN) :: wf_friction_
      INTEGER,          INTENT(IN) :: nit_
      INTEGER,          INTENT(IN) :: nsd_
      INTEGER,          INTENT(IN) :: nsteps_
      REAL(DP),         INTENT(IN) :: tolw_
      LOGICAL,          INTENT(IN) :: adapt_
      INTEGER,          INTENT(IN) :: calwf_
      INTEGER,          INTENT(IN) :: nwf_
      INTEGER,          INTENT(IN) :: wffort_
      INTEGER,          INTENT(IN) :: iplot_(:)
      LOGICAL,          INTENT(IN) :: writev_
      CHARACTER(LEN=*), INTENT(IN) :: restart_mode_
      !
      !
      wf_efield   = wf_efield_
      wf_switch   = wf_switch_
      sw_len      = sw_len_
      efx0        = efx0_
      efy0        = efy0_
      efz0        = efz0_
      efx1        = efx1_
      efy1        = efy1_
      efz1        = efz1_
      wfsd        = wfsd_
      wfdt        = wfdt_
!==================================================================
!exx_wf related
      neigh       = neigh_
      vnbsp       = vnbsp_
      poisson_eps = poisson_eps_
      dis_cutoff  = dis_cutoff_
      exx_ps_rcut_s = exx_ps_rcut_s_
      exx_me_rcut_s = exx_me_rcut_s_
      exx_ps_rcut_p = exx_ps_rcut_p_
      exx_me_rcut_p = exx_me_rcut_p_
      texx_cube = texx_cube_
!==================================================================
      maxwfdt     = maxwfdt_
      wf_q        = wf_q_
      wf_friction = wf_friction_
      nit         = nit_
      nsd         = nsd_
      nsteps      = nsteps_
      tolw        = tolw_
      adapt       = adapt_
      calwf       = calwf_
      nwf         = nwf_
      wffort      = wffort_
      writev      = writev_
      !
      IF ( calwf == 1 .AND. nwf == 0 ) &
         CALL errore( 'wannier_init ', &
                    & 'when calwf = 1, nwf must be larger that 0', 1 )
      !
      IF ( nwf > 0 ) THEN
         !
         ALLOCATE( iplot( nwf ) )
         !
         iplot(:) = iplot_(1:nwf)
         !
      END IF
      !
      IF ( TRIM( restart_mode_ ) == "from_scratch" ) THEN
         !
         IF ( wf_efield ) &
            CALL errore( 'wannier_init', 'electric field not ' // &
                       & 'allowed when starting from scratch', 1 )
         !
      END IF
      !
    END SUBROUTINE wannier_init
    !
    !
    !
    SUBROUTINE deallocate_wannier_base()
       IF( ALLOCATED( iplot ) ) DEALLOCATE( iplot )
       IF( ALLOCATED( wfg1 ) ) DEALLOCATE( wfg1 )
       IF( ALLOCATED( wfg ) ) DEALLOCATE( wfg )
       IF( ALLOCATED( indexplus ) ) DEALLOCATE( indexplus )
       IF( ALLOCATED( indexminus ) ) DEALLOCATE( indexminus )
       IF( ALLOCATED( indexplusz ) ) DEALLOCATE( indexplusz )
       IF( ALLOCATED( indexminusz ) ) DEALLOCATE( indexminusz )
       IF( ALLOCATED( tag ) ) DEALLOCATE( tag )
       IF( ALLOCATED( tagp ) ) DEALLOCATE( tagp )
       IF( ALLOCATED( weight ) ) DEALLOCATE( weight )
       IF( ALLOCATED( gnx ) ) DEALLOCATE( gnx )
       IF( ALLOCATED( gnn ) ) DEALLOCATE( gnn )
       IF( ALLOCATED( expo ) ) DEALLOCATE( expo )
       !
       IF (ALLOCATED(becwf))  DeALLOCATE(becwf)
       IF (ALLOCATED(U2))     DeALLOCATE(U2)
       IF (ALLOCATED(cwf))    DeALLOCATE(cwf)
       IF (ALLOCATED(bec2))   DeALLOCATE(bec2)
       IF (ALLOCATED(bec3))   DeALLOCATE(bec3)
       IF (ALLOCATED(bec2up)) DeALLOCATE(bec2up)
       IF (ALLOCATED(bec3up)) DeALLOCATE(bec3up)
       IF (ALLOCATED(bec2dw)) DEALLOCATE(bec2dw)
       IF (ALLOCATED(bec3dw)) DEALLOCATE(bec3dw)
       IF (ALLOCATED(tagz))   DEALLOCATE(tagz)
       IF (ALLOCATED(O))      DEALLOCATE(O)
       IF (ALLOCATED(X))      DEALLOCATE(X)
       IF (ALLOCATED(Oa))     DEALLOCATE(Oa)
       IF (ALLOCATED(X2))     DEALLOCATE(X2)
       IF (ALLOCATED(X3))     DEALLOCATE(X3)
       IF (ALLOCATED(c_p))    DEALLOCATE(c_p)
       IF (ALLOCATED(c_m))    DEALLOCATE(c_m)
       IF (ALLOCATED(qv))     DEALLOCATE(qv)
       IF (ALLOCATED(fg1))    DEALLOCATE(fg1)
       IF (ALLOCATED(wr))     DEALLOCATE(wr)
       IF (ALLOCATED(W))      DEALLOCATE(W)
       IF (ALLOCATED(gr))     DEALLOCATE(gr)
       IF (ALLOCATED(EW))     DEALLOCATE(EW)
       IF (ALLOCATED(f3))     DEALLOCATE(f3)
       IF (ALLOCATED(f4))     DEALLOCATE(f4)
       IF (ALLOCATED(mt0))    DEALLOCATE(mt0)
       IF (ALLOCATED(mt))     DEALLOCATE(mt)
       IF (ALLOCATED(Xsp))    DEALLOCATE(Xsp)
       IF (ALLOCATED(c_psp))  DEALLOCATE(c_psp)
       IF (ALLOCATED(c_msp))  DEALLOCATE(c_msp)
#if defined (__MPI)
       IF (ALLOCATED(ns))        DEALLOCATE(ns)
       IF (ALLOCATED(psitot))    DEALLOCATE(psitot)
       IF (ALLOCATED(psitot_pl)) DEALLOCATE(psitot_pl)
       IF (ALLOCATED(psitot_mi)) DEALLOCATE(psitot_mi)
#endif
       RETURN
    END SUBROUTINE deallocate_wannier_base
    !
    !
    !
    SUBROUTINE wannier_base_resize_scratch_only_once(nbsp)
       !! This routine is needed due to a potential compiler bug in NVHPC (present on 21.5)
       !! see https://gitlab.com/QEF/q-e/-/issues/346 for more detail
       INTEGER, INTENT(IN) :: nbsp
       IF (.NOT.ALLOCATED(X)) ALLOCATE(X(nbsp, nbsp))
    END SUBROUTINE wannier_base_resize_scratch_only_once
    !
    !
    !
END MODULE wannier_base
