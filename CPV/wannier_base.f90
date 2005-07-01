!
! Copyright (C) 2002-2005 Quantum-ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
MODULE wannier_base
  !----------------------------------------------------------------------------
  !
  USE kinds, ONLY : dbl
  !
  IMPLICIT NONE
  !
  ! ... input variables
  !
  LOGICAL        :: wf_efield
  LOGICAL        :: wf_switch
  INTEGER        :: sw_len
  REAL(KIND=dbl) :: efx0, efy0, efz0
  REAL(KIND=dbl) :: efx1, efy1, efz1
  LOGICAL        :: wfsd
  REAL(KIND=dbl) :: wfdt
  REAL(KIND=dbl) :: maxwfdt
  REAL(KIND=dbl) :: wf_q
  REAL(KIND=dbl) :: wf_dt
  REAL(KIND=dbl) :: wf_friction
  INTEGER        :: nit
  INTEGER        :: nsd
  INTEGER        :: nsteps
  REAL(KIND=dbl) :: tolw
  LOGICAL        :: adapt
  INTEGER        :: calwf
  INTEGER        :: nwf
  INTEGER        :: wffort
  INTEGER        :: iwf
  LOGICAL        :: writev
  !
  ! ... other internal variables
  !
  INTEGER                        :: nw, nwrwf ,jwf
  INTEGER, ALLOCATABLE           :: iplot(:), wfg1(:), wfg(:,:)
  INTEGER, ALLOCATABLE           :: indexplus(:,:), indexminus(:,:)
  INTEGER, ALLOCATABLE           :: indexplusz(:), indexminusz(:)
  INTEGER, ALLOCATABLE           :: tag(:,:), tagp(:,:)
  REAL(KIND=dbl),    ALLOCATABLE :: weight(:)            ! weights of G vectors
  REAL(KIND=dbl),    ALLOCATABLE :: gnx(:,:)
  INTEGER ,          ALLOCATABLE :: gnn(:,:)
  COMPLEX(KIND=dbl), ALLOCATABLE :: expo(:,:)
  !
  CONTAINS
    !
    !------------------------------------------------------------------------
    SUBROUTINE wannier_init( wf_efield_, wf_switch_, sw_len_, efx0_, efy0_, &
                             efz0_, efx1_, efy1_, efz1_, wfsd_, wfdt_,      &
                             maxwfdt_, wf_q_, wf_dt_, wf_friction_, nit_,   &
                             nsd_, nsteps_, tolw_, adapt_, calwf_, nwf_,    &
                             wffort_, iwf_, writev_, restart_mode_ )
      !------------------------------------------------------------------------
      !
      IMPLICIT NONE
      !
      LOGICAL,          INTENT(IN) :: wf_efield_
      LOGICAL,          INTENT(IN) :: wf_switch_
      INTEGER,          INTENT(IN) :: sw_len_
      REAL(KIND=dbl),   INTENT(IN) :: efx0_, efy0_, efz0_
      REAL(KIND=dbl),   INTENT(IN) :: efx1_, efy1_, efz1_
      LOGICAL,          INTENT(IN) :: wfsd_
      REAL(KIND=dbl),   INTENT(IN) :: wfdt_
      REAL(KIND=dbl),   INTENT(IN) :: maxwfdt_
      REAL(KIND=dbl),   INTENT(IN) :: wf_q_
      REAL(KIND=dbl),   INTENT(IN) :: wf_dt_
      REAL(KIND=dbl),   INTENT(IN) :: wf_friction_
      INTEGER,          INTENT(IN) :: nit_
      INTEGER,          INTENT(IN) :: nsd_
      INTEGER,          INTENT(IN) :: nsteps_
      REAL(KIND=dbl),   INTENT(IN) :: tolw_
      LOGICAL,          INTENT(IN) :: adapt_
      INTEGER,          INTENT(IN) :: calwf_
      INTEGER,          INTENT(IN) :: nwf_
      INTEGER,          INTENT(IN) :: wffort_
      INTEGER,          INTENT(IN) :: iwf_
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
      maxwfdt     = maxwfdt_
      wf_q        = wf_q_
      wf_dt       = wf_dt_
      wf_friction = wf_friction_
      nit         = nit_
      nsd         = nsd_
      nsteps      = nsteps_
      tolw        = tolw_
      adapt       = adapt_
      calwf       = calwf_
      nwf         = nwf_
      wffort      = wffort_
      iwf         = iwf_
      writev      = writev_
      !
      IF ( TRIM( restart_mode_ ) == "from_scratch" ) THEN
         !
         IF ( wf_efield == .TRUE.  ) &
            CALL errore( 'wannier_init ', &
                       & 'electric field not allowed from scratch', 1 )
         !
      END IF
      !
    END SUBROUTINE wannier_init
    !
END MODULE wannier_base
