PROGRAM spectral_fn_broadening
!-------------------------------------------------------------------------
!! authors: Marios Zacharias, Feliciano Giustino 
!! acknowledgement: Sabyasachi Tiwari for help packaging this release
!! version: v1.2
!! license: GNU
!
! This routine is used to calculate the spectral function, or momentum resolved
! density of states, from a band structure unfolding calculation. If you are 
! using this routine please cite:
!
! ! https://doi.org/10.1103/PhysRevResearch.2.013357
!
USE kinds,       ONLY : dp
USE mp_global,   ONLY : inter_pool_comm
USE mp_world,    ONLY : world_comm
USE mp,          ONLY : mp_bcast, mp_barrier, mp_sum
USE io_global,   ONLY : ionode, ionode_id, stdout
USE mp_global,   ONLY : mp_startup, mp_global_end
USE environment, ONLY : environment_start, environment_end
USE constants,   ONLY : pi
!
IMPLICIT NONE
!
  CHARACTER(LEN=256)    :: flspfn, flin
  INTEGER               :: ksteps, esteps, steps, ik, ie, i, ios
  INTEGER               :: lower_bnd, upper_bnd
  REAL(DP)              :: kmin, kmax, emin, emax
  REAL(DP)              :: maxv, sigmak, sigmae, argk, arge
  REAL(DP), ALLOCATABLE :: k(:), e(:), k0(:), e0(:), w0(:), spctrlfn(:, :)
  !
  NAMELIST /input/ steps, ksteps, esteps, kmin, kmax, emin, emax, &
  &                flspfn, flin
  !
CALL mp_startup()
CALL environment_start('SPFN_BROAD')
!
  IF (ionode) CALL input_from_file ( )
  !
  ! ... all calculations are done by the first cpu
  !
  ! set namelist default
  !
  steps  = 5000
  ksteps = 200
  esteps = 200
  kmin   = 0.00000
  kmax   = 1.00000
  emin   = -5.00
  emax   = 5.00
  !
  IF (ionode) READ (5,input,IOSTAT=ios)
  CALL mp_bcast(ios, ionode_id, world_comm)
  CALL errore('spctrlfn', 'reading input namelist', ABS(ios))
  CALL mp_bcast(steps, ionode_id, world_comm)
  CALL mp_bcast(ksteps, ionode_id, world_comm)
  CALL mp_bcast(esteps, ionode_id, world_comm)
  CALL mp_bcast(kmin, ionode_id, world_comm)
  CALL mp_bcast(kmax, ionode_id, world_comm)
  CALL mp_bcast(emin, ionode_id, world_comm)
  CALL mp_bcast(emax, ionode_id, world_comm)
  CALL mp_bcast(flspfn, ionode_id, world_comm)
  CALL mp_bcast(flin, ionode_id, world_comm)

ALLOCATE(k(ksteps), e(esteps), k0(steps), e0(steps), w0(steps), spctrlfn(ksteps, esteps))

sigmak = (kmax - kmin) / float(ksteps)
sigmae = (emax - emin) / float(esteps)

DO ik = 1, ksteps
  k(ik) = kmin + (kmax - kmin) / float(ksteps - 1) * float(ik - 1)
ENDDO
!
DO ie = 1, esteps
  e(ie) = emin + (emax - emin) / float(esteps - 1) * float(ie - 1)
ENDDO

IF (ionode) THEN
  OPEN(UNIT=20, FILE=flin)
  !
  DO i = 1, steps
    READ(20,*) k0(i), e0(i), w0(i)
  ENDDO
  !
  CLOSE(20)
ENDIF
!
CALL mp_bcast(k0, ionode_id, world_comm)
CALL mp_bcast(e0, ionode_id, world_comm)
CALL mp_bcast(w0, ionode_id, world_comm)
!
spctrlfn = 0.d0
!
CALL fkbounds( ksteps, lower_bnd, upper_bnd )
!
DO ik = lower_bnd, upper_bnd
  DO ie = 1, esteps
    DO i = 1, steps
      argk = -(k(ik) - k0(i))**2.d0 / 2.d0 / sigmak**2.d0
      arge = -(e(ie) - e0(i))**2.d0 / 2.d0 / sigmae**2.d0
      spctrlfn(ik, ie) = spctrlfn(ik, ie) + w0(i) * exp(argk) * exp(arge)
    ENDDO
    spctrlfn(ik, ie) = spctrlfn(ik, ie) / (2 * pi) / sigmak / sigmae
  ENDDO
ENDDO
!
CALL mp_sum( spctrlfn, inter_pool_comm)
CALL mp_barrier( inter_pool_comm)
!
IF (ionode) THEN
  maxv = maxval(spctrlfn)
  OPEN(10, FILE = flspfn)
  DO ik = 1, ksteps
    DO ie = 1, esteps
      WRITE(10,'(3F16.8)') k(ik), e(ie), spctrlfn(ik, ie) / maxv
    ENDDO
    WRITE(10,*)
  ENDDO
  CLOSE(10)
ENDIF
!
DEALLOCATE(k, e, k0, e0, w0, spctrlfn)
!
CALL environment_end('SPFN_BROAD')
!
CALL mp_global_end()
!
STOP
!
END PROGRAM

SUBROUTINE fkbounds( nktot, lower_bnd, upper_bnd )
  !-----------------------------------------------------------------------
  !!
  !!   Subroutine from EPW finds the lower and upper bounds a k-grid in
  !parallel
  !!
  !! @ Note: 
  !!    If you have 19 kpts and 2 pool, this routine will RETURN
  !!    lower_bnd= 1 and upper_bnd = 10 for the first pool
  !!    lower_bnd= 1 and upper_bnd = 9 for the second pool
  !-----------------------------------------------------------------------
  !
  USE mp_global,    ONLY: my_pool_id, npool
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT (in) :: nktot
  !! nktot k-points splited over pools
  INTEGER, INTENT (out) :: lower_bnd
  !! Lower kpt bounds for that image pool 
  INTEGER, INTENT (out) :: upper_bnd
  !! Upper kpt for that image pool
  !
#if defined(__MPI)
  !
  INTEGER :: nkl, nkr
  !
  ! find the bounds of k-dependent arrays in the parallel case
  ! number of kpoint blocks, kpoints per pool and reminder
  !
  nkl = nktot / npool
  nkr = nktot - nkl * npool
  !
  ! the reminder goes to the first nkr pools (0...nkr-1)
  !
  IF (my_pool_id < nkr ) nkl = nkl + 1
  !
  ! the index of the first k point in this pool
  !
  lower_bnd = my_pool_id * nkl + 1
  IF ( my_pool_id >= nkr ) lower_bnd = my_pool_id * nkl + 1 + nkr
  !
  ! the index of the last k point in this pool
  !
  upper_bnd = lower_bnd + nkl - 1
  !
#else  
  !     
  ! In serial the definitions are much easier 
  !     
  lower_bnd = 1
  upper_bnd = nktot
  !     
#endif 
  !
  RETURN
  !
  END SUBROUTINE fkbounds
