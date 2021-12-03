PROGRAM disca_broadening
!-------------------------------------------------------------------------
!! authors: Marios Zacharias, Feliciano Giustino 
!! acknowledgement: Hyungjun Lee for help packaging this release
!! version: v0.1
!! license: GNU
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
  CHARACTER(LEN=256)     :: flstrfin, flstrfout
  INTEGER                :: kres1, kres2, ik, iky, col1, col2
  INTEGER                :: steps, ii, lower_bnd, upper_bnd, Np, ios
  REAL(DP)               :: kmin, kmax, maxv
  REAL(DP)               :: jump, sf_smearingx, sf_smearingy !, pi 
  REAL(DP), ALLOCATABLE  :: kgridx(:), kgridy(:), structure_fact(:, :), structure_fact_out(:, :)
!
NAMELIST /input/ steps, kres1, kres2, kmin, kmax, &
&                col1, col2, Np, flstrfin, flstrfout 
!
CALL mp_startup()
CALL environment_start('DISCA_BROADENING')
!
  IF (ionode) CALL input_from_file ( )
  !
  ! ... all calculations are done by the first cpu
  !
  ! set namelist default
  !
  steps = 10000
  kres1 = 250
  kres2 = 250 
  kmin = -5
  kmax = 5
  col1 = 1
  col2 = 2
  Np = 400000
  flstrfin = 'input_data.dat'
  flstrfout = 'output_data.dat'
  ! 
  IF (ionode) READ (5,input,IOSTAT=ios)
  CALL mp_bcast(ios, ionode_id, world_comm)
  CALL errore('disca', 'reading input namelist', ABS(ios))
  CALL mp_bcast(steps, ionode_id, world_comm)
  CALL mp_bcast(kres1, ionode_id, world_comm)
  CALL mp_bcast(kres2, ionode_id, world_comm)
  CALL mp_bcast(kmin, ionode_id, world_comm)
  CALL mp_bcast(kmax, ionode_id, world_comm)
  CALL mp_bcast(col1, ionode_id, world_comm)
  CALL mp_bcast(col2, ionode_id, world_comm)
  CALL mp_bcast(Np, ionode_id, world_comm)
  CALL mp_bcast(flstrfin, ionode_id, world_comm)
  CALL mp_bcast(flstrfout, ionode_id, world_comm)


!OPEN(46,FILE=flstrfin)
!    READ(46,*) steps
!    READ(46,*) kres1
!    READ(46,*) kres2
!    READ(46,*) kmin
!    READ(46,*) kmax
!    READ(46,*) col1
!    READ(46,*) col2
!    READ(46,*) Np
!CLOSE(46)

ALLOCATE(structure_fact(steps,4), kgridx(kres1), kgridy(kres2))
ALLOCATE(structure_fact_out(kres1,kres2))
!kmin = -10.0
!kmax = 10.0

jump = (kmax - kmin) / dble(kres1 - 1)
DO ik = 1, kres1
  kgridx(ik) = kmin + (ik - 1) * jump
ENDDO
sf_smearingx = (kmax - kmin)/dble(kres1)
!
jump = (kmax - kmin) / dble(kres2 - 1)
DO ik = 1, kres2
  kgridy(ik) = kmin + (ik- 1)*jump
ENDDO
sf_smearingy = (kmax-kmin)/dble(kres2)
!! 
OPEN(46, FILE = flstrfin)
  DO ii = 1, steps, 1
    READ(46, *)  structure_fact(ii,:)
  ENDDO ! ii
CLOSE(46)

structure_fact_out = 0.d0
!
CALL fkbounds( steps, lower_bnd, upper_bnd )
!
DO ii = lower_bnd, upper_bnd
  DO ik = 1, kres1 !
    DO iky = 1, kres2
    !
    structure_fact_out(ik, iky) =  structure_fact_out(ik, iky) +  &
                          structure_fact(ii, 4) / sf_smearingx / sqrt(2.0d0 * pi) / sf_smearingy / sqrt(2.0d0 * pi)* &
                          (EXP(-(structure_fact(ii, col1)- kgridx(ik))**2.d0 / sf_smearingx**2.d0 / 2.d0))*&
                          (EXP(-(structure_fact(ii, col2)- kgridy(iky))**2.d0 / sf_smearingy**2.d0 / 2.d0))
    !
    ENDDO
  ENDDO
ENDDO
!
CALL mp_sum(structure_fact_out, inter_pool_comm)
CALL mp_barrier(inter_pool_comm)
!
IF (ionode) THEN
  OPEN(46,FILE=flstrfout)
  ! maxv = DBLE(6.162035539820569E+019)
  maxv =  maxval(structure_fact_out)
  WRITE(46,*) "#", maxv, maxval(structure_fact_out)
  DO ik = 1, kres1
    DO iky = 1, kres2
      WRITE(46,'(3F28.12)') kgridx(ik), kgridy(iky), structure_fact_out(ik,iky) * Np**(-2.0d0) !1.169035831194574E+019 
                                         !/ !maxval(abs(structure_fact_out)) !!* Np**(-2.0d0) ! / maxv
    ENDDO
    WRITE(46,*)
  ENDDO
  CLOSE(46)
ENDIF
!
!
!
DEALLOCATE(structure_fact, structure_fact_out, kgridx, kgridy)
CALL environment_end('DISCA_BROADENING')
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
  !!    lower_bnd= 1 and upper_bnd= 10 for the first pool
  !!    lower_bnd= 1 and upper_bnd=9 for the second pool
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
  ! the reminder goes to the first nkr pools (0...nkr- 1)
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
