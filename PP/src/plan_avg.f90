!
! Copyright (C) 2001-2009 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
PROGRAM plan_avg
  !-----------------------------------------------------------------------
  !
  ! calculate planar averages of each wavefunction
  !
  USE kinds,     ONLY : DP
  USE run_info,  ONLY: title
  USE cell_base, ONLY : ibrav, celldm, at
  USE fft_base,  ONLY : dfftp
  USE gvect,     ONLY : gcutm
  USE gvecs,     ONLY : dual
  USE klist,     ONLY : nkstot, xk
  USE ions_base, ONLY : nat, ntyp=>nsp, ityp, tau, atm, zv
  USE io_files,  ONLY : tmp_dir, prefix
  USE io_global, ONLY : ionode, ionode_id
  USE wvfct,     ONLY : nbnd
  USE gvecw,     ONLY : ecutwfc
  USE mp,        ONLY : mp_bcast
  USE mp_world,  ONLY : world_comm
  USE mp_global, ONLY : mp_startup
  USE control_flags, ONLY : gamma_only
  USE environment,   ONLY : environment_start, environment_end
  !
  IMPLICIT NONE
  !
  CHARACTER(LEN=256), EXTERNAL :: trimcheck
  !
  INTEGER :: ninter
  CHARACTER(len=256) :: filplot, outdir
  REAL(DP), ALLOCATABLE :: averag (:,:,:), plan (:,:,:)
  !
  INTEGER :: iunplot = 4, ios, ibnd, ik, ir, nt, na, i
  !
  NAMELIST / inputpp / outdir, prefix, filplot
  !
  ! initialise environment
  !
#if defined(__MPI)
  CALL mp_startup ( )
#endif
  CALL environment_start ( 'plan-avg' )
  !
  IF ( ionode )  CALL input_from_file ( )
  !
  !
  !   set default values for variables in namelist
  !
  prefix = 'pwscf'
  CALL get_environment_variable( 'ESPRESSO_TMPDIR', outdir )
  IF ( trim( outdir ) == ' ' ) outdir = './'
  filplot = 'tmp.pp'
  !
  ios = 0
  !
  IF ( ionode )  THEN
     !
     !     reading the namelist inputpp
     !
     READ (5, inputpp, iostat = ios)
     tmp_dir = trimcheck (outdir)
     !
  ENDIF
  !
  CALL mp_bcast( ios, ionode_id, world_comm )
  IF ( ios /= 0 ) CALL errore ('plan_avg', 'reading inputpp namelist', abs(ios))
  !
  ! ... Broadcast variables
  !
  CALL mp_bcast( tmp_dir, ionode_id, world_comm )
  CALL mp_bcast( prefix, ionode_id, world_comm )
  CALL mp_bcast( filplot, ionode_id, world_comm )
  !
  !   Now allocate space for pwscf variables, read and check them.
  !
  CALL read_file ( )
  !
  IF (gamma_only) CALL errore ('plan_avg', &
       ' planar average with gamma tricks not yet implemented',2)
  !
  CALL openfil_pp ( )
  !
  ALLOCATE (averag( nat, nbnd, nkstot))
  ALLOCATE (plan(dfftp%nr3, nbnd, nkstot))
  !
  CALL do_plan_avg (averag, plan, ninter)
  !
  IF ( ionode ) THEN
     !
     OPEN (UNIT = iunplot, FILE = filplot, FORM = 'formatted', &
          STATUS = 'unknown', err = 100, IOSTAT = ios)
100  CALL errore ('plan_avg', 'opening file '//trim(filplot), abs (ios) )
     WRITE (iunplot, '(a)') title
     WRITE (iunplot, '(8i8)') dfftp%nr1x, dfftp%nr2x, dfftp%nr3x, dfftp%nr1, dfftp%nr2, dfftp%nr3, nat, ntyp
     WRITE (iunplot, '(i6,6f12.8)') ibrav, celldm
     IF  (ibrav == 0) THEN
        WRITE ( iunplot, * ) at(:,1)
        WRITE ( iunplot, * ) at(:,2)
        WRITE ( iunplot, * ) at(:,3)
     ENDIF
     WRITE (iunplot, '(3f20.10,i6)') gcutm, dual, ecutwfc, 9
     WRITE (iunplot, '(i4,3x,a2,3x,f5.2)') &
          (nt, atm (nt), zv (nt), nt=1, ntyp)
     WRITE (iunplot, '(i4,3x,3f15.9,3x,i2)') (na, &
          (tau (i, na), i = 1, 3), ityp (na), na = 1, nat)
     !
     WRITE (iunplot, '(3i8)') ninter, nkstot, nbnd
     DO ik = 1, nkstot
        DO ibnd = 1, nbnd
           WRITE (iunplot, '(3f15.9,i5)') xk (1, ik) , xk (2, ik) , xk (3, &
                ik) , ibnd
           WRITE (iunplot, '(4(1pe17.9))') (averag (ir, ibnd, ik) , ir = 1, &
                ninter)
           DO ir = 1, dfftp%nr3
              WRITE (iunplot, * ) ir, plan (ir, ibnd, ik)
           ENDDO
        ENDDO
     ENDDO
     !
     CLOSE (UNIT = iunplot, STATUS = 'keep')
     !
  ENDIF
  !
  DEALLOCATE (plan)
  DEALLOCATE (averag)
  !
  CALL environment_end ( 'plan-avg' )
  !
  CALL stop_pp ( )

CONTAINS
!
SUBROUTINE do_plan_avg (averag, plan, ninter)
  !
  !    This routine computes the planar average on the xy plane
  !    for the charge density of each state of the system.
  !    The routine should work on parallel machines.
  !    On these machines the results are collected for all
  !    k points and on exit each processor contains the
  !    planar average of all k points (even those of other pools).
  !    In the US case the augmentation part is added only in one
  !    dimension, so that no overload with respect to the NC case
  !    is expected.
  !
  !    Furthermore the amount of charge contained in each plane is
  !    evaluated and given as output. The number of planes is
  !    computed starting from the atomic positions
  !
  USE cell_base, ONLY: celldm, omega, alat
  USE ions_base, ONLY: nat, ntyp=>nsp, ityp, tau
  USE gvect
  USE klist, ONLY: nks, nkstot, xk, ngk, igk_k
  USE lsda_mod, ONLY: lsda, current_spin, isk
  USE uspp, ONLY: vkb, nkb
  USE wvfct, ONLY: npwx, nbnd, wg
  USE wavefunctions_module,  ONLY: evc
  USE noncollin_module, ONLY : noncolin, npol
  USE io_files, ONLY: iunwfc, nwordwfc
  USE becmod, ONLY: bec_type, becp, calbec, allocate_bec_type, deallocate_bec_type

  IMPLICIT NONE
  INTEGER :: ninter
  ! output: the number of planes
  real(DP) :: averag (nat, nbnd, nkstot), plan (dfftp%nr3, nbnd, nkstot)
  ! output: the average charge on ea
  ! output: the planar average
  !
  !      Local variables
  !
  INTEGER :: ik, ibnd, iin, na, ir, ij, ind, i1 (nat), ntau (nat + 1), npw
  ! counter on k points
  ! counter on bands
  ! counter on planes
  ! counter on atoms
  ! counter on points
  ! counter on coordinates and planes
  ! starting point of each plane
  ! the number of tau per plane

  real(DP) :: sp_min, avg (nat), z1 (nat), sum, zdim
  ! minimum plane distance
  ! the average position of each plane
  ! auxiliary for coordinates
  ! length in a.u. of the cell along z

  IF ( celldm(3) == 0.d0 ) celldm(3) = celldm(1)
  zdim = alat * celldm (3)
  sp_min = 2.d0 / alat
  !
  !     Compute the number of planes and the coordinates on the mesh of th
  !     points which define each plane
  !
  avg(:) = 0.d0
  ninter = 1
  z1 (ninter) = tau (3, 1)
  avg (ninter) = tau (3, 1)
  ntau (ninter) = 1
  DO na = 2, nat
     DO iin = 1, ninter
        IF (abs (mod (z1(iin)-tau(3,na), celldm(3)) ) < sp_min) THEN
           avg (iin) = avg (iin) + tau (3, na)
           ntau (iin) = ntau (iin) + 1
           GOTO 100
        ENDIF
     ENDDO
     ninter = ninter + 1
     z1 (ninter) = tau (3, na)
     avg (ninter) = tau (3, na)
     ntau (ninter) = 1
100  CONTINUE
  ENDDO
  !
  !     for each plane compute the average position of the central plane
  !     and first point in the fft mesh
  !
  DO iin = 1, ninter
     z1 (iin) = mod (avg (iin), celldm (3) ) / ntau (iin)
     ind = (z1 (iin) / celldm (3) ) * dfftp%nr3 + 1
     IF (ind<=0) ind = ind+dfftp%nr3
     i1 (iin) = ind
  ENDDO
  !
  !    order the points
  !
  DO iin = 1, ninter
     ntau (iin) = i1 (iin)
     DO ik = iin + 1, ninter
        IF (i1 (ik) <ntau (iin) ) THEN
           ij = ntau (iin)
           ntau (iin) = i1 (ik)
           i1 (ik) = ij
        ENDIF
     ENDDO
  ENDDO
  ntau (ninter + 1) = ntau (1) + dfftp%nr3
  !
  !    and compute the point associated to each plane
  !
  DO iin = 1, ninter
     i1 (iin) = (ntau (iin) + ntau (iin + 1) ) / 2
  ENDDO
  !
  !     for each state compute the planar average
  !
  averag(:,:,:) = 0.d0
  plan(:,:,:) = 0.d0
  CALL allocate_bec_type ( nkb, nbnd, becp )
  DO ik = 1, nks
     IF (lsda) current_spin = isk (ik)
     npw = ngk(ik)
     CALL davcio (evc, 2*nwordwfc, iunwfc, ik, - 1)
     CALL init_us_2 (npw, igk_k(1,ik), xk (1, ik), vkb)

     CALL calbec ( npw, vkb, evc, becp)

     DO ibnd = 1, nbnd
        CALL local_dos1d (ik, ibnd, plan (1, ibnd, ik) )
        !
        !     compute the integrals of the charge
        !
        DO ir = 1, i1 (1) - 1
           averag (1, ibnd, ik) = averag (1, ibnd, ik) + plan (ir, ibnd, ik)
        ENDDO
        DO ir = i1 (ninter), dfftp%nr3
           averag (1, ibnd, ik) = averag (1, ibnd, ik) + plan (ir, ibnd, ik)
        ENDDO
        averag (1, ibnd, ik) = averag (1, ibnd, ik) * zdim / dfftp%nr3
        sum = averag (1, ibnd, ik)
        DO iin = 2, ninter
           DO ir = i1 (iin - 1), i1 (iin) - 1
              averag(iin,ibnd,ik) = averag(iin,ibnd,ik) + plan(ir,ibnd,ik)
           ENDDO
           averag (iin, ibnd, ik) = averag (iin, ibnd, ik) * zdim / dfftp%nr3
           sum = sum + averag (iin, ibnd, ik)
        ENDDO
     ENDDO
  ENDDO
  CALL deallocate_bec_type (becp)
#if defined(__MPI)
  CALL poolrecover (plan, dfftp%nr3 * nbnd, nkstot, nks)
  CALL poolrecover (averag, nat * nbnd, nkstot, nks)
  CALL poolrecover (xk, 3, nkstot, nks)
#endif
  RETURN
END SUBROUTINE do_plan_avg

END PROGRAM plan_avg

