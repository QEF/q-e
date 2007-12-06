!
! Copyright (C) 2001-2007 Quantum-Espresso group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
PROGRAM do_plan_avg
  !-----------------------------------------------------------------------
  !
  ! calculate planar averages of each wavefunction
  !
  USE kinds,     ONLY : DP
  USE char,      ONLY : title
  USE cell_base, ONLY : ibrav, celldm, at
  USE gvect,     ONLY : nrx1, nrx2, nrx3, nr1, nr2, nr3, gcutm, dual, ecutwfc
  USE klist,     ONLY : nkstot, xk
  USE ions_base, ONLY : nat, ntyp=>nsp, ityp, tau, atm, zv
  USE io_files,  ONLY : tmp_dir, prefix, nd_nmbr, trimcheck
  USE io_global, ONLY : ionode, ionode_id
  USE wvfct,     ONLY : nbnd
  USE mp,        ONLY : mp_bcast
  USE control_flags, ONLY : gamma_only
  !
  IMPLICIT NONE
  INTEGER :: ninter
  CHARACTER(len=256) :: filplot, outdir
  REAL(DP), ALLOCATABLE :: averag (:,:,:), plan (:,:,:)
  !
  INTEGER :: iunplot = 4, ios, ibnd, ik, ir, nt, na, i
  !
  NAMELIST / inputpp / outdir, prefix, filplot
  !
  ! initialise parallel environment
  !
  CALL start_postproc (nd_nmbr)
  IF ( ionode )  CALL input_from_file ( )
  !
  !
  !   set default values for variables in namelist
  !
  prefix = 'pwscf'
  CALL get_env( 'ESPRESSO_TMPDIR', outdir )
  IF ( TRIM( outdir ) == ' ' ) outdir = './'
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
  END IF
  !
  CALL mp_bcast( ios, ionode_id )
  IF ( ios /= 0 ) CALL errore ('plan_avg', 'reading inputpp namelist', ABS(ios))
  !
  ! ... Broadcast variables
  !
  CALL mp_bcast( tmp_dir, ionode_id )
  CALL mp_bcast( prefix, ionode_id )
  CALL mp_bcast( filplot, ionode_id )
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
  ALLOCATE (plan(nr3, nbnd, nkstot))    
  !
  CALL plan_avg (averag, plan, ninter)
  !
  IF ( ionode ) THEN
     !
     OPEN (UNIT = iunplot, FILE = filplot, FORM = 'formatted', &
          STATUS = 'unknown', err = 100, IOSTAT = ios)
100  CALL errore ('plan_avg', 'opening file '//TRIM(filplot), abs (ios) )
     WRITE (iunplot, '(a)') title
     WRITE (iunplot, '(8i8)') nrx1, nrx2, nrx3, nr1, nr2, nr3, nat, ntyp
     WRITE (iunplot, '(i6,6f12.8)') ibrav, celldm
     IF  (ibrav == 0) THEN
        WRITE ( iunplot, * ) at(:,1)
        WRITE ( iunplot, * ) at(:,2)
        WRITE ( iunplot, * ) at(:,3)
     END IF
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
           DO ir = 1, nr3
              WRITE (iunplot, * ) ir, plan (ir, ibnd, ik)
           ENDDO
        ENDDO
     ENDDO
     !
  ENDIF
  !
  DEALLOCATE (plan)
  DEALLOCATE (averag)

END PROGRAM do_plan_avg
!
subroutine plan_avg (averag, plan, ninter)
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
#include "f_defs.h"
  USE cell_base, ONLY: celldm, omega, alat, tpiba2
  USE ions_base, ONLY: nat, ntyp=>nsp, ityp, tau
  USE gvect
  USE klist, ONLY: nks, nkstot, xk
  USE lsda_mod, ONLY: lsda, current_spin, isk
  USE uspp, ONLY: vkb, nkb
  USE wvfct, ONLY: npw, npwx, nbnd, wg, igk, g2kin
  USE wavefunctions_module,  ONLY: evc
  USE noncollin_module, ONLY : noncolin, npol
  USE io_files, ONLY: iunwfc, nwordwfc
  USE becmod, ONLY: becp, becp_nc, calbec

  implicit none
  integer :: ninter
  ! output: the number of planes
  real(DP) :: averag (nat, nbnd, nkstot), plan (nr3, nbnd, nkstot)
  ! output: the average charge on ea
  ! output: the planar average
  !
  !      Local variables
  !
  integer :: ik, ibnd, iin, na, ir, ij, ind, i1 (nat), ntau (nat + 1)
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

  if ( celldm(3) == 0.d0 ) celldm(3) = celldm(1)
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
  do na = 2, nat
     do iin = 1, ninter
        if (abs (mod (z1(iin)-tau(3,na), celldm(3)) ) .lt. sp_min) then
           avg (iin) = avg (iin) + tau (3, na)
           ntau (iin) = ntau (iin) + 1
           goto 100
        endif
     enddo
     ninter = ninter + 1
     z1 (ninter) = tau (3, na)
     avg (ninter) = tau (3, na)
     ntau (ninter) = 1
100  continue
  enddo
  !
  !     for each plane compute the average position of the central plane
  !     and first point in the fft mesh
  !
  do iin = 1, ninter
     z1 (iin) = mod (avg (iin), celldm (3) ) / ntau (iin)
     ind = (z1 (iin) / celldm (3) ) * nr3 + 1
     if (ind.le.0) ind = ind+nr3
     i1 (iin) = ind
  enddo
  !
  !    order the points
  !
  do iin = 1, ninter
     ntau (iin) = i1 (iin)
     do ik = iin + 1, ninter
        if (i1 (ik) .lt.ntau (iin) ) then
           ij = ntau (iin)
           ntau (iin) = i1 (ik)
           i1 (ik) = ij
        endif
     enddo
  enddo
  ntau (ninter + 1) = ntau (1) + nr3
  !
  !    and compute the point associated to each plane
  !
  do iin = 1, ninter
     i1 (iin) = (ntau (iin) + ntau (iin + 1) ) / 2
  enddo
  !
  !     for each state compute the planar average
  !
  averag(:,:,:) = 0.d0
  plan(:,:,:) = 0.d0
  if (noncolin) then
     allocate(becp_nc(nkb,npol,nbnd))
  else
     allocate(becp(nkb,nbnd))
  endif
  CALL init_us_1 ( )
  do ik = 1, nks
     if (lsda) current_spin = isk (ik)
     call gk_sort (xk (1, ik), ngm, g, ecutwfc / tpiba2, npw, igk, g2kin)
     call davcio (evc, nwordwfc, iunwfc, ik, - 1)
     call init_us_2 (npw, igk, xk (1, ik), vkb)

     if (noncolin) then
        call calbec ( npw, vkb, evc, becp_nc)
     else
        call calbec ( npw, vkb, evc, becp)
     endif
     do ibnd = 1, nbnd
        call local_dos1d (ik, ibnd, plan (1, ibnd, ik) )
        !
        !     compute the integrals of the charge
        !
        do ir = 1, i1 (1) - 1
           averag (1, ibnd, ik) = averag (1, ibnd, ik) + plan (ir, ibnd, ik)
        enddo
        do ir = i1 (ninter), nr3
           averag (1, ibnd, ik) = averag (1, ibnd, ik) + plan (ir, ibnd, ik)
        enddo
        averag (1, ibnd, ik) = averag (1, ibnd, ik) * zdim / nr3
        sum = averag (1, ibnd, ik)
        do iin = 2, ninter
           do ir = i1 (iin - 1), i1 (iin) - 1
              averag(iin,ibnd,ik) = averag(iin,ibnd,ik) + plan(ir,ibnd,ik)
           enddo
           averag (iin, ibnd, ik) = averag (iin, ibnd, ik) * zdim / nr3
           sum = sum + averag (iin, ibnd, ik)
        enddo
     enddo
  enddo
  if (noncolin) then
     deallocate (becp_nc)
  else
     deallocate (becp)
  endif
#ifdef __PARA
  call poolrecover (plan, nr3 * nbnd, nkstot, nks)
  call poolrecover (averag, nat * nbnd, nkstot, nks)
  call poolrecover (xk, 3, nkstot, nks)
#endif
  return
end subroutine plan_avg
