!
! Copyright (C) 2001-2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "f_defs.h"
!
!-----------------------------------------------------------------------
PROGRAM bands
  !-----------------------------------------------------------------------
  !
  USE io_files,  ONLY : nd_nmbr, prefix, tmp_dir
  USE mp_global, ONLY : npool
  USE io_global, ONLY : ionode, ionode_id
  USE mp,        ONLY : mp_bcast
  !
  IMPLICIT NONE
  !
  CHARACTER (len=256) :: filband
  CHARACTER (len=256) :: outdir
  INTEGER :: ios
  !
  NAMELIST / inputpp / outdir, prefix, filband
  !                                  
  CHARACTER (LEN=256) :: input_file
  INTEGER             :: nargs, iiarg, ierr, ILEN
  INTEGER, EXTERNAL   :: iargc
  !
  !
  CALL start_postproc (nd_nmbr)
  !
  !   set default values for variables in namelist
  !
  prefix = 'pwscf'
  outdir = './'
  filband = 'bands.out'
  !
  IF ( npool > 1 ) CALL errore('bands','pools not implemented',npool)
  !
  IF ( ionode )  THEN
     !
     ! ... Input from file ?
     !
     nargs = iargc()
     !
     DO iiarg = 1, ( nargs - 1 )
        !
        CALL getarg( iiarg, input_file )
        IF ( TRIM( input_file ) == '-input' .OR. &
             TRIM( input_file ) == '-inp'   .OR. &
             TRIM( input_file ) == '-in' ) THEN
           !
           CALL getarg( ( iiarg + 1 ) , input_file )
           OPEN ( UNIT = 5, FILE = input_file, FORM = 'FORMATTED', &
                STATUS = 'OLD', IOSTAT = ierr )
           CALL errore( 'iosys', 'input file ' // TRIM( input_file ) // &
                & ' not found' , ierr )
           !
        END IF
        !
     END DO

     READ (5, inputpp, err = 200, iostat = ios)
200  CALL errore ('do_bands', 'reading inputpp namelist', ABS (ios) )
     !
     tmp_dir = TRIM(outdir)
     !
  END IF
  !
  ! ... Broadcast variables
  !
  CALL mp_bcast( tmp_dir, ionode_id )
  CALL mp_bcast( prefix, ionode_id )
  CALL mp_bcast( filband, ionode_id )
  !
  !   Now allocate space for pwscf variables, read and check them.
  !
  CALL read_file
  CALL openfil_pp
  CALL init_us_1
  !
  CALL punch_band (filband)
  !
  CALL stop_pp
  STOP
END PROGRAM bands
!
!-----------------------------------------------------------------------
SUBROUTINE punch_band (filband)
  !-----------------------------------------------------------------------
  !
  !    This routine writes the band energies on a file. The routine orders
  !    the eigenvalues using the overlap of the eigenvectors to give
  !    an estimate crossing and anticrossing of the bands. This simplified
  !    method works in many, but not in all the cases.
  !
  !
  USE atom
  USE ions_base,            ONLY : nat, ityp, ntyp => nsp
  USE cell_base
  USE constants,            ONLY : rytoev
  USE gvect
  USE lsda_mod,             ONLY : nspin
  USE klist
  USE io_files,             ONLY : iunpun, nwordwfc, iunwfc
  USE wvfct
  USE uspp,                 ONLY : nkb, vkb, qq
  USE uspp_param,           ONLY : tvanp, nh, nhm
  USE wavefunctions_module, ONLY : evc
  USE io_global,            ONLY : ionode 

  IMPLICIT NONE
  CHARACTER (len=*) :: filband
  REAL(kind=DP) :: proold
  ! the best overlap product
  COMPLEX(kind=DP) :: pro
  ! the product of wavefunctions

  COMPLEX(kind=DP), ALLOCATABLE :: psiold (:,:), old (:), NEW (:)
  ! psiold: eigenfunctions at previous k-point, ordered
  ! old, new: contain one band resp. at previous and current k-point
  COMPLEX(kind=DP), ALLOCATABLE :: becp(:,:), becpold (:,:)
  ! becp   : <psi|beta> at current  k-point
  ! becpold: <psi|beta> at previous k-point

  INTEGER :: ibnd, jbnd, ik, ikb, ig, npwold, ios
  ! counters
  INTEGER, ALLOCATABLE :: ok (:), igkold (:), il (:)
  ! ok: keeps track of which bands have been already ordered
  ! igkold: indices of k+G at previous k-point
  ! il: band ordering
  INTEGER, PARAMETER :: maxdeg = 4
  ! maxdeg : max allowed degeneracy
  INTEGER :: ndeg, deg, nd
  ! ndeg : number of degenerate states
  INTEGER, ALLOCATABLE :: degeneracy(:), degbands(:,:), INDEX(:)
  ! degbands keeps track of which states are degenerate
  REAL(kind=DP), ALLOCATABLE:: edeg(:)
  REAL(kind=DP), PARAMETER :: eps = 0.001
  ! threshold (Ry) for degenerate states 
  COMPLEX(kind=DP), EXTERNAL :: cgracsc
 ! scalar product with the S matrix

  IF (filband == ' ') RETURN
  IF (nspin == 2) CALL errore('bands','LSDA bands not implemented',-nspin)
  iunpun = 18
  !
  IF ( ionode ) THEN
     !
     OPEN (unit = iunpun, file = filband, status = 'unknown', form = &
          'formatted', err = 100, iostat = ios)
100  CALL errore ('punch_band', 'Opening filband file', ABS (ios) )
     REWIND (iunpun)
     !
  END IF
  !
  ALLOCATE (psiold( npwx, nbnd))    
  ALLOCATE (old(ngm), NEW(ngm))    
  ALLOCATE (becp(nkb, nbnd), becpold(nkb, nbnd))    
  ALLOCATE (igkold (npwx))    
  ALLOCATE (ok (nbnd), il (nbnd))    
  ALLOCATE (degeneracy(nbnd), edeg(nbnd))
  ALLOCATE (INDEX(maxdeg), degbands(nbnd,maxdeg))
  !
  DO ik = 1, nks
     !
     !    prepare the indices of this k point
     !
     CALL gk_sort (xk (1, ik), ngm, g, ecutwfc / tpiba2, npw, &
          igk, g2kin)
     !
     !   read eigenfunctions
     !
     CALL davcio (evc, nwordwfc, iunwfc, ik, - 1)
     !
     ! calculate becp = <psi|beta> 
     ! 
     CALL init_us_2 (npw, igk, xk (1, ik), vkb)
     CALL ccalbec (nkb, npwx, npw, nbnd, becp, vkb, evc)
     !
     IF (ik == 1) THEN
        !
        !  first k-point in the list:
        !  save eigenfunctions in the current order (increasing energy)
        !
        DO ibnd = 1, nbnd
           il (ibnd) = ibnd
        END DO
     ELSE
        !
        !  following  k-points in the list:
        !  determine eigenfunction order in array il
        !
        DO ibnd = 1, nbnd
           ok (ibnd) = 0
        ENDDO
        DO ibnd = 1, nbnd
           old(:) = (0.d0, 0.d0)
           DO ig = 1, npwold
              old (igkold (ig) ) = psiold (ig, ibnd)
           ENDDO
           proold = 0.d0
           DO jbnd = 1, nbnd
              IF (ok (jbnd) == 0) THEN
                 NEW (:) = (0.d0, 0.d0)
                 DO ig = 1, npw
                    NEW (igk (ig) ) = evc (ig, jbnd)
                 ENDDO
                 pro = cgracsc (nkb, becp (1, jbnd), becpold (1, ibnd), &
                      nhm, ntyp, nh, qq, nat, ityp, ngm, NEW, old, tvanp)
                 IF (ABS (pro) > proold) THEN
                    il (ibnd) = jbnd
                    proold = ABS (pro)
                 ENDIF
              ENDIF
           ENDDO
           ok (il (ibnd) ) = 1
        ENDDO
        !
        !  if there were bands crossing at degenerate eigenvalues
        !  at previous k-point, re-order those bands so as to keep
        !  lower band indices corresponding to lower bands
        !
        DO nd = 1, ndeg
           DO deg = 1, degeneracy (nd)
              INDEX(deg) = il(degbands(nd,deg))
              edeg (deg) = et(il(degbands(nd,deg)), ik)
           END DO
           CALL hpsort(degeneracy (nd), edeg, index)
           DO deg = 1, degeneracy (nd)
              il(degbands(nd,deg)) = INDEX(deg)
           END DO
        END DO
     END IF
     !
     !   Now the order of eigenfunctions has been established
     !   for this k-point -- prepare data for next k point
     !
     DO ibnd = 1, nbnd
        DO ig = 1, npw
           psiold (ig, ibnd) = evc (ig, il (ibnd) )
        ENDDO
        DO ikb = 1, nkb
           becpold (ikb, ibnd) = becp (ikb, il (ibnd) )
        ENDDO
     ENDDO
     DO ig = 1, npw
        igkold (ig) = igk (ig)
     ENDDO
     npwold = npw
     !
     !  find degenerate eigenvalues
     !
     deg  = 0
     ndeg = 0
     DO ibnd = 2, nbnd
        IF ( ABS (et(ibnd, ik) - et(ibnd-1, ik)) < eps ) THEN
           IF ( deg == 0 ) THEN
              ndeg = ndeg + 1
              edeg (ndeg) = et(ibnd, ik)
           END IF
           deg = 1
        ELSE
           deg = 0
        END IF
     END DO
     !
     !  locate band crossings at degenerate eigenvalues
     !
     DO nd = 1, ndeg
        deg = 0
        DO ibnd = 1, nbnd
           IF ( ABS (et(il(ibnd), ik) - edeg (nd)) < eps ) THEN
              deg = deg + 1
              IF (deg > maxdeg) CALL errore ('punch_band', &
                   ' increase maxdeg', deg)
              degbands(nd,deg) = ibnd
           END IF
        END DO
        degeneracy (nd) = deg
     END DO
     !
     IF ( ionode ) THEN
        !
        IF (ik == 1) THEN
           WRITE (iunpun, '(" &plot nbnd=",i4,", nks=",i4," /")') &
                nbnd, nks
        END IF
        WRITE (iunpun, '(14x,3f7.4)') xk(1,ik),xk(2,ik),xk(3,ik)
        WRITE (iunpun, '(10f8.3)') (et (il (ibnd) , ik) &
             * rytoev, ibnd = 1, nbnd)
        !
     END IF
     !
  ENDDO
  !
  DEALLOCATE (index, degbands)
  DEALLOCATE (edeg, degeneracy)
  DEALLOCATE (il, ok)
  DEALLOCATE (igkold)
  DEALLOCATE (becpold, becp)
  DEALLOCATE (NEW, old)
  DEALLOCATE (psiold)
  !
  IF ( ionode ) CLOSE (iunpun)
  !
  RETURN
  !
END SUBROUTINE punch_band
