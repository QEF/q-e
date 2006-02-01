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
  USE io_files,  ONLY : nd_nmbr, prefix, tmp_dir, trimcheck
  USE mp_global, ONLY : npool
  USE io_global, ONLY : ionode, ionode_id
  USE mp,        ONLY : mp_bcast
  !
  IMPLICIT NONE
  !
  CHARACTER (len=256) :: filband, outdir
  LOGICAL :: lsigma(4)
  INTEGER :: spin_component
  INTEGER :: ios
  !
  NAMELIST / inputpp / outdir, prefix, filband, spin_component, lsigma
  !                                  
  !
  CALL start_postproc (nd_nmbr)
  !
  !   set default values for variables in namelist
  !
  prefix = 'pwscf'
  outdir = './'
  filband = 'bands.out'
  lsigma = .false.
  spin_component = 1
  !
  IF ( npool > 1 ) CALL errore('bands','pools not implemented',npool)
  !
  IF ( ionode )  THEN
     !
     CALL input_from_file ( )
     !
     READ (5, inputpp, err = 200, iostat = ios)
200  CALL errore ('do_bands', 'reading inputpp namelist', ABS (ios) )
     !
     lsigma(4)=.false.
     tmp_dir = trimcheck (outdir)
     !
  END IF
  !
  ! ... Broadcast variables
  !
  CALL mp_bcast( tmp_dir, ionode_id )
  CALL mp_bcast( prefix, ionode_id )
  CALL mp_bcast( filband, ionode_id )
  CALL mp_bcast( spin_component, ionode_id )
  CALL mp_bcast( lsigma(:), ionode_id )
  !
  !   Now allocate space for pwscf variables, read and check them.
  !
  CALL read_file
  CALL openfil_pp
  CALL init_us_1
  !
  CALL punch_band (filband, spin_component, lsigma)
  !
  CALL stop_pp
  STOP
END PROGRAM bands
!
!-----------------------------------------------------------------------
SUBROUTINE punch_band (filband, spin_component, lsigma)
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
  USE noncollin_module,     ONLY : noncolin, npol
  USE wavefunctions_module, ONLY : evc, evc_nc
  USE io_global,            ONLY : ionode 

  IMPLICIT NONE
  CHARACTER (len=*) :: filband
  REAL(DP) :: proold
  ! the best overlap product
  COMPLEX(DP) :: pro
  ! the product of wavefunctions
  INTEGER :: spin_component
  LOGICAL :: lsigma(4)

  COMPLEX(DP), ALLOCATABLE :: psiold (:,:), old (:), new (:)
  ! psiold: eigenfunctions at previous k-point, ordered
  ! old, new: contain one band resp. at previous and current k-point
  COMPLEX(DP), ALLOCATABLE :: becp(:,:), becpold (:,:)
  ! becp   : <psi|beta> at current  k-point
  ! becpold: <psi|beta> at previous k-point
  COMPLEX(DP), ALLOCATABLE :: psiold_nc (:,:,:), old_nc(:,:), new_nc(:,:)
  COMPLEX(DP), ALLOCATABLE :: becp_nc(:,:,:), becpold_nc(:,:,:)
  ! as above for the noncolinear case
  INTEGER :: ibnd, jbnd, ik, ikb, ig, npwold, ios, nks1, nks2, ipol, ih, is1
  ! counters
  INTEGER, ALLOCATABLE :: ok (:), igkold (:), il (:)
  ! ok: keeps track of which bands have been already ordered
  ! igkold: indices of k+G at previous k-point
  ! il: band ordering
  INTEGER :: maxdeg 
  ! maxdeg : max allowed degeneracy
  INTEGER :: ndeg, deg, nd
  ! ndeg : number of degenerate states
  INTEGER, ALLOCATABLE :: degeneracy(:), degbands(:,:), INDEX(:)
  ! degbands keeps track of which states are degenerate
  INTEGER :: iunpun_sigma(4)
  CHARACTER(LEN=256) :: nomefile
  REAL(DP), ALLOCATABLE:: edeg(:)
  REAL(DP), ALLOCATABLE:: sigma_avg(:,:,:)
  ! expectation value of sigma
  REAL(DP), PARAMETER :: eps = 0.001
  ! threshold (Ry) for degenerate states 
  COMPLEX(DP), EXTERNAL :: cgracsc, cgracsc_nc
 ! scalar product with the S matrix

  IF (filband == ' ') RETURN
  DO ipol=1,4
     IF (lsigma(ipol).and..not.noncolin) THEN
        CALL errore ('punch_band', 'lsigma requires noncollinear run', &
                    ABS (ios) )
        lsigma=.false.
     ENDIF
  ENDDO
  
  iunpun = 18
  maxdeg = 4 * npol 
  !
  IF ( ionode ) THEN
     !
     OPEN (unit = iunpun, file = filband, status = 'unknown', form = &
          'formatted', err = 100, iostat = ios)
100  CALL errore ('punch_band', 'Opening filband file', ABS (ios) )
     REWIND (iunpun)
     DO ipol=1,4
        IF (lsigma(ipol)) THEN
           iunpun_sigma(ipol)=iunpun+ipol
           WRITE(nomefile,'(".",i1)') ipol
           OPEN (unit = iunpun_sigma(ipol),  &
                 file = TRIM(filband)//TRIM(nomefile), &
                 status = 'unknown', form='formatted', err = 200, iostat = ios)
200        CALL errore ('punch_band', 'Opening filband.1 file', ABS (ios) )
           REWIND (iunpun_sigma(ipol))
        ENDIF
     ENDDO
     !
  END IF
  !
  IF (noncolin) THEN
     ALLOCATE (psiold_nc( npwx, npol, nbnd))
     ALLOCATE (becp_nc(nkb, npol, nbnd), becpold_nc(nkb, npol, nbnd))
     ALLOCATE (old_nc(ngm,npol), new_nc(ngm,npol))
     ALLOCATE (sigma_avg(4,nbnd,nks))
  ELSE
     ALLOCATE (psiold( npwx, nbnd))    
     ALLOCATE (old(ngm), new(ngm))    
     ALLOCATE (becp(nkb, nbnd), becpold(nkb, nbnd))    
  END IF

  ALLOCATE (igkold (npwx))    
  ALLOCATE (ok (nbnd), il (nbnd))    
  ALLOCATE (degeneracy(nbnd), edeg(nbnd))
  ALLOCATE (INDEX(maxdeg), degbands(nbnd,maxdeg))
  !
  IF (nspin==1.OR.nspin==4) THEN
     nks1=1
     nks2=nks
     if (spin_component .ne. 1)  &
        CALL errore('punch_bands','uncorrect spin_component',1)
  ELSE IF (nspin.eq.2) THEN
     IF (spin_component == 1) THEN
        nks1=1
        nks2=nks/2
     ELSE IF (spin_component==2) THEN
        nks1=nks/2 + 1
        nks2=nks
     ELSE
        CALL errore('punch_bands','uncorrect spin_component',1)
     END IF
  ELSE
     CALL errore('punch_bands','nspin not allowed in bands', 1)
  END IF
  DO ik = nks1, nks2
     !
     !    prepare the indices of this k point
     !
     CALL gk_sort (xk (1, ik), ngm, g, ecutwfc / tpiba2, npw, &
          igk, g2kin)
     !
     !   read eigenfunctions
     !
     IF (noncolin) THEN
        CALL davcio (evc_nc, nwordwfc, iunwfc, ik, - 1)
     ELSE
        CALL davcio (evc, nwordwfc, iunwfc, ik, - 1)
     END IF
     !
     ! calculate becp = <psi|beta> 
     ! 
     CALL init_us_2 (npw, igk, xk (1, ik), vkb)
     IF (noncolin) THEN
        CALL ccalbec_nc (nkb, npwx, npw, npol, nbnd, becp_nc, vkb, evc_nc)
        CALL compute_sigma_avg(sigma_avg(1,1,ik),becp_nc,ik,lsigma)
     ELSE
        CALL ccalbec (nkb, npwx, npw, nbnd, becp, vkb, evc)
     END IF
     !
     IF (ik == nks1) THEN
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
           IF (noncolin) THEN
              old_nc = (0.d0, 0.d0)
              DO ipol=1, npol
                 DO ig = 1, npwold
                    old_nc(igkold(ig),ipol)=psiold_nc(ig,ipol,ibnd)
                 END DO
              END DO
              proold = 0.d0
              DO jbnd = 1, nbnd
                 IF (ok (jbnd) == 0) THEN
                    new_nc = (0.d0, 0.d0)
                    DO ipol=1,npol
                       DO ig = 1, npw
                          new_nc (igk (ig),ipol ) = evc_nc (ig,ipol, jbnd)
                       END DO
                    END DO
                    pro = cgracsc_nc (nkb,becp_nc(1,1,jbnd), &
                              becpold_nc(1,1,ibnd), nhm, ntyp, nh, &
                              nat, ityp, ngm, npol, new_nc, old_nc, tvanp)
                    IF (abs (pro) > proold) THEN
                       il (ibnd) = jbnd
                       proold = abs (pro)
                    END IF
                 END IF
              END DO
              ok (il (ibnd) ) = 1
           ELSE
              old = (0.d0, 0.d0)
              DO ig = 1, npwold
                 old (igkold (ig) ) = psiold (ig, ibnd)
              END DO
              proold = 0.d0
              DO jbnd = 1, nbnd
                 IF (ok (jbnd) == 0) THEN
                    new (:) = (0.d0, 0.d0)
                    DO ig = 1, npw
                       new (igk (ig) ) = evc (ig, jbnd)
                    END DO
                    pro = cgracsc (nkb, becp (1, jbnd), becpold (1, ibnd), &
                         nhm, ntyp, nh, qq, nat, ityp, ngm, NEW, old, tvanp)
                    IF (ABS (pro) > proold) THEN
                       il (ibnd) = jbnd
                       proold = ABS (pro)
                    END IF
                 END IF
              END DO
              ok (il (ibnd) ) = 1
           END IF
        END DO
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
        IF (noncolin) THEN
           DO ipol=1,npol
              DO ig = 1, npw
                 psiold_nc(ig, ipol, ibnd) = evc_nc(ig, ipol, il (ibnd))
              END DO
              DO ikb = 1, nkb
                 becpold_nc(ikb, ipol, ibnd) = becp_nc(ikb, ipol, il(ibnd))
              END DO
           END DO
        ELSE
           DO ig = 1, npw
              psiold (ig, ibnd) = evc (ig, il (ibnd) )
           END DO
           DO ikb = 1, nkb
              becpold (ikb, ibnd) = becp (ikb, il (ibnd) )
           END DO
        END IF
     END DO
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
        IF (ik == nks1) THEN
           WRITE (iunpun, '(" &plot nbnd=",i4,", nks=",i4," /")') &
                nbnd, nks2-nks1+1
           DO ipol=1,4
              IF (lsigma(ipol)) WRITE(iunpun_sigma(ipol), &
                            '(" &plot nbnd=",i4,", nks=",i4," /")') &
                             nbnd, nks2-nks1+1
           END DO
        END IF
        WRITE (iunpun, '(10x,3f10.6)') xk(1,ik),xk(2,ik),xk(3,ik)
        WRITE (iunpun, '(10f8.3)') (et (il (ibnd) , ik)             &
             * rytoev, ibnd = 1, nbnd)
        DO ipol=1,4
           IF (lsigma(ipol)) THEN
              WRITE (iunpun_sigma(ipol), '(10x,3f10.6)')            &
                                          xk(1,ik),xk(2,ik),xk(3,ik)
              WRITE (iunpun_sigma(ipol), '(10f8.3)')                &
                            (sigma_avg(ipol, il (ibnd) , ik), ibnd = 1, nbnd)
           END IF
        END DO
        !
     END IF
     !
  END DO
  !
  DEALLOCATE (index, degbands)
  DEALLOCATE (edeg, degeneracy)
  DEALLOCATE (il, ok)
  DEALLOCATE (igkold)
  IF (noncolin) THEN
     DEALLOCATE (sigma_avg)
     DEALLOCATE (becpold_nc, becp_nc)
     DEALLOCATE (new_nc, old_nc)
     DEALLOCATE (psiold_nc)
  ELSE
     DEALLOCATE (becpold, becp)
     DEALLOCATE (new, old)
     DEALLOCATE (psiold)
  END IF
  !
  IF ( ionode ) THEN
     CLOSE (iunpun)
     DO ipol=1,4
        IF (lsigma(ipol)) CLOSE(iunpun_sigma(ipol))
     ENDDO
  ENDIF
  !
  RETURN
  !
END SUBROUTINE punch_band
