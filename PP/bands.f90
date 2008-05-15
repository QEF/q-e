!
! Copyright (C) 2001-2007 Quantum-Espresso group
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
  USE mp_global, ONLY : npool, nproc, nproc_pool, nproc_file, nproc_pool_file
  USE control_flags, ONLY : twfcollect
  USE wvfct,     ONLY : nbnd
  USE klist,     ONLY : nkstot, two_fermi_energies
  USE noncollin_module, ONLY : i_cons
  USE io_global, ONLY : ionode, ionode_id, stdout
  USE mp,        ONLY : mp_bcast

  !
  IMPLICIT NONE
  !
  CHARACTER (len=256) :: filband, filp, outdir
  LOGICAL :: lsigma(4), lsym, lp, no_overlap
  INTEGER :: spin_component, firstk, lastk
  INTEGER :: ios
  !
  NAMELIST / inputpp / outdir, prefix, filband, filp, spin_component, lsigma,&
                       lsym, lp, filp, firstk, lastk, no_overlap
  !                                  
  !
  CALL start_postproc (nd_nmbr)
  !
  !   set default values for variables in namelist
  !
  prefix = 'pwscf'
  CALL get_env( 'ESPRESSO_TMPDIR', outdir )
  IF ( TRIM( outdir ) == ' ' ) outdir = './'
  filband = 'bands.out'
  lsym=.false.
  lsigma=.false.
  filp='p_avg.dat'
  lp=.false.
  firstk=0
  lastk=10000000
  spin_component = 1
  no_overlap=.false.
  !
  ios = 0
  !
  IF ( ionode )  THEN
     !
     CALL input_from_file ( )
     !
     READ (5, inputpp, iostat = ios)
     !
     lsigma(4)=.false.
     tmp_dir = trimcheck (outdir)
     !
  END IF
  !
  !
  CALL mp_bcast( ios, ionode_id )
  IF (ios /= 0) CALL errore ('do_bands', 'reading inputpp namelist', ABS(ios) )
  !
  ! ... Broadcast variables
  !
  CALL mp_bcast( tmp_dir, ionode_id )
  CALL mp_bcast( prefix, ionode_id )
  CALL mp_bcast( filband, ionode_id )
  CALL mp_bcast( filp, ionode_id )
  CALL mp_bcast( spin_component, ionode_id )
  CALL mp_bcast( firstk, ionode_id )
  CALL mp_bcast( lastk, ionode_id )
  CALL mp_bcast( lp, ionode_id )
  CALL mp_bcast( lsym, ionode_id )
  CALL mp_bcast( lsigma, ionode_id )
  CALL mp_bcast( no_overlap, ionode_id )

  IF ( npool > 1 .and..not.lsym) CALL errore('bands', &
                                             'pools not implemented',npool)
  !
  !   Now allocate space for pwscf variables, read and check them.
  !
  CALL read_file()

  IF (nproc /= nproc_file .and. .not. twfcollect)  &
     CALL errore('bands',&
     'pw.x run with a different number of processors. Use wf_collect=.true.',1)

  IF (nproc_pool /= nproc_pool_file .and. .not. twfcollect)  &
     CALL errore('bands',&
     'pw.x run with a different number of pools. Use wf_collect=.true.',1)

  IF (two_fermi_energies.or.i_cons /= 0) &
     CALL errore('bands',&
     'The bands code with constrained magnetization has not been tested',1)

  CALL openfil_pp()
  !
  IF (lsym) no_overlap=.true.
  CALL punch_band(filband,spin_component,lsigma,no_overlap)
  IF (lsym) call sym_band(filband,spin_component,firstk,lastk)
  IF (lp) CALL write_p_avg(filp,spin_component,firstk,lastk)
  !
  CALL stop_pp
  STOP
END PROGRAM bands
!
!-----------------------------------------------------------------------
SUBROUTINE punch_band (filband, spin_component, lsigma, no_overlap)
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
  USE uspp_param,           ONLY : upf, nh, nhm
  USE noncollin_module,     ONLY : noncolin, npol
  USE wavefunctions_module, ONLY : evc
  USE io_global,            ONLY : ionode, ionode_id
  USE mp,                   ONLY : mp_bcast
  USE becmod,               ONLY : calbec

  IMPLICIT NONE
  CHARACTER (len=*) :: filband
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
  COMPLEX(DP), ALLOCATABLE :: psiold_nc (:,:), old_nc(:,:), new_nc(:,:)
  COMPLEX(DP), ALLOCATABLE :: becp_nc(:,:,:), becpold_nc(:,:,:)
  LOGICAL :: no_overlap
  ! as above for the noncolinear case
  INTEGER :: ibnd, jbnd, ik, ikb, ig, npwold, nks1, nks2, ipol, ih, is1
  INTEGER :: nks1tot, nks2tot
  ! counters
  INTEGER, ALLOCATABLE :: ok (:), igkold (:), il (:), ilold(:)
  ! ok: keeps track of which bands have been already ordered
  ! igkold: indices of k+G at previous k-point
  ! il: band ordering
  INTEGER :: maxdeg 
  ! maxdeg : max allowed degeneracy
  INTEGER :: ndeg, deg, nd
  ! ndeg : number of degenerate states
  INTEGER, ALLOCATABLE :: degeneracy(:), degbands(:,:), idx(:)
  ! degbands keeps track of which states are degenerate
  INTEGER :: iunpun_sigma(4), ios(0:4), indjbnd
  CHARACTER(LEN=256) :: nomefile
  REAL(DP), ALLOCATABLE:: edeg(:)
  REAL(DP), ALLOCATABLE:: sigma_avg(:,:,:)
  ! expectation value of sigma
  REAL(DP), PARAMETER :: eps = 0.00001d0
  ! threshold (Ry) for degenerate states 
  REAL(DP) :: minene
  COMPLEX(DP), EXTERNAL :: cgracsc, cgracsc_nc
 ! scalar product with the S matrix

  IF (filband == ' ') RETURN
  DO ipol=1,4
     IF (lsigma(ipol).and..not.noncolin) THEN
        CALL errore ('punch_band', 'lsigma requires noncollinear run', &
                    ipol )
        lsigma=.false.
     ENDIF
  ENDDO
  
  iunpun = 18
  maxdeg = 30 * npol 
  ! 
  ios(:) = 0
  IF ( ionode ) THEN
     !
     OPEN (unit = iunpun, file = filband, status = 'unknown', form = &
          'formatted', iostat = ios(0))
     REWIND (iunpun)
     DO ipol=1,4
        IF (lsigma(ipol)) THEN
           iunpun_sigma(ipol)=iunpun+ipol
           WRITE(nomefile,'(".",i1)') ipol
           OPEN (unit = iunpun_sigma(ipol),  &
                 file = TRIM(filband)//TRIM(nomefile), &
                 status = 'unknown', form='formatted', iostat = ios(ipol))
           REWIND (iunpun_sigma(ipol))
        ENDIF
     ENDDO
     !
  END IF
  !
  CALL mp_bcast( ios, ionode_id )
  IF ( ios(0) /= 0 ) &
     CALL errore ('punch_band', 'Opening filband file', ABS(ios(0)) )
  DO ipol=1,4
     IF ( ios(ipol) /= 0 ) &
        CALL errore ('punch_band', 'Opening filband.N file ', ipol)
  END DO
  !
  IF (noncolin) THEN
     ALLOCATE (psiold_nc( npwx*npol, nbnd))
     ALLOCATE (becp_nc(nkb, npol, nbnd), becpold_nc(nkb, npol, nbnd))
     ALLOCATE (old_nc(ngm,npol), new_nc(ngm,npol))
     ALLOCATE (sigma_avg(4,nbnd,nkstot))
  ELSE
     ALLOCATE (psiold( npwx, nbnd))    
     ALLOCATE (old(ngm), new(ngm))    
     ALLOCATE (becp(nkb, nbnd), becpold(nkb, nbnd))    
  END IF

  ALLOCATE (igkold (npwx))    
  ALLOCATE (ok (nbnd), il (nbnd), ilold(nbnd) )
  ALLOCATE (degeneracy(nbnd), edeg(nbnd))
  ALLOCATE (idx(nbnd), degbands(nbnd,maxdeg))
  !
  IF (spin_component.NE.1.AND.nspin.NE.2) &
     CALL errore('punch_bands','uncorrect spin_component',1)
  IF (spin_component<1.OR.spin_component>2) &
     CALL errore('punch_bands','uncorrect lsda spin_component',1)

  CALL find_nks1nks2(1,nkstot,nks1tot,nks1,nks2tot,nks2,spin_component)

  DO ibnd = 1, nbnd
     il (ibnd) = ibnd
  END DO

  DO ik = nks1, nks2
     !
     !    prepare the indices of this k point
     !
     IF (.NOT.no_overlap.OR.lsigma(1).OR.lsigma(2).OR.lsigma(3).OR.lsigma(4)) THEN
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
        IF (noncolin) THEN
           CALL calbec ( npw, vkb, evc, becp_nc )
           CALL compute_sigma_avg(sigma_avg(1,1,ik),becp_nc,ik,lsigma)
        ELSE
           CALL calbec ( npw, vkb, evc, becp )
        END IF
     END IF
     !
     IF (ik==nks1.or.no_overlap) THEN
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
!
! The bands are checked in order of increasing energy.
!
        DO ibnd=1,nbnd
           idx(ibnd)=ibnd
           edeg(ibnd)=et(il(ibnd),ik-1)
        ENDDO
        CALL hpsort(nbnd, edeg, idx)
        DO ibnd = 1, nbnd
           IF (noncolin) THEN
              old_nc = (0.d0, 0.d0)
              DO ig = 1, npwold
                 old_nc(igkold(ig), 1)=psiold_nc(ig     ,idx(ibnd))
                 old_nc(igkold(ig), 2)=psiold_nc(ig+npwx,idx(ibnd))
              END DO
           ELSE
              old = (0.d0, 0.d0)
              DO ig = 1, npwold
                 old (igkold (ig) ) = psiold (ig, idx(ibnd))
              END DO
           END IF
           DO jbnd = 1, nbnd
              IF (ok (jbnd) == 0) THEN
                 IF (noncolin) THEN
                    new_nc = (0.d0, 0.d0)
                    DO ig = 1, npw
                       new_nc (igk (ig), 1) = evc (ig     , jbnd)
                       new_nc (igk (ig), 2) = evc (ig+npwx, jbnd)
                    END DO
                    pro = cgracsc_nc (nkb,becp_nc(1,1,jbnd), &
                              becpold_nc(1,1,idx(ibnd)), nhm, ntyp, nh, &
                              nat, ityp, ngm, npol, new_nc, old_nc, upf)
                 ELSE
                    new (:) = (0.d0, 0.d0)
                    DO ig = 1, npw
                       new (igk (ig) ) = evc (ig, jbnd)
                    END DO
                    pro=cgracsc(nkb,becp(1,jbnd),becpold(1,idx(ibnd)), &
                         nhm, ntyp, nh, qq, nat, ityp, ngm, NEW, old, upf)
                 END IF
!                 write(6,'(3i5,f15.10)') ik,idx(ibnd), jbnd, abs(pro)
                 IF (abs (pro) > 1.d-2 ) THEN
                    il (idx(ibnd)) = jbnd
                    GOTO 10
                 END IF
              END IF
           END DO
!           WRITE(6,*) '  no band found', ik, ilold(idx(ibnd)), &
!                        et(ilold(idx(ibnd)),ik-1)*rytoev
!
!     no band found. Takes the closest in energy. NB: This should happen only
!     for high energy bands. 
!
           minene=1.d10
           DO jbnd = 1, nbnd
              IF (ok (jbnd) == 0) THEN
                 IF (abs(et(idx(ibnd),ik)-et(jbnd,ik))<minene) THEN
                    indjbnd=jbnd
                    minene=abs(et(idx(ibnd),ik)-et(jbnd,ik))
                 ENDIF
              ENDIF
           ENDDO
           il(idx(ibnd))=indjbnd
10         CONTINUE
           ok (il (idx(ibnd)) ) = 1
        END DO
        !
        !  if there were bands crossing at degenerate eigenvalues
        !  at previous k-point, re-order those bands so as to keep
        !  lower band indices corresponding to lower bands
        !
        DO nd = 1, ndeg
           DO deg = 1, degeneracy (nd)
              idx(deg) = il(degbands(nd,deg))
              edeg (deg) = et(il(degbands(nd,deg)), ik)
           END DO
           CALL hpsort(degeneracy (nd), edeg, idx)
           DO deg = 1, degeneracy (nd)
              il(degbands(nd,deg)) = idx(deg)
           END DO
        END DO
     END IF
     !
     !   Now the order of eigenfunctions has been established
     !   for this k-point -- prepare data for next k point
     !
     IF (.NOT.no_overlap.OR.lsigma(1).OR.lsigma(2).OR.lsigma(3).OR.lsigma(4)) THEN
        DO ibnd = 1, nbnd
           IF (noncolin) THEN
              psiold_nc(:, ibnd) = evc(:, il (ibnd))
              DO ipol=1,npol
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
        ilold=il
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
     END IF
  END DO
#ifdef __PARA
  IF (noncolin) CALL poolrecover(sigma_avg,4*nbnd,nkstot,nks)
#endif
  !
  IF ( ionode ) THEN
     !
     DO ik=nks1tot,nks2tot
        IF (ik == nks1) THEN
           WRITE (iunpun, '(" &plot nbnd=",i4,", nks=",i4," /")') &
             nbnd, nks2tot-nks1tot+1
           DO ipol=1,4
              IF (lsigma(ipol)) WRITE(iunpun_sigma(ipol), &
                            '(" &plot nbnd=",i4,", nks=",i4," /")') &
                             nbnd, nks2tot-nks1tot+1
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
     END DO
  END IF
  !
  DEALLOCATE (idx, degbands)
  DEALLOCATE (edeg, degeneracy)
  DEALLOCATE (ilold, il, ok)
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
