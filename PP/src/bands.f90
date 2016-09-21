!
! Copyright (C) 2001-2016 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
PROGRAM do_bands
  !-----------------------------------------------------------------------
  !
  ! See files INPUT_BANDS.* in Doc/ directory for usage
  ! 
  !
  USE io_files,  ONLY : prefix, tmp_dir
  USE mp_global, ONLY : npool, nproc_pool, nproc_file, &
                        nproc_pool_file, mp_startup
  USE control_flags, ONLY : twfcollect, gamma_only
  USE environment,   ONLY : environment_start, environment_end
  USE wvfct,     ONLY : nbnd
  USE klist,     ONLY : nkstot, two_fermi_energies
  USE noncollin_module, ONLY : noncolin, i_cons
  USE lsda_mod,  ONLY : nspin
  USE io_global, ONLY : ionode, ionode_id, stdout
  USE mp,        ONLY : mp_bcast
  USE mp_world,  ONLY : world_comm
  !
  IMPLICIT NONE
  !
  CHARACTER(LEN=256), EXTERNAL :: trimcheck
  !
  CHARACTER (len=256) :: filband, filp, outdir
  LOGICAL :: lsigma(4), lsym, lp, no_overlap, plot_2d
  INTEGER :: spin_component, firstk, lastk
  INTEGER :: ios
  !
  NAMELIST / bands / outdir, prefix, filband, filp, spin_component, lsigma,&
                       lsym, lp, filp, firstk, lastk, no_overlap, plot_2d
  !
  ! initialise environment
  !
#if defined(__MPI)
  CALL mp_startup ( )
#endif
  CALL environment_start ( 'BANDS' )
  !
  !   set default values for variables in namelist
  !
  prefix = 'pwscf'
  CALL get_environment_variable( 'ESPRESSO_TMPDIR', outdir )
  IF ( trim( outdir ) == ' ' ) outdir = './'
  filband = 'bands.out'
  lsym=.true.
  no_overlap=.true.
  plot_2d=.false.
  lsigma=.false.
  lp=.false.
  filp='p_avg.dat'
  firstk=0
  lastk=10000000
  spin_component = 1
  !
  ios = 0
  !
  IF ( ionode )  THEN
     !
     CALL input_from_file ( )
     !
     READ (5, bands, iostat = ios)
     !
     lsigma(4)=.false.
     tmp_dir = trimcheck (outdir)
     !
  ENDIF
  !
  !
  CALL mp_bcast( ios, ionode_id, world_comm )
  IF (ios /= 0) CALL errore ('bands', 'reading bands namelist', abs(ios) )
  !
  ! ... Broadcast variables
  !
  CALL mp_bcast( tmp_dir, ionode_id, world_comm )
  CALL mp_bcast( prefix, ionode_id, world_comm )
  CALL mp_bcast( filband, ionode_id, world_comm )
  CALL mp_bcast( filp, ionode_id, world_comm )
  CALL mp_bcast( spin_component, ionode_id, world_comm )
  CALL mp_bcast( firstk, ionode_id, world_comm )
  CALL mp_bcast( lastk, ionode_id, world_comm )
  CALL mp_bcast( lp, ionode_id, world_comm )
  CALL mp_bcast( lsym, ionode_id, world_comm )
  CALL mp_bcast( lsigma, ionode_id, world_comm )
  CALL mp_bcast( no_overlap, ionode_id, world_comm )
  CALL mp_bcast( plot_2d, ionode_id, world_comm )

  IF (plot_2d) THEN
     lsym=.false.
     lp=.false.
     no_overlap=.true.
  ENDIF
  IF (lsym) no_overlap=.true.

  IF ( npool > 1 .and..not.(lsym.or.no_overlap)) CALL errore('bands', &
                                             'pools not implemented',npool)
  IF ( spin_component < 1 .OR. spin_component > 2 ) &
     CALL errore('bands','incorrect spin_component',1)
  !
  !   Now allocate space for pwscf variables, read and check them.
  !
  CALL read_file()
  !
  IF (gamma_only) CALL errore('bands','gamma_only case not implemented',1)
  IF (nproc_pool /= nproc_pool_file .and. .not. twfcollect)  &
     CALL errore('bands',&
     'pw.x run with a different number of procs/pools. Use wf_collect=.true.',1)
  IF (two_fermi_energies.or.i_cons /= 0) &
     CALL errore('bands',&
     'The bands code with constrained magnetization has not been tested',1)
  IF ( ANY(lsigma) .AND. .NOT.noncolin ) &
     CALL errore ('punch_band', 'lsigma requires noncollinear run', 1 )
  IF ( spin_component/=1 .and. nspin/=2 ) &
     CALL errore('punch_bands','incorrect spin_component',1)
  !
  CALL openfil_pp()
  !
  IF (plot_2d) THEN
     CALL punch_band_2d(filband,spin_component)
  ELSE
     CALL punch_band(filband,spin_component,lsigma,no_overlap)
     IF (lsym) CALL sym_band(filband,spin_component,firstk,lastk)
     IF (lp) CALL write_p_avg(filp,spin_component,firstk,lastk)
  END IF
  !
  CALL environment_end ( 'BANDS' )
  !
  CALL stop_pp
  STOP
END PROGRAM do_bands
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
  USE kinds,                ONLY : dp
  USE ions_base,            ONLY : nat, ityp, ntyp => nsp
  USE cell_base,            ONLY : at
  USE constants,            ONLY : rytoev
  USE gvect,                ONLY : g, ngm
  USE klist,                ONLY : xk, nks, nkstot, ngk, igk_k
  USE io_files,             ONLY : iunpun, nwordwfc, iunwfc
  USE wvfct,                ONLY : nbnd, et, npwx
  USE uspp,                 ONLY : nkb, vkb, qq
  USE uspp_param,           ONLY : upf, nh, nhm
  USE noncollin_module,     ONLY : noncolin, npol
  USE wavefunctions_module, ONLY : evc
  USE io_global,            ONLY : ionode, ionode_id, stdout
  USE mp,                   ONLY : mp_bcast
  USE mp_world,             ONLY : world_comm
  USE becmod,               ONLY : calbec, bec_type, allocate_bec_type, &
                                   deallocate_bec_type

  IMPLICIT NONE
  CHARACTER (len=*) :: filband
  INTEGER, INTENT(IN) :: spin_component
  LOGICAL, INTENT(IN) :: lsigma(4), no_overlap

  TYPE(bec_type):: becp
  ! becp   : <psi|beta> at current  k-point
  INTEGER :: ibnd, jbnd, i, ik, ig, ig1, ig2, ipol, npw, ngmax, jmax
  INTEGER :: nks1tot, nks2tot, nks1, nks2
  INTEGER :: iunpun_sigma(4), ios(0:4), done(nbnd)
  CHARACTER(len=256) :: nomefile
  REAL(dp):: pscur, psmax, psr(nbnd)
  COMPLEX(dp), ALLOCATABLE :: psi(:,:), spsi(:,:), ps(:,:)
  INTEGER, ALLOCATABLE :: work(:), igg(:)
  INTEGER, ALLOCATABLE :: closest_band(:,:)! index for band ordering
  REAL(DP), ALLOCATABLE:: sigma_avg(:,:,:) ! expectation value of sigma
  REAL(DP), ALLOCATABLE:: et_(:,:) ! reordered eigenvalues in eV

  
  IF (filband == ' ') RETURN

  iunpun = 19
  ios(:) = 0
  !
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
                 file = trim(filband)//trim(nomefile), &
                 status = 'unknown', form='formatted', iostat = ios(ipol))
           REWIND (iunpun_sigma(ipol))
        ENDIF
     ENDDO
     !
  ENDIF
  !
  CALL mp_bcast( ios, ionode_id, world_comm )
  IF ( ios(0) /= 0 ) &
     CALL errore ('punch_band', 'Opening filband file', abs(ios(0)) )
  DO ipol=1,4
     IF ( ios(ipol) /= 0 ) &
        CALL errore ('punch_band', 'Opening filband.N file ', ipol)
  ENDDO
  !
  CALL find_nks1nks2(1,nkstot,nks1tot,nks1,nks2tot,nks2,spin_component)
  !
  ! index of largest G in plane-wave basis set across k-points
  !
  ALLOCATE ( closest_band(nbnd,nkstot)  )
  DO ik=1,nkstot
     !
     ! default ordering: leave bands as they are
     !
     DO ibnd=1,nbnd
        closest_band(ibnd,ik) = ibnd
     END DO
  END DO
  !
  IF ( noncolin ) ALLOCATE ( sigma_avg(4,nbnd,nkstot) )
  ALLOCATE ( psi(npwx*npol,nbnd), spsi(npwx*npol,nbnd), ps(nbnd,nbnd) )
  CALL allocate_bec_type(nkb, nbnd, becp)
  !
  DO ik = nks1, nks2
     !
     !   read eigenfunctions
     !
     CALL davcio (evc, 2*nwordwfc, iunwfc, ik, - 1)
     !
     ! calculate becp = <psi|beta>, needed to compute spsi = S|psi>
     !
     npw = ngk(ik)
     CALL init_us_2 (npw, igk_k(1,ik), xk(1,ik), vkb)
     CALL calbec ( npw, vkb, evc, becp )
     !
     ! calculate average magnetization for noncolinear case
     !
     IF (noncolin) &
           CALL compute_sigma_avg(sigma_avg(1,1,ik),becp%nc,ik,lsigma)
     !
     IF ( ik > nks1 .AND. .NOT. no_overlap ) THEN
        !
        ! compute correspondence between k+G_i indices at current and previous k
        !
        ngmax = MAXVAL ( igk_k(:,ik-1:ik) )
        ALLOCATE ( work(ngmax), igg(ngmax) )
        work(:) = 0
        DO ig1 = 1, ngk(ik-1)
           work( igk_k(ig1,ik-1)) = ig1
        END DO
        igg(:) = 0
        DO ig2 = 1, npw
           ig1 = work( igk_k(ig2,ik)) 
           IF (ig1 > 0) igg(ig2) = ig1
        END DO
        !
        ! compute overlap <\psi_k|S\psi_{k-1}> (without the Bloch factor)
        ! psi(G) = \psi(k+G) (order of G-vectors same as for previous k)
        !
        psi(:,:) = (0.0_dp,0.0_dp)
        DO ig2 = 1, npw
           IF ( igg(ig2) > 0 ) psi(igg(ig2),:) = evc(ig2,:)
        END DO
        IF ( noncolin) THEN
           DO ig2 = 1, npw
              IF ( igg(ig2) > 0 ) psi(npwx+igg(ig2),:) = evc(npwx+ig2,:)
           END DO
        END IF
        DEALLOCATE (igg, work)
        CALL calbec (ngk(ik-1), spsi, psi, ps )
        !
        ! ps(ibnd,jbnd) = <S\psi_{k-1,ibnd} | \psi_{k,jbnd}>
        !
        ! assign bands on the basis of the relative overlap
        ! simple cases first: large or very small overlap
        !
        closest_band(:,ik) = -1
        done(:) =0
        !ndone = 0
        !nlost = 0
        DO ibnd=1,nbnd
           !
           psr(:) = real(ps(ibnd,:))**2+aimag(ps(ibnd,:))**2
           psmax = MAXVAL( psr )
           !
           IF ( psmax > 0.75 ) THEN
              ! simple case: large overlap with one specific band
              closest_band(ibnd,ik) = MAXLOC( psr, 1 )
              ! record that this band at ik has been linked to a band at ik-1 
              done( closest_band(ibnd,ik) ) = 1
              ! ndone = ndone + 1
              !
        !   ELSE IF ( psmax < 0.05 ) THEN
        !      ! simple case: negligible overlap with all bands
        !      closest_band(ibnd,ik) = 0
        !      nlost = nlost + 1
              !
           END IF
        END DO
        !  
        ! assign remaining bands so as to maximise overlap
        !
        DO ibnd=1,nbnd
           !
           ! for unassigned bands ...
           !
           IF ( closest_band(ibnd,ik) == -1 ) THEN
              psmax = 0.0_dp
              jmax  = 0
              DO jbnd = 1, nbnd
                 !
                 ! check if this band was already assigne ...
                 !
                 IF ( done(jbnd) > 0 ) CYCLE
                 pscur = real(ps(ibnd,jbnd))**2+aimag(ps(ibnd,jbnd))**2
                 IF ( pscur > psmax ) THEN
                    psmax = pscur
                    jmax = jbnd
                 END IF
              END DO
              closest_band(ibnd,ik) = jmax
              done(jmax) = 1
           END IF
        END DO
     ENDIF
     !
     IF ( ik < nks2 .AND. .NOT. no_overlap ) THEN
        !
        ! compute S\psi_k to be used at next k-point
        !
        CALL s_psi( npwx, npw, nbnd, evc, spsi )
        !
     END IF
     !
  ENDDO
  !
  IF (noncolin) CALL poolrecover(sigma_avg,4*nbnd,nkstot,nks)
  !
  CALL deallocate_bec_type(becp)
  DEALLOCATE ( psi, spsi, ps )
  !
  IF ( ionode ) THEN
     !
     ! Re-order eigenvalues according to overlap
     ! (ibnd=band index, jbnd=re-ordered index)
     !
     ALLOCATE (et_(nbnd,nkstot))
     DO ibnd=1,nbnd
        jbnd = ibnd
        DO ik=nks1tot,nks2tot
           et_(ibnd,ik) = et(jbnd,ik)* rytoev
           jbnd = closest_band(jbnd,ik)
        END DO
     END DO
     !
     CALL punch_plottable_bands ( filband, nks1tot, nks2tot, nkstot, nbnd, &
                                  xk, et_ )
     !
     DO ik=nks1tot,nks2tot
        IF (ik == nks1) THEN
           WRITE (iunpun, '(" &plot nbnd=",i4,", nks=",i6," /")') &
             nbnd, nks2tot-nks1tot+1
           DO ipol=1,4
              IF (lsigma(ipol)) WRITE(iunpun_sigma(ipol), &
                            '(" &plot nbnd=",i4,", nks=",i6," /")') &
                             nbnd, nks2tot-nks1tot+1
           ENDDO
        ENDIF
        WRITE (iunpun, '(10x,3f10.6)') xk(1,ik),xk(2,ik),xk(3,ik)
        WRITE (iunpun, '(10f9.3)') (et_(ibnd, ik), ibnd = 1, nbnd)
        DO ipol=1,4
           IF (lsigma(ipol)) THEN
              WRITE (iunpun_sigma(ipol), '(10x,3f10.6)')            &
                                          xk(1,ik),xk(2,ik),xk(3,ik)
              WRITE (iunpun_sigma(ipol), '(10f9.3)')                &
                            (sigma_avg(ipol, ibnd, ik), ibnd = 1, nbnd)
           ENDIF
        ENDDO
     ENDDO
     !
     DEALLOCATE ( et_ )
     !
  ENDIF
  !
  IF (noncolin) DEALLOCATE (sigma_avg)
  !
  IF ( ionode ) THEN
     CLOSE (iunpun)
     WRITE ( stdout, &
          '(5x,"Bands written to file ",A)') TRIM(filband)
     DO ipol=1,4
        IF (lsigma(ipol)) CLOSE(iunpun_sigma(ipol))
     ENDDO
  ENDIF
  !
  RETURN
  !
END SUBROUTINE punch_band

SUBROUTINE punch_band_2d(filband,spin_component)
!
!  This routine opens a file for each band and writes on output 
!  kx, ky, energy, 
!  kx, ky, energy
!  .., .., ..
!  where kx and ky are proportional to the length
!  of the vectors k_1 and k_2 specified in the input of the 2d plot.
!
!  The k points are supposed to be in the form
!  xk(i,j) = xk_0 + dkx *(i-1) + dky * (j-1)      1<i<n1, 1<j<n2
!
!  kx(i,j) = (i-1) |dkx|
!  ky(i,j) = (j-1) |dky|
!
   USE kinds, ONLY : DP
   USE constants, ONLY : eps8, rytoev
   USE lsda_mod,  ONLY : nspin
   USE klist, ONLY : xk, nkstot, nks
   USE wvfct, ONLY : et, nbnd
   USE io_files, ONLY : iuntmp
   USE io_global, ONLY : ionode, ionode_id
   USE mp, ONLY : mp_bcast
   USE mp_world, ONLY : world_comm

   IMPLICIT NONE
   CHARACTER(LEN=256),INTENT(IN) :: filband
   INTEGER, INTENT(IN) :: spin_component
   REAL(DP) :: xk0(3), xk1(3), xk2(3), dkx(3), dky(3), xkdum(3), mdkx, mdky
   INTEGER :: n1, n2
   INTEGER :: ik, i, i1, i2, ibnd, ijk, start_k, last_k, nks_eff, j, ios
   CHARACTER(LEN=256) :: filename
   CHARACTER(LEN=6), EXTERNAL :: int_to_char
   REAL(DP), ALLOCATABLE :: xk_collect(:,:), et_collect(:,:)
   
   ALLOCATE(xk_collect(3,nkstot))
   ALLOCATE(et_collect(nbnd,nkstot))
   CALL poolcollect(    3, nks, xk, nkstot, xk_collect)
   CALL poolcollect( nbnd, nks, et, nkstot, et_collect)

   start_k=1
   last_k=nkstot
   nks_eff=nkstot
   IF (nspin==2) THEN
      nks_eff=nkstot/2
      IF (spin_component==1) THEN
         start_k=1
         last_k=nks_eff
      ELSE
         start_k=nks_eff+1
         last_k=nkstot
      ENDIF
   ENDIF
!
!  Determine xk0
!
   xk0(:)=xk_collect(:,start_k)
!
! Determine dkx
!
   dky(:)=xk_collect(:,start_k+1)-xk0(:)
!
! Determine n2 and dky
!

loop_k:  DO j=start_k+2, nkstot
     xkdum(:)=xk0(:)+(j-1)*dky(:)
     IF (ABS(xk_collect(1,j)-xkdum(1))>eps8.OR.   &
         ABS(xk_collect(2,j)-xkdum(2))>eps8.OR.   &
         ABS(xk_collect(3,j)-xkdum(3))>eps8) THEN    
         n2=j-1
         dkx(:)=xk_collect(:,j)-xk0(:)
         EXIT loop_k
     ENDIF
  ENDDO  loop_k
  n1=nks_eff/n2
  IF (n1*n2 /= nks_eff) CALL errore('punch_band_2d',&
                                    'Problems with k points',1)
  mdkx = sqrt( dkx(1)**2 + dkx(2)**2 + dkx(3)**2 )
  mdky = sqrt( dky(1)**2 + dky(2)**2 + dky(3)**2 )
!   
!  write the output, a band per file
!
  DO ibnd=1,nbnd
     filename=TRIM(filband) // '.' // TRIM(int_to_char(ibnd))
     IF (ionode) &
     open(unit=iuntmp,file=filename,status='unknown', err=100, iostat=ios)
     CALL mp_bcast(ios,ionode_id, world_comm)
100  CALL errore('punch_band_2d','Problem opening outputfile',ios)
     ijk=0
     DO i1=1,n1
        DO i2=1,n2
           ijk=ijk+1
           IF (ionode) &
           WRITE(iuntmp,'(3f16.6)') mdkx*(i1-1), mdky*(i2-1), &
                                    et_collect(ibnd,ijk)*rytoev
        ENDDO 
     ENDDO
     IF (ionode) CLOSE(unit=iuntmp,status='KEEP')
  ENDDO

  DEALLOCATE(xk_collect)
  DEALLOCATE(et_collect)

  RETURN
  END
!
!----------------------------------------------------------------------------
SUBROUTINE punch_plottable_bands ( filband, nks1tot, nks2tot, nkstot, nbnd, &
                                   xk, et )
  !---------------------------------------------------------------------------
  !
  USE kinds, ONLY : dp
  IMPLICIT NONE
  CHARACTER(LEN=*), INTENT(IN) :: filband
  INTEGER, INTENT(IN) :: nks1tot, nks2tot, nkstot, nbnd
  REAL(dp), INTENT(IN) :: xk(3,nkstot), et(nbnd,nks1tot)
  !
  INTEGER, PARAMETER :: max_lines = 100, stdout=6, iunpun0=18
  INTEGER:: ios, i, n, nlines, npoints(max_lines), point(max_lines)
  LOGICAL :: high_symmetry(nkstot), opnd
  REAL(dp):: k1(3), k2(3), kx(nkstot), ps, dxmod, dxmod_save
  !
  !
  IF ( nks1tot < 1 .OR. nks2tot > nkstot .OR. nkstot < 1 .OR. nbnd < 1 ) THEN
     CALL infomsg('punch_plottable_bands','incorrect input data, exiting')
     RETURN
  END IF
  ios = 0
  OPEN (unit = iunpun0, file = TRIM(filband)//'.gnu', status = 'unknown',&
       form = 'formatted', iostat = ios)
  IF ( ios /= 0 ) THEN
     WRITE ( stdout, &
          '(/,5x,"Error opening plottable file ",A)') TRIM(filband)//'.gnu'
     RETURN
  END IF
  !
  !  Find high-symmetry points (poor man's algorithm)
  !
  DO n=nks1tot,nks2tot
     IF (n==nks1tot .OR. n==nks2tot) THEN
        high_symmetry(n) = .true.
     ELSE
        k1(:) = xk(:,n) - xk(:,n-1)
        k2(:) = xk(:,n+1) - xk(:,n)
        IF ( k1(1)*k1(1) + k1(2)*k1(2) + k1(3)*k1(3) < 1.0d-8 .OR. &
             k2(1)*k2(1) + k2(2)*k2(2) + k2(3)*k2(3) < 1.0d-8 ) THEN
           CALL infomsg('punch_plottable_bands','two consecutive same k, exiting')
           RETURN
        END IF
        ps = ( k1(1)*k2(1) + k1(2)*k2(2) + k1(3)*k2(3) ) / &
             sqrt( k1(1)*k1(1) + k1(2)*k1(2) + k1(3)*k1(3) ) / &
             sqrt( k2(1)*k2(1) + k2(2)*k2(2) + k2(3)*k2(3) )
        high_symmetry(n) = (ABS(ps-1.d0) >1.0d-4)
        !
        !  The gamma point is a high symmetry point
        !
        IF (xk(1,n)**2+xk(2,n)**2+xk(3,n)**2 < 1.0d-8) high_symmetry(n)=.true.
        !
        !   save the typical length of dk
        !
        IF (n==nks1tot+1) dxmod_save = sqrt( k1(1)**2 + k1(2)**2 + k1(3)**2)
     ENDIF
  ENDDO

  kx(nks1tot) = 0.0_dp
  DO n=nks1tot+1,nks2tot
     dxmod=sqrt ( (xk(1,n)-xk(1,n-1))**2 + &
                  (xk(2,n)-xk(2,n-1))**2 + &
                  (xk(3,n)-xk(3,n-1))**2 )
     IF (dxmod > 5*dxmod_save) THEN
        !
        !   A big jump in dxmod is a sign that the point xk(:,n) and xk(:,n-1)
        !   are quite distant and belong to two different lines. We put them on
        !   the same point in the graph 
        !
        kx(n)=kx(n-1)
     ELSEIF (dxmod > 1.d-4) THEN
        !
        !  This is the usual case. The two points xk(:,n) and xk(:,n-1) are in
        !  the same path.
        !
        kx(n) = kx(n-1) +  dxmod
        dxmod_save = dxmod
     ELSE
        !
        !  This is the case in which dxmod is almost zero. The two points
        !  coincide in the graph, but we do not save dxmod.
        !
        kx(n) = kx(n-1) +  dxmod
        !
     ENDIF
  ENDDO
  !
  !  Now we compute how many paths there are: nlines
  !  The first point of this path: point(iline)
  !  How many points are in each path: npoints(iline)
  !
  DO n=nks1tot,nks2tot
     IF (high_symmetry(n)) THEN
        IF (n==nks1tot) THEN
           !
           !   first point. Initialize the number of lines, and the number of point
           !   and say that this line start at the first point
           !
           nlines=1
           npoints(1)=1
           point(1)=1
        ELSEIF (n==nks2tot) THEN
           !
           !    Last point. Here we save the last point of this line, but
           !    do not increase the number of lines
           !
           npoints(nlines) = npoints(nlines)+1
           point(nlines+1)=n
        ELSE
           !
           !   Middle line. The current line has one more points, and there is
           !   a new line that has to be initialized. It has one point and its
           !   first point is the current k.
           !
           npoints(nlines) = npoints(nlines)+1
           nlines=nlines+1
           IF (nlines>max_lines) THEN
              CALL infomsg('punch_plottable_bands','too many lines, exiting')
              RETURN
           END IF
           npoints(nlines) = 1
           point(nlines)=n
        ENDIF
        !
        WRITE( stdout,'(5x,"high-symmetry point: ",3f7.4,&
                         &"   x coordinate",f9.4)') (xk(i,n),i=1,3), kx(n)
     ELSE
        !
        !   This k is not an high symmetry line so we just increase the number
        !   of points of this line.
        !
        npoints(nlines) = npoints(nlines)+1
     ENDIF
  ENDDO
  !
  DO i=1,nbnd
     WRITE (iunpun0,'(2f10.4)') (kx(n), et(i,n),n=nks1tot,nks2tot)
     WRITE (iunpun0,*)
  ENDDO
  !
  WRITE ( stdout, &
       '(/,5x,"Plottable bands (eV) written to file ",A)') TRIM(filband)//'.gnu'
  CLOSE(unit=iunpun0, STATUS='KEEP')

  RETURN
  !
END SUBROUTINE punch_plottable_bands

