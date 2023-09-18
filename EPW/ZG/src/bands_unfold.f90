! Copyright (C) 2020-2021 Marios Zacharias, Feliciano Giustino 
!                                                                            
! This file is distributed under the terms of the GNU General Public         
! License. See the file `LICENSE' in the root directory of the               
! present distribution, or http://www.gnu.org/copyleft.gpl.txt .
!
! This is a modified a routine of bands.f90. Changes can be traced 
! by "mz_b" and "mz_e". In this routine we perform band structure 
! unfolding based on the wavefunctions calculated from a supercell 
! calculation. Theory follows the work in: 
!
! https://doi.org/10.1103/PhysRevB.85.085201 by Popescu and Zunger. 
! 
! If you are using this routine also cite: 
!
! https://doi.org/10.1103/PhysRevResearch.2.013357
!
!-----------------------------------------------------------------------
PROGRAM do_bands
  !-----------------------------------------------------------------------
  !
  ! See files INPUT_BANDS.* in Doc/ directory for usage
  ! 
  !
  USE io_files,  ONLY : prefix, tmp_dir
  USE mp_global, ONLY : mp_startup
  USE mp_pools,  ONLY : npool
  USE control_flags, ONLY : gamma_only
  USE environment,   ONLY : environment_start, environment_end
  USE wvfct,     ONLY : nbnd
  USE klist,     ONLY : nkstot, two_fermi_energies
  USE noncollin_module, ONLY : noncolin, i_cons
  USE lsda_mod,  ONLY : nspin
  USE io_global, ONLY : ionode, ionode_id, stdout
  USE mp,        ONLY : mp_bcast
  USE mp_images,  ONLY : intra_image_comm
  ! mz_b
  USE kinds,     ONLY : dp
  USE mp_bands,  ONLY : nproc_bgrp
  ! mz_e
  !
  IMPLICIT NONE
  !
  CHARACTER(LEN=256), EXTERNAL :: trimcheck
  !
  CHARACTER (len=256) :: filband, filp, outdir
  LOGICAL :: lsigma(4), lsym, lp, no_overlap, plot_2d
  INTEGER :: spin_component, firstk, lastk
  ! mz_b dim1,dim2,dim3
  INTEGER  :: dim1,dim2,dim3 
  logical  :: poors_man 
  ! mz_e 
  INTEGER :: ios
  !
  ! mz_b adds dim1,dim2,dim3 and , 
  NAMELIST / bands / outdir, prefix, filband, filp, spin_component, lsigma,&
                     lsym, lp, filp, firstk, lastk, no_overlap, plot_2d, &
                     dim1, dim2, dim3, poors_man 
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
  ! mz_b 
  dim1 = 1
  dim2 = 1
  dim3 = 1
  poors_man = .true.
  ! mz_e
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
  CALL mp_bcast( ios, ionode_id, intra_image_comm )
  IF (ios /= 0) CALL errore ('bands', 'reading bands namelist', abs(ios) )
  !
  ! ... Broadcast variables
  !
  CALL mp_bcast( tmp_dir, ionode_id, intra_image_comm )
  CALL mp_bcast( prefix, ionode_id, intra_image_comm )
  CALL mp_bcast( filband, ionode_id, intra_image_comm )
  CALL mp_bcast( filp, ionode_id, intra_image_comm )
  CALL mp_bcast( spin_component, ionode_id, intra_image_comm )
  CALL mp_bcast( firstk, ionode_id, intra_image_comm )
  CALL mp_bcast( lastk, ionode_id, intra_image_comm )
  CALL mp_bcast( lp, ionode_id, intra_image_comm )
  CALL mp_bcast( lsym, ionode_id, intra_image_comm )
  CALL mp_bcast( lsigma, ionode_id, intra_image_comm )
  CALL mp_bcast( no_overlap, ionode_id, intra_image_comm )
  CALL mp_bcast( plot_2d, ionode_id, intra_image_comm )
  !mz_b
  CALL mp_bcast( dim1, ionode_id, intra_image_comm )
  CALL mp_bcast( dim2, ionode_id, intra_image_comm )
  CALL mp_bcast( dim3, ionode_id, intra_image_comm )
  CALL mp_bcast( poors_man, ionode_id, intra_image_comm )
  ! 
  ! 

  IF (plot_2d) THEN
     lsym=.false.
     lp=.false.
     no_overlap=.true.
  ENDIF
  IF (lsym) no_overlap=.true.
  IF ( npool > 1 .and. poors_man ) CALL errore('bands_unfold', &
                                             'pools not implemented',npool)
  IF ( npool > 1 .and..not.(lsym.or.no_overlap)) CALL errore('bands_unfold', &
                                             'pools not implemented',npool)
  IF ( spin_component < 1 .OR. spin_component > 2 ) &
     CALL errore('bands','incorrect spin_component',1)
  !
  !   Now allocate space for pwscf variables, read and check them.
  !
  CALL read_file()
  !
  IF (gamma_only) CALL errore('bands','gamma_only case not implemented',1)
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
     !mz_b
     CALL punch_band(filband,spin_component,lsigma,no_overlap,dim1,dim2,dim3, &
                     poors_man)
 !! uncomment lsym because takes time
!!     IF (lsym) CALL sym_band(filband,spin_component,firstk,lastk)
!!
     !mz_e
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
SUBROUTINE punch_band (filband, spin_component, lsigma, no_overlap,dim1,dim2,dim3, &
                      poors_man)
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
  USE constants,            ONLY : rytoev, pi
  USE gvect,                ONLY : g, ngm
  USE klist,                ONLY : xk, nks, nkstot, ngk, igk_k
  USE io_files,             ONLY : iunpun, nwordwfc, iunwfc
  USE wvfct,                ONLY : nbnd, et, npwx
  USE uspp,                 ONLY : nkb, vkb
  USE uspp_param,           ONLY : upf, nh, nhm
  USE noncollin_module,     ONLY : noncolin, npol
  USE wavefunctions, ONLY : evc
  USE io_global,            ONLY : ionode, ionode_id, stdout
  USE mp,                   ONLY : mp_bcast
  USE mp_images,             ONLY : intra_image_comm
  USE becmod,               ONLY : calbec, bec_type, allocate_bec_type, &
                                   deallocate_bec_type, becp
  !mz_b
  USE cell_base,            ONLY : at, bg
  USE mp_bands,             ONLY : intra_bgrp_comm
  USE mp,                   ONLY : mp_sum 
  USE uspp_init,            ONLY : init_us_2
  !mz_e

  IMPLICIT NONE
  CHARACTER (len=*) :: filband
  INTEGER, INTENT(IN) :: spin_component
  LOGICAL, INTENT(IN) :: lsigma(4), no_overlap

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
  !mz_b
  logical, intent(in)      :: poors_man
  INTEGER, intent(in)      :: dim1,dim2,dim3
  REAL(DP), ALLOCATABLE    :: g_mz(:,:)
  INTEGER                  :: ctr,ctr2, kbnd 
  INTEGER                  :: i_mz, ig_mz, kkx, kky, kkz
  REAL(DP), ALLOCATABLE    :: P_mk(:,:), et_mz(:) !!, 
  COMPLEX(dp), ALLOCATABLE :: psi_mz(:,:)
  CHARACTER(len=256)       :: filename_mz,pointer_mz
  INTEGER                  :: index_mz! those are from file  
  COMPLEX(DP), ALLOCATABLE :: evc_new(:,:), evc_new2(:,:)
  REAL(DP), PARAMETER      :: eps = 0.00001d0
  COMPLEX(DP)              :: pro_mz ! pro_mz for adding becp contribution if noncolinear
  COMPLEX(DP), ALLOCATABLE :: pro_mz_t(:)
  COMPLEX(DP)              :: cgracsc_nc, cgracsc
  ! becp   : <psi|beta> at current  k-point
  TYPE(bec_type)::  becp_mz
  ! mz_adds  pro_mz_t(:,:) ! In pro_mz_t we get \sum_i Q_ik <\phi_nk|\beta_ik><\beta_ik|\phi_nk> 
  !see --> https://doi.org/10.1103/PhysRevB.41.7892
  !mz_e
  

  
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
  CALL mp_bcast( ios, ionode_id, intra_image_comm )
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
  ! mz_b
  CALL allocate_bec_type(nkb, nbnd, becp_mz)
  !
  !
  IF (poors_man) ALLOCATE(P_mk(nbnd,nkstot))
  IF (poors_man) P_mk(:,:) = 0.0d0
  ! mz_e
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
     !mz_b
     call cryst_to_cart( npw, g, at, -1 ) ! here we convert them to crystal coordinates
     ALLOCATE(g_mz(3,ngm))
     !
     g_mz(:,:) = g(:,:) ! save reciprocal lattice vectors G = m1*B1+m2*B2+m3*B3
     ! 
      ! Calculate the poor's man spectral weights
      IF ( poors_man ) THEN
        IF (noncolin) ALLOCATE(evc_new(npol*npwx,nbnd)) ! to compute becp contribution PAW or ultrasof
        ALLOCATE(evc_new2(npol*npwx,nbnd))
        IF (noncolin) evc_new (:,:) = 0
        evc_new2 = (0.d0,0.d0)
        !
        ctr=0              
        pwloop:  DO ig_mz = 1, npw
             IF    ( abs(mod(NINT(g(1,igk_k(ig_mz, ik))), dim1)) .EQ. 0  .AND. & ! to check if integer
                     abs(mod(NINT(g(2,igk_k(ig_mz, ik))), dim2)) .EQ. 0  .AND. & ! to check if integer
                     abs(mod(NINT(g(3,igk_k(ig_mz, ik))), dim3)) .EQ. 0  )  THEN ! to check if integer
                 ! P_mk as in eq.(4) of PRB 89, 041407(R) (2014) 
                 ! Separate treatment if we have spin-orbit coupling
                 ctr=ctr+1
                 IF (noncolin) THEN
                   DO jbnd = 1, nbnd ! bands loop
                      P_mk(jbnd,ik) = P_mk(jbnd,ik) + &
                                         CONJG(evc(ig_mz+npwx,jbnd))*evc(ig_mz+npwx,jbnd) &
                                         + CONJG(evc(ig_mz,jbnd))*evc(ig_mz,jbnd)
                     evc_new(ig_mz,jbnd) =  evc(ig_mz,jbnd) 
                     evc_new(ig_mz+npwx,jbnd) = evc(ig_mz+npwx,jbnd)!
                    END DO
                 ELSE
                   DO jbnd = 1, nbnd ! bands loop
                      P_mk(jbnd,ik) = P_mk(jbnd,ik) + &
                                     CONJG(evc(ig_mz,jbnd))*evc(ig_mz,jbnd)
                      evc_new2(ig_mz,jbnd) = evc(ig_mz,jbnd)
                   END DO
                 ENDIF ! noncol
             ENDIF ! if abs
         END DO pwloop
       call mp_sum( P_mk(:,ik), intra_bgrp_comm ) ! collect P_mk
       ! If PAW pseudos are used
       ! To compute the new scalar product with S matrix. The above subroutine
       ! gives the new <\beta_ik|\phi_nk>.
       ! This is for one q-point. I need to implement correctly the loop for manu p-points
        IF (noncolin) THEN
            ALLOCATE(pro_mz_t(nbnd))
            pro_mz_t(:) = 0.d0
            CALL calbec ( npw, vkb, evc_new, becp )
              DO jbnd = 1, nbnd ! bands loop
         !          CALL compute_sigma_avg(sigma_avg(1,1,ik),becp%nc,ik,lsigma)
         ! now compute new  \sum_i Q_ik <\phi_nk|\beta_ik><\beta_ik|\phi_nk>            
                   pro_mz = cgracsc_nc (nkb,becp%nc(1,1,jbnd), &
                             becp%nc(1,1,jbnd), nhm, ntyp, nh, &
                            nat, ityp, npol,upf)
         ! Note that pro_mz is parallelized over band in subroutine cgracsc
                   pro_mz_t(jbnd) = pro_mz_t(jbnd) + pro_mz
         ! computation of "scal" value in cgracsc.f90 which is equal to <\phi_nk|\phi_nk>. 
         ! We compute only \sum_i Q_ik <\phi_nk|\beta_ik><\beta_ik|\phi_nk> .
         ! This is substracted from the total spectral weight ! See Vanderbilt's paper eq.14
         ! "Soft self-consistent pseudopotentials in a generalized eigenvalue formalism"  
         !
              ENDDO
            P_mk(:,ik) = P_mk(:,ik) + dble(pro_mz_t(:))
            DEALLOCATE(pro_mz_t)
            DEALLOCATE(evc_new)
         ELSE
            ALLOCATE(pro_mz_t(nbnd))
            pro_mz_t(:) = 0.d0
              CALL calbec ( npw, vkb, evc_new2, becp_mz )
              ! it is very important to use evc_new2 and becp_mz
              DO jbnd = 1, nbnd ! bands loop
                  ! now compute new  \sum_i Q_ik <\phi_nk|\beta_ik><\beta_ik|\phi_nk> 
                     pro_mz = cgracsc(nkb,becp_mz%k(:,jbnd),becp_mz%k(:,jbnd), &
                                    nhm, ntyp, nh, nat, ityp, npw, upf)
          !
                     pro_mz_t(jbnd) = pro_mz_t(jbnd) + pro_mz
              ENDDO
            P_mk(:,ik) = P_mk(:,ik) + dble(pro_mz_t(:))
            DEALLOCATE(pro_mz_t)
         ENDIF ! noncolin
         DEALLOCATE(evc_new2)
      !
      END IF ! poors_man 
      !
      !
      ! Find all C(k+q + g + K)
     ALLOCATE ( psi_mz(npwx*npol,nbnd) )
     psi_mz(:,:) = (0.0_dp,0.0_dp)
     g_mz(:,:) = 2000.0 !g(:,:) ! to save and print gvecs
     ctr=0
     DO ig_mz = 1, npw
        ! Separate treatment if we have spin-orbit coupling
        ctr=ctr+1
        IF (noncolin) THEN
           g_mz(:,ctr) = g(:,igk_k(ig_mz,ik))
           DO jbnd = 1, nbnd ! bands loop
              psi_mz(ctr,jbnd) =  evc(ig_mz,jbnd)
              psi_mz(ctr+npwx,jbnd) = evc(ig_mz+npwx,jbnd)
           END DO
        ELSE
           g_mz(:,ctr) = g(:,igk_k(ig_mz,ik))
           DO jbnd = 1, nbnd ! bands loop
              psi_mz(ctr,jbnd) = evc(ig_mz,jbnd)
           END DO
        ENDIF ! noncol
     END DO
     ! Bring them back to cartesian, for the main loops
      call cryst_to_cart( npw, g, bg, 1 ) 
     !
     !
    DEALLOCATE(psi_mz,g_mz)
  ENDDO ! k-loop
  !call mp_sum( P_mk, intra_bgrp_comm ) ! collect P_mk
  !
  !
  IF (noncolin) CALL poolrecover(sigma_avg,4*nbnd,nkstot,nks)
  !
  IF ( ionode ) THEN
     write(pointer_mz,'(i2.2)') 1
     filename_mz = 'bands' // TRIM( pointer_mz ) // '.dat'
     OPEN (unit = 44, file = filename_mz, status = 'unknown', form = &
               'formatted', iostat = ios(0))
     WRITE (44, '(" &plot nbnd=",i4,", nks=",i6," /")') &
           nbnd, nks2tot-nks1tot+1
     DO ik=nks1tot,nks2tot
       !WRITE(*,*) "hiiii",  ik, nks1, nks2, nks1tot,nks2tot
       WRITE (44, '(10x,3f10.6)') dble(xk(1,ik)/dim1), dble(xk(2,ik)/dim2), dble(xk(3,ik)/dim3) !
       WRITE (44,'(10f14.6)') ( et(i_mz,ik)*rytoev, i_mz = 1, nbnd) ! write the energies of the supercell E_mK
     ENDDO
     CLOSE(44)
     ! 
     IF (poors_man) THEN
       filename_mz = 'spectral_weights' // TRIM( pointer_mz ) // '.dat'
       OPEN (unit = 25, file = filename_mz, status = 'unknown', form = &
               'formatted', iostat = ios(0))
       WRITE (25, '(" &plot nbnd=",i4,", nks=",i6," /")') & 
             nbnd, nks2tot-nks1tot+1
       DO ik=nks1tot,nks2tot
         ! We write the header for the output files and the new band structure
         WRITE (25, '(10x,3f10.6)') dble(xk(1, ik) / dim1), dble(xk(2, ik) / dim2),& 
                                  dble(xk(3, ik) / dim3) !xk(1,ik)/dim1, xk(2,ik)/dim2, xk(3,ik)/dim3
         WRITE (25,'(10f12.6)') ( P_mk(i_mz, ik), i_mz = 1, nbnd)
       ENDDO! ik
       CLOSE(25)          
      ! 
     END IF ! poorsman 
  ENDIF
  IF (poors_man) DEALLOCATE(P_mk)
  CALL deallocate_bec_type(becp_mz)
  !mz_e
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
     DEALLOCATE (closest_band)
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
        WRITE (iunpun, '(10f14.8)') (et_(ibnd, ik), ibnd = 1, nbnd)
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
  !
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
   USE mp_images, ONLY : intra_image_comm

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
     CALL mp_bcast(ios,ionode_id, intra_image_comm)
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
  REAL(dp), INTENT(IN) :: xk(3,nkstot), et(nbnd,nkstot)
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
!-----------------------------------------------------------------------
FUNCTION cgracsc_nc (nkb, bec1, bec2, nhm, ntyp, nh, nat, ityp, npol, upf)
  !-----------------------------------------------------------------------
  !
  !     This function computes the scalar product between two wavefunction
  !     and the S matrix of the US pseudopotential: <psi1 | S | psi2 >.
  !     It assumes that the product of psi1 with all the beta functions
  !     is in bec1, and the product of psi2 is in bec2.
  !
  !
  USE kinds
  USE uspp, ONLY: qq_so
  USE noncollin_module, ONLY: lspinorb
  USE pseudo_types, ONLY : pseudo_upf
  USE mp_bands,  ONLY : intra_bgrp_comm
  USE mp,         ONLY : mp_sum
  IMPLICIT NONE
  !
  !     here the dummy variables
  !

  INTEGER :: nkb, npol, nhm, ntyp, nat, ityp (nat), nh (ntyp)
  ! input: the number of beta functions
  ! input: the number of plane waves
  ! input: the maximum number of solid be
  ! input: the number of types of atoms
  ! input: the number of atoms
  ! input: the type of each atom
  ! input: the number of beta for each ty

  COMPLEX(DP) :: bec1 (nkb,npol), bec2 (nkb,npol), cgracsc_nc
  ! input: the first wavefunction
  ! input: the second wavefunction
  ! output: the value of the scalar produ

  TYPE(pseudo_upf) :: upf (ntyp)
  ! input: if true the pseudo is vanderb
  !
  !    Here the local variables
  !

  INTEGER :: ikb, jkb, na, np, ijkb0, ih, jh, ipol, jpol, ijh
  ! counter on total beta functions
  ! counter on total beta functions
  ! counter on atoms
  ! the pseudopotential of each atom
  ! auxiliary variable to compute ikb and jkb
  ! counter on solid beta functions
  ! counter on solid beta functions

  COMPLEX(DP) :: scal, zdotc
  !
  scal = 0.d0
  ijkb0 = 0
  DO np = 1, ntyp
  !   IF (upf(np)%tvanp ) THEN
        DO na = 1, nat
           IF (ityp (na) ==np) THEN
              DO ih = 1, nh (np)
                 ikb = ijkb0 + ih
                 DO jh = 1, nh (np)
                    jkb = ijkb0 + jh
                    IF (lspinorb) THEN
                       ijh=0
                       DO ipol=1,npol
                          DO jpol=1,npol
                            ijh=ijh+1
                          !  scal = scal + qq_so(ih,jh,1,np) * conjg( bec1(jkb, ipol) ) * bec2(ikb, jpol) &
                          !              + qq_so(ih,jh,2,np) * conjg( bec1(jkb, ipol) ) * bec2(ikb, jpol) &
                          !              + qq_so(ih,jh,3,np) * conjg( bec1(jkb, ipol) ) * bec2(ikb, jpol) &
                          !              + qq_so(ih,jh,4,np) * conjg( bec1(jkb, ipol) ) * bec2(ikb, jpol) 
                            scal=scal+qq_so(ih,jh,ijh,np)* &
                                      conjg(bec1(ikb,ipol))*bec2(jkb,jpol)
                          ENDDO
                       ENDDO
                    ENDIF
                 ENDDO
              ENDDO
              ijkb0 = ijkb0 + nh (np)
           ENDIF
        ENDDO
   !  ELSE
  !      DO na = 1, nat
  !         IF (ityp (na) ==np) ijkb0 = ijkb0 + nh (np)
  !      ENDDO
  !   ENDIF
  ENDDO

  cgracsc_nc = scal
  RETURN
END FUNCTION cgracsc_nc


FUNCTION cgracsc (nkb, bec1, bec2, nhm, ntyp, nh, nat, ityp, &
     npw, upf)
  !-----------------------------------------------------------------------
  !
  !     This function computes the scalar product between two wavefunction
  !     and the S matrix of the US pseudopotential: <psi1 | S | psi2 >.
  !     It assumes that the product of psi1 with all the beta functions
  !     is in bec1, and the product of psi2 is in bec2.
  !
  !
  USE kinds
  USE pseudo_types, ONLY : pseudo_upf
  USE mp_global,  ONLY : intra_bgrp_comm
  USE mp,         ONLY : mp_sum
  USE uspp, ONLY: qq_nt


  IMPLICIT NONE
  !
  !     here the dummy variables
  !

  INTEGER :: nkb, npw, nhm, ntyp, nat, ityp (nat), nh (ntyp)
  ! input: the number of beta functions
  ! input: the number of plane waves
  ! input: the maximum number of solid be
  ! input: the number of types of atoms
  ! input: the number of atoms
  ! input: the type of each atom
  ! input: the number of beta for each ty

  COMPLEX(DP) :: bec1 (nkb), bec2 (nkb),cgracsc
  ! input: the first wavefunction
  ! input: the second wavefunction
  ! output: the value of the scalar produ

!  real(DP) :: qq (nhm, nhm, ntyp)
  ! input: the q values defining S
  TYPE(pseudo_upf) :: upf (ntyp)
! input: if true the pseudo is vanderb
  !
  !    Here the local variables
  !

  INTEGER :: ikb, jkb, na, np, ijkb0, ih, jh
  ! counter on total beta functions
  ! counter on total beta functions
  ! counter on atoms
  ! the pseudopotential of each atom
  ! auxiliary variable to compute ikb and jkb
  ! counter on solid beta functions
  ! counter on solid beta functions

  COMPLEX(DP) :: scal
  !
  scal = 0.d0
  ijkb0 = 0
  DO np = 1, ntyp
    ! IF (upf(np)%tvanp ) THEN
        DO na = 1, nat
           IF (ityp (na) ==np) THEN
              DO ih = 1, nh (np)
                 ikb = ijkb0 + ih
                 DO jh = 1, nh (np)
                    jkb = ijkb0 + jh
                    scal = scal + qq_nt(ih,jh,np)*conjg(bec1(ikb))*bec2(jkb)
            !    conjg( becp%k(ikb,m) ) * becp2(jkb,n)
                 ENDDO
              ENDDO
              ijkb0 = ijkb0 + nh (np)
           ENDIF
        ENDDO
     !ELSE
     !   DO na = 1, nat
     !      IF (ityp (na) == np) ijkb0 = ijkb0 + nh (np)
     !   ENDDO
     !ENDIF
  ENDDO

  cgracsc = scal
  RETURN
END FUNCTION cgracsc
