  !                                                                            
  ! Copyright (C) 2010-2016 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino
  ! Copyright (C) 2007-2009 Jesse Noffsinger, Brad Malone, Feliciano Giustino  
  !                                                                            
  ! This file is distributed under the terms of the GNU General Public         
  ! License. See the file `LICENSE' in the root directory of the               
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .             
  !                                                                            
  !----------------------------------------------------------------------
  SUBROUTINE ephwann_shuffle( nqc, xqc )
  !---------------------------------------------------------------------
  !!
  !!  Wannier interpolation of electron-phonon vertex
  !!
  !!  Scalar implementation   Feb 2006
  !!  Parallel version        May 2006
  !!  Disentenglement         Oct 2006
  !!  Compact formalism       Dec 2006
  !!  Phonon irreducible zone Mar 2007
  !!
  !!  No ultrasoft now
  !!  No spin polarization
  !!
  !!  RM - add noncolin case
  !-----------------------------------------------------------------------
  !
  USE kinds,         ONLY : DP
  USE pwcom,         ONLY : nbnd, nks, nkstot, isk, &
                            et, xk, ef,  nelec
  USE cell_base,     ONLY : at, bg, omega, alat
  USE start_k,       ONLY : nk1, nk2, nk3
  USE ions_base,     ONLY : nat, amass, ityp, tau
  USE phcom,         ONLY : nq1, nq2, nq3, nmodes
  USE epwcom,        ONLY : nbndsub, fsthick, epwread, longrange,     &
                            epwwrite, ngaussw, degaussw, lpolar, lifc, lscreen, &
                            nbndskip, parallel_k, parallel_q, etf_mem, scr_typ, &
                            elecselfen, phonselfen, nest_fn, a2f, specfun_ph,   &
                            vme, eig_read, ephwrite, nkf1, nkf2, nkf3,          & 
                            efermi_read, fermi_energy, specfun_el, band_plot,   &
                            scattering, nstemp, int_mob, scissor, carrier,      &
                            iterative_bte, longrange, scatread, nqf1, prtgkk,   &
                            nqf2, nqf3, mp_mesh_k, restart, ncarrier, plselfen, &
                            specfun_pl, lindabs
  USE noncollin_module, ONLY : noncolin
  USE constants_epw, ONLY : ryd2ev, ryd2mev, one, two, eps2, zero, czero, &
                            twopi, ci, kelvin2eV
  USE io_files,      ONLY : prefix, diropn, tmp_dir
  USE io_global,     ONLY : stdout, ionode
  USE io_epw,        ONLY : lambda_phself, linewidth_phself, iunepmatwe,        &
                            iunepmatwp, crystal, iunepmatwp2
  USE elph2,         ONLY : nrr_k, nrr_q, cu, cuq, lwin, lwinq, irvec, ndegen_k,&
                            ndegen_q, wslen, chw, chw_ks, cvmew, cdmew, rdw,    &
                            epmatwp, epmatq, wf, etf, etf_k, etf_ks, xqf, xkf,  &
                            wkf, dynq, nqtotf, nkqf, epf17, nkf, nqf, et_ks,    &
                            ibndmin, ibndmax, lambda_all, dmec, dmef, vmef,     &
                            sigmai_all, sigmai_mode, gamma_all, epsi, zstar,    &
                            efnew, sigmar_all, zi_all, nkqtotf, eps_rpa,   &
                            nkqtotf, sigmar_all, zi_allvb, inv_tau_all, Fi_all, &
                            F_current, F_SERTA, inv_tau_allcb, zi_allcb
  USE transportcom,  ONLY : transp_temp, mobilityh_save, mobilityel_save, lower_bnd, &
                            upper_bnd, ixkqf_tr,  s_BZtoIBZ_full
  USE wan2bloch,     ONLY : dmewan2bloch, hamwan2bloch, dynwan2bloch, &
                            ephwan2blochp, ephwan2bloch, vmewan2bloch, &
                            dynifc2blochf, dynifc2blochc 
  USE bloch2wan,     ONLY : hambloch2wan, dmebloch2wan, dynbloch2wan, &
                            vmebloch2wan, ephbloch2wane, ephbloch2wanp, &
                            ephbloch2wanp_mem
  USE superconductivity, ONLY : write_ephmat, count_kpoints, kmesh_fine, &
                            kqmap_fine
  USE transport,     ONLY : transport_coeffs, iterativebte, scattering_rate_q
#ifdef __NAG
  USE f90_unix_io,   ONLY : flush
#endif
  USE mp,            ONLY : mp_barrier, mp_bcast, mp_sum
  USE io_global,     ONLY : ionode_id
  USE mp_global,     ONLY : inter_pool_comm, intra_pool_comm, root_pool
  USE mp_world,      ONLY : mpime, world_comm
#if defined(__MPI)
  USE parallel_include, ONLY: MPI_MODE_RDONLY, MPI_INFO_NULL
#endif
  !
  implicit none
  !
  INTEGER, INTENT (in) :: nqc
  !! number of qpoints in the coarse grid
  !
  REAL(kind=DP), INTENT (in) :: xqc(3,nqc)
  !! qpoint list, coarse mesh
  ! 
  ! Local  variables
  LOGICAL :: already_skipped
  !! Skipping band during the Wannierization
  LOGICAL :: exst
  !! If the file exist
  LOGICAL :: first_cycle
  !! Check wheter this is the first cycle after a restart. 
  LOGICAL :: first_time
  !! Check wheter this is the first timeafter a restart. 
  !
  CHARACTER (len=256) :: filint
  !! Name of the file to write/read 
  CHARACTER (len=30)  :: myfmt
  !! Variable used for formatting output
  ! 
  INTEGER, ALLOCATABLE :: ndegen_kk(:)
  !! Number of degeneresence with nrr_k bounds. 
  INTEGER, ALLOCATABLE :: ndegen_qq(:)
  !! Number of degeneresence with nrr_q bounds. 
  INTEGER, ALLOCATABLE :: irvec_kk(:,:)
  !! Irevec on the nrr_k bounds
  INTEGER, ALLOCATABLE :: irvec_qq(:,:)
  !! Irevec on the nrr_k bounds
  INTEGER :: ios
  !! integer variable for I/O control
  INTEGER :: iq 
  !! Counter on coarse q-point grid
  INTEGER :: iq_restart
  !! Counter on coarse q-point grid
  INTEGER :: ik
  !! Counter on coarse k-point grid
  INTEGER :: ikk
  !! Counter on k-point when you have paired k and q
  INTEGER :: ikq
  !! Paired counter so that q is adjacent to its k
  INTEGER :: ibnd
  !! Counter on band
  INTEGER :: jbnd
  !! Counter on band
  INTEGER :: imode
  !! Counter on mode
  INTEGER :: na
  !! Counter on atom
  INTEGER :: mu
  !! counter on mode
  INTEGER :: nu
  !! counter on mode
  INTEGER :: fermicount
  !! Number of states at the Fermi level
  INTEGER :: lrepmatw
  !! record length while reading file
  INTEGER :: ikx
  !! Counter on the coase k-grid
  INTEGER :: ikfx 
  !! Counter on the fine k-grid. 
  INTEGER :: xkk1, xkq1
  !! Integer of xkk when multiplied by nkf/nk
  INTEGER :: xkk2, xkq2
  !! Integer of xkk when multiplied by nkf/nk
  INTEGER :: xkk3, xkq3
  !! Integer of xkk when multiplied by nkf/nk
  INTEGER :: ir
  !! Counter for WS loop
  INTEGER :: nrws
  !! Number of real-space Wigner-Seitz
  INTEGER :: valueRSS(2)
  !! Return virtual and resisdent memory from system
  INTEGER :: ierr
  !! Error status
  INTEGER, PARAMETER :: nrwsx=200
  !! Maximum number of real-space Wigner-Seitz
  !  
  REAL(kind=DP) :: rdotk_scal
  !! Real (instead of array) for $r\cdot k$
  REAL(kind=DP) :: xxq(3)
  !! Current q-point 
  REAL(kind=DP) :: xxk(3)
  !! Current k-point on the fine grid
  REAL(kind=DP) :: xkk(3)
  !! Current k-point on the fine grid
  REAL(kind=DP) :: xkq(3)
  !! Current k+q point on the fine grid
  REAL(kind=DP) :: rws(0:3,nrwsx)
  !! Real-space wigner-Seitz vectors
  REAL(kind=DP) :: atws(3,3)
  !! Maximum vector: at*nq
  REAL(kind=DP), EXTERNAL :: efermig
  !! External function to calculate the fermi energy
  REAL(kind=DP), EXTERNAL :: efermig_seq
  !! Same but in sequential
  REAL(kind=DP), PARAMETER :: eps = 0.01/ryd2mev
  !! Tolerence
  REAL(kind=DP), ALLOCATABLE :: w2 (:)
  !! Interpolated phonon frequency
  REAL(kind=DP), ALLOCATABLE :: irvec_r (:,:)
  !! Wigner-Size supercell vectors, store in real instead of integer
  REAL(kind=DP), ALLOCATABLE :: rdotk(:)
  !! $r\cdot k$
  !
  COMPLEX(kind=DP) :: tablex (4*nk1+1,nkf1)
  !! Look-up table for the exponential (speed optimization) in the case of
  !! homogeneous grids.
  COMPLEX(kind=DP) :: tabley (4*nk2+1,nkf2)
  !! Look-up table for the exponential (speed optimization) in the case of
  !! homogeneous grids.
  COMPLEX(kind=DP) :: tablez (4*nk3+1,nkf3)
  !! Look-up table for the exponential (speed optimization) in the case of
  !! homogeneous grids.
  COMPLEX(kind=DP) :: tableqx (4*nk1+1,2*nkf1+1)
  !! Look-up table for the exponential (speed optimization) in the case of
  !! homogeneous grids.
  COMPLEX(kind=DP) :: tableqy (4*nk2+1,2*nkf2+1)
  !! Look-up table for the exponential (speed optimization) in the case of
  !! homogeneous grids.
  COMPLEX(kind=DP) :: tableqz (4*nk3+1,2*nkf3+1)
  !! Look-up table for the exponential (speed optimization) in the case of
  !! homogeneous grids.
  COMPLEX(kind=DP), ALLOCATABLE :: epmatwe  (:,:,:,:,:)
  !! e-p matrix  in wannier basis - electrons
  COMPLEX(kind=DP), ALLOCATABLE :: epmatwe_mem  (:,:,:,:)
  !! e-p matrix  in wannier basis - electrons (written on disk)
  COMPLEX(kind=DP), ALLOCATABLE :: epmatwef (:,:,:,:)
  !! e-p matrix  in el wannier - fine Bloch phonon grid
  COMPLEX(kind=DP), ALLOCATABLE :: epmatf( :, :, :)
  !! e-p matrix  in smooth Bloch basis, fine mesh
  COMPLEX(kind=DP), ALLOCATABLE :: cufkk ( :, :)
  !! Rotation matrix, fine mesh, points k
  COMPLEX(kind=DP), ALLOCATABLE :: cufkq ( :, :)
  !! the same, for points k+q
  COMPLEX(kind=DP), ALLOCATABLE :: uf( :, :)
  !! Rotation matrix for phonons
  COMPLEX(kind=DP), ALLOCATABLE :: bmatf ( :, :)
  !! overlap U_k+q U_k^\dagger in smooth Bloch basis, fine mesh
  COMPLEX(kind=DP), ALLOCATABLE :: cfac(:)
  !! Used to store $e^{2\pi r \cdot k}$ exponential 
  COMPLEX(kind=DP), ALLOCATABLE :: cfacq(:)
  !! Used to store $e^{2\pi r \cdot k+q}$ exponential

  ! Conductivity ------------
  INTEGER, PARAMETER :: maxiter = 300
  !! Maximum number of iteration in the case of iterative BTE
  INTEGER :: iter
  !! Current iteration number
  INTEGER :: itemp
  !! Temperature index
  INTEGER :: icbm
  !! Index of the CBM
  REAL(KIND=DP) :: error_h
  !! Error in the hole iterative BTE
  REAL(KIND=DP) :: error_el
  !! Error in the electron iterative BTE
  REAL(KIND=DP) :: etemp
  !! Temperature in Ry (this includes division by kb)
  REAL(KIND=DP) :: ef0(nstemp)
  !! Fermi level for the temperature itemp  
  REAL(KIND=DP) :: efcb(nstemp)
  !! Second Fermi level for the temperature itemp  
  REAL(kind=DP), ALLOCATABLE :: etf_all(:,:)
  !! Eigen-energies on the fine grid collected from all pools in parallel case
  REAL(KIND=DP), EXTERNAL :: fermicarrier
  !! Function that returns the Fermi level so that n=p (if int_mob = .true.)  
  ! -----------------
  ! 
  IF (nbndsub.ne.nbnd) &
       WRITE(stdout, '(/,14x,a,i4)' ) 'band disentanglement is used:  nbndsub = ', nbndsub
  !
  ALLOCATE ( cu ( nbnd, nbndsub, nks), & 
             cuq ( nbnd, nbndsub, nks), & 
             lwin ( nbnd, nks ), &
             lwinq ( nbnd, nks ), &
             irvec (3, 20*nk1*nk2*nk3), &
             ndegen_k (20*nk1*nk2*nk3), &
             ndegen_q (20*nq1*nq2*nq3), &
             wslen(20*nk1*nk2*nk3)      )
  !
  CALL start_clock ( 'ephwann' )
  !
  ! DBSP
  ! HERE loadkmesh
  IF ( epwread ) THEN
    !
    ! Might have been pre-allocate depending of the restart configuration 
    IF(ALLOCATED(tau))  DEALLOCATE( tau )
    IF(ALLOCATED(ityp)) DEALLOCATE( ityp )
    IF(ALLOCATED(w2))   DEALLOCATE( w2 )
    ! 
    ! We need some crystal info
    IF (mpime.eq.ionode_id) THEN
      !
      OPEN(unit=crystal,file='crystal.fmt',status='old',iostat=ios)
      READ (crystal,*) nat
      READ (crystal,*) nmodes
      READ (crystal,*) nelec
      READ (crystal,*) at
      READ (crystal,*) bg
      READ (crystal,*) omega
      READ (crystal,*) alat
      ALLOCATE( tau( 3, nat ) )
      READ (crystal,*) tau
      READ (crystal,*) amass
      ALLOCATE( ityp( nat ) )
      READ (crystal,*) ityp
      READ (crystal,*) isk
      READ (crystal,*) noncolin
      ! 
    ENDIF
    CALL mp_bcast (nat     , ionode_id, inter_pool_comm)
    CALL mp_bcast (nat     , root_pool, intra_pool_comm)  
    IF (mpime /= ionode_id) ALLOCATE( ityp( nat ) )
    CALL mp_bcast (nmodes  , ionode_id, inter_pool_comm)
    CALL mp_bcast (nmodes  , root_pool, intra_pool_comm)  
    CALL mp_bcast (nelec   , ionode_id, inter_pool_comm)
    CALL mp_bcast (nelec   , root_pool, intra_pool_comm)  
    CALL mp_bcast (at      , ionode_id, inter_pool_comm)
    CALL mp_bcast (at      , root_pool, intra_pool_comm)  
    CALL mp_bcast (bg      , ionode_id, inter_pool_comm)
    CALL mp_bcast (bg      , root_pool, intra_pool_comm)  
    CALL mp_bcast (omega   , ionode_id, inter_pool_comm)
    CALL mp_bcast (omega   , root_pool, intra_pool_comm)  
    CALL mp_bcast (alat    , ionode_id, inter_pool_comm)
    CALL mp_bcast (alat    , root_pool, intra_pool_comm)  
    IF (mpime /= ionode_id) ALLOCATE( tau( 3, nat ) )
    CALL mp_bcast (tau     , ionode_id, inter_pool_comm)
    CALL mp_bcast (tau     , root_pool, intra_pool_comm)  
    CALL mp_bcast (amass   , ionode_id, inter_pool_comm)
    CALL mp_bcast (amass   , root_pool, intra_pool_comm)  
    CALL mp_bcast (ityp    , ionode_id, inter_pool_comm)
    CALL mp_bcast (ityp    , root_pool, intra_pool_comm)  
    CALL mp_bcast (isk     , ionode_id, inter_pool_comm)
    CALL mp_bcast (isk     , root_pool, intra_pool_comm)  
    CALL mp_bcast (noncolin, ionode_id, inter_pool_comm)
    CALL mp_bcast (noncolin, root_pool, intra_pool_comm)  
    IF (mpime.eq.ionode_id) THEN
      CLOSE(crystal)
    ENDIF
    CALL mp_barrier(inter_pool_comm)
    ! 
  ELSE
    continue
  ENDIF
  !
  ALLOCATE( w2( 3*nat) )
  !
  ! determine Wigner-Seitz points
  !
  CALL wigner_seitz2 &
       ( nk1, nk2, nk3, nq1, nq2, nq3, nrr_k, nrr_q, irvec, wslen, ndegen_k, ndegen_q )
  ! 
  ALLOCATE ( ndegen_kk(nrr_k) )
  ALLOCATE ( ndegen_qq(nrr_q) )
  ALLOCATE ( irvec_kk(3,nrr_k) )
  ALLOCATE ( irvec_qq(3,nrr_q) )

  DO ir = 1, nrr_k
    ndegen_kk(ir) = ndegen_k(ir)
    irvec_kk(:,ir) = irvec(:,ir)
  ENDDO
  DO ir = 1, nrr_q
    ndegen_qq(ir) = ndegen_q(ir)
    irvec_qq(:,ir) = irvec(:,ir)
  ENDDO

  DEALLOCATE(ndegen_k) 
  DEALLOCATE(ndegen_q) 
  DEALLOCATE(irvec) 
  !
#ifndef __MPI  
  ! Open like this only in sequential. Otherwize open with MPI-open
  IF ((etf_mem == 1) .AND. (ionode)) THEN
    ! open the .epmatwe file with the proper record length
    lrepmatw   = 2 * nbndsub * nbndsub * nrr_k * nmodes
    filint    = trim(prefix)//'.epmatwp'
    CALL diropn (iunepmatwp, 'epmatwp', lrepmatw, exst)
  ENDIF
#endif
  ! 
  ! At this point, we will interpolate the Wannier rep to the Bloch rep 
  !
  IF ( epwread ) THEN
     !
     !  read all quantities in Wannier representation from file
     !  in parallel case all pools read the same file
     !
     CALL epw_read
     !
  ELSE !if not epwread (i.e. need to calculate fmt file)
     ! 
     IF ((etf_mem == 1) .AND. (ionode)) THEN
       lrepmatw   = 2 * nbndsub * nbndsub * nrr_k * nmodes
       filint    = trim(prefix)//'.epmatwe'
       CALL diropn (iunepmatwe, 'epmatwe', lrepmatw, exst)
#ifdef __MPI       
       filint    = trim(prefix)//'.epmatwp'
       CALL diropn (iunepmatwp, 'epmatwp', lrepmatw, exst)
#endif
     ENDIF
     !
     xxq = 0.d0 
     CALL loadumat &
          ( nbnd, nbndsub, nks, nkstot, xxq, cu, cuq, lwin, lwinq )  
     !
     ! ------------------------------------------------------
     !   Bloch to Wannier transform
     ! ------------------------------------------------------
     !
     ALLOCATE ( chw     ( nbndsub, nbndsub, nrr_k ),        &
          chw_ks  ( nbndsub, nbndsub, nrr_k ),        &
          cdmew   ( 3, nbndsub, nbndsub, nrr_k ),     &
          rdw     ( nmodes,  nmodes,  nrr_q ) )
     IF (vme) ALLOCATE(cvmew   ( 3, nbndsub, nbndsub, nrr_k ) )
     ! 
     ! SP : Let the user chose. If false use files on disk
     IF (etf_mem == 0) THEN
       ALLOCATE(epmatwe ( nbndsub, nbndsub, nrr_k, nmodes, nqc))
       ALLOCATE (epmatwp ( nbndsub, nbndsub, nrr_k, nmodes, nrr_q))
     ELSE
       ALLOCATE(epmatwe_mem ( nbndsub, nbndsub, nrr_k, nmodes))
       epmatwe_mem(:,:,:,:) = czero
     ENDIF
     !
     ! Hamiltonian
     !
     CALL hambloch2wan &
          ( nbnd, nbndsub, nks, nkstot, et, xk, cu, lwin, nrr_k, irvec_kk, wslen, chw )
     !
     ! Kohn-Sham eigenvalues
     !
     IF (eig_read) THEN
       WRITE (stdout,'(5x,a)') "Interpolating MB and KS eigenvalues"
       CALL hambloch2wan &
            ( nbnd, nbndsub, nks, nkstot, et_ks, xk, cu, lwin, nrr_k, irvec_kk, wslen, chw_ks )
     ENDIF
     !
     ! Dipole
     !
    ! CALL dmebloch2wan &
    !      ( nbnd, nbndsub, nks, nkstot, nkstot, dmec, xk, cu, nrr_k, irvec_kk, wslen )
     CALL dmebloch2wan &
          ( nbnd, nbndsub, nks, nkstot, dmec, xk, cu, nrr_k, irvec_kk, wslen, lwin )
     !
     ! Dynamical Matrix 
     !
     IF (.not. lifc) CALL dynbloch2wan &
                          ( nmodes, nqc, xqc, dynq, nrr_q, irvec_qq, wslen )
     !
     ! Transform of position matrix elements
     ! PRB 74 195118  (2006)
     IF (vme) CALL vmebloch2wan &
         ( nbnd, nbndsub, nks, nkstot, xk, cu, nrr_k, irvec_kk, wslen )
     !
     ! Electron-Phonon vertex (Bloch el and Bloch ph -> Wannier el and Bloch ph)
     !
     DO iq = 1, nqc
       !
       xxq = xqc (:, iq)
       !
       ! we need the cu again for the k+q points, we generate the map here
       !
       CALL loadumat ( nbnd, nbndsub, nks, nkstot, xxq, cu, cuq, lwin, lwinq )
       !
       DO imode = 1, nmodes
         !
         IF (etf_mem == 0) THEN 
           CALL ephbloch2wane &
             ( nbnd, nbndsub, nks, nkstot, xk, cu, cuq, &
             epmatq (:,:,:,imode,iq), nrr_k, irvec_kk, wslen, epmatwe(:,:,:,imode,iq) )
         ELSE
           CALL ephbloch2wane &
             ( nbnd, nbndsub, nks, nkstot, xk, cu, cuq, &
             epmatq (:,:,:,imode,iq), nrr_k, irvec_kk, wslen, epmatwe_mem(:,:,:,imode) )
           !
         ENDIF
         !
       ENDDO
       ! Only the master node writes 
       IF ((etf_mem == 1) .AND. (ionode)) THEN
         ! direct write of epmatwe for this iq 
         CALL rwepmatw ( epmatwe_mem, nbndsub, nrr_k, nmodes, iq, iunepmatwe, +1)       
         !   
       ENDIF   
       !
     ENDDO
     !
     ! Electron-Phonon vertex (Wannier el and Bloch ph -> Wannier el and Wannier ph)
     !
     ! Only master perform this task. Need to be parallelize in the future (SP)
     IF (ionode) THEN
       IF (etf_mem == 0) THEN
         CALL ephbloch2wanp &
           ( nbndsub, nmodes, xqc, nqc, irvec_kk, nrr_k, nrr_q, epmatwe )
       ELSE
          CALL ephbloch2wanp_mem &
           ( nbndsub, nmodes, xqc, nqc, irvec_kk, nrr_k, nrr_q, epmatwe_mem )
       ENDIF
     ENDIF
     !
     CALL mp_barrier(inter_pool_comm)
     !
     IF ( epwwrite ) THEN
        CALL epw_write 
        CALL epw_read 
     ENDIF
     !
  ENDIF
  !
  !
  IF ( ALLOCATED (epmatwe) ) DEALLOCATE (epmatwe)
  IF ( ALLOCATED (epmatwe_mem) ) DEALLOCATE (epmatwe_mem)
  IF ( ALLOCATED (epmatq) )  DEALLOCATE (epmatq)
  IF ( ALLOCATED (cu) )      DEALLOCATE (cu)
  IF ( ALLOCATED (cuq) )     DEALLOCATE (cuq)
  IF ( ALLOCATED (lwin) )    DEALLOCATE (lwin)
  IF ( ALLOCATED (lwinq) )   DEALLOCATE (lwinq)
  IF (etf_mem == 1) THEN
    CLOSE(iunepmatwe, status = 'delete')
  ELSE
    CLOSE(iunepmatwe)
  ENDIF
#ifdef __MPI
  CLOSE(iunepmatwp)
#endif
  ! 
  ! Check Memory usage
  CALL system_mem_usage(valueRSS)
  ! 
  WRITE(stdout, '(a)' )             '     ==================================================================='
  WRITE(stdout, '(a,i10,a)' ) '     Memory usage:  VmHWM =',valueRSS(2)/1024,'Mb'
  WRITE(stdout, '(a,i10,a)' ) '                   VmPeak =',valueRSS(1)/1024,'Mb'
  WRITE(stdout, '(a)' )             '     ==================================================================='
  WRITE(stdout, '(a)' )             '     '
  
  !
  !  At this point, we will interpolate the Wannier rep to the Bloch rep 
  !  for electrons, phonons and the ep-matrix
  !
  !  need to add some sort of parallelization (on g-vectors?)  what
  !  else can be done when we don't ever see the wfcs??
  !
  ! SP: k-point parallelization should always be efficient here
  IF (parallel_k) THEN
     CALL loadqmesh_serial
     CALL loadkmesh_para
  ELSEIF(parallel_q) THEN
     CALL loadkmesh_serial
     CALL loadqmesh_para
     ALLOCATE(etf_k ( nbndsub, nkqf))
  ELSE
     CALL errore('ephwann_shuffle', "parallel k and q not (yet) implemented",1)
  ENDIF
  !
  ALLOCATE ( epmatwef( nbndsub, nbndsub, nrr_k, nmodes),             &
       wf ( nmodes,  nqf ), etf ( nbndsub, nkqf),                   &
       etf_ks ( nbndsub, nkqf),                   &
       epmatf( nbndsub, nbndsub, nmodes), cufkk ( nbndsub, nbndsub), &
       cufkq ( nbndsub, nbndsub), uf ( nmodes, nmodes),              &
       bmatf( nbndsub, nbndsub), eps_rpa( nmodes) )
  !
  ! Need to be initialized
  epmatf(:,:,:) = czero
  ! allocate dipole matrix elements after getting grid size
  !
  ALLOCATE ( dmef(3, nbndsub, nbndsub, 2 * nkf) )
  IF (vme) ALLOCATE ( vmef(3, nbndsub, nbndsub, 2 * nkf) )
  !
  ALLOCATE(cfac(nrr_k))
  ALLOCATE(cfacq(nrr_k))
  ALLOCATE(rdotk(nrr_k))
  ! This is simply because dgemv take only real number (not integer)
  ALLOCATE(irvec_r(3,nrr_k))
  irvec_r = REAL(irvec_kk,KIND=dp)
  ! 
  ! SP: Create a look-up table for the exponential of the factor. 
  !     This can only work with homogeneous fine grids.
  IF ( (nkf1 >0) .AND. (nkf2 > 0) .AND. (nkf3 > 0) .AND. &
       (nqf1 >0) .AND. (nqf2 > 0) .AND. (nqf3 > 0) .AND. .NOT. mp_mesh_k .AND. .NOT. lscreen ) THEN
    ! Make a check   
    IF ((nqf1>nkf1) .or. (nqf2>nkf2) .or. (nqf3>nkf3)) &
            CALL errore('The fine q-grid cannot be larger than the fine k-grid',1)   
    ! Along x
    DO ikx = -2*nk1, 2*nk1
      DO ikfx = 0, nkf1-1
        !rdotk = twopi * ( xk(1)*irvec_kk(1,ir))
        rdotk_scal = twopi * ( (REAL(ikfx,kind=DP)/nkf1) * ikx )
        tablex(ikx+2*nk1+1,ikfx+1) = exp( ci*rdotk_scal )  
      ENDDO
    ENDDO
    ! For k+q
    DO ikx = -2*nk1, 2*nk1
      DO ikfx = 0, 2*nkf1
        rdotk_scal = twopi * ( (REAL(ikfx,kind=DP)/nkf1) * ikx )
        tableqx(ikx+2*nk1+1,ikfx+1) = exp( ci*rdotk_scal )
      ENDDO
    ENDDO
    ! Along y
    DO ikx = -2*nk2, 2*nk2
      DO ikfx = 0, nkf2-1
        rdotk_scal = twopi * ( (REAL(ikfx,kind=DP)/nkf2) * ikx )
        tabley(ikx+2*nk2+1,ikfx+1) = exp( ci*rdotk_scal )
      ENDDO
    ENDDO  
    ! For k+q
    DO ikx = -2*nk2, 2*nk2
      DO ikfx = 0, 2*nkf2
        rdotk_scal = twopi * ( (REAL(ikfx,kind=DP)/nkf2) * ikx )
        tableqy(ikx+2*nk2+1,ikfx+1) = exp( ci*rdotk_scal )
      ENDDO
    ENDDO
    ! Along z
    DO ikx = -2*nk3, 2*nk3
      DO ikfx = 0, nkf3-1
        rdotk_scal = twopi * ( (REAL(ikfx,kind=DP)/nkf3) * ikx )
        tablez(ikx+2*nk3+1,ikfx+1) = exp( ci*rdotk_scal )
      ENDDO
    ENDDO
    ! For k+q
    DO ikx = -2*nk3, 2*nk3
      DO ikfx = 0, 2*nkf3
        rdotk_scal = twopi * ( (REAL(ikfx,kind=DP)/nkf3) * ikx )
        tableqz(ikx+2*nk3+1,ikfx+1) = exp( ci*rdotk_scal )
      ENDDO
    ENDDO
  ENDIF
  !
  ! 
  ! ------------------------------------------------------
  ! Hamiltonian : Wannier -> Bloch (preliminary)
  ! ------------------------------------------------------
  !
  ! We here perform a preliminary interpolation of the hamiltonian
  ! in order to determine the fermi window ibndmin:ibndmax for later use.
  ! We will interpolate again afterwards, for each k and k+q separately
  !
  xxq = 0.d0
  !
  ! nkqf is the number of kpoints in the pool
  ! parallel_k case = nkqtotf/npool
  ! parallel_q case = nkqtotf
  !
  DO ik = 1, nkqf
     !
     xxk = xkf (:, ik)
     !
     IF ( 2*(ik/2).eq.ik ) THEN
        !
        !  this is a k+q point : redefine as xkf (:, ik-1) + xxq
        !
        CALL cryst_to_cart ( 1, xxq, at,-1 )
        xxk = xkf (:, ik-1) + xxq
        CALL cryst_to_cart ( 1, xxq, bg, 1 )
        !
     ENDIF
     !
     ! SP: Compute the cfac only once here since the same are use in both hamwan2bloch and dmewan2bloch
     ! + optimize the 2\pi r\cdot k with Blas
     CALL dgemv('t', 3, nrr_k, twopi, irvec_r, 3, xxk, 1, 0.0_DP, rdotk, 1 )
     cfac(:) = exp( ci*rdotk ) / ndegen_kk(:)
     ! 
     CALL hamwan2bloch &
          ( nbndsub, nrr_k, cufkk, etf (:, ik), chw, cfac)
     !
     !
  ENDDO
  !
  ! 27/06/2012 RM
  ! in the case when a random or uniform fine k-mesh is used
  ! calculate the Fermi level corresponding to the fine k-mesh 
  ! this Fermi level is then used as a reference in fermiwindow 
  ! 06/05/2014 CV
  ! calculate the Fermi level corresponding to the fine k-mesh
  ! or read it from input (Fermi level from the coarse grid 
  ! may be wrong or inaccurate)
  !
  WRITE(stdout,'(/5x,a,f10.6,a)') 'Fermi energy coarse grid = ', ef * ryd2ev, ' eV'
  !
  IF( efermi_read ) THEN
     !
     ef = fermi_energy
     WRITE(stdout,'(/5x,a)') repeat('=',67)
     WRITE(stdout, '(/5x,a,f10.6,a)') &
         'Fermi energy is read from the input file: Ef = ', ef * ryd2ev, ' eV'
     WRITE(stdout,'(/5x,a)') repeat('=',67)
     !
     ! SP: even when reading from input the number of electron needs to be correct
     already_skipped = .false.
     IF ( nbndskip .gt. 0 ) THEN
        IF ( .not. already_skipped ) THEN
           IF ( noncolin ) THEN
              nelec = nelec - one * nbndskip
           ELSE
              nelec = nelec - two * nbndskip
           ENDIF
           already_skipped = .true.
           WRITE(stdout,'(/5x,"Skipping the first ",i4," bands:")') nbndskip
           WRITE(stdout,'(/5x,"The Fermi level will be determined with ",f9.5," electrons")') nelec
        ENDIF
     ENDIF
     !      
  ELSEIF( band_plot ) THEN 
     !
     WRITE(stdout,'(/5x,a)') repeat('=',67)
     WRITE(stdout, '(/5x,"Fermi energy corresponds to the coarse k-mesh")')
     WRITE(stdout,'(/5x,a)') repeat('=',67) 
     !
  ELSE 
     ! here we take into account that we may skip bands when we wannierize
     ! (spin-unpolarized)
     ! RM - add the noncolin case
     already_skipped = .false.
     IF ( nbndskip .gt. 0 ) THEN
        IF ( .not. already_skipped ) THEN
           IF ( noncolin ) THEN 
              nelec = nelec - one * nbndskip
           ELSE
              nelec = nelec - two * nbndskip
           ENDIF
           already_skipped = .true.
           WRITE(stdout,'(/5x,"Skipping the first ",i4," bands:")') nbndskip
           WRITE(stdout,'(/5x,"The Fermi level will be determined with ",f9.5," electrons")') nelec
        ENDIF
     ENDIF
     !
     ! Fermi energy
     !  
     ! since wkf(:,ikq) = 0 these bands do not bring any contribution to Fermi level
     !  
     IF (parallel_k) efnew = efermig(etf, nbndsub, nkqf, nelec, wkf, degaussw, ngaussw, 0, isk)
     IF (parallel_q) THEN 
     IF (mpime .eq. ionode_id) THEN
       efnew = efermig_seq(etf, nbndsub, nkqf, nelec, wkf, degaussw, ngaussw, 0, isk)
       ! etf on the full k-grid is later required for selfen_phon_k          
       etf_k = etf
     ENDIF
     CALL mp_bcast (efnew, ionode_id, inter_pool_comm)
     CALL mp_bcast (etf_k, ionode_id, inter_pool_comm)
     ENDIF
     !
     WRITE(stdout, '(/5x,a,f10.6,a)') &
         'Fermi energy is calculated from the fine k-mesh: Ef = ', efnew * ryd2ev, ' eV'
     !
     ! if 'fine' Fermi level differs by more than 250 meV, there is probably something wrong
     ! with the wannier functions, or 'coarse' Fermi level is inaccurate
     IF (abs(efnew - ef) * ryd2eV .gt. 0.250d0 .and. (.not.eig_read) ) &
        WRITE(stdout,'(/5x,a)') 'Warning: check if difference with Fermi level fine grid makes sense'
     WRITE(stdout,'(/5x,a)') repeat('=',67)
     !
     ef=efnew
     !
  ENDIF
  !
  ! identify the bands within fsthick from the Fermi level
  ! (in shuffle mode this actually does not depend on q)
  !
  !
  CALL fermiwindow
  !
  !  xqf must be in crystal coordinates
  !
  ! this loops over the fine mesh of q points.
  ! if parallel_k then this is the entire q-list (nqftot)
  ! if parallel_q then this is nqftot/npool
  !
  ! ---------------------------------------------------------------------------------------
  ! ---------------------------------------------------------------------------------------
  IF (lifc) THEN
    !
    ! build the WS cell corresponding to the force constant grid
    !
    atws(:,1) = at(:,1)*DBLE(nq1)
    atws(:,2) = at(:,2)*DBLE(nq2)
    atws(:,3) = at(:,3)*DBLE(nq3)
    ! initialize WS r-vectors
    CALL wsinit(rws,nrwsx,nrws,atws)
  ENDIF
  !
  ! Open the ephmatwp file here
#if defined(__MPI)
  IF (etf_mem == 1) then
    ! Check for directory given by "outdir"
    !      
    filint = trim(tmp_dir)//trim(prefix)//'.epmatwp1'
    CALL MPI_FILE_OPEN(world_comm,filint,MPI_MODE_RDONLY,MPI_INFO_NULL,iunepmatwp2,ierr)
    IF( ierr /= 0 ) CALL errore( 'ephwann_shuffle', 'error in MPI_FILE_OPEN',1 )
    IF( parallel_q ) CALL errore( 'ephwann_shuffle', 'q-parallel+etf_mem = 1 is not supported',1 )
  ENDIF
#endif
  !
  IF (parallel_k) THEN
    !
    ! get the size of the matrix elements stored in each pool
    ! for informational purposes.  Not necessary
    !
    CALL mem_size(ibndmin, ibndmax, nmodes, nkf)
    !
    ! Fine mesh set of g-matrices.  It is large for memory storage
    ALLOCATE ( epf17 (ibndmax-ibndmin+1, ibndmax-ibndmin+1, nmodes, nkf) )
    ALLOCATE ( etf_all ( nbndsub, nkqtotf ) )
    ALLOCATE ( inv_tau_all (nstemp, ibndmax-ibndmin+1, nkqtotf/2) )
    inv_tau_all(:,:,:) = zero
    ALLOCATE ( zi_allvb (nstemp, ibndmax-ibndmin+1, nkqtotf/2) )
    zi_allvb(:,:,:) = zero
    ! 
    IF (int_mob .AND. carrier) THEN
      ALLOCATE ( inv_tau_allcb (nstemp, ibndmax-ibndmin+1, nkqtotf/2) )
      inv_tau_allcb(:,:,:) = zero
      ALLOCATE ( zi_allcb (nstemp, ibndmax-ibndmin+1, nkqtotf/2) )
      zi_allcb(:,:,:) = zero
    ENDIF
    !
    IF (iterative_bte) THEN
      ALLOCATE(Fi_all(3,ibndmax-ibndmin+1,nkqtotf/2))
      ! Current iterative F(i+1) function
      ALLOCATE(F_current(3,ibndmax-ibndmin+1,nkqtotf/2))
      ALLOCATE(F_SERTA(3,ibndmax-ibndmin+1,nkqtotf/2))
      Fi_all(:,:,:) = zero
      F_current(:,:,:) = zero
      F_SERTA(:,:,:) = zero
    ENDIF 
    ! 
    ! Define it only once for the full run. 
    CALL fkbounds( nkqtotf/2, lower_bnd, upper_bnd )
    ! 
    !  Start iteration index
    iter = 1
    !iter = 0
    ! Eerror for iterative BTE
    IF (int_mob .OR. (ncarrier > 1E5)) THEN
      error_el = 10_DP
    ELSE
      error_el = 0_DP
    ENDIF
    IF (int_mob .OR. (ncarrier < 1E5)) THEN
      error_h = 10_DP
    ELSE
      error_h = 0_DP
    ENDIF
    mobilityh_save = 0.0_DP
    mobilityel_save = 0.0_DP
    !
    ! Restart calculation
    iq_restart = 1
    first_cycle = .FALSE.
    first_time = .TRUE.
    IF (restart) THEN
      ! 
      IF ( elecselfen ) THEN
        IF ( .not. ALLOCATED (sigmar_all) ) ALLOCATE( sigmar_all(ibndmax-ibndmin+1, nkqtotf/2) )
        IF ( .not. ALLOCATED (sigmai_all) ) ALLOCATE( sigmai_all(ibndmax-ibndmin+1, nkqtotf/2) )
        IF ( .not. ALLOCATED (zi_all) )     ALLOCATE( zi_all(ibndmax-ibndmin+1, nkqtotf/2) )
        sigmar_all(:,:) = zero
        sigmai_all(:,:) = zero
        zi_all(:,:) = zero
        !
        CALL electron_read(iq_restart, nqf, nkqtotf/2, sigmar_all, sigmai_all, zi_all)
        !
      ENDIF
      IF ( scattering ) THEN
        ! 
        IF (int_mob .AND. carrier) THEN
          !
          ! Here inv_tau_all and inv_tau_allcb gets updated
          CALL tau_read(iq_restart, nqf, nkqtotf/2, .TRUE.)
          !
        ELSE
          ! Here inv_tau_all gets updated
          CALL tau_read(iq_restart, nqf, nkqtotf/2, .FALSE.)
          !
        ENDIF
        !
      ENDIF
      IF ( iterative_bte ) THEN
        ! 
        CALL F_read(iter, iq_restart, nqf,nkqtotf/2, error_h, error_el)
        ! 
        IF (int_mob .OR. (ncarrier < 1E5)) THEN
          IF ( error_h < eps2 ) WRITE(stdout,'(5x,a)') repeat('=',67)
          IF ( error_h < eps2 ) &
            WRITE( stdout,'(5x,"IBTE is converged with value for hole mobility of",1E18.6," "/)') mobilityh_save
          IF ( error_h < eps2 ) WRITE(stdout,'(5x,a)') repeat('=',67)
        ENDIF
        IF (int_mob .OR. (ncarrier > 1E5)) THEN
          IF ( error_el < eps2 ) WRITE(stdout,'(5x,a)') repeat('=',67)
          IF ( error_el < eps2 ) &
            WRITE( stdout,'(5x,"IBTE is converged with value for electron mobility of",1E18.6," "/)') mobilityel_save
          IF ( error_el < eps2 ) WRITE(stdout,'(5x,a)') repeat('=',67)
        ENDIF
        !
      ENDIF     
      !
      ! If you restart from reading a file. This prevent 
      ! the case were you restart but the file does not exist
      IF (iq_restart > 1) first_cycle = .TRUE.
      ! 
    ENDIF
    ! 
    ! Scatread assumes that you alread have done the full q-integration
    ! We just do one loop to get interpolated eigenenergies.  
    IF(scatread) iq_restart = nqf -1 
    !      
    DO WHILE ( (error_h > eps2) .OR. (error_el > eps2) )
      !
      IF ( iterative_bte ) THEN
        IF (iter==1 .OR. first_cycle) THEN
          WRITE(stdout,'(5x,a)') ' '
          WRITE(stdout,'(5x,a)') repeat('=',67)
          WRITE(stdout,'(5x,"Start solving iterative Boltzmann Transport Equation")')
          WRITE(stdout,'(5x,a/)') repeat('=',67)
        ENDIF
        WRITE(stdout,'(/5x,"Iteration number:", i10," "/)') iter
        ! 
        IF (nstemp > 1) CALL errore('ephwann_shuffle', &
            'Iterative BTE can only be done at 1 temperature, nstemp = 1.',1)  
        ! 
        IF (iter > maxiter) CALL errore('ephwann_shuffle', &
          'The iteration reached the maximum but did not converge. ',1)
        ! 
        ! Reading iter X if the file exist from a restart
        !CALL F_read(Fi_all, ibndmax-ibndmin+1, nkqtotf/2, iter)
        ! 
      ENDIF
      ! 
      !iter = iter +1
      !
      DO iq=iq_restart, nqf
         !   
         CALL start_clock ( 'ep-interp' )
         !
         ! In case of big calculation, show progression of iq (especially usefull when
         ! elecselfen = true as nothing happen during the calculation otherwise. 
         !
         IF ( .not. phonselfen) THEN 
           IF (MOD(iq,50) == 0) THEN
             WRITE(stdout, '(a,i10,a,i10)' ) '     Progression iq (fine) = ',iq,'/',nqf
           ENDIF
         ENDIF
         !
         xxq = xqf (:, iq)
         !
         ! ------------------------------------------------------
         ! dynamical matrix : Wannier -> Bloch
         ! ------------------------------------------------------
         !
         IF (.not. lifc) THEN
           CALL dynwan2bloch &
               ( nmodes, nrr_q, irvec_qq, ndegen_qq, xxq, uf, w2 )
         ELSE
           CALL dynifc2blochf ( nmodes, rws, nrws, xxq, uf, w2 )
         ENDIF
         !
         ! ...then take into account the mass factors and square-root the frequencies...
         !
         DO nu = 1, nmodes
            !
            ! wf are the interpolated eigenfrequencies
            ! (omega on fine grid)
            !
            IF ( w2 (nu) .gt. 0.d0 ) THEN
               wf(nu,iq) =  sqrt(abs( w2 (nu) ))
            ELSE
               wf(nu,iq) = -sqrt(abs( w2 (nu) ))
            ENDIF
            !
            DO mu = 1, nmodes
               na = (mu - 1) / 3 + 1
               uf (mu, nu) = uf (mu, nu) / sqrt(amass(ityp(na)))
            ENDDO
         ENDDO
         !
         ! --------------------------------------------------------------
         ! epmat : Wannier el and Wannier ph -> Wannier el and Bloch ph
         ! --------------------------------------------------------------
         !
         !DBSP
         !CALL start_clock ( 'cl2' )
         IF (.NOT. longrange) THEN
           CALL ephwan2blochp &
               ( nmodes, xxq, irvec_qq, ndegen_qq, nrr_q, uf, epmatwef, nbndsub, nrr_k )
         ENDIF
         !CALL stop_clock ( 'cl2' )
         !
         !
         !  number of k points with a band on the Fermi surface
         fermicount = 0
         !
         IF (lscreen) THEN
           IF (scr_typ == 0) CALL rpa_epsilon (xxq, wf(:,iq), nmodes, epsi, eps_rpa)
           IF (scr_typ == 1) CALL tf_epsilon (xxq, nmodes, epsi, eps_rpa)
         ENDIF
         ! this is a loop over k blocks in the pool
         ! (size of the local k-set)
         DO ik = 1, nkf
           !print*,'ik in ephwann ', ik
           !
           ! xkf is assumed to be in crys coord
           !
           ikk = 2 * ik - 1
           ikq = ikk + 1
           !
           xkk = xkf(:, ikk)
           xkq = xkk + xxq
           !
           ! SP: Compute the cfac only once here since the same are use in both hamwan2bloch and dmewan2bloch
           ! + optimize the 2\pi r\cdot k with Blas
           IF ( (nkf1 >0) .AND. (nkf2 > 0) .AND. (nkf3 > 0) .AND. &
              (nqf1 > 0) .AND. (nqf2 > 0) .AND. (nqf3 > 0) .AND. .NOT. mp_mesh_k .AND. .NOT. lscreen ) THEN
             ! We need to use NINT (nearest integer to x) rather than INT
             xkk1 = NINT(xkk(1)*(nkf1)) + 1
             xkk2 = NINT(xkk(2)*(nkf2)) + 1
             xkk3 = NINT(xkk(3)*(nkf3)) + 1
             xkq1 = NINT(xkq(1)*(nkf1)) + 1
             xkq2 = NINT(xkq(2)*(nkf2)) + 1
             xkq3 = NINT(xkq(3)*(nkf3)) + 1
             ! 
             ! SP: Look-up table is more effecient than calling the exp function.
             DO ir = 1, nrr_k
               cfac(ir) = ( tablex(irvec_kk(1,ir)+2*nk1+1,xkk1) *&
                       tabley(irvec_kk(2,ir)+2*nk2+1,xkk2) * tablez(irvec_kk(3,ir)+2*nk3+1,xkk3) ) / ndegen_kk(ir)
               cfacq(ir) = ( tableqx(irvec_kk(1,ir)+2*nk1+1,xkq1) *&
                       tableqy(irvec_kk(2,ir)+2*nk2+1,xkq2) * tableqz(irvec_kk(3,ir)+2*nk3+1,xkq3) ) /  ndegen_kk(ir)
             ENDDO
             !DBSP
             !IF ( (iq == 1) .and. (ik ==12)) THEN
             !  CALL dgemv('t', 3, nrr_k, twopi, irvec_r, 3, xkk, 1, 0.0_DP, rdotk, 1 )
             !  cfac1(:) = exp( ci*rdotk(:) ) / ndegen_k(:)
             !  CALL dgemv('t', 3, nrr_k, twopi, irvec_r, 3, xkq, 1, 0.0_DP, rdotk, 1 )
             !  cfacq1(:) = exp( ci*rdotk(:) ) / ndegen_k(:)
             !ENDIF
           ELSE
             CALL dgemv('t', 3, nrr_k, twopi, irvec_r, 3, xkk, 1, 0.0_DP, rdotk, 1 )
             cfac(:) = exp( ci*rdotk(:) ) / ndegen_kk(:)
             CALL dgemv('t', 3, nrr_k, twopi, irvec_r, 3, xkq, 1, 0.0_DP, rdotk, 1 )
             cfacq(:) = exp( ci*rdotk(:) ) / ndegen_kk(:)
           ENDIF
           !
           ! ------------------------------------------------------        
           ! hamiltonian : Wannier -> Bloch 
           ! ------------------------------------------------------
           !
           ! Kohn-Sham first, then get the rotation matricies for following interp.
           IF (eig_read) THEN
              CALL hamwan2bloch &
                ( nbndsub, nrr_k, cufkk, etf_ks (:, ikk), chw_ks, cfac)
              CALL hamwan2bloch &
                ( nbndsub, nrr_k, cufkq, etf_ks (:, ikq), chw_ks, cfacq)
           ENDIF
           !
           CALL hamwan2bloch &
                ( nbndsub, nrr_k, cufkk, etf (:, ikk), chw, cfac)
           CALL hamwan2bloch &
                ( nbndsub, nrr_k, cufkq, etf (:, ikq), chw, cfacq)
           !
           ! ------------------------------------------------------        
           !  dipole: Wannier -> Bloch
           ! ------------------------------------------------------        
           !
           CALL dmewan2bloch &
                ( nbndsub, nrr_k, cufkk, dmef (:,:,:, ikk), etf(:,ikk), etf_ks(:,ikk), cfac)
           CALL dmewan2bloch &
               ( nbndsub, nrr_k, cufkq, dmef (:,:,:, ikq), etf(:,ikq), etf_ks(:,ikq), cfacq)
           ! 
           ! ------------------------------------------------------        
           !  velocity: Wannier -> Bloch
           ! ------------------------------------------------------        
           !
           IF (vme) THEN
              IF (eig_read) THEN
                 CALL vmewan2bloch &
                      ( nbndsub, nrr_k, irvec_kk, ndegen_kk, xkk, cufkk, vmef (:,:,:, ikk), etf(:,ikk), etf_ks(:,ikk), chw_ks)
                 CALL vmewan2bloch &
                      ( nbndsub, nrr_k, irvec_kk, ndegen_kk, xkq, cufkq, vmef (:,:,:, ikq), etf(:,ikq), etf_ks(:,ikq), chw_ks)
              ELSE
                 CALL vmewan2bloch &
                      ( nbndsub, nrr_k, irvec_kk, ndegen_kk, xkk, cufkk, vmef (:,:,:, ikk), etf(:,ikk), etf_ks(:,ikk), chw)
                 CALL vmewan2bloch &
                      ( nbndsub, nrr_k, irvec_kk, ndegen_kk, xkq, cufkq, vmef (:,:,:, ikq), etf(:,ikq), etf_ks(:,ikq), chw)
              ENDIF
           ENDIF
           !
           IF (.NOT. scatread) THEN
             ! interpolate ONLY when (k,k+q) both have at least one band 
             ! within a Fermi shell of size fsthick 
             !
             IF ( (( minval ( abs(etf (:, ikk) - ef) ) < fsthick ) .and. ( minval ( abs(etf (:, ikq) - ef) ) < fsthick )) ) THEN
               !
               !  fermicount = fermicount + 1
               !
               ! --------------------------------------------------------------
               ! epmat : Wannier el and Bloch ph -> Bloch el and Bloch ph
               ! --------------------------------------------------------------
               !
               !
               ! SP: Note: In case of polar materials, computing the long-range and short-range term 
               !     separately might help speed up the convergence. Indeed the long-range term should be 
               !     much faster to compute. Note however that the short-range term still contains a linear
               !     long-range part and therefore could still be a bit more difficult to converge than 
               !     non-polar materials. 
               ! 
               IF (longrange) THEN
                 !      
                 epmatf = czero
                 !
               ELSE
                 !
                 CALL ephwan2bloch &
                   ( nbndsub, nrr_k, irvec_kk, ndegen_kk, epmatwef, xkk, cufkk, cufkq, epmatf, nmodes )
                 !
               ENDIF
               !
               IF (lpolar) THEN
                 !
                 CALL compute_umn_f( nbndsub, cufkk, cufkq, bmatf )
                 !
                 IF ( (abs(xxq(1)) > eps) .or. (abs(xxq(2)) > eps) .or. (abs(xxq(3)) > eps) ) THEN
                   !      
                   CALL cryst_to_cart (1, xxq, bg, 1)
                   CALL rgd_blk_epw_fine(nq1, nq2, nq3, xxq, uf, epmatf, &
                                         nmodes, epsi, zstar, bmatf, +1.d0)
                   CALL cryst_to_cart (1, xxq, at, -1)
                   !
                 ENDIF
                 !
               ENDIF
               ! 
               ! Store epmatf in memory
               !
               DO jbnd = ibndmin, ibndmax
                 DO ibnd = ibndmin, ibndmax
                   ! 
                   IF (lscreen) THEN
                      epf17(ibnd-ibndmin+1,jbnd-ibndmin+1,:,ik) = epmatf(ibnd,jbnd,:) / eps_rpa(:)
                   ELSE
                      epf17(ibnd-ibndmin+1,jbnd-ibndmin+1,:,ik) = epmatf(ibnd,jbnd,:)
                   ENDIF
                   !
                 ENDDO
               ENDDO
               !
               !if (ik==2) then
               !  do imode = 1, nmodes
               !    write(*,*) 'epmatf ',SUM((REAL(REAL(epmatf(:,:,imode))))**2)+SUM((REAL(AIMAG(epmatf(:,:,imode))))**2)
               !  enddo
               !endif
               !
             ENDIF
           ENDIF ! scatread 
           !
         ENDDO  ! end loop over k points
         !
         IF (prtgkk     ) CALL print_gkk( iq )
         IF (phonselfen ) CALL selfen_phon_q( iq )
         IF (elecselfen ) CALL selfen_elec_q( iq, first_cycle )
         IF (plselfen   ) CALL selfen_pl_q( iq )
         IF (nest_fn    ) CALL nesting_fn_q( iq )
         IF (specfun_el ) CALL spectral_func_q( iq )
         IF (specfun_ph ) CALL spectral_func_ph( iq )
         IF (specfun_pl ) CALL spectral_func_pl_q( iq )
         IF (ephwrite) THEN
            IF ( iq .eq. 1 ) THEN 
               CALL kmesh_fine
               CALL kqmap_fine
            ENDIF
            CALL write_ephmat( iq ) 
            CALL count_kpoints(iq)
         ENDIF
         ! 
         IF (.NOT. scatread) THEN
           ! 
           ! Indirect absorption ---------------------------------------------------------
           ! If Indirect absortpion, keep unshifted values:
           IF ( lindabs .AND. .NOT. scattering ) etf_ks(:,:) = etf(:,:)
           ! 
           ! Apply a scissor shift to CBM if required by user
           ! The shift is apply to k and k+q
           IF (ABS(scissor) > 0.000001) THEN
             IF ( noncolin ) THEN
               icbm = FLOOR(nelec/1.0d0) +1
             ELSE
               icbm = FLOOR(nelec/2.0d0) +1
             ENDIF
             !
             DO ik = 1, nkf
               ikk = 2 * ik - 1
               ikq = ikk + 1
               DO ibnd = icbm, nbndsub
                 ! 
                 etf (ibnd, ikk) = etf (ibnd, ikk) + scissor
                 etf (ibnd, ikq) = etf (ibnd, ikq) + scissor
               ENDDO
             ENDDO
             IF ( iq .eq. 1 ) THEN
               WRITE(stdout, '(5x,"Applying a scissor shift of ",f9.5," eV to the conduction states")' ) scissor * ryd2ev
             ENDIF
           ENDIF
           !  
           ! Indirect absorption
           IF ( lindabs .AND. .NOT. scattering )  CALL indabs(iq)  
           ! 
           ! Conductivity ---------------------------------------------------------
           IF (scattering) THEN
             !   
             ! If we want to compute intrinsic mobilities, call fermicarrier to 
             ! correctly positionned the ef0 level.
             ! This is only done once for iq = 0 
             IF ( iq == iq_restart ) THEN
               ! 
               DO itemp = 1, nstemp
                 ! 
                 etemp = transp_temp(itemp) 
                 WRITE(stdout, '(/5x,"Temperature ",f8.3," K")' ) etemp * ryd2ev / kelvin2eV
                 ! 
                 ! Small gap semiconductor. Computes intrinsic mobility by placing 
                 ! the Fermi level such that carrier density is equal for electron and holes
                 IF (int_mob .AND. .NOT. carrier) THEN               
                   !
                   ef0(itemp) = fermicarrier( etemp )
                   WRITE(stdout, '(5x,"Mobility Fermi level ",f10.6," eV")' )  ef0(itemp) * ryd2ev  
                   ! We only compute 1 Fermi level so we do not need the other
                   efcb(itemp) = 0
                   !   
                 ENDIF
                 ! 
                 ! Large bandgap semiconductor. Place the gap at the value ncarrier.
                 ! The user want both VB and CB mobilities. 
                 IF (int_mob .AND. carrier) THEN
                   ! 
                   ncarrier = - ABS(ncarrier) 
                   ef0(itemp) = fermicarrier( etemp )
                   WRITE(stdout, '(5x,"Mobility VB Fermi level ",f10.6," eV")' )  ef0(itemp) * ryd2ev 
                   ! 
                   ncarrier = ABS(ncarrier) 
                   efcb(itemp) = fermicarrier( etemp )
                   WRITE(stdout, '(5x,"Mobility CB Fermi level ",f10.6," eV")' )  efcb(itemp) * ryd2ev
                   !  
                 ENDIF   
                 ! 
                 ! User decide the carrier concentration and choose to only look at VB or CB  
                 IF (.NOT. int_mob .AND. carrier) THEN
                   ! SP: Determination of the Fermi level for intrinsic or doped carrier 
                   !     One also need to apply scissor before calling it.
                   ! 
                   ef0(itemp) = fermicarrier( etemp )               
                   WRITE(stdout, '(5x,"Mobility Fermi level ",f10.6," eV")' )  ef0(itemp) * ryd2ev
                   ! We only compute 1 Fermi level so we do not need the other
                   efcb(itemp) = 0
                   ! 
                 ENDIF
                 ! 
                 IF (.NOT. int_mob .AND. .NOT. carrier ) THEN
                   IF ( efermi_read ) THEN
                     !
                     ef0(itemp) = fermi_energy
                     !
                   ELSE !SP: This is added for efficiency reason because the efermig routine is slow
                     ef0(itemp) = efnew
                   ENDIF
                   ! We only compute 1 Fermi level so we do not need the other
                   efcb(itemp) = 0
                   !  
                 ENDIF
                 ! 
               ENDDO
               !
               ! 
             ENDIF ! iq=0
             !   
             IF ( iterative_bte) THEN
               ! First iteration is just SERTA
               IF (iter == 1) THEN 
                 !
                 CALL scattering_rate_q( iq, ef0, efcb, first_cycle )
                 !print*,'SUM(inv_tau_all) after ',SUM(inv_tau_all)
                 !
                 ! Compute the SERTA mobility for the first iteration
                 IF (iq == nqf) CALL transport_coeffs (ef0,efcb)
                 IF (iq == nqf) iter = iter + 1 
                 ! 
               ELSE
                 ! 
                 IF (int_mob .AND. carrier) THEN
                   call errore('ephwann_shuffle','The iterative solution cannot be solved with int_mob AND carrier at the moment',1)
                 ELSE
                   CALL iterativebte(iter, iq, ef0(1), error_h, error_el, first_cycle, first_time)
                 ENDIF   
                 !
                 IF (iq == nqf) iter = iter + 1 
               ENDIF
               !
             ELSE
               !
               CALL scattering_rate_q( iq, ef0, efcb, first_cycle )
               ! 
             ENDIF
             ! 
           ENDIF ! scattering
           ! --------------------------------------       
           !
           CALL stop_clock ( 'ep-interp' )
           !
         ENDIF ! scatread
      ENDDO  ! end loop over q points
      !
      IF (iterative_bte) iq_restart = 1 
      ! If we do not do iterative BTE, then exist the while loop.
      IF (.NOT. iterative_bte) error_h = 0.0_DP
      IF (.NOT. iterative_bte) error_el = 0.0_DP
      ! 
    ENDDO ! End the while loop  
    IF (iterative_bte) DEALLOCATE (Fi_all)
    IF (iterative_bte) DEALLOCATE (F_current)
    IF (iterative_bte) DEALLOCATE (F_SERTA)
    IF (iterative_bte) DEALLOCATE (inv_tau_all)
    IF (iterative_bte) DEALLOCATE (zi_allvb)
    IF (iterative_bte) DEALLOCATE (s_BZtoIBZ_full)
    IF (iterative_bte) DEALLOCATE (ixkqf_tr)
    IF (int_mob .AND. carrier .AND. iterative_bte) DEALLOCATE (inv_tau_allcb)
    IF (int_mob .AND. carrier .AND. iterative_bte) DEALLOCATE (zi_allcb)
    ! 
  ENDIF ! end parallel_k
  ! 
  IF (parallel_q) THEN
    !      
    ! get the size of the matrix elements stored in each pool
    ! for informational purposes.  Not necessary
    !
    CALL mem_size(ibndmin, ibndmax, nmodes, nkf)
    !
    ALLOCATE ( epf17 (ibndmax-ibndmin+1, ibndmax-ibndmin+1, nmodes, nqf) )
    !
    DO ik = 1, nkf
      ! 
      CALL start_clock ( 'ep-interp' )
      !
      DO iq = 1, nqf
        !
        ikq = 2 * iq
        ikk = ikq - 1 
        !   
        ! xkf is assumed to be in crys coord
        !
        xkk = xkf (:, 2*ik-1)
        xxq = xqf (:, iq)
        xkq = xkk + xxq
        !
        ! ------------------------------------------------------
        ! dynamical matrix : Wannier -> Bloch
        ! ------------------------------------------------------
        !
        CALL dynwan2bloch &
             ( nmodes, nrr_q, irvec_qq, ndegen_qq, xxq, uf, w2 )
        !
        ! ...then take into account the mass factors and square-root the frequencies...
        !
        DO nu = 1, nmodes
           !
           ! wf are the interpolated eigenfrequencies
           ! (omega on fine grid)
           !
           IF ( w2 (nu) .gt. 0.d0 ) THEN
              wf(nu,iq) =  sqrt(abs( w2 (nu) ))
           ELSE
              wf(nu,iq) = -sqrt(abs( w2 (nu) ))
           ENDIF
           !
           DO mu = 1, nmodes
              na = (mu - 1) / 3 + 1
              uf (mu, nu) = uf (mu, nu) / sqrt(amass(ityp(na)))
           ENDDO
        ENDDO
        ! ------------------------------------------------------------
        ! epmat : Wannier el and Wannier ph -> Wannier el and Bloch ph
        ! --------------------------------------------------------------
        !
        
        CALL ephwan2blochp &
             ( nmodes, xxq, irvec_qq, ndegen_qq, nrr_q, uf, epmatwef, nbndsub, nrr_k )
        !
        !  number of k points with a band on the Fermi surface
        fermicount = 0
        !
        ! this is a loop over k blocks in the pool
        ! (size of the local k-set)
        !
        ! SP: Compute the cfac only once here since the same are use in both hamwan2bloch and dmewan2bloch
        ! + optimize the 2\pi r\cdot k with Blas
        CALL dgemv('t', 3, nrr_k, twopi, irvec_r, 3, xkk, 1, 0.0_DP, rdotk, 1 )
        cfac(:) = exp( ci*rdotk ) / ndegen_kk(:)
        CALL dgemv('t', 3, nrr_k, twopi, irvec_r, 3, xkq, 1, 0.0_DP, rdotk, 1 )
        cfacq(:) = exp( ci*rdotk ) / ndegen_kk(:)
        ! ------------------------------------------------------        
        ! hamiltonian : Wannier -> Bloch
        ! ------------------------------------------------------
        !
        ! Kohn-Sham first, then get the rotation matricies for following interp.
        IF (eig_read) THEN
           CALL hamwan2bloch &
             ( nbndsub, nrr_k, cufkk, etf_ks (:, ikk), chw_ks, cfac)     
           CALL hamwan2bloch &
             ( nbndsub, nrr_k, cufkq, etf_ks (:, ikq), chw_ks, cfacq)
        ENDIF
        !
        CALL hamwan2bloch &
             ( nbndsub, nrr_k, cufkk, etf (:, ikk), chw, cfac)        
        CALL hamwan2bloch &
             ( nbndsub, nrr_k, cufkq, etf (:, ikq), chw, cfacq)
        !
        ! ------------------------------------------------------        
        !  dipole: Wannier -> Bloch
        ! ------------------------------------------------------        
        !
        CALL dmewan2bloch &
             ( nbndsub, nrr_k, cufkk, dmef (:,:,:, ikk), etf(:,ikk), etf_ks(:,ikk), cfac)
        CALL dmewan2bloch &
             ( nbndsub, nrr_k, cufkq, dmef (:,:,:, ikq), etf(:,ikq), etf_ks(:,ikq), cfac)
        !
        ! ------------------------------------------------------        
        !  velocity: Wannier -> Bloch
        ! ------------------------------------------------------        
        !
        IF (vme) THEN
           IF (eig_read) THEN
              CALL vmewan2bloch &
                   ( nbndsub, nrr_k, irvec_kk, ndegen_kk, xkk, cufkk, vmef (:,:,:, ikk), etf(:,ikk), etf_ks(:,ikk), chw_ks)
              CALL vmewan2bloch &
                   ( nbndsub, nrr_k, irvec_kk, ndegen_kk, xkq, cufkq, vmef (:,:,:, ikq), etf(:,ikq), etf_ks(:,ikq), chw_ks)
           ELSE
              CALL vmewan2bloch &
                   ( nbndsub, nrr_k, irvec_kk, ndegen_kk, xkk, cufkk, vmef (:,:,:, ikk), etf(:,ikk), etf_ks(:,ikk), chw)
              CALL vmewan2bloch &
                   ( nbndsub, nrr_k, irvec_kk, ndegen_kk, xkq, cufkq, vmef (:,:,:, ikq), etf(:,ikq), etf_ks(:,ikq), chw)
           ENDIF
        ENDIF
        !
        ! interpolate ONLY when (k,k+q) both have at least one band 
        ! within a Fermi shell of size fsthick 
        !
        IF ( (( minval ( abs(etf (:, ikk) - ef) ) .lt. fsthick ) .and. &
             ( minval ( abs(etf (:, ikq) - ef) ) .lt. fsthick ))) THEN
           !
           !  fermicount = fermicount + 1
           !
           ! --------------------------------------------------------------
           ! epmat : Wannier el and Bloch ph -> Bloch el and Bloch ph
           ! --------------------------------------------------------------
           !
           ! SP: Note: In case of polar materials, computing the long-range and short-range term 
           !     separately might help speed up the convergence. Indeed the long-range term should be 
           !     much faster to compute. Note however that the short-range term still contains a linear
           !     long-range part and therefore could still be a bit more difficult to converge than 
           !     non-polar materials. 
           ! 
           IF (longrange) THEN
             !      
             epmatf = czero
             !
           ELSE
             !
             CALL ephwan2bloch &
               ( nbndsub, nrr_k, irvec_kk, ndegen_kk, epmatwef, xkk, cufkk, cufkq, epmatf, nmodes )
             !
           ENDIF
           ! 
           IF (lpolar) THEN
             !
             CALL compute_umn_f( nbndsub, cufkk, cufkq, bmatf )
             !
             IF ( (abs(xxq(1)).gt.eps) .or. (abs(xxq(2)).gt.eps) .or. (abs(xxq(3)).gt.eps) ) THEN
                CALL cryst_to_cart (1, xxq, bg, 1)
                DO ibnd = 1, nbndsub
                  DO jbnd = 1, nbndsub
                    CALL rgd_blk_epw(nq1, nq2, nq3, xxq, uf, epmatf(ibnd,jbnd,:), &
                          nmodes, epsi, zstar, bmatf(ibnd,jbnd), +1.d0)
                  ENDDO
                ENDDO
                CALL cryst_to_cart (1, xxq, at, -1)
             ENDIF
             !
           ENDIF
           !            
           ! write epmatf to file / store in memory
           !
           !
           DO jbnd = ibndmin, ibndmax
             DO ibnd = ibndmin, ibndmax
               ! 
               epf17(ibnd-ibndmin+1,jbnd-ibndmin+1,:,iq) = epmatf(ibnd,jbnd,:)
               !
             ENDDO
           ENDDO
           ! 
           !if (ik==2) then
           !  do imode = 1, nmodes
           !    write(*,*) 'epmatf ',SUM((REAL(REAL(epmatf(:,:,imode))))**2)+SUM((REAL(AIMAG(epmatf(:,:,imode))))**2)
           !  enddo
           !endif
           !
        ENDIF
        !
      ENDDO  ! end loop over q points
      !
      IF (phonselfen ) CALL selfen_phon_k( ik )
      IF (elecselfen ) CALL selfen_elec_k( ik )
      IF (nest_fn    ) CALL nesting_fn_k( ik )
      IF (specfun_el ) CALL spectral_func_k( ik )
      IF (ephwrite) THEN
        CALL errore('ephwann_shuffle', "parallel q not (yet) implemented with ephwrite = .true.",1)
      ENDIF
      !
      CALL stop_clock ( 'ep-interp' )
      !
    ENDDO  ! end loop over k points
    ! 
  ENDIF ! end parallel_q
  !
  !  Close th epmatwp file
#if defined(__MPI)
  IF (etf_mem == 1) then
    CALL MPI_FILE_CLOSE(iunepmatwp2,ierr)
    IF( ierr /= 0 ) CALL errore( 'ephwann_shuffle', 'error in MPI_FILE_CLOSE',1 )
  ENDIF
#endif 
  ! 
  ! Check Memory usage
  CALL system_mem_usage(valueRSS)
  ! 
  WRITE(stdout, '(a)' )             '     ==================================================================='
  WRITE(stdout, '(a,i10,a)' ) '     Memory usage:  VmHWM =',valueRSS(2)/1024,'Mb'
  WRITE(stdout, '(a,i10,a)' ) '                   VmPeak =',valueRSS(1)/1024,'Mb'
  WRITE(stdout, '(a)' )             '     ==================================================================='
  WRITE(stdout, '(a)' )
  ! 
  ! ---------------------------------------------------------------------------------------
  ! ---------------------------------------------------------------------------------------  
  !
  ! SP: Added lambda and phonon lifetime writing to file.
  ! 
  CALL mp_barrier(inter_pool_comm)
  IF (mpime.eq.ionode_id) THEN
    !
    IF (phonselfen .and. parallel_k ) THEN
      OPEN(unit=lambda_phself,file='lambda.phself')
      WRITE(lambda_phself, '(/2x,a/)') '#Lambda phonon self-energy'
      WRITE(lambda_phself, *) '#Modes     ',(imode, imode=1,nmodes)
      DO iq = 1, nqtotf
          !
          !myfmt = "(*(3x,E15.5))"  This does not work with PGI
        myfmt = "(1000(3x,E15.5))"
        WRITE(lambda_phself,'(i9,4x)',advance='no') iq
        WRITE(lambda_phself, fmt=myfmt) (REAL(lambda_all(imode,iq,1)),imode=1,nmodes)
          !
      ENDDO
      CLOSE(lambda_phself)
      OPEN(unit=linewidth_phself,file='linewidth.phself')
      WRITE(linewidth_phself, '(a)') '# Phonon frequency and phonon lifetime in meV '
      WRITE(linewidth_phself,'(a)') '# Q-point  Mode   Phonon freq (meV)   Phonon linewidth (meV)'
      DO iq = 1, nqtotf
        !
        DO imode=1, nmodes
          WRITE(linewidth_phself,'(i9,i6,E20.8,E22.10)') iq,imode,&
                                 ryd2mev*wf(imode,iq),ryd2mev*REAL(gamma_all(imode,iq,1))
        ENDDO
        !
      ENDDO
      CLOSE(linewidth_phself)
    ENDIF
  ENDIF
  IF (band_plot) CALL plot_band
  !
  IF (a2f) CALL eliashberg_a2f
  ! 
  ! if scattering is read then Fermi level and scissor have not been computed.
  IF (scatread) THEN
    IF (ABS(scissor) > 0.000001) THEN
      icbm = FLOOR(nelec/2.0d0) + nbndskip + 1
      DO ik = 1, nkf
        ikk = 2 * ik - 1
        ikq = ikk + 1
        DO ibnd = icbm, nbndsub
          etf (ibnd, ikk) = etf (ibnd, ikk) + scissor
          etf (ibnd, ikq) = etf (ibnd, ikq) + scissor
        ENDDO
      ENDDO
      WRITE( stdout, '(5x,"Applying a scissor shift of ",f9.5," eV to the conduction states")' ) scissor * ryd2ev
    ENDIF          
    DO itemp = 1, nstemp
      etemp = transp_temp(itemp)      
      IF (int_mob .OR. carrier) THEN
        ! SP: Determination of the Fermi level for intrinsic or doped carrier 
        !     One also need to apply scissor before calling it.
        ef0(itemp) = fermicarrier( etemp )
      ELSE
        IF ( efermi_read ) THEN
          ef0(itemp) = fermi_energy
        ELSE !SP: This is added for efficiency reason because the efermig routine is slow
          ef0(itemp) = efnew
        ENDIF
      ENDIF
    ENDDO ! itemp
  ENDIF ! if scattering 
  ! 
  IF ( ALLOCATED(lambda_all) )   DEALLOCATE( lambda_all )
  IF ( ALLOCATED(gamma_all) )    DEALLOCATE( gamma_all )
  IF ( ALLOCATED(sigmai_all) )   DEALLOCATE( sigmai_all )
  IF ( ALLOCATED(sigmai_mode) )  DEALLOCATE( sigmai_mode )
  IF ( ALLOCATED(w2) )           DEALLOCATE( w2 )
  DEALLOCATE(cfac)
  DEALLOCATE(cfacq)
  DEALLOCATE(rdotk)
  DEALLOCATE(irvec_r)
  DEALLOCATE(irvec_kk)
  DEALLOCATE(irvec_qq)
  ! 
  IF (.not. iterative_bte) CALL transport_coeffs (ef0,efcb)
  !
  IF ( ALLOCATED(inv_tau_all) )     DEALLOCATE( inv_tau_all )
  IF ( ALLOCATED(inv_tau_allcb) )   DEALLOCATE( inv_tau_allcb )
  IF ( ALLOCATED(zi_allvb) )        DEALLOCATE( zi_allvb )
  IF ( ALLOCATED(zi_allcb) )        DEALLOCATE( zi_allcb )
  !
  CALL stop_clock ( 'ephwann' )
  !
  END SUBROUTINE ephwann_shuffle
  !
  !-------------------------------------------
  SUBROUTINE epw_write
  !-------------------------------------------
  !
  USE epwcom,    ONLY : nbndsub, vme, eig_read, etf_mem
  USE pwcom,     ONLY : ef, nelec, isk
  USE elph2,     ONLY : nrr_k, nrr_q, chw, rdw, cdmew, cvmew, chw_ks, &
                        zstar, epsi, epmatwp
  USE ions_base, ONLY : amass, ityp, nat, tau
  USE cell_base, ONLY : at, bg, omega, alat
  USE phcom,     ONLY : nmodes  
  USE io_epw,    ONLY : epwdata, iundmedata, iunvmedata, iunksdata, iunepmatwp, &
                        crystal
  USE noncollin_module, ONLY : noncolin              
  USE io_files,  ONLY : prefix, diropn
  USE mp,        ONLY : mp_barrier
  USE mp_global, ONLY : inter_pool_comm
  USE mp_world,  ONLY : mpime
  USE io_global, ONLY : ionode_id
  !
  implicit none
  LOGICAL             :: exst
  INTEGER             :: ibnd, jbnd, jmode, imode, irk, irq, ipol, lrepmatw
  CHARACTER (len=256) :: filint
  !
  WRITE(6,'(/5x,"Writing Hamiltonian, Dynamical matrix and EP vertex in Wann rep to file"/)')
  !
  IF (mpime.eq.ionode_id) THEN
    !
    OPEN(unit=epwdata,file='epwdata.fmt')
    OPEN(unit=crystal,file='crystal.fmt')
    OPEN(unit=iundmedata,file='dmedata.fmt')
    IF (vme) OPEN(unit=iunvmedata,file='vmedata.fmt')
    IF (eig_read) OPEN(unit=iunksdata,file='ksdata.fmt')
    WRITE (crystal,*) nat
    WRITE (crystal,*) nmodes
    WRITE (crystal,*) nelec
    WRITE (crystal,*) at
    WRITE (crystal,*) bg
    WRITE (crystal,*) omega
    WRITE (crystal,*) alat
    WRITE (crystal,*) tau
    WRITE (crystal,*) amass
    WRITE (crystal,*) ityp
    WRITE (crystal,*) isk
    WRITE (crystal,*) noncolin
    WRITE (epwdata,*) ef
    WRITE (epwdata,*) nbndsub, nrr_k, nmodes, nrr_q
    WRITE (epwdata,*) zstar, epsi
    !
    DO ibnd = 1, nbndsub
      DO jbnd = 1, nbndsub
        DO irk = 1, nrr_k
          WRITE (epwdata,*) chw(ibnd,jbnd,irk)
          IF (eig_read) WRITE (iunksdata,*) chw_ks(ibnd,jbnd,irk)
          DO ipol = 1,3
            WRITE (iundmedata,*) cdmew(ipol, ibnd,jbnd,irk)
            IF (vme) WRITE (iunvmedata,*) cvmew(ipol, ibnd,jbnd,irk)
          ENDDO
        ENDDO
      ENDDO
    ENDDO
    !
    DO imode = 1, nmodes
      DO jmode = 1, nmodes
        DO irq = 1, nrr_q
          WRITE (epwdata,*) rdw(imode,jmode,irq) 
        ENDDO
      ENDDO
    ENDDO
    !
    IF (etf_mem == 0) THEN
      ! SP: The call to epmatwp is now inside the loop
      !     This is important as otherwise the lrepmatw integer 
      !     could become too large for integer(kind=4).
      !     Note that in Fortran the record length has to be a integer
      !     of kind 4. 
      lrepmatw   = 2 * nbndsub * nbndsub * nrr_k * nmodes
      filint    = trim(prefix)//'.epmatwp'
      CALL diropn (iunepmatwp, 'epmatwp', lrepmatw, exst)
      DO irq = 1, nrr_q
        CALL davcio ( epmatwp(:,:,:,:,irq), lrepmatw, iunepmatwp, irq, +1 )
      ENDDO
      ! 
      CLOSE(iunepmatwp)
      IF (ALLOCATED(epmatwp)) DEALLOCATE ( epmatwp )
    ENDIF 
    !
    CLOSE(epwdata)
    CLOSE(crystal)
    CLOSE(iundmedata)
    IF (vme) CLOSE(iunvmedata)
    IF (eig_read) CLOSE(iunksdata)
    !
  ENDIF
  CALL mp_barrier(inter_pool_comm)
  !---------------------------------
  END SUBROUTINE epw_write
  !---------------------------------
  !---------------------------------
  SUBROUTINE epw_read()
  !---------------------------------
  USE epwcom,    ONLY : nbndsub, vme, eig_read, etf_mem, lifc
  USE pwcom,     ONLY : ef
  USE elph2,     ONLY : nrr_k, nrr_q, chw, rdw, epmatwp, &
                        cdmew, cvmew, chw_ks, zstar, epsi
  USE ions_base, ONLY : nat
  USE phcom,     ONLY : nmodes  
  USE io_global, ONLY : stdout
  USE io_files,  ONLY : prefix, diropn
  USE io_epw,    ONLY : epwdata, iundmedata, iunvmedata, iunksdata, iunepmatwp
  USE constants_epw, ONLY :  czero
#if defined(__NAG)
  USE f90_unix_io,ONLY : flush
#endif
  USE io_global, ONLY : ionode_id
  USE mp,        ONLY : mp_barrier, mp_bcast
  USE mp_global, ONLY : intra_pool_comm, inter_pool_comm, root_pool
  USE mp_world,  ONLY : mpime
  !
  implicit none
  !
  LOGICAL             :: exst
  CHARACTER (len=256) :: filint
  INTEGER             :: ibnd, jbnd, jmode, imode, irk, irq, &
                         ipol, ios, lrepmatw
  !
  WRITE(stdout,'(/5x,"Reading Hamiltonian, Dynamical matrix and EP vertex in Wann rep from file"/)')
  call flush(6)
  ! 
  ! This is important in restart mode as zstar etc has not been allocated
  IF (.NOT. ALLOCATED (zstar) ) ALLOCATE( zstar(3,3,nat) )
  IF (.NOT. ALLOCATED (epsi) ) ALLOCATE( epsi(3,3) )

  IF (mpime.eq.ionode_id) THEN
    !
    OPEN(unit=epwdata,file='epwdata.fmt',status='old',iostat=ios)
    IF (ios /= 0) call errore ('ephwann_shuffle', 'error opening epwdata.fmt',epwdata)
    OPEN(unit=iundmedata,file='dmedata.fmt',status='old',iostat=ios)
    IF (ios /= 0) call errore ('ephwann_shuffle', 'error opening dmedata.fmt',iundmedata)
    IF (eig_read) OPEN(unit=iunksdata,file='ksdata.fmt',status='old',iostat=ios)
    IF (eig_read .AND. ios /= 0) call errore ('ephwann_shuffle', 'error opening ksdata.fmt',iunksdata)
    IF (vme) OPEN(unit=iunvmedata,file='vmedata.fmt',status='old',iostat=ios)
    IF (vme .AND. ios /= 0) call errore ('ephwann_shuffle', 'error opening vmedata.fmt',iunvmedata)
    READ (epwdata,*) ef
    READ (epwdata,*) nbndsub, nrr_k, nmodes, nrr_q
    READ (epwdata,*) zstar, epsi
    ! 
  ENDIF
  CALL mp_bcast (ef, ionode_id, inter_pool_comm)
  CALL mp_bcast (ef, root_pool, intra_pool_comm)
  !
  CALL mp_bcast (nbndsub, ionode_id, inter_pool_comm)
  CALL mp_bcast (nbndsub, root_pool, intra_pool_comm)
  !
  CALL mp_bcast (nrr_k, ionode_id, inter_pool_comm)
  CALL mp_bcast (nrr_k, root_pool, intra_pool_comm)
  !
  CALL mp_bcast (nmodes, ionode_id, inter_pool_comm)
  CALL mp_bcast (nmodes, root_pool, intra_pool_comm)
  !
  CALL mp_bcast (nrr_q, ionode_id, inter_pool_comm)
  CALL mp_bcast (nrr_q, root_pool, intra_pool_comm)
  !
  CALL mp_bcast (zstar, ionode_id, inter_pool_comm)
  CALL mp_bcast (zstar, root_pool, intra_pool_comm)
  !
  CALL mp_bcast (epsi, ionode_id, inter_pool_comm)
  CALL mp_bcast (epsi, root_pool, intra_pool_comm)
  !
  IF (.not. ALLOCATED(chw)    ) ALLOCATE ( chw ( nbndsub, nbndsub, nrr_k )            )
  IF (.not. ALLOCATED(chw_ks) ) ALLOCATE ( chw_ks ( nbndsub, nbndsub, nrr_k )         )
  IF (.not. ALLOCATED(rdw)    ) ALLOCATE ( rdw ( nmodes,  nmodes,  nrr_q )            )
  IF (.not. ALLOCATED(cdmew)  ) ALLOCATE ( cdmew ( 3, nbndsub, nbndsub, nrr_k )       )
  IF (vme .and. (.not.ALLOCATED(cvmew))  ) ALLOCATE ( cvmew   ( 3, nbndsub, nbndsub, nrr_k )     )
  !
  IF (mpime.eq.ionode_id) THEN
    !
    DO ibnd = 1, nbndsub
       DO jbnd = 1, nbndsub
          DO irk = 1, nrr_k
             READ (epwdata,*) chw(ibnd,jbnd,irk)
             IF (eig_read) READ (iunksdata,*) chw_ks(ibnd,jbnd,irk)
             DO ipol = 1,3
                READ (iundmedata,*) cdmew(ipol, ibnd,jbnd,irk)
                IF (vme) READ (iunvmedata,*) cvmew(ipol, ibnd,jbnd,irk)
             ENDDO
          ENDDO
       ENDDO
    ENDDO
    !
    IF (.not. lifc) THEN
       DO imode = 1, nmodes
          DO jmode = 1, nmodes
             DO irq = 1, nrr_q
                READ (epwdata,*) rdw(imode,jmode,irq)
             ENDDO
          ENDDO
       ENDDO
    ENDIF
    !
  ENDIF
  !
  CALL mp_bcast (chw, ionode_id, inter_pool_comm)
  CALL mp_bcast (chw, root_pool, intra_pool_comm)
  !
  IF (eig_read) CALL mp_bcast (chw_ks, ionode_id, inter_pool_comm)
  IF (eig_read) CALL mp_bcast (chw_ks, root_pool, intra_pool_comm)
  !
  IF (.not. lifc) CALL mp_bcast (rdw, ionode_id, inter_pool_comm)
  IF (.not. lifc) CALL mp_bcast (rdw, root_pool, intra_pool_comm)
  !
  CALL mp_bcast (cdmew, ionode_id, inter_pool_comm)
  CALL mp_bcast (cdmew, root_pool, intra_pool_comm)
  !
  IF (vme) CALL mp_bcast (cvmew, ionode_id, inter_pool_comm)
  IF (vme) CALL mp_bcast (cvmew, root_pool, intra_pool_comm)
  
  IF (lifc) CALL read_ifc
  !
  IF (etf_mem == 0) then
    IF (.not. ALLOCATED(epmatwp)) ALLOCATE ( epmatwp ( nbndsub, nbndsub, nrr_k, nmodes, nrr_q) )
    epmatwp = czero
    IF (mpime.eq.ionode_id) THEN
      ! SP: The call to epmatwp is now inside the loop
      !     This is important as otherwise the lrepmatw integer 
      !     could become too large for integer(kind=4).
      !     Note that in Fortran the record length has to be a integer
      !     of kind 4.      
      lrepmatw   = 2 * nbndsub * nbndsub * nrr_k * nmodes
      filint    = trim(prefix)//'.epmatwp'
      CALL diropn (iunepmatwp, 'epmatwp', lrepmatw, exst)
      DO irq = 1, nrr_q
        CALL davcio ( epmatwp(:,:,:,:,irq), lrepmatw, iunepmatwp, irq, -1 )
      ENDDO
      !  
      CLOSE(iunepmatwp)
    ENDIF
    !
    CALL mp_bcast (epmatwp, ionode_id, inter_pool_comm)
    CALL mp_bcast (epmatwp, root_pool, intra_pool_comm)
    !
  ENDIF
  !
  CALL mp_barrier(inter_pool_comm)
  IF (mpime.eq.ionode_id) THEN
    CLOSE(epwdata)
    CLOSE(iundmedata)
    IF (vme) CLOSE(iunvmedata)
  ENDIF
  !
  WRITE(stdout,'(/5x,"Finished reading Wann rep data from file"/)')
  !
  !---------------------------------
  END SUBROUTINE epw_read
  !---------------------------------
  !---------------------------------
  SUBROUTINE mem_size(ibndmin, ibndmax, nmodes, nkf) 
  !---------------------------------
  !!
  !!  SUBROUTINE estimates the amount of memory taken up by 
  !!  the $$<k+q| dV_q,nu |k>$$ on the fine meshes and prints 
  !!  out a useful(?) message   
  !!
  USE io_global,     ONLY : stdout
  USE kinds,         ONLY : DP
  !
  implicit none
  !
  INTEGER, INTENT (in) :: ibndmin
  !! Min band
  INTEGER, INTENT (in) :: ibndmax
  !! Min band
  INTEGER, INTENT (in) :: nmodes
  !! Number of modes
  INTEGER, INTENT (in) :: nkf
  !! Number of k-points in pool 
  !
  ! Work variables
  INTEGER             :: imelt
  REAL(kind=DP)       :: rmelt
  CHARACTER (len=256) :: chunit
  !
  imelt = (ibndmax-ibndmin+1)**2 * nmodes * nkf
  rmelt = imelt * 8 / 1048576.d0 ! 8 bytes per number, value in Mb
  IF (rmelt .lt. 1000.0 ) THEN
     chunit =  ' Mb '
     IF (rmelt .lt. 1.0 ) THEN
        chunit = ' Kb '
        rmelt  = rmelt * 1024.d0
     ENDIF
  ELSE
     rmelt = rmelt / 1024.d0
     chunit = ' Gb '
  ENDIF
  WRITE(stdout,'(/,5x,a, i13, a,f7.2,a,a)') "Number of ep-matrix elements per pool :", &
       imelt, " ~= ", rmelt, trim(chunit), " (@ 8 bytes/ DP)"
  !
  !---------------------------------
  END SUBROUTINE mem_size
  !---------------------------------
  ! 
  !--------------------------------------------------------------------
  FUNCTION efermig_seq (et, nbnd, nks, nelec, wk, Degauss, Ngauss, is, isk)
  !--------------------------------------------------------------------
  !!
  !!     Finds the Fermi energy - Gaussian Broadening
  !!     (see Methfessel and Paxton, PRB 40, 3616 (1989 )
  !!
  USE io_global, ONLY : stdout
  USE kinds,     ONLY : DP
  USE constants, ONLY : rytoev
  !
  implicit none
  !
  INTEGER, INTENT (in) :: nks
  !! Number of k-points per pool
  INTEGER, INTENT (in) :: nbnd
  !! Number of band
  INTEGER, INTENT (in) :: Ngauss
  !! 
  INTEGER, INTENT (in) :: is
  !! 
  INTEGER, INTENT (in) :: isk(nks)
  !! 
  ! 
  REAL (kind=DP), INTENT (in) :: wk (nks)
  !!
  REAL (kind=DP), INTENT (in) :: et (nbnd, nks)
  !!
  REAL (kind=DP), INTENT (in) :: Degauss
  !!
  REAL (kind=DP), INTENT (in) :: nelec
  !! Number of electron (charge)
  ! 
  real(DP) :: efermig_seq
  !
  ! Local variables
  ! 
  real(DP), parameter :: eps= 1.0d-10
  integer, parameter :: maxiter = 300
  real(DP) :: Ef, Eup, Elw, sumkup, sumklw, sumkmid
  real(DP), external::  sumkg_seq
  integer :: i, kpoint
  !
  !  find bounds for the Fermi energy. Very safe choice!
  !
  Elw = et (1, 1)
  Eup = et (nbnd, 1)
  do kpoint = 2, nks
     Elw = min (Elw, et (1, kpoint) )
     Eup = max (Eup, et (nbnd, kpoint) )
  enddo
  Eup = Eup + 2 * Degauss
  Elw = Elw - 2 * Degauss
  !
  !  Bisection method
  !
  sumkup = sumkg_seq (et, nbnd, nks, wk, Degauss, Ngauss, Eup, is, isk)
  sumklw = sumkg_seq (et, nbnd, nks, wk, Degauss, Ngauss, Elw, is, isk)
  if ( (sumkup - nelec) < -eps .or. (sumklw - nelec) > eps )  &
       call errore ('efermig_seq', 'internal error, cannot bracket Ef', 1)
  DO i = 1, maxiter
    Ef = (Eup + Elw) / 2.d0
    sumkmid = sumkg_seq (et, nbnd, nks, wk, Degauss, Ngauss, Ef, is, isk)
    if (abs (sumkmid-nelec) < eps) then
       efermig_seq = Ef
       return
    elseif ( (sumkmid-nelec) < -eps) then
       Elw = Ef
    else
       Eup = Ef
    endif
  ENDDO
  IF (is /= 0) WRITE(stdout, '(5x,"Spin Component #",i3)') is
  WRITE( stdout, '(5x,"Warning: too many iterations in bisection"/ &
       &      5x,"Ef = ",f10.6," sumk = ",f10.6," electrons")' ) &
       Ef * rytoev, sumkmid
  !
  efermig_seq = Ef
  RETURN
  !
  end FUNCTION efermig_seq
  !
  !-----------------------------------------------------------------------
  function sumkg_seq (et, nbnd, nks, wk, degauss, ngauss, e, is, isk)
  !-----------------------------------------------------------------------
  !!
  !!  This function computes the number of states under a given energy e
  !!
  USE kinds, ONLY : DP
  USE mp,    ONLY : mp_sum
  ! 
  implicit none
  ! 
  INTEGER, INTENT (in) :: nks
  !! the total number of K points
  INTEGER, INTENT (in) :: nbnd
  !! the number of bands
  INTEGER, INTENT (in) :: ngauss
  !! the type of smearing
  INTEGER, INTENT (in) :: is
  !!
  INTEGER, INTENT (in) :: isk(nks)
  !!
  !
  REAL(kind=DP), INTENT (in) :: wk (nks)
  !! the weight of the k points
  REAL(kind=DP), INTENT (in) :: et (nbnd, nks) 
  !! the energy eigenvalues
  REAL(kind=DP), INTENT (in) :: degauss
  !! gaussian broadening
  REAL(kind=DP), INTENT (in) :: e
  !! the energy to check
  !
  REAL(kind=DP)  :: sumkg_seq 
  !! 
  !
  ! local variables
  !
  real(DP), external :: wgauss
  ! function which compute the smearing 
  real(DP) ::sum1
  integer :: ik, ibnd
  ! counter on k points
  ! counter on the band energy
  !
  sumkg_seq = 0.d0
  DO ik = 1, nks
    sum1 = 0.d0
    if (is /= 0) then
       if (isk(ik).ne.is) cycle
    end if
    do ibnd = 1, nbnd
       sum1 = sum1 + wgauss ( (e-et (ibnd, ik) ) / degauss, ngauss)
    enddo
    sumkg_seq = sumkg_seq + wk (ik) * sum1
  ENDDO
  RETURN
  !
  end function sumkg_seq
  !
  !-----------------------------------------------------------------
  subroutine rwepmatw ( epmatw, nbnd, np, nmodes, nrec, iun, iop)
  !-----------------------------------------------------------------
  !!
  !! A simple wrapper to the davcio routine to read/write arrays
  !! instead of vectors 
  !!
  !-----------------------------------------------------------------
  USE kinds, ONLY : DP
  USE mp,    ONLY : mp_barrier
  ! 
  implicit none
  ! 
  INTEGER, INTENT (in) :: nbnd
  !! Total number of bands
  INTEGER, INTENT (in) :: np
  !! np is either nrr_k or nq (epmatwe and epmatwp have the same structure)
  INTEGER, INTENT (in) :: nmodes
  !! Number of modes
  INTEGER, INTENT (in) :: nrec
  !! Place where to start reading/writing
  INTEGER, INTENT (in) :: iun
  !! Record number
  INTEGER, INTENT (in) :: iop
  !! If -1, read and if +1 write the matrix
  ! 
  COMPLEX(kind=DP), intent (inout) :: epmatw(nbnd,nbnd,np,nmodes)
  !! El-ph matrix to read or write
  !
  ! Local variables
  integer :: lrec, i, ibnd, jbnd, imode, ip
  complex(kind=DP):: aux ( nbnd*nbnd*np*nmodes )
  !
  lrec = 2 * nbnd * nbnd * np * nmodes
  !
  IF ( iop .eq. -1 ) then
    !
    !  read matrix
    !
    CALL davcio ( aux, lrec, iun, nrec, -1 )
    !
    i = 0
    DO imode = 1, nmodes
     DO ip = 1, np
      DO jbnd = 1, nbnd
       DO ibnd = 1, nbnd
         i = i + 1
         epmatw ( ibnd, jbnd, ip, imode ) = aux (i)
         ! 
       ENDDO
      ENDDO
     ENDDO
    ENDDO
    !
  ELSEif ( iop .eq. 1 ) then 
    !
    !  write matrix
    !
    i = 0
    DO imode = 1, nmodes
     DO ip = 1, np
      DO jbnd = 1, nbnd
       DO ibnd = 1, nbnd
         i = i + 1
         aux (i) = epmatw ( ibnd, jbnd, ip, imode )
       ENDDO
      ENDDO
     ENDDO
    ENDDO
    !
    CALL davcio ( aux, lrec, iun, nrec, +1 )
    !
  ELSE
    !
    CALL errore ('rwepmatw','iop not permitted',1)
    !
  ENDIF
  !
  ! ----------------------------------------------------------------------
  end subroutine rwepmatw
  ! ----------------------------------------------------------------------
  ! 
  !-----------------------------------------------------------------------
  FUNCTION fermicarrier( temp )
  !-----------------------------------------------------------------------
  !!
  !!  This subroutine computes the Fermi energy associated with a given 
  !!  carrier concentration using bissection
  !!
  !-----------------------------------------------------------------------
  USE cell_base, ONLY : omega, alat, at
  USE kinds,     ONLY : DP
  USE io_global, ONLY : stdout
  USE elph2,     ONLY : etf, nkf, wkf
  USE constants_epw, ONLY : ryd2ev, bohr2ang, ang2cm
  USE noncollin_module, ONLY : noncolin
  USE pwcom,     ONLY : nelec
  USE epwcom,    ONLY : int_mob, nbndsub, ncarrier, &
                        system_2d, carrier
  USE mp,        ONLY : mp_barrier, mp_sum, mp_max, mp_min
  USE mp_global, ONLY : inter_pool_comm
  !
  IMPLICIT NONE
  !
  REAL(KIND=DP), intent(in) :: temp
  !! Temperature in kBT [Ry] unit.
  INTEGER :: i
  !! Index for the bisection iteration
  INTEGER :: ik
  !! k-point index per pool
  INTEGER :: ikk
  !! Odd index to read etf
  INTEGER :: ibnd
  !! Local band index
  INTEGER :: ivbm
  !! Index of the VBM
  INTEGER :: icbm
  !! Index of the CBM
  INTEGER, PARAMETER :: maxiter = 500 ! 300
  !! Maximum interation
  REAL(KIND=DP) :: fermicarrier
  !! Fermi level returned
  REAL(KIND=DP) :: fnk
  !! Fermi-Diract occupation
  REAL(KIND=DP) :: ks_exp(nbndsub, nkf)
  !! Exponential of the eigenvalues divided by kBT
  REAL(KIND=DP) :: fermi_exp
  !! Fermi level in exponential format
  REAL(KIND=DP) :: rel_err
  !! Relative error
  REAL(KIND=DP) :: factor
  !! Factor that goes from number of carrier per unit cell to number of
  !! carrier per cm^-3
  REAL(KIND=DP) :: arg
  !! Argument of the exponential
  REAL(KIND=DP) :: inv_cell
  !! Inverse of the volume in [Bohr^{-3}]
  REAL(KIND=DP) :: evbm
  !! Energy of the VBM
  REAL(KIND=DP) :: ecbm
  !! Energy of the CBM
  REAL(KIND=DP) :: Ef
  !! Energy of the current Fermi level for the bisection method
  REAL(KIND=DP) :: Elw
  !! Energy lower bound for the bisection method
  REAL(KIND=DP) :: Eup
  !! Energy upper bound for the bisection method
  REAL(KIND=DP) :: hole_density
  !! Hole carrier density
  REAL(KIND=DP) :: electron_density 
  !! Electron carrier density
  REAL(DP), PARAMETER :: maxarg = 200.d0
  !! Maximum value for the argument of the exponential
  REAL(DP), PARAMETER :: eps= 1.0d-5
  !! Tolerence to be converged [relative]
  ! 
  Ef = 0.0d0
  ! 
  inv_cell = 1.0d0/omega
  ! for 2d system need to divide by area (vacuum in z-direction)
  IF ( system_2d ) &
     inv_cell = inv_cell * at(3,3) * alat

  ! vbm index
  IF ( noncolin ) THEN
    ivbm = FLOOR(nelec/1.0d0)
  ELSE
    ivbm = FLOOR(nelec/2.0d0)
  ENDIF  
  icbm = ivbm + 1 ! Nb of bands
  !
  ! Initialization value. Should be large enough ...
  evbm = -10000
  ecbm = 10000 ! In Ry
  ! 
  DO ik = 1, nkf
    ikk = 2 * ik - 1
    DO ibnd = 1, nbndsub
      IF ( ibnd < ivbm+1) THEN
        IF (etf (ibnd, ikk) > evbm ) THEN
          evbm = etf (ibnd, ikk)
        ENDIF
      ENDIF
      ! Find cbm index 
      IF (ibnd > ivbm) THEN
        IF (etf (ibnd, ikk) < ecbm ) THEN
          ecbm = etf (ibnd, ikk)
        ENDIF
      ENDIF
      ! 
    ENDDO
  ENDDO
  !
  ! Find max and min across pools
  !         
  call mp_max( evbm, inter_pool_comm )
  call mp_min( ecbm, inter_pool_comm )
  !    
  WRITE( stdout, '(5x,"Valence band maximum    = ",f10.6," eV")' )  evbm * ryd2ev
  WRITE( stdout, '(5x,"Conduction band minimum = ",f10.6," eV")' )  ecbm * ryd2ev


  !WRITE(stdout,*),'evbm ',evbm * ryd2ev
  !WRITE(stdout,*),'ecbm ',ecbm * ryd2ev
  !WRITE(stdout,*),'temp ',temp * ryd2ev
  ! 
  ! Store e^(e_nk/kbT) on each core
  DO ik = 1, nkf
    DO ibnd = 1, nbndsub
      ikk = 2 * ik - 1
      ! Because the number are so large. It does lead to instabilities
      ! Therefore we rescale everything to the VBM
      arg = (etf (ibnd, ikk) - evbm )/ temp 
      !
      IF (arg .lt. - maxarg) THEN
        ks_exp(ibnd,ik) = 0.0
      ELSE
        ks_exp(ibnd,ik) = exp (arg)
      ENDIF
    ENDDO
  ENDDO
  ! 
  ! Starting bounds energy for the biscection method. The energies are rescaled
  ! to the VBM
  Elw = 1.0d0
  Eup = exp (- (ecbm-evbm) / temp)
  !
  ! Intrinsic mobilities (electron and hole concentration are the same)   
  IF (int_mob .AND. .NOT. carrier) THEN
    ! Use bisection method
    DO i = 1, maxiter
      !
      !WRITE(stdout,*),'Iteration ',i
      ! We want Ef = (Eup + Elw) / 2.d0 but the variables are exp therefore:
      Ef = SQRT(Eup)*SQRT(Elw)     
      ! 
      !WRITE(stdout,*),'Ef ', - log (Ef) * temp * ryd2ev
      hole_density = 0.0
      electron_density = 0.0
      DO ik = 1, nkf
        ikk = 2 * ik - 1
        ! Compute hole carrier concentration
        DO ibnd = 1, ivbm 
          ! Discard very large numbers
          IF (ks_exp(ibnd,ik) * Ef > 10E30) THEN
            fnk = 0.0d0
          ELSE
            fnk = 1.0d0 / ( ks_exp(ibnd,ik) * Ef  + 1.0d0)  
          ENDIF
          ! The wkf(ikk) already include a factor 2
          hole_density = hole_density + wkf(ikk) * (1.0d0 - fnk)
        ENDDO
        ! Compute electron carrier concentration
        DO ibnd = icbm, nbndsub            
          ! Discard very large numbers
          IF (ks_exp(ibnd,ik) * Ef > 10E30) THEN
            fnk = 0.0d0
          ELSE
            fnk = 1.0d0 / ( ks_exp(ibnd,ik) * Ef  + 1.0d0)
          ENDIF
          ! The wkf(ikk) already include a factor 2
          electron_density = electron_density + wkf(ikk) * fnk
        ENDDO    
        ! 
      ENDDO 
      !
      CALL mp_sum( hole_density, inter_pool_comm )
      CALL mp_sum( electron_density, inter_pool_comm )
      ! 
!      WRITE(stdout,*),'hole_density ',hole_density * (1.0d0/omega) * ( bohr2ang * ang2cm  )**(-3)
!      WRITE(stdout,*),'electron_density ',electron_density * (1.0d0/omega) * (bohr2ang * ang2cm  )**(-3)
      rel_err = (hole_density-electron_density) / hole_density
      !
      IF (abs (rel_err) < eps) THEN
        fermi_exp = Ef
        fermicarrier = (- log (fermi_exp) * temp) + evbm
        return
      ELSEIF( (rel_err) > eps) THEN
        Elw = Ef                           
      ELSE                                   
        Eup = Ef 
      ENDIF
    ENDDO ! iteration
  ENDIF 
  ! 
  factor = inv_cell * ( bohr2ang * ang2cm  )**(-3)
  ! Electron doped mobilities (Carrier concentration should be larger than 1E5 cm^-3)   
  IF (ncarrier > 1E5) THEN
    ! Use bisection method
    DO i = 1, maxiter
      ! We want Ef = (Eup + Elw) / 2.d0 but the variables are exp therefore:
      Ef = SQRT(Eup)*SQRT(Elw)
      ! 
      electron_density = 0.0
      DO ik = 1, nkf
        ikk = 2 * ik - 1
        ! Compute electron carrier concentration
        DO ibnd = icbm, nbndsub
          ! Discard very large numbers
          IF (ks_exp(ibnd,ik) * Ef > 10E30) THEN
            fnk = 0.0d0
          ELSE
            fnk = 1.0d0 / ( ks_exp(ibnd,ik) * Ef  + 1.0d0)
          ENDIF
          ! The wkf(ikk) already include a factor 2
          electron_density = electron_density + wkf(ikk) * fnk * factor
        ENDDO
        ! 
      ENDDO
      !
      CALL mp_sum( electron_density, inter_pool_comm )
      ! 
      rel_err = (electron_density-ncarrier) / electron_density
      !
      IF (abs (rel_err) < eps) THEN
        fermi_exp = Ef
        fermicarrier = (- log (fermi_exp) * temp) + evbm
        return
      ELSEIF( (rel_err) > eps) THEN
        Eup = Ef
      ELSE
        Elw = Ef
      ENDIF
    ENDDO ! iteration
  ENDIF
  ! 
  ! Hole doped mobilities (Carrier concentration should be larger than 1E5 cm^-3)   
  IF (ncarrier < -1E5) THEN
    ! Use bisection method
    DO i = 1, maxiter
      ! We want Ef = (Eup + Elw) / 2.d0 but the variables are exp therefore:
      Ef = SQRT(Eup)*SQRT(Elw)
      ! 
      hole_density = 0.0
      DO ik = 1, nkf
        ikk = 2 * ik - 1
        ! Compute hole carrier concentration
        DO ibnd = 1, ivbm
          ! Discard very large numbers
          IF (ks_exp(ibnd,ik) * Ef > 10E30) THEN
            fnk = 0.0d0
          ELSE
            fnk = 1.0d0 / ( ks_exp(ibnd,ik) * Ef  + 1.0d0)
          ENDIF
          ! The wkf(ikk) already include a factor 2
          hole_density = hole_density + wkf(ikk) * (1.0d0 - fnk) * factor
        ENDDO
        ! 
      ENDDO
      !
      CALL mp_sum( hole_density, inter_pool_comm )
      !
      ! In this case ncarrier is a negative number
      rel_err = (hole_density-ABS(ncarrier)) / hole_density
      !
      IF (abs (rel_err) < eps) THEN
        fermi_exp = Ef
        fermicarrier = (- log (fermi_exp) * temp) + evbm
        return
      ELSEIF( (rel_err) > eps) THEN
        Elw = Ef
      ELSE
        Eup = Ef
      ENDIF
    ENDDO ! iteration
  ENDIF
  ! 
  fermi_exp = Ef
  fermicarrier = (- log (fermi_exp) * temp) + evbm
  ! 
  WRITE( stdout, '(5x,"Warning: too many iterations in bisection"/ &
       &      5x,"Ef = ",f10.6)' ) fermicarrier * ryd2ev
  !
  return   
  END FUNCTION fermicarrier
  !--------------------------------------------------------------------------
