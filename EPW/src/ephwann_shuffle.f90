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
  !
  !  Wannier interpolation of electron-phonon vertex
  !
  !  Scalar implementation   Feb 2006
  !  Parallel version        May 2006
  !  Disentenglement         Oct 2006
  !  Compact formalism       Dec 2006
  !  Phonon irreducible zone Mar 2007
  !
  !  No ultrasoft now
  !  No spin polarization
  !
  !  RM - add noncolin case
  !-----------------------------------------------------------------------
  !
  USE kinds,         ONLY : DP
  USE pwcom,         ONLY : nbnd, nks, nkstot, isk, &
                            et, xk, at, bg, ef,  nelec
  USE start_k,       ONLY : nk1, nk2, nk3
  USE ions_base,     ONLY : amass, ityp
  USE phcom,         ONLY : nq1, nq2, nq3, nmodes, w2
  USE control_lr,    ONLY : lgamma
  USE epwcom,        ONLY : nbndsub, lrepmatf, fsthick, lretf, epwread,   &
                            epwwrite, ngaussw, degaussw, lpolar,          &
                            nbndskip, parallel_k, parallel_q, etf_mem,    &
                            elecselfen, phonselfen, nest_fn, a2f, indabs, &
                            epexst, vme, eig_read, ephwrite, mp_mesh_k,   & 
                            efermi_read, fermi_energy, specfun, band_plot
  USE noncollin_module, ONLY : noncolin
  USE constants_epw, ONLY : ryd2ev, ryd2mev, one, two
  USE control_flags, ONLY : iverbosity
  USE io_files,      ONLY : prefix, diropn
  USE io_global,     ONLY : stdout, ionode
  USE io_epw,        ONLY : lambda_phself,linewidth_phself,linewidth_elself, iunepmatf, &
                            iuetf, iunepmatwe, iunepmatwp
  USE elph2,         ONLY : nrr_k, nrr_q, cu, cuq, lwin, lwinq, irvec, ndegen_k, ndegen_q, &
                            wslen, chw, chw_ks, cvmew, cdmew, rdw, epmatwp, epmatq, &
                            wf, etf, etf_k, etf_ks, xqf, xkf, wkf, wqf, &
                            dynq, nqtotf, nkqf, epf17, nkf, nqf, et_ks, &
                            ibndmin, ibndmax, lambda_all, dmec, dmef, vmef, &
                            sigmai_all, sigmai_mode, gamma_all, nkqtotf, epsi, zstar, efnew
#ifdef __NAG
  USE f90_unix_io,   ONLY : flush
  USE,INTRINSIC :: f90_unix_file, ONLY:fstat, stat_t
#endif
#ifdef __PARA
  USE mp,            ONLY : mp_barrier, mp_bcast, mp_sum
  USE io_global,     ONLY : ionode_id
  USE mp_global,     ONLY : my_pool_id, nproc_pool, intra_image_comm, &
                            inter_pool_comm, me_pool, root_pool, intra_pool_comm, &
                            my_pool_id
  USE mp_world,      ONLY : mpime
#endif
  !
  implicit none
  !
#ifdef __NAG
  TYPE(stat_t) :: statb
#endif
#ifndef __NAG
  integer :: fstat,statb(13)
#endif
  !
  complex(kind=DP), ALLOCATABLE :: &
    epmatwe  (:,:,:,:,:),       &! e-p matrix  in wannier basis - electrons
    epmatwe_mem  (:,:,:,:),     &! e-p matrix  in wannier basis - electrons (written on disk)
    epmatwef (:,:,:,:)           ! e-p matrix  in el wannier - fine Bloch phonon grid
  complex(kind=DP), ALLOCATABLE :: &
    epmatf( :, :, :),           &! e-p matrix  in smooth Bloch basis, fine mesh
    cufkk ( :, :),              &! Rotation matrix, fine mesh, points k
    cufkq ( :, :),              &! the same, for points k+q
    uf    ( :, :),              &! Rotation matrix for phonons
    bmatf ( :, :)                ! overlap U_k+q U_k^\dagger in smooth Bloch basis, fine mesh
  integer :: &
    nqc                          ! number of qpoints in the coarse grid
  real(kind=DP) :: &
    xqc (3, nqc)                 ! qpoint list, coarse mesh
  !
  integer :: iq, ik, ikk, ikq, ibnd, jbnd, imode, ir, na, nu, mu, &
    fermicount, nrec, indnew, indold, lrepmatw, ios, irq
  LOGICAL :: already_skipped, exst
  character (len=256) :: filint
  character (len=30)  :: myfmt
  real(kind=DP) :: xxq(3), xxk(3), xkk(3), xkq(3), size_m
  real(kind=DP), external :: efermig
  real(kind=DP), external :: efermig_seq
  real(kind=DP), parameter :: eps = 0.01/ryd2mev
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
  !
  ! determine Wigner-Seitz points
  !
  CALL wigner_seitz2 &
       ( nk1, nk2, nk3, nq1, nq2, nq3, nrr_k, nrr_q, irvec, wslen, ndegen_k, ndegen_q )
  !
  IF ((.NOT. etf_mem) .AND. (ionode)) THEN
    ! open the .epmatwe file with the proper record length
    lrepmatw   = 2 * nbndsub * nbndsub * nrr_k * nmodes
    filint    = trim(prefix)//'.epmatwe'
    CALL diropn (iunepmatwe, 'epmatwe', lrepmatw, exst)  
    filint    = trim(prefix)//'.epmatwp'
    CALL diropn (iunepmatwp, 'epmatwp', lrepmatw, exst)
    !CALL seqopn (iunepmatwp, 'epmatwp', lrepmatw, exst)
  ENDIF
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
     !print*,'nrr_k',nrr_k
     !print*,'nrr_q',nrr_q
     IF (etf_mem) THEN
       ALLOCATE(epmatwe ( nbndsub, nbndsub, nrr_k, nmodes, nqc))
       ALLOCATE (epmatwp ( nbndsub, nbndsub, nrr_k, nmodes, nrr_q))
     ELSE
       ALLOCATE(epmatwe_mem ( nbndsub, nbndsub, nrr_k, nmodes))
     ENDIF
     !
     ! Hamiltonian
     !
     CALL hambloch2wan &
          ( nbnd, nbndsub, nks, nkstot, lgamma, et, xk, cu, lwin, nrr_k, irvec, wslen, chw )
     !
     ! Kohn-Sham eigenvalues
     !
     IF (eig_read) THEN
       WRITE (6,'(5x,a)') "Interpolating MB and KS eigenvalues"
       CALL hambloch2wan &
            ( nbnd, nbndsub, nks, nkstot, lgamma, et_ks, xk, cu, lwin, nrr_k, irvec, wslen, chw_ks )
     ENDIF
     !
     ! Dipole
     !
    ! CALL dmebloch2wan &
    !      ( nbnd, nbndsub, nks, nkstot, nkstot, lgamma, dmec, xk, cu, nrr_k, irvec, wslen )
     CALL dmebloch2wan &
          ( nbnd, nbndsub, nks, nkstot, dmec, xk, cu, nrr_k, irvec, wslen )
     !
     ! Dynamical Matrix 
     !
     CALL dynbloch2wan &
          ( nmodes, nqc, xqc, dynq, nrr_q, irvec, wslen )
     !
     ! Transform of position matrix elements
     ! PRB 74 195118  (2006)
     IF (vme) CALL vmebloch2wan &
         ( nbnd, nbndsub, nks, nks, nkstot, lgamma, xk, cu, nrr_k, irvec, wslen )
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
         IF (etf_mem) THEN 
           CALL ephbloch2wane &
             ( nbnd, nbndsub, nks, nkstot, lgamma, xk, cu, cuq, lwin, lwinq, &
             epmatq (:,:,:,imode,iq), nrr_k, irvec, wslen, epmatwe(:,:,:,imode,iq) )
         ELSE
           CALL ephbloch2wane &
             ( nbnd, nbndsub, nks, nkstot, lgamma, xk, cu, cuq, lwin, lwinq, &
             epmatq (:,:,:,imode,iq), nrr_k, irvec, wslen, epmatwe_mem(:,:,:,imode) )
           !
         ENDIF
         !
       ENDDO
       ! Only the master node writes 
       IF ((.NOT. etf_mem) .AND. (ionode)) THEN
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
       IF (etf_mem) THEN
         CALL ephbloch2wanp &
           ( nbndsub, nmodes, xqc, nqc, irvec, wslen, nrr_k, nrr_q, epmatwe )
       ELSE
          CALL ephbloch2wanp_mem &
           ( nbndsub, nmodes, xqc, nqc, irvec, wslen, nrr_k, nrr_q, epmatwe_mem )
       ENDIF
     ENDIF
     !
#ifdef __PARA
     CALL mp_barrier(inter_pool_comm)
#endif
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
  !
  ! at this point, we will interpolate the Wannier rep to the Bloch rep 
  ! for electrons, phonons and the ep-matrix
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
       bmatf( nbndsub, nbndsub) )

  ! allocate dipole matrix elements after getting grid size
  !
  ALLOCATE ( dmef(3, nbndsub, nbndsub, 2 * nkf) )
  IF (vme) ALLOCATE ( vmef(3, nbndsub, nbndsub, 2 * nkf) )
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
     CALL hamwan2bloch &
          ( nbndsub, nrr_k, irvec, ndegen_k, xxk, cufkk, etf (:, ik), chw)
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
  WRITE(6,'(/5x,a,f10.6,a)') 'Fermi energy coarse grid = ', ef * ryd2ev, ' eV'
  !
  IF( efermi_read ) THEN
     !
     ef = fermi_energy
     WRITE(6,'(/5x,a)') repeat('=',67)
     WRITE(6, '(/5x,a,f10.6,a)') &
         'Fermi energy is read from the input file: Ef = ', ef * ryd2ev, ' eV'
     WRITE(6,'(/5x,a)') repeat('=',67)
     !
  ELSEIF( band_plot ) THEN 
     !
     WRITE(6,'(/5x,a)') repeat('=',67)
     WRITE(stdout, '(/5x,"Fermi energy corresponds to the coarse k-mesh")')
     WRITE(6,'(/5x,a)') repeat('=',67) 
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
           WRITE(6,'(/5x,"Skipping the first ",i4," bands:")') nbndskip
           WRITE(6,'(/5x,"The Fermi level will be determined with ",f9.5," electrons")') nelec
        ENDIF
     ENDIF
     !
     ! Fermi energy
     !  
     ! since wkf(:,ikq) = 0 these bands do not bring any contribution to Fermi level
     !  
     IF (parallel_k) efnew = efermig(etf, nbndsub, nkqf, nelec, wkf, degaussw, ngaussw, 0, isk)
     IF (parallel_q) THEN 
#ifdef __PARA
       IF (mpime .eq. ionode_id) THEN
#endif
         efnew = efermig_seq(etf, nbndsub, nkqf, nelec, wkf, degaussw, ngaussw, 0, isk)
         ! etf on the full k-grid is later required for selfen_phon_k          
         etf_k = etf
#ifdef __PARA
       ENDIF
       CALL mp_bcast (efnew, ionode_id, inter_pool_comm)
       CALL mp_bcast (etf_k, ionode_id, inter_pool_comm)
#endif     
     ENDIF
     !
     WRITE(6, '(/5x,a,f10.6,a)') &
         'Fermi energy is calculated from the fine k-mesh: Ef = ', efnew * ryd2ev, ' eV'
     !
     ! if 'fine' Fermi level differs by more than 250 meV, there is probably something wrong
     ! with the wannier functions, or 'coarse' Fermi level is inaccurate
     IF (abs(efnew - ef) * ryd2eV .gt. 0.250d0 .and. (.not.eig_read) ) &
        WRITE(6,'(/5x,a)') 'Warning: check if difference with Fermi level fine grid makes sense'
     WRITE(6,'(/5x,a)') repeat('=',67)
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
  IF (parallel_k) THEN
    !
    ! get the size of the matrix elements stored in each pool
    ! for informational purposes.  Not necessary
    !
    CALL mem_size(ibndmin, ibndmax, nmodes, nkf, nqf)
    !
    IF (etf_mem) THEN
       ! Fine mesh set of g-matrices.  It is large for memory storage
       ALLOCATE ( epf17 (nkf, ibndmax-ibndmin+1, ibndmax-ibndmin+1, nmodes) )
       !
    ELSE
       !
       !  open epf and etf files with the correct record length
       !
       lrepmatf  = 2 * (ibndmax-ibndmin+1) * (ibndmax-ibndmin+1)
       CALL diropn (iunepmatf, 'epf', lrepmatf, exst)
       !
       lretf     = (ibndmax-ibndmin+1)
       CALL diropn (iuetf, 'etf', lretf, exst)
       !
    ENDIF
    !      
    DO iq = 1, nqf
       !   
       CALL start_clock ( 'ep-interp' )
       !
       ! In case of big calculation, show progression of iq (especially usefull when
       ! elecselfen = true as nothing happen during the calculation otherwise. 
       !
       IF ((.not. etf_mem) .AND. (.not. phonselfen) ) THEN 
         IF (MOD(iq,100) == 0) THEN
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
       CALL dynwan2bloch &
            ( nmodes, nrr_q, irvec, ndegen_q, xxq, uf, w2 )
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
       CALL ephwan2blochp &
            ( nmodes, xxq, irvec, ndegen_q, nrr_q, uf, epmatwef, nbndsub, nrr_k )
       !
       !
       !  number of k points with a band on the Fermi surface
       fermicount = 0
       !
       ! this is a loop over k blocks in the pool
       ! (size of the local k-set)
       DO ik = 1, nkf
          !
          ! xkf is assumed to be in crys coord
          !
          ikk = 2 * ik - 1
          ikq = ikk + 1
          !
          xkk = xkf(:, ikk)
          xkq = xkk + xxq
          !
          ! ------------------------------------------------------        
          ! hamiltonian : Wannier -> Bloch 
          ! ------------------------------------------------------
          !
          ! Kohn-Sham first, then get the rotation matricies for following interp.
          IF (eig_read) THEN
             CALL hamwan2bloch &
               ( nbndsub, nrr_k, irvec, ndegen_k, xkk, cufkk, etf_ks (:, ikk), chw_ks)
             CALL hamwan2bloch &
               ( nbndsub, nrr_k, irvec, ndegen_k, xkq, cufkq, etf_ks (:, ikq), chw_ks)
          ENDIF
          !
          CALL hamwan2bloch &
               ( nbndsub, nrr_k, irvec, ndegen_k, xkk, cufkk, etf (:, ikk), chw)
          CALL hamwan2bloch &
               ( nbndsub, nrr_k, irvec, ndegen_k, xkq, cufkq, etf (:, ikq), chw)
          !
          ! ------------------------------------------------------        
          !  dipole: Wannier -> Bloch
          ! ------------------------------------------------------        
          !
          CALL dmewan2bloch &
               ( nbndsub, nrr_k, irvec, ndegen_k, xkk, cufkk, dmef (:,:,:, ikk), etf(:,ikk), etf_ks(:,ikk))
          CALL dmewan2bloch &
               ( nbndsub, nrr_k, irvec, ndegen_k, xkq, cufkq, dmef (:,:,:, ikq), etf(:,ikq), etf_ks(:,ikq))
          !
          ! ------------------------------------------------------        
          !  velocity: Wannier -> Bloch
          ! ------------------------------------------------------        
          !
          IF (vme) THEN
             IF (eig_read) THEN
                CALL vmewan2bloch &
                     ( nbndsub, nrr_k, irvec, ndegen_k, xkk, cufkk, vmef (:,:,:, ikk), etf(:,ikk), etf_ks(:,ikk), chw_ks)
                CALL vmewan2bloch &
                     ( nbndsub, nrr_k, irvec, ndegen_k, xkq, cufkq, vmef (:,:,:, ikq), etf(:,ikq), etf_ks(:,ikq), chw_ks)
             ELSE
                CALL vmewan2bloch &
                     ( nbndsub, nrr_k, irvec, ndegen_k, xkk, cufkk, vmef (:,:,:, ikk), etf(:,ikk), etf_ks(:,ikk), chw)
                CALL vmewan2bloch &
                     ( nbndsub, nrr_k, irvec, ndegen_k, xkq, cufkq, vmef (:,:,:, ikq), etf(:,ikq), etf_ks(:,ikq), chw)
             ENDIF
          ENDIF
          !
          !
          IF (.not. etf_mem) THEN
             nrec  = ikk
             CALL davcio ( etf (ibndmin:ibndmax, ikk), ibndmax-ibndmin+1, iuetf, nrec, + 1)
             nrec  = ikq
             CALL davcio ( etf (ibndmin:ibndmax, ikq), ibndmax-ibndmin+1, iuetf, nrec, + 1)
             !
          ENDIF
          !
          ! interpolate ONLY when (k,k+q) both have at least one band 
          ! within a Fermi shell of size fsthick 
          !
          IF ( (( minval ( abs(etf (:, ikk) - ef) ) .lt. fsthick ) .and. &
               ( minval ( abs(etf (:, ikq) - ef) ) .lt. fsthick )) .or. indabs ) THEN
             !
             !  fermicount = fermicount + 1
             !
             ! --------------------------------------------------------------
             ! epmat : Wannier el and Bloch ph -> Bloch el and Bloch ph
             ! --------------------------------------------------------------
             !
             !
             CALL ephwan2bloch &
                  ( nbndsub, nrr_k, irvec, ndegen_k, epmatwef, xkk, cufkk, cufkq, epmatf, nmodes )
             !
             IF (lpolar) THEN
               !
               CALL compute_bmn_para2( nbndsub, nkstot, cufkk, cufkq, bmatf )
               !
               IF ( (abs(xxq(1)).gt.eps) .or. (abs(xxq(2)).gt.eps) .or. (abs(xxq(3)).gt.eps) ) THEN
                  CALL cryst_to_cart (1, xxq, bg, 1)
                  DO ibnd = 1, nbndsub
                    DO jbnd = 1, nbndsub
                      CALL rgd_blk_epw2(nq1, nq2, nq3, xxq, uf, epmatf(ibnd,jbnd,:), &
                            nmodes, epsi, zstar, bmatf(ibnd,jbnd), +1.d0)
                    ENDDO
                  ENDDO
                  CALL cryst_to_cart (1, xxq, at, -1)
               ENDIF
               !
             ENDIF
             !
             ! 
             ! write epmatf to file / store in memory
             !
             !
             DO imode = 1, nmodes
                !
                IF (etf_mem) THEN
                   !
                   DO jbnd = ibndmin, ibndmax
                      DO ibnd = ibndmin, ibndmax
                         !
                         epf17(ik,jbnd-ibndmin+1,ibnd-ibndmin+1,imode) = epmatf(jbnd,ibnd,imode)
                         !
                      ENDDO
                   ENDDO
                   !
                ELSE
                   !
                   nrec = (imode-1) * nkf + ik
                   CALL dasmio ( epmatf(ibndmin:ibndmax,ibndmin:ibndmax,imode), &
                        ibndmax-ibndmin+1, lrepmatf, iunepmatf, nrec, +1)
                   !
                ENDIF
                !
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
       ENDDO  ! end loop over k points
       !
       IF (phonselfen ) CALL selfen_phon_q( iq )
       IF (elecselfen ) CALL selfen_elec_q( iq )
       IF (nest_fn    ) CALL nesting_fn_q( iq )
       IF (specfun    ) CALL spectral_func_q( iq )
  !     IF (indabs    ) CALL indabs (iq)
  !     IF (twophoton ) CALL twophoton (iq)
       IF (ephwrite) THEN
          IF ( iq .eq. 1 ) THEN 
             CALL kmesh_fine
             CALL kqmap_fine
          ENDIF
          CALL write_ephmat( iq ) 
          CALL count_kpoints(iq)
       ENDIF
       !
       CALL stop_clock ( 'ep-interp' )
       !
    ENDDO  ! end loop over q points
    ! 
  ENDIF ! end parallel_k
  ! 
  IF (parallel_q) THEN
    !      
    ! get the size of the matrix elements stored in each pool
    ! for informational purposes.  Not necessary
    !
    CALL mem_size(ibndmin, ibndmax, nmodes, nkf, nqf)
    !
    IF (etf_mem) THEN
       ! Fine mesh set of g-matrices.  It is large for memory storage
       ALLOCATE ( epf17 (nqf, ibndmax-ibndmin+1, ibndmax-ibndmin+1, nmodes) )
       !
    ELSE
       !
       !  open epf and etf files with the correct record length
       !
       lrepmatf  = 2 * (ibndmax-ibndmin+1) * (ibndmax-ibndmin+1)
       CALL diropn (iunepmatf, 'epf', lrepmatf, exst)
       !
       lretf     = (ibndmax-ibndmin+1)
       CALL diropn (iuetf, 'etf', lretf, exst)
       !
    ENDIF
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
             ( nmodes, nrr_q, irvec, ndegen_q, xxq, uf, w2 )
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
             ( nmodes, xxq, irvec, ndegen_q, nrr_q, uf, epmatwef, nbndsub, nrr_k )
        !
        !  number of k points with a band on the Fermi surface
        fermicount = 0
        !
        ! this is a loop over k blocks in the pool
        ! (size of the local k-set)
        !
        ! ------------------------------------------------------        
        ! hamiltonian : Wannier -> Bloch
        ! ------------------------------------------------------
        !
        ! Kohn-Sham first, then get the rotation matricies for following interp.
        IF (eig_read) THEN
           CALL hamwan2bloch &
             ( nbndsub, nrr_k, irvec, ndegen_k, xkk, cufkk, etf_ks (:, ikk), chw_ks)     
           CALL hamwan2bloch &
             ( nbndsub, nrr_k, irvec, ndegen_k, xkq, cufkq, etf_ks (:, ikq), chw_ks)
        ENDIF
        !
        CALL hamwan2bloch &
             ( nbndsub, nrr_k, irvec, ndegen_k, xkk, cufkk, etf (:, ikk), chw)        
        CALL hamwan2bloch &
             ( nbndsub, nrr_k, irvec, ndegen_k, xkq, cufkq, etf (:, ikq), chw)
        !
        ! ------------------------------------------------------        
        !  dipole: Wannier -> Bloch
        ! ------------------------------------------------------        
        !
        CALL dmewan2bloch &
             ( nbndsub, nrr_k, irvec, ndegen_k, xkk, cufkk, dmef (:,:,:, ikk), etf(:,ikk), etf_ks(:,ikk))
        CALL dmewan2bloch &
             ( nbndsub, nrr_k, irvec, ndegen_k, xkq, cufkq, dmef (:,:,:, ikq), etf(:,ikq), etf_ks(:,ikq))
        !
        ! ------------------------------------------------------        
        !  velocity: Wannier -> Bloch
        ! ------------------------------------------------------        
        !
        IF (vme) THEN
           IF (eig_read) THEN
              CALL vmewan2bloch &
                   ( nbndsub, nrr_k, irvec, ndegen_k, xkk, cufkk, vmef (:,:,:, ikk), etf(:,ikk), etf_ks(:,ikk), chw_ks)
              CALL vmewan2bloch &
                   ( nbndsub, nrr_k, irvec, ndegen_k, xkq, cufkq, vmef (:,:,:, ikq), etf(:,ikq), etf_ks(:,ikq), chw_ks)
           ELSE
              CALL vmewan2bloch &
                   ( nbndsub, nrr_k, irvec, ndegen_k, xkk, cufkk, vmef (:,:,:, ikk), etf(:,ikk), etf_ks(:,ikk), chw)
              CALL vmewan2bloch &
                   ( nbndsub, nrr_k, irvec, ndegen_k, xkq, cufkq, vmef (:,:,:, ikq), etf(:,ikq), etf_ks(:,ikq), chw)
           ENDIF
        ENDIF
        !
        !
        IF (.not. etf_mem) THEN
           nrec  = ikk
           CALL davcio ( etf (ibndmin:ibndmax, ikk), ibndmax-ibndmin+1, iuetf, nrec, + 1)
           nrec  = ikq
           CALL davcio ( etf (ibndmin:ibndmax, ikq), ibndmax-ibndmin+1, iuetf, nrec, + 1)
           !
        ENDIF
        !
        ! interpolate ONLY when (k,k+q) both have at least one band 
        ! within a Fermi shell of size fsthick 
        !
        IF ( (( minval ( abs(etf (:, ikk) - ef) ) .lt. fsthick ) .and. &
             ( minval ( abs(etf (:, ikq) - ef) ) .lt. fsthick )) .or. indabs ) THEN
           !
           !  fermicount = fermicount + 1
           !
           ! --------------------------------------------------------------
           ! epmat : Wannier el and Bloch ph -> Bloch el and Bloch ph
           ! --------------------------------------------------------------
           !
           !
           CALL ephwan2bloch &
                ( nbndsub, nrr_k, irvec, ndegen_k, epmatwef, xkk, cufkk, cufkq, epmatf, nmodes )
           !
           ! 
           IF (lpolar) THEN
             !
             CALL compute_bmn_para2( nbndsub, nkstot, cufkk, cufkq, bmatf )
             !
             IF ( (abs(xxq(1)).gt.eps) .or. (abs(xxq(2)).gt.eps) .or. (abs(xxq(3)).gt.eps) ) THEN
                CALL cryst_to_cart (1, xxq, bg, 1)
                DO ibnd = 1, nbndsub
                  DO jbnd = 1, nbndsub
                    CALL rgd_blk_epw2(nq1, nq2, nq3, xxq, uf, epmatf(ibnd,jbnd,:), &
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
           DO imode = 1, nmodes
              !
              IF (etf_mem) THEN
                 !
                 DO jbnd = ibndmin, ibndmax
                    DO ibnd = ibndmin, ibndmax
                       !
                       epf17(iq,jbnd-ibndmin+1,ibnd-ibndmin+1,imode) = epmatf(jbnd,ibnd,imode)
                       !
                    ENDDO
                 ENDDO
                 !
              ELSE
                 !
                 nrec = (imode-1) * nqf + iq
                 CALL dasmio ( epmatf(ibndmin:ibndmax,ibndmin:ibndmax,imode), &
                      ibndmax-ibndmin+1, lrepmatf, iunepmatf, nrec, +1)
                 !
              ENDIF
              !
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
      IF (specfun    ) CALL spectral_func_k( ik )
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
  ! ---------------------------------------------------------------------------------------
  ! ---------------------------------------------------------------------------------------  
  !
  ! SP: Added lambda and phonon lifetime writing to file.
  ! 
#ifdef __PARA
  CALL mp_barrier(inter_pool_comm)
  IF (mpime.eq.ionode_id) THEN
#endif
    !
    IF (phonselfen .and. parallel_k ) THEN
      OPEN(unit=lambda_phself,file='lambda.phself')
      WRITE(lambda_phself, '(/2x,a/)') '#Lambda phonon self-energy'
      WRITE(lambda_phself, *) '#Modes     ',(imode, imode=1,nmodes)
      DO iq = 1, nqtotf
          !
        myfmt = "(*(3x,E15.5))"
        WRITE(lambda_phself,'(i9,4x)',advance='no') iq
        WRITE(lambda_phself, fmt=myfmt) (REAL(lambda_all(imode,iq,1)),imode=1,nmodes)
          !
      ENDDO
      CLOSE(lambda_phself)
      OPEN(unit=linewidth_phself,file='linewidth.phself')
      WRITE(linewidth_phself, '(/2x,a/)') '#Phonon lifetime (meV) '
      WRITE(linewidth_phself,'(2x,a)',advance='no') '#Q-point     '
      Do imode=1, nmodes
        WRITE(linewidth_phself, '(a)',advance='no') '      Mode'
        WRITE(linewidth_phself, '(i3)',advance='no') imode
      enddo
      WRITE(linewidth_phself, '(/2x,a/)') ' '
      DO iq = 1, nqtotf
        !
        myfmt = "(*(3x,E15.5))"
        WRITE(linewidth_phself,'(i9,4x)',advance='no') iq
        WRITE(linewidth_phself, fmt=myfmt) (ryd2mev*REAL(gamma_all(imode,iq,1)), imode=1,nmodes)
        !
      ENDDO
      CLOSE(linewidth_phself)
    ENDIF
#ifdef __PARA
  ENDIF
#endif
  IF (band_plot) CALL plot_band
  !
  IF (a2f) CALL eliashberg_a2f
  ! 
  IF ( ALLOCATED(lambda_all) )   DEALLOCATE( lambda_all )
  IF ( ALLOCATED(gamma_all) )   DEALLOCATE( gamma_all )
  IF ( ALLOCATED(sigmai_all) )   DEALLOCATE( sigmai_all )
  IF ( ALLOCATED(sigmai_mode) )   DEALLOCATE( sigmai_mode )
  !
  !
  CALL stop_clock ( 'ephwann' )
  !
  END SUBROUTINE ephwann_shuffle
  !
!-------------------------------------------
SUBROUTINE epw_write
!-------------------------------------------
  !
  USE kinds,     ONLY : DP
  USE epwcom,    ONLY : nbndsub, vme, eig_read, etf_mem
  USE pwcom,     ONLY : ef
  USE elph2,     ONLY : nrr_k, nrr_q, chw, rdw, cdmew, cvmew, chw_ks, &
                        zstar, epsi, epmatwp
  USE phcom,     ONLY : nmodes  
  USE io_epw,    ONLY : epwdata, iundmedata, iunvmedata, iunksdata, iunepmatwp
  USE io_files,  ONLY : prefix, diropn
#ifdef __PARA
  USE mp,        ONLY : mp_barrier
  USE mp_global, ONLY : my_pool_id,inter_pool_comm
  USE mp_world,  ONLY : mpime
  USE io_global, ONLY : ionode_id
#endif
  !
  implicit none
  LOGICAL             :: exst
  INTEGER             :: ibnd, jbnd, jmode, imode, irk, irq, ipol, i, lrepmatw
  character (len=256) :: filint
  complex(kind=DP)    :: aux ( nbndsub*nbndsub*nrr_k*nmodes*nrr_q ) 
  !
  WRITE(6,'(/5x,"Writing Hamiltonian, Dynamical matrix and EP vertex in Wann rep to file"/)')
  !
#ifdef __PARA
     IF (mpime.eq.ionode_id) THEN
#endif     
       !
       OPEN(unit=epwdata,file='epwdata.fmt')
       OPEN(unit=iundmedata,file='dmedata.fmt')
       IF (vme) OPEN(unit=iunvmedata,file='vmedata.fmt')
       IF (eig_read) OPEN(unit=iunksdata,file='ksdata.fmt')
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
       IF (etf_mem) THEN
         lrepmatw   = 2 * nbndsub * nbndsub * nrr_k * nmodes * nrr_q
         i = 0
         DO irq = 1, nrr_q
           DO imode = 1, nmodes
             DO irk = 1, nrr_k
               DO jbnd = 1, nbndsub
                 DO ibnd = 1, nbndsub
                   i = i + 1
                   aux (i) = epmatwp(ibnd,jbnd,irk,imode,irq)
                 ENDDO
               ENDDO
             ENDDO
           ENDDO
         ENDDO
         filint    = trim(prefix)//'.epmatwp'
         CALL diropn (iunepmatwp, 'epmatwp', lrepmatw, exst)
         CALL davcio ( aux, lrepmatw, iunepmatwp, 1, +1 )
         CLOSE(iunepmatwp)
         IF (ALLOCATED(epmatwp)) DEALLOCATE ( epmatwp )
       ENDIF 
       !
       CLOSE(epwdata)
       CLOSE(iundmedata)
       IF (vme) CLOSE(iunvmedata)
       IF (eig_read) CLOSE(iunksdata)
       !
#ifdef __PARA
    ENDIF
    CALL mp_barrier(inter_pool_comm)
#endif     
     !
!---------------------------------
END SUBROUTINE epw_write
!---------------------------------
!---------------------------------
SUBROUTINE epw_read()
!---------------------------------
  USE kinds,     ONLY : DP
  USE epwcom,    ONLY : nbndsub, vme, eig_read, wepexst, etf_mem
  USE pwcom,     ONLY : ef
  USE elph2,     ONLY : nrr_k, nrr_q, chw, rdw, epmatwp, &
                        cdmew, cvmew, chw_ks, zstar, epsi
  USE phcom,     ONLY : nmodes  
  USE io_global, ONLY : stdout
  USE io_files,  ONLY : prefix, diropn
  USE io_epw,    ONLY : epwdata, iundmedata, iunvmedata, iunksdata, iunepmatwp
  USE constants_epw, ONLY :  czero
#ifdef __NAG
  USE f90_unix_io,ONLY : flush
#endif
#ifdef __PARA
  USE io_global, ONLY : ionode_id
  USE mp,        ONLY : mp_barrier, mp_bcast
  USE mp_global, ONLY : my_pool_id, &
                        intra_pool_comm, inter_pool_comm, root_pool
  USE mp_world,  ONLY : mpime
#endif
  !
  implicit none
  !
  !
  LOGICAL             :: exst
  character (len=256) :: filint
  INTEGER             :: ibnd, jbnd, jmode, imode, irk, irq, &
                         ipol, ios, i, lrepmatw
  complex(kind=DP)    :: aux ( nbndsub*nbndsub*nrr_k*nmodes*nrr_q )
     !
  WRITE(stdout,'(/5x,"Reading Hamiltonian, Dynamical matrix and EP vertex in Wann rep from file"/)')
  call flush(6)
#ifdef __PARA
  IF (mpime.eq.ionode_id) THEN
#endif
    !
    OPEN(unit=epwdata,file='epwdata.fmt',status='old',iostat=ios)
    OPEN(unit=iundmedata,file='dmedata.fmt',status='old',iostat=ios)
    IF (eig_read) OPEN(unit=iunksdata,file='ksdata.fmt',status='old',iostat=ios)
    IF (vme) OPEN(unit=iunvmedata,file='vmedata.fmt',status='old',iostat=ios)
    IF (ios /= 0) call errore ('ephwann_shuffle', 'error opening epwdata.fmt',epwdata)
    READ (epwdata,*) ef
    READ (epwdata,*) nbndsub, nrr_k, nmodes, nrr_q
    READ (epwdata,*) zstar, epsi
    ! 
#ifdef __PARA
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
#endif
  !
  IF (.not. ALLOCATED(epmatwp)) ALLOCATE ( epmatwp ( nbndsub, nbndsub, nrr_k, nmodes, nrr_q) )
  IF (.not. ALLOCATED(chw)    ) ALLOCATE ( chw ( nbndsub, nbndsub, nrr_k )            )
  IF (.not. ALLOCATED(chw_ks) ) ALLOCATE ( chw_ks ( nbndsub, nbndsub, nrr_k )         )
  IF (.not. ALLOCATED(rdw)    ) ALLOCATE ( rdw ( nmodes,  nmodes,  nrr_q )            )
  IF (.not. ALLOCATED(cdmew)  ) ALLOCATE ( cdmew ( 3, nbndsub, nbndsub, nrr_k )       )
  IF (vme .and. (.not.ALLOCATED(cvmew))  ) ALLOCATE ( cvmew   ( 3, nbndsub, nbndsub, nrr_k )     )
  !
#ifdef __PARA
  IF (mpime.eq.ionode_id) THEN
#endif
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
     DO imode = 1, nmodes
        DO jmode = 1, nmodes
           DO irq = 1, nrr_q
              READ (epwdata,*) rdw(imode,jmode,irq)
           ENDDO
        ENDDO
     ENDDO
     !
#ifdef __PARA
  ENDIF
  !
  CALL mp_bcast (chw, ionode_id, inter_pool_comm)
  CALL mp_bcast (chw, root_pool, intra_pool_comm)
  !
  IF (eig_read) CALL mp_bcast (chw_ks, ionode_id, inter_pool_comm)
  IF (eig_read) CALL mp_bcast (chw_ks, root_pool, intra_pool_comm)
  !
  CALL mp_bcast (rdw, ionode_id, inter_pool_comm)
  CALL mp_bcast (rdw, root_pool, intra_pool_comm)
  !
  CALL mp_bcast (cdmew, ionode_id, inter_pool_comm)
  CALL mp_bcast (cdmew, root_pool, intra_pool_comm)
  !
  IF (vme) CALL mp_bcast (cvmew, ionode_id, inter_pool_comm)
  IF (vme) CALL mp_bcast (cvmew, root_pool, intra_pool_comm)
  !
#endif
  !
  !
  IF (etf_mem) then
    epmatwp = czero
#ifdef __PARA
    IF (mpime.eq.ionode_id) THEN
#endif
      !
      lrepmatw   = 2 * nbndsub * nbndsub * nrr_k * nmodes * nrr_q
      filint    = trim(prefix)//'.epmatwp'
      CALL diropn (iunepmatwp, 'epmatwp', lrepmatw, exst)
      CALL davcio ( aux, lrepmatw, iunepmatwp, 1, -1 )
      i = 0
      DO irq = 1, nrr_q
        DO imode = 1, nmodes
          DO irk = 1, nrr_k
            DO jbnd = 1, nbndsub
              DO ibnd = 1, nbndsub
                i = i + 1
                epmatwp(ibnd,jbnd,irk,imode,irq) = aux(i)
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO
#ifdef __PARA
    ENDIF
    !
    CALL mp_bcast (epmatwp, ionode_id, inter_pool_comm)
    CALL mp_bcast (epmatwp, root_pool, intra_pool_comm)
    !
#endif
    !
  ENDIF
  !
#ifdef __PARA
  CALL mp_barrier(inter_pool_comm)
  IF (mpime.eq.ionode_id) THEN
#endif
    CLOSE(epwdata)
    CLOSE(iundmedata)
    IF (vme) CLOSE(iunvmedata)
#ifdef __PARA
  ENDIF
#endif
  !
  WRITE(stdout,'(/5x,"Finished reading Wann rep data from file"/)')
  !
!---------------------------------
END SUBROUTINE epw_read
!---------------------------------
!---------------------------------
SUBROUTINE mem_size(ibndmin, ibndmax, nmodes, nkf, nqf) 
!---------------------------------
!
!  SUBROUTINE estimates the amount of memory taken up by 
!  the <k+q| dV_q,nu |k> on the fine meshes and prints 
!  out a useful(?) message   
!
  USE io_global,     ONLY : stdout
  USE kinds,         ONLY : DP
  !
  implicit none
  !
  integer :: imelt, ibndmin, ibndmax, nmodes, nkf, nqf
  real(kind=DP)    :: rmelt
  character (len=256) :: chunit
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

!--------------------------------------------------------------------
FUNCTION efermig_seq (et, nbnd, nks, nelec, wk, Degauss, Ngauss, is, isk)
  !--------------------------------------------------------------------
  !
  !     Finds the Fermi energy - Gaussian Broadening
  !     (see Methfessel and Paxton, PRB 40, 3616 (1989 )
  !
  USE io_global, ONLY : stdout
  USE kinds,     ONLY : DP
  USE constants, ONLY : rytoev
  implicit none
  !  I/O variables
  integer, intent(in) :: nks, nbnd, Ngauss, is, isk(nks)
  real(DP), intent(in) :: wk (nks), et (nbnd, nks), Degauss, nelec
  real(DP) :: efermig_seq
  !
  real(DP), parameter :: eps= 1.0d-10
  integer, parameter :: maxiter = 300
  ! internal variables
  real(DP) :: Ef, Eup, Elw, sumkup, sumklw, sumkmid
  real(DP), external::  sumkg_seq
  integer :: i, kpoint
  !
  !      find bounds for the Fermi energy. Very safe choice!
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
  !      Bisection method
  !
  sumkup = sumkg_seq (et, nbnd, nks, wk, Degauss, Ngauss, Eup, is, isk)
  sumklw = sumkg_seq (et, nbnd, nks, wk, Degauss, Ngauss, Elw, is, isk)
  if ( (sumkup - nelec) < -eps .or. (sumklw - nelec) > eps )  &
       call errore ('efermig_seq', 'internal error, cannot bracket Ef', 1)
  do i = 1, maxiter
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
  enddo
  if (is /= 0) WRITE(stdout, '(5x,"Spin Component #",i3)') is
  WRITE( stdout, '(5x,"Warning: too many iterations in bisection"/ &
       &      5x,"Ef = ",f10.6," sumk = ",f10.6," electrons")' ) &
       Ef * rytoev, sumkmid
  !
  efermig_seq = Ef
  return
end FUNCTION efermig_seq

!-----------------------------------------------------------------------
function sumkg_seq (et, nbnd, nks, wk, degauss, ngauss, e, is, isk)
  !-----------------------------------------------------------------------
  !
  !     This function computes the number of states under a given energy e
  !
  !
  USE kinds
  USE mp_pools, ONLY : inter_pool_comm
  USE mp,       ONLY : mp_sum
  implicit none
  ! Output variable
  real(DP) :: sumkg_seq
  ! Input variables
  integer, intent(in) :: nks, nbnd, ngauss
  ! input: the total number of K points
  ! input: the number of bands
  ! input: the type of smearing
  real(DP), intent(in) :: wk (nks), et (nbnd, nks), degauss, e
  ! input: the weight of the k points
  ! input: the energy eigenvalues 
  ! input: gaussian broadening
  ! input: the energy to check
  integer, intent(in) :: is, isk(nks)
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
  do ik = 1, nks
     sum1 = 0.d0
     if (is /= 0) then
        if (isk(ik).ne.is) cycle
     end if
     do ibnd = 1, nbnd
        sum1 = sum1 + wgauss ( (e-et (ibnd, ik) ) / degauss, ngauss)
     enddo
     sumkg_seq = sumkg_seq + wk (ik) * sum1
  enddo
  return
end function sumkg_seq
!
  !-----------------------------------------------------------------
  subroutine rwepmatw ( epmatw, nbnd, np, nmodes, nrec, iun, iop)
  !-----------------------------------------------------------------
  !
  ! A simple wrapper to the davcio routine to read/write arrays
  ! instead of vectors 
  !-----------------------------------------------------------------
  USE kinds, only : DP
#ifdef __PARA
  use mp, only : mp_barrier
  use mp_global, only : my_pool_id
#endif
  implicit none
  integer :: lrec, iun, nrec, iop, i, nbnd, np, nmodes, ibnd, jbnd, imode, ip
  !
  ! np is either nrr_k or nq (epmatwe and epmatwp have the same structure)
  !
  complex(kind=DP):: epmatw(nbnd,nbnd,np,nmodes), &
     aux ( nbnd*nbnd*np*nmodes )
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
  end subroutine rwepmatw



