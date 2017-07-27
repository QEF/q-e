  !                                                                            
  ! Copyright (C) 2010-2016 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino
  ! Copyright (C) 2007-2009 Jesse Noffsinger, Brad Malone, Feliciano Giustino  
  !                                                                            
  ! This file is distributed under the terms of the GNU General Public         
  ! License. See the file `LICENSE' in the root directory of the               
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .             
  !                                                                            
  !----------------------------------------------------------------------
  SUBROUTINE ephwann_shuffle_mem( nqc, xqc )
  !---------------------------------------------------------------------
  !!
  !!  Wannier interpolation of electron-phonon vertex (memory optimized version)
  !!
  !!  Scalar implementation   Feb 2006
  !!  Parallel version        May 2006
  !!  Disentenglement         Oct 2006
  !!  Compact formalism       Dec 2006
  !!  Phonon irreducible zone Mar 2007
  !!
  !-----------------------------------------------------------------------
  !
  USE kinds,         ONLY : DP
  USE pwcom,         ONLY : nbnd, nks, nkstot, isk, &
                            et, xk, ef,  nelec
  USE cell_base,     ONLY : at, bg, omega, alat
  USE start_k,       ONLY : nk1, nk2, nk3
  USE ions_base,     ONLY : nat, amass, ityp, tau
  USE phcom,         ONLY : nq1, nq2, nq3, nmodes
  USE epwcom,        ONLY : nbndsub, lrepmatf, fsthick, epwread, longrange,     &
                            epwwrite, ngaussw, degaussw, lpolar, lifc,          &
                            nbndskip, parallel_k, parallel_q, etf_mem,          &
                            elecselfen, phonselfen, nest_fn, a2f,               &
                            vme, eig_read, ephwrite, nkf1, nkf2, nkf3,          & 
                            efermi_read, fermi_energy, specfun, band_plot,      &
                            nqf1, nqf2, nqf3, mp_mesh_k, restart, prtgkk
  USE noncollin_module, ONLY : noncolin
  USE constants_epw, ONLY : ryd2ev, ryd2mev, one, two, czero, twopi, ci, zero
  USE io_files,      ONLY : prefix, diropn
  USE io_global,     ONLY : stdout, ionode
  USE io_epw,        ONLY : lambda_phself, linewidth_phself, iunepmatwe,        &
                            iunepmatwp, crystal
  USE elph2,         ONLY : nrr_k, nrr_q, cu, cuq, lwin, lwinq, irvec, ndegen_k,&
                            ndegen_q,  wslen, chw, chw_ks, cvmew, cdmew, rdw,   &
                            epmatwp, epmatq, wf, etf, etf_k, etf_ks, xqf, xkf,  &
                            wkf, dynq, nqtotf, nkqf, epf17, nkf, nqf, et_ks,    &
                            ibndmin, ibndmax, lambda_all, dmec, dmef, vmef,     &
                            sigmai_all, sigmai_mode, gamma_all, epsi, zstar,    &
                            efnew, ifc, sigmar_all, zi_all, nkqtotf
#if defined(__NAG)
  USE f90_unix_io,   ONLY : flush
#endif
  USE mp,            ONLY : mp_barrier, mp_bcast, mp_sum
  USE io_global,     ONLY : ionode_id
  USE mp_global,     ONLY : inter_pool_comm, intra_pool_comm, root_pool
  USE mp_world,      ONLY : mpime
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
  LOGICAL :: opnd
  !! Check whether the file is open.
  !
  CHARACTER (len=256) :: filint
  !! Name of the file to write/read 
  CHARACTER (len=256) :: nameF
  !! Name of the file
  CHARACTER (len=30)  :: myfmt
  !! Variable used for formatting output
  ! 
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
  INTEGER :: nrec
  !! record index when reading file
  INTEGER :: lrepmatw
  !! record length while reading file
  INTEGER :: i 
  !! Index when writing to file
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
  COMPLEX(kind=DP), ALLOCATABLE :: epmatwe_mem  (:,:,:,:)
  !! e-p matrix  in wannier basis - electrons (written on disk)
  COMPLEX(kind=DP), ALLOCATABLE :: epmatwef (:,:,:)
  !! e-p matrix  in el wannier - fine Bloch phonon grid
  COMPLEX(kind=DP), ALLOCATABLE :: epmatf( :, :)
  !! e-p matrix  in smooth Bloch basis, fine mesh
  COMPLEX(kind=DP), ALLOCATABLE :: cufkk ( :, :, :)
  !! Rotation matrix, fine mesh, points k
  COMPLEX(kind=DP), ALLOCATABLE :: cufkk_dum ( :, :)
  !! Dummy rotation matrix  
  COMPLEX(kind=DP), ALLOCATABLE :: cufkq ( :, :, :)
  !! the same, for points k+q
  COMPLEX(kind=DP), ALLOCATABLE :: uf( :, :)
  !! Rotation matrix for phonons
  COMPLEX(kind=DP), ALLOCATABLE :: bmatf ( :, :)
  !! overlap U_k+q U_k^\dagger in smooth Bloch basis, fine mesh
  !COMPLEX(kind=DP), ALLOCATABLE :: cfac1(:)
  !COMPLEX(kind=DP), ALLOCATABLE :: cfacq1(:)
  COMPLEX(kind=DP), ALLOCATABLE :: cfac(:)
  !! Used to store $e^{2\pi r \cdot k}$ exponential 
  COMPLEX(kind=DP), ALLOCATABLE :: cfacq(:)
  !! Used to store $e^{2\pi r \cdot k+q}$ exponential 
  COMPLEX(kind=DP), ALLOCATABLE :: eptmp(:,:,:,:)
  !! Temporary el-ph matrices. 
  COMPLEX(kind=DP), ALLOCATABLE :: epmatlrT(:,:,:,:)
  !! Long-range temp. save   
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
  IF ( epwread ) THEN
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
    CALL mp_bcast (nat, ionode_id, inter_pool_comm)
    CALL mp_bcast (nat, root_pool, intra_pool_comm)  
    IF (mpime /= ionode_id) ALLOCATE( ityp( nat ) )
    CALL mp_bcast (nmodes, ionode_id, inter_pool_comm)
    CALL mp_bcast (nmodes, root_pool, intra_pool_comm)  
    CALL mp_bcast (nelec, ionode_id, inter_pool_comm)
    CALL mp_bcast (nelec, root_pool, intra_pool_comm)  
    CALL mp_bcast (at, ionode_id, inter_pool_comm)
    CALL mp_bcast (at, root_pool, intra_pool_comm)  
    CALL mp_bcast (bg, ionode_id, inter_pool_comm)
    CALL mp_bcast (bg, root_pool, intra_pool_comm)  
    CALL mp_bcast (omega, ionode_id, inter_pool_comm)
    CALL mp_bcast (omega, root_pool, intra_pool_comm)  
    CALL mp_bcast (alat, ionode_id, inter_pool_comm)
    CALL mp_bcast (alat, root_pool, intra_pool_comm)  
    IF (mpime /= ionode_id) ALLOCATE( tau( 3, nat ) )
    CALL mp_bcast (tau, ionode_id, inter_pool_comm)
    CALL mp_bcast (tau, root_pool, intra_pool_comm)  
    CALL mp_bcast (amass, ionode_id, inter_pool_comm)
    CALL mp_bcast (amass, root_pool, intra_pool_comm)  
    CALL mp_bcast (ityp, ionode_id, inter_pool_comm)
    CALL mp_bcast (ityp, root_pool, intra_pool_comm)  
    CALL mp_bcast (isk, ionode_id, inter_pool_comm)
    CALL mp_bcast (isk, root_pool, intra_pool_comm)    
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
#ifndef __MPI  
  ! Open like this only in sequential. Otherwize open with MPI-open
  IF (ionode) THEN
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
     IF (ionode) THEN
       lrepmatw   = 2 * nbndsub * nbndsub * nrr_k * nmodes
       filint    = trim(prefix)//'.epmatwe'
       CALL diropn (iunepmatwe, 'epmatwe', lrepmatw, exst)
       filint    = trim(prefix)//'.epmatwp'
       CALL diropn (iunepmatwp, 'epmatwp', lrepmatw, exst)    
     ENDIF          
     !
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
     ALLOCATE(epmatwe_mem ( nbndsub, nbndsub, nrr_k, nmodes))
     !
     ! Hamiltonian
     !
     CALL hambloch2wan &
          ( nbnd, nbndsub, nks, nkstot, et, xk, cu, lwin, nrr_k, irvec, wslen, chw )
     !
     ! Kohn-Sham eigenvalues
     !
     IF (eig_read) THEN
       WRITE (6,'(5x,a)') "Interpolating MB and KS eigenvalues"
       CALL hambloch2wan &
            ( nbnd, nbndsub, nks, nkstot, et_ks, xk, cu, lwin, nrr_k, irvec, wslen, chw_ks )
     ENDIF
     !
     ! Dipole
     !
    ! CALL dmebloch2wan &
    !      ( nbnd, nbndsub, nks, nkstot, nkstot, dmec, xk, cu, nrr_k, irvec, wslen )
     CALL dmebloch2wan &
          ( nbnd, nbndsub, nks, nkstot, dmec, xk, cu, nrr_k, irvec, wslen, lwin )
     !
     ! Dynamical Matrix 
     !
     IF (.not. lifc) CALL dynbloch2wan &
                          ( nmodes, nqc, xqc, dynq, nrr_q, irvec, wslen )
     !
     ! Transform of position matrix elements
     ! PRB 74 195118  (2006)
     IF (vme) CALL vmebloch2wan &
         ( nbnd, nbndsub, nks, nkstot, xk, cu, nrr_k, irvec, wslen )
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
         CALL ephbloch2wane &
           ( nbnd, nbndsub, nks, nkstot, xk, cu, cuq, &
           epmatq (:,:,:,imode,iq), nrr_k, irvec, wslen, epmatwe_mem(:,:,:,imode) )
         !
       ENDDO
       ! Only the master node writes 
       IF (ionode) THEN
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
        CALL ephbloch2wanp_mem &
         ( nbndsub, nmodes, xqc, nqc, irvec, nrr_k, nrr_q, epmatwe_mem )
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
  IF ( ALLOCATED (epmatwe_mem) ) DEALLOCATE (epmatwe_mem)
  IF ( ALLOCATED (epmatq) )  DEALLOCATE (epmatq)
  IF ( ALLOCATED (cu) )      DEALLOCATE (cu)
  IF ( ALLOCATED (cuq) )     DEALLOCATE (cuq)
  IF ( ALLOCATED (lwin) )    DEALLOCATE (lwin)
  IF ( ALLOCATED (lwinq) )   DEALLOCATE (lwinq)
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
  !
  ! at this point, we will interpolate the Wannier rep to the Bloch rep 
  ! for electrons, phonons and the ep-matrix
  !
  !  need to add some sort of parallelization (on g-vectors?)  what
  !  else can be done when we don't ever see the wfcs??
  !
  ! SP: Only k-para
  CALL loadqmesh_serial
  CALL loadkmesh_para
  !
  ALLOCATE ( epmatwef( nbndsub, nbndsub, nrr_k),                          &
       wf ( nmodes,  nqf ), etf ( nbndsub, nkqf),                         &
       etf_ks ( nbndsub, nkqf), cufkk_dum(nbndsub,nbndsub),               &
       epmatf( nbndsub, nbndsub), cufkk ( nbndsub, nbndsub, nkf), &
       cufkq ( nbndsub, nbndsub, nkf), uf ( nmodes, nmodes),              &
       bmatf( nbndsub, nbndsub)  )
  !
  ! Need to be initialized
  epmatf(:,:) = czero
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
  irvec_r = REAL(irvec,KIND=dp)
  ! 
  ! SP: Create a look-up table for the exponential of the factor. 
  !     This can only work with homogeneous fine grids.
  IF ( (nkf1 >0) .AND. (nkf2 > 0) .AND. (nkf3 > 0) .AND. &
       (nqf1 >0) .AND. (nqf2 > 0) .AND. (nqf3 > 0) .AND. .NOT. mp_mesh_k ) THEN
    ! Make a check   
    IF ((nqf1>nkf1) .or. (nqf2>nkf2) .or. (nqf3>nkf3)) &
            CALL errore('The fine q-grid cannot be larger than the fine k-grid',1)   
    ! Along x
    DO ikx = -2*nk1, 2*nk1
      DO ikfx = 0, nkf1-1
        !rdotk = twopi * ( xk(1)*irvec(1,ir))
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
     cfac(:) = exp( ci*rdotk ) / ndegen_k(:)
     ! 
     CALL hamwan2bloch &
          ( nbndsub, nrr_k, cufkk, etf (:, ik), chw, cfac)
     !
     !
  ENDDO
  !
  WRITE(6,'(/5x,a,f10.6,a)') 'Fermi energy coarse grid = ', ef * ryd2ev, ' eV'
  !
  IF( efermi_read ) THEN
     !
     ef = fermi_energy
     WRITE(stdout,'(/5x,a)') repeat('=',67)
     WRITE(stdout, '(/5x,a,f10.6,a)') &
         'Fermi energy is read from the input file: Ef = ', ef * ryd2ev, ' eV'
     WRITE(stdout,'(/5x,a)') repeat('=',67)
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
           WRITE(6,'(/5x,"Skipping the first ",i4," bands:")') nbndskip
           WRITE(6,'(/5x,"The Fermi level will be determined with ",f9.5," electrons")') nelec
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
  IF (parallel_k) THEN
    !
    ! get the size of the matrix elements stored in each pool
    ! for informational purposes.  Not necessary
    !
    CALL mem_size(ibndmin, ibndmax, nmodes, nkf)
    !
    ! Fine mesh set of g-matrices.  It is large for memory storage
    ! SP: Should not be a memory problem. If so, can always the number of cores to reduce nkf. 
    ALLOCATE ( epf17 (ibndmax-ibndmin+1, ibndmax-ibndmin+1, nmodes, nkf) )
    ALLOCATE ( eptmp (ibndmax-ibndmin+1, ibndmax-ibndmin+1, nmodes, nkf) )
    ALLOCATE ( epmatlrT (nbndsub, nbndsub, nmodes, nkf) )
    ! 
    ! Restart calculation
    iq_restart = 1
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
        CALL electron_read(iq_restart,nqf,nkqtotf/2,sigmar_all,sigmai_all,zi_all)
      ENDIF
      ! 
    ENDIF
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
    DO iq = iq_restart, nqf
       !   
       CALL start_clock ( 'ep-interp' )
       !
       ! In case of big calculation, show progression of iq (especially usefull when
       ! elecselfen = true as nothing happen during the calculation otherwise. 
       !
       IF (.not. phonselfen) THEN 
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
               ( nmodes, nrr_q, irvec, ndegen_q, xxq, uf, w2 )
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
       epf17(:,:,:,:) = czero
       eptmp(:,:,:,:) = czero
       epmatlrT(:,:,:,:) = czero
       cufkk(:,:,:) = czero
       cufkq(:,:,:) = czero
       ! 
       DO imode = 1, nmodes       
         ! 
         epmatwef(:,:,:) = czero
         !
         IF (.NOT. longrange) THEN
           CALL ephwan2blochp_mem (imode, nmodes, xxq, irvec, ndegen_q, nrr_q, uf, epmatwef, nbndsub, nrr_k )
         ENDIF
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
            ! Only do the following for the first mode because it is mode independent
            IF (imode==1) THEN
              ! 
              ! SP: Compute the cfac only once here since the same are use in both hamwan2bloch and dmewan2bloch
              ! + optimize the 2\pi r\cdot k with Blas
              IF ( (nkf1 >0) .AND. (nkf2 > 0) .AND. (nkf3 > 0) .AND. &
                 (nqf1 > 0) .AND. (nqf2 > 0) .AND. (nqf3 > 0) .AND. .NOT. mp_mesh_k) THEN          
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
                  cfac(ir) = ( tablex(irvec(1,ir)+2*nk1+1,xkk1) *&
                          tabley(irvec(2,ir)+2*nk2+1,xkk2) * tablez(irvec(3,ir)+2*nk3+1,xkk3) ) / ndegen_k(ir)
                  cfacq(ir) = ( tableqx(irvec(1,ir)+2*nk1+1,xkq1) *&
                          tableqy(irvec(2,ir)+2*nk2+1,xkq2) * tableqz(irvec(3,ir)+2*nk3+1,xkq3) ) /  ndegen_k(ir)
                ENDDO
              ELSE
                CALL dgemv('t', 3, nrr_k, twopi, irvec_r, 3, xkk, 1, 0.0_DP, rdotk, 1 )
                cfac(:) = exp( ci*rdotk(:) ) / ndegen_k(:)
                CALL dgemv('t', 3, nrr_k, twopi, irvec_r, 3, xkq, 1, 0.0_DP, rdotk, 1 )
                cfacq(:) = exp( ci*rdotk(:) ) / ndegen_k(:)
              ENDIF
              !
              ! ------------------------------------------------------        
              ! hamiltonian : Wannier -> Bloch 
              ! ------------------------------------------------------
              !
              ! Kohn-Sham first, then get the rotation matricies for following interp.
              IF (eig_read) THEN
                 CALL hamwan2bloch &
                   ( nbndsub, nrr_k, cufkk(:,:,ik), etf_ks (:, ikk), chw_ks, cfac)
                 CALL hamwan2bloch &
                   ( nbndsub, nrr_k, cufkq(:,:,ik), etf_ks (:, ikq), chw_ks, cfacq)
              ENDIF
              !
              CALL hamwan2bloch &
                   ( nbndsub, nrr_k, cufkk(:,:,ik), etf (:, ikk), chw, cfac)
              CALL hamwan2bloch &
                   ( nbndsub, nrr_k, cufkq(:,:,ik), etf (:, ikq), chw, cfacq)
              !
              ! ------------------------------------------------------        
              !  dipole: Wannier -> Bloch
              ! ------------------------------------------------------        
              !
              CALL dmewan2bloch &
                   ( nbndsub, nrr_k, irvec, ndegen_k, xkk, cufkk(:,:,ik), dmef (:,:,:, ikk), etf(:,ikk), etf_ks(:,ikk), cfac)
              CALL dmewan2bloch &
                   ( nbndsub, nrr_k, irvec, ndegen_k, xkq, cufkq(:,:,ik), dmef (:,:,:, ikq), etf(:,ikq), etf_ks(:,ikq), cfacq)
              !
              ! ------------------------------------------------------        
              !  velocity: Wannier -> Bloch
              ! ------------------------------------------------------        
              !
              IF (vme) THEN
                 IF (eig_read) THEN
                    CALL vmewan2bloch( nbndsub, nrr_k, irvec, ndegen_k, xkk, &
                                       cufkk(:,:,ik), vmef (:,:,:, ikk), etf(:,ikk), etf_ks(:,ikk), chw_ks)
                    CALL vmewan2bloch( nbndsub, nrr_k, irvec, ndegen_k, xkq, &
                                       cufkq(:,:,ik), vmef (:,:,:, ikq), etf(:,ikq), etf_ks(:,ikq), chw_ks)
                 ELSE
                    CALL vmewan2bloch( nbndsub, nrr_k, irvec, ndegen_k, xkk, &
                                       cufkk(:,:,ik), vmef (:,:,:, ikk), etf(:,ikk), etf_ks(:,ikk), chw)
                    CALL vmewan2bloch( nbndsub, nrr_k, irvec, ndegen_k, xkq, &
                                       cufkq(:,:,ik), vmef (:,:,:, ikq), etf(:,ikq), etf_ks(:,ikq), chw)
                 ENDIF
              ENDIF
              !
            ENDIF ! imode  
            !
            ! interpolate ONLY when (k,k+q) both have at least one band 
            ! within a Fermi shell of size fsthick 
            !
            IF ( (( minval ( abs(etf (:, ikk) - ef) ) < fsthick ) .and. &
                    ( minval ( abs(etf (:, ikq) - ef) ) < fsthick )) ) THEN
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
                 CALL ephwan2bloch_mem &
                   ( imode,nbndsub, nrr_k, irvec, ndegen_k, epmatwef, xkk, cufkk(:,:,ik), cufkq(:,:,ik), epmatf, nmodes )
                 !
               ENDIF
               !
               IF (lpolar) THEN
                 !
                 CALL compute_umn_f( nbndsub, cufkk(:,:,ik), cufkq(:,:,ik), bmatf )
                 !
                 IF ( (abs(xxq(1)) > eps) .or. (abs(xxq(2)) > eps) .or. (abs(xxq(3)) > eps) ) THEN
                   !      
                   CALL cryst_to_cart (1, xxq, bg, 1)
                   CALL rgd_blk_epw_fine_mem (imode, nq1, nq2, nq3, xxq, uf, epmatlrT(:,:,imode,ik), &
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
                   eptmp(ibnd-ibndmin+1,jbnd-ibndmin+1,imode,ik) = epmatf(ibnd,jbnd)        
                   !
                 ENDDO
               ENDDO
               ! 
            ENDIF
            !
         ENDDO  ! end loop over k points
       ENDDO ! mode
       !
       ! Now do the eigenvector rotation:
       ! epmatf(j) = sum_i eptmp(i) * uf(i,j)
       !
       DO ik=1, nkf
         CALL zgemm( 'n', 'n', (ibndmax-ibndmin+1) * (ibndmax-ibndmin+1), nmodes, nmodes, ( 1.d0, 0.d0 ), eptmp(:,:,:,ik),&
               (ibndmax-ibndmin+1) * (ibndmax-ibndmin+1), uf, nmodes, ( 0.d0, 0.d0 ), &
               epf17(:,:,:,ik), (ibndmax-ibndmin+1) * (ibndmax-ibndmin+1) )
       ENDDO
       ! 
       ! After the rotation, add the long-range that is already rotated
       DO jbnd = ibndmin, ibndmax
         DO ibnd = ibndmin, ibndmax
           epf17(ibnd-ibndmin+1,jbnd-ibndmin+1,:,:) = epf17(ibnd-ibndmin+1,jbnd-ibndmin+1,:,:) + epmatlrT(ibnd,jbnd,:,:)
         ENDDO
       ENDDO
       !
       IF (prtgkk)      CALL print_gkk( iq ) 
       IF (phonselfen ) CALL selfen_phon_q( iq )
       IF (elecselfen ) CALL selfen_elec_q( iq )
       IF (nest_fn    ) CALL nesting_fn_q( iq )
       IF (specfun    ) CALL spectral_func_q( iq )
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
  IF ( ALLOCATED(lambda_all) )    DEALLOCATE( lambda_all )
  IF ( ALLOCATED(gamma_all) )     DEALLOCATE( gamma_all )
  IF ( ALLOCATED(sigmai_all) )    DEALLOCATE( sigmai_all )
  IF ( ALLOCATED(sigmai_mode) )   DEALLOCATE( sigmai_mode )
  IF ( ALLOCATED(w2) )            DEALLOCATE( w2 )
  IF ( ALLOCATED(epf17) )            DEALLOCATE( epf17 )
  DEALLOCATE(cfac)
  DEALLOCATE(cfacq)
  DEALLOCATE(rdotk)
  DEALLOCATE(irvec_r)
  !
  CALL stop_clock ( 'ephwann' )
  !
  END SUBROUTINE ephwann_shuffle_mem
  !
