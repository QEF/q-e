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
  USE epwcom,        ONLY : nbndsub, lrepmatf, fsthick, epwread, longrange,     &
                            epwwrite, ngaussw, degaussw, lpolar, lifc,          &
                            nbndskip, parallel_k, parallel_q, etf_mem,          &
                            elecselfen, phonselfen, nest_fn, a2f,               &
                            vme, eig_read, ephwrite, nkf1, nkf2, nkf3,          & 
                            efermi_read, fermi_energy, specfun, band_plot,      &
                            nqf1, nqf2, nqf3, mp_mesh_k, restart
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
  !COMPLEX(kind=DP), ALLOCATABLE :: cfac1(:)
  !COMPLEX(kind=DP), ALLOCATABLE :: cfacq1(:)
  COMPLEX(kind=DP), ALLOCATABLE :: cfac(:)
  !! Used to store $e^{2\pi r \cdot k}$ exponential 
  COMPLEX(kind=DP), ALLOCATABLE :: cfacq(:)
  !! Used to store $e^{2\pi r \cdot k+q}$ exponential 
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
  IF ((.NOT. etf_mem) .AND. (ionode)) THEN
    ! open the .epmatwe file with the proper record length
    lrepmatw   = 2 * nbndsub * nbndsub * nrr_k * nmodes
    filint    = trim(prefix)//'.epmatwe'
    CALL diropn (iunepmatwe, 'epmatwe', lrepmatw, exst)  
    filint    = trim(prefix)//'.epmatwp'
    CALL diropn (iunepmatwp, 'epmatwp', lrepmatw, exst)
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
         IF (etf_mem) THEN 
           CALL ephbloch2wane &
             ( nbnd, nbndsub, nks, nkstot, xk, cu, cuq, &
             epmatq (:,:,:,imode,iq), nrr_k, irvec, wslen, epmatwe(:,:,:,imode,iq) )
         ELSE
           CALL ephbloch2wane &
             ( nbnd, nbndsub, nks, nkstot, xk, cu, cuq, &
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
           ( nbndsub, nmodes, xqc, nqc, irvec, nrr_k, nrr_q, epmatwe )
       ELSE
          CALL ephbloch2wanp_mem &
           ( nbndsub, nmodes, xqc, nqc, irvec, nrr_k, nrr_q, epmatwe_mem )
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
  !
  ! Need to be initialized
  epmatf(:,:,:) = czero
  ! allocate dipole matrix elements after getting grid size
  !
  ALLOCATE ( dmef(3, nbndsub, nbndsub, 2 * nkf) )
  IF (vme) ALLOCATE ( vmef(3, nbndsub, nbndsub, 2 * nkf) )
  !
  ALLOCATE(cfac(nrr_k))
  !DBSP
  !ALLOCATE(cfac1(nrr_k))
  !ALLOCATE(cfacq1(nrr_k))
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
    DO iq = 1, nqf
       !   
       CALL start_clock ( 'ep-interp' )
       !
       ! In case of big calculation, show progression of iq (especially usefull when
       ! elecselfen = true as nothing happen during the calculation otherwise. 
       !
       IF (.not. phonselfen) THEN 
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
       IF (.NOT. longrange) THEN
         CALL ephwan2blochp &
            ( nmodes, xxq, irvec, ndegen_q, nrr_q, uf, epmatwef, nbndsub, nrr_k )
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
          ! DBSP
          !write(900,*),'iq ik ',iq, ' ',ik
          !write(901,*),'iq ik ',iq, ' ',ik
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
            !DBSP
            !IF ( (iq == 1) .and. (ik ==4)) THEN
            !  CALL dgemv('t', 3, nrr_k, twopi, irvec_r, 3, xkk, 1, 0.0_DP, rdotk, 1 )
            !  cfac1(:) = exp( ci*rdotk(:) ) / ndegen_k(:)
            !  CALL dgemv('t', 3, nrr_k, twopi, irvec_r, 3, xkq, 1, 0.0_DP, rdotk, 1 )
            !  cfacq1(:) = exp( ci*rdotk(:) ) / ndegen_k(:)
            !ENDIF
          ELSE
            CALL dgemv('t', 3, nrr_k, twopi, irvec_r, 3, xkk, 1, 0.0_DP, rdotk, 1 )
            cfac(:) = exp( ci*rdotk(:) ) / ndegen_k(:)
            CALL dgemv('t', 3, nrr_k, twopi, irvec_r, 3, xkq, 1, 0.0_DP, rdotk, 1 )
            cfacq(:) = exp( ci*rdotk(:) ) / ndegen_k(:)
          ENDIF
          !DO ir = 1, nrr_k
          !  write(900,*)ir, ' ', cfacq(ir)
          !  write(901,*)ir, ' ', cfacq1(ir)
          !ENDDO
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
          !DBSP
          !write(900,*)'ikk ',etf(:,ikk)
          !write(900,*)'ikq ',etf(:,ikq)
          !
          ! ------------------------------------------------------        
          !  dipole: Wannier -> Bloch
          ! ------------------------------------------------------        
          !
          CALL dmewan2bloch &
               ( nbndsub, nrr_k, irvec, ndegen_k, xkk, cufkk, dmef (:,:,:, ikk), etf(:,ikk), etf_ks(:,ikk), cfac)
          CALL dmewan2bloch &
               ( nbndsub, nrr_k, irvec, ndegen_k, xkq, cufkq, dmef (:,:,:, ikq), etf(:,ikq), etf_ks(:,ikq), cfacq)
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
               CALL ephwan2bloch &
                 ( nbndsub, nrr_k, irvec, ndegen_k, epmatwef, xkk, cufkk, cufkq, epmatf, nmodes )
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
                 epf17(ibnd-ibndmin+1,jbnd-ibndmin+1,:,ik) = epmatf(ibnd,jbnd,:)        
                 !
               ENDDO
             ENDDO
             ! 
             !DBSP
             ! Debug on the long/short range. Usefull to keep commented for now.  
             !if (ik==2) then
             !  !print*,'iq ',iq
             !  !do imode = 1, nmodes
             !    !write(*,*) 'epmatf ',SUM((REAL(REAL(epmatf(:,:,imode))))**2)+SUM((REAL(AIMAG(epmatf(:,:,imode))))**2)
             !    F = SUM((REAL(REAL(epmatf(:,:,:))))**2)+SUM((REAL(AIMAG(epmatf(:,:,:))))**2)
             !    !S = SUM((REAL(REAL(epmatfs(:,:,:))))**2)+SUM((REAL(AIMAG(epmatfs(:,:,:))))**2)
             !    S = SUM((epmatfs(:,:,:))**2)
             !    L = SUM((REAL(REAL(epmatfl(:,:,:))))**2)+SUM((REAL(AIMAG(epmatfl(:,:,:))))**2)
             !    write(*,*) 'F, S+L', F, S+L
             !    DO ibnd = 1, nbndsub
             !      print*,'ibnd ',ibnd
             !      DO jbnd = 1, nbndsub
             !        print*,'jbnd ',jbnd
             !        F = SUM((REAL(REAL(epmatf(ibnd,jbnd,:))))**2)+SUM((REAL(AIMAG(epmatf(ibnd,jbnd,:))))**2)
             !        !S = SUM((REAL(REAL(epmatfs(ibnd,jbnd,:))))**2)+SUM((REAL(AIMAG(epmatfs(ibnd,jbnd,:))))**2)
             !        S = SUM(epmatfs(ibnd,jbnd,:)**2)
             !        L = SUM((REAL(REAL(epmatfl(ibnd,jbnd,:))))**2)+SUM((REAL(AIMAG(epmatfl(ibnd,jbnd,:))))**2)
             !        write(*,*) 'F, S+L', F, S+L
             !        DO imode = 1, nmodes
             !          print*,'imode ',imode
             !          F = (REAL(REAL(epmatf(ibnd,jbnd,imode))))**2+(REAL(AIMAG(epmatf(ibnd,jbnd,imode))))**2
             !          !S = (REAL(REAL(epmatfs(ibnd,jbnd,imode))))**2+(REAL(AIMAG(epmatfs(ibnd,jbnd,imode))))**2
             !          S = (epmatfs(ibnd,jbnd,imode))**2
             !          L = (REAL(REAL(epmatfl(ibnd,jbnd,imode))))**2+(REAL(AIMAG(epmatfl(ibnd,jbnd,imode))))**2
             !          write(*,*) 'F, S+L', F, S+L
             !        ENDDO
             !      ENDDO
             !    ENDDO
             !  !enddo
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
    CALL mem_size(ibndmin, ibndmax, nmodes, nkf)
    !
    ! Fine mesh set of g-matrices.  It is large for memory storage
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
        ! 
        ! SP: Compute the cfac only once here since the same are use in both hamwan2bloch and dmewan2bloch
        ! + optimize the 2\pi r\cdot k with Blas
        CALL dgemv('t', 3, nrr_k, twopi, irvec_r, 3, xkk, 1, 0.0_DP, rdotk, 1 )
        cfac(:) = exp( ci*rdotk ) / ndegen_k(:)
        CALL dgemv('t', 3, nrr_k, twopi, irvec_r, 3, xkq, 1, 0.0_DP, rdotk, 1 )
        cfacq(:) = exp( ci*rdotk ) / ndegen_k(:)
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
             ( nbndsub, nrr_k, irvec, ndegen_k, xkk, cufkk, dmef (:,:,:, ikk), etf(:,ikk), etf_ks(:,ikk), cfac)
        CALL dmewan2bloch &
             ( nbndsub, nrr_k, irvec, ndegen_k, xkq, cufkq, dmef (:,:,:, ikq), etf(:,ikq), etf_ks(:,ikq), cfac)
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
               ( nbndsub, nrr_k, irvec, ndegen_k, epmatwef, xkk, cufkk, cufkq, epmatf, nmodes )
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
  IF ( ALLOCATED(w2) )   DEALLOCATE( w2 )
  DEALLOCATE(cfac)
  DEALLOCATE(cfacq)
  DEALLOCATE(rdotk)
  DEALLOCATE(irvec_r)

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
  INTEGER             :: ibnd, jbnd, jmode, imode, irk, irq, ipol, i, lrepmatw
  character (len=256) :: filint
  complex(kind=DP)    :: aux ( nbndsub*nbndsub*nrr_k*nmodes*nrr_q ) 
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
    IF (etf_mem) THEN
      ! SP: The call to aux is now inside the loop
      !     This is important as otherwise the lrepmatw integer 
      !     could become too large for integer(kind=4).
      !     Note that in Fortran the record length has to be a integer
      !     of kind 4. 
      lrepmatw   = 2 * nbndsub * nbndsub * nrr_k * nmodes
      filint    = trim(prefix)//'.epmatwp'
      CALL diropn (iunepmatwp, 'epmatwp', lrepmatw, exst)
      DO irq = 1, nrr_q
        i = 0      
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
        CALL davcio ( aux, lrepmatw, iunepmatwp, irq, +1 )
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
  USE kinds,     ONLY : DP
  USE epwcom,    ONLY : nbndsub, vme, eig_read, etf_mem, lifc
  USE pwcom,     ONLY : ef
  USE elph2,     ONLY : nrr_k, nrr_q, chw, rdw, ifc, epmatwp, &
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
  character (len=256) :: filint
  INTEGER             :: ibnd, jbnd, jmode, imode, irk, irq, &
                         ipol, ios, i, lrepmatw
  complex(kind=DP)    :: aux ( nbndsub*nbndsub*nrr_k*nmodes*nrr_q )
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
    OPEN(unit=iundmedata,file='dmedata.fmt',status='old',iostat=ios)
    IF (eig_read) OPEN(unit=iunksdata,file='ksdata.fmt',status='old',iostat=ios)
    IF (vme) OPEN(unit=iunvmedata,file='vmedata.fmt',status='old',iostat=ios)
    IF (ios /= 0) call errore ('ephwann_shuffle', 'error opening epwdata.fmt',epwdata)
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
  IF (.not. ALLOCATED(epmatwp)) ALLOCATE ( epmatwp ( nbndsub, nbndsub, nrr_k, nmodes, nrr_q) )
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
  IF (etf_mem) then
    epmatwp = czero
    IF (mpime.eq.ionode_id) THEN
      ! SP: The call to aux is now inside the loop
      !     This is important as otherwise the lrepmatw integer 
      !     could become too large for integer(kind=4).
      !     Note that in Fortran the record length has to be a integer
      !     of kind 4.      
      lrepmatw   = 2 * nbndsub * nbndsub * nrr_k * nmodes
      filint    = trim(prefix)//'.epmatwp'
      CALL diropn (iunepmatwp, 'epmatwp', lrepmatw, exst)
      DO irq = 1, nrr_q
        i = 0     
        CALL davcio ( aux, lrepmatw, iunepmatwp, irq, -1 )
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
    ENDIF
    !
    CALL mp_bcast (epmatwp, ionode_id, inter_pool_comm)
    CALL mp_bcast (epmatwp, root_pool, intra_pool_comm)
    !
  ENDIF
  !
  CALL mp_barrier(inter_pool_comm)
  IF (mpime.eq.ionode_id) THEN
    CLOSE(iunepmatwp)
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
        

