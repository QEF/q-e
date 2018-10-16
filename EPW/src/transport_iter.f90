  !
  ! Copyright (C) 2016-2019 Samuel Ponce', Roxana Margine, Feliciano Giustino
  ! Copyright (C) 2010-2016 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino  
  ! 
  ! This file is distributed under the terms of the GNU General Public         
  ! License. See the file `LICENSE' in the root directory of the               
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .
  !
  !----------------------------------------------------------------------
  MODULE transport_iter
  !----------------------------------------------------------------------
  !! 
  !! This module contains all the subroutine linked with self-consistent electronic transport  
  !! 
  IMPLICIT NONE
  ! 
  CONTAINS
    ! 
    !-----------------------------------------------------------------------
    SUBROUTINE ibte( nind, etf_all, vkk_all, wkf_all, trans_prob, ef0, &
                     sparse_q, sparse_k, sparse_i, sparse_j, sparse_t ) 
    !-----------------------------------------------------------------------
    !!
    !!  This subroutine computes the scattering rate with the iterative BTE
    !!  (inv_tau).
    !!  The fine k-point and q-point grid have to be commensurate. 
    !!  The k-point grid uses crystal symmetry to decrease computational cost.
    !!
    !-----------------------------------------------------------------------
    USE kinds,         ONLY : DP
    USE io_global,     ONLY : stdout
    USE cell_base,     ONLY : alat, at, omega, bg
    USE phcom,         ONLY : nmodes
    USE epwcom,        ONLY : fsthick, mob_maxiter,  & 
                              eps_acustic, degaussw, nstemp, & 
                              system_2d, int_mob, ncarrier, restart, restart_freq,&
                              mp_mesh_k, nkf1, nkf2, nkf3, vme, broyden_beta
    USE pwcom,         ONLY : ef 
    USE elph2,         ONLY : ibndmax, ibndmin, etf, nkqf, nkf, wkf, dmef, vmef, & 
                              wf, xkf, epf17, nqtotf, nkqtotf, & 
                              map_rebal, xqf, wqf, nqf
    USE transportcom,  ONLY : transp_temp, mobilityh_save, mobilityel_save, lower_bnd, &
                              ixkqf_tr, s_BZtoIBZ_full
    USE constants_epw, ONLY : zero, one, two, pi, kelvin2eV, ryd2ev, & 
                              electron_SI, bohr2ang, ang2cm, hbarJ, eps6, eps8, eps10, &
                              eps2, eps4, eps80, eps160
    USE mp,            ONLY : mp_barrier, mp_sum, mp_bcast
    USE mp_global,     ONLY : inter_pool_comm, world_comm
    USE mp_world,      ONLY : mpime
    USE io_global,     ONLY : ionode_id
    USE symm_base,     ONLY : s, t_rev, time_reversal, set_sym_bl, nrot
    USE io_eliashberg, ONLY : kpmq_map
    USE printing,      ONLY : print_serta, print_serta_sym, print_mob, print_mob_sym
    USE io_scattering, ONLY : Fin_write, Fin_read
    USE noncollin_module, ONLY : noncolin
    USE io_files,      ONLY : diropn
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN) :: nind
    !! Total number of elements per cpu
    INTEGER, INTENT(IN) :: sparse_q(nind)
    !! Q-point mapping index
    INTEGER, INTENT(IN) :: sparse_k(nind)
    !! K-point mapping index
    INTEGER, INTENT(IN) :: sparse_i(nind)
    !! Band mapping index
    INTEGER, INTENT(IN) :: sparse_j(nind)
    !! Band mapping index
    INTEGER, INTENT(IN) :: sparse_t(nind)
    !! Temperature mapping index
    REAL(KIND=DP), INTENT(IN) :: etf_all(ibndmax-ibndmin+1,nkqtotf/2)
    !! Eigenenergies
    REAL(KIND=DP), INTENT(IN) :: vkk_all(3,ibndmax-ibndmin+1,nkqtotf/2)
    !! Velocity of k
    REAL(KIND=DP), INTENT(IN) :: wkf_all(nkqtotf/2)
    !! Weight of k
    REAL(KIND=DP), INTENT(IN) :: trans_prob(nind)
    !! Transition probability
    REAL(KIND=DP), INTENT(IN) :: ef0(nstemp)
    !! The Fermi level 
    ! 
    ! Local variables
    LOGICAL :: exst
    INTEGER :: iter, ind, ierr
    !! Iteration number in the IBTE
    INTEGER :: i, iiq, iq
    !! Cartesian direction index 
    INTEGER :: j
    !! Cartesian direction index 
    INTEGER :: ij
    !! Cartesian direction index 
    INTEGER :: ik
    !! K-point index
    INTEGER :: ikk
    !! Odd index to read etf
    INTEGER :: ikq
    !! Even k+q index to read etf
    INTEGER :: ibnd
    !! Local band index
    INTEGER :: jbnd
    !! Local band index
    INTEGER :: imode
    !! Local mode index
    INTEGER :: itemp
    !! Temperature index
    INTEGER :: ipool
    !! Index of the pool
    INTEGER :: nkq
    !! Index of the pool the the k+q point is
    INTEGER :: nkq_abs
    !! Index of the k+q point from the full grid. 
    INTEGER :: nkqtotf_tmp
    !! Temporary k-q points.
    INTEGER :: ikbz
    !! k-point index that run on the full BZ
    INTEGER :: nb
    !! Number of points in the BZ corresponding to a point in IBZ 
    INTEGER :: BZtoIBZ_tmp(nkf1*nkf2*nkf3)
    !! Temporary mapping
    INTEGER :: BZtoIBZ(nkf1*nkf2*nkf3)
    !! BZ to IBZ mapping
    INTEGER :: s_BZtoIBZ(3,3,nkf1*nkf2*nkf3)
    !! symmetry 
    INTEGER :: n
    !! Use for averaging
    ! 
    REAL(KIND=DP) :: tau
    !! Relaxation time
    REAL(KIND=DP) :: ekk
    !! Energy relative to Fermi level: $$\varepsilon_{n\mathbf{k}}-\varepsilon_F$$
    REAL(KIND=DP) :: ekq
    !! Energy relative to Fermi level: $$\varepsilon_{m\mathbf{k+q}}-\varepsilon_F$$
    REAL(KIND=DP) :: vkk(3,ibndmax-ibndmin+1)
    !! Electronic velocity $$v_{n\mathbf{k}}$$
    REAL(kind=DP) :: xkf_all(3,nkqtotf)
    !! Collect k-point coordinate (and k+q) from all pools in parallel case
    REAL(kind=DP) :: F_SERTA(3, ibndmax-ibndmin+1, nkqtotf/2, nstemp)
    !! SERTA solution
    REAL(kind=DP) :: F_in(3, ibndmax-ibndmin+1, nkqtotf/2, nstemp)
    !! In solution for iteration i
    REAL(kind=DP) :: F_out(3, ibndmax-ibndmin+1, nkqtotf/2, nstemp)
    !! In solution for iteration i
    REAL(kind=DP) :: F_rot(3)
    !! Rotated Fi_in by the symmetry operation 
    REAL(kind=DP) :: error(nstemp)
    !! Error in the hole mobility
    REAL(kind=DP) :: av_mob_old(nstemp)
    !! Average hole mobility from previous iteration
    REAL(kind=DP) :: av_mob(nstemp)
    !! Average hole mobility
    REAL(kind=DP) :: tmp(ibndmax-ibndmin+1, nkqtotf/2, nstemp)
    !REAL(kind=DP) :: tmp2(ibndmax-ibndmin+1, ibndmax-ibndmin+1,nstemp, nkqtotf/2, nqf)
    REAL(kind=DP) :: tmp2
    !! Used for the averaging
    REAL(kind=DP) :: tmp3(ibndmax-ibndmin+1)
    !! Used for the averaging
    REAL(kind=DP) :: ekk2
    !! Use for averaging
  
    !
    REAL(kind=DP) :: xkf_tmp (3, nkqtotf)
    !! Temporary k-point coordinate (dummy variable)
    REAL(kind=DP) :: wkf_tmp(nkqtotf)
    !! Temporary k-weights (dummy variable)
    ! 
    ! Gather all the k-point coordinate from all the pools
    xkf_all(:,:) = zero
    av_mob(:)    = zero
#ifdef __MPI
    CALL poolgather2 ( 3, nkqtotf, nkqf, xkf, xkf_all)
#else
    xkf_all = xkf
#endif
    !
   ! print*,'mp_mesh_k ',mp_mesh_k
   ! print*,'mpime ',mpime
   ! print*,'nind ',nind
   ! print*,'allocated ',ALLOCATED(ixkqf_tr)
   ! print*,'allocated s_BZtoIBZ_full',ALLOCATED(s_BZtoIBZ_full)
   ! FLUSH(6)
  
    ! Deal with symmetries
    IF (mp_mesh_k) THEN
      ALLOCATE(ixkqf_tr(nind), STAT=ierr)
      ALLOCATE(s_BZtoIBZ_full(3,3,nind), STAT=ierr)
      BZtoIBZ(:) = 0
      s_BZtoIBZ(:,:,:) = 0
      ixkqf_tr(:) = 0
      !call move_alloc(test1, s_BZtoIBZ_full)
      s_BZtoIBZ_full(:,:,:) = 0
      ! 
      IF ( mpime .eq. ionode_id ) THEN
        ! 
        CALL set_sym_bl( )
        !
        ! What we get from this call is BZtoIBZ
        CALL kpoint_grid_epw ( nrot, time_reversal, .false., s, t_rev, bg, nkf1*nkf2*nkf3, &
                   nkf1,nkf2,nkf3, nkqtotf_tmp, xkf_tmp, wkf_tmp,BZtoIBZ,s_BZtoIBZ)
        ! 
        BZtoIBZ_tmp(:) = 0
        DO ikbz=1, nkf1*nkf2*nkf3
          BZtoIBZ_tmp(ikbz) = map_rebal( BZtoIBZ( ikbz ) )
        ENDDO
        BZtoIBZ(:) = BZtoIBZ_tmp(:)
        ! 
      ENDIF ! mpime
      CALL mp_bcast( s_BZtoIBZ, ionode_id, inter_pool_comm )
      CALL mp_bcast( BZtoIBZ, ionode_id, inter_pool_comm )
      ! 
      DO ind=1, nind
        iq    = sparse_q( ind )
        ik    = sparse_k( ind )
        !print*,'ind ik ',ind, ik
        ! 
        CALL kpmq_map( xkf_all(:, 2*ik-1 ), xqf (:, iq), +1, nkq_abs )
        s_BZtoIBZ_full(:,:,ind) = s_BZtoIBZ(:,:,nkq_abs)
        ixkqf_tr(ind) = BZtoIBZ(nkq_abs)
        !print*,'ind iq ik ixkqf_tr ',ind, iq, ik, ixkqf_tr(ind), s_BZtoIBZ_full(1,1,ind)
      ENDDO
      ! 
    ENDIF
    !
    ! First computes the SERTA solution as the first step of the IBTE
    F_SERTA(:,:,:,:) = zero
    tmp(:,:,:) = zero
    !tmp2(:,:,:,:,:) = zero
    ! 
    DO ind=1, nind
      iq    = sparse_q( ind )
      ik    = sparse_k( ind )
      ibnd  = sparse_i( ind )
      jbnd  = sparse_j( ind )
      itemp = sparse_t( ind )
      ! 
      tmp(ibnd, ik, itemp) = tmp(ibnd, ik, itemp)  + trans_prob(ind)
      !tmp2(jbnd, ibnd, itemp, ik, iq)  = trans_prob(ind)
  
      !IF (ik==2 .and. ibnd ==2 .and. itemp ==2) print*,'ind tmp ', ind, tmp(ibnd, ik, itemp)
  
      !IF (ik==2) print*,ind, trans_prob(ind)
      !print*,'ind iq ik ibnd jbnd itemp ',ind, iq, ik, ibnd, jbnd, itemp
      !print*,'tmp ',tmp(ibnd, ik, itemp)
    ENDDO
    !print*,'ind=10, iq==1, ik==2, ibnd=2, jbnd=2, itemp==1 ', trans_prob(10)
    ! 
    CALL mp_sum(tmp, world_comm)
    ! 
    ! Average over degenerate eigenstates:
    WRITE(stdout,'(5x,"Average over degenerate eigenstates is performed")')
    ! 
    tmp3(:) = zero
    DO itemp=1, nstemp
      DO ik = 1, nkqtotf/2
        ! 
        DO ibnd = 1, ibndmax-ibndmin+1
          ekk = etf_all (ibndmin-1+ibnd, ik)
          n = 0
          tmp2 = 0.0_DP
          DO jbnd = 1, ibndmax-ibndmin+1
            ekk2 = etf_all (ibndmin-1+jbnd, ik)
            IF ( ABS(ekk2-ekk) < eps6 ) THEN
              n = n + 1
              tmp2 =  tmp2 + tmp(ibnd,ik,itemp)
            ENDIF
            ! 
          ENDDO ! jbnd
          tmp3(ibnd) = tmp2 / float(n)
          !
        ENDDO ! ibnd
         tmp(:,ik,itemp) = tmp3(:)
        ! 
      ENDDO ! nkqtotf  
    ENDDO ! itemp
    ! 
    !
    DO itemp=1, nstemp 
      DO ik=1, nkqtotf/2
        DO ibnd=1, ibndmax-ibndmin+1
          IF ( ABS(tmp(ibnd, ik, itemp)) > eps160 ) THEN
            F_SERTA(:, ibnd, ik, itemp) = vkk_all(:,ibnd,ik) / ( two * tmp(ibnd,ik,itemp) )  
          ENDIF
        ENDDO
        !IF (itemp==2) print*,'ik ',ik, SUM(F_SERTA(:,:,ik,2)), SUM(vkk_all(:,:,ik))
      ENDDO
    ENDDO
    !  
    !print*,'F_SERTA ',SUM(F_SERTA)
    !print*,'trans_prob ',trans_prob(1:5)
    !print*,'F_SERTA before ',F_SERTA(:,1,1,1)
    ! Now compute and print the electron and hole mobility of SERTA
    IF (mp_mesh_k) THEN
      ! Use k-point symmetry
      CALL print_serta_sym(F_SERTA, BZtoIBZ, s_BZtoIBZ, vkk_all, etf_all, wkf_all, ef0)
    ELSE 
      ! No symmetry
      CALL print_serta(F_SERTA, vkk_all, etf_all, wkf_all, ef0)
    ENDIF
    !STOP
    !print*,'F_SERTA after ',F_SERTA(:,1,1,1)
    ! 
    ! NOW solve IBTE
  
    ! Read from file
    iter = 1
    F_in(:,:,:,:) = zero
    IF (ncarrier > 1E5) THEN
      CALL Fin_read(iter, F_in, av_mob_old, .TRUE.)
    ENDIF
    ! 
    IF (ncarrier < -1E5) THEN
      CALL Fin_read(iter, F_in, av_mob_old, .FALSE.)
    ENDIF
    !write(*,*)'F_in ',sum(F_in)
    !write(*,*)'av_mob_old ',av_mob_old
    ! 
    ! If it is the first time, put to SERTA
    IF (iter == 1) THEN
      F_in(:,:,:,:) = F_SERTA(:,:,:,:)
      av_mob_old(:) = 0.0
    ENDIF
    ! 
    F_out(:,:,:,:) = zero
    error(:) = 1000
    ! 
    ! Now compute the Iterative solution for electron or hole
    WRITE(stdout,'(5x,a)') ' '
    WRITE(stdout,'(5x,a)') repeat('=',67)
    WRITE(stdout,'(5x,"Start solving iterative Boltzmann Transport Equation")')
    WRITE(stdout,'(5x,a/)') repeat('=',67)
    !  
    DO WHILE (MAXVAL(error) > eps6)  
      WRITE(stdout,'(/5x,"Iteration number:", i10," "/)') iter
      ! 
      IF (iter > mob_maxiter) THEN
        WRITE(stdout,'(5x,a)') repeat('=',67)
        WRITE(stdout,'(5x,"The iteration reached the maximum but did not converge.")')
        WRITE(stdout,'(5x,a/)') repeat('=',67)
        exit
      ENDIF
      ! 
      IF (mp_mesh_k) THEN
        ! Use k-point symmetry
        DO ind=1, nind
          !  
          F_rot(:) = zero
          iq    = sparse_q( ind )
          ik    = sparse_k( ind )
          ibnd  = sparse_i( ind )
          jbnd  = sparse_j( ind )
          itemp = sparse_t( ind )
          ! 
          !print*,'before F_in ',ind,F_in(:, jbnd, ixkqf_tr(ind), itemp)
          CALL cryst_to_cart(1,F_in(:, jbnd, ixkqf_tr(ind), itemp), at, -1)
          !print*,'after F_in ',ind,F_in(:, jbnd, ixkqf_tr(ind), itemp)
          CALL dgemv( 'n', 3, 3, 1.d0,&
             REAL(s_BZtoIBZ_full(:,:,ind), kind=DP), 3, F_in(:, jbnd, ixkqf_tr(ind), itemp),1 ,0.d0 , F_rot(:), 1 )
          !print*,'before F_rot ',ind, F_rot
          CALL cryst_to_cart(1, F_in(:, jbnd, ixkqf_tr(ind), itemp), bg, 1)
          CALL cryst_to_cart(1,F_rot,bg,1)
     
          !print*,'after F_rot ',ind,F_rot(:) 
          F_out(:, ibnd, ik, itemp) = F_out(:, ibnd, ik, itemp) + two * trans_prob(ind) * F_rot(:)
          ! 
        ENDDO
      ELSE
        DO ind=1, nind
          !  
          iq    = sparse_q( ind )
          ik    = sparse_k( ind )
          ibnd  = sparse_i( ind )
          jbnd  = sparse_j( ind )
          itemp = sparse_t( ind )
          ! We need F_in at k+q point
          CALL kpmq_map( xkf_all(:, 2*ik-1 ), xqf (:, iq), +1, nkq_abs )  
          ! 
          F_out(:, ibnd, ik, itemp) = F_out(:, ibnd, ik, itemp) + two * trans_prob(ind) * F_in(:, jbnd, nkq_abs, itemp)
          !  
        ENDDO
      ENDIF
      ! 
      CALL mp_sum(F_out, world_comm)  
      ! 
      DO itemp = 1, nstemp
        DO ik = 1, nkqtotf/2
          DO ibnd = 1, ibndmax-ibndmin+1
            IF ( ABS(tmp(ibnd, ik, itemp)) > eps160 ) THEN
              F_out(:, ibnd, ik, itemp) = F_SERTA(:, ibnd, ik, itemp) +&
                               F_out(:, ibnd, ik, itemp) / ( two * tmp(ibnd, ik, itemp) )
            ENDIF
          ENDDO
        ENDDO
      ENDDO
      !  
      IF (mp_mesh_k) THEN
        CALL print_mob_sym(F_out, BZtoIBZ, s_BZtoIBZ, vkk_all, etf_all, wkf_all, ef0, av_mob) 
      ELSE 
        CALL print_mob(F_out, vkk_all, etf_all, wkf_all, ef0, av_mob) 
      ENDIF
      ! 
      ! Computes the error
      DO itemp = 1, nstemp
        error(itemp) = ABS( av_mob(itemp) - av_mob_old(itemp) ) 
      ENDDO
      av_mob_old = av_mob
      WRITE(stdout,'(a)')
      WRITE(stdout,'(45x, 1E18.6, a)') MAXVAL(error), '     Err'
      !
      ! Save F_in
      ! Full mixing 
      !F_in = F_out
      ! Linear mixing
      F_in = (1.0 - broyden_beta ) * F_in + broyden_beta * F_out 
      F_out = zero
      ! 
      iter = iter + 1
      ! 
      ! Save F_in to file:
      IF (ncarrier > 1E5) THEN 
        CALL Fin_write(iter, F_in, av_mob_old, .TRUE.) 
      ENDIF
      ! 
      IF (ncarrier < -1E5) THEN 
        CALL Fin_write(iter, F_in, av_mob_old, .FALSE.) 
      ENDIF
      ! 
      ! 
    ENDDO ! end of while loop
    ! 
    RETURN
    !
    ! ---------------------------------------------------------------------------
    END SUBROUTINE ibte
    !----------------------------------------------------------------------------
    !
    !----------------------------------------------------------------------------
    SUBROUTINE iter_restart(etf_all, wkf_all, vkk_all, ind_tot, ind_totcb, ef0, efcb)
    !----------------------------------------------------------------------------
    !  
    ! This subroutine opens all the required files to restart an IBTE calculation
    ! then call the ibte subroutine to perform the iterations. 
    ! This routine requires that the scattering rates have been computed previously. 
    !  
    ! ----------------------------------------------------------------------------
    USE kinds,            ONLY : DP, i4b
    USE elph2,            ONLY : nkqtotf, ibndmin, ibndmax
    USE mp_world,         ONLY : mpime, world_comm
    USE io_global,        ONLY : ionode_id, stdout
    USE io_files,         ONLY : tmp_dir, prefix
    USE epwcom,           ONLY : nstemp, ncarrier
    USE constants_epw,    ONLY : zero
    USE io_epw,           ONLY : iufilibtev_sup, iunepmat, iunsparseq, iunsparsek, &
                                 iunsparsei, iunsparsej, iunsparset, iunsparseqcb, &
                                 iunsparsekcb, iunrestart, iunsparseicb, iunsparsejcb,&
                                 iunsparsetcb, iunepmatcb
    USE transportcom,     ONLY : lower_bnd, upper_bnd
    USE mp,               ONLY : mp_bcast
    USE division,         ONLY : fkbounds2
#if defined(__MPI)
    USE parallel_include, ONLY : MPI_OFFSET, MPI_MODE_RDONLY, MPI_INFO_NULL, &
                                 MPI_SEEK_SET, MPI_DOUBLE_PRECISION, MPI_STATUS_IGNORE, &
                                 MPI_OFFSET_KIND, MPI_INTEGER4
#endif    
    !
    IMPLICIT NONE
    ! 
#if defined(__MPI)
    INTEGER (kind=MPI_OFFSET_KIND), INTENT(INOUT) :: ind_tot
    !! Total number of component for valence band
    INTEGER (kind=MPI_OFFSET_KIND), INTENT(INOUT) :: ind_totcb
    !! Total number of component for the conduction band
#else
    INTEGER(KIND=8), INTENT(INOUT) :: ind_tot
    !! Total number of component for valence band
    INTEGER(KIND=8), INTENT(INOUT) :: ind_totcb
    !! Total number of component for conduction band
#endif    
    !
    REAL(kind=DP), INTENT(INOUT) :: etf_all(ibndmax-ibndmin+1,nkqtotf/2)
    !! Eigen-energies on the fine grid collected from all pools in parallel case
    REAL(kind=DP), INTENT(INOUT) :: wkf_all(nkqtotf/2)
    !! k-point weights from all the cpu
    REAL(kind=DP), INTENT(INOUT) :: vkk_all(3,ibndmax-ibndmin+1,nkqtotf/2)
    !! velocity from all the k-points
    REAL(KIND=DP), INTENT(INOUT) :: ef0(nstemp)
    !! Fermi level for the temperature itemp     
    REAL(KIND=DP), INTENT(INOUT) :: efcb(nstemp)
    !! Fermi level for the temperature itemp for cb band    
    ! 
    ! Local variables
    !
    CHARACTER (len=256) :: filint
    !! Name of the file to write/read    
    ! 
    INTEGER :: ierr
    !! Error status
    INTEGER :: ik
    !! K-point
    INTEGER :: ios
    !! IO error message
    INTEGER :: itemp
    !! Temperature index
    INTEGER :: ibnd
    !! Counter on bandA
    INTEGER :: iktmp
    !! Dummy counter for k-points
    INTEGER :: ibtmp
    !! Dummy counter for bands
    INTEGER :: nind
    !! Number of local elements per cores. 
    INTEGER(kind=i4b), ALLOCATABLE :: sparse_q( : )
    !! Index mapping for q-points
    INTEGER(kind=i4b), ALLOCATABLE :: sparse_k( : )
    !! Index mapping for k-points
    INTEGER(kind=i4b), ALLOCATABLE :: sparse_i( : )
    !! Index mapping for i bands
    INTEGER(kind=i4b), ALLOCATABLE :: sparse_j( : )
    !! Index mapping for j bands
    INTEGER(kind=i4b), ALLOCATABLE :: sparse_t( : )
    !! Index mapping for temperature 
    INTEGER(kind=i4b), ALLOCATABLE :: sparsecb_q( : )
    !! Index mapping for q-points for cb
    INTEGER(kind=i4b), ALLOCATABLE :: sparsecb_k( : )
    !! Index mapping for k-points for cb
    INTEGER(kind=i4b), ALLOCATABLE :: sparsecb_i( : )
    !! Index mapping for i bands for cb
    INTEGER(kind=i4b), ALLOCATABLE :: sparsecb_j( : )
    !! Index mapping for j bands for cb
    INTEGER(kind=i4b), ALLOCATABLE :: sparsecb_t( : )
    !! Index mapping for temperature for cb
#if defined(__MPI)
    INTEGER (kind=MPI_OFFSET_KIND) :: lrepmatw2
    !! Local core offset for reading
    INTEGER (kind=MPI_OFFSET_KIND) :: lrepmatw4
    !! Local core offset for reading
    INTEGER (kind=MPI_OFFSET_KIND) :: lsize
    !! Offset to tell where to start reading the file
#else
    INTEGER (kind=8) :: lrepmatw2
    !! Local core offset for reading
    INTEGER (kind=i4b) :: lrepmatw4
    !! Local core offset for reading
    INTEGER (kind=8) :: lsize
    !! Offset to tell where to start reading the file
#endif
    ! 
    REAL(kind=DP) :: dum1
    !! Dummy variable
    REAL(kind=DP), ALLOCATABLE :: trans_prob(:)
    !! Transition probabilities
    REAL(kind=DP), ALLOCATABLE :: trans_probcb(:)
    !! Transition probabilities for cb    
    ! 
    etf_all(:,:)   = zero
    wkf_all(:)     = zero
    vkk_all(:,:,:) = zero
    ! 
    ! SP - The implementation only works with MPI so far
#ifdef __MPI
    ! Read velocities
    IF (mpime.eq.ionode_id) THEN
      !
      OPEN(unit=iufilibtev_sup,file='IBTEvel_sup.fmt',status='old',iostat=ios)
      READ(iufilibtev_sup,'(a)')
      READ(iufilibtev_sup,*) ind_tot, ind_totcb
      READ(iufilibtev_sup,'(a)')
      DO itemp=1, nstemp
        READ(iufilibtev_sup,*) dum1, ef0(itemp), efcb(itemp)
      ENDDO
      READ(iufilibtev_sup,'(a)')
      ! 
      DO ik = 1, nkqtotf/2
        DO ibnd = 1, ibndmax-ibndmin+1
          READ(iufilibtev_sup,*) iktmp, ibtmp, vkk_all(:,ibnd,ik), etf_all(ibnd,ik), wkf_all(ik)
        ENDDO
      ENDDO
      !  
    ENDIF
    ! 
    CALL MPI_BCAST( ind_tot, 1, MPI_OFFSET, ionode_id, world_comm, ierr)
    CALL MPI_BCAST( ind_totcb, 1, MPI_OFFSET, ionode_id, world_comm, ierr)
    CALL mp_bcast (ef0, ionode_id, world_comm)
    CALL mp_bcast (efcb, ionode_id, world_comm)
    CALL mp_bcast (vkk_all, ionode_id, world_comm)
    CALL mp_bcast (wkf_all, ionode_id, world_comm)
    CALL mp_bcast (etf_all, ionode_id, world_comm)
    ! 
    ! Now choose HOLE OR ELECTRON (the implementation does not support both)
    ! HOLE
    IF (ncarrier < -1E5) THEN    
      ! 
      ! Split all the matrix elements across all cores. 
      CALL fkbounds2( ind_tot, lower_bnd, upper_bnd )
      ! 
      ! Allocate the local size 
      nind = upper_bnd - lower_bnd + 1
      WRITE(stdout,'(5x,a,i10)') 'Number of elements per core ',nind
      ALLOCATE ( trans_prob ( nind ) )
      trans_prob(:) = 0.0d0
      ! 
      ! Open file containing trans_prob 
      filint = trim(tmp_dir)//trim(prefix)//'.epmatkq1'
      CALL MPI_FILE_OPEN(world_comm,filint,MPI_MODE_RDONLY, MPI_INFO_NULL, iunepmat, ierr)
      IF( ierr /= 0 ) CALL errore( 'iter_restart', 'error in MPI_FILE_OPEN X.epmatkq1',1 )
      !
      ! Offset depending on CPU
      lrepmatw2 = INT( lower_bnd -1, kind = MPI_OFFSET_KIND ) * 8_MPI_OFFSET_KIND
      ! 
      ! Size of what we read
      lsize = INT( nind , kind = MPI_OFFSET_KIND )
      !
      CALL MPI_FILE_SEEK(iunepmat, lrepmatw2, MPI_SEEK_SET, ierr)
      IF( ierr /= 0 ) CALL errore( 'iter_restart', 'error in MPI_FILE_SEEK',1 )
      CALL MPI_FILE_READ(iunepmat, trans_prob(:), lsize, MPI_DOUBLE_PRECISION, MPI_STATUS_IGNORE, ierr)
      IF( ierr /= 0 ) CALL errore( 'iter_restart', 'error in MPI_FILE_READ',1 )
      !      
      ! Now read the sparse matrix mapping
      CALL MPI_FILE_OPEN(world_comm,'sparseq',MPI_MODE_RDONLY ,MPI_INFO_NULL, iunsparseq, ierr)
      IF( ierr /= 0 ) CALL errore( 'iter_restart', 'error in MPI_FILE_OPEN sparseq',1 )
      CALL MPI_FILE_OPEN(world_comm,'sparsek',MPI_MODE_RDONLY ,MPI_INFO_NULL, iunsparsek, ierr)
      IF( ierr /= 0 ) CALL errore( 'iter_restart', 'error in MPI_FILE_OPEN sparsek',1 )
      CALL MPI_FILE_OPEN(world_comm,'sparsei',MPI_MODE_RDONLY ,MPI_INFO_NULL, iunsparsei, ierr)
      IF( ierr /= 0 ) CALL errore( 'iter_restart', 'error in MPI_FILE_OPEN sparsei',1 )
      CALL MPI_FILE_OPEN(world_comm,'sparsej',MPI_MODE_RDONLY ,MPI_INFO_NULL, iunsparsej, ierr)
      IF( ierr /= 0 ) CALL errore( 'iter_restart', 'error in MPI_FILE_OPEN sparsej',1 )
      CALL MPI_FILE_OPEN(world_comm,'sparset',MPI_MODE_RDONLY ,MPI_INFO_NULL, iunsparset, ierr)
      IF( ierr /= 0 ) CALL errore( 'iter_restart', 'error in MPI_FILE_OPEN sparset',1 )
      ! 
      ALLOCATE ( sparse_q ( nind ) )
      ALLOCATE ( sparse_k ( nind ) )
      ALLOCATE ( sparse_i ( nind ) )
      ALLOCATE ( sparse_j ( nind ) )
      ALLOCATE ( sparse_t ( nind ) )
      sparse_q(:) = 0.0d0
      sparse_k(:) = 0.0d0
      sparse_i(:) = 0.0d0
      sparse_j(:) = 0.0d0
      sparse_t(:) = 0.0d0
      !        
      lrepmatw4 = INT( lower_bnd - 1, kind = MPI_OFFSET_KIND ) * 4_MPI_OFFSET_KIND
      !
      CALL MPI_FILE_SEEK(iunsparseq, lrepmatw4, MPI_SEEK_SET, ierr)
      IF( ierr /= 0 ) CALL errore( 'iter_restart', 'error in MPI_FILE_SEEK',1 )
      CALL MPI_FILE_READ(iunsparseq, sparse_q(:), lsize, MPI_INTEGER4, MPI_STATUS_IGNORE, ierr)
      IF( ierr /= 0 ) CALL errore( 'iter_restart', 'error in MPI_FILE_READ',1 )
      CALL MPI_FILE_SEEK(iunsparsek, lrepmatw4, MPI_SEEK_SET, ierr)
      IF( ierr /= 0 ) CALL errore( 'iter_restart', 'error in MPI_FILE_SEEK',1 )
      CALL MPI_FILE_READ(iunsparsek, sparse_k(:), lsize, MPI_INTEGER4, MPI_STATUS_IGNORE, ierr)
      IF( ierr /= 0 ) CALL errore( 'iter_restart', 'error in MPI_FILE_READ',1 )
      CALL MPI_FILE_SEEK(iunsparsei, lrepmatw4, MPI_SEEK_SET, ierr)
      IF( ierr /= 0 ) CALL errore( 'iter_restart', 'error in MPI_FILE_SEEK',1 )
      CALL MPI_FILE_READ(iunsparsei, sparse_i(:), lsize, MPI_INTEGER4, MPI_STATUS_IGNORE, ierr)
      IF( ierr /= 0 ) CALL errore( 'iter_restart', 'error in MPI_FILE_READ',1 )
      CALL MPI_FILE_SEEK(iunsparsej, lrepmatw4, MPI_SEEK_SET, ierr)
      IF( ierr /= 0 ) CALL errore( 'iter_restart', 'error in MPI_FILE_SEEK',1 )
      CALL MPI_FILE_READ(iunsparsej, sparse_j(:), lsize, MPI_INTEGER4, MPI_STATUS_IGNORE, ierr)
      IF( ierr /= 0 ) CALL errore( 'iter_restart', 'error in MPI_FILE_READ',1 )
      CALL MPI_FILE_SEEK(iunsparset, lrepmatw4, MPI_SEEK_SET, ierr)
      IF( ierr /= 0 ) CALL errore( 'iter_restart', 'error in MPI_FILE_SEEK',1 )
      CALL MPI_FILE_READ(iunsparset, sparse_t(:), lsize, MPI_INTEGER4, MPI_STATUS_IGNORE, ierr)
      IF( ierr /= 0 ) CALL errore( 'iter_restart', 'error in MPI_FILE_READ',1 )
      ! 
      ! Now call the ibte to solve the BTE iteratively until convergence
      CALL ibte(nind, etf_all, vkk_all, wkf_all, trans_prob, ef0, sparse_q, sparse_k, sparse_i, sparse_j, sparse_t)
      ! 
      CALL MPI_FILE_CLOSE(iunepmat,ierr)
      IF( ierr /= 0 ) CALL errore( 'iter_restart', 'error in MPI_FILE_CLOSE',1)
      CALL MPI_FILE_CLOSE(iunsparseq,ierr)
      IF( ierr /= 0 ) CALL errore( 'iter_restart', 'error in MPI_FILE_CLOSE',1)
      CALL MPI_FILE_CLOSE(iunsparsek,ierr)
      IF( ierr /= 0 ) CALL errore( 'iter_restart', 'error in MPI_FILE_CLOSE',1)
      CALL MPI_FILE_CLOSE(iunsparsei,ierr)
      IF( ierr /= 0 ) CALL errore( 'iter_restart', 'error in MPI_FILE_CLOSE',1)
      CALL MPI_FILE_CLOSE(iunsparsej,ierr)
      IF( ierr /= 0 ) CALL errore( 'iter_restart', 'error in MPI_FILE_CLOSE',1)
      CALL MPI_FILE_CLOSE(iunsparset,ierr)
      IF( ierr /= 0 ) CALL errore( 'iter_restart', 'error in MPI_FILE_CLOSE',1)
      DEALLOCATE(trans_prob)
      DEALLOCATE(sparse_q)
      DEALLOCATE(sparse_k)
      DEALLOCATE(sparse_i)
      DEALLOCATE(sparse_j)
      DEALLOCATE(sparse_t) 
      ! 
    ENDIF
    ! Electrons
    IF (ncarrier > 1E5) THEN
      ! 
      CALL fkbounds2( ind_totcb, lower_bnd, upper_bnd )
      ! Allocate the local size 
      nind = upper_bnd - lower_bnd + 1
      WRITE(stdout,'(5x,a,i10)') 'Number of elements per core ',nind
      ALLOCATE ( trans_probcb ( nind ) )
      trans_probcb(:) = 0.0d0
      ! 
      ! Open file containing trans_prob 
      filint = trim(tmp_dir)//trim(prefix)//'.epmatkqcb1'
      CALL MPI_FILE_OPEN(world_comm, filint, MPI_MODE_RDONLY, MPI_INFO_NULL, iunepmatcb, ierr)
      IF( ierr /= 0 ) CALL errore( 'iter_restart', 'error in MPI_FILE_OPEN X.epmatkq1', 1 )
      !
      ! Offset depending on CPU
      lrepmatw2 = INT( lower_bnd-1, kind = MPI_OFFSET_KIND ) * 8_MPI_OFFSET_KIND
      ! 
      ! Size of what we read
      lsize = INT( nind, kind = MPI_OFFSET_KIND )
      !
      CALL MPI_FILE_SEEK(iunepmatcb, lrepmatw2, MPI_SEEK_SET, ierr)
      IF( ierr /= 0 ) CALL errore( 'iter_restart', 'error in MPI_FILE_SEEK iunepmatcb', 1 )
      CALL MPI_FILE_READ(iunepmatcb, trans_probcb(:), lsize, MPI_DOUBLE_PRECISION, MPI_STATUS_IGNORE, ierr)
      IF( ierr /= 0 ) CALL errore( 'iter_restart', 'error in MPI_FILE_READ iunepmatcb', 1 )
      !      
      ! Now read the sparse matrix mapping
      CALL MPI_FILE_OPEN(world_comm, 'sparseqcb', MPI_MODE_RDONLY, MPI_INFO_NULL, iunsparseqcb, ierr)
      IF( ierr /= 0 ) CALL errore( 'iter_restart', 'error in MPI_FILE_OPEN sparseqcb', 1 )
      CALL MPI_FILE_OPEN(world_comm, 'sparsekcb', MPI_MODE_RDONLY, MPI_INFO_NULL, iunsparsekcb, ierr)
      IF( ierr /= 0 ) CALL errore( 'iter_restart', 'error in MPI_FILE_OPEN sparsekcb', 1 )
      CALL MPI_FILE_OPEN(world_comm, 'sparseicb', MPI_MODE_RDONLY, MPI_INFO_NULL, iunsparseicb, ierr)
      IF( ierr /= 0 ) CALL errore( 'iter_restart', 'error in MPI_FILE_OPEN sparseicb', 1 )
      CALL MPI_FILE_OPEN(world_comm, 'sparsejcb', MPI_MODE_RDONLY, MPI_INFO_NULL, iunsparsejcb, ierr)
      IF( ierr /= 0 ) CALL errore( 'iter_restart', 'error in MPI_FILE_OPEN sparsejcb', 1 )
      CALL MPI_FILE_OPEN(world_comm, 'sparsetcb', MPI_MODE_RDONLY, MPI_INFO_NULL, iunsparsetcb, ierr)
      IF( ierr /= 0 ) CALL errore( 'iter_restart', 'error in MPI_FILE_OPEN sparsetcb', 1 )    
      ! 
      ALLOCATE ( sparsecb_q ( nind ) )
      ALLOCATE ( sparsecb_k ( nind ) )
      ALLOCATE ( sparsecb_i ( nind ) )
      ALLOCATE ( sparsecb_j ( nind ) )
      ALLOCATE ( sparsecb_t ( nind ) )
      sparsecb_q(:) = 0.0d0
      sparsecb_k(:) = 0.0d0
      sparsecb_i(:) = 0.0d0
      sparsecb_j(:) = 0.0d0
      sparsecb_t(:) = 0.0d0
      !        
      lrepmatw4 = INT( lower_bnd - 1, kind = MPI_OFFSET_KIND ) * 4_MPI_OFFSET_KIND
      !
      CALL MPI_FILE_SEEK(iunsparseqcb, lrepmatw4, MPI_SEEK_SET, ierr)
      IF( ierr /= 0 ) CALL errore( 'iter_restart', 'error in MPI_FILE_SEEK iunsparseqcb',1 )
      CALL MPI_FILE_READ(iunsparseqcb, sparsecb_q(:), lsize, MPI_INTEGER4, MPI_STATUS_IGNORE, ierr)
      IF( ierr /= 0 ) CALL errore( 'iter_restart', 'error in MPI_FILE_READ iunsparseqcb',1 )
      CALL MPI_FILE_SEEK(iunsparsekcb, lrepmatw4, MPI_SEEK_SET, ierr)
      IF( ierr /= 0 ) CALL errore( 'iter_restart', 'error in MPI_FILE_SEEK iunsparsekcb',1 )
      CALL MPI_FILE_READ(iunsparsekcb, sparsecb_k(:), lsize, MPI_INTEGER4, MPI_STATUS_IGNORE, ierr)
      IF( ierr /= 0 ) CALL errore( 'iter_restart', 'error in MPI_FILE_READ iunsparsekcb',1 )
      CALL MPI_FILE_SEEK(iunsparseicb, lrepmatw4, MPI_SEEK_SET, ierr)
      IF( ierr /= 0 ) CALL errore( 'iter_restart', 'error in MPI_FILE_SEEK iunsparseicb',1 )
      CALL MPI_FILE_READ(iunsparseicb, sparsecb_i(:), lsize, MPI_INTEGER4, MPI_STATUS_IGNORE, ierr)
      IF( ierr /= 0 ) CALL errore( 'iter_restart', 'error in MPI_FILE_READ iunsparseicb',1 )
      CALL MPI_FILE_SEEK(iunsparsejcb, lrepmatw4, MPI_SEEK_SET, ierr)
      IF( ierr /= 0 ) CALL errore( 'iter_restart', 'error in MPI_FILE_SEEK iunsparsejcb',1 )
      CALL MPI_FILE_READ(iunsparsejcb, sparsecb_j(:), lsize, MPI_INTEGER4, MPI_STATUS_IGNORE, ierr)
      IF( ierr /= 0 ) CALL errore( 'iter_restart', 'error in MPI_FILE_READ iunsparsejcb',1 )
      CALL MPI_FILE_SEEK(iunsparsetcb, lrepmatw4, MPI_SEEK_SET, ierr)
      IF( ierr /= 0 ) CALL errore( 'iter_restart', 'error in MPI_FILE_SEEK iunsparsetcb',1 )
      CALL MPI_FILE_READ(iunsparsetcb, sparsecb_t(:), lsize, MPI_INTEGER4, MPI_STATUS_IGNORE, ierr)
      IF( ierr /= 0 ) CALL errore( 'iter_restart', 'error in MPI_FILE_READ iunsparsetcb',1 )
      !
      CALL ibte(nind, etf_all, vkk_all, wkf_all, trans_probcb, efcb, &
                   sparsecb_q, sparsecb_k, sparsecb_i, sparsecb_j, sparsecb_t)
      ! 
      CALL MPI_FILE_CLOSE(iunepmatcb,ierr)
      IF( ierr /= 0 ) CALL errore( 'iter_restart', 'error in MPI_FILE_CLOSE',1)
      CALL MPI_FILE_CLOSE(iunsparseqcb,ierr)
      IF( ierr /= 0 ) CALL errore( 'iter_restart', 'error in MPI_FILE_CLOSE',1)
      CALL MPI_FILE_CLOSE(iunsparsekcb,ierr)
      IF( ierr /= 0 ) CALL errore( 'iter_restart', 'error in MPI_FILE_CLOSE',1)
      CALL MPI_FILE_CLOSE(iunsparseicb,ierr)
      IF( ierr /= 0 ) CALL errore( 'iter_restart', 'error in MPI_FILE_CLOSE',1)
      CALL MPI_FILE_CLOSE(iunsparsejcb,ierr)
      IF( ierr /= 0 ) CALL errore( 'iter_restart', 'error in MPI_FILE_CLOSE',1)
      CALL MPI_FILE_CLOSE(iunsparsetcb,ierr)
      IF( ierr /= 0 ) CALL errore( 'iter_restart', 'error in MPI_FILE_CLOSE',1)
      DEALLOCATE(trans_probcb)
      DEALLOCATE(sparsecb_q)
      DEALLOCATE(sparsecb_k)
      DEALLOCATE(sparsecb_i)
      DEALLOCATE(sparsecb_j)
      DEALLOCATE(sparsecb_t)
      ! 
    ENDIF
#endif  
    ! 
    !----------------------------------------------------------------------------
    END SUBROUTINE iter_restart
    !----------------------------------------------------------------------------
    ! 
  END MODULE transport_iter
