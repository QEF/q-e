SUBROUTINE exx_gs(nfi, c)
    !=======================================================================================
    ! Code Version 1.0 (Princeton University, September 2014)
    !=======================================================================================
    ! Note:  From this code exx_potential is returned after multiplying mixing parameter exxalfa.
    !        Later the exx_potential is added with GGA potential in forces.f90.
    !        In the future, full exx_potential should be returned and the mixing parameter exxalfa
    !        should be multiplied in forces.f90.
    !=======================================================================================
    !
    USE kinds,                   ONLY  : DP
    USE constants,               ONLY  : fpi
    USE fft_base,                ONLY  : dffts, dfftp
    USE mp,                      ONLY  : mp_barrier, mp_sum
    USE mp_global,               ONLY  : nproc_image, me_image, root_image, intra_image_comm, intra_bgrp_comm, me_bgrp
    USE parallel_include
    USE io_global,               ONLY  : stdout
    USE cell_base,               ONLY  : omega, ainv, h
    USE cell_base,               ONLY  : ibrav
    USE cell_base,               ONLY  : isotropic  !True if volume option is chosen for cell_dofree
    USE electrons_base,          ONLY  : nbsp, nbspx, nspin
    USE gvecw,                   ONLY  : ngw
    USE wannier_module,          ONLY  : wfc
    USE exx_module,              ONLY  : my_nbspx, my_nbsp, my_nxyz, index_my_nbsp,rk_of_obtl, lindex_of_obtl
    USE exx_module,              ONLY  : selfv, pairv, pair_dist
    USE exx_module,              ONLY  : pair_label, pair_status, pair_step
    USE exx_module,              ONLY  : exx_potential
    USE exx_module,              ONLY  : n_exx, sc_fac ! EXX step and the performance scaling factor
    USE exx_module,              ONLY  : coe_1st_derv, coeke, nord1, nord2, fornberg
    USE exx_module,              ONLY  : thdtood, odtothd_in_sp, thdtood_in_sp
    USE exx_module,              ONLY  : np_in_sp_s   , np_in_sp_me_s   , np_in_sp_p   , np_in_sp_me_p
    USE exx_module,              ONLY  : xx_in_sp,yy_in_sp,zz_in_sp,sc_xx_in_sp,sc_yy_in_sp,sc_zz_in_sp 
    USE energies,                ONLY  : exx
    USE printout_base,           ONLY  : printout_base_open, printout_base_unit, printout_base_close
    USE wannier_base,            ONLY  : neigh, dis_cutoff, texx_cube
    USE mp_wave,                 ONLY  : redistwfr
    !
    USE time_step,               ONLY  : tps                !md time in picoseconds
    USE io_global,               ONLY  : ionode             !logical for I/O node
    USE cp_main_variables,       ONLY  : iprint_stdout      !print control
    USE printout_base,           ONLY  : printout_base_close!close print unit
    USE printout_base,           ONLY  : printout_base_open !open print unit
    USE printout_base,           ONLY  : printout_base_unit !printout_base_unit
    USE exx_module,              ONLY  : dexx_dh
    USE exx_module,              ONLY  : exx_energy_cell_derivative
    USE exx_module,              ONLY  : exxalfa
    USE fft_helper_subroutines
    use exx_module, only : psime_pair_recv, psime_pair_send
    ! cubic domain related "use" variables
    USE exx_module,              ONLY  : selfrho, pairrho 
    USE exx_module,              ONLY  : nrg
    USE exx_module,              ONLY  : s_me_r, s_ps_r, p_me_r, p_ps_r
    USE exx_module,              ONLY  : n_s_me, n_s_ps, n_p_me, n_p_ps
    USE exx_module,              ONLY  : lmax
    USE exx_module,              ONLY  : me_cs
    USE exx_module,              ONLY  : me_rs
    USE exx_module,              ONLY  : me_ri
    USE exx_module,              ONLY  : me_rc
    USE exx_module,              ONLY  : coemicf !MCA/HK : dirty hack for std CG
    USE exx_module,              ONLY  : exx_energy_cell_derivative_cube
    !
    IMPLICIT NONE
    COMPLEX(DP)   c(ngw, nbspx)        ! wave functions at time t
#if defined(__MPI)
    !
    INTEGER  :: istatus(MPI_STATUS_SIZE)
#endif
    INTEGER     ir, ip, i, j,nfi, ierr, nnrtot, nr1s,nr2s,nr3s
    INTEGER     nj_max, iunit
    REAl(DP)    sa1,a(3),ha, hb, hc
    REAl(DP)    hcub, centerx, centery, centerz
    REAL(DP)    middle(3)
    REAL(DP)    d_pair            ! pair distance
    !
    REAL(DP),    ALLOCATABLE ::   vpsil(:,:)
    REAL(DP),    ALLOCATABLE ::   rhol(:),rho_in_sp(:),vl(:)
    ! cubic domain related variables
    INTEGER     l,m
    INTEGER     lid,gid
    INTEGER     s_me_r1,s_me_r2,s_me_r3,s_me_r4,s_me_r5,s_me_r6
    REAl(DP)    inv_omega
    REAl(DP)    dist, dxy, costheta, sintheta
    COMPLEX(DP) cxy
    REAl(DP)    plm(0:lmax,0:lmax)
    REAl(DP)    ha_proj(3), hb_proj(3), hc_proj(3)
    REAL(DP)    dq1,dq2,dq3
    REAL(DP)    dqs1,dqs2,dqs3
    REAL(DP)    start_timer, stop_timer
    REAL(DP),   ALLOCATABLE :: rhome(:),rhops(:),potme(:)
    INTEGER                 :: pos, oldest_step, guess_status  
    REAl(DP),   ALLOCATABLE :: psime(:)
    REAL(DP),    ALLOCATABLE ::   psi(:,:)
    !
    INTEGER   iobtl, gindex_of_iobtl, irank, rk_of_obtl_trcv, rk_of_obtl_tbs
    INTEGER   obtl_tbs, lindex_obtl_tbs, obtl_trcv, lindex_obtl_trcv 
    INTEGER   obtl_tbadd
    REAL(DP)  totalenergy, totalenergyg, tot_energy(nbsp)
    REAL(DP)  total_exx_derv(3,3), total_exx_derv_g(3,3)
    REAL(DP)  selfe, paire(neigh/2), &
        self_dexx_dhab(3,3), pair_dexx_dhab(3,3,neigh/2)
    !
    INTEGER,   ALLOCATABLE   :: isendreq(:)
    INTEGER,   ALLOCATABLE   :: irecvreq(:)
    INTEGER                  :: irecv_count
    INTEGER                  :: isend_count
    INTEGER                  :: itr
    INTEGER                  :: jtr
    INTEGER                  :: tran(3)
    INTEGER                  :: proc
    INTEGER                  :: tmp_iobtl
    INTEGER                  :: me
    !
    INTEGER                  :: k,jj,ii,ia,ib,ic,my_var,my_var2,my_var3,i_fac,va,cgstep
    INTEGER                  :: ndim,nogrp 
    INTEGER,    ALLOCATABLE  :: obtl_recv(:,:), obtl_send(:,:), num_recv(:), num_send(:)
    REAl(DP),   ALLOCATABLE  :: wannierc(:,:),wannierc_tmp(:,:)
    REAl(DP),   ALLOCATABLE  :: psil(:)
    REAL(DP),   ALLOCATABLE  :: exx_tmp(:,:),exx_tmp3(:,:)
    INTEGER,    ALLOCATABLE  :: sdispls(:), sendcount(:)
    INTEGER,    ALLOCATABLE  :: rdispls(:), recvcount(:)
    INTEGER,    ALLOCATABLE  :: sdispls1(:), sendcount1(:)
    INTEGER,    ALLOCATABLE  :: rdispls1(:), recvcount1(:)
    !
    REAL(DP) :: Jim(3,3)  ! jacobian [d/d x]        [d/d a]
    !                                |d/d y| = [J]  |d/d b|
    !                                [d/d z]        [d/d c]
    INTEGER            :: psgsn=3 !MCA: number of steps that arrays are stored. Memory overhead from here?
    !
    !=============================================================================================
    ! 
    CALL start_clock('exx_gs_setup')
    CALL exx_gs_setup_common
    if (texx_cube) then
      CALL exx_gs_setup_cube
    else
      CALL exx_gs_setup_sphere
    end if
    !
    !========================================================================
    !
    CALL stop_clock('exx_gs_setup')
    !
    !-------------------------------------------------------------------------
    ! Get the Wannier center and compute the pair overlap matrix 
    !-------------------------------------------------------------------------
    !
    CALL start_clock('exx_pairs')
    !
    ndim=MAX(nproc_image, nbsp)
    ALLOCATE (wannierc(3,ndim)); wannierc=0.0_DP 
    ALLOCATE (wannierc_tmp(3,nbsp)); wannierc_tmp=0.0_DP 
    !
    ! Adjust Cartesian coordinates of wannier centres according to periodic boundary conditions...
    ! N.B.: PBC are imposed here in the range [0,1)... 
    !
    DO iobtl=1,nbsp
      !
      wannierc_tmp(1,iobtl)=ainv(1,1)*wfc(1,iobtl)+ainv(1,2)*wfc(2,iobtl)+ainv(1,3)*wfc(3,iobtl)   ! s = h^-1 r
      wannierc_tmp(2,iobtl)=ainv(2,1)*wfc(1,iobtl)+ainv(2,2)*wfc(2,iobtl)+ainv(2,3)*wfc(3,iobtl)   ! s = h^-1 r
      wannierc_tmp(3,iobtl)=ainv(3,1)*wfc(1,iobtl)+ainv(3,2)*wfc(2,iobtl)+ainv(3,3)*wfc(3,iobtl)   ! s = h^-1 r
      !
      wannierc_tmp(1,iobtl)=wannierc_tmp(1,iobtl)-FLOOR(wannierc_tmp(1,iobtl))   ! impose PBC on s in range: [0,1)
      wannierc_tmp(2,iobtl)=wannierc_tmp(2,iobtl)-FLOOR(wannierc_tmp(2,iobtl))   ! impose PBC on s in range: [0,1)
      wannierc_tmp(3,iobtl)=wannierc_tmp(3,iobtl)-FLOOR(wannierc_tmp(3,iobtl))   ! impose PBC on s in range: [0,1)
      !
      wannierc(1,iobtl)=h(1,1)*wannierc_tmp(1,iobtl)+h(1,2)*wannierc_tmp(2,iobtl)+h(1,3)*wannierc_tmp(3,iobtl)   ! r = h s
      wannierc(2,iobtl)=h(2,1)*wannierc_tmp(1,iobtl)+h(2,2)*wannierc_tmp(2,iobtl)+h(2,3)*wannierc_tmp(3,iobtl)   ! r = h s
      wannierc(3,iobtl)=h(3,1)*wannierc_tmp(1,iobtl)+h(3,2)*wannierc_tmp(2,iobtl)+h(3,3)*wannierc_tmp(3,iobtl)   ! r = h s
      !
      wannierc_tmp(:,iobtl)=wannierc(:,iobtl) ! keep a temporary copy to compute pair indices
      !
    END DO
    ! 
    ! make copy of wannier centres when number of processors > number of bands
    ! 
    IF(nproc_image.GT.nbsp) THEN
      !
      DO iobtl=nbsp+1,nproc_image
        !
        ir=MOD(iobtl,nbsp)
        IF(ir.EQ.0 )THEN
          wannierc(:,iobtl)=wannierc(:,nbsp)
        ELSE
          wannierc(:,iobtl)=wannierc(:,ir)
        END IF
        !
      END DO
      !
    END IF
    !
    ! overlap is the unique neighbor list that for each band or processor image (ndim)
    ! num_recv is the number of unique neighbors for each band or processor image (ndim)
    !---------------------------------------------------------------------------------------------
    ALLOCATE (obtl_recv(neigh/2, ndim)); obtl_recv=0
    ALLOCATE (obtl_send(neigh, ndim));   obtl_send=0
    ALLOCATE (num_recv(ndim)); num_recv=0
    ALLOCATE (num_send(ndim)); num_send=0
    ! generate the unique neighbor list
    !
    CALL exx_index_pair(wannierc_tmp, obtl_recv, num_recv, nj_max, ndim)
    !
    IF (ALLOCATED(wannierc_tmp))            DEALLOCATE(wannierc_tmp)
    !---------------------------------------------------------------------------------------------
    DO itr = 1, ndim
      !
      DO jtr = 1, neigh/2
        !
        IF (obtl_recv(jtr, itr) .NE. 0) THEN
          !
          num_send( obtl_recv(jtr, itr) ) = num_send( obtl_recv(jtr, itr) ) + 1
          obtl_send( num_send( obtl_recv(jtr, itr) ), obtl_recv(jtr, itr) ) = itr
          !
        END IF
        !
      END DO
      !
    END DO
    !---------------------------------------------------------------------------------------------
    !
    CALL stop_clock('exx_pairs')
    !
    !-------------------------------------------------------------------------
    !
    ! Allocate variables to store potentials for 3 steps ....
    !
    IF (n_exx.EQ.0) THEN
      !
      ! the following variables are used in the extrapolation of exx potentials
      if (texx_cube) then
        ALLOCATE( selfv  ( n_s_ps, psgsn, my_nbspx), stat=ierr ); selfv=0.0_DP
        ALLOCATE( pairv  ( n_p_ps, psgsn, neigh, my_nbspx), stat=ierr ); pairv=0.0_DP
        ALLOCATE( selfrho( n_s_ps, psgsn, my_nbspx), stat=ierr ); selfrho=0.0_DP
        ALLOCATE( pairrho( n_p_ps, psgsn, neigh, my_nbspx), stat=ierr ); pairrho=0.0_DP
        ALLOCATE( pair_label( neigh, my_nbspx ), stat=ierr ); pair_label=0.0_DP
        ALLOCATE( pair_step( neigh, my_nbspx ), stat=ierr ); pair_step=0.0_DP
        ALLOCATE( pair_status( neigh, my_nbspx ), stat=ierr ); pair_status=0.0_DP
      else
        ALLOCATE( selfv  ( np_in_sp_s, psgsn, my_nbspx), stat=ierr ); selfv=0.0_DP
        ALLOCATE( pairv  ( np_in_sp_p, psgsn, neigh, my_nbspx), stat=ierr ); pairv=0.0_DP
        ALLOCATE( pair_dist( psgsn, neigh, my_nbspx), stat=ierr ); pair_dist=0.0_DP
      end if
      !
    END IF
    !
    !-------------------------------------------------------------------------
    !
    ! update exx step ...
    !
    n_exx = n_exx + 1
    !
    !=========================================================================
    !
    ! obtain orbitals on each local processor, stored in psi
    !
    ALLOCATE ( psi(  nnrtot, my_nbsp(me ) ) ); psi=0.0_DP
    ALLOCATE ( vpsil(nnrtot, my_nbsp(me ) ) ); vpsil=0.0_DP
    !
    CALL start_clock('r_orbital')
    !
    ! In c variable, the plane wave wavefunction, is distributed in nr3 parts
    ! and stored in each parallel mpitasks (or processors). 
    ! The exx_psi maps c in to psi, which is in real space grid.
    !
    ! c   -- one processor has a part of information from  all bands
    ! psi -- one processor has complete information of one band (or more)
    !
    CALL exx_psi(c, psi, nnrtot, my_nbsp, my_nxyz, nbsp) 
    !
    CALL stop_clock('r_orbital')
    !
    !===============================================================================
    !
    !           PAIR AND SELF POTENTIALS ARE CALCULATED IN ONE BIG LOOP 
    !
    !===============================================================================
    !
    !========================================================================
    !                      THE MOST OUTER LOOP STARTS:
    !========================================================================
    !
    ! flag identifying whether the pair communication is done or not
    !
    ALLOCATE( irecvreq( neigh * my_nbsp(me) ) )
    ALLOCATE( isendreq( neigh * my_nbsp(me) ) )
    !
    ! obtain psi in sphere (psil) for neighbors
    !
    call start_clock('exx_big_alloc')
    if (texx_cube) then
      if (.not.allocated(psime_pair_send)) then
        ALLOCATE ( psime_pair_send(n_p_me, neigh, my_nbsp(me)) ); psime_pair_send=0.0_DP
      end if
      if (.not.allocated(psime_pair_recv)) then
        ALLOCATE ( psime_pair_recv(n_p_me, neigh, my_nbsp(me)) ); psime_pair_recv=0.0_DP
      end if
    else
      if (.not.allocated(psime_pair_send)) then
        ALLOCATE ( psime_pair_send(np_in_sp_me_p, neigh, my_nbsp(me)) ); psime_pair_send=0.0_DP
      end if
      if (.not.allocated(psime_pair_recv)) then
        ALLOCATE ( psime_pair_recv(np_in_sp_me_p, neigh, my_nbsp(me)) ); psime_pair_recv=0.0_DP
      end if
    end if
    call stop_clock('exx_big_alloc')
    !
    ! initialize totalenergy and derivatives
    !
    totalenergy = 0.0_DP
    total_exx_derv(:,:) =0.0_DP
    !
    ! my_var is the maximum of nbsp or nproc_image
    my_var = MAX(nproc_image, nbsp)
    !
#if defined(__MPI)
    irecv_count = 0
    isend_count = 0
    ! ---------- start of Non-blocking communication (1) -----------
    ! In this region of the code, we start with non-blocking send of
    ! local orbitals to form MLWF product potentials on the non-self
    ! multipole expansion domain. This part allows all pairs to be 
    ! send simultaneously instead of one-by-one...
    ! --------------------------------------------------------------
    CALL start_clock('send_psi')
    ! 
    ! we should use my_nbspx (maxval(my_nbsp(me))) here
    !
    DO iobtl = 1, my_nbspx
      ! 
      gindex_of_iobtl = index_my_nbsp(iobtl, me)
      if (gindex_of_iobtl.gt.nbsp) cycle ! index_my_nbsp is set to nbsp+1 when not assigned
      !
      !========================================================================
      !
      !prep receives
      DO itr = 1, neigh/2
        !
        obtl_tbs = obtl_recv( itr, gindex_of_iobtl )
        !
        IF ( obtl_tbs .NE. 0 ) THEN
          !
          rk_of_obtl_tbs  = rk_of_obtl(obtl_tbs)
          !
          irecv_count = irecv_count + 1
          if (texx_cube) then
            CALL MPI_IRECV( psime_pair_recv(1, itr, iobtl), n_p_me, MPI_DOUBLE_PRECISION, rk_of_obtl_tbs, &
              obtl_tbs*ndim+gindex_of_iobtl, intra_image_comm, irecvreq(irecv_count), ierr)
          else
            CALL MPI_IRECV( psime_pair_recv(1, itr, iobtl), np_in_sp_me_p, MPI_DOUBLE_PRECISION, rk_of_obtl_tbs, &
              obtl_tbs*ndim+gindex_of_iobtl, intra_image_comm, irecvreq(irecv_count), ierr)
          end if
          !
        END IF
        !
      END DO
      !
      !prep sends
      DO itr = 1, neigh
        !
        obtl_trcv = obtl_send( itr, gindex_of_iobtl )
        !
        IF ( obtl_trcv .NE. 0 ) THEN
          !
          rk_of_obtl_trcv = rk_of_obtl(obtl_trcv)
          !
          ! calculate mid point of two wannier centers
          CALL getmiddlewc( wannierc(1, gindex_of_iobtl), wannierc(1, obtl_trcv), h, ainv, middle )
          !
          ! calculate translation vector from the center of the box
          CALL getsftv( nr1s, nr2s, nr3s, h, ainv, middle, tran)
          !
          ! get the localized psi around the mid point of two wannier centers
          ! note: the psime is centered at the center of the box
          ! (using the translation vector "tran" from middle of wfc to the center of box)
          if (texx_cube) then
            CALL getpsicb( nrg, p_me_r, psi(1,iobtl), psime_pair_send(1, itr, iobtl), tran)
          else
            CALL getpsil( nnrtot, np_in_sp_me_p, psi(1, iobtl), psime_pair_send(1, itr, iobtl), tran)
          end if
          !
          isend_count = isend_count + 1
          if (texx_cube) then
            CALL MPI_ISEND( psime_pair_send(1, itr, iobtl), n_p_me, MPI_DOUBLE_PRECISION, rk_of_obtl_trcv, &
                            gindex_of_iobtl*ndim+obtl_trcv, intra_image_comm, isendreq(isend_count), ierr)
          else
            CALL MPI_ISEND( psime_pair_send(1, itr, iobtl), np_in_sp_me_p, MPI_DOUBLE_PRECISION, rk_of_obtl_trcv, &
              gindex_of_iobtl*ndim+obtl_trcv, intra_image_comm, isendreq(isend_count), ierr)
          end if
          !
        END IF
        !
      END DO
      !========================================================================
      !
    END DO
    !
    !==========================================================================
    !
    CALL stop_clock('send_psi')
    !
    !==========================================================================
    !
    CALL start_clock('send_psi_wait')
    !
    DO itr = 1, isend_count
      CALL MPI_WAIT(isendreq(itr), istatus, ierr)
    END DO
    !
    DO itr = 1, irecv_count
      CALL MPI_WAIT(irecvreq(itr), istatus, ierr)
    END DO
    !
    CALL stop_clock('send_psi_wait')
    ! ------- end of Non-blocking communication (1): synchronize by wait ------
#endif
    !
    !=========================================================================
    ! after this loop ( do irank ), all the processor got all the overlapping orbitals 
    ! for the i_obtl orbital and ready to calculate pair potential 
    !=========================================================================
    !
    !=========================================================================
    !                              CALCULATION STARTS HERE
    !=========================================================================
    ! printout header for the cgsteps
    !
    IF ((MOD(nfi,iprint_stdout).EQ.0)) THEN
      !   
      IF (ionode) THEN
        !   
        iunit=printout_base_unit("ncg")
        !
        CALL printout_base_open("ncg")
        !   
        WRITE(iunit,'(I8,F16.8)')nfi,tps
        !   
      END IF    
      !   
    END IF    
    !
    CALL start_clock('getpairv')
    !
    if (texx_cube) then
      ! Do some allocations
      IF(.not.ALLOCATED(psime )) ALLOCATE( psime(max(n_p_me,n_s_me)) ); 
      IF(.not.ALLOCATED(rhome )) ALLOCATE( rhome(max(n_p_me,n_s_me)) );
      IF(.not.ALLOCATED(rhops )) ALLOCATE( rhops(max(n_p_ps,n_s_ps)) ); 
      IF(.not.ALLOCATED(potme )) ALLOCATE( potme(max(n_p_me,n_s_me)) ); 
    end if
    DO iobtl = 1, my_nbspx
      ! 
      middle(:)=0.0_DP
      !
      gindex_of_iobtl = index_my_nbsp(iobtl, me)
      if (gindex_of_iobtl.gt.nbsp) cycle ! index_my_nbsp is set to nbsp+1 when not assigned
      !
      IF ( gindex_of_iobtl .LE. my_var) THEN
        !
        !-- second loop starts: calculate overlapping potential with the j_th orbitals --
        !
        ! my_var3 is the unique neighbor number for gindex_of_iobtl
        my_var3=num_recv( gindex_of_iobtl )
        !
        do j = 1, my_var3
          !
          ! my_var2 is the global index of the unique pair to gindex_of_iobtl
          my_var2=obtl_recv(j,gindex_of_iobtl)
          !
          IF (my_var2 .NE. 0) THEN
            !==================================================================================
            !                      find the position of previous v/rho
            !==================================================================================
            ! 
            if (texx_cube) then
              call solve_a_nonself_pair_cube
            else
              call solve_a_nonself_pair_sphere
            end if
            !
            paire(j) = paire(j) * 0.5_DP* hcub             ! volume element hcub and trapezoidal rule prefactor 0.5_DP are included
            totalenergy = totalenergy + 2.0_DP*paire(j)    ! the factor of two comes from the identity of ij and ji pair
            !
            IF (.NOT. (isotropic .AND. (ibrav.EQ.1) )) THEN
              CALL start_clock('exx_cell_derv')
              !
              ! EXX cell derivative (note: exxalfa is included in vofrho.f90 when calculate stress)
              !
              if (texx_cube) then
                CALL exx_energy_cell_derivative_cube(p_me_r, p_ps_r, tran, rhome, potme, &
                  ha_proj, hb_proj, hc_proj, Jim, pair_dexx_dhab(:,:,j))
              else
                CALL exx_energy_cell_derivative(np_in_sp_me_p, np_in_sp_p, tran,&
                  vl, ha, hb, hc, rhol, pair_dexx_dhab(:,:,j))
              end if
              !
              ! volume element hcub and trapezoidal rule prefactor 0.5_DP are included
              !
              pair_dexx_dhab(:,:,j) = pair_dexx_dhab(:,:,j)*0.5_DP*hcub
              !
              ! accumulate the derivative from different pair terms
              !
              total_exx_derv(:,:) = total_exx_derv(:,:) + 2.0_DP*pair_dexx_dhab(:,:,j)
              !
              !
              ! if isotropic => calculate the stress tensor in vofrho.f90
              !
              CALL stop_clock('exx_cell_derv')
            END IF
            !
            if (.not.texx_cube) then
              IF (ALLOCATED(psil))            DEALLOCATE(psil)
              IF (ALLOCATED(rhol))            DEALLOCATE(rhol)
              IF (ALLOCATED(rho_in_sp))       DEALLOCATE(rho_in_sp)
              IF (ALLOCATED(vl))              DEALLOCATE(vl)
            end if
            !
          END IF
          !
        END DO !for j
        !
      END IF !gindex_of_iobtl <= nbsp
      !
    END DO
    CALL stop_clock('getpairv')
    !
    !===============================================================================
    ! After this loop, each processor finished the pair potentials for the 
    ! iobtl orbital, and shall talk to send/recv vpsiforj
    !===============================================================================
    !
    !===============================================================================
    !                INITIALIZE SEND VPSI BEFORE CALCULATE PAIR
    !===============================================================================
    !
#if defined(__MPI)
    irecv_count = 0
    isend_count = 0
    !
    ! ---------- start of Non-blocking communication (2) -----------
    ! In this region of the code, we use non-blocking MPI feature to
    ! send the EXX contributions to orbital force and asynchronously
    ! overlap with the self-exchange computation...
    ! --------------------------------------------------------------
    CALL start_clock('send_v')
    !
    !========================================================================
    !
    DO iobtl = 1, my_nbspx
      !
      !========================================================================
      !
      gindex_of_iobtl = index_my_nbsp(iobtl, me)
      if (gindex_of_iobtl.gt.nbsp) cycle ! index_my_nbsp is set to nbsp+1 when not assigned
      !
      DO itr = 1, neigh/2
        !
        obtl_trcv = obtl_recv( itr, gindex_of_iobtl )
        !
        IF ( obtl_trcv .NE. 0 ) THEN
          !
          rk_of_obtl_trcv  = rk_of_obtl(obtl_trcv)
          !
          isend_count = isend_count + 1
          if (texx_cube) then
            CALL MPI_ISEND( psime_pair_recv(1,itr,iobtl), n_p_me, MPI_DOUBLE_PRECISION, rk_of_obtl_trcv, &
              gindex_of_iobtl*ndim+obtl_trcv, intra_image_comm, isendreq(isend_count), ierr)
          else
            CALL MPI_ISEND( psime_pair_recv(1,itr,iobtl), np_in_sp_me_p, MPI_DOUBLE_PRECISION, rk_of_obtl_trcv, &
              gindex_of_iobtl*ndim+obtl_trcv, intra_image_comm, isendreq(isend_count), ierr)
          end if
          !
          ! WRITE(my_unit,*) "itr = ", itr
          ! WRITE(my_unit,*) "iobtl = ", iobtl
          ! WRITE(my_unit,*) "gindex_of_iobtl = ",gindex_of_iobtl
          ! WRITE(my_unit,*) "obtl_trcv: ", obtl_trcv
          ! WRITE(my_unit,*) "rk_of_obtl_trcv: ", rk_of_obtl_trcv
          ! WRITE(my_unit,*) "tag: ", gindex_of_iobtl*nbsp+obtl_trcv
          !
        END IF
        !
      END DO
      !
      !========================================================================
      !
      !MCA: send and receive now change their roles. Before, we were communicating
      !the orbitals. Now, it is the potential that we got out of Poisson. Send is
      !recieving and recv is sending. For now on, only psime_pair_send is used for 
      !computation.
      !
      DO itr = 1, neigh
        !
        obtl_tbs = obtl_send( itr, gindex_of_iobtl )
        !
        IF ( obtl_tbs .NE. 0 ) THEN
          !
          rk_of_obtl_tbs = rk_of_obtl(obtl_tbs)
          !
          irecv_count = irecv_count + 1
          if (texx_cube) then
            CALL MPI_IRECV( psime_pair_send(1, itr, iobtl), n_p_me, MPI_DOUBLE_PRECISION, rk_of_obtl_tbs, &
              obtl_tbs*ndim+gindex_of_iobtl, intra_image_comm, irecvreq(irecv_count), ierr)
          else
            CALL MPI_IRECV( psime_pair_send(1, itr, iobtl), np_in_sp_me_p, MPI_DOUBLE_PRECISION, rk_of_obtl_tbs, &
              obtl_tbs*ndim+gindex_of_iobtl, intra_image_comm, irecvreq(irecv_count), ierr)
          end if
          !
          ! WRITE(my_unit,*) "itr = ", itr
          ! WRITE(my_unit,*) "iobtl = ", iobtl
          ! WRITE(my_unit,*) "gindex_of_iobtl = ", gindex_of_iobtl
          ! WRITE(my_unit,*) "obtl_tbs: ", obtl_tbs
          ! WRITE(my_unit,*) "rk_of_obtl_tbs: ", rk_of_obtl_tbs
          ! WRITE(my_unit,*) "tag: ", obtl_tbs*nbsp+gindex_of_iobtl
          !
          !
        END IF
        !
      END DO
      !
      !========================================================================
      !
    END DO
    !
    !========================================================================
    !
    CALL stop_clock('send_v')
#endif
    !
    !=========================================================================
    !                                
    !              SELF POTENTIAL FOR EACH ORBITAL STARTS HERE
    !                                
    !=========================================================================
    !
    CALL start_clock('getselfv')
    !
    !=========================================================================
    DO iobtl = 1, my_nbspx
      !
      IF (iobtl.LE.my_nbsp(me)) THEN ! skip when the loop of my_nbspx goes outside of scope
        !
        IF (me.GT.(nbsp*(sc_fac-1))) THEN ! compatible with more processors than nbsp
          !
          gindex_of_iobtl =  index_my_nbsp(iobtl, me)
          if (gindex_of_iobtl.gt.nbsp) cycle ! index_my_nbsp is set to nbsp+1 when not assigned
          !
          if (texx_cube) then
            call solve_a_self_pair_cube
          else
            call solve_a_self_pair_sphere
          end if
          selfe = selfe * 0.5_DP * hcub               ! volume element hcub and trapezoidal rule prefactor 0.5_DP are included
          totalenergy = totalenergy + selfe
          !
          !IF(me .GT. nbsp) THEN
          !  vpsil(:,iobtl) = 0.0_DP
          !END IF
          !
          IF (.NOT. (isotropic .AND. (ibrav.EQ.1))) THEN
            !
            !  EXX cell derivative (note: need to include exxalfa later)
            CALL start_clock('exx_cell_derv')
            !
            if (texx_cube) then
              CALL exx_energy_cell_derivative_cube(s_me_r, s_ps_r, tran, rhome, potme, &
                ha_proj, hb_proj, hc_proj, Jim, self_dexx_dhab(:,:))
            else
              CALL exx_energy_cell_derivative(np_in_sp_me_s, np_in_sp_s, tran,&
                vl, ha, hb, hc, rhol, self_dexx_dhab(:,:))
            end if
            !
            ! volume element hcub and trapezoidal rule prefactor 0.5_DP are included
            !
            self_dexx_dhab(:,:) = self_dexx_dhab(:,:)*0.5_DP*hcub
            !
            ! combine derivative with pair terms
            !
            total_exx_derv(:,:) = total_exx_derv(:,:) + self_dexx_dhab(:,:)
            !
            CALL stop_clock('exx_cell_derv')
            !
            ! if isotropic => calculate the stress tensor in vofrho.f90
            !
          END IF
          !
          if (.not.texx_cube) then
            IF (ALLOCATED(psil))            DEALLOCATE(psil)
            IF (ALLOCATED(rhol))            DEALLOCATE(rhol)
            IF (ALLOCATED(rho_in_sp))       DEALLOCATE(rho_in_sp)
            IF (ALLOCATED(vl))              DEALLOCATE(vl)
          end if
          !
        END IF ! me
        !
      END IF !iobtl 
      !
    END DO ! iobtl
    if (texx_cube) then
      IF (ALLOCATED(psime))           DEALLOCATE(psime)
      IF (ALLOCATED(rhome))           DEALLOCATE(rhome)
      IF (ALLOCATED(rhops))           DEALLOCATE(rhops)
      IF (ALLOCATED(potme))           DEALLOCATE(potme)
    end if
    CALL stop_clock('getselfv')
    !========================================================================
    !
#if defined(__MPI)
    CALL start_clock('send_v_wait')
    !
    DO itr = 1, isend_count
      CALL MPI_WAIT(isendreq(itr), istatus, ierr)
    END DO
    !
    DO itr = 1, irecv_count
      CALL MPI_WAIT(irecvreq(itr), istatus, ierr)
    END DO
    !
    CALL stop_clock('send_v_wait')
    ! ------- end of Non-blocking communication (2): synchronize by wait ------
    !
    !========================================================================
    CALL start_clock('force_rec')
    !
    DO iobtl = 1, my_nbspx
      !
      gindex_of_iobtl =  index_my_nbsp(iobtl, me)
      if (gindex_of_iobtl.gt.nbsp) cycle ! index_my_nbsp is set to nbsp+1 when not assigned
      !
      DO itr = 1, neigh
        !
        obtl_tbadd = obtl_send(itr, gindex_of_iobtl)
        !
        IF ( obtl_tbadd .NE. 0 ) THEN
          !
          CALL getmiddlewc( wannierc(1,gindex_of_iobtl), wannierc(1,obtl_tbadd), h, ainv, middle )
          !
          ! calculate translation vector from the center of the box
          CALL getsftv( nr1s, nr2s, nr3s, h, ainv, middle, tran )
          !
          ! upadate vpsil PBE0 
          !
          if (texx_cube) then
            CALL updateforce_rec(nrg, p_me_r, vpsil(1,iobtl), psime_pair_send(1,itr,iobtl), tran)
          else
            !$omp parallel do private(ir) 
            DO ip = 1, np_in_sp_me_p
              CALL l2goff (ip,ir,tran) ! local is centered at box center; global index is offset by tran
              vpsil(ir,iobtl) = vpsil(ir,iobtl) + psime_pair_send(ip,itr,iobtl)
            END DO
            !$omp end parallel do
          end if
          !
        END IF
        !
      END DO
      !
    END DO ! iobtl
    !
    CALL stop_clock('force_rec')
    !========================================================================
#endif
    !
    CALL start_clock('totalenergy')
    !
    totalenergyg=0.0_DP ! mpi reduction variable initialization
    exx=0.0_DP          ! exx energy (used to handle the open/closed shell energy)
    ! 
#if defined(__MPI)
    ! collect the totalenergy of each mpi task to totalenergyg
    CALL MPI_ALLREDUCE(totalenergy, totalenergyg, 1, MPI_DOUBLE_PRECISION, &
        &                        MPI_SUM, intra_image_comm, ierr)
    !
#else
    totalenergyg = totalenergy
#endif
    exx = totalenergyg
    IF (nspin .EQ. 1) exx = exx + totalenergyg ! if closed shell double the totalenergy
    !
    !WRITE(stdout, '("EXX Energy",2F30.14," step",I7)')exx,totalenergyg*2.0_DP, nfi
    !
    CALL stop_clock('totalenergy')
    !
    CALL start_clock('exx_cell_derv')
    !
    total_exx_derv_g(:,:) = 0.0_DP ! mpi reduction variable initialization
    !
    IF (.NOT. (isotropic .AND. (ibrav.EQ.1))) THEN
#if defined(__MPI)
      ! collect the total_exx_derv of each mpi task to total_exx_derv_g
      CALL MPI_ALLREDUCE(total_exx_derv(:,:), total_exx_derv_g(:,:), 9, &
          MPI_DOUBLE_PRECISION, MPI_SUM, intra_image_comm, ierr)
      !
#else
      total_exx_derv_g(:,:) = total_exx_derv(:,:)
#endif
    END IF
    !
    ! for closed shell case inclued spin factor of 2
    dexx_dh(:,:) = total_exx_derv_g(:,:)
    IF (nspin .EQ. 1) dexx_dh(:,:) = dexx_dh(:,:) + total_exx_derv_g(:,:)
    !
    CALL stop_clock('exx_cell_derv')
    !
    ! Local to global distribution of EXX potential
    ! vpsil (local) --> exx_potential (global)
    !
    CALL start_clock('vl2vg')
    exx_potential=0.0_DP
    !
    IF (nproc_image .LE. nbsp) THEN 
      !
#ifdef __MPI
      CALL redistwfr ( exx_potential, vpsil, my_nxyz, my_nbsp, intra_image_comm, -1 )
#else
      exx_potential = vpsil
#endif
      !
    ELSE
      !
      !-----------Zhaofeng's vpsil (local) to exx_potential (global) -----------
      !
      nogrp = fftx_ntgrp(dffts)
      !
      ALLOCATE( sdispls(nproc_image), sendcount(nproc_image) ); sdispls=0; sendcount=0
      ALLOCATE( rdispls(nproc_image), recvcount(nproc_image) ); rdispls=0; recvcount=0 
      ALLOCATE( sdispls1(nogrp), sendcount1(nogrp) ); sdispls1=0; sendcount1=0
      ALLOCATE( rdispls1(nogrp), recvcount1(nogrp) ); rdispls1=0; recvcount1=0
      !
      DO proc = 1, nproc_image
        !
        IF (me <= nogrp*nr3s) THEN
          sendcount(proc) =nr1s*nr2s/nogrp
        ELSE
          sendcount(proc) = 0
        END IF
        !proc 1 holds  the nr1s*nr2s/nogrp information of 1024 orbital(duplicate)
        !proc 640 as well
        !however, 641 to 1024 idle
        IF (proc <= nogrp*nr3s) THEN
          recvcount(proc)=nr1s*nr2s/nogrp
        ELSE
          recvcount(proc)=0
        END IF
        !
      END DO
      !
      sdispls(1) = 0
      rdispls(1) = 0
      !
      DO proc = 2,  nproc_image
        sdispls(proc)=  sdispls(proc-1) + sendcount(proc-1)
        rdispls(proc) = rdispls(proc-1) + recvcount(proc-1)
      END DO
      !
      ALLOCATE(exx_tmp (dffts%nnr,nproc_image/nogrp)); exx_tmp=0.0_DP
      ALLOCATE(exx_tmp3(dffts%nnr,nproc_image/nogrp)); exx_tmp3=0.0_DP
      !
#if defined(__MPI)
      !
      CALL mp_barrier( intra_image_comm )
      CALL MPI_ALLTOALLV(vpsil(1,1), recvcount,rdispls,MPI_DOUBLE_PRECISION, &
          &           exx_tmp, sendcount,sdispls, MPI_DOUBLE_PRECISION, &
          &           intra_image_comm, ierr)
#endif
      !
      va = dffts%nnr/nogrp
      DO j=1,nproc_image/nogrp,2
        DO i=1,2*nogrp
          ii=((i-1)/2)*va 
          jj=j+mod(i-1,2)
          ia=(i-1-((i-1)/nogrp)*nogrp)*va
          ib=j+(i-1)/nogrp
          !$omp parallel do 
          DO ir=1,va
            exx_tmp3(ii+ir,jj)=exx_tmp(ia+ir,ib)
          END DO
          !$omp end parallel do 
        END DO
      END DO
      !
      DO proc = 1 , nogrp
        sendcount1(proc) = dffts%nnr/nogrp
        recvcount1(proc) = dffts%nnr/nogrp
      END DO
      !
      rdispls1(1) = 0
      sdispls1(1) = 0
      !
      DO proc = 2, nogrp
        sdispls1(proc) = sdispls1(proc-1) + sendcount1(proc-1)
        rdispls1(proc) = rdispls1(proc-1) + recvcount1(proc-1)
      END DO
      !
#if defined(__MPI)
      !
      DO ir=1,nproc_image/nogrp
        CALL mp_barrier( fftx_tgcomm(dffts) )
        CALL MPI_ALLTOALLV(exx_tmp3(1,ir), sendcount1, sdispls1, MPI_DOUBLE_PRECISION, &
            &         exx_potential(1,ir),recvcount1, rdispls1, MPI_DOUBLE_PRECISION, &
            &         fftx_tgcomm(dffts), ierr)
      END DO
#endif
      !
      DO ir=1,nbsp/nogrp
        DO i=1,sc_fac-1
          ii=i*nbsp/nogrp
          !$omp parallel do 
          DO ia=1,dffts%nnr
            exx_potential(ia,ir)=exx_potential(ia,ir)+exx_potential(ia,ir+ii)
          END DO
          !$omp end parallel do 
        END DO
      END DO
      !
      !-----------Zhaofeng's vpsil (local) to exx_potential (global) -----------
      !
    END IF ! vl2vg
    !
    CALL stop_clock('vl2vg')
    !
    !==============================================================================
    IF (ALLOCATED(vpsil))           DEALLOCATE(vpsil)
    IF (ALLOCATED(psi))             DEALLOCATE(psi)
    IF (ALLOCATED(isendreq))        DEALLOCATE(isendreq)
    IF (ALLOCATED(irecvreq))        DEALLOCATE(irecvreq)
    IF (ALLOCATED(wannierc))        DEALLOCATE(wannierc)
    IF (ALLOCATED(obtl_recv))       DEALLOCATE(obtl_recv)
    IF (ALLOCATED(obtl_send))       DEALLOCATE(obtl_send)
    IF (ALLOCATED(num_recv))        DEALLOCATE(num_recv)
    IF (ALLOCATED(num_send))        DEALLOCATE(num_send)
    IF (ALLOCATED(exx_tmp))         DEALLOCATE(exx_tmp)
    IF (ALLOCATED(exx_tmp3))        DEALLOCATE(exx_tmp3)
    IF (ALLOCATED(sdispls))         DEALLOCATE(sdispls)
    IF (ALLOCATED(rdispls))         DEALLOCATE(rdispls)
    IF (ALLOCATED(sdispls1))        DEALLOCATE(sdispls1)
    IF (ALLOCATED(rdispls1))        DEALLOCATE(rdispls1)
    IF (ALLOCATED(sendcount))       DEALLOCATE(sendcount)
    IF (ALLOCATED(recvcount))       DEALLOCATE(recvcount)
    IF (ALLOCATED(sendcount1))      DEALLOCATE(sendcount1)
    IF (ALLOCATED(recvcount1))      DEALLOCATE(recvcount1)
    !
    RETURN
  contains

    SUBROUTINE  exx_gs_setup_common()
      IMPLICIT NONE
      !
      ! make processor index start from 1
      !
      me = me_image+1
      !
      ! number of real space gird along each lattice parameter directions (a1, a2, a3)
      !
      nr1s=dfftp%nr1; nr2s=dfftp%nr2; nr3s=dfftp%nr3 
      !
      ! the length of each lattice parameters
      !
      a(1)=DSQRT(h(1,1)*h(1,1)+h(2,1)*h(2,1)+h(3,1)*h(3,1))   ! lattice 1 
      a(2)=DSQRT(h(1,2)*h(1,2)+h(2,2)*h(2,2)+h(3,2)*h(3,2))   ! lattice 2 
      a(3)=DSQRT(h(1,3)*h(1,3)+h(2,3)*h(2,3)+h(3,3)*h(3,3))   ! lattice 3 
      !
      ! grid spacing in each lattice parameters
      !
      ha = a(1) / DBLE(nr1s)  !grid spacing in Lattice 1 direction
      hb = a(2) / DBLE(nr2s)  !grid spacing in Lattice 2 direction
      hc = a(3) / DBLE(nr3s)  !grid spacing in Lattice 3 direction
      !
      ! total number of real space grid points in the global mesh 
      ! and the corresponding volume elements for each grid point
      !
      nnrtot = nr1s * nr2s * nr3s
      hcub = omega / DBLE(nnrtot) !nnrtot in parallel
      !
      ! the x,y,z coordinates of the center of the box (gird center)
      ! NOTE: center of the box is set to grid point at int(nr1/2), int(nr2/2), int(nr3/2) for every cell
      !
      centerx = h(1,1)*DBLE(INT(nr1s/2))+h(1,2)*DBLE(INT(nr2s/2))+h(1,3)*DBLE(INT(nr3s/2))   ! r = h s
      centery = h(2,1)*DBLE(INT(nr1s/2))+h(2,2)*DBLE(INT(nr2s/2))+h(2,3)*DBLE(INT(nr3s/2))   ! r = h s
      centerz = h(3,1)*DBLE(INT(nr1s/2))+h(3,2)*DBLE(INT(nr2s/2))+h(3,3)*DBLE(INT(nr3s/2))   ! r = h s
      !
      ! inverse volume
      !
      sa1 = 1.0_DP/omega
      !
      ! Compute coeke and renormalize
      ! This part needs to be done once in constant volume simulation and
      ! needs to be done every step in variable cell simulationulations ...
      !
      ! get ha*d/da, ha^2*d^2/da^2 stencil and cross coefficients
      !   1. for the finite difference coefficients, we follow B. Fornberg in 
      !       Math. Comp. 51 (1988), 699-706
      !
      CALL fornberg(nord1, nord2,coe_1st_derv(:,1),coeke(:,1,1),coeke(:,1,2),ierr)
      !
      IF (ierr .ne. 0) THEN
        WRITE(stdout,*) ' ERROR: Wrong parameter in CALL of Fornberg'
        WRITE(stdout,*) ' STOP in exx_gs'
        RETURN
      END IF
      !RENORMALIZE COEKES WITH RESPECT TO THE GRID SPACING
      ! first derivative coefficients
      !
      coe_1st_derv(:,3) = coe_1st_derv(:,1)/hc ! d/dc stencil
      coe_1st_derv(:,2) = coe_1st_derv(:,1)/hb ! d/db stencil
      coe_1st_derv(:,1) = coe_1st_derv(:,1)/ha ! d/da stencil
      !
      ! NOTE: in the second derivatives there is a additional factor of
      !       -4*pi because we merege that from the Poisson equation
      !
      !                \nabla^2 V = -4*\pi \rho
      !
      ! axial derivatives
      !
      coeke(:,3,3) = -coeke(:,1,1)/(hc*hc*fpi) ! -d^2/dc^2/4pi stencil
      coeke(:,2,2) = -coeke(:,1,1)/(hb*hb*fpi) ! -d^2/db^2/4pi stencil
      coeke(:,1,1) = -coeke(:,1,1)/(ha*ha*fpi) ! -d^2/da^2/4pi stencil
      !
      ! cross derivatives
      !
      coeke(:,2,3) = -coeke(:,1,2)/(hb*hc*fpi) ! -d^2/dbdc/4pi stencil
      coeke(:,1,3) = -coeke(:,1,2)/(ha*hc*fpi) ! -d^2/dadc/4pi stencil
      coeke(:,1,2) = -coeke(:,1,2)/(ha*hb*fpi) ! -d^2/dadb/4pi stencil
      !
      ! -- Jacobian for the general (non-orthogonal) --
      ! please see the following reference for details <todo: EXX paper>
      !
      ! J = transpose(ainv).(diag(a))
      !
      Jim(:,1) = ainv(1,:)*a(1)
      Jim(:,2) = ainv(2,:)*a(2) ! i={xyz}, m={abc}
      Jim(:,3) = ainv(3,:)*a(3)
      !
      ! -- weigh coeke with the Jacobian --
      !
      ! axial derivatives
      !
      coeke(:,3,3) = (Jim(1,3)**2+Jim(2,3)**2+Jim(3,3)**2)*coeke(:,3,3)
      coeke(:,2,2) = (Jim(1,2)**2+Jim(2,2)**2+Jim(3,2)**2)*coeke(:,2,2)
      coeke(:,1,1) = (Jim(1,1)**2+Jim(2,1)**2+Jim(3,1)**2)*coeke(:,1,1)
      !
      ! cross derivatives (needed for non-othogonal grids in the second derivatives)
      !
      coeke(:,2,3) = 2.0_DP*(Jim(1,2)*Jim(1,3)+Jim(2,2)*Jim(2,3)+Jim(3,2)*Jim(3,3))*coeke(:,2,3)
      coeke(:,1,3) = 2.0_DP*(Jim(1,1)*Jim(1,3)+Jim(2,1)*Jim(2,3)+Jim(3,1)*Jim(3,3))*coeke(:,1,3)
      coeke(:,1,2) = 2.0_DP*(Jim(1,1)*Jim(1,2)+Jim(2,1)*Jim(2,2)+Jim(3,1)*Jim(3,2))*coeke(:,1,2)
      coeke(:,3,2) = coeke(:,2,3) ! symmetry of coeke
      coeke(:,2,1) = coeke(:,1,2) ! symmetry of coeke
      coeke(:,3,1) = coeke(:,1,3) ! symmetry of coeke
      !
      ! a samall check on the shape of user defined cell (if any)
      !
      IF ((ibrav.EQ.0).AND.(nfi.EQ.1)) THEN
        WRITE(stdout,*) 'EXX info: If you are using an orthogonal cell without its cell vectors&
          & aligned to the xyz directions, the EXX calculation may be twice more expensive.'
      END IF
      RETURN
    END SUBROUTINE exx_gs_setup_common

    subroutine  exx_gs_setup_cube()
      implicit none
      psgsn=1
      ! consider a grid an unit cell, the lattice vectors for the grid
      ha_proj(:) = ha*h(:,1)/a(1)
      hb_proj(:) = hb*h(:,2)/a(2)
      hc_proj(:) = hc*h(:,3)/a(3)
      inv_omega = 1.0_DP/omega
      s_me_r1=s_me_r(1)
      s_me_r2=s_me_r(2)
      s_me_r3=s_me_r(3)
      s_me_r4=s_me_r(4)
      s_me_r5=s_me_r(5)
      s_me_r6=s_me_r(6)
      DO k = s_me_r(3),s_me_r(6)
        DO j = s_me_r(2),s_me_r(5)
          DO i = s_me_r(1),s_me_r(4)
            !---------------------------------------------------------------------------------------
            dqs1 = (DBLE(i)/DBLE(nr1s)) - DBLE(INT(nr1s/2))/DBLE(nr1s)
            dqs2 = (DBLE(j)/DBLE(nr2s)) - DBLE(INT(nr2s/2))/DBLE(nr2s)
            dqs3 = (DBLE(k)/DBLE(nr3s)) - DBLE(INT(nr3s/2))/DBLE(nr3s)
            !
            ! Here we are computing distances between Grid points and center of the simulation cell, so no MIC is needed ...
            ! Compute distance between grid point and the center of the simulation cell in R space 
            !
            dq1=h(1,1)*dqs1+h(1,2)*dqs2+h(1,3)*dqs3   !r_i = h s_i
            dq2=h(2,1)*dqs1+h(2,2)*dqs2+h(2,3)*dqs3   !r_i = h s_i
            dq3=h(3,1)*dqs1+h(3,2)*dqs2+h(3,3)*dqs3   !r_i = h s_i
            !
            dist = DSQRT(dq1*dq1+dq2*dq2+dq3*dq3)
            !-------------------------------------------------
            me_cs(1,i,j,k)=dq1
            me_cs(2,i,j,k)=dq2
            me_cs(3,i,j,k)=dq3
            !-------------------------------------------------
            me_rs(0,i,j,k)=1.0_DP
            me_rs(1,i,j,k)=dist
            !-------------------------------------------------
            me_ri(0,i,j,k)=1.0_DP
            !
            IF(dist.GE.1.0E-10) THEN
              me_ri(1,i,j,k)=1.0_DP/dist
            ELSE
              me_ri(1,i,j,k)=0.0_DP ! JJ: this seems wrong, but is correct
            END IF
            !-------------------------------------------------
            DO l=2,lmax
              me_rs(l,i,j,k)=me_rs(l-1,i,j,k)*me_rs(1,i,j,k)
            END DO
            !-------------------------------------------------
            DO l=2,lmax+1
              me_ri(l,i,j,k)=me_ri(l-1,i,j,k)*me_ri(1,i,j,k)
            END DO
            !-------------------------------------------------
            IF( (i.GE.s_me_r1).AND.(i.LE.s_me_r4).AND. &
                (j.GE.s_me_r2).AND.(j.LE.s_me_r5).AND. &
                (k.GE.s_me_r3).AND.(k.LE.s_me_r6) ) THEN
              !-------------------------------------------------
              dxy      = DSQRT(dq1*dq1+dq2*dq2)
              !-------------------------------------------------
              me_rc(0,i,j,k)     =1.0_DP
              me_rc(1:lmax,i,j,k)=0.0_DP
              !-------------------------------------------------
              IF (dxy .GT. 1.0E-10) THEN
                !-----------------------------------------------
                cxy = CMPLX(dq1,dq2)/dxy
                !-----------------------------------------------
                DO m=1,lmax
                  me_rc(m,i,j,k)=me_rc(m-1,i,j,k)*cxy
                END DO
                !-----------------------------------------------
              END IF
              !-------------------------------------------------
            END IF
            !-------------------------------------------------
          END DO
        END DO
      END DO
      !---------------------------------------------------------------------------------------------
      return
    end subroutine exx_gs_setup_cube

    subroutine  exx_gs_setup_sphere()
      implicit none
      !========================================================================
      ! Compute distances between grid points and the center of the simulation cell in R space 
      ! This part needs to be done once in constant volume simulation and
      ! needs to be done every step in variable cell simulationulations ...
      !
      !========================================================================
      DO i=1,np_in_sp_me_s
        xx_in_sp(i)=h(1,1)*sc_xx_in_sp(i)+h(1,2)*sc_yy_in_sp(i)+h(1,3)*sc_zz_in_sp(i)   ! r = h s
        yy_in_sp(i)=h(2,1)*sc_xx_in_sp(i)+h(2,2)*sc_yy_in_sp(i)+h(2,3)*sc_zz_in_sp(i)   ! r = h s
        zz_in_sp(i)=h(3,1)*sc_xx_in_sp(i)+h(3,2)*sc_yy_in_sp(i)+h(3,3)*sc_zz_in_sp(i)   ! r = h s
      END DO
      return
    end subroutine exx_gs_setup_sphere

    subroutine  solve_a_nonself_pair_cube()
      implicit none
      call generate_pe_guess_cube
      call start_clock('exx_grid_trans')
      CALL getmiddlewc(wannierc(1,gindex_of_iobtl),wannierc(1,my_var2), h, ainv, middle )
      !
      ! calculate translation vector from the center of the box
      CALL getsftv( nr1s, nr2s, nr3s, h, ainv, middle, tran)
      call stop_clock('exx_grid_trans')
      !      
      ! get the localized psi around the mid point of two wannier centers
      ! note: the psime is centered at the center of the box
      ! (using the translation vector "tran" from middle of wfc to the center of box)
      call start_clock('exx_psicb')
      CALL getpsicb( nrg, p_me_r, psi(1,iobtl), psime(1), tran)
#if ! defined(__MPI)
      ! does not need communication so construct it here; my_var2 is simultaneously the global and local index (serial case)
      CALL getpsicb( nrg, p_me_r, psi(1,my_var2), psime_pair_recv(1, j, iobtl), tran)
#endif
      call stop_clock('exx_psicb')
      ! 
      ! the localized density rhome 
      call start_clock('exx_getrhol')
      CALL getrhol_cube(p_me_r, p_ps_r, psime(1), psime_pair_recv(1, j, iobtl), rhome, rhops, inv_omega)
      call stop_clock('exx_getrhol')
      !
      ! calculate the exx potential from the pair density by solving Poisson
      !
      !--------------------------------------------------------------------------------------
      call start_clock('exx_vofr')
      CALL getvofr_cube( p_me_r, p_ps_r, max(n_p_me,n_s_me), max(n_p_ps,n_s_ps), hcub, &
        rhops, potme, pair_status(pos, iobtl), psgsn, pairrho(:,:,pos,iobtl), &
        pairv(:,:,pos,iobtl), cgstep)
      call stop_clock('exx_vofr')
      !--------------------------------------------------------------------------------------
      !
      !--------------------------------------------------------------------------------------
      ! write cgsteps in the suffix.ncg (unit=44)
      !--------------------------------------------------------------------------------------
      IF ((MOD(nfi,iprint_stdout).EQ.0)) THEN
        IF (ionode) THEN ! maybe not needed for ionode (if one want more information)
          WRITE(iunit,'(3X,"(i,j,cgsteps)",3I6)') gindex_of_iobtl, my_var2, cgstep
        END IF    
      END IF    
      !--------------------------------------------------------------------------------------
      ! 
      !--------------------------------------------------------------------------------------
      ! update force and energy
      !--------------------------------------------------------------------------------------
      call start_clock('exx_force_loc')
      CALL updateforce_loc(nrg, p_me_r, vpsil(:,iobtl), potme, psime, psime_pair_recv(1,j,iobtl),tran)
#if ! defined(__MPI)
      ! does not need communication so construct it here; my_var2 is simultaneously the global and local index (serial case)
      CALL updateforce_loc(nrg, p_me_r, vpsil(:,my_var2), potme, psime_pair_recv(1,j,iobtl), psime, tran)
#endif
      call stop_clock('exx_force_loc')
      !
      call start_clock('exx_penergy')
      CALL vvprod_cube(p_me_r, rhome, potme, paire(j))    ! dot product of the rho and potme  
      call stop_clock('exx_penergy')
      return
    end subroutine solve_a_nonself_pair_cube

    subroutine  solve_a_nonself_pair_sphere()
      implicit none
      call start_clock('exx_grid_trans')
      CALL getmiddlewc(wannierc(1,gindex_of_iobtl),wannierc(1,my_var2), h, ainv, middle )
      ! d_pair is used in the extrapolation scheme (sphere only)
      CALL get_pair_dist(wannierc(1,gindex_of_iobtl),wannierc(1,my_var2),d_pair)
      !
      ! calculate translation vector from the center of the box
      CALL getsftv( nr1s, nr2s, nr3s, h, ainv, middle, tran)
      call stop_clock('exx_grid_trans')
      !      
      ! get the localized psi around the mid point of two wannier centers
      ! note: the psil is centered at the center of the box
      ! (using the translation vector "tran" from middle of wfc to the center of box)
      ALLOCATE ( psil(np_in_sp_me_p) ); psil=0.0_DP ! HK/MCA: (TODO) the allocation can to be done in a reusable fashion
      CALL getpsil( nnrtot, np_in_sp_me_p, psi(1, iobtl), psil(1), tran) 
#if ! defined(__MPI)
      ! does not need communication so construct it here; my_var2 is simultaneously the global and local index (serial case)
      CALL getpsil( nnrtot, np_in_sp_me_p, psi(1, my_var2), psime_pair_recv(1, j, iobtl), tran)
#endif
      !
      ! the localized density rhol 
      ! HK/MCA: (TODO) need to make these array allocation in a reusable fashion
      ALLOCATE ( rhol(np_in_sp_me_p) ); rhol=0.0_DP
      ALLOCATE ( rho_in_sp(np_in_sp_p) ); rho_in_sp=0.0_DP
      CALL getrhol_sphere( np_in_sp_me_p, np_in_sp_p, psil(1), psime_pair_recv(1, j, iobtl), rhol, rho_in_sp, tran, sa1)
      !
      ! calculate the exx potential from the pair density by solving Poisson
      !
      ! calculate the exx potential from the pair density by solving Poisson
      !
      ALLOCATE ( vl(np_in_sp_me_p) ); vl=0.0_DP ! compute potential (vl) in ME sphere 
      CALL start_clock('getvofr')
      ! HK/MCA: d_pair is used in the extrapolation scheme (check if still working) see also ``call get_pair_dist''
      CALL getvofr_sphere( np_in_sp_me_p, np_in_sp_p, hcub, rho_in_sp, vl,&
        pairv(1,1,j,iobtl), pairv(1,2,j,iobtl), pairv(1,3,j,iobtl),&
        .FALSE., d_pair, pair_dist(1,j,iobtl), pair_dist(2,j,iobtl),&
        pair_dist(3,j,iobtl),cgstep)
      CALL stop_clock('getvofr')
      !
      ! write cgsteps in the suffix.ncg (unit=44)
      !
      IF ((MOD(nfi,iprint_stdout).EQ.0)) THEN
        !   
        IF (ionode) THEN ! maybe not needed for ionode (if one want more information)
          !   
          WRITE(iunit,'(3X,"(i,j,cgsteps)",3I6)') gindex_of_iobtl, my_var2, cgstep
          !
        END IF    
        !   
      END IF    
      ! 
      ! update vpsil in the global grid (exxalfa is 0.25 for PBE0)
      !$omp parallel do private(ir) 
      DO ip = 1, np_in_sp_me_p
        CALL l2goff (ip,ir,tran) ! local is centered at box center; global index is offset by tran
        vpsil(ir,iobtl) = vpsil(ir,iobtl) - exxalfa*vl(ip)*psime_pair_recv(ip,j,iobtl) ! to remain 
#if defined(__MPI)
        psime_pair_recv(ip,j,iobtl) =            - exxalfa*vl(ip)*psil(ip)             ! to be sent 
#else
        ! does not need communication so construct it here; my_var2 is simultaneously the global and local index (serial case)
        vpsil(ir,my_var2) = vpsil(ir,my_var2)    - exxalfa*vl(ip)*psil(ip)
#endif
      END DO
      !$omp end parallel do
      !
      CALL vvprod_sphere(np_in_sp_me_p, rhol, vl, paire(j)) ! dot product of the rho and vl !HK (todo): do we need to do PS+ME ?? rho_in_sp may be enough
      return
    end subroutine solve_a_nonself_pair_sphere

    subroutine  generate_pe_guess_cube()
      implicit none
      guess_status = 1        ! guess what?
      pos = 0
      oldest_step = 100000000
      DO itr = 1, neigh
        IF ( pair_label(itr, iobtl) .EQ. 0 ) THEN
          pos = itr
          EXIT
        ELSE IF ( pair_label(itr, iobtl) .EQ. my_var2 ) THEN
          guess_status = 1
          pos = itr
          EXIT
        ELSE IF ( pair_step(itr, iobtl) < oldest_step ) THEN
          pos = itr
          oldest_step = pair_step(itr, iobtl)
        END IF
      END DO
      !
      IF (guess_status .EQ. 1) THEN
        IF (pair_step(pos, iobtl) .EQ. n_exx - 1) THEN
          pair_status(pos, iobtl) = pair_status(pos, iobtl) + 1
        ELSE
          pair_status(pos, iobtl) = 1
        END IF
      ELSE
        pair_status(pos, iobtl) = 0
      END IF
      !
      pair_label(pos, iobtl) = my_var2
      pair_step(pos, iobtl)  = n_exx
      return
    end subroutine generate_pe_guess_cube

    subroutine  solve_a_self_pair_cube()
      implicit none
      ! calculate translation vector from the center of the box
      CALL getsftv(nr1s, nr2s, nr3s, h, ainv, wannierc(1, gindex_of_iobtl), tran)
      !
      ! get the localized psi around the wannier centers
      ! note: the psime is centered at the center of the box
      ! (using the translation vector "tran" from the wfc to the center of box)
      !
      CALL getpsicb( nrg, s_me_r, psi(1,iobtl), psime(1), tran)
      ! get the localized density rhome  
      CALL getrhol_cube(s_me_r, s_ps_r, psime(1), psime(1), rhome, rhops, inv_omega)
      ! 
      ! calculate the exx potential from the pair density by solving Poisson
      !
      !--------------------------------------------------------------------------------------
      CALL start_clock('getvofr')
      !
      CALL getvofr_cube( s_me_r, s_ps_r, n_s_me, n_s_ps, hcub, rhops, potme, n_exx-1, psgsn, &
        selfrho(:,:,iobtl), selfv(:,:,iobtl), cgstep)
      !
      CALL stop_clock('getvofr')
      !--------------------------------------------------------------------------------------
      !
      !--------------------------------------------------------------------------------------
      ! write cgsteps in the suffix.ncg (unit=44)
      !--------------------------------------------------------------------------------------
      IF ((MOD(nfi,iprint_stdout).EQ.0)) THEN
        IF (ionode) THEN ! maybe not needed for ionode (if one want more information)
          WRITE(44,'(3X,"(i,i,cgsteps)",3I6)') gindex_of_iobtl, gindex_of_iobtl, cgstep
          CALL printout_base_close("ncg")
        END IF    
      END IF    
      !--------------------------------------------------------------------------------------
      !
      !--------------------------------------------------------------------------------------
      ! update force and energy
      !--------------------------------------------------------------------------------------
      !
      CALL updateforce_slf(nrg, s_me_r, vpsil(1,iobtl), potme, psime, tran)
      CALL vvprod_cube(s_me_r, rhome, potme, selfe)    ! dot product of the rho and potme 
      return
    end subroutine solve_a_self_pair_cube

    subroutine  solve_a_self_pair_sphere()
      implicit none
      ! calculate translation vector from the center of the box
      CALL getsftv(nr1s, nr2s, nr3s, h, ainv, wannierc(1, gindex_of_iobtl), tran)
      !
      ! get the localized psi around the wannier centers
      ! note: the psil is centered at the center of the box
      ! (using the translation vector "tran" from the wfc to the center of box)
      ALLOCATE ( psil(np_in_sp_me_s) ); psil=0.0_DP
      CALL getpsil( nnrtot, np_in_sp_me_s, psi(1, iobtl), psil(1), tran)
      !
      ! get the localized density rhol  
      ALLOCATE ( rhol(np_in_sp_me_s) ); rhol=0.0_DP
      ALLOCATE ( rho_in_sp(np_in_sp_s) ); rho_in_sp=0.0_DP
      CALL getrhol_sphere( np_in_sp_me_s, np_in_sp_s, psil(1), psil(1), rhol, rho_in_sp, tran, sa1)
      ! 
      ! calculate the exx potential from the pair density by solving Poisson
      !
      ALLOCATE ( vl(np_in_sp_me_s) ); vl=0.0_DP ! compute potential (vl) in ME sphere 
      CALL start_clock('getvofr')
      CALL getvofr_sphere( np_in_sp_me_s,np_in_sp_s,&
        hcub, rho_in_sp, vl, selfv(1,1,iobtl), selfv(1,2,iobtl),&
        selfv(1,3,iobtl), .TRUE., 0.0, 0.0, 0.0, 0.0,cgstep)
      !
      CALL stop_clock('getvofr')
      !
      ! write cgsteps in the suffix.ncg (unit=44)
      !
      IF ((MOD(nfi,iprint_stdout).EQ.0)) THEN
        !   
        IF (ionode) THEN ! maybe not needed for ionode (if one want more information)
          !   
          WRITE(44,'(3X,"(i,i,cgsteps)",3I6)') gindex_of_iobtl, gindex_of_iobtl, cgstep
          !
          CALL printout_base_close("ncg")
          !   
        END IF    
        !   
      END IF    
      !
      ! update vpsil in the global grid (exxalfa is 0.25 for PBE0) 
      !$omp parallel do private(ir) 
      DO ip = 1, np_in_sp_me_s 
        CALL l2goff (ip,ir,tran) ! local is centered at box center; global index is offset by tran
        vpsil(ir,iobtl) = vpsil(ir,iobtl) - exxalfa*vl(ip)*psil(ip) ! PBE0
      END DO
      !$omp end parallel do 
      !
      ! compute exchange energy in ME sphere 
      CALL vvprod_sphere(np_in_sp_me_s, rhol, vl, selfe) ! dot product of the rho and vl !HK (todo): do we need to do PS+ME ?? rho_in_sp may be enough
      return
    end subroutine solve_a_self_pair_sphere

END SUBROUTINE exx_gs
!====================================================================================

!==============================================================================
SUBROUTINE getsftv(nr1s, nr2s, nr3s, h, ainv, wc, tran)
    !
    USE kinds, ONLY  : DP
    !
    IMPLICIT NONE
    !
    INTEGER  nr1s, nr2s, nr3s, tran(3)
    REAL(DP) wc(3), wcs(3), h(3,3), ainv(3,3)
    !
    INTEGER  i, bcm(3), wcm(3)
    !
    bcm(1) = INT(nr1s/2)
    bcm(2) = INT(nr2s/2)
    bcm(3) = INT(nr3s/2)
    !
    wcs(1)=ainv(1,1)*wc(1)+ainv(1,2)*wc(2)+ainv(1,3)*wc(3)   ! s_ij = h^-1 r_ij
    wcs(2)=ainv(2,1)*wc(1)+ainv(2,2)*wc(2)+ainv(2,3)*wc(3)   ! s_ij = h^-1 r_ij
    wcs(3)=ainv(3,1)*wc(1)+ainv(3,2)*wc(2)+ainv(3,3)*wc(3)   ! s_ij = h^-1 r_ij
    !
    ! we map to nearest grid point along each lattice vector
    !
    wcm(1) = IDNINT( wcs(1)*DBLE(nr1s) )  
    wcm(2) = IDNINT( wcs(2)*DBLE(nr2s) )
    wcm(3) = IDNINT( wcs(3)*DBLE(nr3s) )
    !
    DO i = 1, 3
      tran(i) = bcm(i) - wcm(i)
    ENDDO
    !
    RETURN
END SUBROUTINE getsftv
!==============================================================================

!==============================================================================
SUBROUTINE getmiddlewc(wc1, wc2, h, ainv, mid)
    !
    USE kinds, ONLY  : DP
    !
    IMPLICIT NONE
    !
    REAL(DP)  wc1(3), wc2(3), mid(3)
    REAL(DP)  dij(3), dij2(3), h(3,3), ainv(3,3), tmp(3)
    !
    INTEGER   i
    !
    dij(1)=wc2(1)-wc1(1)   ! r_ij = r_i - r_j   
    dij(2)=wc2(2)-wc1(2)   ! r_ij = r_i - r_j   
    dij(3)=wc2(3)-wc1(3)   ! r_ij = r_i - r_j   
    !
    dij2(1)=ainv(1,1)*dij(1)+ainv(1,2)*dij(2)+ainv(1,3)*dij(3)   ! s_ij = h^-1 r_ij
    dij2(2)=ainv(2,1)*dij(1)+ainv(2,2)*dij(2)+ainv(2,3)*dij(3)   ! s_ij = h^-1 r_ij
    dij2(3)=ainv(3,1)*dij(1)+ainv(3,2)*dij(2)+ainv(3,3)*dij(3)   ! s_ij = h^-1 r_ij
    !
    dij2(1)=dij2(1)-IDNINT(dij2(1))   ! impose MIC on s_ij in range: [-0.5,+0.5]
    dij2(2)=dij2(2)-IDNINT(dij2(2))   ! impose MIC on s_ij in range: [-0.5,+0.5]
    dij2(3)=dij2(3)-IDNINT(dij2(3))   ! impose MIC on s_ij in range: [-0.5,+0.5]
    !
    dij(1)=h(1,1)*dij2(1)+h(1,2)*dij2(2)+h(1,3)*dij2(3)   ! r_ij = h s_ij (MIC)
    dij(2)=h(2,1)*dij2(1)+h(2,2)*dij2(2)+h(2,3)*dij2(3)   ! r_ij = h s_ij (MIC)
    dij(3)=h(3,1)*dij2(1)+h(3,2)*dij2(2)+h(3,3)*dij2(3)   ! r_ij = h s_ij (MIC)
    !
    mid(1) = wc1(1) + 0.5_DP*dij(1)
    mid(2) = wc1(2) + 0.5_DP*dij(2)
    mid(3) = wc1(3) + 0.5_DP*dij(3)
    !
    RETURN
END SUBROUTINE getmiddlewc
!==============================================================================

!==============================================================================
SUBROUTINE get_pair_dist (wc1, wc2, P_dist)
    !
    USE kinds,                   ONLY  : DP
    USE cell_base,               ONLY  : omega, ainv, h
    !
    IMPLICIT NONE
    ! P_dist is the pair distance with in the minimum image convention
    REAL(DP)  wc1(3), wc2(3), P_dist
    REAL(DP)  dist_vec_in_r(3), dist_vec_in_s(3)
    !
    INTEGER i
    !
    P_dist = 0.D0
    ! dist vector
    DO i = 1,3
      dist_vec_in_r(i) = wc1(i) - wc2(i)
    END DO
    !
    dist_vec_in_s(1)=ainv(1,1)*dist_vec_in_r(1)+ainv(1,2)*dist_vec_in_r(2)+ainv(1,3)*dist_vec_in_r(3)   ! s = h^-1 r
    dist_vec_in_s(2)=ainv(2,1)*dist_vec_in_r(1)+ainv(2,2)*dist_vec_in_r(2)+ainv(2,3)*dist_vec_in_r(3)   ! s = h^-1 r
    dist_vec_in_s(3)=ainv(3,1)*dist_vec_in_r(1)+ainv(3,2)*dist_vec_in_r(2)+ainv(3,3)*dist_vec_in_r(3)   ! s = h^-1 r
    !
    DO i = 1,3
      dist_vec_in_s(i)=dist_vec_in_s(i)-IDNINT(dist_vec_in_s(i))   
    END DO
    !
    dist_vec_in_r(1)=h(1,1)*dist_vec_in_s(1)+h(1,2)*dist_vec_in_s(2)+h(1,3)*dist_vec_in_s(3)   ! r = h s
    dist_vec_in_r(2)=h(2,1)*dist_vec_in_s(1)+h(2,2)*dist_vec_in_s(2)+h(2,3)*dist_vec_in_s(3)   ! r = h s
    dist_vec_in_r(3)=h(3,1)*dist_vec_in_s(1)+h(3,2)*dist_vec_in_s(2)+h(3,3)*dist_vec_in_s(3)   ! r = h s
    !
    DO i = 1, 3
      P_dist = P_dist + dist_vec_in_r(i)*dist_vec_in_r(i)  
    END DO
    !
    P_dist = DSQRT(P_dist)
    !
    RETURN
END SUBROUTINE get_pair_dist
!==============================================================================

!==============================================================================
SUBROUTINE vvprod_sphere(n, v1, v2, prod)
    !
    USE kinds, ONLY  : DP
    !
    IMPLICIT NONE
    !
    INTEGER  n
    REAL(DP) prod, v1(n), v2(n), vp
    !
    INTEGER  i
    !
    prod = 0.0_DP
    vp = 0.0_DP
    !
    !$omp parallel do reduction(+:vp)
    DO i = 1, n
      vp = vp + v1(i) * v2(i)
    END DO
    !$omp end parallel do 
    !
    prod = vp
    !
    RETURN
END SUBROUTINE vvprod_sphere
!==============================================================================

!==============================================================================
SUBROUTINE vvprod_cube(me_r, v1, v2, prod)
    !
    USE kinds, ONLY  : DP
    !
    IMPLICIT NONE
    !
    !----------------------------------------------------------------
    INTEGER      :: me_r(6)
    REAL(DP)     :: v1(me_r(1):me_r(4),me_r(2):me_r(5),me_r(3):me_r(6))
    REAL(DP)     :: v2(me_r(1):me_r(4),me_r(2):me_r(5),me_r(3):me_r(6))
    REAL(DP)     :: prod
    !----------------------------------------------------------------
    INTEGER      :: i,j,k
    REAL(DP)     :: prodp
    !----------------------------------------------------------------
    !
    prodp=0.0D0
    !
    ! WRITE(*,*) "vvprod"
    !
    !$omp parallel do private(i,j,k) reduction(+:prodp)
    DO k=me_r(3),me_r(6)
      DO j=me_r(2),me_r(5)
        DO i=me_r(1),me_r(4)
          !----------------------------------------------------------
          prodp = prodp + v1(i,j,k)*v2(i,j,k)
          !----------------------------------------------------------
          ! WRITE(*,"(I4,I4,I4,F15.11,F15.11,F15.11)") i,j,k,v1(i,j,k),v2(i,j,k),v1(i,j,k)*v2(i,j,k)
          !----------------------------------------------------------
        END DO
      END DO
    END DO
    !----------------------------------------------------------------
    !$omp end parallel do 
    !
    prod = prodp
    !
    RETURN
END SUBROUTINE vvprod_cube
!==============================================================================

!==============================================================================
SUBROUTINE getpsil( ntot, np_in_sp_me, psi, psi2, tran)
    !
    USE kinds, ONLY  : DP
    !
    IMPLICIT NONE
    !
    INTEGER  ntot, tran(3), np_in_sp_me
    REAL(DP) psi(ntot), psi2(np_in_sp_me)
    !
    INTEGER  ir, ip, i, j, k, ii, jj, kk
    !
    !$omp parallel do private(ir) 
    DO ip = 1, np_in_sp_me 
      CALL l2goff (ip,ir,tran)
      psi2(ip) = psi(ir)
    END DO
    !$omp end parallel do 
    !
    RETURN
END SUBROUTINE getpsil
!==============================================================================

!==============================================================================
SUBROUTINE getrhol_cube(me_r, ps_r, psi1, psi2, rhome, rhops, inv_omega)
    !
    USE kinds, ONLY  : DP
    !
    IMPLICIT NONE
    !
    INTEGER  me_r(6), ps_r(6)
    REAL(DP) psi1 (me_r(1):me_r(4),me_r(2):me_r(5),me_r(3):me_r(6))
    REAL(DP) psi2 (me_r(1):me_r(4),me_r(2):me_r(5),me_r(3):me_r(6))
    REAL(DP) rhome(me_r(1):me_r(4),me_r(2):me_r(5),me_r(3):me_r(6))
    REAL(DP) rhops(ps_r(1):ps_r(4),ps_r(2):ps_r(5),ps_r(3):ps_r(6))
    REAL(DP) inv_omega
    !
    INTEGER  i, j, k
    !
    !$omp parallel do private(i,j,k)
    DO k=ps_r(3),ps_r(6)
      DO j=ps_r(2),ps_r(5)
        DO i=ps_r(1),ps_r(4)
          rhops(i,j,k) = psi1(i,j,k) * psi2(i,j,k) * inv_omega
        END DO
      END DO
    END DO
    !$omp end parallel do 
    !--------------------------------------------------------------------------
    !
    !--------------------------------------------------------------------------
    !$omp parallel do private(i,j,k)
    DO k=me_r(3),me_r(6)
      DO j=me_r(2),me_r(5)
        DO i=me_r(1),me_r(4)
          rhome(i,j,k) = psi1(i,j,k)*psi2(i,j,k) * inv_omega
        END DO
      END DO
    END DO
    !$omp end parallel do 
    !--------------------------------------------------------------------------
    !
    RETURN
END SUBROUTINE getrhol_cube
!==============================================================================

!==============================================================================
SUBROUTINE getrhol_sphere( np_in_sp_me, np_in_sp, psi, psi2, rho, rho_in_sp, tran, sa1)
    !
    USE kinds, ONLY  : DP
    !
    IMPLICIT NONE
    !
    INTEGER  np_in_sp_me, tran(3),np_in_sp
    REAL(DP) psi(np_in_sp_me), psi2(np_in_sp_me), rho(np_in_sp_me),sa1, rho_in_sp(np_in_sp)
    !
    INTEGER  ir, ip, i, j, k, ii, jj, kk
    rho_in_sp(:) = 0.D0
    !
    !$omp parallel do 
    DO ip = 1, np_in_sp_me 
      rho(ip) = psi(ip) * psi2(ip) * sa1
#if ! defined(__PGI)
! workaround for a likely PGI bug: for no apparent reasons,
! this loop crashes with PGI v.19.10 and 20.9 - PG
      IF( ip.LE.np_in_sp ) THEN
        rho_in_sp( ip ) = rho(ip)
      END IF
#endif
    ENDDO
    !$omp end parallel do 
#if defined(__PGI)
    DO ip = 1, np_in_sp
      rho_in_sp( ip ) = rho( ip )
    ENDDO
#endif
    !
    RETURN
END SUBROUTINE getrhol_sphere
!==============================================================================

!==============================================================================
SUBROUTINE l2goff (lind,gind,tran)
    !
    USE exx_module,       ONLY  : odtothd_in_sp, thdtood
    USE fft_base,         ONLY  : dfftp
    !
    IMPLICIT NONE
    !
    INTEGER  tran(3),lind,gind
    INTEGER  ir, ip, i, j, k, ii, jj, kk, nr1s, nr2s, nr3s
    !
    nr1s=dfftp%nr1; nr2s=dfftp%nr2; nr3s=dfftp%nr3 
    !
    i  = odtothd_in_sp(1, lind)
    j  = odtothd_in_sp(2, lind)
    k  = odtothd_in_sp(3, lind)
    !
    ii = i - tran(1)
    jj = j - tran(2)
    kk = k - tran(3)
    !
    IF ( ii .GT. nr1s)ii = ii - nr1s
    IF ( jj .GT. nr2s)jj = jj - nr2s
    IF ( kk .GT. nr3s)kk = kk - nr3s
    !
    IF ( ii .LT. 1)ii = ii + nr1s
    IF ( jj .LT. 1)jj = jj + nr2s
    IF ( kk .LT. 1)kk = kk + nr3s
    !
    gind = thdtood(ii, jj, kk)
    !
    RETURN
END SUBROUTINE l2goff
!==============================================================================
SUBROUTINE getpsicb(nrg,nrl,psig,psil,tran)
    !
    USE kinds, ONLY  : DP
    USE fft_base,         ONLY  : dfftp
    !
    IMPLICIT NONE
    !
    INTEGER      :: nrg(3)
    INTEGER      :: nrl(6)
    REAL(DP)     :: psig(nrg(1),nrg(2),nrg(3))
    REAL(DP)     :: psil(nrl(1):nrl(4),nrl(2):nrl(5),nrl(3):nrl(6))
    INTEGER      :: tran(3)
    INTEGER      :: gid(3)
    !INTEGER      :: lid(3)
    INTEGER      :: i,j,k
    INTEGER      :: ti, tj, tk
    INTEGER      :: gi, gj, gk
    integer, external :: l2gcb
    !
    ti = tran(1); tj = tran(2); tk = tran(3)
    !$omp parallel do private(i,j,k,gi,gj,gk)
    DO k = nrl(3),nrl(6)
      DO j = nrl(2),nrl(5)
        DO i = nrl(1),nrl(4)
          !----------------------------------------------------
          gi = l2gcb(dfftp%nr1,i,ti)
          gj = l2gcb(dfftp%nr2,j,tj)
          gk = l2gcb(dfftp%nr3,k,tk)
          psil(i,j,k)=psig(gi,gj,gk)
          !----------------------------------------------------
        END DO
      END DO
    END DO
    !$omp end parallel do 
    !----------------------------------------------------------
    !
    RETURN
END SUBROUTINE getpsicb
!==============================================================================

integer function l2gcb(n,l,t)
  implicit none
  integer :: n, l, t
  l2gcb = MOD(l-t-1+n, n)+1
end function l2gcb
SUBROUTINE updateforce_loc(nrg, me_r, vpsil, potme, psime1, psime2, tran)
    !
    USE kinds,                   ONLY  : DP
    USE exx_module,              ONLY  : exxalfa
    USE fft_base,         ONLY  : dfftp
    !
    IMPLICIT NONE
    !
    INTEGER      :: nrg(3)
    INTEGER      :: me_r(6)
    REAL(DP)     :: vpsil(nrg(1),nrg(2),nrg(3))
    REAL(DP)     :: potme(me_r(1):me_r(4),me_r(2):me_r(5),me_r(3):me_r(6))
    REAL(DP)     :: psime1(me_r(1):me_r(4),me_r(2):me_r(5),me_r(3):me_r(6))
    REAL(DP)     :: psime2(me_r(1):me_r(4),me_r(2):me_r(5),me_r(3):me_r(6))
    INTEGER      :: tran(3)
    !----------------------------------------------------------------
    INTEGER      :: gid(3)
    INTEGER      :: lid(3)
    INTEGER      :: i,j,k
    INTEGER      :: gi,gj,gk,ti,tj,tk
    integer, external :: l2gcb
    !
    ti=tran(1);tj=tran(2);tk=tran(3)
    !----------------------------------------------------------------
    ! update vpsil in the global grid (exxalfa is 0.25 for PBE0)
    !----------------------------------------------------------------
    !$omp parallel do private(i,j,k,gi,gj,gk)
    !----------------------------------------------------------------
    DO k=me_r(3),me_r(6)
      DO j=me_r(2),me_r(5)
        DO i=me_r(1),me_r(4)
          !----------------------------------------------------------
          !CALL l2gcb(lid,gid,tran)
          gi = l2gcb(dfftp%nr1,i,ti)
          gj = l2gcb(dfftp%nr2,j,tj)
          gk = l2gcb(dfftp%nr3,k,tk)
          !----------------------------------------------------------
          vpsil(gi,gj,gk) = vpsil(gi,gj,gk) &
                                         - exxalfa*potme(i,j,k)*psime2(i,j,k)
#if defined(__MPI)
          psime2(i,j,k) = - exxalfa*potme(i,j,k)*psime1(i,j,k)
#endif
          !----------------------------------------------------------
        END DO
      END DO
    END DO
    !----------------------------------------------------------------
    !$omp end parallel do 
    !----------------------------------------------------------------

    !----------------------------------------------------------------
    RETURN
    !----------------------------------------------------------------
END SUBROUTINE updateforce_loc
!==============================================================================

!==============================================================================
SUBROUTINE updateforce_slf(nrg, me_r, vpsil, potme, psime, tran)
    !
    USE kinds,                   ONLY  : DP
    USE exx_module,              ONLY  : exxalfa
    USE fft_base,         ONLY  : dfftp
    !
    IMPLICIT NONE
    !
    INTEGER      :: nrg(3)
    INTEGER      :: me_r(6)
    REAL(DP)     :: vpsil(nrg(1),nrg(2),nrg(3))
    REAL(DP)     :: potme(me_r(1):me_r(4),me_r(2):me_r(5),me_r(3):me_r(6))
    REAL(DP)     :: psime(me_r(1):me_r(4),me_r(2):me_r(5),me_r(3):me_r(6))
    INTEGER      :: tran(3)
    !----------------------------------------------------------------
    INTEGER      :: gid(3)
    INTEGER      :: lid(3)
    INTEGER      :: i,j,k
    !
    INTEGER      :: gi,gj,gk,ti,tj,tk
    integer, external :: l2gcb
    !
    ti=tran(1);tj=tran(2);tk=tran(3)
    !----------------------------------------------------------------
    ! update vpsil in the global grid (exxalfa is 0.25 for PBE0)
    !----------------------------------------------------------------
    !$omp parallel do private(i,j,k,gi,gj,gk)
    DO k=me_r(3),me_r(6)
      DO j=me_r(2),me_r(5)
        DO i=me_r(1),me_r(4)
          !----------------------------------------------------------
          gi = l2gcb(dfftp%nr1,i,ti)
          gj = l2gcb(dfftp%nr2,j,tj)
          gk = l2gcb(dfftp%nr3,k,tk)
          !----------------------------------------------------------
          vpsil(gi,gj,gk) = vpsil(gi,gj,gk) &
                                         - exxalfa*potme(i,j,k)*psime(i,j,k)
          !----------------------------------------------------------
        END DO
      END DO
    END DO
    !----------------------------------------------------------------
    !$omp end parallel do 
    !----------------------------------------------------------------

    !----------------------------------------------------------------
    RETURN
    !----------------------------------------------------------------
END SUBROUTINE updateforce_slf
!==============================================================================

!==============================================================================
SUBROUTINE updateforce_rec(nrg, me_r, vpsil, force, tran)
    !
    USE kinds,                   ONLY  : DP
    USE exx_module,              ONLY  : exxalfa
    USE fft_base,                ONLY  : dfftp
    !
    IMPLICIT NONE
    !
    INTEGER      :: nrg(3)
    INTEGER      :: me_r(6)
    REAL(DP)     :: vpsil(nrg(1),nrg(2),nrg(3))
    REAL(DP)     :: force(me_r(1):me_r(4),me_r(2):me_r(5),me_r(3):me_r(6))
    INTEGER      :: tran(3)
    !----------------------------------------------------------------
    INTEGER      :: gi, gj, gk
    INTEGER      :: i,j,k
    INTEGER      :: ti,tj,tk
    integer, external :: l2gcb
    !
    ti=tran(1);tj=tran(2);tk=tran(3)
    !
    !----------------------------------------------------------------
    ! update vpsil in the global grid (exxalfa is 0.25 for PBE0)
    !----------------------------------------------------------------
    !$omp parallel do private(i,j,k,gi,gj,gk)
    !----------------------------------------------------------------
    DO k=me_r(3),me_r(6)
      DO j=me_r(2),me_r(5)
        DO i=me_r(1),me_r(4)
          !----------------------------------------------------------
          gi = l2gcb(dfftp%nr1,i,ti)
          gj = l2gcb(dfftp%nr2,j,tj)
          gk = l2gcb(dfftp%nr3,k,tk)
          !----------------------------------------------------------
          vpsil(gi,gj,gk) = vpsil(gi,gj,gk) + force(i,j,k)
          !----------------------------------------------------------
        END DO
      END DO
    END DO
    !----------------------------------------------------------------
    !$omp end parallel do 
    !----------------------------------------------------------------

    !----------------------------------------------------------------
    RETURN
    !----------------------------------------------------------------
END SUBROUTINE updateforce_rec
!==============================================================================
