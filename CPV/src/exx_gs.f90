module dummy_exx
  implicit none
contains

#ifdef __CUDA
  attributes(host,device) &
#endif
  integer function l2gcb(n,l,t)
    implicit none
    integer, value :: n, l, t
    l2gcb = MOD(l-t-1+n, n)+1
  end function l2gcb

end module dummy_exx

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
    USE exx_module,              ONLY  : selfv, selfrho, pairv, pairrho, pair_dist
    USE exx_module,              ONLY  : pair_label, pair_status, pair_step
    USE exx_module,              ONLY  : exx_potential
    USE exx_module,              ONLY  : n_exx, sc_fac ! EXX step and the performance scaling factor
    USE exx_module,              ONLY  : nord1, nord2, fornberg
    USE exx_module,              ONLY  : nrg
    USE exx_module,              ONLY  : s_me_r, s_ps_r, p_me_r, p_ps_r
    USE exx_module,              ONLY  : n_s_me, n_s_ps, n_p_me, n_p_ps
    USE exx_module,              ONLY  : prev_obtl_recv
    USE exx_module,              ONLY  : lmax
    USE exx_module,              ONLY  : me_cs
    USE exx_module,              ONLY  : me_rs
    USE exx_module,              ONLY  : me_ri
    USE exx_module,              ONLY  : me_rc
#ifdef __CUDA
    USE exx_module,              ONLY  : me_cs_d
    USE exx_module,              ONLY  : me_rs_d
    USE exx_module,              ONLY  : me_ri_d
    USE exx_module,              ONLY  : me_rc_d
    !
    USE exx_module,              ONLY  : coe_1st_derv_d, coemicf_d, coeke_d  !MCA/HK : dirty hack for std CG
    !
    USE exx_module,              ONLY  : psi_d,potpsi_d  
#endif
    USE exx_module,              ONLY  : coe_1st_derv, coemicf, coeke !MCA/HK : dirty hack for std CG
    USE energies,                ONLY  : exx
    USE printout_base,           ONLY  : printout_base_open, printout_base_unit, printout_base_close
    USE wannier_base,            ONLY  : neigh, dis_cutoff
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
    USE control_flags,           ONLY  : iverbosity
    use exx_module, only : psime_pair_recv, psime_pair_send
#ifdef __CUDA
    use exx_module, only : psime_pair_recv_d, psime_pair_send_d
#endif
    !
    IMPLICIT NONE
    COMPLEX(DP)              :: c(ngw, nbspx)        ! wave functions at time t
#ifdef __MPI
    !
    INTEGER  :: istatus(MPI_STATUS_SIZE)
#endif
    INTEGER     i,j,k,l,m
    INTEGER     lid,gid
    INTEGER     ir, ip, nfi, ierr, nnrtot, nr1s,nr2s,nr3s
    INTEGER     nj_max, iunit
    INTEGER     s_me_r1,s_me_r2,s_me_r3,s_me_r4,s_me_r5,s_me_r6
    REAl(DP)    inv_omega, a(3), ha, hb, hc
    REAl(DP)    dist, dxy, costheta, sintheta
    COMPLEX(DP) cxy
    REAl(DP)    plm(0:lmax,0:lmax)
    REAl(DP)    ha_proj(3), hb_proj(3), hc_proj(3)
    REAl(DP)    hcub, centerx, centery, centerz
    REAL(DP)    middle(3)
    REAL(DP)    dq1,dq2,dq3
    REAL(DP)    dqs1,dqs2,dqs3
    REAL(DP)    d_pair            ! pair distance
    REAL(DP)    start_timer, stop_timer
    !
    REAL(DP),    ALLOCATABLE ::   potpsi(:,:)
#ifdef __CUDA
    attributes (pinned) :: potpsi
#endif
    REAL(DP),    ALLOCATABLE ::   psi(:,:)
    REAL(DP),    ALLOCATABLE ::   rhome(:),rhops(:),potme(:)
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
    INTEGER                  :: jj,ii,ia,ib,ic,my_var,my_var2,my_var3,i_fac,va,cgstep
    INTEGER                  :: pos, oldest_step, guess_status  
    INTEGER                  :: ndim,nogrp 
    INTEGER,    ALLOCATABLE  :: obtl_recv(:,:), obtl_send(:,:), num_recv(:), num_send(:)
    REAl(DP),   ALLOCATABLE  :: wannierc(:,:),wannierc_tmp(:,:)
    REAl(DP),   ALLOCATABLE  :: psime(:)
    REAL(DP),   ALLOCATABLE  :: exx_tmp(:,:),exx_tmp3(:,:)
    INTEGER,    ALLOCATABLE  :: sdispls(:), sendcount(:)
    INTEGER,    ALLOCATABLE  :: rdispls(:), recvcount(:)
    INTEGER,    ALLOCATABLE  :: sdispls1(:), sendcount1(:)
    INTEGER,    ALLOCATABLE  :: rdispls1(:), recvcount1(:)
    !
    REAL(DP) :: Jim(3,3)  ! jacobian [d/d x]        [d/d a]
    !                                |d/d y| = [J]  |d/d b|
    !                                [d/d z]        [d/d c]
    CHARACTER(LEN=256)       :: chr_me
    INTEGER                  :: ifile
    !
    ! CHARACTER (LEN=20) :: my_filename ! debug
    ! INTEGER            :: my_unit ! debug
    !INTEGER            :: psgsn=3 !MCA: number of steps that arrays are stored. Memory overhead from here?
    INTEGER            :: psgsn=1 !MCA: we are not using extrapolation anyways 
#ifdef __CUDA
    REAL(DP), ALLOCATABLE, DEVICE :: h_d(:,:)
    ATTRIBUTES (device) :: rhome,rhops,potme,psime
#endif
    !
     write(*,*) "entering exx_gs" ! debug
    !
    !=============================================================================================
    !
    ! WRITE (my_filename, "(A, I0.4)") "cg", me_image ! debug
    ! my_unit = me_image + 75332 ! debug
    ! OPEN (UNIT=my_unit, FILE=my_filename) ! debug
    ! 
    CALL start_clock('exx_gs_setup')
    !
    !---------------------------------------------------------------------------------------------
    ! make processor index start from 1
    !---------------------------------------------------------------------------------------------
    me = me_image+1
    !---------------------------------------------------------------------------------------------
    !
    !---------------------------------------------------------------------------------------------
    ! number of real space gird along each lattice parameter directions (a1, a2, a3)
    !---------------------------------------------------------------------------------------------
    nr1s=dfftp%nr1                                         ! fft dense grid in Lattice 1 direction
    nr2s=dfftp%nr2                                         ! fft dense grid in Lattice 2 direction
    nr3s=dfftp%nr3                                         ! fft dense grid in Lattice 3 direction
    !---------------------------------------------------------------------------------------------
    !
    !---------------------------------------------------------------------------------------------
    ! the length of each lattice parameters
    !---------------------------------------------------------------------------------------------
    a(1)=DSQRT(h(1,1)*h(1,1)+h(2,1)*h(2,1)+h(3,1)*h(3,1))                              ! lattice 1 
    a(2)=DSQRT(h(1,2)*h(1,2)+h(2,2)*h(2,2)+h(3,2)*h(3,2))                              ! lattice 2 
    a(3)=DSQRT(h(1,3)*h(1,3)+h(2,3)*h(2,3)+h(3,3)*h(3,3))                              ! lattice 3 
    !---------------------------------------------------------------------------------------------
    !
    !---------------------------------------------------------------------------------------------
    ! grid spacing in each lattice parameters
    !---------------------------------------------------------------------------------------------
    ha = a(1) / DBLE(nr1s)                                   ! grid spacing in Lattice 1 direction
    hb = a(2) / DBLE(nr2s)                                   ! grid spacing in Lattice 2 direction
    hc = a(3) / DBLE(nr3s)                                   ! grid spacing in Lattice 3 direction
    !---------------------------------------------------------------------------------------------
    !
    !---------------------------------------------------------------------------------------------
    ! consider a grid an unit cell, the lattice vectors for the grid
    !---------------------------------------------------------------------------------------------
    ha_proj(:) = ha*h(:,1)/a(1)
    hb_proj(:) = hb*h(:,2)/a(2)
    hc_proj(:) = hc*h(:,3)/a(3)
    !---------------------------------------------------------------------------------------------
    !
    !---------------------------------------------------------------------------------------------
    ! total number of real space grid points in the global mesh 
    ! and the corresponding volume elements for each grid point
    !---------------------------------------------------------------------------------------------
    nnrtot = nr1s * nr2s * nr3s
    hcub   = omega / DBLE(nnrtot)                                  ! volume for a single grid cell
    !
    ! laplacian prefactor initialization
    !call laplacian_setup_meta
    !---------------------------------------------------------------------------------------------
    !
    !---------------------------------------------------------------------------------------------
    ! the x,y,z coordinates of the center of the box (gird center)
    ! NOTE: center of the box is set to grid point at int(nr1/2), int(nr2/2), int(nr3/2)
    !---------------------------------------------------------------------------------------------
    centerx = h(1,1)*DBLE(INT(nr1s/2))+h(1,2)*DBLE(INT(nr2s/2))+h(1,3)*DBLE(INT(nr3s/2))
    centery = h(2,1)*DBLE(INT(nr1s/2))+h(2,2)*DBLE(INT(nr2s/2))+h(2,3)*DBLE(INT(nr3s/2))
    centerz = h(3,1)*DBLE(INT(nr1s/2))+h(3,2)*DBLE(INT(nr2s/2))+h(3,3)*DBLE(INT(nr3s/2))
    !---------------------------------------------------------------------------------------------
    !
    !---------------------------------------------------------------------------------------------
    ! inverse volume
    !---------------------------------------------------------------------------------------------
    inv_omega = 1.0_DP/omega
    !---------------------------------------------------------------------------------------------
    !
    !---------------------------------------------------------------------------------------------
    ! Compute distances between grid points and the center of the simulation cell in R space 
    ! This part needs to be done once in constant volume simulation and
    ! needs to be done every step in variable cell simulationulations ...
    !---------------------------------------------------------------------------------------------
    s_me_r1=s_me_r(1)
    s_me_r2=s_me_r(2)
    s_me_r3=s_me_r(3)
    s_me_r4=s_me_r(4)
    s_me_r5=s_me_r(5)
    s_me_r6=s_me_r(6)
#ifdef __CUDA
    ALLOCATE (h_d, source=h)
    associate (me_cs=>me_cs_d, me_rs=>me_rs_d, me_ri=>me_ri_d, me_rc=>me_rc_d, h=>h_d)
    !$cuf kernel do (3)
#endif
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
#ifdef __CUDA
   end associate
   me_cs = me_cs_d 
   me_rs = me_rs_d 
   me_ri = me_ri_d 
   me_rc = me_rc_d 
   DEALLOCATE(h_d)
#endif

    !---------------------------------------------------------------------------------------------
    ! Compute and renormalize numerical derivative coefficients
    ! This part needs to be done once in constant volume simulation and
    ! needs to be done every step in variable cell simulationulations ...
    !---------------------------------------------------------------------------------------------
    ! Get ha*d/da, ha^2*d^2/da^2 stencil and cross coefficients
    ! For the finite difference coefficients, see B. Fornberg in Math. Comp. 51 (1988), 699-706
    !---------------------------------------------------------------------------------------------
    CALL fornberg(nord1, nord2,coe_1st_derv(:,1),coeke(:,1,1),coeke(:,1,2),ierr)
    !---------------------------------------------------------------------------------------------
    IF (ierr .ne. 0) THEN
      WRITE(stdout,*) ' ERROR: Wrong parameter in CALL of Fornberg'
      WRITE(stdout,*) ' STOP in exx_gs'
      RETURN
    END IF
    !---------------------------------------------------------------------------------------------
    ! RENORMALIZE COEKES WITH RESPECT TO THE GRID SPACING
    ! First derivative coefficients
    !---------------------------------------------------------------------------------------------
    coe_1st_derv(:,3) = coe_1st_derv(:,1)/hc                                        ! d/dc stencil
    coe_1st_derv(:,2) = coe_1st_derv(:,1)/hb                                        ! d/db stencil
    coe_1st_derv(:,1) = coe_1st_derv(:,1)/ha                                        ! d/da stencil
    !---------------------------------------------------------------------------------------------
    ! NOTE: second derivatives contains a additional factor of -4*pi from Poisson equation
    !---------------------------------------------------------------------------------------------
    !                            --- \nabla^2 V = -4*\pi \rho ---
    !---------------------------------------------------------------------------------------------
    ! axial derivatives
    !---------------------------------------------------------------------------------------------
    coeke(:,3,3) = -coeke(:,1,1)/(hc*hc*fpi)                               ! -d^2/dc^2/4pi stencil
    coeke(:,2,2) = -coeke(:,1,1)/(hb*hb*fpi)                               ! -d^2/db^2/4pi stencil
    coeke(:,1,1) = -coeke(:,1,1)/(ha*ha*fpi)                               ! -d^2/da^2/4pi stencil
    !---------------------------------------------------------------------------------------------
    ! cross derivatives
    !---------------------------------------------------------------------------------------------
    coeke(:,2,3) = -coeke(:,1,2)/(hb*hc*fpi)                               ! -d^2/dbdc/4pi stencil
    coeke(:,1,3) = -coeke(:,1,2)/(ha*hc*fpi)                               ! -d^2/dadc/4pi stencil
    coeke(:,1,2) = -coeke(:,1,2)/(ha*hb*fpi)                               ! -d^2/dadb/4pi stencil
    !---------------------------------------------------------------------------------------------
    ! Jacobian for the non-orthogonal cells 
    !---------------------------------------------------------------------------------------------
    !                          --- J = transpose(ainv).(diag(a)) ---
    !---------------------------------------------------------------------------------------------
    Jim(:,1) = ainv(1,:)*a(1)                                                   ! i={xyz}, m={abc}
    Jim(:,2) = ainv(2,:)*a(2) 
    Jim(:,3) = ainv(3,:)*a(3)
    !
    !---------------------------------------------------------------------------------------------
    ! --- weigh coeke with the Jacobian ---
    !---------------------------------------------------------------------------------------------
    ! axial derivatives
    !---------------------------------------------------------------------------------------------
    coeke(:,3,3) = (Jim(1,3)**2+Jim(2,3)**2+Jim(3,3)**2)*coeke(:,3,3)
    coeke(:,2,2) = (Jim(1,2)**2+Jim(2,2)**2+Jim(3,2)**2)*coeke(:,2,2)
    coeke(:,1,1) = (Jim(1,1)**2+Jim(2,1)**2+Jim(3,1)**2)*coeke(:,1,1)
    !---------------------------------------------------------------------------------------------
    !
    !---------------------------------------------------------------------------------------------
    ! cross derivatives (needed for non-othogonal grids in the second derivatives)
    !---------------------------------------------------------------------------------------------
    coeke(:,2,3) = 2.0_DP*(Jim(1,2)*Jim(1,3)+Jim(2,2)*Jim(2,3)+Jim(3,2)*Jim(3,3))*coeke(:,2,3)
    coeke(:,1,3) = 2.0_DP*(Jim(1,1)*Jim(1,3)+Jim(2,1)*Jim(2,3)+Jim(3,1)*Jim(3,3))*coeke(:,1,3)
    coeke(:,1,2) = 2.0_DP*(Jim(1,1)*Jim(1,2)+Jim(2,1)*Jim(2,2)+Jim(3,1)*Jim(3,2))*coeke(:,1,2)
    coeke(:,3,2) = coeke(:,2,3)                                                ! symmetry of coeke
    coeke(:,2,1) = coeke(:,1,2)                                                ! symmetry of coeke
    coeke(:,3,1) = coeke(:,1,3)                                                ! symmetry of coeke
    !---------------------------------------------------------------------------------------------
    !CALL coeff_opt()      ! calculate coemicf (preconditioner)
coemicf = 0.d0 ! MCA/HK: dirty hack for std CG...
#ifdef __CUDA
    coeke_d = coeke
    coemicf_d = coemicf
    coe_1st_derv_d = coe_1st_derv
#endif
    !---------------------------------------------------------------------------------------------
    ! a samall check on the shape of user defined cell (if any)
    !---------------------------------------------------------------------------------------------
    IF ((ibrav.EQ.0).AND.(nfi.EQ.1)) THEN
      WRITE(stdout,*) 'EXX info: If you are using an orthogonal cell without its cell vectors&
          & aligned to the xyz directions, the EXX calculation may be twice more expensive.'
    END IF
    !---------------------------------------------------------------------------------------------
    !
    CALL stop_clock('exx_gs_setup')
    !
    !---------------------------------------------------------------------------------------------
    !
    CALL start_clock('exx_pairs')
    !
    !---------------------------------------------------------------------------------------------
    ndim=MAX(nproc_image, nbsp)
    ALLOCATE (wannierc(3,ndim)); wannierc=0.0_DP 
    ALLOCATE (wannierc_tmp(3,nbsp)); wannierc_tmp=0.0_DP 
    !---------------------------------------------------------------------------------------------
    ! Adjust Cartesian coordinates of wannier centres according to periodic boundary conditions...
    ! N.B.: PBC are imposed here in the range [0,1)... 
    !---------------------------------------------------------------------------------------------
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
    !---------------------------------------------------------------------------------------------
    ! 
    !---------------------------------------------------------------------------------------------
    ! make copy of wannier centres when number of processors > number of bands
    !---------------------------------------------------------------------------------------------
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
    !---------------------------------------------------------------------------------------------
    !
    !---------------------------------------------------------------------------------------------
    ! overlap is the unique neighbor list that for each band or processor image (ndim)
    ! num_recv is the number of unique neighbors for each band or processor image (ndim)
    !---------------------------------------------------------------------------------------------
    ALLOCATE (obtl_recv(neigh/2, ndim)); obtl_recv=0
    ALLOCATE (obtl_send(neigh, ndim));   obtl_send=0
    ALLOCATE (num_recv(ndim)); num_recv=0
    ALLOCATE (num_send(ndim)); num_send=0
    !---------------------------------------------------------------------------------------------
    ! generate the unique neighbor list
    !---------------------------------------------------------------------------------------------
    CALL exx_index_pair(wannierc_tmp, obtl_recv, num_recv, nj_max, ndim)
    !---------------------------------------------------------------------------------------------
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
    !---------------------------------------------------------------------------------------------
    ! DO itr = 1, ndim
    !   DO jtr = 1, neigh/2
    !     WRITE(*,*) itr, "   ", obtl_recv(jtr, itr) 
    !   END DO
    ! END DO
    !---------------------------------------------------------------------------------------------
    ! DO itr = 1, ndim
    !   DO jtr = 1, neigh
    !     WRITE(*,*) itr, "   ", obtl_send(jtr, itr) 
    !   END DO
    ! END DO
    !---------------------------------------------------------------------------------------------
    !
    ! Allocate variables to store potentials for 3 steps ....
    !
    IF (n_exx.EQ.0) THEN
      !
      ! the following variables are used in the extrapolation of exx potentials
      ALLOCATE( selfv  ( n_s_ps, psgsn, my_nbspx), stat=ierr ); selfv=0.0_DP
      ALLOCATE( selfrho( n_s_ps, psgsn, my_nbspx), stat=ierr ); selfrho=0.0_DP
      ALLOCATE( pairv  ( n_p_ps, psgsn, neigh, my_nbspx), stat=ierr ); pairv=0.0_DP
      ALLOCATE( pairrho( n_p_ps, psgsn, neigh, my_nbspx), stat=ierr ); pairrho=0.0_DP
      ALLOCATE( pair_dist( psgsn, neigh, my_nbspx), stat=ierr ); pair_dist=0.0_DP
      ALLOCATE( pair_label( neigh, my_nbspx ), stat=ierr ); pair_label=0.0_DP
      ALLOCATE( pair_step( neigh, my_nbspx ), stat=ierr ); pair_step=0.0_DP
      ALLOCATE( pair_status( neigh, my_nbspx ), stat=ierr ); pair_status=0.0_DP
      ALLOCATE( prev_obtl_recv( neigh/2, ndim ), stat=ierr ); prev_obtl_recv=0
      !
    END IF
    !
    prev_obtl_recv = obtl_recv
    !
    !---------------------------------------------------------------------------------------------
    !
    ! update exx step ...
    !
    n_exx = n_exx + 1
    !
    !=========================================================================
    !
    ! obtain orbitals on each local processor, stored in psi
    !
    ALLOCATE ( psi   ( nnrtot, my_nbsp(me) ) ); psi=0.0_DP
    ALLOCATE ( potpsi( nnrtot, my_nbsp(me) ) ); potpsi=0.0_DP
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
#ifdef __CUDA
    IF(.not.ALLOCATED(psi_d ))    ALLOCATE(psi_d(nnrtot, my_nbsp(me)) )
    psi_d = psi
#endif
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
    ! obtain psi in sphere (psime) for neighbors
    !
    call start_clock('exx_big_alloc')
    ! TODO check : if neigh changes
    if (.not.allocated(psime_pair_send)) then
      ALLOCATE ( psime_pair_send(n_p_me, neigh, my_nbsp(me)) ); psime_pair_send=0.0_DP
    end if
    if (.not.allocated(psime_pair_recv)) then
      ALLOCATE ( psime_pair_recv(n_p_me, neigh, my_nbsp(me)) ); psime_pair_recv=0.0_DP
    end if
#ifdef __CUDA
    if (.not.allocated(psime_pair_send_d)) then
      ALLOCATE ( psime_pair_send_d(n_p_me, neigh, my_nbsp(me)) ); psime_pair_send_d=0.0_DP
    end if
    if (.not.allocated(psime_pair_recv_d)) then
      ALLOCATE ( psime_pair_recv_d(n_p_me, neigh, my_nbsp(me)) ); psime_pair_recv_d=0.0_DP
    end if
#endif
    call stop_clock('exx_big_alloc')
    !
    ! initialize totalenergy and derivatives
    !
    totalenergy = 0.0_DP
    total_exx_derv(:,:) = 0.0_DP
    irecv_count = 0
    isend_count = 0
    !
    ! my_var is the maximum of nbsp or nproc_image
    my_var = MAX(nproc_image, nbsp)
    !
    !==========================================================================
    !
    CALL start_clock('send_psi')
    ! 
    !==========================================================================
    !
    ! we should use my_nbspx (maxval(my_nbsp(me))) here
    !
    !MCA: TODO try to communicate everything via GPU
    !
    DO iobtl = 1, my_nbspx
      ! 
      !========================================================================
      !
      gindex_of_iobtl = index_my_nbsp(iobtl, me)
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
          CALL MPI_IRECV( psime_pair_recv(1, itr, iobtl), n_p_me, MPI_DOUBLE_PRECISION, rk_of_obtl_tbs, &
                          obtl_tbs*ndim+gindex_of_iobtl, intra_image_comm, irecvreq(irecv_count), ierr)
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
#ifdef __CUDA
          CALL getpsicb( nrg, p_me_r, psi_d(1,iobtl), psime_pair_send_d(1, itr, iobtl), tran)
          psime_pair_send(:,itr,iobtl) = psime_pair_send_d(:, itr, iobtl)
#else
          CALL getpsicb( nrg, p_me_r, psi(1,iobtl), psime_pair_send(1, itr, iobtl), tran)
#endif
          !
          isend_count = isend_count + 1
          CALL MPI_ISEND( psime_pair_send(1, itr, iobtl), n_p_me, MPI_DOUBLE_PRECISION, rk_of_obtl_trcv, &
                          gindex_of_iobtl*ndim+obtl_trcv, intra_image_comm, isendreq(isend_count), ierr)
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

    DO itr = 1, irecv_count
      CALL MPI_WAIT(irecvreq(itr), istatus, ierr)
    END DO
    !
    CALL stop_clock('send_psi_wait')
    !
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
    !========================================================================
    !
    CALL start_clock('getpairv')
    !
    !========================================================================
    !
    ! Do some allocations
#ifdef __CUDA
    IF(.not.ALLOCATED(potpsi_d )) ALLOCATE(potpsi_d(nnrtot, my_nbsp(me)) )
    potpsi_d = 0._DP
#endif
    IF(.not.ALLOCATED(psime )) ALLOCATE( psime(max(n_p_me,n_s_me)) ); 
    IF(.not.ALLOCATED(rhome )) ALLOCATE( rhome(max(n_p_me,n_s_me)) );
    IF(.not.ALLOCATED(rhops )) ALLOCATE( rhops(max(n_p_ps,n_s_ps)) ); 
    IF(.not.ALLOCATED(potme )) ALLOCATE( potme(max(n_p_me,n_s_me)) ); 
    !
    DO iobtl = 1, my_nbspx
      ! 
      middle(:)=0.0_DP
      !
      gindex_of_iobtl = index_my_nbsp(iobtl, me)
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
            guess_status = 0        ! no guess
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
            !
            !==================================================================================
            !
            call start_clock('exx_grid_trans')
            CALL getmiddlewc(wannierc(1,gindex_of_iobtl),wannierc(1,my_var2), h, ainv, middle )
            ! get pair distance
            !CALL get_pair_dist(wannierc(1,gindex_of_iobtl),wannierc(1,my_var2),d_pair)
            !
            ! calculate translation vector from the center of the box
            CALL getsftv( nr1s, nr2s, nr3s, h, ainv, middle, tran)
            call stop_clock('exx_grid_trans')
            !      
            ! get the localized psi around the mid point of two wannier centers
            ! note: the psime is centered at the center of the box
            ! (using the translation vector "tran" from middle of wfc to the center of box)
            call start_clock('exx_psicb')
#ifdef __CUDA
            associate(psi=>psi_d)
#endif
            CALL getpsicb( nrg, p_me_r, psi(1,iobtl), psime(1), tran)
#ifdef __CUDA
            end associate
#endif

            call stop_clock('exx_psicb')
            ! 
            ! to be continue
            !
            ! the localized density rhome 
            call start_clock('exx_getrhol')
#ifdef __CUDA
            psime_pair_recv_d(:, j, iobtl)=psime_pair_recv(:, j, iobtl)
            CALL getrhol(p_me_r, p_ps_r, psime(1), psime_pair_recv_d(1, j, iobtl), rhome, rhops, inv_omega)
#else
            CALL getrhol(p_me_r, p_ps_r, psime(1), psime_pair_recv(1, j, iobtl), rhome, rhops, inv_omega)
#endif
            call stop_clock('exx_getrhol')
            !
            ! calculate the exx potential from the pair density by solving Poisson
            !
            !--------------------------------------------------------------------------------------
            call start_clock('exx_vofr')
            CALL getvofr( p_me_r, p_ps_r, n_p_me, n_p_ps, hcub, rhops, potme, pair_status(pos, iobtl), psgsn, &
                          pairrho(:,:,pos,iobtl), pairv(:,:,pos,iobtl), cgstep)
!            CALL getvofr( p_me_r, p_ps_r, max(n_p_me,n_s_me), max(n_p_ps,n_s_ps), hcub, rhops, potme, pair_status(pos, iobtl), psgsn, &
!                          pairrho(:,:,pos,iobtl), pairv(:,:,pos,iobtl), cgstep)
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
#ifdef __CUDA
            CALL updateforce_loc(nrg, p_me_r, potpsi_d(:,iobtl), potme, psime, psime_pair_recv_d(1,j,iobtl),tran)
            psime_pair_recv(:,j,iobtl) = psime_pair_recv_d(:,j,iobtl)
#else
            CALL updateforce_loc(nrg, p_me_r, potpsi(:,iobtl), potme, psime, psime_pair_recv(1,j,iobtl),tran)
#endif
            call stop_clock('exx_force_loc')
            !
            call start_clock('exx_penergy')
            !
            CALL vvprod(p_me_r, rhome, potme, paire(j))    ! dot product of the rho and potme  

            call stop_clock('exx_penergy')
            paire(j) = paire(j) * 0.5_DP* hcub             ! volume element hcub and trapezoidal rule prefactor 0.5_DP are included
            totalenergy = totalenergy + 2.0_DP*paire(j)    ! the factor of two comes from the identity of ij and ji pair
            !--------------------------------------------------------------------------------------
            !
            ! write (my_unit, *) "pair = ", gindex_of_iobtl, "+", my_var2, "pos = ", pos, &
            !                    "cgstep = ", cgstep, "guess = ", pair_status(pos, iobtl) ! debug
            !
            IF (.NOT. (isotropic .AND. (ibrav.EQ.1) )) THEN
              call start_clock('exx_stress')
              !
              ! EXX cell derivative (note: exxalfa is included in vofrho.f90 when calculate stress)
              ! CALL start_clock('exx_cell_derv')
              !
              CALL exx_energy_cell_derivative(p_me_r, p_ps_r, tran, rhome, potme, &
                                              ha_proj, hb_proj, hc_proj, Jim, pair_dexx_dhab(:,:,j))
              !
              ! volume element hcub and trapezoidal rule prefactor 0.5_DP are included
              !
              pair_dexx_dhab(:,:,j) = pair_dexx_dhab(:,:,j)*0.5_DP*hcub
              !
              ! accumulate the derivative from different pair terms
              !
              total_exx_derv(:,:) = total_exx_derv(:,:) + 2.0_DP*pair_dexx_dhab(:,:,j)
              !
              ! CALL stop_clock('exx_cell_derv')
              !
              ! if isotropic => calculate the stress tensor in vofrho.f90
              !
              call stop_clock('exx_stress')
            END IF
            
          END IF
          !
        END DO !for j
        !
      END IF !gindex_of_iobtl <= nbsp
      !
    END DO
    !
    !===============================================================================
    !
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
    irecv_count = 0
    isend_count = 0
    !
    CALL start_clock('send_v')
    !
    !========================================================================
    !
    DO iobtl = 1, my_nbspx
      !
      !========================================================================
      !
      gindex_of_iobtl = index_my_nbsp(iobtl, me)
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
          CALL MPI_ISEND( psime_pair_recv(1,itr,iobtl), n_p_me, MPI_DOUBLE_PRECISION, rk_of_obtl_trcv, &
                          gindex_of_iobtl*ndim+obtl_trcv, intra_image_comm, isendreq(isend_count), ierr)
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
          CALL MPI_IRECV( psime_pair_send(1, itr, iobtl), n_p_me, MPI_DOUBLE_PRECISION, rk_of_obtl_tbs, &
                          obtl_tbs*ndim+gindex_of_iobtl, intra_image_comm, irecvreq(irecv_count), ierr)
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
      ! CALL moments(psi(:,:,:,iobtl), nnrtot, 1)
      ! CALL moments(psi(:,:,:,iobtl), nnrtot, 2)
      ! CALL moments(psi(:,:,:,iobtl), nnrtot, 3)
      !
      IF (iobtl.LE.my_nbsp(me)) THEN ! skip when the loop of my_nbspx goes outside of scope
        !
        IF (me.GT.(nbsp*(sc_fac-1))) THEN ! compatible with more processors than nbsp
          !
          gindex_of_iobtl =  index_my_nbsp(iobtl, me)
          !
          ! calculate translation vector from the center of the box
          CALL getsftv(nr1s, nr2s, nr3s, h, ainv, wannierc(1, gindex_of_iobtl), tran)
          !
          ! get the localized psi around the wannier centers
          ! note: the psime is centered at the center of the box
          ! (using the translation vector "tran" from the wfc to the center of box)
          !
#ifdef __CUDA
          associate( psi=>psi_d )
#endif
          CALL getpsicb( nrg, s_me_r, psi(1,iobtl), psime(1), tran)
#ifdef __CUDA
          end associate
#endif
          ! get the localized density rhome  
          CALL getrhol(s_me_r, s_ps_r, psime(1), psime(1), rhome, rhops, inv_omega)
          ! 
          ! calculate the exx potential from the pair density by solving Poisson
          !
          !--------------------------------------------------------------------------------------
          CALL start_clock('getvofr')
          !
          CALL getvofr( s_me_r, s_ps_r, n_s_me, n_s_ps, hcub, rhops, potme, n_exx-1, psgsn, &
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
#ifdef __CUDA
          associate( potpsi=>potpsi_d )
#endif
          CALL updateforce_slf(nrg, s_me_r, potpsi(1,iobtl), potme, psime, tran)
          CALL vvprod(s_me_r, rhome, potme, selfe)    ! dot product of the rho and potme 
#ifdef __CUDA
          end associate
#endif
          selfe = selfe * 0.5_DP* hcub                ! volume element hcub and trapezoidal rule prefactor 0.5_DP are included
          totalenergy = totalenergy + selfe    ! the factor of two comes from the identity of ij and ji pair
          !--------------------------------------------------------------------------------------
          !
          ! write (*,*) "selfenergy = ", selfe ! debug
          !
          !IF(me .GT. nbsp) THEN
          !  potpsi(:,iobtl) = 0.0_DP
          !END IF
          !
          IF (.NOT. (isotropic .AND. (ibrav.EQ.1))) THEN
            !
            !  EXX cell derivative (note: need to include exxalfa later)
            ! CALL start_clock('exx_cell_derv')
            !
            CALL exx_energy_cell_derivative(s_me_r, s_ps_r, tran, rhome, potme, &
                                            ha_proj, hb_proj, hc_proj, Jim, self_dexx_dhab(:,:))
            !
            ! volume element hcub and trapezoidal rule prefactor 0.5_DP are included
            !
            self_dexx_dhab(:,:) = self_dexx_dhab(:,:)*0.5_DP*hcub
            !
            ! combine derivative with pair terms
            !
            total_exx_derv(:,:) = total_exx_derv(:,:) + self_dexx_dhab(:,:)
            !
            ! CALL stop_clock('exx_cell_derv')
            !
            ! if isotropic => calculate the stress tensor in vofrho.f90
            !
          END IF
          !
        END IF ! me
        !
      END IF !iobtl 
      !
    END DO
    !
    ! clean up
    IF (ALLOCATED(psime))           DEALLOCATE(psime)
    IF (ALLOCATED(rhome))           DEALLOCATE(rhome)
    IF (ALLOCATED(rhops))           DEALLOCATE(rhops)
    IF (ALLOCATED(potme))           DEALLOCATE(potme)
    ! clean up
    !=========================================================================
    !
    CALL stop_clock('getselfv')
    !
    !========================================================================
    !
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
    !
    !========================================================================
    !
    CALL start_clock('force_rec')
    !
    DO iobtl = 1, my_nbspx
      !
      gindex_of_iobtl =  index_my_nbsp(iobtl, me)
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
          
          ! upadate potpsi PBE0 
          !
#ifdef __CUDA
          psime_pair_send_d(:,itr,iobtl) = psime_pair_send(:,itr,iobtl)
          CALL updateforce_rec(nrg, p_me_r, potpsi_d(1,iobtl), psime_pair_send_d(1,itr,iobtl), tran)
#else
          CALL updateforce_rec(nrg, p_me_r, potpsi(1,iobtl), psime_pair_send(1,itr,iobtl), tran)
#endif
          !

          !
        END IF
        !
      END DO
      !
    END DO ! iobtl
    !
#ifdef __CUDA
    potpsi=potpsi_d
#endif
    !
    !MCA: end of GPU stuff
    !
    CALL stop_clock('force_rec')
    !
    !============================================================================
    !THE MOST OUTER LOOP FOR PAIR POTENTIAL ENDS HERE
    !============================================================================
    !
    !============================================================================
    CALL start_clock('totalenergy')
    !
    totalenergyg=0.0_DP ! mpi reduction variable initialization
    exx=0.0_DP          ! exx energy (used to handle the open/closed shell energy)
    ! 
#ifdef __MPI
    ! collect the totalenergy of each mpi task to totalenergyg
    CALL MPI_ALLREDUCE(totalenergy, totalenergyg, 1, MPI_DOUBLE_PRECISION, &
        &                        MPI_SUM, intra_image_comm, ierr)
    !
#endif
    exx = totalenergyg
    IF (nspin .EQ. 1) exx = exx + totalenergyg ! if closed shell double the totalenergy
#include "debug_patch0.f90"
    !
    !WRITE(stdout, '("EXX Energy",2F30.14," step",I7)')exx,totalenergyg*2.0_DP, nfi
    !
    CALL stop_clock('totalenergy')
    !
    ! CALL start_clock('exx_cell_derv')
    !
    total_exx_derv_g(:,:) = 0.0_DP ! mpi reduction variable initialization
    !
#ifdef __MPI
    IF (.NOT. (isotropic .AND. (ibrav.EQ.1))) THEN
      ! collect the total_exx_derv of each mpi task to total_exx_derv_g
      CALL MPI_ALLREDUCE(total_exx_derv(:,:), total_exx_derv_g(:,:), 9, &
          MPI_DOUBLE_PRECISION, MPI_SUM, intra_image_comm, ierr)
      !
    END IF
#endif
    !
    ! for closed shell case inclued spin factor of 2
    dexx_dh(:,:) = total_exx_derv_g(:,:)
    IF (nspin .EQ. 1) dexx_dh(:,:) = dexx_dh(:,:) + total_exx_derv_g(:,:)
    !
    ! CALL stop_clock('exx_cell_derv')
    !
    ! Local to global distribution of EXX potential
    ! potpsi (local) --> exx_potential (global)
    !
    CALL start_clock('vl2vg')
    exx_potential=0.0_DP
    !
    IF (nproc_image .LE. nbsp) THEN 
      !
      CALL redistwfr ( exx_potential, potpsi, my_nxyz, my_nbsp, intra_image_comm, -1 )
#include "debug_patch.f90"
      !
    ELSE
      !
      !-----------Zhaofeng's potpsi (local) to exx_potential (global) -----------
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
#ifdef __MPI
      !
      CALL mp_barrier( intra_image_comm )
      CALL MPI_ALLTOALLV(potpsi, recvcount,rdispls,MPI_DOUBLE_PRECISION, &
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
#ifdef __MPI
      !
      DO ir=1,nproc_image/nogrp
        !CALL mp_barrier( dffts%ogrp_comm )
        CALL mp_barrier( fftx_tgcomm(dffts) )
        CALL MPI_ALLTOALLV(exx_tmp3(1,ir), sendcount1, sdispls1, MPI_DOUBLE_PRECISION, &
            &         exx_potential(1,ir),recvcount1, rdispls1, MPI_DOUBLE_PRECISION, &
            &         fftx_tgcomm(dffts), ierr)
        !    &         dffts%ogrp_comm, ierr)
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
      !-----------Zhaofeng's potpsi (local) to exx_potential (global) -----------
      !
    END IF ! vl2vg
    !
    CALL stop_clock('vl2vg')
    !
!     IF ( iverbosity > 1 ) THEN
!       !
!      DO ir=1,nbsp/nogrp
!        DO ia=1,dffts%nnr
!          ifile = 1+10000*ir+me
!          write(chr_me,"(I4.4,'-',I4.4)") me, ir
!          open(unit=ifile,file='exx_potential.'//trim(chr_me))
!          write (ifile,*) exx_potential(ia,ir)
!          close(ifile)
!        END DO
!      END DO
!       !
!     END IF

     write(*,*) "exiting exx_gs" ! debug
    !
    !==============================================================================
    IF (ALLOCATED(potpsi))           DEALLOCATE(potpsi)
    IF (ALLOCATED(psi))             DEALLOCATE(psi)
#ifdef __CUDA
    IF (ALLOCATED(psi_d))             DEALLOCATE(psi_d)
#endif
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
END SUBROUTINE exx_gs
!====================================================================================

!==============================================================================
SUBROUTINE getsftv(nr1s, nr2s, nr3s, h, ainv, wc, tran)
    !
    USE kinds, ONLY : DP
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
SUBROUTINE vvprod(me_r, v1, v2, prod)
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
#ifdef __CUDA
    attributes(device) :: v1,v2
#endif
    !----------------------------------------------------------------
    !
    prodp=0.0D0
    !
    ! WRITE(*,*) "vvprod"
    !
#ifdef __CUDA
    !$cuf kernel do (3)
#else
    !$omp parallel do private(i,j,k) reduction(+:prodp)
#endif    
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
#ifndef __CUDA
    !$omp end parallel do 
#endif
    !
    prod = prodp
    !
    RETURN
END SUBROUTINE vvprod
!==============================================================================

!==============================================================================
SUBROUTINE getpsicb(nrg,nrl,psig,psil,tran)
    !
    USE kinds, ONLY  : DP
    USE dummy_exx, ONLY  : l2gcb
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
#ifdef __CUDA
    attributes(device) :: psil, psig
#endif
    !
    ti = tran(1); tj = tran(2); tk = tran(3)
#ifdef __CUDA
    !$cuf kernel do (3)
#else
    !$omp parallel do private(i,j,k,gi,gj,gk)
#endif
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
#ifndef __CUDA
    !$omp end parallel do 
#endif
    !----------------------------------------------------------
    !
    RETURN
END SUBROUTINE getpsicb
!==============================================================================
SUBROUTINE getpsicb_cpu(nrg,nrl,psig,psil,tran)
    !
    USE kinds, ONLY  : DP
    USE dummy_exx, ONLY  : l2gcb
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
END SUBROUTINE getpsicb_cpu

!====================================================================================
SUBROUTINE moments(wavefunc, ntot, control)
    !
    USE kinds,                   ONLY  : DP
    USE cell_base,               ONLY  : omega, h
    USE fft_base,                ONLY  : dfftp
    USE io_global,               ONLY  : stdout
    !
    IMPLICIT NONE
    INTEGER                 ::  i, j, k, offset, ntot, nr1s, nr2s, nr3s, control
    REAL(DP)                ::  wavefunc(ntot), pos_ratio(3), pos(3)
    REAL(DP)                ::  rho_r
    REAL(DP)                ::  e_x, e_y, e_z
    REAL(DP)                ::  e_xx, e_yy, e_zz, e_xy, e_yz, e_xz
    REAL(DP)                ::  e_xxx, e_xxy, e_xxz, e_xyy, e_xyz, e_xzz, e_yyy, e_yyz, e_yzz, e_zzz
    REAL(DP)                ::  x, y, z, nelec, hcub, norm
    REAL(DP), ALLOCATABLE   ::  xcoeff(:), ycoeff(:), zcoeff(:)
    ! 
    hcub = omega/DBLE(ntot)
    nr1s=dfftp%nr1
    nr2s=dfftp%nr2
    nr3s=dfftp%nr3 
    norm=1.0_DP/DBLE(nr1s*nr2s*nr3s)
    !
    ALLOCATE(xcoeff(nr1s+1))
    ALLOCATE(ycoeff(nr2s+1))
    ALLOCATE(zcoeff(nr3s+1))
    !
    IF (control .EQ. 1) THEN
      CALL get_trapezoid_coeff(xcoeff, nr1s+1)
      CALL get_trapezoid_coeff(ycoeff, nr2s+1)
      CALL get_trapezoid_coeff(zcoeff, nr3s+1)
    ELSE IF (control .EQ. 2) THEN
      CALL get_simpson_coeff(xcoeff, nr1s+1)
      CALL get_simpson_coeff(ycoeff, nr2s+1)
      CALL get_simpson_coeff(zcoeff, nr3s+1)
    ELSE IF (control .EQ. 3) THEN
      CALL get_simpson38_coeff(xcoeff, nr1s+1)
      CALL get_simpson38_coeff(ycoeff, nr2s+1)
      CALL get_simpson38_coeff(zcoeff, nr3s+1)
    END IF
    !
    !--------------------
    e_x=0.0_DP
    e_y=0.0_DP
    e_z=0.0_DP
    !--------------------
    e_xx=0.0_DP
    e_yy=0.0_DP
    e_zz=0.0_DP
    e_xy=0.0_DP
    e_yz=0.0_DP
    e_xz=0.0_DP
    !--------------------
    e_xxx=0.0_DP
    e_xxy=0.0_DP
    e_xxz=0.0_DP
    e_xyy=0.0_DP
    e_xyz=0.0_DP
    e_xzz=0.0_DP
    e_yyy=0.0_DP
    e_yyz=0.0_DP
    e_yzz=0.0_DP
    e_zzz=0.0_DP
    !--------------------
    !
    DO i = 1, nr1s+1
      pos_ratio(1) = DBLE(i-1)/DBLE(nr1s) 
      DO j = 1, nr2s+1
        pos_ratio(2) = DBLE(j-1)/DBLE(nr2s) 
        DO k = 1, nr3s+1
          pos_ratio(3) = DBLE(k-1)/DBLE(nr3s) 
          !------------------------------------------------------------------------------------
          pos(1)=h(1,1)*pos_ratio(1)+h(1,2)*pos_ratio(2)+h(1,3)*pos_ratio(3)
          pos(2)=h(2,1)*pos_ratio(1)+h(2,2)*pos_ratio(2)+h(2,3)*pos_ratio(3)
          pos(3)=h(3,1)*pos_ratio(1)+h(3,2)*pos_ratio(2)+h(3,3)*pos_ratio(3)
          !------------------------------------------------------------------------------------
          offset = MOD(i,nr1s) + MOD(j-1,nr2s)*nr1s + MOD(k-1,nr3s)*nr1s*nr2s
          rho_r  = wavefunc(offset) * wavefunc(offset)
          !------------------------------------------------------------------------------------
          !
          !------------------------------------------------------------------------------------
          nelec  = nelec + rho_r                            * xcoeff(i) * ycoeff(j) * zcoeff(k)  
          !------------------------------------------------------------------------------------
          e_x    = e_x   + rho_r * pos(1)                   * xcoeff(i) * ycoeff(j) * zcoeff(k)
          e_y    = e_y   + rho_r * pos(2)                   * xcoeff(i) * ycoeff(j) * zcoeff(k)
          e_z    = e_z   + rho_r * pos(3)                   * xcoeff(i) * ycoeff(j) * zcoeff(k)  
          !------------------------------------------------------------------------------------
          e_xx   = e_xx  + rho_r * pos(1) * pos(1)          * xcoeff(i) * ycoeff(j) * zcoeff(k)  
          e_yy   = e_yy  + rho_r * pos(2) * pos(2)          * xcoeff(i) * ycoeff(j) * zcoeff(k)  
          e_zz   = e_zz  + rho_r * pos(3) * pos(3)          * xcoeff(i) * ycoeff(j) * zcoeff(k)  
          e_xy   = e_xy  + rho_r * pos(1) * pos(2)          * xcoeff(i) * ycoeff(j) * zcoeff(k)  
          e_yz   = e_yz  + rho_r * pos(2) * pos(3)          * xcoeff(i) * ycoeff(j) * zcoeff(k)  
          e_xz   = e_xz  + rho_r * pos(2) * pos(3)          * xcoeff(i) * ycoeff(j) * zcoeff(k)  
          !------------------------------------------------------------------------------------
          e_xxx  = e_xxx + rho_r * pos(1) * pos(1) * pos(1) * xcoeff(i) * ycoeff(j) * zcoeff(k)  
          e_xxy  = e_xxy + rho_r * pos(1) * pos(1) * pos(2) * xcoeff(i) * ycoeff(j) * zcoeff(k)  
          e_xxz  = e_xxz + rho_r * pos(1) * pos(1) * pos(3) * xcoeff(i) * ycoeff(j) * zcoeff(k)  
          e_xyy  = e_xyy + rho_r * pos(1) * pos(2) * pos(2) * xcoeff(i) * ycoeff(j) * zcoeff(k)  
          e_xyz  = e_xyz + rho_r * pos(1) * pos(2) * pos(3) * xcoeff(i) * ycoeff(j) * zcoeff(k)  
          e_xzz  = e_xzz + rho_r * pos(1) * pos(3) * pos(3) * xcoeff(i) * ycoeff(j) * zcoeff(k)  
          e_yyy  = e_yyy + rho_r * pos(2) * pos(2) * pos(2) * xcoeff(i) * ycoeff(j) * zcoeff(k)  
          e_yyz  = e_yyz + rho_r * pos(2) * pos(2) * pos(3) * xcoeff(i) * ycoeff(j) * zcoeff(k)  
          e_yzz  = e_yzz + rho_r * pos(2) * pos(3) * pos(3) * xcoeff(i) * ycoeff(j) * zcoeff(k)  
          e_zzz  = e_zzz + rho_r * pos(3) * pos(3) * pos(3) * xcoeff(i) * ycoeff(j) * zcoeff(k)  
          !------------------------------------------------------------------------------------
        END DO
      END DO
    END DO
    !
    !--------------------
    nelec  = norm*nelec
    !--------------------
    e_x    = norm*e_x
    e_y    = norm*e_y
    e_z    = norm*e_z
    !--------------------
    e_xx   = norm*e_xx
    e_yy   = norm*e_yy
    e_zz   = norm*e_zz
    e_xy   = norm*e_xy
    e_yz   = norm*e_yz
    e_xz   = norm*e_xz
    !--------------------
    e_xxx  = norm*e_xxx
    e_xxy  = norm*e_xxy
    e_xxz  = norm*e_xxz
    e_xyy  = norm*e_xyy
    e_xyz  = norm*e_xyz
    e_xzz  = norm*e_xzz
    e_yyy  = norm*e_yyy
    e_yyz  = norm*e_yyz
    e_yzz  = norm*e_yzz
    e_zzz  = norm*e_zzz
    !--------------------
    !
    !-----------------------------
    WRITE(stdout,*)  "  x:", e_x
    WRITE(stdout,*)  "  y:", e_y
    WRITE(stdout,*)  "  z:", e_z
    !-----------------------------
    WRITE(stdout,*)  " xx:", e_xx
    WRITE(stdout,*)  " yy:", e_yy
    WRITE(stdout,*)  " zz:", e_zz
    WRITE(stdout,*)  " xy:", e_xy
    WRITE(stdout,*)  " yz:", e_yz
    WRITE(stdout,*)  " xz:", e_xz
    !-----------------------------
    WRITE(stdout,*)  "xxx:", e_xxx
    WRITE(stdout,*)  "xxy:", e_xxy
    WRITE(stdout,*)  "xxz:", e_xxz
    WRITE(stdout,*)  "xyy:", e_xyy
    WRITE(stdout,*)  "xyz:", e_xyz
    WRITE(stdout,*)  "xzz:", e_xzz
    WRITE(stdout,*)  "yyy:", e_yyy
    WRITE(stdout,*)  "yyz:", e_yyz
    WRITE(stdout,*)  "yzz:", e_yzz
    WRITE(stdout,*)  "zzz:", e_zzz
    !-----------------------------
    !
    DEALLOCATE(xcoeff)
    DEALLOCATE(ycoeff)
    DEALLOCATE(zcoeff)
    !
END SUBROUTINE moments
!====================================================================================

!====================================================================================
SUBROUTINE get_trapezoid_coeff(coeff, n)
    !
    USE kinds, ONLY : DP
    !
    IMPLICIT NONE
    !
    INTEGER    :: n, i
    REAL(DP)   :: coeff(n)
    !
    IF (n .GE. 11) THEN
      coeff(1)    = 1.0_DP/2 
      coeff(n)    = 1.0_DP/2 
      !
      DO i = 2, n-1
        coeff(i)  = 1.0_DP
      END DO
    ELSE
      CALL errore( 'exx_gs', 'trapezoid rule do not apply for n < 11', 1 )
    END IF
    !
END SUBROUTINE get_trapezoid_coeff
!===================================================================================

!====================================================================================
SUBROUTINE get_simpson_coeff(coeff, n)
    !
    USE kinds, ONLY : DP
    !
    IMPLICIT NONE
    !
    INTEGER    :: n, i, ep
    REAL(DP)   :: coeff(n)
    !
    ep =2*((n-1)/2)+1
    !
    IF (n .GE. 21) THEN
      coeff(1)    = 1.0_DP/3
      coeff(ep)   = 1.0_DP/3
      !
      DO i = 2, ep-1
        IF (MOD(i,2) .EQ. 0) THEN
          coeff(i) = 4.0_DP/3
        ELSE
          coeff(i) = 2.0_DP/3
        END IF
      END DO
      !
      DO i = ep+1, n
        coeff(i) = 1.0_DP
      END DO
    ELSE
      CALL errore( 'exx_gs', 'simpson rule do not apply for n < 21', 1 )
    END IF
    !
END SUBROUTINE get_simpson_coeff
!====================================================================================

!====================================================================================
SUBROUTINE get_simpson38_coeff(coeff, n)
    !
    USE kinds, ONLY : DP
    !
    IMPLICIT NONE
    !
    INTEGER    :: n, i, ep
    REAL(DP)   :: coeff(n)
    !
    ep =3*((n-1)/3)+1
    !
    IF (n .GE. 22) THEN
      coeff(1)    = 3.0_DP/8 
      coeff(2)    = 9.0_DP/8 
      coeff(ep-1) = 9.0_DP/8 
      coeff(ep)   = 3.0_DP/8 
      !
      DO i = 3, ep-2
        IF (MOD(i,3) .EQ. 1) THEN
          coeff(i) = 6.0_DP/8
        ELSE
          coeff(i) = 9.0_DP/8
        END IF
      END DO
      !
      DO i = ep+1, n
        coeff(i) = 1.0_DP
      END DO
    ELSE
      CALL errore( 'exx_gs', 'simpson38 rule do not apply for n < 22', 1 )
    END IF
    !
END SUBROUTINE get_simpson38_coeff
!====================================================================================

!==============================================================================
SUBROUTINE getrhol(me_r, ps_r, psi1, psi2, rhome, rhops, inv_omega)
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
#ifdef __CUDA
    attributes(device) :: psi1, psi2, rhops, rhome
#endif
    !
    INTEGER  i, j, k
    !
#ifdef __CUDA
    !$cuf kernel do (3)
#else
    !$omp parallel do private(i,j,k)
#endif
    DO k=ps_r(3),ps_r(6)
      DO j=ps_r(2),ps_r(5)
        DO i=ps_r(1),ps_r(4)
          rhops(i,j,k) = psi1(i,j,k) * psi2(i,j,k) * inv_omega
        END DO
      END DO
    END DO
#ifndef __CUDA
    !$omp end parallel do 
#endif
    !--------------------------------------------------------------------------
    !
    !--------------------------------------------------------------------------
#ifdef __CUDA
    !$cuf kernel do (3)
#else
    !$omp parallel do private(i,j,k)
#endif
    DO k=me_r(3),me_r(6)
      DO j=me_r(2),me_r(5)
        DO i=me_r(1),me_r(4)
          rhome(i,j,k) = psi1(i,j,k)*psi2(i,j,k) * inv_omega
        END DO
      END DO
    END DO
#ifndef __CUDA
    !$omp end parallel do 
#endif
    !--------------------------------------------------------------------------
    !
    RETURN
END SUBROUTINE getrhol
!==============================================================================

!==============================================================================
SUBROUTINE updateforce_loc(nrg, me_r, potpsi, potme, psime1, psime2, tran)
    !
    USE kinds,                   ONLY  : DP
    USE exx_module,              ONLY  : exxalfa
    USE dummy_exx, ONLY  : l2gcb
    USE fft_base,         ONLY  : dfftp
    !
    IMPLICIT NONE
    !
    INTEGER      :: nrg(3)
    INTEGER      :: me_r(6)
    REAL(DP)     :: potpsi(nrg(1),nrg(2),nrg(3))
    REAL(DP)     :: potme(me_r(1):me_r(4),me_r(2):me_r(5),me_r(3):me_r(6))
    REAL(DP)     :: psime1(me_r(1):me_r(4),me_r(2):me_r(5),me_r(3):me_r(6))
    REAL(DP)     :: psime2(me_r(1):me_r(4),me_r(2):me_r(5),me_r(3):me_r(6))
    INTEGER      :: tran(3)
    !----------------------------------------------------------------
    INTEGER      :: gid(3)
    INTEGER      :: lid(3)
    INTEGER      :: i,j,k
    INTEGER      :: gi,gj,gk,ti,tj,tk
#ifdef __CUDA
    attributes (device) :: potpsi,potme,psime1,psime2
#endif
    !
    ti=tran(1);tj=tran(2);tk=tran(3)
    !----------------------------------------------------------------
    ! update potpsi in the global grid (exxalfa is 0.25 for PBE0)
    !----------------------------------------------------------------
#ifdef __CUDA
    !$cuf kernel do (3)
#else
    !$omp parallel do private(i,j,k,gi,gj,gk)
#endif
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
          potpsi(gi,gj,gk) = potpsi(gi,gj,gk) &
                                         - exxalfa*potme(i,j,k)*psime2(i,j,k)
          psime2(i,j,k) = - exxalfa*potme(i,j,k)*psime1(i,j,k)
          !----------------------------------------------------------
        END DO
      END DO
    END DO
    !----------------------------------------------------------------
#ifndef __CUDA
    !$omp end parallel do 
#endif
    !----------------------------------------------------------------

    !----------------------------------------------------------------
    RETURN
    !----------------------------------------------------------------
END SUBROUTINE updateforce_loc
!==============================================================================

!==============================================================================
SUBROUTINE updateforce_slf(nrg, me_r, potpsi, potme, psime, tran)
    !
    USE kinds,                   ONLY  : DP
    USE exx_module,              ONLY  : exxalfa
    USE dummy_exx, ONLY  : l2gcb
    USE fft_base,         ONLY  : dfftp
    !
    IMPLICIT NONE
    !
    INTEGER      :: nrg(3)
    INTEGER      :: me_r(6)
    REAL(DP)     :: potpsi(nrg(1),nrg(2),nrg(3))
    REAL(DP)     :: potme(me_r(1):me_r(4),me_r(2):me_r(5),me_r(3):me_r(6))
    REAL(DP)     :: psime(me_r(1):me_r(4),me_r(2):me_r(5),me_r(3):me_r(6))
    INTEGER      :: tran(3)
    !----------------------------------------------------------------
    INTEGER      :: gid(3)
    INTEGER      :: lid(3)
    INTEGER      :: i,j,k
    !
    INTEGER      :: gi,gj,gk,ti,tj,tk
#ifdef __CUDA
    attributes (device) :: potpsi,potme,psime
#endif
    !
    ti=tran(1);tj=tran(2);tk=tran(3)
    !----------------------------------------------------------------
    ! update potpsi in the global grid (exxalfa is 0.25 for PBE0)
    !----------------------------------------------------------------
#ifdef __CUDA
    !$cuf kernel do (3)
#else
    !$omp parallel do private(i,j,k,gi,gj,gk)
#endif
    DO k=me_r(3),me_r(6)
      DO j=me_r(2),me_r(5)
        DO i=me_r(1),me_r(4)
          !----------------------------------------------------------
          gi = l2gcb(dfftp%nr1,i,ti)
          gj = l2gcb(dfftp%nr2,j,tj)
          gk = l2gcb(dfftp%nr3,k,tk)
          !----------------------------------------------------------
          potpsi(gi,gj,gk) = potpsi(gi,gj,gk) &
                                         - exxalfa*potme(i,j,k)*psime(i,j,k)
          !----------------------------------------------------------
        END DO
      END DO
    END DO
    !----------------------------------------------------------------
#ifndef __CUDA
    !$omp end parallel do 
#endif
    !----------------------------------------------------------------

    !----------------------------------------------------------------
    RETURN
    !----------------------------------------------------------------
END SUBROUTINE updateforce_slf
!==============================================================================

!==============================================================================
SUBROUTINE updateforce_rec(nrg, me_r, potpsi, force, tran)
    !
    USE kinds,                   ONLY  : DP
    USE exx_module,              ONLY  : exxalfa
    USE dummy_exx,               ONLY  : l2gcb
    USE fft_base,                ONLY  : dfftp
    !
    IMPLICIT NONE
    !
    INTEGER      :: nrg(3)
    INTEGER      :: me_r(6)
    REAL(DP)     :: potpsi(nrg(1),nrg(2),nrg(3))
    REAL(DP)     :: force(me_r(1):me_r(4),me_r(2):me_r(5),me_r(3):me_r(6))
    INTEGER      :: tran(3)
    !----------------------------------------------------------------
    INTEGER      :: gi, gj, gk
    INTEGER      :: i,j,k
    INTEGER      :: ti,tj,tk
#ifdef __CUDA
    attributes (device) :: potpsi, force
#endif
    !
    ti=tran(1);tj=tran(2);tk=tran(3)
    !
    !----------------------------------------------------------------
    ! update potpsi in the global grid (exxalfa is 0.25 for PBE0)
    !----------------------------------------------------------------
#ifdef __CUDA
    !$cuf kernel do (3)
#else
    !$omp parallel do private(i,j,k,gi,gj,gk)
#endif
    !----------------------------------------------------------------
    DO k=me_r(3),me_r(6)
      DO j=me_r(2),me_r(5)
        DO i=me_r(1),me_r(4)
          !----------------------------------------------------------
          gi = l2gcb(dfftp%nr1,i,ti)
          gj = l2gcb(dfftp%nr2,j,tj)
          gk = l2gcb(dfftp%nr3,k,tk)
          !----------------------------------------------------------
          potpsi(gi,gj,gk) = potpsi(gi,gj,gk) + force(i,j,k)
          !----------------------------------------------------------
        END DO
      END DO
    END DO
    !----------------------------------------------------------------
#ifndef __CUDA
    !$omp end parallel do 
#endif
    !----------------------------------------------------------------

    !----------------------------------------------------------------
    RETURN
    !----------------------------------------------------------------
END SUBROUTINE updateforce_rec
!==============================================================================

!==============================================================================
SUBROUTINE fcn(m, n, x, fvec, iflag)
    USE kinds,              ONLY: DP                 !double-precision kind (selected_real_kind(14,200))
    USE exx_module,         ONLY: PScubeSL_p
    USE exx_module,         ONLY: coeke
    IMPLICIT NONE
    !--------------------------------------------------------------
    INTEGER,  INTENT(IN)  :: m, n, iflag
    REAL(DP), INTENT(IN)  :: x(n)
    REAL(DP), INTENT(OUT) :: fvec(m)
    
    INTEGER               :: i
    REAL(DP)              :: y(4)
    REAL(DP)              :: gm

    ! here we use PScubeSL_p since it is smaller than PScubeSL_s
    gm = (1-2.0D0/PScubeSL_p)**(0.25)

    y = [coeke(0,1,1), coeke(1,1,1), coeke(2,1,1), coeke(3,1,1)]
    
    fvec(1) = (x(1)**2.0D0 + 3*x(2)**2.0D0 + 3*x(3)**2.0D0 + 3*x(4)**2.0D0) +&
      & ((3*x(2)**2.0D0 + 3*x(3)**2.0D0 + 3*x(4)**2.0D0) + &
      & (6*x(2)*x(3) + 6*x(2)*x(4) + 6*x(3)*x(4)))*2*gm - y(1)*3
    fvec(2) = x(1)*x(2) + x(2)*x(3) + x(3)*x(4) - y(2)
    fvec(3) = x(1)*x(3) + x(2)*x(4) - y(3)
    fvec(4) = x(1)*x(4) - y(4)
END SUBROUTINE fcn
!==============================================================================


!==============================================================================
SUBROUTINE l2gcb(lid,gid,tran)
    !
    USE fft_base,         ONLY  : dfftp
    !
    IMPLICIT NONE
    !
    INTEGER  lid(3),gid(3),tran(3)
    INTEGER  nr1s, nr2s, nr3s
    !
    nr1s=dfftp%nr1 
    nr2s=dfftp%nr2 
    nr3s=dfftp%nr3 
    !
    gid(1) = MOD(lid(1)-tran(1)-1+nr1s, nr1s)+1
    gid(2) = MOD(lid(2)-tran(2)-1+nr2s, nr2s)+1
    gid(3) = MOD(lid(3)-tran(3)-1+nr3s, nr3s)+1
    !
    RETURN
END SUBROUTINE l2gcb
!==============================================================================

