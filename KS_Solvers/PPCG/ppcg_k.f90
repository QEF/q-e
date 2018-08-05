!
SUBROUTINE ppcg_k( h_psi, s_psi, overlap, precondition, &
                 npwx, npw, nbnd, npol, psi, e, btype, &
                 ethr, maxter, notconv, avg_iter, sbsize, rr_step, scf_iter)
  !
  !----------------------------------------------------------------------------
  !
  ! E.V. Ignore btype, use ethr as threshold on subspace residual subspace
  ! SdG  restore btype use in the eigenvalue locking procedure
  !
  USE ppcg_param,         ONLY : DP, stdout
  USE mp,                 ONLY : mp_bcast, mp_root_sum, mp_sum
  USE mp_bands_util,      ONLY : intra_bgrp_comm, inter_bgrp_comm, root_bgrp_id, nbgrp, my_bgrp_id
  USE descriptors,        ONLY : la_descriptor, descla_init, descla_local_dims
  USE parallel_toolkit,   ONLY : zsqmher
  USE mp_diag,            ONLY : ortho_comm, np_ortho, me_ortho, ortho_comm_id, leg_ortho, &
                                 ortho_parent_comm, ortho_cntx, do_distr_diag_inside_bgrp
  !
  IMPLICIT NONE
  COMPLEX (DP), PARAMETER :: C_ONE = (1.D0,0.D0), C_ZERO = (0.D0,0.D0)
  !
  ! ... I/O variables
  !
  LOGICAL,      INTENT(IN)    :: overlap ! whether the eigenvalue problem is a generalized one or not
  INTEGER,      INTENT(IN)    :: npwx, npw, nbnd, npol, maxter
  ! maximum number of PW for wavefunctions
  ! the number of plane waves
  ! number of bands
  ! number of independent spin components of the wfc.
  ! maximum number of iterations
  COMPLEX (DP), INTENT(INOUT) :: psi(npwx*npol,nbnd)
  INTEGER,      INTENT(IN)    :: btype(nbnd) ! one if the corresponding state has to be
                                             ! ...converged to full accuracy, zero otherwise (src/pwcom.f90)
  REAL (DP),    INTENT(INOUT) :: e(nbnd)
  REAL (DP),    INTENT(IN)    :: precondition(npw), ethr
  ! the diagonal preconditioner
  ! the convergence threshold for eigenvalues
  INTEGER,      INTENT(OUT)   :: notconv
  REAL(DP),     INTENT(OUT)   :: avg_iter
  ! number of notconverged elements
  ! average number of iterations in PPCG
  INTEGER,      INTENT(IN)   ::  sbsize, rr_step      ! sub-block size (num. of vectors in sub-blocks)
                                                      ! ...to be used in PPCG block splitting. by default, sbsize=1
                                                      ! run the Rayleigh Ritz procedure every rr_step
  INTEGER, INTENT(IN)         :: scf_iter             ! this variable should be removed in future version, used for timing purpose
  !
  ! ... local variables
  !
  COMPLEX(DP), ALLOCATABLE ::  hpsi(:,:), spsi(:,:), w(:,:), hw(:,:), sw(:,:), p(:,:), hp(:,:), sp(:,:)
  COMPLEX(DP)              ::  buffer(npwx*npol,nbnd), buffer1(npwx*npol,nbnd)
  COMPLEX(DP), ALLOCATABLE ::  K(:,:), K_store(:,:), M(:,:), M_store(:,:), cwork(:)
  REAL(DP), ALLOCATABLE    ::  rwork(:)
  INTEGER,  ALLOCATABLE    ::  iwork(:)
  REAL (DP)                ::  trG, trG1, trdif, trtol, lock_tol
  COMPLEX(DP)              ::  coord_psi(sbsize,sbsize), coord_w(sbsize,sbsize), coord_p(sbsize,sbsize),  &
                               G(nbnd,nbnd), G1(nbnd, nbnd)
  REAL (DP)                ::  D(3*sbsize)
  INTEGER                  ::  nsb, sbsize_last, sbsize3, dimp,                        &
                               l, i, j, iter, total_iter, ierr = 0, info = 0,                       &
                               col_idx(sbsize), lcwork = -1, lrwork = -1, liwork = -1
  INTEGER                  ::  nact, nact_old, act_idx(nbnd), idx(nbnd)
                               ! number of active columns
                               ! number of active columns on previous iteration
                               ! indices of active columns
                               ! auxiliary indices
  INTEGER                  ::  n_start, n_end, my_n ! auxiliary indices for band group parallelization
  INTEGER                  ::  print_info     ! If > 0 then iteration information is printed
  INTEGER                  ::  kdim, kdimx, ipol, ibnd
  LOGICAL                  ::  clean
  REAL(DP), EXTERNAL :: ZLANGE

  EXTERNAL h_psi, s_psi
    ! h_psi(npwx,npw,nvec,psi,hpsi)
    !     calculates H|psi>
    ! s_psi(npwx,npw,nvec,psi,spsi)
    !     calculates S|psi> (if needed)
    !     Vectors psi,hpsi,spsi are dimensioned (npwx,nvec)

  COMPLEX (DP), ALLOCATABLE    ::  Gl(:,:)
  !
  TYPE(la_descriptor)  :: desc
  ! descriptor of the current distributed Gram matrix
  LOGICAL :: la_proc
  ! flag to distinguish procs involved in linear algebra
  INTEGER, ALLOCATABLE :: irc_ip( : )
  INTEGER, ALLOCATABLE :: nrc_ip( : )
  INTEGER, ALLOCATABLE :: rank_ip( :, : )
  ! matrix distribution descriptors
  LOGICAL :: force_repmat    ! = .TRUE. to force replication of the Gram matrices for Cholesky
                             ! Needed if the sizes of these  matrices become too small after locking
  REAL                          :: res_array(maxter)

    res_array     = 0.0
  !
  CALL start_clock( 'ppcg_k' )
  !
  !  ... Initialization and validation
!
  print_info = 0 ! 3 
  sbsize3 = sbsize*3
  kdim  = npwx*(npol-1) + npw
  kdimx = npwx*npol
  clean = (npw < npwx) .AND. ( npol == 2 )
  if (npol> 2) CALL errore( 'ppcg ',' wrong npol value: npol > 2 ', npol )
  if (npol<=0) CALL errore( 'ppcg ',' non positive npol value: errcode = 1+abs(npol) ', 1+abs(npol) )
  !
  nact      =  nbnd
  nact_old  =  nbnd
  lock_tol  =  SQRT( ethr )  ! 1.D-4*SQRT( ethr )  ! 1.0D-7
  trdif     =  -1.D0
  ! trdif is used for stopping. It equals either
  ! the difference of traces of psi'*hpsi on two consecutive
  ! itrations or -1.D0 when this difference is undefined
  ! (e.g., initial iteration or the iteration following RR)
  trG    =  0.D0
  trG1   =  0.D0
  ! trG and trG1 are use to evaluate the difference of traces
  ! of psi'*hpsi on two consecutive iterations
  iter      =  1
  force_repmat = .FALSE.
  !
  CALL allocate_all
  !
  ! ... Compute block residual w = hpsi - psi*(psi'hpsi) (psi is orthonormal on input)
  !
  call start_clock('ppcg:hpsi')
                                                         if (clean)  psi(npw+1:npwx,:) = C_ZERO
  CALL h_psi( npwx, npw, nbnd, psi, hpsi )             ; if (clean) hpsi(npw+1:npwx,:) = C_ZERO
  if (overlap) CALL s_psi( npwx, npw, nbnd, psi, spsi) ; if (clean) spsi(npw+1:npwx,:) = C_ZERO

  avg_iter = 1.d0
  call stop_clock('ppcg:hpsi')
  !
  !     G = psi'hpsi
  call start_clock('ppcg:zgemm')
  G = C_ZERO
  CALL divide(inter_bgrp_comm,nbnd,n_start,n_end); my_n = n_end - n_start + 1; !write (*,*) nbnd,n_start,n_end
  if (n_start .le. n_end) &
  CALL ZGEMM('C','N', nbnd, my_n, kdim, C_ONE, psi, kdimx, hpsi(1,n_start), kdimx, C_ZERO, G(1,n_start), nbnd)
  CALL mp_sum( G, inter_bgrp_comm )
  !
  CALL mp_sum( G, intra_bgrp_comm )
  call stop_clock('ppcg:zgemm')
  !
  !    w = hpsi - spsi*G
  call start_clock('ppcg:zgemm')
  CALL threaded_assign( w, hpsi, kdimx, nact, bgrp_root_only=.true. )
  CALL divide(inter_bgrp_comm,nbnd,n_start,n_end); my_n = n_end - n_start + 1; !write (*,*) nbnd,n_start,n_end
  if (overlap) then
     if (n_start .le. n_end) &
     CALL ZGEMM('N','N',kdim, nbnd, my_n, -C_ONE,spsi(1,n_start), kdimx, G(n_start,1), nbnd, C_ONE, w, kdimx)
  else
     if (n_start .le. n_end) &
     CALL ZGEMM('N','N',kdim, nbnd, my_n, -C_ONE, psi(1,n_start), kdimx, G(n_start,1), nbnd, C_ONE, w, kdimx)
  end if
  CALL mp_sum( w, inter_bgrp_comm )
  call stop_clock('ppcg:zgemm')
  !
  !
  ! ... Lock converged eigenpairs (set up act_idx and nact and store current nact in nact_old)
  call start_clock('ppcg:lock')
  nact_old = nact;
  CALL lock_epairs(kdim, nbnd, btype, w, kdimx, lock_tol, nact, act_idx)
  call stop_clock('ppcg:lock')
  !
  ! ... Set up iteration parameters after locking
  CALL setup_param
  !
  G1(1:nact,1:nact) = G(act_idx(1:nact), act_idx(1:nact))
  trG = get_trace( G1, nbnd, nact )
  !
  ! Print initial info ...
  IF (print_info >= 1)  THEN
     WRITE(stdout, '("Ethr: ",1pD9.2,", npw: ", I10, ", nbnd: ", I10, " , "   &
                "maxter: ",I5, ", sbsize:  ", I10,", nsb: ", I10 ,", nact: ", I10, ", trtol: ", 1pD9.2 )')  &
                ethr, npw, nbnd, maxter, sbsize, nsb, nact, trtol
     IF (print_info == 3) THEN
        CALL print_rnrm
        WRITE(stdout,'("Res. norm:  ", 1pD9.2)') res_array(iter)
     END IF
     CALL flush( stdout )
  END IF
  !
  !---Begin the main loop
  !
  DO WHILE ( ((trdif > trtol) .OR. (trdif == -1.D0))  .AND. (iter <= maxter) .AND. (nact > 0) )
     !
     ! ... apply the diagonal preconditioner
     !
     !$omp parallel do collapse(2)
     DO j = 1, nact
     do ipol=0,npol-1
        w(1+npwx*ipol:npw+npwx*ipol, act_idx(j)) = w(1+npwx*ipol:npw+npwx*ipol,act_idx(j)) / precondition(1:npw)
     end do
     END DO
     !$omp end parallel do
     !
     call start_clock('ppcg:zgemm')
     call threaded_assign( buffer, w, kdimx, nact, act_idx )
     G(1:nbnd,1:nact) = C_ZERO
     CALL divide(inter_bgrp_comm,nbnd,n_start,n_end); my_n = n_end - n_start + 1; !write (*,*) nbnd,n_start,n_end
     if (overlap) then
        if (n_start .le. n_end) &
        CALL ZGEMM( 'C','N', my_n, nact, kdim, C_ONE, spsi(1,n_start), kdimx, buffer, kdimx, C_ZERO, G(n_start,1), nbnd )
     else
        if (n_start .le. n_end) &
        CALL ZGEMM( 'C','N', my_n, nact, kdim, C_ONE,psi(1,n_start), kdimx, buffer, kdimx, C_ZERO, G(n_start,1), nbnd )
     end if
     CALL mp_sum( G(1:nbnd,1:nact), inter_bgrp_comm )
     !
     CALL mp_sum( G(1:nbnd,1:nact), intra_bgrp_comm )
     call stop_clock('ppcg:zgemm')
     !
     !     w = w - psi*G
     call start_clock('ppcg:zgemm')
     call threaded_assign( buffer, w, kdimx, nact, act_idx, bgrp_root_only=.true. )
     if (n_start .le. n_end) &
     CALL ZGEMM('N','N', kdim, nact, my_n, -C_ONE, psi(1,n_start), kdimx, G(n_start,1), nbnd, C_ONE, buffer, kdimx)
     CALL mp_sum( buffer(:,1:nact), inter_bgrp_comm )
     w(:,act_idx(1:nact)) = buffer(:,1:nact)
     call stop_clock('ppcg:zgemm')
     !
     ! ... Compute h*w
     call start_clock('ppcg:hpsi')
     call threaded_assign( buffer1, w, kdimx, nact, act_idx )
     CALL h_psi( npwx, npw, nact, buffer1, buffer )       ; if (clean) buffer (npw+1:npwx,1:nact) = C_ZERO
     hw(:,act_idx(1:nact)) = buffer(:,1:nact)
     if (overlap) then ! ... Compute s*w
        CALL s_psi( npwx, npw, nact, buffer1, buffer )    ; if (clean) buffer(npw+1:npwx,1:nact) = C_ZERO
        sw(:,act_idx(1:nact)) = buffer(:,1:nact)
     end if
     avg_iter = avg_iter + nact/dble(nbnd)
     call stop_clock('ppcg:hpsi')
     !
     ! ... orthogonalize p against psi and w
!ev     IF ( MOD(iter, rr_step) /= 1 ) THEN    ! In this case, P is skipped after each RR
     IF ( iter  /=  1 ) THEN
        !  G = spsi'p
        call start_clock('ppcg:zgemm')
        G(1:nact,1:nact) = C_ZERO
        if (overlap) then
           call threaded_assign( buffer, spsi, kdimx, nact, act_idx )
        else
           call threaded_assign( buffer,  psi, kdimx, nact, act_idx )
        end if
        call threaded_assign( buffer1, p, kdimx, nact, act_idx )
        CALL divide(inter_bgrp_comm,nact,n_start,n_end); my_n = n_end - n_start + 1; !write (*,*) nact,n_start,n_end
        if (n_start .le. n_end) &
        CALL ZGEMM('C','N', my_n, nact, kdim, C_ONE, buffer(1,n_start), kdimx, buffer1, kdimx, C_ZERO, G(n_start,1), nbnd)
        CALL mp_sum( G(1:nact,1:nact), inter_bgrp_comm )
        !
        CALL mp_sum( G(1:nact,1:nact), intra_bgrp_comm )
        call stop_clock('ppcg:zgemm')
        !
        ! p = p - psi*G, hp = hp - hpsi*G, sp = sp - spsi*G
        call start_clock('ppcg:zgemm')
        call threaded_assign( buffer, p, kdimx, nact, act_idx, bgrp_root_only = .true. )
        call threaded_assign( buffer1, psi, kdimx, nact, act_idx )
        if (n_start .le. n_end) & ! could be done differently
        CALL ZGEMM('N','N', kdim, nact, my_n, -C_ONE, buffer1(1,n_start), kdimx, G(n_start,1), nbnd, C_ONE, buffer, kdimx)
        CALL mp_sum( buffer(:,1:nact), inter_bgrp_comm )
        p(:,act_idx(1:nact)) = buffer(:,1:nact)
        call stop_clock('ppcg:zgemm')
        !
        call start_clock('ppcg:zgemm')
        call threaded_assign( buffer, hp, kdimx, nact, act_idx, bgrp_root_only = .true. )
        call threaded_assign( buffer1, hpsi, kdimx, nact, act_idx )
        if (n_start .le. n_end) &
        CALL ZGEMM('N','N', kdim, nact, my_n, -C_ONE, buffer1(1,n_start), kdimx, G(n_start,1), nbnd, C_ONE, buffer, kdimx)
        CALL mp_sum( buffer(:,1:nact), inter_bgrp_comm )
        hp(:,act_idx(1:nact)) = buffer(:,1:nact)
        call stop_clock('ppcg:zgemm')
        !
        if (overlap) then
           call start_clock('ppcg:zgemm')
           call threaded_assign( buffer, sp, kdimx, nact, act_idx, bgrp_root_only = .true. )
           call threaded_assign( buffer1, spsi, kdimx, nact, act_idx )
           if (n_start .le. n_end) &
           CALL ZGEMM('N','N', kdim, nact, my_n, -C_ONE, buffer1(1,n_start), kdimx, G(n_start,1), nbnd, C_ONE, buffer, kdimx)
           CALL mp_sum( buffer(:,1:nact), inter_bgrp_comm )
           sp(:,act_idx(1:nact)) = buffer(:,1:nact)
           call stop_clock('ppcg:zgemm')
        end if
     END IF
     !
     !  ... for each sub-block construct the small projected matrices K and M
     !      and store in K_store and M_store
     !
     K_store = C_ZERO
     M_store = C_ZERO
     !
     CALL divide(inter_bgrp_comm,nsb,n_start,n_end); my_n = n_end - n_start + 1; !write (*,*) nsb,n_start,n_end
     DO j = n_start, n_end
        !
        ! Get size of the sub-block and define indices of the corresponding columns
        IF ( j < nsb )  THEN
           l = sbsize
        ELSE
           l = sbsize_last
        END IF
        col_idx(1:l) = act_idx(  (/ (i, i = (j-1)*sbsize + 1, (j-1)*sbsize + l) /) )
        !
        ! ... form the local Gramm matrices (K,M)
        K = C_ZERO
        M = C_ZERO
        !
        call start_clock('ppcg:zgemm')
        call threaded_assign( buffer, psi, kdimx, l, col_idx )
        call threaded_assign( buffer1, hpsi, kdimx, l, col_idx )
        CALL ZGEMM('C','N', l, l, kdim, C_ONE, buffer, kdimx, buffer1, kdimx, C_ZERO, K, sbsize3)
        !
        if (overlap) then
           call threaded_assign( buffer1, spsi, kdimx, l, col_idx ); if (clean) buffer1(npw+1:npwx,1:l) = C_ZERO
        else
           call threaded_assign( buffer1, buffer, kdimx, l )
        end if
        CALL ZGEMM('C','N', l, l, kdim, C_ONE, buffer, kdimx, buffer1, kdimx, C_ZERO, M, sbsize3)
        !
        ! ---
        call threaded_assign( buffer, w, kdimx, l, col_idx )
        call threaded_assign( buffer1, hw, kdimx, l, col_idx )
        CALL ZGEMM('C','N', l, l, kdim, C_ONE, buffer, kdimx, buffer1, kdimx, C_ZERO, K(l+1, l+1), sbsize3)
        !
        if (overlap) then
           call threaded_assign( buffer1, sw, kdimx, l, col_idx )
        else
           call threaded_assign( buffer1, buffer, kdimx, l )
        end if
        CALL ZGEMM('C','N', l, l, kdim, C_ONE, buffer, kdimx, buffer1, kdimx, C_ZERO, M(l+1, l+1 ), sbsize3)
        !
        ! ---
        call threaded_assign( buffer, psi, kdimx, l, col_idx )
        call threaded_assign( buffer1, hw, kdimx, l, col_idx ) ; if (clean) buffer1(npw+1:npwx,1:l) = C_ZERO
        CALL ZGEMM('C','N', l, l, kdim, C_ONE, buffer, kdimx, buffer1, kdimx, C_ZERO, K(1, l+1), sbsize3)
        !
        if (overlap) then
           call threaded_assign( buffer1, sw, kdimx, l, col_idx )
        else
           call threaded_assign( buffer1,  w, kdimx, l, col_idx )
        end if
        CALL ZGEMM('C','N', l, l, kdim, C_ONE, buffer, kdimx, buffer1, kdimx, C_ZERO, M(1, l+1), sbsize3)
        call stop_clock('ppcg:zgemm')
        !
        ! ---
        !
!ev        IF ( MOD(iter,rr_step) /= 1 ) THEN   ! In this case, P is skipped after each RR
        IF ( iter  /= 1 ) THEN
          call start_clock('ppcg:zgemm')
          call threaded_assign( buffer,  p, kdimx, l, col_idx )
          call threaded_assign( buffer1,  hp, kdimx, l, col_idx )
          CALL ZGEMM('C','N', l, l, kdim, C_ONE, buffer, kdimx, buffer1, kdimx, C_ZERO, K(2*l + 1, 2*l+1), sbsize3)
          !
          if (overlap) then
             call threaded_assign( buffer1,  sp, kdimx, l, col_idx )
          else
             call threaded_assign( buffer1,  buffer, kdimx, l )
          end if
          CALL ZGEMM('C','N', l, l, kdim, C_ONE, buffer, kdimx, buffer1, kdimx, C_ZERO, M(2*l + 1, 2*l+1), sbsize3)
          !
          ! ---
          call threaded_assign( buffer,  psi, kdimx, l, col_idx )
          call threaded_assign( buffer1,  hp, kdimx, l, col_idx )
          CALL ZGEMM('C','N', l, l, kdim, C_ONE, buffer, kdimx, buffer1, kdimx, C_ZERO, K(1, 2*l+1), sbsize3)
          !
          if (overlap) then
             call threaded_assign( buffer1,  sp, kdimx, l, col_idx )
          else
             call threaded_assign( buffer1,  p, kdimx, l, col_idx )
          end if
          CALL ZGEMM('C','N', l, l, kdim, C_ONE, buffer, kdimx, buffer1, kdimx, C_ZERO, M(1, 2*l+1), sbsize3)
          call stop_clock('ppcg:zgemm')
          !
          ! ---
          !
          call start_clock('ppcg:zgemm')
          call threaded_assign( buffer,  w, kdimx, l, col_idx )
          call threaded_assign( buffer1,  hp, kdimx, l, col_idx )
          CALL ZGEMM('C','N', l, l, kdim, C_ONE, buffer, kdimx, buffer1, kdimx, C_ZERO, K(l+1, 2*l+1), sbsize3)
          !
          if (overlap) then
             call threaded_assign( buffer1,  sp, kdimx, l, col_idx )
          else
             call threaded_assign( buffer1,  p, kdimx, l, col_idx )
          end if
          CALL ZGEMM('C','N', l, l, kdim, C_ONE, buffer, kdimx, buffer1, kdimx, C_ZERO, M(l+1, 2*l+1), sbsize3)
          call stop_clock('ppcg:zgemm')
          !
        END IF
        !
        ! ... store the projected matrices
        K_store(:, (j-1)*sbsize3 + 1 : j*sbsize3 ) = K
        M_store(:, (j-1)*sbsize3 + 1 : j*sbsize3 ) = M
        !
     END DO
     CALL mp_sum(K_store,inter_bgrp_comm)
     CALL mp_sum(M_store,inter_bgrp_comm)
     !
     CALL mp_sum(K_store,intra_bgrp_comm)
     CALL mp_sum(M_store,intra_bgrp_comm)
     !
     ! ... perform nsb 'separate RQ minimizations' and update approximate subspace

     idx(:) = 0 ! find the inactive columns to be kept by root_bgrp_id only
     if (my_bgrp_id == root_bgrp_id) then
        idx(1:nbnd) = 1 ; idx(act_idx(1:nact)) = 0
     end if
     DO j = n_start, n_end
       !
       ! Get size of the sub-block and define indices of the corresponding columns
       IF ( j < nsb )  THEN
          l = sbsize
       ELSE
          l = sbsize_last
       END IF

       col_idx(1:l) = act_idx( (/ (i, i = (j-1)*sbsize + 1, (j-1)*sbsize + l) /)  )
       !
       K = K_store(:, (j-1)*sbsize3 + 1 : j*sbsize3)
       M = M_store(:, (j-1)*sbsize3 + 1 : j*sbsize3)
       !
       lcwork =       6*sbsize +  9*sbsize**2
       lrwork =  1 + 15*sbsize + 18*sbsize**2
       liwork =  3 + 15*sbsize
!ev       IF ( MOD(iter, rr_step) /= 1) THEN    ! set the dimension of the separate projected eigenproblem
       IF ( iter /= 1 ) THEN    ! set the dimension of the separate projected eigenproblem
          dimp = 3*l
       ELSE
          dimp = 2*l
       END IF
       !
       CALL ZHEGVD(1, 'V','U', dimp, K, sbsize3, M, sbsize3, D, cwork, lcwork, rwork, lrwork, iwork, liwork, info)
       IF (info /= 0) THEN
       ! reset the matrix and try again with psi and w only
          K = K_store(:, (j-1)*sbsize3 + 1 : j*sbsize3)
          M = M_store(:, (j-1)*sbsize3 + 1 : j*sbsize3)
          dimp = 2*l
          CALL ZHEGVD(1, 'V','U', dimp, K, sbsize3, M, sbsize3, D, cwork, lcwork, rwork, lrwork, iwork, liwork, info)
          IF (info /= 0) THEN
             CALL errore( 'ppcg ',' zhegvd failed ', info )
             STOP
          END IF
       END IF
       !
       coord_psi(1 : l, 1 : l) = K(1 : l, 1 : l)
       coord_w(1 : l, 1 : l) = K(l+1 : 2*l, 1 : l)
       !
       ! ... update the sub-block of P and AP
!ev       IF ( MOD(iter, rr_step) /= 1 ) THEN
!sdg      IF ( iter /= 1 ) THEN
       IF ( dimp == 3*l ) THEN
          !
          coord_p(1 : l, 1 : l) = K(2*l+1 : 3*l, 1 : l)
          !
          call start_clock('ppcg:zgemm')
          call threaded_assign( buffer1,  p, kdimx, l, col_idx )
          CALL ZGEMM('N','N', kdim, l, l, C_ONE, buffer1, kdimx, coord_p, sbsize, C_ZERO, buffer, kdimx)
          call threaded_assign( buffer1,  w, kdimx, l, col_idx )
          CALL ZGEMM('N','N', kdim, l, l, C_ONE, buffer1, kdimx, coord_w, sbsize, C_ONE, buffer, kdimx)
          p(:,col_idx(1:l))  = buffer(:,1:l)
          call stop_clock('ppcg:zgemm')
          !
          call start_clock('ppcg:zgemm')
          call threaded_assign( buffer1,  hp, kdimx, l, col_idx )
          CALL ZGEMM('N','N', kdim, l, l, C_ONE, buffer1, kdimx, coord_p, sbsize, C_ZERO, buffer, kdimx)
          call threaded_assign( buffer1,  hw, kdimx, l, col_idx )
          CALL ZGEMM('N','N', kdim, l, l, C_ONE, buffer1, kdimx, coord_w, sbsize, C_ONE, buffer, kdimx)
          hp(:,col_idx(1:l))  = buffer(:,1:l)
          call stop_clock('ppcg:zgemm')
          !
          if (overlap) then
             call start_clock('ppcg:zgemm')
             call threaded_assign( buffer1,  sp, kdimx, l, col_idx )
             CALL ZGEMM('N','N', kdim, l, l, C_ONE, buffer1, kdimx, coord_p, sbsize, C_ZERO, buffer, kdimx)
             call threaded_assign( buffer1,  sw, kdimx, l, col_idx )
             CALL ZGEMM('N','N', kdim, l, l, C_ONE, buffer1, kdimx, coord_w, sbsize, C_ONE, buffer, kdimx)
             sp(:,col_idx(1:l))  = buffer(:,1:l)
             call stop_clock('ppcg:zgemm')
          end if
       ELSE
          !
          call start_clock('ppcg:zgemm')
          call threaded_assign( buffer1,  w, kdimx, l, col_idx )
          CALL ZGEMM('N','N', kdim, l, l, C_ONE, buffer1, kdimx, coord_w, sbsize, C_ZERO, buffer, kdimx)
          p(:,col_idx(1:l)) = buffer(:, 1:l)
          call stop_clock('ppcg:zgemm')
          !
          call start_clock('ppcg:zgemm')
          call threaded_assign( buffer1,  hw, kdimx, l, col_idx )
          CALL ZGEMM('N','N', kdim, l, l, C_ONE, buffer1, kdimx, coord_w, sbsize, C_ZERO, buffer, kdimx)
          hp(:,col_idx(1:l)) = buffer(:, 1:l)
          call stop_clock('ppcg:zgemm')
          !
          if (overlap) then
             call start_clock('ppcg:zgemm')
             call threaded_assign( buffer1,  sw, kdimx, l, col_idx )
             CALL ZGEMM('N','N', kdim, l, l, C_ONE, buffer1, kdimx, coord_w, sbsize, c_ZERO, buffer, kdimx)
             sp(:,col_idx(1:l)) = buffer(:, 1:l)
             call stop_clock('ppcg:zgemm')
          end if
       END IF
       !
       ! Update the sub-blocks of psi and hpsi (and spsi)
       call start_clock('ppcg:zgemm')
       call threaded_assign( buffer1,  psi, kdimx, l, col_idx )
       CALL ZGEMM('N','N', kdim, l, l, C_ONE, buffer1, kdimx, coord_psi, sbsize, C_ZERO, buffer, kdimx)
       psi(:, col_idx(1:l))  = buffer(:,1:l) + p(:,col_idx(1:l))
       call stop_clock('ppcg:zgemm')
       !
       call start_clock('ppcg:zgemm')
       call threaded_assign( buffer1,  hpsi, kdimx, l, col_idx )
       CALL ZGEMM('N','N', kdim, l, l, C_ONE, buffer1, kdimx, coord_psi, sbsize, C_ZERO, buffer, kdimx)
       hpsi(:,col_idx(1:l)) = buffer(:,1:l) + hp(:,col_idx(1:l))
       call stop_clock('ppcg:zgemm')
       !
       if (overlap) then
          call start_clock('ppcg:zgemm')
          call threaded_assign( buffer1,  spsi, kdimx, l, col_idx )
          CALL ZGEMM('N','N', kdim, l, l, C_ONE, buffer1, kdimx, coord_psi, sbsize, C_ZERO, buffer, kdimx)
          spsi(:,col_idx(1:l)) = buffer(:,1:l) + sp(:,col_idx(1:l))
          call stop_clock('ppcg:zgemm')
       end if
       !
       idx(col_idx(1:l)) = 1 ! keep track of which columns this bgrp has acted on
     END DO  ! end 'separate RQ minimizations'
     ! set to zero the columns not assigned to this bgrp, inactive colums are assigned to root_bgrp
     do j=1,nbnd
        if (idx(j)==0) then
           psi (:,j) = C_ZERO ; hpsi (:,j) = C_ZERO ; p(:,j) = C_ZERO ; hp(:,j) = C_ZERO
        end if
     end do
     CALL mp_sum(psi ,inter_bgrp_comm)
     CALL mp_sum(hpsi,inter_bgrp_comm)
     CALL mp_sum(p ,inter_bgrp_comm)
     CALL mp_sum(hp,inter_bgrp_comm)
     if (overlap) then
        do j=1,nbnd
           if (idx(j)==0) then
              spsi (:,j) = C_ZERO ; sp(:,j) = C_ZERO
           end if
        end do
        CALL mp_sum(spsi,inter_bgrp_comm)
        CALL mp_sum(sp,inter_bgrp_comm)
    end if
    !
    !
    ! ... Perform the RR procedure every rr_step
!    IF ( (MOD(iter, rr_step) == 0) .AND. (iter /= maxter) ) THEN
    IF ( MOD(iter, rr_step) == 0 ) THEN
       !
       call start_clock('ppcg:RR')
       CALL extract_epairs_dmat(kdim, nbnd, kdimx, e, psi, hpsi, spsi )
       call stop_clock('ppcg:RR')
       !
       IF (print_info >= 2) WRITE(stdout, *) 'RR has been invoked.' ; !CALL flush( stdout )
       !
       ! ... Compute the new residual vector block by evaluating
       !     residuals for individual eigenpairs in psi and e
       if (overlap) then
          !$omp parallel do collapse(2)
          DO j = 1, nbnd
             DO i = 1, kdim
                w( i ,j ) = hpsi( i, j ) - spsi( i, j )*e( j )
             END DO
          END DO
          !$omp end parallel do
       else
          !$omp parallel do collapse(2)
          DO j = 1, nbnd
             DO i = 1, kdim
                w( i, j ) = hpsi( i, j ) -  psi( i, j )*e( j )
             END DO
          END DO
          !$omp end parallel do
       end if
       !
       ! ... Lock converged eigenpairs (set up act_idx and nact)
       call start_clock('ppcg:lock')
       nact_old = nact;
       CALL lock_epairs(kdim, nbnd, btype, w, kdimx, lock_tol, nact, act_idx)
       call stop_clock('ppcg:lock')
       !
       ! ... Set up iteration parameters after locking
       CALL setup_param
       !
       !
       trG1  =  0.D0
       trdif = -1.D0
       trG   = SUM( e(act_idx(1:nact)) )
       !
    ELSE
       !
!! EV begin
      ! ... orthogonalize psi and update hpsi accordingly
      IF ( .NOT. force_repmat ) THEN
         !
         call threaded_assign( buffer,  psi, kdimx, nact, act_idx )
         if (overlap) then
            call threaded_assign( buffer1,  spsi, kdimx, nact, act_idx )
         else
            call threaded_assign( buffer1,  buffer, kdimx, nact )
         end if
         !
         call start_clock('ppcg:cholQR')
         CALL cholQR_dmat(kdim, nact, buffer, buffer1, kdimx, Gl, desc)
         call stop_clock('ppcg:cholQR')
         !
         psi(:,act_idx(1:nact)) = buffer(:,1:nact)
         !
         call threaded_assign( buffer1,  hpsi, kdimx, nact, act_idx )
         call start_clock('ppcg:ZTRSM')
         CALL zgemm_dmat( kdim, nact, kdimx, desc, C_ONE, buffer1, Gl, C_ZERO, buffer )
         call stop_clock('ppcg:ZTRSM')
         !
         hpsi(:,act_idx(1:nact)) = buffer(:,1:nact)
         !
         if (overlap) then
            call threaded_assign( buffer1,  spsi, kdimx, nact, act_idx )
            call start_clock('ppcg:ZTRSM')
            CALL zgemm_dmat( kdim, nact, kdimx, desc, C_ONE, buffer1, Gl, C_ZERO, buffer )
            call stop_clock('ppcg:ZTRSM')
            !
            spsi(:,act_idx(1:nact)) = buffer(:,1:nact)
         end if
      ELSE
         !
         call threaded_assign( buffer,  psi, kdimx, nact, act_idx )
         if (overlap) then
            call threaded_assign( buffer1,  spsi, kdimx, nact, act_idx )
         else
            call threaded_assign( buffer1,  buffer, kdimx, nact )
         end if
         !
         call start_clock('ppcg:cholQR')
         CALL cholQR(kdim, nact, buffer, buffer1, kdimx, G, nbnd)
         call stop_clock('ppcg:cholQR')
         !
         psi(:,act_idx(1:nact)) = buffer(:,1:nact)
         !
         call threaded_assign( buffer,  hpsi, kdimx, nact, act_idx )
         !
         call start_clock('ppcg:ZTRSM')
         CALL ZTRSM('R', 'U', 'N', 'N', kdim, nact, C_ONE, G, nbnd, buffer, kdimx)
         call stop_clock('ppcg:ZTRSM')
         !
         hpsi(:,act_idx(1:nact)) = buffer(:,1:nact)
         !
         if (overlap) then
            call threaded_assign( buffer,  spsi, kdimx, nact, act_idx )
            !
            call start_clock('ppcg:ZTRSM')
            CALL ZTRSM('R', 'U', 'N', 'N', kdim, nact, C_ONE, G, nbnd, buffer, kdimx)
            call stop_clock('ppcg:ZTRSM')
            !
            spsi(:,act_idx(1:nact)) = buffer(:,1:nact)
         end if
      END IF
!! EV end
       !
       ! ... Compute the new subspace residual for active columns
       !
       !  G = psi'hpsi
       call threaded_assign( buffer,  psi, kdimx, nact, act_idx )
       call threaded_assign( buffer1,  hpsi, kdimx, nact, act_idx )
       call start_clock('ppcg:zgemm')
       G = C_ZERO
       CALL divide(inter_bgrp_comm,nact,n_start,n_end); my_n = n_end - n_start + 1; !write (*,*) nact,n_start,n_end
       if (n_start .le. n_end) &
       CALL ZGEMM('C','N', nact, my_n, kdim, C_ONE, buffer, kdimx, buffer1(1,n_start), kdimx, C_ZERO, G(1,n_start), nbnd)
       CALL mp_sum(G(1:nact,1:nact), inter_bgrp_comm)
       !
       CALL mp_sum(G(1:nact,1:nact), intra_bgrp_comm)
       call stop_clock('ppcg:zgemm')
       !
       ! w = hpsi - spsi*G
       call threaded_assign( buffer,  hpsi, kdimx, nact, act_idx, bgrp_root_only=.true. )
       if (overlap) then
          call threaded_assign( buffer1,  spsi, kdimx, nact, act_idx )
       else
          call threaded_assign( buffer1,   psi, kdimx, nact, act_idx )
       end if
       call start_clock('ppcg:zgemm')
       if (n_start .le. n_end) &
       CALL ZGEMM('N','N', kdim, nact, my_n, -C_ONE, buffer1(1,n_start), kdimx, G(n_start,1), nbnd, C_ONE, buffer, kdimx)
       CALL mp_sum( buffer(:,1:nact), inter_bgrp_comm )
       w(:,act_idx(1:nact)) = buffer(:,1:nact);

       call stop_clock('ppcg:zgemm')
       !
       ! ... Compute trace of the projected matrix on current iteration
       ! trG1  = get_trace( G(1:nact, 1:nact), nact )
       trG1  = get_trace( G, nbnd, nact )
       trdif = ABS(trG1 - trG)
       trG   = trG1
       !
    END IF
    !
    ! Print iteration info ...
    IF (print_info >= 1) THEN
       WRITE(stdout, '("iter: ", I5, " nact = ", I5, ", trdif = ", 1pD9.2, ", trtol = ", 1pD9.2 )') iter, nact, trdif, trtol
       IF (print_info == 3) THEN
          CALL print_rnrm
          WRITE(stdout,'("Res. norm:  ", 1pD9.2)') res_array(iter)
       END IF
       CALL flush( stdout )
    END IF
    !
    total_iter = iter
    iter   = iter + 1
    !
    !
 END DO   !---End the main loop
 !
! IF (nact > 0) THEN
 IF ( MOD(iter-1, rr_step) /= 0 ) THEN        ! if RR has not just been performed
 ! if nact==0 then the RR has just been performed
 ! in the main loop
    call start_clock('ppcg:RR')
    CALL extract_epairs_dmat(kdim, nbnd, kdimx, e, psi, hpsi, spsi )
    call stop_clock('ppcg:RR')
    !
    ! ... Compute residuals
    if (overlap) then
       !$omp parallel do collapse(2)
       DO j = 1, nbnd
          DO i =1, kdim
             w( i, j ) = hpsi( i, j ) - spsi( i, j )*e( j )
          END DO
       END DO
       !$omp end parallel do
    else
       !$omp parallel do collapse(2)
       DO j = 1, nbnd
          DO i =1, kdim
             w( i, j ) = hpsi( i, j ) -  psi( i, j )*e( j )
          END DO
       END DO
       !$omp end parallel do
    end if
    !
    ! ... Get the number of converged eigenpairs and their indices
    !     Note: The tolerance is 10*lock_tol, i.e., weaker tham lock_tol
! E.V. notconv issue should be addressed
    call start_clock('ppcg:lock')
    CALL lock_epairs(kdim, nbnd, btype, w, kdimx, 10*lock_tol, nact, act_idx)
    call stop_clock('ppcg:lock')
    !
 END IF
 !
! E.V. notconv issue comment
 notconv = 0 ! nact
 !
 IF (print_info >= 1) THEN
    WRITE(stdout, *) '-----------PPCG result summary ...  ----------------'
    WRITE(stdout, '("avg_iter: ", f6.2,  ", notconv: ", I5)') avg_iter, notconv
    CALL flush( stdout )
 END IF
 !
 CALL deallocate_all
 !
 CALL stop_clock( 'ppcg_k' )
 !

!!!EV-BANDS
if (print_info == 3) then
    write (stdout,'(1pD9.2)') ( res_array(j), j=1,maxter )
end if

 !
 RETURN
   !
CONTAINS
  !
  SUBROUTINE allocate_all
    !
    ! This subroutine allocates memory for the eigensolver
    !
    INTEGER :: nx
    ! maximum local block dimension
    !
    ALLOCATE ( hpsi(kdimx,nbnd), stat = ierr )
    IF (ierr /= 0) CALL errore( 'ppcg ',' cannot allocate hpsi ', ABS(ierr) )
    if (overlap) ALLOCATE ( spsi(kdimx,nbnd), stat = ierr )
    IF (ierr /= 0) CALL errore( 'ppcg ',' cannot allocate spsi ', ABS(ierr) )
    ALLOCATE ( w(kdimx,nbnd), hw(kdimx,nbnd), stat = ierr )
    IF (ierr /= 0) CALL errore( 'ppcg ',' cannot allocate w and hw ', ABS(ierr) )
    if (overlap) ALLOCATE ( sw(kdimx,nbnd), stat = ierr )
    IF (ierr /= 0) CALL errore( 'ppcg ',' cannot allocate sw ', ABS(ierr) )
    ALLOCATE ( p(kdimx,nbnd), hp(kdimx,nbnd), stat = ierr )
    IF (ierr /= 0) CALL errore( 'ppcg ',' cannot allocate p and hp ', ABS(ierr) )
    if (overlap) ALLOCATE ( sp(kdimx,nbnd), stat = ierr )
    IF (ierr /= 0) CALL errore( 'ppcg ',' cannot allocate sp ', ABS(ierr) )
    ALLOCATE ( K(sbsize3, sbsize3), M(sbsize3,sbsize3), stat = ierr )
    IF (ierr /= 0) CALL errore( 'ppcg ',' cannot allocate K and M ', ABS(ierr) )
    ALLOCATE ( cwork( 1 + 18*sbsize + 18*sbsize**2 ), rwork( 1 + 18*sbsize + 18*sbsize**2 ), iwork(3 + 15*sbsize), stat = ierr )
    IF (ierr /= 0) CALL errore( 'ppcg ',' cannot allocate lapack work arrays ', ABS(ierr) )
    !
    ALLOCATE( irc_ip( np_ortho(1) ), STAT=ierr )
    IF( ierr /= 0 ) CALL errore( 'ppcg ',' cannot allocate irc_ip ', ABS(ierr) )
    ALLOCATE( nrc_ip( np_ortho(1) ), STAT=ierr )
    IF( ierr /= 0 ) CALL errore( 'ppcg ',' cannot allocate nrc_ip ', ABS(ierr) )
    ALLOCATE( rank_ip( np_ortho(1), np_ortho(2) ), STAT=ierr )
    IF( ierr /= 0 ) CALL errore( 'ppcg ',' cannot allocate rank_ip ', ABS(ierr) )
    !
    CALL desc_init( nbnd, desc, irc_ip, nrc_ip  )
    !
    nx = desc%nrcx
    !
    IF ( la_proc ) THEN
       ALLOCATE( Gl( nx, nx ), STAT=ierr )
    ELSE
       ALLOCATE( Gl( 1, 1 ), STAT=ierr )
    END IF
    IF( ierr /= 0 ) CALL errore( 'ppcg ',' cannot allocate Gl ', ABS(ierr) )

  END SUBROUTINE allocate_all
  !
  !
  !
  SUBROUTINE cholQR(n, k, X, SX, ldx,  R, ldr)
    !
    !  This subroutine orthogonalizes X using the Choleski decomposition.
    !  The matrix X'X and its Choleski decomposition is replicated on all processors.
    !  TBD: If the Choleslki orthogonalization fails, a slower but more accuarte
    !  Householder QR is performed.
    !
    IMPLICIT NONE
    !
    ! ... I/O variables
    !
    INTEGER,     INTENT (IN) :: n, k, ldx, ldr
    COMPLEX(DP), INTENT (INOUT) :: X(ldx,*), SX(ldx,*)
    COMPLEX(DP), INTENT(OUT) :: R(ldr,*)
    !
    ! ... local variables
    !
    COMPLEX(DP):: XTX(k,k)
    INTEGER    :: ierr = 0, info = 0
    !
    ! ... do Cholesky QR unless X is rank deficient
    !
    CALL ZGEMM('C','N', k, k, n, C_ONE, X, ldx, SX, ldx, C_ZERO, XTX, k)
    !
    CALL mp_sum( XTX, intra_bgrp_comm )
    !
    CALL ZPOTRF('U', k, XTX, k, info)
    IF ( info == 0 ) THEN
       !
       ! ... triangualar solve
       CALL ZTRSM('R', 'U', 'N', 'N', n, k, C_ONE, XTX, k, X, ldx)
       !
    ELSE
       ! TBD: QR
       WRITE(stdout,*) '[Q, R] = qr(X, 0) failed'
       STOP
!           ![X, R] = qr(X, 0)
!           ! QR without forming Q
!           !
!           call zgeqrf(nn,kk,X,nn,tau,wqr,lwqr,info)
!           !
!           if (info<0) then
!              write(*,*), '[Q, R] = qr(X, 0) failed'
!              stop
!           endif
!           ! QR: constructing Q
!           !
!           call zungqr(nn,kk,kk,X,nn,tau,wqr,lwqr,info)
!           !
    END IF
    !
    ! ... also return R factor if needed
    !
    CALL ZLACPY('U', k, k,  XTX, k,  R, ldr)
    !
    RETURN
    !
  END SUBROUTINE cholQR
  !
  !
  !
  SUBROUTINE cholQR_dmat(kdim, k, X, SX, kdimx,  Rl, desc)
    !
    ! Distributed version of cholQR
    !
    IMPLICIT NONE
    !
    ! ... I/O variables
    !
    INTEGER,     INTENT (IN) :: kdim, k, kdimx
    COMPLEX(DP), INTENT (INOUT) :: X(kdimx,k), SX(kdimx,k)
    TYPE(la_descriptor), INTENT (IN)  :: desc
    COMPLEX(DP), INTENT(OUT) :: Rl(:, :)
    ! inverse of the upper triangular Cholesky factor
    !
    ! ... local variables
    !
    COMPLEX(DP)  :: buffer(kdimx,k)
    COMPLEX(DP), ALLOCATABLE   :: XTXl(:,:)
    INTEGER    :: nx, ierr = 0
#ifdef __SCALAPACK
    INTEGER     :: desc_sca( 16 ), info
#endif
    !
    nx = desc%nrcx
    !
    IF ( la_proc ) THEN
       !
       ALLOCATE( XTXl( nx, nx ), STAT=ierr )
       IF( ierr /= 0 ) &
          CALL errore( 'ppcg ',' cannot allocate XTXl ', ABS(ierr) )
       !
    ELSE
       !
       ALLOCATE( XTXl( 1, 1 ), STAT=ierr )
       IF( ierr /= 0 ) &
          CALL errore( 'ppcg ',' cannot allocate XTXl ', ABS(ierr) )
       !
    END IF
    !
    ! ... Perform Cholesky of X'X
    !
    CALL compute_distmat(XTXl, desc, X, SX, k)
    !
    IF ( la_proc ) THEN
       !
#ifdef __SCALAPACK
       CALL descinit( desc_sca, k, k, nx, nx, 0, 0, ortho_cntx, SIZE( XTXl, 1 ), info )
       IF( info /= 0 ) CALL errore( ' ppcg ', ' descinit ', ABS( info ) )
       !
       CALL PZPOTRF( 'U', k, XTXl, 1, 1, desc_sca, info )
!      IF( info /= 0 ) CALL errore( ' ppcg ', ' problems computing cholesky ', ABS( info ) )
       !
       IF ( info == 0 ) THEN
          !
!          ! set the lower triangular part to zero
!          CALL sqr_zsetmat( 'L', k, C_ZERO, XTXl, size(XTXl,1), desc )
          !
          ! find inverse of the upper triangular Cholesky factor R
          CALL PZTRTRI( 'U', 'N', k, XTXl, 1, 1, desc_sca, info )
          IF( info /= 0 ) CALL errore( ' ppcg ', ' problems computing inverse ', ABS( info ) )
          !
          ! set the lower triangular part to zero
          CALL sqr_zsetmat( 'L', k, C_ZERO, XTXl, size(XTXl,1), desc )
       !
       ELSE
       ! TBD: QR
          WRITE(stdout,*) '[Q, R] = qr(X, 0) failed'
          STOP
!           ![X, R] = qr(X, 0)
!           ! QR without forming Q
!           !
!           call zgeqrf(nn,kk,X,nn,tau,wqr,lwqr,info)
!           !
!           if (info<0) then
!              write(*,*), '[Q, R] = qr(X, 0) failed'
!              stop
!           endif
!           ! QR: constructing Q
!           !
!           call zungqr(nn,kk,kk,X,nn,tau,wqr,lwqr,info)
!           !
       END IF
#else
       CALL qe_pzpotrf( XTXl, nx, k, desc )
       !
       CALL qe_pztrtri ( XTXl, nx, k, desc )
#endif
    !
    !
    END IF
    !
    CALL zgemm_dmat( kdim, k, kdimx, desc, C_ONE, X, XTXl, C_ZERO, buffer )
    !
    X = buffer
    ! ... also return R factor
    Rl = XTXl
    !
    DEALLOCATE(XTXl)
    !
    RETURN
    !
  END SUBROUTINE cholQR_dmat
  !
  !
  !
  SUBROUTINE lock_epairs(kdim, nbnd, btype, w, kdimx, tol, nact, act_idx)
     !
     ! Lock converged eigenpairs: detect "active" columns of w
     ! by checking if each column has norm greater than tol.
     ! Returns the number of active columns and their indices.
     !
     IMPLICIT NONE
     !
     ! ... I/O variables
     !
     INTEGER,     INTENT (IN) :: kdim, kdimx, nbnd, btype(nbnd)
     COMPLEX(DP), INTENT (IN) :: w(kdimx,nbnd)
     REAL(DP),    INTENT (IN) :: tol
     INTEGER,     INTENT(OUT) :: nact, act_idx(nbnd)
     !
     ! ... local variables
     !
     INTEGER         :: j
     REAL(DP)        :: rnrm_store(nbnd), band_tollerance 
     REAL(DP), EXTERNAL :: DDOT

     !
     nact = 0
     ! ... Compute norms of each column of psi
     rnrm_store = 0.D0
     CALL divide(inter_bgrp_comm,nbnd,n_start,n_end); my_n = n_end - n_start + 1; !write (*,*) nbnd,n_start,n_end
     DO j = n_start, n_end
        !
        rnrm_store(j)   =  DDOT(2*kdim, w(:,j), 1, w(:,j), 1)
        !
     END DO
     CALL mp_sum( rnrm_store, inter_bgrp_comm )
     !
     CALL mp_sum( rnrm_store, intra_bgrp_comm )
     !
     DO j = 1, nbnd
        !
        if ( btype(j) == 0 ) then
             band_tollerance = max(2.5*tol,1.d-3)
        else
             band_tollerance = tol
        end if
        !
        rnrm_store(j) = SQRT( rnrm_store(j) )
        !
        IF ( (print_info >= 2) .AND. (iter > 1) )  THEN
          write(stdout, '( "Eigenvalue ", I5, " = ", 1pe12.4, ". Residual norm = ",  1pe9.2)') &
                      j, e(j), rnrm_store(j)
        END IF
        !
        IF ( rnrm_store(j) > band_tollerance ) THEN
           nact = nact + 1
           act_idx(nact) = j
        END IF
        !
     END DO
     !
     !
  END SUBROUTINE lock_epairs
  !
  !
  SUBROUTINE setup_param
     !
     ! Based on the information about active columns of psi,
     ! set up iteration parameters, such as:
     !   - number of subblock
     !   - size of the last sub-block
     !   - tolerance level for the trace difference of reduced hamiltonians
     !     on consecutive iterations
     !   - replicate or re-distribute the Gram matrix G
     !
     IMPLICIT NONE
     !
     ! ... I/O variables
     !
!     INTEGER,     INTENT (IN) :: nbnd, nact, act_idx(nbnd), sbsize
!     INTEGER,     INTENT(OUT) :: nsb, sbsize_last
!     REAL(DP),    INTENT(OUT) :: tol
     !
     ! ... local variables
     !
     INTEGER         :: ierr
     INTEGER :: nx, act_thresh
     ! maximum local block dimension
     ! Threshold on the number of active vectors that determines if
     ! Choleski is distributed or replicated
     !
     !
     ! ... Compute the number of sub-blocks and the size of the last sub-block
     sbsize_last = sbsize
     nsb   =  FLOOR( DBLE(nact)/DBLE(sbsize) )
     IF ( MOD( nact, sbsize ) /= 0) THEN
        sbsize_last =  nact - sbsize*nsb
        nsb         =  nsb + 1
     END IF
     !
     ! ... Compute the current tolerance level for nact pairs
     trtol   = ethr*SQRT(DBLE(nact)) ! MIN( ethr*SQRT(DBLE(nact)), 1.0D-2 )
     !
     ! ... Redistribute matrices for Cholesky because number of active vectors has changed
     !
     ! Don't want to run distributed Cholesky for matrices of size <= 100
     act_thresh = MAX( 100, np_ortho(1) )
     !
     IF ( ( nact > act_thresh ) .AND. (nact /= nact_old) ) THEN
        !
        IF ( ALLOCATED(Gl) ) DEALLOCATE(Gl)
        !
        CALL desc_init( nact, desc, irc_ip, nrc_ip  )
        !
        nx = desc%nrcx
        !
        IF ( la_proc ) THEN
           !
           ALLOCATE( Gl( nx, nx ), STAT=ierr )
           IF( ierr /= 0 ) &
              CALL errore( 'ppcg ',' cannot allocate Gl ', ABS(ierr) )
           !
        ELSE
           !
           ALLOCATE( Gl( 1, 1 ), STAT=ierr )
           IF( ierr /= 0 ) &
           CALL errore( 'ppcg ',' cannot allocate Gl ', ABS(ierr) )
           !
        END IF
        !
        force_repmat = .FALSE.
        !
     ELSE
       ! replicate Cholesky because number of active vectors is small
       !
       IF (nact <= act_thresh) THEN
          !
          force_repmat = .TRUE.
          !
          IF ( ALLOCATED(Gl) ) DEALLOCATE(Gl)
          !
        ELSE
          !
          force_repmat = .FALSE.
          !
        END IF
           !
     END IF
     !
     !
     ! ... Allocate storage for the small matrices in separate minimization loop
     IF ( ALLOCATED(K_store) ) DEALLOCATE(K_store)
     IF ( ALLOCATED(M_store) ) DEALLOCATE(M_store)
     ALLOCATE ( K_store(sbsize3, sbsize3*nsb), M_store(sbsize3,sbsize3*nsb), stat = ierr )
     IF (ierr /= 0) &
        CALL errore( 'ppcg ',' cannot allocate K_store and M_store ', ABS(ierr) )
     !
     RETURN
     !
  END SUBROUTINE setup_param
  !
  !
  !
  !
  SUBROUTINE extract_epairs_dmat(kdim, nbnd, kdimx, e, psi, hpsi, spsi)
     !
     ! Perform the Rayleigh-Ritz to "rotate" psi to eigenvectors
     ! The Gram matrices are distributed over la processor group
     !
     !
     IMPLICIT NONE
     !
     ! ... I/O variables
     !
     INTEGER,     INTENT (IN)  :: kdim, kdimx, nbnd
     REAL(DP),    INTENT (OUT) :: e(nbnd)
     COMPLEX(DP), INTENT (INOUT)  :: psi(kdimx,nbnd), hpsi(kdimx,nbnd), spsi(kdimx,nbnd)
!     LOGICAL,     INTENT(IN), OPTIONAL :: ortho
     ! ortho = .true.(default) then orthogonalization of psi prior to the RR is enforced
     !
     ! ... local variables
     !
     COMPLEX (DP), ALLOCATABLE    :: Hl(:,:), Sl(:,:)
     ! local part of projected Hamiltonian and of the overlap matrix
     TYPE(la_descriptor)       :: desc
     !
     ! Matrix distribution descriptors to temporary store the "global" current descriptor
     LOGICAL :: la_proc_store
     ! flag to distinguish procs involved in linear algebra
     INTEGER, ALLOCATABLE :: irc_ip_store( : )
     INTEGER, ALLOCATABLE :: nrc_ip_store( : )
     INTEGER, ALLOCATABLE :: rank_ip_store( :, : )
     !
     COMPLEX(DP), ALLOCATABLE  :: psi_t(:, :), hpsi_t(:, :), spsi_t(:, :)
     COMPLEX(DP), ALLOCATABLE  :: vl(:,:)
     COMPLEX(DP)               :: buffer(kdimx,nbnd)
!     REAL(DP)                  :: R(nbnd, nbnd)
     INTEGER                   :: nx
!     LOGICAL                   :: do_orth
     !
     ALLOCATE ( psi_t(kdimx,nbnd), hpsi_t(kdimx,nbnd), stat = ierr )
     IF (ierr /= 0) CALL errore( 'ppcg ',' cannot allocate psi_t and hpsi_t ', ABS(ierr) )
     if (overlap) ALLOCATE ( spsi_t(kdimx,nbnd), stat = ierr )
     IF (ierr /= 0) CALL errore( 'ppcg ',' cannot allocate spsi_t ', ABS(ierr) )
     !
     ! Store current distributed matrix descriptor information
     ALLOCATE ( irc_ip_store( np_ortho(1)  ), stat = ierr )
     IF (ierr /= 0) CALL errore( 'ppcg ',' cannot allocate irc_ip_store ', ABS(ierr) )
     !
     ALLOCATE ( nrc_ip_store( np_ortho(1) ), stat = ierr )
     IF (ierr /= 0) CALL errore( 'ppcg ',' cannot allocate nrc_ip_store ', ABS(ierr) )
     !
     ALLOCATE ( rank_ip_store( np_ortho(1), np_ortho(2) ), stat = ierr )
     IF (ierr /= 0) CALL errore( 'ppcg ',' cannot allocate rank_ip_store ', ABS(ierr) )
     !
     irc_ip_store  = irc_ip
     nrc_ip_store  = nrc_ip
     rank_ip_store = rank_ip
     !
     CALL desc_init( nbnd, desc, irc_ip, nrc_ip  )
     !
     nx = desc%nrcx
     !
     IF ( la_proc ) THEN
        !
        ALLOCATE( vl( nx, nx ), STAT=ierr )
        IF( ierr /= 0 ) &
           CALL errore( 'ppcg ',' cannot allocate vl ', ABS(ierr) )
        !
        ALLOCATE( Sl( nx, nx ), STAT=ierr )
        IF( ierr /= 0 ) &
          CALL errore( 'ppcg ',' cannot allocate Sl ', ABS(ierr) )
        !
        ALLOCATE( Hl( nx, nx ), STAT=ierr )
        IF( ierr /= 0 ) &
          CALL errore( 'ppcg ',' cannot allocate Hl ', ABS(ierr) )
        !
     ELSE
        !
        ALLOCATE( vl( 1, 1 ), STAT=ierr )
        IF( ierr /= 0 ) &
           CALL errore( 'pregterg ',' cannot allocate vl ', ABS(ierr) )
        !
        ALLOCATE( Sl( 1, 1 ), STAT=ierr )
        IF( ierr /= 0 ) &
           CALL errore( 'ppcg ',' cannot allocate Sl ', ABS(ierr) )
        !
        ALLOCATE( Hl( 1, 1 ), STAT=ierr )
        IF( ierr /= 0 ) &
           CALL errore( 'ppcg ',' cannot allocate Hl ', ABS(ierr) )
        !
     END IF
     !
     !
     !  G = psi'*hpsi
     CALL compute_distmat(Hl, desc, psi, hpsi, nbnd)
     if (overlap) then
        CALL compute_distmat(Sl, desc, psi, spsi, nbnd)
     else
        CALL compute_distmat(Sl, desc, psi,  psi, nbnd)
     end if
     !
     ! ... diagonalize the reduced hamiltonian
     !     Calling block parallel algorithm
     !
     IF ( do_distr_diag_inside_bgrp ) THEN ! NB on output of pcdiaghg e and vl are the same across ortho_parent_comm
        ! only the first bgrp performs the diagonalization
        IF( my_bgrp_id == root_bgrp_id) CALL pcdiaghg( nbnd, Hl, Sl, nx, e, vl, desc )
        IF( nbgrp > 1 ) THEN ! results must be brodcast to the other band groups
          CALL mp_bcast( vl, root_bgrp_id, inter_bgrp_comm )
          CALL mp_bcast( e,  root_bgrp_id, inter_bgrp_comm )
        ENDIF
     ELSE
        CALL pcdiaghg( nbnd, Hl, Sl, nx, e, vl, desc )
     END IF
     !
     ! "Rotate" psi to eigenvectors
     !
     CALL zgemm_dmat( kdim, nbnd, kdimx, desc, C_ONE,  psi, vl, C_ZERO,  psi_t )
     CALL zgemm_dmat( kdim, nbnd, kdimx, desc, C_ONE, hpsi, vl, C_ZERO, hpsi_t )
     if (overlap) CALL zgemm_dmat( kdim, nbnd, kdimx, desc, C_ONE, spsi, vl, C_ZERO, spsi_t )
     !
     psi   = psi_t ; hpsi  = hpsi_t ;  if (overlap) spsi  = spsi_t
     !
     ! Restore current "global" distributed  matrix descriptor
     !
     irc_ip  = irc_ip_store
     nrc_ip  = nrc_ip_store
     rank_ip = rank_ip_store
     !
     DEALLOCATE ( irc_ip_store, nrc_ip_store, rank_ip_store )
     DEALLOCATE(psi_t, hpsi_t) ; if (overlap) DEALLOCATE(spsi_t)
     DEALLOCATE(Hl, Sl)
     DEALLOCATE(vl)
     !
  END SUBROUTINE extract_epairs_dmat
  !
  !
  !
  REAL(DP) FUNCTION get_trace(G, ld, k)
     !
     !  This function returns trace of a k-by-k matrix G
     !
     IMPLICIT NONE
     !
     ! ... I/O variables
     !
     INTEGER,     INTENT (IN) :: ld, k
     COMPLEX(DP), INTENT (IN) :: G(ld,*)
     !
     ! ... local variables
     !
     INTEGER    :: j
     !
     get_trace = 0.D0
     DO j = 1, k
        get_trace = get_trace + DBLE(G(j, j))
     END DO
     !
     RETURN
     !
  END FUNCTION get_trace
  !
  !
  !
  SUBROUTINE deallocate_all
    !
    ! This subroutine releases the allocated memory
    !
    IF ( ALLOCATED(hpsi) )    DEALLOCATE ( hpsi )
    IF ( ALLOCATED(spsi) )    DEALLOCATE ( spsi )
    IF ( ALLOCATED(w) )       DEALLOCATE ( w )
    IF ( ALLOCATED(hw) )      DEALLOCATE ( hw )
    IF ( ALLOCATED(sw) )      DEALLOCATE ( sw )
    IF ( ALLOCATED(p) )       DEALLOCATE ( p )
    IF ( ALLOCATED(hp) )      DEALLOCATE ( hp )
    IF ( ALLOCATED(sp) )      DEALLOCATE ( sp )
    IF ( ALLOCATED(K) )       DEALLOCATE ( K )
    IF ( ALLOCATED(M) )       DEALLOCATE ( M )
    IF ( ALLOCATED(K_store) ) DEALLOCATE ( K_store )
    IF ( ALLOCATED(M_store) ) DEALLOCATE ( M_store )
    IF ( ALLOCATED(cwork) )   DEALLOCATE ( cwork )
    IF ( ALLOCATED(rwork) )   DEALLOCATE ( rwork )
    IF ( ALLOCATED(iwork) )   DEALLOCATE ( iwork )
    IF ( ALLOCATED(irc_ip) )  DEALLOCATE ( irc_ip )
    IF ( ALLOCATED(nrc_ip) )  DEALLOCATE ( nrc_ip )
    IF ( ALLOCATED(rank_ip) ) DEALLOCATE ( rank_ip )
    IF ( ALLOCATED(Gl) )      DEALLOCATE ( Gl )
!    IF ( ALLOCATED(G) )       DEALLOCATE( G )
!    IF ( ALLOCATED(Sl) )      DEALLOCATE( Sl )
    !
  END SUBROUTINE deallocate_all
  !
  !
! dmat begin
  SUBROUTINE desc_init( nsiz, desc, irc_ip, nrc_ip )
!    copy-paste from pregterg
     !
     INTEGER, INTENT(IN)  :: nsiz
     TYPE(la_descriptor), INTENT(OUT) :: desc
     INTEGER, INTENT(OUT) :: irc_ip(:)
     INTEGER, INTENT(OUT) :: nrc_ip(:)

     INTEGER :: i, j, rank
     !
     CALL descla_init( desc, nsiz, nsiz, np_ortho, me_ortho, ortho_comm, ortho_cntx, ortho_comm_id)
!     !
!     nx = desc%nrcx   ! nx should be initialized outside whenever needed
     !
     DO j = 0, desc%npc - 1
        CALL descla_local_dims( irc_ip( j + 1 ), nrc_ip( j + 1 ), desc%n, desc%nx, np_ortho(1), j )
        DO i = 0, desc%npr - 1
           CALL GRID2D_RANK( 'R', desc%npr, desc%npc, i, j, rank )
           rank_ip( i+1, j+1 ) = rank * leg_ortho
        END DO
     END DO
     !
     la_proc = .FALSE.
     IF( desc%active_node > 0 ) la_proc = .TRUE.
     !
     RETURN
  END SUBROUTINE desc_init
  !
  !
  SUBROUTINE compute_distmat( dm, desc, v, w, k)
!    Copy-paste from pcegterg and desc added as a parameter
     !
     !  This subroutine compute <vi|wj> and store the
     !  result in distributed matrix dm
     !
     IMPLICIT NONE
     !
     ! ... I/O variables
     !
     COMPLEX(DP), INTENT(OUT)   :: dm( :, : )
     TYPE(la_descriptor), INTENT(IN) :: desc
     COMPLEX(DP), INTENT(IN) :: v(:,:), w(:,:)
     INTEGER, INTENT(IN)     :: k
     ! global size of dm = number of vectors in v and w blocks
     !
     ! ... local variables
     !
     INTEGER :: ipc, ipr
     INTEGER :: nr, nc, ir, ic, root
     COMPLEX(DP), ALLOCATABLE :: work( :, : )
     INTEGER :: nx
     ! maximum local block dimension
     !
     nx = desc%nrcx
     !
     ALLOCATE( work( nx, nx ) )
     !
     work = C_ZERO
     !
     DO ipc = 1, desc%npc !  loop on column procs
        !
        nc = nrc_ip( ipc )
        ic = irc_ip( ipc )
        !
        DO ipr = 1, ipc ! use symmetry for the loop on row procs
           !
           nr = nrc_ip( ipr )
           ir = irc_ip( ipr )
           !
           !  rank of the processor for which this block (ipr,ipc) is destinated
           !
           root = rank_ip( ipr, ipc )

           ! use blas subs. on the matrix block

           CALL ZGEMM( 'C','N', nr, nc, kdim, C_ONE, v(1,ir), kdimx, w(1,ic), kdimx, C_ZERO, work, nx )

           ! accumulate result on dm of root proc.

           CALL mp_root_sum( work, dm, root, ortho_parent_comm )

        END DO
        !
     END DO
     !
     if (ortho_parent_comm.ne.intra_bgrp_comm .and. nbgrp > 1) dm = dm/nbgrp
     !
!     CALL zsqmher( nbnd, dm, nx, desc )
     CALL zsqmher( k, dm, nx, desc )
     !
     DEALLOCATE( work )
     !
     RETURN
  END SUBROUTINE compute_distmat
  !
  !
  SUBROUTINE zgemm_dmat( n, k, ld, desc, alpha, X, Gl, beta, Y  )
!    Copy-paste from refresh_evc in pregterg with some modifications
     !
     ! Compute Y = alpha*(X*G) + beta*Y, where G is distributed across la processor group
     ! and is given by local block Gl.
     !
     IMPLICIT NONE
     !
     ! ... I/O variables
     !
     INTEGER, INTENT(IN) :: n, k, ld
     ! number of rows of X and Y
     ! number of columns of X,Y  and size/leading dimension of (global) G
     ! leading dimension of X and Y
     TYPE(la_descriptor), INTENT(IN) :: desc
     ! descriptor of G
     COMPLEX(DP), INTENT(IN)      ::  alpha, beta
     COMPLEX(DP), INTENT (IN)     ::  X(ld, k)
     COMPLEX(DP), INTENT (INOUT)  ::  Y(ld, k)
     COMPLEX(DP), INTENT(IN)      ::  Gl( :, :)
     !
     ! ... local variables
     !
     INTEGER :: ipc, ipr
     INTEGER :: nr, nc, ir, ic, root
     COMPLEX(DP), ALLOCATABLE :: Gltmp( :, : )
     COMPLEX(DP), ALLOCATABLE :: Xtmp( :, : )
     COMPLEX(DP) :: gamm
     INTEGER :: nx
     !
     nx = desc%nrcx
     !
     ALLOCATE( Gltmp( nx, nx ) )
     ALLOCATE( Xtmp( ld, k ) )
     !
     DO ipc = 1, desc%npc
        !
        nc = nrc_ip( ipc )
        ic = irc_ip( ipc )
        !
        IF( ic <= k ) THEN
           !
           nc = min( nc, k - ic + 1 )
           !
           gamm = C_ZERO

           DO ipr = 1, desc%npr
              !
              nr = nrc_ip( ipr )
              ir = irc_ip( ipr )
              !
              root = rank_ip( ipr, ipc )
              IF( ipr-1 == desc%myr .AND. ipc-1 == desc%myc .AND. la_proc ) THEN
                 !
                 !  this proc sends his block
                 !
                 CALL mp_bcast( Gl(:,1:nc), root, ortho_parent_comm )
                 CALL ZGEMM( 'N','N', n, nc, nr, C_ONE, X(1,ir), ld, Gl, nx, gamm, Xtmp(1,ic), ld )
              ELSE
                 !
                 !  all other procs receive
                 !
                 CALL mp_bcast( Gltmp(:,1:nc), root, ortho_parent_comm )
                 CALL ZGEMM( 'N','N', n, nc, nr, C_ONE, X(1,ir), ld, Gltmp, nx, gamm, Xtmp(1,ic), ld )
              END IF
              !
              gamm = C_ONE
              !
           END DO
           !
        END IF
        !
     END DO
     !
     IF (beta /= 0.D0) THEN
        Y = alpha*Xtmp + beta*Y
     ELSE
       Y = alpha*Xtmp
     END IF
     !
     DEALLOCATE( Gltmp )
     DEALLOCATE( Xtmp )
     !
     RETURN
     !
  END SUBROUTINE zgemm_dmat
  !

! dmat end

  SUBROUTINE print_rnrm
     !
     !  Compute the subspce residual
     !
     COMPLEX(DP), ALLOCATABLE :: psi_t(:,:), hpsi_t(:,:), res(:,:)
     COMPLEX(DP), ALLOCATABLE :: G(:,:), work(:)
     REAL(DP) :: rnrm
     REAL(DP), EXTERNAL       :: ZLANGE
!ev new begin
integer   :: nwanted, nguard
nguard = 0 ! 24 ! 50
!ev new  end

     rnrm  = 0.D0
     ALLOCATE(psi_t(kdimx,nbnd), hpsi_t(kdimx,nbnd), res(kdimx,nbnd))
     ALLOCATE(G(nbnd, nbnd))
     !
!!! EV begin comment
!     psi_t(:,1:nbnd)  = psi(:,1:nbnd)
!     hpsi_t(:,1:nbnd) = hpsi(:,1:nbnd)
!     !
!     !     G = psi_t'hpsi_t
!     CALL ZGEMM('C','N', nbnd, nbnd, npw, C_ONE, psi_t, npwx, hpsi_t, npwx, C_ZERO, G, nbnd)
!     CALL mp_sum( G, intra_bgrp_comm )
!     !
!     !    res = hpsi_t - psi*G
!     res = hpsi_t;
!     CALL ZGEMM('N','N',npw, nbnd, nbnd, -C_ONE, psi_t, npwx, G, nbnd, C_ONE, res, npwx)
!     !
!     ! ... get the Frobenius norm of the residual
!     rnrm = ZLANGE('F', npw, nbnd, res, npwx, work)  ! work is not referenced for Frobenius norm
!     rnrm = abs(rnrm)**2
!!! EV end comment
!!!
!!!
!  EV begin new
!   Instead computing the norm of the whole residual, compute it for X(:,nbnd-nguard)
     nwanted = nbnd - nguard
     psi_t(:,1:nwanted)  = psi(:,1:nwanted)
     hpsi_t(:,1:nwanted) = hpsi(:,1:nwanted)
     !
     !     G = psi_t'hpsi_t
     CALL ZGEMM('C','N', nwanted, nwanted, kdim, C_ONE, psi_t, kdimx, hpsi_t, kdimx, C_ZERO, G, nbnd)
     CALL mp_sum( G, intra_bgrp_comm )
     !
     !    res = hpsi_t - psi*G
     res(:,1:nwanted) = hpsi_t(:,1:nwanted);
     CALL ZGEMM('N','N',kdim, nwanted, nwanted, -C_ONE, psi_t, kdimx, G, nbnd, C_ONE, res, kdimx)
     !
     ! ... get the Frobenius norm of the residual
     rnrm = ZLANGE('F', kdim, nwanted, res, kdimx, work)  ! work is not referenced for Frobenius norm
     rnrm = rnrm**2
!   EV end new

     !
     CALL mp_sum( rnrm, intra_bgrp_comm )
     !
     rnrm = SQRT(ABS(rnrm))
     res_array(iter) = rnrm
     !
     DEALLOCATE(psi_t, hpsi_t, res)
     DEALLOCATE(G)
     !

!     !CALL flush( stdout )
     !
  END SUBROUTINE print_rnrm

!- routines to perform threaded assignements

  SUBROUTINE threaded_assign(array_out, array_in, kdimx, nact, act_idx, bgrp_root_only)
  ! 
  !  assign (copy) a complex array in a threaded way
  !
  !  array_out( 1:kdimx, 1:nact ) = array_in( 1:kdimx, 1:nact )       or
  !
  !  array_out( 1:kdimx, 1:nact ) = array_in( 1:kdimx, act_idx(1:nact) )
  !
  !  if the index array act_idx is given
  !
  !  if  bgrp_root_only is present and .true. the assignement is made only by the 
  !  MPI root process of the bgrp and array_out is zeroed otherwise
  !
  USE util_param,   ONLY : DP
  !
  IMPLICIT NONE
  !
  COMPLEX(DP), INTENT(OUT) :: array_out( kdimx, nact )
  COMPLEX(DP), INTENT(IN)  :: array_in ( kdimx, * )
  INTEGER, INTENT(IN)      :: kdimx, nact
  INTEGER, INTENT(IN), OPTIONAL :: act_idx( * )
  LOGICAL, INTENT(IN), OPTIONAL :: bgrp_root_only
  !
  INTEGER, PARAMETER :: blocksz = 256
  INTEGER :: nblock

  INTEGER :: i, j
  !
  IF (kdimx <=0 .OR. nact<= 0) RETURN
  !
  IF (present(bgrp_root_only) ) THEN
     IF (bgrp_root_only .AND. ( my_bgrp_id /= root_bgrp_id ) ) THEN
        call threaded_memset( array_out, 0.d0, 2*kdimx*nact )
        RETURN
     END IF
  END IF

  nblock = (kdimx - 1)/blocksz  + 1

  IF (present(act_idx) ) THEN
     !$omp parallel do collapse(2)
     DO i=1, nact
        DO j=1, nblock
           array_out(1+(j-1)*blocksz:MIN(j*blocksz,kdimx), i ) = array_in(1+(j-1)*blocksz:MIN(j*blocksz,kdimx), act_idx( i ) ) 
        ENDDO
     ENDDO
     !$omp end parallel do
  ELSE
     !$omp parallel do collapse(2)
     DO i=1, nact
        DO j=1, nblock
           array_out(1+(j-1)*blocksz:MIN(j*blocksz,kdimx), i ) = array_in(1+(j-1)*blocksz:MIN(j*blocksz,kdimx), i ) 
        ENDDO
     ENDDO
     !$omp end parallel do
  END IF
  !
  END SUBROUTINE threaded_assign

END SUBROUTINE ppcg_k
