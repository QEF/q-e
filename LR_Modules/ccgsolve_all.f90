!
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------
subroutine ccgsolve_all (ch_psi, ccg_psi, e, d0psi, dpsi, h_diag, &
     ndmx, ndim, ethr, ik, kter, conv_root, anorm, nbnd, npol, freq_c)
  !----------------------------------------------------------------------
  !
  !     iterative solution of the linear systems (i=1,nbnd):
  !
  !                 ( H - e_i + Q ) * dpsi_i = d0psi_i                      (1)
  !
  !     where H is a complex hermitean matrix, e_i is a real scalar, Q is a
  !     projector on occupied states, dpsi_i and d0psi_ are complex vectors
  !
  !     on input:
  !                 ch_psi   EXTERNAL  name of a subroutine:
  !                          Calculates  (H-e+Q)*psi products.
  !                          Vectors psi and psip should be dimensioned
  !                          (ndmx,nbnd)
  !
  !                 cg_psi   EXTERNAL  name of a subroutine:
  !                          which calculates (h-e)^-1 * psi, with
  !                          some approximation, e.g. (diag(h)-e)
  !
  !                 e        real     unperturbed eigenvalue.
  !
  !                 dpsi     contains an estimate of the solution
  !                          vector.
  !
  !                 d0psi    contains the right hand side vector
  !                          of the system.
  !
  !                 ndmx     integer row dimension of dpsi, ecc.
  !
  !                 ndim     integer actual row dimension of dpsi
  !
  !                 ethr     real     convergence threshold. solution
  !                          improvement is stopped when the error in
  !                          eq (1), defined as l.h.s. - r.h.s., becomes
  !                          less than ethr in norm.
  !
  !     on output:  dpsi     contains the refined estimate of the
  !                          solution vector.
  !
  !                 d0psi    is corrupted on exit
  !
  !   revised (extensively)       6 Apr 1997 by A. Dal Corso & F. Mauri
  !   revised (to reduce memory) 29 May 2004 by S. de Gironcoli
  !
  USE kinds,          ONLY : DP
  USE mp_bands,       ONLY : intra_bgrp_comm, inter_bgrp_comm, use_bgrp_in_hpsi
  USE mp,             ONLY : mp_sum, mp_barrier
  USE control_flags,  ONLY : gamma_only
  USE gvect,          ONLY : gstart
  USE eqv,            ONLY : evq
#if defined(__CUDA)
 USE cublas
#endif  

  implicit none
  !
  !   first the I/O variables
  !
  integer :: ndmx, & ! input: the maximum dimension of the vectors
             ndim, & ! input: the actual dimension of the vectors
             kter, & ! output: counter on iterations
             nbnd, & ! input: the number of bands
             npol, & ! input: number of components of the wavefunctions
             ik      ! input: the k point

  real(DP) :: &
             e(nbnd), & ! input: the actual eigenvalue
             anorm,   & ! output: the norm of the error in the solution
             ethr       ! input: the required precision

  complex(DP) :: &
             dpsi (ndmx*npol, nbnd), & ! output: the solution of the linear syst
             h_diag(ndmx*npol,nbnd), & ! input: an estimate of ( H - \epsilon -w)
                                       ! w is complex
             freq_c                , & ! complex frequency
             d0psi (ndmx*npol, nbnd)   ! input: the known term

  logical :: conv_root ! output: if true the root is converged
  external ch_psi      ! input: the routine computing ch_psi
  external ccg_psi      ! input: the routine computing cg_psi
  !
  !  here the local variables
  !
  integer, parameter :: maxter = 1000
  ! the maximum number of iterations
  integer :: iter, ibnd, ibnd_, lbnd
  ! counters on iteration, bands
  integer , allocatable :: conv (:)
  ! if 1 the root is converged

  complex(DP), allocatable :: g (:,:),  t (:,:),  h (:,:),  hold (:,:)
  COMPLEX(DP), ALLOCATABLE :: gs (:,:), ts (:,:), hs (:,:), hsold (:,:)
  !  the gradient of psi
  !  the preconditioned gradient
  !  the delta gradient
  !  the conjugate gradient
  !  work space
  complex(DP) ::  dcgamma, dcgamma1, dclambda, dclambda1
  !  the ratio between rho
  !  step length
!  REAL(kind=dp), EXTERNAL :: ddot
  REAL(kind=dp), EXTERNAL :: myddot
  REAL(kind=dp), EXTERNAL :: myddotv2
  !  the scalar product
  real(DP), allocatable    :: eu (:)
  complex(DP), allocatable :: rho (:), rhoold (:), euc (:), a(:), c(:)
  ! the residue
  ! auxiliary for h_diag
  real(DP) :: kter_eff
  ! account the number of iterations with b
  ! coefficient of quadratic form
  integer, allocatable :: indb(:)
  
  ! bgrp parallelization auxiliary variables
  INTEGER :: n_start, n_end, my_nbnd
  logical :: lsave_use_bgrp_in_hpsi
  !auxiliary for ddot
  real(DP) :: ddotval, addot, cddot
  !
  call start_clock ('ccgsolve')

  call divide (inter_bgrp_comm,nbnd,n_start,n_end)
  my_nbnd = n_end - n_start + 1

  ! allocate workspace (bgrp distributed)
  allocate ( conv(nbnd) )
  allocate ( g(ndmx*npol,my_nbnd), t(ndmx*npol,my_nbnd), h(ndmx*npol,my_nbnd), &
             hold(ndmx*npol,my_nbnd) )
  allocate ( gs(ndmx*npol,my_nbnd), ts(ndmx*npol,my_nbnd), hs(ndmx*npol,my_nbnd), &
             hsold(ndmx*npol,my_nbnd) )
  allocate ( a(my_nbnd), c(my_nbnd) )
  allocate ( rho(my_nbnd), rhoold(my_nbnd) )
  allocate ( eu(my_nbnd) )
  allocate ( euc(my_nbnd) )
  allocate ( indb(my_nbnd) )
  !      WRITE( stdout,*) g,t,h,hold

  kter_eff = 0.d0 ; conv (1:nbnd) = 0

  !!g=(0.d0,0.d0); 
!  t=(0.d0,0.d0); h=(0.d0,0.d0); hold=(0.d0,0.d0)
!  gs=(0.d0,0.d0); ts=(0.d0,0.d0); hs=(0.d0,0.d0); hsold=(0.d0,0.d0)

  rho= (0.0d0, 0.0d0); rhoold=(0.0d0,0.d0)

  ! bgrp parallelization is done outside h_psi/s_psi. set use_bgrp_in_hpsi temporarily to false
  lsave_use_bgrp_in_hpsi = use_bgrp_in_hpsi ; use_bgrp_in_hpsi = .false.
  !$acc enter data create(a(1:my_nbnd),c(1:my_nbnd),eu(1:my_nbnd),euc(1:my_nbnd),t(1:ndmx*npol,1:my_nbnd),ts(1:ndmx*npol,1:my_nbnd),g(1:ndmx*npol,1:my_nbnd),gs(1:ndmx*npol,1:my_nbnd),h(1:ndmx*npol,1:my_nbnd),hs(1:ndmx*npol,1:my_nbnd),hold(1:ndmx*npol,1:my_nbnd),hsold(1:ndmx*npol,1:my_nbnd)) copyin(e(1:nbnd),dpsi(1:ndmx*npol,1:nbnd),evq,h_diag(1:ndmx*npol,1:nbnd),d0psi(1:ndmx*npol,1:nbnd))
  !$acc kernels present(g,t,h,hold,gs,ts,hs,hsold)
  g=(0.d0,0.d0)
  t=(0.d0,0.d0)
  h=(0.d0,0.d0)
  hold=(0.d0,0.d0)
  gs=(0.d0,0.d0)
  ts=(0.d0,0.d0)
  hs=(0.d0,0.d0)
  hsold=(0.d0,0.d0)  
  !$acc end kernels

!!!!  !$acc kernels present(indb)
  do ibnd = 1, nbnd
     indb(ibnd) = ibnd
  enddo
!!!!  !$acc end kernels
  !$acc kernels present(eu,e)
  eu(:) = e(:)
  !$acc end kernels
  do iter = 1, maxter
     !
     !    compute the gradient. can reuse information from previous step
     !
     if (iter == 1) then
!!!        !$acc kernels present(euc,e,indb)
        do ibnd = n_start, n_end
           euc(ibnd) = CMPLX(e(indb(ibnd))+DREAL(freq_c), DIMAG(freq_c))
        ENDDO
!!!        !$acc end kernels
!$acc update device(euc)

        !$acc data present(dpsi,euc,g)
        call ch_psi (ndim, dpsi, g, euc, ik, my_nbnd)
        !$acc end data
        do ibnd = n_start, n_end ; ibnd_ = ibnd - n_start + 1
           !$acc data present(g,d0psi)
           !$acc host_data use_device(d0psi,g)
           call zaxpy (ndim, (-1.d0,0.d0), d0psi(1,ibnd), 1, g(1,ibnd_), 1)
           !$acc end host_data
           !$acc end data
        enddo
        IF (npol==2) THEN
           do ibnd = n_start, n_end ; ibnd_ = ibnd - n_start + 1
              !$acc data present(g,d0psi)
              !$acc host_data use_device(d0psi,g)
              call zaxpy (ndim, (-1.d0,0.d0), d0psi(ndmx+1,ibnd), 1, g(ndmx+1,ibnd_), 1)
              !$acc end host_data
              !$acc end data
           enddo
        END IF
        !$acc kernels present(g,gs)
        gs(:,:) = CONJG(g(:,:))
        !$acc end kernels
     endif
     !
     !    compute preconditioned residual vector and convergence check
     !
     lbnd = 0
     CALL start_clock('loop1')
     do ibnd = n_start, n_end ;  ibnd_ = ibnd - n_start + 1
        if (conv (ibnd) .eq.0) then
           lbnd = lbnd+1
           !$acc data present(g,h,gs,hs)
           !$acc host_data use_device(g,h,gs,hs)
           call zcopy (ndmx*npol, g (1, ibnd_), 1, h (1, ibnd_), 1)
           call zcopy (ndmx*npol, gs (1, ibnd_), 1, hs (1, ibnd_), 1)
           !$acc end host_data
           !$acc end data

           call ccg_psi(ndmx, ndim, 1, h(1,ibnd_), h_diag(1,ibnd), 1 )
           call ccg_psi(ndmx, ndim, 1, hs(1,ibnd_), h_diag(1,ibnd), -1 )

           IF (gamma_only) THEN
              !$acc data present(g,h)
              !$acc host_data use_device(g,h)
              ddotval=2.0d0*myddot(2*ndmx*npol,h(1,ibnd_),1,g(1,ibnd_),1)
              !$acc end host_data
              !$acc end data
!              rho(lbnd)=2.0d0*ddot(2*ndmx*npol,h(1,ibnd_),1,g(1,ibnd_),1)
              rho(lbnd) = ddotval
              IF(gstart==2) THEN
                 !$acc kernels present(h,g) copy(rho(1:my_nbnd))     
                 rho(lbnd)=rho(lbnd)-DBLE(h(1,ibnd_))*DBLE(g(1,ibnd_))
                 !$acc end kernels
              ENDIF
           ELSE
              !$acc data present(g,hs) 
              !$acc host_data use_device(g,hs) 
              ddotval = myddot (2*ndmx*npol, hs(1,ibnd_), 1, g(1,ibnd_), 1)
              !$acc end host_data
              !$acc end data
              rho(lbnd) = ddotval     
!              rho(lbnd) = dot_product (hs(:,ibnd_), g(:,ibnd_))
           ENDIF

        endif
     enddo
     CALL stop_clock('loop1')
     kter_eff = kter_eff + DBLE (lbnd) / DBLE (nbnd)
     call mp_sum( rho(1:lbnd), intra_bgrp_comm )
     do ibnd = n_end, n_start, -1 ; ibnd_ = ibnd - n_start + 1
        if (conv(ibnd).eq.0) then
           rho(ibnd_)=rho(lbnd)
           lbnd = lbnd -1
           anorm = sqrt ( abs (rho (ibnd_)) )
           if (anorm.lt.ethr) conv (ibnd) = 1
        endif
     enddo
     !
     conv_root = .true.
     do ibnd = n_start, n_end
        conv_root = conv_root.and. (conv (ibnd) .eq.1)
     enddo
     if (conv_root) goto 100
     !
     !        compute the step direction h. Conjugate it to previous step
     !
     lbnd = 0
     CALL start_clock('loop2')
     do ibnd = n_start, n_end ; ibnd_ = ibnd - n_start + 1
        if (conv (ibnd) .eq.0) then
           !
           ! change sign to h and hs
           !
           !$acc kernels present(h,hs)
           h(:,ibnd_)=-1.0d0*h(:,ibnd_)
           hs(:,ibnd_)=-1.0d0*hs(:,ibnd_)
           !$acc end kernels          
!           call dscal (2 * ndmx * npol, - 1.d0, h (1, ibnd_), 1)
!           call dscal (2 * ndmx * npol, - 1.d0, hs (1, ibnd_), 1)
           if (iter.ne.1) then
              dcgamma = rho (ibnd_) / rhoold (ibnd_)
              dcgamma1 = CONJG(dcgamma)
              !$acc data present(h,hold,hs,hsold)
              !$acc host_data use_device(hold,h,hsold,hs)
              call zaxpy (ndmx*npol, dcgamma, hold (1, ibnd_), 1, h (1, ibnd_), 1)
              CALL zaxpy (ndmx*npol, dcgamma1, hsold (1, ibnd_), 1, hs (1, ibnd_), 1)
              !$acc end host_data
              !$acc end data
           endif

           !
           ! here hold is used as auxiliary vector in order to efficiently compute t = A*h
           ! it is later set to the current (becoming old) value of h
           !
           lbnd = lbnd+1
           !$acc data present(h,hold,hs,hsold)
           !$acc host_data use_device(hold,h,hsold,hs)
           call zcopy (ndmx*npol, h (1, ibnd_), 1, hold (1, lbnd), 1)
           CALL zcopy (ndmx*npol, hs (1, ibnd_), 1, hsold (1, lbnd), 1)
           !$acc end host_data
           !$acc end data

!!!           !$acc kernels present(indb)
           indb (lbnd) = ibnd
!!!           !$acc end kernels
        endif
     enddo
     CALL stop_clock('loop2')
     !
     !        compute t = A*h and  ts= A^+ * h 
     !
!!!     !$acc kernels present(indb,euc,e)
     DO ibnd=1,lbnd
        euc(ibnd) = CMPLX(e(indb(ibnd))+DREAL(freq_c), DIMAG(freq_c))
     ENDDO
!!     !$acc end kernels
     !$acc update device(euc)

     !$acc data present(euc, hold, t)
     call ch_psi (ndim, hold, t, euc, ik, lbnd)
     !$acc end data
     !$acc data present(euc, hsold, ts)
     call ch_psi (ndim, hsold, ts, conjg(euc), ik, lbnd)
     !$acc end data
     !
     !        compute the coefficients a and c for the line minimization
     !        compute step length lambda
     !

     lbnd=0
     CALL start_clock('loop3')
     do ibnd = n_start, n_end ; ibnd_ = ibnd - n_start + 1
        if (conv (ibnd) .eq.0) then
           lbnd=lbnd+1
           IF (gamma_only) THEN
              !$acc data present(g,hs,t,a,c) 
              !$acc host_data use_device(g,hs,t) 
              addot = 2.0d0*myddot(2*ndmx*npol,hs(1,ibnd_),1,g(1,ibnd_),1)
              cddot = 2.0d0*myddot(2*ndmx*npol,hs(1,ibnd_),1,t(1,lbnd),1)
              !$acc end host_data
              !$acc kernels present(a,c)
              a(lbnd) = addot
              c(lbnd) = cddot
              !$acc end kernels
              !$acc end data
!              a(lbnd) = 2.0d0*ddot(2*ndmx*npol,hs(1,ibnd_),1,g(1,ibnd_),1)
!              c(lbnd) = 2.0d0*ddot(2*ndmx*npol,hs(1,ibnd_),1,t(1,lbnd),1)
              IF (gstart == 2) THEN
                 !$acc kernels present(a,c,hs,g,t)     
                 a(lbnd)=a(lbnd)-DBLE(hs(1,ibnd_))*DBLE(g(1,ibnd_))
                 c(lbnd)=c(lbnd)-DBLE(hs(1,ibnd_))*DBLE(t(1,lbnd))
                 !$acc end kernels
              ENDIF
           ELSE
              !$acc data present(g,hs,t,a,c) 
              !$acc host_data use_device(g,hs,t)
              addot = myddot (2*ndmx*npol, hs(1,ibnd_), 1, g(1,ibnd_), 1)
              cddot = myddot (2*ndmx*npol, hs(1,ibnd_), 1, t(1,lbnd), 1)
              !$acc end host_data
              !$acc serial present(a,c)
              a(lbnd) = addot
              c(lbnd) = cddot
              !$acc end serial
              !$acc end data
!              a(lbnd) = dot_product (hs(:,ibnd_), g(:,ibnd_))
!              c(lbnd) = dot_product (hs(:,ibnd_), t(:,lbnd))
           ENDIF
        end if
     end do
     CALL stop_clock('loop3')
     !
     !$acc host_data use_device(a,c)
     call mp_sum(  a(1:lbnd), intra_bgrp_comm )
     call mp_sum(  c(1:lbnd), intra_bgrp_comm )
     !$acc end host_data
     lbnd=0
     CALL start_clock('loop4')
     do ibnd = n_start, n_end ; ibnd_ = ibnd - n_start + 1
        if (conv (ibnd) .eq.0) then
           lbnd=lbnd+1
           !$acc update host(a,c)
           dclambda = - a(lbnd) / c(lbnd)
           dclambda1 = CONJG(dclambda)
           !
           !    move to new position
           !
           !$acc data present(dpsi,h)
           !$acc host_data use_device(dpsi,h)
           call zaxpy (ndmx*npol, dclambda, h(1,ibnd_), 1, dpsi(1,ibnd), 1)
           !$acc end host_data
           !$acc end data
           !
           !    update to get the gradient
           !
           !g=g+lam
           !$acc data present(t,g,ts,gs)
           !$acc host_data use_device(t,g,ts,gs)           
           call zaxpy (ndmx*npol, dclambda, t(1,lbnd), 1, g(1,ibnd_), 1)
           CALL zaxpy (ndmx*npol, dclambda1, ts(1,lbnd), 1, gs(1,ibnd_), 1)
           !$acc end host_data
           !$acc end data
           !
           !    save current (now old) h and rho for later use
           !
           !$acc data present(hold,h,hsold,hs)
           !$acc host_data use_device(h,hold,hs,hsold)
           call zcopy (ndmx*npol, h(1,ibnd_), 1, hold(1,ibnd_), 1)
           CALL zcopy (ndmx*npol, hs(1,ibnd_), 1, hsold(1,ibnd_), 1)
           !$acc end host_data
           !$acc end data
           rhoold (ibnd_) = rho (ibnd_)
        endif
     enddo
     CALL stop_clock('loop4')
  enddo

100 continue
  !$acc exit data delete(evq,a,c,g,gs,h,hs,h_diag,d0psi,hold,hsold,t,ts,eu,euc,e) copyout(dpsi)
  ! deallocate workspace not needed anymore
  deallocate (eu) ; deallocate (rho, rhoold) ; deallocate (a,c) ; deallocate (g, t, h, hold)
  deallocate (gs, ts, hs, hsold, indb)
  deallocate (euc)

  ! wait for all bgrp to complete their task
  CALL mp_barrier( inter_bgrp_comm )

  ! check if all root converged across all bgrp
  call mp_sum( conv, inter_bgrp_comm )
  conv_root = .true.
  do ibnd = 1, nbnd
     conv_root = conv_root.and. (conv (ibnd) .eq.1)
  enddo
  deallocate (conv)

  ! collect the result
  if (n_start > 1 ) dpsi(:, 1:n_start-1) = (0.d0,0.d0) ; if (n_end < nbnd) dpsi(:, n_end+1:nbnd) = (0.d0,0.d0)
  call mp_sum( dpsi, inter_bgrp_comm )

  call mp_sum( kter_eff, inter_bgrp_comm )
  kter = kter_eff

  ! restore the value of use_bgrp_in_hpsi to its saved value
  use_bgrp_in_hpsi = lsave_use_bgrp_in_hpsi


  call stop_clock ('ccgsolve')
  return
end subroutine ccgsolve_all
