!
! Copyright (C) 2003-2020 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
MODULE bfgs_module
   !----------------------------------------------------------------------------
   !! Ionic relaxation through the Newton-Raphson optimization scheme
   !! based on the Broyden-Fletcher-Goldfarb-Shanno algorithm for the
   !! estimate of the inverse Hessian matrix.  
   !! The ionic relaxation is performed converting cartesian (and cell) 
   !! positions into internal coordinates.  
   !! The algorithm uses a "trust radius" line search based on Wolfe 
   !! conditions. Steps are rejected until the first Wolfe condition
   !! (sufficient energy decrease) is satisfied. Updated step length
   !! is estimated from quadratic interpolation.  
   !! When the step is accepted inverse hessian is updated according to 
   !! BFGS scheme and a new search direction is obtained from NR or GDIIS
   !! method. The corresponding step length is limited by trust_radius_max 
   !! and can't be larger than the previous step multiplied by a certain 
   !! factor determined by Wolfe and other convergence conditions.
   !
   !! Originally written (5/12/2003) and maintained (2003-2007) by 
   !! Carlo Sbraccia.  
   !! Modified for variable-cell-shape relaxation (2007-2008) by 
   !! Javier Antonio Montoya, Lorenzo Paulatto and Stefano de Gironcoli.  
   !! Re-analyzed by Stefano de Gironcoli (2010).  
   !! FCP case added by Minoru Otani and collaborators (2020).
   !
   !! References:  
   !! 1) Roger Fletcher, Practical Methods of Optimization, John Wiley and
   !!    Sons, Chichester, 2nd edn, 1987.  
   !! 2) Salomon R. Billeter, Alexander J. Turner, Walter Thiel,
   !!    Phys. Chem. Chem. Phys. 2, 2177 (2000).  
   !! 3) Salomon R. Billeter, Alessandro Curioni, Wanda Andreoni,
   !!    Comput. Mat. Science 27, 437, (2003).  
   !! 4) Ren Weiqing, PhD Thesis: Numerical Methods for the Study of Energy
   !!    Landscapes and Rare Events.
   !
   USE kinds,     ONLY : DP
   USE constants, ONLY : eps4, eps8, eps16, RYTOEV
   !
   USE basic_algebra_routines
   USE matrix_inversion
   !
   IMPLICIT NONE
   !
   PRIVATE
   !
   ! ... public methods
   !
   PUBLIC :: bfgs, init_bfgs, terminate_bfgs, bfgs_get_n_iter  
   !
   ! ... global module variables
   !
   SAVE
   !
   CHARACTER (len=18) :: fname="energy"
   !! name of the function to be minimized
   CHARACTER (len=320):: bfgs_file=" "
   !! name of file (with path) used to store and retrieve the status
   !
   REAL(DP), ALLOCATABLE :: pos(:)
   !! positions + cell
   REAL(DP), ALLOCATABLE :: grad(:)
   !! gradients + cell_force
   REAL(DP), ALLOCATABLE :: pos_p(:)
   !! positions at the previous accepted iteration
   REAL(DP), ALLOCATABLE :: grad_p(:)
   !! gradients at the previous accepted iteration
   REAL(DP), ALLOCATABLE :: inv_hess(:,:)
   !! inverse hessian matrix (updated using BFGS formula)
   REAL(DP), ALLOCATABLE :: metric(:,:)
   REAL(DP), ALLOCATABLE :: inv_metric(:,:)
   REAL(DP), ALLOCATABLE :: h_block(:,:)
   REAL(DP), ALLOCATABLE :: hinv_block(:,:)
   REAL(DP), ALLOCATABLE :: step(:)
   !! the (new) search direction (normalized NR step)
   REAL(DP), ALLOCATABLE :: step_old(:)
   !! the previous search direction (normalized NR step)
   REAL(DP), ALLOCATABLE :: pos_old(:,:)
   !! list of m old positions - used only by gdiis
   REAL(DP), ALLOCATABLE :: grad_old(:,:)
   !! list of m old gradients - used only by gdiis
   REAL(DP), ALLOCATABLE :: pos_best(:)
   !! best extrapolated positions - used only by gdiis
   !
   REAL(DP) :: nr_step_length
   !! length of (new) Newton-Raphson step
   REAL(DP) :: nr_step_length_old
   !! length of previous Newton-Raphson step
   REAL(DP) :: trust_radius
   !! new displacement along the search direction
   REAL(DP) :: trust_radius_old
   !! old displacement along the search direction
   REAL(DP) :: energy_p
   !! energy at previous accepted iteration
   !
   INTEGER :: scf_iter
   !! number of scf iterations
   INTEGER :: bfgs_iter
   !! number of bfgs iterations
   INTEGER :: gdiis_iter
   !! number of gdiis iterations
   INTEGER :: tr_min_hit = 0
   !! set to 1 if the trust_radius has already been set to the minimum value
   !! at the previous step. Set to 2 if trust_radius is reset again: exit.
   LOGICAL :: conv_bfgs = .FALSE.
   !! .TRUE. when bfgs convergence has been achieved
   !
   ! ... default values for bfgs_ndim, trust_radius_max, trust_radius_min,
   ! ... trust_radius_ini, w_1, w_2, are set in Modules/read_namelist.f90
   ! ... (SUBROUTINE ions_defaults) and can be assigned in the input
   !
   LOGICAL :: use_gdiis_step
   LOGICAL :: bfgs_initialized = .FALSE.
   INTEGER :: stdout
   !! standard output for writing
   INTEGER :: bfgs_ndim
   !! dimension of the subspace for GDIIS for standard BFGS algorithm,
   !! \(\text{bfgs_ndim}=1\)
   REAL(DP) :: trust_radius_ini
   !! suggested initial displacement
   REAL(DP) :: trust_radius_min
   !! minimum allowed displacement
   REAL(DP) :: trust_radius_max
   !! maximum allowed displacement
   !
   ! ... parameters for Wolfe conditions
   !
   REAL(DP) :: w_1
   !! 1st Wolfe condition: sufficient energy decrease
   REAL(DP) :: w_2
   !! 2nd Wolfe condition: sufficient gradient decrease
   !
CONTAINS
   !
   !------------------------------------------------------------------------
   SUBROUTINE init_bfgs( stdout_, bfgs_ndim_, trust_radius_max_, &
                   trust_radius_min_, trust_radius_ini_, w_1_, w_2_, use_gdiis_step_)
     !------------------------------------------------------------------------
     !! set values for several parameters of the algorithm
     !
     INTEGER, INTENT(IN) :: &
          stdout_, &
          bfgs_ndim_
     REAL(DP), INTENT(IN)  :: &
        trust_radius_ini_, &
        trust_radius_min_, &
        trust_radius_max_, &
        w_1_,              &
        w_2_
      LOGICAL,INTENT(IN) :: use_gdiis_step_
     !
     stdout            = stdout_
     bfgs_ndim         = bfgs_ndim_
     trust_radius_max  = trust_radius_max_
     trust_radius_min  = trust_radius_min_
     trust_radius_ini  = trust_radius_ini_
     w_1               = w_1_
     w_2               = w_2_
     use_gdiis_step = use_gdiis_step_
     bfgs_initialized  = .true.
     !
   END SUBROUTINE init_bfgs
   !
   !------------------------------------------------------------------------
   SUBROUTINE bfgs( filebfgs, pos_in, h, nelec, energy, &
                   grad_in, fcell, iforceh, felec, &
                   energy_thr, grad_thr, cell_thr, fcp_thr, &
                   energy_error, grad_error, cell_error, fcp_error, &
                   lmovecell, lfcp, fcp_cap, fcp_hess, step_accepted, &
                   stop_bfgs, failed, istep )
      !------------------------------------------------------------------------
      !! BFGS algorithm.
      !
      USE constants,              ONLY : ry_kbar
      !
      IMPLICIT NONE
      !
      REAL(DP),         INTENT(INOUT) :: pos_in(:)
      !! vector containing 3N coordinates of the system ( x )
      REAL(DP),         INTENT(INOUT) :: h(3,3)
      !! 3x3 matrix of the primitive lattice vectors
      REAL(DP),         INTENT(INOUT) :: nelec
      !! number of electrons (for FCP)
      REAL(DP),         INTENT(INOUT) :: energy
      !! energy of the system ( V(x) )
      REAL(DP),         INTENT(INOUT) :: grad_in(:)
      !! vector containing 3N components of grad( V(x) )
      REAL(DP),         INTENT(INOUT) :: fcell(3,3)
      !! 3x3 matrix containing the stress tensor
      REAL(DP),         INTENT(INOUT) :: felec 
      !! force on FCP
      CHARACTER(LEN=*), INTENT(IN)    :: filebfgs
      !! file name for storing and retrieving data
      REAL(DP),         INTENT(IN)    :: energy_thr
      !! threshold on energy difference for BFGS convergence
      REAL(DP),         INTENT(IN)    :: grad_thr
      !! threshold on grad difference for BFGS convergence
      !! the largest component of grad( V(x) ) is considered
      REAL(DP),         INTENT(IN)    :: cell_thr
      !! threshold on grad of cell for BFGS convergence
      REAL(DP),         INTENT(IN)    :: fcp_thr
      !! treshold on grad of FCP for BFGS convergence
      INTEGER,          INTENT(IN)    :: iforceh(3,3)
      !! 3x3 matrix containing cell constraints
      LOGICAL,          INTENT(IN)    :: lmovecell
      LOGICAL,          INTENT(IN)    :: lfcp
      !! include FCP, or not ?
      REAL(DP),         INTENT(IN)    :: fcp_cap
      !! capacitance for FCP
      REAL(DP),         INTENT(IN)    :: fcp_hess
      !! hessian for FCP
      REAL(DP),         INTENT(OUT)   :: energy_error
      !! energy difference \(\| V(x_i) - V(x_{i-1}) \|\)
      REAL(DP),         INTENT(OUT)   :: grad_error
      !! the largest component of \(\|\text{grad}(V(x_i)) -
      !! \text{grad}(V(x_{i-1}))\|\)
      REAL(DP),         INTENT(OUT)   :: cell_error
      !! the largest component of: omega*(stress-press*I)
      REAL(DP),         INTENT(OUT)   :: fcp_error
      !! \(\| \text{grad}(V(x_{FCP})) \|\)
      LOGICAL,          INTENT(OUT)   :: step_accepted
      !! .TRUE. if a new BFGS step is done
      LOGICAL,          INTENT(OUT)   :: stop_bfgs
      !! .TRUE. if BFGS convergence has been achieved
      LOGICAL,          INTENT(OUT)   :: failed
      !! .TRUE. if BFGS failed
      INTEGER,          INTENT(OUT)   :: istep
      !
      ! ... local variables
      !
      INTEGER  :: n, i, j, k, nat
      LOGICAL  :: lwolfe
      REAL(DP) :: dE0s, den
      REAL(DP), ALLOCATABLE :: step_tmp(:)
      ! ... for scaled coordinates
      REAL(DP) :: hinv(3,3),g(3,3),ginv(3,3), omega
      !
      ! ... additional dimensions of cell and FCP
      INTEGER, PARAMETER :: NADD = 9 + 1
      !
      !
      failed = .FALSE.
      !
      IF ( .NOT.bfgs_initialized ) CALL errore('bfgs',' not initialized',1)
      !
      IF ( bfgs_file == " ") bfgs_file = TRIM(filebfgs)
      lwolfe=.false.
      n = SIZE( pos_in ) + NADD
      nat = size (pos_in) / 3
      if (nat*3 /= size (pos_in)) call errore('bfgs',' strange dimension',1)
      !
      ! ... work-space allocation
      !
      ALLOCATE( pos(    n ) )
      ALLOCATE( grad(   n ) )
      !
      ALLOCATE( grad_old( n, bfgs_ndim ) )
      ALLOCATE( pos_old(  n, bfgs_ndim ) )
      !
      ALLOCATE( inv_hess( n, n ) )
      !
      ALLOCATE( pos_p(    n ) )
      ALLOCATE( grad_p(   n ) )
      ALLOCATE( step(     n ) )
      ALLOCATE( step_old( n ) )
      ALLOCATE( pos_best( n ) )
      ! ... scaled coordinates work-space
      ALLOCATE( hinv_block( n-NADD, n-NADD ) )
      ! ... cell related work-space
      ALLOCATE( metric( n , n  ) )
      ALLOCATE( inv_metric( n , n  ) )
      !
      ! ... the BFGS file read (pos & grad) in scaled coordinates
      !
      call invmat(3, h, hinv, omega)
      ! volume is defined to be positive even for left-handed vector triplet
      omega = abs(omega) 
      !
      hinv_block = 0.d0
      FORALL ( k=0:nat-1, i=1:3, j=1:3 ) hinv_block(i+3*k,j+3*k) = hinv(i,j)
      !
      ! ... generate metric to work with scaled ionic coordinates
      g = MATMUL(TRANSPOSE(h),h)
      call invmat(3,g,ginv)
      metric = 0.d0
      FORALL ( k=0:nat-1,   i=1:3, j=1:3 ) metric(i+3*k,j+3*k) = g(i,j)
      FORALL ( k=nat:nat+2, i=1:3, j=1:3 ) metric(i+3*k,j+3*k) = 0.04d0 * omega * ginv(i,j)
      !
      IF ( lfcp ) THEN
         IF ( fcp_cap > eps4 ) THEN
            metric(n,n) = (0.5d0 / (0.1d0 * fcp_cap)) ** 2
         ELSE
            metric(n,n) = 1.d0
         END IF
      ELSE
         metric(n,n) = 1.d0
      END IF
      !
      ! ... generate inverse of metric
      inv_metric = 0.d0
      FORALL ( k=0:nat-1,   i=1:3, j=1:3 ) inv_metric(i+3*k,j+3*k) = ginv(i,j)
      FORALL ( k=nat:nat+2, i=1:3, j=1:3 ) inv_metric(i+3*k,j+3*k) = g(i,j) / (0.04d0 * omega)
      inv_metric(n,n) = 1.d0 / metric(n,n)
      !
      ! ... generate bfgs vectors for the degrees of freedom and their gradients
      pos = 0.d0
      pos(1:n-NADD) = pos_in
      if (lmovecell) FORALL( i=1:3, j=1:3)  pos( n-NADD + j+3*(i-1) ) = h(i,j)
      if (lfcp) pos( n ) = nelec
      grad = 0.d0
      grad(1:n-NADD) = grad_in
      if (lmovecell) FORALL( i=1:3, j=1:3) grad( n-NADD + j+3*(i-1) ) = fcell(i,j)*iforceh(i,j)
      if (lfcp) grad( n ) = felec
      !
      ! if the cell moves the quantity to be minimized is the enthalpy
      fname = "energy"
      IF ( lmovecell ) fname="enthalpy"
      IF ( lfcp ) fname = "grand-" // TRIM(fname)
      !
      CALL read_bfgs_file( pos, grad, energy, n, lfcp, fcp_hess )
      !
      scf_iter = scf_iter + 1
      istep    = scf_iter
      !
      ! ... convergence is checked here
      !
      IF (scf_iter .EQ. 1) THEN
          ! Compute predicted energy change for the initial iteration, since
          ! there is no previous energy.
          ! See section 6.2 of Nocedal and Wright "Numerical Optimization",
          ! there it is called "predicted reduction".
          ! Note that since it is first iteration, hessian is metric.
          ALLOCATE (step_tmp( n ) )
          !
          step_tmp(:) = inv_metric(:,:) .times. grad(:)
          if (lmovecell) FORALL( i=1:3, j=1:3) step_tmp( n-NADD+j+3*(i-1) ) = &
              step_tmp( n-NADD+j+3*(i-1) )*iforceh(i,j)
          !
          energy_error = abs(grad(:) .dot. step_tmp(:) + &
              0.5_DP * (step_tmp(:) .dot. (metric(:, :) .times. step_tmp(:))))
          !
          DEALLOCATE( step_tmp )
      ELSE
          energy_error = ABS( energy_p - energy )
      END IF
      !
      ! ... obscure PGI bug as of v.17.4
      !
#if defined (__PGI)
      grad_in = MATMUL( TRANSPOSE(hinv_block), grad(1:n-NADD) )
      grad_error = MAXVAL( ABS( grad_in ) )
#else
      grad_error = MAXVAL( ABS( MATMUL( TRANSPOSE(hinv_block), grad(1:n-NADD)) ) )
#endif
      conv_bfgs = energy_error < energy_thr
      conv_bfgs = conv_bfgs .AND. ( grad_error < grad_thr )
      !
      IF( lmovecell) THEN
          cell_error = MAXVAL( ABS( MATMUL ( TRANSPOSE ( RESHAPE( grad(n-NADD+1:n-1), (/ 3, 3 /) ) ),&
                                             TRANSPOSE(h) ) ) ) / omega
          conv_bfgs = conv_bfgs .AND. ( cell_error < cell_thr ) 
#undef DEBUG
#if defined(DEBUG)
           write (*,'(3f15.10)') TRANSPOSE ( RESHAPE( grad(n-NADD+1:n-1), (/ 3, 3 /) ) )
           write (*,*)
           write (*,'(3f15.10)') TRANSPOSE(h)
           write (*,*)
           write (*,'(3f15.10)') MATMUL (TRANSPOSE( RESHAPE( grad(n-NADD+1:n-1), (/ 3, 3 /) ) ),&
                                             TRANSPOSE(h) ) / omega
           write (*,*)
           write (*,*) cell_error/cell_thr*0.5d0
#endif
      END IF
      !
      IF( lfcp ) THEN
          fcp_error = ABS( grad( n ) )
          conv_bfgs = conv_bfgs .AND. ( fcp_error < fcp_thr )
      END IF
      !
      ! ... converged (or useless to go on): quick return
      !
      IF ( .NOT. conv_bfgs .AND. ( tr_min_hit > 1 ) ) THEN
         CALL infomsg( 'bfgs',&
              'history already reset at previous step: exiting' )
         failed = .TRUE.
      END IF
      !
      conv_bfgs = conv_bfgs .OR. ( tr_min_hit > 1 )
      !
      WRITE(stdout, '(5X,"Energy error",T30,"= ",1PE12.1," Ry")') energy_error
      WRITE(stdout, '(5X,"Gradient error",T30,"= ",1PE12.1," Ry/Bohr")') grad_error
      IF( lmovecell ) WRITE(stdout, &
          '(5X,"Cell gradient error",T30,"= ",1PE12.1," kbar")') &
          cell_error * ry_kbar
      IF( lfcp ) WRITE(stdout, &
         '(5X,"FCP gradient error",T30,"= ",1PE12.1,/)') fcp_error
      !
      IF ( conv_bfgs ) GOTO 1000
      !
      ! ... some output is written
      !
      WRITE( UNIT = stdout, &
           & FMT = '(/,5X,"number of scf cycles",T30,"= ",I3)' ) scf_iter
      WRITE( UNIT = stdout, &
           & FMT = '(5X,"number of bfgs steps",T30,"= ",I3,/)' ) bfgs_iter
      IF ( scf_iter > 1 ) WRITE( UNIT = stdout, &
           & FMT = '(5X,A," old",T30,"= ",F18.10," Ry")' ) fname,energy_p
      WRITE( UNIT = stdout, &
           & FMT = '(5X,A," new",T30,"= ",F18.10," Ry",/)' ) fname,energy
      !
      ! ... the bfgs algorithm starts here
      !
      IF ( .NOT. energy_wolfe_condition( energy ) .AND. (scf_iter > 1) .AND. (.NOT. lfcp) ) THEN
         !
         ! ... Wolfe condition not met: the previous step is rejected,
         ! ... line search goes on
         ! ... Note that for the FCP case, the Wolfe condition is ignored
         !
         step_accepted = .FALSE.
         !
         WRITE( UNIT = stdout, &
              & FMT = '(5X,"CASE: ",A,"_new > ",A,"_old",/)' ) fname,fname
         !
         ! ... the new trust radius is obtained by quadratic interpolation
         !
         ! ... E(s) = a*s*s + b*s + c      ( we use E(0), dE(0), E(s') )
         !
         ! ... s_min = - 0.5*( dE(0)*s'*s' ) / ( E(s') - E(0) - dE(0)*s' )
         !
         if (abs(scnorm(step_old(:))-1._DP) > 1.d-10) call errore('bfgs', &
                  ' step_old is NOT normalized ',1)
         ! (normalized) search direction is the same as in previous step
         step(:) = step_old(:)
         !
         dE0s = ( grad_p(:) .dot. step(:) ) * trust_radius_old
         IF (dE0s > 0._DP ) CALL errore( 'bfgs', &
                  'dE0s is positive which should never happen', 1 )
         den = energy - energy_p - dE0s
         !
         ! estimate new trust radius by interpolation
         trust_radius = - 0.5_DP*dE0s*trust_radius_old / den
         !
         WRITE( UNIT = stdout, &
              & FMT = '(5X,"new trust radius",T30,"= ",F18.10," bohr")' ) &
              trust_radius
         !
         ! ... values from the last successful bfgs step are restored
         !
         pos(:)  = pos_p(:)
         energy  = energy_p
         grad(:) = grad_p(:)
         !
         IF ( trust_radius < trust_radius_min ) THEN
            !
            ! ... the history is reset ( the history can be reset at most two
            ! ... consecutive times )
            !
            WRITE( UNIT = stdout, &
                   FMT = '(/,5X,"trust_radius < trust_radius_min")' )
            WRITE( UNIT = stdout, FMT = '(/,5X,"resetting bfgs history",/)' )
            !
            ! ... if tr_min_hit=1 the history has already been reset at the 
            ! ... previous step : something is going wrong
            !
            IF ( tr_min_hit == 1 ) THEN
               CALL infomsg( 'bfgs', &
                            'history already reset at previous step: stopping' )
               tr_min_hit = 2 
            ELSE
               tr_min_hit = 1
            END IF
            !
            CALL reset_bfgs( n, lfcp, fcp_hess )
            !
            step(:) = - ( inv_hess(:,:) .times. grad(:) )
            if (lmovecell) FORALL( i=1:3, j=1:3) step( n-NADD+j+3*(i-1) ) = step( n-NADD+j+3*(i-1) )*iforceh(i,j)
            ! normalize step but remember its length
            nr_step_length = scnorm(step)
            step(:) = step(:) / nr_step_length
            !
            trust_radius = min(trust_radius_ini, nr_step_length)
            !
         ELSE
            !
            tr_min_hit = 0
            !
         END IF
         !
      ELSE
         !
         ! ... a new bfgs step is done
         !
         bfgs_iter = bfgs_iter + 1
         !
         IF ( bfgs_iter == 1 ) THEN
            !
            ! ... first iteration
            !
            step_accepted = .FALSE.
            !
         ELSE
            !
            step_accepted = .TRUE.
            !
            nr_step_length_old = nr_step_length
            !
            IF ( .NOT. lfcp ) THEN
               !
               WRITE( UNIT = stdout, &
                    & FMT = '(5X,"CASE: ",A,"_new < ",A,"_old",/)' ) fname,fname
               !
            END IF
            !
            CALL check_wolfe_conditions( lwolfe, energy, grad )
            !
            CALL update_inverse_hessian( pos, grad, n, lfcp, fcp_hess )
            !
         END IF
         ! compute new search direction and store NR step length
         IF ( bfgs_ndim > 1 ) THEN
            !
            ! ... GDIIS extrapolation
            !
            CALL gdiis_step()
            !
         ELSE
            !
            ! ... standard Newton-Raphson step
            !
            step(:) = - ( inv_hess(:,:) .times. grad(:) )
            if (lmovecell) FORALL( i=1:3, j=1:3) step( n-NADD+j+3*(i-1) ) = step( n-NADD+j+3*(i-1) )*iforceh(i,j)
            !
         END IF
         IF ( ( grad(:) .dot. step(:) ) > 0.0_DP ) THEN
            !
            WRITE( UNIT = stdout, &
                   FMT = '(5X,"uphill step: resetting bfgs history",/)' )
            !
            CALL reset_bfgs( n, lfcp, fcp_hess )
            step(:) = - ( inv_hess(:,:) .times. grad(:) )
            if (lmovecell) FORALL( i=1:3, j=1:3) step( n-NADD+j+3*(i-1) ) = step( n-NADD+j+3*(i-1) )*iforceh(i,j)
            !
         END IF
         !
         ! normalize the step and save the step length
         nr_step_length = scnorm(step)
         step(:) = step(:) / nr_step_length
         !
         ! ... the new trust radius is computed
         !
         IF ( bfgs_iter == 1 ) THEN
            !
            trust_radius = min(trust_radius_ini, nr_step_length)
            tr_min_hit = 0
            !
         ELSE
            !
            CALL compute_trust_radius( lwolfe, energy, grad, n, lfcp, fcp_hess )
            !
         END IF
         !
         WRITE( UNIT = stdout, &
              & FMT = '(5X,"new trust radius",T30,"= ",F18.10," bohr")' ) &
              trust_radius
         !
      END IF
      !
      ! ... step along the bfgs direction
      !
      IF ( nr_step_length < eps16 ) &
         CALL errore( 'bfgs', 'NR step-length unreasonably short', 1 )
      !
      ! ... information required by next iteration is saved here ( this must
      ! ... be done before positions are updated )
      !
      CALL write_bfgs_file( pos, energy, grad )
      !
      ! ... positions and cell are updated
      !
      IF (use_gdiis_step .and. bfgs_ndim .gt. 1 ) THEN 
         pos(:) = pos + nr_step_length * step(:)
      ELSE 
         pos(:) = pos(:) + trust_radius * step(:)
      END IF 
      !
1000  stop_bfgs = conv_bfgs
      ! ... input ions+cell variables
      pos_in = pos(1:n-NADD)
      IF ( lmovecell ) FORALL( i=1:3, j=1:3) h(i,j) = pos( n-NADD + j+3*(i-1) )
      IF ( lfcp ) nelec = pos( n )
      ! ... update forces
      grad_in = grad(1:n-NADD)
      !
      ! ... work-space deallocation
      !
      DEALLOCATE( pos )
      DEALLOCATE( grad )
      DEALLOCATE( pos_p )
      DEALLOCATE( grad_p )
      DEALLOCATE( pos_old )
      DEALLOCATE( grad_old )
      DEALLOCATE( inv_hess )
      DEALLOCATE( step )
      DEALLOCATE( step_old )
      DEALLOCATE( pos_best )
      DEALLOCATE( hinv_block )
      DEALLOCATE( metric )
      DEALLOCATE( inv_metric )
      !
      RETURN
      !
   CONTAINS
      !
      !--------------------------------------------------------------------
      SUBROUTINE gdiis_step()
         !--------------------------------------------------------------------
         !! GDIIS extrapolation
         !
         USE basic_algebra_routines
         IMPLICIT NONE
         !
         REAL(DP), ALLOCATABLE :: res(:,:), overlap(:,:), work(:)
         INTEGER,  ALLOCATABLE :: iwork(:)
         INTEGER               :: k, k_m, info
         REAL(DP)              :: gamma0
         !
         !
         gdiis_iter = gdiis_iter + 1
         !
         k   = MIN( gdiis_iter, bfgs_ndim )
         k_m = k + 1
         !
         ALLOCATE( res( n, k ) )
         ALLOCATE( overlap( k_m, k_m ) )
         ALLOCATE( work( k_m ), iwork( k_m ) )
         !
         work(:)  = 0.0_DP
         iwork(:) = 0
         !
         ! ... the new direction is added to the workspace
         !
         DO i = bfgs_ndim, 2, -1
            !
            pos_old(:,i)  = pos_old(:,i-1)
            grad_old(:,i) = grad_old(:,i-1)
            !
         END DO
         !
         pos_old(:,1)  = pos(:)
         grad_old(:,1) = grad(:)
         !
         ! ... |res_i> = H^-1 \times |g_i>
         !
         CALL DGEMM( 'N', 'N', n, k, n, 1.0_DP, &
                     inv_hess, n, grad_old, n, 0.0_DP, res, n )
         !
         ! ... overlap_ij = <grad_i|res_j>
         !
         CALL DGEMM( 'T', 'N', k, k, n, 1.0_DP, &
                     res, n, res, n, 0.0_DP, overlap, k_m )
         !
         overlap( :, k_m) = 1.0_DP
         overlap(k_m, : ) = 1.0_DP
         overlap(k_m,k_m) = 0.0_DP
         !
         ! make sure the overlap matrix is not singular 
         !
         FORALL( i = 1:k ) overlap(i,i) = overlap(i,i) * ( 1.0_DP + eps8 ) ; overlap(k_m,k_m) = eps8
         !
         ! ... overlap is inverted via Bunch-Kaufman diagonal pivoting method
         !
         CALL DSYTRF( 'U', k_m, overlap, k_m, iwork, work, k_m, info )
         CALL DSYTRI( 'U', k_m, overlap, k_m, iwork, work, info )
         CALL errore( 'gdiis_step', 'error in Bunch-Kaufman inversion', info )
         !
         ! ... overlap is symmetrised
         !
         FORALL( i = 1:k_m, j = 1:k_m, j > i ) overlap(j,i) = overlap(i,j)
         !
         pos_best(:) = 0.0_DP
         step(:)     = 0.0_DP
         !
         DO i = 1, k
            !
            gamma0 = overlap(k_m,i)
            !
            pos_best(:) = pos_best(:) + gamma0*pos_old(:,i)
            !
            step(:) = step(:) - gamma0*res(:,i)
            !
         END DO
         !
         ! ... the step must be consistent with the last positions
         !
         step(:) = step(:) + ( pos_best(:) - pos(:) )
         !
         IF ( ( grad(:) .dot. step(:) ) > 0.0_DP ) THEN
            !
            ! ... if the extrapolated direction is uphill use only the
            ! ... last gradient and reset gdiis history
            !
            step(:) = - ( inv_hess(:,:) .times. grad(:) )
            if (lmovecell) FORALL( i=1:3, j=1:3) step( n-NADD+j+3*(i-1) ) = step( n-NADD+j+3*(i-1) )*iforceh(i,j)
            !
            gdiis_iter = 0
            !
         END IF
         !
         DEALLOCATE( res, overlap, work, iwork )
         !
      END SUBROUTINE gdiis_step
      !
   END SUBROUTINE bfgs
   !
   !------------------------------------------------------------------------
   SUBROUTINE reset_bfgs( n, lfcp, fcp_hess )
      !------------------------------------------------------------------------
      !! \(\text{inv_hess}\) in re-initialized to the initial guess 
      !! defined as the inverse metric.
      !
      INTEGER,  INTENT(IN) :: n
      LOGICAL,  INTENT(IN) :: lfcp
      REAL(DP), INTENT(IN) :: fcp_hess
      !
      inv_hess = inv_metric
      !
      IF ( lfcp ) THEN
         !
         IF ( fcp_hess > eps4 ) THEN
            !
            inv_hess(n,n) = fcp_hess
            !
         END IF
         !
      END IF
      !
      gdiis_iter = 0
      !
   END SUBROUTINE reset_bfgs
   !
   !------------------------------------------------------------------------
   SUBROUTINE read_bfgs_file( pos, grad, energy, n, lfcp, fcp_hess )
      !------------------------------------------------------------------------
      !
      IMPLICIT NONE
      !
      REAL(DP),         INTENT(INOUT) :: pos(:)
      REAL(DP),         INTENT(INOUT) :: grad(:)
      INTEGER,          INTENT(IN)    :: n
      LOGICAL,          INTENT(IN)    :: lfcp
      REAL(DP),         INTENT(IN)    :: fcp_hess
      REAL(DP),         INTENT(INOUT) :: energy
      !
      INTEGER            :: iunbfgs
      LOGICAL            :: file_exists
      !
      !
      INQUIRE( FILE = bfgs_file, EXIST = file_exists )
      !
      IF ( file_exists ) THEN
         !
         ! ... bfgs is restarted from file
         !
         OPEN( NEWUNIT = iunbfgs, FILE = bfgs_file, &
               STATUS = 'UNKNOWN', ACTION = 'READ' )
         !
         READ( iunbfgs, * ) pos_p
         READ( iunbfgs, * ) grad_p
         READ( iunbfgs, * ) scf_iter
         READ( iunbfgs, * ) bfgs_iter
         READ( iunbfgs, * ) gdiis_iter
         READ( iunbfgs, * ) energy_p
         READ( iunbfgs, * ) pos_old
         READ( iunbfgs, * ) grad_old
         READ( iunbfgs, * ) inv_hess
         READ( iunbfgs, * ) tr_min_hit
         READ( iunbfgs, * ) nr_step_length
         !
         CLOSE( UNIT = iunbfgs )
         !
         step_old = ( pos(:) - pos_p(:) ) 
         trust_radius_old = scnorm( step_old )
         step_old = step_old / trust_radius_old
         !
      ELSE
         !
         ! ... bfgs initialization
         !
         WRITE( UNIT = stdout, FMT = '(/,5X,"BFGS Geometry Optimization")' )
         !
         ! initialize the inv_hess to the inverse of the metric
         inv_hess = inv_metric
         !
         IF ( lfcp ) THEN
            !
            IF ( fcp_hess > eps4 ) THEN
               !
               inv_hess(n,n) = fcp_hess
               !
            END IF
            !
         END IF
         !
         pos_p      = 0.0_DP
         grad_p     = 0.0_DP
         scf_iter   = 0
         bfgs_iter  = 0
         gdiis_iter = 0
         energy_p   = energy
         step_old   = 0.0_DP
         nr_step_length = 0.0_DP
         !
         trust_radius_old = trust_radius_ini
         !
         pos_old(:,:)  = 0.0_DP
         grad_old(:,:) = 0.0_DP
         !
         tr_min_hit = 0
         !
      END IF
      !
   END SUBROUTINE read_bfgs_file
   !
   !------------------------------------------------------------------------
   SUBROUTINE write_bfgs_file( pos, energy, grad )
      !------------------------------------------------------------------------
      !
      IMPLICIT NONE
      !
      REAL(DP),         INTENT(IN) :: pos(:)
      REAL(DP),         INTENT(IN) :: energy
      REAL(DP),         INTENT(IN) :: grad(:)
      !
      INTEGER :: iunbfgs
      !
      OPEN( NEWUNIT = iunbfgs, FILE = bfgs_file, STATUS = 'UNKNOWN', ACTION = 'WRITE' )
      !
      WRITE( iunbfgs, * ) pos
      WRITE( iunbfgs, * ) grad
      WRITE( iunbfgs, * ) scf_iter
      WRITE( iunbfgs, * ) bfgs_iter
      WRITE( iunbfgs, * ) gdiis_iter
      WRITE( iunbfgs, * ) energy
      WRITE( iunbfgs, * ) pos_old
      WRITE( iunbfgs, * ) grad_old
      WRITE( iunbfgs, * ) inv_hess
      WRITE( iunbfgs, * ) tr_min_hit
      WRITE( iunbfgs, * ) nr_step_length
      !
      CLOSE( UNIT = iunbfgs )
      !
   END SUBROUTINE write_bfgs_file
   !
   SUBROUTINE update_inverse_hessian( pos, grad, n, lfcp, fcp_hess )
      !------------------------------------------------------------------------
      !
      IMPLICIT NONE
      !
      REAL(DP), INTENT(IN)  :: pos(:)
      REAL(DP), INTENT(IN)  :: grad(:)
      INTEGER,  INTENT(IN)  :: n
      LOGICAL,  INTENT(IN)  :: lfcp
      REAL(DP), INTENT(IN)  :: fcp_hess
      INTEGER               :: info
      !
      REAL(DP), ALLOCATABLE :: y(:), s(:)
      REAL(DP), ALLOCATABLE :: Hy(:), yH(:)
      REAL(DP)              :: sdoty, sBs, Theta
      REAL(DP), ALLOCATABLE :: B(:,:)
      !
      ALLOCATE( y( n ), s( n ), Hy( n ), yH( n ) )
      !
      s(:) = pos(:)  - pos_p(:)
      y(:) = grad(:) - grad_p(:)
      !
      sdoty = ( s(:) .dot. y(:) )
      !
      IF ( ABS( sdoty ) < eps16 ) THEN
         !
         ! ... the history is reset
         !
         WRITE( stdout, '(/,5X,"WARNING: unexpected ", &
                         &     "behaviour in update_inverse_hessian")' )
         WRITE( stdout, '(  5X,"         resetting bfgs history",/)' )
         !
         CALL reset_bfgs( n, lfcp, fcp_hess )
         !
         RETURN
         !
      ELSE
!       Conventional Curvature Trap here
!       See section 18.2 (p538-539 ) of Nocedal and Wright "Numerical
!       Optimization"for instance
!       LDM Addition, April 2011
!
!       While with the Wolfe conditions the Hessian in most cases
!       remains positive definite, if one is far from the minimum
!       and/or "bonds" are being made/broken the curvature condition
!              Hy = s ; or s = By
!       cannot be satisfied if s.y < 0. In addition, if s.y is small
!       compared to s.B.s too greedy a step is taken.
!
!       The trap below is conventional and "OK", and has been around
!       for ~ 30 years but, unfortunately, is rarely mentioned in
!       introductory texts and hence often neglected.
!
!       First, solve for inv_hess*t = s ; i.e. t = B*s
!       Use yH as workspace here

        ALLOCATE (B(n,n) )
        B = inv_hess
        yH= s
        call DPOSV('U',n,1,B,n, yH, n, info)
!       Info .ne. 0 should be trapped ...
        if(info .ne. 0)write( stdout, '(/,5X,"WARNING: info=",i3," for Hessian")' )info
        DEALLOCATE ( B )
!
!       Calculate s.B.s
        sBs = ( s(:) .dot. yH(:) )
!
!       Now the trap itself
        if ( sdoty < 0.20D0*sBs ) then
!               Conventional damping
                Theta = 0.8D0*sBs/(sBs-sdoty)
                WRITE( stdout, '(/,5X,"WARNING: bfgs curvature condition ", &
                &     "failed, Theta=",F6.3)' )theta
                y = Theta*y + (1.D0 - Theta)*yH
        endif
      END IF
      !
      Hy(:) = ( inv_hess .times. y(:) )
      yH(:) = ( y(:) .times. inv_hess )
      !
      ! ... BFGS update
      !
      inv_hess = inv_hess + 1.0_DP / sdoty * &
                 ( ( 1.0_DP + ( y .dot. Hy ) / sdoty ) * matrix( s, s ) - &
                  ( matrix( s, yH ) +  matrix( Hy, s ) ) )
      !
      DEALLOCATE( y, s, Hy, yH )
      !
      RETURN
      !
   END SUBROUTINE update_inverse_hessian
   !
   !------------------------------------------------------------------------
   SUBROUTINE check_wolfe_conditions( lwolfe, energy, grad )
      !------------------------------------------------------------------------
      IMPLICIT NONE
      REAL(DP), INTENT(IN)  :: energy
      REAL(DP), INTENT(IN)  :: grad(:)
      LOGICAL,  INTENT(OUT) :: lwolfe
      !
      lwolfe =  energy_wolfe_condition ( energy ) .AND. &
                gradient_wolfe_condition ( grad )
      !
   END SUBROUTINE check_wolfe_conditions
   !
   !------------------------------------------------------------------------
   LOGICAL FUNCTION energy_wolfe_condition ( energy )
      !------------------------------------------------------------------------
      IMPLICIT NONE
      REAL(DP), INTENT(IN)  :: energy
      !
      energy_wolfe_condition = &
          ( energy-energy_p ) < w_1 * ( grad_p.dot.step_old ) * trust_radius_old
      !
   END FUNCTION energy_wolfe_condition
   !
   !------------------------------------------------------------------------
   LOGICAL FUNCTION gradient_wolfe_condition ( grad )
      !------------------------------------------------------------------------
      IMPLICIT NONE
      REAL(DP), INTENT(IN)  :: grad(:)
      !
      gradient_wolfe_condition = &
          ABS( grad .dot. step_old ) < - w_2 * ( grad_p .dot. step_old )
      !
   END FUNCTION gradient_wolfe_condition
   !
   !------------------------------------------------------------------------
   SUBROUTINE compute_trust_radius( lwolfe, energy, grad, n, lfcp, fcp_hess )
      !-----------------------------------------------------------------------
      !
      IMPLICIT NONE
      !
      LOGICAL,  INTENT(IN)  :: lwolfe
      REAL(DP), INTENT(IN)  :: energy
      REAL(DP), INTENT(IN)  :: grad(:)
      INTEGER,  INTENT(IN)  :: n
      LOGICAL,  INTENT(IN)  :: lfcp
      REAL(DP), INTENT(IN)  :: fcp_hess
      !
      REAL(DP) :: a
      LOGICAL  :: ltest
      !
      ltest = energy_wolfe_condition( energy )
      !
      ! The instruction below replaces the original instruction:
      !    ltest = ltest .AND. ( nr_step_length_old > trust_radius_old )
      ! which gives a random result if trust_radius was set equal to 
      ! nr_step_length at previous step. I am not sure what the best 
      ! action should be in that case, though (PG)
      !
      ltest = ltest .AND. ( nr_step_length_old > trust_radius_old + eps8 )
      !
      IF ( ltest ) THEN
         a = 1.5_DP
      ELSE
         a = 1.1_DP
      END IF
      IF ( lwolfe ) a = 2._DP * a
      !
      trust_radius = MIN( trust_radius_max, a*trust_radius_old, nr_step_length )
      !
      IF ( trust_radius < trust_radius_min ) THEN
         !
         ! ... the history is reset
         !
         ! ... if tr_min_hit the history has already been reset at the 
         ! ... previous step : something is going wrong
         !
         IF ( tr_min_hit == 1 ) THEN
            CALL infomsg( 'bfgs :: compute_trust_radius', &
                          'history already reset at previous step: stopping' )
            tr_min_hit = 2 
         ELSE
            tr_min_hit = 1
         END IF
         !
         WRITE( UNIT = stdout, &
                FMT = '(5X,"small trust_radius: resetting bfgs history",/)' )
         !
         CALL reset_bfgs( n, lfcp, fcp_hess )
         step(:) = - ( inv_hess(:,:) .times. grad(:) )
         !
         nr_step_length = scnorm(step)
         step(:) = step(:) / nr_step_length
         !
         trust_radius = min(trust_radius_min, nr_step_length )
         !
      ELSE
         !
         tr_min_hit = 0
         !
      END IF
      !
   END SUBROUTINE compute_trust_radius
   !
   !----------------------------------------------------------------------- 
   REAL(DP) FUNCTION scnorm1( vect )
      !-----------------------------------------------------------------------
      IMPLICIT NONE
      REAL(DP), INTENT(IN) :: vect(:)
      !
      scnorm1 = SQRT( DOT_PRODUCT( vect  ,  MATMUL( metric, vect ) ) )
      !
   END FUNCTION scnorm1
   !
   !----------------------------------------------------------------------- 
   REAL(DP) FUNCTION scnorm( vect )
      !-----------------------------------------------------------------------
      IMPLICIT NONE
      REAL(DP), INTENT(IN) :: vect(:)
      REAL(DP) :: ss
      INTEGER :: i,k,l,n
      INTEGER :: ndim
      !
      ndim = SIZE (vect)
      !
      scnorm = 0._DP
      n = (ndim - 1) / 3
      do i=1,n
         ss = 0._DP
         do k=1,3
            do l=1,3
               ss = ss + &
                    vect(k+(i-1)*3)*metric(k+(i-1)*3,l+(i-1)*3)*vect(l+(i-1)*3)
            end do
         end do
         scnorm = MAX (scnorm, SQRT (ss) )
      end do
      !
      ss = vect(ndim) * metric(ndim,ndim) * vect(ndim)
      scnorm = MAX (scnorm, SQRT (ss) )
      !
   END FUNCTION scnorm
   !
   !------------------------------------------------------------------------
   SUBROUTINE terminate_bfgs( energy, energy_thr, grad_thr, cell_thr, fcp_thr, &
                              lmovecell, lfcp, failed )
      !------------------------------------------------------------------------
      !
      USE io_files, ONLY : delete_if_present
      !
      IMPLICIT NONE
      REAL(DP),         INTENT(IN) :: energy, energy_thr, grad_thr, cell_thr, fcp_thr
      LOGICAL,          INTENT(IN) :: lmovecell, lfcp, failed
      !
      IF ( conv_bfgs ) THEN
         !
         IF ( failed ) THEN
            WRITE( UNIT = stdout, &
              & FMT = '(/,5X,"bfgs failed after ",I3," scf cycles and ", &
              &         I3," bfgs steps, convergence not achieved")' ) &
              & scf_iter, bfgs_iter
         ELSE
            WRITE( UNIT = stdout, &
              & FMT = '(/,5X,"bfgs converged in ",I3," scf cycles and ", &
              &         I3," bfgs steps")' ) scf_iter, bfgs_iter
         END IF
         IF ( lmovecell ) THEN
            WRITE( UNIT = stdout, &
              & FMT = '(5X,"(criteria: energy < ",ES8.1," Ry, force < ",ES8.1,&
              &       " Ry/Bohr, cell < ",ES8.1," kbar)")') energy_thr, grad_thr, cell_thr
         ELSE
            WRITE( UNIT = stdout, &
              & FMT = '(5X,"(criteria: energy < ",ES8.1," Ry, force < ",ES8.1,&
              &            " Ry/Bohr)")') energy_thr, grad_thr
         END IF
         !
         IF ( lfcp ) THEN
            WRITE( UNIT = stdout, &
              & FMT = '(5X,"(criteria: force on FCP < ",ES8.1," eV)")') fcp_thr * RYTOEV
         END IF
         !
         WRITE( UNIT = stdout, &
              & FMT = '(/,5X,"End of BFGS Geometry Optimization")' )
         WRITE( UNIT = stdout, &
              & FMT = '(/,5X,"Final ",A," = ",F18.10," Ry")' ) fname, energy
         !
         CALL delete_if_present( bfgs_file )
         bfgs_file = " "
         !
      ELSE
         !
         WRITE( UNIT = stdout, &
                FMT = '(/,5X,"The maximum number of steps has been reached.")' )
         WRITE( UNIT = stdout, &
                FMT = '(/,5X,"End of BFGS Geometry Optimization")' )
         !
      END IF
      !
   END SUBROUTINE terminate_bfgs
   !
   FUNCTION bfgs_get_n_iter (what)  RESULT(n_iter)
   !  
   IMPLICIT NONE
   INTEGER                         :: n_iter
   CHARACTER(10),INTENT(IN)        :: what
   SELECT CASE (TRIM(what)) 
      CASE ('bfgs_iter') 
           n_iter = bfgs_iter
      CASE ( 'scf_iter') 
           n_iter = scf_iter
      CASE default 
           n_iter = -1
   END SELECT
   END FUNCTION bfgs_get_n_iter
END MODULE bfgs_module
