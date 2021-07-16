!
! Copyright (C) 2001-2016 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
SUBROUTINE lr_setup_nscf ()
  !---------------------------------------------------------------------
  !
  ! This subroutine initializes variables for the non-scf calculations at k
  ! and k+q required by the linear response calculation at finite q.
  ! Determines k and k+q points in the irreducible BZ.
  ! Needed on input (read from data file):
  ! "nkstot" k-points in the irreducible BZ wrt lattice symmetry.
  ! Inspired by PH/set_defaults_pw.f90
  !
  ! Created by Iurii Timrov (2013)
  !
  USE kinds,              ONLY : DP
  USE lr_variables,       ONLY : magnons
  USE constants,          ONLY : eps8
  USE parameters,         ONLY : npk
  USE constants,          ONLY : pi, degspin
  USE cell_base,          ONLY : at, bg, alat, tpiba, tpiba2, ibrav, omega
  USE ions_base,          ONLY : nat, tau, ityp, zv
  USE force_mod,          ONLY : force
  USE basis,              ONLY : natomwfc
  USE klist,              ONLY : xk, wk, nks, nelec, degauss, lgauss, &
                                 nkstot, qnorm
  USE lsda_mod,           ONLY : lsda, nspin, current_spin, isk
  USE symm_base,          ONLY : s, t_rev, nrot, nsym, time_reversal
  USE wvfct,              ONLY : nbnd, nbndx
  USE control_flags,      ONLY : ethr, isolve, david, use_para_diag, &
                                 & noinv, max_cg_iter
  USE control_lr,         ONLY : ethr_nscf
  USE mp_pools,           ONLY : kunit
  USE spin_orb,           ONLY : domag
  USE noncollin_module,   ONLY : noncolin
  USE start_k,            ONLY : nks_start, xk_start, wk_start, &
                                 nk1, nk2, nk3, k1, k2, k3
  USE upf_ions,           ONLY : n_atom_wfc
  USE lr_symm_base,       ONLY : nsymq, minus_q
  USE qpoint,             ONLY : xq
  !
  USE io_global,          ONLY : ionode
  USE io_files,           ONLY : prefix
  ! 
  IMPLICIT NONE
  !
  LOGICAL :: magnetic_sym 
  !
  CALL start_clock( 'lr_setup_nscf' )
  ! 
  IF ( .NOT. ALLOCATED( force ) ) ALLOCATE( force( 3, nat ) )
  !
  ! ... threshold for diagonalization ethr
  !
  ethr = ethr_nscf
  !
  ! ... variables for iterative diagonalization (Davidson is assumed)
  !
  isolve = 0
  david  = 4
  nbndx  = david*nbnd
  max_cg_iter = 20
  natomwfc = n_atom_wfc( nat, ityp, noncolin )
  !
  CALL set_para_diag( nbnd, use_para_diag )
  !
  ! Symmetry section
  !
  ! time_reversal = use q=>-q symmetry for k-point generation
  !
  magnetic_sym = noncolin .AND. domag
  !
  !time_reversal = .NOT. noinv .AND. .NOT. magnetic_sym
  !
  ! Determine the small group of q : S q = q + G
  !
  CALL lr_smallgq (xq)
  !
  !  K points section   
  !
  ! Input k-points are assumed to be given in the IBZ of the Bravais
  ! lattice, with the full point symmetry of the lattice.
  !
  IF ( magnons ) THEN
     !
     ! Generate a new set of k-points for magnons
     !
     ! In case of magnetic response we need to couple +k and -k
     !
     CALL kpoint_grid_no_t_rev (bg, npk, k1,k2,k3, nk1,nk2,nk3, nkstot, xk, wk)
     !
     ! Add k+Q and k-Q to the list of k-points 
     ! need irreducible_BZ???
     ! 
     CALL set_kplusq_kminusq( xk, wk, xq, nkstot, npk )
     !
  ELSE
     !
     IF ( nks_start > 0 ) THEN
        !
        ! In this case keep the same points of the charge-density calculation
        !
        nkstot = nks_start
        xk(:,1:nkstot) = xk_start(:,1:nkstot)
        wk(1:nkstot)   = wk_start(1:nkstot)
        !
     ELSE
        !
        ! Generate a new set of k-points
        !
        CALL kpoint_grid ( nrot, time_reversal,.false., s, t_rev, &
                         bg, nk1*nk2*nk3, k1,k2,k3, nk1,nk2,nk3, nkstot, xk, wk)
        !
     ENDIF
     !
     ! If some symmetries of the lattice are missing in the crystal,
     ! "irreducible_BZ" computes the missing k-points.
     !
     CALL irreducible_BZ (nrot, s, nsymq, minus_q, magnetic_sym, &
                          at, bg, npk, nkstot, xk, wk, t_rev)
     !
     ! Add k+q to the list of k
     !
     CALL set_kplusq( xk, wk, xq, nkstot, npk )
     !
  ENDIF
  !
  IF ( lsda ) THEN
     !
     ! ... LSDA case: two different spin polarizations,
     ! ...            each with its own kpoints
     !
     IF (nspin /= 2) CALL errore ('lr_setup_nscf','nspin should be 2; check iosys',1)
     !
     CALL set_kup_and_kdw( xk, wk, isk, nkstot, npk )
     !
  ELSEIF ( noncolin ) THEN
     !
     ! ... noncolinear magnetism: potential and charge have dimension 4 (1+3)
     !
     IF (nspin /= 4) CALL errore ('lr_setup_nscf','nspin should be 4; check iosys',1)
     current_spin = 1
     !
  ELSE
     !
     ! ... LDA case: the two spin polarizations are identical
     !
     wk(1:nkstot) = wk(1:nkstot) * degspin
     current_spin = 1
     !
     IF ( nspin /= 1 ) &
        CALL errore( 'lr_setup_nscf', 'nspin should be 1; check iosys', 1 )
     !
  ENDIF
  !
  IF ( nkstot > npk ) CALL errore( 'lr_setup_nscf', 'too many k points', nkstot )
  !
  ! ...notice: qnorm is used by allocate_nlpot to determine
  ! the correct size of the interpolation table "qrad"
  !
  qnorm = sqrt(xq(1)**2 + xq(2)**2 + xq(3)**2) * tpiba
  !
  ! ... set the granularity for k-point distribution
  !
  IF ( ABS( xq(1) ) < eps8 .AND. ABS( xq(2) ) < eps8 .AND. &
       ABS( xq(3) ) < eps8 ) THEN
     !
     kunit = 1
     !
  ELSEIF (magnons) THEN
     !
     ! In the magnons case we need
     !
     !  k,  k+Q,  k-Q
     ! -k, -k+Q, -k-Q
     !
     ! in the same pool, therefore kunit=6
     !
     kunit = 6
     !
  ELSE
     !
     kunit = 2
     !
  ENDIF
  !
  ! ... distribute k-points (and their weights and spin indices)
  !
  CALL divide_et_impera( nkstot, xk, wk, isk, nks )
  !
  CALL stop_clock( 'lr_setup_nscf' )
  !
  RETURN
  CONTAINS
   !
   SUBROUTINE kpoint_grid_no_t_rev ( bg, npk, k1,k2,k3, nk1,nk2,nk3, nks, xk, wk)
     !
     USE kinds, ONLY: DP
     USE io_global,  ONLY : stdout  
     !
     IMPLICIT NONE
     !
     INTEGER, INTENT(in):: npk, k1, k2, k3, nk1, nk2, nk3
     real(DP), INTENT(in):: bg(3,3)
     !
     INTEGER, INTENT(out) :: nks
     real(DP), INTENT(out):: xk(3,npk)
     real(DP), INTENT(out):: wk(npk)
     ! LOCAL:
     real(DP), PARAMETER :: eps=1.0d-5
     real(DP) :: xkr(3), fact, xx, yy, zz
     real(DP), ALLOCATABLE:: xkg(:,:), wkk(:)
     INTEGER :: nkr, i,j,k, ns, n, nk
     INTEGER, ALLOCATABLE :: equiv(:)
     LOGICAL :: in_the_list
     LOGICAL :: is_gamma

     CHARACTER(len=256) :: filename

     !
     nkr=nk1*nk2*nk3
     ALLOCATE (xkg( 3,nkr))
     ALLOCATE (equiv( nkr))
     !
     DO i=1,nk1
        DO j=1,nk2
           DO k=1,nk3
              !  this is nothing but consecutive ordering
              n = (k-1) + (j-1)*nk3 + (i-1)*nk2*nk3 + 1
              !  xkg are the components of the complete grid in crystal axis
              xkg(1,n) = dble(i-1)/nk1 + dble(k1)/2/nk1
              xkg(2,n) = dble(j-1)/nk2 + dble(k2)/2/nk2
              xkg(3,n) = dble(k-1)/nk3 + dble(k3)/2/nk3
           ENDDO
        ENDDO
     ENDDO

     !if (ionode) write(stdout,*) "nkr =", nkr


     DO nk=1,nkr
        equiv(nk)=nk
     ENDDO

     DO nk=1,nkr
     !  check if this k-point has already been found to be minus another
       IF (equiv(nk) == nk) THEN
         !  check if there are equivalent k-point to this in the list
         !  (excepted those previously found to be equivalent to another)
         !  check both k and -k
         DO i = 1, 3
           xkr(i) = -xkg(i,nk) + nint( xkg(i,nk) )
         ENDDO
         !
         xx = xkr(1)*nk1 - 0.5d0*k1
         yy = xkr(2)*nk2 - 0.5d0*k2
         zz = xkr(3)*nk3 - 0.5d0*k3
         in_the_list = abs(xx-nint(xx))<=eps .and. &
                       abs(yy-nint(yy))<=eps .and. &
                       abs(zz-nint(zz))<=eps
         IF (in_the_list) THEN
            i = mod ( nint ( xkr(1)*nk1 - 0.5d0*k1 + 2*nk1), nk1 ) + 1
            j = mod ( nint ( xkr(2)*nk2 - 0.5d0*k2 + 2*nk2), nk2 ) + 1
            k = mod ( nint ( xkr(3)*nk3 - 0.5d0*k3 + 2*nk3), nk3 ) + 1
            n = (k-1) + (j-1)*nk3 + (i-1)*nk2*nk3 + 1
            IF (n>nk .and. equiv(n)==n) THEN
               equiv(n) = nk
            ELSE
               IF (equiv(n)/=nk .or. n<nk ) CALL errore('kpoint_grid', &
                  'something wrong in the checking algorithm',1)
            ENDIF
         ENDIF
       ENDIF
     ENDDO

     !  count irreducible points and order them

     nks=0
     fact=0.0d0
     DO nk=1,nkr
        IF (equiv(nk)==nk) THEN
           nks=nks+2
           is_gamma = abs(xkg(1,nk)-nint(xkg(1,nk)))<=eps .and. &
                      abs(xkg(2,nk)-nint(xkg(2,nk)))<=eps .and. &
                      abs(xkg(3,nk)-nint(xkg(3,nk)))<=eps
           IF (is_gamma) THEN
             wk(nks-1) = 0.5
             wk(nks) = 0.5
           ELSE
             wk(nks-1) = 1.0
             wk(nks) = 1.0 
           ENDIF
           fact = fact+wk(nks-1)+wk(nks)
           !
           !IF (nks>npk) CALL errore('kpoint_grid','too many k-points',1)
           !  bring back into to the first BZ
           DO i=1,3
              xk(i,nks-1) = xkg(i,nk)-nint(xkg(i,nk))
              xk(i,nks) = -xk(i,nks-1) 
           ENDDO
        ENDIF
     ENDDO

     !if (ionode) write(stdout,*) "nks =", nks

     !  go to cartesian axis (in units 2pi/a0)
     CALL cryst_to_cart(nks,xk,bg,1)
     !
     !  normalize weights to one
     DO nk=1,nks
        wk(nk) = wk(nk)/fact
     ENDDO
     !
     !filename = trim(prefix) // '__k_mesh.xyz' 
     !if (ionode) OPEN (159, file = filename, form = 'formatted', status = 'unknown')
     !if (ionode) write(159,*) nks
     !if (ionode) write(159,*) ""
     !  normalize weights to one
     !DO nk=1,nks
     !   if (ionode .and. mod(nk,2) == 0 ) then
     !     write(159,'("C",3f16.9)') xk(1,nk-1), xk(2,nk-1), xk(3,nk-1)
     !     write(159,'("F",3f16.9)')  xk(1,nk), xk(2,nk), xk(3,nk)
     !   endif
     !ENDDO
     !if (ionode) close(159)

     DEALLOCATE(equiv)
     DEALLOCATE(xkg)

     RETURN
   END SUBROUTINE kpoint_grid_no_t_rev
   !
   !
   SUBROUTINE set_kplusq_kminusq (xk, wk, xq, nks, npk)
     !
     USE kinds, only : DP
     implicit none
     !
     !    First the dummy variables
     !
   
     integer :: npk, nks
     ! input-output: maximum allowed number of k
     ! input-output: starting and ending number of
     real(DP) :: xk (3, npk), wk (npk), eps, xq (3)
     ! input-output: coordinates of k points
     ! input-output: weights of k points
     ! the smallest xq
     ! input: coordinates of a q-point
     !
     !    And then the local variables
     !
   
     logical :: lgamma
     ! true if xq is the gamma point
     integer :: ik, j
     ! counter on k
     ! counter
     !
     eps = 1.d-12
     !
     !
     lgamma = abs (xq (1) ) .lt.eps.and.abs (xq (2) ) .lt.eps.and.abs ( &
          xq (3) ) .lt.eps
   
     if (.not.lgamma) then
        !!
        !! in magnons the limit npk should maybe be increased
        !!
        if (3 * nks.gt. npk) call errore ('set_kplusq', 'too many k points', &
             & nks)
        do ik = nks, 1, - 1
           do j = 1, 3
              xk (j, 3 * ik - 2) = xk (j, ik)
              xk (j, 3 * ik - 1) = xk (j, ik) + xq (j)
              xk (j, 3 * ik) = xk (j, ik) - xq (j)
           enddo
           wk (3 * ik - 2) = wk (ik)
           wk (3 * ik - 1) = 0.0d0
           wk (3 * ik) = 0.d0
        enddo
        nks = 3 * nks
   
     endif
     return
   END SUBROUTINE set_kplusq_kminusq
  !
END SUBROUTINE lr_setup_nscf
