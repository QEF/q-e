!
! Copyright (C) 2008-2009 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------------------------------------------
!
MODULE symme
  !
  USE kinds,      ONLY : DP
  !
  ! ... The variables needed to describe the symmetry properties
  !  
  SAVE
  !
  PRIVATE
  !
  PUBLIC :: s, sname, ftau, nrot, nsym, t_rev, time_reversal, irt, &
            invs, invsym, d1, d2, d3
  !
  INTEGER :: &
       s(3,3,48),            &! simmetry matrices, in crystal axis
       invs(48),             &! index of inverse operation: S^{-1}_i=S(invs(i))
       ftau(3,48),           &! fractional translations, in FFT coordinates
       nrot,                 &! number of bravais lattice symmetries 
       nsym                   ! number of crystal symmetries
  INTEGER :: &
       t_rev(48) = 0          ! time reversal flag, for noncolinear magnetisation
  INTEGER, ALLOCATABLE :: &
       irt(:,:)               ! symmetric atom for each atom and sym.op.
  LOGICAL :: &
       time_reversal=.true., &! if .TRUE. the system has time_reversal symmetry
       invsym                 ! if .TRUE. the system has inversion symmetry
  REAL(DP),TARGET :: &
       d1(3,3,48),           &! matrices for rotating spherical
       d2(5,5,48),           &! harmonics (d1 for l=1, ...)
       d3(7,7,48)             !
  CHARACTER(LEN=45) ::  sname(48)   ! name of the symmetries
  !
  ! For symmetrization in reciprocal space (all variables are private)
  !
  PUBLIC :: sym_rho_init, sym_rho, sym_rho_deallocate, inverse_s
  !
  LOGICAL :: &
       no_rho_sym=.true.      ! do not perform symetrization of charge density
  INTEGER :: ngs              ! number of symmetry-related G-vector shells
  TYPE shell_type
     INTEGER, POINTER :: vect(:)
  END TYPE shell_type
  ! shell contains a list of symmetry-related G-vectors for each shell
  TYPE(shell_type), ALLOCATABLE :: shell(:)
  ! Arrays used for parallel symmetrization
  INTEGER, ALLOCATABLE :: sendcnt(:), recvcnt(:), sdispls(:), rdispls(:)
  !
CONTAINS
   !
   LOGICAL FUNCTION rho_sym_needed ( )
      rho_sym_needed = .NOT. no_rho_sym
   END FUNCTION rho_sym_needed
   !
   SUBROUTINE inverse_s ( )
     !
     ! Locate index of S^{-1}
     !
     IMPLICIT NONE
     !
     INTEGER :: isym, jsym, ss (3, 3)
     LOGICAL :: found
     !
     DO isym = 1, nsym
        found = .FALSE.
        DO jsym = 1, nsym
           !
           ss = MATMUL (s(:,:,jsym),s(:,:,isym))
           ! s(:,:,1) is the identity
           IF ( ALL ( s(:,:,1) == ss(:,:) ) ) THEN
              invs (isym) = jsym
              found = .TRUE.
           END IF
        END DO
        IF ( .NOT.found) CALL errore ('inverse_s', ' Not a group', 1)
     END DO
     !
END SUBROUTINE inverse_s 
   !
   SUBROUTINE sym_rho_init ( gamma_only )
    !-----------------------------------------------------------------------
    !
    !  Initialize arrays needed for symmetrization in reciprocal space
    ! 
    USE gvect, ONLY : ngm, g
    !
    LOGICAL, INTENT(IN) :: gamma_only
    !
    no_rho_sym = gamma_only .OR. (nsym==1)
    IF (no_rho_sym) RETURN
#ifdef __PARA
    CALL sym_rho_init_para ( )
#else
    CALL sym_rho_init_shells( ngm, g )
#endif
    !
  END SUBROUTINE sym_rho_init
   !
#ifdef __PARA
  !
  SUBROUTINE sym_rho_init_para ( )
    !-----------------------------------------------------------------------
    !
    !  Initialize arrays needed for parallel symmetrization
    ! 
    USE mp_global, ONLY : nproc_pool, me_pool, intra_pool_comm
    USE parallel_include
    USE gvect, ONLY : ngm, gcutm, g, gg
    !
    IMPLICIT NONE
    !
    REAL(DP), PARAMETER :: twothirds = 0.6666666666666666_dp
    REAL(DP), ALLOCATABLE :: gcut_(:), g_(:,:)
    INTEGER :: np, ig, ngloc, ngpos, ierr, ngm_
    !
    ALLOCATE ( sendcnt(nproc_pool), recvcnt(nproc_pool), &
               sdispls(nproc_pool), rdispls(nproc_pool) )
    ALLOCATE ( gcut_(nproc_pool) )
    !
    ! the gcut_ cutoffs are estimated in such a way that there is an similar
    ! number of G-vectors in each shell gcut_(i) < G^2 < gcut_(i+1)
    !
    DO np = 1, nproc_pool
       gcut_(np) = gcutm * np**twothirds/nproc_pool**twothirds
    END DO
    !
    ! find the number of G-vectors in each shell (defined as above)
    ! beware: will work only if G-vectors are in order of increasing |G|
    !
    ngpos=0
    DO np = 1, nproc_pool
       sdispls(np) = ngpos
       ngloc=0
       DO ig=ngpos+1,ngm
          IF ( gg(ig) > gcut_(np) ) EXIT
          ngloc = ngloc+1
       END DO
       IF ( ngloc < 1 ) CALL infomsg('sym_rho_init', &
            'likely internal error: no G-vectors found')
       sendcnt(np) = ngloc
       ngpos = ngpos + ngloc
       IF ( ngpos > ngm ) &
            CALL errore('sym_rho','internal error: too many G-vectors', ngpos)
    END DO
    IF ( ngpos /= ngm .OR. ngpos /= SUM (sendcnt)) &
         CALL errore('sym_rho_init', &
         'internal error: inconsistent number of G-vectors', ngpos)
    DEALLOCATE ( gcut_ )
  !
  ! sendcnt(i) = n_j(i) = number of G-vectors in shell i for processor j (this)
  ! sdispls(i) = \sum_{k=1}^i n_j(k) = starting position of shell i for proc j
  ! we need the number and positions of G-vector shells for other processors
  !
    CALL mpi_alltoall( sendcnt, 1, MPI_INTEGER, recvcnt, 1, MPI_INTEGER, &
         intra_pool_comm, ierr)
    !
    rdispls(1) = 0
    DO np = 2, nproc_pool
       rdispls(np) = rdispls(np-1)+ recvcnt(np-1)
    END DO
    !
    ! recvcnt(i) = n_i(j) = number of G-vectors in shell j for processor i
    ! rdispls(i) = \sum_{k=1}^i n_k(j) = start.pos. of shell j for proc i
    !
    ! now collect G-vector shells on each processor
    !
    ngm_ = SUM(recvcnt)
    ALLOCATE (g_(3,ngm_))
    ! remember that G-vectors have 3 components
    sendcnt(:) = 3*sendcnt(:)
    recvcnt(:) = 3*recvcnt(:)
    sdispls(:) = 3*sdispls(:)
    rdispls(:) = 3*rdispls(:)
    CALL mpi_alltoallv ( g , sendcnt, sdispls, MPI_DOUBLE_PRECISION, &
                         g_, recvcnt, rdispls, MPI_DOUBLE_PRECISION, &
                         intra_pool_comm, ierr)
    sendcnt(:) = sendcnt(:)/3
    recvcnt(:) = recvcnt(:)/3
    sdispls(:) = sdispls(:)/3
    rdispls(:) = rdispls(:)/3
    !
    ! find shells of symetry-related G-vectors
    !
    CALL sym_rho_init_shells( ngm_, g_ )
    !
    DEALLOCATE (g_)
    !
  END SUBROUTINE sym_rho_init_para
  !
#endif
  !
  SUBROUTINE sym_rho_init_shells ( ngm_, g_ )
    !-----------------------------------------------------------------------
    !
    !  Initialize G-vector shells needed for symmetrization
    ! 
    !
    USE cell_base, ONLY : at
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN) :: ngm_
    REAL(DP), INTENT(IN) :: g_(3,ngm_)
    !
    LOGICAL, ALLOCATABLE :: done(:)
    INTEGER, ALLOCATABLE :: n(:,:)
    INTEGER :: i,j,is,ig, ng, sn(3), gshell(3,48)
    LOGICAL :: found
    !
    ngs = 0
    ! shell should be allocated to the number of symmetry shells
    ! since this is unknown, we use the number of all G-vectors
    ALLOCATE ( shell(ngm_) )
    ALLOCATE ( done(ngm_), n(3,ngm_) )
    DO ig=1,ngm_
       !
       done(ig) = .false.
       ! G-vectors are stored as integer indices in crystallographic axis:
       !    G = n(1)*at(1) + n(2)*at(2) + n(3)*at(3)
       n(:,ig) = nint ( at(1,:)*g_(1,ig) + at(2,:)*g_(2,ig) + at(3,:)*g_(3,ig) )
       !
    END DO
    !
    DO ig=1,ngm_
       !
       IF ( done(ig) ) CYCLE
       !
       ! we start a new shell of symmetry-equivalent G-vectors
       ngs = ngs+1
       ! ng: counter on G-vectors in this shell
       ng  = 0
       DO is=1,nsym
          ! integer indices for rotated G-vector
          sn(:)=s(:,1,is)*n(1,ig)+s(:,2,is)*n(2,ig)+s(:,3,is)*n(3,ig)
          found = .false.
          ! check if this rotated G-vector is equivalent to any other
          ! vector already present in this shell
shelloop: DO i=1,ng
             found = ( sn(1)==gshell(1,i) .and. &
                       sn(2)==gshell(2,i) .and. &
                       sn(3)==gshell(3,i) )
             if (found) exit shelloop
          END DO shelloop
          IF ( .not. found ) THEN
             ! add rotated G-vector to this shell
             ng = ng + 1
             IF (ng > 48) CALL errore('sym_rho_init_shell','internal error',48)
             gshell(:,ng) = sn(:)
          END IF
       END DO
       ! there are ng vectors gshell in shell ngs
       ! now we have to locate them in the list of G-vectors
       ALLOCATE ( shell(ngs)%vect(ng))
       DO i=1,ng
gloop:    DO j=ig,ngm_
             IF (done(j)) CYCLE gloop
                found = ( gshell(1,i)==n(1,j) .and. &
                          gshell(2,i)==n(2,j) .and. &
                          gshell(3,i)==n(3,j) )
             IF ( found ) THEN
                done(j)=.true.
                shell(ngs)%vect(i) = j
                EXIT gloop
             END IF
          END DO gloop
          IF (.not. found) CALL errore('sym_rho_init_shell','lone vector',i)
       END DO
       !
    END DO
    DEALLOCATE ( n, done ) 

  END SUBROUTINE sym_rho_init_shells
  !
  !-----------------------------------------------------------------------
  SUBROUTINE sym_rho (nspin, rhog)
    !-----------------------------------------------------------------------
    !
    !     Symmetrize the charge density rho in reciprocal space
    !     Distributed parallel algorithm: collects entire shells of G-vectors
    !     and corresponding rho(G), calls sym_rho_serial to perform the
    !     symmetrization, re-distributed rho(G) into original ordering
    !     rhog(ngm,nspin) components of rho: rhog(ig) = rho(G(:,ig))
    !                     unsymmetrized on input, symmetrized on output
    !     nspin=1,2,4     unpolarized, LSDA, non-colinear magnetism     
    !
    USE constants,            ONLY : eps8, eps6
    USE gvect,                ONLY : ngm, g
#ifdef __PARA
    USE parallel_include
    USE mp_global,            ONLY : nproc_pool, me_pool, intra_pool_comm
#endif
!	use cell_base
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN) :: nspin
    COMPLEX(DP), INTENT(INOUT) :: rhog(ngm,nspin)
    !
    REAL(DP), allocatable :: g0(:,:), g_(:,:), gg_(:) 
    REAL(DP) :: gg0_, gg1_
    COMPLEX(DP), allocatable :: rhog_(:,:)
    INTEGER :: is, ig, igl, np, ierr, ngm_
    !
    IF ( no_rho_sym) RETURN
#ifndef __PARA
    !
    CALL sym_rho_serial ( ngm, g, nspin, rhog )
    !
#else
    !
    ! we transpose the matrix of G-vectors and their coefficients
    !
    ngm_ = SUM(recvcnt)
    ALLOCATE (rhog_(ngm_,nspin),g_(3,ngm_))
    DO is=1,nspin
       CALL mpi_alltoallv (rhog (1,is) , sendcnt, sdispls, MPI_DOUBLE_COMPLEX,&
            rhog_(1,is), recvcnt, rdispls, MPI_DOUBLE_COMPLEX, &
            intra_pool_comm, ierr)
    END DO
    ! remember that G-vectors have 3 components
    sendcnt(:) = 3*sendcnt(:)
    recvcnt(:) = 3*recvcnt(:)
    sdispls(:) = 3*sdispls(:)
    rdispls(:) = 3*rdispls(:)
    CALL mpi_alltoallv ( g , sendcnt, sdispls, MPI_DOUBLE_PRECISION, &
         g_, recvcnt, rdispls, MPI_DOUBLE_PRECISION, &
         intra_pool_comm, ierr)
    !
    !   Now symmetrize
    !
    CALL sym_rho_serial ( ngm_, g_, nspin, rhog_ )
    !
    DEALLOCATE ( g_ )
    !
    ! bring symmetrized rho(G) back to original distributed form
    !
    sendcnt(:) = sendcnt(:)/3
    recvcnt(:) = recvcnt(:)/3
    sdispls(:) = sdispls(:)/3
    rdispls(:) = rdispls(:)/3
    DO is = 1, nspin
       CALL mpi_alltoallv (rhog_(1,is), recvcnt, rdispls, MPI_DOUBLE_COMPLEX, &
            rhog (1,is), sendcnt, sdispls, MPI_DOUBLE_COMPLEX, &
            intra_pool_comm, ierr)
    END DO
    DEALLOCATE ( rhog_ )
#endif
    !
    RETURN
  END SUBROUTINE sym_rho
  !
  !-----------------------------------------------------------------------
  SUBROUTINE sym_rho_serial ( ngm_, g_, nspin_, rhog_ )
    !-----------------------------------------------------------------------
    !
    !     symmetrize the charge density rho in reciprocal space 
    !     Serial algorithm - requires in input: 
    !     g_(3,ngm_)      list of G-vectors
    !     nspin_          number of spin components to be symmetrized
    !     rhog_(ngm_,nspin_) rho in reciprocal space: rhog_(ig) = rho(G(:,ig))
    !                      unsymmetrized on input, symmetrized on output
    !
    USE kinds
    USE constants,            ONLY : tpi
    USE cell_base,            ONLY : at, bg
    USE gvect,                ONLY : nr1,nr2,nr3
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT (IN) :: ngm_, nspin_
    REAL(DP) , INTENT (IN) :: g_( 3, ngm_ )
    COMPLEX(DP) , INTENT (INOUT) :: rhog_( ngm_, nspin_ )
    !
    REAL(DP), ALLOCATABLE :: g0(:,:)
    REAL(DP) :: sg(3), ft(3,48), arg
    COMPLEX(DP) :: fact, rhosum(2), mag(3), magrot(3), magsum(3)
    INTEGER :: irot(48), ig, isg, igl, ng, ns, nspin_lsda, is
    LOGICAL, ALLOCATABLE :: done(:)
    LOGICAL :: non_symmorphic(48)
    !
    ! convert fractional translations to a.u.
    !
    DO ns=1,nsym
       non_symmorphic(ns) = (ftau(1,ns)/=0 .OR. ftau(2,ns)/=0 .OR. ftau(3,ns)/=0)
       IF ( non_symmorphic(ns) ) ft(:,ns) = at(:,1)*ftau(1,ns)/nr1 + &
                                            at(:,2)*ftau(2,ns)/nr2 + &
                                            at(:,3)*ftau(3,ns)/nr3
    END DO
    !
    IF ( nspin_ == 4 ) THEN
       nspin_lsda = 1
       ! S^{-1} are needed as well
       call inverse_s ( )
       !
    ELSE IF ( nspin_ == 1 .OR. nspin_ == 2 ) THEN
       nspin_lsda = nspin_
    ELSE
       CALL errore('sym_rho_serial','incorrect value of nspin',nspin_)
    END IF
    !
    ! scan shells of G-vectors
    !
    DO igl=1, ngs
       !
       ! symmetrize: \rho_sym(G) = \sum_S rho(SG) for all G-vectors in the star
       !
       ng = SIZE ( shell(igl)%vect )
       allocate ( g0(3,ng), done(ng) )
       IF ( ng < 1 ) CALL errore('sym_rho_serial','internal error',1)
       !
       !  bring G-vectors to crystal axis
       !
       DO ig=1,ng
          g0(:,ig) = g_(:,shell(igl)%vect(ig) )
       END DO
       CALL cryst_to_cart (ng, g0, at,-1)
       !
       !  rotate G-vectors
       !
       done(1:ng) = .false.
       DO ig=1,ng
          IF ( .NOT. done(ig)) THEN
             rhosum(:) = (0.0_dp, 0.0_dp)
             magsum(:) = (0.0_dp, 0.0_dp)
             DO ns=1,nsym
                sg(:)=s(:,1,ns)*g0(1,ig)+s(:,2,ns)*g0(2,ig)+s(:,3,ns)*g0(3,ig)
                irot(ns) = 0
                DO isg=1,ng
                   IF ( ABS ( sg(1)-g0(1,isg) ) < 1.0D-5 .AND. &
                        ABS ( sg(2)-g0(2,isg) ) < 1.0D-5 .AND. &
                        ABS ( sg(3)-g0(3,isg) ) < 1.0D-5 ) THEN
                      irot(ns) = isg
                      EXIT
                   END IF
                END DO
                IF ( irot(ns) < 1 .OR. irot(ns) > ng ) &
                     CALL errore('sym_rho_serial','internal error',2)
                ! isg is the index of rotated G-vector
                isg = shell(igl)%vect(irot(ns))
                !
                ! non-spin-polarized case: component 1 is the charge
                ! LSDA case: components 1,2 are spin-up and spin-down charge
                ! non colinear case: component  1 is the charge density,
                !                    components 2,3,4 are the magnetization
                ! non colinear case: components 2,3,4 are the magnetization
                !
                IF ( nspin_ == 4 ) THEN
                   ! bring magnetization to crystal axis
                   mag(:) = rhog_(isg, 2) * bg(1,:) + &
                            rhog_(isg, 3) * bg(2,:) + &
                            rhog_(isg, 4) * bg(3,:)
                   ! rotate and add magnetization
                   magrot(:) = s(1,:,ns) * mag(1) + &
                               s(2,:,ns) * mag(2) + &
                               s(3,:,ns) * mag(3)
                   IF (sname(ns)(1:3)=='inv') magrot(:)=-magrot(:)
                   IF (t_rev(ns) == 1)        magrot(:)=-magrot(:)
                END IF
                IF ( non_symmorphic (ns) )  THEN
                   arg = tpi * ( g_(1,isg) * ft(1,ns) + &
                                 g_(2,isg) * ft(2,ns) + &
                                 g_(3,isg) * ft(3,ns) )
                   fact = CMPLX ( COS(arg), -SIN(arg), KIND=dp )
                   DO is=1,nspin_lsda
                      rhosum(is) = rhosum(is) + rhog_(isg, is) * fact
                   END DO
                   IF ( nspin_ == 4 ) &
                        magsum(:) = magsum(:) + magrot(:) * fact
                ELSE
                   DO is=1,nspin_lsda
                      rhosum(is) = rhosum(is) + rhog_(isg, is)
                   END DO
                   IF ( nspin_ == 4 ) &
                        magsum(:) = magsum(:) + magrot(:)
                END IF
             END DO
             !
             DO is=1,nspin_lsda
                rhosum(is) = rhosum(is) / nsym
             END DO
             IF ( nspin_ == 4 ) magsum(:) = magsum(:) / nsym
             !
             !  now fill the shell of G-vectors with the symmetrized value
             !
             DO ns=1,nsym
                isg = shell(igl)%vect(irot(ns))
                IF ( nspin_ == 4 ) THEN
                   ! rotate magnetization
                   magrot(:) = s(1,:,invs(ns)) * magsum(1) + &
                               s(2,:,invs(ns)) * magsum(2) + &
                               s(3,:,invs(ns)) * magsum(3)
                   IF (sname(ns)(1:3)=='inv') magrot(:)=-magrot(:)
                   IF (t_rev(ns) == 1)        magrot(:)=-magrot(:)
                   ! back to cartesian coordinates
                   mag(:) = magrot(1) * at(:,1) + &
                            magrot(2) * at(:,2) + &
                            magrot(3) * at(:,3)
                END IF
                IF ( non_symmorphic (ns) )  THEN
                   arg = tpi * ( g_(1,isg) * ft(1,ns) + &
                                 g_(2,isg) * ft(2,ns) + &
                                 g_(3,isg) * ft(3,ns) )
                   fact = CMPLX ( COS(arg), SIN(arg), KIND=dp )
                   DO is=1,nspin_lsda
                      rhog_(isg,is) = rhosum(is) * fact
                   END DO
                   IF ( nspin_ == 4 ) THEN
                      DO is=2,nspin_
                         rhog_(isg, is) = mag(is-1)*fact
                      END DO
                   END IF
                ELSE
                   DO is=1,nspin_lsda
                      rhog_(isg,is) = rhosum(is)
                   END DO
                   IF ( nspin_ == 4 ) THEN
                      DO is=2,nspin_
                         rhog_(isg, is) = mag(is-1)
                      END DO
                   END IF
                END IF
                done(irot(ns)) =.true.
             END DO
          END IF
       END DO
       DEALLOCATE ( done, g0 )
    END DO
    !
    RETURN
  END SUBROUTINE sym_rho_serial

  SUBROUTINE sym_rho_deallocate ( )
    !
    IF ( ALLOCATED (rdispls) ) DEALLOCATE (rdispls) 
    IF ( ALLOCATED (recvcnt) ) DEALLOCATE (recvcnt) 
    IF ( ALLOCATED (sdispls) ) DEALLOCATE (sdispls) 
    IF ( ALLOCATED (sendcnt) ) DEALLOCATE (sendcnt) 
    IF ( ALLOCATED (shell) ) THEN
       DO i=1,SIZE(shell)
          IF ( ASSOCIATED(shell(i)%vect) ) DEALLOCATE (shell(i)%vect)
       END DO
       DEALLOCATE (shell)
    END IF
    !
  END SUBROUTINE sym_rho_deallocate
  !
END MODULE symme
