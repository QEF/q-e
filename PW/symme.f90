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
           invsym, d1, d2, d3
  PUBLIC :: sym_rho_init, sym_rho, sym_rho_deallocate
  !
  INTEGER :: &
       s(3,3,48),            &! simmetry matrices, in crystal axis
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
  ! For parallel symmetrization in reciprocal space
  !
  INTEGER, ALLOCATABLE :: sendcnt(:), recvcnt(:), sdispls(:), rdispls(:)
  !
CONTAINS
  !
  SUBROUTINE sym_rho_init ( )
    !-----------------------------------------------------------------------
    !
    !  Initialize arrays needed for parallel symmetrization in reciprocal space
    ! 
#ifdef __PARA
    USE gvect, ONLY : ngm, gcutm, gg
    USE mp_global, ONLY : nproc_pool, me_pool, intra_pool_comm
    USE parallel_include
    !
    implicit none
    !
    real(DP), parameter :: twothirds = 0.6666666666666666_dp
    real(DP), allocatable :: gcut_(:)
    integer :: np, ig, ngloc, ngpos, ierr
    !
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
#endif
  END SUBROUTINE sym_rho_init
  !
  !-----------------------------------------------------------------------
  subroutine sym_rho (nspin, rhog)
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
    USE gvect,                ONLY : ngm, g, ngl
    USE wavefunctions_module, ONLY : psic
#ifdef __PARA
    USE parallel_include
    USE mp_global,            ONLY : nproc_pool, me_pool, intra_pool_comm
#endif
!	use cell_base
    !
    implicit none
    !
    INTEGER, INTENT(IN) :: nspin
    complex(DP), INTENT(INOUT) :: rhog(ngm,nspin)
    !
    real(DP), allocatable :: g0(:,:), g_(:,:), gg_(:) 
    real(DP) :: gg0_, gg1_
    complex(DP), allocatable :: rhog_(:,:)
    integer :: is, ig, igl, np, ierr, ngm_, ngl_
    integer, allocatable :: igtongl_(:), index_(:)
    !
    !
#ifndef __PARA
    !
    !     The serial algorithm requires in input G-vectors in order of
    !     increasing module, plus a list of shells. All this is available
    !     in fixed-cell calculations, but in variable-cell calculations,
    !     the correct ordering is lost and variables g, gg, igtongl
    !     In the following we re-order the shells of G
    !     Variables ending with "_" refere to "re-ordered" G-vectors
    !
    ngm_=ngm
    ALLOCATE (rhog_(ngm_,nspin),g_(3,ngm_))
    rhog_(:,:) = rhog(:,:)
    g_(:,:) = g(:,:)
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
#endif
    !
    ! the local array of G-vectors is no longer ordered with increasing |G|
    ! sort G-vectors wrt module - index_ contains indices
    !
    ALLOCATE ( index_(ngm_) )
    ALLOCATE ( gg_(ngm_) )
    DO ig = 1, ngm_
       gg_(ig)= g_(1,ig)**2+g_(2,ig)**2+g_(3,ig)**2
    END DO
    ! 
    index_(1) = 0
    CALL hpsort_eps( ngm_, gg_, index_, eps8 )
    !
    DO igl = 1, 3
       gg_(:)= g_(igl,index_(:))
       g_(igl,:) = gg_(:)
    END DO
    DEALLOCATE ( gg_ )
    !
    ! use psic as work space for reordering
    !
    IF ( ngm_ > SIZE(psic) ) CALL errore ('symrho', &
         'Insufficient psic dimensions (should not happen)', ngm_)
    DO is = 1, nspin
       psic(1:ngm_)= rhog_(index_(:),is)
       rhog_(:,is)= psic(1:ngm_)
    END DO
    !
    ! find shells of G-vectors for this process
    !
    ALLOCATE ( igtongl_(ngm_) )
    ngl_ = 1
    igtongl_(1) = 1
    gg0_= g_(1,1)**2+g_(2,1)**2+g_(3,1)**2
    DO ig = 2, ngm_
       gg1_ = g_(1,ig)**2+g_(2,ig)**2+g_(3,ig)**2
       IF (gg1_ > gg0_ + eps6) THEN
          ngl_ = ngl_ + 1
          gg0_=gg1_
       END IF
       igtongl_(ig) = ngl_
    END DO
    !
    !   Now symmetrize
    !
    CALL sym_rho_serial ( ngm_, g_, ngl_, igtongl_, nspin, rhog_ )
    !
    DEALLOCATE ( igtongl_, g_ )
    !
    ! Restore starting ordering - use psic as work space
    !
#ifdef __PARA
    sendcnt(:) = sendcnt(:)/3
    recvcnt(:) = recvcnt(:)/3
    sdispls(:) = sdispls(:)/3
    rdispls(:) = rdispls(:)/3
    DO is = 1, nspin
       psic(1:ngm_) = (0.0_dp, 0.0_dp)
       psic(index_(:))= rhog_(:,is)
       rhog_(:,is)= psic(1:ngm_)
       !
       ! bring symmetrized rho(G) back to original distributed form
       !
       CALL mpi_alltoallv (rhog_(1,is), recvcnt, rdispls, MPI_DOUBLE_COMPLEX, &
            rhog (1,is), sendcnt, sdispls, MPI_DOUBLE_COMPLEX, &
            intra_pool_comm, ierr)
    END DO
#else
    DO is = 1, nspin
       psic(1:ngm_) = (0.0_dp, 0.0_dp)
       psic(index_(:))= rhog_(:,is)
       rhog(:,is)= psic(1:ngm_)
    END DO
#endif
    DEALLOCATE ( index_ )
    DEALLOCATE ( rhog_ )
    !
    RETURN
  END SUBROUTINE sym_rho
  !
  !-----------------------------------------------------------------------
  subroutine sym_rho_serial ( ngm_, g_, ngl_, igtongl_, nspin_, rhog_ )
    !-----------------------------------------------------------------------
    !
    !     symmetrize the charge density rho in reciprocal space 
    !     Serial algorithm - requires in input: 
    !     g_(3,ngm_)      list of G-vectors, ordered for increasing |G|
    !     ngl_            number of shells of G
    !     igtongl_(ngm_)  the ig-th G-vector belongs to igtongl_(ig) shell
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
    INTEGER, INTENT (IN) :: ngm_, ngl_, nspin_, igtongl_(ngm_)
    REAL(DP) , INTENT (IN) :: g_( 3, ngm_ )
    COMPLEX(DP) , INTENT (INOUT) :: rhog_( ngm_, nspin_ )
    !
    real(DP), allocatable :: g0(:,:)
    real(DP) :: sg(3), ft(3,48), arg
    complex(DP) :: fact, rhosum(2), mag(3), magrot(3), magsum(3)
    integer :: irot(48), ig, isg, igl, firstg, lastg, ng, ns, nspin_mag, is
    integer :: table(48, 48), invs(3, 3, 48)
    logical, allocatable :: done(:)
    logical :: non_symmorphic(48)
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
!!!!    nspin_mag = 1
       nspin_mag = 0 !!! TEMPTEMPTEMP 
       ! S^{-1} are needed as well
       call multable (nsym, s, table)
       call inverse_s (nsym, s, table, invs)
       !
    ELSE IF ( nspin_ == 1 .OR. nspin_ == 2 ) THEN
       nspin_mag = nspin_
    ELSE
       CALL errore('sym_rho_serial','incorrect value of nspin',nspin_)
    END IF
    !
    ! scan shells of G-vectors
    !
    firstg=1
    DO igl=1, ngl_
       !
       ! search for first and last G-vector in shell igl
       ! remember: G-vectors are in order of increasing module
       !
       lastg=firstg
       DO ig=firstg+1, ngm_
          IF ( igtongl_(ig) > igl ) EXIT
          lastg=ig
       END DO
       !
       ! symmetrize: \rho_sym(G) = \sum_S rho(SG) for all G-vectors in the star
       !
       ng = lastg-firstg+1
       allocate ( g0(3,ng), done(ng) )
       IF ( ng < 1 ) CALL errore('sym_rho_serial','internal error',1)
       !
       !  bring G-vectors to crystal axis
       !
       g0(:,1:ng) = g_(:,firstg:lastg)
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
                isg = firstg+irot(ns)-1
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
                   DO is=1,nspin_mag
                      rhosum(is) = rhosum(is) + rhog_(isg, is) * fact
                   END DO
                   IF ( nspin_ == 4 ) &
                        magsum(:) = magsum(:) + magrot(:) * fact
                ELSE
                   DO is=1,nspin_mag
                      rhosum(is) = rhosum(is) + rhog_(isg, is)
                   END DO
                   IF ( nspin_ == 4 ) &
                        magsum(:) = magsum(:) + magrot(:)
                END IF
             END DO
             !
             DO is=1,nspin_mag
                rhosum(is) = rhosum(is) / nsym
             END DO
             IF ( nspin_ == 4 ) magsum(:) = magsum(:) / nsym
             !
             DO ns=1,nsym
                isg = firstg+irot(ns)-1
                IF ( nspin_ == 4 ) THEN
                   ! rotate magnetization
                   magrot(:) = invs(1,:,ns) * magsum(1) + &
                               invs(2,:,ns) * magsum(2) + &
                               invs(3,:,ns) * magsum(3)
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
                   DO is=1,nspin_mag
                      rhog_(isg,is) = rhosum(is) * fact
                   END DO
                   IF ( nspin_ == 4 ) THEN
                      DO is=2,nspin_
                         rhog_(isg, is) = mag(is-1)*fact
                      END DO
                   END IF
                ELSE
                   DO is=1,nspin_mag
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
       firstg=lastg+1
    END DO
    ! check that all G-vectors were taken into account
    IF ( lastg /= ngm_.OR. firstg /= ngm_+1) &
         CALL errore('sym_rho_serial','internal error',3)
    !
    RETURN
  END SUBROUTINE sym_rho_serial

  SUBROUTINE sym_rho_deallocate ( )
    IF ( ALLOCATED (rdispls) ) DEALLOCATE (rdispls) 
    IF ( ALLOCATED (recvcnt) ) DEALLOCATE (recvcnt) 
    IF ( ALLOCATED (sdispls) ) DEALLOCATE (sdispls) 
    IF ( ALLOCATED (sendcnt) ) DEALLOCATE (sendcnt) 
  END SUBROUTINE sym_rho_deallocate
  !
END MODULE symme
