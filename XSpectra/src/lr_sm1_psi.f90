!
! Copyright (C) 2001-2006 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE sm1_psi( recalc, lda, n, m, psi, spsi)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !----------------------------------------------------------------------------
  !
  !    This routine applies the S^{-1} matrix to m wavefunctions psi
  !    and puts the results in spsi.
  !    Requires the products of psi with all beta functions
  !    in array becp(nkb,m) (calculated in h_psi or by ccalbec)
  ! input:
  !     recalc decides if the overlap of beta functions is recalculated or not.
  !            this is needed e.g. if ions are moved and the overlap changes accordingly
  !     lda   leading dimension of arrays psi, spsi
  !     n     true dimension of psi, spsi
  !     m     number of states psi
  !     psi
  ! output:
  !     spsi  S^{-1}*psi
  !
  !   Original routine written by Ralph Gebauer
  !        Modified by Christos Gougoussis
  !
  USE kinds,      ONLY : DP
  USE control_flags,      ONLY : gamma_only
  USE uspp,         ONLY : okvan, vkb, nkb, qq
  USE uspp_param, ONLY : upf, nh
  USE ldaU,       ONLY : lda_plus_u
  USE ions_base,   ONLY : nat, ntyp => nsp, ityp  
  use becmod, only : calbec
  !
  IMPLICIT NONE
  !
  ! ... First the dummy variables
  !
  LOGICAL          :: recalc
  INTEGER          :: lda, n, m
  COMPLEX(KIND=DP) :: psi(lda,m), spsi(lda,m)
  !
  CALL start_clock( 'sm1' )  
  !
  IF ( gamma_only ) THEN
     CALL sm1_psi_gamma()
  ELSE
     !
     CALL sm1_psi_k()
     !
  END IF
  !
  CALL stop_clock( 'sm1' )
  !
  RETURN
  !
CONTAINS
  !-----------------------------------------------------------------------
  SUBROUTINE sm1_psi_gamma()
    !-----------------------------------------------------------------------
    !
    USE becmod,  ONLY : becp
    !
    IMPLICIT NONE
    !
    ! ... local variables
    !
    INTEGER :: ikb, jkb, ih, jh, na, nt, ijkb0, ibnd, ii
    ! counters
    real(KIND=DP), ALLOCATABLE :: ps(:,:)
    real(kind=dp), allocatable, save :: BB_(:,:)

    ! the product vkb and psi
    !
    ! ... initialize  spsi
    !
    CALL zcopy( lda * m, psi, 1, spsi, 1 )
    !
    ! ... The product with the beta functions
    !
    IF ( nkb == 0 .OR. .NOT. okvan ) RETURN
    !
    if(.not.allocated(BB_)) then
       allocate(BB_(nkb,nkb))
       recalc = .true.
    endif

    if(recalc) then
       call errore('sm1_psi','recalculating BB_ matrix',-1)

       call pw_gemm('Y',nkb,nkb,n,vkb,lda,vkb,lda,BB_,nkb) 

       ALLOCATE( ps( nkb, nkb ) )    
       ps(:,:) = (0.d0)
       ijkb0 = 0
       do nt=1,ntyp
          if (upf(nt)%tvanp) then
             do na=1,nat
                if(ityp(na).eq.nt) then
                   do ii=1,nkb
                      do jh=1,nh(nt)
                         jkb=ijkb0 + jh
                         do ih=1,nh(nt)
                            ikb = ijkb0 + ih
                            ps(ikb,ii) = ps(ikb,ii) + qq(ih,jh,nt)*BB_(jkb,ii)
                         enddo
                      enddo
                   enddo
                   ijkb0 = ijkb0+nh(nt)
                endif
             enddo
          else
             DO na = 1, nat
                IF ( ityp(na) == nt ) ijkb0 = ijkb0 + nh(nt)
             END DO
          endif
       enddo

       do ii=1,nkb
          ps(ii,ii) = ps(ii,ii) + 1.d0
       enddo

       call dinv_matrix(ps,nkb)
       BB_(:,:) = 0.d0
       ijkb0 = 0
       do nt=1,ntyp
          if (upf(nt)%tvanp) then
             do na=1,nat
                if(ityp(na).eq.nt) then
                   do ii=1,nkb
                      do jh=1,nh(nt)
                         jkb=ijkb0 + jh
                         do ih=1,nh(nt)
                            ikb = ijkb0 + ih
                            BB_(ii,jkb) = BB_(ii,jkb) - ps(ii,ikb)*qq(ih,jh,nt)
                         enddo
                      enddo
                   enddo
                   ijkb0 = ijkb0+nh(nt)
                endif
             enddo
          else
             DO na = 1, nat
                IF ( ityp(na) == nt ) ijkb0 = ijkb0 + nh(nt)
             END DO
          endif
       enddo

       deallocate(ps)
    endif

    call pw_gemm('Y',nkb,m,n,vkb,lda,psi,lda,becp%r,nkb) 

    !
    ALLOCATE( ps( nkb, m ) )    
    ps(:,:) = 0.D0
    !
    do ibnd=1,m
       do jkb=1,nkb
          do ii=1,nkb
             ps(jkb,ibnd) = ps(jkb,ibnd)+BB_(jkb,ii)*becp%r(ii,ibnd)
          enddo
       enddo
    enddo
    !

    do ibnd=1,m
       do ii=1,nkb
          call zaxpy(n,CMPLX(ps(ii,ibnd),0.0d0,dp),vkb(1,ii),1,spsi(1,ibnd),1) 
       enddo
    enddo
    !
    DEALLOCATE( ps )
    !
    RETURN
    !
  END SUBROUTINE sm1_psi_gamma
  !
  !-----------------------------------------------------------------------
  SUBROUTINE sm1_psi_k()
    !-----------------------------------------------------------------------
    !
    ! ... k-points version
    !
    USE becmod,  ONLY : becp
    USE klist, only :         xk
    USE mp,         only :  mp_sum  ! CG
    USE mp_global,  ONLY : intra_pool_comm ! CG
    !
    IMPLICIT NONE
    !
    ! ... local variables
    !
    INTEGER :: ikb, jkb, ih, jh, na, nt, ijkb0, ibnd, ii, ik1
    ! counters
    complex(KIND=DP), ALLOCATABLE :: ps(:,:)
    complex(kind=dp), allocatable, save :: BB_(:,:)

    ! the product vkb and psi
    !
    ! ... initialize  spsi
    !
    CALL zcopy( lda * m, psi, 1, spsi, 1 )
    !
    ! ... The product with the beta functions
    !
    IF ( nkb == 0 .OR. .NOT. okvan ) RETURN
    !
    if(.not.allocated(BB_)) then
       allocate(BB_(nkb,nkb))
       recalc = .true.
    endif

    if(recalc) then
       call errore('sm1_psi','recalculating BB_ matrix',-1)

       ALLOCATE( ps( nkb, nkb ) )    

       call zgemm('C','N',nkb,nkb,n,(1.d0,0.d0),vkb,lda,vkb,lda,(0.d0,0.d0),BB_(1,1),nkb)
       !******          CALL reduce( 2 * nkb * nkb, BB_ )
       CALL mp_sum(  BB_, intra_pool_comm ) !CG

       ps(:,:) = (0.d0,0.d0)
       ijkb0 = 0
       do nt=1,ntyp
          if (upf(nt)%tvanp) then
             do na=1,nat
                if(ityp(na).eq.nt) then
                   do ii=1,nkb
                      do jh=1,nh(nt)
                         jkb=ijkb0 + jh
                         do ih=1,nh(nt)
                            ikb = ijkb0 + ih
                            ps(ikb,ii) = ps(ikb,ii) + BB_(jkb,ii)*qq(ih,jh,nt)
                         enddo
                      enddo
                   enddo
                   ijkb0 = ijkb0+nh(nt)
                endif
             enddo
          else
             DO na = 1, nat
                IF ( ityp(na) == nt ) ijkb0 = ijkb0 + nh(nt)
             END DO
          endif
       enddo

       do ii=1,nkb
          ps(ii,ii) = ps(ii,ii) + (1.d0,0.d0)
       enddo

       call zinv_matrix(ps,nkb)
       BB_(:,:) = (0.d0,0.d0)
       ijkb0 = 0
       do nt=1,ntyp
          if (upf(nt)%tvanp) then
             do na=1,nat
                if(ityp(na).eq.nt) then
                   do ii=1,nkb
                      do jh=1,nh(nt)
                         jkb=ijkb0 + jh
                         do ih=1,nh(nt)
                            ikb = ijkb0 + ih
                            ! BB_(ii,jkb) = BB_(ii,jkb) - ps(ii,jkb)*qq(ih,jh,nt) ! this seems false
                            BB_(ii,jkb) = BB_(ii,jkb) - ps(ii,ikb)*qq(ih,jh,nt) ! modified by CG
                         enddo
                      enddo
                   enddo
                   ijkb0 = ijkb0+nh(nt)
                endif
             enddo
          else
             DO na = 1, nat
                IF ( ityp(na) == nt ) ijkb0 = ijkb0 + nh(nt)
             END DO
          endif
       enddo
       deallocate(ps)
    endif

    !     call calbec ( lda, vkb, psi, becp ) ! erreur ici ?

    call calbec ( n, vkb, psi, becp%k )

    !
    ALLOCATE( ps( nkb, m ) )    
    ps(:,:) = (0.d0,0.d0)
    !
    do ibnd=1,m
       do jkb=1,nkb
          do ii=1,nkb
             ps(jkb,ibnd) = ps(jkb,ibnd)+BB_(jkb,ii)*becp%k(ii,ibnd)
          enddo
       enddo
    enddo
    !
    !
    CALL zgemm( 'N', 'N', n, m, nkb, (1.D0, 0.D0), vkb, &
         lda, ps, nkb, (1.D0, 0.D0), spsi, lda )


    !    CALL zcopy( lda * m, psi, 1, spsi, 1 ) ! remove this !!!
    DEALLOCATE( ps )
    !
    RETURN
    !
  END SUBROUTINE sm1_psi_k
  !
END subroutine sm1_psi
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine dinv_matrix(M,N)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  USE kinds,      ONLY : DP

  implicit none

  integer :: N                              !  matrix dimension
  real(kind=dp), dimension(0:N-1,0:N-1) :: M ! MAtrix to be inverted

  real(kind=dp), dimension(:), allocatable :: work
  integer, dimension(:), allocatable :: ipiv

  integer :: i,lwork,info
  integer, save :: lworkfact

  data lworkfact /64/

  lwork = lworkfact*N

  allocate(ipiv(0:N-1))
  allocate(work(1:lwork))

  ! Factorize Matrix M

  call dgetrf( N, N, M, N, ipiv, info )
  if (info.ne.0) then
     call errore('dinv_matrix','error in dgetrf',info)
  endif

  ! Invert Matrix

  call dgetri( N, M, N, ipiv, work, lwork, info )
  if (info.ne.0) then
     call errore('dinv_matrix','error in dgetri',info)
  else
     lworkfact = int(work(1)/N)
  endif

  deallocate(work)
  deallocate(ipiv)

end subroutine dinv_matrix

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine zinv_matrix(M,N)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  USE kinds,      ONLY : DP

  implicit none

  integer :: N                            !  matrix dimension
  complex(kind=dp), dimension(0:N-1,0:N-1) :: M ! MAtrix to be inverted

  complex(kind=dp), dimension(:), allocatable :: work
  integer, dimension(:), allocatable :: ipiv

  integer :: i,lwork,info
  integer, save :: lworkfact

  data lworkfact /64/

  lwork = lworkfact*N

  allocate(ipiv(0:N-1))
  allocate(work(1:lwork))

  ! Factorize Matrix M

  call zgetrf( N, N, M, N, ipiv, info )
  if (info.ne.0) then
     call errore('zinv_matrix','error in zgetrf',info)
  endif

  ! Invert Matrix

  call zgetri( N, M, N, ipiv, work, lwork, info )
  if (info.ne.0) then
     call errore('zinv_matrix','error in zgetri',info)
  else
     lworkfact = int(work(1)/N)
  endif

  deallocate(work)
  deallocate(ipiv)

end subroutine zinv_matrix

!
! Copyright (C) 2002-2005 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------------
SUBROUTINE pw_gemm( sum_over_nodes, na, nb, n, a, lda, b, ldb, c, ldc )
  !----------------------------------------------------------------------------
  !
  ! ... matrix times matrix with summation index running on G-vectors or PWs
  ! ... c(ij)=real(a(ik)*b(kj)) using half G vectors or half PWs
  !
  ! ccalbec( nkb, npwx, npw, nbnd, vkb, psi, bec ) =>
  !    pw_gemm( 'Y', nkb, nbnd, npw, vkb, npwx, psi, npwx, bec, nkb )
  !
  !
  USE kinds, ONLY : DP
  USE gvect, ONLY : gstart
  USE mp,         only :  mp_sum  ! CG
  USE mp_global,  ONLY : intra_pool_comm ! CG
  !
  IMPLICIT NONE
  !
  ! ... input
  !
  CHARACTER(LEN=1) :: sum_over_nodes
  INTEGER          :: na, nb, n, lda, ldb, ldc
  COMPLEX(DP)      :: a(lda,na), b(ldb,nb)
  !
  ! ... output
  !
  REAL(DP) :: c(ldc,nb)
  !
  !
  IF ( na == 0 .OR. nb == 0 ) RETURN
  !
  CALL start_clock( 'pw_gemm' )
  !
  IF ( nb == 1 ) THEN
     !
     CALL dgemv( 'C', 2*n, na, 2.D0, a, 2*lda, b, 1, 0.D0, c, 1 )
     !
     IF ( gstart == 2 ) c(:,1) = c(:,1) - a(1,:) * b(1,1)
     !
  ELSE
     !
     CALL dgemm( 'C', 'N', na, nb, 2*n, 2.D0, a, 2*lda, b, 2*ldb, 0.D0, c, ldc )
     !
     IF ( gstart == 2 ) &
          CALL DGER( na, nb, -1.D0, a, 2*lda, b, 2*ldb, c, ldc )
     !
  END IF
  !
  !********  IF ( sum_over_nodes == 'y' .OR. &
  !********       sum_over_nodes == 'Y' ) CALL reduce( ldc*nb, c )
  IF ( sum_over_nodes == 'y' .OR. &
       sum_over_nodes == 'Y' )  CALL mp_sum(  c, intra_pool_comm ) !CG

  !
  CALL stop_clock( 'pw_gemm' )
  !
  RETURN
  !
END SUBROUTINE pw_gemm
