!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE sm1_psi( recalculate, ik, lda, n, m, psi, spsi)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !----------------------------------------------------------------------------
  !
  !    This routine applies the S^{-1} matrix to m wavefunctions psi
  !    and puts the results in spsi.
  !    Requires the products of psi with all beta functions
  !    in array becp(nkb,m) (calculated in h_psi or by ccalbec)
  ! input:
  !     recalculate decides if the overlap of beta functions is recalculated or not.
  !            this is needed e.g. if ions are moved and the overlap changes accordingly
  !     lda   leading dimension of arrays psi, spsi
  !     n     true dimension of psi, spsi
  !     m     number of states psi
  !     psi
  ! output:
  !     spsi  S^{-1}*psi
  !
  ! Modified by Osman Baris Malcioglu (2009)

  USE kinds,      ONLY : DP
  USE control_flags,      ONLY : gamma_only
  USE uspp,       ONLY : okvan, vkb, nkb, qq
  USE uspp_param, ONLY : nh, upf
  USE wvfct,      ONLY : igk, g2kin
  USE ions_base,  ONLY : ityp,nat,ntyp=>nsp
  USE mp,         ONLY : mp_sum
  USE mp_global,            ONLY : inter_pool_comm, intra_bgrp_comm
  USE lr_variables,   ONLY : lr_verbosity
  USE io_global,      ONLY : stdout
  !
  IMPLICIT NONE
  !
  ! ... First the dummy variables
  !
  LOGICAL, INTENT(in)           :: recalculate
  INTEGER, INTENT(in)           :: lda, n, m, ik
  COMPLEX(kind=DP), INTENT(in)  :: psi(lda,m)
  COMPLEX(kind=DP), INTENT(out) :: spsi(lda,m)
  !
  LOGICAL ::recalc
  !
  IF (lr_verbosity > 5) THEN
    WRITE(stdout,'("<lr_sm1_psi>")')
  ENDIF
  !
  CALL start_clock( 'lr_sm1_psi' )
  !
  recalc=recalculate
  IF ( gamma_only ) THEN
     CALL sm1_psi_gamma()
  ELSE
     !
     CALL sm1_psi_k()
     !
  ENDIF
  !
  CALL stop_clock( 'lr_sm1_psi' )
  !
  RETURN
  !
CONTAINS
  !-----------------------------------------------------------------------
  SUBROUTINE sm1_psi_gamma()
    !-----------------------------------------------------------------------
    !
    USE becmod,               ONLY : bec_type,becp,calbec
    !use real_beta,            only : ccalbecr_gamma,s_psir,fft_orbital_gamma
    USE realus,              ONLY : real_space, fft_orbital_gamma, initialisation_level, &
                                    bfft_orbital_gamma, calbec_rs_gamma, add_vuspsir_gamma, &
                                    v_loc_psir, s_psir_gamma, real_space_debug

    !
    IMPLICIT NONE
    !
    ! ... local variables
    !
    INTEGER :: ikb, jkb, ih, jh, na, nt, ijkb0, ibnd, ii
    ! counters
    real(kind=DP), ALLOCATABLE :: ps(:,:)
    real(kind=dp), ALLOCATABLE, SAVE :: BB_(:,:)
    LOGICAL, SAVE :: first_entry = .true.
    IF(first_entry) THEN
      IF(allocated(BB_)) DEALLOCATE(BB_)
      first_entry = .false.
      recalc=.true.
    ENDIF


    ! the product vkb and psi
    !
    ! ... initialize  spsi
    !
    CALL ZCOPY( lda * m, psi, 1, spsi, 1 )
    !
    ! ... The product with the beta functions
    !
    IF ( nkb == 0 .or. .not. okvan ) RETURN
    !
    !BB_ = sum
    !if (allocated(BB_)) then
    !  print *, "BB is allocated, ", BB_(1,1)
    !else
    !  print *, "BB is not allocated"
    !endif
    !OBM - For improved restart handling
    IF (.not.allocated(BB_)) recalc = .true.
    IF (recalc .and. allocated(BB_)) DEALLOCATE(BB_)

    IF(recalc) THEN
       ALLOCATE(BB_(nkb,nkb))
       BB_=0.d0
       CALL errore('sm1_psi','recalculating BB_ matrix',-1)
       !print *, "did you see the recalculating message?"
       IF (lr_verbosity > 1) THEN
          WRITE(stdout,'(5X,"Calculating S^-1")')
       ENDIF
       !call pw_gemm('Y',nkb,nkb,n,vkb,lda,vkb,lda,BB_,nkb)
       CALL calbec (n,vkb,vkb,BB_,nkb)
       ALLOCATE( ps( nkb, nkb ) )
       ps(:,:) = (0.d0)
       ijkb0 = 0
       DO nt=1,ntyp
          IF (upf(nt)%tvanp) THEN
             DO na=1,nat
                IF(ityp(na)==nt) THEN
                   DO ii=1,nkb
                      DO jh=1,nh(nt)
                         jkb=ijkb0 + jh
                         DO ih=1,nh(nt)
                            ikb = ijkb0 + ih
                            ps(ikb,ii) = ps(ikb,ii) + qq(ih,jh,nt)*BB_(jkb,ii)
                         ENDDO
                      ENDDO
                   ENDDO
                   ijkb0 = ijkb0+nh(nt)
                ENDIF
             ENDDO
          ELSE
             DO na = 1, nat
                IF ( ityp(na) == nt ) ijkb0 = ijkb0 + nh(nt)
             ENDDO
          ENDIF
       ENDDO

       DO ii=1,nkb
          ps(ii,ii) = ps(ii,ii) + 1.d0
       ENDDO

       CALL dinv_matrix(ps,nkb)
       BB_(:,:) = 0.d0
       ijkb0 = 0
       DO nt=1,ntyp
          IF (upf(nt)%tvanp) THEN
             DO na=1,nat
                IF(ityp(na)==nt) THEN
                   DO ii=1,nkb
                      DO jh=1,nh(nt)
                         jkb=ijkb0 + jh
                         DO ih=1,nh(nt)
                            ikb = ijkb0 + ih
                            BB_(ii,jkb) = BB_(ii,jkb) - ps(ii,ikb)*qq(ih,jh,nt)
                         ENDDO
                      ENDDO
                   ENDDO
                   ijkb0 = ijkb0+nh(nt)
                ENDIF
             ENDDO
          ELSE
             DO na = 1, nat
                IF ( ityp(na) == nt ) ijkb0 = ijkb0 + nh(nt)
             ENDDO
          ENDIF
       ENDDO

       DEALLOCATE(ps)
    ENDIF
    !print *, "BB is now, ", BB_(1,1)

    IF (real_space_debug>3) THEN !was 3
      DO ibnd=1,m,2
       CALL fft_orbital_gamma(psi,ibnd,m)
       CALL calbec_rs_gamma(ibnd,m,becp%r)
      ENDDO
    ELSE
     CALL calbec(n,vkb,psi,becp,m)
    !call pw_gemm('Y',nkb,m,n,vkb,lda,psi,lda,rbecp,nkb)
    ENDIF
    !
    ALLOCATE( ps( nkb, m ) )
!    ps(:,:) = 0.D0
    !
!    do ibnd=1,m
!       do jkb=1,nkb
!          do ii=1,nkb
!             ps(jkb,ibnd) = ps(jkb,ibnd)+BB_(jkb,ii)*rbecp(ii,ibnd)
!          enddo
!       enddo
!    enddo
    !
    CALL DGEMM( 'N','N',nkb,m,nkb,1.d0,BB_,nkb,becp%r,nkb,0.d0,ps,nkb)


!   do ibnd=1,m
!      do ii=1,nkb
!          call ZAXPY(n,cmplx(ps(ii,ibnd),0.0d0,dp),vkb(1,ii),1,spsi(1,ibnd),1)
!       enddo
!    enddo
    CALL DGEMM('N','N',2*n,m,nkb,1.d0,vkb,2*lda,ps,nkb,1.d0,spsi,2*lda)

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
    USE becmod,        ONLY : bec_type,becp,calbec
    !USE lr_variables,    ONLY: igk_k, npw_k
    USE realus,        ONLY : igk_k,npw_k
    USE klist,         ONLY : nks, xk
    !
    IMPLICIT NONE
    !
    ! ... local variables
    !
    INTEGER :: ikb, jkb, ih, jh, na, nt, ijkb0, ibnd, ii, ik1
    ! counters
    COMPLEX(kind=DP), ALLOCATABLE :: ps(:,:)
    COMPLEX(kind=dp), ALLOCATABLE, SAVE :: BB_(:,:,:)

    ! the product vkb and psi
    !
    ! ... initialize  spsi
    !
    CALL ZCOPY( lda * m, psi, 1, spsi, 1 )
    !
    ! ... The product with the beta functions
    !
    IF ( nkb == 0 .or. .not. okvan ) RETURN
    !
    IF (.not.allocated(BB_)) recalc = .true.
    IF (recalc .and. allocated(BB_)) DEALLOCATE(BB_)

    IF(recalc) THEN
       ALLOCATE(BB_(nkb,nkb,nks))
       BB_=(0.d0,0.d0)
       CALL errore('sm1_psi','recalculating BB_ matrix',-1)

       ALLOCATE( ps( nkb, nkb ) )

       DO ik1 = 1,nks
          CALL init_us_2(npw_k(ik1),igk_k(:,ik1),xk(1,ik1),vkb)
          CALL zgemm('C','N',nkb,nkb,npw_k(ik1),(1.d0,0.d0),vkb,lda,vkb,lda,(0.d0,0.d0),BB_(1,1,ik1),nkb)
#ifdef __MPI
          !CALL reduce( 2 * nkb * nkb, BB_(:,:,ik1) )
          CALL mp_sum(BB_(:,:,ik1), intra_bgrp_comm)
#endif

          ps(:,:) = (0.d0,0.d0)
          ijkb0 = 0
          DO nt=1,ntyp
             IF (upf(nt)%tvanp) THEN
                DO na=1,nat
                   IF(ityp(na)==nt) THEN
                      DO ii=1,nkb
                         DO jh=1,nh(nt)
                            jkb=ijkb0 + jh
                            DO ih=1,nh(nt)
                               ikb = ijkb0 + ih
                               ps(ikb,ii) = ps(ikb,ii) + BB_(jkb,ii, ik1)*qq(ih,jh,nt)
                            ENDDO
                         ENDDO
                      ENDDO
                      ijkb0 = ijkb0+nh(nt)
                   ENDIF
                ENDDO
             ELSE
                DO na = 1, nat
                   IF ( ityp(na) == nt ) ijkb0 = ijkb0 + nh(nt)
                ENDDO
             ENDIF
          ENDDO

          DO ii=1,nkb
             ps(ii,ii) = ps(ii,ii) + (1.d0,0.d0)
          ENDDO

          CALL zinv_matrix(ps,nkb)
          BB_(:,:,ik1) = (0.d0,0.d0)
          ijkb0 = 0
          DO nt=1,ntyp
             IF (upf(nt)%tvanp) THEN
                DO na=1,nat
                   IF(ityp(na)==nt) THEN
                      DO ii=1,nkb
                         DO jh=1,nh(nt)
                            jkb=ijkb0 + jh
                            DO ih=1,nh(nt)
                               ikb = ijkb0 + ih
                               BB_(ii,jkb,ik1) = BB_(ii,jkb,ik1) - ps(ii,ikb)*qq(ih,jh,nt)
                            ENDDO
                         ENDDO
                      ENDDO
                      ijkb0 = ijkb0+nh(nt)
                   ENDIF
                ENDDO
             ELSE
                DO na = 1, nat
                   IF ( ityp(na) == nt ) ijkb0 = ijkb0 + nh(nt)
                ENDDO
             ENDIF
          ENDDO
       ENDDO
       DEALLOCATE(ps)
    ENDIF

    CALL init_us_2(npw_k(ik),igk_k(:,ik),xk(1,ik),vkb)
    !call ccalbec( nkb, lda, n, m, becp, vkb, psi )
    CALL calbec(n,vkb,psi,becp,m)

    !
    ALLOCATE( ps( nkb, m ) )
    ps(:,:) = (0.d0,0.d0)
    !
    DO ibnd=1,m
       DO jkb=1,nkb
          DO ii=1,nkb
             ps(jkb,ibnd) = ps(jkb,ibnd)+BB_(jkb,ii,ik)*becp%k(ii,ibnd)
          ENDDO
       ENDDO
    ENDDO
    !
    !
    CALL ZGEMM( 'N', 'N', n, m, nkb, (1.D0, 0.D0), vkb, &
         lda, ps, nkb, (1.D0, 0.D0), spsi, lda )


    DEALLOCATE( ps )
    !
    RETURN
    !
  END SUBROUTINE sm1_psi_k
  !
END SUBROUTINE sm1_psi
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE dinv_matrix(M,N)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  USE kinds,      ONLY : DP

  IMPLICIT NONE

  INTEGER :: N                              !  matrix dimension
  real(kind=dp), DIMENSION(0:N-1,0:N-1) :: M ! MAtrix to be inverted

  real(kind=dp), DIMENSION(:), ALLOCATABLE :: work
  INTEGER, DIMENSION(:), ALLOCATABLE :: ipiv

  INTEGER :: i,lwork,info
  INTEGER, SAVE :: lworkfact

  data lworkfact /64/

  lwork = lworkfact*N

  ALLOCATE(ipiv(0:N-1))
  ALLOCATE(work(1:lwork))

! Factorize Matrix M

  CALL dgetrf( N, N, M, N, ipiv, info )
  IF (info/=0) THEN
     CALL errore('dinv_matrix','error in dgetrf',info)
  ENDIF

! Invert Matrix

  CALL dgetri( N, M, N, ipiv, work, lwork, info )
  IF (info/=0) THEN
     CALL errore('dinv_matrix','error in dgetri',info)
  ELSE
     lworkfact = int(work(1)/N)
  ENDIF

  DEALLOCATE(work)
  DEALLOCATE(ipiv)

END SUBROUTINE dinv_matrix

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE zinv_matrix(M,N)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  USE kinds,      ONLY : DP

  IMPLICIT NONE

  INTEGER :: N                            !  matrix dimension
  COMPLEX(kind=dp), DIMENSION(0:N-1,0:N-1) :: M ! MAtrix to be inverted

  COMPLEX(kind=dp), DIMENSION(:), ALLOCATABLE :: work
  INTEGER, DIMENSION(:), ALLOCATABLE :: ipiv

  INTEGER :: i,lwork,info
  INTEGER, SAVE :: lworkfact

  data lworkfact /64/

  lwork = lworkfact*N

  ALLOCATE(ipiv(0:N-1))
  ALLOCATE(work(1:lwork))

! Factorize Matrix M

  CALL zgetrf( N, N, M, N, ipiv, info )
  IF (info/=0) THEN
     CALL errore('zinv_matrix','error in zgetrf',info)
  ENDIF

! Invert Matrix

  CALL zgetri( N, M, N, ipiv, work, lwork, info )
  IF (info/=0) THEN
     CALL errore('zinv_matrix','error in zgetri',info)
  ELSE
     lworkfact = int(work(1)/N)
  ENDIF

  DEALLOCATE(work)
  DEALLOCATE(ipiv)

END SUBROUTINE zinv_matrix
