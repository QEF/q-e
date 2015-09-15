!
! Copyright (C) 2013-2015 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! Written by Lorenzo Paulatto (2012-2013) 
! Gamma-only tricks by Simon Binnie
! G-space code based on addusdens.f90 and compute_becsum.f90,
! modified by Paolo Giannozzi (2015) for speed
! Real space code based on realus.f90
!-----------------------------------------------------------------------
MODULE us_exx
  !-----------------------------------------------------------------------
  ! Most of the USPP+EXX code is here.
  ! Notes: 
  !   * compute_becxx is still in exx.f90 as it uses plenty of global variables from there
  !   * some tests and loops are done directly in exx.f90        
  !   * PAW specific parts are in paw_exx.f90
  !
  USE kinds,   ONLY : DP
  USE becmod,  ONLY : bec_type, calbec, ALLOCATE_bec_type, DEALLOCATE_bec_type
  !
  IMPLICIT NONE
  SAVE
  !
  TYPE(bec_type),ALLOCATABLE :: becxx(:) ! <beta_I|phi_j,k>, with the wavefunctions from exxbuff
  ! the visible index is k; while I and J are inside bec_type

  COMPLEX(DP),ALLOCATABLE :: becxx_gamma(:,:)  ! gamma only version of becxx%r
                                               ! two bands stored per stripe 
  COMPLEX(DP),ALLOCATABLE :: qgm(:,:)          ! used in addusxx_g and newdxx_g

 CONTAINS ! ~~+~~---//--~~~-+
  !
  FUNCTION bexg_merge( w, m,n, imin, imax, i)
    ! used at Gamma point when number of bands is odd,
    ! especially for band parallelisation when band group is odd
    ! returns w(i)+i w(i+1) if imin<=i<imax
    !         w(i)          if i==imax
    !         iw(i+1)       if i==imin-1
    !         0             otherwise
    COMPLEX(DP) :: bexg_merge(m)
    !
    INTEGER,INTENT(in)  :: m,n, imin, imax, i
    REAL(DP),INTENT(in) :: w(m,n)
    !
    bexg_merge = (0._dp, 0._dp)
    IF(imin<=i .and. i<imax) THEN
      bexg_merge = CMPLX(w(:,i),w(:,i+1), kind=DP)
    ELSE IF(i==imax) THEN
      bexg_merge = CMPLX(w(:,i),   0._dp, kind=DP)
    ELSE IF(i+1==imin) THEN
      bexg_merge = CMPLX( 0._dp,w(:,i+1), kind=DP)
    ELSE
      bexg_merge = CMPLX( 0._dp,   0._dp, kind=DP)
    ENDIF
    RETURN
  END FUNCTION bexg_merge

  !-----------------------------------------------------------------------
  SUBROUTINE qvan_init ( xkq, xk )
  !-----------------------------------------------------------------------
    ! allocate and store augmentation charges in G space Q(G) for USPP
    !
    USE ions_base,           ONLY : ntyp => nsp
    USE uspp_param,          ONLY : upf, nh, lmaxq
    USE gvect,               ONLY : g
    USE gvecs,               ONLY : ngms
    IMPLICIT NONE
    !
    REAL(dp),   INTENT(in)    :: xkq(3), xk(3)
    !
    REAL(dp), ALLOCATABLE :: ylmk0(:,:), qmod(:), q(:,:), qq(:)
    INTEGER  :: nij, ijh, ig, nt, ih, jh
    !
    CALL start_clock('qvan_init')
    !
    ! nij = number of (ih,jh) pairs for all atom types
    !
    nij = 0
    DO nt = 1, ntyp
       IF ( upf(nt)%tvanp ) nij = nij + nh(nt)*(nh(nt)+1)/2
    END DO
    ALLOCATE ( qgm(ngms,nij) )
    !
    ALLOCATE( ylmk0(ngms, lmaxq * lmaxq), qmod(ngms) )
    ALLOCATE( q(3,ngms), qq(ngms) )
    DO ig = 1, ngms
       q(:,ig) = xk(:) - xkq(:) + g(:,ig)
       qq(ig)  = SUM(q(:,ig)**2)
       qmod(ig)= SQRT(qq(ig))
    ENDDO
    CALL ylmr2 (lmaxq * lmaxq, ngms, q, qq, ylmk0)
    DEALLOCATE(qq, q)
    !
    ! ijh = position of (ih,jh) pairs for atom type nt
    !
    ijh = 0
    DO nt = 1, ntyp
       IF ( upf(nt)%tvanp ) THEN
          DO ih = 1, nh(nt)
             DO jh = ih, nh(nt)
                ijh = ijh + 1
                CALL qvan2(ngms, ih, jh, nt, qmod, qgm(1,ijh), ylmk0)
             END DO
          END DO
       END IF
    END DO
    DEALLOCATE(qmod, ylmk0)
    CALL stop_clock('qvan_init')

  END SUBROUTINE qvan_init

  !-----------------------------------------------------------------------
  SUBROUTINE qvan_clean ( )
  !-----------------------------------------------------------------------
    !
    DEALLOCATE (qgm)
    !
  END SUBROUTINE qvan_clean
  !
  !-----------------------------------------------------------------------
  SUBROUTINE addusxx_g(rhoc, xkq, xk, flag, becphi_c, becpsi_c, becphi_r, becpsi_r )
    !-----------------------------------------------------------------------
    ! 
    ! Add US contribution to rhoc for hybrid functionals
    !   flag = 'c': add complex contribution
    !   flag = 'r': add real contribution to the real part of rhoc 
    !   flag = 'i': add real contribution to the imaginary part
    ! The two latter cases are used together with gamma tricks to store contributions
    ! from two bands into the real and the imaginary part separately
    !
    USE constants,           ONLY : tpi
    USE ions_base,           ONLY : nat, ntyp => nsp, ityp, tau
    USE uspp,                ONLY : nkb, vkb,  okvan, indv_ijkb0, ijtoh
    USE uspp_param,          ONLY : upf, nh, nhm, lmaxq
    USE fft_base,            ONLY : dffts
    USE gvect,               ONLY : ngm, nl, nlm, g, &
                                    eigts1, eigts2, eigts3, mill, gstart
    USE gvecs,               ONLY : ngms, nls, nlsm
    USE cell_base,           ONLY : tpiba
    USE control_flags,       ONLY : gamma_only
    IMPLICIT NONE
    !
    ! In input I get a slice of <beta|left> and <beta|right> only for this kpoint and this band
    COMPLEX(DP),INTENT(inout) :: rhoc(dffts%nnr)
    COMPLEX(DP),INTENT(in), OPTIONAL  :: becphi_c(nkb), becpsi_c(nkb)
    REAL(DP),   INTENT(in), OPTIONAL  :: becphi_r(nkb), becpsi_r(nkb)
    REAL(DP),   INTENT(in)    :: xkq(3), xk(3)
    CHARACTER(LEN=1), INTENT(in) :: flag
    !
    ! ... local variables
    !
    COMPLEX(DP),ALLOCATABLE :: aux1(:), aux2(:), eigqts(:)
    INTEGER :: ikb, jkb, ijkb0, ih, jh, na, nt, ig, nij, ijh
    COMPLEX(DP) :: becfac_c
    REAL(DP) :: arg, becfac_r
    LOGICAL :: add_complex, add_real, add_imaginary
    !
    IF(.not.okvan) RETURN
    CALL start_clock( 'addusxx' )
    !
    add_complex = ( flag=='c' .OR. flag=='C' )
    add_real    = ( flag=='r' .OR. flag=='R' )
    add_imaginary=( flag=='i' .OR. flag=='I' )
    IF ( .NOT. (add_complex .OR. add_real .OR. add_imaginary) ) &
       CALL errore('addusxx_g', 'called with incorrect flag: '//flag, 1 )
    IF ( .NOT. gamma_only .AND. ( add_real .OR. add_imaginary) ) &
       CALL errore('addusxx_g', 'need gamma tricks for this flag: '//flag, 2 )
    IF ( gamma_only .AND. add_complex ) &
       CALL errore('addusxx_g', 'gamma trick not good for this flag: '//flag, 3 )
    IF ( ( add_complex .AND. (.NOT. PRESENT(becphi_c) .OR. .NOT. PRESENT(becpsi_c) ) ) .OR. &
         ( add_real    .AND. (.NOT. PRESENT(becphi_r) .OR. .NOT. PRESENT(becpsi_r) ) ) .OR. &
         ( add_imaginary.AND.(.NOT. PRESENT(becphi_r) .OR. .NOT. PRESENT(becpsi_r) ) ) )    &
       CALL errore('addusxx_g', 'called with incorrect arguments', 2 )
    !
    ALLOCATE( aux1(ngms), aux2(ngms) )
    ALLOCATE(eigqts(nat))
    !
    DO na = 1, nat
      arg = tpi* SUM( (xk(:) - xkq(:))*tau(:,na) )
      eigqts(na) = CMPLX( COS(arg), -SIN(arg), kind=DP)
    END DO
    !
    nij = 0
    DO nt = 1, ntyp
       !
       IF ( upf(nt)%tvanp ) THEN
          !
          DO na = 1, nat
             !
             IF (ityp(na)==nt) THEN
                !
                ! ijkb0 points to the manifold of beta functions for atom na
                !
                ijkb0 = indv_ijkb0(na) 
                !
                aux2(:) = (0.0_dp, 0.0_dp)
                DO ih = 1, nh(nt)
                   ikb = ijkb0 + ih
                   aux1(:) = (0.0_dp, 0.0_dp)
                   DO jh = 1, nh(nt)
                      jkb = ijkb0 + jh
                      IF ( add_complex ) THEN
!$omp parallel do default(shared) private(ig)
                         DO ig = 1, ngms
                            aux1(ig) = aux1(ig) + qgm(ig,nij+ijtoh(ih,jh,nt)) * &
                                 becpsi_c(jkb)
                         ENDDO
!$omp end parallel do
                      ELSE
!$omp parallel do default(shared) private(ig)
                         DO ig = 1, ngms
                            aux1(ig) = aux1(ig) + qgm(ig,nij+ijtoh(ih,jh,nt)) * &
                                 becpsi_r(jkb)
                         ENDDO
!$omp end parallel do
                      END IF
                   END DO
                   IF ( add_complex ) THEN
!$omp parallel do default(shared) private(ig)
                      DO ig = 1,ngms
                         aux2(ig) = aux2(ig) + aux1(ig) * CONJG(becphi_c(ikb))
                      ENDDO
!$omp end parallel do
                   ELSE
!$omp parallel do default(shared) private(ig)
                      DO ig = 1,ngms
                         aux2(ig) = aux2(ig) + aux1(ig) * becphi_r(ikb)
                      ENDDO
!$omp end parallel do
                   END IF
                END DO
!$omp parallel do default(shared) private(ig)
                DO ig = 1, ngms
                   aux2(ig) = aux2(ig) * eigqts(na) * &
                                 eigts1 (mill (1,ig), na) * &
                                 eigts2 (mill (2,ig), na) * &
                                 eigts3 (mill (3,ig), na)
                ENDDO
!$omp end parallel do
                IF ( add_complex ) THEN
                   DO ig = 1, ngms
                      rhoc(nls(ig)) = rhoc(nls(ig)) + aux2(ig)
                   END DO
                ELSE IF ( add_real ) THEN
                   DO ig = 1,ngms
                      rhoc(nls(ig)) = rhoc(nls(ig)) + aux2(ig)
                   ENDDO
                   DO ig = gstart,ngms
                      rhoc(nlsm(ig)) = rhoc(nlsm(ig)) + CONJG(aux2(ig))
                   ENDDO
                ELSE IF ( add_imaginary ) THEN
                   DO ig = 1,ngms
                      rhoc(nls(ig)) = rhoc(nls(ig)) + (0.0_dp,1.0_dp) * aux2(ig)
                   ENDDO
                   DO ig = gstart,ngms
                      rhoc(nlsm(ig)) = rhoc(nlsm(ig)) + (0.0_dp,1.0_dp)* &
                           CONJG( aux2(ig) )
                   ENDDO
                ENDIF
             ENDIF
          ENDDO   ! nat
          !
          nij = nij + nh(nt)*(nh(nt)+1)/2
          !
       END IF
       !
    ENDDO   ! nt
    !
    DEALLOCATE( eigqts, aux2, aux1)
    !
    CALL stop_clock( 'addusxx' )
    !
    RETURN
  !
  !-----------------------------------------------------------------------
  END SUBROUTINE addusxx_g
  !-----------------------------------------------------------------------
  !
  !-----------------------------------------------------------------------
  SUBROUTINE newdxx_g(vc, xkq, xk, flag, deexx, becphi_r, becphi_c)
    !-----------------------------------------------------------------------
    !
    ! This subroutine computes some sort of EXX contribution to the non-local 
    ! part of the hamiltonian. 
    !   alpha_Ii = \int \sum_Jj Q_IJ(r) V^{i,j}_Fock <beta_J|phi_j> d^3(r)
    ! The actual contribution will be (summed outside)
    !  H = H+\sum_I |beta_I> alpha_Ii
    !   flag = 'c': V(G) is contained in complex array vc
    !   flag = 'r': V(G)=v_1(G)+i v_2(G): select v_1(G)
    !   flag = 'i': V(G)=v_1(G)+i v_2(G): select v_2(G)
    ! The two latter cases are used together with gamma tricks
    ! 
    USE constants,      ONLY : tpi
    USE ions_base,      ONLY : nat, ntyp => nsp, ityp, tau
    USE uspp,           ONLY : nkb, vkb,  okvan, indv_ijkb0, ijtoh
    USE uspp_param,     ONLY : upf, nh, nhm, lmaxq
    USE fft_base,       ONLY : dffts
    USE gvect,          ONLY : ngm, nl, nlm, gg, g, gstart, &
                               eigts1, eigts2, eigts3, mill
    USE gvecs,          ONLY : ngms, nls, nlsm
    USE cell_base,      ONLY : tpiba, omega
    USE control_flags,  ONLY : gamma_only
    !
    IMPLICIT NONE
    !
    COMPLEX(DP),INTENT(in)    :: vc(dffts%nnr)
    ! In input I get a slice of <beta|left> and <beta|right> only for this kpoint and this band
    COMPLEX(DP),INTENT(in), OPTIONAL :: becphi_c(nkb)
    REAL(DP),   INTENT(in), OPTIONAL :: becphi_r(nkb)
    COMPLEX(DP),INTENT(inout) :: deexx(nkb)
    REAL(DP),INTENT(in)       :: xk(3), xkq(3)
    CHARACTER(LEN=1), INTENT(IN) :: flag
    !
    ! ... local variables
    INTEGER :: ig, ikb, jkb, ijkb0, ih, jh, na, nt, nij
    REAL(DP) :: fact
    COMPLEX(DP), EXTERNAL :: ZDOTC
    !
    COMPLEX(DP),ALLOCATABLE :: auxvc(:), &  ! vc in order of |g|
                               eigqts(:), aux1(:), aux2(:)
    COMPLEX(DP) :: fp, fm
    REAL(DP) :: arg
    LOGICAL :: add_complex, add_real, add_imaginary
    !
    IF(.not.okvan) RETURN
    !
    add_complex = ( flag=='c' .OR. flag=='C' )
    add_real    = ( flag=='r' .OR. flag=='R' )
    add_imaginary=( flag=='i' .OR. flag=='I' )
    IF ( .NOT. (add_complex .OR. add_real .OR. add_imaginary) ) &
       CALL errore('newdxx_g', 'called with incorrect flag: '//flag, 1 )
    IF ( .NOT. gamma_only .AND. ( add_real .OR. add_imaginary) ) &
       CALL errore('newdxx_g', 'need gamma tricks for this flag: '//flag, 2 )
    IF ( gamma_only .AND. add_complex ) &
       CALL errore('newdxx_g', 'gamma trick not good for this flag: '//flag, 3 )
    IF ( ( add_complex .AND. .NOT. PRESENT(becphi_c) ) .OR. &
         ( add_real    .AND. .NOT. PRESENT(becphi_r) ) .OR. &
         ( add_imaginary.AND..NOT. PRESENT(becphi_r) ) )    &
       CALL errore('newdxx_g', 'called with incorrect arguments', 2 )
    !
    CALL start_clock( 'newdxx' )
    !
    ALLOCATE(aux1(ngms), aux2(ngms), auxvc( ngms))
    ALLOCATE(eigqts(nat))
    !
    DO na = 1, nat
      arg = tpi* SUM( (xk(:) - xkq(:))*tau(:,na) )
      eigqts(na) = CMPLX( COS(arg), -SIN(arg), kind=DP)
    END DO
    !
    ! reindex just once at the beginning
    ! select real or imaginary part if so desired
    ! fact=2 to account for G and -G components
    !
    auxvc = (0._dp, 0._dp)
    IF ( add_complex ) THEN
       auxvc(1:ngms) = vc(nls(1:ngms) )
       fact=omega
    ELSE IF ( add_real ) THEN
       DO ig = 1, ngms
          fp = (vc(nls(ig)) + vc(nlsm(ig)))/2.0_dp
          fm = (vc(nls(ig)) - vc(nlsm(ig)))/2.0_dp
          auxvc(ig) = CMPLX( DBLE(fp), AIMAG(fm), KIND=dp)
       END DO
       fact=2.0_dp*omega
    ELSE IF ( add_imaginary ) THEN
       DO ig = 1, ngms
          fp = (vc(nls(ig)) + vc(nlsm(ig)))/2.0_dp
          fm = (vc(nls(ig)) - vc(nlsm(ig)))/2.0_dp
          auxvc(ig) = CMPLX( AIMAG(fp), -DBLE(fm), KIND=dp)
       END DO
       fact=2.0_dp*omega
    END IF 
    !
    nij = 0
    DO nt = 1, ntyp
       !
       IF ( upf(nt)%tvanp ) THEN
          !
          DO na = 1, nat
             !
             IF (ityp(na)==nt) THEN
                !
                ! ijkb0 points to the manifold of beta functions for atom na
                !
                ijkb0 = indv_ijkb0(na) 
                !
!$omp parallel do default(shared) private(ig)
                DO ig = 1, ngms
                   aux2(ig) = CONJG( auxvc(ig) ) * eigqts(na) * &
                              eigts1(mill(1,ig), na) * &
                              eigts2(mill(2,ig), na) * &
                              eigts3(mill(3,ig), na)
                END DO
!$omp end parallel do
                DO ih = 1, nh(nt)
                   ikb = ijkb0 + ih
                   aux1(:) = (0.0_dp, 0.0_dp)
                   DO jh = 1, nh(nt)
                      jkb = ijkb0 + jh
                      IF ( gamma_only ) THEN
!$omp parallel do default(shared) private(ig)
                         DO ig = 1, ngms
                            aux1(ig) = aux1(ig) + becphi_r(jkb) * &
                                 CONJG( qgm(ig,nij+ijtoh(ih,jh,nt)) )
                         ENDDO
!$omp end parallel do
                      ELSE
!$omp parallel do default(shared) private(ig)
                         DO ig = 1, ngms
                            aux1(ig) = aux1(ig) + becphi_c(jkb) * &
                                 CONJG( qgm(ig,nij+ijtoh(ih,jh,nt)) )
                         ENDDO
!$omp end parallel do
                      END IF
                   END DO
                   !
                   deexx(ikb) = deexx(ikb) + fact*ZDOTC(ngms, aux2, 1, aux1, 1)
                   IF( gamma_only .AND. gstart == 2 ) &
                        deexx(ikb) =  deexx(ikb) - fact*CONJG (aux2(1))*aux1(1)
                ENDDO
             ENDIF
             !
          ENDDO  ! nat
          !
          nij = nij + nh(nt)*(nh(nt)+1)/2
          !
       END IF
    ENDDO
    !
    DEALLOCATE( eigqts, auxvc, aux2, aux1)
    CALL stop_clock( 'newdxx' )
    !
    RETURN
    !
    !-----------------------------------------------------------------------
  END SUBROUTINE newdxx_g
  !-----------------------------------------------------------------------
  !
  !-----------------------------------------------------------------------
   SUBROUTINE add_nlxx_pot(lda, hpsi, xkp, npwp, igkp, deexx, eps_occ, exxalfa)
    !-----------------------------------------------------------------------
    !
    ! This subroutine computes some sort of EXX contribution to the non-local 
    ! part of the hamiltonian. 
    !   alpha_Ii = \int \sum_Jj Q_IJ(r) V^{i,j}_Fock <beta_J|phi_j> d^3(r)
    ! The actual contribution will be (summed outside)
    !  H = H+\sum_I |beta_I> alpha_Ii
    ! 
    USE ions_base,           ONLY : nat, ntyp => nsp, ityp
    USE uspp,                ONLY : nkb, okvan,indv_ijkb0
    USE uspp_param,          ONLY : upf, nh
    USE gvecs,               ONLY : nls
    USE wvfct,               ONLY : nbnd, npwx !, ecutwfc
    USE control_flags,       ONLY : gamma_only
    IMPLICIT NONE
    !
    ! In input I get a slice of <beta|left> and <beta|right> only for this kpoint and this band
    INTEGER,INTENT(in)        :: lda              ! leading dimension of hpsi
    COMPLEX(DP),INTENT(inout) :: hpsi(lda)!*npol)   ! the hamiltonian
    COMPLEX(DP),INTENT(in)    :: deexx(nkb)       ! \int \sum_J Q_IJ <beta_J|phi_i> d3r
    REAL(DP),INTENT(in)       :: xkp(3)           ! current k point
    REAL(DP),INTENT(in)       :: exxalfa       ! fraction of ex. exchange to add
    REAL(DP),INTENT(in)       :: eps_occ ! skip band where occupation is less than this
    INTEGER,INTENT(IN)        :: npwp, igkp(npwp)
    !
    ! ... local variables
    INTEGER :: ikb, ijkb0, ih, na, np
    INTEGER :: ig
    !
    COMPLEX(DP),ALLOCATABLE :: vkbp(:,:)  ! the <beta_I| function
    !
    CALL start_clock( 'nlxx_pot' )
    !
    IF(.not.okvan) RETURN
    !
    ALLOCATE(vkbp(npwx,nkb))
    !
    CALL init_us_2(npwp, igkp, xkp, vkbp)
    !
    DO np = 1, ntyp
      ONLY_FOR_USPP : &
      IF ( upf(np)%tvanp ) THEN
          DO na = 1, nat
            IF (ityp(na)==np) THEN
              DO ih = 1, nh(np)
                ikb = indv_ijkb0(na) + ih
                !
                IF(ABS(deexx(ikb))<eps_occ) CYCLE
                !
                IF (gamma_only) THEN
                   DO ig = 1,npwp
                      hpsi(ig) = hpsi(ig)-exxalfa*DBLE(deexx(ikb))*vkbp(ig,ikb)
                   ENDDO
                ELSE
                   DO ig = 1,npwp
                      hpsi(ig) = hpsi(ig) - exxalfa*deexx(ikb)*vkbp(ig,ikb)
                   ENDDO
                ENDIF
              ENDDO
            END IF
          ENDDO ! nat
      END IF &
      ONLY_FOR_USPP 
    ENDDO
    !
    DEALLOCATE(vkbp)
    CALL stop_clock( 'nlxx_pot' )
    !
    RETURN
    !
    !-----------------------------------------------------------------------
  END SUBROUTINE add_nlxx_pot
  !-----------------------------------------------------------------------
  !
  !------------------------------------------------------------------------
  SUBROUTINE addusxx_r(rho,becphi,becpsi)
    !------------------------------------------------------------------------
    ! This routine adds to the two wavefunctions density (in real space) 
    ! the part which is due to the US augmentation.
    ! NOTE: the density in this case is NOT real and NOT normalized to 1, 
    !       except when (bec)phi and (bec)psi are equal, or with gamma tricks.
    ! With gamma tricks: input rho must contain contributions from band 1
    ! in real part, from band 2 in imaginary part; call routine twice, with
    ! becphi=<beta|phi(1)> (real), then with becphi=-i*<beta|phi(2)> (imaginary)
    !
    USE ions_base,        ONLY : nat, ityp
    USE cell_base,        ONLY : omega
    USE fft_base,         ONLY : dffts
    USE uspp,             ONLY : okvan, nkb, ijtoh, indv_ijkb0
    USE uspp_param,       ONLY : upf, nh
    USE spin_orb,         ONLY : domag
    !
    USE realus, ONLY : tabs
    !
    IMPLICIT NONE
    !
    COMPLEX(DP),INTENT(inout) :: rho(dffts%nnr)
    COMPLEX(DP),INTENT(in)    :: becphi(nkb)
    COMPLEX(DP),INTENT(in)    :: becpsi(nkb)
    !
    INTEGER :: ia, nt, ir, irb, ih, jh, mbia
    INTEGER :: ikb, jkb, ijkb0
    !
    IF ( .not. okvan ) RETURN
    CALL start_clock( 'addusxx' )
    !
    DO ia = 1, nat
      !
      mbia = tabs(ia)%maxbox
      IF ( mbia == 0 ) CYCLE
      !
      nt = ityp(ia)
      IF ( .not. upf(nt)%tvanp ) CYCLE
        DO ih = 1, nh(nt)
        DO jh = 1, nh(nt)
            ijkb0 = indv_ijkb0(ia)
            ikb = ijkb0 + ih
            jkb = ijkb0 + jh
            !
            DO ir = 1, mbia
                irb = tabs(ia)%box(ir)
                rho(irb) = rho(irb) + tabs(ia)%qr(ir,ijtoh(ih,jh,nt)) &
                                         *CONJG(becphi(ikb))*becpsi(jkb)
            ENDDO
        ENDDO
        ENDDO
    ENDDO
    !
    CALL stop_clock( 'addusxx' )
    !
    RETURN
    !-----------------------------------------------------------------------
  END SUBROUTINE addusxx_r
  !-----------------------------------------------------------------------
  !
  !------------------------------------------------------------------------
  SUBROUTINE newdxx_r(vr,becphi,deexx)
    !------------------------------------------------------------------------
    !   This routine computes the integral of the perturbed potential with
    !   the Q function in real space
    USE cell_base,        ONLY : omega
    USE fft_base,         ONLY : dffts
    USE ions_base,        ONLY : nat, ityp
    USE uspp_param,       ONLY : upf, nh, nhm
    USE uspp,             ONLY : nkb, ijtoh, indv_ijkb0
    USE control_flags,    ONLY : tqr
    USE noncollin_module, ONLY : nspin_mag
    USE mp,               ONLY : mp_sum

    USE realus, ONLY : tabs

    IMPLICIT NONE
    ! Input: potential , output: contribution to integral
    COMPLEX(DP),INTENT(in)    :: vr(dffts%nnr)
    COMPLEX(DP),INTENT(in)    :: becphi(nkb)
    COMPLEX(DP),INTENT(inout) :: deexx(nkb)
    !Internal
    INTEGER     :: ia, ih, jh, ir, nt
    INTEGER     :: mbia
    INTEGER     :: ikb, jkb, ijkb0
    REAL(DP)    :: domega
    COMPLEX(DP) :: aux
    !
    CALL start_clock( 'newdxx' )
    domega = omega/(dffts%nr1*dffts%nr2*dffts%nr3)
    !
    DO ia = 1, nat
      !
      mbia = tabs(ia)%maxbox
      IF ( mbia == 0 ) CYCLE
      !
      nt = ityp(ia)
      IF ( .not. upf(nt)%tvanp ) CYCLE
      !
      DO ih = 1, nh(nt)
      DO jh = 1, nh(nt)
          ijkb0 = indv_ijkb0(ia)
          ikb = ijkb0 + ih
          jkb = ijkb0 + jh
          !
          aux = 0._dp
          DO ir = 1, mbia
              aux = aux + tabs(ia)%qr(ir,ijtoh(ih,jh,nt))*vr(tabs(ia)%box(ir))
          ENDDO
          deexx(ikb) = deexx(ikb) + becphi(jkb)*domega*aux
          !
      ENDDO
      ENDDO
      !
    ENDDO
    !
    CALL stop_clock( 'newdxx' )
  !------------------------------------------------------------------------
  END SUBROUTINE newdxx_r
  !------------------------------------------------------------------------
  !
!-----------------------------------------------------------------------
END MODULE us_exx
!-----------------------------------------------------------------------
