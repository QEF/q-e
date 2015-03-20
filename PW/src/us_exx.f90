!
! Copyright (C) 2013 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! Written by Lorenzo Paulatto (2012-2013) 
! Gamma-only tricks by Simon Binnie
! G-space code based on addusdens.f90 and compute_becsum.f90 
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
  LOGICAL,PARAMETER :: dovanxx = .true. ! DEBUG option
  !
  TYPE(bec_type),ALLOCATABLE :: becxx(:) ! <beta_I|phi_j,k>, with the wavefunctions from exxbuff
  ! the visible index is k; while I and J are inside bec_type

  COMPLEX(DP),ALLOCATABLE :: becxx_gamma(:,:)  ! gamma only version of becxx%r
                                               ! two bands stored per stripe 
  ! FIXME: put somewhere else (there is a copy in exx)
  REAL(DP),PARAMETER :: eps_occ  = 1.d-8 ! skip band where occupation is less than this
                                          

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
    REAL(DP),INTENT(in) :: w(m,n)
    INTEGER,INTENT(in)  :: m,n, imin, imax, i
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
    USE uspp,                ONLY : nkb, vkb,  okvan, indv_ijkb0
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
    REAL(DP),ALLOCATABLE    :: qmod(:), q(:,:), qq(:),  &! the modulus of G
                               ylmk0(:,:)  ! the spherical harmonics
    COMPLEX(DP),ALLOCATABLE :: qgm(:), aux(:), eigqts(:)
    INTEGER :: ikb, jkb, ijkb0, ih, jh, na, np, ig
    COMPLEX(DP) :: skk, becfac_c
    REAL(DP) :: arg, becfac_r
    LOGICAL :: add_complex, add_real, add_imaginary
    !
    IF(.not.(okvan .and. dovanxx)) RETURN
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
    ALLOCATE(qmod(ngms), qgm(ngms), aux(ngms))
    ALLOCATE(ylmk0(ngms, lmaxq * lmaxq))    
    ALLOCATE(qq(ngms), q(3,ngm))
    !
    DO ig = 1, ngms
      q(:,ig) = xk(:) - xkq(:) + g(:,ig)
      qq(ig)  = SUM(q(:,ig)**2)
      qmod(ig)= SQRT(qq(ig))
    ENDDO
    !
    CALL ylmr2 (lmaxq * lmaxq, ngms, q, qq, ylmk0)
    !
    DEALLOCATE(qq, q)
    ALLOCATE(eigqts(nat))
    DO na = 1, nat
      arg = tpi* SUM( (xk(:) - xkq(:))*tau(:,na) )
      eigqts(na) = CMPLX( COS(arg), -SIN(arg), kind=DP)
    END DO
    !
    DO np = 1, ntyp
      ONLY_FOR_USPP : &
      IF ( upf(np)%tvanp .and. ANY(ityp(1:nat) == np) ) THEN
        !
        DO ih = 1, nh(np)
        DO jh = 1, nh(np)
            !
            CALL qvan2(ngms, ih, jh, np, qmod, qgm, ylmk0)
            !
            ATOMS_LOOP : &
            DO na = 1, nat
              IF (ityp(na)==np) THEN
                !
                ! ijkb0 points to the manifold of beta functions for atom na
                !
                ijkb0 = indv_ijkb0(na) 
                ikb = ijkb0 + ih
                jkb = ijkb0 + jh

                IF ( add_complex ) THEN
                   becfac_c = CONJG(becphi_c(ikb))*becpsi_c(jkb)
                   DO ig = 1, ngms
                      skk = eigts1(mill(1,ig), na) * &
                            eigts2(mill(2,ig), na) * &
                            eigts3(mill(3,ig), na)
                      aux(ig) = qgm(ig)*eigqts(na)*skk*becfac_c
                   ENDDO
                   DO ig = 1,ngms
                     rhoc(nls(ig)) = rhoc(nls(ig)) + aux(ig)
                   ENDDO
                ELSE 
                   becfac_r = becphi_r(ikb)*becpsi_r(jkb)
                   DO ig = 1, ngms
                      skk = eigts1(mill(1,ig), na) * &
                            eigts2(mill(2,ig), na) * &
                            eigts3(mill(3,ig), na)
                      aux(ig) = qgm(ig)*eigqts(na)*skk*becfac_r
                   ENDDO
                   IF ( add_real ) THEN
                      DO ig = 1,ngms
                         rhoc(nls(ig)) = rhoc(nls(ig)) + aux(ig)
                      ENDDO
                      DO ig = gstart,ngms
                         rhoc(nlsm(ig)) = rhoc(nlsm(ig)) + CONJG(aux(ig))
                      ENDDO
                   ELSE IF ( add_imaginary ) THEN
                      DO ig = 1,ngms
                         rhoc(nls(ig)) = rhoc(nls(ig)) + (0.0_dp,1.0_dp)*aux(ig)
                      ENDDO
                      DO ig = gstart,ngms
                         rhoc(nlsm(ig)) = rhoc(nlsm(ig)) + &
                                          (0.0_dp,1.0_dp)*CONJG(aux(ig))
                      ENDDO
                   ENDIF
                ENDIF
                !
              END IF
            ENDDO ATOMS_LOOP  ! nat
            !
        END DO ! jh
        END DO ! ih
      END IF &
      ONLY_FOR_USPP 
    ENDDO
    !
    DEALLOCATE( ylmk0, qmod, qgm, eigqts, aux)
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
    USE uspp,           ONLY : nkb, vkb,  okvan, indv_ijkb0
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
    INTEGER :: ikb, jkb, ijkb0, ih, jh, na, np !, ijh
    INTEGER :: ig, fact
    COMPLEX(DP) :: skk
    !
    REAL(DP),ALLOCATABLE    :: qmod (:), q(:,:), qq(:), &
                               ylmk0 (:,:)  ! the spherical harmonics
    COMPLEX(DP),ALLOCATABLE :: qgm(:),    & ! the Q(r) function
                               auxvc(:), &  ! vc in order of |g|
                               eigqts(:)
    COMPLEX(DP) :: fp, fm
    REAL(DP) :: arg
    LOGICAL :: add_complex, add_real, add_imaginary
    !
    IF(.not.(okvan .and. dovanxx)) RETURN
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
    ALLOCATE(qgm(ngms), auxvc(ngms), qmod( ngms))
    ALLOCATE(ylmk0(ngms, lmaxq**2))    
    ALLOCATE(qq(ngms), q(3,ngm))
    !
    DO ig = 1, ngms
      q(:,ig) = xk(:) - xkq(:) + g(:,ig)
      qq(ig)  = SUM(q(:,ig)**2)
      qmod(ig)= SQRT(qq (ig) )
    ENDDO
    CALL ylmr2 (lmaxq * lmaxq, ngms, q, qq, ylmk0)
    !
    DEALLOCATE(qq, q)
    ALLOCATE(eigqts(nat))
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
       fact=1.0_dp
    ELSE IF ( add_real ) THEN
       DO ig = 1, ngms
          fp = (vc(nls(ig)) + vc(nlsm(ig)))/2.0_dp
          fm = (vc(nls(ig)) - vc(nlsm(ig)))/2.0_dp
          auxvc(ig) = CMPLX( DBLE(fp), AIMAG(fm), KIND=dp)
       END DO
       fact=2.0_dp
    ELSE IF ( add_imaginary ) THEN
       DO ig = 1, ngms
          fp = (vc(nls(ig)) + vc(nlsm(ig)))/2.0_dp
          fm = (vc(nls(ig)) - vc(nlsm(ig)))/2.0_dp
          auxvc(ig) = CMPLX( AIMAG(fp), -DBLE(fm), KIND=dp)
       END DO
       fact=2.0_dp
    END IF 
    !
    DO np = 1, ntyp
      ONLY_FOR_USPP : &
      IF ( upf(np)%tvanp ) THEN
        DO ih = 1, nh(np)
          DO jh = 1, nh(np)
            !
            CALL qvan2(ngms, ih, jh, np, qmod, qgm, ylmk0)
            !
            ATOMS_LOOP : &
            DO na = 1, nat
              IF (ityp(na)==np) THEN
                !
                ijkb0 = indv_ijkb0(na)
                ikb = ijkb0 + ih
                jkb = ijkb0 + jh
                !
                IF(gamma_only) THEN
                   DO ig = 1, ngms
                      skk = eigts1(mill(1,ig), na) * &
                            eigts2(mill(2,ig), na) * &
                            eigts3(mill(3,ig), na)
                      ! \sum_J Q_IJ V_F
                      deexx(ikb) = deexx(ikb) + becphi_r(jkb)*auxvc(ig)*fact &
                                      * omega*CONJG(eigqts(na)*skk*qgm(ig)) 
                   ENDDO
                   !
                   IF(gstart==2) deexx(ikb) = deexx(ikb) - becphi_r(jkb)* &
                               auxvc(1)*omega*CONJG(eigqts(na)*skk*qgm(1))
                ELSE
                   DO ig = 1, ngms
                      skk = eigts1(mill(1,ig), na) * &
                            eigts2(mill(2,ig), na) * &
                            eigts3(mill(3,ig), na)
                      ! \sum_J Q_IJ V_F
                      deexx(ikb) = deexx(ikb) + becphi_c(jkb)*auxvc(ig)*fact &
                                      * omega*CONJG(eigqts(na)*skk*qgm(ig)) 
                   ENDDO
                ENDIF
                !
              END IF
            ENDDO ATOMS_LOOP ! nat
          ENDDO ! jh
        ENDDO ! ih
      END IF &
      ONLY_FOR_USPP 
    ENDDO
    !
    DEALLOCATE( ylmk0, qmod, qgm, auxvc, eigqts)
    CALL stop_clock( 'newdxx' )
    !
    RETURN
    !
    !-----------------------------------------------------------------------
  END SUBROUTINE newdxx_g
  !-----------------------------------------------------------------------
  !
!   !----------------------------------------------------------------------
!   SUBROUTINE addusxx_force(forcenl)
!     !----------------------------------------------------------------------
!     !
!     !   This routine computes the contribution to atomic forces due
!     !   to the dependence of the Q function on the atomic position.
!     !   On output: the contribution is added to forcenl
!     !
!     USE kinds,      ONLY : DP
!     USE ions_base,  ONLY : nat, ntyp => nsp, ityp
!     USE cell_base,  ONLY : omega, tpiba
!     USE fft_base,   ONLY : dfftp
!     USE gvect,      ONLY : ngm, nl, nlm, gg, g, eigts1, eigts2, eigts3, mill
!     USE scf,        ONLY : v, vltot
!     USE uspp,       ONLY : becsum, okvan
!     USE uspp_param, ONLY : upf, lmaxq, nh, nhm
!     USE mp_bands,   ONLY : intra_bgrp_comm
!     USE mp,         ONLY : mp_sum
!     USE noncollin_module,   ONLY : nspin_mag
!     USE control_flags,      ONLY : gamma_only
!     USE fft_interfaces,     ONLY : fwfft
!     !
!     IMPLICIT NONE
!     !
!     REAL(DP) :: forcenl (3, nat)
!     ! local variables
!     INTEGER :: ig, ir, dim, nt, ih, jh, ijh, ipol, is, na
!     COMPLEX(DP):: cfac
!     REAL(DP) :: fact, ddot
!     ! work space
!     COMPLEX(DP),ALLOCATABLE :: aux(:,:), aux1(:,:), vg(:), qgm(:), eigqts(:)
!     REAL(DP),ALLOCATABLE :: ddeeq(:,:,:,:), qmod(:), ylmk0(:,:)
! 
!     !
!     if (.not.okvan) return
!     !
!     DO ig = 1, ngms
!       q(:,ig) = xk(:) - xkq(:) + g(:,ig)
!       qq(ig)  = SUM(q(:,ig)**2)
!       qmod(ig)= SQRT(qq (ig) )
!     ENDDO
!     !
!     ALLOCATE(eigqts(nat))
!     DO na = 1, nat
!       arg = tpi* SUM( (xk(:) - xkq(:))*tau(:,na) )
!       eigqts(na) = CMPLX( COS(arg), -SIN(arg), kind=DP)
!     END DO
!     !
!     IF (gamma_only) THEN
!       fact = 2.d0
!     ELSE
!       fact = 1.d0
!     ENDIF
!     ALLOCATE (aux(ngm,nspin_mag))    
!     !
!     ! fourier transform of the total effective potential
!     !
!     ALLOCATE (vg(dfftp%nnr))    
!     DO is = 1, nspin_mag
!       IF (nspin_mag.eq.4.and.is.ne.1) then
!           vg (:) = v%of_r(:,is)
!       ELSE
!           vg (:) = vltot (:) + v%of_r (:, is)
!       ENDIF
!       CALL fwfft ('Dense', vg, dfftp)
!       aux (:, is) = vg (nl (:) ) * tpiba * (0.d0, -1.d0)
!     ENDDO
!     DEALLOCATE (vg)
!     !
!     ALLOCATE (aux1(ngm,3))    
!     ALLOCATE (ddeeq( 3, (nhm*(nhm+1))/2,nat,nspin_mag))    
!     ALLOCATE (qgm( ngm))
!     ALLOCATE (qmod( ngm))    
!     ALLOCATE (ylmk0(ngm,lmaxq*lmaxq))    
!     !
!     ddeeq(:,:,:,:) = 0.d0
!     !
!     CALL ylmr2 (lmaxq * lmaxq, ngm, g, gg, ylmk0)
!     !
!     qmod (:) = sqrt (gg (:) )
!     !
!     ! here we compute the integral Q*V for each atom,
!     !       I = sum_G i G_a exp(-iR.G) Q_nm v^*
!     !
!     DO nt = 1, ntyp
!       IF ( upf(nt)%tvanp ) then
!           ijh = 1
!           DO ih = 1, nh (nt)
!             DO jh = ih, nh (nt)
!                 call qvan2 (ngm, ih, jh, nt, qmod, qgm, ylmk0)
!                 DO na = 1, nat
!                   IF (ityp (na) == nt) then
!                       !
!                       ! The product of potential, structure factor and iG
!                       !
!                       DO is = 1, nspin_mag
!                         DO ig = 1, ngm
!                             cfac = aux(ig, is) * eigqts(na) * &
!                                          CONJG(eigts1(mill(1,ig), na) *&
!                                                eigts2(mill(2,ig), na) *&
!                                                eigts3(mill(3,ig), na) )
!                             aux1(ig, 1) = g(1, ig) * cfac
!                             aux1(ig, 2) = g(2, ig) * cfac
!                             aux1(ig, 3) = g(3, ig) * cfac
!                         ENDDO
!                         !
!                         !    and the product with the Q functions
!                         !    G=0 term gives no contribution
!                         !
!                         DO ipol = 1, 3
!                             ddeeq (ipol, ijh, na, is) = omega * fact * &
!                                 ddot (2 * ngm, aux1(1, ipol), 1, qgm, 1)
!                         ENDDO
!                       ENDDO
!                   ENDIF
!                 ENDDO
!                 ijh = ijh + 1
!             ENDDO
!           ENDDO
!       ENDIF
! 
!     ENDDO
! 
!     call mp_sum ( ddeeq, intra_bgrp_comm )
!     !
!     DO is = 1, nspin_mag
!       DO na = 1, nat
!           nt = ityp (na)
!           dim = (nh (nt) * (nh (nt) + 1) ) / 2
!           DO ipol = 1, 3
!             DO ir = 1, dim
!                 forcenl(ipol, na) = forcenl(ipol, na) + &
!                     ddeeq(ipol, ir, na, is) * becsum(ir, na, is)
!             ENDDO
!           ENDDO
!       ENDDO
!     ENDDO
!     !
!     DEALLOCATE(ylmk0,qgm,qmod,ddeeq,aux1,aux,eigqts)
!     RETURN
!     !-----------------------------------------------------------------------
!   END SUBROUTINE addusxx_force
!   !-----------------------------------------------------------------------
  !
  !-----------------------------------------------------------------------
   SUBROUTINE add_nlxx_pot(lda, hpsi, xkp, npwp, igkp, deexx, exxalfa)
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
    IF(.not.(okvan .and. dovanxx)) RETURN
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
  !------------------------------------------------------------------------
  END SUBROUTINE newdxx_r
  !------------------------------------------------------------------------
  !
!-----------------------------------------------------------------------
END MODULE us_exx
!-----------------------------------------------------------------------
