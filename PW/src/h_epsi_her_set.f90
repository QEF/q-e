!
! Copyright (C) 2005 Paolo Umari
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

!-----------------------------------------------------------------------
SUBROUTINE h_epsi_her_set( pdir, e_field )
  !-----------------------------------------------------------------------
  !! This subroutine builds the hermitian operators  w_k w_k*, 
  !! (as in Souza,et al.  PRB B 69, 085106 (2004)).
  !
  !! Wavefunctions from previous iteration are read into \(\textrm{evcel}\);
  !! Spin polarized systems are supported only with fixed occupations.
  !
  USE noncollin_module,   ONLY: noncolin, npol
  USE spin_orb,           ONLY: lspinorb
  USE kinds,              ONLY: DP
  USE wvfct,              ONLY: npwx, nbnd
  USE ldaU,               ONLY: lda_plus_u
  USE lsda_mod,           ONLY: current_spin, nspin
  USE scf,                ONLY: vrs  
  USE gvect
  USE fft_base,           ONLY: dfftp
  USE uspp,               ONLY: okvan, nkb, vkb
  USE uspp_param,         ONLY: upf, nh, nhm, nbetam, lmaxq
  USE bp,                 ONLY: nppstr_3d, fact_hepsi, evcel, evcp=>evcelp, &
                                evcm=>evcelm, mapgp_global, mapgm_global, nx_el
  USE klist
  USE cell_base,          ONLY: at, alat, tpiba, omega, bg
  USE ions_base,          ONLY: ityp, tau, nat,ntyp => nsp
  USE io_files,           ONLY: iunwfc, nwordwfc, iunefieldm, iunefieldp
  USE buffers,            ONLY: get_buffer, save_buffer
  USE constants,          ONLY: e2, pi, tpi, fpi
  USE fixed_occ
  USE mp,                 ONLY: mp_sum
  USE mp_bands,           ONLY: intra_bgrp_comm
  USE becmod,             ONLY: bec_type, becp, calbec,ALLOCATE_bec_type, &
                                deALLOCATE_bec_type
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: pdir
  !! the direction on which the polarization is calculated
  REAL(DP) :: e_field
  !! the electric field along pdir
  !
  ! ... local variables
  !
   COMPLEX(DP), ALLOCATABLE  :: evct(:,:) !for temporary wavefunctios
   INTEGER :: i, ipol
   INTEGER :: igk1(npwx)
   INTEGER :: igk0(npwx)
   INTEGER :: ig
   INTEGER :: info
   INTEGER :: is
   INTEGER :: iv
   INTEGER :: ivpt(nbnd)
   INTEGER :: j
   INTEGER :: jkb
   INTEGER :: jkb_bp
   INTEGER :: jkb1
   INTEGER :: jv
   INTEGER :: m
   INTEGER :: mb
   INTEGER :: mk1
   INTEGER :: mk2
   INTEGER :: mk3
   INTEGER :: n1
   INTEGER :: n2
   INTEGER :: n3
   INTEGER :: na
   INTEGER :: nb
   INTEGER :: ng
   INTEGER :: nhjkb
   INTEGER :: nhjkbm
   INTEGER :: nkbtona(nkb)
   INTEGER :: nkbtonh(nkb)
   INTEGER :: np
   INTEGER :: npw1
   INTEGER :: npw0
   !INTEGER :: nstring
   INTEGER :: nt
   INTEGER :: ik_stringa!k-point index inside string
   REAL(DP) :: dk(3)
   REAL(DP) :: dkm(3)! -dk
   REAL(DP) :: dkmod
   REAL(DP) :: eps
   REAL(DP) :: fac
   REAL(DP) :: gpar(3)
   REAL(DP) :: gtr(3)
   !REAL(DP) :: gvec
   REAL(DP), ALLOCATABLE :: ln(:,:,:)
   REAL(DP), ALLOCATABLE  :: ln0(:,:,:)!map g-space global to g-space k-point dependent
   REAL(DP) :: qrad_dk(nbetam,nbetam,lmaxq,ntyp)
   REAL(DP) :: ylm_dk(lmaxq*lmaxq)
   COMPLEX(DP), ALLOCATABLE :: aux(:)
   COMPLEX(DP), ALLOCATABLE  :: aux0(:)
   ! Also for noncollinear calculation
   COMPLEX(DP), ALLOCATABLE :: aux_2(:)
   COMPLEX(DP), ALLOCATABLE :: aux0_2(:), aux0vec(:,:), aux1vec(:,:)
   !
   COMPLEX(DP) :: cdet(2)
   COMPLEX(DP) :: cdwork(nbnd)
   COMPLEX(DP), ALLOCATABLE :: mat(:,:)
   COMPLEX(DP) :: pref
   COMPLEX(DP) :: q_dk(nhm,nhm,ntyp)
   COMPLEX(DP) :: q_dkp(nhm,nhm,ntyp)!to store the terms T^dagger e^(iGx) T
   COMPLEX(DP) :: struc(nat)
   !
   COMPLEX(DP) :: sca,sca1
   COMPLEX(DP) :: ps(nkb,nbnd*npol)
   COMPLEX(DP) :: matbig(nks,nbnd,nbnd)
   INTEGER :: mdone(nks)
   INTEGER :: ijkb0, ibnd, jh, ih, ikb, ik, ikk
   !
   LOGICAL, ALLOCATABLE :: l_cal(:) ! flag for empty/occupied states
   INTEGER, ALLOCATABLE  :: map_g(:)
   TYPE(bec_type) :: becp_bp,becp0
   !
   REAL(DP) :: dkfact
   LOGICAL  :: l_para! if true new parallel treatment
   COMPLEX(kind=DP), ALLOCATABLE :: aux_g(:)
   COMPLEX(kind=DP), ALLOCATABLE :: aux_g_2(:) ! non-collinear case
   !
   COMPLEX(DP), ALLOCATABLE :: q_dk_so(:,:,:,:), q_dkp_so(:,:,:,:)
   !
   !
   !  --- Define a small number ---
   eps = 0.000001_DP
   !
   ALLOCATE( mat(nbnd,nbnd) )
   !
   IF(ABS(e_field)<eps) RETURN
   !
   CALL start_clock('h_epsi_set')
   !
   !
   IF (okvan) THEN
      !
      CALL ALLOCATE_bec_type( nkb, nbnd, becp0 )
      CALL ALLOCATE_bec_type( nkb, nbnd, becp_bp )
      !
      IF (lspinorb) ALLOCATE( q_dk_so(nhm,nhm,4,ntyp) )
      IF (lspinorb) ALLOCATE( q_dkp_so(nhm,nhm,4,ntyp) )
      !
   ENDIF
   !
   !-------------------------------------------------------------------------!
   !                               INITIALIZATIONS                           !
   !-------------------------------------------------------------------------!
   !
   IF (pdir==3) THEN
      l_para = .FALSE.
   ELSE
      l_para = .TRUE.
   ENDIF
   !
   ALLOCATE( evct(npwx*npol,nbnd) )
   ALLOCATE( map_g(npwx) )
   !
   ALLOCATE( ln(-dfftp%nr1:dfftp%nr1,-dfftp%nr2:dfftp%nr2,-dfftp%nr3:dfftp%nr3), &
             ln0(-dfftp%nr1:dfftp%nr1,-dfftp%nr2:dfftp%nr2,-dfftp%nr3:dfftp%nr3) )
   !
   ALLOCATE( aux(ngm), aux0(ngm) )
   IF (noncolin) ALLOCATE( aux_2(ngm), aux0_2(ngm) )
   !
   ALLOCATE( l_cal(nbnd) )
   !
   ! --- Recalculate FFT correspondence (see ggen.f90) ---
   DO ng = 1, ngm
      mk1 = NINT( g(1,ng)*at(1,1) + g(2,ng)*at(2,1) + g(3,ng)*at(3,1) )
      mk2 = NINT( g(1,ng)*at(1,2) + g(2,ng)*at(2,2) + g(3,ng)*at(3,2) )
      mk3 = NINT( g(1,ng)*at(1,3) + g(2,ng)*at(2,3) + g(3,ng)*at(3,3) )
      ln(mk1,mk2,mk3) = ng
   ENDDO
   ! --- Find vector along strings ---
   IF (nppstr_3d(pdir) /= 1) THEN
      !
      gpar(1) = ( xk(1,nx_el(nppstr_3d(pdir),pdir))-xk(1,nx_el(1,pdir)) ) * &
                  DBLE(nppstr_3d(pdir)) / DBLE(nppstr_3d(pdir)-1)
      gpar(2) = ( xk(2,nx_el(nppstr_3d(pdir),pdir))-xk(2,nx_el(1,pdir)) ) * &
                  DBLE(nppstr_3d(pdir)) / DBLE(nppstr_3d(pdir)-1)
      gpar(3) = ( xk(3,nx_el(nppstr_3d(pdir),pdir))-xk(3,nx_el(1,pdir)) ) * &
                  DBLE(nppstr_3d(pdir)) / DBLE(nppstr_3d(pdir)-1)
   ELSE
      !
      gpar = bg(:,pdir)
      !
   ENDIF
   !
   ! gvec = dsqrt(gpar(1)**2+gpar(2)**2+gpar(3)**2)*tpiba
   !
   matbig = (0.0d0,0.0d0)
   mdone = 0
   !
   CALL factor_a( pdir, at, dkfact )
   !
   dkfact = tpiba / dkfact / DBLE(nppstr_3d(pdir))
   !
   ! determines the spin polarization
   DO ik = 1, nks
      !  --- Find vector between consecutive points in strings ---
      dk = gpar/nppstr_3d(pdir)
      dkmod = SQRT(DOT_PRODUCT(dk,dk))*tpiba
      dkm = -dk
      !
      CALL get_buffer( evcel, nwordwfc, iunwfc, nx_el(ik,pdir) )
      !
      IF (nspin==2) THEN
         IF (ik <= nks/2) THEN
            is = 1
         ELSE
            is = 2
         ENDIF
      ELSE
         is = 1
      ENDIF
      !
      ! l_cal(n) = .true./.false. if n-th state is occupied/empty
      DO nb = 1, nbnd
         IF ( nspin == 2 .AND. tfixed_occ) THEN
            l_cal(nb) = ( f_inp(nb,is) /= 0.0_DP )
         ELSE
            IF (noncolin) THEN
                l_cal(nb) = ( nb <= NINT( nelec) )
             ELSE
                l_cal(nb) = ( nb <= NINT( nelec/2.0_DP ) )
             ENDIF
         ENDIF
      ENDDO
      !
      ik_stringa = MOD(ik-1, nppstr_3d(pdir)) + 1
      !
      IF (okvan) THEN
         !  --- Initialize arrays ---
         jkb_bp = 0
         DO nt = 1, ntyp
            DO na = 1, nat
               IF (ityp(na) == nt) THEN
                  DO i = 1, nh(nt)
                     jkb_bp = jkb_bp + 1
                     nkbtona(jkb_bp) = na
                     nkbtonh(jkb_bp) = i
                  ENDDO
               ENDIF
            ENDDO
         ENDDO
      ENDIF
      !
      !-------------------------------------------------------------------------!
      !           electronic polarization: set values for k-points strings      !
      !-------------------------------------------------------------------------!
      !
      ! ... Calculates fact factor.
      ! Electronic charge is sqrt(2.) (Rydberg units).
      ! The factor (-i)/2 comes form operator Im.
      !
      fact_hepsi( nx_el(ik,pdir), pdir ) = (0.d0,-1.d0) * e_field * &
                                           DSQRT(2.d0) / 2.d0 / dkfact
      !
      evcm(:,:,pdir) = (0.d0,0.d0)
      evcp(:,:,pdir) = (0.d0,0.d0)
      !
      IF (okvan) THEN
         !-------------------------------------------------------------------------!
         !                  electronic polarization: structure factor              !
         !-------------------------------------------------------------------------!
         !
         !  --- Calculate structure factor e^{-i dk*R} ---
         DO na = 1, nat
            fac = (dk(1)*tau(1,na)+dk(2)*tau(2,na)+dk(3)*tau(3,na))*tpi 
            struc(na) = CMPLX( COS(fac), -SIN(fac), kind=DP )
         ENDDO
         !
         !-------------------------------------------------------------------------!
         !                     electronic polarization: form factor                !
         !-------------------------------------------------------------------------!

         !  --- Calculate Bessel transform of Q_ij(|r|) at dk [Q_ij^L(|r|)] ---
         CALL calc_btq( dkmod, qrad_dk, 0 )
         !
         !  --- Calculate the q-space real spherical harmonics at dk [Y_LM] --- 
         dkmod = dk(1)**2 + dk(2)**2 + dk(3)**2
         CALL ylmr2( lmaxq*lmaxq, 1, dk, dkmod, ylm_dk )
         !
         !  --- Form factor: 4 pi sum_LM c_ij^LM Y_LM(Omega) Q_ij^L(|r|) ---
         q_dk = (0.d0,0.d0)
         DO np = 1, ntyp
            IF ( upf(np)%tvanp ) THEN
               DO iv = 1, nh(np)
                  DO jv = iv, nh(np)
                     CALL qvan3( iv, jv, np, pref, ylm_dk, qrad_dk )
                     q_dk(iv,jv,np) = omega*pref
                     q_dk(jv,iv,np) = omega*pref
                  ENDDO
               ENDDO
            ENDIF
         ENDDO
         IF (lspinorb) CALL transform_qq_so( q_dk, q_dk_so )
         !
         !  --- Calculate the q-space real spherical harmonics at -dk [Y_LM] --- 
         !
         dkmod = dkm(1)**2 + dkm(2)**2 + dkm(3)**2
         CALL ylmr2( lmaxq*lmaxq, 1, dkm, dkmod, ylm_dk )
         !
         !  --- Form factor: 4 pi sum_LM c_ij^LM Y_LM(Omega) Q_ij^L(|r|) ---
         !
         q_dkp = (0.d0,0.d0)
         DO np = 1, ntyp
            IF ( upf(np)%tvanp ) THEN
               DO iv = 1, nh(np)
                  DO jv = iv, nh(np)
                     CALL qvan3( iv, jv, np, pref, ylm_dk, qrad_dk )
                     q_dkp(iv,jv,np) = omega*pref
                     q_dkp(jv,iv,np) = omega*pref
                  ENDDO
               ENDDO
            ENDIF
         ENDDO
         !
         IF (lspinorb) CALL transform_qq_so( q_dkp, q_dkp_so )
         !
      ENDIF
      !
      !-------------------------------------------------------------------------!
      !                   electronic polarization: strings phases               !
      !-------------------------------------------------------------------------!
      !    
      ! calculate the term  S-1(k,k-1)
      !       
      IF (ik_stringa /= 1) THEN
         !
         ikk = nx_el(ik-1,pdir)
         npw0 = ngk(ikk)
         igk0(:) = igk_k(:,ikk)
         !
         CALL get_buffer( evct, nwordwfc, iunwfc, nx_el(ik-1,pdir) )
         !
         ! --- Calculate dot products between wavefunctions
         !
         ! --- Dot wavefunctions and betas for PREVIOUS k-point ---
         IF (okvan) THEN
            CALL init_us_2( npw0, igk0, xk(1,nx_el(ik-1,pdir)), vkb )
            CALL calbec( npw0, vkb, evct, becp0 )
         ENDIF
         !
         ! --- Dot wavefunctions and betas for CURRENT k-point ---
         ikk = nx_el(ik,pdir)
         npw1 = ngk(ikk)
         igk1(:) = igk_k(:,ikk)
         !
         !  --- Recalculate FFT correspondence (see ggen.f90) ---
         !
         ln0 = 0 !set array to 0
         DO ig = 1, npw1
            mk1 = NINT(g(1,igk1(ig))*at(1,1)+g(2,igk1(ig))*at(2,1)+g(3,igk1(ig))*at(3,1))
            mk2 = NINT(g(1,igk1(ig))*at(1,2)+g(2,igk1(ig))*at(2,2)+g(3,igk1(ig))*at(3,2))
            mk3 = NINT(g(1,igk1(ig))*at(1,3)+g(2,igk1(ig))*at(2,3)+g(3,igk1(ig))*at(3,3))
            ln0(mk1,mk2,mk3) = ig
         ENDDO
         !
         IF (okvan) THEN
            CALL init_us_2( npw1, igk1, xk(1,nx_el(ik,pdir)), vkb )
            CALL calbec( npw1, vkb, evcel, becp_bp )
         ENDIF
         ! --- Matrix elements calculation ---
         !
         IF (mdone(nx_el(ik,pdir)) == 0) THEN
            !
            mat = (0.d0,0.d0)
            !
            DO nb = 1, nbnd
                DO mb = 1, nbnd
                  IF ( .NOT. l_cal(nb) .OR. .NOT. l_cal(mb) ) THEN
                     IF ( nb == mb )  mat(nb,mb)=0.d0
                  ELSE
                     aux = (0.d0,0.d0)
                     aux0 = (0.d0,0.d0)
                     IF (noncolin) THEN
                        aux_2 = (0.d0,0.d0)
                        aux0_2 = (0.d0,0.d0)
                     ENDIF
                     DO ig = 1, npw1
                        aux0(igk1(ig)) = evcel(ig,nb)
                        IF (noncolin) aux0_2(igk1(ig)) = evcel(ig+npwx,nb)
                     ENDDO
                     !
                     DO ig = 1, npw0
                        aux(igk0(ig)) = evct(ig,mb)
                        IF (noncolin) aux_2(igk0(ig)) = evct(ig+npwx,mb)
                     ENDDO
                     !
                     mat(nb,mb) = dot_product(aux0(1:ngm),aux(1:ngm))
                     IF (noncolin) mat(nb,mb) = mat(nb,mb) + dot_product(aux0_2(1:ngm),aux_2(1:ngm))
                     ! --- Calculate the augmented part: ij=KB projectors, ---
                     ! --- R=atom index: SUM_{ijR} q(ijR) <u_nk|beta_iR>   ---
                     ! --- <beta_jR|u_mk'> e^i(k-k')*R =                   ---
                     ! --- also <u_nk|beta_iR>=<psi_nk|beta_iR> = becp^*   ---
                  ENDIF
               ENDDO
            ENDDO
            !
            CALL mp_sum( mat, intra_bgrp_comm )
            !
            DO nb=1,nbnd
               IF(.NOT. l_cal(nb)) mat(nb,nb) = 1.d0
            ENDDO
            !
            DO nb = 1, nbnd
               DO mb = 1, nbnd
                  IF ( l_cal(nb) .AND. l_cal(mb) ) THEN
                     IF (okvan) THEN
                        pref = (0.d0,0.d0)
                        DO jkb = 1, nkb
                           nhjkb = nkbtonh(jkb)
                           na = nkbtona(jkb)
                           np = ityp(na)
                           nhjkbm = nh(np)
                           jkb1 = jkb - nhjkb
                           IF (lspinorb) THEN
                              DO j = 1,nhjkbm
                                 pref = pref + CONJG(becp_bp%nc(jkb,1,nb))* becp0%nc(jkb1+j,1,mb) &
                                               *q_dkp_so(nhjkb,j,1,np)*CONJG(struc(na))   
                                 pref = pref + CONJG(becp_bp%nc(jkb,1,nb))* becp0%nc(jkb1+j,2,mb) &
                                               *q_dkp_so(nhjkb,j,2,np)*CONJG(struc(na)) 
                                 pref = pref + CONJG(becp_bp%nc(jkb,2,nb))* becp0%nc(jkb1+j,1,mb) &
                                               *q_dkp_so(nhjkb,j,3,np)*CONJG(struc(na)) 
                                 pref = pref + CONJG(becp_bp%nc(jkb,2,nb))* becp0%nc(jkb1+j,2,mb) &
                                               *q_dkp_so(nhjkb,j,4,np)*CONJG(struc(na))
                              ENDDO
                           ELSE
                              DO j = 1,nhjkbm
                                 pref = pref + CONJG(becp_bp%k(jkb,nb))*becp0%k(jkb1+j,mb) &
                                               *q_dkp(nhjkb,j,np)*CONJG(struc(na))
                              ENDDO
                           ENDIF
                        ENDDO
                        mat(nb,mb) = mat(nb,mb) + pref
                     ENDIF
                  ENDIF
               ENDDO
            ENDDO
            ! --- Calculate matrix inverse ---
            CALL zgefa( mat, nbnd, nbnd, ivpt, info )
            CALL errore( 'h_epsi_her_set','error in zgefa', ABS(info) )
            CALL zgedi( mat, nbnd, nbnd, ivpt, cdet, cdwork, 1 )
            matbig( nx_el(ik,pdir),:,: ) = mat
            mdone(nx_el(ik,pdir)) = 1
            !
         ELSE
            !
            mat = matbig( nx_el(ik,pdir),:,: )
            !
         ENDIF
         !
         DO nb = 1, nbnd
            DO mb = 1, nbnd
               IF(.NOT.l_cal(nb).OR. .NOT. l_cal(mb)) mat(mb,nb) = (0.d0,0.d0)
            ENDDO
         ENDDO
         !
         ! mat=S^-1(k,k-1)
         DO ig = 1, npw0
            gtr(1) = g(1,igk0(ig))
            gtr(2) = g(2,igk0(ig))         
            gtr(3) = g(3,igk0(ig))         
            !     --- Find crystal coordinates of gtr, n1,n2,n3 ---
            !     --- and the position ng in the ngm array      ---
            IF (gtr(1)**2+gtr(2)**2+gtr(3)**2 <= gcutm) THEN
               n1 = NINT( gtr(1)*at(1,1) + gtr(2)*at(2,1) + &
                          gtr(3)*at(3,1) )
               n2 = NINT( gtr(1)*at(1,2) + gtr(2)*at(2,2) + &
                          gtr(3)*at(3,2) )
               n3 = NINT( gtr(1)*at(1,3) + gtr(2)*at(2,3) + &
                          gtr(3)*at(3,3) )
               ng = ln0(n1,n2,n3)
               IF (ng > 0) THEN
                  DO m = 1, nbnd
                     DO nb = 1, nbnd
                        evcm(ng,m,pdir) = evcm(ng,m,pdir) + mat(nb,m)*evct(ig,nb)
                        IF (noncolin) evcm(ng+npwx,m,pdir) = evcm(ng+npwx,m,pdir) + &
                                                             mat(nb,m)*evct(ig+npwx,nb)
                     ENDDO
                  ENDDO
               ENDIF
            ENDIF
         ENDDO
         !
         ! ... add US terms into evcm
         ! ... calculate |beta_(ik,na,ih)>Q_dkp(na,ih,ij)<|beta_(ik-1,na,ih)| 
         !
         IF (okvan) THEN
            !
            evct(:,:) =  (0.d0, 0.d0)
            ps(:,:) = (0.d0, 0.d0)
            ijkb0 = 0
            DO nt = 1, ntyp
               DO na = 1, nat
                  IF (ityp(na) == nt) THEN
                     DO ibnd = 1, nbnd
                        DO jh = 1, nh(nt)
                           jkb = ijkb0 + jh
                           DO ih = 1, nh(nt)
                              ikb = ijkb0 + ih
                              IF (lspinorb) THEN
                                 ps(ikb, (ibnd-1)*npol+1) = ps(ikb,(ibnd-1)*npol+1 ) + &
                                     q_dkp_so(ih,jh,1,ityp(na))*CONJG(struc(na))* becp0%nc(jkb,1,ibnd)
                                 ps(ikb, (ibnd-1)*npol+1) = ps(ikb,(ibnd-1)*npol+1 ) + &
                                     q_dkp_so(ih,jh,2,ityp(na))*CONJG(struc(na))* becp0%nc(jkb,2,ibnd)
                                 ps(ikb, (ibnd-1)*npol+2) = ps(ikb,(ibnd-1)*npol+2 ) + &
                                     q_dkp_so(ih,jh,3,ityp(na))*CONJG(struc(na))* becp0%nc(jkb,1,ibnd)
                                 ps(ikb, (ibnd-1)*npol+2) = ps(ikb,(ibnd-1)*npol+2 ) + &
                                     q_dkp_so(ih,jh,4,ityp(na))*CONJG(struc(na))* becp0%nc(jkb,2,ibnd)
                              ELSE
                                 ps(ikb, ibnd) = ps(ikb, ibnd) + &
                                          q_dkp(ih,jh,ityp(na))*CONJG(struc(na))* becp0%k(jkb,ibnd)
                              ENDIF
                           ENDDO
                        ENDDO
                     ENDDO
                     ijkb0 = ijkb0 + nh(nt)
                  ENDIF
               ENDDO
            ENDDO
            !
            ! vkb is relative to the last ik read
            CALL ZGEMM( 'N', 'N', npw1, nbnd*npol , nkb, (1.d0, 0.d0) , vkb, &
                        npwx, ps, nkb, (1.d0, 0.d0) , evct, npwx )
            !
            DO m = 1, nbnd
               DO nb = 1, nbnd
                  DO ig = 1, npw1
                     evcm(ig,m,pdir) = evcm(ig,m,pdir) + mat(nb,m)*evct(ig,nb)
                  ENDDO
                  IF (noncolin) THEN
                     DO ig = 1, npw1
                        evcm(ig+npwx,m,pdir)=evcm(ig+npwx,m,pdir) + mat(nb,m)*evct(ig+npwx,nb)
                     ENDDO
                  ENDIF
               ENDDO
            ENDDO
            !
         ENDIF
         !  --- End of dot products between wavefunctions and betas ---
      ELSE !(ik_stringa == 1)
         !
         ikk = nx_el(ik+nppstr_3d(pdir)-1,pdir)
         npw0 = ngk(ikk)
         igk0(:) = igk_k(:,ikk)
         !
         CALL get_buffer( evct, nwordwfc, iunwfc, nx_el(ik+nppstr_3d(pdir)-1,pdir) )
         !        
         ! --- Calculate dot products between wavefunctions
         !
         ! --- Dot wavefunctions and betas for PREVIOUS k-point ---
         !
         IF (okvan) THEN
            CALL init_us_2( npw0, igk0, xk(1,nx_el(ik+nppstr_3d(pdir)-1,pdir)), vkb )
            CALL calbec( npw0, vkb, evct, becp0 )
         ENDIF
         ! --- Dot wavefunctions and betas for CURRENT k-point ---
         !
         ikk = nx_el(ik,pdir)
         npw1 = ngk(ikk)
         igk1(:)= igk_k(:,ikk)
         !
         ! --- Recalculate FFT correspondence (see ggen.f90) ---
         !
         IF (.NOT.l_para) THEN
            ln0 = 0 !set to 0
            DO ig = 1, npw1
               mk1 = NINT(g(1,igk1(ig))*at(1,1)+g(2,igk1(ig))*at(2,1)+g(3,igk1(ig))*at(3,1))
               mk2 = NINT(g(1,igk1(ig))*at(1,2)+g(2,igk1(ig))*at(2,2)+g(3,igk1(ig))*at(3,2))
               mk3 = NINT(g(1,igk1(ig))*at(1,3)+g(2,igk1(ig))*at(2,3)+g(3,igk1(ig))*at(3,3))
               ln0(mk1,mk2,mk3) = ig
            ENDDO
         ENDIF
         !
         IF (okvan) THEN
            CALL init_us_2( npw1, igk1, xk(1,nx_el(ik,pdir)), vkb )
            CALL calbec( npw1, vkb, evcel, becp_bp )
         ENDIF
         !
         ! --- Matrix elements calculation ---
         !
         IF (mdone(nx_el(ik,pdir)) == 0) THEN
            !
            mat = (0.d0,0.d0)
            !
            IF (.NOT. l_para) THEN
               !
               map_g(:) = 0
               !
               DO ig = 1, npw0
                  !  --- If k'=k+G_o, the relation psi_k+G_o (G-G_o) ---
                  !  --- = psi_k(G) is used, gpar=G_o, gtr = G-G_o   ---
                  gtr(1) = g(1,igk0(ig)) + gpar(1)
                  gtr(2) = g(2,igk0(ig)) + gpar(2) 
                  gtr(3) = g(3,igk0(ig)) + gpar(3) 
                  !  --- Find crystal coordinates of gtr, n1,n2,n3 ---
                  !  --- and the position ng in the ngm array      ---
                  IF (gtr(1)**2+gtr(2)**2+gtr(3)**2 <= gcutm) THEN
                     n1 = NINT( gtr(1)*at(1,1) + gtr(2)*at(2,1) + &
                                gtr(3)*at(3,1) )
                     n2 = NINT( gtr(1)*at(1,2) + gtr(2)*at(2,2) + &
                                gtr(3)*at(3,2) )
                     n3 = NINT( gtr(1)*at(1,3) + gtr(2)*at(2,3) + &
                                gtr(3)*at(3,3) )
                     ng = ln(n1,n2,n3) 
                     IF ( (ABS(g(1,ng)-gtr(1)) > eps) .OR. &
                          (ABS(g(2,ng)-gtr(2)) > eps) .OR. &
                          (ABS(g(3,ng)-gtr(3)) > eps) ) THEN
                        !
                        WRITE(6,*) ' error hepsiher: translated G=', gtr(1), gtr(2),  &
                                   gtr(3), ' with crystal coordinates', n1, n2, n3,   &
                                   ' corresponds to ng=', ng, ' but G(ng)=', g(1,ng), &
                                   g(2,ng), g(3,ng)
                        WRITE(6,*) ' probably because G_par is NOT', &
                                   ' a reciprocal lattice vector '
                        WRITE(6,*) ' Possible choices as smallest ', ' G_par:'
                        !
                        DO i = 1, 50
                           WRITE(6,*) ' i=', i,'   G=', g(1,i), g(2,i), g(3,i)
                        ENDDO
                        STOP
                        !
                     ENDIF
                  ELSE
                     WRITE(6,*) ' |gtr| > gcutm  for gtr=', gtr(1), gtr(2), gtr(3) 
                     STOP
                  ENDIF
                  map_g(ig)=ng
                  !
               ENDDO
               !
            ENDIF
            !
            ! ... OPTIMIZATION BY AM.
            ! NOTE: THERE ARE TOO MANY COMMUNICATION CALLS FOR GLOBAL ARRAY.
            !       WE CAN REDUCE THEM SIGNIFICANTLY !
            ! NOTE: CHANGED ORDER OF LOOPS OVER BANDS !
            !
            DO mb = 1, nbnd
               IF (l_para) THEN
                  ! allocate global array
                  ALLOCATE( aux_g(ngm_g) )
                  aux_g = (0.d0,0.d0)
                  IF (noncolin) THEN
                     ALLOCATE( aux_g_2(ngm_g) )
                     aux_g_2 = (0.d0,0.d0)
                  ENDIF
                  ! put psi1 on global array
                  DO ig = 1, npw0
                     aux_g(mapgp_global(ig_l2g(igk0(ig)),pdir)) = evct(ig,mb)
                     IF (noncolin) aux_g_2(mapgp_global(ig_l2g(igk0(ig)),pdir)) = evct(ig+npwx,mb)
                  ENDDO
                  CALL mp_sum( aux_g(:), intra_bgrp_comm )
                  IF (noncolin) CALL mp_sum( aux_g_2(:), intra_bgrp_comm )
               ENDIF
               !
               DO nb = 1, nbnd
                  IF (.NOT. l_cal(nb) .OR. .NOT. l_cal(mb)) THEN
                     IF (nb == mb)  mat(nb,mb) = 0.d0
                  ELSE
                     IF (.NOT. l_para) THEN
                        aux = (0.d0,0.d0)
                        aux0 = (0.d0,0.d0)
                        IF (noncolin) aux_2 = (0.d0,0.d0)
                        IF (noncolin) aux0_2 = (0.d0,0.d0)
                        DO ig = 1, npw1
                           aux0(igk1(ig)) = evcel(ig,nb)
                           IF (noncolin) aux0_2(igk1(ig)) = evcel(ig+npwx,nb)
                        ENDDO
                        !
                        DO ig = 1, npw0
                           aux(map_g(ig)) = evct(ig,mb)
                           IF (noncolin) aux_2(map_g(ig)) = evct(ig+npwx,mb)
                        ENDDO
                        !
                        mat(nb,mb) = dot_product(aux0(1:ngm),aux(1:ngm))
                        IF (noncolin) mat(nb,mb) = mat(nb,mb)+dot_product(aux0_2(1:ngm),aux_2(1:ngm))
                     ELSE
                        sca = (0.d0,0.d0)
                        ! do scalar product
                        DO ig = 1, npw1
                           sca = sca + CONJG( evcel(ig,nb))*aux_g(ig_l2g(igk1(ig)) )
                           IF (noncolin) sca = sca + CONJG( evcel(ig+npwx,nb))*aux_g_2(ig_l2g(igk1(ig)) )
                        ENDDO
                        ! mp_sum is done later
                        mat(nb,mb) = sca
                     ENDIF
                  ENDIF
               ENDDO
               !
               IF (l_para) THEN
                  DEALLOCATE( aux_g )
                  IF (noncolin) DEALLOCATE( aux_g_2 )
               ENDIF
               !
            ENDDO
            !
            CALL mp_sum( mat, intra_bgrp_comm )
            !
            DO nb = 1, nbnd
               IF(.NOT. l_cal(nb)) mat(nb,nb) = 1.d0
            ENDDO
            !
            DO nb = 1, nbnd
               DO mb = 1, nbnd
                  !
                  IF ( l_cal(nb) .AND. l_cal(mb) ) THEN
                     !   --- Calculate the augmented part: ij=KB projectors, ---
                     !   --- R=atom index: SUM_{ijR} q(ijR) <u_nk|beta_iR>   ---
                     !   --- <beta_jR|u_mk'> e^i(k-k')*R =                   ---
                     !   --- also <u_nk|beta_iR>=<psi_nk|beta_iR> = becp^*   ---
                     IF (okvan) THEN
                        !
                        pref = (0.d0,0.d0)
                        DO jkb=1,nkb
                           nhjkb = nkbtonh(jkb)
                           na = nkbtona(jkb)
                           np = ityp(na)
                           nhjkbm = nh(np)
                           jkb1 = jkb - nhjkb
                           DO j = 1,nhjkbm
                              IF (lspinorb) THEN
                                 pref = pref + CONJG(becp_bp%nc(jkb,1,nb))*becp0%nc(jkb1+j,1,mb) &
                                               * q_dkp_so(nhjkb,j,1,np)*CONJG(struc(na))
                                 pref = pref + CONJG(becp_bp%nc(jkb,1,nb))*becp0%nc(jkb1+j,2,mb) &
                                               * q_dkp_so(nhjkb,j,2,np)*CONJG(struc(na))
                                 pref = pref + CONJG(becp_bp%nc(jkb,2,nb))*becp0%nc(jkb1+j,1,mb) &
                                               * q_dkp_so(nhjkb,j,3,np)*CONJG(struc(na))
                                 pref = pref + CONJG(becp_bp%nc(jkb,2,nb))*becp0%nc(jkb1+j,2,mb) &
                                               * q_dkp_so(nhjkb,j,4,np)*CONJG(struc(na))
                              ELSE
                                 pref = pref + CONJG(becp_bp%k(jkb,nb))*becp0%k(jkb1+j,mb) &
                                               * q_dkp(nhjkb,j,np)*CONJG(struc(na))
                              ENDIF
                           ENDDO
                        ENDDO
                        mat(nb,mb) = mat(nb,mb) + pref
                        !
                     ENDIF
                     !
                  ENDIF
                  !
               ENDDO
            ENDDO
            !
            !  --- Calculate matrix inverse ---
            CALL zgefa( mat, nbnd, nbnd, ivpt, info )
            CALL errore( 'h_epsi_her_set','error in zgefa', ABS(info) )
            CALL zgedi( mat, nbnd, nbnd, ivpt, cdet, cdwork, 1 )
            matbig(nx_el(ik,pdir),:,:) = mat
            mdone(nx_el(ik,pdir)) = 1
            !
         ELSE !mdone(nx_el(ik,pdir)) /= 0
            !
            mat = matbig(nx_el(ik,pdir),:,:)
            !
         ENDIF
         !
         DO nb = 1, nbnd
            DO mb = 1, nbnd
               IF ( .NOT. l_cal(nb) .OR. .NOT. l_cal(mb) ) mat(mb,nb) = (0.d0,0.d0)
            ENDDO
         ENDDO
         !
         ! mat=S^-1(k,k-1)
         !
         IF (.NOT. l_para) THEN
            !
            DO ig = 1, npw0
               gtr(1) = g(1,igk0(ig)) + gpar(1)
               gtr(2) = g(2,igk0(ig)) + gpar(2)        
               gtr(3) = g(3,igk0(ig)) + gpar(3) 
               ! 
               ! --- Find crystal coordinates of gtr, n1,n2,n3 ---
               ! --- and the position ng in the ngm array ---
               IF (gtr(1)**2+gtr(2)**2+gtr(3)**2 <= gcutm) THEN
                  !
                  n1 = NINT( gtr(1)*at(1,1) + gtr(2)*at(2,1) + &
                             gtr(3)*at(3,1) )
                  n2 = NINT( gtr(1)*at(1,2) + gtr(2)*at(2,2) + &
                             gtr(3)*at(3,2) )
                  n3 = NINT( gtr(1)*at(1,3) + gtr(2)*at(2,3) + &
                             gtr(3)*at(3,3) )
                  ng = ln0(n1,n2,n3)
                  !
                  IF (ng > 0) THEN
                     DO m = 1, nbnd
                        DO nb = 1, nbnd
                           evcm(ng,m,pdir) = evcm(ng,m,pdir) + mat(nb,m)*evct(ig,nb)
                           IF (noncolin) evcm(ng+npwx,m,pdir) = evcm(ng+npwx,m,pdir) + &
                                                                mat(nb,m)*evct(ig+npwx,nb)
                        ENDDO
                     ENDDO
                  ENDIF
                  !
               ENDIF
            ENDDO
            !
         ELSE
            !
            ALLOCATE( aux_g(ngm_g) )
            IF (noncolin) ALLOCATE( aux_g_2(ngm_g) )
            !
            DO nb = 1, nbnd
               aux_g(:) = (0.d0,0.d0)
               IF (noncolin) aux_g_2(:) = (0.d0,0.d0)
               DO ig = 1, npw0
                  aux_g(mapgp_global(ig_l2g(igk0(ig)),pdir)) = evct(ig,nb)
                  IF (noncolin) aux_g_2(mapgp_global(ig_l2g(igk0(ig)),pdir)) = evct(ig+npwx,nb)
               ENDDO
               ! put evct on global  array
               CALL mp_sum( aux_g(:), intra_bgrp_comm )
               IF (noncolin) CALL mp_sum( aux_g_2(:), intra_bgrp_comm )
               DO m = 1, nbnd
                  DO ig = 1, npw1
                     evcm(ig,m,pdir) = evcm(ig,m,pdir) + mat(nb,m)*aux_g(ig_l2g(igk1(ig)))
                     IF (noncolin) evcm(ig+npwx,m,pdir) = evcm(ig+npwx,m,pdir) + &
                                               mat(nb,m)*aux_g_2(ig_l2g(igk1(ig)))
                  ENDDO
               ENDDO
            ENDDO
            !
            DEALLOCATE( aux_g )
            IF (noncolin) DEALLOCATE( aux_g_2 )
            !
         ENDIF
         !
         IF (okvan) THEN
            evct(:,:) =  (0.d0, 0.d0)
            ps(:,:) = (0.d0, 0.d0)
            ijkb0 = 0
            !
            DO nt = 1, ntyp
               DO na = 1, nat
                  !
                  IF (ityp(na) == nt) THEN
                     DO ibnd = 1, nbnd
                        !
                        DO jh = 1, nh(nt)
                           jkb = ijkb0 + jh
                           DO ih = 1, nh(nt)
                              !
                              ikb = ijkb0 + ih
                              IF (lspinorb) THEN
                                 ps(ikb,(ibnd-1)*npol+1) = ps(ikb,(ibnd-1)*npol+1 ) + &
                                      q_dkp_so(ih,jh,1,ityp(na))*CONJG(struc(na))* becp0%nc(jkb,1,ibnd)
                                 ps(ikb,(ibnd-1)*npol+1) = ps(ikb,(ibnd-1)*npol+1 ) + &
                                      q_dkp_so(ih,jh,2,ityp(na))*CONJG(struc(na))* becp0%nc(jkb,2,ibnd)
                                 ps(ikb,(ibnd-1)*npol+2) = ps(ikb,(ibnd-1)*npol+2 ) + &
                                      q_dkp_so(ih,jh,3,ityp(na))*CONJG(struc(na))* becp0%nc(jkb,1,ibnd)
                                 ps(ikb,(ibnd-1)*npol+2) = ps(ikb,(ibnd-1)*npol+2 ) + &
                                      q_dkp_so(ih,jh,4,ityp(na))*CONJG(struc(na))* becp0%nc(jkb,2,ibnd)
                              ELSE
                                 ps(ikb,ibnd) = ps(ikb,ibnd) + &
                                      q_dkp(ih,jh,ityp(na))*CONJG(struc(na))* becp0%k(jkb,ibnd)
                              ENDIF
                              !
                           ENDDO
                        ENDDO
                        !
                     ENDDO
                     ijkb0 = ijkb0 + nh(nt)
                  ENDIF
                  !
               ENDDO
            ENDDO
            !
            CALL ZGEMM( 'N', 'N', npw1, nbnd*npol , nkb, (1.d0, 0.d0) , vkb, & ! vkb is relative to
                        npwx, ps, nkb, (1.d0, 0.d0) , evct, npwx )             ! the last ik read.
            !
            DO m = 1, nbnd
               DO nb = 1, nbnd
                  DO ig = 1, npw1
                     evcm(ig,m,pdir) = evcm(ig,m,pdir) + mat(nb,m)*evct(ig,nb)
                     IF (noncolin) THEN
                        evcm(ig+npwx,m,pdir) = evcm(ig+npwx,m,pdir) + mat(nb,m)*evct(ig+npwx,nb)
                     ENDIF
                  ENDDO
               ENDDO
            ENDDO
            !
         ENDIF ! okvan
         !
      ENDIF ! ik_stringa
      !
      ! calculate  S-1(k,k+1)
      ! 
      IF (ik_stringa /= nppstr_3d(pdir)) THEN
         !
         ikk = nx_el(ik+1,pdir)
         npw0   = ngk(ikk)
         igk0(:) = igk_k(:,ikk)
         !
         CALL get_buffer( evct, nwordwfc, iunwfc, nx_el(ik+1,pdir) )
         !        
         ! --- Calculate dot products between wavefunctions
         !
         ! --- Dot wavefunctions and betas for PREVIOUS k-point ---
         !
         IF (okvan) THEN
            CALL init_us_2( npw0, igk0, xk(1,nx_el(ik+1,pdir)), vkb )
            CALL calbec( npw0, vkb, evct, becp0 )
         ENDIF
         ! --- Dot wavefunctions and betas for CURRENT k-point ---
         !
         ikk = nx_el(ik,pdir)
         npw1 = ngk(ikk)
         igk1(:) = igk_k(:,ikk)
         !
         !  --- Recalculate FFT correspondence (see ggen.f90) ---
         !
         ln0 = 0 !set to  0
         DO ig = 1, npw1
            mk1 = NINT( g(1,igk1(ig))*at(1,1) + g(2,igk1(ig))*at(2,1) + g(3,igk1(ig))*at(3,1) )
            mk2 = NINT( g(1,igk1(ig))*at(1,2) + g(2,igk1(ig))*at(2,2) + g(3,igk1(ig))*at(3,2) )
            mk3 = NINT( g(1,igk1(ig))*at(1,3) + g(2,igk1(ig))*at(2,3) + g(3,igk1(ig))*at(3,3) )
            ln0(mk1,mk2,mk3) = ig
         ENDDO
         !
         IF (okvan) THEN
            CALL init_us_2( npw1, igk1, xk(1,nx_el(ik,pdir)), vkb )
            CALL calbec( npw1, vkb, evcel, becp_bp )
         ENDIF
         !
         !  --- Matrix elements calculation ---
         !
         IF (mdone(nx_el(ik+1,pdir)) == 0) THEN
            !
            mat = (0.d0,0.d0)
            ALLOCATE( aux0vec(ngm,nbnd), aux1vec(ngm,nbnd) )
            aux0vec = (0.d0,0.d0)
            aux1vec = (0.d0,0.d0)
            !
            DO nb = 1, nbnd
               DO ig = 1, npw1
                  aux0vec(igk1(ig),nb) = evcel(ig,nb)
               ENDDO
            ENDDO
            !
            DO nb=1,nbnd
               DO ig=1,npw0
                  aux1vec(igk0(ig),nb)=evct(ig,nb)
               ENDDO
            ENDDO
            !
            CALL ZGEMM( 'C', 'N', nbnd, nbnd, ngm, (1.d0,0.d0), aux0vec, ngm, &
                        aux1vec, ngm, (0.d0,0.d0), mat, nbnd )
            !
            IF (noncolin) THEN
               aux0vec = (0.d0,0.d0)
               aux1vec = (0.d0,0.d0)
               DO nb = 1, nbnd
                  DO ig = 1, npw1
                     aux0vec(igk1(ig),nb) = evcel(ig+npwx,nb)
                  ENDDO
               ENDDO
               DO nb = 1, nbnd
                  DO ig = 1, npw0
                     aux1vec(igk0(ig),nb) = evct(ig+npwx,nb)
                  ENDDO
               ENDDO
               CALL ZGEMM( 'C', 'N', nbnd, nbnd, ngm, (1.d0,0.d0), aux0vec, ngm, &
                           aux1vec, ngm, (1.d0,0.d0), mat, nbnd )
            ENDIF
            !
            DEALLOCATE( aux0vec, aux1vec )
            !
            CALL mp_sum( mat, intra_bgrp_comm )
            !
            DO nb=1,nbnd
               DO mb=1,nbnd
                  IF (.NOT. l_cal(nb) .OR. .NOT.l_cal(mb)) THEN
                     IF (nb == mb) THEN
                        mat(nb,mb) = (1.d0,0.d0)
                     ELSE
                        mat(nb,mb) = (0.d0,0.d0)
                     ENDIF
                  ENDIF
               ENDDO
            ENDDO
            !
            DO nb = 1, nbnd
               DO mb = 1, nbnd
                  !
                  IF ( l_cal(nb) .AND. l_cal(mb) ) THEN
                     IF (okvan) THEN
                        pref = (0.d0,0.d0)
                        !
                        DO jkb = 1, nkb
                           nhjkb = nkbtonh(jkb)
                           na = nkbtona(jkb)
                           np = ityp(na)
                           nhjkbm = nh(np)
                           jkb1 = jkb - nhjkb
                           DO j = 1,nhjkbm
                              !
                              IF (lspinorb) THEN
                                 pref = pref + CONJG(becp_bp%nc(jkb,1,nb))*becp0%nc(jkb1+j,1,mb) &
                                               * q_dk_so(nhjkb,j,1,np)*struc(na)
                                 pref = pref + CONJG(becp_bp%nc(jkb,1,nb))*becp0%nc(jkb1+j,2,mb) &
                                               * q_dk_so(nhjkb,j,2,np)*struc(na)
                                 pref = pref + CONJG(becp_bp%nc(jkb,2,nb))*becp0%nc(jkb1+j,1,mb) &
                                               * q_dk_so(nhjkb,j,3,np)*struc(na)
                                 pref = pref + CONJG(becp_bp%nc(jkb,2,nb))*becp0%nc(jkb1+j,2,mb) &
                                               * q_dk_so(nhjkb,j,4,np)*struc(na)
                              ELSE
                                 pref = pref + CONJG(becp_bp%k(jkb,nb))*becp0%k(jkb1+j,mb) &
                                               * q_dk(nhjkb,j,np)*struc(na)
                              ENDIF
                              !
                           ENDDO
                        ENDDO
                        !
                        mat(nb,mb) = mat(nb,mb) + pref
                     ENDIF
                  ENDIF
                  !
               ENDDO
            ENDDO
            !
            !  --- Calculate matrix inverse ---
            !
            CALL zgefa( mat, nbnd, nbnd, ivpt, info )
            CALL errore( 'h_epsi_her_set', 'error in zgefa', ABS(info) )
            CALL zgedi( mat, nbnd, nbnd, ivpt, cdet, cdwork, 1 )
            matbig(nx_el(ik+1,pdir),:,:) = TRANSPOSE(CONJG(mat))
            mdone(nx_el(ik+1,pdir)) = 1
         ELSE
            mat = TRANSPOSE(CONJG(matbig(nx_el(ik+1,pdir),:,:)))
         ENDIF
         !
         ! in case of occupations from input set to 0 parts of mat
         !
         DO nb = 1, nbnd
            DO mb = 1, nbnd
               IF (.NOT. l_cal(nb) .OR. .NOT. l_cal(mb)) mat(mb,nb) = (0.d0,0.d0)
            ENDDO
         ENDDO
         !
         ! mat=S^-1(k,k-1)
         !
         DO ig = 1, npw0
            !
            gtr(1) = g(1,igk0(ig))
            gtr(2) = g(2,igk0(ig))         
            gtr(3) = g(3,igk0(ig))
            !
            ! --- Find crystal coordinates of gtr, n1,n2,n3 ---
            ! --- and the position ng in the ngm array      ---
            !
            IF (gtr(1)**2 + gtr(2)**2 + gtr(3)**2 <= gcutm) THEN
               n1 = NINT( gtr(1)*at(1,1) + gtr(2)*at(2,1) + &
                          gtr(3)*at(3,1) )
               n2 = NINT( gtr(1)*at(1,2) + gtr(2)*at(2,2) + &
                          gtr(3)*at(3,2) )
               n3 = NINT( gtr(1)*at(1,3) + gtr(2)*at(2,3) + &
                          gtr(3)*at(3,3))
               ng = ln0(n1,n2,n3)
               !
               IF (ng > 0) THEN
                  DO m = 1, nbnd
                     DO nb = 1, nbnd
                        evcp(ng,m,pdir) = evcp(ng,m,pdir) + mat(nb,m)*evct(ig,nb)
                        IF (noncolin) evcp(ng+npwx,m,pdir) = evcp(ng+npwx,m,pdir) + &
                                                             mat(nb,m)*evct(ig+npwx,nb)
                     ENDDO
                  ENDDO
               ENDIF
               !
            ENDIF
            !
         ENDDO
         !
         IF (okvan) THEN
            evct(:,:) = (0.d0, 0.d0)
            ps(:,:) = (0.d0, 0.d0)
            ijkb0 = 0
            DO nt = 1, ntyp
               DO na = 1, nat
                  IF (ityp(na) == nt) THEN
                     DO ibnd = 1, nbnd
                        !
                        DO jh = 1, nh(nt)
                           jkb = ijkb0 + jh
                           DO ih = 1, nh(nt)
                              ikb = ijkb0 + ih
                              !
                              IF (noncolin) THEN
                                 ps(ikb,(ibnd-1)*npol+1) = ps(ikb,(ibnd-1)*npol+1 ) + &
                                      q_dk_so(ih,jh,1,ityp(na))*struc(na)* becp0%nc(jkb,1,ibnd)
                                 ps(ikb,(ibnd-1)*npol+1) = ps(ikb,(ibnd-1)*npol+1 ) + &
                                      q_dk_so(ih,jh,2,ityp(na))*struc(na)* becp0%nc(jkb,2,ibnd)
                                 ps(ikb,(ibnd-1)*npol+2) = ps(ikb,(ibnd-1)*npol+2 ) + &
                                      q_dk_so(ih,jh,3,ityp(na))*struc(na)* becp0%nc(jkb,1,ibnd)
                                 ps(ikb,(ibnd-1)*npol+2) = ps(ikb,(ibnd-1)*npol+2 ) + &
                                      q_dk_so(ih,jh,4,ityp(na))*struc(na)* becp0%nc(jkb,2,ibnd)
                              ELSE
                                 ps(ikb,ibnd) = ps(ikb,ibnd) + &
                                      q_dk(ih,jh,ityp(na))*struc(na)* becp0%k(jkb,ibnd)
                              ENDIF
                              !
                           ENDDO
                        ENDDO
                        !
                     ENDDO
                     ijkb0 = ijkb0 + nh(nt)
                  ENDIF
               ENDDO
            ENDDO
            !
            ! vkb is relative to the last ik read
            CALL ZGEMM( 'N', 'N', npw1, nbnd*npol, nkb, (1.d0, 0.d0), vkb, &
                        npwx, ps, nkb, (1.d0, 0.d0), evct, npwx )
            !        
            evct(npw1+1:npwx,1:nbnd) = (0.d0,0.d0)
            IF (noncolin) evct(npwx+npw1+1:2*npwx,1:nbnd) = (0.d0,0.d0)
            evcp(npw1+1:npwx,1:nbnd,pdir) = (0.d0,0.d0)
            IF (noncolin) evcp(npwx+npw1+1:2*npwx,1:nbnd,pdir) = (0.d0,0.d0)
            IF (.NOT.noncolin) THEN
               CALL ZGEMM( 'N', 'N', npw1, nbnd, nbnd, (1.d0,0.d0), evct, npwx, mat, nbnd, &
                           (1.d0,0.d0), evcp(1,1,pdir), npwx )
            ELSE
               CALL ZGEMM( 'N', 'N', npwx*npol, nbnd, nbnd, (1.d0,0.d0), evct, npwx*npol,  &
                           mat, nbnd, (1.d0,0.d0), evcp(1,1,pdir), npwx*npol )
            ENDIF
            !
         ENDIF
         ! --- End of dot products between wavefunctions and betas ---
      ELSE
         !      
         ikk = nx_el(ik-nppstr_3d(pdir)+1,pdir)
         npw0 = ngk(ikk)
         igk0(:)= igk_k(:,ikk)
         !
         CALL get_buffer( evct, nwordwfc, iunwfc, nx_el(ik-nppstr_3d(pdir)+1,pdir) )
         !        
         ! --- Calculate dot products between wavefunctions
         !
         ! --- Dot wavefunctions and betas for PREVIOUS k-point ---
         !
         IF (okvan) THEN
            CALL init_us_2( npw0, igk0, xk(1,nx_el(ik-nppstr_3d(pdir)+1,pdir)), vkb )
            CALL calbec( npw0, vkb, evct, becp0 )
         ENDIF
         ! --- Dot wavefunctions and betas for CURRENT k-point ---
         !
         ikk = nx_el(ik,pdir)
         npw1 = ngk(ikk)
         igk1(:) = igk_k(:,ikk)
         !
         !  --- Recalculate FFT correspondence (see ggen.f90) ---
         !
         IF (.NOT. l_para) THEN
            ln0 = 0  ! set to 0
            DO ig = 1, npw1
               mk1 = NINT( g(1,igk1(ig))*at(1,1)+g(2,igk1(ig))*at(2,1)+g(3,igk1(ig))*at(3,1) )
               mk2 = NINT( g(1,igk1(ig))*at(1,2)+g(2,igk1(ig))*at(2,2)+g(3,igk1(ig))*at(3,2) )
               mk3 = NINT( g(1,igk1(ig))*at(1,3)+g(2,igk1(ig))*at(2,3)+g(3,igk1(ig))*at(3,3) )
               ln0(mk1,mk2,mk3) = ig
            ENDDO
         ENDIF
         !
         IF (okvan) THEN
            CALL init_us_2( npw1, igk1, xk(1,nx_el(ik,pdir)), vkb )
            CALL calbec( npw1, vkb, evcel, becp_bp )
         ENDIF
         ! --- Matrix elements calculation ---
         !
         IF (mdone(nx_el(ik-nppstr_3d(pdir)+1,pdir)) == 0) THEN
            !
            IF (.NOT. l_para) THEN
               !
               map_g(:) = 0
               !
               DO ig = 1, npw0
                  !  --- If k'=k+G_o, the relation psi_k+G_o (G-G_o) ---
                  !  --- = psi_k(G) is used, gpar=G_o, gtr = G-G_o   ---
                  gtr(1) = g(1,igk0(ig)) - gpar(1)
                  gtr(2) = g(2,igk0(ig)) - gpar(2) 
                  gtr(3) = g(3,igk0(ig)) - gpar(3) 
                  ! --- Find crystal coordinates of gtr, n1,n2,n3 ---
                  ! --- and the position ng in the ngm array ---
                  IF (gtr(1)**2+gtr(2)**2+gtr(3)**2 <= gcutm) THEN
                     n1 = NINT( gtr(1)*at(1,1)+gtr(2)*at(2,1) + &
                                gtr(3)*at(3,1) )
                     n2 = NINT( gtr(1)*at(1,2)+gtr(2)*at(2,2) + &
                                gtr(3)*at(3,2) )
                     n3 = NINT( gtr(1)*at(1,3)+gtr(2)*at(2,3) + &
                                gtr(3)*at(3,3) )
                     ng = ln(n1,n2,n3) 
                     IF ((ABS(g(1,ng)-gtr(1)) > eps) .OR. &
                          (ABS(g(2,ng)-gtr(2)) > eps) .OR. &
                          (ABS(g(3,ng)-gtr(3)) > eps)) THEN
                        WRITE(6,*) ' error hepsiher: translated G=', &
                             gtr(1),gtr(2),gtr(3),                   &
                             ' with crystal coordinates',n1,n2,n3,   &
                             ' corresponds to ng=',ng,' but G(ng)=', &
                             g(1,ng),g(2,ng),g(3,ng)
                        WRITE(6,*) ' probably because G_par is NOT', &
                             ' a reciprocal lattice vector '
                        WRITE(6,*) ' Possible choices as smallest ', &
                             ' G_par:'
                        DO i = 1, 50
                           WRITE(6,*) ' i=',i,'   G=', g(1,i),g(2,i),g(3,i)
                        ENDDO
                        STOP
                        !
                     ENDIF
                  ELSE 
                     WRITE(6,*) ' |gtr| > gcutm  for gtr=',gtr(1),gtr(2),gtr(3) 
                     STOP
                  ENDIF
                  map_g(ig)=ng 
               ENDDO
               !
            ENDIF
            !
            mat=(0.d0,0.d0)
            !
            DO mb = 1, nbnd
               !
               IF (l_para) THEN
                  !ALLOCATE global array
                  ALLOCATE( aux_g(ngm_g) )
                  aux_g = (0.d0,0.d0)
                  IF (noncolin) THEN
                    ALLOCATE( aux_g_2(ngm_g) )
                    aux_g_2 = (0.d0,0.d0)
                  ENDIF
                  ! put psi1 on global array
                  DO ig = 1, npw0
                     aux_g(mapgm_global(ig_l2g(igk0(ig)),pdir)) = evct(ig,mb)
                     IF (noncolin) aux_g_2(mapgm_global(ig_l2g(igk0(ig)),pdir)) = evct(ig+npwx,mb)
                  ENDDO
                  CALL mp_sum( aux_g(:), intra_bgrp_comm )
                  IF (noncolin) CALL mp_sum( aux_g_2(:), intra_bgrp_comm )
               ENDIF
               !
               DO nb = 1, nbnd
                  IF ( .NOT. l_cal(nb) .OR. .NOT. l_cal(mb) ) THEN
                     IF ( nb == mb )  mat(nb,mb) = 0.d0
                  ELSE
                     !
                     IF (.NOT. l_para) THEN
                        aux = (0.d0,0.d0)
                        aux0 = (0.d0,0.d0)
                        IF (noncolin) aux_2 = (0.d0,0.d0)
                        IF (noncolin) aux0_2 = (0.d0,0.d0)
                        DO ig = 1, npw1
                           aux0(igk1(ig)) = evcel(ig,nb)
                           IF(noncolin) aux0_2(igk1(ig)) = evcel(ig+npwx,nb)
                        ENDDO
                        DO ig = 1, npw0
                           aux(map_g(ig)) = evct(ig,mb)
                           IF(noncolin) aux_2(map_g(ig)) = evct(ig+npwx,mb)
                        ENDDO
                        mat(nb,mb) = dot_product(aux0(1:ngm),aux(1:ngm))
                        IF (noncolin) mat(nb,mb) = mat(nb,mb) + dot_product(aux0_2(1:ngm),aux_2(1:ngm))
                     ELSE
                        sca = (0.d0,0.d0)
                        ! do scalar product
                        DO ig = 1, npw1
                           sca = sca + CONJG(evcel(ig,nb))*aux_g(ig_l2g(igk1(ig)))
                           IF (noncolin) sca = sca + CONJG( evcel(ig+npwx,nb)) * &
                                                     aux_g_2(ig_l2g(igk1(ig)) )
                        ENDDO
                        ! mp_sum is done later
                        mat(nb,mb) = sca
                        !
                     ENDIF
                  ENDIF
               ENDDO
               !
               IF (l_para) THEN
                  DEALLOCATE( aux_g )
                  IF (noncolin) DEALLOCATE( aux_g_2 )
               ENDIF
               !
            ENDDO
            ! 
            CALL mp_sum( mat, intra_bgrp_comm )
            DO nb = 1, nbnd
               IF (.NOT. l_cal(nb)) mat(nb,nb) = 1.d0
            ENDDO
            !
            DO nb = 1, nbnd
               DO mb = 1, nbnd
                  !
                  IF ( l_cal(nb) .AND. l_cal(mb) ) THEN
                     IF (okvan) THEN
                        pref = (0.d0,0.d0)
                        DO jkb = 1, nkb
                           nhjkb = nkbtonh(jkb)
                           na = nkbtona(jkb)
                           np = ityp(na)
                           nhjkbm = nh(np)
                           jkb1 = jkb - nhjkb
                           DO j = 1,nhjkbm
                              IF (lspinorb) THEN
                                 pref = pref + CONJG(becp_bp%nc(jkb,1,nb))*becp0%nc(jkb1+j,1,mb) &
                                               *q_dk_so(nhjkb,j,1,np)*struc(na)
                                 pref = pref + CONJG(becp_bp%nc(jkb,1,nb))*becp0%nc(jkb1+j,2,mb) &
                                               *q_dk_so(nhjkb,j,2,np)*struc(na)
                                 pref = pref + CONJG(becp_bp%nc(jkb,2,nb))*becp0%nc(jkb1+j,1,mb) &
                                               *q_dk_so(nhjkb,j,3,np)*struc(na)
                                 pref = pref + CONJG(becp_bp%nc(jkb,2,nb))*becp0%nc(jkb1+j,2,mb) &
                                               *q_dk_so(nhjkb,j,4,np)*struc(na)
                              ELSE
                                 pref = pref + CONJG(becp_bp%k(jkb,nb))*becp0%k(jkb1+j,mb) &
                                               *q_dk(nhjkb,j,np)*struc(na)
                              ENDIF
                           ENDDO
                        ENDDO
                        mat(nb,mb) = mat(nb,mb) + pref
                     ENDIF
                  ENDIF
                  !
               ENDDO
            ENDDO
            ! --- Calculate matrix inverse ---
            CALL zgefa( mat, nbnd, nbnd, ivpt, info )
            CALL errore( 'h_epsi_her_set','error in zgefa', ABS(info) )
            CALL zgedi( mat, nbnd, nbnd, ivpt, cdet, cdwork, 1 )
            matbig(nx_el(ik-nppstr_3d(pdir)+1,pdir),:,:) = TRANSPOSE( CONJG(mat) )
            mdone(nx_el(ik-nppstr_3d(pdir)+1,pdir)) = 1
            !
         ELSE
            !
            mat = TRANSPOSE( CONJG(matbig(nx_el(ik-nppstr_3d(pdir)+1,pdir),:,:)) )
            !
         ENDIF
         !
         !
         DO nb = 1, nbnd
            DO mb = 1, nbnd
               IF (.NOT. l_cal(nb) .OR. .NOT. l_cal(mb)) mat(mb,nb) = (0.d0,0.d0)
            ENDDO
         ENDDO
         !
         !
         ! mat=S^-1(k,k-1)
         IF (.NOT. l_para) THEN
            DO ig = 1, npw0
               gtr(1) = g(1,igk0(ig)) - gpar(1)
               gtr(2) = g(2,igk0(ig)) - gpar(2)        
               gtr(3) = g(3,igk0(ig)) - gpar(3)        
               ! --- Find crystal coordinates of gtr, n1,n2,n3 ---
               ! --- and the position ng in the ngm array      ---
               IF (gtr(1)**2 + gtr(2)**2 + gtr(3)**2 <= gcutm) THEN
                  n1 = NINT( gtr(1)*at(1,1)+gtr(2)*at(2,1) + &
                             gtr(3)*at(3,1) )
                  n2 = NINT( gtr(1)*at(1,2)+gtr(2)*at(2,2) + &
                             gtr(3)*at(3,2))
                  n3 = NINT( gtr(1)*at(1,3)+gtr(2)*at(2,3) + &
                             gtr(3)*at(3,3))
                  ng = ln0(n1,n2,n3)
                  !
                  IF (ng > 0) THEN
                     DO m = 1, nbnd
                        DO nb = 1, nbnd
                           evcp(ng,m,pdir) = evcp(ng,m,pdir) + mat(nb,m)*evct(ig,nb)
                           IF (noncolin) evcp(ng+npwx,m,pdir) = evcp(ng+npwx,m,pdir) &
                                                         + mat(nb,m)*evct(ig+npwx,nb)
                        ENDDO
                     ENDDO
                  ENDIF
               ENDIF
            ENDDO
            !
         ELSE
            !
            ALLOCATE( aux_g(ngm_g) )
            IF (noncolin) ALLOCATE( aux_g_2(ngm_g) )
            !
            DO nb = 1, nbnd
               aux_g(:)=(0.d0,0.d0)
               IF (noncolin) aux_g_2(:)=(0.d0,0.d0)
               DO ig=1,npw0
                  aux_g(mapgm_global(ig_l2g(igk0(ig)),pdir)) = evct(ig,nb)
                  IF (noncolin) aux_g_2(mapgm_global(ig_l2g(igk0(ig)),pdir)) = evct(ig+npwx,nb)
               ENDDO
               ! put evct on global array
               CALL mp_sum( aux_g(:), intra_bgrp_comm )
               IF (noncolin) CALL mp_sum( aux_g_2(:), intra_bgrp_comm )
               !
               DO m = 1, nbnd
                  DO ig = 1, npw1
                     evcp(ig,m,pdir) = evcp(ig,m,pdir) + mat(nb,m)*aux_g(ig_l2g(igk1(ig)))
                     IF (noncolin) evcp(ig+npwx,m,pdir) = evcp(ig+npwx,m,pdir) &
                                             + mat(nb,m)*aux_g_2(ig_l2g(igk1(ig)))
                  ENDDO
               ENDDO
               !
            ENDDO
            !
            DEALLOCATE( aux_g )
            IF (noncolin) DEALLOCATE( aux_g_2 )
            !
         ENDIF
         !
         IF (okvan) THEN
            evct(:,:) =  (0.d0, 0.d0)
            ps(:,:) = (0.d0, 0.d0)
            ijkb0 = 0
            DO nt = 1, ntyp
               DO na = 1, nat
                  IF (ityp(na) == nt) THEN
                     DO ibnd = 1, nbnd
                        DO jh = 1, nh (nt)
                           jkb = ijkb0 + jh
                           DO ih = 1, nh (nt)
                              ikb = ijkb0 + ih
                              IF (lspinorb) THEN
                                 ps(ikb, (ibnd-1)*npol+1) = ps(ikb,(ibnd-1)*npol+1 ) + &
                                      q_dk_so(ih,jh,1,ityp(na))*struc(na)* becp0%nc(jkb,1,ibnd)
                                 ps(ikb, (ibnd-1)*npol+1) = ps(ikb,(ibnd-1)*npol+1 ) + &
                                      q_dk_so(ih,jh,2,ityp(na))*struc(na)* becp0%nc(jkb,2,ibnd)
                                 ps(ikb, (ibnd-1)*npol+2) = ps(ikb,(ibnd-1)*npol+2 ) + &
                                      q_dk_so(ih,jh,3,ityp(na))*struc(na)* becp0%nc(jkb,1,ibnd)
                                 ps(ikb, (ibnd-1)*npol+2) = ps(ikb,(ibnd-1)*npol+2 ) + &
                                      q_dk_so(ih,jh,4,ityp(na))*struc(na)* becp0%nc(jkb,2,ibnd)
                              ELSE
                                 ps(ikb, ibnd) = ps(ikb, ibnd) + &
                                      q_dk(ih,jh,ityp(na))*struc(na)* becp0%k(jkb,ibnd)
                              ENDIF
                           ENDDO
                        ENDDO
                     ENDDO
                     ijkb0 = ijkb0 + nh (nt)
                  ENDIF
               ENDDO
            ENDDO
            !
            !vkb is relative to the ik read
            CALL ZGEMM( 'N', 'N', npw1, nbnd*npol , nkb, (1.d0, 0.d0), vkb, &
                        npwx, ps, nkb, (1.d0, 0.d0), evct, npwx )
            !
            DO m = 1, nbnd
               DO nb = 1, nbnd
                  DO ig = 1, npw1
                     evcp(ig,m,pdir) = evcp(ig,m,pdir) + mat(nb,m)*evct(ig,nb)
                     IF (noncolin) THEN
                        evcp(ig+npwx,m,pdir) = evcp(ig+npwx,m,pdir) + mat(nb,m)*evct(ig+npwx,nb)
                     ENDIF
                  ENDDO
               ENDDO
            ENDDO
            !
         ENDIF
         !
      ENDIF
      !writes projectors to disk 
      CALL save_buffer( evcm(:,:,pdir), nwordwfc, iunefieldm, nx_el(ik,pdir)+(pdir-1)*nks )
      CALL save_buffer( evcp(:,:,pdir), nwordwfc, iunefieldp, nx_el(ik,pdir)+(pdir-1)*nks )
      !
   ENDDO !on ik
   !
   !
   DEALLOCATE( l_cal )
   DEALLOCATE( evct  )
   DEALLOCATE( map_g )
   DEALLOCATE( ln, ln0 )
   DEALLOCATE( aux, aux0 )
   IF (ALLOCATED( aux_2  )) DEALLOCATE( aux_2  )
   IF (ALLOCATED( aux0_2 )) DEALLOCATE( aux0_2 )
   IF (okvan) CALL deallocate_bec_type( becp0   )
   IF (okvan) CALL deallocate_bec_type( becp_bp )
   IF (okvan .AND. lspinorb) DEALLOCATE( q_dk_so  )
   IF (okvan .AND. lspinorb) DEALLOCATE( q_dkp_so )
   DEALLOCATE( mat )
   !
   CALL stop_clock( 'h_epsi_set' )
   !
   !
   RETURN
   !
END SUBROUTINE h_epsi_her_set
!
!
!--------------------------------------------------------------------
SUBROUTINE factor_a( dir, a, fact )
   !----------------------------------------------------------------
   !
   USE kinds, ONLY: DP
   !
   IMPLICIT NONE
   !
   REAL(kind=DP):: a(3,3), fact
   INTEGER :: dir
   !
   INTEGER :: d1, d2
   REAL(kind=DP) :: v(3), sca
   !
   IF (dir==1) THEN
      d1=2
      d2=3
   ELSEIF (dir==2) THEN
      d1=3
      d2=1
   ELSEIF (dir==3) THEN
      d1=1
      d2=2
   ENDIF
   !
   !calculate vect(a(d1,:) X a(d2,:)
   !
   v(1) =  a(2,d1)*a(3,d2) - a(3,d1)*a(2,d2)
   v(2) = -a(1,d1)*a(3,d2) + a(3,d1)*a(1,d2)
   v(3) =  a(1,d1)*a(2,d2) - a(2,d1)*a(1,d2)
   !
   ! normalize v
   !
   sca = SQRT( v(1)**2.d0 + v(2)**2.d0 + v(3)**2.d0 )
   v(:) = v(:) / sca
   !
   !calculate a(dir:)*v(:)
   !
   fact = v(1)*a(1,dir) + v(2)*a(2,dir) + v(3)*a(3,dir)
   !
   ! ---
   !
   fact = DSQRT( a(1,dir)**2.d0 + a(2,dir)**2.d0 + a(3,dir)**2.d0 )
   !
   fact = ABS(fact)
   !
   !
   RETURN
   !
 END SUBROUTINE factor_a
