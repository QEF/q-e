!
! Copyright (C) 2001-2004 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!============================================================================!
SUBROUTINE c_phase_field( el_pola, ion_pola, fact_pola, pdir )
   !-------------------------------------------------------------------------!
   !! Geometric phase calculation along a strip of \(\text{nppstr_3d(pdir)}\)
   !! k-points averaged over a 2D grid of \(\text{nkort}\) k-points orthogonal
   !! to \(\text{nppstr_3d(pdir)}\).
   !
   !! This routine is used to calculate the electronic polarization when a
   !! finite electric field, described through the modern theory of the 
   !! polarization, is applied.  
   !! It is very similar to the routine \(\texttt{c_phase}\) in 
   !! \(\texttt{bp_c_phase}\), however, the numbering of the k-points in the
   !! strings is different.
   !
   USE kinds,                ONLY: DP
   USE io_global,            ONLY: stdout, ionode, ionode_id
   USE io_files,             ONLY: iunwfc, nwordwfc,prefix,tmp_dir
   USE buffers,              ONLY: get_buffer
   USE ions_base,            ONLY: nat, ntyp => nsp, ityp, tau, zv, atm
   USE cell_base,            ONLY: at, alat, tpiba, omega
   USE constants,            ONLY: pi, tpi
   USE fft_base,             ONLY: dfftp
   USE gvect,                ONLY: ngm, g, gcutm, ngm_g
   USE uspp,                 ONLY: nkb, vkb, okvan
   USE uspp_param,           ONLY: upf, lmaxq, nbetam, nh, nhm
   USE upf_spinorb,          ONLY: transform_qq_so
   USE lsda_mod,             ONLY: nspin
   USE klist,                ONLY: nelec, degauss, nks, xk, wk, ngk, igk_k
   USE wvfct,                ONLY: npwx, nbnd
   USE noncollin_module,     ONLY: noncolin, npol, lspinorb
   USE bp,                   ONLY: nppstr_3d, mapgm_global, nx_el,phase_control
   USE fixed_occ
   USE gvect,                ONLY: ig_l2g
   USE mp,                   ONLY: mp_sum, mp_bcast
   USE mp_bands,             ONLY: intra_bgrp_comm
   USE mp_pools,             ONLY: intra_pool_comm
   USE becmod,               ONLY: calbec,bec_type,allocate_bec_type,deallocate_bec_type
   USE uspp_init,            ONLY : init_us_2
   IMPLICIT NONE
   !
   REAL(DP), INTENT(OUT) :: el_pola
   !! electronic polarization
   REAL(DP), INTENT(OUT) :: ion_pola
   !! ionic polarization
   REAL(DP), INTENT(OUT) :: fact_pola
   !! the prefactor of the polarization
   INTEGER, INTENT(IN) :: pdir
   !! direction on which the polarization is calculated
   !
   !  ... local variables
   !
   INTEGER :: i, ik
   INTEGER :: igk1(npwx)
   INTEGER :: igk0(npwx)
   INTEGER :: ig
   INTEGER :: info
   INTEGER :: is
   INTEGER :: istring
   INTEGER :: iv
   INTEGER :: ivpt(nbnd)
   INTEGER :: j
   INTEGER :: jkb
   INTEGER :: jkb_bp
   INTEGER :: jkb1
   INTEGER :: jv
   INTEGER :: kort
   INTEGER :: kpar
   INTEGER :: kpoint
   INTEGER :: kstart
   INTEGER :: mb
   INTEGER :: mk1
   INTEGER :: mk2
   INTEGER :: mk3
   INTEGER , ALLOCATABLE :: mod_elec(:)
   INTEGER , ALLOCATABLE :: ln(:,:,:)
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
   INTEGER :: nkort
   INTEGER :: np
   INTEGER :: npw1
   INTEGER :: npw0
   INTEGER :: nstring
   INTEGER :: nt
   INTEGER :: nspinnc
   REAL(DP) :: dk(3)
   REAL(DP) :: dk2
   REAL(DP) :: dkmod
   REAL(DP) :: el_loc
   REAL(DP) :: eps
   REAL(DP) :: fac
   REAL(DP) :: gpar(3)
   REAL(DP) :: gtr(3)
   REAL(DP) :: gvec
   REAL(DP), ALLOCATABLE :: loc_k(:)
   REAL(DP), ALLOCATABLE :: pdl_elec(:)
   REAL(DP), ALLOCATABLE :: phik(:)
   REAL(DP) :: weight
   REAL(DP) :: pola, pola_ion
   REAL(DP), ALLOCATABLE :: wstring(:)
   REAL(DP) :: zeta_mod
   COMPLEX(DP), ALLOCATABLE :: aux(:,:)
   COMPLEX(DP), ALLOCATABLE :: aux0(:,:)
   ! For noncollinear calculations
   COMPLEX(DP), ALLOCATABLE :: aux_2(:,:)
   COMPLEX(DP), ALLOCATABLE :: aux0_2(:,:)
   COMPLEX(DP), ALLOCATABLE :: cphik(:)
   COMPLEX(DP) :: det
   COMPLEX(DP) :: mat(nbnd,nbnd)
   COMPLEX(DP) :: pref
   COMPLEX(DP) :: q_dk(nhm,nhm,ntyp)
   COMPLEX(DP) :: struc(nat)
   COMPLEX(DP) :: zeta
   !
   COMPLEX(DP), ALLOCATABLE :: psi(:,:)
   COMPLEX(DP), ALLOCATABLE :: psi1(:,:)
   COMPLEX(DP) :: zeta_loc
   !
   LOGICAL, ALLOCATABLE :: l_cal(:)       ! flag for occupied/empty states
   INTEGER, ALLOCATABLE :: map_g(:)
   !
   REAL(DP) :: dkfact
   COMPLEX(DP) :: zeta_tot
   !
   LOGICAL :: l_para                      ! if true new parallel treatment
   COMPLEX(DP) :: sca
   COMPLEX(DP), ALLOCATABLE :: aux_g(:)
   COMPLEX(DP), ALLOCATABLE :: aux_g_2(:) ! noncollinear case
   TYPE(bec_type) :: becp0, becp_bp
   COMPLEX(DP), ALLOCATABLE :: q_dk_so(:,:,:,:)
   !
   COMPLEX(DP), ALLOCATABLE :: zetas(:,:) ! string data for phase control
   INTEGER, EXTERNAL :: find_free_unit
   INTEGER :: iun_phase
   INTEGER :: idumm1, idumm2
   REAL(DP) :: zetam
   CHARACTER(len=80) :: iun_name
   !
   !  -------------------------------------------------------------------------   !
   !                               INITIALIZATIONS                                !
   !  -------------------------------------------------------------------------   !
   !
   CALL start_clock( 'c_phase_field' )
   !
   SELECT CASE( pdir )
      CASE( 1 )
         iun_name='1'
      CASE( 2 )
         iun_name='2'
      CASE( 3 )
         iun_name='3'
   END SELECT
   !
   ALLOCATE( psi1(npol*npwx,nbnd) )
   ALLOCATE( psi(npol*npwx,nbnd)  )
   ALLOCATE( aux(ngm,nbnd)  )
   ALLOCATE( aux0(ngm,nbnd) )
   nspinnc=nspin
   IF (noncolin) THEN
      nspinnc=1
      ALLOCATE( aux_2(ngm,nbnd)  )
      ALLOCATE( aux0_2(ngm,nbnd) )
   END IF
   ALLOCATE( map_g(npwx) )
   ALLOCATE( l_cal(nbnd) )
   !
   IF ( pdir==3 ) THEN
      l_para=.FALSE.
   ELSE
      l_para=.TRUE.
   ENDIF
   !
   !
   IF (okvan) THEN
      CALL allocate_bec_type( nkb, nbnd, becp0 )
      CALL allocate_bec_type( nkb, nbnd, becp_bp )
      IF (lspinorb) ALLOCATE( q_dk_so(nhm,nhm,4,ntyp) )
   ENDIF
   !
   !
   pola=0.d0 !set to 0 electronic polarization   
   zeta_tot=(1.d0,0.d0)
   !
   !  --- Check that we are working with an insulator with no empty bands ---
   IF ( degauss > 0.0_DP )  CALL errore( 'c_phase_field', &
           'Polarization only for insulators and no empty bands', 1 )
   !
   !  --- Define a small number ---
   eps=1.0E-6_DP
   !
   !  --- Recalculate FFT correspondence (see ggen.f90) ---
   ALLOCATE( ln(-dfftp%nr1:dfftp%nr1, -dfftp%nr2:dfftp%nr2, -dfftp%nr3:dfftp%nr3) )
   DO ng=1,ngm
      mk1=NINT( g(1,ng)*at(1,1)+g(2,ng)*at(2,1)+g(3,ng)*at(3,1) )
      mk2=NINT( g(1,ng)*at(1,2)+g(2,ng)*at(2,2)+g(3,ng)*at(3,2) )
      mk3=NINT( g(1,ng)*at(1,3)+g(2,ng)*at(2,3)+g(3,ng)*at(3,3) )
      ln(mk1,mk2,mk3) = ng
   END DO
   !
   IF (okvan) THEN
      !  --- Initialize arrays ---
      jkb_bp=0
      DO nt=1,ntyp
         DO na=1,nat
            IF (ityp(na) == nt) THEN
               DO i=1, nh(nt)
                  jkb_bp=jkb_bp+1
                  nkbtona(jkb_bp) = na
                  nkbtonh(jkb_bp) = i        
               END DO
            END IF
         END DO
      END DO
   ENDIF
   !  --- Get the number of strings ---
   nstring=nks/nppstr_3d(pdir)
   nkort=nstring/(nspinnc) ! Include noncollinear case
   !
   !  --- Allocate memory for arrays ---
   ALLOCATE( phik(nstring)  )
   ALLOCATE( loc_k(nstring) )
   ALLOCATE( cphik(nstring) )
   ALLOCATE( wstring(nstring)  )
   ALLOCATE( pdl_elec(nstring) )
   ALLOCATE( mod_elec(nstring) )

   ALLOCATE( zetas(nkort,nspinnc) )
   !  -------------------------------------------------------------------------   !
   !           electronic polarization: set values for k-points strings           !
   !  -------------------------------------------------------------------------   !

   !  --- Find vector along strings ---
   IF ( nppstr_3d(pdir) /= 1) THEN
      gpar(1) = (xk(1,nx_el(nppstr_3d(pdir),pdir))-xk(1,nx_el(1,pdir)))*&
                 &DBLE(nppstr_3d(pdir))/DBLE(nppstr_3d(pdir)-1)
      gpar(2) = (xk(2,nx_el(nppstr_3d(pdir),pdir))-xk(2,nx_el(1,pdir)))*&
                 &DBLE(nppstr_3d(pdir))/DBLE(nppstr_3d(pdir)-1)
      gpar(3) = (xk(3,nx_el(nppstr_3d(pdir),pdir))-xk(3,nx_el(1,pdir)))*&
                 &DBLE(nppstr_3d(pdir))/DBLE(nppstr_3d(pdir)-1)
      gvec = DSQRT( gpar(1)**2+gpar(2)**2+gpar(3)**2 )*tpiba
   ELSE
      gpar(1) = 0.d0
      gpar(2) = 0.d0
      gpar(3) = 0.d0
      gpar(pdir) = 1.d0 / at(pdir,pdir)
      gvec = tpiba / SQRT( at(pdir,1)**2.d0 + at(pdir,2)**2.d0 + at(pdir,3)**2.d0 )
   ENDIF      
   !  --- Find vector between consecutive points in strings ---
   IF ( nppstr_3d(pdir) /= 1) THEN  ! orthorhombic cell 
      dk(1) = xk(1,nx_el(2,pdir))-xk(1,nx_el(1,pdir))
      dk(2) = xk(2,nx_el(2,pdir))-xk(2,nx_el(1,pdir)) 
      dk(3) = xk(3,nx_el(2,pdir))-xk(3,nx_el(1,pdir))
      dkmod = SQRT( dk(1)**2 + dk(2)**2 + dk(3)**2 )*tpiba 
   ELSE ! Gamma point case, only cubic cell for now
      dk(1) = 0.d0
      dk(2) = 0.d0
      dk(3) = 0.d0
      dk(pdir) = 1.d0 / at(pdir,pdir)
      dkmod = tpiba / SQRT( at(pdir,1)**2.d0 + at(pdir,2)**2.d0 + at(pdir,3)**2.d0 )
   ENDIF
   !
   !  -------------------------------------------------------------------------   !
   !                   electronic polarization: weight strings                    !
   !  -------------------------------------------------------------------------   !
   !
   !  --- Calculate string weights, normalizing to 1 (no spin) or 1+1 (spin) ---
   DO is=1,nspinnc ! Include noncollinear case
      weight = 0.0_dp
      DO kort=1,nkort
         istring = kort + (is-1)*nkort
         wstring(istring) = wk(nppstr_3d(pdir)*istring)
         weight = weight + wstring(istring)
      END DO
      DO kort=1,nkort
         istring = kort+(is-1)*nkort
         wstring(istring) = wstring(istring) / weight
      END DO
   END DO  
   !
   !  -------------------------------------------------------------------------   !
   !                  electronic polarization: structure factor                   !
   !  -------------------------------------------------------------------------   !
   !   
   !  --- Calculate structure factor e^{-i dk*R} ---
   !
   DO na=1,nat
      fac = (dk(1)*tau(1,na)+dk(2)*tau(2,na)+dk(3)*tau(3,na))*tpi 
      struc(na) = CMPLX( COS(fac),-SIN(fac),KIND=DP)
   END DO
   !
   !  -------------------------------------------------------------------------   !
   !                     electronic polarization: form factor                     !
   !  -------------------------------------------------------------------------   !
   IF (okvan) THEN
      !  --- Bessel transform of Q_ij(|r|) at dk [Q_ij^L(|r|)] in array qrad ---
      ! CALL calc_btq(dkmod,qrad_dk,0) is no longer needed, see bp_c_phase
      CALL compute_qqc ( tpiba, dk, omega, q_dk )
      IF (lspinorb) CALL transform_qq_so(q_dk,q_dk_so)
      !
   ENDIF
   !
   !  -------------------------------------------------------------------------    !
   !                   electronic polarization: strings phases                     !
   !  -------------------------------------------------------------------------    !
   !
   el_loc=0.d0
   kpoint=0
   zeta=(1.d0,0.d0)
   !
   if(ionode .AND. phase_control>0) THEN
      iun_phase=find_free_unit()
      IF ( phase_control==1 ) THEN
         OPEN( iun_phase, file=TRIM(tmp_dir)//TRIM(prefix)//'.phase.data'//TRIM(iun_name),status='unknown')
      ELSE
         OPEN( iun_phase, file=TRIM(tmp_dir)//TRIM(prefix)//'.phase.data'//TRIM(iun_name),status='OLD')
         DO is=1,nspinnc
            DO kort=1,nkort
               read(iun_phase,*) idumm1,idumm2,zetas(kort,is)
               zetam = DBLE( CONJG(zetas(kort,is))*zetas(kort,is) )
               zetam = 1.d0 / DSQRT( zetam )
               zetas(kort,is) = CONJG( zetam*zetas(kort,is) )
            ENDDO
         ENDDO
      ENDIF
   ENDIF
   !
   ! FIXME: not sure which is the proper communicator here
   !
   IF (phase_control==2) &
      CALL mp_bcast(zetas,   ionode_id, intra_pool_comm )
   !
   !  --- Start loop over spin ---
   DO is=1,nspinnc ! Include noncollinear case 
      !
      ! l_cal(n) = .true./.false. if n-th state is occupied/empty
      DO nb = 1, nbnd
         IF ( nspin == 2 .AND. tfixed_occ) THEN
            l_cal(nb) = ( f_inp(nb,is) /= 0.0_dp )
         ELSE
             IF (noncolin) THEN
                l_cal(nb) = ( nb <= NINT( nelec) )
             ELSE
                l_cal(nb) = ( nb <= NINT( nelec/2.0_dp ) )
             ENDIF
         ENDIF
      END DO
      ! 
      !     --- Start loop over orthogonal k-points ---
      !
      DO kort=1,nkort
         zeta_loc=(1.d0,0.d0)
         !        --- Index for this string ---
         istring=kort+(is-1)*nkort
         !
         !   --- Initialize expectation value of the phase operator ---
         !
         zeta_mod = 1.d0
         !
         !   --- Start loop over parallel k-points ---
         DO kpar = 1,nppstr_3d(pdir)+1
            !
            !    --- Set index of k-point ---
            kpoint = kpoint + 1
            !
            !    --- Calculate dot products between wavefunctions and betas ---
            !
            IF (kpar /= 1 ) THEN
               !
               !   --- Dot wavefunctions and betas for PREVIOUS k-point ---
               !
               ik = nx_el(kpoint-1,pdir)
               npw0   = ngk(ik)
               igk0(:)= igk_k(:,ik)
               CALL get_buffer( psi,nwordwfc,iunwfc,nx_el(kpoint-1,pdir) )
               IF (okvan) THEN
                  CALL init_us_2( npw0,igk0,xk(1,nx_el(kpoint-1,pdir)),vkb )
                  CALL calbec( npw0, vkb, psi, becp0)
               ENDIF
               !
               !   --- Dot wavefunctions and betas for CURRENT k-point ---
               !
               IF (kpar /= (nppstr_3d(pdir)+1)) THEN
                  ik = nx_el(kpoint,pdir)
                  npw1   = ngk(ik)
                  igk1(:)= igk_k(:,ik)
                  CALL get_buffer( psi1,nwordwfc,iunwfc,nx_el(kpoint,pdir) )
                  IF (okvan) THEN
                     CALL init_us_2( npw1,igk1,xk(1,nx_el(kpoint,pdir)),vkb )
                     CALL calbec( npw1, vkb, psi1, becp_bp )
                  ENDIF
               ELSE
                  kstart = kpoint-(nppstr_3d(pdir)+1)+1
                  ik = nx_el(kstart,pdir)
                  npw1   = ngk(ik)
                  igk1(:)= igk_k(:,ik)
                  CALL get_buffer( psi1,nwordwfc,iunwfc,nx_el(kstart,pdir) )
                  IF (okvan) THEN
                     CALL init_us_2( npw1,igk1,xk(1,nx_el(kstart,pdir)),vkb )
                     CALL calbec( npw1, vkb, psi1, becp_bp )
                  ENDIF
               ENDIF
               !
               !           --- Matrix elements calculation ---
               !
               IF (kpar == (nppstr_3d(pdir)+1) .AND. .NOT.l_para) THEN
                  map_g(:) = 0
                  DO ig = 1, npw1
                     !           --- If k'=k+G_o, the relation psi_k+G_o (G-G_o) ---
                     !           --- = psi_k(G) is used, gpar=G_o, gtr = G-G_o ---
                     !          
                     gtr(1)=g(1,igk1(ig)) - gpar(1)
                     gtr(2)=g(2,igk1(ig)) - gpar(2) 
                     gtr(3)=g(3,igk1(ig)) - gpar(3) 
                     !           --- Find crystal coordinates of gtr, n1,n2,n3 ---
                     !           --- and the position ng in the ngm array ---
                     IF (gtr(1)**2+gtr(2)**2+gtr(3)**2 <= gcutm) THEN
                        n1 = NINT(gtr(1)*at(1,1)+gtr(2)*at(2,1) &
                              +gtr(3)*at(3,1))
                        n2 = NINT(gtr(1)*at(1,2)+gtr(2)*at(2,2) &
                              +gtr(3)*at(3,2))
                        n3 = NINT(gtr(1)*at(1,3)+gtr(2)*at(2,3) &
                              +gtr(3)*at(3,3))
                        ng = ln(n1,n2,n3) 
                        IF ( (ABS(g(1,ng)-gtr(1)) > eps) .OR. &
                             (ABS(g(2,ng)-gtr(2)) > eps) .OR. &
                             (ABS(g(3,ng)-gtr(3)) > eps) ) THEN
                           WRITE (6,*) ' error: translated G=', &
                                gtr(1),gtr(2),gtr(3), &
                                &     ' with crystal coordinates',n1,n2,n3, &
                                &     ' corresponds to ng=',ng,' but G(ng)=', &
                                &     g(1,ng),g(2,ng),g(3,ng)
                           WRITE (6,*) ' probably because G_par is NOT', &
                                &    ' a reciprocal lattice vector '
                           WRITE (6,*) ' Possible choices as smallest ', &
                                ' G_par:'
                           DO i=1,50
                              WRITE (6,*) ' i=',i,'   G=', &
                                   g(1,i),g(2,i),g(3,i)
                           ENDDO
                           STOP
                        ENDIF
                     ELSE 
                        WRITE (6,*) ' |gtr| > gcutm  for gtr=', &
                             gtr(1),gtr(2),gtr(3) 
                        STOP
                     END IF
                     map_g(ig)=ng
                  ENDDO
               ENDIF

               mat=(0.d0,0.d0)
               aux=(0.d0,0.d0)
               IF ( noncolin ) aux_2=(0.d0,0.d0)
               DO mb=1,nbnd
                  IF ( l_cal(mb) ) THEN
                     IF (kpar /= (nppstr_3d(pdir)+1)) THEN
                        DO ig=1,npw1
                           aux(igk1(ig),mb)=psi1(ig,mb)
                           IF (noncolin) aux_2(igk1(ig),mb)=psi1(ig+npwx,mb)
                        ENDDO
                     ELSE IF( .not. l_para) THEN
                        DO ig=1,npw1
                           aux(map_g(ig),mb)=psi1(ig,mb)
                           IF (noncolin) aux_2(map_g(ig),mb)=psi1(ig+npwx,mb)
                        ENDDO
                     ELSE
                        ! allocate global array
                        ALLOCATE (aux_g(ngm_g))
                        IF(noncolin) ALLOCATE (aux_g_2(ngm_g))
                        aux_g=(0.d0,0.d0)
                        IF(noncolin) aux_g_2=(0.d0,0.d0)
                        ! put psi1 on global array
                        DO ig=1,npw1
                           aux_g(mapgm_global(ig_l2g(igk1(ig)),pdir))=psi1(ig,mb)
                           IF(noncolin) aux_g_2(mapgm_global(ig_l2g(igk1(ig)),pdir))=psi1(ig+npwx,mb)
                        ENDDO
                        CALL mp_sum(aux_g(:),intra_bgrp_comm)
                        IF (noncolin) CALL mp_sum(aux_g_2(:),intra_bgrp_comm) !non-collinear
                        DO ig=1,ngm
                           aux(ig,mb) = aux_g(ig_l2g(ig))
                           IF (noncolin) aux_2(ig,mb) = aux_g_2(ig_l2g(ig))
                        ENDDO
                        DEALLOCATE (aux_g)
                        IF(noncolin) DEALLOCATE (aux_g_2)
                     END IF
                  END IF
               END DO
               aux0=(0.d0,0.d0)
               IF (noncolin) aux0_2 = (0.d0,0.d0)
               DO nb=1,nbnd
                  DO ig=1,npw0
                     aux0(igk0(ig),nb) = psi(ig,nb)
                     IF (noncolin) aux0_2(igk0(ig),nb) = psi(ig+npwx,nb)
                  END DO
               ENDDO
               CALL ZGEMM('C','N',nbnd,nbnd,ngm,(1.d0,0.d0),aux0,ngm,aux,ngm,(0.d0,0.d0),mat,nbnd)
               IF (noncolin) CALL ZGEMM('C','N',nbnd,nbnd,ngm,(1.d0,0.d0),aux0_2,ngm,aux_2,ngm,(1.d0,0.d0),mat,nbnd)
               !
               CALL  mp_sum( mat, intra_bgrp_comm )
               !
               ! --- Calculate the augmented part: ij=KB projectors, ---
               ! --- R=atom index: SUM_{ijR} q(ijR) <u_nk|beta_iR>   ---
               ! --- <beta_jR|u_mk'> e^i(k-k')*R =                   ---
               ! --- also <u_nk|beta_iR>=<psi_nk|beta_iR> = becp^*   ---
               !
               IF (okvan) THEN
                  DO mb=1,nbnd
                     DO nb=1,nbnd
                        IF ( l_cal(mb) .AND. l_cal(nb) ) THEN
                           pref = (0.d0,0.d0)
                           DO jkb=1,nkb
                              nhjkb = nkbtonh(jkb)
                              na = nkbtona(jkb)
                              np = ityp(na)
                              nhjkbm = nh(np)
                              jkb1 = jkb - nhjkb
                              DO j = 1,nhjkbm
                                 IF (lspinorb) THEN
                                    pref = pref+CONJG(becp0%nc(jkb,1,nb))*becp_bp%nc(jkb1+j,1,mb) &
                                         *q_dk_so(nhjkb,j,1,np)*struc(na)
                                    pref = pref+CONJG(becp0%nc(jkb,1,nb))*becp_bp%nc(jkb1+j,2,mb) &
                                         *q_dk_so(nhjkb,j,2,np)*struc(na)
                                    pref = pref+CONJG(becp0%nc(jkb,2,nb))*becp_bp%nc(jkb1+j,1,mb) &
                                         *q_dk_so(nhjkb,j,3,np)*struc(na)
                                    pref = pref+CONJG(becp0%nc(jkb,2,nb))*becp_bp%nc(jkb1+j,2,mb) &
                                         *q_dk_so(nhjkb,j,4,np)*struc(na)
                                    !
                                 ELSE
                                    !
                                    pref = pref+CONJG(becp0%k(jkb,nb))*becp_bp%k(jkb1+j,mb) &
                                         *q_dk(nhjkb,j,np)*struc(na)
                                 ENDIF
                              ENDDO
                           ENDDO
                           !
                           mat(nb,mb) = mat(nb,mb) + pref
                        ENDIF
                     ENDDO
                  ENDDO
               ENDIF
               !
               DO mb=1,nbnd
                  IF ( .NOT.l_cal(mb) ) mat(mb,mb)=(1.d0,0.d0)
               END DO
               !
               !  --- Calculate matrix determinant ---
               !
               CALL ZGETRF( nbnd,nbnd,mat,nbnd,ivpt,info )
               CALL errore('c_phase_field','error in zgetrf', ABS(info))
               det=(1.d0,0.d0)
               DO nb=1,nbnd
                  IF (nb /= ivpt(nb)) det = -det
                  det = det*mat(nb,nb)
               ENDDO
               !  --- Multiply by the already calculated determinants ---
               zeta=zeta*det
               zeta_loc=zeta_loc*det
               !
               !  --- End of dot products between wavefunctions and betas ---
            ENDIF
            !
         !  --- End loop over parallel k-points ---
         END DO 
         !
         IF (phase_control==1) THEN
            IF (ionode) write(iun_phase,*) kort,is,zeta_loc
         ELSE IF (phase_control==2) THEN
            zeta_loc=zeta_loc*zetas(kort,is)
         ENDIF
         !
         zeta_tot=zeta_tot*(zeta_loc**wstring(istring))
         !
         ! ... uncomment the following line for printing string data
         !WRITE (stdout,*) 'String :', kort, zeta_loc
         !
         pola = pola + wstring(istring)*AIMAG(LOG(zeta_loc))
         !
         kpoint = kpoint-1
         ! --- Calculate the phase for this string ---
         phik(istring) = AIMAG(LOG(zeta))
         cphik(istring) = COS(phik(istring))*(1.0_DP,0.0_DP) &
                         +SIN(phik(istring))*(0.0_DP,1.0_DP)
         !
         ! --- Calculate the localization for current kort ---
         zeta_mod = DBLE(CONJG(zeta)*zeta)
         loc_k(istring) = -(nppstr_3d(pdir)-1) / gvec**2 / nbnd *LOG(zeta_mod)

      ! --- End loop over orthogonal k-points ---
      END DO
      !
   !  --- End loop over spin ---
   END DO
   !-----calculate polarization
   !-----the factor 2. is because of spin
   !new system for avoiding phases problem
   !pola=aimag(log(zeta_tot))
   !
   IF (nspin==1) pola=pola*2.d0
   !
   CALL factor_a( pdir, at, dkfact )
   !
   !factor sqrt(2) is the electronic charge in Rydberg units 
   pola = pola * DSQRT(2.d0) / tpiba * dkfact
   !
   !write output
   WRITE (stdout,*)
   WRITE (stdout,*) "    Expectation value of exp(iGx):",zeta_tot,dkfact
   WRITE (stdout,*) "    Electronic Dipole per cell (Ry a.u.)",pola
   !
   !  -------------------------------------------------------------------------   !
   !                              ionic polarization                              !
   !  -------------------------------------------------------------------------   !
   !
   !factor sqrt(2) is the electronic charge in Rydberg units
   pola_ion=0.d0
   DO na = 1, nat
     pola_ion = pola_ion + zv(ityp(na))*tau(pdir,na)*alat*DSQRT(2.d0)
   END DO
   !
   WRITE (stdout,*) "    Ionic Dipole per cell (Ry a.u.)",pola_ion
   !
   !
   el_pola   = pola
   ion_pola  = pola_ion
   fact_pola = DSQRT(2.d0)/tpiba*dkfact
   !
   IF ( ionode .AND. phase_control>0 ) close(iun_phase)
   !
   !  -------------------------------------------------------------------------   !
   !
   !  --- Free memory ---
   DEALLOCATE( l_cal    )
   DEALLOCATE( pdl_elec )
   DEALLOCATE( mod_elec )
   DEALLOCATE( wstring  )
   DEALLOCATE( loc_k )
   DEALLOCATE( phik  )
   DEALLOCATE( cphik )
   DEALLOCATE( ln    )
   DEALLOCATE( map_g )
   DEALLOCATE( aux   )
   DEALLOCATE( aux0  )
   DEALLOCATE( psi   )
   DEALLOCATE( psi1  )
   DEALLOCATE( zetas )
   IF (ALLOCATED(aux_2) ) DEALLOCATE( aux_2  )
   IF (ALLOCATED(aux0_2)) DEALLOCATE( aux0_2 )
   !
   IF (okvan) THEN
      CALL deallocate_bec_type( becp0 )
      CALL deallocate_bec_type( becp_bp )
      IF (lspinorb) DEALLOCATE( q_dk_so )
   ENDIF
   CALL stop_clock( 'c_phase_field' )
   !
   !------------------------------------------------------------------------------!
   !
END SUBROUTINE c_phase_field

!==============================================================================!
