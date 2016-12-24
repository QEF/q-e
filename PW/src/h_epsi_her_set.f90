!
! Copyright (C) 2005 Paolo Umari
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

!-----------------------------------------------------------------------
subroutine h_epsi_her_set(pdir, e_field)
  !-----------------------------------------------------------------------
  !
  ! this subroutine  builds the hermitean operators  w_k w_k*, 
  ! (as in Souza,et al.  PRB B 69, 085106 (2004))
  !
  ! wavefunctions from previous iteration are read into 'evcel'
  ! spin polarized systems supported only with fixed occupations

  USE noncollin_module,     ONLY : noncolin, npol
  USE spin_orb, ONLY: lspinorb
  USE kinds,    ONLY : DP
  USE us
  USE wvfct,    ONLY : npwx, nbnd
  USE ldaU,     ONLY : lda_plus_u
  USE lsda_mod, ONLY : current_spin, nspin
  USE scf,      ONLY : vrs  
  USE gvect
  USE fft_base, ONLY : dfftp
  USE uspp
  USE uspp_param, ONLY: upf, nh, nhm, nbetam, lmaxq
  USE bp,         ONLY : nppstr_3d, fact_hepsi, evcel, evcp=>evcelp, &
                         evcm=>evcelm, mapgp_global, mapgm_global, nx_el
  USE klist
  USE cell_base, ONLY: at, alat, tpiba, omega, bg
  USE ions_base, ONLY: ityp, tau, nat,ntyp => nsp
  USE io_files,  ONLY: iunwfc, nwordwfc, iunefieldm, iunefieldp
  USE buffers,   ONLY: get_buffer, save_buffer
  USE constants, ONLY : e2, pi, tpi, fpi
  USE fixed_occ
  USE mp,        ONLY : mp_sum
  USE mp_bands,  ONLY : intra_bgrp_comm
  USE becmod,    ONLY : bec_type, becp, calbec,allocate_bec_type, deallocate_bec_type
  !
  implicit none
  !
  INTEGER, INTENT(in) :: pdir!direction on which the polarization is calculated
  REAL(DP) :: e_field!electric field along pdir
  !
  !  --- Internal definitions ---

   COMPLEX(DP), ALLOCATABLE  :: evct(:,:)!for temporary wavefunctios
   INTEGER :: i,ipol
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
   REAL(dp) :: dk(3)
   REAL(dp) :: dkm(3)! -dk
   REAL(dp) :: dkmod
   REAL(dp) :: eps
   REAL(dp) :: fac
   REAL(dp) :: gpar(3)
   REAL(dp) :: gtr(3)
   !REAL(dp) :: gvec
   REAL(dp), ALLOCATABLE :: ln(:,:,:)
   REAL(dp), ALLOCATABLE  :: ln0(:,:,:)!map g-space global to g-space k-point dependent
   REAL(dp) :: qrad_dk(nbetam,nbetam,lmaxq,ntyp)
   REAL(dp) :: ylm_dk(lmaxq*lmaxq)
   COMPLEX(dp), ALLOCATABLE :: aux(:)
   COMPLEX(dp), ALLOCATABLE  :: aux0(:)
   ! Also for noncollinear calculation
   COMPLEX(DP), ALLOCATABLE :: aux_2(:)
   COMPLEX(DP), ALLOCATABLE :: aux0_2(:),aux0vec(:,:),aux1vec(:,:)
  
   COMPLEX(dp) :: cdet(2)
   COMPLEX(dp) :: cdwork(nbnd)
   COMPLEX(dp), ALLOCATABLE :: mat(:,:)
   COMPLEX(dp) :: pref
   COMPLEX(dp) :: q_dk(nhm,nhm,ntyp)
   COMPLEX(dp) :: q_dkp(nhm,nhm,ntyp)!to store the terms T^dagger e^(iGx) T
   COMPLEX(dp) :: struc(nat)
   COMPLEX(dp) :: zdotc

   
   COMPLEX(dp) :: sca,sca1
   COMPLEX(dp) :: ps(nkb,nbnd*npol)
   COMPLEX(dp) :: matbig(nks,nbnd,nbnd)
   INTEGER :: mdone(nks)
   INTEGER :: ijkb0,  ibnd,jh, ih, ikb, ik, ikk


   LOGICAL, ALLOCATABLE :: l_cal(:) ! flag for empty/occupied states
   INTEGER, ALLOCATABLE  :: map_g(:)
   TYPE(bec_type) :: becp_bp,becp0
     

   REAL(dp) :: dkfact
   LOGICAL  :: l_para! if true new parallel treatment
   COMPLEX(kind=DP), ALLOCATABLE :: aux_g(:)
   COMPLEX(kind=DP), ALLOCATABLE :: aux_g_2(:) ! non-collinear case

   COMPLEX(DP), ALLOCATABLE :: q_dk_so(:,:,:,:),q_dkp_so(:,:,:,:)



   !  --- Define a small number ---
   eps=0.000001d0

   allocate(mat(nbnd,nbnd))

   if(ABS(e_field)<eps) return
   call start_clock('h_epsi_set')



  if(okvan) then
      CALL allocate_bec_type (nkb,nbnd,becp0)
     CALL allocate_bec_type (nkb,nbnd,becp_bp)
     IF (lspinorb) ALLOCATE(q_dk_so(nhm,nhm,4,ntyp))
     IF (lspinorb) ALLOCATE(q_dkp_so(nhm,nhm,4,ntyp))
  endif


!  -------------------------------------------------------------------------   !
!                               INITIALIZATIONS
!  -------------------------------------------------------------------------   !
  
   if(pdir==3) then
      l_para=.false.
   else
      l_para=.true.
   endif


   ALLOCATE( evct(npwx*npol,nbnd))
   ALLOCATE( map_g(npwx))

   ALLOCATE(ln(-dfftp%nr1:dfftp%nr1,-dfftp%nr2:dfftp%nr2,-dfftp%nr3:dfftp%nr3),&
            ln0(-dfftp%nr1:dfftp%nr1,-dfftp%nr2:dfftp%nr2,-dfftp%nr3:dfftp%nr3))
   ALLOCATE(aux(ngm),aux0(ngm))
   IF(noncolin) ALLOCATE (aux_2(ngm),aux0_2(ngm))

   ALLOCATE (l_cal(nbnd))

!  --- Recalculate FFT correspondence (see ggen.f90) ---
      DO ng=1,ngm
         mk1=nint(g(1,ng)*at(1,1)+g(2,ng)*at(2,1)+g(3,ng)*at(3,1))
         mk2=nint(g(1,ng)*at(1,2)+g(2,ng)*at(2,2)+g(3,ng)*at(3,2))
         mk3=nint(g(1,ng)*at(1,3)+g(2,ng)*at(2,3)+g(3,ng)*at(3,3))
         ln(mk1,mk2,mk3) = ng
      END DO
 
! !  --- Find vector along strings ---
      if(nppstr_3d(pdir) .ne. 1) then
         gpar(1)=(xk(1,nx_el(nppstr_3d(pdir),pdir))-xk(1,nx_el(1,pdir)))*&
              &DBLE(nppstr_3d(pdir))/DBLE(nppstr_3d(pdir)-1)
         gpar(2)=(xk(2,nx_el(nppstr_3d(pdir),pdir))-xk(2,nx_el(1,pdir)))*&
              &DBLE(nppstr_3d(pdir))/DBLE(nppstr_3d(pdir)-1)
         gpar(3)=(xk(3,nx_el(nppstr_3d(pdir),pdir))-xk(3,nx_el(1,pdir)))*&
              &DBLE(nppstr_3d(pdir))/DBLE(nppstr_3d(pdir)-1)

      else
         gpar=bg(:,pdir)
      endif
      !gvec=dsqrt(gpar(1)**2+gpar(2)**2+gpar(3)**2)*tpiba
      
      matbig=(0.0d0,0.0d0)
      mdone=0
     
      call factor_a(pdir,at,dkfact)
      dkfact=tpiba/dkfact/dble(nppstr_3d(pdir))
  
!determines the spin polarization

   DO ik=1,nks
  
!  --- Find vector between consecutive points in strings ---
      dk=gpar/nppstr_3d(pdir)
      dkmod=SQRT(DOT_PRODUCT(dk,dk))*tpiba
      dkm=-dk

      CALL get_buffer ( evcel, nwordwfc, iunwfc, nx_el(ik,pdir) )

      if(nspin==2) then
         if(ik <= nks/2) then
            is = 1
         else
            is = 2
         endif
      else
         is = 1
      end if

      ! l_cal(n) = .true./.false. if n-th state is occupied/empty
      DO nb = 1, nbnd
         IF ( nspin == 2 .AND. tfixed_occ) THEN
            l_cal(nb) = ( f_inp(nb,is) /= 0.0_dp )
         ELSE
            IF (noncolin) THEN
                l_cal(nb) = ( nb <= NINT ( nelec) )
             ELSE
                l_cal(nb) = ( nb <= NINT ( nelec/2.0_dp ) )
             END IF
         ENDIF
      END DO

      ik_stringa=mod(ik-1,nppstr_3d(pdir))+1
     
 
      if(okvan) then
!  --- Initialize arrays ---
         jkb_bp=0
         DO nt=1,ntyp
            DO na=1,nat
               IF (ityp(na).eq.nt) THEN
                  DO i=1, nh(nt)
                     jkb_bp=jkb_bp+1
                     nkbtona(jkb_bp) = na
                     nkbtonh(jkb_bp) = i
                  END DO
               END IF
            END DO
         END DO
      endif

!  -------------------------------------------------------------------------   !
!           electronic polarization: set values for k-points strings           !
!  -------------------------------------------------------------------------   !

!calculates fact factor
!electronic charge is sqrt(2.) (Rydberg units)
!the factor (-i)/2 comes form operator Im

      fact_hepsi(nx_el(ik,pdir),pdir)=(0.d0,-1.d0)*e_field*dsqrt(2.d0)/2.d0/dkfact
  
      evcm(:,:,pdir)=(0.d0,0.d0)
      evcp(:,:,pdir)=(0.d0,0.d0)

      if(okvan) then
!  -------------------------------------------------------------------------   !
!                  electronic polarization: structure factor                   !
!  -------------------------------------------------------------------------   !

!  --- Calculate structure factor e^{-i dk*R} ---
         DO na=1,nat
            fac=(dk(1)*tau(1,na)+dk(2)*tau(2,na)+dk(3)*tau(3,na))*tpi 
            struc(na)=CMPLX(cos(fac),-sin(fac),kind=DP)
         END DO


!  -------------------------------------------------------------------------   !
!                     electronic polarization: form factor                     !
!  -------------------------------------------------------------------------   !

!  --- Calculate Bessel transform of Q_ij(|r|) at dk [Q_ij^L(|r|)] ---
         CALL calc_btq(dkmod,qrad_dk,0)
   
!  --- Calculate the q-space real spherical harmonics at dk [Y_LM] --- 
         dkmod=dk(1)**2+dk(2)**2+dk(3)**2
         CALL ylmr2(lmaxq*lmaxq, 1, dk, dkmod, ylm_dk)
      
!  --- Form factor: 4 pi sum_LM c_ij^LM Y_LM(Omega) Q_ij^L(|r|) ---
         q_dk=(0.d0,0.d0)
         DO np =1, ntyp
            if( upf(np)%tvanp ) then
               DO iv = 1, nh(np)
                  DO jv = iv, nh(np)
                     call qvan3(iv,jv,np,pref,ylm_dk,qrad_dk)
                     q_dk(iv,jv,np) = omega*pref
                     q_dk(jv,iv,np) = omega*pref
                  ENDDO
               ENDDO
            endif
         ENDDO
         IF (lspinorb) CALL transform_qq_so(q_dk,q_dk_so)
!  --- Calculate the q-space real spherical harmonics at -dk [Y_LM] --- 
         dkmod=dkm(1)**2+dkm(2)**2+dkm(3)**2
         CALL ylmr2(lmaxq*lmaxq, 1, dkm, dkmod, ylm_dk)
  
!  --- Form factor: 4 pi sum_LM c_ij^LM Y_LM(Omega) Q_ij^L(|r|) ---

         q_dkp=(0.d0,0.d0)
         DO np =1, ntyp
            if( upf(np)%tvanp ) then
               DO iv = 1, nh(np)
                  DO jv = iv, nh(np)
                     call qvan3(iv,jv,np,pref,ylm_dk,qrad_dk)
                     q_dkp(iv,jv,np) = omega*pref
                     q_dkp(jv,iv,np) = omega*pref
                  ENDDO
               ENDDO
            endif
         ENDDO
         IF (lspinorb) CALL transform_qq_so(q_dkp,q_dkp_so)
      
      endif

!  -------------------------------------------------------------------------   !
!                   electronic polarization: strings phases                    !
!  -------------------------------------------------------------------------   !
  

 
       
!    

!calculate the term  S-1(k,k-1)       


!       
      if(ik_stringa /= 1) then

         ikk = nx_el(ik-1,pdir)
         npw0   = ngk(ikk)
         igk0(:)= igk_k(:,ikk)

         CALL get_buffer (evct,nwordwfc,iunwfc,nx_el(ik-1,pdir))
!        
!           --- Calculate dot products between wavefunctions

!              --- Dot wavefunctions and betas for PREVIOUS k-point ---
         if(okvan) then
            CALL init_us_2 (npw0,igk0,xk(1,nx_el(ik-1,pdir)),vkb)
            CALL calbec( npw0, vkb, evct, becp0 )
         endif
!              --- Dot wavefunctions and betas for CURRENT k-point ---
         
         ikk = nx_el(ik,pdir)
         npw1   = ngk(ikk)
         igk1(:)= igk_k(:,ikk)

         !  --- Recalculate FFT correspondence (see ggen.f90) ---

         ln0=0!set array to 0
         DO ig=1,npw1
            mk1=nint(g(1,igk1(ig))*at(1,1)+g(2,igk1(ig))*at(2,1)+g(3,igk1(ig))*at(3,1))
            mk2=nint(g(1,igk1(ig))*at(1,2)+g(2,igk1(ig))*at(2,2)+g(3,igk1(ig))*at(3,2))
            mk3=nint(g(1,igk1(ig))*at(1,3)+g(2,igk1(ig))*at(2,3)+g(3,igk1(ig))*at(3,3))
            ln0(mk1,mk2,mk3) = ig
         END DO
      
         if(okvan) then
            CALL init_us_2 (npw1,igk1,xk(1,nx_el(ik,pdir)),vkb)
            CALL calbec( npw1, vkb, evcel, becp_bp )
         endif




!              --- Matrix elements calculation ---

        IF(mdone(nx_el(ik,pdir))==0) THEN
         mat=(0.d0,0.d0)
         DO nb=1,nbnd
            DO mb=1,nbnd
               IF ( .NOT. l_cal(nb) .OR. .NOT. l_cal(mb) ) THEN
                  IF ( nb == mb )  mat(nb,mb)=0.d0
               ELSE
                  aux=(0.d0,0.d0)
                  aux0=(0.d0,0.d0)
                  if(noncolin) then
                     aux_2=(0.d0,0.d0)
                     aux0_2=(0.d0,0.d0)
                  end if
                  DO ig=1,npw1
                     aux0(igk1(ig))=evcel(ig,nb)
                     IF (noncolin) aux0_2(igk1(ig))=evcel(ig+npwx,nb)
                  END DO

                  DO ig=1,npw0
                     aux(igk0(ig))=evct(ig,mb)
                     IF (noncolin) aux_2(igk0(ig))=evct(ig+npwx,mb)
                  END DO
                   
                  mat(nb,mb) = zdotc(ngm,aux0,1,aux,1)
                  if (noncolin) mat(nb,mb)=mat(nb,mb)+zdotc(ngm,aux0_2,1,aux_2,1)
!                    --- Calculate the augmented part: ij=KB projectors, ---
!                    --- R=atom index: SUM_{ijR} q(ijR) <u_nk|beta_iR>   ---
!                    --- <beta_jR|u_mk'> e^i(k-k')*R =                   ---
!                    --- also <u_nk|beta_iR>=<psi_nk|beta_iR> = becp^*   ---
               end if
            END DO
         END DO
         !
         call mp_sum( mat, intra_bgrp_comm )
         DO nb=1,nbnd
            IF(.not.l_cal(nb)) mat(nb,nb)=1.d0
         END DO
         !
         DO nb=1,nbnd
            DO mb=1,nbnd
               IF ( l_cal(nb) .AND. l_cal(mb) ) THEN
                  if(okvan) then
                     pref = (0.d0,0.d0)
                     DO jkb=1,nkb
                        nhjkb = nkbtonh(jkb)
                        na = nkbtona(jkb)
                        np = ityp(na)
                        nhjkbm = nh(np)
                        jkb1 = jkb - nhjkb
                        if(lspinorb) THEN
                           DO j = 1,nhjkbm
                              pref = pref+CONJG(becp_bp%nc(jkb,1,nb))* becp0%nc(jkb1+j,1,mb)  &
                                   *q_dkp_so(nhjkb,j,1,np)*CONJG(struc(na))   
                              pref= pref+CONJG(becp_bp%nc(jkb,1,nb))* becp0%nc(jkb1+j,2,mb)  &
                                   *q_dkp_so(nhjkb,j,2,np)*CONJG(struc(na)) 
                              pref= pref +CONJG(becp_bp%nc(jkb,2,nb))* becp0%nc(jkb1+j,1,mb)  &
                                   *q_dkp_so(nhjkb,j,3,np)*CONJG(struc(na)) 
                              pref= pref +CONJG(becp_bp%nc(jkb,2,nb))* becp0%nc(jkb1+j,2,mb)   &
                                   *q_dkp_so(nhjkb,j,4,np)*CONJG(struc(na))

                           ENDDO

                        ELSE
                           DO j = 1,nhjkbm
                              pref = pref+CONJG(becp_bp%k(jkb,nb))*becp0%k(jkb1+j,mb) &
                                   &        *q_dkp(nhjkb,j,np)*CONJG(struc(na))
                           ENDDO

                        END IF
                     ENDDO
                     mat(nb,mb) = mat(nb,mb) + pref
                  endif
               ENDIF
            ENDDO
         ENDDO

!              --- Calculate matrix inverse ---
         CALL zgefa(mat,nbnd,nbnd,ivpt,info)
         CALL errore('h_epsi_her_set','error in zgefa',abs(info))
         CALL zgedi(mat,nbnd,nbnd,ivpt,cdet,cdwork,1)
         matbig(nx_el(ik,pdir),:,:)=mat
         mdone(nx_el(ik,pdir))=1
        ELSE
         mat=matbig(nx_el(ik,pdir),:,:)
        END IF

   
        DO nb = 1, nbnd
           DO mb = 1, nbnd
              if(.not.l_cal(nb).or. .not. l_cal(mb)) mat(mb,nb)=(0.d0,0.d0)
           END DO
        END DO
        

!    mat=S^-1(k,k-1)
         do ig=1,npw0
            gtr(1)=g(1,igk0(ig))
            gtr(2)=g(2,igk0(ig))         
            gtr(3)=g(3,igk0(ig))         
            !                          --- Find crystal coordinates of gtr, n1,n2,n3 ---
            !                          --- and the position ng in the ngm array ---
            IF (gtr(1)**2+gtr(2)**2+gtr(3)**2 <= gcutm) THEN
               n1=NINT(gtr(1)*at(1,1)+gtr(2)*at(2,1) &
                    &    +gtr(3)*at(3,1))
               n2=NINT(gtr(1)*at(1,2)+gtr(2)*at(2,2) &
                    &   +gtr(3)*at(3,2))
               n3=NINT(gtr(1)*at(1,3)+gtr(2)*at(2,3) &
                    &   +gtr(3)*at(3,3))
               ng=ln0(n1,n2,n3)
               if(ng .gt. 0) then
                  do m=1,nbnd
                     do nb=1,nbnd
                        evcm(ng,m,pdir)=evcm(ng,m,pdir) + mat(nb,m)*evct(ig,nb)
                        IF (noncolin) evcm(ng+npwx,m,pdir)=evcm(ng+npwx,m,pdir) &
                          +mat(nb,m)*evct(ig+npwx,nb)
                     enddo
                  enddo
               end if
            ENDIF
         enddo

! add US terms into evcm
! calculate |beta_(ik,na,ih)>Q_dkp(na,ih,ij)<|beta_(ik-1,na,ih)| 
         if(okvan) then
            evct(:,:) =  (0.d0, 0.d0)
            ps (:,:) = (0.d0, 0.d0)
            ijkb0 = 0
            do nt = 1, ntyp
               do na = 1, nat
                  if (ityp (na) .eq.nt) then
                     do ibnd = 1, nbnd
                        do jh = 1, nh (nt)
                           jkb = ijkb0 + jh
                           do ih = 1, nh (nt)
                              ikb = ijkb0 + ih
                              if(lspinorb) then
                                   ps (ikb, (ibnd-1)*npol+1) = ps (ikb,(ibnd-1)*npol+1 ) + &
                                      q_dkp_so(ih,jh,1,ityp(na))*CONJG(struc(na))* becp0%nc(jkb,1,ibnd)
                                   ps (ikb, (ibnd-1)*npol+1) = ps (ikb,(ibnd-1)*npol+1 ) + &
                                        q_dkp_so(ih,jh,2,ityp(na))*CONJG(struc(na))* becp0%nc(jkb,2,ibnd)
                                   ps (ikb, (ibnd-1)*npol+2) = ps (ikb,(ibnd-1)*npol+2 ) + &
                                        q_dkp_so(ih,jh,3,ityp(na))*CONJG(struc(na))* becp0%nc(jkb,1,ibnd)
                                   ps (ikb, (ibnd-1)*npol+2) = ps (ikb,(ibnd-1)*npol+2 ) + &
                                        q_dkp_so(ih,jh,4,ityp(na))*CONJG(struc(na))* becp0%nc(jkb,2,ibnd)

                              else
                                 ps (ikb, ibnd) = ps (ikb, ibnd) + &
                                      q_dkp(ih,jh,ityp(na))*CONJG(struc(na))* becp0%k(jkb,ibnd)
                              endif
                           enddo
                        enddo
                     enddo
                     ijkb0 = ijkb0 + nh (nt)
                  endif
               enddo
            enddo

            call ZGEMM ('N', 'N', npw1, nbnd*npol , nkb, (1.d0, 0.d0) , vkb, &!vkb is relative to the last ik read
                 npwx, ps, nkb, (1.d0, 0.d0) , evct, npwx)
            do m=1,nbnd
               do nb=1,nbnd
                  do ig=1,npw1
                     evcm(ig,m,pdir)=evcm(ig,m,pdir) + mat(nb,m)*evct(ig,nb)
                  enddo
                  if(noncolin) then
                     do ig=1,npw1
                        evcm(ig+npwx,m,pdir)=evcm(ig+npwx,m,pdir) + mat(nb,m)*evct(ig+npwx,nb)
                     enddo
                  endif
               enddo
            enddo
         endif

 
!           --- End of dot products between wavefunctions and betas ---
      ELSE !(ik_stringa == 1)
      
         ikk = nx_el(ik+nppstr_3d(pdir)-1,pdir)
         npw0   = ngk(ikk)
         igk0(:)= igk_k(:,ikk)

         CALL get_buffer (evct,nwordwfc,iunwfc,nx_el(ik+nppstr_3d(pdir)-1,pdir))
!        

!           --- Calculate dot products between wavefunctions

!              --- Dot wavefunctions and betas for PREVIOUS k-point ---
          
         if(okvan) then
            CALL init_us_2 (npw0,igk0,xk(1,nx_el(ik+nppstr_3d(pdir)-1,pdir)),vkb)
            CALL calbec( npw0, vkb, evct, becp0 )
         endif
!              --- Dot wavefunctions and betas for CURRENT k-point ---
         
         ikk = nx_el(ik,pdir)
         npw1   = ngk(ikk)
         igk1(:)= igk_k(:,ikk)

         !  --- Recalculate FFT correspondence (see ggen.f90) ---

         if(.not.l_para) then
            ln0=0!set to 0
            DO ig=1,npw1
               mk1=nint(g(1,igk1(ig))*at(1,1)+g(2,igk1(ig))*at(2,1)+g(3,igk1(ig))*at(3,1))
               mk2=nint(g(1,igk1(ig))*at(1,2)+g(2,igk1(ig))*at(2,2)+g(3,igk1(ig))*at(3,2))
               mk3=nint(g(1,igk1(ig))*at(1,3)+g(2,igk1(ig))*at(2,3)+g(3,igk1(ig))*at(3,3))
               ln0(mk1,mk2,mk3) = ig
            END DO
         endif
         if(okvan) then
            CALL init_us_2 (npw1,igk1,xk(1,nx_el(ik,pdir)),vkb)
            CALL calbec( npw1, vkb, evcel, becp_bp )
         endif
!              --- Matrix elements calculation ---

        IF(mdone(nx_el(ik,pdir))==0) THEN
         mat=(0.d0,0.d0)
         if(.not. l_para) then
            map_g(:) = 0
            do ig=1,npw0
!                          --- If k'=k+G_o, the relation psi_k+G_o (G-G_o) ---
!                          --- = psi_k(G) is used, gpar=G_o, gtr = G-G_o ---
               gtr(1)=g(1,igk0(ig)) + gpar(1)
               gtr(2)=g(2,igk0(ig)) + gpar(2) 
               gtr(3)=g(3,igk0(ig)) + gpar(3) 
               !                          --- Find crystal coordinates of gtr, n1,n2,n3 ---
!                          --- and the position ng in the ngm array ---
               IF (gtr(1)**2+gtr(2)**2+gtr(3)**2 <= gcutm) THEN
                  n1=NINT(gtr(1)*at(1,1)+gtr(2)*at(2,1) &
                       +gtr(3)*at(3,1))
                  n2=NINT(gtr(1)*at(1,2)+gtr(2)*at(2,2) &
                       +gtr(3)*at(3,2))
                  n3=NINT(gtr(1)*at(1,3)+gtr(2)*at(2,3) &
                       +gtr(3)*at(3,3))
                  ng=ln(n1,n2,n3) 
                  IF ((ABS(g(1,ng)-gtr(1)) > eps) .OR. &
                       (ABS(g(2,ng)-gtr(2)) > eps) .OR. &
                       (ABS(g(3,ng)-gtr(3)) > eps)) THEN
                     WRITE(6,*) ' error hepsiher: translated G=', &
                          gtr(1),gtr(2),gtr(3), &
                          ' with crystal coordinates',n1,n2,n3, &
                          ' corresponds to ng=',ng,' but G(ng)=', &
                          g(1,ng),g(2,ng),g(3,ng)
                     WRITE(6,*) ' probably because G_par is NOT', &
                          ' a reciprocal lattice vector '
                     WRITE(6,*) ' Possible choices as smallest ', &
                          ' G_par:'
                     DO i=1,50
                        WRITE(6,*) ' i=',i,'   G=', &
                             g(1,i),g(2,i),g(3,i)
                     ENDDO
                     STOP
                  ENDIF
               ELSE 
                  WRITE(6,*) ' |gtr| > gcutm  for gtr=', &
                       gtr(1),gtr(2),gtr(3) 
                  STOP
               END IF
               map_g(ig)=ng
            enddo
         endif


         ! OPTIMIZATION BY AM
         ! NOTE THERE ARE TOO MANY COMMUNICATION CALLS FOR GLOBAL ARRAY
         ! CAN REDUCE THEM SIGNIFICANTLY  !!!
         ! NOTE CHANGED ORDER OF LOOPS OVER BANDS!!!
            DO mb=1,nbnd
              IF(l_para) THEN
                !allocate global array
                allocate(aux_g(ngm_g))
                aux_g=(0.d0,0.d0)
                IF (noncolin) THEN
                  allocate(aux_g_2(ngm_g))
                  aux_g_2=(0.d0,0.d0)
                END IF
!put psi1 on global array
                do ig=1,npw0
                  aux_g(mapgp_global(ig_l2g(igk0(ig)),pdir))=evct(ig,mb)
                  IF (noncolin) aux_g_2(mapgp_global(ig_l2g(igk0(ig)),pdir))=evct(ig+npwx,mb)
                enddo
                call mp_sum(aux_g(:),intra_bgrp_comm)
                IF (noncolin) call mp_sum(aux_g_2(:),intra_bgrp_comm)
              END IF

         DO nb=1,nbnd
               IF ( .NOT. l_cal(nb) .OR. .NOT. l_cal(mb) ) THEN
                  IF ( nb == mb )  mat(nb,mb)=0.d0
               ELSE
                  if(.not.l_para) then
                     aux=(0.d0,0.d0)
                     aux0=(0.d0,0.d0)
                     IF(noncolin) aux_2=(0.d0,0.d0)
                     IF(noncolin) aux0_2=(0.d0,0.d0)
                     DO ig=1,npw1
                        aux0(igk1(ig))=evcel(ig,nb)
                        IF(noncolin) aux0_2(igk1(ig))=evcel(ig+npwx,nb)
                     END DO
                     
                     do ig=1,npw0
                        aux(map_g(ig))=evct(ig,mb)
                        IF (noncolin) aux_2(map_g(ig))=evct(ig+npwx,mb)
                     ENDDO
              
                     mat(nb,mb) = zdotc(ngm,aux0,1,aux,1)
                     IF (noncolin) mat(nb,mb) = mat(nb,mb)+zdotc(ngm,aux0_2,1,aux_2,1)
                  else
                     sca=(0.d0,0.d0)
!do scalar product
                     do ig=1,npw1
                        sca=sca+conjg(evcel(ig,nb))*aux_g(ig_l2g(igk1(ig)))
                        IF (noncolin) sca=sca+conjg(evcel(ig+npwx,nb))*aux_g_2(ig_l2g(igk1(ig)))
                     enddo
! mp_sum is done later!!!
                     mat(nb,mb)=sca
                  endif
               endif
            END DO
                    IF(l_para) THEN
                     deallocate(aux_g)
                     IF (noncolin) deallocate(aux_g_2)
                    END IF
         END DO
         !
         call mp_sum( mat, intra_bgrp_comm )
         DO nb=1,nbnd
            IF(.not.l_cal(nb)) mat(nb,nb)=1.d0
         END DO

         !
         DO nb=1,nbnd
            DO mb=1,nbnd
               IF ( l_cal(nb) .AND. l_cal(mb) ) THEN
!                    --- Calculate the augmented part: ij=KB projectors, ---
!                    --- R=atom index: SUM_{ijR} q(ijR) <u_nk|beta_iR>   ---
!                    --- <beta_jR|u_mk'> e^i(k-k')*R =                   ---
!                    --- also <u_nk|beta_iR>=<psi_nk|beta_iR> = becp^*   ---
                  if (okvan) then 
                     pref = (0.d0,0.d0)
                     DO jkb=1,nkb
                        nhjkb = nkbtonh(jkb)
                        na = nkbtona(jkb)
                        np = ityp(na)
                        nhjkbm = nh(np)
                        jkb1 = jkb - nhjkb
                        DO j = 1,nhjkbm
                           if(lspinorb) then
                                pref = pref+CONJG(becp_bp%nc(jkb,1,nb))*becp0%nc(jkb1+j,1,mb) &
                                     *q_dkp_so(nhjkb,j,1,np)*CONJG(struc(na))
                                pref = pref+CONJG(becp_bp%nc(jkb,1,nb))*becp0%nc(jkb1+j,2,mb) &
                                     *q_dkp_so(nhjkb,j,2,np)*CONJG(struc(na))
                                pref = pref+CONJG(becp_bp%nc(jkb,2,nb))*becp0%nc(jkb1+j,1,mb) &
                                     *q_dkp_so(nhjkb,j,3,np)*CONJG(struc(na))
                                pref = pref+CONJG(becp_bp%nc(jkb,2,nb))*becp0%nc(jkb1+j,2,mb) &
                                     *q_dkp_so(nhjkb,j,4,np)*CONJG(struc(na))

                           else
                              pref = pref+CONJG(becp_bp%k(jkb,nb))*becp0%k(jkb1+j,mb) &
                                   *q_dkp(nhjkb,j,np)*CONJG(struc(na))
                           endif
                           ENDDO
                     ENDDO
                     mat(nb,mb) = mat(nb,mb) + pref
                  endif
               endif
            ENDDO
         ENDDO

!              --- Calculate matrix inverse ---
         CALL zgefa(mat,nbnd,nbnd,ivpt,info)
         CALL errore('h_epsi_her_set','error in zgefa',abs(info))
         CALL zgedi(mat,nbnd,nbnd,ivpt,cdet,cdwork,1)
         matbig(nx_el(ik,pdir),:,:)=mat
         mdone(nx_el(ik,pdir))=1
        ELSE
         mat=matbig(nx_el(ik,pdir),:,:)
        END IF
        
        DO nb = 1, nbnd
           DO mb = 1, nbnd
              if(.not.l_cal(nb).or. .not. l_cal(mb)) mat(mb,nb)=(0.d0,0.d0)
           END DO
        END DO
    
     
!    mat=S^-1(k,k-1)
        
         if(.not.l_para) then
            do ig=1,npw0
               gtr(1)=g(1,igk0(ig)) + gpar(1)
               gtr(2)=g(2,igk0(ig)) + gpar(2)        
               gtr(3)=g(3,igk0(ig)) + gpar(3) 
               
            !                          --- Find crystal coordinates of gtr, n1,n2,n3 ---
            !                          --- and the position ng in the ngm array ---
               IF (gtr(1)**2+gtr(2)**2+gtr(3)**2 <= gcutm) THEN
                  n1=NINT(gtr(1)*at(1,1)+gtr(2)*at(2,1) &
                       &      +gtr(3)*at(3,1))
                  n2=NINT(gtr(1)*at(1,2)+gtr(2)*at(2,2) &
                       &      +gtr(3)*at(3,2))
                  n3=NINT(gtr(1)*at(1,3)+gtr(2)*at(2,3) &
                       &      +gtr(3)*at(3,3))
                  ng=ln0(n1,n2,n3) 
                  if(ng .gt. 0) then
                     do m=1,nbnd
                        do nb=1,nbnd
                           evcm(ng,m,pdir)=evcm(ng,m,pdir) + mat(nb,m)*evct(ig,nb)
                           IF (noncolin) evcm(ng+npwx,m,pdir)=evcm(ng+npwx,m,pdir) &
                             +mat(nb,m)*evct(ig+npwx,nb)
                        enddo
                     enddo
                  endif
               ENDIF
            enddo
         else
!allocate
            allocate(aux_g(ngm_g))
            IF (noncolin) allocate(aux_g_2(ngm_g))
!loop on nb
            do nb=1,nbnd
               aux_g(:)=(0.d0,0.d0)
               IF (noncolin) aux_g_2(:)=(0.d0,0.d0)
               do ig=1,npw0
                  aux_g(mapgp_global(ig_l2g(igk0(ig)),pdir))=evct(ig,nb)
                  IF (noncolin) aux_g_2(mapgp_global(ig_l2g(igk0(ig)),pdir))=evct(ig+npwx,nb)
               enddo
!put evct on global  array
               call mp_sum(aux_g(:),intra_bgrp_comm)
               IF (noncolin) call mp_sum(aux_g_2(:),intra_bgrp_comm)
               do m=1,nbnd
                  do ig=1,npw1
                     evcm(ig,m,pdir)=evcm(ig,m,pdir)+mat(nb,m)*aux_g(ig_l2g(igk1(ig)))
                     IF (noncolin) evcm(ig+npwx,m,pdir)=evcm(ig+npwx,m,pdir) &
                       +mat(nb,m)*aux_g_2(ig_l2g(igk1(ig)))
                  enddo
               enddo
            enddo
            deallocate(aux_g)
            IF (noncolin) deallocate(aux_g_2)
         endif
         if(okvan) then
            evct(:,:) =  (0.d0, 0.d0)
            ps (:,:) = (0.d0, 0.d0)
            ijkb0 = 0
            do nt = 1, ntyp
               do na = 1, nat
                  if (ityp (na) .eq.nt) then
                     do ibnd = 1, nbnd
                        do jh = 1, nh (nt)
                           jkb = ijkb0 + jh
                           do ih = 1, nh (nt)
                              ikb = ijkb0 + ih
                              if(lspinorb) then
                                 ps (ikb, (ibnd-1)*npol+1) = ps (ikb,(ibnd-1)*npol+1 ) + &
                                      q_dkp_so(ih,jh,1,ityp(na))*CONJG(struc(na))* becp0%nc(jkb,1,ibnd)
                                 ps (ikb, (ibnd-1)*npol+1) = ps (ikb,(ibnd-1)*npol+1 ) + &
                                      q_dkp_so(ih,jh,2,ityp(na))*CONJG(struc(na))* becp0%nc(jkb,2,ibnd)
                                 ps (ikb, (ibnd-1)*npol+2) = ps (ikb,(ibnd-1)*npol+2 ) + &
                                      q_dkp_so(ih,jh,3,ityp(na))*CONJG(struc(na))* becp0%nc(jkb,1,ibnd)
                                 ps (ikb, (ibnd-1)*npol+2) = ps (ikb,(ibnd-1)*npol+2 ) + &
                                        q_dkp_so(ih,jh,4,ityp(na))*CONJG(struc(na))* becp0%nc(jkb,2,ibnd)
                              else
                                 ps (ikb, ibnd) = ps (ikb, ibnd) + &
                                      q_dkp(ih,jh,ityp(na))*CONJG(struc(na))* becp0%k(jkb,ibnd)
                              endif
                           enddo
                        enddo
                     enddo
                     ijkb0 = ijkb0 + nh (nt)
                  endif
               enddo
            enddo
            
            call ZGEMM ('N', 'N', npw1, nbnd*npol , nkb, (1.d0, 0.d0) , vkb, &!vkb is relative to the last ik read
                 npwx, ps, nkb, (1.d0, 0.d0) , evct, npwx)
            do m=1,nbnd
               do nb=1,nbnd
                  do ig=1,npw1
                     evcm(ig,m,pdir)=evcm(ig,m,pdir) + mat(nb,m)*evct(ig,nb)
                     if(noncolin) then
                        evcm(ig+npwx,m,pdir)=evcm(ig+npwx,m,pdir) + mat(nb,m)*evct(ig+npwx,nb)
                     endif
                  enddo
               enddo
            enddo
         endif

      ENDIF

! calculate  S-1(k,k+1)

    
!       
      if(ik_stringa /= nppstr_3d(pdir)) then
         ikk = nx_el(ik+1,pdir)
         npw0   = ngk(ikk)
         igk0(:)= igk_k(:,ikk)

         CALL get_buffer (evct,nwordwfc,iunwfc,nx_el(ik+1,pdir))
!        

!           --- Calculate dot products between wavefunctions

!              --- Dot wavefunctions and betas for PREVIOUS k-point ---
         
         if(okvan) then
            CALL init_us_2 (npw0,igk0,xk(1,nx_el(ik+1,pdir)),vkb)
            CALL calbec( npw0, vkb, evct, becp0)
         endif
!              --- Dot wavefunctions and betas for CURRENT k-point ---
         
         ikk = nx_el(ik,pdir)
         npw1   = ngk(ikk)
         igk1(:)= igk_k(:,ikk)

         !  --- Recalculate FFT correspondence (see ggen.f90) ---

         ln0=0!set to  0
         DO ig=1,npw1
            mk1=nint(g(1,igk1(ig))*at(1,1)+g(2,igk1(ig))*at(2,1)+g(3,igk1(ig))*at(3,1))
            mk2=nint(g(1,igk1(ig))*at(1,2)+g(2,igk1(ig))*at(2,2)+g(3,igk1(ig))*at(3,2))
            mk3=nint(g(1,igk1(ig))*at(1,3)+g(2,igk1(ig))*at(2,3)+g(3,igk1(ig))*at(3,3))
            ln0(mk1,mk2,mk3) = ig
         END DO
        
         if(okvan) then
            CALL init_us_2 (npw1,igk1,xk(1,nx_el(ik,pdir)),vkb)
            CALL calbec( npw1, vkb, evcel, becp_bp )
         endif
         



!              --- Matrix elements calculation ---

        IF(mdone(nx_el(ik+1,pdir))==0) THEN
         mat=(0.d0,0.d0)
         allocate(aux0vec(ngm,nbnd),aux1vec(ngm,nbnd))
         aux0vec=(0.d0,0.d0)
         aux1vec=(0.d0,0.d0)
         do nb=1,nbnd
            DO ig=1,npw1
               aux0vec(igk1(ig),nb)=evcel(ig,nb)
            END DO
         end do
         do nb=1,nbnd
            DO ig=1,npw0
               aux1vec(igk0(ig),nb)=evct(ig,nb)
            END DO            
         enddo
         call ZGEMM('C','N',nbnd,nbnd,ngm,(1.d0,0.d0),aux0vec,ngm,aux1vec,ngm,(0.d0,0.d0),mat,nbnd)
         if(noncolin) then
             aux0vec=(0.d0,0.d0)
             aux1vec=(0.d0,0.d0)
             do nb=1,nbnd
               DO ig=1,npw1
                  aux0vec(igk1(ig),nb)=evcel(ig+npwx,nb)
               END DO
            end do
            do nb=1,nbnd
               DO ig=1,npw0
                  aux1vec(igk0(ig),nb)=evct(ig+npwx,nb)
               END DO
            enddo
            call ZGEMM('C','N',nbnd,nbnd,ngm,(1.d0,0.d0),aux0vec,ngm,aux1vec,ngm,(1.d0,0.d0),mat,nbnd)
         endif
         deallocate(aux0vec,aux1vec)


      !
      call mp_sum( mat, intra_bgrp_comm )
      !
      DO nb=1,nbnd
         DO mb=1,nbnd
            IF(.not.l_cal(nb).or..not.l_cal(mb)) THEN
               IF(nb==mb) THEN
                  mat(nb,mb)=(1.d0,0.d0)
               ELSE
                  mat(nb,mb)=(0.d0,0.d0)
               END IF
            END IF
         END DO
      END DO

      DO nb=1,nbnd
         DO mb=1,nbnd
            IF ( l_cal(nb) .AND. l_cal(mb) ) THEN
               if(okvan) then
                  pref = (0.d0,0.d0)
                  DO jkb=1,nkb
                     nhjkb = nkbtonh(jkb)
                     na = nkbtona(jkb)
                     np = ityp(na)
                     nhjkbm = nh(np)
                     jkb1 = jkb - nhjkb
                     DO j = 1,nhjkbm
                        if(lspinorb) then
                            pref = pref+CONJG(becp_bp%nc(jkb,1,nb))*becp0%nc(jkb1+j,1,mb) &
                                 *q_dk_so(nhjkb,j,1,np)*struc(na)
                            pref = pref+CONJG(becp_bp%nc(jkb,1,nb))*becp0%nc(jkb1+j,2,mb) &
                                 *q_dk_so(nhjkb,j,2,np)*struc(na)
                            pref = pref+CONJG(becp_bp%nc(jkb,2,nb))*becp0%nc(jkb1+j,1,mb) &
                                 *q_dk_so(nhjkb,j,3,np)*struc(na)
                            pref = pref+CONJG(becp_bp%nc(jkb,2,nb))*becp0%nc(jkb1+j,2,mb) &
                                 *q_dk_so(nhjkb,j,4,np)*struc(na)


                        else
                           pref = pref+CONJG(becp_bp%k(jkb,nb))*becp0%k(jkb1+j,mb) &
                                *q_dk(nhjkb,j,np)*struc(na)
                        endif
                     ENDDO
                  ENDDO
                  mat(nb,mb) = mat(nb,mb) + pref
               endif
            ENDIF
         ENDDO
      ENDDO


!              --- Calculate matrix inverse ---
      CALL zgefa(mat,nbnd,nbnd,ivpt,info)
      CALL errore('h_epsi_her_set','error in zgefa',abs(info))
      CALL zgedi(mat,nbnd,nbnd,ivpt,cdet,cdwork,1)
      matbig(nx_el(ik+1,pdir),:,:)=TRANSPOSE(CONJG(mat))
      mdone(nx_el(ik+1,pdir))=1
   ELSE
      mat=TRANSPOSE(CONJG(matbig(nx_el(ik+1,pdir),:,:)))
   END IF
!in case of occupations from input set to 0 parts of mat

   DO nb = 1, nbnd
      DO mb = 1, nbnd
         if(.not.l_cal(nb).or. .not. l_cal(mb)) mat(mb,nb)=(0.d0,0.d0)
      END DO
   END DO

!    mat=S^-1(k,k-1)
      do ig=1,npw0
         gtr(1)=g(1,igk0(ig))
         gtr(2)=g(2,igk0(ig))         
         gtr(3)=g(3,igk0(ig))         
            !                          --- Find crystal coordinates of gtr, n1,n2,n3 ---
            !                          --- and the position ng in the ngm array ---
         IF (gtr(1)**2+gtr(2)**2+gtr(3)**2 <= gcutm) THEN
            n1=NINT(gtr(1)*at(1,1)+gtr(2)*at(2,1) &
             &       +gtr(3)*at(3,1))
            n2=NINT(gtr(1)*at(1,2)+gtr(2)*at(2,2) &
             &       +gtr(3)*at(3,2))
            n3=NINT(gtr(1)*at(1,3)+gtr(2)*at(2,3) &
             &       +gtr(3)*at(3,3))
            ng=ln0(n1,n2,n3)
            if(ng .gt. 0) then
               do m=1,nbnd
                  do nb=1,nbnd
                     evcp(ng,m,pdir)=evcp(ng,m,pdir) + mat(nb,m)*evct(ig,nb)
                     IF (noncolin) evcp(ng+npwx,m,pdir)=evcp(ng+npwx,m,pdir) &
                       + mat(nb,m)*evct(ig+npwx,nb)
                  enddo
               enddo
            endif
         ENDIF
      enddo
      if(okvan) then
         evct(:,:) =  (0.d0, 0.d0)
         ps (:,:) = (0.d0, 0.d0)
         ijkb0 = 0
         do nt = 1, ntyp
            do na = 1, nat
               if (ityp (na) .eq.nt) then
                  do ibnd = 1, nbnd
                     do jh = 1, nh (nt)
                        jkb = ijkb0 + jh
                        do ih = 1, nh (nt)
                           ikb = ijkb0 + ih
                           if(noncolin) then
                              ps (ikb, (ibnd-1)*npol+1) = ps (ikb,(ibnd-1)*npol+1 ) + &
                                   q_dk_so(ih,jh,1,ityp(na))*struc(na)* becp0%nc(jkb,1,ibnd)
                              ps (ikb, (ibnd-1)*npol+1) = ps (ikb,(ibnd-1)*npol+1 ) + &
                                   q_dk_so(ih,jh,2,ityp(na))*struc(na)* becp0%nc(jkb,2,ibnd)
                              ps (ikb, (ibnd-1)*npol+2) = ps (ikb,(ibnd-1)*npol+2 ) + &
                                   q_dk_so(ih,jh,3,ityp(na))*struc(na)* becp0%nc(jkb,1,ibnd)
                              ps (ikb, (ibnd-1)*npol+2) = ps (ikb,(ibnd-1)*npol+2 ) + &
                                   q_dk_so(ih,jh,4,ityp(na))*struc(na)* becp0%nc(jkb,2,ibnd)

                           else
                              ps (ikb, ibnd) = ps (ikb, ibnd) + &
                                   q_dk(ih,jh,ityp(na))*struc(na)* becp0%k(jkb,ibnd)
                           endif
                        enddo
                     enddo
                  enddo
                  ijkb0 = ijkb0 + nh (nt)
               endif
            enddo
         enddo

         call ZGEMM ('N', 'N', npw1, nbnd*npol , nkb, (1.d0, 0.d0) , vkb, &!vkb is relative to the last ik read
                 npwx, ps, nkb, (1.d0, 0.d0) , evct, npwx)
                 
         evct(npw1+1:npwx,1:nbnd)=(0.d0,0.d0)
         if(noncolin) evct(npwx+npw1+1:2*npwx,1:nbnd)=(0.d0,0.d0)
         evcp(npw1+1:npwx,1:nbnd,pdir)=(0.d0,0.d0)
         if(noncolin) evcp(npwx+npw1+1:2*npwx,1:nbnd,pdir)=(0.d0,0.d0)
         if(.not.noncolin) then
            call ZGEMM('N','N',npw1,nbnd,nbnd,(1.d0,0.d0),evct,npwx,mat,nbnd,&
                 &(1.d0,0.d0),evcp(1,1,pdir),npwx)
         else
            call ZGEMM('N','N',npwx*npol,nbnd,nbnd,(1.d0,0.d0),evct,npwx*npol,mat,nbnd,&
              &(1.d0,0.d0),evcp(1,1,pdir),npwx*npol)
         end if




      endif
         
!           --- End of dot products between wavefunctions and betas ---
     
   else
      
      ikk = nx_el(ik-nppstr_3d(pdir)+1,pdir)
      npw0   = ngk(ikk)
      igk0(:)= igk_k(:,ikk)

      CALL get_buffer (evct,nwordwfc,iunwfc,nx_el(ik-nppstr_3d(pdir)+1,pdir))
!        

!           --- Calculate dot products between wavefunctions

!              --- Dot wavefunctions and betas for PREVIOUS k-point ---

      if(okvan) then
         CALL init_us_2 (npw0,igk0,xk(1,nx_el(ik-nppstr_3d(pdir)+1,pdir)),vkb)
         CALL calbec( npw0, vkb, evct, becp0 )
      endif
!              --- Dot wavefunctions and betas for CURRENT k-point ---
         
      ikk = nx_el(ik,pdir)
      npw1   = ngk(ikk)
      igk1(:)= igk_k(:,ikk)

         !  --- Recalculate FFT correspondence (see ggen.f90) ---

      if(.not.l_para) then
         ln0=0! set to 0
         DO ig=1,npw1
            mk1=nint(g(1,igk1(ig))*at(1,1)+g(2,igk1(ig))*at(2,1)+g(3,igk1(ig))*at(3,1))
            mk2=nint(g(1,igk1(ig))*at(1,2)+g(2,igk1(ig))*at(2,2)+g(3,igk1(ig))*at(3,2))
            mk3=nint(g(1,igk1(ig))*at(1,3)+g(2,igk1(ig))*at(2,3)+g(3,igk1(ig))*at(3,3))
            ln0(mk1,mk2,mk3) = ig
         END DO
      endif

      if(okvan) then
         CALL init_us_2 (npw1,igk1,xk(1,nx_el(ik,pdir)),vkb)
         CALL calbec( npw1, vkb, evcel, becp_bp )
      endif
!              --- Matrix elements calculation ---

     IF(mdone(nx_el(ik-nppstr_3d(pdir)+1,pdir))==0) THEN
      if(.not.l_para) then
         map_g(:) = 0
         do ig=1,npw0
!                          --- If k'=k+G_o, the relation psi_k+G_o (G-G_o) ---
!                          --- = psi_k(G) is used, gpar=G_o, gtr = G-G_o ---
            gtr(1)=g(1,igk0(ig)) - gpar(1)
            gtr(2)=g(2,igk0(ig)) - gpar(2) 
            gtr(3)=g(3,igk0(ig)) - gpar(3) 
         !                          --- Find crystal coordinates of gtr, n1,n2,n3 ---
!                          --- and the position ng in the ngm array ---
            IF (gtr(1)**2+gtr(2)**2+gtr(3)**2 <= gcutm) THEN
               n1=NINT(gtr(1)*at(1,1)+gtr(2)*at(2,1) &
                    +gtr(3)*at(3,1))
               n2=NINT(gtr(1)*at(1,2)+gtr(2)*at(2,2) &
                    +gtr(3)*at(3,2))
               n3=NINT(gtr(1)*at(1,3)+gtr(2)*at(2,3) &
                    +gtr(3)*at(3,3))
               ng=ln(n1,n2,n3) 
               IF ((ABS(g(1,ng)-gtr(1)) > eps) .OR. &
                    (ABS(g(2,ng)-gtr(2)) > eps) .OR. &
                    (ABS(g(3,ng)-gtr(3)) > eps)) THEN
                  WRITE(6,*) ' error hepsiher: translated G=', &
                       gtr(1),gtr(2),gtr(3), &
                       ' with crystal coordinates',n1,n2,n3, &
                       ' corresponds to ng=',ng,' but G(ng)=', &
                       g(1,ng),g(2,ng),g(3,ng)
                  WRITE(6,*) ' probably because G_par is NOT', &
                       ' a reciprocal lattice vector '
                  WRITE(6,*) ' Possible choices as smallest ', &
                       ' G_par:'
                  DO i=1,50
                     WRITE(6,*) ' i=',i,'   G=', &
                          g(1,i),g(2,i),g(3,i)
                  ENDDO
                  STOP
               ENDIF
            ELSE 
               WRITE(6,*) ' |gtr| > gcutm  for gtr=', &
                    gtr(1),gtr(2),gtr(3) 
               STOP
            END IF
            map_g(ig)=ng 
         ENDDO
      endif


      mat=(0.d0,0.d0)
         DO mb=1,nbnd
               if(l_para) then
 !allocate global array
                  allocate(aux_g(ngm_g))
                  aux_g=(0.d0,0.d0)
                  IF (noncolin) THEN
                    allocate(aux_g_2(ngm_g))
                    aux_g_2=(0.d0,0.d0)
                  END IF
!put psi1 on global array
                  do ig=1,npw0
                     aux_g(mapgm_global(ig_l2g(igk0(ig)),pdir))=evct(ig,mb)
                     IF (noncolin) aux_g_2(mapgm_global(ig_l2g(igk0(ig)),pdir))=evct(ig+npwx,mb)
                  enddo
                  call mp_sum(aux_g(:),intra_bgrp_comm)
                  IF (noncolin) call mp_sum(aux_g_2(:),intra_bgrp_comm)
               end if

      DO nb=1,nbnd
            IF ( .NOT. l_cal(nb) .OR. .NOT. l_cal(mb) ) THEN
               IF ( nb == mb )  mat(nb,mb)=0.d0
            ELSE
               if(.not.l_para) then
                  aux=(0.d0,0.d0)
                  aux0=(0.d0,0.d0)
                  IF(noncolin) aux_2=(0.d0,0.d0)
                  IF(noncolin) aux0_2=(0.d0,0.d0)
                  DO ig=1,npw1
                     aux0(igk1(ig))=evcel(ig,nb)
                     IF(noncolin) aux0_2(igk1(ig))=evcel(ig+npwx,nb)
                  END DO
                  do ig=1,npw0
                     aux(map_g(ig))=evct(ig,mb)
                     IF(noncolin) aux_2(map_g(ig))=evct(ig+npwx,mb)
                  ENDDO
                  mat(nb,mb) = zdotc(ngm,aux0,1,aux,1)
                  IF (noncolin) mat(nb,mb)=mat(nb,mb)+zdotc(ngm,aux0_2,1,aux_2,1)
               else
                  sca=(0.d0,0.d0)
!do scalar product
                  do ig=1,npw1
                     sca=sca+conjg(evcel(ig,nb))*aux_g(ig_l2g(igk1(ig)))
                     IF (noncolin) sca=sca+conjg(evcel(ig+npwx,nb)) &
                       *aux_g_2(ig_l2g(igk1(ig)))
                  enddo
! mp_sum is done later!!!
                  mat(nb,mb)=sca
                  
               endif
            endif
         END DO
                IF(l_para) THEN
                  deallocate(aux_g)
                  IF (noncolin) deallocate(aux_g_2)
                END IF
      END DO
      ! 
      call mp_sum( mat, intra_bgrp_comm )
      DO nb=1,nbnd
         IF(.not.l_cal(nb)) mat(nb,nb)=1.d0
      END DO

      !
      DO nb=1,nbnd
         DO mb=1,nbnd
            IF ( l_cal(nb) .AND. l_cal(mb) ) THEN
               if(okvan) then 
                  pref = (0.d0,0.d0)
                  DO jkb=1,nkb
                     nhjkb = nkbtonh(jkb)
                     na = nkbtona(jkb)
                     np = ityp(na)
                     nhjkbm = nh(np)
                     jkb1 = jkb - nhjkb
                     DO j = 1,nhjkbm
                        if(lspinorb) then
                           pref = pref+CONJG(becp_bp%nc(jkb,1,nb))*becp0%nc(jkb1+j,1,mb) &
                                *q_dk_so(nhjkb,j,1,np)*struc(na)
                           pref = pref+CONJG(becp_bp%nc(jkb,1,nb))*becp0%nc(jkb1+j,2,mb) &
                                *q_dk_so(nhjkb,j,2,np)*struc(na)
                           pref = pref+CONJG(becp_bp%nc(jkb,2,nb))*becp0%nc(jkb1+j,1,mb) &
                                *q_dk_so(nhjkb,j,3,np)*struc(na)
                           pref = pref+CONJG(becp_bp%nc(jkb,2,nb))*becp0%nc(jkb1+j,2,mb) &
                                *q_dk_so(nhjkb,j,4,np)*struc(na)
                        else
                           pref = pref+CONJG(becp_bp%k(jkb,nb))*becp0%k(jkb1+j,mb) &
                                *q_dk(nhjkb,j,np)*struc(na)
                        endif
                     ENDDO
                  ENDDO
                  mat(nb,mb) = mat(nb,mb) + pref
               endif
            ENDIF
         ENDDO
      ENDDO

!              --- Calculate matrix inverse ---
      CALL zgefa(mat,nbnd,nbnd,ivpt,info)
      CALL errore('h_epsi_her_set','error in zgefa',abs(info))
      CALL zgedi(mat,nbnd,nbnd,ivpt,cdet,cdwork,1)
      matbig(nx_el(ik-nppstr_3d(pdir)+1,pdir),:,:)=TRANSPOSE(CONJG(mat))
      mdone(nx_el(ik-nppstr_3d(pdir)+1,pdir))=1
     ELSE
      mat=TRANSPOSE(CONJG(matbig(nx_el(ik-nppstr_3d(pdir)+1,pdir),:,:)))
   END IF


      DO nb = 1, nbnd
         DO mb = 1, nbnd
            if(.not.l_cal(nb).or. .not. l_cal(mb)) mat(mb,nb)=(0.d0,0.d0)
         END DO
      END DO

       
!    mat=S^-1(k,k-1)
      if(.not.l_para) then
         do ig=1,npw0
            gtr(1)=g(1,igk0(ig)) - gpar(1)
            gtr(2)=g(2,igk0(ig)) - gpar(2)        
            gtr(3)=g(3,igk0(ig)) - gpar(3)        
            !                          --- Find crystal coordinates of gtr, n1,n2,n3 ---
            !                          --- and the position ng in the ngm array ---
            IF (gtr(1)**2+gtr(2)**2+gtr(3)**2 <= gcutm) THEN
               n1=NINT(gtr(1)*at(1,1)+gtr(2)*at(2,1) &
                    &    +gtr(3)*at(3,1))
               n2=NINT(gtr(1)*at(1,2)+gtr(2)*at(2,2) &
                    &    +gtr(3)*at(3,2))
               n3=NINT(gtr(1)*at(1,3)+gtr(2)*at(2,3) &
                    &    +gtr(3)*at(3,3))
               ng=ln0(n1,n2,n3)
               if(ng .gt. 0) then
                  do m=1,nbnd
                     do nb=1,nbnd
                        evcp(ng,m,pdir)=evcp(ng,m,pdir) + mat(nb,m)*evct(ig,nb)
                        IF (noncolin) evcp(ng+npwx,m,pdir)=evcp(ng+npwx,m,pdir) &
                          +mat(nb,m)*evct(ig+npwx,nb)
                     end do
                  enddo
               end if
            ENDIF
         enddo
      else

!allocate
         allocate(aux_g(ngm_g))
         IF (noncolin) allocate(aux_g_2(ngm_g))
!loop on nb
         do nb=1,nbnd
            aux_g(:)=(0.d0,0.d0)
            IF (noncolin) aux_g_2(:)=(0.d0,0.d0)
            do ig=1,npw0
               aux_g(mapgm_global(ig_l2g(igk0(ig)),pdir))=evct(ig,nb)
               IF (noncolin) aux_g_2(mapgm_global(ig_l2g(igk0(ig)),pdir))=evct(ig+npwx,nb)
            enddo
!put evct on global  array
            call mp_sum(aux_g(:),intra_bgrp_comm)
            IF (noncolin) call mp_sum(aux_g_2(:),intra_bgrp_comm)
            do m=1,nbnd
               do ig=1,npw1
                  evcp(ig,m,pdir)=evcp(ig,m,pdir)+mat(nb,m)*aux_g(ig_l2g(igk1(ig)))
                  IF (noncolin) evcp(ig+npwx,m,pdir)=evcp(ig+npwx,m,pdir) &
                    +mat(nb,m)*aux_g_2(ig_l2g(igk1(ig)))
               enddo
            enddo
         enddo
          deallocate(aux_g)
          IF (noncolin) deallocate(aux_g_2)
       endif
      if(okvan) then
         evct(:,:) =  (0.d0, 0.d0)
         ps (:,:) = (0.d0, 0.d0)
         ijkb0 = 0
         do nt = 1, ntyp
            do na = 1, nat
               if (ityp (na) .eq.nt) then
                  do ibnd = 1, nbnd
                     do jh = 1, nh (nt)
                        jkb = ijkb0 + jh
                        do ih = 1, nh (nt)
                           ikb = ijkb0 + ih
                           if(lspinorb) then
                              ps (ikb, (ibnd-1)*npol+1) = ps (ikb,(ibnd-1)*npol+1 ) + &
                                   q_dk_so(ih,jh,1,ityp(na))*struc(na)* becp0%nc(jkb,1,ibnd)
                              ps (ikb, (ibnd-1)*npol+1) = ps (ikb,(ibnd-1)*npol+1 ) + &
                                   q_dk_so(ih,jh,2,ityp(na))*struc(na)* becp0%nc(jkb,2,ibnd)
                              ps (ikb, (ibnd-1)*npol+2) = ps (ikb,(ibnd-1)*npol+2 ) + &
                                   q_dk_so(ih,jh,3,ityp(na))*struc(na)* becp0%nc(jkb,1,ibnd)
                              ps (ikb, (ibnd-1)*npol+2) = ps (ikb,(ibnd-1)*npol+2 ) + &
                                   q_dk_so(ih,jh,4,ityp(na))*struc(na)* becp0%nc(jkb,2,ibnd)

                           else
                              ps (ikb, ibnd) = ps (ikb, ibnd) + &
                                   q_dk(ih,jh,ityp(na))*struc(na)* becp0%k(jkb,ibnd)
                           endif
                        enddo
                     enddo
                  enddo
                  ijkb0 = ijkb0 + nh (nt)
               endif
            enddo
         enddo

         call ZGEMM ('N', 'N', npw1, nbnd*npol , nkb, (1.d0, 0.d0) , vkb, &!vkb is relative to the ik read
                 npwx, ps, nkb, (1.d0, 0.d0) , evct, npwx)
         do m=1,nbnd
            do nb=1,nbnd
               do ig=1,npw1
                  evcp(ig,m,pdir)=evcp(ig,m,pdir) + mat(nb,m)*evct(ig,nb)
                  if(noncolin) then
                     evcp(ig+npwx,m,pdir)=evcp(ig+npwx,m,pdir) + mat(nb,m)*evct(ig+npwx,nb)
                  endif
               enddo
            enddo
         enddo
      endif

   ENDIF
   


!writes projectors to disk 
   call save_buffer (evcm(:,:,pdir), nwordwfc, iunefieldm, nx_el(ik,pdir)+(pdir-1)*nks)
   call save_buffer (evcp(:,:,pdir), nwordwfc, iunefieldp, nx_el(ik,pdir)+(pdir-1)*nks)

  END  DO !on ik


  DEALLOCATE (l_cal)
  DEALLOCATE( evct)
  DEALLOCATE( map_g)
  deallocate(ln,ln0)
  DEALLOCATE(aux,aux0)
  IF (ALLOCATED(aux_2)) DEALLOCATE(aux_2)
  IF (ALLOCATED(aux0_2)) DEALLOCATE(aux0_2)
  if(okvan)CALL deallocate_bec_type (becp0)
  if(okvan)CALL deallocate_bec_type (becp_bp)
  if(okvan.and.lspinorb) deallocate(q_dk_so)
  if(okvan.and.lspinorb) deallocate(q_dkp_so)
  deallocate(mat)

  call stop_clock('h_epsi_set')
!  --
!------------------------------------------------------------------------------!
   return
 END SUBROUTINE h_epsi_her_set

!==============================================================================!

 SUBROUTINE factor_a(dir, a,fact)

   USE kinds, ONLY : DP

   IMPLICIT NONE
   
   REAL(kind=DP):: a(3,3),fact
   INTEGER :: dir

   INTEGER :: d1,d2
   REAL(kind=DP) :: v(3), sca
   
   if(dir==1) then
      d1=2
      d2=3
   else if(dir==2) then
      d1=3
      d2=1
   else if(dir==3) then
      d1=1
      d2=2
   endif
      
   
   !calculate vect(a(d1,:) X a(d2,:)
   
   v(1)=a(2,d1)*a(3,d2)-a(3,d1)*a(2,d2)
   v(2)=-a(1,d1)*a(3,d2)+a(3,d1)*a(1,d2)
   v(3)=a(1,d1)*a(2,d2)-a(2,d1)*a(1,d2)
 


  !normalize v

   sca=sqrt(v(1)**2.d0+v(2)**2.d0+v(3)**2.d0)
   v(:)=v(:)/sca

   !calculate a(dir:)*v(:)
 fact=v(1)*a(1,dir)+v(2)*a(2,dir)+v(3)*a(3,dir)
!!!!!!!!!!!!!!
 fact=dsqrt(a(1,dir)**2.d0+a(2,dir)**2.d0+a(3,dir)**2.d0)





 fact=abs(fact)



   return
 END SUBROUTINE factor_a
