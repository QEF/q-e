!
! Copyright (C) 2001-2012 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

SUBROUTINE forces_ion_efield (forces_bp, pdir, e_field)

!calculate ionic contribution , which is in the 
!a_gdir direction

  USE kinds,                ONLY : dp
  USE cell_base,            ONLY : at
  USE ions_base,            ONLY : nat,zv, ityp

  implicit none

  INTEGER, INTENT(in) :: pdir!direction on which the polarization is calculated
  REAL(DP), INTENT(in) :: e_field!intensity of the field
  REAL(DP), INTENT(inout) :: forces_bp(3,nat)!atomic forces to be update

  INTEGER i
  REAL(DP) :: e!electronic charge (Ry. a.u.)
  REAL(DP) :: a(3),sca

  e=dsqrt(2.d0)
!  a(:)=at(:,pdir)
!  sca=dsqrt(a(1)**2.d0 + a(2)**2.d0 + a(3)**2.d0)
!  a(:)=a(:)/sca

  do i=1,nat
!     forces_bp(:,i)=forces_bp(:,i)+ e*e_field*zv(ityp(i))*a(:) !ATTENZIONE
      forces_bp(pdir,i)=forces_bp(pdir,i)+ e*e_field*zv(ityp(i))
  enddo
  
  return


END SUBROUTINE forces_ion_efield



SUBROUTINE forces_us_efield(forces_bp, pdir, e_field)

!----------------------------------------------------------------------!

!it calculates the US correction to the atomic forces 
!due to Berry's phase electric field


!  --- Make use of the module with common information ---
   USE kinds,                ONLY : DP
   USE io_global,            ONLY : stdout
   USE io_files,             ONLY : iunwfc, nwordwfc
   USE buffers,              ONLY : get_buffer
   USE ions_base,            ONLY : nat, ntyp => nsp, ityp, tau, zv, atm
   USE cell_base,            ONLY : at, alat, tpiba, omega, tpiba2
   USE constants,            ONLY : pi, tpi
   USE gvect,                ONLY : ngm,  g, gcutm, ngm_g
   USE fft_base,             ONLY : dfftp
   USE uspp,                 ONLY : nkb, vkb, okvan
   USE uspp_param,           ONLY : upf, lmaxq, nbetam, nh, nhm
   USE lsda_mod,             ONLY : nspin
   USE klist,                ONLY : nelec, degauss, nks, xk, wk
   USE wvfct,                ONLY : npwx, npw, nbnd, ecutwfc
   USE wavefunctions_module, ONLY : evc
   USE bp,                   ONLY : nppstr_3d, mapgm_global, nx_el
   USE fixed_occ
   USE gvect,   ONLY : ig_l2g
   USE mp,                   ONLY : mp_sum
   USE mp_global,            ONLY : intra_bgrp_comm
   USE becmod,               ONLY : calbec
!  --- Avoid implicit definitions ---
   IMPLICIT NONE

   REAL(DP), INTENT(inout) :: forces_bp(3,nat)!atomic forces to be update 
   INTEGER, INTENT(in) :: pdir!direction of electric field
   REAL(DP), INTENT(in) :: e_field!initensity of the field

!  --- Internal definitions ---
   INTEGER :: i
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
   REAL(dp) :: dk(3)
   REAL(dp) :: dkmod
   REAL(dp) :: el_loc
   REAL(dp) :: eps
   REAL(dp) :: fac
   REAL(dp) :: g2kin_bp(npwx)
   REAL(dp) :: gpar(3)
   REAL(dp) :: gtr(3)
   REAL(dp) :: gvec
   REAL(dp), ALLOCATABLE :: loc_k(:)
   REAL(dp), ALLOCATABLE :: pdl_elec(:)
   REAL(dp), ALLOCATABLE :: phik(:)
   REAL(dp) :: qrad_dk(nbetam,nbetam,lmaxq,ntyp)
   REAL(dp) :: weight
   REAL(dp) :: pola, pola_ion
   REAL(dp), ALLOCATABLE :: wstring(:)
   REAL(dp) :: ylm_dk(lmaxq*lmaxq)
   REAL(dp) :: zeta_mod
   COMPLEX(dp), ALLOCATABLE :: aux(:)
   COMPLEX(dp), ALLOCATABLE :: aux0(:)
   COMPLEX(dp) :: becp0(nkb,nbnd)
   COMPLEX(dp) :: becp_bp(nkb,nbnd)
   COMPLEX(dp) , ALLOCATABLE :: cphik(:)
   COMPLEX(dp) :: det
   COMPLEX(dp), ALLOCATABLE :: mat(:,:)
   COMPLEX(dp) :: cdet(2)
   COMPLEX(dp) :: cdwork(nbnd)

   COMPLEX(dp) :: pref
   COMPLEX(dp) :: q_dk(nhm,nhm,ntyp)
   COMPLEX(dp) :: struc(nat),struc_r(3,nat)
   COMPLEX(dp) :: zdotc
   COMPLEX(dp) :: zeta

   COMPLEX(dp), ALLOCATABLE :: psi(:,:)
   COMPLEX(dp), ALLOCATABLE :: psi1(:,:)
   COMPLEX(dp) :: zeta_loc

   LOGICAL, ALLOCATABLE :: l_cal(:) ! flag for occupied/empty bands
   INTEGER, ALLOCATABLE :: map_g(:)

   REAL(dp) :: dkfact
   COMPLEX(dp) :: zeta_tot

   LOGICAL :: l_para! if true new parallel treatment
   COMPLEX(kind=DP) :: sca
   COMPLEX(kind=DP), ALLOCATABLE :: aux_g(:)
   COMPLEX(DP), ALLOCATABLE :: dbecp0(:,:,:), dbecp_bp(:,:,:),vkb1(:,:)
   INTEGER :: ipol
   COMPLEX(DP) :: forces_tmp(3,nat)
   REAL(DP) :: fact

!  -------------------------------------------------------------------------   !
!                               INITIALIZATIONS
!  -------------------------------------------------------------------------   !

   ALLOCATE (psi1(npwx,nbnd))
   ALLOCATE (psi(npwx,nbnd))
   ALLOCATE (aux(ngm))
   ALLOCATE (aux0(ngm))
   ALLOCATE (map_g(npwx))
   ALLOCATE (mat(nbnd,nbnd))
   ALLOCATE (dbecp0( nkb, nbnd, 3 ) ,dbecp_bp( nkb, nbnd, 3 ))
   ALLOCATE( vkb1( npwx, nkb ) )
   ALLOCATE( l_cal(nbnd) )

   if(pdir==3) then
      l_para=.false.
   else
      l_para=.true.
   endif


   pola=0.d0 !set to 0 electronic polarization   
   zeta_tot=(1.d0,0.d0)

!  --- Check that we are working with an insulator with no empty bands ---
   IF ( degauss > 0.0_dp ) CALL errore('forces_us_efield', &
            'Polarization only for insulators and no empty bands',1)

   !  --- Define a small number ---
   eps=1.0E-6_dp

!  --- Recalculate FFT correspondence (see ggen.f90) ---
   ALLOCATE (ln (-dfftp%nr1:dfftp%nr1, -dfftp%nr2:dfftp%nr2, -dfftp%nr3:dfftp%nr3) )
   DO ng=1,ngm
      mk1=nint(g(1,ng)*at(1,1)+g(2,ng)*at(2,1)+g(3,ng)*at(3,1))
      mk2=nint(g(1,ng)*at(1,2)+g(2,ng)*at(2,2)+g(3,ng)*at(3,2))
      mk3=nint(g(1,ng)*at(1,3)+g(2,ng)*at(2,3)+g(3,ng)*at(3,3))
      ln(mk1,mk2,mk3) = ng
   END DO

   if (okvan) then
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
!  --- Get the number of strings ---
   nstring=nks/nppstr_3d(pdir)
   nkort=nstring/(nspin)

!  --- Allocate memory for arrays ---
   ALLOCATE(phik(nstring))
   ALLOCATE(loc_k(nstring))
   ALLOCATE(cphik(nstring))
   ALLOCATE(wstring(nstring))
   ALLOCATE(pdl_elec(nstring))
   ALLOCATE(mod_elec(nstring))

!  -------------------------------------------------------------------------   !
!           electronic polarization: set values for k-points strings           !
!  -------------------------------------------------------------------------   !

!  --- Find vector along strings ---
   if(nppstr_3d(pdir) .ne. 1) then
      gpar(1)=(xk(1,nx_el(nppstr_3d(pdir),pdir))-xk(1,nx_el(1,pdir)))*&
           &DBLE(nppstr_3d(pdir))/DBLE(nppstr_3d(pdir)-1)
      gpar(2)=(xk(2,nx_el(nppstr_3d(pdir),pdir))-xk(2,nx_el(1,pdir)))*&
           &DBLE(nppstr_3d(pdir))/DBLE(nppstr_3d(pdir)-1)
      gpar(3)=(xk(3,nx_el(nppstr_3d(pdir),pdir))-xk(3,nx_el(1,pdir)))*&
           &DBLE(nppstr_3d(pdir))/DBLE(nppstr_3d(pdir)-1)
      gvec=dsqrt(gpar(1)**2+gpar(2)**2+gpar(3)**2)*tpiba
   else
      gpar(1)=0.d0
      gpar(2)=0.d0
      gpar(3)=0.d0
      gpar(pdir)=1.d0/at(pdir,pdir)!
      gvec=tpiba/sqrt(at(pdir,1)**2.d0+at(pdir,2)**2.d0+at(pdir,3)**2.d0)
   endif      
      
!  --- Find vector between consecutive points in strings ---
   if(nppstr_3d(pdir).ne.1) then  ! orthorhombic cell 
      dk(1)=xk(1,nx_el(2,pdir))-xk(1,nx_el(1,pdir))
      dk(2)=xk(2,nx_el(2,pdir))-xk(2,nx_el(1,pdir))
      dk(3)=xk(3,nx_el(2,pdir))-xk(3,nx_el(1,pdir))
      dkmod=SQRT(dk(1)**2+dk(2)**2+dk(3)**2)*tpiba 
   else ! Gamma point case, only cubic cell for now
      dk(1)=0.d0
      dk(2)=0.d0
      dk(3)=0.d0
      dk(pdir)=1.d0/at(pdir,pdir)
      dkmod=tpiba/sqrt(at(pdir,1)**2.d0+at(pdir,2)**2.d0+at(pdir,3)**2.d0)
   endif

!  -------------------------------------------------------------------------   !
!                   electronic polarization: weight strings                    !
!  -------------------------------------------------------------------------   !

!  --- Calculate string weights, normalizing to 1 (no spin) or 1+1 (spin) ---
   DO is=1,nspin
      weight=0.0_dp
      DO kort=1,nkort
         istring=kort+(is-1)*nkort
         wstring(istring)=wk(nppstr_3d(pdir)*istring)
         weight=weight+wstring(istring)
      END DO
      DO kort=1,nkort
         istring=kort+(is-1)*nkort
         wstring(istring)=wstring(istring)/weight
      END DO
   END DO  
  
!  -------------------------------------------------------------------------   !
!                  electronic polarization: structure factor                   !
!  -------------------------------------------------------------------------   !
   
!  --- Calculate structure factor e^{-i dk*R} ---

   DO na=1,nat
      fac=(dk(1)*tau(1,na)+dk(2)*tau(2,na)+dk(3)*tau(3,na))*tpi 
      struc(na)=CMPLX(cos(fac),-sin(fac),kind=DP)
   END DO

! Calculate derivatives of structure factors
   do na=1,nat
      do ipol=1,3
         struc_r(ipol,na)=struc(na)*CMPLX(0.d0,-1.d0, kind=dp)*dk(ipol)
      enddo
   enddo

!  -------------------------------------------------------------------------   !
!                     electronic polarization: form factor                     !
!  -------------------------------------------------------------------------   !
   if(okvan) then
!  --- Calculate Bessel transform of Q_ij(|r|) at dk [Q_ij^L(|r|)] ---
      CALL calc_btq(dkmod,qrad_dk,0)

!  --- Calculate the q-space real spherical harmonics at dk [Y_LM] --- 
      dkmod = dk(1)**2+dk(2)**2+dk(3)**2
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
   endif

!calculate factor

   
   call factor_a(pdir,at,dkfact)
   fact=dsqrt(2.d0)*e_field*dkfact
   if(nspin==1) fact=fact*2.d0


   
!  -------------------------------------------------------------------------   !
!                   electronic polarization: strings phases                    !
!  -------------------------------------------------------------------------   !

   el_loc=0.d0
   kpoint=0
   zeta=(1.d0,0.d0)
!  --- Start loop over spin ---
   DO is=1,nspin 
      ! l_cal(n) = .true./.false. if n-th state is occupied/empty
      DO nb = 1, nbnd
         IF ( nspin == 2 .AND. tfixed_occ) THEN
            l_cal(nb) = ( f_inp(nb,is) /= 0.0_dp )
         ELSE
            l_cal(nb) = ( nb <= NINT ( nelec/2.0_dp ) )
         ENDIF
      END DO
!     --- Start loop over orthogonal k-points ---
      DO kort=1,nkort
         zeta_loc=(1.d0,0.d0)
!        --- Index for this string ---
         istring=kort+(is-1)*nkort

!        --- Initialize expectation value of the phase operator ---
      
         zeta_mod = 1.d0

!        --- Start loop over parallel k-points ---
         DO kpar = 1,nppstr_3d(pdir)+1

!           --- Set index of k-point ---
            kpoint = kpoint + 1

!           --- Calculate dot products between wavefunctions and betas ---
            IF (kpar /= 1 ) THEN
               
!              --- Dot wavefunctions and betas for PREVIOUS k-point ---
               CALL gk_sort(xk(1,nx_el(kpoint-1,pdir)),ngm,g,ecutwfc/tpiba2, &
                            npw0,igk0,g2kin_bp) 
               CALL get_buffer (psi,nwordwfc,iunwfc,nx_el(kpoint-1,pdir))
               if (okvan) then
                  CALL init_us_2 (npw0,igk0,xk(1,nx_el(kpoint-1,pdir)),vkb)
                  CALL calbec( npw0, vkb, psi, becp0)
                  DO ipol = 1, 3
                     DO jkb = 1, nkb
                        DO ig = 1, npw0
                           vkb1(ig,jkb) = vkb(ig,jkb)*(0.D0,-1.D0)*g(ipol,igk0(ig))
                        END DO
                     END DO
                     IF ( nkb > 0 ) &
                          CALL ZGEMM( 'C', 'N', nkb, nbnd, npw0, ( 1.D0, 0.D0 ),   &
                          vkb1, npwx, psi, npwx, ( 0.D0, 0.D0 ),      &
                          dbecp0(1,1,ipol), nkb )
                  
                  ENDDO
               endif
!              --- Dot wavefunctions and betas for CURRENT k-point ---
               IF (kpar /= (nppstr_3d(pdir)+1)) THEN
                  CALL gk_sort(xk(1,nx_el(kpoint,pdir)),ngm,g,ecutwfc/tpiba2, &
                               npw1,igk1,g2kin_bp)        
                  CALL get_buffer (psi1,nwordwfc,iunwfc,nx_el(kpoint,pdir))
                  if(okvan) then
                     CALL init_us_2 (npw1,igk1,xk(1,nx_el(kpoint,pdir)),vkb)
                     CALL calbec( npw1, vkb, psi1, becp_bp)
                     DO ipol = 1, 3
                        DO jkb = 1, nkb
                           DO ig = 1, npw1
                              vkb1(ig,jkb) = vkb(ig,jkb)*(0.D0,-1.D0)*g(ipol,igk1(ig))
                           END DO
                        END DO
                        IF ( nkb > 0 ) &
                             CALL ZGEMM( 'C', 'N', nkb, nbnd, npw1, ( 1.D0, 0.D0 ),   &
                             vkb1, npwx, psi1, npwx, ( 0.D0, 0.D0 ),      &
                             dbecp_bp(1,1,ipol), nkb )
                  
                     ENDDO
                  endif
               ELSE
                  kstart = kpoint-(nppstr_3d(pdir)+1)+1
                  CALL gk_sort(xk(1,nx_el(kstart,pdir)),ngm,g,ecutwfc/tpiba2, &
                               npw1,igk1,g2kin_bp)  
                  CALL get_buffer (psi1,nwordwfc,iunwfc,nx_el(kstart,pdir))
                  if(okvan) then
                     CALL init_us_2 (npw1,igk1,xk(1,nx_el(kstart,pdir)),vkb)
                     CALL calbec( npw1, vkb, psi1, becp_bp)
                     DO ipol = 1, 3
                        DO jkb = 1, nkb
                           DO ig = 1, npw1
                              vkb1(ig,jkb) = vkb(ig,jkb)*(0.D0,-1.D0)*g(ipol,igk1(ig))
                           END DO
                        END DO
                        IF ( nkb > 0 ) &
                             CALL ZGEMM( 'C', 'N', nkb, nbnd, npw1, ( 1.D0, 0.D0 ),   &
                             vkb1, npwx, psi1, npwx, ( 0.D0, 0.D0 ),      &
                             dbecp_bp(1,1,ipol), nkb )
                  
                     ENDDO
                  endif
               ENDIF
           
!              --- Matrix elements calculation ---

               IF (kpar == (nppstr_3d(pdir)+1) .and. .not.l_para) THEN
                  map_g(:) = 0
                  do ig=1,npw1
!                          --- If k'=k+G_o, the relation psi_k+G_o (G-G_o) ---
!                          --- = psi_k(G) is used, gpar=G_o, gtr = G-G_o ---
                           
                     gtr(1)=g(1,igk1(ig)) - gpar(1)
                     gtr(2)=g(2,igk1(ig)) - gpar(2) 
                     gtr(3)=g(3,igk1(ig)) - gpar(3) 
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
                        IF ( (ABS(g(1,ng)-gtr(1)) > eps) .OR. &
                             (ABS(g(2,ng)-gtr(2)) > eps) .OR. &
                             (ABS(g(3,ng)-gtr(3)) > eps) ) THEN
                           WRITE(6,*) ' error: translated G=', &
                                gtr(1),gtr(2),gtr(3), &
                                &     ' with crystal coordinates',n1,n2,n3, &
                                &     ' corresponds to ng=',ng,' but G(ng)=', &
                                &     g(1,ng),g(2,ng),g(3,ng)
                           WRITE(6,*) ' probably because G_par is NOT', &
                                &    ' a reciprocal lattice vector '
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
               ENDIF

               mat=(0.d0,0.d0)
               DO nb=1,nbnd
                  DO mb=1,nbnd
                     IF ( .NOT. l_cal(nb) .OR. .NOT. l_cal(mb) ) THEN
                        IF ( nb == mb )  mat(nb,mb)=1.d0
                     ELSE
                        aux=(0.d0,0.d0)
                        aux0=(0.d0,0.d0)
                        DO ig=1,npw0
                           aux0(igk0(ig))=psi(ig,nb)
                        END DO
                        IF (kpar /= (nppstr_3d(pdir)+1)) THEN
                           do ig=1,npw1
                              aux(igk1(ig))=psi1(ig,mb)
                           enddo
                        ELSE IF( .not. l_para) THEN
                           do ig=1,npw1
                              aux(map_g(ig))=psi1(ig,mb)
                           enddo
                        ELSE
! allocate global array
                           allocate(aux_g(ngm_g))
                           aux_g=(0.d0,0.d0)
! put psi1 on global array
                           do ig=1,npw1
                              aux_g(mapgm_global(ig_l2g(igk1(ig)),pdir))=psi1(ig,mb)
                           enddo
                           call mp_sum(aux_g(:))
                           sca=(0.d0,0.d0)
! do scalar product
                           do ig=1,ngm
                              sca=sca+conjg(aux0(ig))*aux_g(ig_l2g(ig))
                           enddo
! do mp_sum
                           call mp_sum(sca)
                           mat(nb,mb)=sca
                           deallocate(aux_g)
                        ENDIF
                        
                        if(kpar /= (nppstr_3d(pdir)+1).or..not. l_para) then
                           mat(nb,mb) = zdotc(ngm,aux0,1,aux,1)                           
                           call mp_sum( mat(nb,mb), intra_bgrp_comm )
                        endif
!                    --- Calculate the augmented part: ij=KB projectors, ---
!                    --- R=atom index: SUM_{ijR} q(ijR) <u_nk|beta_iR>   ---
!                    --- <beta_jR|u_mk'> e^i(k-k')*R =                   ---
!                    --- also <u_nk|beta_iR>=<psi_nk|beta_iR> = becp^*   ---
                        if(okvan) then
                           pref = (0.d0,0.d0)
                           DO jkb=1,nkb
                              nhjkb = nkbtonh(jkb)
                              na = nkbtona(jkb)
                              np = ityp(na)
                              nhjkbm = nh(np)
                              jkb1 = jkb - nhjkb
                              DO j = 1,nhjkbm
                                 pref = pref+CONJG(becp0(jkb,nb))*becp_bp(jkb1+j,mb) &
                                      *q_dk(nhjkb,j,np)*struc(na)
                              ENDDO
                           ENDDO
                      
                           mat(nb,mb) = mat(nb,mb) + pref
                        endif
                     endif !on l_cal
                  ENDDO
               ENDDO

!              --- Calculate matrix determinant ---

! calculate inverse
               CALL zgefa(mat,nbnd,nbnd,ivpt,info)
               CALL errore('forces_us_efield','error in zgefa',abs(info))
               CALL zgedi(mat,nbnd,nbnd,ivpt,cdet,cdwork,1)


!calculate terms
               forces_tmp(:,:)=(0.d0,0.d0)
               if(okvan) then
                  do jkb=1,nkb
                     nhjkb = nkbtonh(jkb)
                     na = nkbtona(jkb)
                     np = ityp(na)
                     nhjkbm = nh(np)
                     jkb1 = jkb - nhjkb
                     do j = 1,nhjkbm
                        do nb=1,nbnd
                           do mb=1,nbnd
                              forces_tmp(:,na)= forces_tmp(:,na)+CONJG(becp0(jkb,nb))*becp_bp(jkb1+j,mb) &
                             *q_dk(nhjkb,j,np)*struc_r(:,na)*mat(mb,nb)
                              forces_tmp(:,na)= forces_tmp(:,na)+CONJG(dbecp0(jkb,nb,:))*becp_bp(jkb1+j,mb) &
                             *q_dk(nhjkb,j,np)*struc(na)*mat(mb,nb)
                               forces_tmp(:,na)= forces_tmp(:,na)+CONJG(becp0(jkb,nb))*dbecp_bp(jkb1+j,mb,:) &
                             *q_dk(nhjkb,j,np)*struc(na)*mat(mb,nb)
                           enddo
                        enddo
                        
                                            
                     enddo
                  
                  enddo
               endif

               forces_bp(:,:)=forces_bp(:,:)+fact*aimag(forces_tmp(:,:))*wstring(istring)

!           --- End of dot products between wavefunctions and betas ---
            ENDIF

!        --- End loop over parallel k-points ---
         END DO
         kpoint=kpoint-1
!     --- End loop over orthogonal k-points ---
      END DO

!  --- End loop over spin ---
   END DO

!  -------------------------------------------------------------------------   !

!  --- Free memory ---
   DEALLOCATE(l_cal)
   DEALLOCATE(pdl_elec)
   DEALLOCATE(mod_elec)
   DEALLOCATE(wstring)
   DEALLOCATE(loc_k)
   DEALLOCATE(phik)
   DEALLOCATE(cphik)
   DEALLOCATE(ln)
   DEALLOCATE(map_g)
   DEALLOCATE(aux)
   DEALLOCATE(aux0)
   DEALLOCATE(psi)
   DEALLOCATE(psi1)
   DEALLOCATE(mat)
!------------------------------------------------------------------------------!

 END SUBROUTINE forces_us_efield

 SUBROUTINE stress_bp_efield (sigmael )
!calculate the stress contribution due to the  electric field
!electronic part

   USE kinds,                ONLY : DP
   USE bp,                   ONLY : efield_cart, el_pol, fc_pol,l3dstring
   USE cell_base, ONLY: at, alat, tpiba, omega, tpiba2
   USE constants, ONLY : pi
   implicit none
   
   REAL(DP), INTENT(out) :: sigmael(3,3)!stress contribution to be calculated

   REAL(DP) :: phases(3)
   INTEGER :: i,j,ipol

   sigmael(:,:)=0.d0
   if(.not.l3dstring ) return
   phases(:)=el_pol(:)/fc_pol(:)
   do ipol=1,3
      do i=1,3
         do j=1,3
            sigmael(i,j)=sigmael(i,j)-efield_cart(i)*at(j,ipol)*phases(ipol)
         enddo
      enddo
   enddo
   sigmael(:,:)=sigmael(:,:)*alat*dsqrt(2.d0)/(2.d0*pi)/omega
   return
 END SUBROUTINE stress_bp_efield


 SUBROUTINE stress_ion_efield (sigmaion )
!calculate the stress contribution due to the  electric field 
!ionic part   
   USE kinds,                ONLY : DP
   USE bp,                   ONLY : efield_cart, ion_pol,l3dstring
   USE cell_base, ONLY: at, alat, omega, bg
   USE constants, ONLY : pi
   implicit none

   REAL(DP), INTENT(out) :: sigmaion(3,3)!stress contribution to be calculated   
   REAL(DP) :: pol_cry(3)
   INTEGER :: i,j,ipol

   sigmaion(:,:)=0.d0
   if(.not.l3dstring ) return
   pol_cry(:)=ion_pol(:)
   call cryst_to_cart (1, pol_cry, at, -1)
   do ipol=1,3
      do i=1,3
         do j=1,3
            sigmaion(i,j)=sigmaion(i,j)-efield_cart(i)*at(j,ipol)*pol_cry(ipol)
         enddo
      enddo
   enddo
   sigmaion(:,:)=sigmaion(:,:)/omega

   return

 END SUBROUTINE stress_ion_efield
