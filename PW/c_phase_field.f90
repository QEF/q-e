!
! Copyright (C) 2001-2004 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

! this routine is used to calculate the electronic polarization
! when a finite  electric field describe through the modern
! theory of the polarization is applied.
! is very close to the routine c_phase in bp_c_phase
! however the numbering of the k-points in the strings is different

#include "f_defs.h"
!======================================================================!

SUBROUTINE c_phase_field

!----------------------------------------------------------------------!

!   Geometric phase calculation along a strip of nppstr k-points
!   averaged over a 2D grid of nkort k-points ortogonal to nppstr 


!  --- Make use of the module with common information ---
   USE kinds,                ONLY : DP
   USE parameters,           ONLY : nbrx
   USE io_global,            ONLY : stdout
   USE io_files,             ONLY : iunwfc, nwordwfc
   USE ions_base,            ONLY : nat, ntyp => nsp, ityp, tau, zv, atm
   USE cell_base,            ONLY : at, alat, tpiba, omega, tpiba2
   USE constants,            ONLY : pi, tpi
   USE gvect,                ONLY : ngm, nr1, nr2, nr3, nrx1, nrx2, nrx3, &
                                    ecutwfc, g, gcutm
   USE uspp,                 ONLY : nkb, vkb, okvan
   USE uspp_param,           ONLY : lmaxq, nh, nhm, tvanp
   USE lsda_mod,             ONLY : nspin
   USE klist,                ONLY : nelec, degauss, nks, xk, wk
   USE wvfct,                ONLY : npwx, npw, nbnd
   USE wavefunctions_module, ONLY : evc
   USE bp
   USE fixed_occ


!  --- Make use of the module with common information ---


!  --- Avoid implicit definitions ---
   IMPLICIT NONE

!  --- Internal definitions ---
   INTEGER :: i
   INTEGER :: igk1(npwx)
   INTEGER :: igk0(npwx)
   INTEGER :: ig
   INTEGER :: ind1
   INTEGER :: info
   INTEGER :: is
   INTEGER :: istring
   INTEGER :: iv
   INTEGER :: ivpt(nbnd)
   INTEGER :: j
   INTEGER :: jkb
   INTEGER :: jkb_bp
   INTEGER :: jkb1
   INTEGER :: job
   INTEGER :: jv
   INTEGER :: kindex
   INTEGER :: kort
   INTEGER :: kpar
   INTEGER :: kpoint
   INTEGER :: kstart
   INTEGER :: mb
   INTEGER :: mk1
   INTEGER :: mk2
   INTEGER :: mk3
   INTEGER , ALLOCATABLE :: mod_elec(:)
   INTEGER :: mod_elec_dw
   INTEGER :: mod_elec_tot
   INTEGER :: mod_elec_up
   INTEGER :: mod_ion(nat)
   INTEGER :: mod_ion_dw
   INTEGER :: mod_ion_tot
   INTEGER :: mod_ion_up
   INTEGER :: mod_tot
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
   LOGICAL :: lodd
   REAL(dp) :: dk(3)
   REAL(dp) :: dkmod
   REAL(dp) :: el_loc
   REAL(dp) :: eps
   REAL(dp) :: fac
   REAL(dp) :: g2kin_bp(npwx)
   REAL(dp) :: gpar(3)
   REAL(dp) :: gtr(3)
   REAL(dp) :: gvec
   REAL(dp) :: ln(-nr1:nr1,-nr2:nr2,-nr3:nr3)
   REAL(dp), ALLOCATABLE :: loc_k(:)
   REAL(dp) , ALLOCATABLE :: pdl_elec(:)
   REAL(dp) :: pdl_elec_dw
   REAL(dp) :: pdl_elec_tot
   REAL(dp) :: pdl_elec_up
   REAL(dp) :: pdl_ion(nat)
   REAL(dp) :: pdl_ion_dw
   REAL(dp) :: pdl_ion_tot
   REAL(dp) :: pdl_ion_up
   REAL(dp) :: pdl_tot
   REAL(dp) , ALLOCATABLE :: phik(:)
   REAL(dp) :: phidw
   REAL(dp) :: phiup
   REAL(dp) :: rmod
   REAL(dp) :: qrad_dk(nbrx,nbrx,lmaxq,ntyp)
   REAL(dp) :: upol(3)
   REAL(dp) :: weight
   REAL(dp), ALLOCATABLE :: wstring(:)
   REAL(dp) :: ylm_dk(lmaxq*lmaxq)
   REAL(dp) :: zeta_mod
   REAL(dp) :: chi
   REAL(dp) :: pola, pola_ion
   COMPLEX(dp) :: aux(ngm)
   COMPLEX(dp) :: aux0(ngm)
   COMPLEX(dp) :: becp0(nkb,nbnd)
   COMPLEX(dp) :: becp_bp(nkb,nbnd)
   COMPLEX(dp) :: cdet(2)
   COMPLEX(dp) :: cdwork(nbnd)
   COMPLEX(dp) :: cave
   COMPLEX(dp) :: cave_dw
   COMPLEX(dp) :: cave_up
   COMPLEX(dp) , ALLOCATABLE :: cphik(:)
   COMPLEX(dp) :: det
   COMPLEX(dp) :: dtheta
   COMPLEX(dp) :: mat(nbnd,nbnd)
   COMPLEX(dp) :: pref
   COMPLEX(dp) :: q_dk(nhm,nhm,ntyp)
   COMPLEX(dp) :: struc(nat)
   COMPLEX(dp) :: theta0
   COMPLEX(dp) :: zdotc
   COMPLEX(dp) :: zeta

   COMPLEX(dp) :: psi(npwx,nbnd)
   COMPLEX(dp) :: psi1(npwx,nbnd)
   COMPLEX(dp) :: zeta_loc


   INTEGER ipivi(nbnd,nbnd),ii

   LOGICAL l_cal!flag for doing mat calculation
   INTEGER, ALLOCATABLE :: map_g(:)


!  -------------------------------------------------------------------------   !
!                               INITIALIZATIONS
!  -------------------------------------------------------------------------   !

  allocate(map_g(npwx))

  pola=0.d0!set to 0 electronic polarization   

!  --- Check that we are working with an insulator with no empty bands ---
!   IF ((degauss > 0.01) .OR. (nbnd /= nelec/2)) CALL errore('c_phase', &
!                'Polarization only for insulators and no empty bands',1)

IF ((degauss > 0.01) .OR. (nbnd /= nelec/2)) &
   &write(stdout,*) 'PAY ATTENTION: EL FIELD AND OCCUPATIONS'

!  --- Define a small number ---
   eps=1.0E-6_dp

!  --- Recalculate FFT correspondence (see ggen.f90) ---
   DO ng=1,ngm
      mk1=nint(g(1,ng)*at(1,1)+g(2,ng)*at(2,1)+g(3,ng)*at(3,1))
      mk2=nint(g(1,ng)*at(1,2)+g(2,ng)*at(2,2)+g(3,ng)*at(3,2))
      mk3=nint(g(1,ng)*at(1,3)+g(2,ng)*at(2,3)+g(3,ng)*at(3,3))
      ln(mk1,mk2,mk3) = ng
   END DO

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
!  --- Get the number of strings ---
   nstring=nks/nppstr
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
   if(nppstr .ne. 1) then
      gpar(1)=(xk(1,nppstr)-xk(1,1))*DBLE(nppstr)/DBLE(nppstr-1)
      gpar(2)=(xk(2,nppstr)-xk(2,1))*DBLE(nppstr)/DBLE(nppstr-1)
      gpar(3)=(xk(3,nppstr)-xk(3,1))*DBLE(nppstr)/DBLE(nppstr-1)
      gvec=dsqrt(gpar(1)**2+gpar(2)**2+gpar(3)**2)*tpiba
   else
      gpar(1)=0.
      gpar(2)=0.
      gpar(3)=0.
      gpar(gdir)=1./at(gdir,gdir)!
      gvec=tpiba/sqrt(at(gdir,1)**2.+at(gdir,2)**2.+at(gdir,3)**2.)
   endif
      
      
!  --- Find vector between consecutive points in strings ---
   if(nppstr.ne.1) then
      dk(1)=xk(1,2)-xk(1,1)
      dk(2)=xk(2,2)-xk(2,1) 
      dk(3)=xk(3,2)-xk(3,1)
      dkmod=SQRT(dk(1)**2+dk(2)**2+dk(3)**2)*tpiba!orthorombic cell
     
   else!caso punto gamma, per adesso solo cella cubica
      dk(1)=0.
      dk(2)=0.
      dk(3)=0.
      dk(gdir)=1./at(gdir,gdir)
      dkmod=tpiba/sqrt(at(gdir,1)**2.+at(gdir,2)**2.+at(gdir,3)**2.)
   endif
   
  
 
!  -------------------------------------------------------------------------   !
!                   electronic polarization: weight strings                    !
!  -------------------------------------------------------------------------   !

!  --- Calculate string weights, normalizing to 1 (no spin) or 1+1 (spin) ---
   DO is=1,nspin
      weight=0.0_dp
      DO kort=1,nkort
         istring=kort+(is-1)*nkort
         wstring(istring)=wk(nppstr*istring)
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
      struc(na)=CMPLX(cos(fac),-sin(fac))
   END DO
  
  

!  -------------------------------------------------------------------------   !
!                     electronic polarization: form factor                     !
!  -------------------------------------------------------------------------   !
   if(okvan) then
!  --- Calculate Bessel transform of Q_ij(|r|) at dk [Q_ij^L(|r|)] ---
      CALL calc_btq(dkmod,qrad_dk,0)

!  --- Calculate the q-space real spherical harmonics at dk [Y_LM] --- 
      CALL ylm_q(lmaxq*lmaxq,dk,dkmod,ylm_dk)!questa no funzia perche'??
    
      
!  --- Form factor: 4 pi sum_LM c_ij^LM Y_LM(Omega) Q_ij^L(|r|) ---
!      CALL setv(2*nhm*nhm*ntyp,0.d0,q_dk,1)
      q_dk=(0.d0,0.d0)
      DO np =1, ntyp
         if(tvanp(np)) then
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
 
   
!  -------------------------------------------------------------------------   !
!                   electronic polarization: strings phases                    !
!  -------------------------------------------------------------------------   !

   el_loc   = 0.d0
   kpoint=0
   zeta=(1.d0,0.d0)
!  --- Start loop over spin ---
   DO is=1,nspin 
       
!     --- Start loop over orthogonal k-points ---
      DO kort=1,nkort
         zeta_loc=(1.d0,0.d0)
!        --- Index for this string ---
         istring=kort+(is-1)*nkort

!        --- Initialize expectation value of the phase operator ---
      
         zeta_mod = 1.d0

!        --- Start loop over parallel k-points ---
         DO kpar = 1,nppstr+1

!           --- Set index of k-point ---
            kpoint = kpoint + 1

!           --- Calculate dot products between wavefunctions and betas ---
            IF (kpar /= 1 ) THEN
             
!              --- Dot wavefunctions and betas for PREVIOUS k-point ---
               CALL gk_sort(xk(1,kpoint-1),ngm,g,ecutwfc/tpiba2, &
                            npw0,igk0,g2kin_bp) 
               CALL davcio(psi,nwordwfc,iunwfc,kpoint-1,-1)
               if(okvan) then
                  CALL init_us_2 (npw0,igk0,xk(1,kpoint-1),vkb)
                  CALL ccalbec(nkb, npwx, npw0, nbnd, becp0, vkb, psi)
               endif
              
!              --- Dot wavefunctions and betas for CURRENT k-point ---
               IF (kpar /= (nppstr+1)) THEN
                  CALL gk_sort(xk(1,kpoint),ngm,g,ecutwfc/tpiba2, &
                               npw1,igk1,g2kin_bp)        
                  CALL davcio(psi1,nwordwfc,iunwfc,kpoint,-1)
                  if(okvan) then
                     CALL init_us_2 (npw1,igk1,xk(1,kpoint),vkb)
                     CALL ccalbec(nkb,npwx,npw1,nbnd,becp_bp,vkb,psi1)
                  endif
               ELSE
                  kstart = kpoint-(nppstr+1)+1
                  CALL gk_sort(xk(1,kstart),ngm,g,ecutwfc/tpiba2, &
                               npw1,igk1,g2kin_bp)  
                  CALL davcio(psi1,nwordwfc,iunwfc,kstart,-1)
                  if(okvan) then
                     CALL init_us_2 (npw1,igk1,xk(1,kstart),vkb)
                     CALL ccalbec(nkb,npwx,npw1,nbnd,becp_bp,vkb,psi1)
                  endif
               ENDIF

           
!              --- Matrix elements calculation ---
!               CALL setv(2*nbnd*nbnd,0.d0,mat,1)
! 






               
               IF (kpar == (nppstr+1)) THEN
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
                        IF ((ABS(g(1,ng)-gtr(1)) > eps) .OR. &
                             (ABS(g(2,ng)-gtr(2)) > eps) .OR. &
                             (ABS(g(3,ng)-gtr(3)) > eps)) THEN
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
!                     CALL setv(2*ngm,0.d0,aux,1)
!
!added support for spin polarized case
                     l_cal=.true.
                     if( nspin==2 .and. tfixed_occ) then
                        if(f_inp(nb,is)==0.d0 .or. f_inp(mb,is)==0.d0) then
                           l_cal=.false.
                           if(nb==mb) then
                              mat(nb,mb)=1.d0
                           else
                              mat(nb,mb)=0.d0
                           endif
                        endif
                     endif
                     if(l_cal) then
!                        CALL setv(2*ngm,0.d0,aux0,1)
                        aux=(0.d0,0.d0)
                        aux0=(0.d0,0.d0)
                        DO ig=1,npw0
                           aux0(igk0(ig))=psi(ig,nb)
                        END DO
                 
                   
                        IF (kpar /= (nppstr+1)) THEN
                           do ig=1,npw1
                              aux(igk1(ig))=psi1(ig,mb)
                           enddo
                        ELSE
                           do ig=1,npw1

                              aux(map_g(ig))=psi1(ig,mb)
                           enddo
                           
                        ENDIF
                     
                        mat(nb,mb) = zdotc(ngm,aux0,1,aux,1)

                        call reduce(2,mat(nb,mb))
             
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
!               CALL zgefa(mat,nbnd,nbnd,ivpt,info)
!               CALL errore('c_phase','error in zgefa',abs(info))
!               job=10
!               CALL zgedi(mat,nbnd,nbnd,ivpt,cdet,cdwork,job)
!               det=cdet(1)*10.**cdet(2)

               det=(1.,0.)
               call zgetrf(nbnd,nbnd,mat,nbnd,ipivi,info)
               CALL errore('c_phase','error in zgetrf',abs(info))
               do ii=1,nbnd
                  if(ii.ne.ipivi(ii,1)) det=-det
               enddo
               do ii=1,nbnd
                  det = det*mat(ii,ii)
               enddo

!               write(6,*) "LogDet:", LOG(det)
!              --- Multiply by the already calculated determinants ---
               zeta=zeta*det
               zeta_loc=zeta_loc*det
               

!           --- End of dot products between wavefunctions and betas ---
            ENDIF

!        --- End loop over parallel k-points ---
         END DO 
         pola=pola+wstring(istring)*aimag(log(zeta_loc))

         kpoint=kpoint-1
!        --- Calculate the phase for this string ---
         phik(istring)=AIMAG(LOG(zeta))
         cphik(istring)=COS(phik(istring))*(1.0_dp,0.0_dp) &
                     +SIN(phik(istring))*(0.0_dp,1.0_dp)

!        --- Calculate the localization for current kort ---
         zeta_mod= DBLE(CONJG(zeta)*zeta)
         loc_k(istring)= - (nppstr-1) / gvec**2 / nbnd *log(zeta_mod)

!     --- End loop over orthogonal k-points ---
      END DO

!  --- End loop over spin ---
   END DO
!-----calculate polarization
!-----the factor 2. is because of spin
   if(nspin==1) pola=pola*2.d0
   pola=pola/(gpar(gdir)*tpiba)

!write output
   write(stdout,*)
   write(stdout,*) "    Expectation value of exp(iGx):",zeta
   write(stdout,*) "    Electronic Dipole per cell (a.u.)",pola
 


!  -------------------------------------------------------------------------   !
!                              ionic polarization                              !
!  -------------------------------------------------------------------------   !

  pola_ion=0.d0
  DO na=1,nat
    pola_ion=pola_ion+zv(ityp(na))*tau(gdir,na)
  END DO

   write(stdout,*) "    Ionic Dipole per cell (a.u.)",pola_ion
   write(stdout,*)

!  -------------------------------------------------------------------------   !

!  --- Free memory ---
   DEALLOCATE(pdl_elec)
   DEALLOCATE(mod_elec)
   DEALLOCATE(wstring)
   DEALLOCATE(loc_k)
   DEALLOCATE(phik)
   DEALLOCATE(cphik)
   DEALLOCATE(map_g)

!------------------------------------------------------------------------------!

END SUBROUTINE c_phase_field

!==============================================================================!
