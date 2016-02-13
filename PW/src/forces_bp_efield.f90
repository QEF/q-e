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


  do i=1,nat
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
   USE cell_base,            ONLY : at, alat, tpiba, omega
   USE constants,            ONLY : pi, tpi
   USE gvect,                ONLY : ngm,  g, gcutm, ngm_g, ngmx, ig_l2g
   USE fft_base,             ONLY : dfftp
   USE uspp,                 ONLY : nkb, vkb, okvan
   USE uspp_param,           ONLY : upf, lmaxq, nbetam, nh, nhm
   USE lsda_mod,             ONLY : nspin
   USE klist,                ONLY : nelec, degauss, nks, xk, wk, ngk, igk_k
   USE wvfct,                ONLY : npwx, nbnd
   USE wavefunctions_module, ONLY : evc
   USE bp,                   ONLY : nppstr_3d, mapgm_global, nx_el,mapg_owner
   USE fixed_occ
   USE mp,                   ONLY : mp_sum,mp_barrier
   USE mp_world,             ONLY : world_comm,mpime,nproc
   USE mp_bands,             ONLY : intra_bgrp_comm
   USE becmod,    ONLY : bec_type, becp, calbec,allocate_bec_type, deallocate_bec_type
   USE noncollin_module,     ONLY : noncolin, npol
   USE spin_orb, ONLY: lspinorb
   USE mytime
   USE parallel_include

!  --- Avoid implicit definitions ---
   IMPLICIT NONE

   REAL(DP), INTENT(inout) :: forces_bp(3,nat)!atomic forces to be update 
   INTEGER, INTENT(in) :: pdir!direction of electric field
   REAL(DP), INTENT(in) :: e_field!initensity of the field

!  --- Internal definitions ---
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
   REAL(dp) :: dk(3)
   REAL(dp) :: dkmod
   REAL(dp) :: el_loc
   REAL(dp) :: eps
   REAL(dp) :: fac
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
   COMPLEX(dp), ALLOCATABLE :: aux(:),aux_2(:)
   COMPLEX(dp), ALLOCATABLE :: aux0(:),aux0_2(:)
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

   COMPLEX(kind=DP) :: sca
   COMPLEX(kind=DP), ALLOCATABLE :: aux_g(:),aux_g_mpi(:,:),aux_proc(:,:),aux_rcv(:,:)
   COMPLEX(DP), ALLOCATABLE :: dbecp0(:,:,:), dbecp_bp(:,:,:),vkb1(:,:)
   INTEGER :: ipol
   COMPLEX(DP) :: forces_tmp(3,nat)
   REAL(DP) :: fact
   TYPE(bec_type) :: becp0, becp_bp
   INTEGER :: nspin_eff
   COMPLEX(DP), ALLOCATABLE :: q_dk_so(:,:,:,:)

   COMPLEX(kind=DP), ALLOCATABLE :: fbmata_1(:,:,:),fbmatb_1(:,:,:,:)
   COMPLEX(kind=DP), ALLOCATABLE :: fbmata_2(:,:,:,:),fbmatb_2(:,:,:,:,:)
   COMPLEX(kind=DP), ALLOCATABLE :: fbmata_3(:,:,:),fbmatb_3(:,:,:,:,:)
   COMPLEX(kind=DP), ALLOCATABLE :: dbecp0_ord(:,:,:,:), dbecp_bp_ord(:,:,:,:)
   
   INTEGER :: igg,max_aux,max_aux_proc,iproc
   INTEGER, ALLOCATABLE :: aux_g_mpi_ind(:,:),ind_g(:),aux_proc_ind(:,:),aux_rcv_ind(:,:)
   INTEGER :: req, ierr

!  -------------------------------------------------------------------------   !
!                               INITIALIZATIONS
!  -------------------------------------------------------------------------   !



   allocate(ind_g(nproc))

   nspin_eff=nspin
   if(noncolin) then
      nspin_eff=1
   endif

   ALLOCATE (psi1(npwx*npol,nbnd))
   ALLOCATE (psi(npwx*npol,nbnd))
   ALLOCATE (aux(ngm),aux_2(ngm))
   ALLOCATE (aux0(ngm),aux0_2(ngm))
   ALLOCATE (map_g(npwx))
   ALLOCATE (mat(nbnd,nbnd))
   ALLOCATE (dbecp0( nkb, nbnd*npol, 3 ) ,dbecp_bp( nkb, nbnd*npol, 3 ))
   ALLOCATE (dbecp0_ord( nkb,npol,3, nbnd) ,dbecp_bp_ord( nkb, npol,3,nbnd))

   ALLOCATE( vkb1( npwx, nkb ) )
   ALLOCATE( l_cal(nbnd) )
   if(okvan) then
      CALL allocate_bec_type (nkb,nbnd,becp0)
      CALL allocate_bec_type (nkb,nbnd,becp_bp)
      IF (lspinorb) ALLOCATE(q_dk_so(nhm,nhm,4,ntyp))
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
   ln=0
   DO ng=1,ngm
      mk1=nint(g(1,ng)*at(1,1)+g(2,ng)*at(2,1)+g(3,ng)*at(3,1))
      mk2=nint(g(1,ng)*at(1,2)+g(2,ng)*at(2,2)+g(3,ng)*at(3,2))
      mk3=nint(g(1,ng)*at(1,3)+g(2,ng)*at(2,3)+g(3,ng)*at(3,3))
      ln(mk1,mk2,mk3) = ng
   END DO
   call mp_sum(ln,intra_bgrp_comm)
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
   nkort=nstring/(nspin_eff)

!  --- Allocate memory for arrays ---
   ALLOCATE(phik(nstring))
   ALLOCATE(loc_k(nstring))
   ALLOCATE(cphik(nstring))
   ALLOCATE(wstring(nstring))
   ALLOCATE(pdl_elec(nstring))
   ALLOCATE(mod_elec(nstring))

   FLUSH(stdout)


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
   DO is=1,nspin_eff
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
      IF (lspinorb) CALL transform_qq_so(q_dk,q_dk_so)
   endif

!calculate factor

   
   call factor_a(pdir,at,dkfact)
   fact=dsqrt(2.d0)*e_field*dkfact
   if(nspin_eff==1.and. .not.noncolin) fact=fact*2.d0
   

   
!  -------------------------------------------------------------------------   !
!                   electronic polarization: strings phases                    !
!  -------------------------------------------------------------------------   !

   el_loc=0.d0
   kpoint=0
   zeta=(1.d0,0.d0)
!  --- Start loop over spin ---
   DO is=1,nspin_eff 
      ! l_cal(n) = .true./.false. if n-th state is occupied/empty
      DO nb = 1, nbnd
         IF ( nspin_eff == 2 .AND. tfixed_occ) THEN
            l_cal(nb) = ( f_inp(nb,is) /= 0.0_dp )
         ELSE
            IF(noncolin) THEN
               l_cal(nb) = ( nb <= NINT ( nelec ) )
            ELSE
               l_cal(nb) = ( nb <= NINT ( nelec/2.0_dp ) )
            END IF
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
               ik = nx_el(kpoint-1,pdir)
               npw0   = ngk(ik)
               igk0(:)= igk_k(:,ik)
 
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
                          CALL ZGEMM( 'C', 'N', nkb, nbnd*npol, npw0, ( 1.D0, 0.D0 ),   &
                          vkb1, npwx, psi, npwx, ( 0.D0, 0.D0 ),      &
                          dbecp0(1,1,ipol), nkb )
                          call mp_sum(dbecp0(1:nkb,1:nbnd*npol,ipol),intra_bgrp_comm)
                  ENDDO
               endif

           
!              --- Dot wavefunctions and betas for CURRENT k-point ---
               IF (kpar /= (nppstr_3d(pdir)+1)) THEN

                  ik = nx_el(kpoint,pdir)
                  npw1   = ngk(ik)
                  igk1(:)= igk_k(:,ik)

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
                             CALL ZGEMM( 'C', 'N', nkb, nbnd*npol, npw1, ( 1.D0, 0.D0 ),   &
                             vkb1, npwx, psi1, npwx, ( 0.D0, 0.D0 ),      &
                             dbecp_bp(1,1,ipol), nkb )
                             call mp_sum(dbecp_bp(1:nkb,1:nbnd*npol,ipol),intra_bgrp_comm)
                     ENDDO
                  endif

               ELSE

                  kstart = kpoint-(nppstr_3d(pdir)+1)+1
                  ik = nx_el(kstart,pdir)
                  npw1   = ngk(ik)
                  igk1(:)= igk_k(:,ik)

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
                             CALL ZGEMM( 'C', 'N', nkb, nbnd*npol, npw1, ( 1.D0, 0.D0 ),   &
                             vkb1, npwx, psi1, npwx, ( 0.D0, 0.D0 ),      &
                             dbecp_bp(1,1,ipol), nkb )
                             call mp_sum(dbecp_bp(1:nkb,1:nbnd*npol,ipol),intra_bgrp_comm)
                     ENDDO
                  endif

               ENDIF
           
!              --- Matrix elements calculation ---

               mat=(0.d0,0.d0)
               DO nb=1,nbnd
                  aux=(0.d0,0.d0)
                  aux0=(0.d0,0.d0)
                  IF(noncolin) THEN
                     aux_2=(0.d0,0.d0)
                     aux0_2=(0.d0,0.d0)
                  ENDIF
                  DO ig=1,npw0
                     aux0(igk0(ig))=psi(ig,nb)
                  END DO
                  if(noncolin) then
                     DO ig=1,npw0
                        aux0_2(igk0(ig))=psi(ig+npwx,nb)
                     END DO
                  endif
                  IF (kpar /= (nppstr_3d(pdir)+1)) THEN
                     DO mb=1,nbnd
                        IF ( .NOT. l_cal(nb) .OR. .NOT. l_cal(mb) ) THEN
                           IF ( nb == mb )  mat(nb,mb)=1.d0
                        ELSE
                           
                           do ig=1,npw1
                              aux(igk1(ig))=psi1(ig,mb)
                           enddo
                           IF(noncolin) THEN
                              do ig=1,npw1
                                 aux_2(igk1(ig))=psi1(ig+npwx,mb)
                              enddo
                           END IF
                           
                           mat(nb,mb) = zdotc(ngm,aux0,1,aux,1)
                           
                           if(noncolin) then
                              mat(nb,mb) = mat(nb,mb) + zdotc(ngm,aux0_2,1,aux_2,1)
                           endif
                           call mp_sum( mat(nb,mb), intra_bgrp_comm )
                        END IF
                     END DO
                  END IF
               END DO
           
               IF (kpar == (nppstr_3d(pdir)+1) ) THEN
                     
! allocate global array
                  
                  allocate(aux_g(ngm_g),aux_g_mpi(ngmx,nproc),aux_g_mpi_ind(ngmx,nproc))
                  do ipol=0,npol-1
                     do mb=1,nbnd
                        
                        aux_g_mpi=0.d0
                        aux_g_mpi_ind=0
                        ind_g=0

                        do ig=1,npw1
                           igg=mapgm_global(ig_l2g(igk1(ig)),pdir)
                           ind_g(mapg_owner(1,igg))=ind_g(mapg_owner(1,igg))+1
                           aux_g_mpi(ind_g(mapg_owner(1,igg)),mapg_owner(1,igg))=psi1(ig+npwx*ipol,mb)
                           aux_g_mpi_ind(ind_g(mapg_owner(1,igg)),mapg_owner(1,igg))=mapg_owner(2,igg)
                        enddo
                        max_aux=0

                        do iproc=1,nproc
                           if(iproc/=mpime+1) then
                              max_aux_proc=0
                              do ig=1,ngmx
                                 if(aux_g_mpi_ind(ig,iproc) > 0) then
                                    max_aux_proc=max_aux_proc+1
                                 else
                                    exit
                                 endif
                              enddo
                              if(max_aux_proc>max_aux) max_aux=max_aux_proc
                           endif
                        enddo
                        max_aux_proc=max_aux

#if defined (__MPI)                        
                        CALL MPI_ALLREDUCE( max_aux_proc,max_aux,1,MPI_INTEGER, MPI_MAX,intra_bgrp_comm, req,IERR )
#endif
                        allocate(aux_proc(max_aux,nproc),aux_proc_ind(max_aux,nproc))
                        allocate(aux_rcv(max_aux,nproc),aux_rcv_ind(max_aux,nproc))
                        aux_proc=(0.d0,0.d0)
                        aux_proc_ind=0
                        do iproc=1,nproc
                           if(iproc/=mpime+1) then
                              do ig=1,max_aux
                                 if(aux_g_mpi_ind(ig,iproc) > 0) then
                                    aux_proc(ig,iproc)=aux_g_mpi(ig,iproc)
                                    aux_proc_ind(ig,iproc)=aux_g_mpi_ind(ig,iproc)
                                 else
                                    exit
                                 end if
                              enddo
                           endif
                        enddo

                        
#if defined (__MPI)
                        CALL MPI_ALLTOALL( aux_proc, max_aux, MPI_DOUBLE_COMPLEX,  &
                             aux_rcv, max_aux, MPI_DOUBLE_COMPLEX, intra_bgrp_comm, ierr )
                        CALL MPI_ALLTOALL( aux_proc_ind, max_aux, MPI_INTEGER,  &
                             aux_rcv_ind, max_aux, MPI_INTEGER, intra_bgrp_comm, ierr )
#else
                        aux_rcv(1:max_aux,1)=aux_proc(1:max_aux,1)
                        aux_rcv_ind(1:max_aux,1)=aux_proc_ind(1:max_aux,1)
#endif

                        
                        do nb=1,nbnd
                           IF ( .NOT. l_cal(nb) .OR. .NOT. l_cal(mb) ) THEN
                              IF ( nb == mb )  mat(nb,mb)=1.d0
                           ELSE

                              aux=(0.d0,0.d0)
                              aux0=(0.d0,0.d0)
                              IF(noncolin) THEN
                                 aux_2=(0.d0,0.d0)
                                 aux0_2=(0.d0,0.d0)
                              ENDIF
                              DO ig=1,npw0
                                 aux0(igk0(ig))=psi(ig,nb)
                              END DO
                              if(noncolin) then
                              DO ig=1,npw0
                                 aux0_2(igk0(ig))=psi(ig+npwx,nb)
                              END DO
                           endif
                           sca=0.d0
                           do iproc=1,nproc
                              if(iproc/=mpime+1) then
                                 do ig=1,max_aux
                                    if(aux_rcv_ind(ig,iproc)/=0) then
                                       if(aux_rcv_ind(ig,iproc)<0.or.aux_rcv_ind(ig,iproc)> ngm) then
                                          write(stdout,*) 'OH BOY', aux_rcv_ind(ig,iproc)
                                       else
                                          if(ipol==0) then
                                             sca=sca+conjg(aux0(aux_rcv_ind(ig,iproc)))*aux_rcv(ig,iproc)
                                          else
                                             sca=sca+conjg(aux0_2(aux_rcv_ind(ig,iproc)))*aux_rcv(ig,iproc)
                                          endif
                                       endif
                                    else
                                       exit
                                    endif
                                 enddo
                              endif
                           enddo
                              
                           do ig=1,ngmx
                              if(aux_g_mpi_ind(ig,mpime+1)/=0) then
                                 if(aux_g_mpi_ind(ig,mpime+1)<0.or.aux_g_mpi_ind(ig,mpime+1)>ngm) then 
                                    write(stdout,*) 'OH BOY',aux_g_mpi_ind(ig,mpime+1)
                                 else
                                    if(ipol==0) then
                                       sca=sca+conjg(aux0(aux_g_mpi_ind(ig,mpime+1)))*aux_g_mpi(ig,mpime+1)
                                    else
                                       sca=sca+conjg(aux0_2(aux_g_mpi_ind(ig,mpime+1)))*aux_g_mpi(ig,mpime+1)
                                    endif
                                 endif
                              else
                                 exit
                              endif
                           enddo
                                 
                           call mp_sum(sca,intra_bgrp_comm)
                           mat(nb,mb)=mat(nb,mb)+sca
                        endif
                        enddo
                        deallocate(aux_proc,aux_proc_ind)
                        deallocate(aux_rcv,aux_rcv_ind)
                           
                     enddo
                  enddo
                  deallocate(aux_g,aux_g_mpi,aux_g_mpi_ind)
                  
                         
               ENDIF
               DO nb=1,nbnd
                  do mb=1,nbnd 
                     IF ( l_cal(nb) .AND. l_cal(mb) ) THEN
                      

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
                                 if(lspinorb) then
                                    pref = pref+CONJG(becp0%nc(jkb,1,nb))*becp_bp%nc(jkb1+j,1,mb) &
                                         *q_dk_so(nhjkb,j,1,np)*struc(na)
                                    pref = pref+CONJG(becp0%nc(jkb,1,nb))*becp_bp%nc(jkb1+j,2,mb) &
                                         *q_dk_so(nhjkb,j,2,np)*struc(na)
                                    pref = pref+CONJG(becp0%nc(jkb,2,nb))*becp_bp%nc(jkb1+j,1,mb) &
                                         *q_dk_so(nhjkb,j,3,np)*struc(na)
                                    pref = pref+CONJG(becp0%nc(jkb,2,nb))*becp_bp%nc(jkb1+j,2,mb) &
                                         *q_dk_so(nhjkb,j,4,np)*struc(na)

                                 else
                                    pref = pref+CONJG(becp0%k(jkb,nb))*becp_bp%k(jkb1+j,mb) &
                                         *q_dk(nhjkb,j,np)*struc(na)
                                 endif
                              ENDDO
                           ENDDO
                      
                           mat(nb,mb) = mat(nb,mb) + pref
                        endif
                     endif !on l_cal
                  ENDDO
               ENDDO
               
!              --- Calculate matrix determinant ---

! calculate inverse
!              

               CALL zgefa(mat,nbnd,nbnd,ivpt,info)
               CALL errore('forces_us_efield','error in zgefa',abs(info))
               CALL zgedi(mat,nbnd,nbnd,ivpt,cdet,cdwork,1)

               DO nb = 1, nbnd
                  DO mb = 1, nbnd
                     if(.not.l_cal(nb).or. .not. l_cal(mb)) mat(mb,nb)=(0.d0,0.d0)
                  END DO
               END DO

               
!calculate terms

               forces_tmp(:,:)=(0.d0,0.d0)

               if(okvan) then
                  allocate(fbmatb_1(nkb,npol,nkb,npol),fbmata_1(nbnd,nkb,npol))
                  allocate(fbmatb_2(nkb,npol,nkb,npol,3),fbmata_2(nbnd,nkb,npol,3))
                  allocate(fbmatb_3(nkb,npol,3,nkb,npol),fbmata_3(nbnd,nkb,npol))

                  if(lspinorb) then
                     call ZGEMM('N','C',nbnd,nkb*npol,nbnd,(1.d0,0.d0),&
                          &mat,nbnd,becp0%nc(1,1,1),nkb*npol,(0.d0,0.d0),fbmata_1,nbnd)
                     call ZGEMM('N','N',nkb*npol,nkb*npol,nbnd,(1.d0,0.d0),&
                          &becp_bp%nc(1,1,1),nkb*npol,fbmata_1,nbnd,(0.d0,0.d0),fbmatb_1,nkb*npol)
                     do ipol=1,npol
                        do nb=1,nbnd
                           dbecp0_ord(1:nkb,ipol,1:3,nb)=dbecp0(1:nkb,(nb-1)*npol+ipol,1:3)
                           dbecp_bp_ord(1:nkb,ipol,1:3,nb)=dbecp_bp(1:nkb,(nb-1)*npol+ipol,1:3)
                        enddo
                     enddo

                     call ZGEMM('N','C',nbnd,nkb*npol*3,nbnd,(1.d0,0.d0),&
                          &mat,nbnd,dbecp0_ord,nkb*npol*3,(0.d0,0.d0),fbmata_2,nbnd)
                     call ZGEMM('N','N',nkb*npol,nkb*npol*3,nbnd,(1.d0,0.d0),&
                          &becp_bp%nc(1,1,1),nkb*npol,fbmata_2,nbnd,(0.d0,0.d0),fbmatb_2,nkb*npol)

                     call ZGEMM('N','C',nbnd,nkb*npol,nbnd,(1.d0,0.d0),&
                          &mat,nbnd,becp0%nc(1,1,1),nkb*npol,(0.d0,0.d0),fbmata_3,nbnd)
                     call ZGEMM('N','N',nkb*npol*3,nkb*npol,nbnd,(1.d0,0.d0),&
                          &dbecp_bp_ord,nkb*npol*3,fbmata_3,nbnd,(0.d0,0.d0),fbmatb_3,nkb*npol*3)


                     do jkb=1,nkb
                        nhjkb = nkbtonh(jkb)
                        na = nkbtona(jkb)
                        np = ityp(na)
                        nhjkbm = nh(np)
                        jkb1 = jkb - nhjkb
                        do j = 1,nhjkbm
                              forces_tmp(1:3,na)= forces_tmp(1:3,na)+ &
                                     & q_dk_so(nhjkb,j,1,np)*struc_r(1:3,na)*fbmatb_1(jkb1+j,1,jkb,1)
                              forces_tmp(1:3,na)= forces_tmp(1:3,na)+ &
                                   q_dk_so(nhjkb,j,1,np)*struc(na)*fbmatb_2(jkb1+j,1,jkb,1,1:3)
                              forces_tmp(1:3,na)= forces_tmp(1:3,na)+ &
                                   q_dk_so(nhjkb,j,1,np)*struc(na)*fbmatb_3(jkb1+j,1,1:3,jkb,1)

                              forces_tmp(1:3,na)= forces_tmp(1:3,na)+ &
                                   & q_dk_so(nhjkb,j,2,np)*struc_r(1:3,na)*fbmatb_1(jkb1+j,2,jkb,1)
                              forces_tmp(1:3,na)= forces_tmp(1:3,na)+ &
                                   q_dk_so(nhjkb,j,2,np)*struc(na)*fbmatb_2(jkb1+j,2,jkb,1,1:3)
                              forces_tmp(1:3,na)= forces_tmp(1:3,na)+ &
                                   q_dk_so(nhjkb,j,2,np)*struc(na)*fbmatb_3(jkb1+j,2,1:3,jkb,1)

                              forces_tmp(1:3,na)= forces_tmp(1:3,na)+ &
                                   & q_dk_so(nhjkb,j,3,np)*struc_r(1:3,na)*fbmatb_1(jkb1+j,1,jkb,2)
                              forces_tmp(1:3,na)= forces_tmp(1:3,na)+ &
                                   q_dk_so(nhjkb,j,3,np)*struc(na)*fbmatb_2(jkb1+j,1,jkb,2,1:3)
                              forces_tmp(1:3,na)= forces_tmp(1:3,na)+ &
                                   q_dk_so(nhjkb,j,3,np)*struc(na)*fbmatb_3(jkb1+j,1,1:3,jkb,2)
                              
                              forces_tmp(1:3,na)= forces_tmp(1:3,na)+ &
                                   & q_dk_so(nhjkb,j,4,np)*struc_r(1:3,na)*fbmatb_1(jkb1+j,2,jkb,2)
                              forces_tmp(1:3,na)= forces_tmp(1:3,na)+ &
                                   q_dk_so(nhjkb,j,4,np)*struc(na)*fbmatb_2(jkb1+j,2,jkb,2,1:3)
                              forces_tmp(1:3,na)= forces_tmp(1:3,na)+ &
                                   q_dk_so(nhjkb,j,4,np)*struc(na)*fbmatb_3(jkb1+j,2,1:3,jkb,2)


                        enddo
                     enddo
                  endif
                  if(.not.lspinorb) then
                  do jkb=1,nkb
                     nhjkb = nkbtonh(jkb)
                     na = nkbtona(jkb)
                     np = ityp(na)
                     nhjkbm = nh(np)
                     jkb1 = jkb - nhjkb
                     do j = 1,nhjkbm
                        do nb=1,nbnd
                           do mb=1,nbnd
                              if(lspinorb) then
                                ! forces_tmp(1:3,na)= forces_tmp(1:3,na)+CONJG(becp0%nc(jkb,1,nb))*becp_bp%nc(jkb1+j,1,mb) &
                                !      *q_dk_so(nhjkb,j,1,np)*struc_r(1:3,na)*mat(mb,nb)
                                ! forces_tmp(1:3,na)= forces_tmp(1:3,na)+CONJG(dbecp0(jkb,(nb-1)*npol+1,1:3)) &
                                !      *becp_bp%nc(jkb1+j,1,mb)*q_dk_so(nhjkb,j,1,np)*struc(na)*mat(mb,nb)
                                ! forces_tmp(1:3,na)= forces_tmp(1:3,na)+CONJG(becp0%nc(jkb,1,nb)) &
                                ! *dbecp_bp(jkb1+j,(mb-1)*npol+1,1:3)*q_dk_so(nhjkb,j,1,np)*struc(na)*mat(mb,nb)

                                ! forces_tmp(1:3,na)= forces_tmp(1:3,na)+CONJG(becp0%nc(jkb,1,nb))*becp_bp%nc(jkb1+j,2,mb) &
                                !      *q_dk_so(nhjkb,j,2,np)*struc_r(1:3,na)*mat(mb,nb)
                                ! forces_tmp(1:3,na)= forces_tmp(1:3,na)+CONJG(dbecp0(jkb,(nb-1)*npol+1,1:3))&
                                !      *becp_bp%nc(jkb1+j,2,mb)*q_dk_so(nhjkb,j,2,np)*struc(na)*mat(mb,nb)
                                ! forces_tmp(1:3,na)= forces_tmp(1:3,na)+CONJG(becp0%nc(jkb,1,nb))&
                                !      *dbecp_bp(jkb1+j,(mb-1)*npol+2,1:3)*q_dk_so(nhjkb,j,2,np)*struc(na)*mat(mb,nb)

                                ! forces_tmp(1:3,na)= forces_tmp(1:3,na)+CONJG(becp0%nc(jkb,2,nb))*becp_bp%nc(jkb1+j,1,mb) &
                                !      *q_dk_so(nhjkb,j,3,np)*struc_r(1:3,na)*mat(mb,nb)
                                ! forces_tmp(1:3,na)= forces_tmp(1:3,na)+CONJG(dbecp0(jkb,(nb-1)*npol+2,1:3))&
                                !      *becp_bp%nc(jkb1+j,1,mb)*q_dk_so(nhjkb,j,3,np)*struc(na)*mat(mb,nb)
                                ! forces_tmp(1:3,na)= forces_tmp(1:3,na)+CONJG(becp0%nc(jkb,2,nb)) &
                                !      *dbecp_bp(jkb1+j,(mb-1)*npol+1,1:3)*q_dk_so(nhjkb,j,3,np)*struc(na)*mat(mb,nb)

                                 
                                 !forces_tmp(1:3,na)= forces_tmp(1:3,na)+CONJG(becp0%nc(jkb,2,nb))*becp_bp%nc(jkb1+j,2,mb) &
                                 !     *q_dk_so(nhjkb,j,4,np)*struc_r(1:3,na)*mat(mb,nb)
                                 !forces_tmp(1:3,na)= forces_tmp(1:3,na)+CONJG(dbecp0(jkb,(nb-1)*npol+2,1:3))&
                                 !     *becp_bp%nc(jkb1+j,2,mb)*q_dk_so(nhjkb,j,4,np)*struc(na)*mat(mb,nb)
                                 !forces_tmp(1:3,na)= forces_tmp(1:3,na)+CONJG(becp0%nc(jkb,2,nb))&
                                 !     *dbecp_bp(jkb1+j,(mb-1)*npol+2,1:3)*q_dk_so(nhjkb,j,4,np)*struc(na)*mat(mb,nb)

                              else
                                 forces_tmp(:,na)= forces_tmp(:,na)+CONJG(becp0%k(jkb,nb))*becp_bp%k(jkb1+j,mb) &
                                      *q_dk(nhjkb,j,np)*struc_r(:,na)*mat(mb,nb)
                                 forces_tmp(:,na)= forces_tmp(:,na)+CONJG(dbecp0(jkb,nb,:))*becp_bp%k(jkb1+j,mb) &
                                      *q_dk(nhjkb,j,np)*struc(na)*mat(mb,nb)
                                 forces_tmp(:,na)= forces_tmp(:,na)+CONJG(becp0%k(jkb,nb))*dbecp_bp(jkb1+j,mb,:) &
                                      *q_dk(nhjkb,j,np)*struc(na)*mat(mb,nb)
                              endif
                           enddo
                        enddo
                     
                                            
                     enddo
                  end do
               end if
                  deallocate(fbmata_1,fbmatb_1)
                  deallocate(fbmata_2,fbmatb_2)
                  deallocate(fbmata_3,fbmatb_3)
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
   DEALLOCATE(aux,aux_2)
   DEALLOCATE(aux0,aux0_2)
   DEALLOCATE(psi)
   DEALLOCATE(psi1)
   DEALLOCATE(mat)
   if(okvan) then
      call deallocate_bec_type(becp0)
      call deallocate_bec_type(becp_bp)
      if(lspinorb) deallocate(q_dk_so)
   endif
   DEALLOCATE(dbecp0,dbecp0_ord,dbecp_bp,dbecp_bp_ord)

   DEALLOCATE(ind_g)


!------------------------------------------------------------------------------!

 END SUBROUTINE forces_us_efield

 SUBROUTINE stress_bp_efield (sigmael )
!calculate the stress contribution due to the  electric field
!electronic part

   USE kinds,                ONLY : DP
   USE bp,                   ONLY : efield_cart, el_pol, fc_pol,l3dstring
   USE cell_base, ONLY: at, alat, tpiba, omega
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
