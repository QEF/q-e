!
! Copyright (C) 2001-2012 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
SUBROUTINE forces_ion_efield( forces_bp, pdir, e_field )
  !--------------------------------------------------------------------------
  !! Calculate ionic contribution, which is in the a_gdir direction.
  !
  USE kinds,                ONLY : DP
  USE cell_base,            ONLY : at
  USE ions_base,            ONLY : nat, zv, ityp
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(in) :: pdir
  !! direction on which the polarization is calculated
  REAL(DP), INTENT(in) :: e_field
  !! intensity of the field
  REAL(DP), INTENT(inout) :: forces_bp(3,nat)
  !!atomic forces to be update
  !
  ! ... local variables
  !
  INTEGER i
  REAL(DP) :: e !electronic charge (Ry. a.u.)
  REAL(DP) :: a(3), sca
  !
  e = DSQRT(2.d0)
  !
  DO i = 1, nat
     forces_bp(pdir,i) = forces_bp(pdir,i) + e*e_field*zv(ityp(i))
  ENDDO
  !
  RETURN
  !
END SUBROUTINE forces_ion_efield
!
!
!--------------------------------------------------------------------
SUBROUTINE forces_us_efield( forces_bp, pdir, e_field )
   !---------------------------------------------------------------------
   !! It calculates the US correction to the atomic forces 
   !! due to Berry's phase electric field.
   !
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
   USE upf_spinorb,          ONLY : transform_qq_so
   USE lsda_mod,             ONLY : nspin
   USE klist,                ONLY : nelec, degauss, nks, xk, wk, ngk, igk_k
   USE wvfct,                ONLY : npwx, nbnd
   USE wavefunctions,        ONLY : evc
   USE bp,                   ONLY : nppstr_3d, mapgm_global, nx_el, mapg_owner
   USE fixed_occ
   USE mp,                   ONLY : mp_sum, mp_max, mp_barrier
   USE mp_world,             ONLY : world_comm,mpime,nproc
   USE mp_bands,             ONLY : intra_bgrp_comm
   USE becmod,               ONLY : bec_type, becp, calbec,ALLOCATE_bec_type, &
                                    DEALLOCATE_bec_type
   USE noncollin_module,     ONLY : noncolin, npol, lspinorb
   USE parallel_include
   USE uspp_init,            ONLY : init_us_2
   !
   IMPLICIT NONE
   !
   REAL(DP), INTENT(inout) :: forces_bp(3,nat)
   !! atomic forces to be update 
   INTEGER, INTENT(in) :: pdir
   !! direction of electric field
   REAL(DP), INTENT(in) :: e_field
   !! intensity of the field
   !
   ! ... local variables
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
   REAL(dp) :: dk(3)
   REAL(dp) :: dk2
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
   REAL(dp) :: weight
   REAL(dp) :: pola, pola_ion
   REAL(dp), ALLOCATABLE :: wstring(:)
   REAL(dp) :: zeta_mod
   COMPLEX(dp), ALLOCATABLE :: aux(:),aux_2(:)
   COMPLEX(dp), ALLOCATABLE :: aux0(:),aux0_2(:)
   COMPLEX(dp) , ALLOCATABLE :: cphik(:)
   COMPLEX(dp) :: det
   COMPLEX(dp), ALLOCATABLE :: mat(:,:)
   COMPLEX(dp) :: cdet(2)
   COMPLEX(dp) :: cdwork(nbnd)
   !
   COMPLEX(dp) :: pref
   COMPLEX(dp) :: q_dk(nhm,nhm,ntyp)
   COMPLEX(dp) :: struc(nat),struc_r(3,nat)
   COMPLEX(dp) :: zeta
   !
   COMPLEX(dp), ALLOCATABLE :: psi(:,:)
   COMPLEX(dp), ALLOCATABLE :: psi1(:,:)
   COMPLEX(dp) :: zeta_loc
   !
   LOGICAL, ALLOCATABLE :: l_cal(:) ! flag for occupied/empty bands
   INTEGER, ALLOCATABLE :: map_g(:)
   !
   REAL(dp) :: dkfact
   COMPLEX(dp) :: zeta_tot
   !
   COMPLEX(kind=DP) :: sca
   COMPLEX(kind=DP), ALLOCATABLE :: aux_g(:),aux_g_mpi(:,:),aux_proc(:,:),aux_rcv(:,:)
   COMPLEX(DP), ALLOCATABLE :: dbecp0(:,:,:), dbecp_bp(:,:,:),vkb1(:,:)
   INTEGER :: ipol
   COMPLEX(DP) :: forces_tmp(3,nat)
   REAL(DP) :: fact
   TYPE(bec_type) :: becp0, becp_bp
   INTEGER :: nspin_eff
   COMPLEX(DP), ALLOCATABLE :: q_dk_so(:,:,:,:)
   !
   COMPLEX(kind=DP), ALLOCATABLE :: fbmata_1(:,:,:),fbmatb_1(:,:,:,:)
   COMPLEX(kind=DP), ALLOCATABLE :: fbmata_2(:,:,:,:),fbmatb_2(:,:,:,:,:)
   COMPLEX(kind=DP), ALLOCATABLE :: fbmata_3(:,:,:),fbmatb_3(:,:,:,:,:)
   COMPLEX(kind=DP), ALLOCATABLE :: dbecp0_ord(:,:,:,:), dbecp_bp_ord(:,:,:,:)
   !
   INTEGER :: igg,max_aux,max_aux_proc,iproc
   INTEGER, ALLOCATABLE :: aux_g_mpi_ind(:,:),ind_g(:),aux_proc_ind(:,:),aux_rcv_ind(:,:)
   INTEGER :: req, ierr
   !
   !  -------------------------------------------------------------------------   !
   !                               INITIALIZATIONS
   !  -------------------------------------------------------------------------   !
   !
   ALLOCATE( ind_g(nproc) )
   !
   nspin_eff = nspin
   IF (noncolin) THEN
      nspin_eff=1
   ENDIF
   !
   ALLOCATE( psi1(npwx*npol,nbnd)  )
   ALLOCATE( psi(npwx*npol,nbnd)   )
   ALLOCATE( aux(ngm),aux_2(ngm)   )
   ALLOCATE( aux0(ngm),aux0_2(ngm) )
   ALLOCATE( map_g(npwx) )
   ALLOCATE( mat(nbnd,nbnd) )
   ALLOCATE( dbecp0(nkb,nbnd*npol,3), dbecp_bp(nkb,nbnd*npol,3) )
   ALLOCATE( dbecp0_ord(nkb,npol,3,nbnd), dbecp_bp_ord(nkb,npol,3,nbnd) )
   !
   ALLOCATE( vkb1(npwx,nkb) )
   ALLOCATE( l_cal(nbnd) )
   IF (okvan) THEN
      CALL ALLOCATE_bec_type( nkb,nbnd,becp0 )
      CALL ALLOCATE_bec_type( nkb,nbnd,becp_bp )
      IF (lspinorb) ALLOCATE( q_dk_so(nhm,nhm,4,ntyp) )
   ENDIF
   !
   pola = 0.d0 !set to 0 electronic polarization   
   zeta_tot = (1.d0,0.d0)
   !
   !  --- Check that we are working with an insulator with no empty bands ---
   IF ( degauss > 0.0_dp ) CALL errore( 'forces_us_efield', &
            'Polarization only for insulators and no empty bands', 1 )
   !
   !  --- Define a small number ---
   eps = 1.0E-6_dp
   !
   !  --- Recalculate FFT correspondence (see ggen.f90) ---
   ALLOCATE( ln(-dfftp%nr1:dfftp%nr1,-dfftp%nr2:dfftp%nr2,-dfftp%nr3:dfftp%nr3) )
   ln = 0
   DO ng = 1, ngm
      mk1 = NINT(g(1,ng)*at(1,1)+g(2,ng)*at(2,1)+g(3,ng)*at(3,1))
      mk2 = NINT(g(1,ng)*at(1,2)+g(2,ng)*at(2,2)+g(3,ng)*at(3,2))
      mk3 = NINT(g(1,ng)*at(1,3)+g(2,ng)*at(2,3)+g(3,ng)*at(3,3))
      ln(mk1,mk2,mk3) = ng
   ENDDO
   !
   CALL mp_sum( ln, intra_bgrp_comm )
   !
   IF (okvan) THEN
   !  --- Initialize arrays ---
      jkb_bp = 0
      !
      DO nt = 1, ntyp
         DO na = 1, nat
            IF (ityp(na) == nt) THEN
               DO i = 1, nh(nt)
                  jkb_bp = jkb_bp+1
                  nkbtona(jkb_bp) = na
                  nkbtonh(jkb_bp) = i        
               ENDDO
            ENDIF
         ENDDO
      ENDDO
   ENDIF
   !  --- Get the number of strings ---
   nstring = nks/nppstr_3d(pdir)
   nkort = nstring/(nspin_eff)
   !
   !  --- Allocate memory for arrays ---
   ALLOCATE( phik(nstring)  )
   ALLOCATE( loc_k(nstring) )
   ALLOCATE( cphik(nstring) ) 
   ALLOCATE( wstring(nstring)  )
   ALLOCATE( pdl_elec(nstring) )
   ALLOCATE( mod_elec(nstring) )
   !
   FLUSH( stdout )
   !
   !  -------------------------------------------------------------------------   !
   !           electronic polarization: set values for k-points strings           !
   !  -------------------------------------------------------------------------   !
   !
   !  --- Find vector along strings ---
   IF (nppstr_3d(pdir) /= 1) THEN
      gpar(1) = (xk(1,nx_el(nppstr_3d(pdir),pdir))-xk(1,nx_el(1,pdir)))*&
           &DBLE(nppstr_3d(pdir))/DBLE(nppstr_3d(pdir)-1)
      gpar(2) = (xk(2,nx_el(nppstr_3d(pdir),pdir))-xk(2,nx_el(1,pdir)))*&
           &DBLE(nppstr_3d(pdir))/DBLE(nppstr_3d(pdir)-1)
      gpar(3) = (xk(3,nx_el(nppstr_3d(pdir),pdir))-xk(3,nx_el(1,pdir)))*&
           &DBLE(nppstr_3d(pdir))/DBLE(nppstr_3d(pdir)-1)
      gvec = DSQRT(gpar(1)**2+gpar(2)**2+gpar(3)**2)*tpiba
   ELSE
      gpar(1)=0.d0
      gpar(2)=0.d0
      gpar(3)=0.d0
      gpar(pdir)=1.d0/at(pdir,pdir)!
      gvec=tpiba/SQRT(at(pdir,1)**2.d0+at(pdir,2)**2.d0+at(pdir,3)**2.d0)
   ENDIF      
   !
   !  --- Find vector between consecutive points in strings ---
   IF (nppstr_3d(pdir) /= 1) THEN  ! orthorhombic cell 
      dk(1) = xk(1,nx_el(2,pdir))-xk(1,nx_el(1,pdir))
      dk(2) = xk(2,nx_el(2,pdir))-xk(2,nx_el(1,pdir))
      dk(3) = xk(3,nx_el(2,pdir))-xk(3,nx_el(1,pdir))
      dkmod = SQRT(dk(1)**2+dk(2)**2+dk(3)**2)*tpiba 
   ELSE ! Gamma point case, only cubic cell for now
      dk(1)=0.d0
      dk(2)=0.d0
      dk(3)=0.d0
      dk(pdir)=1.d0/at(pdir,pdir)
      dkmod=tpiba/SQRT(at(pdir,1)**2.d0+at(pdir,2)**2.d0+at(pdir,3)**2.d0)
   ENDIF
   !
   !  -------------------------------------------------------------------------   !
   !                   electronic polarization: weight strings                    !
   !  -------------------------------------------------------------------------   !
   !
   !  --- Calculate string weights, normalizing to 1 (no spin) or 1+1 (spin) ---
   DO is = 1, nspin_eff
      weight = 0.0_dp
      DO kort = 1, nkort
         istring = kort+(is-1)*nkort
         wstring(istring) = wk(nppstr_3d(pdir)*istring)
         weight = weight+wstring(istring)
      ENDDO
      !
      DO kort = 1, nkort
         istring = kort+(is-1)*nkort
         wstring(istring) = wstring(istring)/weight
      ENDDO
      !
   ENDDO  
   !
   !  -------------------------------------------------------------------------   !
   !                  electronic polarization: structure factor                   !
   !  -------------------------------------------------------------------------   !
   !
   !  --- Calculate structure factor e^{-i dk*R} ---
   !
   DO na = 1, nat
      fac = (dk(1)*tau(1,na)+dk(2)*tau(2,na)+dk(3)*tau(3,na))*tpi 
      struc(na) = CMPLX(COS(fac),-SIN(fac),KIND=DP)
   ENDDO
   !
   ! Calculate derivatives of structure factors
   DO na = 1, nat
      DO ipol = 1, 3
         struc_r(ipol,na) = struc(na)*CMPLX(0.d0,-1.d0, KIND=DP)*dk(ipol)
      ENDDO
   ENDDO
   !
   !  -------------------------------------------------------------------------   !
   !                     electronic polarization: form factor                     !
   !  -------------------------------------------------------------------------   !
   IF (okvan) THEN
      !  --- Bessel transform of Q_ij(|r|) at dk [Q_ij^L(|r|)] in array qrad ---
      ! CALL calc_btq( dkmod, qrad_dk, 0 ) no longer needed
      CALL compute_qqc ( tpiba, dk, omega, q_dk )
      IF (lspinorb) CALL transform_qq_so( q_dk, q_dk_so )
   ENDIF
   !
   ! calculate factor
   !
   CALL factor_a( pdir, at, dkfact )
   fact = DSQRT(2.d0)*e_field*dkfact
   !
   IF (nspin_eff==1 .AND. .NOT.noncolin) fact = fact*2.d0
   !
   !  -------------------------------------------------------------------------   !
   !                   electronic polarization: strings phases                    !
   !  -------------------------------------------------------------------------   !
   !
   el_loc = 0.d0
   kpoint = 0
   zeta = (1.d0,0.d0)
   !  --- Start loop over spin ---
   DO is = 1, nspin_eff 
      ! l_cal(n) = .true./.false. if n-th state is occupied/empty
      DO nb = 1, nbnd
         IF (nspin_eff == 2 .AND. tfixed_occ) THEN
            l_cal(nb) = ( f_inp(nb,is) /= 0.0_dp )
         ELSE
            IF(noncolin) THEN
               l_cal(nb) = ( nb <= NINT( nelec ) )
            ELSE
               l_cal(nb) = ( nb <= NINT( nelec/2.0_dp ) )
            ENDIF
         ENDIF
      ENDDO
      !
      ! --- Start loop over orthogonal k-points ---
      DO kort = 1, nkort
         zeta_loc = (1.d0,0.d0)
         !  --- Index for this string ---
         istring = kort+(is-1)*nkort
         !
         !  --- Initialize expectation value of the phase operator ---
         !
         zeta_mod = 1.d0
         !
         !  --- Start loop over parallel k-points ---
         DO kpar = 1,nppstr_3d(pdir)+1
            !
            ! --- Set index of k-point ---
            kpoint = kpoint + 1
            !
            ! --- Calculate dot products between wavefunctions and betas ---
            IF (kpar /= 1 ) THEN
               !
               ! --- Dot wavefunctions and betas for PREVIOUS k-point ---
               ik = nx_el(kpoint-1,pdir)
               npw0   = ngk(ik)
               igk0(:)= igk_k(:,ik)
               !
               CALL get_buffer( psi, nwordwfc, iunwfc, nx_el(kpoint-1,pdir) )
               !
               IF (okvan) THEN
                  CALL init_us_2( npw0, igk0, xk(1,nx_el(kpoint-1,pdir)), vkb )
                  CALL calbec( npw0, vkb, psi, becp0 )
                  DO ipol = 1, 3
                     DO jkb = 1, nkb
                        DO ig = 1, npw0
                           vkb1(ig,jkb) = vkb(ig,jkb)*(0.D0,-1.D0)*g(ipol,igk0(ig))
                        ENDDO
                     ENDDO
                     !
                     IF ( nkb > 0 ) &
                          CALL ZGEMM( 'C', 'N', nkb, nbnd*npol, npw0, ( 1.D0, 0.D0 ), &
                                      vkb1, npwx, psi, npwx, ( 0.D0, 0.D0 ),          &
                                      dbecp0(1,1,ipol), nkb )
                     CALL mp_sum( dbecp0(1:nkb,1:nbnd*npol,ipol), intra_bgrp_comm )
                  ENDDO
               ENDIF
               !
               ! --- Dot wavefunctions and betas for CURRENT k-point ---
               IF (kpar /= (nppstr_3d(pdir)+1)) THEN
                  !
                  ik = nx_el(kpoint,pdir)
                  npw1   = ngk(ik)
                  igk1(:)= igk_k(:,ik)
                  !
                  CALL get_buffer( psi1, nwordwfc, iunwfc, nx_el(kpoint,pdir) )
                  !
                  IF (okvan) THEN
                     CALL init_us_2 (npw1,igk1,xk(1,nx_el(kpoint,pdir)),vkb)
                     CALL calbec( npw1, vkb, psi1, becp_bp)
                     DO ipol = 1, 3
                        DO jkb = 1, nkb
                           DO ig = 1, npw1
                              vkb1(ig,jkb) = vkb(ig,jkb)*(0.D0,-1.D0)*g(ipol,igk1(ig))
                           ENDDO
                        ENDDO
                        !
                        IF ( nkb > 0 ) &
                           CALL ZGEMM( 'C', 'N', nkb, nbnd*npol, npw1, ( 1.D0, 0.D0 ), &
                                       vkb1, npwx, psi1, npwx, ( 0.D0, 0.D0 ),         &
                                       dbecp_bp(1,1,ipol), nkb )
                        CALL mp_sum( dbecp_bp(1:nkb,1:nbnd*npol,ipol), intra_bgrp_comm )
                     ENDDO
                  ENDIF
                  !
               ELSE
                  !
                  kstart = kpoint-(nppstr_3d(pdir)+1)+1
                  ik = nx_el(kstart,pdir)
                  npw1 = ngk(ik)
                  igk1(:) = igk_k(:,ik)
                  !
                  CALL get_buffer( psi1, nwordwfc, iunwfc, nx_el(kstart,pdir) )
                  !
                  IF (okvan) THEN
                     CALL init_us_2( npw1, igk1, xk(1,nx_el(kstart,pdir)), vkb )
                     CALL calbec( npw1, vkb, psi1, becp_bp )
                     DO ipol = 1, 3
                        DO jkb = 1, nkb
                           DO ig = 1, npw1
                              vkb1(ig,jkb) = vkb(ig,jkb)*(0.D0,-1.D0)*g(ipol,igk1(ig))
                           ENDDO
                        ENDDO
                        IF ( nkb > 0 ) &
                             CALL ZGEMM( 'C', 'N', nkb, nbnd*npol, npw1, ( 1.D0, 0.D0 ), &
                                         vkb1, npwx, psi1, npwx, ( 0.D0, 0.D0 ),         &
                                         dbecp_bp(1,1,ipol), nkb )
                        CALL mp_sum( dbecp_bp(1:nkb,1:nbnd*npol,ipol), intra_bgrp_comm )
                     ENDDO
                  ENDIF
                  !
               ENDIF
               !
               ! --- Matrix elements calculation ---
               !
               mat = (0.d0,0.d0)
               DO nb = 1, nbnd
                  aux = (0.d0,0.d0)
                  aux0 = (0.d0,0.d0)
                  IF (noncolin) THEN
                     aux_2 = (0.d0,0.d0)
                     aux0_2 = (0.d0,0.d0)
                  ENDIF
                  !
                  DO ig = 1, npw0
                     aux0(igk0(ig)) = psi(ig,nb)
                  ENDDO
                  !
                  IF (noncolin) THEN
                     DO ig = 1, npw0
                        aux0_2(igk0(ig)) = psi(ig+npwx,nb)
                     ENDDO
                  ENDIF
                  !
                  IF (kpar /= (nppstr_3d(pdir)+1)) THEN
                     DO mb = 1, nbnd
                        IF ( l_cal(nb) .AND. l_cal(mb) ) THEN
                           DO ig = 1, npw1
                              aux(igk1(ig)) = psi1(ig,mb)
                           ENDDO
                           !
                           IF(noncolin) THEN
                              DO ig = 1, npw1
                                 aux_2(igk1(ig)) = psi1(ig+npwx,mb)
                              ENDDO
                           ENDIF
                           !
                           mat(nb,mb) = dot_product(aux0(1:ngm),aux(1:ngm))
                           !
                           IF (noncolin) THEN
                                   mat(nb,mb) = mat(nb,mb) + dot_product(aux0_2(1:ngm),aux_2(1:ngm))
                           ENDIF
                           !
                           CALL mp_sum( mat(nb,mb), intra_bgrp_comm )
                           !
                        ENDIF
                     ENDDO
                  ENDIF
               ENDDO
               !
               IF ( kpar == (nppstr_3d(pdir)+1) ) THEN
                  !
                  ! ALLOCATE global array
                  !   
                  ALLOCATE( aux_g(ngm_g), aux_g_mpi(ngmx,nproc), aux_g_mpi_ind(ngmx,nproc) )
                  !
                  DO ipol = 0, npol-1
                     DO mb = 1, nbnd
                        !
                        aux_g_mpi = 0.d0
                        aux_g_mpi_ind = 0
                        ind_g = 0
                        !
                        DO ig = 1, npw1
                           igg = mapgm_global(ig_l2g(igk1(ig)),pdir)
                           ind_g(mapg_owner(1,igg)) = ind_g(mapg_owner(1,igg))+1
                           aux_g_mpi(ind_g(mapg_owner(1,igg)),mapg_owner(1,igg)) = psi1(ig+npwx*ipol,mb)
                           aux_g_mpi_ind(ind_g(mapg_owner(1,igg)),mapg_owner(1,igg)) = mapg_owner(2,igg)
                        ENDDO
                        !
                        max_aux = 0
                        !
                        DO iproc=1,nproc
                           IF (iproc/=mpime+1) THEN
                              max_aux_proc=0
                              DO ig=1,ngmx
                                 IF (aux_g_mpi_ind(ig,iproc) > 0) THEN
                                    max_aux_proc=max_aux_proc+1
                                 ELSE
                                    exit
                                 ENDIF
                              ENDDO
                              IF (max_aux_proc>max_aux) max_aux = max_aux_proc
                           ENDIF
                        ENDDO
                        !
                        CALL mp_max(max_aux, intra_bgrp_comm )
                        ALLOCATE( aux_proc(max_aux,nproc), aux_proc_ind(max_aux,nproc) )
                        ALLOCATE( aux_rcv(max_aux,nproc), aux_rcv_ind(max_aux,nproc) )
                        aux_proc = (0.d0,0.d0)
                        aux_proc_ind = 0
                        DO iproc = 1, nproc
                           IF (iproc /= mpime+1) THEN
                              DO ig = 1, max_aux
                                 IF (aux_g_mpi_ind(ig,iproc) > 0) THEN
                                    aux_proc(ig,iproc) = aux_g_mpi(ig,iproc)
                                    aux_proc_ind(ig,iproc) = aux_g_mpi_ind(ig,iproc)
                                 ELSE
                                    EXIT
                                 ENDIF
                              ENDDO
                           ENDIF
                        ENDDO
                        !
#if defined (__MPI)
                        CALL MPI_ALLTOALL( aux_proc, max_aux, MPI_DOUBLE_COMPLEX,  &
                             aux_rcv, max_aux, MPI_DOUBLE_COMPLEX, intra_bgrp_comm, ierr )
                        CALL MPI_ALLTOALL( aux_proc_ind, max_aux, MPI_INTEGER,     &
                             aux_rcv_ind, max_aux, MPI_INTEGER, intra_bgrp_comm, ierr )
#else
                        aux_rcv(1:max_aux,1) = aux_proc(1:max_aux,1)
                        aux_rcv_ind(1:max_aux,1) = aux_proc_ind(1:max_aux,1)
#endif
                        DO nb = 1, nbnd
                           IF ( l_cal(nb) .AND. l_cal(mb) ) THEN
                              aux = (0.d0,0.d0)
                              aux0 = (0.d0,0.d0)
                              IF (noncolin) THEN
                                 aux_2 = (0.d0,0.d0)
                                 aux0_2 = (0.d0,0.d0)
                              ENDIF
                              !
                              DO ig = 1, npw0
                                 aux0(igk0(ig)) = psi(ig,nb)
                              ENDDO
                              !
                              IF (noncolin) THEN
                                 DO ig = 1, npw0
                                    aux0_2(igk0(ig)) = psi(ig+npwx,nb)
                                 ENDDO
                              ENDIF
                              !
                              sca = 0.d0
                              !
                              DO iproc = 1, nproc
                                 IF (iproc /= mpime+1) THEN
                                    DO ig = 1, max_aux
                                       IF (aux_rcv_ind(ig,iproc)/=0) THEN
                                          IF (aux_rcv_ind(ig,iproc)<0 .OR. aux_rcv_ind(ig,iproc)> ngm) THEN
                                             WRITE(stdout,*) 'OH BOY', aux_rcv_ind(ig,iproc)
                                          ELSE
                                             IF (ipol == 0) THEN
                                                sca = sca+CONJG(aux0(aux_rcv_ind(ig,iproc)))*aux_rcv(ig,iproc)
                                             ELSE
                                                sca = sca+CONJG(aux0_2(aux_rcv_ind(ig,iproc)))*aux_rcv(ig,iproc)
                                             ENDIF
                                          ENDIF
                                       ELSE
                                          EXIT
                                       ENDIF
                                    ENDDO
                                 ENDIF
                              ENDDO
                              !
                              DO ig = 1, ngmx
                                 IF (aux_g_mpi_ind(ig,mpime+1) /= 0) THEN
                                    IF (aux_g_mpi_ind(ig,mpime+1)<0 .OR. aux_g_mpi_ind(ig,mpime+1)>ngm) THEN 
                                       WRITE(stdout,*) 'OH BOY', aux_g_mpi_ind(ig,mpime+1)
                                    ELSE
                                       IF (ipol==0) THEN
                                          sca = sca+CONJG(aux0(aux_g_mpi_ind(ig,mpime+1)))*aux_g_mpi(ig,mpime+1)
                                       ELSE
                                          sca = sca+CONJG(aux0_2(aux_g_mpi_ind(ig,mpime+1)))*aux_g_mpi(ig,mpime+1)
                                       ENDIF
                                    ENDIF
                                 ELSE
                                    EXIT
                                 ENDIF
                              ENDDO
                              !                                 
                              CALL mp_sum( sca, intra_bgrp_comm )
                              mat(nb,mb) = mat(nb,mb)+sca
                           ENDIF
                        ENDDO
                        !
                        DEALLOCATE( aux_proc, aux_proc_ind )
                        DEALLOCATE( aux_rcv, aux_rcv_ind )
                        !
                     ENDDO
                  ENDDO
                  !
                  DEALLOCATE( aux_g, aux_g_mpi, aux_g_mpi_ind )                  
                  !
               ENDIF
               !
               DO nb = 1, nbnd
                  DO mb = 1, nbnd 
                     !
                     IF ( l_cal(nb) .AND. l_cal(mb) ) THEN
                        !
                        ! --- Calculate the augmented part: ij=KB projectors, ---
                        ! --- R=atom index: SUM_{ijR} q(ijR) <u_nk|beta_iR>   ---
                        ! --- <beta_jR|u_mk'> e^i(k-k')*R =                   ---
                        ! --- also <u_nk|beta_iR>=<psi_nk|beta_iR> = becp^*   ---
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
                                    pref = pref+CONJG(becp0%nc(jkb,1,nb))*becp_bp%nc(jkb1+j,1,mb) &
                                         *q_dk_so(nhjkb,j,1,np)*struc(na)
                                    pref = pref+CONJG(becp0%nc(jkb,1,nb))*becp_bp%nc(jkb1+j,2,mb) &
                                         *q_dk_so(nhjkb,j,2,np)*struc(na)
                                    pref = pref+CONJG(becp0%nc(jkb,2,nb))*becp_bp%nc(jkb1+j,1,mb) &
                                         *q_dk_so(nhjkb,j,3,np)*struc(na)
                                    pref = pref+CONJG(becp0%nc(jkb,2,nb))*becp_bp%nc(jkb1+j,2,mb) &
                                         *q_dk_so(nhjkb,j,4,np)*struc(na)

                                 ELSE
                                    pref = pref+CONJG(becp0%k(jkb,nb))*becp_bp%k(jkb1+j,mb) &
                                         *q_dk(nhjkb,j,np)*struc(na)
                                 ENDIF
                              ENDDO
                           ENDDO
                           !
                           mat(nb,mb) = mat(nb,mb) + pref
                        ENDIF
                     ENDIF !on l_cal
                     !
                  ENDDO
               ENDDO
               DO nb = 1, nbnd
                  IF ( .NOT. l_cal(nb) ) mat(nb,nb) = 1.d0
               END DO
               !
               ! --- Calculate matrix determinant ---
               !
               ! calculate inverse
               !              
               CALL zgefa( mat, nbnd, nbnd, ivpt, info )
               CALL errore( 'forces_us_efield', 'error in zgefa', ABS(info) )
               CALL zgedi( mat, nbnd, nbnd, ivpt, cdet, cdwork, 1 )
               !
               DO nb = 1, nbnd
                  DO mb = 1, nbnd
                     IF (.NOT.l_cal(nb) .OR. .NOT. l_cal(mb)) mat(mb,nb)=(0.d0,0.d0)
                  ENDDO
               ENDDO
               !
               ! calculate terms
               !
               forces_tmp(:,:) = (0.d0,0.d0)
               !
               IF (okvan) THEN
                  ALLOCATE( fbmatb_1(nkb,npol,nkb,npol),   fbmata_1(nbnd,nkb,npol)   )
                  ALLOCATE( fbmatb_2(nkb,npol,nkb,npol,3), fbmata_2(nbnd,nkb,npol,3) )
                  ALLOCATE( fbmatb_3(nkb,npol,3,nkb,npol), fbmata_3(nbnd,nkb,npol)   )
                  !
                  IF (lspinorb) THEN
                     CALL ZGEMM( 'N','C',nbnd,nkb*npol,nbnd,(1.d0,0.d0), &
                                 mat,nbnd,becp0%nc(1,1,1),nkb*npol,(0.d0,0.d0),fbmata_1,nbnd )
                     CALL ZGEMM( 'N','N',nkb*npol,nkb*npol,nbnd,(1.d0,0.d0), &
                                 becp_bp%nc(1,1,1),nkb*npol,fbmata_1,nbnd,(0.d0,0.d0),fbmatb_1,nkb*npol )
                     !
                     DO ipol = 1, npol
                        DO nb = 1, nbnd
                           dbecp0_ord(1:nkb,ipol,1:3,nb) = dbecp0(1:nkb,(nb-1)*npol+ipol,1:3)
                           dbecp_bp_ord(1:nkb,ipol,1:3,nb) = dbecp_bp(1:nkb,(nb-1)*npol+ipol,1:3)
                        ENDDO
                     ENDDO
                     !
                     CALL ZGEMM( 'N','C',nbnd,nkb*npol*3,nbnd,(1.d0,0.d0), &
                                 mat,nbnd,dbecp0_ord,nkb*npol*3,(0.d0,0.d0),fbmata_2,nbnd )
                     CALL ZGEMM( 'N','N',nkb*npol,nkb*npol*3,nbnd,(1.d0,0.d0), &
                                 becp_bp%nc(1,1,1),nkb*npol,fbmata_2,nbnd,(0.d0,0.d0),fbmatb_2,nkb*npol )
                     !
                     CALL ZGEMM( 'N','C',nbnd,nkb*npol,nbnd,(1.d0,0.d0), &
                                 mat,nbnd,becp0%nc(1,1,1),nkb*npol,(0.d0,0.d0),fbmata_3,nbnd )
                     CALL ZGEMM( 'N','N',nkb*npol*3,nkb*npol,nbnd,(1.d0,0.d0), &
                                 dbecp_bp_ord,nkb*npol*3,fbmata_3,nbnd,(0.d0,0.d0),fbmatb_3,nkb*npol*3 )
                     !
                     DO jkb = 1, nkb
                        nhjkb = nkbtonh(jkb)
                        na = nkbtona(jkb)
                        np = ityp(na)
                        nhjkbm = nh(np)
                        jkb1 = jkb - nhjkb
                        DO j = 1,nhjkbm
                           forces_tmp(1:3,na) = forces_tmp(1:3,na)+ &   
                                q_dk_so(nhjkb,j,1,np)*struc_r(1:3,na)*fbmatb_1(jkb1+j,1,jkb,1)   
                           forces_tmp(1:3,na) = forces_tmp(1:3,na)+ &   
                                q_dk_so(nhjkb,j,1,np)*struc(na)*fbmatb_2(jkb1+j,1,jkb,1,1:3)   
                           forces_tmp(1:3,na) = forces_tmp(1:3,na)+ &   
                                q_dk_so(nhjkb,j,1,np)*struc(na)*fbmatb_3(jkb1+j,1,1:3,jkb,1)   
                           !   
                           forces_tmp(1:3,na) = forces_tmp(1:3,na)+ &   
                                q_dk_so(nhjkb,j,2,np)*struc_r(1:3,na)*fbmatb_1(jkb1+j,2,jkb,1)   
                           forces_tmp(1:3,na) = forces_tmp(1:3,na)+ &   
                                q_dk_so(nhjkb,j,2,np)*struc(na)*fbmatb_2(jkb1+j,2,jkb,1,1:3)   
                           forces_tmp(1:3,na) = forces_tmp(1:3,na)+ &   
                                q_dk_so(nhjkb,j,2,np)*struc(na)*fbmatb_3(jkb1+j,2,1:3,jkb,1)   
                           !   
                           forces_tmp(1:3,na) = forces_tmp(1:3,na)+ &   
                                q_dk_so(nhjkb,j,3,np)*struc_r(1:3,na)*fbmatb_1(jkb1+j,1,jkb,2)   
                           forces_tmp(1:3,na) = forces_tmp(1:3,na)+ &   
                                q_dk_so(nhjkb,j,3,np)*struc(na)*fbmatb_2(jkb1+j,1,jkb,2,1:3)   
                           forces_tmp(1:3,na) = forces_tmp(1:3,na)+ &   
                                q_dk_so(nhjkb,j,3,np)*struc(na)*fbmatb_3(jkb1+j,1,1:3,jkb,2)   
                           !   
                           forces_tmp(1:3,na) = forces_tmp(1:3,na)+ &   
                                q_dk_so(nhjkb,j,4,np)*struc_r(1:3,na)*fbmatb_1(jkb1+j,2,jkb,2)   
                           forces_tmp(1:3,na) = forces_tmp(1:3,na)+ &   
                                q_dk_so(nhjkb,j,4,np)*struc(na)*fbmatb_2(jkb1+j,2,jkb,2,1:3)   
                           forces_tmp(1:3,na) = forces_tmp(1:3,na)+ &   
                                q_dk_so(nhjkb,j,4,np)*struc(na)*fbmatb_3(jkb1+j,2,1:3,jkb,2)   
                           !   
                        ENDDO
                     ENDDO
                     !
                  ENDIF
                  !
                  IF (.NOT.lspinorb) THEN
                     !
                     DO jkb = 1, nkb
                        nhjkb = nkbtonh(jkb)
                        na = nkbtona(jkb)
                        np = ityp(na)
                        nhjkbm = nh(np)
                        jkb1 = jkb - nhjkb
                        !
                        DO j = 1,nhjkbm
                           !
                           DO nb = 1, nbnd
                              DO mb = 1, nbnd
                                 IF (lspinorb) THEN
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

                                 ELSE
                                    forces_tmp(:,na) = forces_tmp(:,na)+CONJG(becp0%k(jkb,nb))*becp_bp%k(jkb1+j,mb) &
                                         *q_dk(nhjkb,j,np)*struc_r(:,na)*mat(mb,nb)
                                    forces_tmp(:,na) = forces_tmp(:,na)+CONJG(dbecp0(jkb,nb,:))*becp_bp%k(jkb1+j,mb) &
                                         *q_dk(nhjkb,j,np)*struc(na)*mat(mb,nb)
                                    forces_tmp(:,na) = forces_tmp(:,na)+CONJG(becp0%k(jkb,nb))*dbecp_bp(jkb1+j,mb,:) &
                                         *q_dk(nhjkb,j,np)*struc(na)*mat(mb,nb)
                                 ENDIF
                              ENDDO
                           ENDDO
                           !
                        ENDDO
                        !
                     ENDDO
                     !
                  ENDIF
                  !
                  DEALLOCATE( fbmata_1, fbmatb_1 )
                  DEALLOCATE( fbmata_2, fbmatb_2 )
                  DEALLOCATE( fbmata_3, fbmatb_3 )
                  !
               ENDIF
               !
               forces_bp(:,:) = forces_bp(:,:)+fact*AIMAG(forces_tmp(:,:))*wstring(istring)
               !
               ! --- End of dot products between wavefunctions and betas ---
            ENDIF
            !
            ! --- End loop over parallel k-points ---
            !
         ENDDO
         !
         kpoint = kpoint-1
         ! --- End loop over orthogonal k-points ---
      ENDDO
      !
      ! --- End loop over spin ---
   ENDDO
   !
   !  --- Free memory ---
   DEALLOCATE( l_cal )
   DEALLOCATE( pdl_elec )
   DEALLOCATE( mod_elec )
   DEALLOCATE( wstring )
   DEALLOCATE( loc_k   )
   DEALLOCATE( phik  )
   DEALLOCATE( cphik )
   DEALLOCATE( ln    )
   DEALLOCATE( map_g )
   DEALLOCATE( aux, aux_2   )
   DEALLOCATE( aux0, aux0_2 )
   DEALLOCATE( psi  )
   DEALLOCATE( psi1 )
   DEALLOCATE( mat  )
   IF (okvan) THEN
      CALL deallocate_bec_type( becp0 )
      CALL deallocate_bec_type( becp_bp )
      IF (lspinorb) DEALLOCATE( q_dk_so )
   ENDIF
   DEALLOCATE( dbecp0, dbecp0_ord, dbecp_bp, dbecp_bp_ord )
   !
   DEALLOCATE( ind_g )
   !
END SUBROUTINE forces_us_efield
!
!
!---------------------------------------------------------------------------------
SUBROUTINE stress_bp_efield( sigmael )
   !-----------------------------------------------------------------------------
   !! Calculate the stress contribution due to the electric field
   !! electronic part.
   !
   USE kinds,                ONLY : DP
   USE bp,                   ONLY : efield_cart, el_pol, fc_pol, l3dstring
   USE cell_base,            ONLY : at, alat, tpiba, omega
   USE constants,            ONLY : pi
   !
   IMPLICIT NONE
   !
   REAL(DP), INTENT(out) :: sigmael(3,3)
   !! stress contribution to be calculated
   !
   REAL(DP) :: phases(3)
   INTEGER :: i, j, ipol
   !
   sigmael(:,:) = 0.d0
   IF (.NOT.l3dstring ) RETURN
   phases(:) = el_pol(:)/fc_pol(:)
   DO ipol = 1, 3
      DO i = 1, 3
         DO j = 1, 3
            sigmael(i,j) = sigmael(i,j)-efield_cart(i)*at(j,ipol)*phases(ipol)
         ENDDO
      ENDDO
   ENDDO
   !
   sigmael(:,:) = sigmael(:,:)*alat*DSQRT(2.d0)/(2.d0*pi)/omega
   !
   RETURN
   !
END SUBROUTINE stress_bp_efield
!
!
!----------------------------------------------------------------------------------
SUBROUTINE stress_ion_efield( sigmaion )
   !-------------------------------------------------------------------------------
   !! Calculate the stress contribution due to the  electric field ionic part.
   !
   USE kinds,                ONLY : DP
   USE bp,                   ONLY : efield_cart, ion_pol,l3dstring
   USE cell_base,            ONLY : at, alat, omega, bg
   USE constants,            ONLY : pi
   !
   IMPLICIT NONE
   !
   REAL(DP), INTENT(out) :: sigmaion(3,3)!stress contribution to be calculated   
   REAL(DP) :: pol_cry(3)
   INTEGER :: i,j,ipol
   !
   sigmaion(:,:) = 0.d0
   IF (.NOT.l3dstring ) RETURN
   !
   pol_cry(:) = ion_pol(:)
   !
   CALL cryst_to_cart( 1, pol_cry, at, -1 )
   !
   DO ipol = 1, 3
      DO i = 1, 3
         DO j = 1, 3
            sigmaion(i,j) = sigmaion(i,j)-efield_cart(i)*at(j,ipol)*pol_cry(ipol)
         ENDDO
      ENDDO
   ENDDO
   !
   sigmaion(:,:) = sigmaion(:,:)/omega
   !
   RETURN
   !
END SUBROUTINE stress_ion_efield
