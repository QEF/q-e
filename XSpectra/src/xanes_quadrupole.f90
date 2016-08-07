
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!  XANES calculation in the electric quadrupole approximation
!------------------------------------------------------------------------------
SUBROUTINE xanes_quadrupole(a,b,ncalcv,xnorm,core_wfn,paw_iltonhb,&
                            terminator,verbosity)
  !----------------------------------------------------------------------------
  USE io_global,       ONLY: stdout     ! Modules/io_global.f90
  USE kinds,           ONLY: DP
  USE constants,       ONLY: pi
  USE parameters,      ONLY: ntypx
  USE radial_grids,    ONLY: ndmx
  USE ions_base,       ONLY: nat, ntyp => nsp, ityp
  USE wvfct,           ONLY: npwx, nbnd, et, current_k
  USE gvecw,           ONLY: gcutw
  USE lsda_mod,        ONLY: nspin,lsda,isk,current_spin
  USE cell_base,       ONLY: tpiba2, bg
  USE wavefunctions_module, ONLY: evc
  USE klist,           ONLY: &
       nkstot,               & ! total number of k-points
       nks,                  & ! number of k-points per pool
       xk,                   & ! k-points coordinates
       wk,                   & ! k-points weight
       ngk, igk_k
  USE gvect,           ONLY: g,ngm,ngl
  USE fft_base,        ONLY: dfftp
  USE paw_gipaw,       ONLY: &
       paw_vkb,              & ! |p> projectors
       paw_becp,             & ! product of projectors and wf.
       paw_nkb,              & ! total number of beta functions, with st.fact.
       paw_lmaxkb,           & 
       paw_recon
  USE becmod,          ONLY: becp, allocate_bec_type, deallocate_bec_type ! CG
  USE scf,             ONLY: vltot, v, vrs, kedtau !CG
  USE gvecs,           ONLY: doublegrid
  USE mp_pools,        ONLY: intra_pool_comm, root_pool, npool
  USE mp_world,        ONLY: world_comm
  USE mp,              ONLY: mp_sum,mp_barrier, mp_bcast !CG
  USE xspectra,        ONLY: xiabs, xanes_qua, xang_mom, xniter, xnitermax,&
                             xkvec, xepsilon, save_file_kind,              &
                             calculated, time_limit
  USE atom,            ONLY: rgrid, msh
  !  use atom,        ONLY : &
  !       mesh,     &!mesh(ntypx) number of mesh points
  !       msh ,     &!msh(ntypx)the point at rcut=end of radial integration
  !       r   
  USE radin_mod
  USE uspp,            ONLY: vkb, nkb, okvan !CG
  USE ldaU,            ONLY: lda_plus_u
  USE basis,           ONLY: natomwfc
  !<CG>
  USE xspectra_paw_variables, ONLY: xspectra_paw_nhm
  !</CG>

  IMPLICIT NONE
  !
  REAL(dp), INTENT(INOUT) :: a(xnitermax,1,nks)
  REAL(dp), INTENT(INOUT) :: b(xnitermax,1,nks)     
  REAL(dp), INTENT(INOUT) :: xnorm(1,nks)
  REAL(dp), INTENT(IN)    :: core_wfn(ndmx)
  INTEGER, INTENT(INOUT)  :: ncalcv(1,nks)
  INTEGER, INTENT(IN)     :: paw_iltonhb(0:paw_lmaxkb,xspectra_paw_nhm, ntyp)
  LOGICAL, INTENT(IN)     :: terminator
  !
  !... Local variables
  !
  INTEGER :: is,ik,iabso,nr,ip,jp,l,j,icrd,ip_l,nrc,nt,na
  INTEGER :: ipx,ipx_0,ipy,ipz,nline,nrest,npw, npw_partial
  INTEGER :: nunfinished
  LOGICAL :: recalc
  REAL (dp) :: pref,prefb,v_of_0,xnorm_partial,prefm2,prefm1,prefm0
  REAL (dp) :: norm,normps


  REAL (dp), ALLOCATABLE :: aux(:)
  COMPLEX(KIND=DP), EXTERNAL :: zdotc
  COMPLEX(dp), ALLOCATABLE :: paw_vkb_cplx(:,:)
  COMPLEX(KIND=dp), ALLOCATABLE :: psi(:)
  COMPLEX(dp), ALLOCATABLE :: psiwfc(:), spsiwfc(:)
  CHARACTER(LEN=4) :: verbosity

  REAL(DP), EXTERNAL ::  get_clock
  REAL(dp) :: timenow
  EXTERNAL zdscal

  nunfinished=0
  timenow=0

  ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  !  Constant Definitions
  ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$


  prefm2=SQRT(3.0/40.0)/3.0
  prefm1=prefm2
  prefm0=2.0*SQRT(1.0/40.0)/3.0

  pref=SQRT(2.d0)
  prefb=1.0/pref


  ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  !  Variable allocation and initialization
  ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ 

  ALLOCATE(aux(rgrid(xiabs)%mesh))!overdimensionated, necessary only up to msh
  ALLOCATE(psi(npwx))
  ALLOCATE (paw_vkb_cplx( npwx,  paw_nkb))
  ALLOCATE(xanes_qua(paw_recon(xiabs)%paw_nl(xang_mom)))
  ALLOCATE(psiwfc(npwx))
  xanes_qua(:)=0.d0
  psi(:)=(0.d0)

  ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  ! Radial Quadrupole Matrix element
  ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  !  Radial part :   <\phi_n|r|\psi_i>=int r^3 \phi_n \psi_i dr.
  !  Here only ps and ae wavefunctions from XX.recon
  !

  

  !WRITE( stdout,*) 'Calculation Quadrupole matrix element'
  !WRITE( stdout,*) 'There are ',paw_recon(xiabs)%paw_nl(xang_mom),'
  !projectors/channel'
  !WRITE( stdout,*) 'xang_mom=',xang_mom,' xiabs=',xiabs

  ! ... Checks that the core wf is correctly normalized

  nr=msh(xiabs)  ! extended up to all the NON ZERO points in the mesh.
  IF(TRIM(verbosity).EQ.'high') THEN
     aux(1:nr)=core_wfn(1:nr)*core_wfn(1:nr)
     WRITE (stdout,'(8x,"Norm of core wfc =",f10.6)') &
            SQRT(para_radin(aux(1:nr),rgrid(xiabs)%r(1:nr),nr))
  ENDIF

  ! ... Calculate the radial integral

  ip_l=0

  DO ip=1,paw_recon(xiabs)%paw_nbeta
     IF(paw_recon(xiabs)%psphi(ip)%label%l.EQ.xang_mom) THEN
        nrc=paw_recon(xiabs)%aephi(ip)%label%nrc
        ip_l=ip_l+1
        !    here below, recall that psi is r*psi and you have a Jacobian=r^2
        !        aux(1:nr)=r(1:nr,xiabs)*r(1:nr,xiabs)* &
        !             psphi(xiabs,ip)%psi(1:nr)*core_wfn(1:nr)
        aux(1:nrc)=rgrid(xiabs)%r(1:nrc)*rgrid(xiabs)%r(1:nrc)*&
                   paw_recon(xiabs)%aephi(ip)%psi(1:nrc)*core_wfn(1:nrc)
        !    here we have to integrate only inside the augmentation region.
        xanes_qua(ip_l)=para_radin(aux(1:nrc),rgrid(xiabs)%r(1:nrc),nrc)
     ENDIF
  ENDDO
 
  ! ... Write the radial transition matrix element(s)


  DO ip=1,paw_recon(xiabs)%paw_nl(xang_mom) !
     WRITE(stdout,'(8x,a,i1,a,i1,a,f14.9)') &
            '| For PAW proj. (l=',xang_mom,') #',&
            ip, ': radial matrix element =', xanes_qua(ip)
  ENDDO
  WRITE(stdout,*)
  DEALLOCATE(aux)
  DEALLOCATE(psi)

  ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  !  Determines the index of the first projector of the absorbing atom
  ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$  

  ipx_0=0
  IF(xiabs.NE.1) THEN
     DO nt=1, xiabs-1
        DO na=1, nat
           IF (ityp(na).EQ.nt) THEN
              ipx_0=ipx_0+paw_recon(nt)%paw_nh
           ENDIF
        ENDDO
     ENDDO
  ENDIF


  ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  ! Starts the loop over the k-points
  ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$  


  !*apsi  call set_vrs(vrs,vltot,vr,nrxx,nspin,doublegrid)
  !  CALL set_vrs(vrs,vltot,v%of_r,nrxx,nspin,doublegrid)
  !  CALL newd

  !<CG>
  ! set_vrs (vrs, vltot, vr, kedtau, kedtaur,nrxx, nspin, doublegrid)
  CALL set_vrs(vrs,vltot,v%of_r,kedtau, v%kin_r,dfftp%nnr,nspin,doublegrid)
  !</CG>


  DO ik=1,nks

     WRITE(stdout,'(8x,a)')&
           '|-------------------------------------------------------------'
     WRITE(stdout,'(8x,a ,i5,a,3(f7.4,a))') '! k-point # ',ik, &
        ':  (', xk(1,ik),', ',xk(2,ik),', ',xk(3,ik),') '
     WRITE(stdout,'(8x,a ,f7.4,a,i3)') '! weight:',wk(ik), &
        ' spin state:', isk(ik)
     WRITE(stdout,'(8x,a)')&
           '|-------------------------------------------------------------'
     
     IF (calculated(1,ik).EQ.1) CYCLE

     timenow=get_clock( 'xanes' )
     CALL mp_bcast(timenow,root_pool, intra_pool_comm)

     IF( timenow.GT.time_limit) THEN
       nunfinished=1
       EXIT
     ENDIF

     IF (verbosity.eq.'high') &
     WRITE(stdout,'(8x, "| Total cpu time spent up to now is ",F9.2," s")')&
           timenow

     current_k=ik 
     IF(lsda) current_spin=isk(ik)
     CALL g2_kin(ik)

     npw = ngk(ik)

     npw_partial = npw
     CALL mp_sum( npw_partial, intra_pool_comm )

     IF(xniter.ge.npw_partial) THEN
        xniter = npw_partial

        WRITE(stdout,'(8x,a)') '|   Hilbert space is saturated'
        WRITE(stdout,'(8x,a,i10)') '|   xniter is set equal to ',npw_partial
        WRITE(stdout,'(8x,a)') &
               '|   Increase kinetic-energy cutoff in your SCF calculation!'
        !        CALL stop_pp
     ENDIF

     !<CG>
     CALL init_gipaw_2(npw,igk_k(1,ik),xk(1,ik),paw_vkb)
     !</CG>
     if(.not.lda_plus_u) CALL init_us_2(npw,igk_k(1,ik),xk(1,ik),vkb)
     IF (lda_plus_u) CALL orthoUwfc_k(ik)

     ! Angular Matrix element
     !
     !... Calculates the complex PAW projectors, paw_vkb_cplx, from
     !     LC of paw_vkb, i.e., the real PAW projectors expressed using
     !     real spherical harmonics)
     !
     !*************************************************************************
     ! Here I define human projectors
     !
     ! paw_vkb(1:npw,ipx) are real spherical harmonics. The real spherical 
     ! harmonics are defined as
     ! 
     !     y_{l,2m}  =[(-)^{m} Y_{l,m}+Y_{l,-m}]/SQRT(2)
     !     y_{l,2m+1}=[(-)^{m} Y_{l,m}-Y_{l,-m}]/(i*SQRT(2))
     !
     ! The complex spherical harmonics can be written has a function of the real
     ! ones as :
     !
     !     Y_{l,m}  =       [y_{l,2m}+iy_{l,2m+1}]/SQRT(2)     
     !     Y_{l,-m} = (-)^m [y_{l,2m}-iy_{l,2m+1}]/SQRT(2)     
     !
     !  The paw_vkb_cplx are the Y_{l,m}, so the usual spherical harmonics
     !*************************************************************************



     DO ip=1,paw_recon(xiabs)%paw_nl(xang_mom)   

        ipx=ipx_0+paw_iltonhb(xang_mom,ip,xiabs)
        paw_vkb_cplx(1:npw,ipx)= paw_vkb(1:npw,ipx)   !m=0
        paw_vkb_cplx(1:npw,ipx+1)=       &
             (paw_vkb(1:npw,ipx+1)+(0.d0,1.d0)*paw_vkb(1:npw,ipx+2))/SQRT(2.0)  
             !m=+1
        paw_vkb_cplx(1:npw,ipx+2)=       &
             -(paw_vkb(1:npw,ipx+1)-(0.d0,1.d0)*paw_vkb(1:npw,ipx+2))/SQRT(2.0)
             !m=-1
        paw_vkb_cplx(1:npw,ipx+3)=       &
             (paw_vkb(1:npw,ipx+3)+(0.d0,1.d0)*paw_vkb(1:npw,ipx+4))/SQRT(2.0)
             !m=+2
        paw_vkb_cplx(1:npw,ipx+4)=       &
             (paw_vkb(1:npw,ipx+3)-(0.d0,1.d0)*paw_vkb(1:npw,ipx+4))/SQRT(2.0)
             !m=-2
     ENDDO


     !... Calculates the initial state for the Lanczos procedure,
     !    stored in psiwfc.
     !    It includes the radial matrix element, the angular matrix element,
     !    and the associated PAW projector.
     ! 

     psiwfc(:)=(0.d0,0.d0)

     DO ip=1,paw_recon(xiabs)%paw_nl(xang_mom)   
        ipx=ipx_0+paw_iltonhb(xang_mom,ip,xiabs)

        !      WARNING, storage of spherical harmonics is the following:
        !         given Y_{l,m}=P_{lm}exp(im\phi) one has
        !                                                         counter
        !   l, m=0  -------> Y_{l,0}                              1     z
        !   l, m=+1 -------> Y_{l,1}+Y_{l,-1} or cos(m\phi)       2     
        !   l, m=-1 -------> Y_{l,1}-Y_{l,-1} or sin(m\phi)       3     
        !  .....
        !   l, m=+l -------> Y_{l,1}+Y_{l,-1} or cos(m\phi)      2*l
        !   l, m=-l -------> Y_{l,1}-Y_{l,-1} or sin(m\phi)      2*l+1

        psiwfc(1:npw)=psiwfc(1:npw)+prefm2*&
             (xepsilon(1)-(0.d0,1.d0)*xepsilon(2))*&
             (xkvec(1)-(0.d0,1.d0)*xkvec(2))*&
             xanes_qua(ip)*&
             paw_vkb_cplx(1:npw,ipx+3)

        psiwfc(1:npw)=psiwfc(1:npw)+prefm2*&
             (xepsilon(1)+(0.d0,1.d0)*xepsilon(2))*&
             (xkvec(1)+(0.d0,1.d0)*xkvec(2))*&
             xanes_qua(ip)*&
             paw_vkb_cplx(1:npw,ipx+4)

        psiwfc(1:npw)=psiwfc(1:npw)-prefm1*( &
             (xepsilon(1)-(0.d0,1.d0)*xepsilon(2))* &
             xkvec(3)+ &
             (xkvec(1)-(0.d0,1.d0)*xkvec(2))* &
             xepsilon(3) )* &
             xanes_qua(ip)*&
             paw_vkb_cplx(1:npw,ipx+1)

        psiwfc(1:npw)=psiwfc(1:npw)+prefm1*( &
             (xepsilon(1)+(0.d0,1.d0)*xepsilon(2))* &
             xkvec(3)+ &
             (xkvec(1)+(0.d0,1.d0)*xkvec(2))* &
             xepsilon(3) )* &
             xanes_qua(ip)*&
             paw_vkb_cplx(1:npw,ipx+2)

        psiwfc(1:npw)=psiwfc(1:npw)+prefm0*( &
             pref*xkvec(3)*xepsilon(3)-&
             prefb*(xepsilon(1)*xkvec(1)+xepsilon(2)*xkvec(2)))*&
             xanes_qua(ip)*&
             paw_vkb_cplx(1:npw,ipx)

     ENDDO


      !... Normalizes the wavefunction psiwfc(1:npw)

!     CALL allocate_bec_type(nkb,1, becp) ! CG
     CALL allocate_bec_type(nkb,natomwfc, becp) 
     !<CG>
     IF (okvan) THEN
        ALLOCATE(spsiwfc(npwx))
        spsiwfc(:)=(0.d0,0.d0)
        recalc=.true.
        CALL sm1_psi(recalc,npwx, npw, 1, psiwfc, spsiwfc)
        xnorm_partial=zdotc(npw,psiwfc,1,spsiwfc,1)
        DEALLOCATE(spsiwfc)
     ELSE
        xnorm_partial=zdotc(npw,psiwfc,1,psiwfc,1)
     ENDIF
     !</CG>

     CALL mp_sum( xnorm_partial, intra_pool_comm )

     xnorm(1,ik)=SQRT(xnorm_partial)
     WRITE(stdout,'(8x,a,e15.8)') '|   Norm of the initial Lanczos vector:',&
                                     xnorm(1,ik)
     norm=1.d0/xnorm(1,ik)

     CALL zdscal(npw,norm,psiwfc,1)
    
     !... Starts the Lanczos procedure

     IF (okvan) THEN
        CALL lanczos_uspp(a(:,1,ik),b(:,1,ik),npw,psiwfc,ncalcv(1,ik), terminator)
     ELSE
        CALL lanczos(a(:,1,ik),b(:,1,ik),npw,psiwfc,ncalcv(1,ik), terminator)
     ENDIF
     
     !!      Then I write small report of the lanczos results
     !IF(TRIM(verbosity).EQ.'high') THEN
      ! WRITE( stdout,*) '-----------------------------------------'
      ! WRITE( stdout,*) 'k-point number =',ik
      ! WRITE( stdout,*) 'Norm of the initial vector =',xnorm(1,ik)
      ! WRITE( stdout,*) 'Number of iterations =',ncalcv(1,ik)
      ! WRITE( stdout,*) 'k-point coordinate, isk'
      ! WRITE( stdout,'(3f12.6,1i2)') xk(1,ik),xk(2,ik),xk(3,ik),isk(ik)
        !     nline=ncalcv(1,ik)/6
        !     nrest=ncalcv(1,ik)-nline*6
        !     WRITE( stdout,*) 'a vectors:'
        !     DO ip=1,nline
        !        WRITE( stdout,"(6(f10.6,3x))") (a((ip-1)*6+j,1,ik),j=1,6)
        !     ENDDO
        !     WRITE( stdout,"(6(f10.6,3x))") (a(nline*6+j,1,ik),j=1,nrest)
        !     WRITE( stdout,*) 'b vectors:'
        !     DO ip=1,nline
        !        WRITE( stdout,"(6(f10.6,3x))") (b((ip-1)*6+j,1,ik),j=1,6)
        !     ENDDO
        !     WRITE( stdout,"(6(f10.6,3x))") (b(nline*6+j,1,ik),j=1,nrest)
      ! WRITE( stdout,*) '-----------------------------------------'
    !ENDIF

     CALL deallocate_bec_type (becp) ! CG
     calculated(1,ik)=1

  ENDDO   !LOOP on k-points

  CALL mp_barrier(world_comm)
  CALL mp_sum(nunfinished, world_comm)
  IF (nunfinished >= 1) THEN
    save_file_kind='unfinished'
    write(stdout,'(5x,a)') 'Calculation not finished'
  ENDIF

  DEALLOCATE(psiwfc)
  DEALLOCATE (paw_vkb_cplx)
  DEALLOCATE(xanes_qua)

END SUBROUTINE xanes_quadrupole

