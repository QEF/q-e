
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!  XANES K-edge calculation in the electric dipole approximation
!------------------------------------------------------------------------------
SUBROUTINE xanes_dipole(a,b,ncalcv,xnorm,core_wfn,paw_iltonhb,&
                        terminator,verbosity)
  !----------------------------------------------------------------------------
  USE constants,       ONLY : fpi
  USE io_global,       ONLY : stdout     ! Modules/io_global.f90
  USE kinds,           ONLY : DP
  USE parameters,      ONLY : ntypx
  USE radial_grids,    ONLY : ndmx
  USE ions_base,       ONLY : nat, ntyp => nsp, ityp
  USE wvfct,           ONLY : npwx, nbnd, et, current_k
  USE gvecw,           ONLY : gcutw
  USE symm_base,       ONLY : d1,d2,d3
  USE noncollin_module,ONLY : noncolin
  USE lsda_mod,        ONLY : nspin,lsda,isk,current_spin
  USE cell_base,       ONLY : tpiba2, bg
  USE wavefunctions_module, ONLY: evc
  USE klist,           ONLY : &
       nkstot,                & ! total number of k-points
       nks,                   & ! number of k-points per pool
       xk,                    & ! k-points coordinates
       wk,                    & ! k-points weight
       ngk, igk_k
  USE gvect,           ONLY: g, ngm, ngl
  USE fft_base,        ONLY: dfftp
  USE paw_gipaw,       ONLY : &
       paw_vkb,               & ! |p> projectors
       paw_becp,              & ! product of projectors and wf.
       paw_nkb,               & ! total number of beta functions, with st.fact.
       paw_lmaxkb,paw_recon
  USE becmod,          ONLY : becp, allocate_bec_type, deallocate_bec_type !CG
  USE scf,             ONLY : vltot, vrs, v, kedtau
  USE gvecs,           ONLY : doublegrid
  USE mp_world,        ONLY : world_comm
  USE mp_pools,        ONLY : intra_pool_comm, root_pool, npool
  USE mp,              ONLY : mp_sum, mp_bcast, mp_barrier !CG
  USE io_global,       ONLY : ionode

  USE xspectra,        ONLY : xiabs, xanes_dip, xang_mom, xniter,&
                              xnitermax, xepsilon,time_limit,calculated,&
                              save_file_kind
  USE atom,            ONLY : rgrid, msh
  !  use atom,        ONLY : &
  !       mesh,     &!mesh(ntypx) number of mesh points
  !       msh ,     &!msh(ntypx)the point at rcut=end of radial integration
  !       r   
  USE radin_mod
  USE basis,           ONLY : natomwfc
  USE uspp,            ONLY : vkb, nkb, okvan !CG
  USE uspp_param,      ONLY : upf
  USE ldaU,            ONLY : lda_plus_u, init_lda_plus_u 
  !<CG>
  USE xspectra_paw_variables, ONLY : xspectra_paw_nhm
  !</CG>

  IMPLICIT NONE
  !
  REAL(dp), INTENT(INOUT) :: a(xnitermax,1,nks)
  REAL(dp), INTENT(INOUT) :: b(xnitermax,1,nks)     
  REAL(dp), INTENT(INOUT) :: xnorm(1,nks)
  REAL(dp), INTENT(IN)    :: core_wfn(ndmx)
  INTEGER, INTENT(INOUT)  :: ncalcv(1,nks)
  INTEGER, INTENT(IN)     :: paw_iltonhb(0:paw_lmaxkb,xspectra_paw_nhm,ntyp)
  LOGICAL, INTENT(IN)     :: terminator
  !
  !... Local variables
  !
  INTEGER  :: is,ik,iabso,nr,ip,jp,l,j,icrd,ip_l,nrc,nt,na
  INTEGER  :: ipx,ipx_0,ipy,ipz,nline,nrest,npw, npw_partial
  INTEGER  :: nunfinished
  LOGICAL  :: recalc
  REAL(dp) :: pref,prefb,v_of_0,xnorm_partial
  REAL(dp) :: norm, normps
  REAL(dp), ALLOCATABLE :: aux(:)
  COMPLEX(KIND=DP), EXTERNAL :: zdotc
  COMPLEX(dp), ALLOCATABLE :: paw_vkb_cplx(:,:)
  COMPLEX(dp), ALLOCATABLE :: psiwfc(:), spsiwfc(:)
  CHARACTER(LEN=4) :: verbosity

  REAL(dp) :: timenow 
  REAL(DP), EXTERNAL ::  get_clock
  EXTERNAL :: zdscal

  timenow=0
  nunfinished=0

  ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  !  Constant definitions
  ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  pref=SQRT(3.d0/2.d0)
  prefb=SQRT(3.d0)

  ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  !  Variable allocation and initialization
  ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  ALLOCATE(aux(rgrid(xiabs)%mesh)) !allocation too big, it needs only up to msh
  ALLOCATE (paw_vkb_cplx( npwx,  paw_nkb))
  ALLOCATE(xanes_dip(paw_recon(xiabs)%paw_nl(xang_mom)))
  ALLOCATE(psiwfc(npwx))
  xanes_dip(:)=0.d0

  ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  ! Radial Dipole Matrix element
  ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  !  radial part:  <\phi_n|r|\psi_i>=int r^3 \phi_n \psi_i dr.
  !  Here only ps and ae wavefunctions from XX.recon
  !
  WRITE(stdout,'(8x,a)') &
       'Radial transition matrix element(s) used in the calculation of the'
  WRITE(stdout,'(8x,a)') &
       'initial vector of the Lanczos basis (|tilde{phi}_abs> normalized)'
  !WRITE(stdout,'(5x,a,i2,a,i2,a,i3,a)') 'There are ',&
  !                                      paw_recon(xiabs)%paw_nl(xang_mom),&
  !                                      ' projector(s)/channel for l=',&
  !                                      xang_mom,' and atom type',&
  !                                      xiabs,'.'

  ! ... Checks that the core wf is correctly normalized
  !
  nr=msh(xiabs)  ! extended up to all the NON ZERO points in the mesh.
  IF(TRIM(verbosity).EQ.'high') THEN
     aux(1:nr)=core_wfn(1:nr)*core_wfn(1:nr)
     WRITE (stdout,'(8x,"Norm of core wfc = ",f10.6)') &
           SQRT(para_radin(aux(1:nr),rgrid(xiabs)%r(1:nr),nr))
  ENDIF

  ! ... Calculates the radial integral
  !
  ip_l=0

  DO ip=1,paw_recon(xiabs)%paw_nbeta
     !  IF(psphi(xiabs,ip)%label%l.EQ.xang_mom) THEN
     IF(paw_recon(xiabs)%aephi(ip)%label%l.EQ.xang_mom) THEN
        nrc=paw_recon(xiabs)%aephi(ip)%label%nrc
        ip_l=ip_l+1
        ! here below, recall that psi is r*psi and you have a Jacobian=r^2
        aux(1:nrc) = rgrid(xiabs)%r(1:nrc) * &
                     paw_recon(xiabs)%aephi(ip)%psi(1:nrc) * &
                     core_wfn(1:nrc)
        ! here we have to integrate only inside the augmentation region.
        xanes_dip(ip_l)=para_radin(aux(1:nrc),rgrid(xiabs)%r(1:nrc),nrc)
     ENDIF
  ENDDO

  !... Writes the radial transition matrix element(s)
  ! 
  DO ip=1,paw_recon(xiabs)%paw_nl(xang_mom) !
     WRITE(stdout,'(8x,a,i1,a,i1,a,f14.9)') &
            '| For PAW proj. (l=',xang_mom,') #',&
            ip, ': radial matrix element =', xanes_dip(ip)
  ENDDO
  WRITE(stdout,*)
  DEALLOCATE(aux)

  ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  !  Determines the index of the first projector of the absorbing atom
  ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  !
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
  !<CG>
  CALL set_vrs(vrs,vltot,v%of_r,kedtau, v%kin_r,dfftp%nnr,nspin,doublegrid)
  !</CG>

  IF (npool /= 1) WRITE(stdout,'(a,i5,a,i3)') 'NB: the ', nks,&
     ' k-point are not all listed below because npool=', npool

  DO ik=1,nks

     WRITE(stdout,'(8x,a)')&
           '|-------------------------------------------------------------'
     WRITE(stdout,'(8x,a ,i5,a,3(f7.4,a),f7.4,a,i3)') '! k-point # ',ik, &
        ':  (', xk(1,ik),', ',xk(2,ik),', ',xk(3,ik),'), ',wk(ik),', ', isk(ik)
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
       WRITE(stdout,'(8x,"|   Total cpu time spent up to now: ",F9.2," s")')&
             timenow 

     current_k=ik 
     IF(lsda) current_spin=isk(ik)
     CALL g2_kin (ik)

     npw = ngk(ik)

     npw_partial = npw
     CALL mp_sum( npw_partial, intra_pool_comm )

     IF(xniter.ge.npw_partial) THEN
        xniter = npw_partial
        WRITE(stdout,'(8x,a)') '|   Hilbert space is saturated'
        WRITE(stdout,'(8x,a,i10)') '|   xniter is set equal to ',npw_partial
        WRITE(stdout,'(8x,a)') &
               '|   Increase kinetic-energy cutoff in your SCF calculation!'
     ENDIF

     !<CG>        
     CALL init_gipaw_2(npw,igk_k(1,ik),xk(1,ik),paw_vkb)
     !</CG>
     IF (.NOT.lda_plus_u) CALL init_us_2(npw,igk_k(1,ik),xk(1,ik),vkb)
     IF (lda_plus_u) CALL orthoUwfc_k(ik)

     ! Angular Matrix element
     !
     !... Calculates the complex PAW projectors, paw_vkb_cplx, from
     !     LC of paw_vkb, i.e., the real PAW projectors expressed using
     !     real spherical harmonics)
     !
     !*************************************************************************
     ! Here I define human projectors <CG>
     !
     ! paw_vkb(1:npw,ipx) are real spherical harmonics. 
     ! The real spherical harmonics are defined as
     ! 
     !     y_{l,2m}  =[Y_{l,m}+(-1)^m Y_{l,-m}]/SQRT(2)
     !     y_{l,2m+1}=[Y_{l,m}-(-1)^m Y_{l,-m}]/(i*SQRT(2))
     !
     ! (remember Y_{l,m}=(-1)^m Y_{l,-m)^*  )
     !
     ! The complex spherical harmonics can be written has a function of the real
     ! ones as:
     !
     !     Y_{l,m}  =        [y_{l,2m}+iy_{l,2m+1}]/SQRT(2)     
     !     Y_{l,-m} = (-1)^m [y_{l,2m}-iy_{l,2m+1}]/SQRT(2)     
     !
     ! The paw_vkb_cplx are the Y_{l,m} so the usual spherical harmonics.
     !
     ! Rotational invariance has been checked.
     !*************************************************************************

     

     DO ip=1,paw_recon(xiabs)%paw_nl(xang_mom)

        ipx=ipx_0+paw_iltonhb(xang_mom,ip,xiabs)
        ! m= 0
        paw_vkb_cplx(1:npw,ipx)= paw_vkb(1:npw,ipx)
        ! m=+1 
        paw_vkb_cplx(1:npw,ipx+1)=       &
             (paw_vkb(1:npw,ipx+1)+(0.d0,1.d0)*paw_vkb(1:npw,ipx+2))/SQRT(2.0)
        ! m=-1 
        paw_vkb_cplx(1:npw,ipx+2)=-       &
             (paw_vkb(1:npw,ipx+1)-(0.d0,1.d0)*paw_vkb(1:npw,ipx+2))/SQRT(2.0)

     ENDDO



     !... Calculates the initial state for the Lanczos procedure,
     !    stored in psiwfc.
     !    It includes the radial matrix element, the angular matrix element,
     !    and the associated PAW projector.
     !
     psiwfc(1:npw)=(0.d0,0.d0)
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


        psiwfc(1:npw)=psiwfc(1:npw) + ( pref                                   &
                                        * paw_vkb_cplx(1:npw,ipx+2)            &
                                        * (xepsilon(1)+(0.d0,1.d0)*xepsilon(2))&
                                       +pref                                   &
                                        * paw_vkb_cplx(1:npw,ipx+1)            &
                                        *(-xepsilon(1)+(0.d0,1.d0)*xepsilon(2))&
                                       +prefb                                  &
                                        *paw_vkb_cplx(1:npw,ipx)               &
                                        *xepsilon(3)                           &
                                      )                                        &
                                    * xanes_dip(ip)/SQRT(fpi)

     ENDDO
     psiwfc(1:npw)=psiwfc(1:npw)*SQRT(fpi)/3.0

     !... Normalizes the wavefunction psiwfc(1:npw)
     !

     !<CG>
     CALL allocate_bec_type(nkb,1,becp)

     write(6,*) 'okvan=',okvan
     IF (okvan) THEN
        ALLOCATE(spsiwfc(npwx))
        spsiwfc(:)=(0.d0,0.d0)
        recalc=.true.
        CALL sm1_psi(recalc,npwx, npw, 1, psiwfc, spsiwfc)
        xnorm_partial=zdotc(npw,psiwfc,1,spsiwfc,1)
        DEALLOCATE(spsiwfc)
     ELSE
!        xnorm_partial=0.d0
!        do ip=1,npw
!          xnorm_partial=xnorm_partial+conjg(psiwfc(ip))*psiwfc(ip)
!       enddo
        xnorm_partial=real(zdotc(npw,psiwfc,1,psiwfc,1),dp)

     ENDIF
     !</CG>

     CALL mp_sum( xnorm_partial, intra_pool_comm )

     xnorm(1,ik)=SQRT(xnorm_partial)
     WRITE(stdout,'(8x,a,e15.8)') '|   Norm of the initial Lanczos vector:',&
          xnorm(1,ik)

     norm=1.d0/xnorm(1,ik)
     
     CALL zdscal(npw,norm,psiwfc,1)
    
     !... Starts the Lanczos procedure
     !
     IF (okvan) THEN
        CALL lanczos_uspp(a(:,1,ik),b(:,1,ik),npw,psiwfc,ncalcv(1,ik),terminator)
     ELSE
        CALL lanczos(a(:,1,ik),b(:,1,ik),npw,psiwfc,ncalcv(1,ik),terminator)
     ENDIF

     !!      Then I write small report of the lanczos results
     !!
     !IF(TRIM(verbosity).EQ.'high') THEN
     !   WRITE( stdout,*) '-----------------------------------------'
     !   WRITE( stdout,*) 'k-point number =',ik
     !   WRITE( stdout,*) 'k-point coordinate, isk'
     !   WRITE( stdout,'(3f12.6,1i2)') xk(1,ik),xk(2,ik),xk(3,ik),isk(ik)
     !   WRITE( stdout,*) 'Norm of the initial vector =',xnorm(1,ik)
     !   WRITE( stdout,*) 'Number of iterations =',ncalcv(1,ik)
     !   !        nline=ncalcv(icrd,ik)/6
     !   !        nrest=ncalcv(icrd,ik)-nline*6
     !   !        WRITE( stdout,*) 'a vectors:'
     !   !        DO ip=1,nline
     !   !           WRITE( stdout,"(6(f10.6,3x))") (a((ip-1)*6+j,icrd,ik),j=1,6)
     !   !        ENDDO
     !   !        WRITE( stdout,"(6(f10.6,3x))") (a(nline*6+j,icrd,ik),j=1,nrest)
     !   !        WRITE( stdout,*) 'b vectors:'
     !   !        DO ip=1,nline
     !   !           WRITE( stdout,"(6(f10.6,3x))") (b((ip-1)*6+j,icrd,ik),j=1,6)
     !   !        ENDDO
     !   !        WRITE( stdout,"(6(f10.6,3x))") (b(nline*6+j,icrd,ik),j=1,nrest)
     !   WRITE( stdout,*) '-----------------------------------------'
     !ENDIF

     CALL deallocate_bec_type ( becp ) ! CG
     calculated(1,ik)=1  

  ENDDO  ! LOOP on k-points
 
  CALL mp_barrier(world_comm)
  CALL mp_sum(nunfinished, world_comm)
  IF (nunfinished >= 1) THEN 
    save_file_kind='unfinished'
    WRITE(stdout,'(5x,a)') 'Calculation not finished'
  ENDIF

  DEALLOCATE(psiwfc)
  DEALLOCATE(xanes_dip)
  DEALLOCATE (paw_vkb_cplx)

END SUBROUTINE xanes_dipole
