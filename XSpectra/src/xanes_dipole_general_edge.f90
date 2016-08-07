!check xang_mom in this case, for dipole it should be one..

! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!  General Calculation ! OB
! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$  


SUBROUTINE xanes_dipole_general_edge(a,b,ncalcv,nl_init, xnorm,core_wfn,paw_iltonhb,terminator, verbosity)
  USE kinds,           ONLY : DP
  USE constants,       ONLY : pi, fpi
  USE io_global,       ONLY : stdout     ! Modules/io_global.f90
  USE parameters,      ONLY : ntypx
  USE radial_grids,    ONLY : ndmx
  USE ions_base,       ONLY : nat, ntyp => nsp, ityp
  USE wvfct,           ONLY : npwx, nbnd, et, current_k
  USE gvecw,           ONLY : gcutw
  USE lsda_mod,        ONLY : nspin,lsda,isk,current_spin
  USE cell_base,       ONLY: tpiba2, bg
  USE wavefunctions_module, ONLY: evc
  USE klist,           ONLY : &
       nkstot,            & ! total number of k-points
       nks,               & ! number of k-points per pool
       xk,                & ! k-points coordinates
       wk,                & ! k-points weight
       ngk, igk_k
  USE gvect,           ONLY: g, ngm, ngl
  USE fft_base,        ONLY: dfftp
  USE paw_gipaw,       ONLY : &
       paw_vkb,             & ! |p> projectors
       paw_becp,            & ! product of projectors and wf.
       paw_nkb,             & ! total number of beta functions, with st.fact.
       paw_lmaxkb,paw_recon
  USE becmod,          ONLY : becp, allocate_bec_type, deallocate_bec_type !CG
  USE scf,             ONLY : vltot, vrs, v, kedtau
  USE gvecs,           ONLY : doublegrid
  USE mp_global,       ONLY : intra_pool_comm, root_pool, world_comm
  USE mp,              ONLY : mp_sum, mp_bcast, mp_barrier !CG
  USE mp_pools,        ONLY : npool
  USE io_global,       ONLY : ionode

  USE xspectra,        ONLY : edge, n_lanczos, xiabs, xang_mom, xniter, &
                              xnitermax, xepsilon,time_limit,calculated,&
                              save_file_kind, lplus, lminus,  two_edges
  USE coef_gaunt
  USE coef_CG
  USE atom,            ONLY : rgrid, msh
  USE radin_mod
  USE uspp,   ONLY : vkb, nkb, okvan !CG
  USE ldaU,   ONLY : lda_plus_u
  !<CG>
  USE xspectra_paw_variables, ONLY : xspectra_paw_nhm
  !</CG>

  IMPLICIT NONE
  REAL(dp) core_wfn(ndmx)
  REAL(dp) a(xnitermax,n_lanczos,nks),b(xnitermax,n_lanczos,nks)     
  REAL (dp)  xnorm(n_lanczos,nks)
  INTEGER :: is,ik,iabso,nr,ip,jp,l,icrd,ip_l,nrc,nt,na,lf, &
             lm, llm, m, lmi, lmf, isg, isgf, mu, ll, mf, mpl, mmi, ms,&
             lmbd, lmmax, lfmax, no, noj,jloop
  INTEGER :: ipx,ipx_0,ipy,ipz,nline,nrest,npw, npw_partial
  integer :: l_final_dim, i_lanczos
  integer, dimension(2), intent(in):: nl_init
  integer, allocatable:: l_final(:)
  logical, allocatable:: q_final(:)
  INTEGER :: ncalcv(n_lanczos,nks)
  logical :: nocalc(n_lanczos,nks)
  INTEGER :: paw_iltonhb(0:paw_lmaxkb,xspectra_paw_nhm, ntyp)
  REAL (dp) facxy, facz, facE2, facE1, v_of_0,xnorm_partial
  REAL (dp) norm, mss, s, j, mj, dl, disg, CG
  REAL (dp), ALLOCATABLE :: aux(:)
  REAL (dp), ALLOCATABLE :: Mxanes(:,:)
  COMPLEX(KIND=DP), EXTERNAL :: zdotc
  COMPLEX(dp), ALLOCATABLE :: paw_vkb_cplx(:,:)
  LOGICAL :: terminator
  REAL(dp) :: normps
  COMPLEX(dp),dimension(3):: y_eps
  COMPLEX(dp), ALLOCATABLE :: y_k(:)
  EXTERNAL zdscal

  LOGICAL :: recalc
  COMPLEX(dp), ALLOCATABLE :: psiwfc(:), spsiwfc(:), term(:)
  REAL(DP), EXTERNAL ::  get_clock
  REAL(dp) :: timenow=0 
  INTEGER :: nunfinished=0
  CHARACTER (LEN=4)   :: verbosity

  ! <OB>
  ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  ! Identify the l of the final states, according to the selection rules of the specific calculation 
  !         plus  = .true.  -> \delta l =  1
  !         lminus = .true.  -> \delta l = -1
  !         lplus=.false .and. lminus=.false                       -> \delta l = \pm 1
  ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  !
  ! MCB: the quadrupolar part has never been tested, so only the dipolar is used.
  !



  l_final_dim = 5
  allocate( l_final(l_final_dim) )                         ! max available l for the final states
  allocate( q_final(l_final_dim) )                         ! flag (quadrupole trs or not of the final state)

  ! Initialization:
  l_final(:) = -1
  q_final(:) = .false.
  
  l_final(1) = nl_init(2) + xang_mom
  if( .not. lplus .and. nl_init(2) > xang_mom - 1 ) &    ! to make sure these states exist
       l_final(2) = nl_init(2) - xang_mom

  lfmax = nl_init(2) + xang_mom  
  
  
  lmmax = lfmax**2 + 2*lfmax + 1

  ! <OB>
  ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  !  Variable allocation and initialization
  ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$  


  ALLOCATE(aux(rgrid(xiabs)%mesh)) !allocation too big, it needs only up to msh
  ALLOCATE (paw_vkb_cplx( npwx,  paw_nkb))
  ALLOCATE( Mxanes(paw_recon(xiabs)%paw_nl(lfmax),lmmax) )
  ALLOCATE(psiwfc(npwx),term(npwx))
  Mxanes(:,:)=0.d0


  !<OB>
  if( lplus )   write( stdout,'("Transitions calculated according to the \delta_l = +",i2," rule")') xang_mom
  if( lminus )  write( stdout,'("Transitions calculated according to the \delta_l = -",i2," rule")') xang_mom
  write(stdout,*)

  !<OB>

  ! This has been adapted for the selection rules explained above 
 
  do ll = 1, l_final_dim
     lf = l_final(ll)
     if( lf == -1 ) exit      ! state does not exist
  end do

  
  ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  ! Radial Dipole Matrix element
  ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$  


  !  I compute the radial part, 
  !          <\phi_n|r|\psi_i>=int r^3 \phi_n \psi_i dr.
  !  Here only ps and ae wavefunctions from XX.recon


  !
  WRITE(stdout,'(8x,a)') &
       'Radial transition matrix element(s) used in the calculation of the'
  WRITE(stdout,'(8x,a)') &
       'initial vector of the Lanczos basis (|tilde{phi}_abs> normalized)'

  !
  ! I check that the core wf is correctly normalized
  
  nr=msh(xiabs)  ! extended up to all the NON ZERO points in the mesh.
  IF(TRIM(verbosity).EQ.'high') THEN
     aux(1:nr)=core_wfn(1:nr)*core_wfn(1:nr)
     WRITE (stdout,'(8x,"Norm of core wfc = ",f10.6)') &
          SQRT(para_radin(aux(1:nr),rgrid(xiabs)%r(1:nr),nr))
  ENDIF


  
  ! <OB>
  ! Calculation of the radial integral
  ! this was modified to contain the several radial integrals required by selection rules 
  !                 and the  m quantum number under the unified index lm
  ! Convention: lm = l**2 + 2m       m>0
  !             lm = l**2 + 2m + 1   m<0
  !             lm = l**2            m=0 
  ! In Mxanes, only the l (and m) compatible with the selection rules are filled
  ! At this stage, the m quantitative dependence of the radial wfc is neglected
  ! <OB>
  
  do ll = 1, l_final_dim
     lf = l_final(ll)
     if( lf == -1 ) exit
     ip_l = 0                                ! no of projectors for the current l 
     DO ip=1,paw_recon(xiabs)%paw_nbeta      ! paw_nbeta disregards m, as I need here
        IF(paw_recon(xiabs)%aephi(ip)%label%l == lf) THEN
           ip_l = ip_l + 1
           nrc=paw_recon(xiabs)%aephi(ip)%label%nrc
           ! <OB>
           if( q_final(ll) ) then
              aux(1:nrc)=rgrid(xiabs)%r(1:nrc)*paw_recon(xiabs)%aephi(ip)%psi(1:nrc) * core_wfn(1:nrc)
              aux(1:nrc)=aux(1:nrc)*rgrid(xiabs)%r(1:nrc)
           else
              ! <OB>
              !    psi is r*psi and the Jacobian=r^2
              aux(1:nrc)=rgrid(xiabs)%r(1:nrc)*paw_recon(xiabs)%aephi(ip)%psi(1:nrc) * core_wfn(1:nrc)
           end if
           !    one has to integrate only inside the augmentation region.
           lmi = lf**2 + 1 
           lmf = lf**2 + 2*lf + 1  
           Mxanes(ip_l,lmi:lmf)=para_radin( aux(1:nrc),rgrid(xiabs)%r(1:nrc),nrc )
        ENDIF
     ENDDO
  end do

  
  
  do ll = 1, l_final_dim
     lf = l_final(ll)
     if( lf == -1 ) exit
     DO ip=1,paw_recon(xiabs)%paw_nl(lf) !
        WRITE( stdout,'("l = ",i2,"    Radial matrix element proj. (",i2,")=",f14.8)') & 
             lf, ip, Mxanes(ip,lf**2+1)
     ENDDO
  end do

  DEALLOCATE(aux)
  
  
  !
  !  We count the projectors for all the atoms
  !
  !  
  
  if( xiabs == 1 ) then  ! this is default
     ipx_0=0    
  else
     DO nt=1, xiabs-1
        DO na=1, nat
           IF (ityp(na).EQ.nt) ipx_0=ipx_0+paw_recon(nt)%paw_nh 
        ENDDO
     ENDDO
  end if


  

  IF (npool /= 1) WRITE(stdout,'(a,i5,a,i3)') 'NB: the ', nks,&
     ' k-point are not all listed below because npool=', npool

  
  ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  ! Beginning the loop over the Lanczos procedures   OB
  ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$  
  
  do i_lanczos = 1, n_lanczos
      
    ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    ! Beginning the loop over the k-points
    ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$  

     write(stdout,*)
     write(stdout,*) '              Lanczos Number      ', i_lanczos, ' of ', n_lanczos
     write(stdout,*)

     CALL set_vrs(vrs,vltot,v%of_r,kedtau, v%kin_r,dfftp%nnr,nspin,doublegrid)     
     
     ! Loop on the k points
     DO ik=1,nks
        
        WRITE(stdout,'(8x,a)')&
             '|-------------------------------------------------------------'
        WRITE(stdout,'(8x,a ,i5,a,3(f7.4,a),f7.4,a,i3)') '! k-point # ',ik, &
             ':  (', xk(1,ik),', ',xk(2,ik),', ',xk(3,ik),'), ',wk(ik),', ', isk(ik)
        WRITE(stdout,'(8x,a)')&
             '|-------------------------------------------------------------'
        
        
        IF (calculated(i_lanczos,ik).EQ.1) CYCLE
        
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
        ENDIF
        
        
        !<CG>        
        CALL init_gipaw_2(npw,igk_k(1,ik),xk(1,ik),paw_vkb)
        !</CG>
        IF (.NOT.lda_plus_u) CALL init_us_2(npw,igk_k(1,ik),xk(1,ik),vkb)
        IF (lda_plus_u) CALL orthoUwfc_k(ik)
        
        
        ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
        ! Angular Matrix element
        ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$  
        !
        ! Here I define human projectors   <OB> I changed notations
        !
        ! paw_vkb(1:npw,ipx) are real spherical harmonics. The real spherical harmonics are
        ! defined as
        ! 
        !     y_{l,m}  = [Y_{l,m}+(-1)^m Y_{l,-m}]/SQRT(2)     if m > 0
        !     y_{l,m}  = [Y_{l,m}-(-1)^m Y_{l,-m}]/(i*SQRT(2)) if m < 0
        !
        ! (remember Y_{l,m}=(-1)^m Y_{l,-m)^*  )
        !
        ! The complex spherical harmonics can be written has a function of the real
        ! ones as (m > 0):
        !
        !     Y_{l,m}  =        [y_{l,m}+iy_{l,-m}]/SQRT(2)     
        !     Y_{l,-m} = (-1)^m [y_{l,m}-iy_{l,-m}]/SQRT(2)     
        !
        !  The paw_vkb_cplx are the Y_{l,m} so the usual spherical harmonics
        !
        ! rotational invariance has been checked
        
        
        !      WARNING, storage of spherical harmonics is the following:
        !         given Y_{l,m}=P_{lm}exp(im\phi) one has
        !                                                         counter
        !   l, m=0  -------> Y_{l,0}                              1     z
        !   l, m=+1 -------> Y_{l,1}+Y_{l,-1} or cos(m\phi)       2     
        !   l, m=-1 -------> Y_{l,1}-Y_{l,-1} or sin(m\phi)       3     
        !  .....
        !   l, m=+l -------> Y_{l,1}+Y_{l,-1} or cos(m\phi)      2*l
        !   l, m=-l -------> Y_{l,1}-Y_{l,-1} or sin(m\phi)      2*l+1
        
        
        
        !<OB>
        ! If no spin-orbit on the final states, the radial part of the proj. is independent of m
        ! only the significant elements (according to the selection rules) are filled in paw_vkb_cplx 
        !<OB>
        
        paw_vkb_cplx(:,:) = (0.d0,0.d0)
        do ll = 1, l_final_dim
           lf = l_final(ll)
           if( lf == -1 ) exit
           DO ip=1,paw_recon(xiabs)%paw_nl(lf)
              
              ipx=ipx_0+paw_iltonhb(lf,ip,xiabs)  ! corresponds to lm 
              do mf = 0, lf
                 if( mf /= 0 ) then
                    mpl = 2*abs(mf) - 1            ! m positive, storage convention
                    mmi = 2*abs(mf)                ! m negative
                    paw_vkb_cplx(1:npw,ipx+mpl)=       &
                         (paw_vkb(1:npw,ipx+mpl)+(0.d0,1.d0)*paw_vkb(1:npw,ipx+mmi))/SQRT(2.0)
                    paw_vkb_cplx(1:npw,ipx+mmi)=       &
                         (paw_vkb(1:npw,ipx+mpl)-(0.d0,1.d0)*paw_vkb(1:npw,ipx+mmi))/SQRT(2.0) 
                    if( mod(mf, 2) /= 0 ) & 
                         paw_vkb_cplx(1:npw,ipx+mmi) = -paw_vkb_cplx(1:npw,ipx+mmi)      
                 else                          ! m = 0
                    paw_vkb_cplx(1:npw,ipx)= paw_vkb(1:npw,ipx)
                 end if
              end do
           ENDDO
        end do
 


        !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
        !          Constants and prefactors
        !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
        
        s = 0.5d0                                                                            ! electron spin
        facxy = 0.5d0*dsqrt(6/fpi)                                                           
        facz  = dsqrt(2.0d0)*facxy                                                               
        facE1 = fpi/3.                                                                       ! prefactor E1E1 (from multipolar development)
        facE2 = (fpi/3.) * (fpi/3.) * 0.5d0                                                  !           E2E2 
        
        
        
        psiwfc(1:npw)=(0.d0,0.d0)                                                            ! recalculated each k and i_lanczos
        
        !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
        !                Define j and m_j for incoming i_lanczos
        !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
        
        if( .not. two_edges ) then                                                           
           j = (n_lanczos - 1)/2.                                                             
        else                                                                                 
           if( i_lanczos < n_lanczos/2 ) then                                                 
              select case( edge(2:3) )                                                         
              case('23')                                                                     
                 j = 0.5                                                                      
              case('45')                                                                     
                 j = 1.5                                                                      
              case('67')                                                                     
                 j = 2.5                                                                      
              case default                                                                   
                 write(stdout,*) 'Needs to be extended'                                       
              end select
           else                                                                               
              select case( edge(2:3) )                                                         
              case('23')                                                                     
                 j = 1.5                                                                      
              case('45')                                                                     
                 j = 2.5                                                                      
              case('67')                                                                     
                 j = 3.5                                                                      
              case default                                                                   
                 write(stdout,*) 'Needs to be extended'                                       
              end select
           end if
        end if
                                                                                            
        if( .not. two_edges ) then
           mj =  -j - 1 + i_lanczos
        else
           if( i_lanczos < n_lanczos/2 ) then
              mj =  -j - 1 + i_lanczos
           else
              ! works for *23, *45, *67 only
              mj =  -3*j + i_lanczos                                        
           end if
        end if
        
        !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$  
        ! Define the spherical harmonics of the c.c. polarization vector    Y_1^mu* = Y_1^-mu
        !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$                                                                                        

        
        y_eps(3) = facxy * ( xepsilon(1) - (0.d0,1.d0) * xepsilon(2) )                       ! Y_1^1*  
        y_eps(2) = facz *  xepsilon(3)                                                       ! Y_1^0*
        y_eps(1) = (-1) * facxy * ( xepsilon(1) + (0.d0,1.d0) * xepsilon(2) )                ! Y_1^-1*
        
        
        !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ 
        ! Calculation of \tilde\Phi_0 :     1 Lanczos / jm_j state / polarization 
        !$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
        
        
        if( nspin == 2 ) then                                                            
           isgf = 2                                                                           ! spin of the final state                                                                
        else                                                                             
           isgf = 1                                                                       
        end if
        do mu = -1, 1                                                                       ! polarisation
           do ll = 1, l_final_dim                                                            
              lf = l_final(ll)                                                                ! l of the final state
              if( lf == -1 ) exit                                                             ! dummy element 
              if( lminus .and. lf /= nl_init(2) - xang_mom ) cycle                            ! transitions to s states only
              lmi = 0                                                                     
              lmf = 2*lf                                                                  
              do lm = lmi, lmf                                                                ! m of the final state
                 if( mod(lm,2) == 0 ) then                                                     ! Transform lm to m
                    m = -(lm/2)                                                             
                 else                                                                      
                    m = (lm+1)/2                                                            
                 end if
                 llm = lf**2 + 1 + lm                                                          ! indexation in Mxanes is different
                 !
                 ! thi s could be suppressed
                 !
                 if( nl_init(2) == 0 ) then                                                    ! Matteo prefers l0 = 0 treated separately
                    if( .not. q_final(ll) .and. gaunt(lf,m,0,0,1,mu) == 0. ) cycle
                    if( .not. q_final(ll) ) then                                                ! E1-E1
                       do ip=1, paw_recon(xiabs)%paw_nl(lf)                                      ! no of proj.
                          ipx=ipx_0+paw_iltonhb(lf,ip,xiabs)                             
                          term(1:npw) = facE1     * y_eps(mu+2)                               &   ! mu* 
                               * paw_vkb_cplx(1:npw,ipx+lm)                &   ! projector: contains radial part
                               * gaunt(lf,m,0,0,1,mu)                      &   ! the Gaunt coefficient
                               * Mxanes(ip,llm)                                ! the radial integral
                          psiwfc(1:npw) = psiwfc(1:npw) + term(1:npw)                         
                       end do
                    else                                                                        ! E2-E2
                       do lmbd = -1, 1
                          if( gaunt4Y(lf,m,0,0,1,mu,1,lmbd) == 0. ) cycle
                          do ip= 1, paw_recon(xiabs)%paw_nl(lf)                                   ! no of proj.
                             ipx=ipx_0+paw_iltonhb(lf,ip,xiabs)
                             term(1:npw) =  facE2     * y_k(lmbd+2)                            &   ! orientation of the wavevector (E2 only) 
                                  * y_eps(mu+2)                            &   ! mu*                             
                                  * gaunt4Y(lf,m,0,0,1,mu,1,lmbd)          &   ! the generalized Gaunt coefficient
                                  * paw_vkb_cplx(1:npw,ipx+lm)             &   ! projector: contains radial part
                                  * Mxanes(ip,llm)                             ! the radial integral
                             !if (ionode) print*, "proj", paw_vkb_cplx(50,ipx+lm) 
                             !if (ionode) print*, "rad", Mxanes(ip,llm) 
                             psiwfc(1:npw) = psiwfc(1:npw) + term(1:npw)    
                          end do
                       end do
                    end if
                 else
                    do isg = 1, isgf 
                       if( nspin == 2 ) then                         ! no magnetism: sigma_up = sigma_down
                          if( isg == 1 .and. ik > nkstot/2 ) cycle   ! the up/down is defined by nkstot and  
                          if( isg == 2 .and. ik < nkstot/2+1 ) cycle !             not by nks
                       end if
                       ms = nint( mj - isg + 1.5d0 )             ! m of the initial state (int) 
                       if( .not. q_final(ll) .and. gaunt(lf,m,nl_init(2),ms,1,mu) == 0. ) cycle
                       dl  = dble( nl_init(2) )                       ! ClebschG requires dble input values
                       disg = dble(isg) - 1.5d0                       ! the actual spin
                       mss = mj - disg                                ! m of the initial state (dble) 
                       CG = ClebschG(dl,s,j,mss,disg,mj)
                       if( .not. q_final(ll) ) then                    ! E1-E1
                          do ip=1, paw_recon(xiabs)%paw_nl(lf)          ! number of projectors
                             ipx=ipx_0+paw_iltonhb(lf,ip,xiabs)
                             term(1:npw) = facE1   * y_eps(mu+2)                               &   ! mu* 
                                  * gaunt(lf,m,nl_init(2),ms,1,mu)            &   ! the Gaunt coef
                                  * CG                                        &   ! Clebsch Gordan
                                  * paw_vkb_cplx(1:npw,ipx+lm)                &   ! projector: contains radial part
                                  * Mxanes(ip,llm)                                ! the radial integral
                             ! OB: if I do this, I don't get the same result as in magnetic calc with up = down :
                             !             if( nspin == 1 ) term(1:npw) = 2 * term(1:npw) 
                             ! it has to do with the doubling of k points if magnetic 
                             psiwfc(1:npw) = psiwfc(1:npw) + term(1:npw)
                          end do
                       else                                                                      ! E2-E2
                          do lmbd = -1, 1  
                             if( gaunt4Y(lf,m,nl_init(2),ms,1,mu,1,lmbd) == 0. ) cycle
                             do ip=1, paw_recon(xiabs)%paw_nl(lf)                                  ! no of proj. 
                                ipx=ipx_0+paw_iltonhb(lf,ip,xiabs)
                                term(1:npw) = facE2 * y_k(lmbd+2)                               &   ! orientation of the wavevector 
                                     * y_eps(mu+2)                               &   ! mu*
                                     * gaunt4Y(lf,m,nl_init(2),ms,1,mu,1,lmbd)   &   ! the generalized Gaunt coefficient
                                     * CG                                        &   ! Clebsch Gordan
                                     * paw_vkb_cplx(1:npw,ipx+lm)                &   ! projector: contains radial part
                                     * Mxanes(ip,llm)                                ! the radial integral
                                psiwfc(1:npw) = psiwfc(1:npw) + term(1:npw)
                             end do
                          end do
                       end if
                    end do
                 end if
              end do
           end do
        end do
        
        
        
        ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
        ! Starting Lanczos
        ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$  
        
        
        !   
        !   I normalize the wavefunction psiwfc(1:npw)
        !
        
        !<CG>
        CALL allocate_bec_type(nkb,1,becp)
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
        
        xnorm(i_lanczos,ik)=SQRT(xnorm_partial)

        WRITE(stdout,'(8x,a,e15.8)') '|   Norm of the initial Lanczos vector:',&
             xnorm(1,ik)
        
        ! <OB> in the spin unpolarized case at the L3 (M5) edges, the contribution of one of the initial states is null
        !      at all k, if nspin=2, half of the k points should not give any signal
        !          (for sigma up pure down states do not count)
        !      hence the norm is 0
        ! <OB>
        
        ! <OB> for Matteo:  perhaps the 1.d-9 criterion is too strict
        !  MCB: not needed
        !
        if( abs( xnorm(i_lanczos,ik) ) > 1.d-9 ) then
           norm=1.d0/xnorm(i_lanczos,ik)
           
           CALL zdscal(npw,norm,psiwfc,1)
           !
           !      Then I call the lanczos routine
           !

           IF (okvan) THEN
              CALL lanczos_uspp(a(:,i_lanczos,ik),b(:,i_lanczos,ik),npw,psiwfc,ncalcv(i_lanczos,ik), terminator)
           ELSE
              CALL lanczos(a(:,i_lanczos,ik),b(:,i_lanczos,ik),npw,psiwfc,ncalcv(i_lanczos,ik), terminator)
           ENDIF
        else
           ncalcv(i_lanczos,ik) = 0
           nocalc(i_lanczos,ik) = .true. 
      end if
        
        
        CALL deallocate_bec_type ( becp ) ! CG
        calculated(i_lanczos,ik)=1
        
        timenow=get_clock( 'xanes' ) 
        WRITE( stdout,'(" total cpu time spent 4 is ",F9.2," secs")') timenow
        
     ENDDO  !on k points
     CALL mp_barrier(world_comm)
     CALL mp_sum(nunfinished, world_comm)

!     WRITE(6,'(3f16.8)') ((a(jloop,i_lanczos,ik),jloop=1,ncalcv(i_lanczos,ik)),ik=1,1)
!     write(stdout,*)
!     WRITE(6,*) ((b(jloop,i_lanczos,ik),jloop=1,ncalcv(i_lanczos,ik)),ik=1,nkstot)
!     stop
  end do ! on Lanczos
  
  ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  ! Array deallocation
  ! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$  
  
  CALL mp_barrier(world_comm)
  

  IF (nunfinished >= 1) THEN 
     save_file_kind='unfinished'
     write(stdout,*) 'calculation not finished'
  ENDIF
  
  deallocate( l_final )
  deallocate( q_final )
  DEALLOCATE( psiwfc, term )
  DEALLOCATE (paw_vkb_cplx)
  DEALLOCATE( Mxanes )
END SUBROUTINE xanes_dipole_general_edge


