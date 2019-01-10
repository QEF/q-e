!-----------------------------------------------------------------------
!
! Program written by Yang Jiao,  Oct 2016, GPL, No warranties.
!
!    Oct 2018, adapted to QE6.3
!
!-----------------------------------------------------------------------
PROGRAM do_ppacf
  !-----------------------------------------------------------------------
  !
  ! ... This routine computes the coupling constant dependency of
  !     exchange correlation potential
  ! ... E_{xc,\lambda}, \lambda \in [0:1]
  !     and the spatial distribution of exchange correlation energy 
  !     density and kinetic correlation energy density according to
  !
  !     Y. Jiao, E. Schr\"oder, and P. Hyldgaard, 
  !     Phys. Rev. B 97, 085115 (2018).
  !     
  !     For an illustration of how to use this routine to set hybrid 
  !     mixing parameter, please refer to 
  !     
  !     Y. Jiao, E. Schr\"oder, and P. Hyldgaard, 
  !     J. Chem. Phys. 148, 194115 (2018).     
  ! ...
  !
  USE basis,                ONLY : starting_wfc
  USE constants,            ONLY : e2,pi,fpi
  USE kinds,                ONLY : DP
  USE klist,                ONLY : nks, xk, ngk, igk_k
  USE gvect,                ONLY : ngm, g
  USE gvecw,                ONLY : ecutwfc,gcutw
  USE io_files,             ONLY : pseudo_dir,prefix,tmp_dir
  USE io_global,            ONLY : stdout, ionode, ionode_id
  USE cell_base,            ONLY : omega
  USE mp,                   ONLY : mp_bcast, mp_sum
  USE mp_world,             ONLY : world_comm
  USE mp_global,            ONLY : mp_startup
  USE mp_bands,             ONLY : intra_bgrp_comm
  USE exx,                  ONLY : exxinit, exxenergy2, fock2, ecutfock, & 
                                   use_ace, aceinit
  USE exx_base,             ONLY : exx_grid_init, exx_mp_init, exx_div_check,exxdiv_treatment
  USE exx_base,             ONLY : nq1, nq2, nq3
  USE fft_base,             ONLY : dfftp
  USE fft_interfaces,       ONLY : invfft
  USE lsda_mod,             ONLY : nspin
  USE scf,                  ONLY : scf_type,create_scf_type,destroy_scf_type
  USE scf,                  ONLY : scf_type_COPY
  USE scf,                  ONLY : rho, rho_core, rhog_core,vltot
  USE funct,                ONLY : xc,xc_spin,gcxc,gcx_spin,gcc_spin,dft_is_nonlocc,nlc
  USE funct,                ONLY : get_iexch, get_icorr, get_igcx, get_igcc
  USE funct,                ONLY : set_exx_fraction,set_auxiliary_flags,enforce_input_dft
  USE wvfct,                ONLY : npw, npwx
  USE environment,          ONLY : environment_start, environment_end
  USE kernel_table,         ONLY : Nqs, vdw_table_name, kernel_file_name
  USE vdW_DF,               ONLY : get_potential, vdW_energy
  USE vdW_DF_scale,         ONLY : xc_vdW_DF_ncc, xc_vdW_DF_spin_ncc, &
                                   get_q0cc_on_grid, get_q0cc_on_grid_spin
  USE vasp_xml,             ONLY : readxmlfile_vasp

  ! 
  IMPLICIT NONE
  !
  LOGICAL :: lplot,ltks,lfock,lecnl_qxln,lecnl_qx
  INTEGER :: code_num
  !  From which code to read in the calculation data
  !  1 Quantum ESPRESSO (default)
  !  2 VASP
  INTEGER :: icc, ncc
  INTEGER :: n_lambda
  REAL(DP):: rs,rs3,s,q,qx,qc
  REAL(DP):: etxclambda,etx,etxlda,etxgc,etcnlccc,etcnlcccp,etcnlcccm
  REAL(DP) :: etcldalambda,etcgclambda,etcnlclambda, etcnl_check,ttcnl_check, tcnl_int
  REAL(DP), ALLOCATABLE :: Ec_nl_ngamma(:)
  REAL(DP) :: etc,etclda,etcgc
  REAL(dp), EXTERNAL :: exxenergyace
  ! !
  INTEGER :: is, ir,iq,ig,icar, nnrtot
  INTEGER :: iexch,icorr,igcx,igcc,inlc
  ! counter on mesh points
  ! counter on nspin
  INTEGER  :: ierr,ios
  REAL(DP) :: cc,dcc,ccp,ccm,ccp2,ccm2,ccp3,ccm3,ccp4,ccm4,ccp8,ccm8,cc3
  ! coupling constant
  ! local exchange energy, local correlation energy
  ! local exchange potential, local correlation potential
  REAL(DP) :: rhox,arhox
  ! the charge in each point
  ! the absolute value of the charge
  REAL(DP) :: etxc, vtxc
  REAL(DP) :: ex,ec,vx(2),vc(2),expp,ecpp,exm,ecm
  REAL(DP) :: ec_l,ecgc_l,Ec_nl
  REAL(DP) :: etxccc,etxcccnl,etxcccnlp,etxcccnlm,vtxccc,vtxccc_buf,vtxcccnl 
  REAL(DP) :: grho2(2),sx, sc,scp,scm, v1x, v2x, v1c, v2c, &
              v1xup, v1xdw, v2xup, v2xdw, v1cup, v1cdw,  &
              etxcgc, vtxcgc, segno, fac, zeta, rh, grh2, amag 
  real(dp) :: dq0_dq                              ! The derivative of the saturated
  real(dp) :: grid_cell_volume
  REAL(DP), ALLOCATABLE :: q0(:)
  REAL(DP), ALLOCATABLE :: vofrcc(:,:)
  REAL(DP), PARAMETER :: vanishing_charge = 1.D-10, &
                         vanishing_mag    = 1.D-20
  real(DP), parameter :: small = 1.E-10_DP,  third = 1.0_DP / 3.0_DP, &
       pi34 = 0.6203504908994_DP  ! pi34=(3/4pi)^(1/3)
  REAL(DP), PARAMETER :: epsr = 1.D-6, epsg = 1.D-10
  REAL(DP), ALLOCATABLE :: grho(:,:,:),rhoout(:,:)
  COMPLEX(DP), ALLOCATABLE :: rhogsum(:,:)
  REAL(DP), ALLOCATABLE :: tot_grad_rho(:,:),grad_rho(:,:,:)
  REAL(DP), ALLOCATABLE :: tot_rho(:)

  INTEGER :: ik                 ! counter on k points
  INTEGER, ALLOCATABLE  :: igk_buf(:)
  REAL(dp), ALLOCATABLE :: gk(:) ! work space

  CHARACTER (LEN=256) :: outdir

  CHARACTER(len=256)  :: filplot
  INTEGER  :: plot_num

  REAL(DP) :: ttclda        ! the kinetic energy T_c^LDA=E_c,lambda=1^LDA-E_c^LDA
  REAL(DP) :: gtimesr
  REAL(DP)              :: r_shift(3)
  real(dp), allocatable :: dq0_drho(:,:)
  real(dp), allocatable :: dq0_drho_up(:)       ! The derivative of the saturated q0
  real(dp), allocatable :: dq0_drho_down(:)     ! (equation 5 of SOLER) with respect
                                                ! to the charge density (see
                                                ! get_q0_on_grid subroutine for details).

  real(dp), allocatable :: dq0_dgradrho(:,:) 
  real(dp), allocatable :: dq0_dgradrho_up(:)   ! The derivative of the saturated q0
  real(dp), allocatable :: dq0_dgradrho_down(:) ! (equation 5 of SOLER) with respect
                                                ! to the gradient of the charge density
                                                ! (again, see get_q0_on_grid subroutine).
  real(dp), allocatable :: potential_vdW(:,:)   ! The vdW contribution to the potential
  complex(dp), allocatable :: u_vdW(:,:)        ! 
  COMPLEX(DP), ALLOCATABLE :: up_vdW(:,:), um_vdW(:,:)
  complex(dp), allocatable :: thetas(:,:)       ! These are the functions of equation 8 of
                                                ! SOLER. They will be forward Fourier transformed
                                                ! in place to get theta(k) and worked on in
                                                ! place to get the u_alpha(r) of equation 11
                                                ! in SOLER. They are formatted as follows:
                                                ! thetas(grid_point, theta_i).
  COMPLEX(DP), ALLOCATABLE :: thetasp(:,:), thetasm(:,:)
  COMPLEX(DP), ALLOCATABLE :: ecnl_c(:), tcnl_c(:)

  REAL(DP), ALLOCATABLE    :: kin_r(:,:)        ! the kinetic energy density in R-space

  type (scf_type) :: exlda, eclda
  type (scf_type) :: tclda  ! the kinetic energy per particle
  type (scf_type) :: ecnl   ! the non-local correlation energy per partical
  type (scf_type) :: tcnl
  type (scf_type) :: exgc, ecgc, tcgc

  NAMELIST / ppacf / code_num,outdir,prefix,n_lambda,lplot,ltks,lfock,use_ace, &
                     pseudo_dir,vdw_table_name,lecnl_qxln,lecnl_qx

  !
  ! initialise environment
  !
#ifdef __MPI
  CALL mp_startup ()
#endif
!!!!!!!!!!!!!!! READ IN PREFIX !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  CALL environment_start ('ppacf')
  !
  ! set default values for variables in namelist
  !
  code_num = 1
  outdir = './'
  prefix = 'ppacf'
  n_lambda = 1
  lplot = .False.
  ltks  = .False.
  lfock  = .False.
  use_ace = .True.
  dcc = 1.E-6_DP 
  filplot = 'tmp.pp'
  plot_num = -1
  lecnl_qxln=.False.
  lecnl_qx=.False.
  !
  
  IF (ionode) THEN
     !
     CALL input_from_file ()
     !
     READ (5, ppacf, iostat = ios) 
     !
  ENDIF
!!!!!!!!!!!!!!!!!! READ IN DATA !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! 
  ! Broadcast variables
  !
  CALL mp_bcast(code_num,ionode_id,world_comm)
  CALL mp_bcast(outdir,ionode_id,world_comm)
  CALL mp_bcast(prefix,ionode_id,world_comm)
  CALL mp_bcast(n_lambda,ionode_id,world_comm)
  CALL mp_bcast(lplot,ionode_id,world_comm)
  CALL mp_bcast(ltks,ionode_id,world_comm)
  CALL mp_bcast(lfock,ionode_id,world_comm)
  CALL mp_bcast(lecnl_qxln,ionode_id,world_comm)
  CALL mp_bcast(lecnl_qx,ionode_id,world_comm)
  CALL mp_bcast(dcc,ionode_id,world_comm)
  CALL mp_bcast(pseudo_dir,ionode_id,world_comm)
  CALL mp_bcast(vdw_table_name,ionode_id,world_comm)
  ncc=n_lambda 
  WRITE( stdout, '(//5x,"entering subroutine acf ..."/)')
  ! Write out the ppacf information.
  IF(ionode) CALL ppacf_info()
  !
  CALL start_clock( 'acf_etxclambda' )
  ! 
!  WRITE(stdout,9093) dcc

  IF (code_num == 1) THEN
     !
     tmp_dir=TRIM(outdir) 
!     CALL read_xml_file_internal(.TRUE.)
     CALL  read_file()

!     Check exchange correlation functional
     iexch = get_iexch()
     icorr = get_icorr()
     igcx  = get_igcx()
     igcc  = get_igcc()
  
  ELSEIF (code_num == 2) THEN
     !
     tmp_dir=TRIM(outdir)
     CALL readxmlfile_vasp(iexch,icorr,igcx,igcc,inlc,ierr)
     IF(ionode) WRITE(stdout,'(5X,a)') "Read data from VASP output 'vasprun.xml'"
     ! 
  ELSE
     CALL errore ('ppacf', 'code_num not implemented', 1)
  ENDIF
  
  ALLOCATE(vofrcc(1:dfftp%nnr,1:nspin))
  
  IF(dft_is_nonlocc()) THEN
     WRITE( stdout, '(//5x,"ACF coupling-constant  Exc_lambda (Ry)  E_c,lambda^LDA (Ry)  E_c,lambda^nl (Ry)"/)')
  ELSE
     WRITE( stdout, '(//5x,"ACF coupling-constant  Exc_lambda (Ry)  E_c,lambda^LDA (Ry)  E_c,lambda^GC (Ry)"/)')
  END IF
  !
  IF (dft_is_nonlocc() .AND. lecnl_qxln) THEN
     etxcccnl=0._DP
     vtxcccnl=0._DP
     vofrcc=0._DP
     CALL nlc( rho%of_r, rho_core, nspin, etxcccnl, vtxcccnl, vofrcc )
     CALL mp_sum(  etxcccnl , intra_bgrp_comm )
  END IF
  !
   
  !
  ! ... add gradiend corrections (if any)
  !
  if (nspin==4) CALL errore ('ppacf', 'Noncollinear not implemented', 1)
  fac = 1.D0 / DBLE( nspin )
  !
  ALLOCATE( grho( 3, dfftp%nnr, nspin) )
  ALLOCATE( rhoout( dfftp%nnr, nspin) )
  ALLOCATE( tot_rho( dfftp%nnr) )
  ALLOCATE( rhogsum( ngm, nspin ) )
  !
  ! ... calculate the gradient of rho + rho_core in real space
  !
  !
  !
  DO is = 1, nspin
     !
     rhoout(:,is)  = fac * rho_core(:)  + rho%of_r(:,is)  
     rhogsum(:,is) = fac * rhog_core(:) + rho%of_g(:,is)
     !
     CALL fft_gradient_g2r( dfftp, rhogsum(1,is), g, grho(1,1,is))
     !
  END DO
  !
  DEALLOCATE( rhogsum )
  !
  IF (nspin == 1) THEN
     tot_rho(:)=rhoout(:,1)
  ELSEIF(nspin==2) THEN
     tot_rho(:)=rhoout(:,1)+rhoout(:,2)
  ELSE
     CALL errore ('ppacf','vdW-DF not available for noncollinear spin case',1)
  END IF
  !
  CALL create_scf_type(exlda)
  exlda%of_r(:,:)=0._DP
  CALL create_scf_type(eclda)
  eclda%of_r(:,:)=0._DP
  CALL create_scf_type(tclda)
  tclda%of_r(:,:)=0._DP
  CALL create_scf_type(exgc)
  exgc%of_r(:,:)=0._DP
  IF (dft_is_nonlocc()) THEN
     CALL create_scf_type(ecnl)
     CALL create_scf_type(tcnl)
     ecnl%of_r(:,:)=0._DP
     tcnl%of_r(:,:)=0._DP
     ALLOCATE(Ec_nl_ngamma(1:ncc))
     Ec_nl_ngamma=0._DP
  ELSEIF(igcc .NE. 0) THEN
     CALL create_scf_type(ecgc)
     CALL create_scf_type(tcgc)
     ecgc%of_r(:,:)=0._DP
     tcgc%of_r(:,:)=0._DP
  END IF
  ttclda = 0._DP
  !
  !! coupling constant > 0
  ! 
  DO icc = 0, ncc
     cc=DBLE(icc)/DBLE(ncc)
     etxclambda=0._DP
     etxccc=0._DP
     vofrcc=0._DP
     etx=0._DP
     etxlda=0._DP
     etxgc=0._DP
     etcldalambda=0._DP
     etcgclambda=0._DP
     etc=0._DP
     etclda=0._DP
     etcgc=0._DP
     !
     IF (nspin == 1 ) THEN
     !
     ! ... spin-unpolarized case
     !
     DO ir = 1, dfftp%nnr
        !
        rhox = rho%of_r(ir,1) + rho_core(ir)
        arhox = ABS(rhox)
        IF (arhox > vanishing_charge) THEN
           IF(iexch==1) THEN
              rs = pi34 /arhox**third
              CALL slater(rs,ex,vx(1))     ! \epsilon_x,\lambda[n]=\epsilon_x[n]
           ELSE
              CALL xc( arhox, ex, ec, vx(1), vc(1) )
           ENDIF
           etx=etx+e2*ex*rhox
           etxlda=etxlda+e2*ex*rhox
           grho2(1) = grho(1,ir,1)**2 + grho(2,ir,1)**2 + grho(3,ir,1)**2
           IF(cc>0._DP) THEN
                 ccp=cc+dcc
                 ccm=cc-dcc
                 ccp2=ccp*ccp
                 ccp3=ccp2*ccp
                 ccp4=ccp3*ccp
                 ccp8=ccp4*ccp4
                 ccm2=ccm*ccm
                 ccm3=ccm2*ccm
                 ccm4=ccm3*ccm
                 ccm8=ccm4*ccm4
              IF(icorr==4) THEN
                 CALL pwcc (rs,cc,ec,vc(1),ec_l)
              ELSE
                 CALL xc(arhox/ccp3,expp,ecpp,vx(1),vc(1))
                 CALL xc(arhox/ccm3,exm,ecm,vx(1),vc(1))
                 ec_l=(ccp2*ecpp-ccm2*ecm)/dcc*0.5_DP
              ENDIF
              etcldalambda=etcldalambda+e2*ec_l*rhox
              IF(icc == ncc) THEN
                 IF(icorr.NE.4) CALL xc(arhox,ex,ec,vx(1),vc(1))
                 tclda%of_r(ir,1)=e2*(ec-ec_l)*rhox
                 ttclda=ttclda+e2*(ec-ec_l)*rhox
              ENDIF
              IF(grho2(1) > epsg .AND. igcc .NE. 0) THEN
                 segno = SIGN( 1.D0, rhoout(ir,1) )
                 CALL gcxc(arhox/ccp3,grho2(1)/ccp8,sx,scp,v1x,v2x,v1c,v2c)
                 CALL gcxc(arhox/ccm3,grho2(1)/ccm8,sx,scm,v1x,v2x,v1c,v2c)
                 ecgc_l=(ccp2*scp*ccp3-ccm2*scm*ccm3)/dcc*0.5_DP
                 etcgclambda=etcgclambda+e2*ecgc_l*segno
              ENDIF
           ENDIF
           CALL xc(arhox,ex,ec,vx(1),vc(1))
           etclda=etclda+e2*ec*rhox
           etc=etc+e2*ec*rhox
           IF(icc==ncc) THEN
              exlda%of_r(ir,1)=e2*ex*rhox
              eclda%of_r(ir,1)=e2*ec*rhox
           END IF
           IF ( grho2(1) > epsg ) THEN
              segno = SIGN( 1.D0, rhoout(ir,1) )
              CALL gcxc( arhox, grho2(1), sx, sc, v1x, v2x, v1c, v2c )
              etx=etx+e2*sx*segno
              etxgc=etxgc+e2*sx*segno
              etc=etc+e2*sc*segno
              etcgc=etcgc+e2*sc*segno
              IF(icc==ncc) THEN
                 exgc%of_r(ir,1)=e2*sx*segno
                 IF(igcc.NE.0) THEN
                    ecgc%of_r(ir,1)=e2*sc*segno
                    tcgc%of_r(ir,1)=e2*(sc-ecgc_l)*segno
                 END IF
              ENDIF
           ENDIF
        ENDIF
     ENDDO
     ELSE IF (nspin ==2) THEN
     !
     ! ... spin-polarized case
     !
     DO ir = 1, dfftp%nnr
        rhox = rho%of_r(ir,1) + rho%of_r(ir,2) + rho_core(ir)
        arhox = ABS( rhox )
        IF (arhox > vanishing_charge) THEN
           rs = pi34 /arhox**third
           zeta = (rho%of_r(ir,1)-rho%of_r(ir,2))/arhox
           IF( ABS( zeta ) > 1.D0 ) zeta = SIGN(1.D0, zeta)
           IF(iexch==1) THEN
              CALL slater_spin (arhox, zeta, ex, vx(1), vx(2))
           ELSE
              CALL xc_spin( arhox, zeta, ex, ec, vx(1), vx(2), vc(1), vc(2) )
           ENDIF
           etx=etx+e2*ex*rhox
           etxlda=etxlda+e2*ex*rhox
           grh2 = ( grho(1,ir,1) + grho(1,ir,2) )**2 + &
                  ( grho(2,ir,1) + grho(2,ir,2) )**2 + &
                  ( grho(3,ir,1) + grho(3,ir,2) )**2
           IF(cc>0._DP) THEN
                 ccp=cc+dcc
                 ccm=cc-dcc
                 ccp2=ccp*ccp
                 ccp3=ccp2*ccp
                 ccp4=ccp3*ccp
                 ccp8=ccp4*ccp4
                 ccm2=ccm*ccm
                 ccm3=ccm2*ccm
                 ccm4=ccm3*ccm
                 ccm8=ccm4*ccm4
              IF(icorr==4) THEN
                 CALL pwcc_spin (rs,cc, zeta, ec, vc(1), vc(2), ec_l)
              ELSE
                 CALL xc_spin( arhox/ccp3, zeta, expp, ecpp, vx(1), vx(2), vc(1), vc(2) )
                 CALL xc_spin( arhox/ccm3, zeta, exm, ecm, vx(1), vx(2), vc(1), vc(2) )
                 ec_l=(ccp2*ecpp-ccm2*ecm)/dcc*0.5_DP
              ENDIF
              etcldalambda=etcldalambda+e2*ec_l*rhox
              IF(icc == ncc) THEN
                 tclda%of_r(ir,1)=e2*(ec-ec_l)*rhox
                 ttclda=ttclda+e2*(ec-ec_l)*rhox
              ENDIF
              IF(igcc .NE. 0) THEN
                 CALL gcc_spin(rhox/ccp3,zeta,grh2/ccp8,scp,v1cup,v1cdw,v2c)
                 CALL gcc_spin(rhox/ccm3,zeta,grh2/ccm8,scm,v1cup,v1cdw,v2c)
                 ecgc_l=(ccp2*scp*ccp3-ccm2*scm*ccm3)/dcc*0.5_DP
                 etcgclambda=etcgclambda+e2*ecgc_l
              ENDIF
           ENDIF
           CALL xc_spin( arhox, zeta, ex, ec, vx(1), vx(2), vc(1), vc(2) )
           etclda=etclda+e2*ec*rhox
           etc=etc+e2*ec*rhox
           IF(icc==ncc) THEN
              exlda%of_r(ir,1)=e2*ex*rhox
              eclda%of_r(ir,1)=e2*ec*rhox
           END IF
           grho2(:) = grho(1,ir,:)**2 + grho(2,ir,:)**2 + grho(3,ir,:)**2
           CALL gcx_spin( rhoout(ir,1), rhoout(ir,2), grho2(1), &
                          grho2(2), sx, v1xup, v1xdw, v2xup, v2xdw )
           etx=etx+e2*sx
           etxgc=etxgc+e2*sx
           CALL gcc_spin( rhox, zeta, grh2, sc, v1cup, v1cdw, v2c )
           etcgc=etcgc+e2*sc
           etc=etc+e2*sc
           IF(icc==ncc) THEN
              exgc%of_r(ir,1)=e2*sx
              IF(igcc.NE.0) THEN
                 ecgc%of_r(ir,1)=e2*sc
                 tcgc%of_r(ir,1)=e2*(sc-ecgc_l)
              END IF
           END IF
        END IF
     END DO
     ELSE IF (nspin==4) THEN
     CALL errore ('ppacf', 'Noncollinear not implemented', 1)
     END IF
  ! 
  grid_cell_volume=omega/DBLE(dfftp%nr1*dfftp%nr2*dfftp%nr3)
  etx = omega*etx/DBLE(dfftp%nr1*dfftp%nr2*dfftp%nr3)
  etc = omega*etc/DBLE(dfftp%nr1*dfftp%nr2*dfftp%nr3)
  etxlda = omega*etxlda/DBLE(dfftp%nr1*dfftp%nr2*dfftp%nr3)
  etclda = omega*etclda/DBLE(dfftp%nr1*dfftp%nr2*dfftp%nr3)
  etxgc = omega*etxgc/DBLE(dfftp%nr1*dfftp%nr2*dfftp%nr3)
  etcgc = omega*etcgc/DBLE(dfftp%nr1*dfftp%nr2*dfftp%nr3)
  etcldalambda=omega*etcldalambda/DBLE(dfftp%nr1*dfftp%nr2*dfftp%nr3)
  etcgclambda=omega*etcgclambda/DBLE(dfftp%nr1*dfftp%nr2*dfftp%nr3)
  ttclda=omega*ttclda/DBLE(dfftp%nr1*dfftp%nr2*dfftp%nr3)
  CALL mp_sum(  etx , intra_bgrp_comm )
  CALL mp_sum(  etc , intra_bgrp_comm )
  CALL mp_sum(  etxlda , intra_bgrp_comm )
  CALL mp_sum(  etclda , intra_bgrp_comm )
  CALL mp_sum(  etxgc , intra_bgrp_comm )
  CALL mp_sum(  etcgc , intra_bgrp_comm )
  CALL mp_sum(  etcldalambda , intra_bgrp_comm )
  CALL mp_sum(  etcgclambda, intra_bgrp_comm )
  CALL mp_sum(  ttclda, intra_bgrp_comm )
  !
  !
! Non-local correlation 
  etcnlclambda=0._DP
! Approximation of exchange linear lambda-dependence
  IF (dft_is_nonlocc()) THEN
    IF (lecnl_qxln) THEN
      etcnlclambda = 2._DP*cc*etxcccnl
! Full lambda dependence
    ELSE  
      IF(cc>0._DP) THEN
        ccp=cc+dcc
        ccm=cc-dcc
        IF(nspin .EQ. 1) THEN
          CALL xc_vdW_DF_ncc(cc,lecnl_qx,etcnlccc)
          CALL xc_vdW_DF_ncc(ccp,lecnl_qx,etcnlcccp) 
          CALL xc_vdW_DF_ncc(ccm,lecnl_qx,etcnlcccm) 
          etcnlclambda=2._DP*cc*etcnlccc+0.5_DP*cc*cc*(etcnlcccp-etcnlcccm)/dcc
          Ec_nl_ngamma(icc)=etcnlccc
        ELSE IF(nspin .EQ. 2) THEN
          CALL xc_vdW_DF_spin_ncc(cc,lecnl_qx,etcnlccc) 
          CALL xc_vdW_DF_spin_ncc(ccp,lecnl_qx,etcnlcccp) 
          CALL xc_vdW_DF_spin_ncc(ccm,lecnl_qx,etcnlcccm) 
          etcnlclambda=2._DP*cc*etcnlccc+0.5_DP*cc*cc*(etcnlcccp-etcnlcccm)/dcc
          Ec_nl_ngamma(icc)=etcnlccc
        END IF
      END IF
    END IF
    etc=etc+etcnlccc 
  END IF

  
!  !
  etxclambda=etx+etcldalambda+etcgclambda+etcnlclambda
  IF(dft_is_nonlocc()) THEN
      WRITE(stdout,9091) cc, etxclambda,etcldalambda,etcnlclambda
  ELSE
      WRITE(stdout,9092) cc, etxclambda,etcldalambda,etcgclambda
  END IF
  FLUSH(stdout)
  !

  END DO  ! icc
  
  IF(dft_is_nonlocc()) THEN
    WRITE(stdout,'(5x,a)') 'Ec_nl(n_1/lambda): '
    DO icc = 1, ncc
      cc=DBLE(icc)/DBLE(ncc)
      WRITE(stdout,9095) cc, Ec_nl_ngamma(icc)
    END DO
  ENDIF

  WRITE(stdout,'(a32,0PF17.8,a3)') 'Exchange', etx, 'Ry'   !,etxlda,etxgc
  WRITE(stdout,'(a32,0PF17.8,a3)') 'Correlation', etc, 'Ry'   !,etclda,etcgc
  WRITE(stdout,'(a32,0PF17.8,a3)') 'Exchange + Correlation', etx+etc, 'Ry'
  WRITE(stdout,'(a32,0PF17.8,a3)') 'T_c^LDA', ttclda, 'Ry'
  DEALLOCATE(vofrcc)

! Non-local correlation energy density
  IF(dft_is_nonlocc() .AND. lplot) THEN
     ALLOCATE( q0(dfftp%nnr), dq0_drho(dfftp%nnr,nspin), dq0_dgradrho(dfftp%nnr,nspin))
     ALLOCATE( tot_grad_rho(3,dfftp%nnr), grad_rho(3,dfftp%nnr,nspin))
     ALLOCATE( thetas(dfftp%nnr, Nqs) )
     ALLOCATE( u_vdW(dfftp%nnr,Nqs), potential_vdW(dfftp%nnr,nspin))
     ALLOCATE( thetasp(dfftp%nnr,Nqs), thetasm(dfftp%nnr,Nqs)) 
     ALLOCATE( up_vdW(dfftp%nnr,Nqs), um_vdW(dfftp%nnr,Nqs)) 
     ! Here we calculate the gradient in reciprocal space using FFT.
     !
     CALL fft_gradient_r2r (dfftp, tot_rho, g, tot_grad_rho)
     DO is =1, nspin
        CALL fft_gradient_r2r (dfftp, rhoout(:,is), g, grad_rho(:,:,is))
     ENDDO
     !
     dq0_drho=0._DP
     dq0_dgradrho=0._DP
     ccp=1._DP+dcc
     ccm=1._DP-dcc
     !
     ! thetas
     !
     IF(nspin==1) THEN
        CALL get_q0cc_on_grid (1._DP,lecnl_qx,tot_rho, tot_grad_rho, q0, thetas)
        CALL get_q0cc_on_grid (ccp,lecnl_qx,tot_rho, tot_grad_rho, q0, thetasp)
        CALL get_q0cc_on_grid (ccm,lecnl_qx,tot_rho, tot_grad_rho, q0, thetasm)
     ELSEIF(nspin==2) THEN
        CALL get_q0cc_on_grid_spin (1._DP,lecnl_qx,tot_rho, rhoout(:,1), rhoout(:,2), &
          tot_grad_rho, grad_rho(:,:,1), grad_rho(:,:,2), q0, thetas)
        CALL get_q0cc_on_grid_spin (ccp,lecnl_qx,tot_rho, rhoout(:,1), rhoout(:,2), &
          tot_grad_rho, grad_rho(:,:,1), grad_rho(:,:,2), q0, thetasp)
        CALL get_q0cc_on_grid_spin (ccm,lecnl_qx,tot_rho, rhoout(:,1), rhoout(:,2), &
          tot_grad_rho, grad_rho(:,:,1), grad_rho(:,:,2), q0, thetasm)
     ELSE
     ENDIF
     u_vdW(:,:)=thetas(:,:)
     up_vdW(:,:)=thetasp(:,:)
     um_vdW(:,:)=thetasm(:,:)
     DO iq=1, Nqs
        CALL invfft('Rho',thetas(:,iq),dfftp)
        CALL invfft('Rho',thetasp(:,iq),dfftp)
        CALL invfft('Rho',thetasm(:,iq),dfftp)
     ENDDO
     CALL vdW_energy(u_vdW,Ec_nl)
     CALL mp_sum(Ec_nl,intra_bgrp_comm)
     write(stdout,*) '     Non-local energy : ', Ec_nl
     CALL vdW_energy(up_vdW,Ec_nl)
     CALL vdW_energy(um_vdW,Ec_nl)
     do iq= 1, Nqs
        CALL invfft('Rho', u_vdW(:,iq), dfftp)
        CALL invfft('Rho', up_vdW(:,iq), dfftp)
        CALL invfft('Rho', um_vdW(:,iq), dfftp)
     end do
     IF( nspin == 1 ) THEN
        CALL get_potential (q0,dq0_drho(:,1), dq0_dgradrho(:,1), tot_grad_rho, u_vdW, potential_vdW(:,1)) 
     ELSEIF( nspin == 2 ) THEN
        CALL get_potential (q0,dq0_drho(:,1), dq0_dgradrho(:,1), grad_rho(:,:,1), u_vdW, potential_vdW(:,1)) 
        CALL get_potential (q0,dq0_drho(:,2), dq0_dgradrho(:,2), grad_rho(:,:,2), u_vdW, potential_vdW(:,2)) 
     END IF

     ecnl%of_r(:,1)=0._DP
     tcnl%of_r(:,1)=0._DP
     ALLOCATE(ecnl_c(dfftp%nnr))
     ALLOCATE(tcnl_c(dfftp%nnr))

     ecnl_c=DCMPLX(0._DP,0._DP)
     tcnl_c=DCMPLX(0._DP,0._DP)
     etcnl_check=0._DP
     ttcnl_check=0._DP
     DO ir=1,dfftp%nnr
        arhox = ABS(tot_rho(ir))
        IF (arhox > vanishing_charge) THEN
           DO iq=1,Nqs
              ecnl_c(ir)=ecnl_c(ir)+  &
                thetas(ir,iq)*u_vdW(ir,iq)
              tcnl_c(ir)=tcnl_c(ir) - thetas(ir,iq)*u_vdW(ir,iq)  & 
                -0.5_DP/dcc*(thetasp(ir,iq)*up_vdW(ir,iq)-thetasm(ir,iq)*um_vdW(ir,iq))
           ENDDO
        etcnl_check=etcnl_check+ecnl_c(ir)
        ttcnl_check=ttcnl_check+tcnl_c(ir)
        ENDIF
     ENDDO
     etcnl_check=e2*0.5_DP*omega*etcnl_check/DBLE(dfftp%nr1*dfftp%nr2*dfftp%nr3)
     ttcnl_check=e2*0.5_DP*omega*ttcnl_check/DBLE(dfftp%nr1*dfftp%nr2*dfftp%nr3)
     CALL mp_sum(etcnl_check,intra_bgrp_comm)
     CALL mp_sum(ttcnl_check,intra_bgrp_comm)
     DEALLOCATE( q0,dq0_drho,dq0_dgradrho,thetas)
     DEALLOCATE( tot_grad_rho, grad_rho)
     DEALLOCATE( u_vdW)
     DEALLOCATE( thetasp, thetasm, up_vdW, um_vdW)
     ecnl%of_r(:,1)=e2*0.5_DP*DBLE(ecnl_c(:))
     tcnl%of_r(:,1)=e2*0.5_DP*DBLE(tcnl_c(:))
     DEALLOCATE(ecnl_c,tcnl_c)
     write(stdout,*) '     Summation of ecnl: ', etcnl_check
     write(stdout,*) '     Summation of tcnl: ', ttcnl_check
  ENDIF


  IF(lfock .OR. (lplot .AND. ltks)) THEN
     IF (code_num==1) THEN
        starting_wfc='file'
        CALL wfcinit()
     ELSE IF (code_num==2) THEN
        CALL errore( 'ppacf', 'wavefunction not implemented for VASP postprocessing', 1)
     ELSE
        CALL errore ('ppacf', 'code_num not implemented', 1)
     END IF
  END IF
! Fock exchange energy from readin wavefunctions

  IF(lfock) THEN
  nq1=0
  nq2=0
  nq3=0
  exxdiv_treatment="gygi-baldereschi"
  ecutfock=ecutwfc
  CALL set_exx_fraction(1._DP)
  CALL enforce_input_dft('HF')
  CALL set_auxiliary_flags
  !
  ALLOCATE (igk_buf(npwx), gk(npwx) )
  igk_k(:,:) = 0
  DO ik = 1, nks
     !
     CALL gk_sort( xk(1,ik), ngm, g, gcutw, npw, igk_buf, gk )
     ngk(ik) = npw
     igk_k(1:npw,ik)= igk_buf(1:npw)
     !
  END DO
  DEALLOCATE ( igk_buf, gk )
  !
!  CALL setup()
  CALL exx_grid_init()
  CALL exx_mp_init()
  CALL exx_div_check()
!!  CALL init_run()
  CALL exxinit(.FALSE.)
  
  IF ( use_ace) THEN
     CALL aceinit ( ) 
     fock2 = exxenergyace()
  ELSE
     fock2 = exxenergy2()
  ENDIF
  WRITE(stdout, 9068) 0.5_DP*fock2
9068 FORMAT( '     Fock energy               =',0PF17.8,' Ry' )
  END IF


! output data in 3D
  IF (lplot) THEN
     IF (code_num==2) THEN
        IF(nspin==1) THEN
           filplot=trim(prefix)//'.chg'
           plot_num=2
           CALL dcopy(dfftp%nnr, rho%of_r(:,1), 1, vltot, 1)
           CALL punch_plot (filplot, plot_num, 0.,0.,0.,0.,0.,0,0,0, .False.)
        ELSEIF(nspin==2) THEN
           filplot=trim(prefix)//'.chg1'
           plot_num=2
           CALL dcopy(dfftp%nnr, rho%of_r(:,1), 1, vltot, 1)
           CALL punch_plot (filplot, plot_num, 0.,0.,0.,0.,0.,0,0,0, .False.)
           filplot=trim(prefix)//'.chg2'
           plot_num=2
           CALL dcopy(dfftp%nnr, rho%of_r(:,2), 1, vltot, 1)
           CALL punch_plot (filplot, plot_num, 0.,0.,0.,0.,0.,0,0,0, .False.)
        END IF
     END IF

     IF (ltks) THEN
        ALLOCATE(kin_r(dfftp%nnr,nspin))
!        CALL init_run()
        CALL sum_band_kin(kin_r)

        IF(nspin==1) THEN
           filplot=trim(prefix)//'.tks'
           plot_num=2
           CALL dcopy(dfftp%nnr, kin_r(:,1), 1, vltot, 1)
           CALL punch_plot (filplot, plot_num, 0.,0.,0.,0.,0.,0,0,0, .False.)
        ELSEIF(nspin==2) THEN
           filplot=trim(prefix)//'.tks1'
           plot_num=2
           CALL dcopy(dfftp%nnr, kin_r(:,1), 1, vltot, 1)
           CALL punch_plot (filplot, plot_num, 0.,0.,0.,0.,0.,0,0,0, .False.)
           filplot=trim(prefix)//'.tks2'
           plot_num=2
           CALL dcopy(dfftp%nnr, kin_r(:,2), 1, vltot, 1)
           CALL punch_plot (filplot, plot_num, 0.,0.,0.,0.,0.,0,0,0, .False.)
        END IF
        DEALLOCATE(kin_r)
     END IF

     filplot=trim(prefix)//'.exlda'
     plot_num=2
     CALL dcopy (dfftp%nnr, exlda%of_r(:,1), 1, vltot, 1)
     CALL punch_plot (filplot, plot_num, 0.,0.,0.,0.,0.,0,0,0, .False.)
     !
     filplot=trim(prefix)//'.eclda'
     plot_num=2
     CALL dcopy (dfftp%nnr, eclda%of_r(:,1), 1, vltot, 1)
     CALL punch_plot (filplot, plot_num, 0.,0.,0.,0.,0.,0,0,0, .False.)
     !
     filplot=trim(prefix)//'.tclda'
     plot_num=2
     CALL dcopy (dfftp%nnr, tclda%of_r(:,1), 1, vltot, 1)
     CALL punch_plot (filplot, plot_num, 0.,0.,0.,0.,0.,0,0,0, .False.)
     !
     IF(igcx.NE.0) THEN
        filplot=trim(prefix)//'.exgc'
        plot_num=2
        CALL dcopy (dfftp%nnr, exgc%of_r(:,1), 1, vltot, 1)
        CALL punch_plot (filplot, plot_num, 0.,0.,0.,0.,0.,0,0,0,.False.)
        CALL destroy_scf_type(exgc)
     ENDIF
     IF(dft_is_nonlocc()) THEN
        filplot=trim(prefix)//'.ecnl'
        plot_num=2
        CALL dcopy (dfftp%nnr, ecnl%of_r(:,1), 1, vltot, 1)
        CALL punch_plot (filplot, plot_num, 0.,0.,0.,0.,0.,0,0,0,.False.)
        !
        filplot=trim(prefix)//'.tcnl'
        plot_num=2
        CALL dcopy (dfftp%nnr, tcnl%of_r(:,1), 1, vltot, 1)
        CALL punch_plot (filplot, plot_num, 0.,0.,0.,0.,0.,0,0,0,.False.)
        !
        IF(nspin==1) THEN
           filplot=trim(prefix)//'.vcnl'
           plot_num=2
           CALL dcopy (dfftp%nnr, potential_vdW(:,1), 1, vltot, 1)
           CALL punch_plot (filplot, plot_num, 0.,0.,0.,0.,0.,0,0,0,.False.)
        ELSEIF(nspin==2) THEN
           filplot=trim(prefix)//'.vcnl1'
           plot_num=2
           CALL dcopy (dfftp%nnr, potential_vdW(:,1), 1, vltot, 1)
           CALL punch_plot (filplot, plot_num, 0.,0.,0.,0.,0.,0,0,0,.False.)
           filplot=trim(prefix)//'.vcnl2'
           plot_num=2
           CALL dcopy (dfftp%nnr, potential_vdW(:,2), 1, vltot, 1)
           CALL punch_plot (filplot, plot_num, 0.,0.,0.,0.,0.,0,0,0,.False.)
        END IF
        DEALLOCATE(potential_vdW)
        CALL destroy_scf_type(ecnl)
        CALL destroy_scf_type(tcnl)
     ELSEIF(igcc.NE.0) THEN
        filplot=trim(prefix)//'.ecgc'
        plot_num=2
        CALL dcopy (dfftp%nnr, ecgc%of_r(:,1), 1, vltot, 1)
        CALL punch_plot (filplot, plot_num, 0.,0.,0.,0.,0.,0,0,0,.False.)
        filplot=trim(prefix)//'.tcgc'
        plot_num=2
        CALL dcopy (dfftp%nnr, tcgc%of_r(:,1), 1, vltot, 1)
        CALL punch_plot (filplot, plot_num, 0.,0.,0.,0.,0.,0,0,0,.False.)
        CALL destroy_scf_type(ecgc)
        CALL destroy_scf_type(tcgc)
     END IF
  !
     CALL stop_clock( 'acf_etxclambda' )
     WRITE( stdout, '(//5x,"exiting subroutine acf ..."/)')
  END IF
  !
  
  CALL destroy_scf_type(exlda)
  CALL destroy_scf_type(eclda)
  CALL destroy_scf_type(tclda)
  !
9091 FORMAT( 0PF17.8,0PF17.8,0PF17.8,0PF17.8,0PF17.8,0PF17.8 )
9092 FORMAT( 0PF17.8,0PF17.8,0PF17.8,0PF17.8 )
9095 FORMAT( 0PF17.8,0PF17.8)
!9090 FORMAT(/' ACF coupling-constant         =',0PF17.8,'   ' &
!            /'     Exc_lambda                =',0PF17.8,' Ry' )
!9091 FORMAT(/'     Non-local contribution    =',0PF17.8,' Ry' )
9093 FORMAT(/' delta coupling constant        =',0PE17.4E3,' ')
  CALL environment_end ('ppacf')
  CALL stop_pp
!-----------------------------------------------------------------------


END PROGRAM do_ppacf
!
!-----------------------------------------------------------------------
subroutine pwcc (rs,cc,ec, vc, ec_l)
  !-----------------------------------------------------------------------
  !     This subroutine is adapted from subroutine pw from 
  !     Modules/functionals.f90
  !
  !     iflag=1: J.P. Perdew and Y. Wang, PRB 45, 13244 (1992)
  !
  USE kinds, ONLY : DP
  implicit none
  real(dp), intent(in) :: rs, cc
  real(dp), intent(out):: ec, vc, ec_l
  integer  :: iflag
  real(DP) :: a, b1, b2, c0, c1, c2, c3, d0, d1
  parameter (a = 0.031091d0, b1 = 7.5957d0, b2 = 3.5876d0, c0 = a, &
       c1 = 0.046644d0, c2 = 0.00664d0, c3 = 0.01043d0, d0 = 0.4335d0, &
       d1 = 1.4408d0)
  real(DP) :: lnrs, rs12, rs32, rs2, om, dom, olog
  real(DP) :: lrs,lrs12,lrs32,lrs2,lom,dlom,lolog
  real(DP) :: lec,dlec
  real(DP) :: a1 (2), b3 (2), b4 (2)
  data a1 / 0.21370d0, 0.026481d0 /, b3 / 1.6382d0, -0.46647d0 /, &
       b4 / 0.49294d0, 0.13354d0 /

  iflag = 1
  ! interpolation formula
  rs12 = sqrt (rs)
  rs32 = rs * rs12
  rs2 = rs**2
  lrs = cc*rs
  lrs12 = sqrt(lrs)
  lrs32 = lrs*lrs12
  lrs2 = lrs**2
  om = 2.d0 * a * (b1 * rs12 + b2 * rs + b3 (iflag) * rs32 + b4 ( &
       iflag) * rs2)
  lom = 2.d0 * a * (b1 * lrs12 + b2 * lrs + b3 (iflag) * lrs32 + b4 ( &
       iflag) * lrs2)
  dom = 2.d0 * a * (0.5d0 * b1 * rs12 + b2 * rs + 1.5d0 * b3 ( &
       iflag) * rs32 + 2.d0 * b4 (iflag) * rs2)
  dlom = 2.d0 * a * (0.5d0*b1*rs/lrs12+b2*rs+1.5d0*b3(iflag)*lrs12*rs &
       +2.d0*b4(iflag)*lrs*rs)
  olog = log (1.d0 + 1.0d0 / om)
  lolog = log (1.d0+1.0d0/lom)
  ec = - 2.d0 * a * (1.d0 + a1 (iflag) * rs) * olog
  lec = - 2.d0 * a * (1.d0 + a1 (iflag) * lrs) * lolog
  dlec = -2.d0*a*a1(iflag)*rs*lolog  &
         +2.d0*a*(1.d0+a1(iflag)*lrs)*dlom/(lom*(lom+1.d0))
  ec_l=2.d0*cc*lec+cc*cc*dlec
  vc = - 2.d0 * a * (1.d0 + 2.d0 / 3.d0 * a1 (iflag) * rs) &
       * olog - 2.d0 / 3.d0 * a * (1.d0 + a1 (iflag) * rs) * dom / &
       (om * (om + 1.d0) )
  return
end subroutine pwcc
!
!-----------------------------------------------------------------------
subroutine pwcc_spin (rs,cc, zeta, ec, vcup, vcdw, ec_l)
  !-----------------------------------------------------------------------
  !     This subroutine is adapted from subroutine pw_spin from 
  !     Modules/lsda_functionals.f90
  !
  !     J.P. Perdew and Y. Wang, PRB 45, 13244 (1992)
  !
  USE kinds, ONLY : DP
  implicit none
  real(DP) :: rs,cc, zeta, ec, vcup, vcdw,ec_l
  ! xc parameters, unpolarised
  real(DP) :: a, a1, b1, b2, b3, b4, c0, c1, c2, c3, d0, d1
  parameter (a = 0.031091d0, a1 = 0.21370d0, b1 = 7.5957d0, b2 = &
       3.5876d0, b3 = 1.6382d0, b4 = 0.49294d0, c0 = a, c1 = 0.046644d0, &
       c2 = 0.00664d0, c3 = 0.01043d0, d0 = 0.4335d0, d1 = 1.4408d0)
  ! xc parameters, polarised
  real(DP) :: ap, a1p, b1p, b2p, b3p, b4p, c0p, c1p, c2p, c3p, d0p, &
       d1p
  parameter (ap = 0.015545d0, a1p = 0.20548d0, b1p = 14.1189d0, b2p &
       = 6.1977d0, b3p = 3.3662d0, b4p = 0.62517d0, c0p = ap, c1p = &
       0.025599d0, c2p = 0.00319d0, c3p = 0.00384d0, d0p = 0.3287d0, d1p &
       = 1.7697d0)
  ! xc parameters, antiferro
  real(DP) :: aa, a1a, b1a, b2a, b3a, b4a, c0a, c1a, c2a, c3a, d0a, &
       d1a
  parameter (aa = 0.016887d0, a1a = 0.11125d0, b1a = 10.357d0, b2a = &
       3.6231d0, b3a = 0.88026d0, b4a = 0.49671d0, c0a = aa, c1a = &
       0.035475d0, c2a = 0.00188d0, c3a = 0.00521d0, d0a = 0.2240d0, d1a &
       = 0.3969d0)
  real(DP) :: fz0
  parameter (fz0 = 1.709921d0)
  real(DP) :: rs12, rs32, rs2, zeta2, zeta3, zeta4, fz, dfz
  real(DP) :: om, dom, olog, epwc, vpwc
  real(DP) :: omp, domp, ologp, epwcp, vpwcp
  real(DP) :: oma, doma, ologa, alpha, vpwca
  real(DP) :: lrs,lrs12,lrs32,lrs2
  real(DP) :: lom,dlom,lolog,lepwc,dlepwc
  real(DP) :: lomp,dlomp,lologp,lepwcp,dlepwcp
  real(DP) :: loma,dloma,lologa,lalpha,dlalpha
  real(DP) :: lec,dlec
  !
  !     if(rs.lt.0.5d0) then
  ! high density formula (not implemented)
  !
  !     else if(rs.gt.100.d0) then
  ! low density formula  (not implemented)
  !
  !     else
  ! interpolation formula
  zeta2 = zeta * zeta
  zeta3 = zeta2 * zeta
  zeta4 = zeta3 * zeta
  rs12 = sqrt (rs)
  rs32 = rs * rs12
  rs2 = rs**2
  lrs=cc*rs
  lrs12 = sqrt (lrs)
  lrs32 = lrs*lrs12
  lrs2=lrs**2
  ! unpolarised
  om = 2.d0 * a * (b1 * rs12 + b2 * rs + b3 * rs32 + b4 * rs2)
  dom = 2.d0 * a * (0.5d0 * b1 * rs12 + b2 * rs + 1.5d0 * b3 * rs32 &
       + 2.d0 * b4 * rs2)
  olog = log (1.d0 + 1.0d0 / om)
  epwc = - 2.d0 * a * (1.d0 + a1 * rs) * olog
  vpwc = - 2.d0 * a * (1.d0 + 2.d0 / 3.d0 * a1 * rs) * olog - 2.d0 / &
       3.d0 * a * (1.d0 + a1 * rs) * dom / (om * (om + 1.d0) )
  ! polarized
  omp = 2.d0 * ap * (b1p * rs12 + b2p * rs + b3p * rs32 + b4p * rs2)
  domp = 2.d0 * ap * (0.5d0 * b1p * rs12 + b2p * rs + 1.5d0 * b3p * &
       rs32 + 2.d0 * b4p * rs2)
  ologp = log (1.d0 + 1.0d0 / omp)
  epwcp = - 2.d0 * ap * (1.d0 + a1p * rs) * ologp
  vpwcp = - 2.d0 * ap * (1.d0 + 2.d0 / 3.d0 * a1p * rs) * ologp - &
       2.d0 / 3.d0 * ap * (1.d0 + a1p * rs) * domp / (omp * (omp + 1.d0) &
       )
  ! antiferro
  oma = 2.d0 * aa * (b1a * rs12 + b2a * rs + b3a * rs32 + b4a * rs2)
  doma = 2.d0 * aa * (0.5d0 * b1a * rs12 + b2a * rs + 1.5d0 * b3a * &
       rs32 + 2.d0 * b4a * rs2)
  ologa = log (1.d0 + 1.0d0 / oma)
  alpha = 2.d0 * aa * (1.d0 + a1a * rs) * ologa
  vpwca = + 2.d0 * aa * (1.d0 + 2.d0 / 3.d0 * a1a * rs) * ologa + &
       2.d0 / 3.d0 * aa * (1.d0 + a1a * rs) * doma / (oma * (oma + 1.d0) &
       )
  ! coupling constant depandent
  lom = 2.d0*a*(b1*lrs12+b2*lrs+b3*lrs32+b4*lrs2)
  dlom = 2.d0 * a * (0.5d0*b1*rs/lrs12+b2*rs+1.5d0*b3*lrs12*rs &
       +2.d0*b4*lrs*rs)
  lolog = log (1.d0+1.0d0/lom)
  lepwc=-2.d0*a*(1.d0+a1*lrs)*lolog 
  dlepwc=-2.d0*a*a1*rs*lolog &
         +2.d0*a*(1.d0+a1*lrs)*dlom/(lom*(lom+1.d0))

  lomp=2.d0*ap*(b1p*lrs12+b2p*lrs+b3p*lrs32+b4p*lrs2)
  dlomp=2.d0*ap*(0.5d0*b1p*rs/lrs12+b2p*rs+1.5d0*b3p*lrs12*rs &
       +2.d0*b4p*lrs*rs)
  lologp=log(1.d0+1.d0/lomp)
  lepwcp=-2.d0*ap*(1.d0+a1p*lrs)*lologp
  dlepwcp=-2.d0*ap*a1p*rs*lologp+2.d0*ap*(1.d0+a1p*lrs)*dlomp/(lomp*(lomp+1.d0))

  loma=2.d0*aa*(b1a*lrs12+b2a*lrs+b3a*lrs32+b4a*lrs2)
  dloma=2.d0*aa*(0.5d0*b1a*rs/lrs12+b2a*rs+1.5d0*b3a*lrs12*rs &
       +2.d0*b4a*lrs*rs)
  lologa=log(1.d0+1.d0/loma)
  lalpha=2.d0*aa*(1.d0+a1a*lrs)*lologa
  dlalpha=2.d0*aa*a1a*rs*lologa  &
          -2.d0*aa*(1.d0+a1a*lrs)*dloma/(loma*(loma+1.d0)) 

  !
  fz = ( (1.d0 + zeta) ** (4.d0 / 3.d0) + (1.d0 - zeta) ** (4.d0 / &
       3.d0) - 2.d0) / (2.d0** (4.d0 / 3.d0) - 2.d0)
  dfz = ( (1.d0 + zeta) ** (1.d0 / 3.d0) - (1.d0 - zeta) ** (1.d0 / &
       3.d0) ) * 4.d0 / (3.d0 * (2.d0** (4.d0 / 3.d0) - 2.d0) )
  !
  ec = epwc + alpha * fz * (1.d0 - zeta4) / fz0 + (epwcp - epwc) &
       * fz * zeta4
  lec = lepwc+lalpha*fz*(1.d0-zeta4)/fz0+(lepwcp-lepwc)*fz*zeta4
  dlec=dlepwc+dlalpha*fz*(1.d0-zeta4)/fz0+(dlepwcp-dlepwc)*fz*zeta4
  ec_l=2.d0*cc*lec+cc*cc*dlec
  !
  vcup = vpwc + vpwca * fz * (1.d0 - zeta4) / fz0 + (vpwcp - vpwc) &
       * fz * zeta4 + (alpha / fz0 * (dfz * (1.d0 - zeta4) - 4.d0 * fz * &
       zeta3) + (epwcp - epwc) * (dfz * zeta4 + 4.d0 * fz * zeta3) ) &
       * (1.d0 - zeta)

  vcdw = vpwc + vpwca * fz * (1.d0 - zeta4) / fz0 + (vpwcp - vpwc) &
       * fz * zeta4 - (alpha / fz0 * (dfz * (1.d0 - zeta4) - 4.d0 * fz * &
       zeta3) + (epwcp - epwc) * (dfz * zeta4 + 4.d0 * fz * zeta3) ) &
       * (1.d0 + zeta)
  !      endif
  !
  return
end subroutine pwcc_spin
!
!----------------------------------------
! ####################################################################
!                          |              |
!                          |  PPACF_INFO  |
!                          |____________ _|

SUBROUTINE ppacf_info

  USE io_global,            ONLY : stdout
  implicit none




  WRITE(stdout,'(/)')
  WRITE(stdout,'(5x,"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")')
  WRITE(stdout,'(5x,"%                                                                      %")')
  WRITE(stdout,'(5x,"% You are using PPACF, please cite the following paper:                %")') 
  WRITE(stdout,'(5x,"%                                                                      %")')
  WRITE(stdout,'(5x,"%   Y. Jiao, E. Schr\""oder, and P. Hyldgaard, PRB 97, 085115 (2018).   %")')
  WRITE(stdout,'(5x,"%                                                                      %")')
  WRITE(stdout,'(5x,"% If you are using this code for hybrid mixing value, please also cite:%")')
  WRITE(stdout,'(5x,"%                                                                      %")')
  WRITE(stdout,'(5x,"%   Y. Jiao, E. Schr\""oder, and P. Hyldgaard, JCP 148, 194115 (2018).  %")')
  WRITE(stdout,'(5x,"%                                                                      %")')
  WRITE(stdout,'(5x,"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")')
  WRITE(stdout,'(/)')


END SUBROUTINE

