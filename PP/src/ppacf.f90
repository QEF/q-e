!-----------------------------------------------------------------------
!
! Program written by Yang Jiao,  Oct 2016, GPL, No warranties.
!
!    Oct 2018, adapted to QE6.3
!    Jan 2019, adapted to change in rho from (up,down) to (tot,magn)
!    June 2020, adapted to extended vdW_DF analysis and refined output (PH)
!
!-----------------------------------------------------------------------
PROGRAM do_ppacf
  !-----------------------------------------------------------------------
  !! This routine computes the coupling constant dependency of exchange
  !! correlation potential \( E_{\text{xc},\lambda}, \lambda \in \[0:1\]
  !! and the spatial distribution of exchange correlation energy 
  !! density and kinetic correlation energy density according to:
  !
  !! Y. Jiao, E. Schr\"oder, and P. Hyldgaard, Phys. Rev. B 97, 085115 (2018).
  !
  !! For an illustration of how to use this routine to set hybrid 
  !! mixing parameter, please refer to:
  !     
  !! Y. Jiao, E. Schr\"oder, and P. Hyldgaard, J. Chem. Phys. 148, 194115 (2018).     
  !
  !! Finally, this routine can also be used to set isolate the 
  !  Ashcroft-type pure-dispersion component of E_{c;vdw}^nl 
  !  (or the cumulant reminder, E_{c;alpha}^nl, defining a local-field susceptibility):
  !     
  !! P. Hyldgaard, Y. Jiao, and V. Shukla, J. Phys.: Condens. Matt. 32, 393001 (2020):     
  !! https://iopscience.iop.org/article/10.1088/1361-648X/ab8250
  !
  USE basis,                ONLY : starting_wfc
  USE constants,            ONLY : e2, pi, fpi
  USE kinds,                ONLY : DP
  USE klist,                ONLY : nks, xk, ngk, igk_k
  USE gvect,                ONLY : ngm, g
  USE gvecw,                ONLY : ecutwfc, gcutw
  USE io_files,             ONLY : prefix, tmp_dir
  USE io_global,            ONLY : stdout, ionode, ionode_id
  USE cell_base,            ONLY : omega
  USE mp,                   ONLY : mp_bcast, mp_sum
  USE mp_world,             ONLY : world_comm
  USE mp_global,            ONLY : mp_startup
  USE mp_bands,             ONLY : intra_bgrp_comm
  USE exx,                  ONLY : exxinit, exxenergy2, fock2, ecutfock, & 
                                   use_ace, aceinit, local_thr, nbndproj
  USE exx_base,             ONLY : exx_grid_init, exx_mp_init, exx_div_check, &
                                   exxdiv_treatment
  USE exx_base,             ONLY : nq1, nq2, nq3
  USE fft_base,             ONLY : dfftp
  USE fft_interfaces,       ONLY : invfft
  USE lsda_mod,             ONLY : nspin
  USE scf,                  ONLY : scf_type, create_scf_type, destroy_scf_type
  USE scf,                  ONLY : scf_type_COPY
  USE scf,                  ONLY : rho, rho_core, rhog_core, vltot
  USE funct,                ONLY : dft_is_nonlocc, nlc, enforce_input_dft   
  USE xc_lib,               ONLY : xclib_get_id, xclib_set_auxiliary_flags,    &
                                   xclib_set_exx_fraction, xclib_dft_is_libxc, &
                                   xclib_set_threshold, xc, xc_gcx
  USE wvfct,                ONLY : npw, npwx
  USE environment,          ONLY : environment_start, environment_end
  USE vdW_DF,               ONLY : Nqs, vdW_DF_potential, vdW_DF_energy, vdW_DF_analysis
  USE vdW_DF_scale,         ONLY : xc_vdW_DF_ncc, xc_vdW_DF_spin_ncc, &
                                   get_q0cc_on_grid, get_q0cc_on_grid_spin
  USE vasp_xml,             ONLY : readxmlfile_vasp
  !
  IMPLICIT NONE
  !
  INTEGER :: code_num
  ! From which code to read in the calculation data:  
  ! 1 \(\rightarrow\) Quantum ESPRESSO (default);  
  ! 2 \(\rigtharrow\) VASP.
  INTEGER :: vdW_analysis
  ! Analysis of which component of vdW-DF nonlocal-correlation term?
  ! 0 \(\rightarrow\) Use full-method kernel (default, sum of 1+2 kernels);  
  ! 1 \(\rigtharrow\) Susceptibility-type vdW_DF kernel (postprocessing ONLY!)
  ! 2 \(\rigtharrow\) Pure-vdW (longitudinal) vdW_DF kernel (postprocessing ONLY!)
  ! Kernels defined/explained in IOP JCPM focused review.
  LOGICAL :: is_libxc(2), isnonlocc
  LOGICAL :: lplot, ltks, lfock
  ! lplot: "2 track spatial variation or not 2 track spatial variation?"
  ! if (lplot) ltks: "2 track Kohn-sham kinetic energy?" (requires wavefunctions)
  ! if (lplot) lfock: "Compute Fock exchange value?" (requires wavefunctions)
  LOGICAL :: lecnl_qxln, lecnl_qx
  ! lecnl_qxln: "Linearized coupling-constant variation in Ecnl analysis" (Sanity check)
  ! lecnl_qx: "TESTING ONLY! Ignore qc (but not qx) in setting q_0 in
  !            Ecnl coupling-constant analysis" (Test! STRONGLY DEPRECIATED).
  INTEGER :: icc, ncc
  INTEGER :: n_lambda
  REAL(DP):: rs, rs3, s, q, qx, qc
  REAL(DP):: etxclambda, etx, etxlda, etxgc, etcnlccc, etcnlcccp, etcnlcccm
  REAL(DP) :: etcldalambda, etcgclambda, etcnlclambda, etcnl_check, &
              ttcnl_check, tcnl_int
  REAL(DP), ALLOCATABLE :: Ec_nl_ngamma(:)
  REAL(DP) :: etc, etclda, etcgc,etcnl,etcnlncc
  REAL(DP), EXTERNAL :: exxenergyace
  !
  INTEGER :: is, ir, iq, ig, icar, nnrtot
  ! counter on mesh points
  ! counter on nspin
  INTEGER :: iexch, icorr, igcx, igcc, inlc
  INTEGER :: ierr, ios
  LOGICAL :: needwf = .FALSE.
  REAL(DP) :: cc, dcc, ccp, ccm, ccp2, ccm2, ccp3, ccm3, ccp4, ccm4, ccp8, ccm8, cc3
  ! coupling constant
  ! local exchange energy, local correlation energy
  ! local exchange potential, local correlation potential
  REAL(DP) :: rhox, arhox, zeta, grho2
  ! the charge in each point
  ! the absolute value of the charge
  !
  ! XC input/output arrays
  REAL(DP), ALLOCATABLE :: v1x(:,:), v1c(:,:), v2x(:,:), v2c(:,:)
  REAL(DP), ALLOCATABLE :: ex(:), ec(:), sx(:), sc(:)
  REAL(DP), ALLOCATABLE :: expp(:), ecpp(:), exm(:), ecm(:)
  REAL(DP), ALLOCATABLE :: sxd(:), scp(:), scm(:)
  !
  REAL(DP) :: etxc, vtxc
  REAL(DP) :: ec_l, ecgc_l, Ec_nl, ec_cc
  REAL(DP) :: etxccc, etxcccnl, etxcccnlp, etxcccnlm, vtxccc, vtxccc_buf, vtxcccnl 
  !
  REAL(DP) :: etxcgc, vtxcgc, segno, rh
  !
  REAL(DP) :: dq0_dq  ! The derivative of the saturated
  REAL(DP) :: grid_cell_volume
  !
  REAL(DP), ALLOCATABLE :: q0(:)
  REAL(DP), ALLOCATABLE :: vofrcc(:,:)
  !
  REAL(DP), PARAMETER :: vanishing_charge = 1.D-10, &
                         vanishing_mag = 1.D-20
  REAL(DP), PARAMETER :: small = 1.E-10_DP,  third = 1.0_DP / 3.0_DP, &
                         pi34 = 0.6203504908994_DP  ! pi34=(3/4pi)^(1/3)
  REAL(DP), PARAMETER :: epsr = 1.D-6, epsg = 1.D-10
  REAL(DP), ALLOCATABLE :: grho(:,:,:), rhoout(:,:), rhoout_tz(:,:)
  COMPLEX(DP), ALLOCATABLE :: rhogsum(:,:)
  REAL(DP), ALLOCATABLE :: tot_grad_rho(:,:), grad_rho(:,:,:)
  REAL(DP), ALLOCATABLE :: tot_rho(:)
  !
  INTEGER :: ik ! counter on k points
  INTEGER, ALLOCATABLE  :: igk_buf(:)
  REAL(DP), ALLOCATABLE :: gk(:) ! work space
  !
  CHARACTER(LEN=256) :: outdir
  !
  CHARACTER(LEN=256) :: filplot
  INTEGER  :: plot_num
  !
  REAL(DP) :: ttclda        ! the LDA kinetic-correlation energy T_c^LDA=E_c,lambda=1^LDA-E_c^LDA
  REAL(DP) :: ttcgc         ! the GGA-LDA kinetic-correlation energy T_c^gc=E_c,lambda=1^gc-E_c^gc
  REAL(DP) :: gtimesr
  REAL(DP) :: r_shift(3)
  REAL(DP), ALLOCATABLE :: dq0_drho(:,:)
  REAL(DP), ALLOCATABLE :: dq0_drho_up(:)       ! The derivative of the saturated q0
  REAL(DP), ALLOCATABLE :: dq0_drho_down(:)     ! (equation 5 of SOLER) with respect
                                                ! to the charge density (see
                                                ! get_q0_on_grid subroutine for details).
  !
  REAL(DP), ALLOCATABLE :: dq0_dgradrho(:,:) 
  REAL(DP), ALLOCATABLE :: dq0_dgradrho_up(:)   ! The derivative of the saturated q0
  REAL(DP), ALLOCATABLE :: dq0_dgradrho_down(:) ! (equation 5 of SOLER) with respect
                                                ! to the gradient of the charge density
                                                ! (again, see get_q0_on_grid subroutine).
  REAL(DP), ALLOCATABLE :: potential_vdW(:,:)   ! The vdW contribution to the potential
  COMPLEX(DP), ALLOCATABLE :: u_vdW(:,:)        ! 
  COMPLEX(DP), ALLOCATABLE :: up_vdW(:,:), um_vdW(:,:)
  COMPLEX(DP), ALLOCATABLE :: thetas(:,:)       ! These are the functions of equation 8 of
                                                ! SOLER. They will be forward Fourier transformed
                                                ! in place to get theta(k) and worked on in
                                                ! place to get the u_alpha(r) of equation 11
                                                ! in SOLER. They are formatted as follows:
                                                ! thetas(grid_point, theta_i).
  COMPLEX(DP), ALLOCATABLE :: thetasp(:,:), thetasm(:,:)
  COMPLEX(DP), ALLOCATABLE :: ecnl_c(:), tcnl_c(:)

  REAL(DP), ALLOCATABLE    :: kin_r(:,:) ! the kinetic energy density in R-space
  !
  TYPE(scf_type) :: exlda, eclda
  TYPE(scf_type) :: tclda  ! the LDA kinetic-correlation energy per particle
  TYPE(scf_type) :: ecnl   ! the nonlocal correlation energy per particle
  TYPE(scf_type) :: tcnl   ! the Ecnl kinetic-correlation energy per particle
  TYPE(scf_type) :: exgc, ecgc, tcgc
  !
  NAMELIST / ppacf / code_num,outdir, prefix, n_lambda, vdW_analysis, lplot, ltks, lfock, use_ace
  !
  ! initialise environment
  !
#ifdef __MPI
  CALL mp_startup()
#endif
  !--------------- READ IN PREFIX --------------------------------!
  CALL environment_start( 'ppacf' )
  !
  is_libxc(1) = xclib_dft_is_libxc('GGA','EXCH')
  is_libxc(2) = xclib_dft_is_libxc('GGA','CORR')
  !
  IF ( ANY(.NOT.is_libxc(1:2)) ) CALL xclib_set_threshold( 'gga', 1.E-10_DP, 1.E-10_DP )
  !
  ! ... set default values for variables in namelist
  !
  code_num = 1
  outdir = './'
  prefix = 'ppacf'
  n_lambda = 1
  vdW_analysis = 0 
  lplot = .FALSE.
  ltks  = .FALSE.
  lfock = .FALSE.
  use_ace = .TRUE.
  dcc = 1.E-6_DP 
  filplot = 'tmp.pp'
  plot_num = -1
  lecnl_qxln = .FALSE.
  lecnl_qx = .FALSE.
  !
  IF (ionode) THEN
     !
     CALL input_from_file()
     !
     READ(5, ppacf, iostat = ios) 
     !
  ENDIF
  !
  IF (nspin == 4) CALL errore( 'ppacf', 'noncollinear spin not implemented', 1 )
  !
  !--------------- READ IN DATA ----------------------------------------------!
  ! 
  ! ... Broadcast variables
  CALL mp_bcast( code_num,   ionode_id, world_comm )
  CALL mp_bcast( outdir,     ionode_id, world_comm )
  CALL mp_bcast( prefix,     ionode_id, world_comm )
  CALL mp_bcast( n_lambda,   ionode_id, world_comm ) 
  CALL mp_bcast( lplot,      ionode_id, world_comm )
  CALL mp_bcast( ltks,       ionode_id, world_comm )
  CALL mp_bcast( lfock,      ionode_id, world_comm )
  CALL mp_bcast( vdW_analysis, ionode_id, world_comm )
  CALL mp_bcast( lecnl_qxln, ionode_id, world_comm )
  CALL mp_bcast( lecnl_qx,   ionode_id, world_comm )
  CALL mp_bcast( dcc,        ionode_id, world_comm )
  !
  SELECT CASE ( vdW_analysis )
  CASE ( 0 ) 
     vdW_DF_analysis = vdW_analysis
  CASE ( 1,2 )
     vdW_DF_analysis = vdW_analysis
     IF (ionode) THEN 
       WRITE(stdout,'(/)') 
       WRITE(stdout,'(5X,a)') "NOTE: Coupling-constant analysis of correlation for" 
       IF (vdW_analysis == 1) THEN
         WRITE(stdout,'(5X,a)') "Ecnl_Alpha (susceptibility kernel), not for Ecnl."
       ELSE
         WRITE(stdout,'(5X,a)') "Ecnl_vdW (pure-vdW kernel), not for full Ecnl."
       END IF
     END IF
  CASE DEFAULT
     CALL errore( 'ppacf', 'vdW_DF analysis kernel not implemented', 1 )
  END SELECT
  !
  ncc = n_lambda
  WRITE( stdout, '(//5x,"entering subroutine acf ..."/)')
  !
  ! ... Write out the ppacf information.
  IF (ionode) CALL ppacf_info()
  !
  CALL start_clock( 'acf_etxclambda' )
  !
  ! WRITE(stdout,9093) dcc
  !
  IF (code_num == 1) THEN
     !
     tmp_dir = TRIM(outdir) 
     IF ( lfock .OR. (lplot .AND. ltks) ) THEN
        CALL read_file ( )
     ELSE
        CALL read_file_new ( needwf )
     ENDIF
     !
     ! Check exchange correlation functional
     !
     iexch = xclib_get_id('LDA','EXCH')
     icorr = xclib_get_id('LDA','CORR')
     igcx  = xclib_get_id('GGA','EXCH')
     igcc  = xclib_get_id('GGA','CORR')
     !
  ELSEIF (code_num == 2) THEN
     !
     tmp_dir = TRIM(outdir)
     CALL readxmlfile_vasp( iexch, icorr, igcx, igcc, inlc, ierr ) ! Presently has generate_kernel() 
                             ! if nonlocal dft --- So this call should be after vdW_DF_analysis is set
     IF (ionode) WRITE(stdout,'(5X,a)') "Read data from VASP output 'vasprun.xml'"
     ! 
  ELSE
     CALL errore( 'ppacf', 'code_num not implemented', 1 )
  ENDIF
  !
  ALLOCATE( vofrcc(1:dfftp%nnr,1:nspin) )
  !
  IF(dft_is_nonlocc()) THEN
     WRITE( stdout, '(//5x,"ACF coupling-constant  Exc_lambda (Ry)  &
                    &E_c,lambda^LDA (Ry)  E_c,lambda^nl (Ry)"/)')
  ELSE
     WRITE( stdout, '(//5x,"ACF coupling-constant  Exc_lambda (Ry)  &
                    &E_c,lambda^LDA (Ry)  E_c,lambda^GC (Ry)"/)')
  ENDIF
  !
  IF (dft_is_nonlocc() .AND. lecnl_qxln) THEN
     etxcccnl = 0._DP
     vtxcccnl = 0._DP
     vofrcc = 0._DP
     !
     CALL nlc( rho%of_r, rho_core, nspin, etxcccnl, vtxcccnl, vofrcc )
     !
     CALL mp_sum(  etxcccnl, intra_bgrp_comm )
  ENDIF
  !
  ! ... add gradient corrections (if any)
  !
  ALLOCATE( grho(3,dfftp%nnr,nspin) )
  ALLOCATE( rhoout(dfftp%nnr,nspin) )
  ALLOCATE( rhoout_tz(dfftp%nnr,nspin) ) 
  ALLOCATE( tot_rho(dfftp%nnr) )
  ALLOCATE( rhogsum(ngm,nspin) )
  !
  ! ... calculate the gradient of rho + rho_core in real space
  ! ... note: input rho is (tot,magn), output rhoout, grho are (up,down)
  !
  IF ( nspin == 1 ) THEN
     !
     rhoout (:,1) = rho_core(:)  + rho%of_r(:,1) 
     rhogsum(:,1) = rhog_core(:) + rho%of_g(:,1)
     !
  ELSE IF ( nspin == 2 ) THEN
     !
     rhoout (:,1) = (rho_core(:)  + rho%of_r(:,1) + rho%of_r(:,2) )/2.0_dp
     rhoout (:,2) = (rho_core(:)  + rho%of_r(:,1) - rho%of_r(:,2) )/2.0_dp
     rhogsum(:,1) = (rhog_core(:) + rho%of_g(:,1) + rho%of_g(:,2) )/2.0_dp
     rhogsum(:,2) = (rhog_core(:) + rho%of_g(:,1) - rho%of_g(:,2) )/2.0_dp
     !
  END IF
  !
  DO is = 1, nspin
     CALL fft_gradient_g2r( dfftp, rhogsum(1,is), g, grho(1,1,is))
  END DO
  !
  DEALLOCATE( rhogsum )
  !
  tot_rho (:) = rho_core(:)  + rho%of_r(:,1) 
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
     CALL create_scf_type( ecnl )
     CALL create_scf_type( tcnl )
     ecnl%of_r(:,:) = 0._DP
     tcnl%of_r(:,:) = 0._DP
     ALLOCATE( Ec_nl_ngamma(1:ncc) )
     Ec_nl_ngamma = 0._DP
  ELSEIF(igcc /= 0) THEN
     CALL create_scf_type( ecgc )
     CALL create_scf_type( tcgc )
     ecgc%of_r(:,:) = 0._DP
     tcgc%of_r(:,:) = 0._DP
  ENDIF
  ttclda = 0._DP
  ttcgc = 0._DP
  etcnl = 0._DP
  etcnlncc = 0._DP
  !
  !
  ALLOCATE( v1x(dfftp%nnr,nspin), v1c(dfftp%nnr,nspin) )
  ALLOCATE( v2x(dfftp%nnr,nspin), v2c(dfftp%nnr,nspin) )
  ALLOCATE( ex(dfftp%nnr), ec(dfftp%nnr) )
  ALLOCATE( sx(dfftp%nnr), sc(dfftp%nnr) )
  ALLOCATE( expp(dfftp%nnr), ecpp(dfftp%nnr) )
  ALLOCATE( exm(dfftp%nnr), ecm(dfftp%nnr) )
  ALLOCATE( sxd(dfftp%nnr), scp(dfftp%nnr), scm(dfftp%nnr) )
  !
  rhoout_tz(:,1) = rho%of_r(:,1) + rho_core(:)
  IF (nspin==2) rhoout_tz(:,2) = rho%of_r(:,2)
  !
  CALL xc( dfftp%nnr, nspin, nspin, rhoout_tz, ex, ec, v1x, v1c )
  !
  CALL xc_gcx( dfftp%nnr, nspin, rhoout, grho, sx, sc, v1x, v2x, v1c, v2c )
  !
  !
  ! ... coupling constant > 0
  !
  DO icc = 0, ncc
     cc = DBLE(icc)/DBLE(ncc)
     etxclambda = 0._DP
     etxccc = 0._DP
     vofrcc = 0._DP
     etx = 0._DP
     etxlda = 0._DP
     etxgc  = 0._DP
     etcldalambda = 0._DP
     etcgclambda  = 0._DP
     etc = 0._DP
     etclda = 0._DP
     etcgc  = 0._DP
     !
     IF (cc > 0) THEN
       !
       ccp = cc+dcc
       ccm = cc-dcc
       ccp2 = ccp*ccp
       ccp3 = ccp2*ccp
       ccp4 = ccp3*ccp
       ccp8 = ccp4*ccp4
       ccm2 = ccm*ccm
       ccm3 = ccm2*ccm
       ccm4 = ccm3*ccm
       ccm8 = ccm4*ccm4
       !
       IF (icorr /= 4 .OR. is_libxc(2)) THEN
         CALL xc( dfftp%nnr, nspin, nspin, rhoout_tz/ccp3, expp, ecpp, v1x, v1c )
         CALL xc( dfftp%nnr, nspin, nspin, rhoout_tz/ccm3, exm,  ecm,  v1x, v1c )
       ENDIF
       !
       IF (igcc /= 0) THEN
         CALL xc_gcx( dfftp%nnr, nspin, rhoout/ccp3, grho/ccp4, sxd, scp, v1x,&
                      v2x, v1c, v2c )
         !
         CALL xc_gcx( dfftp%nnr, nspin, rhoout/ccm3, grho/ccm4, sxd, scm, v1x,&
                      v2x, v1c, v2c )
       ENDIF
       !   
     ENDIF
     !
     !
     IF (nspin == 1 ) THEN
       !
       ! ... spin-unpolarized case
       !
       DO ir = 1, dfftp%nnr
          !
          rhox = rhoout_tz(ir,1)
          arhox = ABS(rhox)
          IF (arhox > vanishing_charge) THEN
             !
             etx = etx + e2*ex(ir)*rhox
             etxlda = etxlda + e2*ex(ir)*rhox
             grho2 = grho(1,ir,1)**2 + grho(2,ir,1)**2 + grho(3,ir,1)**2
             !
             IF (cc > 0._DP) THEN
                !
                IF (icorr == 4 .AND. .NOT.is_libxc(2) ) THEN
                   rs = pi34 /arhox**third
                   CALL pwcc( rs, cc, ec_cc, v1c(1,1), ec_l )
                ELSE
                   ec_l = (ccp2*ecpp(ir)-ccm2*ecm(ir))/dcc*0.5_DP
                ENDIF
                !
                etcldalambda = etcldalambda + e2*ec_l*rhox
                !
                IF(icc == ncc) THEN
                   IF (icorr /= 4 .OR. is_libxc(2)) THEN
                     tclda%of_r(ir,1) = e2*(ec(ir)-ec_l)*rhox
                     ttclda = ttclda + e2*(ec(ir)-ec_l)*rhox
                   ELSE
                     tclda%of_r(ir,1) = e2*(ec_cc-ec_l)*rhox
                     ttclda = ttclda + e2*(ec_cc-ec_l)*rhox
                   ENDIF
                ENDIF
                !
                IF (grho2>epsg .AND. igcc/=0) THEN
                   segno = SIGN( 1.D0, rhoout(ir,1) )
                   ecgc_l = (ccp2*scp(ir)*ccp3-ccm2*scm(ir)*ccm3)/dcc*0.5_DP
                   etcgclambda = etcgclambda + e2*ecgc_l*segno
                ENDIF
             ENDIF
             !
             etclda = etclda + e2*ec(ir)*rhox
             etc = etc + e2*ec(ir)*rhox
             !
             IF (icc == ncc) THEN
                exlda%of_r(ir,1) = e2*ex(ir)*rhox
                eclda%of_r(ir,1) = e2*ec(ir)*rhox
             ENDIF
             !
             IF ( grho2 > epsg ) THEN
                segno = SIGN( 1.D0, rhoout(ir,1) )
                !
                etx = etx + e2*sx(ir)*segno
                etxgc = etxgc + e2*sx(ir)*segno
                etc = etc + e2*sc(ir)*segno
                etcgc = etcgc + e2*sc(ir)*segno
                IF (icc == ncc) THEN
                   exgc%of_r(ir,1) = e2*sx(ir)*segno
                   IF (igcc /= 0) THEN
                      ecgc%of_r(ir,1) = e2*sc(ir)*segno
                      tcgc%of_r(ir,1) = e2*(sc(ir)-ecgc_l)*segno 
                      ttcgc=ttcgc+e2*(sc(ir)-ecgc_l)*segno
                   ENDIF
                ENDIF
             ENDIF
             !
          ENDIF
          !
       ENDDO
       !
     ELSEIF (nspin == 2) THEN
       !
       ! ... spin-polarized case
       !
       DO ir = 1, dfftp%nnr
          rhox = rhoout_tz(ir,1)
          arhox = ABS( rhox )
          IF (arhox > vanishing_charge) THEN
             rs = pi34 /arhox**third
             zeta = rhoout_tz(ir,2)/arhox
             IF (ABS(zeta) > 1.D0) zeta = SIGN(1.D0, zeta)
             etx = etx + e2*ex(ir)*rhox
             etxlda = etxlda+e2*ex(ir)*rhox
             !
             IF (cc > 0._DP) THEN
                !
                IF (icorr == 4 .AND. .NOT.is_libxc(2)) THEN
                   CALL pwcc_spin( rs, cc, zeta, ec_cc, v1c(1,1), v1c(1,2), ec_l )
                ELSE
                   ec_l = (ccp2*ecpp(ir)-ccm2*ecm(ir))/dcc*0.5_DP
                ENDIF
                !
                etcldalambda = etcldalambda+e2*ec_l*rhox
                IF (icc == ncc) THEN
                   IF (icorr /= 4 .OR. is_libxc(2) ) THEN
                     tclda%of_r(ir,1) = e2*(ec(ir)-ec_l)*rhox
                     ttclda = ttclda + e2*(ec(ir)-ec_l)*rhox
                   ELSE
                     tclda%of_r(ir,1) = e2*(ec_cc-ec_l)*rhox
                     ttclda = ttclda + e2*(ec_cc-ec_l)*rhox
                   ENDIF
                ENDIF
                !
                IF (igcc /= 0) THEN
                   ecgc_l = (ccp2*scp(ir)*ccp3-ccm2*scm(ir)*ccm3)/dcc*0.5_DP
                   etcgclambda = etcgclambda+e2*ecgc_l
                ENDIF
             ENDIF
             !
             etclda = etclda + e2*ec(ir)*rhox
             etc = etc + e2*ec(ir)*rhox
             !
             IF (icc == ncc) THEN
                exlda%of_r(ir,1) = e2*ex(ir)*rhox
                eclda%of_r(ir,1) = e2*ec(ir)*rhox
             ENDIF
             !
             etx = etx + e2*sx(ir)
             etxgc = etxgc + e2*sx(ir)
             !
             etcgc = etcgc + e2*sc(ir)
             etc = etc + e2*sc(ir)
             IF ( icc == ncc ) THEN
                exgc%of_r(ir,1)=e2*sx(ir)
                IF(igcc /= 0) THEN
                   ecgc%of_r(ir,1)=e2*sc(ir)
                   tcgc%of_r(ir,1)=e2*(sc(ir)-ecgc_l)
                   ttcgc=ttcgc+e2*(sc(ir)-ecgc_l)
                ENDIF
             ENDIF
          ENDIF
       ENDDO
       !
     ELSEIF (nspin == 4) THEN
       !
       CALL errore( 'ppacf', 'Noncollinear not implemented', 1 )
       !
     ENDIF
     ! 
     grid_cell_volume = omega/DBLE(dfftp%nr1*dfftp%nr2*dfftp%nr3)
     etx = etx * grid_cell_volume
     etc = etc * grid_cell_volume
     etxlda = etxlda * grid_cell_volume
     etclda = etclda * grid_cell_volume
     etxgc = etxgc * grid_cell_volume
     etcgc = etcgc * grid_cell_volume
     etcldalambda = etcldalambda * grid_cell_volume
     etcgclambda = etcgclambda * grid_cell_volume
     ttclda = ttclda * grid_cell_volume
     ttcgc = ttcgc * grid_cell_volume
     !
     CALL mp_sum( etx, intra_bgrp_comm )
     CALL mp_sum( etc, intra_bgrp_comm )
     CALL mp_sum( etxlda, intra_bgrp_comm )
     CALL mp_sum( etclda, intra_bgrp_comm )
     CALL mp_sum( etxgc, intra_bgrp_comm )
     CALL mp_sum( etcgc, intra_bgrp_comm )
     CALL mp_sum( etcldalambda, intra_bgrp_comm )
     CALL mp_sum( etcgclambda, intra_bgrp_comm )
     CALL mp_sum( ttclda, intra_bgrp_comm )
     CALL mp_sum( ttcgc, intra_bgrp_comm )
     !
     ! Non-local correlation 
     etcnlclambda = 0._DP
     ! Approximation of exchange linear lambda-dependence
     IF (dft_is_nonlocc()) THEN
       IF (lecnl_qxln) THEN
         etcnlclambda = 2._DP*cc*etxcccnl
         ! Full lambda dependence
       ELSE  
         IF (cc > 0._DP) THEN
           ccp = cc+dcc
           ccm = cc-dcc
           IF (nspin == 1) THEN
             CALL xc_vdW_DF_ncc( cc,  lecnl_qx, etcnlccc  )
             CALL xc_vdW_DF_ncc( ccp, lecnl_qx, etcnlcccp ) 
             CALL xc_vdW_DF_ncc( ccm, lecnl_qx, etcnlcccm ) 
             etcnlclambda = 2._DP*cc*etcnlccc+0.5_DP*cc*cc*(etcnlcccp-etcnlcccm)/dcc
             Ec_nl_ngamma(icc) = etcnlccc
           ELSEIF(nspin == 2) THEN
             CALL xc_vdW_DF_spin_ncc( cc,  lecnl_qx, etcnlccc  )
             CALL xc_vdW_DF_spin_ncc( ccp, lecnl_qx, etcnlcccp ) 
             CALL xc_vdW_DF_spin_ncc( ccm, lecnl_qx, etcnlcccm )
             etcnlclambda = 2._DP*cc*etcnlccc+0.5_DP*cc*cc*(etcnlcccp-etcnlcccm)/dcc
             Ec_nl_ngamma(icc) = etcnlccc
           ENDIF
         ENDIF
       ENDIF
       !
       etc = etc + etcnlccc
       !
     ENDIF
     ! 
     etxclambda = etx+etcldalambda+etcgclambda+etcnlclambda
     IF (dft_is_nonlocc()) THEN
        WRITE(stdout,9091) cc, etxclambda, etcldalambda, etcnlclambda
        IF ( icc == ncc) etcnlncc = etcnlclambda
     ELSE
         WRITE(stdout,9092) cc, etxclambda, etcldalambda, etcgclambda
     ENDIF
     FLUSH(stdout)
     !
  ENDDO  ! icc
  !
  DEALLOCATE( rhoout_tz ) 
  DEALLOCATE( v1x, v1c, v2x, v2c )
  DEALLOCATE( ex, ec, sx, sc )
  DEALLOCATE( expp, ecpp, exm, ecm )
  DEALLOCATE( sxd, scp, scm )
  !
  IF ( dft_is_nonlocc() ) THEN
    WRITE(stdout,'(5x,a)') 'Ec_nl(n_1/lambda): '
    DO icc = 1, ncc
      cc = DBLE(icc)/DBLE(ncc)
      WRITE(stdout,9095) cc, Ec_nl_ngamma(icc)
      IF (icc == ncc) etcnl = Ec_nl_ngamma(icc)
    ENDDO
  ENDIF
  !
  IF (dft_is_nonlocc()) THEN
     IF (vdW_analysis==0) WRITE(stdout,'(a32,0PF17.8,a3)') 'Lieb-Oxford ratio, XC', (etx+etc)/etxlda !,etxc/etxlda
     IF (vdW_analysis==0) WRITE(stdout,'(a32,0PF17.8,a3)') 'Lieb-Oxford ratio, Exchange', etx/etxlda !,etx/etxlda
  ELSE
     WRITE(stdout,'(a32,0PF17.8,a3)') 'Lieb-Oxford ratio, XC', (etx+etc)/etxlda !,etxc/etxlda
     WRITE(stdout,'(a32,0PF17.8,a3)') 'Lieb-Oxford ratio, Exchange', etx/etxlda !,etx/etxlda
  END IF
  !
  WRITE(stdout,'(a32,0PF17.8,a3)') 'Exchange', etx, 'Ry'      !,etxlda,etxgc
  WRITE(stdout,'(a32,0PF17.8,a3)') 'LDA Exchange', etxlda, 'Ry'   !,etxlda
  WRITE(stdout,'(a32,0PF17.8,a3)') 'Correlation', etc, 'Ry'   !,etclda,etcgc
  WRITE(stdout,'(a32,0PF17.8,a3)') 'LDA Correlation', etclda, 'Ry'   !,etclda
  !
  IF (igcc/=0) THEN
     WRITE(stdout,'(a32,0PF17.8,a3)') 'E_c^gc', etcgc, 'Ry'
  END IF
  IF (dft_is_nonlocc()) THEN
     IF (vdW_analysis==0) WRITE(stdout,'(a32,0PF17.8,a3)') 'E_c^nl', etcnl , 'Ry'
     IF (vdW_analysis==1) WRITE(stdout,'(a32,0PF17.8,a3)') 'E_c^nl(Alpha)', etcnl , 'Ry'
     IF (vdW_analysis==2) WRITE(stdout,'(a32,0PF17.8,a3)') 'E_c^nl(vdW)', etcnl , 'Ry'
  END IF
  WRITE(stdout,'(a32,0PF17.8,a3)') 'Exchange + Correlation', etx+etc, 'Ry'
  WRITE(stdout,'(a32,0PF17.8,a3)') 'T_c^LDA', ttclda, 'Ry'
  IF (igcc /= 0) THEN
     WRITE(stdout,'(a32,0PF17.8,a3)') 'T_c^gc', ttcgc, 'Ry'
     WRITE(stdout,'(a32,0PF17.8,a3)') 'Kinetic-correlation Energy', ttclda+ttcgc, 'Ry'
  END IF
  IF (dft_is_nonlocc()) THEN
     IF (vdW_analysis==0) WRITE(stdout,'(a32,0PF17.8,a3)') 'T_c^nl', etcnl - etcnlncc, 'Ry'
     IF (vdW_analysis==1) WRITE(stdout,'(a32,0PF17.8,a3)') 'T_c^nl(Alpha)', etcnl - etcnlncc, 'Ry'
     IF (vdW_analysis==2) WRITE(stdout,'(a32,0PF17.8,a3)') 'T_c^nl(vdW)', etcnl - etcnlncc, 'Ry'
     !
     IF (vdW_analysis==0) WRITE(stdout,'(a32,0PF17.8,a3)') 'Kinetic-correlation Energy', ttclda+(etcnl-etcnlncc), 'Ry'
  END IF
  !
  DEALLOCATE( vofrcc )
  !
  ! ... Non-local correlation energy density
  IF (dft_is_nonlocc() .AND. lplot) THEN
     !
     ALLOCATE( q0(dfftp%nnr), dq0_drho(dfftp%nnr,nspin), dq0_dgradrho(dfftp%nnr,nspin) )
     ALLOCATE( tot_grad_rho(3,dfftp%nnr), grad_rho(3,dfftp%nnr,nspin) )
     ALLOCATE( thetas(dfftp%nnr,Nqs) )
     ALLOCATE( u_vdW(dfftp%nnr,Nqs), potential_vdW(dfftp%nnr,nspin) )
     ALLOCATE( thetasp(dfftp%nnr,Nqs), thetasm(dfftp%nnr,Nqs) )
     ALLOCATE( up_vdW(dfftp%nnr,Nqs), um_vdW(dfftp%nnr,Nqs) )
     !
     ! ... Here we calculate the gradient in reciprocal space using FFT.
     CALL fft_gradient_r2r( dfftp, tot_rho, g, tot_grad_rho )
     DO is = 1, nspin
        CALL fft_gradient_r2r( dfftp, rhoout(:,is), g, grad_rho(:,:,is) )
     ENDDO
     !
     dq0_drho = 0._DP
     dq0_dgradrho = 0._DP
     ccp = 1._DP+dcc
     ccm = 1._DP-dcc
     !
     ! ... thetas
     IF (nspin == 1) THEN
        CALL get_q0cc_on_grid( 1._DP, lecnl_qx, tot_rho, tot_grad_rho, q0, thetas  )
        CALL get_q0cc_on_grid( ccp,   lecnl_qx, tot_rho, tot_grad_rho, q0, thetasp )
        CALL get_q0cc_on_grid( ccm,   lecnl_qx, tot_rho, tot_grad_rho, q0, thetasm )
     ELSEIF (nspin == 2) THEN
        CALL get_q0cc_on_grid_spin( 1._DP, lecnl_qx, tot_rho, rhoout(:,1), rhoout(:,2), &
                            tot_grad_rho, grad_rho(:,:,1), grad_rho(:,:,2), q0, thetas  )
        CALL get_q0cc_on_grid_spin( ccp,   lecnl_qx, tot_rho, rhoout(:,1), rhoout(:,2), &
                            tot_grad_rho, grad_rho(:,:,1), grad_rho(:,:,2), q0, thetasp )
        CALL get_q0cc_on_grid_spin( ccm,   lecnl_qx, tot_rho, rhoout(:,1), rhoout(:,2), &
                            tot_grad_rho, grad_rho(:,:,1), grad_rho(:,:,2), q0, thetasm )
     ENDIF
     !
     u_vdW(:,:)  = thetas(:,:)
     up_vdW(:,:) = thetasp(:,:)
     um_vdW(:,:) = thetasm(:,:)
     !
     DO iq = 1, Nqs
        CALL invfft( 'Rho', thetas(:,iq),  dfftp )
        CALL invfft( 'Rho', thetasp(:,iq), dfftp )
        CALL invfft( 'Rho', thetasm(:,iq), dfftp )
     ENDDO
     !
     CALL vdW_DF_energy( u_vdW, Ec_nl )
     CALL mp_sum( Ec_nl, intra_bgrp_comm )
     WRITE(stdout,*) '     Non-local energy : ', Ec_nl
     CALL vdW_DF_energy( up_vdW, Ec_nl )
     CALL vdW_DF_energy( um_vdW, Ec_nl )
     !
     DO iq = 1, Nqs
        CALL invfft( 'Rho', u_vdW(:,iq),  dfftp )
        CALL invfft( 'Rho', up_vdW(:,iq), dfftp )
        CALL invfft( 'Rho', um_vdW(:,iq), dfftp )
     ENDDO
     !
     IF (nspin == 1) THEN
        CALL vdW_DF_potential( q0,dq0_drho(:,1), dq0_dgradrho(:,1), tot_grad_rho,    &
                                                        u_vdW, potential_vdW(:,1) )
     ELSEIF (nspin == 2) THEN
        CALL vdW_DF_potential( q0,dq0_drho(:,1), dq0_dgradrho(:,1), grad_rho(:,:,1), &
                                                        u_vdW, potential_vdW(:,1) )
        CALL vdW_DF_potential( q0,dq0_drho(:,2), dq0_dgradrho(:,2), grad_rho(:,:,2), &
                                                        u_vdW, potential_vdW(:,2) )
     ENDIF
     !
     ecnl%of_r(:,1) = 0._DP
     tcnl%of_r(:,1) = 0._DP
     !
     ALLOCATE( ecnl_c(dfftp%nnr) )
     ALLOCATE( tcnl_c(dfftp%nnr) )
     !
     ecnl_c = DCMPLX(0._DP,0._DP)
     tcnl_c = DCMPLX(0._DP,0._DP)
     etcnl_check = 0._DP
     ttcnl_check = 0._DP
     !
     DO ir = 1, dfftp%nnr
        arhox = ABS(tot_rho(ir))
        IF (arhox > vanishing_charge) THEN
           DO iq = 1, Nqs
              ecnl_c(ir) = ecnl_c(ir) + thetas(ir,iq)*u_vdW(ir,iq)
              tcnl_c(ir) = tcnl_c(ir) - thetas(ir,iq)*u_vdW(ir,iq) & 
                                      - 0.5_DP/dcc*(thetasp(ir,iq)*up_vdW(ir,iq) &
                                      - thetasm(ir,iq)*um_vdW(ir,iq))
           ENDDO
           !
           etcnl_check = etcnl_check+ecnl_c(ir)
           ttcnl_check = ttcnl_check+tcnl_c(ir)
           !
        ENDIF
     ENDDO
     !
     etcnl_check = e2*0.5_DP*omega*etcnl_check/DBLE(dfftp%nr1*dfftp%nr2*dfftp%nr3)
     ttcnl_check = e2*0.5_DP*omega*ttcnl_check/DBLE(dfftp%nr1*dfftp%nr2*dfftp%nr3)
     !
     CALL mp_sum( etcnl_check, intra_bgrp_comm )
     CALL mp_sum( ttcnl_check, intra_bgrp_comm )
     !
     DEALLOCATE( q0, dq0_drho, dq0_dgradrho, thetas )
     DEALLOCATE( tot_grad_rho, grad_rho )
     DEALLOCATE( u_vdW )
     DEALLOCATE( thetasp, thetasm, up_vdW, um_vdW )
     !
     ecnl%of_r(:,1) = e2*0.5_DP*DBLE(ecnl_c(:))
     tcnl%of_r(:,1) = e2*0.5_DP*DBLE(tcnl_c(:))
     !
     DEALLOCATE( ecnl_c, tcnl_c )
     WRITE(stdout,*) '     Summation of ecnl: ', etcnl_check
     WRITE(stdout,*) '     Summation of tcnl: ', ttcnl_check
  ENDIF
  !
  !
  IF ( lfock .OR. (lplot .AND. ltks) ) THEN
     IF (code_num == 1) THEN
        starting_wfc = 'file'
        CALL wfcinit()
     ELSEIF (code_num == 2) THEN
        CALL errore( 'ppacf', 'wavefunction not implemented for VASP postprocessing', 1)
     ELSE
        CALL errore( 'ppacf', 'code_num not implemented', 1 )
     ENDIF
  ENDIF
  ! Fock exchange energy from readin wavefunctions
  !
  IF (lfock) THEN
    nq1 = 0
    nq2 = 0
    nq3 = 0
    nbndproj = 0
    exxdiv_treatment = "gygi-baldereschi"
    ecutfock = ecutwfc
    !
    CALL xclib_set_exx_fraction( 1._DP )
    CALL enforce_input_dft( 'HF' )
    isnonlocc = dft_is_nonlocc()
    CALL xclib_set_auxiliary_flags( isnonlocc )
    !
    ALLOCATE( igk_buf(npwx), gk(npwx) )
    igk_k(:,:) = 0
    DO ik = 1, nks
       CALL gk_sort( xk(1,ik), ngm, g, gcutw, npw, igk_buf, gk )
       ngk(ik) = npw
       igk_k(1:npw,ik) = igk_buf(1:npw)
    ENDDO
    DEALLOCATE( igk_buf, gk )
    !
    !  CALL setup()
    CALL exx_grid_init()
    CALL exx_mp_init()
    CALL exx_div_check()
    !  CALL init_run()
    CALL exxinit( .FALSE. )
    !
    IF ( use_ace) THEN
       CALL aceinit( DOLOC = (local_thr > 0._DP) )
       fock2 = exxenergyace()
    ELSE
       fock2 = exxenergy2()
    ENDIF
    !
    WRITE(stdout, 9068) 0.5_DP*fock2
9068 FORMAT( '     Fock energy               =',0PF17.8,' Ry' )
  ENDIF
  !
  ! output data in 3D
  IF (lplot) THEN
     IF (code_num == 2) THEN
        IF (nspin == 1) THEN
           filplot = TRIM(prefix)//'.chg'
           plot_num = 2
           CALL dcopy( dfftp%nnr, rho%of_r(:,1), 1, vltot, 1 )
           CALL punch_plot( filplot, plot_num, 0., 0., 0., 0., 0., 0, 0, 0, .FALSE. )
        ELSEIF (nspin == 2) THEN
           filplot = TRIM(prefix)//'.chg1'
           plot_num = 2
           CALL dcopy( dfftp%nnr, rho%of_r(:,1), 1, vltot, 1 )
           CALL punch_plot( filplot, plot_num, 0., 0., 0., 0., 0., 0, 0, 0, .FALSE. )
           filplot = TRIM(prefix)//'.chg2'
           plot_num = 2
           CALL dcopy( dfftp%nnr, rho%of_r(:,2), 1, vltot, 1 )
           CALL punch_plot( filplot, plot_num, 0., 0., 0., 0., 0., 0, 0, 0, .FALSE. )
        ENDIF
     ENDIF
     !
     IF (ltks) THEN
        ALLOCATE( kin_r(dfftp%nnr,nspin) )
        ! CALL init_run()
        CALL sum_band_kin( kin_r )
        !
        IF (nspin == 1) THEN
           filplot = TRIM(prefix)//'.tks'
           plot_num = 2
           CALL dcopy(dfftp%nnr, kin_r(:,1), 1, vltot, 1 )
           CALL punch_plot( filplot, plot_num, 0., 0., 0., 0., 0., 0, 0, 0, .FALSE. )
        ELSEIF (nspin == 2) THEN
           filplot = TRIM(prefix)//'.tks1'
           plot_num = 2
           CALL dcopy( dfftp%nnr, kin_r(:,1), 1, vltot, 1 )
           CALL punch_plot( filplot, plot_num, 0., 0., 0., 0., 0., 0, 0, 0, .FALSE. )
           filplot = TRIM(prefix)//'.tks2'
           plot_num = 2
           CALL dcopy( dfftp%nnr, kin_r(:,2), 1, vltot, 1 )
           CALL punch_plot( filplot, plot_num, 0., 0., 0., 0., 0., 0, 0, 0, .FALSE. )
        ENDIF
        DEALLOCATE( kin_r )
     ENDIF
     !
     filplot = TRIM(prefix)//'.exlda'
     plot_num = 2
     CALL dcopy( dfftp%nnr, exlda%of_r(:,1), 1, vltot, 1 )
     CALL punch_plot( filplot, plot_num, 0., 0., 0., 0., 0., 0, 0, 0, .FALSE. )
     !
     filplot = TRIM(prefix)//'.eclda'
     plot_num = 2
     CALL dcopy( dfftp%nnr, eclda%of_r(:,1), 1, vltot, 1 )
     CALL punch_plot( filplot, plot_num, 0., 0., 0., 0., 0., 0, 0, 0, .FALSE. )
     !
     filplot = TRIM(prefix)//'.tclda'
     plot_num = 2
     CALL dcopy( dfftp%nnr, tclda%of_r(:,1), 1, vltot, 1 )
     CALL punch_plot( filplot, plot_num, 0., 0., 0., 0., 0., 0, 0, 0, .FALSE. )
     !
     IF (igcx /= 0) THEN
        filplot = TRIM(prefix)//'.exgc'
        plot_num = 2
        CALL dcopy( dfftp%nnr, exgc%of_r(:,1), 1, vltot, 1 )
        CALL punch_plot( filplot, plot_num, 0., 0., 0., 0., 0., 0, 0, 0, .FALSE.)
        CALL destroy_scf_type( exgc )
     ENDIF
     IF (dft_is_nonlocc()) THEN
        IF (vdW_analysis==0) filplot = TRIM(prefix)//'.ecnl'
        IF (vdW_analysis==1) filplot = TRIM(prefix)//'.ecnl_Alpha'
        IF (vdW_analysis==2) filplot = TRIM(prefix)//'.ecnl_vdW'
        plot_num = 2
        CALL dcopy( dfftp%nnr, ecnl%of_r(:,1), 1, vltot, 1 )
        CALL punch_plot( filplot, plot_num, 0., 0., 0., 0., 0., 0, 0, 0, .FALSE. )
        !
        IF (vdW_analysis==0) filplot = TRIM(prefix)//'.tcnl'
        IF (vdW_analysis==1) filplot = TRIM(prefix)//'.tcnl_Alpha'
        IF (vdW_analysis==2) filplot = TRIM(prefix)//'.tcnl_vdW'
        plot_num = 2
        CALL dcopy( dfftp%nnr, tcnl%of_r(:,1), 1, vltot, 1 )
        CALL punch_plot( filplot, plot_num, 0., 0., 0., 0., 0., 0, 0, 0, .FALSE. )
        !
        IF (nspin == 1) THEN
           IF (vdW_analysis==0) filplot = TRIM(prefix)//'.vcnl'
           IF (vdW_analysis==1) filplot = TRIM(prefix)//'.vcnl_Alpha'
           IF (vdW_analysis==2) filplot = TRIM(prefix)//'.vcnl_vdW'
           plot_num = 2
           CALL dcopy( dfftp%nnr, potential_vdW(:,1), 1, vltot, 1 )
           CALL punch_plot( filplot, plot_num, 0., 0., 0., 0., 0., 0, 0, 0, .FALSE. )
        ELSEIF (nspin == 2) THEN
           IF (vdW_analysis==0) filplot = TRIM(prefix)//'.vcnl1'
           IF (vdW_analysis==1) filplot = TRIM(prefix)//'.vcnl1_Alpha'
           IF (vdW_analysis==2) filplot = TRIM(prefix)//'.vcnl1_vdW'
           plot_num = 2
           CALL dcopy( dfftp%nnr, potential_vdW(:,1), 1, vltot, 1 )
           CALL punch_plot( filplot, plot_num, 0., 0., 0., 0., 0., 0, 0, 0, .FALSE. )
           IF (vdW_analysis==0) filplot = TRIM(prefix)//'.vcnl2'
           IF (vdW_analysis==1) filplot = TRIM(prefix)//'.vcnl2_Alpha'
           IF (vdW_analysis==2) filplot = TRIM(prefix)//'.vcnl2_vdW'
           plot_num = 2
           CALL dcopy( dfftp%nnr, potential_vdW(:,2), 1, vltot, 1 )
           CALL punch_plot( filplot, plot_num, 0., 0., 0., 0., 0., 0, 0, 0, .FALSE. )
        ENDIF
        !
        DEALLOCATE( potential_vdW )
        !
        CALL destroy_scf_type( ecnl )
        CALL destroy_scf_type( tcnl )
        !
     ELSEIF (igcc /= 0) THEN
        !
        filplot = TRIM(prefix)//'.ecgc'
        plot_num = 2
        CALL dcopy( dfftp%nnr, ecgc%of_r(:,1), 1, vltot, 1 )
        CALL punch_plot( filplot, plot_num, 0., 0., 0., 0., 0., 0, 0, 0, .FALSE. )
        filplot = TRIM(prefix)//'.tcgc'
        plot_num = 2
        CALL dcopy( dfftp%nnr, tcgc%of_r(:,1), 1, vltot, 1 )
        CALL punch_plot( filplot, plot_num, 0., 0., 0., 0., 0., 0, 0, 0, .FALSE. )
        CALL destroy_scf_type( ecgc )
        CALL destroy_scf_type( tcgc )
        !
     ENDIF
     !
     CALL stop_clock( 'acf_etxclambda' )
     WRITE( stdout, '(//5x,"exiting subroutine acf ..."/)')
  ENDIF
  !
  CALL destroy_scf_type( exlda )
  CALL destroy_scf_type( eclda )
  CALL destroy_scf_type( tclda )
  !
9091 FORMAT( 0PF17.8,0PF17.8,0PF17.8,0PF17.8,0PF17.8,0PF17.8 )
9092 FORMAT( 0PF17.8,0PF17.8,0PF17.8,0PF17.8 )
9095 FORMAT( 0PF17.8,0PF17.8)
!9090 FORMAT(/' ACF coupling-constant         =',0PF17.8,'   ' &
!            /'     Exc_lambda                =',0PF17.8,' Ry' )
!9091 FORMAT(/'     Non-local contribution    =',0PF17.8,' Ry' )
9093 FORMAT(/' delta coupling constant        =',0PE17.4E3,' ')
  CALL environment_end('ppacf')
  CALL stop_pp
!
END PROGRAM do_ppacf
!
!
!-----------------------------------------------------------------------
SUBROUTINE pwcc( rs, cc, ec, vc, ec_l )
  !-----------------------------------------------------------------------
  !! This subroutine is adapted from subroutine pw from 
  !! Modules/functionals.f90
  !
  !! iflag=1: J.P. Perdew and Y. Wang, PRB 45, 13244 (1992)
  !
  USE kinds,  ONLY : DP
  !
  IMPLICIT NONE
  !
  REAL(DP), INTENT(IN) :: rs, cc
  REAL(DP), INTENT(OUT):: ec, vc, ec_l
  INTEGER  :: iflag
  REAL(DP), PARAMETER :: a=0.031091d0, b1=7.5957d0, b2=3.5876d0, &
                         c0=a, c1=0.046644d0, c2=0.00664d0, c3=0.01043d0, &
                         d0=0.4335d0, d1=1.4408d0
  REAL(DP) :: lnrs, rs12, rs32, rs2, om, dom, olog
  REAL(DP) :: lrs, lrs12, lrs32, lrs2, lom, dlom, lolog
  REAL(DP) :: lec, dlec
  REAL(DP) :: a1(2), b3(2), b4(2)
  DATA a1 / 0.21370d0, 0.026481d0 /, b3 / 1.6382d0, -0.46647d0 /, &
       b4 / 0.49294d0, 0.13354d0 /
  !
  iflag = 1
  ! interpolation formula
  rs12 = SQRT(rs)
  rs32 = rs * rs12
  rs2  = rs * rs
  !
  lrs   = cc * rs
  lrs12 = SQRT(lrs)
  lrs32 = lrs * lrs12
  lrs2  = lrs * lrs
  !
  om   = 2.d0 * a * (b1 * rs12 + b2 * rs + b3 (iflag) * rs32 + b4 ( &
         iflag) * rs2)
  lom  = 2.d0 * a * (b1 * lrs12 + b2 * lrs + b3 (iflag) * lrs32 + b4 ( &
         iflag) * lrs2)
  dom  = 2.d0 * a * (0.5d0 * b1 * rs12 + b2 * rs + 1.5d0 * b3 ( &
         iflag) * rs32 + 2.d0 * b4 (iflag) * rs2)
  dlom = 2.d0 * a * (0.5d0*b1*rs/lrs12+b2*rs+1.5d0*b3(iflag)*lrs12*rs &
         +2.d0*b4(iflag)*lrs*rs)
  !
  olog  = LOG(1.d0 + 1.0d0 /  om)
  lolog = LOG(1.d0 + 1.0d0 / lom) 
  !
  ec = - 2.d0 * a * (1.d0 + a1(iflag) * rs) * olog
  lec = - 2.d0 * a * (1.d0 + a1(iflag) * lrs) * lolog
  dlec = - 2.d0*a*a1(iflag)*rs*lolog  &
         + 2.d0*a*(1.d0+a1(iflag)*lrs)*dlom/(lom*(lom+1.d0))
  ec_l = 2.d0*cc*lec+cc*cc*dlec
  !
  vc = - 2.d0 * a * (1.d0 + 2.d0 / 3.d0 * a1(iflag) * rs) &
       * olog - 2.d0 / 3.d0 * a * (1.d0 + a1(iflag) * rs) &
       * dom / (om * (om + 1.d0) )
  !
  RETURN
  !
END SUBROUTINE pwcc
!
!
!-----------------------------------------------------------------------
SUBROUTINE pwcc_spin( rs, cc, zeta, ec, vcup, vcdw, ec_l )
  !-----------------------------------------------------------------------
  !! This subroutine is adapted from subroutine pw_spin from 
  !! Modules/lsda_functionals.f90
  !
  !! J.P. Perdew and Y. Wang, PRB 45, 13244 (1992)
  !
  USE kinds, ONLY : DP
  !
  IMPLICIT NONE
  !
  REAL(DP) :: rs,cc, zeta, ec, vcup, vcdw,ec_l
  ! xc PARAMETERs, unpolarised
  REAL(DP) :: a, a1, b1, b2, b3, b4, c0, c1, c2, c3, d0, d1
  PARAMETER (a = 0.031091d0, a1 = 0.21370d0, b1 = 7.5957d0, b2 = &
       3.5876d0, b3 = 1.6382d0, b4 = 0.49294d0, c0 = a, c1 = 0.046644d0, &
       c2 = 0.00664d0, c3 = 0.01043d0, d0 = 0.4335d0, d1 = 1.4408d0)
  ! xc PARAMETERs, polarised
  REAL(DP) :: ap, a1p, b1p, b2p, b3p, b4p, c0p, c1p, c2p, c3p, d0p, &
       d1p
  PARAMETER (ap = 0.015545d0, a1p = 0.20548d0, b1p = 14.1189d0, b2p &
       = 6.1977d0, b3p = 3.3662d0, b4p = 0.62517d0, c0p = ap, c1p = &
       0.025599d0, c2p = 0.00319d0, c3p = 0.00384d0, d0p = 0.3287d0, d1p &
       = 1.7697d0)
  ! xc PARAMETERs, antiferro
  REAL(DP) :: aa, a1a, b1a, b2a, b3a, b4a, c0a, c1a, c2a, c3a, d0a, &
       d1a
  PARAMETER (aa = 0.016887d0, a1a = 0.11125d0, b1a = 10.357d0, b2a = &
       3.6231d0, b3a = 0.88026d0, b4a = 0.49671d0, c0a = aa, c1a = &
       0.035475d0, c2a = 0.00188d0, c3a = 0.00521d0, d0a = 0.2240d0, d1a &
       = 0.3969d0)
  REAL(DP) :: fz0
  PARAMETER (fz0 = 1.709921d0)
  REAL(DP) :: rs12, rs32, rs2, zeta2, zeta3, zeta4, fz, dfz
  REAL(DP) :: om, dom, olog, epwc, vpwc
  REAL(DP) :: omp, domp, ologp, epwcp, vpwcp
  REAL(DP) :: oma, doma, ologa, alpha, vpwca
  REAL(DP) :: lrs,lrs12,lrs32,lrs2
  REAL(DP) :: lom,dlom,lolog,lepwc,dlepwc
  REAL(DP) :: lomp,dlomp,lologp,lepwcp,dlepwcp
  REAL(DP) :: loma,dloma,lologa,lalpha,dlalpha
  REAL(DP) :: lec,dlec
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
  !
  rs12 = SQRT(rs)
  rs32 = rs * rs12
  rs2  = rs * rs
  !
  lrs   = cc*rs
  lrs12 = SQRT(lrs)
  lrs32 = lrs*lrs12
  lrs2  = lrs*lrs
  !
  ! unpolarised
  om   = 2.d0 * a * (b1 * rs12 + b2 * rs + b3 * rs32 + b4 * rs2)
  dom  = 2.d0 * a * (0.5d0 * b1 * rs12 + b2 * rs + 1.5d0 * b3 * rs32 &
        +2.d0 * b4 * rs2)
  olog = LOG(1.d0 + 1.0d0 / om)
  epwc = - 2.d0 * a * (1.d0 + a1 * rs) * olog
  vpwc = - 2.d0 * a * (1.d0 + 2.d0 / 3.d0 * a1 * rs) * olog - 2.d0 / &
           3.d0 * a * (1.d0 + a1 * rs) * dom / (om * (om + 1.d0) )
  ! polarized
  omp  = 2.d0 * ap * (b1p * rs12 + b2p * rs + b3p * rs32 + b4p * rs2)
  domp = 2.d0 * ap * (0.5d0 * b1p * rs12 + b2p * rs + 1.5d0 * b3p * &
         rs32 + 2.d0 * b4p * rs2)
  ologp = LOG(1.d0 + 1.0d0 / omp)
  epwcp = - 2.d0 * ap * (1.d0 + a1p * rs) * ologp
  vpwcp = - 2.d0 * ap * (1.d0 + 2.d0 / 3.d0 * a1p * rs) * ologp - &
            2.d0 / 3.d0 * ap * (1.d0 + a1p * rs) * domp / (omp * (omp + 1.d0) )
  ! antiferro
  oma = 2.d0 * aa * (b1a * rs12 + b2a * rs + b3a * rs32 + b4a * rs2)
  doma = 2.d0 * aa * (0.5d0 * b1a * rs12 + b2a * rs + 1.5d0 * b3a * &
         rs32 + 2.d0 * b4a * rs2)
  ologa = LOG(1.d0 + 1.0d0 / oma)
  alpha = 2.d0 * aa * (1.d0 + a1a * rs) * ologa
  vpwca = + 2.d0 * aa * (1.d0 + 2.d0 / 3.d0 * a1a * rs) * ologa + &
            2.d0 / 3.d0 * aa * (1.d0 + a1a * rs) * doma / (oma * (oma + 1.d0) )
  ! coupling constant depandent
  lom = 2.d0*a*(b1*lrs12+b2*lrs+b3*lrs32+b4*lrs2)
  dlom = 2.d0 * a * (0.5d0*b1*rs/lrs12+b2*rs+1.5d0*b3*lrs12*rs &
        +2.d0*b4*lrs*rs)
  lolog = LOG(1.d0+1.0d0/lom)
  lepwc = -2.d0*a*(1.d0+a1*lrs)*lolog 
  dlepwc = -2.d0*a*a1*rs*lolog &
           +2.d0*a*(1.d0+a1*lrs)*dlom/(lom*(lom+1.d0))

  lomp = 2.d0*ap*(b1p*lrs12+b2p*lrs+b3p*lrs32+b4p*lrs2)
  dlomp = 2.d0*ap*(0.5d0*b1p*rs/lrs12+b2p*rs+1.5d0*b3p*lrs12*rs &
         +2.d0*b4p*lrs*rs)
  lologp = LOG(1.d0+1.d0/lomp)
  lepwcp = -2.d0*ap*(1.d0+a1p*lrs)*lologp
  dlepwcp = -2.d0*ap*a1p*rs*lologp+2.d0*ap*(1.d0+a1p*lrs)*dlomp/(lomp*(lomp+1.d0))

  loma = 2.d0*aa*(b1a*lrs12+b2a*lrs+b3a*lrs32+b4a*lrs2)
  dloma = 2.d0*aa*(0.5d0*b1a*rs/lrs12+b2a*rs+1.5d0*b3a*lrs12*rs &
         +2.d0*b4a*lrs*rs)
  lologa = LOG(1.d0+1.d0/loma)
  lalpha = 2.d0*aa*(1.d0+a1a*lrs)*lologa
  dlalpha = 2.d0*aa*a1a*rs*lologa  &
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
  !
  vcdw = vpwc + vpwca * fz * (1.d0 - zeta4) / fz0 + (vpwcp - vpwc) &
       * fz * zeta4 - (alpha / fz0 * (dfz * (1.d0 - zeta4) - 4.d0 * fz * &
       zeta3) + (epwcp - epwc) * (dfz * zeta4 + 4.d0 * fz * zeta3) ) &
       * (1.d0 + zeta)
  !      endif
  !
  RETURN
  !
END SUBROUTINE pwcc_spin
!
!----------------------------------------
! ####################################################################
!                          |              |
!                          |  PPACF_INFO  |
!                          |____________ _|
!
SUBROUTINE ppacf_info
  !
  USE io_global,  ONLY : stdout
  !
  IMPLICIT NONE
  !
  !
  WRITE(stdout,'(/)')
  WRITE(stdout,'(5x,"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")')
  WRITE(stdout,'(5x,"%                                                                      %")')
  WRITE(stdout,'(5x,"% You are using PPACF, please cite the following paper:                %")') 
  WRITE(stdout,'(5x,"%                                                                      %")')
  WRITE(stdout,'(5x,"%   Y. Jiao, E. Schr\""oder, and P. Hyldgaard, PRB 97, 085115 (2018).   %")')
  WRITE(stdout,'(5x,"%                                                                      %")')
  WRITE(stdout,'(5x,"% Illustrations of use for analysis (vdW-DFs and vdW-DF-hybrids):      %")')
  WRITE(stdout,'(5x,"%                                                                      %")')
  WRITE(stdout,'(5x,"%   Y. Jiao, E. Schr\""oder, and P. Hyldgaard, JCP 148, 194115 (2018).  %")')
  WRITE(stdout,'(5x,"%   P. Hyldgaard, Y. Jiao, and V. Shukla, IOP JPCM 32, 393001 (2020).   %")')
  WRITE(stdout,'(5x,"%                                                                      %")')
  WRITE(stdout,'(5x,"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")')
  WRITE(stdout,'(/)')
  !
  !
END SUBROUTINE

