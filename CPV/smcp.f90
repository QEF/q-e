!
! Copyright (C) 2002-2005 Quantum-ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "f_defs.h"
!
!----------------------------------------------------------------------------
SUBROUTINE smdmain( tau, fion_out, etot_out, nat_out )
  !----------------------------------------------------------------------------
  !
  ! ... main subroutine for SMD by Yosuke Kanai
  !
  USE kinds,                    ONLY : DP
  USE control_flags,            ONLY : iprint, isave, thdyn, tpre, tbuff, &
                                       iprsta, trhor, tfor, tvlocw, trhow,     &
                                       taurdr, tprnfor, tsdc, ndr, ndw, nbeg,  &
                                       nomore, tsde, tortho, tnosee, tnosep,   &
                                       trane, tranp, tsdp, tcp, tcap, ampre,   &
                                       amprp, tnoseh, tolp,ortho_eps, ortho_max
  USE atom,                     ONLY : nlcc
  USE core,                     ONLY : nlcc_any, rhoc
  USE uspp_param,               ONLY : nhm
  USE uspp,                     ONLY : nkb, vkb, rhovan => becsum, deeq
  USE cvan,                     ONLY : nvb
  USE energies,                 ONLY : eht, epseu, exc, etot, eself, enl, ekin, &
                                       esr, print_energies
  USE electrons_base,           ONLY : nbspx, nbsp, f, nspin
  USE electrons_base,           ONLY : nel, iupdwn, nupdwn
  USE gvecp,                    ONLY : ngm
  USE gvecs,                    ONLY : ngs
  USE gvecb,                    ONLY : ngb
  USE gvecw,                    ONLY : ngw, ngwt
  USE reciprocal_vectors,       ONLY : gstart, mill_l, gzero
  USE ions_base,                ONLY : na, nat, pmass, nax, nsp, rcmax
  USE ions_base,                ONLY : ind_srt, ions_vel, ions_cofmass, &
                                       ions_kinene, ions_temp
  USE ions_base,                ONLY : ions_thermal_stress, ions_vrescal, &
                                       fricp, greasp, iforce
  USE grid_dimensions,          ONLY : nr1, nr2, nr3, nnrx
  USE cell_base,                ONLY : ainv, a1, a2, a3, r_to_s, celldm, ibrav
  USE cell_base,                ONLY : omega, alat, frich, greash, press
  USE cell_base,                ONLY : h, hold, hnew, velh,  deth, wmass, &
                                       s_to_r, iforceh, cell_force
  USE cell_base,                ONLY : thdiag, tpiba2
  USE smooth_grid_dimensions,   ONLY : nnrsx, nr1s, nr2s, nr3s
  USE smallbox_grid_dimensions, ONLY : nr1b, nr2b, nr3b
  USE local_pseudo,             ONLY : allocate_local_pseudo
  USE io_global,                ONLY : io_global_start, stdout, ionode
  USE dener,                    ONLY : detot
  USE derho,                    ONLY : drhog, drhor
  USE cdvan,                    ONLY : dbec, drhovan
  USE stre,                     ONLY : stress
  USE gvecw,                    ONLY : ggp, ecutw
  USE gvecp,                    ONLY : ecutp
  USE restart_file,             ONLY : writefile, readfile
  USE parameters,               ONLY : nacx, natx, nsx, nbndxx, nhclm
  USE constants,                ONLY : pi, factem, au_gpa, au_ps, gpa_au
  USE io_files,                 ONLY : psfile, pseudo_dir, smwout, outdir
  USE wave_base,                ONLY : wave_steepest, wave_verlet
  USE wave_base,                ONLY : wave_speed2, frice, grease
  USE control_flags,            ONLY : conv_elec, tconvthrs
  USE check_stop,               ONLY : check_stop_now
  USE cpr_subroutines,          ONLY : print_atomic_var, print_lambda, &
                                       compute_stress, elec_fakekine2
  USE ions_positions,           ONLY : tau0, velsp
  USE ions_positions,           ONLY : ions_hmove, ions_move
  USE cell_base,                ONLY : cell_kinene, cell_move, cell_gamma, &
                                       cell_hmove
  USE ions_nose,                ONLY : gkbt, qnp, tempw, kbt, nhpcl, &
                                       ndega, nhpdim, atm2nhp, ekin2nhp
  USE cell_nose,                ONLY : xnhh0, xnhhm, xnhhp, vnhh, temph, qnh
  USE time_step,                ONLY : delt
  USE cp_electronic_mass,       ONLY : emass, emass_cutoff, emass_precond
  USE electrons_nose,           ONLY : qne, ekincw, electrons_nose_nrg, &
                                       electrons_nose_shiftvar
  USE path_variables,           ONLY : smx, smd_p, smd_ptr, smd_lmfreq, &
                                       smd_tol, smd_codfreq, smd_forfreq, &
                                       smd_cp, smd_lm, smd_opt, smd_linr, &
                                       smd_polm, smd_stcd, smd_kwnp,      &
                                       smd_maxlm, smd_ene_ini, smd_ene_fin
  USE smd_rep,                  ONLY : rep, rep_el, deallocate_smd_rep
  USE smd_ene,                  ONLY : etot_ar, ekin_ar, eht_ar, epseu_ar, &
                                       exc_ar, esr_ar, deallocate_smd_ene
  USE from_restart_module,      ONLY : from_restart
  USE runcp_module,             ONLY : runcp_uspp
  USE phase_factors_module,     ONLY : strucf
  USE cp_main_variables,        ONLY : ei1, ei2, ei3, eigr, sfac, irb, taub, &
                                       eigrb, rhog, rhor, rhos, becdr, bephi, &
                                       becp, ema0bg, allocate_mainvar, nfi
  USE mp_global,                ONLY : mp_global_start
  USE mp,                       ONLY : mp_sum
  USE fft_base,                 ONLY : dfftp
  !
#if ! defined __NOSMD
  !
  IMPLICIT NONE
  !
  ! output variables
  !
  INTEGER :: nat_out
  REAL(DP) :: tau( 3, nat_out, 0:* )
  REAL(DP) :: fion_out( 3, nat_out, 0:* )
  REAL(DP) :: etot_out( 0:* )
  !
  !
  ! control variables
  !
  LOGICAL tbump
  LOGICAL tfirst, tlast
  LOGICAL tstop
  LOGICAL tconv
  REAL(DP) :: delta_etot
  !
  !
  ! work variables
  !
  COMPLEX(DP) :: speed
  REAL(DP) :: vnhp(nhclm,smx), xnhp0(nhclm,smx), xnhpm(nhclm,smx), &
                    tempp(smx), xnhe0(smx), fccc(smx), xnhem(smx), vnhe(smx),  &
                    epot(smx), xnhpp(nhclm,smx), xnhep(smx), epre(smx), enow(smx),   &
                    econs(smx), econt(smx), ccc(smx)
  REAL(DP) :: temps(nsx) 
  REAL(DP) :: verl1, verl2, verl3, anor, saveh, tps, bigr, dt2,            &
                    dt2by2, twodel, gausp, dt2bye, dt2hbe, savee, savep
  REAL(DP) :: ekinc0(smx), ekinp(smx), ekinpr(smx), ekincm(smx),   &
                    ekinc(smx), pre_ekinc(smx), enthal(smx)
  !
  INTEGER  ::is, nacc, ia, j, iter, i, isa, ipos
  !
  !
  ! work variables, 2
  !
  REAL(DP) :: fcell(3,3), hgamma(3,3)
  REAL(DP) :: cdm(3)
  REAL(DP) :: qr(3)
  REAL(DP) :: temphh(3,3)
  REAL(DP) :: thstress(3,3) 
  !
  INTEGER        :: k, ii, l, m
  REAL(DP) :: ekinh, alfar, temphc, alfap, temp1, temp2, randy
  REAL(DP) :: ftmp

  CHARACTER(len=256) :: filename
  CHARACTER(len=256) :: dirname
  INTEGER            :: strlen, dirlen
  CHARACTER(len=2)   :: unitsmw, unitnum 
  !
  REAL(DP) :: t_arc_pre, t_arc_now, t_arc_tot
  INTEGER        :: sm_k,sm_file,sm_ndr,sm_ndw,unico,unifo,unist
  INTEGER        :: smpm,con_ite
  REAL(DP) :: workvec(3,natx,nsx)
  !
  REAL(DP), ALLOCATABLE :: deviation(:)
  REAL(DP), ALLOCATABLE :: maxforce(:)
  REAL(DP), ALLOCATABLE :: arc_now(:)
  REAL(DP), ALLOCATABLE :: arc_tot(:)
  REAL(DP), ALLOCATABLE :: arc_pre(:)
  REAL(DP), ALLOCATABLE :: paraforce(:)
  REAL(DP), ALLOCATABLE :: err_const(:)
  INTEGER,        ALLOCATABLE :: pvvcheck(:)
  !
  TYPE(smd_ptr), ALLOCATABLE :: p_tau0(:)
  TYPE(smd_ptr), ALLOCATABLE :: p_taup(:)
  TYPE(smd_ptr), ALLOCATABLE :: p_tan(:)
  !
  REAL(DP), ALLOCATABLE :: mat_z(:,:,:)
  !
  !
  !     CP loop starts here
  !
  CALL start_clock( 'initialize' )

  tps     = 0.0d0

  !
  !     ==================================================================
  !
  WRITE(stdout,*) '                                                 '
  WRITE(stdout,*) ' ================================================'
  WRITE(stdout,*) '                                                 '
  WRITE(stdout,*) ' CPSMD : replica = 0 .. ', smd_p
  WRITE(stdout,*) '                                                 '
  WRITE(stdout,*) '                                                 '
  IF(smd_opt) THEN
     WRITE(stdout,*) ' CP optimizations for replicas 1, and 2  '
     WRITE(stdout,*) '                                                 '
  ENDIF
  WRITE(stdout,*) ' String Method : '
  WRITE(stdout,*) '                                                 '
  WRITE(stdout,*) '                                                 '
  WRITE(stdout,*) '      smd_lm    = ', smd_lm
  WRITE(stdout,*) '      CPMD    = ', smd_cp
  WRITE(stdout,*) '      smd_linr    = ', smd_linr 
  WRITE(stdout,*) '      smd_polm    = ', smd_polm, 'Pt = ', smd_kwnp
  WRITE(stdout,*) '      smd_stcd    = ', smd_stcd
  WRITE(stdout,*) '      SMfile  = ', smwout
  WRITE(stdout,*) '      TFOR  = ', tfor
  WRITE(stdout,*) '                                                 '
  IF(smd_lm) THEN
     WRITE(stdout,*) '      CONST. smd_tol  = ', smd_tol
     WRITE(stdout,*) '      Frequency   = ', smd_lmfreq
  ENDIF
  WRITE(stdout,*) '                                                 '
  WRITE(stdout,*) ' ================================================'
  !
  !
  !
  !     ===================================================================
  !
  !     general variables
  !
  tfirst = .TRUE.
  tlast  = .FALSE.
  nacc = 5
  !
  !
  !
  gausp = delt * SQRT(tempw/factem)
  twodel = 2.d0 * delt
  dt2 = delt * delt
  dt2by2 = .5d0 * dt2
  dt2bye = dt2/emass
  dt2hbe = dt2by2/emass
  smpm = smd_p -1
  !
  !
  IF(tbuff) THEN
     WRITE(stdout,*) "TBUFF set to .FALSE., Option not implemented with SMD"
     tbuff = .FALSE.
  ENDIF
  !
  !
  !     Opening files
  !
  !
  unico = 41
  unifo = 42
  unist = 43
  !
  IF (ionode) THEN

     dirlen = INDEX(outdir,' ') - 1
     filename = 'fort.8'
     IF( dirlen >= 1 ) THEN
        filename = outdir(1:dirlen) // '/' // filename
     END IF
     strlen  = INDEX(filename,' ') - 1
     WRITE( stdout, * ) ' UNIT8 = ', filename
     OPEN(unit=8, file=filename(1:strlen), status='unknown')

     filename = 'cod.out'
     IF( dirlen >= 1 ) THEN
        filename = outdir(1:dirlen) // '/' // filename
     END IF
     strlen  = INDEX(filename,' ') - 1
     OPEN(unit=unico, file=filename(1:strlen), status='unknown')

     filename = 'for.out'
     IF( dirlen >= 1 ) THEN
        filename = outdir(1:dirlen) // '/' // filename
     END IF
     strlen  = INDEX(filename,' ') - 1
     OPEN(unit=unifo, file=filename(1:strlen), status='unknown')

     filename = 'str.out'
     IF( dirlen >= 1 ) THEN
        filename = outdir(1:dirlen) // '/' // filename
     END IF
     strlen  = INDEX(filename,' ') - 1
     OPEN(unit=unist, file=filename(1:strlen), status='unknown')


     WRITE(unitsmw, '(i2)') smwout

     DO sm_k=1,smpm
        sm_file = smwout + sm_k

        unitnum = "00"

        IF(sm_k < 10) THEN
           WRITE(unitnum(2:2), '(i1)') sm_k
        ELSE
           WRITE(unitnum,'(i2)') sm_k
        ENDIF

        filename = 'rep.'
        IF( dirlen >= 1 ) THEN
           filename = outdir(1:dirlen) // '/' // filename
        END IF
        strlen  = INDEX(filename,' ') - 1

        OPEN(unit=sm_file,file=filename(1:strlen)//unitnum, status='unknown')
     ENDDO

  END IF

  !
  !     ==================================================================
  !     initialize g-vectors, fft grids
  !     ==================================================================
  !
  !      ... taus is calculated here.
  !
  CALL sminit( ibrav, celldm, ecutp, ecutw, ndr, nbeg, tfirst, delt, tps, iforce )
  !
  !
  CALL r_to_s( rep(0)%tau0,    rep(0)%taus,    na, nsp, ainv)
  CALL r_to_s( rep(smd_p)%tau0, rep(smd_p)%taus, na, nsp, ainv)
  !
  !
  WRITE( stdout,*) ' out from init'
  !
  !     more initialization requiring atomic positions
  !
  nax = MAXVAL( na( 1 : nsp ) )

  DO sm_k=1,smpm
     sm_file = smwout + sm_k
     IF( iprsta > 1 ) THEN
        IF(ionode) WRITE( sm_file,*) ' tau0 '
        IF(ionode) WRITE( sm_file,'(3f14.8)') ((rep(sm_k)%tau0(i,ia),i=1,3),ia=1,SUM(na(1:nsp)))
     ENDIF
  ENDDO

  !
  !     ==================================================================
  !     allocate and initialize nonlocal potentials
  !     ==================================================================

  CALL nlinit
  !

  WRITE( stdout,*) ' out from nlinit'

  !     ==================================================================
  !     allocation of all arrays not already allocated in init and nlinit
  !     ==================================================================
  !
  !
  WRITE( stdout,*) ' Allocation begun, smpm = ', smpm
  !
  CALL flush_unit( stdout )
  !
  DO sm_k=1,smpm
     ALLOCATE(rep_el(sm_k)%c0(ngw,nbspx))
     ALLOCATE(rep_el(sm_k)%cm(ngw,nbspx))
     ALLOCATE(rep_el(sm_k)%phi(ngw,nbspx))
     ALLOCATE(rep_el(sm_k)%lambda(nbspx,nbspx))
     ALLOCATE(rep_el(sm_k)%lambdam(nbspx,nbspx))
     ALLOCATE(rep_el(sm_k)%lambdap(nbspx,nbspx))
     ALLOCATE(rep_el(sm_k)%bec  (nkb,nbsp))
     ALLOCATE(rep_el(sm_k)%rhovan(nhm*(nhm+1)/2,nat,nspin))
  ENDDO
  !
  WRITE( stdout,*) " Allocation for W.F. s : successful " 

  CALL allocate_mainvar &
       ( ngw, ngwt, ngb, ngs, ngm, nr1, nr2, nr3, dfftp%nr1x, dfftp%nr1x, dfftp%npl, &
         nnrx, nnrsx, nat, nax, nsp, nspin, nbsp, nbspx, 0, nupdwn, nkb, gzero, 1,   &
         'gamma', smd = .TRUE. )
  !
  !
  CALL allocate_local_pseudo( ngs, nsp )
  !
  ALLOCATE(vkb(ngw,nkb))
  ALLOCATE(deeq(nhm,nhm,nat,nspin))
  ALLOCATE(dbec (nkb,nbsp,3,3))
  ALLOCATE(drhog(ngm,nspin,3,3))
  ALLOCATE(drhor(nnrx,nspin,3,3))
  ALLOCATE(drhovan(nhm*(nhm+1)/2,nat,nspin,3,3))
  ALLOCATE( mat_z( 1, 1, 1 ) )
  !
  WRITE( stdout,*) ' Allocation for CP core : successful '
  !
  ALLOCATE(etot_ar(0:smd_p))
  ALLOCATE(ekin_ar(0:smd_p))
  ALLOCATE(eht_ar(0:smd_p))
  ALLOCATE(epseu_ar(0:smd_p))
  ALLOCATE(exc_ar(0:smd_p))
  ALLOCATE(esr_ar(0:smd_p))
  !
  ALLOCATE(deviation(smpm))
  ALLOCATE(maxforce(smpm))
  ALLOCATE(arc_now(0:smd_p))
  ALLOCATE(arc_tot(0:smd_p))
  ALLOCATE(arc_pre(0:smd_p))
  ALLOCATE(paraforce(smpm))
  ALLOCATE(err_const(smd_p))
  ALLOCATE(pvvcheck(smpm))
  !
  WRITE( stdout,*) ' Allocation for SM variables : successful '
  !
  ALLOCATE(p_tau0(0:smd_p))
  ALLOCATE(p_taup(0:smd_p))
  ALLOCATE(p_tan(0:smd_p))
  !
  WRITE( stdout,*) ' Allocation for pointers : successful '
  !
  !
  con_ite = 0
  deeq(:,:,:,:) = 0.d0
  deviation(:) = 0.d0
  maxforce(:) = 0.d0
  pre_ekinc(:) = 0.d0
  paraforce(:) = 0.d0
  err_const(:) = 0.d0
  ekinc(:) = 0.d0
  arc_now(:) = 0.d0
  arc_tot(:) = 0.d0
  arc_pre(:) = 0.d0
  !
  !
  WRITE( stdout,*) ' Allocation ended '
  !
  CALL flush_unit( stdout )
  !

666 CONTINUE
  !
  !
  DO sm_k=0,smd_p
     p_tau0(sm_k)%d3 => rep(sm_k)%tau0
     IF(tfor) THEN
        p_taup(sm_k)%d3 => rep(sm_k)%taup
     ELSE
        p_taup(sm_k)%d3 => rep(sm_k)%tau0
     ENDIF
     p_tan(sm_k)%d3 => rep(sm_k)%tan
  ENDDO
  !
  !
  temp1=tempw+tolp
  temp2=tempw-tolp
  gkbt = DBLE( ndega ) * tempw / factem
  kbt = tempw / factem

  etot_ar(0   ) = smd_ene_ini
  etot_ar(smd_p) = smd_ene_fin

  DO sm_k=0,smd_p
     rep(sm_k)%tausm=rep(sm_k)%taus
     rep(sm_k)%tausp=0.
     rep(sm_k)%taum=rep(sm_k)%tau0
     rep(sm_k)%taup=0.
     rep(sm_k)%vels  = 0.0d0
     rep(sm_k)%velsm = 0.0d0
  ENDDO
  !
  velsp = 0.
  !
  hnew = h
  !
  DO sm_k=1,smpm
     rep_el(sm_k)%lambda = 0.d0
     rep_el(sm_k)%cm = (0.d0, 0.d0)
     rep_el(sm_k)%c0 = (0.d0, 0.d0)
  ENDDO
  !
  !     mass preconditioning: ema0bg(i) = ratio of emass(g=0) to emass(g)
  !     for g**2>emass_cutoff the electron mass ema0bg(g) rises quadratically
  !
  CALL emass_precond( ema0bg, ggp, ngw, tpiba2, emass_cutoff )
  !
  ! ... calculating tangent for force transformation.
  !
  IF(smd_lm) CALL TANGENT(p_tau0,p_tan)
  !
  !
  INI_REP_LOOP : DO sm_k = 1, smpm 
     !

     sm_file = smwout + sm_k
     sm_ndr  = ndr + sm_k
     !
     !
     IF ( nbeg < 0 ) THEN 

        !======================================================================
        !    nbeg = -1 
        !======================================================================

        ! IF( trdwfc ) THEN  ! add a new flag
        !   CALL readfile                                            &
        !         &     ( 0, sm_ndr,h,hold,nfi,rep_el(sm_k)%cm,rep_el(sm_k)%cm,rep(sm_k)%taus,  &
        !         &       rep(sm_k)%tausm,rep(sm_k)%vels,rep(sm_k)%velsm,rep(sm_k)%acc,         &
        !         &       rep_el(sm_k)%lambda,rep_el(sm_k)%lambdam,                             &
        !         &       xnhe0(sm_k),xnhem(sm_k),vnhe(sm_k),xnhp0(:,sm_k),xnhpm(:,sm_k),vnhp(:,sm_k),&
        !         &       nhpcl, ekincm(sm_k),                           &
        !         &       xnhh0,xnhhm,vnhh,velh,ecutp,ecutw,delt,pmass,ibrav,celldm,rep(sm_k)%fion, &
        !         &       tps, mat_z, f )
        ! ENDIF
        !     

        CALL phfac( rep(sm_k)%tau0, ei1, ei2, ei3, eigr )
        !
        CALL initbox ( rep(sm_k)%tau0, taub, irb )
        !
        CALL phbox( taub, eigrb )
        !
        IF(trane) THEN    
           IF(sm_k == 1) THEN
              !       
              !     random initialization
              !
              CALL randin(1,nbsp,gstart,ngw,ampre,rep_el(sm_k)%cm)
           ELSE
              rep_el(sm_k)%cm = rep_el(1)%cm
           ENDIF
        ELSE 
           ! IF(sm_k == 1) THEN   ! To be checked
           !    !       
           !    !     gaussian initialization
           !    !
           !    CALL gausin(eigr,rep_el(sm_k)%cm)
           ! ELSE
           !    rep_el(sm_k)%cm = rep_el(1)%cm 
           ! ENDIF
        END IF           
        !
        !     prefor calculates vkb (used by gram)
        !
        CALL prefor(eigr,vkb)
        !
        CALL gram( vkb, rep_el(sm_k)%bec, nkb, rep_el(sm_k)%cm, ngw, nbsp )
        !
        IF(iprsta.GE.3) CALL dotcsc(eigr,rep_el(sm_k)%cm)
        !     
        nfi=0
        !
        !     strucf calculates the structure factor sfac
        !
        CALL strucf( sfac, ei1, ei2, ei3, mill_l, ngs )

        IF(sm_k == 1) CALL formf(tfirst,eself) 

        CALL calbec (1,nsp,eigr,rep_el(sm_k)%cm,rep_el(sm_k)%bec)
        IF (tpre) CALL caldbec( ngw, nkb, nbsp, 1,nsp,eigr,rep_el(sm_k)%cm,dbec,.true.)
        !
        !
        WRITE(stdout,*) " "
        WRITE(stdout,*) " --------------------- Replica : ", sm_k
        WRITE(stdout,*) " "

        !

        CALL rhoofr (nfi,rep_el(sm_k)%cm,irb,eigrb,rep_el(sm_k)%bec, & 
             & rep_el(sm_k)%rhovan,rhor,rhog,rhos,enl,ekin_ar(sm_k))

        ekin = ekin_ar(sm_k)

        !
        IF(iprsta.GT.0) WRITE( stdout,*) ' out from rhoofr'
        !
        !     put core charge (if present) in rhoc(r)
        !
        IF ( ANY( nlcc ) ) CALL set_cc(irb,eigrb,rhoc)
        !
        !
        CALL vofrho(nfi,rhor,rhog,rhos,rhoc,tfirst,tlast,             &
             &        ei1,ei2,ei3,irb,eigrb,sfac,rep(sm_k)%tau0,rep(sm_k)%fion)

        IF( ionode .AND. &
            ( ( nfi == 0 ) .or. ( MOD( nfi-1, iprint ) == 0 ) .or. tfirst .or. tlast ) ) &
            CALL print_energies( .FALSE. )

        !
        etot_ar(sm_k) = etot
        eht_ar(sm_k) = eht
        epseu_ar(sm_k) = epseu
        exc_ar(sm_k) = exc
        !
        !
        CALL compute_stress( stress, detot, h, omega )
        !
        !
        IF(iprsta.GT.0 .AND. ionode) WRITE( sm_file,*) ' out from vofrho'
        IF(iprsta.GT.2) CALL print_atomic_var( rep(sm_k)%fion, na, nsp, ' fion ', iunit = sm_file )
        ! 
        !     forces for eigenfunctions
        !
        !     newd calculates deeq and a contribution to fion
        !
        CALL newd(rhor,irb,eigrb,rep_el(sm_k)%rhovan,rep(sm_k)%fion)

        CALL prefor(eigr,vkb)
        !
        !     if nbsp is odd => c(*,nbsp+1)=0
        !
        CALL runcp_uspp( nfi, fccc(sm_k), ccc(sm_k), ema0bg, dt2bye, rhos, &
             rep_el(sm_k)%bec, rep_el(sm_k)%cm, rep_el(sm_k)%c0, fromscra = .TRUE. )
        !
        !     buffer for wavefunctions is unit 21
        !
        IF(tbuff) REWIND 21
        !
        !     nlfq needs deeq calculated in newd
        !
        IF ( tfor .OR. tprnfor ) CALL nlfq(rep_el(sm_k)%cm,eigr, &
             & rep_el(sm_k)%bec,becdr,rep(sm_k)%fion)
        !
        !     imposing the orthogonality
        !     ==========================================================
        !
        CALL calphi(rep_el(sm_k)%cm,ema0bg,rep_el(sm_k)%bec, &
             & vkb,rep_el(sm_k)%phi)
        !
        !
        IF(ionode) WRITE( sm_file,*) ' out from calphi'
        !     ==========================================================
        !
        IF(tortho) THEN
           CALL ortho  (eigr,rep_el(sm_k)%c0,rep_el(sm_k)%phi,rep_el(sm_k)%lambda, &
                &                   bigr,iter,ccc(sm_k),ortho_eps,ortho_max,delt,bephi,becp)
        ELSE
           CALL gram( vkb, rep_el(sm_k)%bec, nkb, rep_el(sm_k)%c0, ngw, nbsp )
           !
           IF(ionode) WRITE( sm_file,*) ' gram  c0 '
        ENDIF
        !
        !     nlfl needs lambda becdr and bec
        !
        IF ( tfor .OR. tprnfor ) CALL nlfl(rep_el(sm_k)%bec,becdr, &
             & rep_el(sm_k)%lambda,rep(sm_k)%fion)
        !
        IF((tfor .OR. tprnfor) .AND. ionode) WRITE( sm_file,*) ' out from nlfl'
        !
        IF(iprsta.GE.3) CALL print_lambda( rep_el(sm_k)%lambda, nbsp, 9, ccc(sm_k), iunit = sm_file )
        !
        IF(tpre) THEN
           CALL nlfh(rep_el(sm_k)%bec,dbec,rep_el(sm_k)%lambda)
           WRITE( stdout,*) ' out from nlfh'
        ENDIF
        !
        IF(tortho) THEN
           CALL updatc(ccc(sm_k),rep_el(sm_k)%lambda,rep_el(sm_k)%phi, &
                & bephi,becp,rep_el(sm_k)%bec,rep_el(sm_k)%c0)
           !
           IF(ionode) WRITE( sm_file,*) ' out from updatc'
        ENDIF
        CALL calbec (nvb+1,nsp,eigr,rep_el(sm_k)%c0,rep_el(sm_k)%bec)
        IF (tpre) CALL caldbec(ngw,nkb,nbsp,1,nsp,eigr,rep_el(sm_k)%cm,dbec,.true.)
        !
        IF(ionode) WRITE( sm_file,*) ' out from calbec'
        !
        !     ==============================================================
        !     cm now orthogonalized
        !     ==============================================================
        IF(iprsta.GE.3) CALL dotcsc(eigr,rep_el(sm_k)%c0)
        !     
        IF(thdyn) THEN
           CALL cell_force( fcell, ainv, stress, omega, press, wmass )
           CALL cell_hmove( h, hold, delt, iforceh, fcell )
           CALL invmat( 3, h, ainv, deth )
        ENDIF
        !
        IF( tfor ) THEN

           IF(smd_lm) CALL PERP(rep(sm_k)%fion,rep(sm_k)%tan,paraforce(sm_k)) 

           CALL ions_hmove( rep(sm_k)%taus, rep(sm_k)%tausm, iforce, pmass, rep(sm_k)%fion, ainv, delt, na, nsp )
           CALL s_to_r( rep(sm_k)%taus, rep(sm_k)%tau0, na, nsp, h )

        ENDIF
        !
        !     
     ELSE   
        !
        !======================================================================
        !       nbeg = 0, nbeg = 1
        !======================================================================

        CALL readfile                                           &
             &     ( 1, sm_ndr,h,hold,nfi,rep_el(sm_k)%c0,rep_el(sm_k)%cm,rep(sm_k)%taus, &
             &       rep(sm_k)%tausm,rep(sm_k)%vels,rep(sm_k)%velsm,rep(sm_k)%acc,         &
             &       rep_el(sm_k)%lambda,rep_el(sm_k)%lambdam,                   &
             &       xnhe0(sm_k),xnhem(sm_k),vnhe(sm_k),xnhp0(:,sm_k),xnhpm(:,sm_k),vnhp(:,sm_k),&
             &       nhpcl, nhpdim, ekincm(sm_k),                                                         &
             &       xnhh0,xnhhm,vnhh,velh,ecutp,ecutw,delt,pmass,ibrav,celldm,rep(sm_k)%fion, &
             &       tps, mat_z, f )


        CALL from_restart( tfirst, rep(sm_k)%taus, rep(sm_k)%tau0, h, eigr, &
             rep_el(sm_k)%bec, rep_el(sm_k)%c0, rep_el(sm_k)%cm, ei1, ei2, ei3, sfac, eself )
        !
        !
     END IF


     !=============== END of NBEG selection ================================

     !
     !
     !     =================================================================
     !     restart with new averages and nfi=0
     !     =================================================================
     IF( nbeg <= 0 ) THEN
        rep(sm_k)%acc = 0.0d0
        nfi=0
     END IF
     !
     IF( ( .NOT. tfor ) .AND. ( .NOT. tprnfor ) ) THEN
        rep(sm_k)%fion = 0.d0
     END IF
     !
     IF( .NOT. tpre ) THEN
        stress = 0.d0
     ENDIF
     !         
     fccc = 1.0d0
     !
     IF(sm_k==1) nomore = nomore + nfi
     !
     !
     !      call cofmass( rep(sm_k)%taus, rep(sm_k)%cdm0 )
     !
     !
  ENDDO  INI_REP_LOOP      ! <<<<<<<<<<<<<<<<<<<<<<< ! 
  !
  !
  !
  IF(smd_lm .AND. nbeg < 0) THEN

     !  .... temp assignment ...

     DO sm_k=0,smd_p
        NULLIFY(p_taup(sm_k)%d3)
        p_taup(sm_k)%d3 => rep(sm_k)%taum
     ENDDO

     !
     !    ... calculate LM for smd .
     !

     CALL SMLAMBDA(p_tau0,p_taup,p_tan,con_ite,err_const)

     IF(smd_maxlm <= con_ite ) WRITE(stdout, *) "Warning ! : ", smd_maxlm, con_ite 
     !
     !    ... to reduced (crystal) coordinates.
     ! 

     DO sm_k =1,smpm
        CALL r_to_s(rep(sm_k)%tau0,rep(sm_k)%taus, na, nsp, ainv)
     ENDDO

     !
     !  .... back to regular assignemnt
     !
     DO sm_k=0,smd_p
        NULLIFY(p_taup(sm_k)%d3)
        p_taup(sm_k)%d3 => rep(sm_k)%taup
     ENDDO

     !
  ENDIF
  !
  !
  INI2_REP_LOOP : DO sm_k=1,smpm   ! >>>>>>>>>>>>>>>>>>>>>> !

     IF(tfor .OR. smd_lm) THEN  
        CALL phfac(rep(sm_k)%tau0,ei1,ei2,ei3,eigr)
        CALL calbec (1,nsp,eigr,rep_el(sm_k)%c0,rep_el(sm_k)%bec)
        IF (tpre) CALL caldbec(ngw,nkb,nbsp,1,nsp,eigr,rep_el(sm_k)%c0,dbec,.true.)
     ENDIF
     !
     xnhp0(:,sm_k)=0.
     xnhpm(:,sm_k)=0.
     vnhp(:,sm_k) =0.
     rep(sm_k)%fionm=0.d0
     CALL ions_vel(  rep(sm_k)%vels,  rep(sm_k)%taus,  rep(sm_k)%tausm, na, nsp, delt )

     velh = ( h - hold ) / delt
     !
     !     ======================================================
     !     kinetic energy of the electrons
     !     ======================================================
     !
     ekincm(sm_k)=0.0
     CALL elec_fakekine2( ekincm(sm_k), ema0bg, emass, rep_el(sm_k)%c0, rep_el(sm_k)%cm, ngw, nbsp, delt )

     xnhe0(sm_k)=0.
     xnhem(sm_k)=0.
     vnhe(sm_k) =0.
     !
     rep_el(sm_k)%lambdam(:,:)=rep_el(sm_k)%lambda(:,:)
     !
  ENDDO INI2_REP_LOOP  ! <<<<<<<<<<<<<<<<<<<  !

  !
  !
  ! ... Copying the init and final state.
  !
  rep(0)%taup = rep(0)%tau0
  rep(smd_p)%taup = rep(smd_p)%tau0

  ! ... Center of mass
  !
  WRITE(stdout,*) " "
  WRITE(stdout,*) " Center of mass "
  WRITE(stdout,*) " "
  DO sm_k=0,smd_p
     CALL ions_cofmass( rep(sm_k)%tau0, pmass, na, nsp, rep(sm_k)%cdm0)
     WRITE(stdout,'(i4,1x,3f8.5)') sm_k, (rep(sm_k)%cdm0(i),i=1,3)
  ENDDO
  WRITE(stdout,*) " "
  !
  ! ... Initial geometry ..
  !
  IF(ionode) THEN
     IF(nbeg < 0) THEN
        WRITE(unico,'(9(1x,f9.5))') &
             & (((rep(j)%taum(i,ia),i=1,3),ia=1,SUM(na(1:nsp))),j=0,smd_p)
     ELSE
        WRITE(unico,'(9(1x,f9.5))') &
             & (((rep(j)%tau0(i,ia),i=1,3),ia=1,SUM(na(1:nsp))),j=0,smd_p)
     ENDIF
     !
     CALL flush_unit( stdout )
     !
  ENDIF
  !
  !   basic loop for molecular dynamics starts here
  !
  WRITE(stdout,*) " _____________________________________________"
  WRITE(stdout,*) " "
  WRITE(stdout,*) " *****    Entering SMD LOOP   ***** "
  WRITE(stdout,*) " _____________________________________________"
  !
  !
  CALL stop_clock( 'initialize' )

  MAIN_LOOP: DO

     !
     !
     EL_REP_LOOP: DO sm_k=1,smpm        ! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !
        !
        sm_file = smwout + sm_k
        sm_ndr = ndr + sm_k
        !
        !     calculation of velocity of nose-hoover variables
        !
        IF(.NOT.tsde) fccc(sm_k)=1./(1.+frice)
        !
        CALL initbox ( rep(sm_k)%tau0, taub, irb )
        CALL phbox(taub,eigrb)
        CALL phfac( rep(sm_k)%tau0,ei1,ei2,ei3,eigr) 
        !
        !     strucf calculates the structure factor sfac
        !
        CALL strucf( sfac, ei1, ei2, ei3, mill_l, ngs )
        IF (thdyn) CALL formf(tfirst,eself)
        !
        IF(sm_k==1) THEN
           nfi=nfi+1
           tlast=(nfi.EQ.nomore)
        ENDIF
        !
        !
        IF((nfi.EQ.0) .OR. (MOD(nfi-1,iprint).EQ.0)) THEN
           WRITE(stdout,*) " "
           WRITE(stdout,*) " -----------  REPLICA : ", sm_k
           WRITE(stdout,*) " "
        ENDIF
        !
        !

        CALL rhoofr (nfi,rep_el(sm_k)%c0,irb,eigrb,rep_el(sm_k)%bec, &
             & rep_el(sm_k)%rhovan,rhor,rhog,rhos,enl,ekin_ar(sm_k))

        ekin = ekin_ar(sm_k)


        IF(trhow) THEN
           WRITE(stdout,*) " TRHOW set to .FASLE., it is not implemented with SMD" 
           trhow = .FALSE.   
        ENDIF
        !
        ! Y.K.
        !#ifdef __PARA     
        !      if(trhow .and. tlast) call write_rho(47,nspin,rhor)
        !#else
        !      if(trhow .and. tlast) write(47) ((rhor(i,is),i=1,nnrx),is=1,nspin)
        !#endif
        !
        !     put core charge (if present) in rhoc(r)
        !
        IF ( ANY( nlcc ) ) CALL set_cc(irb,eigrb,rhoc)
        !
        !
        IF(tlast) THEN
           WRITE(stdout,*) " "
           WRITE(stdout,*) " -----------  REPLICA : ", sm_k
           WRITE(stdout,*) " "
        ENDIF
        !
        IF(.NOT. tfirst) esr = esr_ar(sm_k)
        !

        CALL vofrho(nfi,rhor,rhog,rhos,rhoc,tfirst,tlast,                 &
                    ei1,ei2,ei3,irb,eigrb,sfac,rep(sm_k)%tau0,rep(sm_k)%fion)

        IF( ionode .AND. &
            ( ( nfi == 0 ) .or. ( MOD( nfi-1, iprint ) == 0 ) .or. tfirst .or. tlast ) ) &
            CALL print_energies( .FALSE. )
        !
        etot_ar(sm_k)  = etot
        eht_ar(sm_k)   = eht
        epseu_ar(sm_k) = epseu   
        exc_ar(sm_k)   = exc   
        !
        IF(tfirst) esr_ar(sm_k) = esr
        !
        !
        CALL compute_stress( stress, detot, h, omega )
        !
        enthal(sm_k)=etot+press*omega
        !
        !
        CALL newd(rhor,irb,eigrb,rep_el(sm_k)%rhovan,rep(sm_k)%fion)
        !
        CALL prefor(eigr,vkb)

        CALL runcp_uspp( nfi, fccc(sm_k), ccc(sm_k), ema0bg, dt2bye, rhos, &
             rep_el(sm_k)%bec, rep_el(sm_k)%c0, rep_el(sm_k)%cm )
        !
        !     buffer for wavefunctions is unit 21
        !
        IF(tbuff) REWIND 21
        !
        !----------------------------------------------------------------------
        !                 contribution to fion due to lambda
        !----------------------------------------------------------------------
        !
        !     nlfq needs deeq bec
        !
        IF ( tfor .OR. tprnfor ) CALL nlfq(rep_el(sm_k)%c0,eigr, &
             & rep_el(sm_k)%bec,becdr,rep(sm_k)%fion)
        !
        IF( tfor .OR. thdyn ) THEN
           !
           ! interpolate new lambda at (t+dt) from lambda(t) and lambda(t-dt):
           !
           rep_el(sm_k)%lambdap(:,:) = 2.d0*rep_el(sm_k)%lambda(:,:)-rep_el(sm_k)%lambdam(:,:)
           rep_el(sm_k)%lambdam(:,:)= rep_el(sm_k)%lambda (:,:)
           rep_el(sm_k)%lambda (:,:)= rep_el(sm_k)%lambdap(:,:)
        ENDIF
        !
        !     calphi calculates phi
        !     the electron mass rises with g**2
        !
        CALL calphi(rep_el(sm_k)%c0,ema0bg,rep_el(sm_k)%bec,vkb,rep_el(sm_k)%phi)
        !
        !     begin try and error loop (only one step!)
        !
        !       nlfl and nlfh need: lambda (guessed) becdr
        !
        IF ( tfor .OR. tprnfor ) CALL nlfl(rep_el(sm_k)%bec,becdr,rep_el(sm_k)%lambda, &
             & rep(sm_k)%fion)
        !
        !
        !     This part is not compatible with SMD.
        !
        IF(tpre) THEN
           CALL nlfh(rep_el(sm_k)%bec,dbec,rep_el(sm_k)%lambda)
           CALL ions_thermal_stress( thstress, pmass, omega, h, rep(sm_k)%vels, nsp, na )
           stress = stress + thstress
        ENDIF

        !
     ENDDO EL_REP_LOOP   ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !


     !_________________________________________________________________________!
     !
     !
     IF(smd_lm) THEN ! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !
        !
        !
        ! ... Transforming ionic forces to perp force. 
        !

        CALL TANGENT(p_tau0,p_tan)

        DO sm_k=1,smpm
           CALL PERP(rep(sm_k)%fion,rep(sm_k)%tan,paraforce(sm_k))
        ENDDO


        DO sm_k=1,smpm 

           deviation(sm_k) = 0.d0
           isa = 0
           DO is=1,nsp
              DO ia=1,na(is)
                 isa = isa + 1
                 DO i=1,3

                    deviation(sm_k) = deviation(sm_k) + (rep(sm_k)%fion(i,isa)* &
                         &    iforce(i,isa))**2.d0   

                    workvec(i,ia,is) = rep(sm_k)%fion(i,isa)*iforce(i,isa)

                 ENDDO
              ENDDO
           ENDDO

           deviation(sm_k) = DSQRT(deviation(sm_k))
           maxforce(sm_k) = MAX(ABS(MAXVAL(workvec)),ABS(MINVAL(workvec)))
        ENDDO
        !
     ENDIF ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<!
     !
     !
     !  ... Writing main output ...
     !
     !
     IF(tfirst) THEN
        WRITE(stdout,*) " " 
        WRITE(stdout,1997)  
        WRITE(stdout,*) " " 
     ENDIF
     !
     IF(smd_lm) THEN
        WRITE(stdout,1998) nfi,SUM(etot_ar(0:smd_p)),con_ite,MAXVAL(err_const), &
             & MAX(MAXLOC(err_const),0),MAXVAL(ekinc(1:smpm)), &
             & MAXVAL(ekinc(1:smpm)-pre_ekinc(1:smpm)), & 
             & SUM(deviation)/(smd_p-1),MAXVAL(maxforce)
     ENDIF
     !
1997 FORMAT(/2x,'nfi',3x,'SUM(E)',3x,'con_ite',5x,'con_usat',5x, &
          &     'cons_pls',5x,'max_ekinc',3x,   &
          &     'deviation',3x,'maxforce')
     !
1998 FORMAT(i5,1x,F14.5,1x,i5,1x,E12.5, &
          & 1x,i3,1x,f8.5, &
          & 1x,f8.5, &
          & 1x,E12.5,1x,E12.5)
     !
     !
     !________________________________________________________________________!
     !
     !=======================================================================
     !
     !              verlet algorithm
     !
     !     loop which updates cell parameters and ionic degrees of freedom
     !     hnew=h(t+dt) is obtained from hold=h(t-dt) and h=h(t)
     !     tausp=pos(t+dt) from tausm=pos(t-dt) taus=pos(t) h=h(t)
     !
     !           guessed displacement of ions
     !=======================================================================
     !
     !     CELL_dynamics
     !
     hgamma(:,:) = 0.d0
     IF(thdyn) THEN

        CALL cell_force( fcell, ainv, stress, omega, press, wmass )

        CALL cell_move( hnew, h, hold, delt, iforceh, fcell, frich, tnoseh, vnhh, velh, tsdc )
        !
        velh(:,:) = ( hnew(:,:) - hold(:,:) ) / twodel
        !
        CALL cell_gamma( hgamma, ainv, h, velh )
        !
     ENDIF
     !
     !======================================================================
     !
     TFOR_IF : IF( tfor ) THEN

        ION_REP_LOOP : DO sm_k=1,smpm

           CALL ions_move( rep(sm_k)%tausp, rep(sm_k)%taus, rep(sm_k)%tausm, iforce, pmass, &
                rep(sm_k)%fion, ainv, delt, na, nsp, fricp, hgamma, rep(sm_k)%vels, tsdp, &
                tnosep, rep(sm_k)%fionm, vnhp(:,sm_k), velsp, rep(sm_k)%velsm, 1, 1, atm2nhp )
           !
           !cc   call cofmass(velsp,rep(sm_k)%cdmvel)
           !         call cofmass(rep(sm_k)%tausp,cdm)
           !         do is=1,nsp
           !            do ia=1,na(is)
           !               do i=1,3
           !cc   velsp(i,ia,is)=velsp(i,ia,is)-cdmvel(i)
           !                  tausp(i,ia,is)=tausp(i,ia,is) ! +cdm0(i)-cdm(i)
           !               enddo
           !            enddo
           !         enddo
           !
           !  ... taup is obtained from tausp ...
           !
           CALL  s_to_r( rep(sm_k)%tausp, rep(sm_k)%taup, na, nsp, hnew )

        ENDDO ION_REP_LOOP

     ENDIF TFOR_IF
     !
     !---------------------------------------------------------------------------
     !              String method const  .. done in real coordiantes ...
     !---------------------------------------------------------------------------
     !
     CALL ARC(p_taup,arc_pre,t_arc_pre,1)
     !
     IF( MOD( nfi, smd_lmfreq ) == 0 ) THEN
        IF(smd_lm) THEN
           !
           CALL SMLAMBDA( p_taup, p_tau0, p_tan, con_ite, err_const )
           !
        ENDIF
     ENDIF
     !
     CALL ARC( p_taup, arc_now, t_arc_now, 1 )
     CALL ARC( p_taup, arc_tot, t_arc_tot, 0 )
     !
     !     ... move back to reduced coordiinates
     !     
     DO sm_k=1,smpm 
        CALL r_to_s(rep(sm_k)%taup,rep(sm_k)%tausp, na, nsp, ainv)
     ENDDO
     ! 
     !    
     POST_REP_LOOP : DO sm_k = 1, smpm 
        !
        sm_file  =  smwout + sm_k
        sm_ndw = ndw + sm_k
        !
        !---------------------------------------------------------------------------
        !              initialization with guessed positions of ions
        !---------------------------------------------------------------------------
        !
        !  if thdyn=true g vectors and pseudopotentials are recalculated for 
        !  the new cell parameters
        !
        IF ( tfor .OR. thdyn ) THEN
           IF( thdyn ) THEN
              hold = h
              h = hnew
              CALL newinit( h )
              CALL newnlinit
           ELSE
              hold = h
           ENDIF
           !        ... phfac calculates eigr
           !
           CALL phfac(rep(sm_k)%taup,ei1,ei2,ei3,eigr)
           !
        ELSE 
           !
           CALL phfac(rep(sm_k)%tau0,ei1,ei2,ei3,eigr)
           !
        END IF
        !
        !        ... prefor calculates vkb
        !
        CALL prefor(eigr,vkb)
        !
        !---------------------------------------------------------------------------
        !                    imposing the orthogonality
        !---------------------------------------------------------------------------
        !
        IF(tortho) THEN
           CALL ortho                                                     &
                &         (eigr,rep_el(sm_k)%cm,rep_el(sm_k)%phi,rep_el(sm_k)%lambda, &
                & bigr,iter,ccc(sm_k),ortho_eps,ortho_max,delt,bephi,becp)
        ELSE
           CALL gram( vkb, rep_el(sm_k)%bec, nkb, rep_el(sm_k)%cm, ngw, nbsp )
           IF(iprsta.GT.4) CALL dotcsc(eigr,rep_el(sm_k)%cm)
        ENDIF
        !
        !---------------------------------------------------------------------------
        !                   correction to displacement of ions
        !---------------------------------------------------------------------------
        !
        IF(iprsta.GE.3) CALL print_lambda( rep_el(sm_k)%lambda, nbsp, 9, 1.0d0 )
        !
        IF(tortho) CALL updatc(ccc(sm_k),rep_el(sm_k)%lambda,rep_el(sm_k)%phi,bephi, &
             & becp,rep_el(sm_k)%bec,rep_el(sm_k)%cm)
        !
        CALL calbec (nvb+1,nsp,eigr,rep_el(sm_k)%cm,rep_el(sm_k)%bec)
        IF (tpre) CALL caldbec(ngw,nkb,nbsp,1,nsp,eigr,rep_el(sm_k)%cm,dbec,.true.)
        !
        IF(iprsta.GE.3)  CALL dotcsc(eigr,rep_el(sm_k)%cm)
        !
        !---------------------------------------------------------------------------
        !                  temperature monitored and controlled
        !---------------------------------------------------------------------------
        !
        ekinp(sm_k)  = 0.d0
        ekinpr(sm_k) = 0.d0
        !
        !     ionic kinetic energy 
        !
        IF( tfor ) THEN
           CALL ions_vel( rep(sm_k)%vels, rep(sm_k)%tausp, rep(sm_k)%tausm, na, nsp, delt )
           CALL ions_kinene( ekinp(sm_k), rep(sm_k)%vels, na, nsp, hold, pmass )
        ENDIF
        !
        !     ionic temperature
        !
        IF( tfor ) THEN
           CALL ions_temp( tempp(sm_k), temps, ekinpr(sm_k), rep(sm_k)%vels, &
                na, nsp, hold, pmass, ndega, nhpdim, atm2nhp,     &
                ekin2nhp )
        ENDIF
        !
        !     fake electronic kinetic energy
        !
        CALL elec_fakekine2( ekinc0(sm_k), ema0bg, emass, rep_el(sm_k)%c0, rep_el(sm_k)%cm, ngw, nbsp, delt )
        !
        !     ... previous ekinc
        !
        pre_ekinc(sm_k) = ekinc(sm_k)

        ekinc(sm_k) = 0.5 * ( ekinc0(sm_k) + ekincm(sm_k) )
        !
        !     fake cell-parameters kinetic energy
        !
        ekinh=0.
        IF(thdyn) THEN
           CALL cell_kinene( ekinh, temphh, velh )
        ENDIF
        IF(thdiag) THEN
           temphc=2.*factem*ekinh/3.
        ELSE
           temphc=2.*factem*ekinh/9.
        ENDIF
        !
        ! warning! thdyn and tcp/tcap are not compatible yet!!!
        !
        IF(tcp.OR.tcap.AND.tfor.AND.(.NOT.thdyn)) THEN
           IF(tempp(sm_k).GT.temp1.OR.tempp(sm_k).LT.temp2.AND.tempp(sm_k).NE.0.d0) THEN
              CALL  ions_vrescal( tcap, tempw, tempp(sm_k), rep(sm_k)%taup, rep(sm_k)%tau0, rep(sm_k)%taum, &
                   na, nsp, rep(sm_k)%fion, iforce, pmass, delt )
           END IF
        END IF
        !
        ! ------------------------------
        !
        IF(MOD(nfi-1,iprint).EQ.0 .OR. tlast) THEN
           WRITE(stdout,*) " "
           WRITE(stdout,*) " ---------------------- Replica : ", sm_k
           WRITE(stdout,*) " "
        ENDIF
        !
        !
        IF(MOD(nfi-1,iprint).EQ.0 .OR. (nfi.EQ.(nomore))) THEN
           CALL eigs0(.true.,nspin,nbspx,nupdwn,iupdwn,f,rep_el(sm_k)%lambda)
           WRITE( stdout,*)
        ENDIF
        !
        epot(sm_k)=eht_ar(sm_k)+epseu_ar(sm_k)+exc_ar(sm_k)
        !
        rep(sm_k)%acc(1)=rep(sm_k)%acc(1)+ekinc(sm_k)
        rep(sm_k)%acc(2)=rep(sm_k)%acc(2)+ekin_ar(sm_k)
        rep(sm_k)%acc(3)=rep(sm_k)%acc(3)+epot(sm_k)
        rep(sm_k)%acc(4)=rep(sm_k)%acc(4)+etot_ar(sm_k)
        rep(sm_k)%acc(5)=rep(sm_k)%acc(5)+tempp(sm_k)
        !
        econs(sm_k)=ekinp(sm_k)+ekinh+enthal(sm_k)
        econt(sm_k)=econs(sm_k)+ekinc(sm_k)
        !
        IF(tnosee)THEN
           econt(sm_k)=econt(sm_k)+ electrons_nose_nrg( xnhe0(sm_k), vnhe(sm_k), qne, ekincw )
        ENDIF
        !
        !     ... Writing the smfiles ...
        !
        IF(MOD(nfi-1,iprint).EQ.0.OR.tfirst)  THEN
           IF(ionode) WRITE( sm_file,*)
           IF(ionode) WRITE( sm_file,1949)
        END IF
        !
        tps = nfi * delt * AU_PS
        !
        IF(ionode) WRITE( sm_file,1950) nfi, ekinc(sm_k), INT(tempp(sm_k)), &
             &              etot_ar(sm_k), econs(sm_k), econt(sm_k),              &
             &              arc_now(sm_k),t_arc_now,arc_pre(sm_k),arc_tot(sm_k),  &
             &              deviation(sm_k),maxforce(sm_k),paraforce(sm_k)

        ! Y.K.
        !      write(8,2948) tps,ekinc,temphc,tempp,enthal,econs,      &
        !     &              econt,                                              &
        !     &              vnhh(3,3),xnhh0(3,3),vnhp,xnhp0

        !
        !     c              frice,frich,fricp
        ! 
1949    FORMAT(/2x,'nfi',4x,'ekinc',1x,'tempp',8x,'etot',7x,'econs',7x,'econt',   &
             &     3x,'arc_diff',2x,'real',4x,'ori_arc',4x,'tot_arc',3x,'dev',4x,'maxF',4x,'paraF')
        !
        !cc     f       7x,'econs',7x,'econt',3x,'frice',2x,'frich',2x,'fricp')
        !
1950    FORMAT(i5,1x,f8.5,1x,i5,1x,f11.5,1x,f11.5,1x,f11.5, &
             & 1x,f8.5,1x,f8.5,1x,f8.5,1x,f8.5,1x,f8.5,1x,f8.5,1x,f8.5)
        !
#if defined (FLUSH) 
        CALL flush_unit( sm_file )
#endif
        !
        ! 
2948    FORMAT(f8.5,1x,f8.5,1x,f6.1,1x,f6.1,3(1x,f11.5),4(1x,f7.4))
        !
        IF( tfor ) THEN
           IF ( ionode ) THEN
              IF(tlast .OR. (MOD(nfi,smd_codfreq) == 0)) THEN
                 WRITE(unico,*) "== COORD ==  rep : ", sm_k
                 WRITE(unico,3340) ((h(i,j),i=1,3),j=1,3)
                 WRITE(unico,'(3f12.8)') ((rep(sm_k)%tau0(i,ia),i=1,3),ia=1,SUM(na(1:nsp)))
              ENDIF
              !
              IF(tlast .OR. (MOD(nfi,smd_forfreq) == 0)) THEN
                 WRITE(unifo,*) "== FORCE ==  rep : ", sm_k
                 WRITE(unifo,'(3f12.8)') ((rep(sm_k)%fion(i,ia),i=1,3),ia=1,SUM(na(1:nsp)))
              ENDIF
              !
              IF(tlast) THEN
                 WRITE(unist,3340) ((stress(i,j),i=1,3),j=1,3)
              ENDIF
              !
#if defined (FLUSH) 
              CALL flush_unit( unico )
#endif
3340          FORMAT(9(1x,f9.5))
           ENDIF
           !
           !     new variables for next step
           !
           rep(sm_k)%tausm(:,1:nat)=rep(sm_k)%taus(:,1:nat)
           rep(sm_k)%taus(:,1:nat)=rep(sm_k)%tausp(:,1:nat)
           rep(sm_k)%taum(:,1:nat)=rep(sm_k)%tau0(:,1:nat)
           rep(sm_k)%tau0(:,1:nat)=rep(sm_k)%taup(:,1:nat)
           rep(sm_k)%velsm(:,1:nat) = rep(sm_k)%vels(:,1:nat)
           rep(sm_k)%vels(:,1:nat)  = velsp(:,1:nat)
           IF(tnosep) THEN
              xnhpm(:,sm_k) = xnhp0(:,sm_k)
              xnhp0(:,sm_k) = xnhpp(:,sm_k)
           ENDIF
           IF(tnosee) THEN
              xnhem(sm_k) = xnhe0(sm_k)
              xnhe0(sm_k) = xnhep(sm_k)
           ENDIF
           IF(tnoseh) THEN
              xnhhm(:,:) = xnhh0(:,:)
              xnhh0(:,:) = xnhhp(:,:)
           ENDIF
        END IF
        !
        IF(thdyn)THEN
           CALL emass_precond( ema0bg, ggp, ngw, tpiba2, emass_cutoff )
        ENDIF
        !
        ekincm(sm_k)=ekinc0(sm_k)
        !  
        !     cm=c(t+dt) c0=c(t)
        !
        CALL DSWAP(2*ngw*nbsp,rep_el(sm_k)%c0,1,rep_el(sm_k)%cm,1)
        !
        !     now:  cm=c(t) c0=c(t+dt)
        !
        IF (tfirst) THEN
           epre(sm_k) = etot_ar(sm_k)
           enow(sm_k) = etot_ar(sm_k)
        ENDIF
        !
        IF(sm_k == smpm) tfirst=.FALSE.
        !
        !     write on file ndw each isave
        !
        IF( ( MOD( nfi, isave ) == 0 ) .AND. ( nfi < nomore ) ) THEN

           CALL writefile                                                    &
                &     ( sm_ndw,h,hold,nfi,rep_el(sm_k)%c0,rep_el(sm_k)%cm,rep(sm_k)%taus, &
                &       rep(sm_k)%tausm,rep(sm_k)%vels,rep(sm_k)%velsm,rep(sm_k)%acc,     &
                &       rep_el(sm_k)%lambda,rep_el(sm_k)%lambdam,xnhe0(sm_k),xnhem(sm_k), &
                &       vnhe(sm_k),xnhp0(:,sm_k),xnhpm(:,sm_k),vnhp(:,sm_k),nhpcl,nhpdim, ekincm(sm_k),       &
                &       xnhh0,xnhhm,vnhh,velh,ecutp,ecutw,delt,pmass,ibrav,celldm,         &
                &       rep(sm_k)%fion, tps, mat_z, f, rhor )

        ENDIF
        !
        !
        epre(sm_k) = enow(sm_k)
        enow(sm_k) = etot_ar(sm_k)
        !
        frice = frice * grease
        fricp = fricp * greasp
        frich = frich * greash

        !     =====================================================
        !
        delta_etot = ABS( epre(sm_k) - enow(sm_k) )
        tstop = check_stop_now()
        tconv = .FALSE.

        ! Y.K.
        !      IF( tconvthrs%active ) THEN
        !       tconv = ( delta_etot < tconvthrs%derho ) .AND. ( ekinc(sm_k) < tconvthrs%ekin )
        !      END IF

        !
        !
        IF(tconv) THEN
           WRITE(stdout,*) "TCONV set to .F. : not implemented with calculation = smd "
           tconv = .FALSE. 
        ENDIF
        !
        IF( tconv ) THEN
           IF( ionode ) THEN
              WRITE( stdout,fmt= &
                   "(/,3X,'MAIN:',10X,'EKINC   (thr)',10X,'DETOT   (thr)',7X,'MAXFORCE   (thr)')" )
              WRITE( stdout,fmt="(3X,'MAIN: ',3(D14.6,1X,D8.1))" ) &
                   ekinc, tconvthrs%ekin, delta_etot, tconvthrs%derho, 0.0d0, tconvthrs%force
              WRITE( stdout,fmt="(3X,'MAIN: convergence achieved for system relaxation')")
           END IF
        END IF


        tstop = tstop .OR. tconv

     ENDDO POST_REP_LOOP ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !



     IF( (nfi >= nomore) .OR. tstop ) EXIT MAIN_LOOP

  END DO MAIN_LOOP

  !
  !============================= END of main LOOP ======================
  !


  !
  !  Here copy relevant physical quantities into the output arrays/variables
  !  for 0 and smd_p replicas.
  !

  etot_out(0) = etot_ar(0)
  etot_out(smd_p) = etot_ar(smd_p)

  isa = 0
  DO is = 1, nsp
     DO ia = 1, na(is)
        isa = isa + 1
        ipos = ind_srt( isa )
        tau( 1, ipos , 0 ) = rep(0)%tau0( 1, isa )
        tau( 2, ipos , 0 ) = rep(0)%tau0( 2, isa )
        tau( 3, ipos , 0 ) = rep(0)%tau0( 3, isa )
        tau( 1, ipos , smd_p) = rep(smd_p)%tau0( 1, isa )
        tau( 2, ipos , smd_p) = rep(smd_p)%tau0( 2, isa )
        tau( 3, ipos , smd_p) = rep(smd_p)%tau0( 3, isa )
        fion_out( 1, ipos, 0 ) = rep(0)%fion( 1, isa )
        fion_out( 2, ipos, 0 ) = rep(0)%fion( 2, isa )
        fion_out( 3, ipos, 0 ) = rep(0)%fion( 3, isa )
        fion_out( 1, ipos, smd_p ) = rep(smd_p)%fion( 1, isa )
        fion_out( 2, ipos, smd_p ) = rep(smd_p)%fion( 2, isa )
        fion_out( 3, ipos, smd_p ) = rep(smd_p)%fion( 3, isa )
     END DO
  END DO


  FIN_REP_LOOP : DO sm_k=1,smpm     ! >>>>>>>>>>>>>>>>>>>>>>>>>> !
     !
     sm_file = smwout + sm_k
     sm_ndw = ndw + sm_k
     !
     ! 
     !  Here copy relevant physical quantities into the output arrays/variables
     !
     !
     etot_out(sm_k) = etot_ar(sm_k)
     !
     isa = 0
     DO is = 1, nsp
        DO ia = 1, na(is)
           isa = isa + 1
           ipos = ind_srt( isa )
           tau( 1, ipos , sm_k) = rep(sm_k)%tau0( 1, isa )
           tau( 2, ipos , sm_k) = rep(sm_k)%tau0( 2, isa )
           tau( 3, ipos , sm_k) = rep(sm_k)%tau0( 3, isa )
           fion_out( 1, ipos, sm_k ) = rep(sm_k)%fion( 1, isa )
           fion_out( 2, ipos, sm_k ) = rep(sm_k)%fion( 2, isa )
           fion_out( 3, ipos, sm_k ) = rep(sm_k)%fion( 3, isa )
        END DO
     END DO

     !  Calculate statistics

     anor=1.d0/DBLE(nfi)
     DO i=1,nacc
        rep(sm_k)%acc(i)=rep(sm_k)%acc(i)*anor
     END DO
     !
     !
     IF(ionode) WRITE( sm_file,1951)
1951 FORMAT(//'              averaged quantities :',/,                 &
          &       9x,'ekinc',10x,'ekin',10x,'epot',10x,'etot',5x,'tempp')
     IF(ionode) WRITE( sm_file,1952) (rep(sm_k)%acc(i),i=1,nacc)
1952 FORMAT(4f14.5,f10.1)
     !
     !
#if defined (FLUSH)
     CALL flush_unit(sm_file)
#endif
     !
     !
     IF( sm_k == smpm ) THEN
        CALL print_clock( 'initialize' )
        CALL print_clock( 'formf' )
        CALL print_clock( 'rhoofr' )
        CALL print_clock( 'vofrho' )
        CALL print_clock( 'dforce' )
        CALL print_clock( 'calphi' )
        CALL print_clock( 'ortho' )
        CALL print_clock( 'updatc' )
        CALL print_clock( 'gram' )
        CALL print_clock( 'newd' )
        CALL print_clock( 'calbec' )
        CALL print_clock( 'prefor' )
        CALL print_clock( 'strucf' )
        CALL print_clock( 'nlfl' )
        CALL print_clock( 'nlfq' )
        CALL print_clock( 'set_cc' )
        CALL print_clock( 'rhov' )
        CALL print_clock( 'nlsm1' )
        CALL print_clock( 'nlsm2' )
        CALL print_clock( 'forcecc' )
        CALL print_clock( 'fft' )
        CALL print_clock( 'ffts' )
        CALL print_clock( 'fftw' )
        CALL print_clock( 'fftb' )
        CALL print_clock( 'rsg' )
        CALL print_clock( 'reduce' )
     END IF
     !
     !

     CALL writefile                                                    &
          &     ( sm_ndw,h,hold,nfi,rep_el(sm_k)%c0,rep_el(sm_k)%cm,rep(sm_k)%taus, &
          &       rep(sm_k)%tausm,rep(sm_k)%vels,rep(sm_k)%velsm,rep(sm_k)%acc,     &
          &       rep_el(sm_k)%lambda,rep_el(sm_k)%lambdam,xnhe0(sm_k),xnhem(sm_k), &
          &       vnhe(sm_k),xnhp0(:,sm_k),xnhpm(:,sm_k),vnhp(:,sm_k),nhpcl,nhpdim,ekincm(sm_k),       &
          &       xnhh0,xnhhm,vnhh,velh,ecutp,ecutw,delt,pmass,ibrav,celldm,         &
          &       rep(sm_k)%fion, tps, mat_z, f, rhor )

     !
     !
     !
     IF(iprsta.GT.1) CALL print_lambda( rep_el(sm_k)%lambda, nbsp, nbsp, 1.0d0, iunit = sm_file )
     !
     IF( tfor .OR. tprnfor ) THEN
        IF(ionode) WRITE( sm_file,1970) ibrav, alat
        IF(ionode) WRITE( sm_file,1971)
        DO i=1,3
           IF(ionode) WRITE( sm_file,1972) (h(i,j),j=1,3)
        ENDDO
        IF(ionode) WRITE( sm_file,1973)
        isa = 0
        DO is=1,nsp
           DO ia=1,na(is)
              isa = isa + 1
              IF(ionode) WRITE( sm_file,1974) is,ia,(rep(sm_k)%tau0(i,isa),i=1,3),       &
                   &            ((ainv(j,1)*rep(sm_k)%fion(1,isa)+ainv(j,2)*rep(sm_k)%fion(2,isa)+    &
                   &              ainv(j,3)*rep(sm_k)%fion(3,isa)),j=1,3)
           END DO
        END DO
        IF(ionode) WRITE( sm_file,1975)
        isa = 0
        DO is=1,nsp
           DO ia=1,na(is)
              isa = isa + 1
              IF(ionode) WRITE( sm_file,1976) is,ia,(rep(sm_k)%taus(i,isa),i=1,3)
           END DO
        END DO
     ENDIF
     conv_elec = .TRUE.



1970 FORMAT(1x,'ibrav :',i4,'  alat : ',f10.4,/)
1971 FORMAT(1x,'lattice vectors',/)
1972 FORMAT(1x,3f10.4)
1973 FORMAT(/1x,'Cartesian coordinates (a.u.)              forces (redu. units)' &
          &       /1x,'species',' atom #', &
          &           '   x         y         z      ', &
          &           '   fx        fy        fz'/)
1974 FORMAT(1x,2i5,3f10.4,2x,3f10.4)
1975 FORMAT(/1x,'Scaled coordinates '/1x,'species',' atom #')
1976 FORMAT(1x,2i5,3f10.4)
     IF(ionode) WRITE( sm_file,1977) 

     !
600  CONTINUE
     !
  ENDDO FIN_REP_LOOP               ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<< !
  !
  !
  CALL memory
  !      
1977 FORMAT(5x,//'====================== end cprvan ',                 &
       &            '======================',//)

  IF( ALLOCATED( vkb ) ) DEALLOCATE( vkb )
  IF( ALLOCATED( deeq ) ) DEALLOCATE( deeq )
  IF( ALLOCATED( deviation )) DEALLOCATE( deviation )
  IF( ALLOCATED( maxforce )) DEALLOCATE( maxforce )
  IF( ALLOCATED( arc_now )) DEALLOCATE( arc_now )
  IF( ALLOCATED( arc_tot )) DEALLOCATE( arc_tot )
  IF( ALLOCATED( arc_pre )) DEALLOCATE( arc_pre )
  IF( ALLOCATED( paraforce )) DEALLOCATE( paraforce )
  IF( ALLOCATED( err_const )) DEALLOCATE( err_const )
  IF( ALLOCATED( pvvcheck )) DEALLOCATE( pvvcheck  )
  !
  IF( ALLOCATED( mat_z )) DEALLOCATE( mat_z  )

  DO sm_k=0,smd_p  
     IF( ASSOCIATED( p_tau0(sm_k)%d3 )) NULLIFY( p_tau0(sm_k)%d3 ) 
     IF( ASSOCIATED( p_taup(sm_k)%d3 )) NULLIFY( p_taup(sm_k)%d3 ) 
     IF( ASSOCIATED( p_tan(sm_k)%d3 )) NULLIFY( p_tan(sm_k)%d3 ) 
  ENDDO

  IF( ALLOCATED( p_tau0 )) DEALLOCATE( p_tau0  )
  IF( ALLOCATED( p_taup )) DEALLOCATE( p_taup )
  IF( ALLOCATED( p_tan )) DEALLOCATE( p_tan )

  CALL deallocate_modules_var()

  CALL deallocate_smd_rep()
  CALL deallocate_smd_ene()

  IF( ionode ) THEN
     CLOSE( 8 )
     CLOSE( unico )
     CLOSE( unifo )
     CLOSE( unist )
     DO sm_k=1,smpm
        sm_file =  smwout + sm_k
        CLOSE(sm_file)
     ENDDO
  END IF
  !
#endif
  !
  RETURN
END SUBROUTINE smdmain
