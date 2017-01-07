  !                                                                            
  ! Copyright (C) 2010-2016 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino 
  ! Copyright (C) 2007-2009 Jesse Noffsinger, Brad Malone, Feliciano Giustino  
  !                                                                            
  ! This file is distributed under the terms of the GNU General Public         
  ! License. See the file `LICENSE' in the root directory of the               
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .             
  !                                                                            
  ! Adapted from the code PH/phq_setup - Quantum-ESPRESSO group                
  !-----------------------------------------------------------------------
  SUBROUTINE epw_setup
  !-----------------------------------------------------------------------
  !!
  !! EPW setup.  
  !!
  !! @Note: RM - Nov 2014: Noncolinear case implemented
  !!
  USE kinds,         ONLY : DP
  USE ions_base,     ONLY : tau, nat, ntyp => nsp, ityp
  USE cell_base,     ONLY : at, bg  
  USE io_global,     ONLY : stdout, ionode, ionode_id
  USE io_files,      ONLY : tmp_dir
  USE ener,          ONLY : Ef
  USE klist,         ONLY : xk, lgauss, degauss, ngauss, nks, nelec, nkstot, ltetra
  USE lsda_mod,      ONLY : nspin, lsda, starting_magnetization
  USE scf,           ONLY : v, vrs, vltot, rho, rho_core, kedtau
  USE gvect,         ONLY : ngm
  USE symm_base,     ONLY : nsym, s, irt, t_rev, time_reversal, invs,&
                            sr
  USE uspp_param,    ONLY : upf
  USE spin_orb,      ONLY : domag
  USE constants,     ONLY : degspin, pi
  USE noncollin_module,     ONLY : noncolin, m_loc, angle1, angle2, ux
  USE wvfct,         ONLY : nbnd, et
  USE output,        ONLY : fildrho
  USE eqv,           ONLY : dmuxc
  USE nlcc_ph,       ONLY : drc
  USE uspp,          ONLY : nlcc_any
  USE control_ph,    ONLY : lgamma_gamma, search_sym, start_irr, &
                            last_irr, niter_ph, alpha_mix, all_done,  &
                            flmixdpot, reduce_io, u_from_file
  USE control_lr,    ONLY : lgamma, alpha_pv, nbnd_occ
  USE gamma_gamma,   ONLY : has_equivalent, asr, nasr, n_diff_sites, &
                            equiv_atoms, n_equiv_atoms, with_symmetry
  USE partial,       ONLY : comp_irr, atomo, nat_todo, all_comp, &
                            done_irr
  USE modes,         ONLY : u, npertx, npert, nirr, nmodes, num_rap_mode
  USE lr_symm_base,  ONLY : gi, gimq, irotmq, minus_q, nsymq, invsymq, rtau
  USE qpoint,        ONLY : xq
  USE control_flags, ONLY : modenum, noinv
  USE funct,         ONLY : dmxc, dmxc_spin, dmxc_nc, dft_is_gradient
  USE mp_global,     ONLY : world_comm
  USE mp,            ONLY : mp_bcast
  USE mp,            ONLY : mp_max, mp_min
  USE mp_pools,      ONLY : inter_pool_comm
  USE epwcom,        ONLY : xk_cryst
  USE fft_base,      ONLY : dfftp
  USE gvecs,         ONLY : doublegrid
  USE start_k,       ONLY : nk1, nk2, nk3
  !
  implicit none
  ! 
  INTEGER :: ir
  !! counter on mesh points
  INTEGER :: table (48, 48)
  !! the multiplication table of the point g
  INTEGER :: isym
  !! counter on symmetries
  INTEGER :: jsym
  !! counter on symmetries
  INTEGER :: ik
  !! counter on rotations
  INTEGER :: jk
  !! counter on k points
  INTEGER :: ibnd
  !! counter on k points
  INTEGER :: ipol
  !! counter on bands
  INTEGER :: mu
  !! counter on polarizations
  INTEGER :: imode0
  !! counter on modes
  INTEGER :: irr
  !! the starting mode
  INTEGER :: ipert
  !! counter on representation and perturbat
  INTEGER :: na
  !! counter on atoms
  INTEGER :: it
  !! counter on iterations
  INTEGER :: is
  !! counter on atomic type
  INTEGER :: js
  !! counter on atomic type
  INTEGER :: last_irr_eff
  !! Last effective irr
  INTEGER, ALLOCATABLE :: ifat(:)
  !!  
  REAL(kind=DP) :: rhotot
  !! total charge
  REAL(kind=DP) :: rhoup
  !! total up charge
  REAL(kind=DP) :: rhodw
  !! total down charge
  REAL(kind=DP) :: target
  !! 
  REAL(kind=DP) :: small
  !! 
  REAL(kind=DP) :: fac
  !! 
  REAL(kind=DP) :: xmax
  !! to set nbnd_occ in the metallic case
  REAL(kind=DP) :: emin
  !! minimum band energy
  REAL(kind=DP) ::emax
  !! maximum band energy
  REAL(kind=DP) :: xx_c
  !! 
  REAL(kind=DP) :: yy_c
  !! 
  REAL(kind=DP) :: zz_c
  !! 
  REAL(kind=DP) :: eps
  !! 
  REAL(kind=DP) :: auxdmuxc(4,4)
  !!
  LOGICAL :: magnetic_sym
  !! the symmetry operations
  LOGICAL :: symmorphic_or_nzb
  !!
  !
  CALL start_clock ('epw_setup')
  !
  ! 0) Set up list of kpoints in crystal coordinates
  !
  DO jk = 1, nkstot
     xk_cryst(:,jk) = xk(:,jk)
  END DO
  CALL cryst_to_cart (nkstot, xk_cryst, at, -1)
  CALL mp_bcast(xk_cryst,ionode_id,world_comm)
  !
  !  loosy tolerance: not important 
  eps = 1.d-5
  DO jk = 1, nkstot
     xx_c = xk_cryst(1,jk)*nk1
     yy_c = xk_cryst(2,jk)*nk2
     zz_c = xk_cryst(3,jk)*nk3
     !
     ! check that the k-mesh was defined in the positive region of 1st BZ
     !
     IF ( xx_c .lt. -eps .or. yy_c .lt. -eps .or. zz_c .lt. -eps ) &
        call errore('epw_setup','coarse k-mesh needs to be strictly positive in 1st BZ',1)
     !
  ENDDO
  !
  ! 1) Computes the total local potential (external+scf) on the smooth grid
  !
  CALL set_vrs (vrs, vltot, v%of_r, kedtau, v%kin_r, dfftp%nnr, nspin, doublegrid)
  !
  ! 2) Set non linear core correction stuff
  !
  nlcc_any = ANY ( upf(1:ntyp)%nlcc )
  IF (nlcc_any) allocate (drc( ngm, ntyp))    
  !
  !  3) If necessary calculate the local magnetization. This information is
  !      needed in sgama 
  !
  IF (.not.ALLOCATED(m_loc)) ALLOCATE( m_loc( 3, nat ) )
  IF (noncolin.and.domag) THEN
     DO na = 1, nat
        !
        m_loc(1,na) = starting_magnetization(ityp(na)) * &
                      SIN( angle1(ityp(na)) ) * COS( angle2(ityp(na)) )
        m_loc(2,na) = starting_magnetization(ityp(na)) * &
                      SIN( angle1(ityp(na)) ) * SIN( angle2(ityp(na)) )
        m_loc(3,na) = starting_magnetization(ityp(na)) * &
                      COS( angle1(ityp(na)) )
     ENDDO
     ux=0.0_DP
     IF (dft_is_gradient()) CALL compute_ux(m_loc,ux,nat)
  ENDIF
  !
  ! 3) Computes the derivative of the xc potential
  !
  dmuxc(:,:,:) = 0.d0
  IF (lsda) THEN
     DO ir = 1, dfftp%nnr
        rhoup = rho%of_r (ir, 1) + 0.5d0 * rho_core (ir)
        rhodw = rho%of_r (ir, 2) + 0.5d0 * rho_core (ir)
        CALL dmxc_spin (rhoup, rhodw, dmuxc(ir,1,1), dmuxc(ir,2,1), &
                                      dmuxc(ir,1,2), dmuxc(ir,2,2) )
     ENDDO
  ELSE
     IF (noncolin.and.domag) THEN
        DO ir = 1, dfftp%nnr
           rhotot = rho%of_r (ir, 1) + rho_core (ir)
           CALL dmxc_nc (rhotot, rho%of_r(ir,2), rho%of_r(ir,3), rho%of_r(ir,4), auxdmuxc)
           DO is=1,nspin
              DO js=1,nspin
                 dmuxc(ir,is,js)=auxdmuxc(is,js)
              ENDDO
           ENDDO
        ENDDO
     ELSE
        DO ir = 1, dfftp%nnr
           rhotot = rho%of_r (ir, 1) + rho_core (ir)
           IF (rhotot.gt.1.d-30) dmuxc (ir, 1, 1) = dmxc (rhotot)
           IF (rhotot.lt. - 1.d-30) dmuxc (ir, 1, 1) = - dmxc ( - rhotot)
        ENDDO
     ENDIF
  ENDIF
  !
  ! 3.1) Setup all gradient correction stuff
  !
  CALL setup_dgc
  !
  ! 4) Computes the inverse of each matrix
  !
  CALL multable (nsym, s, table)
  DO isym = 1, nsym
     DO jsym = 1, nsym
        IF (table (isym, jsym) == 1) invs (isym) = jsym
     ENDDO
  ENDDO
  !
  ! 5) Computes the number of occupied bands for each k point
  !
  IF (lgauss) THEN
     !
     ! discard conduction bands such that w0gauss(x,n) < small
     !
     ! hint:
     !   small = 1.0333492677046d-2  ! corresponds to 2 gaussian sigma
     !   small = 6.9626525973374d-5  ! corresponds to 3 gaussian sigma
     !   small = 6.3491173359333d-8  ! corresponds to 4 gaussian sigma
     !
     small = 6.9626525973374d-5
     !
     ! - appropriate limit for gaussian broadening (used for all ngauss)
     !
     xmax = sqrt ( - log (sqrt (pi) * small) )
     !
     ! - appropriate limit for Fermi-Dirac
     !
     IF (ngauss.eq. - 99) THEN
        fac = 1.d0 / sqrt (small)
        xmax = 2.d0 * log (0.5d0 * (fac + sqrt (fac * fac - 4.d0) ) )
     ENDIF
     target = ef + xmax * degauss
     DO ik = 1, nks
        DO ibnd = 1, nbnd
           IF (et (ibnd, ik) .lt.target) nbnd_occ (ik) = ibnd
        ENDDO
        IF (nbnd_occ (ik) .eq.nbnd) WRITE( stdout, '(5x,/,&
             &"Possibly too few bands at point ", i4,3f10.5)') &
             ik,  (xk (ipol, ik) , ipol = 1, 3)
     ENDDO
  ELSEIF (ltetra) THEN
     CALL errore('epw_setup','phonon + tetrahedra not implemented', 1)
  ELSE
     IF (lsda) call infomsg('epw_setup','occupation numbers probably wrong')
     IF (noncolin) THEN
        nbnd_occ = nint (nelec) 
     ELSE
        DO ik = 1, nks
           nbnd_occ (ik) = nint (nelec) / nint(degspin)
        ENDDO
     ENDIF
  ENDIF
  !
  ! 6) Computes alpha_pv
  !
  emin = et (1, 1)
  DO ik = 1, nks
     DO ibnd = 1, nbnd
        emin = min (emin, et (ibnd, ik) )
     ENDDO
  ENDDO
  ! find the minimum across pools
  CALL mp_min( emin, inter_pool_comm )
  IF (lgauss) THEN
     emax = target
     alpha_pv = emax - emin
  ELSE
     emax = et (1, 1)
     DO ik = 1, nks
        DO ibnd = 1, nbnd_occ(ik)
           emax = max (emax, et (ibnd, ik) )
        ENDDO
     ENDDO
     ! find the maximum across pools
     CALL mp_max( emax, inter_pool_comm )
     alpha_pv = 2.d0 * (emax - emin)
  ENDIF
  ! avoid zero value for alpha_pv
  alpha_pv = max (alpha_pv, 1.0d-2)
  !
  ! 7) set all the variables needed to use the pattern representation
  !
  magnetic_sym = noncolin .AND. domag
  time_reversal = .NOT. noinv .AND. .NOT. magnetic_sym

  nmodes = 3 * nat

!  IF (npool > 1) THEN
!     CALL xk_collect( xk_col, xk, nkstot, nks )
!  ELSE
!     xk_col(:,1:nks) = xk(:,1:nks)
!  ENDIF
  !
  !   The small group of q may be known. At a given q it is calculated
  !   by set_nscf, at gamma it coincides with the point group and we
  !   take nsymq=nsym
  !
  IF (lgamma.AND.modenum==0) THEN
     nsymq=nsym
     minus_q=.TRUE.
  ENDIF
  !
  !   If the code arrives here and nsymq is still 0 the small group of q has 
  !   not been calculated by set_nscf because this is a recover run. 
  !   We recalculate here the small group of q.
  !
  IF (nsymq==0) CALL set_small_group_of_q(nsymq, invsymq, minus_q)
  IF ( .NOT. time_reversal ) minus_q = .FALSE.
  !
  !
  IF (modenum > 0) THEN
     search_sym=.FALSE.
     minus_q = .FALSE.
  ENDIF
  !
  ! allocate and calculate rtau, the bravais lattice vector associated
  ! to a rotation
  !
  call sgam_ph_new (at, bg, nsym, s, irt, tau, rtau, nat)
  !
  !    and calculate the vectors G associated to the symmetry Sq = q + G
  !    if minus_q is true calculate also irotmq and the G associated to Sq=-g+G
  !
  CALL set_giq (xq,s,nsymq,nsym,irotmq,minus_q,gi,gimq)

  search_sym = search_sym .AND. symmorphic_or_nzb()
  !
  num_rap_mode=-1
  ! 
  IF (search_sym) CALL prepare_sym_analysis(nsymq,sr,t_rev,magnetic_sym)
  !
  IF (.NOT.u_from_file) THEN
  ! SP: THIS CALLS SET THE u
     CALL find_irrep()
  ENDIF
  CALL find_irrep_sym()
  !
  IF (lgamma_gamma) THEN
     ALLOCATE(has_equivalent(nat))
     ALLOCATE(with_symmetry(3*nat))
     ALLOCATE(n_equiv_atoms(nat))
     ALLOCATE(equiv_atoms(nat,nat))
     CALL find_equiv_sites (nat,nsym,irt,has_equivalent,n_diff_sites, &
                       n_equiv_atoms,equiv_atoms)

     IF (n_diff_sites .LE. 0 .OR. n_diff_sites .GT. nat)            &
          &      CALL errore('phq_setup','problem with n_diff_sites',1)
     !
     ! look if ASR can be exploited to reduce the number of calculations
     ! we need to locate an independent atom with no equivalent atoms
     nasr=0
     IF (asr.AND.n_diff_sites.GT.1) THEN
        DO na = 1, n_diff_sites
           IF (n_equiv_atoms(na).EQ.1 ) THEN
              nasr = equiv_atoms(na, 1)
              GO TO 1
           END IF
        END DO
 1      CONTINUE
     END IF
  END IF


  if (fildrho.ne.' '.and.ionode) call io_pattern(nat,fildrho,nirr,npert,u,xq,tmp_dir,+1)

  if (start_irr < 0) call errore('phq_setup', 'wrong start_irr', 1)
  last_irr_eff=last_irr
  if (last_irr > nirr.or.last_irr<0) last_irr_eff=nirr
  !
  !  set the alpha_mix parameter
  !
  do it = 2, niter_ph
     if (ABS(alpha_mix (it)) < eps ) alpha_mix (it) = alpha_mix (it - 1)
  enddo
  !
  ! Set flmixdpot
  !
  if (reduce_io) then
     flmixdpot = ' '
  else
     flmixdpot = 'mixd'
  endif
  !
  !  8) set the variables needed for the partial computation: nat_todo, 
  !     atomo, comp_irr
  !
  if (lgamma_gamma) then
     with_symmetry=1
     comp_irr = .FALSE.
     comp_irr(0)=.TRUE.
     do na=1,nat
        if (has_equivalent(na)==0) then
            do ipol=1,3
               comp_irr(3*(na-1)+ipol)=.TRUE.
               with_symmetry(3*(na-1)+ipol)=0
            enddo
        endif
     enddo
     if (nasr>0) then
        do ipol=1,3
           comp_irr(3*(nasr-1)+ipol)=.FALSE.
           with_symmetry(3*(nasr-1)+ipol)=0
        enddo
     endif
     IF (start_irr <= last_irr_eff) THEN
        DO irr=1,start_irr-1
           comp_irr(irr) = .FALSE.
        ENDDO
        DO irr=last_irr_eff+1,3*nat
           comp_irr(irr) = .FALSE.
        ENDDO
     ENDIF
  endif
  !
  !  Compute how many atoms moves and set the list atomo
  !
  ALLOCATE(ifat(nat))
  ifat = 0
  imode0 = 0
  DO irr = 1, nirr
     if (comp_irr (irr)) then
        do ipert = 1, npert (irr)
           do na = 1, nat
              do ipol = 1, 3
                 mu = 3 * (na - 1) + ipol
                 if (abs (u (mu, imode0+ipert) ) > 1.d-12) ifat (na) = 1
              enddo
           enddo
        enddo
     endif
     imode0 = imode0 + npert (irr)
  ENDDO
  nat_todo = 0
  DO na = 1, nat
     IF (ifat (na) == 1) THEN
        nat_todo = nat_todo + 1
        atomo (nat_todo) = na
     ENDIF
  ENDDO

  DEALLOCATE(ifat)
  !
  !   Initialize done_irr, find max dimension of the irreps
  !
  all_comp=.true.
  DO irr=1,nirr
     IF (.NOT.comp_irr(irr)) all_comp=.false.
  ENDDO
  all_comp = all_comp.OR.lgamma_gamma
  all_done = .FALSE.
  npertx = 0
  done_irr = .FALSE.
  DO irr = 1, nirr
     npertx = max (npertx, npert (irr) )
  ENDDO

  CALL stop_clock ('epw_setup')
  RETURN
  !
  END SUBROUTINE epw_setup
  !
  !-----------------------------------------------------------------------
  SUBROUTINE epw_setup_restart
  !-----------------------------------------------------------------------
  !!
  !! Setup in the case of a restart
  !! 
  ! ----------------------------------------------------------------------
  USE kinds,         ONLY : DP
  USE constants_epw, ONLY : zero
  USE io_global,     ONLY : ionode_id
  USE mp_global,     ONLY : world_comm
  USE mp,            ONLY : mp_bcast
  USE epwcom,        ONLY : scattering, nstemp, tempsmin, tempsmax, &
                            temps
  USE transportcom,  ONLY : transp_temp
  !
  implicit none
  !
  INTEGER :: itemp
  !! Counter on temperature
  ! 
  CALL start_clock ('epw_setup')
  !
  IF (.NOT. ALLOCATED(transp_temp)) ALLOCATE( transp_temp(nstemp) )
  !
  transp_temp(:) = zero
  ! In case of scattering calculation
  IF ( scattering ) THEN
    ! 
    IF ( maxval(temps(:)) > zero ) THEN
      transp_temp(:) = temps(:)
    ELSE
      IF ( nstemp .eq. 1 ) THEN
        transp_temp(1) = tempsmin
      ELSE
        DO itemp = 1, nstemp
          transp_temp(itemp) = tempsmin + dble(itemp-1) * &
                              ( tempsmax - tempsmin ) / dble(nstemp-1)
        ENDDO
      ENDIF
    ENDIF
  ENDIF
  ! We have to bcast here because before it has not been allocated
  CALL mp_bcast (transp_temp, ionode_id, world_comm)    !  
  ! 
  CALL stop_clock ('epw_setup')
  RETURN
  !-----------------------------------------------------------------------
  END SUBROUTINE epw_setup_restart
  !-----------------------------------------------------------------------



