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
  !! RM - Nov 2018: Updated based on QE 6.3
  !! 
  USE kinds,         ONLY : DP
  USE ions_base,     ONLY : tau, nat, ntyp => nsp, ityp
  USE cell_base,     ONLY : at, bg  
  USE io_global,     ONLY : stdout, ionode, ionode_id
  USE io_files,      ONLY : tmp_dir
  USE klist,         ONLY : xk, nks, nkstot
  USE lsda_mod,      ONLY : nspin, starting_magnetization
  USE scf,           ONLY : v, vrs, vltot, rho, kedtau
  USE gvect,         ONLY : ngm
  USE symm_base,     ONLY : nsym, s, irt, t_rev, time_reversal, invs, sr, &
                            inverse_s
  USE uspp_param,    ONLY : upf
  USE spin_orb,      ONLY : domag
  USE constants_epw, ONLY : zero, eps5
  USE noncollin_module,     ONLY : noncolin, m_loc, angle1, angle2, ux
  USE nlcc_ph,       ONLY : drc
  USE uspp,          ONLY : nlcc_any
  USE control_ph,    ONLY : search_sym, u_from_file
  USE control_lr,    ONLY : alpha_pv, nbnd_occ
  USE modes,         ONLY : u, npertx, npert, nirr, nmodes, num_rap_mode
  USE lr_symm_base,  ONLY : gi, gimq, irotmq, minus_q, nsymq, invsymq, rtau
  USE qpoint,        ONLY : xq
  USE control_flags, ONLY : modenum, noinv
  USE funct,         ONLY : dft_is_gradient
  USE mp_global,     ONLY : world_comm
  USE mp,            ONLY : mp_bcast
  USE mp_pools,      ONLY : inter_pool_comm
  USE epwcom,        ONLY : xk_cryst, scattering, nstemp, tempsmin, tempsmax, &
                            temps
  USE fft_base,      ONLY : dfftp
  USE gvecs,         ONLY : doublegrid
  USE start_k,       ONLY : nk1, nk2, nk3
  USE transportcom,  ONLY : transp_temp
  !
  IMPLICIT NONE
  ! 
  INTEGER :: jk
  !! counter on k points
  INTEGER :: irr
  !! counter on irrepr
  INTEGER :: na
  !! counter on atoms
  INTEGER :: itemp
  !! counter on temperatures 
  REAL(kind=DP) :: xx_c, yy_c, zz_c
  !! k-points in crystal coords. in multiple of nk1, nk2, nk3
  LOGICAL :: magnetic_sym
  !! the symmetry operations
  LOGICAL :: symmorphic_or_nzb
  !!
  !
  CALL start_clock('epw_setup')
  !
  ! 0) Set up list of kpoints in crystal coordinates
  !
  DO jk = 1, nkstot
    xk_cryst(:,jk) = xk(:,jk)
  ENDDO
  !  bring k-points from cartesian to crystal coordinates
  CALL cryst_to_cart(nkstot, xk_cryst, at, -1)
  CALL mp_bcast(xk_cryst,ionode_id,world_comm)
  !
  !  loosy tolerance: not important 
  DO jk = 1, nkstot
    xx_c = xk_cryst(1,jk) * nk1
    yy_c = xk_cryst(2,jk) * nk2
    zz_c = xk_cryst(3,jk) * nk3
    !
    ! check that the k-mesh was defined in the positive region of 1st BZ
    !
    IF ( xx_c .lt. -eps5 .or. yy_c .lt. -eps5 .or. zz_c .lt. -eps5 ) &
      CALL errore('epw_setup','coarse k-mesh needs to be strictly positive in 1st BZ',1)
    !
  ENDDO
  !
  ! 1) Computes the total local potential (external+scf) on the smooth grid
  !
  CALL set_vrs(vrs, vltot, v%of_r, kedtau, v%kin_r, dfftp%nnr, nspin, doublegrid)
  !
  ! Set non linear core correction stuff
  !
  nlcc_any = ANY( upf(1:ntyp)%nlcc )
  IF (nlcc_any) ALLOCATE(drc(ngm, ntyp))    
  !
  !  2) If necessary calculate the local magnetization. This information is
  !      needed in sgama 
  !
  IF (.not.ALLOCATED(m_loc)) ALLOCATE(m_loc(3, nat))
  IF (noncolin .AND. domag) THEN
    DO na = 1, nat
      !
      m_loc(1,na) = starting_magnetization(ityp(na)) * &
                    SIN( angle1(ityp(na)) ) * COS( angle2(ityp(na)) )
      m_loc(2,na) = starting_magnetization(ityp(na)) * &
                    SIN( angle1(ityp(na)) ) * SIN( angle2(ityp(na)) )
      m_loc(3,na) = starting_magnetization(ityp(na)) * &
                    COS( angle1(ityp(na)) )
    ENDDO
    ux = zero
    IF (dft_is_gradient()) CALL compute_ux(m_loc,ux,nat)
  ENDIF
  !
  ! 3) Computes the derivative of the xc potential
  !
  CALL setup_dmuxc()
  !
  ! 3.1) Setup all gradient correction stuff
  !
  CALL setup_dgc
  !  
  ! 4) Computes the inverse of each matrix of the crystal symmetry group
  !
  CALL inverse_s()
  !
  ! 5) Computes the number of occupied bands for each k point
  !
  CALL setup_nbnd_occ()
  !
  ! 6) Computes alpha_pv
  !
  CALL setup_alpha_pv()
  !
  ! 7) set all the variables needed to use the pattern representation
  !
  magnetic_sym = noncolin .AND. domag
  time_reversal = .NOT. noinv .AND. .NOT. magnetic_sym
  !
  nmodes = 3 * nat
  !
  !   If the code arrives here and nsymq is still 0 the small group of q has 
  !   not been calculated by set_nscf because this is a recover run. 
  !   We recalculate here the small group of q.
  !
  IF (nsymq==0) CALL set_small_group_of_q(nsymq, invsymq, minus_q)
  IF ( .NOT. time_reversal ) minus_q = .FALSE.
  !
  IF (modenum > 0) THEN
    search_sym = .FALSE.
    minus_q = .FALSE.
  ENDIF
  !
  ! allocate and calculate rtau, the Bravais lattice vector associated
  ! to a rotation
  !
  CALL sgam_lr(at, bg, nsym, s, irt, tau, rtau, nat)
  !
  !    and calculate the vectors G associated to the symmetry Sq = q + G
  !    if minus_q is true calculate also irotmq and the G associated to Sq=-g+G
  !
  CALL set_giq(xq, s, nsymq, nsym, irotmq, minus_q, gi, gimq)
  !
  search_sym = search_sym .AND. symmorphic_or_nzb()
  !
  num_rap_mode = -1
  IF (search_sym) CALL prepare_sym_analysis(nsymq, sr, t_rev, magnetic_sym)
  !
  IF (.NOT. u_from_file) THEN
  ! SP: These calls set the u
     CALL find_irrep()
  ENDIF
  CALL find_irrep_sym()
  !
  !  8) set max perturbation
  !   
  npertx = 0
  DO irr = 1, nirr
    npertx = max(npertx, npert(irr))
  ENDDO
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
  CALL mp_bcast(transp_temp, ionode_id, world_comm)    !  
  ! 
  CALL stop_clock('epw_setup')
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
  USE constants_epw, ONLY : zero
  USE io_global,     ONLY : ionode_id
  USE mp_global,     ONLY : world_comm
  USE mp,            ONLY : mp_bcast
  USE epwcom,        ONLY : scattering, nstemp, tempsmin, tempsmax, &
                            temps
  USE transportcom,  ONLY : transp_temp
  !
  IMPLICIT NONE
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
  ! 
  ! We have to bcast here because before it has not been allocated
  CALL mp_bcast (transp_temp, ionode_id, world_comm)    !  
  ! 
  CALL stop_clock ('epw_setup')
  !
  RETURN
  !-----------------------------------------------------------------------
  END SUBROUTINE epw_setup_restart
  !-----------------------------------------------------------------------
