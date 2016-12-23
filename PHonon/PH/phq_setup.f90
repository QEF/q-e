!
! Copyright (C) 2001-2016 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
subroutine phq_setup
  !-----------------------------------------------------------------------
  !
  !  This subroutine prepares several variables which are needed in the
  !  phonon program:
  !  1) computes the total local potential (external+scf) on the smooth
  !     grid to be used in h_psi and similia
  !  2) computes the local magnetization (if necessary)
  !  3) computes dmuxc (with GC if needed)
  !  4) set the inverse of every matrix invs
  !  5) for metals sets the occupied bands
  !  6) computes alpha_pv
  !  7) computes the variables needed to pass to the pattern representation
  !     u      the patterns
  !     t      the matrices of the small group of q on the pattern basis
  !     tmq    the matrix of the symmetry which sends q -> -q + G
  !     gi     the G associated to each symmetry operation
  !     gimq   the G of the q -> -q+G symmetry
  !     nsymq  the order of the small group of q
  !     irotmq the index of the q->-q+G symmetry
  !     nirr   the number of irreducible representation
  !     npert  the dimension of each irreducible representation
  !     nmodes the number of modes
  !     minus_q true if there is a symmetry sending q -> -q+G
  !  8) for testing purposes it sets ubar
  !  9) set the variables needed to deal with nlcc
  !  10) set the variables needed for the partial computation
  !       of the dynamical matrix
  !
  !  IMPORTANT NOTE ABOUT SYMMETRIES:
  !  nrot  is the number of sym.ops. of the Bravais lattice
  !        read from data file, only used in set_default_pw
  !  nsym  is the number of sym.ops. of the crystal symmetry group
  !        read from data file, should never be changed
  !  nsymq is the number of sym.ops. of the small group of q
  !        it is calculated in set_defaults_pw for each q
  !  The matrices "s" of sym.ops are ordered as follows:
  !   first the nsymq sym.ops. of the small group of q
  !   (the ordering is done in subroutine copy_sym in set_defaults_pw),
  !   followed by the remaining nsym-nsymq sym.ops. of the crystal group,
  !   followed by the remaining nrot-nsym sym.ops. of the Bravais  group
  !
  !
  USE kinds,         ONLY : DP
  USE ions_base,     ONLY : tau, nat, ntyp => nsp, ityp
  USE cell_base,     ONLY : at, bg
  USE io_global,     ONLY : ionode
  USE io_files,      ONLY : tmp_dir
  USE klist,         ONLY : xk, nks, nkstot
  USE lsda_mod,      ONLY : nspin, starting_magnetization
  USE scf,           ONLY : v, vrs, vltot, kedtau
  USE fft_base,      ONLY : dfftp
  USE gvect,         ONLY : ngm
  USE gvecs,       ONLY : doublegrid
  USE symm_base,     ONLY : nrot, nsym, s, ftau, irt, t_rev, time_reversal, &
                            sr, invs, inverse_s
  USE uspp_param,    ONLY : upf
  USE uspp,          ONLY : nlcc_any
  USE spin_orb,      ONLY : domag
  USE noncollin_module, ONLY : noncolin, m_loc, angle1, angle2, ux
  USE nlcc_ph,       ONLY : drc
  USE control_ph,    ONLY : rec_code, lgamma_gamma, search_sym, start_irr, &
                            last_irr, niter_ph, alpha_mix, all_done,  &
                            trans, epsil, recover, where_rec, &
                            flmixdpot, reduce_io, rec_code_read, &
                            done_epsil, zeu, done_zeu, current_iq, u_from_file
  USE el_phon,       ONLY : elph, comp_elph, done_elph, elph_nbnd_min, elph_nbnd_max
  USE output,        ONLY : fildrho
  USE modes,         ONLY : u, npertx, npert, nirr, t, tmq, nmodes, num_rap_mode
  USE dynmat,        ONLY : dyn, dyn_rec, dyn00
  USE efield_mod,    ONLY : epsilon, zstareu
  USE partial,       ONLY : comp_irr, atomo, nat_todo, all_comp, &
                            done_irr
  USE gamma_gamma,   ONLY : has_equivalent, asr, nasr, n_diff_sites, &
                            equiv_atoms, n_equiv_atoms, with_symmetry
  USE ph_restart,    ONLY : ph_writefile, ph_readfile
  USE control_flags, ONLY : modenum, noinv
  USE grid_irr_iq,   ONLY : comp_irr_iq
  USE funct,         ONLY : dft_is_gradient
  USE ramanm,        ONLY : lraman, elop, ramtns, eloptns, done_lraman, &
                            done_elop

  USE mp_pools,      ONLY : inter_pool_comm, npool
  !
  USE acfdtest,      ONLY : acfdt_is_active, acfdt_num_der
  USE elph_tetra_mod, ONLY : elph_tetra
  USE wvfct,         ONLY : nbnd, et
  USE ener,          ONLY : ef
  USE mp,            ONLY : mp_max, mp_min

  USE lr_symm_base,  ONLY : gi, gimq, irotmq, minus_q, invsymq, nsymq, rtau
  USE qpoint,        ONLY : xq, xk_col
  USE control_lr,    ONLY : lgamma

  implicit none

  real(DP) :: sr_is(3,3,48)

  integer :: isym, jsym, irot, ik, ibnd, ipol, &
       mu, nu, imode0, irr, ipert, na, it, nt, nsym_is, last_irr_eff
  ! counters

  real(DP), allocatable :: wg_up(:,:), wg_dw(:,:)

  logical :: sym (48), magnetic_sym
  LOGICAL :: symmorphic_or_nzb
  ! the symmetry operations
  integer, allocatable :: ifat(:)
  integer :: ierr

  call start_clock ('phq_setup')
  ! 0) A few checks
  !
  IF (dft_is_gradient().and.(lraman.or.elop)) call errore('phq_setup', &
     'third order derivatives not implemented with GGA', 1)

  IF (nsymq==0) CALL errore('phq_setup', &
                 'The small group of q is no more calculated in phq_setup',1)
  !
  !  read the displacement patterns
  !
  IF (u_from_file) THEN
     CALL ph_readfile('data_u',current_iq,0,ierr)
     IF (ierr /= 0) CALL errore('phq_setup', 'problem with modes file',1)
  ENDIF
  !
  ! 1) Computes the total local potential (external+scf) on the smooth grid
  !
!!!!!!!!!!!!!!!!!!!!!!!! ACFDT TEST !!!!!!!!!!!!!!!!
  IF (acfdt_is_active) THEN
     ! discard set_vrs for numerical derivatives
     if (.not.acfdt_num_der) then 
        call set_vrs (vrs, vltot, v%of_r, kedtau, v%kin_r, dfftp%nnr, nspin, doublegrid)
     end if
  ELSE
     call set_vrs (vrs, vltot, v%of_r, kedtau, v%kin_r, dfftp%nnr, nspin, doublegrid)
  ENDIF
!!!!!!!!!!!!!!!!!!!!!!!!END OF  ACFDT TEST !!!!!!!!!!!!!!!!
  !
  ! Set non linear core correction stuff
  !
  nlcc_any = ANY ( upf(1:ntyp)%nlcc )
  if (nlcc_any) allocate (drc( ngm, ntyp))
  !
  ! 2) If necessary calculate the local magnetization. This information is
  !      needed in find_sym
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
     END DO
     ux=0.0_DP
     if (dft_is_gradient()) call compute_ux(m_loc,ux,nat)
  ENDIF
  !
  ! 3) Computes the derivative of the XC potential
  !
  call setup_dmuxc()
  !
  ! Setup all gradient correction stuff
  !
  call setup_dgc()
  !
  ! 4) Computes the inverse of each matrix of the crystal symmetry group
  !
  call inverse_s()
  !
  ! 5) Computes the number of occupied bands for each k point
  !
  call setup_nbnd_occ()
  !
  ! 6) Computes alpha_pv
  !
  call setup_alpha_pv()
  !
  ! 7) set all the variables needed to use the pattern representation
  !
  magnetic_sym = noncolin .AND. domag
  time_reversal = .NOT. noinv .AND. .NOT. magnetic_sym

  nmodes = 3 * nat

  ! FIXME: xk_col is used only to write the output
  IF (npool > 1) THEN
     CALL poolcollect( 3, nks, xk, nkstot, xk_col )
  ELSE
     xk_col(:,1:nks) = xk(:,1:nks)
  ENDIF
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

  num_rap_mode=-1
  IF (search_sym) CALL prepare_sym_analysis(nsymq,sr,t_rev,magnetic_sym)

  IF (.NOT.u_from_file) THEN
     CALL find_irrep()
     CALL ph_writefile('data_u',current_iq,0,ierr)
  ENDIF
  CALL find_irrep_sym()


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


  if (fildrho.ne.' '.and.ionode) call io_pattern (nat,fildrho,nirr,npert,u,xq,tmp_dir,+1)

  if (start_irr < 0) call errore('phq_setup', 'wrong start_irr', 1)
  last_irr_eff=last_irr
  if (last_irr > nirr.or.last_irr<0) last_irr_eff=nirr
  !
  !  set the alpha_mix parameter
  !
  do it = 2, niter_ph
     if (alpha_mix (it) .eq.0.d0) alpha_mix (it) = alpha_mix (it - 1)
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
  !
  ! 8) Set the ubar
  !
  ! ubar removed on 16/02/2012, used only for debugging
  !
  !  9) set the variables needed for the partial computation:
  !     nat_todo, atomo, comp_irr

  DO irr=0,nirr
     comp_irr(irr)=comp_irr_iq(irr,current_iq)
     IF (elph .AND. irr>0) comp_elph(irr)=comp_irr(irr)
  ENDDO
  !
  !  The gamma_gamma case needs a different treatment
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
  IF (elph) done_elph = .FALSE.
  DO irr = 1, nirr
     npertx = max (npertx, npert (irr) )
  ENDDO
!
!  set to zero the variable written on file
!
  dyn=(0.0_DP,0.0_DP)
  dyn00=(0.0_DP,0.0_DP)
  dyn_rec=(0.0_DP,0.0_DP)
  IF (epsil.and..not.done_epsil) epsilon=0.0_DP
  IF (zeu.and..not.done_zeu) zstareu=0.0_DP
  IF (lraman.and..not.done_lraman) ramtns=0.0_DP
  IF (elop.and..not.done_elop)  eloptns=0.0_DP

  where_rec='phq_setup.'
  rec_code=-40
  CALL ph_writefile('status_ph',current_iq,0,ierr)
  !
  ! Bands for the electron-phonon calculation
  !
  IF(elph_tetra == 1 .OR. elph_tetra == 2) THEN
     !
     elph_nbnd_min = nbnd
     elph_nbnd_max = 1
     !     
     DO ibnd = 1, nbnd
        IF(MAXVAL(et(ibnd, 1:nks)) > ef .AND. ibnd < elph_nbnd_min) elph_nbnd_min = ibnd
        IF(MINVAL(et(ibnd, 1:nks)) < ef .AND. ibnd > elph_nbnd_max) elph_nbnd_max = ibnd
     END DO
     !
     call mp_min(elph_nbnd_min, inter_pool_comm)
     call mp_max(elph_nbnd_max, inter_pool_comm)
     !
  END IF
  !
  CALL stop_clock ('phq_setup')
  RETURN
END SUBROUTINE phq_setup
