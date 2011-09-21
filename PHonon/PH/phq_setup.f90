!
! Copyright (C) 2001-2008 Quantum ESPRESSO group
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
  !  2) computes dmuxc 3) with GC if needed
  !  4) set the inverse of every matrix invs
  !  5) for metals sets the occupied bands
  !  6) computes alpha_pv
  !  7) computes the variables needed to pass to the pattern representation
  !     u      the patterns
  !     t      the matrices of the small group of q on the pattern basis
  !     tmq    the matrix of the symmetry which sends q -> -q + G
  !     gi     the G associated to each symmetry operation
  !     gimq   the G of the q -> -q+G symmetry
  !     irgq   the small group indices
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
  USE ions_base,     ONLY : tau, nat, ntyp => nsp, ityp, amass
  USE cell_base,     ONLY : at, bg
  USE io_global,     ONLY : stdout
  USE ener,          ONLY : ef, ef_up, ef_dw
  USE klist,         ONLY : xk, lgauss, degauss, ngauss, nks, nelec, nelup, &
                            neldw, two_fermi_energies, wk
  USE ktetra,        ONLY : ltetra, tetra
  USE lsda_mod,      ONLY : nspin, lsda, starting_magnetization, isk
  USE scf,           ONLY : v, vrs, vltot, rho, rho_core, kedtau
  USE fft_base,      ONLY : dfftp
  USE gvect,         ONLY : ngm
  USE gvecs,       ONLY : doublegrid
  USE symm_base,     ONLY : nrot, nsym, s, ftau, irt, t_rev, time_reversal, &
                            sname, sr, invs, inverse_s, copy_sym
  USE uspp_param,    ONLY : upf
  USE spin_orb,      ONLY : domag
  USE constants,     ONLY : degspin, pi
  USE noncollin_module, ONLY : noncolin, m_loc, angle1, angle2, ux, nspin_mag
  USE wvfct,         ONLY : nbnd, et
  USE nlcc_ph,       ONLY : drc, nlcc_any
  USE eqv,           ONLY : dmuxc
  USE control_ph,    ONLY : rec_code, lgamma_gamma, search_sym, start_irr, &
                            last_irr, niter_ph, alpha_mix, all_done, elph, &
                            trans, epsil, lgamma, recover, where_rec, alpha_pv,&
                            nbnd_occ, flmixdpot, reduce_io, rec_code_read, &
                            done_epsil, zeu, done_zeu, current_iq, u_from_file
  USE output,        ONLY : fildrho
  USE modes,         ONLY : u, ubar, npertx, npert, gi, gimq, nirr, &
                            t, tmq, irotmq, irgq, minus_q, &
                            nsymq, nmodes, rtau, name_rap_mode, num_rap_mode
  USE dynmat,        ONLY : dyn, dyn_rec, dyn00
  USE efield_mod,    ONLY : epsilon, zstareu
  USE qpoint,        ONLY : xq
  USE partial,       ONLY : comp_irr, atomo, nat_todo, all_comp, &
                            done_irr
  USE gamma_gamma,   ONLY : has_equivalent, asr, nasr, n_diff_sites, &
                            equiv_atoms, n_equiv_atoms, with_symmetry
  USE ph_restart,    ONLY : ph_writefile, ph_readfile
  USE control_flags, ONLY : iverbosity, modenum, noinv
  USE disp,          ONLY : comp_irr_iq
  USE funct,         ONLY : dmxc, dmxc_spin, dmxc_nc, dft_is_gradient
  USE ramanm,        ONLY : lraman, elop, ramtns, eloptns, done_lraman, &
                            done_elop

  USE mp,            ONLY : mp_max, mp_min
  USE mp_global,     ONLY : inter_pool_comm, nimage
  !
  USE acfdtest,      ONLY : acfdt_is_active, acfdt_num_der

  implicit none

  real(DP) :: rhotot, rhoup, rhodw, target, small, fac, xmax, emin, emax
  ! total charge
  ! total up charge
  ! total down charge
  ! auxiliary variables used
  ! to set nbnd_occ in the metallic case
  ! minimum band energy
  ! maximum band energy

  real(DP) :: sr_is(3,3,48)

  integer :: ir, isym, jsym, irot, ik, ibnd, ipol, &
       mu, nu, imode0, irr, ipert, na, it, nt, is, js, nsym_is, last_irr_eff
  ! counters

  real(DP) :: auxdmuxc(4,4)
  real(DP), allocatable :: w2(:), wg_up(:,:), wg_dw(:,:)

  logical :: sym (48), magnetic_sym
  ! the symmetry operations
  integer, allocatable :: ifat(:)
  integer :: ierr

  call start_clock ('phq_setup')
  ! 0) A few checks
  !
  IF (dft_is_gradient().and.(lraman.or.elop)) call errore('phq_setup', &
     'third order derivatives not implemented with GGA', 1)
  !
  !  read the displacement patterns
  !
  IF (u_from_file) THEN
     CALL ph_readfile('data_u',ierr)
     IF (ierr /= 0) CALL errore('ph_setup', 'problem with modes file',1)
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
  ! 2) Set non linear core correction stuff
  !
  nlcc_any = ANY ( upf(1:ntyp)%nlcc )
  if (nlcc_any) allocate (drc( ngm, ntyp))
  !
  !  3) If necessary calculate the local magnetization. This information is
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
  ! 3) Computes the derivative of the xc potential
  !
  dmuxc(:,:,:) = 0.d0
  if (lsda) then
     do ir = 1, dfftp%nnr
        rhoup = rho%of_r (ir, 1) + 0.5d0 * rho_core (ir)
        rhodw = rho%of_r (ir, 2) + 0.5d0 * rho_core (ir)
        call dmxc_spin (rhoup, rhodw, dmuxc(ir,1,1), dmuxc(ir,2,1), &
                                      dmuxc(ir,1,2), dmuxc(ir,2,2) )
     enddo
  else
     IF (noncolin.and.domag) THEN
        do ir = 1, dfftp%nnr
           rhotot = rho%of_r (ir, 1) + rho_core (ir)
           call dmxc_nc (rhotot, rho%of_r(ir,2), rho%of_r(ir,3), rho%of_r(ir,4), auxdmuxc)
           DO is=1,nspin_mag
              DO js=1,nspin_mag
                 dmuxc(ir,is,js)=auxdmuxc(is,js)
              END DO
           END DO
        enddo
     ELSE
        do ir = 1, dfftp%nnr
           rhotot = rho%of_r (ir, 1) + rho_core (ir)
           if (rhotot.gt.1.d-30) dmuxc (ir, 1, 1) = dmxc (rhotot)
           if (rhotot.lt. - 1.d-30) dmuxc (ir, 1, 1) = - dmxc ( - rhotot)
        enddo
     END IF
  endif
  !
  ! 3.1) Setup all gradient correction stuff
  !
  call setup_dgc
  !
  ! 4) Computes the inverse of each matrix of the crystal symmetry group
  !
  call inverse_s ( )
  !
  ! 5) Computes the number of occupied bands for each k point
  !
  if (lgauss) then
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
     if (ngauss.eq. - 99) then
        fac = 1.d0 / sqrt (small)
        xmax = 2.d0 * log (0.5d0 * (fac + sqrt (fac * fac - 4.d0) ) )
     endif
     target = ef + xmax * degauss
     do ik = 1, nks
        do ibnd = 1, nbnd
           if (et (ibnd, ik) .lt.target) nbnd_occ (ik) = ibnd
        enddo
        if (nbnd_occ (ik) .eq.nbnd) WRITE( stdout, '(5x,/,&
             &"Possibly too few bands at point ", i4,3f10.5)') &
             ik,  (xk (ipol, ik) , ipol = 1, 3)
     enddo
  else if (ltetra) then
     call errore('phq_setup','phonon + tetrahedra not implemented', 1)
  else
     if (noncolin) then
        nbnd_occ = nint (nelec)
     else
        IF ( two_fermi_energies ) THEN
           !
           ALLOCATE(wg_up(nbnd,nks))
           ALLOCATE(wg_dw(nbnd,nks))
           CALL iweights( nks, wk, nbnd, nelup, et, ef_up, wg_up, 1, isk )
           CALL iweights( nks, wk, nbnd, neldw, et, ef_dw, wg_dw, 2, isk )
           DO ik = 1, nks
              DO ibnd=1,nbnd
                 IF (isk(ik)==1) THEN
                    IF (wg_up(ibnd,ik) > 0.0_DP) nbnd_occ (ik) = nbnd_occ(ik)+1
                 ELSE
                    IF (wg_dw(ibnd,ik) > 0.0_DP) nbnd_occ (ik) = nbnd_occ(ik)+1
                 ENDIF
              ENDDO
           ENDDO
           !
           ! the following line to prevent NaN in Ef
           !
           ef = ( ef_up + ef_dw ) / 2.0_dp
           !
           DEALLOCATE(wg_up)
           DEALLOCATE(wg_dw)
        ELSE
          if (lsda) call infomsg('phq_setup', &
                                 'occupation numbers probably wrong')
           do ik = 1, nks
              nbnd_occ (ik) = nint (nelec) / degspin
           enddo
        ENDIF
     endif
  endif
  !
  ! 6) Computes alpha_pv
  !
  emin = et (1, 1)
  do ik = 1, nks
     do ibnd = 1, nbnd
        emin = min (emin, et (ibnd, ik) )
     enddo
  enddo
#ifdef __PARA
  ! find the minimum across pools
  call mp_min( emin, inter_pool_comm )
#endif
  if (lgauss) then
     emax = target
     alpha_pv = emax - emin
  else
     emax = et (1, 1)
     do ik = 1, nks
        do ibnd = 1, nbnd_occ(ik)
           emax = max (emax, et (ibnd, ik) )
        enddo
     enddo
#ifdef __PARA
     ! find the maximum across pools
     call mp_max( emax, inter_pool_comm )
#endif
     alpha_pv = 2.d0 * (emax - emin)
  endif
  ! avoid zero value for alpha_pv
  alpha_pv = max (alpha_pv, 1.0d-2)
  !
  ! 7) set all the variables needed to use the pattern representation
  !
  magnetic_sym = noncolin .AND. domag
  time_reversal = .NOT. noinv .AND. .NOT. magnetic_sym
  !
  ! allocate and calculate rtau, the rotated position of each atom
  !
  sym (1:nsym) = .true.
  call sgam_ph (at, bg, nsym, s, irt, tau, rtau, nat, sym)
  !
  nmodes = 3 * nat
  minus_q = (modenum .eq. 0)
  ! if minus_q=.t. set_irr will search for Sq=-q+G symmetry.
  ! On output minus_q=.t. if such a symmetry has been found
  ! TEMP: set_irr_* should not find again the small group of q
  IF (lgamma) nsymq=nsym
  isym=nsymq
  allocate( w2(3*nat) )
  if (modenum .ne. 0) then
     ! workaround: isym in the following call should be nsymq,
     ! but nsymq is re-calculated inside, so a copy is needed
     IF (isym==0) THEN
        sym(1:nsym)=.true.
        call smallg_q (xq, modenum, at, bg, nsym, s, ftau, sym, minus_q)
        call sgam_ph (at, bg, nsym, s, irt, tau, rtau, nat, sym)
        call mode_group (modenum, xq, at, bg, nat, nsym, s, irt, &
                         minus_q, rtau, sym)
        isym = copy_sym ( nsym, sym )
     ENDIF
     npertx=1
     CALL allocate_pert()
     call set_irr_mode (nat, at, bg, xq, s, invs, isym, rtau, irt, &
          irgq, nsymq, minus_q, irotmq, t, tmq, npertx, u, npert, &
          nirr, gi, gimq, iverbosity, modenum)
  else
     if (nsym > 1.and..not.lgamma_gamma.and.(trans.or.zeu.or.elph)) then
        call set_irr (nat, at, bg, xq, s, sr, tau, ntyp, ityp, ftau, invs,&
                    nsym, rtau, irt, irgq, nsymq, minus_q, irotmq, u, npert, &
                    nirr, gi, gimq, iverbosity, u_from_file, w2, search_sym, &
                    nspin_mag, t_rev, amass, num_rap_mode, name_rap_mode)
        npertx = 0
        DO irr = 1, nirr
           npertx = max (npertx, npert (irr) )
        ENDDO
        CALL allocate_pert()
        CALL set_irr_sym (nat, at, bg, xq, s, rtau, irt, irgq, nsymq,  &
                          minus_q, irotmq, t, tmq, u, npert, nirr, npertx )
     else
        search_sym=.FALSE.
        npertx=1
        CALL allocate_pert()
        call set_irr_nosym (nat, at, bg, xq, s, invs, nsym, rtau, irt, &
             irgq, nsymq, minus_q, irotmq, t, tmq, npertx, u, npert, &
             nirr, gi, gimq, iverbosity)
        IF (lgamma_gamma) THEN
           search_sym=.TRUE.
           CALL prepare_sym_analysis(nsymq,sr,t_rev,magnetic_sym)
        ENDIF
     endif
  endif
  IF ( isym /= nsymq ) CALL errore('phq_setup',&
             'internal error: mismatch in the order of small group',isym)
  IF (.NOT.time_reversal) minus_q=.false.

  IF (lgamma_gamma) THEN
     ALLOCATE(has_equivalent(nat))
     ALLOCATE(with_symmetry(3*nat))
     ALLOCATE(n_equiv_atoms(nat))
     ALLOCATE(equiv_atoms(nat,nat))
     CALL find_equiv_sites (nat,nat,nsym,irt,has_equivalent,n_diff_sites, &
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


  if (fildrho.ne.' ') call io_pattern (nat,fildrho,nirr,npert,u,+1)

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

  ubar(:) =( 0.d0,0.d0)
  !
  !   NB: the following instructions are for testing purposes of delta rho
  !       the user must know how many atoms there are in the system
  !
  !      ubar(1)=(1.d-3,0.d0)
  !      ubar(5)=(1.d0,0.d0)
  !      ubar(6)=(1.d0,0.d0)
  !
  !  9) set the variables needed for the partial computation:
  !     nat_todo, atomo, comp_irr
  ALLOCATE(ifat(nat))
  comp_irr = 0
  comp_irr(0)=1
  IF (nat_todo==0.AND.modenum==0) THEN
     !
     !  Case 1)  The partial computation option is not used, make all
     !           representation between start_irr and last_irr
     !
     IF (start_irr <= last_irr_eff) comp_irr(start_irr: last_irr_eff) = 1
     !
  ELSEIF (nat_todo /= 0) THEN
     !
     !   Case 2) Sets the atoms which must be computed: the requested
     !   atoms and all the symmetry related atoms
     !
     ifat = 0
     DO na = 1, nat_todo
        ifat (atomo (na) ) = 1
        DO isym = 1, nsymq
           irot = irgq (isym)
           ifat (irt (irot, atomo (na) ) ) = 1
        ENDDO
     ENDDO
     !
     !  Find the irreducible representations where the required atoms moves
     !
     imode0 = 0
     do irr = 1, nirr
        do ipert = 1, npert (irr)
           mu = imode0 + ipert
           do na = 1, nat
              if (ifat (na) == 1 .and. comp_irr (irr) == 0) then
                 do ipol = 1, 3
                    nu = 3 * (na - 1) + ipol
                    if (abs (u (nu, mu) ) > 1.d-6)  comp_irr (irr) = 1
                 enddo
              endif
           enddo
        enddo
        imode0 = imode0 + npert (irr)
     enddo
  ELSEIF (modenum /= 0) THEN
     comp_irr(modenum)=1
  ELSE
     call errore('phq_setup','nat_todo or nrap wrong',1)
  ENDIF
  !
  !  The gamma_gamma case needs a different treatment
  !
  if (lgamma_gamma) then
     with_symmetry=1
     comp_irr = 0
     comp_irr(0)=1
     do na=1,nat
        if (has_equivalent(na)==0) then
            do ipol=1,3
               comp_irr(3*(na-1)+ipol)=1
               with_symmetry(3*(na-1)+ipol)=0
            enddo
        endif
     enddo
     if (nasr>0) then
        do ipol=1,3
           comp_irr(3*(nasr-1)+ipol)=0
           with_symmetry(3*(nasr-1)+ipol)=0
        enddo
     endif
     IF (start_irr <= last_irr_eff) THEN
        DO irr=1,start_irr-1
           comp_irr(irr) = 0
        ENDDO
        DO irr=last_irr_eff+1,3*nat
           comp_irr(irr) = 0
        ENDDO
     ENDIF
  endif
!
!  In this case the number of irreducible representations to compute
!  has been done elsewhere
!
  IF (nimage > 1 ) THEN
     DO irr=0,nirr
        comp_irr(irr)=comp_irr_iq(irr,current_iq)
     ENDDO
  ENDIF

  !
  !  Compute how many atoms moves and set the list atomo
  !
  ifat = 0
  imode0 = 0
  DO irr = 1, nirr
     if (comp_irr (irr) .eq.1) then
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
  all_comp = nat_todo.eq.nat.or.lgamma_gamma
  all_done = .FALSE.
  npertx = 0
  done_irr = 0
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
  CALL ph_writefile('data',0)

  deallocate(w2)
  CALL stop_clock ('phq_setup')
  RETURN
END SUBROUTINE phq_setup
