!
! Copyright (C) 2001-2007 Quantum-Espresso group
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
  !  Revised 1 Oct. 1995 by Andrea Dal Corso
  !  March   1997: nlcc stuff added (SdG)
  !  April   1997: parallel stuff added (SdG)
  !  Oct-Nov 1998: minor stuff added (SdG)
  !
#include "f_defs.h"
  !
  USE kinds,         ONLY : DP
  USE ions_base,     ONLY : tau, nat, ntyp => nsp, ityp
  USE cell_base,     ONLY : at, bg  
  USE io_global,     ONLY : stdout
  USE ener,          ONLY : Ef
  USE klist,         ONLY : xk, lgauss, degauss, ngauss, nks, nelec
  USE ktetra,        ONLY : ltetra, tetra
  USE lsda_mod,      ONLY : nspin, lsda, starting_magnetization
  USE scf,           ONLY : v, vrs, vltot, rho, rho_core, kedtau
  USE gvect,         ONLY : nrxx, ngm
  USE gsmooth,       ONLY : doublegrid
  USE symme,         ONLY : nsym, s, ftau, irt, t_rev, time_reversal
  USE uspp_param,    ONLY : upf
  USE spin_orb,      ONLY : domag
  USE constants,     ONLY : degspin, pi
  USE noncollin_module, ONLY : noncolin, m_loc, angle1, angle2, ux
  USE wvfct,         ONLY : nbnd, et
  USE rap_point_group,      ONLY : code_group, nclass, nelem, elem, which_irr,&
                                  char_mat, name_rap, gname, name_class, ir_ram
  USE rap_point_group_is,   ONLY : code_group_is, gname_is
  use phcom
  USE control_ph,    ONLY : rec_code
  USE ph_restart,    ONLY : ph_writefile
  USE control_flags, ONLY : iverbosity, modenum
  USE funct,         ONLY : dmxc, dmxc_spin, dmxc_nc, dft_is_gradient
  USE mp,            ONLY : mp_max, mp_min
  USE mp_global,     ONLY : inter_pool_comm
  USE ramanm,        ONLY : lraman, elop, ramtns, eloptns

  implicit none

  real(DP) :: rhotot, rhoup, rhodw, target, small, fac, xmax, emin, emax
  ! total charge
  ! total up charge
  ! total down charge
  ! auxiliary variables used
  ! to set nbnd_occ in the metallic case
  ! minimum band energy
  ! maximum band energy

  real(DP) :: sr(3,3,48), sr_is(3,3,48)

  integer :: ir, table (48, 48), isym, jsym, irot, ik, ibnd, ipol, &
       mu, nu, imode0, irr, ipert, na, it, nt, is, js, nsym_is, last_irr_eff
  ! counter on mesh points
  ! the multiplication table of the point g
  ! counter on symmetries
  ! counter on symmetries
  ! counter on rotations
  ! counter on k points
  ! counter on bands
  ! counter on polarizations
  ! counter on modes
  ! the starting mode
  ! counter on representation and perturbat
  ! counter on atoms
  ! counter on iterations
  ! counter on atomic type

  real(DP) :: auxdmuxc(4,4)

  logical :: sym (48), is_symmorphic
  ! the symmetry operations
  integer, allocatable :: ifat(:)

  call start_clock ('phq_setup')
  ! 0) A few checks
  !
  IF (dft_is_gradient().and.(lraman.or.elop)) call errore('phq_setup', &
     'third order derivatives not implemented with GGA', 1) 
  !
  !
  ! 1) Computes the total local potential (external+scf) on the smooth grid
  !
  call set_vrs (vrs, vltot, v%of_r, kedtau, v%kin_r, nrxx, nspin, doublegrid)
  !
  ! 2) Set non linear core correction stuff
  !
  nlcc_any = ANY ( upf(1:ntyp)%nlcc )
  if (nlcc_any) allocate (drc( ngm, ntyp))    
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
     END DO
     ux=0.0_DP
     if (dft_is_gradient()) call compute_ux(m_loc,ux,nat)
  ENDIF
  !
  ! 3) Computes the derivative of the xc potential
  !
  dmuxc(:,:,:) = 0.d0
  if (lsda) then
     do ir = 1, nrxx
        rhoup = rho%of_r (ir, 1) + 0.5d0 * rho_core (ir)
        rhodw = rho%of_r (ir, 2) + 0.5d0 * rho_core (ir)
        call dmxc_spin (rhoup, rhodw, dmuxc(ir,1,1), dmuxc(ir,2,1), &
                                      dmuxc(ir,1,2), dmuxc(ir,2,2) )
     enddo
  else
     IF (noncolin.and.domag) THEN
        do ir = 1, nrxx
           rhotot = rho%of_r (ir, 1) + rho_core (ir)
           call dmxc_nc (rhotot, rho%of_r(ir,2), rho%of_r(ir,3), rho%of_r(ir,4), auxdmuxc)
           DO is=1,nspin
              DO js=1,nspin
                 dmuxc(ir,is,js)=auxdmuxc(is,js)
              END DO
           END DO
        enddo
     ELSE
        do ir = 1, nrxx
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
  ! 4) Computes the inverse of each matrix
  !
  call multable (nsym, s, table)
  do isym = 1, nsym
     do jsym = 1, nsym
        if (table (isym, jsym) == 1) invs (isym) = jsym
     enddo
  enddo
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
     if (lsda) call infomsg('phq_setup','occupation numbers probably wrong')
     if (noncolin) then
        nbnd_occ = nint (nelec) 
     else
        do ik = 1, nks
           nbnd_occ (ik) = nint (nelec) / degspin
        enddo
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
  time_reversal = .NOT. (noncolin .AND. domag)
  ! 
  ! allocate and calculate rtau, the rotated position of each atom
  !
  do isym = 1, nsym
     sym (isym) = .true.
  enddo

  call sgam_ph (at, bg, nsym, s, irt, tau, rtau, nat, sym)
  nmodes = 3 * nat
  ! if minus_q=.t. set_irr will search for
  minus_q = (modenum .eq. 0)
  ! Sq=-q+G symmetry. On output minus_q=.t.
  ! if such a symmetry has been found
  if (modenum .ne. 0) then
     call set_irr_mode (nat, at, bg, xq, s, invs, nsym, rtau, irt, &
          irgq, nsymq, minus_q, irotmq, t, tmq, max_irr_dim, u, npert, &
          nirr, gi, gimq, iverbosity, modenum)
  else
     if (nsym > 1.and..not.lgamma_gamma) then
        call set_irr (nat, at, bg, xq, s, invs, nsym, rtau, irt, &
             irgq, nsymq, minus_q, irotmq, t, tmq, max_irr_dim, u, npert, &
             nirr, gi, gimq, iverbosity,rec_code)
     else
        call set_irr_nosym (nat, at, bg, xq, s, invs, nsym, rtau, irt, &
             irgq, nsymq, minus_q, irotmq, t, tmq, max_irr_dim, u, npert, &
             nirr, gi, gimq, iverbosity)
     endif
  endif
  is_symmorphic=.true.
  DO isym=1,nsymq
     is_symmorphic=( is_symmorphic.and.(ftau(1,irgq(isym))==0).and.  &
                                       (ftau(2,irgq(isym))==0).and.  &
                                       (ftau(3,irgq(isym))==0) )
  
  END DO
  search_sym=.true.
  IF (.not.is_symmorphic) THEN
     DO isym=1,nsymq
        search_sym=( search_sym.and.(abs(gi(1,irgq(isym)))<1.d-8).and.  &
                                    (abs(gi(2,irgq(isym)))<1.d-8).and.  &
                                    (abs(gi(3,irgq(isym)))<1.d-8) )
     END DO
  END IF
  IF (search_sym) THEN
     DO isym=1,nsym
        CALL s_axis_to_cart (s(1,1,isym), sr(1,1,isym), at, bg)
     END DO
     CALL find_group(nsym,sr,gname,code_group)
     CALL set_irr_rap(code_group,nclass,char_mat,name_rap,name_class,ir_ram)
     CALL divide_class(code_group,nsym,sr,nclass,nelem,elem,which_irr)
     IF (noncolin .and. domag) THEN
        nsym_is=0.d0
        DO isym=1,nsym
           IF (t_rev(isym)==0) THEN
              nsym_is=nsym_is+1
              CALL s_axis_to_cart (s(1,1,isym), sr_is(1,1,nsym_is), at, bg)
           ENDIF
        END DO
        CALL find_group(nsym_is,sr_is,gname_is,code_group_is)
     ENDIF
  ENDIF

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


  if (fildrho.ne.' ') call io_pattern (fildrho,nirr,npert,u,+1)

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
  IF (nrapp==0 .AND. nat_todo==0) THEN
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
  ELSEIF (nrapp /= 0) THEN
     !
     !   Case 3) The representation which must be computed are given
     !           as input
     !
     if (nrapp > nirr) call errore ('phq_setup', 'too many representations', 1)
     do irr = 1, nirr
        do mu = 1, nrapp
           if (list (mu) == irr) comp_irr (irr) = 1
        enddo
     enddo
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
  IF (epsil.and.lgamma) THEN
     epsilon=0.0_DP
     zstareu=0.0_DP
  ENDIF
  IF (zue.and.lgamma) zstarue0=(0.0_DP,0.0_DP)
  IF (lraman.and.lgamma) ramtns=0.0_DP
  IF (elop.and.lgamma)  eloptns=0.0_DP
!
!  if this was not a recover run, save the displacement pattern found here 
!
  IF (.NOT. recover) THEN
     where_rec='setup'
     rec_code=0
     CALL ph_writefile('data',0)
  ENDIF

  CALL stop_clock ('phq_setup')
  RETURN
END SUBROUTINE phq_setup
