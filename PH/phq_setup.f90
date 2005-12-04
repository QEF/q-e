!
! Copyright (C) 2001 PWSCF group
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
  !  1) put the correct units on the masses
  !  2) computes the total local potential (external+scf) on the smoot
  !     grid to be used in h_psi and similia
  !  3) computes dmuxc 3.1) with GC if needed
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
  USE ions_base,     ONLY : tau, nat, ntyp => nsp
  USE io_global,     ONLY : stdout
  use pwcom
  USE atom,          ONLY : nlcc
  USE constants,     ONLY : degspin
  USE kinds,         ONLY : DP
  use phcom
  USE control_flags, ONLY : iverbosity, modenum
  USE funct,         ONLY : dmxc, dmxc_spin  
  implicit none

  real(DP) :: rhotot, rhoup, rhodw, target, small, fac, xmax, emin, emax
  ! total charge
  ! total up charge
  ! total down charge
  ! auxiliary variables used
  ! to set nbnd_occ in the metallic case
  ! minimum band energy
  ! maximum band energy

  integer :: ir, table (48, 48), isym, jsym, irot, ik, ibnd, ipol, &
       mu, nu, imode0, irr, ipert, na, it, nt
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

  logical :: sym (48)
  ! the symmetry operations

  call start_clock ('phq_setup')
  !
  ! 1) We start with the mass renormalization
  !
!  call DSCAL (ntyp, amconv, amass, 1)
!     (now done in phq_readin.f90 to avoid to be done more than once)
!
  !
  ! 2) Computes the total local potential (external+scf) on the smooth grid
  !
  call set_vrs (vrs, vltot, vr, nrxx, nspin, doublegrid)
  !
  ! 2.a) Set non linear core correction stuff
  !
  nlcc_any = .false.
  do nt = 1, ntyp
     nlcc_any = nlcc_any.or.nlcc (nt)

  enddo
  if (nlcc_any) allocate (drc( ngm, ntyp))    
  !
  ! 3) Computes the derivative of the xc potential
  !
  dmuxc(:,:,:) = 0.d0
  if (lsda) then
     do ir = 1, nrxx
        rhoup = rho (ir, 1) + 0.5d0 * rho_core (ir)
        rhodw = rho (ir, 2) + 0.5d0 * rho_core (ir)
        call dmxc_spin (rhoup, rhodw, dmuxc(ir,1,1), dmuxc(ir,2,1), &
                                      dmuxc(ir,1,2), dmuxc(ir,2,2) )
     enddo
  else
     do ir = 1, nrxx
        rhotot = rho (ir, nspin) + rho_core (ir)
        if (rhotot.gt.1.d-30) dmuxc (ir, 1, 1) = dmxc (rhotot)
        if (rhotot.lt. - 1.d-30) dmuxc (ir, 1, 1) = - dmxc ( - rhotot)
     enddo
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
        if (table (isym, jsym) .eq.1) invs (isym) = jsym
     enddo
  enddo
  !
  ! 5) Computes the number of occupied bands for each k point
  !
  if (degauss.ne.0.d0) then
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
     if (lsda) call infomsg('phq_setup','occupation numbers probably wrong', -1)
     do ik = 1, nks
        nbnd_occ (ik) = nint (nelec) / degspin
     enddo
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
  call poolextreme (emin, -1)
#endif
  if (degauss.ne.0.d0) then
     emax = target
     alpha_pv = emax - emin
  else
     emax = et (1, 1)
     do ik = 1, nks
        do ibnd = 1, nbnd
           emax = max (emax, et (ibnd, ik) )
        enddo
     enddo
#ifdef __PARA
     ! find the maximum across pools
     call poolextreme (emax, + 1)
#endif
     alpha_pv = 2.d0 * (emax - emin)
  endif
  ! avoid zero value for alpha_pv
  alpha_pv = max (alpha_pv, 1.0d-2)
  !
  ! 7) set all the variables needed to use the pattern representation
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
     if (nsym.gt.1) then
        call set_irr (nat, at, bg, xq, s, invs, nsym, rtau, irt, &
             irgq, nsymq, minus_q, irotmq, t, tmq, max_irr_dim, u, npert, &
             nirr, gi, gimq, iverbosity)
     else
        call set_irr_nosym (nat, at, bg, xq, s, invs, nsym, rtau, irt, &
             irgq, nsymq, minus_q, irotmq, t, tmq, max_irr_dim, u, npert, &
             nirr, gi, gimq, iverbosity)
     endif
  endif

  if (fildrho.ne.' ') call io_pattern (fildrho,nirr,npert,u,+1)

  !
  !  set maxirr if not already set
  !
  if (maxirr.le.0.or.maxirr.gt.nirr) maxirr = nirr + 1
  if (niter_ph.lt.maxter) maxirr = 1
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
  !  10) set the variables needed for the partial computation
  !
  if (nrapp.eq.0) then
     if (nat_todo.eq.0) then
        !
        !    The partial computation option is not used, compute all atoms
        !
        do na = 1, nat
           atomo (na) = na
        enddo
        nat_todo = nat
     endif
     !
     !   Sets the atoms which must be computed: the requested atoms and all
     !   the symmetry related atoms
     !
     do na = 1, nat
        ifat (na) = 0
     enddo
     do na = 1, nat_todo
        ifat (atomo (na) ) = 1
        do isym = 1, nsymq
           irot = irgq (isym)
           ifat (irt (irot, atomo (na) ) ) = 1
        enddo
     enddo
     !
     !    Computes again nat_todo, prepare the list atomo and sets all_comp
     !
     nat_todo = 0
     do na = 1, nat
        if (ifat (na) .eq.1) then
           nat_todo = nat_todo + 1
           atomo (nat_todo) = na
        endif
     enddo
     !
     !     Find the irreducible representations to be computed
     !
     imode0 = 0
     do irr = 1, nirr
        comp_irr (irr) = 0
        do ipert = 1, npert (irr)
           mu = imode0 + ipert
           do na = 1, nat
              if (ifat (na) .eq.1.and.comp_irr (irr) .eq.0) then
                 do ipol = 1, 3
                    nu = 3 * (na - 1) + ipol
                    if (abs (u (nu, mu) ) .gt.1.d-6) comp_irr (irr) = 1
                 enddo
              endif
           enddo
        enddo
        imode0 = imode0 + npert (irr)
     enddo
  else
     if (nrapp.gt.nirr) call errore ('phq_setup', 'too many representati &
          &on', 1)
     do irr = 1, nirr
        comp_irr (irr) = 0
        do mu = 1, nrapp
           if (list (mu) .eq.irr) comp_irr (irr) = 1
        enddo
     enddo
     do na = 1, nat
        ifat (na) = 0
     enddo
     imode0 = 0
     do irr = 1, nirr
        if (comp_irr (irr) .eq.1) then
           do ipert = 1, npert (irr)
              do na = 1, nat
                 do ipol = 1, 3
                    mu = 3 * (na - 1) + ipol
                    if (abs (u (mu, imode0 + ipert) ) .gt.1.d-12) ifat (na) &
                         = 1
                 enddo
              enddo
           enddo
        endif
        imode0 = imode0 + npert (irr)
     enddo
     nat_todo = 0
     do na = 1, nat
        if (ifat (na) .eq.1) then
           nat_todo = nat_todo + 1
           atomo (nat_todo) = na
        endif
     enddo
  endif
  !
  !   Initialize done_irr, find max dimension of the irreps
  !
  all_comp = nat_todo.eq.nat
  npertx = 0
  do irr = 1, nirr
     done_irr (irr) = 0
     npertx = max (npertx, npert (irr) )

  enddo
  call stop_clock ('phq_setup')
  return
end subroutine phq_setup
