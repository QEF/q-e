!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
subroutine d3_setup  
  !-----------------------------------------------------------------------
  !
  !  This subroutine prepares several variables which are needed in the
  !  d3toten program:
  !  1) computes the total local potential (external+scf) on the smoot
  !     grid to be used in h_psi and similia
  !  2) computes dmuxc 3.1) with GC if needed
  !  3) for metals sets the occupated bands
  !  4) computes alpha_pv
  !  5.1) computes the variables needed to pass to the pattern representat
  !       of the small group of q
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
  !  5.2) computes the variables needed to pass to the pattern representat
  !       of the group of the crystal
  !     ug0     the patterns
  !     tg0     the matrices of the group on the pattern basis
  !     nsymg0  the order of the group of the crystal
  !     nirrg0  the number of irreducible representation
  !     npertg0 the dimension of each irreducible representation
  !  6) set the variables needed to deal with nlcc
  !  7) set the variables needed to distribute one loop between pools
  !  8) set the variables needed to calculate only selected q=0 modes
  !
#include "machine.h"

  use pwcom
  use phcom
  use d3com
  use allocate
#ifdef PARA
  use para
#endif
  implicit none
  real (8) :: rhotot, rhoup, rhodw, target, small, fac, xmax, emin, &
       emax, dmxc, wrk
  ! total charge
  ! total up charge
  ! total down charge
  ! auxiliary variables used
  ! to set nbnd_occ in the metallic case
  ! minimum band energy
  ! maximum band energy
  ! computes derivative of xc potential
  ! working array

  integer :: ir, table (48, 48), isym, jsym, iinv, irot, jrot, ik, &
       ibnd, ipol, mu, nu, imode0, irr, ipert, nt, ii, nu_i
  ! counter on mesh points
  ! the multiplication table of the point g
  ! counter on symmetries
  ! counter on symmetries
  ! the index of the inverse
  ! counter on rotations
  ! counter on rotations
  ! counter on k points
  ! counter on bands
  ! counter on polarizations
  ! counter on modes
  ! the starting mode
  ! counter on representation and perturbat
  ! counter on atomic type

  logical :: sym (48)  
  ! the symmetry operations

#ifdef PARA
  integer :: nlength_w, nlength (npool), nresto  
#endif
  call start_clock ('d3_setup')  
  !
  ! 1) Computes the total local potential (external+scf) on the smoot grid
  !
  call set_vrs (vrs, vltot, vr, nrxx, nspin, doublegrid)  
  !
  ! 2) Computes the derivative of the xc potential
  !
  call setv (nrxx * nspin * nspin, 0.d0, dmuxc, 1)  
  if (lsda) then  
     do ir = 1, nrxx  
        rhoup = rho (ir, 1) + 0.5d0 * rho_core (ir)  
        rhodw = rho (ir, 2) + 0.5d0 * rho_core (ir)  
        call dmxc_spin (rhoup, rhodw, dmuxc (ir, 1, 1), dmuxc (ir, 2, &
             1), dmuxc (ir, 1, 2), dmuxc (ir, 2, 2) )
     enddo
  else  
     do ir = 1, nrxx  
        rhotot = rho (ir, nspin) + rho_core (ir)  
        if (rhotot.gt.1.d-30) dmuxc (ir, 1, 1) = dmxc (rhotot)  
        if (rhotot.lt. - 1.d-30) dmuxc (ir, 1, 1) = - dmxc ( - rhotot)  
     enddo
  endif
  !
  ! 3) Computes the number of occupated bands for each k point
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
     ! - limit appropriated for gaussian broadening (used for all ngauss)
     !
     xmax = sqrt ( - log (sqrt (pi) * small) )  
     !
     ! - limit appropriated for Fermi-Dirac
     !
     if (ngauss.eq. - 99) then  
        fac = 1.d0 / sqrt (small)  
        xmax = 2.d0 * log (0.5 * (fac + sqrt (fac * fac - 4.0) ) )  

     endif
     target = ef + xmax * degauss  
     do ik = 1, nks  
        do ibnd = 1, nbnd  
           if (et (ibnd, ik) .lt.target) nbnd_occ (ik) = ibnd  
        enddo
        if (nbnd_occ (ik) .eq.nbnd) &
             write (6, '(5x,/,"Possibly too few bands at point ", &
             & i4,3f10.5)') ik,  (xk (ipol, ik) , ipol = 1, 3)
     enddo
  else  
     if (lsda) call error ('d3_setup', 'occupation numbers probably wro &
          &ng',  - 1)
     do ik = 1, nks  
        nbnd_occ (ik) = nint (nelec) / degspin  
     enddo
  endif
  !
  ! 4) Computes alpha_pv
  !
  emin = et (1, 1)  
  do ik = 1, nks  
     do ibnd = 1, nbnd  
        emin = min (emin, et (ibnd, ik) )  
     enddo
  enddo
  ! find the minimum across pools

  call poolextreme (emin, - 1)  
  emax = et (1, 1)  
  do ik = 1, nks  
     do ibnd = 1, nbnd  
        emax = max (emax, et (ibnd, ik) )  
     enddo
  enddo
  ! find the maximum across pools

  call poolextreme (emax, + 1)  
  alpha_pv = 2.d0 * (emax - emin)  
  ! avoid zero value for alpha_pv
  alpha_pv = max (alpha_pv, 1.0d-2)  
  !
  ! 5) set all the variables needed to use the pattern representation
  !
  !
  ! 5.0) Computes the inverse of each matrix
  !
  call multable (nsym, s, table)  
  do isym = 1, nsym  
     do jsym = 1, nsym  
        if (table (isym, jsym) .eq.1) invs (isym) = jsym  
     enddo
  enddo
  !
  ! 5.1) Finds the variables needeed for the pattern representation
  !      of the small group of q
  !
  do isym = 1, nsym  
     sym (isym) = .true.  
  enddo

  call sgam_ph (at, bg, nsym, s, irt, tau, rtau, nat, sym)  
  nmodes = 3 * nat  
  ! if minus_q=.t. set_irr will search for
  minus_q = (iswitch.gt. - 3)  
  ! Sq=-q+G symmetry. On output minus_q=.t.
  ! if such a symmetry has been found
  if (iswitch.eq. - 4) then  
     call set_irr_mode (nat, at, bg, xq, s, invs, nsym, rtau, irt, &
          irgq, nsymq, minus_q, irotmq, t, tmq, u, npert, nirr, gi, gimq, &
          iverbosity, modenum)
  else  
     if (nsym.gt.1) then  
        call io_pattern(fildrho,nirr,npert,u,-1)
        call set_sym_irr (nat, at, bg, xq, s, invs, nsym, rtau, irt, &
             irgq, nsymq, minus_q, irotmq, t, tmq, u, npert, nirr, gi, &
             gimq, iverbosity)
     else  
        call set_irr_nosym (nat, at, bg, xq, s, invs, nsym, rtau, &
             irt, irgq, nsymq, minus_q, irotmq, t, tmq, u, npert, nirr, &
             gi, gimq, iverbosity)
     endif
  endif
  !
  ! 5.2) Finds the variables needeed for the pattern representation
  !      of the small group of the crystal
  !
  if (lgamma) then  
     nsymg0 = nsymq  
     nirrg0 = nirr  
  else  
     !
     ! It finds which symmetries of the lattice are symmetries of the crystal
     ! it calculates the order of the crystal group:   nsymg0
     ! and reorder the s matrices in this way:
     !  a) the first nsymg0 matrices are symmetries of the crystal
     !  b) the first nsymq matrices are symmetries for the small group of q
     !
     call sgama_d3 (nsymq, nat, s, ityp, nr1, nr2, nr3, nsymg0, irt, &
          ftau, at, bg, tau)
     !
     ! Recalculates the inverse of each rotation
     !
     call multable (nsymg0, s, table)  
     do irot = 1, nsymg0  
        do jrot = 1, nsymg0  
           if (table (irot, jrot) .eq.1) invs (irot) = jrot  
        enddo
     enddo
     !
     ! Calculates rtau
     !
     do isym = 1, nsymg0  
        sym (isym) = .true.  
     enddo
     call sgam_ph (at, bg, nsymg0, s, irt, tau, rtau, nat, sym)  
     !
     ! Calculates the variables need for the pattern representation
     ! for the q=0 symmetries
     !
     call set_d3irr  
  endif
  !
  ! 6) Set non linear core correction stuff
  !
  nlcc_any = .false.  
  do nt = 1, ntyp  
     nlcc_any = nlcc_any.or.nlcc (nt)  

  enddo

  if (nlcc_any) call mallocate(drc, ngm, ntyp)  
  !
  ! 7) Sets up variables needed to distribute one loop between pools
  !
  npert_i = 1  
  npert_f = 3 * nat  
#ifdef PARA
  nlength_w = (3 * nat) / npool  
  nresto = 3 * nat - nlength_w * npool  
  do ii = 1, npool  
     if (ii.le.nresto) then  
        nlength (ii) = nlength_w + 1  
     else  
        nlength (ii) = nlength_w  
     endif
  enddo
  npert_i = 1  
  do ii = 1, mypool - 1  
     npert_i = npert_i + nlength (ii)  
  enddo

  npert_f = npert_i - 1 + nlength (mypool)  
#endif
  !
  ! 8) Sets up variables needed to calculate only selected
  !    modes at q=0 --the first index of the third order matrix--
  !
  if (q0mode_todo (1) .le.0) then  
     do ii = 1, 3 * nat  
        q0mode (ii) = .true.  
     enddo
  else  
     do ii = 1, 3 * nat  
        q0mode (ii) = .false.  
     enddo
     ii = 1  
     do while (q0mode_todo (ii) .gt.0)  
        q0mode (q0mode_todo (ii) ) = .true.  
        ii = ii + 1  
     enddo
  endif
  !
  ! if you want to compute all the modes; and lgamma=.true.
  ! the calculation can be simplyfied, in this case allmodes
  ! is set .true.
  !
  allmodes = lgamma.and.q0mode_todo (1) .le.0  
  !
  ! Sets up variables needed to write only selected
  ! modes at q=0 --the first index of the third order matrix--
  !
  do ii = 1, 3 * nat  
     wrk = 0.d0  
     do nu_i = 1, 3 * nat  
        if (q0mode (nu_i) ) then  
           wrk = wrk + ug0 (ii, nu_i) * conjg (ug0 (ii, nu_i) )  
        endif
     enddo
     wrmode (ii) = .false.  
     if (wrk.gt.1.d-8) wrmode (ii) = .true.  

  enddo
  call stop_clock ('d3_setup')  
  return  
end subroutine d3_setup
