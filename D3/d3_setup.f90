!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "f_defs.h"
!
!-----------------------------------------------------------------------
SUBROUTINE d3_setup()
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
  USE ions_base,     ONLY : nat, ityp, ntyp => nsp, tau
  USE io_global,     ONLY : stdout
  USE kinds,         ONLY : DP
  USE pwcom
  USE scf, only : rho, rho_core, v, vltot, vrs, kedtau
  USE uspp_param,    ONLY : upf
  USE control_flags, ONLY : iverbosity, modenum
  USE constants,     ONLY : degspin
  USE phcom
  USE d3com
  USE mp_global,     ONLY : npool, me_pool
  USE funct,         ONLY : dmxc, dmxc_spin
  !
  IMPLICIT NONE
  !
  REAL (DP) :: rhotot, rhoup, rhodw, TARGET, small, fac, xmax, emin, &
       emax, wrk
  ! total charge
  ! total up charge
  ! total down charge
  ! auxiliary variables used
  ! to set nbnd_occ in the metallic case
  ! minimum band energy
  ! maximum band energy
  ! working array

  INTEGER :: table (48, 48)
  ! the multiplication table of the point group
  INTEGER :: ir, isym, jsym, iinv, irot, jrot, ik, &
       ibnd, ipol, mu, nu, imode0, irr, ipert, nt, ii, nu_i
  ! counters
  LOGICAL :: sym (48)
  ! the symmetry operations

#ifdef __PARA
  INTEGER :: nlength_w, nlength (npool), nresto
#endif
  CALL start_clock ('d3_setup')
  !
  ! 1) Computes the total local potential (external+scf) on the smoot grid
  !
  CALL set_vrs (vrs, vltot, v%of_r, kedtau, v%kin_r, nrxx, nspin, doublegrid)
  !
  ! 2) Computes the derivative of the xc potential
  !
  dmuxc (:,:,:) = 0.d0
  IF (lsda) THEN
     DO ir = 1, nrxx
        rhoup = rho%of_r (ir, 1) + 0.5d0 * rho_core (ir)
        rhodw = rho%of_r (ir, 2) + 0.5d0 * rho_core (ir)
        CALL dmxc_spin (rhoup, rhodw, dmuxc (ir, 1, 1), &
             dmuxc (ir, 2, 1), dmuxc (ir, 1, 2), dmuxc (ir, 2, 2) )
     ENDDO
  ELSE
     DO ir = 1, nrxx
        rhotot = rho%of_r (ir, nspin) + rho_core (ir)
        IF (rhotot > 1.d-30) dmuxc (ir, 1, 1) = dmxc (rhotot)
        IF (rhotot < - 1.d-30) dmuxc (ir, 1, 1) = - dmxc ( - rhotot)
     ENDDO
  ENDIF
  !
  ! 3) Computes the number of occupated bands for each k point
  !
  IF (degauss /= 0.d0) THEN
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
     xmax = SQRT ( - LOG (SQRT (pi) * small) )
     !
     ! - limit appropriated for Fermi-Dirac
     !
     IF (ngauss == - 99) THEN
        fac = 1.d0 / SQRT (small)
        xmax = 2.d0 * LOG (0.5d0 * (fac + SQRT (fac * fac - 4.0d0) ) )
     ENDIF
     TARGET = ef + xmax * degauss
     DO ik = 1, nks
        DO ibnd = 1, nbnd
           IF (et (ibnd, ik) < TARGET) nbnd_occ (ik) = ibnd
        ENDDO
        IF (nbnd_occ (ik) == nbnd) &
             WRITE( stdout, '(5x,/,"Possibly too few bands at point ", &
             & i4,3f10.5)') ik,  (xk (ipol, ik) , ipol = 1, 3)
     ENDDO
  ELSE
     IF (lsda) CALL infomsg ('d3_setup', 'occupation numbers probably wrong')
     DO ik = 1, nks
        nbnd_occ (ik) = NINT (nelec) / degspin
     ENDDO
  ENDIF
  !
  ! 4) Computes alpha_pv
  !
  emin = et (1, 1)
  DO ik = 1, nks
     DO ibnd = 1, nbnd
        emin = MIN (emin, et (ibnd, ik) )
     ENDDO
  ENDDO
  ! find the minimum across pools

  CALL poolextreme (emin, - 1)
  emax = et (1, 1)
  DO ik = 1, nks
     DO ibnd = 1, nbnd
        emax = MAX (emax, et (ibnd, ik) )
     ENDDO
  ENDDO
  ! find the maximum across pools

  CALL poolextreme (emax, + 1)
  alpha_pv = 2.d0 * (emax - emin)
  ! avoid zero value for alpha_pv
  alpha_pv = MAX (alpha_pv, 1.0d-2)
  !
  ! 5) set all the variables needed to use the pattern representation
  !
  !
  ! 5.0) Computes the inverse of each matrix
  !
  CALL multable (nsym, s, table)
  DO isym = 1, nsym
     DO jsym = 1, nsym
        IF (table (isym, jsym) .EQ.1) invs (isym) = jsym
     ENDDO
  ENDDO
  !
  ! 5.1) Finds the variables needeed for the pattern representation
  !      of the small group of q
  !
  DO isym = 1, nsym
     sym (isym) = .TRUE.
  ENDDO

  CALL sgam_ph (at, bg, nsym, s, irt, tau, rtau, nat, sym)
  nmodes = 3 * nat
  ! if minus_q=.t. set_irr will search for
  minus_q = (modenum .eq. 0)
  ! Sq=-q+G symmetry. On output minus_q=.t.
  ! if such a symmetry has been found
  IF (modenum .ne. 0) THEN
     CALL set_irr_mode (nat, at, bg, xq, s, invs, nsym, rtau, irt, &
          irgq, nsymq, minus_q, irotmq, t, tmq, max_irr_dim, u,    &
          npert, nirr, gi, gimq, iverbosity, modenum)
  ELSE
     IF (nsym > 1) THEN
        CALL io_pattern(fildrho,nirr,npert,u,-1)
        CALL set_sym_irr (nat, at, bg, xq, s, invs, nsym, rtau, irt, &
             irgq, nsymq, minus_q, irotmq, t, tmq, max_irr_dim, u,   &
             npert, nirr, gi, gimq, iverbosity)
     ELSE
        CALL set_irr_nosym (nat, at, bg, xq, s, invs, nsym, rtau, &
             irt, irgq, nsymq, minus_q, irotmq, t, tmq, max_irr_dim, u, & 
             npert, nirr, gi, gimq, iverbosity)
     ENDIF
  ENDIF
  !
  ! 5.2) Finds the variables needeed for the pattern representation
  !      of the small group of the crystal
  !
  IF (lgamma) THEN
     nsymg0 = nsymq
     nirrg0 = nirr
  ELSE
     !
     ! It finds which symmetries of the lattice are symmetries of the crystal
     ! it calculates the order of the crystal group:   nsymg0
     ! and reorder the s matrices in this way:
     !  a) the first nsymg0 matrices are symmetries of the crystal
     !  b) the first nsymq matrices are symmetries for the small group of q
     !
     CALL sgama_d3 (nsymq, nat, s, ityp, nr1, nr2, nr3, nsymg0, irt, &
          ftau, at, bg, tau)
     !
     ! Recalculates the inverse of each rotation
     !
     CALL multable (nsymg0, s, table)
     DO irot = 1, nsymg0
        DO jrot = 1, nsymg0
           IF (table (irot, jrot) == 1) invs (irot) = jrot
        ENDDO
     ENDDO
     !
     ! Calculates rtau
     !
     DO isym = 1, nsymg0
        sym (isym) = .TRUE.
     ENDDO
     CALL sgam_ph (at, bg, nsymg0, s, irt, tau, rtau, nat, sym)
     !
     ! Calculates the variables need for the pattern representation
     ! for the q=0 symmetries
     !
     CALL set_d3irr
  ENDIF
  !
  npertx = 0
  do irr = 1, nirr
      npertx = max (npertx, npert (irr) )
  enddo
  do irr = 1, nirrg0
      npertx = max (npertx, npertg0 (irr) )
  enddo
  !
  ! 6) Set non linear core correction stuff
  !
  nlcc_any = ANY ( upf(1:ntyp)%nlcc )
  !
  IF (nlcc_any) ALLOCATE (drc( ngm, ntyp))    
  !
  ! 7) Sets up variables needed to distribute one loop between pools
  !
  npert_i = 1
  npert_f = 3 * nat
#ifdef __PARA
  nlength_w = (3 * nat) / npool
  nresto = 3 * nat - nlength_w * npool
  DO ii = 1, npool
     IF (ii <= nresto) THEN
        nlength (ii) = nlength_w + 1
     ELSE
        nlength (ii) = nlength_w
     ENDIF
  ENDDO
  npert_i = 1
  DO ii = 1, me_pool
     npert_i = npert_i + nlength (ii)
  ENDDO

  npert_f = npert_i - 1 + nlength (me_pool+1)
#endif
  !
  ! 8) Sets up variables needed to calculate only selected
  !    modes at q=0 --the first index of the third order matrix--
  !
  IF (q0mode_todo (1) <= 0) THEN
     DO ii = 1, 3 * nat
        q0mode (ii) = .TRUE.
     ENDDO
  ELSE
     DO ii = 1, 3 * nat
        q0mode (ii) = .FALSE.
     ENDDO
     ii = 1
     DO WHILE (q0mode_todo (ii) > 0)
        q0mode (q0mode_todo (ii) ) = .TRUE.
        ii = ii + 1
     ENDDO
  ENDIF
  !
  ! if you want to compute all the modes; and lgamma=.true.
  ! the calculation can be simplyfied, in this case allmodes
  ! is set .true.
  !
  allmodes = lgamma.AND.q0mode_todo (1) <= 0
  !
  ! Sets up variables needed to write only selected
  ! modes at q=0 --the first index of the third order matrix--
  !
  DO ii = 1, 3 * nat
     wrk = 0.d0
     DO nu_i = 1, 3 * nat
        IF (q0mode (nu_i) ) THEN
           wrk = wrk + ug0 (ii, nu_i) * CONJG (ug0 (ii, nu_i) )
        ENDIF
     ENDDO
     wrmode (ii) = .FALSE.
     IF (wrk > 1.d-8) wrmode (ii) = .TRUE.
  ENDDO
  CALL stop_clock ('d3_setup')
  RETURN
END SUBROUTINE d3_setup
