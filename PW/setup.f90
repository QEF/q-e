!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------

subroutine setup
  !-----------------------------------------------------------------------
  !
  !  This routine computes the main quantities necessary to describe a
  !  periodic solid. In particular it sets the following quantities:
  !  msh       the point which correspond to a cut-off of 10 a.u.
  !  ecutwfc   the cut-off of each iteration
  !  zv        the charge of each atomtype
  !  nelec     the total number of electrons
  !  nbnd      the total number of bands
  !  nbndx     the maximum number of bands
  !  alat      the measure of length
  !  tpiba     the measure unit in reciprocal space
  !  tpiba2    the square of tpiba
  !  gcutm     the cut-off in g space
  !  gcutms    the cut-off in g space for smooth functions
  !  ethr      the convergence limit of iterative diagonalization
  !  at        the direct lattice vectors
  !  omega     the volume of the unit cell
  !  bg        the reciprocal lattice vectors
  !  g, gg     the reciprocal g vectors
  !  ngm       the number of g vectors in reciprocal space
  !  nl        the correspondence between g vectors <--> fft mesh
  !  s         the symmetry matrices in the direct lattice vectors basis
  !  nsym      the total number of symmetry operations
  !  ftau      the fractionary translations
  !  irt       for each atom gives the corresponding symmetric
  !  invsym    if true the system has inversion symmetry
  !
  !  Furthermore it transforms the  BHS parameters to a working form.
  !
  !
  !
#include "machine.h"
  use pwcom
#ifdef PARA
  use para
#endif
  implicit none
  !
  !    local variables
  !
  !    rcut: cut-off radius for radial integrations
  !
  real(kind=DP), parameter :: rcut = 10.d0, eps = 1.0d-12
  !
  integer :: na, ir, nt, input_nks, nrot, iter, ierr, irot, isym, &
       ipol, jpol, tipo
  integer, external ::n_atom_wfc
  !
  ! lchk_tauxk tests that atomic coordinates do not overlap
  !
  logical, external :: lchk_tauxk
  logical :: minus_q

  do nt = 1, ntyp
     do ir = 1, mesh (nt)
        if (r (ir, nt) .gt.rcut) then
           msh (nt) = ir
           goto 5
        endif
     enddo

     msh (nt) = mesh (nt)
     !
     ! force msh to be odd for simpson integration
     !
5    msh (nt) = 2 * ( (msh (nt) + 1) / 2) - 1
  enddo
  !
  !     Compute the ionic charge for each atom type
  !
  do nt = 1, ntyp
     zv (nt) = zp (nt)
  enddo
  !
  !     Set the number of electrons equal to the total ionic charge
  !
  if (nelec.eq.0.d0) then
     do na = 1, nat
        nelec = nelec + zv (ityp (na) )
     enddo
  endif
  !
  ! For metals: check whether Gaussian broadening or Tetrahedron method is
  !
  lgauss = degauss.ne.0.d0
  !
  ! Check: if there is an odd number of electrons, the crystal is a metal
  !
  if (lscf.and.abs (nint (nelec / 2.d0) - nelec / 2.d0) &
       .gt.1.0d-8.and..not.lgauss.and..not.ltetra) call error ('setup', &
       'the crystal is a metal', 1)
  !
  !     Set the number of occupied bands if not given in input
  !
  if (nbnd.eq.0) then
     nbnd = nint (nelec) / 2.d0
     if (lgauss.or.ltetra) then
        !
        !     metallic case: add 20% more bands, with a minimum of 4
        !
        nbnd = max (nint (1.20d0 * nelec / 2.d0), nbnd+4)
     endif
  else
     if (nbnd.lt.nint (nelec) / 2.d0.and.lscf) call error ('setup', &
          'too few bands', 1)
  endif
  !
  ! Here we  set the precision of the diagonalization
  !
  if (startingpot=='file') then
     ! starting potential is expected to be a good one :
     ! do not spoil it with a lousy first diagonalization
     if (imix.ge.0) ethr = 0.1d0 * min (1.d-2, tr2/nelec)
     if (imix.lt.0) ethr = 0.1d0 * min (1.0d-6, sqrt (tr2) )
  else
     ! starting atomic potential is probably far from scf..
     ! do not waste iterations in the first diagonalizations
     if (imix.ge.0) ethr = 1.0d-2
     if (imix.lt.0) ethr = 1.0d-5
  end if

  if (.not.lscf) niter = 1

  starting_diag_threshold = ethr
  !
  !    set number of atomic wavefunctions
  !
  natomwfc = n_atom_wfc (nat, npsx, ityp, newpseudo, nchix, nchi, oc, lchi)
  !
  !    set the maximum number of bands
  !
  nbndx = nbnd
  if (isolve.eq.0) then
     nbndx = david * nbnd
  elseif (isolve.eq.2) then
     if (diis_start_dav.gt.0) then
        nbndx = david * nbnd
     else
        if (diis_wfc_keep) then
           nbndx = 2 * nbnd
        else
           nbndx = nbnd
        endif
     endif
  endif

  if (startingwfc.eq.'atomic'.and..not.restart) nbndx = max (natomwfc, nbndx)
  !
  !     Set the units in real and reciprocal space
  !
  tpiba = 2.d0 * pi / alat
  tpiba2 = tpiba**2
  !
  !     Compute the cut-off of the G vectors
  !
  gcutm = dual * ecutwfc / tpiba2
  doublegrid = dual.gt.4.0d0
  if (doublegrid) then
     gcutms = 4.d0 * ecutwfc / tpiba2
  else
     gcutms = gcutm
  endif
  !
  !     Generate the reciprocal lattice vectors
  !

  call recips ( at(1,1), at(1,2), at(1,3), bg(1,1), bg(1,2), bg(1,3) )
  !
  !     If  lxkcry = .true. , the input k-point components in crystallogra
  !     units are transformed in cartesian coordinates
  !
  if (lxkcry) call cryst_to_cart (nks, xk, bg, 1)
  !
  !     Test that atomic coordinates are different
  !
  if (.not. (lchk_tauxk (nat, tau, bg) ) ) call error ('setup', &
       'Wrong atomic coordinates ', 1)
  !
  ! set dtau_ref for relaxation and dynamics
  ! this is done here because dtau_ref is updated in cg
  !
  dtau_ref = 0.2d0
  !
  !   calculate dimensions of the FFT grid
  !
  call set_fft_dim
  !
  !  generate transformation matrices for the crystal point group
  !  First we generate all the symmetry matrices of the Bravais lattice
  !
  if (ibrav.eq.4.or.ibrav.eq.5) then
     !
     !  here the hexagonal or trigonal bravais lattice
     !
     call hexsym (at, s, sname, nrot)
     tipo = 2
  elseif (ibrav.ge.1.and.ibrav.le.14) then
     !
     !  here for the cubic bravais lattice
     !
     call cubicsym (at, s, sname, nrot)
     tipo = 1
  elseif (ibrav.eq.0) then
     if (symm_type.eq.'cubic') then
        tipo = 1
        call cubicsym (at, s, sname, nrot)
     endif
     if (symm_type.eq.'hexagonal') then
        tipo = 2
        call hexsym (at, s, sname, nrot)
     endif
  else
     call error ('setup', 'wrong ibrav', 1)
  endif
  !
  !   if noinv is true eliminate all symmetries which exchange z with -z
  !
  if (noinv) then
     irot = 0
     do isym = 1, nrot
        if (s (1, 3, isym) .eq.0.and.s (3, 1, isym) .eq.0.and.s (2, 3, &
             isym) .eq.0.and.s (3, 2, isym) .eq.0.and.s (3, 3, isym) .eq.1) &
             then
           irot = irot + 1
           do ipol = 1, 3
              do jpol = 1, 3
                 s (ipol, jpol, irot) = s (ipol, jpol, isym)
              enddo
           enddo
           sname (irot) = sname (isym)
        endif
     enddo
     nrot = irot
  endif

  !
  !   If nosym is true do not use any point-group symmetry
  !

  if (nosym) nrot = 1

  !
  !    Automatic generation of k-points (if required)
  !

  if (nks.lt.0) then
     call setupkpoint (s, nrot, xk, wk, nks, npk, nk1, &
          nk2, nk3, k1, k2, k3, at, bg, tipo)
  else if (nks.eq.0) then
     call kpoint_grid ( nrot, s, bg, npk, k1,k2,k3, nk1,nk2,nk3, nks, xk, wk)
  end if
  !
  input_nks = nks
  !
  !   allocate space for irt
  !
  allocate (irt( 48, nat))    
  !
  !   "sgama" eliminates rotations that are not symmetry operations
  !   Input k-points are assumed to be  given in the IBZ of the Bravais
  !   lattice, with the full point symmetry of the lattice.
  !   If some symmetries are missing in the crystal, "sgama" computes
  !   the missing k-points. If nosym is true (see above) we do not use
  !   any point-group symmetry and leave k-points unchanged.
  !
  !
  call sgama (nrot, nat, s, sname, at, bg, tau, ityp, nsym, nr1, &
       nr2, nr3, irt, ftau, npk, nks, xk, wk, invsym, minus_q, xqq, &
       iswitch, modenum)

  call checkallsym (nsym, s, nat, tau, ityp, at, bg, nr1, nr2, nr3, &
       irt, ftau)
  !
  ! if dynamics is done the system should have no symmetries
  ! (inversion symmetry alone is allowed)
  !

  if (iswitch.gt.2.and. (nsym.eq.2.and..not.invsym.or.nsym.gt.2) &
     .and. .not. ( calc.eq.'mm' .or. calc.eq.'nm' ) ) &
     call error ('setup', 'Dynamics, you should have no symmetries', -1)
  !
  !     Automatic generation of k-points
  !

  if (ltetra) then
     ntetra = 6 * nk1 * nk2 * nk3
     allocate (tetra(4,ntetra))    
     call tetrahedra (nsym,s,minus_q,at,bg,npk,k1,k2,k3, &
          nk1,nk2,nk3, nks,xk,wk,ntetra,tetra)
  end if
  !
  !     In a non scf calculation we check if k-points are compatible
  !     with the symmetry of the crystal but do not change them
  !
  if ( iswitch.eq.-1 .and. input_nks.ne.nks .and. .not.ltetra ) then
     write (6, '(5x,"Warning: missing inequivalent points")')
     nks = input_nks
  endif

  if (iswitch.le.-2) call set_kplusq (xk, wk, xqq, nks, npk)
  if (lsda) then
     !
     ! LSDA case: two different spin polarizations, each with its own kpoints
     !
     nspin = 2
     call set_kup_and_kdw (xk, wk, isk, nks, npk)
  else
     !
     ! LDA case: the two spin polarizations are identical
     !
     nspin = 1
     current_spin = 1
  endif
  if (nks.gt.npk) call error ('setup', 'too many k points', nks)
#ifdef PARA
  call init_pool
  !
  ! set the granularity for k-point distribution
  !
  if ( (abs (xqq (1) ) .lt.eps.and.abs (xqq (2) ) .lt.eps.and.abs ( &
       xqq (3) ) .lt.eps) .or. (iswitch.gt. - 2) ) then
     kunit = 1
  else
     kunit = 2
  endif
  !
  ! distribute the k-points (and their weights and spin indices)
  !

  call divide_et_impera (xk, wk, isk, lsda, nkstot, nks)
#else
  !
  !  set nkstot which is used to write results for all k-points
  !

  nkstot = nks
#endif
  !
  !   allocate memory for G- and R-space fft arrays
  !
  call allocate_fft
  !
  !   generate reciprocal-lattice vectors and fft indices
  !

  call ggen
#ifndef PARA
  !
  !   generates pencils for 3d-fft of psi and related quantities
  !
  call set_pencils (nks, xk, ngms, gg, nls, ecutwfc / tpiba2, &
       nrx1s, nr1s, nr2s, nr3s)
#endif
  !
  !   For Bachelet-Hamann-Schluter pseudopotentials only
  !
  do nt = 1, ntyp
     if (.not.tvanp (nt) ) then
        if (nlc (nt) .eq.2.and.nnl (nt) .eq.3.and.bhstype (nt) ) call &
             bachel (alps (1, 0, nt), aps (1, 0, nt), 1, lmax (nt) )
     endif
  enddo
  !
  !   okvan = .true. : at least one pseudopotential is US
  !
  okvan = .false.
  do nt = 1, ntyp
     okvan = okvan.or.tvanp (nt)
  enddo
  !
  ! initialize parameters for charge density extrapolation during dynamics
  !
  alpha0 = 1.d0
  beta0 = 0.d0
  !
  ! Needed for LDA+U
  !
  ! initialize d1 and d2 to rotate the spherical harmonics
  !

  if (lda_plus_u) call d_matrix (d1, d2)
  !
  return
end subroutine setup
!
!-----------------------------------------------------------------------

integer function n_atom_wfc (nat, npsx, ityp, newpseudo, nchix, nchi, oc, lchi)
  !-----------------------------------------------------------------------
  !
  ! Find max number of bands needed
  !
  use parameters, only : DP
  implicit none
  integer :: nat, npsx, ityp (nat), nchix, nchi (npsx), lchi (nchix, npsx)
  real(kind=DP) :: oc (nchix, npsx)
  logical :: newpseudo (npsx)

  integer :: na, nt, n
  n_atom_wfc = 0
  do na = 1, nat
     nt = ityp (na)
     do n = 1, nchi (nt)
        if (oc (n, nt) .gt.0.d0.or..not.newpseudo (nt) ) n_atom_wfc = &
             n_atom_wfc + 2 * lchi (n, nt) + 1
     enddo
  enddo

  return
end function n_atom_wfc
