!
! Copyright (C) 2001-2003 PWSCF group
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
  !  This routine 
  !  1) determines various parameters of the calculation
  !  2) finds actual crystal symmetry, determine lattice
  !  3) generates k-points corresponding to the crystal symmetry
  !  4) transforms the BHS parameters to a working form
  !
  !  Calculated parameters:
  !    msh       mesh point (atomic grid) for which R(msh) > Rcut = 10a.u.
  !    zv        charge of each atomic type
  !    nelec     total number of electrons
  !    nbnd      total number of bands
  !    nbndx     max number of bands used in iterative diagonalization
  !    tpiba     2 pi / a (a = lattice parameter)
  !    tpiba2    square of tpiba
  !    gcutm     cut-off in g space
  !    gcutms    cut-off in g space for smooth functions
  !    ethr      convergence limit of iterative diagonalization
  !    at        direct lattice vectors
  !    omega     volume of the unit cell
  !    bg        reciprocal lattice vectors
  !    s         symmetry matrices in the direct lattice vectors basis
  !    nsym      total number of symmetry operations
  !    ftau      fractionary translations
  !    irt       for each atom gives the corresponding symmetric
  !    invsym    if true the system has inversion symmetry
  !  + LDA+U-related quantities.
  !
#include "machine.h"
  use pwcom
#ifdef __PARA
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
  logical :: minus_q  
  integer, external ::n_atom_wfc, set_Hubbard_l
  !
  ! lchk_tauxk tests that atomic coordinates do not overlap
  !
  logical, external :: lchk_tauxk

  do nt = 1, ntyp
     do ir = 1, mesh (nt)
        if (r (ir, nt) > rcut) then
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
  zv (1:ntyp) = zp (1:ntyp)
  !
  !     Set the number of electrons equal to the total ionic charge
  !
  if (nelec == 0.d0) nelec = SUM(zv (ityp (1:nat) ))
  !
  ! For metals: check whether Gaussian broadening or Tetrahedron method is
  !
  lgauss = degauss .ne. 0.d0
  !
  ! Check: if there is an odd number of electrons, the crystal is a metal
  !
  if (lscf .and. abs (nint (nelec / 2.d0) - nelec / 2.d0) > 1.0d-8 &
           .and. .not.lgauss .and. .not.ltetra ) call errore ('setup', &
       'the crystal is a metal', 1)
  !
  !     Set the number of occupied bands if not given in input
  !
  if (nbnd == 0) then
     nbnd = nint (nelec) / 2.d0
     if (lgauss .or. ltetra) then
        !
        !     metallic case: add 20% more bands, with a minimum of 4
        !
        nbnd = max (nint (1.20d0 * nelec / 2.d0), nbnd+4)
     endif
  else
     if (nbnd < nint(nelec)/2.d0 .and. lscf) call errore ('setup', &
          'too few bands', 1)
  endif
  !
  ! Here we  set the precision of the diagonalization
  !
  if (ethr == 0.d0) then
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
  end if

  if (.not.lscf) niter = 1

  starting_diag_threshold = ethr
  !
  !    set number of atomic wavefunctions
  !
  natomwfc = n_atom_wfc (nat, npsx, ityp, newpseudo, nchix, nchi, oc, lchi)
  !
  !    set the max number of bands used in iterative diagonalization
  !
  nbndx = nbnd
  if (isolve == 0) then
     nbndx = david * nbnd
  endif

  if (startingwfc == 'atomic' .and. .not.restart) nbndx = max (natomwfc, nbndx)
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
  if (.not. (lchk_tauxk (nat, tau, bg) ) ) call errore ('setup', &
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
  if (ibrav == 4 .or. ibrav == 5) then
     !
     !  here the hexagonal or trigonal bravais lattice
     !
     call hexsym (at, s, sname, nrot)
     tipo = 2
  elseif (ibrav >=1  .and. ibrav <= 14) then
     !
     !  here for the cubic bravais lattice
     !
     call cubicsym (at, s, sname, nrot)
     tipo = 1
  elseif (ibrav == 0) then
     if (symm_type == 'cubic') then
        tipo = 1
        call cubicsym (at, s, sname, nrot)
     endif
     if (symm_type == 'hexagonal') then
        tipo = 2
        call hexsym (at, s, sname, nrot)
     endif
  else
     call errore ('setup', 'wrong ibrav', 1)
  endif
  !
  !   if noinv is true eliminate all symmetries which exchange z with -z
  !
  if (noinv) then
     irot = 0
     do isym = 1, nrot
        if (s (1, 3, isym) == 0 .and. s (3, 1, isym) == 0 .and. &
            s (2, 3, isym) == 0 .and. s (3, 2, isym) == 0 .and. &
            s (3, 3, isym) == 1) then
           irot = irot + 1
           s (:, :, irot) = s (:, :, isym)
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

  if (nks < 0) then
     call setupkpoint (s, nrot, xk, wk, nks, npk, nk1, &
          nk2, nk3, k1, k2, k3, at, bg, tipo)
  else if (nks == 0) then
     if (lberry) then
       call kp_strings &
       ( nppstr, gdir, nrot, s, bg, npk, k1,k2,k3, nk1,nk2,nk3, nks, xk, wk)
     else
       call kpoint_grid ( nrot, s, bg, npk, k1,k2,k3, nk1,nk2,nk3, nks, xk, wk)
     end if
  end if
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
  input_nks = nks
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

  if (iswitch > 2 .and. (nsym == 2 .and. .not.invsym .or. nsym > 2) &
     .and. .not. ( calc.eq.'mm' .or. calc.eq.'nm' ) ) &
     call errore ('setup', 'Dynamics, you should have no symmetries', -1)
  !
  !     Automatic generation of k-points
  !

  if (ltetra) then
     ntetra = 6 * nk1 * nk2 * nk3
     allocate (tetra(4,ntetra))    
     call tetrahedra (nsym,s,minus_q,at,bg,npk,k1,k2,k3, &
          nk1,nk2,nk3, nks,xk,wk,ntetra,tetra)
  else
     ntetra = 0
  end if
  !
  ! Berry phase calculation: do not change the number of k-points
  !
  if ( lberry ) nks = input_nks
  !
  ! phonon calculation: add k+q to the list of k
  !
  if (iswitch <= -2) call set_kplusq (xk, wk, xqq, nks, npk)
  !
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
  if (nks > npk) call errore ('setup', 'too many k points', nks)
#ifdef __PARA
  call init_pool
  !
  ! set the granularity for k-point distribution
  !
  if ( (abs (xqq (1) ) < eps .and. abs (xqq (2) ) < eps .and. &
        abs (xqq (3) ) < eps) .or. (iswitch > - 2) ) then
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
  !   For Bachelet-Hamann-Schluter pseudopotentials only
  !
  do nt = 1, ntyp
     if (.not.tvanp (nt) ) then
        if (nlc(nt) == 2 .and. nnl(nt) == 3 .and. bhstype (nt) ) &
             call bachel (alps (1, 0, nt), aps (1, 0, nt), 1, lmax (nt) )
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

  if (lda_plus_u) then
     Hubbard_lmax = -1
     do nt=1,ntyp
        if (Hubbard_U(nt) /= 0.d0 .or. Hubbard_alpha(nt) /= 0.d0) then
           Hubbard_l(nt) = set_Hubbard_l( psd(nt) )
           Hubbard_lmax = max(Hubbard_lmax,Hubbard_l(nt))
           write (6,*) ' HUBBARD L FOR TYPE ',psd(nt),' IS ', Hubbard_l(nt)
        end if
     end do
     write (6,*) ' MAXIMUM HUBBARD L IS ', Hubbard_lmax
     if (Hubbard_lmax == -1) &
        call errore ('setup','lda_plus_u calculation but Hubbard_l not set',1)
     call d_matrix (d1, d2, d3)  
  end if
  !
  return
end subroutine setup
!
!-----------------------------------------------------------------------
integer function n_atom_wfc &
     (nat, npsx, ityp, newpseudo, nchix, nchi, oc, lchi)
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
  !
  n_atom_wfc = 0
  do na = 1, nat
     nt = ityp (na)
     do n = 1, nchi (nt)
        if (oc (n, nt) > 0.d0 .or. .not.newpseudo (nt) ) &
             n_atom_wfc = n_atom_wfc + 2 * lchi (n, nt) + 1
     enddo
  enddo

  return
end function n_atom_wfc
