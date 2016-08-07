  !
  ! Copyright (C) 2010-2016 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino 
  !
  ! Copyright (C) 2001 PWSCF group
  ! This file is distributed under the terms of the
  ! GNU General Public License. See the file `License'
  ! in the root directory of the present distribution,
  ! or http://www.gnu.org/copyleft/gpl.txt .
  !
  ! Adapted from QE.
  !
  !-----------------------------------------------------------------------
  subroutine sgama2 (nrot, nat, s, sname, t_rev, at, bg, tau, ityp, nsym,&
     nr1, nr2, nr3, irt, ftau, npk, nks, xk, wk, invsym, minus_q, xq, &
     modenum, time_reversal, magnetic_sym, m_loc)
  !-----------------------------------------------------------------------
  !
  !     This routine performs the following tasks:
  !     1)  It finds the point group of the crystal, by eliminating the
  !         symmetries of the Bravais lattice which are not allowed
  !         by the atomic positions.
  !     1a) If xq.ne.0 it restricts the symmetries to those of the small
  !         group of q. In this case the small group of q is determined
  !         seeking all sym.op. such that Sq=q+G (and Sq=-q+G is also
  !         considered) when iswitch=-2, while only sym.op. such that
  !         Sq=q (exactly, without G) when iswitch=-3.
  !     1b) if iswitch.eq.-4 it keep only the symmetries which send
  !         a mode in itself. The mode is given by modenum
  !     2)  It finds the special points in the irreducible wedge of the
  !         true point group (or small group of q) of the crystal starting
  !         from the points in the irreducible wedge of the point group
  !         of the Bravais lattice.
  !     3)  It checks if the point group has the inversion symmetry.
  !
  !     This routine is mainly the driver of separate routines which
  !     perform each single task.
  !
  !     Modified by SdG to include the "small group of q" stuff for the
  !     linear-response preparation run.
  !
  USE kinds, only : DP
  implicit none
  !
  !    First the I/O variables
  !

  integer :: nrot, nat, s (3, 3, 48), ityp (nat), nsym, nr1, nr2, &
       nr3, irt (48, nat), ftau (3, 48), npk, nks, modenum
  ! input: number of symmetries of the original
  ! input: number of atoms in the cell
  ! input: matrices of the symmetry operations
  ! input: type of each atom
  ! output: real number of symmetry operations
  !  input: dimensions of the fft mesh
  ! output: for each atom gives the rotated ato
  ! output: fractionary translation of each sym
  ! input: maximum number of k points
  ! input-output: starting and ending number of
  ! input: the mode to be computed
  ! input: main switch of the program; used whe
  !        xq<>0 to restrict the small group of

  real(DP) :: at (3, 3), bg (3, 3), tau (3, nat), xk (3, npk), &
       wk (npk), xq (3), m_loc(3,nat)
  ! input: direct lattice vectors
  ! input: reciprocal lattice vectors
  ! input: coordinates of atomic positions
  ! input-output: coordinates of k points
  ! input-output: weights of k points
  ! input: coordinates of a q-point
  logical, intent(in) :: time_reversal, magnetic_sym
  ! time_reversal=true : use time-reversal symmetry (q=>-q)
  ! magnetic_sym =true : find symmetries that leave magnetization unchanged
  !
  logical, intent(out) :: invsym, minus_q
  ! output: if true the crystal has inversion
  ! output: if true a symmetry sends q->-q+G
  character :: sname (48) * 45
  ! input: name of the rotation part of each symmetry operation
  !
  !    And then the local variables
  !

  real(DP), allocatable :: rtau (:,:,:)
  ! direct translations of each point
  integer :: table (48, 48), irot, jrot, ipol, jpol, invs (3, 3, 48) &
       , irg (48), temp, na
  ! multiplication table of the group
  ! counter over the rotations
  ! counter over the rotations
  ! counter over the polarizations
  ! counter over the polarizations
  ! contains the inverse of each rotation
  ! gives the correspondence of symmetry
  ! operations forming a n-th coset
  ! auxilary variable
  ! counter on atoms

  integer :: t_rev(48)
  ! for magnetic symmetries: if 1 there is time reversal operation
  logical :: sym (48)
  ! if true the corresponding operation is a symmetry operation

  allocate(rtau (3, 48, nat))
  !
  !    Here we find the true symmetries of the crystal
  !
  IF ( magnetic_sym ) THEN
     CALL sgam_at_mag (nrot, s, nat, tau, ityp, bg, &
                  nr1, nr2, nr3, sym, irt, ftau, m_loc, sname, t_rev)
  ELSE
     CALL sgam_at (nrot, s, nat, tau, ityp, bg, nr1, nr2, nr3, sym, &
       irt, ftau)
  END IF
  !
  !    If xq.ne.(0,0,0) this is a preparatory run for a linear response
  !    calculation at xq. The relevant point group is therefore only the
  !    small group of q. Here we exclude from the list the symmetries
  !    that do not belong to it
  !
  call smallg_q (xq, modenum, at, bg, nrot, s, ftau, sym, minus_q)
  !
  IF ( .not. time_reversal ) THEN
     !
     minus_q=.false.
!     IF ( ABS(DOT_PRODUCT(xq,xq)) > 1.0D-07 ) CALL errore ('sgama', &
!          'phonon not implemented with non collinear magnetism', 1)
     ! If somebody wants to implement phonon calculations in non 
     ! collinear magnetic case he/she has to pay attention to the
     ! fact that in non collinear case the symmetry k -> -k is not
     ! always allowed as in collinear case. Adriano
  ENDIF
  !
  if (modenum .ne. 0) then
     call sgam_ph (at, bg, nrot, s, irt, tau, rtau, nat, sym)
     call mode_group (modenum, xq, at, bg, nat, nrot, s, irt, rtau, &
          sym, minus_q)
  endif
  !
  !    We compute the multiplication table of the group
  !
  call multable (nrot, s, table)
  !
  !   And we set the matrices of the inverse
  !
  call inverse_s (nrot, s, table, invs)
  !
  !    Find the coset in the point group of the Bravais lattice
  !
  call coset (nrot, table, sym, nsym, irg)
  !
  !    here we set the k-points in the irreducible wedge of the point grou
  !    of the crystal
  !
  call irrek (npk, nks, xk, wk, at, bg, nrot, invs, nsym, irg, minus_q)
  !
  ! copy symm. operations in sequential order so that
  ! s(i,j,irot) , irot <= nsym          are the sym.ops. of the crystal
  !               nsym+1 < irot <= nrot are the sym.ops. of the lattice
  !
  jrot = 0
  do irot = 1, nrot
     if (sym (irot) ) then
        jrot = jrot + 1
        do ipol = 1, 3
           do jpol = 1, 3
              temp = s (ipol, jpol, jrot)
              s (ipol, jpol, jrot) = s (ipol, jpol, irot)
              s (ipol, jpol, irot) = temp
           enddo
           ftau (ipol, jrot) = ftau (ipol, irot)
        enddo
        do na = 1, nat
           irt (jrot, na) = irt (irot, na)
        enddo
        sname (jrot) = sname (irot)
        t_rev(jrot) = t_rev(irot)
     endif
  enddo
  if (jrot.ne.nsym) call errore ('sgama', 'unexpected', 1)
  !
  ! Sets to zero the first matrix that is not a symmetry of the crystal.
  ! This will be used by d3toten program.
  !
  if (nrot.lt.48) then
     do ipol = 1, 3
        do jpol = 1, 3
           s (ipol, jpol, nrot + 1) = 0
        enddo
     enddo
  endif
  !
  ! check if inversion (I) is a symmetry.
  ! If so, it should be the (nsym/2+1)-th operation of the group
  !
  invsym = .true.
  irot = nsym / 2 + 1
  do ipol = 1, 3
     do jpol = 1, 3
        invsym = invsym.and.s (ipol, jpol, irot) .eq. - s (ipol, jpol, 1)
     enddo

  enddo

  deallocate (rtau)
  return

end subroutine sgama2
!-----------------------------------------------------------------------
!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
subroutine sgam_at_mag (nrot, s, nat, tau, ityp, bg, nr1, nr2, &
     nr3, sym, irt, ftau, m_loc, sname, t_rev)
  !-----------------------------------------------------------------------
  !
  !     given a point group, this routine finds the subgroup which is
  !     the point group of the crystal under consideration
  !     non symmorphic groups non allowed, provided that fractional
  !     translations are commensurate with the FFT grid
  !
  !     It sets the array sym, which for each operation of the original
  !     point group is true if this operation is also an operation of the
  !     total point group
  !
  USE io_global,   ONLY : stdout
  USE kinds
  implicit none
  !
  !     input variables
  !
  integer :: nrot, s (3, 3, 48), nat, ityp (nat), nr1, nr2, nr3
  REAL(DP) :: m_loc(3,nat), tau (3, nat), bg (3, 3)
  ! nrot : order of the parent group
  ! s    : symmetry operations of parent group
  ! nat  : number of atoms in the unit cell
  ! ityp : species of each atom in the unit cell
  ! nr*  : dimensions of the FFT mesh
  ! tau  : cartesian coordinates of the atoms
  ! at   : basis of the real-space lattice
  ! bg   :  "   "   "  reciprocal-space lattice
  !
  !     output variables
  !
  integer :: irt (48, nat), ftau (3, 48)
  logical :: sym (48)
  ! irt(isym,na) : sym.op. isym sends atom na into atom irt(isym,na)
  ! ftau(:,isym) : fractional translation associated to sym.op. isym
  !                (in FFT coordinates: crystal axis, multiplied by nr*)
  ! sym(isym)    : flag indicating if sym.op. isym in the parent group
  !                is a true symmetry operation of the crystal
  !
  !    local variables
  !
  integer :: na, kpol, nb, irot, i, j
  ! counters
  real(DP) , allocatable :: xau (:,:), rau (:,:), mxau(:,:), mrau(:,:)
  ! atomic coordinates in crystal axis
  logical :: fractional_translations
  real(DP) :: ft (3), ft1, ft2, ft3
  integer :: t_rev(48)
  character :: sname (48) * 45
  !
  external checksym
  !
  allocate(xau(3,nat))
  allocate(rau(3,nat))
  ALLOCATE(mxau(3,nat))
  ALLOCATE(mrau(3,nat))
  !
  !     Compute the coordinates of each atom in the basis of
  !     the direct lattice vectors
  !
  do na = 1, nat
     do kpol = 1, 3
        xau (kpol, na) = bg (1, kpol) * tau (1, na) + &
                         bg (2, kpol) * tau (2, na) + &
                         bg (3, kpol) * tau (3, na)
        mxau (kpol, na)= bg (1, kpol) * m_loc (1, na) + &
                         bg (2, kpol) * m_loc (2, na) + &
                         bg (3, kpol) * m_loc (3, na)
     enddo
  enddo
  !
  !      check if the identity has fractional translations
  !      (this means that the cell is actually a supercell).
  !      When this happens, fractional translations are disabled,
  !      because there is no guarantee that the generated sym.ops.
  !      form a group
  !
  nb = 1
  irot = 1
  fractional_translations = .true.
  do na = 2, nat
     if (ityp (nb) .eq.ityp (na) ) then
        ft (1) = xau(1,na) - xau(1,nb) - nint( xau(1,na) - xau(1,nb) )
        ft (2) = xau(2,na) - xau(2,nb) - nint( xau(2,na) - xau(2,nb) )
        ft (3) = xau(3,na) - xau(3,nb) - nint( xau(3,na) - xau(3,nb) )


        call checksym_mag (irot, nat, ityp, xau, xau, ft, sym, irt, mxau,&
                           mxau, t_rev(irot))

        if (sym (irot) .and. (abs (ft (1) **2 + ft (2) **2 + ft (3) ** &
             2) ) .lt.1.d-8) call errore ('sgam_at', 'overlapping atoms', na)
        if (sym (irot) ) then
           fractional_translations = .false.
           WRITE( stdout, '(5x,"Found additional translation:",3f10.4)') ft
        endif
     endif

  enddo
  do irot = 1, nrot
     !
     ! check that the grid is compatible with the S rotation
     !
     if ( mod (s (2, 1, irot) * nr1, nr2) .ne.0 .or. &
          mod (s (3, 1, irot) * nr1, nr3) .ne.0 .or. &
          mod (s (1, 2, irot) * nr2, nr1) .ne.0 .or. &
          mod (s (3, 2, irot) * nr2, nr3) .ne.0 .or. &
          mod (s (1, 3, irot) * nr3, nr1) .ne.0 .or. &
          mod (s (2, 3, irot) * nr3, nr2) .ne.0 ) then
        sym (irot) = .false.
        WRITE( stdout, '(5x,"warning: symmetry operation # ",i2, &
             &         " not compatible with FFT grid. ")') irot
        WRITE( stdout, '(3i4)') ( (s (i, j, irot) , j = 1, 3) , i = 1, 3)
        goto 100

     endif
     do na = 1, nat
        do kpol = 1, 3
           ! rau = rotated atom coordinates
           rau (kpol, na) = s (1, kpol, irot) * xau (1, na) + &
                            s (2, kpol, irot) * xau (2, na) + &
                            s (3, kpol, irot) * xau (3, na)
           mrau(kpol,na) = s(1,kpol,irot)*mxau(1,na) + & ! magnetic moment of the
                           s(2,kpol,irot)*mxau(2,na) + & ! atom rotated by the
                           s(3,kpol,irot)*mxau(3,na)   ! present symmerty
                                                       ! operation
        enddo
     enddo
     if (sname(irot)(1:3)=='inv') mrau=-mrau
     !
     !      first attempt: no fractional translation
     !
     do kpol = 1, 3
        ftau (kpol, irot) = 0
        ! input for checksym
        ft (kpol) = 0.d0
     enddo

     call checksym_mag (irot, nat, ityp, xau, rau, ft, sym, irt, mxau, &
                        mrau, t_rev(irot))
     if (.not.sym (irot) .and.fractional_translations) then
        nb = 1
        do na = 1, nat
           if (ityp (nb) .eq.ityp (na) ) then
              !
              !      second attempt: check all possible fractional translations
              !
              ft (1) = rau(1,na) - xau(1,nb) - nint( rau(1,na) - xau(1,nb) )
              ft (2) = rau(2,na) - xau(2,nb) - nint( rau(2,na) - xau(2,nb) )
              ft (3) = rau(3,na) - xau(3,nb) - nint( rau(3,na) - xau(3,nb) )

              call checksym_mag (irot,nat,ityp,xau,rau,ft,sym,irt,mxau, &
                                 mrau,t_rev(irot))
              if (sym (irot) ) then
                 ! convert ft to FFT coordinates
                 ! for later use in symmetrization
                 ft1 = ft (1) * nr1
                 ft2 = ft (2) * nr2
                 ft3 = ft (3) * nr3
                 ! check if the fractional translations are commensurate
                 ! with the FFT grid, discard sym.op. if not
                 if (abs (ft1 - nint (ft1) ) / nr1.gt.1.0d-5 .or. &
                     abs (ft2 - nint (ft2) ) / nr2.gt.1.0d-5 .or. &
                     abs (ft3 - nint (ft3) ) / nr3.gt.1.0d-5) then
                    WRITE( stdout, '(5x,"warning: symmetry operation", &
                         &     " # ",i2," not allowed.   fractional ", &
                         &     "translation:"/5x,3f11.7,"  in crystal", &
                         &     " coordinates")') irot, ft
                    sym (irot) = .false.
                 endif
                 ftau (1, irot) = nint (ft1)
                 ftau (2, irot) = nint (ft2)
                 ftau (3, irot) = nint (ft3)
                 goto 100
              endif
           endif
        enddo

     endif
100  continue
  enddo
  !
  !   deallocate work space
  !
  deallocate (rau)
  deallocate (xau)
  !
  return
END SUBROUTINE sgam_at_mag

!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------
subroutine checksym_mag (ir,nat,ityp,xau,rau,ft,sym,irt,mxau,mrau,t_rev)
  !-----------------------------------------------------------------------
  !
  !   This routine receives as input all the atomic positions xau,
  !   and the rotated rau by the symmetry operation ir. It sets to true
  !   sym(ir) if for each atom na, it is possible to find an atom nb
  !   which is of the same type of na, and coincide with it after the
  !   symmetry operation. Fractional translations are allowed.
  !
  !   Revised layout 1 may 1995 by A. Dal Corso
  !
  USE kinds
  implicit none
  !
  !     first the dummy variables
  !
  integer :: nat, ityp (nat), irt (48, nat), ir, t_rev
  ! input: the total number of atoms
  ! input: the type of each atom
  ! output: the rotated of each atom
  ! input: the rotation to be tested
  ! output: time reversal operation is present
  real(DP) :: xau (3, nat), rau (3, nat), ft (3), eps
  ! input: the initial vectors
  ! input: the rotated vectors
  ! input: the possible fractionary translat
  REAL(DP) :: mxau(3,nat), mrau(3,nat)
  ! input: the rotated vectors
  ! input: the possible fractionary translation
  logical :: sym (48)
  ! output: if true this is a symmetry opera
  !
  !  few local variables
  !
  integer :: na, nb, na1, t1, t2
  ! counter on atoms
  ! counter on atoms
  logical :: eqvect
  ! the testing function

  external eqvect
  
! SP: We have to add a tolerence factor because the interface of eqvect
  ! changed
  eps = 1d-6

  t1 = 1
  t2 = 1

  do na = 1, nat
     na1 = 0
     do nb = 1, nat
        if(ityp(na).eq.ityp(nb).and. &
           eqvect(rau (1, na),xau(1,nb),ft,eps)) na1 = nb
     enddo
     !
     IF ( na1 /= 0 ) THEN
        !
        if( abs(mrau(1,na) - mxau(1,na1))+       &
            abs(mrau(2,na) - mxau(2,na1))+       &
            abs(mrau(3,na) - mxau(3,na1)).gt.1.0D-5) t1 = 0
        if( abs(mrau(1,na) + mxau(1,na1))+       &
            abs(mrau(2,na) + mxau(2,na1))+       &
            abs(mrau(3,na) + mxau(3,na1)).gt.1.0D-5) t2 = 0
        !
     END IF
     !
     if(na1.eq.0.or.(t1+t2).eq.0) then
       sym(ir) = .false.
       t_rev = 0
       return
     else
       irt(ir,na) = na1
     endif
  enddo

  if(t1+t2.eq.2) then
    sym(ir) = .true.
    t_rev = 0
  elseif(t1.eq.1) then
    sym(ir) = .true.
    t_rev = 0
  else
    sym(ir) = .true.
    t_rev = 1
  endif

  return
END SUBROUTINE checksym_mag

!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
subroutine sgam_at (nrot, s, nat, tau, ityp, bg, nr1, nr2, &
     nr3, sym, irt, ftau)
  !-----------------------------------------------------------------------
  !
  !     given a point group, this routine finds the subgroup which is
  !     the point group of the crystal under consideration
  !     non symmorphic groups non allowed, provided that fractional
  !     translations are commensurate with the FFT grid
  !
  !     It sets the array sym, which for each operation of the original
  !     point group is true if this operation is also an operation of the
  !     total point group
  !
  USE io_global,  ONLY : stdout
  USE kinds
  implicit none
  !
  !     input variables
  !
  integer :: nrot, s (3, 3, 48), nat, ityp (nat), nr1, nr2, nr3
  real(DP) :: tau (3, nat), bg (3, 3)
  ! nrot : order of the parent group
  ! s    : symmetry operations of parent group
  ! nat  : number of atoms in the unit cell
  ! ityp : species of each atom in the unit cell
  ! nr*  : dimensions of the FFT mesh
  ! tau  : cartesian coordinates of the atoms
  ! bg   :  "   "   "  reciprocal-space lattice
  !
  !     output variables
  !
  integer :: irt (48, nat), ftau (3, 48)
  logical :: sym (48)
  ! irt(isym,na) : sym.op. isym sends atom na into atom irt(isym,na)
  ! ftau(:,isym) : fractional translation associated to sym.op. isym
  !                (in FFT coordinates: crystal axis, multiplied by nr*)
  ! sym(isym)    : flag indicating if sym.op. isym in the parent group
  !                is a true symmetry operation of the crystal
  !
  !    local variables
  !
  integer :: na, kpol, nb, irot, i, j
  ! counters
  real(DP) , allocatable :: xau (:,:), rau (:,:)
  ! atomic coordinates in crystal axis
  logical :: fractional_translations
  real(DP) :: ft (3), ft1, ft2, ft3
  !
  external checksym
  !
  allocate(xau(3,nat))
  allocate(rau(3,nat))
  !
  !     Compute the coordinates of each atom in the basis of
  !     the direct lattice vectors
  !
  do na = 1, nat
     do kpol = 1, 3
        xau (kpol, na) = bg (1, kpol) * tau (1, na) + &
                         bg (2, kpol) * tau (2, na) + &
                         bg (3, kpol) * tau (3, na)
     enddo
  enddo
  !
  !      check if the identity has fractional translations
  !      (this means that the cell is actually a supercell).
  !      When this happens, fractional translations are disabled,
  !      because there is no guarantee that the generated sym.ops.
  !      form a group
  !
  nb = 1
  irot = 1
  fractional_translations = .true.
  do na = 2, nat
     if (ityp (nb) .eq.ityp (na) ) then
        ft (1) = xau(1,na) - xau(1,nb) - nint( xau(1,na) - xau(1,nb) )
        ft (2) = xau(2,na) - xau(2,nb) - nint( xau(2,na) - xau(2,nb) )
        ft (3) = xau(3,na) - xau(3,nb) - nint( xau(3,na) - xau(3,nb) )


        call checksym (irot, nat, ityp, xau, xau, ft, sym, irt)

        if (sym (irot) .and. (abs (ft (1) **2 + ft (2) **2 + ft (3) ** &
             2) ) .lt.1.d-8) call errore ('sgam_at', 'overlapping atoms', na)
        if (sym (irot) ) then
           fractional_translations = .false.
           WRITE( stdout, '(5x,"Found additional translation:",3f10.4)') ft
        endif
     endif

  enddo
  do irot = 1, nrot
     !
     ! check that the grid is compatible with the S rotation
     !
     if ( mod (s (2, 1, irot) * nr1, nr2) .ne.0 .or. &
          mod (s (3, 1, irot) * nr1, nr3) .ne.0 .or. &
          mod (s (1, 2, irot) * nr2, nr1) .ne.0 .or. &
          mod (s (3, 2, irot) * nr2, nr3) .ne.0 .or. &
          mod (s (1, 3, irot) * nr3, nr1) .ne.0 .or. &
          mod (s (2, 3, irot) * nr3, nr2) .ne.0 ) then
        sym (irot) = .false.
        WRITE( stdout, '(5x,"warning: symmetry operation # ",i2, &
             &         " not compatible with FFT grid. ")') irot
        WRITE( stdout, '(3i4)') ( (s (i, j, irot) , j = 1, 3) , i = 1, 3)
        goto 100

     endif
     do na = 1, nat
        do kpol = 1, 3
           ! rau = rotated atom coordinates
           rau (kpol, na) = s (1, kpol, irot) * xau (1, na) + &
                            s (2, kpol, irot) * xau (2, na) + &
                            s (3, kpol, irot) * xau (3, na)
        enddo
     enddo
     !
     !      first attempt: no fractional translation
     !
     do kpol = 1, 3
        ftau (kpol, irot) = 0
        ! input for checksym
        ft (kpol) = 0.d0
     enddo
     call checksym (irot, nat, ityp, xau, rau, ft, sym, irt)
     if (.not.sym (irot) .and.fractional_translations) then
        nb = 1
        do na = 1, nat
           if (ityp (nb) .eq.ityp (na) ) then
              !
              !      second attempt: check all possible fractional translations
              !
              ft (1) = rau(1,na) - xau(1,nb) - nint( rau(1,na) - xau(1,nb) )
              ft (2) = rau(2,na) - xau(2,nb) - nint( rau(2,na) - xau(2,nb) )
              ft (3) = rau(3,na) - xau(3,nb) - nint( rau(3,na) - xau(3,nb) )

              call checksym (irot, nat, ityp, xau, rau, ft, sym, irt)
              if (sym (irot) ) then
                 ! convert ft to FFT coordinates
                 ! for later use in symmetrization
                 ft1 = ft (1) * nr1
                 ft2 = ft (2) * nr2
                 ft3 = ft (3) * nr3
                 ! check if the fractional translations are commensurate
                 ! with the FFT grid, discard sym.op. if not
                 if (abs (ft1 - nint (ft1) ) / nr1.gt.1.0d-5 .or. &
                     abs (ft2 - nint (ft2) ) / nr2.gt.1.0d-5 .or. &
                     abs (ft3 - nint (ft3) ) / nr3.gt.1.0d-5) then
                    WRITE( stdout, '(5x,"warning: symmetry operation", &
                         &     " # ",i2," not allowed.   fractional ", &
                         &     "translation:"/5x,3f11.7,"  in crystal", &
                         &     " coordinates")') irot, ft
                    sym (irot) = .false.
                 endif
                 ftau (1, irot) = nint (ft1)
                 ftau (2, irot) = nint (ft2)
                 ftau (3, irot) = nint (ft3)
                 goto 100
              endif
           endif
        enddo

     endif
100  continue
  enddo
  !
  !   deallocate work space
  !
  deallocate (rau)
  deallocate (xau)
  !
  return
end subroutine sgam_at

!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------

subroutine mode_group (modenum, xq, at, bg, nat, nrot, s, irt, &
     rtau, sym, minus_q)
  !-----------------------------------------------------------------------
  !
  ! This routine selects, among the symmetry matrices of the point group
  ! of a crystal, the symmetry operations which leave a given mode unchang
  ! For the moment it assume that the mode modenum displaces the atom
  ! modenum/3 in the direction mod(modenum,3)+1
  ! Also the minus_q operation is tested.
  !
  !  input-output variables
  !
  USE kinds
  USE constants, ONLY : tpi
  implicit none

  integer :: nat, s (3, 3, 48), irt (48, nat), nrot, modenum
  ! input: the number of atoms of the system
  ! input: the symmetry matrices
  ! input: the rotated atom
  ! input: number of symmetry operations
  ! input: the displacement pattern


  real(DP) :: xq (3), rtau (3, 48, nat), bg (3, 3), at (3, 3)
  ! input: the q point
  ! input: the translations of each atom
  ! input: the reciprocal lattice vectors
  ! input: the direct lattice vectors
  logical :: minus_q, sym (48)
  ! input: if true minus_q symmetry is used
  ! input-output: .true. if symm. op. do not change
  ! mode
  !
  !  local variables
  !

  integer :: isym, nas, ipols, na, sna, ipol, jpol
  ! counters
  ! counter on polarizations
  ! counter on polarizations

  real(DP) :: arg
  ! auxiliary

  complex(DP), allocatable :: u (:,:)
  ! the original pattern
  complex(DP)              :: fase, sum
  ! the phase of the mode
  ! check for orthogonality
  complex(DP), allocatable :: work_u (:,:), work_ru (:,:)
  ! the working pattern
  ! the rotated working pattern


  allocate(u(3, nat), work_u(3, nat), work_ru (3, nat))

  if (modenum.gt.3 * nat.or.modenum.lt.1) call errore ('mode_group', &
       'wrong modenum', 1)
  nas = (modenum - 1) / 3 + 1
  ipols = mod (modenum - 1, 3) + 1
  u (:,:) = (0.d0, 0.d0)
  u (ipols, nas) = (1.d0, 0.d0)
  do na = 1, nat
     call trnvecc (u (1, na), at, bg, - 1)
  enddo
  do isym = 1, nrot
     if (sym (isym) ) then
        do na = 1, nat
           do ipol = 1, 3
              work_u (ipol, na) = u (ipol, na)
           enddo
        enddo
        work_ru (:,:) = (0.d0, 0.d0)
        do na = 1, nat
           sna = irt (isym, na)
           arg = 0.d0
           do ipol = 1, 3
              arg = arg + xq (ipol) * rtau (ipol, isym, na)
           enddo
           arg = arg * tpi
           if (isym.eq.nrot.and.minus_q) then
              fase = CMPLX (cos (arg), sin (arg), kind=DP )
           else
              fase = CMPLX (cos (arg), - sin (arg), kind=DP )
           endif
           do ipol = 1, 3
              do jpol = 1, 3
                 work_ru (ipol, sna) = work_ru (ipol, sna) + s (jpol, ipol, &
                      isym) * work_u (jpol, na) * fase
              enddo
           enddo
        enddo
        !
        !    Transform back the rotated pattern
        !
        do na = 1, nat
           call trnvecc (work_ru (1, na), at, bg, 1)
           call trnvecc (work_u (1, na), at, bg, 1)
        enddo
        !
        !   only if the pattern remain the same ap to a phase we keep
        !   the symmetry
        !
        sum = (0.d0, 0.d0)
        do na = 1, nat
           do ipol = 1, 3
              sum = sum + CONJG(work_u (ipol, na) ) * work_ru (ipol, na)
           enddo
        enddo
        sum = abs (sum)
        if (abs (sum - 1.d0) .gt.1.d-7) sym (isym) = .false.
     endif

  enddo
  deallocate ( work_ru, work_u, u)
  return

end subroutine mode_group

!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
subroutine multable (nsym, s, table)
  !-----------------------------------------------------------------------
  !
  !  sets up the multiplication table for a group represented by 3x3
  !  integer matrices and checks that {s} is a group indeed:
  !
  !  table(n,m) = index( s(n)*s(m) )
  !
  USE kinds
  implicit none
  !
  !    here the dummy variables
  !
  integer :: nsym, s (3, 3, 48), table (48, 48)
  ! input: the number of symmetry of the
  ! input: the symmetry matrices
  ! output: the multiplication table
  !
  !  and here the local variables
  !
  integer :: irot, jrot, krot, ipol, jpol, kpol, ss (3, 3)
  ! \
  !   counter on rotations
  ! /
  ! \
  !   counters on polarizations
  ! /
  ! buffer multiplication matrix

  logical :: found, smn
  ! if true the table has been set
  ! used to check symmetries
  do irot = 1, nsym
     do jrot = 1, nsym
        !
        do ipol = 1, 3
           ! sets up th
           do jpol = 1, 3
              ! product
              ss (ipol, jpol) = 0
              ! matrix
              do kpol = 1, 3
                 !
                 ss (ipol, jpol) = ss (ipol, jpol) + s (ipol, kpol, jrot) * s ( &
                      kpol, jpol, irot)
                 ! ss=s(j)*s(
                 !
              enddo
              !
           enddo
           !
        enddo
        !
        !     here checks that the input matrices really form a group
        !     and sets the multiplication table
        !
        found = .false.
        do krot = 1, nsym
           smn = .true.
           do ipol = 1, 3
              do jpol = 1, 3
                 smn = smn.and. (s (ipol, jpol, krot) .eq.ss (ipol, jpol) )
              enddo
           enddo
           if (smn) then
              if (found) call errore ('Multable', 'Not a group', 1)
              found = .true.
              table (jrot, irot) = krot
           endif
        enddo

     enddo
     if (.not.found) call errore ('Multable', ' Not a group', 2)

  enddo
  return
end subroutine multable
!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------

subroutine inverse_s (nrot, s, table, invs)
  !---------------------------------------------------------------------
  implicit none
  integer :: nrot, s (3, 3, 48), table (48, 48), invs (3, 3, 48)
  ! input: number of symmetries of the original
  ! input: matrices of the symmetry operations
  ! input: multiplication table of the group
  ! output: contains the inverse of each rotati


  integer :: irot, jrot, ipol, jpol
  ! counter over the rotations
  ! counter over the rotations
  ! counter over the polarizations
  ! counter over the polarizations
  do irot = 1, nrot
     do jrot = 1, nrot
        if (table (irot, jrot) .eq.1) then
           do ipol = 1, 3
              do jpol = 1, 3
                 invs (ipol, jpol, irot) = s (ipol, jpol, jrot)
              enddo
           enddo
        endif
     enddo

  enddo
  return
end subroutine inverse_s

!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
subroutine coset (nrot, table, sym, nsym, irg)
  !-----------------------------------------------------------------------
  !
  !  Divides the elements of a given group into left cosets of one
  !  of its subgroups.
  !  The input is the array sym which is true only for the
  !  operations of the subgroup, the output is nsym, and the array irg,
  !  which contains as its first elements the indices of the subgroup,
  !  and then its right cosets.
  !
  !  revised layout 1 may 1995 by A. Dal Corso
  !
  USE kinds
  implicit none
  !
  !    first the dummy variables
  !
  integer :: nrot, table (48, 48), nsym, irg (48)
  ! input: order of the group
  ! input: multiplication table of the group
  ! output: order of the subgroup
  ! output: gives the correspondence of symme
  ! operations forming a n-th coset
  ! input: flag indicating if an operations
  logical :: sym (48)
  ! belongs to the subgroup
  !
  ! here the local variables
  !
  logical :: done (48)
  ! if true the operation has been already ch

  integer :: irot, ncos, isym, nc, nelm
  ! counter on rotations
  ! number of cosets (=nrot/nsym)
  ! counter on symmetries
  ! counter on cosets
  ! counter on the number of elements
  !
  !    here we count the elements of the subgroup and set the first part o
  !    irg which contain the subgroup
  !
  nsym = 0
  do irot = 1, nrot
     done (irot) = sym (irot)
     if (sym (irot) ) then
        nsym = nsym + 1
        irg (nsym) = irot
     endif
  enddo
  !
  !     we check that the order of the subgroup is a divisor of the order
  !     total group. ncos is the number of cosets
  !
  IF ( nsym == 0 ) CALL errore( 'coset', 'nsym == 0', 1 ) 
  !
  ncos = nrot / nsym
  if (ncos * nsym.ne.nrot) call errore ('coset', &
  'The order'//' of the group is not a multiple of that of the subgroup', 1)
  !
  !     here we set the other elements of irg, by using the multiplication
  !
  nelm = nsym
  do nc = 2, ncos
     do irot = 1, nrot
        if (.not.done (irot) ) then
           do isym = 1, nsym
              nelm = nelm + 1
              irg (nelm) = table (irot, irg (isym) )
              done (irg (nelm) ) = .true.
           enddo
        endif
     enddo

  enddo
  return
end subroutine coset
!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
subroutine irrek (npk, nks, xk, wk, at, bg, nrot, invs, nsym, irg, &
     minus_q)
  !-----------------------------------------------------------------------
  !
  !  Given a set of special points in the Irreducible Wedge of some
  !  group, finds the equivalent special points in the IW of one of
  !  its subgroups.
  !
  USE kinds, only : DP
  implicit none
  !
  integer, intent(inout) :: nks
  ! number of special points
  integer, intent(in) :: npk, nrot, nsym, invs (3, 3, 48), irg (nrot)
  ! maximum number of special points
  ! order of the parent point group
  ! order of the subgroup
  ! inverse of the elements of the symmetry group
  ! partition of the elements of the symmetry group into left cosets,
  ! as given by SUBROUTINE COSET
  real(DP), intent(inout) :: xk (3, npk), wk (npk)
  ! special points and weights
  real(DP), intent(in) :: at (3, 3), bg (3, 3)
  ! basis vectors of the Bravais and reciprocal lattice
  logical, intent(in) :: minus_q
  ! .true. if symmetries q = -q+G are acceptable
  !
  !    here the local variables
  !
  integer :: nks0, jk, kpol, irot, jrot, ncos, jc, ic, isym
  ! nks0: used to save the initial number of k-points
  ! ncos: total number of cosets
  real(DP) :: xkg (3), xks (3, 48), w (48), sw, one
  ! coordinates of the k point in crystal axis
  ! coordinates of the rotated k point
  ! weight of each coset
  ! buffer which contains the weight of k points
  ! total weight of k-points
  logical :: latm, satm
  ! true if a k-point is equivalent to a previous one
  ! true if equivalent point found
  real(DP) :: eps=1.0d-8

  nks0 = nks
  do jk = 1, nks0
     !
     !     The k point is first computed in crystal axis
     !
     do kpol = 1, 3
        ! xkg are the components ofx k in the crystal RL base
        xkg (kpol) = at (1, kpol) * xk (1, jk) + &
                     at (2, kpol) * xk (2, jk) + &
                     at (3, kpol) * xk (3, jk)
     enddo
     !
     !   Then it is rotated with each symmetry of the global group. Note that
     !   the irg vector is used to divide all the rotated vector in cosets
     !
     do irot = 1, nrot
        jrot = irg (irot)
        do kpol = 1, 3
           ! the rotated of xkg with respect to the group operations
           xks (kpol, irot) = invs (kpol, 1, jrot) * xkg (1) + &
                              invs (kpol, 2, jrot) * xkg (2) + &
                              invs (kpol, 3, jrot) * xkg (3)
        enddo
     enddo
     !
     !    For each coset one point is tested with all the preceding
     !
     ncos = nrot / nsym
     do ic = 1, ncos
        irot = (ic - 1) * nsym + 1
        latm = .false.
        !
        !  latm = .true. if the present k-vector is equivalent to some previous
        !
        do jc = 1, ic - 1
           do isym = 1, nsym
              !
              !   satm = .true. if the present symmetry operation makes 
              !   the ir and ik k-vectors equivalent ...
              !
              jrot = (jc - 1) * nsym + isym
              satm = abs (xks (1, irot) - xks (1, jrot) - &
                     nint (xks (1, irot) - xks (1, jrot) ) ) < 1.0d-5 .and. &
                     abs (xks (2, irot) - xks (2, jrot) - &
                     nint (xks (2, irot) - xks (2, jrot) ) ) < 1.0d-5 .and. &
                     abs (xks (3, irot) - xks (3, jrot) - &
                     nint (xks (3, irot) - xks (3, jrot) ) ) < 1.0d-5
              !
              !  .... or equivalent to minus each other when minus_q=.t.
              !
              if (minus_q) satm = satm .or. &
                   abs (xks (1, irot) + xks (1, jrot) - &
                   nint (xks (1, irot) + xks (1, jrot) ) ) < 1.0d-5 .and. &
                   abs (xks (2, irot) + xks (2, jrot) - &
                   nint (xks (2, irot) + xks (2, jrot) ) ) < 1.0d-5 .and. &
                   abs (xks (3, irot) + xks (3, jrot) - &
                   nint (xks (3, irot) + xks (3, jrot) ) ) < 1.0d-5
              latm = latm .or. satm
              if (satm .and. ABS(w (jc)) < eps) then
                 w (jc) = w (jc) + 1.d0
                 goto 100
              endif
           enddo

        enddo
100     continue
        if (latm) then
           w (ic) = 0.d0
        else
           w (ic) = 1.d0
        endif
     enddo
     !
     !     here the k-point list is updated
     !
     sw = wk (jk) / SUM (w(1:ncos))
     wk (jk) = sw * w (1)
     do ic = 2, ncos
        irot = (ic - 1) * nsym + 1
        if (ABS(w(ic)) < eps) then
           nks = nks + 1
           if (nks > npk) call errore ('irrek', 'too many k-points', nks)
           wk (nks) = sw * w (ic)
           do kpol = 1, 3
              xk (kpol, nks) = bg (kpol, 1) * xks (1, irot) + &
                               bg (kpol, 2) * xks (2, irot) + &
                               bg (kpol, 3) * xks (3, irot)
           enddo
        endif
     enddo

  enddo
  !
  ! normalize weights to one
  !
  one = SUM (wk(1:nks))
  if ( one > 0.d0 ) wk(1:nks) = wk(1:nks) / one
  !
  return
end subroutine irrek

!-----------------------------------------------------------------------
subroutine checksym (ir, nat, ityp, xau, rau, ft, sym, irt)
  !-----------------------------------------------------------------------
  !
  !   This routine receives as input all the atomic positions xau,
  !   and the rotated rau by the symmetry operation ir. It sets to true
  !   sym(ir) if for each atom na, it is possible to find an atom nb
  !   which is of the same type of na, and coincide with it after the
  !   symmetry operation. Fractional translations are allowed.
  !
  !   Revised layout 1 may 1995 by A. Dal Corso
  !
  USE kinds
!  USE symm_base,        ONLY : accep
  implicit none
  !
  !     first the dummy variables
  !
  integer :: nat, ityp (nat), irt (48, nat), ir
  ! input: the total number of atoms
  ! input: the type of each atom
  ! output: the rotated of each atom
  ! input: the rotation to be tested
  real(DP) :: xau (3, nat), rau (3, nat), ft (3)
  ! input: the initial vectors
  ! input: the rotated vectors
  ! input: the possible fractionary translat
  logical :: sym (48)
  ! output: if true this is a symmetry opera
  !
  !  few local variables
  !
  integer :: na, nb
  ! counter on atoms
  ! counter on atoms
  logical :: eqvect
  ! the testing function

  external eqvect
  do na = 1, nat
     do nb = 1, nat
      !  sym (ir) = ityp (na) .eq.ityp (nb) .and.eqvect (rau (1, na), &
      !       xau (1, nb), ft,accep)
        sym (ir) = ityp (na) .eq.ityp (nb) .and.eqvect (rau (1, na), &
             xau (1, nb), ft,1.0d-5)
        if (sym (ir) ) then
           !
           ! the rotated atom does coincide with one of the like atoms
           ! keep track of which atom the rotated atom coincides with
           !
           irt (ir, na) = nb
           goto 10
        endif
     enddo
     !
     ! the rotated atom does not coincide with any of the like atoms
     ! s(ir) + ft is not a symmetry operation
     !
     return
10   continue
  enddo
  !
  ! s(ir) + ft is a symmetry operation
  !
  return
end subroutine checksym
