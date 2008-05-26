!
! Copyright (C) 2001-2008 Quantum-Espresso group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
subroutine sgam_at (nrot, s, nat, tau, ityp, at, bg, nr1, nr2, &
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
  integer, intent(in) :: nrot, s (3, 3, 48), nat, ityp (nat), nr1, nr2, nr3
  ! nrot : order of the parent group
  ! s    : symmetry operations of parent group
  ! nat  : number of atoms in the unit cell
  ! ityp : species of each atom in the unit cell
  ! nr*  : dimensions of the FFT mesh
  !
  real(DP), intent(in) :: tau (3, nat), at (3, 3), bg (3, 3)
  !
  ! tau  : cartesian coordinates of the atoms
  ! at   : basis of the real-space lattice
  ! bg   :  "   "   "  reciprocal-space lattice
  !
  !     output variables
  !
  integer, intent(out) :: irt (48, nat), ftau (3, 48)
  logical, intent(out) :: sym (48)
  ! irt(isym,na) : sym.op. isym sends atom na into atom irt(isym,na)
  ! ftau(:,isym) : fractional translation associated to sym.op. isym
  !                (in FFT coordinates: crystal axis, multiplied by nr*)
  ! sym(isym)    : flag indicating if sym.op. isym in the parent group
  !                is a true symmetry operation of the crystal
  !
  integer :: na, kpol, nb, irot, i, j
  ! counters
  real(DP) , allocatable :: xau (:,:), rau (:,:)
  ! atomic coordinates in crystal axis
  logical :: fractional_translations
  real(DP) :: ft (3), ft1, ft2, ft3
  !
  allocate(xau(3,nat))
  allocate(rau(3,nat))
  !
  !     Compute the coordinates of each atom in the basis of
  !     the direct lattice vectors
  !
  do na = 1, nat
     xau(:,na) = bg(1,:) * tau(1,na) + bg(2,:) * tau(2,na) + bg(3,:) * tau(3,na)
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
     if (ityp (nb) == ityp (na) ) then
        ft (:) = xau(:,na) - xau(:,nb) - nint( xau(:,na) - xau(:,nb) )
        !
        call checksym (irot, nat, ityp, xau, xau, ft, sym, irt)
        !
        if ( sym (irot) .and. &
            (abs (ft(1) **2 + ft(2) **2 + ft (3) **2) < 1.d-8) ) &
            call errore ('sgam_at', 'overlapping atoms', na)
        if (sym (irot) ) then
           fractional_translations = .false.
           WRITE( stdout, '(5x,"Fractionary translation:",3f10.4, &
          &   "is a symmetry operation:" / 5x, "This is a supercell, ",&
          &   "fractionary translation are disabled:",3f10.4)') ft
        endif
     endif
  enddo
  !
  do irot = 1, nrot
     !
     ! check that the grid is compatible with the S rotation
     !
     if ( mod (s (2, 1, irot) * nr1, nr2) /= 0 .or. &
          mod (s (3, 1, irot) * nr1, nr3) /= 0 .or. &
          mod (s (1, 2, irot) * nr2, nr1) /= 0 .or. &
          mod (s (3, 2, irot) * nr2, nr3) /= 0 .or. &
          mod (s (1, 3, irot) * nr3, nr1) /= 0 .or. &
          mod (s (2, 3, irot) * nr3, nr2) /= 0 ) then
        sym (irot) = .false.
        WRITE( stdout, '(5x,"warning: symmetry operation # ",i2, &
             &         " not compatible with FFT grid. ")') irot
        WRITE( stdout, '(3i4)') ( (s (i, j, irot) , j = 1, 3) , i = 1, 3)
        goto 100

     endif
     do na = 1, nat
        ! rau = rotated atom coordinates
        rau (:, na) = s (1,:, irot) * xau (1, na) + &
                      s (2,:, irot) * xau (2, na) + &
                      s (3,:, irot) * xau (3, na)
     enddo
     !
     !      first attempt: no fractional translation
     !
     ftau (:, irot) = 0
     ft (:) = 0.d0
     !
     call checksym (irot, nat, ityp, xau, rau, ft, sym, irt)
     !
     if (.not.sym (irot) .and. fractional_translations) then
        nb = 1
        do na = 1, nat
           if (ityp (nb) == ityp (na) ) then
              !
              !      second attempt: check all possible fractional translations
              !
              ft (:) = rau(:,na) - xau(:,nb) - nint( rau(:,na) - xau(:,nb) )
              !
              call checksym (irot, nat, ityp, xau, rau, ft, sym, irt)
              !
              if (sym (irot) ) then
                 ! convert ft to FFT coordinates
                 ! for later use in symmetrization
                 ft1 = ft (1) * nr1
                 ft2 = ft (2) * nr2
                 ft3 = ft (3) * nr3
                 ! check if the fractional translations are commensurate
                 ! with the FFT grid, discard sym.op. if not
                 if (abs (ft1 - nint (ft1) ) / nr1 > 1.0d-5 .or. &
                     abs (ft2 - nint (ft2) ) / nr2 > 1.0d-5 .or. &
                     abs (ft3 - nint (ft3) ) / nr3 > 1.0d-5) then
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

