!
! Copyright (C) 2010-2015 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------------------------------------------
!
MODULE symm_base

  USE kinds,      ONLY : DP
  USE io_global,  ONLY : stdout
  USE cell_base,  ONLY : at, bg
  !
  ! ... The variables needed to describe the symmetry properties
  ! ... and the routines to find crystal symmetries
  !
  ! ... these are acceptance criteria
  !
  REAL(DP), PARAMETER :: eps1 = 1.0d-6, eps2 = 1.0d-5
  !
  SAVE
  !
  PRIVATE
  !
  ! ... Exported variables
  !
  PUBLIC :: s, sr, sname, ft, ftau, nrot, nsym, nsym_ns, nsym_na, t_rev, &
            no_t_rev, time_reversal, irt, invs, invsym, d1, d2, d3, &
            allfrac, nofrac, nosym, nosym_evc
  INTEGER :: &
       s(3,3,48),            &! symmetry matrices, in crystal axis
       invs(48),             &! index of inverse operation: S^{-1}_i=S(invs(i))
       ftau(3,48),           &! fractional translations, in FFT coordinates
       nrot,                 &! number of bravais lattice symmetries
       nsym = 1,             &! total number of crystal symmetries
       nsym_ns = 0,          &! nonsymmorphic (fractional translation) symms
       nsym_na = 0            ! excluded nonsymmorphic symmetries because
                              ! fract. transl. is noncommensurate with FFT grid
  REAL (DP) :: &
       ft (3,48),            &! fractional translations, in crystal axis
       sr (3,3,48),          &! symmetry matrices, in cartesian axis
       accep = 1.0d-5         ! initial value of the acceptance threshold
                              ! for position comparison by eqvect in checksym
  !
  ! ... note: ftau are used for symmetrization in real space (phonon, exx)
  ! ... in which case they must be commensurated with the FFT grid
  !
  CHARACTER(len=45) ::  sname(48)   ! name of the symmetries
  INTEGER :: &
       t_rev(48) = 0          ! time reversal flag, for noncolinear magnetism
  INTEGER, ALLOCATABLE :: &
       irt(:,:)               ! symmetric atom for each atom and sym.op.
  LOGICAL :: &
       time_reversal=.true., &! if .TRUE. the system has time reversal symmetry
       invsym,               &! if .TRUE. the system has inversion symmetry
       nofrac= .false.,      &! if .TRUE. fract. translations are not allowed
       allfrac= .false.,     &! if .TRUE. all fractionary translations allowed,
                              ! even those not commensurate with FFT grid
       nosym = .false.,      &! if .TRUE. no symmetry is used
       nosym_evc = .false.,  &! if .TRUE. symmetry is used only to symmetrize
                              ! k points
       no_t_rev=.false.       ! if .TRUE. remove the symmetries that
                              ! require time reversal
  REAL(DP),TARGET :: &
       d1(3,3,48),           &! matrices for rotating spherical
       d2(5,5,48),           &! harmonics (d1 for l=1, ...)
       d3(7,7,48)             !
  !
  ! ... Exported routines
  !
  PUBLIC ::  find_sym, inverse_s, copy_sym, checkallsym, &
             s_axis_to_cart, set_sym, set_sym_bl, find_sym_ifc
  !
CONTAINS
   !
   SUBROUTINE inverse_s ( )
     !-----------------------------------------------------------------------
     !
     ! Locate index of S^{-1}
     !
     IMPLICIT NONE
     !
     INTEGER :: isym, jsym, ss (3, 3)
     LOGICAL :: found
     !
     DO isym = 1, nsym
        found = .false.
        DO jsym = 1, nsym
           !
           ss = matmul (s(:,:,jsym),s(:,:,isym))
           ! s(:,:,1) is the identity
           IF ( all ( s(:,:,1) == ss(:,:) ) ) THEN
              invs (isym) = jsym
              found = .true.
           ENDIF
        ENDDO
        IF ( .not.found) CALL errore ('inverse_s', ' Not a group', 1)
     ENDDO
     !
   END SUBROUTINE inverse_s
   !
!-----------------------------------------------------------------------
SUBROUTINE set_sym_bl ( )
  !-----------------------------------------------------------------------
  !
  ! Provides symmetry operations for all bravais lattices
  ! Tests first the 24 proper rotations for the cubic lattice;
  ! then the 8 rotations specific for the hexagonal axis (special axis c);
  ! then inversion is added
  !
  USE matrix_inversion
  IMPLICIT NONE
  !
  ! sin3 = sin(pi/3), cos3 = cos(pi/3), msin3 = -sin(pi/3), mcos3 = -cos(pi/3)
  !
  REAL(DP), PARAMETER :: sin3 = 0.866025403784438597d0, cos3 = 0.5d0, &
                        msin3 =-0.866025403784438597d0, mcos3 = -0.5d0
  REAL(DP) :: s0(3, 3, 32), overlap (3, 3), rat (3), rot (3, 3), value
  ! s0: the s matrices in cartesian axis
  ! overlap: inverse overlap matrix between direct lattice
  ! rat: the rotated of a direct vector ( cartesian )
  ! rot: the rotated of a direct vector ( crystal axis )
  ! value: component of the s matrix in axis basis
  INTEGER :: jpol, kpol, mpol, irot, imat(32)
  ! counters over the polarizations and the rotations

  CHARACTER (len=45) :: s0name (64)
  ! full name of the rotational part of each symmetry operation

  data s0/ 1.d0,  0.d0,  0.d0,  0.d0,  1.d0,  0.d0,  0.d0,  0.d0,  1.d0, &
          -1.d0,  0.d0,  0.d0,  0.d0, -1.d0,  0.d0,  0.d0,  0.d0,  1.d0, &
          -1.d0,  0.d0,  0.d0,  0.d0,  1.d0,  0.d0,  0.d0,  0.d0, -1.d0, &
           1.d0,  0.d0,  0.d0,  0.d0, -1.d0,  0.d0,  0.d0,  0.d0, -1.d0, &
           0.d0,  1.d0,  0.d0,  1.d0,  0.d0,  0.d0,  0.d0,  0.d0, -1.d0, &
           0.d0, -1.d0,  0.d0, -1.d0,  0.d0,  0.d0,  0.d0,  0.d0, -1.d0, &
           0.d0, -1.d0,  0.d0,  1.d0,  0.d0,  0.d0,  0.d0,  0.d0,  1.d0, &
           0.d0,  1.d0,  0.d0, -1.d0,  0.d0,  0.d0,  0.d0,  0.d0,  1.d0, &
           0.d0,  0.d0,  1.d0,  0.d0, -1.d0,  0.d0,  1.d0,  0.d0,  0.d0, &
           0.d0,  0.d0, -1.d0,  0.d0, -1.d0,  0.d0, -1.d0,  0.d0,  0.d0, &
           0.d0,  0.d0, -1.d0,  0.d0,  1.d0,  0.d0,  1.d0,  0.d0,  0.d0, &
           0.d0,  0.d0,  1.d0,  0.d0,  1.d0,  0.d0, -1.d0,  0.d0,  0.d0, &
          -1.d0,  0.d0,  0.d0,  0.d0,  0.d0,  1.d0,  0.d0,  1.d0,  0.d0, &
          -1.d0,  0.d0,  0.d0,  0.d0,  0.d0, -1.d0,  0.d0, -1.d0,  0.d0, &
           1.d0,  0.d0,  0.d0,  0.d0,  0.d0, -1.d0,  0.d0,  1.d0,  0.d0, &
           1.d0,  0.d0,  0.d0,  0.d0,  0.d0,  1.d0,  0.d0, -1.d0,  0.d0, &
           0.d0,  0.d0,  1.d0,  1.d0,  0.d0,  0.d0,  0.d0,  1.d0,  0.d0, &
           0.d0,  0.d0, -1.d0, -1.d0,  0.d0,  0.d0,  0.d0,  1.d0,  0.d0, &
           0.d0,  0.d0, -1.d0,  1.d0,  0.d0,  0.d0,  0.d0, -1.d0,  0.d0, &
           0.d0,  0.d0,  1.d0, -1.d0,  0.d0,  0.d0,  0.d0, -1.d0,  0.d0, &
           0.d0,  1.d0,  0.d0,  0.d0,  0.d0,  1.d0,  1.d0,  0.d0,  0.d0, &
           0.d0, -1.d0,  0.d0,  0.d0,  0.d0, -1.d0,  1.d0,  0.d0,  0.d0, &
           0.d0, -1.d0,  0.d0,  0.d0,  0.d0,  1.d0, -1.d0,  0.d0,  0.d0, &
           0.d0,  1.d0,  0.d0,  0.d0,  0.d0, -1.d0, -1.d0,  0.d0,  0.d0, &
           cos3,  sin3, 0.d0, msin3,  cos3, 0.d0, 0.d0, 0.d0,  1.d0, &
           cos3, msin3, 0.d0,  sin3,  cos3, 0.d0, 0.d0, 0.d0,  1.d0, &
          mcos3,  sin3, 0.d0, msin3, mcos3, 0.d0, 0.d0, 0.d0,  1.d0, &
          mcos3, msin3, 0.d0,  sin3, mcos3, 0.d0, 0.d0, 0.d0,  1.d0, &
           cos3, msin3, 0.d0, msin3, mcos3, 0.d0, 0.d0, 0.d0, -1.d0, &
           cos3,  sin3, 0.d0,  sin3, mcos3, 0.d0, 0.d0, 0.d0, -1.d0, &
          mcos3, msin3, 0.d0, msin3,  cos3, 0.d0, 0.d0, 0.d0, -1.d0, &
          mcos3,  sin3, 0.d0,  sin3,  cos3, 0.d0, 0.d0, 0.d0, -1.d0 /

  data s0name/  'identity                                     ',&
                '180 deg rotation - cart. axis [0,0,1]        ',&
                '180 deg rotation - cart. axis [0,1,0]        ',&
                '180 deg rotation - cart. axis [1,0,0]        ',&
                '180 deg rotation - cart. axis [1,1,0]        ',&
                '180 deg rotation - cart. axis [1,-1,0]       ',&
                ' 90 deg rotation - cart. axis [0,0,-1]       ',&
                ' 90 deg rotation - cart. axis [0,0,1]        ',&
                '180 deg rotation - cart. axis [1,0,1]        ',&
                '180 deg rotation - cart. axis [-1,0,1]       ',&
                ' 90 deg rotation - cart. axis [0,1,0]        ',&
                ' 90 deg rotation - cart. axis [0,-1,0]       ',&
                '180 deg rotation - cart. axis [0,1,1]        ',&
                '180 deg rotation - cart. axis [0,1,-1]       ',&
                ' 90 deg rotation - cart. axis [-1,0,0]       ',&
                ' 90 deg rotation - cart. axis [1,0,0]        ',&
                '120 deg rotation - cart. axis [-1,-1,-1]     ',&
                '120 deg rotation - cart. axis [-1,1,1]       ',&
                '120 deg rotation - cart. axis [1,1,-1]       ',&
                '120 deg rotation - cart. axis [1,-1,1]       ',&
                '120 deg rotation - cart. axis [1,1,1]        ',&
                '120 deg rotation - cart. axis [-1,1,-1]      ',&
                '120 deg rotation - cart. axis [1,-1,-1]      ',&
                '120 deg rotation - cart. axis [-1,-1,1]      ',&
                ' 60 deg rotation - cryst. axis [0,0,1]       ',&
                ' 60 deg rotation - cryst. axis [0,0,-1]      ',&
                '120 deg rotation - cryst. axis [0,0,1]       ',&
                '120 deg rotation - cryst. axis [0,0,-1]      ',&
                '180 deg rotation - cryst. axis [1,-1,0]      ',&
                '180 deg rotation - cryst. axis [2,1,0]       ',&
                '180 deg rotation - cryst. axis [0,1,0]       ',&
                '180 deg rotation - cryst. axis [1,1,0]       ',&
                'inversion                                    ',&
                'inv. 180 deg rotation - cart. axis [0,0,1]   ',&
                'inv. 180 deg rotation - cart. axis [0,1,0]   ',&
                'inv. 180 deg rotation - cart. axis [1,0,0]   ',&
                'inv. 180 deg rotation - cart. axis [1,1,0]   ',&
                'inv. 180 deg rotation - cart. axis [1,-1,0]  ',&
                'inv.  90 deg rotation - cart. axis [0,0,-1]  ',&
                'inv.  90 deg rotation - cart. axis [0,0,1]   ',&
                'inv. 180 deg rotation - cart. axis [1,0,1]   ',&
                'inv. 180 deg rotation - cart. axis [-1,0,1]  ',&
                'inv.  90 deg rotation - cart. axis [0,1,0]   ',&
                'inv.  90 deg rotation - cart. axis [0,-1,0]  ',&
                'inv. 180 deg rotation - cart. axis [0,1,1]   ',&
                'inv. 180 deg rotation - cart. axis [0,1,-1]  ',&
                'inv.  90 deg rotation - cart. axis [-1,0,0]  ',&
                'inv.  90 deg rotation - cart. axis [1,0,0]   ',&
                'inv. 120 deg rotation - cart. axis [-1,-1,-1]',&
                'inv. 120 deg rotation - cart. axis [-1,1,1]  ',&
                'inv. 120 deg rotation - cart. axis [1,1,-1]  ',&
                'inv. 120 deg rotation - cart. axis [1,-1,1]  ',&
                'inv. 120 deg rotation - cart. axis [1,1,1]   ',&
                'inv. 120 deg rotation - cart. axis [-1,1,-1] ',&
                'inv. 120 deg rotation - cart. axis [1,-1,-1] ',&
                'inv. 120 deg rotation - cart. axis [-1,-1,1] ',&
                'inv.  60 deg rotation - cryst. axis [0,0,1]  ',&
                'inv.  60 deg rotation - cryst. axis [0,0,-1] ',&
                'inv. 120 deg rotation - cryst. axis [0,0,1]  ',&
                'inv. 120 deg rotation - cryst. axis [0,0,-1] ',&
                'inv. 180 deg rotation - cryst. axis [1,-1,0] ',&
                'inv. 180 deg rotation - cryst. axis [2,1,0]  ',&
                'inv. 180 deg rotation - cryst. axis [0,1,0]  ',&
                'inv. 180 deg rotation - cryst. axis [1,1,0]  ' /

  !    compute the overlap matrix for crystal axis

  DO jpol = 1,3
     DO kpol = 1,3
        rot(kpol,jpol) = at(1,kpol)*at(1,jpol) +&
                         at(2,kpol)*at(2,jpol) +&
                         at(3,kpol)*at(3,jpol)
     ENDDO
  ENDDO
  !
  !    then its inverse (rot is used as work space)
  !
  CALL invmat (3, rot, overlap)

  nrot = 1
  DO irot = 1,32
     !
     !   for each possible symmetry
     !
     DO jpol = 1,3
        DO mpol = 1,3
           !
           !   compute, in cartesian coordinates the rotated vector
           !
           rat(mpol) = s0(mpol,1,irot)*at(1,jpol) +&
                       s0(mpol,2,irot)*at(2,jpol) +&
                       s0(mpol,3,irot)*at(3,jpol)
        ENDDO

        DO kpol = 1,3
           !
           !   the rotated vector is projected on the direct lattice
           !
           rot(kpol,jpol) = at(1,kpol)*rat(1) +&
                            at(2,kpol)*rat(2) +&
                            at(3,kpol)*rat(3)
        ENDDO
     ENDDO
     !
     !  and the inverse of the overlap matrix is applied
     !
     DO jpol = 1,3
        DO kpol = 1,3
           value = overlap(jpol,1)*rot(1,kpol) +&
           &       overlap(jpol,2)*rot(2,kpol) +&
           &       overlap(jpol,3)*rot(3,kpol)
           IF ( abs(dble(nint(value))-value) > eps1 ) THEN
              !
              ! if a noninteger is obtained, this implies that this operation
              ! is not a symmetry operation for the given lattice
              !
              GOTO 10
           ENDIF
           s(kpol,jpol,nrot) = nint(value)
        ENDDO
     ENDDO
     sname(nrot)=s0name(irot)
     imat(nrot)=irot
     nrot = nrot+1
10   CONTINUE
  ENDDO
  nrot = nrot-1
  IF ( nrot /= 1 .AND. nrot /= 2 .AND. nrot /= 4 .AND. nrot /= 6 .AND. &
       nrot /= 8 .AND. nrot /=12 .AND. nrot /=24 ) THEN
       WRITE (stdout, '(80("-"),/,"NOTICE: Bravais lattice has wrong number (",&
      & i2,") of symmetries - symmetries are disabled",/,80("-"))' ) nrot
      nrot = 1
  END IF
  !
  !     set the inversion symmetry ( Bravais lattices have always inversion
  !     symmetry )
  !
  DO irot = 1, nrot
     sname(irot+nrot) = s0name(imat(irot)+32)
     DO kpol = 1,3
        DO jpol = 1,3
           s(kpol,jpol,irot+nrot) = -s(kpol,jpol,irot)
        ENDDO
     ENDDO
  ENDDO
  nrot = 2*nrot
  !
  !    reset fractional translations to zero before checking the group
  ! 
  ft(:,:) = 0.0_dp
  IF ( .not. is_group ( nrot ) ) THEN
  !    This happens for instance for an hexagonal lattice with one axis 
  !    oriented at 15 degrees from the x axis, the other along (-1,1,0)
      WRITE (stdout, '(80("-"),/,"NOTICE: Symmetry group for Bravais lattice &
     & is not a group - symmetries are disabled",/,80("-"))' ) nrot
      nrot = 1
  ENDIF
  !
  RETURN
  !
END SUBROUTINE set_sym_bl
!
!-----------------------------------------------------------------------
SUBROUTINE find_sym ( nat, tau, ityp, nr1, nr2, nr3, magnetic_sym, m_loc, &
                      no_z_inv )
  !-----------------------------------------------------------------------
  !
  !     This routine finds the point group of the crystal, by eliminating
  !     the symmetries of the Bravais lattice which are not allowed
  !     by the atomic positions (or by the magnetization if present)
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(in) :: nat, ityp (nat), nr1, nr2, nr3
  REAL(DP), INTENT(in) :: tau (3,nat), m_loc(3,nat)
  LOGICAL, INTENT(in) :: magnetic_sym
  LOGICAL, INTENT(IN), OPTIONAL :: no_z_inv
  ! no_z_inv: if .true., disable symmetries sending z into -z.
  !           Some calculations (e.g. monopole fields) require this
  !
  INTEGER :: i
  LOGICAL :: sym (48)
  ! if true the corresponding operation is a symmetry operation
  !
  IF ( .not. allocated(irt) ) ALLOCATE( irt( 48, nat ) )
  irt( :, : ) = 0
  !
  !    Here we find the true symmetries of the crystal
  !
  symm: DO i=1,3 !emine: if it is not resolved in 3 steps it is sth else?
    IF ( PRESENT(no_z_inv) ) THEN
       CALL sgam_at ( nat, tau, ityp, nr1, nr2, nr3, sym, no_z_inv )
    ELSE
       CALL sgam_at ( nat, tau, ityp, nr1, nr2, nr3, sym )
    ENDIF
    !
    !    Here we check for magnetic symmetries
    !
    IF ( magnetic_sym ) CALL sgam_at_mag ( nat, m_loc, sym )
    !
    !  If nosym_evc is true from now on we do not use the symmetry any more
    !
    IF (nosym_evc) THEN
       sym=.false.
       sym(1)=.true.
    ENDIF
    !
    !    Here we re-order all rotations in such a way that true sym.ops
    !    are the first nsym; rotations that are not sym.ops. follow
    !
    nsym = copy_sym ( nrot, sym )
    !
    IF ( .not. is_group ( nsym ) ) THEN
       IF (i == 1) CALL infomsg ('find_sym', &
                      'Not a group! Trying with lower acceptance parameter...')
       accep = accep * 0.5d0
       IF (i == 3) THEN
         CALL infomsg ('find_sym', 'Still not a group! symmetry disabled')
         nsym = 1
       ENDIF
       CYCLE symm
    ELSE
       IF (i > 1) CALL infomsg ('find_sym', 'Symmetry operations form a group')
       exit symm
    ENDIF
  ENDDO symm
  !
  ! check if inversion (I) is a symmetry.
  ! If so, it should be the (nsym/2+1)-th operation of the group
  !
  invsym = all ( s(:,:,nsym/2+1) == -s(:,:,1) )
  !
  CALL inverse_s ( )
  !
  CALL s_axis_to_cart ( )
  !
  RETURN
  !
END SUBROUTINE find_sym
!
!-----------------------------------------------------------------------
SUBROUTINE sgam_at ( nat, tau, ityp, nr1, nr2, nr3, sym, no_z_inv)
  !-----------------------------------------------------------------------
  !
  !     Given the point group of the Bravais lattice, this routine finds
  !     the subgroup which is the point group of the considered crystal.
  !     Non symmorphic groups are allowed, provided that fractional
  !     translations are allowed (nofrac=.false), that the unit cell is
  !     not a supercell, and that they are commensurate with the FFT grid
  !
  !     On output, the array sym is set to .true.. for each operation
  !     of the original point group that is also a symmetry operation
  !     of the crystal symmetry point group
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(in) :: nat, ityp (nat), nr1, nr2, nr3
  ! nat  : number of atoms in the unit cell
  ! ityp : species of each atom in the unit cell
  ! nr*  : dimensions of the FFT mesh
  !
  REAL(DP), INTENT(in) :: tau (3, nat)
  ! tau  : cartesian coordinates of the atoms
  !
  LOGICAL, INTENT(in), OPTIONAL :: no_z_inv
  ! no_z_inv: if .true., disable symmetry operations sending z into -z.
  !           Some calculations (e.g. monopole fields) require this
  !
  LOGICAL, INTENT(out) :: sym (48)
  ! sym(isym)    : flag indicating if sym.op. isym in the parent group
  !                is a true symmetry operation of the crystal
  !
  INTEGER :: na, kpol, nb, irot, i, j
  ! counters
  REAL(DP) , ALLOCATABLE :: xau (:,:), rau (:,:)
  ! atomic coordinates in crystal axis
  LOGICAL :: fractional_translations, no_z
  REAL(DP) :: ft_(3), ftaux(3)
  !
  ALLOCATE(xau(3,nat))
  ALLOCATE(rau(3,nat))
  !
  !     Compute the coordinates of each atom in the basis of
  !     the direct lattice vectors
  !
  DO na = 1, nat
     xau(:,na) = bg(1,:) * tau(1,na) + bg(2,:) * tau(2,na) + bg(3,:) * tau(3,na)
  ENDDO
  !
  !      check if the identity has fractional translations
  !      (this means that the cell is actually a supercell).
  !      When this happens, fractional translations are disabled,
  !      because there is no guarantee that the generated sym.ops.
  !      form a group
  !
  nb = 1
  irot = 1
  !
  fractional_translations = .not. nofrac
  IF ( fractional_translations ) THEN
     DO na = 2, nat
        IF (ityp (nb) == ityp (na) ) THEN
           !
           ft_(:) = xau(:,na) - xau(:,nb) - nint( xau(:,na) - xau(:,nb) )
           sym(irot) = checksym ( irot, nat, ityp, xau, xau, ft_ )
           IF (sym (irot) ) THEN
              fractional_translations = .false.
              WRITE( stdout, '(5x,"Found symmetry operation: I + (",&
             &   3f8.4, ")",/,5x,"This is a supercell,", &
             &   " fractional translations are disabled")') ft_
              GOTO 10
           ENDIF
           !
        ENDIF
     ENDDO
  ENDIF
  !
10 CONTINUE
  nsym_ns = 0
  DO irot = 1, nrot
     !
     ! check that the grid is compatible with the S rotation
     !
     IF ( mod (s (2, 1, irot) * nr1, nr2) /= 0 .or. &
          mod (s (3, 1, irot) * nr1, nr3) /= 0 .or. &
          mod (s (1, 2, irot) * nr2, nr1) /= 0 .or. &
          mod (s (3, 2, irot) * nr2, nr3) /= 0 .or. &
          mod (s (1, 3, irot) * nr3, nr1) /= 0 .or. &
          mod (s (2, 3, irot) * nr3, nr2) /= 0 ) THEN
        sym (irot) = .false.
        WRITE( stdout, '(5x,"warning: symmetry operation # ",i2, &
             &         " not compatible with FFT grid. ")') irot
        WRITE( stdout, '(3i4)') ( (s (i, j, irot) , j = 1, 3) , i = 1, 3)
        GOTO 20
     ENDIF

     DO na = 1, nat
        ! rau = rotated atom coordinates
        rau (:, na) = s (1,:, irot) * xau (1, na) + &
                      s (2,:, irot) * xau (2, na) + &
                      s (3,:, irot) * xau (3, na)
     ENDDO
     !
     !      first attempt: no fractional translation
     !
     ftau (:, irot) = 0
     ft (:, irot) = 0
     ft_(:) = 0.d0
     !
     sym(irot) = checksym ( irot, nat, ityp, xau, rau, ft_ )
     !
     IF (.not.sym (irot) .and. fractional_translations) THEN
        nb = 1
        DO na = 1, nat
           IF (ityp (nb) == ityp (na) ) THEN
              !
              !    second attempt: check all possible fractional translations
              !
              ft_ (:) = rau(:,na) - xau(:,nb) - nint( rau(:,na) - xau(:,nb) )
              !
              !    ft_ is in crystal axis and is a valid fractional translation
              !    only if ft_(i)=0 or ft_(i)=1/n, with n=2,3,4,6
              !    The check below is less strict: n must be integer
              !
              DO i=1,3
                 IF ( ABS (ft_(i)) > eps2 ) THEN
                    ftaux(i) = ABS (1.0_dp/ft_(i) - NINT(1.0_dp/ft_(i)) ) 
                 ELSE
                    ftaux(i) = 0.0_dp
                 END IF
              END DO
              IF ( ANY ( ftaux(:) > eps2 ) ) CYCLE
              !
              sym(irot) = checksym ( irot, nat, ityp, xau, rau, ft_ )
              !
              IF (sym (irot) ) THEN
                 nsym_ns = nsym_ns + 1
                 ft (:,irot) = ft_(:)
                 GOTO 20
              ENDIF
           ENDIF
        ENDDO

     ENDIF
20   CONTINUE
  ENDDO
  !
  ! convert ft to FFT coordinates, check if compatible with FFT grid
  ! for real-space symmetrization (if done: currently, exx, phonon)
  !
  nsym_na = 0
  DO irot =1, nrot
     IF ( sym(irot) .and. .not. allfrac ) THEN
        ftaux(1) = ft(1,irot) * nr1
        ftaux(2) = ft(2,irot) * nr2
        ftaux(3) = ft(3,irot) * nr3
        ! check if the fractional translations are commensurate
        ! with the FFT grid, discard sym.op. if not
        ! (needed because ph.x symmetrizes in real space)
        IF (abs (ftaux(1) - nint (ftaux(1)) ) / nr1 > eps2 .or. &
            abs (ftaux(2) - nint (ftaux(2)) ) / nr2 > eps2 .or. &
            abs (ftaux(3) - nint (ftaux(3)) ) / nr3 > eps2 ) THEN
            !     WRITE( stdout, '(5x,"warning: symmetry operation", &
            !          &     " # ",i2," not allowed.   fractional ", &
            !          &     "translation:"/5x,3f11.7,"  in crystal", &
            !          &     " coordinates")') irot, ft_
            sym (irot) = .false.
            nsym_na = nsym_na + 1
            nsym_ns = nsym_ns - 1
         ENDIF
         ftau (:, irot) = nint (ftaux(:))
      ENDIF
  ENDDO
  ! disable all symmetries z -> -z
  IF ( PRESENT(no_z_inv) ) THEN
     IF ( no_z_inv ) THEN
        DO irot=1,nrot
           IF (s(3,3,irot) == -1) sym(irot)=.false.
        END DO
     ENDIF
  ENDIF
  !
  !   deallocate work space
  !
  DEALLOCATE (rau)
  DEALLOCATE (xau)
  !
  RETURN
END SUBROUTINE sgam_at
!
!-----------------------------------------------------------------------
SUBROUTINE sgam_at_mag ( nat, m_loc, sym )
  !-----------------------------------------------------------------------
  !
  !   Find magnetic symmetries, i.e. point-group symmetries that are
  !   also symmetries of the local magnetization - including
  !   rotation + time reversal operations
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(in) :: nat
  REAL(DP), INTENT(in) :: m_loc(3, nat)
  !
  ! m_loc: local magnetization, must be invariant under the sym.op.
  !
  LOGICAL, INTENT(inout) :: sym (48)
  !
  ! sym(isym) = .true. if rotation isym is a sym.op. of the crystal
  !                    (i.e. not of the bravais lattice only)
  !
  INTEGER :: na, nb, irot
  LOGICAL :: t1, t2
  REAL(DP) , ALLOCATABLE ::  mxau(:,:), mrau(:,:)
  ! magnetization and rotated magnetization in crystal axis
  !
  ALLOCATE ( mxau(3,nat), mrau(3,nat) )
  !
  !     Compute the local magnetization of each atom in the basis of
  !     the direct lattice vectors
  !
  DO na = 1, nat
     mxau (:, na)= bg (1, :) * m_loc (1, na) + &
                   bg (2, :) * m_loc (2, na) + &
                   bg (3, :) * m_loc (3, na)
  ENDDO
  !
  DO irot = 1, nrot
     !
     t_rev(irot) = 0
     !
     IF ( sym (irot) ) THEN
        !
        ! mrau = rotated local magnetization
        !
        DO na = 1, nat
            mrau(:,na) = s(1,:,irot) * mxau(1,na) + &
                         s(2,:,irot) * mxau(2,na) + &
                         s(3,:,irot) * mxau(3,na)
        ENDDO
        IF (sname(irot)(1:3)=='inv') mrau = -mrau
        !
        ! check if this a magnetic symmetry
        !
        t1 = .true.
        t2 = .true.
        DO na = 1, nat
           !
           nb = irt (irot,na)
           IF ( nb < 1 .or. nb > nat ) CALL errore ('check_mag_sym', &
               'internal error: out-of-bound atomic index', na)
           !
           t1 = ( abs(mrau(1,na) - mxau(1,nb)) +       &
                  abs(mrau(2,na) - mxau(2,nb)) +       &
                  abs(mrau(3,na) - mxau(3,nb)) < eps2 ) .and. t1
           t2 = ( abs(mrau(1,na) + mxau(1,nb))+       &
                  abs(mrau(2,na) + mxau(2,nb))+       &
                  abs(mrau(3,na) + mxau(3,nb)) < eps2 ) .and. t2
           !
        ENDDO
        !
        IF ( .not.t1 .and. .not.t2 ) THEN
           ! not a magnetic symmetry
           sym(irot) = .false.
        ELSEIF( t2 .and. .not. t1 ) THEN
           ! magnetic symmetry with time reversal, if allowed
           IF (no_t_rev) THEN
              sym(irot) = .false.
           ELSE
              t_rev(irot) = 1
           ENDIF
        ENDIF
        IF ((.NOT. sym(irot)) .AND. (ftau(1,irot) /= 0 .OR. ftau(2,irot) /=0 &
                     .OR. ftau(3,irot) /=0)) nsym_ns=nsym_ns-1
        !
     ENDIF
     !
  ENDDO
  !
  !   deallocate work space
  !
  DEALLOCATE ( mrau, mxau )
  !
  RETURN
END SUBROUTINE sgam_at_mag
!
SUBROUTINE set_sym(nat, tau, ityp, nspin_mag, m_loc, nr1, nr2, nr3)
  !
  ! This routine receives as input atomic types and positions, if there
  ! is noncollinear magnetism and the initial magnetic moments, the fft
  ! dimensions nr1, nr2, nr3; it sets the symmetry elements of this module.
  ! Note that at and bg are those in cell_base. It sets nrot, nsym, s,
  ! sname, sr, invs, ftau, irt, t_rev,  time_reversal, and invsym
  !
  !-----------------------------------------------------------------------
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(in)  :: nat, ityp(nat), nspin_mag, nr1, nr2, nr3
  REAL(DP), INTENT(in) :: tau(3,nat)
  REAL(DP), INTENT(in) :: m_loc(3,nat)
  !
  time_reversal = (nspin_mag /= 4)
  t_rev(:) = 0
  CALL set_sym_bl ( )
  CALL find_sym ( nat, tau, ityp, nr1, nr2, nr3, .not.time_reversal, m_loc)
  !
  RETURN
  END SUBROUTINE set_sym
!
!-----------------------------------------------------------------------
INTEGER FUNCTION copy_sym ( nrot_, sym )
!-----------------------------------------------------------------------
  !
  IMPLICIT NONE
  INTEGER, INTENT(in) :: nrot_
  LOGICAL, INTENT(inout) :: sym(48)
  !
  INTEGER :: stemp(3,3), ftemp(3), ttemp, irot, jrot
  REAL(dp) :: ft_(3)
  INTEGER, ALLOCATABLE :: irtemp(:)
  CHARACTER(len=45) :: nametemp
  !
  ! copy symm. operations in sequential order so that
  ! s(i,j,irot) , irot <= nsym          are the sym.ops. of the crystal
  !               nsym+1 < irot <= nrot are the sym.ops. of the lattice
  ! on exit copy_sym returns nsym
  !
  ALLOCATE ( irtemp( size(irt,2) ) )
  jrot = 0
  DO irot = 1, nrot_
     IF (sym (irot) ) THEN
        jrot = jrot + 1
        IF ( irot > jrot ) THEN
           stemp = s(:,:,jrot)
           s (:,:, jrot) = s (:,:, irot)
           s (:,:, irot) = stemp
           ftemp(:) = ftau(:,jrot)
           ftau (:, jrot) = ftau (:, irot)
           ftau (:, irot) = ftemp(:)
           ft_(:) = ft(:,jrot)
           ft (:, jrot) = ft (:, irot)
           ft (:, irot) = ft_(:)
           irtemp (:) = irt (jrot,:)
           irt (jrot,:) = irt (irot,:)
           irt (irot,:) = irtemp (:)
           nametemp = sname (jrot)
           sname (jrot) = sname (irot)
           sname (irot) = nametemp
           ttemp = t_rev(jrot)
           t_rev(jrot) = t_rev(irot)
           t_rev(irot) = ttemp
        ENDIF
     ENDIF
  ENDDO
  sym (1:jrot) = .true.
  sym (jrot+1:nrot_) = .false.
  DEALLOCATE ( irtemp )
  !
  copy_sym = jrot
  RETURN
  !
END FUNCTION copy_sym

!
!-----------------------------------------------------------------------
LOGICAL FUNCTION is_group ( nsym_ )
  !-----------------------------------------------------------------------
  !
  !  Checks that {S} is a group
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(in) :: nsym_
  INTEGER :: isym, jsym, ksym, ss (3, 3)
  REAL (DP) :: st(3), dt(3)
  LOGICAL :: found
  !
  DO isym = 1, nsym_
     DO jsym = 1, nsym_
        !
        ss = matmul (s(:,:,isym),s(:,:,jsym))
        st(:)= ft(:,jsym) + s(1,:,jsym)*ft(1,isym) + &
                            s(2,:,jsym)*ft(2,isym) + &
                            s(3,:,jsym)*ft(3,isym)
        !
        !     here we check that the input matrices really form a group:
        !        S(k)   = S(i)*S(j)
        !        ftau_k = S(j)*ftau_i+ftau_j (modulo a lattice vector)
        !
        found = .false.
        DO ksym = 1, nsym_
           dt(:) = ft(:,ksym) - st(:) - nint( ft(:,ksym) - st(:) )
           IF ( all( s(:,:,ksym) == ss(:,:) ) .and. &
                ( abs ( dt(1) ) < eps2 ) .and. &
                ( abs ( dt(2) ) < eps2 ) .and. &
                ( abs ( dt(3) ) < eps2 ) ) THEN
              IF (found) THEN
                 is_group = .false.
                 RETURN
              ENDIF
              found = .true.
           ENDIF
        ENDDO
        IF ( .NOT. found ) THEN
           is_group = .false.
           RETURN
        ENDIF
     ENDDO
  ENDDO
  is_group=.true.
  RETURN
  !
END FUNCTION is_group
!
!-----------------------------------------------------------------------
LOGICAL FUNCTION checksym ( irot, nat, ityp, xau, rau, ft_ )
  !-----------------------------------------------------------------------
  !
  !   This function receives as input all the atomic positions xau,
  !   and the rotated rau by the symmetry operation ir. It returns
  !   true if for each atom na, it is possible to find an atom nb
  !   which is of the same type of na, and coincide with it after the
  !   symmetry operation. Fractional translations are allowed.
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(in) :: nat, ityp (nat), irot
  ! nat : number of atoms
  ! ityp: the type of each atom
  REAL(DP), INTENT(in) :: xau (3, nat), rau (3, nat), ft_(3)
  ! xau: the initial vectors (in crystal coordinates)
  ! rau: the rotated vectors (as above)
  ! ft_: fractionary translation (as above)
  !
  INTEGER :: na, nb
  LOGICAL, EXTERNAL :: eqvect
  ! the testing function
  !
  DO na = 1, nat
     DO nb = 1, nat
        IF( ityp (nb) == ityp (na) ) THEN
           checksym =  eqvect (rau (1, na), xau (1, nb), ft_ , accep) 
           IF ( checksym ) THEN
              !
              ! the rotated atom does coincide with one of the like atoms
              ! keep track of which atom the rotated atom coincides with
              !
              irt (irot, na) = nb
              GOTO 10
           ENDIF
        ENDIF
     ENDDO
     !
     ! the rotated atom does not coincide with any of the like atoms
     ! s(ir) + ft is not a symmetry operation
     !
     checksym=.FALSE.
     RETURN
10   CONTINUE
  ENDDO
  !
  ! s(ir) + ft is a symmetry operation
  !
  RETURN
END FUNCTION checksym
!
!-----------------------------------------------------------------------
SUBROUTINE checkallsym ( nat, tau, ityp, nr1, nr2, nr3 )
  !-----------------------------------------------------------------------
  !     given a crystal group this routine checks that the actual
  !     atomic positions and bravais lattice vectors are compatible with
  !     it. Used in relaxation/MD runs to check that atomic motion is
  !     consistent with assumed symmetry.
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(in) :: nat, ityp (nat), nr1, nr2, nr3
  REAL(DP), INTENT(in) :: tau (3, nat)
  !
  INTEGER :: na, kpol, isym, i, j, k, l
  LOGICAL :: loksym (48)
  REAL(DP) :: sx (3, 3), sy(3,3)
  REAL(DP) , ALLOCATABLE :: xau(:,:), rau(:,:)
  !
  ALLOCATE (xau( 3 , nat))
  ALLOCATE (rau( 3 , nat))
  !
  !     check that s(i,j, isym) is an orthogonal operation
  !
  DO isym = 1, nsym
     sx = dble( s(:,:,isym) )
     sy = matmul ( bg, sx )
     sx = matmul ( sy, transpose(at) )
     ! sx is s in cartesian axis
     sy = matmul ( transpose ( sx ), sx )
     ! sy = s*transpose(s) = I
     DO i = 1, 3
        sy (i,i) = sy (i,i) - 1.0_dp
     ENDDO
     IF (any (abs (sy) > eps1 ) ) &
         CALL errore ('checkallsym', 'not orthogonal operation', isym)
  ENDDO
  !
  !     Compute the coordinates of each atom in the basis of the lattice
  !
  DO na = 1, nat
     DO kpol = 1, 3
        xau (kpol, na) = bg (1, kpol) * tau (1, na) + &
                         bg (2, kpol) * tau (2, na) + &
                         bg (3, kpol) * tau (3, na)
     ENDDO
  ENDDO
  !
  !     generate the coordinates of the rotated atoms
  !
  DO isym = 1, nsym
     DO na = 1, nat
        DO kpol = 1, 3
           rau (kpol, na) = s (1, kpol, isym) * xau (1, na) + &
                            s (2, kpol, isym) * xau (2, na) + &
                            s (3, kpol, isym) * xau (3, na)
        ENDDO
     ENDDO
     !
     loksym(isym) =  checksym ( isym, nat, ityp, xau, rau, ft(1,isym) )
     !
  ENDDO
  !
  !   deallocate work space
  !
  DEALLOCATE(rau)
  DEALLOCATE(xau)
  !
  DO isym = 1,nsym
     IF (.not.loksym (isym) ) CALL errore ('checkallsym', &
          'the following symmetry operation is not satisfied  ', -isym)
  ENDDO
  IF (any (.not.loksym (1:nsym) ) ) THEN
      !call symmetrize_at (nsym, s, invs, ft, irt, nat, tau, at, bg, &
      !                    alat, omega)
      CALL errore ('checkallsym', &
           'some of the original symmetry operations not satisfied ',1)
  ENDIF
  !
  RETURN
END SUBROUTINE checkallsym

!----------------------------------------------------------------------
SUBROUTINE s_axis_to_cart ( )
  !----------------------------------------------------------------------
  !
  !     This routine transforms symmetry matrices expressed in the
  !     basis of the crystal axis into rotations in cartesian axis
  !
  IMPLICIT NONE
  !
  INTEGER :: isym
  REAL(dp):: sa(3,3), sb(3,3)
  !
  DO isym = 1,nsym
     sa (:,:) = dble ( s(:,:,isym) )
     sb = matmul ( bg, sa )
     sr (:,:, isym) = matmul ( at, transpose (sb) )
  ENDDO
  !
 END SUBROUTINE s_axis_to_cart

!-----------------------------------------------------------------------
SUBROUTINE find_sym_ifc ( nat, tau, ityp)
  !-----------------------------------------------------------------------
  !
  !     This routine finds the point group of the crystal, by eliminating
  !     the symmetries of the Bravais lattice which are not allowed
  !     by the atomic positions (for use in the FD package)
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(in) :: nat, ityp (nat)
  REAL(DP), INTENT(in) :: tau (3,nat)
  !
  INTEGER :: i
  LOGICAL :: sym (48)
  ! if true the corresponding operation is a symmetry operation
  !
  IF ( .not. allocated(irt) ) ALLOCATE( irt( 48, nat ) )
  irt( :, : ) = 0
  !
  !    Here we find the true symmetries of the crystal
  !
  CALL sgam_at_ifc ( nat, tau, ityp, sym )
  !
  !    Here we re-order all rotations in such a way that true sym.ops
  !    are the first nsym; rotations that are not sym.ops. follow
  !
  nsym = copy_sym ( nrot, sym )
  !
  ! check if inversion (I) is a symmetry.
  ! If so, it should be the (nsym/2+1)-th operation of the group
  !
  invsym = all ( s(:,:,nsym/2+1) == -s(:,:,1) )
  !
  CALL inverse_s ( )
  !
  CALL s_axis_to_cart ( )
  !
  RETURN
  !
END SUBROUTINE find_sym_ifc
!
!-----------------------------------------------------------------------
SUBROUTINE sgam_at_ifc ( nat, tau, ityp, sym )
  !-----------------------------------------------------------------------
  !
  !     Given the point group of the Bravais lattice, this routine finds
  !     the subgroup which is the point group of the considered crystal.
  !     Non symmorphic groups are allowed, provided that fractional
  !     translations are allowed (nofrac=.false), that the unit cell is
  !     not a supercell.
  !
  !     On output, the array sym is set to .true.. for each operation
  !     of the original point group that is also a symmetry operation
  !     of the crystal symmetry point group
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(in) :: nat, ityp (nat)
  ! nat  : number of atoms in the unit cell
  ! ityp : species of each atom in the unit cell
  !
  REAL(DP), INTENT(in) :: tau (3, nat)
  !
  ! tau  : cartesian coordinates of the atoms
  !
  !     output variables
  !
  LOGICAL, INTENT(out) :: sym (48)
  ! sym(isym)    : flag indicating if sym.op. isym in the parent group
  !                is a true symmetry operation of the crystal
  !
  INTEGER :: na, kpol, nb, irot, i, j
  ! counters
  REAL(DP) , ALLOCATABLE :: xau (:,:), rau (:,:)
  ! atomic coordinates in crystal axis
  LOGICAL :: fractional_translations
  REAL(DP) :: ft_(3)
  !
  ALLOCATE(xau(3,nat))
  ALLOCATE(rau(3,nat))
  !
  !     Compute the coordinates of each atom in the basis of
  !     the direct lattice vectors
  !
  DO na = 1, nat
     xau(:,na) = bg(1,:) * tau(1,na) + bg(2,:) * tau(2,na) + bg(3,:) * tau(3,na)
  ENDDO
  !
  !      check if the identity has fractional translations
  !      (this means that the cell is actually a supercell).
  !      When this happens, fractional translations are disabled,
  !      because there is no guarantee that the generated sym.ops.
  !      form a group
  !
  nb = 1
  irot = 1
  !
  fractional_translations = .not. nofrac
  DO na = 2, nat
     IF ( fractional_translations ) THEN
        IF (ityp (nb) == ityp (na) ) THEN
           ft_(:) = xau(:,na) - xau(:,nb) - nint( xau(:,na) - xau(:,nb) )
           !
           sym(irot) = checksym ( irot, nat, ityp, xau, xau, ft_ )
           !
           IF ( sym (irot) .and. &
               (abs (ft_(1) **2 + ft_(2) **2 + ft_(3) **2) < 1.d-8) ) &
               CALL errore ('sgam_at_ifc', 'overlapping atoms', na)
        ENDIF
     ENDIF
  ENDDO
  !
  nsym_ns = 0
  DO irot = 1, nrot
     DO na = 1, nat
        ! rau = rotated atom coordinates
        rau (:, na) = s (1,:, irot) * xau (1, na) + &
                      s (2,:, irot) * xau (2, na) + &
                      s (3,:, irot) * xau (3, na)
     ENDDO
     !
     !      first attempt: no fractional translation
     !
     ftau (:, irot) = 0
     ft (:, irot) = 0
     ft_(:) = 0.d0
     !
     sym(irot) = checksym ( irot, nat, ityp, xau, rau, ft_ )
     !
     IF (.not.sym (irot) .and. fractional_translations) THEN
        nb = 1
        DO na = 1, nat
           IF (ityp (nb) == ityp (na) ) THEN
              !
              !      second attempt: check all possible fractional translations
              !
              ft_ (:) = rau(:,na) - xau(:,nb) - nint( rau(:,na) - xau(:,nb) )
              !
              sym(irot) = checksym ( irot, nat, ityp, xau, rau, ft_ )
              !
              IF (sym (irot) ) THEN
                 nsym_ns = nsym_ns + 1
                 ft (:,irot) = ft_(:)
                 GOTO 100
              ENDIF
           ENDIF
        ENDDO

     ENDIF
100  CONTINUE
  ENDDO
  !
  DEALLOCATE (rau)
  DEALLOCATE (xau)
  !
  RETURN
END SUBROUTINE sgam_at_ifc

END MODULE symm_base
