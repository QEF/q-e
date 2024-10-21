!
! Copyright (C) 2010-2017 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!---------------------------------------------------------------------------------
!
MODULE symm_base
  !
  !! This module contains the variables needed to describe the symmetry properties
  !! and the routines to find crystal symmetries.
  !
  USE kinds,      ONLY : DP
  USE io_global,  ONLY : stdout
  USE cell_base,  ONLY : at, bg
  USE ions_base,  ONLY : atm
  USE noncollin_module, ONLY : colin_mag
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
  PUBLIC :: s, sr, sname, ft, nrot, nsym, nsym_ns, nsym_na, t_rev,  &
            no_t_rev, time_reversal, irt, invs, invsym, d1, d2, d3, &
            allfrac, nofrac, nosym, nosym_evc, fft_fact, spacegroup,&
            chem_symb
  !
  INTEGER :: s(3,3,48)
  !! symmetry matrices, in crystal axis
  INTEGER :: invs(48)
  !! index of inverse operation: S^{-1}_i=S(invs(i))
  INTEGER :: fft_fact(3) = 1
  !! FFT dimensions must be multiple of fft_fact
  INTEGER :: nrot
  !! number of bravais lattice symmetries
  INTEGER :: spacegroup = 0
  !! space group index, as read from input
  INTEGER :: nsym = 1
  !! total number of crystal symmetries
  INTEGER :: nsym_ns = 0
  !! nonsymmorphic (fractional translation) symms
  INTEGER :: nsym_na = 0 
  !! excluded nonsymmorphic symmetries because
  !! fract. transl. is noncommensurate with FFT grid
  REAL(DP) :: ft (3,48)
  !! fractional translations, in crystal axis
  REAL(DP) :: sr (3,3,48)
  !! symmetry matrices, in cartesian axis
  REAL(DP) :: accep = 1.0d-5
  !! initial value of the acceptance threshold
  !! for position comparison by eqvect in checksym
  !
  CHARACTER(LEN=45) :: sname(48)
  !! name of the symmetries
  INTEGER :: t_rev(48) = 0 
  !! time reversal flag, for noncolinear magnetism
  INTEGER, ALLOCATABLE :: irt(:,:)
  !! symmetric atom for each atom and sym.op.
  LOGICAL :: time_reversal = .TRUE.
  !! if .TRUE. the system has time reversal symmetry
  LOGICAL :: invsym
  !! if .TRUE. the system has inversion symmetry
  LOGICAL :: nofrac = .FALSE.
  !! if .TRUE. fract. translations are not allowed
  LOGICAL :: allfrac = .FALSE.
  !! if .TRUE. all fractionary translations allowed,
  !! even those not commensurate with FFT grid
  LOGICAL :: nosym = .FALSE.
  !! if .TRUE. no symmetry is used
  LOGICAL :: nosym_evc = .FALSE.
  !! if .TRUE. symmetry is used only to symmetrize
  !! k points
  LOGICAL :: no_t_rev = .FALSE.
  !! if .TRUE. remove the symmetries that
  !! require time reversal
  REAL(DP), TARGET :: d1(3,3,48)
  !! matrix d1 for rotation of spherical harmonics (d1 for l=1, ...)
  REAL(DP), TARGET :: d2(5,5,48)
  !! matrix d2 for rotation of spherical harmonics
  REAL(DP), TARGET :: d3(7,7,48)
  !! matrix d3 for rotation of spherical harmonics
  !
  ! ... Exported routines
  !
  PUBLIC ::  find_sym, inverse_s, copy_sym, checkallsym, &
             s_axis_to_cart, set_sym, set_sym_bl, check_grid_sym 
  PUBLIC ::  find_sym_ifc ! FIXME: should be merged with find_sym
  PUBLIC ::  remove_sym   ! FIXME: is this still useful?
  !
CONTAINS
   !
   !-----------------------------------------------------------------------
   SUBROUTINE inverse_s()
     !-----------------------------------------------------------------------
     !! Locate index of \(S^{-1}\).
     !
     IMPLICIT NONE
     !
     INTEGER :: isym, jsym, ss(3,3)
     LOGICAL :: found
     !
     DO isym = 1, nsym
        found = .FALSE.
        DO jsym = 1, nsym
           !
           ss = MATMUL( s(:,:,jsym), s(:,:,isym) )
           ! s(:,:,1) is the identity
           IF (ALL( s(:,:,1) == ss(:,:) )) THEN
              invs(isym) = jsym
              found = .TRUE.
           ENDIF
           !
        ENDDO
        IF (.NOT.found) CALL errore( 'inverse_s', ' Not a group', 1 )
     ENDDO
     !
   END SUBROUTINE inverse_s
   !
   !
   !-----------------------------------------------------------------------
   SUBROUTINE set_sym_bl()
     !---------------------------------------------------------------------
     !! Provides symmetry operations for all bravais lattices.  
     !! Tests the 24 proper rotations for the cubic lattice first, then
     !! the 8 rotations specific for the hexagonal axis (special axis c), 
     !! then inversion is added.
     !
     USE matrix_inversion
     !
     IMPLICIT NONE
     !
     CHARACTER(LEN=6), EXTERNAL :: int_to_char
     !
     ! sin3 = sin(pi/3), cos3 = cos(pi/3), msin3 = -sin(pi/3), mcos3 = -cos(pi/3)
     !
     REAL(DP), PARAMETER :: sin3 = 0.866025403784438597d0,  cos3 =  0.5d0, &
                           msin3 =-0.866025403784438597d0, mcos3 = -0.5d0
     !
     REAL(DP) :: s0(3,3,32), overlap(3,3), rat(3), rot(3,3), value
     ! s0: the s matrices in cartesian axis
     ! overlap: inverse overlap matrix between direct lattice
     ! rat: the rotated of a direct vector ( cartesian )
     ! rot: the rotated of a direct vector ( crystal axis )
     ! value: component of the s matrix in axis basis
     INTEGER :: jpol, kpol, mpol, irot, imat(32)
     ! counters over the polarizations and the rotations
     !
     CHARACTER(LEN=45) :: s0name(64)
     ! full name of the rotational part of each symmetry operation
     !
     DATA s0/ 1.d0,  0.d0,  0.d0,  0.d0,  1.d0,  0.d0,  0.d0,  0.d0,  1.d0, &
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
     !
     DATA s0name/  'identity                                     ',&
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
     !
     ! ... compute the overlap matrix for crystal axis
     DO jpol = 1, 3
        DO kpol = 1, 3
           rot(kpol,jpol) = at(1,kpol)*at(1,jpol) + &
                            at(2,kpol)*at(2,jpol) + &
                            at(3,kpol)*at(3,jpol)
        ENDDO
     ENDDO
     !
     ! ... then its inverse (rot is used as work space)
     CALL invmat( 3, rot, overlap )
     !
     nrot = 1
     !
     DO irot = 1, 32
        !
        ! ... for each possible symmetry
        DO jpol = 1, 3
           DO mpol = 1, 3
              !
              ! ... compute, in cartesian coordinates the rotated vector
              rat(mpol) = s0(mpol,1,irot)*at(1,jpol) + &
                          s0(mpol,2,irot)*at(2,jpol) + &
                          s0(mpol,3,irot)*at(3,jpol)
           ENDDO

           DO kpol = 1, 3
              !
              ! ... the rotated vector is projected on the direct lattice
              rot(kpol,jpol) = at(1,kpol)*rat(1) + &
                               at(2,kpol)*rat(2) + &
                               at(3,kpol)*rat(3)
           ENDDO
        ENDDO
        !
        ! ... and the inverse of the overlap matrix is applied
        DO jpol = 1,3
           DO kpol = 1,3
              value = overlap(jpol,1)*rot(1,kpol) + &
                      overlap(jpol,2)*rot(2,kpol) + &
                      overlap(jpol,3)*rot(3,kpol)
              !
              IF ( ABS(DBLE(NINT(value))-value) > eps1 ) THEN
                 !
                 ! ... if a noninteger is obtained, this implies that this operation
                 ! is not a symmetry operation for the given lattice
                 !
                 GOTO 10
              ENDIF
              !
              s(kpol,jpol,nrot) = NINT(value)
           ENDDO
        ENDDO
        !
        sname(nrot) = s0name(irot)
        imat(nrot) = irot
        nrot = nrot+1
        !
   10   CONTINUE
        !
     ENDDO
     !
     nrot = nrot-1
     !
     IF ( nrot /= 1 .AND. nrot /= 2 .AND. nrot /= 4 .AND. nrot /= 6 .AND. &
          nrot /= 8 .AND. nrot /=12 .AND. nrot /=24 ) THEN
          WRITE (stdout, '(80("-"),/,"NOTICE: Bravais lattice has wrong number (",&
         & i2,") of symmetries - symmetries are disabled",/,80("-"))' ) nrot
         nrot = 1
     ENDIF
     !
     ! ... set the inversion symmetry (Bravais lattices have always inversion symmetry)
     DO irot = 1, nrot
        sname(irot+nrot) = s0name(imat(irot)+32)
        DO kpol = 1, 3
           DO jpol = 1, 3
              s(kpol,jpol,irot+nrot) = -s(kpol,jpol,irot)
           ENDDO
        ENDDO
     ENDDO
     !
     nrot = 2*nrot
     !
     ! ... reset fractional translations to zero before checking the group
     ft(:,:) = 0.0_dp
     IF ( .NOT. is_group(nrot) ) THEN
        ! ... This happens for instance for an hexagonal lattice with one axis 
        ! oriented at 15 degrees from the x axis, the other along (-1,1,0)
        CALL infomsg( 'set_sym_bl', 'NOTICE: Symmetry group for Bravais lattice &
                     &is not a group (' // TRIM(int_to_char(nrot)) // &
                      ') - symmetries are disabled' )
         nrot = 1
     ENDIF
     !
     RETURN
     !
   END SUBROUTINE set_sym_bl
   !
   !
   !-----------------------------------------------------------------------
   SUBROUTINE find_sym( nat, tau, ityp, magnetic_sym, m_loc, no_z_inv )
     !-----------------------------------------------------------------------
     !! This routine finds the point group of the crystal, by eliminating
     !! the symmetries of the Bravais lattice which are not allowed
     !! by the atomic positions (or by the magnetization if present).
     !
     IMPLICIT NONE
     !
     INTEGER, INTENT(IN) :: nat
     !! total number of atoms of all species
     INTEGER, INTENT(IN) :: ityp(nat)
     !! the type of each i-th atom in stdin
     REAL(DP), INTENT(IN) :: tau(3,nat)
     !! atomic positions
     REAL(DP), INTENT(IN) :: m_loc(3,nat)
     !! local integrated magnetization
     LOGICAL, INTENT(IN) :: magnetic_sym
     !! magnetic_sym = noncolin .AND. domag 
     LOGICAL, INTENT(IN), OPTIONAL :: no_z_inv
     !! if .TRUE., disable symmetries sending z into -z.  
     !! Some calculations (e.g. gate fields) require this.
     !
     ! ... local variables
     !
     INTEGER :: i
     LOGICAL :: sym(48)
     ! if true the corresponding operation is a symmetry operation
     !
     IF ( .NOT. ALLOCATED(irt) ) ALLOCATE( irt(48,nat) )
     irt(:,:) = 0
     !
     !    Here we find the true symmetries of the crystal
     !
     symm: DO i = 1, 3 !emine: if it is not resolved in 3 steps it is sth else?
       IF ( PRESENT(no_z_inv) ) THEN
          CALL sgam_at( nat, tau, ityp, sym, no_z_inv )
       ELSE
          CALL sgam_at( nat, tau, ityp, sym )
       ENDIF
       !
       ! ... Here we check for magnetic symmetries
       IF ( magnetic_sym ) THEN
          CALL sgam_at_mag( nat, m_loc, sym )
       ! ... Here we check for time reversal symmetries for collinear systems
       ! NOTE: This check should be performed in the consistent way as in setup.f90
       ! However, we temporarily use this way not to change the interface
       ! until the structure of the code is fixed.
       ELSE IF (colin_mag >= 1) THEN
          CALL sgam_at_collin( nat, m_loc, sym )
       ! ... If nosym_evc is true from now on we do not use the symmetry any more
       ENDIF
       !
       IF (nosym_evc) THEN
          sym = .FALSE.
          sym(1) = .TRUE.
       ENDIF
       !
       ! ... Here we re-order all rotations in such a way that true sym.ops
       ! are the first nsym; rotations that are not sym.ops. follow
       nsym = copy_sym( nrot, sym )
       !
       IF ( .NOT. is_group( nsym ) ) THEN
          IF (i == 1) CALL infomsg( 'find_sym', &
                         'Not a group! Trying with lower acceptance parameter...' )
          accep = accep * 0.5d0
          IF (i == 3) THEN
            CALL infomsg( 'find_sym', 'Still not a group! symmetry disabled' )
            nsym = 1
          ENDIF
          CYCLE symm
       ELSE
          IF (i > 1) CALL infomsg( 'find_sym', 'Symmetry operations form a group' )
          EXIT symm
       ENDIF
     ENDDO symm
     !
     ! ... check if inversion (I) is a symmetry.
     ! If so, it should be the (nsym/2+1)-th operation of the group
     !
     invsym = ALL( s(:,:,nsym/2+1) == -s(:,:,1) )
     !
     CALL inverse_s()
     !
     CALL s_axis_to_cart()
     !
     RETURN
     !
   END SUBROUTINE find_sym
   !
   !
   !-----------------------------------------------------------------------
   SUBROUTINE sgam_at( nat, tau, ityp, sym, no_z_inv )
     !-----------------------------------------------------------------------
     !! Given the point group of the Bravais lattice, this routine finds
     !! the subgroup which is the point group of the considered crystal.  
     !! Non symmorphic groups are allowed, provided that fractional
     !! translations are allowed (nofrac=.false) and that the unit cell
     !! is not a supercell.
     !
     !! On output, the array sym is set to .TRUE.. for each operation
     !! of the original point group that is also a symmetry operation
     !! of the crystal symmetry point group.
     !
     IMPLICIT NONE
     !
     INTEGER, INTENT(IN) :: nat
     !! number of atoms in the unit cell
     INTEGER, INTENT(IN) :: ityp(nat)
     !! species of each atom in the unit cell
     REAL(DP), INTENT(IN) :: tau(3,nat)
     !! cartesian coordinates of the atoms
     LOGICAL, INTENT(IN), OPTIONAL :: no_z_inv
     !! if .TRUE., disable symmetry operations sending z into -z.  
     !! Some calculations (e.g. gate fields) require this
     LOGICAL, INTENT(OUT) :: sym(48)
     !! flag indicating if sym.op. isym in the parent group
     !! is a true symmetry operation of the crystal.
     !
     ! ... local variables
     !
     INTEGER :: na, kpol, nb, irot, i, j
     ! counters
     REAL(DP) , ALLOCATABLE :: xau(:,:), rau(:,:)
     ! atomic coordinates in crystal axis
     LOGICAL :: fractional_translations, no_z
     INTEGER :: nfrac
     REAL(DP) :: ft_(3), ftaux(3)
     !
     ALLOCATE( xau(3,nat) )
     ALLOCATE( rau(3,nat) )
     !
     ! ... Compute the coordinates of each atom in the basis of
     ! the direct lattice vectors
     DO na = 1, nat
        xau(:,na) = bg(1,:)*tau(1,na) + bg(2,:)*tau(2,na) + bg(3,:)*tau(3,na)
     ENDDO
     !
     ! ... check if the identity has fractional translations (this means
     ! that the cell is actually a supercell). When this happens, fractional
     ! translations are disabled, because there is no guarantee that the 
     ! generated sym.ops. form a group.
     !
     nb = 1
     irot = 1
     !
     fractional_translations = .NOT. nofrac
     !
     IF ( fractional_translations ) THEN
        DO na = 2, nat
            IF ( (colin_mag >= 0 .AND. chem_symb( atm(ityp(nb)) ) == chem_symb( atm(ityp(na)) ) ) & 
                 .OR. (colin_mag < 0 .AND. ityp(nb) == ityp(na) ) )THEN
            !IF ( ityp(nb) == ityp(na) ) THEN
               !
              ft_(:) = xau(:,na) - xau(:,nb) - NINT( xau(:,na) - xau(:,nb) )
              sym(irot) = checksym ( irot, nat, ityp, xau, xau, ft_ )
              IF (sym(irot)) THEN
                 fractional_translations = .FALSE.
                 WRITE( stdout, '(5x,"Found identity + (",&
                &   3f8.4, ") symmetry",/,5x,"This is a supercell,", &
                &   " fractional translations are disabled")') ft_
                 GOTO 10
              ENDIF
              !
           ENDIF
        ENDDO
     ENDIF
     !
   10 CONTINUE
     ! 
     nsym_ns = 0
     fft_fact(:) = 1
     !
     DO irot = 1, nrot
        !
        DO na = 1, nat
           ! rau = rotated atom coordinates
           rau(:,na) = s(1,:,irot) * xau(1,na) + &
                       s(2,:,irot) * xau(2,na) + &
                       s(3,:,irot) * xau(3,na)
        ENDDO
        !
        ! ... first attempt: no fractional translation
        ft(:,irot) = 0
        ft_(:) = 0.d0
        !
        sym(irot) = checksym( irot, nat, ityp, xau, rau, ft_ )
        !
        IF (.NOT.sym(irot) .AND. fractional_translations) THEN
           nb = 1
           DO na = 1, nat
               IF ( (colin_mag >= 0 .AND. chem_symb( atm(ityp(nb)) ) == chem_symb( atm(ityp(na)) ) ) & 
                    .OR. (colin_mag < 0 .AND. ityp(nb) == ityp(na) ) )THEN
               !IF ( ityp(nb) == ityp(na) ) THEN
                  !
                 ! ... second attempt: check all possible fractional translations
                 ft_(:) = rau(:,na) - xau(:,nb) - NINT( rau(:,na) - xau(:,nb) )
                 !
                 ! ... ft_ is in crystal axis and is a valid fractional translation
                 ! only if ft_(i)=0 or ft_(i)=1/n, with n=2,3,4,6
                 !
                 DO i = 1, 3
                    IF ( ABS (ft_(i)) > eps2 ) THEN
                       ftaux(i) = ABS (1.0_dp/ft_(i) - NINT(1.0_dp/ft_(i)) )
                       nfrac = NINT(1.0_dp/ABS(ft_(i)))
                       IF ( ftaux(i) < eps2 .AND. nfrac /= 2 .AND. &
                            nfrac /= 3 .AND. nfrac /= 4 .AND. nfrac /= 6 ) &
                            ftaux(i) = 2*eps2
                    ELSE
                       ftaux(i) = 0.0_dp
                    ENDIF
                 ENDDO
                 !
                 IF ( ANY( ftaux(:) > eps2 ) ) CYCLE
                 !
                 sym(irot) = checksym( irot, nat, ityp, xau, rau, ft_ )
                 !
                 IF ( sym(irot) ) THEN
                    nsym_ns = nsym_ns + 1
                    ft(:,irot) = ft_(:)
                    !
                    ! ... Find factors that must be present in FFT grid dimensions
                    ! in order to ensure that fractional translations are
                    ! commensurate with FFT grids.
                    DO i = 1, 3
                       IF ( ABS (ft_(i)) > eps2 ) THEN
                          nfrac = NINT(1.0_dp/ABS(ft_(i)))
                       ELSE
                          nfrac = 0
                       END IF
                       fft_fact(i) = mcm(fft_fact(i),nfrac)
                    ENDDO
                    !
                    GOTO 20
                 ENDIF
              ENDIF
           ENDDO
           !
        ENDIF
        !
   20   CONTINUE
        !
     ENDDO
     !
     ! ... disable all symmetries z -> -z
     IF ( PRESENT(no_z_inv) ) THEN
        IF ( no_z_inv ) THEN
           DO irot = 1, nrot
              IF (s(3,3,irot) == -1) sym(irot) = .FALSE.
           ENDDO
        ENDIF
     ENDIF
     !
     ! ... deallocate work space
     DEALLOCATE( rau )
     DEALLOCATE( xau )
     !
     RETURN
     !
   END SUBROUTINE sgam_at
   !
   !
   !-----------------------------------------------------------------------
   SUBROUTINE sgam_at_mag( nat, m_loc, sym )
     !-----------------------------------------------------------------------
     !! Find magnetic symmetries, i.e. point-group symmetries that are
     !! also symmetries of the local magnetization - including
     !! rotation + time reversal operations.
     !
     IMPLICIT NONE
     !
     INTEGER, INTENT(IN) :: nat
     !! numbero of atoms of all species
     REAL(DP), INTENT(IN) :: m_loc(3,nat)
     !! local magnetization, must be invariant under the sym.op.
     LOGICAL, INTENT(INOUT) :: sym(48)
     !! .TRUE. if rotation isym is a sym.op. of the crystal
     !! (i.e. not of the bravais lattice only)
     !
     ! ... local variables
     !
     INTEGER :: na, nb, irot
     LOGICAL :: t1, t2
     REAL(DP) , ALLOCATABLE ::  mxau(:,:), mrau(:,:)
     ! magnetization and rotated magnetization in crystal axis
     !
     ALLOCATE( mxau(3,nat), mrau(3,nat) )
     !
     ! ... Compute the local magnetization of each atom in the basis of
     ! the direct lattice vectors
     DO na = 1, nat
        mxau(:,na) = bg(1,:) * m_loc(1,na) + &
                     bg(2,:) * m_loc(2,na) + &
                     bg(3,:) * m_loc(3,na)
     ENDDO
     !
     DO irot = 1, nrot
        !
        t_rev(irot) = 0
        !
        IF ( sym(irot) ) THEN
           !
           ! ... mrau = rotated local magnetization
           DO na = 1, nat
               mrau(:,na) = s(1,:,irot) * mxau(1,na) + &
                            s(2,:,irot) * mxau(2,na) + &
                            s(3,:,irot) * mxau(3,na)
           ENDDO
           !
           IF (sname(irot)(1:3) == 'inv') mrau = -mrau
           !
           ! ... check if this a magnetic symmetry
           t1 = .TRUE.
           t2 = .TRUE.
           !
           DO na = 1, nat
              !
              nb = irt(irot,na)
              IF ( nb < 1 .OR. nb > nat ) CALL errore( 'check_mag_sym', &
                  'internal error: out-of-bound atomic index', na )
              !
              t1 = ( ABS(mrau(1,na) - mxau(1,nb)) +       &
                     ABS(mrau(2,na) - mxau(2,nb)) +       &
                     ABS(mrau(3,na) - mxau(3,nb)) < eps2 ) .AND. t1
              t2 = ( ABS(mrau(1,na) + mxau(1,nb)) +       &
                     ABS(mrau(2,na) + mxau(2,nb)) +       &
                     ABS(mrau(3,na) + mxau(3,nb)) < eps2 ) .AND. t2
              !
           ENDDO
           !
           IF ( .NOT. t1 .AND. .NOT. t2 ) THEN
              ! not a magnetic symmetry
              sym(irot) = .FALSE.
           ELSEIF( t2 .AND. .NOT. t1 ) THEN
              ! magnetic symmetry with time reversal, if allowed
              IF (no_t_rev) THEN
                 sym(irot) = .FALSE.
              ELSE
                 t_rev(irot) = 1
              ENDIF
           ENDIF
           IF ( (.NOT. sym(irot) ) .AND. &
                ( ABS(ft(1,irot)) > eps2 .OR. &
                  ABS(ft(2,irot)) > eps2 .OR. &
                  ABS(ft(3,irot)) > eps2 ) ) nsym_ns = nsym_ns-1
           !
        ENDIF
        !
     ENDDO
     !
     ! ... deallocate work space
     DEALLOCATE( mrau, mxau )
     !
     RETURN
     !
   END SUBROUTINE sgam_at_mag
   !
   !
   !-----------------------------------------------------------------------
   SUBROUTINE sgam_at_collin( nat, m_loc, sym )
      !-----------------------------------------------------------------------
      !! Find spin-space-group symmetries of the collinear system, i.e.
      !! the pair of the point-group symmetries and the spin operations
      !! that are symmetries of the atomic configurations and the local magnetization.
      !! The spin operations include the identity and the time reversal.
      !
      IMPLICIT NONE
      !
      INTEGER, INTENT(IN) :: nat
      !! numbero of atoms of all species
      REAL(DP), INTENT(IN) :: m_loc(3,nat)
      !! local magnetization, must be invariant under the sym.op.
      LOGICAL, INTENT(INOUT) :: sym(48)
      !! .TRUE. if rotation isym is a sym.op. of the crystal
      !! (i.e. not of the bravais lattice only)
      !
      ! ... local variables
      !
      INTEGER :: na, nb, irot
      LOGICAL :: t1, t2
      REAL(DP) , ALLOCATABLE ::  m_org(:), m_op(:)
      ! magnetization and rotated magnetization in crystal axis
      !
      ALLOCATE( m_org(nat), m_op(nat) )
      !
      ! Set original magnetization
      DO na = 1, nat
         m_org(na) = m_loc(3,na)
      ENDDO

      ! Check for time reversal
      DO irot = 1, nrot
         !
         t_rev(irot) = 0
         !
         IF ( sym(irot) ) THEN
            DO na = 1, nat
               nb = irt(irot,na)
               IF ( nb < 1 .OR. nb > nat ) CALL errore( 'check_mag_sym', &
                   'internal error: out-of-bound atomic index', na )

               m_op(nb) = m_org(na)

            ENDDO

            IF (ALL( ABS(m_op - m_org) < 1.0D-6)) THEN
               ! the operation is a symmetry without time-reversal
               t_rev(irot) = 0
            ELSE IF (ALL( ABS(m_op + m_org) < 1.0D-6)) THEN
               IF ( colin_mag == 1) THEN 
                  ! discard symmteries with time-reversal
                  sym(irot) = .FALSE. 
               ELSE ! IF ( colin_mag == 2) THEN
                  ! the operation is a symmetry with time-reversal
                  t_rev(irot) = 1
               ENDIF
            ELSE
               ! the operation is not a symmetry
               sym(irot) = .FALSE.
            ENDIF


            IF ( (.NOT. sym(irot) ) .AND. &
               ( ABS(ft(1,irot)) > eps2 .OR. &
                 ABS(ft(2,irot)) > eps2 .OR. &
                 ABS(ft(3,irot)) > eps2 ) ) nsym_ns = nsym_ns-1

         ENDIF
      ENDDO
      ! ... deallocate work space
      DEALLOCATE( m_op, m_org )
      !
      RETURN
      !
   END SUBROUTINE sgam_at_collin
   !
   !
   !-------------------------------------------------------------------------
   SUBROUTINE set_sym( nat, tau, ityp, nspin_mag, m_loc )
     !-----------------------------------------------------------------------
     !! This routine receives as input atomic types and positions, if there
     !! is noncollinear magnetism and the initial magnetic moments 
     !! it sets the symmetry elements of this module.  
     !! Note that \(at\) and \(bg\) are those in \(\textrm{cell_base}\). It sets nrot, nsym, s,
     !! sname, sr, invs, ft, irt, t_rev,  time_reversal, and invsym.
     !
     IMPLICIT NONE
     !
     INTEGER, INTENT(IN) :: nat
     !! number of atoms in the unit cell
     INTEGER, INTENT(IN) :: ityp(nat)
     !! species of each atom in the unit cell
     INTEGER, INTENT(IN) :: nspin_mag
     !! =1 when nspin=1,4 (domag=.false.), =2 when
     !! nspin=2, =4 nspin=4 (domag=.true.)
     REAL(DP), INTENT(IN) :: tau(3,nat)
     !! cartesian coordinates of the atoms
     REAL(DP), INTENT(IN) :: m_loc(3,nat)
     !! local magnetization, must be invariant under the sym.op.
     !
     time_reversal = (nspin_mag /= 4)
     t_rev(:) = 0
     !
     CALL set_sym_bl()
     CALL find_sym( nat, tau, ityp, .NOT.time_reversal, m_loc )
     !
     RETURN
     !
   END SUBROUTINE set_sym
   !
   !
   !-----------------------------------------------------------------------
   INTEGER FUNCTION copy_sym( nrot_, sym )
     !-----------------------------------------------------------------------
     !! Copy symmetry operations in sequential order so that:
     !
     !! * \(s(i,j,\text{irot})\), with \(\text{irot} \leq \text{nsym}\) are the symmetry
     !!   operations of the crystal;
     !! * \(s(i,j,\text{irot})\), with \(\text{nsym}+1<\text{irot}\leq \text{nrot}\) are 
     !!   the symmetry operations of the lattice.
     !
     !! On exit \(\textrm{copy_sym}\) returns nsym.
     !
     IMPLICIT NONE
     !
     INTEGER, INTENT(IN) :: nrot_
     !! number of rotations
     LOGICAL, INTENT(INOUT) :: sym(48)
     !! .TRUE. if rotation isym is a sym.op. of the crystal
     !! (i.e. not of the bravais lattice only)
     !
     ! ... local variables
     !
     INTEGER :: stemp(3,3), ftemp(3), ttemp, irot, jrot
     REAL(DP) :: ft_(3)
     INTEGER, ALLOCATABLE :: irtemp(:)
     CHARACTER(LEN=45) :: nametemp
     !
     !
     ALLOCATE ( irtemp(SIZE(irt,2)) )
     !
     jrot = 0
     !
     DO irot = 1, nrot_
        IF ( sym(irot) ) THEN
           jrot = jrot + 1
           IF (irot > jrot) THEN
              stemp = s(:,:,jrot)
              s(:,:,jrot) = s(:,:,irot)
              s(:,:,irot) = stemp
              ft_(:) = ft(:,jrot)
              ft(:,jrot) = ft(:,irot)
              ft(:,irot) = ft_(:)
              irtemp(:) = irt(jrot,:)
              irt(jrot,:) = irt(irot,:)
              irt(irot,:) = irtemp (:)
              nametemp = sname(jrot)
              sname(jrot) = sname(irot)
              sname(irot) = nametemp
              ttemp = t_rev(jrot)
              t_rev(jrot) = t_rev(irot)
              t_rev(irot) = ttemp
           ENDIF
        ENDIF
     ENDDO
     !
     sym(1:jrot) = .TRUE.
     sym(jrot+1:nrot_) = .FALSE.
     !
     DEALLOCATE( irtemp )
     !
     copy_sym = jrot
     !
     RETURN
     !
   END FUNCTION copy_sym
   !
   !
   !-----------------------------------------------------------------------
   LOGICAL FUNCTION is_group( nsym_ )
     !-----------------------------------------------------------------------
     !! Checks that {S} is a group.
     !
     IMPLICIT NONE
     !
     INTEGER, INTENT(IN) :: nsym_
     INTEGER :: isym, jsym, ksym, ss (3,3)
     REAL(DP) :: st(3), dt(3)
     LOGICAL :: found
     !
     DO isym = 1, nsym_
        DO jsym = 1, nsym_
           !
           ss = MATMUL(s(:,:,isym), s(:,:,jsym))
           st(:) = ft(:,jsym) + s(1,:,jsym)*ft(1,isym) + &
                                s(2,:,jsym)*ft(2,isym) + &
                                s(3,:,jsym)*ft(3,isym)
           !
           ! ... here we check that the input matrices really form a group:
           ! S(k) = S(i)*S(j)
           ! ftau_k = S(j)*ftau_i+ftau_j (modulo a lattice vector)
           !
           found = .FALSE.
           !
           DO ksym = 1, nsym_
              dt(:) = ft(:,ksym) - st(:) - NINT(ft(:,ksym) - st(:))
              IF ( ALL( s(:,:,ksym) == ss(:,:) ) .AND. &
                   ( ABS( dt(1) ) < eps2 ) .AND. &
                   ( ABS( dt(2) ) < eps2 ) .AND. &
                   ( ABS( dt(3) ) < eps2 ) ) THEN
                 IF (found) THEN
                    is_group = .FALSE.
                    RETURN
                 ENDIF
                 found = .TRUE.
              ENDIF
           ENDDO
           !
           IF ( .NOT. found ) THEN
              is_group = .FALSE.
              RETURN
           ENDIF
           !
        ENDDO
     ENDDO
     !
     is_group=.TRUE.
     !
     RETURN
     !
   END FUNCTION is_group
   !
   !
   !-----------------------------------------------------------------------
   LOGICAL FUNCTION checksym( irot, nat, ityp, xau, rau, ft_ )
     !-----------------------------------------------------------------------
     !! This function receives as input all the atomic positions xau,
     !! and the rotated rau by the symmetry operation ir. It returns
     !! .TRUE. if, for each atom na, it is possible to find an atom nb
     !! which is of the same type of na, and coincides with it after the
     !! symmetry operation. Fractional translations are allowed.
     !
     IMPLICIT NONE
     !
     INTEGER, INTENT(IN) :: nat
     !! number of atoms
     INTEGER, INTENT(IN) :: ityp(nat)
     !! the type of each atom
     INTEGER, INTENT(IN) :: irot
     !! rotation index
     REAL(DP), INTENT(IN) :: xau(3,nat)
     !! the initial vectors (in crystal coordinates)
     REAL(DP), INTENT(IN) :: rau(3,nat)
     !! the rotated vectors (as above)
     REAL(DP), INTENT(IN) :: ft_(3)
     !! fractionary translation (as above)
     !
     ! ... local variables
     !
     INTEGER :: na, nb
     LOGICAL, EXTERNAL :: eqvect
     ! the testing function
     !
     DO na = 1, nat
        DO nb = 1, nat
           !
            IF ( (colin_mag >= 0 .AND. chem_symb( atm(ityp(nb)) ) == chem_symb( atm(ityp(na)) ) ) & 
                 .OR. (colin_mag < 0 .AND. ityp(nb) == ityp(na) ) )THEN
            !IF ( ityp(nb) == ityp(na) ) THEN
               checksym =  eqvect( rau(1,na), xau(1,nb), ft_ , accep )
              IF ( checksym ) THEN
                 !
                 ! ... the rotated atom does coincide with one of the like atoms
                 ! keep track of which atom the rotated atom coincides with
                 irt (irot, na) = nb
                 GOTO 10
                 !
              ENDIF
           ENDIF
           !
        ENDDO
        !
        ! ... the rotated atom does not coincide with any of the like atoms
        ! s(ir) + ft is not a symmetry operation
        checksym = .FALSE.
        RETURN
        !
   10   CONTINUE
     ENDDO
     !
     ! ... s(ir) + ft is a symmetry operation
     !
     RETURN
     !
   END FUNCTION checksym
   !
   !
   !-----------------------------------------------------------------------
   SUBROUTINE checkallsym( nat, tau, ityp )
     !-----------------------------------------------------------------------
     !! Given a crystal group this routine checks that the actual atomic
     !! positions and bravais lattice vectors are compatible with it.  
     !! Used in relaxation/MD runs to check that atomic motion is
     !! consistent with assumed symmetry.
     !
     IMPLICIT NONE
     !
     INTEGER, INTENT(IN) :: nat
     !! number of atoms
     INTEGER, INTENT(IN) :: ityp(nat)
     !! the type of each atom
     REAL(DP), INTENT(IN) :: tau(3,nat)
     !! postion of each atom
     !
     ! ... local variables
     !
     INTEGER :: na, kpol, isym, i, j, k, l
     LOGICAL :: loksym(48)
     REAL(DP) :: sx(3,3), sy(3,3)
     REAL(DP), ALLOCATABLE :: xau(:,:), rau(:,:)
     !
     ALLOCATE( xau(3,nat) )
     ALLOCATE( rau(3,nat) )
     !
     ! ... check that s(i,j, isym) is an orthogonal operation
     DO isym = 1, nsym
        sx = DBLE( s(:,:,isym) )
        sy = MATMUL( bg, sx )
        sx = MATMUL( sy, TRANSPOSE(at) )
        ! sx is s in cartesian axis
        sy = MATMUL( TRANSPOSE(sx), sx )
        ! sy = s*TRANSPOSE(s) = I
        DO i = 1, 3
           sy(i,i) = sy(i,i) - 1.0_dp
        ENDDO
        IF ( ANY(ABS(sy) > eps1) ) &
            CALL errore( 'checkallsym', 'not orthogonal operation', isym )
     ENDDO
     !
     ! ... Compute the coordinates of each atom in the basis of the lattice
     DO na = 1, nat
        DO kpol = 1, 3
           xau(kpol,na) = bg(1,kpol) * tau(1,na) + &
                          bg(2,kpol) * tau(2,na) + &
                          bg(3,kpol) * tau(3,na)
        ENDDO
     ENDDO
     !
     ! ... Generate the coordinates of the rotated atoms
     DO isym = 1, nsym
        DO na = 1, nat
           DO kpol = 1, 3
              rau(kpol,na) = s(1,kpol,isym) * xau(1,na) + &
                             s(2,kpol,isym) * xau(2,na) + &
                             s(3,kpol,isym) * xau(3,na)
           ENDDO
        ENDDO
        !
        loksym(isym) =  checksym( isym, nat, ityp, xau, rau, ft(1,isym) )
     ENDDO
     !
     ! ... deallocate work space
     !
     DEALLOCATE( rau )
     DEALLOCATE( xau )
     !
     DO isym = 1,nsym
        IF ( .NOT.loksym(isym) ) CALL errore( 'checkallsym', &
             'the following symmetry operation is not satisfied  ', -isym )
     ENDDO
     !
     IF (ANY(.NOT.loksym(1:nsym) ) ) THEN
         !call symmetrize_at (nsym, s, invs, ft, irt, nat, tau, at, bg, &
         !                    alat, omega)
         CALL errore( 'checkallsym', &
              'some of the original symmetry operations not satisfied ',1 )
     ENDIF
     !
     RETURN
     !
   END SUBROUTINE checkallsym
   !
   !
   !----------------------------------------------------------------------
   SUBROUTINE s_axis_to_cart()
     !----------------------------------------------------------------------
     !! This routine transforms symmetry matrices expressed in the
     !! basis of the crystal axis into rotations in cartesian axis.
     !
     IMPLICIT NONE
     !
     INTEGER :: isym
     REAL(DP) :: sa(3,3), sb(3,3)
     !
     DO isym = 1,nsym
        sa(:,:) = DBLE( s(:,:,isym) )
        sb = MATMUL( bg, sa )
        sr(:,:,isym) = MATMUL( at, TRANSPOSE(sb) )
     ENDDO
     !
    END SUBROUTINE s_axis_to_cart
    !
    !
    !-----------------------------------------------------------------------
    SUBROUTINE find_sym_ifc( nat, tau, ityp )
      !-----------------------------------------------------------------------
      !! This routine finds the point group of the crystal, by eliminating
      !! the symmetries of the Bravais lattice which are not allowed
      !! by the atomic positions (for use in the FD package).
      !
      IMPLICIT NONE
      !
      INTEGER, INTENT(IN) :: nat
      !! number of atoms
      INTEGER, INTENT(IN) :: ityp(nat)
      !! the type of each atom
      REAL(DP), INTENT(IN) :: tau(3,nat)
      !! postion of each atom
      !
      ! ... local variables
      !
      INTEGER :: i
      LOGICAL :: sym(48)
      ! if true the corresponding operation is a symmetry operation
      !
      IF ( .NOT. ALLOCATED(irt) ) ALLOCATE( irt(48,nat) )
      irt(:,:) = 0
      !
      ! ... Here we find the true symmetries of the crystal
      CALL sgam_at_ifc( nat, tau, ityp, sym )
      !
      ! ... Here we re-order all rotations in such a way that true sym.ops
      ! are the first nsym; rotations that are not sym.ops. follow
      nsym = copy_sym( nrot, sym )
      !
      ! ... check if inversion (I) is a symmetry.
      ! If so, it should be the (nsym/2+1)-th operation of the group
      invsym = ALL( s(:,:,nsym/2+1) == -s(:,:,1) )
      !
      CALL inverse_s()
      !
      CALL s_axis_to_cart()
      !
      RETURN
      !
    END SUBROUTINE find_sym_ifc
    !
    !
    !-----------------------------------------------------------------------
    SUBROUTINE sgam_at_ifc( nat, tau, ityp, sym )
      !-----------------------------------------------------------------------
      !! Given the point group of the Bravais lattice, this routine finds
      !! the subgroup which is the point group of the considered crystal.
      !! Non symmorphic groups are allowed, provided that fractional
      !! translations are allowed (nofrac=.false), that the unit cell is
      !! not a supercell.
      !
      !! On output, the array sym is set to .TRUE.. for each operation
      !! of the original point group that is also a symmetry operation
      !! of the crystal symmetry point group.
      !
      IMPLICIT NONE
      !
      INTEGER, INTENT(IN) :: nat
      !! number of atoms in the unit cell
      INTEGER, INTENT(IN) :: ityp(nat)
      !! species of each atom in the unit cell
      REAL(DP), INTENT(IN) :: tau(3,nat)
      !! cartesian coordinates of the atoms
      LOGICAL, INTENT(OUT) :: sym(48)
      !! flag indicating if sym.op. isym in the parent group
      !! is a true symmetry operation of the crystal
      !
      ! ... local variables
      !
      INTEGER :: na, kpol, nb, irot, i, j
      ! counters
      REAL(DP) , ALLOCATABLE :: xau(:,:), rau(:,:)
      ! atomic coordinates in crystal axis
      LOGICAL :: fractional_translations
      REAL(DP) :: ft_(3)
      !
      ALLOCATE( xau(3,nat) )
      ALLOCATE( rau(3,nat) )
      !
      ! ... Compute the coordinates of each atom in the basis of
      ! the direct lattice vectors.
      !
      DO na = 1, nat
         xau(:,na) = bg(1,:)*tau(1,na) + bg(2,:)*tau(2,na) + bg(3,:)*tau(3,na)
      ENDDO
      !
      ! ... check if the identity has fractional translations
      ! (this means that the cell is actually a supercell).
      ! When this happens, fractional translations are disabled,
      ! because there is no guarantee that the generated sym.ops.
      ! form a group.
      !
      nb = 1
      irot = 1
      !
      fractional_translations = .NOT. nofrac
      !
      DO na = 2, nat
         IF ( fractional_translations ) THEN
            IF ( (colin_mag >= 0 .AND. chem_symb( atm(ityp(nb)) ) == chem_symb( atm(ityp(na)) ) ) & 
                 .OR. (colin_mag < 0 .AND. ityp(nb) == ityp(na) ) )THEN

            !IF ( ityp(nb) == ityp(na) ) THEN
               ft_(:) = xau(:,na) - xau(:,nb) - NINT( xau(:,na) - xau(:,nb) )
               !
               sym(irot) = checksym( irot, nat, ityp, xau, xau, ft_ )
               !
               IF ( sym (irot) .AND. &
                   (ABS(ft_(1)**2 + ft_(2)**2 + ft_(3)**2) < 1.d-8) ) &
                   CALL errore( 'sgam_at_ifc', 'overlapping atoms', na )
            ENDIF
         ENDIF
      ENDDO
      !
      nsym_ns = 0
      !
      DO irot = 1, nrot
         DO na = 1, nat
            ! rau = rotated atom coordinates
            rau(:,na) = s(1,:,irot) * xau(1,na) + &
                        s(2,:,irot) * xau(2,na) + &
                        s(3,:,irot) * xau(3,na)
         ENDDO
         !
         ! ... first attempt: no fractional translation
         ft(:,irot) = 0
         ft_(:) = 0.d0
         !
         sym(irot) = checksym( irot, nat, ityp, xau, rau, ft_ )
         !
         IF (.NOT.sym(irot) .AND. fractional_translations) THEN
            nb = 1
            !
            DO na = 1, nat
               IF ( (colin_mag >= 0 .AND. chem_symb( atm(ityp(nb)) ) == chem_symb( atm(ityp(na)) ) ) & 
                    .OR. (colin_mag < 0 .AND. ityp(nb) == ityp(na) ) )THEN
               !IF ( ityp(nb) == ityp(na) ) THEN
                  !
                  !      second attempt: check all possible fractional translations
                  !
                  ft_(:) = rau(:,na) - xau(:,nb) - NINT( rau(:,na) - xau(:,nb) )
                  !
                  sym(irot) = checksym( irot, nat, ityp, xau, rau, ft_ )
                  !
                  IF (sym(irot) ) THEN
                     nsym_ns = nsym_ns + 1
                     ft(:,irot) = ft_(:)
                     GOTO 100
                  ENDIF
               ENDIF
            ENDDO
            !
         ENDIF
         !
    100  CONTINUE
         !
      ENDDO
      !
      DEALLOCATE( rau )
      DEALLOCATE( xau )
      !
      RETURN
      !
    END SUBROUTINE sgam_at_ifc
    !
    !-----------------------------------------------------------------------
    FUNCTION check_grid_sym ( nr1, nr2, nr3 ) RESULT ( compatible )
      !---------------------------------------------------------------------
      !! Check that symmetry operations and FFT grid are compatible
      !! Needed to prevent trouble with real-space symmetrization
      !
      IMPLICIT NONE
      !
      INTEGER, INTENT(IN) :: nr1, nr2, nr3
      LOGICAL :: compatible, bad
      INTEGER :: isym,i,j
      !
      compatible = .true.
      DO isym = 1, nsym
         !
         bad = ( MOD( s(2,1,isym)*nr1, nr2) /= 0 .OR. &
                 MOD( s(3,1,isym)*nr1, nr3) /= 0 .OR. &
                 MOD( s(1,2,isym)*nr2, nr1) /= 0 .OR. &
                 MOD( s(3,2,isym)*nr2, nr3) /= 0 .OR. &
                 MOD( s(1,3,isym)*nr3, nr1) /= 0 .OR. &
                 MOD( s(2,3,isym)*nr3, nr2) /= 0 ) 
         IF ( bad ) THEN
            WRITE( stdout, '(5x,"warning: symmetry operation # ",i2, &
                 &         " not compatible with FFT grid. ")') isym
            WRITE( stdout, '(3i4)') ( (s(i,j,isym), j=1,3), i=1,3 )
            compatible = .false.
         ENDIF
         !
      ENDDO
      !
    END FUNCTION check_grid_sym
    !
    !-----------------------------------------------------------------------
    SUBROUTINE remove_sym( nr1, nr2, nr3 )
      !---------------------------------------------------------------------
      !! Compute ftau used for symmetrization in real space (phonon, exx)
      !! ensure that they are commensurated with the FFT grid.
      !
      IMPLICIT NONE
      !
      INTEGER, INTENT(IN) :: nr1, nr2, nr3
      !
      ! ... local variables
      !
      LOGICAL :: sym(48)
      INTEGER :: isym, nsym_, i, j
      REAL(dp) :: ftaux(3)
      !
      nsym_ = nsym
      sym(1:nsym_) = .TRUE.
      nsym_na = 0
      !
      DO isym = 1, nsym_
         !
         ! check that the grid is compatible with the S rotation
         !
         IF ( MOD( s(2,1,isym)*nr1, nr2) /= 0 .OR. &
              MOD( s(3,1,isym)*nr1, nr3) /= 0 .OR. &
              MOD( s(1,2,isym)*nr2, nr1) /= 0 .OR. &
              MOD( s(3,2,isym)*nr2, nr3) /= 0 .OR. &
              MOD( s(1,3,isym)*nr3, nr1) /= 0 .OR. &
              MOD( s(2,3,isym)*nr3, nr2) /= 0 ) THEN
            sym(isym) = .FALSE.
            WRITE( stdout, '(5x,"warning: symmetry operation # ",i2, &
                 &         " not compatible with FFT grid. ")') isym
            WRITE( stdout, '(3i4)') ( (s(i,j,isym), j=1,3), i=1,3 )
            sym(isym) = .FALSE.
            IF ( ABS(ft(1,isym)) > eps2 .OR. &
                 ABS(ft(2,isym)) > eps2 .OR. &
                 ABS(ft(3,isym)) > eps2 ) nsym_ns = nsym_ns-1
         ENDIF
         !
         ! convert ft to FFT coordinates, check if compatible with FFT grid
         ! for real-space symmetrization 
         !
         ftaux(1) = ft(1,isym) * nr1
         ftaux(2) = ft(2,isym) * nr2
         ftaux(3) = ft(3,isym) * nr3
         ! check if the fractional translations are commensurate
         ! with the FFT grid, discard sym.op. if not
         ! (needed because ph.x symmetrizes in real space)
         IF ( ABS(ftaux(1) - NINT(ftaux(1)) ) / nr1 > eps2 .OR. &
              ABS(ftaux(2) - NINT(ftaux(2)) ) / nr2 > eps2 .OR. &
              ABS(ftaux(3) - NINT(ftaux(3)) ) / nr3 > eps2 ) THEN
            !     WRITE( stdout, '(5x,"warning: symmetry operation", &
            !          &     " # ",i2," not allowed.   fractional ", &
            !          &     "translation:"/5x,3f11.7,"  in crystal", &
            !          &     " coordinates")') isym, ft_
            sym(isym) = .FALSE.
            nsym_na = nsym_na + 1
            nsym_ns = nsym_ns - 1
         ENDIF
         !
      ENDDO
      !
      ! ... count symmetries, reorder them exactly as in "find_sym" 
      !
      nsym = copy_sym( nsym_, sym )
      invsym = ALL( s(:,:,nsym/2+1) == -s(:,:,1) )
      !
      CALL inverse_s()
      CALL s_axis_to_cart()
      !
    END SUBROUTINE remove_sym
    !
    !
    !--------------------------------------------------------------------------
    INTEGER FUNCTION mcm( i, j )
      !------------------------------------------------------------------------
      !! Returns minimum common multiple of two integers
      !! if i=0, returns j, and vice versa; if i<0 or j<0, returns -1.
      !
      INTEGER, INTENT(IN) :: i,j
      INTEGER :: n1,n2,k
      !
      IF (i < 0 .OR. j < 0) THEN
         mcm = -1
      ELSEIF (i == 0 .AND. j == 0) THEN
         mcm = 0
      ELSE
         n1 = MIN(i,j)
         n2 = MAX(i,j)
         DO k = 1, n1
            mcm = k*n2 
            IF (MOD(mcm,n1) == 0) RETURN
         ENDDO
         mcm = n2
      ENDIF
      !
    END FUNCTION mcm
    !
    !
    !--------------------------------------------------------------------------
    CHARACTER FUNCTION chem_symb( symbol )
      !------------------------------------------------------------------------
      !! Returns the chemical symbol used to identify the symmetry
      !
      IMPLICIT NONE
      !
      CHARACTER(LEN=*), INTENT(IN) :: symbol
      !
      IF ( SCAN( symbol ,"0123456789") == 0 ) THEN
         chem_symb = symbol
      ELSE
         chem_symb = symbol( 1:SCAN( symbol , "0123456789_-" ) -1 )
      ENDIF
      !
    END FUNCTION chem_symb
    !
    !

END MODULE symm_base
