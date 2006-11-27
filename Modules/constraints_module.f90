!
! Copyright (C) 2002-2005 Quantum-ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "f_defs.h"
!
#define __REMOVE_CONSTRAINT_FORCE
!#define __DEBUG_CONSTRAINTS
#define __USE_PBC
!
!----------------------------------------------------------------------------
MODULE constraints_module
  !----------------------------------------------------------------------------
  !
  ! ... variables and methods for constraint Molecular Dynamics and
  ! ... constrained ionic relaxations (the SHAKE algorithm based on 
  ! ... lagrange multipliers) are defined here.
  !
  ! ... most of these variables and methods are also used for meta-dynamics
  ! ... and free-energy smd : indeed the collective variables are implemented
  ! ... as constraints.
  !
  ! ... written by Carlo Sbraccia ( 24/02/2004 )
  !
  ! ... references :
  !
  ! ... 1) M. P. Allen and D. J. Tildesley, Computer Simulations of Liquids,
  ! ...    Clarendon Press - Oxford (1986)
  !
  USE kinds,     ONLY : DP
  USE constants, ONLY : eps8, eps16, eps32, tpi, fpi
  USE io_global, ONLY : stdout
  !
  USE basic_algebra_routines
  !
  IMPLICIT NONE
  !
  SAVE
  !
  PRIVATE
  !
  ! ... public methods
  !
  PUBLIC :: init_constraint,       &
            check_constraint,      &
            remove_constr_force,   &
            deallocate_constraint, &
            compute_dmax,          &
            pbc
  !
  ! ... public variables (assigned in the CONSTRAINTS input card)
  !
  PUBLIC :: nconstr,     &
            constr_tol,  &
            constr_type, &
            constr,      &
            lagrange,    &
            target,      &
            dmax
  !
  ! ... global variables
  !
  INTEGER               :: nconstr 
  REAL(DP)              :: constr_tol
  INTEGER,  ALLOCATABLE :: constr_type(:)
  REAL(DP), ALLOCATABLE :: constr(:,:)
  REAL(DP), ALLOCATABLE :: target(:)
  REAL(DP), ALLOCATABLE :: lagrange(:)
  REAL(DP)              :: dmax
  !
  CONTAINS
     !
     ! ... public methods
     !
     !-----------------------------------------------------------------------
     SUBROUTINE init_constraint( nat, tau, ityp, tau_units )
       !-----------------------------------------------------------------------
       !
       ! ... this routine is used to initialize constraints variables and
       ! ... collective variables (notice that collective variables are
       ! ... implemented as normal constraints but are read using specific
       ! ... input variables)
       !
       USE input_parameters, ONLY : nconstr_inp, constr_tol_inp, &
                                    ncolvar_inp, colvar_tol_inp, &
                                    constr_type_inp, constr_inp, &
                                    colvar_type_inp, colvar_inp, &
                                    constr_target, constr_target_set, &
                                    colvar_target, colvar_target_set, nc_fields
       !
       IMPLICIT NONE
       !
       INTEGER,  INTENT(IN) :: nat
       REAL(DP), INTENT(IN) :: tau(3,nat)
       INTEGER,  INTENT(IN) :: ityp(nat)
       REAL(DP), INTENT(IN) :: tau_units
       !
       INTEGER     :: i, j
       INTEGER     :: n, ia, ia0, ia1, ia2, ia3, n_type_coord1
       REAL(DP)    :: d0(3), d1(3), d2(3)
       REAL(DP)    :: C00, C01, C02, C11, C12, C22
       REAL(DP)    :: D01, D12
       REAL(DP)    :: smoothing, r_c
       INTEGER     :: type_coord1, type_coord2
       REAL(DP)    :: dtau(3), norm_dtau
       REAL(DP)    :: k(3), phase, norm_k
       COMPLEX(DP) :: struc_fac
       !
       CHARACTER(LEN=6), EXTERNAL :: int_to_char
       !
       !
       nconstr    = ncolvar_inp + nconstr_inp
       constr_tol = MAX( constr_tol_inp, colvar_tol_inp )
       !
       ALLOCATE( lagrange(    nconstr ) )
       ALLOCATE( target(      nconstr ) )
       ALLOCATE( constr_type( nconstr ) )
       !
       ALLOCATE( constr( nc_fields, nconstr ) )
       !
       ! ... setting constr to 0 to findout which elements have
       !     been set to an atomic index. This is required for CP.
       !
       constr = 0.0d0
       !
       ! ... NB: the first "ncolvar" constraints are collective variables (used
       ! ...     for meta-dynamics and free-energy smd), the remaining are real
       ! ...     constraints
       !
       constr(:,1:ncolvar_inp)         = colvar_inp(:,1:ncolvar_inp)
       constr(:,ncolvar_inp+1:nconstr) = constr_inp(:,1:nconstr_inp)
       !
       ! ... set the largest possible distance among two atoms within
       ! ... the supercell
       !
       IF ( ncolvar_inp > 0 ) THEN
          !
          IF ( ANY( colvar_type_inp(:) == 'distance' ) ) CALL compute_dmax()
          !
       ELSE IF ( nconstr_inp > 0 ) THEN
          !
          IF ( ANY( constr_type_inp(:) == 'distance' ) ) CALL compute_dmax()
          !
       END IF
       !
       ! ... initializations of target values for the constraints :
       !
       ! ... first the initialization of the collective variables
       !
       DO ia = 1, ncolvar_inp
          !
          SELECT CASE ( colvar_type_inp(ia) )
          CASE( 'type_coord' )
             !
             ! ... constraint on global coordination-number, i.e. the average 
             ! ... number of atoms of type B surrounding the atoms of type A
             !
             constr_type(ia) = 1
             !
             IF ( colvar_target_set(ia) ) THEN
                !
                target(ia) = colvar_target(ia)
                !
                CYCLE
                !
             ELSE
                !
                CALL set_type_coord( ia )
                !
             END IF             
             !
          CASE( 'atom_coord' )
             !
             ! ... constraint on local coordination-number, i.e. the average 
             ! ... number of atoms of type A surrounding a specific atom
             !
             constr_type(ia) = 2
             !
             IF ( colvar_target_set(ia) ) THEN
                !
                target(ia) = colvar_target(ia)
                !
                CYCLE
                !
             ELSE
                !
                CALL set_atom_coord( ia )
                !
             END IF
             !
          CASE( 'distance' )
             !
             constr_type(ia) = 3
             !
             IF ( colvar_target_set(ia) ) THEN
                !
                target(ia) = colvar_target(ia)
                !
             ELSE
                !
                ia1 = ANINT( constr(1,ia) )
                ia2 = ANINT( constr(2,ia) )
                !
                dtau(:) = pbc( ( tau(:,ia1) - tau(:,ia2) ) * tau_units )
                !
                target(ia) = norm( dtau(:) )
                !
             END IF
             !
             IF ( target(ia) > dmax )  THEN
                !
                WRITE( stdout, '(/,5X,"target = ",F12.8,/, &
                                &  5X,"dmax   = ",F12.8)' ) target(ia), dmax
                !
                CALL errore( 'init_constraint', 'the target for coll.var. '  //&
                           & TRIM( int_to_char( ia ) ) // ' is larger than ' //&
                           & 'the largest possible value', 1 )
                !
             END IF
             !
          CASE( 'planar_angle' )
             !
             ! ... constraint on planar angle (for the notation used here see 
             ! ... Appendix C of the Allen-Tildesley book)
             !
             constr_type(ia) = 4
             !
             IF ( colvar_target_set(ia) ) THEN
                !
                ! ... in the input target for the angle (in degrees) is 
                ! ... converted to the cosine of the angle
                !
                target(ia) = COS( ( 180.D0 - colvar_target(ia) )* tpi / 360.D0 )
                !
                CYCLE
                !
             ELSE
                !
                CALL set_planar_angle( ia )
                !
             END IF
             !
          CASE( 'torsional_angle' )
             !
             ! ... constraint on torsional angle (for the notation used here 
             ! ... see Appendix C of the Allen-Tildesley book)
             !
             constr_type(ia) = 5
             !
             IF ( colvar_target_set(ia) ) THEN
                !
                ! ... in the input target for the torsional angle (in degrees)
                ! ... is converted to the cosine of the angle
                !
                target(ia) = COS( colvar_target(ia) * tpi / 360.D0 )
                !
                CYCLE
                !
             ELSE
                !
                CALL set_torsional_angle( ia )
                !
             END IF
             !
          CASE( 'struct_fac' )
             !
             ! ... constraint on structure factor at a given k-vector
             !
             constr_type(ia) = 6
             !
             IF ( colvar_target_set(ia) ) THEN
                !
                target(ia) = colvar_target(ia)
                !
                CYCLE
                !
             ELSE
                !
                CALL set_structure_factor( ia )
                !
             END IF
             !
          CASE( 'sph_struct_fac' )
             !
             ! ... constraint on spherical average of the structure factor for
             ! ... a given k-vector of norm k
             !
             constr_type(ia) = 7
             !
             IF ( colvar_target_set(ia) ) THEN
                !
                target(ia) = colvar_target(ia)
                !
                CYCLE
                !
             ELSE
                !
                CALL set_sph_structure_factor( ia )
                !
             END IF
             !
          CASE( 'bennett_proj' )
             !
             ! ... constraint on the projection onto a given direction of the 
             ! ... vector defined by the position of one atom minus the center
             ! ... of mass of the others
             ! ... ( Ch.H. Bennett in Diffusion in Solids, Recent Developments,
             ! ...   Ed. by A.S. Nowick and J.J. Burton, New York 1975 )
             !
             constr_type(ia) = 8
             !
             IF ( colvar_target_set(ia) ) THEN
                !
                target(ia) = colvar_target(ia)
                !
                CYCLE
                !
             ELSE
                !
                CALL set_bennett_proj( ia )
                !
             END IF
             !
          CASE DEFAULT
             !
             CALL errore( 'init_constraint', &
                          'collective-variable type not implemented', 1 )
             !
          END SELECT
          !
       END DO       
       !
       ! ... then then the initialization of the real constraints
       !
       DO n = 1, nconstr_inp
          !
          ia = ncolvar_inp + n
          !
          SELECT CASE ( constr_type_inp(n) )
          CASE( 'type_coord' )
             !
             ! ... constraint on global coordination-number, i.e. the average 
             ! ... number of atoms of type B surrounding the atoms of type A
             !
             constr_type(ia) = 1
             !
             IF ( constr_target_set(n) ) THEN
                !
                target(ia) = constr_target(n)
                !
                CYCLE
                !
             ELSE
                !
                CALL set_type_coord( ia )
                !
             END IF             
             !
          CASE( 'atom_coord' )
             !
             ! ... constraint on local coordination-number, i.e. the average 
             ! ... number of atoms of type A surrounding a specific atom
             !
             constr_type(ia) = 2
             !
             IF ( constr_target_set(n) ) THEN
                !
                target(ia) = constr_target(n)
                !
                CYCLE
                !
             ELSE
                !
                CALL set_atom_coord( ia )
                !
             END IF
             !
          CASE( 'distance' )
             !
             constr_type(ia) = 3
             !
             IF ( constr_target_set(n) ) THEN
                !
                target(ia) = constr_target(n)
                !
             ELSE
                !
                ia1 = ANINT( constr(1,ia) )
                ia2 = ANINT( constr(2,ia) )
                !
                dtau(:) = pbc( ( tau(:,ia1) - tau(:,ia2) ) * tau_units )
                !
                target(ia) = norm( dtau(:) )
                !
             END IF
             !
             IF ( target(ia) > dmax )  THEN
                !
                WRITE( stdout, '(/,5X,"target = ",F12.8,/, &
                                &  5X,"dmax   = ",F12.8)' ) target(ia), dmax
                !
                CALL errore( 'init_constraint', 'the target for constraint ' //&
                           & TRIM( int_to_char( n ) ) // ' is larger than '  //&
                           & 'the largest possible value', 1 )
                !
             END IF
             !
          CASE( 'planar_angle' )
             !
             ! ... constraint on planar angle (for the notation used here see 
             ! ... Appendix C of the Allen-Tildesley book)
             !
             constr_type(ia) = 4
             !
             IF ( constr_target_set(n) ) THEN
                !
                ! ... in the input target for the angle (in degrees) is 
                ! ... converted to the cosine of the angle
                !
                target(ia) = COS( ( 180.D0 - constr_target(n) )* tpi / 360.D0 )
                !
                CYCLE
                !
             ELSE
                !
                CALL set_planar_angle( ia )
                !
             END IF
             !
          CASE( 'torsional_angle' )
             !
             ! ... constraint on torsional angle (for the notation used here 
             ! ... see Appendix C of the Allen-Tildesley book)
             !
             constr_type(ia) = 5
             !
             IF ( constr_target_set(n) ) THEN
                !
                ! ... in the input target for the torsional angle (in degrees)
                ! ... is converted to the cosine of the angle
                !
                target(ia) = COS( constr_target(n) * tpi / 360.D0 )
                !
                CYCLE
                !
             ELSE
                !
                CALL set_torsional_angle( ia )
                !
             END IF
             !
          CASE( 'struct_fac' )
             !
             ! ... constraint on structure factor at a given k-vector
             !
             constr_type(ia) = 6
             !
             IF ( constr_target_set(n) ) THEN
                !
                target(ia) = constr_target(n)
                !
                CYCLE
                !
             ELSE
                !
                CALL set_structure_factor( ia )
                !
             END IF
             !
          CASE( 'sph_struct_fac' )
             !
             ! ... constraint on spherical average of the structure factor for
             ! ... a given k-vector of norm k
             !
             constr_type(ia) = 7
             !
             IF ( constr_target_set(n) ) THEN
                !
                target(ia) = constr_target(n)
                !
                CYCLE
                !
             ELSE
                !
                CALL set_sph_structure_factor( ia )
                !
             END IF
             !
          CASE( 'bennett_proj' )
             !
             ! ... constraint on the projection onto a given direction of the 
             ! ... vector defined by the position of one atom minus the center
             ! ... of mass of the others
             ! ... ( Ch.H. Bennett in Diffusion in Solids, Recent Developments,
             ! ...   Ed. by A.S. Nowick and J.J. Burton, New York 1975 )
             !
             constr_type(ia) = 8
             !
             IF ( constr_target_set(n) ) THEN
                !
                target(ia) = constr_target(n)
                !
                CYCLE
                !
             ELSE
                !
                CALL set_bennett_proj( ia )
                !
             END IF
             !
          CASE DEFAULT
             !
             CALL errore( 'init_constraint', &
                          'constraint type not implemented', 1 )
             !
          END SELECT
          !
       END DO
       !
       RETURN
       !
       CONTAINS
         !
         !-------------------------------------------------------------------
         SUBROUTINE set_type_coord( ia )
           !-------------------------------------------------------------------
           !
           INTEGER, INTENT(IN) :: ia
           !
           type_coord1 = ANINT( constr(1,ia) )
           type_coord2 = ANINT( constr(2,ia) )
           !
           r_c  = constr(3,ia)
           !
           smoothing = 1.D0 / constr(4,ia)
           !
           target(ia) = 0.D0
           !
           n_type_coord1 = 0
           !
           DO ia1 = 1, nat
              !
              IF ( ityp(ia1) /= type_coord1 ) CYCLE
              !
              DO ia2 = 1, nat
                 !
                 IF ( ia2 == ia1 ) CYCLE
                 !
                 IF ( ityp(ia2) /= type_coord2 ) CYCLE
                 !
                 dtau(:) = pbc( ( tau(:,ia1) - tau(:,ia2) ) * tau_units )
                 !
                 norm_dtau = norm( dtau(:) )
                 !
                 target(ia) = target(ia) + 1.D0 / &
                              ( EXP( smoothing * ( norm_dtau - r_c ) ) + 1.D0 )
                 !
              END DO
              !
              n_type_coord1 = n_type_coord1 + 1
              !
           END DO
           !
           target(ia) = target(ia) / DBLE( n_type_coord1 )
           !
         END SUBROUTINE set_type_coord
         !
         !-------------------------------------------------------------------
         SUBROUTINE set_atom_coord( ia )
           !-------------------------------------------------------------------
           !
           INTEGER, INTENT(IN) :: ia
           !
           ia1         = ANINT( constr(1,ia) )
           type_coord1 = ANINT( constr(2,ia) )
           !
           r_c = constr(3,ia)
           !
           smoothing = 1.D0 / constr(4,ia)
           !
           target(ia) = 0.D0
           !
           DO ia2 = 1, nat
              !
              IF ( ia2 == ia1 ) CYCLE
              !
              IF ( ityp(ia2) /= type_coord1 ) CYCLE
              !
              dtau(:) = pbc( ( tau(:,ia1) - tau(:,ia2) ) * tau_units )
              !
              norm_dtau = norm( dtau(:) )
              !
              target(ia) = target(ia) + 1.D0 / &
                           ( EXP( smoothing * ( norm_dtau - r_c ) ) + 1.D0 )
              !
           END DO
           !
         END SUBROUTINE set_atom_coord
         !
         !-------------------------------------------------------------------
         SUBROUTINE set_planar_angle( ia )
           !-------------------------------------------------------------------
           !
           INTEGER, INTENT(IN) :: ia
           !
           ia0 = ANINT( constr(1,ia) )
           ia1 = ANINT( constr(2,ia) )
           ia2 = ANINT( constr(3,ia) )
           !
           d0(:) = pbc( ( tau(:,ia0) - tau(:,ia1) ) * tau_units )
           d1(:) = pbc( ( tau(:,ia1) - tau(:,ia2) ) * tau_units )
           !
           d0(:) = d0(:) / norm( d0(:) )
           d1(:) = d1(:) / norm( d1(:) )
           !
           target(ia) = d0(:) .dot. d1(:)
           !
         END SUBROUTINE set_planar_angle
         !
         !-------------------------------------------------------------------
         SUBROUTINE set_torsional_angle( ia )
           !-------------------------------------------------------------------
           !
           INTEGER, INTENT(IN) :: ia
           !
           ia0 = ANINT( constr(1,ia) )
           ia1 = ANINT( constr(2,ia) )
           ia2 = ANINT( constr(3,ia) )
           ia3 = ANINT( constr(4,ia) )
           !
           d0(:) = pbc( ( tau(:,ia0) - tau(:,ia1) ) * tau_units )
           d1(:) = pbc( ( tau(:,ia1) - tau(:,ia2) ) * tau_units )
           d2(:) = pbc( ( tau(:,ia2) - tau(:,ia3) ) * tau_units )
           !
           C00 = d0(:) .dot. d0(:)
           C01 = d0(:) .dot. d1(:)
           C11 = d1(:) .dot. d1(:)
           C02 = d0(:) .dot. d2(:)
           C12 = d1(:) .dot. d2(:)
           C22 = d2(:) .dot. d2(:)
           !
           D01 = C00 * C11 - C01 * C01
           D12 = C11 * C22 - C12 * C12
           !
           target(ia) = ( C01 * C12 - C02 * C11 ) / SQRT( D01 * D12 )
           !
         END SUBROUTINE set_torsional_angle
         !
         !-------------------------------------------------------------------
         SUBROUTINE set_structure_factor( ia )
           !-------------------------------------------------------------------
           !
           INTEGER, INTENT(IN) :: ia
           !
           k(1) = constr(1,ia) * tpi / tau_units
           k(2) = constr(2,ia) * tpi / tau_units
           k(3) = constr(3,ia) * tpi / tau_units
           !
           struc_fac = ( 0.D0, 0.D0 )
           !
           DO i = 1, nat
              !
              dtau(:) = pbc( ( tau(:,i) - tau(:,1) ) * tau_units )
              !
              phase = k(:) .dot. dtau(:)
              !
              struc_fac = struc_fac + CMPLX( COS( phase ), SIN( phase ) )
              !
           END DO
           !
           target(ia) = ( CONJG( struc_fac ) * struc_fac ) / DBLE( nat*nat )
           !
         END SUBROUTINE set_structure_factor
         !
         !-------------------------------------------------------------------
         SUBROUTINE set_sph_structure_factor( ia )
           !-------------------------------------------------------------------
           !
           INTEGER, INTENT(IN) :: ia
           !
           norm_k = constr(1,ia) * tpi / tau_units
           !
           target(ia) = 0.D0
           !
           DO i = 1, nat - 1
              !
              DO j = i + 1, nat
                 !
                 dtau(:) = pbc( ( tau(:,i) - tau(:,j) ) * tau_units )
                 !
                 norm_dtau = norm( dtau(:) )
                 !
                 phase = norm_k * norm_dtau
                 !
                 IF ( phase < eps32 ) THEN
                    !
                    target(ia) = target(ia) + 1.D0
                    !
                 ELSE
                    !
                    target(ia) = target(ia) + SIN( phase ) / phase
                    !
                 END IF
                 !
              END DO
              !
           END DO
           !
           target(ia) = 2.D0 * fpi * target(ia) / DBLE( nat )
           !
         END SUBROUTINE set_sph_structure_factor
         !
         !-------------------------------------------------------------------
         SUBROUTINE set_bennett_proj( ia )
           !-------------------------------------------------------------------
           !
           INTEGER, INTENT(IN) :: ia
           !
           ia0 = ANINT( constr(1,ia) )
           !
           d0(:) = tau(:,ia0)
           d1(:) = SUM( tau(:,:), DIM = 2 )
           !
           d1(:) = pbc( ( d1(:) - d0(:) )*tau_units ) / DBLE( nat - 1 ) - &
                   pbc( d0(:)*tau_units )
           !
           d2(:) = constr(2:4,ia)
           !
           target(ia) = ( d1(:) .dot. d2(:) ) / tau_units
           !
         END SUBROUTINE set_bennett_proj
         !
     END SUBROUTINE init_constraint
     !
     !-----------------------------------------------------------------------
     SUBROUTINE constraint_grad( idx, nat, tau, &
                                 if_pos, ityp, tau_units, g, dg )
       !-----------------------------------------------------------------------
       !
       ! ... this routine computes the value of the constraint equation and 
       ! ... the corresponding constraint gradient 
       !
       IMPLICIT NONE
       !
       INTEGER,  INTENT(IN)  :: idx
       INTEGER,  INTENT(IN)  :: nat
       REAL(DP), INTENT(IN)  :: tau(:,:)
       INTEGER,  INTENT(IN)  :: if_pos(:,:)
       INTEGER,  INTENT(IN)  :: ityp(:)
       REAL(DP), INTENT(IN)  :: tau_units
       REAL(DP), INTENT(OUT) :: dg(:,:)
       REAL(DP), INTENT(OUT) :: g
       !
       INTEGER     :: i, j
       INTEGER     :: ia, ia0, ia1, ia2, ia3, n_type_coord1
       REAL(DP)    :: d0(3), d1(3), d2(3)
       REAL(DP)    :: inv_den, fac
       REAL(DP)    :: C00, C01, C02, C11, C12, C22
       REAL(DP)    :: D01, D12, invD01, invD12
       REAL(DP)    :: smoothing, r_c
       INTEGER     :: type_coord1, type_coord2
       REAL(DP)    :: dtau(3), norm_dtau, norm_dtau_sq, expo
       REAL(DP)    :: r0(3), ri(3), k(3), phase, ksin(3), norm_k, sinxx
       COMPLEX(DP) :: struc_fac
       !
       REAL(DP), EXTERNAL :: DDOT
       !
       !
       dg(:,:) = 0.D0
       !
       SELECT CASE ( constr_type(idx) )
       CASE( 1 )
          !
          ! ... constraint on global coordination
          !
          type_coord1 = ANINT( constr(1,idx) )
          type_coord2 = ANINT( constr(2,idx) )
          !
          r_c = constr(3,idx)
          !
          smoothing = 1.D0 / constr(4,idx)
          !
          g = 0.D0
          !
          n_type_coord1 = 0
          !
          DO ia1 = 1, nat
             !
             IF ( ityp(ia1) /= type_coord1 ) CYCLE
             !
             DO ia2 = 1, nat
                !
                IF ( ia2 == ia1 ) CYCLE
                !
                IF ( ityp(ia2) /= type_coord2 ) CYCLE
                !
                dtau(:) = pbc( ( tau(:,ia1) - tau(:,ia2) ) * tau_units )
                !
                norm_dtau = norm( dtau(:) )
                !
                dtau(:) = dtau(:) / norm_dtau
                !
                expo = EXP( smoothing * ( norm_dtau - r_c ) )
                !
                g = g + 1.D0 / ( expo + 1.D0 )
                !
                dtau(:) = dtau(:) * smoothing * expo / ( expo + 1.D0 )**2
                !
                dg(:,ia2) = dg(:,ia2) + dtau(:)
                dg(:,ia1) = dg(:,ia1) - dtau(:)
                !
             END DO
             !
             n_type_coord1 = n_type_coord1 + 1
             !
          END DO
          !
          g  = g  / DBLE( n_type_coord1 )
          dg = dg / DBLE( n_type_coord1 )
          !
          g = ( g - target(idx) )
          !
       CASE( 2 )
          !
          ! ... constraint on local coordination
          !
          ia          = ANINT( constr(1,idx) )
          type_coord1 = ANINT( constr(2,idx) )
          !
          r_c = constr(3,idx)
          !
          smoothing = 1.D0 / constr(4,idx)
          !
          g = 0.D0
          !
          DO ia1 = 1, nat
             !
             IF ( ia1 == ia ) CYCLE
             !
             IF ( ityp(ia1) /= type_coord1 ) CYCLE
             !
             dtau(:) = pbc( ( tau(:,ia) - tau(:,ia1) ) * tau_units )
             !
             norm_dtau = norm( dtau(:) )
             !
             dtau(:) = dtau(:) / norm_dtau
             !
             expo = EXP( smoothing * ( norm_dtau - r_c ) )
             !
             g = g + 1.D0 / ( expo + 1.D0 )
             !
             dtau(:) = dtau(:) * smoothing * expo / ( expo + 1.D0 )**2
             !
             dg(:,ia1) = dg(:,ia1) + dtau(:)
             dg(:,ia)  = dg(:,ia)  - dtau(:)
             !
          END DO
          !
          g = ( g - target(idx) )
          !
       CASE( 3 )
          !
          ! ... constraint on distances
          !
          ia1 = ANINT( constr(1,idx) )
          ia2 = ANINT( constr(2,idx) )
          !
          dtau(:) = pbc( ( tau(:,ia1) - tau(:,ia2) ) * tau_units )
          !
          norm_dtau = norm( dtau(:) )
          !
          g = ( norm_dtau - target(idx) )
          !
          dg(:,ia1) = dtau(:) / norm_dtau
          !
          dg(:,ia2) = - dg(:,ia1)
          !
       CASE( 4 )
          !
          ! ... constraint on planar angles (for the notation used here see
          ! ... Appendix C of the Allen-Tildesley book)
          !
          ia0 = ANINT( constr(1,idx) )
          ia1 = ANINT( constr(2,idx) )
          ia2 = ANINT( constr(3,idx) )
          !
          d0(:) = pbc( ( tau(:,ia0) - tau(:,ia1) ) * tau_units )
          d1(:) = pbc( ( tau(:,ia1) - tau(:,ia2) ) * tau_units )
          !
          C00 = d0(:) .dot. d0(:)
          C01 = d0(:) .dot. d1(:)
          C11 = d1(:) .dot. d1(:)
          !
          inv_den = 1.D0 / SQRT( C00 * C11 )
          !
          g = ( C01 * inv_den - target(idx) )
          !
          dg(:,ia0) = ( d1(:) - C01 / C00 * d0(:) ) * inv_den
          dg(:,ia2) = ( C01 / C11 * d1(:) - d0(:) ) * inv_den
          dg(:,ia1) = - dg(:,ia0) - dg(:,ia2)
          !
       CASE( 5 )
          !
          ! ... constraint on torsional angle (for the notation used here 
          ! ... see Appendix C of the Allen-Tildesley book)
          !
          ia0 = ANINT( constr(1,idx) )
          ia1 = ANINT( constr(2,idx) )
          ia2 = ANINT( constr(3,idx) )
          ia3 = ANINT( constr(4,idx) )
          !
          d0(:) = pbc( ( tau(:,ia0) - tau(:,ia1) ) * tau_units )
          d1(:) = pbc( ( tau(:,ia1) - tau(:,ia2) ) * tau_units )
          d2(:) = pbc( ( tau(:,ia2) - tau(:,ia3) ) * tau_units )
          !
          C00 = d0(:) .dot. d0(:)
          C01 = d0(:) .dot. d1(:)
          C11 = d1(:) .dot. d1(:)
          C02 = d0(:) .dot. d2(:)
          C12 = d1(:) .dot. d2(:)
          C22 = d2(:) .dot. d2(:)
          !
          D01 = C00 * C11 - C01 * C01
          D12 = C11 * C22 - C12 * C12
          !
          IF ( ABS( D01 ) < eps32 .OR. ABS( D12 ) < eps32 ) &
             CALL errore( 'constraint_grad', 'either D01 or D12 is zero', 1 )
          !
          invD01 = 1.D0 / D01
          invD12 = 1.D0 / D12
          !
          fac = C01 * C12 - C02 * C11
          !
          inv_den = 1.D0 / SQRT( D01 * D12 )
          !
          g = ( ( C01 * C12 - C02 * C11 ) * inv_den - target(idx) )
          !
          dg(:,ia0) = ( C12 * d1(:) - C11 * d2(:) - &
                        invD01 * fac * ( C11 * d0(:) - C01 * d1(:) ) ) * inv_den
          !
          dg(:,ia2) = ( C01 * ( d1(:) - d2(:) ) - &
                        ( C11 + C12 ) * d0(:) + 2.D0 * C02 * d1(:) - &
                        invD12 * fac * ( ( C11 + C12 ) * d2(:) - &
                                         ( C12 + C22 ) * d1(:) ) - &
                        invD01 * fac * ( C01 * d0(:) - C00 * d1(:) ) ) * inv_den
          !
          dg(:,ia3) = ( C11 * d0(:) - C01 * d1(:) - &
                        invD12 * fac * ( C12 * d1(:) - C11 * d2(:) ) ) * inv_den
          !
          dg(:,ia1) = - dg(:,ia0) - dg(:,ia2) - dg(:,ia3)
          !
       CASE( 6 )
          !
          ! ... constraint on structure factor at a given k vector
          !
          k(1) = constr(1,idx) * tpi / tau_units
          k(2) = constr(2,idx) * tpi / tau_units
          k(3) = constr(3,idx) * tpi / tau_units
          !
          struc_fac = ( 1.D0, 0.D0 )
          !
          r0(:) = tau(:,1)
          !
          DO i = 1, nat - 1
             !
             dtau(:) = pbc( ( tau(:,i+1) - r0(:) ) * tau_units )
             !
             phase = k(1)*dtau(1) + k(2)*dtau(2) + k(3)*dtau(3)
             !
             struc_fac = struc_fac + CMPLX( COS( phase ), SIN( phase ) )
             !
             ri(:) = tau(:,i)
             !
             DO j = i + 1, nat
                !
                dtau(:) = pbc( ( tau(:,j) - ri(:) ) * tau_units )
                !
                phase = k(1)*dtau(1) + k(2)*dtau(2) + k(3)*dtau(3)
                !
                ksin(:) = k(:) * SIN( phase )
                !
                dg(:,i) = dg(:,i) + ksin(:)
                dg(:,j) = dg(:,j) - ksin(:)
                !
             END DO
             !
          END DO
          !
          g = ( CONJG( struc_fac ) * struc_fac ) / DBLE( nat*nat )
          !
          g = ( g - target(idx) )
          !
          dg(:,:) = dg(:,:) * 2.D0 / DBLE( nat*nat )
          !
       CASE( 7 )
          !
          ! ... constraint on spherical average of the structure factor for
          ! ... a given k-vector of norm k
          !
          norm_k = constr(1,idx) * tpi / tau_units
          !
          g = 0.D0
          !
          DO i = 1, nat - 1
             !
             ri(:) = tau(:,i)
             !
             DO j = i + 1, nat
                !
                dtau(:) = pbc( ( ri(:) - tau(:,j) ) * tau_units )
                !
                norm_dtau_sq = dtau(1)**2 + dtau(2)**2 + dtau(3)**2
                !
                norm_dtau = SQRT( norm_dtau_sq )
                !
                phase = norm_k * norm_dtau
                !
                IF ( phase < eps32 ) THEN
                   !
                   g = g + 1.D0
                   !
                ELSE
                   !
                   sinxx = SIN( phase ) / phase
                   !
                   g = g + sinxx
                   !
                   dtau(:) = dtau(:) / norm_dtau_sq * ( COS( phase ) - sinxx )
                   !
                   dg(:,i) = dg(:,i) + dtau(:)
                   dg(:,j) = dg(:,j) - dtau(:)
                   !
                END IF
                !
             END DO
             !
          END DO
          !
          g = ( 2.D0 * fpi * g / DBLE( nat ) - target(idx) )
          !
          dg(:,:) = 4.D0 * fpi * dg(:,:) / DBLE( nat )
          !
       CASE( 8 )
          !
          ! ... constraint on Bennett projection
          !
          ia0 = ANINT( constr(1,idx) )
          !
          d0(:) = tau(:,ia0)
          d1(:) = SUM( tau(:,:), DIM = 2 )
          !
          d1(:) = pbc( ( d1(:) - d0(:) )*tau_units ) / DBLE( nat - 1 ) - &
                  pbc( d0(:)*tau_units )
          !
          d2(:) = constr(2:4,idx)
          !
          g = ( d1(:) .dot. d2(:) ) / tau_units - target( idx )
          !
          dg = 0.D0
          !
          C00 = ( 1.D0 / DBLE( nat - 1 ) ) / tau_units
          C01 = -1.D0 / tau_units
          !
          DO i = 1, nat
             !
             dg(:,i) = d2(:)*C00
             !
          END DO
          !
          dg(:,ia0) = d2(:)*C01
          !
       END SELECT
       !
       dg(:,:) = dg(:,:) * DBLE( if_pos(:,:) )
       !
       RETURN
       !
     END SUBROUTINE constraint_grad
     !
     !-----------------------------------------------------------------------
     SUBROUTINE check_constraint( nat, taup, tau0, &
                                  force, if_pos, ityp, tau_units, dt, massconv )
       !-----------------------------------------------------------------------
       !
       ! ... update taup (predicted positions) so that the constraint equation 
       ! ... g=0 is satisfied, using the recursion formula:
       !
       ! ...                       g(taup)
       ! ... taup = taup - ----------------------- * dg(tau0)
       ! ...               M^-1<dg(taup)|dg(tau0)>
       !
       ! ... in normal cases the constraint equation should be satisfied at 
       ! ... the very first iteration.
       !
       USE ions_base, ONLY : amass
       !
       IMPLICIT NONE
       !
       INTEGER,  INTENT(IN)    :: nat
       REAL(DP), INTENT(INOUT) :: taup(3,nat)
       REAL(DP), INTENT(IN)    :: tau0(3,nat)
       INTEGER,  INTENT(IN)    :: if_pos(3,nat)
       REAL(DP), INTENT(INOUT) :: force(3,nat)
       INTEGER,  INTENT(IN)    :: ityp(nat)
       REAL(DP), INTENT(IN)    :: tau_units
       REAL(DP), INTENT(IN)    :: dt
       REAL(DP), INTENT(IN)    :: massconv
       !
       INTEGER               :: na, i, idx, dim
       REAL(DP), ALLOCATABLE :: gp(:), dgp(:,:), dg0(:,:,:)
       REAL(DP)              :: g0
       REAL(DP)              :: lambda, fac, invdtsq
       LOGICAL, ALLOCATABLE  :: ltest(:)
       LOGICAL               :: global_test
       INTEGER, PARAMETER    :: maxiter = 100
       !
       REAL(DP), EXTERNAL :: DDOT
       !
       !
       ALLOCATE( dgp( 3, nat ) )
       ALLOCATE( dg0( 3, nat, nconstr ) )
       !
       ALLOCATE( gp(    nconstr ) )
       ALLOCATE( ltest( nconstr ) )
       !
       invdtsq  = 1.D0 / dt**2
       !
       dim = 3*nat
       !
       DO idx = 1, nconstr
          !
          CALL constraint_grad( idx, nat, tau0, &
                                if_pos, ityp, tau_units, g0, dg0(:,:,idx) )
          !
       END DO
       !
       outer_loop: DO i = 1, maxiter
          !
          inner_loop: DO idx = 1, nconstr
             !
             ltest(idx) = .FALSE.
             !
             CALL constraint_grad( idx, nat, taup, &
                                   if_pos, ityp, tau_units, gp(idx), dgp )
             !
             ! ... check if gp = 0
             !
#if defined (__DEBUG_CONSTRAINTS)
             WRITE( stdout, '(2(2X,I3),F12.8)' ) i, idx, ABS( gp(idx) )
#endif
             !
             IF ( ABS( gp(idx) ) < constr_tol ) THEN
                !
                ltest(idx) = .TRUE.
                !
                CYCLE inner_loop
                !
             END IF
             !
             ! ... if  gp <> 0  find new taup and check again 
             ! ... ( gp is in bohr and taup in tau_units )
             !
             DO na = 1, nat
                !
                dgp(:,na) = dgp(:,na) / ( amass(ityp(na)) * massconv )
                !
             END DO
             !
             lambda = gp(idx) / DDOT( dim, dgp, 1, dg0(:,:,idx), 1 )
             !
             DO na = 1, nat
                !
                fac = amass(ityp(na)) * massconv * tau_units
                !
                taup(:,na) = taup(:,na) - lambda * dg0(:,na,idx) / fac
                !
             END DO
             !
             lagrange(idx) = lagrange(idx) + lambda * invdtsq
             !
             force(:,:) = force(:,:) - lambda * dg0(:,:,idx) * invdtsq
	     !
          END DO inner_loop
          !
          global_test = ALL( ltest(:) )
          !
          ! ... all constraints are satisfied
          !
          IF ( global_test ) EXIT outer_loop
          !
       END DO outer_loop
       !
       IF ( .NOT. global_test ) THEN
          !
          ! ... error messages
          !
          WRITE( stdout, '(/,5X,"Number of step(s): ",I3)') MIN( i, maxiter )
          WRITE( stdout, '(/,5X,"target convergence: ")' )
          !
          DO i = 1, nconstr
             !
             WRITE( stdout, '(5X,"constr # ",I3,2X,L1,3(2X,F16.10))' ) &
                 i, ltest(i), ABS( gp(i) ), constr_tol, target(i)
             !
          END DO
          !
          CALL errore( 'check_constraint', &
                       'on some constraint g = 0 is not satisfied', 1 )
          !
       END IF
       !
       DEALLOCATE( dgp )
       DEALLOCATE( dg0 )
       DEALLOCATE( gp )
       DEALLOCATE( ltest )
       !
       RETURN
       !
     END SUBROUTINE check_constraint
     !
     !-----------------------------------------------------------------------
     SUBROUTINE remove_constr_force( nat, tau, &
                                     if_pos, ityp, tau_units, force )
       !-----------------------------------------------------------------------
       !
       ! ... the component of the force that is orthogonal to the
       ! ... ipersurface defined by the constraint equations is removed
       ! ... and the corresponding value of the lagrange multiplier computed
       !
       IMPLICIT NONE
       !
       INTEGER,  INTENT(IN)    :: nat
       REAL(DP), INTENT(IN)    :: tau(:,:)
       INTEGER,  INTENT(IN)    :: if_pos(:,:)
       INTEGER,  INTENT(IN)    :: ityp(:)
       REAL(DP), INTENT(IN)    :: tau_units
       REAL(DP), INTENT(INOUT) :: force(:,:)
       !
       INTEGER               :: i, j, dim
       REAL(DP)              :: g, ndg, dgidgj
       REAL(DP)              :: norm_before, norm_after
       REAL(DP), ALLOCATABLE :: dg(:,:,:)
       REAL(DP), ALLOCATABLE :: dg_matrix(:,:)
       INTEGER,  ALLOCATABLE :: iwork(:)
       !
       REAL(DP), EXTERNAL :: DDOT, DNRM2
       !
       !
       dim = 3 * nat
       !
       lagrange(:) = 0.D0
       !
#if defined (__REMOVE_CONSTRAINT_FORCE)
       !
       norm_before = DNRM2( 3*nat, force, 1 )
       !
       ALLOCATE( dg( 3, nat, nconstr ) )
       !
       IF ( nconstr == 1 ) THEN
          !
          CALL constraint_grad( 1, nat, tau, &
                                if_pos, ityp, tau_units, g, dg(:,:,1) )
          !
          lagrange(1) = DDOT( dim, force, 1, dg(:,:,1), 1 )
          !
          ndg = DDOT( dim, dg(:,:,1), 1, dg(:,:,1), 1 )
          !
          force(:,:) = force(:,:) - lagrange(1) * dg(:,:,1) / ndg
          !
       ELSE
          !
          ALLOCATE( dg_matrix( nconstr, nconstr ) )
          ALLOCATE( iwork( nconstr ) )
          !
          DO i = 1, nconstr
             !
             CALL constraint_grad( i, nat, tau, &
                                   if_pos, ityp, tau_units, g, dg(:,:,i) )
             !
          END DO
          !
          DO i = 1, nconstr
             !
             dg_matrix(i,i) = DDOT( dim, dg(:,:,i), 1, dg(:,:,i), 1 )
             !
             lagrange(i) = DDOT( dim, force, 1, dg(:,:,i), 1 )
             !
             DO j = i + 1, nconstr
                !
                dgidgj = DDOT( dim, dg(:,:,i), 1, dg(:,:,j), 1 )
                !
                dg_matrix(i,j) = dgidgj
                dg_matrix(j,i) = dgidgj  
                !
             END DO
             !
          END DO
          !
          CALL DGESV( nconstr, 1, dg_matrix, &
                      nconstr, iwork, lagrange, nconstr, i )
          !
          IF ( i /= 0 ) &
             CALL errore( 'remove_constr_force', &
                          'error in the solution of the linear system', 1 )
          !
          DO i = 1, nconstr
             !
             force(:,:) = force(:,:) - lagrange(i) * dg(:,:,i)
             !
          END DO
          !
          DEALLOCATE( dg_matrix )
          DEALLOCATE( iwork )
          !
       END IF
       !
#if defined (__DEBUG_CONSTRAINTS)
       !
       WRITE( stdout, '(/,5X,"Intermediate forces (Ry/au):",/)')
       !
       DO i = 1, nat
          !
          WRITE( stdout, '(5X,"atom ",I3," type ",I2,3X,"force = ",3F14.8)' ) &
              i, ityp(i), force(:,i)
          !
       END DO
       !
#endif
       !
       norm_after = DNRM2( dim, force, 1 )
       !
       IF ( norm_before < norm_after ) THEN
          !
          WRITE( stdout, '(/,5X,"norm before = ",F16.10)' ) norm_before
          WRITE( stdout, '(  5X,"norm after  = ",F16.10)' ) norm_after
          !
          CALL errore( 'remove_constr_force', &
                       'norm(F) before < norm(F) after', 1 )
          !
       END IF
       !
       DEALLOCATE( dg )
       !
#endif
       !
     END SUBROUTINE remove_constr_force
     !
     !-----------------------------------------------------------------------
     SUBROUTINE deallocate_constraint()
       !-----------------------------------------------------------------------
       !
       IMPLICIT NONE
       !
       !
       IF ( ALLOCATED( lagrange ) )    DEALLOCATE( lagrange )
       IF ( ALLOCATED( constr ) )      DEALLOCATE( constr )
       IF ( ALLOCATED( constr_type ) ) DEALLOCATE( constr_type )
       IF ( ALLOCATED( target ) )      DEALLOCATE( target )       
       !
       RETURN
       !
     END SUBROUTINE deallocate_constraint
     !
     !-----------------------------------------------------------------------
     FUNCTION pbc( vect )
       !-----------------------------------------------------------------------
       !
       ! ... periodic boundary conditions ( vect is assumed to be given
       ! ... in cartesian coordinates and in atomic units )
       !
       USE cell_base, ONLY : at, bg, alat
       !
       IMPLICIT NONE
       !
       REAL(DP), INTENT(IN) :: vect(3)
       REAL(DP)             :: pbc(3)
       !
       !
#if defined (__USE_PBC)
       !
       pbc(:) = MATMUL( vect(:), bg(:,:) ) / alat
       !
       pbc(:) = pbc(:) - ANINT( pbc(:) )
       !
       pbc(:) = MATMUL( at(:,:), pbc(:) ) * alat
       !
#else
       !
       pbc(:) = vect(:)
       !
#endif
       RETURN
       !
     END FUNCTION pbc
     !
     !-----------------------------------------------------------------------
     SUBROUTINE compute_dmax()
       !-----------------------------------------------------------------------
       !
       ! ... dmax corresponds to one half the shortest edge of the cell
       !
       USE cell_base, ONLY : at, alat
       !
       IMPLICIT NONE
       !
       REAL(DP), PARAMETER :: x(3) = (/ 0.5D0, 0.0D0, 0.0D0 /), &
                              y(3) = (/ 0.0D0, 0.5D0, 0.0D0 /), &
                              z(3) = (/ 0.0D0, 0.0D0, 0.5D0 /)
       !
       dmax = norm( MATMUL( at(:,:), x(:) ) )       
       !
       dmax = MIN( dmax, norm( MATMUL( at(:,:), y(:) ) ) )
       !
       dmax = MIN( dmax, norm( MATMUL( at(:,:), z(:) ) ) )
       !
       dmax = dmax * alat
       !
       RETURN
       !
     END SUBROUTINE compute_dmax
     !
END MODULE constraints_module
