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
!#define __USE_LJ_PBC
!
!----------------------------------------------------------------------------
MODULE constraints_module
  !----------------------------------------------------------------------------
  ! 
  ! ... variables and methods for constraint Molecular Dynamics and
  ! ... constrained ionic relaxations (the SHAKE algorithm based on 
  ! ... lagrange multipliers) are defined here.
  !
  ! ... written by Carlo Sbraccia ( 24/02/2004 )
  !
  ! ... references :
  !
  ! ... 1) M. P. Allen and D. J. Tildesley, Computer Simulations of Liquids,
  ! ...    Clarendon Press - Oxford (1986)
  !
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
  PUBLIC :: init_constraint,     &
            check_constraint,    &
            remove_constr_force, &
            deallocate_constraint
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
  INTEGER, PARAMETER :: tab_dim = 1000
  REAL(DP)           :: sin_tab( 0 : tab_dim - 1 )
  !
  CONTAINS
     !
     ! ... public methods
     !
     !-----------------------------------------------------------------------
     SUBROUTINE init_constraint( nat, tau, tau_units, ityp )
       !-----------------------------------------------------------------------
       !
       USE input_parameters, ONLY : nconstr_inp, constr_tol_inp, &
                                    constr_type_inp, constr_inp, &
                                    constr_target, constr_target_set
       USE parser,           ONLY : int_to_char
       !
       IMPLICIT NONE
       !
       INTEGER,  INTENT(IN) :: nat
       REAL(DP), INTENT(IN) :: tau(3,nat)
       REAL(DP), INTENT(IN) :: tau_units
       INTEGER,  INTENT(IN) :: ityp(nat)
       !
       INTEGER     :: i, i_sin, i_cos
       INTEGER     :: ia, ia0, ia1, ia2, ia3, n_type_coord1
       REAL(DP)    :: d0(3), d1(3), d2(3)
       REAL(DP)    :: C00, C01, C02, C11, C12, C22
       REAL(DP)    :: D01, D12
       REAL(DP)    :: smearing, r_c, r_max
       INTEGER     :: type_coord1, type_coord2
       REAL(DP)    :: dtau(3), norm_dtau
       REAL(DP)    :: k(3), phase, norm_k
       COMPLEX(DP) :: struc_fac
       !
       !
       nconstr    = nconstr_inp
       constr_tol = constr_tol_inp
       !
       ALLOCATE( lagrange(    nconstr ) )
       ALLOCATE( target(      nconstr ) )
       ALLOCATE( constr_type( nconstr ) )
       !
       ALLOCATE( constr( SIZE( constr_inp(:,:), DIM = 1 ), nconstr ) )
       !
       constr_type(:) = constr_type_inp(1:nconstr)
       constr(:,:)    = constr_inp(:,1:nconstr)
       !
       ! ... set the largest possible distance among two atoms within
       ! ... the supercell
       !
       IF ( ANY( constr_type(:) == 3 ) ) CALL compute_dmax()
       !
       ! ... target value of the constrain ( in bohr )
       !
       DO ia = 1, nconstr
          !
          SELECT CASE ( constr_type(ia) )
          CASE( 1 )
             !
             ! ... constraint on global coordination-number, i.e. the average 
             ! ... number of atoms of type B surrounding the atoms of type A
             !
             IF ( constr_target_set(ia) ) THEN
                !
                target(ia) = constr_target(ia)
                !
                CYCLE
                !
             END IF
             !
             type_coord1 = ANINT( constr(1,ia) )
             type_coord2 = ANINT( constr(2,ia) )
             !
             r_c  = constr(3,ia)
             !
             smearing = 1.D0 / constr(4,ia)
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
                   dtau(:) = pbc( tau(:,ia1) - tau(:,ia2) ) * tau_units
                   !
                   norm_dtau = norm( dtau(:) )
                   !
                   target(ia) = target(ia) + 1.D0 / &
                                ( EXP( smearing * ( norm_dtau - r_c ) ) + 1.D0 )
                   !
                END DO
                !
                n_type_coord1 = n_type_coord1 + 1
                !
             END DO
             !
             target(ia) = target(ia) / DBLE( n_type_coord1 )
             !
          CASE( 2 )
             !
             ! ... constraint on local coordination-number, i.e. the average 
             ! ... number of atoms of type A surrounding a specific atom
             !
             IF ( constr_target_set(ia) ) THEN
                !
                target(ia) = constr_target(ia)
                !
                CYCLE
                !
             END IF
             !
             ia1         = ANINT( constr(1,ia) )
             type_coord1 = ANINT( constr(2,ia) )
             !
             r_c = constr(3,ia)
             !
             smearing = 1.D0 / constr(4,ia)
             !
             target(ia) = 0.D0
             !
             DO ia2 = 1, nat
                !
                IF ( ia2 == ia1 ) CYCLE
                !
                IF ( ityp(ia2) /= type_coord1 ) CYCLE
                !
                dtau(:) = pbc( tau(:,ia1) - tau(:,ia2) ) * tau_units
                !
                norm_dtau = norm( dtau(:) )
                !
                target(ia) = target(ia) + 1.D0 / &
                             ( EXP( smearing * ( norm_dtau - r_c ) ) + 1.D0 )
                !
             END DO
             !
          CASE( 3 )
             !
             ! ... constraint on distance
             !
             IF ( constr_target_set(ia) ) THEN
                !
                target(ia) = constr_target(ia)
                !
             ELSE
                !
                ia1 = ANINT( constr(1,ia) )
                ia2 = ANINT( constr(2,ia) )
                !
                dtau(:) = pbc( tau(:,ia1) - tau(:,ia2) ) * tau_units
                !
                target(ia) = norm( dtau(:) )
                !
             END IF
             !
             IF ( target(ia) > dmax ) &
                CALL errore( 'init_constraint', 'the target for constraint ' //&
                           & TRIM( int_to_char( ia ) ) // ' is larger than ' //&
                           & 'the largest possible value', 1 )
             !
          CASE( 4 )
             !
             ! ... constraint on planar angle (for the notation used here see 
             ! ... Appendix C of the Allen-Tildesley book)
             !
             IF ( constr_target_set(ia) ) THEN
                !
                ! ... in the input target for the angle (in degrees) is 
                ! ... converted to the cosine of the angle
                !
                target(ia) = COS( ( 180.D0 - constr_target(ia) )* tpi / 360.D0 )
                !
                CYCLE
                !
             END IF
             !
             ia0 = ANINT( constr(1,ia) )
             ia1 = ANINT( constr(2,ia) )
             ia2 = ANINT( constr(3,ia) )
             !
             d0(:) = pbc( tau(:,ia0) - tau(:,ia1) ) * tau_units
             d1(:) = pbc( tau(:,ia1) - tau(:,ia2) ) * tau_units
             !
             d0(:) = d0(:) / norm( d0(:) )
             d1(:) = d1(:) / norm( d1(:) )
             !
             target(ia) = d0(:) .dot. d1(:)
             !
          CASE( 5 )
             !
             ! ... constraint on torsional angle (for the notation used here 
             ! ... see Appendix C of the Allen-Tildesley book)
             !
             IF ( constr_target_set(ia) ) THEN
                !
                ! ... in the input target for the torsional angle (in degrees)
                ! ... is converted to the cosine of the angle
                !
                target(ia) = COS( constr_target(ia) * tpi / 360.D0 )
                !
                CYCLE
                !
             END IF
             !
             ia0 = ANINT( constr(1,ia) )
             ia1 = ANINT( constr(2,ia) )
             ia2 = ANINT( constr(3,ia) )
             ia3 = ANINT( constr(4,ia) )
             !
             d0(:) = pbc( tau(:,ia0) - tau(:,ia1) ) * tau_units
             d1(:) = pbc( tau(:,ia1) - tau(:,ia2) ) * tau_units
             d2(:) = pbc( tau(:,ia2) - tau(:,ia3) ) * tau_units
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
          CASE( 6 )
             !
             ! ... constraint on structure factor at a given k-vector
             !
             IF ( constr_target_set(ia) ) THEN
                !
                target(ia) = constr_target(ia)
                !
                CYCLE
                !
             END IF
             !
             k(1) = constr(1,ia) * tpi
             k(2) = constr(2,ia) * tpi
             k(3) = constr(3,ia) * tpi
             !
             struc_fac = ( 0.D0, 0.D0 )
             !
             DO i = 1, nat
                !
                dtau(:) = pbc( tau(:,i) - tau(:,1) )
                !
                phase = k(:) .dot. dtau(:)
                !
                struc_fac = struc_fac + CMPLX( COS( phase ), SIN( phase ) )
                !
             END DO
             !
             target(ia) = ( CONJG( struc_fac ) * struc_fac ) / DBLE( nat )**2
             !
          CASE( 7 )
             !
             ! ... constraint on spherical average of the structure factor for
             ! ... a given k-vector of norm k ( the 4*pi coefficient is
             ! ... understood )
             !
             IF ( constr_target_set(ia) ) THEN
                !
                target(ia) = constr_target(ia)
                !
                CYCLE
                !
             END IF
             !
             norm_k = constr(1,ia) * tpi / tau_units
             !
             target(ia) = 0.D0
             !
             DO i = 1, nat
                !
                dtau(:) = pbc( tau(:,i) - tau(:,1) ) * tau_units
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
             target(ia) = target(ia) / DBLE( nat )
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
     END SUBROUTINE init_constraint
     !
     !-----------------------------------------------------------------------
     SUBROUTINE constraint_grad( index, nat, tau, &
                                 if_pos, ityp, tau_units, g, dg )
       !-----------------------------------------------------------------------
       !
       ! ... this routine computes the value of the constraint equation and 
       ! ... the corresponding constraint gradient 
       !
       USE cell_base, ONLY : tpiba
       
       IMPLICIT NONE
       !
       INTEGER,  INTENT(IN)  :: index
       INTEGER,  INTENT(IN)  :: nat
       REAL(DP), INTENT(IN)  :: tau(:,:)
       INTEGER,  INTENT(IN)  :: if_pos(:,:)
       INTEGER,  INTENT(IN)  :: ityp(:)
       REAL(DP), INTENT(IN)  :: tau_units
       REAL(DP), INTENT(OUT) :: dg(:,:)
       REAL(DP), INTENT(OUT) :: g
       !
       ! ... local variables
       !
       INTEGER     :: i, j, i_sin, i_cos
       INTEGER     :: ia, ia0, ia1, ia2, ia3, n_type_coord1
       REAL(DP)    :: d0(3), d1(3), d2(3)
       REAL(DP)    :: inv_den, fac
       REAL(DP)    :: C00, C01, C02, C11, C12, C22
       REAL(DP)    :: D01, D12, invD01, invD12
       REAL(DP)    :: smearing, r_c, r_max
       INTEGER     :: type_coord1, type_coord2
       REAL(DP)    :: dtau(3), norm_dtau, expo
       REAL(DP)    :: k(3), phase, sin_phase, norm_k, sinxx
       COMPLEX(DP) :: struc_fac
       !
       ! ... external function
       !
       REAL(DP), EXTERNAL :: DDOT
       !
       !
       dg(:,:) = 0.D0
       !
       SELECT CASE ( constr_type(index) )
       CASE( 1 )
          !
          ! ... constraint on global coordination
          !
          type_coord1 = ANINT( constr(1,index) )
          type_coord2 = ANINT( constr(2,index) )
          !
          r_c = constr(3,index)
          !
          smearing = 1.D0 / constr(4,index)
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
                dtau(:) = pbc( tau(:,ia1) - tau(:,ia2) ) * tau_units
                !
                norm_dtau = norm( dtau(:) )
                !
                dtau(:) = dtau(:) / norm_dtau
                !
                expo = EXP( smearing * ( norm_dtau - r_c ) )
                !
                g = g + 1.D0 / ( expo + 1.D0 )
                !
                dtau(:) = dtau(:) * smearing * expo / ( expo + 1.D0 )**2
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
          g = ( g - target(index) )
          !
       CASE( 2 )
          !
          ! ... constraint on local coordination
          !
          ia          = ANINT( constr(1,index) )
          type_coord1 = ANINT( constr(2,index) )
          !
          r_c = constr(3,index)
          !
          smearing = 1.D0 / constr(4,index)
          !
          g = 0.D0
          !
          DO ia1 = 1, nat
             !
             IF ( ia1 == ia ) CYCLE
             !
             IF ( ityp(ia1) /= type_coord1 ) CYCLE
             !
             dtau(:) = pbc( tau(:,ia) - tau(:,ia1) ) * tau_units
             !
             norm_dtau = norm( dtau(:) )
             !
             dtau(:) = dtau(:) / norm_dtau
             !
             expo = EXP( smearing * ( norm_dtau - r_c ) )
             !
             g = g + 1.D0 / ( expo + 1.D0 )
             !
             dtau(:) = dtau(:) * smearing * expo / ( expo + 1.D0 )**2
             !
             dg(:,ia1) = dg(:,ia1) + dtau(:)
             dg(:,ia)  = dg(:,ia)  - dtau(:)
             !
          END DO
          !
          g = ( g - target(index) )
          !
       CASE( 3 )
          !
          ! ... constraint on distances
          !
          ia1 = ANINT( constr(1,index) )
          ia2 = ANINT( constr(2,index) )
          !
          dtau(:) = pbc( tau(:,ia1) - tau(:,ia2) ) * tau_units
          !
          norm_dtau = norm( dtau(:) )
          !
          g = ( norm_dtau - target(index) )
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
          ia0 = ANINT( constr(1,index) )
          ia1 = ANINT( constr(2,index) )
          ia2 = ANINT( constr(3,index) )
          !
          d0(:) = pbc( tau(:,ia0) - tau(:,ia1) ) * tau_units
          d1(:) = pbc( tau(:,ia1) - tau(:,ia2) ) * tau_units
          !
          C00 = d0(:) .dot. d0(:)
          C01 = d0(:) .dot. d1(:)
          C11 = d1(:) .dot. d1(:)
          !
          inv_den = 1.D0 / SQRT( C00 * C11 )
          !
          g = ( C01 * inv_den - target(index) )
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
          ia0 = ANINT( constr(1,index) )
          ia1 = ANINT( constr(2,index) )
          ia2 = ANINT( constr(3,index) )
          ia3 = ANINT( constr(4,index) )
          !
          d0(:) = pbc( tau(:,ia0) - tau(:,ia1) ) * tau_units
          d1(:) = pbc( tau(:,ia1) - tau(:,ia2) ) * tau_units
          d2(:) = pbc( tau(:,ia2) - tau(:,ia3) ) * tau_units
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
          g = ( ( C01 * C12 - C02 * C11 ) * inv_den - target(index) )
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
          k(1) = constr(1,index) * tpi
          k(2) = constr(2,index) * tpi
          k(3) = constr(3,index) * tpi
          !
          struc_fac = ( 1.D0, 0.D0 )
          !
          DO i = 1, nat - 1
             !
             dtau(:) = pbc( tau(:,i+1) - tau(:,1) )
             !
             phase = k(:) .dot. dtau(:)
             !
             struc_fac = struc_fac + CMPLX( COS( phase ), SIN( phase ) )
             !
             DO j = i + 1, nat
                !
                dtau(:) = pbc( tau(:,j) - tau(:,i) )
                !
                phase = k(:) .dot. dtau(:)
                !
                sin_phase = SIN( phase )
                !
                dg(:,i) = dg(:,i) + k(:) * sin_phase
                dg(:,j) = dg(:,j) - k(:) * sin_phase
                !
             END DO
             !
          END DO
          !
          g = ( CONJG( struc_fac ) * struc_fac ) / DBLE( nat )**2
          !
          g = ( g - target(index) )
          !
          dg(:,:) = 2.D0 * dg(:,:) / DBLE( nat )**2
          !
       CASE( 7 )
          !
          ! ... constraint on spherical average of the structure factor for
          ! ... a given k-vector of norm k (the 4*pi coefficient is understood)
          !
          norm_k = constr(1,index) * tpi / tau_units
          !
          g = 0.D0
          !
          DO i = 1, nat
             !
             dtau(:) = pbc( tau(:,i) - tau(:,1) ) * tau_units
             !
             norm_dtau = norm( dtau(:) )
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
                dg(:,i) = dtau(:) / norm_dtau**2 * ( COS( phase ) - sinxx )
                !
             END IF
             !
          END DO
          !
          g = ( g / DBLE( nat ) - target(index) )
          !
          dg(:,:) = dg(:,:) / DBLE( nat )
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
       ! ... update tau so that the constraint equation g=0 is satisfied,
       ! ... use the recursion formula:
       !
       ! ...                       g(taup)
       ! ... taup = taup - ----------------------- * dg(tau0)
       ! ...               M^-1<dg(taup)|dg(tau0)>
       !
       ! ... in normal cases the constraint equation should be always 
       ! ... satisfied at the very first iteration.
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
       INTEGER               :: na, i, index, dim
       REAL(DP), ALLOCATABLE :: gp(:), dgp(:,:), dg0(:,:,:)
       REAL(DP), ALLOCATABLE :: norm_dg0(:)
       REAL(DP)              :: norm_dgp
       REAL(DP)              :: g0
       REAL(DP)              :: lambda, fac, invdtsq
       LOGICAL               :: ltest(nconstr), global_test
       INTEGER, PARAMETER    :: maxiter = 100
       !
       REAL(DP), EXTERNAL :: DDOT, DNRM2
       !
       !
       ALLOCATE( dgp( 3, nat ) )
       ALLOCATE( dg0( 3, nat, nconstr ) )
       !
       ALLOCATE( gp(       nconstr ) )
       ALLOCATE( norm_dg0( nconstr ) )
       !
       invdtsq  = 1.D0 / dt**2
       !
       dim = 3 * nat
       !
       DO index = 1, nconstr
          !
          CALL constraint_grad( index, nat, tau0, &
                                if_pos, ityp, tau_units, g0, dg0(:,:,index) )
          !
          norm_dg0(index) = DNRM2( dim, dg0(:,:,index), 1 )
          !
       END DO
       !
       outer_loop: DO i = 1, maxiter
          !
          inner_loop: DO index = 1, nconstr
             !
             ltest(index) = .FALSE.
             !
             CALL constraint_grad( index, nat, taup, &
                                   if_pos, ityp, tau_units, gp(index), dgp )
             !
             norm_dgp = DNRM2( dim, dgp(:,:), 1 )
             !
             ! ... check if gp = 0
             !
#if defined (__DEBUG_CONSTRAINTS)
             WRITE( stdout, '(2(2X,I3),3(2X,F12.8))' ) &
                 i, index, ABS( gp(index) ), norm_dg0(index), norm_dgp
#endif
             !
             IF ( ABS( gp(index) ) < constr_tol .OR. &
                  norm_dg0(index) < eps8 .OR. norm_dgp < eps8 ) THEN
                !
                ltest(index) = .TRUE.
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
             lambda = gp(index) / DDOT( dim, dgp, 1, dg0(:,:,index), 1 )
             !
             DO na = 1, nat
                !
                fac = amass(ityp(na)) * massconv * tau_units
                !
                taup(:,na) = taup(:,na) - lambda * dg0(:,na,index) / fac
                !
             END DO
             !
             force(:,:) = force(:,:) - lambda * dg0(:,:,index) * invdtsq
             !
             lagrange(index) = lagrange(index) + lambda * invdtsq
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
       DEALLOCATE( norm_dg0 )
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
       REAL(DP)              :: g, norm_dg, norm_before, norm_after
       REAL(DP), ALLOCATABLE :: dg(:,:,:)
       REAL(DP), ALLOCATABLE :: dg_matrix(:,:)
       INTEGER,  ALLOCATABLE :: iwork(:)
       LOGICAL,  ALLOCATABLE :: skip(:)
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
       norm_before = DNRM2( 3 * nat, force, 1 )
       !
       ALLOCATE( dg( 3, nat, nconstr ) )
       ALLOCATE( dg_matrix( nconstr, nconstr ) )
       ALLOCATE( iwork( nconstr ) )
       ALLOCATE( skip(  nconstr ) )
       !
       DO i = 1, nconstr
          !
          CALL constraint_grad( i, nat, tau, &
                                if_pos, ityp, tau_units, g, dg(:,:,i) )
          !
          norm_dg = DNRM2( dim, dg(:,:,i), 1 )
          !
          IF ( norm_dg > eps8 ) THEN
             !
             dg(:,:,i) = dg(:,:,i) / norm_dg
             !
          ELSE
             !
             skip(i) = .TRUE.
             !
          END IF
          !
       END DO
       !
       DO i = 1, nconstr
          !
          IF ( skip(i) ) THEN
             !
             ! ... special case:  do not project out the components of the 
             ! ...                force along the constraint gradient if the 
             ! ...                latter is too small
             !
             dg_matrix(i,i) = 1.D0
             !
             lagrange(i) = 0.D0
             !
             DO j = i + 1, nconstr
                !
                dg_matrix(i,j) = 0.D0
                dg_matrix(j,i) = 0.D0
                !
             END DO
             !
          ELSE
             !
             dg_matrix(i,i) = 1.D0
             !
             lagrange(i) = DDOT( dim, force, 1, dg(:,:,i), 1 )
             !
             DO j = i + 1, nconstr
                !
                dg_matrix(i,j) = DDOT( dim, dg(:,:,i), 1, dg(:,:,j), 1 )
                dg_matrix(j,i) = dg_matrix(i,j)
                !
             END DO
             !
          END IF
          !
       END DO
       !
       IF ( nconstr == 1 ) THEN
          !
          force(:,:) = force(:,:) - lagrange(1) * dg(:,:,1)
          !
       ELSE
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
       DEALLOCATE( dg_matrix )
       DEALLOCATE( iwork )
       DEALLOCATE( skip )
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
       ! ... in cartesian units )
       !
       USE cell_base, ONLY : at, bg
       !
       IMPLICIT NONE
       !
       REAL(DP), INTENT(IN) :: vect(3)
       REAL(DP)             :: pbc(3)
       !
       !
#if defined (__USE_PBC)
       !
#if defined (__USE_LJ_PBC)
       !
       pbc(:) = vect(:) - ANINT( vect(:) )
       !
#else
       !
       pbc(:) = MATMUL( vect(:), bg(:,:) )
       !
       pbc(:) = pbc(:) - ANINT( pbc(:) )
       !
       pbc(:) = MATMUL( at(:,:), pbc(:) )
       !
#endif
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
       USE cell_base, ONLY : at, alat
       !
       IMPLICIT NONE
       !
       dmax = norm( MATMUL( at(:,:), (/ 0.5D0, 0.D0, 0.D0 /) ) ) * alat
       !
       dmax = MIN( dmax, &
                   norm( MATMUL( at(:,:), (/ 0.0D0, 0.5D0, 0.0D0 /) ) ) * alat )
       !
       dmax = MIN( dmax, &
                   norm( MATMUL( at(:,:), (/ 0.0D0, 0.0D0, 0.5D0 /) ) ) * alat )
       !
       RETURN
       !
     END SUBROUTINE compute_dmax
     !
END MODULE constraints_module
