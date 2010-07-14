!
! Copyright (C) 2002-2005 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
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
            remove_constr_vec,     &
            deallocate_constraint, &
            compute_dmax,          &
            pbc,                   &
            constraint_grad
  !
  !
  ! ... public variables (assigned in the CONSTRAINTS input card)
  !
  PUBLIC :: nconstr,       &
            constr_tol,    &
            constr_type,   &
            constr,        &
            lagrange,      &
            constr_target, &
            dmax, &
            gp
  !
  ! ... global variables
  !
  INTEGER               :: nconstr=0
  REAL(DP)              :: constr_tol
  INTEGER,  ALLOCATABLE :: constr_type(:)
  REAL(DP), ALLOCATABLE :: constr(:,:)
  REAL(DP), ALLOCATABLE :: constr_target(:)
  REAL(DP), ALLOCATABLE :: lagrange(:)
  REAL(DP), ALLOCATABLE :: gp(:)
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
                                    constr_target_inp, &
                                    constr_target_set, colvar_target_inp, &
                                    colvar_target_set, nc_fields
       !
       IMPLICIT NONE
       !
       INTEGER,  INTENT(in) :: nat
       REAL(DP), INTENT(in) :: tau(3,nat)
       INTEGER,  INTENT(in) :: ityp(nat)
       REAL(DP), INTENT(in) :: tau_units
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
       CHARACTER(len=6), EXTERNAL :: int_to_char
       !
       !
       nconstr    = ncolvar_inp + nconstr_inp
       !
       ! Be careful about what tolerance we want
       !
       IF (( ncolvar_inp > 0 ) .and. ( nconstr_inp == 0)) &
       constr_tol = colvar_tol_inp
       IF (( ncolvar_inp == 0 ) .and. ( nconstr_inp > 0)) &
       constr_tol = constr_tol_inp
       IF (( ncolvar_inp > 0 ) .and. ( nconstr_inp > 0)) &
       constr_tol = max( constr_tol_inp, colvar_tol_inp )
       !
       ALLOCATE( lagrange(      nconstr ) )
       ALLOCATE( constr_target( nconstr ) )
       ALLOCATE( constr_type(   nconstr ) )
       !
       ALLOCATE( constr( nc_fields, nconstr ) )
       ALLOCATE( gp(          nconstr ) )
       !
       ! ... setting constr to 0 to findout which elements have been
       ! ... set to an atomic index. This is required for CP.
       !
       constr(:,:) = 0.0_DP
       !
       ! ... NB: the first "ncolvar" constraints are collective variables (used
       ! ...     for meta-dynamics and free-energy smd), the remaining are real
       ! ...     constraints
       !
       IF (ncolvar_inp > 0) &
          constr(:,1:ncolvar_inp)        = colvar_inp(:,1:ncolvar_inp)
       IF (nconstr_inp > 0) &
          constr(:,ncolvar_inp+1:nconstr) = constr_inp(:,1:nconstr_inp)
       !
       ! ... set the largest possible distance among two atoms within
       ! ... the supercell
       !
       IF ( ncolvar_inp > 0 ) THEN
          !
          IF ( any( colvar_type_inp(:) == 'distance' ) ) CALL compute_dmax()
          !
       ELSEIF ( nconstr_inp > 0 ) THEN
          !
          IF ( any( constr_type_inp(:) == 'distance' ) ) CALL compute_dmax()
          !
       ENDIF
       !
       ! ... initializations of constr_target values for the constraints :
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
                constr_target(ia) = colvar_target_inp(ia)
                !
             ELSE
                !
                CALL set_type_coord( ia )
                !
             ENDIF
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
                constr_target(ia) = colvar_target_inp(ia)
                !
             ELSE
                !
                CALL set_atom_coord( ia )
                !
             ENDIF
             !
          CASE( 'distance' )
             !
             constr_type(ia) = 3
             !
             IF ( colvar_target_set(ia) ) THEN
                !
                constr_target(ia) = colvar_target_inp(ia)
                !
             ELSE
                !
                ia1 = anint( constr(1,ia) )
                ia2 = anint( constr(2,ia) )
                !
                dtau(:) = pbc( ( tau(:,ia1) - tau(:,ia2) )*tau_units )
                !
                constr_target(ia) = norm( dtau(:) )
                !
             ENDIF
             !
             IF ( constr_target(ia) > dmax )  THEN
                !
                WRITE( stdout, '(/,5X,"target = ",F12.8,/, &
                                &  5X,"dmax   = ",F12.8)' ) &
                    constr_target(ia), dmax
                !
                CALL errore( 'init_constraint', 'the target for coll.var. '  //&
                           & trim( int_to_char( ia ) ) // ' is larger than ' //&
                           & 'the largest possible value', 1 )
                !
             ENDIF
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
                ! ... the input value of target for the torsional angle (given
                ! ... in degrees) is converted to the cosine of the angle
                !
                constr_target(ia) = cos( ( 180.0_DP - &
                                           colvar_target_inp(ia) )*tpi/360.0_DP )
                !
             ELSE
                !
                CALL set_planar_angle( ia )
                !
             ENDIF
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
                ! ... the input value of target for the torsional angle (given
                ! ... in degrees) is converted to the cosine of the angle
                !
                constr_target(ia) = cos( colvar_target_inp(ia)*tpi/360.0_DP )
                !
             ELSE
                !
                CALL set_torsional_angle( ia )
                !
             ENDIF
             !
          CASE( 'struct_fac' )
             !
             ! ... constraint on structure factor at a given k-vector
             !
             constr_type(ia) = 6
             !
             IF ( colvar_target_set(ia) ) THEN
                !
                constr_target(ia) = colvar_target_inp(ia)
                !
             ELSE
                !
                CALL set_structure_factor( ia )
                !
             ENDIF
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
                constr_target(ia) = colvar_target_inp(ia)
                !
             ELSE
                !
                CALL set_sph_structure_factor( ia )
                !
             ENDIF
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
                constr_target(ia) = colvar_target_inp(ia)
                !
             ELSE
                !
                CALL set_bennett_proj( ia )
                !
             ENDIF
             !
          CASE DEFAULT
             !
             CALL errore( 'init_constraint', &
                          'collective-variable type not implemented', 1 )
             !
          END SELECT
          !
       ENDDO
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
                constr_target(ia) = constr_target_inp(n)
                !
             ELSE
                !
                CALL set_type_coord( ia )
                !
             ENDIF
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
                constr_target(ia) = constr_target_inp(n)
                !
             ELSE
                !
                CALL set_atom_coord( ia )
                !
             ENDIF
             !
          CASE( 'distance' )
             !
             constr_type(ia) = 3
             !
             IF ( constr_target_set(n) ) THEN
                !
                constr_target(ia) = constr_target_inp(n)
                !
             ELSE
                !
                ia1 = anint( constr(1,ia) )
                ia2 = anint( constr(2,ia) )
                !
                dtau(:) = pbc( ( tau(:,ia1) - tau(:,ia2) )*tau_units )
                !
                constr_target(ia) = norm( dtau(:) )
                !
             ENDIF
             !
             IF ( constr_target(ia) > dmax )  THEN
                !
                WRITE( stdout, '(/,5X,"target = ",F12.8,/, &
                                &  5X,"dmax   = ",F12.8)' ) &
                    constr_target(ia), dmax
                !
                CALL errore( 'init_constraint', 'the target for constraint ' //&
                           & trim( int_to_char( n ) ) // ' is larger than '  //&
                           & 'the largest possible value', 1 )
                !
             ENDIF
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
                ! ... the input value of target for the torsional angle (given
                ! ... in degrees) is converted to the cosine of the angle
                !
                constr_target(ia) = cos( ( 180.0_DP - &
                                           constr_target_inp(n) )*tpi/360.0_DP )
                !
             ELSE
                !
                CALL set_planar_angle( ia )
                !
             ENDIF
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
                ! ... the input value of target for the torsional angle (given
                ! ... in degrees) is converted to the cosine of the angle
                !
                constr_target(ia) = cos( constr_target_inp(n)*tpi/360.0_DP )
                !
             ELSE
                !
                CALL set_torsional_angle( ia )
                !
             ENDIF
             !
          CASE( 'struct_fac' )
             !
             ! ... constraint on structure factor at a given k-vector
             !
             constr_type(ia) = 6
             !
             IF ( constr_target_set(n) ) THEN
                !
                constr_target(ia) = constr_target_inp(n)
                !
             ELSE
                !
                CALL set_structure_factor( ia )
                !
             ENDIF
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
                constr_target(ia) = constr_target_inp(n)
                !
             ELSE
                !
                CALL set_sph_structure_factor( ia )
                !
             ENDIF
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
                constr_target(ia) = constr_target_inp(n)
                !
             ELSE
                !
                CALL set_bennett_proj( ia )
                !
             ENDIF
             !
          CASE DEFAULT
             !
             CALL errore( 'init_constraint', &
                          'constraint type not implemented', 1 )
             !
          END SELECT
          !
       ENDDO
       !
       RETURN
       !
       CONTAINS
         !
         !-------------------------------------------------------------------
         SUBROUTINE set_type_coord( ia )
           !-------------------------------------------------------------------
           !
           INTEGER, INTENT(in) :: ia
           !
           type_coord1 = anint( constr(1,ia) )
           type_coord2 = anint( constr(2,ia) )
           !
           r_c  = constr(3,ia)
           !
           smoothing = 1.0_DP / constr(4,ia)
           !
           constr_target(ia) = 0.0_DP
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
                 dtau(:) = pbc( ( tau(:,ia1) - tau(:,ia2) )*tau_units )
                 !
                 norm_dtau = norm( dtau(:) )
                 !
                 constr_target(ia) = constr_target(ia) + 1.0_DP / &
                                ( exp( smoothing*( norm_dtau - r_c ) ) + 1.0_DP )
                 !
              ENDDO
              !
              n_type_coord1 = n_type_coord1 + 1
              !
           ENDDO
           !
           constr_target(ia) = constr_target(ia) / dble( n_type_coord1 )
           !
         END SUBROUTINE set_type_coord
         !
         !-------------------------------------------------------------------
         SUBROUTINE set_atom_coord( ia )
           !-------------------------------------------------------------------
           !
           INTEGER, INTENT(in) :: ia
           !
           ia1         = anint( constr(1,ia) )
           type_coord1 = anint( constr(2,ia) )
           !
           r_c = constr(3,ia)
           !
           smoothing = 1.0_DP / constr(4,ia)
           !
           constr_target(ia) = 0.0_DP
           !
           DO ia2 = 1, nat
              !
              IF ( ia2 == ia1 ) CYCLE
              !
              IF ( ityp(ia2) /= type_coord1 ) CYCLE
              !
              dtau(:) = pbc( ( tau(:,ia1) - tau(:,ia2) )*tau_units )
              !
              norm_dtau = norm( dtau(:) )
              !
              constr_target(ia) = constr_target(ia) + 1.0_DP / &
                                ( exp( smoothing*( norm_dtau - r_c ) ) + 1.0_DP )
              !
           ENDDO
           !
         END SUBROUTINE set_atom_coord
         !
         !-------------------------------------------------------------------
         SUBROUTINE set_planar_angle( ia )
           !-------------------------------------------------------------------
           !
           INTEGER, INTENT(in) :: ia
           !
           ia0 = anint( constr(1,ia) )
           ia1 = anint( constr(2,ia) )
           ia2 = anint( constr(3,ia) )
           !
           d0(:) = pbc( ( tau(:,ia0) - tau(:,ia1) )*tau_units )
           d1(:) = pbc( ( tau(:,ia1) - tau(:,ia2) )*tau_units )
           !
           d0(:) = d0(:) / norm( d0(:) )
           d1(:) = d1(:) / norm( d1(:) )
           !
           constr_target(ia) = d0(:) .dot. d1(:)
           !
         END SUBROUTINE set_planar_angle
         !
         !-------------------------------------------------------------------
         SUBROUTINE set_torsional_angle( ia )
           !-------------------------------------------------------------------
           !
           INTEGER, INTENT(in) :: ia
           !
           ia0 = anint( constr(1,ia) )
           ia1 = anint( constr(2,ia) )
           ia2 = anint( constr(3,ia) )
           ia3 = anint( constr(4,ia) )
           !
           d0(:) = pbc( ( tau(:,ia0) - tau(:,ia1) )*tau_units )
           d1(:) = pbc( ( tau(:,ia1) - tau(:,ia2) )*tau_units )
           d2(:) = pbc( ( tau(:,ia2) - tau(:,ia3) )*tau_units )
           !
           C00 = d0(:) .dot. d0(:)
           C01 = d0(:) .dot. d1(:)
           C11 = d1(:) .dot. d1(:)
           C02 = d0(:) .dot. d2(:)
           C12 = d1(:) .dot. d2(:)
           C22 = d2(:) .dot. d2(:)
           !
           D01 = C00*C11 - C01*C01
           D12 = C11*C22 - C12*C12
           !
           constr_target(ia) = ( C01*C12 - C02*C11 ) / sqrt( D01*D12 )
           !
         END SUBROUTINE set_torsional_angle
         !
         !-------------------------------------------------------------------
         SUBROUTINE set_structure_factor( ia )
           !-------------------------------------------------------------------
           !
           INTEGER, INTENT(in) :: ia
           !
           k(1) = constr(1,ia) * tpi / tau_units
           k(2) = constr(2,ia) * tpi / tau_units
           k(3) = constr(3,ia) * tpi / tau_units
           !
           struc_fac = ( 0.0_DP, 0.0_DP )
           !
           DO i = 1, nat
              !
              dtau(:) = pbc( ( tau(:,i) - tau(:,1) )*tau_units )
              !
              phase = k(:) .dot. dtau(:)
              !
              struc_fac = struc_fac + cmplx( cos(phase), sin(phase), kind=DP )
              !
           ENDDO
           !
           constr_target(ia) = conjg( struc_fac )*struc_fac / dble( nat*nat )
           !
         END SUBROUTINE set_structure_factor
         !
         !-------------------------------------------------------------------
         SUBROUTINE set_sph_structure_factor( ia )
           !-------------------------------------------------------------------
           !
           INTEGER, INTENT(in) :: ia
           !
           norm_k = constr(1,ia)*tpi/tau_units
           !
           constr_target(ia) = 0.0_DP
           !
           DO i = 1, nat - 1
              !
              DO j = i + 1, nat
                 !
                 dtau(:) = pbc( ( tau(:,i) - tau(:,j) )*tau_units )
                 !
                 norm_dtau = norm( dtau(:) )
                 !
                 phase = norm_k*norm_dtau
                 !
                 IF ( phase < eps32 ) THEN
                    !
                    constr_target(ia) = constr_target(ia) + 1.0_DP
                    !
                 ELSE
                    !
                    constr_target(ia) = constr_target(ia) + sin( phase ) / phase
                    !
                 ENDIF
                 !
              ENDDO
              !
           ENDDO
           !
           constr_target(ia) = 2.0_DP * fpi * constr_target(ia) / dble( nat )
           !
         END SUBROUTINE set_sph_structure_factor
         !
         !-------------------------------------------------------------------
         SUBROUTINE set_bennett_proj( ia )
           !-------------------------------------------------------------------
           !
           INTEGER, INTENT(in) :: ia
           !
           ia0 = anint( constr(1,ia) )
           !
           d0(:) = tau(:,ia0)
           d1(:) = sum( tau(:,:), dim = 2 )
           !
           d1(:) = pbc( ( d1(:) - d0(:) )*tau_units ) / dble( nat - 1 ) - &
                   pbc( d0(:)*tau_units )
           !
           d2(:) = constr(2:4,ia)
           !
           constr_target(ia) = ( d1(:) .dot. d2(:) ) / tau_units
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
       INTEGER,  INTENT(in)  :: idx
       INTEGER,  INTENT(in)  :: nat
       REAL(DP), INTENT(in)  :: tau(:,:)
       INTEGER,  INTENT(in)  :: if_pos(:,:)
       INTEGER,  INTENT(in)  :: ityp(:)
       REAL(DP), INTENT(in)  :: tau_units
       REAL(DP), INTENT(out) :: dg(:,:)
       REAL(DP), INTENT(out) :: g
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
       REAL(DP), EXTERNAL :: ddot
       !
       !
       dg(:,:) = 0.0_DP
       !
       SELECT CASE ( constr_type(idx) )
       CASE( 1 )
          !
          ! ... constraint on global coordination
          !
          type_coord1 = anint( constr(1,idx) )
          type_coord2 = anint( constr(2,idx) )
          !
          r_c = constr(3,idx)
          !
          smoothing = 1.0_DP / constr(4,idx)
          !
          g = 0.0_DP
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
                dtau(:) = pbc( ( tau(:,ia1) - tau(:,ia2) )*tau_units )
                !
                norm_dtau = norm( dtau(:) )
                !
                dtau(:) = dtau(:) / norm_dtau
                !
                expo = exp( smoothing*( norm_dtau - r_c ) )
                !
                g = g + 1.0_DP / ( expo + 1.0_DP )
                !
                dtau(:) = dtau(:) * smoothing*expo / ( expo + 1.0_DP )**2
                !
                dg(:,ia2) = dg(:,ia2) + dtau(:)
                dg(:,ia1) = dg(:,ia1) - dtau(:)
                !
             ENDDO
             !
             n_type_coord1 = n_type_coord1 + 1
             !
          ENDDO
          !
          g  = g  / dble( n_type_coord1 )
          dg = dg / dble( n_type_coord1 )
          !
          g = ( g - constr_target(idx) )
          !
       CASE( 2 )
          !
          ! ... constraint on local coordination
          !
          ia          = anint( constr(1,idx) )
          type_coord1 = anint( constr(2,idx) )
          !
          r_c = constr(3,idx)
          !
          smoothing = 1.0_DP / constr(4,idx)
          !
          g = 0.0_DP
          !
          DO ia1 = 1, nat
             !
             IF ( ia1 == ia ) CYCLE
             !
             IF ( ityp(ia1) /= type_coord1 ) CYCLE
             !
             dtau(:) = pbc( ( tau(:,ia) - tau(:,ia1) )*tau_units )
             !
             norm_dtau = norm( dtau(:) )
             !
             dtau(:) = dtau(:) / norm_dtau
             !
             expo = exp( smoothing*( norm_dtau - r_c ) )
             !
             g = g + 1.0_DP / ( expo + 1.0_DP )
             !
             dtau(:) = dtau(:) * smoothing * expo / ( expo + 1.0_DP )**2
             !
             dg(:,ia1) = dg(:,ia1) + dtau(:)
             dg(:,ia)  = dg(:,ia)  - dtau(:)
             !
          ENDDO
          !
          g = ( g - constr_target(idx) )
          !
       CASE( 3 )
          !
          ! ... constraint on distances
          !
          ia1 = anint( constr(1,idx) )
          ia2 = anint( constr(2,idx) )
          !
          dtau(:) = pbc( ( tau(:,ia1) - tau(:,ia2) )*tau_units )
          !
          norm_dtau = norm( dtau(:) )
          !
          g = ( norm_dtau - constr_target(idx) )
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
          ia0 = anint( constr(1,idx) )
          ia1 = anint( constr(2,idx) )
          ia2 = anint( constr(3,idx) )
          !
          d0(:) = pbc( ( tau(:,ia0) - tau(:,ia1) )*tau_units )
          d1(:) = pbc( ( tau(:,ia1) - tau(:,ia2) )*tau_units )
          !
          C00 = d0(:) .dot. d0(:)
          C01 = d0(:) .dot. d1(:)
          C11 = d1(:) .dot. d1(:)
          !
          inv_den = 1.0_DP / sqrt( C00*C11 )
          !
          g = ( C01 * inv_den - constr_target(idx) )
          !
          dg(:,ia0) = ( d1(:) - C01/C00*d0(:) ) * inv_den
          dg(:,ia2) = ( C01/C11*d1(:) - d0(:) ) * inv_den
          dg(:,ia1) = - dg(:,ia0) - dg(:,ia2)
          !
       CASE( 5 )
          !
          ! ... constraint on torsional angle (for the notation used here
          ! ... see Appendix C of the Allen-Tildesley book)
          !
          ia0 = anint( constr(1,idx) )
          ia1 = anint( constr(2,idx) )
          ia2 = anint( constr(3,idx) )
          ia3 = anint( constr(4,idx) )
          !
          d0(:) = pbc( ( tau(:,ia0) - tau(:,ia1) )*tau_units )
          d1(:) = pbc( ( tau(:,ia1) - tau(:,ia2) )*tau_units )
          d2(:) = pbc( ( tau(:,ia2) - tau(:,ia3) )*tau_units )
          !
          C00 = d0(:) .dot. d0(:)
          C01 = d0(:) .dot. d1(:)
          C11 = d1(:) .dot. d1(:)
          C02 = d0(:) .dot. d2(:)
          C12 = d1(:) .dot. d2(:)
          C22 = d2(:) .dot. d2(:)
          !
          D01 = C00*C11 - C01*C01
          D12 = C11*C22 - C12*C12
          !
          IF ( abs( D01 ) < eps32 .or. abs( D12 ) < eps32 ) &
             CALL errore( 'constraint_grad', 'either D01 or D12 is zero', 1 )
          !
          invD01 = 1.0_DP / D01
          invD12 = 1.0_DP / D12
          !
          fac = C01*C12 - C02*C11
          !
          inv_den = 1.0_DP / sqrt( D01*D12 )
          !
          g = ( ( C01*C12 - C02*C11 )*inv_den - constr_target(idx) )
          !
          dg(:,ia0) = ( C12*d1(:) - C11*d2(:) - &
                        invD01*fac*( C11*d0(:) - C01*d1(:) ) )*inv_den
          !
          dg(:,ia2) = ( C01*( d1(:) - d2(:) ) - &
                        ( C11 + C12 )*d0(:) + 2.0_DP*C02*d1(:) - &
                        invD12*fac*( ( C11 + C12 )*d2(:) - &
                                         ( C12 + C22 )*d1(:) ) - &
                        invD01*fac*( C01*d0(:) - C00*d1(:) ) )*inv_den
          !
          dg(:,ia3) = ( C11*d0(:) - C01*d1(:) - &
                        invD12*fac*( C12*d1(:) - C11*d2(:) ) )*inv_den
          !
          dg(:,ia1) = - dg(:,ia0) - dg(:,ia2) - dg(:,ia3)
          !
       CASE( 6 )
          !
          ! ... constraint on structure factor at a given k vector
          !
          k(1) = constr(1,idx)*tpi/tau_units
          k(2) = constr(2,idx)*tpi/tau_units
          k(3) = constr(3,idx)*tpi/tau_units
          !
          struc_fac = ( 1.0_DP, 0.0_DP )
          !
          r0(:) = tau(:,1)
          !
          DO i = 1, nat - 1
             !
             dtau(:) = pbc( ( tau(:,i+1) - r0(:) )*tau_units )
             !
             phase = k(1)*dtau(1) + k(2)*dtau(2) + k(3)*dtau(3)
             !
             struc_fac = struc_fac + cmplx( cos(phase), sin(phase), kind=DP )
             !
             ri(:) = tau(:,i)
             !
             DO j = i + 1, nat
                !
                dtau(:) = pbc( ( tau(:,j) - ri(:) )*tau_units )
                !
                phase = k(1)*dtau(1) + k(2)*dtau(2) + k(3)*dtau(3)
                !
                ksin(:) = k(:)*sin( phase )
                !
                dg(:,i) = dg(:,i) + ksin(:)
                dg(:,j) = dg(:,j) - ksin(:)
                !
             ENDDO
             !
          ENDDO
          !
          g = ( conjg( struc_fac )*struc_fac ) / dble( nat*nat )
          !
          g = ( g - constr_target(idx) )
          !
          dg(:,:) = dg(:,:)*2.0_DP/dble( nat*nat )
          !
       CASE( 7 )
          !
          ! ... constraint on spherical average of the structure factor for
          ! ... a given k-vector of norm k
          !
          norm_k = constr(1,idx)*tpi/tau_units
          !
          g = 0.0_DP
          !
          DO i = 1, nat - 1
             !
             ri(:) = tau(:,i)
             !
             DO j = i + 1, nat
                !
                dtau(:) = pbc( ( ri(:) - tau(:,j) )*tau_units )
                !
                norm_dtau_sq = dtau(1)**2 + dtau(2)**2 + dtau(3)**2
                !
                norm_dtau = sqrt( norm_dtau_sq )
                !
                phase = norm_k * norm_dtau
                !
                IF ( phase < eps32 ) THEN
                   !
                   g = g + 1.0_DP
                   !
                ELSE
                   !
                   sinxx = sin( phase ) / phase
                   !
                   g = g + sinxx
                   !
                   dtau(:) = dtau(:) / norm_dtau_sq*( cos( phase ) - sinxx )
                   !
                   dg(:,i) = dg(:,i) + dtau(:)
                   dg(:,j) = dg(:,j) - dtau(:)
                   !
                ENDIF
                !
             ENDDO
             !
          ENDDO
          !
          g = ( 2.0_DP*fpi*g / dble( nat ) - constr_target(idx) )
          !
          dg(:,:) = 4.0_DP*fpi*dg(:,:) / dble( nat )
          !
       CASE( 8 )
          !
          ! ... constraint on Bennett projection
          !
          ia0 = anint( constr(1,idx) )
          !
          d0(:) = tau(:,ia0)
          d1(:) = sum( tau(:,:), dim = 2 )
          !
          d1(:) = pbc( ( d1(:) - d0(:) )*tau_units ) / dble( nat - 1 ) - &
                  pbc( d0(:)*tau_units )
          !
          d2(:) = constr(2:4,idx)
          !
          g = ( d1(:) .dot. d2(:) ) / tau_units - constr_target( idx )
          !
          dg = 0.0_DP
          !
          C00 = ( 1.0_DP / dble( nat - 1 ) ) / tau_units
          C01 = -1.0_DP / tau_units
          !
          DO i = 1, nat
             !
             dg(:,i) = d2(:)*C00
             !
          ENDDO
          !
          dg(:,ia0) = d2(:)*C01
          !
       END SELECT
       !
       dg(:,:) = dg(:,:)*dble( if_pos(:,:) )
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
       INTEGER,  INTENT(in)    :: nat
       REAL(DP), INTENT(inout) :: taup(3,nat)
       REAL(DP), INTENT(in)    :: tau0(3,nat)
       INTEGER,  INTENT(in)    :: if_pos(3,nat)
       REAL(DP), INTENT(inout) :: force(3,nat)
       INTEGER,  INTENT(in)    :: ityp(nat)
       REAL(DP), INTENT(in)    :: tau_units
       REAL(DP), INTENT(in)    :: dt
       REAL(DP), INTENT(in)    :: massconv
       !
       INTEGER               :: na, i, idx, dim
       REAL(DP), ALLOCATABLE :: dgp(:,:), dg0(:,:,:)
       REAL(DP)              :: g0
       REAL(DP)              :: lambda, fac, invdtsq
       LOGICAL, ALLOCATABLE  :: ltest(:)
       LOGICAL               :: global_test
       INTEGER, PARAMETER    :: maxiter = 100
       !
       REAL(DP), EXTERNAL :: ddot
       !
       !
       ALLOCATE( dgp( 3, nat ) )
       ALLOCATE( dg0( 3, nat, nconstr ) )
       !
!       ALLOCATE( gp(    nconstr ) )
       ALLOCATE( ltest( nconstr ) )
       !
       invdtsq  = 1.0_DP / dt**2
       !
       dim = 3*nat
       !
       DO idx = 1, nconstr
          !
          CALL constraint_grad( idx, nat, tau0, &
                                if_pos, ityp, tau_units, g0, dg0(:,:,idx) )
          !
       ENDDO
       !
       outer_loop: DO i = 1, maxiter
          !
          inner_loop: DO idx = 1, nconstr
             !
             ltest(idx) = .false.
             !
             CALL constraint_grad( idx, nat, taup, &
                                   if_pos, ityp, tau_units, gp(idx), dgp )
             !
             ! ... check if gp = 0
             !
#if defined (__DEBUG_CONSTRAINTS)
             WRITE( stdout, '(2(2X,I3),F12.8)' ) i, idx, abs( gp(idx) )
#endif
             !
             IF ( abs( gp(idx) ) < constr_tol ) THEN
                !
                ltest(idx) = .true.
                !
                CYCLE inner_loop
                !
             ENDIF
             !
             ! ... if  gp <> 0  find new taup and check again
             ! ... ( gp is in bohr and taup in tau_units )
             !
             DO na = 1, nat
                !
                dgp(:,na) = dgp(:,na) / ( amass(ityp(na))*massconv )
                !
             ENDDO
             !
             lambda = gp(idx) / ddot( dim, dgp, 1, dg0(:,:,idx), 1 )
             !
             DO na = 1, nat
                !
                fac = amass(ityp(na))*massconv*tau_units
                !
                taup(:,na) = taup(:,na) - lambda*dg0(:,na,idx)/fac
                !
             ENDDO
             !
             lagrange(idx) = lagrange(idx) + lambda*invdtsq
             !
             force(:,:) = force(:,:) - lambda*dg0(:,:,idx)*invdtsq
             !
          ENDDO inner_loop
          !
          global_test = all( ltest(:) )
          !
          ! ... all constraints are satisfied
          !
          IF ( global_test ) exit outer_loop
          !
       ENDDO outer_loop
       !
       IF ( .not. global_test ) THEN
          !
          ! ... error messages
          !
          WRITE( stdout, '(/,5X,"Number of step(s): ",I3)') min( i, maxiter )
          WRITE( stdout, '(/,5X,"constr_target convergence: ")' )
          !
          DO i = 1, nconstr
             !
             WRITE( stdout, '(5X,"constr # ",I3,2X,L1,3(2X,F16.10))' ) &
                 i, ltest(i), abs( gp(i) ), constr_tol, constr_target(i)
             !
          ENDDO
          !
          CALL errore( 'check_constraint', &
                       'on some constraint g = 0 is not satisfied', 1 )
          !
       ENDIF
       !
       DEALLOCATE( dgp )
       DEALLOCATE( dg0 )
!       DEALLOCATE( gp )
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
       INTEGER,  INTENT(in)    :: nat
       REAL(DP), INTENT(in)    :: tau(:,:)
       INTEGER,  INTENT(in)    :: if_pos(:,:)
       INTEGER,  INTENT(in)    :: ityp(:)
       REAL(DP), INTENT(in)    :: tau_units
       REAL(DP), INTENT(inout) :: force(:,:)
       !
       INTEGER               :: i, j, dim
       REAL(DP)              :: g, ndg, dgidgj
       REAL(DP)              :: norm_before, norm_after
       REAL(DP), ALLOCATABLE :: dg(:,:,:)
       REAL(DP), ALLOCATABLE :: dg_matrix(:,:)
       INTEGER,  ALLOCATABLE :: iwork(:)
       !
       REAL(DP), EXTERNAL :: ddot, dnrm2
       !
       !
       dim = 3*nat
       !
       lagrange(:) = 0.0_DP
       !
#if defined (__REMOVE_CONSTRAINT_FORCE)
       !
       norm_before = dnrm2( 3*nat, force, 1 )
       !
       ALLOCATE( dg( 3, nat, nconstr ) )
       !
       IF ( nconstr == 1 ) THEN
          !
          CALL constraint_grad( 1, nat, tau, &
                                if_pos, ityp, tau_units, g, dg(:,:,1) )
          !
          lagrange(1) = ddot( dim, force, 1, dg(:,:,1), 1 )
          !
          ndg = ddot( dim, dg(:,:,1), 1, dg(:,:,1), 1 )
          !
          force(:,:) = force(:,:) - lagrange(1)*dg(:,:,1)/ndg
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
          ENDDO
          !
          DO i = 1, nconstr
             !
             dg_matrix(i,i) = ddot( dim, dg(:,:,i), 1, dg(:,:,i), 1 )
             !
             lagrange(i) = ddot( dim, force, 1, dg(:,:,i), 1 )
             !
             DO j = i + 1, nconstr
                !
                dgidgj = ddot( dim, dg(:,:,i), 1, dg(:,:,j), 1 )
                !
                dg_matrix(i,j) = dgidgj
                dg_matrix(j,i) = dgidgj
                !
             ENDDO
             !
          ENDDO
          !
          CALL DGESV( nconstr, 1, dg_matrix, &
                      nconstr, iwork, lagrange, nconstr, i )
          !
          IF ( i /= 0 ) &
             CALL errore( 'remove_constr_force', &
                          'error in the solution of the linear system', i )
          !
          DO i = 1, nconstr
             !
             force(:,:) = force(:,:) - lagrange(i)*dg(:,:,i)
             !
          ENDDO
          !
          DEALLOCATE( dg_matrix )
          DEALLOCATE( iwork )
          !
       ENDIF
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
       ENDDO
       !
#endif
       !
       norm_after = dnrm2( dim, force, 1 )
       !
       IF ( norm_before < norm_after ) THEN
          !
          WRITE( stdout, '(/,5X,"norm before = ",F16.10)' ) norm_before
          WRITE( stdout, '(  5X,"norm after  = ",F16.10)' ) norm_after
          !
          CALL errore( 'remove_constr_force', &
                       'norm(F) before < norm(F) after', 1 )
          !
       ENDIF
       !
       DEALLOCATE( dg )
       !
#endif
       !
     END SUBROUTINE remove_constr_force
     !
     !-----------------------------------------------------------------------
     SUBROUTINE remove_constr_vec( nat, tau, &
                                   if_pos, ityp, tau_units, vec )
       !-----------------------------------------------------------------------
       !
       ! ... the component of a displacement vector that is orthogonal to the
       ! ... ipersurface defined by the constraint equations is removed
       ! ... and the corresponding value of the lagrange multiplier computed
       !
       IMPLICIT NONE
       !
       INTEGER,  INTENT(in)    :: nat
       REAL(DP), INTENT(in)    :: tau(:,:)
       INTEGER,  INTENT(in)    :: if_pos(:,:)
       INTEGER,  INTENT(in)    :: ityp(:)
       REAL(DP), INTENT(in)    :: tau_units
       REAL(DP), INTENT(inout) :: vec(:,:)
       !
       INTEGER               :: i, j, dim
       REAL(DP)              :: g, ndg, dgidgj
       REAL(DP), ALLOCATABLE :: dg(:,:,:), dg_matrix(:,:), lambda(:)
       INTEGER,  ALLOCATABLE :: iwork(:)
       !
       REAL(DP), EXTERNAL :: ddot, dnrm2
       !
       !
       dim = 3*nat
       !
       ALLOCATE( lambda( nconstr ) )
       ALLOCATE( dg( 3, nat, nconstr ) )
       !
       IF ( nconstr == 1 ) THEN
          !
          CALL constraint_grad( 1, nat, tau, &
                                if_pos, ityp, tau_units, g, dg(:,:,1) )
          !
          lambda(1) = ddot( dim, vec, 1, dg(:,:,1), 1 )
          !
          ndg = ddot( dim, dg(:,:,1), 1, dg(:,:,1), 1 )
          !
          vec(:,:) = vec(:,:) - lambda(1)*dg(:,:,1)/ndg
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
          ENDDO
          !
          DO i = 1, nconstr
             !
             dg_matrix(i,i) = ddot( dim, dg(:,:,i), 1, dg(:,:,i), 1 )
             !
             lambda(i) = ddot( dim, vec, 1, dg(:,:,i), 1 )
             !
             DO j = i + 1, nconstr
                !
                dgidgj = ddot( dim, dg(:,:,i), 1, dg(:,:,j), 1 )
                !
                dg_matrix(i,j) = dgidgj
                dg_matrix(j,i) = dgidgj
                !
             ENDDO
             !
          ENDDO
          !
          CALL DGESV( nconstr, 1, dg_matrix, &
                      nconstr, iwork, lambda, nconstr, i )
          !
          IF ( i /= 0 ) &
             CALL errore( 'remove_constr_vec', &
                          'error in the solution of the linear system', i )
          !
          DO i = 1, nconstr
             !
             vec(:,:) = vec(:,:) - lambda(i)*dg(:,:,i)
             !
          ENDDO
          !
          DEALLOCATE( dg_matrix )
          DEALLOCATE( iwork )
          !
       ENDIF
       !
       DEALLOCATE( lambda, dg )
       !
     END SUBROUTINE remove_constr_vec
     !
     !-----------------------------------------------------------------------
     SUBROUTINE deallocate_constraint()
       !-----------------------------------------------------------------------
       !
       IMPLICIT NONE
       !
       !
       IF ( allocated( lagrange ) )      DEALLOCATE( lagrange )
       IF ( allocated( constr ) )        DEALLOCATE( constr )
       IF ( allocated( constr_type ) )   DEALLOCATE( constr_type )
       IF ( allocated( constr_target ) ) DEALLOCATE( constr_target )
       IF ( allocated( gp     ) )      DEALLOCATE( gp )
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
       REAL(DP), INTENT(in) :: vect(3)
       REAL(DP)             :: pbc(3)
       !
       !
#if defined (__USE_PBC)
       !
       pbc(:) = matmul( vect(:), bg(:,:) )/alat
       !
       pbc(:) = pbc(:) - anint( pbc(:) )
       !
       pbc(:) = matmul( at(:,:), pbc(:) )*alat
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
       ! ... dmax corresponds to one half the longest diagonal of the cell
       !
       USE cell_base, ONLY : at, alat
       !
       IMPLICIT NONE
       !
       INTEGER  :: x,y,z
       REAL(DP) :: diago(3)
       !
       dmax = 0._dp !norm(at(:,1)+at(:,2)+at(:,3))
       !
       DO z = -1,1,2
       DO y = -1,1,2
       DO x = -1,1,2
          diago = x*at(:,1) + y*at(:,2) + z*at(:,3)
          dmax = max(dmax, norm(diago))
       ENDDO
       ENDDO
       ENDDO
       !
       dmax= dmax*alat*.5_dp
       !
       RETURN
       !
     END SUBROUTINE compute_dmax
     !
END MODULE constraints_module
