!
! Copyright (C) 2002 FPMD group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!


!  AB INITIO COSTANT PRESSURE MOLECULAR DYNAMICS
!  ----------------------------------------------
!  Car-Parrinello Parallel Program
!  Carlo Cavazzoni - Gerardo Ballabio
!  SISSA, Trieste, Italy - 1997-99
!  Last modified: Fri Nov  5 13:33:57 MET 1999

!  BEGIN manual

!=----------------------------------------------------------------------------=!
      MODULE phase_factors_module
!=----------------------------------------------------------------------------=!

!  (describe briefly what this module does...)
!  ----------------------------------------------
!  routines in this module:
!  SUBROUTINE phfacs(eigr,atoms)
!  SUBROUTINE nlpre(eigr,gv)
!  SUBROUTINE strucf(atoms,eigr,gv)
!  ----------------------------------------------

!  END manual

! ...   declare modules
        USE kinds

        IMPLICIT NONE
        SAVE

        PRIVATE

        PUBLIC :: strucf

! ...   --+ end of module-scope declarations +--

!=----------------------------------------------------------------------------=!
      CONTAINS
!=----------------------------------------------------------------------------=!

!=----------------------------------------------------------------------------=!
!  BEGIN manual

        SUBROUTINE phfacs(eigr, atoms)

!  this routine computes the phase factors
!
!    eigr%x(ix,ia) = exp (i ix G_1 dot R(ia))
!    eigr%y(iy,ia) = exp (i iy G_2 dot R(ia))
!    eigr%z(iz,ia) = exp (i iz G_3 dot R(ia))
!
!    G_1,G_2,G_3 = reciprocal lattice generators
!
!    ia = index of ion
!    ix,iy,iz = Miller indices
!  ----------------------------------------------
 
!  END manual

          USE constants, ONLY: tpi
          USE atoms_type_module, ONLY: atoms_type
          USE cp_types, ONLY: phase_factors

! ...     declare subroutine arguments
          TYPE (phase_factors) :: eigr
          TYPE (atoms_type) :: atoms  ! ionic positions

! ...     declare other variables
          COMPLEX(dbl) :: ctep1, ctep2, ctep3, ctem1, ctem2, ctem3
          REAL(dbl) :: ar1, ar2, ar3
          INTEGER :: i, j, k, isa, nr1, nr2, nr3

! ...     --+ end of declarations +--

          nr1 = UBOUND( eigr%x, 1 )
          nr2 = UBOUND( eigr%y, 1 )
          nr3 = UBOUND( eigr%z, 1 )
          

          DO isa = 1, atoms%nat

! ...         Miller index = 0: exp(i 0 dot R(ia)) = 1

              eigr%x( 0, isa ) = CMPLX( 1.d0, 0.d0 )
              eigr%y( 0, isa ) = CMPLX( 1.d0, 0.d0 )
              eigr%z( 0, isa ) = CMPLX( 1.d0, 0.d0 )

! ...         let R_1,R_2,R_3 be the direct lattice generators,
! ...         G_1,G_2,G_3 the reciprocal lattice generators
! ...         by definition G_i dot R_j = 2 pi delta_{ij}
! ...         ionic coordinates are in units of R_1,R_2,R_3
! ...         then G_i dot R(ia) = 2 pi R(ia)_i

              ar1 = tpi * atoms%taus(1,isa)  ! G_1 dot R(ia)
              ar2 = tpi * atoms%taus(2,isa)  ! G_2 dot R(ia)
              ar3 = tpi * atoms%taus(3,isa)  ! G_3 dot R(ia)

! ...         Miller index = 1: exp(i G_i dot R(ia))

              ctep1 = CMPLX( cos( ar1 ), -sin( ar1 ) )
              ctep2 = CMPLX( cos( ar2 ), -sin( ar2 ) )
              ctep3 = CMPLX( cos( ar3 ), -sin( ar3 ) )

! ...         Miller index = -1: exp(-i G_i dot R(ia))

              ctem1 = CONJG(ctep1)
              ctem2 = CONJG(ctep2)
              ctem3 = CONJG(ctep3)

! ...         Miller index > 0: exp(i N G_i dot R(ia)) =
! ...           = exp(i G_i dot R(ia)) exp(i (N-1) G_i dot R(ia))
! ...         Miller index < 0: exp(-i N G_i dot R(ia)) =
! ...           = exp(-i G_i dot R(ia)) exp(-i (N-1) G_i dot R(ia))

              DO i = 1, nr1 - 1
                eigr%x(  i, isa ) = eigr%x(  i - 1, isa ) * ctep1
                eigr%x( -i, isa ) = eigr%x( -i + 1, isa ) * ctem1
              END DO
              DO j = 1, nr2 - 1
                eigr%y(  j, isa ) = eigr%y(  j - 1, isa ) * ctep2
                eigr%y( -j, isa ) = eigr%y( -j + 1, isa ) * ctem2
              END DO
              DO k = 1, nr3 - 1
                eigr%z(  k, isa ) = eigr%z(  k - 1, isa ) * ctep3
                eigr%z( -k, isa ) = eigr%z( -k + 1, isa ) * ctem3
              END DO

          END DO

          RETURN
        END SUBROUTINE phfacs


!=----------------------------------------------------------------------------=!
!  BEGIN manual

        SUBROUTINE nlpre(eigr, gv)

!  this routine computes the phase factors
!
!    eigr%xyz(ig,ia) = exp(i G dot R(ia)) =
!                    = eigr%x(ix,ia) * eigr%y(iy,ia) * eigr%z(iz,ia)
!
!    eigr%x(ix,ia) = exp (i ix G_1 dot R(ia))
!    eigr%y(iy,ia) = exp (i iy G_2 dot R(ia))
!    eigr%z(iz,ia) = exp (i iz G_3 dot R(ia))
!
!    G_1,G_2,G_3 = reciprocal lattice generators
!
!    ia = index of ion
!    ig = index of G vector
!    ix,iy,iz = Miller indices of G vector
!  ----------------------------------------------
!  END manual

          USE cp_types, ONLY: phase_factors, recvecs

! ...     declare subroutine arguments
          TYPE (phase_factors) :: eigr
          TYPE (recvecs) :: gv

! ...     declare other variables
          INTEGER :: i, ig, ig1, ig2, ig3, nat, ngw

! ...     --+ end of declarations +--

          nat = SIZE( eigr%xyz, 2 )
          ngw = SIZE( eigr%xyz, 1 )
          IF( ngw > SIZE( gv%mill, 2 ) ) THEN
            CALL errore(' nlpre ',' eigr inconsisten size ',ngw)
          END IF

          DO ig = 1, ngw
            ig1 = gv%mill(1,ig)
            ig2 = gv%mill(2,ig)
            ig3 = gv%mill(3,ig)
            DO i = 1, nat
              eigr%xyz( ig, i ) = &
                 eigr%x( ig1, i ) * eigr%y( ig2, i ) * eigr%z( ig3, i )
            END DO
          END DO
          RETURN
        END SUBROUTINE nlpre

!=----------------------------------------------------------------------------=!
!  BEGIN manual

      SUBROUTINE strucf( sfac, atoms, eigr, gv )

!  this routine computes the structure factors
!
!    sfac(ig,is) = (sum over ia) exp(i G dot R(ia)) =
!                  (sum over ia) eigr%x(ix,ia) * eigr%y(iy,ia) * eigr%z(iz,ia)
!
!    eigr%x(ix,ia) = exp (i ix G_1 dot R(ia))
!    eigr%y(iy,ia) = exp (i iy G_2 dot R(ia))
!    eigr%z(iz,ia) = exp (i iz G_3 dot R(ia))
!
!    G_1,G_2,G_3 = reciprocal lattice generators
!
!    ia = index of ion (running over ions of species is)
!    ig = index of G vector
!    is = index of atomic species
!    ix,iy,iz = Miller indices of G vector
!  ----------------------------------------------
!  END manual

      USE atoms_type_module, ONLY: atoms_type
      USE cp_types, ONLY: phase_factors, recvecs

      IMPLICIT NONE

! ... declare subroutine arguments
      TYPE (atoms_type) :: atoms
      TYPE (phase_factors) :: eigr
      TYPE (recvecs) :: gv
      COMPLEX(dbl), INTENT(OUT) :: sfac(:,:)

! ... declare other variables
      INTEGER :: is, ig, ia, ig1, ig2, ig3, isa

! ...     --+ end of declarations +--
      
      CALL phfacs( eigr, atoms )
      CALL nlpre( eigr, gv )

      DO ig = 1, gv%ng_l
        ig1 = gv%mill( 1, ig ) 
        ig2 = gv%mill( 2, ig ) 
        ig3 = gv%mill( 3, ig )
        isa = 1
        DO is = 1, atoms%nsp
          sfac( is, ig ) = CMPLX (0.0_dbl, 0.0_dbl)
          DO ia = 1, atoms%na(is)
            sfac( is, ig ) = sfac( is, ig ) + &
              eigr%x( ig1, isa ) * eigr%y( ig2, isa ) * eigr%z( ig3, isa )
            isa = isa + 1
          END DO
        END DO
      END DO

      RETURN
      END SUBROUTINE strucf

!=----------------------------------------------------------------------------=!
      END MODULE phase_factors_module
!=----------------------------------------------------------------------------=!

