!
! Copyright (C) 2002 FPMD group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

!------------------------------------------------------------------------------!
  MODULE ions_base
!------------------------------------------------------------------------------!

      USE kinds, ONLY: dbl
      USE parameters, ONLY: nsx, natx
!
      IMPLICIT NONE
      SAVE

      !     nsp       = number of species
      !     na(is)    = number of atoms of species is
      !     nax       = max number of atoms of a given species
      !     nat       = total number of atoms of all species

      INTEGER :: nsp     = 0
      INTEGER :: na(nsx) = 0    
      INTEGER :: nax     = 0
      INTEGER :: nat     = 0
      INTEGER :: isort_pos(natx,nsx) = 0

      !     zv(is)    = (pseudo-)atomic charge
      !     pmass(is) = mass (converted to a.u.) of ions
      !     rcmax(is) = Ewald radius (for ion-ion interactions)

      REAL(dbl) :: zv(nsx)    = 0.0d0
      REAL(dbl) :: pmass(nsx) = 0.0d0
      REAL(dbl) :: rcmax(nsx) = 0.0d0

      INTEGER :: ipp(nsx) = 0
!

!------------------------------------------------------------------------------!
  CONTAINS
!------------------------------------------------------------------------------!

    SUBROUTINE packtau( taup, tau, na, nsp )
      IMPLICIT NONE
      REAL(dbl), INTENT(OUT) :: taup( :, : )
      REAL(dbl), INTENT(IN) :: tau( :, :, : )
      INTEGER, INTENT(IN) :: na( : ), nsp
      INTEGER :: is, ia, isa
      isa = 0
      DO is = 1, nsp
        DO ia = 1, na( is )
          isa = isa + 1
          taup( :, isa ) = tau( :, ia, is )
        END DO
      END DO
      RETURN
    END SUBROUTINE

    SUBROUTINE unpacktau( tau, taup, na, nsp )
      IMPLICIT NONE
      REAL(dbl), INTENT(IN) :: taup( :, : )
      REAL(dbl), INTENT(OUT) :: tau( :, :, : )
      INTEGER, INTENT(IN) :: na( : ), nsp
      INTEGER :: is, ia, isa
      isa = 0
      DO is = 1, nsp
        DO ia = 1, na( is )
          isa = isa + 1
          tau( :, ia, is ) = taup( :, isa )
        END DO
      END DO
      RETURN
    END SUBROUTINE

!------------------------------------------------------------------------------!
  END MODULE ions_base
!------------------------------------------------------------------------------!
