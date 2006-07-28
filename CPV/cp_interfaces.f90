!
! Copyright (C) 2002-2005 FPMD-CPV groups
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

!=----------------------------------------------------------------------------=!
   MODULE cp_interfaces
!=----------------------------------------------------------------------------=!

   IMPLICIT NONE
   PRIVATE

   PUBLIC :: bessel2
   PUBLIC :: bessel3

   INTERFACE bessel2
      SUBROUTINE bessel2(XG, RW, FINT, LNL, INDL, MMAX)
         USE kinds,     ONLY: DP
         USE constants, ONLY: eps14
         IMPLICIT NONE
         REAL(DP), INTENT(IN)  :: XG
         REAL(DP), INTENT(IN)  :: RW(:)
         REAL(DP), INTENT(OUT) :: FINT(:,:)
         INTEGER,   INTENT(IN)  :: INDL(:), LNL, MMAX
      END SUBROUTINE bessel2
   END INTERFACE

   INTERFACE bessel3
      SUBROUTINE BESSEL3(XG, RW, FINT, LNL, INDL, MMAX)
         USE kinds,     ONLY: DP
         USE constants, ONLY: eps14
         REAL(DP), INTENT(IN)  ::  XG 
         REAL(DP), INTENT(IN)  ::  RW(:)
         REAL(DP), INTENT(OUT)  ::  FINT(:,:)
         INTEGER, INTENT(IN) ::  INDL(:), LNL, MMAX
      END SUBROUTINE BESSEL3
   END INTERFACE

!=----------------------------------------------------------------------------=!
   END MODULE
!=----------------------------------------------------------------------------=!
