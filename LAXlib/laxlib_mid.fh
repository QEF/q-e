!
! Copyright (C) 2003-2013 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!

INTERFACE laxlib_pdsyevd
SUBROUTINE laxlib_pdsyevd_x( tv, n, idesc, hh, ldh, e )
         IMPLICIT NONE
         include 'laxlib_param.fh'
         include 'laxlib_kinds.fh'
         LOGICAL, INTENT(IN) :: tv
         INTEGER, INTENT(IN) :: n, ldh
         INTEGER, INTENT(IN) :: idesc(LAX_DESC_SIZE)
         REAL(DP) :: hh( ldh, ldh )
         REAL(DP) :: e( n )
END SUBROUTINE
END INTERFACE

INTERFACE laxlib_pzheevd
SUBROUTINE laxlib_pzheevd_x( tv, n, idesc, hh, ldh, e )
   IMPLICIT NONE
   include 'laxlib_param.fh'
   include 'laxlib_kinds.fh'
   LOGICAL, INTENT(IN) :: tv
   INTEGER, INTENT(IN) :: n, ldh
   INTEGER, INTENT(IN) :: idesc(LAX_DESC_SIZE)
   COMPLEX(DP) :: hh( ldh, ldh )
   REAL(DP) :: e( n )
END SUBROUTINE
END INTERFACE

INTERFACE laxlib_pzpotrf
SUBROUTINE laxlib_pzpotrf_x( sll, ldx, n, idesc )
   implicit none
   include 'laxlib_param.fh'
   include 'laxlib_kinds.fh'
   integer :: n, ldx
   integer, INTENT(IN) :: idesc(LAX_DESC_SIZE)
   complex(DP) :: sll( ldx, ldx )
END SUBROUTINE
END INTERFACE

INTERFACE laxlib_pdpotrf
SUBROUTINE laxlib_pdpotrf_x( sll, ldx, n, idesc )
   implicit none
   include 'laxlib_param.fh'
   include 'laxlib_kinds.fh'
   integer  :: n, ldx
   INTEGER, INTENT(IN) :: idesc(LAX_DESC_SIZE)
   REAL(DP) :: sll( ldx, ldx )
END SUBROUTINE
END INTERFACE

INTERFACE laxlib_pztrtri
SUBROUTINE laxlib_pztrtri_x ( sll, ldx, n, idesc )
   implicit none
   include 'laxlib_param.fh'
   include 'laxlib_kinds.fh'
   INTEGER, INTENT( IN ) :: n, ldx
   INTEGER, INTENT(IN) :: idesc(LAX_DESC_SIZE)
   COMPLEX(DP), INTENT( INOUT ) :: sll( ldx, ldx )
END SUBROUTINE
END INTERFACE

INTERFACE laxlib_pdtrtri
SUBROUTINE laxlib_pdtrtri_x ( sll, ldx, n, idesc )
   implicit none
   include 'laxlib_param.fh'
   include 'laxlib_kinds.fh'
   INTEGER, INTENT( IN ) :: n, ldx
   INTEGER, INTENT(IN) :: idesc(LAX_DESC_SIZE)
   REAL(DP), INTENT( INOUT ) :: sll( ldx, ldx )
END SUBROUTINE
END INTERFACE
