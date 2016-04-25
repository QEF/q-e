!
! Copyright (C) 2001-2004 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
SUBROUTINE print_clock_pwcond()
   !---------------------------------------------------------------------------
   !
   ! ... this routine prints out the clocks at the end of the pwcond run
   ! ... it tries to construct the calling tree of the program.
   !
   USE io_global,     ONLY : stdout
   USE mp_world,      ONLY : mpime, root
   USE cond,          ONLY : ikind
   !
   IMPLICIT NONE
   !
   !
   WRITE( stdout, * )
   !
   CALL print_clock( 'PWCOND' )
   CALL print_clock( 'init' )
   CALL print_clock( 'poten' )
   CALL print_clock( 'local' )
   !
   !
   WRITE( stdout, * )
   !
   CALL print_clock( 'scatter_forw' )
   CALL print_clock( 'integrals' )
   CALL print_clock( 'scatter' )
   CALL print_clock( 'rotatef' )
   CALL print_clock( 'rotateb' )
   CALL print_clock( 'scatter_back' )
   !
   WRITE( stdout, * )

   CALL print_clock( 'compbs' )
   CALL print_clock( 'compbs_2' )
   !
   WRITE( stdout, * )

   if (ikind.gt.0) then
      CALL print_clock( 'transmit' )
      CALL print_clock( 'set_ls' )
      CALL print_clock( 'solve_ls' )
   endif
   !
   RETURN
   !
END SUBROUTINE print_clock_pwcond
