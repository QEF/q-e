!
! Copyright (C) 2001-2004 Quantum-ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
SUBROUTINE gk_sort( k, ngm, g, ecut, ngk, igk, gk )
  !----------------------------------------------------------------------------
  !
  ! ... sorts k+g in order of increasing magnitude, up to ecut
  ! ... NB: this version will yield the same ordering for different ecut
  ! ...     and the same ordering in all machines
  !
  USE kinds,     ONLY : DP
  USE constants, ONLY : eps8
  USE wvfct,     ONLY : npwx
  !
  IMPLICIT NONE
  !
  ! ... Here the dummy variables
  !
  INTEGER, INTENT(IN) :: ngm
    ! input        : the number of g vectors
  INTEGER, INTENT(INOUT) :: ngk
    ! input/output : the number of k+G vectors inside the "ecut sphere"
  INTEGER, INTENT(OUT) :: igk(npwx)  
    ! output       : the correspondence k+G <-> G

  REAL(DP), INTENT(IN) :: k(3), g(3,ngm), ecut
    ! input  : the k point
    ! input  : the coordinates of G vectors
    ! input  : the cut-off energy
  REAL(DP), INTENT(OUT) :: gk(npwx)
    ! output : the moduli of k+G
  !
  INTEGER :: ng, nk
    ! counter on   G vectors
    ! counter on k+G vectors
  REAL(DP) :: q, q2x
    ! |k+G|^2
    ! upper bound for |G|
  !
  !
  ! ... first we count the number of k+G vectors inside the cut-off sphere
  !
  q2x = ( SQRT( k(1)**2 + k(2)**2 + k(3)**2 ) + SQRT( ecut ) )**2
  !
  ngk = 0
  !
  DO ng = 1, ngm
     !
     q = ( k(1) + g(1,ng) )**2 + ( k(2) + g(2,ng) )**2 + ( k(3) + g(3,ng) )**2
     !
     ! ... here if |k+G|^2 <= Ecut
     !
     IF ( q <= ecut ) THEN
        !
        ngk = ngk + 1
        !
        ! ... gk is a fake quantity giving the same ordering on all machines
        !
        IF ( ngk > npwx ) &
           CALL errore( 'gk_sort', 'array gk out-of-bounds', 1 )
        !
        IF ( q > eps8 ) THEN 
           !
           gk(ngk) = q 
           !
        ELSE
           !
           gk(ngk) = 0.D0
           !
        END IF
        !
        ! ... set the initial value of index array
        !
        igk(ngk) = ng
        !
     ELSE
        !
        ! ... if |G| > |k| + SQRT( Ecut )  stop search and order vectors
        !
        IF ( ( g(1,ng)**2 + g(2,ng)**2 + g(3,ng)**2 ) > ( q2x + eps8 ) ) EXIT
        !
     END IF
     !
  END DO
  !
  IF ( ng > ngm ) &
     CALL infomsg( 'gk_sort', 'unexpected exit from do-loop')
  !
  ! ... order vector gk keeping initial position in index
  !
  CALL hpsort_eps( ngk, gk, igk, eps8 )
  !
  ! ... now order true |k+G|
  !
  DO nk = 1, ngk
     !
     gk(nk) = ( k(1) + g(1,igk(nk) ) )**2 + &
              ( k(2) + g(2,igk(nk) ) )**2 + &
              ( k(3) + g(3,igk(nk) ) )**2
     !         
  END DO
  !
  RETURN
  !
END SUBROUTINE gk_sort
