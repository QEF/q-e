!
! Copyright (C) 2001-2004 PWSCF group
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
  USE parameters, ONLY : DP
  USE constants,  ONLY : eps8
  !
  IMPLICIT NONE
  !
  ! ... Here the dummy variables
  !
  INTEGER :: ngm, ngk, igk(ngk)
    ! input        : the number of g vectors
    ! input/output : the number of k+G vectors inside the "ecut sphere"
    ! output       : the correspondence k+G <-> G
  REAL(KIND=DP) :: k(3), g(3,ngm), ecut, gk(ngk)
    ! input  : the k point
    ! input  : the coordinates of G vectors
    ! input  : the cut-off energy
    ! output : the moduli of k+G
  !
  ! ... here the local variables
  !
  INTEGER :: ng, nk, ngk_in
    ! counter on   G vectors
    ! counter on k+G vectors
    ! the size of vector gk (ngk is overwritten)
  REAL(KIND=DP) :: q, q2x
    ! |k+G|^2
    ! upper bound for |G|
  !
  !
  ! ... first we count the number of k+G vectors inside the cut-off sphere
  !
  q2x = ( SQRT( k(1)**2 + k(2)**2 + k(3)**2 ) + SQRT( ecut ) )**2
  !
  ngk_in = ngk
  ngk    = 0
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
        IF ( ngk > ngk_in ) &
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
        igk (ngk) = ng
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
  IF( ng > ngm ) CALL errore( 'gk_sort', 'unexpected exit from do-loop', -1 )
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
!
!
!----------------------------------------------------------------------------
SUBROUTINE gk_l2gmap( ngm, ig_l2g, ngk, igk, igk_l2g )
  !----------------------------------------------------------------------------
  !
  ! ... This subroutine maps local G+k index to the global G vector index
  ! ... the mapping is used to collect wavefunctions subsets distributed
  ! ... across processors.
  ! ... Written by Carlo Cavazzoni
  !
  USE parameters, ONLY : DP
  !
  IMPLICIT NONE
  !
  ! ... Here the dummy variables
  !
  INTEGER :: ngm, ngk, igk(ngk), ig_l2g(ngm)   ! input
  INTEGER :: igk_l2g(ngk)                      ! output
  INTEGER :: nk
  !
  ! input: mapping between local and global G vector index
  !
  !
  DO nk = 1, ngk
     !
     igk_l2g(nk) = ig_l2g( igk(nk) )
     !
  END DO
  !
  RETURN
  !
END SUBROUTINE gk_l2gmap
