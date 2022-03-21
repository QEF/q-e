!
! Copyright (C) 2001-2010 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
SUBROUTINE gk_sort( k, ngm, g, ecut, ngk, igk, gk )
   !----------------------------------------------------------------------------
   !! Sorts k+g in order of increasing magnitude, up to ecut.
   !
   !! NB1: this version should yield the same ordering for different ecut
   !!      and the same ordering in all machines AS LONG AS INPUT DATA
   !!      IS EXACTLY THE SAME
   !! NB2: this version assumes that input G-vectors are ordered by the same
   !!      routine, "hpsort_eps", used here. In principle for k=0 this should
   !!      guarantee that the ordering of k+G is the same as for G-vectors.
   !!      In practice, in some special cases (primitive lattice vectors that
   !!      are close to but not exactly equal to a symmetric Bravais lattice)
   !!      this does not hold, presumably due to a limitation of "hpsort_eps".
   !!      This is a source of trouble for Gamma-only calculations, so here
   !!      we explicitly set igk(i)=i for k=0.
   !
   USE kinds,      ONLY: DP
   USE constants,  ONLY: eps8
   USE wvfct,      ONLY: npwx
   !
   IMPLICIT NONE
   !
   REAL(DP), INTENT(IN)  :: k(3)
   !! the k point
   INTEGER,  INTENT(IN)  :: ngm
   !! the number of g vectors
   REAL(DP), INTENT(IN)  :: g(3,ngm)
   !! the coordinates of G vectors
   REAL(DP), INTENT(IN)  :: ecut
   !! the cut-off energy
   INTEGER,  INTENT(OUT) :: ngk
   !! the number of k+G vectors inside the "ecut sphere"
   INTEGER,  INTENT(OUT) :: igk(npwx)
   !! the correspondence k+G <-> G
   REAL(DP), INTENT(OUT) :: gk(npwx)
   !! the moduli of k+G
   !
   !  ... local variables
   !
   INTEGER :: ng   ! counter on   G vectors
   INTEGER :: nk   ! counter on k+G vectors
   REAL(DP) :: q   ! |k+G|^2
   REAL(DP) :: q2x ! upper bound for |G|
   !
   ! ... first we count the number of k+G vectors inside the cut-off sphere
   !
   q2x = ( SQRT( SUM(k(:)**2) ) + SQRT( ecut ) )**2
   !
   ngk = 0
   igk(:) = 0
   gk (:) = 0.0_DP
   !
   DO ng = 1, ngm
      q = SUM( ( k(:) + g(:,ng) )**2 )
      IF ( q <= eps8 ) q = 0.0_DP
      !
      ! ... here if |k+G|^2 <= Ecut
      !
      IF ( q <= ecut ) THEN
         ngk = ngk + 1
         IF ( ngk > npwx ) &
            CALL errore( 'gk_sort', 'array gk out-of-bounds', 1 )
         !
         gk(ngk) = q
         !
         ! set the initial value of index array
         igk(ngk) = ng
      ELSE
         ! if |G| > |k| + SQRT( Ecut )  stop search and order vectors
         IF ( SUM( g(:,ng)**2 ) > ( q2x + eps8 ) ) EXIT
      ENDIF
   ENDDO
   !
   IF ( ng > ngm ) &
      CALL infomsg( 'gk_sort', 'unexpected exit from do-loop' )
   !
   ! ... order vector gk keeping initial position in index
   ! ... see comments above about the k=0 case
   !
   IF ( k(1)**2 + k(2)**2 + k(3)**2 > eps8 ) THEN
      CALL hpsort_eps( ngk, gk, igk, eps8 )
      !
      ! ... now order true |k+G|
      !
      DO nk = 1, ngk
         gk(nk) = SUM( (k(:) + g(:,igk(nk)) )**2 )
      ENDDO
   END IF
   !
END SUBROUTINE gk_sort
