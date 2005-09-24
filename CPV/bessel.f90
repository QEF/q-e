!
! Copyright (C) 2002-2005 FPMD-CPV groups
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

!=----------------------------------------------------------------------------=!
   MODULE bessel_functions
!=----------------------------------------------------------------------------=!

        USE kinds

        IMPLICIT NONE
        SAVE
 
        PRIVATE

        REAL(DP) :: small = 1.0d-14

        PUBLIC :: bessel1, bessel2, bessel3

!=----------------------------------------------------------------------------=!
   CONTAINS
!=----------------------------------------------------------------------------=!

! ----------------------------------------------------------------------------
! BEGIN manual

  SUBROUTINE BESSEL1(XG, RW, JL, DJL, MMAX)

!  This subroutine Compute:
!     JL(x)  =  J_0(x);   
!     DJL(x) =  d J_0(x) / dx ;   
!     x = XG * RW(i);  i = 1, ..., mmax

! END manual
! ----------------------------------------------------------------------------

      IMPLICIT NONE

! ... Argument Variables

      REAL(DP), INTENT(IN)  :: XG
      REAL(DP), INTENT(IN)  :: RW(:)
      REAL(DP), INTENT(OUT) :: JL(:)
      REAL(DP), INTENT(OUT) :: DJL(:)
      INTEGER, INTENT(IN) :: MMAX

! ... Local Variables

      REAL(DP) :: ARG_S(MMAX)
      REAL(DP) :: XRGM1(MMAX)
      INTEGER :: IR

! ... Subroutine Body

      IF( ABS(XG) < small ) THEN 
        CALL errore( ' bessel1',' xg too small ', 1)
      END IF

      DO ir = 1, mmax
        arg_s( ir ) = rw( ir ) * xg 
      END DO
      DO ir = 1, mmax
        xrgm1( ir ) = 1.0d0 / arg_s( ir )
      END DO

      CALL sph_bes( mmax, rw(1), xg,  0, jl(1) )
      CALL sph_bes( mmax, rw(1), xg, -1, djl(1) )
      djl(1:mmax) = jl(1:mmax) * xrgm1(1:mmax) - djl(1:mmax)

      RETURN

  END SUBROUTINE bessel1


! ----------------------------------------------------------------------------
! BEGIN manual

  SUBROUTINE bessel2(XG, RW, FINT, LNL, INDL, MMAX)

!  This subroutine Compute:
!     Fint(x,l) = J_l(x);   l = INDL(j); j = 1, LNL
!     x = XG * RW(i);  i = 1, ..., mmax

! END manual
! ----------------------------------------------------------------------------

      IMPLICIT NONE

! ... Argument Variables

      REAL(DP), INTENT(IN)  :: XG
      REAL(DP), INTENT(IN)  :: RW(:)
      REAL(DP), INTENT(OUT) :: FINT(:,:)
      INTEGER,   INTENT(IN)  :: INDL(:), LNL, MMAX

! ... Local Variables

      REAL(DP) :: ARG_S(MMAX)
      REAL(DP) :: XRGM1(MMAX)
      REAL(DP) :: J0(MMAX)
      REAL(DP) :: J1(MMAX)
      REAL(DP) :: J2(MMAX)
      REAL(DP) :: J3(MMAX)
      INTEGER :: IR, L, LL, LMAX

! ... Subroutine Body

      IF( ABS(XG) < small ) THEN 
        CALL errore( ' bessel2 ',' xg too small ', 2)
      END IF

      ARG_S(1:mmax) = RW(1:mmax) * XG 
      xrgm1 = 1.0d0 / arg_s

      LMAX = MAXVAL( INDL ) + 1

      IF ( LMAX > 0 ) THEN
! ...   Calculate J0(|G||r|) = SIN(|G||r|) / (|G||r|)
        CALL sph_bes( mmax, rw(1), xg,  0, j0(1) )
      END IF

      IF ( LMAX > 1 ) THEN
! ...   Calculate J1(|G||r|) = SIN(|G||r|) / (|G||r|)^2 - COS(|G||r|) / (|G||r|)
        CALL sph_bes( mmax, rw(1), xg,  1, j1(1) )
      END IF

      IF ( LMAX > 2 ) THEN
! ...   Calculate J2(|G||r|) = 3 * J1(|G||r|) / (|G||r|) - J0
        J2(1:mmax) = (J1(1:mmax) * 3.0d0 * XRGM1(1:mmax) - J0(1:mmax))
      END IF

      IF ( LMAX > 3 ) THEN
! ...   Calculate J3(|G||r|) = 5 * J2(|G||r|) / (|G||r|) - J1
        J3(1:mmax) = (J2(1:mmax) * 5.0d0 * XRGM1(1:mmax) - J1(1:mmax))
      END IF

      DO L = 1,LNL
        LL = INDL(L)
        IF(LL == 0) THEN 
! ...     FINT = FUNT * J0
          FINT(1:mmax,L) = J0(1:mmax)
        ELSE IF (LL == 1) THEN
! ...     FINT = FUNT * J1
          FINT(1:mmax,L) = J1(1:mmax)
        ELSE IF (LL == 2) THEN
! ...     FINT = FUNT * J2
          FINT(1:mmax,L) = J2(1:mmax)
        ELSE IF (LL == 3) THEN
! ...     FINT = FUNT * J3
          FINT(1:mmax,L) = J3(1:mmax)
        ELSE
          CALL errore(" bessel2 "," ll value not programmed ", MAX( 1, ABS(ll) ) )
        END IF
      END DO

      RETURN
      END SUBROUTINE bessel2

! ----------------------------------------------------------------------------
! BEGIN manual

  SUBROUTINE BESSEL3(XG, RW, FINT, LNL, INDL, MMAX)

!  This subroutine Compute:
!     Fint(x,l) = f_l(x);   l = INDL(j); j = 1, LNL
!     x = XG * RW(i);  i = 1, ..., mmax
!     f_0(x) = cos(x) 
!     f_l(x) = x * j_(l-1)(x); l > 0
!
! END manual
! ----------------------------------------------------------------------------

      IMPLICIT NONE

! ... Argument Variables

      REAL(DP), INTENT(IN)  ::  XG
      REAL(DP), INTENT(IN)  ::  RW(:)
      REAL(DP), INTENT(OUT)  ::  FINT(:,:)
      INTEGER, INTENT(IN) ::  INDL(:), LNL, MMAX

! ... Local Variables

      REAL(DP) :: XRG(MMAX)
      REAL(DP) :: XRGM1(MMAX)
      REAL(DP) :: F0(MMAX)
      REAL(DP) :: F1(MMAX)
      REAL(DP) :: F2(MMAX)
      REAL(DP) :: F3(MMAX)
      INTEGER :: IR, L, LL, LMAX

! ... Subroutine Body

      LMAX = MAXVAL( INDL ) + 1

      IF( ABS(XG) < small ) THEN 
        CALL errore( ' bessel3 ',' xg too small ', 1)
      END IF

      xrg(1:mmax) = RW(1:mmax) * XG 

      IF( LMAX > 0 ) THEN
        
! ...   Calculate F0(|G||r|) = COS(|G||r|) 

        CALL sph_bes( mmax, rw(1), xg, -1, F0(1) )
        F0 = F0 * xrg

      END IF

      IF( LMAX > 1 ) THEN

! ...   Calculate F1(|G||r|) = SIN(|G||r|)  = |G||r| * J0(|G||r|)

        CALL sph_bes( mmax, rw(1), xg,  0, F1(1) )
        F1 = F1 * xrg

      END IF

      IF( LMAX > 2 ) THEN

! ...   Calculate F2(|G||r|) = SIN(|G||r|) / |G||r| - COS(|G||r|)  = |G||r| * J1(|G||r|)

        XRGM1(1:mmax) = 1.0d0/xrg(1:mmax)

        F2(1:mmax) = (F1(1:mmax) * XRGM1(1:mmax) - F0(1:mmax))

      END IF

      IF( LMAX > 3 ) THEN

! ...   Calculate F3(|G||r|) = 3 F2(|G||r|)/|G||r| - F1(|G||r|) = |G||r| * J2(|G||r|)

        F3(1:mmax) = (3.0d0 * F2(1:mmax) * XRGM1(1:mmax) - F1(1:mmax))

      END IF

      DO L = 1,LNL
        LL = INDL(L)
        IF(LL.EQ.0) THEN 
          FINT(1:mmax, L) = F0(1:mmax)
        ELSE IF (LL.EQ.1) THEN
          FINT(1:mmax, L) = F1(1:mmax)
        ELSE IF (LL.EQ.2) THEN
          FINT(1:mmax, L) = F2(1:mmax)
        ELSE IF (LL.EQ.3) THEN
          FINT(1:mmax, L) = F3(1:mmax)
        ELSE
          CALL errore(" bessel3 "," ll value not programmed ", MAX( 1, ABS(ll) ) )
        END IF
      END DO

      RETURN
      END SUBROUTINE bessel3


!=----------------------------------------------------------------------------=!
   END MODULE bessel_functions
!=----------------------------------------------------------------------------=!
