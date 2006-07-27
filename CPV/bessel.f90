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

        USE kinds,     ONLY: DP
        USE constants, ONLY: eps14

        IMPLICIT NONE
        SAVE
 
        PRIVATE

        PUBLIC :: bessel1, bessel2, bessel3

!=----------------------------------------------------------------------------=!
   CONTAINS
!=----------------------------------------------------------------------------=!

! ----------------------------------------------------------------------------
! BEGIN manual

  SUBROUTINE BESSEL1(XG, RW, JL, DJL, MMAX)

!  This subroutine Compute:
!     JL(x)  =  J_0(x);   
!     DJL(x) =  d J_0(x) / dx ; (I think there is a minus sign - PG)
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

      INTEGER :: IR, mmin

! ... Subroutine Body

      IF( ABS(XG) < eps14 ) THEN 
        CALL errore( ' bessel1',' xg too small ', 1)
      END IF

      CALL sph_bes( mmax, rw(1), xg,  0,  jl(1) )

      ! djl = - d j_0(x) /dx = + j_1(x) 

      CALL sph_bes( mmax, rw(1), xg, +1, djl(1) )

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

      REAL(DP) :: J0(MMAX)
      REAL(DP) :: J1(MMAX)
      REAL(DP) :: J2(MMAX)
      REAL(DP) :: J3(MMAX)
      INTEGER :: IR, L, LL, LMAX

! ... Subroutine Body

      IF( ABS(XG) < eps14 ) THEN 
        CALL errore( ' bessel2 ',' xg too small ', 2)
      END IF
      !
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
        CALL sph_bes( mmax, rw(1), xg,  2, j2(1) )
      END IF

      IF ( LMAX > 3 ) THEN
! ...   Calculate J3(|G||r|) = 5 * J2(|G||r|) / (|G||r|) - J1
        CALL sph_bes( mmax, rw(1), xg,  3, j3(1) )
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
      REAL(DP) :: F0(MMAX)
      REAL(DP) :: F1(MMAX)
      REAL(DP) :: F2(MMAX)
      REAL(DP) :: F3(MMAX)
      INTEGER :: IR, L, LL, LMAX, mmin

! ... Subroutine Body

      LMAX = MAXVAL( INDL ) + 1


      IF( ABS( xg * rw( 1 ) ) < eps14 ) THEN
         mmin = 2
      ELSE
         mmin = 1
      END IF

      xrg(1:mmax) = RW(1:mmax) * XG 

      IF( LMAX > 0 ) THEN
        
        ! ...   Calculate F0(|G||r|) = COS(|G||r|) 

        CALL sph_bes( (mmax-mmin+1), rw(mmin), xg, -1, F0(mmin) )
        !
        F0(mmin:mmax) = F0(mmin:mmax) * xrg(mmin:mmax)

        IF( mmin == 2 ) F0( 1 ) = F0( 2 )

      END IF

      IF( LMAX > 1 ) THEN

! ...   Calculate F1(|G||r|) = SIN(|G||r|)  = |G||r| * J0(|G||r|)

        CALL sph_bes( mmax, rw(1), xg,  0, F1(1) )

        F1 = F1 * xrg

      END IF

      IF( LMAX > 2 ) THEN

! ...   Calculate F2(|G||r|) = SIN(|G||r|) / |G||r| - COS(|G||r|)  = |G||r| * J1(|G||r|)

        F2(mmin:mmax) = (F1(mmin:mmax) / XRG(mmin:mmax) - F0(mmin:mmax))

        IF( mmin == 2 ) F2( 1 ) = F2( 2 )

      END IF

      IF( LMAX > 3 ) THEN

! ...   Calculate F3(|G||r|) = 3 F2(|G||r|)/|G||r| - F1(|G||r|) = |G||r| * J2(|G||r|)

        F3(mmin:mmax) = (3.0d0 * F2(mmin:mmax) / XRG(mmin:mmax) - F1(mmin:mmax))

        IF( mmin == 2 ) F3( 1 ) = F3( 2 )

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
