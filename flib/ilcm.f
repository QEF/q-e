      INTEGER FUNCTION ILCM( M, N )
!
!  -- ScaLAPACK tools routine (version 1.5) --
!     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
!     and University of California, Berkeley.
!     May 1, 1997
!
!     .. Scalar Arguments ..
      INTEGER            M, N
!     ..
!
!  Purpose
!  =======
!
!  ILCM computes and returns the Least Common Multiple (LCM) of two
!  positive integers M and N. In fact the routine computes the greatest
!  common divisor (GCD) and use the fact that M*N = GCD*LCM.
!
!  Arguments
!  =========
!
!  M       (input) INTEGER
!          On entry, M >=0. Unchanged on exit.
!
!  N       (input) INTEGER
!          On entry, N >=0. Unchanged on exit.
!
!  =====================================================================
!
!     .. Local Scalars ..
      INTEGER            IA, IQ, IR
!     ..
!     .. Executable Statements ..
!
      IF( M.GE.N ) THEN
         IA = M
         ILCM = N
      ELSE
         IA = N
         ILCM = M
      ENDIF
!
   10 CONTINUE
         IQ = IA / ILCM
         IR = IA - IQ * ILCM
         IF( IR.EQ.0 ) THEN
            ILCM = ( M * N ) / ILCM
            RETURN
         END IF
         IA = ILCM
         ILCM = IR
      GO TO 10
!
!     End of ILCM
!
      END
