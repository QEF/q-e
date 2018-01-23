!
! Copyright (C) 2011 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!=----------------------------------------------------------------------=
MODULE fft_ggen
!=----------------------------------------------------------------------=

  !  ... subroutines generating variables nl* needed to map G-vector
  !  ... components onto the FFT grid(s) in reciprocal space
  !
   USE fft_param

   PRIVATE
   SAVE

   PUBLIC :: fft_set_nl, fft_set_nlm

!=----------------------------------------------------------------------=
CONTAINS
!=----------------------------------------------------------------------=
!
!-----------------------------------------------------------------------
   SUBROUTINE fft_set_nl ( dfft, at, g, mill  )
!----------------------------------------------------------------------
   !
   ! Input:  FFT descriptor dfft, lattice vectors at, list of G-vectors g
   ! Output: indices nl such that G_fft(nl(i)) = G(i)
   !         optionally, Miller indices: if bg = reciprocal lattice vectors,
   ! G(:,i) = mill(1,i)*bg(:,1) + mill(2,i)*bg(:,2) + mill(3,i)*bg(:,3)
   !  
   USE fft_types,  ONLY : fft_type_descriptor
   !
   IMPLICIT NONE
   !
   TYPE (fft_type_descriptor), INTENT(inout) :: dfft
   REAL(DP), INTENT(IN) :: g(:,:)
   REAL(DP), INTENT(IN) :: at(:,:)
   INTEGER, OPTIONAL, INTENT(OUT) :: mill(:,:)
   INTEGER :: ng, n1, n2, n3
   !
   IF( ALLOCATED( dfft%nl ) ) DEALLOCATE( dfft%nl )
   ALLOCATE( dfft%nl( dfft%ngm ) )
   !
   DO ng = 1, dfft%ngm
      n1 = nint (sum(g (:, ng) * at (:, 1)))
      IF(PRESENT(mill)) mill (1,ng) = n1
      IF (n1<0) n1 = n1 + dfft%nr1

      n2 = nint (sum(g (:, ng) * at (:, 2)))
      IF(PRESENT(mill)) mill (2,ng) = n2
      IF (n2<0) n2 = n2 + dfft%nr2

      n3 = nint (sum(g (:, ng) * at (:, 3)))
      IF(PRESENT(mill)) mill (3,ng) = n3
      IF (n3<0) n3 = n3 + dfft%nr3

      IF (n1>=dfft%nr1 .or. n2>=dfft%nr2 .or. n3>=dfft%nr3) &
         CALL fftx_error__('ggen','Mesh too small?',ng)

      IF ( dfft%lpara) THEN
         dfft%nl (ng) = 1 + n3 + ( dfft%isind ( 1 + n1 + n2*dfft%nr1x) - 1) * dfft%nr3x
      ELSE
         dfft%nl (ng) = 1 + n1 + n2 * dfft%nr1x + n3 * dfft%nr1x * dfft%nr2x
      ENDIF
   ENDDO
   !

   END SUBROUTINE fft_set_nl 
   !
   !
   !-----------------------------------------------------------------------
   SUBROUTINE fft_set_nlm( dfft, mill )
   !----------------------------------------------------------------------
   !
   ! Input:  FFT descriptor dfft, miller indices
   ! Output: indices nlm such that G_fft(nlm(i)) = -G(i)
   !
   USE fft_types,  ONLY : fft_type_descriptor
   !
   IMPLICIT NONE
   !
   TYPE (fft_type_descriptor), INTENT(inout) :: dfft
   INTEGER, INTENT(IN) :: mill(:,:)
   !
   INTEGER :: n1, n2, n3, ng
   !
   IF( ALLOCATED( dfft%nlm ) ) DEALLOCATE( dfft%nlm )
   ALLOCATE( dfft%nlm( dfft%ngm ) )
   !
   DO ng = 1, dfft%ngm
      n1 = -mill (1,ng)
      IF (n1 < 0) THEN
         n1 = n1 + dfft%nr1
      END IF
      n2 = -mill (2,ng)
      IF (n2 < 0) THEN
         n2 = n2 + dfft%nr2
      END IF
      n3 = -mill (3,ng)
      IF (n3 < 0) THEN
         n3 = n3 + dfft%nr3
      END IF

      IF (n1>=dfft%nr1 .or. n2>=dfft%nr2 .or. n3>=dfft%nr3) THEN
         CALL fftx_error__('index_minusg','Mesh too small?',ng)
      ENDIF

      IF ( dfft%lpara ) THEN
         dfft%nlm(ng) = 1 + n3 + ( dfft%isind ( 1 + n1 + n2*dfft%nr1x) - 1) * dfft%nr3x
      ELSE
         dfft%nlm(ng) = 1 + n1 + n2 * dfft%nr1x + n3 * dfft%nr1x * dfft%nr2x
      ENDIF
   ENDDO

   END SUBROUTINE fft_set_nlm
   !
!=----------------------------------------------------------------------=
   END MODULE fft_ggen
!=----------------------------------------------------------------------=
