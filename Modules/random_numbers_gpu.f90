!
! Copyright (C) 2001-2012 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
MODULE random_numbers_gpum
  !----------------------------------------------------------------------------
  !
  USE kinds, ONLY : DP
#if defined(__CUDA)
  USE cudafor
#endif
  !
  IMPLICIT NONE
  !
  CONTAINS
    !
    !------------------------------------------------------------------------
    FUNCTION randy_gpu ( irand )
      !------------------------------------------------------------------------
      REAL(DP) :: randy_gpu
      INTEGER, optional    :: irand
#if defined(__CUDA)
      attributes(DEVICE) :: randy_gpu
#endif
      call errore('randy','use randy_vect_gpu on GPUs',1)
    END FUNCTION randy_gpu
    !------------------------------------------------------------------------
    SUBROUTINE randy_vect_gpu ( r, n, irand )
      !------------------------------------------------------------------------
      !
      ! randy_vect_gpu(r, n, irand): reseed with initial seed idum=irand ( 0 <= irand <= ic, see below)
      !                     if randyv is not explicitly initialized, it will be
      !                     initialized with seed idum=0 the first time it is called
      ! randy_vect_gpu(r, n) : generate uniform real(DP) numbers x in [0,1]
      !
      USE random_numbers, ONLY : randy
#if defined(__CUDA)
      USE curand
#endif
      REAL(DP) :: r(n)
#if defined(__CUDA)
      attributes(DEVICE) :: r
#endif
      INTEGER              :: i, n
      INTEGER, optional    :: irand
      !
      INTEGER              :: ist
      INTEGER, SAVE        :: idum=0
#if defined(__CUDA)
      attributes(DEVICE)          :: istat
      type(curandGenerator), SAVE :: gen

      LOGICAL, SAVE        :: first=.true.
      !
      IF ( present(irand) ) THEN
         idum = MIN( ABS(irand), idum) 
         first=.true.
      END IF
      !
      IF ( first ) THEN
         !
         first = .false.
         ist=curandDestroyGenerator(gen)
         ist=curandCreateGenerator(gen, CURAND_RNG_PSEUDO_XORWOW) 
         ist=curandSetPseudoRandomGeneratorSeed(gen, idum)
         !
      END IF
      !
      ist=curandGenerateUniformDouble(gen,r,n)
      !
#else
      ! randy_vect_gpu is not a GPU array in this case
      !
      ! ist means starting index here
      ist = 1
      IF ( present(irand) ) THEN
         r(1) = randy(irand)
         ist = 2
      END IF
      DO i = ist, n
         r(i) = randy()
      END DO
#endif
      RETURN
      !
    END SUBROUTINE randy_vect_gpu
    !
    !------------------------------------------------------------------------
    SUBROUTINE randy_vect_debug_gpu (r, n, irand )
      !------------------------------------------------------------------------
      !
      ! randy_vect_debug_gpu(r, n, irand): reseed with initial seed idum=irand ( 0 <= irand <= ic, see below)
      !                           if randyv is not explicitly initialized, it will be
      !                           initialized with seed idum=0 the first time it is called
      ! randy_vect_debug_gpu(r, n) : generate uniform real(DP) numbers x in [0,1]
      !
      USE random_numbers, ONLY : randy
      !
      REAL(DP) :: r(n)
      INTEGER, optional    :: irand
#if defined(__CUDA)
      attributes(DEVICE) :: r
#endif
      INTEGER :: n, i, ist
      REAL(DP), ALLOCATABLE :: aux_v(:)
      !
      print *, 'Allocation of ', n
      ALLOCATE(aux_v(n))
      !
      ist = 1
      IF ( present(irand) ) THEN
         aux_v(1) = randy(irand)
         ist = 2
      END IF
      !
      DO i = ist, n
         aux_v(i) = randy()
      END DO
      !
      r(1:n) = aux_v(1:n)
      !
      DEALLOCATE(aux_v)
    END SUBROUTINE randy_vect_debug_gpu
    !
    !------------------------------------------------------------------------
    SUBROUTINE set_random_seed ( )
      !------------------------------------------------------------------------
      !
      ! poor-man random seed for randy
      !
      INTEGER, DIMENSION (8) :: itime
      INTEGER  :: iseed
      REAL(DP) :: drand(1)
#if defined(__CUDA)
      attributes(DEVICE) :: drand
#endif
      !
      CALL date_and_time ( values = itime ) 
      ! itime contains: year, month, day, time difference in minutes, hours,
      !                 minutes, seconds and milliseconds. 
      iseed = ( itime(8) + itime(6) ) * ( itime(7) + itime(4) )
      CALL randy_vect_gpu ( drand, 1, iseed )
      CALL randy_vect_debug_gpu (drand, 1, iseed )
      !
    END SUBROUTINE set_random_seed
    !
END MODULE random_numbers_gpum
