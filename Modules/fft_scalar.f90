!
! Copyright (C) 2002 FPMD group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

!----------------------------------------------------------------------
! FFT scalar driver Module.
! Written by Carlo Cavazzoni 
! Last update 18 February 2002
!----------------------------------------------------------------------

!=----------------------------------------------------------------------=!
   MODULE fft_scalar
!=----------------------------------------------------------------------=!

        USE kinds

        IMPLICIT NONE
        SAVE

        PRIVATE
        PUBLIC :: initialize_tables, fft_x, fft_y, fft_z, tabmesh
        PUBLIC :: cft_1z

! ...   Local Parameter

        INTEGER, PARAMETER :: ndims = 4

#if defined __AIX

        INTEGER, PARAMETER :: nfftx = 2048
        INTEGER, PARAMETER :: lwork = 100000
        INTEGER, PARAMETER :: ltabl = 20000

#else

        INTEGER, PARAMETER :: nfftx = 1024
        INTEGER, PARAMETER :: lwork = 20 * nfftx
        INTEGER, PARAMETER :: ltabl = 4 * nfftx


#endif

#if defined __FFTW

#  if defined __TRU64 || defined __SGI64

        INTEGER*8 :: fw_plan_x(ndims) = 0
        INTEGER*8 :: fw_plan_y(ndims) = 0
        INTEGER*8 :: fw_plan_z(ndims) = 0
        INTEGER*8 :: bw_plan_x(ndims) = 0
        INTEGER*8 :: bw_plan_y(ndims) = 0
        INTEGER*8 :: bw_plan_z(ndims) = 0

#  else

        INTEGER :: fw_plan_x(ndims) = 0
        INTEGER :: fw_plan_y(ndims) = 0
        INTEGER :: fw_plan_z(ndims) = 0
        INTEGER :: bw_plan_x(ndims) = 0
        INTEGER :: bw_plan_y(ndims) = 0
        INTEGER :: bw_plan_z(ndims) = 0

#  endif

#elif defined __AIX

        REAL (dbl) :: work(lwork) 
        REAL (dbl) :: fw_tablez(ltabl,ndims)
        REAL (dbl) :: fw_tablex(ltabl,ndims)
        REAL (dbl) :: fw_tabley(ltabl,ndims)
        REAL (dbl) :: bw_tablez(ltabl,ndims)
        REAL (dbl) :: bw_tablex(ltabl,ndims)
        REAL (dbl) :: bw_tabley(ltabl,ndims)

#else

        REAL (dbl) :: work(lwork) 
        REAL (dbl) :: tablez(ltabl,ndims)
        REAL (dbl) :: tablex(ltabl,ndims)
        REAL (dbl) :: tabley(ltabl,ndims)

#endif

        REAL (dbl) :: scale
        INTEGER :: isys = 0

!=----------------------------------------------------------------------=!
      CONTAINS
!=----------------------------------------------------------------------=!


!=----------------------------------------------------------------------=!
!  initializations for the 1D scalar FFTs
!
        SUBROUTINE initialize_tables (nx, ny, nz)

          INTEGER, INTENT(IN) :: nx, ny, nz

          IF( (nx * ny * nz) == 0 ) THEN
            CALL errore(" fft_scalar: initialize_tables ", " an fft dimension is equal to zero ", 0)
          END IF
          IF( nx < 0 .OR. nx > nfftx ) THEN
            CALL errore(" fft_scalar: initialize_tables ", " nx out of range ", nx)
          END IF
          IF( ny < 0 .OR. ny > nfftx ) THEN
            CALL errore(" fft_scalar: initialize_tables ", " ny out of range ", ny)
          END IF
          IF( nz < 0 .OR. nz > nfftx ) THEN
            CALL errore(" fft_scalar: initialize_tables ", " nz out of range ", nz)
          END IF

          scale = 1.d0 / REAL(nx * ny * nz)
!
          RETURN 
        END SUBROUTINE initialize_tables

!
!=----------------------------------------------------------------------=!
!
!
!
!         FFT along "z" direction
!
!
!
!=----------------------------------------------------------------------=!
!
        SUBROUTINE fft_z(isign, c, ldc, nz, nsl)
          INTEGER, INTENT(IN) :: isign
          INTEGER, INTENT(IN) :: nsl, nz, ldc
          COMPLEX (dbl) :: c(:,:) 
          REAL(dbl)  :: tscale
          INTEGER    :: i, j
          INTEGER    :: err, idir, ip
          INTEGER, SAVE :: zdims( 3, ndims ) = -1
          INTEGER, SAVE :: icurrent = 1

          IF( nsl < 0 ) THEN
            CALL errore(" fft_scalar: fft_z ", " nsl out of range ", nsl)
          END IF

          IF( ( isign /= 0 ) .AND. ( ldc /= SIZE(c,1) ) ) THEN
            WRITE( 6, fmt = "( ' MSG: ', 2I5 )" ) SIZE(c,1), ldc
            CALL errore(" fft_scalar: fft_z ", " wrong ldc size ", ldc)
          END IF

          !
          !   Here initialize table only if necessary
          !

          ip = -1
          DO i = 1, ndims
            IF( ( nz == zdims(1,i) ) .and. ( nsl == zdims(2,i) ) .and. ( ldc == zdims(3,i) ) ) THEN
              ip = i
              EXIT
            END IF
          END DO

          IF( ip == -1 ) THEN

            WRITE(6, fmt="('DEBUG fft_z, reinitializing tables ', I3)" ) icurrent

#if defined __FFTW

            IF( fw_plan_z(icurrent) /= 0 ) CALL DESTROY_PLAN( fw_plan_z(icurrent) )
            IF( bw_plan_z(icurrent) /= 0 ) CALL DESTROY_PLAN( bw_plan_z(icurrent) )
            idir =  1; CALL CREATE_PLAN( fw_plan_z(icurrent), nz, idir) 
            idir = -1; CALL CREATE_PLAN( bw_plan_z(icurrent), nz, idir) 

#elif defined __T3E

            CALL CCFFT (0, nz, 1.0d0, c, c, tablez(1,icurrent), work(1), isys)

#elif defined __SGI

            CALL ZFFT1DI( nz, tablez(1,icurrent) )

#elif defined __AIX

            CALL DCFT ( 1, c(1,1), 1, ldc, c(1,1), 1, ldc, nz, nsl,  1, &
               scale, fw_tablez(1,icurrent), ltabl, work(1), lwork)
            CALL DCFT ( 1, c(1,1), 1, ldc, c(1,1), 1, ldc, nz, nsl, -1, &
               1.0d0, bw_tablez(1,icurrent), ltabl, work(1), lwork)

#else 

            CALL errore(' fft_z ',' no scalar fft driver specified ', 1)

#endif

            zdims(1,icurrent) = nz; zdims(2,icurrent) = nsl; zdims(3,icurrent) = ldc;
            ip = icurrent
            icurrent = MOD( icurrent, ndims ) + 1

          END IF

          !
          !   Now perform the FFTs using machine specific drivers
          !

#if defined __FFTW

          IF (isign > 0) THEN
            CALL FFT_Z_STICK(fw_plan_z(ip), c(1,1), ldc, nsl)
            CALL zdscal(SIZE(c), scale, c(1,1), 1)
          ELSE IF (isign < 0) THEN
            CALL FFT_Z_STICK(bw_plan_z(ip), c(1,1), ldc, nsl)
          END IF

#elif defined __T3E

          IF (isign /= 0) THEN
            DO i = 1, nsl
              CALL CCFFT (isign, nz, 1.0d0, c(1,i), c(1,i), tablez(1,ip), work, isys)
            END DO
            IF( isign > 0) THEN
              CALL csscal(SIZE(c), scale, c(1,1), 1)
            END IF
          END IF

#elif defined __SGI

          IF (isign /= 0) THEN
            CALL zfftm1d( isign, nz, nsl, c(1,1), 1, ldc, tablez(i,ip) )
            IF (isign > 0) THEN
              CALL zdscal(SIZE(c), scale, c(1,1), 1)
            END IF
          END IF

#elif defined __AIX

          IF( isign > 0 ) THEN
            tscale = scale
            CALL DCFT (0, c(1,1), 1, ldc, c(1,1), 1, ldc, nz, nsl, isign, &
               tscale, fw_tablez(1,ip), ltabl, work, lwork)
          ELSE IF( isign < 0 ) THEN
            tscale = 1.0d0
            CALL DCFT (0, c(1,1), 1, ldc, c(1,1), 1, ldc, nz, nsl, isign, &
               tscale, bw_tablez(1,ip), ltabl, work, lwork)
          END IF

#else 
                                                                                                      
          CALL errore(' fft_z ',' no scalar fft driver specified ', 1)

#endif

          RETURN
        END SUBROUTINE fft_z

!
!=----------------------------------------------------------------------=!
!
!
!
!         FFT along "y" direction
!
!
!
!=----------------------------------------------------------------------=!
!
        SUBROUTINE fft_y(isign, r, ldx, ldy, pl2ix, nxl, ny, nzl)

          INTEGER, INTENT(IN) :: isign, pl2ix(:), ldx, ldy, nxl, ny, nzl
          COMPLEX (dbl) :: r(:,:,:)
          COMPLEX (dbl) :: yt(ny)
          integer :: i, k, j, err, idir, ip
          INTEGER, SAVE :: icurrent = 1
          INTEGER, SAVE :: dims(2,ndims) = -1

          IF( ( isign /= 0 ) .AND. ( ldx /= SIZE(r,1) ) ) &
            CALL errore(" fft_scalar: fft_y ", " wrong ldx size ", ldx)
          IF( ( isign /= 0 ) .AND. ( ldy /= SIZE(r,2) ) ) &
            CALL errore(" fft_scalar: fft_y ", " wrong ldy size ", ldy)

          !
          !   Here initialize table only if necessary
          !

          ip = -1
          DO i = 1, ndims
            IF( ( ny == dims(1,i) ) .and. ( ldx == dims(2,i) ) ) THEN
              ip = i
              EXIT
            END IF
          END DO

          IF( ip == -1 ) THEN

#if defined __FFTW

            IF( fw_plan_y(icurrent) /= 0 )   CALL DESTROY_PLAN( fw_plan_y(icurrent) )
            IF( bw_plan_y(icurrent) /= 0 )   CALL DESTROY_PLAN( bw_plan_y(icurrent) )
            idir =  1; CALL CREATE_PLAN( fw_plan_y(icurrent), ny, idir)
            idir = -1; CALL CREATE_PLAN( bw_plan_y(icurrent), ny, idir)

#elif defined __T3E

            CALL CCFFT (0, ny, 1.0d0, yt, yt, tabley(1,icurrent), work(1), isys)

#elif defined __AIX

            CALL DCFT ( 1, r(1,1,1), ldx, 1, r(1,1,1), ldx, 1, ny, 1,  1, 1.0d0, &
               fw_tabley(1,icurrent), ltabl, work(1), lwork)
            CALL DCFT ( 1, r(1,1,1), ldx, 1, r(1,1,1), ldx, 1, ny, 1, -1, 1.0d0, &
               bw_tabley(1,icurrent), ltabl, work(1), lwork)

#elif defined __SGI

            CALL ZFFT1DI( ny, tabley(1, icurrent) )

#else

            CALL errore(' fft_y ',' no scalar fft driver specified ', 1)

#endif

            dims(1,icurrent) = ny; dims(2,icurrent) = ldx
            ip = icurrent
            icurrent = MOD( icurrent, ndims ) + 1

          END IF

          !
          !   Now perform the FFTs using machine specific drivers
          !

#if defined __FFTW

          IF( isign /= 0 ) THEN
            do i = 1, nxl
              do k = 1, nzl
                IF( pl2ix( i ) > 0 ) THEN
                  IF( isign > 0 ) THEN
                    call FFT_Y_STICK(fw_plan_y(ip), r(i,1,k), ny, ldx) 
                  ELSE
                    call FFT_Y_STICK(bw_plan_y(ip), r(i,1,k), ny, ldx) 
                  END IF
                END IF
              end do
            end do
          END IF

#elif defined __AIX

          IF( isign /= 0 ) THEN
            do i=1,nxl
              do k=1,nzl
                IF( pl2ix( i ) > 0 ) THEN
                  if(isign.gt.0) then
                    call DCFT ( 0, r(i,1,k), ldx, 1, r(i,1,k), ldx, 1, ny, 1, &
                      isign, 1.0d0, fw_tabley(1,ip), ltabl, work, lwork)
                  else
                    call DCFT ( 0, r(i,1,k), ldx, 1, r(i,1,k), ldx, 1, ny, 1, &
                      isign, 1.0d0, bw_tabley(1,ip), ltabl, work, lwork)
                  endif
                endif
              end do
            end do
          END IF

#elif defined __T3E

          IF( isign /= 0 ) THEN
            do i=1,nxl
              do k=1,nzl
                IF( pl2ix( i ) > 0 ) THEN
                  do j=1,ny
                    yt(j) = r(i,j,k)
                  end do
                  call CCFFT ( isign, ny, 1.0, yt, yt, tabley(ip), work, isys)
                  do j=1,ny
                    r(i,j,k) = yt(j)
                  end do
                END IF
              end do
            end do
          END IF

#elif defined __SGI

          IF( isign /= 0 ) THEN
            do i=1,nxl
              IF( pl2ix( i ) > 0 ) THEN
              call zfftm1d( isign, ny, nzl, r(i,1,1), ldx, ldx*ldy, tabley(ip) )
              END IF
            end do
          END IF

#else

          CALL errore(' fft_y ',' no scalar fft driver specified ', 1)

#endif

          return
        end subroutine fft_y

!
!=----------------------------------------------------------------------=!
!
!
!
!         FFT along "x" direction
!
!
!
!=----------------------------------------------------------------------=!
!

        SUBROUTINE fft_x(isign, r, ldx, ldy, nx, nyl, nzl)

          INTEGER, INTENT(IN) :: isign, ldx, ldy, nzl, nyl, nx
          COMPLEX (dbl) :: r(:,:,:)
          INTEGER :: i, j, k, err, idir, ip
          INTEGER, SAVE :: dims(1,ndims) = -1
          INTEGER, SAVE  :: icurrent = 1

          IF( ( isign /= 0 ) .AND. ( ldx /= SIZE(r,1) ) ) &
            CALL errore(" fft_scalar: fft_x ", " wrong ldx size ", ldx)
          IF( ( isign /= 0 ) .AND. ( ldy /= SIZE(r,2) ) ) &
            CALL errore(" fft_scalar: fft_x ", " wrong ldy size ", ldy)

          ip = -1
          DO i = 1, ndims
            IF( ( nx == dims(1,i) ) ) THEN
              ip = i
              EXIT
            END IF
          END DO

          IF( ip == -1 ) THEN

#if defined __FFTW

            IF( fw_plan_x(icurrent) /= 0 ) CALL DESTROY_PLAN( fw_plan_x(icurrent) )
            IF( bw_plan_x(icurrent) /= 0 ) CALL DESTROY_PLAN( bw_plan_x(icurrent) )
            idir =  1; CALL CREATE_PLAN( fw_plan_x(icurrent), nx, idir) 
            idir = -1; CALL CREATE_PLAN( bw_plan_x(icurrent), nx, idir) 

#elif defined __T3E

            CALL CCFFT (0, nx, 1.0d0, r(1,1,1), r(1,1,1), tablex(1,icurrent), work(1), isys)

#elif defined __SGI

            CALL ZFFT1DI( nx, tablex(1,icurrent) )

#elif defined __AIX

            CALL DCFT ( 1, r(1,1,1), 1, 1, r(1,1,1), 1, 1, nx, 1,  1, &
               1.0d0, fw_tablex(1,icurrent), ltabl, work(1), lwork)
            CALL DCFT ( 1, r(1,1,1), 1, 1, r(1,1,1), 1, 1, nx, 1, -1, &
               1.0d0, bw_tablex(1,icurrent), ltabl, work(1), lwork)

#else

            CALL errore(' fft_x ',' no scalar fft driver specified ', 1)

#endif

            dims(1,icurrent) = nx
            ip = icurrent
            icurrent = MOD( icurrent, ndims ) + 1

          END IF


#if defined __FFTW

          IF( isign > 0 ) THEN
            CALL FFT_X_STICK( fw_plan_x(ip), r(1,1,1), nx, nyl, nzl, ldx, ldy ) 
          ELSE IF( isign < 0 ) THEN
            CALL FFT_X_STICK( bw_plan_x(ip), r(1,1,1), nx, nyl, nzl, ldx, ldy ) 
          END IF

#elif defined __T3E

          IF( isign /= 0 ) THEN
            DO i = 1, nzl
              DO j = 1, nyl
                call CCFFT (isign, nx, 1.0, r(1,j,i), r(1,j,i), tablex(1,ip), work, isys)
              END DO
            END DO
          END IF

#elif defined __SGI

          IF( isign /= 0 ) THEN
            DO i = 1, nzl
              call zfftm1d( isign, nx, nyl, r(1,1,i), 1, ldx, tablex(1,ip) )
            END DO
          END IF

#elif defined __AIX

          IF( isign /= 0 ) THEN
            DO i = 1, nzl
              DO j = 1, nyl
                IF( isign > 0 ) THEN
                  CALL DCFT ( 0, r(1,j,i), 1, 1, r(1,j,i), 1, 1, nx, 1, isign, &
                    1.0d0, fw_tablex(1,ip), ltabl, work, lwork)
                ELSE
                  CALL DCFT ( 0, r(1,j,i), 1, 1, r(1,j,i), 1, 1, nx, 1, isign, &
                    1.0d0, bw_tablex(1,ip), ltabl, work, lwork)
                END IF
              END DO
            END DO
          END IF

#else

          CALL errore(' fft_x ',' no scalar fft driver specified ', 1)

#endif

          RETURN
        END SUBROUTINE fft_x

!=----------------------------------------------------------------------=!

!     ==================================================================
      INTEGER FUNCTION tabmesh(nr, np)
!     ==--------------------------------------------------------------==
      INTEGER, INTENT(IN) ::  nr  ! input trial radix
      INTEGER, INTENT(IN) ::  np  ! number of processors
      INTEGER  :: I
      LOGICAL  :: exist

#if defined __AIX
!     ==================================================================
!     ==   The following table list of acceptable values for          ==
!     ==   the transform lengths in the FFT is taken from pag. 758    ==
!     ==   of the ESSL manual (vol. 3)                                ==
!     ==================================================================
      INTEGER, PARAMETER :: NMX = 165
      INTEGER  LFT(NMX)
      DATA LFT /   2,   4,   6,   8,  10,  12,  14,  16,  18,  20, &
                  22,  24,  28,  30,  32,  36,  40,  42,  44,  48, &
                  56,  60,  64,  66,  70,  72,  80,  84,  88,  90, &
                  96, 110, 112, 120, 126, 128, 132, 140, 144, 154, &
                 160, 168, 176, 180, 192, 198, 210, 220, 224, 240, &
                 252, 256, 264, 280, 288, 308, 320, 330, 336, 352, &
                 360, 384, 396, 420, 440, 448, 462, 480, 504, 512, &
                 528, 560, 576, 616, 630, 640, 660, 672, 704, 720, &
                 768, 770, 792, 840, 880, 896, 924, 960, 990,1008, &
                1024,1056,1120,1152,1232,1260,1280,1320,1344,1386, &
                1408,1440,1536,1540,1584,1680,1760,1792,1848,1920, &
                1980,2016,2048,2112,2240,2304,2310,2464,2520,2560, &
                2640,2688,2772,2816,2880,3072,3080,3168,3360,3520, &
                3584,3696,3840,3960,4032,4096,4224,4480,4608,4620, &
                4928,5040,5120,5280,5376,5544,5632,5760,6144,6160, &
                6336,6720,6930,7040,7168,7392,7680,7920,8064,8192, &
                8448,8960,9216,9240,9856     /
!     ==--------------------------------------------------------------==

#else


      INTEGER, PARAMETER :: NMX = 173
      INTEGER  LFT(NMX)
      DATA LFT /   2,   3,   4,   5,   6,   8,   9,  10,  12,  15,  &
                  16,  18,  20,  24,  25,  27,  30,  32,  36,  40,  &
                  45,  48,  50,  54,  60,  64,  72,  75,  80,  81,  &
                  90,  96, 100, 108, 120, 125, 128, 135, 144, 150,  &
                 160, 162, 180, 192, 200, 216, 225, 240, 243, 250,  &
                 256, 270, 288, 300, 320, 324, 360, 375, 384, 400,  &
                 405, 432, 450, 480, 486, 500, 512, 540, 576, 600,  &
                 625, 640, 648, 675, 720, 729, 750, 768, 800, 810,  &
                 864, 900, 960, 972,1000,1024,1080,1125,1152,1200,  &
                1215,1250,1280,1296,1350,1440,1458,1500,1536,1600,  &
                1620,1728,1800,1875,1920,1944,2000,2025,2048,2160,  &
                2187,2250,2304,2400,2430,2500,2560,2592,2700,2880,  &
                2916,3000,3072,3125,3200,3240,3375,3456,3600,3645,  &
                3750,3840,3888,4000,4050,4096,4320,4374,4500,4608,  &
                4800,4860,5000,5120,5184,5400,5625,5760,5832,6000,  &
                6075,6144,6250,6400,6480,6561,6750,6912,7200,7290,  &
                7500,7680,7776,8000,8100,8192,8640,8748,9000,9216,  &
                9375,9600,9720 /

#endif

      exist = .FALSE.
      RADIX: DO I = 1, NMX
        tabmesh = lft(i)
        IF( ( tabmesh >= NR) .AND. (MOD( tabmesh, NP ) == 0) ) THEN
          exist = .TRUE.
          EXIT RADIX
        END IF
      END DO RADIX
      IF( .NOT. exist ) THEN
        CALL errore(' tabmesh ', ' radix not found ', nr )
      END IF
      RETURN
      END FUNCTION TABMESH

!
!=----------------------------------------------------------------------=!
!
!
!   Subroutine CPV
!
!
!=----------------------------------------------------------------------=!
!

!
        SUBROUTINE cft_1z(c, nsl, nz, ldc, sgn, cout)

!     driver routine for m 1d complex fft's 
!     nx=n+1 is allowed (in order to avoid memory conflicts)
!     A separate initialization is stored each combination of input sizes
!     NOTA BENE: the output in fout !

          INTEGER, INTENT(IN) :: sgn
          INTEGER, INTENT(IN) :: nsl, nz, ldc
          COMPLEX (dbl) :: c(:), cout(:) 
          REAL(dbl)  :: tscale
          INTEGER    :: i, j
          INTEGER    :: err, idir, ip, isign
          INTEGER, SAVE :: zdims( 3, ndims ) = -1
          INTEGER, SAVE :: icurrent = 1

          IF( nsl < 0 ) THEN
            CALL errore(" fft_scalar: cft_1 ", " nsl out of range ", nsl)
          END IF

          isign = -sgn

          !
          !   Here initialize table only if necessary
          !

          ip = -1
          DO i = 1, ndims
            IF( ( nz == zdims(1,i) ) .and. ( nsl == zdims(2,i) ) .and. ( ldc == zdims(3,i) ) ) THEN
              ip = i
              EXIT
            END IF
          END DO

          IF( ip == -1 ) THEN

            WRITE(6, fmt="('DEBUG fft_z, reinitializing tables ', I3)" ) icurrent

#if defined __FFTW

            IF( fw_plan_z(icurrent) /= 0 ) CALL DESTROY_PLAN( fw_plan_z(icurrent) )
            IF( bw_plan_z(icurrent) /= 0 ) CALL DESTROY_PLAN( bw_plan_z(icurrent) )
            idir = -1; CALL CREATE_PLAN( fw_plan_z(icurrent), nz, idir) 
            idir =  1; CALL CREATE_PLAN( bw_plan_z(icurrent), nz, idir) 

#elif defined __T3E

            CALL CCFFT (0, nz, 1.0d0, c, c, tablez(1,icurrent), work(1), isys)

#elif defined __SGI

            CALL ZFFT1DI( nz, tablez(1,icurrent) )

#elif defined __AIX

            tscale = 1.0d0 / nz
            CALL DCFT ( 1, c(1), 1, ldc, c(1), 1, ldc, nz, nsl,  1, &
               tscale, fw_tablez(1,icurrent), ltabl, work(1), lwork)
            CALL DCFT ( 1, c(1), 1, ldc, c(1), 1, ldc, nz, nsl, -1, &
               1.0d0, bw_tablez(1,icurrent), ltabl, work(1), lwork)

#else 

            CALL errore(' cft_1 ',' no scalar fft driver specified ', 1)

#endif

            zdims(1,icurrent) = nz; zdims(2,icurrent) = nsl; zdims(3,icurrent) = ldc;
            ip = icurrent
            icurrent = MOD( icurrent, ndims ) + 1

          END IF

          !
          !   Now perform the FFTs using machine specific drivers
          !

#if defined __FFTW

          IF (isign > 0) THEN
            tscale = 1.0d0 / nz
            CALL FFT_Z_STICK(fw_plan_z(ip), c(1), ldc, nsl)
            CALL zdscal(SIZE(c), tscale, c(1), 1)
          ELSE IF (isign < 0) THEN
            CALL FFT_Z_STICK(bw_plan_z(ip), c(1), ldc, nsl)
          END IF

#elif defined __T3E

          IF (isign /= 0) THEN
            DO i = 1, nsl
              j = (i-1) * ldc + 1
              CALL CCFFT (isign, nz, 1.0d0, c(j), c(j), tablez(1,ip), work, isys)
            END DO
            IF( isign > 0) THEN
              tscale = 1.0d0 / nz
              CALL csscal(SIZE(c), tscale, c(1), 1)
            END IF
          END IF

#elif defined __SGI

          IF (isign /= 0) THEN
            CALL zfftm1d( isign, nz, nsl, c(1), 1, ldc, tablez(i,ip) )
            IF (isign > 0) THEN
              tscale = 1.0d0 / nz
              CALL zdscal(SIZE(c), tscale, c(1), 1)
            END IF
          END IF

#elif defined __AIX

          IF( isign > 0 ) THEN
            tscale = 1.0d0 / nz
            idir   = 1
            CALL DCFT (0, c(1), 1, ldc, c(1), 1, ldc, nz, nsl, idir, &
               tscale, fw_tablez(1,ip), ltabl, work, lwork)
          ELSE IF( isign < 0 ) THEN
            idir   = -1
            tscale = 1.0d0
            CALL DCFT (0, c(1), 1, ldc, c(1), 1, ldc, nz, nsl, idir, &
               tscale, bw_tablez(1,ip), ltabl, work, lwork)
          END IF

#else 
                                                                                                      
          CALL errore(' cft_1 ',' no scalar fft driver specified ', 1)

#endif

          cout( 1 : ldc * nsl ) = c( 1 : ldc * nsl )

          RETURN
        END SUBROUTINE cft_1z


!=----------------------------------------------------------------------=!
   END MODULE fft_scalar
!=----------------------------------------------------------------------=!
