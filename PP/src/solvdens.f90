!
! Copyright (C) 2016 National Institute of Advanced Industrial Science and Technology (AIST)
! [ This code is written by Satomichi Nishihara. ]
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------------------------------------------
SUBROUTINE solvdens(filplot, lpunch)
  !--------------------------------------------------------------------------
  !
  ! ... Writes the solvent density and potential
  ! ... into a file format suitable for plotting
  !
  USE cell_base,      ONLY : at, bg, celldm, ibrav, alat, omega, tpiba, tpiba2
  USE chdens_module,  ONLY : plot_1d, plot_2d, plot_3d, plot_2ds, plot_fast
  USE constants,      ONLY : eps6, tpi
  USE control_flags,  ONLY : gamma_only
  USE fft_base,       ONLY : dfftp, dffts
  USE fft_types,      ONLY : fft_type_allocate
  USE fft_interfaces, ONLY : fwfft
  USE gvect,          ONLY : ngm, g, gcutm, gg, mill, ig_l2g, ngm_g, gstart
  USE gvecs,          ONLY : gcutms, doublegrid, dual, ngms
  USE gvecw,          ONLY : ecutwfc
  USE io_files,       ONLY : prefix
  USE io_global,      ONLY : stdin, stdout, ionode, ionode_id
  USE ions_base,      ONLY : nat, ityp, atm, ntyp => nsp, tau, zv
  USE kinds,          ONLY : DP
  USE lsda_mod,       ONLY : nspin
  USE mp,             ONLY : mp_bcast
  USE mp_bands,       ONLY : intra_bgrp_comm, nyfft
  USE mp_images,      ONLY : intra_image_comm
  USE recvec_subs,    ONLY : ggen, ggens
  USE run_info,       ONLY : title
  USE scatter_mod,    ONLY : scatter_grid
  !
  IMPLICIT NONE
  !
  CHARACTER(LEN=*), INTENT(IN) :: filplot
  LOGICAL,          INTENT(IN) :: lpunch
  !
  LOGICAL                        :: laue
  LOGICAL                        :: fast3d
  LOGICAL                        :: avoid_fft
  INTEGER                        :: ios
  INTEGER                        :: ounit
  INTEGER                        :: ounit1
  INTEGER                        :: iflag
  INTEGER                        :: nx, ny, nz
  INTEGER                        :: output_format
  INTEGER                        :: isite
  INTEGER                        :: nsite
  INTEGER                        :: ifft
  INTEGER                        :: nfft
  INTEGER                        :: na, ipol
  INTEGER                        :: lebedev
  REAL(DP)                       :: e1(3), e2(3), e3(3)
  REAL(DP)                       :: x0(3)
  REAL(DP)                       :: radius
  REAL(DP)                       :: m1, m2, m3
  REAL(DP)                       :: rhotot
  REAL(DP),          ALLOCATABLE :: rhor(:,:)
  REAL(DP),          ALLOCATABLE :: rhos(:)
  REAL(DP),          ALLOCATABLE :: raux(:)
  COMPLEX(DP),       ALLOCATABLE :: rhog(:)
  COMPLEX(DP),       ALLOCATABLE :: caux(:)
  CHARACTER(LEN=12), ALLOCATABLE :: asite(:)
  CHARACTER(LEN=16)              :: fileext
  CHARACTER(LEN=256)             :: fileout
  CHARACTER(LEN=256)             :: fileout0
  CHARACTER(LEN=20)              :: interpolation
  !
  CHARACTER(LEN=75)              :: title_
  INTEGER                        :: nsite_
  INTEGER                        :: nr1x_
  INTEGER                        :: nr2x_
  INTEGER                        :: nr3x_
  INTEGER                        :: nr1_
  INTEGER                        :: nr2_
  INTEGER                        :: nr3_
  INTEGER                        :: nat_
  INTEGER                        :: ntyp_
  INTEGER                        :: ibrav_
  REAL(DP)                       :: celldm_(6)
  REAL(DP)                       :: at_(3, 3)
  REAL(DP)                       :: gcutm_
  LOGICAL                        :: laue_
  !
  INTEGER,           EXTERNAL    :: find_free_unit
  !
  CHARACTER(LEN=7),  PARAMETER   :: NAME_FOURIER = 'fourier'
  CHARACTER(LEN=7),  PARAMETER   :: NAME_BSPLINE = 'bspline'
  !
  CHARACTER(LEN=13), DIMENSION(0:7) :: formatname = (/ &
        & 'gnuplot      ', &
        & 'obsolete!    ', &
        & 'plotrho.x    ', &
        & 'XCrySDen     ', &
        & 'obsolete!    ', &
        & 'XCrySDen     ', &
        & 'Gaussian cube', &
        & 'gnuplot x,y,f' /)
  !
  CHARACTER(LEN=5), DIMENSION(0:7) :: formatext = (/ &
        & '.gnu ', &
        & '     ', &
        & '     ', &
        & '.xsf ', &
        & '     ', &
        & '.xsf ', &
        & '.cube', &
        & '.gnu ' /)
  !
  CHARACTER(LEN=20), DIMENSION(0:4) :: plotname = (/ &
        & '1D spherical average', &
        & '1D along a line     ', &
        & '2D contour          ', &
        & '3D                  ', &
        & '2D polar on a sphere' /)
  !
  NAMELIST / plot / iflag, output_format, fileout, &
                  & e1, e2, e3, x0, nx, ny, nz, radius, interpolation, &
                  & lebedev
  !
  ! ... check filplot
  IF (LEN_TRIM(filplot) < 1) THEN
    CALL errore('solvdens', 'filplot is empty', 1)
  END IF
  !
  ! ...
  ! ... Read Input
  ! .................................................................................
  !
  ! ... set default values
  iflag         = -1
  output_format = -1
  fileout       = ''
  e1(:)         = 0.0_DP
  e2(:)         = 0.0_DP
  e3(:)         = 0.0_DP
  x0(:)         = 0.0_DP
  nx            = 0
  ny            = 0
  nz            = 0
  radius        = 1.0_DP
  interpolation = NAME_FOURIER
  lebedev       = 302
  !
  ! ... read the namelist 'plot'
  IF (ionode) THEN
    READ(stdin, plot, iostat=ios)
    !
    ! ... set default of fileout
    IF (LEN_TRIM(fileout) < 1) THEN
      fileext = ''
      IF (0 <= output_format .AND. output_format <= 7) THEN
        fileext = TRIM(formatext(output_format))
      END IF
      !
      IF (LEN_TRIM(fileext) < 1) THEN
        fileout = TRIM(ADJUSTL(prefix)) // '.3drism'
      ELSE
        fileout = TRIM(ADJUSTL(prefix)) // '.3drism_*' // TRIM(fileext)
      END IF
    END IF
  END IF
  !
  CALL mp_bcast(ios, ionode_id, intra_image_comm)
  !
  IF (ios /= 0) THEN
    CALL infomsg('solvdens', 'namelist plot not found or invalid, exiting')
    RETURN
  END IF
  !
  CALL mp_bcast(iflag,         ionode_id, intra_image_comm)
  CALL mp_bcast(output_format, ionode_id, intra_image_comm)
  CALL mp_bcast(fileout,       ionode_id, intra_image_comm)
  CALL mp_bcast(e1,            ionode_id, intra_image_comm)
  CALL mp_bcast(e2,            ionode_id, intra_image_comm)
  CALL mp_bcast(e3,            ionode_id, intra_image_comm)
  CALL mp_bcast(x0,            ionode_id, intra_image_comm)
  CALL mp_bcast(nx,            ionode_id, intra_image_comm)
  CALL mp_bcast(ny,            ionode_id, intra_image_comm)
  CALL mp_bcast(nz,            ionode_id, intra_image_comm)
  CALL mp_bcast(radius,        ionode_id, intra_image_comm)
  CALL mp_bcast(interpolation, ionode_id, intra_image_comm)
  CALL mp_bcast(lebedev,       ionode_id, intra_image_comm)
  !
  ! ...
  ! ... Check Input
  ! .................................................................................
  !
  ! ... check output_format
  IF (output_format == -1) THEN
    CALL infomsg('solvdens', 'output format not set, exiting')
    RETURN
  END IF
  !
  IF (output_format == 1 .OR. output_format == 4 .OR. &
      output_format <  0 .OR. output_format >  7) THEN
    CALL errore ('solvdens', 'output_format is wrong or obsolete', 1)
    RETURN
  END IF
  !
  ! ... check iflag
  IF (iflag == -1) THEN
    CALL infomsg('solvdens', 'iflag not set, exiting')
    RETURN
  END IF
  !
  ! ... check fileout
  IF (LEN_TRIM(fileout) < 1) THEN
    CALL infomsg('solvdens', 'fileout not set, exiting')
    RETURN
  END IF
  !
  IF (iflag == 0) THEN
    ! ... 1D plot (spherical average)
    IF (nx <= 0) THEN
      CALL errore('solvdens', 'wrong nx', 1)
    END IF
    !
  ELSE IF (iflag == 1) THEN
    ! ... 1D plot
    IF ((e1(1) ** 2 + e1(2) ** 2 + e1(3) ** 2) < eps6) THEN
      CALL errore('solvdens', 'missing e1 vector', 1)
    END IF
    !
    IF (nx <= 0) THEN
      CALL errore('solvdens', 'wrong nx', 1)
    END IF
    !
  ELSE IF (iflag == 2) THEN
    ! ... 2D plot
    IF ((e1(1) ** 2 + e1(2) ** 2 + e1(3) ** 2) <  eps6 .OR. &
      & (e2(1) ** 2 + e2(2) ** 2 + e2(3) ** 2) <  eps6) THEN
      CALL errore('solvdens', 'missing e1/e2 vectors', 2)
    END IF
    !
    IF (ABS(e1(1) * e2(1) + e1(2) * e2(2) + e1(3) * e2(3)) > eps6) THEN
      CALL errore('solvdens', 'e1 and e2 are not orthogonal', 2)
    END IF
    !
    IF (nx <= 0 .OR. ny <= 0) THEN
      CALL errore('solvdens', 'wrong nx/ny', 2)
    END IF
    !
    IF (output_format /= 0 .AND. &  ! gnuplot
      & output_format /= 2 .AND. &  ! plotrho.x
      & output_format /= 3 .AND. &  ! XCRYSDEN
      & output_format /= 7) THEN    ! gnuplot x,y,f
      CALL errore('solvdens', 'incompatible iflag/output_format', 2)
    END IF
    !
  ELSE IF (iflag == 3) THEN
    ! ... 3D plot
    IF (ABS(e1(1) * e2(1) + e1(2) * e2(2) + e1(3) * e2(3)) > eps6 .OR. &
      & ABS(e1(1) * e3(1) + e1(2) * e3(2) + e1(3) * e3(3)) > eps6 .OR. &
      & ABS(e2(1) * e3(1) + e2(2) * e3(2) + e2(3) * e3(3)) > eps6 ) THEN
      CALL errore('solvdens', 'e1, e2, e3 are not orthogonal', 3)
    END IF
    !
    IF (output_format /= 3 .AND. &  ! XCRYSDEN (user-supplied 3D region)
      & output_format /= 5 .AND. &  ! XCRYSDEN (using entire FFT grid)
      & output_format /= 6) THEN    ! GAUSSIAN CUBE
      CALL errore('solvdens', 'incompatible iflag/output_format', 3)
    END IF
    !
  ELSE IF (iflag  == 4) THEN
    ! ... 2D polar plot on a sphere
    IF (nx <= 0 .OR. ny <= 0) THEN
      CALL errore ('solvdens', 'wrong nx/ny', 4)
    END IF
    !
    IF (output_format /= 0) THEN    ! gnuplot
      CALL errore('solvdens', 'incompatible iflag/output_format', 4)
    END IF
    !
  ELSE
    ! ... incorrect iflag
    CALL errore('solvdens', 'iflag not implemented', MAX(1, ABS(iflag)))
    !
  END IF
  !
  ! ... check interpolation
  IF (TRIM(interpolation) /= NAME_FOURIER .AND. TRIM(interpolation) /= NAME_BSPLINE) THEN
    CALL errore('solvdens', 'wrong interpolation: ' // TRIM(interpolation), 1)
  END IF
  !
  ! ...
  ! ... Read Punched Data
  ! .................................................................................
  !
  ! ... clean data, if punched
  IF (lpunch) THEN
    CALL clean_pw(.FALSE.)
    !
    dfftp%nr1 = 0
    dfftp%nr2 = 0
    dfftp%nr3 = 0
    !
    dffts%nr1 = 0
    dffts%nr2 = 0
    dffts%nr3 = 0
  END IF
  !
  ! ... read header
  IF (ionode) THEN
    CALL read_rism_header(TRIM(ADJUSTL(filplot)), title, nsite, &
    & dfftp%nr1x, dfftp%nr2x, dfftp%nr3x, dfftp%nr1, dfftp%nr2, dfftp%nr3, &
    & nat, ntyp, ibrav, celldm, at, gcutm, laue)
  END IF
  !
  CALL mp_bcast(title,      ionode_id, intra_image_comm)
  CALL mp_bcast(nsite,      ionode_id, intra_image_comm)
  CALL mp_bcast(dfftp%nr1x, ionode_id, intra_image_comm)
  CALL mp_bcast(dfftp%nr2x, ionode_id, intra_image_comm)
  CALL mp_bcast(dfftp%nr3x, ionode_id, intra_image_comm)
  CALL mp_bcast(dfftp%nr1,  ionode_id, intra_image_comm)
  CALL mp_bcast(dfftp%nr2,  ionode_id, intra_image_comm)
  CALL mp_bcast(dfftp%nr3,  ionode_id, intra_image_comm)
  CALL mp_bcast(nat,        ionode_id, intra_image_comm)
  CALL mp_bcast(ntyp,       ionode_id, intra_image_comm)
  CALL mp_bcast(ibrav,      ionode_id, intra_image_comm)
  CALL mp_bcast(celldm,     ionode_id, intra_image_comm)
  CALL mp_bcast(at,         ionode_id, intra_image_comm)
  CALL mp_bcast(gcutm,      ionode_id, intra_image_comm)
  CALL mp_bcast(laue,       ionode_id, intra_image_comm)
  !
  ! ... allocate atomic data, if not punched
  IF (.NOT. lpunch) THEN
    ALLOCATE(tau(3, nat))
    ALLOCATE(ityp(nat))
  END IF
  !
  ! ... initialize cell variables
  CALL latgen(ibrav, celldm, at(1, 1), at(1, 2), at(1, 3), omega)
  alat   = celldm(1)
  at     = at  / alat
  tpiba  = tpi / alat
  tpiba2 = tpiba ** 2
  !
  CALL recips(at(1, 1), at(1, 2), at(1, 3), bg(1, 1), bg(1, 2), bg(1, 3))
  CALL volume(alat, at(1, 1), at(1, 2), at(1, 3), omega)
  !
  ! ... initialize FFT
  dual       = 4.0_DP
  ecutwfc    = gcutm * tpiba2 / dual
  doublegrid = .FALSE.
  gcutms     = gcutm
  nfft       = dfftp%nr1x * dfftp%nr2x * dfftp%nr3x
  !
  CALL fft_type_allocate(dfftp, at, bg, gcutm,  intra_bgrp_comm, nyfft=nyfft)
  CALL fft_type_allocate(dffts, at, bg, gcutms, intra_bgrp_comm, nyfft=nyfft)
  !
  ! ... allocate punched data
  IF (ionode) THEN
    ALLOCATE(asite(nsite))
    ALLOCATE(rhor(nfft, nsite))
  ELSE
    ALLOCATE(asite(1))
    ALLOCATE(rhor(1, 1))
  END IF
  asite = ''
  rhor  = 0.0_DP
  !
  ! ... read punched data
  IF (ionode) THEN
    CALL plot_rism(TRIM(ADJUSTL(filplot)), title_, nsite_, &
    & nr1x_, nr2x_, nr3x_, nr1_, nr2_, nr3_, nat_, ntyp_, &
    & ibrav_, celldm_, at_, gcutm_, laue_, asite, atm, ityp, zv, tau, rhor, -1)
    !
    ! ... remove noise
    DO isite = 1, nsite
      DO ifft = 1, nfft
        IF (ABS(rhor(ifft, isite)) < 0.9999E-99_DP) THEN
          rhor(ifft, isite) = 0.0_DP
        END IF
      END DO
    END DO
  END IF
  !
  CALL mp_bcast(atm,  ionode_id, intra_image_comm)
  CALL mp_bcast(ityp, ionode_id, intra_image_comm)
  CALL mp_bcast(zv,   ionode_id, intra_image_comm)
  CALL mp_bcast(tau,  ionode_id, intra_image_comm)
  !
  ! ...
  ! ... Modify and Set Variables
  ! .................................................................................
  !
  ! ... re-define vectors
  m1 = SQRT(e1(1) ** 2 + e1(2) ** 2 + e1(3) ** 2)
  IF (ABS(m1) < eps6) THEN
    e1(:) = at(:, 1)
    m1 = SQRT(e1(1) ** 2 + e1(2) ** 2 + e1(3) ** 2)
  END IF
  e1(:) = e1(:) / m1
  !
  m2 = SQRT(e2(1) ** 2 + e2(2) ** 2 + e2(3) ** 2)
  IF (ABS(m2) < eps6) THEN
    e2(:) = at(:, 2)
    m2 = SQRT(e2(1) ** 2 + e2(2) ** 2 + e2(3) ** 2)
  END IF
  e2(:) = e2(:) / m2
  !
  m3 = SQRT(e3(1) ** 2 + e3(2) ** 2 + e3(3) ** 2)
  IF (ABS(m3) < eps6) THEN
    e3(:) = at(:, 3)
    m3 = SQRT(e3(1) ** 2 + e3(2) ** 2 + e3(3) ** 2)
  ENDIF
  e3(:) = e3(:) / m3
  !
  ! ... change unit of radius
  radius = radius / alat
  !
  ! ... check fast3d
  fast3d = (TRIM(interpolation) == NAME_FOURIER)
  !
  fast3d = fast3d .AND. &
  & (e1(2) == 0.0_DP .AND. e1(3) == 0.0_DP) .AND. &
  & (e2(1) == 0.0_DP .AND. e2(3) == 0.0_DP) .AND. &
  & (e3(1) == 0.0_DP .AND. e3(2) == 0.0_DP)
  !
  fast3d = fast3d .AND. &
  & (at(2, 1) == 0.0_DP .AND. at(3, 1) == 0.0_DP) .AND. &
  & (at(1, 2) == 0.0_DP .AND. at(3, 2) == 0.0_DP) .AND. &
  & (at(1, 3) == 0.0_DP .AND. at(2, 3) == 0.0_DP)
  !
  ! ... to avoid a priori fft ?
  avoid_fft = (TRIM(interpolation) /= NAME_FOURIER)
  !
  avoid_fft = avoid_fft .OR. &
            & (iflag == 3 .AND. &         ! 3D plot
            & (output_format == 5 .OR. &  ! XCRYSDEN (using entire FFT grid)
            &  output_format == 6 .OR. &  ! GAUSSIAN CUBE
            &  fast3d))                   ! fast 3D-fitting
  !
  ! ... Laue-RISM cannot support FFT
  IF (laue .AND. (.NOT. avoid_fft)) THEN
    avoid_fft     = .TRUE.
    interpolation = NAME_BSPLINE
    CALL infomsg('solvdens', 'FOURIER cannot be use for Laue-RISM, B-SPLINE is substituted.')
  END IF
  !
  ! ... initialize FFT
  IF (.NOT. avoid_fft) THEN
    nspin      = 1
    gamma_only = .FALSE.
    !
    CALL data_structure(gamma_only)
    CALL allocate_fft()
    !
    CALL ggen (dfftp, gamma_only, at, bg, gcutm, ngm_g, ngm, &
               g, gg, mill, ig_l2g, gstart)
    CALL ggens(dffts, gamma_only, at, g, gg, mill, gcutms, ngms)
    !
    ALLOCATE(rhog(ngm))
  END IF
  !
  ! ...
  ! ... Plot to Files
  ! .................................................................................
  DO isite = 1, nsite
    !
    ! ... set output file name
    IF (ionode) THEN
      fileout0 = ADJUSTL(fileout)
      CALL set_filename(fileout0, isite - 4, asite(isite))
    END IF
    !
    CALL mp_bcast(fileout0, ionode_id, intra_image_comm)
    !
    IF (LEN_TRIM(fileout0) < 1) THEN
      CALL errore('solvdens', 'wrong output file name: ' // TRIM(ADJUSTL(fileout)), isite)
    END IF
    !
    ! ... open output file
    IF (ionode) THEN
      ounit = find_free_unit()
      OPEN(unit=ounit, file=TRIM(ADJUSTL(fileout0)), status='unknown', form='formatted', iostat=ios)
      WRITE(stdout, '(/5X,"Writing data to be plotted to file ",A)') TRIM(ADJUSTL(fileout0))
    END IF
    !
    CALL mp_bcast(ios, ionode_id, intra_image_comm)
    !
    IF (ios /= 0) THEN
      CALL errore('solvdens', 'error to open file: ' // TRIM(ADJUSTL(fileout0)), ABS(ios))
    END IF
    !
    ! ... FFT rho(r) -> rho(G)
    IF (.NOT. avoid_fft) THEN
      !
      ALLOCATE(caux(dfftp%nnr))
      !
#if defined(__MPI)
      ALLOCATE(rhos(nfft))
      ALLOCATE(raux(dfftp%nnr))
      !
      IF (ionode) THEN
        rhos = rhor(:, isite)
      END IF
      !
      raux = 0.0_DP
      CALL scatter_grid(dfftp, rhos, raux)
      caux = CMPLX(raux, 0.0_DP, kind=DP)
      !
      DEALLOCATE(rhos)
      DEALLOCATE(raux)
      !
#else
      caux = CMPLX(rhor(:, isite), 0.0_DP, kind=DP)
      !
#endif
      CALL fwfft('Rho', caux, dfftp)
      !
      rhog(:) = caux(dfftp%nl(:))
      !
      DEALLOCATE(caux)
      !
    END IF
    !
    ! ... plot (rhog in G-space with parallel, rhor in R-space with serial)
    !
    IF (iflag == 0 .OR. iflag == 1) THEN
      ! ... 1D plot
      IF (TRIM(interpolation) == NAME_FOURIER) THEN
        CALL plot_1d(nx, m1, x0, e1, ngm, g, rhog, alat, iflag, ounit)
        !
      ELSE IF (TRIM(interpolation) == NAME_BSPLINE) THEN
        IF (iflag == 0) THEN
          ! ... spherical average by Lebedev Quadrature
          IF (ionode) THEN
            CALL plot_sphere_bspline(nx, lebedev, m1, x0, rhor(:, isite), alat, ounit, laue)
          END IF
        ELSE !IF (iflag == 1) THEN
          IF (ionode) THEN
            CALL plot_1d_bspline(nx, m1, x0, e1, rhor(:, isite), alat, iflag, ounit, laue)
          END IF
        END IF
      END IF
      !
    ELSE IF (iflag == 2) THEN
      ! ... 2D plot
      IF (TRIM(interpolation) == NAME_FOURIER) THEN
        CALL plot_2d(nx, ny, m1, m2, x0, e1, e2, ngm, g, rhog, alat, &
                   & at, nat, tau, atm, ityp, output_format, ounit)
        !
      ELSE IF (TRIM(interpolation) == NAME_BSPLINE) THEN
        IF (ionode) THEN
          CALL plot_2d_bspline(nx, ny, m1, m2, x0, e1, e2, rhor(:, isite), alat, &
                             & at, nat, tau, atm, ityp, output_format, ounit, laue)
        END IF
      END IF
      !
      IF (output_format == 2) THEN  ! plotrho.x
        IF (ionode) THEN
          WRITE(ounit, '(I4)') nat
          WRITE(ounit, '(3F8.4,I3)') ((tau(ipol, na), ipol = 1, 3), 1, na = 1, nat)
          WRITE(ounit, '(F10.6)') celldm(1)
          WRITE(ounit, '(3(3F12.6/))') at
        END IF
      END IF
      !
    ELSE IF (iflag == 3) THEN
      ! ... 3D plot
      !
      ! ... XCRYSDEN FORMAT
      IF (output_format == 5) THEN
        IF (ionode) THEN
          CALL xsf_struct(alat, at, nat, tau, atm, ityp, ounit)
          CALL xsf_fast_datagrid_3d(rhor(:, isite), &
          & dfftp%nr1, dfftp%nr2, dfftp%nr3, dfftp%nr1x, dfftp%nr2x, dfftp%nr3x, at, alat, ounit)
          !
          CALL infomsg('solvdens', 'neither FOURIER nor B-SPLINE is done.')
        END IF
        !
      ! ... GAUSSIAN CUBE FORMAT
      ELSE IF (output_format == 6) THEN
        IF (TRIM(interpolation) == NAME_FOURIER) THEN
          IF (ionode) THEN
            CALL write_cubefile(alat, at, bg, nat, tau, atm, ityp, rhor(:, isite), &
            & dfftp%nr1, dfftp%nr2, dfftp%nr3, dfftp%nr1x, dfftp%nr2x, dfftp%nr3x, ounit)
            !
            CALL infomsg('solvdens', 'FOURIER is not done.')
          END IF
          !
        ELSE IF (TRIM(interpolation) == NAME_BSPLINE) THEN
          IF (nx <= 0 .OR. ny <= 0 .OR. nz <= 0) THEN
            CALL errore('solvdens', 'nx,ny,nz are required', 1)
          END IF
          !
          IF (ionode) THEN
            CALL plot_3d_bspline(celldm(1), at, nat, tau, atm, ityp, rhor(:, isite), &
               & nx, ny, nz, m1, m2, m3, x0, e1, e2, e3, output_format, ounit, rhotot, laue)
          END IF
        END IF
        !
      ! ... XCRYSDEN FORMAT (fast-3D)
      ELSE IF (fast3d) THEN
        IF (ionode) THEN
          CALL plot_fast(celldm(1), at, nat, tau, atm, ityp, &
             & dfftp%nr1x, dfftp%nr2x, dfftp%nr3x, dfftp%nr1, dfftp%nr2, dfftp%nr3, &
             & rhor(:, isite), bg, m1, m2, m3, x0, e1, e2, e3, output_format, ounit, rhotot)
          !
          CALL infomsg('solvdens', 'FAST-3D is done.')
        END IF
        !
      ! ... XCRYSDEN FORMAT (standard)
      ELSE
        IF (nx <= 0 .OR. ny <= 0 .OR. nz <= 0) THEN
          CALL errore('solvdens', 'nx,ny,nz are required', 2)
        END IF
        !
        IF (TRIM(interpolation) == NAME_FOURIER) THEN
          CALL plot_3d(celldm(1), at, nat, tau, atm, ityp, ngm, g, rhog, &
             & nx, ny, nz, m1, m2, m3, x0, e1, e2, e3, output_format, ounit, rhotot)
          !
        ELSE IF (TRIM(interpolation) == NAME_BSPLINE) THEN
          IF (ionode) THEN
            CALL plot_3d_bspline(celldm(1), at, nat, tau, atm, ityp, rhor(:, isite), &
               & nx, ny, nz, m1, m2, m3, x0, e1, e2, e3, output_format, ounit, rhotot, laue)
          END IF
        END IF
        !
      END IF
      !
    ELSE IF (iflag == 4) THEN
      ! ... 2D polar plot on a sphere
      CALL plot_2ds(nx, ny, radius, ngm, g, rhog, output_format, ounit)
      !
    ELSE
      !
      CALL errore('solvdens', 'wrong iflag', MAX(1, ABS(iflag)))
      !
    END IF
    !
  END DO ! end of isite
  !
  WRITE(stdout, '()')
  WRITE(stdout, '(5X,"Plot Type     : ",A)') TRIM(plotname(iflag))
  WRITE(stdout, '(5X,"Output format : ",A)') TRIM(formatname(output_format))
  !
  ! ... deallocate punched data
  DEALLOCATE(asite)
  DEALLOCATE(rhor)
  !
  IF (.NOT. avoid_fft) THEN
    DEALLOCATE(rhog)
  END IF
  !
CONTAINS
  !
  !--------------------------------------------------------------------------
  SUBROUTINE set_filename(filename, num, label)
    !--------------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    CHARACTER(LEN=*), INTENT(INOUT) :: filename
    INTEGER,          INTENT(IN)    :: num
    CHARACTER(LEN=*), INTENT(IN)    :: label
    !
    INTEGER            :: i
    INTEGER            :: nsize
    LOGICAL            :: lhead
    CHARACTER(LEN=1)   :: s
    CHARACTER(LEN=256) :: shead
    CHARACTER(LEN=256) :: stail
    CHARACTER(LEN=256) :: sbody
    CHARACTER(LEN=256) :: filename0
    !
    ! ... create filename
    nsize = LEN_TRIM(filename)
    lhead = .TRUE.
    shead = ''
    stail = ''
    sbody = ''
    !
    DO i = 1, nsize
      s = filename(i:i)
      !
      IF (lhead .AND. s == '*') THEN
        lhead = .FALSE.
        !
      ELSE
        IF (lhead) THEN
          shead = TRIM(shead) // s
        ELSE
          stail = TRIM(stail) // s
        END IF
      END IF
    END DO
    !
    IF (LEN_TRIM(shead) < 1 .AND. LEN_TRIM(stail) < 1) THEN
      filename = ''
      RETURN
    END IF
    !
    IF (lhead) THEN
      shead = TRIM(shead) // '_'
    END IF
    !
    IF (num > 0) THEN
      WRITE(sbody, *) num
      sbody = 'rho#' // TRIM(ADJUSTL(sbody)) // '_'
    END IF
    !
    sbody = TRIM(sbody) // TRIM(ADJUSTL(label))
    !
    IF (LEN_TRIM(sbody) < 1) THEN
      filename = ''
      RETURN
    END IF
    !
    filename = TRIM(shead) // TRIM(sbody) // TRIM(stail)
    !
    ! ... modify filename
    nsize = LEN_TRIM(filename)
    filename0 = ''
    !
    DO i = 1, nsize
      s = filename(i:i)
      IF (s == ' ') THEN
        CYCLE
      END IF
      ! IF (s == '/' .OR. s == '\' .OR. s == '|') THEN
      ! PGI/NVHPC doesn't like '\' 
      IF (s == '/' .OR. s == CHAR(92) .OR. s == '|') THEN
        s = '.'
      END IF
      !
      filename0 = TRIM(filename0) // s
    END DO
    !
    filename = TRIM(filename0)
    !
  END SUBROUTINE set_filename
  !
END SUBROUTINE solvdens
