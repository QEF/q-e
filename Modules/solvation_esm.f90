!
! Copyright (C) 2016 National Institute of Advanced Industrial Science and Technology (AIST)
! [ This code is written by Satomichi Nishihara. ]
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!---------------------------------------------------------------------------
SUBROUTINE solvation_esm_potential(rismt, iref, vref, ierr)
  !---------------------------------------------------------------------------
  !
  ! ... calculate solvation potential of ESM(BC1), from Laue-RISM.
  ! ... calculation is performed around the expanded cell.
  !
  ! ... Variables:
  ! ...   iref: reference position of solvation potential
  ! ...   vref: reference value of solvation potential
  !
  USE cell_base,     ONLY : alat, tpiba, tpiba2
  USE constants,     ONLY : tpi, fpi, e2
  USE err_rism,      ONLY : IERR_RISM_NULL, IERR_RISM_INCORRECT_DATA_TYPE
  USE kinds,         ONLY : DP
  USE lauefft,       ONLY : fw_lauefft_1z_exp, inv_lauefft_1z_exp
  USE mp,            ONLY : mp_sum
  USE rism,          ONLY : rism_type, ITYPE_LAUERISM
  USE rism3d_facade, ONLY : IREFERENCE_AVERAGE, IREFERENCE_RIGHT, IREFERENCE_LEFT
  !
  IMPLICIT NONE
  !
  TYPE(rism_type), INTENT(INOUT) :: rismt
  INTEGER,         INTENT(IN)    :: iref
  REAL(DP),        INTENT(OUT)   :: vref
  INTEGER,         INTENT(OUT)   :: ierr
  !
  INTEGER                  :: iz
  INTEGER                  :: igz
  INTEGER                  :: igxy
  INTEGER                  :: jgxy
  REAL(DP)                 :: z
  REAL(DP)                 :: zright
  REAL(DP)                 :: zleft
  REAL(DP)                 :: zstart
  REAL(DP)                 :: dz
  REAL(DP)                 :: gz
  REAL(DP)                 :: gxy
  REAL(DP)                 :: ggxy
  REAL(DP)                 :: fac1  ! factors to    -> Energy/G/G
  REAL(DP)                 :: fac2  ! convert units -> Energy*R/G
  REAL(DP)                 :: fac3  !               -> Energy*R*R
  REAL(DP)                 :: phir
  REAL(DP)                 :: phil
  REAL(DP)                 :: rho0
  REAL(DP)                 :: realr
  REAL(DP)                 :: reall
  REAL(DP)                 :: imager
  REAL(DP)                 :: imagel
  REAL(DP)                 :: vsolv
  REAL(DP)                 :: vsolu
  COMPLEX(DP)              :: coeffr
  COMPLEX(DP)              :: coeffl
  COMPLEX(DP), ALLOCATABLE :: rhogt(:,:)
  COMPLEX(DP), ALLOCATABLE :: rhogz(:)
  COMPLEX(DP), ALLOCATABLE :: vpott(:,:)
  COMPLEX(DP), ALLOCATABLE :: expigzr(:)
  COMPLEX(DP), ALLOCATABLE :: expigzl(:)
  !
  COMPLEX(DP), PARAMETER   :: C_ZERO = CMPLX(0.0_DP, 0.0_DP, kind=DP)
  !
  ! ... check data type
  IF (rismt%itype /= ITYPE_LAUERISM) THEN
    ierr = IERR_RISM_INCORRECT_DATA_TYPE
    RETURN
  END IF
  !
  IF (rismt%nrzl < rismt%lfft%nrz) THEN
    ierr = IERR_RISM_INCORRECT_DATA_TYPE
    RETURN
  END IF
  !
  IF (rismt%ngxy < rismt%lfft%ngxy) THEN
    ierr = IERR_RISM_INCORRECT_DATA_TYPE
    RETURN
  END IF
  !
  ! ... allocate memory
  IF (rismt%lfft%ngz_x * rismt%lfft%ngxy > 0) THEN
    ALLOCATE(rhogt(  rismt%lfft%ngz_x, rismt%lfft%ngxy))
    ALLOCATE(vpott(  rismt%lfft%ngz_x, rismt%lfft%ngxy))
  END IF
  IF (rismt%lfft%ngz_x > 0) THEN
    ALLOCATE(rhogz(  rismt%lfft%ngz_x))
    ALLOCATE(expigzr(rismt%lfft%ngz_x))
    ALLOCATE(expigzl(rismt%lfft%ngz_x))
  END IF
  !
  ! ... set variables
  zright = rismt%lfft%zright
  zleft  = rismt%lfft%zleft
  zstart = rismt%lfft%zleft + rismt%lfft%zoffset
  dz     = rismt%lfft%zstep
  fac1   = e2 * fpi / tpiba2        ! -> Energy/G/G
  fac2   = e2 * fpi * alat / tpiba  ! -> Energy/R/G
  fac3   = e2 * fpi * alat * alat   ! -> Energy*R*R
  !
  ! ... initialize reference potential
  vref = 0.0_DP
  !
  ! ... calculate exp(i*gz*zright) and exp(i*gz*zleft)
!$omp parallel do default(shared) private(igz, gz, phir, phil)
  DO igz = 1, rismt%lfft%ngz_x
    gz   = rismt%lfft%gz_x(igz)
    phir = tpi * gz * zright
    phil = tpi * gz * zleft
    expigzr(igz) = CMPLX(COS(phir), SIN(phir), kind=DP)
    expigzl(igz) = CMPLX(COS(phil), SIN(phil), kind=DP)
  END DO
!$omp end parallel do
  !
  ! ... 1D-FFT of charge: Laue-rep. -> G-space
  IF (rismt%lfft%ngz_x * rismt%lfft%ngxy > 0) THEN
    rhogt = C_ZERO
    CALL fw_lauefft_1z_exp(rismt%lfft, rismt%rhog, rismt%nrzl, rhogt, rismt%lfft%ngz_x)
  END IF
  !
  ! ... Hartree potential: part of 4pi/G^2
  IF (rismt%lfft%ngz_x * rismt%lfft%ngxy > 0) THEN
    vpott = C_ZERO
  END IF
  !
  DO igxy = rismt%lfft%gxystart, rismt%lfft%ngxy
    ggxy = rismt%lfft%ggxy(igxy)
!$omp parallel do default(shared) private(igz, gz)
    DO igz = 1, rismt%lfft%ngz_x
      gz = rismt%lfft%gz_x(igz)
      vpott(igz, igxy) = rhogt(igz, igxy) * (fac1 / (gz * gz + ggxy))
    END DO
!$omp end parallel do
  END DO
  !
  IF (rismt%lfft%gxystart > 1) THEN
    igxy = 1
!$omp parallel do default(shared) private(igz, gz)
    DO igz = 1, rismt%lfft%ngz_x
      IF (igz /= rismt%lfft%gzzero_x) THEN
        gz = rismt%lfft%gz_x(igz)
        vpott(igz, igxy) = rhogt(igz, igxy) * (fac1 / (gz * gz))
      END IF
    END DO
!$omp end parallel do
  END IF
  !
  ! ... 1D-FFT of Hartree potential: G-space -> Laue-rep.
  IF (rismt%nrzl * rismt%ngxy > 0) THEN
    rismt%vpot = C_ZERO
  END IF
  IF (rismt%lfft%ngz_x * rismt%lfft%ngxy > 0) THEN
    CALL inv_lauefft_1z_exp(rismt%lfft, vpott, rismt%lfft%ngz_x, rismt%vpot, rismt%nrzl)
  END IF
  !
  ! ...
  ! ... Hartree potential, when Gxy /= 0
  ! ...
  DO igxy = rismt%lfft%gxystart, rismt%lfft%ngxy
    !
    jgxy = rismt%nrzl * (igxy - 1)
    gxy  = rismt%lfft%gnxy(igxy)
    !
    IF (rismt%lfft%ngz_x > 0) THEN
      rhogz(:) = rhogt(:, igxy)
    END IF
    !
    ! ... coeffr means
    !      ----           exp( i*gz*zright)
    !      >    rho(g) * -------------------
    !      ----              i*gz - gxy
    !       gz
    !
    ! ... coeffl means
    !      ----           exp( i*gz*zleft)
    !      >    rho(g) * -------------------
    !      ----              i*gz + gxy
    !       gz
    !
    coeffr = C_ZERO
    coeffl = C_ZERO
    !
!$omp parallel do default(shared) private(igz, gz) reduction(+:coeffr, coeffl)
    DO igz = 1, rismt%lfft%ngz_x
      gz = rismt%lfft%gz_x(igz)
      coeffr = coeffr + rhogz(igz) * expigzr(igz) / CMPLX(-gxy, gz, kind=DP)
      coeffl = coeffl + rhogz(igz) * expigzl(igz) / CMPLX( gxy, gz, kind=DP)
    END DO
!$omp end parallel do
    !
    ! ... when zleft <= z <= zright, potential is
    !
    !            exp( gxy*(z-zright))    ----           exp( i*gz*zright)
    !  (+4pi) * ---------------------- * >    rho(g) * -------------------
    !                 2 * gxy            ----              i*gz - gxy
    !                                     gz
    !
    !            exp(-gxy*(z+zleft))     ----           exp( i*gz*zleft)
    !  (-4pi) * ---------------------- * >    rho(g) * -------------------
    !                 2 * gxy            ----              i*gz + gxy
    !                                     gz
    !
!$omp parallel do default(shared) private(iz, z)
    DO iz = 1, rismt%lfft%nrz
      z = zstart + dz * DBLE(iz - 1)
      rismt%vpot(iz + jgxy) = rismt%vpot(iz + jgxy) + fac1 * ( &
      & + (0.5_DP / gxy) * EXP( tpi * gxy * (z - zright)) * coeffr &
      & - (0.5_DP / gxy) * EXP(-tpi * gxy * (z - zleft )) * coeffl )
    END DO
!$omp end parallel do
    !
  END DO
  !
  ! ...
  ! ... Hartree potential, when Gxy = 0
  ! ...
  IF (rismt%lfft%gxystart > 1) THEN
    !
    igxy = 1
    jgxy = rismt%nrzl * (igxy - 1)
    gxy  = rismt%lfft%gnxy(igxy)
    !
    IF (rismt%lfft%ngz_x > 0) THEN
      rhogz(:) = rhogt(:, igxy)
      rho0 = rhogz(rismt%lfft%gzzero_x)
    END IF
    !
    ! ... realr means
    !      ----  Re[ rho(gz) * exp( i*gz*zright) ]
    !      >    -----------------------------------
    !      ----                gz^2
    !      gz>0
    !
    ! ... reall means
    !      ----  Re[ rho(gz) * exp( i*gz*zleft) ]
    !      >    -----------------------------------
    !      ----                gz^2
    !      gz>0
    !
    ! ... imager means
    !      ----  Im[ rho(gz) * exp( i*gz*zright) ]
    !      >    -----------------------------------
    !      ----                gz
    !      gz>0
    !
    ! ... imagel means
    !      ----  Im[ rho(gz) * exp( i*gz*zleft) ]
    !      >    -----------------------------------
    !      ----                gz
    !      gz>0
    !
    realr  = 0.0_DP
    reall  = 0.0_DP
    imager = 0.0_DP
    imagel = 0.0_DP
    !
!$omp parallel do default(shared) private(igz, gz) reduction(+:realr, reall, imager, imagel)
    DO igz = (rismt%lfft%gzzero_x + 1), rismt%lfft%ngz_x
      gz = rismt%lfft%gz_x(igz)
      realr  = realr  + DBLE( rhogz(igz) * expigzr(igz)) / gz / gz
      reall  = reall  + DBLE( rhogz(igz) * expigzl(igz)) / gz / gz
      imager = imager + AIMAG(rhogz(igz) * expigzr(igz)) / gz
      imagel = imagel + AIMAG(rhogz(igz) * expigzl(igz)) / gz
    END DO
!$omp end parallel do
    !
    ! ... when -z0 <= z <= z0, potential is
    !
    !           ----  Re[ rho(gz) * exp( i*gz*zright) ]
    !  (-4pi) * >    -----------------------------------
    !           ----                gz^2
    !           gz>0
    !
    !           ----  Re[ rho(gz) * exp( i*gz*zrleft) ]
    !  (-4pi) * >    -----------------------------------
    !           ----                gz^2
    !           gz>0
    !
    !                        ----  Im[ rho(gz) * exp( i*gz*zright) ]
    !  (+4pi) * (z-zright) * >    -----------------------------------
    !                        ----                gz
    !                        gz>0
    !
    !                        ----  Im[ rho(gz) * exp( i*gz*zleft) ]
    !  (+4pi) * (z-zleft)  * >    -----------------------------------
    !                        ----                gz
    !                        gz>0
    !
    !  (-4pi) * ((z-zright)^2 + (z-zleft)^2) * rho(0) / 4
    !
!$omp parallel do default(shared) private(iz, z)
    DO iz = 1, rismt%lfft%nrz
      z = zstart + dz * DBLE(iz - 1)
      rismt%vpot(iz + jgxy) = rismt%vpot(iz + jgxy) + CMPLX( &
      & + fac1 * ( - realr - reall )                &
      & + fac2 * ( + (z - zright) * imager          &
      &            + (z - zleft ) * imagel )        &
      & + fac3 * 0.25_DP * rho0 * (                 &
      &            - (z - zright) * (z - zright)    &
      &            - (z - zleft ) * (z - zleft ) )  &
      & , 0.0_DP, kind=DP)
    END DO
!$omp end parallel do
    !
    ! ... modify reference of solvation potential
    vsolv = 0.0_DP
    vsolu = 0.0_DP
    !
    IF (iref == IREFERENCE_AVERAGE) THEN
      vsolv = 0.0_DP
      vsolu = 0.0_DP
      !
    ELSE IF (iref == IREFERENCE_RIGHT) THEN
      vsolv = fac1 * ( + realr - reall) &
          & + fac2 * ( + zright * imager - zleft * imagel ) &
          & + fac3 * 0.25_DP * rho0 * ( + zright * zright - zleft * zleft )
      vsolu = AIMAG(rismt%vright(igxy))
      !
    ELSE IF (iref == IREFERENCE_LEFT) THEN
      vsolv = fac1 * ( - realr + reall) &
          & + fac2 * ( - zright * imager + zleft * imagel ) &
          & + fac3 * 0.25_DP * rho0 * ( - zright * zright + zleft * zleft )
      vsolu = AIMAG(rismt%vleft(igxy))
    END IF
    !
    vref = vsolv + vsolu
    !
!$omp parallel do default(shared) private(iz)
    DO iz = 1, rismt%lfft%nrz
      rismt%vpot(iz + jgxy) = rismt%vpot(iz + jgxy) - CMPLX(vref, 0.0_DP, kind=DP)
    END DO
!$omp end parallel do
    !
  END IF
  !
  CALL mp_sum(vref, rismt%mp_site%intra_sitg_comm)
  !
  ! ... deallocate memory
  IF (rismt%lfft%ngz_x * rismt%lfft%ngxy > 0) THEN
    DEALLOCATE(rhogt)
    DEALLOCATE(vpott)
  END IF
  IF (rismt%lfft%ngz_x > 0) THEN
    DEALLOCATE(rhogz)
    DEALLOCATE(expigzr)
    DEALLOCATE(expigzl)
  END IF
  !
  ! ... normally done
  ierr = IERR_RISM_NULL
  !
END SUBROUTINE solvation_esm_potential
!
!---------------------------------------------------------------------------
SUBROUTINE solvation_esm_force(rismt, alpha, force, ierr)
  !---------------------------------------------------------------------------
  !
  ! ... calculate solvation force of ESM(BC1), from Laue-RISM.
  ! ... local potential is derived from Gaussian functions:
  ! ...
  ! ...                      1               [   |r - R|^2 ]
  ! ...   rho(r) = -------------------- * exp[- -----------]
  ! ...             pi^(3/2) * alpha^3       [    alpha^2  ]  .
  !
  ! ... Variables:
  ! ...   alpha: gaussian width (in alat units)
  ! ...   force: solvation force from local potential of ESM(BC1)
  !
  USE cell_base,     ONLY : alat
  USE constants,     ONLY : pi, tpi, e2
  USE control_flags, ONLY : gamma_only
  USE err_rism,      ONLY : IERR_RISM_NULL, IERR_RISM_INCORRECT_DATA_TYPE
  USE gvect,         ONLY : eigts1, eigts2
  USE ions_base,     ONLY : nat, tau, ityp, zv
  USE kinds,         ONLY : DP
  USE mp,            ONLY : mp_sum
  USE rism,          ONLY : rism_type, ITYPE_LAUERISM
  !
  IMPLICIT NONE
  !
  TYPE(rism_type), INTENT(IN)  :: rismt
  REAL(DP),        INTENT(IN)  :: alpha
  REAL(DP),        INTENT(OUT) :: force(1:3, 1:*)
  INTEGER,         INTENT(OUT) :: ierr
  !
  INTEGER                  :: ia
  INTEGER                  :: it
  INTEGER                  :: ipol
  INTEGER                  :: iz
  INTEGER                  :: igxy
  INTEGER                  :: jgxy
  INTEGER                  :: mx, my
  REAL(DP)                 :: z
  REAL(DP)                 :: za
  REAL(DP)                 :: zstart
  REAL(DP)                 :: dz
  REAL(DP)                 :: gx, gy
  REAL(DP)                 :: gxy
  REAL(DP)                 :: qa
  REAL(DP)                 :: mult
  REAL(DP)                 :: rterm1
  REAL(DP)                 :: rterm2
  REAL(DP)                 :: rcoeff
  REAL(DP)                 :: rhogr
  REAL(DP)                 :: rhogi
  REAL(DP)                 :: dvlocr
  REAL(DP)                 :: dvloci
  REAL(DP)                 :: forctmp(3)
  REAL(DP),    ALLOCATABLE :: forcesm(:,:)
  COMPLEX(DP)              :: ccoeff
  COMPLEX(DP)              :: strf_xy
  COMPLEX(DP), ALLOCATABLE :: dvloc(:,:)
  COMPLEX(DP), ALLOCATABLE :: rhogz(:)
  !
  COMPLEX(DP), PARAMETER   :: C_ZERO = CMPLX(0.0_DP, 0.0_DP, kind=DP)
  !
  ! ... check data type
  IF (rismt%itype /= ITYPE_LAUERISM) THEN
    ierr = IERR_RISM_INCORRECT_DATA_TYPE
    RETURN
  END IF
  !
  IF (rismt%nrzl < rismt%lfft%nrz) THEN
    ierr = IERR_RISM_INCORRECT_DATA_TYPE
    RETURN
  END IF
  !
  IF (rismt%ngxy < rismt%lfft%ngxy) THEN
    ierr = IERR_RISM_INCORRECT_DATA_TYPE
    RETURN
  END IF
  !
  ! ... allocate memory
  IF (nat > 0) THEN
    ALLOCATE(forcesm(3, nat))
  END IF
  IF (rismt%lfft%nrz > 0) THEN
    ALLOCATE(dvloc(3, rismt%lfft%nrz))
    ALLOCATE(rhogz(   rismt%lfft%nrz))
  END IF
  !
  ! ... set variables
  zstart = rismt%lfft%zleft + rismt%lfft%zoffset
  dz     = rismt%lfft%zstep
  IF (gamma_only) THEN
    mult = 2.0_DP
  ELSE
    mult = 1.0_DP
  END IF
  !
  ! ... initialize force
  forcesm = 0.0_DP
  !
  ! ...
  ! ... local potential, when Gxy /= 0
  ! ...
  DO igxy = rismt%lfft%gxystart, rismt%lfft%ngxy
    !
    jgxy = rismt%nrzl * (igxy - 1)
    gx   = rismt%lfft%gxy(1, igxy)
    gy   = rismt%lfft%gxy(2, igxy)
    gxy  = rismt%lfft%gnxy(igxy)
    mx   = rismt%lfft%millxy(1, igxy)
    my   = rismt%lfft%millxy(2, igxy)
    !
    IF (rismt%lfft%nrz > 0) THEN
      rhogz(:) = rismt%rhog((1 + jgxy):(rismt%lfft%nrz + jgxy))
    END IF
    !
    DO ia = 1, nat
      !
      it = ityp(ia)
      qa = zv(it)
      za = tau(3, ia)
      !
      ! ... ccoeff means
      !
      !      pi*exp(-i(gx*xa+gy*ya))
      !
      strf_xy = eigts1(mx, ia) * eigts2(my, ia)
      rcoeff  = -qa * e2 * pi
      ccoeff  = rcoeff * strf_xy
      !
      dvloc = C_ZERO
      !
!$omp parallel do default(shared) private(iz, z, rterm1, rterm2)
      DO iz = 1, rismt%lfft%nrz
        z = zstart + dz * DBLE(iz - 1)
        !
        ! ... rterm1 means
        !                              gxy*alpha     z-za
        !     exp( gxy*(z-za)) * erfc(----------- + -------)
        !                                  2         alpha
        !
        ! ... rterm2 means
        !                              gxy*alpha     z-za
        !     exp(-gxy*(z-za)) * erfc(----------- - -------)
        !                                  2         alpha
        !
        ! ... NOTE: to avoid overflows,
        ! ...       exp(var1)*erfc(var2) = exp(var1 + log(erfc(var2))) .
        !
        rterm1 = EXP( tpi * gxy * (z - za) + LOG(erfc(0.5_DP * tpi * gxy * alpha + (z - za) / alpha)))
        rterm2 = EXP(-tpi * gxy * (z - za) + LOG(erfc(0.5_DP * tpi * gxy * alpha - (z - za) / alpha)))
        !
        ! ... derive by X
        !
        !      -i*gx                              [                          gxy*alpha     z-za
        !     ------- * pi*exp(-i(gx*xa+gy*ya)) * [ exp( gxy*(z-za)) * erfc(----------- + -------)
        !       gxy                               [                              2         alpha
        !                                                                    gxy*alpha     z-za    ]
        !                                         + exp(-gxy*(z-za)) * erfc(----------- - -------) ]
        !                                                                        2         alpha   ]
        !
        ! ... derive by Y
        !
        !      -i*gx                              [                          gxy*alpha     z-za
        !     ------- * pi*exp(-i(gx*xa+gy*ya)) * [ exp( gxy*(z-za)) * erfc(----------- + -------)
        !       gxy                               [                              2         alpha
        !                                                                    gxy*alpha     z-za    ]
        !                                         + exp(-gxy*(z-za)) * erfc(----------- - -------) ]
        !                                                                        2         alpha   ]
        ! ... derive by Z
        !
        !                                 [                          gxy*alpha     z-za
        !     - pi*exp(-i(gx*xa+gy*ya)) * [ exp( gxy*(z-za)) * erfc(----------- + -------)
        !                                 [                              2         alpha
        !                                                            gxy*alpha     z-za    ]
        !                                 - exp(-gxy*(z-za)) * erfc(----------- - -------) ]
        !                                                                2         alpha   ]
        !
        dvloc(1, iz) = CMPLX(0.0_DP, -gx / gxy, kind=DP) * ccoeff * (rterm1 + rterm2)
        dvloc(2, iz) = CMPLX(0.0_DP, -gy / gxy, kind=DP) * ccoeff * (rterm1 + rterm2)
        dvloc(3, iz) = -ccoeff * (rterm1 - rterm2)
      END DO
!$omp end parallel do
      !
      forctmp = 0.0_DP
!$omp parallel do default(shared) private(iz, ipol, rhogr, rhogi, dvlocr, dvloci) reduction(+:forctmp)
      DO iz = 1, rismt%lfft%izleft_gedge
        rhogr = -DBLE( rhogz(iz))
        rhogi = -AIMAG(rhogz(iz))
        DO ipol = 1, 3
          dvlocr = DBLE( dvloc(ipol, iz))
          dvloci = AIMAG(dvloc(ipol, iz))
          forctmp(ipol) = forctmp(ipol) - mult * (dvlocr * rhogr + dvloci * rhogi)
        END DO
      END DO
!$omp end parallel do
      forcesm(:, ia) = forcesm(:, ia) + forctmp(:)
      !
      forctmp = 0.0_DP
!$omp parallel do default(shared) private(iz, ipol, rhogr, rhogi, dvlocr, dvloci) reduction(+:forctmp)
      DO iz = rismt%lfft%izright_gedge, rismt%lfft%nrz
        rhogr = -DBLE( rhogz(iz))
        rhogi = -AIMAG(rhogz(iz))
        DO ipol = 1, 3
          dvlocr = DBLE( dvloc(ipol, iz))
          dvloci = AIMAG(dvloc(ipol, iz))
          forctmp(ipol) = forctmp(ipol) - mult * (dvlocr * rhogr + dvloci * rhogi)
        END DO
      END DO
!$omp end parallel do
      forcesm(:, ia) = forcesm(:, ia) + forctmp(:)
      !
    END DO
    !
  END DO
  !
  ! ...
  ! ... local potential, when Gxy = 0
  ! ...
  IF (rismt%lfft%gxystart > 1) THEN
    !
    igxy = 1
    jgxy = rismt%nrzl * (igxy - 1)
    gx   = rismt%lfft%gxy(1, igxy)
    gy   = rismt%lfft%gxy(2, igxy)
    gxy  = rismt%lfft%gnxy(igxy)
    !
    IF (rismt%lfft%nrz > 0) THEN
      rhogz(:) = rismt%rhog((1 + jgxy):(rismt%lfft%nrz + jgxy))
    END IF
    !
    DO ia = 1, nat
      !
      it = ityp(ia)
      qa = zv(it)
      za = tau(3, ia)
      !
      dvloc = C_ZERO
      !
!$omp parallel do default(shared) private(iz, z)
      DO iz = 1, rismt%lfft%nrz
        z = zstart + dz * DBLE(iz - 1)
        !
        ! ... derive by Z
        !
        !               z-za
        !    2pi * erf(-------)
        !               alpha
        !
        dvloc(1, iz) = C_ZERO
        dvloc(2, iz) = C_ZERO
        dvloc(3, iz) = CMPLX((-qa * e2 * tpi) * erf((z - za) / alpha), 0.0_DP, kind=DP)
      END DO
!$omp end parallel do
      !
      forctmp = 0.0_DP
!$omp parallel do default(shared) private(iz, ipol, rhogr, dvlocr) reduction(+:forctmp)
      DO iz = 1, rismt%lfft%izleft_gedge
        rhogr = -DBLE(rhogz(iz))
        DO ipol = 1, 3
          dvlocr = DBLE(dvloc(ipol, iz))
          forctmp(ipol) = forctmp(ipol) - dvlocr * rhogr
        END DO
      END DO
!$omp end parallel do
      forcesm(:, ia) = forcesm(:, ia) + forctmp(:)
      !
      forctmp = 0.0_DP
!$omp parallel do default(shared) private(iz, ipol, rhogr, dvlocr) reduction(+:forctmp)
      DO iz = rismt%lfft%izright_gedge, rismt%lfft%nrz
        rhogr = -DBLE(rhogz(iz))
        DO ipol = 1, 3
          dvlocr = DBLE(dvloc(ipol, iz))
          forctmp(ipol) = forctmp(ipol) - dvlocr * rhogr
        END DO
      END DO
!$omp end parallel do
      forcesm(:, ia) = forcesm(:, ia) + forctmp(:)
      !
    END DO
    !
  END IF
  !
  IF (nat > 0) THEN
    CALL mp_sum(forcesm, rismt%mp_site%intra_sitg_comm)
    force(1:3, 1:nat) = forcesm(1:3, 1:nat) * dz * alat
  END IF
  !
  ! ... deallocate memory
  IF (nat > 0) THEN
    DEALLOCATE(forcesm)
  END IF
  IF (rismt%lfft%nrz > 0) THEN
    DEALLOCATE(dvloc)
    DEALLOCATE(rhogz)
  END IF
  !
  ! ... normally done
  ierr = IERR_RISM_NULL
  !
END SUBROUTINE solvation_esm_force
!
!---------------------------------------------------------------------------
SUBROUTINE solvation_esm_stress(rismt, alpha, sigma, ierr)
  !---------------------------------------------------------------------------
  !
  ! ... calculate solvation stress of ESM(BC1), from Laue-RISM.
  ! ... local potential is derived from Gaussian functions:
  ! ...
  ! ...                      1               [   |r - R|^2 ]
  ! ...   rho(r) = -------------------- * exp[- -----------]
  ! ...             pi^(3/2) * alpha^3       [    alpha^2  ]  .
  !
  ! ... Variables:
  ! ...   alpha: gaussian width (in alat units)
  ! ...   sigma: solvation stress from local potential of ESM(BC1)
  !
  USE cell_base,     ONLY : alat
  USE constants,     ONLY : pi, tpi, e2
  USE control_flags, ONLY : gamma_only
  USE err_rism,      ONLY : IERR_RISM_NULL, IERR_RISM_INCORRECT_DATA_TYPE
  USE gvect,         ONLY : eigts1, eigts2
  USE ions_base,     ONLY : nat, tau, ityp, zv
  USE kinds,         ONLY : DP
  USE mp,            ONLY : mp_sum
  USE rism,          ONLY : rism_type, ITYPE_LAUERISM
  !
  IMPLICIT NONE
  !
  TYPE(rism_type), INTENT(IN)  :: rismt
  REAL(DP),        INTENT(IN)  :: alpha
  REAL(DP),        INTENT(OUT) :: sigma(3, 3)
  INTEGER,         INTENT(OUT) :: ierr
  !
  INTEGER                  :: ia
  INTEGER                  :: it
  INTEGER                  :: ipol
  INTEGER                  :: iz
  INTEGER                  :: igxy
  INTEGER                  :: jgxy
  INTEGER                  :: mx, my
  REAL(DP)                 :: z
  REAL(DP)                 :: za
  REAL(DP)                 :: zstart
  REAL(DP)                 :: dz
  REAL(DP)                 :: gx, gy
  REAL(DP)                 :: gxy
  REAL(DP)                 :: qa
  REAL(DP)                 :: mult
  REAL(DP)                 :: rterm1
  REAL(DP)                 :: rterm2
  REAL(DP)                 :: rcoeff
  REAL(DP)                 :: rhogr
  REAL(DP)                 :: rhogi
  REAL(DP)                 :: dvlocr
  REAL(DP)                 :: dvloci
  REAL(DP)                 :: sigmaesm(3, 3)
  COMPLEX(DP)              :: ccoeff
  COMPLEX(DP)              :: strf_xy
  COMPLEX(DP), ALLOCATABLE :: dvloc(:,:)
  COMPLEX(DP), ALLOCATABLE :: rhogz(:)
  !
  COMPLEX(DP), PARAMETER   :: C_ZERO = CMPLX(0.0_DP, 0.0_DP, kind=DP)
  !
  ! ... check data type
  IF (rismt%itype /= ITYPE_LAUERISM) THEN
    ierr = IERR_RISM_INCORRECT_DATA_TYPE
    RETURN
  END IF
  !
  IF (rismt%nrzl < rismt%lfft%nrz) THEN
    ierr = IERR_RISM_INCORRECT_DATA_TYPE
    RETURN
  END IF
  !
  IF (rismt%ngxy < rismt%lfft%ngxy) THEN
    ierr = IERR_RISM_INCORRECT_DATA_TYPE
    RETURN
  END IF
  !
  ! ... allocate memory
  IF (rismt%lfft%nrz > 0) THEN
    ALLOCATE(dvloc(3, rismt%lfft%nrz))
    ALLOCATE(rhogz(   rismt%lfft%nrz))
  END IF
  !
  ! ... set variables
  zstart = rismt%lfft%zleft + rismt%lfft%zoffset
  dz     = rismt%lfft%zstep
  IF (gamma_only) THEN
    mult = 2.0_DP
  ELSE
    mult = 1.0_DP
  END IF
  !
  ! ... initialize stress
  sigmaesm = 0.0_DP
  !
  ! ...
  ! ... local potential, when Gxy /= 0
  ! ...
  DO igxy = rismt%lfft%gxystart, rismt%lfft%ngxy
    !
    jgxy = rismt%nrzl * (igxy - 1)
    gx   = rismt%lfft%gxy(1, igxy)
    gy   = rismt%lfft%gxy(2, igxy)
    gxy  = rismt%lfft%gnxy(igxy)
    mx   = rismt%lfft%millxy(1, igxy)
    my   = rismt%lfft%millxy(2, igxy)
    !
    IF (rismt%lfft%nrz > 0) THEN
      rhogz(:) = rismt%rhog((1 + jgxy):(rismt%lfft%nrz + jgxy))
    END IF
    !
    DO ia = 1, nat
      !
      it = ityp(ia)
      qa = zv(it)
      za = tau(3, ia)
      !
      ! TODO
      ! TODO set sigmaesm
      ! TODO
      !
    END DO
    !
  END DO
  !
  ! ...
  ! ... local potential, when Gxy = 0
  ! ...
  IF (rismt%lfft%gxystart > 1) THEN
    !
    igxy = 1
    jgxy = rismt%nrzl * (igxy - 1)
    gx   = rismt%lfft%gxy(1, igxy)
    gy   = rismt%lfft%gxy(2, igxy)
    gxy  = rismt%lfft%gnxy(igxy)
    !
    IF (rismt%lfft%nrz > 0) THEN
      rhogz(:) = rismt%rhog((1 + jgxy):(rismt%lfft%nrz + jgxy))
    END IF
    !
    DO ia = 1, nat
      !
      it = ityp(ia)
      qa = zv(it)
      za = tau(3, ia)
      !
      ! TODO
      ! TODO set sigmaesm
      ! TODO
      !
    END DO
    !
  END IF
  !
  CALL mp_sum(sigmaesm, rismt%mp_site%intra_sitg_comm)
  sigma = sigmaesm * dz * alat
  !
  ! ... deallocate memory
  IF (rismt%lfft%nrz > 0) THEN
    DEALLOCATE(dvloc)
    DEALLOCATE(rhogz)
  END IF
  !
  ! ... normally done
  ierr = IERR_RISM_NULL
  !
END SUBROUTINE solvation_esm_stress
