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
SUBROUTINE potential_esm_hartree(rismt, rhog, vlaue, vright, vleft, ierr)
  !---------------------------------------------------------------------------
  !
  ! ... calculate Hartree potential of ESM(BC1), for Laue-RISM.
  ! ... however, " 4pi*Rho(G)/G^2 " is NOT included.
  ! ...
  ! ... Variables:
  ! ...   rhog:   electronic density in G-space
  ! ...   vlaue:  coulomb potential in Laue-rep.
  ! ...           its leading dimension is rismt%nrzl.
  ! ...   vright: coefficient of coulomb potential when zright < z
  ! ...   vleft:  coefficient of coulomb potential when z < zleft
  !
  USE cell_base, ONLY : at, alat, tpiba, tpiba2
  USE constants, ONLY : tpi, fpi, e2
  USE err_rism,  ONLY : IERR_RISM_NULL, IERR_RISM_INCORRECT_DATA_TYPE
  USE kinds,     ONLY : DP
  USE rism,      ONLY : rism_type, ITYPE_LAUERISM
  !
  IMPLICIT NONE
  !
  TYPE(rism_type), INTENT(IN)    :: rismt
  COMPLEX(DP),     INTENT(IN)    :: rhog(1:*)
  COMPLEX(DP),     INTENT(INOUT) :: vlaue(1:*)
  COMPLEX(DP),     INTENT(INOUT) :: vright(1:*)
  COMPLEX(DP),     INTENT(INOUT) :: vleft(1:*)
  INTEGER,         INTENT(OUT)   :: ierr
  !
  INTEGER                  :: iz
  INTEGER                  :: ig
  INTEGER                  :: igz
  INTEGER                  :: igxy
  INTEGER                  :: jgxy
  REAL(DP)                 :: z
  REAL(DP)                 :: z0
  REAL(DP)                 :: zright
  REAL(DP)                 :: zleft
  REAL(DP)                 :: zstart
  REAL(DP)                 :: dz
  REAL(DP)                 :: gz
  REAL(DP)                 :: gxy
  REAL(DP)                 :: fac1  ! factors to    -> Energy/G/G
  REAL(DP)                 :: fac2  ! convert units -> Energy*R/G
  REAL(DP)                 :: fac3  !               -> Energy*R*R
  REAL(DP)                 :: phi
  REAL(DP)                 :: cos1, sin1
  REAL(DP)                 :: cos2, sin2
  REAL(DP)                 :: rho0
  COMPLEX(DP)              :: coeff1
  COMPLEX(DP)              :: coeff2
  COMPLEX(DP)              :: coeff3
  COMPLEX(DP)              :: coeff4
  COMPLEX(DP), ALLOCATABLE :: rhogt(:,:)
  COMPLEX(DP), ALLOCATABLE :: rhogz(:)
  COMPLEX(DP), ALLOCATABLE :: expigz(:)
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
  IF (rismt%lfft%ngz * rismt%lfft%ngxy > 0) THEN
    ALLOCATE(rhogt( rismt%lfft%ngz, rismt%lfft%ngxy))
  END IF
  IF (rismt%lfft%ngz > 0) THEN
    ALLOCATE(rhogz( rismt%lfft%ngz))
    ALLOCATE(expigz(rismt%lfft%ngz))
  END IF
  !
  ! ... set variables
  z0     = 0.5_DP * at(3, 3)
  zright = rismt%lfft%zright
  zleft  = rismt%lfft%zleft
  zstart = rismt%lfft%zleft + rismt%lfft%zoffset
  dz     = rismt%lfft%zstep
  fac1   = e2 * fpi / tpiba2        ! -> Energy/G/G
  fac2   = e2 * fpi * alat / tpiba  ! -> Energy/R/G
  fac3   = e2 * fpi * alat * alat   ! -> Energy*R*R
  !
  ! ... calculate exp(i*gz*z0)
  DO igz = 1, rismt%lfft%ngz
    gz   = rismt%lfft%gz(igz)
    phi  = tpi * gz * z0
    expigz(igz) = CMPLX(COS(phi), SIN(phi), kind=DP)
  END DO
  !
  ! ... mapping sphere Rho(G) -> cyrinder Rho(Gz, Gxy)
  rhogt = C_ZERO
  DO ig = 1, rismt%gvec%ngm
    igz  = rismt%lfft%igtoigz(1, ig)
    igxy = rismt%lfft%igtoigxy(ig)
    rhogt(igz, igxy) = rhog(ig)
    !
    ! ... for gamma_only
    igz = rismt%lfft%igtoigz(2, ig)
    IF (igz > 0) THEN
      rhogt(igz, igxy) = CONJG(rhog(ig))
    END IF
  END DO
  !
  ! ...
  ! ... Hartree potential, when Gxy /= 0
  ! ...
  DO igxy = rismt%lfft%gxystart, rismt%lfft%ngxy
    !
    jgxy = rismt%nrzl * (igxy - 1)
    gxy  = rismt%lfft%gnxy(igxy)
    !
    IF (rismt%lfft%ngz > 0) THEN
      rhogz(:) = rhogt(:, igxy)
    END IF
    !
    ! ... coeff1 means
    !      ----           exp( i*gz*z0)
    !      >    rho(g) * ---------------
    !      ----            i*gz + gxy
    !       gz
    !
    ! ... coeff2 means
    !      ----           exp(-i*gz*z0)
    !      >    rho(g) * ---------------
    !      ----            i*gz + gxy
    !       gz
    !
    ! ... coeff3 means
    !      ----           exp( i*gz*z0)
    !      >    rho(g) * ---------------
    !      ----            i*gz - gxy
    !       gz
    !
    ! ... coeff4 means
    !      ----           exp(-i*gz*z0)
    !      >    rho(g) * ---------------
    !      ----            i*gz - gxy
    !       gz
    !
    coeff1 = C_ZERO
    coeff2 = C_ZERO
    coeff3 = C_ZERO
    coeff4 = C_ZERO
    !
!$omp parallel do default(shared) private(igz, gz) reduction(+:coeff1, coeff2, coeff3, coeff4)
    DO igz = 1, rismt%lfft%ngz
      gz = rismt%lfft%gz(igz)
      coeff1 = coeff1 + rhogz(igz) *       expigz(igz)  / CMPLX( gxy, gz, kind=DP)
      coeff2 = coeff2 + rhogz(igz) * CONJG(expigz(igz)) / CMPLX( gxy, gz, kind=DP)
      coeff3 = coeff3 + rhogz(igz) *       expigz(igz)  / CMPLX(-gxy, gz, kind=DP)
      coeff4 = coeff4 + rhogz(igz) * CONJG(expigz(igz)) / CMPLX(-gxy, gz, kind=DP)
    END DO
!$omp end parallel do
    !
    ! ... when z < zleft, potential is
    !
    !                        [           exp( gxy*(zleft-z0))    ----           exp( i*gz*z0)
    !  exp( gxy*(z-zleft)) * [ (+4pi) * ---------------------- * >    rho(g) * ---------------
    !                        [                2 * gxy            ----            i*gz - gxy
    !                                                             gz
    !
    !                                    exp( gxy*(zleft+z0))    ----           exp(-i*gz*z0)  ]
    !                          (-4pi) * ---------------------- * >    rho(g) * --------------- ]
    !                                         2 * gxy            ----            i*gz - gxy    ]
    !                                                             gz
    !
    vleft(igxy) = vleft(igxy) + fac1 * ( &
    & + (0.5_DP / gxy) * EXP( tpi * gxy * (zleft - z0)) * coeff3 &
    & - (0.5_DP / gxy) * EXP( tpi * gxy * (zleft + z0)) * coeff4 )
    !
    ! ... when zleft <= z < -z0, potential is
    !
    !            exp( gxy*(z-z0))    ----           exp( i*gz*z0)
    !  (+4pi) * ------------------ * >    rho(g) * ---------------
    !               2 * gxy          ----            i*gz - gxy
    !                                 gz
    !
    !            exp( gxy*(z+z0))    ----           exp(-i*gz*z0)
    !  (-4pi) * ------------------ * >    rho(g) * ---------------
    !               2 * gxy          ----            i*gz - gxy
    !                                 gz
    !
!$omp parallel do default(shared) private(iz, z)
    DO iz = 1, (rismt%lfft%izcell_start - 1)
      z = zstart + dz * DBLE(iz - 1)
      vlaue(iz + jgxy) = vlaue(iz + jgxy) + fac1 * ( &
      & + (0.5_DP / gxy) * EXP( tpi * gxy * (z - z0)) * coeff3 &
      & - (0.5_DP / gxy) * EXP( tpi * gxy * (z + z0)) * coeff4 )
    END DO
!$omp end parallel do
    !
    ! ... when -z0 <= z <= z0, potential is
    !
    !            exp( gxy*(z-z0))    ----           exp( i*gz*z0)
    !  (+4pi) * ------------------ * >    rho(g) * ---------------
    !               2 * gxy          ----            i*gz - gxy
    !                                 gz
    !
    !            exp(-gxy*(z+z0))    ----           exp(-i*gz*z0)
    !  (-4pi) * ------------------ * >    rho(g) * ---------------
    !               2 * gxy          ----            i*gz + gxy
    !                                 gz
    !
!$omp parallel do default(shared) private(iz, z)
    DO iz = rismt%lfft%izcell_start, rismt%lfft%izcell_end
      z = zstart + dz * DBLE(iz - 1)
      vlaue(iz + jgxy) = vlaue(iz + jgxy) + fac1 * ( &
      & + (0.5_DP / gxy) * EXP( tpi * gxy * (z - z0)) * coeff3 &
      & - (0.5_DP / gxy) * EXP(-tpi * gxy * (z + z0)) * coeff2 )
    END DO
!$omp end parallel do
    !
    ! ... when z0 < z <= zright, potential is
    !
    !            exp(-gxy*(z-z0))    ----           exp( i*gz*z0)
    !  (+4pi) * ------------------ * >    rho(g) * ---------------
    !               2 * gxy          ----            i*gz + gxy
    !                                 gz
    !
    !            exp(-gxy*(z+z0))    ----           exp(-i*gz*z0)
    !  (-4pi) * ------------------ * >    rho(g) * ---------------
    !               2 * gxy          ----            i*gz + gxy
    !                                 gz
    !
!$omp parallel do default(shared) private(iz, z)
    DO iz = (rismt%lfft%izcell_end + 1), rismt%lfft%nrz
      z = zstart + dz * DBLE(iz - 1)
      vlaue(iz + jgxy) = vlaue(iz + jgxy) + fac1 * ( &
      & + (0.5_DP / gxy) * EXP(-tpi * gxy * (z - z0)) * coeff1 &
      & - (0.5_DP / gxy) * EXP(-tpi * gxy * (z + z0)) * coeff2 )
    END DO
!$omp end parallel do
    !
    ! ... when zright < z, potential is
    !
    !                         [           exp(-gxy*(zright-z0))    ----           exp( i*gz*z0)
    !  exp(-gxy*(z-zright)) * [ (+4pi) * ----------------------- * >    rho(g) * ---------------
    !                         [                2 * gxy             ----            i*gz + gxy
    !                                                               gz
    !
    !                                     exp(-gxy*(zright+z0))    ----           exp(-i*gz*z0)  ]
    !                           (-4pi) * ----------------------- * >    rho(g) * --------------- ]
    !                                          2 * gxy             ----            i*gz + gxy    ]
    !                                                               gz
    !
    vright(igxy) = vright(igxy) + fac1 * ( &
    & + (0.5_DP / gxy) * EXP(-tpi * gxy * (zright - z0)) * coeff1 &
    & - (0.5_DP / gxy) * EXP(-tpi * gxy * (zright + z0)) * coeff2 )
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
    IF (rismt%lfft%ngz > 0) THEN
      rhogz(:) = rhogt(:, igxy)
      rho0 = DBLE(rhogz(rismt%lfft%gzzero))
    END IF
    !
    ! ... cos1 means
    !      ----               2*cos(gz*z0)
    !      >    Im[rho(g)] * --------------
    !      ----                   gz
    !      gz>0
    !
    ! ... sin1 means
    !      ----               2*sin(gz*z0)
    !      >    Re[rho(g)] * --------------
    !      ----                   gz
    !      gz>0
    !
    ! ... cos2 means
    !      ----               2*cos(gz*z0)
    !      >    Re[rho(g)] * --------------
    !      ----                   gz^2
    !      gz>0
    !
    ! ... sin2 means
    !      ----               2*sin(gz*z0)
    !      >    Im[rho(g)] * --------------
    !      ----                   gz^2
    !      gz>0
    !
    cos1 = 0.0_DP
    sin1 = 0.0_DP
    cos2 = 0.0_DP
    sin2 = 0.0_DP
    !
!$omp parallel do default(shared) private(igz, gz) reduction(+:cos1, sin1, cos2, sin2)
    DO igz = (rismt%lfft%gzzero + 1), rismt%lfft%ngz
      gz = rismt%lfft%gz(igz)
      cos1 = cos1 + 2.0_DP * AIMAG(rhogz(igz)) * DBLE( expigz(igz)) / gz
      sin1 = sin1 + 2.0_DP * DBLE( rhogz(igz)) * AIMAG(expigz(igz)) / gz
      cos2 = cos2 + 2.0_DP * DBLE( rhogz(igz)) * DBLE( expigz(igz)) / gz / gz
      sin2 = sin2 + 2.0_DP * AIMAG(rhogz(igz)) * AIMAG(expigz(igz)) / gz / gz
    END DO
!$omp end parallel do
    !
    ! ... when z < -z0, potential is
    !
    !           ----               2*sin(gz*z0)
    !  (-4pi) * >    Im[rho(g)] * --------------
    !           ----                   gz^2
    !           gz>0
    !
    !                ----               2*cos(gz*z0)
    !  (-4pi) * z0 * >    Im[rho(g)] * --------------
    !                ----                   gz
    !                gz>0
    !
    !                ----               2*sin(gz*z0)
    !  (+4pi) * z  * >    Re[rho(g)] * --------------
    !                ----                   gz
    !                gz>0
    !
    !  (+4pi) * z*z0 * rho(0)
    !
!$omp parallel do default(shared) private(iz, z)
    DO iz = 1, (rismt%lfft%izcell_start - 1)
      z = zstart + dz * DBLE(iz - 1)
      vlaue(iz + jgxy) = vlaue(iz + jgxy) + CMPLX( &
      & + fac1 * ( - sin2 )               &
      & + fac2 * ( - z0 * cos1            &
      &            + z  * sin1 )          &
      & + fac3 * ( + z  * z0  * rho0 )    &
      & , 0.0_DP, kind=DP)
    END DO
!$omp end parallel do
    !
    vleft(igxy) = vleft(igxy) + &
    CMPLX(fac2 * sin1 + fac3 * z0 * rho0, -fac1 * sin2 - fac2 * z0 * cos1, kind=DP)
    !
    ! ... when -z0 <= z <= z0, potential is
    !
    !           ----               2*cos(gz*z0)
    !  (-4pi) * >    Re[rho(g)] * --------------
    !           ----                   gz^2
    !           gz>0
    !
    !                ----               2*cos(gz*z0)
    !  (+4pi) * z  * >    Im[rho(g)] * --------------
    !                ----                   gz
    !                gz>0
    !
    !                ----               2*sin(gz*z0)
    !  (-4pi) * z0 * >    Re[rho(g)] * --------------
    !                ----                   gz
    !                gz>0
    !
    !  (-4pi) * z0^2 / 2 * rho(0)
    !
    !  (-4pi) * z^2  / 2 * rho(0)
    !
!$omp parallel do default(shared) private(iz, z)
    DO iz = rismt%lfft%izcell_start, rismt%lfft%izcell_end
      z = zstart + dz * DBLE(iz - 1)
      vlaue(iz + jgxy) = vlaue(iz + jgxy) + CMPLX( &
      & + fac1 * ( - cos2 )                        &
      & + fac2 * ( + z  * cos1                     &
      &            - z0 * sin1 )                   &
      & + fac3 * ( - z0 * z0 * 0.5_DP * rho0       &
      &            - z  * z  * 0.5_DP * rho0 )     &
      & , 0.0_DP, kind=DP)
    END DO
!$omp end parallel do
    !
    ! ... when z0 < z, potential is
    !
    !           ----               2*sin(gz*z0)
    !  (+4pi) * >    Im[rho(g)] * --------------
    !           ----                   gz^2
    !           gz>0
    !
    !                ----               2*cos(gz*z0)
    !  (+4pi) * z0 * >    Im[rho(g)] * --------------
    !                ----                   gz
    !                gz>0
    !
    !                ----               2*sin(gz*z0)
    !  (-4pi) * z  * >    Re[rho(g)] * --------------
    !                ----                   gz
    !                gz>0
    !
    !  (-4pi) * z*z0 * rho(0)
    !
!$omp parallel do default(shared) private(iz, z)
    DO iz = (rismt%lfft%izcell_end + 1), rismt%lfft%nrz
      z = zstart + dz * DBLE(iz - 1)
      vlaue(iz + jgxy) = vlaue(iz + jgxy) + CMPLX( &
      & + fac1 * ( + sin2 )               &
      & + fac2 * ( + z0 * cos1            &
      &            - z  * sin1 )          &
      & + fac3 * ( - z  * z0  * rho0 )    &
      & , 0.0_DP, kind=DP)
    END DO
!$omp end parallel do
    !
    vright(igxy) = vright(igxy) + &
    CMPLX(-fac2 * sin1 - fac3 * z0 * rho0, fac1 * sin2 + fac2 * z0 * cos1, kind=DP)
    !
  END IF
  !
  ! ... deallocate memory
  IF (rismt%lfft%ngz * rismt%lfft%ngxy > 0) THEN
    DEALLOCATE(rhogt)
  END IF
  IF (rismt%lfft%ngz > 0) THEN
    DEALLOCATE(rhogz)
    DEALLOCATE(expigz)
  END IF
  !
  ! ... normally done
  ierr = IERR_RISM_NULL
  !
END SUBROUTINE potential_esm_hartree
!
!---------------------------------------------------------------------------
SUBROUTINE potential_esm_local(rismt, alpha, vlaue, vright, vleft, ierr)
  !---------------------------------------------------------------------------
  !
  ! ... calculate local potential of ESM(BC1), for Laue-RISM.
  ! ... local potential is derived from Gaussian functions:
  ! ...
  ! ...                      1               [   |r - R|^2 ]
  ! ...   rho(r) = -------------------- * exp[- -----------]
  ! ...             pi^(3/2) * alpha^3       [    alpha^2  ]  .
  ! ...
  ! ... Variables:
  ! ...   alpha:  gaussian width (in alat units)
  ! ...   vlaue:  coulomb potential in Laue-rep.
  ! ...   vright: coefficient of coulomb potential when zright < z
  ! ...   vleft:  coefficient of coulomb potential when z < zleft
  !
  USE cell_base, ONLY : at, alat, tpiba
  USE constants, ONLY : pi, tpi, sqrtpi, e2
  USE err_rism,  ONLY : IERR_RISM_NULL, IERR_RISM_INCORRECT_DATA_TYPE
  USE gvect,     ONLY : eigts1, eigts2
  USE ions_base, ONLY : nat, tau, ityp, zv
  USE kinds,     ONLY : DP
  USE rism,      ONLY : rism_type, ITYPE_LAUERISM
  !
  IMPLICIT NONE
  !
  TYPE(rism_type), INTENT(IN)    :: rismt
  REAL(DP),        INTENT(IN)    :: alpha
  COMPLEX(DP),     INTENT(INOUT) :: vlaue(1:*)
  COMPLEX(DP),     INTENT(INOUT) :: vright(1:*)
  COMPLEX(DP),     INTENT(INOUT) :: vleft(1:*)
  INTEGER,         INTENT(OUT)   :: ierr
  !
  INTEGER     :: ia
  INTEGER     :: it
  INTEGER     :: iz
  INTEGER     :: igxy
  INTEGER     :: jgxy
  INTEGER     :: mx, my
  REAL(DP)    :: z
  REAL(DP)    :: za
  REAL(DP)    :: zright
  REAL(DP)    :: zleft
  REAL(DP)    :: zstart
  REAL(DP)    :: dz
  REAL(DP)    :: gxy
  REAL(DP)    :: qa
  REAL(DP)    :: area_xy
  REAL(DP)    :: fac1    ! factors to    -> Energy/R/R/G
  REAL(DP)    :: fac2    ! convert units -> Energy/R
  REAL(DP)    :: rcoeff
  COMPLEX(DP) :: ccoeff
  COMPLEX(DP) :: strf_xy
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
  ! ... set variables
  zright  = rismt%lfft%zright
  zleft   = rismt%lfft%zleft
  zstart  = rismt%lfft%zleft + rismt%lfft%zoffset
  dz      = rismt%lfft%zstep
  area_xy = ABS(at(1, 1) * at(2, 2) - at(1, 2) * at(2, 1))
  fac1    = e2 / alat / alat / tpiba  ! -> Energy/R/R/G
  fac2    = e2 / alat                 ! -> Energy/R
  !
  ! ...
  ! ... local potential, when Gxy /= 0
  ! ...
  DO igxy = rismt%lfft%gxystart, rismt%lfft%ngxy
    !
    jgxy = rismt%nrzl * (igxy - 1)
    gxy  = rismt%lfft%gnxy(igxy)
    mx   = rismt%lfft%millxy(1, igxy)
    my   = rismt%lfft%millxy(2, igxy)
    !
    DO ia = 1, nat
      !
      it = ityp(ia)
      qa = zv(it)
      za = tau(3, ia)
      !
      strf_xy = eigts1(mx, ia) * eigts2(my, ia)
      rcoeff  = -qa * fac1 * pi / area_xy / gxy
      ccoeff  = rcoeff * strf_xy
      !
      ! ... when z < zleft, potential is
      !
      !                         [  2pi*exp(-i(gx*xa+gy*ya))                        ]
      !   exp( gxy*(z-zleft)) * [ ------------------------- * exp( gxy*(zleft-za)) ]
      !                         [          S * gxy                                 ]
      !
      vleft(igxy) = vleft(igxy) + 2.0_DP * ccoeff * EXP( tpi * gxy * (zleft - za))
      !
      ! ... when zleft <= z <= zright, potential is
      !
      !    pi*exp(-i(gx*xa+gy*ya))  [                          gxy*alpha     z-za
      !   ------------------------- [ exp( gxy*(z-za)) * erfc(----------- + -------)
      !            S * gxy          [                              2         alpha
      !
      !                                                        gxy*alpha     z-za    ]
      !                             + exp(-gxy*(z-za)) * erfc(----------- - -------) ]
      !                                                            2         alpha   ]
      !
      ! ... NOTE: to avoid overflows,
      ! ...       exp(var1)*erfc(var2) = exp(var1 + log(erfc(var2))) .
      !
!$omp parallel do default(shared) private(iz, z)
      DO iz = 1, rismt%lfft%nrz
        z = zstart + dz * DBLE(iz - 1)
        vlaue(iz + jgxy) = vlaue(iz + jgxy) + ccoeff * ( &
        &   EXP( tpi * gxy * (z - za) + LOG(erfc(0.5_DP * tpi * gxy * alpha + (z - za) / alpha))) &
        & + EXP(-tpi * gxy * (z - za) + LOG(erfc(0.5_DP * tpi * gxy * alpha - (z - za) / alpha))) )
      END DO
!$omp end parallel do
      !
      ! ... when zright < z, potential is
      !
      !                          [  2pi*exp(-i(gx*xa+gy*ya))                         ]
      !   exp(-gxy*(z-zright)) * [ ------------------------- * exp(-gxy*(zright-za)) ]
      !                          [          S * gxy                                  ]
      !
      vright(igxy) = vright(igxy) + 2.0_DP * ccoeff * EXP(-tpi * gxy * (zright - za))
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
    gxy  = rismt%lfft%gnxy(igxy)
    !
    DO ia = 1, nat
      !
      it = ityp(ia)
      qa = zv(it)
      za = tau(3, ia)
      !
      rcoeff = -qa * fac2 * tpi / area_xy
      ccoeff = CMPLX(rcoeff, 0.0_DP, kind=DP)
      !
      ! ... potential is
      !
      !    2pi  [     alpha           (z-za)^2                   z-za    ]
      !   ----- [ - ---------- exp(- ----------) - (z-za) * erf(-------) ]
      !     S   [    pi^(1/2)         alpha^2                    alpha   ]
      !
!$omp parallel do default(shared) private(iz, z)
      DO iz = 1, rismt%lfft%nrz
        z = zstart + dz * DBLE(iz - 1)
        vlaue(iz + jgxy) = vlaue(iz + jgxy) + ccoeff * ( &
        & - (alpha / sqrtpi) * EXP(-(z - za) * (z - za) / alpha / alpha) &
        & - (z - za) * erf((z - za) / alpha)                           )
      END DO
!$omp end parallel do
      !
      vright(igxy) = vright(igxy) + CMPLX(-rcoeff,  rcoeff * za, kind=DP)
      !
      vleft( igxy) = vleft( igxy) + CMPLX( rcoeff, -rcoeff * za, kind=DP)
      !
    END DO
    !
  END IF
  !
  ! ... normally done
  ierr = IERR_RISM_NULL
  !
END SUBROUTINE potential_esm_local
!
!---------------------------------------------------------------------------
SUBROUTINE charge_esm(rhog, charge)
  !---------------------------------------------------------------------------
  !
  ! ... calculate total charge.
  ! ...
  ! ... Variables:
  ! ...   rhog:   electronic density in G-space
  ! ...   charge: total charge
  !
  USE cell_base, ONLY : omega
  USE gvect,     ONLY : gstart
  USE ions_base, ONLY : nat, ityp, zv
  USE kinds,     ONLY : DP
  USE mp,        ONLY : mp_sum
  USE mp_bands,  ONLY : intra_bgrp_comm
  !
  IMPLICIT NONE
  !
  COMPLEX(DP), INTENT(IN)  :: rhog(1:*)
  REAL(DP),    INTENT(OUT) :: charge
  !
  INTEGER  :: ia
  INTEGER  :: it
  REAL(DP) :: qa
  REAL(DP) :: qion
  REAL(DP) :: qele
  !
  qion = 0.0_DP
  !
  DO ia = 1, nat
    it = ityp(ia)
    qa = zv(it)
    qion = qion + qa
  END DO
  !
  IF (gstart > 1) THEN
    qele = omega * DBLE(rhog(1))
  ELSE
    qele = 0.0_DP
  END IF
  CALL mp_sum(qele, intra_bgrp_comm)
  !
  charge = qion - qele
  !
END SUBROUTINE charge_esm
!
!---------------------------------------------------------------------------
SUBROUTINE check_esm_outside(rismt, ierr)
  !---------------------------------------------------------------------------
  !
  ! ... check whether coulomb potential in [ z < zleft, zright < z ] is finite value or negligible.
  !
  USE cell_base, ONLY : at, alat, tpiba
  USE constants, ONLY : tpi, e2, eps6
  USE err_rism,  ONLY : IERR_RISM_NULL, IERR_RISM_INCORRECT_DATA_TYPE
  USE kinds,     ONLY : DP
  USE rism,      ONLY : rism_type, ITYPE_LAUERISM
  !
  IMPLICIT NONE
  !
  TYPE(rism_type), INTENT(INOUT) :: rismt
  INTEGER,         INTENT(OUT)   :: ierr
  !
  INTEGER     :: igxy
  REAL(DP)    :: z0
  REAL(DP)    :: zright
  REAL(DP)    :: zleft
  REAL(DP)    :: gxy
  REAL(DP)    :: area_xy
  REAL(DP)    :: fac  ! factors to convert units -> Energy/R/R/G
  REAL(DP)    :: eright
  REAL(DP)    :: eleft
  REAL(DP)    :: rcoeff
  !
  ! ... check data type
  IF (rismt%itype /= ITYPE_LAUERISM) THEN
    ierr = IERR_RISM_INCORRECT_DATA_TYPE
    RETURN
  END IF
  !
  IF (rismt%ngxy < rismt%lfft%ngxy) THEN
    ierr = IERR_RISM_INCORRECT_DATA_TYPE
    RETURN
  END IF
  !
  ! ... set variables
  z0      = 0.5_DP * at(3, 3)
  zright  = rismt%lfft%zright
  zleft   = rismt%lfft%zleft
  area_xy = ABS(at(1, 1) * at(2, 2) - at(1, 2) * at(2, 1))
  fac     = e2 / alat / alat / tpiba  ! -> Energy/R/R/G
  !
  ! ... initialize value
  rismt%do_vright = .FALSE.
  rismt%do_vleft  = .FALSE.
  !
  ! ...
  ! ... check when Gxy /= 0
  ! ...
  DO igxy = rismt%lfft%gxystart, rismt%lfft%ngxy
    !
    gxy  = rismt%lfft%gnxy(igxy)
    !
    rcoeff = fac * tpi / area_xy / gxy
    !
    ! ... when z < zleft, check
    !
    !    2pi*exp( gxy*(zleft+z0))
    !   -------------------------
    !            S * gxy
    !
    eleft = rcoeff * EXP( tpi * gxy * (zleft  + z0))
    rismt%do_vleft(igxy) = (ABS(eleft) > eps6)
    !
    ! ... when zright < z, check
    !
    !    2pi*exp(-gxy*(zright-z0))
    !   -------------------------
    !            S * gxy
    !
    eright = rcoeff * EXP(-tpi * gxy * (zright - z0))
    rismt%do_vright(igxy) = (ABS(eright) > eps6)
    !
  END DO
  !
  ! ...
  ! ... check when Gxy = 0
  ! ...
  IF (rismt%lfft%gxystart > 1) THEN
    rismt%do_vright(1) = .TRUE.
    rismt%do_vleft( 1) = .TRUE.
  END IF
  !
  ! ... normally done
  ierr = IERR_RISM_NULL
  !
END SUBROUTINE check_esm_outside
