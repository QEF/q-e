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
SUBROUTINE chempot_lauerism(rismt, ierr)
  !---------------------------------------------------------------------------
  !
  ! ... calculate chemical potential of solvation.
  ! ... although a chemical potential for each site (usol) does not have a meaning,
  ! ... only solvation' energy (esol) has a physical meaning.
  ! ...
  ! ... NOTE: the current version support only Gaussian Fluctuation model.
  !
  USE constants, ONLY : K_BOLTZMANN_RY
  USE kinds,     ONLY : DP
  USE mp,        ONLY : mp_sum
  USE err_rism,  ONLY : IERR_RISM_NULL, IERR_RISM_INCORRECT_DATA_TYPE
  USE rism,      ONLY : rism_type, get_chempot_type, ITYPE_LAUERISM, CHEMPOT_GF
  USE solvmol,   ONLY : get_nuniq_in_solVs, solVs, &
                      & iuniq_to_nsite, iuniq_to_isite, isite_to_isolV, isite_to_iatom
  !
  IMPLICIT NONE
  !
  TYPE(rism_type), INTENT(INOUT) :: rismt
  INTEGER,         INTENT(OUT)   :: ierr
  !
  INTEGER  :: ichempot
  INTEGER  :: nq
  INTEGER  :: iq
  INTEGER  :: iiq
  INTEGER  :: iv
  INTEGER  :: nv
  INTEGER  :: isolV
  INTEGER  :: iatom
  REAL(DP) :: rhov_right
  REAL(DP) :: rhov_left
  REAL(DP) :: qv
  REAL(DP) :: beta
  REAL(DP) :: usol
  !
  ! ... number of sites in solvents
  nq = get_nuniq_in_solVs()
  !
  ! ... check data type
  IF (rismt%itype /= ITYPE_LAUERISM) THEN
    ierr = IERR_RISM_INCORRECT_DATA_TYPE
    RETURN
  END IF
  !
  IF (rismt%mp_site%nsite < nq) THEN
    ierr = IERR_RISM_INCORRECT_DATA_TYPE
    RETURN
  END IF
  !
  IF (rismt%nrzs < rismt%dfft%nr3) THEN
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
  ! ... if no data, return as normally done
  IF (rismt%nsite < 1) THEN
    GOTO 100
  END IF
  !
  ! ... set variables
  ichempot = get_chempot_type(rismt)
  beta = 1.0_DP / K_BOLTZMANN_RY / rismt%temp
  !
  ! ... chemical potential for each site
  IF (rismt%nsite > 0) THEN
    rismt%usol    = 0.0_DP
    rismt%usol_GF = 0.0_DP
  END IF
  !
  DO iq = rismt%mp_site%isite_start, rismt%mp_site%isite_end
    iiq        = iq - rismt%mp_site%isite_start + 1
    iv         = iuniq_to_isite(1, iq)
    nv         = iuniq_to_nsite(iq)
    isolV      = isite_to_isolV(iv)
    iatom      = isite_to_iatom(iv)
    rhov_right = DBLE(nv) * solVs(isolV)%density
    rhov_left  = DBLE(nv) * solVs(isolV)%subdensity
    qv         = solVs(isolV)%charge(iatom)
    !
    CALL chempot_laue_for_a_site(rismt, ichempot, &
                               & iiq, rhov_right, rhov_left, qv, beta, usol)
    !
    rismt%usol(iiq) = usol
    !
    CALL chempot_laue_for_a_site(rismt, CHEMPOT_GF, &
                               & iiq, rhov_right, rhov_left, qv, beta, usol)
    !
    rismt%usol_GF(iiq) = usol
    !
  END DO
  !
  CALL mp_sum(rismt%usol,    rismt%mp_site%intra_sitg_comm)
  CALL mp_sum(rismt%usol_GF, rismt%mp_site%intra_sitg_comm)
  !
  ! ... normally done
100 CONTINUE
  ierr = IERR_RISM_NULL
  !
END SUBROUTINE chempot_lauerism
!
!---------------------------------------------------------------------------
SUBROUTINE chempot_laue_for_a_site(rismt, ichempot, isite, rhov1, rhov2, qv, beta, usol)
  !---------------------------------------------------------------------------
  !
  ! ... calculate chemical potential for a site
  !
  USE kinds, ONLY : DP
  USE rism,  ONLY : rism_type, CHEMPOT_HNC, CHEMPOT_KH, CHEMPOT_GF
  !
  IMPLICIT NONE
  !
  TYPE(rism_type), INTENT(IN)  :: rismt
  INTEGER,         INTENT(IN)  :: ichempot
  INTEGER,         INTENT(IN)  :: isite
  REAL(DP),        INTENT(IN)  :: rhov1
  REAL(DP),        INTENT(IN)  :: rhov2
  REAL(DP),        INTENT(IN)  :: qv
  REAL(DP),        INTENT(IN)  :: beta
  REAL(DP),        INTENT(OUT) :: usol
  !
  IF (ichempot == CHEMPOT_HNC) THEN
    CALL chempot_laue_HNC_x(rismt, isite, rhov1, rhov2, qv, beta, usol)
    !
  ELSE IF (ichempot == CHEMPOT_KH) THEN
    CALL chempot_laue_KH_x(rismt, isite, rhov1, rhov2, qv, beta, usol)
    !
  ELSE IF (ichempot == CHEMPOT_GF) THEN
    CALL chempot_laue_GF_x(rismt, isite, rhov1, rhov2, qv, beta, usol)
    !
  ELSE
    usol = 0.0_DP
  END IF
  !
  usol = usol / beta
  !
END SUBROUTINE chempot_laue_for_a_site
!
!---------------------------------------------------------------------------
SUBROUTINE chempot_laue_HNC_x(rismt, isite, rhov1, rhov2, qv, beta, usol)
  !---------------------------------------------------------------------------
  !
  ! ... HyperNetted-Chain model
  ! ... (J.P.Hansen et al., Theory of simple liquids. Academic Press, London, 1990.)
  ! ...
  ! ...            /    [  1                   1              ]
  ! ...   kB * T * | dr [ --- h(r)^2 - c(r) - --- h(r) * c(r) ]
  ! ...            /    [  2                   2              ]
  !
  USE kinds, ONLY : DP
  USE rism,  ONLY : rism_type
  !
  IMPLICIT NONE
  !
  TYPE(rism_type), INTENT(IN)  :: rismt
  INTEGER,         INTENT(IN)  :: isite
  REAL(DP),        INTENT(IN)  :: rhov1
  REAL(DP),        INTENT(IN)  :: rhov2
  REAL(DP),        INTENT(IN)  :: qv
  REAL(DP),        INTENT(IN)  :: beta
  REAL(DP),        INTENT(OUT) :: usol
  !
  usol = 0.0_DP
  !
END SUBROUTINE chempot_laue_HNC_x
!
!---------------------------------------------------------------------------
SUBROUTINE chempot_laue_KH_x(rismt, isite, rhov1, rhov2, qv, beta, usol)
  !---------------------------------------------------------------------------
  !
  ! ... Kovalenko and Hirata's model
  ! ... (A.Kovalenko, F.Hirata, J. Chem. Phys. 1999, 110, 10095-10112)
  ! ...
  ! ...            /    [  1                                  1              ]
  ! ...   kB * T * | dr [ --- h(r)^2 * Theta(-h(r)) - c(r) - --- h(r) * c(r) ]
  ! ...            /    [  2                                  2              ]
  !
  USE kinds, ONLY : DP
  USE rism,  ONLY : rism_type
  !
  IMPLICIT NONE
  !
  TYPE(rism_type), INTENT(IN)  :: rismt
  INTEGER,         INTENT(IN)  :: isite
  REAL(DP),        INTENT(IN)  :: rhov1
  REAL(DP),        INTENT(IN)  :: rhov2
  REAL(DP),        INTENT(IN)  :: qv
  REAL(DP),        INTENT(IN)  :: beta
  REAL(DP),        INTENT(OUT) :: usol
  !
  usol = 0.0_DP
  !
END SUBROUTINE chempot_laue_KH_x
!
!---------------------------------------------------------------------------
SUBROUTINE chempot_laue_GF_x(rismt, isite, rhov1, rhov2, qv, beta, usol)
  !---------------------------------------------------------------------------
  !
  ! ... Gaussian Fluctuation model
  ! ... (D.Chandler et al., J. Chem. Phys. 1984, 81, 1975-1982)
  ! ...
  ! ...            /    [           1              ]
  ! ...   kB * T * | dr [ - c(r) - --- h(r) * c(r) ]
  ! ...            /    [           2              ]
  !
  USE cell_base,     ONLY : at, alat
  USE control_flags, ONLY : gamma_only
  USE kinds,         ONLY : DP
  USE rism,          ONLY : rism_type
  !
  IMPLICIT NONE
  !
  TYPE(rism_type), INTENT(IN)  :: rismt
  INTEGER,         INTENT(IN)  :: isite
  REAL(DP),        INTENT(IN)  :: rhov1
  REAL(DP),        INTENT(IN)  :: rhov2
  REAL(DP),        INTENT(IN)  :: qv
  REAL(DP),        INTENT(IN)  :: beta
  REAL(DP),        INTENT(OUT) :: usol
  !
  INTEGER  :: iz
  INTEGER  :: jz
  INTEGER  :: igxy
  INTEGER  :: jgxy
  INTEGER  :: kgxy
  REAL(DP) :: usol1
  REAL(DP) :: usol2
  REAL(DP) :: usol_
  REAL(DP) :: area_xy
  REAL(DP) :: zstep
  REAL(DP) :: csr, csi
  REAL(DP) :: clr, cli
  REAL(DP) :: cr, ci
  REAL(DP) :: hr, hi
  !
  usol1 = 0.0_DP
  usol2 = 0.0_DP
  !
  ! ... in case Gxy = 0
  IF (rismt%lfft%gxystart > 1) THEN
    !
    usol_ = 0.0_DP
!$omp parallel do default(shared) private(iz, csr, clr, cr, hr) reduction(+:usol_)
    DO iz = 1, rismt%lfft%izleft_gedge
      csr = rismt%csdg0(iz, isite)
      clr = -beta * qv * DBLE(rismt%vlgz(iz))
      cr  = csr + clr
      hr  = DBLE(rismt%hsgz(iz, isite) + rismt%hlgz(iz, isite))
      usol_ = usol_ - rhov2 * (cr + 0.5_DP * hr * cr)
    END DO
!$omp end parallel do
    usol1 = usol1 + usol_
    !
    usol_ = 0.0_DP
!$omp parallel do default(shared) private(iz, csr, clr, cr, hr) reduction(+:usol_)
    DO iz = rismt%lfft%izright_gedge, rismt%lfft%nrz
      csr = rismt%csdg0(iz, isite)
      clr = -beta * qv * DBLE(rismt%vlgz(iz))
      cr  = csr + clr
      hr  = DBLE(rismt%hsgz(iz, isite) + rismt%hlgz(iz, isite))
      usol_ = usol_ - rhov1 * (cr + 0.5_DP * hr * cr)
    END DO
!$omp end parallel do
    usol1 = usol1 + usol_
    !
  END IF
  !
  ! ... in case Gxy /= 0
  DO igxy = rismt%lfft%gxystart, rismt%lfft%ngxy
    jgxy = rismt%nrzs * (igxy - 1)
    kgxy = rismt%nrzl * (igxy - 1)
    !
    usol_ = 0.0_DP
!$omp parallel do default(shared) private(iz, jz, csr, csi, clr, cli, cr, ci, hr, hi) &
!$omp                             reduction(+:usol_)
    DO iz = 1, rismt%lfft%izleft_gedge
      jz = iz - rismt%lfft%izcell_start + 1
      IF (jz > 0) THEN
        csr = DBLE( rismt%csgz(jz + jgxy, isite))
        csi = AIMAG(rismt%csgz(jz + jgxy, isite))
      ELSE
        csr = 0.0_DP
        csi = 0.0_DP
      END IF
      clr = -beta * qv * DBLE( rismt%vlgz(iz + kgxy))
      cli = -beta * qv * AIMAG(rismt%vlgz(iz + kgxy))
      cr  = csr + clr
      ci  = csi + cli
      hr  = DBLE( rismt%hsgz(iz + kgxy, isite) + rismt%hlgz(iz + kgxy, isite))
      hi  = AIMAG(rismt%hsgz(iz + kgxy, isite) + rismt%hlgz(iz + kgxy, isite))
      usol_ = usol_ - rhov2 * 0.5_DP * (hr * cr + hi * ci)
    END DO
!$omp end parallel do
    usol2 = usol2 + usol_
    !
    usol_ = 0.0_DP
!$omp parallel do default(shared) private(iz, jz, csr, csi, clr, cli, cr, ci, hr, hi) &
!$omp                             reduction(+:usol_)
    DO iz = rismt%lfft%izright_gedge, rismt%lfft%nrz
      jz = iz - rismt%lfft%izcell_start + 1
      IF (jz <= rismt%dfft%nr3) THEN
        csr = DBLE( rismt%csgz(jz + jgxy, isite))
        csi = AIMAG(rismt%csgz(jz + jgxy, isite))
      ELSE
        csr = 0.0_DP
        csi = 0.0_DP
      END IF
      clr = -beta * qv * DBLE( rismt%vlgz(iz + kgxy))
      cli = -beta * qv * AIMAG(rismt%vlgz(iz + kgxy))
      cr  = csr + clr
      ci  = csi + cli
      hr  = DBLE( rismt%hsgz(iz + kgxy, isite) + rismt%hlgz(iz + kgxy, isite))
      hi  = AIMAG(rismt%hsgz(iz + kgxy, isite) + rismt%hlgz(iz + kgxy, isite))
      usol_ = usol_ - rhov1 * 0.5_DP * (hr * cr + hi * ci)
    END DO
!$omp end parallel do
    usol2 = usol2 + usol_
    !
  END DO
  !
  IF (gamma_only) THEN
    usol2 = usol2 * 2.0_DP
  END IF
  !
  area_xy = alat * alat * ABS(at(1, 1) * at(2, 2) - at(1, 2) * at(2, 1))
  zstep   = alat * rismt%lfft%zstep
  usol    = (usol1 + usol2) * area_xy * zstep
  !
END SUBROUTINE chempot_laue_GF_x
