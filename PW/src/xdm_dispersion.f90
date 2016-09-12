! Module for the calculation of the XDM dispersion correction. See:
!   A. Otero de la Roza and E. R. Johnson, J. Chem. Phys. 136 (2012) 174109 and 138, 204109 (2013).
!   A. D. Becke and E. R. Johnson, J. Chem. Phys. 127, 154108 (2007) and references therein.

! Copyright (C) 2013 A. Otero-de-la-Roza and E. R. Johnson, University of California-Merced.
! Copyright (C) 2013 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
module xdm_module
  USE kinds, ONLY: DP
  USE constants, ONLY: bohr_radius_angs, pi, fpi
  IMPLICIT NONE

  PRIVATE
  PUBLIC :: a1i, a2i    ! the damping function coefficients (real_dp)
  PUBLIC :: init_xdm    ! initialize XDM: calculate atomic volumes, radial densities,...
  PUBLIC :: energy_xdm  ! compute the XDM dispersion energy and derivatives 
  PUBLIC :: force_xdm   ! fetch the forces calculated by energy_xdm
  PUBLIC :: stress_xdm  ! fetch the stresses calculated by energy_xdm
  PUBLIC :: cleanup_xdm ! deallocate arrays

  ! atomic environments
  INTEGER :: nenv
  REAL(DP), ALLOCATABLE :: xenv(:,:)
  INTEGER, ALLOCATABLE :: ienv(:), lvec(:,:)
  INTEGER :: nvec

  ! moments, polarizabilities, radii, dispersion coefficients
  REAL(DP), ALLOCATABLE :: alpha(:), ml(:,:)
  REAL(DP), ALLOCATABLE :: cx(:,:,:), rvdw(:,:)
  REAL(DP) :: maxc6

  ! energies, forces and stresses
  REAL(DP) :: esave = 0._DP
  REAL(DP) :: esaveold = 0._DP
  REAL(DP), ALLOCATABLE :: fsave(:,:), ssave(:,:)

  ! have moments been computed before?
  LOGICAL :: saved = .FALSE.

  ! a1 and a2 coefficients, with defaults for pw86pbe
  REAL(DP) :: a1i = 0.6836_DP
  REAL(DP) :: a2i = 1.5045_DP

  ! radial atomic densities
  REAL(DP), ALLOCATABLE :: rfree(:,:), w2free(:,:), rmaxg2(:)
  REAL(DP), ALLOCATABLE :: rcore(:,:), w2core(:,:), rmaxcore2(:)

  ! free volumes
  REAL(DP), ALLOCATABLE :: afree(:)

  ! free atomic polarizabilities from CRC handbook, 88th ed. (ang^3->bohr^3)
  REAL(DP), PARAMETER :: alpha_free(1:102) = (/0.6668_DP, 0.2051_DP, 24.3300_DP, 5.6000_DP,&
     3.0300_DP, 1.7600_DP, 1.1000_DP, 0.8020_DP, 0.5570_DP, 0.3956_DP, 24.1100_DP, 10.6000_DP,&
     6.8000_DP, 5.3800_DP, 3.6300_DP, 2.9000_DP, 2.1800_DP, 1.6411_DP, 43.4000_DP, 22.8000_DP,&
     17.8000_DP, 14.6000_DP, 12.4000_DP, 11.6000_DP, 9.4000_DP, 8.4000_DP, 7.5000_DP, 6.8000_DP,&
     6.2000_DP, 5.7500_DP, 8.1200_DP, 6.0700_DP, 4.3100_DP, 3.7700_DP, 3.0500_DP, 2.4844_DP,&
     47.3000_DP, 27.6000_DP, 22.7000_DP, 17.9000_DP, 15.7000_DP, 12.8000_DP, 11.4000_DP, 9.6000_DP,&
     8.6000_DP, 4.8000_DP, 7.2000_DP, 7.3600_DP, 10.2000_DP, 7.7000_DP, 6.6000_DP, 5.5000_DP,&
     5.3500_DP, 4.0440_DP, 59.4200_DP, 39.7000_DP, 31.1000_DP, 29.6000_DP, 28.2000_DP, 31.4000_DP,&
     30.1000_DP, 28.8000_DP, 27.7000_DP, 23.5000_DP, 25.5000_DP, 24.5000_DP, 23.6000_DP, 22.7000_DP,&
     21.8000_DP, 21.0000_DP, 21.9000_DP, 16.2000_DP, 13.1000_DP, 11.1000_DP, 9.7000_DP, 8.5000_DP,&
     7.6000_DP, 6.5000_DP, 5.8000_DP, 5.0200_DP, 7.6000_DP, 6.8000_DP, 7.4000_DP, 6.8000_DP, 6.0000_DP,&
     5.3000_DP, 48.6000_DP, 38.3000_DP, 32.1000_DP, 32.1000_DP, 25.4000_DP, 24.9000_DP, 24.8000_DP,&
     24.5000_DP, 23.3000_DP, 23.0000_DP, 22.7000_DP, 20.5000_DP, 19.7000_DP, 23.8000_DP, 18.2000_DP,&
     17.5000_DP /) / bohr_radius_angs**3

  ! factorials and lambda_l1,l2
  REAL(DP), PARAMETER :: fact(0:8) = REAL((/1,1,2,6,24,120,720,5040,40320/),DP)

CONTAINS

  SUBROUTINE init_xdm()
    ! Initialize storage arrays, calculate the atomic and core radial densities
    ! and integrate the free volumes.
    USE ions_base, ONLY: nat
    USE uspp_param, ONLY : upf
    USE ions_base, ONLY: ntyp => nsp
    USE atom, ONLY : rgrid, msh
    USE splinelib, ONLY : spline

    INTEGER :: i, j, ialloc, nn
    REAL(DP), ALLOCATABLE :: d1y(:), d2y(:)

    IF ( .NOT. ALL (upf(1:ntyp)%tpawp) ) &
       CALL errore("init_xdm","XDM only implemented for PAW",1)

    ! allocate c6, etc.
    ALLOCATE(cx(nat,nat,2:4),rvdw(nat,nat),STAT=ialloc)
    IF (ialloc /= 0) CALL alloc_failed("c6, c8, c10, rvdw")
    ALLOCATE(alpha(nat),ml(3,nat),STAT=ialloc)
    IF (ialloc /= 0) CALL alloc_failed("ml, alpha")
    ALLOCATE(fsave(3,nat),STAT=ialloc)
    IF (ialloc /= 0) CALL alloc_failed("fsave")
    ALLOCATE(ssave(3,3),STAT=ialloc)
    IF (ialloc /= 0) CALL alloc_failed("ssave")

    ! free atomic and core densities
    nn = 0
    DO i = 1, ntyp
       nn = MAX(nn,msh(i))
    END DO
    ALLOCATE(rfree(nn,ntyp),w2free(nn,ntyp),d1y(nn),d2y(nn),rmaxg2(ntyp),STAT=ialloc)
    IF (ialloc /= 0) CALL alloc_failed("rfree")
    ALLOCATE(rcore(nn,ntyp),w2core(nn,ntyp),rmaxcore2(ntyp),STAT=ialloc)
    IF (ialloc /= 0) CALL alloc_failed("rcore")
    DO i = 1, ntyp
       nn = msh(i)
       rfree(1:nn,i) = upf(i)%rho_at(1:nn) / (fpi*rgrid(i)%r(1:nn)**2) + upf(i)%paw%ae_rho_atc(1:nn)
       CALL radial_gradient(rfree(1:nn,i),d1y(1:nn),rgrid(i)%r(1:nn),nn,1)
       CALL radial_gradient(d1y(1:nn),d2y(1:nn),rgrid(i)%r(1:nn),nn,1)
       CALL spline(rgrid(i)%r(1:nn),rfree(1:nn,i),d1y(1),d2y(1),w2free(1:nn,i))
       rmaxg2(i) = rgrid(i)%r(nn)**2

       rcore(1:nn,i) = upf(i)%paw%ae_rho_atc(1:nn)
       CALL radial_gradient(rcore(1:nn,i),d1y(1:nn),rgrid(i)%r(1:nn),nn,1)
       CALL radial_gradient(d1y(1:nn),d2y(1:nn),rgrid(i)%r(1:nn),nn,1)
       CALL spline(rgrid(i)%r(1:nn),rcore(1:nn,i),d1y(1),d2y(1),w2core(1:nn,i))
       if (rcore(1,i) > 1e-8_DP) then
          DO j = nn, 1, -1
             IF (rcore(j,i) > 1e-8_DP) EXIT
          END DO
       else
          j = 1
       end if
       rmaxcore2(i) = rgrid(i)%r(j)**2
    END DO

    ! free volumes
    ALLOCATE(afree(ntyp))
    DO i = 1, ntyp
       nn = msh(i)
       d1y = rfree(1:nn,i) * rgrid(i)%r(1:nn)**5 * fpi
       CALL simpson(nn,d1y,rgrid(i)%rab(1:nn),afree(i))
    END DO
    DEALLOCATE(d1y,d2y)

  END SUBROUTINE init_xdm

  SUBROUTINE cleanup_xdm()
    ! Free all the allocated arrays.

    IF (ALLOCATED(rvdw)) DEALLOCATE(rvdw)
    IF (ALLOCATED(cx)) DEALLOCATE(cx)
    IF (ALLOCATED(alpha)) DEALLOCATE(alpha)
    IF (ALLOCATED(ml)) DEALLOCATE(ml)
    IF (ALLOCATED(fsave)) DEALLOCATE(fsave)
    IF (ALLOCATED(ssave)) DEALLOCATE(ssave)
    IF (ALLOCATED(rfree)) DEALLOCATE(rfree)
    IF (ALLOCATED(w2free)) DEALLOCATE(w2free)
    IF (ALLOCATED(rmaxg2)) DEALLOCATE(rmaxg2)
    IF (ALLOCATED(rcore)) DEALLOCATE(rcore)
    IF (ALLOCATED(w2core)) DEALLOCATE(w2core)
    IF (ALLOCATED(rmaxcore2)) DEALLOCATE(rmaxcore2)
    IF (ALLOCATED(afree)) DEALLOCATE(afree)
    IF (ALLOCATED(xenv)) DEALLOCATE(xenv)
    IF (ALLOCATED(ienv)) DEALLOCATE(ienv)
    IF (ALLOCATED(lvec)) DEALLOCATE(lvec)

  END SUBROUTINE cleanup_xdm

  FUNCTION energy_xdm() RESULT(evdw)
    ! Calculate the XDM dispersion energy correction, forces (Cx coefficients are assumed constant)
    ! and stresses using the electron density and the kinetic energy density to obtain
    ! the dispersion coefficients. The computed coefficients are saved for geometry optimization 
    ! runs. In addition, forces and stresses are saved for subsequent calls to force_xdm
    ! and stress_xdm.
    USE control_flags, ONLY: lbfgs, lmd
    USE scf, ONLY: rho
    USE io_global, ONLY: stdout, ionode
    USE fft_base, ONLY : dfftp
    USE cell_base, ONLY : at, alat, omega
    USE ions_base, ONLY: nat, tau, atm, ityp, ntyp => nsp
    USE constants, ONLY: au_gpa
    USE lsda_mod, ONLY: nspin
    USE atom, ONLY: msh, rgrid
    USE splinelib, ONLY : splint
    USE mp_images, ONLY : me_image, nproc_image, intra_image_comm
    USE mp_pools,  ONLY : me_pool
    USE mp, ONLY : mp_sum

    REAL(DP) :: evdw

    ! energy cutoff for max. interaction distance
    REAL(DP), PARAMETER :: ecut = 1e-11_DP

    INTEGER :: ialloc
    INTEGER :: i, iat, n, ix, iy, iz, j, jj
    REAL(DP), ALLOCATABLE :: gaux(:,:), ggaux(:,:,:), rhoat(:), rhocor(:), rhoae(:)
    REAL(DP), ALLOCATABLE :: lapr(:), gmod(:), avol(:), b(:)
    REAL(DP) :: taus, rhos, ds, qs, rhs, xroot, xshift, xold, expx, gx, fx, ffx
    REAL(DP) :: grho, lap, rhot, rhofree, db2, ri2, rhosf, rhoaf, rc
    REAL(DP) :: x(3), wei, weic, db, ri, atb(3,3), taub(3)
    REAL(DP) :: xij(3), ehadd(6:10), eat, ee
    INTEGER :: l1, l2, ll, m1, m2
    LOGICAL :: docalc
    REAL(DP) :: a1, a2, rmax, rmax2, den, den2
    REAL(DP) :: dij2
    REAL(DP) :: rvdwx, dijx, dijxm2, fxx, cn0
    INTEGER :: i3, nn
    REAL(DP) :: for(3,nat), sigma(3,3), sat(3,3)
    INTEGER :: resto, divid, first, last, it
    INTEGER :: idx0, idx, ispin
    INTEGER, EXTERNAL :: atomic_number

    real*8 :: iix, iiy, iiz

    ! initialize
    IF (nspin > 2) CALL errore('energy_xdm','nspin > 2 not implemented',1)
    evdw = 0._DP
    fsave = 0._DP
    ssave = 0._DP
    atb = at * alat

    ! do we need to recalculate the coefficients?
    docalc = .NOT.saved .OR. .NOT.(lbfgs .OR. lmd)

    ! Define damping coefficients
    a1 = a1i
    a2 = a2i / bohr_radius_angs
    IF (ionode) THEN
       WRITE (stdout,'(/"* XDM dispersion")')
       WRITE (stdout,'("  a1 = ",F12.6)') a1
       WRITE (stdout,'("  a2 (ang) = ",F12.6)') a2i
       WRITE (stdout,'("  a2 (bohr) = ",F12.6)') a2
    END IF

    ! calculate the interaction coefficients
    IF (docalc) THEN

       ! set up the atomic environment for densities
       rmax = SQRT(MAXVAL(rmaxg2))
       CALL set_environ(rmax)

       ! total and core promolecular density
       ALLOCATE(rhoat(dfftp%nnr),rhocor(dfftp%nnr),STAT=ialloc)
       IF (ialloc /= 0) CALL alloc_failed("rhoat/rhocor")
       CALL promolecular_rho(rhoat,rhocor)

       ! all-electron density
       ALLOCATE(rhoae(dfftp%nnr),STAT=ialloc)
       IF (ialloc /= 0) CALL alloc_failed("rhoae")
       CALL PAW_make_ae_charge_xdm(rho,rhoae)
       rhoae = (rhoae + rhocor) / REAL(nspin,DP)

       ! don't need the core anymore
       DEALLOCATE(rhocor)

       ! define the starting index for each processor
       idx0 = dfftp%nr1x * dfftp%nr2x * dfftp%ipp(me_pool+1)

       ! allocate arrays and initialize
       ALLOCATE(b(dfftp%nnr),STAT=ialloc)
       if (ialloc /= 0) CALL alloc_failed("b")
       ALLOCATE(avol(nat),STAT=ialloc)
       IF (ialloc /= 0) CALL alloc_failed("avol")
       ALLOCATE(lapr(dfftp%nnr),STAT=ialloc)
       IF (ialloc /= 0) CALL alloc_failed("lapr")
       ALLOCATE(gmod(dfftp%nnr),STAT=ialloc)
       IF (ialloc /= 0) CALL alloc_failed("gmod")
       avol = 0._DP
       ml = 0._DP
       b = 0._DP

       ! loop over spins
       DO ispin = 1, nspin
          ALLOCATE(gaux(3,dfftp%nnr),ggaux(3,3,dfftp%nnr),STAT=ialloc)
          IF (ialloc /= 0) CALL alloc_failed("gaux, ggaux")

          ! valence gradient and laplacian
          CALL external_ggradient(rho%of_r(:,ispin),gaux,ggaux)
          lapr = ggaux(1,1,:) + ggaux(2,2,:) + ggaux(3,3,:)
          DEALLOCATE(ggaux)
          gmod = sqrt(gaux(1,:)**2 + gaux(2,:)**2 + gaux(3,:)**2)
          DEALLOCATE(gaux)

          ! calculate b on the real-space grid
          DO i = 1, dfftp%nnr
             IF (rho%of_r(i,ispin) < 1e-14_DP) CYCLE
             rhot = MAX(rho%of_r(i,ispin),1e-14_DP)
             rhos = rhot / 2._DP
             grho = gmod(i)
             lap = lapr(i)
             taus = rho%kin_r(i,ispin)
             IF (nspin > 1) THEN
                rhos = rhot
                rhot = MAX(rho%of_r(i,1)+rho%of_r(i,2),1e-14_DP)
             ELSE
                grho = grho / 2._DP
                lap = lap / 2._DP
                taus = taus / 2._DP
             END IF

             ds = taus - 0.25_DP * grho**2 / rhos
             qs = 1._DP/6._DP * (lap - 2._DP * ds)
             rhs = 2._DP/3._DP * pi**(2._DP/3._DP) * (rhos)**(5._DP/3._DP) / qs

             ! newton seed
             IF (rhs > 0._DP) THEN
                xroot = 3._DP
                xshift = 1._DP
                DO WHILE ((xroot * EXP(-2._DP*xroot/3._DP)) / (xroot - 2._DP) < rhs)
                   xshift = xshift * 0.1_DP
                   xroot = 2._DP + xshift
                END DO
             ELSE
                xroot = 1._DP
                xshift = 1._DP
                DO WHILE ( (xroot * EXP(-2._DP*xroot/3._DP)) / (xroot - 2._DP) > rhs)
                   xshift = xshift * 0.1_DP
                   xroot = 2._DP - xshift
                END DO
             END IF

             ! do newton
             xold = 2._DP
             DO WHILE (ABS(xroot - xold) > 1e-10_DP)
                xold = xroot
                expx = EXP(-2._DP * xroot / 3._DP)
                gx = (xroot * expx) / (xroot - 2._DP)
                fx = gx - rhs
                ffx = gx * (1._DP / xroot - 2._DP/3._DP - 1._DP / (xroot - 2._DP))
                xroot = xroot - fx / ffx
             END DO
             b(i) = xroot * (EXP(-xroot) / (8._DP*pi*rhos))**(1._DP/3._DP)
          END DO

          ! integrate atomic volumes and moments
          DO iat = 1, nat
             it = ityp(iat)
             nn = msh(it)
             taub = tau(:,iat) * alat
             DO n = 1, dfftp%nr1x*dfftp%nr2x * dfftp%npl
                idx = idx0 + n - 1

                iz = idx / (dfftp%nr1x*dfftp%nr2x)
                idx = idx - (dfftp%nr1x*dfftp%nr2x)*iz
                iy = idx / dfftp%nr1x
                idx = idx - dfftp%nr1x*iy
                ix = idx
                iix = ix / REAL(dfftp%nr1,DP)
                iiy = iy / REAL(dfftp%nr2,DP)
                iiz = iz / REAL(dfftp%nr3,DP)

                rhosf = rho%of_r(n,ispin) / rhoat(n)
                rhoaf = rhoae(n) / rhoat(n)

                DO ll = 1, nvec
                   x = (lvec(1,ll) + iix) * atb(:,1) + (lvec(2,ll) + iiy) * atb(:,2) + (lvec(3,ll) + iiz) * atb(:,3) - taub
                   ri = x(1)*x(1) + x(2)*x(2) + x(3)*x(3)
                   IF (ri > rmaxg2(it)) CYCLE
                   ri = SQRT(ri)

                   rhofree = splint(rgrid(it)%r(1:nn),rfree(1:nn,it),w2free(1:nn,it),ri)
                   wei = rhofree * rhosf
                   db = MAX(ri-b(n),0._DP)
                   ri2 = 1._DP
                   db2 = 1._DP
                   DO i = 1, 3
                      ri2 = ri2 * ri
                      db2 = db2 * db
                      ml(i,iat) = ml(i,iat) + wei * (ri2 - db2)**2
                   END DO
                   weic = rhofree * rhoaf
                   avol(iat) = avol(iat) + weic * ri2
                END DO ! ll
             END DO ! n
          END DO ! iat
       END DO ! ispin
#if defined(__MPI)
       CALL mp_sum(avol,intra_image_comm)
       CALL mp_sum(ml,intra_image_comm)
#endif
       avol = avol * omega / (dfftp%nr1*dfftp%nr2*dfftp%nr3)
       ml = ml * omega / (dfftp%nr1*dfftp%nr2*dfftp%nr3)

       ! deallocate stuff
       IF (ALLOCATED(b)) DEALLOCATE(b)
       IF (ALLOCATED(rhoat)) DEALLOCATE(rhoat)
       IF (ALLOCATED(rhoae)) DEALLOCATE(rhoae)
       IF (ALLOCATED(lapr)) DEALLOCATE(lapr)
       IF (ALLOCATED(gmod)) DEALLOCATE(gmod)

       ! atom-in-molecule polarizabilities
       DO iat = 1, nat
          it = ityp(iat)
          alpha(iat) = MIN(avol(iat) / afree(it),1._DP) * alpha_free(atomic_number(atm(it)))
       END DO

       ! output the volumes and moments
       IF (ionode) THEN
          WRITE (stdout,*)
          WRITE (stdout,'("+ Volumes and moments")')
          WRITE (stdout,'("# All results in atomic units (Hartree,bohr)")')
          WRITE (stdout,'("# i        V             Vfree           M1             M2             M3")')
          DO iat = 1, nat
             it = ityp(iat)
             WRITE (stdout,'(I3,1p,5(1X,E14.6))') iat, avol(iat), afree(it), ml(1:3,iat)
          END DO
          WRITE (stdout,*)
       END IF

       ! calculate dispersion coefficients and rvdw
       IF (ionode) THEN
          WRITE (stdout,'("+ Dispersion coefficients")')
          WRITE (stdout,'("# All results in atomic units (Hartree,bohr).")')
          WRITE (stdout,'("# i   j      C6             C8             C10            Rc            Rvdw")')
       END IF

       ! critical radii, compute c6, c8, c10, rvdw
       maxc6 = -1._DP
       DO i = 1, nat
          DO j = 1, i
             cx(i,j,2) = alpha(i)*alpha(j)*ml(1,i)*ml(1,j) / (ml(1,i)*alpha(j) + ml(1,j)*alpha(i))
             maxc6 = MAX(cx(i,j,2),maxc6)
             cx(j,i,2) = cx(i,j,2)
             cx(i,j,3) = 3._DP/2._DP * (alpha(i)*alpha(j)*(ml(1,i)*ml(2,j)+ml(2,i)*ml(1,j))) /&
                (ml(1,i)*alpha(j)+ml(1,j)*alpha(i))
             cx(j,i,3) = cx(i,j,3)
             cx(i,j,4) = 2 * alpha(i)*alpha(j) * (ml(1,i)*ml(3,j) + ml(3,i)*ml(1,j)) /&
                (ml(1,i)*alpha(j) + ml(1,j)*alpha(i)) + 21._DP/5._DP * alpha(i)*alpha(j)*&
                ml(2,i)*ml(2,j) / (alpha(j)*ml(1,i)+alpha(i)*ml(1,j))
             cx(j,i,4) = cx(i,j,4)
             rc = (SQRT(cx(i,j,3)/cx(i,j,2)) + SQRT(cx(i,j,4)/cx(i,j,3)) + (cx(i,j,4)/cx(i,j,2))**(0.25_DP)) / 3
             rvdw(i,j) = a1 * rc + a2
             rvdw(j,i) = rvdw(i,j)

             WRITE (stdout,'(I3,1X,I3,1p,5(1X,E14.6))') i, j, cx(i,j,2), cx(i,j,3), cx(i,j,4), rc, rvdw(i,j)
          END DO
       END DO

       ! clean up and mark as done
       IF (ALLOCATED(avol)) DEALLOCATE(avol)
       saved = .TRUE.
    END IF

    ! calculate energy contributions
    IF (ionode) THEN
       WRITE (stdout,*)
       WRITE (stdout,'("+ van der Waals energies, forces and stresses (Ry,bohr)")')
    END IF
    evdw = 0._DP
    for = 0._DP
    sigma = 0._DP

    ehadd = 0._DP

    ! set the atomic environment for the energy sum -> it would be nice to rewrite
    ! this using ewald: Williams, Acta Cryst. A 27 (1971) 452 and some other paper, maybe
    ! in the international tables for crystallography.
    rmax = (maxc6/ecut)**(1._DP/6._DP)
    rmax2 = rmax*rmax
    CALL set_environ(rmax)

    ! parallelize over atoms
#if defined __MPI
    resto = MOD ( nat , nproc_image )
    divid = nat / nproc_image
    IF ( me_image + 1 <= resto ) THEN
       first = ( divid  + 1 ) * me_image + 1
       last  = ( divid  + 1 ) * ( me_image + 1 )
    ELSE
       first = ( ( divid + 1 ) * resto ) + ( divid ) * ( me_image-resto ) + 1
       last  = ( divid  + 1 ) * resto + ( divid ) * ( me_image - resto + 1 )
    END IF
#else
    first = 1
    last  = nat
#endif

    DO i = first, last
       sat = 0._DP
       eat = 0._DP
       taub = tau(:,i) * alat

       ! C6, C8, C10
       DO i3 = 2, 4
          ! order R^nn
          nn = 2 * (i3 + 1)

          DO j = 1, nenv
             jj = ienv(j)
             xij = xenv(:,j) - taub
             dij2 = xij(1)*xij(1) + xij(2)*xij(2) + xij(3)*xij(3)
             IF (dij2 < 1e-15_DP .OR. dij2>rmax2) CYCLE
             dijx = dij2**(i3+1)
             dijxm2 = dijx / dij2

             ! energy contribution
             cn0 = cx(i,jj,i3)
             den = 1 / (rvdw(i,jj)**nn + dijx)
             ee = cn0 * den
             ehadd(nn) = ehadd(nn) + ee
             eat = eat + ee

             ! force and stress contribution
             den2 = den * den
             fxx = nn * cn0 * dijxm2 * den2

             for(:,i) = for(:,i) + fxx * xij
             DO m1 = 1, 3
                sat(m1,m1) = sat(m1,m1) + fxx * xij(m1) * xij(m1)
                DO m2 = 1, m1-1
                   sat(m1,m2) = sat(m1,m2) + fxx * xij(m1) * xij(m2)
                END DO
             END DO
          END DO ! j

       END DO ! i3
       sat(1,2) = sat(2,1)
       sat(1,3) = sat(3,1)
       sat(2,3) = sat(3,2)
       evdw = evdw + eat
       sigma = sigma + sat
    END DO ! i
    evdw= -0.5_DP * evdw
    sigma = -0.5_DP * sigma / omega
    ehadd = -0.5_DP * ehadd

#if defined(__MPI)
    CALL mp_sum(evdw,intra_image_comm)
    CALL mp_sum(for,intra_image_comm)
    CALL mp_sum(sigma,intra_image_comm)
    DO nn = 6, 10
       CALL mp_sum(ehadd(nn),intra_image_comm)
    ENDDO
#endif

    ! Convert to Ry
    evdw = evdw * 2
    for = for * 2
    sigma = sigma * 2
    ehadd = ehadd * 2

    ! save energy, forces and stress tensor
    esaveold = esave
    esave = evdw
    fsave = for(:,1:nat)
    ssave = sigma

    IF (ionode) THEN
       WRITE (stdout,'("  Evdw(total,Ry)   = ",1p,E20.12)') evdw
       WRITE (stdout,'("  Evdw(C6,Ry)      = ",1p,E20.12)') ehadd(6)
       WRITE (stdout,'("  Evdw(C8,Ry)      = ",1p,E20.12)') ehadd(8)
       WRITE (stdout,'("  Evdw(C10,Ry)     = ",1p,E20.12)') ehadd(10)
       DO i = 1, nat
          WRITE (stdout,'("  Fvdw (",I3.3,",Ry/bohr) = ",1p,3(E20.12,1X))') i, for(:,i)
       END DO
       WRITE (stdout,'("  sigma_vdw (Ry/bohr**3) = ",1p,3(E20.12,1X)," ")') sigma(1,:)
       WRITE (stdout,'("                           ",1p,3(E20.12,1X)," ")') sigma(2,:)
       WRITE (stdout,'("                           ",1p,3(E20.12,1X)," ")') sigma(3,:)
       WRITE (stdout,'("  sigma_vdw (GPa) = ",1p,3(E20.12,1X)," ")') 0.5_DP*sigma(1,:)*au_gpa
       WRITE (stdout,'("                    ",1p,3(E20.12,1X)," ")') 0.5_DP*sigma(2,:)*au_gpa
       WRITE (stdout,'("                    ",1p,3(E20.12,1X)," ")') 0.5_DP*sigma(3,:)*au_gpa
       WRITE (stdout,*)
    END IF

  END FUNCTION energy_xdm

  FUNCTION force_xdm(nat) RESULT(fvdw)
    ! Fetch the dispersion contribution to forces from a previous energy_xdm execution.
    INTEGER, INTENT(IN) :: nat
    REAL(DP) :: fvdw(3,nat)

    fvdw = fsave

  END FUNCTION force_xdm

  FUNCTION stress_xdm() RESULT(svdw)
    ! Fetch the dispersion contribution to stress from a previous energy_xdm execution.
    REAL(DP) :: svdw(3,3)

    svdw = ssave

  END FUNCTION stress_xdm

  ! --- private ---
  SUBROUTINE set_environ (rcut)
    ! Calculate an atomic environemnt of the entire unit cell up to a distance rcut. 
    ! This environment is saved in the host module arrays ienv, xenv and lvec.
    USE cell_base, ONLY: at, bg, alat, omega, tpiba2
    USE ions_base, ONLY: nat, tau, ityp, atm
    USE io_global, ONLY: stdout, ionode

    REAL(DP), INTENT(IN) :: rcut

    INTEGER :: nadd, imax, jmax, kmax, ialloc
    REAL(DP) :: rmat(3,3), gtensor(3,3), alp, bet, gam, aa, bb, cc, xx(3)
    INTEGER :: ii, jj, kk, m, nsize, lsize
    INTEGER, ALLOCATABLE :: ienvaux(:), lvecaux(:,:)
    REAL(DP), ALLOCATABLE :: xenvaux(:,:)
    INTEGER, PARAMETER :: menv = 1000, lenv=100

    ! allocate the initial environment
    nenv = 0
    IF (ALLOCATED(ienv)) DEALLOCATE(ienv)
    IF (ALLOCATED(xenv)) DEALLOCATE(xenv)
    ALLOCATE(ienv(menv),STAT=ialloc)
    IF (ialloc /= 0) CALL alloc_failed("ienv")
    ALLOCATE(xenv(3,menv),STAT=ialloc)
    IF (ialloc /= 0) CALL alloc_failed("xenv")
    nsize = menv

    ! allocate the array of lattice vectors
    nvec = 0
    IF (ALLOCATED(lvec)) DEALLOCATE(lvec)
    ALLOCATE(lvec(3,lenv),STAT=ialloc)
    IF (ialloc /= 0) CALL alloc_failed("lenv")
    lsize = lenv

    ! determine number of cells (adapted from gulp, by J. Gale)
    rmat = at * alat
    gtensor = MATMUL(TRANSPOSE(rmat),rmat)
    aa = SQRT(gtensor(1,1))
    bb = SQRT(gtensor(2,2))
    cc = SQRT(gtensor(3,3))
    alp = ACOS(gtensor(2,3) / bb / cc) * 180._dp / pi
    bet = ACOS(gtensor(1,3) / aa / cc) * 180._dp / pi
    gam = ACOS(gtensor(1,2) / aa / bb) * 180._dp / pi
    IF (alp<30 .OR. bet<30 .OR. gam<30 .OR. alp>150 .OR. bet>150 .OR. gam>150) THEN
       nadd = 5
    ELSE IF (alp<50 .OR. bet<50 .OR. gam<50 .OR. alp>130 .OR. bet>130 .OR. gam>130) THEN
       nadd = 4
    ELSE IF (alp<70 .OR. bet<70 .OR. gam<70 .OR. alp>110 .OR. bet>110 .OR. gam>110) THEN
       nadd = 3
    ELSE
       nadd = 2
    END IF
    imax = NINT(rcut / aa) + nadd
    jmax = NINT(rcut / bb) + nadd
    kmax = NINT(rcut / cc) + nadd

    ! build the environment arrays
    DO ii = -imax, imax
       DO jj = -jmax, jmax
          DO kk = -kmax, kmax

             ! run over the ions in the (i,j,k) cell:
             DO m = 1, nat
                xx = tau(:,m) + ii*at(:,1) + jj*at(:,2) + kk*at(:,3)

                ! dynamically increase the array size
                nenv = nenv + 1
                IF (nenv > nsize) THEN
                   ALLOCATE(ienvaux(NINT(1.5*nsize)),STAT=ialloc)
                   IF (ialloc /= 0) CALL alloc_failed("ienvaux")
                   ALLOCATE(xenvaux(3,NINT(1.5*nsize)),STAT=ialloc)
                   IF (ialloc /= 0) CALL alloc_failed("xenvaux")
                   ienvaux(1:nsize) = ienv
                   xenvaux(:,1:nsize) = xenv
                   CALL move_alloc(ienvaux,ienv)
                   CALL move_alloc(xenvaux,xenv)
                   nsize = NINT(1.5*nsize)
                END IF
                xenv(:,nenv) = xx * alat
                ienv(nenv) = m
             ENDDO  ! m

             ! one more lattice vector
             nvec = nvec + 1
             IF (nvec > lsize) THEN
                ALLOCATE(lvecaux(3,NINT(1.5*lsize)),STAT=ialloc)
                IF (ialloc /= 0) CALL alloc_failed("lvecaux")
                lvecaux(:,1:lsize) = lvec
                CALL move_alloc(lvecaux,lvec)
                lsize = NINT(1.5*lsize)
             END IF
             lvec(:,nvec) = (/ii,jj,kk/)
          END DO ! kk
       END DO ! jj
    END DO ! ii

    ! fit memory snugly
    ALLOCATE(ienvaux(nsize),STAT=ialloc)
    IF (ialloc /= 0) CALL alloc_failed("ienvaux")
    ienvaux(1:nsize) = ienv
    CALL move_alloc(ienvaux,ienv)
    ALLOCATE(xenvaux(3,nsize),STAT=ialloc)
    IF (ialloc /= 0) CALL alloc_failed("xenvaux")
    xenvaux(:,1:nsize) = xenv
    CALL move_alloc(xenvaux,xenv)
    ALLOCATE(lvecaux(3,lsize),STAT=ialloc)
    IF (ialloc /= 0) CALL alloc_failed("lvecaux")
    lvecaux(:,1:lsize) = lvec
    CALL move_alloc(lvecaux,lvec)

  END SUBROUTINE set_environ

  SUBROUTINE PAW_make_ae_charge_xdm(rho,rhoout)
    ! Build the true valence electron density from the pseudo-electron density using
    ! the PAW transformation. This is necessary for the calculation of the atom-in-molecule
    ! volumes. Adapted from PP.
    USE paw_variables, ONLY : paw_info
    USE paw_onecenter, ONLY : paw_rho_lm
    USE atom,          ONLY : g => rgrid
    USE ions_base,     ONLY : nat, ityp, tau, ntyp => nsp
    USE lsda_mod,      ONLY : nspin
    USE uspp_param,    ONLY : nh, nhm, upf
    USE scf,           ONLY : scf_type
    USE fft_base,      ONLY : dfftp
    USE mp,            ONLY : mp_bcast, mp_sum
    USE mp_pools,      ONLY : me_pool
    USE mp_images,     ONLY : intra_image_comm
    USE io_global,     ONLY : ionode_id
    USE splinelib,     ONLY : spline, splint
    USE cell_base,     ONLY : at, bg, alat

    TYPE(scf_type), INTENT(IN) :: rho
    REAL(DP), INTENT(OUT) :: rhoout(dfftp%nnr)

    TYPE(paw_info)          :: i
    INTEGER                 :: ipol, ir, is, it, lm
    INTEGER                 :: j, k, l, idx, idx0
    INTEGER                 :: ia, il, im, ml, mm
    REAL(DP),ALLOCATABLE    :: wsp_lm(:,:), ylm_posi(:,:), d1y(:,:), d2y(:,:)
    REAL(DP),ALLOCATABLE    :: rho_lm(:,:,:), rho_lm_ae(:,:,:), rho_lm_ps(:,:,:)
    REAL(DP)                :: posi(3), first, second
    REAL(DP)                :: inv_nr1, inv_nr2, inv_nr3, distsq, g0, g1, g2, r0, r1, rqq
    INTEGER                 :: nkk
    INTEGER, ALLOCATABLE    :: iatom(:)

    ! Some initialization
    inv_nr1 = 1._DP / DBLE(  dfftp%nr1 )
    inv_nr2 = 1._DP / DBLE(  dfftp%nr2 )
    inv_nr3 = 1._DP / DBLE(  dfftp%nr3 )

    ! copy the density to the output
    rhoout = 0._DP
    DO ir = 1, dfftp%nnr
       DO is = 1, nspin
          rhoout(ir) = rhoout(ir) + rho%of_r(ir,is)
       END DO
    END DO

    ! allocate the rho_lm array
    mm = 0
    ml = 0
    DO it = 1, ntyp
       mm = MAX(g(it)%mesh,mm)
       ml = MAX(upf(it)%lmax_rho + 1,ml)
    END DO
    ALLOCATE(rho_lm(mm,ml**2,nat))
    rho_lm = 0._DP

    ! count the number of processors per atom
    ALLOCATE(iatom(nat))
    iatom = 0
    DO ia = 1, nat
       IF (ASSOCIATED(upf(ityp(ia))%paw%pfunc)) iatom(ia) = iatom(ia) + 1
    END DO
    CALL mp_sum(iatom,intra_image_comm)

    ! run over atoms and build rho_lm. Not all atom types are visible to all processors.
    DO ia = 1, nat
       i%a = ia                      ! atom's index
       i%t = ityp(ia)                ! type of atom ia
       i%m = g(i%t)%mesh             ! radial mesh size for atom i%t
       i%b = upf(i%t)%nbeta          ! number of beta functions for i%t
       i%l = upf(i%t)%lmax_rho+1     ! max ang.mom. in augmentation for ia

       IF (.NOT.upf(i%t)%tpawp) call errore('paw_make_ae_charge_xdm','non-paw pseudo',1)
       IF (.NOT.ASSOCIATED(upf(i%t)%paw%pfunc)) CYCLE

       ALLOCATE(rho_lm_ae(i%m,i%l**2,nspin), rho_lm_ps(i%m,i%l**2,nspin))

       CALL PAW_rho_lm(i, rho%bec, upf(i%t)%paw%pfunc,  rho_lm_ae)
       CALL PAW_rho_lm(i, rho%bec, upf(i%t)%paw%ptfunc, rho_lm_ps, upf(i%t)%qfuncl)

       DO is=1,nspin
          DO lm = 1,i%l**2
             DO ir = 1, i%m
                rho_lm(ir,lm,ia) = rho_lm(ir,lm,ia) + (rho_lm_ae(ir,lm,is) - rho_lm_ps(ir,lm,is) ) * g(i%t)%rm2(ir) / iatom(ia)
             ENDDO
          ENDDO
       ENDDO
       DEALLOCATE(rho_lm_ae, rho_lm_ps)
    END DO
    call mp_sum(rho_lm,intra_image_comm)
    DEALLOCATE(iatom)

    ! Not parallelizing over atoms, because it is already parallelized over charge slabs
    atoms: DO ia = 1, nat
       i%a = ia                      ! atom's index
       i%t = ityp(ia)                ! type of atom ia
       i%m = g(i%t)%mesh             ! radial mesh size for atom i%t
       i%b = upf(i%t)%nbeta          ! number of beta functions for i%t
       i%l = upf(i%t)%lmax_rho+1 ! max ang.mom. in augmentation for ia

       ! spline the rho_lm
       ALLOCATE( d1y(upf(i%t)%kkbeta,i%l**2))
       ALLOCATE( d2y(upf(i%t)%kkbeta,i%l**2))
       ALLOCATE(wsp_lm(i%m,i%l**2))
       DO lm = 1, i%l**2
          CALL radial_gradient(rho_lm(1:upf(i%t)%kkbeta,lm,ia),d1y(:,lm),g(i%t)%r,upf(i%t)%kkbeta,1)
          CALL radial_gradient(d1y(:,lm),d2y(:,lm),g(i%t)%r,upf(i%t)%kkbeta,1)
          first  = d1y(1,lm) ! first derivative in first point
          second = d2y(1,lm) ! second derivative in first point
          ! prepare interpolation
          CALL spline( g(i%t)%r(:), rho_lm(1:i%m,lm,ia), first, second, wsp_lm(:,lm) )
       ENDDO
       DEALLOCATE(d1y,d2y)

       ! define the starting index for each processor
       idx0 = dfftp%nr1x * dfftp%nr2x * dfftp%ipp(me_pool+1)

       ALLOCATE(ylm_posi(1,i%l**2))
       rsp_point : DO ir = 1, dfftp%nr1x*dfftp%nr2x * dfftp%npl
          ! three dimensional indices (i,j,k)
          idx   = idx0 + ir - 1
          k     = idx / ( dfftp%nr1x* dfftp%nr2x)
          idx   = idx - ( dfftp%nr1x* dfftp%nr2x)*k
          j     = idx /  dfftp%nr1x
          idx   = idx -  dfftp%nr1x*j
          l     = idx

          ! ... do not include points outside the physical range!
          IF ( l >=  dfftp%nr1 .OR. j >=  dfftp%nr2 .OR. k >=  dfftp%nr3 ) CYCLE rsp_point
          !
          DO ipol = 1, 3
             posi(ipol) = DBLE( l )*inv_nr1*at(ipol,1) + &
                DBLE( j )*inv_nr2*at(ipol,2) + &
                DBLE( k )*inv_nr3*at(ipol,3)
          ENDDO
          !
          ! find the distance of real-space grid's point ir w.r.t
          ! closer periodic image of atom ia
          !
          posi(:) = posi(:) - tau(:,ia)
          CALL cryst_to_cart( 1, posi, bg, -1 )
          posi(:) = posi(:) - ANINT( posi(:) )
          CALL cryst_to_cart( 1, posi, at, 1 )
          !
          posi(:) = posi(:) * alat
          distsq = posi(1)**2 + posi(2)**2 + posi(3)**2
          ! don't consider points too far from the atom:
          IF ( distsq > g(i%t)%r2(upf(i%t)%kkbeta) ) &
             CYCLE rsp_point
          !
          ! generate the atomic charge on point posi(:), which means
          ! sum over l and m components rho_lm_ae-rho_lm_ps
          ! interpolate the radial function at distance |posi(:)|
          !
          ! prepare spherical harmonics
          CALL ylmr2( i%l**2, 1, posi, distsq, ylm_posi )

          rqq = SQRT(distsq)
          DO lm = 1, i%l**2
             rhoout(ir)= rhoout(ir) + ylm_posi(1,lm)  * splint(g(i%t)%r(:),rho_lm(:,lm,ia),wsp_lm(:,lm),rqq)
          ENDDO
       ENDDO rsp_point
       DEALLOCATE(ylm_posi, wsp_lm)
    ENDDO atoms
    DEALLOCATE(rho_lm)

  END SUBROUTINE PAW_make_ae_charge_xdm

  SUBROUTINE promolecular_rho(rhot,rhoc)
    ! Calculate the promolecular density (i.e. the sum of atomic
    ! densitites) and the sum of core densities on the real-space
    ! grid. Unfortunately, aliasing errors prevent using the atomic
    ! form factor trick, so we're stuck with summing over an
    ! environment. 
    USE io_global, ONLY: ionode
    USE kinds,     ONLY : DP
    USE atom,      ONLY : rgrid, msh
    USE ions_base, ONLY : ityp, ntyp => nsp
    USE cell_base, ONLY : at
    USE fft_base,  ONLY : dfftp
    USE splinelib, ONLY : splint
    use cell_base, ONLY : alat
    USE mp_pools,  ONLY : me_pool
    implicit none

    real(DP), intent(out) :: rhoc(dfftp%nnr) ! core density in the real-space grid
    real(DP), intent(out) :: rhot(dfftp%nnr) ! core density in the real-space grid

    integer :: i, it, nn
    integer :: n, idx0, idx, ix, iy, iz
    real(DP) :: x(3), xx(3), r, r2, rrho

    ! define the starting index for each processor
    idx0 = dfftp%nr1x * dfftp%nr2x * dfftp%ipp(me_pool+1)

    rhot = 0._DP
    rhoc = 0._DP

    ! run over the real-space density grid
    DO n = 1, dfftp%nnr
       idx = idx0 + n - 1

       iz = idx / (dfftp%nr1x*dfftp%nr2x)
       idx = idx - (dfftp%nr1x*dfftp%nr2x)*iz
       iy = idx / dfftp%nr1x
       idx = idx - dfftp%nr1x*iy
       ix = idx

       x = ix / REAL(dfftp%nr1,DP) * at(:,1) + iy / REAL(dfftp%nr2,DP) * at(:,2) + &
          iz / REAL(dfftp%nr3,DP) * at(:,3)
       x = x * alat

       ! contributions from the environment
       DO i = 1, nenv
          it = ityp(ienv(i))
          nn = msh(it)
          xx = x - xenv(:,i)
          r2 = xx(1)*xx(1) + xx(2)*xx(2) + xx(3)*xx(3)
          IF (r2 > rmaxg2(it)) CYCLE
          r = SQRT(r2)
          rrho = splint(rgrid(it)%r(1:nn),rfree(1:nn,it),w2free(1:nn,it),r)
          rhot(n) = rhot(n) + rrho

          IF (r2 > rmaxcore2(it)) CYCLE
          rrho = splint(rgrid(it)%r(1:nn),rcore(1:nn,it),w2core(1:nn,it),r)
          rhoc(n) = rhoc(n) + rrho
       END DO
       rhot(n) = MAX(rhot(n),1e-14_DP)
    END DO

  END SUBROUTINE promolecular_rho

  SUBROUTINE alloc_failed(message)
    ! Error message and horrible death
    CHARACTER*(*), INTENT(IN) :: message

    CALL errore('energy_xdm','allocation failed: '//TRIM(ADJUSTL(message)),1)

  END SUBROUTINE alloc_failed

END MODULE xdm_module
