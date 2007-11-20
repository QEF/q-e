!
! Copyright (C) 2007 Quantum-Espresso group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "f_defs.h"
#define OPTIONAL_CALL if(.false.) CALL
!
MODULE paw_init
  !
  USE kinds, ONLY : DP
  !
  IMPLICIT NONE

  PUBLIC :: PAW_init_becsum
  PUBLIC :: PAW_init_onecenter
#ifdef __PARA
  PUBLIC :: PAW_post_init
#endif

  PUBLIC :: allocate_paw_internals, deallocate_paw_internals

!!!=========================================================================
 CONTAINS

  ! Allocate PAW internal variables require for SCF calculation
  SUBROUTINE allocate_paw_internals
    USE lsda_mod,           ONLY : nspin
    USE ions_base,          ONLY : nat, ntyp => nsp
    USE uspp_param,         ONLY : nhm
    USE gvect,              ONLY : ngl
    !
    USE paw_variables
    !
    IMPLICIT NONE
    !
    ALLOCATE(ddd_paw( nhm, nhm, nat, nspin))
    !
  END SUBROUTINE allocate_paw_internals

  ! Called from clean_pw
  SUBROUTINE deallocate_paw_internals
    USE paw_variables
    USE uspp_param, ONLY : upf
    USE ions_base,  ONLY : ntyp => nsp
    !
    IMPLICIT NONE
    INTEGER :: nt
    !
    IF(allocated(ddd_paw)) DEALLOCATE (ddd_paw)
    !
    ! Allocated in read_paw:
    DO nt = 1,ntyp
        IF(associated(upf(nt)%paw%aug)) DEALLOCATE (upf(nt)%paw%aug)
    ENDDO
    !
  END SUBROUTINE deallocate_paw_internals

#ifdef __PARA
! Deallocate variables that are used only at init and then no more necessary.
! This is only useful in parallel, as each node only does a limited number of atoms
SUBROUTINE PAW_post_init()
    ! this routine does nothing at this moment...
    USE ions_base,          ONLY : nat, ntyp=>nsp, ityp
    USE uspp_param,         ONLY : upf
    USE mp_global,          ONLY : mpime, nproc
    USE io_global,          ONLY : stdout, ionode
    !
    INTEGER :: nt, na, np, first_nat, last_nat
    LOGICAL :: info(nproc,ntyp)

    IF(ionode) &
    WRITE(stdout,"(5x,a)") &
        'Checking if some PAW data can be deallocated... '
    info(:,:) = .false.
    CALL divide (nat, first_nat, last_nat)
    !
    types : &
    DO nt = 1,ntyp
        DO na = first_nat, last_nat
            IF (ityp(na) == nt ) CYCLE types
        ENDDO
        ! If I can't find any atom within first_nat and last_nat
        ! which is of type nt, then I can deallocate:
        DEALLOCATE( upf(nt)%paw%aug,        &
                    upf(nt)%paw%ae_rho_atc, &
                    upf(nt)%paw%pfunc,      &
                    upf(nt)%paw%ptfunc,     &
                    upf(nt)%paw%ae_vloc     &
                  )
        DEALLOCATE( upf(nt)%vloc)
!        DEALLOCATE( upf(nt)%rho_atc)
        info(mpime,nt) = .true.
    ENDDO types
    CALL reduce(nproc*ntyp, info)

    IF(ionode) THEN
        DO np = 1,nproc
        DO nt = 1,ntyp
            IF( info(np,nt) ) &
            WRITE(*,"(7x,a,i3,a,10i3)") "node ",np,&
                    ", deallocated PAW data for type:", nt
        ENDDO
        ENDDO
    ENDIF

END SUBROUTINE PAW_post_init
#endif

! Initialize becsum with atomic occupations (for PAW atoms only)
! Notice: requires exact correspondence chi <--> beta in the atom,
! that is that all wavefunctions considered for PAW generation are
! counted in chi (otherwise the array "oc" does not correspond to beta)
SUBROUTINE PAW_init_becsum()
    USE uspp,               ONLY : becsum, nhtol, indv
    USE uspp_param,         ONLY : upf, nh, upf
    USE ions_base,          ONLY : nat, ityp
    USE lsda_mod,           ONLY : nspin, starting_magnetization
    USE paw_variables,      ONLY : okpaw
    USE paw_onecenter,      ONLY : PAW_symmetrize
    IMPLICIT NONE
    INTEGER :: ispin, na, nt, ijh, ih, jh, nb, mb
    !
    IF (.NOT. okpaw) RETURN
    !
    if (nspin.GT.2) &
        CALL errore('PAW_init_becsum', 'Atomic becsum not implemented for nspin>2',1)
    !
    na_loop: DO na = 1, nat
       nt = ityp(na)
       is_paw: IF (upf(nt)%tpawp) THEN
          !
          ijh = 1
          ih_loop: DO ih = 1, nh(nt)
             nb = indv(ih,nt)
             !
             IF (nspin==1) THEN
                !
                becsum(ijh,na,1) = upf(nt)%paw%oc(nb) / REAL(2*nhtol(ih,nt)+1,DP)
                !
             ELSE IF (nspin==2) THEN
                !
                becsum(ijh,na,1) = 0.5d0 * (1.d0 + starting_magnetization(nt))* &
                                   upf(nt)%paw%oc(nb) / REAL(2*nhtol(ih,nt)+1,DP)
                becsum(ijh,na,2) = 0.5d0 * (1.d0 - starting_magnetization(nt))* &
                                   upf(nt)%paw%oc(nb) / REAL(2*nhtol(ih,nt)+1,DP)
                !
             END IF
             ijh = ijh + 1
             !
             jh_loop: &
              DO jh = ( ih + 1 ), nh(nt)
                !mb = indv(jh,nt)
                DO ispin = 1, nspin
                   becsum(ijh,na,ispin) = 0._DP
                END DO
                ijh = ijh + 1
                !
             END DO jh_loop
          END DO ih_loop
       END IF is_paw
    END DO na_loop

    CALL PAW_symmetrize(becsum)

END SUBROUTINE PAW_init_becsum

! This allocates space to store onecenter potential and 
! calls PAW_rad_init to initialize onecenter integration.
SUBROUTINE PAW_init_onecenter()
    USE ions_base,              ONLY : nat, ityp, ntyp => nsp
    USE paw_variables,          ONLY : xlm, saved, rad, paw_is_init
    USE atom,                   ONLY : g => rgrid
    USE radial_grids,           ONLY : do_mesh
    USE uspp_param,             ONLY : lmaxq, upf
    USE lsda_mod,               ONLY : nspin
    USE funct,                  ONLY : dft_is_gradient

    INTEGER :: na, nt, first_nat, last_nat, lmax_safe

    IF( paw_is_init ) THEN
        CALL infomsg('PAW_init_onecenter', 'Already initialized!', 1)
        RETURN
    ENDIF

    ! First a bit of generic initialization:
    IF (allocated(saved)) DEALLOCATE(saved)
    ALLOCATE(saved(nat)) ! allocate space to store the potentials
    !
    ! Parallelizing this loop every node only allocs the potential
    ! for the atoms that it will actually use later.
    CALL divide (nat, first_nat, last_nat)
    DO na = first_nat, last_nat
        nt = ityp(na)
        ! note that if the atom is not paw it is left unallocated
        IF ( upf(nt)%tpawp ) THEN
            IF (associated(saved(na)%v)) DEALLOCATE(saved(na)%v)
            ALLOCATE( saved(na)%v(g(nt)%mesh, (upf(nt)%paw%lmax_rho+1)**2, nspin, 2 ) )
                      !                                     {AE|PS}
        ENDIF
    ENDDO

    ! initialize for integration on angular momentum and gradient, integrating
    ! up to 2*lmaxq (twice the maximum angular momentum of rho) is enough for
    ! H energy and for XC energy. If I have gradient correction I have to go a bit higher
    IF(allocated(rad)) DEALLOCATE( rad )
    ALLOCATE( rad(ntyp) )
    !
    IF ( dft_is_gradient() ) THEN
        lmax_safe = xlm
    ELSE
        lmax_safe = 0
    ENDIF
    !
    types : &
    DO nt = 1,ntyp
        ! only allocate radial grid integrator for atomic species
        ! that are actually present on this parallel node:
        DO na = first_nat, last_nat
        IF (ityp(na) == nt ) THEN
            CALL PAW_rad_init(2*upf(nt)%paw%lmax_rho + lmax_safe, rad(nt))
            !
            CYCLE types
        ENDIF
        ENDDO
    ENDDO types

    paw_is_init = .true.

END SUBROUTINE PAW_init_onecenter

!___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___
!!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!! 
!!! initialize several quantities related to radial integration: spherical harmonics and their 
!!! gradients along a few (depending on lmaxq) directions, weights for spherical integration
!!
!!! IMPORTANT: routine PW/summary.f90 has the initialization parameters hardcoded inside
!!!            remember to update it if you change this!!!
SUBROUTINE PAW_rad_init(l, rad)
    USE constants,              ONLY : pi, fpi, eps8
    USE funct,                  ONLY : dft_is_gradient
    USE paw_variables,          ONLY : paw_radial_integrator
    INTEGER,INTENT(IN)          :: l ! max angular momentum component that will be
                                     ! integrated exactly (to numerical precision)
    TYPE(paw_radial_integrator),INTENT(OUT) :: &
                                   rad ! containt weights and more info to integrate
                                       ! on radial grid up to lmax = l

    REAL(DP),ALLOCATABLE        :: x(:),&       ! nx versors in smart directions
                                   w(:),&       ! temporary integration weights
                                   r(:,:),&     ! integration directions
                                   r2(:),&      ! square modulus of r
                                   ath(:),aph(:)! angles in sph coords for r

    INTEGER                     :: i,ii,n       ! counters
    INTEGER                     :: lm,lm2,m     ! indexes for ang.mom
    REAL(DP)                    :: phi,dphi,rho ! spherical coordinates
    REAL(DP)                    :: z            ! cartesian coordinates
    ! for gradient corrections:
    INTEGER                     :: ipol
    REAL(DP),ALLOCATABLE        :: aux(:,:),&   ! workspace
                                   s(:,:),&     ! integration directions + delta
                                   s2(:)        ! square modulus of s
    REAL(DP)                    :: vth(3), vph(3) !versors for theta and phi
    !
    CHARACTER(len=100)          :: message

    OPTIONAL_CALL start_clock ('PAW_rad_init')

    ! maximum value of l correctly integrated
    rad%lmax = l
    ! volume element for angle phi
    dphi = 2._dp*pi/(rad%lmax+1)
    ! number of samples for theta angle
    n = (rad%lmax+2)/2
    ALLOCATE(x(n),w(n))
    ! compute weights for theta integration
    CALL weights(x,w,n)

    ! number of integration directions
    rad%nx = n*(rad%lmax+1)
!     WRITE(message,"(a,i3,a,i2)") "Setup to integrate on ",&
!           rad%nx," directions; integration exact up to l = ",rad%lmax
!     CALL infomsg('PAW_rad_init', message)
    ALLOCATE(r(3,rad%nx),r2(rad%nx), rad%ww(rad%nx), ath(rad%nx), aph(rad%nx))

    ! compute real weights multiplying theta and phi weights
    ii = 0
    do i=1,n
        z = x(i)
        rho=sqrt(1._dp-z**2)
        do m=0,rad%lmax
            ii= ii+1
            phi = dphi*m
            r(1,ii) = rho*cos(phi)
            r(2,ii) = rho*sin(phi)
            r(3,ii) = z
            rad%ww(ii) = w(i)*2._dp*pi/(rad%lmax+1)
            r2(ii) = r(1,ii)**2+r(2,ii)**2+r(3,ii)**2
            ! these will be used later:
            ath(ii) = acos(z/sqrt(r2(ii)))
            aph(ii) = phi
        end do
    end do
    ! cleanup
    DEALLOCATE (x,w)

    ! initialize spherical harmonics that will be used
    ! to convert rho_lm to radial grid
    rad%lm_max = (rad%lmax+1)**2
    ALLOCATE( rad%ylm(rad%nx, rad%lm_max) )
    CALL ylmr2(rad%lm_max, rad%nx, r,r2,rad%ylm)
    ! As I will mostly use the product ww*ylm I can 
    ! precompute it here:
    ALLOCATE( rad%wwylm(rad%nx, rad%lm_max) )
    DO i = 1,rad%nx
        DO lm = 1, rad%lm_max
            rad%wwylm(i, lm) = rad%ww(i) * rad%ylm(i, lm)
        ENDDO
    ENDDO

    ! if gradient corrections will be used than we need
    ! to initialize the gradient of ylm, as we are working in spherical
    ! coordinates the formula involves \hat{theta} and \hat{phi}
    gradient: IF (dft_is_gradient()) THEN
        ALLOCATE(s(3,rad%nx),s2(rad%nx))
        ALLOCATE( rad%dylmt(rad%nx,rad%lm_max),&
                  rad%dylmp(rad%nx,rad%lm_max),&
                  aux(rad%nx,rad%lm_max)       )
        ALLOCATE(rad%cotg_th(rad%nx) )
        rad%dylmt(:,:) = 0._dp
        rad%dylmp(:,:) = 0._dp

        ! compute derivative along x, y and z => gradient, then compute the
        ! scalar products with \hat{theta} and \hat{phi} and store them in
        ! dylmt and dylmp respectively
        DO ipol = 1,3 !x,y,z
            CALL dylmr2(rad%lm_max, rad%nx, r,r2, aux, ipol)
            DO lm = 1, rad%lm_max
            DO i = 1, rad%nx
                vph = (/-sin(aph(i)), cos(aph(i)), 0._dp/)
                ! this is the explicit form, but the cross product trick (below) is much faster:
                ! vth = (/cos(aph(i))*cos(ath(i)), sin(aph(i))*cos(ath(i)), -sin(ath(i))/)
                vth = (/vph(2)*r(3,i)-vph(3)*r(2,i),&
                        vph(3)*r(1,i)-vph(1)*r(3,i),&
                        vph(1)*r(2,i)-vph(2)*r(1,i)/)
                rad%dylmt(i,lm) = rad%dylmt(i,lm) + aux(i,lm)*vth(ipol)
                ! CHECK: the 1/sin(th) factor should be correct, but deals wrong result, why?
                rad%dylmp(i,lm) = rad%dylmp(i,lm) + aux(i,lm)*vph(ipol) !/sin(ath(i))
                rad%cotg_th(i) = cos(ath(i))/sin(ath(i))
            ENDDO
            ENDDO
        ENDDO
        DEALLOCATE(aux)
    ENDIF gradient
    ! cleanup
    DEALLOCATE (r,r2)

    OPTIONAL_CALL stop_clock ('PAW_rad_init')

 CONTAINS
    ! Computes weights for gaussian integrals,
    ! from numerical recipes
    SUBROUTINE weights(x,w,n)
    implicit none
    integer :: n, i,j,m
    real(8), parameter :: eps=1.d-14
    real(8) :: x(n),w(n), z,z1, p1,p2,p3,pp,pi
    
    pi = 4._dp*atan(1._dp)
    m=(n+1)/2
    do i=1,m
        z1 = 2._dp
        z=cos(pi*(i-0.25_dp)/(n+0.5_dp))
        do while (abs(z-z1).gt.eps)
        p1=1._dp
        p2=0._dp
        do j=1,n
            p3=p2
            p2=p1
            p1=((2._dp*j-1._dp)*z*p2-(j-1._dp)*p3)/j
        end do
        pp = n*(z*p1-p2)/(z*z-1._dp)
        z1=z
        z=z1-p1/pp
        end do
        x(i) = -z
        x(n+1-i) = z
        w(i) = 2._dp/((1._dp-z*z)*pp*pp)
        w(n+1-i) = w(i)
    end do

    END SUBROUTINE weights
END SUBROUTINE PAW_rad_init 


END MODULE paw_init
