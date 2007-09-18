!
! Copyright (C) 2001-2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
MODULE rad_paw_routines
    !
    USE kinds,      ONLY : DP
    !
    IMPLICIT NONE
    PUBLIC
    SAVE
    !
    LOGICAL              :: is_init = .false.
    ! if set to true (in init) we have to do gradient correction
    LOGICAL              :: do_gcxc = .false.

    ! the following variables are used to convert spherical harmonics expansion
    ! to radial sampling, they are initialized for an angular momentum up to
    ! l = max_l and (l+1)**2 = lm_max = nx
    ! see function PAW_rad_init for details
    INTEGER              :: l_max  = 0
    INTEGER              :: lm_max = 0
    INTEGER              :: nx     = 0
    REAL(DP),ALLOCATABLE :: ww(:)
    ! additional variables for gradient correction
    REAL(DP),PARAMETER   :: delta = 1.e-5_dp ! for numerical derivatives
    REAL(DP),ALLOCATABLE :: ylm(:,:),&  ! Y_lm(nx,lm_max)
                            dylmt(:,:),&! |d(ylm)/dtheta|**2
                            dylmp(:,:)  ! |d(ylm)/dphi|**2

CONTAINS
! these has to be modularized too:
!#include "../atomic-fake/vxc_t.f90"
!#include "../atomic-fake/exc_t.f90"
!#include "../atomic-fake/vxcgc.f90"

!___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___
!!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!
!!! This is the main driver of PAW routines, it uses directly or indirectly
!!! all the other routines of the module
!!
SUBROUTINE PAW_energy(becsum,correction)
    USE kinds,                  ONLY : DP
    USE constants,              ONLY : pi
    USE lsda_mod,               ONLY : nspin
    USE radial_grids,           ONLY : ndmx
    USE ions_base,              ONLY : nat, ityp

    USE grid_paw_variables,     ONLY : pfunc, ptfunc, tpawp, &
                                       aerho_atc, psrho_atc, aug
    USE uspp_param,             ONLY : nhm, lmaxq

    REAL(DP), INTENT(IN)    :: becsum(nhm*(nhm+1)/2,nat,nspin)! cross band occupations
    REAL(DP), INTENT(OUT),OPTIONAL :: &
                               correction(nat,2,2) ! {# of atoms}, {XC|H}, {AE|PS}

    INTEGER, PARAMETER      :: AE = 1, PS = 2,& ! All-Electron and Pseudo
                               XC = 1, H  = 2   ! XC and Hartree
    INTEGER                 :: i_what           ! counter on AE and PS
    INTEGER                 :: na,nt,first_nat,last_nat ! atoms counters and indexes

    ! hartree energy scalar fields expanded on Y_lm
    REAL(DP), ALLOCATABLE   :: rho_lm(:,:,:)      ! radial density expanded on Y_lm
    REAL(DP), ALLOCATABLE   :: v_h_lm(:,:)        ! hartree potential, summed on spins
    !
    REAL(DP)                :: e_h(lmaxq**2,nspin)! hartree energy components
    REAL(DP)                :: e,e1,e2            ! placeholders
    ! xc variables:
    REAL(DP)                :: e_xc(lmaxq**2,nspin)! hartree energy components
    REAL(DP), POINTER       :: rho_core(:,:)      ! pointer to AE/PS core charge density 

    ! BEWARE THAT HARTREE ONLY DEPENDS ON THE TOTAL RHO NOT ON RHOUP AND RHODW SEPARATELY...
    ! TREATMENT OF NSPIN>1 MUST BE CHECKED AND CORRECTED

    CALL start_clock ('PAW_energy')

    ! initialize for integration on angular momentum and gradient
    !CALL PAW_rad_init(4*lmaxq)
    CALL PAW_rad_init(6*lmaxq)
    ! nullify energy components:
    IF( present(correction) ) correction(:,:,:) = 0._dp
    ! The following operations will be done, first on AE, then on PS part
    ! furthermore they will be done for one atom at a time (to reduce memory usage)
    ! in the future code will be parallelized on atoms:
    !
    ! 1. build rho_lm (PAW_rho_lm)
    ! 2. compute v_h_lm and use it to compute HARTREE 
    !    energy (PAW_h_energy)
    ! 3. compute XC energy
    !   a. build rho_rad(theta, phi) from rho_lm
    !     + if using gradient correction compute also grad(rho)
    !   b. compute XC energy from rho_rad
    !   c. iterate on enough directions to compute the 
    !      spherical integral

    ! CHECK: maybe we don't need to alloc/dealloc rho_lm every time

    first_nat = 1
    last_nat  = nat
    ! Operations from here on are (will be) parallelized on atoms
    atoms: DO na = first_nat, last_nat
        !
        nt = ityp(na) ! the type of atom na
        ifpaw: IF (tpawp(nt)) THEN
        whattodo: DO i_what = AE, PS
            ! STEP: 1 [ build rho_lm (PAW_rho_lm) ]
            ALLOCATE(rho_lm(ndmx,lmaxq**2,nspin))
            NULLIFY(rho_core)
            IF (i_what == AE) THEN
                ! passing "na" as an argument is dirtyer but faster and
                ! uses less memory than passing only a hyperslice of the array
                CALL PAW_rho_lm(na, becsum, pfunc, rho_lm)
                ! used later for xc energy:
                rho_core => aerho_atc
            ELSE
                CALL PAW_rho_lm(na, becsum, ptfunc, rho_lm, aug) !fun)
                !     optional argument for pseudo part --> ^^^^^^
                ! used later for xc energy:
                rho_core => psrho_atc
            ENDIF
            ! STEP: 2 [ compute Hartree energy ]
            ALLOCATE(v_h_lm(ndmx,lmaxq**2))
            !   2a. use rho_lm to compute hartree potential (PAW_v_h)
            e = PAW_h_energy(na, rho_lm, v_h_lm, e_h)
            IF( present(correction) ) correction(na,H,i_what) = e
!             WRITE(6,*) "========================="
!             WRITE(6,*) "== radial PAW energies =="
!             WRITE(6,*) "== Hartree: ", e
            !WRITE(6,'(a,i1,a,f15.7)') ("==RADIAL PAW ENERGY (LM=",lm,"):",e_h(lm,1),lm=1,lmaxq**2)

            ! STEP: 3 [ compute XC energy ]
            e1 = PAW_xc_energy(na, rho_lm, rho_core, v_h_lm, e_xc,0)
            IF( present(correction) ) correction(na,XC,i_what) = e1
            e2 = PAW_xc_energy(na, rho_lm, rho_core, v_h_lm, e_xc,1)
!             WRITE(6,"(x,a,2i3,3f25.15)") "== XC: ", i_what, na, e1,e2
!             WRITE(6,*) "========================="

            ! check if integration is working
            !e = PAW_sph_integral(rho_lm, v_h_lm)
            !write(6,*) "==radial hartree integral --> ",e
            !
            DEALLOCATE(v_h_lm)
            !
            DEALLOCATE(rho_lm)

        ENDDO whattodo
        ENDIF ifpaw
    ENDDO atoms

    CALL stop_clock ('PAW_energy')


END SUBROUTINE PAW_energy

!___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___
!!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!! 
!!! initialize several quantities related to radial integration: spherical harmonics and their 
!!! gradients along a few (depending on lmaxq) directions, weights for spherical integration
!!
SUBROUTINE PAW_rad_init(l)
    USE constants,              ONLY : pi, fpi, eps8
    USE funct,                  ONLY : dft_is_gradient
    INTEGER,INTENT(IN)          :: l            ! max angular momentum

    REAL(DP),ALLOCATABLE        :: x(:),&       ! nx versors in smart directions
                                   w(:),&       ! temporary integration weights
                                   r(:,:),&     ! integration directions
                                   r2(:),&      ! square modulus of r
                                   ath(:),aph(:)! angles in sph coords for r

    INTEGER                     :: i,ii,n       ! counters
    INTEGER                     :: lm,m         ! indexes for ang.mom
    REAL(DP)                    :: phi,dphi,rho ! spherical coordinates
    REAL(DP)                    :: z            ! cartesian coordinates
    ! for gradient corrections:
    INTEGER                     :: ipol
    REAL(DP),ALLOCATABLE        :: aux(:,:),&   ! workspace
                                   s(:,:),&     ! integration directions + delta
                                   s2(:)        ! square modulus of s
    REAL(DP)                    :: vth(3), vph(3) !versors for theta and phi

    ! reinit if necessary
    IF( is_init ) THEN
        IF ( l /= l_max ) THEN
            CALL infomsg('PAW_rad_init',&
             'PAW radial integration already initialized but for a different l: reinitializing.')
            DEALLOCATE(ww, ylm)
            IF (allocated(dylmt)) DEALLOCATE(dylmt)
            IF (allocated(dylmp)) DEALLOCATE(dylmp)
        ELSE
            ! if already initialized correctly nothing to be done
            RETURN
        ENDIF
    ENDIF

    CALL start_clock ('PAW_rad_init')

    ! maximum value of l correctly integrated
    l_max = l
    ! volume element for angle phi
    dphi = 2.d0*pi/(l_max+1)
    ! number of samples for theta angle
    n = (l_max+2)/2
    ALLOCATE (x(n),w(n))
    ! compute weights for theta integration
    CALL weights(x,w,n)

    ! number of integration directions
    nx = n*(l_max+1)
    ALLOCATE (r(3,nx),r2(nx), ww(nx), ath(nx), aph(nx))

    ! compute real weights multiplying theta and phi weights
    ii = 0
    do i=1,n
        z = x(i)
        rho=sqrt(1.d0-z**2)
        do m=0,l_max
            ii= ii+1
            phi = dphi*m
            r(1,ii) = rho*cos(phi)
            r(2,ii) = rho*sin(phi)
            r(3,ii) = z
            ww(ii) = w(i)*2._dp*pi/(l_max+1)
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
    lm_max = (l_max+1)**2
    ALLOCATE(ylm(nx,lm_max))
    CALL ylmr2(lm_max, nx, r,r2,ylm)

    ! if gradient corrections will be used than we need
    ! to initialize the gradient of ylm, as we are working in spherical
    ! coordinates the formula involves \hat{theta} and \hat{phi}
    do_gcxc = dft_is_gradient()
    gradient: IF (do_gcxc) THEN
        ALLOCATE (s(3,nx),s2(nx))
        ALLOCATE(dylmt(nx,lm_max),dylmp(nx,lm_max),aux(nx,lm_max))
        dylmt(:,:) = 0._dp
        dylmp(:,:) = 0._dp

        ! compute derivative along x, y and z => gradient, then compute the
        ! scalar products with \hat{theta} and \hat{phi} and store them in
        ! dylmt and dylmp respectively
        DO ipol = 1,3 !x,y,z
            CALL dylmr2(lm_max, nx, r,r2, aux, ipol)
            DO lm = 1, lm_max
            DO i = 1,nx
                vph = (/-sin(aph(i)), cos(aph(i)), 0._dp/)
                ! this is the explicit form, but the cross product trick (below) is much faster:
                ! vth = (/cos(aph(i))*cos(ath(i)), sin(aph(i))*cos(ath(i)), -sin(ath(i))/)
                vth = (/vph(2)*r(3,i)-vph(3)*r(2,i),&
                        vph(3)*r(1,i)-vph(1)*r(3,i),&
                        vph(1)*r(2,i)-vph(2)*r(1,i)/)
                dylmt(i,lm) = dylmt(i,lm) + aux(i,lm)*vth(ipol)
                ! CHECK: the 1/cos(th) factor should be correct, but deals wrong result, why?
                dylmp(i,lm) = dylmp(i,lm) + aux(i,lm)*vph(ipol) !/cos(ath(i))
            ENDDO
            ENDDO
        ENDDO
    DEALLOCATE(aux)
    ENDIF gradient

    ! cleanup
    DEALLOCATE (r,r2)

    ! success
    is_init = .true.

    CALL stop_clock ('PAW_rad_init')

CONTAINS
    ! Computes weights for gaussian integrals,
    ! from numerical recipes
    SUBROUTINE weights(x,w,n)
    implicit none
    integer :: n, i,j,m
    real(8), parameter :: eps=1.d-14
    real(8) :: x(n),w(n), z,z1, p1,p2,p3,pp,pi
    
    pi = 4.d0*atan(1.d0)
    m=(n+1)/2
    do i=1,m
        z1 = 2.d0
        z=cos(pi*(i-0.25d0)/(n+0.5d0))
        do while (abs(z-z1).gt.eps)
        p1=1.d0
        p2=0.d0
        do j=1,n
            p3=p2
            p2=p1
            p1=((2.d0*j-1.d0)*z*p2-(j-1.d0)*p3)/j
        end do
        pp = n*(z*p1-p2)/(z*z-1.d0)
        z1=z
        z=z1-p1/pp
        end do
        x(i) = -z
        x(n+1-i) = z
        w(i) = 2.d0/((1.d0-z*z)*pp*pp)
        w(n+1-i) = w(i)
    end do

    END SUBROUTINE weights
END SUBROUTINE PAW_rad_init 
!___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___
!!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!! 
!!! use the density produced by sum_rad_rho to compute xc potential and energy, as
!!! xc functional is not diagonal on angular momentum numerical integartion is performed
FUNCTION PAW_xc_energy(na, rho_lm, rho_core, pot_lm, e_lm, task)
    USE kinds,                  ONLY : DP
    USE constants,              ONLY : fpi
    USE parameters,             ONLY : npsx
    USE radial_grids,           ONLY : ndmx
    USE lsda_mod,               ONLY : nspin
    USE uspp_param,             ONLY : lmaxq
    USE ions_base,              ONLY : ityp
    USE atom,                   ONLY : rgrid

    REAL(DP)              :: PAW_xc_energy      ! total xc energy
    !
    INTEGER,  INTENT(IN)  :: na                         ! the number of the atom
    INTEGER,  INTENT(IN)  :: task!remove me
    REAL(DP), INTENT(IN)  :: rho_lm(ndmx,lmaxq**2,nspin)! charge density as lm components
    REAL(DP), INTENT(IN)  :: rho_core(ndmx,npsx)        ! core charge, radial and spherical
    ! TODO:
    REAL(DP), INTENT(OUT) :: pot_lm(ndmx,lmaxq**2)      ! out: potential as lm components
    REAL(DP), OPTIONAL,INTENT(OUT) :: e_lm(lmaxq**2)    ! out: energy components 
    !
    INTEGER               :: lsd          ! switch to control local spin density
    !
    REAL(DP)              :: rho_loc(2) = (/0._dp, 0._dp/) 
                             ! local density (workspace), up and down
    REAL(DP)              :: e            ! workspace
    REAL(DP)              :: e_rad(ndmx)  ! radial energy (to be integrated)
    REAL(DP)              :: rho_rad(ndmx,nspin) ! workspace (radial slice of rho)
    INTEGER               :: nt, &        ! ityp(na)
                             ix,k         ! counters on directions and radial grid
    ! for gradient correction only:
    REAL(DP),ALLOCATABLE  :: grho_rad(:,:)! workspace (radial slice of grad(rho))

    pot_lm = 0._dp
    e_lm   = 0._dp

    CALL start_clock ('PAW_xc_nrg')
    lsd = nspin-1
    nt = ityp(na)

    ! init for gradient correction
    IF (do_gcxc) ALLOCATE(grho_rad(ndmx,nspin))

    PAW_xc_energy = 0._dp
    DO ix = 1, nx
        ! LDA (and LSDA) part (no gradient correction):
        CALL PAW_lm2rad(ix, rho_lm, rho_rad)
        !
        DO k = 1,rgrid(nt)%mesh
            rho_loc(1:nspin) = rho_rad(k,1:nspin)/rgrid(nt)%r2(k)
            !
            e_rad(k) = exc_t(rho_loc, rho_core(k,nt), lsd)&
                     * (SUM(rho_rad(k,1:nspin))+rho_core(k,nt)*rgrid(nt)%r2(k))
        ENDDO
        gradient_correction:& ! add it!
        IF (do_gcxc) THEN
            IF (task == 1) e_rad(:) = 0._dp ! reset the energy, so that only the correction is displayed
            CALL PAW_grad(na, ix, rho_lm, rho_rad, rho_core, grho_rad)
            !                                      v-------------^
            CALL PAW_gcxc(na, rho_rad,rho_core,grho_rad, e_rad)
        ENDIF gradient_correction
        !
        ! integrate radial slice of xc energy:
        CALL simpson (rgrid(nt)%mesh,e_rad,rgrid(nt)%rab,e)
        ! integrate on sph. surface     v----------------^
        PAW_xc_energy = PAW_xc_energy + e * ww(ix)
    ENDDO

    ! cleanup for gc
    IF (do_gcxc) DEALLOCATE(grho_rad)

    CALL stop_clock ('PAW_xc_nrg')

END FUNCTION PAW_xc_energy
!___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___
!!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!! 
!!! add gradient correction, code mostly adapted from ../atomic/vxcgc.f90
SUBROUTINE PAW_gcxc(na, rho,core,grho,e)
    USE kinds,                  ONLY : DP
    USE ions_base,              ONLY : ityp
    USE radial_grids,           ONLY : ndmx
    USE lsda_mod,               ONLY : nspin
    USE atom,                   ONLY : g => rgrid
    USE parameters,             ONLY : npsx
    USE constants,              ONLY : fpi,pi,e2
    USE funct,                  ONLY : gcxc, gcx_spin, gcc_spin
    !
    INTEGER, INTENT(IN)    :: na              ! atom index
    REAL(DP), INTENT(IN)   :: rho(ndmx,nspin) ! radial density,
    REAL(DP), INTENT(IN)   :: grho(ndmx,nspin)! square gradient of rho
    REAL(DP), INTENT(IN)   :: core(ndmx,npsx) ! spherical core density
    REAL(DP), INTENT(INOUT):: e(ndmx)         ! radial local xc energy
    !
    INTEGER            :: i               ! counter for mesh
    INTEGER            :: nt              ! atom type
    ! workspaces:
    REAL(DP)           :: arho, sgn
    REAL(DP)           :: sx,sc           ! x and c energy from gcxc
    REAL(DP)           :: v1x,v2x,v1c,v2c ! potentials from gcxc (unused)
    REAL(DP),PARAMETER :: eps = 1.e-30_dp ! 1.e-12 may be small enough
    ! for nspin>1
    REAL(DP)           :: rh,grh2,zeta
    REAL(DP)           :: v1cup, v1cdw,v1xup, v1xdw, v2xup, v2xdw !(unused)
    
    nt = ityp(na)
    
    IF (nspin.eq.1) THEN
        DO i=1,g(nt)%mesh
            arho = rho(i,1)/g(nt)%r2(i) + core(i,ityp(na))
            sgn  = sign(1.0_dp,arho)
            arho = abs(arho)
            IF (arho.gt.eps.and.abs(grho(i,1)).gt.eps) THEN
                CALL gcxc(arho,grho(i,1),sx,sc,v1x,v2x,v1c,v2c)
                e(i) = e(i) + sgn *e2 *(sx+sc)*g(nt)%r2(i) !&
            ENDIF
        ENDDO
    ELSE
        !   this is the \sigma-GGA case
        DO i=1,g(nt)%mesh
            CALL gcx_spin (rho(i, 1), rho(i, 2), grho(i,1), grho(i,2), &
                           sx, v1xup, v1xdw, v2xup, v2xdw)
            rh = rho(i, 1) + rho(i, 2)
            IF (rh.gt.eps) THEN
                zeta = (rho (i, 1) - rho (i, 2) ) / rh
                grh2 = (sqrt(grho(i,1)) + sqrt(grho(i,2)) ) **2 
                CALL gcc_spin (rh, zeta, grh2, sc, v1cup, v1cdw, v2c)
            ELSE
                sc = 0.0_dp
!                 v1cup = 0.0_dp
!                 v1cdw = 0.0_dp
!                 v2c = 0.0_dp
            ENDIF
            e(i) = e(i) + e2*(sx+sc)
!             vgc(i,1)= v1xup+v1cup
!             vgc(i,2)= v1xdw+v1cdw
!             h(i,1)  =((v2xup+v2c)*grho(i,1)+v2c*grho(i,2))*r2(i)
!             h(i,2)  =((v2xdw+v2c)*grho(i,2)+v2c*grho(i,1))*r2(i)
        enddo
    ENDIF

END SUBROUTINE PAW_gcxc
!___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___
!!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!
!!! build gradient of radial charge distribution from its spherical harmonics expansion
!!! uses pre-computed rho_rad
SUBROUTINE PAW_grad(na, ix, rho_lm, rho_rad, rho_core, grho_rad)
    USE kinds,                  ONLY : DP
    USE constants,              ONLY : fpi
    USE uspp_param,             ONLY : lmaxq
    USE lsda_mod,               ONLY : nspin
    USE radial_grids,           ONLY : ndmx
    USE ions_base,              ONLY : ityp
    USE parameters,             ONLY : npsx
    USE atom,                   ONLY : g => rgrid

    INTEGER, INTENT(IN)  :: ix ! line of the dylm2 matrix to use
                               ! actually it is one of the nx directions
    INTEGER, INTENT(IN)  :: na ! atom index
    REAL(DP), INTENT(IN) :: rho_lm(ndmx,lmaxq**2,nspin)! Y_lm expansion of rho
    REAL(DP), INTENT(IN) :: rho_rad(ndmx,nspin)        ! radial density along direction ix
    REAL(DP), INTENT(IN) :: rho_core(ndmx,npsx)        ! core density
    REAL(DP), INTENT(OUT):: grho_rad(ndmx,nspin)       ! grad of charge density on rad. grid
    !
    REAL(DP)             :: aux(ndmx),aux2(ndmx)       ! workspace
    INTEGER              :: i, is, lm, nt,k   !

    CALL start_clock ('PAW_grad')
    nt = ityp(na)
    ! 1. build real charge density = rho/r**2 + rho_core
    ! 2. compute the partial derivative of rho_rad
    grho_rad(:,:) = 0._dp
    DO is = 1,nspin
        ! build real charge density
        aux(1:g(nt)%mesh) = rho_rad(1:g(nt)%mesh,is)/g(nt)%r2(1:g(nt)%mesh) &
                          + rho_core(1:g(nt)%mesh,nt)/nspin
        ! numerical derivative by ADC ../atomic/vxcgc.f90
        DO k  = 2,g(nt)%mesh-1
            aux2(k) = (  (g(nt)%r(k+1)-g(nt)%r(k))**2 * (aux(k-1)-aux(k)) &
                        -(g(nt)%r(k-1)-g(nt)%r(k))**2 * (aux(k+1)-aux(k)) &
                      )&
                    / (   (g(nt)%r(k+1)-g(nt)%r(k)  ) &
                        * (g(nt)%r(k-1)-g(nt)%r(k)  ) &
                        * (g(nt)%r(k+1)-g(nt)%r(k-1)) &
                      )
        ENDDO
        ! extremes:
        aux2(g(nt)%mesh)=0.0_dp
        aux2(1)= aux2(2) + (aux2(3)-aux2(2)) &
               * (g(nt)%r(1)-g(nt)%r(2))/(g(nt)%r(3)-g(nt)%r(2))
        ! compute the square
        grho_rad(:,is) = aux2(:)**2
    ENDDO

    aux(:)  = 0._dp
    aux2(:) = 0._dp
    DO is = 1,nspin
    ! Spherical (lm=1) component (that would also include core correction) can be omitted
    ! as its contribution to non-radial derivative is zero
    DO lm = 2,lmaxq**2
        ! 5. [ \sum_{lm} rho(r) (dY_{lm}/dphi /cos(theta))  ]**2
        aux(1:g(nt)%mesh) = aux(1:g(nt)%mesh) + dylmp(ix,lm)* (rho_lm(1:g(nt)%mesh,lm,is))
        ! 6. [ \sum_{lm} rho(r) (dY_{lm}/dtheta)  ]**2
        aux2(1:g(nt)%mesh) = aux2(1:g(nt)%mesh) + dylmt(ix,lm)* (rho_lm(1:g(nt)%mesh,lm,is))
    ENDDO
    ! Square and sum up these 2 components, the (1/r**2)**3 factor come from:
    !  a. 1/r**2 from the derivative in spherical coordinates
    !  b. (1/r**2)**2 from rho_lm being multiplied by r**2 
    grho_rad(1:g(nt)%mesh,is) = grho_rad(1:g(nt)%mesh,is)&
                              + (aux(1:g(nt)%mesh)**2 + aux2(1:g(nt)%mesh)**2)&
                                / g(nt)%r2(1:g(nt)%mesh)**3
    ENDDO

    CALL stop_clock ('PAW_grad')

END SUBROUTINE PAW_grad

!___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___
!!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!! 
!!! use the density produced by sum_rad_rho to compute hartree potential 
!!! the potential is then directly integrated to compute hartree energy
FUNCTION PAW_h_energy(na, rho_lm, pot_lm, e_lm)
    USE kinds,                  ONLY : DP
    USE constants,              ONLY : fpi
    USE parameters,             ONLY : npsx
    USE radial_grids,           ONLY : ndmx, hartree
    USE lsda_mod,               ONLY : nspin
    USE uspp_param,             ONLY : nhm, nh, lmaxq
    USE ions_base,              ONLY : ityp
    USE atom,                   ONLY : rgrid

    REAL(DP)                       :: PAW_h_energy      ! total hartree energy
    !
    INTEGER,  INTENT(IN)  :: na                         ! the number of the atom
    REAL(DP), INTENT(IN)  :: rho_lm(ndmx,lmaxq**2,nspin)! charge density as lm components
    REAL(DP), INTENT(OUT) :: pot_lm(ndmx,lmaxq**2)      ! out: potential as lm components
    REAL(DP), OPTIONAL,INTENT(OUT) :: e_lm(lmaxq**2)    ! out: energy components
    !
    REAL(DP)              :: aux(ndmx)   ! workspace
    REAL(DP)              :: par_energy  ! workspace
    REAL(DP)              :: pref        ! workspace

    INTEGER               :: nt,&        ! atom type
                             ispin, &    ! counter on spins
                             lm,l        ! counter on composite angmom lm = l**2 +m
    CALL start_clock ('PAW_h_energy')

    ! get type of atom
    nt = ityp(na)

    ! init total energy and its lm components
    PAW_h_energy = 0._dp
    IF (present(e_lm)) e_lm(:) = 0._dp

    ! this loop computes the hartree potential using the following formula:
    !               l is the first argument in hartree subroutine
    !               r1 = min(r,r'); r2 = MAX(r,r')
    ! V_h(r) = \sum{lm} Y_{lm}(\hat{r})/(2L+1) \int dr' 4\pi r'^2 \rho^{lm}(r') (r1^l/r2^{l+1})
    !     done here --> ^^^^^^^^^^^^^^^^^^^^^           ^^^^^^^^^^^^^^^^^^^^^^ <-- input to the hartree subroutine
    !                 output from the h.s. --> ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    DO lm = 1, lmaxq**2
        l = INT(sqrt(DBLE(lm-1)))     ! l has to start from *zero*
            ! this should be the more efficient way to threat nspin>1 cases
            ! with minimal overhead for nspin=1 and no duplicated code
            pref = fpi/(2*l+1)
            aux(:) = pref * rho_lm(:,lm,1)
            DO ispin = 2,nspin ! if nspin < 2 it jumps to ENDDO
                aux(:) = aux(:)+fpi/(2*l+1)*rho_lm(:,lm,ispin)
            ENDDO
            CALL hartree(l, 2*l+2, rgrid(nt)%mesh, rgrid(nt), aux(:), pot_lm(:,lm))
            !
            ! now energy is computed integrating v_h^{lm} * \sum_{spin} rho^{lm}
            ! and summing on lm, aux already contains \sum_{spin} rho^{lm}
            ! but I have to redivide by 4pi/(2l+1)
            aux(:) = aux(:) * pot_lm(:,lm)
            CALL simpson (rgrid(nt)%mesh,aux,rgrid(nt)%rab,par_energy)
            !
            PAW_h_energy = PAW_h_energy + par_energy / pref
            IF (present(e_lm)) e_lm(lm) = par_energy / pref
    ENDDO ! lm

    CALL stop_clock ('PAW_h_energy')

END FUNCTION PAW_h_energy
!___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___
!!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!
!!! sum up pfuncs x occupation to build radial density's angular momentum components
SUBROUTINE PAW_rho_lm(na, becsum, pfunc, rho_lm, aug)
    USE kinds,                  ONLY : DP
    USE atom,                   ONLY : msh
    USE ions_base,              ONLY : ityp, ntyp => nsp, nat 
    USE lsda_mod,               ONLY : nspin
    USE uspp_param,             ONLY : nhm, nh, lmaxq
    USE uspp,                   ONLY : indv, ap, nhtolm,lpl,lpx
    USE parameters,             ONLY : nbrx, lqmax
    USE radial_grids,           ONLY : ndmx
    USE grid_paw_variables,     ONLY : augfun_t

    INTEGER,  INTENT(IN)  :: na     ! index of atom to use
    REAL(DP), INTENT(IN)  :: becsum(nhm*(nhm+1)/2,nat,nspin)! cross band occupation
    REAL(DP), INTENT(IN)  :: pfunc(ndmx,nbrx,nbrx,ntyp)     ! psi_i * psi_j
    REAL(DP), INTENT(OUT) :: rho_lm(ndmx,lmaxq**2,nspin)    ! AE charge density on rad. grid
!    REAL(DP), OPTIONAL,INTENT(IN):: augfun(ndmx,nbrx,nbrx,0:lqmax,ntyp)! augmentation functions (only for PS part)
    TYPE(augfun_t), OPTIONAL,INTENT(IN) :: &
                             aug(ntyp) ! augmentation functions (only for PS part)

    REAL(DP)                :: pref ! workspace (ap*becsum)

    INTEGER                 :: nb,mb, &     ! counters for pfunc nb,mb = 1, nh
                               nmb, &       ! composite "triangular" index for pfunc nmb = 1,nh*(nh+1)/2
                               lm,lp,l, &   ! counters for angular momentum lm = l**2+m
                               ispin,&      ! counter for spin (FIXME: may be unnecessary)
                               nt           ! type of atom na

    ! This subroutine computes the angular momentum components of rho
    ! using the following formula:
    !   rho(\vec{r}) = \sum_{LM} Y_{LM} \sum_{i,j} (\hat{r}) a_{LM}^{(lm)_i(lm)_j} becsum_ij pfunc_ij(r)
    !
    ! actually different angular momentum components are stored separately:
    !   rho^{LM}(\vec{r}) = \sum_{i,j} (\hat{r}) a_{LM}^{(lm)_i(lm)_j} becsum_ij pfunc_ij(r)
    !
    ! notice that pfunc's are already multiplied by r^2 and they are indexed on the atom
    ! (they only depends on l, not on m), the augmentation charge depend only on l
    ! but the becsum depend on both l and m

    CALL start_clock ('PAW_rho_lm')

    ! get type of atom
    nt = ityp(na)

    ! initialize density
    rho_lm(:,:,:) = 0._dp

    spins: DO ispin = 1, nspin
    nmb = 0
        ! loop on all pfunc for this kind of pseudo
        DO nb = 1, nh(nt)
        DO mb = nb, nh(nt)
            nmb = nmb+1 ! nmb = 1, nh*(nh+1)/2
            angular_momentum: &
            DO lp = 1, lpx (nhtolm(mb,nt), nhtolm(nb,nt)) !lmaxq**2
                ! the lpl array contains the possible combination of LM,lm_j,lm_j that
                ! have non-zero a_{LM}^{(lm)_i(lm)_j} (it saves some loops)
                lm = lpl (nhtolm(mb,nt), nhtolm(nb,nt), lp)
                ! becsum already contains a factor 2 for off-diagonal pfuncs
                pref = becsum(nmb,na,ispin) * ap(lm, nhtolm(nb,nt), nhtolm(mb,nt))
                !
                rho_lm(1:msh(nt),lm,ispin) = rho_lm(1:msh(nt),lm,ispin) +&
                                pref * pfunc(1:msh(nt), indv(nb,nt), indv(mb,nt), nt)
                IF (present(aug)) THEN
                    ! if I'm doing the pseudo part I have to add the augmentation charge
                    l = INT(sqrt(DBLE(lm-1))) ! l has to start from zero
                    rho_lm(1:msh(nt),lm,ispin) = rho_lm(1:msh(nt),lm,ispin) +&
                                pref * aug(nt)%fun(1:msh(na), indv(nb,nt), indv(mb,nt), l)
                ENDIF ! augfun
            ENDDO angular_momentum 
        ENDDO !mb
        ENDDO !nb
    ENDDO spins

    CALL stop_clock ('PAW_rho_lm')

END SUBROUTINE PAW_rho_lm
!___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___
!!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!
!!! build radial charge distribution from its spherical harmonics expansion
SUBROUTINE PAW_lm2rad(ix, rho_lm, rho_rad)
    USE kinds,                  ONLY : DP
    USE constants,              ONLY : eps8, pi
    USE uspp_param,             ONLY : lmaxq
    USE lsda_mod,               ONLY : nspin
    USE radial_grids,           ONLY : ndmx

    REAL(DP), INTENT(IN)        :: rho_lm(ndmx,lmaxq**2,nspin)! Y_lm expansion of rho
    INTEGER                     :: ix ! line of the ylm matrix to use
                                      ! actually it is one of the nx directions
    REAL(DP), INTENT(OUT)       :: rho_rad(ndmx,nspin) ! charge density on rad. grid
    !
    INTEGER                     :: ispin, lm ! counters on angmom and spin

    CALL start_clock ('PAW_lm2rad')
    rho_rad(:,:) = 0._dp
    ! cycling on spin is a bit less general...
    spins: DO ispin = 1,nspin
        DO lm = 1, lmaxq**2 ! 
            rho_rad(:,ispin) = rho_rad(:,ispin) +&
                    ylm(ix,lm)*rho_lm(:,lm,ispin)
        ENDDO ! lm
    ENDDO spins

    CALL stop_clock ('PAW_lm2rad')

END SUBROUTINE PAW_lm2rad

! REMOVE ME:
! #include "rad_paw_trash.f90"
!---------------------------------------------------------------
function exc_t(rho,rhoc,lsd)
  !---------------------------------------------------------------
  !
  use kinds, only : DP
  use funct
  implicit none
  integer:: lsd
  real(DP) :: exc_t, rho(2),arho,rhot, zeta,rhoc
  real(DP) :: ex, ec, vx(2), vc(2)

  real(DP),parameter:: e2 =2.0_DP

  exc_t=0.0_DP

  if(lsd == 0) then
     !
     !     LDA case
     !
     rhot = rho(1) + rhoc
     arho = abs(rhot)
     if (arho.gt.1.e-30_DP) then      
 !subroutine xc (rho, ex, ec, vx, vc)
        call xc(arho,ex,ec,vx(1),vc(1))
        exc_t=e2*(ex+ec)
     endif
  else
     !
     !     LSDA case
     !
     rhot = rho(1)+rho(2)+rhoc
     arho = abs(rhot)
     if (arho.gt.1.e-30_DP) then      
        zeta = (rho(1)-rho(2)) / arho
        call xc_spin(arho,zeta,ex,ec,vx(1),vx(2),vc(1),vc(2))
        exc_t=e2*(ex+ec)
     endif
  endif

  exc_t=e2*(ex+ec)

  return
end function exc_t
!
! Copyright (C) 2004 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!---------------------------------------------------------------
subroutine vxcgc(ndm,mesh,nspin,r,r2,rho,rhoc,vgc,egc)
  !---------------------------------------------------------------
  !
  !
  !     This routine compute the exchange and correlation potential and
  !     energy to be added to the local density, to have the first
  !     gradient correction.
  !     In input the density is rho(r) (multiplied by 4*pi*r2).
  !
  !     The units of the potential are Ryd.
  !
  use kinds, only : DP
  use funct
  implicit none
  integer :: ndm,mesh,nspin,ndm1
  real(DP) :: r(mesh), r2(mesh), rho(ndm,2), rhoc(ndm), &
       vgc(ndm,2), egc(ndm)

  integer :: i, is, ierr
  real(DP) :: sx,sc,v1x,v2x,v1c,v2c,gaux
  real(DP) :: v1xup, v1xdw, v2xup, v2xdw, v1cup, v1cdw
  real(DP) :: segno, arho, grho2(2)
  real(DP) :: rh, zeta, grh2
  real(DP),parameter :: eps=1.e-12_dp, fourpi=3.14159265358979_DP*4.0_DP

  real(DP), pointer :: grho(:,:), h(:,:), dh(:)
  !
  !      First compute the charge and the charge gradient, assumed  
  !      to have spherical symmetry. The gradient is the derivative of
  !      the charge with respect to the modulus of r. The last point is
  !      assumed to have zero gradient as happens in an atom.
  !
  allocate(grho(mesh,2),stat=ierr)
  allocate(h(mesh,2),stat=ierr)
  allocate(dh(mesh),stat=ierr)

  egc=0.0_dp
  vgc=0.0_dp

  do is=1,nspin
     do i=1, mesh
        rho(i,is)=(rho(i,is)+rhoc(i)/nspin)/fourpi/r2(i)
     enddo
     do i=2, mesh-1
        grho(i,is)=( (r(i+1)-r(i))**2*(rho(i-1,is)-rho(i,is)) &
             -(r(i-1)-r(i))**2*(rho(i+1,is)-rho(i,is)) )   &
             /((r(i+1)-r(i))*(r(i-1)-r(i))*(r(i+1)-r(i-1)))
     enddo
     grho(mesh,is)=0.0_dp
     !     
     !     The gradient in the first point is a linear interpolation of the
     !     gradient at point 2 and 3. The final result is not really sensitive to
     !     the value of these derivatives.
     !     
     grho(1,is)=grho(2,is)+(grho(3,is)-grho(2,is)) &
          *(r(1)-r(2))/(r(3)-r(2))
  enddo

  if (nspin.eq.1) then
     !
     !     GGA case
     !
     do i=1,mesh
        arho=abs(rho(i,1)) 
        segno=sign(1.0_dp,rho(i,1))
        if (arho.gt.eps.and.abs(grho(i,1)).gt.eps) then
           call gcxc(arho,grho(i,1)**2,sx,sc,v1x,v2x,v1c,v2c)
           egc(i)=(sx+sc)*segno
           vgc(i,1)= v1x+v1c
           h(i,1)  =(v2x+v2c)*grho(i,1)*r2(i)
           !            if (i.lt.4) write(6,'(f20.12,e20.12,2f20.12)') &
           !                          rho(i,1), grho(i,1)**2,  &
           !                          vgc(i,1),h(i,1)
        else if (i.gt.mesh/2) then
           !
           ! these are asymptotic formulae (large r) 
           !
           vgc(i,1)=-1.0_dp/r2(i)
           egc(i)=-0.0_dp/(2.0_dp*r(i))
           h(i,1)=h(i-1,1)
        else
           vgc(i,1)=0.0_dp
           egc(i)=0.0_dp
           h(i,1)=0.0_dp
        endif
     end do

  else
     !
     !   this is the \sigma-GGA case
     !       
     do i=1,mesh
        !
        !  NB: the special or wrong cases where one or two charges 
        !      or gradients are zero or negative must
        !      be detected within the gcxc_spin routine
        !
        !            call gcxc_spin(rho(i,1),rho(i,2),grho(i,1),grho(i,2),  &
        !                           sx,sc,v1xup,v1xdw,v2xup,v2xdw,          &
        !                           v1cup,v1cdw,v2c)
        !
        !    spin-polarised case
        !
        do is = 1, nspin
           grho2(is)=grho(i,is)**2
        enddo

        call gcx_spin (rho(i, 1), rho(i, 2), grho2(1), grho2(2), &
             sx, v1xup, v1xdw, v2xup, v2xdw)
        rh = rho(i, 1) + rho(i, 2)
        if (rh.gt.eps) then
           zeta = (rho (i, 1) - rho (i, 2) ) / rh
           grh2 = (grho (i, 1) + grho (i, 2) ) **2 
           call gcc_spin (rh, zeta, grh2, sc, v1cup, v1cdw, v2c)
        else
           sc = 0.0_dp
           v1cup = 0.0_dp
           v1cdw = 0.0_dp
           v2c = 0.0_dp
        endif

        egc(i)=sx+sc
        vgc(i,1)= v1xup+v1cup
        vgc(i,2)= v1xdw+v1cdw
        h(i,1)  =((v2xup+v2c)*grho(i,1)+v2c*grho(i,2))*r2(i)
        h(i,2)  =((v2xdw+v2c)*grho(i,2)+v2c*grho(i,1))*r2(i)
        !            if (i.lt.4) write(6,'(f20.12,e20.12,2f20.12)') &
        !                          rho(i,1)*2.0_dp, grho(i,1)**2*4.0_dp, &
        !                          vgc(i,1),  h(i,2)
     enddo
  endif
  !     
  !     We need the gradient of h to calculate the last part of the exchange
  !     and correlation potential.
  !     
  do is=1,nspin
     do i=2,mesh-1
        dh(i)=( (r(i+1)-r(i))**2*(h(i-1,is)-h(i,is))  &
             -(r(i-1)-r(i))**2*(h(i+1,is)-h(i,is)) ) &
             /( (r(i+1)-r(i))*(r(i-1)-r(i))*(r(i+1)-r(i-1)) )
     enddo

     dh(1)=dh(2)+(dh(3)-dh(2)) &
          *(r(1)-r(2))/(r(3)-r(2))
     dh(mesh)=0.0_dp
     !
     !     Finally we compute the total exchange and correlation energy and
     !     potential. We put the original values on the charge and multiply
     !     by two to have as output Ry units.

     do i=1, mesh
        vgc(i,is)=vgc(i,is)-dh(i)/r2(i)
        rho(i,is)=rho(i,is)*fourpi*r2(i)-rhoc(i)/nspin
        vgc(i,is)=2.0_dp*vgc(i,is)
        if (is.eq.1) egc(i)=2.0_dp*egc(i)
        !            if (is.eq.1.and.i.lt.4) write(6,'(3f20.12)') &
        !                                      vgc(i,1)
     enddo
  enddo

  deallocate(dh)
  deallocate(h)
  deallocate(grho)

  return
end subroutine vxcgc
!
! Copyright (C) 2004 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!---------------------------------------------------------------
subroutine vxc_t(rho,rhoc,lsd,vxc)
  !---------------------------------------------------------------
  !
  !  this function returns the XC potential in LDA or LSDA approximation
  !

  use kinds, only : DP
  use funct
  implicit none
  integer:: lsd
  real(DP):: vxc(2), rho(2),rhoc,arho,zeta
  real(DP):: vx(2), vc(2), ex, ec
  !
  real(DP), parameter :: e2=2.0_dp, eps=1.e-30_dp

  vxc(1)=0.0_dp
  if (lsd.eq.1) vxc(2)=0.0_dp

  if (lsd.eq.0) then
     !
     !     LDA case
     !
     arho=abs(rho(1)+rhoc)
     if (arho.gt.eps) then      
        call xc(arho,ex,ec,vx(1),vc(1))
        vxc(1)=e2*(vx(1)+vc(1))
     endif
  else
     !
     !     LSDA case
     !
     arho = abs(rho(1)+rho(2)+rhoc)
     if (arho.gt.eps) then      
        zeta = (rho(1)-rho(2)) / arho
        if (abs(zeta).gt.1.0_dp) then 
           write(6,*) 'zeta= me', zeta, rho(1),rho(2),rhoc
        else
           call xc_spin(arho,zeta,ex,ec,vx(1),vx(2),vc(1),vc(2))
           vxc(1) = e2*(vx(1)+vc(1))
           vxc(2) = e2*(vx(2)+vc(2))
        endif
     endif
  endif

  return
end subroutine vxc_t

END MODULE rad_paw_routines
