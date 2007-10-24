!     WRITE(6,*) "==PAW RADIAL ENERGIES: "
!     WRITE(6,*) "==AE 1     :", e(1,1)
!     WRITE(6,*) "==AE 2     :", e(2,1)
!     WRITE(6,*) "==PS 1     :", e(1,2)
!     WRITE(6,*) "==PS 2     :", e(2,2)
!     WRITE(6,*) "==AE tot   :", SUM(e(:,1))
!     WRITE(6,*) "==PS tot   :", SUM(e(:,2))
!     WRITE(6,*) "==AE-PS 1  :", e(1,1)-e(1,2)
!     WRITE(6,*) "==AE-PS 2  :", e(2,1)-e(2,2)
!     WRITE(6,"(a,2f15.7)") "==AE-PS tot:", SUM(e(:,1))-SUM(e(:,2))
!     WRITE(6,"(a,2f15.7)") "== LM=1:",ecomps(1,1,1,1),ecomps(1,1,1,2)
!     WRITE(6,"(a,2f15.7)") "== LM=2:",ecomps(2,1,1,1),ecomps(2,1,1,2)
!     WRITE(6,"(a,2f15.7)") "== LM=3:",ecomps(3,1,1,1),ecomps(3,1,1,2)
!     WRITE(6,"(a,2f15.7)") "== LM=4:",ecomps(4,1,1,1),ecomps(4,1,1,2)
!     WRITE(6,"(a,2f15.7)") "== LM=5:",ecomps(5,1,1,1),ecomps(5,1,1,2)
!     WRITE(6,"(a,2f15.7)") "== LM=6:",ecomps(6,1,1,1),ecomps(6,1,1,2)
!     WRITE(6,"(a,2f15.7)") "== LM=7:",ecomps(7,1,1,1),ecomps(7,1,1,2)
!     WRITE(6,"(a,2f15.7)") "== LM=8:",ecomps(8,1,1,1),ecomps(8,1,1,2)
!     WRITE(6,"(a,2f15.7)") "== LM=9:",ecomps(9,1,1,1),ecomps(9,1,1,2)
!     WRITE(6,"(a,2f15.7)") "=============================================="
!___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___
!!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!
!!! This is the main driver of PAW routines, it uses directly or indirectly
!!! all the other routines of the module
!!
SUBROUTINE PAW_energy(becsum,e_atom)
    USE kinds,                  ONLY : DP
    USE lsda_mod,               ONLY : nspin
    USE ions_base,              ONLY : nat, ityp
    USE atom,                   ONLY : g => rgrid
    USE grid_paw_variables,     ONLY : pfunc, ptfunc, tpawp, &
                                       aerho_atc, psrho_atc, aug
    USE uspp_param,             ONLY : nhm, lmaxq

    REAL(DP), INTENT(IN)    :: becsum(nhm*(nhm+1)/2,nat,nspin)! cross band occupations
    REAL(DP), INTENT(OUT),OPTIONAL :: &
                               e_atom(nat,2,2) ! {# of atoms}, {XC|H}, {AE|PS}

    INTEGER, PARAMETER      :: AE = 1, PS = 2,& ! All-Electron and Pseudo
                               XC = 1, H  = 2   ! XC and Hartree
    INTEGER                 :: i_what           ! counter on AE and PS
    INTEGER                 :: na,first_nat,last_nat ! atoms counters and indexes

    ! hartree energy scalar fields expanded on Y_lm
    REAL(DP), ALLOCATABLE   :: rho_lm(:,:,:)      ! radial density expanded on Y_lm
!    REAL(DP), ALLOCATABLE   :: v_h_lm(:,:)        ! hartree potential, summed on spins
    !
    REAL(DP)                :: e_h(lmaxq**2,nspin)! hartree energy components
    REAL(DP)                :: e,e1,e2            ! placeholders
    ! xc variables:
    REAL(DP)                :: e_xc(lmaxq**2,nspin)! UNUSED! XC energy components
    REAL(DP), POINTER       :: rho_core(:,:)       ! pointer to AE/PS core charge density 
    TYPE(paw_info)          :: i                   ! minimal info on atoms

    ! BEWARE THAT HARTREE ONLY DEPENDS ON THE TOTAL RHO NOT ON RHOUP AND RHODW SEPARATELY...
    ! TREATMENT OF NSPIN>1 MUST BE CHECKED AND CORRECTED

    CALL start_clock ('PAW_energy')

    ! nullify energy components:
    IF( present(e_atom) ) e_atom(:,:,:) = 0._dp
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
    ! Operations from here on are (will be) parallelized on atoms, for the
    ! moment this OpenMP directive gives single-core parallelization:
!$OMP PARALLEL DO default(shared) reduction(+:e_atom) &
!$OMP private(rho_lm, e_xc, e_h, e,e1,e2, rho_core, na, i_what, i)
    atoms: DO na = first_nat, last_nat
    !
    i%a = na         ! the index of the atom
    i%t = ityp(na)   ! the type of atom na
    i%m = g(i%t)%mesh! radial mesh size for atom na
    !
    ifpaw: IF (tpawp(i%t)) THEN
        !
        ! Arrays are allocated inside the cycle to allow reduced
        ! memory usage as differnt atoms have different meshes
        ALLOCATE(rho_lm(i%m,lmaxq**2,nspin))
        !
        whattodo: DO i_what = AE, PS
            i%w = i_what     ! spherical_gradient likes to know
            ! STEP: 1 [ build rho_lm (PAW_rho_lm) ]
            NULLIFY(rho_core)
            IF (i_what == AE) THEN
                ! to pass the atom indes is dirtyer but faster and
                ! uses less memory than to pass only a hyperslice of the array
                CALL PAW_rho_lm(i, becsum, pfunc, rho_lm)
                ! used later for xc energy:
                rho_core => aerho_atc
            ELSE
                CALL PAW_rho_lm(i, becsum, ptfunc, rho_lm, aug)
                !     optional argument for pseudo part --> ^^^
                ! used later for xc energy:
                rho_core => psrho_atc
            ENDIF
            ! STEP: 2 [ compute Hartree energy ]
            !   2a. use rho_lm to compute hartree potential (PAW_v_h)
            e = PAW_h_energy(i, rho_lm, e_h)
            IF( present(e_atom) ) e_atom(i%a, H, i_what) = e

            ! STEP: 3 [ compute XC energy ]
            e1 = PAW_xc_energy(i, rho_lm, rho_core, e_xc,0)
            IF( present(e_atom) ) e_atom(i%a, XC, i_what) = e1
            e2 = PAW_xc_energy(i, rho_lm, rho_core, e_xc,1)
            write(*,"(a,2i2,2f10.5,a)") " == ", i_what, na, e, e1,"  (2)"

        ENDDO whattodo


        ! cleanup
        DEALLOCATE(rho_lm)

    ENDIF ifpaw
    ENDDO atoms
!$OMP END PARALLEL DO

    CALL stop_clock ('PAW_energy')


END SUBROUTINE PAW_energy

!___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___
!!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!! 
!!! add gradient correction to energy, code mostly adapted from ../atomic/vxcgc.f90
SUBROUTINE PAW_gcxc_energy(i, rho,core,grho,e)
    USE kinds,                  ONLY : DP
    USE ions_base,              ONLY : ityp
    USE lsda_mod,               ONLY : nspin
    USE atom,                   ONLY : g => rgrid
    USE parameters,             ONLY : ntypx
    USE constants,              ONLY : fpi,pi,e2
    USE funct,                  ONLY : gcxc, gcx_spin, gcc_spin
    !
    TYPE(paw_info)  :: i                     ! atom's minimal info
    REAL(DP), INTENT(IN)   :: rho(i%m,nspin) ! radial density,
    REAL(DP), INTENT(IN)   :: grho(i%m,nspin)! square gradient of rho
    REAL(DP), INTENT(IN)   :: core(i%m,ntypx) ! spherical core density
    REAL(DP), INTENT(INOUT):: e(i%m)         ! radial local xc energy
    !
    INTEGER            :: k               ! counter for mesh
    ! workspaces:
    REAL(DP)           :: arho, sgn
    REAL(DP)           :: sx,sc           ! x and c energy from gcxc
    REAL(DP)           :: v1x,v2x,v1c,v2c ! potentials from gcxc (unused)
    REAL(DP),PARAMETER :: eps = 1.D-12    ! 1.e-12 may be small enough
    ! for nspin>1
    REAL(DP)           :: rh,grh2,zeta
    REAL(DP)           :: v1cup, v1cdw,v1xup, v1xdw, v2xup, v2xdw !(unused)
    
    OPTIONAL_CALL start_clock ('PAW_gcxc_e')

    IF (nspin.eq.1) THEN
        DO k=1,i%m
            arho = rho(k,1)*g(i%t)%rm2(k) + core(k,i%t)
            sgn  = sign(1._dp,arho)
            arho = abs(arho)
            IF (arho.gt.eps.and.abs(grho(k,1)).gt.eps) THEN
                CALL gcxc(arho,grho(k,1),sx,sc,v1x,v2x,v1c,v2c)
                e(k) = e(k) + sgn *e2 *(sx+sc)*g(i%t)%r2(k) !&
            ENDIF
        ENDDO
    ELSE
        !   this is the \sigma-GGA case
        DO k=1,i%m
            CALL gcx_spin (rho(k,1), rho(k,2), grho(k,1), grho(k,2), &
                           sx, v1xup, v1xdw, v2xup, v2xdw)
            rh = rho(k,1) + rho(k,2)
            IF (rh.gt.eps) THEN
                zeta = (rho (k,1) - rho (k,2) ) / rh
                ! FIXME: next line is wrong because it looses sign!!
                grh2 = (sqrt(grho(k,1)) + sqrt(grho(k,2)) ) **2 
                CALL gcc_spin (rh, zeta, grh2, sc, v1cup, v1cdw, v2c)
            ELSE
                sc = 0._dp
!                 v1cup = 0.0_dp
!                 v1cdw = 0.0_dp
!                 v2c = 0.0_dp
            ENDIF
            e(k) = e(k) + e2*(sx+sc)
!             vgc(i,1)= v1xup+v1cup
!             vgc(i,2)= v1xdw+v1cdw
!             h(i,1)  =((v2xup+v2c)*grho(i,1)+v2c*grho(i,2))*r2(i)
!             h(i,2)  =((v2xdw+v2c)*grho(i,2)+v2c*grho(i,1))*r2(i)
        ENDDO
    ENDIF

    OPTIONAL_CALL stop_clock ('PAW_gcxc_e')

END SUBROUTINE PAW_gcxc_energy
!___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___
!!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!! 
!!! use the density produced by sum_rad_rho to compute hartree potential 
!!! the potential is then directly integrated to compute hartree energy
FUNCTION PAW_h_energy(i, rho_lm, e_lm)
    USE kinds,                  ONLY : DP
    USE constants,              ONLY : fpi
    USE parameters,             ONLY : ntypx
    USE lsda_mod,               ONLY : nspin
    USE uspp_param,             ONLY : nh, lmaxq
    USE atom,                   ONLY : g => rgrid

    REAL(DP)                       :: PAW_h_energy      ! total hartree energy
    !
    TYPE(paw_info)  :: i                                ! the number of the atom
    REAL(DP), INTENT(IN)  :: rho_lm(i%m,lmaxq**2,nspin) ! charge density as lm components
    REAL(DP), OPTIONAL,INTENT(OUT) :: e_lm(lmaxq**2)    ! out: energy components
    !
    REAL(DP)              :: aux(i%m)              ! workspace
    REAL(DP)              :: par_energy            ! workspace
    REAL(DP)              :: v_lm(i%m,lmaxq**2)    ! potential as lm components
    REAL(DP)              :: trho_lm(i%m,lmaxq**2) !  total spin-summed charge density

    INTEGER               :: ispin, &    ! counter on spins
                             lm,l        ! counter on composite angmom lm = l**2 +m
    OPTIONAL_CALL start_clock ('PAW_h_energy')

    ! init total energy and its lm components
    PAW_h_energy = 0._dp
    IF (present(e_lm)) e_lm(:) = 0._dp

    ! sum up spin components, respecting column-majorness
    trho_lm(:,:) = rho_lm(:,:,1)
    DO ispin = 2,nspin ! if nspin < 2 it jumps to ENDDO
        trho_lm(:,:) = trho_lm(:,:) + rho_lm(:,:,ispin)
    ENDDO

    ! compute the hartree potential
    CALL PAW_h_potential(i, rho_lm, v_lm)

    DO lm = 1, lmaxq**2
            !
            ! now energy is computed integrating v_h^{lm} * \sum_{spin} rho^{lm}
            ! and summing on lm 
            aux(:) = trho_lm(:,lm) * v_lm(:,lm)
            CALL simpson (i%m, aux, g(i%t)%rab, par_energy)
            !
            PAW_h_energy = PAW_h_energy + par_energy
            IF (present(e_lm)) e_lm(lm) = par_energy
    ENDDO ! lm

    OPTIONAL_CALL stop_clock ('PAW_h_energy')

END FUNCTION PAW_h_energy

!___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___   ___
!!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!!  !!!! 
!!! use the density produced by sum_rad_rho to compute xc potential and energy, as
!!! xc functional is not diagonal on angular momentum numerical integration is performed
FUNCTION PAW_xc_energy(i, rho_lm, rho_core, e_lm, task)
    USE kinds,                  ONLY : DP
    USE constants,              ONLY : fpi
    USE parameters,             ONLY : ntypx
    USE lsda_mod,               ONLY : nspin
    USE uspp_param,             ONLY : lmaxq
    USE ions_base,              ONLY : ityp
    USE radial_grids,           ONLY : ndmx
    USE atom,                   ONLY : g => rgrid

    REAL(DP)              :: PAW_xc_energy      ! total xc energy
    !
    TYPE(paw_info)  :: i                               ! atom's minimal info
    INTEGER,  INTENT(IN)  :: task!remove me
    REAL(DP), INTENT(IN)  :: rho_lm(i%m,lmaxq**2,nspin)! charge density as lm components
    REAL(DP), INTENT(IN)  :: rho_core(ndmx,ntypx)       ! core charge, radial and spherical
    ! TODO:
    REAL(DP), OPTIONAL,INTENT(OUT) :: e_lm(lmaxq**2)      ! out: energy components 
    !
    INTEGER               :: lsd          ! switch to control local spin density
    !
    REAL(DP)              :: rho_loc(2) = (/0._dp, 0._dp/) 
                             ! local density (workspace), up and down
    REAL(DP)              :: e            ! workspace
    REAL(DP)              :: e_rad(i%m)   ! radial energy (to be integrated)
    REAL(DP)              :: rho_rad(i%m,nspin) ! workspace (radial slice of rho)
    INTEGER               :: ix,k         ! counters on directions and radial grid
    ! for gradient correction only:
    REAL(DP),ALLOCATABLE  :: grho_rad(:,:)! workspace (radial slice of grad(rho))
    ! functions from atomic code:
    REAL(DP),EXTERNAL     :: exc_t

    ! to prevent warnings (this quantity is not implemented,
    ! and maybe it never will):
    e_lm   = 0._dp

    OPTIONAL_CALL start_clock ('PAW_xc_nrg')
    lsd = nspin-1

    ! init for gradient correction
    IF (do_gcxc) ALLOCATE(grho_rad(i%m,nspin))

    PAW_xc_energy = 0._dp
    DO ix = 1, nx
        ! LDA (and LSDA) part (no gradient correction):
        CALL PAW_lm2rad(i, ix, rho_lm, rho_rad)
        !
        DO k = 1,i%m
            rho_loc(1:nspin) = rho_rad(k,1:nspin)*g(i%t)%rm2(k)
            !
            e_rad(k) = exc_t(rho_loc, rho_core(k,i%t), lsd)&
                     * (SUM(rho_rad(k,1:nspin))+rho_core(k,i%t)*g(i%t)%r2(k))
        ENDDO
        gradient_correction:& ! add it!
        IF (do_gcxc) THEN
            IF (task == 1) e_rad(:) = 0._dp ! reset the energy, so that only the correction is displayed
            CALL PAW_gradient(i, ix, rho_lm, rho_rad, rho_core, grho_rad)
            !                                          v-------------^
            CALL PAW_gcxc_energy(i, rho_rad,rho_core,grho_rad, e_rad)
        ENDIF gradient_correction
        !
        ! integrate radial slice of xc energy:
        CALL simpson (i%m, e_rad, g(i%t)%rab, e)
        ! integrate on sph. surface     v-----^
        PAW_xc_energy = PAW_xc_energy + e * ww(ix)
    ENDDO

    ! cleanup for gc
    IF (do_gcxc) DEALLOCATE(grho_rad)

    OPTIONAL_CALL stop_clock ('PAW_xc_nrg')

END FUNCTION PAW_xc_energy

SUBROUTINE PAW_brute_radial_ddd(becsum)
    USE kinds,                  ONLY : DP
    USE lsda_mod,               ONLY : nspin
    USE ions_base,              ONLY : nat, ityp
    USE uspp_param,             ONLY : nhm, nh
    USE grid_paw_variables,     ONLY : tpawp

    REAL(DP), INTENT(IN)    :: becsum(nhm*(nhm+1)/2,nat,nspin)! cross band occupations
    REAL(DP)     :: ddd(nhm,nhm,nat,nspin,1:2)! cross band occupations

    INTEGER, PARAMETER      :: AE = 1, PS = 2, AEmPS=0 ! All-Electron and Pseudo
    INTEGER  :: ih, jh, ijh
    INTEGER  :: na, nt, ispin
    !
    REAL(DP), PARAMETER :: eps = 1.D-6
    !
    REAL(DP) :: b1,b2,delta
    REAL(DP) :: bectmp(nhm*(nhm+1)/2,nat,nspin)
    REAL(DP) :: e(nat,2,2), ee ! {# of atoms}, {XC|H}, {AE|PS}
    
    CALL start_clock('PAW_brddd')
    !
    ddd = 0._dp
    ispin = 1 ! FIXME
    bectmp = becsum

    DO na = 1,nat
    nt = ityp(na)
    IF ( tpawp(nt) ) THEN
        !
        DO ih = 1, nh(nt)
        DO jh = ih, nh(nt)
!             write(6,"(a,3i4)") "now doing:",ih,jh,ijh
            !
            ijh = jh * (jh-1) / 2 + ih
            !
            b1 = becsum(ijh,na,ispin) + eps
            b2 = becsum(ijh,na,ispin) - eps
            delta = (b1-b2)
            !
            bectmp(ijh,na,ispin) = b1
            CALL PAW_potential(bectmp, ee, e)
            ddd(ih,jh,na,ispin,AE) = - SUM(e(na,:,AE))
            ddd(ih,jh,na,ispin,PS) = - SUM(e(na,:,PS))
            !
            bectmp(ijh,na,ispin) = b2
            CALL PAW_potential(bectmp,ee, e)
            ddd(ih,jh,na,ispin,AE) = ddd(ih,jh,na,ispin,AE) &
                                    + SUM(e(na,:,AE))
            ddd(ih,jh,na,ispin,PS) = ddd(ih,jh,na,ispin,PS) &
                                    + SUM(e(na,:,PS))
            !
            ddd(ih,jh,na,ispin,:)  = ddd(ih,jh,na,ispin,:) / delta
!             write(6,*) ddd(ih,jh,na,ispin,:)
            ddd(jh,ih,na,ispin,AE) = ddd(ih,jh,na,ispin,AE)
            ddd(jh,ih,na,ispin,PS) = ddd(ih,jh,na,ispin,PS)
            !
            bectmp(ijh,na,ispin) = becsum(ijh,na,ispin)
            !
        END DO ! jh
        END DO ! ih
        !
    ENDIF
    ENDDO
    !
    CALL stop_clock('PAW_brddd')

    PRINT *, "RADIAL VERSION"
    PRINT *, 'D - D1'
    PRINT '(8f15.7)', ((ddd(jh,ih,1,1,AE),jh=1,nh(nt)),ih=1,nh(nt))
    PRINT *, 'D - D1~'
    PRINT '(8f15.7)', ((ddd(jh,ih,1,1,PS),jh=1,nh(nt)),ih=1,nh(nt))


END SUBROUTINE PAW_brute_radial_ddd


SUBROUTINE integrate_pfunc
    !
    USE kinds,      ONLY : DP
    USE parameters, ONLY : lmaxx, nbrx, lqmax
    USE radial_grids, ONLY: ndmx
    USE constants,  ONLY : fpi, eps8, eps4
    USE atom,       ONLY : r, rab, mesh, msh
    USE ions_base,  ONLY : ntyp => nsp
    USE cell_base,  ONLY : omega, tpiba
    USE gvect,      ONLY : g, gg
    USE lsda_mod,   ONLY : nspin
    USE us,         ONLY : nqxq, dq, nqx, tab, qrad
    USE uspp 
    USE uspp_param
    USE spin_orb,   ONLY : lspinorb, rot_ylm, fcoef
    !
    USE grid_paw_variables, ONLY: tpawp, pfunc, ptfunc, pp, ppt, prad, ptrad, okpaw
    ! for FFt method
    USE gvect,         ONLY : gg, gi =>gstart, ngm
    USE grid_paw_routines, ONLY: pvan2
    !
    IMPLICIT NONE
    !, int_pfunc_(:,:,:)
    REAL(DP), POINTER :: pfunc_(:,:,:,:), prad_(:,:,:,:), pp_(:,:,:), int_pfunc_(:,:,:,:,:)
    REAL(DP),TARGET   :: int_pfunc(nbrx,nbrx,nbrx,nbrx,ntyp),&
                         int_ptfunc(nbrx,nbrx,nbrx,nbrx,ntyp)
    REAL(DP)          :: integral, ap2
    !
    INTEGER :: i_what, terms
    REAL(DP) :: aux2(ndmx)
    !
    ! here a few local variables
    !
    INTEGER :: nt, ih, jh, nb, mb, nc, mc, nmb, l, m, lm, ir, iq, is, ndm ! various counters
    REAL(DP), ALLOCATABLE :: aux (:), aux1 (:)
    ! various work space
    INTEGER :: n1, m0, m1, n, li, mi, vi, vj, ijs, is1, is2, &
              lk, mk, vk, kh, lh, sph_ind, nnbrx, ll,j
    COMPLEX(DP) :: coeff, qgm(1)
    REAL(DP) :: ap_tot
    REAL(DP) :: spinor, ji, jk

    ! for FFT method
    REAL(DP), ALLOCATABLE :: qmod (:),  & ! the modulus of G
                             ylmk0 (:,:)  ! the spherical harmonics
    COMPLEX(DP), ALLOCATABLE  :: pft(:,:,:)  ! the AE/PS wfc products
    COMPLEX(DP), TARGET  :: int_pfunc_fft(nbrx,nbrx,nbrx,nbrx, ntyp),&  ! the AE/PS wfc products
                            int_ptfunc_fft(nbrx,nbrx,nbrx,nbrx, ntyp)
    COMPLEX(DP), POINTER :: int_pfunc_fft_(:,:,:,:,:)
    INTEGER :: terms_counter(9,9,9,9,2,1)

WRITE(6,*) "RADIAL PAW ROUTINES: integrate_pfunc (start)"

RETURN
    !
    ! part1: compute P_ij * P_ab on radial, real space, grid
    !--------------------------------------------------------------------------------

whattodo: DO i_what=1, 2
       ! associate a pointer to the AE or PS part
       NULLIFY(pfunc_,int_pfunc_)
       IF (i_what==1) THEN
          pfunc_=> pfunc
          int_pfunc_ => int_pfunc
       ELSE IF (i_what==2) THEN
          pfunc_=> ptfunc
          int_pfunc_ => int_ptfunc
       END IF
       ! Compute the integrals of pfunc
       DO nt = 1, ntyp ! ntype is the # of atomic species (PP's) is .le. than the # of atoms
          IF (tpawp(nt)) THEN
            ! I have to cicle on pfunc TWICE, and each pfunc has 2 indexes => 4 indexes
            ih = 0
            DO nc = 1, nh(nt)
            DO mc = 1, nh(nt)
                DO nb = 1, nh(nt)
                DO mb = 1, nh(nt)
!                    WRITE(6,*) MAXVAL(pfunc_(1:msh(nt), nb, mb, nt))
                    ih = ih+1
                    int_pfunc_(nc,mc,nb,mb,nt) = 0._DP
                    terms = 0
                    DO lm = 1, lmaxq**2 ! FIXME: is this the right upper bound??
                        ap2 = ap(lm,nhtolm(nc,nt),nhtolm(mc,nt))*ap(lm,nhtolm(nb,nt),nhtolm(mb,nt))
                        IF ( ABS(ap2) > eps8 ) THEN
                            terms = terms +1
                            ! if I don't have the augfun the integral have to be computed only once
                            IF ((i_what == 1) .and. (terms == 1)) THEN
                                aux2(1:msh(nt)) = (pfunc_(1:msh(nt), indv(nb,nt), indv(mb,nt), nt)/r(1:msh(nt),nt))*&
                                                  (pfunc_(1:msh(nt), indv(nc,nt), indv(mc,nt), nt)/r(1:msh(nt),nt))
                                CALL simpson (msh(nt),aux2,rab(1,nt),integral)
                            ENDIF
                            ! with the augfun than I have to compute the integral for each value of lm
                            IF ((i_what == 2)) THEN
                                l = INT(sqrt(DBLE(lm-1))) ! the "+1" is not required, as augfun are labelled 0..l
                                aux2(1:msh(nt)) = (pfunc_(1:msh(nt), indv(nb,nt), indv(mb,nt), nt)&
                                                    +augfun(1:msh(nt), indv(nb,nt), indv(mb,nt), l, nt))/r(1:msh(nt),nt)*&
                                                  (pfunc_(1:msh(nt), indv(nc,nt), indv(mc,nt), nt)&
                                                    +augfun(1:msh(nt), indv(nc,nt), indv(mc,nt), l, nt))/r(1:msh(nt),nt)
                                ! the following line is duplicated (better than using two IF..THEN)
                                CALL simpson (msh(nt),aux2,rab(1,nt),integral)
                            ENDIF
                            ! anyway I have to sum
                            int_pfunc_(nc,mc,nb,mb,nt) = int_pfunc_(nc,mc,nb,mb,nt) + integral*ap2
                        ENDIF
                    ENDDO ! l = 1, lmaxq
                    IF (terms > 0.and. i_what==2)&
                       WRITE(1001,"(i4,2i4,3i2,f14.8,i3,2f8.4)"), ih,i_what,nc,mc,nb,mb,  int_pfunc_(nc,mc,nb,mb,nt),&
                                   terms,integral, ap2
                    terms_counter(nc,mc,nb,mb,i_what,nt) = terms
                END DO !mb
                END DO !nb 
            END DO !mc
            END DO !nc 
            !
          END IF ! tpawp
       END DO ! nt
    END DO whattodo

    !
    ! part2: the same in FFT
    !--------------------------------------------------------------------------------
    WRITE(6,*) "done: radial"
    ALLOCATE (qmod(ngm), pft(ngm,nbrx,nbrx), ylmk0(ngm,lmaxq*lmaxq))
    !
    CALL ylmr2 (lmaxq * lmaxq, ngm, g, gg, ylmk0)
    qmod(:) = SQRT(gg(:))
    !
    whattodo2: DO i_what=1, 2
        !
        NULLIFY(prad_,int_pfunc_fft_)
        IF (i_what==1) THEN
            prad_ => prad
            int_pfunc_fft_ => int_pfunc_fft
        ELSE IF (i_what==2) THEN
            ! ***NOTE: ptrad already has the augmentation charge***
            prad_ => ptrad
            int_pfunc_fft_ => int_ptfunc_fft
        END IF
         pft(:,:,:) = 0._DP !probably unnecessary
        !
        int_pfunc_fft_ (:,:,:,:,:)  = (0.d0, 0.d0)
        !
        DO nt = 1, ntyp
            ih = 0
            pft (:,:,:) = (0.d0, 0.d0)
            IF (tpawp(nt)) THEN
                DO nc = 1, nh (nt)
                DO mc = 1, nh (nt)
                    CALL pvan2 (ngm, mc, nc, nt, qmod, pft(1,nc,mc), ylmk0, prad_, &
                        SIZE(prad_,1),SIZE(prad_,2),SIZE(prad_,3),SIZE(prad_,4))
                ENDDO ! jh
                ENDDO ! ih

                DO nc = 1, nh(nt)
                DO mc = 1, nh(nt)
                    DO nb = 1, nh(nt)
                    DO mb = 1, nh(nt)
                    ih = ih+1
                    int_pfunc_fft_ (mb,nb,mc,nc, nt) = OMEGA *& 
                    SUM( DBLE( CONJG(pft(:,mc,nc))*pft(:,mb,nb) ))
                    !
                    !int_pfunc_fft_ (ijh2, ijh, nt) = CONJG( int_pfunc_fft_ (ijh, ijh2, nt) )
                 !
                    IF (ABS(int_pfunc_fft_(nc,mc,nb,mb,nt))>eps8 .and. i_what==2) &
                        WRITE(1002,"(i4,2i4,3i2,2f14.8)"), ih,i_what,nc,mc,nb,mb,  int_pfunc_fft_(nc,mc,nb,mb,nt)
                    ENDDO ! mb
                    ENDDO ! nb
                ENDDO ! mc
                ENDDO ! nc
            ENDIF ! tpawp
        ENDDO ! nt
     !
    ENDDO whattodo2
    !
    DEALLOCATE (qmod, pft, ylmk0) 
    WRITE(6,*) "done: FFT"

    lm = 1
    whattodo3: DO i_what=1, 2
        !
        NULLIFY(prad_,int_pfunc_fft_,int_pfunc_)
        IF (i_what==1) THEN
            int_pfunc_ => int_pfunc
            int_pfunc_fft_ => int_pfunc_fft
        ELSE IF (i_what==2) THEN
            int_pfunc_ => int_ptfunc
            int_pfunc_fft_ => int_ptfunc_fft
        END IF
        DO nt = 1, ntyp
            ih = 0
            IF (tpawp(nt)) THEN
                DO nc = 1, nh(nt)
                DO mc = 1, nh(nt)
                    DO nb = 1, nh(nt)
                    DO mb = 1, nh(nt)
!                 DO nc = 1, nh(nt)
!                 DO mc = nc, nh(nt)
!                     DO nb = mc, nh(nt)
!                     DO mb = nb, nh(nt)
                        ih = ih+1
                        IF( ABS(int_pfunc_(nc,mc,nb,mb,nt) - int_pfunc_fft_(nc,mc,nb,mb,nt)) > eps4) THEN
                            WRITE(1003,"(3i4,3i2,f14.8,f16.8,f14.8,i3,4i2,a)"), ih,i_what,nc,mc,nb,mb, &
                             int_pfunc_(nc,mc,nb,mb,nt), int_pfunc_fft_(nc,mc,nb,mb,nt),terms_counter(nc,mc,nb,mb,i_what,nt),&
                             nhtolm(nc,nt),nhtolm(mc,nt),nhtolm(nb,nt),nhtolm(mb,nt),"yadda"
                        ELSE IF ( ABS(int_pfunc_(nc,mc,nb,mb,nt)) > eps8 ) THEN
                            WRITE(1004,"(3i4,3i2,f14.8,f16.8,f14.8,i3,4i2,a)"),ih,i_what,nc,mc,nb,mb, &
                             int_pfunc_(nc,mc,nb,mb,nt), int_pfunc_fft_(nc,mc,nb,mb,nt), terms_counter(nc,mc,nb,mb,i_what,nt),&
                             nhtolm(nc,nt),nhtolm(mc,nt),nhtolm(nb,nt),nhtolm(mb,nt),"blah"
                        ENDIF
                    ENDDO ! mb
                    ENDDO ! nb
                ENDDO ! mc
                ENDDO ! nc
            ENDIF ! tpawp
        ENDDO ! nt
    ENDDO whattodo3
WRITE(6,*) "RADIAL PAW ROUTINES: integrate_pfunc (end)"

STOP

END SUBROUTINE integrate_pfunc


! analogous to compute_onecenter_charges
!
  SUBROUTINE coc_pwned(becnew, rho1new, rho1tnew,lm)
    !
    USE kinds,                ONLY : DP
    USE constants,            ONLY : eps8
    USE ions_base,            ONLY : nat, ntyp => nsp, ityp, tau, atm
    USE gvect,                ONLY : nr1, nr2, nr3, nrx1, nrx2, nrx3, nrxx, &
                                     ngm, nl, nlm, gg, g
    USE lsda_mod,             ONLY : nspin
    USE uspp_param,           ONLY : lmaxq, nh, nhm
    USE wvfct,                ONLY : gamma_only
    USE wavefunctions_module, ONLY : psic
    !
    USE grid_paw_variables, ONLY: prad, ptrad, pp, tpawp, okpaw
    USE grid_paw_routines
    USE us,                 ONLY: qrad
    USE uspp,                   ONLY : indv, ap, nhtolm
    !
    IMPLICIT NONE
    !
    !first input-output variables
    ! 
    REAL(DP), INTENT(IN) :: becnew (nhm*(nhm+1)/2,nat,nspin)
    REAL(DP), TARGET, INTENT(OUT) :: &
         rho1new(nrxx, nspin, nat), rho1tnew(nrxx,nspin,nat)
    !
    INTEGER :: ig, na, nt, ih, jh, ijh, is ! counters
    INTEGER,INTENT(IN) :: lm!DEBUG
    !
    REAL(DP), ALLOCATABLE :: qmod (:), & ! the modulus of G
                             ylmk0 (:,:) ! the spherical harmonics
    COMPLEX(DP), ALLOCATABLE :: aux (:,:,:), & ! work space for rho(G,nspin)
                                qgm(:)         ! Fourier transform of q

    REAL(DP), POINTER :: rho1_(:,:,:), prad_(:,:,:,:)
    INTEGER :: i_what

    IF (.NOT.okpaw) RETURN
    call start_clock ('one-charge')

    ALLOCATE (aux(ngm,nspin,nat), qmod(ngm), qgm(ngm), ylmk0(ngm,lmaxq*lmaxq))    
    !  
    CALL ylmr2 (lmaxq * lmaxq, ngm, g, gg, ylmk0)
    qmod(:) = SQRT(gg(:))

    !WRITE(20,*) "becsum used in GRID:"
    DO i_what =1,2
        atoms: DO na = 1, nat
        nt = ityp(na)
            spins: DO is = 1, nspin
            ijh = 0
                ! loop on all pfunc for this kind of pseudo
                DO ih = 1, nh(nt)
                DO jh = ih, nh(nt)
                    ijh = ijh+1
                    !WRITE(20,"(a,i3,a,4i3,f12.6)") "-->",ijh,":",ih,jh,na,is,becnew(ijh,na,is)
                ENDDO
                ENDDO
            ENDDO spins
        ENDDO atoms
    ENDDO

    whattodo: DO i_what=1, 2
       NULLIFY(prad_,rho1_)
       IF (i_what==1) THEN
          prad_ => prad
          rho1_ => rho1new
       ELSE IF (i_what==2) THEN
          prad_ => ptrad
          rho1_ => rho1tnew
       END IF
       aux (:,:,:) = (0.d0, 0.d0)

       DO nt = 1, ntyp
          IF (tpawp (nt) ) THEN
             ijh = 0
             DO ih = 1, nh (nt)
                DO jh = ih, nh (nt)
                   !
                   ijh = ijh + 1
                   CALL pvan2_pwned (ngm, ih, jh, nt, qmod, qgm, ylmk0, prad_, &
                        SIZE(prad_,1),SIZE(prad_,2),SIZE(prad_,3),SIZE(prad_,4),lm)
                   DO na = 1, nat
                      !
                      IF (ityp(na).NE.nt) CYCLE
                      DO is = 1, nspin
!________________________________________________________________________
!                        lm: DO lm = 1, lmaxq**2
!                        ap__: IF ( ABS(ap(lm, nhtolm(ih,nt), nhtolm(jh,nt))) > eps8 ) THEN
!________________________________________________________________________

                         DO ig = 1, ngm
                            aux(ig,is,na) = aux(ig,is,na) +         &
                                            !ap(lm, nhtolm(ih,ityp(na)), nhtolm(jh,ityp(na)))*&
                                            qgm(ig)*becnew(ijh,na,is)
                         ENDDO
!________________________________________________________________________
!                      ENDIF ap__
!                      ENDDO lm
!________________________________________________________________________
                      ENDDO
                      !
                   ENDDO
                   !
                ENDDO
             ENDDO
          ENDIF
       ENDDO
       !
       !     convert aux to real space
       !
       DO na = 1, nat
          IF (tpawp(ityp(na))) THEN
             DO is = 1, nspin
                psic(:) = (0.d0, 0.d0)
                psic( nl(:) ) = aux(:,is,na)
                IF (gamma_only) psic( nlm(:) ) = CONJG(aux(:,is,na))
                CALL cft3 (psic, nr1, nr2, nr3, nrx1, nrx2, nrx3, 1)
                rho1_ (:, is, na) = DBLE (psic (:) )
             ENDDO
          END IF
       END DO
       !
    END DO whattodo
    !
    DEALLOCATE (ylmk0, qgm, qmod, aux)
    call stop_clock ('one-charge')

  END SUBROUTINE coc_pwned


  ! Analogous to PW/qvan2.f90
  SUBROUTINE pvan2_pwned (ngy, ih, jh, np, qmod, qg, ylmk0, prad_, s1, s2, s3, s4,lm)
    !
    !#include "f_defs.h"
    USE kinds, ONLY: DP
    USE us, ONLY: dq!, qrad
    USE uspp_param, ONLY: lmaxq, nbrx
    USE uspp, ONLY: nlx, lpl, lpx, ap, indv, nhtolm
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: s1,s2,s3,s4
    REAL(DP), INTENT(IN) :: prad_(s1,s2,s3,s4)
    INTEGER,INTENT(IN) :: lm

    INTEGER :: ngy, & ! input: the number of G vectors to compute
               ih,  & ! input: the first index of Q
               jh,  & ! input: the second index of Q
               np     ! input: the number of the pseudopotential

    REAL(DP) :: ylmk0 (ngy, lmaxq * lmaxq), & ! the spherical harmonics
                qmod (ngy)         ! input:  moduli of the q+g vectors
    COMPLEX(DP) :: qg (ngy)        ! output: the fourier transform of interest
    !
    !     here the local variables
    !
    COMPLEX(DP) :: sig ! (-i)^L

    INTEGER :: nb,          & ! the atomic index corresponding to ih
               mb,          & ! the atomic index corresponding to jh
               nmb,         & ! combined index (nb,mb)
               ivl,         & ! the lm corresponding to ih
               jvl,         & ! the lm corresponding to jh
               ig,          & ! counter on g vectors
               lp,          & ! the actual LM
               l,           & ! the angular momentum L
!               lm,          & ! the possible LM's compatible with ih,j
               i0, i1, i2, i3 ! counters for interpolation table

    REAL(DP) :: sixth,                & ! 1 divided by six
                dqi,                  & ! 1 divided dq
                qm,                   & ! qmod/dq
                px,                   & ! measures for interpolation table
                ux, vx, wx, uvx, pwx, & ! auxiliary variables for intepolation
                work                    ! auxiliary variable
    !
    LOGICAL :: new_qmod
    !
    ! compute the indices which correspond to ih,jh
    !
    call start_clock ('pvan2')
    sixth = 1.d0 / 6.d0
    dqi = 1 / dq
    nb = indv (ih, np)
    mb = indv (jh, np)
    IF (nb.GE.mb) THEN
       nmb = nb * (nb - 1) / 2 + mb
    ELSE
       nmb = mb * (mb - 1) / 2 + nb
    ENDIF
    ivl = nhtolm(ih, np)
    jvl = nhtolm(jh, np)
    IF (nb.GT.nbrx) CALL errore (' pvan2 ', ' nb.gt.nbrx ', nb)
    IF (mb.GT.nbrx) CALL errore (' pvan2 ', ' mb.gt.nbrx ', mb)
    IF (ivl.GT.nlx) CALL errore (' pvan2 ', ' ivl.gt.nlx  ', ivl)
    IF (jvl.GT.nlx) CALL errore (' pvan2 ', ' jvl.gt.nlx  ', jvl)
    qg(:) = (0.d0, 0.d0)
    !
    !    and make the sum over the non zero LM
    !
    !DO lm = 1, lpx (ivl, jvl)
       !lp = lpl (ivl, jvl, lm)
        lp=lm
       !
       ! extraction of angular momentum l from lp:
       !
       if (lp<1)  CALL errore (' qvan ', ' lp < 1 ', lp)
       l = sqrt(DBLE(lp-1)) + 1
       if (lp>49) CALL errore (' qvan ', ' lp > 49 ', lp)
       !
       sig = (0.d0, -1.d0) ** (l - 1)
       sig = sig * ap (lp, ivl, jvl)
       !
       new_qmod = .true.
       DO ig = 1, ngy
          !
          ! calculate quantites depending on the module of G only when needed
          !
          IF ( ig > 1 ) new_qmod = ABS( qmod(ig) - qmod(ig-1) ) > 1.0D-6
          IF ( new_qmod ) THEN
             qm = qmod (ig) * dqi
             px = qm - INT (qm)
             ux = 1.d0 - px
             vx = 2.d0 - px
             wx = 3.d0 - px
             i0 = INT( qm ) + 1
             i1 = i0 + 1
             i2 = i0 + 2
             i3 = i0 + 3
             uvx = ux * vx * sixth
             pwx = px * wx * 0.5d0
             work = prad_ (i0, nmb, l, np) * uvx * wx + &
                    prad_ (i1, nmb, l, np) * pwx * vx - &
                    prad_ (i2, nmb, l, np) * pwx * ux + &
                    prad_ (i3, nmb, l, np) * px * uvx
          ENDIF
          qg (ig) = qg (ig) + sig * ylmk0 (ig, lp) * work
       ENDDO
    !ENDDO
    call stop_clock ('pvan2')

    RETURN
  END SUBROUTINE pvan2_pwned

SUBROUTINE rad_dipole(rho1rad, rho1trad)
    USE kinds,                  ONLY : DP
    USE cell_base,              ONLY : at, alat
    USE constants,              ONLY : fpi, eps8
    USE parameters,             ONLY : npsx
    USE radial_grids,           ONLY : ndmx
    USE lsda_mod,               ONLY : nspin
    USE uspp_param,             ONLY : nhm, nh, lmaxq
    USE uspp,                   ONLY : ap
    USE ions_base,              ONLY : nat, ityp, ntyp => nsp
    USE atom,                   ONLY : r, rab, mesh, msh

    REAL(DP), TARGET, INTENT(IN) :: rho1rad(ndmx,lmaxq**2,nspin,nat) ! AE charge density on rad. grid
    REAL(DP), TARGET, INTENT(IN) :: rho1trad(ndmx,lmaxq**2,nspin,nat)! the same, but pseudo
    REAL(DP), POINTER            :: rho1rad_(:,:,:,:)                ! pointer to both

    INTEGER, PARAMETER           :: AE = 1, PS = 2   ! All-Electron and Pseudo
    INTEGER                      :: i_what, &
                                    na,nt, &     ! counter on atoms and atom types
                                    ispin, &     ! counter on spins
                                    lm,l,&       ! counter on composite angmom lm = l**2 +m
                                    k,& !DEBUG
                                    x,y,z

    REAL(DP)                     :: dipole(3,nspin,nat,2) !DEBUG
    REAL(DP)                     :: monopole(nspin,nat,2) !DEBUG
    INTEGER                      :: lm2c(9) = (/-1, 3, 1, 2, -1, -1, -1, -1, -1 /) ! lm to cartesian
    REAL(DP)                     :: aux(ndmx)
    REAL(DP)                     :: v(3),d,de(nspin,nat,2),me(nspin,nat,2)

    WRITE(6,*) "Compunting radial dipole..."
    
    WRITE(6,*) "=== === === === === === === === === ==="
    dipole(:,:,:,:) = 0._dp
    whattodo: DO i_what = AE, PS
    NULLIFY(rho1rad_)
    IF (i_what == AE) THEN
        rho1rad_ => rho1rad
    ELSE IF (i_what == PS) THEN
        rho1rad_ => rho1trad
    ENDIF
        atoms: DO na = 1, nat
        nt = ityp(na)
        spins: DO ispin = 1,nspin
            DO lm = 2, 4 ! l=1 m=z,y,x
!                 CALL simpson (msh(nt),rho1rad_(:,lm,ispin,na),rab(1,nt),monopole(ispin,na,i_what))
                k = lm2c(lm) ! k = 3,1,2 = z,x,y
                aux(:) = r(:,nt)*rho1rad_(:,lm,ispin,na)
                CALL simpson (msh(nt),aux,rab(1,nt),dipole(k,ispin,na,i_what))
            ENDDO ! lm
            WRITE(6,"(a,3i2,3f12.6)") " ===",ispin,na,i_what, dipole(:,ispin,na,i_what)
        ENDDO spins
        ENDDO atoms
    ENDDO whattodo

    de(:,:,:) = 0._dp
    me(:,:,:) = 0._dp
    dipole(1,1,1,1) = 0.1
    dipole(2,1,1,1) = 0.2
    dipole(3,1,1,1) = 0.4
    k = 128
    whattodo2: DO i_what = AE, PS
        atoms2: DO na = 1, nat
        nt = ityp(na)
        spins2: DO ispin = 1,nspin
            DO z = -k,k
            !WRITE(6,*) "==",z
            DO y = -k,k
            DO x = -k,k
                IF ((x/=0).or.(y/=0).or.(z/=0)) THEN
                ! v is the versor from center of cell [000] to center of cell [xyz]
                ! d is the distance between them
                v(:) = at(:,1)*x + at(:,2)*y + at(:,3)*z
                d = sqrt(SUM(v(:)*v(:)))
                v(:) = v(:) / d
                d = d*alat
!                 WRITE(6,"(a,3i2,4f15.6)")," ===",x,y,z,&
!                    (SUM(dipole(:,ispin,na,i_what)**2)-3._dp*SUM(dipole(:,ispin,na,i_what)*v(:))**2),d**3
!                 (SUM(dipole(:,ispin,na,i_what)**2)),(-3._dp*SUM(dipole(:,ispin,na,i_what)*v(:))**2),(d**3)
                de(ispin,na,i_what) = de(ispin,na,i_what)+&
                        6._dp *((SUM(dipole(:,ispin,na,i_what)**2)-3._dp*SUM(dipole(:,ispin,na,i_what)*v(:))**2)/(d**3))
!                 me(ispin,na,i_what) = me(ispin,na,i_what)+&
!                         monopole(ispin,na,i_what)**2/d
                ENDIF
            ENDDO
            ENDDO
            ENDDO
        WRITE(6,"(a,3i3,f15.8)") "dipole ==>", ispin,na,i_what,de(ispin,na,i_what)
        ENDDO spins2
        ENDDO atoms2
    ENDDO whattodo2

!     DO i_what = AE, PS
!     DO na = 1, nat
!     DO ispin = 1,nspin    
!         WRITE(6,"(a,3i3,f15.8)") "monopole ==>", ispin,na,i_what,me(ispin,na,i_what)
!     ENDDO
!     ENDDO
!     ENDDO
    WRITE(6,*) "=== === === === === === === === === ==="

END SUBROUTINE rad_dipole


