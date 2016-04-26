!#define __DEBUG_V_H_vs_SPHEROPOLE
!
! Copyright (C) 2004 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
MODULE atomic_paw
  !
  !============================================================================
  !
  !   Module for Projector Augmented Wave (PAW) calculations assuming
  !   spherical symmetry. Kresse's notations are adopted.
  !   Contains the type describing a PAW dataset, the routine for
  !   generating it from the ld1 code, and the routines to compute
  !   the hamiltonian and energy.
  !   GGA and LSD are implemented, relativistic calculations are not.
  !
  !   References:
  !   P.E.Blochl, Phys. Rev. B 50, 17953 (1994)
  !   G.Kresse, D.Joubert, Phys. Rev. B 59, 1758 (1999)
  !
  !   WARNINGS:
  !   Still work in progress on many aspects.
  !   The PAW dataset is written in a temporary format which is not supposed
  !   to be used any longer.
  !   Consistency with the input of the ld1 code yet to be checked
  !
  !   Written by Guido Fratesi, february 2005
  !   Modified by Riccardo Mazzarello, july 2006
  !   Fully Relativistic generalization by Andrea Dal Corso, november 2009 
  !   Other people involved: Lorenzo Paulatto and Stefano de Gironcoli
  !
  !============================================================================
  !
  USE kinds,            ONLY: dp
  USE ld1_parameters,   ONLY: nwfsx
  USE parameters,       ONLY: lmaxx
  USE constants,        ONLY: pi, fpi, e2, eps8
  USE radial_grids,     ONLY: ndmx, radial_grid_type
  USE paw_type,         ONLY: paw_t, nullify_pseudo_paw, allocate_pseudo_paw
  !
  IMPLICIT NONE
  PRIVATE
  SAVE
  !
  REAL(dp), PARAMETER :: ZERO=0._dp, ONE=1._dp, TWO=2._dp, HALF=0.5_dp
  !
  ! TEMP, to be substituted by module constants
!   REAL(dp), PARAMETER :: &
!        PI=3.14159265358979323846_dp, FPI=4._dp*PI, E2=2._dp, EPS8=1.0e-8_dp
  !
  !============================================================================
  !
  PUBLIC :: paw_t
  PUBLIC :: us2paw
  PUBLIC :: paw2us
  PUBLIC :: check_multipole
  PUBLIC :: new_paw_hamiltonian
  PUBLIC :: find_bes_qi
  PUBLIC :: compute_nonlocal_coeff_ion
  !
CONTAINS
  !
  !============================================================================
  !                          PUBLIC ROUTINES                                !!!
  !============================================================================
  !
  ! Compute the values of the local pseudopotential and the NL coefficients
  ! Compute also the total energy
  ! 
  SUBROUTINE new_paw_hamiltonian (veffps_, ddd_, etot_, &
       pawset_, nwfc_, l_, j_, nspin_, spin_, oc_, pswfc_, eig_, paw_energy,dddion_)
    IMPLICIT NONE
    REAL(dp), INTENT(OUT) :: veffps_(ndmx,2)
    REAL(dp), INTENT(OUT) :: ddd_(nwfsx,nwfsx,2)
    REAL(dp), INTENT(OUT) :: etot_
    TYPE(paw_t),   INTENT(IN)  :: pawset_
    INTEGER,       INTENT(IN)  :: nwfc_
    INTEGER,       INTENT(IN)  :: l_(nwfsx)
    INTEGER,       INTENT(IN)  :: nspin_
    INTEGER,       INTENT(IN)  :: spin_(nwfsx)
    REAL(dp), INTENT(IN)  :: j_(nwfsx)
    REAL(dp), INTENT(IN)  :: oc_(nwfsx)
    REAL(dp), INTENT(IN)  :: pswfc_(ndmx,nwfsx)
    REAL(dp), INTENT(IN)  :: eig_(nwfsx)
    REAL(dp), OPTIONAL :: dddion_(nwfsx,nwfsx)
    REAL(dp), INTENT(OUT), OPTIONAL :: paw_energy(5,3) 
    !
    REAL(dp) :: &                                        ! one center:
         eps,             e1,             e1ps,             & ! energies;
                          veff1(ndmx,2),   veff1ps(ndmx,2),   & ! eff potentials;
         chargeps(ndmx,2), charge1(ndmx,2), charge1ps(ndmx,2), & ! charges.
         projsum(nwfsx,nwfsx,2), eigsum !  sum of projections, sum of eigenval.
    !
    INTEGER :: ns, is, n
    REAL(dp) :: energy(5,3)
    !
    ! Compute the valence charges
    CALL compute_charges(projsum, chargeps, charge1, charge1ps, &
       pawset_, nwfc_, l_, j_, nspin_, spin_, oc_, pswfc_, 1 )
 !
 !  Check for negative charge
 !
     do is=1,nspin_
        do n=2,pawset_%grid%mesh
!           write(6,*) n, pawset_%grid%r(n), chargeps(n,is)
           if (chargeps(n,is)<-1.d-12) &
                   call  errore('new_paw_hamiltonian','negative rho',1)
        enddo
     enddo

!     write(766,"(4f12.6)") (pawset_%grid%r(ns), chargeps(ns,1), charge1(ns,1), charge1ps(ns,1), ns=1,pawset_%grid%mesh)
!     write(767,"(4f12.6)") (pawset_%grid%r(ns), chargeps(ns,2), charge1(ns,2), charge1ps(ns,2), ns=1,pawset_%grid%mesh)
    !
    ! Compute the one-center energy and effective potential:
    ! E = Eh[n_v] + Exc[n_v+n_c] - Int[veff*n_v],
    ! veff = vh[n_v] + v_xc[n_v+n_c],
    ! where n_v can be n~+n^ or n1 or n1~+n^ and n_c can be nc or nc~
    ! n~+n^ , nc~
    CALL compute_onecenter_energy ( eps,  veffps_, &
       pawset_, chargeps,  pawset_%nlcc, pawset_%psccharge, nspin_,&
             pawset_%grid%mesh, pawset_%psloc, energy(1,1))
    ! n1 , nc
    CALL compute_onecenter_energy ( e1,   veff1, &
       pawset_, charge1,  .TRUE.,        pawset_%aeccharge, nspin_,&
            pawset_%irc, pawset_%aeloc, energy(1,2))
    ! n1~+n^ , nc~
    CALL compute_onecenter_energy ( e1ps, veff1ps, &
       pawset_, charge1ps, pawset_%nlcc, pawset_%psccharge, nspin_,&
            pawset_%irc, pawset_%psloc, energy(1,3))
    ! Add the local part
    DO is=1,nspin_
       veffps_(1:pawset_%grid%mesh,is) = veffps_(1:pawset_%grid%mesh,is) +  &
            pawset_%psloc(1:pawset_%grid%mesh)
       veff1  (1:pawset_%grid%mesh,is) = veff1  (1:pawset_%grid%mesh,is) +  &
            pawset_%aeloc(1:pawset_%grid%mesh)
       veff1ps(1:pawset_%grid%mesh,is) = veff1ps(1:pawset_%grid%mesh,is) +  &
            pawset_%psloc(1:pawset_%grid%mesh)
    END DO
    !
    ! Compute the nonlocal D coefficients
    CALL compute_nonlocal_coeff (ddd_,pawset_,nspin_,veffps_,veff1,veff1ps)
    IF (PRESENT(dddion_)) THEN
       CALL compute_nonlocal_coeff_ion (dddion_,pawset_)
    END IF
    !
    ! Compute total energy
    eigsum=ZERO
    DO ns=1,nwfc_
       IF (oc_(ns)>ZERO) eigsum = eigsum + oc_(ns)*eig_(ns)
    END DO
    etot_ = eigsum + eps + e1 - e1ps

    if (PRESENT(paw_energy)) paw_energy=energy
    !
  END SUBROUTINE new_paw_hamiltonian
  !
  !============================================================================
  !
  ! Convert the USPP to a PAW dataset
  !
  SUBROUTINE us2paw (pawset_,                                        &
       zval, grid, rmatch_augfun, ikk,                               &
       nbeta, lls, jjs, ocs, enls, els, rcutus, psipaw, psipaw_rel,  &
       phis, betas,  qvan, kindiff,                                  &
       nlcc, aerhoc, psrhoc, aevtot, psvtot, which_paw_augfun,rel     )

    USE funct,        ONLY : get_dft_long
    USE ld1inc,       ONLY : zed, file_screen
    USE paw_type,     ONLY : nullify_pseudo_paw, allocate_pseudo_paw
    USE io_global,    ONLY : stdout, ionode, ionode_id
    USE radial_grids, ONLY : allocate_radial_grid
    USE mp,           only : mp_bcast
    USE mp_world,     only : world_comm
    IMPLICIT NONE
    TYPE(paw_t),      INTENT(OUT) :: pawset_
    REAL(dp), INTENT(IN)  :: zval
    type(radial_grid_type), INTENT(IN) :: grid
    REAL(dp), INTENT(IN)  :: rmatch_augfun
    INTEGER,  INTENT(IN)  :: ikk(nwfsx)
    !
    INTEGER,  INTENT(IN)  :: nbeta, rel
    INTEGER,  INTENT(IN)  :: lls(nwfsx)
    REAL(dp), INTENT(IN)  :: ocs(nwfsx)
    CHARACTER(LEN=2), INTENT(IN)  :: els(nwfsx)
    REAL(dp), INTENT(IN)  :: jjs(nwfsx)
    REAL(dp), INTENT(IN)  :: rcutus(nwfsx)
    REAL(dp), INTENT(IN)  :: enls(nwfsx)
    REAL(dp), INTENT(IN)  :: psipaw(ndmx,nwfsx)
    REAL(dp), INTENT(IN)  :: psipaw_rel(ndmx,nwfsx)
    REAL(dp), INTENT(IN)  :: phis(ndmx,nwfsx)
    REAL(dp), INTENT(IN)  :: betas(ndmx,nwfsx)
    !
    REAL(dp), INTENT(IN)  :: qvan(ndmx,nwfsx,nwfsx)
    REAL(dp), INTENT(IN)  :: kindiff(nwfsx,nwfsx)
    LOGICAL,  INTENT(IN)  :: nlcc
    REAL(dp), INTENT(IN)  :: aerhoc(ndmx)
    REAL(dp), INTENT(IN)  :: psrhoc(ndmx)
    REAL(dp), INTENT(IN)  :: aevtot(ndmx)
    REAL(dp), INTENT(IN)  :: psvtot(ndmx)
    CHARACTER(20), INTENT(IN)  :: which_paw_augfun
    !
    REAL(DP),  EXTERNAL :: int_0_inf_dr
    CHARACTER(len=2), EXTERNAL :: atom_name
    REAL(dp) :: vps(ndmx,2), projsum(nwfsx,nwfsx,2), ddd(nwfsx,nwfsx,2), dddion(nwfsx,nwfsx)
    INTEGER  :: irc, ns, ns1, n, leading_power, mesh, ios
    REAL(dp) :: aux(ndmx), aux2(ndmx,2), raux
    REAL(dp) :: aecharge(ndmx,2), pscharge(ndmx,2)
    REAL(dp) :: etot
    INTEGER  :: nspin=1, spin(nwfsx)=1 ! PAW generat. from spin-less calculation
    !
    ! variables for aug. functions generation
    ! 
    REAL(DP) :: energy(5,3), max_aug_cutoff = -1._dp
    INTEGER  :: nc, iok          ! index Bessel function, ...
    INTEGER :: l1,l2,l3, lll, ircm, ir, ir0
    REAL(dp) :: twosigma2, rm                  ! needed for "gaussian" augfun
    REAL(dp) :: zeta, resi,usigma=4._dp        ! needed for "Bloechl" augfun
    REAL(DP),ALLOCATABLE :: gaussian(:)        ! needed for gaussian and Bloechl
    REAL(DP),ALLOCATABLE :: aug_real(:,:,:)    ! needed for PSQ augfun
    REAL(DP) :: qc(2), xc(2), b1(2), b2(2)     ! needed for BESSEL augfun
    REAL(DP), ALLOCATABLE :: j1(:,:)           ! needed for BESSEL augfun
    !
    mesh = grid%mesh
    irc = maxval(ikk(1:nbeta))+8
    CALL nullify_pseudo_paw(pawset_)
    CALL allocate_pseudo_paw(pawset_,ndmx,nwfsx,lmaxx)
    CALL allocate_radial_grid(pawset_%grid, mesh)
    pawset_%symbol = atom_name(nint(zed))
    pawset_%zval = zval
    !
    ! Copy the mesh
    pawset_%grid%mesh  = grid%mesh
    pawset_%grid%r(:)  = grid%r(:)
    pawset_%grid%r2(:) = grid%r2(:)
    pawset_%grid%rab(:)= grid%rab(:)
    pawset_%grid%sqr(:)= grid%sqr(:)
    pawset_%grid%xmin  = grid%xmin
    pawset_%grid%rmax  = grid%rmax
    pawset_%grid%zmesh = grid%zmesh
    pawset_%grid%dx    = grid%dx
    !
    pawset_%rmatch_augfun = rmatch_augfun
    !if (rmatch_augfun <= 0.0_dp) pawset_%rmatch_augfun = grid%r(irc)
    if (rmatch_augfun <= 0.0_dp) pawset_%rmatch_augfun = MAXVAL(rcutus(1:nbeta))
    pawset_%rel = rel
    pawset_%irc = irc
    pawset_%ikk(1:nbeta)=ikk(1:nbeta)
    !
    ! Copy the wavefunctions
    pawset_%nwfc = nbeta
    pawset_%l(1:nbeta) = lls(1:nbeta)
    pawset_%lmax = MAXVAL(pawset_%l(1:nbeta))
    pawset_%oc(1:nbeta) = ocs(1:nbeta)
    pawset_%jj(1:nbeta) = jjs(1:nbeta)
    pawset_%rcutus(1:nbeta) = rcutus(1:nbeta)
    pawset_%els(1:nbeta)= els(1:nbeta)
    pawset_%enl(1:nbeta)= enls(1:nbeta)
    pawset_%aewfc(1:mesh,1:nbeta) = psipaw(1:mesh,1:nbeta)
    pawset_%aewfc_rel(1:mesh,1:nbeta) = psipaw_rel(1:mesh,1:nbeta)
    pawset_%pswfc(1:mesh,1:nbeta) = phis  (1:mesh,1:nbeta)
    pawset_%proj (1:mesh,1:nbeta) = betas (1:mesh,1:nbeta)
    !
    pawset_%augfun = 0._dp
    !
    !
    ! Augmentation functions:
    ! The specific radial part is not important, as long as the
    ! multipole moments of the exact augmentation functions are
    ! reproduced. So we write on the file the exact augmentation
    ! functions and their multipole moments, keeping in mind that
    ! the PW program should not use the radial part as is but
    ! substitute it with some smoothened analogue.
    !
    ! Sdg>> not quite sure this is correct because the shape of augfun,
    ! arbitrary as it may be, determines the shape of vloc since this
    ! is obtained by unscreening vscf with the GIVEN augfun...
    !
    ! moreover in order to give the right electrostatics in the solid the
    ! augmentation functions should not overlap.
    !
    ! Compute the exact augmentation functions and their moments
    !
    write(stdout,'(//,4x,a)') 'multipoles (all-electron charge) - (pseudo charge)'
    write(stdout,'(7x,2a,":",2a,2x,6a)') ' ns',' l1','ns1',' l2',&
                 '  l=0   ','  l=1   ','  l=2   ','  l=3   ','  l=4   ', '  l=5  '
    pawset_%augfun(:,:,:,:) = 0.0_dp
    pawset_%augmom(:,:,:) = 0.0_dp
    DO ns=1,nbeta
       l1=pawset_%l(ns)
       DO ns1=1,ns
          l2=pawset_%l(ns1)
          do l3 = max(l1-l2,l2-l1), l1+l2, 2
             pawset_%augfun(1:mesh,ns,ns1,l3) = &
                 pawset_%aewfc(1:mesh,ns) * pawset_%aewfc(1:mesh,ns1) - &
                 pawset_%pswfc(1:mesh,ns) * pawset_%pswfc(1:mesh,ns1)
             IF (pawset_%rel==2) &
             pawset_%augfun(1:irc,ns,ns1,l3) =pawset_%augfun(1:irc,ns,ns1,l3)&
                 +pawset_%aewfc_rel(1:irc,ns) * pawset_%aewfc_rel(1:irc,ns1)
             pawset_%augfun(1:mesh,ns1,ns,l3) = pawset_%augfun(1:mesh,ns,ns1,l3)
             aux(1:irc) = pawset_%augfun(1:irc,ns,ns1,l3) * pawset_%grid%r(1:irc)**l3
             lll = l1 + l2 + 2 + l3
             pawset_%augmom(ns,ns1,l3)=int_0_inf_dr(aux(1:irc),pawset_%grid,irc,lll)
             pawset_%augmom(ns1,ns,l3)=pawset_%augmom(ns,ns1,l3)
          end do
          write (stdout,'(7x,2i3,a,2i3,10f8.4)') ns,l1,":",ns1,l2,&
                              (pawset_%augmom(ns,ns1,l3), l3=0,l1+l2)
       END DO
    END DO
    !
    !
    IF (which_paw_augfun/='AE') THEN
       ! The following lines enables to test the idependence on the 
       ! specific shape of the radial part of the augmentation functions
       ! in a spherically averager system (as is the case in atoms) only
       ! the zero-th moment contribute to the scf charge
       write(stdout, '(5x,3a,f12.6)') 'Required augmentation: ',TRIM(which_paw_augfun), "   radius:", pawset_%rmatch_augfun
       !
101    pawset_%augfun(:,:,:,:) = 0.0_dp
       DO ns=1,nbeta
          l1 = pawset_%l(ns)
          DO ns1=1,ns
             l2 = pawset_%l(ns1)
             !
             SELECT CASE (which_paw_augfun)
             CASE ('PSQ')
                continue ! the work is done at the end
             CASE ('QVAN')
                IF(ns==1 .and. ns1==1) &
                CALL infomsg('us2paw', 'WARNING: QVAN augmentation function are for testing ONLY: '//&
                                       'they will not work in pw!')
                ! Choose the shape of the augmentation functions: NC Q ...
                pawset_%augfun(1:mesh,ns,ns1,0) = qvan(1:mesh,ns,ns1)
                ! Renormalize the aug. functions so that their integral is correct
                leading_power = l1 + l2 + 2 
                raux=int_0_inf_dr(pawset_%augfun(1:irc,ns,ns1,0),pawset_%grid,irc,leading_power)
                IF (ABS(raux) < eps8) CALL errore &
                   ('ld1_to_paw','norm of augm.func. is too small',ns*100+ns1)
                raux=pawset_%augmom(ns,ns1,0)/raux
                pawset_%augfun(1:mesh,ns,ns1,0)=raux*pawset_%augfun(1:mesh,ns,ns1,0)
                !
             CASE ('BG')
                IF(ns==1 .and. ns1==1) &
                CALL infomsg('us2paw', 'WARNING: using Bloechl style augmentation functions '//&
                                       'is not a good idea, as analytical overlap are not '//&
                                       'implemented in pwscf: use BESSEL or GAUSS instead.')
                IF(.not. allocated(gaussian)) ALLOCATE(gaussian(mesh))
                ! use Bloechl style augmentation functions, as linear combinations of
                ! gaussians functions (this is quite pointless if the the plane-wave
                ! code doesn't use double-augmentation with analytical gaussian overlaps)
                DO ir0 = 1,mesh
                    IF(grid%r(ir0) > pawset_%rmatch_augfun) &
                        exit
                ENDDO
                ! look for a sigma large enough to keep (almost) all the gaussian in the aug sphere
                zeta = (usigma/pawset_%rmatch_augfun)**2
                DO ir = 1,mesh
                    gaussian(ir) = exp( - zeta*(grid%r(ir))**2 ) * grid%r2(ir)
                ENDDO
                DO l3 = max (l1-l2,l2-l1), l1+l2
                    ! Functions has to be normalized later, so I can use a constant factor 
                    ! = rc**l3 to prevent very large numbers when integrating:
                    aux(1:grid%mesh) = gaussian(1:grid%mesh) * grid%r(1:grid%mesh)**l3
                    ! Normalization to have unitary multipole l3
                    ! and check norm of augmentation functions.
                    raux = int_0_inf_dr(aux,pawset_%grid,mesh,l3+2)
                    IF (abs(raux) < eps8) CALL errore &
                        ('ld1_to_paw','norm of augm. func. is too small',ns*100+ns1)
                    ! Check if gaussians are localized enough into the augmentation sphere,
                    ! otherwise try again with a larger sigma.
                    resi = int_0_inf_dr(aux,pawset_%grid,ir0, l3+2)
                    resi = (raux-resi)/raux
                    IF (abs(resi) > 1.e-10_dp) THEN
                        usigma = usigma + .0625_dp
                        GOTO 101
                    ENDIF
                    raux=pawset_%augmom(ns,ns1,l3)/raux

                    pawset_%augfun(1:mesh,ns,ns1,l3) = raux*gaussian(1:mesh)
                    pawset_%augfun(1:mesh,ns1,ns,l3) = raux*gaussian(1:mesh)
                ENDDO
                DEALLOCATE(gaussian)
                !
             CASE ('GAUSS')
                ! use linear combinations of gaussians functions, not the Bloechl style
                ! but it looks a bit alike... (used for testing, now obsolete)
                IF(ns==1 .and. ns1==1) &
                CALL infomsg('us2paw', 'GAUSS augmentation functions are ususally '//&
                                       'harder than BESSEL; use BESSEL instead unless'//&
                                       ' you have discontinuity in local potential')
                ALLOCATE(gaussian(mesh))
                !
                rm = pawset_%rmatch_augfun
                twosigma2 = TWO*(rm/3.0_dp)**2
                do ir=1,mesh
                   if (grid%r(ir) <= rm) then
                      gaussian(ir) = ( EXP(-grid%r(ir)**2/twosigma2) + &
                                       EXP(-(grid%r(ir)-TWO*rm)**2/twosigma2 ) - &
                                       TWO*EXP(-rm**2/twosigma2) ) * grid%r2(ir)
                   else
                      gaussian(ir) = 0.0_dp
                   endif
                end do
                DO l3 = max (l1-l2,l2-l1), l1+l2
                   ! 
                   aux(1:irc) = gaussian(1:irc) * pawset_%grid%r(1:irc)**l3
                   ! calculate the corresponding multipole
                   raux=int_0_inf_dr(aux,pawset_%grid,irc,l3+2)
                   IF (ABS(raux) < eps8) CALL errore &
                      ('ld1_to_paw','norm of augm. func. is too small',ns*100+ns1)
                   raux=pawset_%augmom(ns,ns1,l3)/raux
                   pawset_%augfun(1:mesh,ns,ns1,l3) = raux*gaussian(1:mesh)
                   pawset_%augfun(1:mesh,ns1,ns,l3) = raux*gaussian(1:mesh)
                   ! 
                END DO
                DEALLOCATE(gaussian)
                !
             CASE ('BESSEL')
                ! Use Kresse-Joubert style augmentation functions
                ! Defined as linear combination of Bessel functions.
                ALLOCATE (j1 (pawset_%grid%mesh,2))
                do ir=1,irc
                   !if (grid%r(ir)<pawset_%rmatch_augfun) ircm=ir
                   if (grid%r(ir)<pawset_%rmatch_augfun) ircm=ir
                end do
                do l3 = max(l1-l2,l2-l1), l1+l2 
                   ! 
                   CALL find_bes_qi(qc,pawset_%grid%r(ircm),l3,2,iok)
                   IF (iok.ne.0) CALL errore('compute_augfun', &
                         'problems with find_bes_qi',1)
                   DO nc = 1, 2
                      !
                      CALL sph_bes(irc,grid%r(1),qc(nc),l3,j1(1,nc))
                      aux(1:irc) = j1(1:irc,nc) * grid%r(1:irc)**(l3+2)
                      b1(nc) = j1(ircm,nc)
                      b2(nc) = int_0_inf_dr(aux,pawset_%grid,ircm,l3+2)
                      !
                   ENDDO
                   xc(1) = b1(2) / (b1(2) * b2(1) - b1(1) * b2(2))
                   xc(2) = - b1(1) * xc(1) / b1(2)
                   pawset_%augfun(1:ircm,ns,ns1,l3) =                   &
                          pawset_%augmom(ns,ns1,l3) * grid%r2(1:ircm) * &
                          (xc(1) * j1(1:ircm,1) + xc(2) * j1(1:ircm,2))
                   ! Symmetrize augmentation functions:
                   pawset_%augfun(1:mesh,ns1,ns,l3)=pawset_%augfun(1:mesh,ns,ns1,l3)
                   !
                   ! Save higher bessel coefficient to compute suggested cutoff
                   max_aug_cutoff=MAX( max_aug_cutoff, MAXVAL(qc(1:2))**2)
                   !
                end do 
                DEALLOCATE (j1)
                !
             CASE DEFAULT
                !
                CALL errore ('ld1_to_paw','Specified augmentation functions ('// &
                             TRIM(which_paw_augfun)//') not allowed or coded',1)
                !
             END SELECT
          END DO
       END DO
       IF ( which_paw_augfun == 'BG') &
         write(stdout,"(5x,a,f12.6)") "Gaussians generated with zeta: ", zeta
       IF ( which_paw_augfun == 'BESSEL') &
         write(stdout,'(5x, "Suggested rho cutoff for augmentation:",f7.2," Ry")') max_aug_cutoff
       !
       IF ( which_paw_augfun == 'PSQ') THEN
            ALLOCATE(aug_real(ndmx,nwfsx,nwfsx))
            DO ns=1,nbeta
            DO ns1=ns,nbeta
               aug_real(1:mesh,ns,ns1) = &
                     pawset_%aewfc(1:mesh,ns) * pawset_%aewfc(1:mesh,ns1) - &
                     pawset_%pswfc(1:mesh,ns) * pawset_%pswfc(1:mesh,ns1)
               IF (pawset_%rel==2) &
               aug_real(1:irc,ns,ns1) = aug_real(1:irc,ns,ns1) + &
                     pawset_%aewfc_rel(1:irc,ns)*pawset_%aewfc_rel(1:irc,ns1) 
               aug_real(1:mesh,ns1,ns) = aug_real(1:mesh,ns,ns1)
            ENDDO
            ENDDO
            !
            CALL pseudo_q(aug_real,pawset_%augfun)
            !
            DEALLOCATE(aug_real)
       ENDIF

    END IF
!
!   Outside the PAW spheres augfun should be exactly 0. On some machine
!   it is equal to zero to machine precision and sometimes it is negative, 
!   so as to confuse the check for negative charge. So we set it to zero
!   explicitely.
!
    DO ns = 1, nbeta
       l1 = pawset_%l(ns)
       DO ns1 = ns, nbeta
          l2 = pawset_%l(ns1)
          DO l3 = max (l1-l2,l2-l1), l1+l2
             DO n = pawset_%irc+1, pawset_%grid%mesh
                pawset_%augfun(n,ns,ns1,l3) = 0.0_DP
                pawset_%augfun(n,ns1,ns,l3) = 0.0_DP
             END DO
          END DO
       END DO
    END DO
    !
    !
    ! Copy kinetic energy differences
    pawset_%kdiff(1:nbeta,1:nbeta)=kindiff(1:nbeta,1:nbeta)
    !
    ! Copy the core charge (if not NLCC copy only the AE one)
    pawset_%nlcc=nlcc
    IF (pawset_%nlcc) pawset_%psccharge(1:mesh)=psrhoc(1:mesh)
    pawset_%aeccharge(1:mesh)=aerhoc(1:mesh)
    !
    ! Copy the local potentials
    pawset_%psloc(1:mesh)=psvtot(1:mesh)
    pawset_%aeloc(1:mesh)=aevtot(1:mesh)
    ! and descreen them:
    CALL compute_charges(projsum, pscharge, aecharge, aux2, &
       pawset_, nbeta, lls, jjs, nspin, spin, ocs, phis )
    pawset_%pscharge(1:mesh)=pscharge(1:mesh,1)
    !
    CALL compute_onecenter_energy ( raux,  aux2, &
       pawset_, pscharge, pawset_%nlcc, pawset_%psccharge, nspin, &
       pawset_%grid%mesh )
    pawset_%psloc(1:mesh)=psvtot(1:mesh)-aux2(1:mesh,1)
    !
    CALL compute_onecenter_energy ( raux,  aux2, &
       pawset_, aecharge, .TRUE.,       pawset_%aeccharge, nspin, &
       pawset_%grid%mesh)
    pawset_%aeloc(1:mesh)=aevtot(1:mesh)-aux2(1:mesh,1)
    !
    if (file_screen .ne.' ') then
        if (ionode) &
            open(unit=20,file=file_screen, status='unknown', iostat=ios, err=105 )
105     call mp_bcast(ios, ionode_id, world_comm)
        call errore('descreening','opening file'//file_screen,abs(ios))
        if (ionode) then
            write(20,'(a)') "#   n, r(n),       aeloc(n),   psloc(n),   pscharge(n)"
            do n=1,grid%mesh
            write(20,'(i5,7e12.4)') n,grid%r(n), pawset_%aeloc(n), &
                    pawset_%psloc(n), pawset_%pscharge(n)
            enddo
            close(20)
        endif
    endif
    !
    write(pawset_%dft,'(80x)') !fill it with spaces
    pawset_%dft = get_dft_long ( )
    !
    ! Generate the paw hamiltonian for test (should be equal to the US one)
    CALL new_paw_hamiltonian (vps, ddd, etot, &
       pawset_, pawset_%nwfc, pawset_%l, pawset_%jj, nspin, spin, pawset_%oc, &
       pawset_%pswfc, pawset_%enl, energy, dddion)
    pawset_%dion(1:nbeta,1:nbeta)=dddion(1:nbeta,1:nbeta)
    WRITE(stdout,'(/5x,A,f12.6,A)') 'Estimated PAW energy =',etot,' Ryd'
    WRITE(stdout,'(/5x,A)') 'The PAW screened D coefficients'
    DO ns1=1,pawset_%nwfc
       WRITE(stdout,'(6f12.5)') (ddd(ns1,ns,1),ns=1,pawset_%nwfc)
    END DO
    WRITE(stdout,'(/5x,A)') 'The PAW descreened D coefficients (US)'
    DO ns1=1,pawset_%nwfc
       WRITE(stdout,'(6f12.5)') (dddion(ns1,ns),ns=1,pawset_%nwfc)
    END DO
    !
    !
  END SUBROUTINE us2paw
  !
  !============================================================================
  !
  ! ...
  !
  SUBROUTINE paw2us (pawset_,zval,grid,nbeta,lls,jjs,ikk,betas,qq,qvan,&
                     vpsloc, bmat, rhos, els, rcutus, pseudotype,psipaw_rel)
    USE funct, ONLY : set_dft_from_name
    use radial_grids, only: radial_grid_type, allocate_radial_grid
    IMPLICIT NONE
    TYPE(radial_grid_type), INTENT(OUT) :: grid
    TYPE(paw_t),   INTENT(IN)  :: pawset_
    REAL(dp), INTENT(OUT) :: zval
    INTEGER,       INTENT(OUT) :: nbeta
    INTEGER,       INTENT(OUT) :: lls(nwfsx)
    INTEGER,       INTENT(OUT) :: ikk(nwfsx)
    INTEGER,       INTENT(OUT) :: pseudotype
    CHARACTER(LEN=2), INTENT(OUT) :: els(nwfsx)
    REAL(dp), INTENT(OUT) :: jjs(nwfsx)
    REAL(dp), INTENT(OUT) :: rcutus(nwfsx)
    REAL(dp), INTENT(OUT) :: betas(ndmx,nwfsx)
    REAL(dp), INTENT(OUT) :: psipaw_rel(ndmx,nwfsx)
    REAL(dp), INTENT(OUT) :: qq(nwfsx,nwfsx)
    REAL(dp), INTENT(OUT) :: qvan(ndmx,nwfsx,nwfsx)
    REAL(dp), INTENT(OUT) :: vpsloc(ndmx)     ! the local pseudopotential
    REAL(dp), INTENT(OUT) :: bmat(nwfsx,nwfsx)! the pseudo coefficients (unscreened D)
    REAL(dp), INTENT(OUT) :: rhos(ndmx)     ! the pseudo density

    INTEGER :: ns, ns1, mesh
    zval=pawset_%zval
    !
    mesh=pawset_%grid%mesh
    ! Copy the mesh
    call allocate_radial_grid(grid, mesh)
    grid%mesh  = pawset_%grid%mesh
    grid%r(1:mesh)  = pawset_%grid%r(1:mesh)
    grid%r2(1:mesh) = pawset_%grid%r2(1:mesh)
    grid%rab(1:mesh)= pawset_%grid%rab(1:mesh)
    grid%sqr(1:mesh)= pawset_%grid%sqr(1:mesh)
    grid%xmin  = pawset_%grid%xmin
    grid%rmax  = pawset_%grid%rmax
    grid%zmesh = pawset_%grid%zmesh
    grid%dx    = pawset_%grid%dx

    !
    nbeta=pawset_%nwfc
    lls(1:nbeta)=pawset_%l(1:nbeta)
    jjs(1:nbeta)=pawset_%jj(1:nbeta)
    els(1:nbeta)=pawset_%els(1:nbeta)
    rcutus(1:nbeta)=pawset_%rcutus(1:nbeta)
    ikk(1:nbeta)=pawset_%ikk(1:nbeta)
    !
    DO ns=1,nbeta
       DO ns1=1,nbeta
          qvan(1:mesh,ns,ns1)=pawset_%augfun(1:mesh,ns,ns1,0)
          IF (lls(ns)==lls(ns1)) THEN
             qq(ns,ns1)=pawset_%augmom(ns,ns1,0)
          ELSE ! different spherical harmonic => 0
             qq(ns,ns1)=ZERO
          END IF
       END DO
    END DO
    !
    vpsloc(1:mesh)=pawset_%psloc(1:mesh)
    bmat(1:nbeta,1:nbeta)=pawset_%dion(1:nbeta,1:nbeta)
    !
    rhos(1:mesh)=pawset_%pscharge(1:mesh)
    !
    psipaw_rel(1:mesh,1:nbeta)=pawset_%aewfc_rel(1:mesh,1:nbeta)
    betas(1:mesh,1:nbeta)=pawset_%proj(1:mesh,1:nbeta)
    pseudotype=3
    !
    CALL set_dft_from_name (pawset_%dft)

  END SUBROUTINE paw2us
  !
  !============================================================================
  !
  ! ...
  !
  SUBROUTINE check_multipole (pawset_)
    USE radial_grids, ONLY: hartree
    USE io_global, ONLY : stdout
    IMPLICIT NONE
    TYPE(paw_t), INTENT(IN)  :: pawset_
    INTEGER:: mesh
    REAL(dp) :: r(ndmx), r2(ndmx), sqr(ndmx), dx
    INTEGER :: nbeta
    INTEGER :: lls(nwfsx)
    INTEGER :: ir, ns1, ns2, l1, l2, l3, irc, ir0
    REAL(dp) :: auxpot(ndmx,0:2*lmaxx+2), auxrho(ndmx)
    !
    ! set a few internal variables 
    write (stdout,*) "check_multipole : lmaxx =",lmaxx
    mesh=pawset_%grid%mesh
    r(1:mesh)=pawset_%grid%r(1:mesh)
    r2(1:mesh)=pawset_%grid%r2(1:mesh)
    sqr(1:mesh)=pawset_%grid%sqr(1:mesh)
    dx=pawset_%grid%dx
    irc = pawset_%irc
    !
    nbeta=pawset_%nwfc
    lls(1:nbeta)=pawset_%l(1:nbeta)
    !
    do ns1=1,nbeta
       l1 = lls(ns1)
       do ns2=1,nbeta
          l2 = lls(ns2)
          auxpot(:,:) = 0.0_dp
          do l3 = max(l1-l2,l2-l1), l1+l2
             auxrho(1:mesh) = &
                pawset_%aewfc(1:mesh,ns1) * pawset_%aewfc(1:mesh,ns2) - &
                pawset_%pswfc(1:mesh,ns1) * pawset_%pswfc(1:mesh,ns2) - &
                pawset_%augfun(1:mesh,ns1,ns2,l3) 
             call hartree(l3,l1+l2+2,mesh,pawset_%grid,auxrho,auxpot(1,l3))
          end do
          write (stdout,'(a,2i3,a,2i3)') " MULTIPOLO DI ",ns1,l1,":",ns2, l2
          do ir=1,irc
             if (r(ir) < 1.0_dp) ir0 = ir
          end do
          do ir=ir0,irc+30, 3
             write (stdout,'(10f8.4)') r(ir),(auxpot(ir,l3), l3=0,l1+l2)
          end do
       end do
    end do
    return
  END SUBROUTINE check_multipole
  !
  !============================================================================
  !                          PRIVATE ROUTINES                               !!!
  !============================================================================
  !
  SUBROUTINE compute_charges (projsum_, chargeps_, charge1_, charge1ps_, &
       pawset_, nwfc_, l_, j_, nspin_, spin_, oc_, pswfc_ , iflag, unit_)
    USE io_global, ONLY : ionode
    IMPLICIT NONE
    REAL(dp), INTENT(OUT) :: projsum_(nwfsx,nwfsx,2)
    REAL(dp), INTENT(OUT) :: chargeps_(ndmx,2)
    REAL(dp), INTENT(OUT) :: charge1_(ndmx,2)
    REAL(dp), INTENT(OUT) :: charge1ps_(ndmx,2)
    TYPE(paw_t),   INTENT(IN)  :: pawset_
    INTEGER,       INTENT(IN)  :: nwfc_
    INTEGER,       INTENT(IN)  :: l_(nwfsx)
    INTEGER,       INTENT(IN)  :: nspin_
    INTEGER,       INTENT(IN)  :: spin_(nwfsx)
    REAL(dp), INTENT(IN)  :: j_(nwfsx)
    REAL(dp), INTENT(IN)  :: oc_(nwfsx)
    REAL(dp), INTENT(IN)  :: pswfc_(ndmx,nwfsx)
    INTEGER, OPTIONAL :: unit_, iflag
    REAL(dp) :: augcharge(ndmx,2), chargetot
    INTEGER :: i, n, iflag0

    iflag0=0
    if (present(iflag)) iflag0=iflag
    CALL compute_sumwfc2(chargeps_,pawset_,nwfc_,pswfc_,oc_,spin_)
    CALL compute_projsum(projsum_,pawset_,nwfc_,l_,j_,spin_,pswfc_,oc_)
!    WRITE (6200,'(20e20.10)') ((projsum_(ns,ns1,1),ns=1,ns1),ns1=1,pawset_%nwfc)
    CALL compute_onecenter_charge(charge1ps_,pawset_,projsum_,nspin_,"PS")
    CALL compute_onecenter_charge(charge1_  ,pawset_,projsum_,nspin_,"AE")
    ! add augmentation charges
    CALL compute_augcharge(augcharge,pawset_,projsum_,nspin_)

    if (present(unit_).and.ionode) then
       write(unit_,*) 
       write(unit_,*) "#"
       do i=1,pawset_%grid%mesh
          write (unit_,'(4f12.8)') pawset_%grid%r(i), augcharge(i,1), chargeps_(i,1), charge1ps_(i,1)
       end do
    end if
    chargeps_ (1:pawset_%grid%mesh,1:nspin_) = chargeps_ (1:pawset_%grid%mesh,1:nspin_) &
                                        + augcharge(1:pawset_%grid%mesh,1:nspin_)
    charge1ps_(1:pawset_%grid%mesh,1:nspin_) = charge1ps_(1:pawset_%grid%mesh,1:nspin_) &
                                        + augcharge(1:pawset_%grid%mesh,1:nspin_)
!
!  If there are unbounded scattering states in the pseudopotential generation,
!  n1 and n~1 separately diverge. This makes the hartree and exchange and
!  correlation energies separately diverging. The cancellation becomes
!  very difficult numerically. Outside the sphere however the two charges 
!  should be equal and opposite and we set them to zero. 
!  Is this the right solution? 
!
    do n=pawset_%irc+1,pawset_%grid%mesh
       chargetot=chargeps_(n,1)
       if (nspin_==2) chargetot=chargetot+chargeps_(n,2)
       if (chargetot<1.d-11.or.iflag0==1) then
          charge1_(n,1:nspin_)=0.0_DP
          charge1ps_(n,1:nspin_)=0.0_DP
       endif
    enddo
!    do n=1,pawset_%grid%mesh
!       write(6,'(4f20.15)') pawset_%grid%r(n), chargeps_(n,1), charge1_(n,1),
!                                               charge1ps_(n,1)

  END SUBROUTINE compute_charges
  !
  !============================================================================
  !
  ! Compute the one-center energy and effective potential:
  ! E = Eh[n_v] + Exc[n_v+n_c] - Int[veff*n_v],
  ! veff = vh[n_v] + v_xc[n_v+n_c],
  ! where n_v can be n~+n^ or n1 or n1~+n^ and n_c can be nc or n~c
  !
  SUBROUTINE compute_onecenter_energy ( totenergy_, veff_, &
       pawset_, vcharge_, nlcc_, ccharge_, nspin_, iint, vloc, energies_ , unit_)
    USE funct, ONLY: dft_is_gradient
    USE radial_grids, ONLY: hartree
    USE io_global, ONLY : stdout, ionode
    IMPLICIT NONE
    REAL(dp), INTENT(OUT) :: totenergy_            ! H+XC+DC
    REAL(dp), INTENT(OUT) :: veff_(ndmx,2)         ! effective potential
    TYPE(paw_t),   INTENT(IN)  :: pawset_          ! the PAW dataset
    REAL(dp), INTENT(IN)  :: vcharge_(ndmx,2)      ! valence charge
    LOGICAL,       INTENT(IN)  :: nlcc_            ! non-linear core correction
    REAL(dp), INTENT(IN)  :: ccharge_(ndmx)        ! core charge
    INTEGER,       INTENT(IN)  :: nspin_           ! 1 for LDA, 2 for LSDA
    INTEGER,       INTENT(IN)  :: iint             ! integrals up to iint
    REAL(dp), INTENT(IN), OPTIONAL :: vloc(ndmx)   ! 
    REAL(dp), INTENT(OUT), OPTIONAL :: energies_(5)! Etot,H,XC,DC terms
    INTEGER, OPTIONAL :: unit_
    !
    REAL(dp), PARAMETER :: rho_eq_0(ndmx) = ZERO ! ccharge=0 when nlcc=.f.
    !
    REAL(dp) ::        &
         eh, exc, edc, & ! hartree, xc and double counting energies
         eloc,         & ! local energy
         rhovtot(ndmx), & ! total valence charge
         aux(ndmx),     & ! auxiliary to compute integrals
         vh(ndmx),      & ! hartree potential
         vxc(ndmx,2),   & ! exchange-correlation potential (LDA+GGA)
         vgc(ndmx,2),   & ! exchange-correlation potential (GGA only)
         egc(ndmx),     & ! exchange correlation energy density (GGA only)
         rh(2),        & ! valence charge at a given point without 4 pi r^2
         rhc,          & ! core    charge at a given point without 4 pi r^2
         vxcr(2)         ! exchange-correlation potential at a given point
    REAL(dp) :: & ! compatibility with metaGGA - not yet used
         tau(ndmx) = ZERO, vtau(ndmx) = ZERO  !
    !
    INTEGER :: i, is
    INTEGER :: lsd
    REAL(DP), EXTERNAL :: int_0_inf_dr
#if defined __DEBUG_V_H_vs_SPHEROPOLE
    REAL(DP) :: dummy_charge,aux1(ndmx),aux2(ndmx),res1,res2
#endif
    !
    ! Set up total valence charge
    rhovtot(1:pawset_%grid%mesh) = vcharge_(1:pawset_%grid%mesh,1)
    IF (nspin_==2) rhovtot(1:pawset_%grid%mesh) = rhovtot(1:pawset_%grid%mesh) +   &
         vcharge_(1:pawset_%grid%mesh,2)
    !
    ! Hartree
    CALL hartree(0,2,pawset_%grid%mesh,pawset_%grid,rhovtot,vh)
    if (PRESENT(unit_).and.ionode) then
       write (unit_,*)  " " 
       write (unit_,*)  "#" 
       do i=1,pawset_%grid%mesh
          write (unit_,'(3f12.7)') pawset_%grid%r(i),rhovtot(i),vh(i)
       end do
    end if
#if defined __DEBUG_V_H_vs_SPHEROPOLE
    dummy_charge=int_0_inf_dr(rhovtot,pawset_%grid,pawset_%grid%mesh,2)
    aux1(1:pawset_%grid%mesh) = FPI*pawset_%grid%r2(1:pawset_%grid%mesh)*vh(1:pawset_%grid%mesh) - &
         FPI*dummy_charge * pawset_%grid%r(1:pawset_%grid%mesh)
    aux2(1:pawset_%grid%mesh) = rhovtot(1:pawset_%grid%mesh)*pawset_%grid%r2(1:pawset_%grid%mesh)
    res1 = int_0_inf_dr(aux1,pawset_%grid,pawset_%grid%mesh,1)
    res2 = int_0_inf_dr(aux2,pawset_%grid,pawset_%grid%mesh,4)
    WRITE (stdout,'(4(A,1e15.7))') ' INT rho', dummy_charge,' INT V_H', &
         res1, ' INT r^2*rho', res2, ' ERR:', (1.d0-  res1/ (-res2 * (2.d0*PI/3.d0)))
#endif
    vh(1:pawset_%grid%mesh) = e2 * vh(1:pawset_%grid%mesh)
    aux(1:pawset_%grid%mesh) = vh(1:pawset_%grid%mesh) * rhovtot(1:pawset_%grid%mesh)
    eh = HALF * int_0_inf_dr(aux,pawset_%grid,iint,2)
    !
    ! Exhange-Correlation
    rh=(/ZERO,ZERO/)
    rhc=ZERO
    lsd=nspin_-1
    DO i=1,pawset_%grid%mesh
       DO is=1,nspin_
          rh(is) = vcharge_(i,is)/pawset_%grid%r2(i)/FPI
       ENDDO
       IF (nlcc_) rhc = ccharge_(i)/pawset_%grid%r2(i)/FPI
       CALL vxc_t(lsd,rh,rhc,exc,vxcr)
       vxc(i,1:nspin_)=vxcr(1:nspin_)
       IF (nlcc_) THEN
          aux(i)=exc * (rhovtot(i)+ccharge_(i))
       ELSE
          aux(i)=exc *  rhovtot(i)
       END IF
    END DO
    IF (dft_is_gradient()) THEN
       IF (nlcc_) THEN
          CALL vxcgc(ndmx,pawset_%grid%mesh,nspin_,pawset_%grid%r,&
                     pawset_%grid%r2,vcharge_,ccharge_,vgc,egc, &
                     tau, vtau, 1)
       ELSE
          CALL vxcgc(ndmx,pawset_%grid%mesh,nspin_,pawset_%grid%r,&
                     pawset_%grid%r2,vcharge_,rho_eq_0,vgc,egc, &
                     tau, vtau, 1)
       END IF
       vxc(1:pawset_%grid%mesh,1:nspin_) = vxc(1:pawset_%grid%mesh,1:nspin_) + &
                                      vgc(1:pawset_%grid%mesh,1:nspin_)
       aux(1:pawset_%grid%mesh) = aux(1:pawset_%grid%mesh) + &
           egc(1:pawset_%grid%mesh) * pawset_%grid%r2(1:pawset_%grid%mesh) * FPI
    END IF
    exc = int_0_inf_dr(aux,pawset_%grid,iint,2)
    !
    ! Double counting
    edc=ZERO
    DO is=1,nspin_
       veff_(1:pawset_%grid%mesh,is)=vxc(1:pawset_%grid%mesh,is)+vh(1:pawset_%grid%mesh)
       aux(1:pawset_%grid%mesh)=veff_(1:pawset_%grid%mesh,is)*vcharge_(1:pawset_%grid%mesh,is)
       edc=edc+int_0_inf_dr(aux,pawset_%grid,iint,2)
    END DO
    !
    eloc=ZERO
    !
    IF (present(vloc)) THEN
       DO is=1,nspin_
          aux(1:pawset_%grid%mesh)=vloc(1:pawset_%grid%mesh)        &
                               *vcharge_(1:pawset_%grid%mesh,is)
          eloc=eloc+int_0_inf_dr(aux,pawset_%grid,iint,2)
       ENDDO
    ENDIF
    !
    ! Total
    totenergy_ = eh + exc - edc
    !
    IF (PRESENT(energies_)) THEN
       energies_(1)=totenergy_
       energies_(2)=eh
       energies_(3)=exc
       energies_(4)=edc
       energies_(5)=eloc
    END IF
    !
  END SUBROUTINE compute_onecenter_energy
  !
  !============================================================================
  !
  ! Compute NL 'D' coefficients = D^ + D1 - D1~
  !
  SUBROUTINE compute_nonlocal_coeff(ddd_, pawset_, nspin_, veffps_, veff1_, veff1ps_)
    IMPLICIT NONE
    REAL(dp), INTENT(OUT) :: ddd_(nwfsx,nwfsx,2)
    TYPE(paw_t),   INTENT(IN)  :: pawset_
    INTEGER,       INTENT(IN)  :: nspin_
    REAL(dp), INTENT(IN)  :: veffps_(ndmx,2)
    REAL(dp), INTENT(IN)  :: veff1_(ndmx,2)
    REAL(dp), INTENT(IN)  :: veff1ps_(ndmx,2)
    INTEGER :: is, ns, ns1
    REAL(dp) :: aux(ndmx), dd
    REAL(DP), EXTERNAL :: int_0_inf_dr
!     REAL(dp):: dddd(nwfsx,nwfsx,3) = 0.d0

    !
    ! D^ = Int Q*v~
    ! D1-D1~ = kindiff + Int[ae*v1*ae - ps*v1~*ps - Q*v1~]
    ddd_(:,:,:)=ZERO
    DO is=1,nspin_
       DO ns=1,pawset_%nwfc
          DO ns1=1,ns
             IF (pawset_%l(ns)==pawset_%l(ns1).and.&
                       ABS(pawset_%jj(ns)-pawset_%jj(ns1))<1.d-8) THEN
                ! Int[Q*v~]
                aux(1:pawset_%grid%mesh) =                        &
                   pawset_%augfun(1:pawset_%grid%mesh,ns,ns1,0) * &
                   veffps_(1:pawset_%grid%mesh,is)
                dd = int_0_inf_dr(aux,pawset_%grid,pawset_%irc,(pawset_%l(ns)+1)*2)
!                 dddd(ns,ns1,1) = int_0_inf_dr(aux,pawset_%grid,pawset_%irc,(pawset_%l(ns)+1)*2)
                ! Int[ae*v1*ae]
                aux(1:pawset_%grid%mesh) =                        &
                     pawset_%aewfc(1:pawset_%grid%mesh,ns ) *     &
                     pawset_%aewfc(1:pawset_%grid%mesh,ns1) *     &
                     veff1_(1:pawset_%grid%mesh,is)
                IF (pawset_%rel==2) &
                    aux(1:pawset_%irc) =  aux(1:pawset_%irc) +     &
                     pawset_%aewfc_rel(1:pawset_%irc,ns ) *     &
                     pawset_%aewfc_rel(1:pawset_%irc,ns1) *     &
                     veff1_(1:pawset_%irc,is)
                dd = dd +                                    &
                     int_0_inf_dr(aux,pawset_%grid,pawset_%irc,(pawset_%l(ns)+1)*2)
!                 dddd(ns,ns1,2) = int_0_inf_dr(aux,pawset_%grid,pawset_%irc,(pawset_%l(ns)+1)*2)
                ! Int[ps*v1~*ps + aufun*v1~]
                aux(1:pawset_%grid%mesh) =                        &
                   ( pawset_%pswfc(1:pawset_%grid%mesh,ns ) *     &
                     pawset_%pswfc(1:pawset_%grid%mesh,ns1) +     &
                     pawset_%augfun(1:pawset_%grid%mesh,ns,ns1,0) ) * &
                     veff1ps_(1:pawset_%grid%mesh,is)
                dd = dd -                                    &
                     int_0_inf_dr(aux,pawset_%grid,pawset_%irc,(pawset_%l(ns)+1)*2)
!                 dddd(ns,ns1,3) = int_0_inf_dr(aux,pawset_%grid,pawset_%irc,(pawset_%l(ns)+1)*2)
                ! collect
                ddd_(ns,ns1,is) = pawset_%kdiff(ns,ns1) + dd
                ddd_(ns1,ns,is) = ddd_(ns,ns1,is)
             END IF
          END DO
       END DO
    END DO
!     write(*,*) 'deeq', pawset_%irc
!     write(*,'(4f13.8)') dddd(1:4,1:4,1)
!     write(*,*) 'ddd ae'
!     write(*,'(4f13.8)') dddd(1:4,1:4,2)
!     write(*,*) 'ddd ps'
!     write(*,'(4f13.8)') dddd(1:4,1:4,3)

  END SUBROUTINE compute_nonlocal_coeff
  !
  !============================================================================
  !
  ! 'D_ion' coefficients = D1 - D1~
  !
  SUBROUTINE compute_nonlocal_coeff_ion(ddd_, pawset_)
    IMPLICIT NONE
    REAL(dp), INTENT(OUT) :: ddd_(nwfsx,nwfsx)
    TYPE(paw_t),   INTENT(IN)  :: pawset_
    INTEGER :: ns, ns1
    REAL(dp) :: aux(ndmx), dd
    REAL(DP), EXTERNAL :: int_0_inf_dr
    !
    ! D^ = Int Q*v~
    ! D1-D1~ = kindiff + Int[ae*v1*ae - ps*v1~*ps - Q*v1~]
!     write(666,"(4f12.6)") (pawset_%grid%r(ns), veffps_(ns,1), veff1_(ns,1), veff1ps_(ns,1), ns=1,pawset_%grid%mesh)
!     write(667,"(4f12.6)") (pawset_%grid%r(ns), veffps_(ns,2), veff1_(ns,2), veff1ps_(ns,2), ns=1,pawset_%grid%mesh)
    ddd_(:,:)=ZERO
    DO ns=1,pawset_%nwfc
       DO ns1=1,ns
          IF (pawset_%l(ns)==pawset_%l(ns1).and. &
                      ABS(pawset_%jj(ns)-pawset_%jj(ns1))<1.d-8 ) THEN
             ! Int[ae*v1*ae]
             aux(1:pawset_%grid%mesh) =                        &
                  pawset_%aewfc(1:pawset_%grid%mesh,ns ) *     &
                  pawset_%aewfc(1:pawset_%grid%mesh,ns1) *     &
                  pawset_%aeloc(1:pawset_%grid%mesh)
             IF (pawset_%rel==2) &
             aux(1:pawset_%irc) = aux(1:pawset_%irc) +         &
                  pawset_%aewfc_rel(1:pawset_%irc,ns ) *     &
                  pawset_%aewfc_rel(1:pawset_%irc,ns1) *     &
                  pawset_%aeloc(1:pawset_%irc)
             dd = int_0_inf_dr(aux,pawset_%grid,pawset_%irc,(pawset_%l(ns)+1)*2)
             ! Int[ps*v1~*ps + Q*v1~]
             aux(1:pawset_%grid%mesh) =                        &
                ( pawset_%pswfc(1:pawset_%grid%mesh,ns ) *     &
                  pawset_%pswfc(1:pawset_%grid%mesh,ns1) +     &
                  pawset_%augfun(1:pawset_%grid%mesh,ns,ns1,0) ) * &
                  pawset_%psloc(1:pawset_%grid%mesh)
             dd = dd - &
                  int_0_inf_dr(aux,pawset_%grid,pawset_%irc,(pawset_%l(ns)+1)*2)
             !
             ddd_(ns,ns1) = pawset_%kdiff(ns,ns1) +  dd
             ddd_(ns1,ns)=ddd_(ns,ns1)
          END IF
       END DO
    END DO
  END SUBROUTINE compute_nonlocal_coeff_ion
  !
  !============================================================================
  !
  ! Write PAW dataset wfc and potentials on files
  !
  SUBROUTINE human_write_paw(pawset_)
    USE io_global, ONLY : ionode
    IMPLICIT NONE
    TYPE(paw_t), INTENT(In) :: pawset_
    INTEGER :: n,ns
    IF (.not.ionode) return
    DO ns=1,pawset_%nwfc
       WRITE (5000+ns,'(A)') "# r                 AEwfc               PSwfc               projector"
       DO n=1,pawset_%grid%mesh
          WRITE (5000+ns,'(4f20.12)') pawset_%grid%r(n), pawset_%aewfc(n,ns), pawset_%pswfc(n,ns),pawset_%proj(n,ns)
       END DO
    END DO
    !
    WRITE (6000,'(A)') "# r                 AEpot               PSpot"
    DO n=1,pawset_%grid%mesh
       WRITE (6000,'(3f20.8)') pawset_%grid%r(n), pawset_%aeloc(n), pawset_%psloc(n)
    END DO
  END SUBROUTINE human_write_paw
  !
  !============================================================================
  !                       LOWER-LEVEL ROUTINES                              !!!
  !============================================================================
  !
  ! Compute Sum_i oc_i * wfc_i^2
  !
  SUBROUTINE compute_sumwfc2(charge_, pawset_, nwfc_, wfc_, oc_, spin_)
    IMPLICIT NONE
    REAL(dp), INTENT(OUT) :: charge_(ndmx,2)
    TYPE(paw_t),   INTENT(IN)  :: pawset_
    INTEGER,       INTENT(IN)  :: nwfc_
    REAL(dp), INTENT(IN)  :: wfc_(ndmx,nwfsx)
    REAL(dp), INTENT(IN)  :: oc_(nwfsx)
    INTEGER,       INTENT(IN)  :: spin_(nwfsx)
    INTEGER :: nf
    charge_(1:pawset_%grid%mesh,:)=ZERO
    DO nf=1,nwfc_
       IF (oc_(nf)>ZERO) charge_(1:pawset_%grid%mesh,spin_(nf)) = &
                         charge_(1:pawset_%grid%mesh,spin_(nf)) + oc_(nf) * &
                         wfc_(1:pawset_%grid%mesh,nf) * wfc_(1:pawset_%grid%mesh,nf)
    END DO
  END SUBROUTINE compute_sumwfc2
  !
  !============================================================================
  !
  ! Compute Sum_n oc_n <pswfc_n|proj_i> <proj_j|pswfc_n>
  !
  SUBROUTINE compute_projsum (projsum_, pawset_, nwfc_, l_, j_, spin_, pswfc_, oc_)
    REAL(dp), INTENT(OUT) :: projsum_(nwfsx,nwfsx,2)
    TYPE(paw_t),   INTENT(IN)  :: pawset_
    INTEGER,       INTENT(IN)  :: nwfc_
    INTEGER,       INTENT(IN)  :: l_(nwfsx)
    INTEGER,       INTENT(IN)  :: spin_(nwfsx)
    REAL(dp),      INTENT(IN)  :: j_(nwfsx)
    REAL(dp), INTENT(IN)  :: pswfc_(ndmx,nwfsx)
    REAL(dp), INTENT(IN)  :: oc_(nwfsx)
    REAL(dp) :: proj_dot_wfc(nwfsx,nwfsx), aux(ndmx)
    INTEGER :: ns, ns1, nf, nr, is, nst
    REAL(DP), EXTERNAL :: int_0_inf_dr
    ! Compute <projector|wavefunction>
    DO ns=1,pawset_%nwfc
       DO nf=1,nwfc_
          IF (pawset_%l(ns)==l_(nf).AND.pawset_%jj(ns)==j_(nf)) THEN
             DO nr=1,pawset_%grid%mesh
                aux(nr)=pawset_%proj(nr,ns)*pswfc_(nr,nf)
             END DO
             nst=(l_(nf)+1)*2
             proj_dot_wfc(ns,nf)=int_0_inf_dr(aux,pawset_%grid,pawset_%ikk(ns),nst)
          ELSE
             proj_dot_wfc(ns,nf)=ZERO
          END IF
       END DO
    END DO
    ! Do the sum over the wavefunctions
    projsum_(:,:,:)=ZERO
    DO ns=1,pawset_%nwfc
       DO ns1=1,ns
          DO nf=1,nwfc_
             is=spin_(nf)
             IF (oc_(nf)>ZERO) THEN
                projsum_(ns,ns1,is) = projsum_(ns,ns1,is) + oc_(nf) * &
                     proj_dot_wfc(ns,nf) * proj_dot_wfc(ns1,nf)
             END IF
          END DO
          projsum_(ns1,ns,:)=projsum_(ns,ns1,:)
       END DO
    END DO
  END SUBROUTINE compute_projsum
  !
  !============================================================================
  !
  !
  ! Compute augmentation charge as Sum_ij W_ij * Q_ij
  !
  SUBROUTINE compute_augcharge(augcharge_, pawset_, projsum_, nspin_)
    IMPLICIT NONE
    REAL(dp), INTENT(OUT) :: augcharge_(ndmx,2)
    TYPE(paw_t),   INTENT(IN)  :: pawset_
    REAL(dp), INTENT(IN)  :: projsum_(nwfsx,nwfsx,2)
    INTEGER,       INTENT(IN)  :: nspin_
    INTEGER :: ns, ns1, is
    REAL(dp) :: factor
    augcharge_=ZERO
    DO is=1,nspin_
       DO ns=1,pawset_%nwfc
          DO ns1=1,ns
             ! multiply by TWO off-diagonal terms
             factor=TWO
             IF (ns1==ns) factor=ONE
             !
             ! NB: in a spherically averaged system only the l=0 component 
             !     of the augmentation functions is present
             !
             augcharge_(1:pawset_%grid%mesh,is) = augcharge_(1:pawset_%grid%mesh,is) + &
                    factor * projsum_(ns,ns1,is) * &
                    pawset_%augfun(1:pawset_%grid%mesh,ns,ns1,0)
          END DO
       END DO
    END DO
  END SUBROUTINE compute_augcharge
  !
  !============================================================================
  !
  ! Compute n1 and n1~, as Sum_ij W_ij wfc_i(r)*wfc_j(r) 
  !
  SUBROUTINE compute_onecenter_charge(charge1_, pawset_, projsum_, nspin_, which_wfc)
    IMPLICIT NONE
    REAL(dp),   INTENT(OUT) :: charge1_(ndmx,2)
    TYPE(paw_t),     INTENT(IN)  :: pawset_
    REAL(dp),   INTENT(IN)  :: projsum_(nwfsx,nwfsx,2)
    INTEGER,         INTENT(IN)  :: nspin_
    CHARACTER(LEN=2),INTENT(IN)  :: which_wfc
    INTEGER :: ns, ns1, is
    REAL(dp) :: factor
    charge1_=ZERO
    DO is=1,nspin_
       DO ns=1,pawset_%nwfc
          DO ns1=1,ns
             ! multiply times TWO off-diagonal terms
             IF (ns1==ns) THEN
                factor=ONE
             ELSE
                factor=TWO
             END IF
             SELECT CASE (which_wfc)
             CASE ("AE")
                charge1_(1:pawset_%grid%mesh,is) = charge1_(1:pawset_%grid%mesh,is) + factor * &
                     projsum_(ns,ns1,is) * pawset_%aewfc(1:pawset_%grid%mesh,ns) *     &
                     pawset_%aewfc(1:pawset_%grid%mesh,ns1)
                IF (pawset_%rel==2) &
                charge1_(1:pawset_%irc,is) = charge1_(1:pawset_%irc,is) + factor * &
                     projsum_(ns,ns1,is) * pawset_%aewfc_rel(1:pawset_%irc,ns) *     &
                     pawset_%aewfc_rel(1:pawset_%irc,ns1)
             CASE ("PS")
                charge1_(1:pawset_%grid%mesh,is) = charge1_(1:pawset_%grid%mesh,is) + factor * &
                     projsum_(ns,ns1,is) * pawset_%pswfc(1:pawset_%grid%mesh,ns) *     &
                     pawset_%pswfc(1:pawset_%grid%mesh,ns1)
             CASE DEFAULT
                call errore ('compute_onecenter_charge','specify AE or PS wavefunctions',1)
             END SELECT
          END DO
       END DO
    END DO
  END SUBROUTINE compute_onecenter_charge
!
!--------------------------------------------------------------------------
SUBROUTINE find_bes_qi(qc,rmatch,lam,ncn,iok)
  !--------------------------------------------------------------------------
  !
  !      This routine finds two values of q such that the
  !      functions f_l have a derivative equal to 0 at rmatch
  !  
  IMPLICIT NONE

  INTEGER,     INTENT(IN)  ::      &
       lam,   & ! input: the angular momentum
       ncn      ! input: the number of qi to compute
  INTEGER,     INTENT(INOUT)  ::      &
       iok      ! output: if 0 the calculation in this routine is ok

  REAL (dp),   INTENT(OUT)   :: &
       qc(ncn)  ! output: the values of qi
  REAL (dp),   INTENT(IN)  :: rmatch

  REAL (dp) ::   &
       zeroderjl (2,7) ! first two zeros of the first derivative of 
                       ! spherical Bessel function j_l for l = 0,...,6

  INTEGER ::    &
       nc       ! counter on the q found

  data zeroderjl / 0.0_dp,                 4.4934094579614_dp, &
                   2.0815759780862_dp,     5.9403699890844_dp, &
                   3.3420936578747_dp,     7.2899322987026_dp, &
                   4.5140996477983_dp,     8.5837549433127_dp, &
                   5.6467036213923_dp,     9.8404460168549_dp, &
                   6.7564563311363_dp,    11.0702068269176_dp, &
                   7.8510776799611_dp,    12.2793339053177_dp  /
  iok=0
  IF (ncn.gt.2) &
       CALL errore('find_aug_qi','ncn is too large',1)

  IF (lam.gt.6) &
       CALL errore('find_aug_qi','l not programmed',1)
  !
  !    fix deltaq and the maximum step number
  !
  DO nc = 1, ncn
     !
     qc(nc) = zeroderjl (nc, lam + 1) / rmatch
     !
  ENDDO
  RETURN
END SUBROUTINE find_bes_qi 
  !
  !============================================================================
  !
END MODULE atomic_paw

