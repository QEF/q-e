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
  !
  !============================================================================
  !
  USE io_global, ONLY: stdout
  USE kinds, ONLY: dp
  USE constants, ONLY: pi, fpi, e2, eps8
  USE ld1_parameters, ONLY: ndm, nwfsx
  USE parameters, ONLY: lmaxx
  !
  IMPLICIT NONE
  PRIVATE
  SAVE
  !
  REAL(dp), PARAMETER :: ZERO=0._dp, ONE=1._dp, TWO=2._dp, HALF=0.5_dp
  !
  !============================================================================
  !
  ! Type describing a PAW dataset
  !
  TYPE :: paw_t
     REAL (dp) :: zval
     REAL (dp) :: z
     REAL (dp) :: zmesh
     CHARACTER(LEN=80) :: dft
     INTEGER        :: mesh      ! the size of the mesh
     REAL (dp) :: r (ndm)   ! the mesh
     REAL (dp) :: r2 (ndm)  ! r^2
     REAL (dp) :: sqr (ndm) ! sqrt(r)
     REAL (dp) :: dx        ! log(r(i+1))-log(r(i))
     LOGICAL :: nlcc ! nonlinear core correction
     INTEGER :: nwfc ! number of wavefunctions/projectors
     INTEGER :: lmax ! maximum angular momentum of projectors
     INTEGER :: l(nwfsx) ! angular momentum of projectors
     INTEGER :: ikk(nwfsx) ! cutoff radius for the projectors
     INTEGER :: irc ! r(irc) = radius of the augmentation sphere
     REAL (dp) :: &
          oc (nwfsx), & ! the occupations
          enl (nwfsx), & ! the energy of the wavefunctions
          aewfc (ndm,nwfsx), &  ! all-electron wavefunctions
          pswfc (ndm,nwfsx),        & ! pseudo wavefunctions
          proj (ndm,nwfsx),     & ! projectors
          augfun(ndm,nwfsx,nwfsx),      & ! augmentation functions
          augmom(nwfsx,nwfsx,0:2*lmaxx) , & ! moments of the augmentation functions
          aeccharge (ndm),  & ! AE core charge
          psccharge (ndm),  & ! PS core charge
          aeloc (ndm),     & ! descreened AE potential: v_AE-v_H[n1]-v_XC[n1+nc]
          psloc (ndm),     & ! descreened local PS potential: v_PS-v_H[n~+n^]-v_XC[n~+n^+n~c]
          kdiff (nwfsx,nwfsx)         ! kinetic energy differences
!!!  Notes about screening:
!!!       Without nlcc, the local PSpotential is descreened with n~+n^ only.
!!!       The local AEpotential is descreened ALWAYS with n1+nc. This improves
!!!       the accuracy, and will not cost in the plane wave code (atomic
!!!       contribution only).
  END TYPE paw_t
  !
  !============================================================================
  !
  PUBLIC :: paw_t
  PUBLIC :: us2paw
  PUBLIC :: paw2us
  PUBLIC :: new_paw_hamiltonian
  PUBLIC :: paw_io
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
       pawset_, nwfc_, l_, nspin_, spin_, oc_, pswfc_, eig_)
    IMPLICIT NONE
    REAL(dp), INTENT(OUT) :: veffps_(ndm,2)
    REAL(dp), INTENT(OUT) :: ddd_(nwfsx,nwfsx,2)
    REAL(dp), INTENT(OUT) :: etot_
    TYPE(paw_t),   INTENT(IN)  :: pawset_
    INTEGER,       INTENT(IN)  :: nwfc_
    INTEGER,       INTENT(IN)  :: l_(nwfsx)
    INTEGER,       INTENT(IN)  :: nspin_
    INTEGER,       INTENT(IN)  :: spin_(nwfsx)
    REAL(dp), INTENT(IN)  :: oc_(nwfsx)
    REAL(dp), INTENT(IN)  :: pswfc_(ndm,nwfsx)
    REAL(dp), INTENT(IN)  :: eig_(nwfsx)
    !
    REAL(dp) :: &                                        ! one center:
         eps,             e1,             e1ps,             & ! energies;
                          veff1(ndm,2),   veff1ps(ndm,2),   & ! eff potentials;
         chargeps(ndm,2), charge1(ndm,2), charge1ps(ndm,2), & ! charges.
         projsum(nwfsx,nwfsx,2), eigsum !  sum of projections, sum of eigenval.
    !
    INTEGER :: ns, ns1, is
    REAL(dp) :: aux(ndm)
    !
    ! Compute the valence charges
    CALL compute_charges(projsum, chargeps, charge1, charge1ps, &
       pawset_, nwfc_, l_, nspin_, spin_, oc_, pswfc_ )
    !
    ! Compute the effective potentials (H+XC)
    CALL compute_onecenter_energy ( eps,  veffps_, &
       pawset_, chargeps,  pawset_%nlcc, pawset_%psccharge, nspin_ )
    CALL compute_onecenter_energy ( e1,   veff1, &
       pawset_, charge1,  .TRUE.,        pawset_%aeccharge, nspin_ )
    CALL compute_onecenter_energy ( e1ps, veff1ps, &
       pawset_, charge1ps, pawset_%nlcc, pawset_%psccharge, nspin_ )
    ! Add the local part
    DO is=1,nspin_
       veffps_(1:pawset_%mesh,is) = veffps_(1:pawset_%mesh,is) +  &
            pawset_%psloc(1:pawset_%mesh)
       veff1  (1:pawset_%mesh,is) = veff1  (1:pawset_%mesh,is) +  &
            pawset_%aeloc(1:pawset_%mesh)
       veff1ps(1:pawset_%mesh,is) = veff1ps(1:pawset_%mesh,is) +  &
            pawset_%psloc(1:pawset_%mesh)
    END DO
    !
    ! Compute the nonlocal D coefficients
    CALL compute_nonlocal_coeff (ddd_,pawset_,nspin_,veffps_,veff1,veff1ps)
    !
    ! Compute total energy
    eigsum=ZERO
    DO ns=1,nwfc_
       IF (oc_(ns)>ZERO) eigsum = eigsum + oc_(ns)*eig_(ns)
    END DO
    etot_ = eigsum + eps + e1 - e1ps
    !
  END SUBROUTINE new_paw_hamiltonian
  !
  !============================================================================
  !
  ! Convert the USPP to a PAW dataset
  !
  SUBROUTINE us2paw (pawset_,                                     &
       zval, mesh, r, r2, sqr, dx, irc, ikk,                      &
       nbeta, lls, ocs, enls, psipaw, phis, betas, qvan, kindiff, &
       nlcc, aerhoc, psrhoc, aevtot,psvtot, do_write_dataset)
    USE funct, only : get_iexch,get_icorr,get_igcx,get_igcc,dft_name
    IMPLICIT NONE
    TYPE(paw_t),   INTENT(OUT) :: pawset_
    REAL(dp), INTENT(IN)  :: zval
    INTEGER,       INTENT(IN)  :: mesh
    REAL(dp), INTENT(IN)  :: r(ndm)
    REAL(dp), INTENT(IN)  :: r2(ndm), sqr(ndm), dx
    INTEGER,       INTENT(IN)  :: irc
    INTEGER,       INTENT(IN)  :: ikk(nwfsx)
    INTEGER,       INTENT(IN)  :: nbeta
    INTEGER,       INTENT(IN)  :: lls(nwfsx)
    REAL(dp), INTENT(IN)  :: ocs(nwfsx)
    REAL(dp), INTENT(IN)  :: enls(nwfsx)
    REAL(dp), INTENT(IN)  :: psipaw(ndm,nwfsx)
    REAL(dp), INTENT(IN)  :: phis(ndm,nwfsx)
    REAL(dp), INTENT(IN)  :: betas(ndm,nwfsx)
    REAL(dp), INTENT(IN)  :: qvan(ndm,nwfsx,nwfsx)
    REAL(dp), INTENT(IN)  :: kindiff(nwfsx,nwfsx)
    LOGICAL,       INTENT(IN)  :: nlcc
    REAL(dp), INTENT(IN)  :: aerhoc(ndm)
    REAL(dp), INTENT(IN)  :: psrhoc(ndm)
    REAL(dp), INTENT(IN)  :: aevtot(ndm)
    REAL(dp), INTENT(IN)  :: psvtot(ndm)
    LOGICAL,INTENT(IN),OPTIONAL:: do_write_dataset
    !
    REAL(DP), EXTERNAL :: int_0_inf_dr
    REAL(dp) :: vps(ndm,2), projsum(nwfsx,nwfsx,2), ddd(nwfsx,nwfsx,2)
    INTEGER :: ns, ns1, n, l
    REAL(dp) :: aux(ndm), aux2(ndm,2), raux
    REAL(dp) :: aecharge(ndm,2), pscharge(ndm,2)
    REAL(dp) :: etot
    INTEGER :: nspin=1, spin(nwfsx)=1
    CHARACTER(LEN=4) :: shortname
    INTEGER :: iexch, icorr, igcx, igcc
    !
    pawset_%zval=zval
    !
    ! Copy the mesh
    pawset_%mesh=mesh
    pawset_%r(1:mesh)=r(1:mesh)
    pawset_%r2(1:mesh)=r2(1:mesh)
    pawset_%sqr(1:mesh)=sqr(1:mesh)
    pawset_%dx=dx
    pawset_%irc=irc
    pawset_%ikk(1:nbeta)=ikk(1:nbeta)
    !
    ! Copy the wavefunctions
    pawset_%nwfc=nbeta
    pawset_%l(1:nbeta)=lls(1:nbeta)
    pawset_%lmax=MAXVAL(pawset_%l(1:pawset_%nwfc))
    pawset_%oc(1:nbeta)=ocs(1:nbeta)
    pawset_%enl(1:nbeta)=enls(1:nbeta)
    pawset_%aewfc(1:mesh,1:nbeta) = psipaw(1:mesh,1:nbeta)
    pawset_%pswfc(1:mesh,1:nbeta) = phis  (1:mesh,1:nbeta)
    pawset_%proj (1:mesh,1:nbeta) = betas (1:mesh,1:nbeta)
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
    ! Compute the exact augmentation functions
    DO ns=1,nbeta
       DO ns1=1,ns
          pawset_%augfun(1:mesh,ns,ns1) = &
              pawset_%aewfc(1:mesh,ns) * pawset_%aewfc(1:mesh,ns1) - &
              pawset_%pswfc(1:mesh,ns) * pawset_%pswfc(1:mesh,ns1)
          pawset_%augfun(1:mesh,ns1,ns) = pawset_%augfun(1:mesh,ns,ns1)
       END DO
    END DO
    !
    ! Compute the moments of the exact augmentation functions
    DO l=0,2*pawset_%lmax
       DO ns=1,nbeta
          DO ns1=1,ns
             aux(1:pawset_%irc) = pawset_%augfun(1:pawset_%irc,ns,ns1) * pawset_%r(1:pawset_%irc)**l
             pawset_%augmom(ns,ns1,l) = int_0_inf_dr(aux(1:pawset_%irc), &
                  pawset_%r,pawset_%r2,pawset_%dx,pawset_%irc,(pawset_%l(ns)+1)*2)
             pawset_%augmom(ns1,ns,l)=pawset_%augmom(ns,ns1,l)
          END DO
       END DO
    END DO
    !
!!! Uncomment the following line if you wish to test the invariance with
!!! respect to the specific shape of the radial part of the augmentation
!!! functions (the following implementation of this test is correct for
!!! spherical symmetry only, ie only the zero-th moment is conserved)
!#define __TEST_AUGFUN__
#if defined __TEST_AUGFUN__
    CALL infomsg ('us2paw','Use for tests only!', -1)
    ! a gaussian
    aux(1:mesh)=EXP(-(r(1:mesh)**2)/(TWO*0.25_dp**2))
    DO ns=1,nbeta
       DO ns1=1,ns
          IF ( (pawset_%l(ns1)==pawset_%l(ns)) .AND. (ABS(pawset_%augmom(ns,ns1,0))>eps8)) THEN
             !
             ! Choose the shape of the augmentation functions: NC Q ...
             pawset_%augfun(1:mesh,ns,ns1) = qvan(1:mesh,ns,ns1)
             ! ... or Gaussian Q ?
             !pawset_%augfun(1:mesh,ns,ns1) = aux(1:mesh) * pawset_%r2(1:mesh)
             !
             ! Renormalize the aug. functions so that their integral is correct
             raux=int_0_inf_dr(pawset_%augfun(1:pawset_%irc,ns,ns1),pawset_%r,pawset_%r2,pawset_%dx,pawset_%irc,(pawset_%l(ns)+1)*2)
             IF (ABS(raux) < eps8) THEN
                CALL errore('us2paw','norm of augmentation funciton too small',ns*100+ns1)
             END IF
             raux=pawset_%augmom(ns,ns1,0)/raux
             !pawset_%augfun(1:mesh,ns,ns1)=raux*pawset_%augfun(1:mesh,ns,ns1)
          ELSE
             pawset_%augfun(1:mesh,ns,ns1)=ZERO
          END IF
          ! Set the symmetric part
          pawset_%augfun(1:mesh,ns1,ns)=pawset_%augfun(1:mesh,ns,ns1)
          !WRITE (900+ns*10+ns1,'(2e20.10)')(r(n),pawset_%augfun(n,ns,ns1),n=1,irc)
       END DO
    END DO
#endif
    !
    !
    ! Copy kinetic energy differences
    pawset_%kdiff(1:nbeta,1:nbeta)=kindiff(1:nbeta,1:nbeta)
    !
    ! Copy the core charge (if not NLCC copy only the AE one)
    pawset_%nlcc=nlcc
    pawset_%aeccharge(1:mesh)=aerhoc(1:mesh)
    IF (pawset_%nlcc) pawset_%psccharge(1:mesh)=psrhoc(1:mesh)
    !
    ! Copy the local potentials
    pawset_%aeloc(1:mesh)=aevtot(1:mesh)
    pawset_%psloc(1:mesh)=psvtot(1:mesh)
    ! and descreen them:
    CALL compute_charges(projsum, pscharge, aecharge, aux2, &
       pawset_, nbeta, lls, nspin, spin, ocs, phis )
    CALL compute_onecenter_energy ( raux,  aux2, &
       pawset_, pscharge, pawset_%nlcc, pawset_%psccharge, nspin )
    pawset_%psloc(1:mesh)=psvtot(1:mesh)-aux2(1:mesh,1)
    CALL compute_onecenter_energy ( raux,  aux2, &
       pawset_, aecharge, .TRUE.,       pawset_%aeccharge, nspin )
    pawset_%aeloc(1:mesh)=aevtot(1:mesh)-aux2(1:mesh,1)
    !WRITE(4444,'(5e20.10)')(r(n),aevtot(n),psvtot(n),pawset_%aeloc(n),pawset_%psloc(n),n=1,mesh)
    !
    pawset_%dft="                                                                                "
    iexch = get_iexch()
    icorr = get_icorr()
    igcx  = get_igcx()
    igcc  = get_igcc()
    CALL dft_name (iexch, icorr, igcx, igcc, pawset_%dft, shortname)
    !
    IF (PRESENT(do_write_dataset)) THEN
       IF (do_write_dataset) CALL human_write_paw(pawset_)
    END IF
    !
    ! Generate the paw hamiltonian for test (should be equal to the US one)
    CALL new_paw_hamiltonian (vps, ddd, etot, &
       pawset_, pawset_%nwfc, pawset_%l, nspin, spin, pawset_%oc, pawset_%pswfc, pawset_%enl)
    WRITE(stdout,'(/5x,A,f12.6,A)') 'Estimated PAW energy =',etot,' Ry'
    WRITE(stdout,'(/5x,A)') 'The PAW screened D coefficients'
    DO ns1=1,pawset_%nwfc
       WRITE(stdout,'(6f12.5)') (ddd(ns1,ns,1),ns=1,pawset_%nwfc)
    END DO
    !
  END SUBROUTINE us2paw
  !
  !============================================================================
  !
  ! ...
  !
  SUBROUTINE paw2us (pawset_,zval,mesh,r,r2,sqr,dx,nbeta,lls,ikk, &
       betas,qq,qvan,pseudotype)
    USE funct, ONLY : set_dft_from_name
    IMPLICIT NONE
    TYPE(paw_t),   INTENT(IN)  :: pawset_
    REAL(dp), INTENT(OUT) :: zval
    INTEGER,       INTENT(OUT) :: mesh
    REAL(dp), INTENT(OUT) :: r(ndm)
    REAL(dp), INTENT(OUT) :: r2(ndm), sqr(ndm), dx
    INTEGER,       INTENT(OUT) :: nbeta
    INTEGER,       INTENT(OUT) :: lls(nwfsx)
    INTEGER,       INTENT(OUT) :: ikk(nwfsx)
    INTEGER,       INTENT(OUT) :: pseudotype
    REAL(dp), INTENT(OUT) :: betas(ndm,nwfsx)
    REAL(dp), INTENT(OUT) :: qq(nwfsx,nwfsx)
    REAL(dp), INTENT(OUT) :: qvan(ndm,nwfsx,nwfsx)
    INTEGER :: ns, ns1
    !
    zval=pawset_%zval
    !
    mesh=pawset_%mesh
    r(1:mesh)=pawset_%r(1:mesh)
    r2(1:mesh)=pawset_%r2(1:mesh)
    sqr(1:mesh)=pawset_%sqr(1:mesh)
    dx=pawset_%dx
    !
    nbeta=pawset_%nwfc
    lls(1:nbeta)=pawset_%l(1:nbeta)
    ikk(1:nbeta)=pawset_%ikk(1:nbeta)
    !
    DO ns=1,nbeta
       DO ns1=1,nbeta
          IF (lls(ns)==lls(ns1)) THEN
             qvan(1:mesh,ns,ns1)=pawset_%augfun(1:mesh,ns,ns1)
             qq(ns,ns1)=pawset_%augmom(ns,ns1,0)
          ELSE ! enforce spherical symmetry
             qvan(1:mesh,ns,ns1)=0._dp
             qq(ns,ns1)=0._dp
          END IF
       END DO
    END DO
    !
    betas(1:mesh,1:nbeta)=pawset_%proj(1:mesh,1:nbeta)
    pseudotype=3
    !
    CALL set_dft_from_name (pawset_%dft)
  END SUBROUTINE paw2us
  !
  !============================================================================
  !
  ! Read/write the PAW dataset
  !
  SUBROUTINE paw_io(pawset_,un,what)
    TYPE(paw_t), INTENT(INOUT) :: pawset_
    INTEGER, INTENT(IN) :: un
    CHARACTER(LEN=3), INTENT(IN) :: what
    INTEGER :: n, ns, ns1, l
    SELECT CASE (what)
    CASE ("OUT")
       WRITE(un,*)
       WRITE(un,'(e20.10)') pawset_%zval
       WRITE(un,'(i8)') pawset_%mesh
       WRITE(un,'(e20.10)') (pawset_%r(n), n=1,pawset_%mesh)
       WRITE(un,'(e20.10)') (pawset_%r2(n), n=1,pawset_%mesh)
       WRITE(un,'(e20.10)') (pawset_%sqr(n), n=1,pawset_%mesh)
       WRITE(un,'(e20.10)') pawset_%dx
       WRITE(un,*) pawset_%nlcc
       WRITE(un,'(i8)') pawset_%nwfc
       WRITE(un,'(i8)') (pawset_%l(ns), ns=1,pawset_%nwfc)
       WRITE(un,'(i8)') (pawset_%ikk(ns), ns=1,pawset_%nwfc)
       WRITE(un,'(i8)') pawset_%irc
       WRITE(un,'(e20.10)') (pawset_%oc(ns), ns=1,pawset_%nwfc)
       WRITE(un,'(e20.10)') (pawset_%enl(ns), ns=1,pawset_%nwfc)
       WRITE(un,'(e20.10)') ((pawset_%aewfc(n,ns), n=1,pawset_%mesh), ns=1,pawset_%nwfc)
       WRITE(un,'(e20.10)') ((pawset_%pswfc(n,ns), n=1,pawset_%mesh), ns=1,pawset_%nwfc)
       WRITE(un,'(e20.10)') ((pawset_%proj(n,ns), n=1,pawset_%mesh), ns=1,pawset_%nwfc)
       WRITE(un,'(e20.10)') (((pawset_%augfun(n,ns,ns1), n=1,pawset_%mesh), ns=1,pawset_%nwfc), ns1=1,pawset_%nwfc)
       WRITE(un,'(i8)') pawset_%lmax
       WRITE(un,'(e20.10)') (((pawset_%augmom(ns,ns1,l), ns=1,pawset_%nwfc), ns1=1,pawset_%nwfc),l=0,2*pawset_%lmax)
       WRITE(un,'(e20.10)') (pawset_%aeccharge(n), n=1,pawset_%mesh)
       IF (pawset_%nlcc) WRITE(un,'(e20.10)') (pawset_%psccharge(n), n=1,pawset_%mesh)
       WRITE(un,'(e20.10)') (pawset_%aeloc(n), n=1,pawset_%mesh)
       WRITE(un,'(e20.10)') (pawset_%psloc(n), n=1,pawset_%mesh)
       WRITE(un,'(e20.10)') ((pawset_%kdiff(ns,ns1), ns=1,pawset_%nwfc), ns1=1,pawset_%nwfc)
       WRITE(un,'(A)') TRIM(pawset_%dft)
    CASE ("INP")
       READ(un,*)
       READ(un,'(e20.10)') pawset_%zval
       READ(un,'(i8)') pawset_%mesh
       READ(un,'(e20.10)') (pawset_%r(n), n=1,pawset_%mesh)
       READ(un,'(e20.10)') (pawset_%r2(n), n=1,pawset_%mesh)
       READ(un,'(e20.10)') (pawset_%sqr(n), n=1,pawset_%mesh)
       READ(un,'(e20.10)') pawset_%dx
       READ(un,*) pawset_%nlcc
       READ(un,'(i8)') pawset_%nwfc
       READ(un,'(i8)') (pawset_%l(ns), ns=1,pawset_%nwfc)
       READ(un,'(i8)') (pawset_%ikk(ns), ns=1,pawset_%nwfc)
       READ(un,'(i8)') pawset_%irc
       READ(un,'(e20.10)') (pawset_%oc(ns), ns=1,pawset_%nwfc)
       READ(un,'(e20.10)') (pawset_%enl(ns), ns=1,pawset_%nwfc)
       READ(un,'(e20.10)') ((pawset_%aewfc(n,ns), n=1,pawset_%mesh), ns=1,pawset_%nwfc)
       READ(un,'(e20.10)') ((pawset_%pswfc(n,ns), n=1,pawset_%mesh), ns=1,pawset_%nwfc)
       READ(un,'(e20.10)') ((pawset_%proj(n,ns), n=1,pawset_%mesh), ns=1,pawset_%nwfc)
       READ(un,'(e20.10)') (((pawset_%augfun(n,ns,ns1), n=1,pawset_%mesh), ns=1,pawset_%nwfc), ns1=1,pawset_%nwfc)
       READ(un,'(i8)') pawset_%lmax
       READ(un,'(e20.10)') (((pawset_%augmom(ns,ns1,l), ns=1,pawset_%nwfc), ns1=1,pawset_%nwfc),l=0,2*pawset_%lmax)
       READ(un,'(e20.10)') (pawset_%aeccharge(n), n=1,pawset_%mesh)
       IF (pawset_%nlcc) READ(un,'(e20.10)') (pawset_%psccharge(n), n=1,pawset_%mesh)
       READ(un,'(e20.10)') (pawset_%aeloc(n), n=1,pawset_%mesh)
       READ(un,'(e20.10)') (pawset_%psloc(n), n=1,pawset_%mesh)
       READ(un,'(e20.10)') ((pawset_%kdiff(ns,ns1), ns=1,pawset_%nwfc), ns1=1,pawset_%nwfc)
       pawset_%dft="                                                                                "
       READ(un,'(A)') pawset_%dft
    CASE DEFAULT
       CALL errore ('paw_io','specify (INP)ut or (OUT)put',1)
    END SELECT
    CLOSE(un)
  END SUBROUTINE paw_io
  !
  !============================================================================
  !                          PRIVATE ROUTINES                               !!!
  !============================================================================
  !
  SUBROUTINE compute_charges (projsum_, chargeps_, charge1_, charge1ps_, &
       pawset_, nwfc_, l_, nspin_, spin_, oc_, pswfc_ )
    IMPLICIT NONE
    REAL(dp), INTENT(OUT) :: projsum_(nwfsx,nwfsx,2)
    REAL(dp), INTENT(OUT) :: chargeps_(ndm,2)
    REAL(dp), INTENT(OUT) :: charge1_(ndm,2)
    REAL(dp), INTENT(OUT) :: charge1ps_(ndm,2)
    TYPE(paw_t),   INTENT(IN)  :: pawset_
    INTEGER,       INTENT(IN)  :: nwfc_
    INTEGER,       INTENT(IN)  :: l_(nwfsx)
    INTEGER,       INTENT(IN)  :: nspin_
    INTEGER,       INTENT(IN)  :: spin_(nwfsx)
    REAL(dp), INTENT(IN)  :: oc_(nwfsx)
    REAL(dp), INTENT(IN)  :: pswfc_(ndm,nwfsx)
    REAL(dp) :: augcharge(ndm,2)
    CALL compute_projsum(projsum_,pawset_,nwfc_,l_,spin_,pswfc_,oc_)
    !WRITE (6200,'(20e20.10)') ((projsum(ns,ns1),ns=1,ns1),ns1=1,pawset_%nwfc)
    CALL compute_sumwfc2(chargeps_,pawset_,nwfc_,pswfc_,oc_,spin_)
    CALL compute_onecenter_charge(charge1ps_,pawset_,projsum_,nspin_,"PS")
    CALL compute_onecenter_charge(charge1_  ,pawset_,projsum_,nspin_,"AE")
    ! add augmentation charges
    CALL compute_augcharge(augcharge,pawset_,projsum_,nspin_)
    chargeps_ (1:pawset_%mesh,1:nspin_) = chargeps_ (1:pawset_%mesh,1:nspin_) + &
         augcharge(1:pawset_%mesh,1:nspin_)
    charge1ps_(1:pawset_%mesh,1:nspin_) = charge1ps_(1:pawset_%mesh,1:nspin_) + &
         augcharge(1:pawset_%mesh,1:nspin_)
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
       pawset_, vcharge_, nlcc_, ccharge_, nspin_, energies_ )
    USE funct, ONLY: dft_is_gradient
    IMPLICIT NONE
    REAL(dp), INTENT(OUT) :: totenergy_       ! H+XC+DC
    REAL(dp), INTENT(OUT) :: veff_(ndm,2)     ! effective potential
    TYPE(paw_t),   INTENT(IN)  :: pawset_          ! the PAW dataset
    REAL(dp), INTENT(IN)  :: vcharge_(ndm,2)  ! valence charge
    LOGICAL,       INTENT(IN)  :: nlcc_            ! non-linear core correction
    REAL(dp), INTENT(IN)  :: ccharge_(ndm)    ! core charge
    INTEGER,       INTENT(IN)  :: nspin_           ! 1 for LDA, 2 for LSDA
    REAL(dp), INTENT(OUT), OPTIONAL :: energies_(4) ! Etot,H,XC,DC terms
    !
    REAL(dp), PARAMETER :: rho_eq_0(ndm) = ZERO ! ccharge=0 when nlcc=.f.
    !
    REAL(dp) :: &
         eh, exc, edc, & ! hartree, xc and double counting energies
         rhovtot(ndm), & ! total valence charge
         aux(ndm),     & ! auxiliary to compute integrals
         vh(ndm),      & ! hartree potential
         vxc(ndm,2),   & ! exchange-correlation potential (LDA+GGA)
         vgc(ndm,2),   & ! exchange-correlation potential (GGA only)
         egc(ndm),     & ! exchange correlation energy density (GGA only)
         rh(2),        & ! valence charge at a given point without 4 pi r^2
         rhc,          & ! core    charge at a given point without 4 pi r^2
         vxcr(2)         ! exchange-correlation potential at a given point
    !
    INTEGER :: ns, i, is
    INTEGER :: lsd
    REAL(DP), EXTERNAL :: int_0_inf_dr, exc_t
    !
    ! Set up total valence charge
    rhovtot(1:pawset_%mesh) = vcharge_(1:pawset_%mesh,1)
    IF (nspin_==2) rhovtot(1:pawset_%mesh) = rhovtot(1:pawset_%mesh) +   &
         vcharge_(1:pawset_%mesh,2)
    !
    ! Hartree
    CALL hartree(0,2,pawset_%mesh,pawset_%r,pawset_%r2,pawset_%sqr,pawset_%dx,rhovtot,vh)
    vh(1:pawset_%mesh) = e2 * vh(1:pawset_%mesh)
    aux(1:pawset_%mesh) = vh(1:pawset_%mesh) * rhovtot(1:pawset_%mesh)
    eh = HALF * int_0_inf_dr(aux,pawset_%r,pawset_%r2,pawset_%dx,pawset_%mesh,2)
    !
    ! Exhange-Correlation
    rh=(/ZERO,ZERO/)
    rhc=ZERO
    lsd=nspin_-1
    DO i=1,pawset_%mesh
       DO is=1,nspin_
          rh(is) = vcharge_(i,is)/pawset_%r2(i)/FPI
       ENDDO
       IF (nlcc_) rhc = ccharge_(i)/pawset_%r2(i)/FPI
       CALL vxc_t(rh,rhc,lsd,vxcr)
       vxc(i,1:nspin_)=vxcr(1:nspin_)
       IF (nlcc_) THEN
          aux(i)=exc_t(rh,rhc,lsd) * (rhovtot(i)+ccharge_(i))
       ELSE
          aux(i)=exc_t(rh,rhc,lsd) *  rhovtot(i)
       END IF
    END DO
    IF ( dft_is_gradient() ) THEN
       IF (nlcc_) THEN
          CALL vxcgc(ndm,pawset_%mesh,nspin_,pawset_%r,pawset_%r2,vcharge_,ccharge_,vgc,egc)
       ELSE
          CALL vxcgc(ndm,pawset_%mesh,nspin_,pawset_%r,pawset_%r2,vcharge_,rho_eq_0,vgc,egc)
       END IF
       vxc(1:pawset_%mesh,1:nspin_) = vxc(1:pawset_%mesh,1:nspin_) + &
                                      vgc(1:pawset_%mesh,1:nspin_)
       aux(1:pawset_%mesh) = aux(1:pawset_%mesh) + &
           egc(1:pawset_%mesh) * pawset_%r2(1:pawset_%mesh) * FPI
    END IF
    exc = int_0_inf_dr(aux,pawset_%r,pawset_%r2,pawset_%dx,pawset_%mesh,2)
    !
    ! Double counting
    edc=ZERO
    DO is=1,nspin_
       veff_(1:pawset_%mesh,is)=vxc(1:pawset_%mesh,is)+vh(1:pawset_%mesh)
       aux(1:pawset_%mesh)=veff_(1:pawset_%mesh,is)*vcharge_(1:pawset_%mesh,is)
       edc=edc+int_0_inf_dr(aux,pawset_%r,pawset_%r2,pawset_%dx,pawset_%mesh,2)
    END DO
    !
    ! Total
    totenergy_ = eh + exc - edc
    !
    IF (PRESENT(energies_)) THEN
       energies_(1)=totenergy_
       energies_(2)=eh
       energies_(3)=exc
       energies_(4)=edc
    END IF
    !
  END SUBROUTINE compute_onecenter_energy
  !
  !============================================================================
  !
  ! Compute NL 'D' coefficients = D^ + D1 - D~1
  !
  SUBROUTINE compute_nonlocal_coeff(ddd_, pawset_, nspin_, veffps_, veff1_, veff1ps_)
    IMPLICIT NONE
    REAL(dp), INTENT(OUT) :: ddd_(nwfsx,nwfsx,2)
    TYPE(paw_t),   INTENT(IN)  :: pawset_
    INTEGER,       INTENT(IN)  :: nspin_
    REAL(dp), INTENT(IN)  :: veffps_(ndm,2)
    REAL(dp), INTENT(IN)  :: veff1_(ndm,2)
    REAL(dp), INTENT(IN)  :: veff1ps_(ndm,2)
    INTEGER :: is, ns, ns1
    REAL(dp) :: aux(ndm)
    REAL(DP), EXTERNAL :: int_0_inf_dr
    !
    ! D^ = Int Q*v~
    ! D1-D1~ = kindiff + Int[ae*v1*ae - ps*v1~*ps - augfun*v~1]
    DO is=1,nspin_
       DO ns=1,pawset_%nwfc
          DO ns1=1,ns
             IF (pawset_%l(ns)==pawset_%l(ns1)) THEN
                aux(1:pawset_%mesh) =                        &
                     pawset_%augfun(1:pawset_%mesh,ns,ns1) * &
                     veffps_(1:pawset_%mesh,is)
                aux(1:pawset_%mesh) = aux(1:pawset_%mesh) +  &
                     pawset_%aewfc(1:pawset_%mesh,ns ) *     &
                     pawset_%aewfc(1:pawset_%mesh,ns1) *     &
                     veff1_(1:pawset_%mesh,is)
                aux(1:pawset_%mesh) = aux(1:pawset_%mesh) -  &
                     pawset_%pswfc(1:pawset_%mesh,ns ) *     &
                     pawset_%pswfc(1:pawset_%mesh,ns1) *     &
                     veff1ps_(1:pawset_%mesh,is)
                aux(1:pawset_%mesh) = aux(1:pawset_%mesh) -  &
                     pawset_%augfun(1:pawset_%mesh,ns,ns1) * &
                     veff1ps_(1:pawset_%mesh,is)
                ddd_(ns,ns1,is) = pawset_%kdiff(ns,ns1) +       &
                     int_0_inf_dr(aux,pawset_%r,pawset_%r2,pawset_%dx,pawset_%irc,(pawset_%l(ns)+1)*2)
             ELSE
                ddd_(ns,ns1,is)=ZERO
             END IF
             ddd_(ns1,ns,is)=ddd_(ns,ns1,is)
          END DO
       END DO
    END DO
  END SUBROUTINE compute_nonlocal_coeff
  !
  !============================================================================
  !
  ! Write PAW dataset wfc and potentials on files
  !
  SUBROUTINE human_write_paw(pawset_)
    IMPLICIT NONE
    TYPE(paw_t), INTENT(In) :: pawset_
    INTEGER :: n,ns
    DO ns=1,pawset_%nwfc
       WRITE (5000+ns,'(A)') "# r                 AEwfc               PSwfc               projector"
       DO n=1,pawset_%mesh
          WRITE (5000+ns,'(4f20.12)') pawset_%r(n), pawset_%aewfc(n,ns), pawset_%pswfc(n,ns),pawset_%proj(n,ns)
       END DO
    END DO
    !
    WRITE (6000,'(A)') "# r                 AEpot               PSpot"
    DO n=1,pawset_%mesh
       WRITE (6000,'(3f20.8)') pawset_%r(n), pawset_%aeloc(n), pawset_%psloc(n)
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
    REAL(dp), INTENT(OUT) :: charge_(ndm,2)
    TYPE(paw_t),   INTENT(IN)  :: pawset_
    INTEGER,       INTENT(IN)  :: nwfc_
    REAL(dp), INTENT(IN)  :: wfc_(ndm,nwfsx)
    REAL(dp), INTENT(IN)  :: oc_(nwfsx)
    INTEGER,       INTENT(IN)  :: spin_(nwfsx)
    INTEGER :: ns
    charge_(1:pawset_%mesh,:)=ZERO
    DO ns=1,nwfc_
       IF (oc_(ns)>ZERO) charge_(1:pawset_%mesh,spin_(ns)) = charge_(1:pawset_%mesh,spin_(ns)) + &
            oc_(ns) * wfc_(1:pawset_%mesh,ns) * wfc_(1:pawset_%mesh,ns)
    END DO
  END SUBROUTINE compute_sumwfc2
  !
  !============================================================================
  !
  ! Compute Sum_n oc_n <pswfc_n|proj_i> <proj_j|pswfc_n>
  !
  SUBROUTINE compute_projsum (projsum_, pawset_, nwfc_, l_, spin_, pswfc_, oc_)
    REAL(dp), INTENT(OUT) :: projsum_(nwfsx,nwfsx,2)
    TYPE(paw_t),   INTENT(IN)  :: pawset_
    INTEGER,       INTENT(IN)  :: nwfc_
    INTEGER,       INTENT(IN)  :: l_(nwfsx)
    INTEGER,       INTENT(IN)  :: spin_(nwfsx)
    REAL(dp), INTENT(IN)  :: pswfc_(ndm,nwfsx)
    REAL(dp), INTENT(IN)  :: oc_(nwfsx)
    REAL(dp) :: proj_dot_wfc(nwfsx,nwfsx), aux(ndm)
    INTEGER :: ns, ns1, nf, nr, is, nst
    REAL(DP), EXTERNAL :: int_0_inf_dr
    ! Compute <projector|wavefunction>
    DO ns=1,pawset_%nwfc
       DO nf=1,nwfc_
          IF (pawset_%l(ns)==l_(nf)) THEN
             DO nr=1,pawset_%mesh
                aux(nr)=pawset_%proj(nr,ns)*pswfc_(nr,nf)
             END DO
             nst=(l_(nf)+1)*2
             proj_dot_wfc(ns,nf)=int_0_inf_dr(aux,pawset_%r,pawset_%r2,pawset_%dx,pawset_%irc,nst)
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
    REAL(dp), INTENT(OUT) :: augcharge_(ndm,2)
    TYPE(paw_t),   INTENT(IN)  :: pawset_
    REAL(dp), INTENT(IN)  :: projsum_(nwfsx,nwfsx,2)
    INTEGER,       INTENT(IN)  :: nspin_
    INTEGER :: ns, ns1, is
    REAL(dp) :: factor
    augcharge_=ZERO
    DO is=1,nspin_
       DO ns=1,pawset_%nwfc
          DO ns1=1,ns
             ! multiply times TWO off-diagonal terms
             IF (ns1==ns) THEN
                factor=ONE
             ELSE
                factor=TWO
             END IF
             augcharge_(1:pawset_%mesh,is) = augcharge_(1:pawset_%mesh,is) + factor * &
                  projsum_(ns,ns1,is) * pawset_%augfun(1:pawset_%mesh,ns,ns1)
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
    REAL(dp),   INTENT(OUT) :: charge1_(ndm,2)
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
                charge1_(1:pawset_%mesh,is) = charge1_(1:pawset_%mesh,is) + factor * &
                     projsum_(ns,ns1,is) * pawset_%aewfc(1:pawset_%mesh,ns) *     &
                     pawset_%aewfc(1:pawset_%mesh,ns1)
             CASE ("PS")
                charge1_(1:pawset_%mesh,is) = charge1_(1:pawset_%mesh,is) + factor * &
                     projsum_(ns,ns1,is) * pawset_%pswfc(1:pawset_%mesh,ns) *     &
                     pawset_%pswfc(1:pawset_%mesh,ns1)
             CASE DEFAULT
                call errore ('compute_onecehnter_charge','specify AE or PS wavefunctions',1)
             END SELECT
          END DO
       END DO
    END DO
  END SUBROUTINE compute_onecenter_charge
  !
  !============================================================================
  !
END MODULE atomic_paw
