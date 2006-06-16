!
! Copyright (C) 2002-2005 FPMD-CPV groups
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

#include "f_defs.h"

MODULE rundiis_module

        IMPLICIT NONE
        SAVE

        PRIVATE

        LOGICAL, PARAMETER :: tforce  = .FALSE.
        LOGICAL, PARAMETER :: tstress = .FALSE.

        INTERFACE set_lambda
          MODULE PROCEDURE set_lambda_r
        END INTERFACE

        PUBLIC :: rundiis, runsdiis


CONTAINS


   SUBROUTINE rundiis( tprint, rhoe, atoms, &
                 bec, becdr, eigr, ei1, ei2, ei3, sfac, c0, cm, cgrad, cdesc, tcel, ht0, fi, eig, &
                 vpot, doions, edft )

   !  this routine computes the electronic ground state via diagonalization
   !  by the DIIS method (Wood and Zunger, J.Phys.A 18,1343 (1985))
   !  resulting wave functions are the Kohn-Sham eigenfunctions
   !
   !  overview of the DIIS method:
   !  1) make a starting guess on the ground-state eigenfunctions, |A(0)>
   !  2) iterate:
   !     a) compute a new difference vector, |dA(n+1)>, by solving
   !        the equation
   !
   !          (H-E(n))|A(n)+dA(n+1)> = 0,  E(n) = <A(n)|H|A(n)>
   !
   !        in the "diagonal approximation"
   !
   !                                    <x(i)|H-E(n)|A(n)>
   !          |dA(n+1)> = -(sum over i) ------------------ |x(i)>
   !                                    <x(i)|H-E(n)|x(i)>
   !
   !        where the |x(i)> are suitable basis vectors
   !     b) compute a new approximate eigenvector, |A(n+1)>, as a linear
   !        combination of all |dA>'s (including |dA(0)> = |A(0)>), with
   !        coefficients chosen as to minimize the norm of the vector
   !
   !          |R(n+1)> = (H-E(n+1))|A(n+1)>,  E(n+1) = <A(n+1)|H|A(n+1)>
   !
   !        (which is exactly zero when |A(n+1)> is an eigenvector of H):
   !        they are obtained as the eigenvector of the lowest eigenvalue
   !        of equation
   !
   !          P|c> = lambda Q|c>
   !
   !          P(i,j) = <(H-E(n))dA(i)|(H-E(n))dA(j)>
   !          Q(i,j) = <dA(i)|dA(j)>
   !
   !        this equation has a small dimensionality (n+2) and is easily
   !        solved by standard techniques
   !  3) stop when <R|R> is zero within a given tolerance

   ! ... declare modules

      USE kinds
      USE mp_global, ONLY: me_image, nproc_image, intra_image_comm
      USE mp, ONLY: mp_sum
      USE io_global, ONLY: ionode
      USE io_global, ONLY: stdout
      USE energies, ONLY: dft_energy_type
      USE cp_electronic_mass, ONLY: emass
      USE electrons_base, ONLY: nupdwn, iupdwn
      USE time_step, ONLY: delt
      USE wave_functions, ONLY: proj, crot
      USE phase_factors_module, ONLY: strucf, phfacs
      USE charge_mix
      USE charge_density, ONLY: rhoofr
      USE guess
      USE diis
      USE cell_module, ONLY: boxdimensions
      USE check_stop, ONLY: check_stop_now
      USE nl, ONLY: nlrh_m
      USE potentials, ONLY: vofrhos
      USE forces
      USE wave_types, ONLY: wave_descriptor
      USE atoms_type_module, ONLY: atoms_type
      USE control_flags, ONLY: force_pairing, gamma_only
      use grid_dimensions,    only: nr1, nr2, nr3
      USE reciprocal_vectors, ONLY: mill_l
      USE gvecp, ONLY: ngm
      USE local_pseudo, ONLY: vps
      USE uspp,             ONLY : vkb, nkb

      IMPLICIT NONE

      ! ... declare subroutine arguments
      !
      LOGICAL    :: tcel, tprint, doions
      TYPE (atoms_type) :: atoms
      COMPLEX(DP), INTENT(INOUT) :: c0(:,:,:), cm(:,:,:), cgrad(:,:,:)
      TYPE (wave_descriptor) :: cdesc
      REAL(DP) :: rhoe(:,:)
      REAL(DP) :: bec(:,:)
      REAL(DP) :: becdr(:,:,:)
      COMPLEX(DP) :: sfac(:,:)
      COMPLEX(DP) :: eigr(:,:)
      COMPLEX(DP) :: ei1(:,:)
      COMPLEX(DP) :: ei2(:,:)
      COMPLEX(DP) :: ei3(:,:)
      TYPE (boxdimensions), INTENT(INOUT) ::  ht0
      REAL(DP)  :: fi(:,:)
      TYPE (dft_energy_type) :: edft

      REAL(DP)    :: eig(:,:)
      REAL(DP)    :: vpot(:,:) 

      ! ... declare other variables
      !
      INTEGER ig, ib, j, k, ngw, i, is, nrt, istate, nrl, ndiis, nowv
      INTEGER idiis
      LOGICAL tlimit
      REAL(DP)  fions(3) 
      REAL(DP)  timepre,s1,s2,s3,s4,s5
      REAL(DP)  dene,etot_m,cnorm, drho
      REAL(DP)  ekinc,svar1,svar2,svar3_0
      REAL(DP)  efermi, sume, entk, rhos, eold, dum_kin, nel
      REAL(DP)  wke( cdesc%ldb, cdesc%nkl, cdesc%nspin )
      REAL(DP),    ALLOCATABLE :: lambda(:,:), fs(:,:)
      COMPLEX(DP), ALLOCATABLE :: clambda(:,:,:)
      COMPLEX(DP), ALLOCATABLE :: c0rot(:,:)

      REAL(DP), EXTERNAL :: cclock


! ... end of declarations
!  ----------------------------------------------

      IF( .NOT. gamma_only ) &
         CALL errore( " rundiis", " diis and k-points not allowed ", 1 )

      IF( force_pairing ) &
        CALL errore(' rundiis ', ' force pairing not implemented ', 1 )

 
      nrt   = 0         ! electronic minimization at fixed potential

      ! ... Initialize DIIS

      treset_diis = .TRUE.
      doions      = .FALSE.
      istate      = 0
 
      ! ... Initialize energy values

      etot_m   = 0.0d0
      cnorm    = 1000.0d0
      ekinc    = 100.0d0
      edft%ent = 0.0d0
      eold     = 1.0d10  ! a large number
      nel = 2.0d0 * cdesc%nbt( 1 ) / DBLE( cdesc%nspin )

      IF( SIZE( fi, 2 ) > 1 ) THEN
        CALL errore(' rundiis ', ' nspin > 1 not allowed ', SIZE( fi, 2 ) )
      END IF

      ! ... initialize occupation numbers

      ALLOCATE( fs( cdesc%ldb, cdesc%nspin ) )
      fs = 2.d0


      ! ... distribute lambda's rows across processors with a blocking factor
      ! ... of 1, ( row 1 to PE 1, row 2 to PE 2, .. row nproc_image+1 to PE 1 and
      ! ... so on).

      ! ... compute local number of rows
      nrl = cdesc%nbl( 1 ) / nproc_image
      IF( me_image < MOD( cdesc%nbl( 1 ), nproc_image) ) THEN
        nrl = nrl + 1
      END IF

! ... Allocate and initialize lambda to the identity matrix
      ALLOCATE(lambda(nrl,cdesc%nbl( 1 )))
      CALL set_lambda(lambda)


! ... starting guess on the wavefunctions
!      CALL guessrho(rhoe, cm, c0, fi, ht0)

      CALL phfacs( ei1, ei2, ei3, eigr, mill_l, atoms%taus, nr1, nr2, nr3, atoms%nat )
      CALL strucf( sfac, ei1, ei2, ei3, mill_l, ngm )
      CALL rhoofr( 1, c0, cdesc, fi, rhoe, ht0)
      CALL newrho(rhoe(:,1), drho, 0)  ! memorize density
      CALL phfacs( ei1, ei2, ei3, eigr, mill_l, atoms%taus, nr1, nr2, nr3, atoms%nat )
      CALL strucf( sfac, ei1, ei2, ei3, mill_l, ngm )
      CALL guessc0( .NOT. gamma_only, bec, c0, cm, cdesc)

! ... Initialize the rotation index srot
      srot = srot0

      DIIS_LOOP: DO idiis = 1, maxnstep

        nrt   = nrt   + 1  ! number of steps since last upgrade of density

        s1 = cclock()

! ...   density upgrade every srot steps
        IF ( nrt >= srot .OR. idiis == 1 ) THEN

          nrt    = 0  ! reset the charge rotation counter
          istate = 0  ! reset the diis internal state

          IF( idiis /= 1 ) THEN

             ! ...       bring wave functions onto KS states

             CALL crot( 1, c0(:,:,1), cdesc, lambda, eig(:,1) )

             call adjef_s(eig(1,1),fi(1,1),efermi,nel, cdesc%nbl( 1 ),temp_elec,sume)
             call entropy_s(fi(1,1),temp_elec,cdesc%nbl( 1 ),edft%ent)

          END IF

! ...     self consistent energy
          edft%enl = nlrh_m(c0, cdesc, tforce, atoms, bec, becdr, eigr)
          CALL rhoofr( 1, c0, cdesc, fi, rhoe, ht0)
          CALL vofrhos(.FALSE., tforce, tstress, rhoe, atoms, &
            vpot, bec, c0, cdesc, fi, eigr, ei1, ei2, ei3, sfac, timepre, ht0, edft)

! ...     density upgrade
          CALL newrho(rhoe(:,1), drho, idiis)
          IF (ionode) WRITE( stdout,45) idiis, edft%etot, drho
          dene   = abs(edft%etot - etot_m)
          etot_m = edft%etot

45        FORMAT('etot  drho ',i3,1x,2(1x,f18.10))

! ...     recalculate potential
          edft%enl = nlrh_m(c0, cdesc, tforce, atoms, bec, becdr, eigr)
          CALL vofrhos(.FALSE., tforce, tstress, rhoe, atoms, &
            vpot, bec, c0, cdesc, fi, eigr, ei1, ei2, ei3, sfac, timepre, ht0, edft)

          IF( idiis /= 1 )THEN
            IF( drho < tolrhof .AND. dene < tolene) EXIT DIIS_LOOP
            IF( drho < tolrho)  srot = sroti
            IF( drho < tolrhoi) srot = srotf
          END IF

! ...     check for exit
          IF ( check_stop_now() ) THEN
            cm = c0
            EXIT DIIS_LOOP
          END IF

! ...     calculate lambda_i,j=<c_i| H |c_j>

          edft%enl = nlrh_m(c0, cdesc, tforce, atoms, bec, becdr, eigr)

          CALL dforce_all( 1, c0(:,:,1), fi(:,1), cgrad(:,:,1), vpot(:,1), eigr, bec, nupdwn, iupdwn )

          CALL proj( 1, cgrad(:,:,1), cdesc, c0(:,:,1), cdesc, lambda )
          CALL crot( 1, c0(:,:,1), cdesc, lambda, eig(:,1) )

          call adjef_s(eig(1,1),fi(1,1),efermi,nel, cdesc%nbl( 1 ),temp_elec,sume)
          call entropy_s(fi(1,1),temp_elec,cdesc%nbl(1),edft%ent)

          edft%enl = nlrh_m(c0, cdesc, tforce, atoms, bec, becdr, eigr)
          CALL dforce_all( 1, c0(:,:,1), fi(:,1), cgrad(:,:,1), vpot(:,1), eigr, bec, nupdwn, iupdwn )

          DO ib = 1, cdesc%nbl( 1 )
            cgrad(:,ib,1) = cgrad(:,ib,1) + eig(ib,1)*c0(:,ib,1)
          END DO

        ELSE

! ...     DIIS on c0 at FIXED potential
          edft%enl = nlrh_m(c0, cdesc, tforce, atoms, bec, becdr, eigr)

          CALL dforce_all( 1, c0(:,:,1), fi(:,1), cgrad(:,:,1), vpot(:,1), eigr, bec, nupdwn, iupdwn )

          CALL proj( 1, cgrad(:,:,1), cdesc, c0(:,:,1), cdesc, lambda)

        END IF

        edft%etot = 0.d0
        DO ib=1,nrl
           edft%etot = edft%etot + lambda(ib,(ib-1)*nproc_image+me_image+1)
        END DO
        CALL mp_sum(edft%etot, intra_image_comm)

        IF (ionode) WRITE( stdout,80) idiis, cnorm, edft%etot, edft%ent
80      FORMAT("STEP NORMG ETOT ENT: ",I3,2X,F12.8,2X,F16.6,4(1x,f8.5))

        s4 = cclock()
        svar1   = 2.d0
        svar2   = -1.d0
        svar3_0 = delt * delt / emass
        CALL simupd(ekinc,doions,c0(:,:,1),cgrad(:,:,1),cdesc, svar1,svar2, &
                  svar3_0,edft%etot,fs(:,1),eigr,sfac,vps, &
                  treset_diis,istate,cnorm,eold,ndiis,nowv)

        CALL gram( vkb, bec, nkb, c0(1,1,1), SIZE(c0,1), cdesc%nbt( 1 ) )

      END DO DIIS_LOOP

      IF( idiis > maxnstep )THEN
        IF (ionode) THEN
          WRITE( stdout,fmt = ' (3X, "NOT CONVERGED IN ",I5," STEPS") ' ) maxnstep
          WRITE( stdout,fmt = ' (3X, "DRHO FINAL = ",D14.6) ' ) drho
        END IF
        cm = c0
      END IF

      IF ( check_stop_now() ) THEN
        IF (ionode) THEN
          WRITE( stdout,fmt = ' (3X, "CPU LIMIT REACHED") ' ) 
        END IF
        cm = c0
      END IF

      IF( drho < tolrhof .AND. dene < tolene) THEN
        IF (ionode) THEN
          WRITE( stdout,fmt = ' (3X, "CONVERGENCE ACHIEVED") ' ) 
        END IF
      END IF

      IF( tprint ) THEN
         WHERE( fi(:,1) /= 0.d0 ) eig(:,1) = eig(:,1) / fi(:,1)
      END IF

      DEALLOCATE(lambda)
      deallocate(fs)

      RETURN

      END SUBROUTINE rundiis




      SUBROUTINE runsdiis(tprint, rhoe, atoms, &
                 bec, becdr, eigr, ei1, ei2, ei3, sfac, c0, cm, cgrad, cdesc, tcel, ht0, fi, eig, &
                 vpot, doions, edft )

!  this routine computes the electronic ground state via diagonalization
!  by the DIIS method (Wood and Zunger, J.Phys.A 18,1343 (1985))
!  resulting wave functions are the Kohn-Sham eigenfunctions
!
!  overview of the DIIS method:
!  1) make a starting guess on the ground-state eigenfunctions, |A(0)>
!  2) iterate:
!     a) compute a new difference vector, |dA(n+1)>, by solving
!        the equation
!
!          (H-E(n))|A(n)+dA(n+1)> = 0,  E(n) = <A(n)|H|A(n)>
!
!        in the "diagonal approximation"
!
!                                    <x(i)|H-E(n)|A(n)>
!          |dA(n+1)> = -(sum over i) ------------------ |x(i)>
!                                    <x(i)|H-E(n)|x(i)>
!
!        where the |x(i)> are suitable basis vectors
!     b) compute a new approximate eigenvector, |A(n+1)>, as a linear
!        combination of all |dA>'s (including |dA(0)> = |A(0)>), with
!        coefficients chosen as to minimize the norm of the vector
!
!          |R(n+1)> = (H-E(n+1))|A(n+1)>,  E(n+1) = <A(n+1)|H|A(n+1)>
!
!        (which is exactly zero when |A(n+1)> is an eigenvector of H):
!        they are obtained as the eigenvector of the lowest eigenvalue
!        of equation
!
!          P|c> = lambda Q|c>
!
!          P(i,j) = <(H-E(n))dA(i)|(H-E(n))dA(j)>
!          Q(i,j) = <dA(i)|dA(j)>
!
!        this equation has a small dimensionality (n+2) and is easily
!        solved by standard techniques
!  3) stop when <R|R> is zero within a given tolerance

! ... declare modules
      USE kinds
      USE mp_global, ONLY: me_image, nproc_image
      USE runcp_module, ONLY: runcp
      USE energies, ONLY: dft_energy_type
      USE electrons_module, ONLY: ei
      USE cp_electronic_mass, ONLY: emass
      USE electrons_base, ONLY: nupdwn, iupdwn
      USE time_step, ONLY: delt
      USE wave_functions, ONLY: proj, update_wave_functions
      USE diis
      USE cell_module, ONLY: boxdimensions
      USE check_stop, ONLY: check_stop_now
      USE potentials, ONLY: vofrhos, kspotential
      USE forces
      USE io_global, ONLY: ionode
      USE io_global, ONLY: stdout
      USE control_flags, ONLY: tortho, tsde
      USE wave_types
      USE atoms_type_module, ONLY: atoms_type
      USE local_pseudo, ONLY: vps
      USE uspp,             ONLY : vkb, nkb

      IMPLICIT NONE

! ... declare subroutine arguments
      LOGICAL   :: tprint, tcel, doions
      TYPE (atoms_type) :: atoms
      COMPLEX(DP), INTENT(INOUT) :: c0(:,:,:), cm(:,:,:), cgrad(:,:,:)
      TYPE (wave_descriptor) :: cdesc
      REAL(DP) :: rhoe(:,:)
      COMPLEX(DP) :: eigr(:,:)
      COMPLEX(DP) :: ei1(:,:)
      COMPLEX(DP) :: ei2(:,:)
      COMPLEX(DP) :: ei3(:,:)
      COMPLEX(DP) :: sfac(:,:)
      TYPE (boxdimensions), INTENT(INOUT) ::  ht0
      REAL(DP)  :: fi(:,:)
      REAL(DP)  :: bec(:,:)
      REAL(DP)  :: becdr(:,:,:)
      TYPE (dft_energy_type) :: edft

      REAL(DP)    :: eig(:,:)
      REAL(DP)    :: vpot(:,:)

! ... declare other variables
      LOGICAL :: tlimit, tsteep
      LOGICAL :: ttreset_diis(SIZE(fi, 2))
      LOGICAL :: ddoions(SIZE(fi, 2))
      INTEGER ig, ib, ibg, j, k, ibl, ispin
      INTEGER nfi_l,nrt,istate
      INTEGER nspin
      INTEGER nx,nrl, ndiis,nowv, isteep
      REAL(DP)  ekinc(2), svar1, svar2, svar3_0
      REAL(DP)  timepre, s0, s1, s2, s3, s4, s5
      REAL(DP) :: seconds_per_iter, old_clock_value
      REAL(DP)  dene, etot_m, cnorm,  drho
      REAL(DP)  efermi, sume, entk,rhos
      REAL(DP)  eold, timerd, timeorto
      REAL(DP)  fccc
      REAL(DP), ALLOCATABLE :: lambda(:,:)
      COMPLEX(DP), ALLOCATABLE :: clambda(:,:,:)
      COMPLEX(DP), ALLOCATABLE :: c0rot(:,:)

      REAL(DP), EXTERNAL :: cclock


! ... end of declarations
!  ----------------------------------------------

      nspin       = SIZE(fi, 2)
      istate      = 0
      eold        = 1.0d10  ! a large number
      isteep      = 0
      timerd      = 0
      timeorto    = 0
      svar1       = 2.d0
      svar2       = -1.d0
      svar3_0     = delt * delt / emass
      fccc        = 1.0d0
      isteep      = nreset

      ttreset_diis = .TRUE.

      s1 = cclock()
      old_clock_value = s1

      IF( ionode ) THEN
        WRITE( stdout,'(/,12X,"DIIS Optimizations for electron, starting ...")' )
        WRITE( stdout,'(  12X,"iter     erho          derho       xxxxx      seconds")' )
      END IF

! ... DO WHILE .NOT.doions

      DIIS_LOOP: DO nfi_l = 1, maxnstep

        ddoions      = .FALSE.

! ...   check for exit
        IF (check_stop_now()) THEN
          cm = c0
          EXIT DIIS_LOOP
        END IF

        CALL kspotential( 1, .FALSE., tforce, tstress, rhoe, &
          atoms, bec, becdr, eigr, ei1, ei2, ei3, sfac, c0, cdesc, tcel, ht0, fi, vpot, edft, timepre )

        s0 = cclock()
        seconds_per_iter = (s0 - old_clock_value)
        old_clock_value = s0

        IF( ionode ) THEN
            WRITE( stdout,113) nfi_l, edft%etot, edft%etot-eold, 0.d0, seconds_per_iter
113         FORMAT(10X,I5,2X,F14.6,2X,3D12.4)
        END IF

        IF (isteep .LT. nreset) THEN

          isteep = isteep + 1
          cm = c0
          CALL runcp(.FALSE., tortho, tsde, cm, c0, cgrad, cdesc, vpot, eigr, fi, &
            ekinc, timerd, timeorto, ht0, ei, bec, fccc )
          CALL update_wave_functions(cm, c0, cgrad, cdesc)

        ELSE

          SPIN: DO ispin = 1, nspin

! ...       compute local number of rows
            nx  = cdesc%nbt( ispin )
            nrl = nx / nproc_image
            IF( me_image .LT. MOD(nx, nproc_image) ) THEN
              nrl = nrl + 1
            END IF

! ...       initialize lambda to the identity matrix
            ALLOCATE(lambda(nrl, nx))
            CALL set_lambda(lambda)

! ...       distribute lambda's rows across processors with a blocking factor
! ...       of 1, ( row 1 to PE 1, row 2 to PE 2, .. row nproc_image+1 to PE 1 and
! ...       so on).

            CALL dforce_all( ispin, c0(:,:,ispin), fi(:,ispin), cgrad(:,:,ispin), &
                             vpot(:,ispin), eigr, bec, nupdwn, iupdwn )

            CALL proj( ispin, cgrad(:,:,ispin), cdesc, c0(:,:,ispin), cdesc, lambda)

            s4 = cclock()
            CALL simupd(ekinc(ispin), ddoions(ispin), c0(:,:,ispin), cgrad(:,:,ispin), cdesc, &
                svar1, svar2, svar3_0, edft%etot, fi(:,ispin), eigr, sfac, &
                vps, ttreset_diis(ispin), istate, cnorm, &
                eold, ndiis, nowv)
            CALL gram( vkb, bec, nkb, c0(1,1,ispin), SIZE(c0,1), cdesc%nbt( ispin ) )
            DEALLOCATE(lambda)
          END DO SPIN

          IF( ANY( ttreset_diis ) ) THEN
            isteep = 0
            ttreset_diis = .TRUE.
          END IF

        END IF

        IF ( ALL( ddoions ) ) THEN
          doions = .TRUE.
          EXIT DIIS_LOOP
        END IF

      END DO DIIS_LOOP

      s2 = cclock()
        

      IF ( doions ) THEN
        IF(ionode) WRITE( stdout,fmt="(12X,'runsdiis: convergence achieved successfully, in ',F8.2,' sec.')") (s2-s1)
      ELSE
        IF(ionode) WRITE( stdout,fmt="(12X,'runsdiis: convergence not achieved, in ',F8.2,' sec.')") (s2-s1)
      END IF

      IF ( tprint ) THEN
        CALL diis_eigs(.TRUE., atoms, c0, cdesc, fi, vpot, cgrad, eigr, bec)
      END IF

      cgrad = c0

      RETURN
      END SUBROUTINE runsdiis

!  ----------------------------------------------

      SUBROUTINE set_lambda_r( lambda )
! ...   initialize lambda to the identity matrix
        USE kinds
        USE mp_global, ONLY: me_image, nproc_image
        REAL(DP) :: lambda(:,:)
        INTEGER ib, ibl, nrl
          nrl = SIZE(lambda, 1)
          lambda = 0.d0
          ib = me_image + 1
          DO ibl = 1, nrl
            lambda(ibl,ib) = 1.d0  ! diagonal elements
            ib = ib + nproc_image
          END DO
        RETURN
      END SUBROUTINE set_lambda_r

!  ----------------------------------------------

      SUBROUTINE diis_eigs(tortho, atoms, c, cdesc, fi, vpot, eforce, eigr, bec )

        USE kinds
        USE wave_types
        USE wave_functions, ONLY: crot
        USE wave_constrains, ONLY: update_lambda
        USE wave_base, ONLY: dotp
        USE cell_base, ONLY: tpiba2
        USE electrons_module, ONLY: eigs, ei, nb_l, ib_owner, ib_local
        USE electrons_base, ONLY: iupdwn, nupdwn
        USE forces, ONLY: dforce_all
        USE orthogonalize
        USE pseudopotential,       ONLY: nspnl
        USE mp,                    ONLY: mp_sum
        USE mp_global,             ONLY: me_image, nproc_image, intra_image_comm
        USE atoms_type_module,     ONLY: atoms_type
        USE reciprocal_vectors,    ONLY: g, gx
        USE control_flags,         ONLY: gamma_only

        IMPLICIT NONE

! ...   ARGUMENTS
        COMPLEX(DP), INTENT(inout) ::  c(:,:,:)
        COMPLEX(DP), INTENT(inout) ::  eforce(:,:,:)
        TYPE (wave_descriptor), INTENT(in) :: cdesc
        REAL (DP), INTENT(in) ::  vpot(:,:), fi(:,:)
        REAL (DP) ::  bec(:,:)
        LOGICAL, INTENT(IN) :: TORTHO
        COMPLEX(DP) :: eigr(:,:)
        TYPE(atoms_type), INTENT(INOUT)  :: atoms     

! ...   LOCALS
        INTEGER     kk, i, k, j, iopt, iter, nwh
        INTEGER     ngw, ngw_g, n_occ, n, n_l, nspin, ispin
        INTEGER     ig, iprinte, nrl, jl, ibl
        LOGICAL     gamma_symmetry, gzero

        REAL(DP),    ALLOCATABLE :: gam(:,:)
        COMPLEX(DP), ALLOCATABLE :: cgam(:,:)

! ...   SUBROUTINE BODY
!
        nspin = cdesc%nspin
        ngw   = cdesc%ngwl
        ngw_g = cdesc%ngwt
        n     = cdesc%nbl( 1 )
        gzero = cdesc%gzero
        gamma_symmetry = cdesc%gamma
        n_l  = nb_l(1)

        ALLOCATE(gam(n_l, n), cgam(n_l, n))

! ...   electronic state diagonalization ==
        DO ispin = 1, nspin

          CALL nlsm1( n, 1, nspnl, eigr, c(1,1,ispin), bec )

! ...     Calculate | dH / dpsi(j) >
          CALL dforce_all( ispin, c(:,:,ispin), fi(:,ispin), eforce(:,:,ispin), &
                           vpot(:,ispin), eigr, bec, nupdwn, iupdwn )

! ...       Calculate Eij = < psi(i) | H | psi(j) > = < psi(i) | dH / dpsi(j) >
            DO i = 1, n
              IF( gamma_symmetry ) THEN
                CALL update_lambda( i,  gam, c( :, :, ispin), cdesc, eforce( :, i, ispin) ) 
              ELSE
                CALL update_lambda( i, cgam, c( :, :, ispin), cdesc, eforce( :, i, ispin) ) 
              END IF
            END DO

            CALL eigs( n, gam, tortho, fi(:,ispin), ei(:,ispin) )
          END DO

        DEALLOCATE(gam, cgam)

        RETURN

      END SUBROUTINE diis_eigs



      END MODULE rundiis_module
