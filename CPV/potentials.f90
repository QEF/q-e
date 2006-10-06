!
! Copyright (C) 2002-2005 FPMD-CPV groups
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!  AB INITIO COSTANT PRESSURE MOLECULAR DYNAMICS
!  ----------------------------------------------
!  Car-Parrinello Parallel Program

#include "f_defs.h"



      SUBROUTINE potential_print_info( iunit )

        USE control_flags, ONLY: iesr

        INTEGER, INTENT(IN) :: iunit

        WRITE(iunit,50)
        WRITE(iunit,115) (2*iesr+1),(2*iesr+1),(2*iesr+1)

   50   FORMAT(//,3X,'Potentials Parameters',/,3X,'---------------------')
  115   FORMAT(   3X,'Ewald sum over ',I1,'*',I1,'*',I1,' cells')

        RETURN
      END SUBROUTINE potential_print_info
     

      SUBROUTINE vofmean_x( sfac, rhops, rhoeg )

        USE kinds,              ONLY: DP
        USE control_flags,      ONLY: vhrmin, vhrmax, vhasse
        USE cp_main_variables,  ONLY: nfi
        USE constants,          ONLY: fpi
        USE cell_base,          ONLY: tpiba2, tpiba
        USE mp,                 ONLY: mp_sum
        USE mp_global,          ONLY: nproc_image, me_image, intra_image_comm
        USE io_global,          ONLY: ionode
        USE io_files,           ONLY: opt_unit
        USE gvecp,              ONLY: ngm
        USE reciprocal_vectors, ONLY: gstart, gx, g

        IMPLICIT NONE

        REAL(DP),    INTENT(IN)   :: RHOPS(:,:)
        COMPLEX(DP), INTENT(IN)   :: RHOEG(:)
        COMPLEX(DP), INTENT(IN)   :: sfac(:,:)

        COMPLEX(DP) :: fpi_tpiba2, rp, vcg
        REAL(DP)    :: gxt, dr, r
        REAL(DP), ALLOCATABLE :: vrmean(:) 
        
        INTEGER :: ig, is, iasse, ipiano1, ipiano2
        INTEGER :: ir
        INTEGER :: vhnr = 640

        IF( (vhasse.NE.'X') .AND. (vhasse.NE.'Y') .AND. (vhasse.NE.'Z') ) THEN
          CALL errore( ' vofmean ', ' wrong asse ',0)
        END IF 
        IF( vhrmax .LE. vhrmin ) THEN
          CALL errore( ' vofmean ', ' wrong rmax or rmin ',0)
        END IF 
        IF( vhnr .LE. 0 ) THEN
          CALL errore( ' vofmean ', ' wrong nr ',0)
        END IF 

        fpi_tpiba2 = CMPLX(FPI/TPIBA2,0.0d0)
        dr = ( vhrmax - vhrmin ) / vhnr

        ALLOCATE(vrmean(vhnr))

        IF( vhasse.EQ.'X' ) THEN
          iasse = 1
          ipiano1 = 2
          ipiano2 = 3
        ELSE IF(vhasse.EQ.'Y') THEN
          iasse = 2
          ipiano1 = 1
          ipiano2 = 3
        ELSE
          iasse = 3
          ipiano1 = 1
          ipiano2 = 2
        END IF

        DO ir = 1, vhnr
          vrmean(ir) = 0.0d0
          r = vhrmin + (ir-1)*dr
          DO ig = gstart, ngm
            rp   = (0.D0,0.D0)
            DO is = 1, SIZE( sfac, 2 )
              rp = rp + sfac( ig, is ) * rhops( ig, is )
            END DO
            IF((gx(ipiano1,IG).EQ.0.d0).AND. &
               (gx(ipiano2,IG).EQ.0.d0))THEN
              vcg       = fpi_tpiba2 * (rhoeg(ig) + rp) / g(ig)
              gxt       = gx(iasse, ig) * tpiba
              vrmean(ir) = vrmean(ir) + DBLE(vcg)  * COS(gxt*r) 
              vrmean(ir) = vrmean(ir) - AIMAG(vcg) * SIN(gxt*r)
            END IF
          END DO
          vrmean(ir) = 2.0d0 * vrmean(ir)
        END DO
        CALL mp_sum( vrmean, intra_image_comm )


        IF(ionode) THEN
          OPEN( unit = opt_unit, file = 'Vh_mean.out', position = 'append' )
          WRITE(opt_unit,*) nfi
          DO ir = 1, vhnr
            r = vhrmin + (ir-1)*dr
            WRITE(opt_unit,100) r, vrmean(ir)
          END DO
          CLOSE( unit = opt_unit )
 100      FORMAT(2X,F14.8,F14.8)
        END IF
        
        DEALLOCATE(vrmean)

        RETURN
      END SUBROUTINE vofmean_x

!  -------------------------------------------------------------------------

      SUBROUTINE kspotential_x &
        ( nfi, tprint, tforce, tstress, rhoe, atoms, bec, becdr, eigr, &
          ei1, ei2, ei3, sfac, c0, cdesc, tcel, ht, fi, vpot, edft )

        USE kinds,             ONLY: DP
        USE cp_interfaces,     ONLY: rhoofr, nlrh, vofrhos
        USE energies,          ONLY: dft_energy_type
        USE cell_base,         ONLY: boxdimensions
        USE atoms_type_module, ONLY: atoms_type
        USE wave_types,        ONLY: wave_descriptor
        USE dener,             ONLY: denl6, dekin6

        IMPLICIT NONE

! ...   declare subroutine arguments
        INTEGER,              INTENT(IN)    :: nfi
        LOGICAL, INTENT(IN) :: tforce, tstress, tprint
        REAL(DP) :: rhoe(:,:)
        TYPE (atoms_type),    INTENT(INOUT) :: atoms
        REAL(DP) :: bec(:,:)
        REAL(DP) :: becdr(:,:,:)
        COMPLEX(DP) :: eigr(:,:)
        COMPLEX(DP) :: ei1(:,:)
        COMPLEX(DP) :: ei2(:,:)
        COMPLEX(DP) :: ei3(:,:)
        COMPLEX(DP), INTENT(IN) :: sfac(:,:)
        COMPLEX(DP),         INTENT(INOUT) :: c0(:,:)
        TYPE (wave_descriptor),  INTENT(IN) :: cdesc
        LOGICAL   :: tcel
        TYPE (boxdimensions), INTENT(INOUT) ::  ht
        REAL(DP), INTENT(IN) :: fi(:)
        REAL(DP)    :: vpot(:,:)
        TYPE (dft_energy_type) :: edft

        CALL nlrh( c0, tforce, tstress, atoms%for, bec, becdr, eigr, edft%enl, denl6 )

        CALL rhoofr( nfi, tstress, c0, fi, rhoe, ht%deth, edft%ekin, dekin6 )

        CALL vofrhos( tprint, tforce, tstress, rhoe, atoms, vpot, bec, &
                      c0, cdesc, fi, eigr, ei1, ei2, ei3, sfac, &
                      ht, edft )

        RETURN
      END SUBROUTINE kspotential_x

!=----------------------------------------------------------------------------=!

   SUBROUTINE vofrhos_x &
      ( tprint, tforce, tstress, rhoe, atoms, vpot, bec, c0, cdesc, fi, &
        eigr, ei1, ei2, ei3, sfac, box, edft )

      !  this routine computes:
      !  ekin = dft_kinetic term of the DFT functional (see dft_kinetic_energy)
      !  the one-particle potential v in real space,
      !  the total energy etot,
      !  the forces fion acting on the ions,
      !  the matrix ngdrt used to compute the pair-correlation
      !  function gdr for the ions
      !
      !  Moreover, this routin is re-written also to calculate the self-interaction-correction
      !  has proposed by Mauri et al. (PRB 2005), taking also into account the 'comment' 
      !  proposed by Sprik et al. (ICR 2005). 
      !  Thus, we introduce the parameters sic_alpha and sic_epsilon to correct the 
      !  the exchange-correlation and the electronic hartree potentials, respectively.
      !  They are two empirical parameters, thus to remain in a ab-initio
      !  set them equal to 1.0d0.
      !  Sprik et al. showed that, in same cases, i.e. OH radical, it should be better
      !  to under estimate the correction to ex-ch, since in same way the exch already
      !  corrects the electronic hartree part.

!fran:  My personal considerations: 
  !     the SIC is a way to correct the self-interaction
  !     of ONE and only ONE e- that lives in an unpaired electronic level
  !     we have choosen for it the spin up
  !     the other e- are fictitious calculate in a LSD approach
  !     even if they feel a different potential
  !     we constrain them to have the same force, and the same eigenvalues, the same eigenstate
  !     When you applied this SIC scheme to a molecule or to an atom, which are neutral,
  !     remeber hat you have to consider another correction to the energy level as proposed
  !     by Landau: infact if you start from a neutral system and subtract the self-intereaction
  !     the unpaired e- feels a charge system. Thus remeber a correction term ~2.317(Madelung)/2L_box
 

      ! ... include modules

      USE kinds,          ONLY: DP
      USE control_flags,  ONLY: tscreen, iprsta, iesr, tvhmean
      USE mp_global,      ONLY: nproc_image, me_image, intra_image_comm
      USE mp,             ONLY: mp_sum
      USE cell_base,      ONLY: tpiba2, boxdimensions
      USE ions_base,      ONLY: rcmax, zv, nsp
      USE fft_base,       ONLY: dfftp
      USE energies,       ONLY: total_energy, dft_energy_type, ekin
      USE cp_interfaces,  ONLY: pstress, stress_kin, compute_gagb, stress_nl, &
                                stress_local, add_drhoph, stress_hartree
      USE stress_param,   ONLY: dalbe
      USE funct,          ONLY: dft_is_gradient
      USE vanderwaals,    ONLY: tvdw, vdw
      USE wave_types,     ONLY: wave_descriptor
      USE io_global,      ONLY: ionode, stdout
      USE sic_module,     ONLY: self_interaction, sic_epsilon, sic_alpha !!TO ADD!!!
      USE gvecp,          ONLY: ngm
      USE local_pseudo,   ONLY: vps, rhops
      USE atom,           ONLY: nlcc
      USE core,           ONLY: nlcc_any, rhocg, drhocg
      USE cp_interfaces,  ONLY: fwfft, invfft, add_core_charge, core_charge_forces
      USE electrons_base, ONLY: iupdwn, nupdwn, nspin
      !
      USE reciprocal_vectors, ONLY: gx, g, gstart
      USE atoms_type_module,  ONLY: atoms_type
      USE cp_interfaces,      ONLY: exch_corr_energy, stress_xc, vofmean
      USE cp_interfaces,      ONLY: vofloc, vofps, self_vofhar, force_loc, fillgrad
      use grid_dimensions,    only: nr1, nr2, nr3, nnrx
      use dener,              only: dekin6, denl6

      IMPLICIT NONE

! ... declare subroutine arguments
      LOGICAL, INTENT(IN) :: tprint, tforce, tstress
      REAL(DP)            :: vpot(:,:)
      REAL(DP),    INTENT(IN) :: fi(:)
      REAL(DP)    :: bec(:,:)
      COMPLEX(DP) :: ei1(:,:)
      COMPLEX(DP) :: ei2(:,:)
      COMPLEX(DP) :: ei3(:,:)
      COMPLEX(DP) :: eigr(:,:)
      COMPLEX(DP),   INTENT(IN) :: c0(:,:)
      TYPE (atoms_type), INTENT(INOUT) :: atoms
      TYPE (wave_descriptor), INTENT(IN) :: cdesc
      TYPE (boxdimensions),    INTENT(INOUT) :: box
      TYPE (dft_energy_type) :: edft
      REAL(DP) :: rhoe(:,:)
      COMPLEX(DP), INTENT(IN) :: sfac(:,:)

      TYPE (dft_energy_type) :: edft_self

! ... declare functions
      REAL(DP)  DDOT

! ... declare other variables

      COMPLEX(DP), ALLOCATABLE :: vloc(:), self_vloc(:)
      COMPLEX(DP), ALLOCATABLE :: rhog(:), drhog(:,:), rhoeg(:,:)
      COMPLEX(DP), ALLOCATABLE :: psi(:)
      COMPLEX(DP), ALLOCATABLE :: screen_coul(:)
      !
      REAL(DP),    ALLOCATABLE :: rhoetr(:,:)
      REAL(DP),    ALLOCATABLE :: fion_vdw(:,:)
      REAL(DP),    ALLOCATABLE :: grho(:,:,:)
      REAL(DP),    ALLOCATABLE :: v2xc(:,:,:)
      REAL(DP),    ALLOCATABLE :: fion(:,:)

      REAL(DP),    ALLOCATABLE :: self_rho(:,:)
      REAL(DP),    ALLOCATABLE :: self_vpot(:,:)
      REAL(DP),    ALLOCATABLE :: self_grho(:,:,:)
      REAL(DP),    ALLOCATABLE :: self_v2xc(:,:,:)

      REAL(DP),    ALLOCATABLE :: gagb(:,:)

      COMPLEX(DP) :: ehtep
      REAL(DP)    :: self_exc, self_vxc

      REAL(DP)  :: summing1, summing2

      COMPLEX(DP) :: ehp, eps

      REAL(DP)  :: dum, exc, vxc, ehr, strvxc
      REAL(DP)  :: omega, desr(6), pesum(16)
      REAL(DP), DIMENSION (6) :: deht, deps, dexc, dvdw

      LOGICAL :: ttscreen, ttsic, tgc

      INTEGER ig1, ig2, ig3, is, ia, ig, isc, iflag, iss
      INTEGER ik, i, j, k, isa, idum
      INTEGER :: ierr

      DATA iflag / 0 /
      SAVE iflag, desr

      !  end of declarations
      !  ----------------------------------------------

      CALL start_clock( 'vofrho' )

      edft%evdw = 0.0d0
      !
      ttscreen = .FALSE.   ! .TRUE. to enable cluster boundary conditions
      !
      ttsic    = ( ABS(self_interaction) /= 0  ) 

      omega    = box%deth
      !
      tgc      = dft_is_gradient()
      !
      !   Allocate local array
      !
      IF( tstress ) THEN
         !
         ALLOCATE( gagb( 6, ngm ) )
         ALLOCATE( drhog( ngm, 6 ) )
         !
         CALL compute_gagb( gagb, gx, ngm, tpiba2 )
         !
      END IF
      
      ALLOCATE( fion( 3, atoms%nat ) )

      fion = atoms%for( 1:3, 1:atoms%nat )
      !

      ALLOCATE( rhoeg ( ngm, nspin ) )
      ALLOCATE( rhog( ngm ) )
      ALLOCATE( vloc( ngm ) )
      ALLOCATE( psi( SIZE( rhoe, 1 ) ) )
      !
      IF( tscreen ) THEN
         !
         IF( tprint .AND. ionode ) THEN
            WRITE( stdout,fmt="(3X,'Using screened Coulomb potential for cluster calculation')")
         END IF
         !
         ALLOCATE( screen_coul( ngm ) )
         !
         CALL cluster_bc( screen_coul, g, box%deth, box%hmat )
         !
      END IF

      !
      !
      IF( tstress .OR. tforce .OR. iflag == 0 )  THEN
         !
         CALL vofesr( iesr, edft%esr, desr, fion, atoms%taus, tstress, box%hmat )
         !
         IF( iflag == 0 ) &
            WRITE( stdout, fmt="(/,3X,'ESR (real part of Ewald sum) = ',D16.8,/)" ) edft%esr
         iflag = 1
         !
      END IF

      ! ... Van Der Waals energy and forces
      !
      IF ( tvdw ) THEN
         CALL VdW( edft%evdw, atoms, fion, box )
         IF( tstress ) THEN
            ! CALL vdw_stress(c6, iesr, stau0, dvdw, na, nax, nsp)
         END IF
      END IF


      ! ... FFT: rho(r) --> rho(g)  
      !
      DO iss = 1, nspin

         psi = rhoe(:,iss)

         CALL fwfft(   'Dense', psi, dfftp%nr1, dfftp%nr2, dfftp%nr3, dfftp%nr1x, dfftp%nr2x, dfftp%nr3x )
         CALL psi2rho( 'Dense', psi, dfftp%nnr, rhoeg(:,iss), ngm )
 
      END DO

      rhog( 1:ngm ) = rhoeg( 1:ngm, 1 )
      IF( nspin > 1 ) rhog( 1:ngm ) = rhog( 1:ngm ) + rhoeg( 1:ngm, 2 )

      IF( tstress ) THEN
         !
         ! add drho_e / dh
         !
         DO k = 1, 6
            drhog( 1:ngm, k ) = - rhog( 1:ngm ) * dalbe( k )
         END DO
         !
      END IF

      ! ... Calculate local part of the pseudopotential and its energy contribution (eps)

      CALL vofps( eps, vloc, rhog, vps, sfac, box%deth )

      edft%epseu = DBLE(eps)

      IF( tstress ) THEN
         !
         CALL stress_local( deps, edft%epseu, gagb, sfac, rhog, drhog, box%deth )
         !
      END IF

      ! ... Calculate hartree potential and energy (eh)
      !
      CALL vofloc( tscreen, edft%ehte, edft%ehti, ehp, vloc, rhog, &
                   rhops, vps, sfac, box%deth, screen_coul )

      IF( tforce ) THEN
         !
         CALL force_loc( tscreen, rhog, fion, rhops, vps, eigr, &
                         ei1, ei2, ei3, sfac, box%deth, screen_coul )
         !
      END IF

      edft%self_ehte = 0.d0

      IF( ttsic ) THEN
         !
         ! ... Calculate Self-interaction correction --- Hartree part
         !
         ALLOCATE ( self_vloc( ngm ) )
         CALL self_vofhar( ttscreen, edft%self_ehte, self_vloc, rhoeg, omega, box%hmat )
         !
      END IF 

      edft%eh  = DBLE( ehp ) - edft%self_ehte
      edft%eht = edft%eh + edft%esr - edft%eself 
     
      IF( tprint .AND. tvhmean ) THEN
         !
         CALL vofmean( sfac, rhops, rhog )
         !
      END IF

      IF( tstress ) THEN
         !
         ! add Ionic pseudo charges  rho_I
         !
         DO is = 1, nsp
            DO ig = gstart, ngm
               rhog( ig ) = rhog( ig ) + sfac( ig, is ) * rhops( ig, is )
            END DO
         END DO
         !
         ! add drho_I / dh
         !
         CALL add_drhoph( drhog, sfac, gagb )
         !
         CALL stress_hartree( deht, edft%eh, sfac, rhog, drhog, gagb, box%deth )

         DEALLOCATE( drhog )
         !
      END IF

      DEALLOCATE( rhog )

      ALLOCATE( rhoetr( nnrx, nspin ) )
      !
      rhoetr = 0.0d0

      DO iss = 1, nspin

         ! ... add core contribution to the charge

         CALL DCOPY( nnrx, rhoe(1,iss), 1, rhoetr(1,iss), 1 )

         IF( nlcc_any ) THEN

            ! ...     add core correction:  rhoeg = rhoeg + cc;  rhoetr = rhoe  + cc

            CALL add_core_charge( rhoeg(:,iss), rhoetr(:,iss), sfac, rhocg, atoms%nsp )

         END IF

      END DO

      !
      !
      !  exchange and correlation potential
      !
      !

      IF(tgc) THEN
         ALLOCATE( grho( nnrx, 3, nspin ) )
         ALLOCATE( v2xc( nnrx, nspin, nspin ) )
      ELSE
         ALLOCATE( grho( 1, 1, 1 ) )
         ALLOCATE( v2xc( 1, 1, 1 ) )
      END IF

      grho = 0.0d0
      v2xc = 0.0d0

      IF( tgc ) THEN
         !
         CALL fillgrad( nspin, rhoeg, grho )
         !
      END IF

      CALL exch_corr_energy( rhoetr, grho, vpot, exc, vxc, v2xc )

      self_exc  = 0.d0
      self_vxc  = 0.d0
      !
      IF ( ttsic ) THEN                

         ALLOCATE (self_rho( nnrx, 2), STAT = ierr)
         IF( ierr /= 0 ) CALL errore(' vofrhos ', ' allocating self_rho ', ierr)
         !
         self_rho(:,1) = rhoetr(:,2)
         self_rho(:,2) = rhoetr(:,2)

         IF ( tgc ) THEN
            !
            ALLOCATE(self_grho( nnrx, 3, nspin ), STAT = ierr)
            IF( ierr /= 0 ) CALL errore(' vofrhos ', ' allocating self_grho ', ierr)
            !
            self_grho(:,:,1) = grho(:,:,2)
            self_grho(:,:,2) = grho(:,:,2)
            !
            ALLOCATE(self_v2xc( nnrx, nspin, nspin ), STAT = ierr)
            IF( ierr /= 0 ) CALL errore(' vofrhos ', ' allocating self_v2xc ', ierr)
            !
         ENDIF

         ALLOCATE ( self_vpot( nnrx, 2 ), STAT = ierr )
         IF( ierr /= 0 ) CALL errore(' vofrhos ', ' allocating self_vpot ', ierr)
   
         self_vpot  = 0.D0

         CALL exch_corr_energy( self_rho, self_grho, self_vpot, self_exc, self_vxc, self_v2xc )

         vpot (:,1) = ( 1.0d0 - sic_alpha ) * vpot(:,1)
         vpot (:,2) = ( 1.0d0 - sic_alpha ) * vpot(:,2) + sic_alpha * ( self_vpot(:,2) + self_vpot(:,1) )

         IF (tgc) THEN
           !
           v2xc(:,1,1) = ( 1.0d0 - sic_alpha ) * v2xc(:,1,1)
           v2xc(:,2,2) = ( 1.0d0 - sic_alpha ) * v2xc(:,2,2) + sic_alpha * ( self_v2xc(:,2,2) + self_v2xc(:,1,1) )
           !
         END IF

         self_exc = sic_alpha * ( exc - self_exc )
         !
         exc      = exc - self_exc
         !
         self_vxc = sic_alpha * ( vxc - self_vxc )
         !
         vxc      = vxc - self_vxc
         !
      END IF  

      IF ( tstress ) THEN
         !
         strvxc = ( exc - vxc ) * omega / DBLE( nr1 * nr2 * nr3 )
         ! 
      END IF

      edft%exc       = exc      * omega / DBLE( nr1 * nr2 * nr3 )
      edft%vxc       = vxc      * omega / DBLE( nr1 * nr2 * nr3 )
      edft%self_exc  = self_exc * omega / DBLE( nr1 * nr2 * nr3 )
      edft%self_vxc  = self_vxc * omega / DBLE( nr1 * nr2 * nr3 )

      CALL mp_sum( edft%vxc,      intra_image_comm )
      CALL mp_sum( edft%exc,      intra_image_comm )
      CALL mp_sum( edft%self_exc, intra_image_comm )
      CALL mp_sum( edft%self_vxc, intra_image_comm )


      IF( nlcc_any ) THEN
        !
        ! ...   xc potential (vpot) from real to G space, to compute nlcc forces
        ! ...   rhoeg = fwfft(vpot)
        !
        DO iss = 1, nspin
           !
           psi = vpot(:,iss)
           !
           CALL fwfft(   'Dense', psi, dfftp%nr1, dfftp%nr2, dfftp%nr3, dfftp%nr1x, dfftp%nr2x, dfftp%nr3x )
           CALL psi2rho( 'Dense', psi, dfftp%nnr, rhoeg(:,iss), ngm )
           ! 
        END DO
        !
        ! ...   now rhoeg contains the xc potential
        !
        IF (tforce) THEN
           !
           CALL core_charge_forces( fion, rhoeg, rhocg, nlcc, atoms, box, ei1, ei2, ei3 )
           !
        END IF
        !
      END IF

      !
      ! ... vloc(g): hartree and local part of the pseudo potentials (in
      ! ...          reciprocal space)
      !

      IF ( ttsic ) THEN

        CALL rho2psi( 'Dense', psi, dfftp%nnr, self_vloc, ngm )
        CALL invfft(  'Dense', psi, dfftp%nr1, dfftp%nr2, dfftp%nr3, dfftp%nr1x, dfftp%nr2x, dfftp%nr3x )
        !
        vpot(:,1) =  vpot(:,1) - DBLE( psi(:) )
        vpot(:,2) =  vpot(:,2) + DBLE( psi(:) )

      END IF

      ! ...   add hartree end local pseudo potentials ( invfft(vloc) )
      ! ...   to xc potential (vpot).
      ! ...   vpot = vpot + invfft(vloc)

      CALL rho2psi( 'Dense', psi, dfftp%nnr, vloc, ngm )
      CALL invfft(  'Dense', psi, dfftp%nr1, dfftp%nr2, dfftp%nr3, dfftp%nr1x, dfftp%nr2x, dfftp%nr3x )

      ! ... now potentials are in real space
      ! ... vpot(r) = hartree + xc + local

      DO iss = 1, nspin
         vpot(:,iss) = vpot(:,iss) + DBLE( psi )
      END DO

      IF( ttsic ) THEN
         IF( tgc ) THEN
            DEALLOCATE( self_grho  )
            DEALLOCATE( self_v2xc  )
         END IF
         DEALLOCATE( self_vpot )
         DEALLOCATE( self_rho  )
         DEALLOCATE( self_vloc )
      END IF

      IF( tstress ) THEN
         !
         ! ... compute exchange & correlation energy contribution
         ! 
         CALL stress_xc( dexc, strvxc, sfac, rhoeg, grho, v2xc, gagb, nlcc, drhocg, box )
         !
      END IF

      ! ... sum up forces
      !
      IF (tforce) THEN
         CALL mp_sum(fion, intra_image_comm)
      END IF

      ! ... sum up energy contributions
      !
      CALL total_energy( edft )


      ! ... sum up stress tensor
      !
      IF( tstress ) THEN
         CALL pstress( box%pail, desr, dekin6, denl6, deps, deht, dexc, box )
      END IF


      ! ... Copy new atomic forces on for type member
      !
      atoms%for( 1:3, 1:atoms%nat ) = fion

      DEALLOCATE( rhoeg, rhoetr, grho, v2xc, fion )
      DEALLOCATE( vloc, psi )

      !
      IF( tscreen ) THEN
         DEALLOCATE( screen_coul )
      END IF
      !
      IF( tstress ) THEN
         DEALLOCATE( gagb )
      END IF

      CALL stop_clock( 'vofrho' )

      ! ... Flush stdout

      CALL flush_unit( stdout )

      RETURN
      END SUBROUTINE vofrhos_x

!=----------------------------------------------------------------------------=!

  SUBROUTINE cluster_bc( screen_coul, hg, omega, hmat )

      USE kinds,           ONLY: DP
      USE green_functions, ONLY: greenf
      USE mp_global,       ONLY: me_image
      USE fft_base,        ONLY: dfftp
      USE cp_interfaces,   ONLY: fwfft
      USE gvecp,           ONLY: ngm
      USE constants,       ONLY: gsmall, pi
      USE cell_base,       ONLY: tpiba2, s_to_r, alat
      use grid_dimensions, only: nr1, nr2, nr3, nr1l, nr2l, nr3l, nnrx

      IMPLICIT NONE
      
      REAL(DP), INTENT(IN) :: hg( ngm )
      REAL(DP), INTENT(IN) :: omega, hmat( 3, 3 )
      COMPLEX(DP) :: screen_coul( ngm )

      ! ... declare external function
      !
      REAL(DP) :: erf, erfc
      EXTERNAL erf, erfc

      ! ... Locals
      !
      COMPLEX(DP), ALLOCATABLE :: grr(:)
      COMPLEX(DP), ALLOCATABLE :: grg(:)
      REAL(DP) :: rc, r(3), s(3), rmod, g2, rc2, arg, fact
      INTEGER   :: ig, i, j, k, ir
      INTEGER   :: ir1, ir2, ir3

      ir1 = 1
      ir2 = 1
      ir3 = 1
      DO k = 1, me_image
        ir3 = ir3 + dfftp%npp( k )
      END DO

      ALLOCATE( grr( nnrx ) )
      ALLOCATE( grg( SIZE( screen_coul ) ) )

      grr = 0.0d0

      ! ... Martina and Tuckerman convergence criterium
      !
      rc  = 7.0d0 / alat
      rc2 = rc**2
      fact  = omega / ( nr1 * nr2 * nr3 )
      IF( MOD(nr1 * nr2 * nr3, 2) /= 0 ) fact = -fact

      DO k = 1, nr3l
        s(3) = DBLE ( (k-1) + (ir3 - 1) ) / nr3 - 0.5d0
        DO j = 1, nr2l
          s(2) = DBLE ( (j-1) + (ir2 - 1) ) / nr2 - 0.5d0
          DO i = 1, nr1l
            s(1) = DBLE ( (i-1) + (ir1 - 1) ) / nr1 - 0.5d0
            CALL S_TO_R( S, R, hmat )
            rmod = SQRT( r(1)**2 + r(2)**2 + r(3)**2 )
            ir =  i + (j-1)*dfftp%nr1x + (k-1)*dfftp%nr1x*dfftp%nr2x
            IF( rmod < gsmall ) THEN
              grr( ir ) = fact * 2.0d0 * rc / SQRT( pi )
            ELSE
              grr( ir ) = fact * erf( rc * rmod ) / rmod
            END IF
          END DO
        END DO
      END DO

      ! grg = FFT( grr )

      CALL fwfft(   'Dense', grr, dfftp%nr1, dfftp%nr2, dfftp%nr3, dfftp%nr1x, dfftp%nr2x, dfftp%nr3x )
      CALL psi2rho( 'Dense', grr, dfftp%nnr, grg, ngm )

      DO ig = 1, SIZE( screen_coul )
        IF( hg(ig) < gsmall ) THEN
          screen_coul(ig) = grg(1) - ( - pi / rc2 )
        ELSE
          g2  = tpiba2 * hg(ig)
          arg = - g2 / ( 4.0d0 * rc2 )
          screen_coul(ig) = grg(ig) - ( 4.0d0 * pi * EXP( arg ) / g2 ) 
        END IF
      END DO

      DEALLOCATE( grr, grg )

    RETURN
  END SUBROUTINE cluster_bc


!=----------------------------------------------------------------------------=!


   SUBROUTINE vofps_x( eps, vloc, rhoeg, vps, sfac, omega )

      !  this routine computes:
      !  omega = ht%deth
      !  vloc_ps(ig)  =  (sum over is) sfac(is,ig) * vps(ig,is)
      !
      !  Eps = Fact * omega * (sum over ig) cmplx ( rho_e(ig) ) * vloc_ps(ig)
      !  if Gamma symmetry Fact = 2 else Fact = 1
      !

      USE kinds,              ONLY: DP
      USE io_global,          ONLY: stdout
      USE ions_base,          ONLY: nsp
      USE gvecp,              ONLY: ngm
      USE reciprocal_vectors, ONLY: gstart, gx, g
      USE mp_global,          ONLY: intra_image_comm
      USE mp,                 ONLY: mp_sum

      IMPLICIT NONE

      ! ... Arguments

      REAL(DP),    INTENT(IN)  :: vps(:,:)
      REAL(DP),    INTENT(IN)  :: omega
      COMPLEX(DP), INTENT(OUT) :: vloc(:)
      COMPLEX(DP), INTENT(IN)  :: rhoeg(:)
      COMPLEX(DP), INTENT(IN)  :: sfac(:,:)
      COMPLEX(DP), INTENT(OUT) :: eps

      ! ... Locals

      INTEGER     :: is, ig
      COMPLEX(DP) :: vp

      ! ... Subroutine body ...
      !
      eps   = (0.D0,0.D0)
      !
      DO ig = gstart, ngm 

        vp   = (0.D0,0.D0)
        DO is = 1, nsp
          vp = vp + sfac( ig, is ) * vps( ig, is )
        END DO

        vloc(ig) = vp
        eps      = eps  +     vp * CONJG( rhoeg( ig ) )

      END DO
      ! ... 
      ! ... G = 0 element
      !
      IF ( gstart == 2 ) THEN
        vp = (0.D0,0.D0)
        DO is = 1, nsp
          vp = vp + sfac( 1, is) * vps(1, is)
        END DO
        vloc(1) = VP
        eps     = eps + vp * CONJG( rhoeg(1) ) * 0.5d0
      END IF
      !
      eps = 2.D0 * eps  * omega
      !
      CALL mp_sum( eps, intra_image_comm )

      RETURN
   END SUBROUTINE vofps_x


!=----------------------------------------------------------------------------=!

  SUBROUTINE vofloc_x( tscreen, ehte, ehti, eh, vloc, rhoeg, &
                     rhops, vps, sfac, omega, screen_coul )

      !  this routine computes:
      !  omega = ht%deth
      !  rho_e(ig)    =  (sum over iss) rhoeg(ig,iss) 
      !  rho_I(ig)    =  (sum over is) sfac(is,ig) * rhops(ig,is) 
      !  vloc_h(ig)   =  fpi / ( g(ig) * tpiba2 ) * { rho_e(ig) + rho_I(ig) }
      !
      !  Eh  = Fact * omega * (sum over ig) * fpi / ( g(ig) * tpiba2 ) *
      !        { rho_e(ig) + rho_I(ig) } * conjugate { rho_e(ig) + rho_I(ig) }
      !  if Gamma symmetry Fact = 1 else Fact = 1/2
      !
      !  Hatree potential and local pseudopotential
      !  vloc(ig)     =  vloc_h(ig) + vloc_ps(ig) 
      !

      USE kinds,              ONLY: DP
      USE constants,          ONLY: fpi
      USE cell_base,          ONLY: tpiba2, tpiba
      USE io_global,          ONLY: stdout
      USE reciprocal_vectors, ONLY: gstart, g
      USE ions_base,          ONLY: nsp
      USE gvecp,              ONLY: ngm
      USE mp_global,          ONLY: intra_image_comm
      USE mp,                 ONLY: mp_sum

      IMPLICIT NONE

      ! ... Arguments

      LOGICAL,     INTENT(IN)    :: tscreen
      REAL(DP),    INTENT(IN)    :: rhops(:,:), vps(:,:)
      COMPLEX(DP), INTENT(INOUT) :: vloc(:)
      COMPLEX(DP), INTENT(IN)    :: rhoeg(:)
      COMPLEX(DP), INTENT(IN)    :: sfac(:,:)
      REAL(DP),    INTENT(OUT)   :: ehte, ehti
      REAL(DP),    INTENT(IN)    :: omega
      COMPLEX(DP), INTENT(OUT)   :: eh
      COMPLEX(DP), INTENT(IN)    :: screen_coul(:)

      ! ... Locals

      INTEGER     :: is, ig
      REAL(DP)    :: fpibg, cost
      COMPLEX(DP) :: rhet, rhog, rp, vscreen

      ! ... Subroutine body ...

      eh    = 0.0d0
      ehte  = 0.0d0
      ehti  = 0.0d0

      DO ig = gstart, ngm 

        rp   = (0.D0,0.D0)
        DO is = 1, nsp
          rp = rp + sfac( ig, is ) * rhops( ig, is )
        END DO

        rhet  = rhoeg( ig )
        rhog  = rhet + rp

        IF( tscreen ) THEN
          fpibg     = fpi / ( g(ig) * tpiba2 ) + screen_coul(ig)
        ELSE
          fpibg     = fpi / ( g(ig) * tpiba2 )
        END IF

        vloc(ig) = vloc(ig)  +  fpibg *        rhog 
        eh       = eh        +  fpibg *        rhog * CONJG(rhog)
        ehte     = ehte      +  fpibg *   DBLE(rhet * CONJG(rhet))
        ehti     = ehti      +  fpibg *   DBLE(  rp * CONJG(rp))

      END DO
      ! ... 
      ! ... G = 0 element
      !
      IF ( gstart == 2 ) THEN
        rp = (0.D0,0.D0)
        IF( tscreen ) THEN
          vscreen = screen_coul(1)
        ELSE
          vscreen = 0.0d0
        END IF
        DO IS = 1, nsp
          rp = rp + sfac( 1, is) * rhops(1, is)
        END DO
        rhet    = rhoeg(1)
        rhog    = rhet + rp
        vloc(1) = vloc(1)   +  vscreen *   rhog
        eh      = eh        +  vscreen *        rhog * CONJG(rhog)
        ehte    = ehte      +  vscreen *   DBLE(rhet * CONJG(rhet))
        ehti    = ehti      +  vscreen *   DBLE(  rp * CONJG(rp))
      END IF
      ! ...
      eh   =        eh   * omega
      ehte =        ehte * omega
      ehti =        ehti * omega
      ! ...
      CALL mp_sum(eh  , intra_image_comm)
      CALL mp_sum(ehte, intra_image_comm)
      CALL mp_sum(ehti, intra_image_comm)
      !
      RETURN
  END SUBROUTINE vofloc_x


  SUBROUTINE force_loc_x( tscreen, rhoeg, fion, rhops, vps, eigr, ei1, ei2, ei3, &
                        sfac, omega, screen_coul )

      !  this routine computes:
      !
      !  Local contribution to the forces on the ions
      !  eigrx(ig,isa)   = ei1( mill(1,ig), isa)
      !  eigry(ig,isa)   = ei2( mill(2,ig), isa)
      !  eigrz(ig,isa)   = ei3( mill(3,ig), isa)
      !  fpibg           = fpi / ( g(ig) * tpiba2 )
      !  tx_h(ig,is)     = fpibg * rhops(ig, is) * CONJG( rho_e(ig) + rho_I(ig) )
      !  tx_ps(ig,is)    = vps(ig,is) * CONJG( rho_e(ig) )
      !  gx(ig)          = CMPLX(0.D0, gx(1,ig)) * tpiba
      !  fion(x,isa)     = fion(x,isa) + 
      !      Fact * omega * ( sum over ig, iss) (tx_h(ig,is) + tx_ps(ig,is)) * 
      !      gx(ig) * eigrx(ig,isa) * eigry(ig,isa) * eigrz(ig,isa) 
      !  if Gamma symmetry Fact = 2.0 else Fact = 1
      !

      USE kinds,              ONLY: DP
      USE constants,          ONLY: fpi
      USE cell_base,          ONLY: tpiba2, tpiba
      USE io_global,          ONLY: stdout
      USE grid_dimensions,    ONLY: nr1, nr2, nr3
      USE reciprocal_vectors, ONLY: mill_l, gstart, gx, g
      USE ions_base,          ONLY: nat, nsp, na
      USE gvecp,              ONLY: ngm

      IMPLICIT NONE

      ! ... Arguments

      LOGICAL     :: tscreen
      REAL(DP)    :: fion(:,:)
      REAL(DP)    :: rhops(:,:), vps(:,:)
      COMPLEX(DP) :: rhoeg(:)
      COMPLEX(DP), INTENT(IN) :: sfac(:,:)
      COMPLEX(DP) :: ei1(-nr1:nr1,nat)
      COMPLEX(DP) :: ei2(-nr2:nr2,nat)
      COMPLEX(DP) :: ei3(-nr3:nr3,nat)
      COMPLEX(DP) :: eigr(:,:)
      REAL(DP)    :: omega
      COMPLEX(DP) :: screen_coul(:)

      ! ... Locals

      INTEGER     :: is, ia, isa, ig, ig1, ig2, ig3
      REAL(DP)    :: fpibg
      COMPLEX(DP) :: cxc, rhet, rhog, vp, rp, gxc, gyc, gzc
      COMPLEX(DP) :: teigr, cnvg, cvn, tx, ty, tz
      COMPLEX(DP), ALLOCATABLE :: ftmp(:,:)

      ! ... Subroutine body ...

      ALLOCATE( ftmp( 3, SIZE( fion, 2 ) ) )
      
      ftmp = 0.0d0

      DO IG = gstart, ngm 

        RP   = (0.D0,0.D0)
        DO IS = 1, nsp
          RP = RP + sfac( ig, is ) * rhops( ig, is )
        END DO

        RHET  = RHOEG( ig )
        RHOG  = RHET + RP

        IF( tscreen ) THEN
          FPIBG     = fpi / ( g(ig) * tpiba2 ) + screen_coul(ig)
        ELSE
          FPIBG     = fpi / ( g(ig) * tpiba2 )
        END IF

        ig1  = mill_l(1,IG)
        ig2  = mill_l(2,IG)
        ig3  = mill_l(3,IG)
        GXC  = CMPLX(0.D0,gx(1,IG))
        GYC  = CMPLX(0.D0,gx(2,IG))
        GZC  = CMPLX(0.D0,gx(3,IG))
        isa = 1
        DO IS = 1, nsp
           CNVG  = RHOPS(IG,is) * FPIBG * CONJG(rhog)
           CVN   = VPS(ig, is)  * CONJG(rhet)
           TX = (CNVG+CVN) * GXC
           TY = (CNVG+CVN) * GYC
           TZ = (CNVG+CVN) * GZC
           DO IA = 1, na(is)
              TEIGR = ei1(IG1,ISA) * ei2(IG2,ISA) * ei3(IG3,ISA)
              ftmp(1,ISA) = ftmp(1,ISA) + TEIGR*TX
              ftmp(2,ISA) = ftmp(2,ISA) + TEIGR*TY
              ftmp(3,ISA) = ftmp(3,ISA) + TEIGR*TZ
              isa = isa + 1
           END DO
        END DO

      END DO
      !
      fion = fion + DBLE(ftmp) * 2.D0 * omega * tpiba

      DEALLOCATE( ftmp )
       
      RETURN
      END SUBROUTINE force_loc_x


!
!=----------------------------------------------------------------------------=!
   SUBROUTINE vofesr( iesr, esr, desr, fion, taus, tstress, hmat )
!=----------------------------------------------------------------------------=!

      USE kinds,       ONLY : DP
      USE constants,   ONLY : sqrtpm1
      USE cell_base,   ONLY : s_to_r, pbcs
      USE mp_global,   ONLY : nproc_image, me_image, intra_image_comm
      USE mp,          ONLY : mp_sum
      USE ions_base,   ONLY : rcmax, zv, nsp, na, nat
 
      IMPLICIT NONE

! ... ARGUMENTS 
      
      INTEGER,  INTENT(IN) :: iesr
      REAL(DP), INTENT(IN) :: taus(3,nat)
      REAL(DP) :: ESR
      REAL(DP) :: DESR(6)
      REAL(DP) :: FION(3,nat)
      LOGICAL,  INTENT(IN) :: TSTRESS
      REAL(DP), INTENT(in) :: hmat( 3, 3 )

! ... declare external function
      REAL(DP) :: erf, erfc
      EXTERNAL erf, erfc

      INTEGER :: ldim_block, gind_block
      EXTERNAL ldim_block, gind_block

      
! ... LOCALS 

      INTEGER :: na_loc, ia_s, ia_e, igis
      INTEGER :: k, i, j, l, m, is, ia, infm, ix, iy, iz, ishft
      INTEGER :: npt, isa, me
      INTEGER :: iakl, iajm
      LOGICAL :: split, tzero, tshift
      INTEGER, ALLOCATABLE   :: iatom(:,:)
      REAL(DP), ALLOCATABLE :: zv2(:,:)
      REAL(DP), ALLOCATABLE :: rc(:,:)  
      REAL(DP), ALLOCATABLE :: fionloc(:,:) 
      REAL(DP)  :: RXLM(3), SXLM(3)
      REAL(DP)  :: xlm, ylm, zlm, erre2, rlm, arg, esrtzero
      REAL(DP)  :: addesr, addpre, repand, fxx
      REAL(DP)  :: rckj_m1
      REAL(DP)  :: zvk, zvj, zv2_kj
      REAL(DP)  :: fact_pre
      REAL(DP)  :: iasp( nsp )

      INTEGER, DIMENSION(6), PARAMETER :: ALPHA = (/ 1,2,3,2,3,3 /)
      INTEGER, DIMENSION(6), PARAMETER :: BETA  = (/ 1,1,1,2,2,3 /)

! ... SUBROUTINE BODY 

      me = me_image + 1

      !  get the index of the first atom of each specie

      isa = 1
      DO is = 1, nsp
         iasp( is ) = isa
         isa = isa + na( is )
      END DO

      !  Here count the pairs of atoms

      npt = 0
      DO k = 1, nsp
        DO j = k, nsp
          DO l = 1, na(k)
            IF ( k == j ) THEN
              infm = l             ! If the specie is the same avoid  
            ELSE                   ! atoms double counting
              infm = 1
            END IF
            DO m = infm, na(j)
              npt = npt + 1
            END DO
          END DO
        END DO
      END DO

      ALLOCATE( iatom( 4, npt ) )
      ALLOCATE( rc( nsp, nsp ) )
      ALLOCATE( zv2( nsp, nsp ) )
      ALLOCATE( fionloc( 3, nat ) )
      rc      = 0.0_DP
      zv2     = 0.0_DP
      fionloc = 0.0_DP

      !  Here pre-compute some factors

      DO k = 1, nsp
        DO j = k, nsp
          zv2( k, j ) = zv( k ) * zv( j )
          rc ( k, j ) = SQRT( rcmax(k)**2 + rcmax(j)**2 )
        END DO
      END DO

      !  Here store the indexes of all pairs of atoms

      npt = 0
      DO k = 1, nsp
        DO j = k, nsp
          DO l = 1, na(k)
            IF (k.EQ.j) THEN
              infm = l
            ELSE
              infm = 1
            END IF
            DO m = infm, na(j)
              npt = npt + 1
              iatom(1,npt) = k
              iatom(2,npt) = j
              iatom(3,npt) = l
              iatom(4,npt) = m
            END DO
          END DO
        END DO
      END DO

      xlm     = 1.0_DP
      ylm     = 1.0_DP
      zlm     = 1.0_DP
      ESR     = 0.0_DP
      DESR    = 0.0_DP

      !  Distribute the atoms pairs to processors

      NA_LOC = ldim_block( npt, nproc_image, me_image)
      IA_S   = gind_block( 1, npt, nproc_image, me_image )
      IA_E   = IA_S + NA_LOC - 1

      DO ia = ia_s, ia_e

        k = iatom(1,ia)
        j = iatom(2,ia)
        l = iatom(3,ia)
        m = iatom(4,ia)

        zv2_kj   = zv2(k,j)
        rckj_m1  = 1.0_DP / rc(k,j)
        fact_pre = (2.0_DP * zv2_kj * sqrtpm1) * rckj_m1

        iakl = iasp(k) + l - 1
        iajm = iasp(j) + m - 1

        IF( (l.EQ.m) .AND. (k.EQ.j)) THEN      
           ! ...     same atoms
           xlm=0.0_DP; ylm=0.0_DP; zlm=0.0_DP; 
           tzero=.TRUE.
        ELSE
           ! ...     different atoms
           xlm = taus(1,iakl) - taus(1,iajm)
           ylm = taus(2,iakl) - taus(2,iajm)
           zlm = taus(3,iakl) - taus(3,iajm)
           CALL pbcs(xlm,ylm,zlm,xlm,ylm,zlm,1)
           TZERO=.FALSE.
        END IF

        DO IX=-IESR,IESR
          SXLM(1) = XLM + DBLE(IX)
          DO IY=-IESR,IESR
            SXLM(2) = YLM + DBLE(IY)
            DO IZ=-IESR,IESR
              TSHIFT= IX.EQ.0 .AND. IY.EQ.0 .AND. IZ.EQ.0
              IF( .NOT. ( TZERO .AND. TSHIFT ) ) THEN
                SXLM(3) = ZLM + DBLE(IZ)
                CALL S_TO_R( SXLM, RXLM, hmat )
                ERRE2 = RXLM(1)**2 + RXLM(2)**2 + RXLM(3)**2
                RLM   = SQRT(ERRE2)
                ARG   = RLM * rckj_m1
                IF (TZERO) THEN
                   ESRTZERO=0.5D0
                ELSE
                   ESRTZERO=1.D0
                END IF
                ADDESR = ZV2_KJ * erfc(ARG) / RLM
                ESR    = ESR + ESRTZERO*ADDESR
                ADDPRE = FACT_PRE * EXP(-ARG*ARG)
                REPAND = ESRTZERO*(ADDESR + ADDPRE)/ERRE2
                !
                DO i = 1, 3
                   fxx = repand * rxlm( i )
                   fionloc( i, iakl ) = fionloc( i, iakl ) + fxx
                   fionloc( i, iajm ) = fionloc( i, iajm ) - fxx
                END DO
                !
                IF( tstress ) THEN
                   DO i = 1, 6
                      fxx = repand * rxlm( alpha( i ) ) * rxlm( beta( i ) )
                      desr( i ) = desr( i ) - fxx
                   END DO
                END IF
                !
              END IF
            END DO    ! IZ
          END DO      ! IY
        END DO        ! IX
      END DO

!
!     each processor add its own contribution to the array FION
!
      isa = 0
      DO IS = 1, nsp
        DO IA = 1, na(is)
          isa = isa + 1
          FION(1,ISA) = FION(1,ISA)+FIONLOC(1,ISA)
          FION(2,ISA) = FION(2,ISA)+FIONLOC(2,ISA)
          FION(3,ISA) = FION(3,ISA)+FIONLOC(3,ISA)
        END DO
      END DO

      CALL mp_sum(esr, intra_image_comm)
     
      DEALLOCATE(iatom)
      DEALLOCATE(rc)
      DEALLOCATE(zv2)
      DEALLOCATE(fionloc)
      
      RETURN
!=----------------------------------------------------------------------------=!
   END SUBROUTINE vofesr
!=----------------------------------------------------------------------------=!



!=----------------------------------------------------------------------------=!
   SUBROUTINE self_vofhar_x( tscreen, self_ehte, vloc, rhoeg, omega, hmat )
!=----------------------------------------------------------------------------=!

      !  adds the hartree part of the self interaction

      USE kinds,              ONLY: DP
      USE constants,          ONLY: fpi
      USE control_flags,      ONLY: gamma_only
      USE cell_base,          ONLY: tpiba2, boxdimensions
      USE gvecp,              ONLY: ngm
      USE reciprocal_vectors, ONLY: gstart, g
      USE sic_module,         ONLY: sic_epsilon, sic_alpha
      USE mp_global,          ONLY: intra_image_comm
      USE mp,                 ONLY: mp_sum

      IMPLICIT NONE

      ! ... Arguments
      LOGICAL     :: tscreen
      COMPLEX(DP) :: vloc(:)
      COMPLEX(DP) :: rhoeg(:,:)
      REAL(DP)    :: self_ehte
      REAL(DP), INTENT(IN) :: omega
      REAL(DP), INTENT(IN) :: hmat( 3, 3 )

      ! ... Locals

      INTEGER      :: ig
      REAL(DP)    :: fpibg
      COMPLEX(DP) :: rhog
      COMPLEX(DP) :: ehte
      COMPLEX(DP) :: vscreen
      COMPLEX(DP), ALLOCATABLE :: screen_coul(:)

      ! ... Subroutine body ...

      IF( tscreen ) THEN
        ALLOCATE( screen_coul( ngm ) )
        CALL cluster_bc( screen_coul, g, omega, hmat )
      END IF

      !==  HARTREE ==

      ehte = 0.D0

      DO IG = gstart, ngm

        rhog  = rhoeg(ig,1) - rhoeg(ig,2)

        IF( tscreen ) THEN
          FPIBG     = fpi / ( g(ig) * tpiba2 ) + screen_coul(ig)
        ELSE
          FPIBG     = fpi / ( g(ig) * tpiba2 )
        END IF

        vloc(ig) = fpibg * rhog
        ehte     = ehte   +  fpibg *   rhog * CONJG(rhog)

      END DO
 
      ! ... G = 0 element
      !
      IF ( gstart == 2 ) THEN
        rhog    = rhoeg(1,1) - rhoeg(1,2)
        IF( tscreen ) THEN
          vscreen = screen_coul(1)
        ELSE
          vscreen = 0.0d0
        END IF
        vloc(1) = vscreen * rhog
        ehte    = ehte   +  vscreen *  rhog * CONJG(rhog)
      END IF

      ! ...

      IF( .NOT. gamma_only ) THEN
        ehte  = ehte  * 0.5d0
      END IF
      !
      self_ehte = DBLE(ehte) * omega * sic_epsilon
      vloc = vloc * sic_epsilon

      CALL mp_sum( self_ehte, intra_image_comm )

      IF( ALLOCATED( screen_coul ) ) DEALLOCATE( screen_coul )

      RETURN
!=----------------------------------------------------------------------------=!
   END SUBROUTINE self_vofhar_x
!=----------------------------------------------------------------------------=!



!=----------------------------------------------------------------------------=!
   SUBROUTINE localisation_x( wfc, atoms_m, ht)
!=----------------------------------------------------------------------------=!


      USE kinds,              ONLY: DP
      USE constants, ONLY: fpi
      USE control_flags, ONLY: gamma_only
      USE atoms_type_module, ONLY: atoms_type
      USE sic_module, ONLY: ind_localisation, nat_localisation, print_localisation
      USE sic_module, ONLY: sic_rloc, pos_localisation
      USE ions_base, ONLY: ind_srt
      USE fft_base, ONLY: dfftp, dffts
      USE cell_base, ONLY: tpiba2, boxdimensions, s_to_r
      USE reciprocal_vectors, ONLY: gstart, g
      USE gvecp, ONLY: ngm
      USE gvecw, ONLY: ngw
      use grid_dimensions, only: nr1, nr2, nr3, nr1l, nr2l, nr3l, nnrx
      USE cp_interfaces, ONLY: fwfft, invfft

      IMPLICIT NONE

! ... Arguments

      COMPLEX(DP), INTENT(IN) :: wfc(:)
      TYPE (atoms_type), INTENT(in) :: atoms_m
      TYPE (boxdimensions), INTENT(in) :: ht


! ... Locals

      REAL(DP)    :: ehte
      INTEGER      :: ig, at, ia, is, isa_input, isa_sorted, isa_loc
      REAL(DP)    :: fpibg, omega, aRe, aR2, R(3)
      INTEGER      :: Xmin, Ymin, Zmin, Xmax, Ymax, Zmax, i,j,k, ir
      REAL(DP)    :: work, work2
      COMPLEX(DP) :: rhog
      COMPLEX(DP), ALLOCATABLE :: density(:), psi(:)
      COMPLEX(DP), ALLOCATABLE :: k_density(:)
      COMPLEX(DP) :: vscreen
      COMPLEX(DP), ALLOCATABLE :: screen_coul(:)

! ... Subroutine body ...


      IF( .FALSE. ) THEN
        ALLOCATE( screen_coul( ngm ) )
        CALL cluster_bc( screen_coul, g, ht%deth, ht%hmat )
      END IF


      omega = ht%deth

      ALLOCATE( density( nnrx ) )
      ALLOCATE( psi( nnrx ) )
      ALLOCATE( k_density( ngm ) )

      CALL c2psi(  psi, dffts%nnr, wfc, wfc, ngw, 1 )
      CALL invfft( 'Wave', psi, dffts%nr1, dffts%nr2, dffts%nr3, dffts%nr1x, dffts%nr2x, dffts%nr3x )

      psi = DBLE( psi )

      isa_sorted = 0
      isa_loc    = 0

      DO is = 1, atoms_m%nsp

        DO ia = 1, atoms_m%na( is )

          isa_sorted = isa_sorted + 1         ! index of the atom as is in the sorted %tau atom_type component
          isa_input  = ind_srt( isa_sorted )  ! index of the atom as is in the input card ATOMIC_POSITIONS
          
          IF( ind_localisation( isa_input ) > 0 ) THEN

            isa_loc = isa_loc + 1  ! index of the localised atom ( 1 ... nat_localisation )
 
            IF( isa_loc > SIZE( pos_localisation, 2 ) ) &
              CALL errore( ' localisation ', ' too many localization ', isa_loc )

            ehte = 0.D0
            R( : ) = atoms_m%taus( :, isa_sorted )
            CALL s_to_r ( R, pos_localisation( 1:3 , isa_loc ), ht )

            !WRITE(6,*) 'ATOM ', ind_localisation( isa_input )
            !WRITE(6,*) 'POS  ', atoms_m%taus( :, isa_sorted )

            work  = nr1l
            work2 = sic_rloc * work
            work  = work * R(1) - work2
            Xmin  = FLOOR(work)
            work  = work + 2*work2
            Xmax  = FLOOR(work)
            IF ( Xmax > nr1l ) Xmax = nr1l
            IF ( Xmin < 1 ) Xmin = 1

            work  = nr2l
            work2 = sic_rloc * work
            work  = work * R(2) - work2
            Ymin  = FLOOR(work)
            work  = work + 2*work2
            Ymax  = FLOOR(work)
            IF ( Ymax > nr2l ) Ymax = nr2l
            IF ( Ymin < 1 ) Ymin = 1

            work  = nr3l
            work2 = sic_rloc * work
            work  = work * R(3) - work2
            Zmin  = FLOOR(work)
            work  = work + 2*work2
            Zmax  = FLOOR(work)
            IF ( Zmax > nr3l ) Zmax = nr3l
            IF ( Zmin < 1 ) Zmin = 1

            density = 0.D0

            DO k = Zmin, Zmax
               DO j = Ymin, Ymax
                  DO i = Xmin, Xmax
                     ir =  i + (j-1)*dfftp%nr1x + (k-1)*dfftp%nr1x*dfftp%nr2x
                     density( ir ) = psi( ir ) * psi( ir )
                  END DO
               END DO
            END DO

            CALL fwfft(   'Dense', density, dfftp%nr1, dfftp%nr2, dfftp%nr3, dfftp%nr1x, dfftp%nr2x, dfftp%nr3x )
            CALL psi2rho( 'Dense', density, dfftp%nnr, k_density, ngm )

! ... G /= 0 elements

            DO IG = gstart, ngm

              rhog  = k_density(ig)
              IF( .FALSE. ) THEN
                FPIBG     = fpi / ( g(ig) * tpiba2 ) + screen_coul(ig)
              ELSE
                FPIBG     = fpi / ( g(ig) * tpiba2 )
              END IF

              ehte       = ehte   +  fpibg *   DBLE(rhog * CONJG(rhog))

            END DO

! ... G = 0 element

            IF ( gstart == 2 ) THEN
              IF( .FALSE. ) THEN
                vscreen = screen_coul(1)
              ELSE
                vscreen = 0.0d0
              END IF
              rhog    = k_density(1)
              ehte    = ehte   +  vscreen *  DBLE(rhog * CONJG(rhog))
            END IF
! ...
            IF( .NOT. gamma_only ) THEN
              ehte  = ehte  * 0.5d0
            END IF
            ehte = ehte * omega
            pos_localisation( 4, isa_loc ) = ehte

          END IF  ! ind_localisation

        END DO  ! ia

      END DO  ! is

      ! CALL errore( 'DEBUG', ' qui ', 1 )

! ...
      IF( ALLOCATED(screen_coul) ) DEALLOCATE( screen_coul )
      DEALLOCATE( k_density, density, psi )

      RETURN
      END SUBROUTINE localisation_x

