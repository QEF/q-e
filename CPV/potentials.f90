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

!=----------------------------------------------------------------------------=!
  MODULE potentials
!=----------------------------------------------------------------------------=!

        USE kinds,         ONLY: DP
        USE control_flags, ONLY: timing

        IMPLICIT NONE 
        PRIVATE
        SAVE

        LOGICAL   :: tvhmean
        INTEGER   :: vhnr, vhiunit
        REAL(DP)  :: vhrmin, vhrmax
        CHARACTER :: vhasse

        INTEGER   :: iesr

        REAL(DP)  :: timtot, timfwft, timesr, timsumg, timforc, timinvft, timscat
        REAL(DP)  :: timxc, timhar, timstr
        INTEGER   :: timcnt = 0

        PUBLIC :: vofrhos, potential_init, potential_print_info, &
                  kspotential, print_vofrho_time, localisation, vofesr

        PUBLIC :: self_vofhar

        REAL(DP), EXTERNAL :: cclock

!=----------------------------------------------------------------------------=!
  CONTAINS
!=----------------------------------------------------------------------------=!

      SUBROUTINE potential_init &
          ( tvhmean_inp, vhnr_inp, vhiunit_inp, vhrmin_inp, vhrmax_inp, &
            vhasse_inp, iesr_inp )

          LOGICAL,   INTENT(IN) :: tvhmean_inp
          INTEGER,   INTENT(IN) :: vhnr_inp, vhiunit_inp
          REAL(DP),  INTENT(IN) :: vhrmin_inp, vhrmax_inp
          CHARACTER, INTENT(IN) :: vhasse_inp
          INTEGER,   INTENT(IN) :: iesr_inp

          tvhmean = tvhmean_inp
          vhnr    = vhnr_inp
          vhiunit = vhiunit_inp
          vhrmin  = vhrmin_inp
          vhrmax  = vhrmax_inp
          vhasse  = vhasse_inp
          iesr    = iesr_inp

          timtot = 0.0d0
          timfwft = 0.0d0
          timesr = 0.0d0
          timsumg = 0.0d0
          timforc = 0.0d0
          timinvft = 0.0d0
          timscat  = 0.0d0
          timxc  = 0.0d0 
          timhar = 0.0d0
          timstr = 0.0d0
          timcnt = 0

        RETURN
      END SUBROUTINE potential_init


      SUBROUTINE potential_print_info( iunit )

        INTEGER, INTENT(IN) :: iunit

        WRITE(iunit,50)
        WRITE(iunit,115) (2*iesr+1),(2*iesr+1),(2*iesr+1)

   50   FORMAT(//,3X,'Potentials Parameters',/,3X,'---------------------')
  115   FORMAT(   3X,'Ewald sum over ',I1,'*',I1,'*',I1,' cells')

        RETURN
      END SUBROUTINE potential_print_info
     

      SUBROUTINE vofmean( sfac, rhops, rhoeg )

        USE constants, ONLY: fpi
        USE cell_base, ONLY: tpiba2, tpiba
        USE mp,        ONLY: mp_sum
        USE mp_global, ONLY: nproc_image, me_image, intra_image_comm
        USE io_global, ONLY: ionode
        USE gvecp,     ONLY: ngm
        USE reciprocal_vectors, ONLY: gstart, gx, g

        REAL(DP),    INTENT(IN)   :: RHOPS(:,:)
        COMPLEX(DP), INTENT(IN)   :: RHOEG(:)
        COMPLEX(DP), INTENT(IN)   :: sfac(:,:)

        COMPLEX(DP) :: fpi_tpiba2, rp, vcg
        REAL(DP)    :: gxt, dr, r
        REAL(DP), ALLOCATABLE :: vrmean(:) 
        
        INTEGER :: ig, is, iasse, ipiano1, ipiano2
        INTEGER :: ir

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
          DO ir = 1, vhnr
            r = vhrmin + (ir-1)*dr
            WRITE(vhiunit,100) r, vrmean(ir)
          END DO
 100      FORMAT(2X,F14.8,F14.8)
        END IF
        
        DEALLOCATE(vrmean)

        RETURN
      END SUBROUTINE vofmean

!  -------------------------------------------------------------------------

      SUBROUTINE kspotential &
        ( nfi, tprint, tforce, tstress, rhoe, atoms, bec, becdr, eigr, &
          ei1, ei2, ei3, sfac, c0, cdesc, tcel, ht, fi, vpot, edft, timepre )

        USE charge_density,    ONLY: rhoofr
        USE nl,                ONLY: nlrh_m
        USE energies,          ONLY: dft_energy_type
        USE cell_module,       ONLY: boxdimensions
        USE atoms_type_module, ONLY: atoms_type
        USE wave_types,        ONLY: wave_descriptor

! ...   declare subroutine arguments
        LOGICAL   :: tcel
        INTEGER,              INTENT(IN)    :: nfi
        TYPE (atoms_type),    INTENT(INOUT) :: atoms
        COMPLEX(DP),         INTENT(INOUT) :: c0(:,:,:,:)
        TYPE (wave_descriptor),  INTENT(IN) :: cdesc
        REAL(DP) :: rhoe(:,:)
        COMPLEX(DP) :: ei1(:,:)
        COMPLEX(DP) :: ei2(:,:)
        COMPLEX(DP) :: ei3(:,:)
        COMPLEX(DP) :: eigr(:,:)
        TYPE (boxdimensions), INTENT(INOUT) ::  ht
        REAL(DP), INTENT(IN) :: fi(:,:,:)
        REAL(DP) :: bec(:,:)
        REAL(DP) :: becdr(:,:,:)
        TYPE (dft_energy_type) :: edft
        REAL(DP)    :: vpot(:,:)
        COMPLEX(DP), INTENT(IN) :: sfac(:,:)
        LOGICAL, INTENT(IN) :: tforce, tstress, tprint
        REAL(DP), INTENT(OUT) :: timepre

        edft%enl = nlrh_m( c0, cdesc, tforce, atoms, fi, bec, becdr, eigr )

        CALL rhoofr( nfi, c0, cdesc, fi, rhoe, ht )

        CALL vofrhos( tprint, tforce, tstress, rhoe, atoms, vpot, bec, &
                      c0, cdesc, fi, eigr, ei1, ei2, ei3, sfac, timepre,  &
                      ht, edft )

        RETURN
      END SUBROUTINE kspotential

!=----------------------------------------------------------------------------=!

   SUBROUTINE vofrhos &
      ( tprint, tforce, tstress, rhoe, atoms, vpot, bec, c0, cdesc, fi, &
        eigr, ei1, ei2, ei3, sfac, timepre, box, edft )

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

      USE control_flags,  ONLY: tscreen, tchi2, iprsta, force_pairing
      USE mp_global,      ONLY: nproc_image, me_image, intra_image_comm
      USE mp,             ONLY: mp_sum
      USE cell_module,    ONLY: boxdimensions
      USE cell_base,      ONLY: tpiba2
      USE ions_base,      ONLY: rcmax, zv
      USE fft_base,       ONLY: dfftp
      USE energies,       ONLY: total_energy, dft_energy_type
      USE stress,         ONLY: pstress
      USE funct,          ONLY: dft_is_gradient
      USE charge_density, ONLY: fillgrad
      USE chi2,           ONLY: rhochi, allocate_chi2, deallocate_chi2
      USE vanderwaals,    ONLY: tvdw, vdw
      USE wave_types,     ONLY: wave_descriptor
      USE io_global,      ONLY: ionode, stdout
      USE sic_module,     ONLY: self_interaction, sic_epsilon, sic_alpha !!TO ADD!!!
      USE gvecp,          ONLY: ngm
      USE local_pseudo,   ONLY: vps, rhops
      USE atom,           ONLY: nlcc
      USE core,           ONLY: nlcc_any, rhocg
      USE core,           ONLY: add_core_charge, core_charge_forces
      USE fft_module,     ONLY: fwfft, invfft
      !
      USE reciprocal_vectors,        ONLY: gx
      USE atoms_type_module,         ONLY: atoms_type
      USE exchange_correlation,      ONLY: exch_corr_energy
      use grid_dimensions,           only: nr1, nr2, nr3, nnrx

      IMPLICIT NONE

! ... declare subroutine arguments
      LOGICAL, INTENT(IN) :: tprint, tforce, tstress
      REAL(DP)            :: vpot(:,:)
      REAL(DP),    INTENT(IN) :: fi(:,:,:)
      REAL(DP)    :: bec(:,:)
      COMPLEX(DP) :: ei1(:,:)
      COMPLEX(DP) :: ei2(:,:)
      COMPLEX(DP) :: ei3(:,:)
      COMPLEX(DP) :: eigr(:,:)
      COMPLEX(DP),   INTENT(IN) :: c0(:,:,:,:)
      TYPE (atoms_type), INTENT(INOUT) :: atoms
      TYPE (wave_descriptor), INTENT(IN) :: cdesc
      TYPE (boxdimensions),    INTENT(INOUT) :: box
      TYPE (dft_energy_type) :: edft
      REAL(DP) :: rhoe(:,:)
      COMPLEX(DP), INTENT(IN) :: sfac(:,:)

      TYPE (dft_energy_type) :: edft_self

      REAL(DP) :: timepre

! ... declare functions
      REAL(DP)  DDOT

! ... declare other variables

      COMPLEX(DP), ALLOCATABLE :: vloc(:), self_vloc(:)
      COMPLEX(DP), ALLOCATABLE :: rho12(:), rhoeg(:,:)
      COMPLEX(DP), ALLOCATABLE :: rhoetg(:,:)
      COMPLEX(DP), ALLOCATABLE :: psi(:)
      REAL(DP), ALLOCATABLE :: rhoetr(:,:)
      REAL(DP), ALLOCATABLE :: fion_vdw(:,:)
      REAL(DP), ALLOCATABLE :: grho(:,:,:)
      REAL(DP), ALLOCATABLE :: v2xc(:,:,:)
      REAL(DP), ALLOCATABLE :: fion(:,:)

      REAL(DP), ALLOCATABLE :: self_rho(:,:)
      REAL(DP), ALLOCATABLE :: self_vpot(:,:)
      REAL(DP), ALLOCATABLE :: self_grho(:,:,:)
      REAL(DP), ALLOCATABLE :: self_v2xc(:,:,:)

      REAL(DP) ::  pail(3,3)

      COMPLEX(DP) :: ehtep
      REAL(DP)    :: self_sxcp, self_vxc

      REAL(DP)  :: summing1, summing2

      COMPLEX(DP) :: ehp, eps

      REAL(DP)  :: dum, sxcp, vxc, ehr, strvxc
      REAL(DP)  :: omega, desr(6), pesum(16)
      REAL(DP)  :: s0, s1, s2, s3, s4, s5, s6, s7, s8

      LOGICAL :: ttscreen, ttsic, tgc

      INTEGER ig1, ig2, ig3, is, ia, ig, isc, iflag, iss
      INTEGER ik, i, j, k, isa, idum, nspin, iswfc
      INTEGER :: ierr

      DATA iflag / 0 /
      SAVE iflag, desr

      REAL(DP), EXTERNAL :: enkin

!  end of declarations
!  ----------------------------------------------

      IF(timing) s0 = cclock()

      nspin = SIZE( rhoe, 2 )
      
      edft%evdw = 0.0d0
      !
      ! ttscreen = .TRUE.
      ttscreen = .FALSE.
      !
      ttsic    = ( ABS(self_interaction) /= 0  ) 

      omega    = box%deth
      !
      tgc      = dft_is_gradient()

      IF(tchi2) THEN
        CALL allocate_chi2(ngm)
      END IF

      ALLOCATE( rhoetr( nnrx, nspin ) )
      rhoetr = 0.0d0

      ALLOCATE( fion( 3, atoms%nat ) )

      fion = atoms%for( 1:3, 1:atoms%nat )
      !
      pail = box%pail

      IF(tgc) THEN
        ALLOCATE( grho( nnrx, 3, nspin ) )
        ALLOCATE( v2xc( nnrx, nspin, nspin) )
      ELSE
        ALLOCATE( grho( 1, 1, 1 ) )
        ALLOCATE( v2xc( 1, 1, 1 ) )
      END IF

      grho = 0.0d0
      v2xc = 0.0d0

      ALLOCATE( rhoeg (ngm, nspin) )
      ALLOCATE( rhoetg(ngm, nspin) )

      IF( ttsic ) THEN
      
        IF ( tgc ) THEN

          ALLOCATE(self_grho( nnrx, 3, nspin ), STAT = ierr)
          IF( ierr /= 0 ) CALL errore(' vofrhos ', ' allocating self_grho ', ierr)
    
          ALLOCATE(self_v2xc( nnrx, nspin, nspin ), STAT = ierr)
          IF( ierr /= 0 ) CALL errore(' vofrhos ', ' allocating self_v2xc ', ierr)
  
          self_grho  = 0.D0
          self_v2xc  = 0.D0

        END IF !on tgc

        ALLOCATE (self_vpot( nnrx, 2 ), STAT = ierr)
        IF( ierr /= 0 ) CALL errore(' vofrhos ', ' allocating self_vpot ', ierr)
   
        self_vpot  = 0.D0
     
        ALLOCATE (self_rho( nnrx, 2), STAT = ierr)
        IF( ierr /= 0 ) CALL errore(' vofrhos ', ' allocating self_rho ', ierr)

        self_rho  = 0.D0

      END IF !on self_interaction

      IF(timing) s1 = cclock()

      ! ... compute kinetic energy

      edft%ekin  = 0.0_DP
      edft%emkin = 0.0_DP

      DO iss = 1, nspin
        iswfc = iss
        IF( force_pairing ) iswfc = 1
        edft%ekin  = edft%ekin + enkin( c0(1,1,1,iswfc), SIZE(c0,1), fi(1,1,iss), cdesc%nbl(iss) )
      END DO

      IF(tprint) THEN
        IF( ionode .AND.  ttscreen ) &
           WRITE( stdout,fmt="(3X,'Using screened Coulomb potential for cluster calculation')")
      END IF

      IF( tstress .OR. tforce .OR. iflag == 0 )  THEN
         CALL vofesr( iesr, edft%esr, desr, fion, atoms%taus, tstress, box%hmat )
         IF( iflag == 0 ) &
            WRITE( stdout, fmt="(/,3X,'ESR (real part of Ewald sum) = ',D16.8,/)" ) edft%esr
         iflag = 1
         ! WRITE(6,*) 'DEBUG esr = ', SUM(fion)
      END IF


      IF(timing) s2 = cclock()

      DO iss = 1, nspin

        ! ... FFT: rho(r) --> rho(g)  

        ALLOCATE( psi( SIZE( rhoe, 1 ) ) )

        psi = rhoe(:,iss)

        CALL fwfft(   'Dense', psi, dfftp%nr1, dfftp%nr2, dfftp%nr3, dfftp%nr1x, dfftp%nr2x, dfftp%nr3x )
        CALL psi2rho( 'Dense', psi, dfftp%nnr, rhoeg(:,iss), ngm )
 
        DEALLOCATE( psi )

        ! ... add core contribution to the charge

        CALL ZCOPY( SIZE(rhoeg,1) , rhoeg(1,iss), 1, rhoetg(1,iss), 1 )
        CALL DCOPY( nnrx, rhoe(1,iss), 1, rhoetr(1,iss), 1 )

        IF( nlcc_any ) THEN

          ! ...     add core correction
          ! ...     rhoetg = rhoeg + cc
          ! ...     rhoetr = rhoe  + cc

          CALL add_core_charge( rhoetg(:,iss), rhoetr(:,iss), sfac, rhocg, atoms%nsp )

        ELSE

          ! ...     no core correction
          ! ...     rhoetg = rhoeg
          ! ...     rhoetr = rhoe

          ! ...     chi2

          IF( tchi2 ) THEN
            IF(nspin.GT.1) CALL errore(' vofrho ',' spin + tchi ',nspin)
            rhochi = rhoeg(:,1)
          END IF

        END IF

      END DO

      IF( tgc ) THEN
         CALL fillgrad( nspin, rhoetg, grho )
      END IF

      IF(timing) s3 = cclock()


      CALL exch_corr_energy(rhoetr, grho, vpot, sxcp, vxc, v2xc)

      edft%sxc       = sxcp
      edft%self_sxc  = 0.d0
      self_vxc       = 0.d0
      !
      IF ( ttsic ) THEN                

         self_rho(:,1) = rhoetr(:,2)
         self_rho(:,2) = rhoetr(:,2)

         IF (tgc) THEN
            self_grho(:,:,1) = grho(:,:,2)
            self_grho(:,:,2) = grho(:,:,2)
         ENDIF

         CALL exch_corr_energy( self_rho, self_grho, self_vpot, &
                 self_sxcp, self_vxc, self_v2xc )

         vpot (:,1) = ( 1.0d0 - sic_alpha ) * vpot(:,1)
         ! 
         vpot (:,2) = ( 1.0d0 - sic_alpha ) * vpot(:,2) + &
                      sic_alpha * ( self_vpot(:,2) + self_vpot(:,1) )

         IF (tgc) THEN
           !
           v2xc(:,1,1) = ( 1.0d0 - sic_alpha ) * v2xc(:,1,1)
           !
           v2xc(:,2,2) = ( 1.0d0 - sic_alpha ) * v2xc(:,2,2) + & 
                         sic_alpha * ( self_v2xc(:,2,2) + self_v2xc(:,1,1) )
           !
         END IF

         edft%self_sxc = sic_alpha * ( sxcp - self_sxcp )
         !
         edft%sxc      = sxcp - edft%self_sxc
         !
         self_vxc      = sic_alpha * ( vxc - self_vxc )
         !
         vxc           = vxc - self_vxc
         !
      END IF  

      IF ( tstress ) THEN
        strvxc = ( edft%sxc - vxc ) * omega / DBLE( nr1 * nr2 * nr3 )
      END IF

      IF( nlcc_any ) THEN
        !
        ! ...   xc potential (vpot) from real to G space, to compute nlcc forces
        ! ...   rhoetg = fwfft(vpot)
        !
        ALLOCATE( psi( SIZE( vpot, 1 ) ) )

        DO iss = 1, nspin
           !
           psi = vpot(:,iss)
           !
           CALL fwfft(   'Dense', psi, dfftp%nr1, dfftp%nr2, dfftp%nr3, dfftp%nr1x, dfftp%nr2x, dfftp%nr3x )
           CALL psi2rho( 'Dense', psi, dfftp%nnr, rhoetg(:,iss), ngm )
           ! 
        END DO

        DEALLOCATE( psi )
        !
        ! ...   now rhoetg contains the xc potential
        !
        IF (tforce) THEN
          CALL core_charge_forces( fion, rhoetg, rhocg, nlcc, atoms, box, ei1, ei2, ei3 )
        END IF
        !
      END IF


      ! ... Van Der Waals energy and forces
      !
      IF ( tvdw ) THEN
        CALL VdW( edft%evdw, atoms, fion, box )
      END IF

      IF(timing) s4 = cclock()


      ! ... Calculate hartree potential and energy (eh), and
      ! ... local part of the pseudopotential and its energy contribution (eps)
      ! ... Self-interaction correction --- Hartree part

      ALLOCATE( vloc( ngm ) )
      !
      CALL vofloc(ttscreen, tforce, edft%ehte, edft%ehti, ehp, & 
           eps, vloc, rhoeg, fion, atoms, rhops, vps, eigr, &
           ei1, ei2, ei3, sfac, box )
      !
      edft%self_ehte = 0.d0

      IF ( ttsic ) THEN

        ALLOCATE( self_vloc( ngm ) )
        ALLOCATE( psi( nnrx ) )

        CALL self_vofhar( ttscreen, edft%self_ehte, self_vloc, rhoeg, omega, box%hmat )
        !
        CALL rho2psi( 'Dense', psi, dfftp%nnr, self_vloc, ngm )
        CALL invfft(  'Dense', psi, dfftp%nr1, dfftp%nr2, dfftp%nr3, dfftp%nr1x, dfftp%nr2x, dfftp%nr3x )
        !
        vpot(:,1) =  vpot(:,1) - DBLE( psi(:) )
        vpot(:,2) =  vpot(:,2) + DBLE( psi(:) )

        DEALLOCATE( self_vloc, psi )

      END IF

      edft%eh = DBLE( ehp ) - edft%self_ehte
     
      IF( ALLOCATED( self_grho  ) ) DEALLOCATE( self_grho  )
      IF( ALLOCATED( self_v2xc  ) ) DEALLOCATE( self_v2xc  )
      IF( ALLOCATED( self_vpot  ) ) DEALLOCATE( self_vpot )
      IF( ALLOCATED( self_rho   ) ) DEALLOCATE( self_rho  )

      ! ... vloc(g): hartree and local part of the pseudo potentials (in
      ! ...          reciprocal space
 
      IF(timing) s5 = cclock()


      ALLOCATE( psi( SIZE( vpot, 1 ) ) )

      DO iss = 1, nspin

        ! ...   add hartree end local pseudo potentials ( invfft(vloc) )
        ! ...   to xc potential (vpot).
        ! ...   vpot = vpot + invfft(vloc)

        CALL rho2psi( 'Dense', psi, dfftp%nnr, vloc, ngm )
        CALL invfft(  'Dense', psi, dfftp%nr1, dfftp%nr2, dfftp%nr3, dfftp%nr1x, dfftp%nr2x, dfftp%nr3x )

        vpot(:,iss) = vpot(:,iss) + DBLE( psi )

      END DO

      DEALLOCATE( vloc, psi )

! ... now potentials are in real space
! ... vpot(r) = hartree + xc + local

      IF(timing) s6 = cclock()

! ... sum up forces
      IF (tforce) THEN
        CALL mp_sum(fion, intra_image_comm)
      END IF
      ! WRITE(6,*) 'DEBUG end = ', SUM(fion)

! ... sum up energies
      CALL mp_sum(eps, intra_image_comm)
      CALL mp_sum(edft%sxc, intra_image_comm)
      CALL mp_sum(edft%self_sxc, intra_image_comm)
      CALL mp_sum(vxc, intra_image_comm)
      CALL mp_sum(edft%eh, intra_image_comm)
      CALL mp_sum(edft%ehte, intra_image_comm)
      CALL mp_sum(edft%ehti, intra_image_comm)
      CALL mp_sum(edft%self_ehte, intra_image_comm)
      ! CALL mp_sum(edft%ekin, group)  ! already summed up
      CALL mp_sum(edft%emkin, intra_image_comm)

      CALL total_energy(edft,omega,vxc,eps,self_vxc,nr1*nr2*nr3)

!fran: the output is introduced only in the print_energies.f90
!fran: in this way you print only each print_step

      ! ... compute stress tensor
      !
      IF( tstress ) THEN
        s8 = cclock()
        CALL pstress( strvxc, rhoeg, rhoetg, pail, desr, bec, c0, cdesc, fi,  &
                      eigr, sfac, grho, v2xc, box, edft)
        timepre = cclock() - s8
      END IF


! ... Copy new atomic forces on for type member
!
      atoms%for( 1:3, 1:atoms%nat ) = fion
      box%pail = pail

! ... print out energies, stress and forces on file 6
! ... just a proof to calculate vhmean for sic

      IF(tprint.AND.tvhmean) THEN
        IF(nspin.GT.2) CALL errore(' vofrho ',' spin + vmean ',nspin)
        IF(nspin==1) CALL vofmean( sfac, rhops, rhoeg(:,1) )
        IF(nspin==2) THEN
             ALLOCATE(rho12(ngm))
            rho12 (:) = rhoeg(:,1)+rhoeg(:,2)
          CALL vofmean( sfac, rhops, rho12)
           DEALLOCATE(rho12)
        END IF
      END IF

      DEALLOCATE( rhoeg, rhoetg, rhoetr, grho, v2xc, fion )

      IF(tchi2) THEN
        CALL deallocate_chi2
      END IF

      IF(iprsta>2) THEN
        CALL memstat(me_image)
      END IF

      IF(timing) THEN
        s7 = cclock()
        timtot   = (s7 - s0) + timtot
        timstr   = (s7 - s6) + timstr
        timinvft = (s6 - s5) + timinvft
        timhar   = (s5 - s4) + timhar
        timxc    = (s4 - s3) + timxc
        timfwft  = (s3 - s2) + timfwft
        timesr   = (s2 - s1) + timesr
        timcnt   = timcnt + 1
      END IF

! ... Flush stdout
!
      CALL flush_unit( stdout )

      RETURN
      END SUBROUTINE vofrhos

!=----------------------------------------------------------------------------=!

      SUBROUTINE print_vofrho_time( iunit )
        USE io_global, ONLY: ionode
        INTEGER, INTENT(IN) :: iunit
        IF( timing .AND. timcnt > 0 ) THEN
          timesr  = timesr/timcnt
          timfwft  = timfwft/timcnt
          timxc  = timxc/timcnt
          timhar = timhar/timcnt
          timinvft  = timinvft/timcnt
          timstr  = timstr/timcnt
          timtot   = timtot/timcnt
          IF(ionode) THEN
            WRITE( iunit, 999 ) timesr, timfwft, timxc, timhar, timinvft, timstr, timtot
          END IF
999       FORMAT(1X,7(1X,F9.3))
        END IF
        timtot = 0.0d0
        timfwft = 0.0d0
        timesr = 0.0d0
        timsumg = 0.0d0
        timforc = 0.0d0
        timinvft = 0.0d0
        timscat  = 0.0d0
        timxc  = 0.0d0 
        timhar = 0.0d0
        timstr = 0.0d0
        timcnt = 0
        RETURN
      END SUBROUTINE print_vofrho_time
     

!=----------------------------------------------------------------------------=!

  SUBROUTINE cluster_bc( screen_coul, hg, omega, hmat )

      USE green_functions, ONLY: greenf
      USE mp_global,       ONLY: me_image
      USE fft_base,        ONLY: dfftp
      USE fft_module,      ONLY: fwfft
      USE gvecp,           ONLY: ngm
      USE cell_module,     ONLY: boxdimensions, s_to_r, alat
      USE constants,       ONLY: gsmall, pi
      USE cell_base,       ONLY: tpiba2
      use grid_dimensions, only: nr1, nr2, nr3, nr1l, nr2l, nr3l, nnrx
      
      REAL(DP), INTENT(IN) :: hg(:)
      REAL(DP), INTENT(IN) :: omega, hmat( 3, 3 )
      COMPLEX(DP) :: screen_coul(:)

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

!  ----------------------------------------------
!  BEGIN manual

      SUBROUTINE vofloc(tscreen, tforce, ehte, ehti, eh, eps, vloc, rhoeg, &
                 fion, atoms, rhops, vps, eigr, ei1, ei2, ei3, sfac, ht )

!  this routine computes:
!  omega = ht%deth
!  rho_e(ig)    =  (sum over iss) rhoeg(ig,iss) 
!  rho_I(ig)    =  (sum over is) sfac(is,ig) * rhops(ig,is) 
!  vloc_h(ig)   =  fpi / ( g(ig) * tpiba2 ) * { rho_e(ig) + rho_I(ig) }
!  vloc_ps(ig)  =  (sum over is) sfac(is,ig) * vps(ig,is)
!
!  Eps = Fact * omega * (sum over ig) cmplx ( rho_e(ig) ) * vloc_ps(ig)
!  if Gamma symmetry Fact = 2 else Fact = 1
!
!  Eh  = Fact * omega * (sum over ig) * fpi / ( g(ig) * tpiba2 ) *
!        { rho_e(ig) + rho_I(ig) } * conjugate { rho_e(ig) + rho_I(ig) }
!  if Gamma symmetry Fact = 1 else Fact = 1/2
!
!  Hatree potential and local pseudopotential
!  vloc(ig)     =  vloc_h(ig) + vloc_ps(ig) 
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
!  ----------------------------------------------
!  END manual

      USE constants, ONLY: fpi
      USE control_flags, ONLY: gamma_only
      USE cell_base, ONLY: tpiba2, tpiba
      USE cell_module, ONLY: boxdimensions
      USE atoms_type_module, ONLY: atoms_type
      USE io_global, ONLY: stdout
      USE grid_dimensions, ONLY: nr1, nr2, nr3
      USE reciprocal_vectors, ONLY: mill_l
      USE ions_base, ONLY: nat
      USE gvecp, ONLY: ngm
      USE reciprocal_vectors, ONLY: gstart, gx, g

      IMPLICIT NONE

! ... Arguments

      TYPE (atoms_type) :: atoms
      TYPE (boxdimensions), INTENT(in) :: ht
      LOGICAL      :: tforce
      LOGICAL      :: tscreen
      REAL(DP)    :: fion(:,:)
      REAL(DP)    :: rhops(:,:), vps(:,:)
      COMPLEX(DP) :: vloc(:)
      COMPLEX(DP) :: rhoeg(:,:)
      COMPLEX(DP), INTENT(IN) :: sfac(:,:)
      REAL(DP)    :: ehte, ehti
      COMPLEX(DP) :: eh, eps
      COMPLEX(DP) :: ei1(-nr1:nr1,nat)
      COMPLEX(DP) :: ei2(-nr2:nr2,nat)
      COMPLEX(DP) :: ei3(-nr3:nr3,nat)
      COMPLEX(DP) :: eigr(:,:)

! ... Locals

      INTEGER      :: is, ia, isa, ig, ig1, ig2, ig3, nspin, iss
      REAL(DP)    :: fpibg, cost, omega
      COMPLEX(DP) :: cxc, rhet, rhog, vp, rp, gxc, gyc, gzc
      COMPLEX(DP) :: teigr, cnvg, cvn, tx, ty, tz, vscreen
      COMPLEX(DP), ALLOCATABLE :: ftmp(:,:)
      COMPLEX(DP), ALLOCATABLE :: screen_coul(:)

! ... Subroutine body ...

      nspin = SIZE(rhoeg,2)
      IF(TFORCE) THEN
        ALLOCATE( ftmp(3, SIZE(fion, 2) ) )
        ftmp = CMPLX(0.0_DP,0.0_DP)
      END IF

      IF( tscreen ) THEN
        ALLOCATE( screen_coul( ngm ) )
        CALL cluster_bc( screen_coul, g, ht%deth, ht%hmat )
      END IF

!=======================================================================
!==  HARTREE AND LOCAL PART OF THE PSEUDO ENERGIES                    ==
!=======================================================================

      omega = ht%deth
      EH    = (0.D0,0.D0)
      EPS   = (0.D0,0.D0)
      ehte  = 0.0d0
      ehti  = 0.0d0

      DO IG = gstart, ngm 

        RP   = (0.D0,0.D0)
        VP   = (0.D0,0.D0)
        DO IS = 1, atoms%nsp
          VP = VP + sfac( ig, is ) * vps( ig, is )
          RP = RP + sfac( ig, is ) * rhops( ig, is )
        END DO

        ! WRITE( stdout,*) 'vp ',ig, vp  ! DEBUG
        ! WRITE( stdout,*) 'rp ',ig, rp  ! DEBUG
        ! WRITE( stdout,*) 'rhoeg ',ig, SUM( RHOEG( ig, : ) )  ! DEBUG
        ! WRITE( stdout,*) 'mill ',ig, mill(1,ig),mill(2,ig),mill(3,ig)

        RHET  = SUM( RHOEG( ig, : ) )
        RHOG  = RHET + RP

        ! WRITE( stdout,*) 'rhet ',ig, rhet  ! DEBUG
        ! WRITE( stdout,*) 'rhog ',ig, rhog  ! DEBUG

        IF( tscreen ) THEN
          FPIBG     = fpi / ( g(ig) * tpiba2 ) + screen_coul(ig)
        ELSE
          FPIBG     = fpi / ( g(ig) * tpiba2 )
        END IF

        ! WRITE( stdout,*) 'fpibg ',ig, fpibg  ! DEBUG
        
        vloc(ig) = vp   +  fpibg *        rhog 
        eh       = eh   +  fpibg *        rhog * CONJG(rhog)
        eps      = eps  +     vp * CONJG(rhet)
        ehte     = ehte +  fpibg *   DBLE(rhet * CONJG(rhet))
        ehti     = ehti +  fpibg *   DBLE(  rp * CONJG(rp))

        IF(TFORCE) THEN
          ig1  = mill_l(1,IG)
          ig2  = mill_l(2,IG)
          ig3  = mill_l(3,IG)
          GXC  = CMPLX(0.D0,gx(1,IG))
          GYC  = CMPLX(0.D0,gx(2,IG))
          GZC  = CMPLX(0.D0,gx(3,IG))
          isa = 1
          DO IS = 1, atoms%nsp
            CNVG  = RHOPS(IG,is) * FPIBG * CONJG(rhog)
            CVN   = VPS(ig, is)  * CONJG(rhet)
            TX = (CNVG+CVN) * GXC
            TY = (CNVG+CVN) * GYC
            TZ = (CNVG+CVN) * GZC
            DO IA = 1, atoms%na(is)
              TEIGR = ei1(IG1,ISA) * ei2(IG2,ISA) * ei3(IG3,ISA)
              ftmp(1,ISA) = ftmp(1,ISA) + TEIGR*TX
              ftmp(2,ISA) = ftmp(2,ISA) + TEIGR*TY
              ftmp(3,ISA) = ftmp(3,ISA) + TEIGR*TZ
              isa = isa + 1
            END DO
          END DO
        END IF

      END DO
! ... 
      IF(TFORCE) THEN
! ...   each processor add its own contribution to the array FION
        IF( gamma_only ) THEN
          cost = 2.D0 * omega * tpiba
        ELSE
          cost = omega * tpiba
        END IF
        FION = FION + DBLE(ftmp) * cost
      END IF

! ... G = 0 element
      IF ( gstart == 2 ) THEN
        vp = (0.D0,0.D0)
        rp = (0.D0,0.D0)
        IF( tscreen ) THEN
          vscreen = screen_coul(1)
        ELSE
          vscreen = 0.0d0
        END IF
        DO IS = 1, atoms%nsp
          vp = vp + sfac( 1, is) * vps(1, is)
          rp = rp + sfac( 1, is) * rhops(1, is)
        END DO
        rhet    = SUM( rhoeg(1, :) )
        rhog    = rhet + rp
        vloc(1) = VP   +  vscreen *   rhog
        eh      = eh   +  vscreen *        rhog * CONJG(rhog)
        ehte    = ehte +  vscreen *   DBLE(rhet * CONJG(rhet))
        ehti    = ehti +  vscreen *   DBLE(  rp * CONJG(rp))
        DO iss = 1, nspin
          IF( gamma_only ) THEN
            eps = eps + vp * CONJG(RHOEG(1,iss)) * 0.5d0
          ELSE
            eps = eps + vp * CONJG(RHOEG(1,iss))
          END IF
        END DO
      END IF
! ...
      IF( .NOT. gamma_only ) THEN
        EPS = EPS * 0.5d0
        EH  = EH  * 0.5d0
      END IF
      eh   =        eh   * omega
      eps  = 2.D0 * eps  * omega
      ehte =        ehte * omega
      ehti =        ehti * omega
! ...
      IF(ALLOCATED(ftmp)) DEALLOCATE(ftmp)
      IF(ALLOCATED(screen_coul)) DEALLOCATE(screen_coul)
       
      RETURN
      END SUBROUTINE vofloc

!
!=----------------------------------------------------------------------------=!
   SUBROUTINE vofesr( iesr, esr, desr, fion, taus, tstress, hmat )
!=----------------------------------------------------------------------------=!

      USE constants,   ONLY : sqrtpm1
      USE cell_module, ONLY : s_to_r, pbcs
      USE mp_global,   ONLY : nproc_image, me_image, intra_image_comm
      USE mp,          ONLY : mp_sum
      USE ions_base,   ONLY : rcmax, zv, nsp, na, nat
 
      IMPLICIT NONE

! ... ARGUMENTS 
      
      INTEGER,  INTENT(IN) :: iesr
      REAL(DP), INTENT(IN) :: taus(:,:)
      REAL(DP) :: ESR
      REAL(DP) :: DESR(:)
      REAL(DP) :: FION(:,:)
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
   SUBROUTINE self_vofhar( tscreen, self_ehte, vloc, rhoeg, omega, hmat )
!=----------------------------------------------------------------------------=!

      !  adds the hartree part of the self interaction

      USE constants,          ONLY: fpi
      USE control_flags,      ONLY: gamma_only
      USE cell_module,        ONLY: boxdimensions
      USE cell_base,          ONLY: tpiba2
      USE gvecp,              ONLY: ngm
      USE reciprocal_vectors, ONLY: gstart, g
      USE sic_module,         ONLY: sic_epsilon, sic_alpha

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

      IF( ALLOCATED( screen_coul ) ) DEALLOCATE( screen_coul )

      RETURN
!=----------------------------------------------------------------------------=!
   END SUBROUTINE self_vofhar
!=----------------------------------------------------------------------------=!



      SUBROUTINE localisation( wfc, atoms_m, ht)

!  adds the hartree part of the self interaction
!
!  ----------------------------------------------
!  END manual

      USE constants, ONLY: fpi
      USE control_flags, ONLY: gamma_only
      USE cell_module, ONLY: boxdimensions, s_to_r
      USE atoms_type_module, ONLY: atoms_type
      USE sic_module, ONLY: ind_localisation, nat_localisation, print_localisation
      USE sic_module, ONLY: sic_rloc, pos_localisation
      USE ions_base, ONLY: ind_srt
      USE fft_base, ONLY: dfftp, dffts
      USE cell_base, ONLY: tpiba2
      USE reciprocal_vectors, ONLY: gstart, g
      USE gvecp, ONLY: ngm
      USE gvecw, ONLY: ngw
      use grid_dimensions, only: nr1, nr2, nr3, nr1l, nr2l, nr3l, nnrx
      USE fft_module, ONLY: fwfft, invfft

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
      END SUBROUTINE localisation


 
!=----------------------------------------------------------------------------=!
   END MODULE potentials
!=----------------------------------------------------------------------------=!
