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

        USE kinds
        USE control_flags, ONLY: timing

        IMPLICIT NONE 
        PRIVATE
        SAVE

        LOGICAL   :: tvhmean
        INTEGER   :: vhnr, vhiunit
        REAL(dbl) :: vhrmin, vhrmax
        CHARACTER :: vhasse

        INTEGER   :: iesr

        REAL(dbl)  :: timtot, timfwft, timesr, timsumg, timforc, timinvft, timscat
        REAL(dbl)  :: timxc, timhar, timstr
        INTEGER    :: timcnt = 0

        PUBLIC :: vofrhos, potential_init, potential_print_info, &
                  kspotential, print_vofrho_time, localisation

        REAL(dbl), EXTERNAL :: cclock

!=----------------------------------------------------------------------------=!
  CONTAINS
!=----------------------------------------------------------------------------=!

      SUBROUTINE potential_init(tvhmean_inp,vhnr_inp, vhiunit_inp, &
          vhrmin_inp, vhrmax_inp, vhasse_inp, iesr_inp)

          LOGICAL, INTENT(IN) :: tvhmean_inp
          INTEGER, INTENT(IN) :: vhnr_inp, vhiunit_inp
          REAL(dbl), INTENT(IN)  :: vhrmin_inp, vhrmax_inp
          CHARACTER, INTENT(IN) :: vhasse_inp
          INTEGER, INTENT(IN) :: iesr_inp

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

!! ...................................................................... !!
!! ...................................................................... !!

      SUBROUTINE potential_print_info(iunit)

        INTEGER, INTENT(IN) :: iunit

        WRITE(iunit,50)

        WRITE(iunit,115) (2*iesr+1),(2*iesr+1),(2*iesr+1)

   50   FORMAT(//,3X,'Potentials Parameters',/,3X,'---------------------')
  115   FORMAT(   3X,'Ewald sum over ',I1,'*',I1,'*',I1,' cells')

        RETURN

      END SUBROUTINE potential_print_info
     
!! ...................................................................... !!
!! ...................................................................... !!

      SUBROUTINE vofmean( sfac, rhops, rhoeg )

        USE constants, ONLY: fpi
        USE cell_base, ONLY: tpiba2, tpiba
        USE mp, ONLY: mp_sum
        USE mp_global, ONLY: nproc, mpime, group, root
        USE io_global, ONLY: ionode
        USE reciprocal_vectors, ONLY: gstart, gx, g
        USE gvecp, ONLY: ngm

        REAL(dbl),    INTENT(IN)   :: RHOPS(:,:)
        COMPLEX(dbl), INTENT(IN)   :: RHOEG(:)
        COMPLEX(dbl), INTENT(IN)   :: sfac(:,:)

        COMPLEX(dbl) :: fpi_tpiba2, rp, vcg
        REAL(dbl)    :: gxt, dr, r
        REAL(dbl), ALLOCATABLE :: vrmean(:) 
        
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
              vrmean(ir) = vrmean(ir) + REAL(vcg)  * COS(gxt*r) 
              vrmean(ir) = vrmean(ir) - AIMAG(vcg) * SIN(gxt*r)
            END IF
          END DO
          vrmean(ir) = 2.0d0 * vrmean(ir)
        END DO
        CALL mp_sum(vrmean,group)

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
      SUBROUTINE kspotential( tprint, tforce, tstress, rhoe, desc, &
        atoms, kp, ps, eigr, ei1, ei2, ei3, sfac, c0, cdesc, tcel, ht, fi, fnl, vpot, edft, timepre )

        USE charge_density, ONLY: rhoofr
        USE nl, ONLY: nlrh_m
        USE energies, ONLY: dft_energy_type
        USE cell_module, ONLY: boxdimensions
        USE brillouin, ONLY: kpoints
        USE cp_types, ONLY: pseudo
        USE atoms_type_module, ONLY: atoms_type
        USE wave_types, ONLY: wave_descriptor
        USE pseudo_projector, ONLY: projector
        USE charge_types, ONLY: charge_descriptor

! ...   declare subroutine arguments
        LOGICAL   :: tcel
        TYPE (atoms_type),    INTENT(INOUT) :: atoms
        COMPLEX(dbl),         INTENT(INOUT) :: c0(:,:,:,:)
        TYPE (wave_descriptor),  INTENT(IN) :: cdesc
        TYPE (pseudo),        INTENT(INOUT) :: ps
        REAL(dbl) :: rhoe(:,:,:,:)
        COMPLEX(dbl) :: ei1(:,:)
        COMPLEX(dbl) :: ei2(:,:)
        COMPLEX(dbl) :: ei3(:,:)
        COMPLEX(dbl) :: eigr(:,:)
        TYPE (kpoints),       INTENT(IN)    ::  kp
        TYPE (boxdimensions), INTENT(INOUT) ::  ht
        REAL(dbl) :: fi(:,:,:)
        TYPE (projector) :: fnl(:,:)
        TYPE (dft_energy_type) :: edft
        REAL(dbl)    :: vpot(:,:,:,:)
        COMPLEX(dbl), INTENT(IN) :: sfac(:,:)
        LOGICAL, INTENT(IN) :: tforce, tstress, tprint
        REAL(dbl), INTENT(OUT) :: timepre
        TYPE (charge_descriptor),  INTENT(IN) :: desc

        edft%enl = nlrh_m(c0, cdesc, tforce, atoms, fi, kp, fnl, ps%wsg, ps%wnl, eigr)

        CALL rhoofr(kp, c0, cdesc, fi, rhoe, desc, ht)

        CALL vofrhos(tprint, rhoe, desc, tforce, tstress, tforce, atoms, &
          kp, fnl, vpot, ps, c0, cdesc, fi, eigr, ei1, ei2, ei3, sfac, timepre, ht, edft)

        RETURN
      END SUBROUTINE kspotential

!=----------------------------------------------------------------------------=!

!  ----------------------------------------------
!  BEGIN manual

      SUBROUTINE vofrhos(tprint, rhoe, desc, tfor, thdyn, tforce, atoms, &
              kp, fnl, vpot, ps, c0, cdesc, fi, eigr, ei1, ei2, ei3, sfac, timepre, box, edft)

!  this routine computes:
!  ekin = dft_kinetic term of the DFT functional (see dft_kinetic_energy)
!  the one-particle potential v in real space,
!  the total energy etot,
!  the forces fion acting on the ions,
!  the matrix ngdrt used to compute the pair-correlation
!  function gdr for the ions
!  ----------------------------------------------
!  END manual


! ... include modules
      USE cell_module, ONLY: boxdimensions
      USE fft, ONLY : pfwfft, pinvfft
      USE cell_base, ONLY: tpiba2
      USE energies, ONLY: total_energy, dft_energy_type
      USE stress, ONLY: pstress
      USE exchange_correlation, ONLY: exch_corr_energy
      USE funct, ONLY: igcx, igcc
      USE charge_density, ONLY: gradrho
      USE non_local_core_correction, ONLY: add_core_charge, core_charge_forces
      USE chi2, ONLY: rhochi, allocate_chi2, deallocate_chi2
      USE mp_global, ONLY: nproc, mpime, root, group
      USE vanderwaals, ONLY: tvdw, vdw
      USE charge_density, ONLY: checkrho
      USE wave_functions, ONLY: dft_kinetic_energy
      USE brillouin, ONLY: kpoints
      USE mp, ONLY: mp_sum
      USE wave_types, ONLY: wave_descriptor
      USE pseudo_projector, ONLY: projector
      USE atoms_type_module, ONLY: atoms_type
      USE fft_base, ONLY: dfftp
      USE charge_types, ONLY: charge_descriptor
      USE control_flags, ONLY: tscreen, tchi2, iprsta
      USE io_global, ONLY: ionode
      USE cp_types, ONLY: pseudo
      USE io_global, ONLY: stdout
      USE sic_module, ONLY: self_interaction, si_epsilon
      USE reciprocal_vectors, ONLY: gx
      USE gvecp, ONLY: ngm
      USE pseudo_base, ONLY: compute_eself

      IMPLICIT NONE

! ... declare subroutine arguments
      LOGICAL, INTENT(IN) :: tprint, tfor, thdyn, tforce
      REAL(dbl) :: vpot(:,:,:,:)

      TYPE (atoms_type), INTENT(INOUT) :: atoms
      REAL(dbl),    INTENT(IN) :: fi(:,:,:)
      TYPE (projector) :: fnl(:,:)
      TYPE (pseudo),  INTENT(IN) :: ps
      COMPLEX(dbl) :: ei1(:,:)
      COMPLEX(dbl) :: ei2(:,:)
      COMPLEX(dbl) :: ei3(:,:)
      COMPLEX(dbl) :: eigr(:,:)
      TYPE (kpoints), INTENT(IN) :: kp
      COMPLEX(dbl),    INTENT(IN) :: c0(:,:,:,:)
      TYPE (wave_descriptor), INTENT(IN) :: cdesc
      TYPE (charge_descriptor),    INTENT(IN) :: desc
      TYPE (boxdimensions),    INTENT(INOUT) :: box
      TYPE (dft_energy_type) :: edft
      REAL(dbl) :: rhoe(:,:,:,:)
      COMPLEX(dbl), INTENT(IN) :: sfac(:,:)

      TYPE (dft_energy_type) :: edft_self

      REAL(dbl) :: timepre

! ... declare functions
      REAL(dbl)  DDOT

! ... declare other variables

      COMPLEX(dbl), ALLOCATABLE :: vloc(:), self_vloc(:)
      COMPLEX(dbl), ALLOCATABLE :: rho12(:), rhoeg(:,:), self_rhoeg(:)
      COMPLEX(dbl), ALLOCATABLE :: rhoetg(:,:)
      REAL(dbl), ALLOCATABLE :: rhoetr(:,:,:,:)
      REAL(dbl), ALLOCATABLE :: fion_vdw(:,:)
      REAL(dbl), ALLOCATABLE :: grho(:,:,:,:,:)
      REAL(dbl), ALLOCATABLE :: v2xc(:,:,:,:,:)
      REAL(dbl), ALLOCATABLE :: fion(:,:)

      REAL(dbl), ALLOCATABLE :: self_rho(:,:,:,:)
      REAL(dbl), ALLOCATABLE :: self_vpot(:,:,:,:)
      REAL(dbl), ALLOCATABLE :: self_grho(:,:,:,:,:)
      REAL(dbl), ALLOCATABLE :: self_v2xc(:,:,:,:,:)

      REAL(dbl) ::  pail(3,3)

      COMPLEX(dbl) :: self_ehtep,ehtep
      REAL(dbl)    :: self_sxcp, self_vxc

      REAL(dbl)  :: summing1, summing2


      COMPLEX(dbl) :: ehp, eps

      REAL(dbl)  :: dum, sxcp, vxc, ehr, strvxc
      REAL(dbl)  :: omega, desr(6), pesum(16)
      REAL(dbl)  :: s0, s1, s2, s3, s4, s5, s6, s7, s8
      REAL(dbl)  :: rsum( SIZE( rhoe, 4 ) )
      REAL(dbl)  :: zv( atoms%nsp ), raggio( atoms%nsp )

      LOGICAL :: ttstress, ttforce, ttscreen, ttsic, tgc

      INTEGER ig1, ig2, ig3, is, ia, ig, isc, iflag, ispin
      INTEGER ik, i, j, k, isa, idum, nspin
      INTEGER nr1_l, nr2_l, nr3_l, nr1_g, nr2_g, nr3_g,  nnr_l
      INTEGER :: nr1x, nr2x, nr3x
      INTEGER :: ierr

      DATA iflag / 0 /
      SAVE iflag, desr

!  end of declarations
!  ----------------------------------------------

      IF(timing) s0 = cclock()

      nspin = desc % nspin
      
      nr1_l = desc % nxl
      nr2_l = desc % nyl
      nr3_l = desc % nzl

      nr1_g = desc % nx
      nr2_g = desc % ny
      nr3_g = desc % nz

      nnr_l    = nr1_l * nr2_l * nr3_l

      nr1x = dfftp%nr1x
      nr2x = dfftp%nr2x
      nr3x = dfftp%npl

      edft%evdw = 0.0d0
      ttstress = thdyn.OR.(iprsta>2).OR.tprint
      ttforce  = tfor.OR.(iprsta>2).OR.tprint.OR.tforce
      !ttscreen = .TRUE.
      ttscreen = .FALSE.
      ttsic    = ( ABS(self_interaction) /= 0 )
      omega    = box%deth
      tgc      = ( igcx > 0 ) .OR. ( igcc > 0 )


      IF(tchi2) THEN
        CALL allocate_chi2(ngm)
      END IF

      ALLOCATE( rhoetr( nr1x, nr2x, nr3x, nspin) )
      ALLOCATE( fion( 3, atoms%nat ) )

      fion = atoms%for( 1:3, 1:atoms%nat )
      pail = box%pail

      IF(tgc) THEN
        ALLOCATE( grho( nr1x, nr2x, nr3x, 3, nspin ) )
        ALLOCATE( v2xc( nr1x, nr2x, nr3x, nspin, nspin) )
      ELSE
        ALLOCATE( grho( 1, 1, 1, 1, 1 ) )
        ALLOCATE( v2xc( 1, 1, 1, 1, 1 ) )
      END IF

      ALLOCATE( rhoeg(ngm, nspin) )
      ALLOCATE( rhoetg(ngm, nspin) )

      edft%self_sxc = 0.d0
      edft%sxc = 0.d0
      edft%self_ehte = 0.d0
      edft%eht = 0.d0
    
      IF( ttsic ) THEN
      
        IF ( tgc ) THEN

          ALLOCATE(self_grho( nr1x, nr2x, nr3x, 3, nspin ), STAT = ierr)
          IF( ierr /= 0 ) CALL errore(' vofrhos ', ' allocating self_grho ', ierr)
    
          ALLOCATE(self_v2xc( nr1x, nr2x, nr3x, nspin, nspin ), STAT = ierr)
          IF( ierr /= 0 ) CALL errore(' vofrhos ', ' allocating self_v2xc ', ierr)
  
          self_v2xc  = 0.D0

        END IF !on tgc

        ALLOCATE (self_vpot( nr1x, nr2x, nr3x, 2 ), STAT = ierr)
        IF( ierr /= 0 ) CALL errore(' vofrhos ', ' allocating self_vpot ', ierr)
   
        self_vpot  = 0.D0
     
        ALLOCATE (self_rho( nr1x, nr2x, nr3x, 2), STAT = ierr)
        IF( ierr /= 0 ) CALL errore(' vofrhos ', ' allocating self_rho ', ierr)

      END IF !on self_interaction


      IF(timing) s1 = cclock()

! ... compute kinetic energy

      edft%ekin  = 0.0_dbl
      edft%emkin = 0.0_dbl
      edft%ekin  = dft_kinetic_energy(c0, cdesc, kp, fi, rsum, edft%emkin)

      IF(tprint) THEN
        CALL checkrho(rhoe, desc, rsum, omega)
        IF( ionode .AND.  ttscreen ) &
           WRITE( stdout,fmt="(3X,'Using screened Coulomb potential for cluster calculation')")
      END IF

! ... reciprocal-space vectors are in units of alat/(2 pi) so a
! ... multiplicative factor (2 pi/alat)**2 is required
      edft%ekin  = edft%ekin  * tpiba2
      edft%emkin = edft%emkin * tpiba2

      IF( ttstress .OR. ttforce .OR. iflag == 0 )  THEN
        CALL vofesr( edft%esr, desr, fion, ps, atoms, ttstress, box )
        IF( iflag == 0 ) &
          WRITE( stdout, fmt="(/,3X,'ESR (real part of Ewald sum) = ',D16.8,/)" ) edft%esr
        iflag = 1
      END IF

      IF(timing) s2 = cclock()

! ... FFT: rho(r) --> rho(g)  
      DO ispin = 1, nspin

        CALL pfwfft( rhoeg(:,ispin), rhoe(:,:,:,ispin) )


! ...   add core contribution to the charge

        CALL ZCOPY( SIZE(rhoetg,1), rhoeg(1,ispin), 1, rhoetg(1,ispin), 1 )
        CALL DCOPY( SIZE(rhoe(:,:,:,ispin)), rhoe(1,1,1,ispin), 1, rhoetr(1,1,1,ispin), 1 )

        IF( ANY(ps%tnlcc) ) THEN

          ! ...     add core correction
          ! ...     rhoetg = rhoeg + cc
          ! ...     rhoetr = rhoe  + cc
          CALL add_core_charge( rhoetg(:,ispin), rhoetr(:,:,:,ispin), sfac, ps, atoms%nsp )
        ELSE

          ! ...     no core correction
          ! ...     rhoetg = rhoeg
          ! ...     rhoetr = rhoe

          ! ...     chi2
          IF(tchi2) THEN
            IF(nspin.GT.1) CALL errore(' vofrho ',' spin + tchi ',nspin)
            rhochi = rhoeg(:,1)
          END IF

        END IF

        IF(tgc) THEN
          CALL gradrho( rhoetg(:,ispin), grho(:,:,:,:,ispin), gx )
        END IF

      END DO

      IF(timing) s3 = cclock()



! ... Self-interaction correction --- Excor Part
!In any case calculate the Excor part with rhoetr
      
      CALL exch_corr_energy(rhoetr, rhoetg, grho, vpot, sxcp, vxc, v2xc)

      IF( ttsic ) THEN
        IF( ionode ) THEN
          write(stdout,*) 
          write(stdout,*) '  KIND of SELF_INTERACTION CHOOSEN  == ', self_interaction
          write(stdout,*) '  EXC before SIC corr               == ', sxcp * omega / REAL( nr1_g * nr2_g * nr3_g )
          write(stdout,*) 
        END IF
      END IF

      SELECT CASE( ABS(self_interaction) )

        CASE default 

          !  no sic correction
          !
          edft%sxc = sxcp
          edft%self_sxc  = 0.d0
          self_vxc  = 0.d0
                
        CASE(1) 

          !   Delta_Esic to xc = Exc[rhoup-rhodown, 0]
          !
          self_rho(:,:,:,1) = rhoetr(:,:,:,1) - rhoetr(:,:,:,2)
          self_rho(:,:,:,2) = 0.D0
          IF (tgc) THEN
                 self_grho(:,:,:,:,1) = grho(:,:,:,:,1) - grho(:,:,:,:,2)
                 self_grho(:,:,:,:,2) = 0.D0
          ENDIF
          CALL exch_corr_energy(self_rho, rhoetg, self_grho, self_vpot, &
                   self_sxcp, self_vxc, self_v2xc)
          vpot(:,:,:,1) =  vpot(:,:,:,1) - self_vpot(:,:,:,1)
          vpot(:,:,:,2) =  vpot(:,:,:,2) + self_vpot(:,:,:,1)
          IF (tgc) THEN
            v2xc(:,:,:,1,1) =  v2xc(:,:,:,1,1) - self_v2xc(:,:,:,1,1)
            v2xc(:,:,:,2,2) =  v2xc(:,:,:,2,2) + self_v2xc(:,:,:,1,1)
          ENDIF

          ! write(stdout,*)  'da exc la parte da rhodwn rhodwn',self_sxcp

          edft%sxc = sxcp + self_sxcp !- self_sxcp
          vxc = vxc - self_vxc
          edft%self_sxc = self_sxcp * omega / REAL(nr1_g*nr2_g*nr3_g)
          self_vxc = self_vxc * omega / REAL(nr1_g*nr2_g*nr3_g)


        CASE(2) 

          !   Delta_Esic to xc = Exc[rhoup,rhodown] - Exc[rhopaired, rhopaired] 
          !   where Exc[rhoup,rhodown==sxc prevoiusly calculated

          self_rho(:,:,:,1) = rhoetr(:,:,:,2)
          self_rho(:,:,:,2) = rhoetr(:,:,:,2)

          IF (tgc) THEN
            self_grho(:,:,:,:,1) = grho(:,:,:,:,2)
            self_grho(:,:,:,:,2) = grho(:,:,:,:,2)
          ENDIF

!          write(stdout,*)'DA XC SELF_RHO1',self_rho(:,:,:,1)
!          write(stdout,*)'DA XC SELF_RHO2',self_rho(:,:,:,2)

          CALL exch_corr_energy(self_rho, rhoetg, self_grho, self_vpot, &
                 self_sxcp, self_vxc, self_v2xc)

          vpot (:,:,:,2) = self_vpot(:,:,:,2) + self_vpot(:,:,:,1)
          vpot (:,:,:,1) = 0.d0
          write(stdout,*)  'da XC with tgc: E_XC==',self_sxcp

          IF (tgc) THEN
             v2xc(:,:,:,1,1) = 0.d0
             v2xc(:,:,:,2,2) = self_v2xc(:,:,:,2,2) + self_v2xc(:,:,:,1,1)
          ENDIF
          write(stdout,*)  'da exc la parte da rhodwn rhodwn',self_sxcp
          edft%self_sxc= sxcp - self_sxcp
          edft%sxc = self_sxcp !!sxcp - edft%self_sxc
          vxc = self_vxc
          edft%self_sxc = edft%self_sxc * omega / REAL(nr1_g*nr2_g*nr3_g)
          self_vxc = self_vxc * omega / REAL(nr1_g*nr2_g*nr3_g)

      END SELECT 

      !on value of self_interaction for Exchange-Correlation Part 
      !!sxcp is the exchange-correlation part to the energy from LSD with rhoup e rhodown

      IF( ttsic ) THEN
          write(stdout,*) '  Exchange-correlation Energy introducing the SIC'
          write(stdout,*) '  -----------------------------------------------'
          write(stdout,*) '  SXCP from first call     :: ', sxcp * omega / REAL( nr1_g * nr2_g * nr3_g )
          write(stdout,*) '  SXC after SIC-correction :: ', edft%sxc * omega / REAL( nr1_g * nr2_g * nr3_g )
          write(stdout,*) '  D_SIC SIC correction     :: ', edft%self_sxc
          write(stdout,*) '  -----------------------------------------------'
      END IF
                 
      IF ( ttstress ) THEN
        strvxc = ( edft%sxc - vxc ) * omega / REAL( nr1_g * nr2_g * nr3_g )
      END IF

      IF( ANY( ps%tnlcc ) ) THEN
! ...   xc potential (vpot) from real to G space, to compute nlcc forces
! ...   rhoetg = fwfft(vpot)
        DO ispin = 1, nspin
          CALL pfwfft( rhoetg(:,ispin), vpot(:,:,:,ispin) )
        END DO
! ...   now rhoetg contains the xc potential
        IF (ttforce) THEN
          CALL core_charge_forces(fion, rhoetg, ps%rhoc1, ps%tnlcc, atoms, box, ei1, ei2, ei3, kp )
        END IF
      END IF

! ... Van Der Waals energy and forces
      IF (tvdw) THEN
        CALL VdW(edft%evdw, atoms, fion, box)
      END IF

      IF(timing) s4 = cclock()


! ... Calculate hartree potential and energy (eh), and
! ... local part of the pseudopotential and its energy contribution (eps)
! ... Self-interaction correction --- Hartree part

      ALLOCATE( vloc( ngm ) )
      CALL vofloc(ttscreen, ttforce, edft%ehte, edft%ehti, ehp, & 
           eps, vloc, rhoeg, fion, atoms, ps%rhops, ps%vps, kp, eigr, &
           ei1, ei2, ei3, sfac, box, desc, ps%ap )

      !       edft%ehte = REAL ( ehtep )

      edft%self_ehte = 0.d0

      IF ( self_interaction > 0 .AND. nspin == 2 ) THEN

        !  Delta_Esic to Hartree is ALWAYS = -EH(rhoup-rhodw)

        ALLOCATE(self_vloc(ngm), self_rhoeg(ngm), STAT = ierr)
        IF( ierr /= 0 )&
          CALL errore(' vofrhos ', ' allocating self_vloc, self_rhoeg ', ierr)

        self_rhoeg = rhoeg(:,1) - rhoeg(:,2)
        self_vpot = 0.d0

        !  working on the total charge density

        CALL self_vofloc(ttscreen, self_ehtep, self_vloc, self_rhoeg, &
              kp, box, desc)
        CALL pinvfft(self_vpot(:,:,:,1), self_vloc(:))

        self_vpot(:,:,:,1) = si_epsilon * self_vpot(:,:,:,1)
        edft%self_ehte = si_epsilon * REAL(self_ehtep)
 
        vpot(:,:,:,1) =  vpot(:,:,:,1) - self_vpot(:,:,:,1)
        vpot(:,:,:,2) =  vpot(:,:,:,2) + self_vpot(:,:,:,1)

        DEALLOCATE(self_vloc, self_rhoeg)

      END IF

      edft%eh = REAL( ehp ) - edft%self_ehte
     
      IF ( ttsic ) THEN
        IF ( ionode ) THEN
          write(stdout,*) '  Hartree Energy Contribution when SIC is introduced'
          write(stdout,*) '  --------------------------------------------------'
          write(stdout,*) '  HARTREE Potential == ' , REAL( edft%eh )
          write(stdout,*) '  EH(rhoup+rhodwn)  == ' , REAL( ehp )
          write(stdout,*) '  EH(rhoup-rhodwn)  == ' , REAL( edft%self_ehte )
          write(stdout,*) '  --------------------------------------------------'
        END IF
      END IF

      IF( ALLOCATED( self_grho  ) ) DEALLOCATE( self_grho  )
      IF( ALLOCATED( self_v2xc  ) ) DEALLOCATE( self_v2xc  )
      IF( ALLOCATED( self_vpot  ) ) DEALLOCATE( self_vpot )
      IF( ALLOCATED( self_rho   ) ) DEALLOCATE( self_rho  )

! ... vloc(g): hartree and local part of the pseudo potentials (in
! ...          reciprocal space
 
      IF(timing) s5 = cclock()


      DO ispin = 1, nspin

! ...   add hartree end local pseudo potentials ( invfft(vloc) )
! ...   to xc potential (vpot).
! ...   vpot = vpot + invfft(vloc)

        CALL pinvfft( vpot(:,:,:,ispin), vloc(:), 1.0d0 )

      END DO

      DEALLOCATE(vloc)

! ... now potentials are in real space
! ... vpot(r) = hartree + xc + local

      IF(timing) s6 = cclock()

! ... sum up forces
      IF (ttforce) THEN
        CALL mp_sum(fion, group)
      END IF

! ... sum up energies
      CALL mp_sum(eps, group)
      CALL mp_sum(edft%sxc, group)
      CALL mp_sum(edft%self_sxc, group)
      CALL mp_sum(vxc, group)
      CALL mp_sum(edft%eh, group)
      CALL mp_sum(edft%ehte, group)
      CALL mp_sum(edft%ehti, group)
      CALL mp_sum(edft%self_ehte, group)
      CALL mp_sum(edft%ekin, group)
      CALL mp_sum(edft%emkin, group)

      ! ... self interaction energy of the pseudocharges
      do is = 1, atoms%nsp
        zv( is ) = ps%ap(is)%zv
        raggio( is ) = ps%ap(is)%raggio
      end do
      edft%eself = compute_eself( atoms%na, zv, raggio, atoms%nsp ) 


      IF( ttsic .and. ionode ) THEN
        write(stdout,*) 'ESELF_EL::', edft%eself
      END IF

      CALL total_energy(edft,omega,vxc,eps,self_vxc,nr1_g*nr2_g*nr3_g)

! ... compute stress tensor
      IF( ttstress .AND. kp%gamma_only ) THEN
        s8 = cclock()
        CALL pstress( strvxc, rhoeg, rhoetg, pail, desr, fnl, &
          ps, c0, cdesc, fi, eigr, sfac, grho, v2xc, box, edft)
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
        IF(nspin==1) CALL vofmean( sfac, ps%rhops, rhoeg(:,1) )
        IF(nspin==2) THEN
             ALLOCATE(rho12(ngm))
            rho12 (:) = rhoeg(:,1)+rhoeg(:,2)
          CALL vofmean( sfac, ps%rhops, rho12)
           DEALLOCATE(rho12)
        END IF
      END IF

      DEALLOCATE( rhoeg, rhoetg, rhoetr, grho, v2xc, fion )

      IF(tchi2) THEN
        CALL deallocate_chi2
      END IF

      IF(iprsta>2) THEN
        CALL memstat(mpime)
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
      CALL cpflush

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

  SUBROUTINE cluster_bc( screen_coul, hg, box, desc )

      USE green_functions, ONLY: greenf
      USE mp_global, ONLY: mpime
      USE fft, ONLY : pfwfft
      USE fft_base, ONLY: dfftp
      USE charge_types, ONLY: charge_descriptor
      USE processors_grid_module, ONLY: get_grid_info
      USE cell_module, ONLY: boxdimensions, s_to_r, alat
      USE constants, ONLY: gsmall, pi
      USE cell_base, ONLY: tpiba2
      
      REAL(dbl), INTENT(IN) :: hg(:)
      TYPE (boxdimensions),    INTENT(IN) :: box
      TYPE (charge_descriptor),    INTENT(IN) :: desc
      COMPLEX(dbl) :: screen_coul(:)

! ... declare external function
      REAL(dbl) :: erf, erfc
      EXTERNAL erf, erfc

! ... Locals
      REAL(dbl), ALLOCATABLE :: grr(:,:,:)
      COMPLEX(dbl), ALLOCATABLE :: grg(:)
      REAL(dbl) :: rc, r(3), s(3), rmod, g2, rc2, arg, omega, fact
      INTEGER   :: ig, i, j, k
      INTEGER   :: nr1_l, nr2_l, nr3_l, nr1_g, nr2_g, nr3_g
      INTEGER   :: ir1, ir2, ir3


      nr1_l = desc % nxl
      nr2_l = desc % nyl
      nr3_l = desc % nzl
      nr1_g = desc % nx
      nr2_g = desc % ny
      nr3_g = desc % nz

      ir1 = 1
      ir2 = 1
      ir3 = 1
      DO k = 1, mpime
        ir3 = ir3 + dfftp%npp( k )
      END DO

      ALLOCATE( grr( dfftp%nr1x, dfftp%nr2x, dfftp%npl ) )
      ALLOCATE( grg( SIZE( screen_coul ) ) )

! ... Martina and Tuckerman convergence criterium
      rc  = 7.0d0 / alat
      rc2 = rc**2
      omega = box%deth
      fact  = omega / ( nr1_g * nr2_g * nr3_g )
      IF( MOD(nr1_g * nr2_g * nr3_g, 2) /= 0 ) fact = -fact

      DO k = 1, nr3_l
        s(3) = REAL ( (k-1) + (ir3 - 1) ) / nr3_g - 0.5d0
        DO j = 1, nr2_l
          s(2) = REAL ( (j-1) + (ir2 - 1) ) / nr2_g - 0.5d0
          DO i = 1, nr1_l
            s(1) = REAL ( (i-1) + (ir1 - 1) ) / nr1_g - 0.5d0
            CALL S_TO_R(S, R, box)
            rmod = SQRT( r(1)**2 + r(2)**2 + r(3)**2 )
            IF( rmod < gsmall ) THEN
              grr(i,j,k) = fact * 2.0d0 * rc / SQRT( pi )
            ELSE
              grr(i,j,k) = fact * erf( rc * rmod ) / rmod
            END IF
          END DO
        END DO
      END DO

      CALL pfwfft( grg, grr ) ! grg = FFT( grr )

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

      SUBROUTINE vofloc(tscreen, ttforce, ehte, ehti, eh, eps, vloc, rhoeg, &
                 fion, atoms, rhops, vps, kp, eigr, ei1, ei2, ei3, sfac, ht, desc, ap)

!  this routine computes:
!  omega = ht%deth
!  rho_e(ig)    =  (sum over ispin) rhoeg(ig,ispin) 
!  rho_I(ig)    =  (sum over is) sfac(is,ig) * rhops(ig,is) 
!  vloc_h(ig)   =  fpi / ( g(ig) * tpiba2 ) * { rho_e(ig) + rho_I(ig) }
!  vloc_ps(ig)  =  (sum over is) sfac(is,ig) * vps(ig,is)
!
!  Eps = Fact * omega * (sum over ig) CMPLX( rho_e(ig) ) * vloc_ps(ig)
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
!      Fact * omega * ( sum over ig, ispin) (tx_h(ig,is) + tx_ps(ig,is)) * 
!      gx(ig) * eigrx(ig,isa) * eigry(ig,isa) * eigrz(ig,isa) 
!  if Gamma symmetry Fact = 2.0 else Fact = 1
!
!  ----------------------------------------------
!  END manual

      USE constants, ONLY: fpi
      USE cell_base, ONLY: tpiba2, tpiba
      USE cell_module, ONLY: boxdimensions
      USE brillouin, ONLY: kpoints
      USE charge_types, ONLY: charge_descriptor
      USE atoms_type_module, ONLY: atoms_type
      USE pseudo_types, ONLY: pseudo_ncpp
      USE io_global, ONLY: stdout
      USE grid_dimensions, ONLY: nr1, nr2, nr3
      USE reciprocal_vectors, ONLY: mill_l
      USE ions_base, ONLY: nat
      USE gvecp, ONLY: ngm
      USE reciprocal_vectors, ONLY: gstart, gx, g

      IMPLICIT NONE

! ... Arguments

      TYPE (atoms_type) :: atoms
      TYPE (pseudo_ncpp), INTENT(IN) :: ap(:)
      TYPE (kpoints), INTENT(in) :: kp
      TYPE (boxdimensions), INTENT(in) :: ht
      TYPE (charge_descriptor), INTENT(IN) :: desc
      LOGICAL      :: ttforce
      LOGICAL      :: tscreen
      REAL(dbl)    :: fion(:,:)
      REAL(dbl)    :: rhops(:,:), vps(:,:)
      COMPLEX(dbl) :: vloc(:)
      COMPLEX(dbl) :: rhoeg(:,:)
      COMPLEX(dbl), INTENT(IN) :: sfac(:,:)
      REAL(dbl)    :: ehte, ehti
      COMPLEX(dbl) :: eh, eps
      COMPLEX(dbl) :: ei1(-nr1:nr1,nat)
      COMPLEX(dbl) :: ei2(-nr2:nr2,nat)
      COMPLEX(dbl) :: ei3(-nr3:nr3,nat)
      COMPLEX(dbl) :: eigr(:,:)

! ... Locals

      INTEGER      :: is, ia, isa, ig, ig1, ig2, ig3, nspin, ispin
      REAL(dbl)    :: fpibg, cost, omega
      COMPLEX(dbl) :: cxc, rhet, rhog, vp, rp, gxc, gyc, gzc
      COMPLEX(dbl) :: teigr, cnvg, cvn, tx, ty, tz, vscreen
      COMPLEX(dbl), ALLOCATABLE :: ftmp(:,:)
      COMPLEX(dbl), ALLOCATABLE :: screen_coul(:)

! ... Subroutine body ...

      nspin = SIZE(rhoeg,2)
      IF(TTFORCE) THEN
        ALLOCATE( ftmp(3, SIZE(fion, 2) ) )
        ftmp = CMPLX(0.0_dbl,0.0_dbl)
      END IF

      IF( tscreen ) THEN
        ALLOCATE( screen_coul( ngm ) )
        CALL cluster_bc( screen_coul, g, ht, desc )
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
        ehte     = ehte +  fpibg *   REAL(rhet * CONJG(rhet))
        ehti     = ehti +  fpibg *   REAL(  rp * CONJG(rp))

        IF(TTFORCE) THEN
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
      IF(TTFORCE) THEN
! ...   each processor add its own contribution to the array FION
        IF(kp%gamma_only) THEN
          cost = 2.D0 * ht%deth * tpiba
        ELSE
          cost = ht%deth * tpiba
        END IF
        FION = FION + REAL(ftmp) * cost
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
        ehte    = ehte +  vscreen *   REAL(rhet * CONJG(rhet))
        ehti    = ehti +  vscreen *   REAL(  rp * CONJG(rp))
        DO ispin = 1, nspin
          IF( kp%gamma_only ) THEN
            eps = eps + vp * CONJG(RHOEG(1,ispin)) * 0.5d0
          ELSE
            eps = eps + vp * CONJG(RHOEG(1,ispin))
          END IF
        END DO
      END IF
! ...
      IF( .NOT. kp%gamma_only ) THEN
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
   SUBROUTINE vofesr(esr, desr, fion, ps, atoms, tstress, ht)
!=----------------------------------------------------------------------------=!

      USE constants, ONLY: sqrtpm1
      USE cell_module, ONLY: s_to_r, boxdimensions, pbcs
      USE mp_global, ONLY: nproc, mpime, group
      USE mp, ONLY: mp_sum
      USE parallel_types, ONLY: BLOCK_PARTITION_SHAPE
      USE descriptors_module, ONLY: global_index, local_dimension
      USE atoms_type_module, ONLY: atoms_type
      USE cp_types, ONLY: pseudo
 
      IMPLICIT NONE

! ... ARGUMENTS 
      
      TYPE (atoms_type) :: atoms
      TYPE (pseudo)     :: ps
      REAL(dbl) :: ESR
      REAL(dbl) :: DESR(:)
      REAL(dbl) :: FION(:,:)
      LOGICAL   :: TSTRESS
      TYPE (boxdimensions), INTENT(in) :: ht

! ... declare external function
      REAL(dbl) :: erf, erfc
      EXTERNAL erf, erfc

      
! ... LOCALS 

      INTEGER :: na_loc, ia_s, ia_e, igis
      INTEGER :: k, i, j, l, m, is, ia, infm, ix, iy, iz, ishft
      INTEGER :: npt, isa, me
      INTEGER :: iakl, iajm
      LOGICAL :: split, tzero, tshift
      INTEGER, ALLOCATABLE   :: iatom(:,:)
      REAL(dbl), ALLOCATABLE :: zv2(:,:)
      REAL(dbl), ALLOCATABLE :: rc(:,:)  
      REAL(dbl), ALLOCATABLE :: fionloc(:,:,:) 
      REAL(dbl)  :: RXLM(3), SXLM(3)
      REAL(dbl)  :: xlm, ylm, zlm, erre2, rlm, arg, esrtzero
      REAL(dbl)  :: addesr, addpre, repand, fxx
      REAL(dbl)  :: rckj_m1
      REAL(dbl)  :: zvk, zvj, zv2_kj
      REAL(dbl)  :: fact_pre

      INTEGER, DIMENSION(6), PARAMETER :: ALPHA = (/ 1,2,3,2,3,3 /)
      INTEGER, DIMENSION(6), PARAMETER :: BETA  = (/ 1,1,1,2,2,3 /)

! ... SUBROUTINE BODY 

      me = mpime + 1

      !  Here count the pairs of atoms

      npt = 0
      DO k = 1, atoms%nsp
        DO j = k, atoms%nsp
          DO l = 1, atoms%na(k)
            IF ( k == j ) THEN
              infm = l             ! If the specie is the same avoid  
            ELSE                   ! atoms double counting
              infm = 1
            END IF
            DO m = infm, atoms%na(j)
              npt = npt + 1
            END DO
          END DO
        END DO
      END DO

      ALLOCATE( iatom( 4, npt ) )
      ALLOCATE( rc( atoms%nsp, atoms%nsp ) )
      ALLOCATE( zv2( atoms%nsp, atoms%nsp ) )
      ALLOCATE( fionloc( 3, atoms%nax, atoms%nsp ) )
      rc      = 0.0_dbl
      zv2     = 0.0_dbl
      fionloc = 0.0_dbl

      !  Here pre-compute some factors

      DO k = 1, atoms%nsp
        DO j = k, atoms%nsp
          zv2( k, j ) = ps%ap(k)%zv * ps%ap(j)%zv 
          rc ( k, j ) = SQRT( ps%ap(k)%raggio**2 + ps%ap(j)%raggio**2 )
        END DO
      END DO

      !  Here store the indexes of all pairs of atoms

      npt = 0
      DO k = 1, atoms%nsp
        DO j = k, atoms%nsp
          DO l = 1, atoms%na(k)
            IF (k.EQ.j) THEN
              infm = l
            ELSE
              infm = 1
            END IF
            DO m = infm, atoms%na(j)
              npt = npt + 1
              iatom(1,npt) = k
              iatom(2,npt) = j
              iatom(3,npt) = l
              iatom(4,npt) = m
            END DO
          END DO
        END DO
      END DO

      xlm     = 1.0_dbl
      ylm     = 1.0_dbl
      zlm     = 1.0_dbl
      ESR     = 0.0_dbl
      DESR    = 0.0_dbl

      ! NA_LOC = LOCALDIM(npt,NPROC,ME)
      NA_LOC = local_dimension( npt, 1, mpime, 0, nproc, BLOCK_PARTITION_SHAPE)
      ! IA_S   = GLOBALINDEX(1,npt,NPROC,ME)
      IA_S   = global_index( 1, npt, 1, mpime, 0, nproc, BLOCK_PARTITION_SHAPE )
      IA_E   = IA_S + NA_LOC - 1

      DO ia = ia_s, ia_e

        k = iatom(1,ia)
        j = iatom(2,ia)
        l = iatom(3,ia)
        m = iatom(4,ia)

        zv2_kj   = zv2(k,j)
        rckj_m1  = 1.0_dbl / rc(k,j)
        fact_pre = (2.0_dbl * zv2_kj * sqrtpm1) * rckj_m1

        IF( (l.EQ.m) .AND. (k.EQ.j)) THEN      
! ...     same atoms
          xlm=0.0_dbl; ylm=0.0_dbl; zlm=0.0_dbl; tzero=.TRUE.
        ELSE
! ...     different atoms
          iakl = atoms%isa(k) + l - 1
          iajm = atoms%isa(j) + m - 1
          xlm = atoms%taus(1,iakl) - atoms%taus(1,iajm)
          ylm = atoms%taus(2,iakl) - atoms%taus(2,iajm)
          zlm = atoms%taus(3,iakl) - atoms%taus(3,iajm)
          CALL pbcs(xlm,ylm,zlm,xlm,ylm,zlm,1)
          TZERO=.FALSE.
        END IF

        DO IX=-IESR,IESR
          SXLM(1) = XLM + REAL(IX)
          DO IY=-IESR,IESR
            SXLM(2) = YLM + REAL(IY)
            DO IZ=-IESR,IESR
              TSHIFT= IX.EQ.0 .AND. IY.EQ.0 .AND. IZ.EQ.0
              IF(.NOT.(TZERO.AND.TSHIFT)) THEN
                SXLM(3) = ZLM + REAL(IZ)
                CALL S_TO_R(SXLM,RXLM,ht)
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
                DO I=1,3
                  FXX = REPAND*RXLM(I)
                  FIONLOC(I,L,K) = FIONLOC(I,L,K) + FXX
                  FIONLOC(I,M,J) = FIONLOC(I,M,J) - FXX
                END DO
                IF(TSTRESS) THEN
                  DO I=1,6
                    FXX = REPAND * RXLM(ALPHA(I))*RXLM(BETA(I))
                    DESR(I) = DESR(I) - FXX
                  END DO
                END IF
              END IF
            END DO    ! IZ
          END DO      ! IY
        END DO        ! IX
      END DO

!
!     each processor add its own contribution to the array FION
!
      isa = 0
      DO IS = 1, atoms%nsp
        DO IA = 1, atoms%na(is)
          isa = isa + 1
          FION(1,ISA) = FION(1,ISA)+FIONLOC(1,IA,IS)
          FION(2,ISA) = FION(2,ISA)+FIONLOC(2,IA,IS)
          FION(3,ISA) = FION(3,ISA)+FIONLOC(3,IA,IS)
        END DO
      END DO

      CALL mp_sum(esr, group)
     
      DEALLOCATE(iatom)
      DEALLOCATE(rc)
      DEALLOCATE(zv2)
      DEALLOCATE(fionloc)
      
      RETURN
!=----------------------------------------------------------------------------=!
   END SUBROUTINE vofesr
!=----------------------------------------------------------------------------=!


!  ----------------------------------------------
!  BEGIN manual

      SUBROUTINE self_vofloc(tscreen, ehte, vloc, rhoeg, kp, ht, desc)

!  adds the hartree part of the self interaction
!
!  ----------------------------------------------
!  END manual

      USE constants
      USE cell_module, ONLY: boxdimensions
      USE brillouin, ONLY: kpoints
      USE charge_types, ONLY: charge_descriptor
      USE cell_base, ONLY: tpiba2
      USE gvecp, ONLY: ngm
      USE reciprocal_vectors, ONLY: gstart, g

      IMPLICIT NONE

! ... Arguments
      TYPE (kpoints), INTENT(in) :: kp
      TYPE (boxdimensions), INTENT(in) :: ht
      TYPE (charge_descriptor), INTENT(IN) :: desc
      LOGICAL      :: tscreen
      COMPLEX(dbl) :: vloc(:)
      COMPLEX(dbl) :: rhoeg(:)
      COMPLEX(dbl) :: ehte

! ... Locals

      INTEGER      :: ig
      REAL(dbl)    :: fpibg, omega
      COMPLEX(dbl) :: rhog
      COMPLEX(dbl) :: vscreen
      COMPLEX(dbl), ALLOCATABLE :: screen_coul(:)

! ... Subroutine body ...


      IF( tscreen ) THEN
        ALLOCATE( screen_coul( ngm ) )
        CALL cluster_bc( screen_coul, g, ht, desc )
      END IF

!=======================================================================
!==  HARTREE AND LOCAL PART OF THE PSEUDO ENERGIES                    ==
!=======================================================================

      omega = ht%deth
      ehte = 0.D0

      DO IG = gstart, ngm

        rhog  = rhoeg(ig)
        IF( tscreen ) THEN
          FPIBG     = fpi / ( g(ig) * tpiba2 ) + screen_coul(ig)
        ELSE
          FPIBG     = fpi / ( g(ig) * tpiba2 )
        END IF

        vloc(ig) = fpibg * rhog
        ehte       = ehte   +  fpibg *   rhog * CONJG(rhog)

      END DO
!!!! always for a comparison with the el-Hartree contribution from lda  

 
! ... G = 0 element
      IF ( gstart == 2 ) THEN
        IF( tscreen ) THEN
          vscreen = screen_coul(1)
        ELSE
          vscreen = 0.0d0
        END IF
        rhog    = rhoeg(1)
        vloc(1) = vscreen * rhog
        ehte      = ehte   +  vscreen *  rhog * CONJG(rhog)
      END IF
! ...
      IF( .NOT. kp%gamma_only ) THEN
        ehte  = ehte  * 0.5d0
      END IF
      ehte =        ehte * omega

! ...
      IF(ALLOCATED(screen_coul)) DEALLOCATE(screen_coul)

      RETURN
      END SUBROUTINE self_vofloc



      SUBROUTINE localisation( wfc, atoms_m, kp, ht, desc)

!  adds the hartree part of the self interaction
!
!  ----------------------------------------------
!  END manual

      USE constants
      USE cell_module, ONLY: boxdimensions, s_to_r
      USE brillouin, ONLY: kpoints
      USE charge_types, ONLY: charge_descriptor
      USE atoms_type_module, ONLY: atoms_type
      USE fft, ONLY : pw_invfft, pfwfft, pinvfft
      USE sic_module, ONLY: ind_localisation, nat_localisation, print_localisation
      USE sic_module, ONLY: rad_localisation, pos_localisation
      USE ions_base, ONLY: ind_srt
      USE fft_base, ONLY: dfftp
      USE cell_base, ONLY: tpiba2
      USE reciprocal_vectors, ONLY: gstart, g
      USE gvecp, ONLY: ngm

      IMPLICIT NONE

! ... Arguments

      COMPLEX(dbl), INTENT(IN) :: wfc(:)
      TYPE (atoms_type), INTENT(in) :: atoms_m
      TYPE (kpoints), INTENT(in) :: kp
      TYPE (boxdimensions), INTENT(in) :: ht
      TYPE (charge_descriptor), INTENT(IN) :: desc


! ... Locals

      REAL(dbl)    :: ehte
      INTEGER      :: ig, at, ia, is, isa_input, isa_sorted, isa_loc
      REAL(dbl)    :: fpibg, omega, aRe, aR2, R(3)
      INTEGER      :: Xmin, Ymin, Zmin, Xmax, Ymax, Zmax
      INTEGER      :: nr1_l, nr2_l, nr3_l
      REAL(dbl)    :: work, work2
      COMPLEX(dbl) :: rhog
      REAL(dbl), ALLOCATABLE :: density(:,:,:), psi(:,:,:)
      COMPLEX(dbl), ALLOCATABLE :: k_density(:), cpsi(:,:,:)
      COMPLEX(dbl) :: vscreen
      COMPLEX(dbl), ALLOCATABLE :: screen_coul(:)
      INTEGER :: nr1x, nr2x, nr3x

! ... Subroutine body ...


      IF( .FALSE. ) THEN
        ALLOCATE( screen_coul( ngm ) )
        CALL cluster_bc( screen_coul, g, ht, desc )
      END IF


      nr1x = dfftp%nr1x
      nr2x = dfftp%nr2x
      nr3x = dfftp%npl

      omega = ht%deth

      nr1_l = desc % nxl
      nr2_l = desc % nyl
      nr3_l = desc % nzl

      ALLOCATE( density( nr1x, nr2x, nr3x ) )
      ALLOCATE( psi( nr1x, nr2x, nr3x ) )
      ALLOCATE( cpsi( nr1x, nr2x, nr3x ) )
      ALLOCATE( k_density( ngm ) )

      CALL pw_invfft( cpsi(:,:,:), wfc(:), wfc(:) )
      psi = REAL( cpsi, dbl )
      DEALLOCATE( cpsi )

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

            WRITE(6,*) 'ATOM ', ind_localisation( isa_input )
            WRITE(6,*) 'POS  ', atoms_m%taus( :, isa_sorted )

            work  = nr1_l
            work2 = rad_localisation * work
            work  = work * R(1) - work2
            Xmin  = FLOOR(work)
            work  = work + 2*work2
            Xmax  = FLOOR(work)
            IF ( Xmax > nr1_l ) Xmax = nr1_l
            IF ( Xmin < 1 ) Xmin = 1
            work  = nr2_l
            work2 = rad_localisation * work
            work  = work * R(2) - work2
            Ymin  = FLOOR(work)
            work  = work + 2*work2
            Ymax  = FLOOR(work)
            IF ( Ymax > nr2_l ) Ymax = nr2_l
            IF ( Ymin < 1 ) Ymin = 1
            work  = nr3_l
            work2 = rad_localisation * work
            work  = work * R(3) - work2
            Zmin  = FLOOR(work)
            work  = work + 2*work2
            Zmax  = FLOOR(work)
            IF ( Zmax > nr3_l ) Zmax = nr3_l
            IF ( Zmin < 1 ) Zmin = 1

            density = 0.D0
            density( Xmin:Xmax, Ymin:Ymax, Zmin:Zmax ) = &
            psi( Xmin:Xmax, Ymin:Ymax, Zmin:Zmax ) * psi( Xmin:Xmax, Ymin:Ymax, Zmin:Zmax )
            CALL pfwfft( k_density, density )

! ... G /= 0 elements

            DO IG = gstart, ngm

              rhog  = k_density(ig)
              IF( .FALSE. ) THEN
                FPIBG     = fpi / ( g(ig) * tpiba2 ) + screen_coul(ig)
              ELSE
                FPIBG     = fpi / ( g(ig) * tpiba2 )
              END IF

              ehte       = ehte   +  fpibg *   REAL(rhog * CONJG(rhog))

            END DO

! ... G = 0 element

            IF ( gstart == 2 ) THEN
              IF( .FALSE. ) THEN
                vscreen = screen_coul(1)
              ELSE
                vscreen = 0.0d0
              END IF
              rhog    = k_density(1)
              ehte    = ehte   +  vscreen *  REAL(rhog * CONJG(rhog))
            END IF
! ...
            IF( .NOT. kp%gamma_only ) THEN
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
