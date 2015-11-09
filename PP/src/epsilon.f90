
! Copyright (C) 2004-2009 Andrea Benassi and Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .

!------------------------------
 MODULE grid_module
!------------------------------
  USE kinds,        ONLY : DP
  IMPLICIT NONE
  PRIVATE

  !
  ! general purpose vars
  !
  REAL(DP), ALLOCATABLE  :: focc(:,:), wgrid(:)
  REAL(DP)               :: alpha
  !
  !
  PUBLIC :: grid_build, grid_destroy
  PUBLIC :: focc, wgrid, alpha
  !
CONTAINS

!---------------------------------------------
  SUBROUTINE grid_build(nw, wmax, wmin)
  !-------------------------------------------
  !
  USE kinds,     ONLY : DP
  USE wvfct,     ONLY : nbnd, wg
  USE klist,     ONLY : nks, wk, nelec
  USE lsda_mod,  ONLY : nspin
  USE uspp,      ONLY : okvan
  !
  IMPLICIT NONE
  !
  ! input vars
  INTEGER,  INTENT(in) :: nw
  REAL(DP), INTENT(in) :: wmax ,wmin
  !
  ! local vars
  INTEGER         :: iw,ik,i,ierr

  !
  ! check on the number of bands: we need to include empty bands in order to allow
  ! to write the transitions
  !
  IF ( REAL(nbnd, DP)  <= nelec / 2.0_DP ) CALL errore('epsilon', 'bad band number', 1)

  !
  ! spin is not implemented
  !
  IF( nspin > 2 ) CALL errore('grid_build','Non collinear spin  calculation not implemented',1)

  !
  ! USPP are not implemented (dipole matrix elements are not trivial at all)
  !
  IF ( okvan ) CALL errore('grid_build','USPP are not implemented',1)

  ALLOCATE ( focc( nbnd, nks), STAT=ierr )
  IF (ierr/=0) CALL errore('grid_build','allocating focc', abs(ierr))
  !
  ALLOCATE( wgrid( nw ), STAT=ierr )
  IF (ierr/=0) CALL errore('grid_build','allocating wgrid', abs(ierr))

  !
  ! check on k point weights, no symmetry operations are allowed
  !
  DO ik = 2, nks
     !
     IF ( abs( wk(1) - wk(ik) ) > 1.0d-8 ) &
        CALL errore('grid_build','non unifrom kpt grid', ik )
     !
  ENDDO
  !
  ! occupation numbers, to be normalized differently
  ! whether we are spin resolved or not
  !
  IF(nspin==1) THEN
    DO ik = 1,nks
    DO i  = 1,nbnd
         focc(i,ik)= wg(i, ik ) * 2.0_DP / wk( ik )
    ENDDO
    ENDDO
  ELSEIF(nspin==2) THEN
    DO ik = 1,nks
    DO i  = 1,nbnd
         focc(i,ik)= wg(i, ik ) * 1.0_DP / wk( ik )
    ENDDO
    ENDDO
  ENDIF
  !
  ! set the energy grid
  !
  alpha = (wmax - wmin) / REAL(nw-1, KIND=DP)
  !
  DO iw = 1, nw
      wgrid(iw) = wmin + (iw-1) * alpha
  ENDDO
  !
END SUBROUTINE grid_build
!
!
!----------------------------------
  SUBROUTINE grid_destroy
  !----------------------------------
  IMPLICIT NONE
  INTEGER :: ierr
  !
  IF ( allocated( focc) ) THEN
      !
      DEALLOCATE ( focc, wgrid, STAT=ierr)
      CALL errore('grid_destroy','deallocating grid stuff',abs(ierr))
      !
  ENDIF
  !
END SUBROUTINE grid_destroy

END MODULE grid_module


!------------------------------
PROGRAM epsilon
!------------------------------
  !
  ! Compute the complex macroscopic dielectric function,
  ! at the RPA level, neglecting local field effects.
  ! Eps is computed both on the real or immaginary axis
  !
  ! Authors: Andrea Benassi, Andrea Ferretti, Carlo Cavazzoni
  !
  ! NOTE: Part of the basic implementation is taken from pw2gw.f90;
  !
  !
  USE kinds,       ONLY : DP
  USE io_global,   ONLY : stdout, ionode, ionode_id
  USE mp,          ONLY : mp_bcast
  USE mp_world,    ONLY : world_comm
  USE iotk_module
  USE xml_io_base
  USE io_files,    ONLY : tmp_dir, prefix, outdir
  USE constants,   ONLY : RYTOEV
  USE ener,        ONLY : ef
  USE klist,       ONLY : lgauss
  USE ktetra,      ONLY : ltetra
  USE wvfct,       ONLY : nbnd
  USE lsda_mod,    ONLY : nspin
  USE mp_global,   ONLY : mp_startup
  USE environment, ONLY : environment_start, environment_end
  !
  IMPLICIT NONE
  !
  CHARACTER(LEN=256), EXTERNAL :: trimcheck
  !
  ! input variables
  !
  INTEGER                 :: nw,nbndmin,nbndmax
  REAL(DP)                :: intersmear,intrasmear,wmax,wmin,shift
  CHARACTER(10)           :: calculation,smeartype
  LOGICAL                 :: metalcalc
  !
  NAMELIST / inputpp / prefix, outdir, calculation
  NAMELIST / energy_grid / smeartype,intersmear,intrasmear,wmax,wmin,nbndmin,nbndmax,nw,shift
  !
  ! local variables
  !
  INTEGER :: ios

!---------------------------------------------
! program body
!---------------------------------------------
!
  ! initialise environment
  !
#ifdef __MPI
  CALL mp_startup ( )
#endif
  CALL environment_start ( 'epsilon' )
  !
  ! Set default values for variables in namelist
  !
  calculation  = 'eps'
  prefix       = 'pwscf'
  shift        = 0.0d0
  outdir       = './'
  intersmear   = 0.136
  wmin         = 0.0d0
  wmax         = 30.0d0
  nbndmin      = 1
  nbndmax      = 0
  nw           = 600
  smeartype    = 'gauss'
  intrasmear   = 0.0d0
  metalcalc    = .false.

  !
  ! this routine allows the user to redirect the input using -input
  ! instead of <
  !
  CALL input_from_file( )

  !
  ! read input file
  !
  IF (ionode) WRITE( stdout, "( 2/, 5x, 'Reading input file...' ) " )
  ios = 0
  !
  IF ( ionode ) READ (5, inputpp, IOSTAT=ios)
  !
  CALL mp_bcast ( ios, ionode_id, world_comm )
  IF (ios/=0) CALL errore('epsilon', 'reading namelist INPUTPP', abs(ios))
  !
  IF ( ionode ) THEN
     !
     READ (5, energy_grid, IOSTAT=ios)
     !
     tmp_dir = trimcheck(outdir)
     !
  ENDIF
  !
  CALL mp_bcast ( ios, ionode_id, world_comm )
  IF (ios/=0) CALL errore('epsilon', 'reading namelist ENERGY_GRID', abs(ios))
  !
  ! ... Broadcast variables
  !
  IF (ionode) WRITE( stdout, "( 5x, 'Broadcasting variables...' ) " )

  CALL mp_bcast( smeartype, ionode_id, world_comm )
  CALL mp_bcast( calculation, ionode_id, world_comm )
  CALL mp_bcast( prefix, ionode_id, world_comm )
  CALL mp_bcast( tmp_dir, ionode_id, world_comm )
  CALL mp_bcast( shift, ionode_id, world_comm )
  CALL mp_bcast( outdir, ionode_id, world_comm )
  CALL mp_bcast( intrasmear, ionode_id, world_comm )
  CALL mp_bcast( intersmear, ionode_id, world_comm)
  CALL mp_bcast( wmax, ionode_id, world_comm )
  CALL mp_bcast( wmin, ionode_id, world_comm )
  CALL mp_bcast( nw, ionode_id, world_comm )
  CALL mp_bcast( nbndmin, ionode_id, world_comm )
  CALL mp_bcast( nbndmax, ionode_id, world_comm )

  !
  ! read PW simulation parameters from prefix.save/data-file.xml
  !
  IF (ionode) WRITE( stdout, "( 5x, 'Reading PW restart file...' ) " )

  CALL read_file
  CALL openfil_pp

  !
  ! few conversions
  !

  IF (ionode) WRITE(stdout,"(2/, 5x, 'Fermi energy [eV] is: ',f8.5)") ef *RYTOEV

  IF (lgauss .or. ltetra) THEN
      metalcalc=.true.
      IF (ionode) WRITE( stdout, "( 5x, 'The system is a metal...' ) " )
  ELSE
      IF (ionode) WRITE( stdout, "( 5x, 'The system is a dielectric...' ) " )
  ENDIF

  IF (nbndmax == 0) nbndmax = nbnd

  !
  ! ... run the specific pp calculation
  !
  IF (ionode) WRITE(stdout,"(/, 5x, 'Performing ',a,' calculation...')") trim(calculation)
  CALL start_clock(trim(calculation))
  SELECT CASE ( trim(calculation) )
  !
  CASE ( 'eps' )
      !
      CALL eps_calc ( intersmear,intrasmear,nw,wmax,wmin,nbndmin,nbndmax,shift,metalcalc,nspin )
      !
  CASE ( 'jdos' )
      !
      CALL jdos_calc ( smeartype,intersmear,nw,wmax,wmin,nbndmin,nbndmax,shift,nspin )
      !
  CASE ( 'offdiag' )
      !
      CALL offdiag_calc ( intersmear,intrasmear,nw,wmax,wmin,nbndmin,nbndmax,shift,metalcalc,nspin )
      !
  CASE ( 'occ' )
      !
      CALL occ_calc ()
      !
  CASE DEFAULT
      !
      CALL errore('epsilon','invalid CALCULATION = '//trim(calculation),1)
      !
  END SELECT
  !
  CALL stop_clock(trim(calculation))
  IF ( ionode ) WRITE( stdout , "(/)" )
  CALL print_clock( trim(calculation) )
  CALL print_clock( 'dipole_calc' )
  IF ( ionode ) WRITE( stdout, *  )
  !
  CALL environment_end ( 'epsilon' )
  !
  CALL stop_pp ()

END PROGRAM epsilon


!-----------------------------------------------------------------------------
SUBROUTINE eps_calc ( intersmear,intrasmear, nw, wmax, wmin, nbndmin, nbndmax, shift, &
                      metalcalc , nspin)
  !-----------------------------------------------------------------------------
  !
  USE kinds,                ONLY : DP
  USE constants,            ONLY : PI, RYTOEV
  USE cell_base,            ONLY : tpiba2, omega
  USE wvfct,                ONLY : nbnd, et
  USE ener,                 ONLY : efermi => ef
  USE klist,                ONLY : nks, nkstot, degauss
  USE io_global,            ONLY : ionode, stdout
  !
  USE grid_module,          ONLY : alpha, focc, wgrid, grid_build, grid_destroy
  USE mp_global,            ONLY : inter_pool_comm
  USE mp,                   ONLY : mp_sum
  !
  IMPLICIT NONE

  !
  ! input variables
  !
  INTEGER,         INTENT(in) :: nw,nbndmin,nbndmax,nspin
  REAL(DP),        INTENT(in) :: wmax, wmin, intersmear,intrasmear, shift
  LOGICAL,         INTENT(in) :: metalcalc
  !
  ! local variables
  !
  INTEGER       :: i, ik, iband1, iband2,is
  INTEGER       :: iw, iwp, ierr
  REAL(DP)      :: etrans, const, w, renorm(3)
  !
  REAL(DP), ALLOCATABLE    :: epsr(:,:), epsi(:,:), epsrc(:,:,:), epsic(:,:,:)
  REAL(DP), ALLOCATABLE    :: ieps(:,:), eels(:,:), iepsc(:,:,:), eelsc(:,:,:)
  REAL(DP), ALLOCATABLE    :: dipole(:,:,:)
  COMPLEX(DP),ALLOCATABLE  :: dipole_aux(:,:,:)
!
!--------------------------
! main routine body
!--------------------------
!
  !
  ! perform some consistency checks, calculate occupation numbers and setup w grid
  !
  CALL grid_build(nw, wmax, wmin)
  !
  ! allocate main spectral and auxiliary quantities
  !
  ALLOCATE( dipole(3, nbnd, nbnd), STAT=ierr )
  IF (ierr/=0) CALL errore('epsilon','allocating dipole', abs(ierr) )
  !
  ALLOCATE( dipole_aux(3, nbnd, nbnd), STAT=ierr )
  IF (ierr/=0) CALL errore('epsilon','allocating dipole_aux', abs(ierr) )

!
! spin unresolved calculation
!
IF (nspin == 1) THEN
  !
  ALLOCATE( epsr( 3, nw), epsi( 3, nw), eels( 3, nw), ieps(3,nw ), STAT=ierr )
  IF (ierr/=0) CALL errore('epsilon','allocating eps', abs(ierr))

  !
  ! initialize response functions
  !
  epsr(:,:)  = 0.0_DP
  epsi(:,:)  = 0.0_DP
  ieps(:,:)  = 0.0_DP

  !
  ! main kpt loop
  !
  kpt_loop: &
  DO ik = 1, nks
     !
     ! For every single k-point: order k+G for
     !                           read and distribute wavefunctions
     !                           compute dipole matrix 3 x nbnd x nbnd parallel over g
     !                           recover g parallelism getting the total dipole matrix
     !
     CALL dipole_calc( ik, dipole_aux, metalcalc , nbndmin, nbndmax)
     !
     dipole(:,:,:)= tpiba2 * REAL( dipole_aux(:,:,:) * conjg(dipole_aux(:,:,:)), DP )

     !
     ! Calculation of real and immaginary parts
     ! of the macroscopic dielettric function from dipole
     ! approximation.
     ! 'intersmear' is the brodening parameter
     !
     !Interband
     !
     DO iband2 = nbndmin,nbndmax
         !
         IF ( focc(iband2,ik) < 2.0d0) THEN
     DO iband1 = nbndmin,nbndmax
         !
         IF (iband1==iband2) CYCLE
         IF ( focc(iband1,ik) >= 1e-4 ) THEN
         IF (abs(focc(iband2,ik)-focc(iband1,ik))< 1e-3) CYCLE
               !
               ! transition energy
               !
               etrans = ( et(iband2,ik) -et(iband1,ik) ) * RYTOEV + shift
               !
               ! loop over frequencies
               !
               DO iw = 1, nw
                   !
                   w = wgrid(iw)
                   !
                   epsi(:,iw) = epsi(:,iw) + dipole(:,iband1,iband2) * intersmear * w* &
                                             RYTOEV**3 * (focc(iband1,ik))/  &
                                  (( (etrans**2 -w**2 )**2 + intersmear**2 * w**2 )* etrans )

                   epsr(:,iw) = epsr(:,iw) + dipole(:,iband1,iband2) * RYTOEV**3 * &
                                             (focc(iband1,ik)) * &
                                             (etrans**2 - w**2 ) / &
                                  (( (etrans**2 -w**2 )**2 + intersmear**2 * w**2 )* etrans )
               ENDDO

         ENDIF
     ENDDO
         ENDIF
     ENDDO
     !
     !Intraband (only if metalcalc is true)
     !
     IF (metalcalc) THEN
     DO iband1 = nbndmin,nbndmax
         !
         IF ( focc(iband1,ik) < 2.0d0) THEN
         IF ( focc(iband1,ik) >= 1e-4 ) THEN
               !
               ! loop over frequencies
               !
               DO iw = 1, nw
                   !
                   w = wgrid(iw)
                   !
                  epsi(:,iw) = epsi(:,iw) +  dipole(:,iband1,iband1) * intrasmear * w* &
                                RYTOEV**2 * (exp((et(iband1,ik)-efermi)/degauss ))/  &
                    (( w**4 + intrasmear**2 * w**2 )*(1+exp((et(iband1,ik)-efermi)/ &
                    degauss))**2*degauss )

                  epsr(:,iw) = epsr(:,iw) - dipole(:,iband1,iband1) * RYTOEV**2 * &
                                            (exp((et(iband1,ik)-efermi)/degauss )) * w**2 / &
                    (( w**4 + intrasmear**2 * w**2 )*(1+exp((et(iband1,ik)-efermi)/ &
                    degauss))**2*degauss )
               ENDDO

         ENDIF
         ENDIF

     ENDDO
     ENDIF
  ENDDO kpt_loop

  !
  ! recover over kpt parallelization (inter_pool)
  !
  CALL mp_sum( epsr, inter_pool_comm )
  CALL mp_sum( epsi, inter_pool_comm )

  !
  ! impose the correct normalization
  !
  const = 64.0d0 * PI / ( omega * REAL(nkstot, DP) )
  epsr(:,:) = 1.0_DP + epsr(:,:) * const
  epsi(:,:) =          epsi(:,:) * const

  !
  ! Calculation of eels spectrum
  !
  DO iw = 1, nw
      !
      eels(:,iw) = epsi(:,iw) / ( epsr(:,iw)**2 + epsi(:,iw)**2 )
      !
  ENDDO

  !
  !  calculation of dielectric function on the immaginary frequency axe
  !

  DO iw = 1, nw
  DO iwp = 2, nw
      !
      ieps(:,iw) = ieps(:,iw) + wgrid(iwp) * epsi(:,iwp) / ( wgrid(iwp)**2 + wgrid(iw)**2)
      !
  ENDDO
  ENDDO

  ieps(:,:) = 1.0d0 + 2 / PI * ieps(:,:) * alpha

  !
  ! check  dielectric function  normalizzation via sumrule
  !
 DO i=1,3
     renorm(i) = alpha * sum( epsi(i,:) * wgrid(:) )
 ENDDO
  !
  IF ( ionode ) THEN
      !
      WRITE(stdout,"(/,5x, 'The bulk xx plasmon frequency [eV] is: ',f15.9 )")  sqrt(renorm(1) * 2.0d0 / PI)
      WRITE(stdout,"(5x, 'The bulk yy plasmon frequency [eV] is: ',f15.9 )")  sqrt(renorm(2) * 2.0d0 / PI)
      WRITE(stdout,"(5x, 'The bulk zz plasmon frequency [eV] is: ',f15.9 )")  sqrt(renorm(3) * 2.0d0 / PI)
      WRITE(stdout,"(/,5x, 'Writing output on file...' )")
      !
      ! write results on data files
      !

      OPEN (30, FILE='epsr.dat', FORM='FORMATTED' )
      OPEN (40, FILE='epsi.dat', FORM='FORMATTED' )
      OPEN (41, FILE='eels.dat', FORM='FORMATTED' )
      OPEN (42, FILE='ieps.dat', FORM='FORMATTED' )
      !
      WRITE(30, "(2x,'# energy grid [eV]     epsr_x  epsr_y  epsr_z')" )
      WRITE(40, "(2x,'# energy grid [eV]     epsi_x  epsi_y  epsi_z')" )
      WRITE(41, "(2x,'# energy grid [eV]  eels components [arbitrary units]')" )
      WRITE(42, "(2x,'# energy grid [eV]     ieps_x  ieps_y  ieps_z ')" )
      !
      DO iw =1, nw
          !
          WRITE(30,"(4f15.6)") wgrid(iw), epsr(1:3, iw)
          WRITE(40,"(4f15.6)") wgrid(iw), epsi(1:3, iw)
          WRITE(41,"(4f15.6)") wgrid(iw), eels(1:3, iw)
          WRITE(42,"(4f15.6)") wgrid(iw), ieps(1:3, iw)
          !
      ENDDO
      !
      CLOSE(30)
      CLOSE(40)
      CLOSE(41)
      CLOSE(42)
      !
  ENDIF

  DEALLOCATE ( epsr, epsi, eels, ieps)
!
! collinear spin calculation
!
ELSEIF (nspin == 2 ) THEN
  !
  ALLOCATE( epsrc( 0:1, 3, nw), epsic( 0:1,3, nw), eelsc( 0:1,3, nw), iepsc(0:1,3,nw ), STAT=ierr )
  IF (ierr/=0) CALL errore('epsilon','allocating eps', abs(ierr))

  !
  ! initialize response functions
  !
  epsrc(:,:,:)  = 0.0_DP
  epsic(:,:,:)  = 0.0_DP
  iepsc(:,:,:)  = 0.0_DP

  !
  ! main kpt loop
  !

spin_loop: &
DO is=0,1
  kpt_loopspin: &
! if nspin=2 the number of nks must be even (even if the calculation
! is performed at gamma point only), so nks must be always a multiple of 2
  DO ik = 1 + is * int(nks/2), int(nks/2) +  is * int(nks/2)
     !
     ! For every single k-point: order k+G for
     !                           read and distribute wavefunctions
     !                           compute dipole matrix 3 x nbnd x nbnd parallel over g
     !                           recover g parallelism getting the total dipole matrix
     !
     CALL dipole_calc( ik, dipole_aux, metalcalc , nbndmin, nbndmax)
     !
     dipole(:,:,:)= tpiba2 * REAL( dipole_aux(:,:,:) * conjg(dipole_aux(:,:,:)), DP )

     !
     ! Calculation of real and immaginary parts
     ! of the macroscopic dielettric function from dipole
     ! approximation.
     ! 'intersmear' is the brodening parameter
     !
     !Interband
     !
     DO iband2 = nbndmin,nbndmax
         !
         IF ( focc(iband2,ik) < 1.0d0) THEN
     DO iband1 = nbndmin,nbndmax
         !
         IF (iband1==iband2) CYCLE
         IF ( focc(iband1,ik) >= 1e-4 ) THEN
         IF (abs(focc(iband2,ik)-focc(iband1,ik))< 1e-3) CYCLE
               !
               ! transition energy
               !
               etrans = ( et(iband2,ik) -et(iband1,ik) ) * RYTOEV + shift
               !
               ! loop over frequencies
               !
               DO iw = 1, nw
                   !
                   w = wgrid(iw)
                   !
                   epsic(is,:,iw) = epsic(is,:,iw) + dipole(:,iband1,iband2) * intersmear * w* &
                                             RYTOEV**3 * (focc(iband1,ik))/  &
                                  (( (etrans**2 -w**2 )**2 + intersmear**2 * w**2 )* etrans )

                   epsrc(is,:,iw) = epsrc(is,:,iw) + dipole(:,iband1,iband2) * RYTOEV**3 * &
                                             (focc(iband1,ik)) * &
                                             (etrans**2 - w**2 ) / &
                                  (( (etrans**2 -w**2 )**2 + intersmear**2 * w**2 )* etrans )
               ENDDO

         ENDIF
     ENDDO
         ENDIF
     ENDDO
     !
     !Intraband (only if metalcalc is true)
     !
     IF (metalcalc) THEN
     DO iband1 = nbndmin,nbndmax
         !
         IF ( focc(iband1,ik) < 1.0d0) THEN
         IF ( focc(iband1,ik) >= 1e-4 ) THEN
               !
               ! loop over frequencies
               !
               DO iw = 1, nw
                   !
                   w = wgrid(iw)
                   !
                  epsic(is,:,iw) = epsic(is,:,iw) +  dipole(:,iband1,iband1) * intrasmear * w* &
                                RYTOEV**2 * (exp((et(iband1,ik)-efermi)/degauss ))/  &
                    (( w**4 + intrasmear**2 * w**2 )*(1+exp((et(iband1,ik)-efermi)/ &
                    degauss))**2*degauss )

                  epsrc(is,:,iw) = epsrc(is,:,iw) - dipole(:,iband1,iband1) * RYTOEV**2 * &
                                            (exp((et(iband1,ik)-efermi)/degauss )) * w**2 / &
                    (( w**4 + intrasmear**2 * w**2 )*(1+exp((et(iband1,ik)-efermi)/ &
                    degauss))**2*degauss )
               ENDDO

         ENDIF
         ENDIF

     ENDDO
     ENDIF
  ENDDO kpt_loopspin
ENDDO spin_loop
  !
  ! recover over kpt parallelization (inter_pool)
  !
  CALL mp_sum( epsr, inter_pool_comm )
  CALL mp_sum( epsi, inter_pool_comm )

  !
  ! impose the correct normalization
  !
  const = 128.0d0 * PI / ( omega * REAL(nkstot, DP) )
  epsrc(:,:,:) = 1.0_DP + epsrc(:,:,:) * const
  epsic(:,:,:) =          epsic(:,:,:) * const

  !
  ! Calculation of eels spectrum
  !
  DO iw = 1, nw
      !
      eelsc(:,:,iw) = epsic(:,:,iw) / ( epsrc(:,:,iw)**2 + epsic(:,:,iw)**2 )
      !
  ENDDO

  !
  !  calculation of dielectric function on the immaginary frequency axe
  !

  DO iw = 1, nw
  DO iwp = 2, nw
      !
      iepsc(:,:,iw) = iepsc(:,:,iw) + wgrid(iwp) * epsic(:,:,iwp) / ( wgrid(iwp)**2 + wgrid(iw)**2)
      !
  ENDDO
  ENDDO

  iepsc(:,:,:) = 1.0d0 + 2.0_DP / PI * iepsc(:,:,:) * alpha

  IF (ionode) THEN
      WRITE(stdout,"(/,5x, 'Writing output on file...' )")
      !
      ! write results on data files
      !

      OPEN (30, FILE='uepsr.dat', FORM='FORMATTED' )
      OPEN (40, FILE='uepsi.dat', FORM='FORMATTED' )
      OPEN (41, FILE='ueels.dat', FORM='FORMATTED' )
      OPEN (42, FILE='uieps.dat', FORM='FORMATTED' )
      OPEN (43, FILE='depsr.dat', FORM='FORMATTED' )
      OPEN (44, FILE='depsi.dat', FORM='FORMATTED' )
      OPEN (45, FILE='deels.dat', FORM='FORMATTED' )
      OPEN (46, FILE='dieps.dat', FORM='FORMATTED' )
      OPEN (47, FILE='epsr.dat', FORM='FORMATTED' )
      OPEN (48, FILE='epsi.dat', FORM='FORMATTED' )
      OPEN (49, FILE='eels.dat', FORM='FORMATTED' )
      OPEN (50, FILE='ieps.dat', FORM='FORMATTED' )
      !
      WRITE(30, "(2x,'# energy grid [eV]     epsr_x  epsr_y  epsr_z')" )
      WRITE(40, "(2x,'# energy grid [eV]     epsi_x  epsi_y  epsi_z')" )
      WRITE(41, "(2x,'# energy grid [eV]  eels components [arbitrary units]')" )
      WRITE(42, "(2x,'# energy grid [eV]     ieps_x  ieps_y  ieps_z ')" )
      WRITE(43, "(2x,'# energy grid [eV]     epsr_x  epsr_y  epsr_z')" )
      WRITE(44, "(2x,'# energy grid [eV]     epsi_x  epsi_y  epsi_z')" )
      WRITE(45, "(2x,'# energy grid [eV]  eels components [arbitrary units]')" )
      WRITE(46, "(2x,'# energy grid [eV]     ieps_x  ieps_y  ieps_z ')" )
      WRITE(47, "(2x,'# energy grid [eV]     epsr_x  epsr_y  epsr_z')" )
      WRITE(48, "(2x,'# energy grid [eV]     epsi_x  epsi_y  epsi_z')" )
      WRITE(49, "(2x,'# energy grid [eV]  eels components [arbitrary units]')" )
      WRITE(50, "(2x,'# energy grid [eV]     ieps_x  ieps_y  ieps_z ')" )
      !
      DO iw =1, nw
          !
          WRITE(30,"(4f15.6)") wgrid(iw), epsrc(0,1:3, iw)
          WRITE(40,"(4f15.6)") wgrid(iw), epsic(0,1:3, iw)
          WRITE(41,"(4f15.6)") wgrid(iw), eelsc(0,1:3, iw)
          WRITE(42,"(4f15.6)") wgrid(iw), iepsc(0,1:3, iw)
          WRITE(43,"(4f15.6)") wgrid(iw), epsrc(1,1:3, iw)
          WRITE(44,"(4f15.6)") wgrid(iw), epsic(1,1:3, iw)
          WRITE(45,"(4f15.6)") wgrid(iw), eelsc(1,1:3, iw)
          WRITE(46,"(4f15.6)") wgrid(iw), iepsc(1,1:3, iw)
          WRITE(47,"(4f15.6)") wgrid(iw), epsrc(1,1:3, iw)+epsrc(0,1:3, iw)
          WRITE(48,"(4f15.6)") wgrid(iw), epsic(1,1:3, iw)+epsic(0,1:3, iw)
          WRITE(49,"(4f15.6)") wgrid(iw), eelsc(1,1:3, iw)+eelsc(0,1:3, iw)
          WRITE(50,"(4f15.6)") wgrid(iw), iepsc(1,1:3, iw)+iepsc(0,1:3, iw)
          !
      ENDDO
      !
      CLOSE(30)
      CLOSE(40)
      CLOSE(41)
      CLOSE(42)
      CLOSE(43)
      CLOSE(44)
      CLOSE(45)
      CLOSE(46)
      CLOSE(47)
      CLOSE(48)
      CLOSE(49)
      CLOSE(50)
      !
  ENDIF
  DEALLOCATE ( epsrc, epsic, eelsc, iepsc)
ENDIF
  !
  ! local cleaning
  !
  CALL grid_destroy()
  !
  DEALLOCATE (  dipole, dipole_aux )

END SUBROUTINE eps_calc

!----------------------------------------------------------------------------------------
SUBROUTINE jdos_calc ( smeartype,intersmear,nw,wmax,wmin,nbndmin,nbndmax,shift,nspin )
  !--------------------------------------------------------------------------------------
  !
  USE kinds,                ONLY : DP
  USE constants,            ONLY : PI, RYTOEV
  USE wvfct,                ONLY : nbnd, et
  USE klist,                ONLY : nks
  USE io_global,            ONLY : ionode, stdout
  USE grid_module,          ONLY : alpha, focc, wgrid, grid_build, grid_destroy
  !
  IMPLICIT NONE

  !
  ! input variables
  !
  INTEGER,         INTENT(in) :: nw,nbndmin,nbndmax,nspin
  REAL(DP),        INTENT(in) :: wmax, wmin, intersmear, shift
  CHARACTER(*),    INTENT(in) :: smeartype
  !
  ! local variables
  !
  INTEGER       :: ik, is, iband1, iband2
  INTEGER       :: iw, ierr
  REAL(DP)      :: etrans, w, renorm, count, srcount(0:1), renormzero,renormuno
  !
  REAL(DP), ALLOCATABLE    :: jdos(:),srjdos(:,:)
!
!--------------------------
! main routine body
!--------------------------
!
! No wavefunctions are needed in order to compute jdos, only eigenvalues,
! they are distributed to each task so
! no mpi calls are necessary in this routine
  !
  ! perform some consistency checks, calculate occupation numbers and setup w grid
  !
  CALL grid_build(nw, wmax, wmin )

!
! spin unresolved calculation
!
IF (nspin == 1) THEN
  !
  ! allocate main spectral and auxiliary quantities
  !
  ALLOCATE( jdos(nw), STAT=ierr )
      IF (ierr/=0) CALL errore('epsilon','allocating jdos',abs(ierr))
  !
  ! initialize jdos
  !
  jdos(:)=0.0_DP

  ! Initialising a counter for the number of transition
  count=0.0_DP

  !
  ! main kpt loop
  !

  IF (smeartype=='lorentz') THEN

    kpt_lor: &
    DO ik = 1, nks
       !
       ! Calculation of joint density of states
       ! 'intersmear' is the brodening parameter
       !
       DO iband2 = 1,nbnd
           IF ( focc(iband2,ik) <  2.0d0) THEN
       DO iband1 = 1,nbnd
           !
           IF ( focc(iband1,ik) >= 1.0d-4 ) THEN
                 !
                 ! transition energy
                 !
                 etrans = ( et(iband2,ik) -et(iband1,ik) ) * RYTOEV  + shift
                 !
                 IF( etrans < 1.0d-10 ) CYCLE

                 count = count + (focc(iband1,ik)-focc(iband2,ik))
                 !
                 ! loop over frequencies
                 !
                 DO iw = 1, nw
                     !
                     w = wgrid(iw)
                     !
                     jdos(iw) = jdos(iw) + intersmear * (focc(iband1,ik)-focc(iband2,ik)) &
                                  / ( PI * ( (etrans -w )**2 + (intersmear)**2 ) )

                 ENDDO

           ENDIF
       ENDDO
           ENDIF
       ENDDO

    ENDDO kpt_lor

  ELSEIF (smeartype=='gauss') THEN

    kpt_gauss: &
    DO ik = 1, nks

       !
       ! Calculation of joint density of states
       ! 'intersmear' is the brodening parameter
       !
       DO iband2 = 1,nbnd
       DO iband1 = 1,nbnd
           !
           IF ( focc(iband2,ik) <  2.0d0) THEN
           IF ( focc(iband1,ik) >= 1.0d-4 ) THEN
                 !
                 ! transition energy
                 !
                 etrans = ( et(iband2,ik) -et(iband1,ik) ) * RYTOEV  + shift
                 !
                 IF( etrans < 1.0d-10 ) CYCLE

                 ! loop over frequencies
                 !

                 count=count+ (focc(iband1,ik)-focc(iband2,ik))

                 DO iw = 1, nw
                     !
                     w = wgrid(iw)
                     !
                     jdos(iw) = jdos(iw) + (focc(iband1,ik)-focc(iband2,ik)) * &
                                exp(-(etrans-w)**2/intersmear**2) &
                                  / (intersmear * sqrt(PI))

                 ENDDO

           ENDIF
           ENDIF
       ENDDO
       ENDDO

    ENDDO kpt_gauss

  ELSE

    CALL errore('epsilon', 'invalid SMEARTYPE = '//trim(smeartype), 1)

  ENDIF

  !
  ! jdos normalizzation
  !

  jdos(:)=jdos(:)/count

  !
  ! check jdos normalization
  !

  renorm = alpha * sum( jdos(:) )
  !
  ! write results on data files
  !
  IF (ionode) THEN
     WRITE(stdout,"(/,5x, 'Integration over JDOS gives: ',f15.9,' instead of 1.0d0' )") renorm
     WRITE(stdout,"(/,5x, 'Writing output on file...' )")

     OPEN (30, FILE='jdos.dat', FORM='FORMATTED' )
     !
     WRITE(30, "(2x,'# energy grid [eV]     JDOS [1/eV] ')" )
     !
     DO iw =1, nw
         !
         WRITE(30,"(4f15.6)") wgrid(iw), jdos(iw)
         !
     ENDDO
     !
     CLOSE(30)
  ENDIF
  !
  ! local cleaning
  !
  DEALLOCATE ( jdos )

!
! collinear spin calculation
!
ELSEIF(nspin==2) THEN
  !
  ! allocate main spectral and auxiliary quantities
  !
  ALLOCATE( srjdos(0:1,nw), STAT=ierr )
      IF (ierr/=0) CALL errore('epsilon','allocating spin resolved jdos',abs(ierr))
  !
  ! initialize jdos
  !
  srjdos(:,:)=0.0_DP

  ! Initialising a counter for the number of transition
  srcount(:)=0.0_DP

  !
  ! main kpt loop
  !

  IF (smeartype=='lorentz') THEN

  DO is=0,1
    ! if nspin=2 the number of nks must be even (even if the calculation
    ! is performed at gamma point only), so nks must be always a multiple of 2
    DO ik = 1 + is * int(nks/2), int(nks/2) +  is * int(nks/2)
       !
       ! Calculation of joint density of states
       ! 'intersmear' is the brodening parameter
       !
       DO iband2 = 1,nbnd
           IF ( focc(iband2,ik) <  2.0d0) THEN
       DO iband1 = 1,nbnd
           !
           IF ( focc(iband1,ik) >= 1.0d-4 ) THEN
                 !
                 ! transition energy
                 !
                 etrans = ( et(iband2,ik) -et(iband1,ik) ) * RYTOEV  + shift
                 !
                 IF( etrans < 1.0d-10 ) CYCLE

                 ! loop over frequencies
                 !
                 srcount(is)=srcount(is)+ (focc(iband1,ik)-focc(iband2,ik))

                 DO iw = 1, nw
                     !
                     w = wgrid(iw)
                     !
                     srjdos(is,iw) = srjdos(is,iw) + intersmear * (focc(iband1,ik)-focc(iband2,ik)) &
                                  / ( PI * ( (etrans -w )**2 + (intersmear)**2 ) )

                 ENDDO

           ENDIF
       ENDDO
           ENDIF
       ENDDO

    ENDDO
 ENDDO

  ELSEIF (smeartype=='gauss') THEN

  DO is=0,1
    ! if nspin=2 the number of nks must be even (even if the calculation
    ! is performed at gamma point only), so nks must be always a multiple of 2
    DO ik = 1 + is * int(nks/2), int(nks/2) +  is * int(nks/2)
       !
       ! Calculation of joint density of states
       ! 'intersmear' is the brodening parameter
       !
       DO iband2 = 1,nbnd
       DO iband1 = 1,nbnd
           !
           IF ( focc(iband2,ik) <  2.0d0) THEN
           IF ( focc(iband1,ik) >= 1.0d-4 ) THEN
                 !
                 ! transition energy
                 !
                 etrans = ( et(iband2,ik) -et(iband1,ik) ) * RYTOEV  + shift
                 !
                 IF( etrans < 1.0d-10 ) CYCLE

                 ! loop over frequencies
                 !

                 srcount(is)=srcount(is)+ (focc(iband1,ik)-focc(iband2,ik))

                 DO iw = 1, nw
                     !
                     w = wgrid(iw)
                     !
                     srjdos(is,iw) = srjdos(is,iw) + (focc(iband1,ik)-focc(iband2,ik)) * &
                                exp(-(etrans-w)**2/intersmear**2) &
                                  / (intersmear * sqrt(PI))

                 ENDDO

           ENDIF
           ENDIF
       ENDDO
       ENDDO

    ENDDO
 ENDDO

  ELSE

    CALL errore('epsilon', 'invalid SMEARTYPE = '//trim(smeartype), 1)

  ENDIF

  !
  ! jdos normalizzation
  !
  DO is = 0,1
    srjdos(is,:)=srjdos(is,:)/srcount(is)
  ENDDO
  !
  ! check jdos normalization
  !

  renormzero = alpha * sum( srjdos(0,:) )
  renormuno = alpha * sum( srjdos(1,:) )
  !
  ! write results on data files
  !
  IF (ionode) THEN
     WRITE(stdout,"(/,5x, 'Integration over spin UP JDOS gives: ',f15.9,' instead of 1.0d0' )") renormzero
     WRITE(stdout,"(/,5x, 'Integration over spin DOWN JDOS gives: ',f15.9,' instead of 1.0d0' )") renormuno
     WRITE(stdout,"(/,5x, 'Writing output on file...' )")

     OPEN (30, FILE='jdos.dat', FORM='FORMATTED' )
     !
     WRITE(30, "(2x,'# energy grid [eV]     UJDOS [1/eV]      DJDOS[1:eV]')" )
     !
     DO iw =1, nw
         !
         WRITE(30,"(4f15.6)") wgrid(iw), srjdos(0,iw), srjdos(1,iw)
         !
     ENDDO
     !
     CLOSE(30)
  ENDIF

  DEALLOCATE ( srjdos )
ENDIF
  !
  ! local cleaning
  !
  CALL grid_destroy()

END SUBROUTINE jdos_calc

!-----------------------------------------------------------------------------
SUBROUTINE offdiag_calc ( intersmear,intrasmear, nw, wmax, wmin, nbndmin, nbndmax,&
                          shift, metalcalc, nspin )
  !-----------------------------------------------------------------------------
  !
  USE kinds,                ONLY : DP
  USE constants,            ONLY : PI, RYTOEV
  USE cell_base,            ONLY : tpiba2, omega
  USE wvfct,                ONLY : nbnd, et
  USE ener,                 ONLY : efermi => ef
  USE klist,                ONLY : nks, nkstot, degauss
  USE grid_module,          ONLY : focc, wgrid, grid_build, grid_destroy
  USE io_global,            ONLY : ionode, stdout
  USE mp_global,            ONLY : inter_pool_comm
  USE mp,                   ONLY : mp_sum

  !
  IMPLICIT NONE

  !
  ! input variables
  !
  INTEGER,         INTENT(in) :: nw,nbndmin,nbndmax,nspin
  REAL(DP),        INTENT(in) :: wmax, wmin, intersmear,intrasmear, shift
  LOGICAL,         INTENT(in) :: metalcalc
  !
  ! local variables
  !
  INTEGER       :: ik, iband1, iband2
  INTEGER       :: iw, ierr, it1, it2
  REAL(DP)      :: etrans, const, w
  !
  COMPLEX(DP), ALLOCATABLE :: dipole_aux(:,:,:)
  COMPLEX(DP), ALLOCATABLE :: epstot(:,:,:),dipoletot(:,:,:,:)
  !
  !--------------------------
  ! main routine body
  !--------------------------
  !
  ! perform some consistency checks, calculate occupation numbers and setup w grid
  !
  CALL grid_build(nw, wmax, wmin )
  !
  ! allocate main spectral and auxiliary quantities
  !
  ALLOCATE( dipoletot(3,3, nbnd, nbnd), STAT=ierr )
  IF (ierr/=0) CALL errore('epsilon','allocating dipoletot', abs(ierr) )
  !
  ALLOCATE( dipole_aux(3, nbnd, nbnd), STAT=ierr )
  IF (ierr/=0) CALL errore('epsilon','allocating dipole_aux', abs(ierr) )
  !
  ALLOCATE(epstot( 3,3, nw),STAT=ierr )
  IF (ierr/=0) CALL errore('epsilon','allocating epstot', abs(ierr))

   !
   ! initialize response functions
   !
   epstot  = (0.0_DP,0.0_DP)
   !
   ! main kpt loop
   !
   DO ik = 1, nks
     !
     ! For every single k-point: order k+G for
     !                           read and distribute wavefunctions
     !                           compute dipole matrix 3 x nbnd x nbnd parallel over g
     !                           recover g parallelism getting the total dipole matrix
     !
     CALL dipole_calc( ik, dipole_aux, metalcalc, nbndmin, nbndmax)
     !
     DO it2 = 1, 3
        DO it1 = 1, 3
           dipoletot(it1,it2,:,:) = tpiba2 * dipole_aux(it1,:,:) * conjg( dipole_aux(it2,:,:) )
        ENDDO
     ENDDO
     !
     ! Calculation of real and immaginary parts
     ! of the macroscopic dielettric function from dipole
     ! approximation.
     ! 'intersmear' is the brodening parameter
     !
     DO iband2 = 1,nbnd
         IF ( focc(iband2,ik) <  2.0d0) THEN
     DO iband1 = 1,nbnd
         !
         IF ( focc(iband1,ik) >= 1e-4 ) THEN
             !
             ! transition energy
             !
             etrans = ( et(iband2,ik) -et(iband1,ik) ) * RYTOEV + shift
             !
             IF (abs(focc(iband2,ik)-focc(iband1,ik))< 1e-4) CYCLE
             !
             ! loop over frequencies
             !
             DO iw = 1, nw
                  !
                  w = wgrid(iw)
                  !
                  epstot(:,:,iw) = epstot(:,:,iw) + dipoletot(:,:,iband1,iband2)*RYTOEV**3/(etrans) *&
                                   focc(iband1,ik)/(etrans**2 - w**2 - (0,1)*intersmear*w)
             ENDDO
             !
         ENDIF
     ENDDO
         ENDIF
     ENDDO
     !
     !Intraband (only if metalcalc is true)
     !
     IF (metalcalc) THEN
     DO iband1 = 1,nbnd
         !
         IF ( focc(iband1,ik) < 2.0d0) THEN
         IF ( focc(iband1,ik) >= 1e-4 ) THEN
               !
               ! loop over frequencies
               !
               DO iw = 1, nw
                   !
                   w = wgrid(iw)
                   !
                   epstot(:,:,iw) = epstot(:,:,iw) - dipoletot(:,:,iband1,iband1)* &
                                RYTOEV**2 * (exp((et(iband1,ik)-efermi)/degauss ))/  &
                    (( w**2 + (0,1)*intrasmear*w)*(1+exp((et(iband1,ik)-efermi)/ &
                    degauss))**2*degauss )
               ENDDO

         ENDIF
         ENDIF

     ENDDO
     ENDIF
  ENDDO

  !
  ! recover over kpt parallelization (inter_pool)
  !
  CALL mp_sum( epstot, inter_pool_comm )
  !
  ! impose the correct normalization
  !
  const = 64.0d0 * PI / ( omega * REAL(nkstot, DP) )
  epstot(:,:,:) = epstot(:,:,:) * const
  !
  ! add diagonal term
  !
  epstot(1,1,:) = 1.0_DP + epstot(1,1,:)
  epstot(2,2,:) = 1.0_DP + epstot(2,2,:)
  epstot(3,3,:) = 1.0_DP + epstot(3,3,:)
  !
  ! write results on data files
  !
  IF (ionode) THEN
      !
      WRITE(stdout,"(/,5x, 'Writing output on file...' )")
      !
      OPEN (41, FILE='epsxx.dat', FORM='FORMATTED' )
      OPEN (42, FILE='epsxy.dat', FORM='FORMATTED' )
      OPEN (43, FILE='epsxz.dat', FORM='FORMATTED' )
      OPEN (44, FILE='epsyx.dat', FORM='FORMATTED' )
      OPEN (45, FILE='epsyy.dat', FORM='FORMATTED' )
      OPEN (46, FILE='epsyz.dat', FORM='FORMATTED' )
      OPEN (47, FILE='epszx.dat', FORM='FORMATTED' )
      OPEN (48, FILE='epszy.dat', FORM='FORMATTED' )
      OPEN (49, FILE='epszz.dat', FORM='FORMATTED' )
      !
      WRITE(41, "(2x,'# energy grid [eV]     epsr     epsi')" )
      WRITE(42, "(2x,'# energy grid [eV]     epsr     epsi')" )
      WRITE(43, "(2x,'# energy grid [eV]     epsr     epsi')" )
      WRITE(44, "(2x,'# energy grid [eV]     epsr     epsi')" )
      WRITE(45, "(2x,'# energy grid [eV]     epsr     epsi')" )
      WRITE(46, "(2x,'# energy grid [eV]     epsr     epsi')" )
      WRITE(47, "(2x,'# energy grid [eV]     epsr     epsi')" )
      WRITE(48, "(2x,'# energy grid [eV]     epsr     epsi')" )
      WRITE(49, "(2x,'# energy grid [eV]     epsr     epsi')" )
      !
      DO iw =1, nw
         !
         WRITE(41,"(4f15.6)") wgrid(iw), REAL(epstot(1,1, iw)), aimag(epstot(1,1, iw))
         WRITE(42,"(4f15.6)") wgrid(iw), REAL(epstot(1,2, iw)), aimag(epstot(1,2, iw))
         WRITE(43,"(4f15.6)") wgrid(iw), REAL(epstot(1,3, iw)), aimag(epstot(1,3, iw))
         WRITE(44,"(4f15.6)") wgrid(iw), REAL(epstot(2,1, iw)), aimag(epstot(2,1, iw))
         WRITE(45,"(4f15.6)") wgrid(iw), REAL(epstot(2,2, iw)), aimag(epstot(2,2, iw))
         WRITE(46,"(4f15.6)") wgrid(iw), REAL(epstot(2,3, iw)), aimag(epstot(2,3, iw))
         WRITE(47,"(4f15.6)") wgrid(iw), REAL(epstot(3,1, iw)), aimag(epstot(3,1, iw))
         WRITE(48,"(4f15.6)") wgrid(iw), REAL(epstot(3,2, iw)), aimag(epstot(3,2, iw))
         WRITE(49,"(4f15.6)") wgrid(iw), REAL(epstot(3,3, iw)), aimag(epstot(3,3, iw))
         !
      ENDDO
      !
      CLOSE(30)
      CLOSE(40)
      CLOSE(41)
      CLOSE(42)
      !
  ENDIF

  !
  ! local cleaning
  !
  CALL grid_destroy()
  DEALLOCATE ( dipoletot, dipole_aux, epstot )

END SUBROUTINE offdiag_calc


!--------------------------------------------------------------------
SUBROUTINE dipole_calc( ik, dipole_aux, metalcalc, nbndmin, nbndmax )
  !------------------------------------------------------------------
  USE kinds,                ONLY : DP
  USE wvfct,                ONLY : npw, nbnd, igk, g2kin, ecutwfc
  USE wavefunctions_module, ONLY : evc
  USE klist,                ONLY : xk
  USE cell_base,            ONLY : tpiba2
  USE gvect,                ONLY : ngm, g
  USE io_files,             ONLY : nwordwfc, iunwfc
  USE grid_module,          ONLY : focc
  USE mp_global,            ONLY : intra_bgrp_comm
  USE mp,                   ONLY : mp_sum

IMPLICIT NONE
  !
  ! global variables
  INTEGER, INTENT(in)        :: ik,nbndmin,nbndmax
  COMPLEX(DP), INTENT(inout) :: dipole_aux(3,nbnd,nbnd)
  LOGICAL, INTENT(in)        :: metalcalc
  !
  ! local variables
  INTEGER :: iband1,iband2,ig
  COMPLEX(DP)   :: caux

  !
  ! Routine Body
  !
  CALL start_clock( 'dipole_calc' )

  !
  ! setup k+G grids for each kpt
  !
  CALL gk_sort (xk (1, ik), ngm, g, ecutwfc / tpiba2, npw, igk, g2kin)
  !
  ! read wfc for the given kpt
  !
  CALL davcio (evc, 2*nwordwfc, iunwfc, ik, - 1)
  !
  ! compute matrix elements
  !
  dipole_aux(:,:,:) = (0.0_DP,0.0_DP)
  !
  DO iband2 = nbndmin,nbndmax
      IF ( focc(iband2,ik) <  2.0d0) THEN
  DO iband1 = nbndmin,nbndmax
      !
      IF ( iband1==iband2 ) CYCLE
      IF ( focc(iband1,ik) >= 1e-4 ) THEN
            !
            DO  ig=1,npw
                 !
                 caux= conjg(evc(ig,iband1))*evc(ig,iband2)
                 !
                 dipole_aux(:,iband1,iband2) = dipole_aux(:,iband1,iband2) + &
                       ( g(:,igk(ig)) ) * caux
                 !
            ENDDO
      ENDIF
      !
  ENDDO
      ENDIF
  ENDDO
  !
  ! The diagonal terms are taken into account only if the system is treated like a metal, not
  ! in the intraband therm. Because of this we can recalculate the diagonal component of the dipole
  ! tensor directly as we need it for the intraband therm, without interference with interband one.
  !
  IF (metalcalc) THEN
     !
     DO iband1 = nbndmin,nbndmax
        DO  ig=1,npw
          !
          caux= conjg(evc(ig,iband1))*evc(ig,iband1)
          !
          dipole_aux(:,iband1,iband1) = dipole_aux(:,iband1,iband1) + &
                                        ( g(:,igk(ig))+ xk(:,ik) ) * caux
          !
        ENDDO
     ENDDO
     !
  ENDIF
  !
  ! recover over G parallelization (intra_bgrp)
  !
  CALL mp_sum( dipole_aux, intra_bgrp_comm )
  !
  CALL stop_clock( 'dipole_calc' )
  !
END SUBROUTINE dipole_calc


!-------------------------------------------------
SUBROUTINE occ_calc ()
  !-------------------------------------------------
  !
  USE kinds,     ONLY : DP
  USE klist,     ONLY : nkstot, wk, degauss
  USE wvfct,     ONLY : nbnd, wg, et
  USE ener,      ONLY : ef
  USE mp_global, ONLY : me_pool
  !
  IMPLICIT NONE
  !
  REAL(DP), ALLOCATABLE  :: focc(:,:),foccp(:,:)
  CHARACTER(25)          :: filename
  INTEGER                :: ierr, i, ik
  !
  ALLOCATE ( focc( nbnd, nkstot), STAT=ierr )
  IF (ierr/=0) CALL errore('grid_build','allocating focc', abs(ierr))
  !
  ALLOCATE ( foccp( nbnd, nkstot), STAT=ierr )
  IF (ierr/=0) CALL errore('grid_build','allocating foccp', abs(ierr))

  IF (me_pool==0) THEN
      !
      filename = 'occupations.dat'
!      WRITE(filename,"(I3,'.occupation.dat')")me_pool
      OPEN (unit=50, file=trim(filename))
      WRITE(50,*) '#energy (Ry)      occupation factor       derivative'

      DO ik = 1,nkstot
      DO i  = 1,nbnd
           focc(i,ik)= wg(i, ik ) * 2.0_DP/wk( ik )
           foccp(i,ik)= 2* exp((et(i,ik)-ef)/degauss)/((1+exp((et(i,ik)-ef)/degauss))**2*degauss)
           WRITE(50,*)et(i,ik),focc(i,ik),foccp(i,ik)
      ENDDO
      ENDDO

      CLOSE (50)
      !
  ENDIF
  !
  DEALLOCATE ( focc, STAT=ierr)
  CALL errore('grid_destroy','deallocating grid stuff',abs(ierr))
  !
  DEALLOCATE ( foccp, STAT=ierr)
  CALL errore('grid_destroy','deallocating grid stuff',abs(ierr))

END SUBROUTINE occ_calc
