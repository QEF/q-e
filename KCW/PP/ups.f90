!
! Copyright (C) 2004-2009 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!------------------------------
 MODULE grid_module
!------------------------------
  USE kinds,        ONLY : DP
  IMPLICIT NONE
  PRIVATE

  !
  ! general purpose vars
  !
  INTEGER                :: nw
  REAL(DP)               :: wmax, wmin
  REAL(DP)               :: alpha, full_occ
  REAL(DP), ALLOCATABLE  :: focc(:,:), wgrid(:)
  !
  PUBLIC :: grid_build, grid_destroy
  PUBLIC :: nw, wmax, wmin
  PUBLIC :: focc, wgrid, alpha, full_occ
  !
CONTAINS

!---------------------------------------------
  SUBROUTINE grid_build(nw_, wmax_, wmin_, metalcalc)
  !-------------------------------------------
  !
  USE kinds,     ONLY : DP
  USE io_global, ONLY : stdout, ionode
  USE wvfct,     ONLY : nbnd, wg
  USE klist,     ONLY : nks, wk, nelec
  USE lsda_mod,  ONLY : nspin
  USE uspp,      ONLY : okvan
  !
  IMPLICIT NONE
  !
  ! input vars
  INTEGER,  INTENT(IN) :: nw_
  REAL(DP), INTENT(IN) :: wmax_ ,wmin_
  LOGICAL,  OPTIONAL, INTENT(IN) :: metalcalc
  !
  ! local vars
  INTEGER         :: iw,ik,i,ierr

  !
  ! check on the number of bands: we need to include empty bands in order
  ! to compute the transitions
  !
  IF ( nspin == 1) full_occ = 2.0d0
  IF ( nspin == 2 .OR. nspin == 4) full_occ = 1.0d0
  !
  IF ( nspin == 2 ) THEN
     IF ( nbnd*full_occ <= nelec/2.d0 ) CALL errore('ups', 'bad band number', 2)
  ELSE
     IF ( nbnd*full_occ <= nelec ) CALL errore('ups', 'bad band number', 1)
  ENDIF
  !
  ! USPP are not implemented (dipole matrix elements are not trivial at all)
  !
  IF ( okvan ) CALL errore('grid_build','USPP are not implemented',1)

  !
  ! store data in module
  !
  nw = nw_
  wmax = wmax_
  wmin = wmin_

  !
  ! workspace
  !
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
        CALL errore('grid_build','non uniform kpt grid', ik )
     !
  ENDDO
  !
  ! occupation numbers, to be normalized differently
  ! whether we are spin resolved or not
  !
  DO ik = 1, nks
    DO i = 1, nbnd
        focc(i, ik) = wg(i, ik) * full_occ / wk(ik)
    ENDDO
  ENDDO

  !
  ! set the energy grid
  !
  IF ( metalcalc .AND. ABS(wmin) <= 0.001d0 ) wmin=0.001d0
  IF ( ionode ) WRITE(stdout,"(5x,a,f12.6)") "metallic system: redefining wmin = ", wmin  
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
  IF ( ALLOCATED( focc) ) THEN
      !
      DEALLOCATE ( focc, wgrid, STAT=ierr)
      CALL errore('grid_destroy','deallocating grid stuff',abs(ierr))
      !
  ENDIF
  !
END SUBROUTINE grid_destroy

END MODULE grid_module
!
!------------------------------
PROGRAM ups
!------------------------------
  !
  ! Compute the Ultraviolet photoemission spectrum.
  ! See Phys. Rev. Lett. 114, 166405 for details
  ! Adapted from epsilon.f90
  !
  ! Authors: 
  !     2014    Linh Nguyen, Andrea Ferrett:   basic implementation (using routine from epsilon.f90)
  !     2022    Nicola Colonna:                Ported to the latest QE and restructured
  !
  USE kinds,       ONLY : DP
  USE io_global,   ONLY : stdout, ionode, ionode_id
  USE mp,          ONLY : mp_bcast
  USE mp_global,   ONLY : mp_startup, mp_global_end
  USE mp_images,   ONLY : intra_image_comm
  USE io_files,    ONLY : tmp_dir, prefix
  USE constants,   ONLY : RYTOEV
  USE ener,        ONLY : ef
  USE klist,       ONLY : lgauss, ltetra
  USE wvfct,       ONLY : nbnd
  USE lsda_mod,    ONLY : nspin
  USE environment, ONLY : environment_start, environment_end
  USE grid_module, ONLY : grid_build, grid_destroy
  !
  IMPLICIT NONE
  !
  CHARACTER(LEN=256), EXTERNAL :: trimcheck
  CHARACTER(LEN=256) :: outdir
  !
  ! input variables
  !
  INTEGER                 :: nw,nbndmin,nbndmax
  REAL(DP)                :: intersmear,intrasmear,wmax,wmin,shift, & 
                             ekin_eout, ekin_error , photon_ener, & 
                             polar_angle, azimuthal_angle, & 
                             photon_angle, e_fermi
  CHARACTER(10)           :: calculation,smeartype
  LOGICAL                 :: metalcalc, homo_gas, wfc_real, & 
                             modified_pw, othor_pw
  !
  NAMELIST / inputpp / prefix, outdir, calculation
  NAMELIST / energy_grid / smeartype, intersmear, intrasmear, nw, wmax, wmin, &
                           nbndmin,nbndmax,shift,ekin_eout,ekin_error, &
                           polar_angle, azimuthal_angle, photon_angle, &
                           photon_ener, homo_gas, wfc_real, e_fermi,   &
                           modified_pw, othor_pw                          
  
  !
  ! local variables
  !
  INTEGER :: ios
  LOGICAL :: needwf = .TRUE.

!---------------------------------------------
! program body
!---------------------------------------------
!
  ! initialise environment
  !
#if defined(__MPI)
  CALL mp_startup ( )
#endif
  CALL environment_start ( 'ups' )
  !
  ! Set default values for variables in namelist
  !
  calculation  = 'eps'
  prefix       = 'pwscf'
  shift        = 0.0d0
  CALL get_environment_variable( 'ESPRESSO_TMPDIR', outdir )
  IF ( trim( outdir ) == ' ' ) outdir = './'
  intersmear   = 0.136
  wmin         = 0.0d0
  wmax         = 30.0d0
  nbndmin      = 1
  nbndmax      = 0
  nw           = 600
  smeartype    = 'gauss'
  intrasmear   = 0.0d0
  metalcalc    = .FALSE.
  !
  ! PHOTOSPEC additional variables
  ekin_eout    = 30.0d0 
  ekin_error   = 1.0d0 
  polar_angle  = 30
  azimuthal_angle  = 30
  photon_angle  = 50
  photon_ener  = 80
  e_fermi      = 4.0
  wfc_real     = .true.
  homo_gas     = .true.
  modified_pw  = .true.
  othor_pw     = .false.
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
  CALL mp_bcast ( ios, ionode_id, intra_image_comm )
  IF (ios/=0) CALL errore('ups', 'reading namelist INPUTPP', abs(ios))
  !
  IF ( ionode ) THEN
     !
     READ (5, energy_grid, IOSTAT=ios)
     !
     tmp_dir = trimcheck(outdir)
     !
  ENDIF
  !
  CALL mp_bcast ( ios, ionode_id, intra_image_comm )
  IF (ios/=0) CALL errore('ups', 'reading namelist ENERGY_GRID', abs(ios))
  !
  ! ... Broadcast variables
  !
  IF (ionode) WRITE( stdout, "( 5x, 'Broadcasting variables...' ) " )

  CALL mp_bcast( smeartype, ionode_id, intra_image_comm )
  CALL mp_bcast( calculation, ionode_id, intra_image_comm )
  CALL mp_bcast( prefix, ionode_id, intra_image_comm )
  CALL mp_bcast( tmp_dir, ionode_id, intra_image_comm )
  CALL mp_bcast( shift, ionode_id, intra_image_comm )
  CALL mp_bcast( intrasmear, ionode_id, intra_image_comm )
  CALL mp_bcast( intersmear, ionode_id, intra_image_comm)
  CALL mp_bcast( wmax, ionode_id, intra_image_comm )
  CALL mp_bcast( wmin, ionode_id, intra_image_comm )
  CALL mp_bcast( nw, ionode_id, intra_image_comm )
  CALL mp_bcast( nbndmin, ionode_id, intra_image_comm )
  CALL mp_bcast( nbndmax, ionode_id, intra_image_comm )
  !
  CALL mp_bcast( ekin_eout,       ionode_id, intra_image_comm ) 
  CALL mp_bcast( ekin_error,      ionode_id, intra_image_comm) 
  CALL mp_bcast( polar_angle,     ionode_id, intra_image_comm) 
  CALL mp_bcast( azimuthal_angle, ionode_id, intra_image_comm ) 
  CALL mp_bcast( photon_angle,    ionode_id, intra_image_comm ) 
  CALL mp_bcast( photon_ener,     ionode_id, intra_image_comm ) 
  CALL mp_bcast( homo_gas,        ionode_id, intra_image_comm ) 
  CALL mp_bcast( wfc_real,        ionode_id, intra_image_comm ) 
  CALL mp_bcast( e_fermi,         ionode_id, intra_image_comm ) 
  CALL mp_bcast( modified_pw,     ionode_id, intra_image_comm ) 
  CALL mp_bcast( othor_pw,        ionode_id, intra_image_comm ) 
  !
  ! read PW simulation parameters from prefix.save/data-file.xml
  !
  IF (ionode) WRITE( stdout, "( 5x, 'Reading PW restart file...' ) " )

  CALL read_file_new( needwf )
  !
  ! few conversions
  !

  IF (ionode) WRITE(stdout,"(2/, 5x, 'Fermi energy [eV] is: ',f8.5)") ef *RYTOEV

  IF (lgauss .or. ltetra) THEN
      metalcalc=.TRUE.
      IF (ionode) WRITE( stdout, "( 5x, 'The system is a metal (occupations are not fixed)...' ) " )
  ELSE
      IF (ionode) WRITE( stdout, "( 5x, 'The system is a dielectric...' ) " )
  ENDIF

  IF (nbndmax == 0) nbndmax = nbnd

  !
  ! perform some consistency checks, 
  ! setup w-grid and occupation numbers
  !
  CALL grid_build(nw, wmax, wmin, metalcalc)
  !
  ! ... run the specific pp calculation
  !
  IF (ionode) WRITE(stdout,"(/, 5x, 'Performing ',a,' calculation...')") trim(calculation)
  CALL start_clock(trim(calculation))
  SELECT CASE ( trim(calculation) )
  !
  CASE ( 'photospec' )
      !
      IF (nspin > 2) CALL errore ('ups', 'photospec NOT implemented for non-collinear spin', 1)
      CALL photoemission_spectr_pw ( intersmear,intrasmear,nbndmin,nbndmax, &
                                     shift,nspin,metalcalc,polar_angle, azimuthal_angle,&
                                     photon_angle,photon_ener,homo_gas,wfc_real,e_fermi, &
                                     modified_pw, othor_pw)
  CASE DEFAULT
      !
      CALL errore('ups','invalid CALCULATION = '//trim(calculation),1)
      !
  END SELECT
  !
  CALL stop_clock(trim(calculation))
  IF ( ionode ) WRITE( stdout , "(/)" )
  !
  CALL print_clock( trim(calculation) )
  CALL print_clock( 'dipole_calc' )
  IF ( ionode ) WRITE( stdout, *  )
  !
  ! cleaning
  !
  CALL grid_destroy()
  !
  CALL environment_end ( 'ups' )
  !
#if defined(__MPI)
  CALL mp_global_end()
#endif

END PROGRAM ups

!----------------------------------------------------------------------------------------
SUBROUTINE photoemission_spectr_pw ( intersmear,intrasmear,nbndmin,nbndmax, &
                                     shift,nspin,metalcalc,polar_angle, azimuthal_angle,&
                                     photon_angle,photon_ener,homo_gas,wfc_real,e_fermi,&
                                     modified_pw,othor_pw)
  !--------------------------------------------------------------------------------------
  !
  USE kinds,                ONLY : DP
  USE constants,            ONLY : PI, RYTOEV
  USE wvfct,                ONLY : npwx
  USE cell_base,            ONLY : tpiba2
  USE wvfct,                ONLY : et
  USE gvect,                ONLY : g
  USE klist,                ONLY : nks, xk, igk_k, ngk
  USE io_global,            ONLY : ionode, stdout
  USE grid_module,          ONLY : nw, wgrid, grid_build, grid_destroy  
  USE mp_global,            ONLY : intra_pool_comm
  USE mp,                   ONLY : mp_sum
  !
  USE lsda_mod,             ONLY : lsda, current_spin, isk
  ! 
  IMPLICIT NONE

  !
  ! input variables
  !
  INTEGER,         INTENT(IN) :: nbndmin,nbndmax,nspin
  REAL(DP),        INTENT(IN) :: intersmear, shift, photon_ener,  & 
                                 polar_angle, azimuthal_angle, photon_angle,  &
                                 e_fermi, intrasmear
  LOGICAL,         INTENT(IN) :: metalcalc, homo_gas, wfc_real, modified_pw,  &
                                 othor_pw
  !
  ! local variables
  !
  INTEGER       :: ik, iband1, ipol, ig
  INTEGER       :: iw, ierr
  INTEGER       :: npw
  REAL(DP)      :: etrans, w
  REAL(DP)      :: polar_angle_radial, azimuthal_angle_radial, photon_angle_radial
  REAL(DP)      :: ekin_eout
  REAL(DP)      :: module_k, kx, ky, kz, delta_ecut_G, &
                   delta_kx_Gx, delta_ky_Gy, delta_kz_Gz, constant,  & 
                   max_sigma_tot

  REAL(DP)      :: ssigma_tot (2, nbndmax)
  REAL(DP)      :: ssigma_2   (2,nbndmax)
  REAL(DP)      :: sbeta      (2,nbndmax)
  REAL(DP)      :: seigen_mol (2,nbndmax)
  !
  REAL(DP), ALLOCATABLE    :: srphotospec(:,:)
  REAL(DP), ALLOCATABLE    :: g2kin (:)
  REAL(DP), ALLOCATABLE    :: dipole_2(:,:)
  REAL(DP), ALLOCATABLE    :: gamma_2_opw(:,:)
  REAL(DP), ALLOCATABLE    :: lambda_2_opw(:,:)
  COMPLEX(DP),ALLOCATABLE  :: dipole_aux(:,:,:)
  COMPLEX(DP),ALLOCATABLE  :: dipole_opw_gamma(:,:,:)
  COMPLEX(DP),ALLOCATABLE  :: scala_opw_lambda(:,:)
  !
  INTEGER :: iss
  !
  !--------------------------
  ! main routine body
  !--------------------------
  !
  ! change angles from degree to radial unit
  !
  polar_angle_radial = (polar_angle/180.0)*PI
  azimuthal_angle_radial = (azimuthal_angle/180.0)*PI
  photon_angle_radial = (photon_angle/180.0)*PI
  !
  ! allocate main spectral and auxiliary quantities
  !
  ALLOCATE( g2kin(npwx), STAT=ierr )
  IF (ierr/=0) CALL errore('ups','allocating photoemission_spectr',ABS(ierr))
  !
  ALLOCATE( dipole_2(nbndmax, npwx), STAT=ierr ) 
  IF (ierr/=0) CALL errore('ups','allocating dipole', ABS(ierr) )
  !
  ALLOCATE( dipole_aux(3, nbndmax, npwx), STAT=ierr )
  IF (ierr/=0) CALL errore('ups','allocating dipole_aux', ABS(ierr) )
  !
  IF (othor_pw) THEN
    ALLOCATE( dipole_opw_gamma(3, nbndmax, npwx) )
    ALLOCATE( scala_opw_lambda(nbndmax, npwx) )
    ALLOCATE( gamma_2_opw(nbndmax, npwx) )
    ALLOCATE( lambda_2_opw(nbndmax, npwx) )
  ENDIF
  !
  ! allocate main spectral and auxiliary quantities
  !
  ALLOCATE( srphotospec(nspin,nw), STAT=ierr )
      IF (ierr/=0) CALL errore('ups','allocating photospectra',ABS(ierr))
  !
  ! initialize 
  !
  srphotospec(:,:)=0.0_DP
  !
  constant = 1.0_DP
  !
  IF (homo_gas) THEN
    ssigma_tot(:, :) = 0.0_DP
    ssigma_2(:, :)   = 0.0_DP
  ENDIF
  !
  ! main kpt loop
  !
  !
  DO ik = 1, nks
    !
    current_spin = 1
    IF (lsda) current_spin = isk(ik)
    !
    ! For every single k-point: order k+G for
    !                           read and distribute wavefunctions
    !                           compute dipole matrix 3 x nbnd x nbnd
    !                           parallel over g 
    !                           recover g parallelism getting the total
    !                           dipole matrix
    !
    npw = ngk(ik)
    !
    CALL dipole_calc_pw( ik, dipole_aux, dipole_opw_gamma, scala_opw_lambda, metalcalc, nbndmin, nbndmax, &
                         photon_angle_radial, azimuthal_angle_radial, &
                         homo_gas, othor_pw)
    !
    dipole_2(:,:) = 0.0_DP
    DO ipol = 1, 3
        dipole_2(:,:) = dipole_2(:,:) + tpiba2 * ( (REAL(dipole_aux(ipol,:,:)))**2 + (AIMAG(dipole_aux(ipol,:,:)))**2 )
    ENDDO
    !
    IF (othor_pw) THEN
      ! 
      gamma_2_opw(:,:) = 0.0_DP
      lambda_2_opw(:,:) = 0.0_DP
      DO ipol = 1, 3
        gamma_2_opw(:,:) = gamma_2_opw(:,:) + tpiba2 * & 
                   ( (REAL(dipole_opw_gamma(ipol,:,:)))**2 + (AIMAG(dipole_opw_gamma(ipol,:,:)))**2 )
      ENDDO
      !
      lambda_2_opw(:,:) = lambda_2_opw(:,:)  + tpiba2 * & 
                   scala_opw_lambda(:,:) * conjg(scala_opw_lambda(:,:))
      !
    ENDIF
    !
    ! Calculation of photoemission spectra
    ! 'intersmear' is the brodening parameter
    !
    DO iband1 = nbndmin, nbndmax 
      !
      IF (homo_gas) THEN
        !
        seigen_mol(current_spin, iband1) = (0.0_DP - et(iband1, ik) )*RYTOEV
        ! 
        IF (modified_pw ) THEN
           ekin_eout = photon_ener
        ELSE
           ekin_eout = photon_ener - (0.0_DP - et(iband1, ik) )*RYTOEV
        ENDIF
        !  
      ELSE
        ! 
        ekin_eout = photon_ener -((0.0_DP - et(iband1, ik) )*RYTOEV + e_fermi)
        !
      ENDIF
      ! 
      IF (ekin_eout > photon_ener .or. ekin_eout <=0.0_DP ) CYCLE
      !  
      module_k = sqrt ((ekin_eout/13.6056923)/tpiba2)
      ! 
      kx = module_k * cos(azimuthal_angle_radial) * sin(polar_angle_radial)
      ky = module_k * sin(azimuthal_angle_radial) * sin(polar_angle_radial)
      kz = module_k * cos(polar_angle_radial)
      !
      ! compute the total cross-section
      ! and ansymmetry parameter
      !
      IF (homo_gas) THEN
        !
        DO ig = 1, npw
          !
          g2kin(ig) = tpiba2*((xk(1,ik)+g(1,igk_k(ig,ik)))**2 + (xk(2,ik)+g(2,igk_k(ig,ik))) **2 + (xk(3,ik)+g(3,igk_k(ig,ik)))**2)
          delta_ecut_G =  intrasmear/( PI * (( g2kin(ig)*13.6056923 - ekin_eout )**2 + (intrasmear)**2 ))
          !
          ssigma_tot(current_spin, iband1) = ssigma_tot(current_spin, iband1) + constant*module_k*dipole_2(iband1, ig)*delta_ecut_G
          !
          IF (othor_pw ) THEN
             ssigma_2(current_spin, iband1) = ssigma_2(current_spin, iband1)   + 1.5*constant*module_k &
                         * (gamma_2_opw(iband1,ig) - lambda_2_opw(iband1,ig))*delta_ecut_G 
          ENDIF
          !
        ENDDO
        !
      ENDIF
      ! 
      DO ig = 1, npw
          g2kin(ig) = tpiba2*((xk(1,ik)+g(1,igk_k(ig,ik)))**2 + (xk(2,ik)+g(2,igk_k(ig,ik))) **2 + (xk(3,ik)+g(3,igk_k(ig,ik)))**2)
          delta_ecut_G =  intrasmear/( PI * (( g2kin(ig)*13.6056923 - ekin_eout  )**2 + (intrasmear)**2 ))
          ! 
          delta_kx_Gx  =  intrasmear/( PI * (( ( kx - ( xk(1,ik) + g(1,igk_k(ig,ik)))) *sqrt(tpiba2))**2 + (intrasmear)**2 ) )
          delta_ky_Gy  =  intrasmear/( PI * (( ( ky - ( xk(2,ik) + g(2,igk_k(ig,ik)))) *sqrt(tpiba2))**2 + (intrasmear)**2 ) )
          delta_kz_Gz  =  intrasmear/( PI * (( ( kz - ( xk(3,ik) + g(3,igk_k(ig,ik)))) *sqrt(tpiba2))**2 + (intrasmear)**2 ) )
          !
          ! transition energy 
          !
          etrans = (0.0d0 - et(iband1, ik) ) * RYTOEV  + shift ! g2kin were called in the dipole_pw routine
          !
          IF( etrans < 1.0d-10 ) CYCLE
          !
          ! loop over frequencies
          !
          DO iw = 1, nw
            !
            w = wgrid(iw)
            !
            IF (homo_gas) THEN
            !
              IF (wfc_real) THEN  
                !
                srphotospec(current_spin, iw) = srphotospec(current_spin, iw) + 2.0D0 * module_k* intersmear &
                                                * dipole_2(iband1, ig) * delta_ecut_G &
                                                / ( PI * ( (etrans - w )**2 + (intersmear)**2 )  ) 
                !
              ELSE
                !
                srphotospec(current_spin, iw) = srphotospec(current_spin, iw) + module_k * intersmear &
                                                * dipole_2(iband1, ig) * delta_ecut_G &
                                                / ( PI * ( (etrans - w )**2 + (intersmear)**2 )  ) 
                !
              ENDIF
              !
            ELSE
              !
              srphotospec(current_spin, iw) = srphotospec(current_spin, iw) + dipole_2(iband1, ig) * delta_ecut_G &
                   * delta_kx_Gx * delta_ky_Gy * delta_kz_Gz * intersmear & 
                   / ( PI * ( (etrans - w )**2 + (intersmear)**2 )  )
              !
            ENDIF
            !
          ENDDO
          ! 
      ENDDO
      !  
    ENDDO
    !
  ENDDO ! kpt_loop
  !
  ! recover over G parallelization (intra_pool)
  !
  DO iss = 1, nspin
    !
    CALL mp_sum( srphotospec(iss,:), intra_pool_comm )
    !
    IF(homo_gas) THEN
      !  
      CALL mp_sum( ssigma_tot(iss,:), intra_pool_comm )
      CALL mp_sum( ssigma_2(iss,:), intra_pool_comm )
      !
      max_sigma_tot = maxval(ssigma_tot(iss, :))
      ! 
      IF (othor_pw) THEN
        !
        DO iband1 = nbndmin, nbndmax
           sbeta(iss, iband1)= 2.0_DP*(1.0_DP-ssigma_2(iss, iband1)/ssigma_tot(iss, iband1) )
        ENDDO
        ! 
      ELSE
        ! 
        sbeta(iss,:) = 2.0_DP
        !
      ENDIF
      !
      IF (ionode) THEN
        IF (nspin == 1) THEN 
          WRITE(stdout,"(/,5x, 'Writing the molecule gas properties' )")
        ELSE
          IF (iss == 0) WRITE(stdout,"(/,5x, 'Writing the molecule gas properties with spin up' )") 
          IF (iss == 1) WRITE(stdout,"(/,5x, 'Writing the molecule gas properties with spin down' )") 
        ENDIF
        DO iband1 = nbndmin, nbndmax
          WRITE(stdout,"(4f15.6)") seigen_mol(iss, iband1), ssigma_tot(iss, iband1)/max_sigma_tot, sbeta(iss, iband1)
        ENDDO
      ENDIF
      !
    ENDIF
    !
  ENDDO
  !
  ! write results on data files
  !
  IF (ionode) THEN
     WRITE(stdout,"(/,5x, 'Writing output on file...' )")

     OPEN (30, FILE='photospectra.dat', FORM='FORMATTED' )
     !
     IF( nspin == 1) WRITE(30, "(2x,'# energy grid [eV]   photospectra ')" )
     IF( nspin == 2) WRITE(30, "(2x,'# energy grid [eV]   photospectra spin_up [ab.]   photospectra spin_dw [ab.] ')" )
     !
     DO iw =1, nw
         ! 
         IF (nspin == 1 ) WRITE(30,"(4f15.10)") wgrid(iw), srphotospec(1, iw)
         IF (nspin == 2 ) WRITE(30,"(4f15.10)") wgrid(iw), srphotospec(1, iw), srphotospec(2,iw) 
         !
     ENDDO
     !
     CLOSE(30)
     !
  ENDIF
  !
  ! local cleaning
  !
  DEALLOCATE (g2kin,STAT=ierr)
  DEALLOCATE (srphotospec,STAT=ierr)
  DEALLOCATE (dipole_2, STAT=ierr)
  DEALLOCATE (dipole_aux, STAT=ierr)
  ! 
  IF (othor_pw) THEN
    DEALLOCATE( dipole_opw_gamma )
    DEALLOCATE( scala_opw_lambda )
    DEALLOCATE( gamma_2_opw )
    DEALLOCATE( lambda_2_opw )
  ENDIF
  !
  CALL grid_destroy()
  !
  return
END SUBROUTINE photoemission_spectr_pw

!--------------------------------------------------------------------
SUBROUTINE dipole_calc_pw ( ik, dipole_aux, dipole_opw_gamma, scala_opw_lambda, metalcalc, &
                            nbndmin, nbndmax, photon_angle, azimuthal_angle, homo_gas, othor_pw) 
  !------------------------------------------------------------------
  USE kinds,                ONLY : DP
  USE wvfct,                ONLY : npwx
  USE wavefunctions,        ONLY : evc
  USE klist,                ONLY : xk, igk_k, ngk 
  USE gvect,                ONLY : g
  USE io_files,             ONLY : restart_dir
  USE mp_global,            ONLY : intra_pool_comm
  USE mp,                   ONLY : mp_sum
  USE pw_restart_new,       ONLY : read_collected_wfc

IMPLICIT NONE
  !
  ! global variables
  INTEGER, INTENT(IN)        :: ik,nbndmin, nbndmax
  REAL(DP),INTENT(IN)        :: photon_angle, azimuthal_angle 
  COMPLEX(DP), INTENT(INOUT) :: dipole_aux(3, nbndmax, npwx)
  COMPLEX(DP), INTENT(INOUT) :: dipole_opw_gamma(3, nbndmax, npwx)
  COMPLEX(DP), INTENT(INOUT) :: scala_opw_lambda(nbndmax, npwx)
  LOGICAL, INTENT(IN)        :: metalcalc, homo_gas, othor_pw
  !
  ! local variables
  INTEGER :: npw
  INTEGER :: iband1, iband2, ig
  REAL(DP):: sqrtk2
  COMPLEX(DP) :: caux1, caux2
  COMPLEX(DP) :: sumx, sumy, sumz  
  COMPLEX(DP) :: ax(npwx), ay(npwx), az(npwx) 
  COMPLEX(DP) :: ax_add(npwx), ay_add(npwx), az_add(npwx) 
  !
  ! Routine Body
  ! 
  CALL start_clock( 'dipole_calc' )
  !
  ! setup k+G grids for each kpt
  !
  !CALL gk_sort (xk (1, ik), ngm, g, ecutwfc / tpiba2, npw, igk, g2kin)  
  !
  ! read wfc for the given kpt
  !
  !CALL davcio (evc, nwordwfc, iunwfc, ik, - 1) 
  CALL read_collected_wfc ( restart_dir(), ik, evc )
  !
  ! compute matrix elements
  !
  dipole_aux(:,:,:) = (0.0_DP,0.0_DP)
  npw = ngk(ik)
  !
  IF (othor_pw) THEN
     dipole_opw_gamma(:,:,:) = (0.0_DP,0.0_DP)
     scala_opw_lambda(:,:) = (0.0_DP,0.0_DP)
  ENDIF
  !
  DO iband1 = nbndmin, nbndmax 
     !
     DO ig = 1, npw
        !
        caux1 = evc(ig,iband1)
        !
        IF (homo_gas) then 
           dipole_aux(1,iband1,ig) = dipole_aux(1,iband1,ig) + &
                                     ( g(1,igk_k(ig,ik)) ) * caux1
           dipole_aux(2,iband1,ig) = dipole_aux(2,iband1,ig) + &
                                     ( g(2,igk_k(ig,ik)) ) * caux1
           dipole_aux(3,iband1,ig) = dipole_aux(3,iband1,ig) + &
                                     ( g(3,igk_k(ig,ik)) ) * caux1
        ELSE
           dipole_aux(1,iband1,ig) = dipole_aux(1,iband1,ig) + &
                                     ( g(1,igk_k(ig, ik)) )* cos(photon_angle) &
                                      * cos(azimuthal_angle) * caux1
           dipole_aux(2,iband1,ig) = dipole_aux(2,iband1,ig) + &
                                     ( g(2,igk_k(ig,ik)) )* cos(photon_angle) &
                                      * sin(azimuthal_angle) * caux1
           dipole_aux(3,iband1,ig) = dipole_aux(3,iband1,ig) + &
                                     ( g(3,igk_k(ig, ik)) )* sin(photon_angle) * caux1
        ENDIF
        ! 
     ENDDO
     !
     IF (othor_pw) THEN
        !
        ax(:) = (0.0_DP,0.0_DP) 
        ay(:) = (0.0_DP,0.0_DP) 
        az(:) = (0.0_DP,0.0_DP)
        ax_add(:) = (0.0_DP,0.0_DP)
        ay_add(:) = (0.0_DP,0.0_DP)
        az_add(:) = (0.0_DP,0.0_DP) 
        !
        DO iband2 = nbndmin, nbndmax
          !
          sumx = (0.0_DP,0.0_DP)
          sumy = (0.0_DP,0.0_DP)
          sumz = (0.0_DP,0.0_DP)
          ! 
          DO ig = 1, npw
            caux1 = evc(ig,iband1) 
            caux2 = evc(ig,iband2) 
            if (homo_gas) then
               sumx = sumx +  g(1,igk_k(ig,ik)) * conjg(caux1)*caux2    
               sumy = sumy +  g(2,igk_k(ig,ik)) * conjg(caux1)*caux2    
               sumz = sumz +  g(3,igk_k(ig,ik)) * conjg(caux1)*caux2    
            else
               sumx = sumx +  g(1,igk_k(ig,ik)) * conjg(caux1)*caux2 * &
                           cos(photon_angle) * cos(azimuthal_angle)
               sumy = sumy +  g(2,igk_k(ig,ik)) * conjg(caux1)*caux2 * &
                           cos(photon_angle) * sin(azimuthal_angle) 
               sumz = sumz +  g(3,igk_k(ig,ik)) * conjg(caux1)*caux2 * &
                           sin(photon_angle)
            endif 
          ENDDO
          !
          CALL mp_sum( sumx, intra_pool_comm )
          CALL mp_sum( sumy, intra_pool_comm )
          CALL mp_sum( sumz, intra_pool_comm )
          ! 
          DO ig = 1, npw
            caux2  = evc(ig,iband2)
            ax(ig) = ax(ig) + conjg(caux2)*sumx 
            ay(ig) = ay(ig) + conjg(caux2)*sumy 
            az(ig) = az(ig) + conjg(caux2)*sumz 
          ENDDO
          !
          DO ig = 1, npw
            caux2  = evc(ig,iband2)
            ax_add(ig) = ax_add(ig) + conjg(caux2)*sumx*g(1,igk_k(ig,ik))
            ay_add(ig) = ay_add(ig) + conjg(caux2)*sumy*g(2,igk_k(ig,ik))
            az_add(ig) = az_add(ig) + conjg(caux2)*sumz*g(3,igk_k(ig,ik))
          ENDDO
          !
        ENDDO
        !
        DO ig = 1, npw
          !
          ! gradient component of total sigma
          ! 
          dipole_aux(1,iband1,ig) = dipole_aux(1,iband1,ig) - ax(ig)
          dipole_aux(2,iband1,ig) = dipole_aux(2,iband1,ig) - ay(ig)
          dipole_aux(3,iband1,ig) = dipole_aux(3,iband1,ig) - az(ig)
          !
          ! gradient component of Gamma
          !   
          dipole_opw_gamma(1,iband1,ig) = dipole_opw_gamma(1,iband1,ig) + (ax(ig))
          dipole_opw_gamma(2,iband1,ig) = dipole_opw_gamma(2,iband1,ig) + (ay(ig))
          dipole_opw_gamma(3,iband1,ig) = dipole_opw_gamma(3,iband1,ig) + (az(ig))
          !   
          sqrtk2 = sqrt((xk(1,ik)+g(1,igk_k(ig,ik)))**2 + (xk(2,ik)+g(2,igk_k(ig,ik))) **2 + (xk(3,ik)+g(3,igk_k(ig,ik)))**2)
          !
          IF (sqrtk2 > 1.0E-05) THEN
            !
            ! scala component of Lambda
            !
            scala_opw_lambda(iband1,ig) = scala_opw_lambda(iband1,ig) + & 
                                          (ax_add(ig) + ay_add(ig) + az_add(ig))/sqrtk2      
            !    
          ENDIF
          !
        ENDDO
        !
     ENDIF
     ! 
  ENDDO 
  !
  CALL stop_clock( 'dipole_calc' )
  !
  return
  !
END SUBROUTINE dipole_calc_pw

