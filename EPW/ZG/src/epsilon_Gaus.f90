!
! Copyright (C) 2004-2009 Andrea Benassi and Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! This is the routine epsilon.f90. 

! A.D. and TYK (02/18/2026): We have added symmetrization of the matrix
! elements and the addition of the nonlocal pseudopotential contribution 
! in eps_calc/dipole_calc2. 

! mz: We make a small modification so that 
! Gaussian broadening is applied on the spectra and avoid numerical artefacts
! at the tails of indirect optical absorption. Changes can be traced 
! by "mz_b" and "mz_e". 
! 
! To compile this routine, it requires to have pp compiled. 
! i.e in main directory type "make pp". 
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
     IF ( nbnd*full_occ <= nelec/2.d0 ) CALL errore('epsilon', 'bad band number', 2)
  ELSE
     IF ( nbnd*full_occ <= nelec ) CALL errore('epsilon', 'bad band number', 1)
  ENDIF


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
  IF ( metalcalc .AND. ABS(wmin) <= 0.001d0 ) THEN
     wmin=0.001d0
     IF ( ionode ) WRITE(stdout,"(5x,a,f18.6)") "metallic system: redefining wmin = ", wmin  
  ENDIF
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
!---------------------------------------------------------------------------
END MODULE grid_module
!
MODULE eps_writer
!------------------------------
! A wrapper for the function that writes the dielectric function
  IMPLICIT NONE
  !
  PRIVATE
  !
  PUBLIC :: eps_writetofile
  !
CONTAINS
!--------------------------------------------------------------------
SUBROUTINE eps_writetofile(namein,desc,nw,wgrid,ncol,var,desc2)
  !------------------------------------------------------------------
  ! Writes the dielectric function to file.
  ! 
  USE kinds,          ONLY : DP
  USE io_files,       ONLY : prefix, tmp_dir
  !
  IMPLICIT NONE
  !
  CHARACTER(LEN=*),   INTENT(IN)           :: namein
  ! Prefix of the output file
  CHARACTER(LEN=*),   INTENT(IN)           :: desc
  ! Header line for the output file
  INTEGER,            INTENT(IN)           :: nw
  ! Number of frequencies 
  INTEGER,            INTENT(IN)           :: ncol
  ! Number of columns in the output
  REAL(DP),           INTENT(IN)           :: wgrid(nw)
  ! The values of the frequencies
  REAL(DP),           INTENT(IN)           :: var(ncol,nw)
  ! The dielectric function component of choice
  CHARACTER(LEN=*),   INTENT(IN), OPTIONAL :: desc2
  ! Optional second header line
  !
  CHARACTER(256) :: str
  ! Temporary string for the output file name
  INTEGER        :: iw
  ! Loop variable for the frequencies
  !
  str = TRIM(namein) // "_" // TRIM(prefix) // ".dat"
  OPEN(40,FILE=TRIM(str))
  !
  WRITE(40,"(a)") "# "// TRIM(desc)
  !
  IF (PRESENT(desc2)) THEN
    WRITE(40, "(a)") "# "// TRIM(desc2)
  ELSE
    WRITE(40,"(a)") "#"
  END IF
  !
  DO iw = 1, nw
     !
     WRITE(40,"(10f30.9)") wgrid(iw), var(1:ncol,iw)
     !
  ENDDO
  !
  CLOSE(40)
  !
END SUBROUTINE eps_writetofile
!
END MODULE eps_writer
!
!------------------------------
PROGRAM epsilon
!------------------------------
  !
  ! Compute the complex macroscopic dielectric function,
  ! at the RPA level, neglecting local field effects.
  ! Eps is computed both on the real or imaginary axis
  !
  ! Authors: 
  !     2006    Andrea Benassi, Andrea Ferretti, Carlo Cavazzoni:   basic implementation (partly taken from pw2gw.f90)
  !     2007    Andrea Benassi:                                     intraband contribution, nspin=2
  !     2016    Tae-Yun Kim, Cheol-Hwan Park:                       bugs fixed
  !     2016    Tae-Yun Kim, Cheol-Hwan Park, Andrea Ferretti:      non-collinear magnetism implemented
  !                                                                 code significantly restructured
  !     2026    Tae-Yun Kim, Adam Denchfield:                       symmetry reduction and USPP implemented for epsilon
  USE kinds,       ONLY : DP
  USE io_global,   ONLY : stdout, ionode, ionode_id
  USE mp,          ONLY : mp_bcast
  USE mp_global,   ONLY : mp_startup
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
  REAL(DP)                :: intersmear,intrasmear,wmax,wmin,shift
  CHARACTER(10)           :: calculation,smeartype
  LOGICAL                 :: metalcalc
  !
  NAMELIST / inputpp / prefix, outdir, calculation
  NAMELIST / energy_grid / smeartype, intersmear, intrasmear, nw, wmax, wmin, &
                           nbndmin, nbndmax, shift
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
  CALL environment_start ( 'epsilon' )
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
  CALL mp_bcast ( ios, ionode_id, intra_image_comm )
  IF (ios/=0) CALL errore('epsilon', 'reading namelist ENERGY_GRID', abs(ios))
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
  CALL mp_bcast( metalcalc, ionode_id, intra_image_comm )

  !
  ! read PW simulation parameters from prefix.save/data-file.xml
  !
  IF (ionode) WRITE( stdout, "( 5x, 'Reading PW restart file...' ) " )

  CALL read_file_new( needwf )
  !
  ! few conversions
  !

  IF (ionode) WRITE(stdout,"(2/, 5x, 'Fermi energy [eV] is: ',f10.5)") ef *RYTOEV

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
  CASE ( 'eps' )
      !
      CALL eps_calc ( intersmear, intrasmear, nbndmin, nbndmax, shift, metalcalc, nspin )
      !
  CASE ( 'jdos' )
      !
      CALL jdos_calc ( smeartype, intersmear, nbndmin, nbndmax, shift, nspin )
      !
  CASE ( 'offdiag' )
      !
      CALL offdiag_calc ( intersmear, intrasmear, nbndmin, nbndmax, shift, metalcalc, nspin )
      !
  CASE DEFAULT
      !
      CALL errore('epsilon','invalid CALCULATION = '//trim(calculation),1)
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
  CALL environment_end ( )
  !
  CALL stop_pp ()

END PROGRAM epsilon


!-----------------------------------------------------------------------------
SUBROUTINE eps_calc ( intersmear,intrasmear, nbndmin, nbndmax, shift, metalcalc , nspin)
  !-----------------------------------------------------------------------------
  ! 
  USE kinds,                ONLY : DP
  USE constants,            ONLY : PI, RYTOEV
  USE cell_base,            ONLY : tpiba2, omega
  USE wvfct,                ONLY : nbnd, et
  USE ener,                 ONLY : efermi => ef
  USE klist,                ONLY : nks, nkstot, degauss, ngauss
  USE io_global,            ONLY : ionode, stdout
  !
  USE grid_module,          ONLY : alpha, focc, full_occ, nw, wgrid, grid_destroy
  USE eps_writer,           ONLY : eps_writetofile
  USE mp_pools,             ONLY : inter_pool_comm
  USE mp,                   ONLY : mp_sum
  USE uspp,                 ONLY : nkb
  USE becmod,               ONLY : becp, allocate_bec_type, deallocate_bec_type
  !
  IMPLICIT NONE

  !
  ! input variables
  !
  INTEGER,         INTENT(in) :: nbndmin, nbndmax, nspin
  REAL(DP),        INTENT(in) :: intersmear, intrasmear, shift
  LOGICAL,         INTENT(in) :: metalcalc
  !
  ! local variables
  !
  INTEGER       :: i, ik, iband1, iband2,is
  INTEGER       :: iw, iwp, ierr
  REAL(DP)      :: etrans, const, w, renorm(3)
  CHARACTER(128):: desc(4)
  CHARACTER(128):: desc2
  !
  REAL(DP),    ALLOCATABLE :: epsr(:,:), epsi(:,:), epsrc(:,:,:), epsic(:,:,:)
  REAL(DP),    ALLOCATABLE :: ieps(:,:), eels(:,:), iepsc(:,:,:), eelsc(:,:,:)
  REAL(DP),    ALLOCATABLE :: dipole(:,:,:)
  COMPLEX(DP), ALLOCATABLE :: dipole_aux(:,:,:)
  !
  REAL(DP) , EXTERNAL :: w0gauss
!
!--------------------------
! main routine body
!--------------------------
!
    !
    ! allocate main spectral and auxiliary quantities
    !
    ALLOCATE( dipole(3, nbnd, nbnd), STAT=ierr )
    IF (ierr/=0) CALL errore('epsilon','allocating dipole', abs(ierr) )
    !
    ALLOCATE( dipole_aux(3, nbnd, nbnd), STAT=ierr )
    IF (ierr/=0) CALL errore('epsilon','allocating dipole_aux', abs(ierr) )
    !
    ALLOCATE( epsr( 3, nw), epsi( 3, nw), eels( 3, nw), ieps(3,nw ), STAT=ierr )
    IF (ierr/=0) CALL errore('epsilon','allocating eps', abs(ierr))
    ! allocate beta projectors for calculation of the non-local contributions
    CALL allocate_bec_type(nkb, nbnd, becp)
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
        ! NOTE: For debug purpose
        !
        ! Option 1: no nonlocal commutatator + symmetrization
        ! CALL dipole_calc( ik, dipole_aux, metalcalc, nbndmin, nbndmax)
        ! CALL dipole_symmetrize( ik, dipole_aux, dipole, nbndmin, nbndmax)
        ! dipole = dipole * tpiba2

        ! Option 2: nonlocal commutator + symmetrization
        !           * commutator_Hx_psi includes tpiba, so no need to multiply tpiba2 at the end
        CALL dipole_calc_2( ik, dipole_aux, metalcalc, nbndmin, nbndmax)


!        IF (nkstot > 1) THEN 
        CALL dipole_symmetrize( ik, dipole_aux, dipole, nbndmin, nbndmax)
!        ELSE ! no need to symmetrize
!           dipole(:,:,:)= tpiba2 * REAL( dipole_aux(:,:,:) * conjg(dipole_aux(:,:,:)), DP )
!        ENDIF

        ! Option 3: no nonlocal commutator + no symmetrization
        ! CALL dipole_calc( ik, dipole_aux, metalcalc, nbndmin, nbndmax)
        ! dipole(:,:,:)= tpiba2 * REAL( dipole_aux(:,:,:) * conjg(dipole_aux(:,:,:)), DP )

        ! Calculation of real and imaginary parts
        ! of the macroscopic dielettric function from dipole
        ! approximation.
        ! 'intersmear' is the brodening parameter
        !
        ! Interband
        !
        DO iband2 = nbndmin, nbndmax
          IF ( focc(iband2, ik) < full_occ ) THEN
              DO iband1 = nbndmin,nbndmax
                IF ( iband1 == iband2 ) CYCLE
                IF ( focc(iband1, ik) >= 0.5d-4*full_occ ) THEN
                    IF ( abs(focc(iband2,ik) - focc(iband1,ik)) < 1.0d-3*full_occ ) CYCLE
                    ! transition energy
                    etrans = ( et(iband2, ik) - et(iband1, ik) ) * RYTOEV + shift
                    ! loop over frequencies
                    DO iw = 1, nw
                      w = wgrid(iw)
                      ! mz_modifies to have Gaus broad for epsi2; A.D. modified to include +w and -w terms, and factor of pi/2 for normalization
                   epsi(:, iw) = epsi(:, iw) + dipole(:, iband1, iband2) / etrans**2   &
                        * RYTOEV**3 * (focc(iband1, ik))   &
                        / intersmear * sqrt(3.141592653589793) * 0.5*(EXP(-(etrans - w)**2 &
                        / intersmear**2)+EXP(-(etrans + w)**2 / intersmear**2))
                   ! mz_ends
                      epsr(:, iw) = epsr(:, iw) + dipole(:, iband1, iband2) * RYTOEV**3 * &
                            (focc(iband1, ik)) * &
                            (etrans**2 - w**2 ) / &
                            (( (etrans**2 -w**2 )**2 + intersmear**2 * w**2 )* etrans )
                    ENDDO
                ENDIF
              ENDDO
          ENDIF
        ENDDO

        !
        ! Intraband (only if metalcalc is true)
        !
        IF (metalcalc) THEN
          DO iband1 = nbndmin,nbndmax
              ! loop over frequencies
              DO iw = 1, nw
                w = wgrid(iw)
                epsi(:, iw) = epsi(:, iw) +  dipole(:, iband1, iband1) * intrasmear * w * &
                      RYTOEV**2 * w0gauss((et(iband1,ik)-efermi)/degauss, ngauss) / &
                      (( w**4 + intrasmear**2 * w**2 )*degauss ) * (0.5d0 * full_occ)
                epsr(:,iw) = epsr(:,iw) - dipole(:,iband1,iband1) * RYTOEV**2 * &
                      w0gauss((et(iband1,ik)-efermi)/degauss, ngauss) * w**2 / &
                      (( w**4 + intrasmear**2 * w**2 )*degauss ) * (0.5d0 * full_occ)
              ENDDO
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
    const = 64.0d0 * PI / ( omega ) * full_occ
    !
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
            ieps(:, iw) = ieps(:, iw) + wgrid(iwp) * epsi(:, iwp) / ( wgrid(iwp)**2 + wgrid(iw)**2)
            !
        ENDDO
    ENDDO
    !
    ieps(:,:) = 1.0d0 + 2.0d0/PI * ieps(:,:) * alpha

    !
    ! check  dielectric function  normalizzation via sumrule
    !
    DO i=1,3
        renorm(i) = alpha * SUM( epsi(i, :) * wgrid(:) )
    ENDDO
    renorm(:) = SQRT( renorm(:) * 2.0d0/PI) 
    !
    IF ( ionode ) THEN
        !
        WRITE(stdout,"(/,5x, 'xx,yy,zz plasmon frequences [eV] are: ',3f30.9 )")  renorm(:)
        WRITE(stdout,"(/,5x, 'Writing output on file...' )")

        !
        ! write results on data files
        !
        desc(1) = "energy grid [eV]     epsr_x  epsr_y  epsr_z"
        WRITE(desc2, "('plasmon frequences [eV]: ',3f30.9)") renorm (:)
        !
        desc(2) = "energy grid [eV]     epsi_x  epsi_y  epsi_z"
        desc(3) = "energy grid [eV]  eels components [arbitrary units]"
        desc(4) = "energy grid [eV]     ieps_x  ieps_y  ieps_z"
        !
        CALL eps_writetofile("epsr",desc(1),nw,wgrid,3,epsr,desc2)
        CALL eps_writetofile("epsi",desc(2),nw,wgrid,3,epsi)
        CALL eps_writetofile("eels",desc(3),nw,wgrid,3,eels)
        CALL eps_writetofile("ieps",desc(4),nw,wgrid,3,ieps)
        !
    ENDIF

    DEALLOCATE ( epsr, epsi, eels, ieps)
    !
    ! local cleaning
    !
    DEALLOCATE (  dipole, dipole_aux )

    CALL deallocate_bec_type(becp)
    
END SUBROUTINE eps_calc


!----------------------------------------------------------------------------------------
SUBROUTINE jdos_calc ( smeartype, intersmear, nbndmin, nbndmax, shift, nspin )
  !--------------------------------------------------------------------------------------
  ! Computes the joint density of states 
  ! 
  USE kinds,                ONLY : DP
  USE constants,            ONLY : PI, RYTOEV
  USE wvfct,                ONLY : nbnd, et
  USE klist,                ONLY : nks
  USE io_global,            ONLY : ionode, stdout
  USE grid_module,          ONLY : alpha, focc, nw, wgrid
  USE eps_writer,           ONLY : eps_writetofile
  !
  IMPLICIT NONE

  !
  ! input variables
  !
  INTEGER,      INTENT(IN) :: nbndmin, nbndmax, nspin
  REAL(DP),     INTENT(IN) :: intersmear, shift
  CHARACTER(*), INTENT(IN) :: smeartype
  !
  ! local variables
  !
  INTEGER  :: ik, is, iband1, iband2
  INTEGER  :: iw, ierr
  REAL(DP) :: etrans, w, renorm, count, srcount(0:1), renormzero,renormuno
  !
  CHARACTER(128)        :: desc
  REAL(DP), ALLOCATABLE :: jdos(:),srjdos(:,:)

  !
  !--------------------------
  ! main routine body
  !--------------------------
  !
  ! No wavefunctions are needed in order to compute jdos, only eigenvalues,
  ! they are distributed to each task so
  ! no mpi calls are necessary in this routine
  !

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
  renorm = alpha * sum( jdos(:) )

  !
  ! write results on data files
  !
  IF (ionode) THEN
      WRITE(stdout,"(/,5x, 'Integration over JDOS gives: ',f30.9,' instead of 1.0d0' )") renorm
      WRITE(stdout,"(/,5x, 'Writing output on file...' )")
      !
      desc = "energy grid [eV]     JDOS [1/eV]"
      CALL eps_writetofile('jdos',desc,nw,wgrid,1,jdos)
      !
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
  renormzero = alpha * sum( srjdos(0,:) )
  renormuno  = alpha * sum( srjdos(1,:) )

  !
  ! write results on data files
  !
  IF (ionode) THEN
      !
      WRITE(stdout,"(/,5x, 'Integration over spin UP JDOS gives: ',f30.9,' instead of 1.0d0' )") renormzero
      WRITE(stdout,"(/,5x, 'Integration over spin DOWN JDOS gives: ',f30.9,' instead of 1.0d0' )") renormuno
      WRITE(stdout,"(/,5x, 'Writing output on file...' )")
      !
      desc = "energy grid [eV]     UJDOS [1/eV]      DJDOS[1/eV]"
      CALL eps_writetofile('jdos',desc,nw,wgrid,2,srjdos(0:1,:))
      !
  ENDIF

  DEALLOCATE ( srjdos )
ENDIF

END SUBROUTINE jdos_calc

!-----------------------------------------------------------------------------
SUBROUTINE offdiag_calc ( intersmear, intrasmear, nbndmin, nbndmax, shift, metalcalc, nspin )
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
  USE mp_pools,             ONLY : inter_pool_comm
  USE mp,                   ONLY : mp_sum
  USE grid_module,          ONLY : focc, nw, wgrid

  !
  IMPLICIT NONE

  !
  ! input variables
  !
  INTEGER,      INTENT(IN) :: nbndmin, nbndmax, nspin
  REAL(DP),     INTENT(IN) :: intersmear, intrasmear, shift
  LOGICAL,      INTENT(IN) :: metalcalc
  !
  ! local variables
  !
  INTEGER  :: ik, iband1, iband2
  INTEGER  :: iw, ierr, it1, it2
  REAL(DP) :: etrans, const, w
  !
  COMPLEX(DP), ALLOCATABLE :: dipole_aux(:,:,:)
  COMPLEX(DP), ALLOCATABLE :: epstot(:,:,:),dipoletot(:,:,:,:)

  !
  !--------------------------
  ! main routine body
  !--------------------------
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
     CALL dipole_calc_2( ik, dipole_aux, metalcalc, nbndmin, nbndmax)
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
         WRITE(41,"(4f30.6)") wgrid(iw), REAL(epstot(1,1, iw)), aimag(epstot(1,1, iw))
         WRITE(42,"(4f30.6)") wgrid(iw), REAL(epstot(1,2, iw)), aimag(epstot(1,2, iw))
         WRITE(43,"(4f30.6)") wgrid(iw), REAL(epstot(1,3, iw)), aimag(epstot(1,3, iw))
         WRITE(44,"(4f30.6)") wgrid(iw), REAL(epstot(2,1, iw)), aimag(epstot(2,1, iw))
         WRITE(45,"(4f30.6)") wgrid(iw), REAL(epstot(2,2, iw)), aimag(epstot(2,2, iw))
         WRITE(46,"(4f30.6)") wgrid(iw), REAL(epstot(2,3, iw)), aimag(epstot(2,3, iw))
         WRITE(47,"(4f30.6)") wgrid(iw), REAL(epstot(3,1, iw)), aimag(epstot(3,1, iw))
         WRITE(48,"(4f30.6)") wgrid(iw), REAL(epstot(3,2, iw)), aimag(epstot(3,2, iw))
         WRITE(49,"(4f30.6)") wgrid(iw), REAL(epstot(3,3, iw)), aimag(epstot(3,3, iw))
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
  DEALLOCATE ( dipoletot, dipole_aux, epstot )

END SUBROUTINE offdiag_calc


!--------------------------------------------------------------------
SUBROUTINE dipole_calc( ik, dipole_aux, metalcalc, nbndmin, nbndmax )
  !------------------------------------------------------------------
  !
  USE kinds,                ONLY : DP
  USE wvfct,                ONLY : nbnd, npwx
  USE wavefunctions,        ONLY : evc
  USE klist,                ONLY : xk, ngk, igk_k
  USE gvect,                ONLY : ngm, g
  USE io_files,             ONLY : restart_dir
  USE pw_restart_new,       ONLY : read_collected_wfc
  USE grid_module,          ONLY : focc, full_occ
  USE mp_bands,             ONLY : intra_bgrp_comm
  USE mp,                   ONLY : mp_sum
  USE lsda_mod,             ONLY : nspin
  !
  IMPLICIT NONE
  !
  ! global variables
  INTEGER,     INTENT(IN)    :: ik,nbndmin,nbndmax
  COMPLEX(DP), INTENT(INOUT) :: dipole_aux(3,nbnd,nbnd)
  LOGICAL,     INTENT(IN)    :: metalcalc
  !
  ! local variables
  INTEGER     :: iband1,iband2,ig,npw
  COMPLEX(DP) :: caux

  !
  ! Routine Body
  !
  CALL start_clock( 'dipole_calc' )
  !
  ! read wfc for the given kpt
  !
  CALL read_collected_wfc ( restart_dir(), ik, evc )
  !
  ! compute matrix elements
  !
  dipole_aux(:,:,:) = (0.0_DP,0.0_DP)
  !
  npw = ngk(ik)
  !
  DO iband2 = nbndmin,nbndmax
      IF ( focc(iband2,ik) <  full_occ) THEN
          DO iband1 = nbndmin,nbndmax
              !
              IF ( iband1==iband2 ) CYCLE
              IF ( focc(iband1,ik) >= 0.5e-4*full_occ ) THEN
                  !
                  DO ig=1,npw
                      !
                      caux= conjg(evc(ig,iband1))*evc(ig,iband2)
                      !
                      ! Non collinear case
                      IF ( nspin == 4 ) THEN
                          caux = caux + conjg(evc(ig+npwx,iband1))*evc(ig+npwx,iband2)
                      ENDIF
                      !
                      dipole_aux(:,iband1,iband2) = dipole_aux(:,iband1,iband2) + &
                            ( g(:,igk_k(ig,ik)) ) * caux
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
  ! tensor directly as we need it for the intraband term, without interference with interband one.
  !
  IF (metalcalc) THEN
     !
     DO iband1 = nbndmin,nbndmax
        DO  ig=1,npw
          !
          caux= conjg(evc(ig,iband1))*evc(ig,iband1)
          !
          ! Non collinear case
          IF ( nspin == 4 ) THEN
              caux = caux + conjg(evc(ig+npwx,iband1))*evc(ig+npwx,iband1)
          ENDIF
          !
          dipole_aux(:,iband1,iband1) = dipole_aux(:,iband1,iband1) + &
                                        ( g(:,igk_k(ig,ik))+ xk(:,ik) ) * caux
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


!--------------------------------------------------------------------
SUBROUTINE dipole_calc_2( ik, dipole_aux, metalcalc, nbndmin, nbndmax )
  !------------------------------------------------------------------
  !
  ! Calculates the dipole matrix elements including nonlocal contribution
  ! of pseudopotential
  ! 
  USE kinds,                ONLY : DP
  USE wvfct,                ONLY : nbnd, npwx, et
  USE wavefunctions,        ONLY : evc
  USE klist,                ONLY : xk, ngk, igk_k, wk
  USE gvect,                ONLY : ngm, g
  USE io_global,            ONLY : stdout, ionode
  USE io_files,             ONLY : restart_dir
  USE pw_restart_new,       ONLY : read_collected_wfc
  USE grid_module,          ONLY : focc, full_occ
  USE mp_bands,             ONLY : intra_bgrp_comm
  USE mp,                   ONLY : mp_sum

  USE uspp,                 ONLY : nkb, vkb, okvan
  USE becmod,               ONLY : bec_type, calbec, allocate_bec_type, &
                                   becp,  deallocate_bec_type  
  USE noncollin_module,     ONLY : npol, noncolin
  USE lsda_mod,             ONLY : nspin, isk, current_spin
  USE uspp_init,            ONLY : init_us_2  

  !
  IMPLICIT NONE
  !
  ! global variables
  INTEGER,        INTENT(IN)      :: ik
  !! kpt index
  INTEGER,        INTENT(IN)      :: nbndmin
  !! lowest specified band for calculation
  INTEGER,        INTENT(IN)      :: nbndmax
  !! highest specified band for calculation
  COMPLEX(DP), INTENT(INOUT) :: dipole_aux(3,nbnd,nbnd)
  !! auxiliary variable containing the dipole matrix elements
  LOGICAL,     INTENT(IN)    :: metalcalc
  !! flag for if we are doing a metal calculation
  !
  ! local variables
  INTEGER     :: iband1,iband2,ig,npw,ipol
  COMPLEX(DP) :: caux
  COMPLEX(DP), ALLOCATABLE    :: ppsi(:,:), ppsi_us(:,:)
  !
  ! Routine Body
  !
  CALL start_clock( 'dipole_calc' )
  !
  ! read wfc for the given kpt
  !
  CALL read_collected_wfc ( restart_dir(), ik, evc )
  !
  ! compute matrix elements
  !
  dipole_aux(:,:,:) = (0.0_DP,0.0_DP)
  !
  npw = ngk(ik)
  !
  IF (nspin == 2) THEN
      current_spin = isk(ik)
  ELSE
      current_spin = 1
  END IF
  !
  ALLOCATE(ppsi(npwx*npol,nbnd))
  IF (okvan) ALLOCATE(ppsi_us(npwx*npol,nbnd))
  !
  CALL init_us_2 (npw, igk_k(1,ik), xk(1,ik), vkb)
  ! 
  CALL calbec (npw, vkb, evc, becp)

  DO ipol=1,3
    CALL compute_ppsi(ppsi, ppsi_us, ik, ipol, nbnd, current_spin)
    !
    !    DO iband2 = 1, nbnd
    DO iband2 = nbndmin, nbndmax
      !
      IF ( focc(iband2,ik) <  full_occ) THEN
         !        DO iband1 = 1, nbnd
         DO iband1 = nbndmin, nbndmax
          IF ( iband1==iband2 ) CYCLE
          !
          IF ( focc(iband1,ik) >= 0.5e-4*full_occ ) THEN
            DO ig = 1, npw
              dipole_aux(ipol,iband1,iband2) = dipole_aux(ipol,iband1,iband2) +&
                                             & conjg(evc(ig,iband1))*ppsi(ig,iband2)
              IF (okvan) THEN
                dipole_aux(ipol,iband1,iband2) = dipole_aux(ipol,iband1,iband2) +&
                     & conjg(evc(ig,iband1))*(0.d0,0.5d0)*(et(iband1,ik) - et(iband2,ik))*ppsi_us(ig,iband2)
              END IF

              !
              ! Non collinear case
              !
              IF ( noncolin ) THEN
                dipole_aux(ipol,iband1,iband2) = dipole_aux(ipol,iband1,iband2) +&
                     & conjg(evc(ig+npwx,iband1))*ppsi(ig+npwx,iband2)
                IF (okvan) THEN
                  dipole_aux(ipol,iband1,iband2) = dipole_aux(ipol,iband1,iband2) +&
                       & conjg(evc(ig+npwx,iband1))*(0.d0,0.5d0)*(et(iband1,ik) -&
                       & et(iband2,ik))*ppsi_us(ig+npwx,iband2)
                END IF
              ENDIF
              !
            ENDDO
          ENDIF
          !
        ENDDO
      ENDIF
      !
    ENDDO

    !
    ! The diagonal terms are taken into account only if the system is treated like a metal, not
    ! in the intraband therm. Because of this we can recalculate the diagonal component of the dipole
    ! tensor directly as we need it for the intraband term, without interference with interband one.
    !
    IF (metalcalc) THEN
      DO iband1 = nbndmin, nbndmax
        DO ig = 1, npw
          dipole_aux(ipol,iband1,iband1) = dipole_aux(ipol,iband1,iband1) +&
                                         & conjg(evc(ig,iband1))*ppsi(ig,iband1)
          !
          ! Non collinear case
          !
          IF ( noncolin ) THEN
            dipole_aux(ipol,iband1,iband1) = dipole_aux(ipol,iband1,iband1) +&
                                           & conjg(evc(ig+npwx,iband1))*ppsi(ig+npwx,iband1)
          ENDIF
          !
        ENDDO
      ENDDO
    ENDIF
  ENDDO

  !
  DEALLOCATE(ppsi)
  IF (okvan) DEALLOCATE(ppsi_us)  
  !
  ! recover over G parallelization (intra_bgrp)
  !
  CALL mp_sum( dipole_aux, intra_bgrp_comm )
  !
  CALL stop_clock( 'dipole_calc' )
  !
END SUBROUTINE dipole_calc_2

SUBROUTINE dipole_symmetrize (ik, dipole_aux, dipole_sqr,nbndmin,nbndmax)
  !
  ! calculates the symmetrized sum of the squared dipole matrix elements
  !
  USE kinds,                ONLY : DP
  USE klist,                ONLY : wk, nkstot
  USE wvfct,                ONLY : nbnd
  USE cell_base,            ONLY : at, bg
  USE symm_base,            ONLY : s, nsym
  USE grid_module,          ONLY : full_occ

  IMPLICIT NONE

  INTEGER,        INTENT(IN)      :: ik
  !! kpt index
  INTEGER,        INTENT(IN)      :: nbndmin
  !! lowest specified band for calculation
  INTEGER,        INTENT(IN)      :: nbndmax
  !! highest specified band for calculation
  COMPLEX(DP),    INTENT(INOUT)   :: dipole_aux(3,nbnd,nbnd)
  !! auxiliary variable containing the dipole matrix elements
  REAL(DP),       INTENT(INOUT)   :: dipole_sqr(3,nbnd,nbnd)
  !! contains the symmetrized and squared dipole matrix elements
  !
  INTEGER     :: iband1
  !! loop band index
  INTEGER     :: iband2
  !! second loop band index
  INTEGER     :: isym
  !! loop index over symmetries
  COMPLEX(DP) :: work(3)
  !! temporary variable for dipole matrix element
  COMPLEX(DP) :: dipole_rot(3)
  !! rotated dipole matrix el via symmetry op
  dipole_sqr(:,:,:) = 0.0_DP
  !
  DO iband2 = nbndmin, nbndmax
    DO iband1 = nbndmin, nbndmax
      !
      DO isym = 1, nsym
        ! bring dipole matrix elements to crystal axis
        !
        work(:) = dipole_aux(1,iband1,iband2) * at(1,:) + &
                  dipole_aux(2,iband1,iband2) * at(2,:) + &
                  dipole_aux(3,iband1,iband2) * at(3,:)
        !
        dipole_rot(:) = s ( :, 1, isym ) * work(1) + &
                        s ( :, 2, isym ) * work(2) + &
                        s ( :, 3, isym ) * work(3) 
        !
        ! bring back to cartesian axis
        !
        work(:) = dipole_rot(1) * bg(:,1) + &
                  dipole_rot(2) * bg(:,2) + &
                  dipole_rot(3) * bg(:,3)
        !
        dipole_sqr(:,iband1,iband2) = dipole_sqr(:,iband1,iband2) + &
                                      REAL(work(:)*conjg(work(:)), DP) / &
                                      REAL(nsym, DP) * wk(ik) 
      ENDDO
    ENDDO
  ENDDO
  !
END SUBROUTINE dipole_symmetrize
