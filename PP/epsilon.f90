
! Copyright (C) 2004 PWSCF group 
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

CONTAINS

!------------------------------------
  SUBROUTINE grid_build(nw,wmax,wmin)
  !----------------------------------
  !
  USE kinds,    ONLY : DP
  USE wvfct,    ONLY : nbnd, wg
  USE klist,    ONLY : nks, wk, nelec
  USE lsda_mod, ONLY : nspin
  USE uspp,     ONLY : okvan
  !
  IMPLICIT NONE
  !
  ! input vars
  INTEGER,  INTENT(IN) :: nw
  REAL(DP), INTENT(IN) :: wmax ,wmin
  !
  ! local vars
  INTEGER :: iw,ik,i,ierr 
 
  !
  ! check on the number of bands: we need to include empty bands in order to allow
  ! to write the transitions
  !
  IF ( REAL(nbnd, DP)  <= nelec / 2.0_DP ) CALL errore('epsilon', 'ban band number', 1)

  !
  ! spin is not implemented
  !
  IF( nspin > 1 ) CALL errore('grid_build','Spin polarization not implemented',1)

  !
  ! USPP are not implemented (dipole matrix elements are not trivial at all)
  !
  IF ( okvan ) CALL errore('grid_build','USPP are not implemented',1)

  ALLOCATE ( focc( nbnd, nks), STAT=ierr )
  IF (ierr/=0) CALL errore('grid_build','allocating focc', ABS(ierr))
  !
  ALLOCATE( wgrid( nw ), STAT=ierr )
  IF (ierr/=0) CALL errore('grid_build','allocating wgrid', ABS(ierr))

  !
  ! check on k point weights, no symmetry operations are allowed 
  !
  DO ik = 2, nks
     !
     IF ( ABS( wk(1) - wk(ik) ) > 1.0d-8 ) &
        CALL errore('grid_build','non unifrom kpt grid', ik )
     !
  ENDDO
  !
  ! occupation numbers
  !
  DO ik = 1,nks
  DO i  = 1,nbnd
       focc(i,ik)= wg(i, ik ) * 2.0_DP/wk( ik )
  ENDDO
  ENDDO
  !
  ! set the energy grid
  !
  alpha = (wmax - wmin) / REAL(nw, DP)
  !
  DO iw = 1, nw 
      wgrid(iw) = wmin + iw * alpha
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
      CALL errore('grid_destroy','deallocating grid stuff',ABS(ierr))
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
  USE iotk_module
  USE xml_io_base
  USE io_files,    ONLY : nd_nmbr, tmp_dir, prefix, outdir, trimcheck
  USE constants,   ONLY : RYTOEV
  USE ener,        ONLY : ef
  USE klist,       ONLY : lgauss  
  USE ktetra,      ONLY : ltetra
  ! 
  IMPLICIT NONE

  !
  ! input variables
  !
  INTEGER                 :: nw,ierr
  REAL(DP)                :: intersmear,intrasmear,wmax,wmin,shift
  CHARACTER(10)           :: calculation,smeartype
  LOGICAL                 :: metalcalc
  !
  NAMELIST / inputpp / prefix, outdir, calculation
  NAMELIST / energy_grid / smeartype, intersmear ,intrasmear, wmax, wmin, nw, shift  
  
  !
  ! local variables
  !
  INTEGER :: ios

!--------------------------------------------- 
! program body
!--------------------------------------------- 
!

  CALL init_clocks( .TRUE. )
  CALL start_clock( 'epsilon' )
  !
  !
  CALL start_postproc(nd_nmbr)

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
  
  IF ( ionode )  THEN 
     !
     READ (5, inputpp, IOSTAT=ios)
     IF (ios/=0) CALL errore('epsilon', 'reading namelist INPUTPP', ABS(ios))
     !
     READ (5, energy_grid, IOSTAT=ios)
     IF (ios/=0) CALL errore('epsilon', 'reading namelist ENERGY_GRID', ABS(ios))
     !
     tmp_dir = trimcheck(outdir)
     ! 
  ENDIF
  
  ! 
  ! ... Broadcast variables 
  !
  IF (ionode) WRITE( stdout, "( 5x, 'Broadcasting variables...' ) " )

  CALL mp_bcast( smeartype, ionode_id ) 
  CALL mp_bcast( calculation, ionode_id )
  CALL mp_bcast( prefix, ionode_id ) 
  CALL mp_bcast( tmp_dir, ionode_id ) 
  CALL mp_bcast( shift, ionode_id ) 
  CALL mp_bcast( outdir, ionode_id ) 
  CALL mp_bcast( intrasmear, ionode_id )
  CALL mp_bcast( intersmear, ionode_id) 
  CALL mp_bcast( wmax, ionode_id ) 
  CALL mp_bcast( wmin, ionode_id ) 
  CALL mp_bcast( nw, ionode_id ) 

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

  IF (lgauss .OR. ltetra) THEN 
      metalcalc=.true. 
      IF (ionode) WRITE( stdout, "( 5x, 'The system is a metal...' ) " )
  ELSE
      IF (ionode) WRITE( stdout, "( 5x, 'The system is a dielectric...' ) " )
  ENDIF

  !
  ! ... run the specific pp calculation
  !
  IF (ionode) WRITE(stdout,"(/, 5x, 'Performing ',a,' calculation...')") TRIM(calculation)

  CALL start_clock( 'calculation' )
  !
  SELECT CASE ( TRIM(calculation) )
  !
  CASE ( 'eps' )
      !
      CALL eps_calc ( intersmear, intrasmear, nw, wmax, wmin, shift, metalcalc )
      !
  CASE ( 'jdos' )
      !
      CALL jdos_calc ( smeartype, intersmear, intrasmear, nw, wmax, wmin, shift )
      !
  CASE ( 'offdiag' )
      !
      CALL offdiag_calc ( intersmear, intrasmear, nw, wmax, wmin, shift, metalcalc )
      !
  CASE DEFAULT
      !
      CALL errore('epsilon','invalid CALCULATION = '//TRIM(calculation),1)
      !
  END SELECT
  !
  CALL stop_clock( 'calculation' )
  
  !
  ! few info about timing
  !
  CALL stop_clock( 'epsilon' )
  !
  IF ( ionode ) WRITE( stdout , "(/)" )
  !
  CALL print_clock( 'epsilon' )
  CALL print_clock( 'calculation' )
  CALL print_clock( 'dipole_calc' )
  !
  IF ( ionode ) WRITE( stdout, *  )

  !
  !
  CALL stop_pp ()

END PROGRAM epsilon 


!-----------------------------------------------------------------------------
SUBROUTINE eps_calc ( intersmear,intrasmear, nw, wmax, wmin, shift, &
                      metalcalc )
  !-----------------------------------------------------------------------------
  !
  USE kinds,                ONLY : DP
  USE constants,            ONLY : PI, RYTOEV
  USE cell_base,            ONLY : tpiba2, omega
  USE wvfct,                ONLY : nbnd, wg, et
  USE ener,                 ONLY : efermi => ef
  USE klist,                ONLY : nks, nkstot, nelec, degauss
  USE io_global,            ONLY : ionode, ionode_id, stdout
  !
  USE grid_module,          ONLY : alpha, focc, wgrid, grid_build, grid_destroy                
  ! 
  IMPLICIT NONE

  !
  ! input variables
  !
  INTEGER,         INTENT(IN) :: nw
  REAL(DP),        INTENT(IN) :: wmax, wmin, intersmear,intrasmear, shift
  LOGICAL,         INTENT(IN) :: metalcalc 
  !
  ! local variables
  !
  INTEGER       :: i, ig, ik, iband1, iband2 
  INTEGER       :: iw, iwp, ierr
  REAL(DP)      :: etrans, const, w, renorm(3)
  COMPLEX(DP)   :: caux 
  !
  REAL(DP), ALLOCATABLE    :: epsr(:,:), epsi(:,:)
  REAL(DP), ALLOCATABLE    :: ieps(:,:), eels(:,:)
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
  CALL grid_build(nw,wmax,wmin)
  !   
  ! allocate main spectral and auxiliary quantities   
  !   
  ALLOCATE( dipole(3, nbnd, nbnd), STAT=ierr )
  IF (ierr/=0) CALL errore('epsilon','allocating dipole', ABS(ierr) )
  !
  ALLOCATE( dipole_aux(3, nbnd, nbnd), STAT=ierr )
  IF (ierr/=0) CALL errore('epsilon','allocating dipole_aux', ABS(ierr) )
  !
  ALLOCATE( epsr( 3, nw), epsi( 3, nw), eels( 3, nw), ieps(3,nw ), STAT=ierr )
  IF (ierr/=0) CALL errore('epsilon','allocating eps', ABS(ierr))

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
     CALL dipole_calc( ik, dipole_aux)
     !
     dipole(:,:,:)= tpiba2 * REAL( dipole_aux(:,:,:) * CONJG(dipole_aux(:,:,:)), DP ) 

     ! 
     ! Calculation of real and immaginary parts 
     ! of the macroscopic dielettric function from dipole
     ! approximation. 
     ! 'intersmear' is the brodening parameter  
     !
     !Interband
     ! 
     DO iband2 = 1,nbnd
     DO iband1 = 1,nbnd
         !
         IF ( focc(iband2,ik) < 2.0d0) THEN
         IF ( focc(iband1,ik) >= 1e-4 ) THEN

               !
               ! transition energy
               !
               etrans = ( et(iband2,ik) -et(iband1,ik) ) * RYTOEV + shift
               !
               IF( etrans < 1d-10 ) CYCLE

               ! loop over frequencies
               !
               DO iw = 1, nw 
                   !
                   w = wgrid(iw)
                   !

                   epsi(:,iw) = epsi(:,iw) + dipole(:,iband1,iband2) * intersmear * w* &
                                             RYTOEV**3 * (focc(iband1,ik)-focc(iband2,ik))/  &
                                  (( (etrans**2 -w**2 )**2 + intersmear**2 * w**2 )* etrans ) 
                                   
                   epsr(:,iw) = epsr(:,iw) + dipole(:,iband1,iband2) * RYTOEV**3 * &
                                             (focc(iband1,ik)-focc(iband2,ik)) * &
                                             (etrans**2 - w**2 ) / &
                                  (( (etrans**2 -w**2 )**2 + intersmear**2 * w**2 )* etrans ) 
                                    
               ENDDO

         ENDIF
         ENDIF

     ENDDO
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

                   epsi(:,iw) = epsi(:,iw) - 1/2 * dipole(:,iband1,iband1) * intrasmear * w* &
                                RYTOEV**3 * (2*EXP((et(iband1,ik)-efermi)/degauss ))/  &
                    (( w**4 + intrasmear**2 * w**2 )*(1+EXP((et(iband1,ik)-efermi)/ & 
                    degauss))**2*degauss * RYTOEV) 
                                   
                   epsr(:,iw) = epsr(:,iw) + 1/2 * dipole(:,iband1,iband1) * RYTOEV**3 * &
                                             (2*EXP((et(iband1,ik)-efermi)/degauss )) * w**2 / & 
                    (( w**4 + intrasmear**2 * w**2 )*(1+EXP((et(iband1,ik)-efermi)/ &
                    degauss))**2*degauss *RYTOEV )
               ENDDO

         ENDIF
         ENDIF

     ENDDO
     ENDIF
  ENDDO kpt_loop

  !
  ! recover over kpt parallelization (inter_pool)
  !
  CALL poolreduce( 3 * (nw), epsr ) 
  CALL poolreduce( 3 * (nw), epsi ) 
  
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
     renorm(i) = alpha * SUM( epsi(i,:) * wgrid(:) ) 
 ENDDO
  !
  IF ( ionode ) THEN
      !
      WRITE(stdout,"(/,5x, 'The bulk xx plasmon frequency [eV] is: ',f15.9 )")  SQRT(renorm(1) * 2.0d0 / PI)
      WRITE(stdout,"(/,5x, 'The bulk yy plasmon frequency [eV] is: ',f15.9 )")  SQRT(renorm(2) * 2.0d0 / PI)
      WRITE(stdout,"(/,5x, 'The bulk zz plasmon frequency [eV] is: ',f15.9 )")  SQRT(renorm(3) * 2.0d0 / PI)
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

  !
  ! local cleaning
  !
  CALL grid_destroy()
  !
  DEALLOCATE (  dipole, dipole_aux )
  DEALLOCATE ( epsr, epsi, eels, ieps)

END SUBROUTINE eps_calc
 
!-----------------------------------------------------------------------------
SUBROUTINE jdos_calc ( smeartype, intersmear,intrasmear, nw, wmax, wmin, shift )
  !-----------------------------------------------------------------------------
  !
  USE kinds,                ONLY : DP
  USE constants,            ONLY : PI, RYTOEV
  USE cell_base,            ONLY : tpiba2, omega
  USE wvfct,                ONLY : nbnd, et
  USE klist,                ONLY : nks, nkstot, nelec
  USE io_global,            ONLY : ionode, ionode_id, stdout
  USE grid_module,          ONLY : alpha, focc, wgrid, grid_build, grid_destroy                
  ! 
  IMPLICIT NONE

  !
  ! input variables
  !
  INTEGER,         INTENT(IN) :: nw
  REAL(DP),        INTENT(IN) :: wmax, wmin, intersmear,intrasmear, shift
  CHARACTER(*),    INTENT(IN) :: smeartype
  !
  ! local variables
  !
  INTEGER       :: i, ik, iband1, iband2 
  INTEGER       :: iw, iwp, ierr, count
  REAL(DP)      :: etrans, const, w, renorm
  !
  REAL(DP), ALLOCATABLE    :: jdos(:)
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
  CALL grid_build(nw,wmax,wmin)
  !
  ! allocate main spectral and auxiliary quantities
  !
  ALLOCATE( jdos(nw), STAT=ierr )
      IF (ierr/=0) CALL errore('epsilon','allocating jdos',ABS(ierr))
  !
  ! initialize jdos
  !
  jdos(:)=0.0_DP

  ! Initialising a counter for the number of transition
  count=0

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
       DO iband1 = 1,nbnd
           !
           IF ( focc(iband2,ik) <  2.0d0) THEN
           IF ( focc(iband1,ik) >= 1.0d-4 ) THEN
                 !
                 ! transition energy
                 !
                 etrans = ( et(iband2,ik) -et(iband1,ik) ) * RYTOEV
                 !
                 IF( etrans < 1.0d-10 ) CYCLE

                 ! loop over frequencies
                 !
                 count=count+ (focc(iband1,ik)-focc(iband2,ik))

                 DO iw = 1, nw
                     !
                     w = wgrid(iw)
                     !
                     jdos(iw) = jdos(iw) + intersmear * (focc(iband1,ik)-focc(iband2,ik)) &
                                  / ( PI * ( (etrans -w )**2 + (intersmear)**2 ) )

                 ENDDO

           ENDIF
           ENDIF

       ENDDO
       ENDDO

    ENDDO kpt_lor

  ELSE IF (smeartype=='gauss') THEN

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
                 etrans = ( et(iband2,ik) -et(iband1,ik) ) * RYTOEV
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
                                EXP(-(etrans-w)**2/intersmear**2) &
                                  / (intersmear * SQRT(PI))

                 ENDDO

           ENDIF
           ENDIF
       ENDDO
       ENDDO

    ENDDO kpt_gauss

  ELSE

    CALL errore('epsilon', 'invalid SMEARTYPE = '//TRIM(smeartype), 1)

  ENDIF

  !
  ! jdos normalizzation
  !

  jdos(:)=jdos(:)/count

  !
  ! check jdos normalization
  !

  renorm = alpha * SUM( jdos(:) )
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
  CALL grid_destroy()
  DEALLOCATE ( jdos )

END SUBROUTINE jdos_calc

!-----------------------------------------------------------------------------
SUBROUTINE offdiag_calc ( intersmear,intrasmear, nw, wmax, wmin, shift, &
                          metalcalc )
  !-----------------------------------------------------------------------------
  !
  USE kinds,                ONLY : DP
  USE constants,            ONLY : PI, RYTOEV
  USE cell_base,            ONLY : tpiba2, omega
  USE wvfct,                ONLY : nbnd, et
  USE ener,                 ONLY : efermi => ef
  USE klist,                ONLY : nks, nkstot, nelec, degauss
  USE io_global,            ONLY : ionode, ionode_id, stdout
  USE grid_module,          ONLY : alpha, focc, wgrid, grid_build, grid_destroy                
  ! 
  IMPLICIT NONE

  !
  ! input variables
  !
  INTEGER,         INTENT(IN) :: nw
  REAL(DP),        INTENT(IN) :: wmax, wmin, intersmear,intrasmear, shift
  LOGICAL,         INTENT(IN) :: metalcalc 
  !
  ! local variables
  !
  INTEGER       :: i, ig, ik, iband1, iband2 
  INTEGER       :: iw, iwp, ierr,it1,it2
  REAL(DP)      :: etrans, const, w, renorm
  COMPLEX(DP)   :: caux 
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
  CALL grid_build(nw,wmax,wmin)
  !
  ! allocate main spectral and auxiliary quantities
  !
  ALLOCATE( dipoletot(3,3, nbnd, nbnd), STAT=ierr )
  IF (ierr/=0) CALL errore('epsilon','allocating dipoletot', ABS(ierr) )
  !
  ALLOCATE( dipole_aux(3, nbnd, nbnd), STAT=ierr )
  IF (ierr/=0) CALL errore('epsilon','allocating dipole_aux', ABS(ierr) )
  !
  ALLOCATE(epstot( 3,3, nw),STAT=ierr )
  IF (ierr/=0) CALL errore('epsilon','allocating epstot', ABS(ierr))

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
     CALL dipole_calc(ik,dipole_aux)
     !
     DO it2 = 1, 3
        DO it1 = 1, 3
           dipoletot(it1,it2,:,:) = tpiba2 * dipole_aux(it1,:,:) * CONJG( dipole_aux(it2,:,:) )
        ENDDO
     ENDDO
     !
     ! Calculation of real and immaginary parts
     ! of the macroscopic dielettric function from dipole
     ! approximation.
     ! 'intersmear' is the brodening parameter
     !
     DO iband2 = 1,nbnd
     DO iband1 = 1,nbnd
         !
         IF ( focc(iband2,ik) <  2.0d0) THEN
         IF ( focc(iband1,ik) >= 1e-4 ) THEN
             !
             ! transition energy
             !
             etrans = ( et(iband2,ik) -et(iband1,ik) ) * RYTOEV + shift
             !
             if( etrans < 1d-10 ) cycle
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
         ENDIF

     ENDDO
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
                   epstot(:,:,iw) = epstot(:,:,iw) + 1/2 * dipoletot(:,:,iband1,iband1)* &
                                RYTOEV**3 * (2*EXP((et(iband1,ik)-efermi)/degauss ))/  &
                    (( w**2 + (0,1)*intrasmear*w)*(1+EXP((et(iband1,ik)-efermi)/ &
                    degauss))**2*degauss *RYTOEV)
               ENDDO

         ENDIF
         ENDIF

     ENDDO
     ENDIF
  ENDDO

  !
  ! recover over kpt parallelization (inter_pool)
  !
  CALL poolreduce( 2 * 3 * 3 * (nw), epstot )

  !
  ! impose the correct normalization
  !
  const = 64.0d0 * PI / ( omega * REAL(nkstot, DP) )
  !
  epstot(:,:,:) = 1.0_DP + epstot(:,:,:) * const

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
         WRITE(41,"(4f15.6)") wgrid(iw), REAL(epstot(1,1, iw)), AIMAG(epstot(1,1, iw))
         WRITE(42,"(4f15.6)") wgrid(iw), REAL(epstot(1,2, iw)), AIMAG(epstot(1,2, iw))
         WRITE(43,"(4f15.6)") wgrid(iw), REAL(epstot(1,3, iw)), AIMAG(epstot(1,3, iw))
         WRITE(44,"(4f15.6)") wgrid(iw), REAL(epstot(2,1, iw)), AIMAG(epstot(2,1, iw))
         WRITE(45,"(4f15.6)") wgrid(iw), REAL(epstot(2,2, iw)), AIMAG(epstot(2,2, iw))
         WRITE(46,"(4f15.6)") wgrid(iw), REAL(epstot(2,3, iw)), AIMAG(epstot(2,3, iw))
         WRITE(47,"(4f15.6)") wgrid(iw), REAL(epstot(3,1, iw)), AIMAG(epstot(3,1, iw))
         WRITE(48,"(4f15.6)") wgrid(iw), REAL(epstot(3,2, iw)), AIMAG(epstot(3,2, iw))
         WRITE(49,"(4f15.6)") wgrid(iw), REAL(epstot(3,3, iw)), AIMAG(epstot(3,3, iw))
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


!------------------------------------
SUBROUTINE dipole_calc( ik, dipole_aux)
  !----------------------------------
  USE kinds,                ONLY : DP
  USE wvfct,                ONLY : npw, nbnd, igk, g2kin
  USE wavefunctions_module, ONLY : evc
  USE klist,                ONLY : xk
  USE cell_base,            ONLY : tpiba2
  USE gvect,                ONLY : ngm, g, ecutwfc
  USE io_files,             ONLY : nwordwfc, iunwfc
  USE grid_module,          ONLY : focc

IMPLICIT NONE
  !
  ! global variables
  INTEGER, INTENT(IN)        :: ik
  COMPLEX(DP), INTENT(INOUT) :: dipole_aux(3,nbnd,nbnd)
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
  CALL davcio (evc, nwordwfc, iunwfc, ik, - 1) 
  !
  ! compute matrix elements
  !
  dipole_aux(:,:,:) = (0.0_DP,0.0_DP)
  !
  DO iband2 = 1,nbnd
  DO iband1 = 1,nbnd
      !
      IF ( focc(iband2,ik) <  2.0d0) THEN
      IF ( focc(iband1,ik) >= 1e-4 ) THEN
            !
            DO  ig=1,npw
                 !
                 caux= CONJG(evc(ig,iband1))*evc(ig,iband2) 
                 !
                 dipole_aux(:,iband1,iband2) = dipole_aux(:,iband1,iband2) + &
                       ( g(:,igk(ig)) ) * caux
                 !
            ENDDO
      ENDIF
      ENDIF
      !
  ENDDO
  ENDDO
  !
  ! recover over G parallelization (intra_pool)
  !
  CALL reduce( 2 * 3 * nbnd * nbnd, dipole_aux ) 
  !
  CALL stop_clock( 'dipole_calc' )
  !
END SUBROUTINE dipole_calc
