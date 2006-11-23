
! Copyright (C) 2004 PWSCF group 
! This file is distributed under the terms of the 
! GNU General Public License. See the file `License' 
! in the root directory of the present distribution, 
! or http://www.gnu.org/copyleft/gpl.txt . 

!------------------------------
PROGRAM epsilon
!------------------------------
  !
  ! Compute the complex macroscopic dielectric function,
  ! at the RPA level, neglecting local field effects.
  ! Eps is computed both on the real or immaginary axis
  !
  ! Authors: Andrea Benassi, Andrea Ferretti
  !
  ! NOTE: Part of the basic implementation is taken from pw2gw.f90;
  !
  !
  USE kinds,     ONLY : DP
  USE io_files,  ONLY : nd_nmbr, prefix, outdir, tmp_dir
  USE io_global, ONLY : stdout, ionode, ionode_id
  USE mp,        ONLY : mp_bcast
  ! 
  IMPLICIT NONE

  !
  ! input variables
  !
  INTEGER            :: nw
  REAL(DP)           :: smear,wmax,shift
  CHARACTER(10)           :: calculation,smeartype
  !
  NAMELIST / inputpp / prefix, outdir, calculation
  NAMELIST / energy_grid / smeartype, smear , wmax, nw, shift  
  
  !
  ! local variables
  !
  INTEGER :: ios



!--------------------------------------------- 
! program body
!--------------------------------------------- 
!

  CALL start_postproc(nd_nmbr)

  !
  ! Set default values for variables in namelist 
  !
  calculation  = 'eps'
  prefix       = 'pwscf'
  shift        = 0.0d0
  outdir       = './'
  smear        = 0.02
  wmax         = 30.0d0
  nw           = 600
  smeartype    = 'drude'

  !
  ! read input
  !
  CALL input_from_file( )
  IF ( ionode )  THEN 
     !
     READ (5, inputpp, IOSTAT=ios)
     IF (ios/=0) CALL errore('epsilon', 'reading namelist INPUTPP', ABS(ios))
     !
     READ (5, energy_grid, IOSTAT=ios)
     IF (ios/=0) CALL errore('epsilon', 'reading namelist ENERGY_GRID', ABS(ios))
     !
     tmp_dir = TRIM(outdir)
  ENDIF

  ! 
  ! ... Broadcast variables 
  !
  CALL mp_bcast( smeartype, ionode_id ) 
  CALL mp_bcast( calculation, ionode_id )
  CALL mp_bcast( prefix, ionode_id ) 
  CALL mp_bcast( tmp_dir, ionode_id ) 
  CALL mp_bcast( shift, ionode_id ) 
  CALL mp_bcast( outdir, ionode_id ) 
  CALL mp_bcast( smear, ionode_id ) 
  CALL mp_bcast( wmax, ionode_id ) 
  CALL mp_bcast( nw, ionode_id ) 
  
  !
  ! ... start postproc
  !
  CALL read_file
  CALL openfil_pp

  !
  ! ... run the specific pp calculation
  !
  CALL eps_calc ( smeartype, smear, nw, wmax, shift, calculation )
  

  CALL stop_pp 

END PROGRAM epsilon 


!--------------------------------------------------------------
SUBROUTINE eps_calc ( smeartype, smear, nw, wmax, shift, calculation )
  !--------------------------------------------------------------
  !
  USE kinds,      ONLY : DP
  USE constants,  ONLY : PI, RYTOEV
  USE cell_base,  ONLY : alat, tpiba2, at, omega
  USE wvfct,      ONLY : npw, npwx, nbnd, gamma_only, igk, g2kin, wg, et
  USE gvect,      ONLY : ngm, g, gg, ecutwfc
  USE klist,      ONLY : nks, nkstot, xk, wk, nelec
  USE lsda_mod,   ONLY : nspin
  USE io_files,   ONLY : nwordwfc, iunwfc
  USE io_global,  ONLY : ionode, ionode_id, stdout
  USE uspp,       ONLY : okvan
  USE wavefunctions_module, ONLY : evc
  ! 
  IMPLICIT NONE

  !
  ! input variables
  !
  INTEGER,       INTENT(IN) :: nw
  REAL(DP),      INTENT(IN) :: wmax, smear, shift
  CHARACTER (10),          INTENT(IN) :: calculation, smeartype
  !
  ! local variables
  !
  INTEGER       :: i, ig, ik, iband1, iband2,it1,it2
  INTEGER       :: iw, iwp, ierr,count
  REAL(DP)      :: alpha, etrans, const, w, renorm 
  COMPLEX(DP)   :: caux, dipole_aux(3),dipole_aux1,dipole_aux2
  !
  INTEGER,  ALLOCATABLE    :: igk_l2g(:)
  REAL(DP), ALLOCATABLE    :: focc(:,:), wgrid(:)
  REAL(DP), ALLOCATABLE    :: epsr(:,:), epsi(:,:)
  REAL(DP), ALLOCATABLE    :: jdos(:), ieps(:,:), eels(:,:)
  REAL(DP), ALLOCATABLE    :: dipole(:,:,:)
  COMPLEX(DP), ALLOCATABLE :: epstot(:,:,:),dipoletot(:,:,:,:) 
!
!--------------------------
! main routine body
!--------------------------
!

  !
  ! check on the number of bands: we need to include empty bands in order to allow
  ! to write the transitions
  !
  IF (nbnd  <= nelec/2 ) CALL errore('epsilon', 'ban band number', 1)

  !
  ! spin is not implemented
  !
  IF( nspin > 1 ) CALL errore('epsilon','Spin polarization not implemented',1)

  !
  ! USPP are not implemented (dipole matrix elements are not trivial at all)
  !
  IF ( okvan ) CALL errore('epsilon','USPP are not implemented',1)
 
  !
  ! pool parallelization not implemented
  !
  IF ( nkstot /= nks ) CALL errore('epsilon','pool are not implemented',1)

!
! start different kinds of calculations
!

 IF (calculation == 'jdos') THEN

!
! initializations
!
  !   
  ! allocate main spectral and auxiliary quantities   
  !   
  ALLOCATE (focc(nbnd,nks), STAT=ierr)
      IF (ierr/=0) CALL errore('epsilon','allocating focc',ABS(ierr))
  ALLOCATE(wgrid(nw+1),STAT=ierr)
      IF (ierr/=0) CALL errore('epsilon','allocating wgrid',ABS(ierr))
  ALLOCATE( jdos(nw+1), STAT=ierr )
      IF (ierr/=0) CALL errore('epsilon','allocating eps',ABS(ierr))

  !
  ! occupation numbers
  !
  DO ik = 1,nks
  DO i  = 1,nbnd
       focc(i,ik)= wg(i,ik) * 2.0_DP/wk(ik)
  ENDDO
  ENDDO

  !
  ! set the energy grid
  !
  alpha = wmax/REAL(nw, DP)
  !
  DO iw = 1, nw + 1
      wgrid(iw) = -wmax/2.0_DP + (iw-1)*alpha
  ENDDO
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

 kpt_loop1_lor: &
 DO ik = 1, nks  

     ! 
     ! Calculation of joint density of states  
     ! 'smear' is the brodening parameter  
     ! 
     DO iband2 = 1,nbnd
     DO iband1 = 1,nbnd
         !
         IF ( focc(iband2,ik) <  2.0d0) THEN
         IF ( focc(iband1,ik) >= 1e-4 ) THEN
               !
               ! transition energy
               !
               etrans = ( et(iband2,ik) -et(iband1,ik) ) * RYTOEV 
               !
               if( etrans < 1d-10 ) cycle

               ! loop over frequencies
               !

               count=count+ (focc(iband1,ik)-focc(iband2,ik))

               DO iw = 1, nw+1 
                   !
                   w = wgrid(iw)
                   !
                   jdos(iw) = jdos(iw) + smear * (focc(iband1,ik)-focc(iband2,ik)) &
                                / ( PI * ( (etrans -w )**2 + (smear)**2 ) )
               
               ENDDO
     
         ENDIF
         ENDIF

     ENDDO
     ENDDO

  ENDDO kpt_loop1_lor

 ELSE IF (smeartype=='gauss') THEN

 kpt_loop1_gauss: &
 DO ik = 1, nks  

     ! 
     ! Calculation of joint density of states  
     ! 'smear' is the brodening parameter  
     ! 
     DO iband2 = 1,nbnd
     DO iband1 = 1,nbnd
         !
         IF ( focc(iband2,ik) <  2.0d0) THEN
         IF ( focc(iband1,ik) >= 1e-4 ) THEN
               !
               ! transition energy
               !
               etrans = ( et(iband2,ik) -et(iband1,ik) ) * RYTOEV 
               !
               if( etrans < 1d-10 ) cycle

               ! loop over frequencies
               !

               count=count+ (focc(iband1,ik)-focc(iband2,ik))

               DO iw = 1, nw+1 
                   !
                   w = wgrid(iw)
                   !
                   jdos(iw) = jdos(iw) + (focc(iband1,ik)-focc(iband2,ik)) * &
                              EXP(-(etrans-w)**2/smear**2) &
                                / (smear * SQRT(PI)) 
               
               ENDDO
     
         ENDIF
         ENDIF

     ENDDO
     ENDDO

  ENDDO kpt_loop1_gauss

 ELSE 

 CALL errore('epsilon', 'reading SMEARTYPE in namelist ENERGY_GRID', 1)
 
 ENDIF
  
  !
  ! jdos normalizzation
  ! 

  jdos(:)=jdos(:)/count

  !
  ! check jdos normalization
  !

  renorm = alpha * SUM( jdos(:) )

IF (ionode) THEN
     
      WRITE(stdout,*)
      WRITE(stdout,*)'-----------------------------------------------------------------'
      WRITE(stdout,*)
      WRITE(stdout,*)'Integration over JDOS gives: ', renorm ,' instead of 1.0d0'   
      WRITE(stdout,*)
      WRITE(stdout,*)'-----------------------------------------------------------------'
      WRITE(stdout,*)

     !
     ! write results on data files
     !

     OPEN (30, FILE='jdos.dat', FORM='FORMATTED' )
     !
     WRITE(30, "(2x,'# energy grid [eV]     JDOS [1/eV] ')" )
     !
     DO iw =1, nw+1
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
  DEALLOCATE ( focc, wgrid )
  DEALLOCATE ( jdos )
 
 ELSE IF (calculation == 'eps' ) THEN

!
! initializations
!
  !   
  ! allocate main spectral and auxiliary quantities   
  !   
  ALLOCATE( dipole(3, nbnd, nbnd), STAT=ierr )
      IF (ierr/=0) CALL errore('epsilon','allocating dipole', ABS(ierr) )
  ALLOCATE (focc(nbnd,nks), STAT=ierr)
      IF (ierr/=0) CALL errore('epsilon','allocating focc', ABS(ierr))
  ALLOCATE(wgrid(nw+1),STAT=ierr)
      IF (ierr/=0) CALL errore('epsilon','allocating wgrid', ABS(ierr))
  ALLOCATE( epsr( 3, nw+1), epsi( 3, nw+1), eels( 3, nw+1), ieps(3,nw+1 ), STAT=ierr )
      IF (ierr/=0) CALL errore('epsilon','allocating eps', ABS(ierr))
 !
  ! occupation numbers
  !
  DO ik = 1,nks
  DO i  = 1,nbnd
       focc(i,ik)= wg(i,ik) * 2.0_DP/wk(ik)
  ENDDO
  ENDDO

  !
  ! set the energy grid
  !
  alpha = wmax/REAL(nw, DP)
  !
  DO iw = 1, nw + 1
      wgrid(iw) = (iw-1)*alpha
  ENDDO
  !
  ! initialize response functions
  !
  epsr(:,:)  = 0.0_DP
  epsi(:,:)  = 0.0_DP
  ieps(:,:)  = 0.0_DP

  !
  ! main kpt loop
  !
  kpt_loop2: &
  DO ik = 1, nks  
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
     dipole   = 0.0_DP
     !
     DO iband2 = 1,nbnd
     DO iband1 = 1,nbnd
         !
         IF ( focc(iband2,ik) <  2.0d0) THEN
         IF ( focc(iband1,ik) >= 1e-4 ) THEN
               !
               dipole_aux = 0.0_DP
               !
               DO  ig=1,npw
                    !
                    caux= CONJG(evc(ig,iband1))*evc(ig,iband2) 
                    !
                    dipole_aux(:) = dipole_aux(:) + &
                          ( g(:,igk(ig)) ) * caux
                   !
               ENDDO
               !
               dipole(:,iband1,iband2)= tpiba2* REAL( dipole_aux(:) * CONJG( dipole_aux(:) ), DP )
               !
         ENDIF
         ENDIF
         !
     ENDDO
     ENDDO

     !
     ! recover G parallelism
     !
     CALL reduce( 3 * nbnd * nbnd, dipole ) 


     ! 
     ! Calculation of real and immaginary parts 
     ! of the macroscopic dielettric function from dipole
     ! approximation. 
     ! 'smear' is the brodening parameter  
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

               ! loop over frequencies
               !
               DO iw = 1, nw+1 
                   !
                   w = wgrid(iw)
                   !

                   epsi(:,iw) = epsi(:,iw) + dipole(:,iband1,iband2) * smear * w* &
                                             RYTOEV**3 * (focc(iband1,ik)-focc(iband2,ik))/  &
                                  (( (etrans**2 -w**2 )**2 + smear**2 * w**2 )* etrans )
                    
                   epsr(:,iw) = epsr(:,iw) + dipole(:,iband1,iband2) * RYTOEV**3 * &
                                             (focc(iband1,ik)-focc(iband2,ik)) * &
                                             (etrans**2 - w**2 ) / &
                                  (( (etrans**2 -w**2 )**2 + smear**2 * w**2 )* etrans )
                   
               ENDDO

         ENDIF
         ENDIF

     ENDDO
     ENDDO

  ENDDO kpt_loop2

  !
  ! impose the correct normalization
  !
  const = 64.0d0 * PI / ( omega * REAL(nkstot, DP) )

  epsr(:,:) = 1.0_DP + epsr(:,:) * const  
  epsi(:,:) =          epsi(:,:) * const 
  !
  ! Calculation of the eels 
  !

  DO iw = 1, nw+1

       eels(:,iw) = epsi(:,iw) / ( epsr(:,iw)**2 + epsi(:,iw)**2 )

  ENDDO
  
  !
  !  calculation of dielectric function on the immaginary frequency axe
  !                   

 DO iw = 1, nw+1
         DO iwp = 2, nw+1 

               ieps(:,iw) = ieps(:,iw) + wgrid(iwp) * epsi(:,iwp) / ( wgrid(iwp)**2 + wgrid(iw)**2) 

         ENDDO
 ENDDO               

  ieps(:,:) = 1.0d0 + 2 / PI * ieps(:,:) * alpha 

  !
  ! check  dielectric function  normalizzation via sumrule  
  !
  
  renorm = alpha* SUM( epsi(1,:) * wgrid(:)) 

  !
  IF ( ionode ) THEN
       WRITE(stdout,*)
       WRITE(stdout,*)'---------------------------------------------------------------------------'
       WRITE(stdout,*)
       WRITE(stdout,*)'Electron density is: ', REAL(nelec, DP) / (omega),' 1/(a.u.)^3'
       WRITE(stdout,*)
       WRITE(stdout,*)'The bulk plasmon frequency is: ', SQRT(renorm * 2.0d0 / PI),' eV'
       WRITE(stdout,*)
       WRITE(stdout,*)'---------------------------------------------------------------------------'
       WRITE(stdout,*) 
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
      DO iw =1, nw+1
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
  DEALLOCATE ( focc, wgrid )
  DEALLOCATE ( epsr, epsi, eels, ieps)

ELSE IF (calculation == 'offdiag' ) THEN
!
! initializations
!
!
! allocate main spectral and auxiliary quantities
!
   ALLOCATE( dipoletot(3,3, nbnd, nbnd), STAT=ierr )
       IF (ierr/=0) CALL errore('epsilon','allocating dipole', ABS(ierr) )
   ALLOCATE (focc(nbnd,nks), STAT=ierr)
       IF (ierr/=0) CALL errore('epsilon','allocating focc', ABS(ierr))
   ALLOCATE(wgrid(nw+1),STAT=ierr)
       IF (ierr/=0) CALL errore('epsilon','allocating wgrid', ABS(ierr))
   ALLOCATE(epstot( 3,3, nw+1),STAT=ierr )
       IF (ierr/=0) CALL errore('epsilon','allocating eps', ABS(ierr))
  !
  ! occupation numbers
  !
   DO ik = 1,nks
   DO i  = 1,nbnd
        focc(i,ik)= wg(i,ik) * 2.0_DP/wk(ik)
   ENDDO
   ENDDO

   !
   ! set the energy grid
   !
   alpha = wmax/REAL(nw, DP)
   !
   DO iw = 1, nw + 1
       wgrid(iw) = (iw-1)*alpha
   ENDDO
   !
   ! initialize response functions
   !
   epstot  = 0.0_DP
   !
   ! main kpt loop
   !
   kpt_loop3: &
   DO ik = 1, nks
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
      dipoletot   = 0.0_DP
      !
      DO it1 = 1, 3
         DO it2 = 1, 3
            DO iband2 = 1,nbnd
               DO iband1 = 1,nbnd
                  !
                  IF ( focc(iband2,ik) <  2.0d0) THEN
                  IF ( focc(iband1,ik) >= 1e-4 ) THEN
                  !
                  dipole_aux1 = 0.0_DP
                  dipole_aux2 = 0.0_DP
                  !
                  DO  ig=1,npw

                     caux= CONJG(evc(ig,iband1))*evc(ig,iband2)

                     dipole_aux1 = dipole_aux1+ &
                           ( g(it1,igk(ig)) ) * caux

                     dipole_aux2 = dipole_aux2+ &
                           ( g(it2,igk(ig)) ) * caux

                 ENDDO
                 !
                 dipoletot(it1,it2,iband1,iband2)= tpiba2* &
                                   dipole_aux1 * CONJG( dipole_aux2 )
                 !
                 ENDIF
                 ENDIF
                 !
              ENDDO
           ENDDO
        ENDDO
     ENDDO
     !
     ! recover G parallelism
     !
     CALL reduce( 2 * 3 * 3 * nbnd * nbnd, dipoletot )
     
     !
     ! Calculation of real and immaginary parts
     ! of the macroscopic dielettric function from dipole
     ! approximation.
     ! 'smear' is the brodening parameter
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

               ! loop over frequencies
               !
              DO iw = 1, nw+1
                   !
                   w = wgrid(iw)
                   !
                   epstot(:,:,iw) = epstot(:,:,iw) + dipoletot(:,:,iband1,iband2)*RYTOEV**3/(etrans) *&
                                  focc(iband1,ik)/(etrans**2 - w**2 - (0,1)*smear*w) 
              
              ENDDO

         ENDIF
         ENDIF

     ENDDO
     ENDDO
  ENDDO kpt_loop3

  !
  ! impose the correct normalization
  !
  const = 64.0d0 * PI / ( omega * REAL(nkstot, DP) )

  epstot(:,:,:) = 1.0_DP + epstot(:,:,:) * const

  !
  ! write results on data files
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
      DO iw =1, nw+1
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
  ! local cleaning
  !
  DEALLOCATE ( focc, wgrid )
  DEALLOCATE ( epstot )
  DEALLOCATE ( dipoletot)

 ELSE

 CALL errore('epsilon', 'reading CALCULATION in namelist INPUTPP', 1)
 
 END IF

END SUBROUTINE eps_calc 

