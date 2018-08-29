!---------------------------------------------------------------------
!
! Copyright (C) 2001-2002 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! Generate a pseudopotential in the Virtual Crystal Approximation:
!
!   V^{(vca)} = V_{loc)^{(vca)} + V_{nl}^{(vca)}
! where
!   V_{loc)^{(vca)} = x V_{loc}^{(1)} + (1-x) V_{loc}^{(2)}
! and
!   V_{nl)^{(vca)} = \sum_{ij} |\beta^{(1)}_i> x D^{(1)}_{ij} <\beta^{(1)}_j|
!                  + \sum_{ij} |\beta^{(2)}_i> (1-x)D^{(2)}_{ij} <\beta^{{2)}_j|
! where
!   V_{loc}^{(n)}(r) is the local part of pseudopot n
!   \beta^{{n)}_i(r) are the projectors for pseudopot n
!   D^{(n))_{ij} are the (bare) components of matrix D for pseudopot n
!
!
! virtual_v2: supports reading of UPF v2 format
! Author: Jingyang Wang (jw598@cornell.edu)
! 
!
PROGRAM virtual_test

  USE pseudo_types, ONLY : pseudo_upf, nullify_pseudo_upf, &
                           deallocate_pseudo_upf
  USE upf_module, ONLY : read_upf
  USE write_upf_module, ONLY : write_upf
  USE radial_grids, ONLY : radial_grid_type, nullify_radial_grid
  USE environment, ONLY: environment_start, environment_end
  USE mp_global, ONLY: mp_startup, mp_global_end
  USE io_global, ONLY: ionode, stdout


  !
  IMPLICIT NONE
  !
  INTEGER :: is, iunps, ierr, ounps
  real(8) :: x
  CHARACTER(256) :: filein(2), fileout
  !
  !  Local variables
  !
  INTEGER :: ios
  TYPE (pseudo_upf) :: upf(2), upf_vca
  TYPE (radial_grid_type) :: grid(2)
#if defined(__MPI)
  CALL mp_startup()
#endif
  CALL environment_start('VIRTUAL_V2.X') 
  IF (ionode) THEN 
     !
     PRINT '('' '')'                                                                     
     PRINT '('' Generate the UPF pseudopotential for a virtual atom '')'                 
     PRINT '('' combining two pseudopootentials in UPF format '')'                       
     PRINT '('' '')'                                                                     
     !                                                                                   
                                                                                         
     ! Read pseudopotentials                                                             
     !                                                                                   
     DO is=1,2                                                                           
                                                                                         
       PRINT '(''  Input PP file # '',i2,'' in UPF format > '',$)', is                   
       READ (5, '(a)', end = 20, err = 20) filein(is)                                    
                                                                                         
       !  nullify objects as soon as they are instantiated                               
                                                                                         
       CALL nullify_pseudo_upf(upf(is))                                                  
       CALL nullify_radial_grid(grid(is))                                                
                                                                                         
                                                                                         
                                                                                         
                                                                                         
       CALL read_upf(upf(is), GRID = grid(is), IERR = ierr,FILENAME =  TRIM(filein(is))) 
       !                                                                                 
       IF (ierr/=0 .AND. ierr/=-1) THEN                                                  
          print *, ierr                                                                  
          CALL errore('virtual_test', 'reading pseudo upf', ierr)                        
       END IF                                                                            
       CLOSE (unit=iunps)                                                                
       PRINT '('' '')'                                                                   
                                                                                         
     ENDDO                                                                               
                                                                                         
     ! Choose mixing parameter x                                                         
     !                                                                                   
     PRINT '('' New Pseudo = x '',a,'' + (1-x) '',a)', (trim(filein(is)), is=1,2)        
   10 CONTINUE                                                                           
     WRITE(stdout,'('' mixing parameter x [0<x<1] = '')', advance="NO")                       
     READ (5,*) x                                                                        
     IF (x<0.d0 .or. x>1)  GOTO 10                                                       
                                                                                         
     ! compute virtual crystal approximation                                             
     !                                                                                   
     CALL compute_virtual(x, filein, upf(1), upf_vca)                                    
                                                                                         
     ! write VCA pseudopotential to file                                                 
     !                                                                                   
     fileout='NewPseudo.UPF'                                                             
     PRINT '("Output PP file in UPF format :  ",a)', fileout                             
     !                                                                                   
     CALL write_upf ( TRIM(fileout), upf_vca, SCHEMA='v2')                               
     !                                                                                   
     CLOSE(ounps)                                                                        
     CALL deallocate_pseudo_upf(upf(1))                                                  
     CALL deallocate_pseudo_upf(upf(2))                                                  
     !     ----------------------------------------------------------                    
     WRITE (stdout,"('Pseudopotential successfully written')")                                
     WRITE (stdout,"('Please review the content of the PP_INFO fields')")                     
     WRITE (stdout,"('*** Please TEST BEFORE USING !!! ***')")                                
     !     ----------------------------------------------------------                    
     !
   END IF
  CALL environment_end('VIRTUAL_V2.X')
#if defined(__MPI) 
  CALL mp_global_end()
#endif 

   STOP
20 CALL errore ('virtual.x', 'error reading pseudo file', 1)   

END PROGRAM virtual_test

!
!---------------------------------------------------------------------
SUBROUTINE compute_virtual(x, filein, upf, upf_vca)

  USE pseudo_types, ONLY : pseudo_upf
  USE splinelib
  USE funct, ONLY : set_dft_from_name, get_iexch, get_icorr, get_igcx, get_igcc

  IMPLICIT NONE

  INTEGER :: i, j, ib, ijv, ijv2, prova
  real(8), INTENT(IN) :: x
  CHARACTER (len=256) :: filein(2)
  TYPE (pseudo_upf), INTENT(IN) :: upf(2)
  TYPE (pseudo_upf), INTENT(INOUT) :: upf_vca

  CHARACTER(len=9) :: day, hour
  CHARACTER (len=5) :: xlabel

  LOGICAL, EXTERNAL :: matches

  !
  ! All variables to be written into the UPF file
  ! (UPF = unified pseudopotential format, v.2)
  !

  ! pp_info
  INTEGER :: upf_rel
  real(8) :: upf_rcloc

  ! pp_header
  INTEGER :: upf_iexch, upf_icorr, upf_igcx, upf_igcc
  INTEGER :: upf_lmax, upf_mesh, upf_nbeta, upf_ntwfc
  real(8) :: upf_zp, upf_ecutrho, upf_ecutwfc, upf_etotps
  LOGICAL :: upf_nlcc
  real(8), ALLOCATABLE :: upf_ocw(:)
  CHARACTER(len=2), ALLOCATABLE :: upf_elsw(:)
  INTEGER, ALLOCATABLE :: upf_lchiw(:)

  ! pp_mesh
  real(8), ALLOCATABLE :: upf_r(:), upf_rab(:)
  real(8), ALLOCATABLE :: upf_rho_atc(:)
  real(8) :: capel
  real(8), ALLOCATABLE :: aux1(:,:), aux2(:,:)
  LOGICAL :: interpolate

  ! pp_local
  real(8), ALLOCATABLE :: upf_vloc0(:)

  ! pp_nonlocal
  ! pp_beta
  real(8), ALLOCATABLE :: upf_betar(:,:)
  INTEGER, ALLOCATABLE :: upf_lll(:), upf_ikk2(:)
  ! pp_dij
  real(8), ALLOCATABLE :: upf_dion(:,:)
  ! pp_qij
  INTEGER :: l, l1, l2
  INTEGER :: upf_nqf, upf_nqlc
  real(8), ALLOCATABLE :: upf_rinner(:), upf_qqq(:,:), upf_qfunc(:,:), upf_qfuncl(:,:,:)
  ! pp_qfcoef
  real(8), ALLOCATABLE :: upf_qfcoef(:,:,:,:)
  !
  ! pp_pswfc
  real(8), ALLOCATABLE :: upf_chi(:,:)
  !
  ! pp_rhoatom
  real(8), ALLOCATABLE :: upf_rho_at(:)
  !
  ! pp_spin_orb
  real(8), ALLOCATABLE :: upf_jjj(:)

  interpolate = .false.

  upf_vca = upf(1)

  !pp_info
  upf_rel = -1
  upf_rcloc = 0.d0
  !

  !pp_header
  upf_vca%generated  = 'Generated using virtual.x code'
  upf_vca%author = 'virtual_v2'

  call date_and_tim(day, hour)
  upf_vca%date = trim(day)

  WRITE( xlabel, '(f5.3)' ) x
  upf_vca%comment    = 'Pseudo = x '//trim(filein(1))//&
                   ' + (1-x) '//trim(filein(2))//', with x='//xlabel
  upf_vca%psd = "Xx"
  upf_vca%typ = "NC"
  IF (matches(upf(1)%typ, "USPP") .or. matches(upf(2)%typ, "USPP")) &
     upf_vca%typ = "USPP"

  CALL set_dft_from_name(upf(1)%dft)
  upf_iexch = get_iexch()
  upf_icorr = get_icorr()
  upf_igcx  = get_igcx()
  upf_igcc  = get_igcc()
  CALL set_dft_from_name(upf(2)%dft)
  IF (get_iexch()/=upf_iexch .or. get_icorr()/=upf_icorr .or. &
      get_igcx()/=upf_igcx .or. get_igcc()/=upf_igcc) &
      CALL errore ('virtual','conflicting DFT functionals',1)

  upf_lmax = max(upf(1)%lmax, upf(2)%lmax)

  IF ( upf(1)%mesh/=upf(2)%mesh ) THEN
     WRITE (*,*) " pseudopotentials have different mesh "
     WRITE (*,*) upf(1)%mesh, upf(2)%mesh
     WRITE (*,*) upf(1)%r(1), upf(2)%r(1)
     WRITE (*,*) upf(1)%r(upf(1)%mesh), upf(2)%r(upf(2)%mesh)
     interpolate = .true.
  ENDIF

  upf_mesh = upf(1)%mesh
  upf_nbeta = upf(1)%nbeta+upf(2)%nbeta
  upf_ntwfc = upf(1)%nwfc
  upf_nlcc  = upf(1)%nlcc.or.upf(2)%nlcc
  upf_ecutrho = upf(1)%ecutrho
  upf_ecutwfc = upf(1)%ecutwfc
  upf_etotps  = upf(1)%etotps


  ALLOCATE( upf_ocw(upf_ntwfc), upf_elsw(upf_ntwfc), upf_lchiw(upf_ntwfc) )
  upf_ocw(1:upf_ntwfc)  = upf(1)%oc(1:upf_ntwfc)
  upf_elsw(1:upf_ntwfc) = upf(1)%els(1:upf_ntwfc)
  upf_lchiw(1:upf_ntwfc) = upf(1)%lchi(1:upf_ntwfc)
  upf_zp    =  x * upf(1)%zp + (1.d0-x) * upf(2)%zp
  !

  !pp_mesh
  capel = 0.d0
  DO i=1,upf_mesh
     IF (i<=upf(2)%mesh) THEN
        capel = capel + abs(upf(1)%r(i)-upf(2)%r(i)) + abs(upf(1)%rab(i)-upf(2)%rab(i))
     ELSE
        capel = capel + abs(upf(1)%r(i)) + abs(upf(1)%rab(i))
     ENDIF
  ENDDO
  IF (capel>1.d-6) THEN
     WRITE (*,*) " pseudopotentials have different mesh "
     WRITE (*,*) "capel = ", capel
     interpolate = .true.
  ENDIF

  WRITE (*,*) "INTERPOLATE = ", interpolate
  IF (interpolate) ALLOCATE ( aux1(1,upf(1)%mesh), aux2(1,upf(2)%mesh) )

  ALLOCATE( upf_r(upf_mesh), upf_rab(upf_mesh) )
  upf_r(1:upf_mesh)   = upf(1)%r(1:upf_mesh)
  upf_rab(1:upf_mesh) = upf(1)%rab(1:upf_mesh)
  !

  !pp_nlcc
  ALLOCATE( upf_rho_atc(upf_mesh) )
  IF (interpolate) THEN
     WRITE (*,*) "interpolate rho_atc"
     aux2(1,1:upf(2)%mesh) = upf(2)%rho_atc(1:upf(2)%mesh)
     CALL dosplineint( upf(2)%r(1:upf(2)%mesh), aux2, upf_r(1:upf_mesh), aux1 )
     ! upf(2)%rho_atc(1:upf_mesh) = aux1(1,1:upf_mesh)
     upf_rho_atc(1:upf_mesh) =    x     * upf(1)%rho_atc(1:upf_mesh) + &
                            (1.d0-x) * aux1(1,1:upf_mesh)
     WRITE (*,*) " done"
  ELSE
     upf_rho_atc(1:upf_mesh) =    x     * upf(1)%rho_atc(1:upf_mesh) + &
                            (1.d0-x) * upf(2)%rho_atc(1:upf_mesh)
  ENDIF
  !

  !pp_local
  ALLOCATE( upf_vloc0(upf_mesh) )
  IF (interpolate) THEN
     WRITE (*,*) " interpolate vloc0"
     aux2(1,1:upf(2)%mesh) =  upf(2)%vloc(1:upf(2)%mesh)

     CALL dosplineint( upf(2)%r(1:upf(2)%mesh), aux2, upf_r(1:upf_mesh), aux1 )

     ! upf(2)%vloc(1:upf_mesh) = aux1(1,1:upf_mesh)

     ! Jivtesh - if the mesh of the first atom extends to a larger radius
     ! than the mesh of the second atom, then, for those radii that are
     ! greater than the maximum radius of the second atom, the local potential
     ! of the second atom is calculated using the expression
     ! v_local = (-2)*Z/r instead of using the extrapolated value.
     ! This is because, typically extrapolation leads to positive potentials.
     ! This is implemented in lines 240-242

     DO i=1,upf(1)%mesh
        IF ( upf(1)%r(i) > upf(2)%r(upf(2)%mesh) ) &
           ! upf(2)%vloc(i) = -(2.0*upf(2)%zp)/upf(1)%r(i)
           aux1(1,i) = -(2.0*upf(2)%zp)/upf(1)%r(i)
     ENDDO
     upf_vloc0(1:upf_mesh) =      x     * upf(1)%vloc(1:upf_mesh) + &
                               (1.d0-x) * aux1(1,1:upf_mesh)
  ELSE
     upf_vloc0(1:upf_mesh) =      x     * upf(1)%vloc(1:upf_mesh) + &
                               (1.d0-x) * upf(2)%vloc(1:upf_mesh)    
  ENDIF
  !

  !pp_nonlocal
  !pp_beta
  ALLOCATE( upf_betar(upf_mesh, upf_nbeta), &
            upf_lll(upf_nbeta), upf_ikk2(upf_nbeta) )
  ib = 0
  DO i=1,upf(1)%nbeta
     ib  = ib + 1
     upf_betar(1:upf_mesh,ib) = upf(1)%beta(1:upf_mesh,i)
     upf_lll(ib)              = upf(1)%lll(i)
     upf_ikk2(ib)             = upf(1)%kbeta(i)
  ENDDO

  DO i=1,upf(2)%nbeta
     ib  = ib + 1
     IF (interpolate) THEN
        WRITE (*,*) " interpolate betar"
        aux2(1,1:upf(2)%mesh) = upf(2)%beta(1:upf(2)%mesh,i)
        CALL dosplineint( upf(2)%r(1:upf(2)%mesh), aux2, upf_r(1:upf_mesh), aux1 )
        ! upf(2)%beta(1:upf_mesh,i) = aux1(1,1:upf_mesh)
        upf_betar(1:upf_mesh,ib) = aux1(1,1:upf_mesh)
     ELSE
        upf_betar(1:upf_mesh,ib) = upf(2)%beta(1:upf_mesh,i)        
     ENDIF

     upf_lll(ib)              = upf(2)%lll(i)

     ! SdG - when the meshes of the two pseudo are different the ikk2 limits
     ! for the beta functions of the second one must be set properly
     ! This is done in lines 273-277
     IF (interpolate) THEN
        j = 1
        DO WHILE ( upf_r(j) < upf(2)%r(upf(2)%kbeta(i)) )
           j = j + 1
        ENDDO
        upf_ikk2(ib) = j
     ELSE
        upf_ikk2(ib) = upf(2)%kbeta(i)
     ENDIF

  ENDDO
  !

  WRITE (*,*) "upf(1)%lll = ", upf(1)%lll
  WRITE (*,*) "upf(2)%lll = ", upf(2)%lll






  !pp_dij
  ALLOCATE( upf_dion(upf_nbeta, upf_nbeta) )
  upf_dion(:,:) = 0.d0

  DO i=1,upf(1)%nbeta
     DO j=1,upf(1)%nbeta
        upf_dion(i,j) = x * upf(1)%dion(i,j)
     ENDDO
  ENDDO

  DO i=1,upf(2)%nbeta
     DO j=1,upf(2)%nbeta
        upf_dion(upf(1)%nbeta+i, upf(1)%nbeta+j) = (1.d0-x) * upf(2)%dion(i,j)
     ENDDO
  ENDDO
  !

  WRITE (*,*) "pp_dij completed."


  IF (matches(upf_vca%typ, "USPP")) THEN

    !pp_qij
    IF (upf(1)%nqf/=upf(2)%nqf) &
        CALL errore ("Virtual","different nqf are not implemented (yet)", 1)
  
    IF (upf(1)%nqlc/=upf(2)%nqlc) &
        CALL errore ("Virtual","different nqlc are not implemented (yet)", 1)
  
    upf_nqf = upf(1)%nqf
    upf_nqlc = upf(1)%nqlc
  
    ALLOCATE( upf_rinner(upf_nqlc) )
    DO i=1,upf_nqlc
      IF(upf(1)%rinner(i)/=upf(2)%rinner(i)) &
         CALL errore("Virtual","different rinner are not implemented (yet)",i)
    ENDDO
  
    upf_rinner(1:upf_nqlc) = upf(1)%rinner(1:upf_nqlc)
  
    ALLOCATE( upf_qqq(upf_nbeta,upf_nbeta) )
    upf_qqq(:,:) = 0.d0
    IF ( upf(1)%q_with_l .NEQV. upf(2)%q_with_l) &
       CALL errore ( 'virtual.x: ', 'Augmentation charges are written using uncompatible formats',1) 
    IF( upf(1)%q_with_l ) THEN
       upf_vca%q_with_l = .TRUE. 
       ALLOCATE( upf_qfuncl(upf_mesh, upf_nbeta*(upf_nbeta+1)/2, 0:2*upf(1)%lmax) )
       upf_qfuncl(:,:,:) = 0.d0
    ELSE
       upf_vca%q_with_l = .FALSE. 
       ALLOCATE( upf_qfunc(upf_mesh, upf_nbeta*(upf_nbeta+1)/2) )
       upf_qfunc(:,:) = 0.d0
    ENDIF
  
    WRITE (*,*) "pp_qij"
  
  
    DO i=1,upf(1)%nbeta
       DO j=i,upf(1)%nbeta
          ijv = j * (j-1)/2 + i
          upf_qqq(i,j) = x * upf(1)%qqq(i,j)
  
  
          IF( ALLOCATED(upf_qfuncl) ) THEN
             l1=upf(1)%lll(i)
             l2=upf(1)%lll(j)
             DO l=abs(l1-l2), l1+l2
                upf_qfuncl(1:upf_mesh,ijv,l) = x * upf(1)%qfuncl(1:upf_mesh,ijv,l)
             ENDDO
          ELSE
             upf_qfunc(1:upf_mesh,ijv) = x * upf(1)%qfunc(1:upf_mesh,ijv)
          ENDIF
  
       ENDDO
    ENDDO
  
  
    DO i=1,upf(2)%nbeta
       DO j=i,upf(2)%nbeta
          ijv = j * (j-1)/2 + i
          !ijv  = (upf(1)%nbeta+j)*(upf(1)%nbeta+j-1)/2 + i + upf(1)%nbeta
          upf_qqq(upf(1)%nbeta+i,upf(1)%nbeta+j) = (1.d0-x) * upf(2)%qqq(i,j)
  
          ijv2 = (upf(1)%nbeta+j) * (upf(1)%nbeta+j-1) / 2 + (upf(1)%nbeta+i)
  
  
          IF( ALLOCATED(upf_qfuncl) ) THEN
             l1=upf(2)%lll(i)
             l2=upf(2)%lll(j)
             DO l=abs(l1-l2), l1+l2
                IF (interpolate) THEN
                   WRITE (*,*) " interpolate qfunc"
                   aux2(1,1:upf(2)%mesh) = upf(2)%qfuncl(1:upf(2)%mesh,ijv,l)
                   CALL dosplineint( upf(2)%r(1:upf(2)%mesh), aux2, upf_r(1:upf_mesh), aux1 )
                   upf_qfuncl(1:upf_mesh, ijv2, l) = (1.d0-x) * aux1(1,1:upf_mesh)
                   WRITE (*,*) " done"  
                ELSE
                   upf_qfuncl(1:upf_mesh, ijv2, l) = (1.d0-x) * upf(2)%qfuncl(1:upf_mesh,ijv,l)
  
                   IF ((i==1) .AND. (j==1)) THEN
                   ENDIF
  
                ENDIF
             ENDDO
          ELSE
             IF (interpolate) THEN
                 aux2(1,1:upf(2)%mesh) = upf(2)%qfunc(1:upf(2)%mesh,ijv)
                CALL dosplineint( upf(2)%r(1:upf(2)%mesh), aux2, upf_r(1:upf_mesh), aux1 )
                upf_qfunc(1:upf_mesh, ijv2) = (1.d0-x) * aux1(1,1:upf_mesh)
             ELSE
                upf_qfunc(1:upf_mesh, ijv2) = (1.d0-x) * upf(2)%qfunc(1:upf_mesh,ijv)
             ENDIF
          ENDIF
          ! ! upf_qfunc(1:upf_mesh, upf(1)%nbeta+i, upf(1)%nbeta+j) = (1.d0-x) * qfunc(1:upf_mesh,i,j,2)
       ENDDO
    ENDDO
    !
  
    !pp_qfcoef
    IF (upf_nqf/=0) THEN
       ALLOCATE( upf_qfcoef(upf_nqf,upf_nqlc,upf_nbeta,upf_nbeta) )
       upf_qfcoef(:,:,:,:) = 0.d0
       DO i=1,upf(1)%nbeta
          DO j=1,upf(1)%nbeta
             upf_qfcoef(1:upf_nqf,1:upf_nqlc,i,j) = &
                 x * upf(1)%qfcoef(1:upf_nqf,1:upf_nqlc,i,j)
          ENDDO
       ENDDO
       DO i=1,upf(2)%nbeta
          DO j=1,upf(2)%nbeta
             upf_qfcoef(1:upf_nqf,1:upf_nqlc,upf(1)%nbeta+i,upf(1)%nbeta+j) = &
                 (1.d0-x) * upf(2)%qfcoef(1:upf_nqf,1:upf_nqlc,i,j)
          ENDDO
       ENDDO
    ENDIF
    !
  ENDIF

  !pp_pswfc
  ALLOCATE ( upf_chi(upf_mesh,upf_ntwfc) )
  IF (upf(1)%nwfc==upf(2)%nwfc) THEN
     DO i=1,upf(2)%nwfc
        IF (interpolate) THEN
           WRITE (*,*) " interpolate chi"
           aux2(1,1:upf(2)%mesh) = upf(2)%chi(1:upf(2)%mesh,i)
           CALL dosplineint( upf(2)%r(1:upf(2)%mesh), aux2, upf_r(1:upf_mesh), aux1 )
           ! chi(1:upf_mesh,i,2) = aux1(1,1:upf_mesh)
           upf_chi(1:upf_mesh,i) =    x     * upf(1)%chi(1:upf_mesh,i) + &
              (1.d0-x) * aux1 (1,1:upf_mesh) 
        ELSE
        upf_chi(1:upf_mesh,i) =    x     * upf(1)%chi(1:upf_mesh,i) + &
                                (1.d0-x) * upf(2)%chi(1:upf_mesh,i)          
        ENDIF
        ! Jivtesh - The wavefunctions are calculated to be the average of the
        ! wavefunctions of the two atoms - lines 365-366
     ENDDO
  ELSE
     WRITE (*,*) "Number of wavefunctions not the same for the two pseudopotentials"
  ENDIF
  !upf_chi(1:upf_mesh,1:upf_ntwfc) = chi(1:upf_mesh,1:upf_ntwfc,1)
  !
  !pp_rhoatm
  ALLOCATE ( upf_rho_at(upf_mesh) )
  IF (interpolate) THEN
     WRITE (*,*) " interpolate rho_at"
     aux2(1,1:upf(2)%mesh) = upf(2)%rho_at(1:upf(2)%mesh)
     CALL dosplineint( upf(2)%r(1:upf(2)%mesh), aux2, upf_r(1:upf_mesh), aux1 )
     ! rho_at(1:upf_mesh,2) = aux1(1,1:upf_mesh)
     upf_rho_at(1:upf_mesh) =    x     * upf(1)%rho_at(1:upf_mesh) + &
                              (1.d0-x) * aux1(1,1:upf_mesh)
     WRITE (*,*) " done"
  ELSE
     upf_rho_at(1:upf_mesh) =    x     * upf(1)%rho_at(1:upf_mesh) + &
                              (1.d0-x) * upf(2)%rho_at(1:upf_mesh)    
  ENDIF

  !pp_spin_orb
  IF (upf_vca%has_so) THEN
    WRITE (*,*) " this pseudopotential has spin orbit coupling!"
    ALLOCATE (upf_jjj(upf_nbeta))
    upf_jjj(1:upf(1)%nbeta) = upf(1)%jjj
    upf_jjj(upf(1)%nbeta+1:upf_nbeta) = upf(2)%jjj
  ENDIF


  ! pp_info
  ! pp_header
  upf_vca%lmax = upf_lmax
  upf_vca%nlcc = upf_nlcc
  upf_vca%zp = upf_zp
  ! pp_nlcc
  upf_vca%rho_atc = upf_rho_atc
  ! pp_local
  upf_vca%vloc = upf_vloc0

  ! pp_pswfc
  upf_vca%chi = upf_chi
  ! pp_rhoatom
  upf_vca%rho_at = upf_rho_at

  ! pp_nonlocal
  ! pp_beta
  upf_vca%nbeta = upf_nbeta

  NULLIFY( upf_vca%kbeta,          &
           upf_vca%lll,            &
           upf_vca%beta,           &
           upf_vca%els_beta,       &
           upf_vca%dion,           &
           upf_vca%qqq,            &
           upf_vca%rcut,           &
           upf_vca%rcutus          &
  )

  ALLOCATE( upf_vca%kbeta(upf_nbeta),          &
            upf_vca%lll(upf_nbeta),            &
            upf_vca%beta(upf_mesh, upf_nbeta), &
            upf_vca%els_beta(upf_nbeta),       &
            upf_vca%dion(upf_nbeta, upf_nbeta),&
            upf_vca%qqq(upf_nbeta, upf_nbeta), &
            upf_vca%rcut(upf_nbeta),           &
            upf_vca%rcutus(upf_nbeta)          &
  )

  upf_vca%kbeta = upf_ikk2
  upf_vca%beta = upf_betar
  upf_vca%lll = upf_lll
  ! pp_dij
  upf_vca%dion = upf_dion

  IF (matches(upf_vca%typ, "USPP")) THEN
    ! pp_qij
    upf_vca%qqq = upf_qqq
    IF( upf_vca%q_with_l ) THEN
       NULLIFY( upf_vca%qfuncl )
       ALLOCATE( upf_vca%qfuncl(upf_mesh, upf_nbeta*(upf_nbeta+1)/2, 0:2*upf_lmax) )
       upf_vca%qfuncl = upf_qfuncl
    ELSE
       NULLIFY( upf_vca%qfunc )
       ALLOCATE( upf_vca%qfunc(upf_mesh, upf_nbeta*(upf_nbeta+1)/2) )
       upf_vca%qfunc = upf_qfunc
    ENDIF
    ! pp_qfcoef
    IF ( ALLOCATED(upf_qfcoef) ) THEN
       NULLIFY( upf_vca%qfcoef )
       ALLOCATE( upf_vca%qfcoef(upf_nqf, upf_nqlc, upf_nbeta, upf_nbeta) )
       upf_vca%qfcoef = upf_qfcoef
    ENDIF
  ENDIF

  upf_vca%kkbeta = maxval(upf_vca%kbeta)

  upf_vca%els_beta(1:upf(1)%nbeta) = upf(1)%els_beta
  upf_vca%els_beta(upf(1)%nbeta+1:upf_nbeta) = upf(2)%els_beta

  upf_vca%rcut(1:upf(1)%nbeta) = upf(1)%rcut
  upf_vca%rcut(upf(1)%nbeta+1:upf_nbeta) = upf(2)%rcut

  upf_vca%rcutus(1:upf(1)%nbeta) = upf(1)%rcutus
  upf_vca%rcutus(upf(1)%nbeta+1:upf_nbeta) = upf(2)%rcutus

  IF (upf_vca%has_so) THEN
     NULLIFY( upf_vca%jjj )
     ALLOCATE( upf_vca%jjj(upf_nbeta) )
     upf_vca%jjj = upf_jjj
  ENDIF


  ! !! DEBUG
  ! !! upf(1)
  ! WRITE (*,*) "upf(1)%nbeta = ", upf(1)%nbeta
  ! WRITE (*,*) "shape of upf(1)%kbeta = ", shape(upf(1)%kbeta)
  ! WRITE (*,*) "upf(1)%kbeta = ", upf(1)%kbeta
  ! WRITE (*,*) "upf(1)%kkbeta = ", upf(1)%kkbeta
  ! WRITE (*,*) "shape of upf(1)%lll = ", shape(upf(1)%lll)
  ! WRITE (*,*) "shape of upf(1)%beta = ", shape(upf(1)%beta)
  ! WRITE (*,*) "shape of upf(1)%els_beta = ", shape(upf(1)%els_beta)
  ! WRITE (*,*) "shape of upf(1)%dion = ", shape(upf(1)%dion)
  ! WRITE (*,*) "shape of upf(1)%qqq = ", shape(upf(1)%qqq)
  ! WRITE (*,*) "shape of upf(1)%qfuncl = ", shape(upf(1)%qfuncl)
  ! WRITE (*,*) "shape of upf(1)%qfcoef = ", shape(upf(1)%qfcoef)
  ! WRITE (*,*) "upf(1)%aewfc = ", upf(1)%aewfc
  ! WRITE (*,*) "upf(1)%pswfc = ", upf(1)%pswfc
  ! WRITE (*,*) "shape of upf(1)%rcut = ", shape(upf(1)%rcut)
  ! WRITE (*,*) "shape of upf(1)%rcutus = ", shape(upf(1)%rcutus)
  ! WRITE (*,*) ""
  ! !!! upf(2)
  ! WRITE (*,*) "upf(2)%nbeta = ", upf(2)%nbeta
  ! WRITE (*,*) "shape of upf(2)%kbeta = ", shape(upf(2)%kbeta)
  ! WRITE (*,*) "upf(2)%kbeta = ", upf(2)%kbeta
  ! WRITE (*,*) "upf(2)%kkbeta = ", upf(2)%kkbeta
  ! WRITE (*,*) "shape of upf(2)%lll = ", shape(upf(2)%lll)
  ! WRITE (*,*) "shape of upf(2)%beta = ", shape(upf(2)%beta)
  ! WRITE (*,*) "shape of upf(2)%els_beta = ", shape(upf(2)%els_beta)
  ! WRITE (*,*) "shape of upf(2)%dion = ", shape(upf(2)%dion)
  ! WRITE (*,*) "shape of upf(2)%qqq = ", shape(upf(2)%qqq)
  ! WRITE (*,*) "shape of upf(2)%qfuncl = ", shape(upf(2)%qfuncl)
  ! WRITE (*,*) "shape of upf(2)%qfcoef = ", shape(upf(2)%qfcoef)
  ! WRITE (*,*) "upf(2)%aewfc = ", upf(2)%aewfc
  ! WRITE (*,*) "upf(2)%pswfc = ", upf(2)%pswfc
  ! WRITE (*,*) "shape of upf(2)%rcut = ", shape(upf(2)%rcut)
  ! WRITE (*,*) "shape of upf(2)%rcutus = ", shape(upf(2)%rcutus)
  ! WRITE (*,*) ""
  ! !!! upf_vca
  ! WRITE (*,*) "upf_vca%nbeta = ", upf_vca%nbeta
  ! WRITE (*,*) "shape of upf_vca%kbeta = ", shape(upf_vca%kbeta)
  ! WRITE (*,*) "upf_vca%kbeta = ", upf_vca%kbeta
  ! WRITE (*,*) "upf_vca%kkbeta = ", upf_vca%kkbeta
  ! WRITE (*,*) "shape of upf_vca%lll = ", shape(upf_vca%lll)
  ! WRITE (*,*) "shape of upf_vca%beta = ", shape(upf_vca%beta)
  ! WRITE (*,*) "shape of upf_vca%els_beta = ", shape(upf_vca%els_beta)
  ! WRITE (*,*) "shape of upf_vca%dion = ", shape(upf_vca%dion)
  ! WRITE (*,*) "shape of upf_vca%qqq = ", shape(upf_vca%qqq)
  ! WRITE (*,*) "shape of upf_vca%qfuncl = ", shape(upf_vca%qfuncl)
  ! WRITE (*,*) "shape of upf_vca%qfcoef = ", shape(upf_vca%qfcoef)
  ! WRITE (*,*) "shape of upf_vca%rcut = ", shape(upf_vca%rcut)
  ! WRITE (*,*) "shape of upf_vca%rcutus = ", shape(upf_vca%rcutus)
  ! !! DEBUG


END SUBROUTINE compute_virtual
