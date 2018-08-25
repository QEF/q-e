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
!   V_{nl)^{(vca)} = \sum_{ij} |\beta^{(1)}_i>   x D^{(1)}_{ij} <\beta^{(1)}_j|
!                  + \sum_{ij} |\beta^{(2)}_i>(1-x)D^{(2)}_{ij} <\beta^{{2)}_j|
! where
!   V_{loc}^{(n)}(r) is the local part of pseudopot n
!   \beta^{{n)}_i(r) are the projectors for pseudopot n
!   D^{(n))_{ij} are the (bare) components of matrix D for pseudopot n
!
!
PROGRAM virtual
  !---------------------------------------------------------------------
  !
  !  Read pseudopotentials in the Unified Pseudopotential Format (UPF)
  !
  USE pseudo_types, ONLY: pseudo_upf, nullify_pseudo_upf
  USE upf_module,   ONLY: read_upf 
  USE write_upf_module, ONLY: write_upf
  USE environment, ONLY: environment_start, environment_end
  USE mp_global, ONLY: mp_startup, mp_global_end
  USE io_global, ONLY: ionode, stdout
  IMPLICIT NONE
  TYPE (pseudo_upf)     :: upf_input(2), upf_out
  INTEGER :: is, ios, iunps = 4
  real (8) :: x
  CHARACTER (len=256) :: filein(2), fileout
  INTERFACE 
     SUBROUTINE compute_virtual(x, upf_in1, upf_in2, upf_out, filein)
        import pseudo_upf 
        implicit none 
        real(8),intent(in)  :: x
        type(pseudo_upf),intent(in)         :: upf_in1, upf_in2 
        type(pseudo_upf),intent(out)        :: upf_out
        character(len=*),intent(in)         :: filein(2)
     END SUBROUTINE
   END INTERFACE
#if defined(__MPI)
  CALL mp_startup()
#endif
  CALL environment_start('VIRTUAL.X') 
  IF (ionode) THEN 
     CALL nullify_pseudo_upf(upf_input(1))
     CALL nullify_pseudo_upf(upf_input(2))
     CALL nullify_pseudo_upf(upf_out) 
     PRINT '(" ")'
     PRINT '(" Generate the UPF pseudopotential for a virtual atom ")'
     PRINT '(" combining two pseudopootentials in UPF format ")'
     PRINT '(" ")'
     !
     DO is=1,2
        WRITE(stdout,'("  Input PP file # ",i2," in UPF format > ")', advance="NO") is
        FLUSH(stdout) 
        READ (5, '(a)', end = 20, err = 20) filein(is)
        !OPEN(unit=iunps,file=filein(is),status='old',form='formatted',iostat=ios)
        CALL read_upf( upf_input(is), IERR = ios, FILENAME = TRIM(filein(is)))  
        IF (ios > 0) CALL errore('virtual.x:' , 'error parsing UPF FILE', ios) 
        WRITE (stdout,*) " IOS= ", ios, is, TRIM(filein(is)) 
        FLUSH (stdout) 
        PRINT '(" ")'
     ENDDO
     PRINT '(" New Pseudo = x ",a," + (1-x) ",a)', (trim(filein(is)), is=1,2)
10   CONTINUE
       WRITE(stdout,'(" mixing parameter x [0<x<1] = ")', advance="NO")
       FLUSH(stdout)
       READ (5,*) x
       IF (x<0.d0 .or. x>1)  GOTO 10

     CALL compute_virtual(x,upf_input(1), upf_input(2), upf_out, filein)

     fileout='NewPseudo.UPF'
     PRINT '("Output PP file in UPF format :  ",a)', fileout

     CALL write_upf ( TRIM(fileout), upf_out, SCHEMA='v2') 
  END IF  
  CALL environment_end('VIRTUAL.X')
#if defined(__MPI) 
  CALL mp_global_end()
#endif 
  STOP
20 CALL errore ('virtual.x', 'error reading pseudo file', 1)   
END PROGRAM virtual
!
!---------------------------------------------------------------------
SUBROUTINE compute_virtual(x_,upf_in1, upf_in2, upf_out,filein)
  USE kinds, ONLY: DP
  USE pseudo_types, ONLY: pseudo_upf , nullify_pseudo_upf 
  USE splinelib
  USE funct, ONLY : set_dft_from_name, get_iexch, get_icorr, get_igcx, get_igcc
  IMPLICIT NONE
  REAL(DP),INTENT(IN) :: x_
  TYPE(pseudo_upf),INTENT(IN) :: upf_in1, upf_in2
  CHARACTER (len=256),INTENT(IN) :: filein(2)
  TYPE(pseudo_upf),INTENT(OUT):: upf_out
  !
  TYPE (pseudo_upf)      :: upf_in(2)  
  INTEGER :: i, j, ib, il, ijv
  CHARACTER (len=5) :: xlabel
  real (DP) :: x, capel,   temp_iexch, temp_icorr, temp_igcc, temp_igcx
  real (DP), ALLOCATABLE :: aux1(:,:), aux2(:,:)
  LOGICAL :: interpolate
  interpolate = .false.
  call nullify_pseudo_upf(upf_out)
  x = x_
  !
  IF (upf_in1%mesh /= upf_in2%mesh ) THEN
     WRITE (*,*) " pseudopotentials have different mesh "
     WRITE (*,*) upf_in1%mesh , upf_in2%mesh
     WRITE (*,*) upf_in1%r(1), upf_in2%r(1)
     WRITE (*,*) upf_in1%r(upf_in1%mesh),upf_in2%r(upf_in2%mesh)
     interpolate = .true.
  ENDIF
  IF ( interpolate .AND. upf_in1%mesh .GT. upf_in2%mesh) THEN 
     upf_in(1) = upf_in1 
     upf_in(2) = upf_in2
  ELSE IF (interpolate) THEN
     upf_in(1) = upf_in2
     upf_in(2) = upf_in1 
     x = 1 -x 
  ELSE 
     upf_in(1) = upf_in1 
     upf_in(2) = upf_in2
  END IF

  !pp_info
  upf_out%nv = upf_in(1)%nv
  upf_out%rel = "Scalar"
  upf_out%rcloc = 0.d0
  !
  !pp_header
  upf_out%generated  = 'Generated using virtual.x code '
  upf_out%author=      'Author unknown. '//&
                   'Refer to original pseudopotential files'
  upf_out%date  = 'Date unknown '//&
                  'Refer to original pseudopotential files'
  WRITE( xlabel, '(f5.3)' ) x_
  !upf_out%comment    =  trim(filein(1)) // &
  !                      "+ " // trim(filein(2)) //". x parm =" // trim(xlabel) 
  upf_out%comment = ""
  upf_out%psd = "Xx"
  upf_out%typ  = "NC"
  IF (upf_in(1)%tvanp .OR. upf_in(2)%tvanp ) THEN 
     upf_out%typ = "US"
     upf_out%tvanp = .TRUE.
  END IF
  IF (upf_in(1)%typ == "PAW" .or. upf_in(2)%typ == "PAW") THEN 
     PRINT '("PAW mixing not implemented")'
     STOP
  ENDIF
  CALL set_dft_from_name(upf_in(1)%dft)
  upf_out%dft = upf_in(1)%dft 
  temp_iexch = get_iexch()
  temp_icorr = get_icorr()
  temp_igcx  = get_igcx()
  temp_igcc  = get_igcc()
  CALL set_dft_from_name(upf_in(2)%dft)
  IF (get_iexch()/=temp_iexch .or. get_icorr()/=temp_icorr .or. &
      get_igcx()/=temp_igcx .or. get_igcc()/=temp_igcc) &
      CALL errore ('virtual','conflicting DFT functionals',1)
  upf_out%lmax = max(upf_in(1)%lmax, upf_in(2)%lmax)
  upf_out%mesh = upf_in(1)%mesh
  upf_out%nbeta = upf_in(1)%nbeta + upf_in(2)%nbeta
  upf_out%nwfc = upf_in(1)%nwfc
  upf_out%nlcc  = upf_in(1)%nlcc .or. upf_in(2)%nlcc
  upf_out%ecutrho = max(upf_in(1)%ecutrho, upf_in(2)%ecutrho)
  upf_out%ecutwfc = max(upf_in(1)%ecutwfc, upf_in(2)%ecutwfc)
  upf_out%etotps  = x*upf_in(1)%etotps + (1.d0 -x )*upf_in(2)%etotps
  ALLOCATE( upf_out%oc(upf_out%nwfc), upf_out%els(upf_out%nwfc), upf_out%lchi(upf_out%nwfc) )
  upf_out%oc(1:upf_out%nwfc)  = upf_in(1)%oc(1:upf_out%nwfc)
  upf_out%els(1:upf_out%nwfc) = upf_in(1)%els(1:upf_out%nwfc)
  upf_out%lchi(1:upf_out%nwfc) = upf_in(1)%lchi(1:upf_out%nwfc)
  ALLOCATE(upf_out%nchi(upf_out%nwfc))                  
  upf_out%nchi(1:upf_out%nwfc)       = upf_in(1)%nchi(1:upf_out%nwfc)
  ALLOCATE(upf_out%rcut_chi(upf_out%nwfc), upf_out%rcutus_chi(upf_out%nwfc), upf_out%epseu(upf_out%nwfc))
  upf_out%rcut_chi = upf_in(1)%rcut_chi
  upf_out%rcutus_chi =upf_in(1)%rcutus_chi
  upf_out%epseu = upf_in(1)%epseu
  upf_out%zp    =  x * upf_in(1)%zp + (1.d0-x) * upf_in(2)%zp
  upf_out%lloc = upf_in(1)%lloc
  upf_out%has_so = .FALSE.
  upf_out%has_gipaw = .FALSE.
  !
  !pp_mesh
  capel = 0.d0
  DO i=1,upf_out%mesh
     capel = capel + abs(upf_in(1)%r(i)-upf_in(2)%r(i)) + abs(upf_in(1)%rab(i)-upf_in(2)%rab(i))
  ENDDO
  IF (capel>1.d-6) THEN
     WRITE (*,*) " pseudopotentials have different mesh "
     interpolate = .true.
  ENDIF
  WRITE (*,*) "INTERPOLATE =", interpolate
  !if (interpolate) call errore ("virtual", &
  !                "grid interpolation is not working yet",1)

  IF (interpolate) ALLOCATE ( aux1(1,upf_in(1)%mesh), aux2(1,upf_in(2)%mesh) )

  ALLOCATE( upf_out%r(upf_out%mesh), upf_out%rab(upf_out%mesh) )
  upf_out%r(1:upf_out%mesh)   = upf_in(1)%r(1:upf_out%mesh) 
  upf_out%rab(1:upf_out%mesh) = upf_in(1)%rab(1:upf_out%mesh)
  !
  !pp_nlcc
  ALLOCATE( upf_out%rho_atc(upf_out%mesh) )
  IF (interpolate) THEN
     WRITE (*,*) "interpolate rho_atc"
     aux2(1,1:upf_in(2)%mesh ) = upf_in(2)%rho_atc(1:upf_in(2)%mesh )
     CALL dosplineint( upf_in(2)%r(1:upf_in(2)%mesh ), aux2, upf_out%r(1:upf_out%mesh), aux1 )
     !rho_atc(1:upf_mesh,2) = aux1(1,1:upf_mesh):
     WRITE (*,*) " done"
     upf_out%rho_atc (1:upf_out%mesh) = x * upf_in(1)%rho_atc(1:upf_out%mesh) + &
                                       (1.d0-x) * aux1(1,1:upf_out%mesh) 
  ELSE
     upf_out%rho_atc(1:upf_out%mesh) =    x     * upf_in(1)%rho_atc(1:upf_out%mesh) + &
                              (1.d0-x) *   upf_in(2)%rho_atc(1:upf_out%mesh)
  END IF
  !
  !pp_local
  ALLOCATE( upf_out%vloc(upf_out%mesh) )
  IF (interpolate) THEN
     WRITE (*,*) " interpolate vloc0"
     aux2(1,1:upf_in(2)%mesh) =  upf_in(2)%vloc(1:upf_in(2)%mesh )

     CALL dosplineint( upf_in(2)%r(1:upf_in(2)%mesh ), aux2, upf_out%r(1:upf_out%mesh), aux1 )

     ! PD. vloc0(1:upf_mesh,2) = aux1(1,1:upf_mesh)

     ! Jivtesh - if the mesh of the first atom extends to a larger radius
     ! than the mesh of the second atom, then, for those radii that are
     ! greater than the maximum radius of the second atom, the local potential
     ! of the second atom is calculated using the expression
     ! v_local = (-2)*Z/r instead of using the extrapolated value.
     ! This is because, typically extrapolation leads to positive potentials.
     ! This is implemented in lines 240-242

     DO i=1,upf_in(1)%mesh
        IF(upf_in(1)%r(i)>upf_in(2)%r(upf_in(2)%mesh)) aux1(1,i) = -(2.0*upf_in(2)%zp)/upf_in(1)%r(i)
     ENDDO
     upf_out%vloc(1:upf_out%mesh) =  x  * upf_in(1)%vloc(1:upf_out%mesh) + &
                                   (1.d0 - x) * aux1(1,1:upf_out%mesh)
  ELSE
     upf_out%vloc(1:upf_out%mesh) =      x     * upf_in(1)%vloc(1:upf_out%mesh) +  &
                            (1.d0-x) * upf_in(2)%vloc(1:upf_out%mesh)
  END IF
  !
  !pp_nonlocal
  !pp_beta
  ALLOCATE( upf_out%beta(upf_out%mesh,upf_out%nbeta), &
            upf_out%lll(upf_out%nbeta), upf_out%kbeta(upf_out%nbeta) , upf_out%els_beta(upf_out%nbeta),& 
            upf_out%rcut(upf_out%nbeta), upf_out%rcutus(upf_out%nbeta)) 
  ib = 0
  DO i=1,upf_in(1)%nbeta
     ib  = ib + 1
     upf_out%beta(1:upf_out%mesh,ib) = upf_in(1)%beta(1:upf_out%mesh,i)
     upf_out%lll(ib)              = upf_in(1)%lll(i)
     upf_out%kbeta(ib)             = upf_in(1)%kbeta(i)
     upf_out%els_beta(ib)          = upf_in(1)%els_beta(i)
     upf_out%rcut(ib)              = upf_in(1)%rcut(i)
     upf_out%rcutus(ib)            = upf_in(1)%rcutus(i)
  ENDDO
  DO i=1,upf_in(2)%nbeta
     ib  = ib + 1
     IF (interpolate) THEN
     WRITE (*,*) " interpolate betar"
        aux2(1,1:upf_in(2)%mesh ) = upf_in(2)%beta(1:upf_in(2)%mesh,i)
        CALL dosplineint( upf_in(2)%r(1:upf_in(2)%mesh), aux2, upf_out%r(1:upf_out%mesh), aux1 )
        ! betar(1:upf_mesh,i,2) = aux1(1,1:upf_mesh)
        upf_out%beta(1:upf_out%mesh,ib) = aux1(1,1:upf_out%mesh) 
     ELSE 
        upf_out%beta(1:upf_out%mesh,ib) = upf_in(2)%beta(1:upf_out%mesh,i)
     ENDIF
     upf_out%lll(ib)              = upf_in(2)%lll(i)
     ! SdG - when the meshes of the two pseudo are different the ikk2 limits
     ! for the beta functions of the second one must be set properly
     ! This is done in lines 273-277
     IF (interpolate) THEN
        j = 1
        DO WHILE ( upf_out%r(j) < upf_in(2)%r( upf_in(2)%kbeta(i)) )
           j = j + 1
        ENDDO
        upf_out%kbeta(ib)             = j
     ELSE
        upf_out%kbeta(ib)             = upf_in(2)%kbeta(i)
     ENDIF
     upf_out%els_beta(ib)             = upf_in(2)%els_beta(i)
     upf_out%rcut(ib)                 = upf_in(2)%rcut(i)
     upf_out%rcutus(ib)               = upf_in(2)%rcut(i)
  ENDDO
  !
  !pp_dij
  ALLOCATE( upf_out%dion(upf_out%nbeta, upf_out%nbeta) )
  upf_out%dion(:,:) = 0.d0
  IF (ASSOCIATED(upf_in(1)%dion)) THEN
     DO i=1,upf_in(1)%nbeta
         DO j=1,upf_in(1)%nbeta
            upf_out%dion(i,j) = x * upf_in(1)%dion(i,j)
         ENDDO
     ENDDO
  END IF 
  IF (ASSOCIATED(upf_in(2)%dion)) THEN
      DO i=1,upf_in(2)%nbeta
         DO j=1,upf_in(2)%nbeta
            upf_out%dion(upf_in(1)%nbeta+i,upf_in(1)%nbeta+j) = (1.d0-x) * upf_in(2)%dion(i,j)
         ENDDO
      ENDDO
  END IF
  !
  !pp_qij
   IF ( upf_out%typ == "US") THEN 
      IF (upf_in(1)%nqf/=upf_in(2)%nqf ) &
         CALL errore ("Virtual","different nqf are not implemented (yet)", 1)
      IF (upf_in(1)%nqlc/=upf_in(2)%nqlc ) &
         CALL errore ("Virtual","different nqlc are not implemented (yet)", 1)
      upf_out%nqf = max(upf_in(1)%nqf,upf_in(2)%nqf)
      upf_out%nqlc = upf_in(1)%nqlc
      IF (  upf_in(1)%q_with_l .NEQV. upf_in(2)%q_with_l) &
                  CALL errore ( "Virtual", "Q_ij in the two files are not compatible for mixing",1)  
      upf_out%q_with_l = upf_in(1)%q_with_l 
      ALLOCATE( upf_out%rinner(upf_out%nqlc), upf_out%qqq(upf_out%nbeta,upf_out%nbeta)) 
      IF (upf_out%q_with_l) THEN
         ALLOCATE (upf_out%qfuncl(upf_out%mesh,upf_out%nbeta*(upf_out%nbeta+1)/2,0:2*upf_out%lmax), &
                   upf_out%qfunc(1,1))
      ELSE 
         ALLOCATE (upf_out%qfunc(upf_out%mesh,upf_out%nbeta*(upf_out%nbeta+1)/2),&
                  upf_out%qfuncl(1,1,1))  
      END IF 
      DO i=1,upf_out%nqlc
         IF(upf_in(1)%rinner(i)/=upf_in(2)%rinner(i)) &
            CALL errore("Virtual","different rinner are not implemented (yet)",i)
      ENDDO
      upf_out%rinner(1:upf_out%nqlc) = upf_in(1)%rinner(1:upf_out%nqlc)

      upf_out%qqq(:,:) = 0.d0
      upf_out%qfunc(:,:) = 0.d0
      upf_out%qfuncl(:,:,:) = 0.d0
      DO i=1, upf_in(1)%nbeta
         DO j= 1,upf_in(1)%nbeta
            print *, upf_in(1)%nbeta,   size(upf_in(1)%qfunc,2)
            upf_out%qqq(i,j) = x * upf_in(1)%qqq(i, j)
            IF ( j .GE. i ) THEN 
               IF ( .NOT. upf_out%q_with_l) THEN 
                  upf_out%qfunc(1:upf_out%mesh, j*(j-1)/2+i ) = &
                                             x * upf_in(1)%qfunc(1:upf_out%mesh, j*(j-1)/2+i)
               ELSE
                  DO il = 0, 2*upf_in(1)%lmax
                     upf_out%qfuncl(1:upf_out%mesh,j*(j-1)/2+i,il) = &
                        x * upf_in(1)%qfuncl(1:upf_out%mesh, j*(j-1)/2+i,il)
                  END DO 
               END IF
            END IF
         ENDDO
      ENDDO
      DO i=1,upf_in(2)%nbeta
         DO j=1,upf_in(2)%nbeta
            ijv  = (upf_in(1)%nbeta+j)*(upf_in(1)%nbeta+j-1)/2 + i + upf_in(1)%nbeta
            upf_out%qqq(upf_in(1)%nbeta +i, upf_in(1)%nbeta +j) = (1.d0-x) * upf_in(2)%qqq(i, j)
            IF (interpolate .AND. (j .GE. i)  ) THEN
               WRITE (*,*) " interpolate qfunc"
               IF (.NOT. upf_out%q_with_l) THEN 
                  aux2(1,1:upf_in(2)%mesh ) = upf_in(2)%qfunc(1:upf_in(2)%mesh,j*(j-1)/2+i )
                  CALL dosplineint( upf_in(2)%r(1:upf_in(2)%mesh), aux2, upf_out%r(1:upf_out%mesh), aux1 )
                  WRITE (*,*) " done" 
                  upf_out%qfunc(1:upf_out%mesh, ijv) = (1.d0-x)*aux1(1,1:upf_out%mesh)
               ELSE 
                  DO il = 0,2*upf_in(2)%lmax
                     PRINT '("Interpolate qfuncl for l = ",I5)',il 
                     aux2(1,1:upf_in(2)%mesh ) = upf_in(2)%qfuncl(1:upf_in(2)%mesh,j*(j-1)/2+i,il)
                     CALL dosplineint( upf_in(2)%r(1:upf_in(2)%mesh), aux2, upf_out%r(1:upf_out%mesh), aux1 )
                     PRINT '("Done")' 
                     upf_out%qfuncl(1:upf_out%mesh, ijv,il) = (1.d0-x)*aux1(1,1:upf_out%mesh)
                  END DO
               END IF
            ELSE IF (j .GE. i ) THEN 
               IF (.NOT. upf_out%q_with_l ) THEN 
                  upf_out%qfunc(1:upf_out%mesh, ijv) = (1.d0-x) * upf_in(2)%qfunc(1:upf_out%mesh,ijv)
               ELSE 
                  DO il = 0, 2*upf_in(2)%lmax 
                     upf_out%qfuncl(1:upf_out%mesh, ijv,il) = (1.d0-x) * upf_in(2)%qfuncl(1:upf_out%mesh, ijv, il)
                  END DO
               END IF
            END IF
         ENDDO
      ENDDO
      !
      !pp_qfcoef
      ALLOCATE( upf_out%qfcoef(upf_out%nqf,upf_out%nqlc,upf_out%nbeta,upf_out%nbeta) )
      upf_out%qfcoef(:,:,:,:) = 0.d0
      DO i=1, upf_in(1)%nbeta
         DO j=1, upf_in(1)%nbeta
            upf_out%qfcoef(1:upf_out%nqf,1:upf_out%nqlc,i,j) = &
               x * upf_in(1)%qfcoef(1:upf_out%nqf,1:upf_out%nqlc,i,j)
         ENDDO
      ENDDO
      DO i=1, upf_in(2)%nbeta
         DO j=1, upf_in(2)%nbeta
            upf_out%qfcoef(1:upf_out%nqf,1:upf_out%nqlc,upf_in(1)%nbeta+i, upf_in(1)%nbeta +j ) = &
                                       (1.d0-x) * upf_in(2)%qfcoef(1:upf_out%nqf,1:upf_out%nqlc,i,j )
         ENDDO
      ENDDO
      !
   END IF
   !pp_pswfc

   ALLOCATE (upf_out%chi(upf_out%mesh,upf_out%nwfc) )

   IF (upf_in(1)%nwfc  == upf_in(2)%nwfc ) THEN
      DO i=1, upf_in(2)%nwfc
         IF (interpolate) THEN
            WRITE (*,*) " interpolate chi"
            aux2(1,1:upf_in(2)%mesh ) = upf_in(2)%chi(1:upf_in(2)%mesh , i)
            CALL dosplineint( upf_in(2)%r(1:upf_in(2)%mesh ), aux2, upf_out%r(1:upf_out%mesh), aux1 )
            WRITE (*,*) " done"
            upf_out%chi(1:upf_out%mesh,i) = x * upf_in(1)%chi(1:upf_out%mesh,i) + & 
                                          (1.d0 -x ) * aux1(1,1:upf_out%mesh) 
         ELSE 
            ! Jivtesh - The wavefunctions are calcuated to be the average of the
            ! wavefunctions of the two atoms - lines 330-331
            upf_out%chi(1:upf_out%mesh,i) =    x     * upf_in(1)%chi(1:upf_out%mesh,i) + &
                                (1.d0-x) * upf_in(2)%chi(1:upf_out%mesh,i)

         ENDIF
      ENDDO
   ELSE
      WRITE (*,*) "Number of wavefunctions not the same for the two pseudopotentials"
      upf_out%nwfc=0
   ENDIF
   !upf_chi(1:upf_mesh,1:upf_ntwfc) = chi(1:upf_mesh,1:upf_ntwfc,1)
   !
   !pp_rhoatm

   ALLOCATE (upf_out%rho_at(upf_out%mesh) )
   IF (interpolate) THEN
      WRITE (*,*) " interpolate rho_at"
      aux2(1,1:upf_in(2)%mesh ) = upf_in(2)%rho_at(1:upf_in(2)%mesh )
      CALL dosplineint( upf_in(2)%r(1:upf_in(2)%mesh), aux2, upf_out%r(1:upf_out%mesh), aux1 )
      WRITE (*,*) " done"
      upf_out%rho_at(1:upf_out%mesh) = x   * upf_in(1)%rho_at(1:upf_out%mesh) + &
                                    (1.d0 - x) * aux1(1,1:upf_out%mesh) 
   ELSE 
      upf_out%rho_at(1:upf_out%mesh) =    x     * upf_in(1)%rho_at(1:upf_out%mesh) + &
                           (1.d0-x) * upf_in(2)%rho_at(1:upf_out%mesh)

   END IF
END SUBROUTINE compute_virtual
