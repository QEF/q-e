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
MODULE pseudo
  !
  ! All variables to be read from the UPF file
  ! (UPF = unified pseudopotential format)
  !
  INTEGER ,PARAMETER :: npsx = 2
  ! npsx  : maximum number of different pseudopotentials
  INTEGER, PARAMETER :: lmaxx  = 3, nchix  = 6, ndm = 2000
  ! lmaxx : maximum non local angular momentum in PP
  ! nchix : maximum number of atomic wavefunctions per PP
  ! ndm   : maximum number of points in the radial mesh
  INTEGER, PARAMETER :: nbrx = 8, lqmax = 5, nqfx = 8
  ! nbrx  : maximum number of beta functions
  ! lqmax : maximum number of angular momentum of Q
  ! nqfx  : maximum number of coefficients in Q smoothing
  !
  ! pp_header
  CHARACTER (len=80):: generated, date_author, comment
  CHARACTER (len=2) :: psd(npsx), pseudotype
  CHARACTER (len=20):: dft(npsx)
  INTEGER :: lmax(npsx), mesh(npsx), nbeta(npsx), ntwfc(npsx)
  LOGICAL :: nlcc(npsx), isus(npsx)
  real(8) :: zp(npsx), ecutrho, ecutwfc, etotps
  real(8) :: oc(nchix,npsx)
  CHARACTER(len=2) :: els(nchix,npsx)
  INTEGER :: lchi(nchix,npsx)
  !
  ! pp_mesh
  real(8) :: r(ndm,npsx), rab(ndm,npsx)
  !   pp_nlcc
  real(8) :: rho_atc(ndm,npsx)
  !
  ! pp_local
  real(8) ::  vloc0(ndm,npsx)
  !
  ! pp_nonlocal
  ! pp_beta
  real(8) :: betar(ndm, nbrx, npsx)
  INTEGER :: lll(nbrx,npsx), ikk2(nbrx,npsx)
  ! pp_dij
  real(8) :: dion(nbrx,nbrx,npsx)
  ! pp_qij
  INTEGER ::  nqf(npsx), nqlc(npsx)
  real(8) :: rinner(lqmax,npsx), qqq(nbrx,nbrx,npsx), &
       qfunc(ndm,nbrx,nbrx,npsx)
  ! pp_qfcoef
  real(8) :: qfcoef(nqfx,lqmax,nbrx,nbrx,npsx)
  !
  ! pp_pswfc
  real(8) :: chi(ndm,nchix,npsx)
  !
  ! pp_rhoatom
  real(8) :: rho_at(ndm,npsx)
END MODULE pseudo
!
PROGRAM virtual
  !---------------------------------------------------------------------
  !
  !  Read pseudopotentials in the Unified Pseudopotential Format (UPF)
  !
  USE pseudo_mod, ONLY : read_pseudo
  IMPLICIT NONE
  INTEGER :: is, ios, iunps = 4
  real (8) :: x
  CHARACTER (len=256) :: filein(2), fileout
  PRINT '(" ")'
  PRINT '(" Generate the UPF pseudopotential for a virtual atom ")'
  PRINT '(" combining two pseudopootentials in UPF format ")'
  PRINT '(" ")'
  !
  DO is=1,2
     WRITE(*,'("  Input PP file # ",i2," in UPF format > ")', advance="NO") is
     READ (5, '(a)', end = 20, err = 20) filein(is)
     OPEN(unit=iunps,file=filein(is),status='old',form='formatted',iostat=ios)
     IF (ios/=0) STOP
     WRITE (*,*) " IOS= ", ios, is, iunps
     CALL read_pseudo(is, iunps)
     CLOSE (unit=iunps)
     PRINT '(" ")'
  ENDDO
  PRINT '(" New Pseudo = x ",a," + (1-x) ",a)', (trim(filein(is)), is=1,2)
10 CONTINUE
  WRITE(*,'(" mixing parameter x [0<x<1] = ")', advance="NO")
  READ (5,*) x
  IF (x<0.d0 .or. x>1)  GOTO 10

  CALL compute_virtual(x,filein)

  fileout='NewPseudo.UPF'
  PRINT '("Output PP file in UPF format :  ",a)', fileout

  OPEN(unit=2,file=fileout,status='unknown',form='formatted')
  CALL write_upf_v1(2)
  CLOSE (unit=2)

20 STOP
END PROGRAM virtual
!
!---------------------------------------------------------------------
SUBROUTINE compute_virtual(x,filein)
  USE pseudo
  USE upf, ONLY : &
           upf_rel => rel, upf_rcloc => rcloc, upf_nwfs => nwfs, &
           upf_oc => oc, upf_rcut => rcut, upf_rcutus => rcutus, &
           upf_epseu => epseu, upf_els => els, &
           upf_lchi => lchi, upf_nns => nns, &
           upf_generated => generated, upf_date_author => date_author, &
           upf_comment => comment, &
           upf_psd => psd, upf_pseudotype => pseudotype, &
           upf_iexch => iexch, &
           upf_icorr => icorr, &
           upf_igcx  => igcx, &
           upf_igcc => igcc, &
           upf_lmax => lmax, upf_mesh => mesh, &
           upf_nbeta => nbeta, upf_ntwfc => ntwfc, upf_nlcc => nlcc, &
           upf_zp => zp, upf_ecutrho => ecutrho, upf_ecutwfc => ecutwfc, &
           upf_etotps => etotps, upf_ocw => ocw, &
           upf_elsw => elsw, upf_lchiw =>lchiw, &
           upf_r => r, upf_rab => rab, &
           upf_rho_atc => rho_atc, &
           upf_vloc0   => vloc0, &
           upf_betar => betar, upf_lll => lll,  upf_ikk2 => ikk2, &
           upf_dion => dion, &
           upf_nqf => nqf, upf_nqlc => nqlc, &
           upf_rinner => rinner, upf_qqq => qqq, upf_qfunc => qfunc, &
           upf_qfcoef => qfcoef, &
           upf_chi => chi, &
           upf_rho_at  => rho_at
  USE splinelib
  USE funct, ONLY : set_dft_from_name, get_iexch, get_icorr, get_igcx, get_igcc
  IMPLICIT NONE
  INTEGER :: i, j, ib
  CHARACTER (len=256) :: filein(2)
  CHARACTER (len=5) :: xlabel
  real (8) :: x, capel
  real (8), ALLOCATABLE :: aux1(:,:), aux2(:,:)
  LOGICAL :: interpolate
  interpolate = .false.
  !
  !pp_info
  upf_rel = -1
  upf_rcloc = 0.d0
  !
  !pp_header
  upf_generated  = 'Generated using virtual.x code '
  upf_date_author= 'Author and generation date: unknown. '//&
                   'Refer to original pseudopotential files'
  WRITE( xlabel, '(f5.3)' ) x
  upf_comment    = 'Pseudo = x '//trim(filein(1))//&
                   ' + (1-x) '//trim(filein(2))//', with x='//xlabel
  upf_psd = "Xx"
  upf_pseudotype = "NC"
  IF (isus(1) .or. isus(2)) upf_pseudotype = "US"
  CALL set_dft_from_name(dft(1))
  upf_iexch = get_iexch()
  upf_icorr = get_icorr()
  upf_igcx  = get_igcx()
  upf_igcc  = get_igcc()
  CALL set_dft_from_name(dft(2))
  IF (get_iexch()/=upf_iexch .or. get_icorr()/=upf_icorr .or. &
      get_igcx()/=upf_igcx .or. get_igcc()/=upf_igcc) &
      CALL errore ('virtual','conflicting DFT functionals',1)
  upf_lmax = max(lmax(1), lmax(2))
  IF (mesh(1)/=mesh(2) ) THEN
     WRITE (*,*) " pseudopotentials have different mesh "
     WRITE (*,*) mesh(1),mesh(2)
     WRITE (*,*) r(1,1), r(1,2)
     WRITE (*,*) r(mesh(1),1),r(mesh(2),2)
     interpolate = .true.
  ENDIF
  upf_mesh = mesh(1)
  upf_nbeta = nbeta(1)+nbeta(2)
  upf_ntwfc = ntwfc(1)
  upf_nlcc  = nlcc(1).or.nlcc(2)
  upf_ecutrho = ecutrho
  upf_ecutwfc = ecutwfc
  upf_etotps  = etotps
  ALLOCATE( upf_ocw(upf_ntwfc), upf_elsw(upf_ntwfc), upf_lchiw(upf_ntwfc) )
  upf_ocw(1:upf_ntwfc)  = oc(1:upf_ntwfc,1)
  upf_elsw(1:upf_ntwfc) = els(1:upf_ntwfc,1)
  upf_lchiw(1:upf_ntwfc) = lchi(1:upf_ntwfc,1)
  upf_zp    =  x * zp(1) + (1.d0-x) * zp(2)
  !
  !pp_mesh
  capel = 0.d0
  DO i=1,upf_mesh
     capel = capel + abs(r(i,1)-r(i,2)) + abs(rab(i,1)-rab(i,2))
  ENDDO
  IF (capel>1.d-6) THEN
     WRITE (*,*) " pseudopotentials have different mesh "
     interpolate = .true.
  ENDIF
  WRITE (*,*) "INTERPOLATE =", interpolate
  !if (interpolate) call errore ("virtual", &
  !                "grid interpolation is not working yet",1)

  IF (interpolate) ALLOCATE ( aux1(1,mesh(1)), aux2(1,mesh(2)) )

  ALLOCATE( upf_r(upf_mesh), upf_rab(upf_mesh) )
  upf_r(1:upf_mesh)   = r(1:upf_mesh,1)
  upf_rab(1:upf_mesh) = rab(1:upf_mesh,1)
  !
  !pp_nlcc
  ALLOCATE( upf_rho_atc(upf_mesh) )
  IF (interpolate) THEN
     WRITE (*,*) "interpolate rho_atc"
     aux2(1,1:mesh(2)) = rho_atc(1:mesh(2),2)
     CALL dosplineint( r(1:mesh(2),2), aux2, upf_r(1:upf_mesh), aux1 )
     rho_atc(1:upf_mesh,2) = aux1(1,1:upf_mesh)
     WRITE (*,*) " done"
  ENDIF
  upf_rho_atc(1:upf_mesh) =    x     * rho_atc(1:upf_mesh,1) + &
                            (1.d0-x) * rho_atc(1:upf_mesh,2)
  !
  !pp_local
  ALLOCATE( upf_vloc0(upf_mesh) )
  IF (interpolate) THEN
     WRITE (*,*) " interpolate vloc0"
     aux2(1,1:mesh(2)) =  vloc0(1:mesh(2),2)

     CALL dosplineint( r(1:mesh(2),2), aux2, upf_r(1:upf_mesh), aux1 )

     vloc0(1:upf_mesh,2) = aux1(1,1:upf_mesh)

     ! Jivtesh - if the mesh of the first atom extends to a larger radius
     ! than the mesh of the second atom, then, for those radii that are
     ! greater than the maximum radius of the second atom, the local potential
     ! of the second atom is calculated using the expression
     ! v_local = (-2)*Z/r instead of using the extrapolated value.
     ! This is because, typically extrapolation leads to positive potentials.
     ! This is implemented in lines 240-242

     DO i=1,mesh(1)
        IF(r(i,1)>r(mesh(2),2)) vloc0(i,2) = -(2.0*zp(2))/r(i,1)
     ENDDO

  ENDIF
  upf_vloc0(1:upf_mesh) =      x     * vloc0(1:upf_mesh,1) +  &
                            (1.d0-x) * vloc0(1:upf_mesh,2)
  !
  !pp_nonlocal
  !pp_beta
  ALLOCATE( upf_betar(upf_mesh,upf_nbeta), &
            upf_lll(upf_nbeta), upf_ikk2(upf_nbeta) )
  ib = 0
  DO i=1,nbeta(1)
     ib  = ib + 1
     upf_betar(1:upf_mesh,ib) = betar(1:upf_mesh,i,1)
     upf_lll(ib)              = lll(i,1)
     upf_ikk2(ib)             = ikk2(i,1)
  ENDDO
  DO i=1,nbeta(2)
     ib  = ib + 1
     IF (interpolate) THEN
     WRITE (*,*) " interpolate betar"
        aux2(1,1:mesh(2)) = betar(1:mesh(2),i,2)
        CALL dosplineint( r(1:mesh(2),2), aux2, upf_r(1:upf_mesh), aux1 )
        betar(1:upf_mesh,i,2) = aux1(1,1:upf_mesh)
     ENDIF
     upf_betar(1:upf_mesh,ib) = betar(1:upf_mesh,i,2)
     upf_lll(ib)              = lll(i,2)
     ! SdG - when the meshes of the two pseudo are different the ikk2 limits
     ! for the beta functions of the second one must be set properly
     ! This is done in lines 273-277
     IF (interpolate) THEN
        j = 1
        DO WHILE ( upf_r(j) < r( ikk2(i,2), 2) )
           j = j + 1
        ENDDO
        upf_ikk2(ib)             = j
     ELSE
        upf_ikk2(ib)             = ikk2(i,2)
     ENDIF
  ENDDO
  !
  !pp_dij
  ALLOCATE( upf_dion(upf_nbeta, upf_nbeta) )
  upf_dion(:,:) = 0.d0
  DO i=1,nbeta(1)
     DO j=1,nbeta(1)
        upf_dion(i,j) = x * dion(i,j,1)
     ENDDO
  ENDDO
  DO i=1,nbeta(2)
     DO j=1,nbeta(2)
        upf_dion(nbeta(1)+i,nbeta(1)+j) = (1.d0-x) * dion(i,j,2)
     ENDDO
  ENDDO
  !
  !pp_qij
  IF (nqf(1)/=nqf(2)) &
      CALL errore ("Virtual","different nqf are not implemented (yet)", 1)
  IF (nqlc(1)/=nqlc(2)) &
      CALL errore ("Virtual","different nqlc are not implemented (yet)", 1)
  upf_nqf = nqf(1)
  upf_nqlc = nqlc(1)
  ALLOCATE( upf_rinner(upf_nqlc), upf_qqq(upf_nbeta,upf_nbeta), &
            upf_qfunc(upf_mesh,upf_nbeta,upf_nbeta) )
  DO i=1,upf_nqlc
    IF(rinner(i,1)/=rinner(i,2)) &
       CALL errore("Virtual","different rinner are not implemented (yet)",i)
  ENDDO
  upf_rinner(1:upf_nqlc) = rinner(1:upf_nqlc,1)

  upf_qqq(:,:) = 0.d0
  upf_qfunc(:,:,:) = 0.d0
  DO i=1,nbeta(1)
     DO j=1,nbeta(1)
        upf_qqq(i,j) = x * qqq(i, j,1)
        upf_qfunc(1:upf_mesh,i,j) = x * qfunc(1:upf_mesh,i,j,1)
     ENDDO
  ENDDO
  DO i=1,nbeta(2)
     DO j=1,nbeta(2)
        upf_qqq(nbeta(1)+i,nbeta(1)+j) = (1.d0-x) * qqq(i, j, 2)
        IF (interpolate) THEN
     WRITE (*,*) " interpolate qfunc"
           aux2(1,1:mesh(2) ) = qfunc(1:mesh(2),i,j,2)
           CALL dosplineint( r(1:mesh(2),2), aux2, upf_r(1:upf_mesh), aux1 )
           qfunc(1:upf_mesh,i,j,2) = aux1(1,1:upf_mesh)
     WRITE (*,*) " done"
        ENDIF
        upf_qfunc(1:upf_mesh,nbeta(1)+i,nbeta(1)+j) = (1.d0-x) * qfunc(1:upf_mesh,i,j,2)
     ENDDO
  ENDDO
  !
  !pp_qfcoef
  ALLOCATE( upf_qfcoef(upf_nqf,upf_nqlc,upf_nbeta,upf_nbeta) )
  upf_qfcoef(:,:,:,:) = 0.d0
  DO i=1,nbeta(1)
     DO j=1,nbeta(1)
        upf_qfcoef(1:upf_nqf,1:upf_nqlc,i,j) = &
            x * qfcoef(1:upf_nqf,1:upf_nqlc,i,j, 1)
     ENDDO
  ENDDO
  DO i=1,nbeta(2)
     DO j=1,nbeta(2)
        upf_qfcoef(1:upf_nqf,1:upf_nqlc,nbeta(1)+i,nbeta(1)+j) = &
            (1.d0-x) * qfcoef(1:upf_nqf,1:upf_nqlc,i,j, 2)
     ENDDO
  ENDDO
  !
  !pp_pswfc

  ALLOCATE (upf_chi(upf_mesh,upf_ntwfc) )

  IF (ntwfc(1)==ntwfc(2)) THEN

     DO i=1,ntwfc(2)
        IF (interpolate) THEN
           WRITE (*,*) " interpolate chi"
           aux2(1,1:mesh(2)) = chi(1:mesh(2),i,2)
           CALL dosplineint( r(1:mesh(2),2), aux2, upf_r(1:upf_mesh), aux1 )
           chi(1:upf_mesh,i,2) = aux1(1,1:upf_mesh)
           WRITE (*,*) " done"
        ENDIF
        ! Jivtesh - The wavefunctions are calcuated to be the average of the
        ! wavefunctions of the two atoms - lines 365-366
        upf_chi(1:upf_mesh,i) =    x     * chi(1:upf_mesh,i,1) + &
                                (1.d0-x) * chi(1:upf_mesh,i,2)
     ENDDO
  ELSE
     WRITE (*,*) "Number of wavefunctions not the same for the two pseudopotentials"
  ENDIF
  !upf_chi(1:upf_mesh,1:upf_ntwfc) = chi(1:upf_mesh,1:upf_ntwfc,1)
  !
  !pp_rhoatm

  ALLOCATE (upf_rho_at(upf_mesh) )
  IF (interpolate) THEN
     WRITE (*,*) " interpolate rho_at"
     aux2(1,1:mesh(2)) = rho_at(1:mesh(2),2)
     CALL dosplineint( r(1:mesh(2),2), aux2, upf_r(1:upf_mesh), aux1 )
     rho_at(1:upf_mesh,2) = aux1(1,1:upf_mesh)
     WRITE (*,*) " done"
  ENDIF
  upf_rho_at(1:upf_mesh) =    x     * rho_at(1:upf_mesh,1) + &
                           (1.d0-x) * rho_at(1:upf_mesh,2)

END SUBROUTINE compute_virtual
