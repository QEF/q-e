!
! Copyright (C) 2001-2009 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file 'License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

! This utility can be used to convert norm-conserving
! pseudopotentials from FPMD old format to UPF format
!
! Usage:
!   fpmd2upf.x < input.namelist
!
!   input.namelist should contain the namelist fpmd_pseudo
!   fpmd_pseudo parameter are:
!
!   psfile   pseudopotential filename in FPMD format
!   nwfs     Number of wavefunction
!   wfl(i)   i = 1, nwfs  Wavefunction label
!   wfoc(i)  i = 1, nwfs  Wavefunction occupation
!   psd      element name
!   zp       valence charge
!   iexch    exchange functional
!   icorr    correlation functional
!   igcx     exchange gradient correction
!   igcc     correlation gradient correction
!
! Example:


MODULE fpmd2upf_module

  USE kinds, ONLY: dp
  USE parameters
  USE radial_grids, ONLY: ndmx

  IMPLICIT NONE
  SAVE

  REAL(dp), PRIVATE :: TOLMESH = 1.d-5

  TYPE pseudo_ncpp
     CHARACTER(len=4) :: psd         ! Element label
     CHARACTER(len=20) :: pottyp     ! Potential type
     LOGICAL :: tmix
     LOGICAL :: tnlcc
     INTEGER :: igau
     INTEGER :: lloc
     INTEGER :: nbeta
     INTEGER :: lll(lmaxx+1)
     INTEGER :: nchan
     INTEGER :: mesh
     REAL(dp) ::  zv
     REAL(dp) ::  dx            ! r(i) = cost * EXP( xmin + dx * (i-1) )
     REAL(dp) ::  rab(ndmx)
     REAL(dp) ::  rw(ndmx)
     REAL(dp) ::  vnl(ndmx, lmaxx+1)
     REAL(dp) ::  vloc(ndmx)
     REAL(dp) ::  vrps(ndmx, lmaxx+1)
     REAL(dp) ::  wgv(lmaxx+1)
     REAL(dp) ::  rc(2)
     REAL(dp) ::  wrc(2)
     REAL(dp) ::  rcl(3,3)
     REAL(dp) ::  al(3,3)
     REAL(dp) ::  bl(3,3)
     INTEGER :: nrps                     ! number of atomic wave function
     INTEGER :: lrps(lmaxx+1)            ! angular momentum
     REAL(dp) :: oc(lmaxx+1)            ! occupation for each rps
     REAL(dp) :: rps(ndmx, lmaxx+1)  ! atomic pseudo wave function
     REAL(dp) :: rhoc(ndmx)          ! core charge
   END TYPE pseudo_ncpp


CONTAINS


  SUBROUTINE read_pseudo_fpmd( ap, psfile )
    TYPE(pseudo_ncpp) :: ap
    CHARACTER(len=256) :: psfile
    CHARACTER(len=80) :: error_msg
    INTEGER :: info, iunit

    iunit = 11
    OPEN(UNIT=iunit,FILE=psfile,STATUS='OLD')
    REWIND( iunit )

    CALL read_head_pp( iunit, ap, error_msg, info)
    IF( info /= 0 ) GOTO 200

    IF( ap%pottyp == 'GIANNOZ' ) THEN

      CALL read_giannoz(iunit, ap, info)
      IF( info /= 0 ) GOTO 200

    ELSEIF( ap%pottyp == 'NUMERIC' ) THEN

      CALL read_numeric_pp( iunit, ap, error_msg, info)
      IF( info /= 0 ) GOTO 200

    ELSEIF( ap%pottyp == 'ANALYTIC' ) THEN

      CALL read_analytic_pp( iunit, ap, error_msg, info)
      IF( info /= 0 ) GOTO 200

    ELSE

      info = 1
      error_msg = ' Pseudopotential type '//trim(ap%pottyp)//' not implemented '
      GOTO 200

    ENDIF
200 CONTINUE
    IF( info /= 0 ) THEN
      CALL errore(' readpseudo ', error_msg, abs(info) )
    ENDIF

    CLOSE(iunit)

    RETURN
  END SUBROUTINE read_pseudo_fpmd


!=----------------------------------------------------------------------------=!

      SUBROUTINE analytic_to_numeric(ap)
        TYPE (pseudo_ncpp), INTENT(inout) :: ap
        INTEGER :: ir, mesh, lmax, l, n, il, ib, ll
        REAL(dp) :: xmin, zmesh, dx, x
!        REAL(dp) :: pi        = 3.14159265358979323846_dp

! ...   declare external function
        REAL(dp), EXTERNAL :: qe_erf

        IF( ap%mesh == 0 ) THEN
! ...     Local pseudopotential, define a logaritmic grid
          mesh  = size( ap%rw )
          xmin  = -5.0d0
          zmesh = 6.0d0
          dx    =  0.025d0
          DO ir = 1, mesh
            x = xmin + dble(ir-1) * dx
            ap%rw(ir)  = exp(x) / zmesh
            IF( ap%rw(ir) > 1000.0d0 ) exit
          ENDDO
          ap%mesh = mesh
          ap%dx   = dx
          ap%rab  = ap%dx * ap%rw
        ENDIF

        ap%vnl = 0.0d0
        ap%vloc = 0.0d0
        ap%vrps = 0.0d0
        DO l = 1, 3
          DO ir = 1, ap%mesh
            ap%vnl(ir,l)= - ( ap%wrc(1) * qe_erf(sqrt(ap%rc(1))*ap%rw(ir)) +  &
                              ap%wrc(2) * qe_erf(sqrt(ap%rc(2))*ap%rw(ir)) ) *&
                            ap%zv / ap%rw(ir)
          ENDDO
          DO ir = 1, ap%mesh
            DO n = 1, ap%igau
              ap%vnl(ir,l)= ap%vnl(ir,l)+ (ap%al(n,l)+ ap%bl(n,l)*ap%rw(ir)**2 )* &
                   exp(-ap%rcl(n,l)*ap%rw(ir)**2)
            ENDDO
          ENDDO
        ENDDO

! ...   Copy local component to a separate array
        ap%vloc(:) = ap%vnl(:,ap%lloc)
        DO l = 1, ap%nbeta
          ll=ap%lll(l) + 1  ! find out the angular momentum (ll-1) of the component stored
                         ! in position l
          ap%vrps(:,l) = ( ap%vnl(:,ll) - ap%vloc(:) ) * ap%rps(:,ll)
        ENDDO

        RETURN
      END SUBROUTINE analytic_to_numeric

!=----------------------------------------------------------------------------=!

      SUBROUTINE read_giannoz(uni, ap, ierr)
!        USE constants, ONLY : fpi
        IMPLICIT NONE
        TYPE (pseudo_ncpp), INTENT(inout) :: ap
        INTEGER, INTENT(in) :: uni
        INTEGER, INTENT(out) :: ierr
        REAL(dp) :: chi( size(ap%rps, 1), size(ap%rps, 2) )
        REAL(dp) :: vnl( size(ap%vnl, 1), size(ap%vnl, 2) )
        REAL(dp) :: rho_core( size(ap%rhoc, 1) )
        REAL(dp) :: r, ra, rb, fac
        REAL(dp) :: oc( size(ap%rps, 2) )
        REAL(dp) :: enl( size(ap%rps, 2) )
        REAL(dp) :: zmesh, xmin, dx, etot
        REAL(dp) :: zval
        INTEGER   :: nn(size(ap%rps, 2)), ll(size(ap%rps, 2))
        INTEGER   :: nwf, mesh, i, j, in1, in2, in3, in4, m
        INTEGER   :: lmax, nlc, nnl, lloc, l, il
        LOGICAL   :: nlcc
        CHARACTER(len=80) :: dft
        CHARACTER(len=4)  :: atom
        CHARACTER(len=2)  :: el( size(ap%rps, 2) )
        CHARACTER(len=80) :: ppinfo
        CHARACTER(len=80) :: strdum
        CHARACTER(len=2) :: sdum1, sdum2

!
        ierr = 0

        READ(uni,fmt='(a)') dft
        READ(uni,fmt='(a4,f5.1,3i2,a2,l1,a2,i2,a)') &
          atom, zval, lmax, nlc, nnl, sdum1, nlcc, sdum2, lloc, ppinfo

        ! WRITE(6,*) ' DEBUG ', atom, zval,lmax, nlc, nnl, nlcc, lloc, ppinfo

        IF( (lmax+1) > size(ap%vnl, 2) ) THEN
          ierr = 1
          RETURN
        ENDIF
        IF( (nlcc .and. .not.ap%tnlcc) .or. (.not.nlcc .and. ap%tnlcc) ) THEN
          ierr = 2
          RETURN
        ENDIF

        READ(uni,fmt='(f8.2,f8.4,f10.6,2i6)') zmesh, xmin, dx, mesh, nwf

        IF( mesh > size(ap%rps, 1) ) THEN
          ierr = 3
          RETURN
        ENDIF
        IF( nwf > size(ap%rps, 2) ) THEN
          ierr = 4
          RETURN
        ENDIF

        DO j = 0, lmax
           READ(uni,fmt="(A16,i1)") strdum, l
           READ(uni,'(4e16.8)') (vnl(i,j+1), i=1,mesh)
        ENDDO
        IF (nlcc) THEN
          READ(uni,fmt='(4e16.8)') (rho_core(i), i=1,mesh)
        ENDIF
        DO j = 1, nwf
          READ(uni,fmt="(A16,a2)") strdum,el(j)
          READ(uni,fmt='(i5,f6.2)') ll(j),oc(j)
          READ(uni,fmt='(4e16.8)') (chi(i,j), i=1,mesh)
        ENDDO

        ap%zv = zval
        ap%nchan = lmax+1
        ap%mesh = mesh
        ap%rw = 0.0d0
        ap%vnl = 0.0d0
        ap%vrps = 0.0d0
        fac = 0.5d0

        ! WRITE(6,*) ' DEBUG ', ap%lloc, ap%numeric, ap%nbeta, ap%raggio, ap%zv

        DO i = 1, mesh
          r = exp(xmin+dble(i-1)*dx)/zmesh
          ap%rw(i) = r
          DO j = 1, lmax+1
            ap%vnl(i,j) = vnl(i,j) * fac
          ENDDO
        ENDDO
        IF( minval( ap%rw(1:mesh) ) <= 0.0d0 ) THEN
           ierr = 5
           RETURN
        ENDIF
        ap%dx  = dx
        ap%rab = ap%dx * ap%rw
        ap%vloc(:) = ap%vnl(:,ap%lloc)

        ap%lrps(1:nwf) = ll(1:nwf)
        ap%oc = 0.0d0
        ap%nrps = nwf
        ap%mesh = mesh
        ap%rps = 0.0d0
!        fac = 1.0d0/SQRT(fpi)
        fac = 1.0d0
        DO i = 1, mesh
          r = exp(xmin+dble(i-1)*dx)/zmesh
          DO j = 1, nwf
            ap%rps(i,j) = chi(i,j) * fac
          ENDDO
        ENDDO

        DO l = 1, ap%nbeta
          il=ap%lll(l) + 1  ! find out the angular momentum (il-1) of the component stored
                         ! in position l
          DO i = 1, mesh
            ap%vrps(i,l) = ( ap%vnl(i,il) - ap%vloc(i) ) * ap%rps(i,il)
          ENDDO
        ENDDO

        IF( nlcc ) THEN
          ap%rhoc = 0.0d0
          DO i = 1, mesh
            r = exp(xmin+dble(i-1)*dx)/zmesh
            ap%rhoc(i) = rho_core(i)
          ENDDO
        ENDIF

        RETURN
      END SUBROUTINE read_giannoz

!=----------------------------------------------------------------------------=!


      SUBROUTINE ap_info( ap )
        TYPE (pseudo_ncpp), INTENT(in) :: ap
        INTEGER   :: in1, in2, in3, in4, m, il, ib, l, i

        IF (ap%nbeta > 0) THEN
          WRITE(6,10) ap%pottyp
          IF (ap%tmix) THEN
            WRITE(6,107)
            WRITE(6,106)  (ap%lll(l),l=1,ap%nbeta)
            WRITE(6,105)  (ap%wgv(l),l=1,ap%nbeta)
          ELSE
            WRITE(6,50) ap%lloc
          ENDIF
          WRITE(6,60) (ap%lll(l),l=1,ap%nbeta)
        ELSE
! ...     A local pseudopotential has been read.
          WRITE(6,11) ap%pottyp
          WRITE(6,50) ap%lloc
        ENDIF

   10   FORMAT(   3X,'Type is ',A10,' and NONLOCAL. ')
  107   FORMAT(   3X,'Mixed reference potential:')
  106   FORMAT(   3X,'  L     :',3(9X,i1))
  105   FORMAT(   3X,'  Weight:',3(2X,F8.5))
   50   FORMAT(   3X,'Local component is ..... : ',I3)
   60   FORMAT(   3X,'Non local components are : ',4I3)
   11   FORMAT(   3X,'Type is ',A10,' and LOCAL. ')
   20   FORMAT(   3X,'Pseudo charge : ',F8.3,', pseudo radius : ',F8.3)

        WRITE(6,20) ap%zv

        IF( ap%pottyp /= 'ANALYTIC' ) THEN

          WRITE(6,131) ap%nchan, ap%mesh, ap%dx
          in1=1
          in2=ap%mesh/4
          in3=ap%mesh/2
          in4=ap%mesh
          WRITE(6,132)
          WRITE(6,120) in1,ap%rw(in1),(ap%vnl(in1,m),m=1,ap%nchan)
          WRITE(6,120) in2,ap%rw(in2),(ap%vnl(in2,m),m=1,ap%nchan)
          WRITE(6,120) in3,ap%rw(in3),(ap%vnl(in3,m),m=1,ap%nchan)
          WRITE(6,120) in4,ap%rw(in4),(ap%vnl(in4,m),m=1,ap%nchan)
  131     FORMAT(/, 3X,'Pseudopotentials Grid    : Channels = ',I2,&
                   ', Mesh = ',I5,/,30X,'dx   = ',F16.14)
  132     FORMAT(   3X,'point      radius        pseudopotential')
  120     FORMAT(I8,F15.10,5F10.6)

        ELSE

          WRITE(6,25) ap%igau
          WRITE(6,30)
          WRITE(6,104) ap%wrc(1),ap%rc(1),ap%wrc(2),ap%rc(2)
   25     FORMAT(/, 3X,'Gaussians used : ',I2,'. Parameters are : ')
   30     FORMAT(   3X,'C (core), Alfa(core) : ')
  104     FORMAT(4(3X,F8.4))

          WRITE(6,40)
          DO il=1,3
            DO ib=1,ap%igau
              WRITE(6,103) ap%rcl(ib,il),ap%al(ib,il),ap%bl(ib,il)
            ENDDO
          ENDDO
   40     FORMAT(   3X,'Hsc radii and coeff. A and B :')
  103     FORMAT(3X,F8.4,2(3X,F15.7))


        ENDIF

        IF( ap%nrps > 0 .and. ap%mesh > 0 ) THEN
          WRITE(6,141) ap%nrps, ap%mesh, ap%dx
          in1=1
          in2=ap%mesh/4
          in3=ap%mesh/2
          in4=ap%mesh
          WRITE(6,145) (ap%oc(i),i=1,ap%nrps)
          WRITE(6,142)
          WRITE(6,120) in1,ap%rw(in1),(ap%rps(in1,m),m=1,ap%nrps)
          WRITE(6,120) in2,ap%rw(in2),(ap%rps(in2,m),m=1,ap%nrps)
          WRITE(6,120) in3,ap%rw(in3),(ap%rps(in3,m),m=1,ap%nrps)
          WRITE(6,120) in4,ap%rw(in4),(ap%rps(in4,m),m=1,ap%nrps)
        ENDIF

  141   FORMAT(/, 3X,'Atomic wavefunction Grid : Channels = ',I2,&
                   ', Mesh = ',I5,/,30X,'dx   = ',F16.14)
  142   FORMAT(   3X,'point      radius        wavefunction')
  145   FORMAT(   3X,'Channels occupation number : ',5F10.4)

        IF( ap%tnlcc ) THEN
          WRITE(6,151) ap%mesh, ap%dx
          in1 = 1
          in2 = ap%mesh / 4
          in3 = ap%mesh / 2
          in4 = ap%mesh
          WRITE(6,152)
          WRITE(6,120) in1,ap%rw(in1),ap%rhoc(in1)
          WRITE(6,120) in2,ap%rw(in2),ap%rhoc(in2)
          WRITE(6,120) in3,ap%rw(in3),ap%rhoc(in3)
          WRITE(6,120) in4,ap%rw(in4),ap%rhoc(in4)
        ENDIF

  151   FORMAT(/, 3X,'Core correction Grid     : Mesh = ',I5, &
             ', dx   = ',F16.14)
  152   FORMAT(   3X,'point      radius        rho core')

        RETURN
      END SUBROUTINE ap_info

!=----------------------------------------------------------------------------=!

      REAL(dp) FUNCTION calculate_dx( a, m )
        REAL(dp), INTENT(in) :: a(:)
        INTEGER, INTENT(in) :: m
        INTEGER :: n
        REAL(dp) :: ra, rb
          n = min( size( a ), m )
          ra = a(1)
          rb = a(n)
          calculate_dx = log( rb / ra ) / dble( n - 1 )
          WRITE(6,*) 'amesh (dx) = ', calculate_dx
        RETURN
      END FUNCTION calculate_dx


SUBROUTINE read_atomic_wf( iunit, ap, err_msg, ierr)
  USE parser, ONLY: field_count
  IMPLICIT NONE
  INTEGER, INTENT(in) :: iunit
  TYPE (pseudo_ncpp), INTENT(inout) :: ap
  CHARACTER(len=*) :: err_msg
  INTEGER, INTENT(out) :: ierr
!
  CHARACTER(len=80) :: input_line
  INTEGER :: i, j, m, strlen, info, nf, mesh
  REAL(dp) :: rdum

! ... read atomic wave functions
! ... nchan : indicate number of atomic wave functions ( s p d )

  ierr = 0
  err_msg = ' error while reading atomic wf '

  ap%rps  = 0.0_dp
  ap%nrps = 0
  ap%oc   = 0.0d0
  ap%lrps = 0

  ! this is for local pseudopotentials
  IF( ap%nbeta == 0 ) RETURN

  READ(iunit,'(A80)',end=100) input_line
  CALL field_count(nf, input_line)

  strlen = len_trim(input_line)

  IF( nf == 2 ) THEN
    READ(input_line(1:strlen),*,IOSTAT=ierr) mesh, ap%nrps
  ELSE
    READ(input_line(1:strlen),*,IOSTAT=ierr) mesh, ap%nrps, ( ap%oc(j), j=1, min(ap%nrps,size(ap%oc)) )
  ENDIF
  IF( ap%nrps > size(ap%rps,2) ) THEN
    ierr = 2
    WRITE( 6, * ) ' nchan = (wf) ', ap%nrps
    err_msg = ' NCHAN NOT PROGRAMMED '
    GOTO 110
  ENDIF
  IF( mesh > size(ap%rw) .or. mesh < 0) THEN
    ierr = 4
    err_msg = ' WAVMESH OUT OF RANGE '
    GOTO 110
  ENDIF

  DO j = 1, mesh
    READ(iunit,*,IOSTAT=ierr) rdum, (ap%rps(j,m),m=1,ap%nrps)
    IF( ap%mesh == 0 ) ap%rw(j) = rdum
    IF( abs(rdum - ap%rw(j))/(rdum+ap%rw(j)) > TOLMESH ) THEN
      ierr = 5
      err_msg = ' radial meshes do not match '
      GOTO 110
    ENDIF
  ENDDO

  IF( ap%mesh == 0 ) THEN
    ap%mesh = mesh
    ap%dx = calculate_dx( ap%rw, ap%mesh )
    ap%rab  = ap%dx * ap%rw
  ENDIF

  GOTO 110
100 ierr = 1
110 CONTINUE

  RETURN
END SUBROUTINE read_atomic_wf

!=----------------------------------------------------------------------------=!

SUBROUTINE read_numeric_pp( iunit, ap, err_msg, ierr)
  IMPLICIT NONE
  INTEGER, INTENT(in) :: iunit
  TYPE (pseudo_ncpp), INTENT(inout) :: ap
  CHARACTER(len=*) :: err_msg
  INTEGER, INTENT(out) :: ierr
!
  CHARACTER(len=80) :: input_line
  INTEGER :: i, j, m, strlen, info, nf, l, ll

! ... read numeric atomic pseudopotential
! ... nchan : indicate number of atomic wave functions ( s p d )

  ierr = 0
  err_msg = ' error while reading atomic numeric pseudo '

  IF(ap%tmix) THEN
    READ(iunit,*) (ap%wgv(l),l=1,ap%nbeta)
  ENDIF

  READ(iunit,*,IOSTAT=ierr) ap%zv
  READ(iunit,*,IOSTAT=ierr) ap%mesh, ap%nchan

  IF((ap%nchan > size(ap%vnl,2) ) .or. (ap%nchan < 1)) THEN
    ierr = 1
    WRITE( 6, * ) ' nchan (pp) = ', ap%nchan
    err_msg = ' NCHAN NOT PROGRAMMED '
    GOTO 110
  ENDIF
  IF((ap%mesh > size(ap%rw) ) .or. (ap%mesh < 0)) THEN
    info = 2
    err_msg = ' NPOTMESH OUT OF RANGE '
    GOTO 110
  ENDIF

  ap%rw = 0.0d0
  ap%vnl = 0.0d0
  ap%vloc = 0.0d0
  ap%vrps = 0.0d0
  DO j = 1, ap%mesh
    READ(iunit,*,IOSTAT=ierr) ap%rw(j), (ap%vnl(j,l),l=1,ap%nchan)
  ENDDO

  IF( minval( ap%rw(1:ap%mesh) ) <= 0.0d0 ) THEN
    info = 30
    err_msg = ' ap rw too small '
    GOTO 110
  ENDIF

! ...  mixed reference potential is in vr(lloc)
  IF(ap%tmix) THEN
    DO j=1,ap%mesh
      ap%vnl(j,ap%lloc)= 0.d0
      DO l=1,ap%nchan
        IF(l /= ap%lloc) THEN
          ap%vnl(j,ap%lloc)=  ap%vnl(j,ap%lloc) + ap%wgv(l) * ap%vnl(j,l)
        ENDIF
      ENDDO
    ENDDO
  ENDIF
  ap%vloc(:) = ap%vnl(:,ap%lloc)
  ap%dx = calculate_dx( ap%rw, ap%mesh )
  ap%rab  = ap%dx * ap%rw

  CALL read_atomic_wf( iunit, ap, err_msg, ierr)
  IF( ierr /= 0 ) GOTO 110

  DO l = 1, ap%nbeta
    ll=ap%lll(l) + 1
    ap%vrps(:,l) = ( ap%vnl(:,ll) - ap%vloc(:) ) * ap%rps(:,ll)
  ENDDO

  IF(ap%tnlcc) THEN
    CALL read_atomic_cc( iunit, ap,  err_msg, ierr)
    IF( ierr /= 0 ) GOTO 110
  ENDIF

  GOTO 110
100 ierr = 1
110 CONTINUE

  RETURN
END SUBROUTINE read_numeric_pp

!

SUBROUTINE read_head_pp( iunit, ap, err_msg, ierr)
  IMPLICIT NONE
  INTEGER, INTENT(in) :: iunit
  TYPE (pseudo_ncpp), INTENT(inout) :: ap
  CHARACTER(len=*) :: err_msg
  INTEGER, INTENT(out) :: ierr
!
  INTEGER :: i, l

! ... read pseudo header

  ierr = 0
  err_msg = ' error while reading header pseudo '

  ap%lll = 0
  READ(iunit, *) ap%tnlcc, ap%tmix
  READ(iunit, *) ap%pottyp, ap%lloc, ap%nbeta, (ap%lll(l), l = 1, min(ap%nbeta, size(ap%lll)) )

  ap%lll = ap%lll - 1

  IF( ap%nbeta > size(ap%lll) .or. ap%nbeta < 0 ) THEN
    ierr = 1
    err_msg = 'LNL out of range'
    GOTO 110
  ENDIF
  IF( ap%lloc < 0 .or. ap%lloc > size(ap%vnl,2) ) THEN
    ierr = 3
    err_msg = 'LLOC out of range'
    GOTO 110
  ENDIF
  IF( ap%tmix .and. ap%pottyp /= 'NUMERIC' ) THEN
    ierr = 4
    err_msg = 'tmix not implemented for pseudo ' // ap%pottyp
    GOTO 110
  ENDIF
  DO l = 2, ap%nbeta
    IF( ap%lll(l) <= ap%lll(l-1)) THEN
      ierr = 5
      err_msg =' NONLOCAL COMPONENTS MUST BE GIVEN IN ASCENDING ORDER'
      GOTO 110
    ENDIF
  ENDDO
  DO l = 1, ap%nbeta
    IF( ap%lll(l)+1 == ap%lloc) THEN
      ierr = 6
      err_msg = ' LLOC.EQ.L NON LOCAL!!'
      GOTO 110
    ENDIF
  ENDDO

  GOTO 110
100 ierr = 1
110 CONTINUE

  RETURN
END SUBROUTINE read_head_pp

!=----------------------------------------------------------------------------=!

SUBROUTINE read_analytic_pp( iunit, ap, err_msg, ierr)
  IMPLICIT NONE
  INTEGER, INTENT(in) :: iunit
  TYPE (pseudo_ncpp), INTENT(inout) :: ap
  CHARACTER(len=*) :: err_msg
  INTEGER, INTENT(out) :: ierr
!
  INTEGER :: i, l

! ... read analytic pseudo gaussians

  ierr = 0
  err_msg = ' error while reading atomic analytic pseudo '

  READ(iunit,*,IOSTAT=ierr) ap%zv, ap%igau

  ap%mesh = 0
  ap%nchan = 0
  ap%dx = 0.0d0
  ap%rab  = 0.0d0
  ap%rw   = 0.0d0
  ap%vnl   = 0.0d0
  ap%vloc   = 0.0d0
  ap%vrps   = 0.0d0

  SELECT CASE (ap%igau)
    CASE ( 1 )
      READ(iunit,*,IOSTAT=ierr) ap%rc(1)
      ap%wrc(1) = 1.d0
      ap%wrc(2) = 0.d0
      ap%rc(2)  = 0.d0
    CASE ( 3 )
      READ(iunit,*,IOSTAT=ierr) ap%wrc(1), ap%rc(1), ap%wrc(2), ap%rc(2)
    CASE DEFAULT
      ierr = 1
      err_msg = ' IGAU NOT PROGRAMMED '
      GOTO 110
  END SELECT

  DO l=1,3
    DO i=1,ap%igau
      READ(iunit,*,IOSTAT=ierr) ap%rcl(i,l), ap%al(i,l), ap%bl(i,l)
    ENDDO
  ENDDO

  CALL read_atomic_wf( iunit, ap, err_msg, ierr)
  IF( ierr /= 0 ) GOTO 110

  IF(ap%tnlcc) THEN
    CALL read_atomic_cc( iunit, ap, err_msg, ierr)
    IF( ierr /= 0 ) GOTO 110
  ENDIF

! ... Analytic pseudo are not supported anymore, conversion
! ... to numeric form is forced
  CALL analytic_to_numeric( ap )

  GOTO 110
100 ierr = 1
110 CONTINUE

  RETURN
END SUBROUTINE read_analytic_pp

!=----------------------------------------------------------------------------=!


SUBROUTINE read_atomic_cc( iunit, ap, err_msg, ierr)
  IMPLICIT NONE
  INTEGER, INTENT(in) :: iunit
  TYPE (pseudo_ncpp), INTENT(inout) :: ap
  CHARACTER(len=*) :: err_msg
  INTEGER, INTENT(out) :: ierr
!
  CHARACTER(len=80) :: input_line
  INTEGER :: j, mesh
  REAL(dp) :: rdum

! ... read atomic core

  ierr = 0
  err_msg = ' error while reading atomic core pseudo '

  ap%rhoc = 0.0d0

  READ(iunit,*,IOSTAT=ierr) mesh
  IF(mesh > size(ap%rw) .or. mesh < 0 ) THEN
    ierr = 17
    err_msg = '  CORE CORRECTION MESH OUT OF RANGE '
    GOTO 110
  ENDIF
  DO j = 1, mesh
    READ(iunit,*,IOSTAT=ierr) rdum, ap%rhoc(j)
    IF( ap%mesh == 0 ) ap%rw(j) = rdum
    IF( abs(rdum - ap%rw(j))/(rdum+ap%rw(j)) > TOLMESH ) THEN
      ierr = 5
      err_msg = ' core cor. radial mesh does not match '
      GOTO 110
    ENDIF
  ENDDO

  IF( ap%mesh == 0 ) THEN
    ap%mesh = mesh
    ap%dx = calculate_dx( ap%rw, ap%mesh )
    ap%rab  = ap%dx * ap%rw
  ENDIF

  GOTO 110
100 ierr = 1
110 CONTINUE

  RETURN
END SUBROUTINE read_atomic_cc


END MODULE fpmd2upf_module




PROGRAM fpmd2upf

  !
  !     Convert a pseudopotential written in the FPMD format
  !     to unified pseudopotential format
  !

  USE kinds
  USE fpmd2upf_module
  USE parameters
  USE upf

  IMPLICIT NONE

  TYPE (pseudo_ncpp) :: ap
  CHARACTER(len=256) :: psfile
  CHARACTER(len=2)   :: wfl( 10 )
  REAL(8)       :: wfoc( 10 )
  INTEGER :: nsp, nspnl, i, lloc, l, ir, iv, kkbeta
  REAL(8) :: rmax = 10
  REAL(8) :: vll
  REAL(8), ALLOCATABLE :: aux(:)

  NAMELIST / fpmd_pseudo / psfile, nwfs, wfl, wfoc, psd, &
                           iexch, icorr, igcx, igcc, zp

  ! ... end of declarations

  CALL input_from_file()

  READ( 5, fpmd_pseudo )

  nsp = 1
  CALL read_pseudo_fpmd(ap, psfile)

  WRITE(generated, '("Generated using unknown code")')
  WRITE(date_author,'("Author: unknown    Generation date: as well")')
  comment = 'Info: automatically converted from CPMD format'

  rcloc = 0.0d0

  ALLOCATE( els(nwfs), oc(nwfs), epseu(nwfs) )
  ALLOCATE( lchi(nwfs), nns(nwfs) )
  ALLOCATE( rcut (nwfs), rcutus (nwfs) )

  els = '?'
  oc  = 0.0d0
  DO i = 1, nwfs
     els(i)   = wfl(i)
     oc(i)    = wfoc(i)
     lchi(i)  = i - 1
     nns (i)  = 0
     rcut(i)  = 0.0d0
     rcutus(i)= 0.0d0
     epseu(i) = 0.0d0
  ENDDO

  pseudotype = 'NC'
  nlcc       = ap%tnlcc
  IF( ap%zv > 0.0d0 ) zp = ap%zv
  etotps     = 0.0d0

  lloc  = ap%lloc
  lmax  = max( maxval( ap%lll( 1:ap%nbeta ) ), ap%lloc - 1 )

  nbeta = ap%nbeta
  mesh  = ap%mesh
  ntwfc = nwfs
  ALLOCATE( elsw(ntwfc), ocw(ntwfc), lchiw(ntwfc) )
  DO i = 1, nwfs
     lchiw(i) = lchi(i)
     ocw(i)   = oc(i)
     elsw(i)  = els(i)
  ENDDO

  ALLOCATE(rab(mesh))
  ALLOCATE(  r(mesh))
  r   = ap%rw
  ap%dx = calculate_dx( ap%rw, ap%mesh )
  rab = ap%rw * ap%dx

  WRITE(6,*) ap%lloc, ap%lll( 1:ap%nbeta ) , ap%nbeta, ap%dx

  ALLOCATE (rho_atc(mesh))
  IF (nlcc) rho_atc = ap%rhoc

  ALLOCATE (vloc0(mesh))
  ! the factor 2 converts from Hartree to Rydberg
  vloc0(:) = ap%vloc * 2.0d0

  IF (nbeta > 0) THEN

     ALLOCATE(ikk2(nbeta), lll(nbeta))
     kkbeta = mesh
     DO ir = 1,mesh
        IF ( r(ir) > rmax ) THEN
           kkbeta=ir
           exit
        ENDIF
     ENDDO
     ikk2(:) = kkbeta
     ALLOCATE(aux(kkbeta))
     ALLOCATE(betar(mesh,nbeta))
     ALLOCATE(qfunc(mesh,nbeta,nbeta))
     ALLOCATE(dion(nbeta,nbeta))
     ALLOCATE(qqq (nbeta,nbeta))
     qfunc(:,:,:)=0.0d0
     dion(:,:) =0.d0
     qqq(:,:)  =0.d0
     iv = 0
     DO i = 1, nwfs
        l = lchi(i)
        IF ( l /= (lloc-1) ) THEN
           iv = iv + 1
           lll( iv ) = l
           DO ir = 1, kkbeta
              ! the factor 2 converts from Hartree to Rydberg
              betar(ir, iv) = 2.d0 * ap%vrps( ir, iv )
              aux(ir) = ap%rps(ir, (l+1) ) * betar(ir, iv)
           ENDDO
           CALL simpson2(kkbeta, aux(1), rab(1), vll)
           dion(iv,iv) = 1.0d0/vll
           WRITE(6,*) aux(2), rab(2), kkbeta, vll
        ENDIF
     ENDDO

  ENDIF

  ALLOCATE (rho_at(mesh))
  rho_at = 0.d0
  DO i = 1, nwfs
     rho_at(:) = rho_at(:) + ocw(i) * ap%rps(:, i) ** 2
  ENDDO

  ALLOCATE (chi(mesh,ntwfc))
  chi = ap%rps
  !     ----------------------------------------------------------
  WRITE (6,'(a)') 'Pseudopotential successfully converted'
  !     ----------------------------------------------------------

  CALL write_upf_v1( 10 )

100 CONTINUE

END PROGRAM fpmd2upf



!----------------------------------------------------------------------
SUBROUTINE simpson2(mesh,func,rab,asum)
  !-----------------------------------------------------------------------
!
  !     simpson's rule integrator for function stored on the
  !     radial logarithmic mesh
  !

  IMPLICIT NONE

  INTEGER :: i, mesh
  real(8) ::  rab(mesh), func(mesh), f1, f2, f3, r12, asum

      !     routine assumes that mesh is an odd number so run check
      !     if ( mesh+1 - ( (mesh+1) / 2 ) * 2 .ne. 1 ) then
      !       write(*,*) '***error in subroutine radlg'
!       write(*,*) 'routine assumes mesh is odd but mesh =',mesh+1
!       stop
!     endif

  asum = 0.0d0
  r12 = 1.0d0 / 12.0d0
  f3  = func(1) * rab(1) * r12

  DO i = 2,mesh-1,2
     f1 = f3
     f2 = func(i) * rab(i) * r12
     f3 = func(i+1) * rab(i+1) * r12
     asum = asum + 4.0d0*f1 + 16.0d0*f2 + 4.0d0*f3
  ENDDO

  RETURN
END SUBROUTINE simpson2



!
!  Description of the Native FPMD pseudopotential format
!
!  The format of the file must be as follows
!  (lowercase text and }'s are comments):
!
!  When POTTYP = 'ANALYTIC' the layout is:
!
!    TCC      TMIX                additional stuff on each line is ignored
!    POTTYP   LLOC LNL ( INDL(i), i = 1, LNL )
!    ( WGV(i), i = 1, LNL )       this line only if tmix(is) is true
!    ZV       IGAU                igau must be 1 or 3     }
!    WRC(1) RC(1) WRC(2) RC(2)    this line if igau = 3   }
!    RC(1)                        this one if igau = 1    }
!    RCL(1,1)    AL(1,1)    BL(1,1)         }             }  this
!     ...         ...        ...            }  l = 0      }  section
!    RCL(IGAU,1) AL(IGAU,1) BL(IGAU,1)      }             }  only if
!    RCL(1,2)    AL(1,2)    BL(1,2)      }                }  pottyp is
!     ...         ...        ...         }     l = 1      }  'ANALYTIC'
!    RCL(IGAU,2) AL(IGAU,2) BL(IGAU,2)   }                }
!    RCL(1,3)    AL(1,3)    BL(1,3)         }             }
!     ...         ...        ...            }  l = 2      }
!    RCL(IGAU,3) AL(IGAU,3) BL(IGAU,3)      }             }
!    NMESH NCHAN                                       }
!    RW( 1 )     ( RPS( 1, j ), j = 1, NCHAN )         }  pseudowave
!     ...         ...              ...                 }
!    RW( NMESH ) ( RPS( NMESH, j ), j = 1, NCHAN )     }
!
!
!  When POTTYP = 'NUMERIC' the layout is:
!
!    TCC      TMIX             additional stuff on each line is ignored
!    POTTYP   LLOC LNL  ( INDL(i), i = 1, LNL )
!    ( WGV(i), i = 1, LNL )       this line only if tmix(is) is true
!    ZV                                             }
!    NMESH NCHAN                                    }    this if
!    RW( 1 )     ( VR( 1, j ), j = 1, NCHAN )       }    pottyp is
!     ...       ...             ...                 }    'NUMERIC'
!    RW( NMESH ) ( VR( NMESH, j ), j = 1, NCHAN )   }
!    NMESH NCHAN                                       }
!    RW( 1 )     ( RPS( 1, j ), j = 1, NCHAN )         }  pseudowave
!     ...         ...              ...                 }
!    RW( NMESH ) ( RPS( NMESH, j ), j = 1, NCHAN )     }
!
!  DETAILED DESCRIPTION OF INPUT PARAMETERS:
!
!    TCC      (logical)   True if Core Correction are required for this
!                         pseudo
!
!    TMIX     (logical)   True if we want to mix nonlocal pseudopotential
!                         components
!
!    WGV(i)   (real)      wheight of the nonlocal components in the
!                         pseudopotential mixing scheme
!                         These parameters are present only if TMIX = .TRUE.
!                         1 <= i <= LNL
!
!    POTTYP   (character) pseudopotential type
!                         pottyp = 'ANALYTIC' : use an analytic expression
!                         pottyp = 'NUMERIC'  : read values from a table
!
!    ZV       (integer)   valence for each species
!
!    IGAU     (integer)   number of Gaussians in the pseudopotentials
!                         expression used only if pottyp='ANALYTIC'
!
!  parameters from Bachelet-Hamann-Schluter's table:
!
!    WRC(2)   (real)      c1, c2 (core)  parameters
!    RC(2)    (real)      alpha1, alpha2 parameters
!
!    RCL(i,3) (real)      alpha1, alpha2, alpha3 for each angular momentum
!                         1 <= i <= IGAU
!    AL(i,3)  (real)      parameters for each angular momentum
!                         1 <= i <= IGAU
!    BL(i,3)  (real)      parameters for each angular momentum
!                         1 <= i <= IGAU
!
!  nonlocality
!    IGAU     (integer)   number of Gaussians for analytic pseudopotentials
!    LLOC     (integer)   index of the angular momentum component added to
!                          the local part  ( s = 1, p = 2, d = 3 )
!    LNL      (integer)   number of non local component
!    INDL(i)  (integer)   indices of non local components
!                         1 <= i <= LNL
!                         ( 1 3 means s and d taken as non local )
!
!  pseudo grids
!    NMESH    (integer)   number of points in the mesh mesh
!    NCHAN    (integer)   numbero of colums, radial components
!    RW(i)    (real)      distance from the core in A.U. (radial mesh)
!                         1 <= i <= NMESH
!    RPS(i,j) (real)      Atomic pseudo - wavefunctions
!                         1 <= i <= NMESH ; 1 <= j <= NCHAN
!    VP(i,j)  (real)      Atomic pseudo - potential
!                         1 <= i <= NMESH ; 1 <= j <= NCHAN
!
!  ----------------------------------------------
!  END manual

