!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file 'License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!
module upf
  !
  ! All variables to be written into the UPF file
  ! (UPF = unified pseudopotential format)
  !
  ! pp_info
  integer :: rel
  real(kind=8) :: rcloc
  integer :: nwfs
  real(kind=8), allocatable :: oc(:), rcut(:), rcutus(:), epseu(:)
  character(len=2), allocatable :: els(:)
  integer, allocatable:: lchi (:), nns (:)
  !
  ! pp_header
  character (len=80):: generated, date_author, comment
  character (len=2) :: psd, pseudotype
  integer :: nv = 0
  integer :: iexch, icorr, igcx, igcc
  integer :: lmax, mesh, nbeta, ntwfc
  logical :: nlcc
  real(kind=8) :: zp, ecutrho, ecutwfc, etotps
  real(kind=8), allocatable :: ocw(:)
  character(len=2), allocatable :: elsw(:)
  integer, allocatable:: lchiw(:)
  !
  ! pp_mesh
  real(kind=8), allocatable :: r(:), rab(:)
  !
  ! pp_nlcc
  real(kind=8), allocatable :: rho_atc(:)
  !
  ! pp_local
  real(kind=8), allocatable ::  vloc0(:)
  !
  ! pp_nonlocal
  ! pp_beta
  real(kind=8), allocatable :: betar(:,:)
  integer, allocatable:: lll(:), ikk2(:)  
  ! pp_dij
  real(kind=8), allocatable :: dion(:,:)
  ! pp_qij
  integer ::  nqf, nqlc
  real(kind=8), allocatable :: rinner(:), qqq(:,:), qfunc(:,:,:)
  ! pp_qfcoef
  real(kind=8), allocatable :: qfcoef(:,:,:,:)
  !
  ! pp_pswfc
  real(kind=8), allocatable :: chi(:,:)
  !
  ! pp_rhoatom
  real(kind=8), allocatable :: rho_at(:)
end module upf
!

module fpmd2upf_module

  USE kinds, ONLY: dbl

  implicit none
  save

  REAL(dbl), PRIVATE :: TOLMESH = 1.d-5

contains

  subroutine read_pseudo_fpmd( ap, psfile )
    USE pseudo_types, ONLY: pseudo_ncpp
    type(pseudo_ncpp) :: ap
    character(len=80) :: psfile
    character(len=80) :: error_msg
    integer :: info, iunit

    iunit = 11
    OPEN(UNIT=iunit,FILE=psfile,STATUS='OLD')
    REWIND( iunit )

    CALL read_head_pp( iunit, ap, error_msg, info)
    IF( info /= 0 ) GO TO 200
    
    IF( ap%pottyp == 'GIANNOZ' ) THEN

      CALL read_giannoz(iunit, ap, info)
      IF( info /= 0 ) GO TO 200

    ELSE IF( ap%pottyp == 'NUMERIC' ) THEN

      CALL read_numeric_pp( iunit, ap, error_msg, info)
      IF( info /= 0 ) GO TO 200

    ELSE IF( ap%pottyp == 'ANALYTIC' ) THEN

      CALL read_analytic_pp( iunit, ap, error_msg, info)
      IF( info /= 0 ) GO TO 200

    ELSE

      info = 1
      error_msg = ' Pseudopotential type '//TRIM(ap%pottyp)//' not implemented '
      GO TO 200

    END IF
200 CONTINUE
    IF( info /= 0 ) THEN
      CALL error(' readpseudo ', error_msg, ABS(info) )
    END IF

    CLOSE(iunit)

    return 
  end subroutine


!=----------------------------------------------------------------------------=!

      SUBROUTINE check_file_type( iunit, info )

! ... This sub. check if a given fortran unit 'iunit' contains a UPF pseudopot.

        use parser, only: matches

        INTEGER, INTENT(IN) :: iunit
        INTEGER, INTENT(OUT) :: info
        CHARACTER(LEN=80) :: dummy
        INTEGER :: ios
        info = 0
        ios  = 0
        header_loop: do while (ios == 0)
          read (iunit, *, iostat = ios, err = 200) dummy  
          if (matches ("<PP_HEADER>", dummy) ) then
            info = 1
            exit header_loop
          endif
        enddo header_loop
200     continue
        RETURN
      END SUBROUTINE check_file_type

!=----------------------------------------------------------------------------=!

      SUBROUTINE analytic_to_numeric(ap)
        USE pseudo_types, ONLY: pseudo_ncpp
        TYPE (pseudo_ncpp), INTENT(INOUT) :: ap
        INTEGER :: ir, mesh, lmax, l, n, il, ib, ll
        REAL(dbl) :: xmin, zmesh, dx, x
        REAL(dbl) :: pi        = 3.14159265358979323846_dbl

! ...   declare external function
        REAL(dbl) :: erf, erfc
        EXTERNAL erf, erfc

        IF( ap%mesh == 0 ) THEN
! ...     Local pseudopotential, define a logaritmic grid
          mesh  = SIZE( ap%rw )
          xmin  = -5.0d0
          zmesh = 6.0d0
          dx    =  0.025d0
          DO ir = 1, mesh
            x = xmin + REAL(ir-1) * dx
            ap%rw(ir)  = EXP(x) / zmesh
          END DO
          ap%mesh = mesh
          ap%dx   = dx
          ap%rab  = ap%dx * ap%rw
        END IF

        ap%vnl = 0.0d0
        ap%vloc = 0.0d0
        ap%vrps = 0.0d0
        do l = 1, 3
          do ir = 1, ap%mesh
            ap%vnl(ir,l)= - ( ap%wrc(1) * erf(SQRT(ap%rc(1))*ap%rw(ir)) +     &
                             ap%wrc(2) * erf(SQRT(ap%rc(2))*ap%rw(ir)) ) * ap%zv / ap%rw(ir)
          end do
          do ir = 1, ap%mesh
            do n = 1, ap%igau
              ap%vnl(ir,l)= ap%vnl(ir,l)+ (ap%al(n,l)+ ap%bl(n,l)*ap%rw(ir)**2 )* &
                   EXP(-ap%rcl(n,l)*ap%rw(ir)**2)
            end do
          end do
        end do

! ...   Copy local component to a separate array
        ap%vloc(:) = ap%vnl(:,ap%lloc)
        DO l = 1, ap%lnl
          ll=ap%indl(l)  ! find out the angular momentum (ll-1) of the component stored
                         ! in position l
          ap%vrps(:,l) = ( ap%vnl(:,ll) - ap%vloc(:) ) * ap%rps(:,ll)
        END DO

        RETURN
      END SUBROUTINE

!=----------------------------------------------------------------------------=!

      SUBROUTINE read_giannoz(uni, ap, ierr)
        USE pseudo_types, ONLY: pseudo_ncpp
        IMPLICIT NONE
        TYPE (pseudo_ncpp), INTENT(INOUT) :: ap
        INTEGER, INTENT(IN) :: uni
        INTEGER, INTENT(OUT) :: ierr
        REAL(dbl) :: pi        = 3.14159265358979323846_dbl
        REAL(dbl) :: chi( SIZE(ap%rps, 1), SIZE(ap%rps, 2) )
        REAL(dbl) :: vnl( SIZE(ap%vnl, 1), SIZE(ap%vnl, 2) )
        REAL(dbl) :: rho_core( SIZE(ap%rhoc, 1) )
        REAL(dbl) :: r, ra, rb, fac
        REAL(dbl) :: oc( SIZE(ap%rps, 2) )
        REAL(dbl) :: enl( SIZE(ap%rps, 2) )
        REAL(dbl) :: zmesh, xmin, dx, etot
        REAL(dbl) :: zval
        INTEGER   :: nn(SIZE(ap%rps, 2)), ll(SIZE(ap%rps, 2))
        INTEGER   :: nwf, mesh, i, j, in1, in2, in3, in4, m
        INTEGER   :: lmax, nlc, nnl, lloc, l, il
        LOGICAL   :: nlcc
        CHARACTER(len=80) :: dft
        CHARACTER(len=4)  :: atom
        CHARACTER(len=2)  :: el( SIZE(ap%rps, 2) )
        CHARACTER(len=80) :: ppinfo
        CHARACTER(len=80) :: strdum
        CHARACTER(len=2) :: sdum1, sdum2

!
        ierr = 0

        READ(uni,fmt='(a)') dft
        READ(uni,fmt='(a4,f5.1,3i2,a2,l1,a2,i2,a)') &
          atom, zval, lmax, nlc, nnl, sdum1, nlcc, sdum2, lloc, ppinfo

        ! WRITE(6,*) ' DEBUG ', atom, zval,lmax, nlc, nnl, nlcc, lloc, ppinfo

        IF( (lmax+1) > SIZE(ap%vnl, 2) ) THEN
          ierr = 1
          RETURN
        END IF
        IF( (nlcc .AND. .NOT.ap%tnlcc) .OR. (.NOT.nlcc .AND. ap%tnlcc) ) THEN
          ierr = 2
          RETURN
        END IF

        READ(uni,fmt='(f8.2,f8.4,f10.6,2i6)') zmesh, xmin, dx, mesh, nwf

        IF( mesh > SIZE(ap%rps, 1) ) THEN
          ierr = 3
          RETURN
        END IF
        IF( nwf > SIZE(ap%rps, 2) ) THEN
          ierr = 4
          RETURN
        END IF

        DO j = 0, lmax
           READ(uni,fmt="(A16,i1)") strdum, l
           READ(uni,'(4e16.8)') (vnl(i,j+1), i=1,mesh)
        END DO
        IF (nlcc) THEN
          READ(uni,fmt='(4e16.8)') (rho_core(i), i=1,mesh)
        END IF   
        DO j = 1, nwf
          READ(uni,fmt="(A16,a2)") strdum,el(j)
          READ(uni,fmt='(i5,f6.2)') ll(j),oc(j)
          READ(uni,fmt='(4e16.8)') (chi(i,j), i=1,mesh)
        END DO

        ap%zv = zval
        ap%nchan = lmax+1
        ap%mesh = mesh
        ap%rw = 0.0d0
        ap%vnl = 0.0d0
        ap%vrps = 0.0d0
        fac = 0.5d0

        ! WRITE(6,*) ' DEBUG ', ap%lloc, ap%numeric, ap%lnl, ap%raggio, ap%zv

        DO i = 1, mesh
          r = EXP(xmin+REAL(i-1)*dx)/zmesh
          ap%rw(i) = r
          DO j = 1, lmax+1
            ap%vnl(i,j) = vnl(i,j) * fac
          END DO
        END DO
        IF( MINVAL( ap%rw(1:mesh) ) <= 0.0d0 ) THEN
           ierr = 5
           RETURN
        END IF
        ap%dx  = dx
        ap%rab = ap%dx * ap%rw
        ap%vloc(:) = ap%vnl(:,ap%lloc)

        ap%lrps(1:nwf) = ll(1:nwf)
        ap%oc = 0.0d0
        ap%nrps = nwf
        ap%mesh = mesh
        ap%rps = 0.0d0
        fac = 1.0d0/SQRT(4.0d0*pi)
        fac = 1.0d0
        DO i = 1, mesh
          r = EXP(xmin+REAL(i-1)*dx)/zmesh
          DO j = 1, nwf
            ap%rps(i,j) = chi(i,j) * fac
          END DO
        END DO

        DO l = 1, ap%lnl
          il=ap%indl(l)  ! find out the angular momentum (il-1) of the component stored
                         ! in position l
          DO i = 1, mesh
            ap%vrps(i,l) = ( ap%vnl(i,il) - ap%vloc(i) ) * ap%rps(i,il)
          END DO
        END DO

        IF( nlcc ) THEN
          ap%rhoc = 0.0d0
          DO i = 1, mesh
            r = EXP(xmin+REAL(i-1)*dx)/zmesh
            ap%rhoc(i) = rho_core(i)
          END DO
        END IF

        RETURN
      END SUBROUTINE 

!=----------------------------------------------------------------------------=!


      SUBROUTINE ap_info( ap )
        USE pseudo_types, ONLY: pseudo_ncpp
        TYPE (pseudo_ncpp), INTENT(IN) :: ap
        INTEGER   :: in1, in2, in3, in4, m, il, ib, l, i

        IF (ap%lnl > 0) THEN
          WRITE(6,10) ap%pottyp
          IF (ap%tmix) THEN
            WRITE(6,107) 
            WRITE(6,106)  (ap%indl(l),l=1,ap%lnl)
            WRITE(6,105)  (ap%wgv(l),l=1,ap%lnl)
          ELSE
            WRITE(6,50) ap%lloc
          END IF
          WRITE(6,60) (ap%indl(l),l=1,ap%lnl)
        ELSE
! ...     A local pseudopotential has been read.
          WRITE(6,11) ap%pottyp
          WRITE(6,50) ap%lloc
        END IF

   10   FORMAT(   3X,'Type is ',A10,' and NONLOCAL. ')
  107   FORMAT(   3X,'Mixed reference potential:')
  106   FORMAT(   3X,'  L     :',3(9X,i1))
  105   FORMAT(   3X,'  Weight:',3(2X,F8.5))
   50   FORMAT(   3X,'Local component is ..... : ',I3)
   60   FORMAT(   3X,'Non local components are : ',4I3)
   11   FORMAT(   3X,'Type is ',A10,' and LOCAL. ')
   20   FORMAT(   3X,'Pseudo charge : ',F8.3,', pseudo radius : ',F8.3)

        WRITE(6,20) ap%zv, ap%raggio

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
            END DO
          END DO
   40     FORMAT(   3X,'Hsc radii and coeff. A and B :')
  103     FORMAT(3X,F8.4,2(3X,F15.7))


        END IF

        IF( ap%nrps > 0 .AND. ap%mesh > 0 ) THEN
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
        END IF

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
        END IF

  151   FORMAT(/, 3X,'Core correction Grid     : Mesh = ',I5, &
             ', dx   = ',F16.14)
  152   FORMAT(   3X,'point      radius        rho core')

        RETURN
      END SUBROUTINE 

!=----------------------------------------------------------------------------=!

      REAL(dbl) FUNCTION calculate_dx( a, m )
        REAL(dbl), INTENT(IN) :: a(:)
        INTEGER, INTENT(IN) :: m 
        INTEGER :: n
        REAL(dbl) :: ra, rb 
          n = MIN( SIZE( a ), m )
          ra = a(1)
          rb = a(n)
          calculate_dx = LOG( rb / ra ) / REAL( n - 1 )
        RETURN
      END FUNCTION 


SUBROUTINE read_atomic_wf( iunit, ap, err_msg, ierr)
  USE pseudo_types, ONLY: pseudo_ncpp
  USE parser, ONLY: field_count
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: iunit
  TYPE (pseudo_ncpp), INTENT(INOUT) :: ap
  CHARACTER(LEN=*) :: err_msg
  INTEGER, INTENT(OUT) :: ierr
!
  CHARACTER(LEN=80) :: input_line
  INTEGER :: i, j, m, strlen, info, nf, mesh
  REAL(dbl) :: rdum

! ... read atomic wave functions
! ... nchan : indicate number of atomic wave functions ( s p d )

  ierr = 0
  err_msg = ' error while reading atomic wf '

  ap%rps  = 0.0_dbl
  ap%nrps = 0
  ap%oc   = 0.0d0
  ap%lrps = 0

  ! this is for local pseudopotentials
  IF( ap%lnl == 0 ) RETURN
              
  READ(iunit,'(A80)',end=100) input_line
  CALL field_count(nf, input_line)

  strlen = len_trim(input_line)

  IF( nf == 2 ) THEN
    READ(input_line(1:strlen),*,IOSTAT=ierr) mesh, ap%nrps
  ELSE
    READ(input_line(1:strlen),*,IOSTAT=ierr) mesh, ap%nrps, ( ap%oc(j), j=1, MIN(ap%nrps,SIZE(ap%oc)) )
  END IF
  IF( ap%nrps > SIZE(ap%rps,2) ) THEN
    ierr = 2   
    err_msg = ' NCHAN NOT PROGRAMMED '
    GO TO 110
  END IF
  IF( mesh > SIZE(ap%rw) .OR. mesh < 0) THEN
    ierr = 4
    err_msg = ' WAVMESH OUT OF RANGE '
    GO TO 110
  END IF

  DO j = 1, mesh
    READ(iunit,*,IOSTAT=ierr) rdum, (ap%rps(j,m),m=1,ap%nrps)
    IF( ap%mesh == 0 ) ap%rw(j) = rdum
    IF( ABS(rdum - ap%rw(j))/(rdum+ap%rw(j)) > TOLMESH ) THEN
      ierr = 5
      err_msg = ' radial meshes do not match '
      GO TO 110
    END IF
  END DO

  IF( ap%mesh == 0 ) THEN
    ap%mesh = mesh
    ap%dx = calculate_dx( ap%rw, ap%mesh )
    ap%rab  = ap%dx * ap%rw
  END IF

  GOTO 110
100 ierr = 1
110 CONTINUE
  
  RETURN
END SUBROUTINE

!=----------------------------------------------------------------------------=!

SUBROUTINE read_numeric_pp( iunit, ap, err_msg, ierr)
  USE pseudo_types, ONLY: pseudo_ncpp
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: iunit
  TYPE (pseudo_ncpp), INTENT(INOUT) :: ap
  CHARACTER(LEN=*) :: err_msg
  INTEGER, INTENT(OUT) :: ierr
!
  CHARACTER(LEN=80) :: input_line
  INTEGER :: i, j, m, strlen, info, nf, l, ll

! ... read numeric atomic pseudopotential
! ... nchan : indicate number of atomic wave functions ( s p d )

  ierr = 0
  err_msg = ' error while reading atomic numeric pseudo '

  IF(ap%tmix) THEN
    READ(iunit,*) (ap%wgv(l),l=1,ap%lnl)
  END IF

  READ(iunit,*,IOSTAT=ierr) ap%zv
  READ(iunit,*,IOSTAT=ierr) ap%mesh, ap%nchan

  IF((ap%nchan > SIZE(ap%vnl,2) ) .OR. (ap%nchan < 1)) THEN
    ierr = 1
    err_msg = ' NCHAN NOT PROGRAMMED '
    GO TO 110
  END IF
  IF((ap%mesh > SIZE(ap%rw) ) .OR. (ap%mesh < 0)) THEN
    info = 2
    err_msg = ' NPOTMESH OUT OF RANGE '
    GO TO 110
  END IF

  ap%rw = 0.0d0
  ap%vnl = 0.0d0
  ap%vloc = 0.0d0
  ap%vrps = 0.0d0
  DO j = 1, ap%mesh
    READ(iunit,*,IOSTAT=ierr) ap%rw(j), (ap%vnl(j,l),l=1,ap%nchan)
  END DO

  IF( MINVAL( ap%rw(1:ap%mesh) ) <= 0.0d0 ) THEN
    info = 30
    err_msg = ' ap rw too small '
    GO TO 110
  END IF

! ...  mixed reference potential is in vr(lloc)
  IF(ap%tmix) THEN
    DO j=1,ap%mesh
      ap%vnl(j,ap%lloc)= 0.d0
      DO l=1,ap%nchan
        IF(l /= ap%lloc) THEN
          ap%vnl(j,ap%lloc)=  ap%vnl(j,ap%lloc) + ap%wgv(l) * ap%vnl(j,l)
        END IF
      END DO
    END DO
  END IF
  ap%vloc(:) = ap%vnl(:,ap%lloc)
  ap%dx = calculate_dx( ap%rw, ap%mesh )
  ap%rab  = ap%dx * ap%rw

  CALL read_atomic_wf( iunit, ap, err_msg, ierr)
  IF( ierr /= 0 ) GO TO 110

  DO l = 1, ap%lnl
    ll=ap%indl(l) 
    ap%vrps(:,l) = ( ap%vnl(:,ll) - ap%vloc(:) ) * ap%rps(:,ll)
  END DO

  IF(ap%tnlcc) THEN
    CALL read_atomic_cc( iunit, ap,  err_msg, ierr)
    IF( ierr /= 0 ) GO TO 110
  END IF

  GOTO 110
100 ierr = 1
110 CONTINUE
  
  RETURN
END SUBROUTINE

!

SUBROUTINE read_head_pp( iunit, ap, err_msg, ierr)
  USE pseudo_types, ONLY: pseudo_ncpp
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: iunit
  TYPE (pseudo_ncpp), INTENT(INOUT) :: ap
  CHARACTER(LEN=*) :: err_msg
  INTEGER, INTENT(OUT) :: ierr
!
  INTEGER :: i, l

! ... read pseudo header

  ierr = 0
  err_msg = ' error while reading header pseudo '

  ap%indl = 0
  READ(iunit, *) ap%tnlcc, ap%tmix
  READ(iunit, *) ap%pottyp, ap%lloc, ap%lnl, (ap%indl(l), l = 1, MIN(ap%lnl, SIZE(ap%indl)) ) 

  IF( ap%lnl > SIZE(ap%indl) .OR. ap%lnl < 0 ) THEN
    ierr = 1
    err_msg = 'LNL out of range'
    GO TO 110
  END IF
  IF( ap%lloc < 0 .OR. ap%lloc > SIZE(ap%vnl,2) ) THEN
    ierr = 3
    err_msg = 'LLOC out of range'
    GO TO 110
  END IF
  IF( ap%tmix .AND. ap%pottyp /= 'NUMERIC' ) THEN
    ierr = 4
    err_msg = 'tmix not implemented for pseudo ' // ap%pottyp
    GO TO 110
  END IF
  DO l = 2, ap%lnl
    IF( ap%indl(l) <= ap%indl(l-1)) THEN
      ierr = 5
      err_msg =' NONLOCAL COMPONENTS MUST BE GIVEN IN ASCENDING ORDER'
      GO TO 110
    END IF
  END DO
  DO l = 1, ap%lnl
    IF( ap%indl(l) == ap%lloc) THEN
      ierr = 6
      err_msg = ' LLOC.EQ.L NON LOCAL!!' 
      GO TO 110
    END IF
  END DO

  GOTO 110
100 ierr = 1
110 CONTINUE
  
  RETURN
END SUBROUTINE

!=----------------------------------------------------------------------------=!

SUBROUTINE read_analytic_pp( iunit, ap, err_msg, ierr)
  USE pseudo_types, ONLY: pseudo_ncpp
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: iunit
  TYPE (pseudo_ncpp), INTENT(INOUT) :: ap
  CHARACTER(LEN=*) :: err_msg
  INTEGER, INTENT(OUT) :: ierr
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
      GO TO 110
  END SELECT

  DO l=1,3
    DO i=1,ap%igau
      READ(iunit,*,IOSTAT=ierr) ap%rcl(i,l), ap%al(i,l), ap%bl(i,l)
    END DO
  END DO

  CALL read_atomic_wf( iunit, ap, err_msg, ierr)
  IF( ierr /= 0 ) GO TO 110

  IF(ap%tnlcc) THEN
    CALL read_atomic_cc( iunit, ap, err_msg, ierr)
    IF( ierr /= 0 ) GO TO 110
  END IF

! ... Analytic pseudo are not supported anymore, conversion
! ... to numeric form is forced
  CALL analytic_to_numeric( ap )

  GOTO 110
100 ierr = 1
110 CONTINUE
  
  RETURN
END SUBROUTINE

!=----------------------------------------------------------------------------=!


SUBROUTINE read_atomic_cc( iunit, ap, err_msg, ierr)
  USE pseudo_types, ONLY: pseudo_ncpp
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: iunit
  TYPE (pseudo_ncpp), INTENT(INOUT) :: ap
  CHARACTER(LEN=*) :: err_msg
  INTEGER, INTENT(OUT) :: ierr
!
  CHARACTER(LEN=80) :: input_line
  INTEGER :: j, mesh
  REAL(dbl) :: rdum

! ... read atomic core

  ierr = 0
  err_msg = ' error while reading atomic core pseudo '

  ap%rhoc = 0.0d0

  READ(iunit,*,IOSTAT=ierr) mesh
  IF(mesh > SIZE(ap%rw) .OR. mesh < 0 ) THEN
    ierr = 17
    err_msg = '  CORE CORRECTION MESH OUT OF RANGE '
    GO TO 110
  END IF
  DO j = 1, mesh
    READ(iunit,*,IOSTAT=ierr) rdum, ap%rhoc(j)
    IF( ap%mesh == 0 ) ap%rw(j) = rdum
    IF( ABS(rdum - ap%rw(j))/(rdum+ap%rw(j)) > TOLMESH ) THEN
      ierr = 5
      err_msg = ' core cor. radial mesh does not match '
      GO TO 110
    END IF
  END DO

  IF( ap%mesh == 0 ) THEN
    ap%mesh = mesh
    ap%dx = calculate_dx( ap%rw, ap%mesh )
    ap%rab  = ap%dx * ap%rw
  END IF

  GOTO 110
100 ierr = 1
110 CONTINUE
  
  RETURN
END SUBROUTINE


end module




program fpmd2upf

  !
  !     Convert a pseudopotential written in the FPMD format
  !     to unified pseudopotential format
  !

  USE kinds
  USE fpmd2upf_module, ONLY: read_pseudo_fpmd, calculate_dx
  USE pseudo_types, ONLY: pseudo_ncpp
  USE parameters
  USE upf

  IMPLICIT NONE

  TYPE (pseudo_ncpp) :: ap
  CHARACTER(LEN=80) :: psfile
  INTEGER :: nsp, nspnl, i, lloc, l, ir, iv, kkbeta
  REAL(kind=8) :: rmax = 10
  REAL(kind=8) :: vll
  REAL(kind=8), allocatable :: aux(:)

! ... end of declarations
!  ----------------------------------------------

      WRITE(6,*) ' Enter the pseudopotential filename in FPMD format : '
      READ(5,'(a)') psfile

      nsp = 1
      CALL read_pseudo_fpmd(ap, psfile)

  write(generated, '("Generated using unknown code")')
  write(date_author,'("Author: unknown    Generation date: as well")')
  comment = 'Info: automatically converted from CPMD format'

  rcloc = 0.0

  write(6, * ) 'Number of wavefunction  > '
  read (5,*) nwfs

  allocate( els(nwfs), oc(nwfs), epseu(nwfs) )
  allocate( lchi(nwfs), nns(nwfs) )
  allocate( rcut (nwfs), rcutus (nwfs) )

  do i = 1, nwfs
     write(6, * ) 'Wavefunction ',i,' label, occupation > '
     read (5,*) els(i), oc(i)
     lchi(i)  = i - 1
     nns (i)  = 0
     rcut(i)  = 0.0
     rcutus(i)= 0.0
     epseu(i) = 0.0
  end do

  psd        = 'XX'
  pseudotype = 'NC'
  nlcc       = ap%tnlcc
  zp         = ap%zv
  etotps     = 0.0

  lloc  = ap%lloc 
  lmax  = MAX( MAXVAL( ap%indl( 1:ap%lnl ) ) - 1, ap%lloc - 1 )

!  go to 100

  nbeta = ap%lnl
  mesh  = ap%mesh
  ntwfc = nwfs
  allocate( elsw(ntwfc), ocw(ntwfc), lchiw(ntwfc) )
  do i = 1, nwfs
     lchiw(i) = lchi(i)
     ocw(i)   = oc(i)
     elsw(i)  = els(i)
  end do

  iexch =  1
  icorr =  1
  igcx  =  1
  igcc  =  1
  !
  ! We have igcc=2 (PW91) and 3 (LYP) exchanged wrt CPMD conventions
  !
  if (igcc.eq.3) then
     igcc=2
  else if (igcc.eq.3) then
     igcc=2
  end if

  allocate(rab(mesh))
  allocate(  r(mesh))
  r   = ap%rw
  ap%dx = calculate_dx( ap%rw, ap%mesh )
  rab = ap%rw * ap%dx

  write(6,*) ap%lloc, ap%indl( 1:ap%lnl ) , ap%lnl, ap%dx


  allocate (rho_atc(mesh))
  if (nlcc) rho_atc = ap%rhoc

  allocate (vloc0(mesh))
  ! the factor 2 converts from Hartree to Rydberg
  vloc0(:) = ap%vloc * 2.0d0

  if (nbeta > 0) then

     allocate(ikk2(nbeta), lll(nbeta))
     kkbeta = mesh
     do ir = 1,mesh
        if ( r(ir) > rmax ) then
           kkbeta=ir
           exit
        end if
     end do
     ikk2(:) = kkbeta
     allocate(aux(kkbeta))
     allocate(betar(mesh,nbeta))
     allocate(qfunc(mesh,nbeta,nbeta))
     allocate(dion(nbeta,nbeta))
     allocate(qqq (nbeta,nbeta))
     qfunc(:,:,:)=0.0d0
     dion(:,:) =0.d0
     qqq(:,:)  =0.d0
     iv = 0
     do i = 1, nwfs
        l = lchi(i)
        if ( l .ne. (lloc-1) ) then
           iv = iv + 1
           lll( iv ) = l
           do ir = 1, kkbeta
              ! the factor 2 converts from Hartree to Rydberg
              betar(ir, iv) = 2.d0 * ap%vrps( ir, iv )
              aux(ir) = ap%rps(ir, (l+1) ) * betar(ir, iv)
           end do
           call simpson2(kkbeta, aux(1), rab(1), vll)
           dion(iv,iv) = 1.0d0/vll
           write(6,*) aux(2), rab(2), kkbeta, vll
        end if
     enddo

  end if

  allocate (rho_at(mesh))
  rho_at = 0.d0
  do i = 1, nwfs
     rho_at(:) = rho_at(:) + ocw(i) * ap%rps(:, i) ** 2
  end do

  allocate (chi(mesh,ntwfc))
  chi = ap%rps
  !     ----------------------------------------------------------
  write (6,'(a)') 'Pseudopotential successfully converted'
  !     ----------------------------------------------------------

  call write_upf( 10 )

100 continue

end program


subroutine write_upf(ounps)

  use upf, only: nlcc

  integer :: ounps

  call write_pseudo_comment(ounps)  
  call write_pseudo_header(ounps)  
  call write_pseudo_mesh(ounps)
  if (nlcc)  call write_pseudo_nlcc(ounps)  
  call write_pseudo_local(ounps)  
  call write_pseudo_nl(ounps)  
  call write_pseudo_pswfc(ounps)
  call write_pseudo_rhoatom(ounps)  
  !
  print '("*** PLEASE TEST BEFORE USING!!! ***")'
  print '("review the content of the PP_INFO fields")'
  !
end subroutine write_upf

  !
  !---------------------------------------------------------------------
  subroutine write_pseudo_comment (ounps)  
    !---------------------------------------------------------------------
    !
    !
    !     This routine writes the comments of the new UPF file
    !
    use upf
    implicit none
    integer :: ounps  

    integer :: nb, ios  

    write (ounps, '(a9)', err = 100, iostat = ios) "<PP_INFO>"  

    write (ounps, '(a)', err = 100, iostat = ios) generated
    write (ounps, '(a)', err = 100, iostat = ios) date_author
    write (ounps, '(a)', err = 100, iostat = ios) comment
    if (rel==2) then  
       write (ounps, '(i5,t14,a)', err = 100, iostat = ios) rel,& 
            &"The Pseudo was generated with a Full-Relativistic Calculation"
    else if (rel==1) then  
       write (ounps, '(i5,t14,a)', err = 100, iostat = ios) rel,& 
            &"The Pseudo was generated with a Scalar-Relativistic Calculation"
    else  
       write (ounps, '(i5,t14,a)', err = 100, iostat = ios) rel, &
            & "The Pseudo was generated with a Non-Relativistic Calculation"
    endif

    write (ounps, '(1pe19.11,t24,a)', err = 100, iostat = ios) &
         rcloc, "Local Potential cutoff radius"

    write (ounps, '(a2,2a3,a6,3a19)', err = 100, iostat = ios) "nl", &
         &" pn", "l", "occ", "Rcut", "Rcut US", "E pseu"
    do nb = 1, nwfs  
       write (ounps, '(a2,2i3,f6.2,3f19.11)') els (nb) , nns (nb) , &
            lchi (nb) , oc (nb) , rcut (nb) , rcutus (nb) , epseu(nb)

    enddo

    write (ounps, '(a10)', err = 100, iostat = ios) "</PP_INFO>"  
    return
100 call error ('write_pseudo_comment', 'Writing pseudo file', abs ( &
         ios))   
  end subroutine write_pseudo_comment

  !
  !---------------------------------------------------------------------
  subroutine write_pseudo_header (ounps)  
    !---------------------------------------------------------------------
    !
    !
    !     This routine writes the header of the new UPF file
    !
    use upf
    implicit none
    integer :: ounps  
    !
    character (len=4) :: shortname
    character (len=20):: dft  
    integer :: nb, ios  
    !
    !
    write (ounps, '(//a11)', err = 100, iostat = ios) "<PP_HEADER>"  

    write (ounps, '(t3,i2,t24,a)', err = 100, iostat = ios) nv, &
         "Version Number"
    write (ounps, '(t3,a,t24,a)', err = 100, iostat = ios) psd , &
         "Element"
    if (pseudotype == 'NC') then  
       write (ounps, '(a5,t24,a)', err = 100, iostat = ios) "NC", &
            "Norm - Conserving pseudopotential"
    else if (pseudotype == 'US') then
       write (ounps, '(a5,t24,a)', err = 100, iostat = ios) "US", &
            "Ultrasoft pseudopotential"
    else
       call error ('write_pseudo_header',&
            'Unknown PP type: '//pseudotype, 1)
    endif
    write (ounps, '(l5,t24,a)', err = 100, iostat = ios) nlcc , &
         "Nonlinear Core Correction"
    call dftname (iexch, icorr, igcx, igcc, dft, shortname)
    write (ounps, '(a,t24,a4,a)', err = 100, iostat = ios) &
         dft, shortname," Exchange-Correlation functional"
    write (ounps, '(f17.11,t24,a)') zp , "Z valence"  
    write (ounps, '(f17.11,t24,a)') etotps, "Total energy"  
    write (ounps, '(2f11.7,t24,a)') ecutrho, ecutwfc, &
         "Suggested cutoff for wfc and rho"  

    write (ounps, '(i5,t24,a)') lmax, "Max angular momentum component"  
    write (ounps, '(i5,t24,a)') mesh, "Number of points in mesh"
    write (ounps, '(2i5,t24,a)', err = 100, iostat = ios) ntwfc, &
         nbeta  , "Number of Wavefunctions, Number of Projectors"
    write (ounps, '(a,t24,a2,a3,a6)', err = 100, iostat = ios) &
         " Wavefunctions", "nl", "l", "occ"
    do nb = 1, ntwfc
       write (ounps, '(t24,a2,i3,f6.2)') elsw(nb), lchiw(nb), ocw(nb)
    enddo
    !---> End header writing

    write (ounps, '(a12)', err = 100, iostat = ios) "</PP_HEADER>"
    return   
100 call error ('write_pseudo_header','Writing pseudo file', abs(ios) )

  end subroutine write_pseudo_header

  !
  !---------------------------------------------------------------------
  subroutine write_pseudo_mesh (ounps)  
    !---------------------------------------------------------------------
    !
    !
    !     This routine writes the atomic charge density to the new UPF file
    !
    use upf
    implicit none
    integer :: ounps  
    !
    integer :: ir, ios  
    !
    write (ounps, '(//a9)', err = 100, iostat = ios) "<PP_MESH>"  

    write (ounps, '(t3,a6)', err = 100, iostat = ios) "<PP_R>"  
    write (ounps, '(1p4e19.11)', err=100, iostat=ios) (r(ir),  ir=1,mesh )
    write (ounps, '(t3,a7)', err = 100, iostat = ios) "</PP_R>"  
    write (ounps, '(t3,a8)', err = 100, iostat = ios) "<PP_RAB>"  
    write (ounps, '(1p4e19.11)', err=100, iostat=ios) (rab(ir), ir=1,mesh )
    write (ounps, '(t3,a9)', err = 100, iostat = ios) "</PP_RAB>"  

    write (ounps, '(a10)', err = 100, iostat = ios) "</PP_MESH>"  

    return

100 call error ('write_pseudo_rhoatom','Writing pseudo file',abs(ios))

  end subroutine write_pseudo_mesh

  !
  !---------------------------------------------------------------------
  subroutine write_pseudo_nlcc (ounps)  
    !---------------------------------------------------------------------
    !
    !
    !     This routine writes the core charge for the nonlinear core
    !     correction of the new UPF file
    !
    use upf
    implicit none
    integer :: ounps  
    !
    integer :: ir, ios  

    write (ounps, '(//a9)', err = 100, iostat = ios) "<PP_NLCC>"  

    write (ounps, '(1p4e19.11)', err=100, iostat=ios) &
                                 ( rho_atc(ir), ir = 1, mesh )
    write (ounps, '(a10)', err = 100, iostat = ios) "</PP_NLCC>"  
    return

100 call error ('write_pseudo_nlcc', 'Writing pseudo file', abs (ios))

  end subroutine write_pseudo_nlcc
  !
  !---------------------------------------------------------------------
  subroutine write_pseudo_local (ounps)  
    !---------------------------------------------------------------------
    !
    !
    !     This routine writes the local part of the new UPF file
    !
    use upf
    implicit none
    integer :: ounps  
    !
    integer :: ir, ios  

    write (ounps, '(//a10)', err = 100, iostat = ios) "<PP_LOCAL>"  
    write (ounps, '(1p4e19.11)', err=100, iostat=ios) &
                                ( vloc0(ir), ir = 1, mesh )
    write (ounps, '(a11)', err = 100, iostat = ios) "</PP_LOCAL>"  
    return
100 call error ('write_pseudo_local', 'Writing pseudo file', abs(ios) )  
  end subroutine write_pseudo_local

  !
  !---------------------------------------------------------------------
  subroutine write_pseudo_nl (ounps)  
    !---------------------------------------------------------------------
    !
    !
    !     This routine writes the non local part of the new UPF file
    !
    use upf
    implicit none
    integer :: ounps  
    !
    integer :: nb, mb, n, ir, nd, i, lp, ios  

    write (ounps, '(//a13)', err = 100, iostat = ios) "<PP_NONLOCAL>"  
    do nb = 1, nbeta  
       write (ounps, '(t3,a9)', err = 100, iostat = ios) "<PP_BETA>"  
       write (ounps, '(2i5,t24,a)', err=100, iostat=ios) &
                                    nb, lll(nb), "Beta    L"
       write (ounps, '(i6)', err=100, iostat=ios) ikk2 (nb)  
       write (ounps, '(1p4e19.11)', err=100, iostat=ios) &
                                    ( betar(ir,nb), ir=1,ikk2(nb) )
       write (ounps, '(t3,a10)', err = 100, iostat = ios) "</PP_BETA>"  
    enddo

    write (ounps, '(t3,a8)', err = 100, iostat = ios) "<PP_DIJ>"  
    nd = 0  
    do nb = 1, nbeta  
       do mb = nb, nbeta  
          if ( abs(dion(nb,mb)) .gt. 1.0d-12 )  nd = nd + 1 
       enddo
    enddo
    write (ounps, '(1p,i5,t24,a)', err=100, iostat=ios) &
                                   nd, "Number of nonzero Dij"
    do nb = 1, nbeta
       do mb = nb, nbeta  
          if ( abs(dion(nb,mb)) .gt. 1.0d-12 ) &
             write(ounps,'(1p,2i5,e19.11)', err=100, iostat=ios) &
                                   nb, mb, dion(nb,mb)
       enddo
    enddo
    write (ounps, '(t3,a9)', err=100, iostat=ios) "</PP_DIJ>"  

    if (pseudotype == 'US') then  
       write (ounps, '(t3,a8)', err = 100, iostat = ios) "<PP_QIJ>"  
       write (ounps, '(i5,a)',err=100, iostat=ios) nqf,"     nqf.&
          & If not zero, Qij's inside rinner are computed using qfcoef's"
       if (nqf.gt.0) then
          write (ounps, '(t5,a11)', err=100, iostat=ios) "<PP_RINNER>"  
          write (ounps,'(i5,1pe19.11)', err=100, iostat=ios) &
                                        (i, rinner(i), i = 1, nqlc)
          write (ounps, '(t5,a12)', err=100, iostat=ios) "</PP_RINNER>"  
       end if
       do nb = 1, nbeta 
          do mb = nb, nbeta
             write (ounps, '(3i5,t24,a)', err=100, iostat=ios) &
                                          nb, mb, lll(mb) , "i  j  (l(j))"
             write (ounps, '(1pe19.11,t24,a)', err=100, iostat=ios) &
                                          qqq(nb,mb), "Q_int"
             write (ounps, '(1p4e19.11)', err=100, iostat=ios) &
                                          ( qfunc (n,nb,mb), n=1,mesh )
             if (nqf.gt.0) then
                write (ounps, '(t5,a11)', err=100, iostat=ios) &
                                          "<PP_QFCOEF>"  
                write(ounps,'(1p4e19.11)', err=100, iostat=ios) &
                                 ((qfcoef(i,lp,nb,mb),i=1,nqf),lp=1,nqlc)
                write (ounps, '(t5,a12)', err=100, iostat=ios) &
                                          "</PP_QFCOEF>"
             end if
          enddo
       enddo
       write (ounps, '(t3,a9)', err = 100, iostat = ios) "</PP_QIJ>"  

    endif
    write (ounps, '(a14)', err = 100, iostat = ios) "</PP_NONLOCAL>"  
    return

100 call error ('write_pseudo_nl', 'Writing pseudo file', abs (ios) )  

  end subroutine write_pseudo_nl

  !
  !---------------------------------------------------------------------
  subroutine write_pseudo_pswfc (ounps)  
    !---------------------------------------------------------------------
    !
    !
    !     This routine writes the pseudo atomic functions
    !     of the new UPF file
    !
    use upf
    implicit none
    integer :: ounps  
    !
    integer :: nb, ir, ios  

    write (ounps, '(//a10)', err = 100, iostat = ios) "<PP_PSWFC>"  
    do nb = 1, ntwfc
       write (ounps,'(a2,i5,f6.2,t24,a)', err=100, iostat=ios) &
            elsw(nb), lchiw(nb), ocw(nb), "Wavefunction"
       write (ounps, '(1p4e19.11)', err=100, iostat=ios) &
            ( chi(ir,nb), ir=1,mesh )
    enddo
    write (ounps, '(a11)', err = 100, iostat = ios) "</PP_PSWFC>"  
    return

100 call error ('write_pseudo_pswfc', 'Writing pseudo file', abs(ios) )  
  end subroutine write_pseudo_pswfc
  !
  !---------------------------------------------------------------------
  subroutine write_pseudo_rhoatom (ounps)  
    !---------------------------------------------------------------------
    !
    !
    !     This routine writes the atomic charge density to the new UPF file
    !
    use upf
    implicit none
    integer :: ounps  
    !
    integer :: ir, ios  

    write (ounps, '(//a12)', err = 100, iostat = ios) "<PP_RHOATOM>"  
    write (ounps, '(1p4e19.11)', err = 100, iostat = ios) &
                               ( rho_at(ir), ir=1,mesh )
    write (ounps, '(a13)', err = 100, iostat = ios) "</PP_RHOATOM>"  
    return

100 call error('write_pseudo_rhoatom','Writing pseudo file',abs(ios))
  end subroutine write_pseudo_rhoatom

  !---------------------------------------------------------------------
  subroutine dftname(iexch, icorr, igcx, igcc, longname, shortname)
  !---------------------------------------------------------------------
  implicit none
  integer iexch, icorr, igcx, igcc
  character (len=4) :: shortname
  character (len=20):: longname
  !
  ! The data used to convert iexch, icorr, igcx, igcc
  ! into a user-readable string
  !
  integer, parameter :: nxc = 1, ncc = 9, ngcx = 3, ngcc = 4 
  character (len=4) :: exc, corr, gradx, gradc  
  dimension exc (0:nxc), corr (0:ncc), gradx (0:ngcx), gradc (0:ngcc)
  data exc / 'NOX ', 'SLA ' /  
  data corr / 'NOC ', 'PZ  ', 'VWN ', 'LYP ', 'PW  ', 'WIG ', 'HL  ',&
              'OBZ ', 'OBW ', 'GL  ' /
  data gradx / 'NOGX', 'B88 ', 'GGX ', 'PBE ' /  
  data gradc / 'NOGC', 'P86 ', 'GGC ', 'BLYP', 'PBE ' /  

  if (iexch==1.and.igcx==0.and.igcc==0) then
     shortname = corr(icorr)
  else if (iexch==1.and.icorr==3.and.igcx==1.and.igcc==3) then
     shortname = 'BLYP'
  else if (iexch==1.and.icorr==1.and.igcx==1.and.igcc==0) then
     shortname = 'B88'
  else if (iexch==1.and.icorr==1.and.igcx==1.and.igcc==1) then
     shortname = 'BP'
  else if (iexch==1.and.icorr==4.and.igcx==2.and.igcc==2) then
     shortname = 'PW91'
  else if (iexch==1.and.icorr==4.and.igcx==3.and.igcc==4) then
     shortname = 'PBE'
  else
     shortname = ' '
  end if
  write(longname,'(4a5)') exc(iexch),corr(icorr),gradx(igcx),gradc(igcc)

  return
end subroutine dftname

subroutine error(a,b,n)
  character(len=*) :: a,b

  write(6,'(//'' program '',a,'':'',a,''.'',8x,i8,8x,''stop'')') a,b,n
  stop
end subroutine error

!----------------------------------------------------------------------
subroutine simpson2(mesh,func,rab,asum)
  !-----------------------------------------------------------------------
!
  !     simpson's rule integrator for function stored on the
  !     radial logarithmic mesh
  !

  implicit none

  integer :: i, mesh
  real(kind=8) ::  rab(mesh), func(mesh), f1, f2, f3, r12, asum

      !     routine assumes that mesh is an odd number so run check
      !     if ( mesh+1 - ( (mesh+1) / 2 ) * 2 .ne. 1 ) then
      !       write(*,*) '***error in subroutine radlg'
!       write(*,*) 'routine assumes mesh is odd but mesh =',mesh+1
!       stop
!     endif

  asum = 0.0d0
  r12 = 1.0d0 / 12.0d0
  f3  = func(1) * rab(1) * r12

  do i = 2,mesh-1,2
     f1 = f3
     f2 = func(i) * rab(i) * r12
     f3 = func(i+1) * rab(i+1) * r12
     asum = asum + 4.0d0*f1 + 16.0d0*f2 + 4.0d0*f3
  enddo

  return
end subroutine simpson2

