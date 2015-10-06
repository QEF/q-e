!
! Copyright (C) 2001-2011 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file 'License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!---------------------------------------------------------------------
PROGRAM cpmd2upf
  !---------------------------------------------------------------------
  !
  !     Convert a pseudopotential written in the CPMD format, TYPE:
  !     NORMCONSERVING [ NUMERIC, LOGARITHMIC, CAR ], single radial grid
  !     to unified pseudopotential format (v.2)
  !
  USE pseudo_types, ONLY : pseudo_upf, nullify_pseudo_upf, &
                           deallocate_pseudo_upf
  USE write_upf_v2_module, ONLY :  write_upf_v2
  !
  IMPLICIT NONE
  TYPE(pseudo_upf) :: upf
  CHARACTER(len=256) filein, fileout
  INTEGER :: ios
  !
  IF ( trim(filein) == ' ') &
       CALL errore ('cpmd2upf', 'usage: cpmd2upf "file-to-be-converted"', 1)
  CALL get_file ( filein )
  OPEN ( unit=1, file=filein, status = 'old', form='formatted', iostat=ios )
  IF ( ios /= 0) CALL errore ('cpmd2upf', 'file: '//trim(filein)//' not found', 2)
  !
  CALL read_cpmd(1)
  !
  CLOSE (unit=1)
  !
  ! convert variables read from CPMD format into those needed
  ! by the upf format - add missing quantities
  !
  CALL nullify_pseudo_upf ( upf )
  !
  CALL convert_cpmd (upf)
  !
  ! write to file
  !
  fileout=trim(filein)//'.UPF'
  PRINT '(''Output PP file in UPF format :  '',a)', fileout
  OPEN(unit=2,file=fileout,status='unknown',form='formatted')
  !
  CALL write_upf_v2 (2, upf )
  !
  CLOSE (unit=2)
  CALL deallocate_pseudo_upf ( upf )
  !     ----------------------------------------------------------
  WRITE (6,"('Pseudopotential successfully written')")
  WRITE (6,"('Please review the content of the PP_INFO fields')")
  WRITE (6,"('*** Please TEST BEFORE USING !!! ***')")
  !     ----------------------------------------------------------
  !
  STOP

END PROGRAM cpmd2upf

MODULE cpmd
  !
  ! All variables read from CPMD file format
  !
  CHARACTER (len=80) title
  !
  INTEGER :: ixc, pstype = 0
  real(8) :: alphaxc
  REAL(8) :: z, zv
  ! Grid variables
  INTEGER :: mesh
  real(8) :: amesh, rmax, xmin
  real(8), ALLOCATABLE :: r(:)
  ! PP variables
  INTEGER, PARAMETER :: lmaxx=3
  INTEGER ::lmax, nwfc=0
  ! Car PP variables
  real(8) :: alphaloc, alpha(0:lmaxx), a(0:lmaxx), b(0:lmaxx)
  ! Goedecker PP variables
  INTEGER, PARAMETER :: ncmax=4, nlmax=3
  INTEGER :: nc, nl(0:lmaxx)
  real(8) :: rc, rl(0:lmaxx), c(ncmax), h(0:lmaxx, nlmax*(nlmax+1)/2 )
  ! Numeric PP variables
  real(8), ALLOCATABLE :: vnl(:,:)
  real(8), ALLOCATABLE :: chi(:,:)
  ! Core correction variables
  LOGICAL :: nlcc=.false.
  real(8), ALLOCATABLE :: rho_atc(:)
  ! Variables used for reading and for checks
  INTEGER :: maxinfo, info_lines
  PARAMETER (maxinfo = 100)
  CHARACTER (len=80), ALLOCATABLE :: info_sect(:)
  !------------------------------

END MODULE cpmd
!
!     ----------------------------------------------------------
SUBROUTINE read_cpmd(iunps)
  !     ----------------------------------------------------------
  !
  USE cpmd
  IMPLICIT NONE
  INTEGER :: iunps
  !
  INTEGER :: found = 0, closed = 0, unknown = 0
  INTEGER :: i, l, dum, ios
  CHARACTER (len=256) line
  CHARACTER (len=4) token
  real (8) :: amesh_, vnl0(0:3)
  LOGICAL :: grid_read = .false., wfc_read=.false.
  LOGICAL, EXTERNAL :: matches
  INTEGER, EXTERNAL :: locate
  REAL(8), EXTERNAL :: qe_erf
  !
  info_lines = 0
10 READ (iunps,'(A)',end=20,err=20) line
  IF (matches ("&ATOM", trim(line)) ) THEN
     found = found + 1
     ! Z
     READ (iunps,'(a)',end=200,err=200) line
     l = len_trim(line)
     i = locate('=',line)
     READ (line(i+1:l),*) z
     ! Zv
     READ (iunps,'(a)',end=200,err=200) line
     l = len_trim(line)
     i = locate('=',line)
     READ (line(i+1:l),*) zv
     ! XC
     READ (iunps,'(a)',end=200,err=200) line
     l = len_trim(line)
     i = locate('=',line)
     READ (line(i+1:l),*) ixc, alphaxc
     ! TYPE
     READ (iunps,'(a)',end=200,err=200) line
     IF ( matches("NORMCONSERVING",line) ) THEN
        IF ( matches("NUMERIC",line) .or. matches("LOGARITHMIC",line)  ) THEN
           pstype = 1
        ELSEIF ( matches("CAR",line) ) THEN
           pstype = 2
        ELSEIF ( matches("GOEDECKER",line) ) THEN
           pstype = 3
        ENDIF
     ENDIF
     IF (pstype == 0 ) CALL errore('read_cpmd','unknown type: '//line,1)
  ELSEIF (matches ("&INFO", trim(line)) ) THEN
     found = found + 1
     ! read (iunps,'(a)') title
     ! store info section for later perusal
     ALLOCATE (info_sect(maxinfo))
     DO i=1,maxinfo
        READ (iunps,'(a)',end=20,err=20) title
        IF (matches ("&END", trim(title)) ) THEN
           closed = closed + 1
           GOTO 10
        ELSE
           info_sect(i) = trim(title)
           info_lines = i
        ENDIF
     ENDDO
  ELSEIF (matches ("&POTENTIAL", trim(line)) ) THEN
     found = found + 1
     READ (iunps,'(a)') line
     IF ( pstype == 1 ) THEN
        !
        ! NORMCONSERVING NUMERIC
        !
        READ (line,*,iostat=ios) mesh, amesh_
        IF ( ios /= 0) THEN
           READ (line,*,iostat=ios) mesh
           amesh_ = -1.0d0
        ENDIF
        IF ( .not. grid_read ) ALLOCATE (r(mesh))
        !
        ! determine the number of angular momenta
        !
        READ (iunps, '(a)') line
        ios = 1
        lmax= 4
        DO WHILE (ios /= 0)
           lmax = lmax - 1
           READ(line,*,iostat=ios) r(1),(vnl0(l),l=0,lmax)
        ENDDO
        ALLOCATE (vnl(mesh,0:lmax))
        vnl(1,0:lmax) = vnl0(0:lmax)
        DO i=2,mesh
           READ(iunps, *) r(i),(vnl(i,l),l=0,lmax)
        ENDDO
        IF ( .not.grid_read ) THEN
           CALL check_radial_grid ( amesh_, mesh, r, amesh )
           grid_read = .true.
        ENDIF
     ELSEIF ( pstype == 2 ) THEN
     !
     ! NORMCONSERVING CAR
     !
        READ(iunps, *) alphaloc
        ! convert r_c's written in file to alpha's: alpha = 1/r_c^2
        alphaloc = 1.d0/alphaloc**2
        DO lmax=-1,2
           READ(iunps, '(A)') line
           IF (matches ("&END", trim(line)) ) THEN
              closed = closed + 1
              exit
           ENDIF
           READ(line, *) alpha(lmax+1), a(lmax+1), b(lmax+1)
           alpha(lmax+1) = 1.d0/alpha(lmax+1)**2
        ENDDO
     ELSEIF ( pstype == 3 ) THEN
     !
     ! NORMCONSERVING GOEDECKER
     !
        c(:) = 0.d0
        rl(:) = 0.d0
        nl(:) = 0
        h(:,:) = 0.d0
        READ(iunps, *) lmax
        lmax = lmax - 1
        IF ( lmax > lmaxx ) &
          CALL errore('read_cpmd','incorrect parameter read',1)
        READ(iunps, *) rc
        READ(iunps, '(A)') line
        READ(line, *) nc
        IF ( nc > ncmax ) &
          CALL errore('read_cpmd','incorrect parameter read',2)
        ! I am not sure if it is possible to use nc in the same line
        ! where it is read. Just in case, better to read twice
        READ(line, *) dum, (c(i), i=1,nc)
        DO l=0,lmax+1
           READ(iunps, '(A)') line
           IF ( matches ("&END", trim(line)) ) THEN
              closed = closed + 1
              exit
           ENDIF
           READ(line, *) rl(l), nl(l)
           IF ( nl(l) > nlmax ) &
             CALL errore('read_cpmd','incorrect parameter read',3)
           IF ( nl(l) > 0 ) &
              READ(line, *) rl(l), dum, ( h(l,i), i=1,nl(l)*(nl(l)+1)/2)
        ENDDO
     ENDIF
  ELSEIF (matches ("&WAVEFUNCTION", trim(line)) ) THEN
     wfc_read=.true.
     found = found + 1
     READ (iunps,'(A)') line
     READ (line,*,iostat=ios) mesh, amesh_
     IF ( ios /= 0) THEN
        READ (line,*,iostat=ios) mesh
        amesh_ = -1.0d0
     ENDIF
     IF ( .not. grid_read )  ALLOCATE(r(mesh))
     ! find number of atomic wavefunctions
     READ (iunps,'(A)') line
     DO nwfc = lmax+1,1,-1
        READ(line,*,iostat=ios) r(1),(vnl0(l),l=0,nwfc-1)
        IF ( ios == 0 ) exit
     ENDDO
     IF ( ios /= 0 ) &
        CALL errore('read_cpmd','at least one atomic wvfct should be present',1)
     ALLOCATE(chi(mesh,nwfc))
     chi(1,1:nwfc) = vnl0(0:nwfc-1)
     DO i=2,mesh
        READ(iunps, *) r(i),(chi(i,l),l=1,nwfc)
     ENDDO
     IF ( .not.grid_read ) THEN
        CALL check_radial_grid ( amesh_, mesh, r, amesh )
        grid_read = .true.
     ENDIF
  ELSEIF (matches ("&NLCC", trim(line)) ) THEN
     found = found + 1
     READ (iunps, '(a)') line
     READ(iunps, *) mesh
     nlcc = ( mesh > 0 )
     IF (nlcc) THEN
        IF ( .not. matches ("NUMERIC", trim(line)) ) &
          CALL errore('read_cpmd',' only NUMERIC core-correction supported',1)
        ALLOCATE (rho_atc(mesh))
        READ(iunps, * ) (r(i), rho_atc(i), i=1,mesh)
     ENDIF
  ELSEIF (matches ("&ATDENS", trim(line)) ) THEN
     ! skip over &ATDENS section, add others here, if there are more.
     DO WHILE(.not. matches("&END", trim(line)))
        READ (iunps,'(a)') line
     ENDDO
  ELSEIF (matches ("&END", trim(line)) ) THEN
     closed = closed + 1
  ELSE
     PRINT*, 'line ignored: ', line
     unknown = unknown + 1
  ENDIF
  GOTO 10

20 CONTINUE

  IF ( pstype /= 3 ) THEN
     IF (nlcc .and. found /= 5 .or. .not.nlcc .and. found /= 4) &
         CALL errore('read_cpmd','some &FIELD card missing',found)
  ELSE
     IF (found /= 3) &
         CALL errore('read_cpmd','some &FIELD card missing',found)
  ENDIF
  IF (closed /= found) &
       CALL errore('read_cpmd','some &END card missing',closed)
  IF (unknown /= 0 ) PRINT '("WARNING: ",i3," cards not read")', unknown
  !
  IF ( .not. grid_read ) THEN
     xmin = -7.0d0
     amesh=0.0125d0
     rmax =100.0d0
     PRINT '("A radial grid must be provided. We use the following one:")'
     PRINT '("r_i = e^{xmin+(i-1)*dx}/Z, i=1,mesh, with parameters:")'
     PRINT '("Z=",f6.2,", xmin=",f6.2," dx=",f8.4," rmax=",f6.1)', &
           z, xmin, amesh, rmax
     mesh = 1 + (log(z*rmax)-xmin)/amesh
     mesh = (mesh/2)*2+1 ! mesh is odd (for historical reasons?)
     ALLOCATE (r(mesh))
     DO i=1, mesh
        r(i) = exp (xmin+(i-1)*amesh)/z
     ENDDO
     PRINT '(I4," grid points, rmax=",f8.4)', mesh, r(mesh)
     grid_read = .true.
  ENDIF
  rmax = r(mesh)
  xmin = log(z*r(1))
  !
  IF ( .not. wfc_read ) PRINT '("Notice: atomic wfcs not found")'
  !
  IF ( pstype == 2 ) THEN
     ALLOCATE (vnl(mesh,0:lmax))
     DO l=0, lmax
        DO i=1, mesh
           vnl(i,l) = ( a(l) + b(l)*r(i)**2 ) * exp (-alpha(l)*r(i)**2) - &
                       zv * qe_erf (sqrt(alphaloc)*r(i))/r(i)
        ENDDO
     ENDDO
  ENDIF

  RETURN
200 CALL errore('read_cpmd','error in reading file',1)

END SUBROUTINE read_cpmd

!     ----------------------------------------------------------
SUBROUTINE convert_cpmd(upf)
  !     ----------------------------------------------------------
  !
  USE cpmd
  USE pseudo_types, ONLY : pseudo_upf
  USE constants, ONLY : e2
  !
  IMPLICIT NONE
  !
  TYPE(pseudo_upf) :: upf
  !
  REAL(8), ALLOCATABLE :: aux(:)
  REAL(8) :: x, x2, vll, rcloc, fac
  REAL(8), EXTERNAL :: mygamma, qe_erf
  CHARACTER (len=20):: dft
  CHARACTER (len=2):: label
  CHARACTER (len=1):: spdf(0:3) = ['S','P','D','F']
  CHARACTER (len=2), EXTERNAL :: atom_name
  INTEGER :: lloc, my_lmax, l, i, j, ij, ir, iv, jv
  !
  ! NOTE: many CPMD pseudopotentials created with the 'Hamann' code
  ! from Juerg Hutter's homepage have additional (bogus) entries for
  ! pseudo-potential and wavefunction. In the 'report' they have
  ! the same rc and energy eigenvalue than the previous angular momentum.
  ! we need to be able to ignore that part or the resulting UPF file
  ! will be useless. so we first print the info section and ask
  ! for the LMAX to really use. AK 2005/03/30.
  !
  DO i=1,info_lines
     PRINT '(A)', info_sect(i)
  ENDDO
  IF ( pstype == 3 ) THEN
     ! not actually used, except by write_upf to write a meaningful message
     lloc = -3
     rcloc=0.0
  ELSE
     WRITE(*,'("max L to use ( <= ",I1," ) > ")', advance="NO") lmax
     READ (5,*) my_lmax
     IF ((my_lmax <= lmax) .and. (my_lmax >= 0)) lmax = my_lmax
     WRITE(*,'("local L ( <= ",I1," ), Rc for local pot (au) > ")', advance="NO") lmax
     READ (5,*) lloc, rcloc
  ENDIF
  !
  IF ( pstype == 3 ) THEN
     upf%generated= "Generated in analytical, separable form"
     upf%author   = "Goedecker/Hartwigsen/Hutter/Teter"
     upf%date     = "Phys.Rev.B58, 3641 (1998); B54, 1703 (1996)"
  ELSE
     upf%generated= "Generated using unknown code"
     upf%author   = "unknown"
     upf%date     = "unknown"
  ENDIF
  upf%nv       = "2.0.1"
  upf%comment  = "Info: automatically converted from CPMD format"
  upf%psd      = atom_name ( nint(z) )
  ! reasonable assumption
  IF (z > 18) THEN
     upf%rel = 'no'
  ELSE
     upf%rel = 'scalar'
  ENDIF
  IF ( pstype == 3 ) THEN
     upf%typ = 'NC'
  ELSE
     upf%typ = 'SL'
  ENDIF
  upf%tvanp = .false.
  upf%tpawp = .false.
  upf%tcoulombp=.false.
  upf%nlcc = nlcc
  !
  IF (ixc==900) THEN
     PRINT '("Pade approx. not implemented! assuming Perdew-Zunger LDA")'
     upf%dft='SLA-PZ-NOGX-NOGC'
  ELSEIF (ixc==1100) THEN
     upf%dft='SLA-PZ-NOGX-NOGC'
  ELSEIF (ixc==1111) THEN
     upf%dft='SLA-PZ-B86-P88'
  ELSEIF (ixc==1134 .or. ixc==1434) THEN
     upf%dft='SLA-PW-PBX-PBC'
  ELSEIF (ixc==1134) THEN
     upf%dft='revPBE'
  ELSEIF (ixc==1197) THEN
     upf%dft='PBESOL'
  ELSEIF (ixc==1312) THEN
     upf%dft='BLYP'
  ELSEIF (ixc==362) THEN
     upf%dft='OLYP'
  ELSEIF (ixc==1372) THEN
     upf%dft='XLYP'
  ELSEIF (ixc==55) THEN
     upf%dft='HCTH'
  ELSE
     WRITE(*,'("Unknown DFT ixc=",i4,". Please provide a DFT name > ")', advance="NO") ixc
     READ *, upf%dft
  ENDIF
  PRINT '("Assuming DFT: ",A," . Please check this is what you want")', &
          trim(upf%dft)
  !
  upf%zp = zv
  upf%etotps =0.0d0
  upf%ecutrho=0.0d0
  upf%ecutwfc=0.0d0
  IF ( lmax == lloc) THEN
     upf%lmax = lmax-1
  ELSE
     upf%lmax = lmax
  ENDIF
  upf%lloc = lloc
  upf%lmax_rho = 0
  upf%nwfc = nwfc
  !
  ALLOCATE( upf%els(upf%nwfc) )
  ALLOCATE( upf%oc(upf%nwfc) )
  ALLOCATE( upf%epseu(upf%nwfc) )
  ALLOCATE( upf%lchi(upf%nwfc) )
  ALLOCATE( upf%nchi(upf%nwfc) )
  ALLOCATE( upf%rcut_chi (upf%nwfc) )
  ALLOCATE( upf%rcutus_chi(upf%nwfc) )

  DO i=1, upf%nwfc
10   WRITE(*,'("Wavefunction # ",i1,": label (e.g. 4s), occupancy > ")', advance="NO") i
     READ (5,*) label, upf%oc(i)
     READ (label(1:1),*, err=10) l
     upf%els(i)  = label
     upf%nchi(i)  = l
     IF ( label(2:2) == 's' .or. label(2:2) == 'S') THEN
        l=0
     ELSEIF ( label(2:2) == 'p' .or. label(2:2) == 'P') THEN
        l=1
     ELSEIF ( label(2:2) == 'd' .or. label(2:2) == 'D') THEN
        l=2
     ELSEIF ( label(2:2) == 'f' .or. label(2:2) == 'F') THEN
        l=3
     ELSE
        l=i-1
     ENDIF
     upf%lchi(i)  = l
     upf%rcut_chi(i)  = 0.0d0
     upf%rcutus_chi(i)= 0.0d0
     upf%epseu(i) = 0.0d0
  ENDDO

  upf%mesh = mesh
  upf%dx   = amesh
  upf%rmax = rmax
  upf%xmin = xmin
  upf%zmesh= z
  ALLOCATE(upf%rab(upf%mesh))
  ALLOCATE(upf%r(upf%mesh))
  upf%r(:) = r(1:upf%mesh)
  upf%rab(:)=upf%r(:)*amesh

  ALLOCATE (upf%rho_atc(upf%mesh))
  IF (upf%nlcc) upf%rho_atc(:) = rho_atc(1:upf%mesh)

  upf%rcloc = rcloc
  ALLOCATE (upf%vloc(upf%mesh))
  !
  ! the factor e2=2 converts from Hartree to Rydberg
  !
  IF ( upf%typ == "SL" ) THEN
     upf%vloc(:) = vnl(1:upf%mesh,upf%lloc)*e2
     ALLOCATE(upf%vnl(upf%mesh,0:upf%lmax,1))
     upf%vnl(:,:,1) = vnl(1:upf%mesh,0:upf%lmax)*e2
     upf%nbeta= lmax
  ELSE
     DO i=1,upf%mesh
        x = upf%r(i)/rc
        x2=x**2
        upf%vloc(i) = e2 * ( -upf%zp*qe_erf(x/sqrt(2.d0))/upf%r(i) + &
             exp ( -0.5d0*x2 ) * (c(1) + x2*( c(2) + x2*( c(3) + x2*c(4) ) ) ) )
     ENDDO
     upf%nbeta=0
     DO l=0,upf%lmax
        upf%nbeta = upf%nbeta + nl(l)
     ENDDO
  ENDIF

  IF (upf%nbeta > 0) THEN
     ALLOCATE(upf%els_beta(upf%nbeta) )
     ALLOCATE(upf%lll(upf%nbeta))
     ALLOCATE(upf%kbeta(upf%nbeta))
     IF ( pstype == 3 ) THEN
        upf%kbeta(:) = upf%mesh
     ELSE
        iv=0  ! counter on beta functions
        DO i=1,upf%nwfc
           l=upf%lchi(i)
           IF (l/=lloc) THEN
              iv=iv+1
              upf%kbeta(iv)=upf%mesh
              DO ir = upf%mesh,1,-1
                 IF ( abs ( vnl(ir,l) - vnl(ir,lloc) ) > 1.0E-6 ) THEN
                    ! include points up to the last with nonzero value
                    upf%kbeta(iv)=ir+1
                    exit
                 ENDIF
              ENDDO
           ENDIF
        ENDDO
        ! the number of points used in the evaluation of integrals
        ! should be even (for simpson integration)
        DO i=1,upf%nbeta
           IF ( mod (upf%kbeta(i),2) == 0 ) upf%kbeta(i)=upf%kbeta(i)+1
           upf%kbeta(i)=MIN(upf%mesh,upf%kbeta(i))
        ENDDO
        upf%kkbeta = maxval(upf%kbeta(:))
     ENDIF
     ALLOCATE(upf%beta(upf%mesh,upf%nbeta))
     ALLOCATE(upf%dion(upf%nbeta,upf%nbeta))
     upf%beta(:,:) =0.d0
     upf%dion(:,:) =0.d0
     ALLOCATE(upf%rcut  (upf%nbeta))
     ALLOCATE(upf%rcutus(upf%nbeta))

     IF ( pstype == 3 ) THEN
        iv=0  ! counter on beta functions
        DO l=0,upf%lmax
           ij = 0
           DO i=1, nl(l)
              iv = iv+1
              upf%lll(iv)=l
              WRITE (upf%els_beta(iv), '(I1,A1)' ) i, spdf(l)
              DO j=i, nl(l)
                 jv = iv+j-i
                 ij=ij+1
                 upf%dion(iv,jv) = h(l,ij)/e2
                 IF ( j > i ) upf%dion(jv,iv) = upf%dion(iv,jv)
              ENDDO
              fac= sqrt(2d0*rl(l)) / ( rl(l)**(l+2*i) * sqrt(mygamma(l+2*i)) )
              DO ir=1,upf%mesh
                 x2 = (upf%r(ir)/rl(l))**2
                 upf%beta(ir,iv) = upf%r(ir)**(l+2*(i-1)) * &
                                     exp ( -0.5d0*x2 ) * fac * e2
                 ! ...remember: the beta functions in the UPF format
                 ! ...have to be multiplied by a factor r !!!
                 upf%beta(ir,iv) = upf%beta(ir,iv)*upf%r(ir)
                 !
              ENDDO
              ! look for index kbeta such that v(i)=0 if i>kbeta
              DO ir=upf%mesh,1,-1
                 IF ( abs(upf%beta(ir,iv)) > 1.D-12 ) exit
              ENDDO
              IF ( ir < 2 ) THEN
                 CALL errore('cpmd2upf','zero beta function?!?',iv)
              ELSEIF ( mod(ir,2) /= 0 ) THEN
                 ! even index
                 upf%kbeta(iv) = ir
              ELSEIF ( ir < upf%mesh .and. mod(ir,2) == 0 ) THEN
                 ! odd index
                 upf%kbeta(iv) = ir+1
              ELSE
                 upf%kbeta(iv) = upf%mesh
              ENDIF
              ! not really the same thing as rc in PP generation
              upf%rcut  (iv) = upf%r(upf%kbeta(iv))
              upf%rcutus(iv) = 0.0
           ENDDO
        ENDDO
        upf%kkbeta = maxval(upf%kbeta(:))
     ELSE
        ALLOCATE(aux(upf%kkbeta))
        iv=0  ! counter on beta functions
        DO i=1,upf%nwfc
           l=upf%lchi(i)
           IF (l/=lloc) THEN
              iv=iv+1
              upf%lll(iv)=l
              upf%els_beta(iv)=upf%els(i)
              DO ir=1,upf%kbeta(iv)
                 ! the factor e2 converts from Hartree to Rydberg
                 upf%beta(ir,iv) = e2 * chi(ir,l+1) * &
                      ( vnl(ir,l) - vnl(ir,lloc) )
                 aux(ir) = chi(ir,l+1) * upf%beta(ir,iv)
              ENDDO
              upf%rcut  (iv) = upf%r(upf%kbeta(iv))
              upf%rcutus(iv) = 0.0
              CALL simpson(upf%kbeta(iv),aux,upf%rab,vll)
              upf%dion(iv,iv) = 1.0d0/vll
           ENDIF
        ENDDO
        DEALLOCATE(aux)
     ENDIF
  ELSE
     ! prevents funny errors when writing file
     ALLOCATE(upf%dion(upf%nbeta,upf%nbeta))
  ENDIF

  ALLOCATE (upf%chi(upf%mesh,upf%nwfc))
  upf%chi(:,:) = chi(1:upf%mesh,1:upf%nwfc)

  ALLOCATE (upf%rho_at(upf%mesh))
  upf%rho_at(:) = 0.d0
  DO i=1,upf%nwfc
     upf%rho_at(:) = upf%rho_at(:) + upf%oc(i) * upf%chi(:,i) ** 2
  ENDDO

  !     ----------------------------------------------------------
  WRITE (6,'(a)') 'Pseudopotential successfully converted'
  !     ----------------------------------------------------------

  RETURN
END SUBROUTINE convert_cpmd
!
! ------------------------------------------------------------------
SUBROUTINE check_radial_grid ( amesh_, mesh, r, amesh )
! ------------------------------------------------------------------
!
   IMPLICIT NONE
   INTEGER, INTENT (in) :: mesh
   REAL(8), INTENT (in) :: amesh_, r(mesh)
   REAL(8), INTENT (out) :: amesh
   INTEGER :: i
   !
   ! get amesh if not available directly, check its value otherwise
   PRINT  "('Radial grid r(i) has ',i4,' points')", mesh
   PRINT  "('Assuming log radial grid: r(i)=exp[(i-1)*amesh]*r(1), with:')"
   IF (amesh_ < 0.0d0) THEN
      amesh = log (r(mesh)/r(1))/(mesh-1)
      PRINT  "('amesh = log (r(mesh)/r(1))/(mesh-1) = ',f10.6)",amesh
   ELSE
   ! not clear whether the value of amesh read from file
   ! matches the above definition, or if it is exp(amesh) ...
      amesh = log (r(mesh)/r(1))/(mesh-1)
      IF ( abs ( amesh - amesh_ ) > 1.0d-5 ) THEN
         IF ( abs ( amesh - exp(amesh_) ) < 1.0d-5 ) THEN
            amesh = log(amesh_)
            PRINT  "('amesh = log (value read from file) = ',f10.6)",amesh
         ELSE
            CALL errore ('cpmd2upf', 'unknown real-space grid',2)
         ENDIF
      ELSE
         amesh = amesh_
         PRINT  "('amesh = value read from file = ',f10.6)",amesh
      ENDIF
   ENDIF
   ! check if the grid is what we expect
   DO i=2,mesh
      IF ( abs(r(i) - exp((i-1)*amesh)*r(1)) > 1.0d-5) THEN
         PRINT  "('grid point ',i4,': found ',f10.6,', expected ',f10.6)",&
              i, r(i),  exp((i-1)*amesh)*r(1)
         CALL errore ('cpmd2upf', 'unknown real-space grid',1)
      ENDIF
   ENDDO
   RETURN
END
! ------------------------------------------------------------------
REAL(8) FUNCTION mygamma ( n )
  !------------------------------------------------------------------
  !
  ! mygamma(n)  = \Gamma(n-1/2) = sqrt(pi)*(2n-3)!!/2**(n-1)
  !
  USE constants, ONLY : pi
  IMPLICIT NONE
  INTEGER, INTENT(in) :: n
  !
  REAL(8) :: x
  INTEGER, EXTERNAL :: semifact
  !
  IF ( n < 2 ) CALL errore('mygamma','unexpected input argument',1)
  mygamma = sqrt(pi) * semifact(2*n-3) / 2.d0**(n-1)
  !
  RETURN
END FUNCTION mygamma
! ------------------------------------------------------------------
INTEGER FUNCTION locate(onechar,string)
  ! ------------------------------------------------------------------
  !
  CHARACTER(len=1) :: onechar
  CHARACTER(len=*) :: string
  !
  INTEGER:: i
  !
  DO i=1,len_trim(string)
     IF (string(i:i) == "=") THEN
        locate = i
        RETURN
     ENDIF
  ENDDO
  locate = 0
  RETURN
END FUNCTION locate
