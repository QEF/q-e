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
  CALL get_file ( filein )
  IF ( TRIM(filein) == ' ') &
       CALL errore ('cpmd2upf', 'usage: cpmd2upf "file-to-be-converted"', 1)
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
  INTEGER :: z, zv
  !
  real(8) :: alphaloc, alpha(0:3), a(0:3), b(0:3)
  INTEGER :: mesh
  real(8) :: amesh, rmax, xmin
  real(8), ALLOCATABLE :: r(:)
  !
  INTEGER ::lmax
  real(8), ALLOCATABLE :: vnl(:,:)
  real(8), ALLOCATABLE :: chi(:,:)
  !
  LOGICAL :: nlcc
  real(8), ALLOCATABLE :: rho_atc(:)
  !
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
  INTEGER :: i, l, ios
  CHARACTER (len=80) line
  CHARACTER (len=4) token
  real (8) :: amesh_, vnl0(0:3)
  LOGICAL :: grid_read = .FALSE.
  LOGICAL, EXTERNAL :: matches
  INTEGER, EXTERNAL :: locate
  REAL(8), EXTERNAL :: qe_erf
  !
  nlcc = .false.
  info_lines = 0
10 READ (iunps,'(A)',end=20,err=20) line
  IF (matches ("&ATOM", trim(line)) ) THEN
     found = found + 1
     ! Z
     READ (iunps,'(a)',end=200,err=200) line
     l = len_trim(line)
     i = locate('=',line)
     READ (line(i+1:l),*) z
     ! ZV
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
     IF ( matches("NORMCONSERVING",line) .AND. ( matches("NUMERIC",line) &
                                     .OR. matches("LOGARITHMIC",line)  ) )THEN
        pstype = 1
     ELSE IF ( matches("NORMCONSERVING",line) .AND. matches("CAR",line) ) THEN
        pstype = 2
     END IF
     IF (pstype == 0 ) CALL errore('read_cpmd','unknown type: '//line,1)
  ELSEIF (matches ("&INFO", trim(line)) ) THEN
     found = found + 1
     ! read (iunps,'(a)') title
     ! store info section for later perusal (FIXME: not yet implemented. 2004/10/12, AK)
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
        IF ( .NOT. grid_read ) ALLOCATE (r(mesh))
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
        IF ( .NOT.grid_read ) THEN
           CALL check_radial_grid ( amesh_, mesh, r, amesh )
           grid_read = .TRUE.
        END IF
     ELSE
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
              EXIT
           END IF
           READ(line, *) alpha(lmax+1), a(lmax+1), b(lmax+1)
           alpha(lmax+1) = 1.d0/alpha(lmax+1)**2
        END DO
     END IF
  ELSEIF (matches ("&WAVEFUNCTION", trim(line)) ) THEN
     found = found + 1
     READ (iunps,'(A)') line
     READ (line,*,iostat=ios) mesh, amesh_
     IF ( ios /= 0) THEN
        READ (line,*,iostat=ios) mesh
        amesh_ = -1.0d0
     ENDIF
     IF ( .NOT. grid_read )  ALLOCATE(r(mesh))
     ALLOCATE(chi(mesh,lmax+1))
     DO i=1,mesh
        READ(iunps, *) r(i),(chi(i,l+1),l=0,lmax)
     ENDDO
     IF ( .NOT.grid_read ) THEN
        CALL check_radial_grid ( amesh_, mesh, r, amesh )
        grid_read = .TRUE.
     END IF
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
     END IF
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
   rmax = r(mesh)
   xmin = log(z*r(1))
  IF (nlcc .and. found /= 5 .or. .not.nlcc .and. found /= 4) &
       CALL errore('read_cpmd','some &FIELD card missing',found)
  IF (closed /= found) &
       CALL errore('read_cpmd','some &END card missing',closed)
  IF (unknown /= 0 ) PRINT '("WARNING: ",i3," cards not read")', unknown
  !
  IF ( pstype == 2 ) THEN
     ALLOCATE (vnl(mesh,0:lmax))
     DO l=0, lmax
        DO i=1, mesh
           vnl(i,l) = ( a(l) + b(l)*r(i)**2 ) * exp (-alpha(l)*r(i)**2) - &
                       zv * qe_erf (sqrt(alphaloc)*r(i))/r(i)
        END DO
     END DO
  END IF

  RETURN
200 CALL errore('read_cpmd','error in reading file',1)

END SUBROUTINE read_cpmd

!     ----------------------------------------------------------
SUBROUTINE convert_cpmd(upf)
  !     ----------------------------------------------------------
  !
  USE cpmd
  USE pseudo_types, ONLY : pseudo_upf
  USE funct, ONLY :  dft_name
  !
  IMPLICIT NONE
  !
  TYPE(pseudo_upf) :: upf 
  !
  real(8), ALLOCATABLE :: aux(:)
  real(8) :: vll, rcloc
  CHARACTER (len=20):: dft
  CHARACTER (len=2):: label
  CHARACTER (len=2), EXTERNAL :: atom_name
  INTEGER :: lloc, my_lmax
  INTEGER :: iexch, icorr, igcx, igcc, inlc, l, i, ir, iv
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
  PRINT '("max L to use ( <= ",I1," ) > ",$)', lmax
  READ (5,*) my_lmax
  IF ((my_lmax <= lmax) .and. (my_lmax >= 0)) lmax = my_lmax
  PRINT '("local L ( <= ",I1," ), Rc for local pot (au) > ",$)', lmax
  READ (5,*) lloc, rcloc
  !
  upf%nv       = "2.0.1"
  upf%generated= "Generated using unknown code"
  upf%author   = "Author: unknown"
  upf%date     = "Generation date: as well"
  upf%comment  = "Info: automatically converted from CPMD format"
  upf%psd      = atom_name (z)
  ! reasonable assumption
  IF (z > 18) THEN
     upf%rel = 'no'
  ELSE
     upf%rel = 'scalar'
  ENDIF
  upf%typ = 'SL'
  upf%tvanp = .FALSE.
  upf%tpawp = .FALSE.
  upf%tcoulombp=.FALSE.
  upf%nlcc = nlcc
  !
  iexch = ixc/1000
  icorr = (ixc-1000*iexch)/100
  igcx  = (ixc-1000*iexch-100*icorr)/10
  igcc  = (ixc-1000*iexch-100*icorr-10*igcx)
  ! We have igcc=2 (PW91) and 3 (LYP) exchanged wrt CPMD conventions
  IF (igcc==3) THEN
     igcc=2
  ELSEIF (igcc==2) THEN
     igcc=3
  ENDIF
  inlc = 0
  CALL dft_name (iexch, icorr, igcx, igcc, inlc, upf%dft, dft)
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
  upf%nwfc = lmax+1
  !
  ALLOCATE( upf%els(upf%nwfc) )
  ALLOCATE( upf%oc(upf%nwfc) )
  ALLOCATE( upf%epseu(upf%nwfc) )
  ALLOCATE( upf%lchi(upf%nwfc) )
  ALLOCATE( upf%nchi(upf%nwfc) )
  ALLOCATE( upf%rcut_chi (upf%nwfc) )
  ALLOCATE( upf%rcutus_chi(upf%nwfc) )

  DO i=1, upf%nwfc
10   PRINT '("Wavefunction # ",i1,": label (e.g. 4s), occupancy > ",$)', i
     READ (5,*) label, upf%oc(i)
     READ (label(1:1),*, err=10) l
     upf%els(i)  = label
     upf%nchi(i)  = l
     IF ( label(2:2) == 's' .OR. label(2:2) == 'S') then
        l=0
     ELSE IF ( label(2:2) == 'p' .OR. label(2:2) == 'P') then
        l=1
     ELSE IF ( label(2:2) == 'd' .OR. label(2:2) == 'D') then
        l=2
     ELSE IF ( label(2:2) == 'f' .OR. label(2:2) == 'F') then
        l=3
     ELSE
        l=i-1
     END IF
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

  ALLOCATE (upf%vloc(upf%mesh))
  ! the factor 2 converts from Hartree to Rydberg
  upf%vloc(:) = vnl(1:upf%mesh,upf%lloc)*2.d0
  upf%rcloc = rcloc

  ALLOCATE(upf%vnl(upf%mesh,0:upf%lmax,1))
  upf%vnl(:,:,1) = vnl(1:upf%mesh,0:upf%lmax)

  upf%nbeta= lmax
  IF (upf%nbeta > 0) THEN

     ALLOCATE(upf%els_beta(upf%nbeta) )
     ALLOCATE(upf%lll(upf%nbeta))
     ALLOCATE(upf%kbeta(upf%nbeta))
     iv=0  ! counter on beta functions
     DO i=1,upf%nwfc
        l=upf%lchi(i)
        IF (l/=lloc) THEN
           iv=iv+1
           DO ir = upf%mesh,1,-1
              IF ( ABS ( vnl(ir,l) - vnl(ir,lloc) ) > 1.0E-6 ) THEN
                 upf%kbeta(iv)=ir
                 exit
              ENDIF
           ENDDO
        ENDIF
     ENDDO
     ! the number of points used in the evaluation of integrals
     ! should be even (for simpson integration)
     DO i=1,upf%nbeta
        IF ( MOD (upf%kbeta(i),2) == 0 .AND. upf%kbeta(i) < upf%mesh) &
           upf%kbeta(i)=upf%kbeta(i)+1
     END DO
     upf%kkbeta = MAXVAL(upf%kbeta(:))
     ALLOCATE(upf%beta(upf%mesh,upf%nbeta))
     ALLOCATE(upf%dion(upf%nbeta,upf%nbeta))
     upf%beta(:,:) =0.d0
     upf%dion(:,:) =0.d0
     ALLOCATE(upf%rcut  (upf%nbeta))
     ALLOCATE(upf%rcutus(upf%nbeta))
     ALLOCATE(aux(upf%kkbeta))
     iv=0  ! counter on beta functions
     DO i=1,upf%nwfc
        l=upf%lchi(i)
        IF (l/=lloc) THEN
           iv=iv+1
           upf%lll(iv)=l
           upf%els_beta(iv)=upf%els(i)
           DO ir=1,upf%kbeta(iv)
              ! the factor 2 converts from Hartree to Rydberg
              upf%beta(ir,iv) = 2.d0 * chi(ir,l+1) * &
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
   INTEGER, INTENT (IN) :: mesh   
   REAL(8), INTENT (IN) :: amesh_, r(mesh)
   REAL(8), INTENT (OUT) :: amesh
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
