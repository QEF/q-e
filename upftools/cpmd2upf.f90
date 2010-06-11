!
! Copyright (C) 2001 PWSCF group
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
  !     Convert a pseudopotential written in the CPMD format
  !     (TYPE=NORMCONSERVING NUMERIC only, single radial grid)
  !     to unified pseudopotential format
  !
  IMPLICIT NONE
  CHARACTER(len=256) filein, fileout
  !
  !
  CALL get_file ( filein )
  OPEN (unit = 1, file = filein, status = 'old', form = 'formatted')
  CALL read_cpmd(1)
  CLOSE (1)

  ! convert variables read from CPMD format into those needed
  ! by the upf format - add missing quantities

  CALL convert_cpmd

  fileout=trim(filein)//'.UPF'
  PRINT '(''Output PP file in UPF format :  '',a)', fileout

  OPEN(unit=2,file=fileout,status='unknown',form='formatted')
  CALL write_upf(2)
  CLOSE (unit=2)

STOP
20 CALL errore ('cpmd2upf', 'Reading pseudo file name ', 1)

END PROGRAM cpmd2upf

MODULE cpmd
  !
  ! All variables read from CPMD file format
  !
  CHARACTER (len=80) title
  !
  INTEGER :: ixc
  real(8) :: alphaxc
  INTEGER :: z, zv
  !
  INTEGER :: mesh_
  real(8) :: amesh, amesh_
  real(8), ALLOCATABLE :: r_(:)
  !
  INTEGER ::lmax_
  real(8), ALLOCATABLE :: vnl(:,:)
  real(8), ALLOCATABLE :: chi_(:,:)
  !
  LOGICAL :: nlcc_
  real(8), ALLOCATABLE :: rho_atc_(:)
  !
  INTEGER :: maxinfo_, info_lines_
  PARAMETER (maxinfo_ = 100)
  CHARACTER (len=80), ALLOCATABLE :: info_sect_(:)
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
  real (8) :: vnl0(0:3)
  LOGICAL, EXTERNAL :: matches
  INTEGER, EXTERNAL :: locate
  !
  nlcc_ = .false.
  info_lines_ = 0
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
     IF (.not. matches("NORMCONSERVING",line) .or. &
         .not. matches("NUMERIC",line) ) &
             CALL errore('read_cpmd','unknown type: '//line,1)
  ELSEIF (matches ("&INFO", trim(line)) ) THEN
     found = found + 1
     ! read (iunps,'(a)') title
     ! store info section for later perusal (FIXME: not yet implemented. 2004/10/12, AK)
     ALLOCATE (info_sect_(maxinfo_))
     DO i=1,maxinfo_
        READ (iunps,'(a)',end=20,err=20) title
        IF (matches ("&END", trim(title)) ) THEN
           closed = closed + 1
           GOTO 10
        ELSE
           info_sect_(i) = trim(title)
           info_lines_ = i
        ENDIF
     ENDDO
  ELSEIF (matches ("&POTENTIAL", trim(line)) ) THEN
     found = found + 1
     READ (iunps,'(a)') line
     READ (line,*,iostat=ios) mesh_, amesh
     IF ( ios /= 0) THEN
        READ (line,*,iostat=ios) mesh_
        amesh = -1.0d0
     ENDIF
     ALLOCATE (r_(mesh_))
     !
     ! determine the number of angular momenta
     !
     READ (iunps, '(a)') line
     ios = 1
     lmax_=4
     DO WHILE (ios /= 0)
        lmax_ = lmax_ - 1
        READ(line,*,iostat=ios) r_(1),(vnl0(l),l=0,lmax_)
     ENDDO
     ALLOCATE (vnl(mesh_,0:lmax_))
     vnl(1,0:lmax_) = vnl0(0:lmax_)
     DO i=2,mesh_
        READ(iunps, *) r_(i),(vnl(i,l),l=0,lmax_)
     ENDDO
     ! get amesh if not available directly, check its value otherwise
     PRINT  "('Radial grid r(i) has ',i4,' points')", mesh_
     PRINT  "('Assuming log radial grid: r(i)=exp[(i-1)*amesh]*r(1), with:')"
     IF (amesh < 0.0d0) THEN
        amesh = log (r_(mesh_)/r_(1))/(mesh_-1)
        PRINT  "('amesh = log (r(mesh)/r(1))/(mesh-1) = ',f10.6)",amesh
     ELSE
        ! not clear whether the value of amesh read from file
        ! matches the above definition, or if it is exp(amesh) ...
        amesh_ = log (r_(mesh_)/r_(1))/(mesh_-1)
        IF ( abs ( amesh - amesh_ ) > 1.0d-5 ) THEN
           IF ( abs ( amesh - exp(amesh_) ) < 1.0d-5 ) THEN
               amesh = log(amesh)
               PRINT  "('amesh = log (value read from file) = ',f10.6)",amesh
           ELSE
               CALL errore ('cpmd2upf', 'unknown real-space grid',2)
           ENDIF
        ELSE
           PRINT  "('amesh = value read from file = ',f10.6)",amesh
        ENDIF
     ENDIF
     ! check if the grid is what we expect
     DO i=2,mesh_
        IF ( abs(r_(i) - exp((i-1)*amesh)*r_(1)) > 1.0d-5) THEN
            PRINT  "('grid point ',i4,': found ',f10.6,', expected ',f10.6)",&
                     i, r_(i),  exp((i-1)*amesh)*r_(1)
            CALL errore ('cpmd2upf', 'unknown real-space grid',1)
        ENDIF
     ENDDO
  ELSEIF (matches ("&WAVEFUNCTION", trim(line)) ) THEN
     found = found + 1
     ! read (iunps,*) mesh_, amesh
     READ (iunps,'(a)') line
     READ (line,*,iostat=ios) mesh_
     ALLOCATE(chi_(mesh_,lmax_+1))
     DO i=1,mesh_
        READ(iunps, *) r_(i),(chi_(i,l+1),l=0,lmax_)
     ENDDO
  ELSEIF (matches ("&NLCC", trim(line)) ) THEN
     found = found + 1
     nlcc_ = .true.
     READ (iunps, '(a)') line
     IF (.not. matches ("NUMERIC", trim(line)) ) &
          CALL errore('read_cpmd',' only NUMERIC core-correction supported',1)
     READ(iunps, *) mesh_
     ALLOCATE (rho_atc_(mesh_))
     READ(iunps, * ) (r_(i), rho_atc_(i), i=1,mesh_)
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
  IF (nlcc_ .and. found /= 5 .or. .not.nlcc_ .and. found /= 4) &
       CALL errore('read_cpmd','some &FIELD card missing',found)
  IF (closed /= found) &
       CALL errore('read_cpmd','some &END card missing',closed)
  IF (unknown /= 0 ) PRINT '("WARNING: ",i3," cards not read")', unknown

  RETURN
200 CALL errore('read_cpmd','error in reading file',1)

END SUBROUTINE read_cpmd

!     ----------------------------------------------------------
SUBROUTINE convert_cpmd
  !     ----------------------------------------------------------
  !
  USE cpmd
  USE upf
  IMPLICIT NONE
  real(8), PARAMETER :: rmax = 10.0d0
  real(8), ALLOCATABLE :: aux(:)
  real(8) :: vll
  CHARACTER (len=20):: dft
  CHARACTER (len=2), EXTERNAL :: atom_name
  INTEGER :: lloc, kkbeta, my_lmax
  INTEGER :: l, i, ir, iv
  !
  WRITE(generated, '("Generated using unknown code")')
  WRITE(date_author,'("Author: unknown    Generation date: as well")')
  comment = 'Info: automatically converted from CPMD format'

  ! NOTE: many CPMD pseudopotentials created with the 'Hamann' code
  ! from Juerg Hutter's homepage have additional (bogus) entries for
  ! pseudo-potential and wavefunction. In the 'report' they have
  ! the same rc and energy eigenvalue than the previous angular momentum.
  ! we need to be able to ignore that part or the resulting UPF file
  ! will be useless. so we first print the info section and ask
  ! for the LMAX to really use. AK 2005/03/30.
  DO i=1,info_lines_
        PRINT '(A)', info_sect_(i)
  ENDDO
  PRINT '("lmax to use. (max.",I2,") > ",$)', lmax_
  READ (5,*) my_lmax
  IF ((my_lmax <= lmax_) .and. (my_lmax >= 0)) lmax_ = my_lmax
  PRINT '("l local (max.",I2,") > ",$)', lmax_
  READ (5,*) lloc
  ! reasonable assumption
  IF (z > 18) THEN
     rel = 1
  ELSE
     rel = 0
  ENDIF
  rcloc = 0.0d0
  nwfs  = lmax_+1
  ALLOCATE( els(nwfs), oc(nwfs), epseu(nwfs))
  ALLOCATE(lchi(nwfs), nns(nwfs) )
  ALLOCATE(rcut (nwfs), rcutus (nwfs))
  DO i=1, nwfs
     PRINT '("Wavefunction # ",i1,": label, occupancy > ",$)', i
     READ (5,*) els(i), oc(i)
     nns (i)  = 0
     lchi(i)  = i-1
     rcut(i)  = 0.0d0
     rcutus(i)= 0.0d0
     epseu(i) = 0.0d0
  ENDDO
  psd   = atom_name (z)
  pseudotype = 'NC'
  nlcc = nlcc_
  zp = zv
  etotps =0.0d0
  ecutrho=0.0d0
  ecutwfc=0.0d0
  IF ( lmax_ == lloc) THEN
     lmax = lmax_-1
  ELSE
     lmax = lmax_
  ENDIF
  nbeta= lmax_
  mesh = mesh_
  ntwfc= nwfs
  ALLOCATE( elsw(ntwfc), ocw(ntwfc), lchiw(ntwfc) )
  DO i=1, nwfs
     lchiw(i) = lchi(i)
     ocw(i)   = oc(i)
     elsw(i)  = els(i)
  ENDDO
  iexch = ixc/1000
  icorr = (ixc-1000*iexch)/100
  igcx  = (ixc-1000*iexch-100*icorr)/10
  igcc  = (ixc-1000*iexch-100*icorr-10*igcx)
  !
  ! We have igcc=2 (PW91) and 3 (LYP) exchanged wrt CPMD conventions
  !
  IF (igcc==3) THEN
     igcc=2
  ELSEIF (igcc==2) THEN
     igcc=3
  ENDIF

  ALLOCATE(rab(mesh))
  ALLOCATE(  r(mesh))
  r = r_
  rab = r * amesh

  ALLOCATE (rho_atc(mesh))
  IF (nlcc) rho_atc = rho_atc_

  ALLOCATE (vloc0(mesh))
  ! the factor 2 converts from Hartree to Rydberg
  vloc0(:) = vnl(:,lloc)*2.d0

  IF (nbeta > 0) THEN

     ALLOCATE(ikk2(nbeta), lll(nbeta))
     kkbeta=mesh
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
     iv=0
     DO i=1,nwfs
        l=lchi(i)
        IF (l/=lloc) THEN
           iv=iv+1
           lll(iv)=l
           DO ir=1,kkbeta
              ! the factor 2 converts from Hartree to Rydberg
              betar(ir,iv) = 2.d0 * chi_(ir,l+1) * &
                   ( vnl(ir,l) - vnl(ir,lloc) )
              aux(ir) = chi_(ir,l+1) * betar(ir,iv)
           ENDDO
           CALL simpson(kkbeta,aux,rab,vll)
           dion(iv,iv) = 1.0d0/vll
        ENDIF
     ENDDO

  ENDIF

  ALLOCATE (rho_at(mesh))
  rho_at = 0.d0
  DO i=1,nwfs
     rho_at(:) = rho_at(:) + ocw(i) * chi_(:,i) ** 2
  ENDDO

  ALLOCATE (chi(mesh,ntwfc))
  chi = chi_
  !     ----------------------------------------------------------
  WRITE (6,'(a)') 'Pseudopotential successfully converted'
  !     ----------------------------------------------------------
  RETURN
END SUBROUTINE convert_cpmd
!
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
