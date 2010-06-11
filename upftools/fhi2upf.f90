!
! Copyright (C) 2001-2006 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file 'License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!---------------------------------------------------------------------
PROGRAM fhi2upf
  !---------------------------------------------------------------------
  !
  !     Convert a pseudopotential file in Fritz-Haber numerical format
  !     either ".cpi" (fhi88pp) or ".fhi" (abinit)
  !     to unified pseudopotential format
  !     May or may not work: carefully check what you get
  !     Adapted from the converter written by Andrea Ferretti
  !
  IMPLICIT NONE
  CHARACTER(len=256) filein, fileout
  !
  !
  CALL get_file ( filein )
  OPEN (unit = 1, file = filein, status = 'old', form = 'formatted')
  CALL read_fhi(1)
  CLOSE (1)

  ! convert variables read from FHI format into those needed
  ! by the upf format - add missing quantities
  CALL convert_fhi

  fileout=trim(filein)//'.UPF'
  PRINT '(''Output PP file in UPF format :  '',a)', fileout

  OPEN(unit=2,file=fileout,status='unknown',form='formatted')
  CALL write_upf(2)
  CLOSE (unit=2)

STOP
20 WRITE (6,'("fhi2upf: error reading pseudopotential file name")')
   STOP
END PROGRAM fhi2upf

MODULE fhi
  !
  ! All variables read from FHI file format
  !

  TYPE angular_comp
     real(8), POINTER     :: pot(:)
     real(8), POINTER     :: wfc(:)
     real(8), POINTER     :: grid(:)
     real(8)              :: amesh
     INTEGER             :: nmesh
     INTEGER             :: lcomp
  END TYPE angular_comp

  !------------------------------

  real(8) :: Zval           ! valence charge
  INTEGER      :: lmax_          ! max l-component used

  LOGICAL      :: nlcc_
  real(8), ALLOCATABLE :: rho_atc_(:) ! core  charge

  TYPE (angular_comp), POINTER :: comp(:)  ! PP numerical info
                                           ! (wfc, grid, potentials...)
  !------------------------------

  ! variables for the abinit header

  real(8) :: Zatom, Zion, r2well, rchrg, fchrg, qchrg
  INTEGER :: pspdat = 0, pspcod = 0 , pspxc = 0, lloc_ = -1, mmax = 0
  CHARACTER(len=256) :: info

END MODULE fhi
!
!     ----------------------------------------------------------
SUBROUTINE read_fhi(iunps)
  !     ----------------------------------------------------------
  !
  USE fhi
  IMPLICIT NONE
  INTEGER, PARAMETER    :: Nl=7  ! max number of l-components
  INTEGER :: iunps
  real(8) :: r, rhoc, drhoc, d2rhoc
  !

  INTEGER               :: l, i, idum, mesh

  ! Start reading file

  READ(iunps,'(a)') info
  READ(info,*,iostat=i) Zval, l
  IF ( i /= 0 .or. zval <= 0.0 .or. zval > 100.0 ) THEN
     WRITE (6,'("read_fhi: assuming abinit format")')
     READ(iunps,*) Zatom, Zion, pspdat
     READ(iunps,*) pspcod, pspxc, lmax_,lloc_, mmax, r2well
     IF (pspcod /= 6) THEN
        WRITE (6,'("read_fhi: unknown PP type ",i1,"...stopping")') pspcod
        STOP
     ENDIF
     READ(iunps,*) rchrg, fchrg, qchrg
     !
     READ(iunps,*)
     READ(iunps,*)
     READ(iunps,*)
     !
     READ(iunps,*) Zval, l
     IF (abs(Zion-Zval) > 1.0d-8) THEN
        WRITE (6,'("read_fhi: Zval/Zion mismatch...stopping")')
        STOP
     ENDIF
     IF (l-1 /= lmax_) THEN
        WRITE (6,'("read_fhi: lmax mismatch...stopping")')
        STOP
     ENDIF
  ELSE
     info = ' '
  ENDIF
  lmax_ = l - 1

  IF (lmax_+1 > Nl) THEN
     WRITE (6,'("read_fhi: too many l-components...stopping")')
     STOP
  ENDIF

  DO i=1,10
     READ(iunps,*)     ! skipping 11 lines
  ENDDO

  ALLOCATE( comp(0:lmax_) )

  DO l=0,lmax_
     comp(l)%lcomp = l
     READ(iunps,*) comp(l)%nmesh, comp(l)%amesh
     IF (mmax > 0 .and. mmax /= comp(l)%nmesh) THEN
        WRITE (6,'("read_fhi: mismatched number of grid points...stopping")')
        STOP
     ENDIF
     IF ( l > 0) THEN
        IF (comp(l)%nmesh /= comp(0)%nmesh .or.   &
            comp(l)%amesh /= comp(0)%amesh )      THEN
           WRITE(6,'("read_fhi: different radial grids not allowed...stopping")')
           STOP
        ENDIF
     ENDIF
     mesh = comp(l)%nmesh
     ALLOCATE( comp(l)%wfc(mesh),            &      ! wave-functions
               comp(l)%pot(mesh),            &      ! potentials
               comp(l)%grid(mesh)            )      ! real space radial grid
     ! read the above quantities
     DO i=1,mesh
        READ(iunps,*) idum, comp(l)%grid(i),   &
                            comp(l)%wfc(i),    &
                            comp(l)%pot(i)
     ENDDO
  ENDDO

  nlcc_ =.false.
  ALLOCATE(rho_atc_(comp(0)%nmesh))
  mesh = comp(0)%nmesh
  DO i=1,mesh
     READ(iunps,*,end=10, err=20) r, rho_atc_(i), drhoc, d2rhoc
     IF ( abs( r - comp(0)%grid(i) ) > 1.d-6 ) THEN
        WRITE(6,'("read_fhi: radial grid for core charge? stopping")')
        STOP
     ENDIF
  ENDDO
  nlcc_ = .true.
  !     ----------------------------------------------------------
  WRITE (6,'(a)') 'Pseudopotential with NLCC successfully read'
  !     ----------------------------------------------------------
  RETURN
10 CONTINUE
  !     ----------------------------------------------------------
  WRITE (6,'(a)') 'Pseudopotential without NLCC successfully read'
  !     ----------------------------------------------------------
  RETURN
  !
20 WRITE(6,'("read_fhi: error reading core charge")')
  STOP
  !
100  WRITE(6,'("read_fhi: error reading pseudopotential file")')
  STOP

END SUBROUTINE read_fhi

!     ----------------------------------------------------------
SUBROUTINE convert_fhi
  !     ----------------------------------------------------------
  !
  USE fhi
  USE upf
  USE funct, ONLY : set_dft_from_name, get_iexch, get_icorr, get_igcx, get_igcc
  USE constants, ONLY : fpi
  IMPLICIT NONE
  real(8), PARAMETER :: rmax = 10.0d0
  real(8), ALLOCATABLE :: aux(:)
  real(8) :: vll
  CHARACTER (len=20):: dft
  CHARACTER (len=2), EXTERNAL:: atom_name
  INTEGER :: lloc, kkbeta
  INTEGER :: l, i, ir, iv
  !
  IF (nint(Zatom) > 0) THEN
     psd = atom_name(nint(Zatom))
  ELSE
     PRINT '("Atom name > ",$)'
     READ (5,'(a)') psd
  ENDIF
  IF ( lloc_ < 0 ) THEN
     PRINT '("l local (max: ",i1,") > ",$)', lmax_
     READ (5,*) lloc
  ELSE
     lloc = lloc_
  ENDIF
  IF (pspxc == 7) THEN
     dft = 'PW'
  ELSE
     IF (pspxc > 0) THEN
        PRINT '("DFT read from abinit file: ",i1)', pspxc
     ENDIF
     PRINT '("DFT > ",$)'
     READ (5,'(a)') dft
  ENDIF
  WRITE(generated, '("Generated using Fritz-Haber code")')
  WRITE(date_author,'("Author: unknown    Generation date: as well")')
  IF (trim(info) /= ' ') THEN
     comment = trim(info)
  ELSE
     comment = 'Info: automatically converted from FHI format'
  ENDIF
  ! reasonable assumption
  rel = 1
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

  pseudotype = 'NC'
  nlcc = nlcc_
  zp   = Zval
  etotps = 0.0d0
  ecutrho=0.0d0
  ecutwfc=0.0d0
  IF ( lmax_ == lloc) THEN
     lmax = lmax_-1
  ELSE
     lmax = lmax_
  ENDIF
  nbeta= lmax_
  mesh = comp(0)%nmesh
  ntwfc= nwfs
  ALLOCATE( elsw(ntwfc), ocw(ntwfc), lchiw(ntwfc) )
  DO i=1, nwfs
     lchiw(i) = lchi(i)
     ocw(i)   = oc(i)
     elsw(i)  = els(i)
  ENDDO
  CALL set_dft_from_name(dft)
  iexch = get_iexch()
  icorr = get_icorr()
  igcx  = get_igcx()
  igcc  = get_igcc()

  ALLOCATE(rab(mesh))
  ALLOCATE(  r(mesh))
  r = comp(0)%grid
  rab = r * log( comp(0)%amesh )

  IF (nlcc) THEN
     ALLOCATE (rho_atc(mesh))
     rho_atc(:) = rho_atc_(:) / fpi
  ENDIF

  ALLOCATE (vloc0(mesh))
  ! the factor 2 converts from Hartree to Rydberg
  vloc0 = 2.d0*comp(lloc)%pot

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
              ! FHI potentials are in Hartree
              betar(ir,iv) = 2.d0 * comp(l)%wfc(ir) * &
                   ( comp(l)%pot(ir) - comp(lloc)%pot(ir) )
              aux(ir) = comp(l)%wfc(ir) * betar(ir,iv)
           ENDDO
           CALL simpson(kkbeta,aux,rab,vll)
           dion(iv,iv) = 1.0d0/vll
        ENDIF
     ENDDO

  ENDIF

  ALLOCATE (rho_at(mesh))
  rho_at = 0.d0
  DO i=1,nwfs
     l=lchi(i)
     rho_at = rho_at + ocw(i) * comp(l)%wfc ** 2
  ENDDO

  ALLOCATE (chi(mesh,ntwfc))
  DO i=1,ntwfc
     chi(:,i) = comp(i-1)%wfc(:)
  ENDDO
  !     ----------------------------------------------------------
  WRITE (6,'(a)') 'Pseudopotential successfully converted'
  !     ----------------------------------------------------------
  RETURN
END SUBROUTINE convert_fhi
