!
! Copyright (C) 2001-2011 Quantum ESPRESSO group
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
  !     to unified pseudopotential format (v.2)
  !     Adapted from the converter written by Andrea Ferretti
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
  IF ( trim(filein) == ' ') &
       CALL errore ('fhi2upf', 'usage: fhi2upf "file-to-be-converted"', 1)
  OPEN ( unit=1, file=filein, status = 'old', form='formatted', iostat=ios )
  IF ( ios /= 0) CALL errore ('fhi2upf', 'file: '//trim(filein)//' not found', 2)
  !
  CALL read_fhi(1)
  !
  CLOSE (1)

  ! convert variables read from FHI format into those needed
  ! by the upf format - add missing quantities
  !
  CALL nullify_pseudo_upf ( upf )
  !
  CALL convert_fhi (upf)
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
  INTEGER      :: lmax          ! max l-component used

  LOGICAL      :: nlcc_
  real(8), ALLOCATABLE :: rho_atc(:) ! core  charge

  TYPE (angular_comp), POINTER :: comp(:)  ! PP numerical info
                                           ! (wfc, grid, potentials...)
  !------------------------------

  ! variables for the abinit header

  real(8) :: Zatom, Zion, r2well, rchrg, fchrg, qchrg
  INTEGER :: pspdat = 0, pspcod = 0 , pspxc = 0, lloc = -1, mmax = 0
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
  real(8) :: r, drhoc, d2rhoc
  !

  INTEGER               :: l, i, idum, mesh

  ! Start reading file

  READ(iunps,'(a)') info
  READ(info,*,iostat=i) Zval, l
  IF ( i /= 0 .or. zval <= 0.0 .or. zval > 100.0 ) THEN
     WRITE (6,'("Assuming abinit format. First line:",/,A)') trim(info)
     READ(iunps,*) Zatom, Zion, pspdat
     READ(iunps,*) pspcod, pspxc, lmax,lloc, mmax, r2well
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
     IF (l-1 /= lmax) THEN
        WRITE (6,'("read_fhi: lmax mismatch...stopping")')
        STOP
     ENDIF
  ELSE
     info = ' '
  ENDIF
  lmax = l - 1

  IF (lmax+1 > Nl) THEN
     WRITE (6,'("read_fhi: too many l-components...stopping")')
     STOP
  ENDIF

  DO i=1,10
     READ(iunps,*)     ! skipping 11 lines
  ENDDO

  ALLOCATE( comp(0:lmax) )

  DO l=0,lmax
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
  ALLOCATE(rho_atc(comp(0)%nmesh))
  mesh = comp(0)%nmesh
  DO i=1,mesh
     READ(iunps,*,end=10, err=20) r, rho_atc(i), drhoc, d2rhoc
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
20 WRITE(6,'("read_fhi: error reading core charge, assuming no core charge")')
   WRITE(6,'("this error may be due to the presence of additional", &
 &           " lines at the end of file")')
10 CONTINUE
  !     ----------------------------------------------------------
  WRITE (6,'(a)') 'Pseudopotential without NLCC successfully read'
  !     ----------------------------------------------------------
  RETURN
  !
  STOP

END SUBROUTINE read_fhi

!     ----------------------------------------------------------
SUBROUTINE convert_fhi (upf)
  !     ----------------------------------------------------------
  !
  USE fhi
  USE pseudo_types, ONLY : pseudo_upf
  USE funct, ONLY : set_dft_from_name, get_iexch, get_icorr, get_igcx, get_igcc
  USE constants, ONLY : fpi
  !
  IMPLICIT NONE
  !
  TYPE(pseudo_upf) :: upf
  !
  real(8), ALLOCATABLE :: aux(:)
  real(8) :: vll
  CHARACTER (len=2):: label
  CHARACTER (len=2), EXTERNAL:: atom_name
  INTEGER :: l, i, ir, iv
  !
  upf%nv       = "2.0.1"
  upf%generated= "Generated using FHI98PP, converted with fhi2upf.x v.5.0.2"
  upf%author   = "unknown"
  upf%date     = "unknown"
  IF (trim(info) /= ' ') THEN
     upf%comment = trim(info)
  ELSE
     upf%comment = 'Info: automatically converted from FHI format'
  ENDIF
  upf%rel = 'scalar'  ! just guessing
  IF (nint(Zatom) > 0) THEN
     upf%psd = atom_name(nint(Zatom))
     IF (nint(Zatom) < 18) upf%rel = 'no' ! just guessing
  ELSE
     WRITE(*,'("Atom name > ")', advance="NO")
     READ (5,'(a)') upf%psd
  ENDIF
  upf%typ = 'SL'
  upf%tvanp = .false.
  upf%tpawp = .false.
  upf%tcoulombp=.false.
  upf%nlcc = nlcc_
  !
  IF (pspxc == 7) THEN
     upf%dft = 'SLA-PW'
  ELSEIF (pspxc == 11) THEN
     upf%dft = 'PBE'
  ELSE
     IF (pspxc > 0) THEN
        PRINT '("DFT read from abinit file: ",i1)', pspxc
     ENDIF
     WRITE(*,'("DFT > ")', advance="NO")
     READ (5,'(a)') upf%dft
  ENDIF
  !
  upf%zp   = Zval
  upf%etotps =0.0d0
  upf%ecutrho=0.0d0
  upf%ecutwfc=0.0d0
  !
! 2014/11/11 JM
! Use lloc (from fhi module) here, otherwise the user input has not effect on the vloc assignment below.
! In case of a .cpi file (i.e. direct output from Martin Fuchs's FHI98PP code), which does not include information about lloc,
! a proper conversion could actually never have been achieved in the past.
  WRITE(*,'("Confirm or modify l max, l loc (read:",2i3,") > ")', advance="NO") lmax, lloc
  READ (5,*) lmax, lloc
  !
  IF ( lmax == lloc) THEN
     upf%lmax = lmax-1
  ELSE
     upf%lmax = lmax
  ENDIF
  upf%lloc = lloc
  upf%lmax_rho = 0
  upf%nwfc  = lmax+1
  !
  ALLOCATE( upf%els(upf%nwfc) )
  ALLOCATE( upf%oc(upf%nwfc) )
  ALLOCATE( upf%epseu(upf%nwfc) )
  ALLOCATE( upf%lchi(upf%nwfc) )
  ALLOCATE( upf%nchi(upf%nwfc) )
  ALLOCATE( upf%rcut_chi (upf%nwfc) )
  ALLOCATE( upf%rcutus_chi(upf%nwfc) )

  PRINT '("PPs in FHI format do not contain information on atomic valence (pseudo-)wavefunctions")'
  PRINT '("Provide the label and the occupancy for each atomic wavefunction used in the PP generation")'
  PRINT '("If unknown: list valence wfcts and occupancies for the atomic ground state ", &
         &"in increasing l order: s,p,d,f")'
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

  upf%mesh = comp(0)%nmesh
  upf%dx   = log( comp(0)%amesh )
  upf%rmax = comp(0)%grid(upf%mesh)
  upf%xmin = log( comp(0)%grid(1)*Zatom )
  upf%zmesh= Zatom
  ALLOCATE(upf%rab(upf%mesh))
  ALLOCATE(upf%r(upf%mesh))
  upf%r(:) = comp(0)%grid
  upf%rab(:)=upf%r(:)*upf%dx

  ALLOCATE (upf%rho_atc(upf%mesh))
  IF (upf%nlcc) upf%rho_atc(:) = rho_atc(1:upf%mesh) / fpi

  ALLOCATE (upf%vloc(upf%mesh))
  ! the factor 2 converts from Hartree to Rydberg
  upf%vloc(:) = 2.d0*comp(lloc)%pot
  upf%rcloc = 0.0d0

  ALLOCATE(upf%vnl(upf%mesh,0:upf%lmax,1))
  DO l=0, upf%lmax
     upf%vnl(:,l,1) = 2.d0*comp(l)%pot(:)
  ENDDO

  ! calculate number of nonlocal projectors
  IF ( upf%lloc >= 0 .and. upf%lloc <= upf%lmax ) THEN
     upf%nbeta= upf%lmax
  ELSE
     upf%nbeta= upf%lmax+1
  ENDIF

  IF (upf%nbeta > 0) THEN

     ALLOCATE(upf%els_beta(upf%nbeta) )
     ALLOCATE(upf%lll(upf%nbeta))
     ALLOCATE(upf%kbeta(upf%nbeta))
     iv=0  ! counter on beta functions
     DO i=1,upf%nwfc
        l=upf%lchi(i)
        IF (l/=upf%lloc) THEN
           iv=iv+1
           upf%kbeta(iv)=upf%mesh
           DO ir = upf%mesh,1,-1
! 2014/11/11 JM
! Explicitly use local potential here:
! Otherwise the IF (lmax == upf%lloc) clause above leads to unnecessary ("maximum") large radial grids for the beta projectors.
!              IF ( abs ( upf%vnl(ir,l,1) - upf%vnl(ir,upf%lloc,1) ) > 1.0E-6 ) THEN
              IF ( abs ( upf%vnl(ir,l,1) - upf%vloc(ir) ) > 1.0E-6 ) THEN
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
        IF (l/=upf%lloc) THEN
           iv=iv+1
           upf%lll(iv)=l
           upf%els_beta(iv)=upf%els(i)
           DO ir=1,upf%kbeta(iv)
              ! the factor 2 converts from Hartree to Rydberg
              upf%beta(ir,iv) = 2.d0 * comp(l)%wfc(ir) * &
                   ( comp(l)%pot(ir) - comp(upf%lloc)%pot(ir) )
              aux(ir) = comp(l)%wfc(ir) * upf%beta(ir,iv)
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
  DO i=1,upf%nwfc
     upf%chi(:,i) = comp(i-1)%wfc(:)
  ENDDO

  ALLOCATE (upf%rho_at(upf%mesh))
  upf%rho_at(:) = 0.d0
  DO i=1,upf%nwfc
     upf%rho_at(:) = upf%rho_at(:) + upf%oc(i) * upf%chi(:,i) ** 2
  ENDDO
  !     ----------------------------------------------------------
  WRITE (6,'(a)') 'Pseudopotential successfully converted'
  !     ----------------------------------------------------------
  RETURN
END SUBROUTINE convert_fhi
