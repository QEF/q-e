!
! Copyright (C) 2023 Quantum ESPRESSO Foundation
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------------------------
subroutine read_psml ( filename, upf )
  !-----------------------------------------------------
  !! Read pseudopotential files in PSML format using "xmltools"
  !! stores data into the "upf" structure
  !! Main differences between the two formats:
  !! - in PSML arrays: local potential, projectors, charges, ... 
  !!   may be shorter than the full length of the radial grid,
  !!   as specified in tag argument "npts";
  !!   in UPF they all have the same full grid size "upf%mesh" so they
  !!   must be completed (local potential) or set to zero beyond npts
  !! - PSML pseudo-core charge must be divided by 4\pi
  !! - PSML projectors must be multiplied by r
  !! Works only for a subset of PSML files
  !! Written by P. Giannozzi, April 2023
  !
  USE xmltools
  USE upf_kinds, ONLY : dp
  USE upf_const, ONLY : e2, fpi
  USE pseudo_types, ONLY: pseudo_upf
  !
  IMPLICIT NONE
  !
  CHARACTER(len=*), INTENT(IN) :: filename
  !! input : name of file in psml format
  TYPE(pseudo_upf), INTENT(OUT) :: upf
  !! the derived type storing the pseudo data
  INTEGER  :: ierr
  !! error code (0 if correctly read)
  CHARACTER(len=30) :: tag
  !! tag where error (ierr != 0) was detected
  INTEGER :: iun
  !! unit for reading data
  INTEGER :: n
  !
  ierr = 0
  iun = xml_open_file ( filename )
  IF ( iun == -1 ) THEN
     ierr = 1
     tag = 'file'
     GO TO 10 
  END IF
  !
  tag = 'psml'
  call xmlr_opentag ( trim(tag), IERR = ierr )
  IF ( ierr /= 0 ) RETURN
  call get_attr ( 'version', upf%nv )
  ! print *, 'version=',upf%nv
  call get_attr ( 'uuid', upf%author )
  !
  tag = 'provenance'
  call xmlr_opentag ( trim(tag), IERR = ierr  )
  IF ( ierr /= 0 ) GO TO 10
  call get_attr ( 'creator', upf%generated )
  upf%author='UUID: '//trim(upf%author)
  call get_attr ( 'date', upf%date )
  upf%comment='PSML file v. '//trim(upf%nv)
  call xmlr_closetag ( ) ! provenance'
  !
  tag = 'pseudo-atom-spec'
  call read_psml_pseudo_atom_spec ( ierr )
  IF ( ierr /= 0 ) GO TO 10
  !
  tag = 'grid'
  call read_psml_grid ( ierr )
  IF ( ierr /= 0 ) GO TO 10
  !
  tag = 'valence-charge'
  call read_psml_radialfunc ( tag, upf%rho_at,ierr )
  IF ( ierr /= 0 ) GO TO 10
  !
  IF ( upf%nlcc ) THEN
     tag = 'pseudocore-charge'
     call read_psml_radialfunc ( tag, upf%rho_atc,ierr )
     IF ( ierr /= 0 ) GO TO 10
     upf%rho_atc(:) = upf%rho_atc(:) / fpi
  END IF
  !
  tag = 'local-potential'
  call read_psml_radialfunc ( tag, upf%vloc, ierr )
  IF ( ierr /= 0 ) GO TO 10
  ! size(vloc) < size(mesh) so vloc is filled with zeros
  ! fill vloc with asymptotic behavior until size(mesh) 
  do n=upf%mesh,1,-1
     if ( upf%vloc(n) == 0.0_dp ) then
        upf%vloc(n) = -upf%zp/upf%r(n)
     else
        exit
     end if
  end do
  !
  tag = 'nonlocal-projectors'
  call read_psml_nonlocal_projectors ( ierr )
  IF ( ierr /= 0 ) GO TO 10
  !
  tag = 'pseudo-wave-functions'
  call read_psml_pseudo_wave_functions ( ierr )
  IF ( ierr /= 0 ) THEN
     print *, 'read_psml: tag ',trim(tag),' not present'
     upf%nwfc = 0 
  END IF
  !
  call xmlr_closetag ( ) ! psml
  call xml_closefile( )
  !
  ! Convert from Hartree (PSML) to Ry
  !
  upf%vloc(:)   = e2*upf%vloc(:)
  upf%dion(:,:) = e2*upf%dion(:,:)
  !
  RETURN
10 print *, 'read_psml: error reading tag ',trim(tag)
   stop
  !
CONTAINS
  !
  SUBROUTINE read_psml_pseudo_atom_spec ( ierr )
    !
    INTEGER :: ierr
    INTEGER :: n, nxc, ndum
    INTEGER :: xc(6)
    CHARACTER(len=3) :: cc
    CHARACTER(len=25), external :: libxc_to_qe
    
    !
    upf%tvanp = .false.
    upf%tpawp = .false.
    upf%has_so = .false.
    upf%has_gipaw = .false.
    upf%paw_as_gipaw = .false.
    upf%tcoulombp = .false.
    upf%is_gth = .false.
    upf%is_multiproj = .false.
    upf%typ = 'NC'
    call xmlr_opentag ( 'pseudo-atom-spec', IERR = ierr )
    if (ierr /= 0) return
    call get_attr ( 'atomic-label', upf%psd )
    call get_attr ( 'atomic-number', upf%zmesh )
    call get_attr ( 'relativity', upf%rel )
    if ( upf%rel(1:5) == 'dirac' ) then
       upf%rel='full'
       upf%has_so = .true.
    end if
    call get_attr ( 'core-corrections', cc )
    upf%nlcc = (cc == 'yes')
    call xmlr_opentag ( 'exchange-correlation', IERR = ierr )
    if (ierr /= 0) return
    call xmlr_opentag ( 'libxc-info', IERR = ierr )
    if (ierr /= 0) return
    call get_attr ( 'number-of-functionals', nxc )
    do n=1,nxc
       call xmlr_readtag ( 'functional', xc(n), IERR = ierr )
       if (ierr > 0) return
       call get_attr ( 'id', xc(n) )
    end do
    call xmlr_closetag ( ) ! libxc-info
    call xmlr_closetag ( ) ! exchange-correlation
    !
    upf%dft = libxc_to_qe (nxc, xc)
    !
    call xmlr_opentag ( 'valence-configuration', IERR = ierr )
    if (ierr /= 0) return
    call get_attr ( 'total-valence-charge', upf%zp )
    ! here just count the number of valence wavefunctions
    n=0
    do 
       call xmlr_readtag ( 'shell', cc, IERR = ierr )
       if (ierr /= -1) exit
       n = n+1
       !call get_attr ( 'n', ndum )
       !call get_attr ( 'l', cc )
       !call get_attr ( 'occupation', ndum )
    end do
    upf%nwfc = n
    call xmlr_closetag ( ) ! valence-configuration
    call xmlr_closetag ( ) ! pseudo-atom-spec
    !
    ierr = 0
    !  
  END SUBROUTINE read_psml_pseudo_atom_spec
  !
  SUBROUTINE read_psml_grid ( ierr )
    !
    INTEGER :: ierr, i, n
    CHARACTER(LEN=1) :: dum
    REAL(dp) :: step, delta, r0, r1, r2
    !
    call xmlr_opentag ( 'grid', IERR = ierr )
    if (ierr /= 0) return
    call get_attr ( 'npts', upf%mesh )
    allocate (upf%r(upf%mesh))
    call xmlr_readtag ( 'annotation', dum, IERR = ierr )
    if (ierr > 0) return
    call get_attr ( 'step',  step )
    call get_attr ( 'scale', r0 )
    call get_attr ( 'delta', delta )
    call xmlr_readtag ( 'grid-data', upf%r, IERR = ierr )
    if (ierr /= 0) return
    call xmlr_closetag ( ) ! grid
    ! compute, or more exactly guess, suitable dr/dx for integration
    allocate (upf%rab(upf%mesh))
    upf%rab(:) = 0.0_dp
    ! i  is the index of the original logarithmic grid
    ! n  is the index of the PSML grid
    i = 0
    n = 0
    r1= 0.0_dp
    do
       i = i + 1
       r2 = r0*exp((i-1)*step)
       if ( r2-r1 > delta) then
          ! points at distance > delta are those selected
          n = n+1
          ! check: is this the correct grid?
          ! write (100,*) n, i, upf%r(n+1), r2
          if ( abs(r2-upf%r(n+1)) > 1.e-8 ) then
             print *, 'Unknown grid !'
             return
          end if
          upf%rab(n+1) = (r2-r1)
          r1 = r2
       end if
       if ( n == upf%mesh -1 ) exit
    end do
    upf%rab(1) = delta
    !
  END SUBROUTINE read_psml_grid
  !
  SUBROUTINE read_psml_radialfunc ( tag, rho, ierr )
    !
    CHARACTER(len=*) :: tag
    REAL(dp), allocatable :: rho(:)
    INTEGER :: ierr
    INTEGER :: npt
    !
    call xmlr_opentag ( trim(tag), IERR = ierr )
    if (ierr /= 0) return
    call xmlr_opentag ( 'radfunc', IERR = ierr )
    if (ierr /= 0) return
    call xmlr_opentag ( 'data', IERR = ierr )
    if (ierr /= 0) return
    call get_attr ( 'npts', npt )
    ! may differ
    if ( npt > upf%mesh ) then
       ierr = 1
       return
    end if
    allocate ( rho(upf%mesh) )
    read (iun,*) rho(1:npt)
    if ( npt < upf%mesh ) rho(npt+1:) = 0.0_dp
    call xmlr_closetag ( ) ! data
    call xmlr_closetag ( ) ! radfunc
    call xmlr_closetag ( ) ! tag
    !
  END SUBROUTINE read_psml_radialfunc
  !
  SUBROUTINE read_psml_nonlocal_projectors ( ierr )
    !
    use upf_utils, only: spdf_to_l
    INTEGER :: ierr
    INTEGER :: n, nb, npt, ndum
    real(dp) :: ekb, j
    CHARACTER(len=1) :: spdf

    !
    call xmlr_opentag('nonlocal-projectors', IERR = ierr )
    if (ierr /= 0) return
    upf%nbeta = 0
    upf%kkbeta= 0
    nb=0
    do 
       call xmlr_opentag ( 'proj', IERR = ierr )
       call get_attr ( 'l', spdf )
       call get_attr ( 'j', j )
       ! call get_attr ( 'seq',ndum )
       call get_attr ( 'ekb',ekb )
       call xmlr_opentag ( 'radfunc' )
       if ( ierr == -10 .and. upf%nbeta == 0 ) then
          ! first scan of the file completed, file has been rewound:
          ! set number of projectors, allocate and read  arrays
          upf%nbeta = nb
          ALLOCATE (upf%els_beta(nb), &
               upf%lll(nb),      &
               upf%kbeta(nb),    &
               upf%rcut(nb),     &
               upf%rcutus(nb),   &
               upf%dion(nb,nb),  &
               upf%qqq(nb,nb)    )
          allocate (upf%beta(upf%mesh,nb))
          IF (upf%has_so) ALLOCATE( upf%jjj(nb) ) 
          upf%rcut(:)  = 0.0_dp
          upf%rcutus(:)= 0.0_dp
          upf%dion(:,:)= 0.0_dp
          upf%qqq(:,:) = 0.0_dp
          ! reset counter
          nb = 0
       end if
       nb = nb+1
       call xmlr_opentag ( 'data', IERR = ierr  )
       call get_attr ( 'npts', npt )
       if ( upf%nbeta > 0 ) then
          ! actual read is done here during the second scan
          read (iun,*) upf%beta(1:npt,nb)
          do n=1,upf%mesh
             if ( n > 1 .and. n <= npt+1 ) then
                upf%beta(n,nb) = upf%beta(n,nb) * upf%r(n)
             else
                upf%beta(n,nb) = 0.0_dp
             end if
          end do
          upf%dion(nb,nb) = ekb
          upf%els_beta(nb) = '*'//spdf
          upf%lll(nb) = spdf_to_l(spdf)
          if (upf%has_so) upf%jjj(nb) = j
          upf%kbeta(nb) = npt
       else
          ! set max length of projectors during the first scan
          upf%kkbeta = max ( upf%kkbeta, npt )
       end if
       call xmlr_closetag () ! data
       call xmlr_closetag () ! radfun
       call xmlr_closetag () ! proj
       if ( nb == upf%nbeta ) exit
    end do
    call xmlr_closetag ( ) ! nonlocal-projectors
    !
  END SUBROUTINE read_psml_nonlocal_projectors
  !
  SUBROUTINE read_psml_pseudo_wave_functions ( ierr )
    !
    use upf_utils, only : spdf_to_l
    INTEGER :: ierr
    INTEGER :: n, npt, ndum
    CHARACTER(len=1) :: spdf
    !
    call xmlr_opentag('pseudo-wave-functions', IERR = ierr )
    if (ierr /= 0) return
    allocate ( upf%chi(upf%mesh,upf%nwfc) )
    allocate ( upf%els(upf%nwfc), &
               upf%oc(upf%nwfc), &
               upf%lchi(upf%nwfc), &
               upf%nchi(upf%nwfc), &
               upf%rcut_chi(upf%nwfc), &
               upf%rcutus_chi(upf%nwfc), &
               upf%epseu(upf%nwfc) )
    upf%rcut_chi(:)  = 0.0_dp
    upf%rcutus_chi(:)= 0.0_dp
    upf%oc(upf%nwfc) = 0.0_dp
    IF ( upf%has_so ) THEN
       allocate ( upf%nn(upf%nwfc) )
       upf%nn(:) = 0       
       allocate ( upf%jchi(upf%nwfc) )
    END IF
    do n=1,upf%nwfc
       call xmlr_opentag ( 'pswf', IERR = ierr )
       if ( ierr /= 0 ) return
       call get_attr ( 'l', spdf )
       call get_attr ( 'n', upf%nchi(n) )
       if ( upf%has_so ) call get_attr ( 'j', upf%jchi(n) )
       call get_attr ( 'energy_level', upf%epseu(n) )
       upf%lchi(n) = spdf_to_l(spdf)
       write(upf%els(n),'(i1,a1)') upf%nchi(n), spdf
       call xmlr_opentag ( 'radfunc' )
       call xmlr_opentag ( 'data', IERR = ierr  )
       call get_attr ( 'npts', npt )
       if ( npt > upf%mesh ) then
          ierr = 1
          return
       end if
       read (iun,*) upf%chi(1:npt,n)
       if ( npt < upf%mesh ) upf%chi(npt+1:upf%mesh,n) = 0.0_dp
       call xmlr_closetag () ! data
       call xmlr_closetag () ! radfun
       call xmlr_closetag () ! pswf
    end do
    call xmlr_closetag ( ) ! pseudo-wave-functions
    !
  END SUBROUTINE read_psml_pseudo_wave_functions
  !
END subroutine read_psml

function libxc_to_qe (nxc, xc)
  integer :: nxc
  integer :: xc(nxc)
  character(len=25) :: libxc_to_qe
  !
  libxc_to_qe = 'Not Recognized'
  print *, 'nxc, nc = ', nxc,xc
  if ( nxc < 2 ) return
  if ( xc(1) == 1 .and. xc(2) == 9 ) then
     libxc_to_qe = 'SLA-PZ' ! Perdew-Zunger
  else if ( xc(1) == 1 .and. xc(2) == 12 ) then
     libxc_to_qe = 'SLA-PW' ! Perdew-Wang
  else if ( xc(1) == 101 .and. xc(2) == 130 ) then
     libxc_to_qe = 'SLA-PW-PBX-PBC' ! PBE
  else if ( xc(1) == 116 .and. xc(2) == 133 ) then
     libxc_to_qe = 'SLA-PW-PSX-PSC' ! PBESOL
  end if
  !
end function libxc_to_qe
