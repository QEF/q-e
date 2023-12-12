!
! Copyright (C) 2023 Quantum ESPRESSO Foundation
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!---------------------------------------------------------------------
MODULE read_psml_module
  !---------------------------------------------------------------------
  !
  PUBLIC :: read_psml
  !
CONTAINS
  !--------------------------------------------------------
  subroutine read_psml ( filename, upf, ierr )
  !-----------------------------------------------------
  !! Read pseudopotential files in PSML format using "xmltools"
  !! stores data into the "upf" structure. Note that:
  !! - PSML uses a nonstandard radial grid that is not well suited 
  !!   for QE integration methods, so all variables are interpolated
  !!   to a uniform grid (grid parameters: dr=0.01, rmax=5 a.u.,
  !!   the latter being the typical rmax used in PSML)
  !! - in PSML, most arrays: local potential, projectors, charges,
  !!   may be shorter than the full length of the radial grid
  !!   (as specified in tag argument "npts"), while in the UPF 
  !!  format and upf% structure they all have the same upf%mesh size
  !! To get QE variables from PSML ones:
  !! - PSML pseudo-core charge must be divided by 4\pi (why? don't know)
  !! - PSML projectors must be multiplied by r (why? don't know)
  !! - PSML potentials must be multiplied by e^2=2 to bring them to Ry
  !! Tested only for a small  subset of PSML files
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
  TYPE(pseudo_upf), INTENT(INOUT) :: upf
  !! the derived type storing the pseudo data
  !! INOUT because many variables are reset to default values in input
  INTEGER, INTENT(OUT) :: ierr
  !! error code (0 if correctly read)
  CHARACTER(len=30) :: tag
  !! tag where error (ierr != 0) was detected
  INTEGER :: iun
  !! unit for reading data
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
  IF ( ierr /= 0 ) GO TO 10
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
  call read_psml_pseudo_atom_spec ( tag, ierr )
  IF ( ierr /= 0 ) GO TO 10
  !
  tag = 'grid'
  call read_psml_grid ( ierr )
  IF ( ierr /= 0 ) GO TO 10
  !
  tag = 'valence-charge'
  call read_psml_radialfunc ( tag, upf%rho_at, ierr )
  IF ( ierr /= 0 ) GO TO 10
  !
  IF ( upf%nlcc ) THEN
     tag = 'pseudocore-charge'
     call read_psml_radialfunc ( tag, upf%rho_atc, ierr )
     IF ( ierr /= 0 ) GO TO 10
     upf%rho_atc(:) = upf%rho_atc(:) / fpi
  END IF
  !
  tag = 'local-potential'
  call read_psml_radialfunc ( tag, upf%vloc, ierr )
  IF ( ierr /= 0 ) GO TO 10
  !
  tag = 'nonlocal-projectors'
  call read_psml_nonlocal_projectors ( ierr )
  IF ( ierr /= 0 ) GO TO 10
  !
  tag = 'pseudo-wave-functions'
  call read_psml_pseudo_wave_functions ( ierr )
  IF ( ierr /= 0 ) THEN
     !! optional tag, may or may not be present
     print *, 'read_psml: tag ',trim(tag),' not present'
     upf%nwfc = 0 
     ierr = 0
  END IF
  !
  call xmlr_closetag ( ) ! psml
  call xml_closefile( )
  !
  ! upf%r now contains the original grid, upf%rab the uniform grid 
  ! all arrays are read in the former and interpolated on the latter
  ! now set the correct grid and grid derivative (upf%rab = (dr/dx)*dx)
  !
  deallocate(upf%r)
  allocate(upf%r(upf%mesh))
  upf%r = upf%rab
  upf%rab = upf%dx
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
  SUBROUTINE read_psml_pseudo_atom_spec ( tag, ierr )
    !
    CHARACTER(len=*), INTENT(inout) :: tag
    INTEGER, INTENT(out)  :: ierr
    INTEGER :: n, nxc, ndum
    INTEGER :: xc(6)
    CHARACTER(len=3) :: cc
    
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
    tag = 'exchange-correlation'
    call xmlr_opentag ( tag, IERR = ierr )
    if (ierr /= 0) return
    tag = 'libxc-info'
    call xmlr_opentag ( tag, IERR = ierr )
    if (ierr /= 0) return
    call get_attr ( 'number-of-functionals', nxc )
    do n=1,nxc
       tag = 'functional'
       call xmlr_readtag ( tag, xc(n), IERR = ierr )
       if (ierr > 0) return
       call get_attr ( 'id', xc(n) )
    end do
    call xmlr_closetag ( ) ! libxc-info
    call xmlr_closetag ( ) ! exchange-correlation
    !
    upf%dft = libxc_to_qe (nxc, xc)
    !
    tag = 'valence-configuration'
    call xmlr_opentag ( tag, IERR = ierr )
    if (ierr /= 0) return
    tag = 'total-valence-charge'
    call get_attr ( tag, upf%zp )
    ! here just count the number of valence wavefunctions
    n=0
    do 
       tag = 'shell'
       call xmlr_readtag ( tag, cc, IERR = ierr )
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
    INTEGER, INTENT(OUT) :: ierr
    INTEGER :: npt, n
    REAL(dp):: step, r0, delta
    CHARACTER(LEN=1) :: dum
    !
    call xmlr_opentag ( 'grid', IERR = ierr )
    if (ierr /= 0) return
    call get_attr ( 'npts', npt )
    allocate (upf%r(npt))
    call xmlr_readtag ( 'annotation', dum, IERR = ierr )
    if (ierr > 0) return
    call get_attr ( 'step',  step )
    call get_attr ( 'scale', r0 )
    call get_attr ( 'delta', delta )
    call xmlr_readtag ( 'grid-data', upf%r, IERR = ierr )
    if (ierr /= 0) return
    call xmlr_closetag ( ) ! grid
    ! Now store in upf%rab a uniform grid with dx=0.01 up to rmax=5 
    upf%dx = 0.01_dp
    upf%rmax=5.0_dp
    upf%mesh = upf%rmax/upf%dx+1
    allocate (upf%rab(upf%mesh))
    do n = 1, upf%mesh
       upf%rab(n) = (n-1)*upf%dx
    end do
    !
  END SUBROUTINE read_psml_grid
  !
  SUBROUTINE read_psml_radialfunc ( tag, rho, ierr )
    !
    USE splinelib,  ONLY : dosplineint
    CHARACTER(len=*) :: tag
    REAL(dp), allocatable :: rho(:)
    REAL(DP), allocatable :: rint(:)
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
    if ( npt > size(upf%r) ) then
       ierr = 1
       return
    end if
    allocate ( rint(npt) )
    allocate ( rho (upf%mesh) )
    read (iun,*) rint
    call dosplineint( upf%r(1:npt), rint, upf%rab, rho )
    call xmlr_closetag ( ) ! data
    call xmlr_closetag ( ) ! radfunc
    call xmlr_closetag ( ) ! tag
    deallocate (rint)
    !
  END SUBROUTINE read_psml_radialfunc
  !
  SUBROUTINE read_psml_nonlocal_projectors ( ierr )
    !
    USE splinelib,  ONLY : dosplineint
    USE upf_utils, only: spdf_to_l
    INTEGER :: ierr
    INTEGER :: n, nb, npt, ndum
    REAL(DP), allocatable :: betaint(:)
    REAL(dp) :: ekb, j
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
          allocate (betaint(npt))
          ! actual read is done here during the second scan
          read (iun,*) betaint
          call dosplineint( upf%r(1:npt), betaint, upf%rab, upf%beta(:,nb))
          do n=1,upf%mesh
             upf%beta(n,nb) = upf%beta(n,nb) * upf%rab(n)
          end do
          upf%dion(nb,nb) = ekb
          upf%els_beta(nb) = '*'//spdf
          upf%lll(nb) = spdf_to_l(spdf)
          if (upf%has_so) upf%jjj(nb) = j
          upf%kbeta(nb) = npt
          deallocate (betaint)
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
    USE upf_utils, only : spdf_to_l
    USE splinelib,  ONLY : dosplineint
    INTEGER :: ierr
    INTEGER :: n, npt, ndum
    REAL(dp), ALLOCATABLE :: chint(:)
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
    IF ( upf%has_so ) allocate ( upf%jchi(upf%nwfc) )
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
       if ( npt > size(upf%r) ) then
          ierr = 1
          return
       end if
       allocate ( chint(npt) )
       read (iun,*) chint
       call dosplineint( upf%r(1:npt), chint, upf%rab, upf%chi(:,n))
       deallocate ( chint) 
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
  ! print *, 'nxc, nc = ', nxc,xc
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

END MODULE read_psml_module
