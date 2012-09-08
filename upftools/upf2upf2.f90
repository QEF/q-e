!
! Copyright (C) 2011 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file 'License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!---------------------------------------------------------------------
PROGRAM upf2upf2
  !---------------------------------------------------------------------
  !
  !  Convert a pseudopotential written in UPF v.1 format to UPF v.2 format
  !
  USE pseudo_types, ONLY : pseudo_upf, nullify_pseudo_upf, &
                           deallocate_pseudo_upf
  USE radial_grids, ONLY: radial_grid_type, nullify_radial_grid
  USE read_upf_v1_module, ONLY : read_upf_v1
  USE write_upf_v2_module, ONLY: write_upf_v2
  !
  IMPLICIT NONE
  TYPE(pseudo_upf) :: upf
  TYPE (radial_grid_type), TARGET :: grid
  CHARACTER(len=256) filein, fileout
  INTEGER :: ios
  INTEGER, EXTERNAL :: atomic_number
  !
  CALL get_file ( filein )
  IF ( trim(filein) == ' ') &
       CALL errore ('upf2upf2', 'usage: upf2upf2 "file-to-be-converted"', 1)
  OPEN ( unit=1, file=filein, status = 'old', form='formatted', iostat=ios )
  IF ( ios /= 0) &
     CALL errore ('upf2upf2', 'file: '//trim(filein)//' not found', 2)
  !
  CALL nullify_pseudo_upf ( upf )
  CALL nullify_radial_grid ( grid )
  upf%grid => grid
  CALL read_upf_v1 (1, upf, grid, ios)
  IF ( ios /= 0) &
     CALL errore ('upf2upf2', 'file '//trim(filein)//' not UPF v.1', 3)
  !
  CLOSE (unit=1)
  !
  ! convert a few variables
  !
  upf%nv       = "2.0.1"
  IF ( .not. associated (upf%epseu) ) THEN
     ALLOCATE ( upf%epseu( upf%nwfc) )
     upf%epseu=0
  ENDIF
  ALLOCATE ( upf%nchi( upf%nwfc) )
  IF ( .not. associated(upf%nn) ) THEN
     upf%nchi=0
  ELSE
     upf%nchi=upf%nn(1:upf%nwfc)
  ENDIF
  ALLOCATE ( upf%rcut_chi( upf%nwfc ) )
  ALLOCATE ( upf%rcutus_chi( upf%nwfc ) )
  upf%rcut_chi=upf%rcut(1:upf%nwfc)
  upf%rcutus_chi=upf%rcutus(1:upf%nwfc)
  !
  upf%rmax = upf%r(upf%mesh)
  upf%dx = log(upf%rmax/upf%r(1))/(upf%mesh-1)
  upf%zmesh = atomic_number( upf%psd )
  upf%xmin = log(upf%r(1)*upf%zmesh )
  IF ( upf%has_so) THEN
     upf%rel="full"
  ELSEIF ( upf%zmesh > 18 ) THEN
     upf%rel="scalar"
  ELSE
     upf%rel="no"
  ENDIF
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

END PROGRAM upf2upf2

