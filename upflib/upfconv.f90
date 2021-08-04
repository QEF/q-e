!
! Copyright (C) 2020 Quantum ESPRESSO Foundation
! This file is distributed under the terms of the
! GNU General Public License. See the file 'License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! Original code "upf2casino" from Simon Binnie, 2011
!---------------------------------------------------------------------
PROGRAM upfconv
  !---------------------------------------------------------------------
  !
  !     Pseudopotential conversion utility, can
  !     - convert from:
  !          UPF v.1 or v.2 containing "&" characters,
  !          old PWSCF norm-conserving and Ultrasoft formats
  !          Vanderbilt ultrasoft PP generation code format
  !          CPMD format (TYPE=NUMERIC, LOGARITHMIC, CAR, GOEDECKER)
  !       to:
  !         UPF v.2 clean, or UPF with schema (xml)
  !     - extract and write to separate files:
  !          wavefunctions
  !          projectors ("beta" functions)
  !          potential
  !          if available, core wavefunctions from GIPAW section
  !     - convert to CASINO tabulated format (obsolete?)
  !
  USE pseudo_types, ONLY : pseudo_upf, deallocate_pseudo_upf
  USE casino_pp,    ONLY : conv_upf2casino, write_casino_tab
  USE write_upf_new,ONLY : write_upf
  !
  IMPLICIT NONE
  TYPE(pseudo_upf) :: upf_in
  INTEGER :: prefix_len, nargs, i,j 
  CHARACTER(LEN=256) :: filein, fileout
  CHARACTER(LEN=2)   :: conversion=' '
  CHARACTER(LEN=5)   :: schema='none'
  !
  nargs = command_argument_count()
  IF ( nargs < 1 .OR. nargs > 2 ) THEN
     WRITE(*,*) 'Usage: upfconv -c|-u|-x|-e "pseudopotential file"'
     WRITE(*,*) '       upfconv -h'
     STOP
  END IF
  !
  CALL get_command_argument(1, conversion)
  !
  IF ( conversion == "-h" ) THEN
     WRITE(*,*) 'Converts a pseudopotential file to either upf v.2, xml or CASINO formats.'
     WRITE(*,*) 'Usage: upfconv -c|-u|-x|-e "pseudopotential file"'
     WRITE(*,*) 'Options: -c     convert to CASINO format. MAY OR MAY NOT WORK (LIKELY IT DOES NOT)'
     WRITE(*,*) 'Options: -c     Make sure that the local channel chosen in the CASINO pp file'
     WRITE(*,*) '                is what you expect'
     WRITE(*,*) '         -u     convert to upf v.2 format'
     WRITE(*,*) '         -x     convert to xml format'
     WRITE(*,*) '         -e     extract GIPAW core wavefunctions if available,'
     WRITE(*,*) '                wavefunctions, projectors, potential'
     WRITE(*,*) '         -h     print this message'
     WRITE(*,*) 'The pseudopotential file can be any of the following:'
     WRITE(*,*) '*.upf or *.UPF      upf (v.1 or 2)'
     WRITE(*,*) '*.vdb or *.van      Vanderbilt US pseudopotential code format'
     WRITE(*,*) '*.rrkj3 or *.RRKJ3  Old US pseudopotential format of atomic code'
     WRITE(*,*) '*.cpi or *.fhi      FHI/abinit formats'
     WRITE(*,*) '*.gth DOES NOT WORK Goedecker-Teter-Hutter NC pseudo format'
     WRITE(*,*) '*.cpmd              CPMD format (TYPE=NUMERIC, LOGARITHMIC, CAR, GOEDECKER)'
     WRITE(*,*) 'none of the above   Old PWSCF norm-conserving format'
     STOP
  END IF
  !
  CALL get_command_argument(2, filein)  
  IF ( INDEX(TRIM(filein),'.UPF' ) > 0) THEN 
     prefix_len = INDEX(TRIM(filein),'.UPF') - 1
  ELSE IF (INDEX(TRIM(filein),'.upf') > 0 ) THEN
     prefix_len = INDEX(TRIM(filein),'.upf') - 1
  ELSE 
     prefix_len = LEN_TRIM(filein)
  ENDIF

  IF ( conversion == "-c" ) THEN
     fileout = filein(1:prefix_len) //'.out'
     WRITE(*,*) 'UPF to CASINO conversion'
     WRITE(*,*) 'All pseudopotential files generated should be &
          &thoroughly checked.'
     WRITE(*,*) 'In particular make sure the local channel chosen&
             & in the CASINO pp file is what you expected.'
  ELSE IF ( conversion == "-x" ) THEN
     fileout = filein(1:prefix_len) //'.xml'
     WRITE(*,*) 'UPF to xml format conversion'
  ELSE IF ( conversion == "-u" ) THEN
     fileout = filein(1:prefix_len) //'.UPF2'
     WRITE(*,*) 'UPF v.1 to UPF v.2 format conversion'
  ELSE IF ( conversion == "-e" .or. conversion == "-E" ) THEN
     fileout = filein(1:prefix_len)
  ELSE
     WRITE(*,*) 'Invalid option ' // conversion
     WRITE(*,*) 'Usage: upfconv -c|-u|-x|-e|-E "pseudopotential file"'
     STOP
  END IF
  IF ( prefix_len < 1 ) THEN
     WRITE(*,*) 'Empty file name, stopping'
     STOP
  END IF
  WRITE(*,*) 'input file: ' // trim(filein), ', output file: ' // trim(fileout)
 
  CALL read_ps ( filein, upf_in )

  IF ( conversion == "-c" ) THEN
     !
     CALL conv_upf2casino(upf_in)
     CALL write_casino_tab(upf_in, fileout)
     !
  ELSE IF ( conversion == "-e" ) THEN
     !
     CALL write_files ( upf_in, fileout )
     !
  ELSE
     !
     CALL conv_upf2xml(upf_in)
     IF ( conversion == "-x" ) THEN
        schema = 'qe_pp'
     ELSE IF ( conversion == "-u" ) THEN
        schema = 'v2'
     END IF
     CALL write_upf (FILENAME = fileout, UPF = upf_in, SCHEMA = schema)
     !
  ENDIF
  
  STOP
END PROGRAM upfconv

SUBROUTINE write_files( upf_in, fileout )
  !
  USE upf_const, ONLY : fpi
  USE pseudo_types, ONLY : pseudo_upf
  !
  IMPLICIT NONE
  TYPE(pseudo_upf), INTENT(inout) :: upf_in
  CHARACTER (LEN=*) :: fileout
  INTEGER :: i, j, n, ios, iunps
  !
  iunps=999
  !
  OPEN ( UNIT=iunps, FILE=TRIM(fileout)//'.wfc', STATUS='unknown', &
       FORM='formatted', IOSTAT=ios)
  IF ( ios /= 0 ) THEN
     WRITE(*,"('cannot write wfc file, stopping')" )
     STOP
  END IF
  WRITE(*,"('writing: ',a)") TRIM(fileout)//'.wfc'
  DO n=1,upf_in%mesh
     WRITE(iunps,'(30f12.6)') upf_in%r(n), (upf_in%chi(n,j), j=1,upf_in%nwfc)
  ENDDO
  CLOSE(iunps)

  OPEN ( UNIT=iunps, FILE=TRIM(fileout)//'.beta', STATUS='unknown', &
       FORM='formatted', IOSTAT=ios)
  IF ( ios /= 0 ) THEN
     WRITE(*,"('cannot write beta file, stopping')" )
     STOP
  END IF
  WRITE(*,"('writing: ',a)") TRIM(fileout)//'.beta'
  DO n=1,upf_in%mesh
     WRITE(iunps,'(30f12.6)') upf_in%r(n), (upf_in%beta(n,j), j=1,upf_in%nbeta)
  ENDDO
  CLOSE(iunps)

  OPEN ( UNIT=iunps, FILE=TRIM(fileout)//'.pot', STATUS='unknown', &
       FORM='formatted', IOSTAT=ios)
  IF ( ios /= 0 ) THEN
     WRITE(*,"('cannot write beta file, stopping')" )
     STOP
  END IF
  WRITE(*,"('writing: ',a)") TRIM(fileout)//'.pot'
  DO n=1,upf_in%mesh
     WRITE(iunps,'(4f12.6)') upf_in%r(n), upf_in%vloc(n), &
          upf_in%rho_at(n), upf_in%rho_atc(n)*fpi*upf_in%r(n)**2
  ENDDO
  CLOSE(iunps)

  IF(upf_in%has_gipaw) THEN 
     DO j = 1, upf_in%gipaw_ncore_orbitals
        OPEN(unit=iunps, file = TRIM(fileout)//TRIM(upf_in%gipaw_core_orbital_el(j))//".out")
        WRITE(*,"('writing: ',a)") TRIM(fileout)//TRIM(upf_in%gipaw_core_orbital_el(j))//".out"
        DO n = 1, upf_in%mesh
           WRITE(iunps,*) upf_in%r(n), upf_in%gipaw_core_orbital(n,j)
           WRITE(iunps+1,*) upf_in%r(n), upf_in%gipaw_core_orbital(n,j)
        ENDDO
        CLOSE(iunps)
     ENDDO

     ! write the same, but all in one file for xspectra
     OPEN(unit=iunps, file = TRIM(fileout)//".xspectra")
     WRITE(*,"('writing: ',a)") TRIM(fileout)//".xspectra"
     WRITE(iunps,"('# mesh size',i6,'; core orbitals:',99(a3,', '))")  upf_in%mesh, upf_in%gipaw_core_orbital_el(:)
     DO j = 1, upf_in%gipaw_ncore_orbitals
         DO n = 1, upf_in%mesh
           WRITE(iunps,*) upf_in%r(n), upf_in%gipaw_core_orbital(n,j)
        ENDDO
     ENDDO
     CLOSE(iunps)

  ELSE
     WRITE(*,*) "Core charge not written: this pseudopotential does not contain gipaw data"
  ENDIF

END SUBROUTINE write_files

SUBROUTINE conv_upf2xml( upf )
  !
  USE pseudo_types, ONLY : pseudo_upf
  USE upf_utils,    ONLY : version_compare
  !
  IMPLICIT NONE
  TYPE(pseudo_upf), INTENT(inout) :: upf
  INTEGER, EXTERNAL :: atomic_number
  !
  ! convert a few variables from UPF v.1 to UPF v.2/xml
  !
  IF ( version_compare(upf%nv,"2.0.1") == 'equal') RETURN
  upf%nv="2.0.1"
  !
  IF ( .NOT. ALLOCATED(upf%nchi) ) THEN
     ALLOCATE(upf%nchi(upf%nwfc))
     upf%nchi(:) = 0
  END IF
  IF ( .NOT. ALLOCATED(upf%rcut_chi) ) THEN
     ALLOCATE(upf%rcut_chi(upf%nwfc))
     upf%rcut_chi(:) = upf%rcut(:)
  END IF
  IF ( .NOT. ALLOCATED(upf%rcutus_chi) ) THEN
     ALLOCATE(upf%rcutus_chi(upf%nwfc))
     upf%rcutus_chi(:) = upf%rcutus(:)
  END IF
  IF ( .NOT. ALLOCATED(upf%epseu) ) THEN
     ALLOCATE(upf%epseu(upf%nwfc))
     upf%epseu(:) = 0.0
  END IF
  IF ( TRIM(upf%rel) == '' ) THEN 
     IF (upf%has_so) THEN
        upf%rel="full"
     ELSE IF ( upf%zmesh > 18 ) THEN
        upf%rel="scalar"
     ELSE
        upf%rel="no"
     ENDIF
  ENDIF
  !
  IF ( .not. upf%has_so) THEN
     upf%rmax = upf%r(upf%mesh)
     upf%zmesh = atomic_number( upf%psd )
     IF (upf%r(1) .GT. 1.d-16) THEN 
        upf%dx = log(upf%rmax/upf%r(1))/(upf%mesh-1)
        upf%xmin = log(upf%r(1)*upf%zmesh )
     END IF
  END IF
  !
END SUBROUTINE conv_upf2xml
