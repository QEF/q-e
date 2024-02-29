!
! Copyright (C) 2020-2023 Quantum ESPRESSO Foundation
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
  !          norm-conserving PPs in .psml format (experimental),
  !          old PWSCF norm-conserving and Ultrasoft formats,
  !          Vanderbilt ultrasoft PP generation code format,
  !          CPMD format (TYPE=NUMERIC, LOGARITHMIC, CAR, GOEDECKER)
  !       to:
  !         UPF v.2, or new UPF with schema (xml)
  !     - extract and write to separate files:
  !          wavefunctions
  !          projectors ("beta" functions)
  !          potential
  !          if available, core wavefunctions from GIPAW section
  !     - convert to CASINO tabulated format (obsolete?)
  !
  USE pseudo_types, ONLY : pseudo_upf, reset_upf, deallocate_pseudo_upf
  USE casino_pp,    ONLY : conv_upf2casino, write_casino_tab
  USE write_upf_new,ONLY : write_upf
  !
  IMPLICIT NONE
  TYPE(pseudo_upf) :: upf
  INTEGER :: prefix_len, nargs, i,j, ierr 
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
     WRITE(*,*) '         -p     write Vloc, betas in plottable format to file,'
     WRITE(*,*) '                both on the radial grid and Bessel-transformed'
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
  ELSE IF (INDEX(TRIM(filein),'.psml') > 0 ) THEN
     prefix_len = INDEX(TRIM(filein),'.psml') - 1
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
     WRITE(*,*) 'conversion to xml format'
  ELSE IF ( conversion == "-u" ) THEN
     fileout = filein(1:prefix_len) //'.UPF2'
     WRITE(*,*) 'conversion to UPF v.2'
  ELSE IF ( conversion == "-e" .or. conversion == "-p" ) THEN
     fileout = filein(1:prefix_len)
  ELSE
     WRITE(*,*) 'Invalid option ' // conversion
     WRITE(*,*) 'Usage: upfconv -c|-u|-x|-e|-p "pseudopotential file"'
     STOP
  END IF
  IF ( prefix_len < 1 ) THEN
     WRITE(*,*) 'Empty file, stopping'
     STOP
  END IF
  WRITE(*,*) 'input file: ' // trim(filein), ', output file: ' // trim(fileout)
  !
  CALL reset_upf( upf )
  !
  CALL read_ps_new ( filein, upf, .false., ierr )
  IF ( ierr > 0  ) THEN
     WRITE(*,*) 'Cannot read file, stopping'
     STOP
  END IF

  IF ( conversion == "-c" ) THEN
     !
     CALL conv_upf2casino(upf)
     CALL write_casino_tab(upf, fileout)
     !
  ELSE IF ( conversion == "-e" ) THEN
     !
     CALL write_files ( upf, fileout )
     !
  ELSE IF ( conversion == "-p" ) THEN
     !
     CALL plot_ps_loc ( upf, fileout )
     CALL plot_ps_beta( upf, fileout )
     !
  ELSE
     !
     CALL conv_upf2xml(upf)
     IF ( conversion == "-x" ) THEN
        schema = 'qe_pp'
     ELSE IF ( conversion == "-u" ) THEN
        schema = 'v2'
     END IF
     CALL write_upf (FILENAME = fileout, UPF = upf, SCHEMA = schema)
     !
  ENDIF
  !
  CALL deallocate_pseudo_upf (upf)  
  !
  STOP
END PROGRAM upfconv

SUBROUTINE write_files( upf, fileout )
  !
  USE upf_const, ONLY : fpi
  USE pseudo_types, ONLY : pseudo_upf
  !
  IMPLICIT NONE
  TYPE(pseudo_upf), INTENT(inout) :: upf
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
  DO n=1,upf%mesh
     WRITE(iunps,'(30f12.6)') upf%r(n), (upf%chi(n,j), j=1,upf%nwfc)
  ENDDO
  CLOSE(iunps)

  OPEN ( UNIT=iunps, FILE=TRIM(fileout)//'.beta', STATUS='unknown', &
       FORM='formatted', IOSTAT=ios)
  IF ( ios /= 0 ) THEN
     WRITE(*,"('cannot write beta file, stopping')" )
     STOP
  END IF
  WRITE(*,"('writing: ',a)") TRIM(fileout)//'.beta'
  DO n=1,upf%mesh
     WRITE(iunps,'(30f12.6)') upf%r(n), (upf%beta(n,j), j=1,upf%nbeta)
  ENDDO
  CLOSE(iunps)

  OPEN ( UNIT=iunps, FILE=TRIM(fileout)//'.pot', STATUS='unknown', &
       FORM='formatted', IOSTAT=ios)
  IF ( ios /= 0 ) THEN
     WRITE(*,"('cannot write beta file, stopping')" )
     STOP
  END IF
  WRITE(*,"('writing: ',a)") TRIM(fileout)//'.pot'
  DO n=1,upf%mesh
     WRITE(iunps,'(4f12.6)') upf%r(n), upf%vloc(n), &
          upf%rho_at(n), upf%rho_atc(n)*fpi*upf%r(n)**2
  ENDDO
  CLOSE(iunps)

  IF(upf%has_gipaw) THEN 
     DO j = 1, upf%gipaw_ncore_orbitals
        OPEN(unit=iunps, file = TRIM(fileout)//TRIM(upf%gipaw_core_orbital_el(j))//".out")
        WRITE(*,"('writing: ',a)") TRIM(fileout)//TRIM(upf%gipaw_core_orbital_el(j))//".out"
        DO n = 1, upf%mesh
           WRITE(iunps,*) upf%r(n), upf%gipaw_core_orbital(n,j)
           WRITE(iunps+1,*) upf%r(n), upf%gipaw_core_orbital(n,j)
        ENDDO
        CLOSE(iunps)
     ENDDO

     ! write the same, but all in one file for xspectra
     OPEN(unit=iunps, file = TRIM(fileout)//".xspectra")
     WRITE(*,"('writing: ',a)") TRIM(fileout)//".xspectra"
     WRITE(iunps,"('# mesh size',i6,'; core orbitals:',99(a3,', '))")  upf%mesh, upf%gipaw_core_orbital_el(:)
     DO j = 1, upf%gipaw_ncore_orbitals
         DO n = 1, upf%mesh
           WRITE(iunps,*) upf%r(n), upf%gipaw_core_orbital(n,j)
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
  IF ( upf%nwfc > 0 ) THEN
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
  !
SUBROUTINE plot_ps_loc ( upf, fileout )
  !
  !  Plots Vloc in both real and reciprocal space and in plottable format,
  !  for testing purposes.
  !        File                    contains
  !  "fileout"-vofr.dat         r, Vloc(r), Vloc(r)+Ze^2/r
  !  "fileout"-vofq.dat         q, q^2*Vloc(q)
  !  where f(q) = \int f(r) j_l(qr) r^2 dr
  !
  USE pseudo_types, ONLY : pseudo_upf
  !
  USE upf_kinds,  ONLY : dp
  USE upf_const,  ONLY : fpi, e2
  !
  IMPLICIT NONE
  !
  TYPE(pseudo_upf), intent(in) :: upf
  CHARACTER(LEN=256) :: fileout
  !
  REAL(DP), parameter :: dq = 0.01_dp
  INTEGER, PARAMETER :: nq = 4000
  REAL(DP) :: q(nq), vq(nq), vq2(nq), aux(upf%mesh), aux1(upf%mesh)
  REAL(DP) :: vq0, vq1, vlq
  INTEGER :: iq, ir, msh
  !
  open( unit = 100, file=trim(fileout) // '-vofr.dat', status='unknown', &
       form='formatted' )
  write(100,'("#   r,  vloc(r), vloc+Ze^2/r")') 
  do ir=1,upf%mesh
     write(100,'(f10.4,2f25.12)') upf%r(ir), upf%vloc(ir), &
          upf%vloc(ir) + upf%zp*e2/max(0.1_dp,upf%r(ir))
  end do
  close(unit=100)

  do iq = 1, nq
     q(iq) = (iq-1)*dq
  end do
  !
  ! Compute vloc(q=0) for the local potential without the divergence
  !
  do ir = 1, upf%mesh
     aux (ir) = upf%r(ir) * (upf%r(ir)*upf%vloc(ir) + upf%zp * e2)
     if ( upf%r(ir) < 10 ) msh = ir
  enddo
  !
  print '("msh, mesh = ",2i5)',msh,upf%mesh
  call simpson (upf%mesh, aux, upf%rab, vq0)
  call simpson (msh, aux, upf%rab, vq1)
  !
  ! Compute q^2*vloc(q). We first store the part of the integrand 
  !   function independent of q in real space
  !
  do ir = 1, upf%mesh
     aux1 (ir) = upf%r(ir) * upf%vloc(ir) + upf%zp * e2 * erf (upf%r(ir) )
  enddo
  !
  !    and here we perform the integral, after multiplication
  !    with the q-dependent part
  !
  do iq = 1, nq
     do ir = 1, upf%mesh
        aux(ir) = aux1(ir) * sin (q(iq)*upf%r(ir)) * q(iq)
     enddo
     call simpson (upf%mesh, aux, upf%rab, vlq)
     !
     !   here we re-add the analytic fourier transform of the erf function
     !
     vq(iq) = vlq  - e2 * upf%zp * exp ( - q(iq)**2/4.0_dp) !!! / q(iq)**2
     !
     call simpson (msh, aux, upf%rab, vlq)
     vq2(iq)= vlq  - e2 * upf%zp * exp ( - q(iq)**2/4.0_dp) !!! / q(iq)**2
  enddo
  !!! vq(:) = vq(:) * fpi / omega
  !
  open( unit = 100, file=trim(fileout) // '-vofq.dat', status='unknown', &
       form='formatted' )
  write(100,'("#   q,  q^2*v(q) (mesh)  q^2*v(q) (msh)  v(0) = ",2f15.8)') &
       vq0,vq1
  do iq=1,nq
     write(100,'(f10.4,2f25.12)') q(iq), vq(iq), vq2(iq)
  end do
  close(unit=100)
  !  
END SUBROUTINE plot_ps_loc
!
SUBROUTINE plot_ps_beta ( upf, fileout )
  !
  !  Plots beta projectors in both real and reciprocal space
  !  and in plottable format for testing purposes.
  !        File                    contains
  ! "fileout"-betar-n.dat        r, beta_n(r)
  ! "fileout"-betaq-n.dat        q, q^2 beta_n(q)
  ! where f(q) = \int f(r) j_l(qr) r^2 dr
  !
  USE pseudo_types, ONLY : pseudo_upf
  !
  USE upf_kinds,  ONLY : dp
  USE upf_const,  ONLY : fpi, e2
  !
  IMPLICIT NONE
  !
  TYPE(pseudo_upf), intent(in) :: upf
  CHARACTER(LEN=256) :: fileout
  !
  REAL(DP) :: dq = 0.01_dp
  INTEGER, PARAMETER :: nq = 4000
  REAL(DP) :: q(nq), aux(upf%mesh), aux1(upf%mesh)
  REAL(DP) :: vq0, vq1
  INTEGER :: iq, ir, nb, iun
  character(len=256) :: filename
  character(len=2) :: number
  !
  print '("kkbeta = ",1i5)',upf%kkbeta
  do nb=1,upf%nbeta
     if ( nb < 10 ) then
        write(filename,'(A,"-betar-",i1,".dat")') trim(fileout),nb
     else
        write(filename,'(A,"-betar-",i2,".dat")') trim(fileout),nb
     end if
     iun = 199 + nb
     open( unit=iun, file=filename, status='unknown', form='formatted' )
     write(iun,'("#   r,  vbeta(r,lm)")') 
     do ir=1,upf%kkbeta
        write(iun,'(f10.4,f25.12)') upf%r(ir), upf%beta(ir,nb)
     end do
     close(unit=iun)
  end do
  
  do iq = 1, nq
     q(iq) = (iq-1)*dq
  end do
  !
  ! Compute q^2*beta(q)
  !
  do nb=1,upf%nbeta
     if ( nb < 10 ) then
        write(filename,'(A,"-betaq-",i1,".dat")') trim(fileout),nb
     else
        write(filename,'(A,"-betaq-",i2,".dat")') trim(fileout),nb
     end if
     iun = 299 + nb
     open( unit=iun, file=filename, status='unknown', form='formatted' )
     write(iun,'("#   q,  q^2*v(q) ((old)  a^2*v(q) (cp90)")')
     do iq = 1, nq
        do ir = 1, upf%kkbeta
           aux(ir) = upf%beta(ir,nb) * sin (q(iq)*upf%r(ir)) * q(iq)
        enddo
        call simpson (upf%kkbeta, aux, upf%rab, vq0)
        call simpson_cp90 (upf%kkbeta, aux, upf%rab, vq1)
        write(iun,'(f10.4,2f25.12)') q(iq), vq0, vq1
     end do
     close(unit=iun)
  enddo
  !
END SUBROUTINE plot_ps_beta
