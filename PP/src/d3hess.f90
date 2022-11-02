!
! Copyright (C) 2001-2022 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
! Author: Ivan Carnimeo (October 2022)
!----------------------------------------------------------------------------
!
!----------------------------------------------------------------------------
program d3hess
  USE io_global,        ONLY: stdout, ionode, ionode_id
  USE io_files,         ONLY: prefix, tmp_dir
  USE kinds,            ONLY: DP
  USE mp,               ONLY: mp_bcast
  USE mp_global,        ONLY: mp_startup
  USE mp_world,         ONLY: world_comm
  USE environment,      ONLY: environment_start, environment_end
  !
  USE cell_base,        ONLY: alat, at
  USE ions_base,        ONLY: nat, tau, ityp, atm
  USE input_parameters, ONLY: dftd3_version, dftd3_threebody
  USE funct,            ONLY: get_dft_short
  USE dftd3_api,        ONLY: dftd3_init, dftd3_set_functional, get_atomic_number, dftd3_pbc_dispersion
  USE dftd3_qe,         ONLY: dftd3_xc, dftd3, dftd3_pbc_gdisp, dftd3_pbc_hdisp, dftd3_in
  !
  IMPLICIT NONE
  INTEGER :: ios
  CHARACTER(LEN=256), EXTERNAL :: trimcheck
  CHARACTER(len=256) :: filhess, outdir
  REAL(DP) :: step
  LOGICAL :: needwf = .FALSE.
  LOGICAL :: q_gamma
  !
  INTEGER :: iat, jat, ixyz, jxyz, irep, jrep, krep, i,j
  INTEGER :: nnat, nrep, nhess, nsize
  INTEGER :: rep_cn(3), rep_vdw(3), rep_hes(3)
  REAL(DP) :: latvecs(3,3)
  REAL(DP) :: stress_d3(3,3)
  CHARACTER(LEN=256):: dft_ , formt
  INTEGER, ALLOCATABLE  :: atnum(:)
  REAL(DP), ALLOCATABLE :: xyz(:,:), buffer(:)
  REAL(DP), ALLOCATABLE :: force_d3(:,:), hess_d3(:,:,:,:,:,:,:)
  !
  NAMELIST /input/ prefix, outdir, step, q_gamma, filhess
  !
9078 FORMAT( '     DFT-D3 Dispersion         =',F17.8,' Ry' )
9035 FORMAT(5X,'atom ',I4,' type ',I2,'   force = ',3F14.8)
  !
  write(stdout,'(A)') 'PROGRAM: d3hess '
  !
  ! initialise environment
  !
#if defined(__MPI)
  CALL mp_startup ( )
#endif
  CALL environment_start ( 'd3hess' )
  !
  ios = 0
  !
  IF ( ionode ) THEN
     !
     !   set default values for variables in namelist
     !
     CALL get_environment_variable( 'ESPRESSO_TMPDIR', outdir )
     IF ( trim( outdir ) == ' ' ) outdir = './'
     prefix ='pwscf'
     filhess=' ' 
     q_gamma=.false. ! whether to use a much cheaper algorithm when q=0,0,0
     step=2.d-5      ! step for numerical differentiation
     !
     CALL input_from_file ( )
     READ (5,input,IOSTAT=ios)
     !
     tmp_dir = trimcheck (outdir)
     IF ( filhess == ' ' ) filhess = trim(prefix)//'.hess'
     filhess = TRIM(tmp_dir)//TRIM(filhess)
     !
  ENDIF
  !
  ! ... Broadcast variables
  !
  CALL mp_bcast(ios, ionode_id, world_comm)
  CALL errore('d3hess', 'reading input namelist', ABS(ios))

  CALL mp_bcast( tmp_dir, ionode_id, world_comm )
  CALL mp_bcast( prefix, ionode_id, world_comm )
  CALL mp_bcast( step, ionode_id, world_comm )
  CALL mp_bcast( q_gamma, ionode_id, world_comm )
  !
  CALL read_file_new ( needwf )
  !
  ! Setting DFT-D3 functional dependent parameters
  !
!civn 
dftd3_threebody = .false.
!
  if (dftd3_version==2) dftd3_threebody=.false.
  dftd3_in%threebody = dftd3_threebody
  CALL dftd3_init(dftd3, dftd3_in)
  dft_ = get_dft_short( )
  dft_ = dftd3_xc ( dft_ )
  CALL dftd3_set_functional(dftd3, func=dft_, version=dftd3_version,tz=.false.)
  WRITE( stdout, '(/,5x,A,f24.12)') 'Differentiation step: ',  step 
!civn 
write(*,'(/,A,/)') '!!!WARNING: FIX DFT-D3 XML FILE READING!!!!'
  WRITE( stdout, '(5x,A,3I4)') 'DFT-D3 version: ',  dftd3_version  
  WRITE( stdout, '(5x,A,L)') 'DFT-D3 threebody: ',  dftd3_threebody  
  IF(q_gamma) THEN
    WRITE( stdout, '(5x,A)') 'Using a cheap algorithm for q=0,0,0 only'
  ELSE
    WRITE( stdout, '(5x,A)') 'Using a general (and memory consuming) algorithm for any q '
  END IF
  !
  ! Computing DFT-D3 forces to get rep_vdw and check consistency with the scf 
  !
  CALL start_clock('force_dftd3')
  ALLOCATE( xyz(3,nat), atnum(nat), force_d3(3,nat) )
  force_d3(:,:) = 0.0_DP
  latvecs(:,:)=at(:,:)*alat
  xyz(:,:)=tau(:,:)*alat
  DO iat = 1, nat
     atnum(iat) = get_atomic_number(TRIM(atm(ityp(iat))))
  ENDDO
  !
  call dftd3_pbc_gdisp(dftd3, xyz, atnum, latvecs, force_d3, stress_d3, rep_cn, rep_vdw)
  force_d3 = -2.d0*force_d3
  !
  WRITE( stdout, '(/,5x,A,3I4)') 'Number of images for CN: ',  rep_cn(:)
  WRITE( stdout, '(/,5x,A,3I4)') 'Number of images for VdW: ', rep_vdw(:)
  WRITE( stdout, '(/,5x,"DFT-D3 dispersion forces in the unit cell:")')
  DO iat = 1, nat
     WRITE( stdout, 9035) iat, ityp(iat), (force_d3(ixyz,iat), ixyz = 1, 3)
  ENDDO
  !
  CALL stop_clock('force_dftd3')
  !
  ! Computing DFT-D3 hessian 
  !
  IF(q_gamma) THEN
    rep_hes(:) = 0
  ELSE
    rep_hes(:) = rep_vdw(:) 
  END IF 
  !
  WRITE( stdout, '(/,5x,A,3I4)') 'Number of cells replicated along each semiaxis: ', rep_vdw(1), rep_vdw(2), rep_vdw(3)
  WRITE( stdout, '(/,5x,A,3I4)') 'Number of cells used for Hessian allocations: ',   rep_hes(1), rep_hes(2), rep_hes(3)
  nrep = (2*rep_hes(1)+1) * (2*rep_hes(2)+1) * (2*rep_hes(3)+1)       
  WRITE( stdout, '(5x,A,I9)') 'Number of cells in the supercell (Hessian): ', nrep 
  nnat = nat * nrep                                                   
  WRITE( stdout, '(5x,A,I9)') 'Number of atoms in the supercell (Hessian): ', nnat
  nhess = (3*nat)**2  * nrep 
  ! note that we are allocating one Hessian for each replicated cell (nhess), 
  ! that is much smaller than the full Hessian of the supercell ((3*nat*nrep)**2), 
  ! because we ultimately want the Hessian only for the unit cell
  WRITE( stdout, '(5x,A,I9)') 'Hessian allocation dimensions: ', nhess 
  !
  ALLOCATE( hess_d3(-rep_hes(1):rep_hes(1),-rep_hes(2):rep_hes(2),-rep_hes(3):rep_hes(3),3,nat,3,nat) )
  !
  nsize = int(size(hess_d3))
  IF( nsize .ne. nhess ) Call errore('d3hess', "Wrong Hessian dimensions", 1)
  !
  CALL start_clock('hessian_dftd3')
  !
  Call dftd3_pbc_hdisp(dftd3, stdout, step, xyz, atnum, latvecs, rep_cn, rep_vdw, hess_d3, q_gamma )
  !
  CALL stop_clock('hessian_dftd3')  
  !
  ! Writing DFT-D3 hessian on file
  !
  WRITE( stdout, '(/,5x,2A)') 'Writing Hessian on file: ',  filhess
  ! 
  WRITE(formt,'(A,I9,A)') '(', 3*nat, 'f24.16)'
  ALLOCATE( buffer(3*nat) )
  !
  OPEN (unit = 1, file = filhess, status = 'unknown')
  !
  WRITE(1,'(A)') 'Hessian matrix of the Grimme-D3 dispersion term'
  WRITE(1,'(A,5I9,3x,L)') 'System: ', rep_hes(:), nat, nnat, q_gamma
  !
  DO irep = -rep_hes(1), rep_hes(1)
    DO jrep = -rep_hes(2), rep_hes(2)
      DO krep = -rep_hes(3), rep_hes(3)
        !
        IF(irep.eq.0 .and. jrep.eq.0 .and. krep.eq.0) THEN
          WRITE(1, '(3(A,I4))') 'Unit cell: irep= ',irep, ' jrep= ',jrep, ' krep= ',krep
        ELSE
          WRITE(1, '(3(A,I4))') 'Replica cell: irep= ',irep, ' jrep= ',jrep, ' krep= ',krep
        END IF
        !
        DO i = 1, 3*nat         ! 1 2 3 4 5 6 7 8 9 ... 3*nat
          iat  = (i+2)/3        ! 1 1 1 2 2 2 3 3 3 ... nat 
          ixyz = i - 3* (iat-1) ! 1 2 3 1 2 3 1 2 3 ... 3 
          !
          DO j = 1, 3*nat
            jat  = (j+2)/3 
            jxyz = j - 3* (jat-1) 
            buffer(j) = hess_d3(irep,jrep,krep,ixyz,iat,jxyz,jat)
          END DO 
          !  
          WRITE(1, formt ) buffer(1:3*nat)
          !
        END DO 
        !
      END DO 
    END DO 
  END DO 
  !
  CLOSE (1)
  !
  DEALLOCATE( buffer, xyz, atnum, force_d3, hess_d3 )
  !
  WRITE( stdout, * ) 
  CALL print_clock('force_dftd3')
  CALL print_clock('hessian_dftd3')
  !
  CALL environment_end ( 'd3hess' )
  !
  CALL stop_pp
  !
  return
  !
end program
