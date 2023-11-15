!
! Copyright (C) 2001-2022 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
MODULE d3hess_mod
  !
  USE kinds,            ONLY: DP
  !
  IMPLICIT NONE
  !
  SAVE
  !
  ! whether to use a much cheaper algorithm when q=0,0,0
  LOGICAL :: q_gamma = .FALSE.
  ! whether to check consistency between hessian, forces and energies
  LOGICAL :: debug = .FALSE.
  ! step for numerical differentiation
  REAL(DP) :: step = 2.d-5
  !
  CHARACTER(LEN=14), PARAMETER :: AUTOMATIC_NAME = 'automatic.hess'
  !
  CONTAINS
!
SUBROUTINE d3hess_sub(filhess)
  !
  USE io_global,        ONLY: stdout, ionode
  !
  USE cell_base,        ONLY: alat, at, bg
  USE ions_base,        ONLY: nat, tau, ityp, atm
  USE dftd3_api,        ONLY: get_atomic_number
  USE dftd3_qe,         ONLY: dftd3, dftd3_pbc_gdisp_new, dftd3_pbc_hdisp
  !
  IMPLICIT NONE
  !
  CHARACTER(len=*), INTENT(IN) :: filhess
  !
  INTEGER :: ios
  !
  INTEGER :: iat, jat, ixyz, jxyz, irep, jrep, krep, i,j
  INTEGER :: nnat, nrep, nhess, nsize
  INTEGER :: rep_cn(3), rep_vdw(3), rep_hes(3)
  REAL(DP) :: latvecs(3,3)
  CHARACTER(LEN=256) :: formt
  INTEGER, ALLOCATABLE  :: atnum(:)
  REAL(DP), ALLOCATABLE :: xyz(:,:), buffer(:)
  REAL(DP), ALLOCATABLE :: force_d3(:,:), hess_d3(:,:,:,:,:,:,:)
  !
9035 FORMAT(5X,'atom ',I4,' type ',I2,'   force = ',3F14.8)
  !
  ! DFT-D3 functional dependent parameters have been set in read_file_new
  !
  WRITE( stdout, '(/,5x,A,f24.12)') 'Differentiation step: ',  step 
  WRITE( stdout, '(5x,A,3I4)') 'DFT-D3 version: ',  dftd3%version  
  WRITE( stdout, '(5x,A,L)') 'DFT-D3 threebody: ',  .not.dftd3%noabc
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
  ! xyz are atomic positions centered around r=0 (in bohr units)
  xyz(:,:) = tau(:,:)
  CALL cryst_to_cart( nat, xyz, bg, -1 )
  xyz(:,:) = xyz(:,:) - NINT(xyz(:,:))
  CALL cryst_to_cart( nat, xyz, at,  1 )
  xyz(:,:) = xyz(:,:)*alat
  DO iat = 1, nat
     atnum(iat) = get_atomic_number(TRIM(atm(ityp(iat))))
  ENDDO
  !
  call dftd3_pbc_gdisp_new(dftd3, xyz, atnum, latvecs, force_d3, rep_cn, rep_vdw)
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
  if(debug) Call d2ionq_dispd3_debug( alat, nat, ityp, at, tau )
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
  ALLOCATE( hess_d3(-rep_hes(3):rep_hes(3),-rep_hes(2):rep_hes(2),-rep_hes(1):rep_hes(1),3,nat,3,nat) )
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
  IF( ionode ) THEN 
    !
    WRITE( stdout, '(/,5x,2A)') 'Writing Hessian on file: ', TRIM(filhess)
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
              buffer(j) = hess_d3(krep,jrep,irep,ixyz,iat,jxyz,jat)
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
  END IF 
  !
  CALL print_clock('force_dftd3')
  CALL print_clock('hessian_dftd3')
  !
END SUBROUTINE d3hess_sub
!
!---------------------------------------------------------------------------
SUBROUTINE d2ionq_dispd3_debug( alat, nat, ityp, at, tau )
  !------------------------------------------------------------------------
  !! This routine calculates the Grimme-D3 contribution to the dynamical matrix.
  !
  USE io_global,        ONLY: stdout
  USE symme,            ONLY: symvector
  USE mp,               ONLY: mp_stop
  USE funct,            ONLY: get_dft_short
  USE dftd3_api,        ONLY: dftd3_init, dftd3_set_functional, dftd3_pbc_dispersion, get_atomic_number
  USE dftd3_qe,         ONLY: dftd3, dftd3_pbc_gdisp_new, print_dftd3_hessian
  USE ions_base,        ONLY: atm
  USE ener,             ONLY: edftd3
  
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: nat
  !! number of atoms in the unit cell
  REAL(DP), INTENT(IN) :: alat
  !! cell parameter (celldm(1))
  INTEGER, INTENT(IN) :: ityp(nat)
  !! atomic types for atoms in the unit cell
  REAL(DP), INTENT(IN) :: at(3,3)
  !! at(:,i) is lattice vector i in alat units
  REAL(DP), INTENT(IN) :: tau(3,nat)
  !! atomic positions in alat units
  !
  !  Local variables  
  CHARACTER(LEN=256):: dft_ , formt
  REAL(DP) :: latvecs(3,3)
  REAL(DP) :: step, eerr, eerl, eelr, eell
  INTEGER:: atnum(1:nat), i, j, iat, ixyz, jat, jxyz
  REAL(DP) :: xyz(3,nat)
  COMPLEX(DP), ALLOCATABLE :: mat(:,:,:,:)
  REAL(DP), ALLOCATABLE :: force_d3(:,:), force_num(:,:), buffer(:)
  REAL(DP), ALLOCATABLE :: der2disp_ene(:,:,:,:), der2disp_frc(:,:,:,:) 
  
9078 FORMAT( '     DFT-D3 Dispersion         =',F17.8,' Ry' )
9035 FORMAT(5X,'atom ',I4,' type ',I2,'   force = ',3F14.8)

  write(stdout,'(A)') 'Adding Grimme-D3 contribution to the dynamical matrix'
  WRITE ( stdout , 9078 ) edftd3 
  !
  ! Setting DFT-D3 functional dependent parameters
  !

  CALL start_clock('force_dftd3')
  ALLOCATE( force_d3(3, nat), force_num(3, nat) )
  force_num = 0._dp
  force_d3(:,:) = 0.0_DP
  latvecs(:,:)=at(:,:)*alat
  xyz(:,:)=tau(:,:)*alat
  DO iat = 1, nat
     atnum(iat) = get_atomic_number(TRIM(atm(ityp(iat))))
  ENDDO
  call dftd3_pbc_gdisp_new(dftd3, xyz, atnum, latvecs, force_d3)
  edftd3=edftd3*2.d0
  force_d3 = -2.d0*force_d3
  
  step=1.d-6
  !step=2.d-5

  do iat = 1, nat
    do ixyz = 1, 3 
      xyz(ixyz,iat)=xyz(ixyz,iat)+step
      call dftd3_pbc_dispersion(dftd3, xyz, atnum, latvecs, eerr)
      eerr = eerr*2.0d0
      xyz(ixyz,iat)=xyz(ixyz,iat)-2*step
      call dftd3_pbc_dispersion(dftd3, xyz, atnum, latvecs, eell)
      eell = eell*2.0d0
      force_num(ixyz,iat)=-0.5*(eerr-eell)/step
      xyz(ixyz,iat)=xyz(ixyz,iat)+step
    end do 
  end do 
  !CALL symvector( nat, force_num )
  CALL stop_clock('force_dftd3')

  WRITE ( stdout , 9078 ) edftd3  
  WRITE( stdout, '(/,5x,"DFT-D3 dispersion contribution to forces (analytical):")')
  DO iat = 1, nat
     WRITE( stdout, 9035) iat, ityp(iat), (force_d3(ixyz,iat), ixyz = 1, 3)
  ENDDO

  WRITE( stdout, '(/,5x,"DFT-D3 dispersion contribution to forces (numerical):")')
  DO iat = 1, nat
     WRITE( stdout, 9035) iat, ityp(iat), (force_num(ixyz,iat), ixyz = 1, 3)
  ENDDO

  WRITE( stdout, '(/,5x,"DFT-D3 analytical vs numerical err:")')
  DO iat = 1, nat
     WRITE( stdout, 9035) iat, ityp(iat), (force_d3(ixyz,iat)-force_num(ixyz,iat), ixyz = 1, 3)
  ENDDO

  step=1.d-3
  CALL start_clock('dftd3')
  ALLOCATE( der2disp_ene(3,nat,3,nat), der2disp_frc(3,nat,3,nat) )
  der2disp_ene = 0._dp
  der2disp_frc = 0._dp

  CALL start_clock('dftd3:frc')

  do iat = 1, nat
    do ixyz = 1, 3 
      write(*,*) 'displacing forces: ', iat, ixyz

      xyz(ixyz,iat)=xyz(ixyz,iat)+step
      call dftd3_pbc_gdisp_new(dftd3, xyz, atnum, latvecs, force_d3)
      force_d3 = -2.d0*force_d3
      xyz(ixyz,iat)=xyz(ixyz,iat)-2*step
      call dftd3_pbc_gdisp_new(dftd3, xyz, atnum, latvecs, force_num)
      force_num = -2.d0*force_num

      der2disp_frc(ixyz,iat,1:3,1:nat) = -0.5 * (force_d3(1:3,1:nat) - force_num(1:3,1:nat) ) / step
      xyz(ixyz,iat)=xyz(ixyz,iat)+step
    end do 
  end do 

  CALL stop_clock('dftd3:frc')

  allocate( mat(3,nat,3,nat) )

  mat(:,:,:,:) = cmplx( der2disp_frc(:,:,:,:), kind=dp )

  CALL print_dftd3_hessian( mat, nat, 'debug' )

  deallocate( mat )
 
  CALL start_clock('dftd3:ene')

  do iat = 1, nat
    do ixyz = 1, 3 
      do jat = 1, nat
        do jxyz = 1, 3 

          write(*,*) 'displacing energy: ', iat, ixyz, jat, jxyz

          xyz(ixyz,iat)=xyz(ixyz,iat)+step
          xyz(jxyz,jat)=xyz(jxyz,jat)+step
          call dftd3_pbc_dispersion(dftd3, xyz, atnum, latvecs, eerr)  ! i+s  j+s
          eerr = eerr*2.0d0

          xyz(ixyz,iat)=xyz(ixyz,iat)-2*step
          call dftd3_pbc_dispersion(dftd3, xyz, atnum, latvecs, eelr)  ! i-s  j+s
          eelr = eelr*2.0d0

          xyz(jxyz,jat)=xyz(jxyz,jat)-2*step
          call dftd3_pbc_dispersion(dftd3, xyz, atnum, latvecs, eell)  ! i-s  j-s
          eell = eell*2.0d0

          xyz(ixyz,iat)=xyz(ixyz,iat)+2*step
          call dftd3_pbc_dispersion(dftd3, xyz, atnum, latvecs, eerl)  ! i+s  j-s
          eerl = eerl*2.0d0

          xyz(ixyz,iat)=xyz(ixyz,iat)-step
          xyz(jxyz,jat)=xyz(jxyz,jat)+step

          der2disp_ene(ixyz,iat,jxyz,jat) = (eerr - eerl - eelr + eell) / 4.0d0 / step / step

        end do 
      end do 
    end do 
  end do 

  CALL stop_clock('dftd3:ene')

  !der2disp = der2disp_ene
  !der2disp = der2disp_frc

  CALL stop_clock('dftd3')

  write(formt,'(A,I5,A)') '(',3*nat,'F12.8)'
  WRITE( stdout, '(/,5x,"DFT-D3 numerical hessian (from Forces):")')
  DO iat = 1, nat
     DO ixyz = 1, 3
       WRITE( stdout, formt) der2disp_frc(ixyz,iat,1:3,1:nat)
     END DO 
  ENDDO

  WRITE( stdout, '(/,5x,"DFT-D3 numerical hessian (from Energies):")')
  DO iat = 1, nat
     DO ixyz = 1, 3
       WRITE( stdout, formt) der2disp_ene(ixyz,iat,1:3,1:nat)
     END DO 
  ENDDO


  WRITE( stdout, '(/,5x,"DFT-D3 numerical hessian (from Forces vs from Energies):")')
  DO iat = 1, nat
     DO ixyz = 1, 3
       WRITE( stdout, formt) der2disp_frc(ixyz,iat,1:3,1:nat)-der2disp_ene(ixyz,iat,1:3,1:nat)
     END DO 
  ENDDO


  DEALLOCATE( force_d3, force_num, der2disp_ene, der2disp_frc)

! Call mp_stop(555)

  RETURN

END SUBROUTINE d2ionq_dispd3_debug
!
END MODULE d3hess_mod
