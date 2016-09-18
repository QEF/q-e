! Copyright (C) 2006-2008 Dmitry Korotin dmitry@korotin.name
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#define ZERO (0.d0,0.d0)
#define ONE (1.d0,0.d0)

!-----------------------------------------------------------------------
SUBROUTINE write_hamiltonian_default(nwan,hamk,iunhamilt)
!-----------------------------------------------------------------------

  USE io_global, ONLY: stdout
  USE kinds, ONLY: DP
  USE constants,  ONLY : rytoev
  USE klist, ONLY: nks, wk
  
  IMPLICIT NONE
  INTEGER, INTENT(in) :: nwan, iunhamilt
  COMPLEX(DP), INTENT(in) :: hamk(nwan,nwan,nks)

  INTEGER :: i,j, ik, ios
  COMPLEX(DP), ALLOCATABLE :: hamk2(:,:)
  real(DP) :: eps = 1.d-8, hr,hi

  WRITE(stdout,'(/5x,a36,i5,a9)') 'Hamiltonian is in the default format,', nks, 'k-points'
  WRITE(stdout,'(5x,a48/)') 'ATTENTION: All k-points weights are real weights'

  OPEN (iunhamilt, file = 'hamilt', status = 'unknown', form = 'formatted', err = 300, iostat = ios)
300 CALL errore ('HMLT', 'Opening hamilt', abs (ios) )

  WRITE(iunhamilt,*) nks,nwan

  ALLOCATE(hamk2(nwan,nwan))

  DO ik = 1, nks

     WRITE(iunhamilt,'(f15.12)') wk(ik)

     hamk2 = hamk(:,:,ik) * rytoev
     
     DO i=1, nwan
        DO j=1, nwan

           hr = abs(dreal(hamk2(i,j)))
           hi = abs(dimag(hamk2(i,j)))
           IF((hr>=eps).and.(hi>=eps)) WRITE(iunhamilt,'(2f13.8)') dreal(hamk2(i,j)), aimag(hamk2(i,j))
           IF ((hr<eps).and.(hi>=eps)) WRITE(iunhamilt,'(f3.0,f13.8)') 0., aimag(hamk2(i,j))
           IF ((hr>=eps).and.(hi<eps)) WRITE(iunhamilt,'(f13.8,f3.0)') dreal(hamk2(i,j)), 0.
           IF ((hr<eps).and.(hi<eps)) WRITE(iunhamilt,'(2f3.0)') 0., 0.

        ENDDO
     ENDDO

  ENDDO

  DEALLOCATE(hamk2)

  CLOSE(iunhamilt)

END SUBROUTINE write_hamiltonian_default

!-----------------------------------------------------------------------
SUBROUTINE write_hamiltonian_amulet(nwan,hamk,hash,iunhamilt)
!-----------------------------------------------------------------------
! Special hamiltonian format for the AMULET code instegration
! www.amulet-code.org

  USE io_global, ONLY: stdout
  USE kinds, ONLY: DP
  USE constants,  ONLY : rytoev
  USE klist, ONLY: nks, wk, xk
  USE lsda_mod, ONLY : nspin
  USE input_parameters, ONLY : title
  USE global_version, ONLY : version_number
  
  IMPLICIT NONE
  INTEGER, INTENT(in) :: nwan, hash, iunhamilt
  COMPLEX(DP), INTENT(in) :: hamk(nwan,nwan,nks)

  INTEGER :: i,j, ik, ios
  COMPLEX(DP), ALLOCATABLE :: hamk2(:,:)
  real(DP) :: eps = 1.d-8, hr,hi
  CHARACTER(LEN=9)  :: cdate, ctime

  CALL date_and_tim( cdate, ctime )

  WRITE(stdout,'(/5x,a36,i5,a9)') 'Hamiltonian is in the AMULET format,', nks/nspin, 'k-points'
  WRITE(stdout,'(5x,a48/)') 'ATTENTION: All k-points weights are real weights'

  OPEN (iunhamilt, file = 'hamilt.am', status = 'unknown', form = 'formatted', err = 300, iostat = ios)
300 CALL errore ('HMLT', 'Opening hamilt', abs (ios) )

  write(iunhamilt,'(a30,2a10/)') '# This file was generated on: ', cdate, ctime
  IF( trim(title) .NE. '' ) write(iunhamilt,'(a2,a80/)') '# ', title
  
  WRITE(iunhamilt,'(a10)') '&codestamp'
  WRITE(iunhamilt,'(a3,a6)') 'QE_', version_number
  WRITE(iunhamilt,*) 

  WRITE(iunhamilt,'(a5)') '&hash'
  WRITE(iunhamilt,*) hash
  WRITE(iunhamilt,*) 

  WRITE(iunhamilt,'(a6)') '&nspin'
  WRITE(iunhamilt,'(i1)') nspin
  WRITE(iunhamilt,*) 

  WRITE(iunhamilt,'(a4)') '&nkp'
  WRITE(iunhamilt,'(i5)') nks/nspin
  WRITE(iunhamilt,*) 

  WRITE(iunhamilt,'(a4)') '&dim'
  WRITE(iunhamilt,'(i3)') nwan
  WRITE(iunhamilt,*) 

  WRITE(iunhamilt,'(a8)') '&kpoints'
  DO ik=1, nks/nspin
    WRITE(iunhamilt,'(f15.12,3f9.5)') wk(ik), xk(:,ik)
  END DO
  WRITE(iunhamilt,*) 

  WRITE(iunhamilt,'(a12)') '&hamiltonian'
  
  ALLOCATE(hamk2(nwan,nwan))

  DO ik = 1, nks

     hamk2 = hamk(:,:,ik) * rytoev
     
     DO i=1, nwan
        DO j=i, nwan

           hr = abs(dreal(hamk2(i,j)))
           hi = abs(dimag(hamk2(i,j)))
           IF((hr>=eps).and.(hi>=eps)) WRITE(iunhamilt,'(2f13.8)') dreal(hamk2(i,j)), aimag(hamk2(i,j))
           IF ((hr<eps).and.(hi>=eps)) WRITE(iunhamilt,'(f3.0,f13.8)') 0., aimag(hamk2(i,j))
           IF ((hr>=eps).and.(hi<eps)) WRITE(iunhamilt,'(f13.8,f3.0)') dreal(hamk2(i,j)), 0.
           IF ((hr<eps).and.(hi<eps)) WRITE(iunhamilt,'(2f3.0)') 0., 0.

        ENDDO
     ENDDO

  ENDDO

  DEALLOCATE(hamk2)

  CLOSE(iunhamilt)

END SUBROUTINE write_hamiltonian_amulet

!-----------------------------------------------------------------------
SUBROUTINE write_systemdata_amulet(hash,nelec,iunsystem)
!-----------------------------------------------------------------------
! Damp of the system data for the AMULET code instegration
! www.amulet-code.org

  USE io_global, ONLY: stdout
  USE kinds, ONLY: DP
  USE constants,  ONLY : rytoev
  USE klist, ONLY: nks, wk, xk
  USE lsda_mod, ONLY : nspin
  USE input_parameters, ONLY : title
  USE ions_base, ONLY : nat, atm, tau, ityp
  USE cell_base, ONLY : alat, at
  USE ener, ONLY : ef
  USE wannier_new, ONLY : nwan, wan_in
  USE global_version, ONLY : version_number

  IMPLICIT NONE
  INTEGER, INTENT(in) :: hash,iunsystem
  REAL(DP), INTENT(in) :: nelec
  INTEGER :: ios, i, j
  CHARACTER(LEN=9)  :: cdate, ctime
  CHARACTER :: l_symb(4)
  DATA l_symb/'s','p','d','f'/
  INTEGER :: orbitals(7,4)
  DATA orbitals/ 1, 0, 0, 0, 0, 0, 0, &
                 3, 4, 2, 0, 0, 0, 0, &
                 7, 8, 6, 9, 5, 0, 0, &
                 13, 14, 12, 15, 11, 16, 10 /
  INTEGER, PARAMETER :: nwannierblocksmax = 25

  INTEGER :: nblocks, block_dim(nwannierblocksmax), block_l(nwannierblocksmax), &
             block_atom(nwannierblocksmax), block_wannier(nwannierblocksmax,9), block_start(nwannierblocksmax)

  interface 
    SUBROUTINE split_basis_into_blocks(nblocks,block_dim,block_l,block_atom,block_wannier,block_start)
      INTEGER, INTENT(OUT) :: nblocks, block_dim(:), block_atom(:), block_l(:), block_wannier(:,:), block_start(:)
    END SUBROUTINE split_basis_into_blocks
  end interface

  CALL date_and_tim( cdate, ctime )

  OPEN (iunsystem, file = 'system.am', status = 'unknown', form = 'formatted', err = 300, iostat = ios)
300 CALL errore ('HMLT', 'Opening system.am', abs (ios) )

  write(iunsystem,'(a30,2a10/)') '# This file was generated on: ', cdate, ctime
  IF( trim(title) .NE. '' ) write(iunsystem,'(a2,a80/)') '# ', title
  
  WRITE(iunsystem,'(a5)') '&hash'
  WRITE(iunsystem,*) hash
  WRITE(iunsystem,*) 

  WRITE(iunsystem,'(a10)') '&codestamp'
  WRITE(iunsystem,'(a3,a6)') 'QE_',version_number
  WRITE(iunsystem,*) 

  WRITE(iunsystem,'(a5)') '&cell'
  WRITE(iunsystem,'(f12.9)') alat
  DO i=1,3
    WRITE(iunsystem,'(3f9.5)') at(:,i)
  END DO
  WRITE(iunsystem,*) 

  WRITE(iunsystem,'(a6)') '&atoms'
  WRITE(iunsystem,'(i5)') nat
  DO i=1, nat
    WRITE(iunsystem,'(a4,1x,3f9.5)') atm(ityp(i)), tau(:,i)
  END DO
  WRITE(iunsystem,*)

  WRITE(iunsystem,'(a6)') '&nelec'
  WRITE(iunsystem,'(f7.2)') nelec
  WRITE(iunsystem,*)

  WRITE(iunsystem,'(a7)') '&efermi'
  WRITE(iunsystem,'(f8.4)') ef*rytoev
  WRITE(iunsystem,*)

  CALL split_basis_into_blocks(nblocks,block_dim,block_l,block_atom,block_wannier,block_start)
  WRITE(iunsystem,'(a20)') '# Basis description:'
  WRITE(iunsystem,'(a14)') '# dim, nblocks'
  WRITE(iunsystem,'(a74)') '# atom_sym, atom_num, l_sym, block_dim, block_start, orbitals(1:block_dim)'
  WRITE(iunsystem,'(a6)') '&basis'
  WRITE(iunsystem,'(i2,i4)') nwan, nblocks
  DO i=1, nblocks
    WRITE(iunsystem,'(a3,i3,a2,i2,i4,4x,7i2)') atm(ityp(block_atom(i))), block_atom(i), l_symb(block_l(i)+1), &
                       block_dim(i), block_start(i), ( orbitals(wan_in(j,1)%ing(1)%m,block_l(i)+1), &
                        j=block_wannier(i,1), block_wannier(i,block_dim(i)) )
  END DO
  WRITE(iunsystem,*)

  CLOSE(iunsystem)

END SUBROUTINE write_systemdata_amulet

SUBROUTINE split_basis_into_blocks(nblocks,block_dim,block_l,block_atom,block_wannier,block_start)
  USE kinds, only : DP
  USE wannier_new, only : nwan, wan_in
  USE io_global, ONLY: stdout

  implicit none
  INTEGER :: i,j, iblock
  INTEGER, INTENT(OUT) :: nblocks, block_dim(:), block_atom(:), block_l(:), block_wannier(:,:), block_start(:)
  
  nblocks = 0
  block_dim = 0
  block_l = -1
  block_atom = 0
  block_wannier = 0

  iblock = 1
  block_start(iblock) = 1
  
  j=1
  DO i = 1, nwan-1
    block_wannier(iblock,j) = i
    j = j+1
    IF( ((wan_in(i,1)%iatom .ne. wan_in(i+1,1)%iatom) .OR. wan_in(i,1)%ing(1)%l .ne. wan_in(i+1,1)%ing(1)%l) &
        ) THEN
      block_dim(iblock) = i - block_start(iblock) + 1
      block_atom(iblock) = wan_in(i,1)%iatom
      block_l(iblock) = wan_in(i,1)%ing(1)%l
      iblock = iblock + 1
      block_start(iblock) = i+1
      j=1
    END IF
  END DO

  ! the last wannier
  block_dim(iblock) = nwan - block_start(iblock) + 1
  block_atom(iblock) = wan_in(nwan,1)%iatom
  block_l(iblock) = wan_in(nwan,1)%ing(1)%l
  block_wannier(iblock,j) = i

  nblocks = iblock
  
!   write(stdout,'(/5x,a19)') 'Blocks of orbitals:'
!   DO iblock = 1, nblocks
!     write(stdout,*) 'Block', iblock, block_dim(iblock), block_atom(iblock), & 
!                             block_l(iblock), block_wannier(iblock,1:block_dim(iblock))
!   END DO

END SUBROUTINE split_basis_into_blocks
