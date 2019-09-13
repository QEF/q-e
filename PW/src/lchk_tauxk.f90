!
! Copyright (C) 2001-2008 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
SUBROUTINE check_atoms( nvec, vec, trmat )
  !-----------------------------------------------------------------------
  !! This routine tests that the atomic coordinates (or k-points)
  !! are different and not related by a lattice translation.
  !
  USE kinds
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: nvec
  !! number of atomic positions (or k-points)
  REAL(DP), INTENT(IN) :: vec(3,nvec)
  !! cartesian coordinates of atomic positions (or k-points)
  REAL(DP), INTENT(IN) :: trmat(3,3)
  !! transformation matrix to crystal axis
  ! ( = bg , basis of the real-space lattice, for atoms
  !   = at , basis of the rec.-space lattice, for k-points )
  !
  ! ... local variables
  !
  INTEGER :: nv1, nv2
  REAL(DP), ALLOCATABLE :: vaux(:,:)
  REAL(DP) :: zero(3) = 0.0_DP
  CHARACTER(LEN=80) :: message
  REAL(DP), PARAMETER :: accep=1.d-5
  LOGICAL, EXTERNAL :: eqvect
  !
  ! ... Copy input positions and transform them to crystal units
  !
  ALLOCATE( vaux(3,nvec) )
  vaux = vec
  !
  CALL cryst_to_cart( nvec, vaux, trmat, -1 )
  !
  ! ... Test that all the atomic positions (or k-points) are different
  !
  DO nv1 = 1, nvec-1
     DO nv2 = nv1+1, nvec
        !
        IF ( eqvect( vaux(1,nv1), vaux(1,nv2), zero, accep ) ) THEN
           zero(:) = vaux(:,nv1) - vaux(:,nv2)
           IF ( ALL ( ABS(zero) < accep ) ) THEN
              WRITE (message,'("atoms #",i4," and #",i4," overlap!")') nv1, nv2
           ELSE
              WRITE (message,'("atoms #",i4," and #",i4," differ by lattice &
                    &vector (",i2,",",i2,",",i2,") in crystal axis")') &
                    nv1, nv2, NINT(zero)
           ENDIF
           CALL errore( 'check_atoms', TRIM(message), 1 )
        ENDIF
        !
     ENDDO
  ENDDO
  !
  DEALLOCATE( vaux )
  !
  RETURN
  !
END SUBROUTINE check_atoms
