!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "f_defs.h"
!
!--------------------------------------------------------------------
PROGRAM voronoy
  !--------------------------------------------------------------------
  !
  !  Calculates charges on atoms by dividing the space into Voronoy
  !  polyhedra. A Voronoy polyhedron around a given atom is defined
  !  as the region of space that is closer to that atom than to the others
  !
  !  Note that this is a very rough way to associate charges to atoms
  !  and that it is well defined only if all atoms are of the same type!
  !
  !  On input: nr1big, nr2bug, nr3big are FFT grid dimensions larger than
  !  the original ones (the largest, the more accurate the integration)
  !
  USE constants,  ONLY : pi
  USE io_global,  ONLY : stdout
  USE kinds
  USE char,       ONLY : title
  USE cell_base,  ONLY : ibrav, celldm, at, bg, alat, omega, tpiba, tpiba2
  USE ions_base,  ONLY : atm, zv, nat, tau, ityp, ntyp=>nsp
  USE lsda_mod,   ONLY : nspin
  USE gvect,      ONLY : ecutwfc, nr1, nr2, nr3, nrxx, gcutm, ngm, g, nl, &
                         dual, nrx1, nrx2, nrx3
  USE scf,        ONLY : rho
  USE fft_scalar, ONLY : good_fft_dimension
  USE io_files,   ONLY : nd_nmbr
  USE mp_global,  ONLY : mpime, root

  IMPLICIT NONE
  INTEGER :: nr1big, nr2big, nr3big, nrx1big
  INTEGER :: n, i, j, ng, na, plot_num
  REAL(DP) :: total_charge, rhodum
  INTEGER, ALLOCATABLE :: nlbig (:)
  REAL(DP), ALLOCATABLE :: partial_charge (:)
  COMPLEX(DP), ALLOCATABLE :: rhobig (:)
  CHARACTER(len=256) :: filename
  !
  CALL start_postproc (nd_nmbr)
  !
  ! Works for parallel machines but only for one processor !!!
  !
  IF ( mpime == root ) THEN
     !
     PRINT '(" Input file > ",$)'
     READ (5, '(a)') filename
     !
     ! read file header and allocate objects
     !
     CALL read_io_header (filename, title, nrx1, nrx2, nrx3, nr1, nr2, nr3, nat, &
          ntyp, ibrav, celldm, at, gcutm, dual, ecutwfc, plot_num)
     !
     ALLOCATE(tau(3, nat))
     ALLOCATE(ityp(nat))

     CALL latgen (ibrav, celldm, at (1, 1), at (1, 2), at (1, 3), omega )
     alat = celldm (1) ! define alat
     at = at / alat    ! bring at in units of alat

     tpiba = 2.d0 * pi / alat
     tpiba2 = tpiba**2
     CALL recips (at (1, 1), at (1, 2), at (1, 3), bg (1, 1), bg (1, 2) &
          , bg (1, 3) )
     CALL volume (alat, at (1, 1), at (1, 2), at (1, 3), omega)
     CALL set_fft_dim
     nrxx = nrx1 * nrx2 * nrx3
     nspin = 1
     CALL allocate_fft
     !
     ! read data from file
     !
     CALL plot_io (filename, title, nrx1, nrx2, nrx3, nr1, nr2, nr3, &
          nat, ntyp, ibrav, celldm, at, gcutm, dual, ecutwfc, plot_num, &
          atm, ityp, zv, tau, rho, - 1)
     !
     ! calculate g-vectors
     !
     CALL ggen
     !
     ! interpolate on a larger grid
     !
     PRINT '(" nr1big, nr2big, nr3big (min: ",3i4," ) > ",$)', nr1, &
          nr2, nr3

     READ (5, * ) nr1big, nr2big, nr3big

     nrx1big = good_fft_dimension (nr1big)

     ALLOCATE (nlbig( ngm))    

     CALL get_fftindex (g, ngm, nr1big, nr2big, nr3big, nrx1big, &
          nr2big, nr3big, at, nlbig)

     ALLOCATE (rhobig(nrx1big * nr2big * nr3big))    

     CALL rhor_to_rhobig (ngm, nr1, nr2, nr3, nrx1, nl, rho, nr1big, &
          nr2big, nr3big, nrx1big, nlbig, rhobig)

     ALLOCATE (partial_charge(nat))

     CALL calculate_partial_charges (nat, tau, at, bg, nrx1big, nr1big, &
          nr2big, nr3big, rhobig, partial_charge)

     WRITE( stdout, '(" nr1, nr2, nr3 = ",3i4)') nr1big, nr2big, nr3big
     total_charge = 0.0
     DO na = 1, nat
        partial_charge (na) = partial_charge (na) * omega / (nr1big * &
             nr2big * nr3big)
        WRITE( stdout, '(" Atom # ",i2," tau = ",3f8.4, &
             &       "  charge: ",f8.4) ') na,  (tau (i, na) , i = 1, 3) , &
             partial_charge (na)
        total_charge = total_charge+partial_charge (na)
     ENDDO

     WRITE( stdout, '(" Check: total charge = ",f8.4)') total_charge

  END IF

  CALL stop_pp
  !
END PROGRAM voronoy
!
!-----------------------------------------------------------------------
SUBROUTINE get_fftindex (g, ngm, nr1, nr2, nr3, nrx1, nrx2, nrx3, at, nl)
  !-----------------------------------------------------------------------
  USE kinds, ONLY: DP
  IMPLICIT NONE
  INTEGER :: nr1, nr2, nr3, nrx1, nrx2, nrx3, ngm
  REAL(DP) :: g (3, ngm), at (3, 3)
  INTEGER :: nl (ngm)
  INTEGER :: n1, n2, n3, ng
  !
  DO ng = 1, ngm
     n1 = NINT (g (1, ng) * at (1, 1) + g (2, ng) * at (2, 1) + g (3, &
          ng) * at (3, 1) ) + 1
     IF (n1.LT.1) n1 = n1 + nr1
     n2 = NINT (g (1, ng) * at (1, 2) + g (2, ng) * at (2, 2) + g (3, &
          ng) * at (3, 2) ) + 1
     IF (n2.LT.1) n2 = n2 + nr2
     n3 = NINT (g (1, ng) * at (1, 3) + g (2, ng) * at (2, 3) + g (3, &
          ng) * at (3, 3) ) + 1
     IF (n3.LT.1) n3 = n3 + nr3
     IF (n1.LE.nr1.AND.n2.LE.nr2.AND.n3.LE.nr3) THEN
        nl (ng) = n1 + (n2 - 1) * nrx1 + (n3 - 1) * nrx1 * nrx2
     ELSE
        STOP ' Mesh too small '
     ENDIF
  ENDDO
  !
  RETURN

END SUBROUTINE get_fftindex
!-----------------------------------------------------------------------
SUBROUTINE rhor_to_rhobig (ngm, nr1, nr2, nr3, nrx1, nl, rho, &
     nr1big, nr2big, nr3big, nrx1big, nlbig, rhobig)
  !-----------------------------------------------------------------------
  USE kinds, ONLY: DP
  IMPLICIT NONE
  INTEGER :: ngm, nr1, nr2, nr3, nrx1, nl (ngm), nr1big, nr2big, &
       nr3big, nrx1big, nlbig (ngm)
  REAL(DP) :: rho (nrx1 * nr2 * nr3)
  COMPLEX(DP) :: rhobig (nrx1big * nr2big * nr3big)
  COMPLEX(DP), ALLOCATABLE :: aux (:)
  INTEGER :: nr, ng

  ALLOCATE (aux(nrx1 * nr2 * nr3))    
  DO nr = 1, nrx1 * nr2 * nr3
     aux (nr) = rho (nr)
  ENDDO

  CALL cft3 (aux, nr1, nr2, nr3, nrx1, nr2, nr3, - 1)

  DO nr = 1, nrx1big * nr2big * nr3big
     rhobig (nr) = 0.0
  ENDDO
  DO ng = 1, ngm
     rhobig (nlbig (ng) ) = aux (nl (ng) )
  ENDDO

  DEALLOCATE(aux)

  CALL cft3 (rhobig, nr1big, nr2big, nr3big, nrx1big, nr2big, &
       nr3big, 1)
  RETURN

END SUBROUTINE rhor_to_rhobig
!-----------------------------------------------------------------------
SUBROUTINE calculate_partial_charges (nat, tau, at, bg, nrx1big, &
     nr1big, nr2big, nr3big, rhobig, partial_charge)
  !-----------------------------------------------------------------------
  USE kinds, ONLY: DP
  IMPLICIT NONE
  INTEGER :: nat, nrx1big, nr1big, nr2big, nr3big
  REAL(DP) :: at (3, 3), bg (3, 3), tau (3, nat), partial_charge(nat)

  COMPLEX(DP) :: rhobig (nrx1big, nr2big, nr3big)
  INTEGER :: atom (nat), equidistant, n1, n2, n3, na, i
  REAL(DP) :: r (3), dr (3), distance (nat), relative_distance
  !
  !
  DO na = 1, nat
     partial_charge (na) = 0.0
  ENDDO
  DO n3 = 1, nr3big
     DO n2 = 1, nr2big
        DO n1 = 1, nr1big
           ! r is the position of point (n1,n2,n3) on the grid
           DO i = 1, 3
              r (i) = (n1 - 1) * at (i, 1) / nr1big + (n2 - 1) * at (i, 2) &
                   / nr2big + (n3 - 1) * at (i, 3) / nr3big
           ENDDO
           DO na = 1, nat
              ! dr is the distance between r and this atom, in crystal axis
              DO i = 1, 3
                 dr (i) = (r (1) - tau (1, na) ) * bg (1, i) + &
                          (r (2) - tau (2, na) ) * bg (2, i) + &
                          (r (3) - tau (3, na) ) * bg (3, i)
                 ! this brings dr back into the unit cell
                 dr (i) = dr (i) - NINT (dr (i) )
              ENDDO
              ! distance is in cartesioan axis
              distance (na) = (dr (1) * at (1, 1) + dr (2) * at (1, 2) + dr (3) &
                   * at (1, 3) ) **2 + (dr (1) * at (2, 1) + dr (2) * at (2, 2) &
                   + dr (3) * at (2, 3) ) **2 + (dr (1) * at (3, 1) + dr (2) * at (3, &
                   2) + dr (3) * at (3, 3) ) **2
           ENDDO
           ! initialization needed by hpsort
           atom (1) = 0
           ! sort distances in increasing order
           CALL hpsort (nat, distance, atom)
           ! find if some atoms are equidistant
           equidistant = 1
           DO na = 2, nat
              relative_distance = distance (na) - distance (1)
              IF (relative_distance.GT.1.0e-4) GOTO 10
              equidistant = equidistant + 1
           ENDDO
10         CONTINUE
           ! the charge is assigned to the closest atom or shared among equidistant
           DO na = 1, equidistant
              partial_charge (atom (na) ) = partial_charge (atom (na) ) + &
                 DBLE ( rhobig (n1, n2, n3) ) / equidistant
           ENDDO
        ENDDO
     ENDDO
  ENDDO
  RETURN
END SUBROUTINE calculate_partial_charges
