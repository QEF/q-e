!
! Copyright (C) 20102011 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!=----------------------------------------------------------------------------=!
   MODULE smallbox_gvec
!=----------------------------------------------------------------------------=!
     !! Variables describing G-vectors for the small box grid.
     !! Basically the same meaning as for the corresponding 
     !! quantities for the true lattice.
     
     USE kinds, ONLY: DP

     IMPLICIT NONE
     SAVE
     !
     INTEGER :: ngb  = 0
     !! local number of G vectors
     INTEGER :: ngbl = 0
     !! number of G-vector shells up to ngw

     REAL(DP), ALLOCATABLE :: gb(:)
     !! \(G(i)^2\) in \((\text{tpi}/\text{alatb})^2\) units
     REAL(DP), ALLOCATABLE :: gxb(:,:)
     !! \(G(:,i)\) in  \(\text{tpi}/\text{alatb}\) units
     REAL(DP), ALLOCATABLE :: glb(:)
     !! shells of \(G(i)^2\)
     INTEGER, ALLOCATABLE :: npb(:)
     !! FFT indices
     INTEGER, ALLOCATABLE :: nmb(:)
     !! FFT indices
     INTEGER, ALLOCATABLE :: mill_b(:,:)
     !! Miller indices

     REAL(DP) :: gcutb = 0.0_DP  ! effective cut-off in (tpi/alatb)**2 units

   CONTAINS

     SUBROUTINE deallocate_smallbox_gvec()
       !! Deallocate small box G-vectors related arrays.
       IF( ALLOCATED( gb ) ) DEALLOCATE( gb )
       IF( ALLOCATED( gxb ) ) DEALLOCATE( gxb )
       IF( ALLOCATED( glb ) ) DEALLOCATE( glb )
       IF( ALLOCATED( npb ) ) DEALLOCATE( npb )
       IF( ALLOCATED( nmb ) ) DEALLOCATE( nmb )
       IF( ALLOCATED( mill_b ) ) DEALLOCATE( mill_b )
     END SUBROUTINE deallocate_smallbox_gvec

!=----------------------------------------------------------------------------=!
   END MODULE smallbox_gvec
!=----------------------------------------------------------------------------=!
