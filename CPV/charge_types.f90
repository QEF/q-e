!
! Copyright (C) 2002 FPMD group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

!  ----------------------------------------------
!  BEGIN manual

!  ----------------------------------------------  !
      MODULE charge_types
!  ----------------------------------------------  !

!  ----------------------------------------------
!  END manual


        USE kinds
        USE parameters, ONLY: nspinx
        IMPLICIT NONE
        PRIVATE
        SAVE

!  BEGIN manual
!  TYPE DEFINITIONS
     
        TYPE charge_descriptor

          INTEGER :: ldx  ! leading dimension for x direction
          INTEGER :: ldy  ! leading dimension for y
          INTEGER :: ldz  ! leading dimension for z
          INTEGER :: lds  ! leading dimension for spin array dimension
       
          INTEGER :: nxl ! local number of cell in x direction
          INTEGER :: nx  ! global number of cell in x direction
          INTEGER :: nyl ! local number of cell in y direction
          INTEGER :: ny  ! global number of cell in y direction
          INTEGER :: nzl ! local number of cell in z direction
          INTEGER :: nz  ! global number of cell in z direction

          INTEGER :: nspin ! number of spin

        END TYPE charge_descriptor

!  ----------------------------------------------
!  END manual

        PUBLIC :: charge_descriptor, charge_descriptor_init, charge_descriptor_info

!  end of module-scope declarations
!  ----------------------------------------------

      CONTAINS

!  ----------------------------------------------  !
!  subroutines


      SUBROUTINE charge_descriptor_init( desc, nx, ny, nz, nxl, nyl, nzl, ldx, ldy, ldz, nspin )

        IMPLICIT NONE

        TYPE (charge_descriptor), INTENT(OUT) :: desc
        INTEGER, INTENT(IN) :: nx
        INTEGER, INTENT(IN) :: ny
        INTEGER, INTENT(IN) :: nz
        INTEGER, INTENT(IN) :: nxl
        INTEGER, INTENT(IN) :: nyl
        INTEGER, INTENT(IN) :: nzl
        INTEGER, INTENT(IN) :: ldx
        INTEGER, INTENT(IN) :: ldy
        INTEGER, INTENT(IN) :: ldz
        INTEGER, INTENT(IN) :: nspin

        INTEGER :: is

        !  g vectors

        IF( nx <= 0 ) &
          CALL errore( ' charge_descriptor_init ', ' arg no. 2 out of range ', 1 ) 
        IF( ny <= 0 ) &
          CALL errore( ' charge_descriptor_init ', ' arg no. 3 out of range ', 1 ) 
        IF( nz <= 0 ) &
          CALL errore( ' charge_descriptor_init ', ' arg no. 4 out of range ', 1 ) 
        IF( nxl > nx ) &
          CALL errore( ' charge_descriptor_init ', ' arg no. 5 out of range ', 1 ) 
        IF( nyl > ny ) &
          CALL errore( ' charge_descriptor_init ', ' arg no. 6 out of range ', 1 ) 
        IF( nzl > nz ) &
          CALL errore( ' charge_descriptor_init ', ' arg no. 7 out of range ', 1 ) 
        IF( ldx < nxl ) &
          CALL errore( ' charge_descriptor_init ', ' arg no. 8 out of range ', 1 ) 
        IF( ldy < nyl ) &
          CALL errore( ' charge_descriptor_init ', ' arg no. 9 out of range ', 1 ) 
        IF( ldz < nzl ) &
          CALL errore( ' charge_descriptor_init ', ' arg no. 10 out of range ', 1 ) 
        IF( nspin < 1 .OR. nspin > 2 ) &
          CALL errore( ' charge_descriptor_init ', ' arg no. 11 out of range ', 1 ) 

        desc % nx = nx
        desc % ny = ny
        desc % nz = nz
        desc % nxl = nxl
        desc % nyl = nyl
        desc % nzl = nzl
        desc % nspin = nspin
        desc % ldx = ldx
        desc % ldy = ldy
        desc % ldz = ldz
        desc % lds = nspin

        RETURN
      END SUBROUTINE


      SUBROUTINE charge_descriptor_info( desc, nam, iunit )

        IMPLICIT NONE

        TYPE (charge_descriptor), INTENT(IN) :: desc
        INTEGER, INTENT(IN) :: iunit
        CHARACTER(LEN=*) :: nam

        WRITE( iunit, 10 ) nam, desc % nx, desc % ny, desc % nz, desc % nspin, &
          desc % nxl, desc % nyl, desc % nzl, desc % ldx, desc % ldy,          &
          desc % ldz, desc % lds


10      FORMAT( 3X, 'Charge density descriptor  . . . . : ',A20,/, &
                3X, 'global dimensions (x,y,z,s)  . . . : ',4I5,/, &
                3X, 'local dimensions (x,y,z) . . . . . : ',3I5,/, &
                3X, 'leading dimensions (lx,ly,lz,ls) . : ',4I5)

        RETURN
      END SUBROUTINE

!  ----------------------------------------------  !
      END MODULE
!  ----------------------------------------------  !
