!
! Copyright (C) 2002-2005 FPMD-CPV groups
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!  ----------------------------------------------
!  BEGIN manual

!  ----------------------------------------------  !
      MODULE wave_types
!  ----------------------------------------------  !

!  ----------------------------------------------
!  END manual


        USE kinds
        IMPLICIT NONE
        PRIVATE
        SAVE

!  BEGIN manual
!  TYPE DEFINITIONS
     
        TYPE wave_descriptor

          INTEGER :: ldg  ! leading dimension for pw array dimension
          INTEGER :: ldb  ! leading dimension for band array dimension
          INTEGER :: lds  ! leading dimension for spin array dimension
          INTEGER :: ldk  ! leading dimension for k-points array dimension
       
          INTEGER :: ngwl  ! local number of pw
          INTEGER :: ngwt  ! global number of pw
          INTEGER :: nbl( 2 )  ! local number of bands
          INTEGER :: nbt( 2 )  ! global number of bands
          INTEGER :: nkl  ! local number of k-points
          INTEGER :: nkt  ! global number of k-points

          INTEGER :: nspin ! number of spin

          INTEGER :: isym  ! symmetry of the wave function 
                           ! ( gamma symmetry: isym == 0 )

          LOGICAL :: gamma ! true if wave functions have gamma symmetry
   
          LOGICAL :: gzero ! true if the first plane wave is the one 
                           ! with |G| == 0

        END TYPE wave_descriptor

!  ----------------------------------------------
!  END manual

        PUBLIC :: wave_descriptor, wave_descriptor_init, wave_descriptor_info

!  end of module-scope declarations
!  ----------------------------------------------

      CONTAINS

!  ----------------------------------------------  !
!  subroutines


      SUBROUTINE wave_descriptor_init( desc, ngwl, ngwt, nbl, nbt, nkl, nkt, &
        nspin, scheme, lgz )

        IMPLICIT NONE

        TYPE (wave_descriptor), INTENT(OUT) :: desc
        INTEGER, INTENT(IN) :: ngwl
        INTEGER, INTENT(IN) :: ngwt
        INTEGER, INTENT(IN) :: nbl( : )
        INTEGER, INTENT(IN) :: nbt( : )
        INTEGER, INTENT(IN) :: nkl
        INTEGER, INTENT(IN) :: nkt
        INTEGER, INTENT(IN) :: nspin
        CHARACTER(LEN=*), INTENT(IN) :: scheme
        LOGICAL, INTENT(IN) :: lgz

        INTEGER :: is

        !  g vectors

        IF( ngwt <  0 ) &
          CALL errore( ' wave_descriptor_init ', ' arg no. 3 out of range ', 1 ) 

        desc % ngwt = ngwt

        IF( ngwl <= 0 ) THEN
          desc % ngwl = ngwt
        ELSE IF( ngwl > ngwt ) THEN
          CALL errore( ' wave_descriptor_init ', ' arg no. 2 incompatible with arg no. 3 ', 1 ) 
        ELSE
          desc % ngwl = ngwl
        END IF

        !  bands

        desc % nbt = 0
        DO is = 1, nspin
          IF( nbt( is ) <  0 ) &
            CALL errore( ' wave_descriptor_init ', ' arg no. 5 out of range ', 1 ) 
          desc % nbt( is ) = nbt( is )
        END DO


        desc % nbl = 0
        DO is = 1, nspin
          IF( nbl( is ) <= 0 ) THEN
            desc % nbl( is ) = nbt( is )
          ELSE IF( nbl( is ) > nbt( is ) ) THEN
            CALL errore( ' wave_descriptor_init ', ' arg no. 4 incompatible with arg no. 5 ', 1 ) 
          ELSE
            desc % nbl( is ) = nbl( is )
          END IF
        END DO
        !  k - points

        IF( nkt <  0 ) &
          CALL errore( ' wave_descriptor_init ', ' arg no. 7 out of range ', 1 ) 

        desc % nkt = nkt

        IF( nkl <= 0 ) THEN
          desc % nkl = nkt
        ELSE IF( nkl > nkt ) THEN
          CALL errore( ' wave_descriptor_init ', ' arg no. 6 incompatible with arg no. 7 ', 1 ) 
        ELSE
          desc % nkl = nkl
        END IF

        ! spin

        IF( nspin < 0 .OR. nspin > 2 ) &
          CALL errore( ' wave_descriptor_init ', ' arg no. 8 out of range ', 1 ) 

        desc % nspin = nspin

        ! other

        IF( TRIM( scheme ) == 'gamma' ) THEN
          desc % isym = 0
          desc % gamma = .TRUE.
        ELSE
          desc % isym = 1
          desc % gamma = .FALSE.
        END IF

        desc % gzero = lgz

        desc % ldg = MAX( 1, desc % ngwl  )
        desc % ldb = MAX( 1, MAXVAL( desc % nbl ) )
        desc % ldk = MAX( 1, desc % nkl   )
        desc % lds = MAX( 1, desc % nspin )

        RETURN
      END SUBROUTINE wave_descriptor_init


      SUBROUTINE wave_descriptor_info( desc, nam, iunit )

        IMPLICIT NONE

        TYPE (wave_descriptor), INTENT(IN) :: desc
        INTEGER, INTENT(IN) :: iunit
        CHARACTER(LEN=*) :: nam

        WRITE( iunit, 10 ) nam, desc%ldg, desc%ldb, desc%ldk, desc%lds, &
          desc%ngwl, desc%ngwt, desc%nbl, desc%nbt, desc%nkl, desc%nkt, &
          desc%nspin, desc%isym, desc%gzero

10      FORMAT( 3X, 'Wave function descriptor . . . . . : ',A20,/ &
               ,3X, 'leading dimensions (g,b,k,s) . . . : ',4I8,/ &
               ,3X, 'num. of plane wave (Local, Global) : ',2I8,/&
               ,3X, 'num. of bands (Local, Global). . . : ',4I5,/&
               ,3X, 'num. of k points (Local, Global) . : ',2I5,/&
               ,3X, 'num. of spin . . . . . . . . . . . : ',I4,/&
               ,3X, 'symmetry . . . . . . . . . . . . . : ',I4,/&
               ,3X, 'has G == 0 vector. . . . . . . . . : ',L7)

        RETURN
      END SUBROUTINE wave_descriptor_info

!  ----------------------------------------------  !
      END MODULE wave_types
!  ----------------------------------------------  !
