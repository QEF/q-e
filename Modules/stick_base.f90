!
! Copyright (C) 2002 FPMD group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

!=----------------------------------------------------------------------=
   MODULE stick_base
!=----------------------------------------------------------------------=

      USE kinds
      USE io_global, ONLY: ionode

        IMPLICIT NONE
        PRIVATE
        SAVE

        PUBLIC :: sticks_maps, sticks_sort, sticks_countg, sticks_dist, sticks_pairup
        PUBLIC :: sticks_owner, sticks_deallocate
 
! ...   sticks_owner :   stick owner, sticks_owner( i, j ) is the index of the processor
! ...     (starting from 1) owning the stick whose x and y coordinate  are i and j.

        INTEGER, ALLOCATABLE, TARGET :: sticks_owner( : , : )       

        INTERFACE sticks_dist
          MODULE PROCEDURE sticks_dist1
        END INTERFACE

!=----------------------------------------------------------------------=
   CONTAINS
!=----------------------------------------------------------------------=

      SUBROUTINE  sticks_maps( tk, ub, lb, b1, b2, b3, gcut, gcutw, gcuts, st, stw, sts )

          USE mp, ONLY: mp_sum
          USE mp_global, ONLY: me_pool, nproc_pool, intra_pool_comm

          LOGICAL, INTENT(IN) :: tk    !  if true use the full space grid
          INTEGER, INTENT(IN) :: ub(:) !  upper bounds for i-th grid dimension
          INTEGER, INTENT(IN) :: lb(:) !  lower bounds for i-th grid dimension
          REAL(dbl) , INTENT(IN) :: b1(:), b2(:), b3(:) ! reciprocal space base vectors
          REAL(dbl) , INTENT(IN) :: gcut   ! cut-off for potentials
          REAL(dbl) , INTENT(IN) :: gcutw  ! cut-off for plane waves
          REAL(dbl) , INTENT(IN) :: gcuts  ! cut-off for smooth mesh
          INTEGER, INTENT(OUT) :: st( lb(1): ub(1), lb(2):ub(2) ) ! stick map for potential
          INTEGER, INTENT(OUT) :: stw(lb(1): ub(1), lb(2):ub(2) ) ! stick map for wave functions
          INTEGER, INTENT(OUT) :: sts(lb(1): ub(1), lb(2):ub(2) ) ! stick map for smooth mesh

          INTEGER :: i, j, k, kip
          REAL(dbl) :: gsq

          stw  = 0
          st   = 0
          sts  = 0

! ...       Here find the basic maps of sticks st, stw and sts for the potential
! ...       cut-off gcut, wavefunction cut-off gcutw, and smooth mesh cut-off gcuts

! ...       st(i,j) will contain the number of G vectors of the stick whose
! ...       indices are (i,j). 
 
          IF( .NOT. tk ) THEN

            kip = 0 + ABS(lb(3)) + 1
            IF( MOD( kip, nproc_pool ) == me_pool ) THEN
              st (0,0) = st (0,0) + 1
              stw(0,0) = stw(0,0) + 1
              sts(0,0) = sts(0,0) + 1
            END IF

            DO i= 0, 0 
              DO j= 0, 0 
                DO k= 1, ub(3)
                  kip = k + ABS(lb(3)) + 1
                  IF( MOD( kip, nproc_pool ) == me_pool ) THEN
                    gsq=    (REAL(i)*b1(1)+REAL(j)*b2(1)+REAL(k)*b3(1) )**2
                    gsq=gsq+(REAL(i)*b1(2)+REAL(j)*b2(2)+REAL(k)*b3(2) )**2
                    gsq=gsq+(REAL(i)*b1(3)+REAL(j)*b2(3)+REAL(k)*b3(3) )**2
                    IF(gsq.LE.gcut ) THEN
                      st(i,j) = st(i,j) + 1
                      IF(gsq.LE.gcutw) THEN
                        stw(i,j) = stw(i,j) + 1
                      END IF
                      IF(gsq.LE.gcuts) THEN
                        sts(i,j) = sts(i,j) + 1
                      END IF
                    END IF
                  END IF
                END DO
              END DO
            END DO

            DO i = 0, 0
              DO j = 1, ub(2)
                DO k = lb(3), ub(3)
                  kip = k + ABS(lb(3)) + 1
                  IF( MOD( kip, nproc_pool) == me_pool ) THEN
                    gsq=    (REAL(i)*b1(1)+REAL(j)*b2(1)+REAL(k)*b3(1) )**2
                    gsq=gsq+(REAL(i)*b1(2)+REAL(j)*b2(2)+REAL(k)*b3(2) )**2
                    gsq=gsq+(REAL(i)*b1(3)+REAL(j)*b2(3)+REAL(k)*b3(3) )**2
                    IF(gsq.LE.gcut ) THEN
                      st(i,j) = st(i,j) + 1
                      IF(gsq.LE.gcutw) THEN
                        stw(i,j) = stw(i,j) + 1
                      END IF
                      IF(gsq.LE.gcuts) THEN
                        sts(i,j) = sts(i,j) + 1
                      END IF
                    END IF
                  END IF
                END DO
              END DO
            END DO

            DO i = 1, ub(1)
              DO j = lb(2), ub(2)
                DO k = lb(3), ub(3)
                  kip = k + ABS(lb(3)) + 1
                  IF( MOD( kip, nproc_pool) == me_pool ) THEN
                    gsq=    (REAL(i)*b1(1)+REAL(j)*b2(1)+REAL(k)*b3(1) )**2
                    gsq=gsq+(REAL(i)*b1(2)+REAL(j)*b2(2)+REAL(k)*b3(2) )**2
                    gsq=gsq+(REAL(i)*b1(3)+REAL(j)*b2(3)+REAL(k)*b3(3) )**2
                    IF(gsq.LE.gcut ) THEN
                      st(i,j) = st(i,j) + 1
                      IF(gsq.LE.gcutw) THEN
                        stw(i,j) = stw(i,j) + 1
                      END IF
                      IF(gsq.LE.gcuts) THEN
                        sts(i,j) = sts(i,j) + 1
                      END IF
                    END IF
                  END IF
                END DO
              END DO
            END DO

          ELSE

            DO i= lb(1), ub(1)
              DO j= lb(2), ub(2)
                DO k= lb(3), ub(3)
                  kip = k + ABS(lb(3)) + 1
                  IF( MOD( kip, nproc_pool ) == me_pool ) THEN
                    gsq=    (REAL(i)*b1(1)+REAL(j)*b2(1)+REAL(k)*b3(1) )**2
                    gsq=gsq+(REAL(i)*b1(2)+REAL(j)*b2(2)+REAL(k)*b3(2) )**2
                    gsq=gsq+(REAL(i)*b1(3)+REAL(j)*b2(3)+REAL(k)*b3(3) )**2
                    IF(gsq.LE.gcut ) THEN
                      st(i,j) = st(i,j) + 1
                    END IF
                    IF(gsq.LE.gcutw) THEN
                      stw(i,j) = stw(i,j) + 1
                    END IF
                    IF(gsq.LE.gcuts) THEN
                      sts(i,j) = sts(i,j) + 1
                    END IF
                  END IF
                END DO
              END DO
            END DO

          END IF

          CALL mp_sum(st  ,intra_pool_comm)
          CALL mp_sum(stw ,intra_pool_comm)
          CALL mp_sum(sts ,intra_pool_comm)

! Test sticks
!          write(6,*) 'testtesttesttesttesttesttesttesttesttest'
!          DO i = lb(1), ub(1)
!            DO j = lb(2), ub(2)
!              write(6,'(2I4,I6)') i,j,stw(i,j)
!            END DO
!          END DO
!          write(6,*) 'testtesttesttesttesttesttesttesttesttest'
! Test sticks

        RETURN
      END SUBROUTINE

!=----------------------------------------------------------------------=

      SUBROUTINE sticks_sort( ngc, ngcw, ngcs, nct, index )

! ...     This subroutine sorts the sticks indexes, according to 
! ...     the lenght and type of the sticks, wave functions sticks
! ...     first, then smooth mesh sticks, and finally potential
! ...     sticks

        USE mp_global, ONLY: nproc_pool

        ! lenghts of sticks, ngc for potential mesh, ngcw for wave functions mesh
        ! and ngcs for smooth mesh

        INTEGER, INTENT(IN) :: ngc(:), ngcw(:), ngcs(:) 

        ! nct, total number of sticks

        INTEGER, INTENT(IN) :: nct                      

        ! index, on output, new sticks indexes

        INTEGER, INTENT(OUT) :: index(:)                

        INTEGER :: mc, nr3x, ic
        REAL(dbl), ALLOCATABLE :: aux(:)

        nr3x = MAXVAL( ngc(1:nct) ) + 1

        IF( nproc_pool > 1 ) THEN
          ALLOCATE( aux( nct ) )
          DO mc = 1, nct
            aux(mc) = -(ngcw(mc)*nr3x**2 + ngcs(mc)*nr3x + ngc(mc))
            index(mc) = 0
          END DO
          CALL hpsort( nct, aux(1), index(1))
          DEALLOCATE( aux )
        ELSE
          ic = 0
          do mc = 1, nct
            if( ngcw(mc) > 0 ) then
              ic = ic + 1
              index(ic) = mc
            endif
          end do
          do mc = 1, nct
            if( ngcs(mc) > 0 .AND. ngcw(mc) == 0 ) then
              ic = ic + 1
              index(ic) = mc
            endif
          end do
          do mc = 1, nct
            if( ngc(mc) > 0 .AND. ngcs(mc) == 0 .AND. ngcw(mc) == 0 ) then
              ic = ic + 1
              index(ic) = mc
            endif
          end do
        END IF

        ! WRITE(6,*) '-----------------'
        ! WRITE(6,*) 'STICKS_SORT DEBUG'
        ! DO mc = 1, nct
        !   WRITE(6, fmt="(4I10)" ) index(mc), ngcw( index(mc) ), ngcs( index(mc) ), ngc( index(mc) )
        ! END DO
        ! WRITE(6,*) '-----------------'

        RETURN
      END SUBROUTINE

!=----------------------------------------------------------------------=

    SUBROUTINE sticks_countg( tk, ub, lb, st, stw, sts, in1, in2, ngc, ngcw, ngcs )

      INTEGER, INTENT(IN) :: ub(:), lb(:)
      INTEGER, INTENT(IN) :: st( lb(1): ub(1), lb(2):ub(2) ) ! stick map for potential
      INTEGER, INTENT(IN) :: stw(lb(1): ub(1), lb(2):ub(2) ) ! stick map for wave functions
      INTEGER, INTENT(IN) :: sts(lb(1): ub(1), lb(2):ub(2) ) ! stick map for smooth mesh
      LOGICAL, INTENT(IN) :: tk

      INTEGER, INTENT(OUT) :: in1(:), in2(:)
      INTEGER, INTENT(OUT) :: ngc(:), ngcw(:), ngcs(:)

      INTEGER :: j1, j2, i1, i2, nct, min_size

!
! ...     initialize the sticks indexes array ist
! ...     nct counts columns containing G-vectors for the dense grid
! ...     ncts counts columns contaning G-vectors for the smooth grid
!
      nct   = 0

      ngc   = 0
      ngcs  = 0
      ngcw  = 0

      min_size = MIN( SIZE( in1 ), SIZE( in2 ), SIZE( ngc ), SIZE( ngcw ), SIZE( ngcs ) )

      DO j2 = 0, ( ub(2) - lb(2) )
        DO j1 = 0, ( ub(1) - lb(1) )

          i1 = j1
          if( i1 > ub(1) ) i1 = lb(1) + ( i1 - ub(1) ) - 1

          i2 = j2
          if( i2 > ub(2) ) i2 = lb(2) + ( i2 - ub(2) ) - 1

          IF( st( i1, i2 ) > 0 ) THEN

            ! this sticks contains G-vectors

            nct = nct + 1
            IF( nct > min_size ) &
              CALL errore(' sticks_countg ',' too many sticks ', nct )

            in1(nct) = i1
            in2(nct) = i2

            ngc(nct) = st( i1 , i2)
            IF( stw( i1, i2 ) .GT. 0 ) ngcw(nct) = stw( i1 , i2)
            IF( sts( i1, i2 ) .GT. 0 ) ngcs(nct) = sts( i1 , i2)

          END IF

          ! WRITE(7,fmt="(5I5)") i1, i2, nct, ngc(nct), ngcw( nct )

        END DO
      END DO

      RETURN
    END SUBROUTINE

!=----------------------------------------------------------------------=

    SUBROUTINE sticks_dist1( tk, ub, lb, index, in1, in2, ngc, ngcw, ngcs, nct, &
                             ncp, ncpw, ncps, ngp, ngpw, ngps, stown, stownw, stowns )

      USE mp_global, ONLY: nproc_pool

      LOGICAL, INTENT(IN) :: tk

      INTEGER, INTENT(IN) :: ub(:), lb(:), index(:)
      INTEGER, INTENT(OUT) :: stown( lb(1): ub(1), lb(2):ub(2) ) ! stick map for potential
      INTEGER, INTENT(OUT) :: stownw(lb(1): ub(1), lb(2):ub(2) ) ! stick map for wave functions
      INTEGER, INTENT(OUT) :: stowns(lb(1): ub(1), lb(2):ub(2) ) ! stick map for smooth mesh

      INTEGER, INTENT(IN) :: in1(:), in2(:)
      INTEGER, INTENT(IN) :: ngc(:), ngcw(:), ngcs(:)
      INTEGER, INTENT(IN) :: nct
      INTEGER, INTENT(OUT) :: ncp(:), ncpw(:), ncps(:)
      INTEGER, INTENT(OUT) :: ngp(:), ngpw(:), ngps(:)

      INTEGER :: mc, i1, i2, i, j, jj

      ncp  = 0
      ncps = 0
      ncpw = 0
      ngp  = 0
      ngps = 0
      ngpw = 0

      stown  = 0
      stownw = 0
      stowns = 0

      DO mc = 1, nct

         i = index( mc )
!
! index contains the desired ordering of sticks (see above)
!
         i1 = in1( i )
         i2 = in2( i )
!
         if ( ( .NOT. tk ) .AND. ( (i1 < 0) .or. ( (i1 == 0) .and. (i2 < 0) ) ) ) go to 30
!
         jj = 1

         if ( ngcw(i) > 0 ) then
!
! this is an active sticks: find which processor has currently
! the smallest number of plane waves
!
            do j = 1, nproc_pool
               if ( ngpw(j) < ngpw(jj) ) then
                 jj = j
               else if ( ( ngpw(j) == ngpw(jj) ) .AND. ( ncpw(j) < ncpw(jj) ) ) then
                 jj = j
               end if
            end do

         else
!
! this is an inactive sticks: find which processor has currently
! the smallest number of G-vectors
!
            do j = 1, nproc_pool
               if ( ngp(j) < ngp(jj) ) jj = j
            end do

         end if
!
         ! potential mesh

         ncp(jj) = ncp(jj) + 1
         ngp(jj) = ngp(jj) + ngc(i)
         stown(i1,i2) = jj

         ! smooth mesh

         if ( ngcs(i) > 0 ) then
            ncps(jj) = ncps(jj) + 1
            ngps(jj) = ngps(jj) + ngcs(i)
            stowns(i1,i2) = jj
         endif

         ! wave functions mesh

         if ( ngcw(i) > 0 ) then
            ncpw(jj) = ncpw(jj) + 1
            ngpw(jj) = ngpw(jj) + ngcw(i)
            stownw(i1,i2) = jj
         endif

 30      continue

      END DO

      RETURN
    END SUBROUTINE

!=----------------------------------------------------------------------=
    
    SUBROUTINE sticks_pairup( tk, ub, lb, index, in1, in2, ngc, ngcw, ngcs, nct, &
                             ncp, ncpw, ncps, ngp, ngpw, ngps, stown, stownw, stowns )

      USE mp_global, ONLY: nproc_pool

      LOGICAL, INTENT(IN) :: tk

      INTEGER, INTENT(IN) :: ub(:), lb(:), index(:)
      INTEGER, INTENT(INOUT) :: stown( lb(1): ub(1), lb(2):ub(2) ) ! stick map for potential
      INTEGER, INTENT(INOUT) :: stownw(lb(1): ub(1), lb(2):ub(2) ) ! stick map for wave functions
      INTEGER, INTENT(INOUT) :: stowns(lb(1): ub(1), lb(2):ub(2) ) ! stick map for wave functions

      INTEGER, INTENT(IN) :: in1(:), in2(:)
      INTEGER, INTENT(IN) :: ngc(:), ngcw(:), ngcs(:)
      INTEGER, INTENT(IN) :: nct
      INTEGER, INTENT(OUT) :: ncp(:), ncpw(:), ncps(:)
      INTEGER, INTENT(OUT) :: ngp(:), ngpw(:), ngps(:)

      INTEGER :: mc, i1, i2, i, j, jj, iss, is, ip

      IF ( .NOT. tk ) THEN

        !  when gamma symmetry is used only the sticks of half reciprocal space
        !  are generated, then here we pair-up the sticks with those of the other
        !  half of the space, using the gamma symmetry relation
        !  Note that the total numero of stick "nct" is not modified

        DO mc = 1, nct
           i = index(mc)
           i1 = in1(i)
           i2 = in2(i)
           IF( i1 == 0 .and. i2 == 0 ) THEN
             jj = stown( i1, i2 )
             if( jj > 0 ) ngp( jj ) = ngp( jj ) + ngc( i ) - 1
             jj = stowns( i1, i2 )
             if( jj > 0 ) ngps( jj ) = ngps( jj ) + ngcs( i ) - 1
             jj = stownw( i1, i2 )
             if( jj > 0 ) ngpw( jj ) = ngpw( jj ) + ngcw( i ) - 1
           ELSE
             jj = stown( i1, i2 )
             if( jj > 0 ) then
               stown( -i1, -i2 ) = jj
               ncp( jj ) = ncp( jj ) + 1
               ngp( jj ) = ngp( jj ) + ngc( i )
             end if
             jj = stowns( i1, i2 )
             if( jj > 0 ) then
               stowns( -i1, -i2 ) = jj
               ncps( jj ) = ncps( jj ) + 1
               ngps( jj ) = ngps( jj ) + ngcs( i )
             end if
             jj = stownw( i1, i2 )
             if( jj > 0 ) then
               stownw( -i1, -i2 ) = jj
               ncpw( jj ) = ncpw( jj ) + 1
               ngpw( jj ) = ngpw( jj ) + ngcw( i )
             end if
           END IF
        END DO

      END IF

      IF( ALLOCATED( sticks_owner ) ) DEALLOCATE( sticks_owner )
      ALLOCATE( sticks_owner( lb(1): ub(1), lb(2):ub(2) ) )

      sticks_owner( :, : ) = ABS( stown( :, :) )

      RETURN
    END SUBROUTINE

!=----------------------------------------------------------------------=


    SUBROUTINE sticks_deallocate
      IF( ALLOCATED( sticks_owner ) ) DEALLOCATE( sticks_owner )
      RETURN
    END SUBROUTINE
   

!=----------------------------------------------------------------------=
   END MODULE stick_base
!=----------------------------------------------------------------------=
