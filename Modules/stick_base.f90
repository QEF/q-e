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
        PUBLIC :: sticks_owner, sticks_deallocate, pstickset, sticks_maps_scalar
 
! ...   sticks_owner :   stick owner, sticks_owner( i, j ) is the index of the processor
! ...     (starting from 1) owning the stick whose x and y coordinate  are i and j.

        INTEGER, ALLOCATABLE, TARGET :: sticks_owner( : , : )       

        INTERFACE sticks_dist
          MODULE PROCEDURE sticks_dist1
        END INTERFACE

!=----------------------------------------------------------------------=
   CONTAINS
!=----------------------------------------------------------------------=

      SUBROUTINE sticks_maps( tk, ub, lb, b1, b2, b3, gcut, gcutw, gcuts, st, stw, sts )

          USE mp, ONLY: mp_sum
          USE mp_global, ONLY: me_pool, nproc_pool, intra_pool_comm

          LOGICAL, INTENT(IN) :: tk    !  if true use the full space grid
          INTEGER, INTENT(IN) :: ub(:) !  upper bounds for i-th grid dimension
          INTEGER, INTENT(IN) :: lb(:) !  lower bounds for i-th grid dimension
          REAL(DP) , INTENT(IN) :: b1(:), b2(:), b3(:) ! reciprocal space base vectors
          REAL(DP) , INTENT(IN) :: gcut   ! cut-off for potentials
          REAL(DP) , INTENT(IN) :: gcutw  ! cut-off for plane waves
          REAL(DP) , INTENT(IN) :: gcuts  ! cut-off for smooth mesh
          INTEGER, INTENT(OUT) :: st( lb(1): ub(1), lb(2):ub(2) ) ! stick map for potential
          INTEGER, INTENT(OUT) :: stw(lb(1): ub(1), lb(2):ub(2) ) ! stick map for wave functions
          INTEGER, INTENT(OUT) :: sts(lb(1): ub(1), lb(2):ub(2) ) ! stick map for smooth mesh

          INTEGER :: i, j, k, kip
          REAL(DP) :: gsq

          stw  = 0
          st   = 0
          sts  = 0

! ...       Here find the basic maps of sticks st, stw and sts for the potential
! ...       cut-off gcut, wavefunction cut-off gcutw, and smooth mesh cut-off gcuts

! ...       st(i,j) will contain the number of G vectors of the stick whose
! ...       indices are (i,j). 
 
#if defined (__EKO)
          write(*,*) ! Workaround for EKOPath compiler bug
#endif
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
                    gsq=    (DBLE(i)*b1(1)+DBLE(j)*b2(1)+DBLE(k)*b3(1) )**2
                    gsq=gsq+(DBLE(i)*b1(2)+DBLE(j)*b2(2)+DBLE(k)*b3(2) )**2
                    gsq=gsq+(DBLE(i)*b1(3)+DBLE(j)*b2(3)+DBLE(k)*b3(3) )**2
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
                    gsq=    (DBLE(i)*b1(1)+DBLE(j)*b2(1)+DBLE(k)*b3(1) )**2
                    gsq=gsq+(DBLE(i)*b1(2)+DBLE(j)*b2(2)+DBLE(k)*b3(2) )**2
                    gsq=gsq+(DBLE(i)*b1(3)+DBLE(j)*b2(3)+DBLE(k)*b3(3) )**2
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
                    gsq=    (DBLE(i)*b1(1)+DBLE(j)*b2(1)+DBLE(k)*b3(1) )**2
                    gsq=gsq+(DBLE(i)*b1(2)+DBLE(j)*b2(2)+DBLE(k)*b3(2) )**2
                    gsq=gsq+(DBLE(i)*b1(3)+DBLE(j)*b2(3)+DBLE(k)*b3(3) )**2
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
                    gsq=    (DBLE(i)*b1(1)+DBLE(j)*b2(1)+DBLE(k)*b3(1) )**2
                    gsq=gsq+(DBLE(i)*b1(2)+DBLE(j)*b2(2)+DBLE(k)*b3(2) )**2
                    gsq=gsq+(DBLE(i)*b1(3)+DBLE(j)*b2(3)+DBLE(k)*b3(3) )**2
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

          CALL mp_sum(st  ,intra_pool_comm )
          CALL mp_sum(stw ,intra_pool_comm )
          CALL mp_sum(sts ,intra_pool_comm )

! Test sticks
!          WRITE( 6,*) 'testtesttesttesttesttesttesttesttesttest'
!          WRITE( 6,*) 'lb = ', lb(1), lb(2)
!          WRITE( 6,*) 'ub = ', ub(1), ub(2)
!          WRITE( 6,*) 'counts    = ', COUNT( st > 0 ), COUNT( stw > 0 ), COUNT( sts > 0 )
!          WRITE( 6,*) 'cut-offs  = ', gcut, gcutw, gcuts
!          WRITE( 6,*) 'b1  = ', b1(1:3)
!          WRITE( 6,*) 'b2  = ', b2(1:3)
!          WRITE( 6,*) 'b3  = ', b3(1:3)
!          DO i = lb(1), ub(1)
!            DO j = lb(2), ub(2)
!              WRITE( 6,'(2I4,I6)') i,j,stw(i,j)
!            END DO
!          END DO
!          WRITE( 6,*) 'testtesttesttesttesttesttesttesttesttest'
! Test sticks

        RETURN
      END SUBROUTINE sticks_maps

!=----------------------------------------------------------------------=

  SUBROUTINE sticks_maps_scalar( lgamma, ub, lb, b1, b2, b3, gcutm, gkcut, gcutms, stw, ngm, ngms )

    LOGICAL, INTENT(IN) :: lgamma !  if true use gamma point simmetry
    INTEGER, INTENT(IN) :: ub(:)  !  upper bounds for i-th grid dimension
    INTEGER, INTENT(IN) :: lb(:)  !  lower bounds for i-th grid dimension
    REAL(DP) , INTENT(IN) :: b1(:), b2(:), b3(:) ! reciprocal space base vectors
    REAL(DP) , INTENT(IN) :: gcutm  ! cut-off for potentials
    REAL(DP) , INTENT(IN) :: gkcut  ! cut-off for plane waves
    REAL(DP) , INTENT(IN) :: gcutms  ! cut-off for smooth mesh
    !
    INTEGER, INTENT(OUT) :: ngm, ngms
    !
    !     stick map for wave functions, note that map is taken in YZ plane
    !
    INTEGER, INTENT(OUT) :: stw( lb(2) : ub(2), lb(3) : ub(3) ) 

    INTEGER :: i1, i2, i3, n1, n2, n3
    REAL(DP) :: amod

    ngm = 0
    ngms = 0

    n1 = MAX( ABS( lb(1) ), ABS( ub(1) ) )
    n2 = MAX( ABS( lb(2) ), ABS( ub(2) ) )
    n3 = MAX( ABS( lb(3) ), ABS( ub(3) ) )

    loop1: do i1 = - n1, n1
       !
       ! Gamma-only: exclude space with x<0
       !
       if (lgamma .and. i1 < 0) cycle loop1
       !
       loop2: do i2 = - n2, n2
          !
          ! Gamma-only: exclude plane with x=0, y<0
          !
          if(lgamma .and. i1 == 0.and. i2 < 0) cycle loop2
          !
          loop3: do i3 = - n3, n3
             !
             ! Gamma-only: exclude line with x=0, y=0, z<0
             !
             if(lgamma .and. i1 == 0 .and. i2 == 0 .and. i3 < 0) cycle loop3
             !
             amod = (i1 * b1 (1) + i2 * b2 (1) + i3 * b3 (1) ) **2 + &
                    (i1 * b1 (2) + i2 * b2 (2) + i3 * b3 (2) ) **2 + &
                    (i1 * b1 (3) + i2 * b2 (3) + i3 * b3 (3) ) **2
             if (amod <= gcutm)  ngm  = ngm  + 1
             if (amod <= gcutms) ngms = ngms + 1
             if (amod <= gkcut ) then
                stw( i2, i3 ) = 1
                if (lgamma) stw( -i2, -i3 ) = 1
             end if
          enddo loop3
       enddo loop2
    enddo loop1

    RETURN
  END SUBROUTINE sticks_maps_scalar


!=----------------------------------------------------------------------=

      SUBROUTINE sticks_sort( ngc, ngcw, ngcs, nct, idx )

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

        INTEGER, INTENT(OUT) :: idx(:)                

        INTEGER :: mc, nr3x, ic
        REAL(DP), ALLOCATABLE :: aux(:)

        nr3x = MAXVAL( ngc(1:nct) ) + 1

        IF( nproc_pool > 1 ) THEN
          ALLOCATE( aux( nct ) )
          DO mc = 1, nct
            aux(mc) = -(ngcw(mc)*nr3x**2 + ngcs(mc)*nr3x + ngc(mc))
            idx(mc) = 0
          END DO
          CALL hpsort( nct, aux(1), idx(1))
          DEALLOCATE( aux )
        ELSE
          ic = 0
          do mc = 1, nct
            if( ngcw(mc) > 0 ) then
              ic = ic + 1
              idx(ic) = mc
            endif
          end do
          do mc = 1, nct
            if( ngcs(mc) > 0 .AND. ngcw(mc) == 0 ) then
              ic = ic + 1
              idx(ic) = mc
            endif
          end do
          do mc = 1, nct
            if( ngc(mc) > 0 .AND. ngcs(mc) == 0 .AND. ngcw(mc) == 0 ) then
              ic = ic + 1
              idx(ic) = mc
            endif
          end do
        END IF

        ! WRITE( stdout,*) '-----------------'
        ! WRITE( stdout,*) 'STICKS_SORT DEBUG'
        ! DO mc = 1, nct
        !   WRITE( stdout, fmt="(4I10)" ) idx(mc), ngcw( idx(mc) ), ngcs( idx(mc) ), ngc( idx(mc) )
        ! END DO
        ! WRITE( stdout,*) '-----------------'

        RETURN
      END SUBROUTINE sticks_sort

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
    END SUBROUTINE sticks_countg

!=----------------------------------------------------------------------=

    SUBROUTINE sticks_dist1( tk, ub, lb, idx, in1, in2, ngc, ngcw, ngcs, nct, &
                             ncp, ncpw, ncps, ngp, ngpw, ngps, stown, stownw, stowns )

      USE mp_global, ONLY: nproc_pool

      LOGICAL, INTENT(IN) :: tk

      INTEGER, INTENT(IN) :: ub(:), lb(:), idx(:)
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

         i = idx( mc )
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
    END SUBROUTINE sticks_dist1

!=----------------------------------------------------------------------=
    
    SUBROUTINE sticks_pairup( tk, ub, lb, idx, in1, in2, ngc, ngcw, ngcs, nct, &
                             ncp, ncpw, ncps, ngp, ngpw, ngps, stown, stownw, stowns )

      USE mp_global, ONLY: nproc_pool

      LOGICAL, INTENT(IN) :: tk

      INTEGER, INTENT(IN) :: ub(:), lb(:), idx(:)
      INTEGER, INTENT(INOUT) :: stown( lb(1): ub(1), lb(2):ub(2) ) ! stick map for potential
      INTEGER, INTENT(INOUT) :: stownw(lb(1): ub(1), lb(2):ub(2) ) ! stick map for wave functions
      INTEGER, INTENT(INOUT) :: stowns(lb(1): ub(1), lb(2):ub(2) ) ! stick map for wave functions

      INTEGER, INTENT(IN) :: in1(:), in2(:)
      INTEGER, INTENT(IN) :: ngc(:), ngcw(:), ngcs(:)
      INTEGER, INTENT(IN) :: nct
      INTEGER, INTENT(OUT) :: ncp(:), ncpw(:), ncps(:)
      INTEGER, INTENT(OUT) :: ngp(:), ngpw(:), ngps(:)

      INTEGER :: mc, i1, i2, i, jj

      IF ( .NOT. tk ) THEN

        !  when gamma symmetry is used only the sticks of half reciprocal space
        !  are generated, then here we pair-up the sticks with those of the other
        !  half of the space, using the gamma symmetry relation
        !  Note that the total numero of stick "nct" is not modified

        DO mc = 1, nct
           i = idx(mc)
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
    END SUBROUTINE sticks_pairup

!=----------------------------------------------------------------------=


        SUBROUTINE pstickset( dfftp, dffts, alat, a1, a2, a3, gcut, gkcut, gcuts, &
          nr1, nr2, nr3, nr1x, nr2x, nr3x, nr1s, nr2s, nr3s, nr1sx, nr2sx, nr3sx, &
          ngw, ngm, ngs )

          USE kinds, ONLY: DP
          USE mp_global, ONLY: me_pool, nproc_pool, intra_pool_comm, nogrp
          USE control_flags, ONLY: gamma_only
          USE io_global, ONLY: ionode
          USE io_global, ONLY: stdout
          USE fft_types, ONLY: fft_dlay_descriptor, fft_dlay_allocate, fft_dlay_set, &
               fft_dlay_scalar


          TYPE(fft_dlay_descriptor), INTENT(INOUT) :: dfftp, dffts
          REAL(DP), INTENT(IN) :: a1(3), a2(3), a3(3), alat
          REAL(DP), INTENT(IN) :: gcut, gkcut, gcuts
          INTEGER, INTENT(IN) :: nr1, nr2, nr3, nr1x, nr2x, nr3x
          INTEGER, INTENT(IN) :: nr1s, nr2s, nr3s, nr1sx, nr2sx, nr3sx
          INTEGER, INTENT(OUT) :: ngw, ngm, ngs

          LOGICAL :: tk
! ...     tk  logical flag, TRUE if the symulation does not have the
! ...         GAMMA symmetry

          REAL(DP) :: b1(3), b2(3), b3(3)
! ...     b1, b2, b3 reciprocal space base vectors.

          INTEGER :: ub(3), lb(3)
! ...     ub(i), lb(i) upper and lower miller indexes

!
! ...     Plane Waves
!

        INTEGER, ALLOCATABLE :: stw(:,:)
! ...   stick map (wave functions), stw(i,j) = number of G-vector in the
! ...     stick whose x and y miller index are i and j

        INTEGER, ALLOCATABLE :: nstpw(:)
! ...   number of sticks (wave functions), nstpw(ip) = number of stick
! ...     for processor ip

        INTEGER, ALLOCATABLE :: sstpw(:)
! ...   number of G-vectors (wave functions), sstpw(ip) = sum of the
! ...     sticks lenght for processor ip = number of G-vectors
! ...     owned by the processor ip

        INTEGER :: nstw, nstpwx
! ...   nstw     local number of sticks (wave functions)
! ...   nstpwx   maximum among all processors of nstw

!
! ...     Potentials
!

        INTEGER, ALLOCATABLE :: st(:,:)
! ...   stick map (potentials), st(i,j) = number of G-vector in the
! ...     stick whose x and y miller index are i and j

        INTEGER, ALLOCATABLE :: nstp(:)
! ...   number of sticks (potentials), nstp(ip) = number of stick
! ...     for processor ip

        INTEGER, ALLOCATABLE :: sstp(:)
! ...   number of G-vectors (potentials), sstp(ip) = sum of the
! ...     sticks lenght for processor ip = number of G-vectors
! ...     owned by the processor ip

        INTEGER :: nst, nstpx
! ...   nst      local number of sticks (potentials)
! ...   nstpx    maximum among all processors of nst

!
! ...     Smooth Mesh
!

        INTEGER, ALLOCATABLE :: sts(:,:)
! ...   stick map (smooth mesh), sts(i,j) = number of G-vector in the
! ...     stick whose x and y miller index are i and j

        INTEGER, ALLOCATABLE :: nstps(:)
! ...   number of sticks (smooth mesh), nstp(ip) = number of stick
! ...     for processor ip

        INTEGER, ALLOCATABLE :: sstps(:)
! ...   number of G-vectors (smooth mesh), sstps(ip) = sum of the
! ...     sticks lenght for processor ip = number of G-vectors
! ...     owned by the processor ip

        INTEGER :: nsts
! ...   nsts      local number of sticks (smooth mesh)


        INTEGER, ALLOCATABLE :: ist(:,:)    ! sticks indexes ordered



          INTEGER :: ip, ngm_ , ngs_
          INTEGER, ALLOCATABLE :: idx(:)

          tk    = .NOT. gamma_only
          ub(1) = ( nr1 - 1 ) / 2
          ub(2) = ( nr2 - 1 ) / 2
          ub(3) = ( nr3 - 1 ) / 2
          lb    = - ub

          ! ... reciprocal lattice generators

          CALL recips( a1, a2, a3, b1, b2, b3 )
          b1 = b1 * alat
          b2 = b2 * alat
          b3 = b3 * alat

          ! ...       Allocate maps

          ALLOCATE( stw ( lb(1):ub(1), lb(2):ub(2) ) )
          ALLOCATE( st  ( lb(1):ub(1), lb(2):ub(2) ) )
          ALLOCATE( sts ( lb(1):ub(1), lb(2):ub(2) ) )

          st  = 0
          stw = 0
          sts = 0

! ...       Fill in the stick maps, for given g-space base (b1,b2,b3)  and cut-off

          CALL sticks_maps( tk, ub, lb, b1, b2, b3, gcut, gkcut, gcuts, st, stw, sts )

! ...       Now count the number of stick nst and nstw

          nst  = COUNT( st  > 0 )
          nstw = COUNT( stw > 0 )
          nsts = COUNT( sts > 0 )

          IF (ionode) THEN
            WRITE( stdout,*)
            WRITE( stdout,10)
 10         FORMAT(3X,'Stick Mesh',/, &
                   3X,'----------')
            WRITE( stdout,15) nst, nstw, nsts
 15         FORMAT( 3X, 'nst =', I6, ',  nstw =', I6, ', nsts =', I6 )
          END IF

          ALLOCATE(ist(nst,5))

          ALLOCATE(nstp(nproc_pool))
          ALLOCATE(sstp(nproc_pool))

          ALLOCATE(nstpw(nproc_pool))
          ALLOCATE(sstpw(nproc_pool))

          ALLOCATE(nstps(nproc_pool))
          ALLOCATE(sstps(nproc_pool))

! ...       initialize the sticks indexes array ist

          CALL sticks_countg( tk, ub, lb, st, stw, sts, &
            ist(:,1), ist(:,2), ist(:,4), ist(:,3), ist(:,5) )

! ...       Sorts the sticks according to their lenght

          ALLOCATE( idx( nst ) )

          CALL sticks_sort( ist(:,4), ist(:,3), ist(:,5), nst, idx )

          ! ... Set as first stick the stick containing the G=0
          !
          !  DO iss = 1, nst
          !    IF( ist( idx( iss ), 1 ) == 0 .AND. ist( idx( iss ), 2 ) == 0 )  EXIT
          !  END DO
          !  itmp         = idx( 1 )
          !  idx( 1 )   = idx( iss )
          !  idx( iss ) = itmp

          CALL sticks_dist( tk, ub, lb, idx, ist(:,1), ist(:,2), ist(:,4), ist(:,3), ist(:,5), &
             nst, nstp, nstpw, nstps, sstp, sstpw, sstps, st, stw, sts )

          ngw = sstpw( me_pool + 1 )
          ngm = sstp( me_pool + 1 )
          ngs = sstps( me_pool + 1 )

          CALL sticks_pairup( tk, ub, lb, idx, ist(:,1), ist(:,2), ist(:,4), ist(:,3), ist(:,5), &
             nst, nstp, nstpw, nstps, sstp, sstpw, sstps, st, stw, sts )

          ! ...   Allocate and Set fft data layout descriptors

#if defined __PARA

          CALL fft_dlay_allocate( dfftp, nproc_pool, nr1x,  nr2x )
          CALL fft_dlay_allocate( dffts, nproc_pool, nr1sx, nr2sx )

          CALL fft_dlay_set( dfftp, tk, nst, nr1, nr2, nr3, nr1x, nr2x, nr3x, (me_pool+1), &
            nproc_pool, nogrp, ub, lb, idx, ist(:,1), ist(:,2), nstp, nstpw, sstp, sstpw, st, stw )
          CALL fft_dlay_set( dffts, tk, nsts, nr1s, nr2s, nr3s, nr1sx, nr2sx, nr3sx, (me_pool+1), &
            nproc_pool, nogrp, ub, lb, idx, ist(:,1), ist(:,2), nstps, nstpw, sstps, sstpw, sts, stw )

#else

          DEALLOCATE( stw )
          ALLOCATE( stw( lb(2) : ub(2), lb(3) : ub(3) ) )

          CALL sticks_maps_scalar( (.not.tk), ub, lb, b1, b2, b3, gcut, gkcut, gcuts, stw, ngm_ , ngs_ )

          IF( ngm_ /= ngm ) CALL errore( ' pstickset ', ' inconsistent ngm ', ABS( ngm - ngm_ ) )
          IF( ngs_ /= ngs ) CALL errore( ' pstickset ', ' inconsistent ngs ', ABS( ngs - ngs_ ) )

          CALL fft_dlay_allocate( dfftp, nproc_pool, MAX(nr1x, nr3x),  nr2x  )
          CALL fft_dlay_allocate( dffts, nproc_pool, MAX(nr1sx, nr3sx), nr2sx )

          CALL fft_dlay_scalar( dfftp, ub, lb, nr1, nr2, nr3, nr1x, nr2x, nr3x, stw )
          CALL fft_dlay_scalar( dffts, ub, lb, nr1s, nr2s, nr3s, nr1sx, nr2sx, nr3sx, stw )

#endif

! ...     Maximum number of sticks (potentials)
          nstpx  = MAXVAL( nstp )
! ...     Maximum number of sticks (wave func.)
          nstpwx = MAXVAL( nstpw  )

          IF (ionode) WRITE( stdout,119)
 119      FORMAT(3X,'     PEs    n.st   n.stw   n.sts    n.g    n.gw   n.gs')
          DO ip = 1, nproc_pool
            IF (ionode) THEN
              WRITE( stdout,120) ip, nstp(ip),  nstpw(ip), nstps(ip), sstp(ip), sstpw(ip), sstps(ip)
            END IF
          END DO
          IF (ionode) THEN
            WRITE( stdout,120) 0, SUM(nstp),  SUM(nstpw), SUM(nstps), SUM(sstp), SUM(sstpw), SUM(sstps)
          END IF
 120      FORMAT(3X,7I8)


          DEALLOCATE( ist )
          DEALLOCATE( idx )

          DEALLOCATE( st, stw, sts )
          DEALLOCATE( sstp )
          DEALLOCATE( nstp )
          DEALLOCATE( sstpw )
          DEALLOCATE( nstpw )
          DEALLOCATE( sstps )
          DEALLOCATE( nstps )

          IF(ionode) WRITE( stdout,*)

          RETURN
        END SUBROUTINE pstickset


!=----------------------------------------------------------------------=

    SUBROUTINE sticks_deallocate
      IF( ALLOCATED( sticks_owner ) ) DEALLOCATE( sticks_owner )
      RETURN
    END SUBROUTINE sticks_deallocate
   

!=----------------------------------------------------------------------=
   END MODULE stick_base
!=----------------------------------------------------------------------=
