!=----------------------------------------------------------------------=
MODULE stick_base
!=----------------------------------------------------------------------=

   USE fft_param
   IMPLICIT NONE
   PRIVATE
   SAVE

   PUBLIC :: sticks_map_set
   PUBLIC :: sticks_map_index, sticks_sort_new, sticks_dist_new
   PUBLIC :: sticks_map, sticks_map_allocate
   PUBLIC :: sticks_map_deallocate, get_sticks

   TYPE sticks_map
      LOGICAL :: lgamma=.false. ! if .true. the map has gamma symmetry
      LOGICAL :: lpara=.false.  ! if .true. the map is set for parallel and serial, if .false. only serial 
      INTEGER :: mype=0   ! my task id (starting from 0)
      INTEGER :: nproc=1  ! number of task (as nproc in fft_type_descriptor)
      INTEGER :: nyfft=1  ! number processors in y-direction (as nproc2 in fft_type_descriptor)
      INTEGER, ALLOCATABLE :: iproc(:,:) ! the processor index (as in fft_type_descriptor)
      INTEGER, ALLOCATABLE :: iproc2(:)  ! the Y group processor index (as in fft_type_descriptor)
#if defined(__MPI)
      INTEGER :: comm     = MPI_COMM_NULL
#else
      INTEGER :: comm     = 0          ! communicator of the fft gruop 
#endif
      INTEGER :: nstx=0   ! a safe maximum number of sticks on the map
      INTEGER :: lb(3)=0  ! map's lower bounds
      INTEGER :: ub(3)=0  ! map's upper bounds
      INTEGER, ALLOCATABLE :: idx(:)   ! the index of each stick
      INTEGER, ALLOCATABLE :: ist(:,:) ! the cartesian coordinates of each stick
      INTEGER, ALLOCATABLE :: stown(:,:) ! the owner of each stick
      INTEGER, ALLOCATABLE :: indmap(:,:) ! the index of each stick (represented on the map)
      REAL(DP) :: bg(3,3) ! base vectors, the generators of the mapped space
   END TYPE

!=----------------------------------------------------------------------=
CONTAINS
!=----------------------------------------------------------------------=
   SUBROUTINE sticks_map_deallocate( smap )
      IMPLICIT NONE
      TYPE( sticks_map ) :: smap
      !write (6,*) ' inside sticks_map_deallocate'; FLUSH(6)
      IF( ALLOCATED( smap%iproc ) ) DEALLOCATE( smap%iproc )
      IF( ALLOCATED( smap%iproc2 ) ) DEALLOCATE( smap%iproc2 )
      IF( ALLOCATED( smap%idx ) ) DEALLOCATE( smap%idx )
      IF( ALLOCATED( smap%ist ) ) DEALLOCATE( smap%ist )
      IF( ALLOCATED( smap%stown ) ) DEALLOCATE( smap%stown )
      IF( ALLOCATED( smap%indmap ) ) DEALLOCATE( smap%indmap )
      smap%ub = 0
      smap%lb = 0
      smap%nstx = 0
   END SUBROUTINE sticks_map_deallocate

   SUBROUTINE sticks_map_allocate( smap, lgamma, lpara, nyfft, iproc, iproc2, nr1, nr2, nr3, bg, comm )
   ! assign 
      IMPLICIT NONE
      TYPE( sticks_map ) :: smap
      LOGICAL, INTENT(IN) :: lgamma
      LOGICAL, INTENT(IN) :: lpara
      INTEGER, INTENT(IN) :: nyfft
      INTEGER, INTENT(IN) :: iproc(:,:)
      INTEGER, INTENT(IN) :: iproc2(:)
      INTEGER, INTENT(IN) :: nr1, nr2, nr3
      INTEGER, INTENT(IN) :: comm
      REAL(DP), INTENT(IN) :: bg(3,3)
      INTEGER :: lb(3), ub(3)
      INTEGER :: nzfft, nstx, ierr
      INTEGER, ALLOCATABLE :: indmap(:,:), stown(:,:), idx(:), ist(:,:)
      !write (6,*) ' inside sticks_map_allocate'; FLUSH(6)
      ub(1) = ( nr1 - 1 ) / 2
      ub(2) = ( nr2 - 1 ) / 2
      ub(3) = ( nr3 - 1 ) / 2
      lb    = - ub
      !write(*,*) 'ub,lb',ub,lb
      nstx = (ub(1)-lb(1)+1)*(ub(2)-lb(2)+1) ! we stay very large indeed
      !write(*,*) smap%nstx,nstx
      IF( smap%nstx == 0 ) THEN
         ! this map is clean, allocate
         !write (*,*) ' ! this map is clean, allocate '
         smap%mype = 0
         smap%nproc = 1
         smap%comm = comm
#if defined(__MPI)
         CALL MPI_COMM_RANK( smap%comm, smap%mype, ierr )
         CALL MPI_COMM_SIZE( smap%comm, smap%nproc, ierr )
#endif
         smap%lgamma = lgamma
         smap%lpara = lpara
         smap%comm = comm
         smap%nstx = nstx
         smap%ub = ub
         smap%lb = lb
         smap%bg = bg
         smap%nyfft = nyfft
         nzfft = smap%nproc / nyfft
         ALLOCATE ( smap%iproc( nyfft, nzfft ), smap%iproc2 ( smap%nproc ) )
         smap%iproc = iproc
         smap%iproc2= iproc2

         IF( ALLOCATED( smap%indmap ) ) THEN
            CALL fftx_error__(' sticks_map_allocate ',' indmap already allocated ', 1 )
         END IF
         IF( ALLOCATED( smap%stown ) ) THEN
            CALL fftx_error__(' sticks_map_allocate ',' stown already allocated ', 1 )
         END IF
         IF( ALLOCATED( smap%idx ) ) THEN
            CALL fftx_error__(' sticks_map_allocate ',' idx already allocated ', 1 )
         END IF
         IF( ALLOCATED( smap%ist ) ) THEN
            CALL fftx_error__(' sticks_map_allocate ',' ist already allocated ', 1 )
         END IF
         ALLOCATE( smap%indmap ( lb(1):ub(1), lb(2):ub(2) ) )
         ALLOCATE( smap%stown ( lb(1):ub(1), lb(2):ub(2) ) )
         ALLOCATE( smap%idx( nstx ) )
         ALLOCATE( smap%ist( nstx , 2) )
         smap%stown = 0
         smap%indmap = 0
         smap%idx = 0
         smap%ist = 0
      ELSE IF( smap%nstx < nstx .OR. smap%ub(3) < ub(3) ) THEN
         !  change the size of the map, but keep the data already there
         !write (*,*) ' !  change the size of the map, but keep the data already there '
         IF( smap%lgamma .neqv.  lgamma ) THEN
            CALL fftx_error__(' sticks_map_allocate ',' changing gamma symmetry not allowed ', 1 )
         END IF
         IF( smap%comm /= comm ) THEN
            CALL fftx_error__(' sticks_map_allocate ',' changing communicator not allowed ', 1 )
         END IF
         ALLOCATE( indmap ( lb(1):ub(1), lb(2):ub(2) ) )
         ALLOCATE( stown ( lb(1):ub(1), lb(2):ub(2) ) )
         ALLOCATE( idx( nstx ) )
         ALLOCATE( ist( nstx , 2) )
         idx = 0
         ist = 0 
         indmap = 0
         stown  = 0
         idx( 1:smap%nstx )      = smap%idx
         ist( 1:smap%nstx, : ) = smap%ist 
         indmap( smap%lb(1):smap%ub(1), smap%lb(2):smap%ub(2) ) = smap%indmap( smap%lb(1):smap%ub(1), smap%lb(2):smap%ub(2) )
         stown( smap%lb(1):smap%ub(1), smap%lb(2):smap%ub(2) ) = smap%stown( smap%lb(1):smap%ub(1), smap%lb(2):smap%ub(2) )
         DEALLOCATE( smap%indmap )
         DEALLOCATE( smap%stown )
         DEALLOCATE( smap%idx )
         DEALLOCATE( smap%ist )
         ALLOCATE( smap%indmap ( lb(1):ub(1), lb(2):ub(2) ) )
         ALLOCATE( smap%stown ( lb(1):ub(1), lb(2):ub(2) ) )
         ALLOCATE( smap%idx( nstx ) )
         ALLOCATE( smap%ist( nstx , 2) )
         smap%indmap = indmap
         smap%stown = stown
         smap%idx = idx
         smap%ist = ist
         DEALLOCATE( indmap )
         DEALLOCATE( stown )
         DEALLOCATE( idx )
         DEALLOCATE( ist )
         smap%nstx = nstx
         smap%ub = ub
         smap%lb = lb
         smap%bg = bg
         smap%nyfft = nyfft
         smap%iproc = iproc
         smap%iproc2= iproc2
      ELSE
         !write(*,*) 'map is the same ... should we update ub,lb ?'
         !write(*,*) smap%ub,smap%lb
         IF( smap%lgamma .neqv. lgamma ) THEN
            CALL fftx_error__(' sticks_map_allocate ',' changing gamma symmetry not allowed ', 2 )
         END IF
         IF( smap%comm /= comm ) THEN
            CALL fftx_error__(' sticks_map_allocate ',' changing communicator not allowed ', 1 )
         END IF
      END IF
      RETURN
   END SUBROUTINE sticks_map_allocate
 
   SUBROUTINE sticks_map_set( lgamma, ub, lb, bg, gcut, st, comm )

      ! .. Compute the basic maps of sticks
      ! .. st(i,j) will contain the number of G vectors of the stick whose indices are (i,j).

      IMPLICIT NONE

      LOGICAL, INTENT(in) :: lgamma !  if true use gamma point simmetry
      INTEGER, INTENT(in) :: ub(:)  !  upper bounds for i-th grid dimension
      INTEGER, INTENT(in) :: lb(:)  !  lower bounds for i-th grid dimension
      REAL(DP) , INTENT(in) :: bg(:,:) ! reciprocal space base vectors
      REAL(DP) , INTENT(in) :: gcut  ! cut-off for potentials
      INTEGER, OPTIONAL, INTENT(in) :: comm ! communicator of the g-vec group
      !
      !     stick map for wave functions, note that map is taken in YZ plane
      !
      INTEGER, INTENT(out) :: st( lb(1) : ub(1), lb(2) : ub(2) )
      REAL(DP) :: b1(3), b2(3), b3(3) 
      INTEGER :: i1, i2, i3, n1, n2, n3, mype, nproc, ierr
      REAL(DP) :: amod
      INTEGER :: ngm
      !write (6,*) 'inside sticks_map_set gcut=',gcut; FLUSH(6)
      !write (6,*) ub,lb

      st = 0
      b1(:) = bg(:,1)
      b2(:) = bg(:,2)
      b3(:) = bg(:,3)

      n1 = max( abs( lb(1) ), abs( ub(1) ) )
      n2 = max( abs( lb(2) ), abs( ub(2) ) )
      n3 = max( abs( lb(3) ), abs( ub(3) ) )

      mype = 0
      nproc = 1
#if defined(__MPI)
      IF( PRESENT( comm ) ) THEN
         CALL MPI_COMM_RANK( comm, mype, ierr )
         CALL MPI_COMM_SIZE( comm, nproc, ierr )
      END IF
#endif
      !write (6,*) ' N1,N2,N3 ', n1,n2,n3
      ngm=0
      loop1: DO i1 = - n1, n1
         !
         ! Gamma-only: exclude space with x<0
         !
         IF ( (lgamma .and. i1 < 0) .OR. ( MOD( i1 + n1, nproc ) /= mype )) CYCLE loop1
         !
         loop2: DO i2 = - n2, n2
            !
            ! Gamma-only: exclude plane with x=0, y<0
            !
            IF(lgamma .and. i1 == 0.and. i2 < 0) CYCLE loop2
            !
            loop3: DO i3 = - n3, n3
               !
               ! Gamma-only: exclude line with x=0, y=0, z<0
               !
               IF(lgamma .and. i1 == 0 .and. i2 == 0 .and. i3 < 0) CYCLE loop3
               !
               amod = (i1 * b1 (1) + i2 * b2 (1) + i3 * b3 (1) ) **2 + &
                      (i1 * b1 (2) + i2 * b2 (2) + i3 * b3 (2) ) **2 + &
                      (i1 * b1 (3) + i2 * b2 (3) + i3 * b3 (3) ) **2
               IF (amod <= gcut ) then
                  st( i1, i2 ) = st( i1, i2 ) + 1
                  ngm =ngm + 1
               endif
            ENDDO loop3
         ENDDO loop2
      ENDDO loop1

#if defined(__MPI)
      IF( PRESENT( comm ) ) THEN
         CALL MPI_ALLREDUCE(MPI_IN_PLACE, st, size(st), MPI_INTEGER, MPI_SUM, comm, ierr)
         CALL MPI_ALLREDUCE(MPI_IN_PLACE, ngm, 1, MPI_INTEGER, MPI_SUM, comm, ierr)
      END IF
#endif

      !write (6,*) ' NGM ', ngm


      RETURN
   END SUBROUTINE sticks_map_set

!=----------------------------------------------------------------------=

   SUBROUTINE sticks_map_index( ub, lb, st, in1, in2, ngc, index_map )

      IMPLICIT NONE

      INTEGER, INTENT(in) :: ub(:), lb(:)
      INTEGER, INTENT(in) :: st( lb(1): ub(1), lb(2):ub(2) ) ! stick map for potential
      INTEGER, INTENT(inout) :: index_map( lb(1): ub(1), lb(2):ub(2) ) ! keep track of sticks index

      INTEGER, INTENT(out) :: in1(:), in2(:)
      INTEGER, INTENT(out) :: ngc(:)

      INTEGER :: j1, j2, i1, i2, i3, nct, min_size, ind
      LOGICAL :: ok
      !write (6,*) ' inside sticks_map_index'; FLUSH(6)

!
! ...     initialize the sticks indexes array ist
! ...     nct counts columns containing G-vectors for the dense grid
! ...     ncts counts columns contaning G-vectors for the smooth grid
!
      nct   = MAXVAL( index_map )
      ngc   = 0

      min_size = min( size( in1 ), size( in2 ), size( ngc ) )

      DO j2 = 0, ( ub(2) - lb(2) )
        DO j1 = 0, ( ub(1) - lb(1) )
          i1 = j1
          IF( i1 > ub(1) ) i1 = lb(1) + ( i1 - ub(1) ) - 1
          i2 = j2
          IF( i2 > ub(2) ) i2 = lb(2) + ( i2 - ub(2) ) - 1
          IF( st( i1, i2 ) > 0 ) THEN
            IF( index_map( i1, i2 ) == 0 ) THEN
              nct = nct + 1
              index_map( i1, i2 ) = nct
            END IF
            ind = index_map( i1, i2 )
            IF( nct > min_size ) &
              CALL fftx_error__(' sticks_map_index ',' too many sticks ', nct )
            in1(ind) = i1
            in2(ind) = i2
            ngc(ind) = st( i1 , i2)
          ENDIF
        ENDDO
      ENDDO

      RETURN
   END SUBROUTINE sticks_map_index

!=----------------------------------------------------------------------=

   SUBROUTINE sticks_sort_new( parallel, ng, nct, idx )

! ...     This subroutine sorts the sticks indexes, according to
! ...     the length and type of the sticks, wave functions sticks
! ...     first, then smooth mesh sticks, and finally potential
! ...     sticks

      ! lengths of sticks, ngc for potential mesh, ngcw for wave functions mesh
      ! and ngcs for smooth mesh

      IMPLICIT NONE

      LOGICAL, INTENT(in) :: parallel
      INTEGER, INTENT(in) :: ng(:)

      ! nct, total number of sticks

      INTEGER, INTENT(in) :: nct

      ! index, on output, new sticks indexes

      INTEGER, INTENT(inout) :: idx(:)

      INTEGER  :: mc, ic, nc
      INTEGER, ALLOCATABLE :: iaux(:)
      INTEGER, ALLOCATABLE :: itmp(:)
      REAL(DP), ALLOCATABLE :: aux(:)
      !write (6,*) ' inside sticks_map_sort_new'; FLUSH(6)

      !  we need to avoid sorting elements already sorted previously
      !  build inverse indexes
      ALLOCATE( iaux( nct ) )
      iaux = 0
      DO mc = 1, nct
         IF( idx( mc ) > 0 ) iaux( idx( mc ) ) = mc
      END DO
      !
      !  check idx has no "hole"
      !
      IF( idx( 1 ) == 0 ) THEN
         ic = 0
         DO mc = 2, nct
            IF( idx( mc ) /= 0 ) THEN
               CALL fftx_error__(' sticks_sort ',' non contiguous indexes 1 ', nct )
            END IF
         END DO
      ELSE
         ic = 1
         DO mc = 2, nct
            IF( idx( mc ) == 0 ) EXIT 
            ic = ic + 1
         END DO
         DO mc = ic+1, nct
            IF( idx( mc ) /= 0 ) THEN
               CALL fftx_error__(' sticks_sort ',' non contiguous indexes 2 ', nct )
            END IF
         END DO
      END IF

      IF( parallel ) THEN
         ALLOCATE( aux( nct ) )
         ALLOCATE( itmp( nct ) )
         itmp = 0
         nc = 0
         DO mc = 1, nct
            IF( ng( mc ) > 0 .AND. iaux( mc ) == 0 ) THEN
               nc = nc + 1
               aux( nc ) = -ng(mc) 
               itmp( nc ) = mc
            END IF
         ENDDO
         CALL hpsort( nc, aux, itmp)
         DO mc = 1, nc
            idx( ic + mc ) = itmp( mc )
         END DO
         DEALLOCATE( itmp )
         DEALLOCATE( aux )
      ELSE
         DO mc = 1, nct
            IF( ng(mc) > 0 .AND. iaux(mc) == 0 ) THEN
               ic = ic + 1
               idx(ic) = mc
            ENDIF
         ENDDO
      ENDIF

      DEALLOCATE( iaux )

      RETURN
   END SUBROUTINE sticks_sort_new

!=----------------------------------------------------------------------=

   SUBROUTINE sticks_dist_new( lgamma, mype, nproc, nyfft, iproc, iproc2, ub, lb, idx, in1, in2, ngc, nct, ncp, ngp, stown, ng )

      IMPLICIT NONE

      LOGICAL, INTENT(in) :: lgamma
      INTEGER, INTENT(in) :: mype
      INTEGER, INTENT(in) :: nproc
      INTEGER, INTENT(in) :: nyfft
      INTEGER, INTENT(in) :: iproc(:,:),iproc2(:)

      INTEGER, INTENT(in) :: ub(:), lb(:), idx(:) ! map's limits
      INTEGER, INTENT(inout) :: stown(lb(1): ub(1), lb(2):ub(2) ) ! stick owner's map, including previously assigned ones
      INTEGER, INTENT(in) :: in1(:), in2(:)   ! cartesian coordinate of each column, ordered according to idx index
      INTEGER, INTENT(in) :: ngc(:)           ! number of G-vector of each column, ordered according to idx index
      INTEGER, INTENT(in) :: nct              ! total number of relevant sticks
      INTEGER, INTENT(out):: ncp(:)           ! number of sticks (columns) per processor
      INTEGER, INTENT(out):: ngp(:)           ! number of G-vectors per processor
      INTEGER, INTENT(out):: ng               ! number of G-vector of this processor. that is ng = ngp(mype+1)

      INTEGER :: mc, i1, i2, j, jj, icnt, gr, j2, j3
      INTEGER, ALLOCATABLE :: yc(:), yg(:) ! number of relevant columns and G-vectors in a given yz plane
      INTEGER, ALLOCATABLE :: ygr(:) ! element in the nyfft group to which a YZ plane belong
      INTEGER, ALLOCATABLE :: ygrp(:), ygrc(:), ygrg(:) ! number of yz planes, relevant columns and G-vectors belonging to a given element in the nyfft group
      !write (6,*) ' inside sticks_map_dist_new'; FLUSH(6)

      ! distribute X values first
      !write (6,*) 'nyfft inside sticks_dist_new', nyfft
      ALLOCATE ( yc(lb(1):ub(1)), yg(lb(1):ub(1)), ygr(lb(1):ub(1)) ) ; yc = 0; yg = 0 ; ygr = 0
      ALLOCATE ( ygrp(nyfft), ygrc(nyfft), ygrg(nyfft) ) ; ygrp = 0; ygrc = 0; ygrg = 0 
      DO mc =1, nct
         if( idx( mc ) < 1 ) CYCLE
         i1 = in1( idx( mc ) )
         i2 = in2( idx( mc ) )
         IF ( ngc( idx(mc) ) > 0 ) THEN
             yc(i1) = yc(i1) + 1
             yg(i1) = yg(i1) + ngc(idx(mc))
         END IF
         IF (stown(i1,i2) > 0) THEN
            gr = iproc2(stown(i1,i2))
            IF (ygr(i1) == 0 ) ygr(i1) = gr
            if (ygr(i1) .ne. gr ) CALL fftx_error__(' sticks_dist ',' ygroups are not compatible ', 1 )
         END IF
      END DO
      do i1 = lb(1),ub(1)
         if (ygr(i1) == 0 ) cycle
         ygrp(ygr(i1)) = ygrp(ygr(i1)) + 1
         ygrc(ygr(i1)) = ygrc(ygr(i1)) + yc(i1)
         ygrg(ygr(i1)) = ygrg(ygr(i1)) + yg(i1)
      end do

      !write (6,*) 'inside sticks_dist_new', lb(1),ub(1)
      !write (6,'(4i4)') (i1,yc(i1),yg(i1),ygr(i1),i1=lb(1),ub(1))
      !write (6,*) 'initial nyfft distribution'
      !write (6,'(4i5)') (i1,ygrp(i1),ygrc(i1),ygrg(i1),i1=1,nyfft)

      ncp = 0
      ngp = 0
      icnt = 0

      DO mc = 1, nct

         IF( idx( mc ) < 1 ) CYCLE

         i1 = in1( idx( mc ) )
         i2 = in2( idx( mc ) )
!
         IF ( lgamma .and. ( (i1 < 0) .or. ( (i1 == 0) .and. (i2 < 0) ) ) ) GOTO 30
!
         IF (ygr(i1) == 0 ) THEN
            j2=1 
            DO j=1,nyfft
               IF ( ygrg(j) < ygrg(j2) ) THEN
                  j2 = j
               ELSEIF ( ( ygrg(j) == ygrg(j2) ) .and. ( ygrc(j) < ygrc(j2) ) ) THEN
                  j2 = j
               END IF
            END DO
            ygr(i1) = j2
            ygrp(j2) = ygrp(j2) + 1
            ygrc(j2) = ygrc(j2) + yc(i1)
            ygrg(j2) = ygrg(j2) + yg(i1)
         ELSE
            j2=ygr(i1)
         END IF

         IF ( ngc( idx(mc) ) > 0 .AND. stown(i1,i2) == 0 ) THEN
            jj = iproc(j2,1)
            DO j3 = 1, nproc / nyfft
               j = iproc(j2,j3)
               IF ( ngp(j) < ngp(jj) ) THEN
                  jj = j
               ELSEIF ( ( ngp(j) == ngp(jj) ) .and. ( ncp(j) < ncp(jj) ) ) THEN
                  jj = j
               ENDIF
            ENDDO
            stown(i1,i2) = jj
         END IF
         IF ( ngc( idx(mc) ) > 0 ) THEN
            ncp( stown(i1,i2) ) = ncp( stown(i1,i2) ) + 1
            ngp( stown(i1,i2) ) = ngp( stown(i1,i2) ) + ngc( idx(mc) )
         ENDIF
 30      CONTINUE
      ENDDO
      !
      ng = ngp( mype + 1 )
      !
      IF ( lgamma ) THEN
         !  when gamma symmetry is used only the sticks of half reciprocal space
         !  are generated, then here we pair-up the sticks with those of the other
         !  half of the space, using the gamma symmetry relation
         !  Note that the total numero of stick "nct" is not modified
         DO mc = 1, nct
            IF( idx( mc ) < 1 ) CYCLE
            IF( ngc( idx(mc) ) < 1 ) CYCLE
            i1 = in1( idx(mc) )
            i2 = in2( idx(mc) )
            IF( i1 == 0 .and. i2 == 0 ) THEN
               jj = stown( i1, i2 )
               IF( jj > 0 ) ngp( jj ) = ngp( jj ) + ngc( idx(mc) ) - 1
            ELSE
               jj = stown( i1, i2 )
               IF( jj > 0 ) THEN
                  stown( -i1, -i2 ) = jj
                  ncp( jj ) = ncp( jj ) + 1
                  ngp( jj ) = ngp( jj ) + ngc( idx(mc) )
               ENDIF
            ENDIF
         ENDDO
      ENDIF

      !write (6,*) 'final nyfft distribution'
      !write (6,'(4i5)') (i1,ygrp(i1),ygrc(i1),ygrg(i1),i1=1,nyfft)
      DEALLOCATE ( ygrp, ygrc, ygrg ) 
      DEALLOCATE ( yc, yg, ygr ) 

      RETURN
   END SUBROUTINE sticks_dist_new

!---------------------------------------------------------------------

   SUBROUTINE get_sticks( smap, gcut, nstp, sstp, st, nst, ng )

      IMPLICIT NONE

      TYPE( sticks_map ), INTENT(INOUT) :: smap
      REAL(DP) , INTENT(in) :: gcut    ! kinetic energy cut-off for this stick distribution

      INTEGER, INTENT(out) :: st(smap%lb(1): smap%ub(1), smap%lb(2):smap%ub(2) )
                                       ! in output it contains (if>0) the processor_id+1 of the owner of a given stick
                                       ! internally it's used to contain the number of G-vectors of a given stick
      INTEGER, INTENT(out) :: nstp(:)  ! number of sticks per processor
      INTEGER, INTENT(out) :: sstp(:)  ! number of G-vectors per processor
      INTEGER, INTENT(out) :: nst      ! total number of sticks
      INTEGER, INTENT(out) :: ng       ! number of G-vector of this processor. that is ng = sstp(mype+1)

      INTEGER, ALLOCATABLE :: ngc(:)
      INTEGER :: ic
      !write (6,*) ' inside get_sticks gcut=',gcut; FLUSH(6)
      !write (6,*) smap%lb, smap%ub

      IF( .NOT. ALLOCATED( smap%stown ) ) THEN
         CALL fftx_error__(' get_sticks ',' sticks map, not allocated ', 1 )
      END IF

      st = 0
      CALL sticks_map_set( smap%lgamma, smap%ub, smap%lb, smap%bg, gcut, st, smap%comm )

      ALLOCATE( ngc ( SIZE( smap%idx ) ) )
      ngc = 0
      CALL sticks_map_index( smap%ub, smap%lb, st, smap%ist(:,1), smap%ist(:,2), ngc, smap%indmap )
      nst = count( st > 0 )
      CALL sticks_sort_new( smap%nproc>1, ngc, SIZE(smap%idx), smap%idx )

      CALL sticks_dist_new( smap%lgamma, smap%mype, smap%nproc, smap%nyfft, smap%iproc, smap%iproc2, &
                            smap%ub, smap%lb, smap%idx, &
                            smap%ist(:,1), smap%ist(:,2), ngc, SIZE(smap%idx), nstp, sstp, smap%stown, ng )

      ! assign the owner of each (relavant) stick 
      st = 0
      DO ic = 1, SIZE( smap%idx )
         IF( smap%idx( ic ) > 0 ) THEN
            IF( ngc( smap%idx( ic ) ) > 0 ) THEN
                st( smap%ist(smap%idx( ic ),1), smap%ist(smap%idx( ic ),2) ) = &
                    smap%stown( smap%ist(smap%idx( ic ),1),smap%ist(smap%idx( ic ),2))
                IF(smap%lgamma) st(-smap%ist(smap%idx( ic),1),-smap%ist(smap%idx( ic ),2)) = &
                   smap%stown( smap%ist(smap%idx( ic ),1),smap%ist(smap%idx( ic ),2))
            END IF
         END IF
      END DO
      DEALLOCATE( ngc )
      RETURN
   END SUBROUTINE get_sticks

!---------------------------------------------------------------------
   SUBROUTINE hpsort (n, ra, ind)
      !---------------------------------------------------------------------
      ! sort an array ra(1:n) into ascending order using heapsort algorithm.
      ! n is input, ra is replaced on output by its sorted rearrangement.
      ! create an index table (ind) by making an exchange in the index array
      ! whenever an exchange is made on the sorted data array (ra).
      ! in case of equal values in the data array (ra) the values in the
      ! index array (ind) are used to order the entries.
      ! if on input ind(1)  = 0 then indices are initialized in the routine,
      ! if on input ind(1) != 0 then indices are assumed to have been
      !                initialized before entering the routine and these
      !                indices are carried around during the sorting process
      !
      ! no work space needed !
      ! free us from machine-dependent sorting-routines !
      !
      ! adapted from Numerical Recipes pg. 329 (new edition)
      !
      IMPLICIT NONE
      !-input/output variables
      INTEGER :: n
      INTEGER :: ind (n)
      REAL(DP) :: ra (n)
      !-local variables
      INTEGER :: i, ir, j, l, iind
      REAL(DP) :: rra
      !
      IF (n < 1 ) RETURN
      ! initialize index array
      IF (ind (1) ==0) THEN
         DO i = 1, n
            ind (i) = i
         ENDDO
      ENDIF
      ! nothing to order
      IF (n < 2) RETURN
      ! initialize indices for hiring and retirement-promotion phase
      l = n / 2 + 1
      ir = n
10    CONTINUE
      ! still in hiring phase
      IF (l>1) THEN
         l = l - 1
         rra = ra (l)
         iind = ind (l)
         ! in retirement-promotion phase.
      ELSE
         ! clear a space at the end of the array
         rra = ra (ir)
         !
         iind = ind (ir)
         ! retire the top of the heap into it
         ra (ir) = ra (1)
         !
         ind (ir) = ind (1)
         ! decrease the size of the corporation
         ir = ir - 1
         ! done with the last promotion
         IF (ir==1) THEN
            ! the least competent worker at all !
            ra (1) = rra
            !
            ind (1) = iind
            RETURN
         ENDIF
      ENDIF
      ! wheter in hiring or promotion phase, we
      i = l
      ! set up to place rra in its proper level
      j = l + l
      !
      DO WHILE (j<=ir)
         IF (j<ir) THEN
            ! compare to better underling
            IF (ra (j) <ra (j + 1) ) THEN
               j = j + 1
            ELSEIF (ra (j) ==ra (j + 1) ) THEN
               IF (ind (j) <ind (j + 1) ) j = j + 1
            ENDIF
         ENDIF
         ! demote rra
         IF (rra<ra (j) ) THEN
            ra (i) = ra (j)
            ind (i) = ind (j)
            i = j
            j = j + j
         ELSEIF (rra==ra (j) ) THEN
            ! demote rra
            IF (iind<ind (j) ) THEN
               ra (i) = ra (j)
               ind (i) = ind (j)
               i = j
               j = j + j
            ELSE
               ! set j to terminate do-while loop
               j = ir + 1
            ENDIF
            ! this is the right place for rra
         ELSE
            ! set j to terminate do-while loop
            j = ir + 1
         ENDIF
      ENDDO
      ra (i) = rra
      ind (i) = iind
      GOTO 10
      !
   END SUBROUTINE hpsort

!=----------------------------------------------------------------------=
END MODULE stick_base
!=----------------------------------------------------------------------=
