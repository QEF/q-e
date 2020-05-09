SUBROUTINE exx_index_pair(wannierc, overlap_final, nj_final, nj_max, ndim )
    !=======================================================================================
    ! Code Version 1.0 (Princeton University, September 2014)
    !=======================================================================================
    !==========================================================================
    !
    !This subroutine computes unique pair list for each Wannier orbital
    !
    !==========================================================================
    USE kinds,                   ONLY  : DP
    USE electrons_base,          ONLY  : nbsp, nspin, nupdwn, iupdwn
    USE mp_global,               ONLY  : me_image, nproc_image, intra_image_comm
    USE cell_base,               ONLY  : h,ainv 
    USE parallel_include
    USE mp,                      ONLY  : mp_barrier 
    USE wannier_base,            ONLY  : neigh, dis_cutoff
    USE io_global,               ONLY  : stdout
    USE exx_module,              ONLY  : prev_obtl_overlap
    USE wrappers                
    !
    IMPLICIT NONE
    !
    REAl(DP),INTENT(IN)    :: wannierc(3, nbsp)
    INTEGER, INTENT(INOUT) :: nj_max, ndim
    INTEGER, INTENT(INOUT) :: overlap_final(neigh/2,ndim), nj_final(ndim) 
    !
    REAl(DP), ALLOCATABLE  :: distance(:)
    INTEGER,  ALLOCATABLE  :: overlap(:,:),overlap2(:),overlap3(:,:),overlap4(:,:),nj(:),my_nj(:)
    INTEGER,  ALLOCATABLE  :: njl(:),njl_tmp(:),njl_position(:),njl_index(:)
    INTEGER                :: i, j, k, jj, ip, ir, ierr, num, iobtl,my_nj_proc,sc_fac,least_est_obtl,pos,temp,target_obtl
    INTEGER                :: least_k_obtl
    INTEGER                :: prev_recv_factor
    REAl(DP)               :: dij(3),dij2(3),least_est_job,my_est_job,least_factor
    !
    INTEGER   gindx_of_i, jj_neib_of_i
    INTEGER   spin_loop, nupdwn_(2),iupdwn_(2)
    INTEGER             :: munit = 1992
    CHARACTER (len=300) :: print_str
    !------------------------------------------------------------------
    !
    !
    ALLOCATE( distance(nbsp) ); distance=0.0_DP
    ALLOCATE( overlap(neigh, nbsp)); overlap=0
    ALLOCATE( overlap2(neigh) ); overlap2=0
    ALLOCATE( nj(nbsp)); nj=0
    !
    IF (nspin == 1) THEN
      nupdwn_(1)=nbsp
      nupdwn_(2)=nbsp
      iupdwn_(1)=1
      iupdwn_(2)=1
      spin_loop=1
    ELSE
      nupdwn_(:)=nupdwn(:)
      iupdwn_(:)=iupdwn(:)
    END IF
    !
    DO i = 1, nbsp
      !
      nj(i) = 0
      distance(:) = 0.0_DP
      !
      IF ( i < iupdwn(2) ) THEN
        spin_loop = 1
      ELSE
        spin_loop = 2
      END IF
      !
      DO j = iupdwn_(spin_loop), iupdwn_(spin_loop) + nupdwn_(spin_loop)-1
        !
        IF ( j .NE. i ) THEN
          !
          ! Compute distance between wfc i and wfc j (according to the minimum image convention)...
          !
          dij(1)=wannierc(1,i)-wannierc(1,j)   ! r_ij = r_i - r_j   
          dij(2)=wannierc(2,i)-wannierc(2,j)   ! r_ij = r_i - r_j   
          dij(3)=wannierc(3,i)-wannierc(3,j)   ! r_ij = r_i - r_j   
          !
          dij2(1)=ainv(1,1)*dij(1)+ainv(1,2)*dij(2)+ainv(1,3)*dij(3)   ! s_ij = h^-1 r_ij
          dij2(2)=ainv(2,1)*dij(1)+ainv(2,2)*dij(2)+ainv(2,3)*dij(3)   ! s_ij = h^-1 r_ij
          dij2(3)=ainv(3,1)*dij(1)+ainv(3,2)*dij(2)+ainv(3,3)*dij(3)   ! s_ij = h^-1 r_ij
          !
          dij2(1)=dij2(1)-IDNINT(dij2(1))   ! impose MIC on s_ij in range: [-0.5,+0.5]
          dij2(2)=dij2(2)-IDNINT(dij2(2))   ! impose MIC on s_ij in range: [-0.5,+0.5]
          dij2(3)=dij2(3)-IDNINT(dij2(3))   ! impose MIC on s_ij in range: [-0.5,+0.5]
          !
          dij(1)=h(1,1)*dij2(1)+h(1,2)*dij2(2)+h(1,3)*dij2(3)   ! r_ij = h s_ij (MIC)
          dij(2)=h(2,1)*dij2(1)+h(2,2)*dij2(2)+h(2,3)*dij2(3)   ! r_ij = h s_ij (MIC)
          dij(3)=h(3,1)*dij2(1)+h(3,2)*dij2(2)+h(3,3)*dij2(3)   ! r_ij = h s_ij (MIC)
          !
          distance(j)=DSQRT(dij(1)*dij(1)+dij(2)*dij(2)+dij(3)*dij(3))   ! |r_i - r_j| (MIC)
          !
          IF ( distance(j) .LT.  dis_cutoff ) THEN
            !
            nj(i) = nj(i) + 1
            !
            IF (nj(i) .EQ. 1 ) THEN
              !
              overlap(nj(i),i) = j
              !
            ELSE IF (distance(overlap(nj(i)-1,i)) .LE. distance(j)) THEN
              !
              overlap(nj(i),i)=j
              !
            ELSE
              !
              overlap2(:)=0
              DO ir=1,nj(i)-1
                !
                IF (distance(overlap(ir,i)) < distance(j)) THEN
                  overlap2(ir)=overlap(ir,i)
                ELSE
                  overlap2(ir)=j
                  DO ip=ir+1,nj(i)
                    overlap2(ip)=overlap(ip-1,i)
                  END DO
                  GO TO 555
                END IF
                !
              END DO
              !
              555      CONTINUE
              !
              DO ir = 1, nj(i)
                overlap(ir,i)=overlap2(ir)
              END DO
              !
            END IF
            !
          END IF !if for distance(j)<5.0d0
          !
        END IF ! j=/i
        !
      END DO  ! j=1,nbsp
      !
      IF (nj(i) > neigh) THEN
        WRITE(print_str,'(3X,"EXX calculation error : maximum non-unique pairs found",1X,I5,&
            & "increase the value of exx_neigh at least to",1X,I5)') nj(i),nj(i) 
        !
        CALL errore('exx_pair',print_str,1)
        !
      END IF
      !
    END DO !i = 1, nbsp
    !
    IF (ALLOCATED(distance))         DEALLOCATE(distance)
    !------------------------------------------------------------------------------------------


    !------------------------------------------------------------------------------------------
    ! Now form unique wfc pair overlap matrix...
    !------------------------------------------------------------------------------------------
    ALLOCATE( overlap3(neigh/2, nbsp)); overlap3=0
    ALLOCATE( njl(nbsp) );              njl=0
    ALLOCATE( njl_tmp(nbsp) );          njl_tmp=0
    ALLOCATE( njl_index(nbsp) );        njl_index=0
    ALLOCATE( njl_position(nbsp) );     njl_position=0
    !------------------------------------------------------------------------------------------
    IF ( .NOT.ALLOCATED(prev_obtl_overlap) ) THEN
      CALL unique_pairlist_scratch(neigh, nbsp, overlap, overlap3, nj, njl, njl_tmp, njl_index, njl_position)
    ELSE
      !debug
      !--------------------------------------------
      write(print_str,fmt='(A,I0.4,A)') "overlap", me_image, ".dat"
      !--------------------------------------------
!     open(unit=munit,file=print_str)
!     !--------------------------------------------
!     write(munit,fmt='(I4)') neigh
!     !--------------------------------------------
!     write(munit,fmt='(I4)') nbsp
!     !--------------------------------------------
!     do j=1,nbsp
!       do i=1,neigh/2
!         write(munit, fmt='(I4)',advance='no') prev_obtl_overlap(i,j)
!       end do
!       write(munit,fmt='(A)') ''
!     end do
!     !--------------------------------------------
!     do j=1,nbsp
!       do i=1,neigh
!         write(munit,fmt='(I4)',advance='no') overlap(i,j)
!       end do
!       write(munit,fmt='(A)') ''
!     end do
!     !--------------------------------------------
!     do j=1,nbsp
!       do i=1,neigh/2
!         write(munit, fmt='(I4)',advance='no') overlap3(i,j)
!       end do
!       write(munit,fmt='(A)') ''
!     end do
!     !--------------------------------------------
!     close(munit)
      !--------------------------------------------
      !debug
      write(print_str,fmt='(A,I0.4,A)') "cout", me_image, ".dat"
      CALL load_balancing(neigh, nbsp, prev_obtl_overlap, overlap, overlap3, TRIM(print_str))
    END IF
    !------------------------------------------------------------------------------------------
    IF ( .NOT.ALLOCATED(prev_obtl_overlap) ) THEN
      ALLOCATE( prev_obtl_overlap(neigh/2, nbsp) ); prev_obtl_overlap=0
    END IF
    !------------------------------------------------------------------------------------------
    prev_obtl_overlap = overlap3
    !------------------------------------------------------------------------------------------
    DO i = 1, nbsp
      num = 0
      DO j = 1, neigh/2
        IF (overlap3(j,i) .NE. 0) num = num + 1
      END DO
      nj(i) = num
    END DO
    !------------------------------------------------------------------------------------------


    !---------------------------------------------------------------------------------------------------------------------------
    !  for using more processors than number of bands a new overlap matrix, overlap4, is made
    !  overlap4 contains indices of distributed orbital pairs when number of processors are larger than number of bands (nbsp)  
    !  an example, if orbital 1 is paired with 15 orbitals with following indices 
    !  12   72   94   17  149   87  183  220   11  180   83  223  115  154   92
    !  then the pairs are redistributed, 8 orbitals in proc 1 and 7 orbitals in proc (1+nbsp) 
    !  proc 1      :       12   72   94   17  149   87  183  220 
    !  proc 1+nbsp :       11  180   83  223  115  154   92
    !
    !  likewise, this can be further distributed among larger number of processors
    !  at present nproc_image has to be integer multiple of nbsp 
    !---------------------------------------------------------------------------------------------------------------------------
    IF (nproc_image.GT.nbsp) THEN
      !
      ALLOCATE(my_nj(nproc_image)); my_nj=0
      ALLOCATE(overlap4(neigh/2,nproc_image)); overlap4=0
      !
      sc_fac = nproc_image/nbsp
      !
      DO i = 1, nbsp
        !
        DO j = 1, neigh/2
          !
          IF ( overlap3(j, i) .NE. 0 ) THEN
            !
            DO k = 1, sc_fac
              !
              IF ( prev_recv_factor( (k-1)*nbsp+i, overlap3(j, i), neigh ) .EQ. 1 ) THEN
                my_nj( (k-1)*nbsp+i )                           = my_nj( (k-1)*nbsp+i ) + 1
                overlap4( my_nj( (k-1)*nbsp+i ), (k-1)*nbsp+i ) = overlap3(j, i)
                overlap3(j, i)                                  = 0
                EXIT
              END IF
              !
            END DO
            !
          END IF
          !
        END DO
        !
        DO j = 1, neigh/2
          !
          IF ( overlap3(j, i) .NE. 0 ) THEN
            !
            least_factor = 10000000.0D0
            least_k_obtl = 0
            !
            DO k = 1, sc_fac
              !
              IF (my_nj( (k-1)*nbsp+i ) + DBLE(k-1)/sc_fac < least_factor) THEN
                least_factor = my_nj( (k-1)*nbsp+i ) + DBLE(k-1)/sc_fac
                least_k_obtl = (k-1)*nbsp+i
              END IF
              !
            END DO
            !
            my_nj(least_k_obtl) = my_nj(least_k_obtl) + 1
            overlap4(my_nj(least_k_obtl), least_k_obtl) = overlap3(j, i)
            overlap3(j, i) = 0
            !
          END IF
          !
        END DO
        !
      END DO
      !
      ! Update final overlap matrix and neighbor list 
      !
      overlap_final=overlap4
      nj_final=my_nj
      nj_max=MAXVAL(my_nj)
      !
    ELSE
      !
      ! Update final overlap matrix and neighbor list 
      !
      overlap_final=overlap3
      nj_final=nj
      nj_max=MAXVAL(nj)
      !
    END IF
    !------------------------------------------------------------------------------------------



    !------------------------------------------------------------------
    IF (ALLOCATED(overlap))          DEALLOCATE(overlap)
    IF (ALLOCATED(overlap2))         DEALLOCATE(overlap2)
    IF (ALLOCATED(overlap3))         DEALLOCATE(overlap3)
    IF (ALLOCATED(overlap4))         DEALLOCATE(overlap4)
    IF (ALLOCATED(nj))               DEALLOCATE(nj)
    IF (ALLOCATED(my_nj))            DEALLOCATE(my_nj)
    IF (ALLOCATED(njl))              DEALLOCATE(njl)
    IF (ALLOCATED(njl_tmp))          DEALLOCATE(njl_tmp)
    IF (ALLOCATED(njl_index))        DEALLOCATE(njl_index)
    IF (ALLOCATED(njl_position))     DEALLOCATE(njl_position)
    !------------------------------------------------------------------
    
    !------------------------------------------------------------------------------------------
    WRITE(stdout,'(3X,"nj_max nj_min nj_avg :",2I6,F6.2)')nj_max,MINVAL(nj_final),DBLE(SUM(nj_final))/DBLE(ndim)
    !------------------------------------------------------------------------------------------

    RETURN
END SUBROUTINE exx_index_pair

!==========================================================================================


!==========================================================================================
SUBROUTINE unique_pairlist_scratch(neigh, nbsp, overlap, overlap3, nj, njl, njl_tmp, njl_index, njl_position)
    !------------------------------------------------------------------
    USE kinds,                 ONLY  : DP
    !------------------------------------------------------------------
    IMPLICIT NONE
    !------------------------------------------------------------------
    ! Pass in variables
    !------------------------------------------------------------------
    INTEGER                 :: neigh, nbsp
    INTEGER                 :: overlap(neigh, nbsp), overlap3(neigh/2, nbsp)
    INTEGER                 :: nj(nbsp), njl(nbsp), njl_tmp(nbsp), njl_index(nbsp), njl_position(nbsp)
    !------------------------------------------------------------------

    !------------------------------------------------------------------
    ! Local variables
    !------------------------------------------------------------------
    INTEGER                 :: i, j
    INTEGER                 :: num, least_est_obtl, pos, target_obtl, temp
    REAL(DP)                :: least_est_job, my_est_job
    !------------------------------------------------------------------


    !------------------------------------------------------------------
    ! Functions
    !------------------------------------------------------------------
    INTEGER                :: prev_factor
    !------------------------------------------------------------------


    !------------------------------------------------------------------
    nj=0
    !
    DO i = 1,nbsp
      num = 0
      DO j = 1,neigh
        IF (overlap(j,i) .NE. 0) num = num+1
      END DO
      njl(i) = num
    END DO
    !
    njl_tmp = njl
    !
    CALL quick_sort(nbsp, njl_tmp, njl_index)
    !
    DO i = 1,nbsp
      njl_position( njl_index(i) ) = i
    END DO
    !
    DO
      IF (njl(njl_index(nbsp)) .EQ. 0) THEN
        EXIT
      ELSE                        
        !------------------------------------------------------------------
        ! Pick where to do the job
        !------------------------------------------------------------------
        least_est_obtl = 100000000
        least_est_job  = 100000000.0D0
        pos            = 0
        !
        DO i = 1, njl(njl_index(nbsp))
          !
          target_obtl = overlap(i,njl_index(nbsp))
          !
          !-----------------------------------------------------------------------
          ! The estimate function is $(number_of_jobs) + 0.5*$(num_of_jobs_left)
          !-----------------------------------------------------------------------
          my_est_job = nj(target_obtl) + 0.5*njl(target_obtl) + (-0.25) * prev_factor(target_obtl, njl_index(nbsp), neigh)
          !
          IF (my_est_job .LT. least_est_job) THEN
            !
            pos = i
            least_est_obtl = target_obtl
            least_est_job  = my_est_job
            !
          END IF
          !
        END DO
        !
        nj(least_est_obtl)                           = nj(least_est_obtl) + 1
        overlap3(nj(least_est_obtl), least_est_obtl) = njl_index(nbsp)
        !
        njl(least_est_obtl)   = njl(least_est_obtl)  - 1
        njl(njl_index(nbsp))  = njl(njl_index(nbsp)) - 1
        !
        DO i = pos, njl(njl_index(nbsp))
          overlap(i, njl_index(nbsp)) = overlap(i+1, njl_index(nbsp))
        END DO
        overlap(njl(njl_index(nbsp))+1, njl_index(nbsp)) = 0
        !
        DO i = 1, neigh
          IF ( overlap(i, least_est_obtl) .EQ. njl_index(nbsp) ) THEN
            pos = i
            EXIT
          END IF
        END DO
        DO i = pos, neigh-1
          overlap(i, least_est_obtl) = overlap(i+1, least_est_obtl)
        END DO
        overlap(neigh, least_est_obtl) = 0
        !
        !=============================================================
        !
        IF ( njl(least_est_obtl) .EQ. 0 ) THEN
          !
          DO j = njl_position(least_est_obtl), 2, -1
              !
              temp                         = njl_position(njl_index(j))
              njl_position(njl_index(j))   = njl_position(njl_index(j-1))
              njl_position(njl_index(j-1)) = temp
              !
              temp                         = njl_index(j)
              njl_index(j)                 = njl_index(j-1)
              njl_index(j-1)               = temp
              !
          END DO
          !
        ELSE
          !
          DO j = njl_position(least_est_obtl), nbsp - 2
            IF (nj(njl_index(j)) + 0.5 * njl(njl_index(j)) > nj(njl_index(j+1)) + 0.5 * njl(njl_index(j+1))) THEN
              !
              temp                         = njl_position(njl_index(j))
              njl_position(njl_index(j))   = njl_position(njl_index(j+1))
              njl_position(njl_index(j+1)) = temp
              !
              temp                = njl_index(j)
              njl_index(j)        = njl_index(j+1)
              njl_index(j+1)      = temp
              !
            ELSE
              EXIT
            END IF
          END DO
          !
        END IF
        !
        !=============================================================
        !
        IF ( njl(njl_index(nbsp)) .EQ. 0 ) THEN
          !
          DO j = nbsp, 2, -1
              !
              temp                = njl_position(njl_index(j))
              njl_position(njl_index(j)) = njl_position(njl_index(j-1))
              njl_position(njl_index(j-1)) = temp
              !
              temp                = njl_index(j)
              njl_index(j)        = njl_index(j-1)
              njl_index(j-1)      = temp
              !
          END DO
          !
        ELSE
          !
          DO j = nbsp, 2, -1
            IF((nj(njl_index(j))+0.5*njl(njl_index(j))<=nj(njl_index(j-1))+0.5*njl(njl_index(j-1))).AND.&
               (njl(njl_index(j-1)).NE.0))THEN
              !
              temp                = njl_position(njl_index(j))
              njl_position(njl_index(j)) = njl_position(njl_index(j-1))
              njl_position(njl_index(j-1)) = temp
              !
              temp                = njl_index(j)
              njl_index(j)        = njl_index(j-1)
              njl_index(j-1)      = temp
              !
            ELSE
              EXIT
            END IF
          END DO
          !
        END IF
        !
        !=============================================================
        !
      END IF
      !
    END DO
    !
    !==========================================================================================
    !
    DO i = 1,nbsp
      num = 0
      DO j = 1,neigh/2
        IF (overlap3(j,i) .NE. 0) num = num+1
      END DO
      njl(i) = num
    END DO
    !
    !==========================================================================================
END SUBROUTINE unique_pairlist_scratch
!==========================================================================================


!==========================================================================================
SUBROUTINE exx_index_pair_nv(wc, overlap, nj, nj_max)
    !====================================================================
    !  This subroutine finds for each state the valence-state neighbors 
    !  Lingzhu Kong
    !====================================================================
    USE kinds,                 ONLY  : DP
    USE electrons_base,        ONLY  : nbsp
    USE exx_module,            ONLY  : vwc
    USE mp_global,             ONLY  : me_image, intra_image_comm, nproc_image
    USE cell_base,             ONLY  : h 
    USE parallel_include
    USE wannier_base,          ONLY  : neigh, dis_cutoff, vnbsp
    !
    IMPLICIT NONE
    !
    ! wc is the wannier centers of the initial quasi-particle states
    REAl(DP),INTENT(IN)    ::    wc(3, nbsp)
    ! number and index of overlapping orbitals for all conduction states
    INTEGER, INTENT(INOUT) ::    overlap(neigh,nbsp), nj(nbsp) , nj_max
    !
    INTEGER     i, j, ierr
    REAl(DP)    centerx, centery, centerz, xi, yi, zi, xj, yj, zj, xij, yij, zij, distance
    !
    !==================================================================
    !
    PRINT *, 'entering exx_index_pair_nv', dis_cutoff, neigh, vnbsp
    centerx = 0.5 * h(1,1)
    centery = 0.5 * h(2,2)
    centerz = 0.5 * h(3,3)
    !
    overlap(:,:) = 0
    !
    DO i = 1, nbsp
      !
      nj(i) = 0
      xi = wc( 1, i )
      yi = wc( 2, i )
      zi = wc( 3, i )
      !
      DO j = 1, vnbsp
        xj = vwc(1,j)
        yj = vwc(2,j)
        zj = vwc(3,j)
        xij = xj - xi - INT( (xj-xi)/centerx )*h(1,1) 
        yij = yj - yi - INT( (yj-yi)/centery )*h(2,2)  
        zij = zj - zi - INT( (zj-zi)/centerz )*h(3,3) 
        distance = sqrt( xij*xij + yij*yij + zij*zij)
        !
        IF ( distance .LT. dis_cutoff ) THEN
          nj(i) = nj(i) + 1
          IF (nj(i) > neigh) THEN
            PRINT *, 'increase neigh, stop in exx_pair', nj(i), neigh
            RETURN
          END IF
          overlap(nj(i),i) = j
        END IF
      END DO  ! j=1,vnbsp
    END DO  
    !
    nj_max = nj(1)
    DO i = 2, nbsp
      IF(nj(i) > nj_max) nj_max = nj(i)
    ENDDO
    !
    PRINT *, 'leave exx_index_pair_nv', nj
    RETURN
    !
END SUBROUTINE exx_index_pair_nv
!==================================================================

!==================================================================
INTEGER FUNCTION prev_factor(target_obtl, obtl_tbc, neigh)
    !
    USE exx_module,     ONLY : prev_obtl_overlap
    !
    INTEGER target_obtl, obtl_tbc, itr
    !
    prev_factor = 0
    IF (ALLOCATED(prev_obtl_overlap)) THEN
      !
      DO itr = 1, neigh/2 
        !
        IF (prev_obtl_overlap(itr, target_obtl) .EQ. obtl_tbc) THEN
          prev_factor = 1
          EXIT
        END IF
        !
      END DO
      !
    END IF
    !
    RETURN
END
!==================================================================

!==================================================================
INTEGER FUNCTION prev_recv_factor(target_obtl, obtl_tbc, neigh)
    !
    USE exx_module,     ONLY : prev_obtl_recv
    !
    INTEGER target_obtl, obtl_tbc, itr
    !
    prev_recv_factor = 0
    IF (ALLOCATED(prev_obtl_recv)) THEN
      !
      DO itr = 1, neigh/2 
        !
        IF (prev_obtl_recv(itr, target_obtl) .EQ. obtl_tbc) THEN
          prev_recv_factor = 1
          EXIT
        END IF
        !
      END DO
      !
    END IF
    !
    RETURN
END
!==================================================================

!==================================================================
RECURSIVE SUBROUTINE quick_sort(array_size, list, order)

    ! Quick sort routine from:
    ! Brainerd, W.S., Goldberg, C.H. & Adams, J.C. (1990) "Programmer's Guide to
    ! Fortran 90", McGraw-Hill  ISBN 0-07-000248-7, pages 149-150.
    ! Modified by Alan Miller to include an associated integer array which gives
    ! the positions of the elements in the original order.

    IMPLICIT NONE
    INTEGER                                 :: array_size
    INTEGER, DIMENSION (:), INTENT(IN OUT)  :: list(array_size)
    INTEGER, DIMENSION (:), INTENT(OUT)     :: order(array_size)

    ! Local variable
    INTEGER :: i

    DO i = 1, SIZE(list)
      order(i) = i
    END DO

    CALL quick_sort_1(1, SIZE(list))

    CONTAINS

    RECURSIVE SUBROUTINE quick_sort_1(left_end, right_end)

    INTEGER, INTENT(IN) :: left_end, right_end

    !     Local variables
    INTEGER             :: i, j, itemp
    INTEGER             :: reference, temp
    INTEGER, PARAMETER  :: max_simple_sort_size = 6

    IF (right_end < left_end + max_simple_sort_size) THEN
      ! Use interchange sort for small lists
      CALL interchange_sort(left_end, right_end)
    
    ELSE
      ! Use partition ("quick") sort
      reference = list((left_end + right_end)/2)
      i = left_end - 1; j = right_end + 1
    
      DO
        ! Scan list from left end until element >= reference is found
        DO
          i = i + 1
          IF (list(i) >= reference) EXIT
        END DO
        ! Scan list from right end until element <= reference is found
        DO
          j = j - 1
          IF (list(j) <= reference) EXIT
        END DO
    
    
        IF (i < j) THEN
          ! Swap two out-of-order elements
          temp = list(i); list(i) = list(j); list(j) = temp
          itemp = order(i); order(i) = order(j); order(j) = itemp
        ELSE IF (i == j) THEN
          i = i + 1
          EXIT
        ELSE
          EXIT
        END IF
      END DO
    
      IF (left_end < j) CALL quick_sort_1(left_end, j)
      IF (i < right_end) CALL quick_sort_1(i, right_end)
    END IF
    
    END SUBROUTINE quick_sort_1
    
    
    SUBROUTINE interchange_sort(left_end, right_end)
    
    INTEGER, INTENT(IN) :: left_end, right_end
    
    !     Local variablesl
    INTEGER             :: i, j, itemp
    INTEGER             :: temp
    
    DO i = left_end, right_end - 1
      DO j = i+1, right_end
        IF (list(i) > list(j)) THEN
          temp = list(i); list(i) = list(j); list(j) = temp
          itemp = order(i); order(i) = order(j); order(j) = itemp
        END IF
      END DO
    END DO
    
    END SUBROUTINE interchange_sort
    
END SUBROUTINE quick_sort
!==================================================================
