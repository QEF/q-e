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
    !
    IMPLICIT NONE
    !
    REAl(DP),INTENT(IN)    :: wannierc(3, nbsp)
    INTEGER, INTENT(INOUT) :: nj_max, ndim
    INTEGER, INTENT(INOUT) :: overlap_final(neigh/2,ndim), nj_final(ndim) 
    !
    REAl(DP), ALLOCATABLE  :: distance(:)
    INTEGER,  ALLOCATABLE  :: overlap(:, :),overlap2(:),overlap3(:,:),overlap4(:,:),nj(:),my_nj(:)
    INTEGER                :: i, j, k,jj, ip, ir, ierr, num, num1, iobtl,my_nj_proc,sc_fac
    REAl(DP)               :: dij(3),dij2(3)
    !
    INTEGER   gindx_of_i, jj_neib_of_i
    INTEGER   spin_loop, nupdwn_(2),iupdwn_(2)
    CHARACTER (len=300) :: print_str
    !------------------------------------------------------------------
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
    !
    ! Now form unique wfc pair overlap matrix...
    !
    ALLOCATE( overlap3(neigh/2, nbsp)); overlap3=0
    !
    num=0; num1=0
    !
    DO j=1, neigh/2
      !
      DO i=1, nbsp
        !
        DO jj=1, nj(i)
          !
          jj_neib_of_i = overlap(jj,i)
          !
          IF (jj_neib_of_i .GT. 0) THEN
            !
            overlap3(j, i) = jj_neib_of_i
            overlap(jj,i) = 0
            num = num + 1
            DO k = 1, nj(jj_neib_of_i)
              IF (overlap(k,jj_neib_of_i) .EQ. i) THEN
                overlap(k,jj_neib_of_i) = 0
                num1 = num1 + 1
                GO TO 666
              END IF
            END DO
            !
          END IF
          !
        END DO !jj
        !
        666 CONTINUE
        !
      END DO !i
      !
    END DO  !j
    !
    IF (num .NE. num1) CALL errore("exx", "sth wrong with overlap",1)
    !
    DO i = 1,nbsp
      num = 0
      DO j = 1,neigh/2
        IF (overlap3(j,i) .NE. 0) num = num+1
      END DO
      nj(i) = num
    END DO
    !
    nj_max=MAXVAL(nj)
    !
    IF (ALLOCATED(overlap))          DEALLOCATE(overlap)
    IF (ALLOCATED(overlap2))         DEALLOCATE(overlap2)
    !
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
    !
    IF (nproc_image.GT.nbsp) THEN
      !
      ALLOCATE(my_nj(nproc_image)); my_nj=0
      ALLOCATE(overlap4(nj_max,nproc_image)); overlap4=0
      !
      sc_fac = nproc_image/nbsp
      !
      DO i = 1, nbsp
        !
        IF( MOD(nj(i), sc_fac) .NE. 0) THEN
          my_nj_proc = nj(i)/sc_fac + 1
        ELSE
          my_nj_proc = nj(i)/sc_fac
        END IF
        !
        jj=0; iobtl=i
        !
        DO j = 1, nj(i)
          jj = jj+1
          IF(jj .GT. my_nj_proc)THEN
            jj=1; iobtl=iobtl+nbsp
          END IF
          overlap4(jj,iobtl)=overlap3(j,i)
          my_nj(iobtl)=my_nj(iobtl)+1
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
    !
    IF (ALLOCATED(overlap3))         DEALLOCATE(overlap3)
    IF (ALLOCATED(overlap4))         DEALLOCATE(overlap4)
    IF (ALLOCATED(nj))               DEALLOCATE(nj)
    IF (ALLOCATED(my_nj))            DEALLOCATE(my_nj)
    !
    WRITE(stdout,'(3X,"nj_max nj_min nj_avg :",2I6,F6.2)')nj_max,MINVAL(nj_final),DBLE(SUM(nj_final))/DBLE(ndim)
    !
    !IF (me_image .eq. 0) THEN
    !  open(unit=20,file='pair.dat',status='unknown',form='formatted')
    !  do i = 1, ndim
    !     write(20, '(I5," :",I5," :",30I5)')i,nj_final(i),(overlap_final(j,i),j=1, nj_final(i))
    !  enddo
    !  close(20)
    !END IF
    !
    !DEBUG: print matrix overlap_final
    !DO i = 1, ndim 
    !  WRITE(stdout,'(I5," :",I5,":",16I5)')i,nj_final(i),(overlap_final(k,i),k=1,nj_final(i))
    !END DO
    !
    RETURN
    !
END SUBROUTINE exx_index_pair

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
