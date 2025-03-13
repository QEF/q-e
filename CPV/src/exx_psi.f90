SUBROUTINE exx_psi(c, psitot2,nnrtot,my_nbsp, my_nxyz, nbsp)
    !=======================================================================================
    ! Code Version 1.0 (Princeton University, September 2014)
    !=======================================================================================
    !====================================================================
    !   get the wave function in real space defined by the FFT grids
    !====================================================================
    !
    USE kinds,                   ONLY  : DP
    USE fft_interfaces,          ONLY  : invfft
    USE fft_base,                ONLY  : dffts, dfftp
    USE gvecw,                   ONLY  : ngw
    USE mp_global,               ONLY  : nproc_image, me_image,intra_image_comm
    USE cell_base,               ONLY  : omega
    USE parallel_include
    USE mp,                      ONLY  : mp_barrier
    USE mp_wave,                 ONLY  : redistwfr
    USE io_global,               ONLY  : stdout         !print/write argument for standard output (to output file)
    USE fft_helper_subroutines
    !
    IMPLICIT NONE
    !
    INTEGER,     INTENT(IN)    :: nnrtot, nbsp
    INTEGER,     INTENT(IN)    :: my_nbsp(nproc_image), my_nxyz(nproc_image)
    COMPLEX(DP), INTENT(INOUT) :: c(ngw,nbsp)
    REAL(DP),    INTENT(INOUT) :: psitot2(nnrtot, my_nbsp(me_image+1) )
    !
    INTEGER i, j, ir,proc, me, ierr, ig, irank, iobtl
    !
    COMPLEX(DP), ALLOCATABLE   ::    psis(:), psis1(:)
    REAL(DP),    ALLOCATABLE   ::    psis2(:,:), psitot(:)
    COMPLEX(DP)  ca(ngw)       
    complex(DP), parameter     ::    ci=(0.0d0,1.0d0)
    !
    INTEGER      sizefft       
    !
    INTEGER,    ALLOCATABLE    :: my_nnr(:)
    INTEGER                    :: sc_fac,ii,jj,eig_offset,eig_index,ngpww1,va,nogrp,nr1s,nr2s,nr3s,nnr2
    INTEGER,    ALLOCATABLE    :: sdispls(:), sendcount(:)
    INTEGER,    ALLOCATABLE    :: rdispls(:), recvcount(:)
    INTEGER,    ALLOCATABLE    :: sdispls1(:), sendcount1(:)
    INTEGER,    ALLOCATABLE    :: rdispls1(:), recvcount1(:)
    INTEGER                    :: igoff, idx
#ifdef __CUDA
    COMPLEX(DP), ALLOCATABLE, DEVICE   ::    psis_d(:), ca_d(:), c_d(:,:)
    REAL(DP),    ALLOCATABLE, DEVICE   ::    psis2_d(:,:)
#endif
    !
    me=me_image + 1
    !
    nr1s=dfftp%nr1; nr2s=dfftp%nr2; nr3s=dfftp%nr3 
    !
    IF (nproc_image .GE. nr3s) THEN
      sizefft=MAX(nr1s*nr2s,nr1s*nr2s*dffts%my_nr3p)
    ELSE
      sizefft=nr1s*nr2s*dffts%my_nr3p
    END IF 
    !
    !write(*,'(10I7)') me, dffts%npp(me), dfftp%npp(me), dffts%tg_npp(me), dfftp%nr1x, dfftp%nr2x, dfftp%nr3x, dffts%nnr, dfftp%nnr, dffts%tg_nnr
    !
    IF (nproc_image .LE. nbsp) THEN
      sc_fac = 1
    ELSE
      sc_fac = nproc_image/nbsp
    END IF
    !
    IF (nproc_image .LE. nbsp) THEN
      !
      ! invfft may encounter overflow when nr3l is not the same throughout all MPI images using psis(sizefft);
      ! therefore, here we add some buffer as dffts%nnr (which is maxval_images(sizefft))
      !
      ALLOCATE ( psis(dffts%nnr)    ); psis=0.0_DP
      ALLOCATE ( psis2(sizefft,nbsp)); psis2=0.0_DP
#ifdef __CUDA
      ALLOCATE ( psis_d(dffts%nnr)    ); psis_d=0.0_DP
      ALLOCATE ( psis2_d(sizefft,nbsp)); psis2_d=0.0_DP
      ALLOCATE ( ca_d(ngw) ); ca_d = (0.0_DP, 0.0_DP)
      ALLOCATE ( c_d, source = c);
#endif
      !
      ca(:)=(0., 0.)
      !
      DO i = 1, nbsp, 2
        !
        IF ( (MOD(nbsp, 2) .NE. 0) .AND. (i .EQ. nbsp ) ) THEN     
#ifdef __CUDA
          CALL fftx_c2psi_gamma_gpu( dffts, psis_d, c_d(:,i), ca_d)
#else
          CALL fftx_c2psi_gamma( dffts, psis, c(:,i:i), ca)
#endif
        ELSE
#ifdef __CUDA
          CALL fftx_c2psi_gamma_gpu( dffts, psis_d, c_d(:,i), c_d(:, i+1))
#else
          CALL fftx_c2psi_gamma( dffts, psis, c(:,i:i), c(:,i+1))
#endif
        END IF 
        !
#ifdef __CUDA
        associate(psis=>psis_d, psis2=>psis2_d)  
#endif
        CALL invfft( 'Wave', psis, dffts )
        !
#ifdef __CUDA
        !$cuf kernel do (1)
#endif
        DO ir = 1,sizefft
          psis2(ir,i) = DBLE(psis(ir))
        END DO
        !
        IF ( (mod(nbsp, 2) .EQ. 0) .OR. (i .NE. nbsp ) ) THEN
#ifdef __CUDA
          !$cuf kernel do (1)
#endif
          DO ir=1,sizefft
            psis2(ir,i+1) = AIMAG(psis(ir))
          END DO
        END IF
        !
#ifdef __CUDA
        end associate
#endif
      END DO !loop over state i

#ifdef __CUDA
      psis2 = psis2_d
#endif
      !
#ifdef __MPI
      CALL redistwfr( psis2, psitot2, my_nxyz, my_nbsp, intra_image_comm, 1 )
#else
      psitot2=psis2
#endif
      !
    ELSE ! (nproc_image .GT. nbsp) 
      !
      !==================================================================
      !
      ! Zhaofeng's code
      !
      !write(stdout,'("entering exx_psi")')
      !write(stdout,'("nogrp:",I5)'),nogrp
      !write(stdout,'("dffts%nnr*nogrp:",I10)'), dffts%nnr*nogrp
      !write(stdout,'("nogrp*nr3 should be smaller or equal to nproc_image:")')
      !
      nogrp = fftx_ntgrp(dffts)
      !
      ALLOCATE( sdispls(nproc_image), sendcount(nproc_image) ); sdispls=0; sendcount=0
      ALLOCATE( rdispls(nproc_image), recvcount(nproc_image) ); rdispls=0; recvcount=0 
      ALLOCATE( sdispls1(nogrp), sendcount1(nogrp) ); sdispls1=0; sendcount1=0
      ALLOCATE( rdispls1(nogrp), recvcount1(nogrp) ); rdispls1=0; recvcount1=0
      !
      DO proc = 1 , nogrp
        sendcount1(proc) = dffts%nnr/nogrp
        recvcount1(proc) = dffts%nnr/nogrp
      END DO
      !
      rdispls1(1) = 0
      sdispls1(1) = 0
      !
      DO proc = 2, nogrp
        sdispls1(proc) = sdispls1(proc-1) + sendcount1(proc-1)
        rdispls1(proc) = rdispls1(proc-1) + recvcount1(proc-1)
      END DO
      !
      ALLOCATE ( psis(dffts%nnr*nogrp) ); psis=0.0_DP
      ALLOCATE ( psis1(dffts%nnr*nogrp) ); psis1=0.0_DP
      ALLOCATE ( psis2(dffts%nnr,nproc_image/nogrp)); psis2=0.0_DP
      !
      jj = 1 ! ??
      !
      va = dffts%nnr/nogrp
      nnr2 = dffts%nnr/2
      !
      !**** nbsp has to be divisible by (2*nogrp) ***
      !
      DO i = 1, nbsp, 2*nogrp
        !
        CALL fftx_c2psi_gamma_tg( dffts, psis, c(:,i:nbsp), ngw, nbsp-i+1 )
        !
        CALL invfft( 'tgWave', psis, dffts )
        !
#if defined(__MPI)
        !
        CALL mp_barrier( fftx_tgcomm(dffts) )
        CALL MPI_ALLTOALLV(psis, sendcount1, sdispls1, MPI_DOUBLE_COMPLEX, &
            &         psis1, recvcount1, rdispls1, MPI_DOUBLE_COMPLEX, &
            &         fftx_tgcomm(dffts), ierr)
#endif
        !
        ngpww1 = 0
        !
        DO ii = 1,2
          DO j = 1,nogrp/2
            ig=(ii-1)*nnr2+(j-1)*va
            !$omp parallel do 
            DO ir = 1,va
              psis2(ir+ngpww1,jj) = DBLE(psis1(ir+ig))
            END DO
            !$omp end parallel do 
            ngpww1=ngpww1+va
            !$omp parallel do 
            DO ir = 1,va
              psis2(ir+ngpww1,jj) = AIMAG(psis1(ir+ig))
            END DO
            !$omp end parallel do 
            ngpww1 = ngpww1 + va
          END DO!loop over j
          jj = jj + 1
          ngpww1 = 0
        END DO!loop over ii
        !
      END DO !loop over state i
      !
      ! ... the wavefunction is duplicated over number of processors that is integer multiple of bands
      ! ... (this part can be improved to distribute wavefunction over any number of processors)
      !
      DO jj=1,nbsp/nogrp
        DO i=1,sc_fac-1
          ig=jj+(i*nbsp/nogrp)
          !$omp parallel do
          DO ir = 1,dffts%nnr
            psis2(ir,ig)=psis2(ir,jj)
          END DO
          !$omp end parallel do 
        END DO
      END DO
      !
      !******************
      !
      ALLOCATE ( psitot(nnrtot*my_nbsp(me)) ); psitot=0.0_DP
      !
      DO proc = 1, nproc_image
        !
        IF (me <= nogrp*nr3s) THEN
          sendcount(proc) =nr1s*nr2s/nogrp
        ELSE
          sendcount(proc) = 0
        END IF
        !
        IF (proc <= nogrp*nr3s) THEN
          recvcount(proc)=nr1s*nr2s/nogrp
        ELSE
          recvcount(proc)=0
        END IF
        !
      END DO
      !
      sdispls(1) = 0
      rdispls(1) = 0
      !
      DO proc = 2,  nproc_image
        sdispls(proc)=  sdispls(proc-1) + sendcount(proc-1)
        rdispls(proc) = rdispls(proc-1) + recvcount(proc-1)
      END DO
      !
#if defined(__MPI)
      !
      CALL mp_barrier( intra_image_comm )
      CALL MPI_ALLTOALLV(psis2, sendcount,sdispls,MPI_DOUBLE_PRECISION, &
          &            psitot, recvcount,rdispls, MPI_DOUBLE_PRECISION, &
          &            intra_image_comm, ierr)
#endif
      !
      !$omp parallel do 
      DO ir = 1, nnrtot
        psitot2(ir,1)=psitot(ir)
      END DO
      !$omp end parallel do 
      !
    END IF ! (nproc_image .LE. nbsp)
    !
    IF (ALLOCATED(psis))         DEALLOCATE(psis)
    IF (ALLOCATED(psis1))        DEALLOCATE(psis1)
    IF (ALLOCATED(psis2))        DEALLOCATE(psis2)
    IF (ALLOCATED(psitot))       DEALLOCATE(psitot)
    !
    IF (ALLOCATED(sdispls))      DEALLOCATE(sdispls)
    IF (ALLOCATED(rdispls))      DEALLOCATE(rdispls)
    IF (ALLOCATED(sdispls1))     DEALLOCATE(sdispls1)
    IF (ALLOCATED(rdispls1))     DEALLOCATE(rdispls1)
    IF (ALLOCATED(sendcount))    DEALLOCATE(sendcount)
    IF (ALLOCATED(recvcount))    DEALLOCATE(recvcount)
    IF (ALLOCATED(sendcount1))   DEALLOCATE(sendcount1)
    IF (ALLOCATED(recvcount1))   DEALLOCATE(recvcount1)
#ifdef __CUDA
    DEALLOCATE (psis_d)
    DEALLOCATE (psis2_d)
    DEALLOCATE (ca_d)
    DEALLOCATE (c_d)
#endif
    !
    !=========================================================================
    !
    RETURN
    !
END SUBROUTINE exx_psi  
