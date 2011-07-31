
!====================================================================
!   get the wave function in real space defined by the FFT grids
!====================================================================

       SUBROUTINE exx_psi(c, psitot2,nnrtot,my_nbsp, my_nxyz, nbsp)

       USE kinds,                   ONLY  : DP
       USE fft_interfaces,           ONLY  : invfft
       USE fft_base,                ONLY  : dffts
       USE gvecw,                   ONLY  : ngw
       USE mp_global,               ONLY  : nproc_image, me_image,intra_image_comm
       USE cell_base,               ONLY  : omega
       USE parallel_include
       USE mp,                      ONLY  : mp_barrier
       USE mp_wave,                  ONLY : redistwfr

       implicit none

       INTEGER,     INTENT(IN)    :: nnrtot, nbsp
       INTEGER,     INTENT(IN)    :: my_nbsp(nproc_image), my_nxyz(nproc_image)
       COMPLEX(DP), INTENT(IN)    :: c(ngw,nbsp)
       REAL(DP),    INTENT(INOUT) :: psitot2(nnrtot, my_nbsp(me_image+1) )
       
       INTEGER i, ir,proc, me, ierr, ig, irank, iobtl

       COMPLEX(DP), ALLOCATABLE ::    psis(:)
       REAL(DP),    ALLOCATABLE ::    psis2(:,:)
       COMPLEX(DP)  ca(ngw)
       complex(DP), parameter   ::    ci=(0.0d0,1.0d0)

       INTEGER      sizefft

       me=me_image + 1
       IF(nproc_image >= dffts%nr3x) THEN
         sizefft=dffts%nnr
       ELSE
         sizefft=dffts%npp(me)*dffts%nr1x*dffts%nr2x
       END IF 

       allocate ( psis(sizefft) )
       allocate ( psis2(sizefft,nbsp))

       print *, "dffts%nnr is", dffts%nnr
       print *, "nbsp and ngw is",nbsp,ngw
       print *, "nrs is", dffts%nr1x, dffts%nr2x, dffts%nr3x

!==================================================================

       ca(:)=(0., 0.)
       DO  i = 1, nbsp, 2
           if( (mod(nbsp, 2) .ne. 0) .and. (i .eq. nbsp ) )then     
              call c2psi( psis, dffts%nnr, c(1,i), ca(1), ngw, 2)
           else
              call c2psi( psis, dffts%nnr, c(1,i), c(1, i+1), ngw, 2)
           end if 

           call invfft( 'Wave', psis, dffts )

           do ir = 1,sizefft
              psis2(ir,i)  =DBLE(psis(ir))
           end do

           if( (mod(nbsp, 2) .eq. 0) .or. (i .ne. nbsp ) )then
              do ir=1,sizefft
                 psis2(ir,i+1)=AIMAG(psis(ir))
              end do
           end if
        END DO !loop over state i

       call redistwfr( psis2, psitot2, my_nxyz, my_nbsp, intra_image_comm, 1 )

       deallocate(psis2, psis )
!=========================================================================
      end subroutine exx_psi  
