MODULE fft_helper_subroutines

  IMPLICIT NONE
  SAVE

  INTERFACE tg_reduce_rho
    MODULE PROCEDURE tg_reduce_rho_1,tg_reduce_rho_2,tg_reduce_rho_3,tg_reduce_rho_4, &
&                    tg_reduce_rho_5
  END INTERFACE

CONTAINS

  SUBROUTINE tg_reduce_rho_1( rhos, tg_rho_nc, tg_rho, ispin, noncolin, domag, desc )
     USE fft_param
     USE fft_types,      ONLY : fft_type_descriptor

     TYPE(fft_type_descriptor), INTENT(in) :: desc
     INTEGER, INTENT(IN) :: ispin
     LOGICAL, INTENT(IN) :: noncolin, domag
     REAL(DP), INTENT(INOUT)  :: tg_rho(:)
     REAL(DP), INTENT(INOUT)  :: tg_rho_nc(:,:)
     REAL(DP), INTENT(OUT) :: rhos(:,:)

     INTEGER :: ierr, ioff, ir3, ir, ipol, ioff_tg, nxyp, npol_
!     write (*,*) ' enter tg_reduce_rho_1'

#if defined(__MPI)
     IF( noncolin) THEN
        CALL MPI_ALLREDUCE( MPI_IN_PLACE, tg_rho_nc, SIZE(tg_rho_nc), MPI_DOUBLE_PRECISION, MPI_SUM, desc%comm2, ierr )
     ELSE
        CALL MPI_ALLREDUCE( MPI_IN_PLACE, tg_rho, SIZE(tg_rho), MPI_DOUBLE_PRECISION, MPI_SUM, desc%comm2, ierr )
     END IF
#endif
     !
     ! copy the charge back to the proper processor location
     !
     nxyp = desc%nr1x * desc%my_nr2p
     IF (noncolin) THEN
         npol_ = 1 ; if (domag) npol_ = 4
     endif
     DO ir3=1,desc%my_nr3p
        ioff    = desc%nr1x * desc%my_nr2p * (ir3-1)
        ioff_tg = desc%nr1x * desc%nr2x    * (ir3-1) + desc%nr1x * desc%my_i0r2p
        IF (noncolin) THEN
!$omp parallel do
           DO ipol=1, npol_
              DO ir = 1, nxyp
                 rhos(ir+ioff,ipol) = rhos(ir+ioff,ipol) + tg_rho_nc(ir+ioff_tg,ipol)
              END DO
           END DO
!$omp end parallel do
        ELSE
!$omp parallel do
           DO ir = 1, nxyp
              rhos(ir+ioff,ispin) = rhos(ir+ioff,ispin) + tg_rho(ir+ioff_tg)
           END DO
!$omp end parallel do
        END IF
     END DO

  END SUBROUTINE



  SUBROUTINE tg_reduce_rho_2( rhos, tmp_rhos, ispin, desc )
     USE fft_param
     USE fft_types,      ONLY : fft_type_descriptor

     TYPE(fft_type_descriptor), INTENT(in) :: desc
     INTEGER, INTENT(IN) :: ispin
     REAL(DP), INTENT(INOUT)  :: tmp_rhos(:)
     REAL(DP), INTENT(OUT) :: rhos(:,:)

     INTEGER :: ierr, ioff, ir3, nxyp, ioff_tg
!     write (*,*) ' enter tg_reduce_rho_2'

     IF ( desc%nproc2 > 1 ) THEN
#if defined(__MPI)
        CALL MPI_ALLREDUCE( MPI_IN_PLACE, tmp_rhos, SIZE(tmp_rhos), MPI_DOUBLE_PRECISION, MPI_SUM, desc%comm2, ierr )
#endif
     ENDIF
     !
     !BRING CHARGE DENSITY BACK TO ITS ORIGINAL POSITION
     !
     nxyp = desc%nr1x * desc%my_nr2p
     DO ir3 = 1, desc%my_nr3p
        ioff    = desc%nr1x * desc%my_nr2p * (ir3-1)
        ioff_tg = desc%nr1x * desc%nr2x    * (ir3-1) + desc%nr1x * desc%my_i0r2p
        rhos(ioff+1:ioff+nxyp,ispin) = rhos(ioff+1:ioff+nxyp,ispin) + tmp_rhos(ioff_tg+1:ioff_tg+nxyp)
     END DO
 
  END SUBROUTINE

  SUBROUTINE tg_reduce_rho_3( rhos, tmp_rhos, desc )
     USE fft_param
     USE fft_types,      ONLY : fft_type_descriptor

     TYPE(fft_type_descriptor), INTENT(in) :: desc
     REAL(DP), INTENT(INOUT)  :: tmp_rhos(:,:)
     REAL(DP), INTENT(OUT) :: rhos(:,:)

     INTEGER :: ierr, from, ir3, ioff, nxyp, ioff_tg
!     write (*,*) ' enter tg_reduce_rho_3'

     IF ( desc%nproc2 > 1 ) THEN
#if defined(__MPI)
        CALL MPI_ALLREDUCE( MPI_IN_PLACE, tmp_rhos, SIZE(tmp_rhos), MPI_DOUBLE_PRECISION, MPI_SUM, desc%comm2, ierr )
#endif
     ENDIF
     !
     !BRING CHARGE DENSITY BACK TO ITS ORIGINAL POSITION
     !
     !If the current processor is not the "first" processor in its
     !orbital group then does a local copy (reshuffling) of its data
     !
     nxyp = desc%nr1x * desc%my_nr2p
     DO ir3 = 1, desc%my_nr3p
        ioff    = desc%nr1x * desc%my_nr2p * (ir3-1)
        ioff_tg = desc%nr1x * desc%nr2x    * (ir3-1) + desc%nr1x * desc%my_i0r2p
        rhos(ioff+1:ioff+nxyp,:) = rhos(ioff+1:ioff+nxyp,:) + tmp_rhos(ioff_tg+1:ioff_tg+nxyp,:)
     END DO
  END SUBROUTINE


  SUBROUTINE tg_reduce_rho_4( rhos, tmp_rhos, desc )
     USE fft_param
     USE fft_types,      ONLY : fft_type_descriptor

     TYPE(fft_type_descriptor), INTENT(in) :: desc
     COMPLEX(DP), INTENT(INOUT)  :: tmp_rhos(:)
     COMPLEX(DP), INTENT(OUT) :: rhos(:)

     INTEGER :: ierr, from, ir3, ioff, nxyp, ioff_tg
!     write (*,*) ' enter tg_reduce_rho_4'

     IF ( desc%nproc2 > 1 ) THEN
#if defined(__MPI)
        CALL MPI_ALLREDUCE( MPI_IN_PLACE, tmp_rhos, 2*SIZE(tmp_rhos), MPI_DOUBLE_PRECISION, MPI_SUM, desc%comm2, ierr )
#endif
     ENDIF
     !
     !BRING CHARGE DENSITY BACK TO ITS ORIGINAL POSITION
     !
     !If the current processor is not the "first" processor in its
     !orbital group then does a local copy (reshuffling) of its data
     !
     nxyp = desc%nr1x * desc%my_nr2p
     DO ir3 = 1, desc%my_nr3p
        ioff    = desc%nr1x * desc%my_nr2p * (ir3-1)
        ioff_tg = desc%nr1x * desc%nr2x    * (ir3-1) + desc%nr1x * desc%my_i0r2p
        rhos(ioff+1:ioff+nxyp) = rhos(ioff+1:ioff+nxyp) + tmp_rhos(ioff_tg+1:ioff_tg+nxyp)
     END DO
  END SUBROUTINE


  SUBROUTINE tg_reduce_rho_5( rhos, tmp_rhos, desc )
     USE fft_param
     USE fft_types,      ONLY : fft_type_descriptor

     TYPE(fft_type_descriptor), INTENT(in) :: desc
     COMPLEX(DP), INTENT(INOUT)  :: tmp_rhos(:,:)
     COMPLEX(DP), INTENT(OUT) :: rhos(:,:)

     INTEGER :: ierr, from, ir3, ioff, nxyp, ioff_tg
!     write (*,*) ' enter tg_reduce_rho_5'

     IF ( desc%nproc2 > 1 ) THEN
#if defined(__MPI)
        CALL MPI_ALLREDUCE( MPI_IN_PLACE, tmp_rhos, 2*SIZE(tmp_rhos), MPI_DOUBLE_PRECISION, MPI_SUM, desc%comm2, ierr )
#endif
     ENDIF
     !
     !BRING CHARGE DENSITY BACK TO ITS ORIGINAL POSITION
     !
     !If the current processor is not the "first" processor in its
     !orbital group then does a local copy (reshuffling) of its data
     !
     nxyp = desc%nr1x * desc%my_nr2p
     DO ir3 = 1, desc%my_nr3p
        ioff    = desc%nr1x * desc%my_nr2p * (ir3-1)
        ioff_tg = desc%nr1x * desc%nr2x    * (ir3-1) + desc%nr1x * desc%my_i0r2p
        rhos(ioff+1:ioff+nxyp,:) = rhos(ioff+1:ioff+nxyp,:) + tmp_rhos(ioff_tg+1:ioff_tg+nxyp,:)
     END DO
  END SUBROUTINE



  SUBROUTINE tg_get_nnr( desc, right_nnr )
     USE fft_param
     USE fft_types,      ONLY : fft_type_descriptor
     TYPE(fft_type_descriptor), INTENT(in) :: desc
     INTEGER, INTENT(OUT) :: right_nnr
     right_nnr = desc%nnr
  END SUBROUTINE

  SUBROUTINE tg_get_local_nr3( desc, val )
     USE fft_param
     USE fft_types,      ONLY : fft_type_descriptor
     TYPE(fft_type_descriptor), INTENT(in) :: desc
     INTEGER, INTENT(OUT) :: val
     val = desc%my_nr3p
  END SUBROUTINE

  SUBROUTINE tg_get_group_nr3( desc, val )
     USE fft_param
     USE fft_types,      ONLY : fft_type_descriptor
     TYPE(fft_type_descriptor), INTENT(in) :: desc
     INTEGER, INTENT(OUT) :: val
     val = desc%my_nr3p
  END SUBROUTINE

  SUBROUTINE tg_get_recip_inc( desc, val )
     USE fft_param
     USE fft_types,      ONLY : fft_type_descriptor
     TYPE(fft_type_descriptor), INTENT(in) :: desc
     INTEGER, INTENT(OUT) :: val
     val = desc%nnr
  END SUBROUTINE

  PURE FUNCTION fftx_ntgrp( desc )
     USE fft_param
     USE fft_types,      ONLY : fft_type_descriptor
     INTEGER :: fftx_ntgrp
     TYPE(fft_type_descriptor), INTENT(in) :: desc
     fftx_ntgrp = desc%nproc2
  END FUNCTION

  PURE FUNCTION fftx_tgpe( desc )
     USE fft_param
     USE fft_types,      ONLY : fft_type_descriptor
     INTEGER :: fftx_tgpe
     TYPE(fft_type_descriptor), INTENT(in) :: desc
     fftx_tgpe = desc%mype2
  END FUNCTION

  PURE FUNCTION fftx_tgcomm( desc )
     USE fft_param
     USE fft_types,      ONLY : fft_type_descriptor
     INTEGER :: fftx_tgcomm
     TYPE(fft_type_descriptor), INTENT(in) :: desc
     fftx_tgcomm = desc%comm2
  END FUNCTION

  SUBROUTINE fftx_add_field( r, f, desc )
     USE fft_param
     USE fft_types,      ONLY : fft_type_descriptor
     TYPE(fft_type_descriptor), INTENT(in) :: desc
     REAL(DP), INTENT(INOUT)  :: r(:,:)
     REAL(DP), INTENT(IN) :: f(:)
     INTEGER :: nspin, ir
     nspin = SIZE( r, 2 )
     IF (nspin.EQ.1) THEN
!$omp parallel do
        DO ir=1,desc%nr1*desc%nr2*desc%my_nr3p
           r(ir,1)=r(ir,1)+f(ir)
        END DO
!$omp end parallel do
     ELSE IF (nspin.EQ.2) THEN
!$omp parallel do
        DO ir=1,desc%nr1*desc%nr2*desc%my_nr3p
           r(ir,1)=r(ir,1)+f(ir)
        END DO
!$omp end parallel do
!$omp parallel do
        DO ir=1,desc%nr1*desc%nr2*desc%my_nr3p
           r(ir,2)=r(ir,2)+f(ir)
        END DO
!$omp end parallel do
     END IF
  END SUBROUTINE


  SUBROUTINE c2psi_gamma( desc, psi, c, ca )
     !
     !  Copy wave-functions from 1D array (c_bgrp) to 3D array (psi) in Fourier space
     !
     USE fft_param
     USE fft_types,      ONLY : fft_type_descriptor
     TYPE(fft_type_descriptor), INTENT(in) :: desc
     complex(DP), INTENT(OUT) :: psi(:)
     complex(DP), INTENT(IN) :: c(:)
     complex(DP), OPTIONAL, INTENT(IN) :: ca(:)
     complex(DP), parameter :: ci=(0.0d0,1.0d0)
     integer :: ig
     !
     psi = 0.0d0
     !
     !  nlm and nl array: hold conversion indices form 3D to
     !     1-D vectors. Columns along the z-direction are stored
     !     contigiously
     !  c array: stores the Fourier expansion coefficients
     !     Loop for all local g-vectors (ngw)
     IF( PRESENT(ca) ) THEN
        do ig = 1, desc%ngw
           psi( desc%nlm( ig ) ) = CONJG( c( ig ) ) + ci * conjg( ca( ig ))
           psi( desc%nl( ig ) ) = c( ig ) + ci * ca( ig )
        end do
     ELSE
        do ig = 1, desc%ngw
           psi( desc%nlm( ig ) ) = CONJG( c( ig ) )
           psi( desc%nl( ig ) ) = c( ig )
        end do
     END IF
  END SUBROUTINE

  SUBROUTINE c2psi_k( desc, psi, c, igk, ngk)
     !
     !  Copy wave-functions from 1D array (c/evc) ordered according (k+G) index igk 
     !  to 3D array (psi) in Fourier space
     !
     USE fft_param
     USE fft_types,      ONLY : fft_type_descriptor
     TYPE(fft_type_descriptor), INTENT(in) :: desc
     complex(DP), INTENT(OUT) :: psi(:)
     complex(DP), INTENT(IN) :: c(:)
     INTEGER, INTENT(IN) :: igk(:), ngk
     ! local variables
     integer :: ig
     !
     !  nl array: hold conversion indices form 3D to 1-D vectors. 
     !     Columns along the z-direction are stored contigiously
     !  c array: stores the Fourier expansion coefficients of the wave function
     !     Loop for all local g-vectors (npw
     psi = 0.0d0
     do ig = 1, ngk
        psi( desc%nl( igk( ig ) ) ) = c( ig )
     end do
     !
  END SUBROUTINE

  SUBROUTINE fftx_oned2threed( desc, psi, c, ca )
     !
     !  Copy charge density from 1D array (c) to 3D array (psi) in Fourier
     !  space
     !
     USE fft_param
     USE fft_types,      ONLY : fft_type_descriptor
     TYPE(fft_type_descriptor), INTENT(in) :: desc
     complex(DP), INTENT(OUT) :: psi(:)
     complex(DP), INTENT(IN) :: c(:)
     complex(DP), OPTIONAL, INTENT(IN) :: ca(:)
     complex(DP), parameter :: ci=(0.0d0,1.0d0)
     integer :: ig
     !
     psi = 0.0d0
     !
     !  nlm and nl array: hold conversion indices form 3D to
     !     1-D vectors. Columns along the z-direction are stored
     !     contigiously
     !  c array: stores the Fourier expansion coefficients
     !     Loop for all local g-vectors (ngw)
     IF( PRESENT(ca) ) THEN
        IF( desc%lgamma ) THEN
           do ig = 1, desc%ngm
              psi( desc%nlm( ig ) ) = CONJG( c( ig ) ) + ci * conjg( ca( ig ))
              psi( desc%nl( ig ) ) = c( ig ) + ci * ca( ig )
           end do
        ELSE
           do ig = 1, desc%ngm
              psi( desc%nl( ig ) ) = c( ig ) + ci * ca( ig )
           end do
        END IF
     ELSE
        IF( desc%lgamma ) THEN
           do ig = 1, desc%ngm
              psi( desc%nlm( ig ) ) = CONJG( c( ig ) )
              psi( desc%nl( ig ) ) = c( ig )
           end do
        ELSE
           DO ig = 1, desc%ngm
              psi( desc%nl( ig ) ) = c( ig )
           END DO
        END IF
     END IF
  END SUBROUTINE

  SUBROUTINE fftx_add_threed2oned_gamma( desc, vin, vout1, vout2 )
     USE fft_param
     USE fft_types,      ONLY : fft_type_descriptor
     TYPE(fft_type_descriptor), INTENT(in) :: desc
     complex(DP), INTENT(INOUT) :: vout1(:)
     complex(DP), OPTIONAL, INTENT(INOUT) :: vout2(:)
     complex(DP), INTENT(IN) :: vin(:)
     COMPLEX(DP) :: fp, fm
     INTEGER :: ig
     IF( PRESENT( vout2 ) ) THEN
        DO ig=1,desc%ngm
           fp=vin(desc%nl(ig))+vin(desc%nlm(ig))
           fm=vin(desc%nl(ig))-vin(desc%nlm(ig))
           vout1(ig) = vout1(ig) + 0.5d0*CMPLX( DBLE(fp),AIMAG(fm),kind=DP)
           vout2(ig) = vout2(ig) + 0.5d0*CMPLX(AIMAG(fp),-DBLE(fm),kind=DP)
        END DO
     ELSE
        DO ig=1,desc%ngm
           vout1(ig) = vout1(ig) + vin(desc%nl(ig))
        END DO
     END IF
  END SUBROUTINE

  SUBROUTINE fftx_threed2oned( desc, vin, vout1, vout2 )
     USE fft_param
     USE fft_types,      ONLY : fft_type_descriptor
     TYPE(fft_type_descriptor), INTENT(in) :: desc
     complex(DP), INTENT(OUT) :: vout1(:)
     complex(DP), OPTIONAL, INTENT(OUT) :: vout2(:)
     complex(DP), INTENT(IN) :: vin(:)
     COMPLEX(DP) :: fp, fm
     INTEGER :: ig
     IF( PRESENT( vout2 ) ) THEN
        DO ig=1,desc%ngm
           fp=vin(desc%nl(ig))+vin(desc%nlm(ig))
           fm=vin(desc%nl(ig))-vin(desc%nlm(ig))
           vout1(ig) = 0.5d0*CMPLX( DBLE(fp),AIMAG(fm),kind=DP)
           vout2(ig) = 0.5d0*CMPLX(AIMAG(fp),-DBLE(fm),kind=DP)
        END DO
     ELSE
        DO ig=1,desc%ngm
           vout1(ig) = vin(desc%nl(ig))
        END DO
     END IF
  END SUBROUTINE


  SUBROUTINE fftx_psi2c_gamma( desc, vin, vout1, vout2 )
     USE fft_param
     USE fft_types,      ONLY : fft_type_descriptor
     TYPE(fft_type_descriptor), INTENT(in) :: desc
     complex(DP), INTENT(OUT) :: vout1(:)
     complex(DP), OPTIONAL, INTENT(OUT) :: vout2(:)
     complex(DP), INTENT(IN) :: vin(:)
     COMPLEX(DP) :: fp, fm
     INTEGER :: ig
     IF( PRESENT( vout2 ) ) THEN
        DO ig=1,desc%ngw
           fp=vin(desc%nl(ig))+vin(desc%nlm(ig))
           fm=vin(desc%nl(ig))-vin(desc%nlm(ig))
           vout1(ig) = CMPLX( DBLE(fp),AIMAG(fm),kind=DP)
           vout2(ig) = CMPLX(AIMAG(fp),-DBLE(fm),kind=DP)
        END DO
     ELSE
        DO ig=1,desc%ngw
           vout1(ig) = vin(desc%nl(ig))
        END DO
     END IF
  END SUBROUTINE


  SUBROUTINE c2psi_gamma_tg(desc, psis, c_bgrp, i, nbsp_bgrp )
     !
     !  Copy all wave-functions of an orbital group 
     !  from 1D array (c_bgrp) to 3D array (psi) in Fourier space
     !
     USE fft_param
     USE fft_types,      ONLY : fft_type_descriptor
     TYPE(fft_type_descriptor), INTENT(in) :: desc
     complex(DP), INTENT(OUT) :: psis(:)
     complex(DP), INTENT(INOUT) :: c_bgrp(:,:)
     INTEGER, INTENT(IN) :: i, nbsp_bgrp
     INTEGER :: eig_offset, eig_index, right_nnr
     !
     !  the i-th column of c_bgrp corresponds to the i-th state (in this band group)
     !
     !  The outer loop goes through i : i + 2*NOGRP to cover
     !  2*NOGRP eigenstates at each iteration
     !
     CALL tg_get_nnr( desc, right_nnr )

!$omp  parallel
!$omp  single

     do eig_index = 1, 2 * fftx_ntgrp(desc), 2

!$omp task default(none) &
!$omp          firstprivate( eig_index, i, nbsp_bgrp, right_nnr ) &
!$omp          private( eig_offset ) &
!$omp          shared( c_bgrp, desc, psis )
        !
        !  here we pack 2*nogrp electronic states in the psis array
        !  note that if nogrp == nproc_bgrp each proc perform a full 3D
        !  fft and the scatter phase is local (without communication)
        !
        !  important: if n is odd => c(*,n+1)=0.
        !
        IF ( ( eig_index + i - 1 ) == nbsp_bgrp ) c_bgrp( : , eig_index + i ) = 0.0d0
        !
        eig_offset = ( eig_index - 1 )/2
        !
        IF ( ( i + eig_index - 1 ) <= nbsp_bgrp ) THEN
           !
           !  The  eig_index loop is executed only ONCE when NOGRP=1.
           !
           CALL c2psi_gamma( desc, psis( eig_offset * right_nnr + 1 : eig_offset * right_nnr + right_nnr ), &
                       c_bgrp( :, i+eig_index-1 ), c_bgrp( :, i+eig_index ) )
           !
        ENDIF
!$omp end task
        !
     end do
!$omp  end single
!$omp  end parallel
     RETURN
  END SUBROUTINE


  SUBROUTINE fft_dist_info( desc, unit )
     USE fft_param
     USE fft_types,      ONLY : fft_type_descriptor
     INTEGER, INTENT(IN) :: unit
     TYPE(fft_type_descriptor), INTENT(in) :: desc
     INTEGER :: i, j, nr3l
     CALL tg_get_local_nr3( desc, nr3l )
     WRITE( stdout,1000) desc%nr1, desc%nr2, desc%nr3, &
                         desc%nr1, desc%my_nr2p, desc%my_nr3p, &
                         1, desc%nproc2, desc%nproc3
     WRITE( stdout,1010) desc%nr1x, desc%nr2x, desc%nr3x
     WRITE( stdout,1020) desc%nnr
     WRITE( stdout,*) '  Number of x-y planes for each processors: '
     WRITE( stdout, fmt = '( 5("  |",I4,",",I4) )' ) ( ( desc%nr2p(j), &
             desc%nr3p(i), i = 1, desc%nproc3 ), j=1,desc%nproc2 )
1000  FORMAT(3X, &
         'Global Dimensions   Local  Dimensions   Processor Grid',/,3X, &
         '.X.   .Y.   .Z.     .X.   .Y.   .Z.     .X.   .Y.   .Z.',/, &
         3(1X,I5),2X,3(1X,I5),2X,3(1X,I5) )
1010  FORMAT(3X, 'Array leading dimensions ( nr1x, nr2x, nr3x )   = ', 3(1X,I5))
1020  FORMAT(3X, 'Local number of cell to store the grid ( nrxx ) = ', 1X, I9 )
     RETURN
  END SUBROUTINE

END MODULE fft_helper_subroutines
