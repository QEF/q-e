MODULE direct_www

USE kinds, ONLY : DP
USE bse_basic_structures

SAVE

type(ii_mat)  :: iimat_direct
type(www_mat) :: wwwmat_direct

CONTAINS

SUBROUTINE free_memory_direct_www
   IMPLICIT NONE
   CALL free_www_mat(wwwmat_direct)
   CALL free_imat(iimat_direct)
   RETURN
END SUBROUTINE

SUBROUTINE initialize_direct_www(fc)
   USE fft_custom_gwl
   USE gvect
   USE lsda_mod, ONLY: nspin
   USE wvfct,    ONLY: npw
   USE wavefunctions, ONLY: psic
   USE mp_world, ONLY: mpime, nproc
   USE mp_pools, ONLY: intra_pool_comm
   
   IMPLICIT NONE
   TYPE(fft_cus), INTENT(IN)  :: fc

   INTEGER :: is, ii, loop_iimat
   COMPLEX(KIND=DP), ALLOCATABLE :: www_t(:,:)
   LOGICAL, PARAMETER :: debug=.FALSE.

   CALL initialize_imat(iimat_direct) 
   CALL initialize_www_mat(wwwmat_direct) 

   DO is = 1, nspin
      CALL read_iimat(iimat_direct, is)
      IF (debug) THEN
         write(*,*) 'I have read iimat into iimat_direct:'
         DO loop_iimat = 1, iimat_direct%numb_v
            write(*,*) iimat_direct%iimat(1:iimat_direct%np_max,loop_iimat)
         END DO
      END IF
   END DO

   CALL read_www_mat(iimat_direct, wwwmat_direct)
   IF (debug) THEN
      write(*,*) 'I have compiled wwwmat_direct:'
      DO loop_iimat = 1, wwwmat_direct%numb_v
         write(*,*) wwwmat_direct%ii_www_mat(1:wwwmat_direct%np_max,loop_iimat)
      END DO
      write(*,*) 'Sanity Check (iimat_direct): '
      DO loop_iimat = 1, iimat_direct%numb_v
         write(*,*) iimat_direct%iimat(1:iimat_direct%np_max,loop_iimat)
      END DO
   END IF

! FFT the Www products
   ALLOCATE(www_t(fc%npwt,wwwmat_direct%ww_tot))

   CALL reorderwfp_col(wwwmat_direct%ww_tot,npw,fc%npwt,wwwmat_direct%www,www_t,npw,fc%npwt, &
           & ig_l2g,fc%ig_l2gt,fc%ngmt_g,mpime, nproc,intra_pool_comm )

   ALLOCATE(wwwmat_direct%www_r(fc%nrxxt,wwwmat_direct%ww_tot))

   do ii=1,wwwmat_direct%ww_tot,2
      psic(1:fc%nrxxt)=(0.d0,0.d0)
      if (ii==wwwmat_direct%ww_tot) then
         psic(fc%nlt(1:fc%npwt))  = www_t(1:fc%npwt,ii)
         psic(fc%nltm(1:fc%npwt)) = CONJG( www_t(1:fc%npwt,ii) )
      else
         psic(fc%nlt(1:fc%npwt))=www_t(1:fc%npwt,ii)+(0.d0,1.d0)*www_t(1:fc%npwt,ii+1)
         psic(fc%nltm(1:fc%npwt))=CONJG(www_t(1:fc%npwt,ii))+(0.d0,1.d0)*CONJG(www_t(1:fc%npwt,ii+1))
      endif
         CALL cft3t( fc, psic, fc%nr1t, fc%nr2t, fc%nr3t, fc%nrx1t, fc%nrx2t, fc%nrx3t, 2 )
          wwwmat_direct%www_r(1:fc%nrxxt,ii)= DBLE(psic(1:fc%nrxxt))
         if(ii/=wwwmat_direct%ww_tot) wwwmat_direct%www_r(1:fc%nrxxt,ii+1)= DIMAG(psic(1:fc%nrxxt))
   enddo

   DEALLOCATE(www_t)
   RETURN
END SUBROUTINE

SUBROUTINE direct_www_exc(a_in, fc, a_out)
   USE bse_basic_structures
   USE exciton
   USE fft_custom_gwl

   IMPLICIT NONE

   TYPE(exc),     INTENT(IN)  :: a_in
   TYPE(fft_cus), INTENT(IN)  :: fc
   TYPE(exc),     INTENT(OUT) :: a_out

   TYPE(exc_r) :: a_in_rt, a_tmp_rt 
   INTEGER, EXTERNAL :: find_free_unit
   INTEGER :: iv,jv, loop_iimat
   LOGICAL, PARAMETER :: debug=.FALSE.

   IF (debug) write(*,*) 'direct_www_exc: DEBUG MODE'

   IF (debug) THEN
      write(*,*) 'Checking wwwmat_direct:'
      DO loop_iimat = 1, wwwmat_direct%numb_v
         write(*,*) wwwmat_direct%ii_www_mat(1:wwwmat_direct%np_max,loop_iimat)
      END DO

      write(*,*) 'Checking iimat_direct: '
      DO loop_iimat = 1, iimat_direct%numb_v
         write(*,*) iimat_direct%iimat(1:iimat_direct%np_max,loop_iimat)
      END DO
   END IF

   ! FFT exciton to real space
   CALL initialize_exc_r(a_in_rt)
   CALL fft_a_exc(a_in, fc, a_in_rt)

   ! Allocate tmp arrays
   CALL initialize_exc_r(a_tmp_rt)
   a_tmp_rt%nrxxt= fc%nrxxt
   a_tmp_rt%numb_v=a_in%numb_v
   a_tmp_rt%label=find_free_unit()
   ALLOCATE(a_tmp_rt%ar(a_tmp_rt%nrxxt, a_tmp_rt%numb_v))
   a_tmp_rt%ar=0.d0

   ! Main loop on iv (v^prime)
   DO iv = 1, a_in%numb_v
      ! ii_mat contains products
      DO jv = 1, iimat_direct%np_max
         ! The product between exc and wannier products
         IF (iimat_direct%iimat(jv,iv) > 0) THEN
            a_tmp_rt%ar(1:a_tmp_rt%nrxxt,iv) = a_tmp_rt%ar(1:a_tmp_rt%nrxxt,iv) &
             + a_in_rt%ar(1:a_tmp_rt%nrxxt,iimat_direct%iimat(jv,iv)) &
             * wwwmat_direct%www_r(1:a_tmp_rt%nrxxt,wwwmat_direct%ii_www_mat(jv,iv))
         ELSE
            
            EXIT
         END IF
      END DO
   END DO

   ! FFT back to G space and clean-up
   CALL fftback_a_exc(a_tmp_rt, fc, a_out)
   CALL free_memory_exc_a_r(a_in_rt)
   CALL free_memory_exc_a_r(a_tmp_rt)

   IF (debug) write(*,*) 'direct_www_exc: EXIT'
   RETURN
END SUBROUTINE

END MODULE
