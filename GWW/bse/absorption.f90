subroutine absorption(vstate_r,psibar,fc,ieig,ampl,ipol)
!this subroutine handles the computation of the absorption spectrum

USE bse_wannier, ONLY : l_finite,r_pola,num_nbndv
USE cell_base,   ONLY : at,alat
USE bse_basic_structures
USE fft_custom_gwl
USE wvfct,       ONLY : npw


implicit none
TYPE(v_state_r),INTENT(in) :: vstate_r
TYPE(fft_cus),  INTENT(in) :: fc
REAL(DP),       INTENT(out):: ampl
INTEGER,        INTENT(in) :: ieig,ipol
COMPLEX(DP),    INTENT(in) :: psibar(npw,num_nbndv(1))

REAL(DP) ::  imod_rpola
REAL(DP) :: upol(3,3)

data upol /1.0d0,0.0d0,0.0d0, 0.0d0,1.0d0,0.0d0, 0.0d0,0.0d0,1.0d0/

call start_clock('absorption')
if(l_finite) then
   r_pola(1:3)=upol(1:3,ipol)
   call amplitude_finite(vstate_r,fc,ieig,ampl)
else
   call amplitude(psibar(1,1),fc,ieig,ampl)
endif

call stop_clock('absorption')
return
end subroutine

subroutine amplitude_finite(vstate_r,fc,ieig,ampl)
!this subroutine computes the amplitude of each exciton for finite systems using 
!the matrix element of the position operator
!WARNING: for now, it should be used only for calculation where the molecule is centered around the
!origin of the supercell 

USE exciton
USE bse_wannier, ONLY : num_nbndv,l_finite,r_pola
USE cell_base,   ONLY : at,alat
USE fft_custom_gwl
USE bse_basic_structures
USE mp_world,   ONLY : mpime, nproc
USE mp_pools,    ONLY: nproc_pool
USE mp,          ONLY : mp_barrier
USE io_global,   ONLY : stdout
USE wvfct,       ONLY : npw
USE mp_world,             ONLY : world_comm


implicit none
TYPE(v_state_r), INTENT(in) :: vstate_r
TYPE(fft_cus),   INTENT(in) :: fc
REAL(DP),        INTENT(out):: ampl
INTEGER,         INTENT(in) :: ieig

TYPE(exc_r) :: rpsiv_r
TYPE(exc)   :: rpsiv

REAL(DP) :: r(3),rr(3),rdote
INTEGER  :: ikstart,iklocal,ij,ii,ik,ifft,iktotal
INTEGER  :: iv
LOGICAL  :: debug
INTEGER  :: iz,iy,ix,iqq
REAL(DP) :: prod,prod2

call start_clock('amplitude_finite')
debug=.true.

! each processor finds the starting index of its collection of FFT planes  
ikstart=1
do ii=1,mpime
   ikstart=ikstart+fc%dfftt%npp(ii)
enddo

iktotal=0
do ii=1,nproc_pool
   iktotal=iktotal+fc%dfftt%npp(ii)
enddo

!if(debug) then
!   write(stdout,*) 'mpime, iktotal=',mpime,iktotal
!   write(stdout,*) 'mpime,fc%nrx3t',mpime,fc%nrx3t
!   write(stdout,*) 'ikstart=',ikstart
!endif

!call flush_unit( stdout )
!call mp_barrier


! create the rpsiv_r excitonic vector (in real space)

call initialize_exc_r(rpsiv_r)
rpsiv_r%nrxxt=fc%nrxxt 
rpsiv_r%numb_v=num_nbndv(1)
rpsiv_r%label=12
allocate(rpsiv_r%ar(rpsiv_r%nrxxt,rpsiv_r%numb_v))

rpsiv_r%ar(1:rpsiv_r%nrxxt,1:rpsiv_r%numb_v)=0.d0

do iz=1,fc%dfftt%npp(mpime+1) 
   do iy=1,fc%nr2t
      do ix=1,fc%nr1t
          iqq=(iz-1)*(fc%nrx1t*fc%nrx2t)+(iy-1)*fc%nrx1t+ix 
          r(:)= (dble(ix-1)/dble(fc%nr1t)-int(2.d0*dble(ix-1)/dble(fc%nr1t)))*at(:,1)*alat+&
            &(dble(iy-1)/dble(fc%nr2t)-int(2.d0*dble(iy-1)/dble(fc%nr2t)))*at(:,2)*alat+&
            &(dble(iz-1+ikstart-1)/dble(fc%nr3t)-int(2.d0*dble(iz-1+ikstart-1)/dble(fc%nr3t)))*at(:,3)*alat

!          r(:)=dble(ix-1)/dble(fc%nr1t)*at(:,1)*alat+&
!               dble(iy-1)/dble(fc%nr2t)*at(:,2)*alat+&
!               dble(iz-1+ikstart-1)/dble(fc%nr3t)*at(:,3)*alat

!          if(debug) then
!             rr(1)=dble(ix-1)/dble(fc%nr1t)-int(2.d0*dble((ix)-1)/dble(fc%nr1t))
!             rr(2)=dble(iy-1)/dble(fc%nr2t)-int(2.d0*dble((iy)-1)/dble(fc%nr2t))
!             rr(3)= dble(iz-1+ikstart-1)/dble(fc%nr3t)-int(2.d0*dble(iz-1+ikstart-1)/dble(fc%nr3t))
!             write(stdout,*) 'rr',rr(1),rr(2),rr(3) 
!             write(stdout,*) 'rc',r(1),r(2),r(3)
!             CALL flush_unit( stdout )   
!          endif

          rdote=r(1)*r_pola(1)+r(2)*r_pola(2)+r(3)*r_pola(3)   

          rpsiv_r%ar(iqq,1:num_nbndv(1))=rdote*vstate_r%wfnrt(iqq,1:num_nbndv(1),1)
  
      enddo
   enddo
enddo

!do ifft=0,fc%nrx1t*fc%nrx2t*fc%dfftt%npp(mpime+1)-1
!
!   iklocal=ifft/(fc%nrx1t*fc%nrx2t)+1
!   ik=ikstart+iklocal-1
!   ij=(ifft-(fc%nrx1t*fc%nrx2t)*(iklocal-1))/fc%nrx1t+1
!   ii=ifft-(fc%nrx1t*fc%nrx2t)*(iklocal-1)-fc%nrx1t*(ij-1)
!!   
!   r(:)= (dble(ii-1)/dble(fc%nrx1t)-int(2.d0*dble((ii)-1)/dble(fc%nrx1t)))*at(:,1)*alat+&
!        &(dble(ij-1)/dble(fc%nrx2t)-int(2.d0*dble((ij)-1)/dble(fc%nrx2t)))*at(:,2)*alat+&
!        &(dble(ik-1)/dble(iktotal)-int(2.d0*dble((ik)-1)/dble(iktotal)))*at(:,3)*alat
!!
!   if(debug) then
!      rr(1)=(dble((ii)-1)/dble(fc%nrx1t)-int(2.d0*dble((ii)-1)/dble(fc%nrx1t)))
!      rr(2)=(dble((ij)-1)/dble(fc%nrx2t)-int(2.d0*dble((ij)-1)/dble(fc%nrx2t)))
!      rr(3)=(dble((ik)-1)/dble(iktotal)-int(2.d0*dble((ik)-1)/dble(iktotal)))
!      write(stdout,*) 'rr',rr(1),rr(2),rr(3) 
!      write(stdout,*) 'rc',r(1),r(2),r(3)
!      CALL flush_unit( stdout )   
!   endif
!
!   rdote=r(1)*r_pola(1)+r(2)*r_pola(2)+r(3)*r_pola(3)   
!
!   rpsiv_r%ar(ifft+1,1:num_nbndv(1))=rdote*vstate_r%wfnrt(ifft+1,1:num_nbndv(1),1)
!
!enddo



!fft rpsiv_r to reciprocal space
call initialize_exc(rpsiv)
rpsiv%label=100
rpsiv%npw=npw
rpsiv%numb_v=num_nbndv(1)
allocate(rpsiv%a(rpsiv%npw,rpsiv%numb_v)) 

!if (debug) then
!   call mp_barrier
!   write(stdout,*) 'rpsiv allocated'
!   CALL flush_unit( stdout )   
!endif

call fftback_a_exc(rpsiv_r,fc,rpsiv)

!if (debug) then
!   call mp_barrier
!   write(stdout,*) 'fft_performed'
!   CALL flush_unit( stdout )   
!endif


!compute the exciton amplitude
if(debug) then
 call sproduct_exc(rpsiv,rpsiv,prod)
 call sproduct_exc(bse_spectrum(ieig),bse_spectrum(ieig),prod2)
 write(*,*) 'ieig, prod1', ieig, prod
 write(*,*) 'ieig, prod2', ieig, prod2
endif

call sproduct_exc(bse_spectrum(ieig),rpsiv,ampl)
ampl=ampl*ampl


!if (debug) then
!   call mp_barrier
!   write(stdout,*) 'amplitude computed'
!   CALL flush_unit( stdout )   
!endif

FLUSH( stdout )   
call free_memory_exc_a_r(rpsiv_r)
call free_memory_exc_a(rpsiv)

call stop_clock('amplitude_finite')
return
end subroutine

subroutine amplitude(psibar,fc,ieig,ampl)
!this subroutine computes the amplitude of each exciton 

USE exciton
USE bse_wannier, ONLY : num_nbndv,l_finite,r_pola
USE cell_base,   ONLY : at,alat
USE fft_custom_gwl
USE bse_basic_structures
USE mp_world,   ONLY : mpime, nproc
USE mp_pools,    ONLY:nproc_pool
USE mp,          ONLY : mp_barrier
USE mp_world,             ONLY : world_comm
USE io_global,   ONLY : stdout
USE wvfct,       ONLY : npw


implicit none
TYPE(fft_cus),   INTENT(in) :: fc
REAL(DP),        INTENT(out):: ampl
INTEGER,         INTENT(in) :: ieig
COMPLEX(DP),    INTENT(in) :: psibar(npw,num_nbndv(1))


REAL(DP) :: r(3),rr(3),rdote,prod,prod2
INTEGER  :: ikstart,iklocal,ij,ii,ik,ifft,iktotal
INTEGER  :: iv
LOGICAL  :: debug
INTEGER  :: iz,iy,ix,iqq
TYPE(exc)   :: rpsiv

call start_clock('amplitude')
debug=.false.

call initialize_exc(rpsiv)
rpsiv%label=100
rpsiv%npw=npw
rpsiv%numb_v=num_nbndv(1)
allocate(rpsiv%a(rpsiv%npw,rpsiv%numb_v)) 

do iv=1,num_nbndv(1)
   rpsiv%a(1:rpsiv%npw,iv)=psibar(1:npw,iv) 
enddo

!check if there is something in the psibar vector
! and in the bse_spectrum_vector
if(debug) then
 call sproduct_exc(rpsiv,rpsiv,prod)
 call sproduct_exc(bse_spectrum(ieig),bse_spectrum(ieig),prod2)
 write(*,*) 'ieig, prod1', ieig, prod
 write(*,*) 'ieig, prod2', ieig, prod2
endif


!compute the exciton amplitude

call sproduct_exc(bse_spectrum(ieig),rpsiv,ampl)
ampl=ampl*ampl


!if (debug) then
!   call mp_barrier
!   write(stdout,*) 'amplitude computed'
!endif

call free_memory_exc_a(rpsiv)

FLUSH( stdout )   
call mp_barrier(world_comm)

call stop_clock('amplitude')
return
end subroutine
