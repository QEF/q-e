module exciton
! this module cointains variables and subroutines related to 
! the excitonic wave functions and energies

USE kinds, ONLY : DP

type exc
!excitonic wavefunction vector in G space
    integer npw ! number of plane waves per processor
    integer numb_v ! number of valnce state
!    real(kind=dp), dimension (:), pointer :: e_ks ! Kohn-Sham energy of valence states  
    complex(kind=dp), dimension (:,:), pointer :: a ! vector on which Hexc can be applied  a(npw,nbnd_v)
    INTEGER :: label!label to read/write to disk
    real(kind=DP) :: e ! energy of the excitonic eigenstate
end type exc

type exc_r
!excitonic wavefunction vector in real space (double grid)
    integer nrxxt ! number of plane waves per processor
    integer numb_v ! number of valnce state
!    real(kind=dp), dimension (:), pointer :: e_ks ! Kohn-Sham energy of valence states  
    real(kind=dp), dimension (:,:), pointer :: ar ! vector on which Hexc can be applied  a(npw,nbnd_v)
    INTEGER :: label!label to read/write to disk
end type exc_r

type spectrum
!collection of eigenvalues and optical amplitudes 
    integer neig
    real(kind=dp), dimension(:), pointer:: en 
    real(kind=dp), dimension(:,:), pointer:: a
end type spectrum


type(exc), allocatable :: bse_spectrum(:)
!type(exc), dimension(:), pointer :: bse_spectrum
!this variable contains the excitonic eigenvectors and eigenvalues
!to be deleted probably for memory issues



contains

   SUBROUTINE initialize_exc(a)
   !this subroutine initializes exc
      implicit none
      TYPE(exc) :: a
      nullify(a%a)
!      nullify(a%e_ks)
      return
   END SUBROUTINE 

   SUBROUTINE initialize_exc_r(a_r)
   !this subroutine initializes exc_r
      implicit none
      TYPE(exc_r) :: a_r
      nullify(a_r%ar)
!      nullify(a%e_ks)
      return
   END SUBROUTINE 

   SUBROUTINE initialize_spectrum(s)
   !this subroutine initializes spectrum
      implicit none
      TYPE(spectrum) :: s
      nullify(s%en)
      nullify(s%a)
      return
   END SUBROUTINE 

   SUBROUTINE free_memory_exc_a(a)
   !this subroutine deallocates exc
      implicit none
      TYPE(exc) a
      if(associated(a%a)) deallocate(a%a)
      nullify(a%a)
!      if(associated(a%e_ks)) deallocate(a%e_ks)
!      nullify(a%e_ks)
      return
   END SUBROUTINE

   SUBROUTINE free_memory_exc_a_r(a_r)
   !this subroutine deallocates exc_r 
      implicit none
      TYPE(exc_r) a_r
      if(associated(a_r%ar)) deallocate(a_r%ar)
      nullify(a_r%ar)
      return
   END SUBROUTINE

   SUBROUTINE free_memory_spectrum(s)
   !this subroutine deallocates exc
      implicit none
      TYPE(spectrum) s
      if(associated(s%en)) deallocate(s%en)
      nullify(s%en)
      if(associated(s%a)) deallocate(s%a)
      nullify(s%a)
      return
   END SUBROUTINE


    SUBROUTINE write_exc(a)
    !this subroutine writes the excitonic vectors on disk
    !the file name is taken from the label
!    USE io_files,             ONLY : find_free_unit, prefix
    USE io_files,             ONLY : prefix,tmp_dir
    USE mp_world,  ONLY : mpime
    implicit none
    INTEGER, EXTERNAL :: find_free_unit
    TYPE(exc) :: a!the exc wavefunction to be written

    INTEGER :: iw, jw, iuna
    CHARACTER(5) :: nfile
    CHARACTER(5) :: nproc

    if(a%label >= 0 ) then
      write(nfile,'(5i1)') &
        & a%label/10000,mod(a%label,10000)/1000,mod(a%label,1000)/100,mod(a%label,100)/10,mod(a%label,10)
      write(nproc,'(5i1)') &
        & mpime/10000,mod(mpime,10000)/1000,mod(mpime,1000)/100,mod(mpime,100)/10,mod(mpime,10)
      iuna = find_free_unit()
      open( unit=iuna, file=trim(tmp_dir)//trim(prefix)//'-exc_a.'// nfile //'.'// nproc , &
           &status='unknown',form='unformatted')
    else
      write(nfile,'(5i1)') &
        & -a%label/10000,mod(-a%label,10000)/1000,mod(-a%label,1000)/100,mod(-a%label,100)/10,mod(-a%label,10)
      write(nproc,'(5i1)') &
        & mpime/10000,mod(mpime,10000)/1000,mod(mpime,1000)/100,mod(mpime,100)/10,mod(mpime,10)
      iuna = find_free_unit()
      open( unit=iuna, file=trim(tmp_dir)//trim(prefix)//'-exc_a.-'// nfile //'.'// nproc , &
           &status='unknown',form='formatted')
    endif
    write(iuna) a%label
    write(iuna) a%npw
    write(iuna) a%numb_v
    write(iuna) a%e
    do iw=1,a%numb_v
       write(iuna)  a%a(1:a%npw,iw)
    enddo
    close(iuna)

    return
    END SUBROUTINE

    SUBROUTINE  read_exc(label, a,l_verbose)
    !this subroutine reads the excitonic vectors from disk
    !the file name is taken from the label


!    USE io_files,             ONLY : find_free_unit, prefix
    USE io_files,             ONLY :  prefix,tmp_dir
    USE io_global,            ONLY : stdout
    USE mp_world,  ONLY : mpime
    implicit none
    INTEGER, EXTERNAL :: find_free_unit
    TYPE(exc) :: a!the excitonic wave function to be read
    INTEGER :: label! the label identifing the required excitonic wavefunction
    LOGICAL, INTENT(in) :: l_verbose
    INTEGER :: iw, jw, iuna
    CHARACTER(5) :: nfile
    CHARACTER(5) :: nproc



!first deallocate
    call free_memory_exc_a(a)


    if(label >= 0 ) then
      write(nfile,'(5i1)') label/10000,mod(label,10000)/1000,mod(label,1000)/100,mod(label,100)/10,mod(label,10)
      write(nproc,'(5i1)') mpime/10000,mod(mpime,10000)/1000,mod(mpime,1000)/100,mod(mpime,100)/10,mod(mpime,10)
      iuna = find_free_unit()
      open( unit=iuna, file=trim(tmp_dir)//trim(prefix)//'-exc_a.'//nfile//'.'//nproc, status='old',form='unformatted')
    else
      write(nfile,'(5i1)') -label/10000,mod(-label,10000)/1000,mod(-label,1000)/100,mod(-label,100)/10,mod(-label,10)
      write(nproc,'(5i1)') mpime/10000,mod(mpime,10000)/1000,mod(mpime,1000)/100,mod(mpime,100)/10,mod(mpime,10)
      iuna = find_free_unit()
      open( unit=iuna, file=trim(tmp_dir)//trim(prefix)//'-exc_a.-'//nfile//'.'//nproc, status='old',form='unformatted')
    endif
    read(iuna) a%label
    read(iuna) a%npw
    read(iuna) a%numb_v
    read(iuna) a%e

!now allocate
    allocate(a%a(a%npw,a%numb_v))
    do iw=1,a%numb_v
       read(iuna)  a%a(1:a%npw,iw)
    enddo
    close(iuna)

    return
    END SUBROUTINE

    SUBROUTINE c_times_exc(a, c )
! this subroutine multiplies each iv line of the excitonic wave-function
! matrix (a%a) with the iv element of the real number vector (c).

    USE kinds, ONLY : DP

    implicit none
    type(exc) :: a
    real(kind=dp) :: c(a%numb_v)

    integer iv
    
    call start_clock('c_times_exc')

    do iv=1, a%numb_v
       a%a(1:a%npw,iv)=cmplx(c(iv),0.d0)*a%a(1:a%npw,iv)  
    enddo

    call stop_clock('c_times_exc')

    return
    END SUBROUTINE

    SUBROUTINE random_exc(a)
! this subroutine returns a random excitonic-wavefunction vector (a%a) 

    USE random_numbers, ONLY : randy
    USE kinds, ONLY : DP
    USE gvect,          ONLY : gstart

    implicit none
    type(exc) ::a

    real(kind=DP):: r1,r2
    integer      ::  iv,ig

    do iv=1,a%numb_v
       do ig=1,a%npw
          r1=randy()
          r2=randy()
          a%a(ig,iv)=cmplx(r1,r2)
          if (gstart==2) a%a(1,iv)=cmplx(r1,0.d0)
       enddo
    enddo  
    return
    END SUBROUTINE


    SUBROUTINE pc_operator_exc(a,v,is)
! this subroutine projects the excitonic wave-function vector into the conduction states manifold    
    USE bse_basic_structures
    USE mp, ONLY : mp_sum
    USE mp_world,             ONLY : world_comm
    USE wvfct,    ONLY : npw,npwx

    implicit none
    type(exc) :: a ! excitonic wfns vector to be projected on the conduction state manifold 
    type(v_state) :: v ! valence states vector in G space  

    REAL(kind=DP), ALLOCATABLE :: prod(:)
    integer :: iv, iiv,is

    call start_clock('pc_operator_exc')

    allocate(prod(a%numb_v))

    do iv=1,a%numb_v
       call dgemm('T','N', a%numb_v,1,2*a%npw,2.d0,v%wfn(:,:,is),2*npw,a%a(:,iv),2*a%npw,&
            & 0.d0,prod,a%numb_v)
        do iiv=1,a%numb_v
           if(v%gstart==2) prod(iiv)=prod(iiv)-dble(conjg(v%wfn(1,iiv,is))*a%a(1,iv))
        enddo
        call mp_sum(prod(:), world_comm )
        call dgemm('N','N',2*a%npw,1,a%numb_v,-1.d0,v%wfn(:,:,is),2*npw,prod,&
          &a%numb_v,1.d0,a%a(:,iv),2*a%npw)
    enddo

    deallocate(prod)

    call stop_clock('pc_operator_exc')

    return
    END SUBROUTINE

    
    SUBROUTINE poutcstate_exc(a_in,a_out,cstate,wcstate)
! this subroutine projects out from the excitonic vectors the component along the conduction states for which the 
! the QP are known, and it multiplies it by a weighted vector (used to avoid scissor)
    USE bse_basic_structures
!    USE qpe, ONLY : qpc
    USE mp, ONLY : mp_sum
    USE mp_world,             ONLY : world_comm
    USE wvfct,    ONLY : npw,npwx
   

    implicit none
    type(exc), intent(in) :: a_in ! excitonic wfns vector to be projected on the conduction state manifold 
    type(exc), intent(out) :: a_out ! final
    
    type(c_state) :: cstate ! GW corrected conduction states vector in G space  
    type(c_state) :: wcstate ! weighted GW corrected conduction states vector in G space  

    REAL(kind=DP), ALLOCATABLE :: prod(:)
    integer :: iv, ic

    call start_clock('poutcstate_exc')

    allocate(prod(cstate%numb_c))

    do iv=1,a_in%numb_v
       call dgemm('T','N', cstate%numb_c,1,2*a_in%npw,2.d0,cstate%wfn(:,:),2*npw,a_in%a(:,iv),2*a_in%npw,&
            & 0.d0,prod,cstate%numb_c)
!       call dgemm('T','N', a%numb_v,1,2*a%npw,2.d0,v%wfn(:,:,is),2*npwx,a%a(:,iv),2*a%npw,&
!            & 0.d0,prod,a%numb_v)
        do ic=1,cstate%numb_c
           if(cstate%gstart==2) prod(ic)=prod(ic)-dble(conjg(cstate%wfn(1,ic))*a_in%a(1,iv))
        enddo
        call mp_sum(prod(:), world_comm )
        call dgemm('N','N',2*a_in%npw,1,cstate%numb_c,1.d0,wcstate%wfn(:,:),2*npw,prod,&
          &cstate%numb_c,0.d0,a_out%a(:,iv),2*a_in%npw)
    enddo

    deallocate(prod)

    call stop_clock('poutcstate_exc')

    return
    END SUBROUTINE


    SUBROUTINE sproduct_exc(a1,a2,prod)
! this subroutine returns the scalar product (prod) between two
! excitonic-wavefunctions, input wave functions are given in G-space
    use io_global, ONLY : stdout, ionode 
    USE kinds, ONLY : DP
    USE mp,             ONLY : mp_sum, mp_barrier
    USE mp_world,             ONLY : world_comm
    use mp_world, ONLY : mpime
    USE gvect,          ONLY : gstart,ngm_g

    implicit none
    REAL(kind=DP), EXTERNAL :: ddot
    type(exc) :: a1,a2
    real(kind=DP) :: prod
    integer :: ii

    logical :: debug
    debug=.false.

!compute the dot product
!    if(debug) then
!      write(*,*) 'sproduct_exc in, mpime=',mpime
!    endif
    call start_clock('sproduct_exc')
    prod=0.d0
    do ii=1,a1%numb_v
       prod=prod+2.d0*ddot(2*a1%npw,a1%a(:,ii),1,a2%a(:,ii),1)
       if (gstart==2) prod=prod-a1%a(1,ii)*a2%a(1,ii)
    enddo
! sum over processor
    call mp_sum(prod, world_comm)



!    if(debug) then
!      write(*,*) 'sproduct_exc out, mpime=',mpime
!      write(*,*) 'sproduct_exc out, ngm_g= ',ngm_g
!    endif

    call stop_clock('sproduct_exc')
    return
    END SUBROUTINE


    SUBROUTINE normalize_exc(a)
! this subroutine normalizes the excitonic wave function to 1   
    USE kinds, ONLY : DP
    use io_global, ONLY : stdout, ionode 
    use mp_world, ONLY : mpime


    implicit none
    type(exc) :: a
    real(kind=DP) :: prod

    logical :: debug
    debug=.false. 

    call start_clock('normalize_exc')

    if(debug) then
      write(*,*) 'normalize_exc in, mpime=',mpime
    endif

    call sproduct_exc(a,a,prod)

    prod=1/sqrt(prod)

    a%a(1:a%npw,1:a%numb_v)= a%a(1:a%npw,1:a%numb_v)*prod
 
    if(debug) then
    ! check normalization'
       call sproduct_exc(a,a,prod)
       if(ionode) write(stdout,*) 'normalize exc check, prod=',prod
    endif 

    if(debug) then
      write(*,*) 'normalize_exc out, mpime=',mpime
    endif
    
    call stop_clock('normalize_exc')

    return
    END SUBROUTINE
    
    SUBROUTINE pout_operator_exc(a,i_state)
!   this subroutine projects out a component |b> from an excitonic
!   wavevector |a> 

    implicit none
    type(exc) :: a
    integer   :: i_state, i
    real(kind=DP), allocatable :: prod(:)

    call start_clock('pout_operator_exc')
    allocate(prod(i_state-1))
    
    do i=1,(i_state-1) 
       call sproduct_exc(a,bse_spectrum(i),prod(i))
    enddo
    
    do i=1,(i_state-1)
       a%a(1:a%npw,1:a%numb_v)=a%a(1:a%npw,1:a%numb_v)-prod(i)*&
                      bse_spectrum(i)%a(1:bse_spectrum(i)%npw,1:bse_spectrum(i)%numb_v)
    enddo

!    call normalize_exc(a)
    deallocate(prod)
    call stop_clock('pout_operator_exc')
    return
    END SUBROUTINE

    SUBROUTINE fft_a_exc(ag,fc,ar)
!   this subroutine performs an FFT to real space of the excitonic wavefunction
!   vector using the dual grid
    USE kinds, ONLY : DP
    USE gvect,                 ONLY : ig_l2g
    USE fft_custom_gwl
    USE bse_wannier, ONLY : dual_bse 
    USE wvfct,    ONLY : npwx
    USE gvecw,              ONLY : gcutw, ecutwfc
    USE io_global, ONLY : stdout, ionode, ionode_id
    USE mp_world, ONLY : mpime, nproc
    USE mp_pools, ONLY : intra_pool_comm
!    USE mp_wave, ONLY : mergewf,splitwf
    USE wavefunctions_module, ONLY :  psic


    implicit none
 
    type(exc) ag
    type(exc_r) ar
    type(fft_cus) :: fc

    COMPLEX(kind=DP), allocatable :: ag_t(:,:)
    COMPLEX(kind=DP), ALLOCATABLE :: evc_g(:)

    integer :: ii

! FFT the wannier function to r-space (dual grid)
    
   call start_clock('fft_a_exc')

    allocate(ag_t(fc%npwt,ag%numb_v))
    
    ar%nrxxt=fc%nrxxt 
    ar%numb_v=ag%numb_v
    ar%label=ag%label

    allocate(ar%ar(ar%nrxxt,ar%numb_v))
    
    allocate(evc_g(fc%ngmt_g )) 

    if(fc%dual_t==4.d0) then
       ag_t(1:fc%npwt,1:ag%numb_v)= ag%a(1:fc%npwt,1:ag%numb_v)
    else
       call reorderwfp_col(ag%numb_v,ag%npw,fc%npwt,ag%a(1,1),ag_t(1,1),ag%npw ,fc%npwt, &
           & ig_l2g,fc%ig_l2gt,fc%ngmt_g,mpime, nproc,intra_pool_comm )

!       do ii=1,ag%numb_v
!          call mergewf(ag%a(:,ii),evc_g,ag%npw,ig_l2g,mpime,nproc,ionode_id,intra_pool_comm)
!          call splitwf(ag_t(:,ii),evc_g,fc%npwt,fc%ig_l2gt,mpime,nproc,ionode_id,intra_pool_comm)
!       enddo
    endif

    do ii=1,ag%numb_v,2
       psic(:)=(0.d0,0.d0)
       if (ii==ag%numb_v) then
          psic(fc%nlt(1:fc%npwt))  = ag_t(1:fc%npwt,ii)
          psic(fc%nltm(1:fc%npwt)) = CONJG( ag_t(1:fc%npwt,ii) )
       else
          psic(fc%nlt(1:fc%npwt))=ag_t(1:fc%npwt,ii)+(0.d0,1.d0)*ag_t(1:fc%npwt,ii+1)
          psic(fc%nltm(1:fc%npwt))=CONJG(ag_t(1:fc%npwt,ii))+(0.d0,1.d0)*CONJG(ag_t(1:fc%npwt,ii+1))
       endif
       CALL cft3t( fc, psic, fc%nr1t, fc%nr2t, fc%nr3t, fc%nrx1t, fc%nrx2t, fc%nrx3t, 2 )
       ar%ar(1:fc%nrxxt,ii)= DBLE(psic(1:fc%nrxxt))
       if(ii/=ag%numb_v) ar%ar(1:fc%nrxxt,ii+1)= DIMAG(psic(1:fc%nrxxt))
    enddo

    deallocate(evc_g)
    call stop_clock('fft_a_exc')
    return
    END SUBROUTINE 

    SUBROUTINE fftback_a_exc(ar,fc,ag)
!   this subroutine performs an FFT G space of the excitonic wavefunction
!   from the dual R-grid and reorders the wavefunction in the pw grid order 
    USE kinds, ONLY : DP
    USE gvect,                 ONLY : ig_l2g
    USE fft_custom_gwl
    USE bse_wannier, ONLY : dual_bse 
    USE wvfct,    ONLY : npwx
    USE gvecw,              ONLY : gcutw, ecutwfc
    USE io_global, ONLY : stdout, ionode, ionode_id
    USE mp_world, ONLY : mpime, nproc
    USE mp_pools, ONLY : intra_pool_comm
!    USE mp_wave, ONLY : mergewf,splitwf
    USE wavefunctions_module, ONLY :  psic


    implicit none
 
    type(exc) ag
    type(exc_r) ar
    type(fft_cus) :: fc
    integer :: ii,iv
    COMPLEX(kind=DP), ALLOCATABLE :: evc_g(:)
    COMPLEX(kind=DP), allocatable :: ag_t(:,:)


! FFT the wannier function to g-space ( from dual grid), and put back in the pw
! grid order 

    call start_clock('fftback_a_exc')
    allocate(ag_t(fc%npwt,ar%numb_v))
    allocate(evc_g(fc%ngmt_g ))


    do iv=1, ag%numb_v,2
       if (iv==ag%numb_v) then 
          psic(1:fc%nrxxt)=dcmplx(ar%ar(1:ar%nrxxt,iv),0.d0)
       else
          psic(1:fc%nrxxt)=dcmplx(ar%ar(1:fc%nrxxt,iv),ar%ar(1:fc%nrxxt,iv+1))
       endif
       CALL cft3t(fc, psic, fc%nr1t, fc%nr2t, fc%nr3t, fc%nrx1t, fc%nrx2t, fc%nrx3t, -2 )
       if (iv==ag%numb_v) then 
          ag_t(1:fc%npwt,iv)=psic(fc%nlt(1:fc%npwt))
       else
          ag_t(1:fc%npwt,iv)=0.5d0*(psic(fc%nlt(1:fc%npwt))+conjg( psic(fc%nltm(1:fc%npwt))))
          ag_t(1:fc%npwt,iv+1)=(0.d0,-0.5d0)*(psic(fc%nlt(1:fc%npwt)) - conjg(psic(fc%nltm(1:fc%npwt))))
       endif
    enddo
 
    if(fc%dual_t==4.d0) then
       ag%a(1:fc%npwt,1:ag%numb_v)=ag_t(1:fc%npwt,1:ag%numb_v)
    else
        call reorderwfp_col(ag%numb_v,fc%npwt,ag%npw,ag_t(1,1),ag%a(1,1),fc%npwt,ag%npw, &
           & fc%ig_l2gt,ig_l2g,fc%ngmt_g,mpime, nproc,intra_pool_comm )

!       do iv=1, ag%numb_v
!          call mergewf(ag_t(:,iv),evc_g,fc%npwt,fc%ig_l2gt,mpime,nproc,ionode_id,intra_pool_comm)
!          call splitwf(ag%a(:,iv),evc_g,ag%npw,ig_l2g,mpime,nproc,ionode_id,intra_pool_comm)
!       enddo
    endif 

    deallocate(evc_g)
    call stop_clock('fftback_a_exc')
    return
  END SUBROUTINE fftback_a_exc

    SUBROUTINE urot_a(a_in,a_out,itrasp)
    USE wvfct,     ONLY : nbnd,npwx
    USE bse_basic_structures, ONLY : u_trans
    USE lsda_mod,             ONLY : nspin
    USE io_global, ONLY : stdout 



    implicit none 
    type(exc):: a_in
    type(exc):: a_out
    integer  :: itrasp ! if 1 takes U^T
    REAL(kind=DP), ALLOCATABLE :: tmp_rot(:,:)
    logical debug
    integer :: ii


    call start_clock('urot_a')

    debug=.false.
    allocate(u_trans(nbnd,nbnd,nspin))
    call read_wannier_matrix


    allocate(tmp_rot(a_in%numb_v,a_in%numb_v))
    tmp_rot(1:a_in%numb_v,1:a_in%numb_v)=dble(u_trans(1:a_in%numb_v,1:a_in%numb_v,1))

!DEBUG 
!    tmp_rot=0.d0
!    do ii=1, a_in%numb_v
!       tmp_rot(ii,ii)=1.d0
!    enddo
!fine DEBUG

   
    if (itrasp==0) call rotate_wannier_gamma_bse(tmp_rot,a_in,a_out,1,0)
    if (itrasp==1) call rotate_wannier_gamma_bse(tmp_rot,a_in,a_out,1,1)

    deallocate(u_trans)
    deallocate(tmp_rot)

    call stop_clock('urot_a')
    return
    END SUBROUTINE 
     

end module  exciton
