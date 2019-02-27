subroutine lanczos(data_input)
  USE input_simple_exc
  USE simple_objects
  USE derived_objects
  USE kinds, ONLY : DP
  USE io_global, ONLY : stdout, ionode
  USE constants, ONLY : rytoev
  USE mp, ONLY : mp_barrier, mp_sum
  USE mp_world,             ONLY : world_comm
  use io_files, only : prefix, tmp_dir
  USE constants, ONLY : rytoev

  implicit none

  TYPE(input_options):: data_input
  TYPE(exc):: x,Hx,xm1,x1,etmp1,etmp2,etmp3,etmp4
  TYPE(epe):: element
  TYPE(bands):: bd
  TYPE(product) :: pd
  TYPE(prod_proj) :: pp
  TYPE(potential) :: pt
  TYPE(prod_mix) :: pm
  COMPLEX(KIND=DP) :: fraz,zeta
  REAL(KIND=DP) :: omin,omax,step, cost, omega
  COMPLEX(KIND=DP),ALLOCATABLE:: cabs(:,:),eemat(:,:,:,:),a(:),b(:)
  INTEGER :: i,dir,v,c,k,io,iun,inizio,fine,count_rate
  INTEGER, EXTERNAL :: find_free_unit
  REAL(KIND=DP),ALLOCATABLE::diel(:,:)
  COMPLEX(kind=DP) :: csca
  REAL(kind=DP) :: norm(3)
  INTEGER :: idir

  call initialize_product(pd)
  call initialize_potential(pt)
  call initialize_prod_mix(pm)
  call initialize_prod_proj(pp)
  call system_clock(inizio,count_rate) 
  call read_bands(data_input,bd)
  call read_product(data_input,pd)
  call read_potential(data_input,pt)
  write(stdout,*)'building derived objects mix'
  call build_prod_mix(data_input,bd,pd,pm,pt)
  write(stdout,*)'building derived objects  proj'
  call build_prod_proj(bd,pd,pp)

  write(stdout,*)'setup exc'
  call setup_exc(bd,x)
  call setup_exc(bd,xm1)
  call setup_exc(bd,Hx)
  call setup_exc(bd,x1)
  call setup_exc(bd,etmp1)
  call setup_exc(bd,etmp2)
  call setup_exc(bd,etmp3)
  call setup_exc(bd,etmp4)
  
  allocate(eemat(bd%numv,bd%numc,bd%nk_loc,3))
  allocate(diel(data_input%spectrum_points,3))
  allocate(cabs(data_input%spectrum_points,3))
  allocate(a(data_input%lanczos_step))
  allocate(b(data_input%lanczos_step-1))
  

  diel=0.d0
  !haydoch recursive method
!  call system_clock(inizio,count_rate)
  do dir=1,3
!     dir=1
     write(stdout,*)'Dir',dir,'constructing starting vector'
     !constructing starting vector
     call build_eemat(data_input,bd,eemat)
     write(stdout,*)'eemat done'
     if ( x%nk_loc >  0 ) then
        do k=1,bd%nk_loc
           do c=1,bd%numc
              do v=1,bd%numv
                 x%avc(v,c,k) = conjg(eemat(v,c,k,dir))/(bd%en_v(v,k)-bd%en_c(c,k))
                  end do
           end do
        end do
     end if
     norm(dir)=x*x
     write(stdout,*) 'Normalizzazione', norm(dir)
     call normalize_exc(x)

     !iterations
     write(stdout,*)'Dir',dir,'matrix construction'
     do i =1,data_input%lanczos_step
        write(stdout,*) 'Step ', i
        call hamiltonian(data_input, 1, bd, pp, pt, pm, x, Hx,0)
        csca=x*Hx
        a(i) = cmplx(dble(csca),0.d0)
        if (i==1) then
           etmp1=a(i)*(-1.d0,0.d0)*x
           call sum_exc_sub(etmp2,Hx,etmp1)
           !etmp2=Hx+etmp1
           csca=etmp2*etmp2
           b(i)=cmplx(sqrt(dble(csca)),0.d0)

         !  b(i) =cmplx(dble( sqrt((Hx+(-1.d0,0.d0)*(a(i)*x))*(Hx+(-1.d0,0.d0)*(a(i)*x)))),0.d0)
           xm1 = x
           etmp1=((-1.d0,0.d0)*a(i))*x
           call sum_exc_sub(etmp2,Hx,etmp1)
           x=(1/b(i))*etmp2
!           x =(1/b(i))*(Hx+(-1.d0,0.d0)*(a(i)*x) )
        else if (i ==data_input%lanczos_step) then
           exit
        else
           
           etmp1=(-1.d0,0.d0)*a(i)*x
           
           etmp2=(-1.d0,0.d0)*b(i-1)*xm1
           
           call sum_exc_sub(etmp3,etmp1,etmp2)

!           etmp3=etmp1+etmp2
           
           call sum_exc_sub(etmp4,Hx,etmp3)
!           etmp4=Hx+etmp3
           
           !b(i)=cmplx(dble(sqrt(etmp4*etmp4)),0.d0)
           b(i)=cmplx(sqrt(dble(etmp4*etmp4)),0.d0)
           

!           b(i) =cmplx(dble(  sqrt((Hx+(-1.d0,0.d0)*(a(i)*x+b(i-1)*xm1))*(Hx+(-1.d0,0.d0)*(a(i)*x+b(i-1)*xm1)))),0.d0)
!           write(*,*)dir,i,'b',b(i)
           etmp1=(-1.d0,0.d0)*a(i)*x
           etmp2=(-1.d0,0.d0)*b(i-1)*xm1
           call sum_exc_sub(etmp3,etmp1,etmp2)
           call sum_exc_sub(etmp4,Hx,etmp3)
           x1= (1/b(i))*etmp4
          ! x1 = (1/b(i))*(Hx+(-1.d0,0.d0)*(a(i)*x+b(i-1)*xm1))! (|i+1>= H|i>-a_i|i>-b_{i-1}|i-1>)/b_i
           
           xm1 = x!|i-1> = |i>
                      
           x = x1!|i> = |i+1>
           
        end if
        write(stdout,*)dir,i
     end do
     write(stdout,*)'writing matrix',dir
    
    
  end do
  call system_clock(fine)
  write(stdout,*)'elapsed time',real(fine-inizio)/real(count_rate)
  if(ionode) then
     iun=find_free_unit()
     open( unit= iun, file='ab_coeff.dat',status='unknown')
     write(iun,*) norm(1:3)
     do idir=1,3
        do i=1,data_input%lanczos_step-1
           write(iun,*)  idir, i, real(a(i))
        enddo
     enddo
     do idir=1,3
        do i=1,data_input%lanczos_step-1
           write(iun,*)  idir, i, real(b(i))
        enddo
     end do
     close(iun)
  endif

    
  call deallocate_bands(bd)
  call deallocate_product(pd)
  call deallocate_potential(pt)
  call deallocate_prod_mix(pm)
  call deallocate_prod_proj(pp)
  call deallocate_exc(x1)
  call deallocate_exc(x)
  call deallocate_exc(xm1)
  call deallocate_exc(Hx)
  call deallocate_exc(etmp1)
  call deallocate_exc(etmp2)
  call deallocate_exc(etmp3)
  call deallocate_exc(etmp4)
  
  deallocate(a)
  deallocate(b)
  deallocate(cabs)
  deallocate(diel)
  deallocate(eemat)
  return
  
end subroutine lanczos

