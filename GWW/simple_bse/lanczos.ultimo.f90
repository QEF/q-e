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
  TYPE(exc):: x,Hx,xm1,x1
  TYPE(epe):: element
  TYPE(bands):: bd
  TYPE(product) :: pd
  TYPE(prod_proj) :: pp
  TYPE(potential) :: pt
  TYPE(prod_mix) :: pm
  COMPLEX(KIND=DP) :: fraz,zeta
  REAL(KIND=DP) :: omin,omax,step, cost, omega
  COMPLEX(KIND=DP),ALLOCATABLE:: cabs(:,:),eemat(:,:,:,:),a(:),b(:)
  INTEGER :: i,dir,v,c,k,io,iun,inizio,fine,count_rate,idir
  INTEGER, EXTERNAL :: find_free_unit
  REAL(KIND=DP),ALLOCATABLE::diel(:,:)
  REAL(kind=DP) :: norm(3)

  call initialize_product(pd)
  call initialize_potential(pt)
  call initialize_prod_mix(pm)
  call initialize_prod_proj(pp)
  call system_clock(inizio,count_rate) 
  call read_bands(data_input,bd)
  call read_product(data_input,pd)
  call read_potential(data_input,pt)

  write(stdout ,*)'building derived objects'

  call build_prod_mix(data_input,bd,pd,pm,pt)
  call build_prod_proj(bd,pd,pp)

  call setup_exc(bd,x)
  call setup_exc(bd,xm1)
  call setup_exc(bd,Hx)
  call setup_exc(bd,x1)
  
  data_input%lanczos_step=data_input%lanczos_step+1
  
  allocate(eemat(bd%numv,bd%numc,bd%nk_loc,3))
  allocate(diel(data_input%spectrum_points,3))
  allocate(cabs(data_input%spectrum_points,3))
  allocate(a(data_input%lanczos_step))
  allocate(b(data_input%lanczos_step-1))
  
  !haydoch recursive method
!  call system_clock(inizio,count_rate)
!  do dir=1,3
     dir=1
     write(stdout,*)'Dir',dir,'constructing starting vector'
     !constructing starting vector
     call build_eemat(data_input,bd,eemat)
     write(stdout,*)'eemat done'
     if ( x%nk_loc >  0 ) then
        do k=1,bd%nk_loc
           do c=1,bd%numc
              do v=1,bd%numv
                 x%avc(v,c,k) = eemat(v,c,k,dir)/(bd%en_v(v,k)-bd%en_c(c,k))
                 !x%avc(v,c,k) = conjg(eemat(v,c,k,dir))/(bd%en_v(v,k)-bd%en_c(c,k))
              end do
           end do
        end do
     end if
     norm(dir)=x*x
     write(stdout,*) 'Normalizzazione', norm(dir)
     call normalize_exc(x)
     !iterations
     write(*,*)'Dir',dir,'matrix construction'
     do i =1,data_input%lanczos_step
        call mp_barrier(world_comm)
        call hamiltonian(data_input, 1, bd, pp, pt, pm, x, Hx)
        call mp_barrier(world_comm)
        a(i) = cmplx(dble(x*Hx),0.d0)
 !       write(*,*)dir,i,'a',a(i)
        if (i==1) then
           b(i) =cmplx(dble( sqrt((Hx+(-1.d0,0.d0)*(a(i)*x))*(Hx+(-1.d0,0.d0)*(a(i)*x)))),0.d0)
 !          write(*,*)dir,i,'b',b(i)
           xm1 = x
           x =(1/b(i))*(Hx+(-1.d0,0.d0)*(a(i)*x))
        else if (i ==data_input%lanczos_step) then
           exit
        else
           b(i) =cmplx(dble(  sqrt((Hx+(-1.d0,0.d0)*(a(i)*x+b(i-1)*xm1))*(Hx+(-1.d0,0.d0)*(a(i)*x+b(i-1)*xm1)))),0.d0)
!           write(*,*)dir,i,'b',b(i)
           x1 = (1/b(i))*(Hx+(-1.d0,0.d0)*(a(i)*x+b(i-1)*xm1))! (|i+1>= H|i>-a_i|i>-b_{i-1}|i-1>)/b_i
           xm1 = x!|i-1> = |i>
           x = x1!|i> = |i+1>
        end if
        write(*,*)dir,i
     end do
     !construction of the continued fraction
     omin = data_input%omega_min/rytoev!Ry
     omax = data_input%omega_max/rytoev!Ry
     step =  (omax-omin)/(data_input%spectrum_points-1)
     cost = - 8/3.14159/10.26!general case for V?
     write(*,*)'Dir',dir,'continued fraction'
     do io=1,data_input%spectrum_points
        fraz=(0.d0,0.d0)
        zeta =cmplx( omin + (io-1)*step, data_input%eta ) !Ry
        fraz = (b(data_input%lanczos_step-1)*b(data_input%lanczos_step-1))/(zeta-a(data_input%lanczos_step))
        do i=2,data_input%lanczos_step-1
           fraz = (b(data_input%lanczos_step-i)*b(data_input%lanczos_step-i))/&
                &(zeta - a(data_input%lanczos_step-i+1)-fraz)
           diel(io,dir)=cost*aimag(1/(zeta-a(1)-fraz))!imaginary part
           
        end do
     end do
  call system_clock(fine)
  write(stdout,*)'elapsed time',real(fine-inizio)/real(count_rate)
  !relation for the total dielectric function
  !saving spectrum
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
  if(ionode) then
     iun=find_free_unit()
     open( unit= iun, file='spectrum_lanczos.dat',status='unknown')
     do io=1,data_input%spectrum_points
        omega = omin + (io-1)*step
!        cabs(io,dir)=omega*aimag(diel(io))/(sqrt(4*3.14159*diel(io)))/(137)
        !        write(iun,*) omega*rytoev, cabs(io)
        write(iun,*) omega*rytoev, (diel(io,1)+diel(io,2)+diel(io,3))/3
     end do
     close(iun)
  end if
  
  call deallocate_bands(bd)
  call deallocate_product(pd)
  call deallocate_potential(pt)
  call deallocate_prod_mix(pm)
  call deallocate_prod_proj(pp)
  call deallocate_exc(x1)
  call deallocate_exc(x)
  call deallocate_exc(xm1)
  call deallocate_exc(Hx)
  deallocate(a)
  deallocate(b)
  deallocate(cabs)
  deallocate(diel)
  deallocate(eemat)
  return
  
end subroutine lanczos
