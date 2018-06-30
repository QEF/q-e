subroutine diago_exc_cg(data_input, bd, pp, pt, pm, x)

  USE input_simple_exc
  USE simple_objects
  USE derived_objects
  USE kinds, ONLY : DP
  USE io_global, ONLY : stdout, ionode
  USE constants, ONLY : rytoev
  USE mp, ONLY : mp_barrier
  USE mp_world,             ONLY : world_comm
  
  implicit none

  TYPE(input_options), INTENT(in) :: data_input
  TYPE(bands), INTENT(in) :: bd
  TYPE(prod_proj), INTENT(in) :: pp
  TYPE(potential), INTENT(in) :: pt
  TYPE(prod_mix), INTENT(in) :: pm
  TYPE(exc) :: xd, x2d, Hxd, Hx2d, h, Hx, Hx1, Hh, xtest
  TYPE(exc) :: x(data_input%nvec)
  INTEGER :: i, j, k, count=1, beginning, rate, end,kf
  COMPLEX(kind=DP) :: d=(0.1d0,0.0d0), a, b, c, fd, f2d, alpha,stim, minus=(-1.0d0,0.0d0), echeck
  COMPLEX(kind=DP), allocatable :: energy(:)
  LOGICAL :: test
  
  TYPE(exc) :: etmp1,etmp2,etmp3,etmp4
  REAL(kind=DP) :: ene0,dene0,ene1,ene2,passop,passo,stima,enew,gamma,esse,essenew,ene
  COMPLEX(kind=DP) :: passot=(1.d0,0.d0),passos,sca, passot2=(0.4d0,0.d0)
  LOGICAL :: l_ok
  REAL(kind=DP) :: norm
  INTEGER :: icon


  call setup_exc(bd,h)
  call setup_exc(bd,Hx)
  call setup_exc(bd,Hx1)
  call setup_exc(bd,Hh)
  call setup_exc(bd,xtest)

  call setup_exc(bd,etmp1)
  call setup_exc(bd,etmp2)
  call setup_exc(bd,etmp3)
  call setup_exc(bd,etmp4)



  call mp_barrier(world_comm)
  allocate(energy(data_input%nvec))
  write(stdout,*)data_input%nvec
  vectors: do j=1,data_input%nvec
     icon=0
     call randomize_exc(x(j))
     call normalize_exc(x(j))
     !x(j)%avc=0.d0
     !x(j)%avc(1,1,1)=(1.d0,0.d0)
     
     call hamiltonian(data_input, 1, bd, pp, pt, pm, x(j), Hx,0)
     call mp_barrier(world_comm)
     energy(j)=x(j)*Hx*rytoev!energy initialization
     h=(-1.d0,0.d0)*Hx!search direction initialization
     !call normalize_exc(h)
     esse=dble(h*h)
!     h=minus*h
     l_ok=.true.
     iterate:do i=1,data_input%max_nstep
            echeck=energy(j)
           call mp_barrier(world_comm)
           call hamiltonian(data_input, 1, bd, pp, pt, pm, h, Hh,0)
           
           ene0=dble(x(j)*Hx)
           dene0=2.d0*dble(x(j)*Hh)-2.d0*ene0*dble(x(j)*h)
           norm=sqrt(dble(h*h))
           if(.true.) then
              if(dene0>0.d0) then
                 !dene0=-dene0
                 passot=-(1.d0,0.d0)
                 passot2=-(2.d0,0.d0)
              else
                 passot=(1.d0,0.d0)
                 passot2=(2.,0.d0)
              endif
           else
              if(dene0>0.d0) then
                 !dene0=-dene0                                                                                                               
                 passot=-(0.5d0,0.d0)/norm
                 passot2=-(1.d0,0.d0)/norm
              else
                 passot=(0.5d0,0.d0)/norm
                 passot2=(1.,0.d0)/norm
              endif
           endif
           etmp1=passot*h
           call sum_exc_sub(etmp2,x(j),etmp1)
           call normalize_exc(etmp2)
           call hamiltonian(data_input, 1, bd, pp, pt, pm, etmp2, etmp3,0)
           ene1=dble(etmp2*etmp3)
        
           etmp1=passot2*h
           call sum_exc_sub(etmp2,x(j),etmp1)
           call normalize_exc(etmp2)
           call hamiltonian(data_input, 1, bd, pp, pt, pm, etmp2, etmp3,0)
           ene2=dble(etmp2*etmp3)

           !call minparabola(ene0,dene0,ene1,dble(passot),passo,stima,l_ok)
          call minparabola2(ene0,ene1,ene2,dble(passot),dble(passot2),passo,stima,l_ok)  
          ! do k=-1000,1000
          !    passos=cmplx(dble(k)/20.d0,0.d0)
          !    etmp1=passos*h
          !    call sum_exc_sub(etmp2,x(j),etmp1)
          !    call normalize_exc(etmp2)
          !    call hamiltonian(data_input, 1, bd, pp, pt, pm, etmp2, etmp3,0)
          !    ene=dble(etmp2*etmp3)
          !    write(stdout,*) dble(k)/20.d0, ene*rytoev
          ! enddo


           if(l_ok) then
              passos=cmplx(passo,0.d0)
              etmp1=passos*h
              call sum_exc_sub(xtest,x(j),etmp1)
           else
              !write(stdout,*) l_ok,norm
              h=minus*Hx
              !ene2=ene0
              !kf=0
              !do k=1,40
              !   sca=cmplx(dble(k)/10.d0,0.d0)
              !   etmp1=sca*h
              !   call sum_exc_sub(xtest,x(j),etmp1)
              !   call normalize_exc(xtest)
              !   call hamiltonian(data_input, 1, bd, pp, pt, pm, xtest, Hx1,0)
              !   enew=dble(xtest*Hx1)
              !   
              !   if(enew<ene2) then
              !      ene2=enew
              !      kf=k
              !   endif
              !enddo
              !write(stdout,*) 'KF:', kf
              !sca=cmplx(dble(kf)/10.d0,0.d0)
              sca=(2.d0,0.d0)
              etmp1=sca*h
              call sum_exc_sub(xtest,x(j),etmp1)
           endif

           call normalize_exc(xtest)
          
           call hamiltonian(data_input, 1, bd, pp, pt, pm, xtest, Hx1,0)
           call mp_barrier(world_comm)
              !     write(stdout,*)i,'check:',real(echeck),real(xtest*Hx1*rytoev)
           enew=dble(xtest*Hx1)
           ! write(stdout,*) 'LIN MIN', dble(echeck),stima*rytoev,enew*rytoev,passo
           write(stdout,*) i,enew*rytoev
        
           if (dble(echeck)-enew*rytoev<0.d0 .or. mod(i,20)==0 .or. .not.l_ok) then
              !write(stdout,*) 'NO GOOD'

              
              h=minus*Hx
              sca=cmplx(2.d0,0.d0)
              etmp1=sca*h
              call sum_exc_sub(xtest,x(j),etmp1)
              

              x(j) = xtest
              call normalize_exc(x(j))
         
              call mp_barrier(world_comm)
              call hamiltonian(data_input, 1, bd, pp, pt, pm, x(j), Hx,0)
              call mp_barrier(world_comm)
              h=Hx
              esse=dble(h*h)
             ! write(stdout,*)j,i,'sd',real(echeck),real(x(j)*Hx*rytoev)
              !write(stdout,*)j,i,real(x(j)*Hx*rytoev)
             
              energy(j)=x(j)*Hx*rytoev
              
           else
              
              x(j)=xtest
              Hx=Hx1
              gamma=dble(Hx*Hx)

              essenew=gamma
              gamma=gamma/esse
              esse=essenew
              sca=cmplx(gamma,0.d0)
              etmp1=sca*h
              etmp2=(1.d0,0.d0)*Hx
              call sum_exc_sub(h,etmp2,etmp1)
             


              !write(stdout,*)j,i,'cg', real(echeck),real(x(j)*Hx*rytoev)
             ! write(stdout,*)j,i,real(x(j)*Hx*rytoev)
           
              energy(j)=x(j)*Hx*rytoev
              if (i==data_input%max_nstep .or. (real(-energy(j)+echeck)<data_input%thr_evc .and. i/=1)) then
                 write(stdout,*)j,'converged'
                 exit 
              end if
           end if
           if(abs(-energy(j)+echeck)<data_input%thr_evc) then
              icon=icon+1
           else
              icon=0
           endif
           if(icon==10) exit
           

     end do iterate
     x(j)%ene(1)=dble(energy(j))
     call hamiltonian(data_input, 1, bd, pp, pt, pm, x(j), Hx,1)
     x(j)%ene(2)=dble(x(j)*Hx*rytoev)
     call hamiltonian(data_input, 1, bd, pp, pt, pm, x(j), Hx,2)
     x(j)%ene(3)=dble(x(j)*Hx*rytoev)
     call hamiltonian(data_input, 1, bd, pp, pt, pm, x(j), Hx,3)
     x(j)%ene(4)=dble(x(j)*Hx*rytoev)



     energy(j)=energy(j)/rytoev!Ry   
  end do vectors
  write(stdout,*) 'calling routine spectrum'
  !spectrum
  !call spectrum(data_input,x,bd,energy)
  !lanczos
  deallocate(energy)
  call deallocate_exc(Hx)
  call deallocate_exc(h)
  call deallocate_exc(Hx1)
  call deallocate_exc(Hh)
  call deallocate_exc(xtest)
  
  call deallocate_exc(etmp1)
  call deallocate_exc(etmp2)
  call deallocate_exc(etmp3)
  call deallocate_exc(etmp4)

  return
  
end subroutine diago_exc_cg



 subroutine minparabola(ene0,dene0,ene1,passop,passo,stima,l_ok)
!this subroutines finds the minimum of a quadratic real function                                                                                          \
                                                                                                                                                           

      use kinds, only : dp
      use io_global, only :stdout

      implicit none
      real(dp) ene0,dene0,ene1,passop,passo,stima
      real(dp) a,b,c!a*x^2+b*x+c                               
      logical :: l_ok
                                                                                                                                                           

      c=ene0
      b=dene0
      a=(ene1-b*passop-c)/(passop**2.d0)

      passo = -b/(2.d0*a)
      if( a.lt.0.d0) then
         l_ok=.false.
       !  write(stdout,*) 'CAZZI'
         if(ene1.lt.ene0) then
            passo=passop
         else
            passo=0.5d0*passop
         endif
      else
         l_ok=.true.
      endif


      stima=a*passo**2.d0+b*passo+c


      return
    end subroutine minparabola

    subroutine minparabola2(ene0,ene1,ene2,x1,x2,x,stima,l_ok)
      
      use kinds, only : DP
      use io_global, only : stdout

      implicit none

      real(kind=dp) :: ene0,ene1,ene2,x1,x2,x,stima
      logical :: l_ok

      real(kind=dp) :: a,b,c
      c=ene0
      a=(ene2-ene1*x2/x1+c*x2/x1-c)/(x2**2.d0-x1*x2)
      b=(ene1-c-a*x1**2.d0)/x1


      x = -b/(2.d0*a)
      if( a.lt.0.d0) then
        ! write(stdout,*) 'CAZZI'
         l_ok=.false.
         if(ene1<=ene0) then
            x=x1
         else
            x=0.5d0*x1
         endif
      else
         l_ok=.true.
      endif


      stima=a*x**2.d0+b*x+c



      
      return
    end subroutine minparabola2
