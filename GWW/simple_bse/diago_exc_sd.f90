SUBROUTINE diago_exc_sd(sin,bd,pp,pt,pm,a)
!this subroutine find the lowest (AT THE MOMENT)
!eigen values /vector of the excitonic Hamiltonia
!with the simple steepest descent scheme
  USE input_simple_exc
  USE simple_objects
  USE derived_objects
  USE kinds, ONLY : DP
  USE io_global, ONLY : stdout, ionode
  USE constants, ONLY : rytoev
  USE mp, ONLY : mp_barrier
  USE mp_world,             ONLY : world_comm
  
  
  implicit none
  
  TYPE(input_options), INTENT(in) :: sin
  TYPE(bands), INTENT(in) :: bd
  TYPE(prod_proj), INTENT(in) :: pp
  TYPE(potential), INTENT(in) :: pt
  TYPE(prod_mix), INTENT(in) :: pm
  TYPE(exc) :: a(sin%nvec)
  
  INTEGER :: ii,jj,kk,is
  TYPE(exc) :: ha,b
  COMPLEX(kind=DP) ::passop,ene

  call setup_exc(bd,ha)
  call setup_exc(bd,b)

  if(ionode) then
     write(stdout,*) 'Routine diago_exc_sd'
    
  endif

!now just one vector
!randoimize vectors 
  call mp_barrier(world_comm)
  write(stdout,*) 'ATT-1'
  do ii=1,sin%nvec
     call randomize_exc(a(ii))
  end do
  !orthonormalize them
  call mp_barrier(world_comm)
  write(stdout,*) 'ATT-2'
  do ii=1,sin%nvec
     call normalize_exc(a(ii))
  enddo
  passop=cmplx(sin%l_step,0)
!lopp
  do is=1,sin%max_nstep
!!find gradient
     call mp_barrier(world_comm)
     write(stdout,*) 'ATT6'   
     call hamiltonian(sin,1,bd,pp,pt,pm,a(1),ha,0)
     call mp_barrier(world_comm)
     write(stdout,*) 'ATT7'   
     ene=a(1)*ha
     if(ionode) then
        write(stdout,*) 'SD step energy :',is, ene*rytoev
     endif
     write(stdout,*) 'ATT1'
     b=(-passop)*ha
     write(stdout,*) 'ATT2'   
     ha=a(1)+b
     write(stdout,*) 'ATT3'   
     call normalize_exc(ha)
     write(stdout,*) 'ATT4'   
     a(1)=ha
     write(stdout,*) 'ATT5'   
!!find gradient
!!find new vector
!!normalize
!!calculate energy
!!print it
  end do


  call deallocate_exc(ha)
  call deallocate_exc(b)
  return

END SUBROUTINE diago_exc_sd
