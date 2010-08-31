! Author: P. Umari
! Modified by G. Stenuit
!
 SUBROUTINE create_hf(options, qp)
!this subroutine creates the perturbative HF energies
!and allocates and sets relevant stuff

   USE basic_structures, ONLY : wannier_u, free_memory
   USE input_gw, ONLY : input_options
   USE constants, ONLY : RYTOEV
   USE kinds, ONLY : DP
   USE energies_gww, ONLY : quasi_particles
   USE io_global, ONLY : stdout

   implicit none

   TYPE(input_options) :: options! for prefix
   TYPE(quasi_particles) :: qp!the descriptor to be build
   TYPE(wannier_u) :: uu
   REAL(kind=DP), ALLOCATABLE :: ene_x(:), ene_h(:)
   INTEGER :: ii


   call read_data_pw_u(uu, options%prefix)


!allocates

   qp%max_i=options%max_i
   allocate(qp%ene_dft_ks(qp%max_i))
   allocate(qp%ene_dft_xc(qp%max_i))
   allocate(qp%ene_dft_h(qp%max_i))
   allocate(qp%ene_gw(qp%max_i))
   allocate(qp%ene_gw_pert(qp%max_i))
   allocate(qp%ene_hf(uu%nums))
   if(options%l_hf_energies) then
      allocate(qp%ene_x(uu%nums))
      allocate(ene_x(uu%nums))
      allocate(ene_h(uu%nums))
   else
      allocate(qp%ene_x(qp%max_i))
      allocate(ene_x(qp%max_i))
      allocate(ene_h(qp%max_i))
   endif
   allocate(qp%ene_h(qp%max_i))

   qp%ene_dft_ks(1:qp%max_i) =  uu%ene(1:qp%max_i)
   qp%ene_dft_xc(1:qp%max_i) =  uu%ene_xc(1:qp%max_i)
   qp%ene_dft_h(1:qp%max_i) =  uu%ene_lda_h(1:qp%max_i)

   if(options%l_lda_hartree) then
       qp%ene_h(1:qp%max_i) =  cmplx(uu%ene_lda_h(1:qp%max_i),0.d0)
    else
!here calculate hartree parte
    endif

!calculate exchange part
   if(options%l_hf_energies) then
      call go_exchange(options,ene_x,ene_h,uu%nums)
      do ii=1,uu%nums
         qp%ene_x(ii)=ene_x(ii)
         qp%ene_hf(ii)=uu%ene(ii)-uu%ene_xc(ii)+ene_x(ii)
          write(stdout,*) 'ENE HF',ii, qp%ene_hf(ii)*RYTOEV
      enddo
   else
      write(stdout,*) 'ENE H', qp%ene_h(1:qp%max_i)
      if(.not.options%l_lda_exchange) then
         call go_exchange(options,ene_x,ene_h,qp%max_i)
      else
!read from file
         write(stdout,*) 'read ene_x from files'
         call read_data_pw_exchange(ene_x,qp%max_i,options%prefix)
         write(stdout,*) 'DONE: read ene_x from files'
         ene_h(1:qp%max_i)=uu%ene_lda_h(1:qp%max_i)
         write(stdout,*) 'done'
      endif
      do ii=1,qp%max_i
         qp%ene_x(ii)=ene_x(ii)
         if(options%l_lda_hartree) then
            qp%ene_hf(ii)=uu%ene(ii)-uu%ene_xc(ii)+ene_x(ii)
         else
            qp%ene_h(ii)=cmplx(ene_h(ii),0.d0)
            qp%ene_hf(ii)=uu%ene(ii)-uu%ene_xc(ii)+ene_x(ii)-uu%ene_lda_h(ii)+ene_h(ii)
         endif
      enddo
   endif

   call free_memory(uu)
   deallocate(ene_x,ene_h)
   return

 END SUBROUTINE create_hf
