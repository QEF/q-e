!
! Copyright (C) 2001-2013 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
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
   REAL(kind=DP), ALLOCATABLE :: ene_x(:,:), ene_h(:,:)
   INTEGER :: ii,is


   call read_data_pw_u(uu, options%prefix)


!allocates

   qp%max_i=options%max_i
   qp%nspin=uu%nspin
   qp%whole_s=.false.!if required whole matrix stuff is set elsewhere
   allocate(qp%ene_dft_ks(qp%max_i,qp%nspin))
   allocate(qp%ene_dft_xc(qp%max_i,qp%nspin))
   allocate(qp%ene_dft_h(qp%max_i,qp%nspin))
   allocate(qp%ene_gw(qp%max_i,qp%nspin))
   allocate(qp%ene_gw_pert(qp%max_i,qp%nspin))
   allocate(qp%ene_hf(uu%nums,qp%nspin))
   if(options%l_hf_energies) then
      allocate(qp%ene_x(uu%nums,qp%nspin))
      allocate(ene_x(uu%nums,qp%nspin))
      allocate(ene_h(uu%nums,qp%nspin))
   else
      allocate(qp%ene_x(qp%max_i,qp%nspin))
      allocate(ene_x(qp%max_i,qp%nspin))
      allocate(ene_h(qp%max_i,qp%nspin))
   endif
   allocate(qp%ene_h(qp%max_i,qp%nspin))

   qp%ene_dft_ks(1:qp%max_i,1:qp%nspin) =  uu%ene(1:qp%max_i,1:qp%nspin)
   qp%ene_dft_xc(1:qp%max_i,1:qp%nspin) =  uu%ene_xc(1:qp%max_i,1:qp%nspin)
   qp%ene_dft_h(1:qp%max_i,1:qp%nspin) =  uu%ene_lda_h(1:qp%max_i,1:qp%nspin)

   if(options%l_lda_hartree) then
       qp%ene_h(1:qp%max_i,1:qp%nspin) =  cmplx(uu%ene_lda_h(1:qp%max_i,1:qp%nspin),0.d0)
    else
!here calculate hartree parte
    endif

!calculate exchange part
   if(options%l_hf_energies) then
      call go_exchange(options,ene_x(:,1),ene_h(:,1),uu%nums)
      do ii=1,uu%nums
         qp%ene_x(ii,1)=ene_x(ii,1)
         qp%ene_hf(ii,1)=uu%ene(ii,1)-uu%ene_xc(ii,1)+ene_x(ii,1)
          write(stdout,*) 'ENE HF',ii, qp%ene_hf(ii,1)*RYTOEV
      enddo
   else
      write(stdout,*) 'ENE H', qp%ene_h(1:qp%max_i,1:qp%nspin)
      if(.not.options%l_lda_exchange) then
         call go_exchange(options,ene_x(:,1),ene_h(:,1),qp%max_i)
      else
!read from file
         call read_data_pw_exchange(ene_x,qp%max_i,options%prefix,qp%nspin)
         ene_h(1:qp%max_i,1:qp%nspin)=uu%ene_lda_h(1:qp%max_i,1:qp%nspin)
      endif
      do is=1,qp%nspin
         do ii=1,qp%max_i
            qp%ene_x(ii,is)=ene_x(ii,is)
            if(options%l_lda_hartree) then
               qp%ene_hf(ii,is)=uu%ene(ii,is)-uu%ene_xc(ii,is)+ene_x(ii,is)
            else
               qp%ene_h(ii,is)=cmplx(ene_h(ii,is),0.d0)
               qp%ene_hf(ii,is)=uu%ene(ii,is)-uu%ene_xc(ii,is)+ene_x(ii,is)-uu%ene_lda_h(ii,is)+ene_h(ii,is)
            endif
         enddo
      enddo
   endif

   call free_memory(uu)
   deallocate(ene_x,ene_h)
   return

 END SUBROUTINE create_hf
