!
! Copyright (C) 2001-2013 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!

 SUBROUTINE create_quasi_particles(options,qp,se)
!given the expansion coeffcients, calculates in a perturbative
!way without self-consistency correction the quasi-particles energies
!relavant arrays are already allocates and set by subroutine create


   USE io_global,        ONLY : stdout
   USE basic_structures, ONLY : wannier_u, free_memory
   USE expansion,        ONLY : self_expansion, value_on_frequency, derivative_on_frequency,value_on_frequency_complex
   USE input_gw,         ONLY : input_options
   USE constants,        ONLY : tpi, RYTOEV
   USE energies_gww,         ONLY : quasi_particles
   USE kinds,            ONLY : DP

   implicit none

   TYPE(input_options) :: options! for prefix
   TYPE(quasi_particles) :: qp!the descriptor to be build
   TYPE(self_expansion) :: se!the descriptor for the multipole expansion


   INTEGER :: ii,jj, it,is
   TYPE(wannier_u) :: uu
   COMPLEX(kind=DP) :: zz, sigmac,dsigmac
   REAL(kind=DP)  :: offset
   REAL(kind=DP), ALLOCATABLE  :: delta_ene(:)



!read in DFT energies
   call read_data_pw_u(uu,options%prefix)

!loop on spin

   do is=1,uu%nspin

      if(.not. options%l_hf_energies) then
         if(uu%nums_occ(is) == 0) then
            offset=-2.d0
         else
            if(uu%nums > uu%nums_occ(is)) then
               if(options%l_lda_hartree) then
                  offset=-(uu%ene(uu%nums_occ(is)+1,is)+uu%ene(uu%nums_occ(is),is))/2.d0
               else
                  offset=-(uu%ene(uu%nums_occ(is)+1,is)+dble(qp%ene_h(uu%nums_occ(is)+1,is))-qp%ene_dft_h(uu%nums_occ(is)+1,is)&
                       & +uu%ene(uu%nums_occ(is),is) +dble(qp%ene_h(uu%nums_occ(is),is))-qp%ene_dft_h(uu%nums_occ(is),is))/2.d0
               endif
            else
               if(options%l_lda_hartree) then
                  offset=-uu%ene(uu%nums_occ(is),is)
               else
                  offset=-(uu%ene(uu%nums_occ(is),is)+dble(qp%ene_h(uu%nums_occ(is),is))-qp%ene_dft_h(uu%nums_occ(is),is))
               endif
            endif
         endif
      else
         if(uu%nums > uu%nums_occ(is)) then
            offset=-(qp%ene_hf(uu%nums_occ(is)+1,is)+qp%ene_hf(uu%nums_occ(is),is))/2.d0
         else
            offset=-qp%ene_hf(uu%nums_occ(is),is)
         endif
      endif
 !  call free_memory(uu)

!set scissor if required

      allocate(delta_ene(options%max_i))
      delta_ene(:)=0.d0
      if(options%l_scissor) then
         do ii=1,uu%nums_occ(is)
            delta_ene(ii)=-options%scissor(1)/RYTOEV
         enddo
         do ii=uu%nums_occ(is)+1,uu%nums
            delta_ene(ii)=-options%scissor(2)/RYTOEV
         enddo
      endif



      do ii=1,qp%max_i
         if(.not. options%l_hf_energies) then
            if(options%l_lda_hartree) then
               call value_on_frequency(se,ii,qp%ene_dft_ks(ii,is)+offset,sigmac,is)
               call derivative_on_frequency(se,ii,qp%ene_dft_ks(ii,is)+offset,dsigmac,is)
            else
               call value_on_frequency(se,ii,qp%ene_dft_ks(ii,is)+offset+dble(qp%ene_h(ii,is))-qp%ene_dft_h(ii,is),sigmac,is)
               call derivative_on_frequency(se,ii,qp%ene_dft_ks(ii,is)+offset+dble(qp%ene_h(ii,is))-qp%ene_dft_h(ii,is),dsigmac,is)
            endif
         else
            call value_on_frequency(se,ii,qp%ene_hf(ii,is)+offset,sigmac,is)
            call derivative_on_frequency(se,ii,qp%ene_hf(ii,is)+offset,dsigmac,is)
         endif
         write(stdout,*) 'value, zeta:',ii,sigmac,dsigmac,is
         zz=(1.d0,0.d0)-dsigmac
         if(.not. options%l_hf_energies) then
            qp%ene_gw(ii,is)=qp%ene_dft_ks(ii,is)+offset +qp%ene_h(ii,is)-qp%ene_dft_h(ii,is)+&
                 & (sigmac+delta_ene(ii)+qp%ene_x(ii,is)-qp%ene_dft_xc(ii,is) )/zz
         else
            qp%ene_gw(ii,is)=qp%ene_hf(ii,is)+offset+(sigmac+delta_ene(ii))/zz
         endif
         write(stdout,*) 'XC-DFT energy',ii,qp%ene_dft_xc(ii,is)
         write(stdout,*) 'H-DFT energy',ii,qp%ene_dft_h(ii,is)*RYTOEV,qp%ene_h(ii,is)*RYTOEV
         write(stdout,*) 'GW-PERT energy', ii,real(qp%ene_gw(ii,is)-offset)*RYTOEV
         qp%ene_gw_pert(ii,is)=qp%ene_gw(ii,is)-offset
!self-consistency loop
         do it=1,10
            call value_on_frequency_complex(se,ii,qp%ene_gw(ii,is),sigmac,is)
            sigmac=sigmac+delta_ene(ii)
            write(stdout,*) 'Iteration energy',it,sigmac
            if(.not. options%l_hf_energies) then
               qp%ene_gw(ii,is)=qp%ene_dft_ks(ii,is)+offset+sigmac+qp%ene_x(ii,is)-qp%ene_dft_xc(ii,is) &
                         &   +qp%ene_h(ii,is)-qp%ene_dft_h(ii,is)
            else
               qp%ene_gw(ii,is)=qp%ene_hf(ii,is)+offset+sigmac
            endif
         enddo
         qp%ene_gw(ii,is)= qp%ene_gw(ii,is)-offset
         FLUSH(stdout)
      enddo


      deallocate(delta_ene)

   enddo!spin
   call free_memory(uu)
   return
 END SUBROUTINE create_quasi_particles


 SUBROUTINE create_quasi_particle_on_real(options,qp,sr)
!given the self-energy on real axis calculate the GW levels
   USE io_global,        ONLY : stdout
   USE basic_structures, ONLY : wannier_u, free_memory
   USE self_energy_storage
   USE input_gw,         ONLY : input_options
   USE constants,        ONLY : tpi, RYTOEV
   USE energies_gww,         ONLY : quasi_particles
   USE kinds,            ONLY : DP

   implicit none

   TYPE(input_options) :: options! for prefix  
   TYPE(quasi_particles) :: qp!the descriptor to be build
   TYPE(self_on_real) :: sr!the descriptor for the self_energy                                                                                   


   INTEGER :: ii,jj, it,is,ierr
   TYPE(wannier_u) :: uu
   COMPLEX(kind=DP) :: zz, sigmac,dsigmac,energy
  
!read in DFT energies
   call read_data_pw_u(uu,options%prefix)

!loop on spin

   do is=1,uu%nspin
      do ii=sr%i_min,sr%i_max
         energy=dcmplx(qp%ene_dft_ks(ii,is),0.d0)
         energy=qp%ene_gw(ii,is)!ATTENZIONE
         call  self_on_real_value(sr,ii,is,energy,sigmac,ierr)
         if(ierr/=0) then
            write(stdout,*) 'OUT OF RANGE:self_on_real_value',energy
            FLUSH(stdout)
            !stop!ATTENZIONE
         endif
         write(stdout,*) 'Iteration energy 0', dble(qp%ene_gw(ii,is))
         qp%ene_gw(ii,is)=qp%ene_dft_ks(ii,is)+sigmac+qp%ene_x(ii,is)-qp%ene_dft_xc(ii,is) 
         write(stdout,*) 'Iteration energy 1', dble(qp%ene_gw(ii,is))
         
!self-consistency loop 
         do it=1,1000
            call  self_on_real_value(sr,ii,is, qp%ene_gw(ii,is),sigmac,ierr)
            if(ierr/=0) then
               write(stdout,*) 'OUT OF RANGE:self_on_real_value',it,qp%ene_gw(ii,is)
               FLUSH(stdout)
               !stop!ATTENZIONE
            endif
            write(stdout,*) 'Iteration energy',it,sigmac,dble(qp%ene_gw(ii,is))
            FLUSH(stdout)
            qp%ene_gw(ii,is)=qp%ene_dft_ks(ii,is)+sigmac+qp%ene_x(ii,is)-qp%ene_dft_xc(ii,is) 
         enddo
      enddo
   
   enddo

   return

 END SUBROUTINE create_quasi_particle_on_real
