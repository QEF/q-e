! Author: P. Umari
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


   INTEGER :: ii,jj, it
   TYPE(wannier_u) :: uu
   COMPLEX(kind=DP) :: zz, sigmac,dsigmac
   REAL(kind=DP)  :: offset
   REAL(kind=DP), ALLOCATABLE  :: remainder(:)



!read in DFT energies
   call read_data_pw_u(uu,options%prefix)


   if(.not. options%l_hf_energies) then
      if(uu%nums > uu%nums_occ) then
         if(options%l_lda_hartree) then
            offset=-(uu%ene(uu%nums_occ+1)+uu%ene(uu%nums_occ))/2.d0
         else
            offset=-(uu%ene(uu%nums_occ+1)+dble(qp%ene_h(uu%nums_occ+1))-qp%ene_dft_h(uu%nums_occ+1)&
               & +uu%ene(uu%nums_occ) +dble(qp%ene_h(uu%nums_occ))-qp%ene_dft_h(uu%nums_occ))/2.d0
         endif
      else
         if(options%l_lda_hartree) then
            offset=-uu%ene(uu%nums_occ)
         else
            offset=-(uu%ene(uu%nums_occ)+dble(qp%ene_h(uu%nums_occ))-qp%ene_dft_h(uu%nums_occ))
         endif
      endif
   else
      if(uu%nums > uu%nums_occ) then
         offset=-(qp%ene_hf(uu%nums_occ+1)+qp%ene_hf(uu%nums_occ))/2.d0
      else
         offset=-qp%ene_hf(uu%nums_occ)
      endif
   endif
   call free_memory(uu)

!set remainders

   allocate(remainder(options%max_i))
   if(options%remainder==3 .or.options%remainder==4) then
      remainder(:)=qp%ene_remainder(:)
   else
      remainder(:)=0.d0
   endif



   do ii=1,qp%max_i
      if(.not. options%l_hf_energies) then
         if(options%l_lda_hartree) then
            call value_on_frequency(se,ii,qp%ene_dft_ks(ii)+offset,sigmac)
            call derivative_on_frequency(se,ii,qp%ene_dft_ks(ii)+offset,dsigmac)
         else
            call value_on_frequency(se,ii,qp%ene_dft_ks(ii)+offset+dble(qp%ene_h(ii))-qp%ene_dft_h(ii),sigmac)
            call derivative_on_frequency(se,ii,qp%ene_dft_ks(ii)+offset+dble(qp%ene_h(ii))-qp%ene_dft_h(ii),dsigmac)
         endif
      else
         call value_on_frequency(se,ii,qp%ene_hf(ii)+offset,sigmac)
         call derivative_on_frequency(se,ii,qp%ene_hf(ii)+offset,dsigmac)
      endif
      write(stdout,*) 'value, zeta:',ii,sigmac,dsigmac
      zz=(1.d0,0.d0)-dsigmac
      if(.not. options%l_hf_energies) then
         qp%ene_gw(ii)=qp%ene_dft_ks(ii)+offset +qp%ene_h(ii)-qp%ene_dft_h(ii)+(sigmac+remainder(ii)+qp%ene_x(ii)-qp%ene_dft_xc(ii)&
                     &  )/zz
      else
         qp%ene_gw(ii)=qp%ene_hf(ii)+offset+(sigmac+remainder(ii))/zz
      endif
      write(stdout,*) 'XC-DFT energy',ii,qp%ene_dft_xc(ii)
      write(stdout,*) 'H-DFT energy',ii,qp%ene_dft_h(ii)*RYTOEV,qp%ene_h(ii)*RYTOEV
      write(stdout,*) 'GW-PERT energy', ii,real(qp%ene_gw(ii)-offset)*RYTOEV
      qp%ene_gw_pert(ii)=qp%ene_gw(ii)-offset
!self-consistency loop
      do it=1,10
         call value_on_frequency_complex(se,ii,qp%ene_gw(ii),sigmac)
         sigmac=sigmac+remainder(ii)
         write(stdout,*) 'Iteration energy',it,sigmac
         if(.not. options%l_hf_energies) then
            qp%ene_gw(ii)=qp%ene_dft_ks(ii)+offset+sigmac+qp%ene_x(ii)-qp%ene_dft_xc(ii) &
                         &   +qp%ene_h(ii)-qp%ene_dft_h(ii)
         else
            qp%ene_gw(ii)=qp%ene_hf(ii)+offset+sigmac
         endif
      enddo
      qp%ene_gw(ii)= qp%ene_gw(ii)-offset
   enddo


   deallocate(remainder)


   return
 END SUBROUTINE create_quasi_particles
