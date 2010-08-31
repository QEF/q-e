!
! P.Umari Program GWW
!
  SUBROUTINE go_dressed_w(options)
!this subroutine reads the polarization on imaginary frequency
!and calculate the dressed interaction on imaginary frequency
!in case of an non-orthogonal basis set transform forth and back
!to the corresponding orthonormalized basis set
!if required uses the symmetrized dielectric matrix


   USE kinds,              ONLY : DP
   USE io_global,          ONLY : stdout
   USE input_gw,           ONLY : input_options
   USE basic_structures,   ONLY : v_pot,ortho_polaw,free_memory, head_epsilon
   USE polarization
   USE para_gww,           ONLY : is_my_pola
   USE mp,                 ONLY : mp_barrier, mp_sum
   USE w_divergence

   implicit none

   TYPE(input_options) :: options! for imaginary time range,number of samples

   TYPE(v_pot) :: vp!bare coulomb potential
   TYPE(polaw) :: pp,ww!polarization and dressed interaction
   TYPE(ortho_polaw) :: op!orthonormalization matrices
   TYPE(head_epsilon) :: he!the head (G=0,G=0) of the symmetrized dielectric matrix
   REAL(kind=DP), ALLOCATABLE :: agz(:)!elements  A_ij<G=0|\tilde{w^P_j> of the head
   REAL(kind=DP), ALLOCATABLE :: awing(:) !elements A_ij wing_j
   REAL(kind=DP), ALLOCATABLE :: awing_c(:) !elements A_ij wing_c_j
   INTEGER :: iw, label, ii, jj
   REAL(kind=DP), ALLOCATABLE :: inv_epsi(:)!for the heads of inverse dielectric matrices
   LOGICAL :: l_divergence
   TYPE(gv_time) :: gt!for handling the W(G=0,G=0) divergence


   call initialize_polaw(pp)
   call initialize_polaw(ww)

   allocate(inv_epsi(options%n+1))

!read coulomb potential
   write(stdout,*) 'ATTEZNIONE1'
   call flush_unit(stdout)
   if(options%w_divergence==2) then
      call read_data_pw_v(vp,options%prefix,options%debug,0,.true.)
   else
      call read_data_pw_v(vp,options%prefix,options%debug,0,.false.)
   endif

    write(stdout,*) 'ATTEZNIONE2'
    call flush_unit(stdout)

!read in orthonormalization matrix

   if(options%lnonorthogonal) then
      call read_data_pw_ortho_polaw(op,options%prefix)
      call orthonormalize_vpot_para(op,vp)
   endif
   write(stdout,*) 'ATTEZNIONE2.5'
   call flush_unit(stdout)
!if symmetric do symmetrize
      if(options%l_symm_epsilon) call square_root_polaw(vp%vmat,vp%numpw)



    write(stdout,*) 'ATTEZNIONE3'
    call flush_unit(stdout)

   allocate(agz(vp%numpw))
   allocate(awing(vp%numpw))
   allocate(awing_c(vp%numpw))
!if required read the head
   if(options%l_symm_epsilon .and. options%l_head_epsilon) then
      call read_data_pw_head_epsilon(he, options%prefix, options%l_wing_epsilon)
   endif

    write(stdout,*) 'ATTEZNIONE4'
    call flush_unit(stdout)

   if(options%w_divergence == 2) then
      l_divergence=.true.
   else
      l_divergence=.false.
   endif

   inv_epsi(:)=0.d0

!loop on imaginary frequencies samples

   do iw=0,options%n
      if(is_my_pola(iw)) then
         write(stdout,*) iw!ATTENZIONE
         call flush_unit(stdout)
         call read_polaw(iw,pp,options%debug)
         if(options%lnonorthogonal) then
            call orthonormalize(op,pp)
         endif
          write(stdout,*) 'call calculate_w'!ATTENZIONE
          call flush_unit(stdout)

          if(options%l_symm_epsilon .and. options%l_head_epsilon) then

!             agz(:)=0.d0
!             do ii=1,op%numpw
!                do jj=1,op%numpw
!                   agz(ii)=agz(ii)+op%on_mat(ii,jj)*he%gzero(jj)
!                enddo
!             enddo

             if(options%lnonorthogonal) then
                call dgemv('N',op%numpw,op%numpw,1.d0,op%on_mat,op%numpw,he%gzero,1,0.d0,agz,1)
             else
                agz(:)= he%gzero(:)
             endif
          endif

          if(options%l_symm_epsilon .and. options%l_wing_epsilon) then

!             awing(:)=0.d0
!             do ii=1,op%numpw
!                do jj=1,op%numpw
!                   awing(ii)=awing(ii)+op%on_mat(ii,jj)*he%wing(jj,iw+1)
!                enddo
!             enddo


             if(options%lnonorthogonal) then
                call dgemv('N',op%numpw,op%numpw,1.d0,op%on_mat,op%numpw,he%wing(:,iw+1),1,0.d0,awing,1)
                call dgemv('N',op%numpw,op%numpw,1.d0,op%on_mat,op%numpw,he%wing_c(:,iw+1),1,0.d0,awing_c,1)
             else
                awing(:)=he%wing(:,iw+1)
                awing_c(:)=he%wing_c(:,iw+1)

             endif
          else
             awing(:)=0.d0
             awing_c(:)=0.d0
          endif
          write(stdout,*) 'call calculate_w2'!ATTENZIONE
          call flush_unit(stdout)
         if(.not.options%l_head_epsilon) then
            call calculate_w(vp,pp,ww,options%xc_together,options%l_symm_epsilon,options%l_head_epsilon, &
                            agz, he%head(iw+1),l_divergence,inv_epsi(iw+1), options%l_wing_epsilon,awing)
          else
             call calculate_w_g(vp,pp,ww,options%xc_together,options%l_symm_epsilon,options%l_head_epsilon, &
               agz, he%head(iw+1),l_divergence,inv_epsi(iw+1), options%l_wing_epsilon,awing,awing_c)
          endif

         write(stdout,*) 'calculated w'!ATTENZIONE
         call flush_unit(stdout)
         if(options%lnonorthogonal) then
            call orthonormalize_inverse(op,ww)
         endif
         call write_polaw(ww,options%debug)

      endif
   enddo

   call mp_sum(inv_epsi)

   call free_memory(vp)
   call free_memory_polaw(pp)
   call free_memory_polaw(ww)
   if(options%lnonorthogonal) then
      call free_memory(op)
   endif
   if(options%l_symm_epsilon .and. options%l_head_epsilon) then
      call free_memory(he)
   endif
   deallocate(agz)
   deallocate(awing,awing_c)

!if required set up gv_time structure and write it on file
   if(options%w_divergence == 2) then
      call initialize_gv_time(gt)
      call read_data_pw_gv_time(gt, options%prefix)

      gt%ontime = .false.
      gt%omega=options%omega
      ii=gt%n+1
      do iw=1,gt%n+1
         gt%inv_epsi(iw)=inv_epsi(ii)
         ii=ii-1
      enddo
      ii=2
      do iw=gt%n+2,2*gt%n+1
         gt%inv_epsi(iw)=inv_epsi(ii)
         ii=ii+1
      enddo

      call write_gv_time(gt)
      call free_memory_gv_time(gt)
   endif
   deallocate(inv_epsi)

   return

 END SUBROUTINE go_dressed_w

  SUBROUTINE control_polarization(options)

   USE kinds,              ONLY : DP
   USE io_global,          ONLY : stdout
   USE input_gw,           ONLY : input_options
   USE basic_structures,   ONLY : v_pot,ortho_polaw,free_memory
   USE polarization
   USE para_gww,           ONLY : is_my_time
   USE mp,                 ONLY : mp_barrier

   implicit none

   TYPE(input_options) :: options! for imaginary time range,number of samples

   TYPE(v_pot) :: vp!bare coulomb potential
   TYPE(polaw) :: pp,ww!polarization and dressed interaction
   TYPE(ortho_polaw) :: op,opi!orthonormalization matrices
   INTEGER :: iw, label

   INTEGER :: i,j
   COMPLEX(kind=DP) :: sum


!read coulomb potential



!read in orthonormalization matrix

   if(options%lnonorthogonal) then
     call read_data_pw_ortho_polaw(op,options%prefix)
   endif

!loop on imaginary frequencies samples
   iw=15
   if(is_my_time(iw)) then
      write(stdout,*) iw!ATTENZIONE
      call read_polaw(iw,pp,options%debug)
      if(options%lnonorthogonal) then
         call orthonormalize(op,pp)
      endif
      sum = (0.d0, 0.d0)
      do i= 1, pp%numpw
         sum=sum+pp%pw(i,i)
      enddo
      write(stdout,*) 'SUM VCVC:', sum
   endif

   call mp_barrier

   call free_memory_polaw(pp)
   call free_memory_polaw(ww)
   if(options%lnonorthogonal) then
      call free_memory(op)
      call free_memory(opi)
   endif
   return

 END SUBROUTINE




