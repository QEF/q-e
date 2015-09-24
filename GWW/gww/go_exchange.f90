!
! Copyright (C) 2001-2013 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!

 SUBROUTINE go_exchange_main( options, qp)
!this subroutines does:
!1)creates and writes green function a t=0+
!2)if required creates and writes contractions
!3)setup qp and its exchange and hf arrays

   USE energies_gww,           ONLY : quasi_particles
   USE compact_product,    ONLY : contraction, free_memory_contraction,do_contraction, write_contraction,&
                                                   &do_contraction_index_state
   USE basic_structures,   ONLY : q_mat, wannier_u,free_memory
   USE green_function,     ONLY : green, free_memory_green, write_green, create_green_part,initialize_green
   USE para_gww,           ONLY : is_my_time
   USE mp,                 ONLY : mp_barrier
   USE mp_world,           ONLY : world_comm
   USE input_gw,           ONLY : input_options
   USE io_global,          ONLY : stdout, ionode
   USE kinds,              ONLY : DP
   USE constants,          ONLY : RYTOEV


   implicit none

   TYPE(input_options), INTENT(in)    ::  options! program options
   TYPE(quasi_particles), INTENT(out) ::  qp!quasi particle structure to be initialized with HF stuff


   TYPE(contraction) :: cr!to speed up calculation
   TYPE(green)       :: gg!green function
   TYPE(q_mat)       :: qm!overlap of wannier products
   TYPE(wannier_u)   :: uu!transformation matrix ks to wannier
   
   REAL(kind=DP)     :: dumm(1)

   nullify(cr%numl)
   nullify(cr%l)
   nullify(cr%q)

   call initialize_green(gg)


   if(options%l_verbose) write(stdout,*) 'Routine go_exchange main1'
   FLUSH(stdout)

   !read U matrix
   call read_data_pw_u(uu,options%prefix)
!read overlap matrix Q
  
 !  call read_data_pw_q(qm,options%prefix,options%l_self_from_pola)

     if(options%l_verbose) write(stdout,*) 'Routine go_exchange main2'
     FLUSH(stdout)

  if(is_my_time(0)) then
     call create_green_part(gg,uu,0.d0, options%debug,.false.,.false.,dumm)
     gg%label=0
     call write_green(gg, options%debug)
  endif
  call free_memory_green(gg)
  if(options%use_contractions) then
     call read_data_pw_q(qm,options%prefix,options%l_self_from_pola)
!in the following max_i defines the max number of KS states not appropriate for HF
!the contraction qr can be defined HERE in another way wp:w_V w_i
     if(.not.options%l_contraction_single_state) then
        call  do_contraction(qm,uu,cr, options%max_i)
        if(options%l_verbose) write(stdout,*) 'Routine go_exchange main2.2'
        call free_memory(uu)
        call free_memory(qm)
        if(options%l_verbose) write(stdout,*) 'Routine go_exchange main2.3'
        call  write_contraction(cr,options)
     else
        call  do_contraction_index_state(qm,uu,options%max_i, options)
        call free_memory(uu)
        call free_memory(qm)
    endif
  endif
  call mp_barrier( world_comm )

  if(.not.options%use_contractions) then
     call free_memory(uu)
  !    call free_memory(qm)
  else
     if(.not.options%l_contraction_single_state) call free_memory_contraction(cr)
  endif



  if(options%l_verbose) write(stdout,*) 'Routine go_exchange main3'

  call create_hf(options, qp)

  if(options%l_verbose) write(*,*) 'go_exchange main hf_ene', qp%ene_hf(2,1)*RYTOEV!ATTENZIONE

  return
   
END SUBROUTINE go_exchange_main









  SUBROUTINE go_exchange(options, ene_x, ene_h, n_max)
!this subroutine calculates the terms, in imaginary time=0
!<\Psi_i|\Sigma(it)_HF|\Psi_j> 
!=O^{P}_n,kl G_{lm}V_{n,o} O^{P}_o,mp U_ki U^{+}_j,p
!for n_max states
!if required calculates also the hartree terms
!<\Psi_i|V_H|\Psi_i>
!and displays result on screen


   USE kinds,                ONLY : DP
   USE io_global,            ONLY : stdout 
   USE basic_structures,     ONLY : wannier_u, q_mat, v_pot, ortho_polaw,free_memory, v_pot_prim, wannier_u_prim
   USE green_function,       ONLY : green, read_green, free_memory_green, initialize_green
   USE input_gw,             ONLY : input_options
   USE compact_product
   USE polarization
   USE para_gww,             ONLY : is_my_state
   USE mp,                   ONLY : mp_sum
   USE mp_world,             ONLY : world_comm
   USE constants, ONLY : RYTOEV

   implicit none

   TYPE(input_options) :: options
   REAL(kind=DP) :: ene_x(n_max)!where to store calculated diagonal values
   REAL(kind=DP) :: ene_h(n_max)!where to store calculated diagonal values
   INTEGER :: n_max!max number of states to be considered

   INTEGER  :: i,j !which element of self enrgy to be calculated
   REAL(kind=DP) :: sene!self energy element
   TYPE(q_mat)  :: qm!descriptors of overlaps of othonormalized wannier producs with wannier products
   TYPE(wannier_u) :: uu!descriptor of transformation matrix from KS states to wanniers
   TYPE(green) :: gf!descriptor of green function 
   TYPE(v_pot) :: vp!bare interaction
   TYPE(contraction) :: cr!for contracted products scheme
   TYPE(contraction_index) :: cri!for contracted products scheme, index part
   TYPE(contraction_state) :: crs!for contracted products scheme, state part
   TYPE(ortho_polaw) :: op!orthonormalization matrix
   REAL(kind=DP), ALLOCATABLE :: qg(:,:)!for the product Q^{P}_{n,l,i}G{l,m} SUPPOSED TO BE REAL
   REAL(kind=DP), ALLOCATABLE :: qu(:)!for the product Q^{P}_{m,v,v'}U{v,v'}SUPPOSED TO BE REAL
   REAL(kind=DP), ALLOCATABLE :: ju(:)!for the product I_{j,m}Q^{P}_{m,v,v'}U{v,v'}SUPPOSED TO BE REAL
   TYPE(v_pot_prim) :: vp_prim
   TYPE(wannier_u_prim) :: wu
   TYPE(v_pot) :: ident

   INTEGER :: k,l,m,n,o,p,v,vv
   INTEGER :: nw,ow
   REAL(kind=DP) :: o_n,o_o

   call initialize_green(gf)


   call read_green(0, gf, options%debug,.false.)

   write(stdout,*) 'GF PART',gf%l_part
   call read_data_pw_v(vp,options%prefix,options%debug,0,.false.)
   if(options%use_contractions) then
      if(.not.options%l_contraction_single_state) then
         call read_contraction(cr,options)
      else
         call read_contraction_index(cri,options)
      endif
   endif
   write(stdout,*) 'contraction read'
   
   if(.not.options%use_contractions) then
      call read_data_pw_u(uu,options%prefix)
      call read_data_pw_q(qm,options%prefix,.false.)
      if(gf%nums /= uu%nums) then
         write(stdout,*) 'Routine self_energy: same nums required'
         stop
      endif
      if(qm%numpw /= vp%numpw) then
         write(stdout,*) 'Routine self_energy: same numpw required',qm%numpw,vp%numpw
         stop
      endif
   endif

   write(stdout,*) 'invert potential'
  if(options%lnonorthogonal) then
     call read_data_pw_ortho_polaw(op,options%prefix)
     call orthonormalize_vpot_para(op,vp)
     call orthonormalize_vpot_inverse_para(op,vp)
     if(options%l_lda_hartree) call free_memory(op)
   endif
  
    write(stdout,*) 'invert potential inverted'
   
   ene_x(:)=0.d0
   ene_h(:)=0.d0

   if(.not.options%l_lda_hartree) then
      call read_data_pw_u(uu,options%prefix)
      call read_data_pw_v_pot_prim(vp_prim,options%prefix,.false.)
      call read_data_pw_u_prim(wu,options%prefix)
      allocate(qu(cri%numpw))
      qu(:)=0.d0
      do v=1,uu%nums_occ(1)
         crs%state=v
         call read_contraction_state(cri,crs,options)
         do m=1,cri%numpw
            do vv=1,cri%numl(m)
               qu(m)=qu(m)+crs%q(m,vv)*dble(uu%umat(v,cri%l(vv,m),1))*2.d0!2. is for spin multiplicity
            enddo
         enddo
         call free_memory_contraction_state(crs)
      enddo

      ident%numpw=cri%numpw
      allocate(ident%vmat(ident%numpw,ident%numpw))
      ident%vmat(:,:)=0.d0
      do m=1,ident%numpw
         ident%vmat(m,m)=1.d0
      enddo
      call orthonormalize_vpot_inverse_para(op,ident)
      allocate(ju(vp_prim%numpw))
      ju(:)=0.d0
      do m=1,vp_prim%numpw
         do n=1,vp_prim%numpw
            ju(m)=ju(m)+ident%vmat(m,n)*qu(n)
         enddo
         write(*,*) 'JU',ju(m)!ATTENZIONE
      enddo
      call free_memory(ident)
      call free_memory(op)
   endif
   if(options%whole_s) then
     write(stdout,*)'Routine go_exchange: whole_s not implemented yet'
     stop
   else
     do i=1,n_max
        if(is_my_state(i)) then
           j=i
!check consistency

           sene=0.d0
           if(.not.options%use_contractions) then
              write(stdout,*) 'ONLY CONTRACTIONS IMPLEMENTED'
              stop
           else
              if(.not.options%l_contraction_single_state) then
                 write(*,*) 'Interno', cr%numpw,cr%nums
                 allocate(qg(cr%numpw,cr%nums))
                 qg(:,:)=0.d0
                 
                 do n=1,cr%numpw!loop on orthonormalized wannier products
                    do m=1,cr%nums
                       do l=1,cr%numl(n)
                          qg(n,m)=qg(n,m)+dble(cr%q(n,l,i))*gf%gf_p(cr%l(l,n),m,1)
                       enddo
                    enddo
                 enddo
                 
                 sene=0.d0
                 do n=1,cr%numpw!loop on orthonormalized wannier products
                    do o=1,cr%numpw!loop on orthonormalized wannier products
                       do m=1,cr%numl(o)
                          sene=sene+qg(n,cr%l(m,o))*vp%vmat(n,o)*dble(cr%q(o,m,j))
                       enddo
                    enddo
                 enddo
                 deallocate(qg)
              else
                 crs%state=i
                 write(stdout,*) 'read state', i
                 call read_contraction_state(cri,crs,options)
                 write(stdout,*) 'Interno state', cri%numpw,cri%nums
                 allocate(qg(cri%numpw,cri%nums))
                 qg(:,:)=0.d0

                 do n=1,cri%numpw!loop on orthonormalized wannier products
                    do m=1,cri%nums
                       do l=1,cri%numl(n)
                          qg(n,m)=qg(n,m)+crs%q(n,l)*gf%gf_p(cri%l(l,n),m,1)
                       enddo
                    enddo
                 enddo

                 sene=0.d0
                 do n=1,cri%numpw!loop on orthonormalized wannier products
                    do o=1,cri%numpw!loop on orthonormalized wannier products
                       do m=1,cri%numl(o)
                          sene=sene+qg(n,cri%l(m,o))*vp%vmat(n,o)*crs%q(o,m)
                       enddo
                    enddo
                 enddo
                 deallocate(qg)

                 if(.not.options%l_lda_hartree) then
!calculate hartree term
                    ene_h(i)=0.d0
                    if(i<=uu%nums_occ(1)) then
                       do l=1,cri%numpw
                          do k=1,cri%numl(l)
                             do m=1,cri%numpw
                                ene_h(i)=ene_h(i)+uu%umat(i,cri%l(k,l),1)*crs%q(l,k)*qu(m)*vp%vmat(l,m)
                             enddo
                          enddo
                       enddo
                    else
                       do l=1,vp_prim%numpw_prim
                          do m=1,cri%numpw
                             ene_h(i)=ene_h(i)+wu%umat(i-wu%nums_occ,vp_prim%ij(1,l))*&
          & uu%umat(i,vp_prim%ij(2,l),1)*vp_prim%vmat(l,m)*ju(m)
                          enddo
                       enddo
                    endif
                    WRITE(STDOUT,*) 'Hartree Energy (ryd) :',i,ene_h(i)*RYTOEV
                 endif
                 call free_memory_contraction_state(crs)
              endif
           endif
           WRITE(STDOUT,*) 'SENE :' ,SENE,gf%factor!ATTENZIONE
           sene=sene*dble(gf%factor*(0.d0,1.d0))
           ene_x(i)=sene
           write(stdout,*) 'Exchange  energies', i,sene
        endif
     enddo
     call mp_sum(ene_x(1:options%max_i),world_comm)
     if(.not.options%l_lda_hartree)  call mp_sum(ene_h(1:options%max_i),world_comm)
  endif

    if(.not.options%use_contractions) then
       call free_memory(qm)
       call free_memory(uu)
    else
       if(.not.options%l_contraction_single_state) then
          call free_memory_contraction(cr)
       else
          call free_memory_contraction_index(cri)
       endif
       if(.not.options%l_lda_hartree) then
          call free_memory(uu)
          call free_memory(wu)
          call free_memory(vp_prim)
          deallocate(qu)
          deallocate(ju)
       endif
    endif
    call free_memory(vp)
    call free_memory_green(gf)
    return
  END SUBROUTINE go_exchange
  

























 
        

