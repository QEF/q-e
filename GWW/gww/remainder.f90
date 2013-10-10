!
! Copyright (C) 2001-2013 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!

!Program GWW P. Umari

SUBROUTINE remainder(options, qp)
!this subroutine calculates the remainder

   USE  constants,         ONLY : eps8
   USE io_global,          ONLY : stdout, ionode, ionode_id
   USE input_gw,           ONLY : input_options
   USE basic_structures,   ONLY : q_mat, wannier_u, wp_psi, wp_psi_cutoff_index,wp_psi_cutoff_data,free_memory
   USE green_function,     ONLY : green,read_green,free_memory_green, initialize_green
   USE polarization,       ONLY : polaw,free_memory_polaw,read_polaw, initialize_polaw
   USE compact_product
   USE mp,                 ONLY : mp_sum, mp_bcast
   USE mp_world,           ONLY : world_comm
   USE para_gww,           ONLY : is_my_time, is_my_pola, is_my_state
   USE energies_gww,           ONLY : quasi_particles
   USE constants,          ONLY : RYTOEV
   USE energies_gww,           ONLY : quasi_particles


   implicit none


   INTERFACE
      SUBROUTINE set_data_wp_psi_cutoff(pw_red,pw,wpi)
        USE kinds, ONLY : DP
        USE basic_structures,   ONLY : wp_psi_cutoff_index
        USE polarization, ONLY : polaw

        COMPLEX(kind=DP), DIMENSION(:), POINTER :: pw_red
        TYPE(polaw) :: pw!data to be contracted
        TYPE(wp_psi_cutoff_index) :: wpi !indices
      END SUBROUTINE set_data_wp_psi_cutoff

      SUBROUTINE self_energy_remainder_cutoff(state,rem,wp,pw_red)
        USE kinds,                ONLY : DP
        USE basic_structures,     ONLY : wp_psi_cutoff_data
        INTEGER  :: state
        COMPLEX(kind=DP) :: rem
        COMPLEX(kind=DP), DIMENSION(:), POINTER :: pw_red
        TYPE(wp_psi_cutoff_data) :: wp
      END SUBROUTINE self_energy_remainder_cutoff

   END INTERFACE



   TYPE(input_options), INTENT(in) :: options
   TYPE(quasi_particles), INTENT(inout) :: qp


   TYPE(green)     :: gg,gm!green function
   TYPE(q_mat)     :: qm!overlap of orthonormalized wannier products with wannier products
   TYPE(polaw)     :: ww!dressed interaction
   TYPE(wannier_u) :: uu!transformation matrix ks to wannier
   TYPE(contraction) :: cr!to speed up calculation
   TYPE(wp_psi) :: wp!for remainder calculations
   TYPE(wp_psi_cutoff_data) :: wpc!for remainder calculations with cutoff
   TYPE(wp_psi_cutoff_index) :: wpci!for remainder calculations with cutoff, index
   TYPE(contraction_index) :: cri! index of contraction
   TYPE(contraction_state) :: crs!state contraction data

   REAL(kind=DP) :: time
   INTEGER       :: iw,ii,jj
   REAL(kind=DP) :: offset
   COMPLEX(kind=DP) :: sca
   COMPLEX(kind=DP), DIMENSION(:), POINTER :: pw_red


   write(stdout,*) 'enter remainder COH'

!allocates
   allocate(qp%ene_remainder(options%max_i,qp%nspin))
   
   call initialize_green(gg)
   call initialize_green(gm)
   call initialize_polaw(ww)


      
!read U matrix
  call read_data_pw_u(uu,options%prefix)
!read overlap matrix Q
  call read_data_pw_q(qm,options%prefix,.false.)
  if(options%use_contractions) then
     if(.not.options%l_contraction_single_state) then
        write(stdout,*) 'call do_contraction'!ATTENZIONE
        call  do_contraction(qm,uu,cr, options%max_i)
        write(stdout,*) 'done do_contraction'!ATTENZIONE
        call  write_contraction(cr,options)
        write(stdout,*) 'done do_contraction'!ATTENZIONE
     else
!contraction index and states already available on disk
        call read_contraction_index(cri, options)
     endif
  endif



  write(stdout,*) 'enter remainder COH'
  

  if(options%remainder /= 4) then!not needed if calculated through pw
     if(options%l_remainder_cutoff) then
        call read_data_pw_wp_psi_cutoff_index(wpci,options%prefix)
        call read_data_pw_wp_psi_cutoff_data(wpci,wpc, options%prefix)
     else
        call read_data_pw_wp_psi(wp,options%prefix)
     endif
  endif

         
  call read_green(0,gg,options%debug,.false.)
  call read_green(0,gm,options%debug,.true.)
  

  call read_polaw(0,ww,options%debug,options%l_verbose)

  write(stdout,*) 'POLAW FACTOR', ww%factor!ATTENZIONE

  if(options%remainder /= 4) then
     if(options%l_remainder_cutoff) then
        call set_data_wp_psi_cutoff(pw_red,ww,wpci)
        WRITE(*,*) 'PW_RED OUT', pw_red(1)
     endif
  endif
      

  time=0.d0

  qp%ene_remainder(:,:) =(0.d0,0.d0)
  do ii=1,options%max_i
     if(is_my_state(ii)) then
        if(.not.options%use_contractions) then
           call self_energy(ii,ii,sca,time,qm,uu,gg,ww)
        else
           if(.not.options%l_contraction_single_state) then
              call self_energy_contraction(ii,ii,sca,time,cr,gg,ww)
           else
              crs%state=ii
               write(stdout,*) 'Call read_contraction_state'
               call read_contraction_state(cri,crs,options)
               call self_energy_contraction_state(ii,ii,sca,time,cri,crs,gg,ww)
            endif
        endif

!sene changes sign because we are on the negative axes!!
        qp%ene_remainder(ii,1)=qp%ene_remainder(ii,1)+0.5d0*dble(sca)

        write(*,*) 'REMAINDER SENE 1', ii, 0.5d0*sca

        if(.not.options%use_contractions) then
           call self_energy(ii,ii,sca,time,qm,uu,gm,ww)
        else
           if(.not.options%l_contraction_single_state) then
              call self_energy_contraction(ii,ii,sca,time,cr,gm,ww)
           else
              call self_energy_contraction_state(ii,ii,sca,time,cri,crs,gm,ww)
              write(stdout,*) 'Call free_memory_contraction_state'
              call free_memory_contraction_state(crs)
           endif
         endif
         qp%ene_remainder(ii,1)=qp%ene_remainder(ii,1)-0.5d0*dble(sca)

         write(*,*) 'REMAINDER SENE 2', ii, 0.5d0*sca


         if(options%remainder /= 4) then
            if(options%l_remainder_cutoff) then
               call self_energy_remainder_cutoff(ii,sca,wpc,pw_red)
            else
               call self_energy_remainder(ii,sca,time,wp,ww)
            endif
            
            qp%ene_remainder(ii,1)=qp%ene_remainder(ii,1)+0.5d0*dble(sca)
            write(*,*) 'REMAINDER SENE 3', ii, 0.5d0*sca
         endif
      endif
   enddo
      


 

   call free_memory_polaw(ww)
   call free_memory_green(gg)
   call free_memory_green(gm)
   write(*,*)  'in of cycle'!ATTENZIONE
   
   
   if(options%remainder /= 4) then
      if(options%l_remainder_cutoff) then
         deallocate(pw_red)
         call free_memory(wpc)
         call free_memory(wpci)
      else
         call free_memory(wp)
      endif
   endif
   call free_memory(uu)
   call free_memory(qm)
   if(.not.options%l_contraction_single_state) then
      call free_memory_contraction(cr)
   else
      call free_memory_contraction_index(cri)
   endif

  call mp_sum(qp%ene_remainder(:,1),world_comm)

  if(options%lconduction) call addconduction_remainder(qp, options)

  if(ionode) then
     do ii=1,options%max_i
        write(stdout,*) 'CORRECTION COH', ii, qp%ene_remainder(ii,1)*RYTOEV
     enddo
  endif

  return
end SUBROUTINE remainder


SUBROUTINE addconduction_remainder(qp, options)
!this subroutine adds to the self_energy of conduction states
!on negative imaginary times, the part due to terms \Psi_c'\Psic\w_P

    USE kinds,         ONLY : DP
    USE io_global,     ONLY : stdout, ionode, ionode_id
    USE input_gw,      ONLY : input_options
    USE basic_structures,  ONLY : v_pot,wannier_u_prim, v_pot_prim,free_memory, ortho_polaw
    USE green_function,    ONLY : green, read_green, free_memory_green,initialize_green
    USE polarization,      ONLY : polaw, free_memory_polaw, read_polaw, invert_v_pot, invert_ortho_polaw,&
         & orthonormalize_inverse, orthonormalize_vpot_para
    USE mp,                ONLY : mp_bcast
    USE mp_world,          ONLY : world_comm
    USE para_gww,          ONLY : is_my_pola
    USE energies_gww,          ONLY : quasi_particles

    implicit none

    TYPE(input_options) :: options
    TYPE(quasi_particles)  :: qp

    TYPE(v_pot) :: vp,vpi
    TYPE(ortho_polaw) :: op,opi
    TYPE(polaw) :: ww!dressed interaction
    TYPE(wannier_u_prim) :: wup
    TYPE(v_pot_prim) :: vpp
    TYPE(green) :: gg

    INTEGER iw,jw,kw,it,ii
    REAL(kind=DP), ALLOCATABLE :: wtemp(:,:)
    REAL(kind=DP), ALLOCATABLE :: cp(:,:,:) !arrys for contraction c',c, numpw
    REAL(kind=DP), ALLOCATABLE :: qg(:,:)
    COMPLEX(kind=DP), ALLOCATABLE :: sene(:)

    REAL(kind=DP), ALLOCATABLE :: gf_t(:,:)
    REAL(kind=DP), ALLOCATABLE :: pwcp_t(:,:)
    REAL(kind=DP), EXTERNAL :: ddot

    nullify(vp%vmat)
    nullify(vpi%vmat)
    nullify(op%on_mat)
    nullify(opi%on_mat)
    nullify(ww%pw)
    nullify(wup%umat)
    nullify(vpp%ij)
    nullify(vpp%vmat)
    

    call initialize_green(gg)





!read coulombian potential and calculate inverse



    write(stdout,*) 'Routine add_cunduction_remainder'

    call read_data_pw_u_prim(wup,options%prefix)
    call read_data_pw_v_pot_prim(vpp, options%prefix,.false.)


    allocate(sene(options%max_i-wup%nums_occ))
    sene(:)=(0.d0,0.d0)


!set up contraction array \sum_j U^{C'}_ij Vjkl

    allocate(cp(vpp%numpw, wup%nums-wup%nums_occ,options%max_i-wup%nums_occ))
    cp(:,:,:)=0.d0


    do iw=1,vpp%numpw_prim
       do ii=1,options%max_i-wup%nums_occ
          do kw=1,vpp%numpw
             cp(kw,vpp%ij(2,iw)-wup%nums_occ,ii)=cp(kw,vpp%ij(2,iw)-wup%nums_occ,ii)+&
                  &dble(wup%umat(ii,vpp%ij(1,iw)))*vpp%vmat(iw,kw)
          enddo
       enddo
    enddo
 
    call free_memory(vpp)

    call read_data_pw_v(vp,options%prefix,options%debug,0,.false.)
    call read_data_pw_ortho_polaw(op,options%prefix)
    call orthonormalize_vpot_para(op,vp)
    call invert_v_pot(vp,vpi)
    call free_memory(vp)
    call invert_ortho_polaw(op,opi)

!loop on negative imaginary times
    if(ionode)  then
       nullify(ww%pw)
       call read_polaw(0,ww,options%debug,options%l_verbose)
       call orthonormalize_inverse(opi,ww)
       allocate(wtemp(ww%numpw,ww%numpw))


      call dgemm('N','N',ww%numpw,ww%numpw,ww%numpw,1.d0,&
                & vpi%vmat,ww%numpw,ww%pw,ww%numpw,0.d0,wtemp,ww%numpw)




      call dgemm('N','N',ww%numpw,ww%numpw,ww%numpw,1.d0,&
                & wtemp,ww%numpw,vpi%vmat,ww%numpw,0.d0,ww%pw,ww%numpw)

       deallocate(wtemp)
       call orthonormalize_inverse(op,ww)
       it=0
       call read_green(it,gg,options%debug,.true.)


       allocate(gf_t(wup%nums-wup%nums_occ,wup%nums-wup%nums_occ))
       do iw=1,(wup%nums-wup%nums_occ)
          do jw=1,(wup%nums-wup%nums_occ)
             gf_t(jw,iw) = gg%gf_p(jw+wup%nums_occ, iw+wup%nums_occ,1)
          enddo
       enddo


       do ii=1,options%max_i-wup%nums_occ




          allocate(qg(op%numpw,wup%nums-wup%nums_occ))
          call dgemm('N','N',op%numpw,wup%nums-wup%nums_occ,wup%nums-wup%nums_occ,1.d0,cp(1,1,ii),&
                 &op%numpw,gf_t,wup%nums-wup%nums_occ,0.d0, qg, op%numpw)

          allocate(pwcp_t(op%numpw,wup%nums-wup%nums_occ))



          call dgemm('N','N',op%numpw,wup%nums-wup%nums_occ,op%numpw,1.d0,ww%pw,op%numpw,&
              &cp(1,1,ii),op%numpw,0.d0, pwcp_t,op%numpw)

          do iw=1,(wup%nums-wup%nums_occ)
             sene(ii) = sene(ii) + ddot(op%numpw,qg(:,iw),1,pwcp_t(:,iw),1)*gg%factor*ww%factor
          enddo

         
          deallocate(pwcp_t)
          deallocate(qg)
          sene(ii)=sene(ii)*(0.d0,1.d0)
       enddo
       deallocate(gf_t)
    endif
    call mp_bcast(sene, ionode_id, world_comm)
    do ii=1,options%max_i-wup%nums_occ
       qp%ene_remainder(ii+wup%nums_occ,1)=qp%ene_remainder(ii+wup%nums_occ,1)-0.5d0*dble(sene(ii))
       write(*,*) 'REMAINDER CONDUCTION', ii, 0.5d0*sene(ii)
    enddo

    call free_memory(vpi)
    call free_memory(op)
    call free_memory(opi)
    call free_memory_polaw(ww)
    call free_memory(wup)
    call free_memory_green(gg)
    deallocate(cp)
    deallocate(sene)
    return
  END SUBROUTINE addconduction_remainder


  SUBROUTINE create_dressed_polarization( options)
!this subroutine calculates the dressed polarization and saves it
!as a polarization with frequency -99999


    USE kinds,         ONLY : DP
    USE io_global,     ONLY : stdout, ionode, ionode_id
    USE input_gw,      ONLY : input_options
    USE basic_structures,  ONLY : v_pot,free_memory, ortho_polaw
    USE green_function,    ONLY : green, read_green, free_memory_green
    USE polarization,      ONLY : polaw, free_memory_polaw, read_polaw, invert_v_pot, invert_ortho_polaw,&
         & orthonormalize_inverse, write_polaw, orthonormalize_vpot
    USE mp,                ONLY : mp_bcast
    USE mp_world,          ONLY : world_comm
    USE para_gww,          ONLY : is_my_pola


    implicit none

    TYPE(input_options) :: options

    TYPE(v_pot) :: vp,vpi
    TYPE(ortho_polaw) :: op,opi
    TYPE(polaw) :: ww!dressed interaction

    INTEGER iw,jw,kw,it,ii
    REAL(kind=DP), ALLOCATABLE :: wtemp(:,:)
    

    nullify(vp%vmat)
    nullify(vpi%vmat)
    nullify(op%on_mat)
    nullify(opi%on_mat)
    nullify(ww%pw)






!read coulombian potential and calculate inverse

    write(stdout,*) 'Routine create_dressed_polarization'


    call read_data_pw_v(vp,options%prefix,options%debug,0,.false.)
    call read_data_pw_ortho_polaw(op,options%prefix)
    call orthonormalize_vpot(op,vp)
    call invert_v_pot(vp,vpi)
    call free_memory(vp)
  
    call invert_ortho_polaw(op,opi)


    if(ionode)  then
       nullify(ww%pw)
       call read_polaw(0,ww,options%debug,options%l_verbose)
       call orthonormalize_inverse(opi,ww)
       allocate(wtemp(ww%numpw,ww%numpw))


      call dgemm('N','N',ww%numpw,ww%numpw,ww%numpw,1.d0,&
                & vpi%vmat,ww%numpw,ww%pw,ww%numpw,0.d0,wtemp,ww%numpw)

      call dgemm('N','N',ww%numpw,ww%numpw,ww%numpw,1.d0,&
                & wtemp,ww%numpw,vpi%vmat,ww%numpw,0.d0,ww%pw,ww%numpw)



       deallocate(wtemp)
       call orthonormalize_inverse(op,ww)
       ww%label=-99999
       call write_polaw(ww,options%debug)
    endif
     
!!!!!!!!!!!
    call free_memory(op)
    call free_memory(opi)
    call free_memory_polaw(ww)

    return
  END SUBROUTINE create_dressed_polarization
