!
! Copyright (C) 2001-2013 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!

!this modules contain the routine for the contour integration

MODULE contour

  USE kinds, ONLY  : DP


  TYPE w_expectation
!descriptor for <Psi_i|Psi_j(r)Psi_j(r')W(r,r')|Psi_i>
     INTEGER ::n!number of steps
     INTEGER :: max_i!number of states considered  
     INTEGER :: i_min!minimum state to be calculated
     INTEGER :: i_max!maximum state to be calculated
     INTEGER :: nspin!spin multiplicity
     COMPLEX(kind=DP), DIMENSION(:), POINTER :: grid
     COMPLEX(kind=DP), DIMENSION(:,:,:,:), POINTER :: diag(:,:,:,:)
  END TYPE w_expectation

  TYPE w_poles
!descriptor for the analytic continuation of <Psi_i|Psi_j(r)Psi_j(r')W(r,r')|Psi_i> 
     INTEGER :: max_i!number of states considered  
     INTEGER :: i_min!minimum state to be calculated
     INTEGER :: i_max!maximum state to be calculated
     INTEGER :: nspin!spin multiplicity                
     INTEGER :: n_multipoles!number of multipoles considered
     COMPLEX(kind=DP), DIMENSION(:,:,:), POINTER :: a_0!parameters a_0 
     COMPLEX(kind=DP), DIMENSION(:,:,:,:), POINTER :: a!parameters a (n_multipoles,max_i(poles),max_i(KS states),nspin)
     COMPLEX(kind=DP), DIMENSION(:,:,:,:), POINTER :: b!parameters b (n_multipoles,max_i(poles),max_i(KS states),nspin)     
  END  TYPE w_poles

  CONTAINS
    SUBROUTINE initialize_w_expectation(we)
      implicit none
      TYPE(w_expectation) :: we
      nullify(we%grid)
      nullify(we%diag)
      return
    END SUBROUTINE initialize_w_expectation

    SUBROUTINE initialize_w_poles(wp)
      implicit none
      TYPE(w_poles) :: wp
      nullify(wp%a_0)
      nullify(wp%a)
      nullify(wp%b)
      return
    END SUBROUTINE initialize_w_poles

    SUBROUTINE free_memory_w_poles(wp)
      implicit none
      TYPE(w_poles) :: wp
      if(associated(wp%a_0)) deallocate (wp%a_0)
      nullify(wp%a_0)
      if(associated(wp%a)) deallocate (wp%a)
      nullify(wp%a)
      if(associated(wp%b)) deallocate (wp%b)
      nullify(wp%b)
      
      return
    END SUBROUTINE free_memory_w_poles

    SUBROUTINE free_memory_w_expectation(we)
      implicit none
      TYPE(w_expectation) :: we
      if(associated(we%grid)) deallocate(we%grid)
      nullify(we%grid)
      if(associated(we%diag)) deallocate(we%diag)
      nullify(we%diag)
      return
    END SUBROUTINE free_memory_w_expectation

  SUBROUTINE write_w_poles(wp)
      USE io_global,          ONLY : stdout, ionode
      USE input_gw,           ONLY : input_options
      USE mp,                 ONLY : mp_barrier
     USE io_files,             ONLY : prefix,tmp_dir
      implicit none
      INTEGER, EXTERNAL :: find_free_unit
      TYPE(w_poles) :: wp!the structure to be written 
      INTEGER :: iun,is

      if(ionode) then
         iun = find_free_unit()
         open( unit=iun, file=trim(tmp_dir)//trim(prefix)//'-'//'wpoles', status='unknown',form='unformatted')
         write(iun) wp%max_i
         write(iun) wp%i_min
         write(iun) wp%i_max
         write(iun) wp%nspin
         write(iun) wp%n_multipoles
         write(iun) wp%a_0(1:wp%max_i,1:wp%max_i,1:wp%nspin)
         write(iun) wp%a(1:wp%n_multipoles,1:wp%max_i,1:wp%max_i,1:wp%nspin)
         write(iun) wp%b(1:wp%n_multipoles,1:wp%max_i,1:wp%max_i,1:wp%nspin)
         close(iun)
        endif
      return
    END SUBROUTINE write_w_poles

    SUBROUTINE read_w_poles(wp)
      USE io_global,          ONLY : stdout, ionode,ionode_id
      USE input_gw,           ONLY : input_options
      USE mp,                 ONLY : mp_bcast
      USE mp_world,           ONLY : world_comm
     USE io_files,             ONLY : prefix,tmp_dir
      implicit none
      INTEGER, EXTERNAL :: find_free_unit
      TYPE(w_poles) :: wp!the structure to be read
      INTEGER :: iun,is

      if(ionode) then
         iun = find_free_unit()
         open( unit=iun, file=trim(tmp_dir)//trim(prefix)//'-'//'wpoles', status='old',form='unformatted')
         read(iun) wp%max_i
         read(iun) wp%i_min
         read(iun) wp%i_max
         read(iun) wp%nspin
         read(iun) wp%n_multipoles
      endif
      call mp_bcast(wp%max_i,ionode_id,world_comm)
      call mp_bcast(wp%i_min,ionode_id,world_comm)
      call mp_bcast(wp%i_max,ionode_id,world_comm)
      call mp_bcast(wp%nspin,ionode_id,world_comm)
      call mp_bcast(wp%n_multipoles,ionode_id,world_comm)
      allocate(wp%a_0(wp%max_i,wp%max_i,wp%nspin))
      allocate(wp%a(wp%n_multipoles,wp%max_i,wp%max_i,wp%nspin))
      allocate(wp%b(wp%n_multipoles,wp%max_i,wp%max_i,wp%nspin))
      if(ionode) then
         read(iun) wp%a_0(1:wp%max_i,1:wp%max_i,1:wp%nspin)
         read(iun) wp%a(1:wp%n_multipoles,1:wp%max_i,1:wp%max_i,1:wp%nspin)
         read(iun) wp%b(1:wp%n_multipoles,1:wp%max_i,1:wp%max_i,1:wp%nspin)
         close(iun)
      endif
      call mp_bcast(wp%a_0,ionode_id,world_comm)
      call mp_bcast(wp%a,ionode_id,world_comm)
      call mp_bcast(wp%b,ionode_id,world_comm)

      return
    END SUBROUTINE read_w_poles

    SUBROUTINE write_w_expectation(we)
      USE io_global,          ONLY : stdout, ionode
      USE input_gw,           ONLY : input_options
      USE mp,                 ONLY : mp_barrier
      USE io_files,  ONLY : prefix, tmp_dir
      
      implicit none
      INTEGER, EXTERNAL :: find_free_unit
      TYPE(w_expectation) :: we!the structure to be written              
      
      INTEGER :: iun,is

      if(ionode) then
         iun = find_free_unit()
         open( unit=iun, file=trim(tmp_dir)//trim(prefix)//'-'//'wexpectation', status='unknown',form='unformatted')
         write(iun) we%n
         write(iun) we%max_i
         write(iun) we%i_min
         write(iun) we%i_max
         write(iun) we%nspin
         write(iun) we%grid(1:we%n)
         do is=1,we%nspin
            write(iun) we%diag(1:we%n,1:we%max_i, 1:we%max_i,is)
         enddo
         close(iun)
      endif
      return
    END SUBROUTINE write_w_expectation

    SUBROUTINE read_w_expectation(we)
      USE io_global,          ONLY : stdout, ionode,ionode_id
      USE input_gw,           ONLY : input_options
      USE mp,                 ONLY : mp_bcast
      USE mp_world,           ONLY : world_comm
      USE io_files,           ONLY : prefix, tmp_dir

      implicit none
      INTEGER, EXTERNAL :: find_free_unit
      TYPE(w_expectation),INTENT(out) :: we!the structure to be written                                                                                                 

      INTEGER :: iun,is

      if(ionode) then
         iun = find_free_unit()
         open( unit=iun, file=trim(tmp_dir)//trim(prefix)//'-'//'wexpectation', status='old',form='unformatted')
         read(iun) we%n
         read(iun) we%max_i
         read(iun) we%i_min
         read(iun) we%i_max
         read(iun) we%nspin
      endif
      call mp_bcast(we%n,ionode_id,world_comm)
      call mp_bcast(we%max_i,ionode_id,world_comm)
      call mp_bcast(we%i_min,ionode_id,world_comm)
      call mp_bcast(we%i_max,ionode_id,world_comm)
      call mp_bcast(we%nspin,ionode_id,world_comm)
      allocate(we%grid(we%n),we%diag(we%n,we%max_i,we%max_i,we%nspin))
      if(ionode) then
         read(iun) we%grid(1:we%n)
         do is=1,we%nspin
            read(iun) we%diag(1:we%n,1:we%max_i, 1:we%max_i,is)
         enddo
         close(iun)
      endif
      call mp_bcast(we%grid,ionode_id,world_comm)
      call mp_bcast(we%diag,ionode_id,world_comm)
      return


    END SUBROUTINE read_w_expectation
    
    SUBROUTINE create_w_expectation(we, tf, options)
!this subroutine create the diagonal elements of W starting from Pgreek read from disk

      USE kinds,             ONLY : DP
      USE io_global,         ONLY : stdout, ionode, ionode_id
      USE input_gw,          ONLY : input_options
      USE basic_structures,  ONLY : v_pot,wannier_u,free_memory, initialize_memory,lanczos_chain, vt_mat_lanczos,tt_mat_lanczos,&
           & contour_terms
      USE green_function,    ONLY : green, read_green, free_memory_green, initialize_green
      USE polarization,      ONLY : polaw, free_memory_polaw, read_polaw, write_polaw,invert_v_pot, initialize_polaw, &
           & read_polaw_global
      USE mp,                ONLY : mp_sum, mp_bcast
      USE mp_world,          ONLY : nproc,mpime,world_comm
      USE times_gw,          ONLY : times_freqs
      USE self_energy_storage, ONLY : self_storage,write_self_storage_ondisk,free_memory_self_storage
      USE lanczos
      USE constants,          ONLY : tpi,pi
      USE para_gww,           ONLY : is_my_pola
      
  

      implicit none
    
      TYPE(w_expectation) :: we
      TYPE(times_freqs), INTENT(in) :: tf!for time frequency grids                                                   
      TYPE(input_options) :: options
    

      TYPE(polaw) :: pp!dressed polarization       
      TYPE(tt_mat_lanczos) :: sg!overlap <s_global|s_local>
      TYPE(vt_mat_lanczos) :: sl!overlap <s_local|v(\Phi)>
      TYPE(contour_terms) :: ct

      INTEGER :: iw, ii,is,jj
      REAL(kind=DP), ALLOCATABLE :: cs_mat(:,:),tmp_mat(:,:),tmp_mat2(:,:),tmp_mat3(:,:)

      INTEGER :: jmin, jmax

      call initialize_memory(ct)
      call initialize_memory(sg)
      call initialize_memory(sl)
      call initialize_polaw(pp)
!allocate
      we%nspin=options%nspin
      we%n=options%n+1
      we%max_i=options%max_i
      we%i_min=options%i_min
      we%i_max=options%i_max
      allocate(we%grid(we%n))
      allocate(we%diag(we%n,we%max_i,we%max_i,we%nspin))
      we%diag(:,:,:,:)=(0.d0,0.d0)
!loop on spin
      if(.not.options%l_big_system) then
         jmin=1
         jmax=options%i_max
      else
         jmin=options%i_min
         jmax=options%i_max
      endif
      do is=1,we%nspin
!read in contour terms
         if(.not.options%l_big_system) call read_data_pw_contour(ct,options%prefix,is,1)
!loop on KS states , poles
         do jj=jmin,jmax
            if(options%l_big_system) call read_data_pw_contour(ct,options%prefix,is,jj)
            call  read_data_pw_tt_mat_lanczos(sg, jj, options%prefix,.false.,is)
            call  read_data_pw_vt_mat_lanczos(sl, jj, options%prefix,.false., is)

            allocate(cs_mat(ct%nums,sl%numpw))
            allocate(tmp_mat(ct%nums,sg%numl)) 
            call dgemm('T','N',ct%nums,sg%numl,sg%numt,1.d0,ct%cmat,ct%numt,sg%tt_mat,sg%numt,0.d0,tmp_mat,ct%nums)
            call dgemm('N','T',ct%nums,sl%numpw,sl%numl,1.d0,tmp_mat,ct%nums,sl%vt_mat,sl%numpw,0.d0,cs_mat,ct%nums)
            allocate(tmp_mat2(ct%nums,sl%numpw),tmp_mat3(ct%nums,ct%nums))
!loop on frequencies
            do iw=0,options%n
               if(is_my_pola(iw)) then
                  call read_polaw(iw,pp,options%debug,options%l_verbose)
                  call dgemm('N','N',ct%nums,pp%numpw,pp%numpw,1.d0,cs_mat,ct%nums,pp%pw,pp%numpw,0.d0,tmp_mat2,ct%nums)
!also off-diagonal elements are calculated although not necessary
                  call dgemm('N','T',ct%nums,ct%nums,pp%numpw,1.d0,tmp_mat2,ct%nums,cs_mat,ct%nums,0.d0,tmp_mat3,ct%nums)
                  do ii=1,ct%nums
                     we%diag(iw+1,jj,ii,is)=tmp_mat3(ii,ii) !GIUSTO CUSSI
                  enddo
                  call free_memory_polaw(pp)
               endif
            enddo
            call free_memory(sg)
            call free_memory(sl)
            deallocate(cs_mat,tmp_mat)
            deallocate(tmp_mat2,tmp_mat3)
            if(options%l_big_system) call free_memory(ct)
         enddo
         if(.not.options%l_big_system) call free_memory(ct)
      enddo
      call mp_sum(we%diag,world_comm)
!now set up frequency grid
      we%grid(1:we%n)=dcmplx(0.d0,tf%freqs(0:tf%n))

      call free_memory_polaw(pp)
      call free_memory(ct)
      call free_memory(sg)
      call free_memory(sl)
      return
    END SUBROUTINE create_w_expectation

    SUBROUTINE create_w_poles(we,wp,options)
!this subroutine perform non-linar fits for finding the expansion
!of the terms for the contour integration
!the sum over poles is distributed among processors 
     
      USE io_global,  ONLY : stdout
      USE input_gw,   ONLY : input_options
      USE para_gww,   ONLY : is_my_state_range
      USE mp,         ONLY : mp_sum
      USE mp_world,   ONLY : world_comm
 
    implicit none

    TYPE(w_expectation), INTENT(in) :: we!data on imaginary frequency
    TYPE(w_poles), INTENT(out) :: wp!poles to be found
    TYPE(input_options), INTENT(in) :: options

    INTEGER :: ii,jj, kk,is,mm
    COMPLEX(kind=DP), ALLOCATABLE :: z(:),s(:)
    REAL(kind=DP) :: chi, chi0
    COMPLEX(kind=DP) :: a_0_old
    COMPLEX(kind=DP), ALLOCATABLE :: a_old(:), b_old(:)


!set up wp 
    call initialize_w_poles(wp)

    wp%nspin=we%nspin
    wp%max_i=we%max_i
    wp%i_max=we%i_max
    wp%i_min=we%i_min
    wp%n_multipoles=options%n_multipoles

!allocate
    allocate(wp%a_0(wp%max_i,wp%max_i,wp%nspin))
    allocate(wp%a(wp%n_multipoles,wp%max_i,wp%max_i,wp%nspin))     
    allocate(wp%b(wp%n_multipoles,wp%max_i,wp%max_i,wp%nspin))
    allocate(a_old(wp%n_multipoles),b_old(wp%n_multipoles))
    wp%a_0(:,:,:)=(0.d0,0.d0)
    wp%a(:,:,:,:)=(0.d0,0.d0)
    wp%b(:,:,:,:)=(0.d0,0.d0)
    
    allocate(z(we%n),s(we%n))
    z(1:we%n)=we%grid(1:we%n)

    do is=1,wp%nspin
       do ii=wp%i_min,wp%i_max
         ! do jj=1,wp%max_i
           do jj=wp%i_min,wp%i_max!ATTENZIONE
             if(is_my_state_range(jj)) then
                s(1:we%n)=we%diag(1:we%n,jj,ii,is)
                wp%a_0(jj,ii,is)=(0.0,0.0d0)
                do mm=1,options%n_multipoles
                   wp%a(mm,jj,ii,is)=cmplx(real(mm)*(0.01d0),0.d0)
                   wp%b(mm,jj,ii,is)=cmplx((0.5d0)*real(mm)*(-1.d0)**real(mm),-0.01d0)
                enddo
                write(stdout,*) 'Call fit_multipole'
                FLUSH(stdout)
                
                call fit_multipole(we%n,wp%n_multipoles,z,s,wp%a_0(jj,ii,is),&
                     &wp%a(:,jj,ii,is),wp%b(:,jj,ii,is),1.d0,options%fit_thres,options%fit_maxiter)
                write(stdout,*) 'Done'
                FLUSH(stdout)

                a_0_old=wp%a_0(jj,ii,is)
                do mm=1,wp%n_multipoles
                   a_old(mm)=wp%a(mm,jj,ii,is)
                   b_old(mm)=wp%b(mm,jj,ii,is)
                enddo


                FLUSH(stdout)
                call fit_multipole_minpack(we%n,wp%n_multipoles,z,s,wp%a_0(jj,ii,is),&
                        &wp%a(:,jj,ii,is),wp%b(:,jj,ii,is),options%fit_thres, options%n_max_minpack, chi)
     
                write(stdout,*) 'FIT pole :', jj, ii,is
                write(stdout,*) 'FIT    a_0:', wp%a_0(jj,ii,is)
                do mm=1,wp%n_multipoles
                   write(stdout,*) 'FIT    a:',mm,wp%a(mm,jj,ii,is)
                   write(stdout,*) 'FIT    b:',mm,wp%b(mm,jj,ii,is)
                enddo
                FLUSH(stdout)
          
                
             endif
          enddo
          call mp_sum(wp%a_0(:,ii,is),world_comm)
          call mp_sum(wp%a(:,:,ii,is),world_comm)
          call mp_sum(wp%b(:,:,ii,is),world_comm)

       enddo
    enddo
    deallocate(a_old,b_old)
    deallocate(z,s)
    return
  END SUBROUTINE create_w_poles

  FUNCTION w_poles_value(energy,wp,jj,ii,ispin)

    implicit none
    COMPLEX(kind=DP) :: w_poles_value!the value of the pole jj, for the KS states ii with spin is
    COMPLEX(kind=DP), INTENT(in) :: energy !frequency considered
    TYPE(w_poles), INTENT(in) :: wp
    INTEGER, INTENT(in) :: jj
    INTEGER, INTENT(in) :: ii
    INTEGER, INTENT(in) :: ispin


    COMPLEX(kind=DP) :: fz
    INTEGER :: ip
  
   
    
    if(dble(energy) >= 0.d0 ) then
       fz=wp%a_0(jj,ii,ispin)
       do ip=1,wp%n_multipoles
          fz=fz+wp%a(ip,jj,ii,ispin)/(energy-wp%b(ip,jj,ii,ispin))
       enddo
    else
       fz=conjg(wp%a_0(jj,ii,ispin))
       do ip=1,wp%n_multipoles
          fz=fz+conjg(wp%a(ip,jj,ii,ispin))/(energy-conjg(wp%b(ip,jj,ii,ispin)))
       enddo
    endif
    
    w_poles_value=fz

    return

  END FUNCTION w_poles_value

  END MODULE contour
