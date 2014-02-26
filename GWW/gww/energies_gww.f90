!
! Copyright (C) 2001-2013 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!


  MODULE energies_gww
!this module contains descriptions and subroutine for final quasi particle energy
! calculations
   USE kinds, ONLY : DP


   TYPE quasi_particles
!energies relative to quasi particles
    INTEGER :: max_i !number of states considered
    INTEGER :: nspin!spin multiplicity
    LOGICAL :: whole_s!if true consider also off diagonal elements
    REAL(kind=DP), POINTER, DIMENSION(:,:) :: ene_dft_ks!kohn sham eigenvalues
    REAL(kind=DP), POINTER, DIMENSION(:,:) :: ene_dft_xc!dft/lda values for exchange and correlation
    REAL(kind=DP), POINTER, DIMENSION(:,:) :: ene_dft_h!dft/lda values for hartree
    COMPLEX(kind=DP), POINTER, DIMENSION(:,:) :: ene_x!effective exchange part
     COMPLEX(kind=DP), POINTER, DIMENSION(:,:) :: ene_h!effective hartree part
    COMPLEX(kind=DP), POINTER, DIMENSION(:,:) :: ene_gw!quasi particle GW eigenenergies
    COMPLEX(kind=DP), POINTER, DIMENSION(:,:) :: ene_gw_pert!perturbative quasi particle GW eigenenergies 
    REAL(kind=DP), POINTER, DIMENSION(:,:) :: ene_hf!perturbative hf energies
    REAL(kind=DP), POINTER, DIMENSION(:,:) :: ene_remainder!for storing remainders
    COMPLEX(kind=DP), POINTER, DIMENSION(:,:) :: ene_gw_off!quasi-particle energies considering the whole Self-energy
    COMPLEX(kind=DP), POINTER, DIMENSION(:,:,:) :: eigen_gw_off!quasi-particle amplitudes considering the whole Self-energy 
    REAL(kind=DP), POINTER, DIMENSION(:,:,:) :: ene_dft_xc_off!off diagonal elements of KS potential
    REAL(kind=DP), POINTER, DIMENSION(:,:,:) :: ene_x_off!off diagonal elements of Fock potential 
 END TYPE quasi_particles

  CONTAINS

    SUBROUTINE write_quasi_particles( qp, options,l_remainder)
!this subroutine write quasi-particles on disk
!ATTENZIONE HF energies not implemented YET
      USE io_global,            ONLY : stdout, ionode
      USE input_gw,             ONLY : input_options
      USE io_files,             ONLY : prefix,tmp_dir

      implicit none

      INTEGER, EXTERNAL :: find_free_unit
      TYPE(quasi_particles) :: qp!object to be written
      TYPE(input_options) :: options!for i/o purposes
      LOGICAL :: l_remainder!if true write also the remainder part

      INTEGER iun
      INTEGER i,j,is

      if(ionode) then
         iun = find_free_unit()
         if(.not. options%debug) then
            open(unit=iun, file=trim(tmp_dir)//trim(prefix)//'-'//'quasi_particles', status='unknown',form='unformatted')
         else
            open(unit=iun, file=trim(tmp_dir)//trim(prefix)//'-'//'quasi_particles', status='unknown',form='formatted')
         endif
         
         if(.not. options%debug) then
            write(iun) qp%max_i
            write(iun) qp%nspin
            write(iun) qp%whole_s
            do is=1,qp%nspin
               write(iun) qp%ene_dft_ks(1:qp%max_i,is)
               write(iun) qp%ene_dft_xc(1:qp%max_i,is)
               write(iun) qp%ene_dft_h(1:qp%max_i,is)
               write(iun) qp%ene_x(1:qp%max_i,is)
               write(iun) qp%ene_h(1:qp%max_i,is)
               write(iun) qp%ene_gw(1:qp%max_i,is)
               write(iun) qp%ene_gw_pert(1:qp%max_i,is)
               write(iun) qp%ene_hf(1:qp%max_i,is)
               if(l_remainder) write(iun) qp%ene_remainder(1:qp%max_i,is)
            enddo
         else
             write(iun,*) qp%max_i
             write(iun,*) qp%nspin
             write(iun,*) qp%whole_s
             do is=1,qp%nspin
                do i=1,qp%max_i
                   write(iun,*) qp%ene_dft_ks(i,is)
                   write(iun,*) qp%ene_dft_xc(i,is)
                   write(iun,*) qp%ene_dft_h(i,is)
                   write(iun,*) qp%ene_x(i,is)
                   write(iun,*) qp%ene_h(i,is)
                   write(iun,*) qp%ene_gw(i,is)
                   write(iun,*) qp%ene_gw_pert(i,is)
                   write(iun,*) qp%ene_hf(i,is)
                   if(l_remainder)write(iun,*) qp%ene_remainder(i,is)
                enddo
             enddo
          endif
          close(iun)
       endif
       return
     END SUBROUTINE write_quasi_particles


    SUBROUTINE read_quasi_particles( qp, options, l_remainder)
!this subroutine write quasi-particles on disk
!HF energies not implemented YET
      USE io_global,            ONLY : stdout, ionode, ionode_id
      USE input_gw,             ONLY : input_options
      USE io_files,             ONLY : prefix,tmp_dir
      USE mp,                   ONLY : mp_bcast
      USE mp_world,             ONLY : world_comm

      implicit none

      INTEGER, EXTERNAL :: find_free_unit
      TYPE(quasi_particles) :: qp!object to be read
      TYPE(input_options) :: options!for i/o purposes
      LOGICAL :: l_remainder !if true read also remainder part

      REAL(kind=DP), ALLOCATABLE :: energies_x(:,:)
      INTEGER iun
      INTEGER i,j,is

      if(ionode) then
         iun = find_free_unit()
         if(.not. options%debug) then
            open(unit=iun, file=trim(tmp_dir)//trim(prefix)//'-'//'quasi_particles', status='old',form='unformatted')
         else
            open(unit=iun, file=trim(tmp_dir)//trim(prefix)//'-'//'quasi_particles', status='old',form='formatted')
         endif

         if(.not. options%debug) then
            read(iun) qp%max_i
            read(iun) qp%nspin
            read(iun) qp%whole_s
         else
            read(iun,*) qp%max_i
            read(iun,*) qp%nspin
            read(iun,*) qp%whole_s
         endif
      endif
      call mp_bcast(qp%max_i, ionode_id,world_comm)
      call mp_bcast(qp%nspin, ionode_id,world_comm)
      call mp_bcast(qp%whole_s, ionode_id,world_comm)
      allocate(qp%ene_dft_ks(qp%max_i,qp%nspin))
      allocate(qp%ene_dft_xc(qp%max_i,qp%nspin))
      allocate(qp%ene_dft_h(qp%max_i,qp%nspin))
      allocate(qp%ene_x(qp%max_i,qp%nspin))
      allocate(qp%ene_h(qp%max_i,qp%nspin))
      allocate(qp%ene_gw(qp%max_i,qp%nspin))
      allocate(qp%ene_gw_pert(qp%max_i,qp%nspin))
      allocate(qp%ene_hf(qp%max_i,qp%nspin))
      if(l_remainder) then
         allocate(qp%ene_remainder(qp%max_i,qp%nspin))
      else
         nullify(qp%ene_remainder)
      endif
      if(ionode) then
         do is=1,qp%nspin
            if(.not.options%debug) then
               read(iun) qp%ene_dft_ks(1:qp%max_i,is)
               read(iun) qp%ene_dft_xc(1:qp%max_i,is)
               read(iun) qp%ene_dft_h(1:qp%max_i,is)
               read(iun) qp%ene_x(1:qp%max_i,is)
               read(iun) qp%ene_h(1:qp%max_i,is)
               read(iun) qp%ene_gw(1:qp%max_i,is)
               read(iun) qp%ene_gw_pert(1:qp%max_i,is)
               read(iun) qp%ene_hf(1:qp%max_i,is)
               if(l_remainder) read(iun) qp%ene_remainder(1:qp%max_i,is)
            else
               do i=1,qp%max_i
                  read(iun,*) qp%ene_dft_ks(i,is)
                  read(iun,*) qp%ene_dft_xc(i,is)
                  read(iun,*) qp%ene_dft_h(i,is)
                  read(iun,*) qp%ene_x(i,is)
                  read(iun,*) qp%ene_h(i,is)
                  read(iun,*) qp%ene_gw(i,is)
                  read(iun,*) qp%ene_gw_pert(i,is)
                  read(iun,*) qp%ene_hf(i,is)
                  if(l_remainder) read(iun,*) qp%ene_remainder(i,is)
               enddo
            endif
         enddo
         close(iun)
      endif
      call mp_bcast(qp%ene_dft_ks(:,:),ionode_id,world_comm)
      call mp_bcast(qp%ene_dft_xc(:,:),ionode_id,world_comm)
      call mp_bcast(qp%ene_dft_h(:,:),ionode_id,world_comm)
      call mp_bcast(qp%ene_x(:,:),ionode_id,world_comm)
      call mp_bcast(qp%ene_h(:,:),ionode_id,world_comm)
      call mp_bcast(qp%ene_gw(:,:),ionode_id,world_comm)
      call mp_bcast(qp%ene_gw_pert(:,:),ionode_id,world_comm)
      call mp_bcast(qp%ene_hf(:,:),ionode_id,world_comm)
      if(l_remainder) call mp_bcast(qp%ene_remainder(:,:),ionode_id,world_comm)

!if required re-read exchange energies

      if(options%l_read_exchange) then
          allocate(energies_x(qp%max_i,qp%nspin))
          call read_data_pw_exchange(energies_x,qp%max_i,options%prefix,qp%nspin)

          qp%ene_hf(1:qp%max_i,1:qp%nspin)=qp%ene_hf(1:qp%max_i,1:qp%nspin)-qp%ene_x(1:qp%max_i,1:qp%nspin)&
               &+energies_x(1:qp%max_i,1:qp%nspin)
          qp%ene_gw(1:qp%max_i,1:qp%nspin)=qp%ene_gw(1:qp%max_i,1:qp%nspin)-qp%ene_x(1:qp%max_i,1:qp%nspin)&
               &+energies_x(1:qp%max_i,1:qp%nspin)
          qp%ene_gw_pert(1:qp%max_i,1:qp%nspin)=qp%ene_gw_pert(1:qp%max_i,1:qp%nspin)&
               &-qp%ene_x(1:qp%max_i,1:qp%nspin)+energies_x(1:qp%max_i,1:qp%nspin)
          qp%ene_x(1:qp%max_i,1:qp%nspin)=cmplx(energies_x(1:qp%max_i,1:qp%nspin),0.d0)

          deallocate(energies_x)
       endif
!if required re-read XC LDA energies
       if(options%l_dft_xc_file) then
          call read_data_pw_dft_xc(qp%ene_dft_xc(:,1),qp%max_i,options%prefix)
       endif

      return
    END SUBROUTINE read_quasi_particles


    SUBROUTINE initialize_quasi_particle(qp)
!this subroutine nullify all arrays

      implicit none
      
      TYPE(quasi_particles) :: qp
      
      nullify(qp%ene_dft_ks)
      nullify(qp%ene_dft_xc)
      nullify(qp%ene_dft_h)
      nullify(qp%ene_x)
      nullify(qp%ene_h)
      nullify(qp%ene_gw)
      nullify(qp%ene_gw_pert)
      nullify(qp%ene_hf)
      nullify(qp%ene_remainder)
      nullify(qp%ene_gw_off)
      nullify(qp%eigen_gw_off)
      nullify(qp%ene_dft_xc_off)
      nullify(qp%ene_x_off)

      return
    END SUBROUTINE initialize_quasi_particle



  SUBROUTINE free_memory_quasi_particles(qp)
!deallocates if allocated
    implicit none

    TYPE(quasi_particles) :: qp

    if(associated(qp%ene_dft_ks)) then
      deallocate(qp%ene_dft_ks)
      nullify(qp%ene_dft_ks)
    endif

   if(associated(qp%ene_dft_xc)) then
      deallocate(qp%ene_dft_xc)
      nullify(qp%ene_dft_xc)
    endif

    if(associated(qp%ene_dft_h)) then
      deallocate(qp%ene_dft_h)
      nullify(qp%ene_dft_h)
    endif


   if(associated(qp%ene_x)) then
      deallocate(qp%ene_x)
      nullify(qp%ene_x)
    endif

   if(associated(qp%ene_h)) then
      deallocate(qp%ene_h)
      nullify(qp%ene_h)
    endif


   if(associated(qp%ene_gw)) then
      deallocate(qp%ene_gw)
      nullify(qp%ene_gw)
    endif

    if(associated(qp%ene_gw_pert)) then
      deallocate(qp%ene_gw_pert)
      nullify(qp%ene_gw_pert)
    endif


    if(associated(qp%ene_hf)) then
      deallocate(qp%ene_hf)
      nullify(qp%ene_hf)
    endif

    if(associated(qp%ene_remainder)) then
       deallocate(qp%ene_remainder)
       nullify(qp%ene_remainder)
    endif

    if(associated(qp%ene_gw_off)) then
       deallocate(qp%ene_gw_off)
       nullify(qp%ene_gw_off)
    endif

    if(associated(qp%eigen_gw_off)) then
       deallocate(qp%eigen_gw_off)
       nullify(qp%eigen_gw_off)
    endif

  if(associated(qp%ene_dft_xc_off)) then
       deallocate(qp%ene_dft_xc_off)
       nullify(qp%ene_dft_xc_off)
    endif

  if(associated(qp%ene_x_off)) then
       deallocate(qp%ene_x_off)
       nullify(qp%ene_x_off)
    endif

   return
 END SUBROUTINE free_memory_quasi_particles


 SUBROUTINE printout_quasi_particles(qp)
!this subroutine prints out the lda and gw energies

  USE io_global, ONLY : stdout, ionode
  USE constants, ONLY : RYTOEV
  USE mp,        ONLY : mp_barrier
    USE io_files,             ONLY : prefix,tmp_dir
 

  implicit none
  INTEGER, EXTERNAL :: find_free_unit
  TYPE(quasi_particles) :: qp

  INTEGER ::  ii,iun,is

  if(ionode) then
     do is=1,qp%nspin
        write(stdout,*) 'QUASI-PARTICLES ENERGIES IN Ev, Spin:', is, qp%nspin
        do ii=1,qp%max_i
           write(stdout,'(''State:'',i5,''DFT  :'',f10.5,'' GW-PERT  :'',f10.5,'' GW  :'',f10.5, '' HF-pert :'',f10.5)') &
                & ii,qp%ene_dft_ks(ii,is)*RYTOEV, real(qp%ene_gw_pert(ii,is))*RYTOEV, &
                & real(qp%ene_gw(ii,is))*RYTOEV,qp%ene_hf(ii,is)*RYTOEV
        enddo
        write(stdout,*) 'IMAGINARY ENERGIES IN Ev:'
        do ii=1,qp%max_i
           write(stdout,'(''State:'',i5,'' GW (Im) :'',f10.5)') ii,aimag(qp%ene_gw(ii,is))*RYTOEV
        enddo
     enddo
!write bands.dat file

     iun =  find_free_unit()
     open( unit=iun, file=trim(tmp_dir)//trim(prefix)//'-'//'bands.dat', status='unknown',form='formatted')
     write(iun,'(i8)') qp%max_i
     write(iun,'(i8)') qp%nspin
     do is=1,qp%nspin
        do ii=1,qp%max_i
           write(iun,'(i5,4f10.5)') ii,qp%ene_dft_ks(ii,is)*RYTOEV, real(qp%ene_gw_pert(ii,is))*RYTOEV, &
                & real(qp%ene_gw(ii,is))*RYTOEV,qp%ene_hf(ii,is)*RYTOEV
        enddo
     enddo
     close(iun)
       


  endif
  



  return

  END SUBROUTINE
  

 SUBROUTINE printout_quasi_particles_off(qp)
!this subroutine prints out the lda and gw energies
!where the whole self energy matrix has been considered

  USE io_global, ONLY : stdout, ionode
  USE constants, ONLY : RYTOEV
  USE mp,        ONLY : mp_barrier
 
  implicit none
  INTEGER, EXTERNAL :: find_free_unit
  TYPE(quasi_particles) :: qp

  INTEGER ::  ii,iun,is


  if(ionode) then
     if(qp%whole_s) then
        write(stdout,*) 'RESULTS FROM WHOLE SE MATRIX:'
        do is=1,qp%nspin
           write(stdout,*) 'QUASI-PARTICLES ENERGIES IN Ev, Spin:', is, qp%nspin
           do ii=1,qp%max_i
              write(stdout,'(''State:'',i5,''DFT  :'',f10.5,'' GW  :'',f10.5, '' HF-pert :'',f10.5)') &
                   & ii,qp%ene_dft_ks(ii,is)*RYTOEV, &
                   & real(qp%ene_gw_off(ii,is))*RYTOEV,qp%ene_hf(ii,is)*RYTOEV
           enddo
           write(stdout,*) 'IMAGINARY ENERGIES IN Ev:'
           do ii=1,qp%max_i
              write(stdout,'(''State:'',i5,'' GW (Im) :'',f10.5)') ii,aimag(qp%ene_gw_off(ii,is))*RYTOEV
           enddo
        enddo
     else
        write(stdout,*) 'OFF DIAGONAL ELEMENTS OF SE NOT AVAILABLE'
     endif
  endif

  return

END SUBROUTINE printout_quasi_particles_off





 END MODULE energies_gww


 
