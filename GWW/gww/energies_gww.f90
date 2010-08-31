!P.Umari Program GWW

  MODULE energies_gww
!this module contains descriptions and subroutine for final quasi particle energy
! calculations
   USE kinds, ONLY : DP


   TYPE quasi_particles
!energies realtive to quasi particles
    INTEGER :: max_i !number of states considered
    REAL(kind=DP), POINTER, DIMENSION(:) :: ene_dft_ks!kohn sham eigenvalues
    REAL(kind=DP), POINTER, DIMENSION(:) :: ene_dft_xc!dft/lda values for exchange and correlation
    REAL(kind=DP), POINTER, DIMENSION(:) :: ene_dft_h!dft/lda values for hartree
    COMPLEX(kind=DP), POINTER, DIMENSION(:) :: ene_x!effective exchange part
     COMPLEX(kind=DP), POINTER, DIMENSION(:) :: ene_h!effective hartree part
    COMPLEX(kind=DP), POINTER, DIMENSION(:) :: ene_gw!quasi particle GW eigenenergies
    COMPLEX(kind=DP), POINTER, DIMENSION(:) :: ene_gw_pert!perturbative quasi particle GW eigenenergies
    REAL(kind=DP), POINTER, DIMENSION(:) :: ene_hf!perturbative hf energies
    REAL(kind=DP), POINTER, DIMENSION(:) :: ene_remainder!for storing remainders
 END TYPE quasi_particles

  CONTAINS

    SUBROUTINE write_quasi_particles( qp, options,l_remainder)
!this subroutine write quasi-particles on disk
!ATTENZIONE HF energies not implemented YET
      USE io_global,            ONLY : stdout, ionode
      USE input_gw,             ONLY : input_options
      USE io_files,             ONLY : find_free_unit, prefix

      implicit none

      TYPE(quasi_particles) :: qp!object to be written
      TYPE(input_options) :: options!for i/o purposes
      LOGICAL :: l_remainder!if true write also the remainder part

      INTEGER iun
      INTEGER i,j

      if(ionode) then
         iun = find_free_unit()
         if(.not. options%debug) then
            open(unit=iun, file='quasi_particles', status='unknown',form='unformatted')
         else
            open(unit=iun, file='quasi_particles', status='unknown',form='formatted')
         endif

         if(.not. options%debug) then
            write(iun) qp%max_i
            write(iun) qp%ene_dft_ks(1:qp%max_i)
            write(iun) qp%ene_dft_xc(1:qp%max_i)
            write(iun) qp%ene_dft_h(1:qp%max_i)
            write(iun) qp%ene_x(1:qp%max_i)
            write(iun) qp%ene_h(1:qp%max_i)
            write(iun) qp%ene_gw(1:qp%max_i)
            write(iun) qp%ene_gw_pert(1:qp%max_i)
            write(iun) qp%ene_hf(1:qp%max_i)
            if(l_remainder) write(iun) qp%ene_remainder(1:qp%max_i)
         else
             write(iun,*) qp%max_i
             do i=1,qp%max_i
                write(iun,*) qp%ene_dft_ks(i)
                write(iun,*) qp%ene_dft_xc(i)
                write(iun,*) qp%ene_dft_h(i)
                write(iun,*) qp%ene_x(i)
                write(iun,*) qp%ene_h(i)
                write(iun,*) qp%ene_gw(i)
                write(iun,*) qp%ene_gw_pert(i)
                write(iun,*) qp%ene_hf(i)
                if(l_remainder)write(iun,*) qp%ene_remainder(i)
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
      USE io_files,             ONLY : find_free_unit, prefix
      USE mp,                   ONLY : mp_bcast

      implicit none

      TYPE(quasi_particles) :: qp!object to be read
      TYPE(input_options) :: options!for i/o purposes
      LOGICAL :: l_remainder !if true read also remainder part

      REAL(kind=DP), ALLOCATABLE :: energies_x(:)
      INTEGER iun
      INTEGER i,j

      if(ionode) then
         iun = find_free_unit()
         if(.not. options%debug) then
            open(unit=iun, file='quasi_particles', status='old',form='unformatted')
         else
            open(unit=iun, file='quasi_particles', status='old',form='formatted')
         endif

         if(.not. options%debug) then
            read(iun) qp%max_i
         else
            read(iun,*) qp%max_i
         endif
      endif
      call mp_bcast(qp%max_i, ionode_id)
      allocate(qp%ene_dft_ks(qp%max_i))
      allocate(qp%ene_dft_xc(qp%max_i))
      allocate(qp%ene_dft_h(qp%max_i))
      allocate(qp%ene_x(qp%max_i))
      allocate(qp%ene_h(qp%max_i))
      allocate(qp%ene_gw(qp%max_i))
      allocate(qp%ene_gw_pert(qp%max_i))
      allocate(qp%ene_hf(qp%max_i))
      if(l_remainder) then
         allocate(qp%ene_remainder(qp%max_i))
      else
         nullify(qp%ene_remainder)
      endif
      if(ionode) then
         if(.not.options%debug) then
            read(iun) qp%ene_dft_ks(1:qp%max_i)
            read(iun) qp%ene_dft_xc(1:qp%max_i)
            read(iun) qp%ene_dft_h(1:qp%max_i)
            read(iun) qp%ene_x(1:qp%max_i)
            read(iun) qp%ene_h(1:qp%max_i)
            read(iun) qp%ene_gw(1:qp%max_i)
            read(iun) qp%ene_gw_pert(1:qp%max_i)
            read(iun) qp%ene_hf(1:qp%max_i)
            if(l_remainder) read(iun) qp%ene_remainder(1:qp%max_i)
         else
            do i=1,qp%max_i
               read(iun,*) qp%ene_dft_ks(i)
               read(iun,*) qp%ene_dft_xc(i)
               read(iun,*) qp%ene_dft_h(i)
               read(iun,*) qp%ene_x(i)
                read(iun,*) qp%ene_h(i)
               read(iun,*) qp%ene_gw(i)
               read(iun,*) qp%ene_gw_pert(i)
               read(iun,*) qp%ene_hf(i)
               if(l_remainder) read(iun,*) qp%ene_remainder(i)
            enddo
         endif
         close(iun)
      endif
      call mp_bcast(qp%ene_dft_ks(:),ionode_id)
      call mp_bcast(qp%ene_dft_xc(:),ionode_id)
      call mp_bcast(qp%ene_dft_h(:),ionode_id)
      call mp_bcast(qp%ene_x(:),ionode_id)
      call mp_bcast(qp%ene_h(:),ionode_id)
      call mp_bcast(qp%ene_gw(:),ionode_id)
      call mp_bcast(qp%ene_gw_pert(:),ionode_id)
      call mp_bcast(qp%ene_hf(:),ionode_id)
      if(l_remainder) call mp_bcast(qp%ene_remainder(:),ionode_id)

!if required re-read exchange energies

      if(options%l_read_exchange) then
          allocate(energies_x(qp%max_i))
          call read_data_pw_exchange(energies_x,qp%max_i,options%prefix)

          qp%ene_hf(:)=qp%ene_hf(:)-qp%ene_x(:)+energies_x(:)
          qp%ene_gw(:)=qp%ene_gw(:)-qp%ene_x(:)+energies_x(:)
          qp%ene_gw_pert(:)=qp%ene_gw_pert(:)-qp%ene_x(:)+energies_x(:)
          qp%ene_x(:)=cmplx(energies_x(:),0.d0)

          deallocate(energies_x)
       endif
!if required re-read XC LDA energies
       if(options%l_dft_xc_file) then
          call read_data_pw_dft_xc(qp%ene_dft_xc,qp%max_i,options%prefix)
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


   return
 END SUBROUTINE free_memory_quasi_particles


 SUBROUTINE printout_quasi_particles(qp)
!this subroutine prints out the lda and gw energies

  USE io_global, ONLY : stdout, ionode
  USE constants, ONLY : RYTOEV
  USE mp,        ONLY : mp_barrier
  USE io_files,  ONLY : find_free_unit, prefix


  implicit none

  TYPE(quasi_particles) :: qp

  INTEGER ::  ii, iun

  if(ionode) then
     write(stdout,*) 'QUASI-PARTICLES ENERGIES IN Ev:'
     do ii=1,qp%max_i
        write(stdout,'(''State:'',i5,''LDA  :'',f10.5,'' GW-PERT  :'',f10.5,'' GW  :'',f10.5, '' HF-pert :'',f10.5)') &
             & ii,qp%ene_dft_ks(ii)*RYTOEV, real(qp%ene_gw_pert(ii))*RYTOEV, &
             & real(qp%ene_gw(ii))*RYTOEV,qp%ene_hf(ii)*RYTOEV
     enddo
     write(stdout,*) 'IMAGINARY ENERGIES IN Ev:'
     do ii=1,qp%max_i
        write(stdout,'(''State:'',i5,'' GW (Im) :'',f10.5)') ii,aimag(qp%ene_gw(ii))*RYTOEV
     enddo
  endif
  !
  ! create the file bands.dat
  ! which can be used for gww_fit.x or projwfc.x
  !
  if(ionode) then
     !
     iun = find_free_unit()
     open(unit=iun, file='bands.dat', status='unknown',form='formatted')
     !
     write(iun,*) qp%max_i
     do ii=1,qp%max_i
        write(iun,'(i5,f10.5,f10.5,f10.5,f10.5)') ii, qp%ene_dft_ks(ii)*RYTOEV, &
                     & real(qp%ene_gw_pert(ii))*RYTOEV, &
                     & real(qp%ene_gw(ii))*RYTOEV, qp%ene_hf(ii)*RYTOEV
     enddo
     close(iun)
  endif
  !
  return

  END SUBROUTINE







 END MODULE energies_gww



