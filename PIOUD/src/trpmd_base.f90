!
! Copyright (C) 2025 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! Written by Aadhityan A, Lorenzo Paulatto, Michele Casula, Tommaso Morresi
!
!---------------------------------------------------------------------------
MODULE trpmd_base
  !---------------------------------------------------------------------------
  !
  ! ... This module contains most of the subroutines and functions needed by
  ! ... the implementation of "PIOUD" method into Quantum ESPRESSO
  !
  ! ... Other relevant files are:
  !
  ! ... ring_variables.f90
  ! ... trpmd_io_routines.f90
  ! ... ring_formats.f90
  ! ... compute_scf_pioud.f90
  !
  USE kinds,     ONLY : DP
  USE constants, ONLY : eps32, pi, autoev, bohr_radius_angs, eV_to_kelvin
  USE ring_io_units_module,  ONLY : iunpath
  USE io_global, ONLY : meta_ionode, meta_ionode_id
  USE mp,        ONLY : mp_bcast
  USE mp_world,  ONLY : world_comm
  !
  USE basic_algebra_routines
  !
  PRIVATE
  !
  PUBLIC :: initialize_polymer
  PUBLIC :: explore_phasespace
  !
  CONTAINS
    !
    ! ... module procedures
    !
    !-----------------------------------------------------------------------
    SUBROUTINE initialize_polymer()
      !-----------------------------------------------------------------------
      !
      USE control_flags,    ONLY : conv_elec
      USE ions_base,        ONLY : amass, ityp
      USE io_files,         ONLY : prefix, tmp_dir
      USE mp_images,        ONLY : nimage
      USE ring_input_parameters_module, ONLY : pos_      => pos, &
                                               input_images, nstep_path_ => nstep_path
      USE ring_input_parameters_module, ONLY : restart_mode
      USE ring_input_parameters_module, ONLY : nat
      USE ring_variables,   ONLY : pos, istep_path, nstep_path,    &
                                   dim1, & 
                                    path_length,  &
                                   tune_load_balance,  &  
                                   pending_image
      USE ring_variables,   ONLY : path_allocation
      USE fcp_variables,        ONLY : lfcpopt
      USE fcp_opt_routines,     ONLY : fcp_opt_allocation
      use pimd_variables, ONLY : nbeadMD, pes
      use mp_world
      !
      IMPLICIT NONE
      !
      INTEGER :: i, fii, lii
      !
      ! ... istep_path is initialised to 1
      !
      istep_path    = 0    !!! <----my mod.
      pending_image = 0
      conv_elec     = .TRUE.
      !
      ! ... the dimension of all "path" arrays (dim1) is set here
      ! ... ( it corresponds to the dimension of the configurational space )
      !
      !
      dim1 = 3*nat
      !
      !
      ! IF ( nimage > 1 ) THEN
      !    !
      !    ! ... the automatic tuning of the load balance in
      !    ! ... image-parallelisation is switched off by default
      !    !
      !    tune_load_balance = .FALSE.
      !    !
      !    ! ... in the case of image-parallelisation the number of images
      !    ! ... to be optimised must be larger than nimage
      !    !
      !    IF ( first_last_opt ) THEN
      !       !
      !       fii = 1
      !       lii = num_of_images
      !       !
      !    ELSE
      !       !
      !       fii = 2
      !       lii = num_of_images - 1
      !       !
      !    END IF
      !    !
      !    IF ( nimage > ( lii - fii + 1 ) ) &
      !       CALL errore( 'initialize_polymer', 'nimage is ' // &
      !                  & 'larger than the available number of images', 1 )
      !    !
      ! END IF
      ! num_of_images=nbeadMD
      ! write(12000+mpime,*) nbeadMD,nimage
      if( nimage > nbeadMD ) &
             CALL errore( 'initialize_polymer', 'nimage is ' // &
                        & 'more computing images than polymer images', 1 )
      fii = 1
      lii = nbeadMD      
      !
      ! ... dynamical allocation of arrays
      !
      CALL path_allocation()
      if ( lfcpopt ) CALL fcp_opt_allocation()
      !
      ! use_masses = .TRUE.

      ! IF ( use_masses ) THEN
      !    !
      !    ! ... mass weighted coordinates are used
      !    !
      !    DO i = 1, nat
      !       !
      !       mass(3*i-2) = amass(ityp(i))
      !       mass(3*i-1) = amass(ityp(i))
      !       mass(3*i-0) = amass(ityp(i))
      !       !
      !    END DO
      !    !
      ! ELSE
      !    !
      !    mass = 1.0_DP
         !
      ! END IF
      !
      ! ... initialization of the allocatable arrays
      !
      !pos(:,1:input_images) = pos_(1:dim1,1:input_images)
      pos(:,1) = pos_(:,1)
      !
      pes          = 0.0_DP
      ! grad_pes     = 0.0_DP
      ! grad         = 0.0_DP
      ! error        = 0.0_DP
      ! frozen       = .FALSE.
      !
      ! k = k_min
      !
      CALL initial_guess()
      !
      ! posold(:,:) = pos(:,:)
      !
      ! ... the actual number of degrees of freedom is computed
      !
      ! deg_of_freedom = 0
      ! !
      ! DO i = 1, nat
      !    !
      !    deg_of_freedom = deg_of_freedom + 3
      !    !
      ! END DO
      !
      RETURN
      !
    END SUBROUTINE initialize_polymer
    !
    !--------------------------------------------------------------------
    SUBROUTINE initial_guess()
      !--------------------------------------------------------------------

      USE ring_input_parameters_module, ONLY : input_images
      USE ring_variables,   ONLY : pos, dim1, path_length ! num_of_images, 
      USE ring_io_units_module,         ONLY : iunpath
      USE pimd_variables, ONLY : nbeadMD
      !
      IMPLICIT NONE
      !
      REAL(DP) :: s
      INTEGER  :: i, j
      LOGICAL  :: tooclose
      REAL(DP), ALLOCATABLE :: pos_n(:,:), dr(:,:), image_spacing(:)
      !
      !
      CALL mp_bcast( pos(:,1),         meta_ionode_id, world_comm )
      if( nbeadMD > 1) then  !!! <----my mod.
          do i=2,nbeadMD    !!! <----my mod.
        
            pos(:,i) = pos(:,1)   !!! <----my mod.
        
          end do    !!! <----my mod.
        end if   !!! <----my mod.
      path_length = 0._dp !
      RETURN
      !
    END SUBROUTINE initial_guess
    !
    ! ... shared routines
    !
    !------------------------------------------------------------------------
    SUBROUTINE born_oppenheimer_pes( stat )
      !------------------------------------------------------------------------
      !
      USE ring_variables, ONLY : pending_image, istep_path !pes, num_of_images
      USE pimd_variables, ONLY : nbeadMD, pes
      !
      IMPLICIT NONE
      !
      LOGICAL, INTENT(OUT) :: stat
      !
      INTEGER  :: fii, lii
      !
      fii = 1
      lii = nbeadMD
      !
      IF ( pending_image /= 0 ) fii = pending_image
      !
      CALL compute_scf_pioud( fii, lii, stat )
      !
      IF ( .NOT. stat ) RETURN
      !
      ! Emin       = MINVAL( pes(1:num_of_images) )
      ! Emax       = MAXVAL( pes(1:num_of_images) )
      ! Emax_index = MAXLOC( pes(1:num_of_images), 1 )
      !
      RETURN
      !
    END SUBROUTINE born_oppenheimer_pes
    !
    !-----------------------------------------------------------------------
    SUBROUTINE explore_phasespace()
      !-----------------------------------------------------------------------
      !
      USE ring_variables,    ONLY :  pos, nstep_path ! ,lneb, lsmd
      USE ring_variables,   ONLY :  istep_path,   &
                                   pending_image !, conv_path &
                                   !pes !CI_scheme
                                  !  Emax_index
      USE trpmd_io_routines, ONLY : write_output
      USE ring_formats,     ONLY : scf_iter_fmt
      USE fcp_variables,    ONLY : lfcpopt
      use mp_world
      !
      USE pimd_variables,   ONLY : forceMD,rpos,unit_dot_out,irun,restart_pimd,nblocks,nstep_block    !!! <----my mod.
      !
      IMPLICIT NONE
      !
      LOGICAL :: stat
      REAL(DP) :: fcp_err_max = 0.0_DP
      REAL(DP) :: potenergy  !!! <----my mod.
      !
      REAL(DP), EXTERNAL :: get_clock
      !
      ! conv_path = .FALSE.
      !
      CALL explore_phasespace_init()
      !
      
      
      ! write(13000+mpime,*)istep_path,nstep_path,pending_image
      IF ( istep_path == nstep_path ) THEN
         !
         CALL write_output()
         !
         pending_image = 0
         !
         RETURN
         !
      END IF
      ! write(13000+mpime,*)istep_path,nstep_path,pending_image
      !IF ( meta_ionode ) 
      CALL pimd_allocation !!! <----my mod.
      
      IF ( meta_ionode .and. restart_pimd) CALL pimd_restart_traj  !!! <----my mod.
      CALL mp_bcast( pos,  meta_ionode_id, world_comm )   !!! <----my mod.
      
      ! write(13000+mpime,*)istep_path,nstep_path,pending_image

      if ( meta_ionode) then 
        write(*,*)
        write(*,*) '-------- computes initial forces and energy --------------' 
      end if
      CALL born_oppenheimer_pes( stat )  !!! <----my mod.
      if ( meta_ionode) then
        write(*,*) '----------------------------------------------------------'
        write(*,*) 
      end if   

      IF (meta_ionode) CALL pimd_get_pot_from_pw(potenergy) !!! <----my mod.
      ! IF (meta_ionode) CALL pimd_get_force_from_pw(forceMD) !!! <----my mod.
      IF (meta_ionode) CALL pimd_pw_convert_pos('pw_to_md')    !!! <----my mod.
      IF (meta_ionode) CALL pimdnvt_init(potenergy)!,forceMD) !!! <----my mod.
      !
      ! ... path optimisation loop
      !
      ! write(13000+mpime,*)istep_path,nstep_path,pending_image

      optimisation: DO
         !
         ! ... new positions are saved on file:  it has to be done here
         ! ... because, in the event of an unexpected crash the new positions
         ! ... would be lost. At this stage the forces and the energies are
         ! ... not yet known (but are not necessary for restarting); the
         ! ... restart file is written again as soon as the energies and
         ! ... forces have been computed.
         !
         !

         IF ( meta_ionode ) &
            WRITE( UNIT = iunpath, FMT = scf_iter_fmt ) istep_path + 1
         
         !
         ! ... energies and gradients acting on each image of the path (in real
         ! ... space) are computed calling a driver for the scf calculations
         !
         !
         IF (meta_ionode .and. irun.eq.0) THEN   !!! <----my mod. 
           
             call pimd_pw_convert_pos('pw_to_md')    !!! <----my mod. 
             CALL vel_verlet_1half_irun0(istep_path)   !!! <----my mod.
             call pimd_pw_convert_pos('md_to_pw')     !!! <----my mod.
             
         ELSE IF ( meta_ionode .and. irun.eq.3)  THEN   !!! <----my mod. 

             call pimd_pw_convert_pos('pw_to_md')    !!! <----my mod.        
             CALL prop_ceriotti_1half_irun3(istep_path)    !!! <----my mod.
             call pimd_pw_convert_pos('md_to_pw')     !!! <----my mod.    
             
         END IF   !!! <----my mod. 
         
         CALL mp_bcast( pos,  meta_ionode_id, world_comm )   !!! <----my mod.
         
         CALL born_oppenheimer_pes( stat )
         
         IF (meta_ionode) CALL pimd_get_pot_from_pw(potenergy) !!! <----my mod.
        !  IF (meta_ionode) CALL pimd_get_force_from_pw(forceMD) !!! <----my mod.        
         
         !
        !  IF ( .NOT. stat ) THEN
        !     !
        !     conv_path = .FALSE.
        !     !
        !     EXIT optimisation
        !     !
        !  END IF
         !
         ! ... istep_path is updated after a self-consistency step has been
         ! ... completed
         !
         istep_path = istep_path + 1
         !
         ! ... information is written on the standard output
         !
         CALL write_output()   
         !
         ! ... the restart file is written
         !
         IF (meta_ionode)    call pimd_pw_convert_pos('pw_to_md')    !!! <----my mod. 
         IF (meta_ionode)    CALL pimdnvt(potenergy)!,forceMD,)     !!! <----my mod.
         IF (meta_ionode)    call pimd_pw_convert_pos('md_to_pw')     !!! <----my mod.
         CALL mp_bcast( pos,  meta_ionode_id, world_comm )
     
         IF ( istep_path .eq. nstep_path) THEN   !!! <----my mod.
         
           EXIT optimisation                     !!! <----my mod.
         
         END IF                                  !!! <----my mod.
         !
      END DO optimisation
      !
      RETURN
      !
    END SUBROUTINE explore_phasespace
    !
    !------------------------------------------------------------------------
    SUBROUTINE explore_phasespace_init()
      !------------------------------------------------------------------------
      !
      USE ring_variables, ONLY : pending_image
      !
      IMPLICIT NONE
      !
      IF ( pending_image /= 0 ) RETURN
      !
      RETURN
      !
    END SUBROUTINE explore_phasespace_init
    !
END MODULE trpmd_base
