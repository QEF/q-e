MODULE input_simple_ip
!this module provides input file routines 

  USE kinds, ONLY: DP

   TYPE input_options_simple_ip
      CHARACTER(len=256) :: prefix = 'prefix'!prefix to designate the files same as in PW        
      CHARACTER(len=256) :: outdir = './'!outdir to designate the files same as in PW   
      INTEGER :: h_level = 2 ! terms included in the k-hamiltonian (for debug) 0=kinetic, 1=kinetic+local >1=kinetic+local+non-local
      INTEGER :: interp_grid(3)    !  k-grid used for the calculation of the optical spectra (Shirley interpolation)
      LOGICAL :: nonlocal_commutator=.true.  ! if true it includes the non-local part of [r,H] in the calculation of the optical spectra
      LOGICAL :: nonlocal_interpolation = .false.   ! if true it interpolates the non-local part of the psp (WARNING: experimental feature)
      LOGICAL :: tetrahedron_method = .false.    ! Use tetrahedron method to calculate Drude plasma frequency, Fermi energy and DOS
      REAL(kind=DP) :: fermi_degauss   ! degauss (in Ry) for Fermi level calculation
      REAL(kind=DP) :: fermi_energy = -1   ! input Fermi energy (in Ry) (if not given in input it is calculated using the Shirley interpolated bands)
      INTEGER  :: fermi_ngauss         ! ngauss for Fermi level calculation (n=-99 --> Fermi-Dirac ; n=-1 --> cold-smearing ; n>=0 --> Methfessel-Paxton)
      REAL(kind=DP) :: drude_degauss   ! degauss (in Ry) for Drude plasma frequency calculation
      INTEGER  :: drude_ngauss = -99        ! ngauss for Drude plasma frequency calculation (n=-99 --> Fermi-Dirac) (WARNING: only Fermi-Dirac smearing is implemented now)
      REAL(kind=DP) :: elec_temp = 0.0018375      ! electronic temperature (in Ry) for the occupation of the electronic states. Default is room temperature
      REAL(kind=DP) :: wmin, wmax         ! Minimun and maximum energy of the energy interval [wmin,wmax] (in Ry)
      INTEGER :: nw                 ! Number of points in the energy interval [wmin,wmax]
      REAL(kind=DP) ::  inter_broadening   ! broadening used for the optical spectra (interband part)
      REAL(kind=DP) ::  intra_broadening   ! broadening used for the optical spectra (intraband part)
      REAL(kind=DP) ::  delta_energy_dos = 0.000735  ! Energy spacing for DOS calculation (in Ry) (default is 0.01 eV)

   END TYPE input_options_simple_ip

   CONTAINS

     SUBROUTINE  read_input_simple_ip( simpleip_in )
       USE io_global,            ONLY : stdout, ionode, ionode_id
       USE mp,                   ONLY : mp_bcast
       USE mp_world,             ONLY : world_comm
       USE io_files,             ONLY : tmp_dir, prefix
       

       implicit none

       CHARACTER(LEN=256), EXTERNAL :: trimcheck
       TYPE(input_options_simple_ip) :: simpleip_in          !in output the input parameters

       NAMELIST/inputsimpleip/simpleip_in

	CHARACTER(LEN=256) :: outdir

       if(ionode) then
          read(*,NML=inputsimpleip)
          outdir = trimcheck(simpleip_in%outdir)
          tmp_dir = outdir
          prefix = trim(simpleip_in%prefix)
       endif

       call mp_bcast( outdir,ionode_id, world_comm )
       call mp_bcast( tmp_dir,ionode_id, world_comm )
       call mp_bcast( prefix,ionode_id, world_comm )
       call mp_bcast( simpleip_in%interp_grid, ionode_id, world_comm)
       call mp_bcast( simpleip_in%h_level, ionode_id, world_comm)
       call mp_bcast( simpleip_in%nonlocal_commutator, ionode_id, world_comm)
       call mp_bcast( simpleip_in%nonlocal_interpolation, ionode_id, world_comm)
       call mp_bcast( simpleip_in%fermi_degauss, ionode_id, world_comm)
       call mp_bcast( simpleip_in%fermi_energy, ionode_id, world_comm)
       call mp_bcast( simpleip_in%fermi_ngauss, ionode_id, world_comm)
       call mp_bcast( simpleip_in%drude_degauss, ionode_id, world_comm)
       call mp_bcast( simpleip_in%drude_ngauss, ionode_id, world_comm)
       call mp_bcast( simpleip_in%elec_temp, ionode_id, world_comm)
       call mp_bcast( simpleip_in%wmin, ionode_id, world_comm)
       call mp_bcast( simpleip_in%wmax, ionode_id, world_comm)
       call mp_bcast( simpleip_in%nw, ionode_id, world_comm)
       call mp_bcast( simpleip_in%inter_broadening, ionode_id, world_comm)
       call mp_bcast( simpleip_in%intra_broadening, ionode_id, world_comm)
       call mp_bcast( simpleip_in%tetrahedron_method, ionode_id, world_comm)
       call mp_bcast( simpleip_in%delta_energy_dos, ionode_id, world_comm)

       return

     END SUBROUTINE read_input_simple_ip

END MODULE input_simple_ip

