MODULE simple_ip_objects
  ! this module describes the most important objects
  !
  USE kinds, ONLY : DP
  !
  TYPE energies
     INTEGER :: nk!total number of k points
     INTEGER :: nk_loc!local number of k points
     INTEGER :: ik_first!first local k point
     INTEGER :: ik_last!last local k point  
     INTEGER :: num_bands   ! number of states
     REAL(kind=DP), DIMENSION(:,:), POINTER :: energy ! energies (num_bands,nk_loc)
     REAL(kind=DP), DIMENSION(:,:,:), POINTER :: energy_der ! derivatives of the energy  (3,num_bands,nk_loc)
  END TYPE
  !
  TYPE kpoints
     INTEGER :: nkgrid(3)         !  k-grid for interpolation (interp_grid)
     REAL(kind=DP), DIMENSION(:,:), POINTER :: qk ! kpoints coordinates (3,nkgrid(1)*nkgrid(2)*nkgrid(3))
     REAL(kind=DP) :: alat        ! lattice paramater 
     REAL(kind=DP) :: bg(3,3)     ! reciprocal basis vectors
     INTEGER :: nk_loc            ! local number of k points
     INTEGER :: ik_first          ! first local k point
     INTEGER :: ik_last           ! last local k point  
     INTEGER :: nk_smooth_loc     ! number of local k-points
     INTEGER :: nk                ! total number of k-points
     REAL(kind=DP), DIMENSION(:,:), POINTER  ::  pos_cube   !   (3,nk_loc) 
     INTEGER, DIMENSION(:,:), POINTER :: coord_cube   ! (8,nk_loc)
  END TYPE
  !
  TYPE shirley
      INTEGER :: ntot_e ! number of optimal basis states
      LOGICAL :: noncolin ! if it is non-collinear
      INTEGER :: nat ! number of atoms
      INTEGER :: ntyp ! number of types of atoms
      INTEGER :: nhm  ! max number of different beta functions per atom
      INTEGER :: nspin  ! number of spins (1=no spin, 2=LSDA)
      INTEGER :: nkb  ! total number of beta functions
      INTEGER :: npol  !
      INTEGER :: nks ! number of k-points of the smooth grid
      INTEGER :: nk_smooth_loc ! number of local k-points
      INTEGER :: num_val , num_cond  ! number of interpolated valence bands, number of interpolated conduction bands
      INTEGER :: num_nbndv(2)       ! total number of occupied bands
      INTEGER :: num_bands       ! total number of considered bands (num_bands = num_val + num_cond)
      LOGICAL :: nonlocal_commutator  !
      INTEGER, DIMENSION(3) :: nkpoints     ! smooth k-points grid on which H(k) is calculated
      INTEGER, DIMENSION(:), POINTER :: ityp ! (nat)
      INTEGER, DIMENSION(:), POINTER :: nh   ! (ntyp)
      INTEGER, DIMENSION(:), POINTER :: indv_ijkb0 ! (nat)
      REAL(kind=8) :: alat, nelec, omega        ! lattice paramater, number of electrons, volume prim. cell
      REAL(kind=8) :: bg(3,3)     ! reciprocal basis vectors
      REAL(kind=8) :: at(3,3)     ! direct basis vectors
      REAL(kind=DP) :: s_bands     ! threshold for the construction of the optimal basis
      COMPLEX(kind=DP), DIMENSION(:,:), POINTER ::  h0 , Vloc ! k-dependent Hamiltonian  (ntot_e,ntot_e)
      COMPLEX(kind=DP), DIMENSION(:,:,:), POINTER :: h1 ! (ntot_e,ntot_e,3)
      COMPLEX(kind=DP), DIMENSION(:,:,:,:), POINTER :: deeq_nc  ! (nhm,nhm,nat,nspin)
      COMPLEX(kind=DP), DIMENSION(:,:,:), POINTER :: deeqc !(nhm,nhm,nat)
      COMPLEX(kind=DP), DIMENSION(:,:,:,:), POINTER :: beck_nc ! (nkb,npol,ntot_e,nk)
      COMPLEX(kind=DP), DIMENSION(:,:,:), POINTER :: beckc ! (nkb,ntot_e,nk)
      COMPLEX(kind=DP),DIMENSION(:,:,:,:), POINTER :: commut ! (ntot_e,ntot_e,3,nk)
  END TYPE
  !
  TYPE eigen
    INTEGER :: num_bands       ! number of states
    REAL(kind=DP) :: q(3)        ! k-point
    REAL(kind=DP), DIMENSION(:), POINTER  :: energy     ! eigenenergies (num_bands)
    COMPLEX(kind=DP), DIMENSION(:,:), POINTER  :: wave_func      ! eigenfunctions (ntot_e,num_bands)
  END TYPE

   CONTAINS

   subroutine initialize_energies(element)
      implicit none
      TYPE(energies) :: element
      nullify(element%energy)
      nullify(element%energy_der)
      return
    end subroutine initialize_energies
    
    subroutine deallocate_energies(element)
      implicit none
      TYPE(energies) :: element
      if(associated(element%energy)) deallocate(element%energy)
      nullify(element%energy)
      if(associated(element%energy_der)) deallocate(element%energy_der)
      nullify(element%energy_der)
      return
    end subroutine deallocate_energies

 
   subroutine initialize_kpoints(element)
      implicit none
      TYPE(kpoints) :: element
      nullify(element%qk)
      nullify(element%pos_cube)
      nullify(element%coord_cube)
      return
    end subroutine initialize_kpoints
    
    SUBROUTINE deallocate_kpoints(element)
      implicit none
      TYPE(kpoints) :: element
      if(associated(element%qk)) deallocate(element%qk)
      nullify(element%qk)
      if(associated(element%pos_cube)) deallocate(element%pos_cube)
      nullify(element%pos_cube)
      if(associated(element%coord_cube)) deallocate(element%coord_cube)
      nullify(element%coord_cube)
      return
    END SUBROUTINE deallocate_kpoints

    SUBROUTINE initialize_shirley(element)
      implicit none
      TYPE(shirley) :: element
      nullify(element%ityp)
      nullify(element%nh)
      nullify(element%indv_ijkb0)
      nullify(element%h0)
      nullify(element%h1)
      nullify(element%Vloc)
      nullify(element%deeq_nc)
      nullify(element%deeqc)
      nullify(element%beck_nc)
      nullify(element%beckc)
      nullify(element%commut)
      return
    END SUBROUTINE initialize_shirley
    
    SUBROUTINE deallocate_shirley(element)
      implicit none
      TYPE(shirley) :: element
      if(associated(element%ityp)) deallocate(element%ityp)
      nullify(element%ityp)
      if(associated(element%nh)) deallocate(element%nh)
      nullify(element%nh)
      if(associated(element%indv_ijkb0)) deallocate(element%indv_ijkb0)
      nullify(element%indv_ijkb0)
      if(associated(element%h0)) deallocate(element%h0)
      nullify(element%h0)
      if(associated(element%h1)) deallocate(element%h1)
      nullify(element%h1)
      if(associated(element%Vloc)) deallocate(element%Vloc)
      nullify(element%Vloc)
      if(associated(element%deeq_nc)) deallocate(element%deeq_nc)
      nullify(element%deeq_nc)
      if(associated(element%deeqc)) deallocate(element%deeqc)
      nullify(element%deeqc)
      if(associated(element%beck_nc)) deallocate(element%beck_nc)
      nullify(element%beck_nc)
      if(associated(element%beckc)) deallocate(element%beckc)
      nullify(element%beckc)
      if(associated(element%commut)) deallocate(element%commut)
      nullify(element%commut)
      return
    END SUBROUTINE deallocate_shirley

   SUBROUTINE initialize_eigen(element)
      implicit none
      TYPE(eigen) :: element
      nullify(element%energy)
      nullify(element%wave_func)
      return
    END SUBROUTINE initialize_eigen
    
    SUBROUTINE deallocate_eigen(element)
      implicit none
      TYPE(eigen) :: element
      if(associated(element%energy)) deallocate(element%energy)
      nullify(element%energy)
      if(associated(element%wave_func)) deallocate(element%wave_func)
      nullify(element%wave_func)
      return
    END SUBROUTINE deallocate_eigen



    SUBROUTINE read_shirley(simpleip_in,sh)
      USE input_simple_ip, ONLY : input_options_simple_ip
      USE mp,                   ONLY : mp_bcast
      USE mp_world,             ONLY : world_comm
      USE io_files,  ONLY : tmp_dir
      USE io_global, ONLY : ionode_id, ionode, stdout

      implicit none

      TYPE(input_options_simple_ip) :: simpleip_in
      TYPE(shirley) :: sh
      INTEGER, EXTERNAL :: find_free_unit
      INTEGER :: iun, idir, nk
      write(stdout,*)'simple_ip: opening file "hamiltonian"'
      if(ionode) then
         iun = find_free_unit()
         open( unit=iun, file=trim(tmp_dir)//trim(simpleip_in%prefix)//'.hamiltonian', status='old',form='unformatted')
         write(stdout,*)'File opened'
         read(iun) sh%ntot_e
      endif   
      call mp_bcast(sh%ntot_e,ionode_id,world_comm)
      allocate(sh%h0(sh%ntot_e,sh%ntot_e), sh%h1(sh%ntot_e,sh%ntot_e,3), sh%Vloc(sh%ntot_e,sh%ntot_e))

      if (ionode) read(iun) sh%h0(1:sh%ntot_e,1:sh%ntot_e) 

      do idir = 1,3
        if (ionode) read(iun) sh%h1(1:sh%ntot_e,1:sh%ntot_e,idir) 
      enddo

      if (ionode) read(iun) sh%Vloc(1:sh%ntot_e,1:sh%ntot_e) 

      call mp_bcast(sh%h0,ionode_id,world_comm)
      call mp_bcast(sh%h1,ionode_id,world_comm)
      call mp_bcast(sh%Vloc,ionode_id,world_comm)
      
      if(ionode) then
        read(iun) sh%noncolin
        read(iun) sh%nat 
        read(iun) sh%ntyp 
        read(iun) sh%nhm
        read(iun) sh%nspin
        read(iun) sh%nkb 
        read(iun) sh%npol
      endif
      call mp_bcast(sh%noncolin,ionode_id,world_comm)
      call mp_bcast(sh%nat,ionode_id,world_comm)
      call mp_bcast(sh%ntyp,ionode_id,world_comm)
      call mp_bcast(sh%nhm,ionode_id,world_comm)
      call mp_bcast(sh%nspin,ionode_id,world_comm)
      call mp_bcast(sh%nkb,ionode_id,world_comm)
      call mp_bcast(sh%npol,ionode_id,world_comm)
      

      allocate(sh%ityp(sh%nat), sh%nh(sh%ntyp), sh%indv_ijkb0(sh%nat))
      if(ionode) then
        read(iun) sh%ityp(1:sh%nat)
        read(iun) sh%nh(1:sh%ntyp)
        read(iun) sh%indv_ijkb0(1:sh%nat)
        read(iun) sh%nkpoints
      endif
      call mp_bcast(sh%ityp,ionode_id,world_comm)
      call mp_bcast(sh%nh,ionode_id,world_comm)
      call mp_bcast(sh%indv_ijkb0,ionode_id,world_comm)
      call mp_bcast(sh%nkpoints,ionode_id,world_comm)
      
      nk = (sh%nkpoints(1))*(sh%nkpoints(2))*(sh%nkpoints(3))

      allocate(sh%deeqc(sh%nhm,sh%nhm,sh%nat), sh%deeq_nc(sh%nhm,sh%nhm,sh%nat,sh%nspin))
      if (sh%noncolin) then
         if (ionode) read(iun) sh%deeq_nc(1:sh%nhm,1:sh%nhm,1:sh%nat,1:sh%nspin)
      else
         if (ionode) read(iun) sh%deeqc(1:sh%nhm,1:sh%nhm,1:sh%nat)
      endif
      
      if (sh%noncolin) then
         call mp_bcast(sh%deeq_nc,ionode_id,world_comm)
      else
         call mp_bcast(sh%deeqc,ionode_id,world_comm)
      endif
            
      if(ionode) then
         read(iun) sh%alat
         read(iun) sh%bg(1:3,1:3)
         read(iun) sh%at(1:3,1:3)
         read(iun) sh%nelec
         read(iun) sh%omega
         read(iun) sh%num_val
         read(iun) sh%num_cond
         read(iun) sh%num_nbndv
         read(iun) sh%nonlocal_commutator
         read(iun) sh%s_bands
      endif
      sh%num_bands = sh%num_val + sh%num_cond
      call mp_bcast(sh%alat,ionode_id,world_comm)
      call mp_bcast(sh%bg,ionode_id,world_comm)
      call mp_bcast(sh%at,ionode_id,world_comm)
      call mp_bcast(sh%nelec,ionode_id,world_comm)
      call mp_bcast(sh%omega,ionode_id,world_comm)
      call mp_bcast(sh%num_val,ionode_id,world_comm)
      call mp_bcast(sh%num_cond,ionode_id,world_comm)
      call mp_bcast(sh%num_nbndv,ionode_id,world_comm)
      call mp_bcast(sh%num_bands,ionode_id,world_comm) 
      call mp_bcast(sh%nonlocal_commutator,ionode_id,world_comm)
      call mp_bcast(sh%s_bands,ionode_id,world_comm)
      !
      if(ionode) then
         close(iun)
         write(stdout,*)'File closed'
      endif
      !
      if (.not. simpleip_in%nonlocal_interpolation .and. ( sh%nkpoints(1) /= simpleip_in%interp_grid(1) .or. &
      & sh%nkpoints(2) /= simpleip_in%interp_grid(2) .or. sh%nkpoints(3) /= simpleip_in%interp_grid(3) ) ) then
        call errore('SIMPLE_IP', 'W/o trilinear interpolation, k-grids from simple and simple_ip have to be equal',1)
      elseif ( simpleip_in%nonlocal_interpolation .and. ( 2*sh%nkpoints(1) /= simpleip_in%interp_grid(1) .or. &
      & 2*sh%nkpoints(2) /= simpleip_in%interp_grid(2) .or. 2*sh%nkpoints(3) /= simpleip_in%interp_grid(3) ) ) then
        call errore('SIMPLE_IP', 'W/ trilinear interpolation, simple_ip k-grid has to be the double of the simple k-grid',1)
      endif
      if( (sh%num_val .ne. sh%num_nbndv(1)) ) then
        call errore('SIMPLE_IP', 'All the occupied bands must be included for Shirley interpolation (num_val=num_nbndv)',1)
      endif
      !
    end subroutine read_shirley


    SUBROUTINE kgrid_creation(simpleip_in,kgrid, sh)
      USE input_simple_ip, ONLY : input_options_simple_ip
      USE mp_world,  ONLY : mpime, nproc
      !
      implicit none
      !
      TYPE(input_options_simple_ip) :: simpleip_in
      TYPE(kpoints) :: kgrid
      TYPE(shirley) :: sh
      INTEGER ::  i, j, k, ii
      INTEGER :: l_blk
      !
      ! WARNING: this must be the same as in create_energies
      kgrid%nkgrid(1:3) = simpleip_in%interp_grid(1:3)
      kgrid%nk = (kgrid%nkgrid(1))*(kgrid%nkgrid(2))*(kgrid%nkgrid(3))
      kgrid%alat = sh%alat
      kgrid%bg(1:3,1:3) = sh%bg(1:3,1:3)
      
      l_blk=kgrid%nk/nproc
      if(l_blk*nproc<kgrid%nk) l_blk=l_blk+1

      if(l_blk*mpime+1 <= kgrid%nk) then
         kgrid%ik_first=l_blk*mpime+1 
         kgrid%ik_last=kgrid%ik_first+l_blk-1
         if(kgrid%ik_last>kgrid%nk) kgrid%ik_last=kgrid%nk
         kgrid%nk_loc=kgrid%ik_last-kgrid%ik_first+1
      else
         kgrid%nk_loc=0
         kgrid%ik_first=0
         kgrid%ik_last=-1
      endif     
      
      allocate(kgrid%qk(3,kgrid%nk))
      ! Uniform grid in [0,1)*bg --> it goes outside first BZ
      ii = 0
      do i=0,kgrid%nkgrid(1)-1
       do j=0,kgrid%nkgrid(2)-1
        do k=0,kgrid%nkgrid(3)-1
           ii = ii + 1
           kgrid%qk(1:3,ii) = kgrid%bg(1:3,1)*dble(i)/dble(kgrid%nkgrid(1)) + &
        & kgrid%bg(1:3,2)*dble(j)/dble(kgrid%nkgrid(2)) + kgrid%bg(1:3,3)*dble(k)/dble(kgrid%nkgrid(3))
        enddo
       enddo
      enddo
      !
    END SUBROUTINE kgrid_creation



    SUBROUTINE create_energies(sh,kgrid,ene)
      USE mp_world,  ONLY : mpime, nproc
      implicit none
      
      TYPE(shirley) :: sh
      TYPE(kpoints) :: kgrid
      TYPE(energies) :: ene
      INTEGER :: l_blk
      
      ene%nk = (kgrid%nkgrid(1))*(kgrid%nkgrid(2))*(kgrid%nkgrid(3))
      ene%num_bands = sh%num_bands

      l_blk=ene%nk/nproc
      if(l_blk*nproc<ene%nk) l_blk=l_blk+1

      if(l_blk*mpime+1 <= ene%nk) then
         ene%ik_first=l_blk*mpime+1 
         ene%ik_last=ene%ik_first+l_blk-1
         if(ene%ik_last>ene%nk) ene%ik_last=ene%nk
         ene%nk_loc=ene%ik_last-ene%ik_first+1
      else
         ene%nk_loc=0
         ene%ik_first=0
         ene%ik_last=-1
      endif
      
      allocate(ene%energy(ene%num_bands,ene%nk_loc), ene%energy_der(3,ene%num_bands,ene%nk_loc))

    END SUBROUTINE create_energies




    SUBROUTINE read_shirley_k(simpleip_in,sh,ene)
      USE input_simple_ip, ONLY : input_options_simple_ip
      USE mp,                   ONLY : mp_bcast, mp_barrier
      USE mp_world,             ONLY : world_comm
      USE io_files,  ONLY : tmp_dir
      USE io_global, ONLY : ionode_id, ionode, stdout
      !
      IMPLICIT NONE
      !
      TYPE(input_options_simple_ip) :: simpleip_in
      TYPE(shirley) :: sh
      TYPE(energies) :: ene
      !
      INTEGER, EXTERNAL :: find_free_unit
      INTEGER :: iun, ll, ii , jj , kk
      COMPLEX(kind=DP), DIMENSION(:,:,:), ALLOCATABLE :: sum_beckc_tmp 
      COMPLEX(kind=DP), DIMENSION(:,:,:,:), ALLOCATABLE :: sum_beck_nc_tmp , sum_commut_tmp 
      !
      write(stdout,*)'simple_ip: opening file "hamiltonian_k"'
      if(ionode) then
         iun = find_free_unit()
         open( unit=iun, file=trim(tmp_dir)//trim(simpleip_in%prefix)//'.hamiltonian_k', status='old',form='unformatted')
         write(stdout,*)'File opened'
      endif   
      !
      if (sh%noncolin) then
        allocate(sum_beck_nc_tmp(sh%nkb,sh%npol,sh%ntot_e,sh%nkpoints(3)))
        allocate(sh%beck_nc(sh%nkb,sh%npol,sh%ntot_e,ene%nk_loc))
      else
        allocate(sum_beckc_tmp(sh%nkb,sh%ntot_e,sh%nkpoints(3)))
        allocate(sh%beckc(sh%nkb,sh%ntot_e,ene%nk_loc))
      endif
      !
      if (sh%nonlocal_commutator) then
        allocate(sum_commut_tmp(sh%ntot_e,sh%ntot_e,3,sh%nkpoints(3)))
        allocate(sh%commut(sh%ntot_e,sh%ntot_e,3,ene%nk_loc))
      endif
      ! 
      write(stdout,*) 'Total number of k-point blocks (for reading): ' , sh%nkpoints(1)*sh%nkpoints(2)
      do ii=1,sh%nkpoints(1)
        do jj=1,sh%nkpoints(2)
        ! 
        write(stdout,*) 'k-point block: ' , jj  + (ii - 1)*sh%nkpoints(2)
        if (sh%noncolin) then
          if (ionode) read(iun) sum_beck_nc_tmp(1:sh%nkb,1:sh%npol,1:sh%ntot_e,1:sh%nkpoints(3))
        else
          if (ionode) read(iun) sum_beckc_tmp(1:sh%nkb,1:sh%ntot_e,1:sh%nkpoints(3))
        endif
        !
        if (sh%nonlocal_commutator) then
          if (ionode) read(iun) sum_commut_tmp(1:sh%ntot_e,1:sh%ntot_e,1:3,1:sh%nkpoints(3))  ! read commutator matrix
        endif
        !
        if (sh%noncolin) then
          call mp_bcast(sum_beck_nc_tmp,ionode_id,world_comm)
        else
          call mp_bcast(sum_beckc_tmp,ionode_id,world_comm)
        endif
        !
        if (sh%nonlocal_commutator) then
          call mp_bcast(sum_commut_tmp,ionode_id,world_comm)
        endif
        !
        do kk=1,sh%nkpoints(3)
          !
          ll = kk + (jj - 1)*sh%nkpoints(3) + (ii - 1)*sh%nkpoints(2)*sh%nkpoints(3)   ! global kindex
          !
          if (ene%ik_first <= ll .and. ll <= ene%ik_last) then
            !
            if (sh%noncolin) then
              sh%beck_nc(1:sh%nkb, 1:sh%npol, 1:sh%ntot_e,ll - ene%ik_first + 1) = &
              & sum_beck_nc_tmp(1:sh%nkb, 1:sh%npol, 1:sh%ntot_e, kk)
            else
              sh%beckc(1:sh%nkb, 1:sh%ntot_e, ll - ene%ik_first + 1) = &
              & sum_beckc_tmp(1:sh%nkb,1:sh%ntot_e, kk)
            endif
            !
            if (sh%nonlocal_commutator) then
              sh%commut(1:sh%ntot_e, 1:sh%ntot_e, 1:3, ll - ene%ik_first + 1) = &
              & sum_commut_tmp(1:sh%ntot_e, 1:sh%ntot_e, 1:3, kk)
            endif
          endif
          ! 
        enddo
        !
        enddo
      enddo
      !
      if(ionode) then 
        close(iun)
        write(stdout,*)'File closed'
      endif
      !
      if (sh%noncolin) then
        deallocate(sum_beck_nc_tmp)
      else
        deallocate(sum_beckc_tmp)
      endif
      !
      if (sh%nonlocal_commutator) then
        deallocate(sum_commut_tmp)
      endif
      !
    END SUBROUTINE read_shirley_k




    !
    SUBROUTINE read_shirley_k_interp(simpleip_in,sh,ene,k)
      USE input_simple_ip, ONLY : input_options_simple_ip
      USE mp,                   ONLY : mp_bcast, mp_barrier
      USE mp_world,             ONLY : world_comm
      USE io_files,  ONLY : tmp_dir
      USE io_global, ONLY : ionode_id, ionode, stdout
   
      implicit none

      TYPE(input_options_simple_ip) :: simpleip_in
      TYPE(shirley) :: sh
      TYPE(energies) :: ene
      TYPE(kpoints)  :: k
      REAL(kind=DP) :: qproj(3), q(3)
      INTEGER :: i,j,ik, n(3), np(3), ii, jj , kk , ll , counter, nk_smooth, iun
      INTEGER, EXTERNAL :: find_free_unit
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: coord_cube_global
      INTEGER, DIMENSION(:), ALLOCATABLE  :: diff_kpoints    
      COMPLEX(kind=DP), DIMENSION(:,:,:), ALLOCATABLE :: sum_beckc_tmp 
      COMPLEX(kind=DP), DIMENSION(:,:,:,:), ALLOCATABLE :: sum_beck_nc_tmp , sum_commut_tmp   

      nk_smooth = sh%nkpoints(1)*sh%nkpoints(2)*sh%nkpoints(3)  ! initial k-grid given to simple.x

      allocate(k%pos_cube(3,k%nk_loc),k%coord_cube(8,k%nk_loc))
      allocate(coord_cube_global(8,k%nk_loc),diff_kpoints(nk_smooth))

      diff_kpoints = 0
      coord_cube_global = 0
      do ik = 1,k%nk_loc
         q(1:3) = k%qk(1:3, ik + k%ik_first - 1)

        ! Project q on the bg basis (it is obtained doing the scalar product of q with the direct lattice vectors at)
        qproj(1:3) = 0
        do i=1,3
          do j=1,3
            qproj(i) = qproj(i) + q(j)*sh%at(j,i) 
          enddo
        enddo
  
        n(1:3) = 0
        do i=1,3
          n(i) = int(qproj(i)*sh%nkpoints(i))
          np(i) = n(i) + 1
          !if ( np(i) >= sh%nkpoints(i) )  np(i) = 0
          ! WARNING: When I reach the border of the grid I take the projector at the previous k-point.
          if ( np(i) >= sh%nkpoints(i) ) np(i) = sh%nkpoints(i) - 1  
        enddo   
  
        coord_cube_global(1,ik) = n(1)*sh%nkpoints(2)*sh%nkpoints(3) + n(2)*sh%nkpoints(3) + n(3) + 1
        coord_cube_global(2,ik) = np(1)*sh%nkpoints(2)*sh%nkpoints(3) + n(2)*sh%nkpoints(3) + n(3) + 1
        coord_cube_global(3,ik) = n(1)*sh%nkpoints(2)*sh%nkpoints(3) + np(2)*sh%nkpoints(3) + n(3) + 1
        coord_cube_global(4,ik) = np(1)*sh%nkpoints(2)*sh%nkpoints(3) + np(2)*sh%nkpoints(3) + n(3) + 1
        coord_cube_global(5,ik) = n(1)*sh%nkpoints(2)*sh%nkpoints(3) + n(2)*sh%nkpoints(3) + np(3) + 1
        coord_cube_global(6,ik) = np(1)*sh%nkpoints(2)*sh%nkpoints(3) + n(2)*sh%nkpoints(3) + np(3) + 1
        coord_cube_global(7,ik) = n(1)*sh%nkpoints(2)*sh%nkpoints(3) + np(2)*sh%nkpoints(3) + np(3) + 1
        coord_cube_global(8,ik) = np(1)*sh%nkpoints(2)*sh%nkpoints(3) + np(2)*sh%nkpoints(3) + np(3) + 1

        do i=1,3
          k%pos_cube(i,ik) = ( qproj(i) - dble(n(i))/dble(sh%nkpoints(i)) ) * dble(sh%nkpoints(i)) !  qproj(i) --> q(i)
        enddo

        do ii =1,8
          diff_kpoints(coord_cube_global(ii,ik)) = 1
        enddo

      enddo

      k%nk_smooth_loc = sum(diff_kpoints)

      if (sh%noncolin) then
        allocate(sh%beck_nc(sh%nkb,sh%npol,sh%ntot_e,k%nk_smooth_loc)) ! (nkb,npol,ntot_e,nk)
        allocate(sum_beck_nc_tmp(sh%nkb,sh%npol,sh%ntot_e,sh%nkpoints(3)))
      else
        allocate(sh%beckc(sh%nkb,sh%ntot_e,k%nk_smooth_loc))
        allocate(sum_beckc_tmp(sh%nkb,sh%ntot_e,sh%nkpoints(3)))
      endif
      if (sh%nonlocal_commutator) then
        allocate(sh%commut(sh%ntot_e,sh%ntot_e,3,k%nk_smooth_loc))
        allocate(sum_commut_tmp(sh%ntot_e,sh%ntot_e,3,sh%nkpoints(3)))
      endif

      counter = 0
      do ii = 1,nk_smooth
        if (diff_kpoints(ii) /= 0) then
          counter = counter + 1
          diff_kpoints(ii) = counter
        endif  
      enddo

      do ik = 1, k%nk_loc
        do ii =1,8
          k%coord_cube(ii,ik) = diff_kpoints(coord_cube_global(ii,ik))
        enddo
      enddo

      write(stdout,*)'simple_ip: opening file "hamiltonian_k"'
      if(ionode) then
         iun = find_free_unit()
         open( unit=iun, file=trim(tmp_dir)//trim(simpleip_in%prefix)//'.hamiltonian_k', status='old',form='unformatted')
         write(stdout,*)'File opened'
      endif   
      !
      write(stdout,*) 'Total number of k-point blocks: ' , sh%nkpoints(1)*sh%nkpoints(2)
      do ii=1,sh%nkpoints(1)
        do jj=1,sh%nkpoints(2)
          ! 
          write(stdout,*) 'k-point block: ' , jj  + (ii - 1)*sh%nkpoints(2)
          if (sh%noncolin) then
            if (ionode) read(iun) sum_beck_nc_tmp(1:sh%nkb,1:sh%npol,1:sh%ntot_e,1:sh%nkpoints(3))
          else
            if (ionode) read(iun) sum_beckc_tmp(1:sh%nkb,1:sh%ntot_e,1:sh%nkpoints(3))
          endif
          !
          if (sh%nonlocal_commutator) then
            if (ionode) read(iun) sum_commut_tmp(1:sh%ntot_e,1:sh%ntot_e,1:3,1:sh%nkpoints(3))  ! read commutator matrix
          endif
          !
          if (sh%noncolin) then
            call mp_bcast(sum_beck_nc_tmp,ionode_id,world_comm)
          else
            call mp_bcast(sum_beckc_tmp,ionode_id,world_comm)
          endif
          !
          if (sh%nonlocal_commutator) then
            call mp_bcast(sum_commut_tmp,ionode_id,world_comm)
          endif
          !
          do kk=1,sh%nkpoints(3)
            !
            ll = kk + (jj - 1)*sh%nkpoints(3) + (ii - 1)*sh%nkpoints(2)*sh%nkpoints(3)   ! global kindex
            !
            if (diff_kpoints(ll) /= 0) then
              if (sh%noncolin) then
                sh%beck_nc(1:sh%nkb,1:sh%npol,1:sh%ntot_e,diff_kpoints(ll)) = sum_beck_nc_tmp(1:sh%nkb,1:sh%npol,1:sh%ntot_e,kk)
              else
                sh%beckc(1:sh%nkb,1:sh%ntot_e,diff_kpoints(ll)) = sum_beckc_tmp(1:sh%nkb,1:sh%ntot_e,kk)
              endif
              if (sh%nonlocal_commutator) then
                  sh%commut(1:sh%ntot_e,1:sh%ntot_e,1:3,diff_kpoints(ll)) = sum_commut_tmp(1:sh%ntot_e,1:sh%ntot_e,1:3,kk) 
              endif
            endif      
            ! 
          enddo
          !
        enddo
      enddo
      
      if(ionode) then 
        close(iun)
        write(stdout,*)'File closed'
      endif

      if (sh%noncolin) then
        deallocate(sum_beck_nc_tmp) 
      else
        deallocate(sum_beckc_tmp)
      endif
      if (sh%nonlocal_commutator) then
        deallocate(sum_commut_tmp)
      endif
      deallocate(coord_cube_global,diff_kpoints)

    END SUBROUTINE read_shirley_k_interp


END MODULE simple_ip_objects
