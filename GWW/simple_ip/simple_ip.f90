program simple_ip
  !
  USE start_end
  USE simple_ip_objects
  USE input_simple_ip
  !
  implicit none
  !
  TYPE(input_options_simple_ip) :: din
  TYPE(shirley) :: sh
  TYPE(kpoints) :: k
  TYPE(energies) :: e
  !
  call initialize_shirley(sh)
  call initialize_kpoints(k)
  call initialize_energies(e)
  !
  !setup MPI environment
  call startup
  !
  call start_clock('simple_ip')
  call start_clock('init (read)')
  !
  call read_input_simple_ip( din )  
  call read_shirley(din,sh) 
  call kgrid_creation(din,k,sh)
  call create_energies(sh,k,e)
  !
  if (din%nonlocal_interpolation) then
    call read_shirley_k_interp(din,sh,e,k)  ! trilinear interpolation
  else
    call read_shirley_k(din,sh,e)
  endif
  !
  call stop_clock('init (read)')
  !
  ! calculate IP dielectric function (interband + intraband)
  call dielectric(sh,din,k,e)
  !
  call stop_run
  call deallocate_shirley(sh)
  call deallocate_kpoints(k)
  call deallocate_energies(e)
  !
  call stop_clock('simple_ip')
  call print_clock('init (read)')
  call print_clock('diagonalization')
  call print_clock('diago_vnloc')
  call print_clock('diago_zheevx')
  call print_clock('optic_elements')
  call print_clock('dielectric')
  call print_clock('simple_ip')
  stop
  !
end program simple_ip

