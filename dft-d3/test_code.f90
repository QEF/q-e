!! using dftd3 API to get pressure by changing the alat


program test_code
 use dftd3_api
 use mp_global, only: mp_startup, mp_global_end
 use io_global, only: ionode, ionode_id, stdout
 use mp_images, only: intra_image_comm
 use mp,        only: mp_bcast
 
 implicit none

    integer,parameter :: wp = kind(1.0d0)

    integer :: nat
    integer :: i, j
    real(wp) :: a0
    real(wp), allocatable :: coords(:,:)
    real(wp), allocatable :: tau(:,:)
    real(wp), allocatable :: test_force(:,:)
    real(wp), dimension(3,3) :: latvecs
    character*2 :: typatom
    type(dftd3_input) :: input
    type(dftd3_calc) :: dftd3
    integer,allocatable :: atnum(:)
    real(wp) :: edisp, edisp1, edisp2
    real(wp) :: omega0
    real(wp), allocatable:: grads(:,:)
    real(wp), dimension(3,3) :: stress
    real(wp), parameter :: delta1 = 0.001_wp, delta2 = 0.0001_wp
    real(wp) :: stress_num(3,3)
    real(wp) :: at(3,3)
    real(wp) :: bg(3,3)

    call mp_startup ()
    !
    if ( ionode ) then
       ! ----- read: lattice parameter a0 (a.u.), lattice vectors (a0 units)
       ! ----- read: number of atoms, atomic tyoe and positions (a0 units)
       read(*,*) a0
       read(*,*) at(1,1), at(2,1), at(3,1)
       read(*,*) at(1,2), at(2,2), at(3,2)
       read(*,*) at(1,3), at(2,3), at(3,3)  
       !
       read(*,*) nat
       !
    else
       open(unit=stdout,file='/dev/null')
    end if
    !
    call mp_bcast( nat,  ionode_id, intra_image_comm )
    !
    allocate(coords(1:3, 1:nat))
    allocate(atnum(1:nat))
    allocate(tau(1:3,1:nat))
    allocate(test_force(1:3,1:nat))
    allocate(grads(1:3, 1:nat))
    !       
    if ( ionode ) then
       do i=1, nat
          read(*,*) typatom, tau(1,i), tau(2,i), tau(3,i)
          atnum(i) = get_atomic_number(typatom)
       end do
    end if
    
    call mp_bcast( a0,      ionode_id, intra_image_comm )
    call mp_bcast( at, ionode_id, intra_image_comm )
    call mp_bcast( atnum,   ionode_id, intra_image_comm )
    call mp_bcast( tau,  ionode_id, intra_image_comm )

    stress(:,:)=0.0
    stress_num(:,:)=0.0

    input%threebody = .true.
    input%numgrad = .false.

! initialize dftd3
    call dftd3_init(dftd3, input)

! choose functional
    call dftd3_set_functional(dftd3, func='pbe',version=3, tz=.false.)

    latvecs(:,:) = at(:,:) * a0
    coords (:,:) = tau(:,:)* a0
    
    write(stdout,*) 'atomic numbers'
    write(stdout,'(20i4)') atnum
    write(stdout,*)    
    write(stdout,*) 'latvecs used for dftd3-API in a.u.'
    write(stdout,"(3ES20.12)") latvecs
    write(stdout,*)  
    write(stdout,*) 'coords used for dftd3-API in a.u.'
    write(stdout,"(3ES20.12)") coords

    call dftd3_pbc_dispersion( dftd3, coords, atnum, latvecs, edisp, grads, stress )
    
    write(stdout,*)
    write(stdout,*) 'edisp [Ry]',edisp*2.00
    write(stdout,*)
    write(stdout,*) 'forces from dftd3-API [Ha/au]'
    write(stdout,"(3ES20.12)") grads
    write(stdout,*)
    write(stdout,*)"forces calculated numerically: [Ha/au]"
    do j=1,nat
        do i=1,3 ! x,y,z coordinates
            coords(i,j) = (tau(i,j) - delta1)*a0
            call dftd3_pbc_dispersion( dftd3, coords, atnum, latvecs, edisp1 )
            coords(i,j) = (tau(i,j) +  delta1)*a0
            call dftd3_pbc_dispersion( dftd3, coords, atnum, latvecs, edisp2 )
            coords(i,j) = tau(i,j)*a0
            test_force(i,j) = (edisp2 - edisp1)/(2.0_wp*delta1*a0)
        end do 
        write(stdout,'(3es20.12)') test_force(1,j),test_force(2,j),test_force(3,j)
    end do 
    write(stdout,*)
    write(stdout,*) 'absolute difference'
    do j=1,nat
       write(stdout,'(3es15.6)') abs( grads(1,j) - test_force(1,j) ), &
            abs( grads(2,j) - test_force(2,j) ), abs( grads(3,j) - test_force(3,j) )
    end do
    !    
    write(stdout,*)
    write(stdout,*) 'stress from dftd3-API [Ha/au^3]'
    write(stdout,"(3ES20.12)") stress
    write(stdout,*)
    !
    ! bring atomic positions tau to crystal axis
    !
    call recips(at(1,1), at(1,2), at(1,3), bg(1,1), bg(1,2), bg(1,3) )
    call cryst_to_cart (nat, tau, bg, -1)
    call volume(a0, at(1,1), at(1,2), at(1,3), omega0 )
    !
    do i=1,3
       do j=1,3
          !
          ! bring atomic positions to cartesian axis with new crystal axis
          !
          latvecs(j,i) = (at(j,i) - delta2)*a0
          coords(:,:) = tau(:,:)
          call cryst_to_cart (nat, coords, latvecs, +1)
          ! latvecs is in a.u. above so on output coords is also in a.u.
          call dftd3_pbc_dispersion( dftd3, coords, atnum, latvecs, edisp1 )
          !
          latvecs(j,i) = (at(j,i) + delta2)*a0
          coords(:,:) = tau(:,:)
          call cryst_to_cart (nat, coords, latvecs, +1)
          ! as above
          call dftd3_pbc_dispersion( dftd3, coords, atnum, latvecs, edisp2 )
          !
          stress_num(j,i) = (edisp2 - edisp1) / (2.0_wp*delta2*a0)
          !
          ! reset latvecs(j,i) to its original value
          latvecs(j,i) = at(j,i)*a0
          !
       end do
    end do
    
    write(stdout,*) 'numerical stress [Ha/au^3]'
    ! write(stdout,"(3ES20.12)") -stress_num
    write(stdout,"(3ES20.12)") matmul(stress_num,transpose(latvecs))/omega0
    write(stdout,*)

10  call mp_global_end ()

end program test_code

