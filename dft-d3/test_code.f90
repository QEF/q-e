!! using dftd3 API to get pressure by changing the alat


program test_code
 use dftd3_api
 use cell_base

 implicit none

    integer,parameter :: wp = kind(1.0d0)

    integer :: natoms, i, j
    real(wp) :: a0,ap,am,edisp1,edisp2
    real(wp), allocatable :: coords(:,:), coords0(:,:)
    real(wp), allocatable :: tau(:,:)
    real(wp), allocatable :: test_force(:,:)
    real(wp), dimension(3,3) :: latvecs, latvecs0
    character*2, allocatable :: typatom(:)
    type(dftd3_input) :: input
    type(dftd3_calc) :: dftd3
    integer,allocatable :: atnum(:)
    real(wp) :: edisp, Ep, Em
    real(wp) :: omega0, omega1, omega2
    real(wp), allocatable:: grads(:,:)
    real(wp), dimension(3,3) :: stress
    real(wp), parameter :: delta = 0.001
     real(wp) :: stress_num(3,3)


! ----- need to read: coordinates, atoms, lattice vectors
    read(*,*) a0
    read(*,*) latvecs(1,1), latvecs(2,1), latvecs(3,1)
    read(*,*) latvecs(1,2), latvecs(2,2), latvecs(3,2)
    read(*,*) latvecs(1,3), latvecs(2,3), latvecs(3,3)  
    latvecs0(:,:) = latvecs(:,:)
    latvecs(:,:) = latvecs(:,:) * a0
	a0=a0 ! * 1.8897261 !! API expects it in Ang, and auto-converts it to bohr.
                            !! The idea here was to have same input as for dftd3.f in
                            !! Grimme's original code. Now not used.
	latvecs(:,:) = latvecs(:,:) ! * 1.8897261
	call volume(a0, latvecs0(1,1), latvecs0(1,2), latvecs0(1,3), omega0 )

    read(*,*) natoms

    allocate(typatom(1:natoms))
    typatom(:) = 'none'
    allocate(coords(1:3, 1:natoms))
    allocate(coords0(1:3, 1:natoms))
    allocate(tau(1:3,1:natoms))
    allocate(test_force(1:3,1:natoms))
    allocate(atnum(1:natoms))
    allocate(grads(1:3, 1:natoms))

    do i=1, natoms
        read(*,*) typatom(i), coords(1,i), coords(2,i), coords(3,i)
    end do
    
	do i=1,natoms
		do j=1,3
			coords(j,i)=coords(1,i)*latvecs0(j,1)+coords(2,i)*latvecs0(j,2)+coords(3,i)*latvecs0(j,3)
		end do
	end do
	coords0(:,:) = coords(:,:) 	
	coords(:,:)=coords(:,:)*a0
write(*,*)'coords in units of a0'
write(*,"(3ES20.12)") coords

    stress(:,:)=0.0
	stress_num(:,:)=0.0

    input%threebody = .true.
	input%numgrad = .false.

    atnum(:) = get_atomic_number(typatom)
  
! initialize dftd3
    call dftd3_init(dftd3, input)

! choose functional
    call dftd3_set_functional(dftd3, func='pbe',version=3, tz=.false.)

    tau(:,:) = coords(:,:)  

write(*,*)'tau'
write(*,"(3ES20.12)") tau
write(*,*)
write(*,*) 'atomic numbers',atnum
write(*,*)
write(*,*) 'latvecs used for dftd3-API in units of a0'
write(*,"(3ES20.12)") latvecs
write(*,*)  
    call dftd3_pbc_dispersion(dftd3,tau,atnum,latvecs,edisp,grads,stress)
write(*,*)
write(*,*) 'edisp [Ry]',edisp*2.00
write(*,*)
write(*,*) 'forces from dftd3-API [Ha/au]'
write(*,"(3ES20.12)") grads

write(*,*)
write(*,*)"Forces calculated numerically: [Ha/au]"
    do j=1,natoms
        do i=1,3 ! x,y,z coordinates
            tau(:,:) = coords(:,:)         
            tau(i,j) = coords(i,j) - (delta/a0)      
            call dftd3_pbc_dispersion(dftd3,tau,atnum,latvecs,edisp1)
            tau(i,j) = coords(i,j) + (delta/a0)
            call dftd3_pbc_dispersion(dftd3,tau,atnum,latvecs,edisp2)
            test_force(i,j) = (edisp2 - edisp1)/(2.000000000000*delta/a0)
        end do 
write(*,*) test_force(1,j),test_force(2,j),test_force(3,j)
    end do 



write(*,*)
write(*,*) 'stress from dftd3-API [Ha/au^3]'
write(*,"(3ES20.12)") stress
write(*,*)


   do i=1,3
	do j=1,3
	tau(:,:)=coords0(:,:)*a0
	latvecs(:,:)=latvecs0(:,:) * a0
   	latvecs(j,i)=latvecs(j,i) + 0.0001
   	call dftd3_pbc_dispersion(dftd3,tau,atnum,latvecs,Ep)
   	latvecs(j,i)=latvecs(j,i) - 0.0002
   	call dftd3_pbc_dispersion(dftd3,tau,atnum,latvecs,Em)
	stress_num(j,i)=(Ep - Em) / (0.0002)
	end do
    end do
write(*,*) 'numerical stress [Ha/au^3]'
write(*,"(3ES20.12)") 1.0/omega0 *(-matmul(stress_num,transpose(latvecs)))
write(*,*)

    ap=a0+0.0001
    latvecs(:,:)=latvecs0(:,:)*ap
    tau(:,:) = coords0(:,:) * ap
    call dftd3_pbc_dispersion(dftd3,tau,atnum,latvecs,Ep)
    call volume(ap, latvecs0(1,1), latvecs0(1,2), latvecs0(1,3), omega2 )

    am=a0-0.0001
    latvecs(:,:)=latvecs0(:,:)*am
    tau(:,:) = coords0(:,:) * am
    call dftd3_pbc_dispersion(dftd3,tau,atnum,latvecs,Em)
    call volume(am, latvecs0(1,1), latvecs0(1,2), latvecs0(1,3), omega1 )
    
write(*,*) 'omega =',omega0, omega2, omega1
write(*,*) 'dE / dV =', 1.0 * (Ep-Em)/(omega2-omega1)
write(*,*)
    CALL recips(latvecs0(1,1),latvecs0(1,2),latvecs0(1,3),bg(1,1),bg(1,2),bg(1,3))
    bg(:,:)=bg(:,:)

    latvecs(:,:)=latvecs0(:,:) * a0
    latvecs(1,1)=latvecs(1,1) + 0.0001
    tau(:,:) = coords0(:,:) * a0
    call dftd3_pbc_dispersion(dftd3,tau,atnum,latvecs,Ep)
    latvecs(:,:)=latvecs0(:,:) * a0
    latvecs(1,1)=latvecs(1,1)-0.0001
    call dftd3_pbc_dispersion(dftd3,tau,atnum,latvecs,Em)
write(*,*) 'dE / dx =', (Ep - Em) / (0.0002)

    latvecs(:,:)=latvecs0(:,:) *a0
    latvecs(2,2)=latvecs(2,2) + 0.0001
    call dftd3_pbc_dispersion(dftd3,tau,atnum,latvecs,Ep)
    latvecs(:,:)=latvecs0(:,:) * a0
    latvecs(2,2)=latvecs(2,2)-0.0001
    call dftd3_pbc_dispersion(dftd3,tau,atnum,latvecs,Em)
write(*,*) 'dE / dy =', (Ep - Em) / (0.0002)

    latvecs(:,:)=latvecs0(:,:) * a0
    latvecs(3,3)=latvecs(3,3) + 0.0001
    call dftd3_pbc_dispersion(dftd3,tau,atnum,latvecs,Ep)
    latvecs(:,:)=latvecs0(:,:) * a0
    latvecs(3,3)=latvecs(3,3) - 0.0001
    call dftd3_pbc_dispersion(dftd3,tau,atnum,latvecs,Em)
write(*,*) 'dE / dz =', (Ep - Em) / (0.0002)

 contains
  function determinant(aa) result(det)
    real*8, intent(in) :: aa(:,:)
    real*8 :: det

    det = aa(1,1) * (aa(2,2) * aa(3,3) - aa(3,2) * aa(2,3))&
        & - aa(1,2) * (aa(2,1) * aa(3,3) - aa(3,1) * aa(2,3))&
        & + aa(1,3) * (aa(2,1) * aa(3,2) - aa(3,1) * aa(2,2))
    
  end function determinant


end program test_code



      subroutine inv(x,a) !x is normal lat, a is lat^(-1)
      IMPLICIT NONE
      real*8, intent(in)   :: x(3,3) !unitcell vectors in direct space
      real*8, intent(out)  :: a(3,3) !unitcell vectors in reciprocal space
      integer i
      real*8 det
      
      a=0.0
      det=x(1,1)*x(2,2)*x(3,3)+x(1,2)*x(2,3)*x(3,1)+x(1,3)*x(2,1)*&
          x(3,2)-x(1,3)*x(2,2)*x(3,1)-x(1,2)*x(2,1)*x(3,3)-x(1,1)*&
          x(2,3)*x(3,2)

      a(1,1)=x(2,2)*x(3,3)-x(2,3)*x(3,2)
      a(2,1)=x(2,3)*x(3,1)-x(2,1)*x(3,3)
      a(3,1)=x(2,1)*x(3,2)-x(2,2)*x(3,1)
      a(1,2)=x(1,3)*x(3,2)-x(1,2)*x(3,3)
      a(2,2)=x(1,1)*x(3,3)-x(1,3)*x(3,1)
      a(3,2)=x(1,2)*x(3,1)-x(1,1)*x(3,2)
      a(1,3)=x(1,2)*x(2,3)-x(1,3)*x(2,2)
      a(2,3)=x(1,3)*x(2,1)-x(1,1)*x(2,3)
      a(3,3)=x(1,1)*x(2,2)-x(1,2)*x(2,1)
      a=a/det
      end subroutine inv
  

