!this subroutines calls solvers for eigen-value problem of the excitions

  SUBROUTINE simple_eigen(sin)

    USE input_simple_exc
    USE simple_objects
    USE kinds, ONLY : DP
    USE io_global, ONLY : stdout
    USE derived_objects

    implicit none

    TYPE(input_options) :: sin
    TYPE(exc), POINTER :: a(:)!excitons to be found
    TYPE(bands) :: bd
    TYPE(product) :: pd
    TYPE(potential) :: pt
    TYPE(prod_proj) :: pp
    TYPE(prod_mix) :: pm

    COMPLEX(kind=DP), ALLOCATABLE :: ene(:)!their energies
    INTEGER :: i

!read in excitonic Hamiltonian stuff
!    write(stdout,*) 'DEBUG 1'
   
    call read_bands(sin,bd)
!     write(stdout,*) 'DEBUG 2'
   
!read in product stuff
    call initialize_product(pd)
    call read_product(sin,pd)
!  write(stdout,*) 'DEBUG 3'
   

!read in potential stuff
    call initialize_potential(pt)
    call read_potential(sin,pt)
!  write(stdout,*) 'DEBUG 4'
   

!set up product contractions
    call initialize_prod_proj(pp)
    call build_prod_proj(bd,pd,pp)
!    write(stdout,*) 'DEBUG 5'
   
!set up mixed contractions

    call initialize_prod_mix(pm)
    call build_prod_mix(sin,bd,pd,pm,pt)
!    write(stdout,*) 'DEBUG 6'
   
    


    allocate(a(sin%nvec))
    do i=1,sin%nvec
       call setup_exc(bd,a(i))
    enddo

    allocate(ene(sin%nvec))
!call solver
    select case(sin%diago)
       case(0)!steepest descent
          call diago_exc_sd(sin,bd,pp,pt,pm,a)
       case(1)!CG
          call diago_exc_cg(sin,bd,pp,pt,pm,a)
    end select

    do i=1,sin%nvec
       call nice_write_exc(bd,sin,a(i),i)
    enddo

    do i=1,sin%nvec
       call deallocate_exc(a(i))
    enddo
    deallocate(a)
    deallocate(ene)
    call deallocate_bands(bd)
    call deallocate_product(pd)
    call deallocate_potential(pt)
    call deallocate_prod_proj(pp)
    call deallocate_prod_mix(pm)

    return
  END SUBROUTINE simple_eigen
