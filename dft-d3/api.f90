! Library interface for the dftd3 program.
!
! Copyright (C) 2016, Bálint Aradi
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 1, or (at your option)
! any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! GNU General Public License for more details.
!
! For the GNU General Public License, see <http://www.gnu.org/licenses/>
!

!> This module contains the API to access DFT-D3 functionality.
!!
module dftd3_api
  use dftd3_sizes
  use dftd3_common
  use dftd3_core
  implicit none
  private
  
  public :: dftd3_input, dftd3_calc
  public :: dftd3_init, dftd3_set_params, dftd3_set_functional
  public :: dftd3_dispersion, dftd3_pbc_dispersion
  public :: get_atomic_number


  !> Input for a dftd3 calculator.
  !!
  type :: dftd3_input
    !> Whether three body term should be calculated
    logical :: threebody = .false.

    !> Whether numerical gradients instead of analytical ones
    logical :: numgrad = .false.

    !> C6 min flags (or unallocated if not needed)
    logical, allocatable :: minc6list(:)

    !> C6 max flags (or unallocated if not needed)
    logical, allocatable :: maxc6list(:)
    
    !> Real space cutoff in atomic units.
    real(wp) :: cutoff = sqrt(9000.0_wp)

    !> Real space cutoff for coordination numbers in atomic units
    real(wp) :: cutoff_cn = sqrt(1600.0_wp)
  end type dftd3_input


  !> State of a dftd3 calculator.
  !!
  type :: dftd3_calc
!
! commented private attribute (QE 2016)
!    private
    logical :: noabc, numgrad
    integer :: version
    real(wp) :: s6, rs6, s18, rs18, alp
    real(wp) :: rthr, cn_thr
    integer  :: rep_vdw(3), rep_cn(3)
    real(wp), allocatable :: r0ab(:,:), c6ab(:,:,:,:,:)
    integer, allocatable :: mxc(:)
  end type dftd3_calc


contains

  !> Initializes a dftd3 calculator.
  !!
  !! \note You also need to call dftd3_set_functional() or dftd3_set_params()
  !!     before you can make an actual calculation.
  !!
  !! \param input  Input parameters for the calculator.
  !!
  subroutine dftd3_init(this, input)
! changed attribute out of "this" structure with inout (QE 2016)
    type(dftd3_calc), intent(inout) :: this
    type(dftd3_input), intent(in) :: input

    logical, allocatable :: minc6list(:), maxc6list(:)
    logical :: minc6, maxc6

    this%noabc = .not. input%threebody
    this%numgrad = input%numgrad

    allocate(minc6list(max_elem))
    if (allocated(input%minc6list)) then
      minc6list(:) = input%minc6list
    else
      minc6list(:) = .false.
    end if

    minc6 = any(minc6list)
    allocate(maxc6list(max_elem))
    if (allocated(input%maxc6list)) then
      maxc6list(:) = input%maxc6list
    else
      maxc6list(:) = .false.
    end if
    maxc6 = any(maxc6list)
    
    allocate(this%c6ab(max_elem, max_elem, maxc, maxc, 3))
    allocate(this%mxc(max_elem))
    call copyc6("", maxc, max_elem, this%c6ab, this%mxc, minc6, minc6list, &
        & maxc6, maxc6list)
    ! local variables deallocated (QE 2018)
    deallocate(maxc6list)
    deallocate(minc6list)
    this%rthr = input%cutoff**2
    this%cn_thr = input%cutoff_cn**2
    allocate(this%r0ab(max_elem, max_elem))
    call setr0ab(max_elem, autoang, this%r0ab)

  end subroutine dftd3_init
    
  !> Sets the parameter for the dftd3 calculator by choosing a functional.
  !!
  !! \param func  Name of the functional.
  !! \param version  Version to use.
  !! \param tz  Whether special TZ-parameters should be used.
  !!
  subroutine dftd3_set_functional(this, func, version, tz)
    type(dftd3_calc), intent(inout) :: this
    character(*), intent(in) :: func
    integer, intent(in) :: version
    logical, intent(in) :: tz

    this%version = version

    call setfuncpar(func, this%version, tz, this%s6, this%rs6, this%s18, &
        & this%rs18, this%alp)

  end subroutine dftd3_set_functional


  !> Sets the parameter for the dftd3 calculator directly.
  !!
  !! \param pars  Parameter to use. The 5 parameters must follow the same 
  !!     order as when specified in the dftd3.local file for the dftd3 program.
  !!     (see the documentation of the dftd3 program for details)
  !! \param version  Version to use. Note, that depending on the version the
  !!     five parameters may have different (or no) meaning.
  !!
  subroutine dftd3_set_params(this, pars, version)
    type(dftd3_calc), intent(inout) :: this
    real(wp), intent(in) :: pars(:)
    integer, intent(in) :: version

    if (size(pars) /= 5) then
      write(*,*) 'Invalid number of custom parameters'
      stop 1
    end if

    this%s6 = pars(1)
    this%rs6 = pars(2)
    this%s18 = pars(3)
    this%rs18 = pars(4)
    this%alp = pars(5)
    this%version = version
    
  end subroutine dftd3_set_params


  !> Calculates the dispersion for a given non-periodic configuration.
  !!
  !! \param coords  Coordinates of the atoms in atomic units. Shape: [3, nAtom].
  !! \param izp  Atomic number of each atom. Shape: [nAtom]. You can determine
  !!    the atomic number using the get_atomic_number() function.
  !! \param disp  Calculated dispersion energy in atomic units.
  !! \param grads  Calculated gradients in atomic units, if present.
  !!
  subroutine dftd3_dispersion(this, coords, izp, disp, grads)
    type(dftd3_calc), intent(in) :: this
    real(wp), intent(in) :: coords(:,:)
    integer, intent(in) :: izp(:)
    real(wp), intent(out) :: disp
    real(wp), optional, intent(out) :: grads(:,:)

    logical, allocatable :: fix(:)
    integer :: natom
    real(wp) :: s6, s18, rs6, rs8, rs10, alp6, alp8, alp10
    real(wp) :: e6, e8, e10, e12, e6abc, gdsp, gnorm

    natom = size(coords, dim=2)
    s6 = this%s6
    s18 = this%s18
    rs6 = this%rs6
    rs8 = this%rs18
    rs10 = this%rs18
    alp6 = this%alp
    alp8 = alp6 + 2.0_wp
    alp10 = alp8 + 2.0_wp
    call edisp(max_elem, maxc, natom, coords, izp, this%c6ab, this%mxc, &
        & r2r4, this%r0ab, rcov, rs6, rs8, rs10, alp6, alp8, alp10, &
        & this%version, this%noabc, this%rthr, this%cn_thr, e6, e8, e10, e12, &
        & e6abc)
    disp = -e6 * this%s6 - e8 * this%s18 - e6abc

    if (.not. present(grads)) then
      return
    end if

    allocate(fix(natom))
    fix(:) = .false.
    grads(:,:) = 0.0_wp
    call gdisp(max_elem, maxc, natom, coords, izp, this%c6ab, this%mxc, r2r4, &
        & this%r0ab, rcov, s6, s18, rs6, rs8, rs10, alp6, alp8, alp10, &
        & this%noabc, this%rthr, this%numgrad, this%version, .false., grads, &
        & gdsp, gnorm, this%cn_thr, fix)
    
  end subroutine dftd3_dispersion


  !> Calculates the dispersion for a given periodic configuration.
  !!
  !! \param coords  Coordinates of the atoms in atomic units. Shape: [3, nAtom].
  !! \param izp  Atomic number of each atom. Shape: [nAtom]. You can determine
  !!     the atomic number using the get_atomic_number() function.
  !! \param latvecs  Lattice vectors in atomic units. Shape: [3, 3].
  !! \param disp  Calculated dispersion energy in atomic units.
  !! \param grads  Calculated gradiens in atomic units, if present.
  !! \param stress  Calculated stress tensor in atomic units, if present.
  !!
  subroutine dftd3_pbc_dispersion(this, coords, izp, latvecs, disp, grads, &
      & stress)
    type(dftd3_calc), intent(in) :: this
    real(wp), intent(in) :: coords(:,:)
    integer, intent(in) :: izp(:)
    real(wp), intent(in) :: latvecs(:,:)
    real(wp), intent(out) :: disp
    real(wp), optional, intent(out) :: grads(:,:), stress(:,:)

    integer :: natom
    real(wp) :: s6, s18, rs6, rs8, rs10, alp6, alp8, alp10
    real(wp) :: e6, e8, e10, e12, e6abc, gnorm, disp2
    real(wp) :: rtmp3(3)
    integer :: rep_cn(3), rep_vdw(3)

    if (present(grads) .neqv. present(stress)) then
      write(*,*) "!!! Error in dftd3_pbc_dispersion"
      write(*,*) "Either both grads and stress must be present or none of them"
      stop
    end if

    natom = size(coords, dim=2)
    s6 = this%s6
    s18 = this%s18
    rs6 = this%rs6
    rs8 = this%rs18
    rs10 = this%rs18
    alp6 = this%alp
    alp8 = alp6 + 2.0_wp
    alp10 = alp8 + 2.0_wp

    call set_criteria(this%rthr, latvecs, rtmp3)
    rep_vdw(:) = int(rtmp3) + 1
    call set_criteria(this%cn_thr, latvecs, rtmp3)
    rep_cn(:) = int(rtmp3) + 1
    call pbcedisp(max_elem, maxc, natom, coords, izp, this%c6ab, this%mxc, &
        & r2r4, this%r0ab, rcov, rs6, rs8, rs10, alp6, alp8, alp10, &
        & this%version, this%noabc, e6, e8, e10, e12, e6abc, latvecs, &
        & this%rthr, rep_vdw, this%cn_thr, rep_cn)
    disp = -e6 * this%s6 - e8 * this%s18 - e6abc

    if (.not. present(grads)) then
      return
    end if

    grads(:,:) = 0.0_wp
    call pbcgdisp(max_elem, maxc, natom, coords, izp, this%c6ab, this%mxc, &
        & r2r4, this%r0ab, rcov, s6, s18, rs6, rs8, rs10, alp6, alp8, alp10, &
        & this%noabc, this%numgrad, this%version, grads, disp2, gnorm, &
        & stress, latvecs, rep_vdw, rep_cn, this%rthr, .false., this%cn_thr)
    ! Note, the stress variable in pbcgdisp contains the *lattice derivatives*
    ! on return, so it needs to be converted to obtain the stress tensor.
    stress(:,:) = -matmul(stress, transpose(latvecs))&
        & / abs(determinant(latvecs))
    
  end subroutine dftd3_pbc_dispersion


  !> Returns the atomic number for a given species.
  !!
  !! \param species  Chemical symbol of the species.
  !! \return  Atomic number.
  !!
  elemental function get_atomic_number(species) result(izp)
    character(*), intent(in) :: species
    integer :: izp

    call elem(trim(species), izp)

  end function get_atomic_number

  
end module dftd3_api
