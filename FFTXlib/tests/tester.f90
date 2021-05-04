! This file is part of fortran_tester
! Copyright 2015 Pierre de Buyl and Peter Colberg
!           2016 Pierre de Buyl and Stefano Szaghi
!           2018 Pierre de Buyl and Pietro Bonfa
! License: BSD

!> Routines to test Fortran programs
!!
!! fortran_tester is a pure-Fortran module. It provides a datatype to hold test results and
!! routines to test for equality, closeness, and positivity of variables. The routines are
!! overloaded and the resulting interface consists of a small number of names.

module tester
  use, intrinsic :: iso_fortran_env, only : int8, int16, int32, int64, real32, real64

  implicit none
  private
  public :: tester_t

  !> The main **tester** class.
  type :: tester_t
     integer(int32) :: n_errors=0_int32                         !< Number of errors.
     integer(int32) :: n_tests=0_int32                          !< Number of tests.
     real(real32)   :: tolerance32=2._real32*epsilon(1._real32) !< Real tolerance, 32 bits.
     real(real64)   :: tolerance64=2._real64*epsilon(1._real64) !< Real tolerance, 64 bits.
   contains
     procedure :: init                           !< Initialize the tester.
     procedure :: print                          !< Print tests results.
     generic, public :: assert_equal =>     &
                        assert_equal_i8,    &
                        assert_equal_i16,   &
                        assert_equal_i32,   &
                        assert_equal_i64,   &
                        assert_equal_r32,   &
                        assert_equal_r64,   &
                        assert_equal_c32,   &
                        assert_equal_c64,   &
                        assert_equal_l,     &
                        assert_equal_i8_1,  &
                        assert_equal_i16_1, &
                        assert_equal_i32_1, &
                        assert_equal_i64_1, &
                        assert_equal_r32_1, &
                        assert_equal_r64_1, &
                        assert_equal_c32_1, &
                        assert_equal_c64_1, &
                        assert_equal_l_1         !< Check if two values (integer, real, complex or logical) are equal.
     procedure, private :: assert_equal_i8       !< Check if two integers (8  bits) are equal.
     procedure, private :: assert_equal_i16      !< Check if two integers (16 bits) are equal.
     procedure, private :: assert_equal_i32      !< Check if two integers (32 bits) are equal.
     procedure, private :: assert_equal_i64      !< Check if two integers (64 bits) are equal.
     procedure, private :: assert_equal_r32      !< Check if two reals (32 bits) are equal.
     procedure, private :: assert_equal_r64      !< Check if two reals (64 bits) are equal.
     procedure, private :: assert_equal_c32      !< Check if two complex numbers (32 bits) are equal.
     procedure, private :: assert_equal_c64      !< Check if two complex numbers (64 bits) are equal.
     procedure, private :: assert_equal_l        !< Check if two logicals are equal.
     procedure, private :: assert_equal_i8_1     !< Check if two integer (8  bits) arrays (rank 1) are equal.
     procedure, private :: assert_equal_i16_1    !< Check if two integer (16 bits) arrays (rank 1) are equal.
     procedure, private :: assert_equal_i32_1    !< Check if two integer (32 bits) arrays (rank 1) are equal.
     procedure, private :: assert_equal_i64_1    !< Check if two integer (64 bits) arrays (rank 1) are equal.
     procedure, private :: assert_equal_r32_1    !< Check if two real (32 bits) arrays (rank 1) are equal.
     procedure, private :: assert_equal_r64_1    !< Check if two real (64 bits) arrays (rank 1) are equal.
     procedure, private :: assert_equal_c32_1    !< Check if two complex (32 bits) arrays (rank 1) are equal.
     procedure, private :: assert_equal_c64_1    !< Check if two complex (64 bits) arrays (rank 1) are equal.
     procedure, private :: assert_equal_l_1      !< Check if two logical arrays (rank 1) are equal.
     generic, public :: assert_positive =>     &
                        assert_positive_i8,    &
                        assert_positive_i16,   &
                        assert_positive_i32,   &
                        assert_positive_i64,   &
                        assert_positive_r32,   &
                        assert_positive_r64,   &
                        assert_positive_i8_1,  &
                        assert_positive_i16_1, &
                        assert_positive_i32_1, &
                        assert_positive_i64_1, &
                        assert_positive_r32_1, &
                        assert_positive_r64_1    !< Check if a number (integer or real) is positive.
     procedure, private :: assert_positive_i8    !< Check if a integer (8  bits) is positive.
     procedure, private :: assert_positive_i16   !< Check if a integer (16 bits) is positive.
     procedure, private :: assert_positive_i32   !< Check if a integer (32 bits) is positive.
     procedure, private :: assert_positive_i64   !< Check if a integer (64 bits) is positive.
     procedure, private :: assert_positive_r32   !< Check if a real (32 bits) is positive.
     procedure, private :: assert_positive_r64   !< Check if a real (64 bits) is positive.
     procedure, private :: assert_positive_i8_1  !< Check if a integer (8  bits) array (rank 1) is positive.
     procedure, private :: assert_positive_i16_1 !< Check if a integer (16 bits) array (rank 1) is positive.
     procedure, private :: assert_positive_i32_1 !< Check if a integer (32 bits) array (rank 1) is positive.
     procedure, private :: assert_positive_i64_1 !< Check if a integer (64 bits) array (rank 1) is positive.
     procedure, private :: assert_positive_r32_1 !< Check if a real (32 bits) array (rank 1) is positive.
     procedure, private :: assert_positive_r64_1 !< Check if a real (64 bits) array (rank 1) is positive.
     generic, public :: assert_close =>     &
                        assert_close_r32,   &
                        assert_close_r64,   &
                        assert_close_c32,   &
                        assert_close_c64,   &
                        assert_close_r32_1, &
                        assert_close_r64_1, &
                        assert_close_c32_1, &
                        assert_close_c64_1       !< Check if two values (real or complex) are close with respect a tolerance.
     procedure, private :: assert_close_r32      !< Check if two reals (32 bits) are close with respect a tolerance.
     procedure, private :: assert_close_r64      !< Check if two reals (64 bits) are close with respect a tolerance.
     procedure, private :: assert_close_c32      !< Check if two complex numbers (32 bits) are close with respect a tolerance.
     procedure, private :: assert_close_c64      !< Check if two complex numbers (64 bits) are close with respect a tolerance.
     procedure, private :: assert_close_r32_1    !< Check if two real (32 bits) arrays (rank 1) are close with respect a tolerance.
     procedure, private :: assert_close_r64_1    !< Check if two real (64 bits) arrays (rank 1) are close with respect a tolerance.
     procedure, private :: assert_close_c32_1    !< Check if two complex (32 bits) arrays (rank 1) are close with respect a tolerance.
     procedure, private :: assert_close_c64_1    !< Check if two complex (64 bits) arrays (rank 1) are close with respect a tolerance.
  end type tester_t

contains

  !> Initialize the tester.
  subroutine init(this, tolerance32, tolerance64)
    class(tester_t), intent(out)          :: this        !< The tester.
    real(real32),    intent(in), optional :: tolerance32 !< Real tolerance, 32 bits.
    real(real64),    intent(in), optional :: tolerance64 !< Real tolerance, 64 bits.

    this% n_errors = 0
    this% n_tests = 0

    if (present(tolerance64)) then
       this% tolerance64 = tolerance64
    else
       this% tolerance64 = 2._real64*epsilon(1._real64)
    end if

    if (present(tolerance32)) then
       this% tolerance32 = tolerance32
    else
       this% tolerance32 = 2._real32*epsilon(1._real32)
    end if

  end subroutine init

  !> Print tests results.
  subroutine print(this, errorstop)
    class(tester_t), intent(in)           :: this      !< The tester.
    logical,         intent(in), optional :: errorstop !< Flag to activate error stop if one test fails.

    logical :: do_errorstop
    if (present(errorstop)) then
       do_errorstop = errorstop
    else
       do_errorstop = .true.
    end if

    write(*,*) 'fortran_tester:', this% n_errors, ' error(s) for', this% n_tests, 'test(s)'

    if (this% n_errors == 0) then
       write(*,*) 'fortran_tester: all tests succeeded'
    else
       write(*,*) 'fortran_tester: tests failed'
       if (do_errorstop) then
          stop 1
       end if
    end if

  end subroutine print

  !> Check if two integers (8 bits) are equal.
  subroutine assert_equal_i8(this, i1, i2, fail)
    class(tester_t), intent(inout)        :: this !< The tester.
    integer(int8),   intent(in)           :: i1   !< Value to compare.
    integer(int8),   intent(in)           :: i2   !< Value to compare.
    logical,         intent(in), optional :: fail !< Fail flag.

    this% n_tests = this% n_tests + 1
    if (i1 .ne. i2) then
       if (.not. present(fail) .or. (present(fail) .and. fail .eqv. .false.)) then
          this% n_errors = this% n_errors + 1
       end if
    end if

  end subroutine assert_equal_i8

  !> Check if two integers (16 bits) are equal.
  subroutine assert_equal_i16(this, i1, i2, fail)
    class(tester_t), intent(inout)        :: this   !< The tester.
    integer(int16),  intent(in)           :: i1     !< Value to compare.
    integer(int16),  intent(in)           :: i2     !< Value to compare.
    logical,         intent(in), optional :: fail   !< Fail flag.

    this% n_tests = this% n_tests + 1
    if (i1 .ne. i2) then
       if (.not. present(fail) .or. (present(fail) .and. fail .eqv. .false.)) then
          this% n_errors = this% n_errors + 1
       end if
    end if

  end subroutine assert_equal_i16

  !> Check if two integers (32 bits) are equal.
  subroutine assert_equal_i32(this, i1, i2, fail)
    class(tester_t), intent(inout)        :: this !< The tester.
    integer(int32),  intent(in)           :: i1   !< Value to compare.
    integer(int32),  intent(in)           :: i2   !< Value to compare.
    logical,         intent(in), optional :: fail !< Fail flag.

    this% n_tests = this% n_tests + 1
    if (i1 .ne. i2) then
       if (.not. present(fail) .or. (present(fail) .and. fail .eqv. .false.)) then
          this% n_errors = this% n_errors + 1
       end if
    end if

  end subroutine assert_equal_i32

  !> Check if two integers (64 bits) are equal.
  subroutine assert_equal_i64(this, i1, i2, fail)
    class(tester_t), intent(inout)        :: this !< The tester.
    integer(int64),  intent(in)           :: i1   !< Value to compare.
    integer(int64),  intent(in)           :: i2   !< Value to compare.
    logical,         intent(in), optional :: fail !< Fail flag.

    this% n_tests = this% n_tests + 1
    if (i1 .ne. i2) then
       if (.not. present(fail) .or. (present(fail) .and. fail .eqv. .false.)) then
          this% n_errors = this% n_errors + 1
       end if
    end if

  end subroutine assert_equal_i64

  !> Check if two reals (32 bits) are equal.
  subroutine assert_equal_r32(this, r1, r2, fail)
    class(tester_t), intent(inout)        :: this !< The tester.
    real(real32),    intent(in)           :: r1   !< Value to compare.
    real(real32),    intent(in)           :: r2   !< Value to compare.
    logical,         intent(in), optional :: fail !< Fail flag.

    this% n_tests = this% n_tests + 1
    if (r1 .ne. r2) then
       if (.not. present(fail) .or. (present(fail) .and. fail .eqv. .false.)) then
          this% n_errors = this% n_errors + 1
       end if
    end if

  end subroutine assert_equal_r32

  !> Check if two reals (64 bits) are equal.
  subroutine assert_equal_r64(this, r1, r2, fail)
    class(tester_t), intent(inout)        :: this !< The tester.
    real(real64),    intent(in)           :: r1   !< Value to compare.
    real(real64),    intent(in)           :: r2   !< Value to compare.
    logical,         intent(in), optional :: fail !< Fail flag.

    this% n_tests = this% n_tests + 1
    if (r1 .ne. r2) then
       if (.not. present(fail) .or. (present(fail) .and. fail .eqv. .false.)) then
          this% n_errors = this% n_errors + 1
       end if
    end if

  end subroutine assert_equal_r64

  !> Check if two complex numbers (32 bits) are equal.
  subroutine assert_equal_c32(this, c1, c2, fail)
    class(tester_t), intent(inout)           :: this !< The tester.
    complex(real32),    intent(in)           :: c1   !< Value to compare.
    complex(real32),    intent(in)           :: c2   !< Value to compare.
    logical,            intent(in), optional :: fail !< Fail flag.

    this% n_tests = this% n_tests + 1
    if (c1 .ne. c2) then
       if (.not. present(fail) .or. (present(fail) .and. fail .eqv. .false.)) then
          this% n_errors = this% n_errors + 1
       end if
    end if

  end subroutine assert_equal_c32

  !> Check if two complex numbers (64 bits) are equal.
  subroutine assert_equal_c64(this, c1, c2, fail)
    class(tester_t), intent(inout)           :: this !< The tester.
    complex(real64),    intent(in)           :: c1   !< Value to compare.
    complex(real64),    intent(in)           :: c2   !< Value to compare.
    logical,            intent(in), optional :: fail !< Fail flag.

    this% n_tests = this% n_tests + 1
    if (c1 .ne. c2) then
       if (.not. present(fail) .or. (present(fail) .and. fail .eqv. .false.)) then
          this% n_errors = this% n_errors + 1
       end if
    end if

  end subroutine assert_equal_c64

  !> Check if two logicals are equal.
 subroutine assert_equal_l(this, l1, l2, fail)
    class(tester_t), intent(inout)        :: this !< The tester.
    logical,         intent(in)           :: l1   !< Value to compare.
    logical,         intent(in)           :: l2   !< Value to compare.
    logical,         intent(in), optional :: fail !< Fail flag.

    this% n_tests = this% n_tests + 1
    if (l1 .neqv. l2) then
       if (.not. present(fail) .or. (present(fail) .and. fail .eqv. .false.)) then
          this% n_errors = this% n_errors + 1
       end if
    end if

  end subroutine assert_equal_l

  !> Check if two integer (8 bits) arrays (rank 1) are equal.
  subroutine assert_equal_i8_1(this, i1, i2, fail)
    class(tester_t),             intent(inout)        :: this !< The tester.
    integer(int8), dimension(:), intent(in)           :: i1   !< Value to compare.
    integer(int8), dimension(:), intent(in)           :: i2   !< Value to compare.
    logical,                     intent(in), optional :: fail !< Fail flag.

    this% n_tests = this% n_tests + 1

    if ( size(i1) .ne. size(i2) ) then
       if (.not. present(fail) .or. (present(fail) .and. fail .eqv. .false.)) then
          this% n_errors = this% n_errors + 1
       end if
    else
       if ( maxval(abs(i1-i2)) > 0 ) then
          if (.not. present(fail) .or. (present(fail) .and. fail .eqv. .false.)) then
             this% n_errors = this% n_errors + 1
          end if
       end if
    end if

  end subroutine assert_equal_i8_1

  !> Check if two integer (16 bits) arrays (rank 1) are equal.
  subroutine assert_equal_i16_1(this, i1, i2, fail)
    class(tester_t),              intent(inout)        :: this !< The tester.
    integer(int16), dimension(:), intent(in)           :: i1   !< Value to compare.
    integer(int16), dimension(:), intent(in)           :: i2   !< Value to compare.
    logical,                      intent(in), optional :: fail !< Fail flag.

    this% n_tests = this% n_tests + 1

    if ( size(i1) .ne. size(i2) ) then
       if (.not. present(fail) .or. (present(fail) .and. fail .eqv. .false.)) then
          this% n_errors = this% n_errors + 1
       end if
    else
       if ( maxval(abs(i1-i2)) > 0 ) then
          if (.not. present(fail) .or. (present(fail) .and. fail .eqv. .false.)) then
             this% n_errors = this% n_errors + 1
          end if
       end if
    end if

  end subroutine assert_equal_i16_1

  !> Check if two integer (32 bits) arrays (rank 1) are equal.
  subroutine assert_equal_i32_1(this, i1, i2, fail)
    class(tester_t),              intent(inout)        :: this !< The tester.
    integer(int32), dimension(:), intent(in)           :: i1   !< Value to compare.
    integer(int32), dimension(:), intent(in)           :: i2   !< Value to compare.
    logical,                      intent(in), optional :: fail !< Fail flag.

    this% n_tests = this% n_tests + 1

    if ( size(i1) .ne. size(i2) ) then
       if (.not. present(fail) .or. (present(fail) .and. fail .eqv. .false.)) then
          this% n_errors = this% n_errors + 1
       end if
    else
       if ( maxval(abs(i1-i2)) > 0 ) then
          if (.not. present(fail) .or. (present(fail) .and. fail .eqv. .false.)) then
             this% n_errors = this% n_errors + 1
          end if
       end if
    end if

  end subroutine assert_equal_i32_1

  !> Check if two integer (64 bits) arrays (rank 1) are equal.
  subroutine assert_equal_i64_1(this, i1, i2, fail)
    class(tester_t),              intent(inout)        :: this !< The tester.
    integer(int64), dimension(:), intent(in)           :: i1   !< Value to compare.
    integer(int64), dimension(:), intent(in)           :: i2   !< Value to compare.
    logical,                      intent(in), optional :: fail !< Fail flag.

    this% n_tests = this% n_tests + 1

    if ( size(i1) .ne. size(i2) ) then
       if (.not. present(fail) .or. (present(fail) .and. fail .eqv. .false.)) then
          this% n_errors = this% n_errors + 1
       end if
    else
       if ( maxval(abs(i1-i2)) > 0 ) then
          if (.not. present(fail) .or. (present(fail) .and. fail .eqv. .false.)) then
             this% n_errors = this% n_errors + 1
          end if
       end if
    end if

  end subroutine assert_equal_i64_1

  !> Check if two real (32 bits) arrays (rank 1) are equal.
  subroutine assert_equal_r32_1(this, r1, r2, fail)
    class(tester_t),            intent(inout)        :: this !< The tester.
    real(real32), dimension(:), intent(in)           :: r1   !< Value to compare.
    real(real32), dimension(:), intent(in)           :: r2   !< Value to compare.
    logical,                    intent(in), optional :: fail !< Fail flag.

    this% n_tests = this% n_tests + 1

    if ( size(r1) .ne. size(r2) ) then
       if (.not. present(fail) .or. (present(fail) .and. fail .eqv. .false.)) then
          this% n_errors = this% n_errors + 1
       end if
    else
       if ( maxval(abs(r1-r2)) > 0 ) then
          if (.not. present(fail) .or. (present(fail) .and. fail .eqv. .false.)) then
             this% n_errors = this% n_errors + 1
          end if
       end if
    end if

  end subroutine assert_equal_r32_1

  !> Check if two real (64 bits) arrays (rank 1) are equal.
  subroutine assert_equal_r64_1(this, r1, r2, fail)
    class(tester_t),            intent(inout)        :: this !< The tester.
    real(real64), dimension(:), intent(in)           :: r1   !< Value to compare.
    real(real64), dimension(:), intent(in)           :: r2   !< Value to compare.
    logical,                    intent(in), optional :: fail !< Fail flag.

    this% n_tests = this% n_tests + 1

    if ( size(r1) .ne. size(r2) ) then
       if (.not. present(fail) .or. (present(fail) .and. fail .eqv. .false.)) then
          this% n_errors = this% n_errors + 1
       end if
    else
       if ( maxval(abs(r1-r2)) > 0 ) then
          if (.not. present(fail) .or. (present(fail) .and. fail .eqv. .false.)) then
             this% n_errors = this% n_errors + 1
          end if
       end if
    end if

  end subroutine assert_equal_r64_1
  
  !> Check if two complex (32 bits) arrays (rank 1) are equal.
  subroutine assert_equal_c32_1(this, c1, c2, fail)
    class(tester_t),               intent(inout)        :: this !< The tester.
    complex(real32), dimension(:), intent(in)           :: c1   !< Value to compare.
    complex(real32), dimension(:), intent(in)           :: c2   !< Value to compare.
    logical,                       intent(in), optional :: fail !< Fail flag.

    this% n_tests = this% n_tests + 1

    if ( size(c1) .ne. size(c2) ) then
       if (.not. present(fail) .or. (present(fail) .and. fail .eqv. .false.)) then
          this% n_errors = this% n_errors + 1
       end if
    else
       if ( maxval(abs(c1-c2)) > 0 ) then
          if (.not. present(fail) .or. (present(fail) .and. fail .eqv. .false.)) then
             this% n_errors = this% n_errors + 1
          end if
       end if
    end if

  end subroutine assert_equal_c32_1

  !> Check if two complex (64 bits) arrays (rank 1) are equal.
  subroutine assert_equal_c64_1(this, c1, c2, fail)
    class(tester_t),               intent(inout)        :: this !< The tester.
    complex(real64), dimension(:), intent(in)           :: c1   !< Value to compare.
    complex(real64), dimension(:), intent(in)           :: c2   !< Value to compare.
    logical,                       intent(in), optional :: fail !< Fail flag.

    this% n_tests = this% n_tests + 1

    if ( size(c1) .ne. size(c2) ) then
       if (.not. present(fail) .or. (present(fail) .and. fail .eqv. .false.)) then
          this% n_errors = this% n_errors + 1
       end if
    else
       if ( maxval(abs(c1-c2)) > 0 ) then
          if (.not. present(fail) .or. (present(fail) .and. fail .eqv. .false.)) then
             this% n_errors = this% n_errors + 1
          end if
       end if
    end if

  end subroutine assert_equal_c64_1

  !> Check if two logical arrays (rank 1) are equal.
  subroutine assert_equal_l_1(this, l1, l2, fail)
    class(tester_t), intent(inout)            :: this !< The tester.
    logical,         intent(in), dimension(:) :: l1   !< Value to compare.
    logical,         intent(in), dimension(:) :: l2   !< Value to compare.
    logical,         intent(in), optional     :: fail !< Fail flag.

    integer :: k

    this% n_tests = this% n_tests + 1

    if ( size(l1) .ne. size(l2) ) then
       if (.not. present(fail) .or. (present(fail) .and. fail .eqv. .false.)) then
          this% n_errors = this% n_errors + 1
       end if
    else
       do k = 1, size(l1)
          if (l1(k) .neqv. l2(k)) then
             if (.not. present(fail) .or. (present(fail) .and. fail .eqv. .false.)) then
                this% n_errors = this% n_errors + 1
             end if
             exit
          end if
       end do
    end if

  end subroutine assert_equal_l_1

  !> Check if a integer (32 bits) is positive.
  subroutine assert_positive_i8(this, i, fail)
    class(tester_t), intent(inout)        :: this !< The tester.
    integer(int8),   intent(in)           :: i    !< Value to check.
    logical,         intent(in), optional :: fail !< Fail flag.

    this% n_tests = this% n_tests + 1
    if (i < 0) then
       if (.not. present(fail) .or. (present(fail) .and. fail .eqv. .false.)) then
          this% n_errors = this% n_errors + 1
       end if
    end if

  end subroutine assert_positive_i8

  !> Check if a integer (16 bits) is positive.
  subroutine assert_positive_i16(this, i, fail)
    class(tester_t), intent(inout)        :: this !< The tester.
    integer(int16),  intent(in)           :: i    !< Value to check.
    logical,         intent(in), optional :: fail !< Fail flag.

    this% n_tests = this% n_tests + 1
    if (i < 0) then
       if (.not. present(fail) .or. (present(fail) .and. fail .eqv. .false.)) then
          this% n_errors = this% n_errors + 1
       end if
    end if

  end subroutine assert_positive_i16

  !> Check if a integer (32 bits) is positive.
  subroutine assert_positive_i32(this, i, fail)
    class(tester_t), intent(inout)        :: this !< The tester.
    integer(int32),  intent(in)           :: i    !< Value to check.
    logical,         intent(in), optional :: fail !< Fail flag.

    this% n_tests = this% n_tests + 1
    if (i < 0) then
       if (.not. present(fail) .or. (present(fail) .and. fail .eqv. .false.)) then
          this% n_errors = this% n_errors + 1
       end if
    end if

  end subroutine assert_positive_i32

  !> Check if a integer (32 bits) is positive.
  subroutine assert_positive_i64(this, i, fail)
    class(tester_t), intent(inout)        :: this !< The tester.
    integer(int64),  intent(in)           :: i    !< Value to check.
    logical,         intent(in), optional :: fail !< Fail flag.

    this% n_tests = this% n_tests + 1
    if (i < 0) then
       if (.not. present(fail) .or. (present(fail) .and. fail .eqv. .false.)) then
          this% n_errors = this% n_errors + 1
       end if
    end if

  end subroutine assert_positive_i64

  !> Check if a real (32 bits) is positive.
  subroutine assert_positive_r32(this, r, fail)
    class(tester_t), intent(inout)        :: this !< The tester.
    real(real32),    intent(in)           :: r    !< Value to check.
    logical,         intent(in), optional :: fail !< Fail flag.

    this% n_tests = this% n_tests + 1
    if (r < 0) then
       if (.not. present(fail) .or. (present(fail) .and. fail .eqv. .false.)) then
          this% n_errors = this% n_errors + 1
       end if
    end if

  end subroutine assert_positive_r32

  !> Check if a real (64 bits) is positive.
  subroutine assert_positive_r64(this, r, fail)
    class(tester_t), intent(inout)        :: this !< The tester.
    real(real64),    intent(in)           :: r    !< Value to check.
    logical,         intent(in), optional :: fail !< Fail flag.

    this% n_tests = this% n_tests + 1
    if (r < 0) then
       if (.not. present(fail) .or. (present(fail) .and. fail .eqv. .false.)) then
          this% n_errors = this% n_errors + 1
       end if
    end if

  end subroutine assert_positive_r64

  !> Check if a integer (8 bits) array (rank 1) is positive.
  subroutine assert_positive_i8_1(this, i, fail)
    class(tester_t),             intent(inout)        :: this !< The tester.
    integer(int8), dimension(:), intent(in)           :: i    !< Value to check.
    logical,                     intent(in), optional :: fail !< Fail flag.

    this% n_tests = this% n_tests + 1

    if ( minval(i) < 0 ) then
       if (.not. present(fail) .or. (present(fail) .and. fail .eqv. .false.)) then
          this% n_errors = this% n_errors + 1
       end if
    end if

  end subroutine assert_positive_i8_1

  !> Check if a integer (16 bits) array (rank 1) is positive.
  subroutine assert_positive_i16_1(this, i, fail)
    class(tester_t),              intent(inout)        :: this !< The tester.
    integer(int16), dimension(:), intent(in)           :: i    !< Value to check.
    logical,                      intent(in), optional :: fail !< Fail flag.

    this% n_tests = this% n_tests + 1

    if ( minval(i) < 0 ) then
       if (.not. present(fail) .or. (present(fail) .and. fail .eqv. .false.)) then
          this% n_errors = this% n_errors + 1
       end if
    end if

  end subroutine assert_positive_i16_1

  !> Check if a integer (32 bits) array (rank 1) is positive.
  subroutine assert_positive_i32_1(this, i, fail)
    class(tester_t),              intent(inout)        :: this !< The tester.
    integer(int32), dimension(:), intent(in)           :: i    !< Value to check.
    logical,                      intent(in), optional :: fail !< Fail flag.

    this% n_tests = this% n_tests + 1

    if ( minval(i) < 0 ) then
       if (.not. present(fail) .or. (present(fail) .and. fail .eqv. .false.)) then
          this% n_errors = this% n_errors + 1
       end if
    end if

  end subroutine assert_positive_i32_1

  !> Check if a integer (64 bits) array (rank 1) is positive.
  subroutine assert_positive_i64_1(this, i, fail)
    class(tester_t),              intent(inout)        :: this !< The tester.
    integer(int64), dimension(:), intent(in)           :: i    !< Value to check.
    logical,                      intent(in), optional :: fail !< Fail flag.

    this% n_tests = this% n_tests + 1

    if ( minval(i) < 0 ) then
       if (.not. present(fail) .or. (present(fail) .and. fail .eqv. .false.)) then
          this% n_errors = this% n_errors + 1
       end if
    end if

  end subroutine assert_positive_i64_1

  !> Check if a real (32 bits) array (rank 1) is positive.
  subroutine assert_positive_r32_1(this, r, fail)
    class(tester_t),            intent(inout)        :: this !< The tester.
    real(real32), dimension(:), intent(in)           :: r    !< Value to check.
    logical,                    intent(in), optional :: fail !< Fail flag.

    this% n_tests = this% n_tests + 1

    if ( minval(r) < 0 ) then
       if (.not. present(fail) .or. (present(fail) .and. fail .eqv. .false.)) then
          this% n_errors = this% n_errors + 1
       end if
    end if

  end subroutine assert_positive_r32_1

  !> Check if a real (64 bits) array (rank 1) is positive.
  subroutine assert_positive_r64_1(this, r, fail)
    class(tester_t),            intent(inout)        :: this !< The tester.
    real(real64), dimension(:), intent(in)           :: r    !< Value to check.
    logical,                    intent(in), optional :: fail !< Fail flag.

    this% n_tests = this% n_tests + 1

    if ( minval(r) < 0 ) then
       if (.not. present(fail) .or. (present(fail) .and. fail .eqv. .false.)) then
          this% n_errors = this% n_errors + 1
       end if
    end if

  end subroutine assert_positive_r64_1

  !> Check if two reals (32 bits) are close with respect a tolerance.
  subroutine assert_close_r32(this, r1, r2, fail)
    class(tester_t), intent(inout)        :: this !< The tester.
    real(real32),    intent(in)           :: r1   !< Value to compare.
    real(real32),    intent(in)           :: r2   !< Value to compare.
    logical,         intent(in), optional :: fail !< Fail flag.

    this% n_tests = this% n_tests + 1

    if ( abs(r1-r2) > this% tolerance32 ) then
       if (.not. present(fail) .or. (present(fail) .and. fail .eqv. .false.)) then
          this% n_errors = this% n_errors + 1
       end if
    end if

  end subroutine assert_close_r32

  !> Check if two reals (64 bits) are close with respect a tolerance.
  subroutine assert_close_r64(this, r1, r2, fail)
    class(tester_t),  intent(inout)        :: this !< The tester.
    real(real64),     intent(in)           :: r1   !< Value to compare.
    real(real64),     intent(in)           :: r2   !< Value to compare.
    logical,          intent(in), optional :: fail !< Fail flag.

    this% n_tests = this% n_tests + 1

    if ( abs(r1-r2) > this% tolerance64 ) then
       if (.not. present(fail) .or. (present(fail) .and. fail .eqv. .false.)) then
          this% n_errors = this% n_errors + 1
       end if
    end if

  end subroutine assert_close_r64

  !> Check if two real (32 bits) arrays (rank 1) are close with respect a tolerance.
  subroutine assert_close_r32_1(this, r1, r2, fail)
    class(tester_t), intent(inout)            :: this !< The tester.
    real(real32),    intent(in), dimension(:) :: r1   !< Value to compare.
    real(real32),    intent(in), dimension(:) :: r2   !< Value to compare.
    logical,         intent(in), optional     :: fail !< Fail flag.

    this% n_tests = this% n_tests + 1

    if ( size(r1) .ne. size(r2) ) then
       if (.not. present(fail) .or. (present(fail) .and. fail .eqv. .false.)) then
          this% n_errors = this% n_errors + 1
       end if
    else
       if ( maxval(abs(r1-r2)) > this% tolerance64 ) then
          if (.not. present(fail) .or. (present(fail) .and. fail .eqv. .false.)) then
             this% n_errors = this% n_errors + 1
          end if
       end if
    end if

  end subroutine assert_close_r32_1

  !> Check if two real (64 bits) arrays (rank 1) are close with respect a tolerance.
  subroutine assert_close_r64_1(this, r1, r2, fail)
    class(tester_t), intent(inout)            :: this !< The tester.
    real(real64),    intent(in), dimension(:) :: r1   !< Value to compare.
    real(real64),    intent(in), dimension(:) :: r2   !< Value to compare.
    logical,         intent(in), optional     :: fail !< Fail flag.

    this% n_tests = this% n_tests + 1

    if ( size(r1) .ne. size(r2) ) then
       if (.not. present(fail) .or. (present(fail) .and. fail .eqv. .false.)) then
          this% n_errors = this% n_errors + 1
       end if
    else
       if ( maxval(abs(r1-r2)) > this% tolerance64 ) then
          if (.not. present(fail) .or. (present(fail) .and. fail .eqv. .false.)) then
             this% n_errors = this% n_errors + 1
          end if
       end if
    end if

  end subroutine assert_close_r64_1
  
  !> Check if two complex numbers (32 bits) are close with respect a tolerance.
  subroutine assert_close_c32(this, c1, c2, fail)
    class(tester_t), intent(inout)        :: this !< The tester.
    complex(real32), intent(in)           :: c1   !< Value to compare.
    complex(real32), intent(in)           :: c2   !< Value to compare.
    logical,         intent(in), optional :: fail !< Fail flag.

    this% n_tests = this% n_tests + 1

    if ( abs(c1-c2) > this% tolerance32 ) then
       if (.not. present(fail) .or. (present(fail) .and. fail .eqv. .false.)) then
          this% n_errors = this% n_errors + 1
       end if
    end if

  end subroutine assert_close_c32

  !> Check if two complex numbers (64 bits) are close with respect a tolerance.
  subroutine assert_close_c64(this, r1, c2, fail)
    class(tester_t),  intent(inout)        :: this !< The tester.
    complex(real64),  intent(in)           :: r1   !< Value to compare.
    complex(real64),  intent(in)           :: c2   !< Value to compare.
    logical,          intent(in), optional :: fail !< Fail flag.

    this% n_tests = this% n_tests + 1

    if ( abs(r1-c2) > this% tolerance64 ) then
       if (.not. present(fail) .or. (present(fail) .and. fail .eqv. .false.)) then
          this% n_errors = this% n_errors + 1
       end if
    end if

  end subroutine assert_close_c64

  !> Check if two complex (32 bits) arrays (rank 1) are close with respect a tolerance.
  subroutine assert_close_c32_1(this, c1, c2, fail)
    class(tester_t), intent(inout)            :: this !< The tester.
    complex(real32), intent(in), dimension(:) :: c1   !< Value to compare.
    complex(real32), intent(in), dimension(:) :: c2   !< Value to compare.
    logical,         intent(in), optional     :: fail !< Fail flag.

    this% n_tests = this% n_tests + 1

    if ( size(c1) .ne. size(c2) ) then
       if (.not. present(fail) .or. (present(fail) .and. fail .eqv. .false.)) then
          this% n_errors = this% n_errors + 1
       end if
    else
       if ( maxval(abs(c1-c2)) > this% tolerance32 ) then
          if (.not. present(fail) .or. (present(fail) .and. fail .eqv. .false.)) then
             this% n_errors = this% n_errors + 1
          end if
       end if
    end if

  end subroutine assert_close_c32_1

  !> Check if two real (64 bits) arrays (rank 1) are close with respect a tolerance.
  subroutine assert_close_c64_1(this, c1, c2, fail)
    class(tester_t), intent(inout)            :: this !< The tester.
    complex(real64), intent(in), dimension(:) :: c1   !< Value to compare.
    complex(real64), intent(in), dimension(:) :: c2   !< Value to compare.
    logical,         intent(in), optional     :: fail !< Fail flag.

    this% n_tests = this% n_tests + 1

    if ( size(c1) .ne. size(c2) ) then
       if (.not. present(fail) .or. (present(fail) .and. fail .eqv. .false.)) then
          this% n_errors = this% n_errors + 1
       end if
    else
       if ( maxval(abs(c1-c2)) > this% tolerance64 ) then
          if (.not. present(fail) .or. (present(fail) .and. fail .eqv. .false.)) then
             this% n_errors = this% n_errors + 1
          end if
       end if
    end if

  end subroutine assert_close_c64_1

end module tester
