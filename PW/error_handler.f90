!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
module error_handler
  implicit none
  private

  public :: init_error, add_name, chop_name, error_mem, warning

  type chain
   character (len=35)   :: routine_name
   type(chain), pointer :: previous_link
  end type chain

  type(chain), pointer :: routine_chain

contains

  subroutine init_error(routine_name)
    implicit none
    character (len=*), intent(in) :: routine_name

    allocate(routine_chain)

    routine_chain%routine_name  =  routine_name
    nullify(routine_chain%previous_link)

    return
  end subroutine init_error

  subroutine add_name(routine_name)
    implicit none
    character (len=*), intent(in) :: routine_name
    type(chain), pointer          :: new_link

    allocate(new_link)
    new_link%routine_name  =  routine_name
    new_link%previous_link => routine_chain
    routine_chain          => new_link

    return
  end subroutine add_name

  subroutine chop_name
    implicit none
    type(chain), pointer :: chopped_chain

    chopped_chain => routine_chain%previous_link
    deallocate(routine_chain)
    routine_chain => chopped_chain

    return
  end subroutine chop_name

  recursive subroutine trace_back(error_code)

    implicit none
    integer :: error_code

    write(unit=*,fmt=*) "   Called by ", routine_chain%routine_name
    if (.not.associated(routine_chain%previous_link)) then
       write(unit=*,fmt=*) &
            " +++ +++ +++ +++ +++ +++ +++ +++ +++ +++ +++ +++ +++ +++ +++"
       write(unit=*,fmt=*) " "
       if( error_code > 0 ) then
          stop
       else
          return
       end if
    end if

    routine_chain => routine_chain%previous_link
    call trace_back(error_code)

  end subroutine trace_back

  subroutine error_mem(message,error_code)
    character (len=*), intent(in) :: message
    integer, intent(in), optional :: error_code
    integer                       :: action_code
    type(chain), pointer          :: save_chain

    if (present(error_code)) then
       action_code = error_code
    else
       action_code = 1
    end if

    if( action_code /= 0 ) then
       write(unit=*,fmt=*) " "
       write(unit=*,fmt=*) &
            " +++ +++ +++ +++ +++ +++ +++ +++ +++ +++ +++ +++ +++ +++ +++"

       if( action_code > 0 ) then
          write(unit=*,fmt=*) "   Fatal error in routine `", &
               trim(routine_chain%routine_name),"': ",message
       else
          write(unit=*,fmt=*) "   Warning from routine `", &
               trim(routine_chain%routine_name),"': ",message
          save_chain => routine_chain
       end if
       write(unit=*,fmt=*) &
            " +++ +++ +++ +++ +++ +++ +++ +++ +++ +++ +++ +++ +++ +++ +++"
       routine_chain => routine_chain%previous_link
       call trace_back(action_code)
       routine_chain => save_chain
    end if

    return
  end subroutine error_mem

  subroutine warning(message)
    character (len=*), intent(in) :: message
    call error_mem(message,-1)
    return
  end subroutine warning

end module error_handler

