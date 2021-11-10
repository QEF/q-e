module debug_utils

use kinds, only: dp, i8b, i4b
USE io_global, only : ionode
USE mp, ONLY: mp_max, mp_sum
implicit none
save
private
public :: init_debug_utils, checkpoint

integer :: counter = 0 !incremented at each call
character(len=64) :: fname_prefix = 'debug.' !prefix of the files to read/write
logical :: testing = .false. !if true read the file and compare with the provided array
logical :: checkpoint_debug_active = .false. !writing/reading of file is active
integer :: communicator

character(len=64) :: msg_pref = '@DEBUG checkpoint: ' 

include 'debug_proc.fh'

contains

        subroutine init_debug_utils(active, prefix, read_file)
        use mp_world, ONLY: mpime
        logical, intent(in) :: active, read_file
        character(len=64), intent(in) :: prefix
        CHARACTER(LEN=6), EXTERNAL :: int_to_char

        counter = 0
        fname_prefix = trim(prefix)//trim(int_to_char(mpime)) // '.'
        checkpoint_debug_active = active
        testing = read_file

        if (ionode .and.  active ) then
           write (*,*) '++++++++++ DEBUG MODE ON +++++++++++'
           write (*,*) 'file prefix = ', prefix
           write (*,*) 'read from file and test = ', read_file
           write (*,*) '++++++++++ ============= +++++++++++'
        end if

        end subroutine

        include 'debug_sub.fh'



end module
