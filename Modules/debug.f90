module debug_utils

use kinds, only: dp
USE io_global, only : ionode
implicit none
save
private
public :: init_debug_utils, checkpoint

integer :: counter = 0 !incremented at each call
character(len=64) :: fname_prefix = 'debug.' !prefix of the files to read/write
logical :: testing = .false. !if true read the file and compare with the provided array
logical :: checkpoint_debug_active = .false. !writing/reading of file is active

character(len=64) :: msg_pref = '@DEBUG checkpoint: ' 

include 'debug_proc.fh'

contains

        subroutine init_debug_utils(active, prefix, read_file)
        logical, intent(in) :: active, read_file
        character(len=64), intent(in) :: prefix

        counter = 0
        fname_prefix = prefix
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
