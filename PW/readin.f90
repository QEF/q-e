!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
subroutine readpp
  !-----------------------------------------------------------------------
  !
  !    Read pseudopotentials
  !
#include "machine.h"
  use pwcom
  use io, only: pseudo_dir, pseudop
  !
  character :: file_pseudo * 256
  ! file name complete with path
  integer :: iunps, isupf, l, nt, ios, pseudo_type
  external pseudo_type
  !
  iunps = 4
  l = len_trim (pseudo_dir)
  do nt = 1, ntyp
     ! add / if needed
     if (pseudo_dir (l:l) .ne.'/') then
        file_pseudo = pseudo_dir (1:l) //'/'//pseudop (nt)
     else
        file_pseudo = pseudo_dir (1:l) //pseudop (nt)
     endif
     !
     open (unit = iunps, file = file_pseudo, status = 'old', form = &
          'formatted', iostat = ios)
     call errore ('readin', 'file '//trim(file_pseudo)//' not found', ios)
     call read_pseudo (nt, iunps, isupf)
     !   The new pseudopotential UPF format is detected via the presence
     !   of the keyword '<PP_HEADER>' at the beginning of the file
     if (isupf /= 0) then
        rewind (unit = iunps)
        !
        !     The type of the pseudopotential is determined by the file name:
        !    *.vdb or *.van  Vanderbilt US pseudopotential code  pseudo_type=1
        !    *.RRKJ3         Andrea's   US new code              pseudo_type=2
        !    none of the above: PWSCF norm-conserving format     pseudo_type=0
        !
        if (pseudo_type (pseudop (nt) ) .eq.1 &
             .or.pseudo_type (pseudop (nt) ) .eq.2) then
           !
           !    The vanderbilt pseudopotential is always in numeric form
           !
           numeric (nt) = .true.
           open (unit = iunps, file = file_pseudo, status = 'old', &
                form = 'formatted', iostat = ios)

           !
           !    newpseudo distinguishes beteween US pseudopotentials
           !    produced by Vanderbilt code and those produced
           !    by Andrea's atomic code.
           !
           if (pseudo_type (pseudop (nt) ) .eq.1) then
              newpseudo (nt) = .false.
              tvanp (nt) = .true.
              call readvan (nt, iunps)
           endif
           if (pseudo_type (pseudop (nt) ) .eq.2) then
              newpseudo (nt) = .true.
              ! tvanp is read inside readnewvan
              call readnewvan (nt, iunps)
           endif
           close (iunps)
        else
           tvanp (nt) = .false.
           newpseudo (nt) = .false.
           open (unit = iunps, file = file_pseudo, status = 'old', &
                form='formatted', iostat = ios)
           ! numeric is read inside read_ncpp
           call read_ncpp (nt, iunps)
           close (iunps)
        endif
     else
        ! UPF is always numeric
        numeric (nt) = .true.
        ! UPF is RRKJ3-like
        newpseudo (nt) = .true.
        close (iunps)
     endif

  enddo
  return
end subroutine readpp
!-----------------------------------------------------------------------
integer function pseudo_type (pseudop)
  !-----------------------------------------------------------------------
  implicit none
  character (len=*) :: pseudop
  integer :: l
  !
  l = len_trim (pseudop)
  pseudo_type = 0
  if (pseudop (l - 3:l) .eq.'.vdb'.or.pseudop (l - 3:l) .eq.'.van') &
       pseudo_type = 1
  if (l > 5) then
     if (pseudop (l - 5:l) .eq.'.RRKJ3') pseudo_type = 2
  end if
  !
  return

end function pseudo_type

