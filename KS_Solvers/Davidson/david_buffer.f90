!
! Copyright (C) Quantum ESPRESSO group
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
MODULE david_buffer
#if defined(__CUDA)
  USE fbuf_david, ONLY : buf_t
  !
  IMPLICIT NONE
  !
  TYPE(buf_t) :: buffer
#endif
END MODULE david_buffer
