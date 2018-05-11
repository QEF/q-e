!
! Copyright (C) Quantum ESPRESSO group
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
MODULE qe_buffers
#if defined(__CUDA)
  USE fbuf_qe, ONLY : bufqe_t
  USE fbuf_qe_cpu, ONLY : buf_pinned_t
  !
  IMPLICIT NONE
  !
  TYPE(bufqe_t) :: qe_buffer
#endif
END MODULE qe_buffers
