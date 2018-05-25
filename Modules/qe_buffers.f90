!
! Copyright (C) Quantum ESPRESSO group
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
MODULE qe_buffers
  USE fbuf_dev, ONLY : fbuf_dev_t
  USE fbuf_pin, ONLY : fbuf_pin_t
  !
  IMPLICIT NONE
  !
  TYPE(fbuf_dev_t) :: qe_buffer
END MODULE qe_buffers
