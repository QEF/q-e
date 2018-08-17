!
! Copyright (C) Quantum ESPRESSO group
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
MODULE gbuffers
  USE tb_dev, ONLY : tb_dev_t
  USE tb_pin, ONLY : tb_pin_t
  !
  IMPLICIT NONE
  !
  TYPE(tb_dev_t) :: dev_buf
  TYPE(tb_pin_t) :: pin_buf
END MODULE gbuffers
