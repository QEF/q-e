!
! Copyright (C) Quantum ESPRESSO group
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
MODULE data_buffer
    USE util_param,  ONLY : DP
#ifdef __CUDA
    USE cudafor
#endif
    !
    IMPLICIT NONE
    !
    REAL(DP), ALLOCATABLE, dimension(:)  :: mp_buff_r
    INTEGER, ALLOCATABLE, dimension(:)  :: mp_buff_i
    PUBLIC :: mp_buff_r, mp_buff_i
    !
#ifdef __CUDA
    REAL(DP), ALLOCATABLE, dimension(:)  :: mp_buff_r_d
    INTEGER, ALLOCATABLE, dimension(:)  :: mp_buff_i_d
    ATTRIBUTES( DEVICE ) :: mp_buff_r_d, mp_buff_i_d
    PUBLIC :: mp_buff_r_d, mp_buff_i_d
#endif

END MODULE data_buffer
