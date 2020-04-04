
! Copyright (C) 2002-2017 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! This module contains interfaces to call the c-routine clib/fletcher32.c
! implementing the Fletcher-32 checksum algorithm as reported in
! https://en.wikipedia.org/wiki/Fletcher%27s_checksum#Optimizations
!
! SdG September 3rd 2017
!
!------------------------------------------------------------------------------!
    MODULE fletcher32_mod
!------------------------------------------------------------------------------!
    USE util_param,     ONLY : DP
    !
    IMPLICIT NONE
    PRIVATE
    integer(2) :: dat(1)

    PUBLIC :: fletcher32_cksum, fletcher32
!
    INTERFACE fletcher32_cksum
       MODULE PROCEDURE fletcher32_i1, fletcher32_r1, fletcher32_c1, fletcher32_z,  fletcher32_l,  &
                        fletcher32_iv, fletcher32_rv, fletcher32_cv, fletcher32_zv, fletcher32_lv, &
                        fletcher32_im, fletcher32_rm, fletcher32_cm,                fletcher32_lm, &
                        fletcher32_it, fletcher32_rt, fletcher32_ct, &
                        fletcher32_i4, fletcher32_r4, fletcher32_c4, &
                                       fletcher32_r5, fletcher32_c5
    END INTERFACE

    INTERFACE
       FUNCTION fletcher32( dat, dat_size ) BIND(C,name="fletcher32") RESULT(t)
          USE ISO_C_BINDING
          integer(kind=c_int16_t) :: dat(*)
          integer(kind=c_int32_t) :: dat_size
          integer(kind=c_int32_t) :: t
       END FUNCTION fletcher32
    END INTERFACE
!
!------------------------------------------------------------------------------!
!
    CONTAINS
!
!------------------------------------------------------------------------------!

!..fletcher32_cksum
!------------------------------------------------------------------------------!
      SUBROUTINE fletcher32_i1(msg, cksum)
         IMPLICIT NONE
         INTEGER, INTENT(IN) :: msg
         INTEGER, INTENT(OUT) :: cksum

         cksum = fletcher32(transfer(msg,dat),size(transfer(msg,dat)))

      END SUBROUTINE fletcher32_i1
!
!------------------------------------------------------------------------------!
      SUBROUTINE fletcher32_iv(msg, cksum)
         IMPLICIT NONE
         INTEGER, INTENT(IN) :: msg(:)
         INTEGER, INTENT(OUT) :: cksum

         cksum = fletcher32(transfer(msg,dat),size(transfer(msg,dat)))

      END SUBROUTINE fletcher32_iv
!
!------------------------------------------------------------------------------!
      SUBROUTINE fletcher32_im( msg, cksum )
         IMPLICIT NONE
         INTEGER, INTENT(IN) :: msg(:,:)
         INTEGER, INTENT(OUT) :: cksum

         cksum = fletcher32(transfer(msg,dat),size(transfer(msg,dat)))

      END SUBROUTINE fletcher32_im
!
!------------------------------------------------------------------------------!
      SUBROUTINE fletcher32_it( msg, cksum )
         IMPLICIT NONE
         INTEGER, INTENT(IN) :: msg(:,:,:)
         INTEGER, INTENT(OUT) :: cksum

         cksum = fletcher32(transfer(msg,dat),size(transfer(msg,dat)))

      END SUBROUTINE fletcher32_it
!
!------------------------------------------------------------------------------!
      SUBROUTINE fletcher32_i4(msg, cksum )
         IMPLICIT NONE
         INTEGER, INTENT(IN) :: msg(:,:,:,:)
         INTEGER, INTENT(OUT) :: cksum

         cksum = fletcher32(transfer(msg,dat),size(transfer(msg,dat)))

      END SUBROUTINE fletcher32_i4
!
!------------------------------------------------------------------------------!
      SUBROUTINE fletcher32_r1( msg, cksum  )
         IMPLICIT NONE
         REAL(DP), INTENT(IN) :: msg
         INTEGER, INTENT(OUT) :: cksum

         cksum = fletcher32(transfer(msg,dat),size(transfer(msg,dat)))

      END SUBROUTINE fletcher32_r1
!
!------------------------------------------------------------------------------!
      SUBROUTINE fletcher32_rv(msg, cksum )
         IMPLICIT NONE
         REAL(DP), INTENT(IN) :: msg(:)
         INTEGER, INTENT(OUT) :: cksum

         cksum = fletcher32(transfer(msg,dat),size(transfer(msg,dat)))

      END SUBROUTINE fletcher32_rv
!
!------------------------------------------------------------------------------!
      SUBROUTINE fletcher32_rm(msg, cksum )
         IMPLICIT NONE
         REAL(DP), INTENT(IN) :: msg(:,:)
         INTEGER, INTENT(OUT) :: cksum

         cksum = fletcher32(transfer(msg,dat),size(transfer(msg,dat)))

      END SUBROUTINE fletcher32_rm
!
!------------------------------------------------------------------------------!
      SUBROUTINE fletcher32_rt(msg, cksum )
         IMPLICIT NONE
         REAL(DP), INTENT(IN) :: msg(:,:,:)
         INTEGER, INTENT(OUT) :: cksum

         cksum = fletcher32(transfer(msg,dat),size(transfer(msg,dat)))

      END SUBROUTINE fletcher32_rt
!
!------------------------------------------------------------------------------!
      SUBROUTINE fletcher32_r4(msg, cksum )
         IMPLICIT NONE
         REAL(DP), INTENT(IN) :: msg(:,:,:,:)
         INTEGER, INTENT(OUT) :: cksum

         cksum = fletcher32(transfer(msg,dat),size(transfer(msg,dat)))

      END SUBROUTINE fletcher32_r4
!
!------------------------------------------------------------------------------!
      SUBROUTINE fletcher32_r5(msg, cksum )
         IMPLICIT NONE
         REAL(DP), INTENT(IN) :: msg(:,:,:,:,:)
         INTEGER, INTENT(OUT) :: cksum

         cksum = fletcher32(transfer(msg,dat),size(transfer(msg,dat)))

      END SUBROUTINE fletcher32_r5
!
!------------------------------------------------------------------------------!
      SUBROUTINE fletcher32_c1(msg, cksum )
         IMPLICIT NONE
         COMPLEX(DP), INTENT(IN) :: msg
         INTEGER, INTENT(OUT) :: cksum

         cksum = fletcher32(transfer(msg,dat),size(transfer(msg,dat)))

      END SUBROUTINE fletcher32_c1
!
!------------------------------------------------------------------------------!
      SUBROUTINE fletcher32_cv(msg, cksum )
         IMPLICIT NONE
         COMPLEX(DP), INTENT(IN) :: msg(:)
         INTEGER, INTENT(OUT) :: cksum

         cksum = fletcher32(transfer(msg,dat),size(transfer(msg,dat)))

      END SUBROUTINE fletcher32_cv
!
!------------------------------------------------------------------------------!
      SUBROUTINE fletcher32_cm(msg, cksum )
         IMPLICIT NONE
         COMPLEX(DP), INTENT(IN) :: msg(:,:)
         INTEGER, INTENT(OUT) :: cksum

         cksum = fletcher32(transfer(msg,dat),size(transfer(msg,dat)))

      END SUBROUTINE fletcher32_cm
!
!------------------------------------------------------------------------------!
      SUBROUTINE fletcher32_ct(msg, cksum )
         IMPLICIT NONE
         COMPLEX(DP), INTENT(IN) :: msg(:,:,:)
         INTEGER, INTENT(OUT) :: cksum

         cksum = fletcher32(transfer(msg,dat),size(transfer(msg,dat)))

      END SUBROUTINE fletcher32_ct
!
!------------------------------------------------------------------------------!
      SUBROUTINE fletcher32_c4(msg, cksum )
         IMPLICIT NONE
         COMPLEX(DP), INTENT(IN) :: msg(:,:,:,:)
         INTEGER, INTENT(OUT) :: cksum

         cksum = fletcher32(transfer(msg,dat),size(transfer(msg,dat)))

      END SUBROUTINE fletcher32_c4
!
!------------------------------------------------------------------------------!
      SUBROUTINE fletcher32_c5(msg, cksum )
         IMPLICIT NONE
         COMPLEX(DP), INTENT(IN) :: msg(:,:,:,:,:)
         INTEGER, INTENT(OUT) :: cksum

         cksum = fletcher32(transfer(msg,dat),size(transfer(msg,dat)))

      END SUBROUTINE fletcher32_c5
!
!------------------------------------------------------------------------------!
      SUBROUTINE fletcher32_l(msg, cksum )
         IMPLICIT NONE
         LOGICAL, INTENT(IN) :: msg
         INTEGER, INTENT(OUT) :: cksum

         cksum = fletcher32(transfer(msg,dat),size(transfer(msg,dat)))

      END SUBROUTINE fletcher32_l
!
!------------------------------------------------------------------------------!
      SUBROUTINE fletcher32_lv(msg, cksum )
         IMPLICIT NONE
         LOGICAL, INTENT(IN) :: msg(:)
         INTEGER, INTENT(OUT) :: cksum

         cksum = fletcher32(transfer(msg,dat),size(transfer(msg,dat)))

      END SUBROUTINE fletcher32_lv
!
!------------------------------------------------------------------------------!
      SUBROUTINE fletcher32_lm(msg, cksum )
         IMPLICIT NONE
         LOGICAL, INTENT(IN) :: msg(:,:)
         INTEGER, INTENT(OUT) :: cksum

         cksum = fletcher32(transfer(msg,dat),size(transfer(msg,dat)))

      END SUBROUTINE fletcher32_lm
!
!------------------------------------------------------------------------------!
      SUBROUTINE fletcher32_z(msg, cksum )
         IMPLICIT NONE
         CHARACTER(len=*), INTENT(IN) :: msg
         INTEGER, INTENT(OUT) :: cksum

         cksum = fletcher32(transfer(msg,dat),size(transfer(msg,dat)))

      END SUBROUTINE fletcher32_z
!
!------------------------------------------------------------------------------!
      SUBROUTINE fletcher32_zv(msg, cksum )
         IMPLICIT NONE
         CHARACTER(len=*), INTENT(IN) :: msg(:)
         INTEGER, INTENT(OUT) :: cksum

         cksum = fletcher32(transfer(msg,dat),size(transfer(msg,dat)))

      END SUBROUTINE fletcher32_zv
!
!------------------------------------------------------------------------------!
    END MODULE fletcher32_mod
!------------------------------------------------------------------------------!
