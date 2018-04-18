!
! Copyright (C) 2002-2011 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

!=----------------------------------------------------------------------------=!
   MODULE wavefunctions_module_gpum
!=----------------------------------------------------------------------------=!
     USE kinds, ONLY :  DP

     IMPLICIT NONE
     SAVE

     !
     COMPLEX(DP), DEVICE, ALLOCATABLE, TARGET :: &
       evc_d(:,:)     ! wavefunctions in the PW basis set
                    ! noncolinear case: first index
                    ! is a combined PW + spin index
     !
     COMPLEX(DP) , DEVICE, ALLOCATABLE, TARGET :: &
       psic_d(:), &      ! additional memory for FFT
       psic_nc_d(:,:)    ! as above for the noncolinear case
     !
     
  CONTAINS
     !
     ! ==== SUBS FOR EVC
     !
     SUBROUTINE start_sync_evc_dev()
        use wavefunctions_module, only: evc
        implicit none
        if (.not. allocated(evc_d)) then 
           allocate(evc_d, source=evc)
           return
        end if
        evc_d = evc
        !
     END SUBROUTINE start_sync_evc_dev
     !
     SUBROUTINE start_sync_evc_hst()
        use wavefunctions_module, only: evc
        implicit none
        evc = evc_d
     END SUBROUTINE start_sync_evc_hst

     SUBROUTINE wait_sync_evc_dev()
        RETURN
     END SUBROUTINE wait_sync_evc_dev
     SUBROUTINE wait_sync_evc_hst()
        RETURN
     END SUBROUTINE wait_sync_evc_hst
     
     
     
!=----------------------------------------------------------------------------=!
   END MODULE wavefunctions_module_gpum
!=----------------------------------------------------------------------------=!
