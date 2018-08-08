!
! Copyright (C) 2010 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!

SUBROUTINE add_paw_to_deeq_gpu(deeq_d)
     ! Add paw contributions to deeq (computed in paw_potential)
  USE kinds,                ONLY : DP
  USE ions_base,     ONLY : nat, ntyp => nsp, ityp
  USE uspp_param,    ONLY : upf, nh, nhm
  USE paw_variables, ONLY : okpaw, ddd_paw
  USE lsda_mod,      ONLY : nspin
  IMPLICIT NONE
  integer :: na, nb, nab, nt, ih, jh, ijh, nhnt, is
  REAL(kind=dp), intent(inout) :: deeq_d( nhm, nhm, nat, nspin )
  REAL(DP), ALLOCATABLE :: ddd_paw_d(:,:,:)
#if defined(__CUDA)
  attributes(DEVICE) :: deeq_d, ddd_paw_d
#endif

! OPTIMIZE HERE: squeeze loop on atoms having PAW pseudo
! OPTIMIZE HERE: use buffers

  if (okpaw) then
     ALLOCATE(ddd_paw_d, SOURCE=ddd_paw)
     do na=1,nat
        nt = ityp(na)
        IF (.not.upf(nt)%tpawp) cycle
        nhnt = nh(nt)
!$cuf kernel do(3)
        do is=1,nspin
           do ih=1,nhnt
              do jh=1,nhnt
                 if (jh >= ih) then
                    ijh = jh + ((ih-1)*(2*nhnt-ih))/2
                    deeq_d(ih,jh,na,is) = deeq_d(ih,jh,na,is) &
                                           + ddd_paw_d(ijh,na,is)
                    deeq_d(jh,ih,na,is) = deeq_d(ih,jh,na,is) 
                 end if
              end do
           end do
        end do
     end do
     deallocate(ddd_paw_d)
  end IF 
  !
  RETURN
  
END SUBROUTINE add_paw_to_deeq_gpu
