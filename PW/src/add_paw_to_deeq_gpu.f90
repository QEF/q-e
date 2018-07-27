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
  INTEGER, ALLOCATABLE :: paw_na_h(:), paw_na_d(:)
#if defined(__CUDA)
  attributes(DEVICE) :: deeq_d, ddd_paw_d, paw_na_d
#endif

! OPTIMIZE HERE: squeeze loop on atoms having PAW pseudo
! OPTIMIZE HERE: use buffers

  if (okpaw) then
     ALLOCATE(ddd_paw_d, SOURCE=ddd_paw)
     ALLOCATE(paw_na_h(nat))
     paw_na_h = -1
     nab=0
     do na=1,nat
        nt = ityp(na)
        IF (upf(nt)%tpawp) nab = nab + 1
        IF (upf(nt)%tpawp) paw_na_h(nab) = na
     end do
     ALLOCATE(paw_na_d, source=paw_na_h)
     
!$cuf kernel do(4)
     do is=1,nspin
        do nb=1,nab
           do ih=1,nhm
              do jh=1,nhm
                 if (jh >= ih) then
                    na = paw_na_d(nb)
                    ijh = jh + ((ih-1)*(2*nhm-ih))/2
                    deeq_d(ih,jh,na,is) = deeq_d(ih,jh,na,is) &
                                           + ddd_paw_d(ijh,na,is)
                    if (jh > ih) deeq_d(jh,ih,na,is) = deeq_d(ih,jh,na,is)
                 end if
              end do
           end do
        end do
     end do
     DEALLOCATE(ddd_paw_d, paw_na_d, paw_na_h)
  end IF 
  
  RETURN
  
END SUBROUTINE add_paw_to_deeq_gpu
