!
! Copyright (C) 2010 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
SUBROUTINE add_paw_to_deeq_gpu(deeq_d)
  !
  !! Add paw contributions to the integral of the perturbed potential 
  !! with the Q function (computed in paw_potential).
  !
  USE kinds,         ONLY : DP
  USE ions_base,     ONLY : nat, ntyp => nsp, ityp
  USE uspp_param,    ONLY : upf, nh, nhm
  USE paw_variables, ONLY : okpaw, ddd_paw
  USE lsda_mod,      ONLY : nspin
  !
  IMPLICIT NONE
  !
  REAL(DP), INTENT(INOUT) :: deeq_d(nhm,nhm,nat,nspin)
  !! integral of the perturbed potential
  !
  ! ... local variables
  !
  INTEGER :: na, nb, nab, nt, ih, jh, ijh, nhnt, is
  REAL(DP), ALLOCATABLE :: ddd_paw_d(:,:,:)
#if defined(__CUDA)
  attributes(DEVICE) :: ddd_paw_d, deeq_d
#endif

! OPTIMIZE HERE: squeeze loop on atoms having PAW pseudo
! OPTIMIZE HERE: use buffers

  IF (okpaw) THEN
     ALLOCATE(ddd_paw_d, SOURCE=ddd_paw)
     DO na=1,nat
        nt = ityp(na)
        IF (.not.upf(nt)%tpawp) CYCLE
        nhnt = nh(nt)
        !$acc parallel loop collapse(3)
        DO is=1,nspin
           DO ih=1,nhnt
              DO jh=1,nhnt
                 IF (jh >= ih) then
                    ijh = jh + ((ih-1)*(2*nhnt-ih))/2
                    deeq_d(ih,jh,na,is) = deeq_d(ih,jh,na,is) &
                                           + ddd_paw_d(ijh,na,is)
                    deeq_d(jh,ih,na,is) = deeq_d(ih,jh,na,is) 
                 ENDIF
              ENDDO
           ENDDO
        ENDDO
     ENDDO
     DEALLOCATE(ddd_paw_d)
  ENDIF
  !
  RETURN
  !
END SUBROUTINE add_paw_to_deeq_gpu
