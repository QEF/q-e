!
! Copyright (C) 2010-2013 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
SUBROUTINE add_vhub_to_deeq_gpu(deeq_d)
  !
  ! Add Hubbard contributions to deeq when U_projection is pseudo
  !
  USE kinds,         ONLY : DP
  USE ions_base,     ONLY : nat, ntyp => nsp, ityp
  USE uspp_param,    ONLY : nh, nhm
  USE lsda_mod,      ONLY : nspin
  USE scf,           ONLY : v
  USE ldaU,          ONLY : is_hubbard, Hubbard_l, offsetU, q_ae
  IMPLICIT NONE
  REAL(KIND=DP), INTENT(INOUT) :: deeq_d( nhm, nhm, nat, nspin )
#if defined(__CUDA)
  attributes(DEVICE) :: deeq_d
#endif
  ! Local variables
  REAL(KIND=DP), ALLOCATABLE :: deeq_aux_h( :, :, : )
  REAL(KIND=DP), ALLOCATABLE :: deeq_aux_d( :, :, : )
#if defined(__CUDA)
  attributes(DEVICE) :: deeq_aux_d
#endif
  INTEGER :: na, nt, ih, jh, ijh, m1, m2, ow1, ow2, is, nhnt
  !
  ! (maybe) OPTIMIZE here ... reorder the loop ?
  !
  ALLOCATE(deeq_aux_h( nhm, nhm, nspin ))
  ALLOCATE(deeq_aux_d( nhm, nhm, nspin ))
  DO na = 1, nat
     !
     nt = ityp(na)
     !
     ! skip atoms without Hubbard U
     IF ( .NOT. is_hubbard(nt) ) CYCLE
     !
     deeq_aux_h = 0.d0
     DO ih = 1, nh(nt)
        DO jh = ih, nh(nt)
           !
           DO m1 = 1, 2 * Hubbard_l(nt) + 1
              DO m2 = 1, 2 * Hubbard_l(nt) + 1
                 !
                 ow1 = offsetU(na)+m1
                 ow2 = offsetU(na)+m2
                 deeq_aux_h(ih,jh,1:nspin) = deeq_aux_h(ih,jh,1:nspin) + &
                    v%ns(m1,m2,1:nspin,na)*q_ae(ow1,ih,na)*q_ae(ow2,jh,na)
                 !
              ENDDO 
           ENDDO
           !
           deeq_aux_h(jh,ih,1:nspin) = deeq_aux_h(ih,jh,1:nspin)
           !
        ENDDO
     ENDDO
     deeq_aux_d = deeq_aux_h
     nhnt = nh(nt)
     !$cuf kernel do(3)
     DO is=1, nspin
        DO ih = 1, nhnt
           DO jh = 1, nhnt
              deeq_d(jh,ih,na,is) = deeq_d(jh,ih,na,is) + deeq_aux_d(jh,ih,is)
           END DO
        END DO
     END DO
     !
  ENDDO
  
  DEALLOCATE(deeq_aux_h, deeq_aux_d)
  !
END SUBROUTINE add_vhub_to_deeq_gpu
