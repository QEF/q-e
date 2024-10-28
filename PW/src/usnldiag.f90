!
! Copyright (C) 2001-2024 QUantum ESPRESSO Foundation
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
SUBROUTINE usnldiag (npw, h_diag, s_diag)
  !-----------------------------------------------------------------------
  !
  !    add nonlocal pseudopotential term to diagonal part of Hamiltonian,
  !    compute the diagonal part of the S matrix
  !
  USE kinds,            ONLY: DP
  USE ions_base,        ONLY: nat, ityp, ntyp => nsp
  USE wvfct,            ONLY: npwx
  USE lsda_mod,         ONLY: current_spin
  USE uspp,             ONLY: deeq, vkb, qq_at, qq_so, deeq_nc, ofsbeta
  USE uspp_param,       ONLY: upf, nh, nhm
  USE noncollin_module, ONLY: noncolin, npol, lspinorb
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(in) :: npw
  ! number of plane waves
  REAL(dp), INTENT(inout) :: h_diag (npwx,npol)
  ! the diagonal part of the hamiltonian
  REAL(dp), INTENT(out)   :: s_diag (npwx,npol)
  ! the diagonal part of the S matrix
  !
  INTEGER :: ikb, jkb, ih, jh, na, nt, ig, ipol, nhnt, offset
  REAL(DP), ALLOCATABLE :: pr1(:,:), pr2(:,:)
  COMPLEX(DP), ALLOCATABLE :: pc1(:,:,:), pc2(:,:,:)
  !$acc declare device_resident(pr1, pr2, pc1, pc2)
  REAL(DP)    :: ar
  COMPLEX(DP) :: ac
  !
  !$acc data present ( h_diag, s_diag, vkb, deeq, deeq_nc, qq_so, qq_at ) 
  !
  ! initialise s_diag
  !
  !$acc kernels
  s_diag(:,:) = 1.d0
  !$acc end kernels
  !
  IF (noncolin) THEN
  !
  ALLOCATE( pc1(nhm,nhm,npol), pc2(nhm,nhm,npol) )
  DO nt = 1, ntyp
     nhnt = nh(nt)
     DO na = 1, nat
        IF (ityp (na) == nt) THEN
           offset = ofsbeta(na)
           !$acc parallel loop collapse(2)
           DO ih = 1, nhnt
              DO jh = 1, nhnt
                 ! Note: in order to keep the code simple, here we deal
                 !       with the generale case of nondiagonal "deeq" only,
                 !       but for non-US and non-multiprojector case, the
                 !       diagonal ih = jh terms are the only nonzero terms
                 pc1(ih,jh,1) = deeq_nc (ih, jh, na, 1)
                 pc1(ih,jh,2) = deeq_nc (ih, jh, na, 4)
                 IF (lspinorb) THEN
                    pc2(ih,jh,1) = qq_so(ih, jh, 1, nt)
                    pc2(ih,jh,2) = qq_so(ih, jh, 4, nt)
                 ELSE
                    pc2(ih,jh,1) = qq_at (ih, jh, na)
                    pc2(ih,jh,2) = qq_at (ih, jh, na)
                 ENDIF
              ENDDO
           ENDDO
           !$acc parallel loop collapse(2)
           DO ig = 1, npw
              DO ipol = 1, npol
                 !$acc loop seq collapse(2)
                 DO ih = 1, nhnt
                    DO jh = 1, nhnt
                       ikb = offset + ih
                       jkb = offset + jh
                       ac = vkb (ig, ikb) *conjg( vkb (ig, jkb))
                       h_diag (ig,ipol) = h_diag (ig,ipol) + pc1(ih,jh,ipol) * ac
                       s_diag (ig,ipol) = s_diag (ig,ipol) + pc2(ih,jh,ipol) * ac
                    ENDDO
                 ENDDO
              ENDDO
           ENDDO
        ENDIF
     ENDDO
  ENDDO
  DEALLOCATE( pc2, pc1 )
  !
  ELSE
  ! 
  ALLOCATE( pr1(nhm,nhm), pr2(nhm,nhm) )
  DO nt = 1, ntyp
     nhnt = nh(nt)
     DO na = 1, nat
        IF (ityp (na) == nt) THEN
           offset = ofsbeta(na)
           !$acc parallel loop collapse(2)
           DO ih = 1, nhnt
              DO jh = 1, nhnt
                 ! Note: in order to keep the code simple, here we deal
                 !       with the generale case of nondiagonal "deeq" only,
                 !       but for non-US and non-multiprojector case, the
                 !       diagonal ih = jh terms are the only nonzero terms
                 pr1(ih,jh) = deeq (ih, jh, na, current_spin)
                 pr2(ih,jh) = qq_at (ih, jh, na)
              ENDDO
           ENDDO
           !$acc parallel loop 
           DO ig = 1, npw
              !$acc loop seq collapse(2)
              DO ih = 1, nhnt
                 DO jh = 1, nhnt
                    ikb = offset + ih
                    jkb = offset + jh
                    ! The diagonal is real
                    ar = vkb (ig, ikb) *conjg( vkb (ig, jkb))
                    h_diag (ig,1) = h_diag (ig,1) + pr1(ih,jh) * ar
                    s_diag (ig,1) = s_diag (ig,1) + pr2(ih,jh) * ar
                 ENDDO
              ENDDO
           ENDDO
        ENDIF
     ENDDO
  ENDDO
  DEALLOCATE( pr2, pr1 )
  !
  ENDIF
  !$acc end data

  RETURN
END SUBROUTINE usnldiag
