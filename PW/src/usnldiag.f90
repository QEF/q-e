!
! Copyright (C) 2001 PWSCF group
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
  !    add nonlocal pseudopotential term to diagonal part of Hamiltonian
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
  COMPLEX(DP), ALLOCATABLE :: ps1(:,:,:), ps2(:,:,:)
  !$acc declare device_resident(ps1, ps2)
  COMPLEX(DP) :: ar
  !
  !$acc data present ( h_diag, s_diag, vkb, deeq, deeq_nc, qq_so, qq_at ) 
  !
  ! initialise s_diag
  !
  !$acc kernels
  s_diag(:,:) = 1.d0
  !$acc end kernels
  !
  !    multiply on projectors
  !
  ALLOCATE( ps1(nhm,nhm,npol), ps2(nhm,nhm,npol) )
  DO nt = 1, ntyp
     nhnt = nh(nt)
     DO na = 1, nat
        IF (ityp (na) == nt) THEN
           offset = ofsbeta(na)
           !$acc parallel loop collapse(2)
           DO ih = 1, nhnt
              DO jh = 1, nhnt
                 IF (lspinorb) THEN
                    ps1(ih,jh,1) = deeq_nc (ih, jh, na, 1)
                    ps1(ih,jh,2) = deeq_nc (ih, jh, na, 4)
                    ps2(ih,jh,1) = qq_so(ih, jh, 1, nt)
                    ps2(ih,jh,2) = qq_so(ih, jh, 4, nt)
                 ELSEIF (noncolin) THEN
                    ps1(ih,jh,1) = deeq_nc (ih, jh, na, 1)
                    ps1(ih,jh,2) = deeq_nc (ih, jh, na, 4)
                    ps2(ih,jh,1) = qq_at (ih, jh, na)
                    ps2(ih,jh,2) = qq_at (ih, jh, na)
                 ELSE
                    ps1(ih,jh,1) = deeq (ih, jh, na, current_spin)
                    ps2(ih,jh,1) = qq_at (ih, jh, na)
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
                       ar = vkb (ig, ikb) *conjg( vkb (ig, jkb))
                       h_diag (ig,ipol) = h_diag (ig,ipol) + ps1(ih,jh,ipol) * ar
                       s_diag (ig,ipol) = s_diag (ig,ipol) + ps2(ih,jh,ipol) * ar
                    ENDDO
                 ENDDO
              ENDDO
           ENDDO
        ENDIF
     ENDDO
  ENDDO
  DEALLOCATE( ps2, ps1 )
  !$acc end data

  RETURN
END SUBROUTINE usnldiag
