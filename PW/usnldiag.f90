!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
subroutine usnldiag (h_diag, s_diag)
  !-----------------------------------------------------------------------
  !
  !    add nonlocal pseudopotential term to diagonal part of Hamiltonian
  !    compute the diagonal part of the S matrix
  !
  USE kinds, ONLY: DP
  USE basis,  ONLY: ntyp, nat, ityp
  USE wvfct, ONLY: npw
  USE lsda_mod, ONLY: current_spin 
  USE uspp,  ONLY: deeq, vkb, qq
  USE uspp_param,  ONLY: nh, tvanp, newpseudo
  implicit none
  !
  !    here the dummy variables
  !
  real(kind=DP) :: h_diag (npw), s_diag (npw)
  ! input/output: the diagonal part of the hamiltonian
  ! output: the diagonal part of the S matrix
  !
  !   and here the local variables
  !
  integer :: ikb, jkb, ih, jh, na, nt, ig, ijkb0
  ! counters
  real(kind=DP) :: ps1, ps2, ar
  !
  ! initialise s_diag
  !
  do ig = 1, npw
     s_diag (ig) = 1.d0
  enddo
  !
  !    multiply on projectors
  !
  ijkb0 = 0
  do nt = 1, ntyp
     do na = 1, nat
        if (ityp (na) .eq.nt) then
           do ih = 1, nh (nt)
              ikb = ijkb0 + ih
              ps1 = deeq (ih, ih, na, current_spin)
              ps2 = qq (ih, ih, nt)
              do ig = 1, npw
                 ar = vkb (ig, ikb)*conjg(vkb (ig, ikb))
                 h_diag (ig) = h_diag (ig) + ps1 * ar
                 s_diag (ig) = s_diag (ig) + ps2 * ar
              enddo
              if (tvanp (nt) .or.newpseudo (nt) ) then
                 do jh = ih + 1, nh (nt)
                    jkb = ijkb0 + jh
                    ps1 = 2.d0 * deeq (ih, jh, na, current_spin)
                    ps2 = 2.d0 * qq (ih, jh, nt)
                    do ig = 1, npw
                       ar = vkb (ig, ikb) *conjg( vkb (ig, jkb))
                       h_diag (ig) = h_diag (ig) + ps1 * ar
                       s_diag (ig) = s_diag (ig) + ps2 * ar
                    enddo
                 enddo
              endif
           enddo
           ijkb0 = ijkb0 + nh (nt)
        endif
     enddo
  enddo
  return
end subroutine usnldiag

