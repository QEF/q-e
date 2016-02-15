!
! Copyright (C) 2001-2007 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
subroutine addcore(mode, drhoc)
  use kinds,     only: DP
  use modes,     only: u
  use qpoint,    only: xq
  use nlcc_ph,   only: drc
  use fft_base,  only: dfftp
  use gvect,     only: ngm
  use ions_base, only: ntyp => nsp
  implicit none
  integer, intent (IN) :: mode
  !  input: the mode
  complex(DP), intent(OUT) :: drhoc (dfftp%nnr)
  ! output: the change of the core charge
  CALL addcore_ofq (xq, mode, u, drc, drhoc, +1._dp)
end subroutine
!
subroutine addcore_ofq (xq, mode, u, drc, drhoc, csign)
  !
  !    This routine computes the change of the core charge
  !    when the atoms moves along the given mode
  !
  !
  USE kinds, only : DP
  use uspp_param, only: upf
  use ions_base, only: nat, ityp, ntyp => nsp
  use cell_base, only: tpiba
  use fft_base,  only: dfftp
  use fft_interfaces, only: invfft
  use gvect, only: ngm, nl, mill, eigts1, eigts2, eigts3, g
  use qpoint, only: eigqts
  use nlcc_ph, only: nlcc_any
  implicit none

  !
  !   The dummy variables
  !
  integer, intent (IN) :: mode              ! the mode
  REAL(DP),INTENT(in)    :: xq(3)           ! the q vector
  COMPLEX(DP),INTENT(in) :: u(3*nat, 3*nat) ! dispalcement patterns
  COMPLEX(DP),INTENT(in) :: drc(ngm, ntyp)  ! core rho without atomic factor
  REAL(DP),INTENT(in) :: csign ! add core with +1, subtract with -1, or anything in between
  !
  complex(DP), intent(OUT) :: drhoc (dfftp%nnr)
  ! output: the change of the core charge
  !
  !   Local variables
  !
  integer :: nt, ig, mu, na
  complex(DP) :: fact, gu, gu0, u1, u2, u3, gtau
  !
  !
  if (.not.nlcc_any) return
  !
  ! compute the derivative of the core charge  along the given mode
  !
  drhoc(:) = (0.d0, 0.d0)
  do na = 1, nat
     nt = ityp (na)
     if (upf(nt)%nlcc) then
        fact = csign*tpiba * (0.d0, -1.d0) * eigqts (na)
        mu = 3 * (na - 1)
        if ( abs (u (mu + 1, mode) ) + &
             abs (u (mu + 2, mode) ) + &
             abs (u (mu + 3, mode) ) > 1.0d-12) then
           u1 = u (mu + 1, mode)
           u2 = u (mu + 2, mode)
           u3 = u (mu + 3, mode)
           gu0 = xq (1) * u1 + xq (2) * u2 + xq (3) * u3
           do ig = 1, ngm
              gtau = eigts1 (mill (1,ig), na) * eigts2 (mill (2,ig), na) &
                   * eigts3 (mill (3,ig), na)
              gu = gu0 + g (1, ig) * u1 + g (2, ig) * u2 + g (3, ig) &
                   * u3
              drhoc (nl (ig) ) = drhoc (nl (ig) ) + drc (ig, nt) * gu * &
                   fact * gtau
           enddo
        endif
     endif
  enddo
  !
  !   transform to real space
  !
  CALL invfft ('Dense', drhoc, dfftp)
  !
  return

end subroutine addcore_ofq
