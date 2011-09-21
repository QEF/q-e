!
! Copyright (C) 2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!---------------------------------------------------------------------
SUBROUTINE dvb_cc (nlcc,npseu,ngm,nrxx,  &
     nl,igtongl,rho_core,dmuxc,ga,aux,dvb_nlcc)
  !---------------------------------------------------------------------
  ! calculate the core-correction contribution to Delta V bare
  !
  USE kinds, ONLY : dp
  USE fft_base, ONLY : dffts, dfftp
  USE fft_interfaces, ONLY : fwfft, invfft
  IMPLICIT NONE
  INTEGER:: npseu,ngm,nrxx,np,ng,i
  LOGICAL :: nlcc(npseu)
  INTEGER :: nl(ngm), igtongl(ngm)
  real(dp) :: rho_core(*), dmuxc(nrxx)
  COMPLEX(dp) :: ga(ngm), dvb_nlcc(ngm), aux(nrxx)
  !
  DO np=1,npseu
     IF(nlcc(np)) GOTO 10
  ENDDO
  RETURN
10 CONTINUE
  !
  aux(:) = (0.d0, 0.d0)
  DO ng=1,ngm
     aux(nl(ng)) = ga(ng) * rho_core(igtongl(ng))
  ENDDO
  CALL invfft ('Dense', aux, dfftp)
  !
  aux(:) = aux(:) * dmuxc(:)
  !
  CALL fwfft ('Dense', aux, dfftp)
  DO ng=1,ngm
     dvb_nlcc(ng) = aux(nl(ng))
  ENDDO
  !
  RETURN
END SUBROUTINE dvb_cc

