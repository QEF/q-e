!
! Copyright (C) 2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!---------------------------------------------------------------------
SUBROUTINE dvb_cc (nlcc,npseu,ngm,nr1,nr2,nr3,nrx1,nrx2,nrx3,  &
     nl,igtongl,rho_core,dmuxc,ga,aux,dvb_nlcc)
  !---------------------------------------------------------------------
  ! calculate the core-correction contribution to Delta V bare
  !
  IMPLICIT NONE
  INTEGER:: npseu,ngm,nr1,nr2,nr3,nrx1,nrx2,nrx3,np,ng,i
  LOGICAL :: nlcc(npseu)
  INTEGER :: nl(ngm), igtongl(ngm)
  real(8) :: rho_core(*), dmuxc(nrx1*nrx2*nrx3)
  COMPLEX(8) :: ga(ngm), dvb_nlcc(ngm), aux(nrx1*nrx2*nrx3)
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
  CALL cft3(aux,nr1,nr2,nr3,nrx1,nr2,nr3,1)
  !
  aux(:) = aux(:) * dmuxc(:)
  !
  CALL cft3(aux,nr1,nr2,nr3,nrx1,nr2,nr3,-1)
  DO ng=1,ngm
     dvb_nlcc(ng) = aux(nl(ng))
  ENDDO
  !
  RETURN
END SUBROUTINE dvb_cc

