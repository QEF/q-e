!
! Copyright (C) 2001-2019 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------------
SUBROUTINE lr_addus_dvpsi (npwq, ik, psi, dvpsi)
  !----------------------------------------------------------------------
  !---------------------------------------------------------------------------
  !
  ! ... Calculate the ultrasoft term of the perturbation exp(iq*r),
  ! ... and then sum up the input wavefunction and the ultrasoft term.
  !
  ! ... input:
  ! ...    ik     given k point
  ! ...    npwq   true dimension of psi
  ! ...    psi    input array
  !
  ! Written by Iurii Timrov (2015)
  ! Generalized to the relativistic case by Andrea Dal Corso (2018)
  !
  USE kinds,            ONLY : DP
  USE uspp_param,       ONLY : upf, nh
  USE uspp,             ONLY : vkb, okvan
  USE lsda_mod,         ONLY : lsda, current_spin, isk
  USE ions_base,        ONLY : ntyp => nsp, nat, ityp
  USE wvfct,            ONLY : nbnd, npwx
  USE noncollin_module, ONLY : noncolin, npol
  USE qpoint,           ONLY : ikks
  USE lr_variables,     ONLY : intq, intq_nc
  USE lrus,             ONLY : becp1
  
  IMPLICIT NONE
  !
  !   The dummy variables
  !
  COMPLEX(DP), INTENT(in) :: psi(npwx*npol,nbnd)
  ! input: wavefunction u_n,k
  COMPLEX(DP), INTENT(out) :: dvpsi(npwx*npol,nbnd)
  ! output: sum of the input wavefunction and the USPP term
  INTEGER :: ik, npwq
  ! input: the k point
  ! input: number of plane waves at the q point
  
  !   And the local variables
  !
  INTEGER :: na, nt, ibnd, ih, jh, ijkb0, ikk, ikb, jkb, is, js, ijs
  ! counter on atoms
  ! counter on atomic types
  ! counter on bands
  ! counter on beta functions
  ! counter on beta functions
  ! auxiliary variable for indexing
  ! counter on the k points
  ! counter on vkb
  ! counter on vkb
  COMPLEX(DP) :: sum0, sum_nc(npol)
  ! auxiliary variable

  IF (.NOT.okvan) RETURN
  CALL start_clock ('lr_addus_dvpsi')
  dvpsi=psi
  ikk = ikks(ik)
  IF (lsda) current_spin = isk (ikk)
  ijkb0 = 0
  DO nt = 1, ntyp
     IF (upf(nt)%tvanp  ) THEN
        DO na = 1, nat
           IF (ityp (na)==nt) THEN
              !
              !   we multiply the integral for the becp term and the beta_n
              !
              DO ibnd = 1, nbnd
                 DO ih = 1, nh (nt)
                    ikb = ijkb0 + ih
                    IF (noncolin) THEN
                       sum_nc = (0.d0, 0.d0)
                    ELSE
                       sum0 = (0.d0, 0.d0)
                    END IF
                    DO jh = 1, nh (nt)
                       jkb = ijkb0 + jh
                       IF (noncolin) THEN
                          ijs=0
                          DO is=1,npol
                             DO js=1,npol
                                ijs=ijs+1
                                sum_nc(is)=sum_nc(is)+         &
                                     intq_nc(ih,jh,na,ijs)*    &
                                     becp1(ik)%nc(jkb, js, ibnd)
                             ENDDO
                          ENDDO
                       ELSE
                          sum0 = sum0 + intq (ih, jh, na)*     &
                                   becp1(ik)%k(jkb, ibnd)
                       ENDIF
                    ENDDO
                    IF (noncolin) THEN
                       CALL zaxpy(npwq,sum_nc(1),vkb(1,ikb),1,dvpsi(1,ibnd),1)
                       CALL zaxpy(npwq,sum_nc(2),vkb(1,ikb),1, &
                                                 dvpsi(1+npwx,ibnd),1)
                    ELSE
                       CALL zaxpy(npwq,sum0,vkb(1,ikb),1,dvpsi(1,ibnd),1)
                    ENDIF
                 ENDDO
              ENDDO
              ijkb0 = ijkb0 + nh (nt)
           ENDIF
        ENDDO
     ELSE
        DO na = 1, nat
           IF (ityp (na)==nt) ijkb0 = ijkb0 + nh (nt)
        ENDDO
     ENDIF
  ENDDO

  CALL stop_clock ('lr_addus_dvpsi')
  RETURN
END SUBROUTINE lr_addus_dvpsi
