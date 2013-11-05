!
! Copyright (C) 2001-2009 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
subroutine ewald_dipole (tens,dipole)
  !-----------------------------------------------------------------------
  !
  ! Calculates the ewald field on each atom due to the presence of dipole, or
  ! the electic field on each atom due to the ionic charge of other atoms,
  ! with both G- and R-space terms.
  ! Determines optimal alpha. Should hopefully work for any structure.
  !
  !
  USE kinds ,     ONLY : dp
  USE gvect ,     ONLY : gcutm, gstart, ngm, g, gg
  USE constants , ONLY : tpi, e2, fpi, pi
  USE cell_base , ONLY : tpiba2, omega, alat, at, bg
  USE ions_base,  ONLY : nat, ntyp => nsp, ityp, tau
  USE vlocal ,    ONLY : strf
  USE mp_bands,   ONLY : intra_bgrp_comm
  USE mp,         ONLY : mp_sum
  !
  implicit none
  !
  real(DP) :: dipole(ntyp),charge, eta, arg, upperbound, temp
  complex(DP) :: tens(nat,3,3)
  complex(DP) :: rhon
  real(DP), external :: qe_erfc
  complex(DP), allocatable:: ewaldg(:,:,:), ewaldr(:,:,:)
  integer :: alpha, beta, na, ng, nt, ipol, nb, nrm, nr

  integer, parameter :: mxr = 50
  real (DP) :: r(3,mxr), r2(mxr), rmax, rr, dtau(3)
  real (DP) :: expcoeff
  complex(DP) :: carg, recarg, recarg_dgg

  allocate (ewaldg(nat,3,3))
  allocate (ewaldr(nat,3,3))

  ewaldg=(0.d0,0.d0)  
  ewaldr=(0.d0,0.d0)


!  e2=1.d0 !hartree
  charge = 0.d0
  do na = 1, nat
     charge = charge+dipole (ityp (na) )
  enddo
  eta = 2.9d0
  do
    eta = eta - 0.1d0
    !
    ! choose alpha in order to have convergence in the sum over G
    ! upperbound is a safe upper bound for the error in the sum over G
    !
    if (eta.le.0.d0) call errore ('ewald_dipole', 'optimal eta not found', 1)
    upperbound = 2.d0 * charge**2 * sqrt (2.d0 * eta / tpi) &
                      * qe_erfc ( sqrt (tpiba2 * gcutm / 4.d0 / eta) )
    if (upperbound.le.1.0d-7) exit
  enddo
  !
  ! G-space sum here.

  do ng = gstart, ngm
     rhon = (0.d0, 0.d0)
     expcoeff = exp ( - gg (ng) * tpiba2 * 0.25d0 / eta )
     do nt = 1, ntyp
        rhon = rhon + dipole (nt) * CONJG(strf (ng, nt) )
     enddo
     do na=1, nat
        arg = (g (1, ng) * tau (1, na) + g (2, ng) * tau (2, na) &
             + g (3, ng) * tau (3, na) ) * tpi
        carg = CMPLX(cos(arg), -sin(arg),kind=DP)
        recarg = rhon*expcoeff*carg
        recarg_dgg = recarg / gg(ng)
        do alpha = 1,3
           do beta=1,3
              ewaldg(na , alpha, beta) = ewaldg(na, alpha, beta) &
                - recarg_dgg * g(alpha,ng) * g(beta,ng)
           enddo

           ewaldg(na , alpha, alpha) = ewaldg(na, alpha, alpha) &
                                      + 1.d0/3.d0 * recarg
        enddo
     enddo
  enddo
  ewaldg = e2 / 2.d0 * fpi / omega * ewaldg !Temp to compare with paratec
!  ewaldg = e2 * fpi / omega * ewaldg
  !
  call mp_sum(  ewaldg, intra_bgrp_comm )
  !
  ! R-space sum here (only for the processor that contains G=0)
  !
  ewaldr = 0.d0
  if (gstart.eq.2) then
     rmax = 4.d0 / sqrt (eta) / alat
     !
     ! with this choice terms up to ZiZj*erfc(4) are counted (erfc(4)=2x10^-8
     !
    do na = 1, nat
       do nb = 1, nat
          do ipol = 1, 3
             dtau (ipol) = tau (ipol, na) - tau (ipol, nb)
          enddo
           !
           ! generates nearest-neighbors shells
           !
          call rgen (dtau, rmax, mxr, at, bg, r, r2, nrm)
           !
           ! and sum to the real space part
           !
          r = r * alat
          do nr = 1, nrm
             rr = sqrt (r2 (nr) ) * alat
             temp= dipole (ityp (na))  * ( 3.d0 / rr**3 *  qe_erfc ( sqrt (eta) * rr) &
                                          + (6.d0 * sqrt (eta/pi) * 1.d0 / rr*2 + 4.d0 * sqrt (eta**3/pi)) &
                                          * exp(-eta* rr**2))
             do alpha=1,3
                do beta=1,3
                   ewaldr(na, alpha,beta) = ewaldr(na, alpha,beta) &
                                           + temp*r(alpha,nr)*r(beta,nr) / rr**2
                enddo
                ewaldr(na, alpha,alpha)= ewaldr(na, alpha,alpha) &
                                        - 1.d0/3.d0 * temp
             enddo
          enddo
       enddo
    enddo
 endif
 ewaldr = e2 *  ewaldr
 !
 call mp_sum(  ewaldr, intra_bgrp_comm )
 ! 
 tens=ewaldg+ewaldr

end subroutine ewald_dipole
