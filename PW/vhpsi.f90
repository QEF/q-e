!
! Copyright (C) 2001-2007 Quantum-Espresso group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "f_defs.h"
!
!-----------------------------------------------------------------------
subroutine vhpsi (ldap, np, mp, psip, hpsi)
  !-----------------------------------------------------------------------
  !
  ! This routine computes the Hubbard potential applied to the electronic
  ! of the current k-point, the result is added to hpsi
  !
  USE kinds,     ONLY : DP
  USE ldaU,      ONLY : Hubbard_lmax, Hubbard_l, HUbbard_U, Hubbard_alpha, &
                        swfcatom
  USE lsda_mod,  ONLY : nspin, current_spin
  USE scf,       ONLY : v
  USE ions_base, ONLY : nat, ntyp => nsp, ityp
  USE basis,     ONLY : natomwfc
  USE wvfct,     ONLY : gamma_only
  USE gvect,     ONLY : gstart
  USE uspp_param,ONLY : upf
  !
  implicit none
  !
  integer :: ldap, np, mp
  complex(DP) :: psip (ldap, mp), hpsi (ldap, mp)
  !
  integer :: ibnd, i, na, nt, n, counter, m1, m2, l
  integer, allocatable ::  offset (:)
  ! offset of localized electrons of atom na in the natomwfc ordering
  complex(DP) :: ZDOTC, temp
  real(DP), external :: DDOT
  complex(DP), allocatable ::  proj (:,:)
  !
  allocate ( offset(nat), proj(natomwfc,mp) ) 
  counter = 0  
  do na = 1, nat  
     nt = ityp (na)  
     do n = 1, upf(nt)%nwfc
        if (upf(nt)%oc (n) >= 0.d0) then  
           l = upf(nt)%lchi (n)  
           if (l.eq.Hubbard_l(nt)) offset (na) = counter  
           counter = counter + 2 * l + 1  
        endif
     enddo
  enddo
  !
  if (counter.ne.natomwfc) call errore ('vhpsi', 'nstart<>counter', 1)
  do ibnd = 1, mp
     do i = 1, natomwfc
        if (gamma_only) then
           proj (i, ibnd) = 2.d0 * &
                DDOT(2*np, swfcatom (1, i), 1, psip (1, ibnd), 1) 
           if (gstart.eq.2) proj (i, ibnd) = proj (i, ibnd) - &
                swfcatom (1, i) * psip (1, ibnd)
        else
           proj (i, ibnd) = ZDOTC (np, swfcatom (1, i), 1, psip (1, ibnd), 1)
        endif
     enddo
  enddo
#ifdef __PARA
  call reduce (2 * natomwfc * mp, proj)
#endif
  do ibnd = 1, mp  
     do na = 1, nat  
        nt = ityp (na)  
        if (Hubbard_U(nt).ne.0.d0 .or. Hubbard_alpha(nt).ne.0.d0) then  
           do m1 = 1, 2 * Hubbard_l(nt) + 1 
              temp = 0.d0
              do m2 = 1, 2 * Hubbard_l(nt) + 1 
                 temp = temp + v%ns( m1, m2, current_spin, na) * &
                                     proj (offset(na)+m2, ibnd)
              enddo
              if (gamma_only) then
                 call DAXPY (2*np, temp, swfcatom(1,offset(na)+m1), 1, &
                                    hpsi(1,ibnd),              1)
              else
                 call ZAXPY (np, temp, swfcatom(1,offset(na)+m1), 1, &
                                    hpsi(1,ibnd),              1)
              endif
           enddo
        endif
     enddo
  enddo
  deallocate (offset, proj)
  return

end subroutine vhpsi

