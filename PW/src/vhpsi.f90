!
! Copyright (C) 2001-2007 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
subroutine vhpsi (ldap, np, mps, psip, hpsi)
  !-----------------------------------------------------------------------
  !
  ! This routine computes the Hubbard potential applied to the electronic
  ! of the current k-point, the result is added to hpsi
  !
  USE kinds,     ONLY : DP
  USE becmod,    ONLY : bec_type, calbec, allocate_bec_type, deallocate_bec_type
  USE ldaU,      ONLY : Hubbard_lmax, Hubbard_l, HUbbard_U, Hubbard_alpha, &
                        swfcatom, oatwfc
  USE lsda_mod,  ONLY : nspin, current_spin
  USE scf,       ONLY : v
  USE ions_base, ONLY : nat, ntyp => nsp, ityp
  USE basis,     ONLY : natomwfc
  USE gvect,     ONLY : gstart
  USE control_flags, ONLY : gamma_only
  USE mp_global, ONLY: intra_pool_comm
  USE mp,        ONLY: mp_sum
  !
  implicit none
  !
  integer, intent (in) :: ldap, np, mps
  complex(DP), intent(in) :: psip (ldap, mps)
  complex(DP), intent(inout) :: hpsi (ldap, mps)
  !
  integer :: ibnd, na, nt, m1, m2
  complex(DP) :: temp
  type (bec_type) :: proj
  !
  ! Offset of atomic wavefunctions initialized in setup and stored in oatwfc
  !
  call allocate_bec_type ( natomwfc,mps, proj )
  CALL calbec (np, swfcatom, psip, proj)
  do ibnd = 1, mps  
     do na = 1, nat  
        nt = ityp (na)  
        if (Hubbard_U(nt).ne.0.d0 .or. Hubbard_alpha(nt).ne.0.d0) then  
           do m1 = 1, 2 * Hubbard_l(nt) + 1 
              temp = 0.d0
              if (gamma_only) then
                 do m2 = 1, 2 * Hubbard_l(nt) + 1 
                    temp = temp + v%ns( m1, m2, current_spin, na) * &
                                    proj%r(oatwfc(na)+m2, ibnd)
                 enddo
                 call daxpy (2*np, temp, swfcatom(1,oatwfc(na)+m1), 1, &
                                    hpsi(1,ibnd),              1)
              else
                 do m2 = 1, 2 * Hubbard_l(nt) + 1 
                    temp = temp + v%ns( m1, m2, current_spin, na) * &
                                    proj%k(oatwfc(na)+m2, ibnd)
                 enddo
                 call zaxpy (np, temp, swfcatom(1,oatwfc(na)+m1), 1, &
                                    hpsi(1,ibnd),              1)
              endif
           enddo
        endif
     enddo
  enddo
  call deallocate_bec_type (proj)
  return

end subroutine vhpsi

