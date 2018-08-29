!
! Copyright (C) 2001-2018 Quantum_ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
SUBROUTINE hp_calc_chi
!-----------------------------------------------------------------------
  !
  ! This routine computes the non-interacting and interaction response
  ! matrices chi0 and chi, respectively, from dns0_tot and dnsscf_tot. 
  ! See Eq. (41) in Ref. [1].
  ! [1] Phys. Rev. B 98, 085127 (2018)
  !
  USE ions_base,    ONLY : nat
  USE io_global,    ONLY : stdout
  USE ldaU_hp,      ONLY : nqsh, nath_sc, nah_pert, chi0, chi, &
                           dnsscf_tot, dns0_tot
  !
  IMPLICIT NONE
  !
  CALL start_clock('hp_calc_chi') 
  !
  ! Compute and write to file chi0 and chi
  ! for a given perturbation of atom nah_pert
  !
  CALL calcchi (dns0_tot,   chi0, 'chi0')
  CALL calcchi (dnsscf_tot, chi,  'chi')
  !
  CALL stop_clock('hp_calc_chi')
  !
  RETURN
  !
CONTAINS
  !
SUBROUTINE calcchi (dns_, chi_, name_)
  !
  ! Compute the trace of the response occupation matrices
  !
  USE kinds,        ONLY : DP
  USE io_files,     ONLY : prefix
  USE ions_base,    ONLY : ntyp => nsp, ityp
  USE ldaU,         ONLY : Hubbard_l, is_hubbard, Hubbard_lmax
  USE lsda_mod,     ONLY : nspin
  USE constants,    ONLY : rytoev
  !
  IMPLICIT NONE
  !
  CHARACTER(len=4), INTENT(IN)    :: name_
  COMPLEX(DP),      INTENT(IN)    :: dns_(2*Hubbard_lmax+1, 2*Hubbard_lmax+1, nspin, nat, nqsh)
  REAL(DP),         INTENT(INOUT) :: chi_(nath_sc, nat)
  !
  ! Local variables
  !
  COMPLEX(DP) :: trace_dns(2), trace_dns_tot
  INTEGER  :: na, nt, is, m, icell, na_sc
  !
  na_sc = 0
  !
  DO icell = 1, nqsh
     !
     DO na = 1, nat
        !
        nt = ityp(na)
        !
        IF ( is_hubbard(nt) ) THEN
           !
           na_sc = na_sc + 1
           !
           trace_dns(:)  = 0.d0
           trace_dns_tot = 0.d0
           !
           ! Divide by rytoev -> conversion of units from Ry to eV
           !
           DO is = 1, nspin
              DO m = 1, 2 * Hubbard_l(nt) + 1
                 trace_dns(is) = trace_dns(is) + dns_(m,m,is,na,icell)/rytoev
              ENDDO
              trace_dns_tot = trace_dns_tot + trace_dns(is)
           ENDDO
           !
           ! If nspin=1, multiply by a factor of 2 due to spin degeneracy
           !
           IF (nspin.EQ.1) trace_dns_tot = 2.0d0 * trace_dns_tot
           !
           chi_(na_sc, nah_pert) = DBLE(trace_dns_tot)
           !
        ENDIF
        !
     ENDDO
     !
  ENDDO
  !
  IF (na_sc.NE.nath_sc) CALL errore( 'hp_calc_chi', "Mismatch in the number of atoms", 1)
  !
  RETURN
  !
END SUBROUTINE calcchi

END SUBROUTINE hp_calc_chi
