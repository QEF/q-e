!
! Copyright (C) 2001-2018 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!------------------------------------------------------------
SUBROUTINE hp_write_dnsq(iq)
!------------------------------------------------------------
  !
  ! Write dns0 and dnsscf for the current q to file.
  !
  USE kinds,         ONLY : DP
  USE io_global,     ONLY : ionode
  USE ions_base,     ONLY : nat, ityp, atm
  USE lsda_mod,      ONLY : nspin
  USE io_files,      ONLY : prefix, tmp_dir
  USE ldaU,          ONLY : Hubbard_lmax, Hubbard_l, is_hubbard
  USE ldaU_hp,       ONLY : nah_pert, dns0, dnsscf, x_q
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: iq
  ! index of the q point
  !
  INTEGER :: iunitdnsq ! unit number
  CHARACTER(len=50) :: filenamednsq
  CHARACTER(len=256) :: tempfile
  CHARACTER(len=6), EXTERNAL :: int_to_char
  INTEGER, EXTERNAL :: find_free_unit
  !
  IF (ionode) THEN
     !
     iunitdnsq = find_free_unit()
     filenamednsq = TRIM(prefix) // TRIM(".dns.") // TRIM("pert_") // TRIM(int_to_char(nah_pert)) // &
                    TRIM(".q_") // TRIM(int_to_char(iq)) // TRIM(".dat")
     tempfile = trim(tmp_dir) // trim(filenamednsq)
     !
     OPEN(iunitdnsq, file = tempfile, form = 'formatted', status = 'unknown')
     !
     WRITE(iunitdnsq, '(1x,"Perturbed atom #",1x,i3,1x,", type : ",a4)') nah_pert, atm(ityp(nah_pert))
     WRITE(iunitdnsq, '(1x,"q #",i4," = (", 3F12.7, " )")') iq, x_q(:,iq)
     WRITE(iunitdnsq, *)
     !
     CALL write_dnsq (dns0(:,:,:,:,iq),   'dns0  ')
     CALL write_dnsq (dnsscf(:,:,:,:,iq), 'dnsscf')
     !
     CLOSE(iunitdnsq)
     !
  ENDIF   
  !
  RETURN
  !
CONTAINS
  !
SUBROUTINE write_dnsq (dns, name_)
  !  
  IMPLICIT NONE
  !
  CHARACTER(len=6), INTENT(IN) :: name_
  COMPLEX(DP),      INTENT(IN) :: dns(2*Hubbard_lmax+1, 2*Hubbard_lmax+1, nspin, nat)
  !
  INTEGER :: na, nt, is, m1, m2
  !
  WRITE(iunitdnsq,'(1x,"Response occupation matrix ", a6, " :")') TRIM(name_)
  DO na = 1, nat
     nt = ityp(na)
     IF (is_hubbard(nt)) THEN
        DO is = 1, nspin
           WRITE(iunitdnsq,'(1x,"Hubbard atom",1x,i2,2x,"spin",1x,i2)') na, is
           WRITE(iunitdnsq,'(1x,"row #",2x,"column #",6x,"Re(",a6,")",12x,"Im(",a6,")")') &
                                                              TRIM(name_), TRIM(name_)
           DO m1 = 1, 2*Hubbard_l(nt)+1
              DO m2 = 1, 2*Hubbard_l(nt)+1
                 WRITE(iunitdnsq,'(1x,i2,6x,i2,4x,f20.15,2x,f20.15)') &
                      m1, m2, DBLE(dns(m1,m2,is,na)), AIMAG(dns(m1,m2,is,na))
              ENDDO
           ENDDO
        ENDDO
     ENDIF
  ENDDO
  WRITE(iunitdnsq,*)
  !
  RETURN
  !
END SUBROUTINE write_dnsq
  !
END SUBROUTINE hp_write_dnsq
