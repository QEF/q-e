!
! Copyright (C) 2001-2018 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!------------------------------------------------------------
SUBROUTINE hp_read_dnsq()
!------------------------------------------------------------
  !
  ! Read dns0 and dnsscf for all q points from file.
  !
  USE kinds,         ONLY : DP
  USE ions_base,     ONLY : nat, ityp
  USE io_global,     ONLY : stdout
  USE io_files,      ONLY : prefix, tmp_dir
  USE lsda_mod,      ONLY : nspin
  USE ldaU,          ONLY : Hubbard_lmax, Hubbard_l, is_hubbard
  USE ldaU_hp,       ONLY : nah_pert, dns0, dnsscf, nqs, tmp_dir_hp
  !
  IMPLICIT NONE
  !
  INTEGER :: iq
  INTEGER :: iunitdnsq ! unit number
  CHARACTER(len=50) :: filenamednsq
  CHARACTER(len=256) :: tempfile
  CHARACTER(len=6), EXTERNAL :: int_to_char
  LOGICAL :: exst
  INTEGER, EXTERNAL :: find_free_unit
  !
  tmp_dir = tmp_dir_hp
  !
  DO iq = 1, nqs
     !
     iunitdnsq = find_free_unit()
     filenamednsq = TRIM(prefix) // TRIM(".dns.") // TRIM("pert_") // TRIM(int_to_char(nah_pert)) // &
                    TRIM(".q_") // TRIM(int_to_char(iq)) // TRIM(".dat")
     tempfile = trim(tmp_dir) // trim(filenamednsq)
     !
     INQUIRE (file = tempfile, exist = exst)
     IF (.NOT.exst) THEN
         WRITE( stdout, * ) "    WARNING: " // TRIM(filenamednsq) // " does not exist !!!"
         WRITE( stdout, * ) "    Check the folder: ", TRIM(tmp_dir)
         CALL errore('hp_read_dnsq','Missing file',1)
     ENDIF
     !
     OPEN(iunitdnsq, file = tempfile, form = 'formatted', status = 'unknown')
     !
     CALL read_dns(dns0(:,:,:,:,iq), dnsscf(:,:,:,:,iq))
     !
     CLOSE(iunitdnsq)
     !
  ENDDO
  !
  RETURN
  !
CONTAINS
  !
SUBROUTINE read_dns(dns0_, dnsscf_)
  ! 
  ! Read dns0 and dnsscf for a specific q point 
  ! 
  IMPLICIT NONE
  !
  COMPLEX(DP), INTENT(OUT) :: dns0_  (2*Hubbard_lmax+1, 2*Hubbard_lmax+1, nspin, nat), &
                              dnsscf_(2*Hubbard_lmax+1, 2*Hubbard_lmax+1, nspin, nat)
  !
  INTEGER :: na, nt, is, m1, m2, m1_, m2_
  REAL(DP) :: re_dns, im_dns
  !
  READ(iunitdnsq, *)
  READ(iunitdnsq, *)
  READ(iunitdnsq, *)
  !
  ! Read dns0
  READ(iunitdnsq, *)
  DO na = 1, nat
     nt = ityp(na)
     IF (is_hubbard(nt)) THEN
        DO is = 1, nspin
           READ(iunitdnsq, *)
           READ(iunitdnsq, *)
           DO m1 = 1, 2*Hubbard_l(nt)+1
              DO m2 = 1, 2*Hubbard_l(nt)+1
                 READ(iunitdnsq, *) m1_, m2_, re_dns, im_dns
                 dns0_(m1,m2,is,na) = CMPLX(re_dns, im_dns, kind=DP)
              ENDDO
           ENDDO
        ENDDO
     ENDIF
  ENDDO
  !
  ! Read dnsscf
  READ(iunitdnsq, *)
  READ(iunitdnsq, *)
  DO na = 1, nat
     nt = ityp(na)
     IF (is_hubbard(nt)) THEN
        DO is = 1, nspin
           READ(iunitdnsq, *)
           READ(iunitdnsq, *)
           DO m1 = 1, 2*Hubbard_l(nt)+1
              DO m2 = 1, 2*Hubbard_l(nt)+1
                 READ(iunitdnsq, *) m1_, m2_, re_dns, im_dns
                 dnsscf_(m1,m2,is,na) = CMPLX(re_dns, im_dns, kind=DP)
              ENDDO
           ENDDO
        ENDDO
     ENDIF
  ENDDO
  !
  RETURN
  !
END SUBROUTINE read_dns
  !
END SUBROUTINE hp_read_dnsq
