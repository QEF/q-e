!
! Copyright (C) 2001-2013 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
SUBROUTINE write_io_header(filplot, iunplot, title, nr1x, nr2x, nr3x, &
           nr1, nr2, nr3, nat, ntyp, ibrav, celldm, at, gcutm, dual, ecutwfc, &
           nkstot,nbnd,natomwfc)
   !-----------------------------------------------------------------------

   USE kinds, ONLY: DP
   USE ions_base, ONLY : zv, atm, tau, ityp
   USE noncollin_module, ONLY: noncolin
   USE spin_orb, ONLY: lspinorb

   IMPLICIT NONE
   CHARACTER (len=*) :: filplot
   CHARACTER (len=*) :: title
   INTEGER :: nr1x, nr2x, nr3x, nr1, nr2, nr3, nat, ntyp, ibrav
   REAL(DP) :: celldm (6), gcutm, dual, ecutwfc, at(3,3)
   INTEGER :: iunplot, ios, na, nt, i
   INTEGER :: nkstot,nbnd,natomwfc
   !
   IF (filplot == ' ') CALL errore ('write_io_h', 'filename missing', 1)

   OPEN (UNIT = iunplot, FILE = filplot, FORM = 'formatted', &
         STATUS = 'unknown', ERR = 101, IOSTAT = ios)
   101     CALL errore ('write_io_header', 'opening file '//trim(filplot), abs (ios) )
   WRITE (iunplot, '(a)') title
   WRITE (iunplot, '(8i8)') nr1x, nr2x, nr3x, nr1, nr2, nr3, nat, ntyp
   WRITE (iunplot, '(i6,6f12.8)') ibrav, celldm
   IF  (ibrav == 0) THEN
       WRITE ( iunplot, * ) at(:,1)
       WRITE ( iunplot, * ) at(:,2)
       WRITE ( iunplot, * ) at(:,3)
   ENDIF
   WRITE (iunplot, '(3f20.10,i6)') gcutm, dual, ecutwfc, 9
   WRITE (iunplot, '(i4,3x,a2,3x,f5.2)') &
         (nt, atm (nt), zv (nt), nt=1, ntyp)
   WRITE (iunplot, '(i4,3x,3f15.9,3x,i2)') (na, &
         (tau (i, na), i = 1, 3), ityp (na), na = 1, nat)
   WRITE (iunplot, '(3i8)') natomwfc, nkstot, nbnd
   WRITE (iunplot, '(2l5)') noncolin, lspinorb

   RETURN
END SUBROUTINE write_io_header
