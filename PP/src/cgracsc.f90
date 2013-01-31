!
! Copyright (C) 2001-2005 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
FUNCTION cgracsc (nkb, bec1, bec2, nhm, ntyp, nh, qq, nat, ityp, &
     npw, psi1, psi2, upf)
  !-----------------------------------------------------------------------
  !
  !     This function computes the scalar product between two wavefunction
  !     and the S matrix of the US pseudopotential: <psi1 | S | psi2 >.
  !     It assumes that the product of psi1 with all the beta functions
  !     is in bec1, and the product of psi2 is in bec2.
  !
  !
  USE kinds
  USE pseudo_types, ONLY : pseudo_upf
  USE mp_global,  ONLY : intra_bgrp_comm
  USE mp,         ONLY : mp_sum


  IMPLICIT NONE
  !
  !     here the dummy variables
  !

  INTEGER :: nkb, npw, nhm, ntyp, nat, ityp (nat), nh (ntyp)
  ! input: the number of beta functions
  ! input: the number of plane waves
  ! input: the maximum number of solid be
  ! input: the number of types of atoms
  ! input: the number of atoms
  ! input: the type of each atom
  ! input: the number of beta for each ty

  COMPLEX(DP) :: bec1 (nkb), bec2 (nkb), psi1 (npw), psi2 (npw), &
       cgracsc
  ! input: the product of beta and psi1
  ! input: the product of beta and psi2
  ! input: the first wavefunction
  ! input: the second wavefunction
  ! output: the value of the scalar produ

  real(DP) :: qq (nhm, nhm, ntyp)
  ! input: the q values defining S
  TYPE(pseudo_upf) :: upf (ntyp)
  ! input: if true the pseudo is vanderb
  !
  !    Here the local variables
  !

  INTEGER :: ikb, jkb, na, np, ijkb0, ih, jh
  ! counter on total beta functions
  ! counter on total beta functions
  ! counter on atoms
  ! the pseudopotential of each atom
  ! auxiliary variable to compute ikb and jkb
  ! counter on solid beta functions
  ! counter on solid beta functions

  COMPLEX(DP) :: scal, zdotc
  !
  scal = zdotc (npw, psi1, 1, psi2, 1)
  CALL mp_sum(  scal, intra_bgrp_comm )
  ijkb0 = 0
  DO np = 1, ntyp
     IF (upf(np)%tvanp ) THEN
        DO na = 1, nat
           IF (ityp (na) ==np) THEN
              DO ih = 1, nh (np)
                 ikb = ijkb0 + ih
                 DO jh = 1, nh (np)
                    jkb = ijkb0 + jh
                    scal = scal + qq (ih,jh,np)*conjg(bec1(ikb))*bec2(jkb)
                 ENDDO
              ENDDO
              ijkb0 = ijkb0 + nh (np)
           ENDIF
        ENDDO
     ELSE
        DO na = 1, nat
           IF (ityp (na) ==np) ijkb0 = ijkb0 + nh (np)
        ENDDO
     ENDIF

  ENDDO

  cgracsc = scal
  RETURN
END FUNCTION cgracsc

!
!-----------------------------------------------------------------------
FUNCTION cgracsc_nc (nkb, bec1, bec2, nhm, ntyp, nh, nat, ityp, &
     npw, npol, psi1, psi2, upf)
  !-----------------------------------------------------------------------
  !
  !     This function computes the scalar product between two wavefunction
  !     and the S matrix of the US pseudopotential: <psi1 | S | psi2 >.
  !     It assumes that the product of psi1 with all the beta functions
  !     is in bec1, and the product of psi2 is in bec2.
  !
  !
  USE kinds
  USE uspp, ONLY: qq, qq_so
  USE spin_orb, ONLY: lspinorb
  USE pseudo_types, ONLY : pseudo_upf
  USE mp_global,  ONLY : intra_bgrp_comm
  USE mp,         ONLY : mp_sum
  IMPLICIT NONE
  !
  !     here the dummy variables
  !

  INTEGER :: nkb, npw, npol, nhm, ntyp, nat, ityp (nat), nh (ntyp)
  ! input: the number of beta functions
  ! input: the number of plane waves
  ! input: the maximum number of solid be
  ! input: the number of types of atoms
  ! input: the number of atoms
  ! input: the type of each atom
  ! input: the number of beta for each ty

  COMPLEX(DP) :: bec1 (nkb,npol), bec2 (nkb,npol), &
                      psi1 (npw,npol), psi2 (npw,npol), cgracsc_nc
  ! input: the product of beta and psi1
  ! input: the product of beta and psi2
  ! input: the first wavefunction
  ! input: the second wavefunction
  ! output: the value of the scalar produ

  TYPE(pseudo_upf) :: upf (ntyp)
  ! input: if true the pseudo is vanderb
  !
  !    Here the local variables
  !

  INTEGER :: ikb, jkb, na, np, ijkb0, ih, jh, ipol, jpol, ijh
  ! counter on total beta functions
  ! counter on total beta functions
  ! counter on atoms
  ! the pseudopotential of each atom
  ! auxiliary variable to compute ikb and jkb
  ! counter on solid beta functions
  ! counter on solid beta functions

  COMPLEX(DP) :: scal, zdotc
  !
  scal = zdotc (npw*npol, psi1, 1, psi2, 1)
  CALL mp_sum(  scal, intra_bgrp_comm )
  ijkb0 = 0
  DO np = 1, ntyp
     IF (upf(np)%tvanp ) THEN
        DO na = 1, nat
           IF (ityp (na) ==np) THEN
              DO ih = 1, nh (np)
                 ikb = ijkb0 + ih
                 DO jh = 1, nh (np)
                    jkb = ijkb0 + jh
                    IF (lspinorb) THEN
                       ijh=0
                       DO ipol=1,npol
                          DO jpol=1,npol
                            ijh=ijh+1
                            scal=scal+qq_so(ih,jh,ijh,np)* &
                                      conjg(bec1(ikb,ipol))*bec2(jkb,jpol)
                          ENDDO
                       ENDDO
                    ELSE
                       DO ipol=1,npol
                          scal=scal+qq(ih,jh,np)* &
                               conjg(bec1(ikb,ipol))*bec2(jkb,ipol)
                       ENDDO
                    ENDIF
                 ENDDO
              ENDDO
              ijkb0 = ijkb0 + nh (np)
           ENDIF
        ENDDO
     ELSE
        DO na = 1, nat
           IF (ityp (na) ==np) ijkb0 = ijkb0 + nh (np)
        ENDDO
     ENDIF
  ENDDO

  cgracsc_nc = scal
  RETURN
END FUNCTION cgracsc_nc
