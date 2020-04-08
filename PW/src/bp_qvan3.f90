!
! Copyright (C) 2004 Vanderbilt's group at Rutgers University, NJ
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! Modified by PG - Oct.2007: removed obsolete comments
!--------------------------------------------------------------------------
SUBROUTINE qvan3( iv, jv, is, qg, ylm_k, qr )
      !---------------------------------------------------------------------
      !! It calculates:
      !! \[ \text{qg} = \sum_{LM} (-I)^L \text{AP}(LM,\text{iv},
      !!                \text{jv}) \text{YR}_{LM} \text{QRAD}(\text{iv},
      !!                \text{jv},L,\text{is}) \]
      !
      ! It calculates: qg = SUM_LM (-I)^L AP(LM,iv,jv) YR_LM QRAD(iv,jv,L,is)
      !
      USE kinds,       ONLY: DP
      USE ions_base,   ONLY: ntyp => nsp
      USE uspp_param,  ONLY: lmaxq, nbetam
      USE uspp,        ONLY: nlx, lpl, lpx, ap, indv, nhtol, nhtolm
      !
      IMPLICIT NONE
      !
      INTEGER :: iv
      !! beta function index
      INTEGER :: jv
      !! beta function index
      INTEGER :: is
      !! atomic type
      COMPLEX(DP) :: qg
      !! output: see routine comments
      REAL(DP) :: ylm_k(lmaxq*lmaxq)
      !! q-space real spherical harmonics at dk [\(Y_{LM}\)]
      REAL(DP) :: qr(nbetam,nbetam,lmaxq,ntyp)
      !! Bessel transform of \(Q_{ij}(|r|)\) at dk [\(Q_{ij}^L(|r|)\)]
      !
      ! ... local variables
      !
      COMPLEX(DP) :: sig
      INTEGER :: ivs, jvs, ivl, jvl, lp, l, i
      !
      !
      ivs = indv(iv,is)
      jvs = indv(jv,is)
      ivl = nhtolm(iv,is)
      jvl = nhtolm(jv,is)
      !
      IF (ivs > nbetam .OR. jvs > nbetam) &
           CALL errore( ' qvan3 ', ' wrong dimensions (1)', MAX(ivs,jvs) )
      IF (ivl > nlx .OR. jvl > nlx) &
           CALL errore( ' qvan3 ', ' wrong dimensions (2)', MAX(ivl,jvl) )
      !
      qg = (0.0_DP,0.0_DP)
      !
!odl                  Write(*,*) 'QVAN3  --  ivs jvs = ',ivs,jvs
!odl                  Write(*,*) 'QVAN3  --  ivl jvl = ',ivl,jvl
      DO i = 1, lpx(ivl,jvl)
!odl                  Write(*,*) 'QVAN3  --  i = ',i
        lp = lpl(ivl,jvl,i)
!odl                  Write(*,*) 'QVAN3  --  lp = ',lp
        !
        ! ... EXTRACTION OF ANGULAR MOMENT L FROM LP:
        !
        IF (lp == 1) THEN
          l = 1
        ELSEIF ((lp >=  2).AND.(lp <=  4)) THEN
          l = 2
        ELSEIF ((lp >=  5).AND.(lp <=  9)) THEN
          l = 3
        ELSEIF ((lp >= 10).AND.(lp <= 16)) THEN
          l = 4
        ELSEIF ((lp >= 17).AND.(lp <= 25)) THEN
          l = 5
        ELSEIF ((lp >= 26).AND.(lp <= 36)) THEN
          l = 6
        ELSEIF ((lp >= 37).AND.(lp <= 49)) THEN
          l = 7
        ELSEIF (lp > 49) THEN
          CALL errore( ' qvan3 ',' l not programmed ', lp )
        ENDIF
        !
        sig = (0.0_DP,-1.0_DP)**(l-1)
        sig = sig * ap(lp,ivl,jvl)
        !
!odl                  Write(*,*) 'QVAN3  --  sig = ',sig
        !
        ! WRITE( stdout,*) 'qvan3',ng1,LP,L,ivs,jvs
        !
        qg = qg + sig * ylm_k(lp) * qr(ivs,jvs,l,is)
        !
      ENDDO
      !
      !
      RETURN
      !
END SUBROUTINE qvan3
