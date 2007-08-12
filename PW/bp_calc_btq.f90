!
! Copyright (C) 2004 Vanderbilt's group at Rutgers University, NJ
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------
SUBROUTINE calc_btq(ql,qr_k,idbes)
  !----------------------------------------------------------------------
  !
  !   Calculates the Bessel-transform (or its derivative if idbes=1) 
  !   of the augmented qrad charges at a given ql point.
  !   Rydberg atomic units are  used.
  !
  USE kinds, ONLY: DP
  USE atom, ONLY: rgrid
  USE ions_base, ONLY : ntyp => nsp
  USE cell_base, ONLY: omega
  USE radial_grids, ONLY:  ndmx
  USE constants, ONLY: fpi
  USE uspp_param, ONLY: lmaxq, qfunc, qfcoef, nqf, rinner, lll, &
       nbeta, nbetam, kkbeta, tvanp
  !
  IMPLICIT NONE
  !
  REAL(DP)  :: ql, qr_k(nbetam,nbetam,lmaxq,ntyp)
  INTEGER :: idbes
  !
  INTEGER :: msh_bp, i, np, l
  INTEGER :: ilmin, ilmax, iv, jv, ijv
  REAL(DP)  :: jl(ndmx), jlp1(ndmx), aux(ndmx), sum
  !
  DO np=1,ntyp
     msh_bp=kkbeta(np)
     IF (tvanp(np)) THEN
        DO iv =1, nbeta(np)
           DO jv =iv, nbeta(np)
              ijv = jv * (jv-1) / 2 + iv
              ilmin = iabs(lll(iv,np)-lll(jv,np))
              ilmax = iabs(lll(iv,np)+lll(jv,np))
              !       only need to calculate for for lmin,lmin+2 ...lmax-2,lmax
              DO l = ilmin,ilmax,2
                 DO i =  msh_bp,2,-1
                    IF (rgrid(np)%r(i) .LT. rinner(l+1,np)) GOTO 100
                    aux(i) = qfunc(i,ijv,np)
                 ENDDO
100              CALL setqf(qfcoef(1,l+1,iv,jv,np),aux(1),rgrid(np)%r &
                      ,nqf(np),l,i)

                 IF (idbes .EQ. 1) THEN
                    !
                    CALL sph_dbes( msh_bp, rgrid(np)%r, ql, l, jl )
                    !
                    ! ... this is the old call
                    !
                    ! CALL dbess( ql, l+1, msh_bp, r(1,np), jl )
                    !
                 ELSE
                    !
                    CALL sph_bes( msh_bp, rgrid(np)%r, ql, l, jl )
                    !
                    ! ... this is the old call
                    !
                    ! CALL bess( ql, l+1, msh_bp, r(1,np), jl )
                    !
                 ENDIF

                 ! jl is the Bessel function (or its derivative) calculated at ql
                 ! now integrate qfunc*jl*r^2 = Bessel transform of qfunc

                 DO i=1, msh_bp
                    jlp1(i) = jl(i)*aux(i)
                 ENDDO
                 !                        if (tlog(np)) then
                 CALL radlg1(msh_bp,jlp1,rgrid(np)%rab,sum) 

                 qr_k(iv,jv,l+1,np) = sum*fpi/omega
                 qr_k(jv,iv,l+1,np) = qr_k(iv,jv,l+1,np)

              END DO
           END DO
        ENDDO
     ENDIF
  ENDDO
  !
  RETURN
END SUBROUTINE calc_btq
