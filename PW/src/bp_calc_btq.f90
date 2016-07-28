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
  USE constants, ONLY: fpi
  USE uspp_param, ONLY: upf, nbetam, lmaxq
  !
  IMPLICIT NONE
  !
  REAL(DP)  :: ql, qr_k(nbetam,nbetam,lmaxq,ntyp)
  INTEGER :: idbes
  !
  INTEGER :: i, np, l, ilmin, ilmax, iv, jv, ijv, ilast
  REAL(DP) :: qrk
  REAL(DP), ALLOCATABLE :: jl(:), aux(:)
  !
  DO np=1,ntyp
     !
     IF ( upf(np)%tvanp ) THEN
        !
        ALLOCATE ( jl(upf(np)%kkbeta), aux(upf(np)%kkbeta) ) 
        DO iv =1, upf(np)%nbeta
           DO jv =iv, upf(np)%nbeta
              ijv = jv * (jv-1) / 2 + iv
              ilmin = abs ( upf(np)%lll(iv) - upf(np)%lll(jv) )
              ilmax =       upf(np)%lll(iv) + upf(np)%lll(jv)
              !       only need to calculate for l=lmin,lmin+2 ...lmax-2,lmax
              DO l = ilmin,ilmax,2
                 aux(:) = 0.0_DP
                 aux(1:upf(np)%kkbeta) =  upf(np)%qfuncl(1:upf(np)%kkbeta,ijv,l)
                 IF (idbes == 1) THEN
                    !
                    CALL sph_dbes( upf(np)%kkbeta, rgrid(np)%r, ql, l, jl )
                    !
                 ELSE
                    !
                    CALL sph_bes( upf(np)%kkbeta, rgrid(np)%r, ql, l, jl )
                    !
                 ENDIF

                 ! jl is the Bessel function (or its derivative) calculated at ql
                 ! now integrate qfunc*jl*r^2 = Bessel transform of qfunc

                 DO i=1, upf(np)%kkbeta
                    aux(i) = jl(i)*aux(i)
                 ENDDO
                 !                        if (tlog(np)) then
                 CALL simpson(upf(np)%kkbeta,aux,rgrid(np)%rab,qrk) 

                 qr_k(iv,jv,l+1,np) = qrk*fpi/omega
                 qr_k(jv,iv,l+1,np) = qr_k(iv,jv,l+1,np)

              END DO
           END DO
        ENDDO
        DEALLOCATE ( aux, jl )
     ENDIF
  ENDDO
  !
  RETURN
END SUBROUTINE calc_btq
