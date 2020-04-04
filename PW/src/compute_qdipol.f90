!
! Copyright (C) 2001-2004 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
SUBROUTINE compute_qdipol( dpqq )
  !! This routine computes the term dpqq, i.e. the dipole moment of the
  !! augmentation charge. The output is given on cartesian coordinates.
  !
  USE kinds,       ONLY: DP
  USE constants,   ONLY: fpi
  USE atom,        ONLY: rgrid
  USE ions_base,   ONLY: ntyp => nsp
  USE uspp,        ONLY: nhtol, nhtolm, indv, nlx, ap
  USE uspp_param,  ONLY: upf, nbetam, nh, nhm
  !
  IMPLICIT NONE
  !
  REAL(DP) :: dpqq(nhm, nhm, 3, ntyp)
  !! Dipole moment of augmentation charge
  REAL(DP), ALLOCATABLE :: qrad2(:,:,:), qtot(:,:,:), aux(:)
  REAL(DP) :: fact
  INTEGER :: nt, l, ir, nb, mb, ijv, ilast, ipol, ih, ivl, jh, jvl, lp, ndm
  !
  CALL start_clock( 'cmpt_qdipol' )
  !
  ndm = MAXVAL( upf(1:ntyp)%kkbeta )
  ALLOCATE( qrad2( nbetam , nbetam, ntyp) )
  ALLOCATE( aux( ndm)                     )
  ALLOCATE( qtot( ndm, nbetam, nbetam)    )

  qrad2(:,:,:)=0.d0
  dpqq=0.d0

  DO nt = 1, ntyp
     IF ( upf(nt)%tvanp ) THEN
        l=1
        !
        !   Only l=1 terms enter in the dipole of Q
        !
        DO nb = 1, upf(nt)%nbeta
           DO mb = nb, upf(nt)%nbeta
              ijv = mb * (mb-1) /2 + nb
              IF ( ( l >= ABS(upf(nt)%lll(nb) - upf(nt)%lll(mb)) ) .AND. &
                   ( l <=     upf(nt)%lll(nb) + upf(nt)%lll(mb)  ) .AND. &
                   (MOD(l+upf(nt)%lll(nb)+upf(nt)%lll(mb), 2) == 0) ) THEN
                 qtot(1:upf(nt)%kkbeta,nb,mb) = upf(nt)%qfuncl(1:upf(nt)%kkbeta,ijv,l)
              ENDIF
           ENDDO
        ENDDO
        do nb=1, upf(nt)%nbeta
           !
           !    the Q are symmetric with respect to indices
           !
           do mb=nb, upf(nt)%nbeta
              if ( ( l >= ABS(upf(nt)%lll(nb) - upf(nt)%lll(mb)) ) .AND. &
                   ( l <=     upf(nt)%lll(nb) + upf(nt)%lll(mb)  ) .AND. &
                   (MOD(l+upf(nt)%lll(nb)+upf(nt)%lll(mb), 2) == 0) ) THEN
                 do ir = 1, upf(nt)%kkbeta
                    aux(ir)=rgrid(nt)%r(ir)*qtot(ir, nb, mb)
                 ENDDO
                 call simpson ( upf(nt)%kkbeta, aux, rgrid(nt)%rab, &
                                qrad2(nb,mb,nt) )
              ENDIF
           ENDDO
        ENDDO
     ENDIF
     ! ntyp
  ENDDO
  
  DO ipol = 1, 3
     fact = -SQRT(fpi/3.d0)
     IF (ipol == 1) lp=3
     IF (ipol == 2) lp=4
     IF (ipol == 3) THEN
        lp=2
        fact=-fact
     ENDIF
     DO nt = 1, ntyp
        IF ( upf(nt)%tvanp ) THEN
           DO ih = 1, nh(nt)
              ivl = nhtolm(ih, nt)
              mb = indv(ih, nt)
              DO jh = ih, nh(nt)
                 jvl = nhtolm(jh, nt)
                 nb = indv(jh,nt)
                 IF (ivl > nlx) CALL errore( 'compute_qdipol',' ivl > nlx', ivl )
                 IF (jvl > nlx) CALL errore( 'compute_qdipol',' jvl > nlx', jvl )
                 IF (nb > nbetam) &
                    CALL errore( 'compute_qdipol',' nb out of bounds', nb )
                 IF (mb > nbetam) &
                    CALL errore( 'compute_qdipol',' mb out of bounds', mb )
                 IF (mb > nb) CALL errore( 'compute_qdipol',' mb > nb', 1 )
                 dpqq(ih,jh,ipol,nt) = fact * ap(lp,ivl,jvl) * qrad2(mb,nb,nt)
                 dpqq(jh,ih,ipol,nt) = dpqq(ih,jh,ipol,nt)
                 ! WRITE( stdout,'(3i5,2f15.9)') ih,jh,ipol,dpqq(ih,jh,ipol,nt)
              ENDDO
           ENDDO
        ENDIF
     ENDDO
  ENDDO
  !
  DEALLOCATE( qtot  )
  DEALLOCATE( aux   )
  DEALLOCATE( qrad2 )
  !
  CALL stop_clock( 'cmpt_qdipol' )
  !
  RETURN
  !
END SUBROUTINE compute_qdipol
