!-----------------------------------------------------------------------
      SUBROUTINE spinsq (c,bec,rhor)
!-----------------------------------------------------------------------
!
!     estimate of <S^2>=s(s+1) in two different ways.
!     1) using as many-body wavefunction a single Slater determinant
!        constructed with Kohn-Sham orbitals:
!
!        <S^2> = (Nup-Ndw)/2 * (Nup-Ndw)/2+1) + Ndw - 
!                \sum_up\sum_dw < psi_up | psi_dw >
!
!        where Nup, Ndw = number of up and down states, the sum is over 
!        occupied states. Not suitable for fractionary occupancy.
!        In the ultrasoft scheme (c is the smooth part of \psi): 
!
!        < psi_up | psi_dw > = \sum_G c*_up(G) c_dw(G) +
!                              \int Q_ij <c_up|beta_i><beta_j|c_dw>
!
!        This is the usual formula, unsuitable for fractionary occupancy.
!     2) using the "LSD model" of Wang, Becke, Smith, JCP 102, 3477 (1995):
!
!        <S^2> = (Nup-Ndw)/2 * (Nup-Ndw)/2+1) + Ndw -
!                \int max(rhoup(r),rhodw(r)) dr
!
!     Requires on input: c=psi, bec=<c|beta>, rhoup(r), rhodw(r)
!     Assumes real psi, with only half G vectors.
!
      USE kinds, only: dp
      USE electrons_base, ONLY: nx => nbspx, n => nbsp, iupdwn, nupdwn, f, nel, nspin
      USE io_global, ONLY: stdout
      USE mp_global, ONLY: intra_bgrp_comm
      USE mp, ONLY: mp_sum
      USE fft_base, ONLY: dfftp
      USE gvecw, ONLY: ngw
      USE gvect, ONLY: gstart
      USE cell_base, ONLY: omega
      USE uspp, ONLY: nhsa => nkb, nhsavb=>nkbus, qq
      USE uspp_param, ONLY: nvb, ish, nh
      USE ions_base, ONLY: na
!
      IMPLICIT NONE
! input
      REAL(dp) :: bec(nhsa,n), rhor(dfftp%nnr,nspin)
      COMPLEX(dp):: c(ngw,nx)
! local variables
      INTEGER :: nup, ndw, ir, i, j, jj, ig, ia, is, iv, jv, inl, jnl
      REAL(dp) :: spin0, spin1, spin2, fup, fdw
      REAL(dp), ALLOCATABLE:: overlap(:,:), temp(:)
      LOGICAL :: frac
!
!
      IF (nspin.EQ.1) RETURN
!
! find spin-up and spin-down states
!
      fup = 0.0d0
      DO i=iupdwn(1),nupdwn(1)
         fup = fup + f(i)
      END DO
      nup = NINT(fup)
      ndw = nel(1)+nel(2) - nup
!
! paranoid checks
!
      frac= ABS(fup-nup).GT.1.0d-6
      fup = 0.0d0
      DO i=1,nup
         fup = fup + f(i)
      END DO
      frac=frac.OR.ABS(fup-nup).GT.1.0d-6
      fdw = 0.0d0
      DO j=iupdwn(2),iupdwn(2)-1+ndw
         fdw = fdw + f(j)
      END DO
      frac=frac.OR.ABS(fdw-ndw).GT.1.0d-6
!
      spin0 = ABS(fup-fdw)/2.d0 * ( ABS(fup-fdw)/2.d0 + 1.d0 ) + fdw
!
!     Becke's formula for spin polarization
!
      spin1 = 0.0d0
      DO ir=1,dfftp%nnr
         spin1 = spin1 - MIN(rhor(ir,1),rhor(ir,2))
      END DO
      CALL mp_sum( spin1, intra_bgrp_comm )
      spin1 = spin0 + omega/(dfftp%nr1*dfftp%nr2*dfftp%nr3)*spin1
      IF (frac) THEN
         WRITE( stdout,'(/" Spin contamination: s(s+1)=",f5.2," (Becke) ",&
     &                             f5.2," (expected)")')              &
     &          spin1, ABS(fup-fdw)/2.d0*(ABS(fup-fdw)/2.d0+1.d0)
         RETURN
      END IF
!
!     Slater formula, smooth contribution to  < psi_up | psi_dw >
!
      ALLOCATE (overlap(nup,ndw))
      ALLOCATE (temp(ngw))
      DO j=1,ndw
         jj=j+iupdwn(2)-1
         DO i=1,nup
            overlap(i,j)=0.d0
            DO ig=1,ngw
               temp(ig)=DBLE(CONJG(c(ig,i))*c(ig,jj))
            END DO
            overlap(i,j) = 2.d0*SUM(temp)
            IF (gstart == 2) overlap(i,j) = overlap(i,j) - temp(1)
         END DO
      END DO
      DEALLOCATE (temp)
      CALL mp_sum( overlap, intra_bgrp_comm )
      DO j=1,ndw
         jj=j+iupdwn(2)-1
         DO i=1,nup
!
!     vanderbilt contribution to  < psi_up | psi_dw >
!
            DO is=1,nvb
               DO iv=1,nh(is)
                  DO jv=1,nh(is)
                     IF(ABS(qq(iv,jv,is)).GT.1.e-5) THEN 
                        DO ia=1,na(is)
                           inl=ish(is)+(iv-1)*na(is)+ia
                           jnl=ish(is)+(jv-1)*na(is)+ia
                           overlap(i,j) = overlap(i,j) +                &
     &                          qq(iv,jv,is)*bec(inl,i)*bec(jnl,jj)
                        END DO
                     ENDIF
                  END DO
               END DO
            END DO
         END DO
      END DO
!
      spin2 = spin0
      DO j=1,ndw
         DO i=1,nup
            spin2 = spin2 - overlap(i,j)**2
         END DO
      END DO
!
      DEALLOCATE (overlap)
!
      WRITE( stdout,'(/" Spin contamination: s(s+1)=",f5.2," (Slater) ",  &
     &          f5.2," (Becke) ",f5.2," (expected)")')              &
     &     spin2,spin1, ABS(fup-fdw)/2.d0*(ABS(fup-fdw)/2.d0+1.d0)
!
      RETURN
      END SUBROUTINE spinsq


