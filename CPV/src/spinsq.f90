!-----------------------------------------------------------------------
   SUBROUTINE spinsq( c, bec, rhor )
!-----------------------------------------------------------------------
      !! Estimate of \( \langle S^2 \rangle=s(s+1) \) in two different 
      !! ways:
      !
      !! * using as many-body wavefunction a single Slater determinant
      !!   constructed with Kohn-Sham orbitals:
      !
      !!   \[ \langle S^2 \rangle = (N_{up}-N_{dw})/2 \cdot (N_{up}-N_{dw})/2+1)  
      !!       + N_{dw} - \sum_{up}\sum_{dw} \langle \psi_{up} | \psi_{dw} \rangle \]
      !
      !!   where \(N_{up}\), \(N_{dw}\) are the number of up and down states and 
      !!   the sum is over occupied states. Not suitable for fractionary occupancy.
      !!   In the ultrasoft scheme (c is the smooth part of \psi): 
      !
      !!   \[ \langle \psi_{up} | \psi_{dw} \rangle = \sum_G c^*_{up}(G) c_{dw}(G) +
      !!     \int Q_{ij} \langle c_{up}|\beta_i\rangle\langle\beta_j|c_{dw}\rangle \]
      !
      !!   This is the usual formula, unsuitable for fractionary occupancy;
      !
      !! * using the "LSD model" of Wang, Becke, Smith, JCP 102, 3477 (1995):
      !!
      !!   \[ \langle S^2 \rangle = (N_{up}-N_{dw})/2 \cdot (N_{up}-N_{dw})/2+1) + N_{dw} -
      !!           \int \text{max}(\rho_{up}(r),\rho_{dw}(r)) dr \]
      !
      !! Requires on input: \(c=\psi\), \(\text{bec}=\langle c|\beta\rangle\), 
      !! \(\rho_{up}(r)\), \(\rho_{dw}(r)\).  
      !! Assumes real \(\psi\), with only half G vectors.
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
      USE uspp, ONLY: nkb, nkbus, qq_nt, ofsbeta 
      USE uspp_param, ONLY: nh, upf
      USE ions_base, ONLY: na, nat, ityp
!
      IMPLICIT NONE
! input
      REAL(dp) :: bec(nkb,n), rhor(dfftp%nnr,nspin)
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
            DO ia=1,nat
               is=ityp(ia)
               IF( upf(is)%tvanp ) THEN
                  DO iv=1,nh(is)
                     DO jv=1,nh(is)
                        IF(ABS(qq_nt(iv,jv,is)).GT.1.e-5) THEN 
                           inl = ofsbeta(ia) + iv
                           jnl = ofsbeta(ia) + jv
                           overlap(i,j) = overlap(i,j) + qq_nt(iv,jv,is)*bec(inl,i)*bec(jnl,jj)
                        ENDIF
                     END DO
                  END DO
               END IF
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


