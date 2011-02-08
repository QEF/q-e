!
! Copyright (C) 2002-2005 FPMD-CPV groups
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!


  MODULE polarization

        USE kinds
        USE berry_phase, only: indi_l, sour_indi, dest_indi, n_indi_rcv, n_indi_snd, icntix

        IMPLICIT NONE

        SAVE

        PRIVATE

        !  variables used for the dipole moment

        REAL(DP) :: p0( 3 ), p( 3 ), pdipole( 3 ), pdipolt( 3 ), pdipole0( 3 )
        REAL(DP) :: cost1, cost2, cost3, fac
        REAL(DP) :: d1old, d2old, d3old
        LOGICAL  :: first = .TRUE.

        PUBLIC :: deallocate_polarization, ddipole
        PUBLIC :: print_dipole


 CONTAINS


   SUBROUTINE deallocate_polarization
      use berry_phase, only: berry_closeup
      call berry_closeup()
      return
   END SUBROUTINE deallocate_polarization



   SUBROUTINE print_dipole( uni, tfile, nfi, tps )
      !
      USE io_global,  ONLY : stdout, ionode
      !
      INTEGER,  INTENT(IN) :: uni
      LOGICAL,  INTENT(IN) :: tfile
      INTEGER,  INTENT(IN) :: nfi
      REAL(DP), INTENT(IN) :: tps

      INTEGER :: i

      WRITE( stdout, 19 )
      WRITE( stdout, 20 ) 'P0_e',  (pdipole0(i),i=1,3)
      WRITE( stdout, 20 ) 'P_e',   (pdipole(i),i=1,3)
      WRITE( stdout, 20 ) 'P0_I',  (p0(i),i=1,3)
      WRITE( stdout, 20 ) 'P_I',   (p(i),i=1,3)
      WRITE( stdout, 20 ) 'P_tot', (pdipolt(i),i=1,3)
      !
      IF (tfile) THEN
         WRITE( uni, 30 ) nfi, tps
         WRITE( uni, 20 ) 'P_e',   (pdipole(i),i=1,3)
         WRITE( uni, 20 ) 'P_I',   (p(i),i=1,3)
         WRITE( uni, 20 ) 'P_tot', (pdipolt(i),i=1,3)
      END IF

19    FORMAT(/,3X,'Dipole moment (AU)')
20    FORMAT(3X,A10,3(F18.8,2X))
30    FORMAT(I7,1X,F11.8)

      RETURN
   END SUBROUTINE print_dipole



   SUBROUTINE ddipole( istep, c2, ngwx, taus, tfor, ngw, n, ht )

      USE mp,          ONLY: mp_sum
      USE constants,   ONLY: pi
      USE cell_base,   ONLY: tpiba, alat, s_to_r
      USE ions_base,   ONLY: zv, nat, nsp, na
      USE mp_global,   ONLY: me_bgrp, nproc_bgrp, intra_bgrp_comm
      USE mp_wave,     ONLY: pwscatter

      IMPLICIT NONE 

      COMPLEX(DP) :: ZDOTU, zdotc
      EXTERNAL ZDOTU, zdotc
      !
      ! ... ARGUMENTS
      !
      INTEGER,  INTENT(IN) :: istep
      INTEGER,  INTENT(IN) :: ngw, n, ngwx
      REAL(DP), INTENT(IN) :: ht( 3, 3 )  !  Transpose of matrix h
      REAL(DP), INTENT(IN) :: taus( 3, nat )  !  Transpose of matrix h
      LOGICAL,  INTENT(IN) :: tfor
      COMPLEX(DP)          :: c2( ngwx, n )

      !
      ! ... LOCALS
      !
      REAL(DP)    :: taup( 3, nat )
      REAL(DP)    :: d1, d2, d3
      REAL(DP)    :: rb1, rb2, rb3
      REAL(DP)    :: rb1m1, rb2m1, rb3m1
      REAL(DP)    :: rdummy
      REAL(DP)    :: bg(3,3), bgm1(3,3), b1( 3 ), b2( 3 ), b3( 3 )
      COMPLEX(DP) :: dumm( n, n ), det, aux( 2*n ), ptemp( ngwx )
      COMPLEX(DP) :: detc( 2 ), ztmp
      INTEGER     :: ipiv( n ), info
      REAL(DP)    :: omega
      REAL(DP)    :: htm1( 3, 3 ), h( 3, 3 )
      INTEGER     :: i, j, is, in2, in1, me, isa

      !
      ! ... Subroutine body
      !

      me = me_bgrp + 1

      h = TRANSPOSE( ht )

      CALL invmat( 3, ht, htm1, omega )

      do i=1,3
        b1(i) = alat * htm1(i,1)
        b2(i) = alat * htm1(i,2)
        b3(i) = alat * htm1(i,3)
      enddo

      CALL s_to_r( taus, taup, na, nsp, h )

      IF(FIRST) THEN
        FAC=2.D0
        RB1=B1(1)*B1(1) + B1(2)*B1(2) + B1(3)*B1(3)
        RB2=B2(1)*B2(1) + B2(2)*B2(2) + B2(3)*B2(3)
        RB3=B3(1)*B3(1) + B3(2)*B3(2) + B3(3)*B3(3)

        RB1M1=1.d0/SQRT(RB1)
        RB2M1=1.d0/SQRT(RB2)
        RB3M1=1.d0/SQRT(RB3)
        COST1=FAC/omega/TPIBA*RB1M1
        COST2=FAC/omega/TPIBA*RB2M1
        COST3=FAC/omega/TPIBA*RB3M1
        DO I=1,9
          BG(I,1)=0.D0                           
        ENDDO
        CALL daxpy(3,RB1M1,B1,1,BG(1,1),1)
        CALL daxpy(3,RB2M1,B2,1,BG(1,2),1)
        CALL daxpy(3,RB3M1,B3,1,BG(1,3),1)
        CALL invmat (3, BG, BGM1, rdummy)

        !
        ! t=0 initial ionic polarization, only if the atoms move.
        !
        IF( tfor ) THEN
          DO J = 1, 3
            P0(J) = 0.D0
            isa = 0
            DO IS = 1, nsp
              DO I = 1, na(is)
                isa = isa + 1
                P0(J) = P0(J) + ZV(is) * TAUP( j, isa )
              ENDDO
            ENDDO
            P0(J) = P0(J) / omega
          ENDDO
        ENDIF
!
      ENDIF
      !
      !..ionic contribution
      !
      DO J = 1, 3
        P(J) = 0.D0
        isa = 0
        DO IS = 1, nsp
          DO I = 1, na(is)
            isa = isa + 1
            P(J) = P(J) + ZV(is) * TAUP(J,isa)
          ENDDO
        ENDDO
        P(J) = P(J) / omega
      ENDDO
!
!..set vectors
!.
!..P(1)  Polarizability along x
!
      dumm = 0.0d0
      DO IN2 = 1, N
        call pwscatter( C2(:,in2), PTEMP, ngw, indi_l(:,1), sour_indi(:,1), &
          dest_indi(:,1), n_indi_rcv(1), n_indi_snd(1), icntix(1), me_bgrp, nproc_bgrp, intra_bgrp_comm )
        DO IN1 = IN2, N
          ztmp = zdotc( NGW, C2(1,IN1), 1, PTEMP(1), 1 )
          call mp_sum( ztmp, intra_bgrp_comm )
          DUMM(IN1,IN2)=ztmp
        ENDDO
        call pwscatter( C2(:,in2), PTEMP, ngw, indi_l(:,3), sour_indi(:,3), &
          dest_indi(:,3), n_indi_rcv(3), n_indi_snd(3), icntix(3), me_bgrp, nproc_bgrp, intra_bgrp_comm )
        DO IN1=IN2,N
          ztmp = ZDOTU( NGW, C2(1,IN1), 1, PTEMP(1), 1 )
          call mp_sum( ztmp, intra_bgrp_comm )
          DUMM(IN1,IN2)=DUMM(IN1,IN2)+ztmp
        ENDDO
        call pwscatter( C2(:,in2), PTEMP, ngw, indi_l(:,2), sour_indi(:,2), &
          dest_indi(:,2), n_indi_rcv(2), n_indi_snd(2), icntix(2), me_bgrp, nproc_bgrp, intra_bgrp_comm )
        DO IN1=IN2,N
          ztmp = zdotc(NGW,PTEMP(1),1,C2(1,IN1),1)
          call mp_sum( ztmp, intra_bgrp_comm )
          DUMM(IN1,IN2)=DUMM(IN1,IN2) + ztmp
        ENDDO
        DO IN1=1,IN2-1
          DUMM(IN1,IN2)=DUMM(IN2,IN1)
        ENDDO
      ENDDO

!
! Compute determinant and then log(det) for P(1)
!
      CALL ZGEFA(DUMM,n,N,IPIV,INFO)
      CALL ZGEDI(DUMM,n,N,IPIV,DETC,AUX,10)
      DET=DETC(1)*10.D0**DETC(2)
      D1= ATAN2 (AIMAG(DET),DBLE(DET))
      IF(.NOT.FIRST) THEN
        IF(ABS(D1-D1OLD).GT.PI) THEN
          D1 = D1 - SIGN(2*PI,D1-D1OLD)
        END IF
      END IF
      D1OLD = D1

!
!..P(2)
      dumm = 0.0d0
      DO IN2=1,N
        call pwscatter( C2(:,in2), PTEMP, ngw, indi_l(:,4), sour_indi(:,4), &
          dest_indi(:,4), n_indi_rcv(4), n_indi_snd(4), icntix(4), me_bgrp, nproc_bgrp, intra_bgrp_comm )
!. contiene il termine ig=0
        DO IN1=IN2,N
          ztmp = zdotc(NGW,C2(1,IN1),1,PTEMP(1),1)
          call mp_sum( ztmp, intra_bgrp_comm )
          DUMM(IN1,IN2)=ztmp
        ENDDO
        call pwscatter( C2(:,in2), PTEMP, ngw, indi_l(:,6), sour_indi(:,6), &
          dest_indi(:,6), n_indi_rcv(6), n_indi_snd(6), icntix(6), me_bgrp, nproc_bgrp, intra_bgrp_comm )
        DO IN1=IN2,N
          ztmp = ZDOTU(NGW,C2(1,IN1),1,PTEMP(1),1)
          call mp_sum( ztmp, intra_bgrp_comm )
          DUMM(IN1,IN2)=DUMM(IN1,IN2) + ztmp
        ENDDO
        call pwscatter( C2(:,in2), PTEMP, ngw, indi_l(:,5), sour_indi(:,5), &
          dest_indi(:,5), n_indi_rcv(5), n_indi_snd(5), icntix(5), me_bgrp, nproc_bgrp, intra_bgrp_comm )
        DO IN1=IN2,N
          ztmp = zdotc(NGW,PTEMP(1),1,C2(1,IN1),1)
          call mp_sum( ztmp, intra_bgrp_comm ) 
          DUMM(IN1,IN2)=DUMM(IN1,IN2) + ztmp
        ENDDO
! simmetrizzo
        DO IN1=1,IN2-1
          DUMM(IN1,IN2)=DUMM(IN2,IN1)
        ENDDO
      ENDDO
!
! Compute determinant and then log(det) for P(2)
!
      CALL ZGEFA(DUMM,n,N,IPIV,INFO)
      CALL ZGEDI(DUMM,n,N,IPIV,DETC,AUX,10)
      DET=DETC(1)*10.D0**DETC(2)
      D2= ATAN2 (AIMAG(DET),DBLE(DET))
      IF(.NOT.FIRST) THEN
        IF(ABS(D2-D2OLD).GT.PI) THEN
          D2 = D2 - SIGN(2*PI,D2-D2OLD)
        END IF
      END IF
      D2OLD = D2
!
!..P(3)
!
      dumm = 0.0d0
      DO IN2=1,N
        call pwscatter( C2(:,in2), PTEMP, ngw, indi_l(:,7), sour_indi(:,7), &
          dest_indi(:,7), n_indi_rcv(7), n_indi_snd(7), icntix(7), me_bgrp, nproc_bgrp, intra_bgrp_comm )
!. contiene il termine ig=0
        DO IN1=IN2,N
          ztmp = zdotc(NGW,C2(1,IN1),1,PTEMP(1),1)
          call mp_sum( ztmp, intra_bgrp_comm ) 
          DUMM(IN1,IN2)=ztmp
        ENDDO
        call pwscatter( C2(:,in2), PTEMP, ngw, indi_l(:,8), sour_indi(:,8), &
          dest_indi(:,8), n_indi_rcv(8), n_indi_snd(8), icntix(8), me_bgrp, nproc_bgrp, intra_bgrp_comm )
        DO IN1=IN2,N
          ztmp = zdotc(NGW,PTEMP(1),1,C2(1,IN1),1)
          call mp_sum( ztmp, intra_bgrp_comm )
          DUMM(IN1,IN2)=DUMM(IN1,IN2)+ztmp
        ENDDO
! simmetrizzo
        DO IN1=1,IN2-1
          DUMM(IN1,IN2)=DUMM(IN2,IN1)
        ENDDO
      ENDDO
!
! Compute determinant and then log(det) for P(3)
!
      CALL ZGEFA(DUMM,n,N,IPIV,INFO)
      CALL ZGEDI(DUMM,n,N,IPIV,DETC,AUX,10)
      DET=DETC(1)*10.D0**DETC(2)
      D3= ATAN2 (AIMAG(DET),DBLE(DET))
      IF(.NOT.FIRST) THEN
        IF(ABS(D3-D3OLD).GT.PI) THEN
          D3 = D3 - SIGN(2*PI,D3-D3OLD)
        END IF
      END IF
      D3OLD = D3

!
! pdipole has the polarization due to the electronic component,
! p has the ionic component, and pdipolt the total polarization.
!
      DO I=1,3
        PDIPOLE(I) = D1*COST1*BGM1(1,I) + D2*COST2*BGM1(2,I) + D3*COST3*BGM1(3,I) 
        PDIPOLT(I) = PDIPOLE(I) + ( P(I) - P0(I) )
      ENDDO
      IF(FIRST.AND.TFOR) THEN
        PDIPOLE0 = PDIPOLE
      ENDIF
!
      FIRST=.false.
!
  100 FORMAT(3F10.5,3(2X,F10.5))
  20  FORMAT(6X,3(F18.10,2X))
  10  FORMAT(6X,I6)
      RETURN
      END subroutine ddipole


   END MODULE POLARIZATION
