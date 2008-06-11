!
! Copyright (C) 2007-2008 Quantum-ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! original version by I. Dabo and N. Marzari (MIT)
!
! contributions by E. Lamas and S. de Gironcoli (SISSA/DEMOCRITOS)
!
!----------------------------------------------------------------------
      SUBROUTINE setlocalcoul
!----------------------------------------------------------------------
       ! 
      ! ... Calculates the Coulomb potential of the ionic cores 
      ! ... in real space 
      !
      USE kinds,         ONLY : DP
      USE ions_base,     ONLY : nat,                                   &
                                ityp,                                  &
                                tau,                                   &
                                ntyp => nsp,                           &
                                zv
!! Eduardo changed everywhere zp by upf(nt)%zp      USE pseud,         ONLY : zp
      USE cell_base,     ONLY : alat,                                  &
                                at
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! DILDCC
!! commented out some variables not used in the subroutine Eduardo
!!
!      USE atom,          ONLY : numeric,                               &
!      USE atom,          ONLY : msh,                                   &
!                                mesh,                                  &
!                                r,                                     &
!                                rab
      USE atom,          ONLY :  msh,                                    &
                                 rgrid
!! notice change from r(ir,nt) to to rgrid(nt)%r(ir) in all ocurrences of r
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Also changing vloc_ati(ir,nt) for upf(nt)%vloc(ir)
!! Eduardo 
!! DILDCC
!!
!      USE uspp_param,    ONLY : vloc_at
      USE uspp_param, ONLY : upf
      USE ee_mod,        ONLY : vloccoul, mr1, mr2, mr3
!!    USE parser,        ONLY : int_to_char   
      USE constants,     ONLY : e2, fpi
      !
      IMPLICIT NONE
      !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! int_to_char moved (see 5 lines above)
      CHARACTER(LEN=6), EXTERNAL :: int_to_char

      INTEGER                :: nt
      INTEGER                :: ir
      INTEGER                :: na
      INTEGER                :: ir1
      INTEGER                :: ir2
      INTEGER                :: ir3
      INTEGER                :: irr
      !
      REAL( DP )             :: delta1
      REAL( DP )             :: delta2
      REAL( DP )             :: delta3
      REAL( DP )             :: dist
      REAL( DP )             :: adist
      REAL( DP )             :: temp
      !
      REAL( DP ), EXTERNAL   :: erf              
      !
      CHARACTER( LEN = 6 )   :: nts
      !
      INTEGER, EXTERNAL      :: COMPINDEX
      !
      DO nt = 1, ntyp
        nts = INT_TO_CHAR( nt )
        OPEN( 200, file = 'vloc' // TRIM( nts )                        &
              // '.dat', status = 'unknown' )
        WRITE( 200, '(3F30.10)' ) ( rgrid(nt)%r(ir), upf(nt)%vloc(ir),    &
                               - upf(nt)%zp * e2 / rgrid(nt)%r(ir),          &
                                 ir = 1, msh (nt) )
      END DO
      !
      delta1 = alat * at( 1, 1 ) / DBLE( mr1 )
      delta2 = alat * at( 2, 2 ) / DBLE( mr2 )
      delta3 = alat * at( 3, 3 ) / DBLE( mr3 )
      
      !
      vloccoul = 0.D0
      !
      DO na = 1, nat
        nt = ityp( na )
        DO ir1 = 1, mr1
          DO ir2 = 1, mr2
            DO ir3 = 1, mr3
              !IF( ( ir1 == 1 ) .OR. ( ir1 == mr1 ) .OR.                &
              !    ( ir2 == 1 ) .OR. ( ir2 == mr2 ) .OR.                &
              !    ( ir3 == 1 ) .OR. ( ir3 == mr3 ) ) THEN
              irr = COMPINDEX( ir1, ir2, ir3, mr1, mr2, mr3 )
              dist =                                                   &
                ( tau( 1, na ) * alat - DBLE( ir1 - 1 ) * delta1 )** 2 &
              + ( tau( 2, na ) * alat - DBLE( ir2 - 1 ) * delta2 )** 2 &
              + ( tau( 3, na ) * alat - DBLE( ir3 - 1 ) * delta3 )** 2
              dist = SQRT( dist )
              ir = 1
              DO WHILE( ( rgrid(nt)%r(ir+1)  < dist )                   &
                        .AND. ( ir + 1 < msh( nt ) ) ) 
                ir = ir + 1
              END DO
              IF( ir + 1 < msh( nt ) ) THEN
                temp = ( upf(nt)%vloc(ir+1) - upf(nt)%vloc(ir) )   &
                  / ( LOG( rgrid(nt)%r( ir + 1 ) ) - LOG( rgrid(nt)%r(ir) ) )
                temp = upf(nt)%vloc(ir)                            &
                  + temp * ( LOG( dist ) - LOG( rgrid(nt)%r(ir) ) )
              ELSE
                temp = - upf(nt)%zp * e2 / dist
              END IF
              vloccoul( irr ) = vloccoul( irr ) + temp    
              !END IF
            END DO
          END DO
        END DO
      END DO
      !
      RETURN
      !
!----------------------------------------------------------------------
      END SUBROUTINE setlocalcoul
!----------------------------------------------------------------------
