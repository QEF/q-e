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
!--------------------------------------------------------------------
      SUBROUTINE write_ee_summary( )
!--------------------------------------------------------------------
      !
      USE io_global,      ONLY : stdout,                               &
                                 ionode
      USE ee_mod,         ONLY : do_comp,       &
                                 which_compensation,                   &
                                 n_self_interaction, poisson_maxiter,  &
                                 poisson_thr,                          &
                                 n_charge_compensation,                &
                                 mixing_charge_compensation,           &
                                 mr1,                                  &
                                 mr2,                                  &
                                 mr3,                                  &
                                 ecutcoarse,                           &
                                 ncompx,                               &
                                 ncompy,                               &
                                 ncompz,                               &
                                 rhocut,                               &
                                 rhocution,                            &
                                 comp_thr,                             &
                                 n_smoothing,                          &
                                 nlev,                                 &
                                 do_mltgrid,                           &
                                 rhoionmax,                            &
                                 which_smoothing                      
      !
      IMPLICIT NONE
      !

      IF( ionode .AND. do_comp ) THEN
        !
        WRITE( UNIT = stdout,                                          &
               FMT  = '(/,5x, "Electrostatic Correction",             &
                       &/,5x, "========================")' )
        !
        WRITE( UNIT = stdout, FMT = 9026 ) comp_thr
        !
        WRITE( UNIT = stdout, FMT = 9005 ) n_charge_compensation,      &
                                       TRIM(which_compensation),       &
                                       mixing_charge_compensation
        !
        WRITE( stdout, * )
        !
      END IF
      !
      IF( ionode ) THEN
        !
        SELECT CASE( TRIM( which_compensation ) )
        !
        CASE( 'dcc' )
          !
          WRITE( UNIT = stdout,                                       &
                 FMT = '(/,5x, "Poisson Solver",                     &
                        &/,5x, "==============")' )
          !
          IF( mr1 .NE. 0 .AND. mr2 .NE. 0 .AND. mr3 .NE. 0 ) THEN
            WRITE( UNIT = stdout, FMT = 9024 ) mr1, mr2, mr3
          ELSE 
            WRITE( UNIT = stdout, FMT = 9025 ) ecutcoarse
          END IF
          !
          IF ( ionode .AND. do_mltgrid ) THEN
            !
            WRITE( UNIT = stdout, FMT = 9030 ) nlev
            !
          END IF
          !
        END SELECT
        !
        WRITE( stdout, * )
        !
      END IF
      !
9005 FORMAT( '     charge-compensation periodicity   = ',  I24,' '     &
            /'     electrostatic correction          = ',  A24,' '     &
            /'     charge-compensation mixing        = ',  F24.2,' ')
9024 FORMAT( '     coarse Poisson-solver grid        = ',  3I8,' ' )
9025 FORMAT( '     coarse-grid energy cutoff (Ry)    = ',  F24.2,' ')
9026 FORMAT( '     compensation onset threshold      = ',  E24.8,' ')
9030 FORMAT( '     number multigrid levels           = ',  I24,' ' )

!--------------------------------------------------------------------
      END SUBROUTINE write_ee_summary
!--------------------------------------------------------------------
