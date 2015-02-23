!
! Copyright (C) 2002 FPMD group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

!------------------------------------------------------------------------------!
  MODULE sic_module
!------------------------------------------------------------------------------!
!
!  The versions after 3.0 contain also the self-interaction-correction method
!  has proposed by Mauri et al. (PRB 2005), taking also into account the 'comment'
!  proposed by Sprik et al. (ICR 2005).
!  Thus, we introduce the parameters sic_alpha and sic_epsilon to correct the
!  exchange-correlation and the electronic hartree potentials, respectively.
!  They are two empirical parameters, thus to remain in a ab-initio
!  set them equal to 1.0_DP.
!  Sprik et al. showed that, in same cases, i.e. OH radical, it should be better
!  to under estimate the correction to ex-ch, since in same way the exch already
!  corrects the electronic hartree part.
! HOW AND WHEN USE THE SIC::
! Fran's personal considerations:
!     the SIC is a way to correct the self-interaction WHEN
!     ONE and only ONE e- lives in an unpaired electronic level
!     we have choosen for it the spin up
!     Remember to select nspin == 2 and nelup = neldw + 1
!     the other e- are fictitious calculate in a LSD approach:
!     infact, even if the paired e- feel a different potential (for spin up and spin dw)
!     we constrain them to have the same force, and the same eigenvalues, the same eigenstates
!     When you applied this SIC scheme to a molecule or to an atom, which are neutral,
!     remember that you have to consider another correction to the energy level as proposed
!     by Landau: infact if you start from a neutral system and subtract the self-intereaction
!     the unpaired e- feels a charge system. Thus remeber a correction term ~2.317(Madelung)/2L_box


      USE kinds, ONLY: DP
!
      IMPLICIT NONE
      SAVE

      INTEGER :: self_interaction = 0 
      REAL(DP) :: sic_epsilon = 0.0_DP
      REAL(DP) :: sic_alpha = 0.0_DP

!------------------------------------------------------------------------------!
  CONTAINS
!------------------------------------------------------------------------------!

    SUBROUTINE sic_initval( nat_ ,  sic_ ,  sic_epsilon_ , sic_alpha_ )

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: nat_
      CHARACTER(LEN=80), INTENT(IN) :: sic_
      REAL(DP), INTENT(IN) :: sic_epsilon_
      REAL(DP), INTENT(IN) :: sic_alpha_

      select case ( TRIM( sic_ ) )
        case ( 'sic_mac' )
          self_interaction = 2
        case default
          self_interaction = 0
      end select
      sic_epsilon     = sic_epsilon_
      sic_alpha       = sic_alpha_
      !
      RETURN
    END SUBROUTINE sic_initval
    
!------------------------------------------------------------------------------!

  SUBROUTINE sic_info( )

    USE io_global, ONLY: stdout
    IMPLICIT NONE

    !
    ! prints the type of USIC we will do :
    !

    IF( self_interaction == 0 ) THEN
      RETURN
    END IF

        WRITE(stdout, 591)
        WRITE(stdout, 592) self_interaction
        WRITE(stdout, 593)
        !!select case (self_interaction)

       IF ( self_interaction /= 0 ) THEN

          write(stdout,*) &
            '  Unpaired-electron self-interaction correction by Mauri', self_interaction
          write(stdout,*) &
            '  E_USIC_EHTE = U_hartree[rho_up + rho_dw]- sic_espilon * U_hartree[rho_up-rhp_down]'
          write(stdout,*) &
            '  E_USIC_XC   = E_xc[rho_up,rho_dw] - sic_alpha( E_xc[rho_up,rho_dw] + E_xc[rho_dw, rho_dw]) '

       END IF !!select

  591 FORMAT(   3X,' ')
  592 FORMAT(   3X,'Introducing a Mauri Avezac Calandra Self_Interaction Correction: ', I3)
  593 FORMAT(   3X,'----------------------------------------')

    RETURN
  END SUBROUTINE sic_info


!------------------------------------------------------------------------------!
  END MODULE sic_module
!------------------------------------------------------------------------------!
