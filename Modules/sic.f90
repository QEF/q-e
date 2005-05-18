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

      USE kinds, ONLY: dbl
      USE parameters, ONLY: natx
!
      IMPLICIT NONE
      SAVE

      INTEGER :: ind_localisation(natx) = 0   ! true if we want to know the localization arount the atom
      INTEGER :: nat_localisation = 0 
      LOGICAL :: print_localisation = .FALSE. ! Calculates hartree energy around specified atoms
      INTEGER :: self_interaction = 0 
      REAL(dbl) :: si_epsilon = 0.0d0
      REAL(dbl) :: rad_localisation = 0.0d0
      REAL(dbl), ALLOCATABLE :: pos_localisation(:,:)

!------------------------------------------------------------------------------!
  CONTAINS
!------------------------------------------------------------------------------!

    SUBROUTINE sic_initval( nat_ , id_loc_ , sic_ ,  sic_epsilon_ , sic_rloc_ )

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: nat_
      INTEGER, INTENT(IN) :: id_loc_ (:)
      CHARACTER(LEN=*), INTENT(IN) :: sic_
      REAL(dbl), INTENT(IN) :: sic_epsilon_
      REAL(dbl), INTENT(IN) :: sic_rloc_

      select case ( TRIM( sic_ ) )
        case ( 'sic_pz' ) 
          self_interaction = 1
        case ( 'sic_mac' )
          self_interaction = 2
        case ( 'only_sich' )
          self_interaction = 3
        case ( 'only_sicxc_pz' )
          self_interaction = -1
        case ( 'only_sicxc_mac' )
          self_interaction = -2
        case default
          self_interaction = 0
      end select
      si_epsilon       = sic_epsilon_
      rad_localisation = sic_rloc_
      ! counting the atoms around which i want to calculate the charge localization
      ind_localisation( 1 : nat_ ) = id_loc_ ( 1 : nat_ )
      nat_localisation = COUNT( ind_localisation > 0 ) 
      IF( ALLOCATED( pos_localisation ) ) DEALLOCATE( pos_localisation )
      ALLOCATE( pos_localisation( 4, MAX( nat_localisation, 1 ) ) )
      !
      IF( nat_localisation > 0 ) print_localisation = .TRUE.
      !
      RETURN
    END SUBROUTINE sic_initval
    
!------------------------------------------------------------------------------!

    SUBROUTINE deallocate_sic()
      IMPLICIT NONE
      IF( ALLOCATED( pos_localisation ) ) DEALLOCATE( pos_localisation )
      RETURN
    END SUBROUTINE deallocate_sic

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
        select case (self_interaction)
        case (1)
          write(stdout,*) &
            '  Unpaired-electron self-interaction correction of type ', self_interaction
          write(stdout,*) &
            '  E_USIC_PZ = - U_hartree[rho_up-rhp_down] - E_excor[rho_up-rho_down,0] '
        case (2)
          write(stdout,*) &
            '  Unpaired-electron self-interaction correction of type ', self_interaction
          write(stdout,*) &
            '  E_USIC_MAC = - U_hartree[rho_up-rhp_down] - E_excor[rho_up,rho_down] + E[rho_down, rho_down] '
        case (3)
          write(stdout,*) &
            '  Unpaired-electron self-interaction correction of type ', self_interaction
          write(stdout,*) &
            '  E_USIC_basis = - U_hartree[rho_up-rhp_down] '
        case (-1)
          write(stdout,*) &
            '  Unpaired-electron self-interaction correction of type ', self_interaction
          write(stdout,*) &
            '  E_USIC_nohartree_PZ =  - E_excor[rho_up-rho_down,0] '
        case (-2)
          write(stdout,*) &
            '  Unpaired-electron self-interaction correction of type ', self_interaction
          write(stdout,*) &
            '  E_USIC_nohartree_MAC =  - E_excor[rho_up,rho_down] + E[rho_down, rho_down] '
        case default
          write(stdout,*) &
            '  Unpaired-electron self-interaction correction of type ', self_interaction
          write(stdout,*) &
            '  No unpaired-electron self-interaction correction '
        end select
  591 FORMAT(   3X,'')
  592 FORMAT(   3X,'Introducing a Self_Interaction Correction case: ', I3)
  593 FORMAT(   3X,'----------------------------------------')

    RETURN
  END SUBROUTINE sic_info


!------------------------------------------------------------------------------!
  END MODULE sic_module
!------------------------------------------------------------------------------!
