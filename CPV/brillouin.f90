!
! Copyright (C) 2002-2005 FPMD-CPV groups
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

!------------------------------------------------------------------------------!
! 
    MODULE brillouin
! 
!------------------------------------------------------------------------------!
      USE kinds, ONLY : DP
! ... 
!
      IMPLICIT NONE
      SAVE
! 
      PRIVATE
!
! ... CP2K Type ... 
      TYPE kpoints
        CHARACTER (len=20) :: scheme
        LOGICAL :: gamma_only
        INTEGER :: nk1, nk2, nk3
        INTEGER :: k1, k2, k3
        REAL (DP) :: shift(3)
        LOGICAL :: symmetry
        INTEGER :: wfn_type
        INTEGER :: nkpt
        REAL (DP), DIMENSION (:), POINTER :: weight
        REAL (DP), DIMENSION (:,:), POINTER :: xk
      END TYPE kpoints
!------------------------------------------------------------------------------!

      TYPE (kpoints) :: kp
      REAL (DP), ALLOCATABLE , TARGET :: weight(:)
      REAL (DP), ALLOCATABLE , TARGET :: xk(:,:)


      PUBLIC :: kpoints, kpoint_info, kpoint_setup, kp
      PUBLIC :: get_kpoints_number
! 
!
    CONTAINS
!
! CP2K input section
!!>----------------------------------------------------------------------------!
!!  SECTION: &kpoint... &end                                                   !
!!                                                                             !
!!  scheme           [Gamma, Monkhorst-Pack, MacDonald, General]               !
!!    { nx ny nz }                                                             !
!!    { nx ny nz  sx sy sz }                                                   !
!!    { nkpt  x1 y1 z1 w1 ... xn yn zn wn }                                    !
!!  symmetry         [on, off]                                                 !
!!  wavefunction     [real, complex]                                           !
!!                                                                             !
!!<----------------------------------------------------------------------------!

      SUBROUTINE kpoint_setup(k_points, nkpt_in, nk1, nk2, nk3, k1, k2, k3, xk_in, weight_in)
        IMPLICIT NONE
        CHARACTER (len=80) :: k_points
        INTEGER :: nk1, nk2, nk3
        INTEGER :: k1, k2, k3
        INTEGER :: nkpt_in
        REAL (DP) :: weight_in(:)
        REAL (DP) :: xk_in(:,:)
        REAL (DP) :: weight_sum

        kp%scheme = 'gamma'
        kp%symmetry = .FALSE.
        kp%wfn_type = 0

        IF( ALLOCATED( xk ) )     DEALLOCATE( xk )
        IF( ALLOCATED( weight ) ) DEALLOCATE( weight )

        IF( TRIM( k_points ) /= 'gamma' ) &
          CALL errore( ' kpoint_setup ', ' only gamma is allowed for CP MD, use PW instead ', 1 )

! ... Kpoint type
        SELECT CASE ( TRIM(k_points) )
          CASE ( 'gamma', 'default' )
          CASE ( 'automatic' )
            CALL errore(' kpoint_setup ',' k_points = '//TRIM(k_points)//' not yet implemented ', 1 )
          CASE ( 'tpiba' )
            kp%scheme = 'general'
            kp%symmetry = .FALSE.
            kp%wfn_type = 1
          CASE ( 'crystal' )
            CALL errore(' kpoint_setup ',' k_points = '//TRIM(k_points)//' not yet implemented ', 1 )
          CASE DEFAULT
            CALL errore(' kpoint_setup ',' unknown k_points '//TRIM(k_points), 1 )
        END SELECT

        kp%nkpt = nkpt_in
        kp%nk1 = nk1
        kp%nk2 = nk2
        kp%nk3 = nk3
        kp%nk1 = k1
        kp%nk2 = k2
        kp%nk3 = k3
        kp%shift = 0.0d0
        kp%gamma_only = .FALSE.

        SELECT CASE (kp%scheme)
        CASE DEFAULT
          CALL errore(' kpoint_setup ',' unknown Scheme '//TRIM(kp%scheme), 1)
        CASE ('gamma')
          kp%nkpt = 1
          ALLOCATE( xk(3,1), weight(1) )
          kp%xk => xk
          kp%weight => weight
          kp%xk = 0.0_DP
          kp%weight = 1.0_DP
          kp%gamma_only = .TRUE.
        CASE ('monkhorst-pack')
          kp%nk1 = nk1
          kp%nk2 = nk2
          kp%nk3 = nk3
        CASE ('macdonald')
          kp%nk1 = nk1
          kp%nk2 = nk2
          kp%nk3 = nk3
          kp%shift = 0.0d0
        CASE ('general')
          kp%nkpt = nkpt_in
          ALLOCATE( xk(3,SIZE(xk_in,2)), weight(SIZE(xk_in,2)) )
          kp%xk => xk
          kp%weight => weight
          kp%xk = xk_in
! ...     normalize and set k points weights
          kp%weight = weight_in
          weight_sum = sum(kp%weight)
          kp%weight = kp%weight / weight_sum
        END SELECT
        RETURN
      END SUBROUTINE kpoint_setup

!------------------------------------------------------------------------------!

      SUBROUTINE kpoint_info(punit)
        IMPLICIT NONE
        INTEGER, INTENT (IN) :: punit
        INTEGER :: i

        WRITE (punit,*)
        WRITE (punit,'(3X,A)') 'K points'
        WRITE (punit,'(3X,A)') '--------'
        IF (kp%scheme=='gamma') THEN
          WRITE (punit,'(3X,A)') 'Gamma-point calculation'
          WRITE (punit,'(3X,A)') 'Wavefunction type: REAL'
        ELSE
          WRITE (punit,'(3X,A,1X,A)') 'K-point scheme: ', adjustr(kp%scheme)
          IF (kp%scheme=='monkhorst-pack') THEN
            WRITE (punit,'(3X,A,3I5)') 'K-Point grid  : ', kp%nk1, kp%nk2, kp%nk3
          ELSE IF (kp%scheme=='macdonald') THEN
            WRITE (punit,'(3X,A,3I5)') 'K-Point grid  : ', kp%nk1, kp%nk2, kp%nk3
            WRITE (punit,'(3X,A,3F10.4)') 'K-Point shift : ', kp%shift
          END IF
          IF (kp%symmetry) THEN
            WRITE (punit,'(3X,A)') 'K-Point symmetry:  ON'
          ELSE
            WRITE (punit,'(3X,A)') 'K-Point symmetry:  OFF'
          END IF
          IF (kp%wfn_type==0) THEN
            WRITE (punit,'(3X,A)') 'Wavefunction type: REAL'
          ELSE
            WRITE (punit,'(3X,A)') 'Wavefunction type: COMPLEX'
          END IF
          WRITE (punit,'(3X,A,I3)') 'Number of K-points: ', kp%nkpt
          WRITE (punit,'(3X,A,T19,A,T37,A,T52,A,T67,A)') &
            ' Number ', 'Weight', 'X', 'Y', 'Z'
          DO i = 1, kp%nkpt
            WRITE (punit,'(3X,A,I5,3X,4F15.5)') &
            ' ', i, kp%weight(i), kp%xk(1,i), kp%xk(2,i), kp%xk(3,i)
          END DO
        END IF
      END SUBROUTINE kpoint_info

!------------------------------------------------------------------------------!

      SUBROUTINE brillouin_info(kp,punit)
        IMPLICIT NONE
        TYPE (kpoints), INTENT (IN) :: kp
        INTEGER, INTENT (IN) :: punit
        INTEGER :: i

        IF (kp%scheme=='gamma') THEN
          WRITE (punit,*)
          WRITE (punit,'(A,T57,A)') ' BRILLOUIN|', ' Gamma-point calculation'
          WRITE (punit,'(A,T76,A)') ' BRILLOUIN| Wavefunction type', ' REAL'
        ELSE
          WRITE (punit,*)
          WRITE (punit,'(A,T61,A)') ' BRILLOUIN| K-point scheme ', &
            adjustr(kp%scheme)
          IF (kp%scheme=='monkhorst-pack') THEN
            WRITE (punit,'(A,T66,3I5)') ' BRILLOUIN| K-Point grid', kp%nk1, kp%nk2, kp%nk3
          ELSE IF (kp%scheme=='macdonald') THEN
            WRITE (punit,'(A,T66,3I5)') ' BRILLOUIN| K-Point grid', kp%nk1, kp%nk2, kp%nk3
            WRITE (punit,'(A,T51,3F10.4)') ' BRILLOUIN| K-Point shift', &
              kp%shift
          END IF
          IF (kp%symmetry) THEN
            WRITE (punit,'(A,T76,A)') ' BRILLOUIN| K-Point symmetry', '   ON'
          ELSE
            WRITE (punit,'(A,T76,A)') ' BRILLOUIN| K-Point symmetry', '  OFF'
          END IF
          IF (kp%wfn_type==0) THEN
            WRITE (punit,'(A,T76,A)') ' BRILLOUIN| Wavefunction type', ' REAL'
          ELSE
            WRITE (punit,'(A,T73,A)') ' BRILLOUIN| Wavefunction type', &
              ' COMPLEX'
          END IF
          WRITE (punit,'(A,T71,I10)') ' BRILLOUIN| Number of K-points ', &
            kp%nkpt
          WRITE (punit,'(A,T30,A,T48,A,T63,A,T78,A)') ' BRILLOUIN| Number ', &
            'Weight', 'X', 'Y', 'Z'
          DO i = 1, kp%nkpt
            WRITE (punit,'(A,I5,3X,4F15.5)') ' BRILLOUIN| ', i, kp%weight(i), &
              kp%xk(1,i), kp%xk(2,i), kp%xk(3,i)
          END DO
        END IF
      END SUBROUTINE brillouin_info

      INTEGER FUNCTION get_kpoints_number()
        get_kpoints_number = kp%nkpt
        RETURN
      END FUNCTION get_kpoints_number

!------------------------------------------------------------------------------!
    END MODULE brillouin
!------------------------------------------------------------------------------!
