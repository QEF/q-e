!----------------------------------------------------------------------------------------
!>
!!
!----------------------------------------------------------------------------------------
MODULE environ_cp_module
    !------------------------------------------------------------------------------------
    !! Contains environment contribution related routines.
#if defined (__ENVIRON)
    !
    USE environ_base_module, ONLY: environ
    !
    USE kinds, ONLY: DP
    USE global_version, ONLY: version_number
    USE io_global, ONLY: stdout
    !
    USE fft_base, ONLY: dfftp
    !
    USE electrons_base, ONLY: nspin, nelt
    !
    !------------------------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    PRIVATE
    !
    PUBLIC :: add_environ_potential, calc_environ_potential
    !
    !------------------------------------------------------------------------------------
CONTAINS
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE add_environ_potential(v)
        !--------------------------------------------------------------------------------
        !! Add environment contribution to local potential.
        !
        IMPLICIT NONE
        !
        REAL(DP), INTENT(INOUT) :: v(dfftp%nnr, nspin)
        !
        INTEGER :: i
        REAL(DP) :: dvtot(dfftp%nnr)
        !
        CHARACTER(LEN=80) :: routine = 'add_environ_potential'
        !
        !--------------------------------------------------------------------------------
        !
        dvtot = 0.5D0 * environ%main%get_dvtot(dfftp%nnr) ! Rydberg to Hartree
        !
        IF (nspin .EQ. 1) THEN
            !
!$omp parallel do
            DO i = 1, dfftp%nnr
                v(i, 1) = v(i, 1) + dvtot(i)
            END DO
!$omp end parallel do
            !
        ELSE IF (nspin .EQ. 2) THEN
            !
!$omp parallel do
            DO i = 1, dfftp%nnr
                v(i, 1) = v(i, 1) + dvtot(i)
                v(i, 2) = v(i, 2) + dvtot(i)
            END DO
!$omp end parallel do
            !
        END IF
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE add_environ_potential
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE calc_environ_potential(rhoin, nfi)
        !--------------------------------------------------------------------------------
        !! update electrons-related quantities in environment and calculate the environment 
        !! contribution to the local potential.
        !
        IMPLICIT NONE
        !
        REAL(DP), INTENT(IN) :: rhoin(dfftp%nnr, nspin)
        INTEGER, INTENT(IN) :: nfi
        !
        LOGICAL :: update_venviron = .FALSE.
        REAL(DP) :: rhoaux(dfftp%nnr)
        !
        CHARACTER(LEN=80) :: routine = 'calc_environ_potential'
        !
        !--------------------------------------------------------------------------------
        ! update electrons-related quantities in environ
        !
        rhoaux = rhoin(:, 1)
        !
        IF (version_number == '6.3') THEN
            IF (nspin == 2) rhoaux = rhoaux + rhoin(:, 2)
        END IF
        !
        CALL environ%main%update_electrons(dfftp%nnr, rhoaux, REAL(nelt, DP))
        !
        !--------------------------------------------------------------------------------
        ! Environ contribution to the local potential
        !
        update_venviron = nfi > environ%setup%get_nskip() .OR. environ%setup%is_restart()
        !
        IF (update_venviron .AND. environ%get_verbosity() > 1) WRITE (stdout, 1000)
        !
        CALL environ%calc%potential(update_venviron)
        !
        !--------------------------------------------------------------------------------
        !
1000    FORMAT(/, 5X, "add environment contribution to local potential")
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE calc_environ_potential
    !------------------------------------------------------------------------------------
    !
#endif
    !------------------------------------------------------------------------------------
END MODULE environ_cp_module
!----------------------------------------------------------------------------------------
