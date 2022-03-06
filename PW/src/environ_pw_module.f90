#if defined (__ENVIRON)
!----------------------------------------------------------------------------------------
!>
!!
!----------------------------------------------------------------------------------------
MODULE environ_pw_module
    !------------------------------------------------------------------------------------
    !
    USE environ_base_module, ONLY: environ
    !
    USE kinds, ONLY: DP
    USE global_version, ONLY: version_number
    USE io_global, ONLY: stdout
    !
    USE control_flags, ONLY: lscf
    !
    USE fft_base, ONLY: dfftp
    !
    USE scf, ONLY: scf_type
    USE klist, ONLY: nelec
    USE lsda_mod, ONLY: nspin
    !
    !------------------------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    PRIVATE
    !
    PUBLIC :: update_environ_potential, calc_environ_potential
    !
    !------------------------------------------------------------------------------------
CONTAINS
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE update_environ_potential(vltot)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        REAL(DP), INTENT(IN) :: vltot(dfftp%nnr)
        !
        CHARACTER(LEN=80) :: sub_name = 'update_environ_potential'
        !
        !--------------------------------------------------------------------------------
        !
        CALL environ%main%update_potential(dfftp%nnr, vltot)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE update_environ_potential
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE calc_environ_potential(rhoin, converged, dr2, vltot)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        TYPE(scf_type), INTENT(IN) :: rhoin
        LOGICAL, INTENT(IN) :: converged
        REAL(DP), INTENT(IN) :: dr2
        !
        REAL(DP), INTENT(INOUT) :: vltot(dfftp%nnr)
        !
        LOGICAL :: update_venviron = .FALSE.
        INTEGER :: local_verbose = 0
        !
        REAL(DP) :: rhoaux(dfftp%nnr)
        !
        CHARACTER(LEN=80) :: sub_name = 'calc_environ_potential'
        !
        !--------------------------------------------------------------------------------
        ! Reduce output at each scf iteration
        !
        IF (.NOT. lscf .OR. converged) local_verbose = 1
        !
        !--------------------------------------------------------------------------------
        ! update electrons-related quantities in environ
        !
        rhoaux = rhoin%of_r(:, 1)
        !
        IF (version_number == '6.3') THEN
            IF (nspin == 2) rhoaux = rhoaux + rhoin%of_r(:, 2)
        END IF
        !
        CALL environ%main%update_electrons(dfftp%nnr, rhoaux, nelec)
        !
        !--------------------------------------------------------------------------------
        ! Environ contribution to the local potential
        !
        IF (dr2 > 0.0_DP) THEN
            update_venviron = .NOT. converged .AND. dr2 < environ%setup%get_threshold()
        ELSE
            update_venviron = environ%setup%is_restart() .OR. environ%setup%is_tddfpt()
            !
            !----------------------------------------------------------------------------
            ! For subsequent steps of optimization or dynamics, compute Environ
            ! contribution during initialization
            !
            CALL environ%setup%set_restart(.TRUE.)
            !
        END IF
        !
        IF (update_venviron) WRITE (stdout, 1000)
        !
        CALL environ%calc%potential(update_venviron, local_verbose)
        !
        vltot = environ%main%get_vzero(dfftp%nnr) + environ%main%get_dvtot(dfftp%nnr)
        !
        IF (.NOT. lscf .OR. converged) CALL environ%main%print_potential_shift()
        !
        !--------------------------------------------------------------------------------
        !
1000    FORMAT(/, 5X, "add environment contribution to local potential")
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE calc_environ_potential
    !------------------------------------------------------------------------------------
    !
    !------------------------------------------------------------------------------------
END MODULE environ_pw_module
!----------------------------------------------------------------------------------------
#endif
