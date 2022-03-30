!----------------------------------------------------------------------------------------
!>
!!
!----------------------------------------------------------------------------------------
MODULE environ_td_module
    !------------------------------------------------------------------------------------
#if defined (__ENVIRON)
    !
    USE environ_base_module, ONLY: environ
    !
    USE kinds, ONLY: DP
    USE io_global, ONLY: stdout
    !
    USE fft_base, ONLY: dfftp
    !
    USE lsda_mod, ONLY: nspin
    USE lr_variables, ONLY: davidson
    !
    !------------------------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    PRIVATE
    !
    PUBLIC :: calc_environ_dpotential
    !
    !------------------------------------------------------------------------------------
CONTAINS
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE calc_environ_dpotential(drho, dv)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        REAL(DP), INTENT(IN) :: drho(dfftp%nnr, nspin)
        !
        REAL(DP), INTENT(INOUT) :: dv(dfftp%nnr, nspin)
        !
        !--------------------------------------------------------------------------------
        !
        IF (.NOT. davidson) WRITE (stdout, 1000)
        !
        IF (environ%setup%optical_permittivity == 1.D0) WRITE (stdout, 1002)
        !
        CALL environ%main%update_response(dfftp%nnr, drho(:, 1))
        !
        CALL environ%calc%dpotential(dfftp%nnr, dv(:, 1))
        !
        !--------------------------------------------------------------------------------
        !
1000    FORMAT(5X, "Calculate Environ contribution to response potential")
        !
1002    FORMAT("Warning: permittivity is set to 1.0 - no Environ contribution")
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE calc_environ_dpotential
    !------------------------------------------------------------------------------------
    !
#endif
    !------------------------------------------------------------------------------------
END MODULE environ_td_module
!----------------------------------------------------------------------------------------
