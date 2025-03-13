!----------------------------------------------------------------------------------------
!>
!!
!----------------------------------------------------------------------------------------
MODULE environ_base_module
    !------------------------------------------------------------------------------------
#if defined (__ENVIRON)
    !
    USE environ_api, ONLY: environ_interface
    !
    USE kinds, ONLY: DP
    USE io_global, ONLY: ionode, ionode_id, stdout
    !
    USE mp_images, ONLY: intra_image_comm
    USE mp_bands, ONLY: intra_bgrp_comm
    !
    USE control_flags, ONLY: conv_elec, iverbosity
    !
    USE fft_base, ONLY: dfftp
    !
    USE ions_base, ONLY: nat, nsp, ityp, zv
    !
    USE uspp_param, ONLY: upf
    !
    USE input_parameters, ONLY : lgcscf
    !
    !------------------------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    PRIVATE
    !
    PUBLIC :: read_environ_input
    !
    PUBLIC :: init_environ_setup, init_environ_base
    !
    PUBLIC :: clean_environ, check_environ_compatibility
    !
    PUBLIC :: update_environ_ions, update_environ_cell
    !
    PUBLIC :: calc_environ_energy, calc_environ_force
    !
    PUBLIC :: print_environ_summary, print_environ_energies, print_environ_clocks
    !
    !------------------------------------------------------------------------------------
    !
    TYPE(environ_interface), PUBLIC :: environ
    !
    !------------------------------------------------------------------------------------
CONTAINS
    !------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------
    !
    !                                   ADMIN ROUTINES
    !
    !------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE read_environ_input()
        !--------------------------------------------------------------------------------
        !
        CALL environ%init_interface()
        !
        CALL environ%init_io(ionode, ionode_id, intra_image_comm, stdout, ionode)
        !
        CALL environ%read_input()
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE read_environ_input
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE init_environ_setup(prog)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: prog
        !
        CHARACTER(LEN=80) :: routine = 'init_environ_setup'
        !
        !--------------------------------------------------------------------------------
        !
        CALL environ%setup%init()
        !
        IF (PRESENT(prog)) THEN
            IF (prog == 'TD') CALL environ%setup%set_tddfpt(.TRUE.)
        END IF
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE init_environ_setup
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE init_environ_base(at, gcutm, do_comp_mt)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        REAL(DP), INTENT(IN) :: at(3, 3)
        REAL(DP), INTENT(IN) :: gcutm
        LOGICAL, OPTIONAL, INTENT(IN) :: do_comp_mt
        !
        INTEGER :: nr(3)
        !
        CHARACTER(LEN=80) :: routine = 'init_environ_base'
        !
        !--------------------------------------------------------------------------------
        !
        nr(1) = dfftp%nr1
        nr(2) = dfftp%nr2
        nr(3) = dfftp%nr3
        !
        CALL environ%setup%init_cell(intra_bgrp_comm, at, gcutm=gcutm, nr=nr)
        !
        CALL environ%setup%init_numerical(do_comp_mt)
        !
        CALL environ%main%init(nat, nsp, ityp, zv, label=upf%psd, lgcscf=lgcscf)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE init_environ_base
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE clean_environ(prog, lflag)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: prog
        LOGICAL, OPTIONAL, INTENT(IN) :: lflag
        !
        LOGICAL :: local_flag
        !
        CHARACTER(LEN=80) :: routine = 'clean_environ'
        !
        !--------------------------------------------------------------------------------
        !
        IF (PRESENT(lflag)) THEN
            local_flag = lflag
        ELSE
            local_flag = .FALSE.
        END IF
        !
        IF (PRESENT(prog)) THEN
            !
            SELECT CASE (prog)
                !
            CASE ('PW')
                !
                !------------------------------------------------------------------------
                ! When called by PW, but inside a TD calculation, do not clean Environ
                ! variables - they have already been cleaned by TD. The lflag input is
                ! used to fully clean the variable or to only clean variables initialized
                ! during the PW run and not the ones initialized while processing the
                ! input. This allows NEB simulations
                !
                IF (local_flag) THEN
                    CALL environ%destroy(1) ! NEB image reading phase
                ELSE IF (.NOT. environ%setup%is_tddfpt()) THEN
                    CALL environ%destroy(2)
                END IF
                !
            CASE ('TD')
                !
                !------------------------------------------------------------------------
                ! When called by TD, use the flag input variable to specify whether
                ! to clean the PW variables or the TD variables. In both cases, the
                ! variables are fully cleaned (no NEB with TD)
                !
                IF (.NOT. local_flag) THEN
                    CALL environ%destroy(3)
                ELSE
                    CALL environ%destroy(4)
                END IF
                !
            CASE DEFAULT
                CALL errore(routine, "Unexpected calling program", 1)
                !
            END SELECT
            !
        ELSE
            CALL environ%destroy()
        END IF
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE clean_environ
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE check_environ_compatibility(calling_routine)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CHARACTER(LEN=*), INTENT(IN) :: calling_routine
        !
        CHARACTER(LEN=80) :: routine = 'check_environ_compatibility'
        !
        !--------------------------------------------------------------------------------
        !
        CALL errore(routine, "Calculation not compatible with Environ embedding", 1)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE check_environ_compatibility
    !------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------
    !
    !                                  UPDATE ROUTINES
    !
    !------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE update_environ_ions(tau)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        REAL(DP), OPTIONAL, INTENT(IN) :: tau(3, nat)
        !
        CHARACTER(LEN=80) :: routine = 'update_environ_ions'
        !
        !--------------------------------------------------------------------------------
        !
        CALL environ%main%update_ions(nat, tau)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE update_environ_ions
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE update_environ_cell(at)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        REAL(DP), INTENT(IN) :: at(3, 3)
        !
        CHARACTER(LEN=80) :: routine = 'update_environ_cell'
        !
        !--------------------------------------------------------------------------------
        !
        CALL environ%update_cell(at)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE update_environ_cell
    !------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------
    !
    !                                COMPUTATION ROUTINES
    !
    !------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------
    !>
    !! Compute environ contributions to total energy
    !! Note: plugin_etot is set to 0.0_dp right before this routine is called
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE calc_environ_energy(plugin_etot, de_flag)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        LOGICAL, OPTIONAL, INTENT(IN) :: de_flag
        !
        REAL(DP), INTENT(INOUT) :: plugin_etot
        !
        CHARACTER(LEN=80) :: routine = 'calc_environ_energy'
        !
        !--------------------------------------------------------------------------------
        !
        IF (PRESENT(de_flag)) THEN
            IF (de_flag) CALL environ%calc%denergy(plugin_etot)
        END IF
        !
        CALL environ%calc%energy(plugin_etot)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE calc_environ_energy
    !------------------------------------------------------------------------------------
    !>
    !! Compute environ contributions to total energy
    !! Note: plugin_etot is set to 0.0_dp right before this routine is called
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE calc_environ_force(force)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        REAL(DP), INTENT(OUT) :: force(3, nat)
        !
        INTEGER :: i, j
        !
        CHARACTER(LEN=80) :: routine = 'calc_environ_force'
        !
        !--------------------------------------------------------------------------------
        !
        force = 0.D0
        !
        CALL environ%calc%force(nat, force)
        !
        IF (iverbosity > 0) THEN
            WRITE (stdout, 1001)
            !
            DO j = 1, nat
                WRITE (stdout, 1002) j, ityp(j), (force(i, j), i=1, 3)
            END DO
            !
            WRITE (stdout, *)
        END IF
        !
        !--------------------------------------------------------------------------------
        !
1001    FORMAT(5X, "The dielectric solvent contribution to forces")
1002    FORMAT(5X, "atom ", I4, " type ", I2, "   force = ", 3F14.8)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE calc_environ_force
    !------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------
    !
    !                                  OUTPUT ROUTINES
    !
    !------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE print_environ_summary()
        !--------------------------------------------------------------------------------
        !
        CALL environ%setup%print_summary(stdout)
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE print_environ_summary
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE print_environ_energies(prog)
        !--------------------------------------------------------------------------------
        !
        IMPLICIT NONE
        !
        CHARACTER(LEN=*), INTENT(IN) :: prog
        !
        CHARACTER(LEN=80) :: routine = 'print_environ_energies'
        !
        !--------------------------------------------------------------------------------
        !
        CALL environ%main%print_energies(prog)
        !
        IF (prog == 'PW' .AND. conv_elec) CALL environ%setup%print_potential_warning()
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE print_environ_energies
    !------------------------------------------------------------------------------------
    !>
    !!
    !------------------------------------------------------------------------------------
    SUBROUTINE print_environ_clocks()
        !--------------------------------------------------------------------------------
        !
        CALL environ%setup%print_clocks()
        !
        !--------------------------------------------------------------------------------
    END SUBROUTINE print_environ_clocks
    !------------------------------------------------------------------------------------
    !
#endif
    !------------------------------------------------------------------------------------
END MODULE environ_base_module
!----------------------------------------------------------------------------------------
