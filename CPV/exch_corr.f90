!
! Copyright (C) 2002 FPMD group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

      MODULE exchange_correlation

        USE kinds
        USE local_density_approximation
        USE local_spin_density

        IMPLICIT NONE
        SAVE

        PRIVATE

        LOGICAL :: use_table

        PUBLIC :: v2gc, vofxc_lsd, exch_corr_energy, vofxc_lda
        PUBLIC :: exch_corr_setup, exch_corr_print_info
        PUBLIC :: deallocate_lda
        PUBLIC :: tgc, narray

      CONTAINS

        SUBROUTINE exch_corr_setup(xc_type, narray_i, rmxxc_i)
          CHARACTER(LEN = 20) xc_type
          INTEGER, INTENT(IN) :: narray_i
          REAL(dbl), INTENT(IN) :: rmxxc_i
            CALL lsd_setup(xc_type)
            CALL lda_setup(narray_i, rmxxc_i)
            IF( xc_type == 'LDA' ) THEN
              use_table = .TRUE.
            ELSE
              use_table = .FALSE.
            END IF
          RETURN
        END SUBROUTINE exch_corr_setup

        SUBROUTINE exch_corr_print_info(iunit)
          INTEGER, INTENT(IN) :: iunit
          WRITE(iunit,800)
          IF( use_table ) THEN
            CALL lda_print_info(iunit)
          ELSE
            CALL lsd_print_info(iunit)
          END IF
800       FORMAT(//,3X,'Exchange and correlations functionals',/ &
                   ,3X,'-------------------------------------')
          RETURN
        END SUBROUTINE exch_corr_print_info


        SUBROUTINE exch_corr_energy(rhoetr, rhoetg, grho, vpot, &
          sxc, vxc, v2xc, v2xc2, gv)
          USE kinds, ONLY: dbl
          USE cp_types, ONLY: recvecs
          REAL (dbl) :: rhoetr(:,:,:,:)
          COMPLEX(dbl) :: rhoetg(:,:)
          REAL (dbl) :: grho(:,:,:,:,:)
          REAL (dbl) :: vpot(:,:,:,:)
          REAL (dbl) :: sxc, vxc
          REAL (dbl) :: v2xc(:,:,:,:)
          REAL (dbl) :: v2xc2(:,:,:)
          TYPE (recvecs), INTENT(IN) :: gv
          INTEGER :: nspin, ispin

          nspin = SIZE(rhoetr, 4)
          IF(.NOT.tgc) THEN
! ...       no gradient correction
! ...       vpot = vxc(rhoetr); vpot(r) <-- u(r)
            IF( (nspin .EQ. 1) .AND. use_table ) THEN
              CALL vofxc_lda(rhoetr(:,:,:,1), vpot(:,:,:,1), sxc, vxc)
            ELSE
              CALL vofxc_lsd(rhoetr, vpot, sxc, vxc)
            END IF
          ELSE
! ...       gradient correction
            CALL vofxc_lsd(rhoetr, vpot, grho, v2xc, v2xc2, sxc, vxc)
            CALL v2gc(v2xc, v2xc2, grho, rhoetr, vpot, gv)
          END If
          RETURN
        END SUBROUTINE

      END MODULE exchange_correlation
