!
! Copyright (C) 2002 FPMD group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

#include "f_defs.h"

!  AB INITIO COSTANT PRESSURE MOLECULAR DYNAMICS
!  ----------------------------------------------

!=----------------------------------------------------------------------------=!
      MODULE from_scratch_module
!=----------------------------------------------------------------------------=!

        IMPLICIT NONE
        SAVE

        PRIVATE

        PUBLIC :: fromscratch

!=----------------------------------------------------------------------------=!
      CONTAINS
!=----------------------------------------------------------------------------=!

!  ----------------------------------------------
!  BEGIN manual

      SUBROUTINE fromscratch(gv, kp, ps, rhoe, desc, cm, c0, cdesc, eigr, sfac, &
        fi, ht, atoms, fnl, vpot, edft)

!  (describe briefly what this routine does...)
!  ----------------------------------------------
!  END manual


! ... declare modules
      USE kinds
      USE cp_types, ONLY: recvecs, pseudo, phase_factors
      USE wave_types, ONLY: wave_descriptor
      USE atoms_type_module, ONLY: atoms_type
      USE wave_functions, ONLY: moveelect, gram, fixwave
      USE wave_base, ONLY: wave_steepest
      USE charge_density, ONLY: rhoofr
      USE phase_factors_module, ONLY: strucf
      USE cell_module, only: boxdimensions
      USE electrons_module, ONLY: nspin, pmss
      USE cp_electronic_mass, ONLY: emass
      USE ions_module, ONLY: taui, cdmi, set_reference_positions, &
          constraints_setup
      USE mp, ONLY: mp_end
      USE nl, ONLY: nlrh_m
      USE energies, ONLY: dft_energy_type, debug_energies
      USE potentials, ONLY: vofrhos
      USE forces, ONLY: dforce_all
      USE orthogonalize, ONLY: ortho
      USE brillouin, ONLY: kpoints
      USE pseudo_projector, ONLY: projector
      USE control_flags, ONLY: tcarpar, tfor, thdyn, tortho, prn, force_pairing
      USE charge_types, ONLY: charge_descriptor
      USE time_step, ONLY: delt

      IMPLICIT NONE

! ... declare subroutine arguments

      TYPE (atoms_type) :: atoms
      TYPE (phase_factors) :: eigr
      TYPE (recvecs) :: gv 
      TYPE (kpoints) :: kp 
      REAL(dbl) :: rhoe(:,:,:,:)
      COMPLEX(dbl) :: sfac(:,:)
      TYPE (charge_descriptor), INTENT(IN) :: desc
      TYPE (wave_descriptor), INTENT(IN) :: cdesc
      COMPLEX(dbl), INTENT(INOUT) :: cm(:,:,:,:), c0(:,:,:,:)
      TYPE (boxdimensions) :: ht 
      TYPE (pseudo) :: ps
      REAL(dbl) :: fi(:,:,:)
      TYPE (projector) :: fnl(:,:)
      REAL (dbl)    ::  vpot(:,:,:,:)
      TYPE (dft_energy_type) :: edft


! ... declare other variables

      LOGICAL, PARAMETER :: ttsde = .TRUE.
      LOGICAL, PARAMETER :: ttforce = .TRUE.
      LOGICAL, PARAMETER :: ttprint = .TRUE.
      REAL(dbl), PARAMETER :: svar1 = 1.d0
      REAL(dbl), PARAMETER :: svar2 = 0.d0
      INTEGER, PARAMETER :: nfi = 0

      INTEGER :: ib, ik, is, ierr
      REAL(dbl)  timepre, dum_kin

      COMPLEX (dbl), ALLOCATABLE :: eforce(:,:,:)
      REAL (dbl), ALLOCATABLE :: dt2bye( : )


!  end of declarations
!  ----------------------------------------------

      atoms%for = 0.0d0
      CALL strucf(sfac, atoms, eigr, gv)
      edft%enl = nlrh_m(cm, cdesc, ttforce, atoms, fi, gv, kp, fnl, ps%wsg, ps%wnl, eigr)
      CALL rhoofr(gv, kp, cm, cdesc, fi, rhoe, desc, ht)
      CALL vofrhos(ttprint, prn, rhoe, desc, tfor, thdyn, ttforce, atoms, &
           gv, kp, fnl, vpot, ps, cm, cdesc, fi, eigr, sfac, timepre, ht, edft)

      ! CALL debug_energies( edft ) ! DEBUG

      ! IF( .FALSE. ) THEN ! DEBUG

      IF( tcarpar ) THEN

        ALLOCATE( dt2bye( SIZE( pmss ) ) )
        dt2bye = delt * delt / pmss

        IF( .NOT. force_pairing ) THEN

          DO is = 1, cdesc%nspin

            ALLOCATE( eforce( SIZE( cm, 1), SIZE( cm, 2), SIZE( cm, 3 ) ), STAT=ierr )
            IF( ierr /= 0 ) CALL errore(' fromscra ', ' allocating eforce ', ierr )

            CALL dforce_all( is, cm(:,:,:,is), cdesc, fi(:,:,is), eforce(:,:,:), &
              gv, vpot(:,:,:,is), fnl(:,is), eigr, ps)    

            DO ik = 1, cdesc%nkl
              DO ib = 1, cdesc%nbl( is )
                CALL wave_steepest( c0(:,ib,ik,is), cm(:,ib,ik,is), dt2bye, eforce(:,ib,ik) )
              END DO
              CALL fixwave( is, c0(:,:,ik,is), cdesc, gv%kg_mask_l(:,ik) )
            END DO


            DEALLOCATE( eforce, STAT=ierr )
            IF( ierr /= 0 ) CALL errore(' fromscra ', ' deallocating eforce ', ierr )
          
          END DO

        ELSE

          c0 = cm

        END IF

        IF( tortho .AND. ( .NOT. force_pairing ) ) THEN
          CALL ortho( cm, c0, cdesc, pmss, emass )
        ELSE
          CALL gram( c0, cdesc )
        END IF

        DEALLOCATE( dt2bye )

      ELSE

        c0 = cm

      END IF

      CALL set_reference_positions(cdmi, taui, atoms, ht)
      CALL constraints_setup(ht, atoms)

      RETURN
      END SUBROUTINE fromscratch

!=----------------------------------------------------------------------------=!
      END MODULE
!=----------------------------------------------------------------------------=!
