!
! Copyright (C) 2002-2005 FPMD-CPV groups
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! ---------------------------------------------------------------------- !
      MODULE runsd_module
! ---------------------------------------------------------------------- !

        USE kinds, ONLY: DP

        IMPLICIT NONE
        SAVE

        PRIVATE

        REAL(DP), PRIVATE :: old_clock_value = 0.0d0

        PUBLIC :: runsd

! ---------------------------------------------------------------------- !
      CONTAINS
! ---------------------------------------------------------------------- !


!  -----------------------------------------------------------------------
!  BEGIN manual

      SUBROUTINE runsd &
                 ( tortho, tprint, tforce, rhoe, atoms_0, bec, becdr, eigr, &
                   vkb, ei1, ei2, ei3, sfac, c0, cm, cp, cdesc, tcel, ht0,  &
                   fi, ei, vpot, doions, edft, maxnstep, sdthr )

!  this routine computes the electronic ground state via steepest descent
!  END manual

! ... declare modules
      USE energies,             ONLY: dft_energy_type, print_energies
      USE check_stop,           ONLY: check_stop_now
      USE io_global,            ONLY: ionode
      USE io_global,            ONLY: stdout
      USE cell_module,          ONLY: boxdimensions
      USE wave_types,           ONLY: wave_descriptor
      USE potentials,           ONLY: kspotential
      USE atoms_type_module,    ONLY: atoms_type
      USE cp_interfaces,        ONLY: runcp, update_wave_functions, strucf, phfacs
      USE control_flags,        ONLY: force_pairing
      use grid_dimensions,      only: nr1, nr2, nr3
      USE reciprocal_vectors,   ONLY: mill_l
      USE gvecp,                ONLY: ngm


      IMPLICIT NONE

! ... declare subroutine arguments
      LOGICAL   :: tortho, tprint, tforce, tcel, doions
      TYPE (atoms_type), INTENT(INOUT) :: atoms_0
      COMPLEX(DP), INTENT(INOUT) :: c0(:,:), cm(:,:), cp(:,:)
      TYPE (wave_descriptor) :: cdesc
      REAL(DP) :: rhoe(:,:)
      COMPLEX(DP) :: sfac(:,:)
      COMPLEX(DP) :: eigr(:,:)
      COMPLEX(DP) :: vkb(:,:)
      COMPLEX(DP) :: ei1(:,:)
      COMPLEX(DP) :: ei2(:,:)
      COMPLEX(DP) :: ei3(:,:)
      TYPE (boxdimensions), INTENT(INOUT) ::  ht0
      REAL(DP)  :: fi(:)
      REAL(DP) :: bec(:,:)
      REAL(DP) :: becdr(:,:,:)
      TYPE (dft_energy_type) :: edft

      REAL(DP)    :: ei(:,:)
      REAL(DP)    :: vpot(:,:)

      INTEGER   :: maxnstep   !  maximum number of iteration
      REAL(DP) :: sdthr      !  threshold for convergence 

! ... declare other variables
      LOGICAL :: ttsde, ttprint, ttforce, ttstress, gzero, ttortho
      LOGICAL :: gammasym

      REAL(DP) :: s0, s1, s2, s3, s4, s5, s6, seconds_per_iter
      REAL(DP) :: eold, ekinc, fccc
      REAL(DP) :: ekinc_old, emin, demin

      INTEGER :: ispin, nspin, iter, ierr

      REAL(DP), EXTERNAL :: cclock

! ... end of declarations
!  ----------------------------------------------

      nspin       = cdesc%nspin
      doions      = .FALSE.
      eold        = 1.0d10  ! a large number
      ttsde       = .TRUE.
      ttprint     = .FALSE.
      ttforce     = .FALSE.
      ttstress    = .FALSE.
      ttortho     = .TRUE.
      fccc        = 1.0d0
      gzero       = cdesc%gzero
      gammasym    = cdesc%gamma

      IF( force_pairing ) &
        CALL errore( ' runsd ', ' force pairing not implemented ', 1 )

      IF( ionode ) THEN
        WRITE( stdout,'(/,12X,"Steepest Descent Optimizations for electron, starting ...")' )
        WRITE( stdout,'(  12X,"iter     erho          derho       ekinc      seconds")' )
      END IF

      old_clock_value = cclock()

      CALL phfacs( ei1, ei2, ei3, eigr, mill_l, atoms_0%taus, nr1, nr2, nr3, atoms_0%nat )
      CALL strucf( sfac, ei1, ei2, ei3, mill_l, ngm )
      CALL prefor( eigr, vkb )

      STEEPEST_DESCENT: DO iter = 1, maxnstep

        s1 = cclock()

        CALL kspotential( 1, ttprint, ttforce, ttstress, rhoe, atoms_0, &
                          bec, becdr, eigr, ei1, ei2, ei3, sfac, c0, cdesc, tcel, ht0,  &
                          fi, vpot, edft )

        s2 = cclock()

        CALL runcp( ttprint, ttortho, ttsde, cm, c0, cp, vpot, vkb, fi, ekinc, ht0, ei, bec, fccc)

        emin  = edft%etot
        demin = eold - emin
        eold  = emin

        s0 = cclock()
        seconds_per_iter = ( s0 - old_clock_value )
        old_clock_value = s0

        IF( ionode ) THEN
            WRITE( stdout,113) iter, emin, demin, ekinc, seconds_per_iter
113         FORMAT(10X,I5,2X,F14.6,2X,3D12.4)
        END IF
        
        CALL update_wave_functions( cm, c0, cp )

        s6 = cclock()

! ...   check for exit
        IF (check_stop_now()) THEN
          EXIT STEEPEST_DESCENT
        END IF
        IF( ekinc .LT. sdthr ) THEN
          IF(ionode) WRITE( stdout,fmt="(12X,'runsd: convergence achieved succesfully')")
          doions = .TRUE.
          EXIT STEEPEST_DESCENT
        END IF
        ekinc_old = ekinc
      END DO STEEPEST_DESCENT

! ... set wave functions velocity to 0
      cm = c0

      IF( tforce ) THEN
        atoms_0%for = 0.0d0
        CALL kspotential( 1, ttprint, tforce, ttstress, rhoe, &
          atoms_0, bec, becdr, eigr, ei1, ei2, ei3, sfac, c0, cdesc, tcel, ht0, fi, vpot, edft )
        IF(ionode ) THEN
          WRITE( stdout,fmt="(12X,'runsd: fion and edft calculated = ',F14.6)") edft%etot
        END IF
      END IF

      IF( (iter .GT. maxnstep) .AND. ionode) THEN
        WRITE( stdout,fmt= &
        "(12X,'runsd: convergence not achieved, maximum number of iteration exceeded')")
      END IF

      RETURN
      END SUBROUTINE runsd

! ---------------------------------------------------------------------- !
      END MODULE runsd_module
! ---------------------------------------------------------------------- !
