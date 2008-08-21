!
! Copyright (C) 2002-2005 Quantum-ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
MODULE from_restart_module
  !----------------------------------------------------------------------------
  !
  IMPLICIT NONE
  SAVE
  !
  PRIVATE
  !
  PUBLIC :: from_restart
  !
  CONTAINS
  !
  !--------------------------------------------------------------------------
  SUBROUTINE from_restart( sfac, eigr, ei1, ei2, ei3, bec, becdr,         &
                           taub, irb, eigrb, b1, b2, b3, nfi, rhog,       &
                           rhor, rhos, rhoc,                              &
                           lambda, lambdam, lambdap, ema0bg, dbec, bephi, &
                           becp, htm, ht0, vpot, atoms0, edft )

    !--------------------------------------------------------------------------
    !
    USE kinds,                ONLY : DP
    USE control_flags,        ONLY : tzeroc, tzerop, tzeroe, dt_old
    USE ions_positions,       ONLY : taus, tausm
    USE cell_base,            ONLY : boxdimensions
    USE uspp,                 ONLY : vkb, nkb
    USE wavefunctions_module, ONLY : c0, cm, phi => cp
    USE io_global,            ONLY : stdout, ionode
    USE ions_nose,            ONLY : xnhp0, xnhpm
    USE energies,             ONLY : eself, dft_energy_type
    USE time_step,            ONLY : delt
    USE atoms_type_module,    ONLY : atoms_type
    !
    COMPLEX(DP) :: eigr(:,:), ei1(:,:), ei2(:,:), ei3(:,:)
    COMPLEX(DP) :: eigrb(:,:)
    INTEGER     :: irb(:,:)
    REAL(DP)    :: bec(:,:), becdr(:,:,:)
    REAL(DP)    :: taub(:,:)
    REAL(DP)    :: b1(:), b2(:), b3(:)
    INTEGER     :: nfi
    COMPLEX(DP) :: sfac(:,:)
    COMPLEX(DP) :: rhog(:,:)
    REAL(DP)    :: rhor(:,:), rhos(:,:), rhoc(:)
    REAL(DP)    :: lambda(:,:,:), lambdam(:,:,:), lambdap(:,:,:)
    REAL(DP)    :: ema0bg(:)
    REAL(DP)    :: dbec(:,:,:,:)
    REAL(DP)    :: bephi(:,:), becp(:,:)
    TYPE(boxdimensions)        :: htm, ht0
    REAL(DP)    :: vpot(:,:)
    TYPE(atoms_type) :: atoms0
    TYPE(dft_energy_type) :: edft
    !
    CALL from_restart_x( ht0, htm, phi, c0, cm, lambdap, lambda, lambdam, &
            ei1, ei2, ei3, eigr, sfac, vkb, nkb, bec, dbec, taub, irb, eigrb )
    !
    edft%eself = eself
    !
    IF( tzerop .or. tzeroe .or. tzeroc ) THEN
       IF( .not. ( tzerop .and. tzeroe .and. tzeroc ) ) THEN
          IF( ionode ) THEN
             WRITE( stdout, * ) 'WARNING setting to ZERO ions, electrons and cell velocities without '
             WRITE( stdout, * ) 'setting to ZERO all velocities could generate meaningles trajectories '
          END IF
       END IF
    END IF
    !
    ! dt_old should be -1.0 here if untouched ...
    !
    if ( dt_old > 0.0d0 ) then
       tausm = taus - (taus-tausm)*delt/dt_old
       xnhpm = xnhp0 - (xnhp0-xnhpm)*delt/dt_old
       WRITE( stdout, '(" tausm & xnhpm were rescaled ")' )
    endif
    !
    RETURN
    !
END SUBROUTINE from_restart




SUBROUTINE from_restart_x( &
   ht0, htm, phi, c0, cm, lambdap, lambda, lambdam, ei1, ei2, ei3, eigr, sfac, vkb, nkb, &
   bec, dbec, taub, irb, eigrb )
   !
   USE kinds,                 ONLY : DP
   USE control_flags,         ONLY : program_name, tbeg, taurdr, tfor, tsdp, tv0rd, &
                                     iprsta, tsde, tzeroe, tzerop, nbeg, tranp, amprp, &
                                     tzeroc, force_pairing, trhor, ampre, trane, tpre
   USE electrons_module,      ONLY : occn_info
   USE electrons_base,        ONLY : nspin, iupdwn, nupdwn, f, nbsp
   USE io_global,             ONLY : ionode, ionode_id, stdout
   USE cell_base,             ONLY : ainv, h, hold, deth, r_to_s, s_to_r, boxdimensions, &
                                     velh, a1, a2, a3
   USE ions_base,             ONLY : na, nsp, iforce, vel_srt, nat, randpos
   USE time_step,             ONLY : tps, delt
   USE ions_positions,        ONLY : taus, tau0, tausm, taum, vels, fion, fionm, set_velocities
   USE grid_dimensions,       ONLY : nr1, nr2, nr3
   USE reciprocal_vectors,    ONLY : mill_l
   USE printout_base,         ONLY : printout_pos
   USE gvecs,                 ONLY : ngs
   USE gvecw,                 ONLY : ngw
   USE cp_interfaces,         ONLY : phfacs, strucf
   USE energies,              ONLY : eself
   USE wave_base,             ONLY : rande_base
   USE mp_global,             ONLY : me_image, mpime
   USE efield_module,         ONLY : efield_berry_setup,  tefield, &
                                     efield_berry_setup2, tefield2
   USE small_box,             ONLY : ainvb
   USE uspp,                  ONLY : okvan
   USE core,                  ONLY : nlcc_any
   !
   IMPLICIT NONE
   !
   TYPE (boxdimensions) :: ht0, htm
   COMPLEX(DP)          :: phi(:,:)
   COMPLEX(DP)          :: c0(:,:)
   COMPLEX(DP)          :: cm(:,:)
   REAL(DP)             :: lambda(:,:,:), lambdam(:,:,:), lambdap(:,:,:)
   COMPLEX(DP)          :: eigr(:,:)
   COMPLEX(DP)          :: ei1(:,:)
   COMPLEX(DP)          :: ei2(:,:)
   COMPLEX(DP)          :: ei3(:,:)
   COMPLEX(DP)          :: sfac(:,:)
   COMPLEX(DP)          :: vkb(:,:)
   INTEGER, INTENT(IN)  :: nkb
   REAL(DP)             :: bec(:,:)
   REAL(DP)             :: dbec(:,:,:,:)
   REAL(DP)             :: taub(:,:)
   INTEGER              :: irb(:,:)
   COMPLEX(DP)          :: eigrb(:,:)
   !
   ! ... We are restarting from file recompute ainv
   !
   CALL invmat( 3, h, ainv, deth )
   !
   ! ... Reset total time counter if the run is not strictly 'restart'
   !
   IF ( nbeg < 1 ) tps = 0.D0
   !
   IF ( taurdr ) THEN
      !
      ! ... Input positions read from input file and stored in tau0
      ! ... in readfile, only scaled positions are read
      !
      CALL r_to_s( tau0, taus, na, nsp, ainv )
      !
   END IF
   !
   IF ( ANY( tranp(1:nsp) ) ) THEN
      !
      ! ... Input positions are randomized
      !
      CALL randpos( taus, na, nsp, tranp, amprp, ainv, iforce )
      !
   END IF
   !
   IF ( tzerop .AND. tfor ) THEN
      !
      CALL r_to_s( vel_srt, vels, na, nsp, ainv )
      !
      CALL set_velocities( tausm, taus, vels, iforce, nat, delt )
      !
   END IF
   !
   CALL s_to_r( taus,  tau0, na, nsp, h )
   !
   CALL s_to_r( tausm, taum, na, nsp, h )
   !
   IF ( tzeroc ) THEN
      !
      hold = h
      velh = 0.D0
      !
      htm      = ht0
      ht0%hvel = 0.D0
      !
   END IF
   !
   fion = 0.D0
   !
   IF( force_pairing ) THEN
      cm(:,iupdwn(2):nbsp) = cm(:,1:nupdwn(2))
      c0(:,iupdwn(2):nbsp) = c0(:,1:nupdwn(2))
      lambdap( :, :, 2) =  lambdap( :, :, 1)
      lambda( :, :, 2) =  lambda( :, :, 1)
      lambdam( :, :, 2) = lambdam( :, :, 1)
   END IF 
   !
   IF ( tzeroe ) THEN
      !
      lambdam = lambda
      !
      cm = c0
      !
      WRITE( stdout, '(" Electronic velocities set to zero")' )
      !
   END IF
   !
   IF( program_name == "FPMD" ) THEN
      !
      CALL occn_info( f )
      !
      ! ... diagnostics
      !
      IF ( ionode ) THEN
         !
         WRITE( stdout, 100 )
         !
         IF ( .NOT. tbeg ) THEN
            WRITE( stdout, 110 )
         ELSE
            WRITE( stdout, 120 )
         END IF
         !
         IF ( .NOT. taurdr ) THEN
            WRITE( stdout, 130 )
         ELSE
            WRITE( stdout, 140 )
         END IF
         !
         IF ( tfor .AND. ( .NOT. tsdp ) ) THEN
            !
            IF ( .NOT. tv0rd ) THEN
               !
               IF ( .NOT. tzerop ) THEN
                  WRITE( stdout, 150 )
               ELSE
                  WRITE( stdout, 155 )
               END IF
               !
            ELSE
               !
               WRITE( stdout, 160 )
               !
            END IF
          !
         END IF
         !
         IF ( iprsta > 1 ) &
            CALL printout_pos( stdout, taus, nat, head = 'Scaled positions from restart module' )
         !
         IF ( .NOT. tsde ) THEN
            !
            IF ( .NOT. tzeroe ) THEN
               WRITE( stdout, 170 )
            ELSE
               WRITE( stdout, 180 )
            END IF
            !
         END IF
         !
         WRITE( stdout, * )
         !
      END IF
      !
   END IF
   !
   ! ... computes form factors and initializes nl-pseudop. according
   ! ... to starting cell (from ndr or again standard input)
   !
   IF ( okvan .or. nlcc_any ) THEN
      CALL initbox( tau0, taub, irb, ainv, a1, a2, a3 )
      CALL phbox( taub, eigrb, ainvb )
   END IF
   !
   CALL phfacs( ei1, ei2, ei3, eigr, mill_l, taus, nr1, nr2, nr3, nat )
   !
   CALL strucf( sfac, ei1, ei2, ei3, mill_l, ngs )
   !
   CALL prefor( eigr, vkb )
   !
   CALL formf( .TRUE. , eself )
   !
   IF ( trane ) THEN
      !
      WRITE( stdout, 515 ) ampre
      !
515   FORMAT(   3X,'Initial random displacement of el. coordinates',/ &
                3X,'Amplitude = ',F10.6 )
      !
      CALL rande_base( c0, ampre )

      CALL gram( vkb, bec, nkb, c0, ngw, nbsp )
      !
      IF( force_pairing ) c0(:,iupdwn(2):nbsp) = c0(:,1:nupdwn(2))
      !
      cm = c0
      !
   END IF
   !
   CALL calbec( 1, nsp, eigr, c0, bec )
   !
   IF ( tpre     ) CALL caldbec( ngw, nkb, nbsp, 1, nsp, eigr, c0, dbec )
   !
   IF ( tefield  ) CALL efield_berry_setup( eigr, tau0 )
   IF ( tefield2 ) CALL efield_berry_setup2( eigr, tau0 )
   !
100 FORMAT( /,3X,'MD PARAMETERS READ FROM RESTART FILE',/ &
             ,3X,'------------------------------------' )
110 FORMAT(   3X,'Cell variables From RESTART file' )
120 FORMAT(   3X,'Cell variables From INPUT file' )
130 FORMAT(   3X,'Ions positions From RESTART file' )
140 FORMAT(   3X,'Ions positions From INPUT file' )
150 FORMAT(   3X,'Ions Velocities From RESTART file' )
155 FORMAT(   3X,'Ions Velocities set to ZERO' )
160 FORMAT(   3X,'Ions Velocities From STANDARD INPUT' )
170 FORMAT(   3X,'Electronic Velocities From RESTART file' )
180 FORMAT(   3X,'Electronic Velocities set to ZERO' )
    !
   !
   RETURN
END SUBROUTINE
  !
  !
  !
END MODULE from_restart_module
