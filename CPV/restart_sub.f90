!
! Copyright (C) 2002 FPMD group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

MODULE from_restart_module

  IMPLICIT NONE
  SAVE

  PRIVATE

  PUBLIC :: from_restart

  INTERFACE from_restart
    MODULE PROCEDURE from_restart_cp, from_restart_fpmd
  END INTERFACE

CONTAINS

!=----------------------------------------------------------------------------=!

  SUBROUTINE from_restart_cp( sfac, eigr, ei1, ei2, ei3, bec, becdr, tfirst, eself, fion, &
      taub, irb, eigrb, b1, b2, b3, nfi, rhog, rhor, rhos, rhoc, enl, ekin, stress,  &
      detot, enthal, etot, lambda, lambdam, lambdap, ema0bg, dbec, delt,  &
      bephi, becp, velh, dt2bye, iforce, fionm, nbeg, xnhe0, xnhem, vnhe, ekincm )

    USE kinds, ONLY: dbl
    USE control_flags, ONLY: tranp, trane, trhor, iprsta, tpre, tzeroc 
    USE control_flags, ONLY: tzerop, tzeroe, tfor, thdyn, lwf, tprnfor, tortho
    USE control_flags, ONLY: amprp, taurdr, ortho_eps, ortho_max
    USE ions_positions, ONLY: taus, tau0, tausm, vels, velsm, ions_hmove
    USE ions_base, ONLY: na, nsp, randpos, zv, ions_vel, pmass
    USE cell_base, ONLY: ainv, h, s_to_r, ibrav, omega, press, hold, r_to_s, deth
    USE cell_base, ONLY: wmass, iforceh, cell_force, cell_hmove
    USE efcalc, ONLY: ef_force
    USE elct, ONLY: n
    USE uspp, ONLY: betae => vkb, rhovan => becsum, deeq
    USE wavefunctions_module, ONLY: c0, cm, phi => cp
    USE io_global, ONLY: stdout
    USE cpr_subroutines
    USE wannier_subroutines
    USE core, ONLY: nlcc_any
    USE gvecw, ONLY: ngw
    USE reciprocal_vectors, ONLY: gstart
    USE wave_base, ONLY: wave_verlet
    USE cvan, ONLY: nvb
    USE ions_nose, ONLY: xnhp0,  xnhpm,  vnhp
    USE cell_nose, ONLY: xnhh0, xnhhm,  vnhh
    USE cp_electronic_mass, ONLY: emass, emaec => emass_cutoff

    COMPLEX(kind=8) :: eigr(:,:,:), ei1(:,:,:),  ei2(:,:,:),  ei3(:,:,:)
    COMPLEX(kind=8) :: eigrb(:,:,:)
    REAL(kind=8) :: bec(:,:), fion(:,:), becdr(:,:,:), fionm(:,:)
    REAL(kind=8) :: eself
    REAL(kind=8) :: taub(:,:)
    REAL(kind=8) :: b1(:), b2(:), b3(:)
    INTEGER :: irb(:,:,:)
    INTEGER :: nfi, iforce(:,:), nbeg
    LOGICAL :: tfirst
    COMPLEX(kind=8) :: sfac(:,:)
    COMPLEX(kind=8) :: rhog(:,:)
    REAL(kind=8) :: rhor(:,:), rhos(:,:), rhoc(:), enl, ekin
    REAL(kind=8) :: stress(:,:), detot(:,:), enthal, etot
    REAL(kind=8) :: lambda(:,:), lambdam(:,:), lambdap(:,:)
    REAL(kind=8) :: ema0bg(:)
    REAL(kind=8) :: dbec(:,:,:,:)
    REAL(kind=8) :: delt
    REAL(kind=8) :: bephi(:,:), becp(:,:)
    REAL(kind=8) :: velh(:,:)
    REAL(kind=8) :: dt2bye, xnhe0, xnhem, vnhe, ekincm


    REAL(kind=8), ALLOCATABLE :: emadt2(:), emaver(:)
    COMPLEX(kind=8), ALLOCATABLE :: c2(:), c3(:)
    REAL(kind=8) :: verl1, verl2
    REAL(kind=8) :: bigr
    INTEGER :: i, j, iter
    LOGICAL :: tlast = .FALSE.
    REAL(kind=8) :: fcell(3,3)


    ! We are restarting from file re compute ainv
    CALL invmat( 3, h, ainv, deth )

    IF( taurdr ) THEN
      ! Input positions read from input file and stored in tau0
      call r_to_s( tau0, taus, na, nsp, h )
    END IF

    IF( ANY( tranp( 1:nsp ) ) ) THEN
      ! Input positions are randomized
      call randpos(taus, na, nsp, tranp, amprp, ainv, iforce )
    END IF

    call s_to_r( taus, tau0, na, nsp, h )
!
    if( trane .and. trhor ) then
      call prefor( eigr, betae )
      call graham( betae, bec, c0 )
      cm(:, 1:n,1,1)=c0(:, 1:n,1,1)
    endif

    IF( tzeroc ) THEN
       hold = h
       velh = 0.0d0
    END IF

    fion = 0.0d0
    if( tzerop ) then
      tausm = taus
      vels  = 0.0d0
    end if

    if( tzeroe ) then
      lambdam = lambda
      WRITE( stdout, * ) 'Electronic velocities set to zero'  
    end if
!     
    call phfac( tau0, ei1, ei2, ei3, eigr )
    call strucf( ei1, ei2, ei3, sfac )
    call formf( tfirst, eself )
    call calbec ( 1, nsp, eigr, c0, bec )
    if (tpre) call caldbec( 1, nsp, eigr, c0 )


    if ( tzerop .or. tzeroe .or. taurdr ) then

      call initbox ( tau0, taub, irb )
      call phbox( taub, eigrb )
!
      if( lwf ) then
        call get_wannier_center( tfirst, c0, bec, becdr, eigr, eigrb, taub, irb, ibrav, b1, b2, b3 )
      end if
!
      call rhoofr ( nfi, c0, irb, eigrb, bec, rhovan, rhor, rhog, rhos, enl, ekin )
!
!     put core charge (if present) in rhoc(r)
!
      if ( nlcc_any ) call set_cc( irb, eigrb, rhoc )

      if( lwf ) then
        call write_charge_and_exit( rhog )
        call ef_tune( rhog, tau0 )
      end if

      call vofrho( nfi, rhor, rhog, rhos, rhoc, tfirst, tlast,                 &
     &  ei1, ei2, ei3, irb, eigrb, sfac, tau0, fion )

      call compute_stress( stress, detot, h, omega )

      if( lwf ) then
        call wf_options( tfirst, nfi, c0, rhovan, bec, becdr, eigr, eigrb, taub, irb, &
             ibrav, b1, b2, b3, rhor, rhog, rhos, enl, ekin  )
      end if

      call newd( rhor, irb, eigrb, rhovan, fion )
      call prefor( eigr, betae )
!

      if( tzeroe ) then

        verl1 = 1.0d0
        verl2 = 0.0d0

        ALLOCATE( emadt2( ngw ) )
        ALLOCATE( emaver( ngw ) )
        ALLOCATE( c2( ngw ) )
        ALLOCATE( c3( ngw ) )
        emadt2 = dt2bye * ema0bg
        emaver = emadt2 * 0.5d0

        if( lwf ) then
          call ef_potential( nfi, rhos, bec, deeq, betae, c0, cm, emadt2, emaver, verl1, verl2, c2, c3 )
        else
          do i = 1, n, 2
             call dforce(bec,betae,i,c0(1,i,1,1),c0(1,i+1,1,1),c2,c3,rhos)
             cm(:,i,1,1) = c0(:,i,1,1)
             CALL wave_verlet( cm(:, i  , 1, 1), c0(:, i  , 1, 1 ), &
                   verl1, verl2, emaver, c2 )
             CALL wave_verlet( cm(:, i+1, 1, 1), c0(:, i+1, 1, 1 ), &
                   verl1, verl2, emaver, c3 )
             if ( gstart == 2 ) then
                cm(1,  i,1,1)=cmplx(real(cm(1,  i,1,1)),0.0)
                cm(1,i+1,1,1)=cmplx(real(cm(1,i+1,1,1)),0.0)
             end if
          end do
        end if

        DEALLOCATE( emadt2 )
        DEALLOCATE( emaver )
        DEALLOCATE( c2 )
        DEALLOCATE( c3 )

      end if
!
!     nlfq needs deeq bec
!
      if( tfor .or. tprnfor ) call nlfq( c0, eigr, bec, becdr, fion )
!
      if( tfor .or. thdyn ) then
!
! interpolate new lambda at (t+dt) from lambda(t) and lambda(t-dt):
!
         lambdap(:,:) = 2.d0*lambda(:,:)-lambdam(:,:)
         lambdam(:,:)=lambda (:,:)
         lambda (:,:)=lambdap(:,:)

      endif
!
!     calphi calculates phi
!     the electron mass rises with g**2
!
      call calphi( c0, ema0bg, bec, betae, phi )
!
!     begin try and error loop (only one step!)
!
!       nlfl and nlfh need: lambda (guessed) becdr
!
      if ( tfor .or. tprnfor ) call nlfl( bec, becdr, lambda, fion )
      if( tpre ) then
         call nlfh( bec, dbec, lambda )
      endif

      if( tortho ) then
         call ortho( eigr, cm, phi, lambda, bigr, iter, dt2bye, ortho_eps, ortho_max, delt, bephi, becp )
      else
         call graham( betae, bec, cm )
      endif
!
      if ( tortho ) call updatc( dt2bye, lambda, phi, bephi, becp, bec, cm )
      call calbec ( nvb+1, nsp, eigr, cm, bec )
      if ( tpre ) call caldbec( 1, nsp, eigr, cm )
!
      if( thdyn ) then
        call cell_force( fcell, ainv, stress, omega, press, wmass )
        call cell_hmove( h, hold, delt, iforceh, fcell )
        call invmat( 3, h, ainv, deth )
      endif

      if( tfor ) then
        if( lwf ) then
          call ef_force( fion, na, nsp, zv )
        end if
        call ions_hmove( taus, tausm, iforce, pmass, fion, ainv, delt, na, nsp )
        CALL s_to_r( taus, tau0, na, nsp, h )
        call phfac( tau0, ei1, ei2, ei3, eigr )
        call calbec ( 1, nsp, eigr, c0, bec)
        if (tpre) call caldbec( 1, nsp, eigr, c0 )
      endif

      xnhp0=0.
      xnhpm=0.
      vnhp =0.
      fionm=0.
      CALL ions_vel( vels, taus, tausm, na, nsp, delt )

      xnhh0(:,:)=0.
      xnhhm(:,:)=0.
      vnhh (:,:) =0.
      velh (:,:)=(h(:,:)-hold(:,:))/delt
!
!     ======================================================
!     kinetic energy of the electrons
!     ======================================================

      call elec_fakekine( ekincm, ema0bg, emass, c0, cm, ngw, n, delt )

      xnhe0=0.
      xnhem=0.
      vnhe =0.

      call DSWAP( 2*ngw*n, c0, 1, cm, 1 )

      lambdam(:,:)=lambda(:,:)

    end if

    return
  end subroutine


!=----------------------------------------------------------------------------=!

      SUBROUTINE from_restart_fpmd( nfi, acc, gv, kp, ps, rhoe, desc, cm, c0, cdesc, &
          eigr, sfac, fi, ht_m, ht_0, atoms_m, atoms_0, fnl, vpot, edft)

!  this routine recreates the starting configuration from a restart file

! ... declare modules
      USE kinds, ONLY: dbl
      USE phase_factors_module, ONLY: strucf
      USE time_step, ONLY: delt
      USE reciprocal_space_mesh, ONLY: newg
      USE charge_density, ONLY: rhoofr
      USE wave_functions, ONLY: moveelect, gram, rande, fixwave
      USE wave_base, ONLY: wave_verlet
      USE electrons_module, ONLY: pmss,emass, nspin
      USE ions_base, ONLY: na, nsp, nax, randpos
      USE ions_module, ONLY: taui, cdmi, &
          set_reference_positions, print_scaled_positions, constraints_setup, set_velocities
      USE energies, ONLY: dft_energy_type
      USE cp_types, ONLY: recvecs, pseudo, phase_factors
      USE pseudopotential, ONLY: formf
      USE cell_module, ONLY: metric_print_info
      USE cell_module, only: boxdimensions, gethinv, alat
      USE cell_base, ONLY: r_to_s, s_to_r
      USE print_out_module, ONLY: printmain
      USE nl, ONLY: nlrh_m
      USE potentials, ONLY: vofrhos
      USE forces, ONLY: dforce_all
      USE orthogonalize, ONLY: ortho
      USE mp_global, ONLY: mpime, root, nproc, group
      USE io_global, ONLY: ionode, ionode_id
      USE io_global, ONLY: stdout
      USE mp, ONLY: mp_bcast
      USE brillouin, ONLY: kpoints
      USE wave_types, ONLY: wave_descriptor
      USE pseudo_projector, ONLY: projector
      USE control_flags, ONLY: tcarpar, nbeg, tranp, amprp, &
          tfor, tsdp, thdyn, tsdc, tbeg, tsde, tortho, prn, tzeroe, &
          tzerop, tzeroc, taurdr, tv0rd, nv0rd, trane, ampre, &
          force_pairing
      USE parameters, ONLY: nacx
      USE atoms_type_module, ONLY: atoms_type
      USE charge_types, ONLY: charge_descriptor
      USE ions_base, ONLY: vel_srt, tau_units

      IMPLICIT NONE

! ... declare subroutine arguments
      INTEGER   :: nfi
      REAL(dbl) :: acc(nacx)
      COMPLEX(dbl) :: sfac(:,:)
      TYPE (atoms_type) :: atoms_0, atoms_m
      TYPE (pseudo) :: ps
      TYPE (phase_factors) :: eigr
      TYPE (recvecs) :: gv
      TYPE (kpoints) :: kp
      COMPLEX(dbl), INTENT(INOUT) :: cm(:,:,:,:), c0(:,:,:,:)
      REAL(dbl) :: fi(:,:,:)
      TYPE (boxdimensions) :: ht_m, ht_0
      REAL (dbl) :: rhoe(:,:,:,:)
      TYPE (charge_descriptor) :: desc
      TYPE (wave_descriptor) :: cdesc
      TYPE (projector) :: fnl(:,:)
      REAL (dbl) :: vpot(:,:,:,:)
      TYPE (dft_energy_type) :: edft


! ... declare other variables

      LOGICAL, PARAMETER :: ttsde = .FALSE.
      REAL(dbl), PARAMETER :: svar1 = 1.d0
      REAL(dbl), PARAMETER :: svar2 = 0.d0

      REAL(dbl), ALLOCATABLE :: svar3(:)
      INTEGER ig, ib, i, j, k, ik, nb, is, ia, ierr, isa
      LOGICAL ttforce
      REAL(dbl)  timepre
      REAL(dbl) :: stau( 3 ), rtau( 3 ), hinv(3,3)

      COMPLEX (dbl), ALLOCATABLE :: eforce(:,:,:)

!  end of declarations
!  ----------------------------------------------

! ... if tbeg .eq. true ht_0 is not read from the restart file, and
! ... has been already been initialized in subroutine init1 togheter
! ... with the g vectors modules
! ... if tbeg is false, ht_0 is read from the restart file and now
! ... we have to compute the inverse and the volume of the cell,
! ... together with the new reciprocal vectors (gv)

      IF ( .NOT. tbeg ) THEN

        CALL newg( gv, kp, ht_0%m1 )

        CALL newgb( ht_0%hmat(:,1), ht_0%hmat(:,2), ht_0%hmat(:,3), ht_0%omega, alat )
!
        IF ( taurdr ) THEN

! ...       positions are read from stdin and not read from restart file, 
! ...       while the cell is read from the restart file, then real 
! ...       position do not correspond 

        END IF

      END IF

! ...   diagnostics
      IF(ionode) THEN

          WRITE( stdout,100)
          IF(.NOT.tbeg) THEN
            WRITE( stdout,110)
          ELSE
            WRITE( stdout,120)
          END IF
          IF(.NOT.taurdr) THEN
            WRITE( stdout,130)
          ELSE
            WRITE( stdout,140)
          END IF
          IF(tfor.AND.(.NOT.tsdp)) THEN
            IF( .NOT. tv0rd) THEN
              IF( .NOT. tzerop) THEN
                WRITE( stdout,150)
              ELSE
                WRITE( stdout,155)
              END IF
            ELSE
              WRITE( stdout,160)
            END IF
          END IF
          IF ( prn ) THEN
            CALL print_scaled_positions(atoms_0, 6, 'from restart module' )
          END IF
          IF(.NOT.tsde) THEN
            IF(.NOT.tzeroe) THEN
              WRITE( stdout,170)
            ELSE
              WRITE( stdout,180)
            END IF
          END IF
          WRITE( stdout,*)
      END IF

      IF(trane) THEN
          WRITE( stdout, 515) ampre
          CALL rande(c0, cdesc, ampre)
          CALL rande(cm, cdesc, ampre)
      END IF

      IF(tzeroc) THEN
          ht_m = ht_0
          ht_0%hvel = 0.0d0
      END IF

      IF( ANY(tranp) ) THEN
        hinv = TRANSPOSE( ht_0%m1 )
        CALL randpos( atoms_0%taus, atoms_0%na, atoms_0%nsp, tranp, amprp, hinv, atoms_0%mobile )
        CALL s_to_r( atoms_0%taus, atoms_0%taur, atoms_0%na, atoms_0%nsp, ht_0%hmat )
      END IF

      IF( tzerop .AND. tfor ) THEN
! ...     set initial velocities
          CALL set_velocities( atoms_m, atoms_0, vel_srt, ht_0, delt )
      END IF


! ...   computes form factors and initializes nl-pseudop. according
! ...   to starting cell (from ndr or again fort.10)
      CALL strucf(sfac, atoms_0, eigr, gv)
      CALL formf(ht_0, gv, kp, ps)

      IF( tzeroe .OR. tzerop ) THEN

! ...     set velocities to zero
! ...     set right initial conditions when c0=cm or stau0=staum
! ...     (the cell is kept fixed)

        ttforce = (tfor .OR. prn)


          atoms_0%for = 0.0d0
          edft%enl = nlrh_m(c0, cdesc, ttforce, atoms_0, fi, gv, kp, fnl, ps%wsg, ps%wnl, eigr)
          CALL rhoofr(gv, kp, c0, cdesc, fi, rhoe, desc, ht_0)
          CALL vofrhos(prn, prn, rhoe, desc, tfor, thdyn, ttforce, atoms_0, gv, kp, fnl, vpot, ps, &
            c0, cdesc, fi, eigr, sfac, timepre, ht_0, edft)

          IF( tzeroe ) THEN

            IF( tcarpar .AND. ( .NOT. force_pairing ) ) THEN

              ALLOCATE( svar3( SIZE( pmss ) ), STAT = ierr )
              IF( ierr /= 0 ) CALL errore(' restart ',' allocating svar3 ', ierr)

              svar3 = delt * delt / pmss * 0.5d0

              DO is = 1, cdesc%nspin

                ALLOCATE( eforce( SIZE( c0, 1 ), SIZE( c0, 2 ), SIZE( c0, 3 ) ), STAT=ierr )
                IF( ierr /= 0 ) CALL errore(' restart ',' allocating eforce ', ierr)

                CALL dforce_all( is, cm(:,:,:,is), cdesc, fi(:,:,is), eforce, gv, vpot(:,:,:,is), &
                   fnl(:,is), eigr, ps)

                DO ik = 1,  cdesc%nkl
                  DO ib = 1,  cdesc%nbl( is )
                    cm(:,ib,ik,is) = c0(:,ib,ik,is)
                    CALL wave_verlet( c0(:,ib,ik,is), cm(:,ib,ik,is), &
                      svar1, svar2, svar3, eforce(:,ib,ik) )
                  END DO
                  CALL fixwave(  is, c0(:,:,ik,is), cdesc, gv%kg_mask_l(:,ik) )
                END DO

                DEALLOCATE( eforce, STAT=ierr )
                IF( ierr /= 0 ) CALL errore(' restart ',' deallocating eforce ', ierr)

              END DO

              DEALLOCATE(svar3, STAT=ierr)
              IF( ierr /= 0 ) CALL errore(' restart ',' deallocating svar3 ', ierr)

              IF(tortho) THEN
                CALL ortho(cm, c0, cdesc, pmss, emass)
              ELSE
                CALL gram( c0, cdesc )
              END IF

            ELSE

              cm = c0

            END IF

          END IF

      END IF


! ... reset some variables if nbeg .LE. 0 (new simulation or step counter
! ... reset to 0)
      IF (nbeg .LE. 0) THEN
        acc   = 0.0d0
        nfi   = 0
        CALL set_reference_positions(cdmi, taui, atoms_0, ht_0)
      END IF

      CALL constraints_setup(ht_0, atoms_0)


  100 FORMAT( /,3X,'MD PARAMETERS READ FROM RESTART FILE',/ &
               ,3X,'------------------------------------')
  110 FORMAT(   3X,'Cell variables From RESTART file')
  120 FORMAT(   3X,'Cell variables From INPUT file')
  130 FORMAT(   3X,'Ions positions From RESTART file')
  140 FORMAT(   3X,'Ions positions From INPUT file')
  150 FORMAT(   3X,'Ions Velocities From RESTART file')
  155 FORMAT(   3X,'Ions Velocities set to ZERO')
  160 FORMAT(   3X,'Ions Velocities From STANDARD INPUT')
  170 FORMAT(   3X,'Electronic Velocities From RESTART file')
  180 FORMAT(   3X,'Electronic Velocities set to ZERO')
  515 FORMAT(   3X,'Initial random displacement of el. coordinates',/ &
                3X,'Amplitude = ',F10.6)

      RETURN
      END SUBROUTINE


END MODULE
