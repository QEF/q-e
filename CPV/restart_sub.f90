MODULE restart_subroutines

  IMPLICIT NONE
  SAVE

CONTAINS

  SUBROUTINE fromscra_sub( sfac, eigr, ei1, ei2, ei3, bec, becdr, tfirst, eself, fion, &
      taub, irb, eigrb, b1, b2, b3, nfi, rhog, rhor, rhos, rhoc, enl, ekin, stress,  &
      detot, enthal, etot, lambda, lambdam, lambdap, ema0bg, dbec, delt,  &
      bephi, becp, velh, dt2bye, iforce, fionm, nbeg, xnhe0, xnhem, vnhe, ekincm )

    USE control_flags, ONLY: tranp, trane, trhor, iprsta, tpre, tzeroc 
    USE control_flags, ONLY: tzerop, tzeroe, tfor, thdyn, lwf, tprnfor, tortho
    USE control_flags, ONLY: amprp, taurdr, ampre, tsde, ortho_eps, ortho_max
    USE ions_positions, ONLY: taus, tau0, tausm, vels, velsm
    USE ions_base, ONLY: na, nsp, randpos, zv, ions_vel, pmass
    USE cell_base, ONLY: ainv, h, s_to_r, ibrav, omega, press, hold, r_to_s, deth
    USE cell_base, ONLY: wmass, iforceh, cell_force
    USE elct, ONLY: n
    USE uspp, ONLY: betae => vkb, rhovan => becsum, deeq
    USE wavefunctions_module, ONLY: c0, cm, phi => cp
    USE io_global, ONLY: stdout
    USE cpr_subroutines
    USE wannier_subroutines
    USE core, ONLY: nlcc_any
    USE gvecw, ONLY: ngw
    USE reciprocal_vectors, ONLY: gstart
    USE wave_base, ONLY: wave_steepest
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
    REAL(kind=8) :: fcell(3,3), ccc, dt2hbe

    dt2hbe = 0.5d0 * dt2bye

    ! Input positions read from input file and stored in tau0

    IF( taurdr ) THEN
      call r_to_s( tau0, taus, na, nsp, h )
    END IF

    IF( ANY( tranp( 1:nsp ) ) ) THEN
      call invmat( 3, h, ainv, deth )
      call randpos(taus, na, nsp, tranp, amprp, ainv, iforce )
      call s_to_r( taus, tau0, na, nsp, h )
    END IF

!
    call phfac( tau0, ei1, ei2, ei3, eigr )
!     
    call initbox ( tau0, taub, irb )
!     
    call phbox( taub, eigrb )
!     
    if( trane ) then
!       
!     random initialization
!     
      call randin(1,n,gstart,ngw,ampre,cm)
      
    else if( nbeg == -3 ) then
!       
!     gaussian initialization
!     
      call gausin(eigr,cm)
      
    end if

    call prefor( eigr, betae )   !     prefor calculates betae (used by graham)
    call graham( betae, bec, cm )

    if( iprsta .ge. 3 ) call dotcsc( eigr, cm )

    hold = h
    velh = 0.0d0
    fion = 0.0d0
    tausm = taus
    vels  = 0.0d0
    lambdam = lambda
!     
    call phfac( tau0, ei1, ei2, ei3, eigr )
    call strucf( ei1, ei2, ei3, sfac )
    call formf( tfirst, eself )
    call calbec ( 1, nsp, eigr, cm, bec )
    if (tpre) call caldbec( 1, nsp, eigr, cm )

    call initbox ( tau0, taub, irb )
    call phbox( taub, eigrb )
!
    call rhoofr ( nfi, cm, irb, eigrb, bec, rhovan, rhor, rhog, rhos, enl, ekin )
!
!     put core charge (if present) in rhoc(r)
!
    if ( nlcc_any ) call set_cc( irb, eigrb, rhoc )

    call vofrho( nfi, rhor(1,1), rhog(1,1), rhos(1,1), rhoc(1), tfirst, tlast,                 &
     &  ei1(1,1,1), ei2(1,1,1), ei3(1,1,1), irb(1,1,1), eigrb(1,1,1), sfac(1,1), & 
     &  tau0(1,1), fion(1,1) )

    call compute_stress( stress, detot, h, omega )

    if(iprsta.gt.2) call print_atomic_var( fion, na, nsp, ' fion ' )

    call newd( rhor, irb, eigrb, rhovan, fion )
    call prefor( eigr, betae )
!
    ALLOCATE( emadt2( ngw ) )
    ALLOCATE( emaver( ngw ) )
    ALLOCATE( c2( ngw ) )
    ALLOCATE( c3( ngw ) )

    ccc = dt2hbe
    if(tsde) ccc = dt2bye
    emadt2 = ccc * ema0bg
    emaver = emadt2

    do i = 1, n, 2
       call dforce(bec,betae,i,cm(1,i,1,1),cm(1,i+1,1,1),c2,c3,rhos)
       CALL wave_steepest( c0(:, i  , 1, 1), cm(:, i  , 1, 1 ), emaver, c2 )
       CALL wave_steepest( c0(:, i+1, 1, 1), cm(:, i+1, 1, 1 ), emaver, c3 )
       if ( gstart == 2 ) then
          c0(1,  i,1,1)=cmplx(real(c0(1,  i,1,1)),0.0)
          c0(1,i+1,1,1)=cmplx(real(c0(1,i+1,1,1)),0.0)
       end if
    end do

    DEALLOCATE( emadt2 )
    DEALLOCATE( emaver )
    DEALLOCATE( c2 )
    DEALLOCATE( c3 )
!
!     nlfq needs deeq bec
!
    if( tfor .or. tprnfor ) call nlfq( cm, eigr, bec, becdr, fion )
!
!     calphi calculates phi
!     the electron mass rises with g**2
!
    call calphi( cm, ema0bg, bec, betae, phi )
!

    if( tortho ) then
       call ortho( eigr, c0, phi, lambda, bigr, iter, ccc, ortho_eps, ortho_max, delt, bephi, becp )
    else
       call graham( betae, bec, c0 )
    endif
!
    if ( tfor .or. tprnfor ) call nlfl( bec, becdr, lambda, fion )

    if ( iprsta >= 3 ) call print_lambda( lambda, n, 9, ccc )

    if ( tpre ) call nlfh( bec, dbec, lambda )
!
    if ( tortho ) call updatc( ccc, lambda, phi, bephi, becp, bec, c0 )
    call calbec ( nvb+1, nsp, eigr, c0, bec )
    if ( tpre ) call caldbec( 1, nsp, eigr, cm )

    if(iprsta.ge.3) call dotcsc(eigr,c0)
!
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

    lambdam(:,:)=lambda(:,:)

    return
  end subroutine

!=----------------------------------------------------------------------------=!

  SUBROUTINE restart_sub( sfac, eigr, ei1, ei2, ei3, bec, becdr, tfirst, eself, fion, &
      taub, irb, eigrb, b1, b2, b3, nfi, rhog, rhor, rhos, rhoc, enl, ekin, stress,  &
      detot, enthal, etot, lambda, lambdam, lambdap, ema0bg, dbec, delt,  &
      bephi, becp, velh, dt2bye, iforce, fionm, nbeg, xnhe0, xnhem, vnhe, ekincm )

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

END MODULE

