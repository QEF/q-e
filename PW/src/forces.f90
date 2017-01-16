!
! Copyright (C) 2001-2011 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
! TB
! included monopole related forces
!----------------------------------------------------------------------------
!
!----------------------------------------------------------------------------
SUBROUTINE forces()
  !----------------------------------------------------------------------------
  !
  ! ... This routine is a driver routine which computes the forces
  ! ... acting on the atoms. The complete expression of the forces
  ! ... contains four parts which are computed by different routines:
  !
  ! ...  a)  force_lc,    local contribution to the forces
  ! ...  b)  force_cc,    contribution due to NLCC
  ! ...  c)  force_ew,    contribution due to the electrostatic ewald term
  ! ...  d)  force_us,    contribution due to the non-local potential
  ! ...  e)  force_corr,  correction term for incomplete self-consistency
  ! ...  f)  force_hub,   contribution due to the Hubbard term
  ! ...  g)  force_london, semi-empirical correction for dispersion forces
  !
  !
  USE kinds,         ONLY : DP
  USE io_global,     ONLY : stdout
  USE cell_base,     ONLY : at, bg, alat, omega  
  USE ions_base,     ONLY : nat, ntyp => nsp, ityp, tau, zv, amass, extfor
  USE fft_base,      ONLY : dfftp
  USE gvect,         ONLY : ngm, gstart, ngl, nl, igtongl, g, gg, gcutm
  USE lsda_mod,      ONLY : nspin
  USE symme,         ONLY : symvector
  USE vlocal,        ONLY : strf, vloc
  USE force_mod,     ONLY : force, lforce, sumfor
  USE scf,           ONLY : rho
  USE ions_base,     ONLY : if_pos
  USE ldaU,          ONLY : lda_plus_u, U_projection
  USE extfield,      ONLY : tefield, forcefield, monopole, forcemono, relaxz
  USE control_flags, ONLY : gamma_only, remove_rigid_rot, textfor, &
                            iverbosity, llondon, lxdm, ts_vdw
  USE plugin_flags
  USE bp,            ONLY : lelfield, gdir, l3dstring, efield_cart, &
                            efield_cry,efield
  USE uspp,          ONLY : okvan
  USE martyna_tuckerman, ONLY: do_comp_mt, wg_corr_force
  USE london_module, ONLY : force_london
  USE xdm_module,    ONLY : force_xdm
  USE tsvdw_module,  ONLY : FtsvdW
  USE esm,           ONLY : do_comp_esm, esm_bc, esm_force_ew
  USE qmmm,          ONLY : qmmm_mode
  !
  IMPLICIT NONE
  !
  REAL(DP), ALLOCATABLE :: forcenl(:,:), &
                           forcelc(:,:), &
                           forcecc(:,:), &
                           forceion(:,:), &
                           force_disp(:,:),&
                           force_disp_xdm(:,:),&
                           force_mt(:,:), &
                           forcescc(:,:), &
                           forces_bp_efield(:,:), &
                           forceh(:,:)
    ! nonlocal, local, core-correction, ewald, scf correction terms, and hubbard
!
! aux is used to store a possible additional density
! now defined in real space
!
  COMPLEX(DP), ALLOCATABLE :: auxg(:), auxr(:)
!
  REAL(DP) :: sumscf, sum_mm
  REAL(DP), PARAMETER :: eps = 1.e-12_dp
  INTEGER  :: ipol, na
    ! counter on polarization
    ! counter on atoms
  !
  !
  CALL start_clock( 'forces' )
  !
  ALLOCATE( forcenl( 3, nat ), forcelc( 3, nat ), forcecc( 3, nat ), &
            forceh( 3, nat ), forceion( 3, nat ), forcescc( 3, nat ) )
  !    
  forcescc(:,:) = 0.D0
  forceh(:,:)   = 0.D0
  force (:,:)   = 0.D0
  !
  ! ... The nonlocal contribution is computed here
  !
  CALL force_us( forcenl )
  !
  ! ... The local contribution
  !
  CALL force_lc( nat, tau, ityp, alat, omega, ngm, ngl, igtongl, &
                 g, rho%of_r, nl, nspin, gstart, gamma_only, vloc, &
                 forcelc )
  !
  ! ... The NLCC contribution
  !
  CALL force_cc( forcecc )
  !
  ! ... The Hubbard contribution
  !     (included by force_us if using beta as local projectors)
  !
  IF ( lda_plus_u .AND. U_projection.NE.'pseudo' ) CALL force_hub( forceh )
  !
  ! ... The ionic contribution is computed here
  !
  IF( do_comp_esm ) THEN
     CALL esm_force_ew( forceion )
  ELSE
     CALL force_ew( alat, nat, ntyp, ityp, zv, at, bg, tau, omega, g, &
                    gg, ngm, gstart, gamma_only, gcutm, strf, forceion )
  END IF
  !
  ! ... the semi-empirical dispersion correction
  !
  IF ( llondon ) THEN
    !
    ALLOCATE ( force_disp ( 3 , nat ) )
    force_disp ( : , : ) = 0.0_DP
    force_disp = force_london( alat , nat , ityp , at , bg , tau )
    !
  END IF
  IF (lxdm) THEN
     ALLOCATE (force_disp_xdm(3,nat))
     force_disp_xdm = 0._dp
     force_disp_xdm = force_xdm(nat)
  end if
     
  !
  ! ... The SCF contribution
  !
  CALL force_corr( forcescc )
  !
  IF (do_comp_mt) THEN
    !
    ALLOCATE ( force_mt ( 3 , nat ) )
    CALL wg_corr_force( .true.,omega, nat, ntyp, ityp, ngm, g, tau, zv, strf, &
                        nspin, rho%of_g, force_mt )
  END IF
  !
  ! ... call void routine for user define/ plugin patches on internal forces
  !
  call plugin_int_forces()
  !
  ! Berry's phase electric field terms
  !
  if(lelfield) then
     ALLOCATE ( forces_bp_efield (3,nat) )
     forces_bp_efield(:,:)=0.d0
     if(.not.l3dstring) then
        if(okvan) call  forces_us_efield(forces_bp_efield,gdir,efield)
        call forces_ion_efield(forces_bp_efield,gdir,efield)
     else
        if(okvan)then
           do ipol=1,3
              call  forces_us_efield(forces_bp_efield,ipol,efield_cry(ipol))
           enddo
        endif
        do ipol=1,3
           call  forces_ion_efield(forces_bp_efield,ipol,efield_cart(ipol))
        enddo
     endif
  endif
  !
  ! ... here we sum all the contributions and compute the total force acting
  ! ... on the crystal
  !
  DO ipol = 1, 3
     !
     sumfor = 0.D0
     !
     DO na = 1, nat
        !
        force(ipol,na) = force(ipol,na)    + &
                         forcenl(ipol,na)  + &
                         forceion(ipol,na) + &
                         forcelc(ipol,na)  + &
                         forcecc(ipol,na)  + &
                         forceh(ipol,na)   + &
                         forcescc(ipol,na)
        !
        IF ( llondon ) force(ipol,na) = force(ipol,na) + force_disp(ipol,na)
        IF ( lxdm )    force(ipol,na) = force(ipol,na) + force_disp_xdm(ipol,na)
        ! factor 2 converts from Ha to Ry a.u.
        IF ( ts_vdw )  force(ipol,na) = force(ipol,na) + 2.0_dp*FtsvdW(ipol,na)
        IF ( tefield ) force(ipol,na) = force(ipol,na) + forcefield(ipol,na)
        IF ( monopole ) force(ipol,na) = force(ipol,na) + forcemono(ipol,na) ! TB
        IF (lelfield)  force(ipol,na) = force(ipol,na) + forces_bp_efield(ipol,na)
        IF (do_comp_mt)force(ipol,na) = force(ipol,na) + force_mt(ipol,na) 

        sumfor = sumfor + force(ipol,na)
        !
     END DO
     !
     !TB
     IF ((monopole.and.relaxz).AND.(ipol==3)) WRITE( stdout, '("Total force in z direction = 0 disabled")')
     !
     IF ( (do_comp_esm .and. ( esm_bc .ne. 'pbc' )).or.(monopole.and.relaxz) ) THEN
        !
        ! ... impose total force along xy = 0
        !
        DO na = 1, nat
           IF ( ipol .ne. 3) force(ipol,na) = force(ipol,na)  &
                                            - sumfor / DBLE ( nat )
        END DO
        !
     ELSE IF ( qmmm_mode < 0 ) THEN
        !
        ! ... impose total force = 0 except in a QM-MM calculation
        !
        DO na = 1, nat
           force(ipol,na) = force(ipol,na) - sumfor / DBLE( nat ) 
        END DO
        !
     ENDIF
     !
  END DO
  !
  ! ... resymmetrize (should not be needed, but ...)
  !
  CALL symvector ( nat, force )
  !
  IF ( remove_rigid_rot ) &
     CALL remove_tot_torque( nat, tau, amass(ityp(:)), force  )
  !
  IF( textfor ) force(:,:) = force(:,:) + extfor(:,:)
  !
  ! ... call void routine for user define/ plugin patches on external forces
  !
  CALL plugin_ext_forces()
  !
  ! ... write on output the forces
  !
  WRITE( stdout, '(/,5x,"Forces acting on atoms (cartesian axes, Ry/au):", / )')
  DO na = 1, nat
     WRITE( stdout, 9035) na, ityp(na), force(:,na)
  END DO
  !
  ! ... forces on fixed coordinates are set to zero ( C.S. 15/10/2003 )
  !
  force(:,:)    = force(:,:)    * DBLE( if_pos )
  forcescc(:,:) = forcescc(:,:) * DBLE( if_pos )
  !
  IF ( iverbosity > 0 ) THEN
     IF ( do_comp_mt ) THEN
        WRITE( stdout, '(5x,"The Martyna-Tuckerman correction term to forces")')
        DO na = 1, nat
           WRITE( stdout, 9035) na, ityp(na), ( force_mt(ipol,na), ipol = 1, 3 )
        END DO
     END IF
     !
     WRITE( stdout, '(5x,"The non-local contrib.  to forces")')
     DO na = 1, nat
        WRITE( stdout, 9035) na, ityp(na), ( forcenl(ipol,na), ipol = 1, 3 )
     END DO
     WRITE( stdout, '(5x,"The ionic contribution  to forces")')
     DO na = 1, nat
        WRITE( stdout, 9035) na, ityp(na), ( forceion(ipol,na), ipol = 1, 3 )
     END DO
     WRITE( stdout, '(5x,"The local contribution  to forces")')
     DO na = 1, nat
        WRITE( stdout, 9035) na, ityp(na), ( forcelc(ipol,na), ipol = 1, 3 )
     END DO
     WRITE( stdout, '(5x,"The core correction contribution to forces")')
     DO na = 1, nat
        WRITE( stdout, 9035) na, ityp(na), ( forcecc(ipol,na), ipol = 1, 3 )
     END DO
     WRITE( stdout, '(5x,"The Hubbard contrib.    to forces")')
     DO na = 1, nat
        WRITE( stdout, 9035) na, ityp(na), ( forceh(ipol,na), ipol = 1, 3 )
     END DO
     WRITE( stdout, '(5x,"The SCF correction term to forces")')
     DO na = 1, nat
        WRITE( stdout, 9035) na, ityp(na), ( forcescc(ipol,na), ipol = 1, 3 )
     END DO
     !
     IF ( llondon) THEN
        WRITE( stdout, '(/,5x,"Dispersion contribution to forces:")')
        DO na = 1, nat
           WRITE( stdout, 9035) na, ityp(na), (force_disp(ipol,na), ipol = 1, 3)
        END DO
     END IF
     !
     IF (lxdm) THEN
        WRITE( stdout, '(/,5x,"XDM contribution to forces:")')
        DO na = 1, nat
           WRITE( stdout, 9035) na, ityp(na), (force_disp_xdm(ipol,na), ipol = 1, 3)
        END DO
     END IF
     !
     IF ( ts_vdw) THEN
        WRITE( stdout, '(/,5x,"TS-VDW contribution to forces:")')
        DO na = 1, nat
           WRITE( stdout, 9035) na, ityp(na), (2.0d0*FtsvdW(ipol,na), ipol=1,3)
        END DO
     END IF
     !
     ! TB monopole forces
     IF ( monopole) THEN
        WRITE( stdout, '(/,5x,"Monopole contribution to forces:")')
        DO na = 1, nat
           WRITE( stdout, 9035) na, ityp(na), (forcemono(ipol,na), ipol = 1, 3)
        END DO
     END IF
     !
  END IF
  !
  sumfor = 0.D0
  sumscf = 0.D0
  !
  DO na = 1, nat
     !
     sumfor = sumfor + force(1,na)**2 + force(2,na)**2 + force(3,na)**2
     sumscf = sumscf + forcescc(1,na)**2 + forcescc(2,na)**2+ forcescc(3,na)**2
     !
  END DO
  !
  sumfor = SQRT( sumfor )
  sumscf = SQRT( sumscf )
  !
  WRITE( stdout, '(/5x,"Total force = ",F12.6,5X, &
              &  "Total SCF correction = ",F12.6)') sumfor, sumscf
  !
  IF ( llondon .AND. iverbosity > 0 ) THEN
     !
     sum_mm = 0.D0
     DO na = 1, nat
        sum_mm = sum_mm + &
                 force_disp(1,na)**2 + force_disp(2,na)**2 + force_disp(3,na)**2
     END DO
     sum_mm = SQRT( sum_mm )
     WRITE ( stdout, '(/,5x, "Total Dispersion Force = ",F12.6)') sum_mm
     !
  END IF
  !
  IF ( lxdm .AND. iverbosity > 0 ) THEN
     !
     sum_mm = 0.D0
     DO na = 1, nat
        sum_mm = sum_mm + &
                 force_disp_xdm(1,na)**2 + force_disp_xdm(2,na)**2 + force_disp_xdm(3,na)**2
     END DO
     sum_mm = SQRT( sum_mm )
     WRITE ( stdout, '(/,5x, "Total XDM Force = ",F12.6)') sum_mm
     !
  END IF
  !
  DEALLOCATE( forcenl, forcelc, forcecc, forceh, forceion, forcescc )
  IF ( llondon )  DEALLOCATE ( force_disp )
  IF ( lxdm ) DEALLOCATE( force_disp_xdm ) 
  IF ( lelfield ) DEALLOCATE ( forces_bp_efield )
  !
  lforce = .TRUE.
  !
  CALL stop_clock( 'forces' )
  !
  IF ( ( sumfor < 10.D0*sumscf ) .AND. ( sumfor > nat*eps ) ) &
  WRITE( stdout,'(5x,"SCF correction compared to forces is large: ", &
                   &  "reduce conv_thr to get better values")')
  !
  IF(ALLOCATED(force_mt))   DEALLOCATE( force_mt )

  RETURN
  !
9035 FORMAT(5X,'atom ',I4,' type ',I2,'   force = ',3F14.8)
  !
END SUBROUTINE forces
