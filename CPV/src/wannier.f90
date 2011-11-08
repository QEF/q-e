!
! Copyright (C) 2002-2005 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! ... wannier function dynamics and electric field
!                                            - M.S
!
!----------------------------------------------------------------------------
MODULE efcalc
  !----------------------------------------------------------------------------
  !
  USE kinds,        ONLY : DP
  USE io_global,    ONLY : stdout
  USE wannier_base, ONLY : wf_efield, wf_switch
  USE wannier_base, ONLY : efx0, efy0, efz0, efx1, efy1, efz1, sw_len
  !
  IMPLICIT NONE
  !
  REAL(DP)              :: efx, efy, efz, sw_step
  REAL(DP), ALLOCATABLE :: xdist(:), ydist(:), zdist(:)
  !
  CONTAINS
  !
  !--------------------------------------------------------------------------
  SUBROUTINE clear_nbeg( nbeg )
    !--------------------------------------------------------------------------
    !
    ! ... some more electric field stuff
    !                              - M.S
    !
    INTEGER, INTENT(INOUT) :: nbeg
    !
    !
    IF ( wf_efield ) THEN
       !
       IF ( wf_switch ) THEN
          !
          WRITE( stdout, '(/,5X,"!----------------------------------!")' )
          WRITE( stdout, '(  5X,"!                                  !")' )
          WRITE( stdout, '(  5X,"! ADIABATIC SWITCHING OF THE FIELD !")' )
          WRITE( stdout, '(  5X,"!                                  !")' )
          WRITE( stdout, '(  5X,"!----------------------------------!",/)' )
          !
          nbeg=0
          !
       END IF
       !
    END IF
    !
    RETURN
    !
  END SUBROUTINE clear_nbeg
  !
  !--------------------------------------------------------------------------
  SUBROUTINE ef_force( fion, na, nsp, zv )
    !--------------------------------------------------------------------------
    !
    ! ... Electric Feild for ions here
    !
    IMPLICIT NONE
    !
    REAL(DP) :: fion(:,:), zv(:)
    INTEGER        :: na(:), nsp
    INTEGER        :: is, ia, isa
    !
    IF ( wf_efield ) THEN
       !
       isa = 0
       !
       DO is =1, nsp
          !
          DO ia = 1, na(is)
             !
             isa = isa + 1
             !
             fion(1,isa) = fion(1,isa) + efx * zv(is)
             fion(2,isa) = fion(2,isa) + efy * zv(is)
             fion(3,isa) = fion(3,isa) + efz * zv(is)
             !
          END DO
          !
       END DO
       !
    END IF
    !
    RETURN
    !
  END SUBROUTINE ef_force
  !
  !
  SUBROUTINE deallocate_efcalc()
     IF( ALLOCATED( xdist ) ) DEALLOCATE( xdist )
     IF( ALLOCATED( ydist ) ) DEALLOCATE( ydist )
     IF( ALLOCATED( zdist ) ) DEALLOCATE( zdist )
  END SUBROUTINE deallocate_efcalc
  !
END MODULE efcalc
!
!--------------------------------------------------------------------------
MODULE tune
  !--------------------------------------------------------------------------
  !
  USE kinds, ONLY : DP
  !
  LOGICAL        :: shift
  INTEGER        :: npts, av0, av1, xdir, ydir, zdir, start
  REAL(DP) :: alpha, b
  !
END MODULE tune
!
!--------------------------------------------------------------------------
MODULE wannier_module
  !--------------------------------------------------------------------------
  !
  ! ... In the presence of an electric field every wannier state feels a 
  ! ... different potantial, which depends on the position of its center. 
  ! ... RHOS is read in as the charge density in subrouting vofrho and 
  ! ... overwritten to be the potential.
  ! ...                                                             -M.S
  !
  USE kinds, ONLY : DP
  !
  IMPLICIT NONE
  !
  SAVE
  !
  LOGICAL                        :: what1, wann_calc
  REAL(DP)                 :: wfx, wfy, wfz, ionx, iony, ionz
  REAL(DP),    ALLOCATABLE :: utwf(:,:)
  REAL(DP),    ALLOCATABLE :: wfc(:,:)
  REAL(DP),    ALLOCATABLE :: rhos1(:,:), rhos2(:,:)
  COMPLEX(DP), ALLOCATABLE :: rhogdum(:,:)
  !
  CONTAINS
  !
  !------------------------------------------------------------------------
  SUBROUTINE allocate_wannier( nbsp, nrxxs, nspin, ng )
    !------------------------------------------------------------------------
    !
    INTEGER, INTENT(in) :: nbsp, nrxxs, nspin, ng
    !
    ALLOCATE( utwf( nbsp, nbsp ) )
    ALLOCATE( wfc( 3, nbsp ) )
    ALLOCATE( rhos1( nrxxs, nspin) )
    ALLOCATE( rhos2( nrxxs, nspin) )
    ALLOCATE( rhogdum( ng, nspin ) )
    !
    RETURN
    !
  END SUBROUTINE allocate_wannier
  !
  !------------------------------------------------------------------------
  SUBROUTINE deallocate_wannier()
    !------------------------------------------------------------------------
    !
    IF ( ALLOCATED( utwf ) )    DEALLOCATE( utwf )
    IF ( ALLOCATED( wfc ) )     DEALLOCATE( wfc )
    IF ( ALLOCATED( rhos1 ) )   DEALLOCATE( rhos1 )
    IF ( ALLOCATED( rhos2 ) )   DEALLOCATE( rhos2 )
    IF ( ALLOCATED( rhogdum ) ) DEALLOCATE( rhogdum )
    !
    RETURN
    !
  END SUBROUTINE deallocate_wannier
  !
END MODULE wannier_module
!
!--------------------------------------------------------------------------
MODULE electric_field_module
  !--------------------------------------------------------------------------
  !
  ! ... 1 Volt / meter  = 1/(5.1412*1.e+11) a.u.
  !
  USE kinds, ONLY : DP
  !
  IMPLICIT NONE
  !
  SAVE
  !
  LOGICAL        :: field_tune, ft
  REAL(DP) :: efe_elec, efe_ion, prefactor, e_tuned(3)
  REAL(DP) :: tt(3), tt2(3)
  REAL(DP) :: par, alen, blen, clen, rel1(3), rel2(3)
  !
END MODULE electric_field_module
!
!--------------------------------------------------------------------------
MODULE wannier_subroutines
  !--------------------------------------------------------------------------
  !
  USE kinds,     ONLY : DP
  USE io_global, ONLY : stdout, ionode
  !
  IMPLICIT NONE
  SAVE
  !
  CONTAINS
  !
  !------------------------------------------------------------------------
  SUBROUTINE wannier_startup( ibrav, alat, a1, a2, a3, b1, b2, b3 )
    !------------------------------------------------------------------------
    !
    USE wannier_module,        ONLY : utwf
    USE efcalc,                ONLY : wf_efield, efx0, efy0, efz0, &
                                      efx1, efy1, efz1, wf_switch, sw_len
    USE wannier_base,          ONLY : calwf, wfsd, wfdt, maxwfdt, nsd, nit, &
                                      wf_q, wf_friction, nsteps
    USE printout_base,         ONLY : printout_base_name
    !
    IMPLICIT NONE
    !
    INTEGER        :: ibrav
    REAL(DP) :: a1(3), a2(3), a3(3)
    REAL(DP) :: b1(3), b2(3), b3(3)
    REAL(DP) :: alat
    CHARACTER(LEN=256) :: fname
    !
    INTEGER :: i
    !
    ! ... More Wannier and Field Initialization
    !
    IF (calwf.GT.1) THEN
       IF (calwf.EQ.3 .AND. ionode ) THEN
          WRITE( stdout, * ) "------------------------DYNAMICS IN THE WANNIER BASIS--------------------------"
          WRITE( stdout, * ) "                             DYNAMICS PARAMETERS "
          IF (wfsd == 1) THEN
             WRITE( stdout, 12125) wf_q
             WRITE( stdout, 12126) wfdt
             WRITE( stdout, 12124) wf_friction
             WRITE( stdout, * ) nsteps,"STEPS OF DAMPED MOLECULAR DYNAMICS FOR OPTIMIZATION OF THE SPREAD"
          ELSE IF (wfsd == 2) THEN
             WRITE( stdout, 12132) wfdt
             WRITE( stdout, 12133) maxwfdt
             WRITE( stdout, * ) nsd,"STEPS OF STEEPEST DESCENT FOR OPTIMIZATION OF THE SPREAD"
             WRITE( stdout, * ) nit-nsd,"STEPS OF CONJUGATE GRADIENT FOR OPTIMIZATION OF THE SPREAD"
          ELSE
             WRITE( stdout, * ) "USING JACOBI ROTATIONS FOR OPTIMIZATION OF THE SPREAD"
          END IF
          WRITE( stdout, * ) "AVERAGE WANNIER FUNCTION SPREAD WRITTEN TO     FORT.24"
          fname = printout_base_name( "spr" )
          WRITE( stdout, * ) "INDIVIDUAL WANNIER FUNCTION SPREAD WRITTEN TO  "//TRIM(fname)
          fname = printout_base_name( "wfc" )
          WRITE( stdout, * ) "WANNIER CENTERS WRITTEN TO                     "//TRIM(fname)
          WRITE( stdout, * ) "SOME PERTINENT RUN-TIME INFORMATION WRITTEN TO FORT.27"
          WRITE( stdout, * ) "-------------------------------------------------------------------------------"
          WRITE( stdout, * )
12124     FORMAT(' DAMPING COEFFICIENT USED FOR WANNIER FUNCTION SPREAD OPTIMIZATION = ',f10.7)
12125     FORMAT(' FICTITIOUS MASS PARAMETER USED FOR SPREAD OPTIMIZATION            = ',f7.1)
12126     FORMAT(' TIME STEP USED FOR DAMPED DYNAMICS                                = ',f10.7)
          !
12132     FORMAT(' SMALLEST TIMESTEP IN THE SD / CG DIRECTION FOR SPREAD OPTIMIZATION= ',f10.7)
12133     FORMAT(' LARGEST TIMESTEP IN THE SD / CG DIRECTION FOR SPREAD OPTIMIZATION = ',f10.7)
       END IF
       WRITE( stdout, * ) "wannier_startup IBRAV SELECTED:",ibrav
       !
       CALL recips( a1, a2, a3, b1, b2, b3 )
       b1 = b1 * alat
       b2 = b2 * alat
       b3 = b3 * alat
       !
       CALL wfunc_init( calwf, b1, b2, b3, ibrav)
       !
       WRITE( stdout, * )
       utwf=0.d0
       DO i=1, SIZE( utwf, 1 )
          utwf(i, i)=1.d0
       END DO
    END IF
    IF(wf_efield) THEN

       CALL grid_map

       IF( ionode ) THEN
          WRITE( stdout, * ) "GRID MAPPING DONE"
          WRITE( stdout, * ) "DYNAMICS IN THE PRESENCE OF AN EXTERNAL ELECTRIC FIELD"
          WRITE( stdout, * )
          WRITE( stdout, * ) "POLARIZATION CONTRIBUTION OUTPUT TO FORT.28 IN THE FOLLOWING FORMAT"
          WRITE( stdout, * )
          WRITE( stdout, * ) "EFX, EFY, EFZ, ELECTRIC ENTHALPY(ELECTRONIC), ELECTRIC ENTHALPY(IONIC)"
          WRITE( stdout, * )
          WRITE( stdout, '(" E0(x) = ",F10.7)' ) efx0
          WRITE( stdout, '(" E0(y) = ",F10.7)' ) efy0
          WRITE( stdout, '(" E0(z) = ",F10.7)' ) efz0
          WRITE( stdout, '(" E1(x) = ",F10.7)' ) efx1
          WRITE( stdout, '(" E1(y) = ",F10.7)' ) efy1
          WRITE( stdout, '(" E1(z) = ",F10.7)' ) efz1
          !
          IF ( wf_switch ) WRITE( stdout, 12127) sw_len
          !
          WRITE( stdout, * )
          !
       END IF
       !
12127  FORMAT(' FIELD WILL BE TURNED ON ADIBATICALLY OVER ',i5,' STEPS')
    END IF
    !
    RETURN
    !
  END SUBROUTINE wannier_startup
  !
  !--------------------------------------------------------------------------
  SUBROUTINE get_wannier_center( tfirst, cm, bec, eigr, &
                                 eigrb, taub, irb, ibrav, b1, b2, b3 )
    !--------------------------------------------------------------------------
    !
    USE efcalc,         ONLY: wf_efield  
    USE wannier_base,   ONLY: calwf, jwf
    USE wannier_module, ONLY: what1, wfc, utwf
    !
    IMPLICIT NONE
    !
    LOGICAL, INTENT(in) :: tfirst
    COMPLEX(DP)   :: cm(:,:)
    REAL(DP)      :: bec(:,:)
    COMPLEX(DP)   :: eigrb(:,:), eigr(:,:)
    INTEGER             :: irb(:,:)
    REAL(DP)      :: taub(:,:)
    INTEGER             :: ibrav
    REAL(DP)      :: b1(:), b2(:), b3(:)
    !
    ! ... Get Wannier centers for the first step if wf_efield=true
    !
    IF ( wf_efield ) THEN
       !
       IF ( tfirst ) THEN
          !
          what1 = .TRUE.
          !
          jwf = 1
          !
          CALL wf( calwf,cm, bec, eigr, eigrb, taub, irb, &
                   b1, b2, b3, utwf, what1, wfc, jwf, ibrav )
          !
          what1 = .FALSE.
          !
       END IF
    END IF
    !
    RETURN
    !
  END SUBROUTINE get_wannier_center
  !
  !--------------------------------------------------------------------------
  SUBROUTINE ef_tune( rhog, tau0 )
    !--------------------------------------------------------------------------
    !
    USE electric_field_module, ONLY: field_tune, e_tuned
    USE wannier_module, ONLY: rhogdum
    !
    IMPLICIT NONE
    !
    COMPLEX(DP) :: rhog(:,:)
    REAL(DP)    :: tau0(:,:)
    !
    ! ... Tune the Electric field
    !
    IF ( field_tune ) THEN
       !
       rhogdum = rhog
       !
       CALL macroscopic_average( rhogdum, tau0, e_tuned )
       !
    END IF
    !
    RETURN
    !
  END SUBROUTINE ef_tune
  !
  !--------------------------------------------------------------------------
  SUBROUTINE write_charge_and_exit( rhog )
    !--------------------------------------------------------------------------
    !
    USE wannier_base, ONLY : writev
    !
    IMPLICIT NONE
    !
    COMPLEX(DP) :: rhog(:,:)
    !
    ! ... Write chargedensity in g-space
    !
    IF ( writev ) THEN
       !
       CALL write_rho_g( rhog )
       !
       CALL stop_run( .TRUE. )
       !
    END IF
    !
    RETURN
    !
  END SUBROUTINE write_charge_and_exit
  !
  !--------------------------------------------------------------------------
  SUBROUTINE wf_options( tfirst, nfi, cm, rhovan, bec, dbec, eigr, eigrb, &
                         taub, irb, ibrav, b1, b2, b3, rhor, drhor, rhog, &
                         drhog ,rhos, enl, ekin  )
    !--------------------------------------------------------------------------
    !
    USE efcalc,         ONLY : wf_efield
    USE wannier_base,   ONLY : nwf, calwf, jwf, wffort, iplot, iwf
    USE wannier_module, ONLY : what1, wfc, utwf
    USE cp_interfaces,  ONLY : rhoofr
    USE dener,          ONLY : denl, dekin6
    !
    IMPLICIT NONE
    !
    LOGICAL, INTENT(IN) :: tfirst
    INTEGER             :: nfi
    COMPLEX(DP)   :: cm(:,:)
    REAL(DP)      :: bec(:,:)
    REAL(DP)      :: dbec(:,:,:,:)
    REAL(DP)      :: rhovan(:,:,:)
    COMPLEX(DP)   :: eigrb(:,:), eigr(:,:)
    INTEGER             :: irb(:,:)
    REAL(DP)      :: taub(:,:)
    INTEGER             :: ibrav
    REAL(DP)      :: b1(:), b2(:), b3(:)
    COMPLEX(DP)   :: rhog(:,:)
    COMPLEX(DP)   :: drhog(:,:,:,:)
    REAL(DP)      :: drhor(:,:,:,:), rhor(:,:), rhos(:,:)
    REAL(DP)      :: enl, ekin 
    !
    INTEGER :: i, j
    !
    !
    ! ... Wannier Function options            - M.S
    !
    jwf=1
    IF (calwf.EQ.1) THEN
       DO i=1, nwf
          iwf=iplot(i)
          j=wffort+i-1
          CALL rhoofr (nfi,cm, irb, eigrb,bec,dbec,rhovan,rhor,drhor,rhog,drhog,rhos,enl,denl,ekin,dekin6,.false.,j)
       END DO
       !
       CALL stop_run( .TRUE. )
       !
    END IF
    !
    IF ( calwf == 2 ) THEN
       !
       ! ... calculate the overlap matrix
       !
       jwf=1
       !
       CALL wf (calwf,cm,bec,eigr,eigrb,taub,irb,b1,b2,b3,utwf,what1,wfc,jwf,ibrav)
       !
       CALL stop_run( .TRUE. )
       !
    END IF
    !
    IF (calwf.EQ.5) THEN
       !
       jwf=iplot(1)
       CALL wf (calwf,cm,bec,eigr,eigrb,taub,irb,b1,b2,b3,utwf,what1,wfc,jwf,ibrav)
       !
       CALL stop_run( .TRUE. )
       !
    END IF
    !
    ! ... End Wannier Function options - M.S
    !
    RETURN
  END SUBROUTINE wf_options
  !
  !--------------------------------------------------------------------------
  SUBROUTINE ef_potential( nfi, rhos, bec, deeq, betae, c0, cm, emadt2, emaver, verl1, verl2 )
    !--------------------------------------------------------------------------
    !
    USE efcalc,                 ONLY : wf_efield, efx, efy, efz, &
                                       efx0, efy0, efz0, efx1, efy1, efz1, &
                                       wf_switch, sw_len, sw_step, xdist,  &
                                       ydist, zdist
    USE electric_field_module,  ONLY : field_tune, e_tuned, par, rel1, rel2
    USE wannier_module,         ONLY : rhos1, rhos2, wfc
    USE electrons_base,         ONLY : nbsp, nspin, nupdwn, f, ispin
    USE cell_base,              ONLY : ainv, alat, at
    USE gvect,                  ONLY : gstart
    USE control_flags,          ONLY : tsde
    USE wave_base,              ONLY : wave_steepest, wave_verlet
    USE cp_interfaces,          ONLY : dforce
    USE fft_base,               ONLY : dffts
    !
    IMPLICIT NONE
    !
    INTEGER :: nfi
    REAL(DP) :: rhos(:,:)
    REAL(DP) :: bec(:,:)
    REAL(DP) :: deeq(:,:,:,:)
    COMPLEX(DP) :: betae(:,:)
    COMPLEX(DP) :: c0( :, : )
    COMPLEX(DP) :: cm( :, : )
    REAL(DP) :: emadt2(:)
    REAL(DP) :: emaver(:)
    REAL(DP) :: verl1, verl2
    REAL(DP) :: a1(3), a2(3), a3(3)
    COMPLEX(DP), ALLOCATABLE :: c2( : ), c3( : )
    INTEGER :: i, ir
    !
    ! ... Potential for electric field
    !
    ALLOCATE( c2( SIZE( c0, 1 )))
    ALLOCATE( c3( SIZE( c0, 1 )))

    a1(:) = at(:,1)/alat ; a2(:) = at(:,2)/alat ; a3(:) = at(:,3)/alat

    IF(wf_efield) THEN
       IF(field_tune) THEN
          efx=e_tuned(1)
          efy=e_tuned(2)
          efz=e_tuned(3)
          WRITE( stdout, '(" wf_efield Now ",3(F12.8,1X))' ) efx, efy,efz
          !
       ELSE
          IF(wf_switch) THEN
             par=0.d0
             IF(nfi.LE.sw_len) THEN
                sw_step=1.0d0/DBLE(sw_len)
                par=nfi*sw_step
                IF(efx1.LT.efx0) THEN
                   efx=efx0-(efx0-efx1)*par**5*(70*par**4-315*par**3+540*par**2-420*par+126)
                ELSE
                   efx=efx0+(efx1-efx0)*par**5*(70*par**4-315*par**3+540*par**2-420*par+126)
                END IF
                IF(efy1.LT.efy0) THEN
                   efy=efy0-(efy0-efy1)*par**5*(70*par**4-315*par**3+540*par**2-420*par+126)
                ELSE
                   efy=efy0+(efy1-efy0)*par**5*(70*par**4-315*par**3+540*par**2-420*par+126)
                END IF
                IF(efz1.LT.efz0) THEN
                   efz=efz0-(efz0-efz1)*par**5*(70*par**4-315*par**3+540*par**2-420*par+126)
                ELSE
                   efz=efz0+(efz1-efz0)*par**5*(70*par**4-315*par**3+540*par**2-420*par+126)
                END IF
             END IF
          ELSE
             efx=efx1
             efy=efy1
             efz=efz1
          END IF
       END IF
    END IF
    DO i=1,nbsp,2
       IF(wf_efield) THEN
          rhos1=0.d0
          rhos2=0.d0
          DO ir=1,dffts%nnr
             rel1(1)=xdist(ir)*a1(1)+ydist(ir)*a2(1)+zdist(ir)*a3(1)-wfc(1,i)
             rel1(2)=xdist(ir)*a1(2)+ydist(ir)*a2(2)+zdist(ir)*a3(2)-wfc(2,i)
             rel1(3)=xdist(ir)*a1(3)+ydist(ir)*a2(3)+zdist(ir)*a3(3)-wfc(3,i)
             !  minimum image convention
             CALL pbc(rel1,a1,a2,a3,ainv,rel1)
             IF(nspin.EQ.2) THEN
                IF(i.LE.nupdwn(1)) THEN
                   rhos1(ir,1)=rhos(ir,1)+efx*rel1(1)+efy*rel1(2)+efz*rel1(3)
                ELSE
                   rhos1(ir,2)=rhos(ir,2)+efx*rel1(1)+efy*rel1(2)+efz*rel1(3)
                END IF
             ELSE
                rhos1(ir,1)=rhos(ir,1)+efx*rel1(1)+efy*rel1(2)+efz*rel1(3)
             END IF
             IF(i.NE.nbsp) THEN
                rel2(1)=xdist(ir)*a1(1)+ydist(ir)*a2(1)+zdist(ir)*a3(1)-wfc(1,i+1)
                rel2(2)=xdist(ir)*a1(2)+ydist(ir)*a2(2)+zdist(ir)*a3(2)-wfc(2,i+1)
                rel2(3)=xdist(ir)*a1(3)+ydist(ir)*a2(3)+zdist(ir)*a3(3)-wfc(3,i+1)
                !  minimum image convention
                CALL pbc(rel2,a1,a2,a3,ainv,rel2)
                IF(nspin.EQ.2) THEN
                   IF(i+1.LE.nupdwn(1)) THEN
                      rhos2(ir,1)=rhos(ir,1)+efx*rel2(1)+efy*rel2(2)+efz*rel2(3)
                   ELSE
                      rhos2(ir,2)=rhos(ir,2)+efx*rel2(1)+efy*rel2(2)+efz*rel2(3)
                   END IF
                ELSE
                   rhos2(ir,1)=rhos(ir,1)+efx*rel2(1)+efy*rel2(2)+efz*rel2(3)
                END IF
             ELSE
                rhos2(ir,:)=rhos1(ir,:)
             END IF
          END DO
          CALL dforce(i,bec,betae,c0,c2,c3,rhos1,dffts%nnr,ispin,f,nbsp,nspin,rhos2)
       ELSE
          CALL dforce(i,bec,betae,c0,c2,c3,rhos,dffts%nnr,ispin,f,nbsp,nspin)
       END IF
       IF(tsde) THEN
          CALL wave_steepest( cm(:, i  ), c0(:, i  ), emadt2, c2 )
          CALL wave_steepest( cm(:, i+1), c0(:, i+1), emadt2, c3 )
       ELSE
          CALL wave_verlet( cm(:, i  ), c0(:, i  ), verl1, verl2, emaver, c2 )
          CALL wave_verlet( cm(:, i+1), c0(:, i+1), verl1, verl2, emaver, c3 )
       ENDIF
       IF (gstart.EQ.2) THEN
          cm(1,  i)=CMPLX(DBLE(cm(1,  i)),0.d0,kind=DP)
          cm(1,i+1)=CMPLX(DBLE(cm(1,i+1)),0.d0,kind=DP)
       END IF
    END DO

    DEALLOCATE( c2 )
    DEALLOCATE( c3 )

    RETURN
  END SUBROUTINE ef_potential
  !
  !--------------------------------------------------------------------
  ! ... Electric Field Implementation for Electric Enthalpy
  ! ...                                              - M.S
  !--------------------------------------------------------------------
  !
  !--------------------------------------------------------------------------
  SUBROUTINE ef_enthalpy( enthal, tau0 )
    !--------------------------------------------------------------------------
    !
    USE efcalc,                ONLY : wf_efield, efx, efy, efz
    USE electric_field_module, ONLY : efe_elec, efe_ion, tt2, tt
    USE wannier_module,        ONLY : wfx, wfy, wfz, ionx, iony, ionz, wfc
    USE electrons_base,        ONLY : nbsp, f
    USE cell_base,             ONLY : ainv, alat, at
    USE ions_base,             ONLY : na, nsp, zv
    USE io_global,             ONLY : ionode
    !
    IMPLICIT NONE
    !
    REAL(DP) :: enthal, tau0(:,:)
    REAL(DP) :: a1(3), a2(3), a3(3)
    INTEGER        :: i, is, ia, isa
    !
    a1(:) = at(:,1)/alat ; a2(:) = at(:,2)/alat ; a3(:) = at(:,3)/alat
    IF(wf_efield) THEN
       !  Electronic Contribution First
       wfx=0.d0
       wfy=0.d0
       wfz=0.d0
       efe_elec=0.d0
       DO i=1,nbsp
          tt2(1)=wfc(1,i)
          tt2(2)=wfc(2,i)
          tt2(3)=wfc(3,i)
          CALL pbc(tt2,a1,a2,a3,ainv,tt2)
          wfx=wfx+f(i)*tt2(1)
          wfy=wfy+f(i)*tt2(2)
          wfz=wfz+f(i)*tt2(3)
       END DO
       efe_elec=efe_elec+efx*wfx+efy*wfy+efz*wfz
       !Then Ionic Contribution
       ionx=0.d0
       iony=0.d0
       ionz=0.d0
       efe_ion=0.d0
       isa = 0
       DO is=1,nsp
          DO ia=1,na(is)
             isa = isa + 1
             tt(1)=tau0(1,isa)
             tt(2)=tau0(2,isa)
             tt(3)=tau0(3,isa)
             CALL pbc(tt,a1,a2,a3,ainv,tt)
             ionx=ionx+zv(is)*tt(1)
             iony=iony+zv(is)*tt(2)
             ionz=ionz+zv(is)*tt(3)
          END DO
       END DO
       efe_ion=efe_ion+efx*ionx+efy*iony+efz*ionz
       IF( ionode ) THEN
          WRITE(28,'(f12.9,1x,f12.9,1x,f12.9,1x,f20.15,1x,f20.15)') efx, efy, efz, efe_elec,-efe_ion
       END IF
    ELSE
       efe_elec = 0.0_DP      
       efe_ion  = 0.0_DP      
    END IF
    enthal=enthal+efe_elec-efe_ion

    RETURN
  END SUBROUTINE ef_enthalpy
  !
  !--------------------------------------------------------------------------
  SUBROUTINE wf_closing_options( nfi, c0, cm, bec, eigr, eigrb, taub,  &
                                 irb, ibrav, b1, b2, b3, taus, tausm, vels,   &
                                 velsm, acc, lambda, lambdam, descla, xnhe0, xnhem,   &
                                 vnhe, xnhp0, xnhpm, vnhp, nhpcl,nhpdim,ekincm,&
                                 xnhh0, xnhhm, vnhh, velh, ecut, ecutw, delt, &
                                 celldm, fion, tps, mat_z, occ_f, rho )
    !--------------------------------------------------------------------------
    !
    USE efcalc,         ONLY : wf_efield
    USE wannier_base,   ONLY : nwf, calwf, jwf, wffort, iplot, iwf
    USE wannier_module, ONLY : what1, wfc, utwf
    USE electrons_base, ONLY : nbsp
    USE gvecw,          ONLY : ngw
    USE control_flags,  ONLY : ndw
    USE cell_base,      ONLY : h, hold
    USE uspp_param,     ONLY : nvb
    USE cp_interfaces,  ONLY : writefile
    USE descriptors,    ONLY : la_descriptor
    !
    IMPLICIT NONE
    !
    INTEGER           :: nfi
    COMPLEX(DP) :: c0(:,:)
    COMPLEX(DP) :: cm(:,:)
    REAL(DP)    :: bec(:,:)
    COMPLEX(DP) :: eigrb(:,:), eigr(:,:)
    INTEGER           :: irb(:,:)
    REAL(DP)    :: taub(:,:)
    INTEGER           :: ibrav
    REAL(DP)    :: b1(:), b2(:), b3(:)
    REAL(DP)    :: taus(:,:), tausm(:,:), vels(:,:), velsm(:,:)
    REAL(DP)    :: acc(:)
    REAL(DP)    :: lambda(:,:,:), lambdam(:,:,:)
    TYPE(la_descriptor), INTENT(IN) :: descla(:)
    REAL(DP)    :: xnhe0, xnhem, vnhe, xnhp0(:), xnhpm(:), vnhp(:), ekincm
    INTEGER           :: nhpcl, nhpdim
    REAL(DP)    :: velh(:,:)
    REAL(DP)    :: xnhh0(:,:), xnhhm(:,:), vnhh(:,:)
    REAL(DP)    :: ecut, ecutw, delt, celldm(:)
    REAL(DP)    :: fion(:,:), tps
    REAL(DP)    :: mat_z(:,:,:), occ_f(:), rho(:,:)
    !
    CALL start_clock('wf_close_opt')
    !
    ! ... More Wannier Function Options
    !
    IF ( calwf == 4 ) THEN
       !
       jwf = 1
       !
       CALL wf( calwf, c0, bec, eigr, eigrb, taub, irb, &
                b1, b2, b3, utwf, what1, wfc, jwf, ibrav )
       !
       IF ( nvb == 0 ) THEN
          !
          CALL wf( calwf, cm, bec, eigr, eigrb, taub, irb, &
                   b1, b2, b3, utwf, what1, wfc, jwf, ibrav )
          !
       ELSE
          !
          cm = c0
          !
       END IF
       !
       CALL writefile( h, hold, nfi, c0, cm, taus, &
                       tausm, vels, velsm,acc, lambda, lambdam, descla, xnhe0, xnhem, &
                       vnhe, xnhp0, xnhpm, vnhp,nhpcl,nhpdim,ekincm, xnhh0, xnhhm,&
                       vnhh, velh, fion, tps, mat_z, occ_f, rho )
       !
       CALL stop_clock('wf_close_opt')
       CALL stop_run( .TRUE. )
       !
    END IF
    !
    IF ( calwf == 3 ) THEN
       !
       ! ... construct overlap matrix and calculate spreads and do Localization
       !
       jwf = 1
       !
       CALL wf( calwf, c0, bec, eigr, eigrb, taub, irb, &
                b1, b2, b3, utwf, what1, wfc, jwf, ibrav )
       !
       CALL stop_clock('wf_close_opt')
       !
    END IF
    !
    RETURN
    !
  END SUBROUTINE wf_closing_options
  !
END MODULE wannier_subroutines
