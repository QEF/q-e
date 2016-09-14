!
! Copyright (C) 2003-2009 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
PROGRAM phcg
  !-----------------------------------------------------------------------
  !
  USE ions_base,     ONLY: nat, tau
  USE io_global,     ONLY: ionode
  USE io_files,      ONLY: seqopn
  USE check_stop,    ONLY: check_stop_init
  USE mp_global,     ONLY: mp_startup, mp_global_end
  USE environment,   ONLY: environment_start
  ! The following instruction is just to make it clear that all modules
  ! from PWscf are needed sooner or later
  USE pwcom          
  USE cgcom

  IMPLICIT NONE

  REAL(DP), ALLOCATABLE :: dchi_dtau(:,:,:,:), dynout(:,:)
  REAL(DP), ALLOCATABLE :: w2(:)
  CHARACTER(len=9) :: cdate, ctime, code = 'PHCG'
  LOGICAL :: exst
  INTEGER :: i
  !
  CALL check_stop_init ()
  !
  ! Initialize MPI, clocks, print initial messages
  !
#if defined(__MPI)
  CALL mp_startup ( )
#endif
  CALL environment_start ( code )
  !
  CALL cg_readin
  !
  CALL cg_setup
  !
  !  calculate eps0, zstar, dynamical matrix for the unperturbed system
  !
  ALLOCATE  ( dynout(3*nat,3*nat))
  ALLOCATE  ( zstar( 3, 3, nat))
  ALLOCATE  ( w2( 3*nat))
  !
  CALL cg_eps0dyn(w2,dynout)
  !
  IF (raman) THEN
     IF (first==1.and.last==nmodes) THEN
        !
        !  calculate dX/dtau (X=polarizability) with finite differences
        !
        ALLOCATE ( dchi_dtau( 3, 3, 3, nat))
        CALL cg_dchi(dchi_dtau)
        !
        !  calculate nonresonant raman intensities for all modes
        !
        IF (trans) CALL raman_cs(dynout,dchi_dtau)
     ELSE
        !
        !  calculate nonresonant raman intensities for selected modes
        !
        CALL raman_cs2(w2,dynout)
     ENDIF
  ENDIF
  !
  CALL stop_clock('PHCG')
  CALL print_clock(' ')
  !
  !  close and delete temporary files, stop
  !
  IF (epsil .and. ionode) THEN
     iubar=1
     CALL seqopn (iubar,'filbar1','unformatted',exst)
     CLOSE(unit=iubar,status='delete')
     CALL seqopn (iubar,'filbar2','unformatted',exst)
     CLOSE(unit=iubar,status='delete')
     CALL seqopn (iubar,'filbar3','unformatted',exst)
     CLOSE(unit=iubar,status='delete')
     iudwf=10
     CALL seqopn (iudwf,'fildwx1','unformatted',exst)
     CLOSE(unit=iudwf,status='delete')
     CALL seqopn (iudwf,'fildwx2','unformatted',exst)
     CLOSE(unit=iudwf,status='delete')
     CALL seqopn (iudwf,'fildwx3','unformatted',exst)
     CLOSE(unit=iudwf,status='delete')
  ENDIF
  !
  CALL mp_global_end ()
  STOP
  !
9000 FORMAT (/5x,'Program ',a12,' starts ...',/5x,                     &
       &            'Today is ',a9,' at ',a9)
END PROGRAM phcg
!
!-----------------------------------------------------------------------
SUBROUTINE cg_dchi(dchi_dtau)
  !-----------------------------------------------------------------------
  !
  !  calculate dX/dtau with finite differences
  !
  USE constants,  ONLY : bohr_radius_angs
  USE ions_base,  ONLY : nat, tau
  USE io_global,  ONLY : stdout, ionode, ionode_id
  USE io_files,   ONLY : iunres, seqopn
  USE mp_world,   ONLY : world_comm
  USE mp,         ONLY : mp_bcast
  USE cell_base,  ONLY : omega, alat
  USE constants,  ONLY : fpi
  USE cgcom

  IMPLICIT NONE
  REAL(DP) :: dchi_dtau(3,3,3,nat)
  !
  REAL(DP) :: delta4(4), coeff4(4), delta2(2), coeff2(2), &
       delta, coeff, convfact
  INTEGER iudyn, nd, na, ipol, nd_, na_, ipol_, jpol, kpol
  LOGICAL :: exst
  DATA delta2/-1.d0, 1.d0/, coeff2/-0.5d0, 0.5d0/
  DATA delta4/-2.d0, -1.d0, 1.d0, 2.d0/
  DATA coeff4/ 0.08333333333333d0,-0.66666666666666d0,              &
       &       0.66666666666667d0,-0.08333333333337d0 /
  !
  CALL start_clock('cg_dchi')
  !
  !  Read partial results (if any)
  !
  na_  =1
  ipol_=1
  nd_  =1
  dchi_dtau(:,:,:,:) = 0.d0
  IF (recover) THEN
     IF (ionode) CALL seqopn( iunres, 'restart_d', 'FORMATTED', exst )
     CALL mp_bcast(exst,ionode_id,world_comm)
     IF ( .not. exst) GOTO 1
     IF (ionode) THEN
        READ(iunres,*,err=1,END=1) na_,ipol_,nd_
        READ(iunres,*,err=1,END=1) dchi_dtau
        CLOSE(unit=iunres)
     END IF
     CALL mp_bcast(na_,ionode_id,world_comm)
     CALL mp_bcast(ipol_,ionode_id,world_comm)
     CALL mp_bcast(nd,ionode_id,world_comm)
     CALL mp_bcast(dchi_dtau,ionode_id,world_comm)
     IF (na_<=nat) THEN
        WRITE( stdout,'(5x,"Restarting from atom ",i2,",  pol ",i1,      &
          &        ", nd=",i1)') na_,ipol_,nd_
     ELSE
        WRITE( stdout,'(5x,"Reading saved data")')
     ENDIF
     CLOSE(unit=iunres)
     GOTO 2
1    WRITE( stdout,'(/5x,"Restart failed, starting new calculation")')
     CLOSE(unit=iunres, status='delete' )
  ELSE
     WRITE( stdout,'(5x,"Starting calculation of Raman coefficients")')
  ENDIF
  !
2 CONTINUE
  !
  convfact = bohr_radius_angs**2
  !
  DO na=na_,nat
     DO ipol=1,3
        IF (na==na_.and.ipol<ipol_) GOTO 11
        DO nd=1,nderiv
           !
           !  Skip results from previous run (if any)
           !
           IF (na==na_.and.ipol==ipol_.and.nd<nd_) GOTO 12
           ! choose type of derivative formula (2- or 4-point)
           IF (nderiv==2) THEN
              delta=delta2(nd)
              coeff=coeff2(nd)
           ELSE
              delta=delta4(nd)
              coeff=coeff4(nd)
           ENDIF
           !
           ! Displaced atomic positions (remember that tau are in units of a0)
           !
           tau(ipol,na) =  tau(ipol,na) + delta*deltatau/alat
           !
           !
           CALL cg_neweps
           !
           tau(ipol,na) =  tau(ipol,na) - delta*deltatau/alat
           !
           ! convfact converts d chi/d tau to A^2 units
           !
           DO kpol=1,3
              DO jpol=1,3
                 dchi_dtau(kpol,jpol,ipol,na) =  &
                      dchi_dtau(kpol,jpol,ipol,na) + &
                      epsilon0(kpol,jpol)*coeff/deltatau * &
                      omega/fpi * convfact
              ENDDO
           ENDDO
           !
           !  Save partial results
           !
           !  parallel case: write only once !
           !
           IF ( ionode ) THEN
              !
              CALL seqopn( iunres, 'restart_d', 'FORMATTED', exst )
              IF (nd/=nderiv) THEN
                 WRITE(iunres,*) na,ipol,nd+1
              ELSEIF(ipol/=3) THEN
                 WRITE(iunres,*) na,ipol+1,1
              ELSE
                 WRITE(iunres,*) na+1,1,1
              ENDIF
              WRITE(iunres,*) dchi_dtau
              CLOSE(unit=iunres)
              !
           ENDIF

12         CONTINUE
        ENDDO
11      CONTINUE
     ENDDO
  ENDDO
  !
  WRITE( stdout,'(/5x, "Raman tensor (A^2)"/)')
  !
  DO na=1,nat
     DO ipol=1,3
        WRITE( stdout,'(/5x,"D X(i,j)",7x,3e14.6     &
             &    /5x,"----------  =  ",3e14.6   &
             &    /5x,"D tau(",i2,")_",i1,4x,3e14.6)')             &
             &  (( dchi_dtau(kpol,jpol,ipol,na), jpol=1,3), kpol=1,2),&
             &     na,ipol, (dchi_dtau(3,jpol,ipol,na), jpol=1,3)
     ENDDO
  ENDDO
  WRITE( stdout,*)
  !
  IF (ionode) THEN
     iudyn = 20
     INQUIRE (FILE = fildyn, EXIST = exst)
     IF (exst) THEN
        OPEN( unit=iudyn, file=fildyn, form='formatted', status='old', &
          position='append')
     ELSE
        OPEN( unit=iudyn, file=fildyn, form='formatted', status='new')
     ENDIF
     WRITE (iudyn,'(/5x,"Raman: D X_{alpha,beta}/D tau_{s,gamma} (units: A^2)"/)')
     DO na=1,nat
        DO ipol=1,3
           WRITE (iudyn,'("atom # ",i4,"   pol. ",i2)') na, ipol
           WRITE (iudyn,'(3e24.12)') &
                ( (dchi_dtau(kpol,jpol,ipol,na), jpol=1,3), kpol=1,3)
        ENDDO
     ENDDO
     CLOSE (unit=iudyn)
  ENDIF
  !
  RETURN
END SUBROUTINE cg_dchi
!
!-----------------------------------------------------------------------
SUBROUTINE cg_eps0dyn(w2,dynout)
  !-----------------------------------------------------------------------
  !
  USE constants,  ONLY : bohr_radius_angs, fpi
  USE cell_base,  ONLY : at, bg, omega
  USE ions_base,  ONLY : nat, tau, ityp, amass
  USE io_global,  ONLY : stdout, ionode, ionode_id
  USE io_files,   ONLY : iunres, seqopn
  USE symm_base,  ONLY : nsym, s, invs, irt
  USE mp_world,   ONLY : world_comm
  USE mp,         ONLY : mp_bcast
  USE cgcom
  !
  IMPLICIT NONE
  !
  REAL(DP) :: w2(3*nat), dynout(3*nat,3*nat), chi(3,3)
  !
  LOGICAL :: exst
  INTEGER :: na, i,j,  nt, iudyn, mode_done
  !
  !   calculate linear response to macroscopic fields
  !
  IF (epsil) THEN
     !
     !   verify if already calculated
     !
     IF (recover) THEN
        IF (ionode) CALL seqopn( iunres, 'restart_e', 'FORMATTED', exst )
        CALL mp_bcast(exst,ionode_id,world_comm)
        IF (.not. exst) GOTO 1
        IF (ionode) THEN
           READ(iunres,*,END=1,err=1) epsilon0
           READ(iunres,*,END=1,err=1) zstar
           CLOSE(unit=iunres)
        END IF
        CALL mp_bcast(epsilon0,ionode_id,world_comm)
        CALL mp_bcast(zstar,ionode_id,world_comm)
        GOTO 2
        !
1       CLOSE(unit=iunres, status='delete')
     ENDIF
     !
     CALL macro
     !
     CALL solve_e
     !
     !   calculate the dielectric tensor and effective charges
     !
     CALL dielec(.true.)
     CALL generate_effective_charges  (nat, nsym, s, invs, irt, at, bg, &
          n_diff_sites, equiv_atoms, has_equivalent, zstar)
     !
     !   save on file results
     !
     IF ( ionode ) THEN
        !
        CALL seqopn( iunres, 'restart_e', 'FORMATTED', exst )
        WRITE(iunres,*) epsilon0
        WRITE(iunres,*) zstar
        CLOSE(unit=iunres)
        !
     ENDIF
     !
     CALL output_tau (.false., .false.)
     !
     DO i=1,3
        DO j=1,3
           IF (i == j) THEN
              chi(i,j) = (epsilon0(i,j)-1.0_dp)*3.0_dp*omega/fpi &
                         /(epsilon0(i,j)+2.0_dp)
           ELSE
              chi(i,j) = epsilon0(i,j)*omega/fpi
           ENDIF
        ENDDO
     ENDDO
     WRITE(stdout,'(/5x,"dielectric constant",20x,"polarizability (A^3)")')
     WRITE(stdout,'(3f10.6,5x,3e14.6)') ( (epsilon0(i,j), j=1,3), &
                           (chi(i,j)*bohr_radius_angs**3, j=1,3), i=1,3)
     WRITE( stdout,'(/5x,"z*(",i2,")",3f10.4,/11x,3f10.4/11x,3f10.4)') &
             (na, ((zstar(i,j,na),j=1,3),i=1,3), na=1,nat)
  ENDIF
  !
2 CONTINUE
  !
  !   calculate linear response to lattice distorsions
  !
  IF (trans) THEN
     !
     !   verify if already calculated
     !
     IF (recover) THEN
        IF (ionode) CALL seqopn( iunres, 'restartph', 'FORMATTED', exst )
        CALL mp_bcast(exst,ionode_id,world_comm)
        IF (.not. exst) GOTO 10
        IF (ionode) READ (iunres,*,err=10,END=10) mode_done
        CALL mp_bcast(mode_done,ionode_id,world_comm)
        IF (mode_done==nmodes+1) THEN
           IF (ionode) THEN
              READ(iunres,*,END=10,err=10) dyn
              READ(iunres,*,END=10,err=10) w2
              CLOSE(unit=iunres)
           ENDIF
           CALL mp_bcast(dyn,ionode_id,world_comm)
           CALL mp_bcast(w2,ionode_id,world_comm)
           GOTO 20
        ELSE
           IF (ionode) CLOSE(unit=iunres)
        ENDIF
        !
10      CLOSE(unit=iunres, status='delete')
     ENDIF
     !
     CALL solve_ph ( )
     !
     !   get the complete dynamical matrix from the irreducible part
     !
     IF (nmodes==3*nat) CALL generate_dynamical_matrix &
          (nat, nsym, s, invs, irt, at, bg, n_diff_sites, equiv_atoms, &
           has_equivalent,dyn)
     !
     !   impose asr on the dynamical matrix
     !
     IF (asr) CALL set_asr(nat,nasr,dyn)
     !
     ! diagonalize the dynamical matrix
     !
     CALL dyndiar(dyn,3*nat,nmodes,u,nat,ityp,amass,w2,dynout)
     !
     !  find new equilibrium positions (in the harmonic approximation)
     !         if (lforce)
     !     &        call equilib(nat,tau,force,nmodes,w2,dyn,3*nat,dtau,alat)
  ENDIF
  !
  IF ( ionode ) THEN
     !
     IF (trans) CALL writedyn ( )
     !
     CALL seqopn( iunres, 'restartph', 'FORMATTED', exst )
     WRITE(iunres,*) nmodes+1
     WRITE(iunres,*) dyn
     WRITE(iunres,*) w2
     CLOSE(unit=iunres)
     !
  ENDIF
  !
20 CONTINUE
  !
  RETURN
  !
END SUBROUTINE cg_eps0dyn
!
!-----------------------------------------------------------------------
SUBROUTINE cg_neweps
  !-----------------------------------------------------------------------
  !
  USE constants, ONLY : bohr_radius_angs, fpi
  USE io_global, ONLY : stdout
  USE cell_base, ONLY : omega
  USE ions_base, ONLY : nat, tau
  USE fft_base,  ONLY : dfftp
  USE scf,       ONLY : rho, rho_core
  USE lsda_mod,  ONLY : current_spin
  USE funct,     ONLY : dmxc
  USE cgcom
  !
  IMPLICIT NONE

  INTEGER :: i, j
  REAL(DP) :: rhotot, chi(3,3)
  !
  !  recalculate self-consistent potential etc
  !
  CALL newscf
  !
  !  new derivative of the xc potential
  !
  dmuxc(:) = 0.d0
  DO i = 1,dfftp%nnr
     rhotot = rho%of_r(i,current_spin)+rho_core(i)
     IF ( rhotot> 1.d-30 ) dmuxc(i)= dmxc( rhotot)
     IF ( rhotot<-1.d-30 ) dmuxc(i)=-dmxc(-rhotot)
  ENDDO
  !
  !  re-initialize data needed for gradient corrections
  !
  CALL cg_setupdgc
  !
  !   calculate linear response to macroscopic fields
  !
  CALL macro
  !
  CALL solve_e
  !
  CALL dielec(.false.)
  !
  CALL output_tau (.false., .false.)
  !
  DO i=1,3
     DO j=1,3
        IF (i == j) THEN
           chi(i,j) = (epsilon0(i,j)-1.0_dp)*3.0_dp*omega/fpi &
                      /(epsilon0(i,j)+2.0_dp)
        ELSE
           chi(i,j) = epsilon0(i,j)*omega/fpi
        ENDIF
     ENDDO
  ENDDO
  !
  WRITE(stdout,'(/5x,"dielectric constant",20x,"polarizability (A^3)")')
  WRITE(stdout,'(3f10.6,5x,3e14.6)') ( (epsilon0(i,j), j=1,3), &
                                       (chi(i,j),j=1,3), i=1,3)
  WRITE(stdout,*)
  !
END SUBROUTINE cg_neweps
!
!-----------------------------------------------------------------------
SUBROUTINE newscf
  !-----------------------------------------------------------------------
  !
  USE basis, ONLY: starting_wfc 
  USE cellmd,ONLY: lmovecell
  USE gvecs, ONLY: doublegrid
  USE wvfct, ONLY: btype
  USE klist, ONLY: nkstot
  USE wvfct, ONLY: nbnd, nbndx
  USE noncollin_module, ONLY: report
  USE symm_base,     ONLY : nsym
  USE io_files,      ONLY : iunwfc, input_drho, output_drho
  USE ldaU,          ONLY : lda_plus_u
  USE control_flags, ONLY : restart, io_level, lscf, iprint, &
                            david, max_cg_iter, &
                            isolve, tr2, ethr, mixing_beta, nmix, niter
  USE extrapolation, ONLY : extrapolate_charge
  !
  IMPLICIT NONE
  INTEGER :: iter
  !
  CALL start_clock('PWSCF')
  !
  !  set all kind of stuff needed by self-consistent (re-)calculation
  !
!  dft='Same as Before'
  restart  =.false.
  io_level = 0
  lscf=.true.
  lda_plus_u=.false.
  doublegrid=.false.
  lmovecell=.false.
  iprint=10000
  input_drho=' '
  output_drho=' '
  starting_wfc='file'
  report=1
  if ( .not. allocated (btype) ) then
     allocate( btype( nbnd, nkstot ) )
     btype(:,:) = 1
  end if
  !
  !  since we use only Gamma we don't need symmetries
  !
  nsym=1
  !
  ! these must be tuned for fast convergence
  !
  david = 4
  nbndx = max (nbndx, david*nbnd)
  max_cg_iter=20
  isolve=0
  tr2 =1.d-8
  ethr=1.d-8
  mixing_beta=0.7d0
  nmix=4
  niter=50
  !
  CALL openfil
  !
  CALL extrapolate_charge( 1 )
  CALL hinit1
  CALL electrons ( )
  !
  CLOSE(unit=iunwfc, status='keep')
  !
  CALL stop_clock('PWSCF')
  !
  RETURN
END SUBROUTINE newscf
!
!-----------------------------------------------------------------------
SUBROUTINE raman_cs(dynout,dchi_dtau)
  !-----------------------------------------------------------------------
  !
  !  calculate Raman cross section
  !
  USE kinds,     ONLY : DP
  USE constants, ONLY : amu_ry
  USE ions_base, ONLY : nat
  USE io_global, ONLY : stdout
  USE cgcom,     ONLY : nmodes
  !
  IMPLICIT NONE
  REAL(DP) :: dynout(3*nat,3*nat), dchi_dtau(3,3,3,nat)
  !
  INTEGER :: nu, na, ipol, jpol, lpol
  REAL(DP), ALLOCATABLE :: raman_activity(:,:,:)
  !
  ALLOCATE  ( raman_activity( 3, 3, nmodes))
  WRITE( stdout,'(/5x, "Raman tensor for mode nu : dX_{alpha,beta}/d nu"/)')
  DO nu=1,nmodes
     !
     DO jpol=1,3
        DO ipol=1,3
           raman_activity(ipol,jpol,nu) = 0.0d0
           DO na=1,nat
              DO lpol=1,3
                 raman_activity(ipol,jpol,nu) = raman_activity(ipol,jpol,nu) +&
                      dchi_dtau(ipol,jpol,lpol,na) * dynout((na-1)*3+lpol,nu)
              ENDDO
           ENDDO
        ENDDO
     ENDDO
     !
     !   conversion factor from (Ry au for mass)^(-1) to amu(-1)
     !
     WRITE( stdout,'(i5,3x,3e14.6,2(/8x,3e14.6))') &
             nu,( ( raman_activity(ipol,jpol,nu)*amu_ry,jpol=1,3), ipol=1,3)
  ENDDO
  DEALLOCATE(raman_activity)
  !
  RETURN
END SUBROUTINE raman_cs
!
!-----------------------------------------------------------------------
SUBROUTINE raman_cs2(w2,dynout)
  !-----------------------------------------------------------------------
  !
  !  calculate d X/d u  (u=phonon mode) with finite differences
  !
  USE constants,  ONLY : bohr_radius_angs, ry_to_thz, ry_to_cmm1, amu_ry,&
                         fpi
  USE ions_base,  ONLY : nat, tau
  USE io_global,  ONLY : stdout, ionode
  USE io_files,   ONLY : iunres, seqopn
  USE cell_base,  ONLY : omega, alat
  USE cgcom
  !
  IMPLICIT NONE
  !
  REAL(DP) :: dynout(3*nat,3*nat), w2(3*nat)
  !
  REAL(DP), ALLOCATABLE :: raman_activity(:,:,:), infrared(:)
  REAL(DP) :: delta4(4), coeff4(4), delta2(2), coeff2(2), &
       delta, norm, coeff, convfact
  LOGICAL :: exst
  INTEGER iudyn, nd, nu, nd_, nu_, na, ipol, jpol
  DATA delta2/-1.d0, 1.d0/, coeff2/-0.5d0, 0.5d0/
  DATA delta4/-2.d0, -1.d0, 1.d0, 2.d0/
  DATA coeff4/ 0.08333333333333d0,-0.66666666666666d0,              &
       &       0.66666666666667d0,-0.08333333333337d0 /
  REAL(8):: polar(3), freq, cmfac, irfac
  REAL(8):: alpha, beta2
  !
  CALL start_clock('raman_cs2')
  !
  !  Read partial results (if any)
  !
  ALLOCATE ( raman_activity( 3, 3, last-first+1))
  nu_=1
  nd_=1
  raman_activity(:,:,:) = 0.d0
  IF (recover) THEN
     CALL seqopn( iunres, 'restart_d', 'FORMATTED', exst )
     IF (.not. exst) GOTO 1
     READ(iunres,*,err=1,END=1) nu_,nd_
     READ(iunres,*,err=1,END=1) raman_activity
     CLOSE(unit=iunres)
     IF (nu_<=last) THEN
        WRITE( stdout,'(5x,"Restarting from mode ",i3,", nd=",i1)') &
             nu_,nd_
     ELSE
        WRITE( stdout,'(5x,"Reading saved data")')
     ENDIF
     CLOSE(unit=iunres)
     GOTO 2
1    WRITE( stdout,'(/5x,"Restart failed, starting new calculation")')
     CLOSE(unit=iunres)
  ELSE
     WRITE( stdout,'(5x,"Starting calculation of Raman coeficients")')
  ENDIF
  !
2 CONTINUE
  !
  !   conversion factor from bohr^2*(Ry au for mass)^(-1/2) to A^2 amu(-1/2)
  !
  convfact = bohr_radius_angs**2*sqrt(amu_ry)
  !
  DO nu=first,last
     IF (nu<nu_) GOTO 11
     !
     ! eigendisplacements are normalized as <u|M|u>=1
     ! we want to add a normalized eigendisplacement instead: <u|u>=1
     !
     norm = 0
     DO na=1,nat
        DO ipol=1,3
           norm = norm + dynout(3*(na-1)+ipol,nu)**2
        ENDDO
     ENDDO
     norm = sqrt(norm)
     !
     DO nd=1,nderiv
        !
        !  Skip results from previous run (if any)
        !
        IF (nu==nu_.and.nd<nd_) GOTO 12
        ! choose type of derivative formula (2- or 4-point)
        IF (nderiv==2) THEN
           delta=delta2(nd)
           coeff=coeff2(nd)
        ELSE
           delta=delta4(nd)
           coeff=coeff4(nd)
        ENDIF
        !
        ! Displaced atomic positions (remember that tau are in units of a0)
        !
        DO na=1,nat
           DO ipol=1,3
              tau(ipol,na) =  tau(ipol,na) + delta * deltatau * &
                   dynout(3*(na-1)+ipol,nu) / norm / alat
           ENDDO
        ENDDO
        !
        CALL cg_neweps
        !
        ! reset atomic positions to equilibrium value
        !
        DO na=1,nat
           DO ipol=1,3
              tau(ipol,na) =  tau(ipol,na) - delta * deltatau * &
                   dynout(3*(na-1)+ipol,nu) / norm / alat
           ENDDO
        ENDDO
        !
        ! calculate derivative, multiply by norm in order to have
        ! raman cross section for a mode normalized as <u|M|u>=1
        !
        DO ipol=1,3
           DO jpol=1,3
              raman_activity(ipol,jpol,nu-first+1) =  &
                   raman_activity(ipol,jpol,nu-first+1) + &
                   epsilon0(ipol,jpol)*coeff/deltatau * norm * &
                   omega/fpi * convfact
           ENDDO
        ENDDO
        !
        !  Save partial results
        !
        !  parallel case: write only once !
        !
        IF ( ionode ) THEN
           !
           CALL seqopn( iunres, 'restart_d', 'FORMATTED', exst )
           IF (nd/=nderiv) THEN
              WRITE(iunres,*) nu,nd+1
           ELSE
              WRITE(iunres,*) nu+1,1
           ENDIF
           WRITE(iunres,*) raman_activity
           CLOSE(unit=iunres)
           !
        ENDIF

12      CONTINUE
     ENDDO
11   CONTINUE
  ENDDO
  !
  WRITE( stdout,'(/5x, "Raman tensor dX_{alpha,beta}/dQ_nu (A^2/amu^1/2)"/)')
  DO nu=first,last
     WRITE( stdout,'(i5,3x,3e14.6,2(/8x,3e14.6))') &
          nu,( ( raman_activity(ipol,jpol,nu-first+1),jpol=1,3), ipol=1,3)
  ENDDO
  !
  !  derivatives of epsilon are translated into derivatives of molecular
  !  polarizabilities by assuming a Clausius-Mossotti behavior
  !  (for anisotropic systems epsilon is replaced by its trace)
  !
  cmfac = 3.d0/( 2.d0 + (epsilon0(1,1) + epsilon0(2,2) + epsilon0(3,3))/3.d0 )
  !
  !   conversion factor for IR cross sections from
  !   (Ry atomic units * e^2)  to  (Debye/A)^2/amu
  !   1 Ry mass unit = 2 * mass of one electron = 2 amu
  !   1 e = 4.80324x10^(-10) esu = 4.80324 Debye/A
  !     (1 Debye = 10^(-18) esu*cm = 0.2081928 e*A)
  !
  irfac = 4.80324d0**2/2.d0*amu_ry
  !
  ALLOCATE (infrared(3*nat))
  !
  DO nu = 1,3*nat
     DO ipol=1,3
        polar(ipol)=0.0d0
     ENDDO
     DO na=1,nat
        DO ipol=1,3
           DO jpol=1,3
              polar(ipol) = polar(ipol) +  &
                   zstar(ipol,jpol,na)*dynout((na-1)*3+jpol,nu)
           ENDDO
        ENDDO
     ENDDO
     !
     ! the factor two is e^2 in Ry atomic units
     !
     infrared(nu) = 2.d0*(polar(1)**2+polar(2)**2+polar(3)**2)*irfac
     !
  ENDDO
  !
  WRITE( stdout,'(/5x,"IR cross sections are in (D/A)^2/amu units")')
  WRITE( stdout,'(5x,"Raman cross sections are in A^4/amu units")')
  WRITE( stdout,'(5x,"multiply by",f9.6," for Clausius-Mossotti correction")')&
        cmfac**2
  WRITE( stdout,'(/"#  mode   [cm-1]     [THz]       IR      Raman")')
  !
  DO nu = 1,3*nat
     !
     freq = sqrt(abs(w2(nu)))
     IF (w2(nu)<0.0) freq = -freq
     !
     ! alpha, beta2: see PRB 54, 7830 (1996) and refs quoted therein
     !
     IF( nu >= first .and. nu<= last ) THEN
        nu_ = nu-first+1
        alpha = (raman_activity(1,1,nu_) + &
                 raman_activity(2,2,nu_) + &
                 raman_activity(3,3,nu_))/3.d0
        beta2 = ( (raman_activity(1,1,nu_) - raman_activity(2,2,nu_))**2 + &
                  (raman_activity(1,1,nu_) - raman_activity(3,3,nu_))**2 + &
                  (raman_activity(2,2,nu_) - raman_activity(3,3,nu_))**2 + &
                  6.d0 * (raman_activity(1,2,nu_)**2 + &
                          raman_activity(1,3,nu_)**2 + &
                          raman_activity(2,3,nu_)**2) )/2.d0
     ELSE
        alpha = 0
        beta2 = 0
     ENDIF
     WRITE( stdout,'(i5,f10.2,f12.4,2f10.4)') &
          nu, freq*ry_to_cmm1, freq*ry_to_thz, infrared(nu), &
          (45.d0*alpha**2 + 7.0d0*beta2)
  ENDDO
  !
  DEALLOCATE (infrared)
  DEALLOCATE (raman_activity)
  RETURN
END SUBROUTINE raman_cs2

