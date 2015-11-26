!
! Copyright (C) 2001-2011 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------------
SUBROUTINE vcsmd( conv_ions )
  !----------------------------------------------------------------------------
  !
  ! Main (interface) routine between PWSCF and the variable-cell shape
  ! molecular dynamics code by R.M. Wentzcovitch, PRB 44, 2358 (1991).
  !
  ! Molecular and/or cell dynamics is performed according to the value of
  ! the switch variable calc:
  !
  !  calc  = 'md'   : standard molecular dynamics
  !  calc  = 'mm'   : structural minimization by damped dynamics
  !  calc  = 'cd'   : Parrinello-Rahman cell dynamics
  !  calc  = 'cm'   : Parrinello-Rahman cell minimization by damped dynamics
  !  calc  = 'nd'   : Wentzcovitch's new cell dynamics
  !  calc  = 'nm'   : Wentzcovitch's new cell minimization by damped dynamics
  !
  ! Dynamics performed using Beeman algorithm, J. Comp. Phys. 20, 130 (1976))
  !
  ! Contraints with vcsmd have been implemented by Vivek Ranjan in 2012
  ! from the Department of Physics, North Carolina State University
  ! Raleigh, North Carolina, USA
  !
  USE kinds,           ONLY : DP
  USE io_global,       ONLY : stdout
  USE constants,       ONLY : e2, ry_kbar, amu_ry
  USE cell_base,       ONLY : omega, alat, at, bg, iforceh, fix_volume, fix_area
  USE ions_base,       ONLY : tau, nat, ntyp => nsp, ityp, atm, if_pos
  USE cellmd,          ONLY : nzero, ntimes, calc, press, at_old, omega_old, &
                              cmass, ntcheck, lmovecell
  USE dynamics_module, ONLY : dt, temperature
  USE ions_base,       ONLY : amass, if_pos 
  USE relax,           ONLY : epse, epsf, epsp
  USE force_mod,       ONLY : force, sigma
  USE control_flags,   ONLY : nstep, istep, tolp, lconstrain
  USE parameters,      ONLY : ntypx
  USE ener,            ONLY : etot
  USE io_files,        ONLY : prefix, delete_if_present, seqopn

  USE constraints_module, ONLY : nconstr
  USE constraints_module, ONLY : remove_constr_force, check_constraint  
  !
  !
  IMPLICIT NONE
  LOGICAL, INTENT (OUT) :: conv_ions
  !
  ! ... I/O variable first
  !
  ! PWSCF variables
  !  nat  = total number of atoms
  !  ntyp = total number of atomic types
  !  ityp(na)  = atomic type for na-th atom
  !  tau(i,na) =  position of the na-th atom
  !  at (icar,ivec) = direct Bravais lattice vectors
  !  bg (icar,ivec) = reciprocal lattice vectors
  !  amass_(nt) = mass (in atomic ryd units) for atom of nt-th type
  !  cmass = cell mass in ryd units.
  !  press = target pressure in ryd/(a.u.)^3
  !
  ! ... local variables
#if ! defined (__REDUCE_OUTPUT)
! for vcsmd with constraints
  REAL(DP), EXTERNAL :: DNRM2
! 
#endif
  !
  REAL(DP) :: p,            & ! virial pressure
                   vcell,        & ! cell volume
                   avec(3,3),    & ! at(3,3) * alat
                   aveci(3,3),   & ! avec at t-dt
                   avecd(3,3),   & ! d(avec)/dt
                   avec2d(3,3),  & ! d2(avec)/dt2
                   avec2di(3,3), & ! d2(avec)/dt2 at t-dt
                   avec0(3,3),   & ! avec at t = 0
                   sig0(3,3),    & ! sigma at t=0
                   v0              ! volume at t=0
  REAL(DP), ALLOCATABLE ::      &
                   amass_(:),    & ! scaled atomic masses
                   rat(:,:),     & ! atomic positions (lattice coord)
                   rati(:,:),    & ! rat at previous step
                   ratd(:,:),    & ! rat derivatives at current step
                   rat2d(:,:),   & ! rat 2nd derivatives at current step
                   rat2di(:,:),  & ! rat 2nd derivatives at previous step
                   tauold(:,:,:)   ! additional  history variables
  REAL(DP) :: &
           avmod(3), theta(3,3), & ! used to monitor cell dynamics
           enew, e_start,        & ! DFT energy at current and first step
           eold,                 & ! DFT energy at previous step
           uta, eka, eta, ekla, utl, etl, ut, ekint, edyn,  & ! other energies
           acu, ack, acp, acpv, avu, avk, avp, avpv,        & ! acc.& avrg. ener
           tnew=0.0_dp, pv,      & ! instantaneous temperature and p*vcell
           sigmamet(3,3),        & ! sigma = avec^-1 * vcell = bg/alat*omega
           vx2(ntypx), vy2(ntypx), vz2(ntypx),     & ! work vectors
           vmean(ntypx), rms(ntypx), ekin(ntypx),  & ! work vectors
           tempo, time_au
  CHARACTER(LEN=3) :: ios          ! status (old or new) for I/O files
  CHARACTER(LEN=6) :: ipos         ! status ('append' or 'asis') for I/O files
  CHARACTER(LEN=80):: calc_long    ! Verbose description of type of calculation
  LOGICAL :: exst
  INTEGER :: na, nst, ipol, i, j, k
    ! counters
  !
  ! ... I/O units
  !
  INTEGER, PARAMETER :: iun_e    = 21, &
                        iun_eal  = 22, &
                        iun_ave  = 23, &
                        iun_p    = 24, &
                        iun_avec = 25, &
                        iun_tv   = 26
  !
  !
  ! ... Allocate work arrays
  !
  ALLOCATE( amass_(ntyp) )
  amass_(1:ntyp) = amass(1:ntyp) * amu_ry
  ALLOCATE( rat(3,nat) )
  ALLOCATE( rati(3,nat) )
  ALLOCATE( ratd(3,nat) )
  ALLOCATE( rat2d(3,nat) )
  ALLOCATE( rat2di(3,nat) )
  ALLOCATE( tauold(3,nat,3) )
  !
  ! ... open MD history file (if not present this is a new run!)
  !
  CALL seqopn( 4, 'md', 'FORMATTED', exst )
  !
  IF ( .NOT. exst ) THEN
     !
     CLOSE( UNIT = 4, STATUS = 'DELETE' )
     !
     IF ( istep /= 0 ) &
        CALL errore( 'vcsmd', 'previous MD history got lost', 1 )
     !
     tnew  = 0.D0
     acu   = 0.D0
     ack   = 0.D0
     acp   = 0.D0
     acpv  = 0.D0
     avu   = 0.D0
     avk   = 0.D0
     avp   = 0.D0
     avpv  = 0.D0
     nzero = 0
     tauold(:,:,:) = 0.D0
     !
     ! ... set value for eold at first iteration
     !
     eold = etot + 2.D0 * epse 
     !
  ELSE
     !
     ! ... read MD run history
     !
     READ( 4, * ) rati, ratd, rat2d, rat2di, tauold
     READ( 4, * ) aveci, avecd, avec2d, avec2di
     READ( 4, * ) avec0, sig0, v0, e_start, eold
     READ( 4, * ) acu, ack, acp, acpv, avu, avk, avp, avpv, sigmamet
     READ( 4, * ) istep, nzero, ntimes
     !
     CLOSE( UNIT = 4, STATUS = 'KEEP' )
     !
     tauold(:,:,3) = tauold(:,:,2)
     tauold(:,:,2) = tauold(:,:,1)     
     !
  END IF
  !
  istep = istep + 1
  !
  IF ( calc == 'cm' ) THEN
     calc_long="Parrinello-Rahman Damped Cell Dynamics Minimization: "
  ELSE IF ( calc == 'nm' ) THEN
     calc_long="Wentzcovitch Damped Cell Dynamics Minimization: "
  ELSE IF ( calc == 'mm' ) THEN
     calc_long="Beeman Damped Dynamics Minimization: "
  ELSE IF ( calc == 'cd' ) THEN
     calc_long="Parrinello-Rahman Cell Dynamics: "
  ELSE IF ( calc == 'nd' ) THEN
     calc_long="Wentzcovitch Cell Dynamics: "
  ELSE IF ( calc == 'md' ) THEN
     calc_long="Beeman Dynamics: "
  END IF
  !
  conv_ions = .FALSE.
  IF ( calc(2:2) == 'm' ) THEN
     !
     ! ... check if convergence for structural minimization is achieved
     !
     conv_ions = ( (eold - etot) < epse ) .AND. ALL(ABS(force(:,1:nat)) < epsf)
     !
     IF ( lmovecell ) THEN
        DO i = 1, 3
           conv_ions = conv_ions .AND. &
              ( ABS( sigma(i,i) - press ) * ry_kbar * iforceh(i,i) < epsp )
           DO j = ( i + 1 ), 3
              conv_ions = conv_ions .AND. &
              ( ABS( sigma(i,j) ) * ry_kbar * iforceh(i,j) < epsp )
           END DO
        END DO
     END IF
     !
     IF ( conv_ions ) THEN
        !
        WRITE( UNIT = stdout, FMT = '(/,5X,A,/,5X,"convergence achieved, ",&
                &  "Efinal=", F15.8)' ) TRIM(calc_long), etot
        !
        IF ( lmovecell ) THEN
           WRITE( UNIT = stdout, &
               FMT = '(/72("-")//5X,"Final estimate of lattice vectors ", &
                               & "(input alat units)")' )
           WRITE( UNIT = stdout, &
               FMT =  '(3F14.9)') ( ( at(i,k) , i = 1, 3 ) , k = 1, 3 )
           WRITE( UNIT = stdout, &
               FMT =  '("  final unit-cell volume =",F12.4," (a.u.)^3")') omega
           WRITE( UNIT = stdout, &
               FMT =  '("  input alat = ",F12.4," (a.u.)")') alat
        END IF
        !
        CALL output_tau( lmovecell, .TRUE. )
        !
        RETURN
        !
     END IF
     !
  END IF
  ! 
  tauold(:,:,1) = tau(:,:)
  !
  time_au = 0.0000242d0 * e2
  !
  tempo = ( istep - 1 ) * dt * time_au
  !
  IF ( istep == 1 .AND. ( calc(2:2) == 'm' ) ) THEN
        WRITE( stdout,'(/5X,A,/,5x,"convergence thresholds EPSE = ",ES8.2, &
             &  "  EPSF = ",ES8.2)' ) TRIM(calc_long), epse, epsf
  END IF
  !
  WRITE( stdout, '(/5X,"Entering Dynamics;  it = ",I5,"   time = ", &
       &                          F8.5," pico-seconds"/)' ) istep, tempo
  !
  IF ( lconstrain ) THEN
     !
     ! ... we first remove the component of the force along the
     ! ... constraint gradient ( this constitutes the initial
     ! ... guess for the calculation of the lagrange multipliers )
     !
     CALL remove_constr_force( nat, tau, if_pos, ityp, alat, force )
     !
  END IF
  !
  ! ... save cell shape of previous step
  !
  at_old = at
  !
  omega_old = omega
  !
  ! ... Translate
  !
  ! ... define rat as the atomic positions in lattice coordinates
  !
  rat = tau
  !
  CALL cryst_to_cart( nat, rat, bg, -1 )
  !
  avec = alat * at
  !
  ! ... convert forces to lattice coordinates
  !
  CALL cryst_to_cart( nat, force, bg, -1 )
  !
  force = force / alat
  !
  ! ... scale stress to stress*omega
  !
  sigma = sigma * omega
  !
  vcell = omega
  !
  IF ( istep == 1 ) THEN
     !   
     e_start = etot
     !
     enew = etot - e_start
     !
     CALL vcinit( ntyp, nat, ntyp, nat, rat, ityp, avec, vcell, force, if_pos, &
                sigma, calc, temperature, vx2, vy2, vz2, rms, vmean, ekin,     &
                avmod, theta, amass_,cmass, press, p, dt, aveci, avecd, avec2d,&
                avec2di, sigmamet, sig0, avec0, v0, rati, ratd, rat2d, rat2di, &
                enew, uta, eka, eta, ekla, utl, etl, ut, ekint, edyn, iforceh )
     !
  ELSE
     !
     enew = etot - e_start
     !
     CALL vcmove( ntyp, nat, ntyp, ityp, rat, avec, vcell, force, if_pos,      &
                sigma, calc, avmod, theta, amass_,cmass, press, p, dt, avecd,  &
                avec2d, aveci, avec2di, sigmamet, sig0, avec0, v0, ratd, rat2d,&
                rati, rat2di, enew, uta, eka, eta, ekla, utl, etl, ut, ekint,  &
                edyn, temperature, tolp, ntcheck, ntimes, istep, tnew, nzero,  &
                nat, acu, ack, acp, acpv, avu, avk, avp, avpv, iforceh)
     !
  END IF
  !
  pv = p * omega
  !
  IF ( calc(2:2) == 'd' ) THEN
     !
     ! ... Dynamics: write to output files several control quantities
     !
     ! ... NB: at the first iteration files should not be present,
     ! ...     for subsequent iterations they should.
     !
     IF ( istep == 1 ) THEN
        !
        CALL delete_if_present( 'e' )
        CALL delete_if_present( 'eal' )
        CALL delete_if_present( 'ave' )
        CALL delete_if_present( 'p' )
        CALL delete_if_present( 'avec' )
        CALL delete_if_present( 'tv' )
        !
        ios  = 'NEW'
        ipos = 'ASIS'
        !
     ELSE
        !
        ios  = 'OLD'
        ipos = 'APPEND'
        !
     END IF
     !
     OPEN( UNIT = iun_e,    FILE = 'e',    STATUS = ios, &
           FORM = 'FORMATTED', POSITION = ipos )
     OPEN( UNIT = iun_eal,  FILE = 'eal',  STATUS = ios, &
           FORM = 'FORMATTED', POSITION = ipos )
     OPEN( UNIT = iun_ave,  FILE = 'ave',  STATUS = ios, &
           FORM = 'FORMATTED', POSITION = ipos )
     OPEN( UNIT = iun_p,    FILE = 'p',    STATUS = ios, &
           FORM = 'FORMATTED', POSITION = ipos )
     OPEN( UNIT = iun_avec, FILE = 'avec', STATUS = ios, &
           FORM = 'FORMATTED', POSITION = ipos )
     OPEN( UNIT = iun_tv,   FILE = 'tv',   STATUS = ios, &
           FORM = 'FORMATTED', POSITION = ipos )
     !
     nst = istep - 1
     !
     WRITE( iun_e,   101 ) ut, ekint, edyn, pv, nst
     WRITE( iun_eal, 103 ) uta, eka, eta, utl, ekla, etl, nst
     WRITE( iun_ave, 104 ) avu, avk, nst
     WRITE( iun_p,   105 ) press, p, avp, nst
     !
     IF ( calc(1:1) /= 'm' ) &
        WRITE( iun_avec, 103 ) &
            avmod(:), theta(1,2), theta(2,3), theta(3,1), nst
     !
     WRITE( iun_tv, 104 ) vcell, tnew, nst
     !
     CLOSE( UNIT = iun_e,    STATUS = 'KEEP' )
     CLOSE( UNIT = iun_eal,  STATUS = 'KEEP' )
     CLOSE( UNIT = iun_ave,  STATUS = 'KEEP' )
     CLOSE( UNIT = iun_p,    STATUS = 'KEEP' )
     CLOSE( UNIT = iun_avec, STATUS = 'KEEP' )
     CLOSE( UNIT = iun_tv,   STATUS = 'KEEP' )
     !
  END IF
  !
  ! ... update configuration in PWSCF variables
  !
  if (fix_volume) call impose_deviatoric_strain(alat*at, avec)
  !
  if (fix_area) call impose_deviatoric_strain_2d(alat*at, avec)
  !
  at = avec / alat
  !
  CALL volume( alat, at(1,1), at(1,2), at(1,3), omega )
  !
  CALL recips( at(1,1), at(1,2), at(1,3), bg(1,1), bg(1,2) , bg(1,3) )
  !
  tau = rat
  !
  IF ( lmovecell ) THEN
     !
     WRITE( stdout, * ) ' new lattice vectors (alat unit) :'
     WRITE( stdout, '(3F14.9)') ( ( at(i,k) , i = 1, 3 ) , k = 1, 3 )
     WRITE( stdout,'(A,F12.4,A)') '  new unit-cell volume =', omega, ' (a.u.)^3'
     !
  END IF
  !
  WRITE( stdout, * ) ' new positions in cryst coord'
  WRITE( stdout,'(A3,3X,3F14.9)') ( atm(ityp(na)), tau(:,na), na = 1, nat )
  WRITE( stdout, * ) ' new positions in cart coord (alat unit)'
  !
  CALL cryst_to_cart( nat, tau, at, 1 )
  !
  WRITE( stdout,'(A3,3X,3F14.9)') ( atm(ityp(na)), tau(:,na), na = 1, nat )
  WRITE( stdout, '(/5X,"Ekin = ",F14.8," Ry    T = ",F6.1," K ", &
       &       " Etot = ",F14.8)') ekint, tnew, edyn + e_start
  !
  CALL cryst_to_cart( nat, force, at, 1 )
  force = force*alat
  !
  CALL output_tau( lmovecell, .FALSE. )
  !
  ! ... for vcsmd with constraints
  !
  IF ( lconstrain ) THEN
     !
     ! ... check if the new positions satisfy the constrain equation
     !
     CALL check_constraint( nat, tau, tauold(:,:,1), &
                            force, if_pos, ityp, alat, dt**2, amu_ry )
     !
#if ! defined (__REDUCE_OUTPUT)
     !
     WRITE( stdout, '(/,5X,"Constrained forces (Ry/au):",/)')
     !
     DO na = 1, nat
        !
        WRITE( stdout, &
               '(5X,"atom ",I3," type ",I2,3X,"force = ",3F14.8)' ) &
            na, ityp(na), force(:,na)
        !
     END DO
     !
     WRITE( stdout, '(/5X,"Total force = ",F12.6)') DNRM2( 3*nat, force, 1 )
     !
#endif
     !
  END IF
  !
  ! ... save MD history to file
  !
  CALL seqopn( 4, 'md', 'FORMATTED', exst )
  !
  WRITE(4,*) rati, ratd, rat2d, rat2di, tauold
  WRITE(4,*) aveci, avecd, avec2d, avec2di
  WRITE(4,*) avec0, sig0, v0, e_start, etot
  WRITE(4,*) acu, ack, acp, acpv, avu, avk, avp, avpv, sigmamet
  WRITE(4,*) istep, nzero, ntimes
  !
  CLOSE( UNIT = 4, STATUS = 'KEEP' )
  !
  DEALLOCATE( amass_, rat, rati, ratd, rat2d, rat2di, tauold )
  !
  RETURN
  !
101 FORMAT(1X,4D12.5,I6)
103 FORMAT(1X,6D12.5,I6)
104 FORMAT(1X,2D12.5,I6)
105 FORMAT(1X,3D12.5,I6)
  !
END SUBROUTINE vcsmd
