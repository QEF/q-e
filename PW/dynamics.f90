!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------

subroutine dynamics
  !-----------------------------------------------------------------------
  !
  !   This routine performs one step of molecular dynamics evolution using
  !   the Verlet algorithm. Parameters:
  !   mass         mass of the atoms
  !   dt           time step
  !   temperature  starting temperature
  !                The starting velocities of atoms are set accordingly
  !                to the starting temperature, in random directions.
  !                The initial velocity distribution is therefore a constant
  !   delta_T, nraise are used to change the temperature as follows:
  !   delta_T = 1 :                 nothing is done.
  !   delta_T <> 1 and delta_T >0 : every step the actual temperature is
  !                                 multiplied by delta_T, this is done
  !                                 rescaling all the velocities
  !   delta_T < 0 :                 every 'nraise' step the temperature
  !                                 reduced by -delta_T
  !
  !     DA 1997
  !
#include "machine.h"
  USE io_global, ONLY : stdout
  USE kinds, ONLY: DP
  USE constants, ONLY: amconv
  USE basis, ONLY: nat, ntyp, tau, ityp, atm
  USE brilz, ONLY: alat
  USE dynam, ONLY: amass, temperature, dt, delta_t, nraise
  USE ener, ONLY: etot
  USE force_mod, ONLY: force
  USE klist, ONLY: nelec
  USE relax, ONLY: fixatom, dtau_ref, starting_diag_threshold
  USE control_flags, ONLY: ethr, upscale, tr2, imix, alpha0, beta0, istep
  use io_files, only : prefix
#ifdef __PARA
  use para
#endif
  implicit none
  integer :: natoms                ! number of moving atoms
  real(kind=DP), allocatable ::  &
             tauold (:,:,:)        ! previous positions of atoms
  real(kind=DP), allocatable ::  &
             a (:,:),  mass (:)    ! accelerations and masses of atoms
  real(kind=DP)  ::  ekin          ! ionic kinetic energy
  real(kind=DP)  ::  total_mass, temp_new, tempo, norm_of_dtau
  real(kind=DP)  ::  ml (3)        ! total linear momentum
  integer :: na, ipol, it ! counters
  logical :: exst
  real(kind=DP), parameter ::  convert_E_to_temp= 315642.28d0 * 0.5d0, &
                               eps = 1.d-6

  allocate (mass(   nat))    
  allocate (a(  3, nat))    
  allocate (tauold(   3,  nat, 3))    
  tauold(:,:,:) = 0.d0

  dtau_ref = 0.2
  natoms = nat - fixatom
  !
  !  one Ryd a.u. of time is 4.84*10^-17 seconds, i.e. 0.0484  femtosecond
  !
  call seqopn (4, trim(prefix)//'.md', 'formatted', exst)
  if (.not.exst) then
     close (unit = 4, status = 'delete')
     WRITE( stdout, '(/5x,"Starting temperature = ",f8.2," K")') &
          temperature
     do na = 1, ntyp
        WRITE( stdout, '(5x,"amass(",i1,") = ",f6.2)') na, amass (na)
     enddo

     WRITE( stdout, '(5x,"Time step = ",f6.2," a.u.,   ",f6.4, &
          &       " femto-seconds")') dt, dt * 0.0484
     !
     ! masses in atomic rydberg units
     !
     total_mass = 0.d0
     do na = 1, nat
        mass (na) = amass (ityp (na) ) * amconv
        total_mass = total_mass + mass (na)
     enddo
     !
     !   initial thermalization. N.B. tau is in units of alat
     !
     call start_therm (mass, tauold)
     tempo = 0.d0
     temp_new = temperature
     it = 0
  else
     read (4, * ) temp_new, mass, total_mass, tauold, tempo, it
     close (unit = 4, status = 'keep')
     istep = it + 1
  endif
  tempo = tempo + dt * 0.0000484
  it = it + 1
  if (mod (it, nraise) .eq.0.and.delta_T.lt.0) then
     WRITE( stdout, '(/5x,"Thermalization: delta_T = ",f6.3, &
          &  ", T = ",f6.1)')  - delta_T, temp_new - delta_T
     call thermalize (temp_new, temp_new - delta_T, tauold)
  endif
  if (delta_T.ne.1.d0.and.delta_T.ge.0) then
     WRITE( stdout, '(/5x,"Thermalization: delta_T = ",f6.3, &
          &  ", T = ",f6.1)') delta_T, temp_new * delta_T
     call thermalize (temp_new, temp_new * delta_T, tauold)
  endif
  WRITE( stdout, '(/5x,"Entering Dynamics;  it = ",i5,"   time = ", &
       &                          f8.5," pico-seconds"/)') it, tempo
  !
  ! calculate accelerations in a.u. units / alat
  !
  do na = 1, natoms
     a (1, na) = force (1, na) / mass (na) / alat
     a (2, na) = force (2, na) / mass (na) / alat
     a (3, na) = force (3, na) / mass (na) / alat
  enddo
  !
  !     save the previous two steps ( a total of three steps is saved)
  !
  tauold (:, :, 3) = tauold (:, :, 2)
  tauold (:, :, 2) = tauold (:, :, 1)
  !
  ! move the atoms accordingly to the classical equation of motion
  !
  call verlet (tau, tauold, a, natoms, ekin, mass, ml, dt)
  !
  !     find the best coefficients for the extrapolation of the potential
  !
  call find_alpha_and_beta (nat, tau, tauold, alpha0, beta0)
#ifdef __PARA
  if (me.eq.1) call poolbcast (1, alpha0)
  call broadcast (1, alpha0)
  if (me.eq.1) call poolbcast (1, beta0)
  call broadcast (1, beta0)
#endif
  !
  ! calculate the "norm" of the step (used to update the convergence threshold)
  !
  norm_of_dtau = 0.d0
  do na = 1, nat
     do ipol = 1, 3
        norm_of_dtau = norm_of_dtau + (tau(ipol,na) - tauold(ipol,na,1)) **2
     enddo
  enddo
  norm_of_dtau = sqrt (norm_of_dtau)
  !
  ! ... ethr is now updated in electrons with a different procedure
  ! ... this value of ethr is overwritten apart when the old style
  ! ... update is used (see OLDSTYLE precompiler variable in electrons)
  !       
  if (imix.lt.0) then
     ethr = starting_diag_threshold * &
         max (1.d0 / upscale, min (1.d0, norm_of_dtau / dtau_ref) ) **2
  else
     ethr = tr2 / nelec
  end if
  !
  ! find the new temperature
  !
  temp_new = 2.d0 / 3.d0 * ekin * alat**2 / natoms * convert_E_to_temp
  !
  ! save on file needed quantity
  !
  call seqopn (4, trim(prefix)//'.md', 'formatted', exst)
  write (4, * ) temp_new, mass, total_mass, tauold, tempo, it
  close (unit = 4, status = 'keep')
  !
  call output_tau( .FALSE. )
  WRITE( stdout, '(/5x,"Ekin = ",f14.8," Ryd   T = ",f6.1," K ", &
       &       " Etot = ",f14.8)') ekin*alat**2, temp_new, ekin*alat**2+etot
  !
  !  total linear momentum must be zero if all atoms move
  !
  if (fixatom == 0) then
     if ( abs(ml(1)) > eps .OR. abs(ml(2)) > eps .OR. abs(ml(3)) > eps ) then
        call errore ('dynamics', 'Total linear momentum <> 0', - 1)
        WRITE( stdout, '(5x,"Linear momentum: ",3f12.8)') ml
     end if
  endif


  deallocate (tauold)
  deallocate (a)
  deallocate (mass)

  return

end subroutine dynamics
!---------------------------------------------------------------------
subroutine verlet (rnew, rold, a, n, ec, mass, ml, dt)
  !---------------------------------------------------------------------
  !
  ! Verlet algorithm to update atomic position
  !
  USE kinds, ONLY: DP
  implicit none
  ! INPUT
  integer :: n  ! number of particles
  real(kind=DP) :: dt, & ! time step
             a (3, n), & ! accelerations
             mass (n)    ! atom masses
  ! OUTPUT
  real(kind=DP) :: ec, ml (3)  ! kinetic energy and total linear momentum
  ! INPUT/OUTPUT
  real(kind=DP) :: rold (3, n), & ! in: previous, out: present atomic positions
                   rnew (3, n)    ! in: present,  out: new     atomic positions
  ! LOCAL
  integer :: i
  real(kind=DP) :: dtsq, dt2, r (3), v (3)
  ! dtsq=dt**2, dt2=2*dt
  !
  ml(:) = 0.d0
  ec = 0.d0
  dtsq = dt**2
  !
  dt2 = dt * 2.0
  do i = 1, n
     r(:) = 2.0 * rnew(:,i) - rold(:,i) + a(:,i) * dtsq
     v(:) = (r(:) - rold(:,i) ) / dt2
     rold(:,i) = rnew(:,i)
     rnew(:,i) = r(:)
     ml(:) = ml(:) + v(:) * mass(i)
     ec = ec + 0.5 * mass(i) * (v(1)**2 + v(2)**2 + v(3)**2)
  enddo
  return

end subroutine verlet
!-----------------------------------------------------------------------
subroutine start_therm (mass, tauold)
  !-----------------------------------------------------------------------
  !
  !     Starting thermalization of the system
  !
  USE kinds, ONLY: DP
  USE basis, ONLY: nat, tau
  USE brilz, ONLY: alat
  USE dynam, ONLY: temperature, dt!, delta_t, nraise
  USE relax, ONLY: fixatom
  USE symme, only: invsym, nsym, irt
#ifdef __PARA
  use para
#endif
  implicit none
  real(kind=DP)  :: mass (nat),  tauold (3, nat)
  !
  integer :: na, nb, natoms
  real(kind=DP)  ::  total_mass, temp_new, aux, convert_E_to_temp, velox,&
       ek, ml(3), direzione_x, direzione_y, direzione_z, modulo, rndm
  ! ek = kinetic energy
  real(kind=DP), allocatable :: step(:,:)
  external rndm
  parameter (convert_E_to_temp = 315642.28d0 * 0.5d0)
  !
  allocate (step(3,nat))    
  aux = temperature / convert_E_to_temp
  natoms = nat - fixatom
  ml(:) = 0.d0
  !
  ! velocity in random direction, with modulus accordingly to mass and
  ! temperature: 3/2KT = 1/2mv^2
  !
#ifdef __PARA
  !
  ! only the first processor calculates ...
  !
  if (me.eq.1.and.mypool.eq.1) then
#endif
     do na = 1, natoms
        !
        ! N.B. velox is in a.u. units /alat
        !
        velox = sqrt (3.d0 * aux / mass (na) ) / alat
        direzione_x = rndm () - .5d0
        direzione_y = rndm () - .5d0
        direzione_z = rndm () - .5d0
        modulo = sqrt (direzione_x**2 + direzione_y**2 + direzione_z**2)
        step (1, na) = velox / modulo * direzione_x
        step (2, na) = velox / modulo * direzione_y
        step (3, na) = velox / modulo * direzione_z
     enddo
#ifdef __PARA
     !
     ! ... and distributes the velocities
     !
  endif
  if (me.eq.1) call poolbcast (3 * natoms, step)

  call broadcast (3 * natoms, step)
#endif
  !
  ! if there is inversion symmetry equivalent atoms have opposite velocities
  !
  if (invsym) then
     do na = 1, natoms
        nb = irt (nsym / 2 + 1, na)
        if (nb.gt.na) then
           step (1, nb) = - step (1, na)
           step (2, nb) = - step (2, na)
           step (3, nb) = - step (3, na)
        endif
        !
        ! the atom on the inversion center is kept fixed
        !
        if (na.eq.nb) then
           step (1, na) = 0.d0
           step (2, na) = 0.d0
           step (3, na) = 0.d0
        endif
     enddo
  else
     !
     ! put total linear momentum equal zero if all atoms move
     !
     if (fixatom.eq.0) then
        total_mass = 0.d0
        do na = 1, natoms
           total_mass = total_mass + mass (na)
           ml (1) = ml (1) + mass (na) * step (1, na)
           ml (2) = ml (2) + mass (na) * step (2, na)
           ml (3) = ml (3) + mass (na) * step (3, na)
        enddo
        ml (1) = ml (1) / total_mass
        ml (2) = ml (2) / total_mass
        ml (3) = ml (3) / total_mass
     endif
  endif
  !
  !     -step is the velocity
  !
  ek = 0.d0
  do na = 1, natoms
     tauold (1, na) = (step (1, na) - ml (1) ) * dt + tau (1, na)
     tauold (2, na) = (step (2, na) - ml (2) ) * dt + tau (2, na)
     tauold (3, na) = (step (3, na) - ml (3) ) * dt + tau (3, na)
     ek = ek + 0.5d0 * mass (na) * ( (step (1, na) - ml (1) ) **2 + &
          (step (2, na) - ml (2) ) **2 + (step (3, na) - ml (3) ) **2)
  enddo
  !
  !     after the velocity of the center of mass has been subtracted the
  !     temperature is usually changed. Set again the temperature to the
  !     right value.
  !
  temp_new = 2.d0 * ek / (3.d0 * natoms) * alat**2 * convert_E_to_temp

  call thermalize (temp_new, temperature, tauold)

  deallocate(step)
  return

end subroutine start_therm
!-------------------------------------------------------------------

subroutine thermalize (temp_old, temp_new, tauold)
  !-------------------------------------------------------------------
  USE kinds, ONLY: DP
  USE basis, ONLY: nat, tau
  USE relax, ONLY: fixatom
  USE dynam, ONLY: dt
  implicit none
  real(kind=DP) :: tauold (3, nat), temp_new, temp_old
  !
  integer :: na, natoms
  real(kind=DP) ::  velox, aux
  !
  natoms = nat - fixatom
  !
  !     rescale the velocities by a factor 3/2KT/Ek
  !
  if (temp_new.gt.0.d0.and.temp_old.gt.0.d0) then
     aux = sqrt (temp_new / temp_old)
  else
     aux = 0.d0
  endif
  do na = 1, natoms
     velox = (tau (1, na) - tauold (1, na) ) / dt
     tauold (1, na) = tau (1, na) - dt * aux * velox
     velox = (tau (2, na) - tauold (2, na) ) / dt
     tauold (2, na) = tau (2, na) - dt * aux * velox
     velox = (tau (3, na) - tauold (3, na) ) / dt
     tauold (3, na) = tau (3, na) - dt * aux * velox
  enddo

  !
  return

end subroutine thermalize
!-----------------------------------------------------------------------
subroutine find_alpha_and_beta (nat, tau, tauold, alpha0, beta0)
  !-----------------------------------------------------------------------
  !
  !     This routine finds the best coefficients alpha0 and beta0 so that
  !
  !     | tau(t+dt) - tau' | is minimum, where
  !
  !     tau' = alpha0 * ( tau(t) - tau(t-dt) ) +
  !             beta0 * ( tau(t-dt) - tau(t-2*dt) )
  !
  USE kinds, ONLY: DP
  implicit none

  integer :: nat, na, ipol

  real(kind=DP) :: chi, alpha0, beta0, tau (3, nat), tauold (3, nat, 3)
  real(kind=DP) :: a11, a12, a21, a22, b1, b2, c, det
  !
  ! solution of the linear system
  !
  a11 = 0.d0
  a12 = 0.d0
  a21 = 0.d0
  a22 = 0.d0 + 1.d-12
  b1 = 0.d0
  b2 = 0.d0
  c = 0.d0
  do na = 1, nat
     do ipol = 1, 3
        a11 = a11 + (tauold (ipol, na, 1) - tauold (ipol, na, 2) ) **2
        a12 = a12 + (tauold (ipol, na, 1) - tauold (ipol, na, 2) ) &
             * (tauold (ipol, na, 2) - tauold (ipol, na, 3) )
        a22 = a22 + (tauold (ipol, na, 2) - tauold (ipol, na, 3) ) **2
        b1 = b1 - (tauold (ipol, na, 1) - tau (ipol, na) ) * (tauold ( &
             ipol, na, 1) - tauold (ipol, na, 2) )
        b2 = b2 - (tauold (ipol, na, 1) - tau (ipol, na) ) * (tauold ( &
             ipol, na, 2) - tauold (ipol, na, 3) )
        c = c + (tauold (ipol, na, 1) - tau (ipol, na) ) **2
     enddo
  enddo
  a21 = a12
  !
  det = a11 * a22 - a12 * a21
  if (det.lt.0d0) call errore ('find_alpha_and_beta', ' det.le.0', 1)
  !
  ! det > 0 case:  a well defined minimum exists
  !
  if (det.gt.0d0) then
     alpha0 = (b1 * a22 - b2 * a12) / det
     beta0 = (a11 * b2 - a21 * b1) / det
  else
     !
     ! det = 0 case: the two increments are linearly dependent, chose
     !               solution with beta=0 (discard oldest configuration)
     !
     alpha0 = 1.d0
     if (a11.gt.0.d0) alpha0 = b1 / a11
     beta0 = 0.d0

  endif
  chi = 0.d0
  do na = 1, nat
     do ipol = 1, 3
        chi = chi + ( (1 + alpha0) * tauold (ipol, na, 1) + (beta0 - &
             alpha0) * tauold (ipol, na, 2) - beta0 * tauold (ipol, na, 3) &
             - tau (ipol, na) ) **2
     enddo
  enddo

!  WRITE( stdout, * ) chi, alpha0, beta0
  return

end subroutine find_alpha_and_beta

