!
! Copyright (C) 2001-2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------------------------------------
FUNCTION efermig( et, nbnd, nks, nelec, wk, Degauss, Ngauss, is, isk )
  !--------------------------------------------------------------------
  !! Finds the Fermi energy - Gaussian Broadening. 
  !! (see Methfessel and Paxton, PRB 40, 3616 (1989 )
  !
  USE io_global, ONLY : stdout
  USE kinds,     ONLY : DP
  USE constants, ONLY : rytoev
  USE mp,        ONLY : mp_max, mp_min
  USE mp_pools,  ONLY : inter_pool_comm
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: nks
  !! number of k points in this pool
  INTEGER, INTENT(IN) :: nbnd
  !! number of bands
  INTEGER, INTENT(IN) :: Ngauss
  !! type of smearing technique
  INTEGER, INTENT(IN) :: is
  !! spin label (0 or 1,2)
  INTEGER, INTENT(IN) :: isk(nks)
  !! for each k-point: 1=spin up, 2=spin down
  REAL(DP), INTENT(IN) :: wk(nks)
  !! weight of k points
  REAL(DP), INTENT(IN) :: et(nbnd,nks)
  !! eigenvalues of the hamiltonian
  REAL(DP), INTENT(IN) :: Degauss
  !! smearing parameter
  REAL(DP), INTENT(IN) :: nelec
  !! number of electrons
  REAL(DP) :: efermig
  !! the Fermi energy
  !
  ! ... local variables
  !
  REAL(DP), PARAMETER :: eps = 1.0d-10, eps_cold = 1.0d-2
  !! tolerance for the number of electrons, important for bisection 
  !! smaller tolerance for the number of electrons, important for M-P and Cold smearings
  INTEGER, PARAMETER :: maxiter = 300
  !
  REAL(DP) :: Ef, Eup, Elw, sumkup, sumklw, sumkmid
  REAL(DP), EXTERNAL :: sumkg, sumkg1, sumkg2
  !! Function to compute the number of electrons for a given energy
  !! Function to compute the first derivative of the number of electrons
  !! Function to compute the second derivative of the number of electrons
  REAL(DP), EXTERNAL :: wgauss, w0gauss, w1gauss
  !! Function to compute the occupation
  !! Function to compute the distribution function ( smearing delta)
  !! Function to compute the first derivative of the distribution function
  INTEGER :: i, kpoint, Ngauss_
  INTEGER :: info, maxiter_aux
  REAL(DP) :: Ef_initial_guess, nelec_ef
  !
  !  ... find (very safe) bounds for the Fermi energy:
  !  Elw = lowest, Eup = highest energy among all k-points.
  !  Works with distributed k-points, also if nks=0 on some processor
  !
  Elw = 1.0E+8
  Eup =-1.0E+8
  DO kpoint = 1, nks
    Elw = MIN( Elw, et(1,kpoint) )
    Eup = MAX( Eup, et(nbnd,kpoint) )
  ENDDO
  Eup = Eup + 2 * Degauss
  Elw = Elw - 2 * Degauss
  !
  ! ... find min and max across pools
  !
  CALL mp_max( eup, inter_pool_comm )
  CALL mp_min( elw, inter_pool_comm )
  !

  open(unit=9909, status="replace", file="minimum_bisection.dat")
  open(unit=66990, status="replace", file="efermig_out.dat")

  ! Perform a preliminary determination with the Gaussian broadening
  ! to safely locate Ef mid-gap in the insulating case
  Ngauss_ = 0
  maxiter_aux = maxiter

  call bisection_find_efermi(num_electrons_minus_nelec, Elw, Eup, ef, eps, maxiter_aux, info)
  efermig = ef
  WRITE( 66990, * ) NEW_LINE('a'), "     Bisection Fermi energy:", efermig*rytoev, " Num. electrons:", num_electrons(efermig)

  ! Error handling
  select case( info )
    case( 1 )
      IF (is /= 0) WRITE(66990, '(5x,"Spin Component #",i3)') is
      WRITE( 66990, '(5x,"Warning: too many iterations in bisection" &
        &      5x,"Ef = ",f15.6," N. electons = ",f10.6)' ) &
        Ef * rytoev, num_electrons_minus_nelec(Ef) + nelec
    case( 2 )
      call errore( 'efermig', 'internal error, cannot bracket Ef', 1 )
  end select

  ! If this initial guess already corresponds to the correct number of electron for the actual occupation function, the solution we are done.
  Ngauss_ = Ngauss

  ! In case Ngauss = 0, the function returns here too.
  if( abs_num_electrons_minus_nelec(ef) < eps .or. Ngauss == 0) then 
    
    efermig = ef
    nelec_ef = num_electrons(efermig)
    
    WRITE( 66990, * ) NEW_LINE('a'), "     Final Fermi energy from Bisection using Gaussian smearing:", efermig*rytoev,&
                                     " Num. electrons:", nelec_ef
    goto 98765
  end if

  ! If the initial prospection Ef did not provide the correct number of electrons, use Newton's methods to improve.
  ! Use the prospected Ef as initial guess:

  WRITE( 66990, * ) NEW_LINE('a'), "     Initial guess for Newton methods:", efermig*rytoev ,&
                                  " Num. electrons:", num_electrons(efermig)

  ! Save the initial prospection
  Ef_initial_guess = ef

  maxiter_aux = maxiter

  if( Ngauss_ > 0  .or.  Ngauss_ == -1 ) then ! If methfessel-paxton method or Cold smearing method

    call newton_minimization(sq_num_electrons_minus_nelec, dev1_sq_num_electrons, dev2_sq_num_electrons, &
                                                                                ef, eps, maxiter_aux, info)
  end if

  ! Error handling
  select case( info )
    case( 0 )

      if( abs_num_electrons_minus_nelec(ef) < eps ) then
        WRITE( 66990, * ) &
           "    Newton's Success: mininum or root reached in ", maxiter_aux, " steps."
      else 
        WRITE( 66990, '(5x,"Warning: Newtons finished well, but the number of electron did not reach the required precision. " &
            &      5x,"Ef = ",f15.6," Num. electrons = ",f10.6,"  Num. steps = ",i0)' ) &
            Ef * rytoev, num_electrons(Ef), maxiter_aux
      end if

    case( 1 )

      IF (is /= 0) WRITE(66990, '(5x,"Spin Component #",i3)') is
      WRITE( 66990, '(5x,"Warning: too many iterations in Newtons"/ &
         &      5x,"Ef = ",f15.6," Num. electrons = ",f10.6)' ) &
         Ef * rytoev, num_electrons(Ef)

    case( 2 )
      ! In case the second derivatives go to zero, one should use bisection
      WRITE( 66990, '(5x,"Warning: second derivative went zero."/ &
         &      5x,"Ef = ",f15.6," Num. electrons = ",f10.6,"  Num. steps = ",i0)' ) &
         Ef * rytoev, num_electrons(Ef), maxiter_aux

  end select

  if( (Ngauss_ == -1 .and. abs_num_electrons_minus_nelec(ef) < eps_cold ) .or. &
      (Ngauss_ >=  1 .and. abs_num_electrons_minus_nelec(ef) < eps_cold )       ) then
    efermig = ef
    nelec_ef = num_electrons(efermig)
    WRITE( 66990, * ) NEW_LINE('a'), "     Final Fermi energy from Newton's method:", efermig*rytoev,&
                                     " Num. electrons:", nelec_ef
  else
    Ngauss_ = Ngauss
    maxiter_aux = maxiter

    call bisection_find_efermi(num_electrons_minus_nelec, Elw, Eup, ef, eps, maxiter_aux, info)

    efermig = ef
    nelec_ef = num_electrons(efermig)
    WRITE( 66990, * ) NEW_LINE('a'), "     Final Fermi energy from Bisection using Ngauss smearing:", efermig*rytoev,&
                                     " Num. electrons:", nelec_ef
    WRITE( 66990, '(5x, a)' ) "Warning: Your 'Degauss' is probably too high!"
  end if

  98765 continue

  close(66990)
  close(9909)

  open(unit=9901, status="replace", file="num_electrons_gauss.dat")
  open(unit=9902, status="replace", file="num_electrons.dat")
  open(unit=9903, status="replace", file="dev1_num_electrons.dat")
  open(unit=9904, status="replace", file="dev2_num_electrons.dat")
  Eup = efermig+2.d0
  Elw = efermig-2.d0
  sumkup = (Eup - Elw)/1000.d0
  do i = 0, 1000
    Ef = Elw + sumkup*i
    
    Ngauss_ = 0
    write(unit=9901, fmt="(2f30.16)") Ef, num_electrons_minus_nelec(Ef)

    Ngauss_ = Ngauss
    if( Ngauss_ > 0 ) then
      write(unit=9902, fmt="(2f30.16)") Ef, num_electrons_minus_nelec(Ef)
      write(unit=9903, fmt="(2f30.16)") Ef, dev1_sq_num_electrons(Ef)
      write(unit=9904, fmt="(2f30.16)") Ef, dev2_sq_num_electrons(Ef)
    else 
      write(unit=9902, fmt="(2f30.16)") Ef, num_electrons_minus_nelec(Ef)
      write(unit=9903, fmt="(2f30.16)") Ef, dev1_sq_num_electrons(Ef)
      write(unit=9904, fmt="(2f30.16)") Ef, dev2_sq_num_electrons(Ef)
      ! write(unit=9903, fmt="(2f30.16)") Ef, dev1_num_electrons(Ef)
      ! write(unit=9904, fmt="(2f30.16)") Ef, dev2_num_electrons(Ef)
    end if   
  end do
  close(9901)
  close(9902)
  close(9903)
  close(9904)

  open(unit=9907, status="replace", file="occupation_fuction_M-P.dat")
  write(unit=9907, fmt="(a)") "# Ef      f(x)        f'(x)           f''(x) "

  do i = 0, 1000
    Ef = -10.d0 + (20.d0/1000.d0)*i
    write(unit=9907, fmt="(4f30.16)") Ef, wgauss(Ef,1), w0gauss(Ef,1), w1gauss(Ef,1)
  end do
  close(9907)

  return

  contains

  function num_electrons_minus_nelec(x)
    real(DP), intent(in) :: x
    real(DP) :: num_electrons_minus_nelec

    num_electrons_minus_nelec = num_electrons(x) - nelec
  end function num_electrons_minus_nelec

  function num_electrons(ef)
    real(DP), intent(in) :: ef
    real(DP) :: num_electrons

    num_electrons = sumkg( et, nbnd, nks, wk, Degauss, Ngauss_, ef, is, isk )
  end function num_electrons

  function abs_num_electrons_minus_nelec(ef)
    real(DP), intent(in) :: ef
    real(DP) :: abs_num_electrons_minus_nelec

    abs_num_electrons_minus_nelec = abs(num_electrons_minus_nelec(ef))
  end function abs_num_electrons_minus_nelec

  function sq_num_electrons_minus_nelec(ef)
    real(DP), intent(in) :: ef
    real(DP) :: sq_num_electrons_minus_nelec

    sq_num_electrons_minus_nelec = (num_electrons_minus_nelec(ef))**2
  end function sq_num_electrons_minus_nelec

  function dev1_num_electrons(ef)
    real(DP), intent(in) :: ef
    real(DP) :: dev1_num_electrons

    dev1_num_electrons = sumkg1( et, nbnd, nks, wk, Degauss, Ngauss_, ef, is, isk )
  end function dev1_num_electrons

  function dev2_num_electrons(ef)
    real(DP), intent(in) :: ef
    real(DP) :: dev2_num_electrons

    dev2_num_electrons = sumkg2( et, nbnd, nks, wk, Degauss, Ngauss_, ef, is, isk )
  end function dev2_num_electrons

  function dev1_sq_num_electrons(ef)
    real(DP), intent(in) :: ef
    real(DP) :: dev1_sq_num_electrons

    dev1_sq_num_electrons = 2.d0 * num_electrons_minus_nelec(ef) * dev1_num_electrons(ef)
  end function dev1_sq_num_electrons

  function dev2_sq_num_electrons(ef)
    real(DP), intent(in) :: ef
    real(DP) :: dev2_sq_num_electrons

    dev2_sq_num_electrons = 2.d0 * ( (dev1_num_electrons(ef))**2 + num_electrons_minus_nelec(ef) * dev2_num_electrons(ef) )
  end function dev2_sq_num_electrons

  subroutine newton_minimization(f, f1, f2, x, tol, Nmax, info)
    real(DP),          intent(inout) :: x
    !! Initial guess in the entry. Solution in the exit
    real(DP),          intent(in)    :: tol
    integer,           intent(inout) :: Nmax
    integer,           intent(out)   :: info
    !! 0 = solution found; 1 = max number of step (Nmax) reached; 2 = second derivative is zero

    real(DP)                         :: abstol, x0, denominator, numerator, factor
    integer                          :: i
    real(DP)                         :: f, f1, f2

    abstol = abs(tol)

    open(unit=9905, status="replace", file="minimum_newton_minimization.dat")
    write(unit=9905, fmt=*) "#      ef (Ry)                N(ef)-N0  "
    write(unit=9905, fmt=*) x, f(x)

    write(66990, *) NEW_LINE('a'), "    -Newton's minimization method"

    x0 = x

    factor = 1.0d0

    do i = 1, Nmax
       
       numerator   = f1(x)
       denominator = abs(f2(x))

       ! Checking if the denominator is zero
       if( denominator > abstol ) then
          x = x0 - factor*numerator/denominator
          write(unit=9905, fmt=*) x, f(x)

          ! Checking if a stationary point was achieved
          if( abs(x0-x) < abstol .or.abs_num_electrons_minus_nelec(x) < abstol) then
             info = 0
             Nmax = i
             close(9905)
             return
          ! If a stationary point was not achieved, continue
          else
             x0 = x
          end if

       ! If denominator is zero, return an error
       else 
          info = 2
          close(9905)
          return
       end if
    end do

    ! Checking if max number of steps was reached
    if( i > Nmax ) then
      info = 1
      close(9905)
      return
    end if 
  end subroutine newton_minimization

  subroutine bisection_find_efermi(f, energy_lower_bound, energy_upper_bound, x, tol, Nmax, info)
    real(DP),      intent(in)    :: energy_lower_bound
    real(DP),      intent(in)    :: energy_upper_bound
    real(DP),      intent(out)   :: x
    !! Found Fermi energy at exit
    real(DP),      intent(in)    :: tol
    integer,       intent(inout) :: Nmax
    !! In entry: Max number of steps. In exit: number of step taken.
    integer,       intent(out)   :: info
    !! 0 = solution found; 1 = max number of step (Nmax) reached; 2 = cannot bracket root

    real(DP)                     :: abs_tol, fx, Elw_local, Eup_local
    integer                      :: i
    real(DP)                     :: f

    abs_tol = abs(tol)

    write(9909, *) "#      ef (Ry)                N(ef)-N0 6"

    write(66990,*) NEW_LINE('a'), "    -Bisection root finding method"   

    Elw_local = energy_lower_bound
    Eup_local = energy_upper_bound

    if( f(Elw_local) > abs_tol .or. f(Eup_local) < -abs_tol ) then
      info = 2
      return
    end if

    do i = 1, Nmax

      x = ( Eup_local + Elw_local ) * 0.5d0
      fx = f(x)
      write(unit=9909, fmt=*) x, fx

      ! Was the root found?
      if( abs(fx) < abs_tol ) then
        info = 0
        Nmax = i
        return
      else
        ! Choosing new boundaries
        if( fx < -abs_tol ) then
          Elw_local = x
        else
          Eup_local = x
        end if
      end if
    end do

    ! Checking if max number of steps was reached
    if( i > Nmax ) then
      info = 1
      return
    end if 
  end subroutine bisection_find_efermi
END FUNCTION efermig

