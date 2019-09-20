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
  REAL(DP), PARAMETER :: eps = 1.0d-10
  INTEGER, PARAMETER :: maxiter = 300
  !
  REAL(DP) :: Ef, Eup, Elw, sumkup, sumklw, sumkmid
  REAL(DP), EXTERNAL :: sumkg
  INTEGER :: i, kpoint, Ngauss_
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
  ! ... Bisection method
  !
  ! ... perform a preliminary determination with the Gaussian broadening
  ! to safely locate Ef mid-gap in the insulating case
  !
  !   Ngauss_ = 0 ! currently disabled
  Ngauss_ = Ngauss
  !
1 CONTINUE
  !
  sumkup = sumkg( et, nbnd, nks, wk, Degauss, Ngauss_, Eup, is, isk )
  sumklw = sumkg( et, nbnd, nks, wk, Degauss, Ngauss_, Elw, is, isk )
  IF ( (sumkup - nelec) < -eps .or. (sumklw - nelec) > eps )  &
       CALL errore( 'efermig', 'internal error, cannot bracket Ef', 1 )
  DO i = 1, maxiter
     Ef = (Eup + Elw) / 2.d0
     sumkmid = sumkg( et, nbnd, nks, wk, Degauss, Ngauss_, Ef, is, isk )
     IF ( ABS( sumkmid-nelec ) < eps) THEN
        efermig = Ef
        ! refine the search with the input Ngauss value if not already done
        IF (Ngauss /= Ngauss_) THEN
           Elw = Ef - Degauss ; Eup = Ef + Degauss ; Ngauss_ = Ngauss
           GOTO 1
        ENDIF
        RETURN
     ELSEIF ( (sumkmid-nelec) < -eps) THEN
        Elw = Ef
     ELSE
        Eup = Ef
     ENDIF
  ENDDO
  IF (is /= 0) WRITE(stdout, '(5x,"Spin Component #",i3)') is
  WRITE( stdout, '(5x,"Warning: too many iterations in bisection"/ &
       &      5x,"Ef = ",f10.6," sumk = ",f10.6," electrons")' ) &
       Ef * rytoev, sumkmid
  !
  efermig = Ef
  !
  !
  RETURN
  !
END FUNCTION efermig
