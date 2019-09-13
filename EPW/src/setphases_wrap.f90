  !                                                                            
  ! Copyright (C) 2010-2016 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino 
  ! Copyright (C) 2007-2009 Jesse Noffsinger, Brad Malone, Feliciano Giustino  
  !                                                                            
  ! This file is distributed under the terms of the GNU General Public         
  ! License. See the file `LICENSE' in the root directory of the               
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .             
  !                                                                            
  !
  !---------------------------------------------------------------------
  SUBROUTINE setphases_wrap()
  !---------------------------------------------------------------------
  !!
  !!  This is the wrapper which is used to set the phases of the wavefunctions  
  !!  at k and k+q on the coarse mesh.  It should only be called once.
  !!  Note that the phases at k+q are for the input 'q' vector, not
  !!  the one in the dynamical coarse list.
  !!  
  !!  The phases for all wavefunctions are now stored in umat_all
  !!  
  !---------------------------------------------------------------------
  USE kinds,           ONLY : DP
  use klist,           only : nkstot
  use io_global,       only : ionode, stdout
  use mp_global,       only : inter_pool_comm, nproc_pool
  use mp,              only : mp_sum
  use elph2,           only : umat, umat_all
  use pwcom,           only : nbnd, nks
  !
  IMPLICIT NONE
  !
  INTEGER :: ik
  !! K-point 
  INTEGER :: ibnd
  !! Band-index
  INTEGER :: jbnd
  !! Band index
  INTEGER :: ierr
  !! Error status
  REAL(KIND = DP) :: zero_vect(3)
  !! Real vector
  !
  IF (nproc_pool>1) CALL errore('setphases_wrap', 'only one proc per pool', 1)
  !
  ALLOCATE(umat_all(nbnd, nbnd, nkstot), STAT = ierr)
  IF (ierr /= 0) CALL errore('setphases_wrap', 'Error allocating umat_all', 1)
  ALLOCATE(umat(nbnd, nbnd, nks), STAT = ierr)
  IF (ierr /= 0) CALL errore('setphases_wrap', 'Error allocating umat', 1)
  umat_all = (0.d0, 0.d0)
  zero_vect = 0.d0
  !
  WRITE(stdout, '(5x,a)') 'No wavefunction gauge setting applied' 
  !
  IF (ionode) THEN
    DO ik = 1, nkstot
      DO ibnd = 1, nbnd
        DO jbnd = 1, nbnd
          IF (ibnd == jbnd) then
            umat_all(ibnd, jbnd, ik) = (1.d0, 0.d0)
          ELSE
            umat_all(ibnd, jbnd, ik) = (0.d0, 0.d0)
          ENDIF
        ENDDO
      ENDDO
    ENDDO
  ENDIF
  DO ik = 1, nks
    DO ibnd = 1, nbnd
      DO jbnd = 1, nbnd
        IF (ibnd == jbnd) then
          umat(ibnd, jbnd, ik) = (1.d0, 0.d0)
        ELSE
          umat(ibnd, jbnd, ik) = (0.d0, 0.d0)
        ENDIF
      ENDDO
    ENDDO
  ENDDO
  !
  ! collect the global phase-setting matrix
  !
  CALL mp_sum(umat_all, inter_pool_comm)
  !
  !---------------------------------------------------------------------
  END SUBROUTINE setphases_wrap
  !---------------------------------------------------------------------
