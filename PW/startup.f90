!
! Copyright (C) 2001-2008 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
SUBROUTINE startup( nd_nmbr, code, version )
  !----------------------------------------------------------------------------
  !
  ! ... This subroutine initializes MPI and various other things
  !
  USE io_global,  ONLY : stdout, meta_ionode
  USE mp_global,  ONLY : mp_startup, nproc, nogrp, nimage, npool, &
                         nproc_pool, me_image, nproc_image, root_image
  USE control_flags, ONLY : use_task_groups, ortho_para
  !
  IMPLICIT NONE
  !
  CHARACTER (LEN=6)  :: nd_nmbr
  CHARACTER (LEN=6)  :: version
  CHARACTER (LEN=9)  :: code
#if defined __OPENMP
  INTEGER, EXTERNAL  :: omp_get_max_threads
#endif
  !
#if defined (__PARA)
  !
  ! ... parallel case setup :  MPI environment is initialized
  !  
  CALL mp_startup ( use_task_groups, ortho_para )
  !
  ! ... set the processor label for files ( remember that 
  ! ... me_image = 0 : ( nproc_image - 1 ) )
  !
  CALL set_nd_nmbr( nd_nmbr, me_image+1, nproc_image )
  !
  ! ... stdout is printed only by the root_image ( set in init_pool() )
  !
#  if defined (DEBUG)
  !
  IF ( me_image /= root_image ) &
     OPEN( UNIT = stdout, FILE = './out_'//nd_nmbr, STATUS = 'UNKNOWN' )
  !   
#  else
  !
  IF ( me_image /= root_image ) &
     OPEN( UNIT = stdout, FILE = '/dev/null', STATUS = 'UNKNOWN' )
  !   
#  endif
  !
  ! ... information printout
  !  
  IF ( meta_ionode ) THEN
     !
     CALL print_message ( code, version )
     !
#if defined __OPENMP
     WRITE( stdout, '(/5X,"Parallel version (MPI & OpenMP), running on ",&
          &I5," processor cores")' ) nproc * omp_get_max_threads()
     !
     WRITE( stdout, '(5X,"Number of MPI processes:           ",I5)' ) nproc
     !
     WRITE( stdout, '(5X,"Threads/MPI process:               ",I4)' ) omp_get_max_threads()
#else
     WRITE( stdout, '(/5X,"Parallel version (MPI), running on ",&
          &I5," processors")' ) nproc 
#endif
     !
     IF ( nimage > 1 ) &
        WRITE( stdout, &
               '(5X,"path-images division:  nimage    = ",I4)' ) nimage
     IF ( npool > 1 ) &
        WRITE( stdout, &
               '(5X,"K-points division:     npool     = ",I4)' ) npool
     IF ( nproc_pool > 1 ) &
        WRITE( stdout, &
               '(5X,"R & G space division:  proc/pool = ",I4)' ) nproc_pool
     IF ( nogrp > 1 ) &
        WRITE( stdout, &
               '(5X,"wavefunctions fft division:  fft/group = ",I4)' ) nogrp
     !
  END IF   
  !
#else
  !
  ! ... serial case setup :  only information printout
  !
  nd_nmbr = '   '
  !
  CALL print_message ( code, version )
  !
#endif
  !
  RETURN
  !
CONTAINS
  !
  SUBROUTINE print_message ( code, version )
    !
    CHARACTER (LEN=9), INTENT(IN)  :: code
    CHARACTER (LEN=6), INTENT(IN)  :: version
    CHARACTER (LEN=9)  :: cdate, ctime
    !
    CALL date_and_tim( cdate, ctime )
    !
    WRITE( stdout, '(/5X,"Program ",A9," v.",A6," starts on ",A9," at ",A9)' )&
         code, version, cdate, ctime
    !
    WRITE( stdout, '(/5X,"This program is part of the open-source Quantum ",&
         &    "ESPRESSO suite", &
         &/5X,"for quantum simulation of materials; please acknowledge",   &
         &/9X,"""P. Giannozzi et al., J. Phys.:Condens. Matter 21 ",&
         &    "395502 (2009);", &
         &/9X," URL http://www.quantum-espresso.org"", ", &
         &/5X,"in publications or presentations arising from this work. More details at",&
         &/5x,"http://www.quantum-espresso.org/wiki/index.php/Citing_Quantum-ESPRESSO")' )
  END SUBROUTINE PRINT_MESSAGE
  !
   SUBROUTINE set_nd_nmbr( nd_nmbr, node_number, nproc_image )
     !
     IMPLICIT NONE
     !
     CHARACTER(LEN=6), INTENT(OUT) :: nd_nmbr
     INTEGER, INTENT(IN) :: node_number
     INTEGER, INTENT(IN) :: nproc_image
     !
     INTEGER :: nmax, nleft, nfact, n
     !
     nd_nmbr = '      '
     nmax = INT ( LOG10 ( nproc_image + 1.0D-8 ) )
     !
     ! nmax+1=number of digits of nproc_image (number of processors)
     ! 1.0D-8 protects from rounding error if nproc_image is a power of 10
     !
     IF ( nmax+1 > LEN (nd_nmbr) ) &
        CALL errore ( ' startup ', 'insufficient size for nd_nmbr', nmax)
     IF ( nmax < 0) &
        CALL errore ( ' startup ', 'incorrect value for nproc_image', nmax)
     !
     nleft = node_number
     !
     DO n = nmax, 0, -1
        !
        ! decompose node_number (index of this process) into powers of 10:
        !    node_number = i*10^nmax+j*10^(nmax-1)+k*10^(nmax-2)...
        ! i,j,k,... can be equal to 0
        !
        nfact = INT ( nleft/10**n )
        IF ( nfact > 9 ) CALL errore ( ' startup ', 'internal error', 1 )
        nleft = nleft - nfact*10**n
        !
        WRITE( nd_nmbr(nmax-n+1:nmax-n+1), '(I1)' ) nfact
        !
     END DO
     !
     IF ( nleft > 0 ) CALL errore ( ' startup ', 'internal error', 2 )
     !
     RETURN
     !
  END SUBROUTINE set_nd_nmbr
  !     
END SUBROUTINE startup
