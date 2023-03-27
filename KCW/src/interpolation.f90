!
! Copyright (C) 2003-2021 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
MODULE interpolation
  !---------------------------------------------------------------------
  !
  !! This module contains all the routines important for the
  !! interpolation of the Hamiltonian matrix
  !
  USE kinds,                ONLY : DP
  USE io_global,            ONLY : stdout
  USE lsda_mod,             ONLY : nspin
  USE klist,                ONLY : nks, nkstot, xk
  USE control_kcw,          ONLY : num_wann, Hamlt, Hamlt_R
  !
  !
  IMPLICIT NONE
  !
  SAVE
  !
  !PRIVATE
  PUBLIC
  !
  PUBLIC :: interpolate_ham
  PUBLIC :: real_ham
  PUBLIC :: dealloc_interpolation
  !
  INTEGER :: ios
  COMPLEX(DP) :: imag = (0.D0,1.D0)
  !
  ! ...  end of module-scope declarations
  !
  !---------------------------------------------------------------------
  !
CONTAINS
  !
  !
  SUBROUTINE interpolate_ham()
    !---------------------------------------------------------------------
    !
    ! ...  This routine takes H(k) on the original mesh of k-points
    ! ...  and interpolates it along the path given in input.
    !
    USE control_kcw,          ONLY : nks_bands, xk_bands, centers, use_ws_distance
    USE constants,            ONLY : rytoev
    !
    !
    IMPLICIT NONE
    !
    COMPLEX(DP) :: ham_int(num_wann,num_wann)                   ! interpolated H(k)
    COMPLEX(DP) :: eigvc(num_wann,num_wann)
    REAL(DP)    :: eigvl(num_wann,nks_bands)
    INTEGER     :: ik
    INTEGER     :: k_to_R     ! FT type: (+1) from k- to R-space, (-1) from R- to k-space
    !
    !
    ALLOCATE( Hamlt_R(nkstot/nspin,num_wann,num_wann) )
    ALLOCATE( centers(3,num_wann) )
    !
    CALL real_ham( Hamlt_R )
    !
    WRITE( stdout, '(/,5x,36("="))')
    WRITE( stdout, '(5x, "STARTING BAND STRUCTURE INTERPOLATION")' )
    WRITE( stdout, '(5x,36("="))')
    !
    IF (use_ws_distance) CALL read_wannier_centers()
    !
    k_to_R = -1
    DO ik = 1, nks_bands
      !
      WRITE( stdout, '(/,8x, "KC interpolated eigenvalues at k=", 3f12.4,2x,/)' ) xk_bands(:,ik)
      !
      CALL FT_ham( Hamlt_R, num_wann, ham_int, ik, k_to_R )
      !
      CALL cdiagh( num_wann, ham_int, num_wann, eigvl(:,ik), eigvc )
      !
      WRITE( stdout, '(6x,8F11.4)' ) eigvl(:,ik)*rytoev
      !
    ENDDO
    !
    CALL print_bands_to_file( eigvl )
    !
    WRITE( stdout, '(/,5x, "ENDING BAND STRUCTURE INTERPOLATION",/)' )
    !
    !
  END SUBROUTINE interpolate_ham
  !
  !
  SUBROUTINE real_ham( ham )
    !---------------------------------------------------------------------
    !
    ! ...  This routine defines H(R) starting from H(k) defined on 
    ! ...  the original mesh from the PWscf calculation
    !
    !
    IMPLICIT NONE
    !
    COMPLEX(DP), INTENT(OUT) :: ham(:,:,:)     ! H(R)
    !
    COMPLEX(DP) :: ham_aux(num_wann,num_wann)
    INTEGER :: ir
    INTEGER :: k_to_R
    !
    !
    k_to_R = +1
    DO ir = 1, nkstot/nspin
      !
      CALL FT_ham( Hamlt(1:nkstot/nspin,:,:), num_wann, ham_aux, ir, k_to_R )
      ham(ir,:,:) = ham_aux 
      !
    ENDDO
    !
    !
  END SUBROUTINE real_ham
  !
  !
  SUBROUTINE FT_ham(ham, h_dim, ham_t, ir, k_to_R)
    !---------------------------------------------------------------------
    !
    ! ...  This routine calculates the Fourier Transform of the input
    ! ...  matrix 'ham' from k-space to R-space ( k_to_R=+1 ) or
    ! ...  from R-space to k-space ( k_to_R=-1 ).
    ! ...  The output matrix H(R(ir)) or H(k(ir)) is stored in 'ham_t'.
    !
    USE mp,                   ONLY : mp_sum
    USE mp_global,            ONLY : inter_pool_comm
    USE constants,            ONLY : tpi
    USE lsda_mod,             ONLY : lsda, isk, nspin
    USE control_kcw,          ONLY : spin_component, irvect, xk_bands, centers, &
                                     use_ws_distance
    USE cell_base,            ONLY : at
    !
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN)      :: h_dim
    INTEGER, INTENT(IN)      :: ir         ! k-point (or R-point) index
    INTEGER, INTENT(IN)      :: k_to_R     ! FT type: (+1) from k- to R-space, (-1) from R- to k-space
    COMPLEX(DP), INTENT(IN)  :: ham(nkstot/nspin,h_dim,h_dim)
    !
    COMPLEX(DP), INTENT(OUT) :: ham_t(h_dim,h_dim)
    !
    INTEGER :: iband, jband, ik, ik_eff
    COMPLEX(DP) :: corr_phase
    REAL(DP) :: xq(3)
    REAL(DP) :: dist_ij(3)       ! distance between WFs |0i> and |0j>
    !
    !
    ham_t(:,:) = (0.D0,0D0)
    !
    IF ( k_to_R == 1 ) THEN            ! build H(R) from H(k)
      !
      DO iband = 1, h_dim
        DO jband = 1, h_dim
          DO ik = 1, nks
            !
            IF ( lsda .AND. isk(ik) /= spin_component) CYCLE
            ik_eff = ik - (spin_component -1)*nkstot/nspin
            !
            xq = xk(:,ik)
            CALL cryst_to_cart( 1, xq, at, -1 )
            !
            ham_t(iband,jband) = ham_t(iband,jband) + &
               EXP( - imag * tpi * DOT_PRODUCT( xq, irvect(:,ir) ) ) * ham(ik_eff,iband,jband)
            !
          ENDDO
          !
          CALL mp_sum(ham_t, inter_pool_comm)
          !
        ENDDO
      ENDDO
      !
      ham_t = ham_t / (nkstot/nspin)   ! 1/Nk factor for the FT from H(k) to H(R)
      !
    ELSE IF ( k_to_R == -1 ) THEN      ! build H(k) from H(R). The k-vector here is taken from 
      !                                ! the list xk_bands coming form the card K_POINTS
      DO iband = 1, h_dim
        DO jband = 1, h_dim
          !
          IF ( use_ws_distance ) THEN
            !
            dist_ij = centers(:,jband) - centers(:,iband)
            !
          ELSE
            !
            dist_ij = (/0.,0.,0./)
            !
          ENDIF
          !
          DO ik = 1, nkstot/nspin
            !
            xq = xk_bands(:,ir)
            !
            CALL correct_phase( dist_ij, irvect(:,ik), xq, corr_phase )
            !
            ham_t(iband,jband) = ham_t(iband,jband) + &
               EXP( + imag * tpi * DOT_PRODUCT( xq, irvect(:,ik) ) ) * corr_phase * ham(ik,iband,jband) 
            !
          ENDDO
        ENDDO
      ENDDO
      !
    ELSE
      !
      CALL errore( 'FT_ham', 'argument k_to_R must be either +1 or -1', ABS(ios) )
      !
    ENDIF
    !
    !
  END SUBROUTINE FT_ham
  !
  !
  SUBROUTINE correct_phase( wf_dist, Rvec, qvec, phase )
    !---------------------------------------------------------------------
    !
    ! ...  This routine accounts for the PBC and calculates the right R-vector
    ! ...  entering in the Fourier transfrom H(R)-->H(k). In the event of multiple
    ! ...  equidistant R-vectors, it considers an average of all of them.
    ! ...  Finally it gives in output the correct phase factor to enter in the FT.
    !
    USE cell_base,            ONLY : at
    USE control_kcw,          ONLY : mp1, mp2, mp3
    USE constants,            ONLY : tpi
    !
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN) :: Rvec(3)
    REAL(DP), INTENT(IN) :: qvec(3)
    REAL(DP), INTENT(IN) :: wf_dist(3)
    !
    COMPLEX(DP), INTENT(OUT) :: phase
    !
    INTEGER :: i, j, k
    INTEGER :: counter           ! counts the equidistant R-vectors
    INTEGER :: Tvec(3)       ! primitive lattice vector of the supercell
    REAL(DP) :: eff_dist(3), eff_dist_aux(3)
    REAL(DP) :: dist_min, dist
    !
    !
    eff_dist = wf_dist + Rvec
    CALL cryst_to_cart( 1, eff_dist, at, 1 )
    dist_min = SQRT( eff_dist(1)**2 + eff_dist(2)**2 + eff_dist(3)**2 )
    !
    phase = 0.
    counter = 0
    !
    DO i = -1, 1
      DO j = -1, 1
        DO k = -1, 1
          !
          Tvec(:) = (/i*mp1, j*mp2, k*mp3/)     ! defines a kind of WS cell for the supercell
          !
          eff_dist_aux = wf_dist + Rvec + Tvec
          !
          CALL cryst_to_cart( 1, eff_dist_aux, at, 1 )
          dist = SQRT( eff_dist_aux(1)**2 + &
                       eff_dist_aux(2)**2 + &
                       eff_dist_aux(3)**2 )
          !
          IF ( ABS( dist - dist_min ) .lt. 1.e-3 ) THEN
            !
            ! an equidistant replica is found
            !
            phase = phase + EXP( + imag * tpi * DOT_PRODUCT( qvec, Tvec ) )
            counter = counter + 1
            !
          ELSE IF ( dist .lt. dist_min ) THEN
            !
            ! a new nearest replica is found
            ! reset dist_min, phase and counter
            !
            dist_min = dist
            phase = EXP( + imag * tpi * DOT_PRODUCT( qvec, Tvec ) )
            counter = 1
            !
          ELSE
            !
            ! this replica is rejected
            !
            CYCLE
            !
            !
          ENDIF
          !
        ENDDO
      ENDDO
    ENDDO
    !
    phase = phase / counter 
    !
    !
  END SUBROUTINE correct_phase
  !
  !
  SUBROUTINE read_wannier_centers( )
    !---------------------------------------------------------------------
    !
    ! ...  This routine reads the centers of Wannier functions from .xyz
    ! ...  file print out by Wannier90, fold them into the R=0 primitive
    ! ...  cell and gives them in output (in crystal units)
    !
    USE io_global,            ONLY : ionode
    USE cell_base,            ONLY : alat, bg
    USE constants,            ONLY : BOHR_RADIUS_ANGS
    USE control_kcw,          ONLY : seedname, have_empty, num_wann_occ, num_wann_emp, &
                                     centers, use_ws_distance
    !
    !
    IMPLICIT NONE
    !
    LOGICAL :: exst
    LOGICAL :: check_emp=.FALSE.
    INTEGER :: n, line, nlines
    CHARACTER(LEN=268) :: filename
    CHARACTER(LEN=256) :: input_line
    !
    !
    IF ( ionode ) THEN
      !
      filename = trim(seedname)//'_centres.xyz'
      nlines = num_wann_occ
      !
50    INQUIRE( file=filename, exist=exst )
      !
      IF ( .not. exst .and. .not. check_emp ) THEN 
        CALL infomsg('read_wannier_centers','WARNING: centres.xyz NOT FOUND, disabling WS distance')
        use_ws_distance = .false.
        RETURN
      ENDIF
      !
      IF ( .not. exst .and. check_emp ) THEN 
        CALL infomsg('read_wannier_centers','WARNING: emp_centres.xyz NOT FOUND, disabling WS distance')
        use_ws_distance = .false.
        RETURN
      ENDIF
      !
      !
      OPEN( 100, file=filename, form='formatted', status='old' )
      !
      READ( 100, *, end=10, err=20 )    ! skip 1st line
      READ( 100, *, end=10, err=20 )    ! skip 2nd line
      !
      DO n = 1, nlines
        !
        READ( 100, '(a256)', end=10, err=20 ) input_line
        !
        IF ( input_line(1:1) .ne. 'X' ) CALL errore( 'read_wannier_centers', &
                'X must precede each Wannier center line', 1 )
        !
        line = n
        IF ( check_emp ) line = n + num_wann_occ
        READ( input_line(2:), *, end=10, err=20 ) centers(:,line)
        !
      ENDDO
      !
      READ( 100, * ) input_line
      IF ( input_line(1:1) == 'X' ) CALL errore( 'read_wannier_centers', &
              'Missing some center! Check num_wann', 1 )
      !
      CLOSE( 100 )
      !
      IF ( have_empty .and. .not. check_emp ) THEN
        !
        filename = trim(seedname)//'_emp_centres.xyz'
        nlines = num_wann_emp
        check_emp = .TRUE.
        GO TO 50
        !
      ENDIF
      !      
    ENDIF
    !
    centers = centers / ( alat * BOHR_RADIUS_ANGS )     ! alat always in BOHR ??? (TO BE CHECKED)
    !
    DO n = 1, num_wann_occ+num_wann_emp
      !
      CALL cryst_to_cart( 1, centers(:,n), bg, -1 )
      !
    ENDDO
    !
    RETURN
    !
10  CALL errore ( 'read_wannier_centers', 'end of file while reading', 1 ) 
20  CALL errore ( 'read_wannier_centers', 'error while reading', 1 ) 
    !
    !
  END SUBROUTINE read_wannier_centers 
  !
  !
  SUBROUTINE print_bands_to_file( eigvl, filename)
    !---------------------------------------------------------------------
    !
    ! ...  This routine prints the interpolated band energies in the  
    ! ...  same format of bands.x program
    !
    USE cell_base,            ONLY : bg
    USE control_kcw,          ONLY : xk_bands, nks_bands, num_wann
    USE io_files,             ONLY : prefix
    USE io_global,            ONLY : ionode
    USE constants,            ONLY : rytoev
    !
    !
    IMPLICIT NONE
    !
    INTEGER :: ik, iband
    REAL(DP), INTENT(IN) :: eigvl(:,:)
    REAL(DP) :: path_coor(nks_bands)     ! linearized coordinate of the path
    REAL(DP) :: kvec(3)
    REAL(DP) :: kvec_aux(3)
    REAL(DP) :: k_dist, k_dist_save 
    CHARACTER(268), INTENT(IN), OPTIONAL :: filename
    CHARACTER(268) :: filename_
    !
    !
    IF ( ionode ) THEN
      !
      DO ik = 1, nks_bands
        !
        kvec = xk_bands(:,ik)
        CALL cryst_to_cart( 1, kvec, bg, 1 )
        !
        IF ( ik == 1 ) THEN
          !
          path_coor(ik) = 0.
          !
          kvec_aux = xk_bands(:,ik+1)
          CALL cryst_to_cart( 1, kvec_aux, bg, 1 )
          k_dist_save = SQRT( (kvec_aux(1) - kvec(1))**2 + &
                              (kvec_aux(2) - kvec(2))**2 + &
                              (kvec_aux(3) - kvec(3))**2 )
          !
          CYCLE
          !
        ENDIF
        !
        kvec_aux = xk_bands(:,ik-1)
        CALL cryst_to_cart( 1, kvec_aux, bg, 1 )
        k_dist = SQRT( (kvec_aux(1) - kvec(1))**2 + &
                       (kvec_aux(2) - kvec(2))**2 + &
                       (kvec_aux(3) - kvec(3))**2 )

        !
        IF ( k_dist > 5*k_dist_save ) THEN
          !
          path_coor(ik) = path_coor(ik-1)
          !
        ELSE IF ( k_dist > 1.e-4 ) THEN
          !
          path_coor(ik) = path_coor(ik-1) + k_dist
          k_dist_save = k_dist
          !
        ELSE
          !
          path_coor(ik) = path_coor(ik-1) + k_dist
          !
        ENDIF
        !
      ENDDO
      !
      !
      filename_ = trim(prefix)//'.kcw_bands.dat'
      IF (PRESENT (filename) ) filename_=filename
      OPEN( 100, file=filename_, status='unknown' )
      !
      DO iband = 1, num_wann
        !
        DO ik = 1, nks_bands
          !
          WRITE( 100, '(2f10.4)' ) path_coor(ik), eigvl(iband,ik)*rytoev
          !
        ENDDO
        !
        WRITE( 100, * )
        !
      ENDDO
      !
      CLOSE( 100 )
      !
    ENDIF
    !
    !
  END SUBROUTINE print_bands_to_file
  !
  !
  SUBROUTINE dealloc_interpolation()
    !---------------------------------------------------------------------
    !
    ! ...  This routine dealloctes all the interpolation-related arrays 
    !
    USE control_kcw,          ONLY : xk_bands, wk_bands, centers
    !
    !
    IMPLICIT NONE
    !
    !
    DEALLOCATE( centers )
    DEALLOCATE( xk_bands )
    DEALLOCATE( wk_bands )
    DEALLOCATE( Hamlt_R )
    !
    !
  END SUBROUTINE dealloc_interpolation
  !
  !
END MODULE interpolation
