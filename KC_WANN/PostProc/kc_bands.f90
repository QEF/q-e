  SUBROUTINE kc_bands()
    !---------------------------------------------------------------------
    !
    ! ...  This routine takes H(k) on the original mesh of k-points
    ! ...  and interpolates it along the path given in input.
    !
    USE control_kc_wann,      ONLY : nks_bands, xk_bands, centers, use_ws_distance
    USE constants,            ONLY : rytoev
    USE kinds,                ONLY : DP
    USE io_global,            ONLY : stdout
    USE control_kc_wann,      ONLY : num_wann, Hamlt_R
    USE interpolation,        ONLY : read_wannier_centers, ft_ham, print_bands_to_file
    !
    !
    IMPLICIT NONE
    !
    COMPLEX(DP) :: ham_int(num_wann,num_wann)                   ! interpolated H(k)
    COMPLEX(DP) :: eigvc(num_wann,num_wann)
    REAL(DP)    :: eigvl(num_wann,nks_bands)
    INTEGER     :: ik
    !
    !
    ALLOCATE( centers(3,num_wann) )
    !
    WRITE( stdout, '(5x, "STARTING BAND STRUCTURE INTERPOLATION")' )
    !
    IF (use_ws_distance) CALL read_wannier_centers()
    !
    DO ik = 1, nks_bands
      !
      WRITE( stdout, '(/,8x, "KC interpolated eigenvalues at k=", 3f12.4,2x,/)' ) xk_bands(:,ik)
      !
      CALL FT_ham( Hamlt_R, num_wann, ham_int, ik, -1 )
      !
      CALL cdiagh( num_wann, ham_int, num_wann, eigvl(:,ik), eigvc )
      !
      WRITE( stdout, '(6x,8F9.4)' ) eigvl(:,ik)*rytoev
      !
    ENDDO
    !
    CALL print_bands_to_file( eigvl )
    !
    WRITE( stdout, '(/,5x, "ENDING BAND STRUCTURE INTERPOLATION",/)' )
    !
    !
  END SUBROUTINE kc_bands

