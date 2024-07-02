  !-----------------------------------------------------------------------
  !copy from exx_base
  SUBROUTINE kcw_set_symm( nr1, nr2, nr3, nr1x, nr2x, nr3x )
    !-----------------------------------------------------------------------
    !! Uses \(\text{nkqs}\) and \(\text{index_sym}\) from module \(\texttt{exx}\),
    !! computes \(\text{rir}\).
    !
    USE symm_base,  ONLY : nsym, s, ft
    USE control_kcw, ONLY: rir
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN) :: nr1, nr2, nr3, nr1x, nr2x, nr3x 
    !
    ! ... local variables
    !
    INTEGER :: ikq, isym, i,j,k, ri,rj,rk, ir, nxxs
    INTEGER, allocatable :: ftau(:,:), s_scaled(:,:,:)
    !
    nxxs = nr1x*nr2x*nr3x
    !
    IF (.NOT. ALLOCATED(rir)) THEN   
        ALLOCATE( rir(nxxs,nsym) )
    ELSEIF ((SIZE(rir,1) /= nxxs) ) THEN 
        DEALLOCATE( rir )
        ALLOCATE( rir(nxxs,nsym) )
    ENDIF
    !
    rir = 0
    ALLOCATE ( ftau(3,nsym), s_scaled(3,3,nsym) )
    CALL scale_sym_ops (nsym, s, ft, nr1, nr2, nr3, s_scaled, ftau)
    DO isym = 1, nsym
       DO k = 1, nr3
          DO j = 1, nr2
             DO i = 1, nr1
                CALL rotate_grid_point( s_scaled(1,1,isym), ftau(1,isym), &
                     i, j, k, nr1, nr2, nr3, ri, rj, rk )
                ir = i + (j-1)*nr1x + (k-1)*nr1x*nr2x
                rir(ir,isym) = ri + (rj-1)*nr1x + (rk-1)*nr1x*nr2x
             ENDDO
          ENDDO
       ENDDO
    ENDDO
    !
    DEALLOCATE ( s_scaled, ftau )
    !
  END SUBROUTINE kcw_set_symm
