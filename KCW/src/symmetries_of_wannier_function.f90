!
!Giovanni Cistaro
!
SUBROUTINE symmetries_of_wannier_function()
! For all symmetry operations, we check if the following equation
! is satisfied:
!        \rho^{0n}_{Rq}(r) = \rho^{0n}_{q}(R^{-1}.r-f) e^{-iqf}
! Those are also the symmetries satisfied by the perturbation
! \Delta V^{0n}, thus the symmetries to satisfy to reduce the q-mesh. 
!
  USE kinds,                 ONLY : DP
  USE constants,             ONLY : tpi
  USE control_kcw,           ONLY : io_real_space, r
  USE control_kcw,           ONLY : tmp_dir_kcwq, tmp_dir_kcw
  USE control_kcw,           ONLY : ir_end, num_wann, Rvect, x_q
  USE control_kcw,           ONLY : spin_component, nqstot, kcw_iverbosity
  USE mp_bands,              ONLY : root_bgrp, intra_bgrp_comm
  USE gvect,                 ONLY : ig_l2g, mill, g
  USE gvecs,                 ONLY : doublegrid, ngms
  USE fft_base,              ONLY : dffts
  USE control_flags,         ONLY : gamma_only
  USE klist,                 ONLY : xk     
  USE symm_base,             ONLY : ft, nsym, s, sr
  USE cell_base,             ONLY : omega
  USE io_kcw,                ONLY : read_rhowann, read_rhowann_g
  USE fft_interfaces,        ONLY : invfft, fwfft
  USE klist,                 ONLY : nkstot
  USE lsda_mod,              ONLY : lsda, isk, nspin, current_spin
  USE cell_base,             ONLY : bg, at
  USE control_kcw,           ONLY : nsym_w_k, nsym_w_q, s_w, ft_w, & 
                                    centers, check_rvect, &
                                    sym_only_for_q !, shift_centers
  USE interpolation,         ONLY : read_wannier_centers
  USE io_global,             ONLY : stdout
  USE mp,                    ONLY : mp_sum  
  IMPLICIT NONE 
  !
  COMPLEX(DP)              ::IMAG
  INTEGER                  :: iq, iwann, isym, ir, iq_, ig
  INTEGER                  :: iRq
  COMPLEX(DP), ALLOCATABLE :: rhowann(:,:,:)
  COMPLEX(DP), ALLOCATABLE :: rhowann_aux(:)
  COMPLEX(DP), ALLOCATABLE :: rho_rotated(:)
  COMPLEX(DP), ALLOCATABLE :: rhog(:)
  CHARACTER (LEN=256)      :: file_base
  REAL(DP)                 :: Gvector(3), Gvector_cryst(3)
  CHARACTER (LEN=6), EXTERNAL :: int_to_char
  COMPLEX(DP), ALLOCATABLE :: phase(:)
  COMPLEX(DP), ALLOCATABLE :: rhowann_(:,:,:)
  REAL(DP), ALLOCATABLE    :: cx(:)
  REAL(DP)                 :: x_qG_cryst(3)
  REAL(DP)                 :: x_q_cryst(3)
  REAL(DP)                 :: ft_cart(3)
  LOGICAL                  :: is_satisfied
  COMPLEX(DP)                 :: delta_rho
  COMPLEX(DP)                 :: delta_rho_R
  COMPLEX(DP) :: sh
  COMPLEX (DP) :: int_rho_Rq, eiRqR
  !
  !
  CALL start_clock ( 'check_symm' )
  !
  IMAG = CMPLX(0.D0, 1.D0, kind=DP)
  ALLOCATE ( rhog (ngms) )
  ALLOCATE (rhowann ( dffts%nnr, nkstot/nspin, num_wann) )
  ALLOCATE ( rhowann_aux(dffts%nnr) )
  ALLOCATE ( rhowann_(dffts%nnr,nqstot,num_wann) )
  ALLOCATE ( rho_rotated(dffts%nnr) )
  ALLOCATE( phase (dffts%nnr) )
  ALLOCATE ( nsym_w_k(num_wann) )
  ALLOCATE ( nsym_w_q(num_wann) )
  ALLOCATE ( sym_only_for_q(48, num_wann) )
  ALLOCATE (s_w(3,3,48,num_wann))
  ALLOCATE (ft_w(3,48,num_wann))
  ALLOCATE( centers(3,num_wann) )
  ALLOCATE(cx(3))
  !
  WRITE( stdout, '(5X, "SYM : Checking Symmetry of the WFs")')
  WRITE( stdout, '(7X, "SYM : nkstot=", I5, 3X, "nsym tot=", I5, 3X, "num_wann=", I5)') nkstot, nsym,num_wann
  !
  !get wannier centres in lattice coordinates
  WRITE(stdout, '(7X, "SYM : read_wannier_centers ...")', advance='no')
  CALL read_wannier_centers()
  WRITE(stdout, '(" DONE")')
  !
  !go to crystal coordinates
  !
  !CALL cryst_to_cart( num_wann, centers, at, +1 )
  WRITE(stdout, '(13X, "Centers of wannier functions (crys)...")')
  DO iwann=1, num_wann
    WRITE(stdout, '(13X, "iwann=", I5, 3X, "centers = (", 3(F20.12, ","), ")" )') &
           iwann, centers(1:3, iwann) 
  END DO
  !
  ! construct rir
  !
  CALL kcw_set_symm( dffts%nr1,  dffts%nr2,  dffts%nr3, &
  dffts%nr1x, dffts%nr2x, dffts%nr3x )  
  nsym_w_k = 0
  nsym_w_q = 0
  !
  ! read all the wannier densities for all the q. store it in rhowann
  ! This is not ported yet to noncollinear case (only one component for now)
  DO iq = 1, nqstot
    !IF ( lsda .AND. isk(iq_) /= spin_component) CYCLE
    !iq = iq_-(spin_component-1)*nkstot/nspin
    !
    !name of the folder with q point iq
    !
    tmp_dir_kcwq= TRIM (tmp_dir_kcw) // 'q' &
        & // TRIM(int_to_char(iq))//'/'
    !
    DO iwann=1, num_wann
      !
      !read density rhowann from file, store it in rho_iwann
      !
      IF ( .NOT. io_real_space ) THEN 
        !
        file_base=TRIM(tmp_dir_kcwq)//'rhowann_g_iwann_'//TRIM(int_to_char(iwann))
        CALL read_rhowann_g( file_base, &
             root_bgrp, intra_bgrp_comm, &
             ig_l2g, 1, rhog(:), .FALSE., gamma_only )
        rhowann_aux=(0.d0,0.d0)
        rhowann_aux(dffts%nl(:)) = rhog(:)
        CALL invfft ('Rho', rhowann_aux, dffts)
        rhowann(:,iq, iwann) = rhowann_aux(:)*omega
        !
      ELSE 
        !
        file_base=TRIM(tmp_dir_kcwq)//'rhowann_iwann_'//TRIM(int_to_char(iwann))
        CALL read_rhowann( file_base, dffts, rhowann_aux )
        rhowann(:,iq, iwann) = rhowann_aux(:)
        !
      ENDIF
      !
      !end of storing rhowann
      !
    END DO!iwann
  END DO !iq
  !copy variable rhowann so that we can shift the center
  !
  rhowann_(:,:,:) = rhowann(:,:,:)
  !
  !density we will use to check symmetries
  !
!  ! Shift center to detect more symmetries. NOT WORKING PROPERLY 
!  ! Kept here for reference 
!  !
!  IF (shift_centers) THEN
!    DO iwann = 1, num_wann
!      DO iq = 1, nqstot
!        x_q_cryst(:) = xk(:,iq)
!        CALL cryst_to_cart(1, x_q_cryst, at, -1)
!        !
!        !go to G space
!        !
!        rhowann_(:,iq,iwann) = rhowann_(:,iq,iwann)!/omega
!        CALL fwfft ('Rho', rhowann_(:,iq,iwann), dffts)  
!        rhog(:) = rhowann_(dffts%nl(:),iq,iwann)
!        !
!        !apply shift
!        !
!        DO ig = 1, ngms
!          x_qG_cryst(:) = g(:,ig)
!          CALL cryst_to_cart(1, x_qG_cryst, at, -1)
!          x_qG_cryst(:) = x_qG_cryst(:) + x_q_cryst(:)
!          rhog(ig) = rhog(ig)*EXP( -IMAG*tpi*DOT_PRODUCT(x_qG_cryst(:),centers(:,iwann)) )
!        END DO
!        !
!        !back to r space
!        !
!        rhowann_(:,iq,iwann) = 0.D0
!        rhowann_(dffts%nl(:),iq,iwann) = rhog(:)
!        CALL invfft ('Rho', rhowann_(:,iq,iwann), dffts)
!      END DO!iq
!    END DO!iwann
!  ENDIF!shift_centers
  !
  ! check which symmetries are satisfied by rhowann(:,:, iwann)
  !
  DO iwann=1, num_wann
    !
    WRITE(stdout, '(/, 7X, "SYM : Checking WF #", I5)') iwann
    !
    DO isym=1, nsym
      !
      sym_only_for_q(nsym_w_q(iwann) + 1, iwann) = .FALSE.
      !
      DO iq = 1, nqstot
        ir = 1
        rhowann_aux = 0.D0
        !WRITE(*,*) "isym", isym, "ft", ft(:, isym)
        !
        ! calculate rho_rotated = rho_q(R^{-1}.r-f)*EXP(-i k.f)
        !
        CALL rotate_evc(isym, rhowann_(:,iq,iwann), rho_rotated)
        x_q_cryst(:)=xk(:,iq)
        CALL cryst_to_cart(1,x_q_cryst,at,-1)
        rho_rotated(:) = rho_rotated(:)*EXP(-IMAG*tpi*dot_product(x_q_cryst(:),ft(:,isym)))
        !
        ! rotate q point
        !
        CALL rotate_xk(iq, isym, iRq, Gvector, Gvector_cryst)
        IF ( ANY( Gvector_cryst .EQ. 0.5 ) ) THEN 
          EXIT
        END IF
      
        !WRITE(*,*) "k=", xk(:,iq), "Rk=", xk(:, iRq), "Gvector", Gvector, "isym=", isym
        !
        ! compare the two rho, in rho_Rq we apply the phase factor:
        !             rho_{Rq}(r) = rho_{q'+G}(r) EXP(-iGr)
        !
        CALL calculate_phase(Gvector, phase)
        !
        !
        rhowann_aux(:) = rho_rotated(:) - phase(:)*rhowann_(:,iRq,iwann) 
        rhowann_aux(:) = rhowann_aux(:)
        !
        ! integrate difference and normalize with respect to number of r points in the grid
        !delta_rho = SUM( ABS(rhowann_aux(:)) )/(dffts%nr1*dffts%nr2*dffts%nr3)
        delta_rho = SUM( (rhowann_aux(:)) )/(dffts%nr1*dffts%nr2*dffts%nr3)
        CALL mp_sum (delta_rho, intra_bgrp_comm)
        !
        int_rho_Rq = SUM( phase(:)*rhowann_(:,iRq,iwann)  ) / (dffts%nr1*dffts%nr2*dffts%nr3)
        CALL mp_sum (int_rho_Rq, intra_bgrp_comm)
        !
        !
        IF (check_rvect .AND. ABS(delta_rho) .gt. 1D-02) THEN 
          ! Try with the same Wannier in different cells:
          ! Each q contribution to the self-Hxc or to the screened self-Hxc (i.e. the alpha coeff)
          ! does not depend on the center of the Wannier density contribution at q. 
          ! This means we can use also the symmetries that send 
          ! \rho_q^{0n}(R^-1r -f) in \rho_Rq^{Ln}(r) = e^{-i(Rq).L}\rho_Rq^{0n}(r)
          ! with L any lattice vector in the SC to reduce the number of q points. 
          !
          ! delta_rho_R = int [\rho_q^{0n}(R^-1r -f) - e^{-i(Rq).L}\rho_Rq^{0n}(r)] =
          !             = int [\rho_q^{0n}(R^-1r -f) - \rho_Rq^{0n}(r)] + (1- e^{-i(Rq).L}) * \int [\rho_Rq^{0n}(r)]
          ! 
          ! NB: this is not true for the density response at q. For the symmetrization of the density response we must
          ! use only the "real" symmetry of the wannier density.
          ! sym_only_for_q store information of wheter the symmetry is a real one (FALSE) or if its an "extra" one to 
          ! be used only for the reduction of the q points (TRUE)
          !
           DO ir = 1, nkstot/nspin
             eiRqR=EXP( -IMAG*tpi*DOT_PRODUCT(x_q(:,iRq),Rvect(:,ir)) )
             delta_rho_R = delta_rho + (CMPLX(1.D0,0.D0, kind=DP) - eiRqR)*int_rho_Rq
             !WRITE (stdout, *) "          ir  =", ir,  "Rvect  =", Rvect(1:3,ir)
             !WRITE (stdout, *) "          iRq =", iRq, "x(iRq) =", x_q(1:3,iRq)
             !WRITE (stdout, *) "          \int rho_Rq(r) dr = ", int_rho_Rq, "exp(-iRq*Rvec) =", eiRqR
             !WRITE (stdout,'(10X, "ir=", I5, 3X, "SUM =", 2F20.12)')&
             !    ir, delta_rho 
             !WRITE (stdout,'(10X, "ir=", I5, 3X, "SUM =", 2F20.12)')&
             !    ir, delta_rho_R
             !WRITE (stdout,*)
             IF (ABS(delta_rho_R) .lt. 1D-02) THEN 
                delta_rho = delta_rho_R
                sym_only_for_q(nsym_w_q(iwann) + 1, iwann) = .TRUE.
                EXIT 
             ENDIF 
           ENDDO
        ENDIF
        ! 
        IF (kcw_iverbosity .gt. 2 ) & 
           !WRITE(stdout,'(7X, "iq=", I5, 3X, "isym =", I5, 3X, "iwann =", I5, 3X, "SUM =", F20.12)')&
           !   iq, isym, iwann, delta_rho
           WRITE(stdout,'(7X, "iq=", I5, 3X, "isym =", I5, 3X, "iwann =", I5, 3X, "SUM =", 2F20.12, 3X, "(ir =", I5 " )")')&
                 iq, isym, iwann, REAL(delta_rho), AIMAG(delta_rho), ir
        !
        IF ( ABS(REAL(delta_rho)) .gt. 1.D-02 .OR. ABS(AIMAG(delta_rho)) .gt. 1D-02)  THEN 
           IF (kcw_iverbosity .gt. 2) WRITE(stdout, '(13X, "isym =", I5, 3X, "NOT RESPECTED, skipping")') isym
           EXIT
        ENDIF
        !
        IF( iq .eq. nqstot ) THEN
           IF(.NOT. sym_only_for_q(nsym_w_q(iwann) + 1, iwann) ) THEN
             nsym_w_k(iwann) = nsym_w_k(iwann) + 1
           END IF 
           nsym_w_q(iwann) = nsym_w_q(iwann) + 1
           s_w(:,:,nsym_w_q(iwann),iwann) = s(:,:,isym)
           ft_w(:, nsym_w_q(iwann), iwann) = ft(:, isym)
           IF (kcw_iverbosity .gt. 1 ) THEN 
              IF (ir /= 1) THEN 
                 WRITE(stdout, '(13X, "isym =", I5, 3X,"    RESPECTED - only q")') isym
              ELSE
                 WRITE(stdout, '(13X, "isym =", I5, 3X,"    RESPECTED")') isym
              ENDIF
           ENDIF
        ENDIF
      END DO!isym
    END DO!iq
    WRITE(stdout,'(/, 13X, "TOTAL NUMBER OF RESPECTED SYMMETRIES (k and q)= ", I5)') nsym_w_k(iwann)
    IF (check_rvect) WRITE(stdout,'(/, 13X, "TOTAL NUMBER OF RESPECTED SYMMETRIES ( only q)= ", I5)') nsym_w_q(iwann)
  END DO !iwann 
  !
  DEALLOCATE(rhowann_)
  !
  CALL stop_clock ( 'check_symm' )
  !
END SUBROUTINE symmetries_of_wannier_function




SUBROUTINE rotate_xk(ik_ToRotate, isym, iRk, Gvector, Gvector_cryst)
! This subroutine rotates the vector with index ik_ToRotate with symmetry isym.
! Then, it finds which of the k vectors in the mesh matches the rotated one, 
! in the sense that they difference is a reciprocal lattice vector. 
! Equation to satisfy:
!         xk(:, iRk) = sr(:,:,isym).xk(:, ik_ToRotate) + Gvector(:)
  USE kinds,                 ONLY : DP
  USE io_global,             ONLY : stdout
  USE symm_base,             ONLY:  sr                !symmetry operation in cartesian coordinates
  USE symm_base,             ONLY:  invs              !index of inverse of symmetry op
  USE klist,                 ONLY:  xk, nkstot        !k points in cartesian coordinates
  USE cell_base,             ONLY : at
  USE control_kcw,           ONLY : x_q
  !
  IMPLICIT NONE 
  !
  INTEGER,  INTENT(IN)  :: ik_torotate
  !k point to rotate
  INTEGER,  INTENT(IN)  :: isym
  !symmetry to apply
  INTEGER,  INTENT(OUT) :: iRk
  !index of Rxk in FBZ
  REAL(DP), INTENT(OUT) :: Gvector(3), Gvector_cryst(3)
  ! G vector connecting xk and Rxk (in cartesian!)
  INTEGER               :: ik
  REAL(DP)              :: Rxk(3)
  ! Rotated k vector Rxk = R.xk(ik) in cartesian
  !
  ! rotate k point
  !
  Rxk(:) = MATMUL(x_q(:, ik_torotate), sr(:,:,isym))
  !
  ! Find k point that differs a reciprocal lattice vector from the rotated one
  !
  CALL find_kpoint(Rxk, iRk, Gvector, Gvector_cryst)
  !
  ! WRITE(*,*) "Rxk:", Rxk, "xk:", xk(:, iRk), "Gvector", Gvector
  !WRITE(*,*) "ik", ik, "xk", xk(:, ik)
END SUBROUTINE 
    
SUBROUTINE find_kpoint(xk_, ik_, Gvector, Gvector_cryst)
  USE kinds,                 ONLY : DP
  USE io_global,             ONLY : stdout
  USE klist,                 ONLY:  xk, nkstot        !k points in cartesian coordinates
  USE cell_base,             ONLY : at
  USE lsda_mod,              ONLY : nspin
  USE control_kcw,           ONLY : x_q
  !
  IMPLICIT NONE
  !
  INTEGER                  :: ik
  !indices
  REAL(DP), INTENT(IN)     :: xk_(3)
  ! k point to find (in cartesian coordinates)
  INTEGER,  INTENT(OUT)    :: ik_
  !index of the k point xk_ in xk
  REAL(DP), INTENT(INOUT)    :: Gvector(3)
  ! G connecting xk_ and xk(ik_):
  !      G = xk_ - xk(ik_)
  REAL(DP), INTENT(INOUT)    :: Gvector_cryst(3)
  ! Same as above but in crystal coordinates
  REAL(DP)              :: delta_xk(3)
  ! differenze of Rxk and xk in cartesian
  REAL(DP)              :: delta_xk_cryst(3)
  ! difference of Rxk and xk in crystal
  LOGICAL               :: isGvec
  !
  Gvector_cryst(:) = 0.5 ! Impossible condition. If not changed at the end of this routine, something wrong
  DO ik = 1, nkstot/nspin
    !
    !delta_xk in cartesian
    delta_xk(:) = xk_(:) - x_q(:, ik)
    !
    !go to crystal
    delta_xk_cryst(:)=delta_xk(:)
    CALL cryst_to_cart(1, delta_xk_cryst, at, -1)
    !
    IF( isGvec(delta_xk_cryst) ) THEN 
      Gvector(:)       = delta_xk(:)
      Gvector_cryst(:) = delta_xk_cryst(:)
      ik_ = ik
      EXIT
    END IF
  END DO!ik
  !
END SUBROUTINE
      
    
    
    
    
!copy from polaron.f90 of epw
!-----------------------------------------------------------------------
FUNCTION isGVec(xxk)
!-----------------------------------------------------------------------
!! Return true if xxk integer times of the reciprocal vector
!! if xxk is the difference of two vectors, then return true if these
!! two vector are the same
!-----------------------------------------------------------------------
!
  USE kinds,                 ONLY : DP
  IMPLICIT NONE
  !
  REAL(KIND = DP), INTENT(in) :: xxk(3)
  !! FIXME
  !
  ! Local variable
  LOGICAL :: isGVec
  !! FIXME
  !
  isGVec = &
    ABS(xxk(1) - NINT(xxk(1))) < 1.e-6 .AND. &
    ABS(xxk(2) - NINT(xxk(2))) < 1.e-6 .AND. &
    ABS(xxk(3) - NINT(xxk(3))) < 1.e-6
  !-----------------------------------------------------------------------
END FUNCTION isGVec
!-----------------------------------------------------------------------
      
