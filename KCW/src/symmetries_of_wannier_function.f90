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
  USE control_kcw,           ONLY : ir_end, num_wann
  USE control_kcw,           ONLY : spin_component
  USE mp_bands,              ONLY : root_bgrp, intra_bgrp_comm
  USE gvect,                 ONLY : ig_l2g, mill
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
  USE control_kcw,           ONLY : nsym_w, s_w, ft_w, centers
  USE interpolation,         ONLY : read_wannier_centers
  !
  IMPLICIT NONE 
  !
  COMPLEX(DP)              ::IMAG
  INTEGER                  :: iq, iwann, isym, ir, iq_, ig
  INTEGER                  :: iRq
  COMPLEX(DP), ALLOCATABLE :: rhowann(:,:,:)
  COMPLEX(DP), ALLOCATABLE :: rhowann_aux(:)
  COMPLEX(DP), ALLOCATABLE :: rhowann_aux2(:)
  COMPLEX(DP), ALLOCATABLE :: rho_rotated(:)
  COMPLEX(DP), ALLOCATABLE :: rhog(:)
  CHARACTER (LEN=256)      :: file_base
  REAL(DP)                 :: Gvector(3)
  REAL(DP)                 :: xk_cryst(3,1)
  CHARACTER (LEN=6), EXTERNAL :: int_to_char
  COMPLEX(DP), ALLOCATABLE :: phase(:)
  LOGICAL                  :: onlyone
  REAL(DP), ALLOCATABLE    :: cx(:)
  REAL(DP)                 :: ggg(3)
  REAL(DP)                 :: ft_cart(3)
  LOGICAL                  :: is_satisfied
  !
  !
  IMAG = CMPLX(0.D0, 1.D0, kind=DP)
  ALLOCATE ( rhog (ngms) )
  ALLOCATE (rhowann ( dffts%nnr, nkstot/nspin, num_wann) )
  ALLOCATE ( rhowann_aux(dffts%nnr) )
  ALLOCATE ( rhowann_aux2(dffts%nnr) )
  ALLOCATE ( rho_rotated(dffts%nnr) )
  ALLOCATE( phase (dffts%nnr) )
  ALLOCATE ( nsym_w(num_wann) )
  ALLOCATE (s_w(3,3,48,num_wann))
  ALLOCATE (ft_w(3,48,num_wann))
  ALLOCATE( centers(3,num_wann) )
  ALLOCATE(cx(3))
  !
  WRITE(*,*) "Inside symmetry function"
  WRITE(*,*) "nkstot= ", nkstot, "nsym= ", nsym, "num_wann= ", num_wann
  !
  WRITE(*,*) "read_wannier_centers"
  !
  !get wannier centres in lattice coordinates
  !
  CALL read_wannier_centers()
  DO iwann=1, num_wann
    WRITE(*,*) centers(:, iwann)
  END DO
  WRITE(*,*) "read_wannier_centers"
  CALL cryst_to_cart( num_wann, centers, at, +1 )
  
  WRITE(*,*) " SYMMETRIES OF WANNIER FUNCTION"
  DO isym=1, nsym
    iwann=5
    ft_cart(:) = ft(:, isym)
    CALL cryst_to_cart(1, ft_cart, at, +1)
    cx = MATMUL(centers(:,iwann), sr(:,:,isym))-ft_cart(:)
    is_satisfied = ( (centers(1,iwann)-cx(1))**2 + (centers(2,iwann)-cx(2))**2 &
                  + (centers(3,iwann)-cx(3))**2 .lt. 1.D-02 )
    IF( is_satisfied)  WRITE(*,*) isym
  END DO
 
  ! construct rir
  !
  CALL kcw_set_symm( dffts%nr1,  dffts%nr2,  dffts%nr3, &
  dffts%nr1x, dffts%nr2x, dffts%nr3x )  
  nsym_w = 0
  !
  ! read all the wannier densities for all the q. store it in rhowann
  DO iq_ = 1, nkstot
    IF ( lsda .AND. isk(iq_) /= spin_component) CYCLE
    iq = iq_-(spin_component-1)*nkstot/nspin
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
  !
  ! check which symmetries are satisfied by rhowann(:,:, iwann)
  !
  DO iwann=1, num_wann
    WRITE(*,*)
    WRITE(*,*)
    WRITE(*,*)
    WRITE(*,*)
    WRITE(*,*)
    WRITE(*,*)
    WRITE(*,*)
    
    DO isym=1, nsym
      DO iq_=1, nkstot
        IF ( lsda .AND. isk(iq_) /= spin_component) CYCLE
        iq = iq_-(spin_component-1)*nkstot/nspin
        rhowann_aux = 0.D0
        !WRITE(*,*) "isym", isym, "ft", ft(:, isym)
        !
        ! calculate rho_rotated = rho_q(R^{-1}.r-f)*EXP(-i k.f)
        !
        CALL rotate_evc(isym, rhowann(:,iq,iwann), rho_rotated)
        xk_cryst(:,1)=xk(:,iq)
        CALL cryst_to_cart(1,xk_cryst,at,-1)
        rho_rotated(:) = rho_rotated(:)*EXP(-IMAG*tpi*dot_product(xk_cryst(:,1),ft(:,isym)))
        !
        ! rotate q point
        !
        CALL rotate_xk(iq, isym, iRq, Gvector)
        !WRITE(*,*) "k=", xk(:,iq), "Rk=", xk(:, iRq), "Gvector", Gvector, "isym=", isym
        !
        ! compare the two rho, in rho_Rq we apply the phase factor:
        !             rho_{Rq}(r) = rho_{q'+G}(r) EXP(-iGr)
        !
        CALL calculate_phase(Gvector, phase)
        
        
        onlyone = .false.
        rhowann_aux2(:) = rho_rotated(:)/rhowann(:,iRq,iwann)
        !
        !trial calculating fourier transform and checking which component is nonzero
           
        CALL fwfft('Rho', rhowann_aux2(:), dffts)
        ig = 1
        do ir=1, dffts%nnr
          if ( ABS(rhowann_aux2(ir)) .gt. 1.d-02 ) THEN
            if ( onlyone ) then
              onlyone = .false.
              exit
            else 
              onlyone = .true.
              ig = ir
            end if 
          end if
        end do
        !
        !trial calculating G from the exponential

        !DO ir = 1, dffts%nnr
        !  phase(ir) = EXP(-IMAG*tpi*dot_product(DBLE(mill(:,ig)),r(ir,:)))
        !END DO

        rhowann_aux(:) = rho_rotated(:) - phase(:)*rhowann(:,iRq,iwann) 
        rhowann_aux(:) = rhowann_aux(:)/dffts%nnr

        !WRITE(*,*) "iwann, isym, iq", iwann, isym, iq, SUM( ABS(rhowann_aux(:)) )
        !
        !WRITE(*,*) "iwann: ", iwann, "iq", iq, "iRq", iRq, SUM (ABS(rhowann_aux(:)))
        !IF (SUM (ABS(rhowann_aux(:))) .lt. 1.D-01) THEN
        !  WRITE(*,*) "YES Gvector", Gvector, "isym=", isym, "onlyone ", onlyone, &
        !  SUM( ABS(rhowann_aux(:)) )
        !  !"ggg", ggg,&
        !  !"g", mill(:,ig)!dffts%nl(ig))
        !  !WRITE(*,*) "iwann", iwann, "iq", iq, "respected symmetry", isym,"."
        !ELSE 
        !  WRITE(*,*) "NO  Gvector", Gvector, "isym=", isym, "onlyone ", onlyone, &
        !  SUM( ABS(rhowann_aux(:)) )   !"ggg", ggg,&
        !     !"g", mill(:,ig)!dffts%nl(ig))
        !END IF
        WRITE(*,*) "iq", iq, "isym", isym, "iwann", iwann, "SUM.", SUM( ABS(rhowann_aux(:)))
        IF ( SUM( ABS(rhowann_aux(:)) ) .gt. 1.D-02 )  THEN
          EXIT
        END IF
        IF( iq .eq. nkstot/nspin ) THEN
          nsym_w(iwann) = nsym_w(iwann) + 1
          s_w(:,:,nsym_w(iwann),iwann) = s(:,:,isym)
          ft_w(:, nsym_w(iwann), iwann) = ft(:, isym)
          WRITE(*,*) "iwann = ", iwann, " respected symmetry ", isym
        ENDIF
      END DO!isym
    END DO!iq
    WRITE(*,*)
    WRITE(*,*) "TOTAL NUMBER OF SYMMETRIES OF WF ", iwann, " = ", nsym_w(iwann)
  END DO !iwann 
END SUBROUTINE




SUBROUTINE rotate_xk(ik_ToRotate, isym, iRk, Gvector)
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
  !
  IMPLICIT NONE 
  !
  INTEGER,  INTENT(IN)  :: ik_torotate
  !k point to rotate
  INTEGER,  INTENT(IN)  :: isym
  !symmetry to apply
  INTEGER,  INTENT(OUT) :: iRk
  !index of Rxk in FBZ
  REAL(DP), INTENT(OUT) :: Gvector(3)
  ! G vector connecting xk and Rxk (in cartesian!)
  INTEGER               :: ik
  REAL(DP)              :: Rxk(3)
  ! Rotated k vector Rxk = R.xk(ik) in cartesian
  !
  ! rotate k point
  !
  Rxk(:) = MATMUL(xk(:, ik_torotate), sr(:,:,isym))
  !
  ! Find k point that differs a reciprocal lattice vector from the rotated one
  !
  CALL find_kpoint(Rxk, iRk, Gvector)
  !
  ! WRITE(*,*) "Rxk:", Rxk, "xk:", xk(:, iRk), "Gvector", delta_xk
  !WRITE(*,*) "ik", ik, "xk", xk(:, ik)
END SUBROUTINE 
    
SUBROUTINE find_kpoint(xk_, ik_, Gvector)
  USE kinds,                 ONLY : DP
  USE io_global,             ONLY : stdout
  USE klist,                 ONLY:  xk, nkstot        !k points in cartesian coordinates
  USE cell_base,             ONLY : at
  !
  IMPLICIT NONE
  !
  INTEGER                  :: ik
  !indices
  REAL(DP), INTENT(IN)     :: xk_(3)
  ! k point to find (in cartesian coordinates)
  INTEGER,  INTENT(OUT)    :: ik_
  !index of the k point xk_ in xk
  REAL(DP), INTENT(OUT)    :: Gvector(3)
  ! G connecting xk_ and xk(ik_):
  !      G = xk_ - xk(ik_)
  REAL(DP)              :: delta_xk(3)
  ! differenze of Rxk and xk in cartesian
  REAL(DP)              :: delta_xk_cryst(3)
  ! difference of Rxk and xk in crystal
  LOGICAL               :: isGvec
  !
  Gvector(:) = 0.5 ! Impossible condition, this to check in the end this value changed
  DO ik = 1, nkstot
    !
    !delta_xk in cartesian
    !
    delta_xk(:) = xk_(:) - xk(:, ik)
    !
    !go to crystal
    !
    delta_xk_cryst(:) = at(1,:)*delta_xk(1) + at(2,:)*delta_xk(2) + &
                        at(3,:)*delta_xk(3)
    IF( isGvec(delta_xk_cryst) ) THEN 
      Gvector(:) = delta_xk(:)!delta_xk_cryst(:)    
      ik_ = ik
      EXIT
    END IF
  END DO!ik
  !
  IF ( ANY( Gvector .EQ. 0.5 ) ) THEN 
    WRITE(stdout, *) "ERROR! Could not find  k point", xk_, "in k mesh grid."
  END IF
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
      
