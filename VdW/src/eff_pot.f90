!--------------------------------------------------------------------
SUBROUTINE eff_pot (rho, nspin, alat, omega, charge, vstart, thresh_veff)
  !--------------------------------------------------------------------
  !
  !     Effective  potential (V_eff) in TF+vW scheme
  !
  USE kinds,                ONLY : DP
  USE constants,            ONLY : fpi, e2
  USE io_global,            ONLY : stdout
  USE io_files,             ONLY : prefix, nwordwfc, iunwfc
  USE klist,                ONLY : nelec
  USE fft_base,             ONLY : dfftp
  USE fft_interfaces,       ONLY : fwfft, invfft
  USE gvect,                ONLY : nlm, g, nl, ngm, gg, gstart
  USE wvfct,                ONLY : g2kin, wg, nbndx, et, nbnd, npwx, &
                                   ecutwfc, igk, npw
  USE uspp,                 ONLY : nkb
  USE scf,                  ONLY : v, vltot, vrs, rho_core
  USE gvecs,              ONLY : nls, nlsm, doublegrid
  USE eff_v,                ONLY : rho_fft, veff
  USE mp_global,            ONLY : intra_pool_comm
  USE mp,                   ONLY : mp_sum

  USE control_flags,        ONLY : gamma_only
!  USE control_vdw,          ONLY : thresh_veff
!  USE wavefunctions_module, ONLY : evc, psic
  IMPLICIT NONE
  !
  !    input
  !
  INTEGER :: nspin 
  real(kind=DP) :: rho(dfftp%nnr, nspin), alat, omega, charge, charge_old
  !
  !    output
  !
  real(kind=DP), ALLOCATABLE :: vv (:,:), rho_in(:,:), k_gamma(:)
  !
  !    local variables
  !
  real(kind=DP) :: tpiba2, fac
  real(kind=DP), ALLOCATABLE :: S(:), psi_smooth(:)
  COMPLEX(kind=DP), ALLOCATABLE ::  aux (:), aux1 (:), psi (:)
  COMPLEX(kind=DP), ALLOCATABLE :: ws_psi(:,:), ws_hpsi(:,:), ws_psic(:)

  INTEGER :: ir, is, ig
  real (kind=DP) :: avg1, avg2, eps
  INTEGER :: nnn, nite, ite = 0
  real (kind=DP) :: vstart, thresh_veff
  CHARACTER (len=10):: str_ite
  real(kind=DP), EXTERNAL :: qe_erf
  !
  CALL start_clock('eff_pot')
  !
  ALLOCATE ( aux(dfftp%nnr), aux1(ngm), psi(dfftp%nnr), psi_smooth(dfftp%nnr), S(dfftp%nnr) )
  ALLOCATE ( vv(dfftp%nnr, nspin) )
  ALLOCATE ( rho_in(dfftp%nnr, nspin) )
  !
  tpiba2 = (fpi / 2.d0 / alat) **2
  !
  rho_in = abs(rho)
  !
  !  set value for psi in real space
  !
  psi(:) = CMPLX ( sqrt( rho_in (:,1) ), 0.0_dp, KIND=dp )
  !
  !  bring psi to G space
  !
  CALL fwfft ('Dense', psi, dfftp)
  !
  !  extract the smooth part of psi in G-space
  !
  DO ig = 1, ngm
    IF ( (tpiba2 * gg(ig)) > ecutwfc ) THEN
      psi(nl(ig)) = (0.0_dp, 0.0_dp)
    ENDIF
  ENDDO
  !
  aux = psi
  !
  !  bring psi_smooth to real space (approximation of psi)
  !
  CALL invfft ('Dense', psi, dfftp)
  psi_smooth(:) = DBLE ( psi(:) )
  deallocate ( psi )
  !
  !  check the difference of the total charge density
  !
  charge = 0.d0
  DO is = 1, nspin
     DO ir = 1, dfftp%nnr
        !  NB: the following will work only if nspin = 1 !
        charge = charge + abs( rho_in(ir,is)-psi_smooth(ir)**2 )
     END DO
  ENDDO
  charge = charge * omega / (dfftp%nr1*dfftp%nr2*dfftp%nr3) / nelec
#ifdef __MPI
  CALL mp_sum( charge, intra_pool_comm )
#endif
  WRITE( stdout, '(/,10x,"Charge difference due to FFT   ",f10.8)' ) charge
  !
  ! compute charge density using smooth wfc
  !
  DO is = 1, nspin
     DO ir = 1, dfftp%nnr
        rho_fft(ir,1) = psi_smooth(ir)**2
     ENDDO
  ENDDO
  !
  !  calculate P^2 |psi> in G-space (NB: V(G=0)=0 )
  !
  aux1(:) = 0.d0
  DO ig = 1, ngm
     fac =  gg(ig) * tpiba2
     aux1(ig) = fac * aux(nl(ig))
  ENDDO
  !
  aux(:) = 0.d0
  DO ig = 1, ngm
     aux(nl(ig)) = aux1(ig)
  ENDDO
  !
  IF (gamma_only) THEN
     DO ig = 1, ngm
        aux(nlm(ig)) = CONJG( aux1(ig) )
     ENDDO
  ENDIF
  !
  !      bring P^2 |psi>  to real space, kinetic term is kept
  !      in aux
  !
  CALL invfft ('Dense', aux, dfftp)
  !
  !
  !  compute V_eff  potential by FT and then use it as initial
  !  potential for iteration
  !
  avg1 = 0.d0
  avg2 = 0.d0
  nnn  = 0
  eps = 0.04d0
  is = 1
  IF (.false.) THEN
     DO ir = 1, dfftp%nnr
        vv(ir,is) = -DBLE(aux(ir))
     ENDDO
  ELSE
     DO ir = 1, dfftp%nnr
        IF (abs(psi_smooth(ir)) > eps ) THEN
           vv(ir,is) = -DBLE(aux(ir)) / psi_smooth(ir)
        ELSE
           avg1 = avg1 - DBLE(aux(ir))
           avg2 = avg2 + psi_smooth(ir)
           nnn = nnn + 1
        ENDIF
     ENDDO
#ifdef __MPI
     CALL mp_sum( avg1, intra_pool_comm )
     CALL mp_sum( avg2, intra_pool_comm )
     CALL mp_sum(nnn, intra_pool_comm)
#endif
     IF (nnn > 0 ) THEN
        DO ir = 1, dfftp%nnr
           IF (abs(psi_smooth(ir)) <= eps ) THEN
              vv(ir,is) = avg1 / avg2
           ENDIF
        ENDDO
     ENDIF
  ENDIF
  !
  ! uncomment the following loop will set local pot. as initial pot.
  !
vstart=20000
  DO ir = 1, dfftp%nnr
     vrs(ir,1) = v%of_r(ir,1) + vltot(ir)
     vv(ir,is) = qe_erf( abs(psi_smooth(ir)*dble(vstart)))*vv(ir,is) + &
           (1.d0-qe_erf( abs(psi_smooth(ir)*dble(vstart))))*vrs(ir,is)
  ENDDO
  !
  !  check the quality of trial potential
  !
  CALL check_v_eff(vv(1,1), charge)
  WRITE( stdout, '(/,10x,"Charge difference due to V_eff (initial)   ",f10.8,/)' ) charge
  !
  ! iterative procedure of V_eff if necessary
  !
  nite = 0
  charge_old = 1000.d0
  DO WHILE ( charge > thresh_veff )
     !
     CALL ite_veff_nhpsi( 1 )
     !
     !  check the quality of veff by solving equation
     !      p^2|phi> + (veff-mu)|phi> = 0
     !  whose the GS solution should be the square root
     !  of the charge density
     !
     CALL check_v_eff(vv, charge)
     nite = nite + 1
!#ifdef __MPI
!     call ireduce(1, nite)
!#endif
     WRITE( stdout, '(10x,"iter #   ", i3, "   charge diff.   ", f10.8, &
                     & "   thresh_veff   ", f10.8,/)' ) nite, charge, thresh_veff
     !
     IF ( charge_old < charge ) THEN
!        CALL io_pot( 1, TRIM( prefix )//'.veff', v, nspin )
!        CALL io_pot( 1, TRIM( prefix )//'.rho-coreff', rho_veff, nspin )
        WRITE(stdout, '( 10x, 10("*"), "unstability happens", 10("*") )' )
!        goto 100
     ENDIF
     !
     charge_old = charge
     !
  ENDDO
  !
  ! set the optmized eff. potential to veff
  !
100 CONTINUE
  veff(:,:) = vv(:,:)
  !
  DEALLOCATE ( vv, rho_in )
  DEALLOCATE (aux,aux1,psi_smooth,S)
  !
  CALL stop_clock('eff_pot')
  !
  RETURN
  !
  CONTAINS
     !
     !------------------------------------------------------------------------
     !
     SUBROUTINE ite_veff_nhpsi( nstep )
     !
     IMPLICIT NONE
     !
     INTEGER :: nstep, nveff
     !
     real (kind=DP) :: alp, beta, s2r2, sr2, s2r, sr, r2, s2, D, Da, Db, w1
     !
     ! Compute S(r) at first step
     !
     CALL start_clock ('ite_veff')
     !
!write (stdout,*) ' enter ite_veff_nhpsi'
!FLUSH( stdout )
     s2 = 0.d0
     DO ir = 1, dfftp%nnr
        S(ir) = psi_smooth(ir) * DBLE (aux(ir)) + vv(ir,1)*psi_smooth(ir)
        s2 = s2 + S(ir)**2
     ENDDO
#ifdef __MPI
     CALL mp_sum( s2, intra_pool_comm )
#endif
     !
     DO nnn = 1, nstep
        !
        !
        ! Compute alpha & beta in Veff = Veff + alp*S(r) + beta
        !
        s2r2 = 0.d0
        s2r  = 0.d0
        sr2  = 0.d0
        sr   = 0.d0
        r2   = 0.d0
        !
        DO ir = 1, dfftp%nnr
           r2   = r2   + psi_smooth(ir)**4
           s2r2 = s2r2 + ( S(ir) * psi_smooth(ir)**2 )**2
           sr2  = sr2  +   S(ir) * psi_smooth(ir)**4
           s2r  = s2r  + ( S(ir)**2) * psi_smooth(ir)**2
           sr   = sr   +   S(ir) * psi_smooth(ir)**2
        ENDDO
#ifdef __MPI
        CALL mp_sum( r2, intra_pool_comm )
        CALL mp_sum( s2r2, intra_pool_comm )
        CALL mp_sum( sr2, intra_pool_comm )
        CALL mp_sum( s2r, intra_pool_comm )
        CALL mp_sum( sr, intra_pool_comm )
#endif
        !
        D  = r2*s2r2 - sr2*sr2
        Da = sr*sr2  - s2r*r2
        Db = sr2*s2r - s2r2*sr
        !
        IF (D>0.d0) THEN
           alp = Da/D
           beta = Db/D
        ELSE
           WRITE(*,*) 'Det. of Hessian matrix is negative'
           STOP
        ENDIF
        !
!        if (mod(nnn,100) .eq. 0) then
!           write(*,*)'iteration ',nnn
!           write(*,*) 's2 = ',s2
!           write(*,*) 'D = ' , D
!           write(*,*) 'Da = ', Da
!           write(*,*) 'Db = ', Db
!           write(*,*) 'alp = ',alp
!           write(*,*) 'beta = ',beta
!           write(*,*)
!        endif
        !
        ! Update V-eff
        !
        DO ir = 1, dfftp%nnr
           vv(ir,1)= vv(ir,1) + alp*S(ir) + beta
           S(ir)   = S(ir) * (1.d0 + alp*psi_smooth(ir)**2) + &
                     beta*psi_smooth(ir)**2
        ENDDO
        !
        s2 = 0.d0
        DO ir = 1, dfftp%nnr
           s2 = s2 + S(ir)**2
        ENDDO
#ifdef __MPI
        CALL mp_sum( s2, intra_pool_comm )
#endif
        !
     ENDDO
     !
     CALL stop_clock ('ite_veff')
     !
     RETURN
     !
     END SUBROUTINE ite_veff_nhpsi
     !
END SUBROUTINE eff_pot

