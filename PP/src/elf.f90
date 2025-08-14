!
! Copyright (C) 2001-2021 Quantum ESPRESSO Foundation
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
SUBROUTINE do_elf_kin(fout,iself,spin_component)
  !-----------------------------------------------------------------------
  !
  !  Calculation of the electron localization function (iself = true)
  !  or the kinetic energy density (iself = false). In the case of the
  !  KED calculation, spin_component = 0 gives the total KED, and
  !  spin_component = 1 or 2 gives the contribution from the
  !  corresponding spin channel.
  !
  !     elf = 1/(1+d**2)
  !  where
  !     d = ( t(r) - t_von_Weizacker(r) ) / t_Thomas-Fermi(r)
  !  and
  !     t (r) = (hbar**2/2m) * \sum_{k,i} |grad psi_{k,i}|**2
  !             (kinetic energy density)
  !     t_von_Weizaecker(r) = (hbar**2/2m) * 0.25 * |grad rho(r)|**2/rho
  !             (non-interacting boson)
  !     t_Thomas-Fermi (r) = (hbar**2/2m) * 3/5 * (3*pi**2)**(2/3) * rho**(5/3)
  !             (free electron gas)
  !
  !  see also http://en.wikipedia.org/wiki/Electron_localization_function
  !
  USE kinds, ONLY: DP
  USE constants, ONLY: pi
  USE cell_base, ONLY: omega, tpiba
  USE fft_base,  ONLY: dffts, dfftp
  USE fft_interfaces, ONLY : fwfft, invfft, fft_interpolate
  USE gvect, ONLY: gcutm, g, ngm
  USE gvecs, ONLY : ngms, doublegrid, dual
  USE io_files, ONLY: restart_dir
  USE klist, ONLY: nks, xk, ngk, igk_k
  USE lsda_mod, ONLY: nspin, isk
  USE scf, ONLY: rho, rhoz_or_updw
  USE symme, ONLY: sym_rho, sym_rho_init
  USE wvfct, ONLY: nbnd, wg
  USE control_flags, ONLY: gamma_only
  USE wavefunctions,  ONLY: evc
  USE mp_pools, ONLY: inter_pool_comm, intra_pool_comm
  USE mp, ONLY: mp_sum
  USE pw_restart_new, ONLY : read_collected_wfc
  !
  ! I/O variables
  !
  IMPLICIT NONE
  REAL(DP) :: fout (dfftp%nnr)
  LOGICAL :: iself
  INTEGER :: spin_component
  !
  ! local variables
  !
  INTEGER :: i, j, k, ibnd, ik, is
  REAL(DP) :: gv(3), w1, d, fac
  REAL(DP), ALLOCATABLE :: kkin (:), tbos (:,:)
  COMPLEX(DP), ALLOCATABLE :: aux (:), aux2 (:)
  !
  CALL infomsg ('do_elf_kin', 'elf/kin + US not fully implemented')
  !
  ALLOCATE (kkin(dfftp%nnr))
  ALLOCATE (aux (dffts%nnr))
  aux(:) = (0.d0,0.d0)
  kkin(:) = 0.d0
  !
  ! Calculates local kinetic energy, stored in kkin
  !
  DO ik = 1, nks
     !
     !   reads the eigenfunctions
     !
     CALL read_collected_wfc ( restart_dir(), ik, evc )

     ! skip if this is a KED calculation with requested spin_component
     IF (.NOT.iself .AND. spin_component > 0 .AND. isk(ik) /= spin_component) CYCLE
     !
     DO ibnd = 1, nbnd
        DO j = 1, 3
           aux(:) = (0.d0,0.d0)
           w1 = wg (ibnd, ik) / omega
           DO i = 1, ngk(ik)
              gv (j) = (xk (j, ik) + g (j, igk_k (i,ik) ) ) * tpiba
              aux (dffts%nl(igk_k(i,ik) ) ) = cmplx(0d0, gv (j) ,kind=DP) * evc(i,ibnd)
              IF (gamma_only) THEN
                 aux (dffts%nlm(igk_k(i,ik) ) ) = cmplx(0d0, -gv (j) ,kind=DP) * &
                      conjg ( evc (i, ibnd) )
              ENDIF
           ENDDO
           CALL invfft ('Wave', aux, dffts)
           DO i = 1, dffts%nnr
              kkin(i) = kkin(i) + w1 * (dble(aux(i))**2 + aimag(aux(i))**2)
           ENDDO
           ! j
        ENDDO
        ! ibnd
     ENDDO
     ! ik
  ENDDO
#if defined(__MPI)
  !
  ! reduce local kinetic energy across pools
  !
  CALL mp_sum( kkin, inter_pool_comm )
#endif
  !
  ! interpolate the local kinetic energy to the dense grid
  ! Note that for US PP this term is incomplete: it contains
  ! only the contribution from the smooth part of the wavefunction
  !
  IF (doublegrid) THEN
     DEALLOCATE (aux)
     ALLOCATE(aux(dfftp%nnr))
     CALL fft_interpolate (dffts, kkin, dfftp, kkin)
  ENDIF
  !
  ! symmetrize the local kinetic energy if needed
  !
  IF ( .not. gamma_only) THEN
     !
     CALL sym_rho_init ( gamma_only )
     !
     aux(:) =  cmplx ( kkin (:), 0.0_dp, kind=dp)
     CALL fwfft ('Rho', aux, dfftp)
     ALLOCATE (aux2(ngm))
     aux2(:) = aux(dfftp%nl(:))
     !
     ! aux2 contains the local kinetic energy in G-space to be symmetrized
     !
     CALL sym_rho ( 1, aux2 )
     !
     aux(:) = (0.0_dp, 0.0_dp)
     aux(dfftp%nl(:)) = aux2(:)
     DEALLOCATE (aux2)
     CALL invfft ('Rho', aux, dfftp)
     kkin (:) = dble(aux(:))
     !
  ENDIF
  IF (iself) then
     !! calculation of the ELF !!
     !
     ! Calculate the bosonic kinetic density, stored in tbos
     !          aux --> charge density in Fourier space
     !         aux2 --> iG * rho(G)
     !
     ALLOCATE ( tbos(dfftp%nnr,nspin), aux2(dfftp%nnr) )
     !
     ! set rho%of_r(:,1) = spin up, rho%of_r(:,2) = spin down charge
     !
     IF (nspin == 2) CALL rhoz_or_updw( rho, 'r_and_g', '->updw' )
     !
     ! FIXME: the following gradient calculation could be performed
     !        in a more straightforward way using "fft_gradient_r2r"
     !
     DO k = 1, nspin
        tbos(:,k) = 0.d0
        aux(:) = cmplx( rho%of_r(:, k), 0.d0 ,kind=DP)
        CALL fwfft ('Rho', aux, dfftp)
        !
        DO j = 1, 3
           aux2(:) = (0.d0,0.d0)
           DO i = 1, ngm
              aux2(dfftp%nl(i)) = aux(dfftp%nl(i)) * cmplx(0.0d0, g(j,i)*tpiba,kind=DP)
           ENDDO
           IF (gamma_only) THEN
              DO i = 1, ngm
                 aux2(dfftp%nlm(i)) = aux(dfftp%nlm(i)) * cmplx(0.0d0,-g(j,i)*tpiba,kind=DP)
              ENDDO
           ENDIF

           CALL invfft ('Rho', aux2, dfftp)
           DO i = 1, dfftp%nnr
              tbos (i, k) = tbos (i, k) + dble(aux2(i))**2
           ENDDO
        ENDDO
     ENDDO
     !
     ! Calculate ELF
     !
     fac = 5.d0 / (3.d0 * (6.d0 * pi**2) ** (2.d0 / 3.d0) )
     fout(:) = 0.d0
     DO i = 1, dfftp%nnr
        IF ( nspin == 2 ) THEN
           IF ( (rho%of_r (i,1) > 1.d-30).and.(rho%of_r (i,2) > 1.d-30) ) THEN
              d = fac / ( rho%of_r(i,1)**(5d0/3d0) + rho%of_r(i,2)**(5d0/3d0) ) &
                 *(kkin(i) - 0.25d0*tbos(i,1)/rho%of_r(i,1) - &
                 0.25d0*tbos(i,2)/rho%of_r(i,2) + 1.d-5)
              fout (i) = 1.0d0 / (1.0d0 + d**2)
           ENDIF
        ELSE
           IF ( rho%of_r (i,1) > 1.d-30 ) THEN
              d = fac / ( rho%of_r(i,1)**(5d0/3d0) ) &
                 *(kkin(i) - 0.25d0*tbos(i,1)/rho%of_r(i,1) + 1.d-5)
              fout (i) = 1.0d0 / (1.0d0 + d**2)
           END IF
        ENDIF
     END DO
     DEALLOCATE (aux, aux2, tbos, kkin)
  ELSE
     !! calculation of the kinetic energy density !!
     fout = kkin
  ENDIF

END SUBROUTINE do_elf_kin
!-----------------------------------------------------------------------
SUBROUTINE do_rdg (rdg)
  !-----------------------------------------------------------------------
  !
  !  reduced density gradient
  !     rdg(r) = (1/2) (1/(3*pi**2))**(1/3) * |\nabla rho(r)|/rho(r)**(4/3)
  !
  USE kinds,                ONLY: DP
  USE constants,            ONLY: pi
  USE fft_base,             ONLY: dfftp
  USE scf,                  ONLY: rho
  USE gvect,                ONLY: g, ngm
  USE lsda_mod,             ONLY: nspin
  IMPLICIT NONE
  REAL(DP), INTENT(OUT) :: rdg (dfftp%nnr)
  REAL(DP), ALLOCATABLE :: grho(:,:)
  REAL(DP) :: fac
  REAL(DP), PARAMETER :: rho_cut = 0.05d0
  INTEGER :: is, i

  fac = (1.d0/2.d0) * 1.d0/(3.d0*pi**2)**(1.d0/3.d0)
  ! gradient of rho
  ALLOCATE( grho(3,dfftp%nnr) )

  CALL fft_gradient_g2r(dfftp, rho%of_g(1,1), g, grho)

  ! calculate rdg
  DO i = 1, dfftp%nnr
     IF (rho%of_r(i,1) > rho_cut) THEN
        rdg(i) = fac * 100.d0 / abs(rho%of_r(i,1))**(4.d0/3.d0)
     ELSE
        rdg(i) = fac * sqrt(grho(1,i)**2 + grho(2,i)**2 + grho(3,i)**2) &
                     / abs(rho%of_r(i,1))**(4.d0/3.d0)
     ENDIF
  ENDDO
  
  DEALLOCATE( grho )
  
  RETURN

END SUBROUTINE do_rdg

!-----------------------------------------------------------------------
SUBROUTINE do_sl2rho (sl2rho)
  !-----------------------------------------------------------------------
  !
  !  Computes sign(l2)*rho(r), where l2 is the second largest eigenvalue
  !  of the electron-density Hessian matrix
  !
  USE kinds,                ONLY: DP
  USE constants,            ONLY: pi
  USE fft_base,             ONLY: dfftp
  USE scf,                  ONLY: rho
  USE gvect,                ONLY: g, ngm
  USE lsda_mod,             ONLY: nspin
  IMPLICIT NONE
  REAL(DP), INTENT(OUT) :: sl2rho (dfftp%nnr)
  REAL(DP), ALLOCATABLE :: grho(:,:), hrho(:,:,:)
  INTEGER :: is, i
  !
  REAL(DP), PARAMETER :: eps = 1.d-14
  REAL(DP) :: vl, vu, work(24), hloc(3,3), e(3), v(3,3)
  INTEGER :: mo, info, iwork(15), ifail(3)
  !

  ! gradient and hessian of rho
  ALLOCATE( grho(3,dfftp%nnr), hrho(3,3,dfftp%nnr) )

  ! calculate hessian of rho (gradient is discarded)
  CALL fft_hessian( dfftp, rho%of_r(:,1), g, grho, hrho )

  ! find eigenvalues of the hessian
  DO i = 1, dfftp%nnr
     !
     IF (     abs(hrho(1,2,i) - hrho(2,1,i)) > eps &
         .OR. abs(hrho(1,3,i) - hrho(3,1,i)) > eps &
         .OR. abs(hrho(2,3,i) - hrho(3,2,i)) > eps   ) THEN
         CALL errore('do_sl2rho', 'hessian not symmetric', i)
     ENDIF
     !
     hloc = hrho(:,:,i)
     v (:,:) = 0.0_dp
     CALL DSYEVX ( 'V', 'I', 'U', 3, hloc, 3, vl, vu, 1, 3, 0.0_dp, mo, e,&
                  v, 3, work, 24, iwork, ifail, info )
     !
     IF ( info > 0) THEN
        CALL errore('do_sl2rho','failed to diagonalize',info)
     ELSEIF (info < 0) THEN
        call errore('do_sl2rho','illegal arguments in DSYEVX',-info)
     ENDIF
    
     sl2rho(i) = sign(1.d0,e(2))*rho%of_r(i,1)
  ENDDO


  DEALLOCATE( grho, hrho )
  
  RETURN

END SUBROUTINE do_sl2rho

!-----------------------------------------------------------------------
SUBROUTINE do_dori (dori)
  !-----------------------------------------------------------------------
  ! D. Yang & Q.Liu
  ! density overlap regions indicator（DOI: 10.1021/ct500490b）
  ! theta(r) = 4 * (laplacian(rho(r)) * grad(rho(r)) * rho(r) 
  !            + | grad(rho(r)) |**2 * grad(rho(r)))
  !            / (| grad(rho(r)) |**2)**3
  ! DORI(r) = theta(r) / (1 + theta(r))

  USE kinds,                ONLY: DP
  USE constants,            ONLY: pi
  USE fft_base,             ONLY: dfftp
  USE scf,                  ONLY: rho
  USE gvect,                ONLY: g, ngm
  USE lsda_mod,             ONLY: nspin
  IMPLICIT NONE
  REAL(DP), INTENT(OUT) :: dori (dfftp%nnr)
  REAL(DP), ALLOCATABLE :: grho(:,:), hrho(:,:,:), temp(:,:,:), sum_grho(:)
  INTEGER :: is, i, j

  ! gradient and hessian of rho
  ALLOCATE( grho(3,dfftp%nnr), hrho(3,3,dfftp%nnr))
  ALLOCATE(sum_grho(dfftp%nnr), temp(2,3,dfftp%nnr))


  ! calculate hessian of rho (gradient is discarded)
  CALL fft_hessian( dfftp, rho%of_r(:,1), g, grho, hrho )

  ! calculate theta(r)
  sum_grho(:) = grho(1,:)**2 + grho(2,:)**2 + grho(3,:)**2
  DO i = 1, 3
     temp(1,i,:) = 0.0d0
     DO j = 1,3
        temp(1,i,:) = temp(1,i,:) + grho(j,:)*hrho(i,j,:)
     ENDDO
     temp(1,i,:) = rho%of_r(:,1)*temp(1,i,:)
     temp(2,i,:) = grho(i,:)*sum_grho(:)
  ENDDO
  dori(:) = 0.d0 
  DO i = 1, 3
     dori(:) = dori(:) + (temp(1,i,:)-temp(2,i,:))**2
  ENDDO
  ! calculate dori(r) 
  dori(:) = 4/(sum_grho(:)+1.d-5)**3*dori(:)
  dori(:) = dori(:) / (1 + dori(:))
  
  DEALLOCATE( grho, hrho, temp )
  DEALLOCATE( sum_grho )

  RETURN

END SUBROUTINE do_dori

