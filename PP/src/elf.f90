   !
! Copyright (C) 2001-2005 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
SUBROUTINE do_elf (elf)
  !-----------------------------------------------------------------------
  !
  !  calculation of the electron localization function;
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
  USE cell_base, ONLY: omega, tpiba, tpiba2
  USE fft_base,  ONLY: dffts, dfftp
  USE fft_interfaces, ONLY : fwfft, invfft
  USE gvect, ONLY: gcutm, g, ngm, nl, nlm
  USE gvecs, ONLY : nls, nlsm, ngms, doublegrid, dual
  USE io_files, ONLY: iunwfc, nwordwfc
  USE klist, ONLY: nks, xk
  USE lsda_mod, ONLY: nspin
  USE scf, ONLY: rho
  USE symme, ONLY: sym_rho, sym_rho_init
  USE wvfct, ONLY: npw, igk, g2kin, nbnd, wg, ecutwfc
  USE control_flags, ONLY: gamma_only
  USE wavefunctions_module,  ONLY: evc
  USE mp_global,            ONLY: inter_pool_comm, intra_pool_comm
  USE mp,                   ONLY: mp_sum
  !
  ! I/O variables
  !
  IMPLICIT NONE
  real(DP) :: elf (dfftp%nnr)
  !
  ! local variables
  !
  INTEGER :: i, j, k, ibnd, ik, is
  real(DP) :: gv(3), w1, d, fac
  real(DP), ALLOCATABLE :: kkin (:), tbos (:)
  COMPLEX(DP), ALLOCATABLE :: aux (:), aux2 (:)
  !
  CALL infomsg ('do_elf', 'elf + US not fully implemented')
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
     !    prepare the indices of this k point
     !
     CALL gk_sort (xk (1, ik), ngm, g, ecutwfc / tpiba2, npw, igk, g2kin)
     !
     !   reads the eigenfunctions
     !
     CALL davcio (evc, nwordwfc, iunwfc, ik, - 1)
     !
     DO ibnd = 1, nbnd
        DO j = 1, 3
           aux(:) = (0.d0,0.d0)
           w1 = wg (ibnd, ik) / omega
           DO i = 1, npw
              gv (j) = (xk (j, ik) + g (j, igk (i) ) ) * tpiba
              aux (nls(igk (i) ) ) = cmplx(0d0, gv (j) ,kind=DP) * evc (i, ibnd)
              IF (gamma_only) THEN
                 aux (nlsm(igk (i) ) ) = cmplx(0d0, -gv (j) ,kind=DP) * &
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
#ifdef __MPI
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
     CALL interpolate (kkin, kkin, 1)
  ENDIF
  !
  ! symmetrize the local kinetic energy if needed
  !
  IF ( .not. gamma_only) THEN
     !
     CALL sym_rho_init ( gamma_only )
     !
     aux(:) =  cmplx ( kkin (:), 0.0_dp, kind=dp)
     CALL fwfft ('Smooth', aux, dffts)
     ALLOCATE (aux2(ngm))
     aux2(:) = aux(nl(:))
     !
     ! aux2 contains the local kinetic energy in G-space to be symmetrized
     !
     CALL sym_rho ( 1, aux2 )
     !
     aux(:) = (0.0_dp, 0.0_dp)
     aux(nl(:)) = aux2(:)
     DEALLOCATE (aux2)
     CALL invfft ('Dense', aux, dfftp)
     kkin (:) = dble(aux(:))
     !
  ENDIF
  !
  ! Calculate the bosonic kinetic density, stored in tbos
  !          aux --> charge density in Fourier space
  !         aux2 --> iG * rho(G)
  !
  ALLOCATE ( tbos(dfftp%nnr), aux2(dfftp%nnr) )
  tbos(:) = 0.d0
  !
  ! put the total (up+down) charge density in rho%of_r(*,1)
  !
  DO is = 2, nspin
     rho%of_r (:, 1) =  rho%of_r (:, 1) + rho%of_r (:, is)
  ENDDO
  !
  aux(:) = cmplx( rho%of_r(:, 1), 0.d0 ,kind=DP)
  CALL fwfft ('Dense', aux, dfftp)
  !
  DO j = 1, 3
     aux2(:) = (0.d0,0.d0)
     DO i = 1, ngm
        aux2(nl(i)) = aux(nl(i)) * cmplx(0.0d0, g(j,i)*tpiba,kind=DP)
     ENDDO
     IF (gamma_only) THEN
        DO i = 1, ngm
           aux2(nlm(i)) = aux(nlm(i)) * cmplx(0.0d0,-g(j,i)*tpiba,kind=DP)
        ENDDO
     ENDIF

     CALL invfft ('Dense', aux2, dffts)
     DO i = 1, dfftp%nnr
        tbos (i) = tbos (i) + dble(aux2(i))**2
     ENDDO
  ENDDO
  !
  ! Calculates ELF
  !
  fac = 5.d0 / (3.d0 * (3.d0 * pi**2) ** (2.d0 / 3.d0) )
  elf(:) = 0.d0
  DO i = 1, dfftp%nnr
     IF (rho%of_r (i,1) > 1.d-30) THEN
        d = fac / rho%of_r(i,1)**(5d0/3d0) * (kkin(i)-0.25d0*tbos(i)/rho%of_r(i,1))
        elf (i) = 1.0d0 / (1.0d0 + d**2)
     ENDIF
  ENDDO
  DEALLOCATE (aux, aux2, tbos, kkin)
  RETURN
END SUBROUTINE do_elf



!-----------------------------------------------------------------------
SUBROUTINE do_rdg (rdg)
  !-----------------------------------------------------------------------
  !
  !  reduced density gradient
  !     rdg(r) = (1/2) (1/(3*pi**2))**(1/3) * |\nabla rho(r)|/rho(r)**(4/3)
  !
  USE kinds,                ONLY: DP
  USE constants,            ONLY: pi
  USE cell_base,            ONLY: omega, tpiba, tpiba2
  USE fft_base,             ONLY: dfftp
  USE scf,                  ONLY: rho
  USE gvect,                ONLY: g, ngm, nl
  USE lsda_mod,             ONLY: nspin
  IMPLICIT NONE
  real(dp), intent(out) :: rdg (dfftp%nnr)
  real(dp), allocatable :: grho(:,:)
  real(dp):: fac
  real(dp), parameter :: rho_cut = 0.05d0
  integer :: is, i

  fac = (1.d0/2.d0) * 1.d0/(3.d0*pi**2)**(1.d0/3.d0)
  ! gradient of rho
  allocate( grho(3,dfftp%nnr) )

  ! put the total (up+down) charge density in rho%of_r(*,1)
  do is = 2, nspin
     rho%of_g(:,1) =  rho%of_g(:,1) + rho%of_g(:,is)
     rho%of_r(:,1) =  rho%of_r(:,1) + rho%of_r(:,is)
  enddo

  ! gradient of rho
  call gradrho(dfftp%nnr, rho%of_g(1,1), ngm, g, nl, grho)

  ! calculate rdg
  do i = 1, dfftp%nnr
    if (rho%of_r(i,1) > rho_cut) then
      rdg(i) = fac * 100.d0 / abs(rho%of_r(i,1))**(4.d0/3.d0)
    else
      rdg(i) = fac * sqrt(grho(1,i)**2 + grho(2,i)**2 + grho(3,i)**2) / abs(rho%of_r(i,1))**(4.d0/3.d0)
    endif
  enddo

  deallocate( grho )
  
  return

END SUBROUTINE do_rdg
