!
! Copyright (C) 2001-2008 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------
subroutine compute_dvloc (uact, addnlcc, dvlocin)
  !----------------------------------------------------------------------
  !! This routine calculates \(dV_\text{bare}/\text{dtau}\) in real space
  !! for one perturbation with a given q. The displacements are described
  !! by a vector u. The result is stored in \(\text{dvlocin}\).
  !!
  !
  !
  USE kinds,            ONLY : DP
  USE funct,            ONLY : dft_is_nonlocc
  USE xc_lib,           ONLY : xclib_dft_is
  USE ions_base,        ONLY : nat, ityp
  USE cell_base,        ONLY : tpiba
  USE fft_base,         ONLY : dfftp, dffts
  USE fft_interfaces,   ONLY : fwfft, invfft
  USE gvect,            ONLY : eigts1, eigts2, eigts3, mill, g, ngm
  USE gvecs,            ONLY : ngms
  USE lsda_mod,         ONLY : nspin, lsda, current_spin
  USE scf,              ONLY : rho, rho_core
  USE noncollin_module, ONLY : nspin_gga
  use uspp_param,       ONLY : upf
  USE nlcc_ph,          ONLY : drc
  USE uspp,             ONLY : nlcc_any
  USE eqv,              ONLY : dmuxc, vlocq
  USE qpoint,           ONLY : xq, eigqts
  USE gc_lr,            ONLY : grho, dvxc_rr, dvxc_sr, dvxc_ss, dvxc_s
  USE modes,            ONLY : nmodes
  USE Coul_cut_2D,      ONLY : do_cutoff_2D
  USE Coul_cut_2D_ph,   ONLY : cutoff_localq
  !
  IMPLICIT NONE
  !
  COMPLEX(DP), INTENT(IN) :: uact(nmodes)
  !! input: the pattern of displacements
  LOGICAL, INTENT(IN) :: addnlcc
  !! input: if .true. the nonlinear core correction contribution is added
  COMPLEX(DP), INTENT(INOUT) :: dvlocin(dffts%nnr)
  !! output: the change of the local potential
  !
  ! ... local variables
  !
  INTEGER :: na
  !! counter on atoms
  INTEGER :: mu
  !! counter on modes
  INTEGER :: ig
  !! counter on G vectors
  INTEGER :: nt
  !! the type of atom
  INTEGER :: ir
  !! counter on real mesh
  !!
  INTEGER :: nnr, nnp, itmp, itmpp
  !!
  complex(DP) :: gtau, gu, fact, u1, u2, u3, gu0
  complex(DP), allocatable :: aux (:,:)
  complex(DP), pointer :: auxs (:)
  COMPLEX(DP), ALLOCATABLE :: drhoc(:,:)
  !
#if defined(__CUDA)
  INTEGER, POINTER, DEVICE :: nl_d(:), nlp_d(:)
  !
  nl_d  => dffts%nl_d
  nlp_d  => dfftp%nl_d
#else
  INTEGER, ALLOCATABLE :: nl_d(:)
  INTEGER, ALLOCATABLE :: nlp_d(:)
  !
  ALLOCATE( nl_d(dffts%ngm) )
  ALLOCATE( nlp_d(dfftp%ngm) )
  nl_d  = dffts%nl
  nlp_d  = dfftp%nl
#endif
  !
  call start_clock_gpu ('com_dvloc')
  !
  !    We start by computing the contribution of the local potential.
  !    The computation of the derivative of the local potential is done in
  !    reciprocal space while the product with the wavefunction is done in
  !    real space
  !
  nnr = dffts%nnr
  !
  !$acc data create(dvlocin(1:nnr)) copyin(vlocq,drc,dmuxc) deviceptr(nl_d, nlp_d)
  !$acc kernels present(dvlocin)
  dvlocin(:) = (0.d0, 0.d0)
  !$acc end kernels
  do na = 1, nat
     fact = tpiba * (0.d0, -1.d0) * eigqts (na)
     mu = 3 * (na - 1)
     u1 = uact (mu + 1)
     u2 = uact (mu + 2)
     u3 = uact (mu + 3)
     if ( abs(u1) + abs(u2) + abs(u3) > 1.0d-12) then
        nt = ityp (na)
        gu0 = xq (1) * u1 + xq (2) * u2 + xq (3) * u3
        !$acc parallel loop present(eigts1, eigts2, eigts3, mill, g, dvlocin)
        do ig = 1, ngms
           gtau = eigts1 (mill(1,ig), na) * eigts2 (mill(2,ig), na) * &
                  eigts3 (mill(3,ig), na)
           gu = gu0 + g (1, ig) * u1 + g (2, ig) * u2 + g (3, ig) * u3
           itmp = nl_d (ig)
           dvlocin (itmp) = dvlocin (itmp) + vlocq (ig, nt) * gu * fact * gtau
        enddo
        IF (do_cutoff_2D) then
           !$acc update host(dvlocin)
           call cutoff_localq( dvlocin, fact, u1, u2, u3, gu0, nt, na)
           !$acc update device(dvlocin)
        ENDIF
        !
     endif
  enddo
  !
  ! add NLCC when present
  !
  if (nlcc_any.and.addnlcc) then
     !CALL errore ('dvqpsi_us', 'openacc fpr ncll_any to be checked', 1)
     allocate (drhoc( dfftp%nnr,nspin))
     allocate (aux( dfftp%nnr,nspin))
     nnp=dfftp%nnr
     !$acc enter data create(drhoc(1:nnp,1:nspin),aux(1:nnp,1:nspin))
     !$acc kernels present(drhoc,aux)
     drhoc(:,:) = (0.d0, 0.d0)
     aux(:,:) = (0.0_dp, 0.0_dp)
     !$acc end kernels
     do na = 1,nat
        fact = tpiba*(0.d0,-1.d0)*eigqts(na)
        mu = 3*(na-1)
        u1 = uact(mu+1)
        u2 = uact(mu+2)
        u3 = uact(mu+3)
        if (abs(u1) + abs(u2) + abs(u3) > 1.0d-12) then
           nt=ityp(na)
           gu0 = xq(1)*u1 +xq(2)*u2+xq(3)*u3
           if (upf(nt)%nlcc) then
              !$acc parallel loop present(eigts1, eigts2, eigts3, g, mill,drhoc)
              do ig = 1,ngm
                 gtau = eigts1(mill(1,ig),na)*   &
                        eigts2(mill(2,ig),na)*   &
                        eigts3(mill(3,ig),na)
                 gu = gu0+g(1,ig)*u1+g(2,ig)*u2+g(3,ig)*u3
                 itmp = nlp_d(ig)
                 drhoc(itmp,1)=drhoc(itmp,1)+drc(ig,nt)*gu*fact*gtau
              enddo
           endif
        endif
     enddo
     !$acc host_data use_device(drhoc)
     CALL invfft ('Rho', drhoc(:,1), dfftp)
     !$acc end host_data
     if (.not.lsda) then
        !$acc parallel loop present(aux,drhoc)
        do ir=1,nnp
           aux(ir,1) = drhoc(ir,1) * dmuxc(ir,1,1)
        end do
     else
        !$acc parallel loop present(drhoc,aux) copyin(current_spin)
        do ir=1,nnp
           drhoc(ir,1) = 0.5d0 * drhoc(ir,1)
           drhoc(ir,2) = drhoc(ir,1)
           aux(ir,1) = drhoc(ir,1) * ( dmuxc(ir,current_spin,1) + &
                                       dmuxc(ir,current_spin,2) )
        enddo
     endif
     rho%of_r(:,1) = rho%of_r(:,1) + rho_core(:)
     !$acc exit data copyout(drhoc)

     IF ( xclib_dft_is('gradient') ) THEN
                    !$acc update host(aux)
                    CALL dgradcorr (dfftp, rho%of_r, grho, dvxc_rr, &
                    dvxc_sr, dvxc_ss, dvxc_s, xq, drhoc, nspin, nspin_gga, g, aux)
                    !$acc update device(aux)
     END IF
     IF (dft_is_nonlocc()) THEN
             !$acc update host(aux)               ! to fix double update
             CALL dnonloccorr(rho%of_r, drhoc, xq, aux)
             !$acc update device(aux)
     END IF
     deallocate (drhoc)

     rho%of_r(:,1) = rho%of_r(:,1) - rho_core(:)

     !$acc host_data use_device(aux)
     CALL fwfft ('Rho', aux(:,1), dfftp)
     !$acc end host_data
!
!  This is needed also when the smooth and the thick grids coincide to
!  cut the potential at the cut-off
!
     allocate (auxs(dffts%nnr))
     !$acc enter data create(auxs(1:nnr))
     !$acc kernels present(auxs)
     auxs(:) = (0.d0, 0.d0)
     !$acc end kernels
     !$acc parallel loop present(auxs,aux)
     do ig=1,ngms
        itmp = nl_d(ig)
        itmpp = nlp_d(ig)
        auxs(itmp) = aux(itmpp,1)
     enddo
     !$acc kernels present(dvlocin,auxs)
     dvlocin(:) = dvlocin(:) + auxs(:)
     !$acc end kernels
     !$acc exit data delete(aux, auxs)
     deallocate (aux)
     deallocate (auxs)
  endif
  !
  ! Now we compute dV_loc/dtau in real space
  !
  !$acc host_data use_device(dvlocin)
  CALL invfft ('Rho', dvlocin, dffts)
  !$acc end host_data
  !
#if !defined(__CUDA)
  DEALLOCATE(nl_d)
  DEALLOCATE(nlp_d)
#endif
  !
  call stop_clock_gpu ('com_dvloc')
  !
end subroutine compute_dvloc
