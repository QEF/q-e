!
! Copyright (C) 2001-2015 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------
subroutine addusstres (sigmanlc)
  !----------------------------------------------------------------------
  !
  !   This routine computes the part of the atomic force which is due
  !   to the dependence of the Q function on the atomic position.
  !   Adds contribution to input sigmanlc, does not sum contributions
  !   from various processors (sum is performed by calling routine)
  !
  USE kinds,      ONLY : DP
  USE ions_base,  ONLY : nat, ntyp => nsp, ityp
  USE cell_base,  ONLY : omega, tpiba
  USE fft_base,   ONLY : dfftp
  USE gvect,      ONLY : ngm, nl, nlm, gg, g, eigts1, eigts2, eigts3, mill
  USE lsda_mod,   ONLY : nspin
  USE scf,        ONLY : v, vltot
  USE uspp,       ONLY : becsum, okvan
  USE uspp_param, ONLY : upf, lmaxq, nh, nhm
  USE control_flags, ONLY : gamma_only
  USE fft_interfaces,ONLY : fwfft
  !
  implicit none
  !
  real(DP), INTENT(INOUT) :: sigmanlc (3, 3)
  ! the nonlocal stress
  integer :: ig, nt, ih, jh, ijh, ipol, jpol, is, na, nij
  ! counters
  complex(DP), allocatable :: aux(:), aux1(:), aux2(:,:), vg(:,:), qgm(:,:)
  ! work space (complex)
  complex(DP)              :: cfac
  real(DP)               :: ps, ddot, sus(3,3)
  ! auxiliary variables
  real(DP) , allocatable :: qmod(:), ylmk0(:,:), dylmk0(:,:), tbecsum(:,:)
  ! work space (real)
  !
  !
  sus(:,:) = 0.d0
  !
  allocate ( aux1(ngm), aux2(ngm,nspin), qmod(ngm) )
  allocate ( ylmk0(ngm,lmaxq*lmaxq), dylmk0(ngm,lmaxq*lmaxq) )
  !
  call ylmr2 (lmaxq * lmaxq, ngm, g, gg, ylmk0)
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ig)
  do ig = 1, ngm
     qmod (ig) = sqrt (gg (ig) )
  enddo
!$OMP END PARALLEL DO
  !
  ! fourier transform of the total effective potential
  !
  allocate ( vg(ngm,nspin))
  allocate ( aux(dfftp%nnr) )
  do is = 1, nspin
     if ( nspin == 4 .and. is /= 1 ) then
        aux(:) = v%of_r(:,is)
     ELSE
        aux(:) = vltot(:) + v%of_r(:,is)
     END IF
     CALL fwfft ('Dense', aux, dfftp)
     do ig = 1, ngm
        vg (ig, is) = aux (nl (ig) )
     enddo
  enddo
  deallocate ( aux )
  !
  ! here we compute the integral Q*V for each atom,
  !       I = sum_G i G_a exp(-iR.G) Q_nm v^*
  ! (no contribution from G=0)
  !
  do ipol = 1, 3
     call dylmr2 (lmaxq * lmaxq, ngm, g, gg, dylmk0, ipol)
     do nt = 1, ntyp
        if ( upf(nt)%tvanp ) then
           nij = nh(nt)*(nh(nt)+1)/2
           allocate (qgm(ngm,nij), tbecsum(nij,nspin) )
           ijh = 0
           do ih = 1, nh (nt)
              do jh = ih, nh (nt)
                 ijh = ijh + 1
                 call dqvan2 (ngm, ih, jh, nt, qmod, qgm(1,ijh), ylmk0, &
                      dylmk0, ipol)
              end do
           end do
           !
           do na = 1, nat
              if (ityp (na) == nt) then
                 !
                 tbecsum(:,:) = becsum(1:nij,na,1:nspin)
                 !
                 CALL dgemm( 'N', 'N', 2*ngm, nspin, nij, 1.0_dp, &
                      qgm, 2*ngm, tbecsum, nij, 0.0_dp, aux2, 2*ngm )
                 do is = 1, nspin
                    do jpol = 1, ipol
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ig, cfac)
                       do ig = 1, ngm
                          cfac = vg (ig, is) * &
                               CONJG( eigts1 (mill (1,ig), na) * &
                               eigts2 (mill (2,ig), na) * &
                               eigts3 (mill (3,ig), na) )
                          aux1 (ig) = cfac * g (jpol, ig)
                       enddo
!$OMP END PARALLEL DO
                       !
                       !    and the product with the Q functions
                       !
                       ps = omega * ddot (2 * ngm, aux1, 1, aux2(1,is), 1)
                       sus (ipol, jpol) = sus (ipol, jpol) - ps
                    enddo
                 enddo
              endif
           enddo
           deallocate ( tbecsum, qgm )
        endif
     enddo

  enddo

  if (gamma_only) then
     sigmanlc(:,:) = sigmanlc(:,:) + 2.d0*sus(:,:)
  else
     sigmanlc(:,:) = sigmanlc(:,:) + sus(:,:)
  end if
  deallocate (ylmk0, dylmk0)
  deallocate (aux1, aux2, vg, qmod)

  return

end subroutine addusstres

