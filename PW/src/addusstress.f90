!
! Copyright (C) 2001-2012 Quantum ESPRESSO group
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
  real(DP) :: sigmanlc (3, 3)
  ! the nonlocal stress

  integer :: ig, nt, ih, jh, ijh, ipol, jpol, is, na
  ! counter on g vectors
  ! counter on mesh points
  ! number of composite nm components
  ! the atom type
  ! counter on atomic beta functions
  ! counter on atomic beta functions
  ! composite index for beta function
  ! counter on polarizations
  ! counter on polarizations
  ! counter on spin polarizations
  ! counter on atoms
  complex(DP), allocatable :: aux(:,:), aux1(:), vg(:), qgm(:)
  complex(DP)              :: cfac
  ! used to contain the potential
  ! used to compute a product
  ! used to contain the structure fac

  real(DP)               :: ps, ddot, sus(3,3)
  real(DP) , allocatable :: qmod(:), ylmk0(:,:), dylmk0(:,:)
  ! the integral
  ! the ultrasoft part of the stress
  ! the modulus of G
  ! the spherical harmonics
  ! the spherical harmonics derivativ
  !  of V_eff and dQ
  ! function which compute the scal.

  allocate ( aux(ngm,nspin), aux1(ngm), vg(dfftp%nnr), qgm(ngm), qmod(ngm) )
  allocate ( ylmk0(ngm,lmaxq*lmaxq), dylmk0(ngm,lmaxq*lmaxq) )

  !
  sus(:,:) = 0.d0
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
  do is = 1, nspin
     if ( nspin == 4 .and. is /= 1 ) then
        !
        vg(:) = v%of_r(:,is)
        !
     ELSE
        !
        vg(:) = vltot(:) + v%of_r(:,is)
        !
     END IF
     CALL fwfft ('Dense', vg, dfftp)
     do ig = 1, ngm
        aux (ig, is) = vg (nl (ig) )
     enddo
  enddo
  !
  ! here we compute the integral Q*V for each atom,
  !       I = sum_G i G_a exp(-iR.G) Q_nm v^*
  ! (no contribution from G=0)
  !
  do ipol = 1, 3
     call dylmr2 (lmaxq * lmaxq, ngm, g, gg, dylmk0, ipol)
     do nt = 1, ntyp
        if ( upf(nt)%tvanp ) then
           ijh = 1
           do ih = 1, nh (nt)
              do jh = ih, nh (nt)
                 call dqvan2 (ngm, ih, jh, nt, qmod, qgm, ylmk0, dylmk0, ipol)
                 do na = 1, nat
                    if (ityp (na) == nt) then
                       !
                       do is = 1, nspin
                          do jpol = 1, ipol
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ig, cfac)
                             do ig = 1, ngm
                                cfac = aux (ig, is) * &
                                       CONJG( eigts1 (mill (1,ig), na) * &
                                              eigts2 (mill (2,ig), na) * &
                                              eigts3 (mill (3,ig), na) )
                                aux1 (ig) = cfac * g (jpol, ig)
                             enddo
!$OMP END PARALLEL DO
                             !
                             !    and the product with the Q functions
                             !
                             ps = omega * ddot (2 * ngm, aux1, 1, qgm, 1)
                             sus (ipol, jpol) = sus (ipol, jpol) - &
                                                ps * becsum (ijh, na, is)
                          enddo
                       enddo
                    endif
                 enddo
                 ijh = ijh + 1
              enddo
           enddo
        endif
     enddo

  enddo

  if (gamma_only) then
     sigmanlc(:,:) = sigmanlc(:,:) + 2.d0*sus(:,:)
  else
     sigmanlc(:,:) = sigmanlc(:,:) + sus(:,:)
  end if
  deallocate (ylmk0, dylmk0)
  deallocate (aux, aux1, vg, qgm, qmod)

  return

end subroutine addusstres

