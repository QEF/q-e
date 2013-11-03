!
! Copyright (C) 2001-2007 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------
subroutine addusforce (forcenl)
  !----------------------------------------------------------------------
  !
  !   This routine computes the contribution to atomic forces due
  !   to the dependence of the Q function on the atomic position.
  !   On output: the contribution is added to forcenl
  !
  USE kinds,      ONLY : DP
  USE ions_base,  ONLY : nat, ntyp => nsp, ityp
  USE cell_base,  ONLY : omega, tpiba
  USE fft_base,   ONLY : dfftp
  USE gvect,      ONLY : ngm, nl, nlm, gg, g, eigts1, eigts2, eigts3, mill
  USE noncollin_module,   ONLY : nspin_mag
  USE scf,        ONLY : v, vltot
  USE uspp,       ONLY : becsum, okvan
  USE uspp_param, ONLY : upf, lmaxq, nh, nhm
  USE mp_bands,   ONLY : intra_bgrp_comm
  USE mp,         ONLY : mp_sum
  USE control_flags, ONLY : gamma_only
  USE fft_interfaces,ONLY : fwfft
  !
  implicit none
  !
  real(DP) :: forcenl (3, nat)
  ! local variables
  integer :: ig, ir, dim, nt, ih, jh, ijh, ipol, is, na
  complex(DP):: cfac
  real(DP) :: fact, ddot
  ! work space
  complex(DP), allocatable :: aux(:,:), aux1(:,:), vg(:), qgm(:)
  real(DP) , allocatable :: ddeeq(:,:,:,:), qmod(:), ylmk0(:,:)

  !
  if (.not.okvan) return
  !
  if (gamma_only) then
     fact = 2.d0
  else
     fact = 1.d0
  end if
  allocate (aux(ngm,nspin_mag))    
  !
  ! fourier transform of the total effective potential
  !
  allocate (vg(dfftp%nnr))    
  do is = 1, nspin_mag
     if (nspin_mag.eq.4.and.is.ne.1) then
        vg (:) = v%of_r(:,is)
     else
        vg (:) = vltot (:) + v%of_r (:, is)
     endif
     CALL fwfft ('Dense', vg, dfftp)
     aux (:, is) = vg (nl (:) ) * tpiba * (0.d0, -1.d0)
  enddo
  deallocate (vg)
  !
  allocate (aux1(ngm,3))    
  allocate (ddeeq( 3, (nhm*(nhm+1))/2,nat,nspin_mag))    
  allocate (qgm( ngm))
  allocate (qmod( ngm))    
  allocate (ylmk0(ngm,lmaxq*lmaxq))    
  !
  ddeeq(:,:,:,:) = 0.d0
  !
  call ylmr2 (lmaxq * lmaxq, ngm, g, gg, ylmk0)
  !
  !qmod (:) = sqrt (gg (:) )
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ig)
  do ig = 1, ngm
     qmod (ig) = sqrt (gg (ig) )
  enddo
!$OMP END PARALLEL DO
  !
  ! here we compute the integral Q*V for each atom,
  !       I = sum_G i G_a exp(-iR.G) Q_nm v^*
  !
  do nt = 1, ntyp
     if ( upf(nt)%tvanp ) then
        ijh = 1
        do ih = 1, nh (nt)
           do jh = ih, nh (nt)
              call qvan2 (ngm, ih, jh, nt, qmod, qgm, ylmk0)
              do na = 1, nat
                 if (ityp (na) == nt) then
                    !
                    ! The product of potential, structure factor and iG
                    !
                    do is = 1, nspin_mag
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ig, cfac)
                       do ig = 1, ngm
                          cfac = aux (ig, is) * CONJG(eigts1 (mill(1,ig), na) *&
                                                      eigts2 (mill(2,ig), na) *&
                                                      eigts3 (mill(3,ig), na) )
                          aux1 (ig, 1) = g (1, ig) * cfac
                          aux1 (ig, 2) = g (2, ig) * cfac
                          aux1 (ig, 3) = g (3, ig) * cfac
                       enddo
!$OMP END PARALLEL DO
                       !
                       !    and the product with the Q functions
                       !    G=0 term gives no contribution
                       !
                       do ipol = 1, 3
                          ddeeq (ipol, ijh, na, is) = omega * fact * &
                               ddot (2 * ngm, aux1 (1, ipol), 1, qgm, 1)
                       enddo
                    enddo
                 endif
              enddo
              ijh = ijh + 1
           enddo
        enddo
     endif

  enddo
  call mp_sum ( ddeeq, intra_bgrp_comm )
  !            WRITE( stdout,'( "dmatrix atom ",i4)') na
  !            do ih = 1, nh(nt)
  !               WRITE( stdout,'(8f9.4)') (ddeeq(ipol,ih,jh,na),jh=1,nh(nt))
  !            end do
  !            WRITE( stdout,'( "dion pseudo ",i4)') nt
  !            do ih = 1, nh(nt)
  !               WRITE( stdout,'(8f9.4)') (dvan(ih,jh,nt),jh=1,nh(nt))
  !            end do
  do is = 1, nspin_mag
     do na = 1, nat
        nt = ityp (na)
        dim = (nh (nt) * (nh (nt) + 1) ) / 2
        do ipol = 1, 3
           do ir = 1, dim
              forcenl (ipol, na) = forcenl (ipol, na) + &
                   ddeeq (ipol, ir, na, is) * becsum (ir, na, is)
           enddo
        enddo
     enddo
  enddo
  deallocate (ylmk0)
  deallocate (qgm)
  deallocate (qmod)
  deallocate (ddeeq)
  deallocate (aux1)
  deallocate (aux)

  return
end subroutine addusforce

