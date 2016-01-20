!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------
subroutine addusddens (drhoscf, dbecsum, mode0, npe, iflag)
  !----------------------------------------------------------------------
  !
  !  This routine adds to the change of the charge and of the
  !  magnetization densities the part due to the US augmentation.
  !  It assumes that the array dbecsum has already accumulated the
  !  change of the becsum term. It calculates Eq. B31 of Ref [1].
  !  If called from drho (iflag=1), dbecsum and drhoscf contain the
  !  orthogonalization contribution to the change of the wavefunctions
  !  and the terms with alphasum and becsum are added. If called
  !  from solve_* (iflag=0) drhoscf and dbecsum contain the contribution
  !  of the solution of the linear system and the terms due to alphasum
  !  and becsum are not added. In this case the change of the charge
  !  calculated by drho (called \Delta \rho in [1]) is read from file
  !  and added. The contribution of the change of
  !  the Fermi energy is not calculated here but added later by ef_shift.
  !  [1] PRB 64, 235118 (2001).
  !
  !
  USE kinds, only : DP
  use fft_base,  only: dfftp
  use fft_interfaces, only: invfft
  USE gvect,  ONLY : gg, ngm, nl, g, eigts1, eigts2, eigts3, mill
  USE uspp,     ONLY : okvan, becsum
  USE cell_base, ONLY : tpiba
  USE ions_base, ONLY : nat, ityp, ntyp => nsp
  USE wavefunctions_module,  ONLY: psic
  USE buffers,    ONLY : get_buffer
  USE uspp_param, ONLY: upf, lmaxq, nh, nhm
  USE paw_variables, ONLY : okpaw
  USE modes,     ONLY : u
  USE phus,    ONLY : becsumort, alphasum
  USE units_ph,  ONLY : iudrhous, lrdrhous
  USE noncollin_module, ONLY : nspin_mag

  USE qpoint,     ONLY : xq, eigqts
  USE control_lr, ONLY : lgamma

  implicit none
  !
  !   the dummy variables
  !

  integer :: iflag, npe
  ! input: if zero does not compute drho
  ! input: the number of perturbations

  complex(DP) :: drhoscf (dfftp%nnr, nspin_mag, npe), &
                      dbecsum (nhm*(nhm+1)/2, nat, nspin_mag, npe)
  ! inp/out: change of the charge density
  !input: sum over kv of bec
  integer ::  mode0
  ! input:the mode of the representation
  !
  !     here the local variables
  !

  integer :: ig, na, nt, ih, jh, mu, mode, ipert, is, ijh
  ! counter on G vectors
  ! counter on atoms
  ! counter on atomic type
  ! counter on beta functions
  ! counter on beta functions
  ! counter on r vectors
  ! pointer on modes
  ! pointer on the mode
  ! counter on perturbations
  ! counter on spin
  ! counter on combined beta functions

  real(DP), allocatable  :: qmod (:), qpg (:,:), ylmk0 (:,:)
  ! the modulus of q+G
  ! the values of q+G
  ! the spherical harmonics

  complex(DP) :: fact, zsum, bb, alpha, alpha_0, u1, u2, u3
  ! auxiliary variables
  complex(DP), allocatable ::  sk (:), qgm(:), drhous (:,:), aux (:,:,:)
  ! the structure factor
  ! q_lm(G)
  ! contain the charge of drho
  ! auxiliary variable for drho(G)

  if (.not.okvan) return
  call start_clock ('addusddens')
  allocate (aux(  ngm , nspin_mag , npe))
  allocate (sk (  ngm))
  allocate (ylmk0(ngm , lmaxq * lmaxq))
  allocate (qgm(  ngm))
  allocate (qmod( ngm))
  if (.not.lgamma) allocate (qpg( 3  , ngm))
  !      WRITE( stdout,*) aux, ylmk0, qmod
  !
  !  And then we compute the additional charge in reciprocal space
  !
  if (.not.lgamma) then
     call setqmod (ngm, xq, g, qmod, qpg)
     call ylmr2 (lmaxq * lmaxq, ngm, qpg, qmod, ylmk0)
     do ig = 1, ngm
        qmod (ig) = sqrt (qmod (ig) )
     enddo
  else
     call ylmr2 (lmaxq * lmaxq, ngm, g, gg, ylmk0)
     do ig = 1, ngm
        qmod (ig) = sqrt (gg (ig) )
     enddo
  endif
  fact = cmplx (0.d0, - tpiba, kind=DP)
  aux(:,:,:) = (0.d0, 0.d0)
  do nt = 1, ntyp
     if (upf(nt)%tvanp  ) then
        ijh = 0
        do ih = 1, nh (nt)
           do jh = ih, nh (nt)
              call qvan2 (ngm, ih, jh, nt, qmod, qgm, ylmk0)
              ijh = ijh + 1
              do na = 1, nat
                 if (ityp (na) .eq.nt) then
                    mu = 3 * (na - 1)
                    !
                    ! calculate the structure factor
                    !
                    do ig = 1, ngm
                       sk (ig) = eigts1 (mill(1,ig), na) * &
                                 eigts2 (mill(2,ig), na) * &
                                 eigts3 (mill(3,ig), na) * &
                                 eigqts (na) * qgm (ig)
                    enddo
                    !
                    !  And qgmq and becp and dbecq
                    !
                    do ipert = 1, npe
                       do is = 1, nspin_mag
                          mode = mode0 + ipert
                          if (iflag==1) then
                             zsum = dbecsum (ijh, na, is, ipert)
                          else
                             zsum = 2.0_DP*dbecsum (ijh, na, is, ipert)
                          endif
                          u1 = u (mu + 1, mode)
                          u2 = u (mu + 2, mode)
                          u3 = u (mu + 3, mode)
                          if (abs(u1) + abs(u2) + abs(u3) .gt.1d-12 .and. &
                              iflag.eq.1) then
                             bb = becsum (ijh, na, is)
                             zsum = zsum + &
                                  ( alphasum (ijh, 1, na, is) * u1 &
                                  + alphasum (ijh, 2, na, is) * u2 &
                                  + alphasum (ijh, 3, na, is) * u3)
                             IF (okpaw) becsumort(ijh,na,is,mode) =  zsum
                             u1 = u1 * fact
                             u2 = u2 * fact
                             u3 = u3 * fact
                             alpha_0 = xq(1)*u1 + xq(2)*u2 + xq(3)*u3
                             do ig = 1, ngm
                                alpha = alpha_0 + &
                                        g(1,ig)*u1 + g(2,ig)*u2 + g(3,ig)*u3
                                aux(ig,is,ipert) = aux(ig,is,ipert) + &
                                                   (zsum + alpha*bb) * sk(ig)
                             enddo
                          else
                             call zaxpy (ngm, zsum, sk, 1, aux(1,is,ipert), 1)
                             IF (okpaw.and.iflag==1) &
                                    becsumort(ijh,na,is,mode) = zsum
                          endif
                       enddo
                    enddo
                 endif
              enddo
           enddo
        enddo
     endif
  enddo
  !
  !     convert aux to real space
  !
  do ipert = 1, npe
     mu = mode0 + ipert
     do is = 1, nspin_mag
        psic(:) = (0.d0, 0.d0)
        do ig = 1, ngm
           psic (nl (ig) ) = aux (ig, is, ipert)
        enddo
        CALL invfft ('Dense', psic, dfftp)
        call daxpy (2*dfftp%nnr, 1.0_DP, psic, 1, drhoscf(1,is,ipert), 1)
     enddo
  enddo
  if (.not.lgamma) deallocate (qpg)
  deallocate (qmod)
  deallocate (qgm)
  deallocate (ylmk0)
  deallocate (sk)
  deallocate (aux)

  if (iflag == 0) then
     allocate (drhous( dfftp%nnr, nspin_mag))
     do ipert = 1, npe
        mu = mode0 + ipert
        call get_buffer (drhous, lrdrhous, iudrhous, mu)
        call daxpy (2*dfftp%nnr*nspin_mag, 1.d0, drhous, 1, drhoscf(1,1,ipert), 1)
     end do
     deallocate (drhous)
  end if

  call stop_clock ('addusddens')
  return
end subroutine addusddens
