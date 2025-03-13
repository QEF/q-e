!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------
subroutine dvanqq
  !----------------------------------------------------------------------
  !! This routine calculates four integrals of the Q functions and
  !! its derivatives with \(\text{Vloc}\) and \(\text{Veff}\) which are used
  !! to compute term \(dV_\text{bare}/d\tau \cdot \psi\) in 
  !! \(\texttt{addusdvqpsi}\) and in \(\texttt{addusdynmat}\).
  !! The result is stored in int1, int2, int4, int5. The routine is called
  !! only once. int4 and int5 are deallocated after use in \(\texttt{addusdynmat}\).
  !
  !! int1: Eq.(B20) of Ref.[1];  
  !! int2: Eq.(B21) of Ref.[1];  
  !! int4: Eq.(B23) of Ref.[1];  
  !! int5: Eq.(B24) of Ref.[1].
  !
  !! [1] PRB 64, 235118 (2001).
  !
  !
  USE kinds, only : DP
  USE cell_base, ONLY : omega, tpiba2, tpiba
  USE ions_base, ONLY : nat, ityp, ntyp => nsp
  USE fft_base,   ONLY: dfftp
  USE fft_interfaces, ONLY: fwfft
  USE gvect, ONLY : ngm, gg, g, mill, eigts1, eigts2, eigts3
  USE scf, ONLY : v, vltot
  USE noncollin_module, ONLY : noncolin, nspin_mag, lspinorb
  USE uspp, ONLY: okvan, ijtoh
  USE uspp_param, ONLY: upf, lmaxq, nh

  USE mp_bands, ONLY: intra_bgrp_comm
  USE mp,        ONLY: mp_sum

  USE phus, ONLY : int1, int2, int4, int4_nc, int5, int5_so
  USE partial, ONLY :  nat_todo, nat_todo_input, atomo, set_local_atomo
  USE lr_symm_base, ONLY  :  nsymq 
  USE symm_base,    ONLY  :  irt    

  USE eqv,        ONLY : vlocq
  USE qpoint,     ONLY : eigqts, xq
  USE control_lr, ONLY : lgamma, rec_code_read

  implicit none
  !
  !   And the local variables
  !
  integer  ::  nat_l
  ! effective number of atoms to do when nat_todo_input > 0 
  integer,allocatable :: atomo_l(:)
  integer :: nt, na, nb, ig, nta, ntb, ir, ih, jh, ijh, ipol, jpol, is, na_l, nb_l 
  ! counters

  real(DP), allocatable :: qmod (:), qmodg (:), qpg (:,:), &
       ylmkq (:,:), ylmk0 (:,:)
  ! the modulus of q+G
  ! the modulus of G
  ! the  q+G vectors
  ! the spherical harmonics

  complex(DP) :: fact, fact1,z9aux(9)
  complex(DP), allocatable :: aux1 (:), aux2 (:),&
       aux3 (:), aux5 (:), aux35(:,:), veff (:,:), sk(:)
  ! work space
  complex(DP), allocatable, target :: qgm(:)
  ! the augmentation function at G
  complex(DP), pointer :: qgmq (:)
  ! the augmentation function at q+G

  if (.not.okvan) return

  if (rec_code_read >= -20) return

  call start_clock ('dvanqq')
  int1(:,:,:,:,:) = (0.d0, 0.d0)
  int2(:,:,:,:,:) = (0.d0, 0.d0)
  int4(:,:,:,:,:) = (0.d0, 0.d0)
  int5(:,:,:,:,:) = (0.d0, 0.d0)
  allocate (sk  (  ngm))
  allocate (aux1(  ngm))
  allocate(aux35(9,ngm))
  allocate (qmodg( ngm))
  allocate (ylmk0( ngm , lmaxq * lmaxq))
  allocate (qgm  ( ngm))
  if (.not.lgamma) then
     allocate (ylmkq(ngm , lmaxq * lmaxq))
     allocate (qmod( ngm))
     allocate (qgmq( ngm))
  else
     qgmq =>qgm
  endif
  !
  !     compute spherical harmonics
  !
  call ylmr2 (lmaxq * lmaxq, ngm, g, gg, ylmk0)
  do ig = 1, ngm
     qmodg (ig) = sqrt (gg (ig) ) * tpiba
  enddo
  if (.not.lgamma) then
     allocate (qpg (3, ngm))
     call setqmod (ngm, xq, g, qmod, qpg)
     call ylmr2 (lmaxq * lmaxq, ngm, qpg, qmod, ylmkq)
     deallocate (qpg)
     do ig = 1, ngm
        qmod (ig) = sqrt (qmod (ig) ) * tpiba
     enddo
  endif
  if (nat_todo_input > 0 ) then 
      call set_local_atomo(nat, nat_todo, atomo, nsymq, irt, nat_l, atomo_l)
  else 
     nat_l = nat 
  end if 
  !
  !   we start by computing the FT of the effective potential
  !
  allocate (veff ( dfftp%nnr , nspin_mag))
  do is = 1, nspin_mag
     if (nspin_mag.ne.4.or.is==1) then
        do ir = 1, dfftp%nnr
           veff (ir, is) = CMPLX(vltot (ir) + v%of_r (ir, is), 0.d0,kind=DP)
        enddo
     else
        do ir = 1, dfftp%nnr
           veff (ir, is) = CMPLX(v%of_r (ir, is), 0.d0,kind=DP)
        enddo
     endif
     CALL fwfft ('Rho', veff (:, is), dfftp)
  enddo
  !
  !     We compute here four of the five integrals needed in the phonon
  !
  fact1 = CMPLX(0.d0, - tpiba * omega,kind=DP)
  !
  do ntb = 1, ntyp
     if (upf(ntb)%tvanp ) then
        do ih = 1, nh (ntb)
           do jh = ih, nh (ntb)
              ijh = ijtoh(ih,jh,ntb) 
              !
              !    compute the augmentation function
              !
              call qvan2 (ngm, ih, jh, ntb, qmodg, qgm, ylmk0)
              !
              if (.not.lgamma) call qvan2 (ngm, ih, jh, ntb, qmod, qgmq, ylmkq)
              !
              !     NB: for this integral the moving atom and the atom of Q
              !     do not necessarily coincide
              !
              do nb = 1, nat
                 if (ityp (nb) == ntb) then
                    do ig = 1, ngm
                       aux1 (ig) = qgmq (ig) * eigts1 (mill(1,ig), nb) &
                                             * eigts2 (mill(2,ig), nb) &
                                             * eigts3 (mill(3,ig), nb)
                    enddo
                    do na_l = 1, nat_l
                       if (nat_l < nat ) then 
                          na = atomo_l(na_l) 
                       else 
                          na = na_l
                       end if 
                       fact = eigqts (na) * CONJG(eigqts (nb) )
                       !
                       !    nb is the atom of the augmentation function
                       !
                       nta = ityp (na)
                       !  
                       do ig=1, ngm
                          sk(ig)=  vlocq(ig,nta) * eigts1(mill(1,ig), na) &
                                                 * eigts2(mill(2,ig), na) &
                                                 * eigts3(mill(3,ig), na)
                       enddo
                       ! 
                       ! FIXME: replace zgemv with zgemm
                       !
                       !$omp parallel do default(shared) private(ig)
                       do ig =1, ngm 
                          aux35(1:3,ig) = conjg(sk(ig)) * (g(1:3,ig) + xq (1:3))
                          aux35(4:6,ig) = aux35(1,ig) *  (g(1:3,ig) + xq (1:3))
                          aux35(7:8,ig) = aux35(2,ig) *  (g(2:3,ig) + xq (2:3))
                          aux35(9,ig)  =  aux35(3,ig) *  (g(3,ig) + xq (3))
                       end do 
                       call zgemv('N', 9,ngm,cmplx(1._dp, 0._dp,kind=dp),aux35,9,aux1,1,(0._dp,0._dp),z9aux,1)
                       z9aux(4:9) = conjg(fact)*tpiba2*omega*z9aux(4:9)
                       !
                       int2(ih,jh,1:3,na,nb) = conjg(z9aux(1:3)) * fact * fact1
                       ! 
                       z9aux(1:3) = z9aux(4:6)
                       z9aux( 4)  = z9aux (2)
                       z9aux(5:6) = z9aux(7:8)
                       z9aux(7)   = z9aux(3)
                       z9aux(8)   = z9aux(6)
                       !
                       call zcopy(9, z9aux(1),1, int5(ijh,1,1,na,nb), size(int5,1)) 
                       !
                    enddo
                    if (.not.lgamma) then
                       do ig = 1, ngm
                          aux1 (ig) = qgm (ig) * eigts1 (mill(1,ig), nb) &
                                               * eigts2 (mill(2,ig), nb) &
                                               * eigts3 (mill(3,ig), nb)
                       enddo
                    endif
                    do is = 1, nspin_mag
                       !$omp parallel do default(shared) private(ig)
                       do ig = 1, ngm 
                          aux35(1:3,ig) = conjg(veff (dfftp%nl (ig), is)) * g (1:3, ig)
                          aux35(4:6,ig) = aux35(1,ig) * g(1:3,ig)
                          aux35(7:8,ig) = aux35(2,ig) * g(2:3,ig)
                          aux35(9,ig)   = aux35(3,ig) * g(3,ig) 
                       end do 
                       call zgemv('N',9,ngm,cmplx(1._dp, 0._dp,kind=dp), aux35,9,aux1,1,(0._dp,0._dp),z9aux,1)
                       !
                       z9aux(4:9) = -tpiba2 * omega * z9aux(4:9)
                       int1(ih,jh,1:3,nb,is) = -fact1  * conjg(z9aux(1:3))
                       ! 
                       z9aux(1:3) = z9aux(4:6)
                       z9aux( 4)  = z9aux (2)
                       z9aux(5:6) = z9aux(7:8)
                       z9aux(7)   = z9aux(3)
                       z9aux(8)   = z9aux(6)
                       !
                       call zcopy(9, z9aux(1),1, int4(ijh,1,1,nb,is), size(int4,1)) 
                       !
                    enddo
                 endif
              enddo
           enddo
        enddo
        do ih = 1, nh (ntb)
           do jh = ih + 1, nh (ntb)
              !
              !    We use the symmetry properties of the integral factor
              !
              do nb =1, nat 
                 if (ityp (nb) == ntb) then
                    do ipol = 1, 3
                       do is = 1, nspin_mag
                          int1(jh,ih,ipol,nb,is) = int1(ih,jh,ipol,nb,is)
                       enddo
                       do na_l = 1, nat_l
                          if (nat_l < nat ) then 
                             na = atomo_l(na_l)
                          else 
                             na = na_l
                          end if 
                          int2(jh,ih,ipol,na,nb) = int2(ih,jh,ipol,na,nb)
                       enddo
                    enddo
                 endif
              enddo
           enddo
        enddo
     endif
  enddo

  call mp_sum(  int1, intra_bgrp_comm )
  call mp_sum(  int2, intra_bgrp_comm )
  call mp_sum(  int4, intra_bgrp_comm )
  call mp_sum(  int5, intra_bgrp_comm )

  IF (noncolin) THEN
     CALL set_int12_nc(0)
     int4_nc = (0.d0, 0.d0)
     IF (lspinorb) int5_so = (0.d0, 0.d0)
     DO nt = 1, ntyp
        IF ( upf(nt)%tvanp ) THEN
           DO na_l = 1, nat_l
              IF (nat_l < nat)  THEN 
                na = atomo_l(na_l)  
              ELSE 
                na = na_l 
              END IF 
              IF (ityp(na)==nt) THEN
                 IF (upf(nt)%has_so) THEN
                    CALL transform_int4_so(int4,na)
                    CALL transform_int5_so(int5,na)
                 ELSE
                    CALL transform_int4_nc(int4,na)
                    IF (lspinorb) CALL transform_int5_nc(int5,na)
                 END IF
              END IF
           END DO
        END IF
     END DO
  END IF


  !      do ih=1,nh(1)
  !         do jh=1,nh(1)
  !            do ipol=1,3
  !            WRITE( stdout,'(3i5,2f20.10)') ipol,ih,jh,int2(ih,jh,ipol,1,1)
  !            enddo
  !         enddo
  !      enddo
  !      call stop_ph(.true.)
  deallocate (veff)
  if (.not.lgamma) then
     deallocate(qgmq)
     deallocate (qmod)
     deallocate (ylmkq)
  endif
  deallocate (qgm)
  deallocate (ylmk0)
  deallocate (qmodg)
  deallocate(aux35)
  deallocate (aux1)
  deallocate (sk)

  call stop_clock ('dvanqq')
  return
end subroutine dvanqq
