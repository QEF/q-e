!
! Copyright (C) 2001-2007 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
subroutine raman_mat
  !-----------------------------------------------------------------------
  !
  ! Reads on the disk all the necessary wavefunctions and computes
  !                      the raman tensor
  !

  USE kinds,    ONLY : DP
  USE becmod,   ONLY : calbec
  USE constants,ONLY : e2, fpi
  USE cell_base,ONLY : at, bg, omega, tpiba
  USE gvect,    ONLY : g
  USE klist,    ONLY : wk, xk, ngk, igk_k
  USE buffers,  ONLY : get_buffer
  USE ions_base,ONLY : nat
  USE symme,    ONLY : symtensor3
  USE uspp,     ONLY : nkb, vkb
  USE wvfct,    ONLY : npwx, nbnd
  USE wavefunctions_module,  ONLY: evc
  USE phus,     ONLY : alphap
  USE units_ph, ONLY : lrdwf, iudwf, lrwfc, iuwfc
  USE ramanm,   ONLY : ramtns, jab, a1j, a2j, lrd2w, iud2w

  USE lrus,     ONLY : becp1
  USE control_lr, ONLY : nbnd_occ
  USE qpoint,   ONLY : nksq
  USE eqv,      ONLY : dvpsi

  USE mp_pools, ONLY : inter_pool_comm
  USE mp_bands, ONLY : intra_bgrp_comm
  USE mp,       ONLY : mp_sum
  implicit none

  logical :: wr_all
  integer :: ik, ig, ipa, ipb, icr, jcr, iat, ibnd, jbnd, imod, nrec, &
             il, ntm, ipol, npw, npwq
  ! counter on k-points
  ! counter on electric field polarizations
  ! counter on electric field polarizations
  ! counter on cartesian coordinates
  ! counter on cartesian coordinates
  ! counter on atoms
  ! counter on bands
  ! counter on atomic displacement mode
  ! record number

  real(DP) , allocatable :: wrk (:,:,:), matram(:,:,:,:), matw(:,:,:,:,:)
  ! working array
  ! the Raman-tensor: the first two indexes referes to the electric fields,
  !                   the last two to the atomic displacemts
  ! components of the Raman-tensor: used only for testing purposes
  real(DP) :: weight, tmp
  ! weight in the summation over k-points
  ! working space

  complex(DP) , allocatable :: uact(:,:), chif(:,:,:),            &
                     depsi (:,:,:), auxg(:), evc_sw (:,:), aux1 (:,:), &
                     ps (:,:,:,:), becp1_sw (:,:), alphap_sw (:,:,:)
  ! pattern of atomic displacements
  ! array of wavefunctions
  ! swap space
  complex(DP) :: tmpc
  ! the scalar product function
  complex(DP), EXTERNAL :: zdotc

  allocate (wrk       (6,3*nat,2)   )
  allocate (matram    (3,3,3,nat)   )
  allocate (matw      (3,3,3,nat,4) )
  allocate (uact      (3*nat,3*nat) )
  allocate (chif      (npwx,nbnd,6) )
  allocate (depsi     (npwx,nbnd,3) )
  allocate (auxg      (npwx)        )
  allocate (evc_sw    (npwx,nbnd)   )
  allocate (aux1      (npwx,nbnd)   )
  allocate (ps       (nbnd,nbnd,3,3))
  allocate (becp1_sw  (nkb,nbnd)    )
  allocate (alphap_sw (nkb,nbnd,3)  )

  !
  ! Set the atomic displacement pattern ( crystal coordinates )
  !
  uact (:,:) = (0.d0, 0.d0)
  do iat = 0, nat - 1
     do icr = 1, 3
        do jcr = 1, 3
           uact (3*iat + jcr, 3*iat + icr) = CMPLX(at (jcr, icr), 0.d0,kind=DP)
        enddo
     enddo
  enddo

  wrk (:,:,:) = 0.d0

  !
  ! The raman tensor is computed as the sum of three different contribution
  ! These contributions are calculated in the following loop and stored
  ! in the two different arrays wrk(:,:,i),i=1,2 ( this may be usefull while
  ! testing ).
  !
  do ik = 1, nksq
  !
  ! Using weight = 2.d0*wk(ik)*e2, calculates the third derivative of the
  !   energy with respect to atomic displacemements and with respect to two
  !   electric fields (units are Bohr^2).
  ! Using weight = -2.d0*wk(ik)*e2*fpi/omega, calculates the derivative
  !   of the dielectric constants with respect to atomic-displacem
  !   (units are Bohr^-1 ).
  !       weight = -2.d0*wk(ik)*e2
     weight = - 2.d0 * wk (ik) * e2 * fpi / omega
     npw = ngk(ik)
     npwq = npw
     if (nksq.gt.1) call get_buffer (evc, lrwfc, iuwfc, ik)
     call init_us_2 (npw, igk_k(1,ik), xk (1,ik), vkb)

     do ipa = 1, 3
        nrec = (ipa - 1) * nksq + ik
        call get_buffer(depsi (1, 1, ipa), lrdwf, iudwf, nrec)
     enddo
     do ipa = 1, 3
        do ipb = 1, 3
           do ibnd = 1, nbnd_occ (ik)
              do jbnd = 1, nbnd_occ (ik)
                 ps (ibnd, jbnd, ipa, ipb) =                &
                      zdotc (npwq, depsi (1, ibnd, ipa), 1, &
                             depsi (1, jbnd, ipb), 1)
              enddo
           enddo
        enddo
     enddo

     call mp_sum ( ps, intra_bgrp_comm )

     do ipa = 1, 6
        nrec = (ipa - 1) * nksq + ik
        call davcio (chif (1, 1, ipa), lrd2w, iud2w, nrec, -1)
     enddo

     do ipa = 1, 6
        do ibnd = 1, nbnd_occ (ik)
           auxg (:) = (0.d0, 0.d0)
           do jbnd = 1, nbnd_occ (ik)
              tmpc = ps (jbnd, ibnd, a1j (ipa), a2j (ipa))
              call zaxpy (npwq, tmpc, evc (1, jbnd), 1, auxg, 1)
           enddo
           call daxpy (2 * npwq, -1.d0, auxg, 1, chif (1, ibnd, ipa), 1)
        enddo
     enddo

     do imod = 1, 3 * nat
        call dvqpsi_us (ik, uact (1, imod),.false. )
        do ipa = 1, 6
           tmp = 0.d0
           do ibnd = 1, nbnd_occ (ik)
              tmp = tmp + weight *  DBLE( zdotc(npwq,             &
                    chif (1, ibnd, ipa), 1, dvpsi (1, ibnd), 1) )
           enddo
           wrk (ipa, imod, 1) = wrk (ipa, imod, 1) + tmp
        enddo
     enddo

     !
     ! evc, becp1, alphap  are written into a swap space
     !
     if (nksq.eq.1) call zcopy (npwx * nbnd, evc, 1, evc_sw, 1)
     call zcopy (nkb * nbnd, becp1(ik)%k, 1, becp1_sw,  1)
     DO ipol=1,3
        call zcopy (nkb * nbnd, alphap (ipol, ik)%k, 1, alphap_sw(1,1,ipol), 1)
     ENDDO

     do ipa = 1, 3
        nrec = (ipa - 1) * nksq + ik
        call get_buffer(chif (1, 1, ipa), lrdwf, iudwf, nrec)
     enddo

     do imod = 1, 3 * nat
        do ipa = 1, 3
           !
           ! initializes some variables used by dvqpsi_us
           !
           call zcopy (npwx * nbnd, chif (1, 1, ipa), 1, evc, 1)
           call calbec (npw, vkb, evc, becp1(ik) )
           do ipb = 1, 3
              do ibnd = 1, nbnd
                 do ig = 1, npw
                    aux1 (ig, ibnd) = evc(ig,ibnd) *               &
                                      tpiba * (0.d0,1.d0) *        &
                                ( xk(ipb,ik) + g(ipb,igk_k(ig,ik)) )
                 enddo
              enddo
              call calbec (npw, vkb, aux1, alphap (ipb,ik) )
           enddo

           call dvqpsi_us (ik, uact (1, imod),.false. )
           do ipb = 1, ipa
              tmp = 0.d0
              do ibnd = 1, nbnd_occ (ik)
                 tmp = tmp + weight *  DBLE(zdotc (npwq, &
                       chif (1, ibnd, ipb), 1, dvpsi (1, ibnd), 1) )
              enddo
              wrk (jab (ipa, ipb), imod, 2) = &
              wrk (jab (ipa, ipb), imod, 2) + tmp
           enddo
        enddo
     enddo
     !
     ! evc, becp1, alphap are restored to their original value
     !
     if (nksq.eq.1) call zcopy (npwx * nbnd, evc_sw, 1, evc, 1)
     call zcopy (nkb * nbnd, becp1_sw,  1, becp1(ik)%k, 1)
     do ipol=1,3
        call zcopy (nkb * nbnd,  alphap_sw(1,1,ipol), 1, alphap(ipol, ik)%k, 1)
     enddo

  enddo

  call mp_sum( wrk, intra_bgrp_comm )
  call mp_sum( wrk, inter_pool_comm )

  do iat = 1, nat
     do icr = 1, 3
        imod = icr + (iat - 1) * 3
        do ipa = 1, 3
           do ipb = 1, ipa
              tmp = wrk (jab (ipa, ipb), imod, 1) + &
                    wrk (jab (ipa, ipb), imod, 2)
              matw (ipa, ipb, icr, iat, 1) = tmp
              matw (ipb, ipa, icr, iat, 1) = tmp
              do il = 1, 2
                 matw (ipa, ipb, icr, iat, il + 1) = &
                       wrk (jab (ipa, ipb), imod, il)
                 matw (ipb, ipa, icr, iat, il + 1) = &
                       wrk (jab (ipa, ipb), imod, il)
              enddo
           enddo
        enddo
     enddo
  enddo
  !
  ! wr_all =.true. ==> writes the two contributions before and after
  !  symmetrization (used for testing purposes only )
  !
  wr_all = .false.

  ntm = 1
  if (wr_all ) ntm = 3

  do il = 1, ntm
     call dcopy(27*nat,matw(1,1,1,1,il),1,matram,1)
     if (wr_all ) then
        if (il.eq.1) then
           write(6,'(/,10x,''Raman tensor: Total '',/)')
        else
           write(6,'(/,10x,''Raman tensor: contribution # '',i3,/)') &
              il - 1
        endif
        write(6,'(/,10x,''Unsymmetrized in crystal axis '',/)')
        call write_raman(matram)
     endif
     !
     ! Symmetrizes the Raman tensor
     ! NOte that the output matrix is in cartesian axis
     !
     call symtensor3 ( nat, matram )
     if (wr_all ) then
        write(6,'(/,10x,''Symmetrized in cartesian axis '',/)')
        call write_raman(matram)
     endif
     !
     write(6,'(/,10x,''Raman tensor (au^-1) in cartesian axis '',/)')
     !
     if (il == 1) ramtns(:,:,:,:) = matram(:,:,:,:)
     if (wr_all ) call write_raman(matram)

     do iat = 1, nat
        write(6,'(10x,'' atom '',i6)') iat
        do icr = 1, 3
           do ipb = 1, 3
              write(6,'(10x,''('',3f18.9,'' )'')')      &
                      (matram(ipa,ipb,icr,iat),ipa=1,3)
           enddo
           write(6,'(10x)')
        enddo
     enddo
  enddo
  !
  ! write Raman tensor dchi/du = (omega/4pi)*deps/du in A^2
  ! it may not be written to file fildyn if trans=.false.
  !
  call write_ramtns (6, ramtns)
  !
  deallocate (wrk       )
  deallocate (matram    )
  deallocate (matw      )
  deallocate (uact      )
  deallocate (chif      )
  deallocate (depsi     )
  deallocate (auxg      )
  deallocate (evc_sw    )
  deallocate (aux1      )
  deallocate (ps        )
  deallocate (becp1_sw  )
  deallocate (alphap_sw )

  return
end subroutine raman_mat
!
!-----------------------------------------------------------------------
subroutine write_raman (matram)
  !-----------------------------------------------------------------------
  !

  use kinds, only : DP
  USE ions_base, ONLY: nat
  USE ramanm, ONLY : a1j, a2j
  implicit none

  real(DP) :: matram(3,3,3,nat)

  integer :: icr, iat, ipa
  character (len=2) :: ch(3), ch2(6)

  data ch /'X','Y','Z'/
  data ch2 /'XX','YY','ZZ','XY','YZ','ZX'/

  write(6,'(''  at'',7x,3(a2,10x),3x,3(a2,10x)  )') &
           ( ch2 (icr), icr = 1, 6)
  do iat = 1, nat
     write(6,'(1x)')
     do icr = 1, 3
        write(6,'(1x,i3,1x,a1,'':'',3f12.6,3x,3f12.6)')      &
                 iat, ch (icr),                              &
                 (matram (a1j (ipa), a2j (ipa), icr, iat), ipa=1,6)
     enddo
  enddo

  return
end subroutine write_raman
