!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------
subroutine d0rhod2v (ipert, drhoscf)
!-----------------------------------------------------------------------
! calculates the term containing the second variation of the potential
! and the first variation of the charge density with respect to a
! perturbation at q=0
!
#include "machine.h"
  use pwcom
  USE wavefunctions,  ONLY: evc
  use phcom
  use d3com
#ifdef __PARA
  use para
#endif
  implicit none
  integer :: ipert              ! index of the perturbation associated with drho
  complex (8) :: drhoscf (nrxx) ! the variation of the charge density
!
  integer :: icart,           & ! counter on polarizations
             jcart,           & ! counter on polarizations
             na_icart,        & ! counter on modes
             na_jcart,        & ! counter on modes
             na,              & ! counter on atoms
             ng,              & ! counter on G vectors
             nt,              & ! counter on atomic types
             ik,              & ! counter on k points
             ikk,             & ! counter on k points
             ig,              & ! counter on G vectors
             ibnd,            & ! counter on bands
             nu_i,            & ! counter on modes
             nu_j,            & ! counter on modes
             nu_k,            & ! counter on modes
             ikb, jkb,        & ! counter on beta functions
             nrec,            & ! record position of dwfc
             ios                ! integer variable for I/O control

  real (8) :: gtau,           & ! the product G*\tau_s
              wgg               ! the weight of a K point

  complex (8) :: ZDOTC, d3dywrk (3*nat,3*nat), fac, alpha(8), work
  complex (8), allocatable :: work0 (:), work1 (:), work2 (:), work3 (:), &
                          work4 (:), work5 (:), work6 (:)
  ! auxiliary space

  allocate (work0(nrxx))    
  allocate (work1(npwx))    
  allocate (work2(npwx))    
  allocate (work3(npwx))    
  allocate (work4(npwx))    
  allocate (work5(npwx))    
  allocate (work6(npwx))    

  call setv (2*9*nat*nat,0.0d0,d3dywrk,1)
!
! Here the contribution deriving from the local part of the potential
#ifdef __PARA
!   ... computed only by the first pool (no sum over k needed)
!
  if (mypool.ne.1) goto 100
#endif
!
  call ZCOPY (nrxx, drhoscf, 1, work0, 1)
  call cft3 (work0, nr1, nr2, nr3, nrx1, nrx2, nrx3, -1)
  do na = 1, nat
     do icart = 1,3
        na_icart = 3*(na-1)+icart
        do jcart = 1,3
           na_jcart = 3*(na-1)+jcart
           do ng = 1, ngm
              gtau = tpi * ( g(1,ng)*tau(1,na) + &
                             g(2,ng)*tau(2,na) + &
                             g(3,ng)*tau(3,na) )

              fac = DCMPLX(cos(gtau),sin(gtau))

              d3dywrk(na_icart,na_jcart) = &
                      d3dywrk(na_icart,na_jcart) - &
                      tpiba2 * g(icart,ng) * g(jcart,ng) * &
                      omega * vloc(igtongl(ng),ityp(na)) * &
                      fac*work0(nl(ng))
           enddo
        enddo
     enddo
     write (*,*) na
     write (*,'(3(2f10.6,2x))') &
           ((d3dywrk(3*(na-1)+icart,3*(na-1)+jcart), &
           jcart=1,3),icart=1,3)

  enddo
#ifdef __PARA
  call reduce(2*9*nat*nat,d3dywrk)
!
! each pool contributes to next term
!
100 continue
#endif
!
! Here we compute the nonlocal (Kleinman-Bylander) contribution.
!
  rewind (unit=iunigk)

  do ik = 1, nksq
     read (iunigk, err = 200, iostat = ios) npw, igk
200  call errore ('d0rhod2v', 'reading igk', abs (ios) )
     if (lgamma) then
        ikk = ik
        npwq = npw
     else
        ikk = 2 * ik - 1
        read (iunigk, err = 300, iostat = ios) npwq, igkq
300     call errore ('d0rhod2v', 'reading igkq', abs (ios) )
        npwq = npw
     endif
     wgg = wk (ikk)
     call davcio (evc, lrwfc, iuwfc, ikk, - 1)

     call init_us_2 (npw, igk, xk (1, ikk), vkb0)
     !
     ! Reads the first variation of the wavefunction projected on conduction
     !
     nrec = (ipert - 1) * nksq + ik
     call davcio (dpsi, lrdwf, iudwf, nrec, - 1)
     !
     ! In the metallic case corrects dpsi so as that the density matrix
     ! will be:   Sum_{k,nu} 2 * | dpsi > < psi |
     !
     if (degauss.ne.0.d0) then
        nrec = ipert + (ik - 1) * 3 * nat
        call davcio (psidqvpsi, lrpdqvp, iupd0vp, nrec, - 1)
        call dpsi_corr (evc, psidqvpsi, ikk, ikk, ipert)
     endif
     do icart = 1, 3
        do jcart = 1, 3
           do ibnd = 1, nbnd
              do ig = 1, npw
                 work1(ig)= evc(ig,ibnd)*tpiba*(xk(icart,ikk)+g(icart,igk(ig)))
                 work2(ig)= evc(ig,ibnd)*tpiba*(xk(jcart,ikk)+g(jcart,igk(ig)))
                 work3(ig)=dpsi(ig,ibnd)*tpiba*(xk(icart,ikk)+g(icart,igk(ig)))
                 work4(ig)=dpsi(ig,ibnd)*tpiba*(xk(jcart,ikk)+g(jcart,igk(ig)))
                 work5(ig)=    work1(ig)*tpiba*(xk(jcart,ikk)+g(jcart,igk(ig)))
                 work6(ig)=    work3(ig)*tpiba*(xk(jcart,ikk)+g(jcart,igk(ig)))
              enddo
              jkb=0
              do nt = 1, ntyp
                 do na = 1, nat
                    if (ityp (na).eq.nt) then
                       na_icart = 3 * (na - 1) + icart
                       na_jcart = 3 * (na - 1) + jcart
                       do ikb = 1, nh (nt)
                          jkb=jkb+1
                          alpha (1) = ZDOTC (npw, work1, 1, vkb0(1,jkb), 1)
                          alpha (2) = ZDOTC (npw, vkb0(1,jkb), 1, work4, 1)
                          alpha (3) = ZDOTC (npw, work2, 1, vkb0(1,jkb), 1)
                          alpha (4) = ZDOTC (npw, vkb0(1,jkb), 1, work3, 1)
                          alpha (5) = ZDOTC (npw, work5, 1, vkb0(1,jkb), 1)
                          alpha (6) = ZDOTC (npw, vkb0(1,jkb), 1, dpsi (1,ibnd), 1)
                          alpha (7) = ZDOTC (npw,  evc (1,ibnd), 1, vkb0(1,jkb), 1)
                          alpha (8) = ZDOTC (npw, vkb0(1,jkb), 1, work6, 1)
#ifdef __PARA
                          call reduce (16, alpha)
#endif
                          d3dywrk (na_icart, na_jcart) = d3dywrk (na_icart, na_jcart) &
                               + (alpha(1)*alpha(2) + alpha(3)*alpha(4) &
                                - alpha(5)*alpha(6) - alpha(7)*alpha(8)) * &
                                  dvan (ikb,ikb,nt) * wgg * 2.0d0
                       enddo
                    end if
                 enddo
              enddo
           enddo
        enddo
     enddo
  enddo
#ifdef __PARA
  call poolreduce (2*9*nat*nat, d3dywrk)
#endif
!
!  Rotate the dynamical matrix on the basis of patterns
!  first index does not need to be rotated
!
  nu_k = ipert
  do nu_i = 1, 3 * nat
     do nu_j = 1, 3 * nat
        work = (0.0d0, 0.0d0)
        do na = 1, nat
           do icart = 1, 3
              na_icart = 3 * (na-1) + icart
              do jcart = 1, 3
                 na_jcart = 3 * (na-1) + jcart
                 work = work + conjg(u(na_icart,nu_i)) * &
                               d3dywrk(na_icart,na_jcart) * &
                               u(na_jcart,nu_j)
              enddo
           enddo
        enddo
        d3dyn(nu_k,nu_i,nu_j) = d3dyn(nu_k,nu_i,nu_j) + work
        if (allmodes) then
           d3dyn(nu_j,nu_k,nu_i) = d3dyn(nu_j,nu_k,nu_i) + work
           d3dyn(nu_i,nu_j,nu_k) = d3dyn(nu_i,nu_j,nu_k) + work
        endif
     enddo

  enddo
  deallocate (work6)
  deallocate (work5)
  deallocate (work4)
  deallocate (work3)
  deallocate (work2)
  deallocate (work1)
  deallocate (work0)

  return
end subroutine d0rhod2v
