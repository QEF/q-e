!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------

subroutine phq_init
  !-----------------------------------------------------------------------
  !
  !     This subroutine computes the quantities necessary to describe the
  !     local and nonlocal pseudopotential in the phononq program.
  !     In detail it computes:
  !     0) initialize the structure factors
  !     a0) compute rhocore for each atomic-type if needed for nlcc
  !     a) The local potential at G-G'. Needed for the part of the dynamic
  !        matrix independent of deltapsi.
  !     b) The local potential at q+G-G'. Needed for the second
  !        second part of the dynamical matrix.
  !     c) The D coefficients for the US pseudopotential or the E_l parame
  !        of the KB pseudo. In the US case it prepares also the integrals
  !        qrad and qradq which are needed for computing Q_nm(G) and
  !        Q_nm(q+G)
  !     d) The functions vkb(k+G) needed for the part of the dynamical mat
  !        independent of deltapsi.
  !     e) The becp functions for the k points
  !     e') The derivative of the becp term with respect to a displacement
  !     f) The functions vkb(k+q+G), needed for the linear sysetm and the
  !        second part of the dynamical matrix.
  !
#include"machine.h"
  USE io_global,      ONLY : stdout
  use pwcom
  USE wavefunctions,  ONLY: evc
  use parameters, only : DP
  use phcom
  implicit none

  integer :: nt, ik, ikq, ipol, ibnd, ikk, na, ig
  ! counter on atom types
  ! counter on k points
  ! counter on k+q points
  ! counter on polarizations
  ! counter on bands
  ! index for wavefunctions at k
  ! counter on atoms
  ! counter on G vectors

  real(kind=DP) :: arg
  ! the argument of the phase


  complex(kind=DP), allocatable :: aux1 (:,:)
  ! used to compute alphap

  call start_clock ('phq_init')
  allocate (aux1( npwx , nbnd))    
  !
  !  initialize structure factor array
  !
  call struc_fact (nat, tau, ntyp, ityp, ngm, g, bg, nr1, nr2, nr3, &
       strf, eigts1, eigts2, eigts3)
  do na = 1, nat
     arg = (xq (1) * tau (1, na) + &
            xq (2) * tau (2, na) + &
            xq (3) * tau (3, na) ) * tpi
     eigqts (na) = DCMPLX (cos (arg), - sin (arg) )
  enddo
  !
  !   a0) compute rhocore for each atomic-type if needed for nlcc
  !

  if (nlcc_any) call set_drhoc (xq)
  !
  !   a) the fourier components of the local potential for each |G|
  !
  call init_vloc
  !
  !   b) the fourier components of the local potential at q+G
  !
  vlocq(:,:) = 0.d0
  do nt = 1, ntyp
     call setlocq (xq, lloc (nt), lmax (nt), numeric (nt), mesh (nt), &
          msh (nt), rab (1, nt), r (1, nt), vnl (1, lloc (nt), nt), &
          cc (1, nt), alpc (1, nt), nlc (nt), nnl (nt), zp (nt), &
          aps(1, 0, nt), alps(1, 0, nt), tpiba2, ngm, g, omega, vlocq(1, nt) )
  enddo
  !
  !   c) the parameters defining the pseudopotential
  !
  !     for the analytic potentials we need to convert in a radial mesh
  !
  call convert_to_num (ntyp, numeric, ndm, mesh, r, lmaxx, lmax, &
       lloc, nnl, aps, alps, vnl)
  !
  !   then we compute the denominators of the KB types, or the
  !   parameters which define the non-local pseudopotential and
  !   which are independent of the k point for the US case
  !

  call init_us_1

  if (nksq.gt.1) rewind (iunigk)
  do ik = 1, nksq
     if (lgamma) then
        ikk = ik
        ikq = ik
        npwq = npw
     else
        ikk = 2 * ik - 1
        ikq = ikk + 1
     endif
     !
     if (lsda) current_spin = isk (ikk)
     !
     !  g2kin is used here as work space
     !
     call gk_sort (xk (1, ikk), ngm, g, ecutwfc / tpiba2, npw, igk, g2kin)
     !
     !  if there is only one k-point evc, evq, npw, igk stay in memory
     !
     if (nksq.gt.1) write (iunigk) npw, igk
     if (.not.lgamma) then
        call gk_sort (xk (1, ikq), ngm, g, ecutwfc / tpiba2, npwq, igkq, g2kin)
        if (nksq.gt.1) write (iunigk) npwq, igkq
        if (abs (xq (1) - (xk (1, ikq) - xk (1, ikk) ) ) .gt.1.d-8 .or. &
            abs (xq (2) - (xk (2, ikq) - xk (2, ikk) ) ) .gt.1.d-8 .or. &
            abs (xq (3) - (xk (3, ikq) - xk (3, ikk) ) ) .gt.1.d-8) then
           WRITE( stdout, * ) ikk, ikq, nksq
           WRITE( stdout, * ) (xq (ipol), ipol = 1, 3)
           WRITE( stdout, * ) (xk (ipol, ikq), ipol = 1, 3)
           WRITE( stdout, * ) (xk (ipol, ikk), ipol = 1, 3)
           call errore ('phq_init', 'wrong order of k points', 1)
        endif
     endif
     !
     !   d) The functions vkb(k+G)
     !
     call init_us_2 (npw, igk, xk (1, ikk), vkb)
     !
     ! read the wavefunctions at k
     !
     call davcio (evc, lrwfc, iuwfc, ikk, - 1)
     !
     !     e) we compute the becp terms which are used in the rest of
     !     the code
     call ccalbec (nkb, npwx, npw, nbnd, becp1 (1, 1, ik), vkb, evc)
     !
     !     e') we compute the derivative of the becp term with respect to an
     !         atomic displacement
     !
     do ipol = 1, 3
        do ibnd = 1, nbnd
           do ig = 1, npw
              aux1 (ig, ibnd) = evc(ig,ibnd) * tpiba * (0.d0,1.d0) * & 
                                ( xk(ipol,ikk) + g(ipol,igk(ig)) )
           enddo
        enddo
        call ccalbec (nkb, npwx, npw, nbnd, alphap(1,1,ipol,ik), vkb, aux1)
     enddo

     !
     ! if there is only one k-point the k+q wavefunctions are read once here
     !
     if (nksq.eq.1.and..not.lgamma) call davcio (evq, lrwfc, iuwfc, ikq, -1)

  enddo

  deallocate (aux1)

  call newd
  call dvanqq
  call drho
  if ((epsil.or.zue).and.okvan) then
     call compute_qdipol
  endif
  call stop_clock ('phq_init')
  return
end subroutine phq_init
