!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------
subroutine dvpsi_e2
  !-----------------------------------------------------------------------
  !
  ! This routine shold be called before the self-consistent cycle used to
  ! compute the second derivative of the wavefunctions with respect to
  ! electric-fields. It computes that part of the potential that remains
  ! constant during the cycle.
  !
#include "f_defs.h"

  use kinds, only : DP
  use pwcom
  use scf, only : rho
  USE io_files, ONLY: iunigk
  USE wavefunctions_module,  ONLY: evc
  use becmod
  use phcom
  USE ramanm
#ifdef __PARA
  USE mp_global, ONLY: my_pool_id
#endif
  implicit none

  integer :: ik, ipa, ipb, ir, ibnd, jbnd, nrec
  ! counter on k-points
  ! counter on polarizations
  ! counter on points of the real-space mesh
  ! counter on bands
  ! the record number
  real(DP), allocatable  :: raux6 (:,:), d2muxc (:)
  ! function on the real space smooth-mesh
  ! second derivative of the XC-potential 
  real(DP) ::  d2mxc, rhotot
  ! external function
  ! total charge on a point
  complex(DP), allocatable :: depsi (:,:,:), auxg (:,:), auxs1 (:), &
               auxs2 (:), aux3s (:,:), aux3 (:,:), ps (:,:,:,:)
  ! d |psi> / dE  (E=electric field)
  ! chi-wavefunction
  ! function on the real space smooth-mesh
  ! function on the real space smooth-mesh
  ! function on the real space smooth-mesh
  ! function on the real space thick-mesh
  complex(DP), pointer :: aux6s (:,:), aux6 (:,:)
  ! function on the real space smooth-mesh
  ! function on the real space thick-mesh
  complex(DP) :: tmp, weight, ZDOTC
  ! working space              
  ! weight in k-point summation
  ! the scalar product function
  !
  call start_clock('dvpsi_e2')
  !
  ! First, calculates the second derivative of the charge-density
  ! -only the part that does not depend on the self-consistent cycle-
  !
  allocate (raux6  (nrxxs,6))
  allocate (depsi  (npwx,nbnd,3))
  allocate (aux3s  (nrxxs,3))
  allocate (ps     (nbnd,nbnd,3,3))

  raux6 (:,:) = 0.d0
  if (nksq.gt.1) rewind (iunigk)
  do ik = 1, nksq
     if (nksq.gt.1) then
        read (iunigk) npw, igk
        npwq = npw
        call davcio (evc, lrwfc, iuwfc, ik, -1)
     endif
     weight = 2.d0 * wk(ik) / omega

     do ipa = 1, 3
        nrec = (ipa - 1) * nksq + ik
        call davcio (depsi (1, 1, ipa), lrdwf, iudwf, nrec, -1)
     enddo

     do ibnd = 1, nbnd_occ (ik)
        do ipa = 1, 3
           call cft_wave (depsi (1, ibnd, ipa), aux3s (1, ipa), +1)
        enddo
        do ipa = 1, 6
           do ir = 1, nrxxs
              tmp = CONJG(aux3s (ir, a1j (ipa))) *    &
                           aux3s (ir, a2j (ipa))
              raux6 (ir, ipa) = raux6 (ir, ipa) + weight *  DBLE (tmp)
           enddo
        enddo
     enddo

     do ipa = 1, 3
        do ipb = 1, 3
           CALL ZGEMM( 'C', 'N', nbnd_occ (ik), nbnd_occ (ik), npwq, &
                (1.d0,0.d0), depsi(1,1, ipa), npwx, depsi(1,1,ipb), npwx, &
                (0.d0,0.d0), ps(1,1,ipa,ipb), nbnd )
        enddo
     enddo
#ifdef __PARA
     call reduce (2 * nbnd * nbnd * 9, ps)
#endif

     do ibnd = 1, nbnd_occ (ik)
        call cft_wave (evc (1, ibnd), aux3s (1,1), +1)
        do jbnd = 1, nbnd_occ (ik)
           call cft_wave (evc (1, jbnd), aux3s (1,2), +1)
           do ipa = 1, 6
              do ir = 1, nrxxs
                 tmp =  aux3s (ir,1) *                           &
                        ps(ibnd, jbnd, a1j (ipa), a2j (ipa)) *   &
                        CONJG(aux3s (ir,2))
                 raux6 (ir, ipa) = raux6 (ir, ipa) - weight *  DBLE (tmp)
              enddo
           enddo
        enddo
     enddo

  enddo

  deallocate (depsi)
  deallocate (aux3s)
  deallocate (ps)

  !
  ! Multiplies the charge with the potential
  !
  if (doublegrid) then
     allocate (auxs1  (nrxxs))
     allocate (aux6   (nrxx,6))
  else
     allocate (aux6s  (nrxxs,6))
     aux6 => aux6s
  endif

  do ipa = 1, 6
     if (doublegrid) then
        do ir = 1, nrxxs
           auxs1 (ir) = CMPLX (raux6 (ir, ipa), 0.d0)
        enddo
        call cinterpolate (aux6 (1, ipa), auxs1, +1)
     else
        do ir = 1, nrxxs
           aux6 (ir, ipa) = CMPLX (raux6 (ir, ipa), 0.d0)
        enddo
     endif
     call dv_of_drho (0, aux6(1, ipa), .false.)
  enddo

  if (doublegrid) deallocate (auxs1)
  deallocate (raux6)

  !
  ! Calculates the term depending on the third derivative of the
  !                     Exchange-correlation energy
  !
  allocate (d2muxc (nrxx))
  allocate (aux3   (nrxx,3))
  do ipa = 1, 3
     call davcio_drho (aux3 (1, ipa), lrdrho, iudrho, ipa, -1)
  enddo

#ifdef __PARA
  if (my_pool_id .ne. 0) goto 100
#endif
  d2muxc (:) = 0.d0
  do ir = 1, nrxx
!     rhotot = rho%of_r(ir,1) + rho_core(ir)
     rhotot = rho%of_r(ir,1)
     if ( rhotot.gt. 1.d-30 ) d2muxc(ir)= d2mxc( rhotot)
     if ( rhotot.lt.-1.d-30 ) d2muxc(ir)=-d2mxc(-rhotot)
  enddo

  do ipa = 1, 6
     do ir = 1, nrxx
        aux6 (ir, ipa) = aux6 (ir, ipa) + d2muxc (ir) * &
                   aux3 (ir, a1j (ipa)) * aux3 (ir, a2j (ipa))
     enddo
  enddo
#ifdef __PARA
 100  continue
  call poolreduce (2 * 6 * nrxx, aux6)
  call psyme2 (aux6)
#else
  call syme2 (aux6)
#endif
  deallocate (d2muxc) 
  deallocate (aux3)


  if (doublegrid) then
     allocate (aux6s  (nrxxs,6))
     do ipa = 1, 6
        call cinterpolate (aux6 (1, ipa), aux6s (1, ipa), -1)
     enddo
     deallocate (aux6)
  endif
  !
  ! Multiplies the obtained potential with the wavefunctions and
  ! writes the results on iuba2; a faster way of proceeding would
  ! be that of keeping the potential in memory and use it directly in
  ! solve_e2
  !
  allocate (auxg   (npwx,nbnd))
  allocate (auxs1  (nrxxs))
  allocate (auxs2  (nrxxs))

  if (nksq.gt.1) rewind (iunigk)
  do ik = 1, nksq
     if (nksq.gt.1) then
        read (iunigk) npw, igk
        npwq = npw
        call davcio (evc, lrwfc, iuwfc, ik, -1)
     endif
     do ipa = 1, 6
        nrec = (ipa - 1) * nksq + ik
        call davcio (auxg, lrchf, iuchf, nrec, -1)
        do ibnd = 1, nbnd_occ (ik)
           call cft_wave (evc (1, ibnd), auxs1, +1)
           do ir = 1, nrxxs
              auxs2 (ir) = auxs1 (ir) * aux6s (ir, ipa)
           enddo
           call cft_wave (auxg (1, ibnd), auxs2, -1)
        enddo
        nrec = (ipa - 1) * nksq + ik
        call davcio (auxg, lrba2, iuba2, nrec, +1)
     enddo
  enddo

  deallocate (auxg)
  deallocate (auxs1)
  deallocate (auxs2)

  deallocate (aux6s)
  call stop_clock('dvpsi_e2')

  return
end subroutine dvpsi_e2
