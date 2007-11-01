!
! Copyright (C) 2001-2007 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
subroutine el_opt
  !-----------------------------------------------------------------------
  !
  ! Calculates electro-optic tensor
  !
#include "f_defs.h"
  use kinds, only : DP
  USE cell_base,  ONLY : omega, at, bg
  USE constants,  ONLY : e2, fpi
  USE gvect,      ONLY : nr1, nr2, nr3, nrxx
  USE klist,      ONLY : wk, ngk
  USE ions_base,  ONLY : nat
  USE scf,        ONLY : rho, rho_core
  USE symme,      ONLY : nsym, irt, s
  USE wvfct,      ONLY : nbnd, npw, npwx
  use phcom
  USE ramanm
  USE io_global, ONLY: ionode_id
#ifdef __PARA
  USE mp, ONLY: mp_bcast
  USE mp_global, ONLY: my_pool_id
#endif
  implicit none

  logical wr_all
  integer :: ik, ir, ipa, ipb, ipc, nrec, ibnd, jbnd, il, ntm
  real(DP) :: weight, fac, elop_ (3, 3, 3, 3), ps3 (3, 3, 3)
  real(DP) :: d2mxc, rhotot
  ! external function
  ! total charge on a point
  real(DP), allocatable :: d2muxc (:)
  complex(DP) :: ps(3, 6), tmp
  complex(DP) , allocatable :: chif(:,:,:), depsi (:,:,:), aux3 (:,:)

  call start_clock('el_opt')
  elop_(:,:,:,:) = 0.0_dp

  allocate (depsi(npwx, nbnd, 3) )
  allocate (chif (npwx, nbnd, 6) )

  do ik = 1, nksq
     weight = wk(ik)
     npw = ngk(ik)
     do ipa = 1, 3
        nrec = (ipa - 1) * nksq + ik
        call davcio (depsi(1,1,ipa), lrdwf, iudwf, nrec, -1)
     enddo
     do ipb = 1, 6
        nrec = (ipb - 1) * nksq + ik
        call davcio (chif(1,1,ipb), lrchf, iuchf, nrec, -1)
     enddo

     ! ps (ipa,ipb) = \sum_i < depsi_i(ipa) | chif_i(ipb) >
     !       do ibnd = 1, nbnd_occ (ik)
     !             ps (ipa, ipb) =  ps (ipa, ipb) +           &
     !                  ZDOTC (npw, depsi (1, ibnd, ipa), 1,  &
     !                               chif (1, ibnd, ipb), 1 )
     !       end do

     CALL ZGEMM( 'C', 'N', 3, 6, npwx*nbnd_occ(ik), (1.0_dp,0.0_dp),  &
                 depsi, npwx*nbnd, chif, npwx*nbnd, &
                 (0.0_dp,0.0_dp), ps, 3 )

     do ipa = 1, 3
        do ipb = 1, 3
           do ipc = 1, 3
              elop_ (ipa, ipb, ipc, 1) = elop_ (ipa, ipb, ipc, 1) + &
                  weight *  DBLE( ps(ipa, jab (ipb, ipc))   + &
                                  ps(ipb, jab (ipc, ipa))   + &
                                  ps(ipc, jab (ipa, ipb)) )
           enddo
        enddo
     enddo
  enddo

#ifdef __PARA
  call     reduce(27, elop_ )
  call poolreduce(27, elop_ )
#endif

  deallocate (chif      )
  deallocate (depsi     )

  !
  ! Calculates the term depending on the third derivative of the
  !                     Exchange-correlation energy
  !
#ifdef __PARA
  if (my_pool_id .ne. 0) goto 100
#endif
  allocate (d2muxc (nrxx))
  allocate (aux3   (nrxx,3))
  do ipa = 1, 3
     call davcio_drho (aux3 (1, ipa), lrdrho, iudrho, ipa, -1)
  enddo

  d2muxc (:) = 0.0_dp
  do ir = 1, nrxx
     rhotot = rho%of_r(ir,1) + rho_core(ir)
     if ( rhotot.gt. 1.d-30 ) d2muxc(ir)= d2mxc( rhotot)
     if ( rhotot.lt.-1.d-30 ) d2muxc(ir)=-d2mxc(-rhotot)
  enddo

  do ipa = 1, 3
     do ipb = 1, 3
        do ipc = 1, 3
           ps3 (ipa, ipb, ipc) = SUM ( DBLE ( d2muxc(:)   * &
                                              aux3(:,ipa) * &
                                              aux3(:,ipb) * &
                                              aux3(:,ipc) ) ) * &
                                 omega / (nr1*nr2*nr3) 
        enddo
     enddo
  enddo

  deallocate (d2muxc )
  deallocate (aux3   )
#ifdef __PARA
  call reduce (27, ps3)
100 continue
  call mp_bcast(ps3, ionode_id)
#endif

  elop_(:,:,:,2) = elop_(:,:,:,1)
  elop_(:,:,:,3) =   ps3(:,:,:)
  elop_(:,:,:,1) = elop_(:,:,:,2) + elop_(:,:,:,3)

  !
  ! Using fac=e2**1.5, calculates the third derivative of the
  !   energy with respect to electric fields.
  !
  ! Using fac=e2**1.5*fpi/omega, calculates the derivative
  !   of the dielectric constants with respect to electric fields.
  ! NB: The result written in output is in Rydberg units, to convert
  !   to pico-meters/Volt you have to multiply per 2.7502
  ! To obtain the static chi^2 multiply by 1/2
  fac = -e2**1.5_dp * fpi / omega
  elop_(:,:,:,:) =  elop_(:,:,:,:) * fac
  !
  ! wr_all =.true. ==> writes separately the two contributions
  !
  wr_all = .true.

  ntm = 1
  if (wr_all ) ntm = 3

  do il = 1, ntm
     !
     ! Symmetrizes the Electro-optic tensor
     !
     call sym_elop(elop_ (1, 1, 1, il), nsym, s, nat, irt)
     !
     ! Transforms from crystal to cartesians axis
     !
     call trntnsr_3 (elop_ (1, 1, 1, il), at, bg, 1)

     if (il.eq.1) then
        write(6,'(/,10x,''    Electro-optic tensor is defined as '')')
        write(6,'(10x  ,''  the derivative of the dielectric tensor '')')
        write(6,'(10x  ,''    with respect to one electric field '')')
        write(6,'(10x  ,''       units are Rydberg a.u. '',/)')
        write(6,'(10x  ,''  to obtain the static chi^2 multiply by 1/2  '',/)')
        write(6,'(10x  ,''  to convert to pm/Volt multiply per 2.7502  '',/)')
        write(6,'(/,10x,''Electro-optic tensor in cartesian axis: '',/)')
        call DCOPY (27, elop_, 1, eloptns, 1)
     else
        write(6,'(/,10x,''Electro-optic tensor: contribution # '',i3,/)') &
           il - 1
     endif

     do ipc = 1, 3
        do ipb = 1, 3
           write(6,'(10x,''('',3f18.9,'' )'')')      &
                   (elop_ (ipa, ipb, ipc, il), ipa = 1, 3)
        enddo
        write(6,'(10x)')
     enddo
  enddo
  call stop_clock('el_opt')
  return
end subroutine el_opt

