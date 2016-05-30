!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
subroutine chi_test (dvscfs, chif, ik, depsi, auxr, auxg)
  !-----------------------------------------------------------------------
  !
  ! This routine is just a debugging tool. Calculates the scalar product
  ! between the chi-wavefunction and Pc DH |psi> in two different ways.
  ! The results of the two procedures should be the same.
  !

  USE kinds, ONLY : DP
  USE wvfct, ONLY : npwx, nbnd
  USE fft_base, ONLY : dffts
  use ramanm,  ONLY : lrd2w, iud2w, jab
  USE units_ph, ONLY : iuwfc, lrdwf, iudwf
  USE buffers, ONLY : get_buffer
  USE qpoint, ONLY : npwq, nksq
  USE eqv, ONLY : dpsi, dvpsi
  USE control_lr, ONLY : nbnd_occ
  USE wavefunctions_module,  ONLY: evc
  implicit none

  integer :: ik
  complex(DP) :: dvscfs(dffts%nnr,3), chif(npwx,nbnd,6), depsi(npwx,nbnd,3), &
            auxr(dffts%nnr), auxg(npwx)

  complex(DP) :: tmp
  complex(DP), EXTERNAL :: zdotc
  complex(DP) , allocatable :: ps1(:,:,:), ps2(:,:,:), &
                     ps3(:,:,:,:), ps4(:,:,:), au2r(:)
  integer :: ip, jp, ib, jb, ipa, ipb, nrec, ir

  allocate (ps1 (nbnd,3,6)      )
  allocate (ps2 (nbnd,3,6)      )
  allocate (ps3 (nbnd,3,nbnd,3) )
  allocate (ps4 (nbnd,3,nbnd)   )
  allocate (au2r (dffts%nnr)        )

  !
  !----------------------------------------------------------
  !
  do ip = 1, 3
     nrec = (ip - 1) * nksq + ik
     call get_buffer (depsi (1, 1, ip), lrdwf, iudwf, nrec)
  enddo

  do jp = 1, 6
     nrec = (jp - 1) * nksq + ik
     call davcio (dpsi, lrd2w, iud2w, nrec, -1)
     do ip = 1, 3
        do ib = 1, nbnd_occ (ik)
              ps1 (ib, ip, jp) =                        &
                    -zdotc (npwq, depsi (1, ib, ip), 1, &
                                  dpsi (1, ib), 1)
        enddo
     enddo
  enddo

  do ip = 1, 3
     do ib = 1, nbnd_occ (ik)
        do jp = 1, 3
           do jb = 1, nbnd_occ (ik)
              ps3(ib, ip, jb, jp) =                              &
        zdotc (npwq, depsi (1, ib, ip), 1, depsi (1, jb, jp), 1)
           enddo
        enddo
     enddo
  enddo

  do ib = 1, nbnd_occ (ik)
     call cft_wave ( ik, evc (1, ib), au2r, +1 )
     do ip = 1, 3
        do ir = 1, dffts%nnr
           auxr (ir) = au2r (ir) * dvscfs (ir, ip)
        end do
        auxg (:) = (0.d0, 0.d0)
        call cft_wave (ik, auxg, auxr, -1 )
        do jb = 1, nbnd_occ (ik)
           ps4 (ib, ip, jb) =                          &
                 zdotc (npwq, auxg, 1, evc (1, jb), 1)
        enddo
     enddo
  enddo

  do ip = 1, 3
     do ib = 1, nbnd_occ (ik)
        do ipa = 1, 3
           do ipb = 1, 3
              tmp = CMPLX(0.d0, 0.d0,kind=DP)
              do jb = 1, nbnd_occ (ik)
                 tmp = tmp +                                 &
                   ps3 (ib, ip, jb, ipa) * ps4 (jb, ipb, ib)
              enddo
              if (ipa.eq.ipb) tmp = tmp * 2.d0
              ps1 (ib, ip, jab (ipa, ipb)) =        &
              ps1 (ib, ip, jab (ipa, ipb)) - tmp
           enddo
        enddo
     enddo
  enddo

  do ip = 1, 3
     do ib = 1, nbnd_occ (ik)
        call cft_wave (ik, depsi (1, ib, ip), au2r, +1 )
        do ipa = 1, 3
           do ir = 1, dffts%nnr
              auxr (ir) = au2r (ir) * dvscfs (ir, ipa)
           enddo
           auxg (:) = (0.d0, 0.d0)
           call cft_wave (ik, auxg, auxr, -1 )
           do ipb = 1, 3
              tmp = zdotc (npwq, auxg, 1, depsi (1, ib, ipb), 1)
              if (ipa.eq.ipb) tmp = tmp * 2.d0
              ps1 (ib, ip, jab (ipa, ipb)) =       &
              ps1 (ib, ip, jab (ipa, ipb)) + tmp
           enddo
        enddo
     enddo
  enddo
!
!----------------------------------------------------------
!
  do ip = 1, 3
     dpsi (:,:) = (0.d0, 0.d0)
     call dvpsi_e (ik, ip)

     do ib = 1, nbnd_occ (ik)
        auxg (:) = (0.d0, 0.d0)
        call daxpy (2 * npwq, -1.d0, dvpsi (1,ib), 1, auxg, 1)
        call cft_wave (ik, evc (1, ib), auxr, +1 )

        do ir = 1, dffts%nnr
           auxr (ir) = auxr (ir) * dvscfs (ir, ip)
        enddo
        call cft_wave (ik, auxg, auxr, -1 )
        do jp = 1, 6
              ps2 (ib, ip, jp) =                             &
                  zdotc (npwq, auxg, 1, chif (1, ib, jp), 1)
        enddo
     enddo
  enddo
!
!----------------------------------------------------------
! If everything is ok, ps1 should be equal to ps2; if not
! there is a problem
!
  do ib = 1, nbnd_occ (ik)
     do jp = 1, 6
!        write(6,9030) ib,jp, (ps1(ib,ip,jp),ip=1,3)
!        write(6,9030) ib,jp, (ps2(ib,ip,jp),ip=1,3)
!        write(6,'(/)')
        write(6,9031) ib,jp, (                               &
             ps1 (ib, ip, jp) / ps2 (ib, ip, jp), ip = 1, 3)
     enddo
  enddo

 9030 format(' bnd:',i5,' ip:',i5,6e12.6)
 9031 format(' bnd:',i5,' ip:',i5,6f12.6)

  deallocate (ps1  )
  deallocate (ps2  )
  deallocate (ps3  )
  deallocate (ps4  )
  deallocate (au2r )

  return
end subroutine chi_test
