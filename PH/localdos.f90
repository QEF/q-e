!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
subroutine localdos (ldos, ldoss, dos_ef)
  !-----------------------------------------------------------------------
  !
  !    This routine compute the local and total density of state at Ef
  !
  !    Note: this routine use psic as auxiliary variable. it should alread
  !          be defined
  !
  !    NB: this routine works only with gamma
  !
#include "machine.h"

  use pwcom
  USE wavefunctions_module,  ONLY: evc, psic
  USE kinds, only : DP
  use phcom
  USE io_files, ONLY: iunigk
  implicit none

  complex(kind=DP) :: ldos (nrxx, nspin), ldoss (nrxxs, nspin)
  ! output: the local density of states at Ef
  ! output: the local density of states at Ef without augmentation
  real(kind=DP) :: dos_ef
  ! output: the density of states at Ef
  !
  !    local variables for Ultrasoft PP's
  !
  integer :: ikb, jkb, ijkb0, ih, jh, na, ijh, nt
  ! counters
  real(kind=DP), allocatable :: becsum1 (:,:,:)
  complex(kind=DP), allocatable :: becp(:,:)
  !
  ! local variables
  !
  real(kind=DP) :: weight, w1, wdelta
  ! weights
  real(kind=DP), external :: w0gauss
  !
  integer :: ik, is, ig, ibnd, j
  ! counters
  integer :: ios
  ! status flag for i/o
  !
  !  initialize ldos and dos_ef
  !
  call start_clock ('localdos')
  allocate (becsum1( (nhm * (nhm + 1)) / 2, nat, nspin), becp(nkb,nbnd) )  

  becsum1 (:,:,:) = 0.d0
  ldos (:,:) = (0d0, 0.0d0)
  ldoss(:,:) = (0d0, 0.0d0)
  dos_ef = 0.d0
  !
  !  loop over kpoints
  !
  if (nksq > 1) rewind (unit = iunigk)
  do ik = 1, nksq
     if (lsda) current_spin = isk (ik)
     if (nksq > 1) then
        read (iunigk, err = 100, iostat = ios) npw, igk
100     call errore ('solve_linter', 'reading igk', abs (ios) )
     endif
     weight = wk (ik)
     !
     ! unperturbed wfs in reciprocal space read from unit iuwfc
     !
     if (nksq > 1) call davcio (evc, lrwfc, iuwfc, ik, - 1)
     call init_us_2 (npw, igk, xk (1, ik), vkb)
     !
     call ccalbec (nkb, npwx, npw, nbnd, becp, vkb, evc)
     do ibnd = 1, nbnd_occ (ik)
        wdelta = w0gauss ( (ef-et(ibnd,ik)) / degauss, ngauss) / degauss
        !
        ! unperturbed wf from reciprocal to real space
        !
        psic (:) = (0.d0, 0.d0)
        do ig = 1, npw
           psic (nls (igk (ig) ) ) = evc (ig, ibnd)
        enddo
        call cft3s (psic, nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, + 1)
        w1 = weight * wdelta / omega
        do j = 1, nrxxs
           ldoss (j, current_spin) = ldoss (j, current_spin) + &
                 w1 * (DREAL ( psic (j) ) **2 + DIMAG (psic (j) ) **2)
        enddo
        !
        !    If we have a US pseudopotential we compute here the sumbec term
        !
        w1 = weight * wdelta
        ijkb0 = 0
        do nt = 1, ntyp
           if (tvanp (nt) ) then
              do na = 1, nat
                 if (ityp (na) == nt) then
                    ijh = 1
                    do ih = 1, nh (nt)
                       ikb = ijkb0 + ih
                       becsum1 (ijh, na, current_spin) = &
                            becsum1 (ijh, na, current_spin) + w1 * &
                            DREAL (conjg(becp(ikb,ibnd))*becp(ikb,ibnd) )
                       ijh = ijh + 1
                       do jh = ih + 1, nh (nt)
                          jkb = ijkb0 + jh
                          becsum1 (ijh, na, current_spin) = &
                               becsum1 (ijh, na, current_spin) + w1 * 2.d0 * &
                               DREAL(conjg(becp(ikb,ibnd))*becp(jkb,ibnd) )
                          ijh = ijh + 1
                       enddo
                    enddo
                    ijkb0 = ijkb0 + nh (nt)
                 endif
              enddo
           else
              do na = 1, nat
                 if (ityp (na) == nt) ijkb0 = ijkb0 + nh (nt)
              enddo
           endif
        enddo
        dos_ef = dos_ef + weight * wdelta
     enddo

  enddo
  if (doublegrid) then
     do is = 1, nspin
        call cinterpolate (ldos (1, is), ldoss (1, is), 1)
     enddo
  else
     ldos (:,:) = ldoss (:,:) 
  endif

  call addusldos (ldos, becsum1)
#ifdef __PARA
  !
  !    Collects partial sums on k-points from all pools
  !
  call poolreduce (2 * nrxxs * nspin, ldoss)
  call poolreduce (2 * nrxx * nspin, ldos)
  call poolreduce (1, dos_ef)
#endif
  !check
  !      check =0.d0
  !      do is=1,nspin
  !         call cft3(ldos(1,is),nr1,nr2,nr3,nrx1,nrx2,nrx3,-1)
  !         check = check + omega*DREAL(ldos(nl(1),is))
  !         call cft3(ldos(1,is),nr1,nr2,nr3,nrx1,nrx2,nrx3,+1)
  !      end do
  !      WRITE( stdout,*) ' check ', check, dos_ef
  !check
  !
  deallocate(becsum1, becp)
  call stop_clock ('localdos')
  return
end subroutine localdos
