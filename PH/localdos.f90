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
  use becmod
  use parameters, only : DP
  use phcom
  implicit none

  complex(kind=DP) :: ldos (nrxx, nspin), ldoss (nrxxs, nspin)
  ! output: the local density of states at
  ! output: the local density of states at
  !         without augmentation
  real(kind=DP) :: dos_ef
  ! output: the density of states at Ef
  !
  !    local variables for Ultrasoft PP's
  !

  integer :: ikb, jkb, ijkb0, ih, jh, na, ijh, nt
  ! counter on beta functions
  ! counter on beta functions
  ! auxiliary variable for ijkb0
  ! counter on solid beta functions
  ! counter on solid beta functions
  ! counter on atoms
  ! counter on composite beta functions
  ! counter on atomic types
  real(kind=DP), allocatable :: becsum1 (:,:,:)
  ! the becsum1 terms
  !
  ! local variables
  !

  real(kind=DP) :: weight, w1, wdelta, w0gauss
  ! kpoint wheight
  ! weight
  ! delta function weight
  ! delta function

  real(kind=DP) :: check


  integer :: ik, is, ig, ibnd, j
  ! kpoint index
  ! counter on spin polarizations
  ! g vector index
  ! band index
  ! fft index

  integer :: ios
  ! status flag for i/o

  !
  !  initialize ldos and dos_ef
  !
  call start_clock ('localdos')
  allocate (becsum1( (nhm * (nhm + 1)) / 2, nat, nspin))    

  call setv ( (nhm * (nhm + 1) ) / 2 * nat * nspin, 0.d0, becsum1, 1)
  call setv (2 * nrxx * nspin, 0.d0, ldos, 1)
  call setv (2 * nrxxs * nspin, 0.d0, ldoss, 1)
  dos_ef = 0.d0
  !
  !  loop over kpoints
  !
  if (nksq.gt.1) rewind (unit = iunigk)
  do ik = 1, nksq
     if (lsda) current_spin = isk (ik)
     if (nksq.gt.1) then
        read (iunigk, err = 100, iostat = ios) npw, igk
100     call errore ('solve_linter', 'reading igk', abs (ios) )
     endif
     weight = wk (ik)
     !
     ! unperturbed wfs in reciprocal space read from unit iuwfc
     !

     if (nksq.gt.1) call davcio (evc, lrwfc, iuwfc, ik, - 1)
     call init_us_2 (npw, igk, xk (1, ik), vkb)

     call ccalbec (nkb, npwx, npw, nbnd, becp, vkb, evc)
     do ibnd = 1, nbnd_occ (ik)
        wdelta = w0gauss ( (ef-et(ibnd,ik)) / degauss, ngauss) / degauss
        !
        ! unperturbed wf from reciprocal to real space
        !
        call setv (2 * nrxxs, 0.d0, psic, 1)
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
                 if (ityp (na) .eq.nt) then
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
                 if (ityp (na) .eq.nt) ijkb0 = ijkb0 + nh (nt)
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
     call ZCOPY (nrxx * nspin, ldoss, 1, ldos, 1)
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
  !      write (*,*) ' check ', check, dos_ef
  !check
  !
  deallocate(becsum1)
  call stop_clock ('localdos')
  return
end subroutine localdos
