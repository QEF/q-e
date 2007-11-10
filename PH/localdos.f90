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
#include "f_defs.h"
  !
  USE ions_base, ONLY : nat, ityp, ntyp => nsp
  use pwcom
  USE noncollin_module, ONLY : noncolin, npol
  USE wavefunctions_module,  ONLY: evc, psic, psic_nc
  USE kinds, only : DP
  USE uspp_param, ONLY: upf, nh, nhm
  use phcom
  USE io_files, ONLY: iunigk
  implicit none

  complex(DP) :: ldos (nrxx, nspin), ldoss (nrxxs, nspin)
  ! output: the local density of states at Ef
  ! output: the local density of states at Ef without augmentation
  real(DP) :: dos_ef
  ! output: the density of states at Ef
  !
  !    local variables for Ultrasoft PP's
  !
  integer :: ikb, jkb, ijkb0, ih, jh, na, ijh, nt
  ! counters
  real(DP), allocatable :: becsum1 (:,:,:)
  complex(DP), allocatable :: becp(:,:), becp_nc(:,:,:), becsum1_nc(:,:,:,:)
  !
  ! local variables
  !
  real(DP) :: weight, w1, wdelta
  ! weights
  real(DP), external :: w0gauss
  !
  integer :: ik, is, ig, ibnd, j, is1, is2
  ! counters
  integer :: ios
  ! status flag for i/o
  !
  !  initialize ldos and dos_ef
  !
  call start_clock ('localdos')
  allocate (becsum1( (nhm * (nhm + 1)) / 2, nat, nspin))
  IF (noncolin) THEN
     allocate (becp_nc(nkb,npol,nbnd) )  
     allocate (becsum1_nc( (nhm * (nhm + 1)) / 2, nat, npol, npol))
     becsum1_nc=(0.d0,0.d0)
  ELSE
     allocate (becp(nkb,nbnd) )  
  ENDIF

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
     IF (noncolin) THEN
        call ccalbec_nc (nkb, npwx, npw, npol, nbnd, becp_nc, vkb, evc)
     ELSE
        call ccalbec (nkb, npwx, npw, nbnd, becp, vkb, evc)
     END IF
     do ibnd = 1, nbnd_occ (ik)
        wdelta = w0gauss ( (ef-et(ibnd,ik)) / degauss, ngauss) / degauss
        w1 = weight * wdelta / omega
        !
        ! unperturbed wf from reciprocal to real space
        !
        IF (noncolin) THEN
           psic_nc = (0.d0, 0.d0)
           do ig = 1, npw
              psic_nc (nls (igk (ig)), 1 ) = evc (ig, ibnd)
              psic_nc (nls (igk (ig)), 2 ) = evc (ig+npwx, ibnd)
           enddo
           call cft3s (psic_nc, nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, + 1)
           call cft3s (psic_nc(1,2), nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, + 1)
           do j = 1, nrxxs
              ldoss (j, current_spin) = ldoss (j, current_spin) + &
                    w1 * ( DBLE(psic_nc(j,1))**2+AIMAG(psic_nc(j,1))**2 + &
                           DBLE(psic_nc(j,2))**2+AIMAG(psic_nc(j,2))**2)
           enddo
        ELSE
           psic (:) = (0.d0, 0.d0)
           do ig = 1, npw
              psic (nls (igk (ig) ) ) = evc (ig, ibnd)
           enddo
           call cft3s (psic, nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, + 1)
           do j = 1, nrxxs
              ldoss (j, current_spin) = ldoss (j, current_spin) + &
                    w1 * ( DBLE ( psic (j) ) **2 + AIMAG (psic (j) ) **2)
           enddo
        END IF
        !
        !    If we have a US pseudopotential we compute here the becsum term
        !
        w1 = weight * wdelta
        ijkb0 = 0
        do nt = 1, ntyp
           if (upf(nt)%tvanp ) then
              do na = 1, nat
                 if (ityp (na) == nt) then
                    ijh = 1
                    do ih = 1, nh (nt)
                       ikb = ijkb0 + ih
                       IF (noncolin) THEN
                          DO is1=1,npol
                             DO is2=1,npol
                                becsum1_nc (ijh, na, is1, is2) = &
                                becsum1_nc (ijh, na, is1, is2) + w1 * &
                                 (CONJG(becp_nc(ikb,is1,ibnd))* &
                                        becp_nc(ikb,is2,ibnd))
                             END DO
                          END DO
                       ELSE
                          becsum1 (ijh, na, current_spin) = &
                            becsum1 (ijh, na, current_spin) + w1 * &
                             DBLE (CONJG(becp(ikb,ibnd))*becp(ikb,ibnd) )
                       ENDIF
                       ijh = ijh + 1
                       do jh = ih + 1, nh (nt)
                          jkb = ijkb0 + jh
                          IF (noncolin) THEN
                             DO is1=1,npol
                                DO is2=1,npol
                                   becsum1_nc(ijh,na,is1,is2) = &
                                      becsum1_nc(ijh,na,is1,is2) + w1* &
                                      (CONJG(becp_nc(ikb,is1,ibnd))* &
                                             becp_nc(jkb,is2,ibnd) )
                                END DO
                             END DO
                          ELSE
                             becsum1 (ijh, na, current_spin) = &
                               becsum1 (ijh, na, current_spin) + w1 * 2.d0 * &
                                DBLE(CONJG(becp(ikb,ibnd))*becp(jkb,ibnd) )
                          END IF
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

  IF (noncolin.and.okvan) THEN
     DO nt = 1, ntyp
        IF ( upf(nt)%tvanp ) THEN
           DO na = 1, nat
              IF (ityp(na)==nt) THEN
                 IF (so(nt)) THEN
                    CALL transform_becsum_so(becsum1_nc,becsum1,na)
                 ELSE
                    CALL transform_becsum_nc(becsum1_nc,becsum1,na)
                 END IF
              END IF
           END DO
        END IF
     END DO
  END IF

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
  !         check = check + omega* DBLE(ldos(nl(1),is))
  !         call cft3(ldos(1,is),nr1,nr2,nr3,nrx1,nrx2,nrx3,+1)
  !      end do
  !      WRITE( stdout,*) ' check ', check, dos_ef
  !check
  !
  deallocate(becsum1)
  IF (noncolin) THEN
     deallocate(becp_nc)
     deallocate(becsum1_nc)
  ELSE
     deallocate(becp)
  ENDIF
  call stop_clock ('localdos')
  return
end subroutine localdos
