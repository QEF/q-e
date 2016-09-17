!
! Copyright (C) 2001-2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
SUBROUTINE work_function (wf)
  !
  ! Print out the workfunction, calculated as the difference between the
  ! potential energy and the fermi energy.
  ! Written for supercells with the main axis along z.
  !
  USE constants, ONLY : rytoev, e2
  USE io_global, ONLY : stdout, ionode, ionode_id
  USE io_files,  ONLY : seqopn
  USE ener,      ONLY : ef
  USE lsda_mod,  ONLY : nspin, current_spin
  USE scf,       ONLY : rho, vltot, v, rho_core, rhog_core
  USE gvect
  USE cell_base, ONLY : omega, alat
  USE fft_base,  ONLY : dfftp
  USE scatter_mod, ONLY : gather_grid
  USE mp,        ONLY : mp_bcast
  USE mp_world,  ONLY : world_comm

  IMPLICIT NONE

  REAL(DP) :: wmean1, wmean2, meancharge, wx1, wx2, wxm, wf, etxc, vtxc
  INTEGER :: n1, n2, ni, nmean, nspin0
  LOGICAL :: exst
  REAL(DP), ALLOCATABLE :: raux1 (:), vaux1 (:), vaux2(:), aux (:)
  REAL(DP), ALLOCATABLE :: vxc(:,:)
  ! auxiliary vectors for charge and potential

  ALLOCATE (raux1( dfftp%nr1x * dfftp%nr2x * dfftp%nr3x))
  ALLOCATE (vaux1( dfftp%nr1x * dfftp%nr2x * dfftp%nr3x))
  ALLOCATE (vaux2( dfftp%nr1x * dfftp%nr2x * dfftp%nr3x))

  nspin0=nspin
  IF (nspin==4) nspin0=1

  ALLOCATE (vxc(dfftp%nnr,nspin))
  CALL v_xc (rho, rho_core, rhog_core, etxc, vtxc, vxc)

  IF ( ionode ) THEN
     !
     CALL seqopn (17, 'workf', 'formatted', exst)
     CALL seqopn (19, 'charge', 'formatted', exst)
     !
  ENDIF

  wf = 0.d0

  DO current_spin=1,nspin0

#if defined(__MPI)
     ALLOCATE (aux  ( dfftp%nnr))
     aux(:) = rho%of_r(:,current_spin) + rho_core(:)/nspin0
     CALL gather_grid (dfftp, aux, raux1)
#else
     raux1(1:dfftp%nnr) = rho%of_r(1:dfftp%nnr,current_spin) + rho_core(1:dfftp%nnr)/nspin0
#endif
     !
#if defined(__MPI)
     aux(:) = vltot(:) + v%of_r(:,current_spin)
     CALL gather_grid (dfftp, aux, vaux1)
     aux(:) = aux(:) - vxc(:,current_spin)
     CALL gather_grid (dfftp, aux, vaux2)
#else
     vaux1(1:dfftp%nnr) = vltot(1:dfftp%nnr) + v%of_r(1:dfftp%nnr,current_spin)
     vaux2(1:dfftp%nnr) = vaux1(1:dfftp%nnr) -vxc(1:dfftp%nnr,current_spin)
#endif
     !
#if defined(__MPI)
     DEALLOCATE(aux)
#endif
     IF ( ionode ) THEN
        !
        IF (nspin == 2) THEN
           IF (current_spin==1) THEN
              WRITE(17,*) " SPIN UP "
              WRITE(19,*) " SPIN UP "
           ELSE
              WRITE(17,*) " SPIN DOWN "
              WRITE(19,*) " SPIN DOWN "
           ENDIF
        ENDIF
        DO nmean = 1, dfftp%nr3
           wmean1 = 0.d0
           wmean2 = 0.d0
           meancharge = 0.d0
           wx1 = 0.d0
           wx2 = 0.d0
           wxm = 0.d0
           DO n2 = 1, dfftp%nr2
              DO n1 = 1, dfftp%nr1
                 ni = n1 + (n2 - 1) * dfftp%nr1x + (nmean - 1) * dfftp%nr1x * dfftp%nr2x
                 meancharge = meancharge+raux1 (ni)
                 wxm = wxm + raux1 (ni) **2
                 wmean1 = wmean1 + vaux1 (ni)
                 wx1 = wx1 + vaux1 (ni) **2
                 wmean2 = wmean2 + vaux2 (ni)
                 wx2 = wx2 + vaux2 (ni) **2
              ENDDO
           ENDDO
           wmean1 = wmean1 / dble (dfftp%nr1 * dfftp%nr2)
           wmean2 = wmean2 / dble (dfftp%nr1 * dfftp%nr2)
           meancharge = meancharge / dble (dfftp%nr1 * dfftp%nr2)
           wx1 = dsqrt (wx1 / dble (dfftp%nr1 * dfftp%nr2) - wmean1 * wmean1)
           wx2 = dsqrt (wx2 / dble (dfftp%nr1 * dfftp%nr2) - wmean2 * wmean2)
           wxm = dsqrt (wxm / dble (dfftp%nr1 * dfftp%nr2) - meancharge**2)
           IF (nmean== (dfftp%nr3 + 1) / 2) THEN
              wf = wf + (wmean2 - ef)
              IF (nspin == 2) THEN
                 IF (current_spin==1) THEN
                    WRITE( stdout,*) " SPIN UP "
                 ELSE
                    WRITE( stdout,*) " SPIN DOWN "
                 ENDIF
              ENDIF
              WRITE( stdout, 9130) rytoev * (wmean1 - ef), wx1 * rytoev, &
                   rytoev * (wmean2 - ef), wx2 * rytoev
           ENDIF
           WRITE (17, * ) nmean, (wmean1 - ef) * rytoev, wx1 * rytoev, &
                (wmean2 - ef) * rytoev, wx2 * rytoev
           WRITE (19, * ) nmean, meancharge, wxm
        ENDDO
        !
     ENDIF
  !
  ENDDO
  wf = wf / nspin0
  !
  CALL mp_bcast( wf, ionode_id, world_comm )

  WRITE( stdout, '(/5x,"Work function written on file workf")')
  WRITE( stdout, '( 5x,"Planar mean charge written on file charge")')

9130 FORMAT (/'     workfunction     = ',f10.4,' +- ',f10.4,' eV', &
    &        /'     without exchcorr = ',f10.4,' +- ',f10.4,' eV')

  CLOSE (17)
  CLOSE (19)

  DEALLOCATE(raux1)
  DEALLOCATE(vaux1)
  DEALLOCATE(vaux2)
  DEALLOCATE(vxc)

  RETURN

END SUBROUTINE work_function
