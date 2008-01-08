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
  USE ener,      ONLY : ef
  USE lsda_mod,  ONLY : nspin, current_spin
  USE scf,       ONLY : rho, vltot, v, rho_core, rhog_core
  USE gvect
  USE cell_base, ONLY : omega, alat
  USE fft_base,  ONLY : grid_gather
  USE mp,        ONLY : mp_bcast

  IMPLICIT NONE

  REAL(DP) :: wmean1, wmean2, meancharge, wx1, wx2, wxm, vx, vc, ex, &
                   ec, rhox, rs, vcca, wf, etxc, vtxc
  INTEGER :: n1, n2, ni, nmean, nspin0
  LOGICAL :: exst
  REAL(DP), ALLOCATABLE :: raux1 (:), vaux1 (:), vaux2(:), aux (:)
  REAL(DP), ALLOCATABLE :: vxc(:,:)
  ! auxiliary vectors for charge and potential

  ALLOCATE (raux1( nrx1 * nrx2 * nrx3))    
  ALLOCATE (vaux1( nrx1 * nrx2 * nrx3))    
  ALLOCATE (vaux2( nrx1 * nrx2 * nrx3))    

  nspin0=nspin
  if (nspin==4) nspin0=1

  ALLOCATE (vxc(nrxx,nspin))
  CALL v_xc (rho, rho_core, rhog_core, etxc, vtxc, vxc)

  IF ( ionode ) THEN
     !
     CALL seqopn (17, 'workf', 'formatted', exst)
     CALL seqopn (19, 'charge', 'formatted', exst)
     !
  END IF

  wf = 0.d0

  DO current_spin=1,nspin0

#ifdef __PARA
     ALLOCATE (aux  ( nrxx))    
     aux(:) = rho%of_r(:,current_spin) + rho_core(:)/nspin0
     CALL grid_gather (aux, raux1)
#else
     raux1(1:nrxx) = rho%of_r(1:nrxx,current_spin) + rho_core(1:nrxx)/nspin0
#endif
     !
#ifdef __PARA
     aux(:) = vltot(:) + v%of_r(:,current_spin)
     CALL grid_gather (aux, vaux1)
     aux(:) = aux(:) - vxc(:,current_spin)
     CALL grid_gather (aux, vaux2)
#else
     vaux1(1:nrxx) = vltot(1:nrxx) + v%of_r(1:nrxx,current_spin)
     vaux2(1:nrxx) = vaux1(1:nrxx) -vxc(1:nrxx,current_spin)
#endif
     !
#ifdef __PARA
     DEALLOCATE(aux)
#endif
     IF ( ionode ) THEN
        !
        IF (nspin == 2) THEN
           IF (current_spin.EQ.1) THEN
              WRITE(17,*) " SPIN UP "
              WRITE(19,*) " SPIN UP "
           ELSE
              WRITE(17,*) " SPIN DOWN "
              WRITE(19,*) " SPIN DOWN "
           END IF
        ENDIF
        DO nmean = 1, nr3
           wmean1 = 0.d0
           wmean2 = 0.d0
           meancharge = 0.d0
           wx1 = 0.d0
           wx2 = 0.d0
           wxm = 0.d0
           DO n2 = 1, nr2
              DO n1 = 1, nr1
                 ni = n1 + (n2 - 1) * nrx1 + (nmean - 1) * nrx1 * nrx2
                 meancharge = meancharge+raux1 (ni)
                 wxm = wxm + raux1 (ni) **2
                 wmean1 = wmean1 + vaux1 (ni)
                 wx1 = wx1 + vaux1 (ni) **2
                 wmean2 = wmean2 + vaux2 (ni) 
                 wx2 = wx2 + vaux2 (ni) **2
              ENDDO
           ENDDO
           wmean1 = wmean1 / DBLE (nr1 * nr2)
           wmean2 = wmean2 / DBLE (nr1 * nr2)
           meancharge = meancharge / DBLE (nr1 * nr2)
           wx1 = dsqrt (wx1 / DBLE (nr1 * nr2) - wmean1 * wmean1)
           wx2 = dsqrt (wx2 / DBLE (nr1 * nr2) - wmean2 * wmean2)
           wxm = dsqrt (wxm / DBLE (nr1 * nr2) - meancharge**2)
           IF (nmean.EQ. (nr3 + 1) / 2) THEN
              wf = wf + (wmean2 - ef)
              IF (nspin == 2) THEN
                 IF (current_spin.EQ.1) THEN
                    WRITE( stdout,*) " SPIN UP "
                 ELSE
                    WRITE( stdout,*) " SPIN DOWN "
                 END IF
              ENDIF
              WRITE( stdout, 9130) rytoev * (wmean1 - ef), wx1 * rytoev, &
                   rytoev * (wmean2 - ef), wx2 * rytoev
           END IF
           WRITE (17, * ) nmean, (wmean1 - ef) * rytoev, wx1 * rytoev, &
                (wmean2 - ef) * rytoev, wx2 * rytoev
           WRITE (19, * ) nmean, meancharge, wxm
        ENDDO
        !
     ENDIF
  !
  END DO
  wf = wf / nspin0
  !
  CALL mp_bcast( wf, ionode_id )

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
