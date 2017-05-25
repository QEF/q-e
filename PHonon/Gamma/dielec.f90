!
! Copyright (C) 2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
SUBROUTINE dielec(do_zstar)
  !-----------------------------------------------------------------------
  !
  !      calculates the dielectric tensor and effective charges
  !
  USE constants, ONLY : fpi
  USE cell_base, ONLY : omega
  USE ions_base, ONLY : nat, zv, ityp
  USE mp_global, ONLY : intra_pool_comm
  USE mp,        ONLY : mp_sum
  USE io_files,  ONLY : seqopn
  USE klist,     ONLY : wk, ngk
  USE wvfct,     ONLY: nbnd, npwx
  USE cgcom

  IMPLICIT NONE
  LOGICAL :: do_zstar
  !
  INTEGER :: npw,ibnd,ipol,jpol,na,nu,ik
  CHARACTER(len=7) :: filbar, fildwf
  real(DP) ::  w, weight
  real(DP), ALLOCATABLE ::  work(:,:)
  COMPLEX(DP), ALLOCATABLE :: dpsi2(:,:), dpsi3(:,:)
  LOGICAL :: done
  !
  CALL start_clock('dielec')
  !
  ALLOCATE (dpsi2( npwx, nbnd))
  ALLOCATE (dpsi3( npwx, nbnd))
  ALLOCATE (work( nbnd, 3))
  !
  epsilon0(:,:) = 0.d0
  IF (do_zstar) zstar (:,:,:) = 0.d0
  !  do ik=1,nks
  ik = 1
  npw= ngk(ik)
  weight = wk(ik)
  w = fpi/omega * weight
  !
  !** calculate Effective Charges (<DeltaV*psi(ion)|DeltaPsi(E)>)
  !
  ! read DeltaPsi(E)
  ! pol. 1
  ipol=1
  iudwf=10+ipol
  WRITE(fildwf,'("fildwx",i1)') ipol
  CALL  seqopn (iudwf,fildwf,'unformatted',done)
  READ (iudwf) dpsi
  CLOSE(unit=iudwf)
  ! pol. 2
  ipol=2
  iudwf=10+ipol
  WRITE(fildwf,'("fildwx",i1)') ipol
  CALL  seqopn (iudwf,fildwf,'unformatted',done)
  READ (iudwf) dpsi2
  CLOSE(unit=iudwf)
  ! pol. 3
  ipol=3
  iudwf=10+ipol
  WRITE(fildwf,'("fildwx",i1)') ipol
  CALL  seqopn (iudwf,fildwf,'unformatted',done)
  READ (iudwf) dpsi3
  CLOSE(unit=iudwf)
  !
  IF (.not.do_zstar) GOTO 10
  !
  DO nu = 1,nmodes
     na  = (nu-1)/3+1
     IF (has_equivalent(na)==0) THEN
        !     DeltaV*psi(ion) for mode nu is recalculated
        CALL dvpsi_kb(ik,nu)
        !
        jpol= mod(nu-1,3)+1
        ! work is the real part of <DeltaV*Psi(ion)|DeltaPsi(E)>
        CALL pw_dot('N',npw,nbnd,dvpsi,npwx,dpsi ,npwx,work(1,1))
        CALL pw_dot('N',npw,nbnd,dvpsi,npwx,dpsi2,npwx,work(1,2))
        CALL pw_dot('N',npw,nbnd,dvpsi,npwx,dpsi3,npwx,work(1,3))
        DO ipol = 1,3
           DO ibnd = 1,nbnd
              zstar(ipol,jpol,na) = zstar(ipol,jpol,na) + 2.0d0*weight*work(ibnd,ipol)
           ENDDO
        ENDDO
     ENDIF
  ENDDO
10 CONTINUE
  !** calculate Dielectric Tensor (<DeltaV*psi(E)\DeltaPsi(E)>)
  !
  DO jpol=1,3
     ! read DeltaV*Psi(elec) for polarization jpol
     iubar=jpol
     WRITE(filbar,'("filbar",i1)') iubar
     CALL  seqopn (iubar,filbar,'unformatted',done)
     READ (iubar) dvpsi
     CLOSE(iubar)
     ! now work is the real part of <DeltaV*psi(E)|DeltaPsi(E)>
     CALL pw_dot('N',npw,nbnd,dvpsi,npwx,dpsi ,npwx,work(1,1))
     CALL pw_dot('N',npw,nbnd,dvpsi,npwx,dpsi2,npwx,work(1,2))
     CALL pw_dot('N',npw,nbnd,dvpsi,npwx,dpsi3,npwx,work(1,3))
     DO ipol = 1,3
        DO ibnd = 1,nbnd
           epsilon0(ipol,jpol) = epsilon0(ipol,jpol) + 4.0d0*w*work(ibnd,ipol)
        ENDDO
     ENDDO
  ENDDO
  !     end do
#if defined(__MPI)
  IF (do_zstar) CALL mp_sum( zstar, intra_pool_comm )
  CALL mp_sum( epsilon0, intra_pool_comm )
#endif
  DEALLOCATE(work)
  DEALLOCATE(dpsi3)
  DEALLOCATE(dpsi2)
  !
  ! add the diagonal part
  !
  DO ipol=1,3
     epsilon0(ipol,ipol) = epsilon0(ipol,ipol) + 1.0d0
     IF (do_zstar) THEN
        DO na=1,nat
           zstar(ipol,ipol,na) = zstar(ipol,ipol,na) + zv(ityp(na))
        ENDDO
     ENDIF
  ENDDO
  !
  CALL stop_clock('dielec')
  !
  RETURN
END SUBROUTINE dielec
