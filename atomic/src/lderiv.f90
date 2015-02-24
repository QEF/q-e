!
! Copyright (C) 2004-201 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!---------------------------------------------------------------
SUBROUTINE lderiv
  !---------------------------------------------------------------
  !
  !  numerical integration of the radial schroedinger equation
  !  computing logarithmic derivatives for Coulomb potential
  !
  !
  USE kinds,     ONLY : dp
  USE radial_grids, ONLY : ndmx
  USE io_global, ONLY : stdout
  USE mp,        ONLY : mp_bcast
  USE ld1_parameters, ONLY : nwfsx
  USE ld1inc,    ONLY : file_logder, grid, vpot, rel, nspin, nld, zed, &
                        npte, deld, eminld, emaxld, rlderiv

  IMPLICIT NONE

  INTEGER ::        &
       lam,    &   ! the angular momentum
       ikrld,  &   ! index of matching radius
       nc,     &   ! counter on logarithmic derivatives
       nin,    &   ! integer variable for lschps
       is,     &   ! counter on spin
       nstop,  &   ! integer to monitor errors
       ios,    &   ! used for I/O control
       n,ie        ! generic counter

  real(DP) ::           &
       aux(ndmx),         & ! the square of the wavefunction
       aux_dir(ndmx,2),   & ! the square of the wavefunction
       ze2,              & ! the nuclear charge in Ry units
       e,                & ! the eigenvalue
       j,                & ! total angular momentum for log_der
       thrdum = 0.0_dp     ! real variable (not used) for lschps

  real(DP), EXTERNAL :: compute_log

  real(DP), ALLOCATABLE ::        &
       ene(:),        &    ! the energy grid
       dlchi(:, :)         ! the logarithmic derivative

  CHARACTER(len=256) :: flld

  IF (nld == 0 .or. file_logder == ' ') RETURN
  IF (nld > nwfsx) CALL errore('lderiv','nld is too large',1)

  ze2=-zed*2.0_dp

  DO n=1,grid%mesh
     IF (grid%r(n) > rlderiv) GOTO 10
  ENDDO
  CALL errore('lderiv','wrong rlderiv?',1)
10 ikrld = n-1
  WRITE(stdout,'(5x,''Computing logarithmic derivative in'',f10.5)') &
       (grid%r(ikrld)+grid%r(ikrld+1))*0.5_dp

  npte= (emaxld-eminld)/deld + 1
  ALLOCATE ( dlchi(npte, nld) )
  ALLOCATE ( ene(npte) )
  DO ie=1,npte
     ene(ie)= eminld+deld*(ie-1)
  ENDDO

  DO is=1,nspin
     DO nc=1,nld
        IF (rel < 2) THEN
           lam=nc-1
           j=0.0_dp
        ELSE
           lam=nc/2
           IF (mod(nc,2)==0) j=lam-0.5_dp
           IF (mod(nc,2)==1) j=lam+0.5_dp
        ENDIF
        DO ie=1,npte
           e=ene(ie)
           !
           !    integrate outward up to ikrld+1
           !
           IF (rel == 1) THEN
              nin = ikrld+5
              CALL lschps (3, zed, thrdum, grid, nin, 1, lam, e, &
                   vpot(1,is), aux, nstop )
           ELSEIF (rel == 2) THEN
              CALL dir_outward(ndmx,ikrld+5,lam,j,e,grid%dx,&
                   aux_dir,grid%r,grid%rab,vpot(1,is))
              aux(:)=aux_dir(:,1)
           ELSE
              CALL intref(lam,e,ikrld+5,grid,vpot(1,is),ze2,aux)
           ENDIF
           !
           !    compute the logarithmic derivative and save in dlchi
           !
           dlchi(ie, nc) = compute_log(aux(ikrld-3),grid%r(ikrld),grid%dx)
        ENDDO
     ENDDO

     IF (nspin == 2 .and. is == 1) THEN
        flld = trim(file_logder)//'up'
     ELSEIF (nspin == 2 .and. is == 2) THEN
        flld = trim(file_logder)//'dw'
     ELSE
        flld = trim(file_logder)
     ENDIF

     CALL write_efun(flld,dlchi,ene,npte,nld)
     !
  ENDDO

  DEALLOCATE (ene)
  DEALLOCATE (dlchi)
  RETURN
END SUBROUTINE lderiv

