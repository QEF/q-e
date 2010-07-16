!
! Copyright (C) 2004 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!---------------------------------------------------------------
SUBROUTINE scf(ic)
  !---------------------------------------------------------------
  !
  !   this routine performs the atomic self-consistent procedure
  !   self-interaction-correction allowed
  !
  USE kinds, ONLY : dp
  USE radial_grids, ONLY : ndmx
  USE constants, ONLY: e2
  USE ld1inc, ONLY : grid, zed, psi, isic, vpot, vh, vxt, rho, iter, &
                     lsd, rel, latt, enne, beta, nspin, tr2, eps0, &
                     nwf, nn, ll, jj, enl, oc, isw, core_state, frozen_core, &
                     vtau, vsic, vsicnew, vhn1, egc, relpert, noscf
  IMPLICIT NONE

  INTEGER, INTENT(in) :: ic

  LOGICAL:: conv
  INTEGER:: nerr, nstop, n, i, is, id, nin, mch
  real(DP) ::  vnew(ndmx,2), rhoc1(ndmx), ze2
  INTEGER, PARAMETER :: maxter=200
  real(DP), PARAMETER :: thresh=1.0e-10_dp
  !
  !
  ze2 = - zed * e2
  rhoc1=0.0_dp
  id=3
  IF (.not.frozen_core.or.ic==1) psi=0.0_dp
  !!!
  IF (isic /= 0 .and. relpert)  id=1 !
  !!!
  DO iter=1,maxter
     nerr=0
     vnew=vpot
     DO n=1,nwf
        IF (oc(n) >= 0.0_dp) THEN
           IF (ic==1.or..not.frozen_core.or..not.core_state(n)) THEN
              is=isw(n)
              IF (isic /= 0 .and. iter > 1) vnew(:,is)=vpot(:,is)-vsic(:,n)
              IF (rel == 0) THEN
                 CALL ascheq (nn(n),ll(n),enl(n),grid%mesh,grid,vnew(1,is),&
                      ze2,thresh,psi(1,1,n),nstop)
              ELSEIF (rel == 1) THEN
                 CALL lschps (1,zed,grid,nin,mch,nn(n),ll(n),enl(n),&
                             vnew(1,is),vtau,psi(1,1,n),nstop)
                 IF (nstop>0.and.oc(n)<1.e-10_DP) nstop=0
              ELSEIF (rel == 2) THEN
                 CALL dirsol (ndmx,grid%mesh,nn(n),ll(n),jj(n),iter,enl(n), &
                      thresh,grid,psi(1,1,n),vnew(1,is),nstop)
              ELSE
                 CALL errore('scf','relativistic not programmed',1)
              ENDIF
              !      write(6,*) el(n),enl(n)
              ! if (nstop /= 0) write(6,'(4i6)') iter,nn(n),ll(n),nstop
              nerr=nerr+nstop
           ENDIF
        ELSE
           enl(n)=0.0_dp
           psi(:,:,n)=0.0_dp
        ENDIF
     ENDDO
     !
     ! calculate charge density (spherical approximation)
     !
     rho=0.0_dp
     IF (noscf) GOTO 500
     DO n=1,nwf
        DO i=1,grid%mesh
           rho(i,isw(n))=rho(i,isw(n))+oc(n)*(psi(i,1,n)**2+psi(i,2,n)**2)
        ENDDO
     ENDDO
     !
     ! calculate new potential
     !
     CALL new_potential(ndmx,grid%mesh,grid,zed,vxt,&
          lsd,.false.,latt,enne,rhoc1,rho,vh,vnew,1)
     !
     ! calculate SIC correction potential (if present)
     !
     IF (isic /= 0) THEN
        DO n=1,nwf
           IF (oc(n) >= 0.0_dp) THEN
              is=isw(n)
              CALL sic_correction(n,vhn1,vsicnew,egc)
              !
              ! use simple mixing for SIC correction
              !
              vsic(:,n) = (1.0_dp-beta)*vsic(:,n)+beta*vsicnew(:)
           ENDIF
        ENDDO
     ENDIF
     !
     ! mix old and new potential
     !
     CALL vpack(grid%mesh,ndmx,nspin,vnew,vpot,1)
     CALL dmixp(grid%mesh*nspin,vnew,vpot,beta,tr2,iter,id,eps0,conv,maxter)
     CALL vpack(grid%mesh,ndmx,nspin,vnew,vpot,-1)
!        write(6,*) iter, eps0
     !
500  IF (noscf) THEN
        conv=.true.
        eps0=0.0_DP
     ENDIF
     IF (conv) THEN
        IF (nerr /= 0) CALL infomsg ('scf','errors in KS equations')
        GOTO 45
     ENDIF
  ENDDO
  CALL infomsg('scf','warning: convergence not achieved')
45 RETURN

END SUBROUTINE scf
