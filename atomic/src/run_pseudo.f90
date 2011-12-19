!
! Copyright (C) 2004 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!---------------------------------------------------------------
SUBROUTINE run_pseudo
  !---------------------------------------------------------------
  !
  !     this routine is a driver to a pseudopotential calculation
  !     with the parameters given in input
  !
  !
  USE kinds, ONLY : dp
  USE radial_grids, ONLY : ndmx
  USE ld1_parameters, ONLY : nwfsx
  USE ld1inc, ONLY : enl, lpaw, nlcc, lsd, latt, pawsetup, &
                     nstoaets, grid, nspin, iter, rhos, rhoc, &
                     nwfts, enlts, llts, jjts, iswts, octs, phits, &
                     vxt, enne, vh, vpsloc, file_potscf, beta, tr2,  &
                     eps0, file_recon, deld, vpstot, nbeta, ddd, etots, &
                     paw_energy, iswitch, lgipaw_reconstruction, use_paw_as_gipaw
  USE atomic_paw, ONLY : new_paw_hamiltonian
  IMPLICIT NONE

  INTEGER :: &
       ns, &    ! counter on pseudowavefunctions
       n,  &    ! counter on mesh
       is       ! counter on spin

  real(DP) :: &
       vnew(ndmx,2)   ! the potential

  INTEGER :: &
       iunps, ios

  LOGICAL :: &
       conv       ! if true convergence reached

  real(DP) :: &
       dddnew(nwfsx,nwfsx,2),   & ! the new D coefficients
       vd(2*(ndmx+nwfsx+nwfsx)), & ! Vloc and D in one array for mixing
       vdnew(2*(ndmx+nwfsx+nwfsx)) ! the new vd array
  INTEGER :: &
       nerr                       ! error message

  real(DP), PARAMETER :: thresh=1.e-10_dp
  INTEGER, PARAMETER :: itmax=200
  CHARACTER(len=256) :: nomefile

  !
  !     initial estimate of the eigenvalues
  !
  DO ns=1,nwfts
     enlts(ns)=enl(nstoaets(ns))
  ENDDO
  !
  !    compute an initial estimate of the potential
  !
  CALL guess_initial_wfc()
  IF (.not.lpaw) THEN
     CALL start_potps ( )
  ELSE
     CALL new_paw_hamiltonian (vpstot, ddd, etots,pawsetup, nwfts, &
                               llts, jjts, nspin, iswts, octs, phits, enlts)
     DO is=1,nspin
        vpstot(1:grid%mesh,is)=vpstot(1:grid%mesh,is)-pawsetup%psloc(1:grid%mesh)
     ENDDO
     CALL vdpack (grid%mesh, ndmx, nbeta, nwfsx, nspin, vpstot, ddd, vd, "PACK")
     DO is=1,nspin
        vpstot(1:grid%mesh,is)=vpstot(1:grid%mesh,is)+pawsetup%psloc(1:grid%mesh)
     ENDDO
  ENDIF
  !
  !     iterate to self-consistency
  !
  DO iter=1,itmax
     CALL ascheqps_drv(vpstot, nspin, thresh, .false., nerr)

     IF (.not.lpaw) THEN
        !
        CALL chargeps(rhos,phits,nwfts,llts,jjts,octs,iswts)
        CALL new_potential(ndmx,grid%mesh,grid,0.0_dp,vxt,lsd,&
             nlcc,latt,enne,rhoc,rhos,vh,vnew,1)

        DO is=1,nspin
           vpstot(:,is)=vpstot(:,is)-vpsloc(:)
        ENDDO

        IF (file_potscf/=' ') THEN
           IF (iter<10) THEN
              WRITE(nomefile,'(a,"_",i1)') trim(file_potscf), iter
           ELSEIF(iter<100) THEN
              WRITE(nomefile,'(a,"_",i2)') trim(file_potscf), iter
           ELSEIF(iter<1000) THEN
              WRITE(nomefile,'(a,"_",i3)') trim(file_potscf), iter
           ELSE
              CALL errore('run_pseudo','problem with iteration',1)
           ENDIF
           OPEN(unit=18,file=trim(nomefile), status='unknown', &
                                             err=100, iostat=ios)
100        CALL errore('run_pseudo','opening file' // nomefile,abs(ios))
           IF (lsd==1) THEN
              DO n=1,grid%mesh
                 WRITE(18,'(5e20.12)') grid%r(n),vnew(n,1)-vpstot(n,1), &
                         vnew(n,1), vnew(n,2)-vpstot(n,2), vnew(n,2)
              ENDDO
           ELSE
              DO n=1,grid%mesh
                 WRITE(18,'(3e26.15)') grid%r(n),vnew(n,1)-vpstot(n,1), &
                          vnew(n,1)
              ENDDO
           ENDIF
           CLOSE(18)
        ENDIF

        CALL vpack(grid%mesh,ndmx,nspin,vnew,vpstot,1)
        CALL dmixp(grid%mesh*nspin,vnew,vpstot,beta,tr2,iter,3,eps0,conv,itmax)
        CALL vpack(grid%mesh,ndmx,nspin,vnew,vpstot,-1)

        DO is=1,nspin
           DO n=1,grid%mesh
              vpstot(n,is)=vpstot(n,is)+vpsloc(n)
           ENDDO
        ENDDO
        CALL newd_at
        !
     ELSE
        !
        CALL new_paw_hamiltonian (vnew, dddnew, etots, &
             pawsetup, nwfts, llts, jjts, nspin, iswts, octs, phits, enlts)
        DO is=1,nspin
           vnew(1:grid%mesh,is)=vnew(1:grid%mesh,is)-pawsetup%psloc(1:grid%mesh)
        ENDDO
        CALL vdpack (grid%mesh, ndmx, nbeta, nwfsx, nspin, vnew, dddnew, vdnew, "PACK")
        CALL dmixp((grid%mesh+nbeta*nbeta)*nspin,vdnew,vd,beta,tr2,iter,3,eps0,conv,itmax)
        CALL vdpack (grid%mesh, ndmx, nbeta, nwfsx, nspin, vpstot, ddd, vd, "UNDO")
        DO is=1,nspin
           vpstot(1:grid%mesh,is)=vpstot(1:grid%mesh,is)+pawsetup%psloc(1:grid%mesh)
        ENDDO
        !
     ENDIF

     IF (conv) THEN
        IF (nerr /= 0) THEN
           IF (iswitch==2) THEN
              CALL infomsg ('run_pseudo','BEWARE! Errors in PS-KS equations')
           ELSE
              CALL errore ('run_pseudo','Errors in PS-KS equation', 1)
           ENDIF
        ENDIF
        GOTO 900
     ENDIF
  ENDDO
  CALL infomsg('run_pseudo','Warning: convergence not achieved')
  !
  !    final calculation with all states
  !
900 CONTINUE

  CALL ascheqps_drv(vpstot, nspin, thresh, .true., nerr)

  IF (.not.lpaw) THEN
     CALL elsdps ( )
  ELSE
     CALL new_paw_hamiltonian (vnew, dddnew, etots, pawsetup, nwfts, &
                llts, jjts, nspin, iswts, octs, phits, enlts, paw_energy)
     CALL elsdps_paw()
  ENDIF

  IF ( lgipaw_reconstruction.and.(.not.use_paw_as_gipaw) ) &
       CALL calculate_gipaw_orbitals()

  IF (file_recon/=' ') CALL write_paw_recon ( )
  !
  !    compute logarithmic derivatives
  !
  IF ( deld > 0.0_dp) CALL lderivps ( )
  !
  !    compute expansion in partial waves
  !
  IF ( deld > 0.0_dp) CALL partial_wave_expansion ( )

  RETURN
END SUBROUTINE run_pseudo
