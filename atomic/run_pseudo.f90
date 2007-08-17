!
! Copyright (C) 2004 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!---------------------------------------------------------------
subroutine run_pseudo
  !---------------------------------------------------------------
  !
  !     this routine is a driver to a pseudopotential calculation
  !     with the parameters given in input
  !
  !
  use kinds, only : dp
  use radial_grids, only : ndmx
  use ld1_parameters, only : nwfsx
  use ld1inc, only : enl, lpaw, nlcc, lsd, latt, pawsetup, &
                     nstoaets, grid, nspin, iter, rhos, rhoc, &
                     nwfts, enlts, llts, jjts, iswts, octs, phits, &
                     vxt, enne, vh, vpsloc, file_potscf, beta, tr2,  &
                     eps0, file_recon, deld, vpstot, nbeta, ddd, etots
  use atomic_paw, only : new_paw_hamiltonian, paw2us
  implicit none

  integer :: &
       ns, &    ! counter on pseudowavefunctions
       n,  &    ! counter on mesh
       is, &    ! counter on spin
       nbf      ! number of beta functions

  real(DP) :: &
       vaux(ndmx),  &   ! auxiliary variable
       vnew(ndmx,2)   ! the potential

  integer :: &
       n1,n2,nst,ikl,ind,ios

  logical :: &
       conv       ! if true convergence reached

  real(DP) :: &
       nvalts,                  & ! number of valence electrons for this conf.
       dddnew(nwfsx,nwfsx,2),   & ! the new D coefficients
       ocstart(nwfsx),          & ! guess for the occupations
       vd(2*(ndmx+nwfsx+nwfsx)), & ! Vloc and D in one array for mixing
       vdnew(2*(ndmx+nwfsx+nwfsx)) ! the new vd array
  integer :: &
       iswstart(nwfsx)            ! guess for the starting spins

  real(DP), parameter :: thresh=1.e-10_dp
  integer, parameter :: itmax=200
  character(len=256) :: nomefile

  !
  !     initial estimate of the eigenvalues
  !
  do ns=1,nwfts
     enlts(ns)=enl(nstoaets(ns))
  enddo
  !
  !    compute an initial estimate of the potential
  !
  if (.not.lpaw) then
     call guess_initial_wfc()
     call start_potps ( )
  else
     ! Set starting occupations by rescaling those of the generating configuration
     nvalts=0._dp
     do ns=1,nbeta
        if (octs(ns)>0._dp) nvalts=nvalts+octs(ns)
     end do
     ocstart(1:pawsetup%nwfc) = pawsetup%oc(1:pawsetup%nwfc) * nvalts &
                               / SUM(pawsetup%oc(1:pawsetup%nwfc))
     iswstart=1
     ! Generate the corresponding local and nonlocal potentials
     CALL new_paw_hamiltonian (vpstot, ddd, etots, &
          pawsetup, pawsetup%nwfc, pawsetup%l, 1,iswstart,ocstart, &
          pawsetup%pswfc, pawsetup%enl)
     vpstot(1:grid%mesh,2)=vpstot(1:grid%mesh,1)
     ddd(1:nbeta,1:nbeta,2)=ddd(1:nbeta,1:nbeta,1)
     do is=1,nspin
        vpstot(1:grid%mesh,is)=vpstot(1:grid%mesh,is)-pawsetup%psloc(1:grid%mesh)
     enddo
     call vdpack (grid%mesh, ndmx, nbeta, nwfsx, nspin, vpstot, ddd, vd, "PACK")
     do is=1,nspin
        vpstot(1:grid%mesh,is)=vpstot(1:grid%mesh,is)+pawsetup%psloc(1:grid%mesh)
     enddo
  endif
  !
  !     iterate to self-consistency
  !
  do iter=1,itmax
     call ascheqps_drv(vpstot, nspin, thresh, .false.)

     if (.not.lpaw) then
        !
        call chargeps(rhos,phits,nwfts,llts,jjts,octs,iswts)
        call new_potential(ndmx,grid%mesh,grid,0.0_dp,vxt,lsd,&
             nlcc,latt,enne,rhoc,rhos,vh,vnew,1)

        do is=1,nspin
           vpstot(:,is)=vpstot(:,is)-vpsloc(:)
        enddo

        if (file_potscf.ne.' ') then
           if (iter<10) then
              write(nomefile,'(a,"_",i1)') trim(file_potscf), iter
           elseif(iter<100) then
              write(nomefile,'(a,"_",i2)') trim(file_potscf), iter
           elseif(iter<1000) then
              write(nomefile,'(a,"_",i3)') trim(file_potscf), iter
           else
              call errore('run_pseudo','problem with iteration',1)
           endif
           open(unit=18,file=trim(nomefile), status='unknown', &
                                             err=100, iostat=ios)
100        call errore('run_pseudo','opening file' // nomefile,abs(ios))
           if (lsd==1) then
              do n=1,grid%mesh
                 write(18,'(5e20.12)') grid%r(n),vnew(n,1)-vpstot(n,1), &
                         vnew(n,1), vnew(n,2)-vpstot(n,2), vnew(n,2)
              enddo
           else
              do n=1,grid%mesh
                 write(18,'(3e26.15)') grid%r(n),vnew(n,1)-vpstot(n,1), &
                          vnew(n,1)
              enddo
           endif
           close(18)
        endif

        call vpack(grid%mesh,ndmx,nspin,vnew,vpstot,1)
        call dmixp(grid%mesh*nspin,vnew,vpstot,beta,tr2,iter,3,eps0,conv)
        call vpack(grid%mesh,ndmx,nspin,vnew,vpstot,-1)

        do is=1,nspin
           do n=1,grid%mesh
              vpstot(n,is)=vpstot(n,is)+vpsloc(n)
           enddo
        enddo
        call newd_at
        !
     else
        !
        call new_paw_hamiltonian (vnew, dddnew, etots, &
             pawsetup, nwfts, llts, nspin, iswts, octs, phits, enlts)
        do is=1,nspin
           vnew(1:grid%mesh,is)=vnew(1:grid%mesh,is)-pawsetup%psloc(1:grid%mesh)
        enddo
        call vdpack (grid%mesh, ndmx, nbeta, nwfsx, nspin, vnew, dddnew, vdnew, "PACK")
        call dmixp((grid%mesh+nbeta*nbeta)*nspin,vdnew,vd,beta,tr2,iter,3,eps0,conv)
        call vdpack (grid%mesh, ndmx, nbeta, nwfsx, nspin, vpstot, ddd, vd, "UNDO")
        do is=1,nspin
           vpstot(1:grid%mesh,is)=vpstot(1:grid%mesh,is)+pawsetup%psloc(1:grid%mesh)
        enddo
        !
     endif

!            write(6,*) 'iteration number',iter, eps0
     if (conv) goto 900
  enddo
  call infomsg('run_pseudo','convergence not achieved')

  !
  !    final calculation with all states
  !
900 continue

  call ascheqps_drv(vpstot, nspin, thresh, .true.)

  if (.not.lpaw) then
     call elsdps ( )
  else
     call new_paw_hamiltonian (vnew, dddnew, etots, &
          pawsetup, nwfts, llts, nspin, iswts, octs, phits, enlts)
  endif

  if (file_recon.ne.' ')  call write_paw_recon ( )

  !
  !    compute logarithmic derivatives
  !
  if ( deld > 0.0_dp) call lderivps ( )

  return
end subroutine run_pseudo
