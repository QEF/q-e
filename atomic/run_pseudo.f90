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
  use ld1inc
  use atomic_paw, only : new_paw_hamiltonian, paw2us
  implicit none

  integer :: &
       ns, &    ! counter on pseudowavefunctions
       n,  &    ! counter on mesh
       is, &    ! counter on spin
       nbf      ! number of beta functions

  real(DP) :: &
       vaux(ndm),  &   ! auxiliary variable
       vnew(ndm,2)   ! the potential

  integer :: &
       n1,n2,nst,ikl,ind

  logical :: &
       conv       ! if true convergence reached

  real(DP) :: &
       nvalts,                  & ! number of valence electrons for this conf.
       dddnew(nwfsx,nwfsx,2),   & ! the new D coefficients
       ocstart(nwfsx),          & ! guess for the occupations
       vd(2*(ndm+nwfsx+nwfsx)), & ! Vloc and D in one array for mixing
       vdnew(2*(ndm+nwfsx+nwfsx)) ! the new vd array
  integer :: &
       iswstart(nwfsx)            ! guess for the starting spins

  real(DP), parameter :: thresh=1.e-10_dp
  integer, parameter :: itmax=200

  !
  !     initial estimate of the eigenvalues
  !
  do ns=1,nwfts
     enls(ns)=enl(nstoae(ns))
  enddo
  if (pseudotype.eq.1) then
     nbf=0
  else
     nbf=nbeta
  endif
  !
  !    compute an initial estimate of the potential
  !
  if (.not.lpaw) then
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
     vpstot(1:mesh,2)=vpstot(1:mesh,1)
     ddd(1:nbeta,1:nbeta,2)=ddd(1:nbeta,1:nbeta,1)
     do is=1,nspin
        vpstot(1:mesh,is)=vpstot(1:mesh,is)-pawsetup%psloc(1:mesh)
     enddo
     call vdpack (mesh, ndm, nbeta, nwfsx, nspin, vpstot, ddd, vd, "PACK")
     do is=1,nspin
        vpstot(1:mesh,is)=vpstot(1:mesh,is)+pawsetup%psloc(1:mesh)
     enddo
  endif
  !
  !     iterate to self-consistency
  !
  do iter=1,itmax
     do ns=1,nwfts
        if (octs(ns).gt.0.0_dp) then
           is=iswts(ns)
           if (pseudotype == 1) then
              if ( rel < 2 .or. llts(ns) == 0 .or. &
                   abs(jjts(ns)-llts(ns)+0.5_dp) < 0.001_dp) then
                 ind=1
              else if ( rel==2 .and. llts(ns)>0 .and. &
                   abs(jjts(ns)-llts(ns)-0.5_dp) < 0.001_dp) then
                 ind=2
              else
                 call errore('run-pseudo','spin-orbit?!?',3)
              endif
              do n=1,mesh
                 vaux(n)=vpstot(n,is)+vnl (n,llts(ns),ind)
              enddo
           else
              do n=1,mesh
                 vaux(n)=vpstot(n,is)
              enddo
           endif
           !
           call ascheqps(nnts(ns), llts(ns), jjts(ns), enls(ns),  &
                mesh, ndm, dx, r, r2, sqr, vaux(1), thresh, phis(1,ns), &
                betas, ddd(1,1,is), qq, nbf, nwfsx, lls, jjs, ikk)
           ! write(6,*) 'run_pseu',ns, nwfts, nnts(ns),llts(ns), jjts(ns),enls(ns)  
        endif
     enddo

     call normalize ( )

     if (.not.lpaw) then
        !
        call chargeps(nwfts,llts,jjts,octs,iswts)
        call new_potential(ndm,mesh,r,r2,sqr,dx,0.0_dp,vxt,lsd, &
             nlcc,latt,enne,rhoc,rhos,vh,vnew)

        do is=1,nspin
           vpstot(:,is)=vpstot(:,is)-vpsloc(:)
        enddo

        !         do n=1,mesh
        !            write(6,'(5f15.4)') r(n),vnew(n,1)-vpstot(n,1),
        !     +               vnew(n,1), vnew(n,2)-vpstot(n,2), vnew(n,2)
        !         enddo
        call vpack(mesh,ndm,nspin,vnew,vpstot,1)
        call dmixp(mesh*nspin,vnew,vpstot,beta,tr2,iter,3,eps0,conv)
        call vpack(mesh,ndm,nspin,vnew,vpstot,-1)

        do is=1,nspin
           do n=1,mesh
              vpstot(n,is)=vpstot(n,is)+vpsloc(n)
           enddo
        enddo
        call newd_at
        !
     else
        !
        call new_paw_hamiltonian (vnew, dddnew, etots, &
             pawsetup, nwfts, llts, nspin, iswts, octs, phis, enls)
        do is=1,nspin
           vnew(1:mesh,is)=vnew(1:mesh,is)-pawsetup%psloc(1:mesh)
        enddo
        call vdpack (mesh, ndm, nbeta, nwfsx, nspin, vnew, dddnew, vdnew, "PACK")
        call dmixp((mesh+nbeta*nbeta)*nspin,vdnew,vd,beta,tr2,iter,3,eps0,conv)
        call vdpack (mesh, ndm, nbeta, nwfsx, nspin, vpstot, ddd, vd, "UNDO")
        do is=1,nspin
           vpstot(1:mesh,is)=vpstot(1:mesh,is)+pawsetup%psloc(1:mesh)
        enddo
        !
     endif

     !         write(6,*) 'iteration number',iter, eps0
     if (conv) goto 900
  enddo
  call infomsg('run_pseudo','convergence not achieved', -1)

  !
  !    final calculation with all states
  !
900 continue

  do ns=1,nwfts
     is=iswts(ns) 
     if (octs(ns).gt.-1.0_dp) then
        if (pseudotype == 1) then
           if ( rel < 2 .or. llts(ns) == 0 .or. &
                abs(jjts(ns)-llts(ns)+0.5_dp) < 0.001_dp) then
              ind=1
           else if ( rel == 2 .and. llts(ns) > 0 .and. &
                abs(jjts(ns)-llts(ns)-0.5_dp) < 0.001_dp) then
              ind=2
           endif
           do n=1,mesh
              vaux(n)=vpstot(n,is)+vnl (n,llts(ns),ind)
           enddo
        else
           do n=1,mesh
              vaux(n)=vpstot(n,is)
           enddo
        endif
        !
        call ascheqps(nnts(ns), llts(ns), jjts(ns), enls(ns),  &
             mesh, ndm, dx, r, r2, sqr, vaux(1), thresh, phis(1,ns), &
             betas, ddd(1,1,is), qq, nbf, nwfsx, lls, jjs, ikk)

        !           write(6,*) ns, nnts(ns),llts(ns), enls(ns)  
     endif
  enddo

  call normalize ( )
  if (.not.lpaw) then
     call elsdps ( )
  else
     call new_paw_hamiltonian (vnew, dddnew, etots, &
          pawsetup, nwfts, llts, nspin, iswts, octs, phis, enls)
  endif
  !
  !   if iswitch=3 we write on the pseudopotential file the calculated 
  !   selfconsistent wavefunctions
  !
  if (iswitch.eq.3) then
     phits=0.0_dp
     do nst=1,nwfts
        if (octs(nst).ge.0.0_dp) then
           phits(:,nst)=phis(:,nst)  
        endif
     enddo
  else
     phits=phis
  endif

  if (file_recon.ne.' ')  call write_paw_recon ( )

  !
  !    compute logarithmic derivatives
  !
  if ( deld > 0.0_dp) call lderivps ( )

  return
end subroutine run_pseudo
