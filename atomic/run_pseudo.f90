!
!---------------------------------------------------------------
subroutine run_pseudo
  !---------------------------------------------------------------
  !
  !     this routine is a driver to an pseudopotential calculation
  !     with the parameters given in input
  !
  !
  use ld1inc
  implicit none

  integer :: &
       ns, &    ! counter on pseudowavefunctions
       n,  &    ! counter on mesh
       is, &    ! counter on spin
       nbf,&    ! number of beta functions
       lam,nwf0,&  ! initial phi
       nstop    ! check if the routine stop

  real(kind=dp) :: &
       xc(6),    & ! coefficients of bessel
       vaux(ndm),  &   ! auxiliary variable
       vnew(ndm,2)   ! the potential

  integer :: &
       n1,n2,nst,ikl,ind

  logical :: &
       conv       ! if true convergence reached

  real(kind=dp), parameter :: thresh=1.e-10_dp
  integer, parameter :: itmax=200

  write(6,110)
110 format (/,5x,14('-'),' Testing the pseudopotential ',24('-'),/)
  !
  !    compute an initial estimate of the potential
  !
  call start_potps
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
  !     iterate to self-consistency
  !
  do iter=1,itmax
     do ns=1,nwfts
        if (octs(ns).gt.0.0_dp) then
           is=iswts(ns)
           if (pseudotype.eq.1) then
              if (rel.eq.2) then
                 if (abs(jjts(ns)-llts(ns)+0.5_dp).lt.1.e-2_dp  &
                      &          .or. llts(ns)==0 ) then
                    ind=1
                 else
                    ind=2
                 endif
                 do n=1,mesh
                    vpstot(n,is)=vpstot(n,is)+vnlo(n,llts(ns),ind)
                    vaux(n)=vnlo(n,llts(ns),ind)
                 enddo
              else
                 do n=1,mesh
                    vpstot(n,is)=vpstot(n,is)+vnl(n,llts(ns))
                    vaux(n)=vnl(n,llts(ns))
                 enddo
              endif
           endif
           call ascheqps(nnts(ns),llts(ns),jjts(ns),enls(ns),  &
                mesh,ndm,dx,r,r2,sqr,vpstot(1,is),     &
                thresh,phis(1,ns),                     &
                betas,ddd(1,1,is),qq,nbf,nwfsx,lls,jjs,ikk)
           !               write(6,*) 'run_pseu',ns, nwfts, nnts(ns),llts(ns), jjts(ns),enls(ns)  
           if (pseudotype.eq.1) then
              do n=1,mesh
                 vpstot(n,is)=vpstot(n,is)-vaux(n)
              enddo
           endif
        endif
     enddo


     call normalize
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

     !         write(6,*) 'iteration number',iter, eps0
     do is=1,nspin
        do n=1,mesh
           vpstot(n,is)=vpstot(n,is)+vpsloc(n)
        enddo
     enddo
     call newd_at
     if (conv) goto 900
  enddo
  call errore('run_pseudo','convergence not achieved',-1)

  !
  !    final calculation with all states
  !
900 continue

  do ns=1,nwfts
     is=iswts(ns) 
     if (octs(ns).gt.-1.0_dp) then
        if (pseudotype.eq.1) then
           if (rel.eq.2) then
              if (abs(jjts(ns)-llts(ns)+0.5_dp).lt.1.e-2_dp  &
                   &          .or. llts(ns)==0 ) then
                 ind=1
              else
                 ind=2
              endif
              do n=1,mesh
                 vpstot(n,is)=vpstot(n,is)+vnlo(n,llts(ns),ind)
                 vaux(n)=vnlo(n,llts(ns),ind)
              enddo
           else
              do n=1,mesh
                 vpstot(n,is)=vpstot(n,is)+vnl(n,llts(ns))
                 vaux(n)=vnl(n,llts(ns))
              enddo
           endif
        endif
        call ascheqps(nnts(ns),llts(ns),jjts(ns),enls(ns), &
             mesh,ndm,dx,r,r2,sqr,vpstot(1,is),    &
             thresh,phis(1,ns),                    &
             betas,ddd(1,1,is),qq,nbf,nwfsx,lls,jjs,ikk)

        !           write(6,*) ns, nnts(ns),llts(ns), enls(ns)  
        if (pseudotype.eq.1) then
           do n=1,mesh
              vpstot(n,is)=vpstot(n,is)-vaux(n)
           enddo
        endif
     endif
  enddo

  call normalize
  call elsdps
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
  !
  !   print results
  !
  call write_resultsps 

  if (file_recon.ne.' ')  call write_paw_recon

  write(6,120)
120 format (/,5x,14('-'), &
       ' End of pseudopotential test ',24('-'),/)
  !
  !    compute logarithmic derivatives
  !
  if (deld.gt.0.0_dp) call lderivps

  return
end subroutine run_pseudo
