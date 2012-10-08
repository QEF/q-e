!
! Copyright (C) 2004 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!---------------------------------------------------------------
subroutine partial_wave_expansion
  !---------------------------------------------------------------
  !
  !  numerical integration of the radial schroedinger equation 
  !  computing logarithmic derivatives for pseudo-potentials
  !  multiple nonlocal projectors are allowed
  !
  use kinds,     only : DP
  use radial_grids, only : ndmx
  use io_global, only : stdout
  use mp,        only : mp_bcast
  use radial_grids, only: series
  use ld1_parameters, only : nwfsx
  use ld1inc,    only : grid, nld, nbeta, nspin, rel, ikk, file_pawexp, &
                        betas, ddd, qq, lls, jjs, pseudotype, vpstot, vnl, &
                        rpwe, npte, emaxld, eminld, deld, phis, rcutus, rcloc
  implicit none

  integer  ::       &
       lam,   &      ! the angular momentum
       ikrld, &      ! index of matching radius
       nc,    &      ! counter on logarithmic derivatives
       nbf,   &      ! number of b functions
       n,ie          ! generic counters

  real(DP) ::  &
       norm, norm2,& ! auxiliary variables
       ze2,     &    ! the nuclear charge in Ry units
       jam,     &    ! the total angular momentum
       e,       &    ! the eigenvalue
       lamsq,   &    ! combined angular momentum
       b(0:3),  &    ! used for starting guess of the solution 
       ddx12     

  real(DP),allocatable :: &
       ene(:),      &  ! the energy grid
       nnn(:,:),    &  ! the expansion fraction
       vaux(:),     &  ! auxiliary: the potential 
       aux(:),      &  ! the square of the wavefunction
       aux2(:),     &  ! auxiliary wavefunction
       al(:)           ! the known part of the differential equation

  real(DP), external :: compute_log
  real(DP), external :: int_0_inf_dr

  integer :: &
       ik, jb,        &  ! auxiliary variables
       ikmin,         &  ! minimum value of ik
       nst,           &  ! auxiliary for integrals
       ios,           &  ! used for I/O control
       is, ind           ! counters on index

  logical :: found

  character(len=256) :: flld


  if (nld == 0 .or. file_pawexp == ' ') return
  if (nld > nwfsx) call errore('partial_wave_expansion','nld is too large',1)

  allocate( al(grid%mesh), aux(grid%mesh), aux2(grid%mesh),vaux(grid%mesh) )

  ze2=0.0_dp

  do n=1,grid%mesh
     if (grid%r(n) > rpwe) go to 10
  enddo
  call errore('partial_wave_expansion','wrong rpwe?',1)
10 ikrld = n-1
  write(stdout,'(5x,''Computing the partial wave expansion '')') 
  npte= (emaxld-eminld)/deld + 1
  allocate ( nnn(npte,nld) )
  allocate ( ene(npte) )
  do ie=1,npte
     ene(ie)=eminld+deld*(ie-1)
  enddo
  ikmin=ikrld+5
  if (nbeta>0) then
     do nbf=1,nbeta
        ikmin=max(ikmin,ikk(nbf))
     enddo
  endif
  spinloop : &
  do is=1,nspin
     nlogdloop : &
     do nc=1,nld
        if (rel < 2) then
           lam=nc-1
           jam=0.0_dp
        else
           lam=nc/2
           if (mod(nc,2)==0) jam=lam-0.5_dp
           if (mod(nc,2)==1) jam=lam+0.5_dp
        endif
        ddx12=grid%dx*grid%dx/12.0_dp
        nst=(lam+1)*2  
        nbf=nbeta
        if (pseudotype == 1) then
           if (rel < 2 .or. lam == 0 .or. abs(jam-lam+0.5_dp) < 0.001_dp) then
              ind=1
           else if (rel==2 .and. lam>0 .and. abs(jam-lam-0.5_dp)<0.001_dp) then
              ind=2
           endif
           do n=1,grid%mesh
              vaux(n) = vpstot(n,is) + vnl(n,lam,ind)
           enddo
           nbf=0
        else
           do n=1,grid%mesh
              vaux(n) = vpstot(n,is)
           enddo
        endif

        do n=1,4
           al(n)=vaux(n)-ze2/grid%r(n)
        enddo
        call series(al,grid%r,grid%r2,b)

        energiesloop : &
        do ie=1,npte
           e=eminld+deld*(ie-1.0_dp)
           lamsq=(lam+0.5_dp)**2
           !
           !     b) find the value of solution s in the first two points
           !
           call start_scheq( lam, e, b, grid, ze2, aux )

           do n=1,grid%mesh
              al(n)=( (vaux(n)-e)*grid%r2(n) + lamsq )*ddx12
              al(n)=1.0_dp-al(n)
           enddo

           !
           ! integrate outward at fixed energy
           !
           call integrate_outward (lam,jam,e,grid%mesh,ndmx,grid,al,b,aux,betas,ddd,&
                qq,nbf,nwfsx,lls,jjs,ikk,ikmin)
           !            
           aux(1:ikmin) = aux(1:ikmin)*grid%sqr(1:ikmin)
           !
           ! calculate the accuracy of the expasion of the solution at 
           ! a given energy in terms of the partial waves at the chosen
           ! reference energies.
           ! a vanishing value of nnn indicates a perfect expansion
           !
           aux2(:) = 0.d0
           ik=min(ikrld,ikmin)
           if (mod(ik,2)==0) ik=ik+1
           found = .false.
           do jb=1,nbeta
              if (lls(jb).eq.lam.and.jjs(jb).eq.jam) then
                 found = .true.
                 IF (grid%r(ik)>max(rcutus(jb),rcloc).and.ie==1) &
                    write(stdout, &
                 '(5x,"R is outside the sphere for l=",i5," j=",f5.2)') lam,jam
                 al(1:ikk(jb))=betas(1:ikk(jb),jb)*aux(1:ikk(jb))
                 norm = int_0_inf_dr(al,grid,ikk(jb),nst)
                 aux2(1:ik) = aux2(1:ik) + phis(1:ik,jb)*norm
              endif
           enddo
           if( .not. found) then
               nnn(:,nc) = 0._dp
               write(stdout, '(7x,a,i3)') "no projector for channel: ", lam
               cycle nlogdloop
           endif

           aux2(ik+1:ikmin) = 0._dp
           al(1:ik)= aux(1:ik)**2
           !
           norm=int_0_inf_dr(al,grid,ik,nst)
           !
           al(1:ik)= (aux2(1:ik)-aux(1:ik))**2
           norm2=int_0_inf_dr(al,grid,ik,nst)
           nnn(ie,nc) = norm2/norm
           !
           !
        enddo energiesloop
     enddo nlogdloop

     flld = trim(file_pawexp)
     nnn=100.0_DP*nnn
     call write_efun(flld,nnn,ene,npte,nld)
  enddo spinloop

  deallocate(ene)
  deallocate(nnn)
  deallocate(vaux, aux, al)

  return
end subroutine partial_wave_expansion
