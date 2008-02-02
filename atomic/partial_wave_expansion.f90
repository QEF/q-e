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
  use io_global, only : stdout, ionode, ionode_id
  use mp,        only : mp_bcast
  use radial_grids, only: series
  use ld1_parameters, only : nwfsx
  use ld1inc,    only : grid, nld, nbeta, nspin, rel, ikk, file_pawexp, &
                        betas, ddd, qq, lls, jjs, pseudotype, vpstot, vnl, &
                        rlderiv, npte, emaxld, eminld, deld, phis, rcutus, rcloc
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
       lamsq,            & ! combined angular momentum
       b(0:3),c(4),      & ! used for starting guess of the solution 
       b0e, rr1,rr2, r1,r2, & ! auxiliary
       xl1, x4l6, ddx12, &
       x6l12, x8l20

  real(DP),allocatable :: &
       nnn(:,:),    &  ! the expansion fraction
       vaux(:),     &  ! auxiliary: the potential 
       aux(:),      &  ! the square of the wavefunction
       aux2(:),     &  ! auxiliary wavefunction
       al(:)           ! the known part of the differential equation

  real(DP), external ::           &
       compute_log, &
       int_0_inf_dr

  integer :: &
       ib,jb,iib,jjb, &  ! counters on beta functions
       ik, ikb,       &  ! auxiliary variables
       ikmin,         &  ! minimum value of ik
       nst,nstop,     &  ! auxiliary for integrals
       ios,           &  ! used for I/O control
       is, ind           ! counters on index

  character(len=256) :: flld


  if (nld == 0 .or. file_pawexp == ' ') return
  if (nld > nwfsx) call errore('lderivps','nld is too large',1)

  allocate( al(grid%mesh), aux(grid%mesh), aux2(grid%mesh),vaux(grid%mesh) )

  ze2=0.0_dp

  do n=1,grid%mesh
     if (grid%r(n) > rlderiv) go to 10
  enddo
  call errore('lderivps','wrong rlderiv?',1)
10 ikrld = n-1
  write(stdout,'(5x,''Computing the partial wave expansion '')') 
  npte= (emaxld-eminld)/deld + 1
  allocate ( nnn(npte,nld) )
  ikmin=ikrld+5
  if (nbeta>0) then
     do nbf=1,nbeta
        ikmin=max(ikmin,ikk(nbf))
     enddo
  endif
  do is=1,nspin
     do nc=1,nld
        if (rel < 2) then
           lam=nc-1
           jam=0.0_dp
        else
           lam=nc/2
           if (mod(nc,2)==0) jam=lam-0.5_dp
           if (mod(nc,2)==1) jam=lam+0.5_dp
        endif
        xl1=lam+1
        x4l6=4*lam+6
        x6l12=6*lam+12
        x8l20=8*lam+20
        ddx12=grid%dx*grid%dx/12.0_dp
        nst=(lam+1)**2  
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

        do ie=1,npte
           e=eminld+deld*(ie-1.0_dp)
           lamsq=(lam+0.5_dp)**2
           !
           !     b) find the value of solution s in the first two points
           !
           b0e=b(0)-e
           c(1)=0.5_dp*ze2/xl1
           c(2)=(c(1)*ze2+b0e)/x4l6
           c(3)=(c(2)*ze2+c(1)*b0e+b(1))/x6l12
           c(4)=(c(3)*ze2+c(2)*b0e+c(1)*b(1)+b(2))/x8l20
           r1=grid%r(1)
           r2=grid%r(2)
           rr1=(1.0_dp+r1*(c(1)+r1*(c(2)+r1*(c(3)+r1*c(4)))))*r1**(lam+1)
           rr2=(1.0_dp+r2*(c(1)+r2*(c(2)+r2*(c(3)+r2*c(4)))))*r2**(lam+1)
           aux(1)=rr1/grid%sqr(1)
           aux(2)=rr2/grid%sqr(2)

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
           aux(1:ikrld+3) = aux(1:ikrld+3)*grid%sqr(1:ikrld+3)
           !
           ! calculate the accuracy of the expasion of the solution at 
           ! a given energy in terms of the partial waves at the chosen
           ! reference energies.
           ! a vanishing value of nnn indicates a perfect expansion
           !
           aux2(:) = 0.d0
           ik = 0
           do jb=1,nbeta
              if (lls(jb).eq.lam.and.jjs(jb).eq.jam) then
                 ikb = 0
                 do while (grid%r(ikb+1) < max(rcutus(jb),rcloc) )
                    ikb=ikb+1
                 end do
                 if (mod(ikb,2) == 0) ikb=ikb+1
                 if (ik==0) ik = ikb
                 if (ikb/=ik) stop "bho"
                 al(1:ik)=betas(1:ik,jb)*aux(1:ik)
                 norm=int_0_inf_dr(al,grid,ik,nst)
                 aux2(:) = aux2(:) + phis(:,jb)*norm
              endif
           enddo
           aux2(ik+1:ikmin) = 0.0
           al(1:ik)= aux(1:ik)**2
           norm=int_0_inf_dr(al,grid,ik,nst)
           al(1:ik)= (aux2(1:ik)-aux(1:ik))**2
           norm2=int_0_inf_dr(al,grid,ik,nst)
           nnn(ie,nc) = norm2/norm
!           if (mod (ie,200) ==0) &
!              write(20000+nc+ie,'(3f20.8)') (grid%r(n),aux(n),aux2(n), n=1,ikmin)
           !
           !
        enddo
     enddo

     flld = trim(file_pawexp)
     if (ionode) &
        open(unit=25, file=flld, status='unknown', iostat=ios, err=350 )
350  call mp_bcast(ios, ionode_id)
     call errore('lderivps','opening file '//flld, abs(ios))

     if (ionode) then
        do ie=1,npte
           e= eminld+deld*(ie-1)
           write(25,'(10f14.6)') e, (max(min(100*nnn(ie,nc),9d4),-9d4),nc=1,nld)
        enddo
        close(unit=25)
     endif
  enddo

  deallocate(nnn)
  deallocate(vaux, aux, al)

  return
end subroutine partial_wave_expansion
