!
! Copyright (C) 2004 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!---------------------------------------------------------------
subroutine lderivps
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
  use ld1inc,    only : grid, nld, nbeta, nspin, rel, ikk, file_logderps, &
                        betas, ddd, qq, lls, jjs, pseudotype, vpstot, vnl, &
                        rlderiv, npte, emaxld, eminld, deld
  implicit none

  integer  ::       &
       lam,   &      ! the angular momentum
       ikrld, &      ! index of matching radius
       nc,    &      ! counter on logarithmic derivatives
       nbf,   &      ! number of b functions
       n,ie          ! generic counters

  real(DP) ::  &
       ze2,     &    ! the nuclear charge in Ry units
       jam,     &    ! the total angular momentum
       e,       &    ! the eigenvalue
       lamsq,            & ! combined angular momentum
       b(0:3),c(4),      & ! used for starting guess of the solution 
       b0e, rr1,rr2, r1,r2, & ! auxiliary
       xl1, x4l6, ddx12, &
       x6l12, x8l20

  real(DP),allocatable :: &
       dlchis(:,:), &  ! the logarithmic derivatives
       vaux(:),     &  ! auxiliary: the potential 
       aux(:),      &  ! the square of the wavefunction
       al(:)           ! the known part of the differential equation

  real(DP), external ::           &
       compute_log, &
       int_0_inf_dr

  integer :: &
       ib,jb,iib,jjb, &  ! counters on beta functions
       ikmin,         &  ! minimum value of ik
       nst,nstop,     &  ! auxiliary for integrals
       ios,           &  ! used for I/O control
       is, ind           ! counters on index

  character(len=256) :: flld


  if (nld == 0 .or. file_logderps == ' ') return
  if (nld > nwfsx) call errore('lderivps','nld is too large',1)

  allocate( al(grid%mesh), aux(grid%mesh), vaux(grid%mesh) )

  ze2=0.0_dp

  do n=1,grid%mesh
     if (grid%r(n) > rlderiv) go to 10
  enddo
  call errore('lderivps','wrong rlderiv?',1)
10 ikrld = n-1
  write(stdout,'(5x,''Computing logarithmic derivative in'',f10.5)') &
       (grid%r(ikrld)+grid%r(ikrld+1))*0.5_dp
  npte= (emaxld-eminld)/deld + 1
  allocate ( dlchis(npte,nld) )
  ikmin=ikrld+5
  if (nbeta>0) then
     do nbf=1,nbeta
        ikmin=max(ikmin,ikk(nbeta))
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

           call integrate_outward (lam,jam,e,grid%mesh,ndmx,grid,al,b,aux,betas,ddd,&
                qq,nbf,nwfsx,lls,jjs,ikmin)

           !
           !    compute the logarithmic derivative and save in dlchi
           !            
           do n=-3,3
              aux(ikrld+n)= aux(ikrld+n)*grid%sqr(ikrld+n)
           enddo

           dlchis(ie,nc)=compute_log(aux(ikrld-3),grid%r(ikrld),grid%dx)
        enddo
     enddo

     if (nspin == 2 .and. is == 1) then
        flld = trim(file_logderps)//'up'
     else if (nspin == 2 .and. is == 2) then
        flld = trim(file_logderps)//'dw'
     else
        flld = trim(file_logderps)
     end if
     if (ionode) &
        open(unit=25, file=flld, status='unknown', iostat=ios, err=300 )
300  call mp_bcast(ios, ionode_id)
     call errore('lderivps','opening file '//flld, abs(ios))

     if (ionode) then
        do ie=1,npte
           e= eminld+deld*(ie-1)
           write(25,'(10f14.6)') e, (dlchis(ie,nc),nc=1,nld)
        enddo
        close(unit=25)
     endif
  enddo

  deallocate(dlchis)
  deallocate(vaux, aux, al)

  return
end subroutine lderivps
