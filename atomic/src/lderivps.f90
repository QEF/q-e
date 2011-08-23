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
  use io_global, only : stdout
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
       lamsq,   &    ! combined angular momentum
       ddx12,   &    !
       b(0:3)          ! used for starting guess of the solution 

  real(DP),allocatable :: &
       ene(:),      &  ! the energy mesh
       dlchis(:,:), &  ! the logarithmic derivatives
       vaux(:),     &  ! auxiliary: the potential 
       aux(:),      &  ! the square of the wavefunction
       al(:)           ! the known part of the differential equation

  real(DP), external :: compute_log
  real(DP), external :: int_0_inf_dr

  integer :: &
       ikmin,         &  ! minimum value of ik
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
  allocate ( ene(npte) )
  do ie=1,npte
     ene(ie)= eminld+deld*(ie-1)
  enddo

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

        ddx12=grid%dx*grid%dx/12.0_dp
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
           e=ene(ie)
           lamsq=(lam+0.5_dp)**2
           !
           !     b) find the value of solution s in the first two points
           !
           call start_scheq( lam, e, b, grid, ze2, aux )

           do n=1,grid%mesh
              al(n)=( (vaux(n)-e)*grid%r2(n) + lamsq )*ddx12
              al(n)=1.0_dp-al(n)
           enddo

           call integrate_outward (lam,jam,e,grid%mesh,ndmx,grid,al,b,aux,betas,ddd,&
                qq,nbf,nwfsx,lls,jjs,ikk,ikmin)

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

     call write_efun(flld,dlchis,ene,npte,nld)
     !
  enddo

  deallocate(ene)
  deallocate(dlchis)
  deallocate(vaux, aux, al)

  return
end subroutine lderivps
