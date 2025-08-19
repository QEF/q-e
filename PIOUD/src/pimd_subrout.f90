
!
! Copyright (C) 2025 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! Written by Aadhityan A, Lorenzo Paulatto, Michele Casula, Tommaso Morresi
!

module random_pioud
  USE kinds, ONLY : DP

  IMPLICIT NONE

contains


    !
    !------------------------------------------------------------------------
FUNCTION pioud_randy ( irand )
  !------------------------------------------------------------------------
  !! * x=pioud_randy(n): reseed with initial seed \(\text{idum}=n\) ( 0 <= n <= ic, see below)
  !!               if pioud_randy is not explicitly initialized, it will be
  !!               initialized with seed \(\text{idum}=0\) the first time it is called;
  !! * x=pioud_randy() : generate uniform REAL(DP) numbers x in [0,1].
  !
  REAL(DP) :: pioud_randy
  INTEGER, optional    :: irand
  !
  INTEGER , PARAMETER  :: m    = 714025, &
                          ia   = 1366, &
                          ic   = 150889, &
                          ntab = 97
  REAL(DP), PARAMETER  :: rm = 1.0_DP / m
  INTEGER              :: j
  INTEGER, SAVE        :: ir(ntab), iy, idum=0
  LOGICAL, SAVE        :: first=.true.
  !
  ! return 0.5_dp

  IF ( present(irand) ) THEN
     idum = MIN( ABS(irand), ic) 
     first=.true.
  END IF

  IF ( first ) THEN
     !
     first = .false.
     idum = MOD( ic - idum, m )
     !
     DO j=1,ntab
        idum=mod(ia*idum+ic,m)
        ir(j)=idum
     END DO
     idum=mod(ia*idum+ic,m)
     iy=idum
  END IF
  j=1+(ntab*iy)/m
  IF( j > ntab .OR. j <  1 ) call errore('prandy','j out of range',ABS(j)+1)
  iy=ir(j)
  pioud_randy=iy*rm
  idum=mod(ia*idum+ic,m)
  ir(j)=idum
  !
  RETURN
  !
END FUNCTION pioud_randy
!
!------------------------------------------------------------------------
SUBROUTINE set_random_seed (iseed )
  !------------------------------------------------------------------------
  !! poor-man random seed for \(\texttt{randy}\).
  !
  INTEGER, DIMENSION (8) :: itime
  INTEGER :: iseed, irand
  !
  irand = pioud_randy ( iseed )
  !
END SUBROUTINE set_random_seed
end module  

! module random_pioud
!   USE kinds, ONLY : DP

!   IMPLICIT NONE

! contains

!   FUNCTION pioud_randy (irand)
!     !! This function is simplified for testing purposes.
!     !! It always returns the same number (e.g., 0.123_DP).
!     REAL(DP) :: pioud_randy
!     INTEGER, optional :: irand

!     ! Return a fixed value for testing
!     pioud_randy = 0.123_DP
!   END FUNCTION pioud_randy

! end module random_pioud

subroutine vel_verlet_1half_irun0(intinp)
  
  use pimd_variables, only : delt,nh
  implicit none
  integer intinp
  real(8) :: enk,enh
  
  if (intinp.eq.0) call pimd_setvel(0)
    
  if(nh) call propNH(0.5d0*delt, enh)
  
  call propP(0.5d0*delt,enk)

  call propX(delt)
  
end subroutine vel_verlet_1half_irun0

subroutine vel_verlet_2half_irun0(enk,enh)
  
  use pimd_variables, only : delt,nh
  implicit none
  real(8) :: enk,enh
  
  call propP(0.5d0*delt,enk)
 
  if(nh) call propNH(0.5d0*delt, enh)
  
end subroutine vel_verlet_2half_irun0

subroutine propX(tauMD)

  use pimd_variables, only : nbeadMD,natMD,indx,ndimMD,rpos,&
                             pimp,amas,fixcm,rcm,mtot
  implicit none

  integer :: ind,i,l,k
  real(8) :: tauMD

  do k=1,nbeadMD 
    do i = 1, natMD
      ind=indx(i)
        do l=1,ndimMD 
          rpos(l,i,k) = rpos(l,i,k) + tauMD*pimp(l,i,k)/amas(ind) 
        enddo 
    end do
  end do 

  if(fixcm) then
    do k=1,nbeadMD
      do l=1,ndimMD
        rcm(l,k)=0.d0
        do i=1,natMD
          rcm(l,k)=rcm(l,k)+amas(indx(i))*rpos(l,i,k)
        enddo
        rcm(l,k)=rcm(l,k)/mtot
!       do i=1,n
!        rpos(l,i)=rpos(l,i)-rcm(l)
!       enddo
      enddo
    enddo
  endif

  return

end subroutine propX

subroutine propP( tauMD, enkin)

  use pimd_variables, ONLY : nh,PSMD,nbeadMD,natMD,indx,ndimMD,&
                             pimp,forceMD,vcm,vel,amas,fixcm,mtot
  implicit none
  real(8) :: tauMD,scaleV,scaleF,enkin
  real(8) :: dts
  real(8), parameter :: eps=1.d-5
  integer :: i,l,ind,k
        
  if(nh) then
    dts=tauMD*PSMD
    scaleV=exp(-dts)
    if(abs(dts).gt.eps) then
      scaleF=(exp(dts)-1.d0)/PSMD
    else
      scaleF=tauMD*(1.d0+0.5d0*tauMD*PSMD*(1.0d0+tauMD*PSMD/3.d0))
    endif

    do k=1,nbeadMD
      do i = 1, natMD
        ind = indx(i)
        do l=1,ndimMD
          pimp(l,i,k) = scaleV*(pimp(l,i,k) + forceMD(l,i,k)*scaleF)
        enddo
      enddo
    enddo
  else
    do k=1,nbeadMD
      do i = 1, natMD
        ind=indx(i)
        do l=1,ndimMD
          pimp(l,i,k) = pimp(l,i,k)+tauMD*forceMD(l,i,k)
        enddo
      enddo
    enddo
  endif

  enkin=0.d0
  vcm=0.d0

  do k=1,nbeadMD
    do i=1,natMD
      ind=indx(i)
      do l=1,ndimMD
        vel(l,i,k) = pimp(l,i,k)/amas(ind)
        enkin=enkin+amas(ind)*vel(l,i,k)**2
      enddo
    enddo
  enddo
  enkin=0.5d0*enkin/nbeadMD
  
  if(fixcm) then
    do k=1,nbeadMD
      do l=1,ndimMD
        vcm(l,k)=0.d0
        do i=1,natMD
          vcm(l,k)=vcm(l,k)+amas(indx(i))*vel(l,i,k)
        enddo
        vcm(l,k)=vcm(l,k)/mtot
      enddo

    enddo
  endif

  return

end subroutine propP

subroutine propNH( tauMD, enh)

  use pimd_variables, ONLY : nbeadMD,natMD,indx,ndimMD,amas, &
                             vel,SMD,PSMD,tempMD,QMD,gMD
  
  implicit none
  real(8) :: tauMD, enh,pps,sumv2
  integer :: i,l,ind,k 
       
  sumv2 = 0.d0
       
  do k=1,nbeadMD 
    do i = 1, natMD
      ind = indx(i)
      do l=1,ndimMD
        sumv2 = sumv2 + amas(ind)*vel(l,i,k)**2
      enddo
    enddo
  enddo 

  sumv2=sumv2/nbeadMD
  SMD = SMD + 0.5d0*tauMD*PSMD
  pps=(sumv2-gMD*tempMD)/QMD
  PSMD = PSMD + tauMD * pps
  SMD = SMD + 0.5d0*tauMD*PSMD
  enh=QMD*0.5d0*PSMD**2 +gMD*SMD*tempMD

  return

end subroutine propNH

subroutine prop_ceriotti_1half_irun3(intinp)
   
  use pimd_variables, ONLY : nbeadMD,tmes_bead,mass_ion,vel,gMD,&
                             cost,fbead,rcentroid,rpos,forceMD,&
                             ekinq,tempMD,natMD,ekinqp,tfakeMD,&
                             ptilde,pimp,delt,rtilde_mode,&
                             yesquantum,ndimMD,sigma_true, rpos_old
  implicit none
  integer :: i,l,ind,k,kk,jj,ii,intinp
  real(8) :: ekinetic
  real(8), allocatable :: rpostmp(:,:,:)
  logical :: stat
 
  if (intinp.eq.0) call pimd_setvel(0)
 
! Instantaneous temperature of each bead
  do k=1,nbeadMD
    tmes_bead(k)=sum(mass_ion(:,:)*vel(:,:,k)**2)/gMD  ! gMD= # degrees of freedom
  enddo
  
  rpos_old=rpos
! potential on the configuration (averaged over quantum images)

 !! call pot(epot,epot_centroid)

  sigma_true=0.d0

  if(yesquantum) then

! Quantum kinetic energy

    allocate(rpostmp(ndimMD,natMD,nbeadMD))
    rpostmp=rpos
    do k=2,nbeadMD
      call my_refold(k,rpostmp(:,:,k-1),rpostmp(:,:,k))
    end do
    
  ! virial
    cost=0.d0
    do ii=1,nbeadMD
      fbead=0.d0
      do jj=1,natMD
        call my_mimage(rcentroid(:,jj)-rpostmp(:,jj,ii),fbead(:,jj))
      end do
      do jj=1,natMD
        do kk=1,ndimMD
          cost=cost+0.5d0*fbead(kk,jj)*forceMD(kk,jj,ii) ! Ry units!!
        enddo
      enddo
    enddo
    ekinq=cost/nbeadMD+1.5d0*tempMD*natMD

  ! primitive
    ekinqp=0.d0
    do k=1,nbeadMD
       if(k.gt.1 .and. k.lt.nbeadMD) then                 
          ekinqp=ekinqp+sum(mass_ion(:,:)*(rpostmp(:,:,k)-rpostmp(:,:,k-1))**2)
       elseif(k.eq.1) then
          ekinqp=ekinqp+sum(mass_ion(:,:)*(rpostmp(:,:,k)-rpostmp(:,:,nbeadMD))**2)
       else
          ekinqp=ekinqp+sum(mass_ion(:,:)*(rpostmp(:,:,k)-rpostmp(:,:,k-1))**2)
       endif
    enddo
    ekinqp=-ekinqp/nbeadMD*tempMD**2*0.5d0+1.5d0*natMD*tempMD ! Ha units!!!
    deallocate(rpostmp)

 ! Propagation in the quantum case
    call normal_modes_fw(ptilde,pimp,.true.) ! half step for the friction in normal modes representation
  
    call normal_modes_bw(ptilde,pimp,.false.) ! back to the real space
  
  
  else 
    call propGamma
  endif ! for yesquantum

  call propP(delt/2.d0,ekinetic) ! half step for the momenta
  if (yesquantum) then 
     call normal_modes_fw(ptilde,pimp,.false.) ! normal modes for momenta
     call normal_modes_fw(rtilde_mode,rpos,.false.) ! normal modes for positions
     call propHarm(delt) ! free ring polymer propagation in normal modes space
     call normal_modes_bw(ptilde,pimp,.false.) ! back to the real space for momenta
     call normal_modes_bw(rtilde_mode,rpos,.false.) ! back to the real space for positions
  else
     call propX(delt)
  endif

end subroutine prop_ceriotti_1half_irun3
 
subroutine prop_ceriotti_2half_irun3(ekinetic,epot)

  use pimd_variables, only : delt,yesquantum,ptilde,&
                             pimp,nbeadMD,natMD,indx,ndimMD,&
                             vel,amas,deltahtilde,forceMD,forceMD_old,&
                             rpos,rpos_old,tmes,tmes_bead,ekinq,&
                             unit_dot_sigma,sigma_true,ekinqp,epot_old,&
                             ener_true
  implicit none
  integer :: i,l,ind,k,kk,jj,ii
  real(8) :: ekin,epot,energyq,ekinetic
  
  ener_true=epot
   
  call propP(delt/2.d0,ekinetic)       ! second half step for the momenta
   
  if (yesquantum) then 
     call normal_modes_fw(ptilde,pimp,.true.) ! second half step for the friction in normal modes representation
     
     call normal_modes_bw(ptilde,pimp,.false.) ! back to the real spaces 

  else
     call propGamma
  endif

! Updating velocities (only momenta are updated during the propagation)
  do k=1,nbeadMD
    do i=1,natMD
      ind=indx(i)
      do l=1,ndimMD
         vel(l,i,k)=pimp(l,i,k)/amas(ind)
      enddo
    enddo
  enddo

! Evaluation of deltaHtilde --> conserved quantity 
! (eq. 19 of the paper 'Accurate sampling using Langevin dynamics' - Bussi)
  do k=1,nbeadMD
    do i=1,natMD
      ind=indx(i)
      do l=1,ndimMD
        deltahtilde=deltahtilde+(rpos(l,i,k)-rpos_old(l,i,k))*(forceMD(l,i,k) + forceMD_old(l,i,k))*0.5d0 &
                    + delt**2/8.d0/amas(ind)*(forceMD(l,i,k)**2 &
                    - forceMD_old(l,i,k)**2)  
      enddo
    enddo 
  enddo 
  deltahtilde=deltahtilde+epot-epot_old

  if(yesquantum) then
  !  tmes already in Ha!!!!
!    true target temperature
    tmes=sum(tmes_bead(:))/nbeadMD**2
    energyq=ener_true+ekinq
    write(unit_dot_sigma,123) energyq,sigma_true,tmes,ekinq,ekinqp 
  else
    tmes=tmes_bead(1)   
    write(unit_dot_sigma,123) ener_true,sigma_true,tmes!,((forceMD_old(l,i,1),l=1,ndimMD),i=1,n)
  endif
  flush(unit_dot_sigma)

123 format(1000000e15.7)
   
  return

end subroutine prop_ceriotti_2half_irun3

subroutine normal_modes_fw(coord_mode,coord,first)
  
  use pimd_variables, only : yesglobal,natMD,ndimMD,nbeadMD,&
                             indx,amas,cov,cmatrix,cost1,friction_mode,&
                             gamma_eigen,gammaMD,delt,cost2,tfakeMD,pi
  USE random_numbers, ONLY: randy
  USE random_pioud, ONLY: pioud_randy
  
  implicit none
  integer :: i,k,l,kk,ind
  real(8) :: drand1,drand2,xi,alpha2,alpha,kin,sumxi2,xi1,testsign
  real(8), dimension(ndimMD,natMD,nbeadMD) :: coord_mode,coord
  real(8), dimension(:,:), allocatable :: sov5,sov6
  real(8), dimension(:), allocatable :: pimpion,sov4
  logical :: first


   
  coord_mode=0.d0
       
  if (yesglobal) then
     kin=0.d0
     sumxi2=0.d0
     xi1=0.d0
  endif

     do k=1,nbeadMD
       do i=1,natMD
         ind=indx(i)
         do l=1,ndimMD
           do kk=1,nbeadMD
             coord_mode(l,i,k) = coord_mode(l,i,k) + coord(l,i,kk)*cmatrix(kk,k)
           enddo
           if (first) then
             if (.not. yesglobal .or. k .gt. 1) then ! PILE_L + PILE_G for modes > 0
                drand1 = pioud_randy()  
                drand2 = pioud_randy()
                ! call random_number(drand1)
                ! call random_number(drand2)
                xi=dsqrt(-2.d0*dlog(1.d0-drand1))*dcos(2.d0*pi*drand2)
                coord_mode(l,i,k) = cost1(k)*coord_mode(l,i,k) + sqrt(amas(ind)*tfakeMD)*cost2(k)*xi
             endif
           endif
         enddo
       enddo
     enddo

     if (first) then
       if (yesglobal) then ! PILE_G for mode = 0 => To be checked because I suspect it does not work properly
          do i=1,natMD
            ind=indx(i)
            do l=1,ndimMD
              kin=kin+coord_mode(l,i,1)**2/amas(ind)
              drand1 = pioud_randy()  
              drand2 = pioud_randy()
              ! call random_number(drand1)
              ! call random_number(drand2)
              xi=dsqrt(-2.d0*dlog(1.d0-drand1))*dcos(2.d0*pi*drand2)
              sumxi2=sumxi2+xi**2
              if (i .eq. 1) xi1=xi1+xi
            enddo 
          enddo
          kin=kin/2.d0 
          alpha2=cost1(1)+cost2(1)*sumxi2*tfakeMD/(2.d0*kin)
          alpha2=alpha2+2.d0*xi1*dsqrt(cost1(1)*cost2(1)*tfakeMD/(2.d0*kin))
          testsign=xi1+dsqrt(2.d0*kin*cost1(1)/(cost2(1)*tfakeMD))
          if (testsign .ge. 0.d0) then 
             alpha=dsqrt(alpha2)
          else
             alpha=-dsqrt(alpha2)
          endif
          coord_mode(l,i,1)=alpha*coord_mode(l,i,1)
        endif
     endif



  return

end subroutine normal_modes_fw

subroutine normal_modes_bw(coord_mode,coord,first)
    
  use pimd_variables, only : ndimMD,natMD,nbeadMD,cmatrix,&
                             amas,cov,indx
  
  implicit none

  integer :: i,k,l,kk,ind
  real(8),dimension(ndimMD,natMD,nbeadMD) :: coord,coord_mode
  real(8),dimension(:,:), allocatable :: sov5,sov6
  real(8),dimension(:), allocatable :: sov4,pimpion
  logical :: first

  coord=0.d0
 
     do k=1,nbeadMD
       do i=1,natMD
         do l=1,ndimMD
           do kk=1,nbeadMD
             coord(l,i,k) = coord(l,i,k) + cmatrix(k,kk)*coord_mode(l,i,kk)
           enddo
         enddo
       enddo
     enddo


  return

end subroutine normal_modes_bw   

subroutine propHarm(timestep)
      
  use pimd_variables, only : ptilde,rtilde_mode,omega_mode,&
                             nbeadMD,indx,ndimMD,natMD,amas
  
  implicit none

  integer :: i,k,l,ind
  real(8) :: cosit,sinut,timestep
  real(8),dimension(ndimMD,natMD,nbeadMD) :: rtilde0,ptilde0

  ptilde0=ptilde
  rtilde0=rtilde_mode
      
  do k=1,nbeadMD
    cosit=cos(omega_mode(k)*timestep)
    sinut=sin(omega_mode(k)*timestep)
    do i=1,natMD
      ind=indx(i)
      do l=1,ndimMD
        ptilde(l,i,k) = cosit*ptilde0(l,i,k) - amas(ind)*omega_mode(k)*sinut*rtilde0(l,i,k)
        if (k .eq. 1) then 
           rtilde_mode(l,i,k) = cosit*rtilde0(l,i,k) + 1.d0/amas(ind)*timestep*ptilde0(l,i,k) 
        else
           rtilde_mode(l,i,k) = cosit*rtilde0(l,i,k) + 1.d0/(amas(ind)*omega_mode(k))*sinut*ptilde0(l,i,k)
        endif      
      enddo 
    enddo
  enddo
    
  return

end subroutine propHarm

subroutine propGamma
     
  use pimd_variables, only : natMD,ndimMD,indx,pimp,amas,&
                             cov,alphaqmc_eig,delt,tempMD,&
                             gamma_eigen,cost1,cost2,pi

  USE random_numbers, ONLY: randy
  USE random_pioud, ONLY: pioud_randy


  implicit none
   
  integer :: i,l,ind
  real(8) :: drand1,drand2,xi,correct_noise,cost0
  real(8), dimension(:), allocatable :: pimpion,sov4
      
     do i=1,natMD
       ind=indx(i)
       do l=1,ndimMD
         drand1 = pioud_randy()  
         drand2 = pioud_randy()
        !  call random_number(drand1)
        !  call random_number(drand2)
         xi=dsqrt(-2.d0*dlog(1.d0-drand1))*dcos(2.d0*pi*drand2)
         pimp(l,i,1) = cost1(1)*pimp(l,i,1) + sqrt(amas(ind)*tempMD)*cost2(1)*xi
       enddo
     enddo
 
  return
 
end subroutine propGamma

subroutine prop_pioud_irun4(ekin,epot)

! Hartree units

  use pimd_variables, only : nbead, natMD, indx, ndimMD, amas, vel, ener_true,&
                             sigma_true, fk, yesquantum, fbead, rcentroid,&
                             rpos ,ekinq, cost, temp, nion, ekinqp, mass_ion,&
                             cov, cov_pimd, velocity, forcedyn, forceMD,&
                             ieskin, averaged_cov, psip, &
                             scalpar,iflagerr, rankMD, friction, &
                             dt, tmes_bead, &
                             normcorr, ener_true, tmes, unit_dot_sigma
  
  implicit none 

  integer :: i,l,ind,jj,kk,ii,k
  real(8) :: ekin,epot,amassq,energyq
  real(8) :: delta_noise
  real(8), allocatable :: rpostmp(:,:,:)
  logical :: yeswrite_loc,stat
  
  ener_true=epot
! Force fluctuations
  fk(:,:)=0.d0

! Force calculation
  !call force0(forceMD)

! kinetic energy averaged over quantum images
  ekin=0.d0
  do k=1,nbead
     do i=1,natMD
        ind=indx(i)
        do l=1,ndimMD
           ekin=ekin+amas(ind)*vel(l,i,k)**2
        enddo
     enddo
  enddo
  ekin=0.5d0*ekin/nbead

! potential on the configuration (averaged over quantum images)

  ener_true=epot
  sigma_true=0.d0

  if(yesquantum) then

! Quantum kinetic energy
     allocate(rpostmp(ndimMD,natMD,nbead))
     rpostmp=rpos
     do k=2,nbead
       call my_refold(k,rpostmp(:,:,k-1),rpostmp(:,:,k))
     end do
    
  ! virial
     cost=0.d0
     do ii=1,nbead
        fbead=0.d0
        do jj=1,natMD
          call my_mimage(rcentroid(:,jj)-rpostmp(:,jj,ii),fbead(:,jj))
        end do
        do jj=1,natMD
           do kk=1,ndimMD
              cost=cost+0.5d0*fbead(kk,jj)*forceMD(kk,jj,ii) 
           enddo
        enddo
     enddo
     ekinq=cost/nbead+1.5d0*temp*nion/nbead

! primitive
    ekinqp=0.d0
    do k=1,nbead
       if(k.gt.1 .and. k.lt.nbead) then                 
          ekinqp=ekinqp+sum(mass_ion(:,:)*(rpostmp(:,:,k)-rpostmp(:,:,k-1))**2)
       elseif(k.eq.1) then
          ekinqp=ekinqp+sum(mass_ion(:,:)*(rpostmp(:,:,k)-rpostmp(:,:,nbead))**2)
       else
          ekinqp=ekinqp+sum(mass_ion(:,:)*(rpostmp(:,:,k)-rpostmp(:,:,k-1))**2)
       endif
    enddo
    ekinqp=-ekinqp/nbead*temp**2*0.5d0+1.5d0*nion*temp ! Ha units!!!
    deallocate(rpostmp) 

  endif

! upload forces, positions, and velocities in the proper TurboRVB vector!!!!!
! ion positions in a_0
! forces must be converted in Ry
! velocities in Ry* (scaled with the masses)

  cov_pimd=0.d0

  do k=1,nbead

     do i=1,natMD
        ind=indx(i)
        amassq=sqrt(amas(ind))
        do l=1,ndimMD 
           velocity(3,l+(i-1)*ndimMD,k)=vel(l,i,k)*amassq  !! qui moltiplico v0_i*sqrt(m_i) ===> è come se passo ai momenti nuovi!!
                                                           !! infatti ora il vettore contiene p0_i/sqrt(m_i)
           forcedyn(l+(i-1)*ndimMD,k)=forceMD(l,i,k)
        enddo
     enddo
  enddo

     
  call reweight_pioud_irun4(psip,ieskin &
                           ,scalpar,iflagerr,rankMD,temp,friction &
                           ,dt,tmes_bead,cov_pimd &
                           ,nion,normcorr,forcedyn,rpos,velocity)

! rescale velocity into atomic units 
  do i=1,natMD
     ind = indx(i)
     amassq=sqrt(amas(ind))
     do k=1,nbead
        do l=1,ndimMD
           vel(l,i,k)=velocity(3,l+(i-1)*ndimMD,k)/amassq
        enddo
     enddo
  enddo

  
  if(yesquantum) then
! tmes already in Ha!!!!
! true target temperature
     tmes=sum(tmes_bead(:))/nbead**2
! fake temperature (T_fake = T_physical * nbead)
!    tmes=sum(tmes_bead(:))/nbead
!    ekin=ieskin*tmes/2.d0
     energyq=ener_true+ekinq
     write(unit_dot_sigma,123) energyq,sigma_true,tmes,ekinq,ekinqp
  else
     tmes=tmes_bead(1)   
     write(unit_dot_sigma,123) ener_true,sigma_true,tmes
!     ekin=ieskin*tmes/2.d0  
  endif
  
  flush(unit_dot_sigma)

123 format(1000000e15.7)

  return

end subroutine prop_pioud_irun4

subroutine reweight_pioud_irun4(psip,ieskin,scalpar,iflagerr,rank &
                              ,temp,friction &
                              ,dtr,tmes,cov &
                              ,nion,normcorr,forza,rion,velocity)

 use pimd_variables, only: kdyn_eig,yessecond &
                         ,write_cov,kdyn,nbead,ndimMD
 USE ring_io_units_module,         ONLY : iunpath
 USE random_numbers, ONLY: randy
 USE random_pioud, ONLY: pioud_randy

! three vectors as fundamental input (bead dependent)
! rion, velion, forza
! forza(ieskin,*),velocity(3,ieskin,*),rion(3,nion,*)
! velocity(1,ieskin,*) and velocity(3,ieskin,*) are "input/output sensitive"
! velocity(1,ieskin,*) "intermediate" positions stored (to be used in the next iteration, updated in output)
! velocity(3,ieskin,*) actual velocities (updated in output)
! velocity(2,ieskin,*) scratch vector
! forza(ieskin,*) ion forces (only input parameter)
! rion(3,nion,*)  ion positions (updated in output)

! cov matrix bead specific!!!

  implicit none 
  integer nion,i,info &
          ,iflagerr,ieskin,ind  &
          ,n2,n3,n4,n5,maxnm,jj,ii   &
          ,indin,ieskin2,kk,kkf,j,ibead,ieskin2_pimd
  real*8 normcorr,psip(*),cost,zeta &
         ,scalpar(*),forza(ieskin,*),velocity(3,ieskin,*) &
         ,temp,pi,dt,friction,tmes(*),dtr &
         ,cov(ieskin*ieskin,*) &
         ,rion(ndimMD,nion,*),cov_sav(ieskin*ieskin)
  real*8 mnoise(2,2),zetan(2),costh,eta_v,eta_r &
         ,Tn,Gn,Gni(2,2),alphaqmc,alphaall &
         ,dth,Gnh,Gnih
  real*8 drand1,drand2,dnrm2
  integer  errnoise
  real(8), dimension(:,:), allocatable:: sov4
  real(8), dimension(:,:,:), allocatable:: sov5
  real(8), dimension(:,:), allocatable:: velion
  real(8), dimension(:,:,:), allocatable:: vel0
  parameter(pi=3.14159265358979323846D0) 
  integer rank
  logical info_check
                    
  info_check=.false.
                           
  if(iflagerr.ne.0) return 
  
  errnoise=0 
  maxnm=ieskin
  n2=maxnm+1
  n3=maxnm+n2 
  n4=maxnm+n3 
  n5=maxnm+n4 
  ieskin2=ieskin*ieskin 
  ieskin2_pimd=ieskin2*nbead

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! bead independent: GLOBAL variables 

! NOTE : psip(1:ieskin) and psip(n3:n3+ieskin-1) are filled before calling ion_dynamics
! they must not be overwritten !!!!
! psip(1:ieskin) ---> eigenvalues of gamma matrix
! psip(n3:n3+ieskin-1) ---> \sqrt{m_0/m_i} with m_0 mass of the first atom in the geometry

  allocate(vel0(3,ieskin,nbead))
  allocate(velion(3,ieskin))
  allocate(sov4(ieskin,2))
  sov4=0.d0

! time step
  dt=dtr
  if(yessecond) then 
     dth=dt/2.d0
  else
     dth=dt
  endif

! mass scaling factor
  do j=1,ieskin 
     psip(n3+j-1)=dsqrt(scalpar(j))
  enddo

!! gamma deve avere dimensioni di (1/tempo); temp ha dimensioni di energia ==> cov deve avere dimensioni di [energia/tempo]
! construct gamma matrix ---> cov/2T
! covariance computed outside (cov given in input)!!!

  if(temp.ne.0.d0) then 
     cost=0.5d0/temp !!! cost ha le dimensioni di un'azione 
     cost=0.d0 
  endif
  call dscal(ieskin2_pimd,cost,cov,1) !!! cov ora ha le dimensioni di (gamma=alfa/2T)


  do ibead=1,nbead
 
! forces stored in sov4(1:ieskin,2) 
     call dcopy(ieskin,forza(1,ibead),1,sov4(1,2),1)

! scale the force by the Mass (in realtà per sqrt(mass_ion(1,1)/mass_ion(i,j)) )                               
     do i=1,ieskin 
        sov4(i,2)=sov4(i,2)*psip(n3+i-1)   
     enddo

! velocities stored in velion(1:ieskin,3)
     vel0(:,:,ibead)=velocity(:,:,ibead)
     velion(:,:)=velocity(:,:,ibead)

! each bead has its own covariance matrix

! scale covariance by masses   (in realtà per mass_ion(1,1)/mass_ion(i,j) )                                                                   
     do j=1,ieskin 
        do i=1,ieskin 
           cov(ieskin*(j-1)+i,ibead)=cov(ieskin*(j-1)+i,ibead)*psip(n3+i-1)*psip(n3+j-1) 
        enddo
     enddo

! add diagonal contribution to gamma matrix
     do i=1,ieskin 
        cov(ieskin*(i-1)+i,ibead)=cov(ieskin*(i-1)+i,ibead)
     enddo
     
! save cov for output before it is destroyed by covariance diagonalization
     if(write_cov) then
        cov_sav(1:ieskin*ieskin)=cov(1:ieskin*ieskin,ibead)
     endif
  
!  do j=1,ieskin
!     write(6,*) j,psip(n3+j-1)
!     do i=1,ieskin
!        write(6,*) i,j,cov(ieskin*(j-1)+i)
!     enddo
!  enddo

! diagonalize gamma matrix
     call dsyev_my('V','L',ieskin,cov(1,ibead),ieskin,psip,info)  !! psip(1:ieskin) contiene gli autovalori di gamma ma 
     ! c'è il fattore di prima psip(n3:4*ieskin)
 
    ! if(ibead.eq.1) then 
    !    write(6,*) ' Eigenvalues covariance'
    !    do i=1,ieskin
    !       write(6,*) i,psip(i)
    !    enddo
!       write(6,*) ' Eigenvectors covariance'
!       do i=1,ieskin
!          write(6,*) i,(cov(ieskin*(i-1)+j),j=1,ieskin)
!       enddo
    ! endif
  
     if(info.ne.0) then 
        write(iunpath,*) ' Error in lapack dsyev dynamic !!! ' 
        errnoise=4 
        info_check=.true.
     endif

!  if(rank.eq.0) write(6,*) ' Ratio dyn =',(psip(ieskin)-friction)/sqrt(temp*157.8873306d0)
     if(ibead.eq.1 .and. rank.eq.0) write(6,*) ' Ratio dyn =',(psip(ieskin))*dt !!!Aadhityan      

!  write(6,*) psip(ieskin),friction,temp

     do i=1,ieskin 
        if(psip(i).gt.0.d0) then 
           psip(n5+i-1)=(1.d0-dexp(-psip(i)*dth))/psip(i) 
        else 
           if(rank.eq.0) write(iunpath,*) ' Refused eigenvalue ',i,ibead 
           psip(n5+i-1)=dth 
        endif
     enddo
!         compute the noise correction                                  
     alphaqmc=0.d0
                                                                        
     do i=1,ieskin 

        if(psip(i).gt.0.d0) then 
!         the following is protected for overflow                       
           cost=exp(-dth*psip(i)) 
           alphaqmc=normcorr*2.d0*temp*(psip(i))
           psip(n4+i-1)=temp*psip(i)**2*(1.d0+cost)/(1.d0-cost)-alphaqmc                          
!          subtracting the noise already present in the forces          
        else 
           psip(n4+i-1)=0.d0 
        endif

        if(psip(n4+i-1).gt.0.d0) then 
           psip(n4+i-1)=dsqrt(psip(n4+i-1)) 
        else 
           psip(n4+i-1)=0.d0 
        
              if(rank.eq.0) write(iunpath,*) 'There should be some error in reweight0',i,ibead,psip(n4+i-1)
              errnoise=5 
           
        endif

     enddo
                                                                        
! now we are able to update the velocity                        
! first change basis actual velocity and force                  
     call dgemv('T',ieskin,ieskin,1.d0,cov(1,ibead),ieskin,velion(3,1),3,0.d0,sov4,1)                               
     call dgemv('T',ieskin,ieskin,1.d0,cov(1,ibead),ieskin,sov4(1,2),1,0.d0,velion(2,1),3)                        

     if(yessecond) then
! Compute the temperature at half time intervals.
        do i=1,ieskin 
           drand1 = pioud_randy()  
           drand2 = pioud_randy()
          !  call random_number(drand1)
          !  call random_number(drand2)
           zeta=dsqrt(-2.d0*dlog(1.d0-drand1))*dcos(2.d0*pi*drand2) 
           sov4(i,1)=sov4(i,1)*dexp(-psip(i)*dth)                    &
                     +psip(n5+i-1)*(velion(2,i)+psip(n4+i-1)*zeta)                     
        enddo
        tmes(ibead)=dnrm2(ieskin,sov4,1)**2/ieskin 
! Second half time interval.
        do i=1,ieskin 
           drand1 = pioud_randy()  
           drand2 = pioud_randy()
          !  call random_number(drand1)
          !  call random_number(drand2)
           zeta=dsqrt(-2.d0*dlog(1.d0-drand1))*dcos(2.d0*pi*drand2) 
           velion(2,i)=sov4(i,1)*dexp(-psip(i)*dth)                    &
                       +psip(n5+i-1)*(velion(2,i)+psip(n4+i-1)*zeta)                     
        enddo
     else
        do i=1,ieskin
           drand1 = pioud_randy()  
           drand2 = pioud_randy()
          !  call random_number(drand1)
          !  call random_number(drand2)
           zeta=dsqrt(-2.d0*dlog(1.d0-drand1))*dcos(2.d0*pi*drand2) 
           velion(2,i)=sov4(i,1)*dexp(-psip(i)*dth)                    &
                       +psip(n5+i-1)*(velion(2,i)+psip(n4+i-1)*zeta)                     
        enddo
     endif
                                                                        
! go back to the original basis                                 
                                                                        
     call dgemv('N',ieskin,ieskin,1.d0,cov(1,ibead),ieskin,velion(2,1),3,0.d0,velion(3,1),3)                                
                                                                        
     if(.not.yessecond) tmes(ibead)=dnrm2(ieskin,velion(3,1),3)**2/ieskin 

! prepare for the next (harmonic) step

     velocity(3,1:ieskin,ibead)=velion(3,1:ieskin) ! velocities
     velocity(2,1:ieskin,ibead)=0.d0 ! forces 
     do j=1,nion
        do i=1,ndimMD
           ind=i+(j-1)*ndimMD
           velocity(1,ind,ibead)=rion(i,j,ibead)/psip(n3+ind-1) ! positions
        enddo
     enddo

!     write(6,*) 'velocity before'
!     do j=1,nion
!        do i=1,ndimMD
!           ind=i+(j-1)*ndimMD
!           write(6,*) i,j,ibead,velocity(1,ind,ibead)
!           write(6,*) i,j,ibead,velocity(2,ind,ibead)
!        enddo
!     enddo     

  enddo

! harmonic forces step
  allocate(sov5(3,ieskin,nbead))
  sov5=0.d0

! rotation in the kdyn basis
! write the above operations with rank-3 blas!!!!
  call dgemm('N','N',3*ieskin,nbead,nbead,1.d0,velocity,3*ieskin,kdyn,nbead,0.d0,sov5,3*ieskin)

! sov5(1,1:ieskin,ibead)  ! positions
! sov5(2,1:ieskin,ibead)  ! forces
! sov5(3,1:ieskin,ibead)  ! velocities

  do ibead=1,nbead

!!!! START equations of motion integration 

     do i=1,ieskin 
        drand1 = pioud_randy()  
        drand2 = pioud_randy()

        ! call random_number(drand1)
        ! call random_number(drand2)
        zetan(1)=dsqrt(-2.d0*dlog(1.d0-drand1))*dcos(2.d0*pi*drand2) 

        drand1 = pioud_randy()  
        drand2 = pioud_randy()
        ! call random_number(drand1)
        ! call random_number(drand2)
        zetan(2)=dsqrt(-2.d0*dlog(1.d0-drand1))*dcos(2.d0*pi*drand2) 
 
! TEST: freezing the random number generator         
!          zetan(1)=0.55
!          zetan(2)=0.75

! the QMC noisy part has been integrated in the first part
        alphaqmc=0.d0

!!!! START equations of motion integration for harmonic part
! gammaall
        psip(i)=2.d0*sqrt(abs(kdyn_eig(ibead)))
        if(psip(i).lt.friction) psip(i)=friction ! cutoff on the lowest eigenvalue. Ceriotti's choice.

        alphaall=2.d0*temp*psip(i)

        call set_turboq(psip(i),kdyn_eig(ibead),dt &
             ,alphaall,alphaqmc,mnoise,Gn,Tn,Gni,Gnh,Gnih)

        call root2mat(mnoise,errnoise)

        if(errnoise.ne.0) write(iunpath,*) ' Error negative definite matrix noise '
                                                                        
        eta_v=mnoise(1,1)*zetan(1)+mnoise(1,2)*zetan(2)                                   
        eta_r=mnoise(2,1)*zetan(1)+mnoise(2,2)*zetan(2)
          
        sov4(i,2)=sov5(3,i,ibead)*Gni(2,1)+Gni(2,2)*sov5(1,i,ibead)+Tn*eta_r
        velion(2,i)=sov5(3,i,ibead)*Gni(1,1)+Gni(1,2)*sov5(1,i,ibead)+Gn*eta_v
       
     enddo

     velocity(2,1:ieskin,ibead)=velion(2,1:ieskin) ! velocities
     velocity(1,1:ieskin,ibead)=sov4(1:ieskin,2) ! coordinatMDes

!     write(6,*) 'velocity'
!     do j=1,nion
!        do i=1,ndimMD
!           ind=i+(j-1)*ndimMD
!           write(6,*) i,j,ibead,velocity(1,ind,ibead)
!           write(6,*) i,j,ibead,velocity(2,ind,ibead)
!        enddo
!     enddo     

  enddo
                                                                        
! Now go back to the original basis                             

! write the above operations with rank-3 blas!!!!
  call dgemm('N','T',3*ieskin,nbead,nbead,1.d0,velocity,3*ieskin,kdyn,nbead,0.d0,sov5,3*ieskin)

! sov5(1,1:ieskin,ibead) coordinatMDes
! sov5(2,1:ieskin,ibead) velocities
     
  do ibead=1,nbead

     call dcopy(ieskin,sov5(2,1,ibead),3,velion(3,1),3)
     call dcopy(ieskin,sov5(1,1,ibead),3,psip(n2),1)

!     write(6,*) 'sov5'
!     do j=1,nion
!        do i=1,ndimMD
!           ind=i+(j-1)*ndimMD
!           write(6,*) i,j,ibead,sov5(1,ind,ibead)
!           write(6,*) i,j,ibead,sov5(2,ind,ibead)
!        enddo
!     enddo     

!       back to the scale of coordinatMDes                                
     do i=1,ieskin 
        psip(n2+i-1)=psip(n2+i-1)*psip(n3+i-1)
        if(psip(n3+i-1).eq.0.d0) then 
           velion(3,i)=0.d0
        endif
     enddo

     if(ibead.eq.1 .and. rank.eq.0) then
   !!     write(iunpath,*)
   !!     write(iunpath,*) ' Temperature (Ha) = ',tmes(ibead) 
        cost=dnrm2(ieskin,psip(n2),1)
   !!     write(iunpath,*) ' Norm change ions =',cost
   !!     write(iunpath,*) ' Coordinates variation for bead',ibead                              
   !!     do i=1,ieskin 
   !!        write(iunpath,*) i,psip(n2+i-1)   
    !!    enddo
    !!    write(iunpath,*)
     endif

! update velocities
     velocity(:,:,ibead)=velion(:,:)

     if(info_check) then 
        write(iunpath,*) ' Warning info =',info 
        write(iunpath,*) 'Continuing with previous param. !!!' 
! do not move ions
        call dscalzero(ieskin,0.d0,psip(n2),1) 
! restore previous velocities
        velocity(:,:,ibead)=vel0(:,:,ibead)
     endif

     kkf=0
     do kk=1,nion
        do jj=1,ndimMD
           kkf=kkf+1
           rion(jj,kk,ibead)=rion(jj,kk,ibead)+psip(n2+kkf-1)
        enddo
     enddo
     
  enddo

! over beads
  
  !!if(write_cov) then
  !!   write(18,123) ((cov_sav(ii+(jj-1)*ieskin),ii=jj,ieskin),jj=1,ieskin)
  !!endif

  if(errnoise.ne.0) then                                                           
     write(iunpath,*) ' Error in reweight0  =',errnoise 
     stop                                                                         
  endif

  deallocate(vel0,velion,sov4,sov5)
                                                                              
  return
     
end subroutine reweight_pioud_irun4

subroutine set_turboq(gamma_eig,kdyn_eig,dt,alphaall,alphaqmc,mnoise,Gn,Tn,Gni,Gnh,Gnih)

  implicit none
  
  real*8, intent(in):: gamma_eig,kdyn_eig,alphaall,alphaqmc,dt
  real*8, intent(out):: mnoise(2,2),Gn,Tn,Gni(2,2),Gnh,Gnih
  real*8 cost,costm1,xbar,sinhx,coshx,gammag,gammagx,gammagxm
  real*8 dergamma,der2gamma,costh,sinhxh,coshxh,costhm1
  real*8, external:: mygamma
  complex*16 gammagxc
  complex*16, external:: mygammac
  logical xbarpos

  if(kdyn_eig.gt.1d-6) then 
    cost=exp(-dt*gamma_eig/2.d0) 
    costh=exp(-dt*gamma_eig/4.d0) 
    xbar=(gamma_eig/2.d0)**2-kdyn_eig   !! determinante 
    
    if(xbar.ge.0) then
      xbar=sqrt(xbar)
      coshx=cosh(xbar*dt)
      sinhx=sinh(xbar*dt)
      coshxh=cosh(xbar*dt/2.d0)
      sinhxh=sinh(xbar*dt/2.d0)
      xbarpos=.true.
    else
      xbar=sqrt(-xbar)
      coshx=cos(xbar*dt)
      sinhx=sin(xbar*dt)
      coshxh=cos(xbar*dt/2.d0)
      sinhxh=sin(xbar*dt/2.d0)
      xbarpos=.false.
    endif

    if(xbar*dt.gt.1d-6) then
      Gn=cost/xbar*sinhx
      Gnh=costh/xbar*sinhxh
      Tn=(1-cost*coshx)/kdyn_eig-0.5d0*cost*gamma_eig*sinhx/(xbar*kdyn_eig)
      Gni(1,1)=cost*(coshx-0.5d0*gamma_eig/xbar*sinhx)
      Gni(1,2)=-cost*kdyn_eig*sinhx/xbar
      Gni(2,1)=Gn
      Gni(2,2)=cost*(coshx+0.5d0*gamma_eig/xbar*sinhx)-1.d0 ! -1 as we remove the identity
      Gnih=costh*(coshxh+0.5d0*gamma_eig/xbar*sinhxh)-1.d0
    else
      Gn=cost*dt
      Gnh=costh*dt/2.d0
      Tn=(1-cost*coshx)/kdyn_eig-0.5d0*cost*gamma_eig*dt/kdyn_eig
      Gni(1,1)=cost*(coshx-0.5d0*gamma_eig*dt)
      Gni(1,2)=-cost*kdyn_eig*dt
      Gni(2,1)=Gn
      Gni(2,2)=cost*(coshx+0.5d0*gamma_eig*dt)-1.d0  ! -1 as we remove the identity
      Gnih=costh*(coshxh+0.25d0*gamma_eig*dt)-1.d0
    endif


!  Protection from roundoff
    if(xbar*dt.le.1d-6) then
      xbar=1d-6/dt
    endif
    if(xbarpos) then 
      gammag=mygamma(gamma_eig,dt)
      gammagx=mygamma(gamma_eig+2*xbar,dt)
      gammagxm=mygamma(gamma_eig-2*xbar,dt)
      dergamma=(gammagx-gammagxm)/xbar
      der2gamma=(gammagx+gammagxm-2*gammag)/xbar**2
      mnoise(1,1)=alphaall*0.25d0*(gammagx+gammagxm+2*gammag)+0.0625d0*alphaall*gamma_eig**2*der2gamma+&
                  & 0.25d0*alphaall*gamma_eig*dergamma
      mnoise(2,2)=0.25d0*alphaall*der2gamma
      mnoise(1,2)=-0.125d0*alphaall*gamma_eig*der2gamma-0.25d0*alphaall*dergamma
    else
      gammag=mygamma(gamma_eig,dt)
      gammagxc=mygammac(dcmplx(gamma_eig,2*xbar),dt)
      dergamma=2*aimag(gammagxc)/xbar
      der2gamma=(2*gammagxc-2*gammag)/xbar**2
      mnoise(1,1)=alphaall*0.25d0*(2*gammagxc+2*gammag)-0.0625d0*alphaall*gamma_eig**2*der2gamma+&
                  & 0.25d0*alphaall*gamma_eig*dergamma
      mnoise(2,2)=-0.25d0*alphaall*der2gamma
      mnoise(1,2)=0.125d0*alphaall*gamma_eig*der2gamma-0.25d0*alphaall*dergamma
    endif

!     Now the bar noise to be added to the force with noise  correction

    mnoise(1,1)=mnoise(1,1)/Gn**2-alphaqmc
    mnoise(2,2)=mnoise(2,2)/Tn**2-alphaqmc
    mnoise(1,2)=mnoise(1,2)/Tn/Gn-alphaqmc
    mnoise(2,1)=mnoise(1,2)

  else

    costh=exp(-dt*gamma_eig/2.d0) 
    if(costh.lt.0.99999d0) then 
      costhm1=costh-1.d0
    else
      costhm1=-dt/2.d0*gamma_eig+(dt*gamma_eig/2.d0)**2/2.d0
    endif
    cost=exp(-dt*gamma_eig) 
    if(cost.lt.0.99999d0) then 
      costm1=cost-1.d0
    else
      costm1=-dt*gamma_eig+(dt*gamma_eig)**2/2.d0
    endif

    Gn=-costm1/gamma_eig
    Gnh=-costhm1/gamma_eig
    Gnih=0.d0
    Tn=1.d0/gamma_eig*(dt-Gn)
    Gni(1,1)=cost
    Gni(2,1)=Gn
    Gni(2,2)=0.d0
    Gni(1,2)=0.d0
!        Tn=1.d0/gamma_eig**2*(dt*gamma_eig+costm1) 
    mnoise(1,1)=-0.5d0*alphaall/gamma_eig*(1.d0+cost)*costm1

!         mnoise(1,1)=-0.5d0*alphaall*gamma_eig*(1.d0+cost)/costm1-alphaqmc
    mnoise(2,2)=alphaall/gamma_eig**2*(dt+1.d0/gamma_eig*costm1*(2.d0-0.5d0*(1.d0+cost)))

!         mnoise(1,2)=alphaall/gamma_eig*&
!    & (-costm1/gamma_eig*(1.-0.5d0*(1+cost)))

    mnoise(1,2)=0.5d0*alphaall*(costm1/gamma_eig)**2
    mnoise(1,1)=mnoise(1,1)/Gn**2-alphaqmc
    mnoise(2,2)=mnoise(2,2)/Tn**2-alphaqmc
    mnoise(1,2)=mnoise(1,2)/Tn/Gn-alphaqmc
    mnoise(2,1)=mnoise(1,2)

  endif

  return

end subroutine
               
function mygamma(x,dt)
  implicit none
  real*8 x,mygamma,dt,xdt
  xdt=x*dt
  if(abs(xdt).gt.1d-6) then
    mygamma=(1.d0-exp(-xdt))/x 
  else
    mygamma=dt-0.5d0*xdt*dt
  endif
  return
end function mygamma
        
function mygammac(x,dt)
  implicit none
  complex*16 x,xdt,mygammac
  real*8 dt
  xdt=x*dt
  if(abs(xdt).gt.1d-6) then
    mygammac=(1.d0-exp(-xdt))/x 
  else
    mygammac=dt-0.5d0*xdt*dt
  endif
  return
end function mygammac

subroutine root2mat(mat,errnoise) 

  real*8 mat(2,2),z,a,b,c,aa,bb,cc
  integer errnoise
!   This subroutine computes the square root of a positive definite
!   symmetric 2x2 matrix A
!   matrix. The output is the unique 2x2 symmetric positive definite matrix B such
!   that B^2 =A. If the input is positive definite the output replaces
!   the input, otherwise the input matrix is unchanged.
  aa=mat(1,1)
  bb=mat(2,1)
  cc=mat(2,2)
  errnoise=0
  if(aa.lt.0.d0 .or. cc.lt.0.d0 .or. aa*cc-bb**2 .lt.0.d0) then
!         check if it is positive definite
     errnoise=2
     return
  endif 

  z=aa+cc+dsqrt(4.d0*(aa*cc-bb**2))

  if(z.le.0.d0) then
!         square root of zero is zero.
     return
  endif
  z=dsqrt(z)

  a=0.5d0*(z+(aa-cc)/z)
  b=bb/z
  c=0.5d0*(z-(aa-cc)/z)
         
  mat(1,1)=a
  mat(2,1)=b
  mat(1,2)=b
  mat(2,2)=c 

  return
end subroutine root2mat

SUBROUTINE DSYEV_MY( JOBZ, UPLO, N, A, LDA, W, INFO)
   implicit none
!     .. Scalar Arguments ..
      CHARACTER          JOBZ, UPLO
      INTEGER            INFO, LDA, N,M,I
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, *), W( * )
      DOUBLE PRECISION, dimension(:,:), allocatable:: z
      double PRECISION, dimension(:), allocatable:: workp
      integer, dimension(:), allocatable:: iwork,ifail

      allocate(z(LDA,N),iwork(5*N),ifail(N),workp(8*N))

        CALL DSYEVX( JOBZ, 'A', UPLO, N, A, LDA, 0.d0, 0.d0, 1, 1,&
     &                   0.d0, M, W, Z, LDA, WORKP, 8*N, IWORK,&
     &                   IFAIL, INFO )


        do i=1,N
        A(1:N,i)=Z(1:N,i)
        enddo


      deallocate(z,iwork,ifail,workp)

      return
END SUBROUTINE DSYEV_MY

subroutine dscalzero(n,zero,vet,m)
  implicit none
  integer n,m,i
  real*8 zero,vet(m,*)
  do i=1,n
    vet(1,i)=zero
  enddo
  return
end subroutine dscalzero

SUBROUTINE pimdnvt_init(epot)!,forcetmp)!(qui devo mettere forze e potenziale)
      !-----------------------------------------------------------------------
      !
      USE pimd_variables,   ONLY : nbeadMD,natMD,indx,ndimMD,avp,ipot,av,ikin,anorm, &
                                   amas,tempMD,ekinq,ekinqp,yesquantum, irun, &
                                   unit_dot_out,iblockMD,kbm1,gMD,fbead,rpos, &
                                   rcentroid,cost,forceMD,mass_ion,forceMD_old,vel, &
                                   ekinq,ekinqp,unit_dot_out,restart_pimd

      USE ring_io_units_module,         ONLY : iunpath
      IMPLICIT NONE
      
      integer :: i,nseed
      integer, allocatable :: seed(:)
      integer :: ii,ind,jj,k,l,kk
      real(8) :: tt,enh,h,en,vs,ep,elr,enk,ttk
      REAL(8) :: epot,ekin
      real(8), allocatable :: rpostmp(:,:,:)
      
      call random_seed(size=nseed)
      allocate(seed(nseed))
      ! call system_clock(count=clock)
      do i=1,nseed
         seed(i) = 3 ! clock + 17*i
      end do
      call random_seed(put=seed)
      
      if(.not.restart_pimd) call pimd_setvel(0)  !  initial velocities of ions
      
      call zeroav(0)         ! set to zero the cumulators
      
     ! forceMD=forcetmp

!!!! **** compute the initial potential & forces... **** !!!!
   
     !! call force0(forceMD)
     !! call pot(epot)

      if (irun .eq. 3) forceMD_old=forceMD
      
      if(.not.restart_pimd) then
         write(iunpath,*)
         write(iunpath,*) '**********************************************************'
         write(iunpath,*) 'initial potential energy:',epot
      
!!!! **** ...and kinetic energy averaged over quantum images **** !!!!
         ekin=0.d0
         do k=1,nbeadMD
           do i=1,natMD
             ind=indx(i)
             do l=1,ndimMD
               ekin=ekin+amas(ind)*vel(l,i,k)**2
             enddo
           enddo
         enddo
         ekin=0.5d0*ekin/nbeadMD
         write(iunpath,*) 'initial kinetic energy:',ekin
      
         en=ekin+epot
         h=en

         ttk=2.d0*ekin/gMD*kbm1 ! Temperature in Kelvin
         tt=2.d0*ekin/gMD    ! Temperature in Hartree

         ttk=ttk/nbeadMD
         tt=tt/nbeadMD

         if(yesquantum) then
! initial EXPECTED value of Quantum kinetic energy calculated in two ways

           allocate(rpostmp(ndimMD,natMD,nbeadMD))
           rpostmp=rpos
           do k=2,nbeadMD
             call my_refold(k,rpostmp(:,:,k-1),rpostmp(:,:,k))
           end do
    
      ! virial
           cost=0.d0
           do ii=1,nbeadMD
              fbead=0.d0
              do jj=1,natMD
               call my_mimage(rcentroid(:,jj)-rpostmp(:,jj,ii),fbead(:,jj))
              end do
              do jj=1,natMD
                do kk=1,ndimMD
                  cost=cost+0.5d0*fbead(kk,jj)*forceMD(kk,jj,ii) ! Ha units!!
                enddo
              enddo
           enddo
           ekinq=cost/nbeadMD+1.5d0*tempMD*natMD/nbeadMD

! primitive
           ekinqp=0.d0
           do k=1,nbeadMD
             if(k.gt.1 .and. k.lt.nbeadMD) then                 
               ekinqp=ekinqp+sum(mass_ion(:,:)*(rpostmp(:,:,k)-rpostmp(:,:,k-1))**2)
             elseif(k.eq.1) then
               ekinqp=ekinqp+sum(mass_ion(:,:)*(rpostmp(:,:,k)-rpostmp(:,:,nbeadMD))**2)
             else
               ekinqp=ekinqp+sum(mass_ion(:,:)*(rpostmp(:,:,k)-rpostmp(:,:,k-1))**2)
             endif
           enddo
           ekinqp=-ekinqp/nbeadMD*tempMD**2*0.5d0+1.5d0*natMD*tempMD ! Ha units!!!
           deallocate(rpostmp)
        
         endif
 
         vs = epot/natMD
         write(iunpath,'('' initial potential energy per atom =  '', g20.10 )' ) vs
         write(iunpath,*) 'quantities after initialization'
         if(nbeadMD.eq.1) then
           write(iunpath,'(//''1)block   2)nmove   3)H   4)H_Nosé   5)E_pot'' &
            ''   6)E_kin   7)Temp(Ha)   8)Temp(K)''//)')
           write(iunpath,'(i2,i2,g15.7,g15.7,g15.7,g15.7,g15.7,g15.7)') &
              0,0,h,en,epot,ekin,tt,ttk
         else
           write(iunpath,'(//''1)block   2)nmove   3)H   4)H_Nosé   5)E_pot'' &
            ''   6)E_kin   7)Temp(Ha)   8)Temp(K)'' &
            ''   9)quantum_kin_virial   10)quantum_kin_primitive  ''//)')
           write(iunpath,'(i2,i2,g15.7,g15.7,g15.7,g15.7,g15.7,g15.7,g15.7,g15.7)') & 
               0,0,h,en,epot,ekin,tt,ttk,ekinq,ekinqp
         endif
         write(iunpath,*) 
         write(iunpath,*)
         write(iunpath,*)'***********************************************************'
         write(iunpath,*)'**************** start of dynamics ************************'
  
         if (.not.restart_pimd) write(unit_dot_out,'('' initial potential =  '', f10.4 )' ) epot
         if (.not.restart_pimd) write(unit_dot_out,'(//'' start of dynamics ''//)')
     
         if (.not.restart_pimd) then 
           if(nbeadMD.eq.1) then
             write(unit_dot_out,'(//''1)block   2)nmove   3)H   4)H_Nosé   5)E_pot'' &
              ''   6)E_kin   7)Temp(Ha)   8)Temp(K)''//)')
           else
             write(unit_dot_out,'(//''1)block   2)nmove   3)H   4)H_Nosé   5)E_pot'' &
              ''   6)E_kin   7)Temp(Ha)   8)Temp(K)'' &
              ''   9)quantum_kin_virial   10)quantum_kin_primitive  ''//)')
           endif
         end if
      end if !! endif restart_pimd
     
      iblockMD=1
     
      return
     
END SUBROUTINE pimdnvt_init

SUBROUTINE pimdnvt( ep)!,forcetmp)
     
     USE ring_variables,   ONLY : istep_path
     USE pimd_variables,   ONLY : iblockMD,nstep_block,irun,deltahtilde,natMD,avp, &
                                  ipot,anormp,iprint,gMD,kbm1,nbeadMD,av,anorm, &
                                  rpos_old,rpos,forceMD_old,forceMD,epot_old,&
                                  ikin,unit_dot_ek,unit_dot_ep,unit_dot_ep2,&
                                  unit_dot_eplr,unit_dot_epsr,ekinq,ekinqp,&
                                  unit_dot_out,ndimMD
     USE ring_io_units_module,         ONLY : iunpath
     IMPLICIT NONE
      
     integer :: i,ii,step_blockMD
     real(8) :: tt,enh,h,en,vs,ep,elr,enk,ttk
  !   real(8) :: forcetmp(ndimMD,natMD,nbeadMD)

   !  forceMD=forcetmp
     enh=0.d0
     enk=0.d0
     
     step_blockMD=istep_path-nstep_block*(iblockMD-1)
          
     if (irun == 0) then 
        call vel_verlet_2half_irun0(enk,enh)
     elseif (irun == 3) then
        call prop_ceriotti_2half_irun3(enk,ep)  
     elseif (irun == 4) then
        call prop_pioud_irun4(enk,ep)
     endif

     en = ep + enk  

     if (irun .eq. 0) then  
        h = en + enh
     elseif (irun .eq. 3) then 
        h = deltahtilde
     else 
        h = en
     endif
     
!  ** calculate potential energy per dumbell **
     vs=ep/natMD !! questa energia in teoria è l'energia potenziale media per atomo (atomo vero si intende)
      
!  ** accumulate step averages **
     call cumul1(vs,avp(ipot(1)),anormp(ipot(1)),1)
     call cumul1(vs**2,avp(ipot(2)),anormp(ipot(2)),1)
     call cumul1(ep,avp(ipot(3)),anormp(ipot(3)),1)
     call cumul1(elr,avp(ipot(4)),anormp(ipot(4)),1)
     call cumul1(enk/natMD,avp(ikin),anormp(ikin),1)
 
!  ** perform periodic operations  **
     if ( mod ( istep_path, iprint ) .eq. 0 ) then 
           
        ttk=2.d0*enk/(gMD*nbeadMD)*kbm1 ! Temperature in Kelvin
        tt=2.d0*enk/(gMD*nbeadMD)  ! Temperature in Hartree
        call checkpoint(ttk)
        if(nbeadMD.gt.1) then
!  ** write out runtime information in quantum case**
!           write(iunpath,*) '**********************************************************'
           write(iunpath,*)
           write(iunpath,'(2i6,8g15.7)') & 
                iblockMD,step_blockMD,h,enh,ep,enk,tt,ttk,ekinq,ekinqp
           write(iunpath,*) '**********************************************************'
           write(unit_dot_out,'(2i6,8g15.7)') & 
                iblockMD,step_blockMD,h,enh,ep,enk,tt,ttk,ekinq,ekinqp
           flush(unit_dot_out)
        else
!  ** write out runtime information in classical case**
!           write(iunpath,*) '**********************************************************'
           write(iunpath,*)
           write(iunpath,'(2i6,6g15.7)') &
                iblockMD,step_blockMD,h,enh,ep,enk,tt,ttk
           write(iunpath,*) '**********************************************************'
           write(unit_dot_out,'(2i6,6g15.7)') &
                     iblockMD,step_blockMD,h,enh,ep,enk,tt,ttk
           flush(unit_dot_out)
        endif

      endif

      if(irun.eq.3) then
! htilde is computed only for Ceriotti integrator
! Saving potential energy, forces and positions to evaluate deltaHtilde
        rpos_old=rpos
        forceMD_old=forceMD
        epot_old=ep
      endif 
     
 
! ** block analysis  
      if(mod(istep_path,nstep_block).eq.0) then
        call cumul(avp(ipot(1)),av(ipot(1)),anorm(ipot(1)),1)
        call cumul(avp(ipot(2)),av(ipot(2)),anorm(ipot(2)),1)
        call cumul(avp(ipot(3)),av(ipot(3)),anorm(ipot(3)),1)
        call cumul(avp(ipot(4)),av(ipot(4)),anorm(ipot(4)),1)
        call cumul(avp(ikin),av(ikin),anorm(ikin),1)
        call sprint(unit_dot_ep,ipot(1),1,'epot')
        flush(unit_dot_ep)
        call sprint(unit_dot_ep2,ipot(2),1,'ept2')
        flush(unit_dot_ep2)
        call sprint(unit_dot_epsr,ipot(3),1,'ept2')
        flush(unit_dot_epsr)
        call sprint(unit_dot_eplr,ipot(4),1,'ept2')
        flush(unit_dot_eplr)
        call sprint(unit_dot_ek,ikin,1,'ekin')
        flush(unit_dot_ek)
      end if                

     if(mod(istep_path,nstep_block).eq.0) then
        call zeroav(1)
        iblockMD=iblockMD+1                    ! zeroed the block cumulators
     end if 

      RETURN
      !
END SUBROUTINE pimdnvt
  
subroutine pimd_allocation
   
   USE pimd_variables
   USE ring_io_units_module,         ONLY : iunpath
   implicit none
   integer :: i,j,ii,jj,ind
   
! WRITE(9999,*) "Alloc dims", ndimMD, natMD, nbeadMD
!FLUSH(9999)

   ! Allocation of variables
   allocate(indx(natMD))
   allocate(el(ndimMD,nbeadMD))
   allocate(rpos(ndimMD,natMD,nbeadMD))
   allocate(rpos_init(ndimMD,natMD,nbeadMD))
   allocate(forceMD(ndimMD,natMD,nbeadMD))
   allocate(forcedyn(ndimMD*natMD,nbeadMD))
   allocate(vel(ndimMD,natMD,nbeadMD))
   allocate(velocity(3,ndimMD*natMD,nbeadMD))
   allocate(pimp(ndimMD,natMD,nbeadMD))
   allocate(rtilde(natMD,ndimMD))
   allocate(vcm(ndimMD,nbeadMD))
   allocate(rcm(ndimMD,nbeadMD))
   ALLOCATE( pes(nbeadMD) )
   ALLOCATE( stress_pes_md( 6, nbeadMD ) )
   
   rpos=0.0
  



! set indexes of the different atoms for interactions
  do i=1,natMD
    indx(i)=i 
  enddo
   
  if(irun.eq.3) then
  
      allocate(cost1(nbeadMD),cost2(nbeadMD))
      allocate(tmes_bead(nbeadMD))
      allocate(mass_ion(ndimMD,natMD))
      allocate(rpos_old(ndimMD,natMD,nbeadMD),forceMD_old(ndimMD,natMD,nbeadMD))
      

      if(yesquantum) then 
        allocate(omega_mode(nbeadMD),friction_mode(nbeadMD)) 
        allocate(cmatrix(nbeadMD,nbeadMD))  
        allocate(ptilde(ndimMD,natMD,nbeadMD))
        allocate(rtilde_mode(ndimMD,natMD,nbeadMD)) 
        allocate(fbead(ndimMD,natMD),rcentroid(ndimMD,natMD))

        do i=1,nbeadMD
          omega_mode(i) = 2.d0*tfakeMD*sin((i-1)*pi/nbeadMD)
          if(i .eq. 1) then 
            friction_mode(i) = gammaMD
          else
            friction_mode(i) = 2.d0*omega_mode(i)
          endif
          if(i .eq. 1 .and. yesglobal) then
            cost1(i) = exp(-friction_mode(i)*delt) ! PILE_G thermostat
            cost2(i) = 1.d0-cost1(i)
          else
            cost1(i) = exp(-friction_mode(i)*delt*0.5d0)
            cost2(i) = sqrt(1.d0-cost1(i)**2)
          endif
          do j=1,nbeadMD
            if(j.eq. 1) then 
              cmatrix(i,j)=1.d0
            elseif (j .ge. 2 .and. j .le. nbeadMD/2) then
              cmatrix(i,j) = sqrt(2.d0)*cos(2*pi*(i-1)*(j-1)/nbeadMD)
            elseif (j .eq. nbeadMD/2+1) then
              cmatrix(i,j) = (-1.d0)**(i-1)
            else 
              cmatrix(i,j) = sqrt(2.d0)*sin(2*pi*(i-1)*(j-1)/nbeadMD)
            endif
          enddo
        enddo
  
        cmatrix=cmatrix*sqrt(1.d0/nbeadMD)

      else ! classical  
        cost1(1) = exp(-gammaMD*delt/2.d0)
        cost2(1) = sqrt(1.d0-cost1(1)**2)
      endif ! yesquantum 

      do ii=1,ndimMD
        do jj=1,natMD
          mass_ion(ii,jj)=amas(jj)
        enddo
      enddo
      
      deltahtilde=0.d0 
      
  else if(irun.eq.4) then
    
      ieskin=ndimMD*natMD
      iflagerr=0
      dt=delt
      nbead=nbeadMD
      temp=tempMD
      nion=natMD
      iscramax=6*ieskin
      iscramax=max(3*nbead,iscramax)
      normcorr=1.d0
      scalecov=1.d0

      if(nbeadMD.ne.1) then
        friction=gammaMD
      endif

! set output control in Turbo Langevin dynamics
      rankMD=1
      if(verbose) rankMD=0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
      if (.not. restart_pimd) then
         write(unit_dot_out,*) '*************************************'
         write(unit_dot_out,*) 'PIOUD integrator parameters'
         write(unit_dot_out,*) 'time step =',dt
         write(unit_dot_out,*) 'friction (gamma) =',gammaMD
         if(yessecond) write(unit_dot_out,*) 'second order Trotter beakup'
         if(nbeadMD.gt.1) then
           write(unit_dot_out,*) 'Quantum MD: Path integral Langevin!!!!'
           write(unit_dot_out,*) 'Generalized dynamics in the bead eigenmodes'
         endif
      end if
       
      allocate(psip(iscramax))
      allocate(scalpar(ieskin))
      allocate(cov_pimd(ieskin*ieskin,nbead))
      allocate(fk(natMD*ndimMD,natMD*ndimMD))

     
      if(.not.restart_pimd)  write(unit_dot_out,*) &
           '*************************************'

      if(yesquantum) then 
        temp=temp*nbead
        allocate(tmes_bead(nbead))
        tmes_bead=0.d0
        allocate(fbead(ndimMD,nion),rcentroid(ndimMD,nion))
        allocate(mass_ion(ndimMD,nion))
        write(iunpath,*) ' Mass particles (unit m_e) '
        do ii=1,ndimMD
          do jj=1,nion
             mass_ion(ii,jj)=amas(jj) 
          enddo
        enddo
        allocate(kdyn(nbead,nbead),kdyn_eig(nbead))
        kdyn=0.d0
        do ii=1,nbead
          kdyn(ii,mod(ii,nbead)+1)=-1.d0
          kdyn(mod(ii,nbead)+1,ii)=-1.d0
          kdyn(ii,ii)=2.d0
        enddo
!     kdyn=kdyn*temp**2
        kdyn=kdyn*temp**2*mass_ion(1,1)  
        do ii=1,ndimMD
          do jj=1,nion
             ind=ii+(jj-1)*ndimMD
!           scalpar(ind)=1.d0/mass_ion(ii,jj)
             scalpar(ind)=mass_ion(1,1)/mass_ion(ii,jj)
          enddo
        enddo
     
        call dsyev('V','L',nbead,kdyn,nbead,kdyn_eig,psip,3*nbead,info)                          

        write(iunpath,*) ' Eigenvalues elastic quantum term '
        do ii=1,nbead
          write(iunpath,*) ii,kdyn_eig(ii)
        enddo

      elseif(nbead.eq.1) then ! classical molecular dynamic case 

        write(iunpath,*) 'classical dynamics not implemented yet!!'
        stop
     
      endif
    
    end if
   
end subroutine pimd_allocation

subroutine pimd_deallocation
   
   use pimd_variables
   implicit none
   
   deallocate(indx)
   deallocate(amas)
   deallocate(rpos)
   deallocate(rpos_init)
   deallocate(forceMD)
   deallocate(vel)
   deallocate(pimp)
   deallocate(rtilde)
   deallocate(el)
   deallocate(vcm)
   deallocate(rcm)
   deallocate(ion_name)
   deallocate(pes)
   deallocate(stress_pes_md)
 
 
  
  if (irun .eq. 3) then
     deallocate(cost1,cost2)
     deallocate(tmes_bead)
     deallocate(mass_ion)
     deallocate(rpos_old,forceMD_old)
     
     if (yesquantum) then  
       deallocate(omega_mode,friction_mode) 
       deallocate(cmatrix)  
       deallocate(ptilde)
       deallocate(rtilde_mode) 
       deallocate(fbead,rcentroid)
     endif
  endif
  
  
  if (irun .eq. 4) then
     deallocate(psip)
     deallocate(scalpar)
     deallocate(cov_pimd)
     deallocate(fk)
     deallocate(tmes_bead)
     deallocate(mass_ion)
     if (yesquantum) then 
       deallocate(rcentroid)
       deallocate(fbead)
       deallocate(kdyn,kdyn_eig)
     endif
  end if
  
  return
  
end subroutine pimd_deallocation

SUBROUTINE pimd_gen_inputs(parse_file_name,engine_prefix,root,comm)
  !
  USE mp,            ONLY : mp_rank
  !
  IMPLICIT NONE
  !
  CHARACTER(len=*), INTENT(in) :: parse_file_name
  CHARACTER(len=*), INTENT(in) :: engine_prefix
  INTEGER, INTENT(in) :: root
  INTEGER, INTENT(in) :: comm
  !
  CHARACTER(len=512) :: dummy
  INTEGER :: i, j
  INTEGER :: parse_unit, md_unit
  INTEGER :: unit_tmp_i
  CHARACTER(len=10) :: a_tmp
  INTEGER :: myrank
  INTEGER, EXTERNAL :: myfind_free_unit
  !
  myrank =  mp_rank(comm)

  parse_unit = myfind_free_unit()
  OPEN(unit=parse_unit,file=trim(parse_file_name),status="old")
  ! ---------------------------------------------------
  ! MD INPUT PART
  ! --------------------------------------------------- 
  md_unit = myfind_free_unit()
  OPEN(unit=md_unit,file='pimd.dat',status="unknown")
  dummy=""
  DO WHILE (len_trim(dummy)<1)
     READ(parse_unit,fmt='(A512)',END=10) dummy
  ENDDO
  
  IF(trim(ADJUSTL(dummy)) == "BEGIN") THEN
     DO WHILE (trim(ADJUSTL(dummy)) /= "END")
        READ(parse_unit,*) dummy
        IF(trim(ADJUSTL(dummy)) == "BEGIN_PIMD_INPUT") THEN
           
           READ(parse_unit,'(A512)') dummy
           
           DO WHILE (trim(ADJUSTL(dummy)) /= "END_PIMD_INPUT")
              IF(myrank==root) WRITE(md_unit,'(A512)') trim(ADJUSTL(dummy))
              READ(parse_unit,'(A512)') dummy
           ENDDO
           
        ENDIF
     ENDDO
  ELSE
     CALL infomsg('path_gen_inputs','key word BEGIN missing')
  ENDIF
  CLOSE(md_unit)
  CLOSE(parse_unit)

10 CONTINUE
  !

END SUBROUTINE pimd_gen_inputs 

subroutine pimd_read_input(unit)

!!! *********************************************************************!!!
!!! initialize all the global variables - not the specific (to irun) ones!!!
!!! *********************************************************************!!!

  use pimd_variables
  USE ring_io_units_module,         ONLY : iunpath
  implicit none   
  integer :: i,j,l,k,ii,jj,ind,ntest
  integer, intent(in) :: unit
  logical :: ifex
  integer :: unit_tmp
  double precision :: dxvar
  character(len=20) :: dummy1,varname1,format_string1
  INTEGER, EXTERNAL :: myfind_free_unit
  INTEGER :: ios
  character(len=64) :: run
  logical, external :: imatches
   

! Namelists for input file
  namelist /dynamics/ restart_pimd,nbeadMD,nunitcells,nblocks,nstep_block  &
                     ,iprint,run   &
                     ,delt,tempMD,gammaMD     &
                     ,delta_force,delta_harm     &
                     ,yessecond         &
                     ,yesturboq,yesglobal,averaged_cov &
                     ,dynmat_cutoff
  
  write(iunpath,*)
  write(iunpath,'(''----------------------------------------------------------'')')
  write(iunpath,'(''      *** Ab Initio PIMD Program ***                      '')')
  write(iunpath,'('' Classical/Quantum (N,V,T) or (N,V,E) molecular dynamics  '')')
  write(iunpath,'('' with Path Integral Langevin approach                     '')')
  write(iunpath,'(''----------------------------------------------------------'')')
  write(iunpath,*)
  write(iunpath,'(''This code can be cited as follows:'')')
  write(iunpath,'(''Mouhat, F., Sorella, S., Vuilleumier, R., Saitta, A. M., & Casula, M. JCTC, 13(6), 2400-2417. (2017).'')')
  write(iunpath,*)

!    ***   System parameters    ***

! Default values defined
  restart_pimd=.false.
  ndimMD=3
  nunitcells=1
  nbeadMD=1

! set verbosity in input
  verbose=.false.
! Default values defined
  nh = .false.
  fixcm = .false.
  averaged_cov=.false.
  dynmat_cutoff=0.d0
  yessecond = .true.
  yesturboq = .true.
  yesglobal = .false. ! PILE_L is the default thermostat for irun=3
  delta_force = 1.d-4       
  delta_harm = 5.d-3

! Namelist 'dynamic' call  
  read(unit,nml=dynamics)

  delt = delt * 41.341374576 ! Convert femtoseconds to atomic unit (time)

!  ios = 0
!  IF ( ionode ) THEN
          !
!     READ( unit, dynamics, iostat = ios )
          !
!  END IF
!  CALL mp_bcast( ios, ionode_id, world_comm )

  if(imatches("verlet",run)) irun = verlet
  if(imatches("ceriotti",run)) irun = ceriotti
  if(imatches("pioud",run)) irun = pioud


  if(irun .eq. 0) then   
    nbeadMD=1
    PSMD=0.5 ! set initial value for NH
  endif 
  
  if(nbeadMD .eq. 1) then  
    yesquantum=.false.
    yesturboq=.false.   
  else
    yesquantum=.true.
  endif
 
  tempMD=tempMD/kbm1
  tfakeMD=tempMD*nbeadMD
  
  if(.not. restart_pimd) then
  
     unit_dot_out  = myfind_free_unit()
     open(unit_dot_out,file='pimd.out',form='formatted')
     unit_dot_ep   = myfind_free_unit()
     open(unit_dot_ep,file='pimd.ep',form='formatted')
     unit_dot_ep2  = myfind_free_unit()
     open(unit_dot_ep2,file='pimd.ep2',form='formatted')
     unit_dot_epsr = myfind_free_unit()
     open(unit_dot_epsr,file='pimd.epsr',form='formatted')
     unit_dot_eplr = myfind_free_unit()
     open(unit_dot_eplr,file='pimd.eplr',form='formatted')
     unit_dot_ek   = myfind_free_unit()
     open(unit_dot_ek,file='pimd.ek',form='formatted')
  
     write(unit_dot_out,'('' number of blocks          '',i10   )') nblocks
     write(unit_dot_out,'('' number of per block       '',i10   )') nstep_block
     write(unit_dot_out,'('' output frequency          '',i10   )') iprint
     write(unit_dot_out,'('' temperature               '',f10.8 )') tempMD
     write(unit_dot_out,'('' Nose-Hoover               '',l )') nh
     write(unit_dot_out,'('' Maximum single part. disp.'',f10.6 )') delt

     write(unit_dot_out,*)
     
     if(irun == 0) then
       write(unit_dot_out,*) 'Classical MD with Nosé-Hoover thermostat'
     elseif(irun == 3) then
       write(unit_dot_out,*) 'Classical/Quantum Langevin dyn. with Ceriotti integrator'
     elseif(irun == 4) then
       write(unit_dot_out,*) 'Quantum Langevin dyn. with PIOUD integrator'

     else
       write(unit_dot_out,*) ' type of run unknown: irun = ',irun
       stop
     endif
     write(unit_dot_out,*)  
     flush(unit_dot_out)
  
     unit_dot_xyz  = myfind_free_unit()
     open(unit_dot_xyz,file='pimd.xyz',form='formatted')
     unit_dot_positions   = myfind_free_unit()
     open(unit_dot_positions,file='positions.dat',form='formatted')
     unit_dot_positions_cen   = myfind_free_unit()
     open(unit_dot_positions_cen,file='positions_cen.dat',form='formatted')
     unit_dot_velocities  = myfind_free_unit()
     open(unit_dot_velocities,file='velocities.dat',form='formatted')
     unit_dot_velocities_cen  = myfind_free_unit()
     open(unit_dot_velocities_cen,file='velocities_cen.dat',form='formatted')
     unit_dot_forces = myfind_free_unit()
     open(unit_dot_forces,file='forces.dat',form='formatted')
     unit_dot_forces_cen = myfind_free_unit()
     open(unit_dot_forces_cen,file='forces_cen.dat',form='formatted')
     unit_dot_stress = myfind_free_unit()
     open(unit_dot_stress,file='stress.dat',form='formatted')
     unit_dot_stress_cen = myfind_free_unit()
     open(unit_dot_stress_cen,file='stress_cen.dat',form='formatted')
     unit_dot_localtemp = myfind_free_unit()
     open(unit_dot_localtemp,file='local_temp.dat',form='formatted')
     unit_dot_sigma   = myfind_free_unit()
     open(unit_dot_sigma,file='sigma.dat',form='formatted')
     rewind(unit_dot_xyz)
     rewind(unit_dot_positions)
     rewind(unit_dot_positions_cen)
     rewind(unit_dot_velocities)
     rewind(unit_dot_velocities_cen)
     rewind(unit_dot_forces)
     rewind(unit_dot_forces_cen)
     rewind(unit_dot_stress) !Not sure what it does Aadhityan A
     rewind(unit_dot_stress_cen)
     rewind(unit_dot_localtemp)
     rewind(unit_dot_sigma)
  
  else
     
     unit_dot_out  = myfind_free_unit()
     open(unit_dot_out,file='pimd.out',position='APPEND',form='formatted')
     unit_dot_ep  = myfind_free_unit()
     open(unit_dot_ep,file='pimd.ep',position='APPEND',form='formatted')
     unit_dot_ep2  = myfind_free_unit()
     open(unit_dot_ep2,file='pimd.ep2',position='APPEND',form='formatted')
     unit_dot_epsr  = myfind_free_unit()
     open(unit_dot_epsr,file='pimd.epsr',position='APPEND',form='formatted')
     unit_dot_eplr  = myfind_free_unit()
     open(unit_dot_eplr,file='pimd.eplr',position='APPEND',form='formatted')
     unit_dot_ek  = myfind_free_unit()
     open(unit_dot_ek,file='pimd.ek',position='APPEND',form='formatted')
     unit_dot_xyz  = myfind_free_unit()
     open(unit_dot_xyz,file='pimd.xyz',position='APPEND',form='formatted')
     
     unit_dot_positions  = myfind_free_unit()
     open(unit_dot_positions,file='positions.dat',form='formatted') !!! not append because I need to read last pos
     unit_dot_velocities  = myfind_free_unit()
     open(unit_dot_velocities,file='velocities.dat',form='formatted') !!! same as pos

     unit_dot_positions_cen  = myfind_free_unit()
     open(unit_dot_positions_cen,file='positions_cen.dat',position='APPEND',form='formatted')      
     unit_dot_velocities_cen  = myfind_free_unit()
     open(unit_dot_velocities_cen,file='velocities_cen.dat',position='APPEND',form='formatted') 
     unit_dot_forces  = myfind_free_unit()
     open(unit_dot_forces,file='forces.dat',position='APPEND',form='formatted')
     unit_dot_forces_cen  = myfind_free_unit()
     open(unit_dot_forces_cen,file='forces_cen.dat',position='APPEND',form='formatted')
     unit_dot_localtemp  = myfind_free_unit()
     open(unit_dot_localtemp,file='local_temp.dat',position='APPEND',form='formatted')
     unit_dot_sigma  = myfind_free_unit()
     open(unit_dot_sigma,file='sigma.dat',position='APPEND',form='formatted')
     unit_dot_stress = myfind_free_unit()
     open(unit_dot_stress,file='stress.dat',position='APPEND',form='formatted')
     unit_dot_stress_cen = myfind_free_unit()
     open(unit_dot_stress_cen,file='stress_cen.dat',position='APPEND',form='formatted')
  
  endif
  
!  close(unit)
     
     
  return

end subroutine pimd_read_input

subroutine pimd_setvel(iff)

  use pimd_variables, only : vcm,nbeadMD,natMD,ndimMD,vel,velocity,&
                             amas,mtot,pimp,indx,fixcm,gMD,tfakeMD

  USE random_numbers, ONLY: randy
  USE random_pioud, ONLY: pioud_randy



  implicit none
  integer :: i,k,l,iff,ind
  real(8), dimension(:), allocatable :: v2,f
  real(8) :: dummy

  allocate(v2(nbeadMD),f(nbeadMD))
      
  vcm=0.d0
  velocity=0.d0
  vel=0.d0
  v2=0.d0
  

  do k=1,nbeadMD
    do i=1,natMD
      if (iff==0) then
        do l=1,ndimMD
          dummy = pioud_randy()  
          vel(l,i,k)=dummy-0.5d0
        enddo
      endif
    enddo
  enddo

  do k=1,nbeadMD
    do i=1,natMD
      do l=1,ndimMD
        vcm(l,k)=vcm(l,k)+(amas(indx(i))*vel(l,i,k))/mtot 
      enddo
    enddo
    
    if(fixcm) then
      do i=1,natMD    
        do l=1,ndimMD
          vel(l,i,k)=vel(l,i,k)-vcm(l,k) 
        enddo
      enddo    
    endif
    
    
    do l=1,ndimMD
      do i=1,natMD
        v2(k)=v2(k)+amas(indx(i))*vel(l,i,k)**2
      enddo
    enddo    
    f(k)=sqrt(gMD*tfakeMD/v2(k))
    do i=1,natMD
      do l=1,ndimMD
        vel(l,i,k)=vel(l,i,k)*f(k)
      enddo
    enddo
    
    !!! check for the center of mass velocity
    do l=1,ndimMD
      vcm(l,k)=0.d0
      do i=1,natMD
        vcm(l,k)=vcm(l,k)+amas(indx(i))*vel(l,i,k)/mtot 
      enddo
    enddo
     
  enddo
  
  do k=1,nbeadMD
    do i=1,natMD
      ind=indx(i)
      do l=1,ndimMD
        pimp(l,i,k)=amas(ind)*vel(l,i,k)
      enddo
    enddo
  enddo
    
  deallocate(v2,f)

  return

end subroutine pimd_setvel

subroutine cumul(a,b,c,n)
      
  implicit none       
  integer :: k,i,n
  real(8),dimension(n) :: a
  real(8),dimension(2*n) :: b,c
 
  k=0
  do i=1,n
    k=k+2
  
    b(k)=(b(k)+b(k-1)**2)*c(k)+a(i)**2
    b(k-1)=b(k-1)*c(k-1)+a(i)
  
    c(k-1)=c(k-1)+1.d0
    c(k)=c(k)+1.d0
  
    b(k-1)=b(k-1)/c(k-1)
    b(k)=b(k)/c(k)-(b(k-1)**2)
  
  enddo

  return

end subroutine cumul

subroutine cumul1(a,b,c,n)

  implicit none

  integer :: n,i
  real(8) :: a
  real(8),dimension(n) :: b,c

  do i=1,n
    b(i)=b(i)*c(i)+a!(i)
    c(i)=c(i)+1.d0
    b(i)=b(i)/c(i)
  enddo

  return

end subroutine cumul1

subroutine scnd(t)

  implicit none
  real(8) :: t
  real(4) :: etime
  real(4), dimension (:), allocatable :: tarray

  allocate(tarray(2))

  t = etime(tarray)
  t = tarray(1)
   
  deallocate(tarray)

  return
end subroutine scnd

subroutine sprint(unit,jndex,dim,string)
      
  use pimd_variables, only : av,avp,anorm

  integer :: jndex,dim,k, unit, l
  real(8) :: cnorm,a1
  character(*) :: string

  cnorm=max(anorm(jndex)-1.d0,1.d0)
  k=jndex-1
 !do l=1,dim
    k=k+2
    if (av(k) == 0.d0) then
      a1=0.d0
    else
      a1=sqrt(av(k)/cnorm)
    endif
    write(unit,'(3g20.10)') avp(k-1),av(k-1),a1
  !enddo
  call flush(unit)

  return
end subroutine sprint

subroutine checkpoint(ttk)
        
  use pimd_variables, only : nbeadMD,nunitcells,unit_dot_xyz,&
                             forceMD,vel,natMD,ion_name,rcentroid,&
                             ndimMD,stress_pes_md,unit_dot_positions,&
                             unit_dot_velocities,unit_dot_forces,&
                             unit_dot_stress,&
                             unit_dot_localtemp,rpos,&
                             unit_dot_positions_cen,&
                             unit_dot_forces_cen,&
                             unit_dot_velocities_cen,&
                             unit_dot_stress_cen
  USE constants,        ONLY : ry_kbar
  implicit none
!    *******************************************************************
!    ** subroutine to write out the configurations and trajectories   **
!    *******************************************************************

  integer :: i,l,k
  real(8) :: ttk
  real(8), allocatable :: fcentroid(:,:),vcentroid(:,:),scentroid(:)
!    ********************************************************************
  allocate(fcentroid(ndimMD,natMD),vcentroid(ndimMD,natMD),scentroid(6))
  fcentroid=0.0; vcentroid=0.0 ; scentroid=0.0

  if(nbeadMD.gt.1) then


      do k=1,nbeadMD
        fcentroid(:,:)=fcentroid(:,:)+forceMD(:,:,k)
        vcentroid(:,:)=vcentroid(:,:)+vel(:,:,k)
        scentroid(:)=scentroid(:)+stress_pes_md(:,k)
        ! write(90000,*)scentroid
        ! write(90001,*)stress_pes_md
      end do
      fcentroid=fcentroid/nbeadMD; vcentroid=vcentroid/nbeadMD; scentroid=scentroid/nbeadMD


!!! -------------------write the xyz file for visualizing trajectories------     
     write(unit_dot_xyz,* ) natMD
     write(unit_dot_xyz,* ) 
     do i=1, natMD
        write(unit_dot_xyz,'(a2,3g15.5)' ) ion_name(i),(rcentroid(l,i), l=1,ndimMD)
     enddo
     flush(unit_dot_xyz)
!!!--------------------writes velocity,force and position files for postprocessing----------
!!!--------------------if the number of unit cells > 1 (crystal system)---------------------
!!!--------------------writes only the centroid coordinatMDes---------------------------------
    !  write(90002,*)scentroid
     write(unit_dot_positions_cen,'(400e15.5)') ((rcentroid(l,i),l=1,ndimMD),i=1,natMD)
     flush(unit_dot_positions_cen)
     write(unit_dot_velocities_cen,'(400e15.5)') ((vcentroid(l,i),l=1,ndimMD),i=1,natMD)
     flush(unit_dot_velocities_cen)
     write(unit_dot_forces_cen,'(400e15.5)') ((fcentroid(l,i),l=1,ndimMD),i=1,natMD)
     flush(unit_dot_forces_cen)
     write(unit_dot_stress_cen,'(400e15.5)') scentroid(:)
     flush(unit_dot_stress_cen)

    !  write(unit_dot_stress_cen,'(400e15.5)') ((scentroid(l,i),l=1,6),i=1,natMD)
    !  flush(unit_dot_stress_cen)
     
     do k=1,nbeadMD
           write(unit_dot_positions,'(400e15.5)') ((rpos(l,i,k),l=1,ndimMD),i=1,natMD)
           flush(unit_dot_positions)
           write(unit_dot_velocities,'(400e15.5)') ((vel(l,i,k),l=1,ndimMD),i=1,natMD)
           flush(unit_dot_velocities)
           write(unit_dot_forces,'(400e15.5)') ((forceMD(l,i,k),l=1,ndimMD),i=1,natMD)
           flush(unit_dot_forces)
           write(unit_dot_stress,'(400e15.5)') (stress_pes_md(l,k),l=1,6)
           flush(unit_dot_stress)
           if(k.eq.1) write(unit_dot_localtemp,'(e15.5)') ttk
           if(k.eq.1) flush(unit_dot_localtemp)
     enddo
     
  else !!! nbeadMD=1 => classical particles
     
     write (unit_dot_xyz,*) natMD
     write (unit_dot_xyz,*) 
     do i=1, natMD
        write (unit_dot_xyz,'(a2,3g15.5)' ) ion_name(i),(rpos(l,i,1), l=1,ndimMD)
     enddo
     flush(unit_dot_xyz)
     
     write(unit_dot_positions,'(1500e15.5)') ((rpos(l,i,1),l=1,ndimMD),i=1,natMD)
     flush(unit_dot_positions)
     write(unit_dot_velocities,'(1500e15.5)') ((vel(l,i,1),l=1,ndimMD),i=1,natMD)
     flush(unit_dot_velocities)
     write(unit_dot_forces,'(1500e15.5)') ((forceMD(l,i,1),l=1,ndimMD),i=1,natMD)
     flush(unit_dot_forces)
     write(unit_dot_stress,'(400e15.5)') (stress_pes_md(l,1),l=1,6)
     flush(unit_dot_stress)
     write(unit_dot_localtemp,'(e15.5)') ttk
     flush(unit_dot_localtemp)

  endif

  deallocate(fcentroid,vcentroid) 

  return
end subroutine checkpoint

subroutine zeroav(iblock)
      
  use pimd_variables, only : ipot,ikin,nav,mav,av,avp,&
                             anorm,anormp
      
  integer :: i, iblock
      
  if(iblock == 0) then
! initialization
    allocate(ipot(4))
 
    ipot(1)=1
    ipot(2)=ipot(1)+2
    ipot(3)=ipot(2)+2
    ipot(4)=ipot(3)+2
    ikin=ipot(4)+2
    nav=ikin+2
    if(nav .gt. mav) stop 'zeroav: nav.gt.mav'

    allocate(av(nav))
    allocate(avp(nav))
    allocate(anorm(nav))
    allocate(anormp(nav))
    do i=1,nav
      av(i)=0.d0
      anorm(i)=0.d0
    enddo
  else
    do i=1,nav
      avp(i)=0.d0
      anormp(i)=0.d0
    enddo
  endif

  return

end subroutine zeroav

subroutine pimd_close_files()
  
  use pimd_variables, only : unit_dot_ek, unit_dot_ep,&
                             unit_dot_ep2, unit_dot_eplr,&
                             unit_dot_epsr, unit_dot_forces,&
                             unit_dot_stress,&
                             unit_dot_localtemp, unit_dot_out,&
                             unit_dot_positions, unit_dot_sigma,&
                             unit_dot_velocities, unit_dot_xyz,&
                             unit_dot_velocities_cen,&
                             unit_dot_positions_cen,&
                             unit_dot_forces_cen

  implicit none
  
  close(unit_dot_ek)
  close(unit_dot_ep)
  close(unit_dot_ep2)
  close(unit_dot_eplr)
  close(unit_dot_epsr)
  close(unit_dot_forces)
  close(unit_dot_stress)
  close(unit_dot_forces_cen)
  close(unit_dot_localtemp)
  close(unit_dot_out)
  close(unit_dot_positions)
  close(unit_dot_positions_cen)
  close(unit_dot_sigma)
  close(unit_dot_velocities)
  close(unit_dot_velocities_cen)
  close(unit_dot_xyz)
   
end subroutine

  FUNCTION myfind_free_unit()
    !--------------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    INTEGER :: myfind_free_unit
    INTEGER :: iunit
    LOGICAL :: opnd
    !
    !
    myfind_free_unit = -1
    unit_loop: DO iunit = 999, 1, -1
       !
       INQUIRE( UNIT = iunit, OPENED = opnd )
       !
       IF ( .NOT. opnd ) THEN
          !
          myfind_free_unit = iunit
          !
          RETURN
          !
       END IF
       !
    END DO unit_loop
    !
    RETURN
    !
  END FUNCTION myfind_free_unit
