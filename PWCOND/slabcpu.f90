!
! Copyright (C) 2003 A. Smogunov 
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
subroutine slabcpu(nrz, nrzp, nkofz, bdl1, bdl2, bdr1, bdr2, z)
!
! To set up the division of slabs among CPU
!
USE kinds, only : DP
implicit none  
  integer :: k, nrz, nrzp, kin, kfin, startk, lastk 
  integer :: nkofz(nrz)  
  real(kind=DP), parameter :: eps=1.d-8
  real(kind=DP) :: bdl1, bdl2, bdr1, bdr2, z(nrz) 

  do k=1, nrz
    nkofz(k)=0
  enddo       

! for the left tip
  do k=1, nrz
    if (z(k).le.bdl1+eps)  kin=k
    if (z(k).le.bdl2-eps)  kfin=k
  enddo            
  if (kfin-kin.gt.0) then
    call divide(kfin-kin+1,startk,lastk)
    startk=startk+kin-1
    lastk=lastk+kin-1
    nrzp=lastk-startk+1
    do k=startk, lastk
      nkofz(k)=k-startk+1
    enddo
  endif

! for the scattering region
  do k=1, nrz
    if (z(k).le.bdl2+eps)  kin=k
    if (z(k).le.bdr1-eps)  kfin=k
  enddo     
  if (kfin-kin.gt.0) then
    call divide(kfin-kin+1,startk,lastk)
    startk=startk+kin-1
    lastk=lastk+kin-1
    do k=startk, lastk
       nkofz(k)=nrzp+k-startk+1
    enddo                           
    nrzp=nrzp+lastk-startk+1
  endif

! for the right tip
  do k=1, nrz
    if (z(k).le.bdr1+eps)  kin=k
    if (z(k).le.bdr2-eps)  kfin=k
  enddo
  if (kfin-kin.gt.0) then
    call divide(kfin-kin+1,startk,lastk)
    startk=startk+kin-1
    lastk=lastk+kin-1
    do k=startk, lastk
      nkofz(k)=nrzp+k-startk+1
    enddo
    nrzp=nrzp+lastk-startk+1    
  endif

  return
end subroutine slabcpu

