!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
subroutine aainit(lli,lqmax,mx,nlx,ap,lpx,lpl)
  !-----------------------------------------------------------------------
  !
  ! this routine computes the coefficients of the expansion of the product
  ! of two real spherical harmonics into real spherical harmonics.
  !
  !     Y_limi(r) * Y_ljmj(r) = \sum_LM  ap(LM,limi,ljmj)  Y_LM(r)
  !
  ! On output:
  ! ap     the expansion coefficients
  ! lpx    for each input limi,ljmj is the number of LM in the sum
  ! lpl    for each input limi,ljmj points to the allowed LM
  !
  ! The indices limi,ljmj and LM assume the order for real spherical harmonics
  ! given in routine ylmr2
  !
  use parameters, only : DP
  implicit none
  !
  ! first the I/O variables
  !
  integer :: &
       lli,            &! input: the maximum li considered
       lqmax,          &! input: array dimension
       mx,             &! input: array dimension
       nlx,            &! input: array dimension, must be >= lli**2
       lpx(nlx,nlx),   &! output: maximum number of LM for limi,ljmj
       lpl(nlx,nlx,mx)  ! output: counter on LM couples

  real(kind=DP) :: &
       ap(lqmax*lqmax,nlx,nlx)    !  output: the expansion coefficients
  !
  ! here the local variables
  !
  integer :: llx, l, li, lj

  real(kind=DP) , allocatable :: r(:,:), rr(:), ylm(:,:), mly(:,:)
  ! an array of random vectors: r(3,llx)
  ! the norm of r: rr(llx)
  ! the real spherical harmonics for array r: ylm(llx,llx)
  ! the inverse of ylm considered as a matrix: mly(llx,llx)

  real(kind=DP) :: compute_ap ! a function computing ap

  if (lli < 0) call errore('aainit','lli not allowed',lli)

  if (lli*lli .gt. nlx) call errore('aainit','nlx is too small ',lli*lli)

  llx = (2*lli-1)**2
  if (2*lli-1 > lqmax) &
      call errore('aainit','ap leading dimension is too small',llx)

  allocate (r( 3, llx ))    
  allocate (rr( llx ))    
  allocate (ylm( llx, llx ))    
  allocate (mly( llx, llx ))    

  r(:,:)   = 0.d0
  ylm(:,:) = 0.d0
  mly(:,:) = 0.d0
  ap(:,:,:)= 0.d0

  ! - generate an array of random vectors (uniform deviate on unitary sphere)
  call gen_rndm_r(llx,r,rr)
  ! - generate the real spherical harmonics for the array: ylm(ir,lm)
  call ylmr2(llx,llx,r,rr,ylm)
  !-  store the inverse of ylm(ir,lm) in mly(lm,ir)
  call invmat(ylm,mly,llx)
  !-  for each li,lj compute ap(l,li,lj) and the indices, lpx and lpl
  do li = 1, lli*lli
     do lj = 1, lli*lli
        lpx(li,lj)=0
        do l = 1, llx
           ap(l,li,lj) = compute_ap(l,li,lj,llx,ylm,mly)
           if (abs(ap(l,li,lj)).gt.1.d-3) then
              lpx(li,lj) = lpx(li,lj) + 1
              if (lpx(li,lj).gt.mx) &
                   call errore('aainit','mx dimension too small', lpx(li,lj))
              lpl(li,lj,lpx(li,lj)) = l
           end if
        end do
     end do
  end do

  deallocate(mly)
  deallocate(ylm)
  deallocate(rr)
  deallocate(r)

  return
end subroutine aainit
!
!-----------------------------------------------------------------------
subroutine gen_rndm_r(llx,r,rr)
  !-----------------------------------------------------------------------
  ! - generate an array of random vectors (uniform deviate on unitary sphere)
  use parameters, only : DP
  implicit none
  !
  ! first the I/O variables
  !
  integer :: llx         ! input: the dimension of r and rr

  real(kind=DP) :: &
       r(3,llx),  &! output: an array of random vectors
       rr(llx)    ! output: the norm of r
  !
  ! here the local variables
  !
  integer :: ir
  real(kind=DP) :: costheta, sintheta, phi, rndm
  !
  ! a parameter
  !
  real(kind=DP) :: tpi ! two times pi
  parameter (tpi = 2.d0 * 3.14159265358979d0)

  external rndm

  do ir = 1, llx
     costheta = 2.d0 * rndm() - 1.d0
     sintheta = sqrt ( 1.d0 - costheta*costheta)
     phi = tpi * rndm()
     r (1,ir) = sintheta * cos(phi)
     r (2,ir) = sintheta * sin(phi)
     r (3,ir) = costheta
     rr(ir)   = 1.d0
  end do

  return
end subroutine gen_rndm_r
!
!-----------------------------------------------------------------------
function compute_ap(l,li,lj,llx,ylm,mly)
  !-----------------------------------------------------------------------
  !-  given an l and a li,lj pair compute ap(l,li,lj)
  use parameters, only : DP
  implicit none
  !
  ! first the I/O variables
  !
  integer :: &
       llx,         &! the dimension of ylm and mly
       l,li,lj       ! the arguments of the array ap

  real(kind=DP) :: &
       compute_ap,  &! this function
       ylm(llx,llx),&! the real spherical harmonics for array r
       mly(llx,llx)  ! the inverse of ylm considered as a matrix
  !
  ! here the local variables
  !
  integer :: ir

  compute_ap = 0.d0
  do ir = 1,llx
     compute_ap = compute_ap + mly(l,ir)*ylm(ir,li)*ylm(ir,lj)
  end do

  return
end function compute_ap
