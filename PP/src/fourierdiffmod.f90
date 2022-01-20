!
! Copyright (C) 2001-2022 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
! Author: Ivan Carnimeo (September 2021)
!----------------------------------------------------------------------------
!
!----------------------------------------------------------------------------
MODULE fourierdiffmod
use fouriermod 
implicit none
CONTAINS
!----------------------------------------------------------------------------
subroutine fourierdiff(Nb, Nq, q, eq, Nk, k, ek, Nsym, at, bg, Op)
!
! compute the band structure with Fourier interpolation from
! Pickett W. E., Krakauer H., Allen P. B., Phys. Rev. B, vol. 38, issue 4, page 2721, 1988 
!
implicit none
  integer,  intent(in) :: Nq, Nk, Nb, Nsym
  real(dp), intent(in) :: q(3,Nq), k(3,Nk), eq(Nq,Nb)
  real(dp), intent(in) :: at(3,3) 
  real(dp), intent(in) :: bg(3,3) 
  real(dp), intent(in) :: Op(1:3,1:3,1:Nsym)
  real(dp), intent(out) :: ek(Nk,Nb)
  !
  ! local variables
  !
  integer :: Na, ib, ik
  complex(dp), allocatable :: fStarsOnQ(:,:) ! Star functions at uniform q-points
  complex(dp), allocatable :: fStarsOnK(:,:) ! Star functions at path of k-points 
  complex(dp), allocatable :: matQQ(:,:)     ! this is exactly H in the reference article
  complex(dp), allocatable :: matKQ(:,:)     ! this is an intermediate quantity S_m(q)*S_m(k)/rho(R_m) to construct J
  complex(dp), allocatable :: matJ(:,:)      ! this is exactly J in the reference article
  complex(dp), allocatable :: ek_c(:,:), eq_c(:,:) ! complex band energies (for ZGEMM)
  real(dp) :: vec(3)
  complex(dp) :: fStar
  !
  ! for matrix inversion
  complex(dp), allocatable :: matQQ_(:,:)     
  real(dp), allocatable :: rmatQQ(:,:)     
  real(dp), allocatable :: rmatQQ_(:,:)     
  real(dp), allocatable :: rmatJ(:,:)     
  ! 
  ! for linear system solution
  integer :: INFO, istar, NCoeff
  integer, allocatable :: IPIV(:)
  real(dp) :: sqrtrho 
  complex(dp), allocatable :: matA(:,:), matX(:,:), matB(:,:)
  complex(dp), allocatable :: matS(:)
  complex(dp), allocatable :: matC(:,:)  ! coefficients for m= 2 ,..., NStars
  complex(dp), allocatable :: matC1(:)   ! coefficient for m=1 are treated separately 
  !
  write(*,'(A)') '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
  write(*,'(A)') 'Fourier difference interpolation method'
  if(check_periodicity) write(*,*) 'Checking Star functions periodicity (WARNING: time consuming)' 
  !
  Na = Nq - 1  ! dimension of the linear system 
  !
  write(*,'(A)') 'Creating b = e(q_i) - e(q_Nq)   i = 1, ... , Nq-1 '
  allocate( matB(Na,Nb), matX(Na,Nb) )
  do ib = 1, Nb
    matB(1:Na,ib) = (One, Zero) * ( eq(1:Na,ib) - eq(Nq,ib) )
  end do 
  matX = matB ! matX will be overwritten with the solution of Ax=b
  !
  write(*,'(A)') 'Creating A'
  write(*,'(A)') 'Creating Star functions...'
  Call find_stars(NSym, Op, at, .true.) 
  !
  if(check_periodicity) then 
    write(*,*) 'Checking Star functions periodicity...'
    Call check_stars(Nq, q, NSym, Op, bg) 
  end if 
  !
  ! fStarsOnQ = [S_m(q_i)-S_m(q_Nq)] / sqrt(rho_m)
  write(*,*) 'Computing fStarsOnQ...'
  allocate( fStarsOnQ(Na,NStars), matS(NStars) )
  fStarsOnQ = (Zero, Zero)
  Call compute_stars(fStarsOnQ, Na, Nq, q, NSym, Op, 2, .true., matS) 
  !
  ! matQQ = fStarsOnQ * fStarsOnQ^T  = sum_m [S_m(q_i)-S_m(q_Nq)]*[S_m(q_j)-S_m(q_Nq)] / rho_m
  write(*,*) 'Computing fStarsOnQ * fStarsOnQ*...'
  allocate(matQQ(Na,Na), matA(Na,Na))
  matQQ = (Zero, Zero)
  Call ZGEMM('N', 'C', Na, Na, NStars, (One,Zero), fStarsOnQ, Na, fStarsOnQ, Na, (Zero,Zero), matQQ, Na)
  matA = matQQ
  !
  write(*,'(A)') 'Solving Ax = b'
  allocate( IPIV(Na) )
  Call ZGESV(Na, Nb, matQQ, Na, IPIV, matX, Na, INFO) 
  deallocate(IPIV) 
  write(*,'(A)') 'Checking A*x - b = 0...'
  Call ZGEMM('N', 'N', Na, Nb, Na, (One,Zero), matA, Na, matX, Na, -(One,Zero), matB, Na)
  Call MatCheck_k('A*x - b = 0', matB, Na, Nb)
  !  
  ! C_m,ib = rho^(-1)_m sum_m=2  lambda_iq,ib * [S_m(q_i)-S_m(q_Nq)] m = 2, ... NStars
  write(*,*) 'Computing coefficients...'  
  allocate( matC(NStars,Nb), matC1(Nb) ) 
  matC = (Zero, Zero)
  Call ZGEMM('C', 'N', NStars, Nb, Na, (One, Zero), fStarsOnQ, Na, matX, Na, (Zero, Zero), matC, NStars)
  do istar = 1, NStars 
    sqrtrho = sqrt_rho(VecStars(:,istar))
    matC(istar,1:Nb) = matC(istar,1:Nb)/sqrtrho ! now matC has the right coefficients for m = 2, ... 
  end do
  do ib = 1, Nb
    matC1(ib) = eq(Nq,ib) - dot_product(matC(:,ib),matS(:))
  end do 
  !  
  write(*,*) 'Computing bands...'  
  ! fStarsOnK = S_m(k) 
  write(*,*) 'Computing fStarsOnK...'
  allocate( fStarsOnK(Nk,NStars), ek_c(Nk,Nb) )
  fStarsOnK = (Zero, Zero)
  Call compute_stars(fStarsOnK, Nk, Nk, k, NSym, Op, 0) 
  ek_c = (Zero, Zero)
  Call ZGEMM('N', 'N', Nk, Nb, NStars, (One, Zero), fStarsOnK, Nk, matC, NStars, (Zero, Zero), ek_c, Nk)
  !
  ! now add the C1 coefficient
  vec(1:3) = Zero
  do ik = 1, Nk
    fStar = star_function(0, k(1:3,ik), vec, NSym, Op)
    do ib = 1, Nb
      ek_c(ik, ib) = ek_c(ik, ib) + matC1(ib) * fStar  
    end do 
  end do 
  !
  ek(:,:) = dble(ek_c(:,:))
  !
  deallocate( matX, matB, matA, matC, matC1, ek_c )
  deallocate( RoughC ) 
  deallocate( VecStars )
  deallocate( fStarsOnQ )
  deallocate( fStarsOnK  )
  !
  return
  !
end subroutine fourierdiff
!----------------------------------------------------------------------------
END MODULE
!----------------------------------------------------------------------------
