! FOR GWW
!
! Author: P. Umari
!
!-------------------------
SUBROUTINE go_wannier_product(numpw,pmat,rot_u, tresh, maxiter)
!-------------------------
!this routine find the localization rotation rot_u
!given the exp(igX) factors pmat
! #ifdef __GWW

  USE kinds,    ONLY : DP
  USE cell_base, ONLY: at, alat, tpiba, omega, tpiba2
  USE constants, ONLY : pi


  implicit none

  INTEGER, INTENT(in) :: numpw ! dimension of space
  REAL(kind=DP), INTENT(in) :: tresh! treshold on wannier wfcs'spread
  INTEGER, INTENT(in) :: maxiter!max number of iterations
  COMPLEX(kind=DP),INTENT(in) :: pmat(3,numpw,numpw)!exp(iGX) factors
  REAL(kind=DP), INTENT(out) :: rot_u(numpw,numpw)!rotation matrix to be found

  !  --- Internal definitions ---
  INTEGER :: i,j,k,l,it
  REAL(kind=DP), ALLOCATABLE  :: matsincos(:,:,:)  !respective sin and cos matrices
  REAL(kind=DP) :: w(6) ! weights for orthorombic cells

  REAL(kind=DP) :: theta, theta4, d2, aa(2), cc, ss, tempi, tempj
  REAL(kind=DP) :: omg0,omg1!omega for testi convergence

  ALLOCATE( matsincos(numpw,numpw,6))


! set weights

  do i=1,3
     w(i)=((at(i,i)*alat)/pi)
     w(i+3)=((at(i,i)*alat)/pi)
  enddo

 ! calculates matsincos

  do i=1,3
     matsincos(1:numpw,1:numpw,i)=w(i)*real(pmat(i,1:numpw,1:numpw))
     matsincos(1:numpw,1:numpw,i+3)=w(i)*aimag(pmat(i,1:numpw,1:numpw))
  enddo

!----------------Young'su stuff-----

! Scale thr according to # of Wannier orbitals and cell length

!  set initial rotation matrix

  rot_u(:,:)=0.d0
  do i=1,numpw
     rot_u(i,i)=1.d0
  end do

! now  valence subspace


! calculate omega
  omg0  = 0.d0
  do k=1,6
     do i=1,numpw
        omg0=omg0 + matsincos(i,i,k)*matsincos(i,i,k)
     enddo
  enddo

  write(6,*) 'LOCALIZING WANNIER FUNCTIONS:'


! Start Iteration =====================================================


  do it=1,maxiter
     do i=1,numpw
        do j=i+1,numpw

! Construct aa
           aa(:)=0.d0
           do k=1,6
              aa(1)=aa(1)+matsincos(i,j,k)*(matsincos(i,i,k)-matsincos(j,j,k))
              aa(2)=aa(2)+matsincos(i,j,k)*matsincos(i,j,k)-0.25d0*(matsincos(i,i,k)-&
              matsincos(j,j,k))*(matsincos(i,i,k)-matsincos(j,j,k))
           end do
           if(abs(aa(2)).gt.1.d-10) then
              theta4=-aa(1)/aa(2)
              theta=0.25*atan(theta4)
           elseif (abs(aa(1)).lt.1.d-10) then
              theta=0.d0
              aa(2)=0.d0
           else
              theta=pi/4.d0
           endif
           d2=aa(1)*sin(4.*theta)-aa(2)*cos(4.*theta)
           if(d2.le.0.d0) theta=theta+pi/4.d0
!
           cc=cos(theta)
           ss=sin(theta)

! update overlap matrices
           do l=1,6
              ! AR
              do k=1,numpw
                 tempi=matsincos(k,i,l)*cc+matsincos(k,j,l)*ss
                 tempj=-matsincos(k,i,l)*ss+matsincos(k,j,l)*cc
                 matsincos(k,i,l)=tempi
                 matsincos(k,j,l)=tempj
              end do
! R^+ A R
              do k=1,numpw
                 tempi=cc*matsincos(i,k,l)+ss*matsincos(j,k,l)
                 tempj=-ss*matsincos(i,k,l)+cc*matsincos(j,k,l)
                 matsincos(i,k,l)=tempi
                 matsincos(j,k,l)=tempj
              end do
           end do

! update U : U=UR
           do k=1,numpw
              tempi=rot_u(k,i)*cc+rot_u(k,j)*ss
              tempj=-rot_u(k,i)*ss+rot_u(k,j)*cc
              rot_u(k,i)=tempi
              rot_u(k,j)=tempj
           end do
!
        end do
     end do

     omg1=0.d0
     do k=1,6
        do i=1,numpw
           omg1=omg1+matsincos(i,i,k)*matsincos(i,i,k)
        end do
     end do


     if(abs(omg1-omg0).lt.tresh ) EXIT
     write(6,*) 'go_wannier_product spread:', omg1,theta
     omg0=omg1

  end do

  deallocate(matsincos)
! #endif __GWW
  return
END SUBROUTINE go_wannier_product

