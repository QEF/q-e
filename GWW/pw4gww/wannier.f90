!
! Copyright (C) 2001-2013 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!

!-------------------------
subroutine go_wannier( iun_wannier, tresh, maxiter,nbndv, itask, ispin)
!-------------------------
! this routine read the wfcs from iun_wannier then
! transfrom to real wfcs, transform to wannier functions
! with treshold tresh and maxiter 
! using Gygi scheme
! and writes wfcs on iun_wannier and writes wannier_centers file
! it's possible to localized two different subspaces


  USE kinds,    ONLY : DP
  USE us
  USE wvfct,    ONLY : nbnd
  USE wavefunctions_module, ONLY : evc
  USE gvect
  USE basis
  USE klist
  USE constants, ONLY : e2, pi, tpi, fpi
  USE io_files, ONLY: nwordwfc
  USE io_global, ONLY : stdout
  USE cell_base, ONLY: at, alat, tpiba, omega, tpiba2
  USE wannier_gw
  USE mp, ONLY : mp_sum, mp_bcast, mp_barrier
  USE mp_world,  ONLY : mpime, nproc
  USE control_flags,  ONLY : gamma_only
  
  implicit none

  INTEGER, INTENT(in) :: iun_wannier !units for reading wfc
  REAL(kind=DP), INTENT(in) :: tresh! treshold on wannier wfcs'spread
  INTEGER, INTENT(in) :: maxiter!max number of iterations
  INTEGER, INTENT(in) :: nbndv !number of first bands wich are localized separately
  INTEGER,INTENT(in)  :: itask! if == 1 calculate {C'} subspace
  INTEGER, INTENT(in) :: ispin!spin channel


  !  --- Internal definitions ---
  INTEGER :: i,j,k,l,it,iw, ii
  COMPLEX(kind=DP), ALLOCATABLE  :: matgp(:,:,:)! for matrix <Psi_i|{exp(iG{x,y,z})|Psi_j>
  REAL(kind=DP), ALLOCATABLE  :: matsincos(:,:,:)  !respective sin and cos matrices
  REAL(kind=DP), ALLOCATABLE  ::rot_u(:,:) ! unitarian matrix which transorm to wannier wfcs
  REAL(kind=DP) :: w(6) ! weights for orthorombic cells
  COMPLEX(kind=DP), ALLOCATABLE :: crot_u(:,:),crot_u_tmp(:,:)
  INTEGER nbnd_start, nbnd_end!for conduction states manifolds

  REAL(kind=DP) :: theta, theta4, d2, aa(2), cc, ss, tempi, tempj
  REAL(kind=DP) :: omg0,omg1!omega for testi convergence
  COMPLEX(kind=DP) :: sca
 
  REAL(kind=DP), ALLOCATABLE :: tmp_mat(:,:)
  INTEGER :: nbnd_second

  INTEGER :: n_oper

  INTEGER :: na,nb
  REAL(kind=DP) :: delta
  REAL(kind=DP), ALLOCATABLE :: vtmpi(:),vtmpj(:)
  REAL(Kind=DP), ALLOCATABLE :: ene_tmp(:)
  INTEGER :: nbnd_par,nbnd_start_me,nbnd_end_me,nbnd_start_proc,nbnd_end_proc,ip
  REAL(kind=DP), ALLOCATABLE :: matsincos_par(:,:,:)
  INTEGER :: n_oper_par,oper_start,oper_end
  REAL(kind=DP) :: sca_mat(6)

!  ALLOCATE( matgp(nbnd,nbnd,3))
  
  ALLOCATE( matsincos(nbnd,nbnd,6))
     

  



!transfrom to real wfcs:
  
  write(stdout,*) 'Transform to real wfcs'
  FLUSH(stdout)

  if(nbndv<0 .OR. nbndv > nbnd) call errore('go_wannier','nbndv: illegal value',1)
  !if(mod(nbndv,2) /= 0 ) call errore('go_wannier','nbndv, odd',1) 

  if(.not.gamma_only) then
    !call real_wfc(u_trans, 1, iun_wannier,nbndv)
  else
  !writes KS wavefunctions on file
    u_trans(:,:,ispin)=(0.d0,0.d0)!it could also be defined as real
    do iw=1,nbnd
      u_trans(iw,iw,ispin)=(1.d0,0.d0)
    enddo
  endif

! set matg p,m matrices
  if(gamma_only) then
    !call matrix_wannier_gamma_big(matgp,1,nset,itask)
     call matrix_wannier_gamma_big(matsincos,ispin,nset,itask)
  else
     write(stdout,*)'ONLY GAMMA POINT IMPLEMENTED'
     stop
  endif

  write(stdout,*) 'Out of matrix_wannier_gamma_big'
  FLUSH(stdout)


! set weights 

  do i=1,3
     w(i)=((at(i,i)*alat)/pi)
     w(i+3)=((at(i,i)*alat)/pi)
  enddo

 ! calculates matsincos

  do i=1,3
     matsincos(1:nbnd,1:nbnd,i)=w(i)*matsincos(1:nbnd,1:nbnd,i)
     matsincos(1:nbnd,1:nbnd,i+3)=w(i)*matsincos(1:nbnd,1:nbnd,i+3)
   !  matsincos(1:nbnd,1:nbnd,i)=w(i)*real(matgp(1:nbnd,1:nbnd,i))
   !  matsincos(1:nbnd,1:nbnd,i+3)=w(i)*aimag(matgp(1:nbnd,1:nbnd,i))
  enddo

 ! deallocate(matgp)


  


  n_oper=6
  

!----------------Young'su stuff-----




! Scale thr according to # of Wannier orbitals and cell length

!  set initial rotation matrix

  ALLOCATE (rot_u(nbnd,nbnd))
  rot_u(:,:)=0.d0
  do i=1,nbnd
     rot_u(i,i)=1.d0
  end do

! now  valence subspace


! calculate omega
  omg0  = 0.d0
  do k=1,6
     do i=1,nbndv
        omg0=omg0 + matsincos(i,i,k)*matsincos(i,i,k)
     enddo
  enddo
  
  write(stdout ,*) 'LOCALIZING WANNIER FUNCTIONS:'
  FLUSH(stdout)
  

! Start Iteration =====================================================


  do it=1,maxiter
     do i=1,nbndv
        do j=i+1,nbndv
           
! Construct aa
           aa(:)=0.d0
           do k=1,n_oper
              aa(1)=aa(1)+matsincos(i,j,k)*(matsincos(i,i,k)-matsincos(j,j,k))
              aa(2)=aa(2)+matsincos(i,j,k)*matsincos(i,j,k)-&
    &0.25d0*(matsincos(i,i,k)-matsincos(j,j,k))*(matsincos(i,i,k)-matsincos(j,j,k))
           end do
           if(abs(aa(2)).gt.1.d-14) then
              theta4=-aa(1)/aa(2)
              theta=0.25*atan(theta4)
           elseif (abs(aa(1)).lt.1.d-14) then
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
           do l=1,n_oper
              ! AR
              do k=1,nbndv
                 tempi=matsincos(k,i,l)*cc+matsincos(k,j,l)*ss
                 tempj=-matsincos(k,i,l)*ss+matsincos(k,j,l)*cc
                 matsincos(k,i,l)=tempi
                 matsincos(k,j,l)=tempj
              end do
! R^+ A R
              do k=1,nbndv
                 tempi=cc*matsincos(i,k,l)+ss*matsincos(j,k,l)
                 tempj=-ss*matsincos(i,k,l)+cc*matsincos(j,k,l)
                 matsincos(i,k,l)=tempi
                 matsincos(j,k,l)=tempj
              end do
           end do

! update U : U=UR
           do k=1,nbndv
              tempi=rot_u(k,i)*cc+rot_u(k,j)*ss
              tempj=-rot_u(k,i)*ss+rot_u(k,j)*cc
              rot_u(k,i)=tempi
              rot_u(k,j)=tempj
           end do
!
        end do
     end do
     
     omg1=0.d0
     do k=1,n_oper
        do i=1,nbndv
           omg1=omg1+matsincos(i,i,k)*matsincos(i,i,k)
        end do
     end do
   
 
     write(stdout,*) 'Spread', omg1,omg0
     FLUSH(stdout)
     if(abs(omg1-omg0).lt.tresh ) EXIT
     omg0=omg1

  end do

!centers of wanniers

  do i=1,nbndv
     do k=1,3
        wannier_centers(k,i,1)=aimag(log(matsincos(i,i,k)+(0.d0,1.d0)*matsincos(i,i,k+3)))
        wannier_centers(k,i,1)=wannier_centers(k,i,1)*at(k,k)/2.d0/pi
        if(wannier_centers(k,i,1) < 0.d0) wannier_centers(k,i,1)=at(k,k)+wannier_centers(k,i,1)
     end do
  enddo
  do i=1,nbndv
     write(stdout,*) 'Center Wannier:', wannier_centers(1,i,1)*alat,wannier_centers(3,i,1)*alat,wannier_centers(3,i,1)*alat
  enddo
  FLUSH(stdout)
!---now conduction subspace



    
 
  deallocate(matsincos)




 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!------ rotate wfc

  if(.not.gamma_only) then
    !call rotate_wannier(rot_u,1,iun_wannier)
  else
    call rotate_wannier_gamma(rot_u,1,0)
  endif



 




!------------update u_trans

!crot_u and u_trans are complex, rot_u is real

  allocate(crot_u(nbnd,nbnd),crot_u_tmp(nbnd,nbnd))
 

  crot_u_tmp(:,:)=dcmplx(rot_u(:,:),0.d0)
  call zgemm('T','N',nbnd,nbnd,nbnd,(1.d0,0.d0),crot_u_tmp,nbnd,u_trans(1,1,ispin),nbnd,&
       &(0.d0,0.d0),crot_u,nbnd)


  u_trans(1:nbnd,1:nbnd,ispin)=crot_u(1:nbnd,1:nbnd)


! 
  deallocate(rot_u)
  deallocate(crot_u,crot_u_tmp)


 return

end subroutine go_wannier



function fast_sin(theta,na,nb,table_sin_a, table_sin_b, table_cos_a, table_cos_b)
!this routine calculates fast the sinus
!using sin(alpha+beta)=sin(alpha)*cos(beta)+sin(beta)cos(alpha)

   USE kinds,       ONLY : DP
   USE constants,   ONLY : pi,tpi

   implicit none

   REAL(kind=DP) :: fast_sin!the calculated value for the sin
   REAL(kind=DP), INTENT(in) :: theta!the input angle in radiants
   INTEGER :: na!number of elements in the table (pi/2)/Na  from 0 to pi/2
   INTEGER :: nb!number of elements in the table (pi/2)/(Na*Nb) for 0 to (pi/2)/(Na)
   REAL(kind=DP) :: table_sin_a(na)!tabel of sin (pi/2)/Na  from 0 to pi/2
   REAL(kind=DP) :: table_sin_b(nb)!tabel of sin (pi/2)/Na  from 0 to pi/2
   REAL(kind=DP) :: table_cos_a(na)!tabel of cos (pi/2)/Na  from 0 to pi/2
   REAL(kind=DP) :: table_cos_b(nb)!tabel of cos (pi/2)/Na  from 0 to pi/2

   REAL(kind=DP) :: angle, angle_p
   INTEGER :: ia,ib
   REAL(kind=DP) :: sign, da, db

!find angle from 0 to pi/2

   angle_p= theta - real(floor(theta/tpi))*tpi
   if( angle_p <= (pi/2.d0)) then
      angle=angle_p
      sign=1.d0
   else if( angle_p <= pi) then
      sign=1.d0
      angle=pi-angle_p
   else if(angle_p <= 3.d0*pi/2.d0) then
      sign=-1.d0
      angle=angle_p-pi
   else
      sign=-1.d0
      angle=tpi-angle_p
   endif

!determines ia, ib

   da=pi/(2.d0*real(na))
   db=da/real(nb)
   ia=floor(angle/(da))
   ib=floor((angle-real(ia)*da)/db)
   ia=ia+1
   ib=ib+1

!sin(a+b)
   fast_sin=sign*(table_sin_a(ia)*table_cos_b(ib)+table_sin_b(ib)*table_cos_a(ia))

   return



 end function fast_sin



function fast_cos(theta,na,nb,table_sin_a, table_sin_b, table_cos_a, table_cos_b)
!this routine calculates fast the sinus
!using cos(alpha+beta)=cos(alpha)*cos(beta)-sin(alpha)*sin(beta)

   USE kinds,       ONLY : DP
   USE constants,   ONLY : pi,tpi

   implicit none

   REAL(kind=DP) :: fast_cos!the calculated value for the cos
   REAL(kind=DP), INTENT(in) :: theta!the input angle in radiants
   INTEGER :: na!number of elements in the table (pi/2)/Na  from 0 to pi/2
   INTEGER :: nb!number of elements in the table (pi/2)/(Na*Nb) for 0 to (pi/2)/(Na)
   REAL(kind=DP) :: table_sin_a(na)!tabel of sin (pi/2)/Na  from 0 to pi/2
   REAL(kind=DP) :: table_sin_b(nb)!tabel of sin (pi/2)/Na  from 0 to pi/2
   REAL(kind=DP) :: table_cos_a(na)!tabel of cos (pi/2)/Na  from 0 to pi/2
   REAL(kind=DP) :: table_cos_b(nb)!tabel of cos (pi/2)/Na  from 0 to pi/2

   REAL(kind=DP) :: angle, angle_p
   INTEGER :: ia,ib
   REAL(kind=DP) :: sign, da, db

!find angle from 0 to pi/2

   angle_p= theta - real(floor(theta/tpi))*tpi
   if( angle_p <= (pi/2.d0)) then
      angle=angle_p
      sign=1.d0
   else if( angle_p <= pi) then
      sign=-1.d0
      angle=pi-angle_p
   else if(angle_p <= 3.d0*pi/2.d0) then
      sign=-1.d0
      angle=angle_p-pi
   else
      sign=1.d0
      angle=tpi-angle_p
   endif

!determines ia, ib

   da=pi/(2.d0*real(na))
   db=da/real(nb)
   ia=floor(angle/(da))
   ib=floor((angle-real(ia)*da)/db)

   ia=ia+1
   ib=ib+1

!cos(a+b)
   fast_cos=sign*(table_cos_a(ia)*table_cos_b(ib)-table_sin_a(ia)*table_sin_b(ib))

   return



 end function fast_cos


 function fast_atan(tg,na,nb,table_sin_a, table_sin_b, table_cos_a, table_cos_b)
!this routine calculates fast the arctan using a bisection algorithm

   USE kinds,       ONLY : DP
   USE constants,   ONLY : pi

   implicit none

   REAL(kind=DP) :: fast_atan!the calculated value for the atan
   REAL(kind=DP), INTENT(in) :: tg!the input tan
   INTEGER :: na!number of elements in the table (pi/2)/Na  from 0 to pi/2
   INTEGER :: nb!number of elements in the table (pi/2)/(Na*Nb) for 0 to (pi/2)/(Na)
   REAL(kind=DP) :: table_sin_a(na)!tabel of sin (pi/2)/Na  from 0 to pi/2
   REAL(kind=DP) :: table_sin_b(nb)!tabel of sin (pi/2)/Na  from 0 to pi/2
   REAL(kind=DP) :: table_cos_a(na)!tabel of cos (pi/2)/Na  from 0 to pi/2
   REAL(kind=DP) :: table_cos_b(nb)!tabel of cos (pi/2)/Na  from 0 to pi/2
   
   INTEGER, PARAMETER :: n=20!number of bisections

   REAL(kind=DP) :: sign, tang
   INTEGER :: i
   REAL(kind=DP) :: ang, delta
   REAL(kind=DP) :: tan_try
   REAL(kind=DP), EXTERNAL :: fast_sin, fast_cos


   if(tg >= 0.d0) then
      sign=1.d0
      tang=tg 
   else
      sign=-1.d0
      tang=-tg
   endif

   
   ang=pi/4.d0
   delta=pi/4.d0
   do i=1,n
      delta=delta/2.d0
      tan_try=fast_sin(ang,na,nb,table_sin_a, table_sin_b, table_cos_a, table_cos_b)/&
               &fast_cos(ang,na,nb,table_sin_a, table_sin_b, table_cos_a, table_cos_b)
      if(tang >= tan_try) then
         ang=ang+delta
      else
         ang=ang-delta
      endif
   enddo
   
   fast_atan=sign*ang

   return
 end function fast_atan
