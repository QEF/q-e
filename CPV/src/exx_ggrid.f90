!===============================================================
!
!  Copyright (C) 2006 
!  This file is part of parsec, http://www.ices.utexas.edu/parsec/
!
!  This routine calculate the coefficients for the laplacian
!  in a general grid that can be non-orthogonal. in such a
!  grid, the laplacian can in general include mixed derivatives
!  (see for example: Modine et al. PRB 55, 10289 (1997) ).
!
!  the way to avoid that in a parallelpiped grid is to use 
!  derivatives in additional auxilary directions along nearest neighbors 
!  that has only pure 2nd derivative terms. A similar idea for 2d problems can 
!  be found in: Brandt et al. SIAM J. Sci. Comp. vol 21, 473-501, 
!  (1999) equation 4.2 and figure near it.       
!
!  instead of: f_xx+f_yy+f_zz we get for the laplacian on a non-cartesian grid:
!
!  a*f_xx+b*f_yy+c*f_zz+d*f_w1+e*f_w2+f*f_w3
!
!  where f_wi - are 2nd derivatives along new auxilary directions.
!
!  the coefficients a-f are calculated according to the
!  lattice vector matrix. 
!
!  2nd derivation along w=x+y direction for example translates 
!  in difference equations to:
!
!  sum(c(k)*grid_value(x+k,y+k,z)
!
!
!  according to how non-orthogonal the grid is, the following 
!  is done:
!
!  1. the grid angles are being checked, if the grid is 
!     orthogonal - nothing is done. if the grid has only 
!     one pair of non-orthogonal axes - exx_ggrid_2d is called.
!
!  2. finding of all the nearest neighbors on the grid is found
!
!  3. an iterative process that checks whether it is possible
!     to find an expression for the Laplacian is done.
!
!  4. the output of the subroutine is:
!
!     grid_lap_neig(3,3) - an array of neighbor directions
!     grid_lap_dir_num   - number of new directions
!     grid_b_lap(6)      - coefficients for the each direction
!                          1:3 - are the original directions
!                          4:6 - are the auxilary directions
!     grid_lap_dir(3)    - 1 or 0 , 1 means that the direction is active.
!
!   the exact equations for the procedure are described in:
!   Natan et al. PRB (2008, to be published).  
!
!  
!  Written by Amir Natan.     
!
!---------------------------------------------------------------
       subroutine exx_ggrid(latt_vec,grid_step,tol)
 
       use kinds, only  : dp     
       use cp_main_variables,    ONLY  : grid_lap_neig=>lap_neig, grid_lap_dir_num=>lap_dir_num, &
                                         grid_lap_dir_step=>lap_dir_step, grid_b_lap=>b_lap,     &
                                         grid_lap_dir=>lap_dir
       use constants

       implicit none
!
!  periodic boundary conditions data
       real(dp), intent(in) :: latt_vec(3,3)

       real(dp), intent(in) :: grid_step(3)
!  tolerance
       real(dp) , intent(in)  :: tol
!  
!  Work variables:     
!
       real(dp)  :: avec_tmp(3,3),ainv(3,3),ainvdot(3,3), adot(3,3)

       real(dp)  :: tmptr, tmpdet, tmpval(3), kval, tmpnorm

!  coeffecients of the laplacian in u,v,w coordinates
!  first 3 are pure derivatives coefficients: u, v, w
!  next are mixed derivatives uv, uw, and vw.
!    
       real(dp)  :: a_lap(6)

!  coefficients of the laplacian after adding the new 
!  first 3 are as a_lap, second are coefficients for
!  w1(u,v), w2(u,w), w3(v,w).      
       real(dp)  :: b_lap(6) 

!  flags for directions w1(u,v), w2(u,w), w3(v,w)
!  0 - direction is not needed
!  1 - direction is active
       integer   :: lap_dir(3)

!  step size for selected direction. 
       real(dp)  :: lap_dir_step(3)

!  temporary duplicates of grid step size. 
!  held twice to support easy calculation of new step size
       real(dp)  :: astep(3),bstep(3), rstep(3)

!  for the case of 2D - laplacian direction that is orthogonal
!  to the other 2, if there is no such direction value is 0.
       integer   :: lap_orth_dir

       integer i,j,k,idx1,op(200,4), disti(200), neig(20,3),flag, ngflag
       integer ndir ! number of mixed derivatives in laplacian of grid.
       real(dp) opg(200,3), opp(200,3), distp(200)
       real(dp) pneig(20,3), neig_un(20,3)
       real(dp) kv,kw
       real(dp) uu(3,3), mt(3,3), mtinv(3,3), bb(3), ff(3),vf

    
!---------------------------------------------------------------
!  First, calculation of normalized lattice vectors matrix is done
!  this is because we continue to work with usual distances on the
!  grid - so we need to normalize the lattice vectors to calculate
!  the correct metric.      

      avec_tmp = latt_vec  
    
      do i = 1, 3
         tmpnorm = sqrt(sum(avec_tmp(:,i)**2))
         avec_tmp(:,i) = avec_tmp(:,i)/tmpnorm
      enddo
    
      ainv = avec_tmp

!  Calculation of the inverse matrix.
      call mtrxin(ainv,tmpdet,tmptr)

!  Calculation of the matrix elements for the laplacian with mixed derivatives.

       ainvdot = matmul(ainv,transpose(ainv))
    
       a_lap(1) = ainvdot(1,1)
       a_lap(2) = ainvdot(2,2)
       a_lap(3) = ainvdot(3,3)
    
       a_lap(4) = ainvdot(1,2) + ainvdot(2,1)
       a_lap(5) = ainvdot(1,3) + ainvdot(3,1)
       a_lap(6) = ainvdot(2,3) + ainvdot(3,2)
    
       write(6,*) 'Original Laplacian coefficients are: ',a_lap

! checking if the grid is orthogonal or if one of the axes 
! is orthogonal to the other 2.

       ndir=0
       do i=1,3
          if(abs(a_lap(3+i))<tol) ndir=ndir+1
       enddo
      
       write(6,*) 'number of zero mixed derivatives is: ',ndir
      
       lap_orth_dir=0
      
       if(ndir==3) then ! means that the grid is orthogonal
          grid_b_lap(1:3)=1.d0
          grid_b_lap(4:6)=0.d0
          grid_lap_dir_num=0
          grid_lap_dir=0
          write(6,*) 'orthogonal grid found'
          return
       else
          if(ndir==2) then
             do i=1,3
                if(abs(a_lap(3+i))>tol) then
                   lap_orth_dir=3-i+1  ! note that this is possible because the way the 
                                   ! laplacian is arrange: f_uv, f_uw, f_vw
                   continue
                endif
             enddo
             call exx_ggrid_2d(latt_vec, grid_step, lap_orth_dir) 
             write(6,*) '2d non-orthogonal grid not implemented/debug yet'
             return
          endif
       endif
     
       write(6,*) 'laplacian orthogonal dir is :',lap_orth_dir
     
! now we calculate the 3 nearest neighbor directions - this 
! should be optimized sometime - it is now done in a very
! lazy/stupid way... - Amir. 

       kv=grid_step(2)/grid_step(1)
       kw=grid_step(3)/grid_step(2)
     
       op=0
       opg=0
       opp=0
       distp=0
       disti=0
     
       idx1=1
       do i=-1,1
        do j=-1,1
         do k=-1,1
           op(idx1,1)=i
           op(idx1,2)=j
           op(idx1,3)=k
           op(idx1,4)=idx1
           opg(idx1,1)=i
           opg(idx1,2)=kv*j
           opg(idx1,3)=kw*k
           call matvec3('N',avec_tmp,opg(idx1,1:3),opp(idx1,1:3))
           distp(idx1)=sqrt(dot_product(opp(idx1,1:3),opp(idx1,1:3)))
           disti(idx1)=idx1
           idx1=idx1+1
         enddo
        enddo
       enddo     
     
       call bubble_sort(27,distp(1:27),disti(1:27))
     
       i=1
       idx1=1

! the following loop assumes that neighbor displacement on each 
! coordinate is +/-1 at most. if this is changed it must be changed.

       ngflag=1
      
       do while((i<28) .and. (ngflag>0))
         i=i+1
         if(sum(abs(op(disti(i),1:3)))>1) then ! check that it is not one of the axis
           if(idx1==1) then
             neig(idx1,1:3)=op(disti(i),1:3)
             neig_un(idx1,1:3)=opg(disti(i),1:3)
             pneig(idx1,1:3)=opp(disti(i),1:3)
             idx1=idx1+1 
           else
             flag=0
             do j=1,idx1-1
              vf = sum((neig(j,1:3) + op(disti(i),1:3))**2)
              if (vf < 0.0001) flag = flag + 1
             enddo
             if(flag==0) then
               neig(idx1,1:3)=op(disti(i),1:3)
               neig_un(idx1,1:3)=opg(disti(i),1:3)
               pneig(idx1,1:3)=opp(disti(i),1:3)
               idx1=idx1+1
             endif
      
           endif
         endif
      
         if(idx1>3) then
     
          write(6,*) 'found the following 3 nearest neighbors'
          do k=1,3
             write(6,*) k, neig(k,1:3)
          enddo


       do k=1,3
        uu(1:3,k)=neig_un(k,1:3)/sqrt(dot_product(pneig(k,1:3),pneig(k,1:3)))
        lap_dir_step(k)=sqrt(dot_product(pneig(k,1:3),pneig(k,1:3)))* &
                       grid_step(1)
       enddo
   
       mt = transpose(uu)
       uu = mt
  
       write(6,*) 'uu matrix is:'
       write(6,*) transpose(uu)
   
       mt(1,:)=2*uu(:,1)*uu(:,2)
       mt(2,:)=2*uu(:,1)*uu(:,3)
       mt(3,:)=2*uu(:,2)*uu(:,3)
   
       write(6,*) 'mt matrix is:'
       write(6,*) transpose(mt)
   
       mtinv=mt
       
       call mtrxin1(mtinv,tmpdet,tmptr,flag)
       if(flag<0) then
         write(6,*) 'mt matrix is singular, replacing last neighbor'
         idx1=idx1-1
       else
         ngflag=0
         call matvec3('N',mtinv,a_lap(4:6),bb) ! bb are the coefficients
       endif
      endif
   
     enddo


!now building the matrix to find the coefficients for the different neighbors

       write(6,*) 'coeffs vec is:'
       write(6,*) bb
    
    
       do i=1,3
         b_lap(i)=a_lap(i)-sum(bb(1:3)*(uu(1:3,i)**2))
         lap_dir(i)=-1
       enddo
    
       b_lap(4:6)=bb



!  Now we have the final laplacian and also the number of 
!  directions we need to use. what is left is to arrange it in
!  a compact way (not taking into account directions that are
!  not needed) so we can save both computation and memory.

!  Calculating number of directions that are different from zero.
  
       do i=1,3
         if(abs(bb(i))<0.0000001) then
           lap_dir(i)=0
         else
           lap_dir(i)=1
         endif
       enddo
      
       write(6,*) 'bb is:' , bb
       write(6,*) 'lap_dir is:', lap_dir

       k = sum(abs(lap_dir(:)))
       grid_lap_dir_num = k
       grid_lap_dir = lap_dir
       grid_lap_dir_step = lap_dir_step
       grid_b_lap = b_lap
      
       do i=1,3
         grid_lap_neig(i,1:3)=neig(i,1:3)
       enddo
      
       write(*,'(1x,a,1x,i2)') 'Laplacian number of directions:', &
            grid_lap_dir_num
       write(*,*) 'Laplacian directions:'
       do i=1,3
          write(*,*) grid_lap_neig(i,1:3), '  weight is: ', grid_b_lap(3+i)
       enddo
       write(*,*) 'Laplacian coefficients data:'
       write(*,'(6(f5.2, 1x))') grid_b_lap
       write(*,*) 'lap_dir_step :'
       write(*,*) grid_lap_dir_step
  
     end subroutine exx_ggrid
!===============================================================

      subroutine bubble_sort(n, valvec,idxvec)
 
        use constants
        use kinds, only   :   DP
        implicit none
        integer, intent(in) :: n
        real(dp), intent(inout) :: valvec(1:n)
        integer, intent(inout) :: idxvec(1:n)

        integer :: i, j, k, itmp
        real(dp) :: minval, tmpval

        minval=valvec(1)
        itmp=1

        do j=1,n-1
          minval=valvec(j)
          itmp=idxvec(j)
          do i=j+1,n
            if(valvec(i)<minval) then
             itmp=idxvec(i)
             tmpval=valvec(i)
             idxvec(i)=idxvec(j)
             valvec(i)=valvec(j)
             idxvec(j)=itmp
                   valvec(j)=tmpval
             minval=tmpval
            endif
          enddo
        enddo

      end subroutine bubble_sort
  
      subroutine mtrxin1(m,det,tr,flag)

        use kinds, only   :   DP
        use constants
        implicit none

! Input/Output variables:
!
        real(dp), intent(inout) :: m(3,3)
        real(dp), intent(out) :: det,tr
        integer, intent(out) :: flag
!
! Work variables:
!
        real(dp) :: a(3,3),del,x
        integer i,j
!---------------------------------------------------------------
!
! compute matrix of cofactors
!

        flag=0

        a(1,1) = m(2,2)*m(3,3) - m(2,3)*m(3,2)
        a(2,1) = -m(2,1)*m(3,3) + m(2,3)*m(3,1)
        a(3,1) = m(2,1)*m(3,2) - m(2,2)*m(3,1)
        a(1,2) = -m(1,2)*m(3,3) + m(1,3)*m(3,2)
        a(2,2) = m(1,1)*m(3,3) - m(1,3)*m(3,1)
        a(3,2) = -m(1,1)*m(3,2) + m(1,2)*m(3,1)
        a(1,3) = m(1,2)*m(2,3) - m(1,3)*m(2,2)
        a(2,3) = -m(1,1)*m(2,3) + m(1,3)*m(2,1)
        a(3,3) = m(1,1)*m(2,2) - m(1,2)*m(2,1)
      
! compute determinant
!
        det = m(1,1)*a(1,1) + m(1,2)*a(2,1) + m(1,3)*a(3,1)
        tr = m(1,1) + m(2,2) + m(3,3)
        del = 1.0d-05
        if (abs(det) < del) then
          flag=-1
          return
        endif
!
! form mi
!
        do i=1,3
           do j=1,3
              x = a(i,j)/det
              m(i,j) = x
           enddo
        enddo
      
      end subroutine mtrxin1
!===============================================================

      subroutine mtrxin(m,det,tr)

      use kinds,    ONLY  : dp
      implicit none
!
! Input/Output variables:
!
      real(dp), intent(inout) :: m(3,3)
      real(dp), intent(out) :: det,tr
!
! Work variables:
!
      real(dp) :: a(3,3),del,x
      integer i,j
!---------------------------------------------------------------
!
! compute matrix of cofactors
!
      a(1,1) = m(2,2)*m(3,3) - m(2,3)*m(3,2)
      a(2,1) = -m(2,1)*m(3,3) + m(2,3)*m(3,1)
      a(3,1) = m(2,1)*m(3,2) - m(2,2)*m(3,1)
      a(1,2) = -m(1,2)*m(3,3) + m(1,3)*m(3,2)
      a(2,2) = m(1,1)*m(3,3) - m(1,3)*m(3,1)
      a(3,2) = -m(1,1)*m(3,2) + m(1,2)*m(3,1)
      a(1,3) = m(1,2)*m(2,3) - m(1,3)*m(2,2)
      a(2,3) = -m(1,1)*m(2,3) + m(1,3)*m(2,1)
      a(3,3) = m(1,1)*m(2,2) - m(1,2)*m(2,1)
!
! compute determinant
!
       det = m(1,1)*a(1,1) + m(1,2)*a(2,1) + m(1,3)*a(3,1)
       tr = m(1,1) + m(2,2) + m(3,3)
       del = 1.0d-05
       if (abs(det) < del) stop 501
!     
!      form mi
!
       do i=1,3
          do j=1,3
             x = a(i,j)/det
             m(i,j) = x
          enddo
       enddo

       end subroutine mtrxin
!===============================================================

!---------------------------------------------------------------
       subroutine matvec3(op,mat,vec,mvec)

       use kinds,     ONLY  : dp
       implicit none
!
!  Input/Output variables:
!
!  operation: transpose(M) or M
       character (len=1), intent(in) :: op
!  input matrix
       real(dp), intent(in) :: mat(3,3)
!  input vector V
       real(dp), intent(in) :: vec(3)
!  output vector MV = op(M)*V
       real(dp), intent(out) :: mvec(3)

!  Work variables:
!
       integer :: ii
       real(dp) :: vtmp(3)
!---------------------------------------------------------------
       do ii = 1, 3
          if (op == 'N') then
             vtmp(ii) = dot_product(mat(ii,:),vec)
          elseif (op == 'T') then
             vtmp(ii) = dot_product(mat(:,ii),vec)
          endif
       enddo
       mvec = vtmp
     
       return
       end subroutine matvec3
!===============================================================

!===============================================================
!
!  Copyright (C) 2006 
!  This file is part of parsec, http://www.ices.utexas.edu/parsec/
!
!  This routine calculate the coefficients for the laplacian
!  in a general grid that can be non-orthogonal on one plane. in such a
!  grid, the laplacian can in general include mixed derivatives
!  (see for example: Modine et al. PRB 55, 10289 (1997) ).
!
!  the way to avoid that in a parallelpiped grid is to use 
!  derivatives in directions such as: x+y, x+z, y+z or x-y, x-z
!  y-z. after some math it is possible to get to a laplacian
!  that has only pure 2nd derivative terms. A similar idea can 
!  be found in: Brandt et al. SIAM J. Sci. Comp. vol 21, 473-501, 
!  (1999) equation 4.2 and figure near it.       
!
!  this subroutine is called when there is only one plane where
!  the grid is not orthogonal, or equivalently, when one of the 
!  lattice directions is orthogonal to the other 2. 
!  In such a case only 1 additional auxilary direction is needed. 
!  
!  Written by Amir Natan.     
!
!---------------------------------------------------------------
       subroutine exx_ggrid_2d(latt_vec,grid_step,plane_num)
      
         use constants
         use kinds, only  : dp
         use cp_main_variables,    ONLY  : grid_lap_neig=>lap_neig, grid_lap_dir_num=>lap_dir_num, &
                                           grid_lap_dir_step=>lap_dir_step, grid_b_lap=>b_lap,     &
                                           grid_lap_dir=>lap_dir


         implicit none
!
!  Input/Output variables:
!
!  periodic boundary conditions data
         real(dp), intent(in) :: latt_vec(3,3)
!  grid structure
         real(dp), intent(in) :: grid_step(3)

!  plane number: 1=uv, 2=uw, 3=vw
!  the plane that includes the pair that is non-orthognal 
!  
         integer    :: plane_num
!  Work variables:     
!
         real(dp)  :: avec_tmp(3,3),ainv(3,3),ainvdot(3,3), adot(3,3)

         real(dp)  :: tmptr, tmpdet, tmpval(3), kval, tmpnorm

!  coeffecients of the laplacian in u,v,w coordinates
!  first 3 are pure derivatives coefficients: u, v, w
!  next are mixed derivatives uv, uw, and vw.
!    
         real(dp)  :: a_lap(6)

!  coefficients of the laplacian after adding the new 
!  first 3 are as a_lap, second are coefficients for
!  w1(u,v), w2(u,w), w3(v,w).      
         real(dp)  :: b_lap(6) 

!  flags for directions w1(u,v), w2(u,w), w3(v,w)
!  0 - direction is not needed
!  1 - u+v is used
!  -1 - u-v is used.      
         integer   :: lap_dir(3)

!  step size for selected direction. 
         real(dp)  :: lap_dir_step(3)

!  temporary duplicates of grid step size. 
!  held twice to support easy calculation of new step size
         real(dp)  :: astep(3),bstep(3), rstep(3)

         integer i,j,k, dflag(3), mflag(3)
    
!---------------------------------------------------------------
!  First, calculation of normalized lattice vectors matrix is done
!  this is because we continue to work with usual distances on the
!  grid - so we need to normalize the lattice vectors to calculate
!  the correct metric.      

         lap_dir=0
         dflag=0
         mflag=0

         avec_tmp = latt_vec

         do i = 1, 3
            tmpnorm = sqrt(sum(avec_tmp(:,i)**2))
            avec_tmp(:,i) = avec_tmp(:,i)/tmpnorm
         enddo

         ainv = avec_tmp

!  Calculation of the inverse matrix.
         call mtrxin(ainv,tmpdet,tmptr)

!  Calculation of the matrix elements for the laplacian with 
!  mixed derivatives.

         ainvdot = matmul(ainv,transpose(ainv))

         a_lap(1) = ainvdot(1,1)
         a_lap(2) = ainvdot(2,2)
         a_lap(3) = ainvdot(3,3)

         a_lap(4) = ainvdot(1,2) + ainvdot(2,1)
         a_lap(5) = ainvdot(1,3) + ainvdot(3,1)
         a_lap(6) = ainvdot(2,3) + ainvdot(3,2)

         select case(plane_num)
          case (3)
           dflag(1)=1
           dflag(2)=1
           mflag(1)=1
           mflag(2)=-1
          case (2)
           dflag(1)=1
           dflag(3)=1
           mflag(1)=1
           mflag(3)=-1
          case (1)
           dflag(2)=1
           dflag(3)=1
           mflag(2)=1
           mflag(3)=-1
         end select 
   

!  Next a check is done for each direction for the following:
!
!  1. if the 2 coordinates are orthogonal (for example u and v)
!  No mixed derivative element exist and hence no need
!  for derivation in a combined direction (like u+v or u-v).      
!  
!  2. a check is done whether it is better to use u+v or u-v 
!     this according to if the angle between u and v is bigger than
!     90 (u+v is used) or smaller than 90 (u-v is used).
!   
!  3. calculation of grid step in the selected direction
!     at the moment - it is assumed that the grid steps
!     on the main axis (u,v,w) are EQUAL!!!!!!!!!!!!!!
!     this will be fixed in the future.
!

         adot = matmul(transpose(avec_tmp),avec_tmp)

!  Check the u-v direction.

!  adot(1,2) is actually equal to the dot product of u and v, 
!  since they were normalized - it is exactly cos(tetha)!       

         tmpval(1) = adot(1,2) ! dot product of u and v
         tmpval(2) = adot(1,3) ! dot product of u and w
         tmpval(3) = adot(2,3) ! dot product of v and w

         astep(1) = grid_step(1)
         astep(2) = grid_step(1)
         astep(3) = grid_step(2)
       
         bstep(1) = grid_step(2)
         bstep(2) = grid_step(3)
         bstep(3) = grid_step(3)

!  rstep definition is to calculate the ratio
!  between grid%step in the different laplacians directions
!  if this is different than 1, special care should be taken
!  because u+v by itself is not good because it is not an 
!  actual grid direction. instead one should take
!  u+r*v or u-r*v where r is the ratio between the 
!  step in v to the step in u. 

         rstep(:) = bstep(:)/astep(:) 

!  all the remarks here refer to u and v as an example.

        i=3-plane_num+1

        if(tmpval(i) > 0) then   ! cos(theta)>0 --> theta < 90
           !           choosing to do derivation along u-v.               
           lap_dir(1) = -1
           !           calculating the grid step for the new direction
           tmpnorm = astep(i)**2+bstep(i)**2- 2*astep(i)*bstep(i)*tmpval(i)
           lap_dir_step(1) = sqrt(tmpnorm)
           !           kval=a_uv*|u-v|^2=a_uv*2*(1-cos(theta))               
           kval = a_lap(3+i)*(1+rstep(i)**2-2*rstep(i)*tmpval(i))
           b_lap(4) = -kval/2/rstep(i)
        else                ! cos(theta)<0 --> theta > 90
           !           choosing to do derivation along u+v               
           lap_dir(1) = 1
           !           calculating the grid step for the new direction
           tmpnorm = astep(i)**2+bstep(i)**2+ &
                2*astep(i)*bstep(i)*tmpval(i)
           lap_dir_step(1) = sqrt(tmpnorm)
           !           kval=a_uv*|u+v|^2=a_uv*2*(1+cos(theta))     
           kval = a_lap(3+i)*(1+rstep(i)**2+2*rstep(i)*tmpval(i))
           b_lap(4) = kval/2/rstep(i)
        end if

        b_lap(1) = a_lap(1)- &
          dflag(1)*(dflag(2)+dflag(3))*lap_dir(1)*a_lap(3+i)/2/rstep(i)

        b_lap(2) = a_lap(2)- &
          dflag(2)*dflag(1)*rstep(i)*lap_dir(1)*a_lap(3+i)/2- &
          dflag(2)*dflag(3)*lap_dir(1)*a_lap(3+i)/2/rstep(i)

        b_lap(3) = a_lap(3)- &
          dflag(3)*(dflag(1)+dflag(2))*rstep(i)*lap_dir(1)*a_lap(3+i)/2

!  Now we have the final laplacian and also the number of 
!  directions we need to use. what is left is to arrange it in
!  a compact way (not taking into account directions that are
!  not needed) so we can save both computation and memory.

!  Calculating number of directions that are different
!  from zero.

        k = sum(abs(lap_dir(:)))

        grid_lap_dir_num = k
      
        grid_lap_dir = lap_dir
       
        grid_lap_dir_step = lap_dir_step
      
        grid_b_lap = b_lap
      

        if(lap_dir(1)>0) then 
          grid_lap_neig(1,:)=dflag
        else 
          grid_lap_neig(1,:)=mflag
        endif
      
        grid_lap_neig(2,:)=0
        grid_lap_neig(3,:)=0

     write(7,'(1x,a,1x,i2)') 'Laplacian number of directions:', &
          grid_lap_dir_num
     write(7,*) grid_lap_neig(1,1:3), '  weight is: ', grid_b_lap(4)

     write(7,*) 'Laplacian coefficients data:'
     write(7,'(6(f5.2, 1x))') grid_b_lap

      end subroutine exx_ggrid_2d
!===============================================================
