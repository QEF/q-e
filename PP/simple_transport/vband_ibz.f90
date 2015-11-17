      subroutine vband_ibz ( nk1,nk2,nk3, nbnd, nksfit, etk, eqk, at, vk, dfk )
      !
      ! Written by Burak Himmetoglu (burakhmmtgl@gmail.com)
      ! Uses some parts of the PW distribution
      !
      ! Computes band velocities and forward derivatives using finite
      ! differences.
      !

!$    use omp_lib

      implicit none

      integer, intent(in) :: nk1,nk2,nk3, nbnd, nksfit, eqk(nk1*nk2*nk3)

      double precision, intent(in) :: etk(nbnd,nksfit), at(3,3)

      double precision, intent(out) :: vk(nbnd,nk1*nk2*nk3,3), dfk(nbnd,nk1*nk2*nk3,3)
                                       ! vk is band velocity
                                       ! dfk is forward derivative (used in adaptive smearing)

      integer :: i,j,k,n,ik,ibnd,nktot,nm,np

      double precision :: temp, temp2, vaux(3,nk1*nk2*nk3)

      ! Create the k-points, then compute velocities   

      nktot = nk1*nk2*nk3
    
      !$omp parallel do default(none) &
      !$omp collapse(4) &
      !$omp private(ibnd,i,j,k,n,np,nm) &
      !$omp shared(nbnd,nk1,nk2,nk3,vk,dfk,etk,eqk)
      do ibnd=1,nbnd
         do i=1,nk1
            do j=1,nk2
               do k=1,nk3
                  n = (k-1) + (j-1)*nk3 + (i-1)*nk2*nk3 + 1
                  !   
                  ! velocity along k1
                  !
                  if ( i .eq. 1 ) then
                     np = (k-1) + (j-1)*nk3 + (i)*nk2*nk3 + 1
                     nm = (k-1) + (j-1)*nk3 + (i-2+nk1)*nk2*nk3 + 1
                  else if (i .eq. nk1) then
                     np = (k-1) + (j-1)*nk3 + (i-nk1)*nk2*nk3 + 1
                     nm = (k-1) + (j-1)*nk3 + (i-2)*nk2*nk3 + 1
                  else
                     np = (k-1) + (j-1)*nk3 + (i)*nk2*nk3 + 1
                     nm = (k-1) + (j-1)*nk3 + (i-2)*nk2*nk3 + 1
                  end if
                  vk(ibnd,n,1) = ( etk(ibnd,eqk(np))-etk(ibnd,eqk(nm)) )/(2.0/nk1)
                  dfk(ibnd,n,1) = ( etk(ibnd,eqk(np))-etk(ibnd,eqk(n)) )/(2.0/nk1)
                  !
                  ! velocity along k2
                  !
                  if ( j .eq. 1 ) then
                     np = (k-1) + (j)*nk3 + (i-1)*nk2*nk3 + 1
                     nm = (k-1) + (j-2+nk2)*nk3 + (i-1)*nk2*nk3 + 1
                  else if (j .eq. nk2) then
                     np = (k-1) + (j-nk2)*nk3 + (i-1)*nk2*nk3 + 1
                     nm = (k-1) + (j-2)*nk3 + (i-1)*nk2*nk3 + 1
                  else
                     np = (k-1) + (j)*nk3 + (i-1)*nk2*nk3 + 1
                     nm = (k-1) + (j-2)*nk3 + (i-1)*nk2*nk3 + 1
                  end if
                  vk(ibnd,n,2) = ( etk(ibnd,eqk(np))-etk(ibnd,eqk(nm)) )/(2.0/nk2)
                  dfk(ibnd,n,2) = ( etk(ibnd,eqk(np))-etk(ibnd,eqk(n)) )/(2.0/nk2)
                  !
                  ! velocity along k3
                  !
                  if ( k .eq. 1 ) then
                     np = (k) + (j-1)*nk3 + (i-1)*nk2*nk3 + 1
                     nm = (k-2+nk3) + (j-1)*nk3 + (i-1)*nk2*nk3 + 1
                  else if (k .eq. nk3) then
                     np = (k-nk3) + (j-1)*nk3 + (i-1)*nk2*nk3 + 1
                     nm = (k-2) + (j-1)*nk3 + (i-1)*nk2*nk3 + 1
                  else
                     np = (k) + (j-1)*nk3 + (i-1)*nk2*nk3 + 1
                     nm = (k-2) + (j-1)*nk3 + (i-1)*nk2*nk3 + 1
                  end if
                  vk(ibnd,n,3) = ( etk(ibnd,eqk(np))-etk(ibnd,eqk(nm)) )/(2.0/nk3)
                  dfk(ibnd,n,3) = ( etk(ibnd,eqk(np))-etk(ibnd,eqk(n)) )/(2.0/nk3)
                  !
               end do  ! k
            end do ! j
         end do ! i
      end do ! ibnd
      !$omp end parallel do
      !
      ! Convert the derivatives in crystal coordinates to cartesian coordinates
      !
      do ibnd=1,nbnd
         !
         do j=1,3
            vaux(j,:) = vk(ibnd,:,j) ! temporary vector
         end do
         !
         ! vaux is updated after this call
         !
         call cryst_to_cart (nktot,vaux,at,1)
         !
         do j=1,3
            vk(ibnd,:,j) = vaux(j,:) 
         end do
         !  
      end do
      !
      end subroutine vband_ibz
