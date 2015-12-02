      subroutine lint ( nsym, s, minus_q, at, bg, npk, k1,k2,k3,        &
     &            nk1,nk2,nk3, nks, xk, kunit, nkBZ, eqBZ, sBZ)

      implicit none
      integer, intent (IN) :: nks, nsym, s(3,3,48), npk, k1, k2, k3,    &
                              nk1, nk2, nk3, kunit, nkBZ
      logical, intent (IN) :: minus_q ! use .true.
      !
      double precision, intent(IN):: at(3,3), bg(3,3), xk(3,npk)
      integer, intent(OUT) :: eqBZ(nkBZ), sBZ(nkBZ)
      !
      double precision :: xkr(3), deltap(3), deltam(3)
      double precision, parameter:: eps=1.0d-4
      !
      double precision, allocatable :: xkg(:,:), xp(:,:)
      integer ::  i,j,k, ns, n, nk
      integer :: nkh
      !
      ! Re-generate a uniform grid of k-points xkg
      !
      allocate (xkg( 3,nkBZ))
      !
      ! if(kunit < 1 .or. kunit > 2) call errore('lint','bad kunit value',kunit)
      !
      ! kunit=2: get only "true" k points, not k+q points, from the list
      !
      !
      nkh = nks/kunit
      !
      allocate (xp(3,nkh))
      !
      if (kunit == 1) then
         xp(:,1:nkh) = xk(:,1:nkh)
      else
         do j=1,nkh
            xp(:,j) = xk(:,2*j-1)
         enddo
      end if
      !write(6,*) 'nkh  ', 'nkBZ  ', 'nsym'
      !write(6,*) nkh, nkBZ, nsym
      ! 
      ! Generate a uniform mesh
      !
      do i=1,nk1
         do j=1,nk2
            do k=1,nk3
               n = (k-1) + (j-1)*nk3 + (i-1)*nk2*nk3 + 1
               xkg(1,n) = dble(i-1)/nk1 + dble(k1)/2/nk1
               xkg(2,n) = dble(j-1)/nk2 + dble(k2)/2/nk2
               xkg(3,n) = dble(k-1)/nk3 + dble(k3)/2/nk3
            end do
         end do
      end do
      !
      !
      call cryst_to_cart (nkh,xp,at,-1)
      !
      !Maybe add here some OpenMP clauses
      do nk=1,nkBZ
         do n=1,nkh
            do ns=1,nsym
               do i=1,3
                  !
                  xkr(i) = s(i,1,ns) * xp(1,n) +  s(i,2,ns) * xp(2,n) + &
     &                     s(i,3,ns) * xp(3,n)
                  !
               end do
               !
               do i=1,3
                  deltap(i) = xkr(i)-xkg(i,nk) - nint (xkr(i)-xkg(i,nk))
                  deltam(i) = xkr(i)+xkg(i,nk) - nint (xkr(i)+xkg(i,nk))
               end do
               !
               if ( sqrt ( deltap(1)**2 + deltap(2)**2 +                &
     &                      deltap(3)**2 ) < eps .or. ( minus_q .and.   &
     &              sqrt ( deltam(1)**2 + deltam(2)**2 +                &
     &                      deltam(3)**2 ) < eps ) ) then
                  ! 
                  eqBZ(nk) = n
                  sBZ(nk) = ns
                  !
                  go to 15
               end if
               ! 
            end do
         end do
         !
         !call errore('lint','cannot locate  k point  xk',nk)
         write(6,*) 'somethings wrong', nk
         ! 
15       continue
      end do
      !
      do n=1,nkh
         do nk=1,nkBZ
            if (eqBZ(nk) == n) go to 20
         end do
         !  this failure of the algorithm may indicate that the displaced grid
         !  (with k1,k2,k3.ne.0) does not have the full symmetry of the lattice
         !
         !call errore('lint','cannot remap grid on k-point list',n)
         write(6,*) 'somethings wrong-2', nk, n
20       continue
      end do

      deallocate(xkg)
      deallocate(xp)

      return
      end subroutine lint

