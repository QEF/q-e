
!====================================================================
!  This subroutine finds for each state the valence-state neighbors 
!====================================================================
      SUBROUTINE exx_index_pair_nv(wc, overlap, nj, nj_max)
      USE kinds,                 ONLY  : DP
      USE electrons_base,        ONLY  : nbsp
      USE cp_main_variables,     ONLY  : vwc
      USE mp_global,             ONLY  : me_image, intra_image_comm, nproc_image
      USE cell_base,             ONLY  : h, ainv, r_to_s, s_to_r
      USE parallel_include
      USE wannier_base,          ONLY  : neigh, dis_cutoff, vnbsp
      IMPLICIT NONE
      
! wc is the wannier centers of the initial quasi-particle states
      REAl(DP),INTENT(IN)    ::    wc(3, nbsp)
! number and index of overlapping orbitals for all conduction states
      INTEGER, INTENT(INOUT) ::    overlap(neigh,nbsp), nj(nbsp) , nj_max
      
      INTEGER     i, j, ii, ierr
      REAl(DP)    ri(3), rj(3), ris(3), rjs(3) , rijs(3), rij(3), distance

!==================================================================

!     print *, 'entering exx_index_pair_nv', dis_cutoff, neigh, vnbsp
      overlap(:,:) = 0

      do i = 1, nbsp

         nj(i) = 0

         do ii = 1, 3
            ri(ii) = wc( ii, i )
         enddo

         call r_to_s(ri, ris, ainv)

         do j = 1, vnbsp

            do ii = 1, 3
               rj(ii) = vwc(ii,j)
            enddo

            call r_to_s(rj, rjs, ainv)

            do ii = 1, 3
               rijs(ii) = rjs(ii) - ris(ii) - ANINT(rjs(ii) - ris(ii))
            enddo

            call s_to_r(rijs, rij, h)

            distance = sqrt( rij(1)*rij(1) + rij(2)*rij(2) + rij(3) *rij(3) )

            if( distance .lt. dis_cutoff )then
                nj(i) = nj(i) + 1
                if(nj(i) > neigh)then
                   print *, 'increase neigh, stop in exx_pair', nj(i), neigh
                   return
                endif
                overlap(nj(i),i) = j
            end if
         enddo  ! j=1,vnbsp
      enddo  

      nj_max = nj(1)
      do i = 2, nbsp
         if(nj(i) > nj_max) nj_max = nj(i)
      enddo

      return
!==================================================================

      END

!==========================================================================
! 
!whether this algorithm works for every system needs tested -- Zhaofeng
!
!==========================================================================

     SUBROUTINE exx_index_pair(wannierc, overlap3, nj, nj_max )
      USE kinds,                   ONLY  : DP
      USE electrons_base,          ONLY  : nbsp, nspin, nupdwn, iupdwn
      USE mp_global,               ONLY  : me_image, intra_image_comm
      USE cell_base,               ONLY  : ainv, h, r_to_s, s_to_r
      USE parallel_include
      USE mp,                      ONLY  : mp_barrier 
      USE wannier_base,            ONLY  : neigh, dis_cutoff
      IMPLICIT NONE
      
      REAl(DP),INTENT(IN)    ::    wannierc(3, nbsp)
      INTEGER, INTENT(INOUT) ::    overlap3(neigh/2,nbsp), nj(nbsp) 
      INTEGER, INTENT(INOUT) ::    nj_max
      
      REAl(DP), ALLOCATABLE  ::    distance(:)
      INTEGER,  ALLOCATABLE  ::    overlap(:, :), overlap2(:)
      INTEGER     i, j, k, ii, jj, ip, ir, ierr, num, num1, nj_avg
      REAl(DP)    ri(3), rj(3), rij(3), ris(3), rjs(3), rijs(3)

      INTEGER   gindx_of_i, jj_neib_of_i
      INTEGER   spin_loop, nupdwn_(2),iupdwn_(2)

!==================================================================

      ALLOCATE( distance(nbsp) )
      ALLOCATE( overlap(neigh, nbsp), overlap2(neigh) )

      if(nspin == 1)then
         nupdwn_(1)=nbsp
         nupdwn_(2)=nbsp
         iupdwn_(1)=1
         iupdwn_(2)=1
         spin_loop=1
      else
         nupdwn_(:)=nupdwn(:)
         iupdwn_(:)=iupdwn(:)
         print *, "nbsp is",nbsp
         print *, "nupdwn, iupdwn",nupdwn_(1), nupdwn_(2),iupdwn_(1), iupdwn_(2)
      end if

      overlap(:,:) = 0

      do i = 1, nbsp

         do ii = 1, 3
            ri(ii) = wannierc( ii, i )
         enddo
 
         call r_to_s(ri, ris, ainv)

         nj(i) = 0
         distance(:) = 0.d0

         if( i < iupdwn(2) ) then
             spin_loop = 1
         else
             spin_loop = 2
         end if

         do j = iupdwn_(spin_loop), iupdwn_(spin_loop) + nupdwn_(spin_loop)-1

            if( j .NE. i )then
               do ii = 1, 3
                  rj(ii) = wannierc(ii, j)
               enddo
               call r_to_s(rj, rjs, ainv)

               do ii = 1, 3
                  rijs(ii) = rjs(ii) - ris(ii) - ANINT(rjs(ii) - ris(ii))
               enddo

               call s_to_r(rijs, rij, h)

               distance(j) = sqrt( rij(1)*rij(1) + rij(2)*rij(2) + rij(3) *rij(3) )

               if( distance(j) .lt.  dis_cutoff )then
                   nj(i) = nj(i) + 1
                   if(nj(i) .eq. 1 )then
                      overlap(nj(i),i) = j
                   else if(distance(overlap(nj(i)-1,i)) .le. distance(j))then
                      overlap(nj(i),i)=j
                   else
                      overlap2(:)=0
                      do ir=1,nj(i)-1
                         if(distance(overlap(ir,i)) < distance(j))then
                            overlap2(ir)=overlap(ir,i)
                         else
                            overlap2(ir)=j
                            do ip=ir+1,nj(i)
                               overlap2(ip)=overlap(ip-1,i)
                            end do
                            GO TO 555
                         end if
                      end do
       555            continue
                      do ir = 1, nj(i)
                         overlap(ir,i)=overlap2(ir)
                      enddo
                   end if
                end  if !if for distance(j)<5.0d0
             end if
          end do  ! j=1,nbsp
       enddo !i = 1, nbsp

!==================================================================

     num  = 0
     num1 = 0
     overlap3(:,:) = 0

     do j = 1, neigh/2
        do i = 1, nbsp
           do jj = 1, nj(i)

              jj_neib_of_i = overlap(jj,i)
              if(jj_neib_of_i .gt. 0)then
                 overlap3(j, i) = jj_neib_of_i
                 overlap(jj,i) = 0
                 num = num + 1
                 do k = 1, nj(jj_neib_of_i)
                    if(overlap(k,jj_neib_of_i) .eq. i)then
                       overlap(k,jj_neib_of_i) = 0
                       num1 = num1 + 1
                       go to   666
                    end if
                 end do

              end if
           end do
   666     continue
        end do !i=1,nbsp
     end do  ! j=1,neigh/2

!==========================================================

     print *,"pair num  and num1 is",num,num1
     if(num .ne. num1)call errore("exx", "sth wrong with overlap",1)
!--------------------------------------------
! njj
     do i = 1,nbsp
        num = 0
        do j = 1,neigh/2
           if(overlap3(j,i) .ne. 0)num = num+1
        end do
        nj(i) = num
     end do

     if(me_image .eq. 0)then
!      open(unit=20,file='pair.dat',status='unknown',form='formatted')
       do i = 1, nbsp
          write(*, '(8I7)')(overlap3(j,i),j=1, nj(i))
       enddo
!      close(20)
     endif

!---------------------------------------------------------------------
!    do j=1,neigh/2
!       num=0
!       do i=1,nbsp
!          if(overlap3(j,i) .ne. 0)num=num+1
!       end do
!       if(num .gt. 0)  print *,"the j low num is",num,j
!    end do
!---------------------------------------------------------------------
     nj_max=1
     nj_avg=0
     do i = 1 , nbsp
        if(nj(i) .GT. nj_max) nj_max = nj(i)
        nj_avg = nj_avg+nj(i)
     end do
     print *,"nj_max and nj_avg are ",nj_max, nj_avg/nbsp
     write(*,'(10I5)')(nj(i),i=1, nbsp)

     DEALLOCATE( distance)
     DEALLOCATE( overlap, overlap2 )

     return
     END
